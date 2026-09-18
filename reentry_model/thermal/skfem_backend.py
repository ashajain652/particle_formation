"""scikit-fem backend: P1 tetrahedra, backward Euler, Picard on k(T)/c_p(T) with Newton on the radiation term.

Assembly. The P1 element stiffness and mass matrices depend on T only through one coefficient per element
(k(T_e) and rho c_p(T_e) at the element-mean temperature), so the unit-coefficient element matrices are computed
once and rescaled into a fixed CSR pattern on every iteration (~20 ms for 50 k tets). scikit-fem's generic `asm`
on the same MeshTet/ElementTetP1 costs 0.2 s per iteration at 12 k nodes (measured 2026-09-18) and is kept as the
reference in `reference_operators` for the conformance test.

Radiation. The boundary functional uses the facet-mean temperature: eps sigma (T_f^4 - T_amb^4) A_f/3 to each of
the facet's three nodes, with its exact Jacobian 4 eps sigma T_f^3 A_f/9 on every node pair of the facet. Each
Newton iterate solves the symmetric positive definite system
    (M/dt + K + B_k) T = M T_old/dt + F_conv - F_rad(T_k) + B_k T_k
with M, K frozen at the previous iterate; because the rows of K sum to zero the discrete energy balance
1^T M (T - T_old) = dt (Q_conv - Q_rad) holds to the Newton tolerance.

Linear solver. `direct`: SciPy SuperLU; `amg` (default): CG preconditioned by a pyamg smoothed-aggregation
hierarchy rebuilt every `amg_rebuild_every` solves (the operator changes slowly). Measured on the 100 mm
2 mm/8 mm mesh (18.9 k nodes): AMG 0.03 s per solve vs SuperLU 0.6 s."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import SIGMA_SB, StepResult


class _Pattern:
    """Fixed CSR sparsity for repeated (rows, cols) entries; `assemble` sums values into it with one bincount."""

    def __init__(self, rows, cols, n):
        key = rows.astype(np.int64) * n + cols.astype(np.int64)
        unique, self._map = np.unique(key, return_inverse=True)
        self._indptr = np.searchsorted(unique // n, np.arange(n + 1))
        self._indices = (unique % n).astype(np.int32)
        self.n, self.nnz = n, unique.size

    def assemble(self, values):
        data = np.bincount(self._map, weights=np.asarray(values).ravel(), minlength=self.nnz)
        return sp.csr_matrix((data, self._indices, self._indptr), shape=(self.n, self.n))


def element_matrices(points, tets):
    """Per-element volumes, unit stiffness (grad phi_a . grad phi_b V_e) and consistent mass (V_e/20 (1 + delta_ab))."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    grad_ref = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    grads = np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), grad_ref)
    Ke = np.einsum("eai,ebi->eab", grads, grads) * vol[:, None, None]
    Me = (np.ones((4, 4)) + np.eye(4))[None] * (vol / 20.0)[:, None, None]
    return vol, Ke, Me


class SkfemThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=False, newton_tol=1e-6, max_iterations=30,
                 amg_rebuild_every=30, cg_tol=1e-10):
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        self.linear_solver, self.lumped_mass = linear_solver, lumped_mass
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.amg_rebuild_every, self.cg_tol = amg_rebuild_every, cg_tol
        self._ml, self._solves, self.last_cg_iterations = None, 0, 0

    def setup(self, mesh, material, emissivity):
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        self.points, self.tets = mesh.points, mesh.tets
        surface = mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        self.vol, self.Ke, self.Me = element_matrices(self.points, self.tets)
        n = len(self.points)
        self.pattern = _Pattern(np.repeat(self.tets, 4, axis=1).ravel(), np.tile(self.tets, (1, 4)).ravel(), n)
        self.facet_pattern = _Pattern(np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel(), n)
        self.Bf = np.ones((3, 3))[None] * (self.areas / 9.0)[:, None, None]
        self.T, self._T_prev = np.full(n, 300.0), None

    def set_temperature(self, T):
        self.T = np.full(len(self.points), float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    def facet_temperature(self, T=None):
        return (self.T if T is None else T)[self.faces].mean(axis=1)

    def facet_load(self, q):
        """Nodal load vector of a per-facet flux q [W/m2]: A_f/3 to each of the facet's nodes."""
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=len(self.points))

    def operators(self, T, T_old=None):
        """Stiffness K with k at the element-mean temperature, and mass M with the secant heat capacity
        [h(T_e) - h(T_old,e)] / (T_e - T_old,e) (c_p_eff where the element has not moved), so that at convergence
        1^T M (T - T_old) equals the enthalpy increment exactly, whatever h(T) is (Step 3's latent heat included)."""
        Te = T[self.tets].mean(axis=1)
        K = self.pattern.assemble(self.Ke * self.material.k(Te)[:, None, None])
        c = self.material.cp_eff(Te)
        if T_old is not None:
            Te_old = T_old[self.tets].mean(axis=1)
            dT = Te - Te_old
            moved = np.abs(dT) > 1e-9
            c = np.where(moved, (self.material.enthalpy(Te) - self.material.enthalpy(Te_old)) / np.where(moved, dT, 1.0), c)
        M = self.pattern.assemble(self.Me * (self.material.rho * c)[:, None, None])
        if self.lumped_mass:
            M = sp.diags(np.asarray(M.sum(axis=1)).ravel(), format="csr")
        return K, M

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def energy(self):
        """Stored enthalpy above T_REF: sum_e rho V_e h(T_e); equals 1^T M T for constant c_p (consistent or lumped mass)."""
        Te = self.T[self.tets].mean(axis=1)
        return float((self.material.rho * self.vol * self.material.enthalpy(Te)).sum())

    def step(self, dt, q_conv, T_amb, dirichlet=None):
        T_old = self.T
        # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        T_new, iteration = T_k, 0
        for iteration in range(1, self.max_iterations + 1):
            K, M = self.operators(T_k, T_old)
            Tf = self.facet_temperature(T_k)
            B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            A = (M / dt + K + B).tocsr()
            b = M @ T_old / dt + F_conv - F_rad + B @ T_k
            T_new = self._solve_with_dirichlet(A, b, T_k, dirichlet)
            if np.linalg.norm(T_new - T_k) <= self.newton_tol * np.linalg.norm(T_new):
                break
            T_k = T_new
        self._T_prev, self.T = T_old, T_new
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration)

    def _solve_with_dirichlet(self, A, b, x0, dirichlet):
        if dirichlet is None:
            return self._solve(A, b, x0)
        nodes, values = dirichlet
        free = np.ones(A.shape[0], dtype=bool)
        free[nodes] = False
        x = np.zeros(A.shape[0])
        x[nodes] = values
        x[free] = self._solve(A[free][:, free], b[free] - A[free][:, ~free] @ x[~free], x0[free], fresh=True)
        return x

    def _solve(self, A, b, x0, fresh=False):
        if self.linear_solver == "direct":
            return spla.spsolve(A.tocsc(), b)
        if fresh or self._ml is None or self._solves % self.amg_rebuild_every == 0:
            import pyamg
            self._ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
        self._solves += 1
        counter = []
        x, info = spla.cg(A, b, x0=x0, rtol=self.cg_tol, maxiter=500, M=self._ml.aspreconditioner(cycle="V"),
                          callback=lambda _: counter.append(1))
        self.last_cg_iterations = len(counter)
        if info != 0:
            raise RuntimeError("CG did not converge (info {})".format(info))
        return x

    def reference_operators(self, T):
        """K and M assembled by scikit-fem with the same element-mean coefficients (conformance test only)."""
        from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP1, MeshTet, asm
        from skfem.helpers import dot, grad
        basis = Basis(MeshTet(self.points.T.copy(), self.tets.T.copy()), ElementTetP1())
        basis0 = basis.with_element(ElementTetP0())
        Te = T[self.tets].mean(axis=1)

        @BilinearForm
        def stiffness(u, v, w):
            return w.k * dot(grad(u), grad(v))

        @BilinearForm
        def mass(u, v, w):
            return w.c * u * v

        K = asm(stiffness, basis, k=basis0.interpolate(self.material.k(Te)))
        M = asm(mass, basis, c=basis0.interpolate(self.material.rho * self.material.cp_eff(Te)))
        return K, M
