"""scikit-fem backend: P1 tetrahedra, backward Euler, Newton on the nodal enthalpy with the radiation term.

Assembly. The P1 stiffness depends on T through one coefficient per element (k at the element-mean temperature)
and the mass through a nodal coefficient (rho c at the nodes, interpolated linearly), so the unit-coefficient
element matrices are computed once and rescaled into a fixed CSR pattern on every iteration (~20 ms for 50 k
tets). scikit-fem's generic `asm` on the same MeshTet/ElementTetP1 costs 0.2 s per iteration at 12 k nodes
(measured 2026-09-18) and is kept as the reference in `reference_operators` for the conformance test. Element
fractions phi_e (Step 3) scale both coefficients; nodes without material are pinned at their temperature.

Time stepping. Each iterate solves the tangent system of the residual R(T) = E(T)/dt + K T - F_conv + F_rad(T) -
F_extra, where E = M(rho c_sec)(T - T_old) is the exact nodal enthalpy increment (secant heat capacity per node, so
1^T E equals the change of rho int h dV with latent heat included); the tangent uses c_p,eff(T_i) per node and the
radiation Jacobian. The update is mapped through the true h(T) per node (enthalpy-consistent update) so that a node
cannot jump across the melting range, with damping as a fallback. With a smooth c_p the scheme is the Step 2
secant-capacity iteration unchanged; measured 2026-09-20: 2 iterations per step without melting, 3-5 with.

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


def nodal_mass_weights():
    """W[k, a, b] = int phi_a phi_b phi_k dV / V_e on a tetrahedron: 1/20 (a = b = k), 1/60 (two equal), 1/120 (all
    distinct); sum_k W[k] is the consistent mass V_e/20 (1 + delta_ab)."""
    W = np.full((4, 4, 4), 1.0 / 120.0)
    for k in range(4):
        for a in range(4):
            W[k, a, a] = 1.0 / 60.0
            W[k, k, a] = W[k, a, k] = 1.0 / 60.0
        W[k, k, k] = 1.0 / 20.0
    return W


def element_matrices(points, tets):
    """Per-element volumes, unit stiffness (grad phi_a . grad phi_b V_e) and the mass matrices for a nodal coefficient:
    M_e(c) = sum_k c_k Mk[e, k] (Mk = V_e W), so that M_e(1) is the consistent mass V_e/20 (1 + delta_ab)."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    grad_ref = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    grads = np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), grad_ref)
    Ke = np.einsum("eai,ebi->eab", grads, grads) * vol[:, None, None]
    Mk = nodal_mass_weights()[None] * vol[:, None, None, None]
    return vol, Ke, Mk


class SkfemThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=True, newton_tol=1e-6, max_iterations=30,
                 amg_rebuild_every=30, cg_tol=1e-10):
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.lumped_mass = linear_solver, lumped_mass
        self.consistent_mass = not lumped_mass
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.amg_rebuild_every, self.cg_tol = amg_rebuild_every, cg_tol
        self._ml, self._solves, self.last_cg_iterations, self.direct_fallbacks = None, 0, 0, 0

    def setup(self, mesh, material, emissivity):
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        self.points, self.tets = mesh.points, mesh.tets
        self.vol, self.Ke, self.Mk = element_matrices(self.points, self.tets)
        n = len(self.points)
        self.pattern = _Pattern(np.repeat(self.tets, 4, axis=1).ravel(), np.tile(self.tets, (1, 4)).ravel(), n)
        self.T, self._T_prev = np.full(n, 300.0), None
        self.phi = np.where(mesh.active, 1.0, 0.0)
        self.film_mass = np.zeros(n)              # kg per node: the melt film riding on the boundary (Step 3)
        self.pinned = np.zeros(n, dtype=bool)
        self._refresh_surface()

    def _refresh_surface(self):
        """Facet data (loads, radiation pattern) of the mesh's current boundary and the pinned (material-free) nodes."""
        surface = self.mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        n = len(self.points)
        self.facet_pattern = _Pattern(np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel(), n)
        self.Bf = np.ones((3, 3))[None] * (self.areas / 9.0)[:, None, None]
        self._update_pinned()

    def _update_pinned(self):
        """Nodes with neither material nor film: they keep their temperature (identity rows). Recomputed whenever
        either the fractions or the film capacity change -- a node that had film and lost it must be pinned again in
        the same breath, or its row empties and the system goes singular (measured 2026-09-22)."""
        pinned = np.ones(len(self.points), dtype=bool)
        pinned[np.unique(self.tets[self.phi > 0.0])] = False
        pinned &= self.film_mass <= 0.0                          # a node carrying film still has a heat capacity
        if not np.array_equal(pinned, self.pinned):
            self._ml = None                                   # the operator's structure changed: fresh AMG hierarchy
        self.pinned = pinned

    def set_film_mass(self, mass):
        """Nodal mass [kg] of the melt film riding on the boundary (Step 3). The film is thermally thin -- q b / k_l =
        0.22 K across a 10 um film at 2 MW/m2, and b^2/alpha = 0.3 ms against a 0.5 s macro step -- so it is given no
        temperature of its own: its mass joins the boundary nodes' capacity with the *same* heat capacity the solid
        uses (tangent c_p,eff in the operator, secant [h(T) - h(T_old)]/(T - T_old) in the enthalpy rate), so the film's
        sensible heat is part of 1^T M dT by construction and `MeltingBody.film_energy` closes the balance exactly."""
        mass = np.asarray(mass, dtype=float)
        if mass.shape != (len(self.points),) or (mass < 0.0).any() or not np.isfinite(mass).all():
            raise ValueError("film mass must be one finite non-negative value per node")
        self.film_mass = mass
        self._update_pinned()

    def set_fractions(self, phi):
        """Element material fractions phi_e in [0, 1] (0 = dead): they scale the heat capacity of every element (the
        conductivity stays that of the full element while phi > 0, see operators()); the boundary follows the mesh's
        active set and nodes without material are pinned at their temperature. Call after the mesh's active set or
        the fractions change."""
        self.phi = np.asarray(phi, dtype=float)
        if self.phi.shape != (len(self.tets),) or (self.phi < 0.0).any() or (self.phi > 1.0).any():
            raise ValueError("fractions must be one value in [0, 1] per element")
        self._refresh_surface()

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

    def mass_matrix(self, c_nodal, c_film=None):
        """Capacity matrix for the nodal coefficient c [J/(m3 K)] (times phi_e per element). Lumped by default:
        diag(sum_e phi_e V_e/4 c_i), so that 1^T M dT = sum_e phi_e V_e/4 sum_i c_i dT_i is exactly the increment of
        the nodal enthalpy integral (the consistent form sum_k c_k Mk integrates the product of the interpolants of c
        and dT, which is not the increment of any energy functional and left a 3e-4 balance error, measured
        2026-09-20). `consistent_mass=True` keeps the consistent form (analytic tests only)."""
        # the film is liquid: its capacity is the liquid c_p (tangent) or the secant of the liquid enthalpy, never
        # the mixture's c_p,eff -- it has already paid its latent heat and does not pay it again on the ramp
        film = self.film_mass * (c_nodal if c_film is None else c_film) / self.material.rho
        if self.consistent_mass:
            M = self.pattern.assemble(np.einsum("ek,ekab->eab", self.phi[:, None] * c_nodal[self.tets], self.Mk))
            return M + sp.diags(film, format="csr") if self.film_mass.any() else M
        return sp.diags(np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * c_nodal[self.tets].ravel(),
                                    minlength=len(self.points)) + film, format="csr")

    def operators(self, T, T_old=None):
        """Stiffness K with k(T_e) at the element-mean temperature (active elements); the tangent mass M_tan with the nodal
        coefficient rho c_p,eff(T_i); and the enthalpy-rate vector E = M(rho c_sec) (T - T_old) with the nodal secant
        heat capacity c_sec,i = [h(T_i) - h(T_old,i)] / (T_i - T_old,i) (c_p,eff where a node has not moved). The
        enthalpy is nodal (spec Step 3 section 6, decided 2026-09-20: the element-mean enthalpy of Step 2 released the
        latent heat over a 4 K window of the element mean while the nodal temperatures span 30 K across a surface
        element, and its Newton iteration cycled): 1^T E is exactly the increment of rho int h dV with h interpolated
        linearly, whatever h(T) is. Without T_old, E is None and M_tan carries c_p,eff (the reference-operator convention)."""
        Te = T[self.tets].mean(axis=1)
        # conductivity is NOT scaled by phi_e: a partly consumed element is a thinner sliver of the same material,
        # which conducts better, not worse; scaling k with phi isolated the surface nodes of nearly consumed
        # elements and drove them to 5000 K (measured 2026-09-20). Dead elements (phi = 0) drop out.
        K = self.pattern.assemble(self.Ke * ((self.phi > 0.0) * self.material.k(Te))[:, None, None])
        c_tan, c_liq = self.material.cp_eff(T), self.material.cp(T)
        M = self.mass_matrix(self.material.rho * c_tan, self.material.rho * c_liq)
        if T_old is None:
            return K, M
        dT = T - T_old
        moved = np.abs(dT) > 1e-9
        c_sec = np.where(moved, (self.material.enthalpy(T) - self.material.enthalpy(T_old)) / np.where(moved, dT, 1.0), c_tan)
        c_sec_l = np.where(moved, (self.material.enthalpy_liquid(T) - self.material.enthalpy_liquid(T_old)) / np.where(moved, dT, 1.0), c_liq)
        # the film's mass enters both M and E through mass_matrix with its own liquid capacity, so 1^T E is exactly
        # the increment of (solid nodal enthalpy + film liquid enthalpy) -- which is what MeltingBody.energy sums
        return K, M, self.mass_matrix(self.material.rho * c_sec, self.material.rho * c_sec_l) @ dT

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def film_weight(self):
        """Share of each node's mass that is melt film, 0 to 1. The film is liquid: its enthalpy has no latent
        plateau, so the enthalpy-consistent Newton update -- which inverts the *mixture* h(T) -- must not be applied
        to a node the film owns, or the node is pushed onto a plateau it is not on and the iteration cycles (Newton
        stopped converging at 30 iterations once a patch was mostly film, measured 2026-09-22). The update is blended
        with the plain linear one by this weight, which is exact in both limits; it changes only the iteration, never
        the residual it converges to."""
        solid = np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * self.material.rho,
                            minlength=len(self.points))
        total = solid + self.film_mass
        return np.divide(self.film_mass, total, out=np.zeros_like(total), where=total > 0.0)

    def nodal_capacity(self):
        """Heat capacity carried by each node [J/K], material and film together (lumped, tangent c_p,eff). The body
        uses it to bound the melt loads it defers: energy booked onto a node that has melted away has nothing to
        heat, and pushing it in anyway moved drained surface nodes by 1e4-1e5 K in a step (measured 2026-09-22)."""
        c = self.material.cp_eff(self.T) * self.material.rho
        return np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * c[self.tets].ravel(),
                           minlength=len(self.points)) + self.film_mass * self.material.cp(self.T)

    def energy(self):
        """Stored enthalpy above T_REF: rho int h dV with the nodal h(T_i) interpolated linearly, i.e.
        sum_e phi_e rho V_e mean_i h(T_i); equals 1^T M T for constant c_p (consistent or lumped mass)."""
        return float(self.element_energies().sum())

    def element_energies(self, T=None):
        """phi_e rho V_e mean_i h(T_i) per element [J]."""
        h = self.material.enthalpy(self.T if T is None else T)
        return self.phi * self.material.rho * self.vol * h[self.tets].mean(axis=1)

    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
        """One backward-Euler step with the per-facet convective flux q_conv [W/m2] and an optional nodal load vector
        [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped)."""
        T_old = self.T
        # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        F_extra = np.zeros(len(self.points)) if nodal_load is None else np.asarray(nodal_load, dtype=float)
        Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
        T_new, iteration = T_k, 0
        converged = False
        last_relative_change, previous_change, damping, decreases = None, None, 1.0, 0
        mat = self.material
        film_w = self.film_weight() if self.film_mass.any() else np.zeros(len(self.points))
        for iteration in range(1, self.max_iterations + 1):
            K, M, E = self.operators(T_k, T_old)
            Tf = self.facet_temperature(T_k)
            B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            # Newton on R(T) = E(T)/dt + K T - F_conv + F_rad(T) - F_extra with the tangent M/dt + K + B:
            # A T_new = A T_k - R(T_k); with a smooth c_p (c_tan = c_sec) this is the Step 2 iteration unchanged
            A = (M / dt + K + B).tocsr()
            b = M @ T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
            if self.pinned.any():                              # material-free nodes keep their temperature (identity rows)
                A = A + sp.diags(self.pinned.astype(float), format="csr")
                b[self.pinned] = T_old[self.pinned]
            T_new = self._solve_with_dirichlet(A, b, T_k, dirichlet)
            if mat.melts:
                # enthalpy-consistent update: the linearised step is an enthalpy increment c_p,eff(T_k) (T_new - T_k)
                # per node; inverting the true h(T) puts a node that would overshoot the melting range where the
                # latent heat actually leaves it (identity where h is linear). Without it nodes jump across the +-2 K
                # ramp of a single-temperature material and the iteration cycles (measured 2026-09-20). The enthalpy
                # inverted is the node's own mixture of material and film (see film_weight), not the material's.
                free, w = ~self.pinned, film_w[~self.pinned]
                T_new[free] = mat.temperature_from_enthalpy_mixed(
                    mat.enthalpy_mixed(T_k[free], w) + mat.cp_mixed(T_k[free], w) * (T_new[free] - T_k[free]), w)
            last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            # damping when the change grows (a fallback for cycling iterates), released again after two decreases
            if previous_change is not None and last_relative_change > previous_change:
                damping, decreases = max(0.25, 0.5 * damping), 0
            elif damping < 1.0:
                decreases += 1
                if decreases >= 2:
                    damping, decreases = 1.0, 0
            previous_change = last_relative_change
            T_k = T_k + damping * (T_new - T_k) if damping < 1.0 else T_new
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev, self.T = T_old, T_new
        self.last_damping = damping
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)

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
        for attempt in (0, 1):
            if fresh or attempt or self._ml is None or self._solves % self.amg_rebuild_every == 0:
                import pyamg
                self._ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
            self._solves += 1
            counter = []
            x, info = spla.cg(A, b, x0=x0, rtol=self.cg_tol, maxiter=500, M=self._ml.aspreconditioner(cycle="V"),
                              callback=lambda _: counter.append(1))
            self.last_cg_iterations = len(counter)
            if info == 0:
                return x
        # Melting drains interior elements to phi ~ 1e-3 while their conduction stays unscaled (spec section 9), which
        # raises the condition number by 1/phi; a hierarchy built before the drain can stall on the emptied region.
        # A fresh hierarchy clears it in every case measured except the step after the heating is switched off, where
        # the operator changes by orders of magnitude in one step -- that one system is solved directly.
        self.direct_fallbacks += 1
        return spla.spsolve(A.tocsc(), b)

    def reference_operators(self, T):
        """K (element-mean k) and M (nodal c_p) assembled by scikit-fem with the same coefficients (conformance test only)."""
        from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP1, MeshTet, asm
        from skfem.helpers import dot, grad
        basis = Basis(MeshTet(self.points.T.copy(), self.tets.T.copy()), ElementTetP1(), intorder=3)   # exact for phi_i phi_j c (cubic)
        basis0 = basis.with_element(ElementTetP0())
        Te = T[self.tets].mean(axis=1)

        @BilinearForm
        def stiffness(u, v, w):
            return w.k * dot(grad(u), grad(v))

        @BilinearForm
        def mass(u, v, w):
            return w.c * u * v

        K = asm(stiffness, basis, k=basis0.interpolate(self.material.k(Te)))
        M = asm(mass, basis, c=basis.interpolate(self.material.rho * self.material.cp_eff(T)))    # nodal (P1) coefficient
        return K, M
