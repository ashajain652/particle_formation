"""The thermal solve on cut cells (design section 7, stage 2 of section 20.11): P1 temperature on the adaptive
tetrahedral mesh, integrated over the metal {phi < 0} only, with the boundary conditions applied on the zero
contour's patches. Level sets are negative inside the metal, patch normals point out of the metal, heat flux into the
metal is positive, SI units throughout (K, W, J, m).

Capacity. Each active cell's metal part is split into sub-tetrahedra (`cutcell.inside_subtets`); its consistent mass
(the integral of N_a N_b over the metal part, unit coefficient) is integrated exactly on them with a four-point rule
exact for quadratics. The cell's temperature for its coefficients is T* = sum_a w_a T_a, w_a the mean of the shape
function a over the metal part, so that 1^T M_e dT = rho c V_in dT* exactly; with the secant capacity
c = [h(T*) - h(T*_old)] / (T* - T*_old) the books close to the Newton tolerance whatever the enthalpy h(T) is (latent
heat included), and the stored enthalpy is sum_e rho V_in,e h(T*_e).

Conduction. P1 gradients are constant per cell, so the stiffness is k(T*) grad N_a . grad N_b V_in.

Boundary. Per patch p (area A_p, weights w_pa = the shape functions at its centroid, exact for a flat patch):
    q_p = q0_p + q1_p T_p - rad_p eps sigma (T_p^4 - T_amb^4),     T_p = sum_a w_pa T_a,
heat into the metal in W/m2; q1 <= 0 is the skins' base law (implicit), rad_p is 0 or 1. The boundary flux is
linearised by Newton; k and the secant capacity are lagged at the iterate (a Picard iteration), so the step is damped
when the change stops falling (a kink of the effective capacity).

The material is duck-typed: `rho` (kg/m3) and vectorised `k(T)` (W/m/K), `cp_eff(T)` (J/kg/K), `enthalpy(T)` (J/kg).
The unknowns are the active nodes (those of a cell with metal in it) only; the other nodes keep their values.

Slivers. A ghost penalty (Burman 2010) on every interior face between two active cells of which at least one is cut,
gamma k h |F| [dN/dn][dN/dn]^T with k the larger of the two cells' conductivities, h the larger cell's size and |F| the
face area, keeps the system conditioned when a cell holds almost no metal. The jump vector [dN/dn in cell 1, -dN/dn in
cell 2] annihilates any function that is linear across the face, so constants are untouched. `ghost` is gamma."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from reentry_model.levelset.cutcell import inside_subtets, subtet_volumes
from reentry_model.levelset.cutcell import patches as cut_patches
from reentry_model.levelset.transfer import _geom, barycentric

SIGMA_SB = 5.670374419e-8          # Stefan-Boltzmann constant, W/m2/K4

# Four-point rule on a tetrahedron, exact for quadratics (barycentric coordinates of the points, weights 1/4 each).
_A, _B = 0.5854101966249685, 0.1381966011250105
_QP = np.array([[_A, _B, _B, _B], [_B, _A, _B, _B], [_B, _B, _A, _B], [_B, _B, _B, _A]])
_FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])     # face i of a tetrahedron is opposite vertex i


class _Pattern:
    """A fixed CSR sparsity pattern for (row, col) pairs; `assemble(values)` sums duplicates with bincount."""

    def __init__(self, rows, cols, n):
        key = rows.astype(np.int64) * n + cols.astype(np.int64)
        unique, self._map = np.unique(key, return_inverse=True)
        self._indptr = np.searchsorted(unique // n, np.arange(n + 1))
        self._indices = (unique % n).astype(np.int32)
        self.n, self.nnz = n, unique.size

    def assemble(self, values):
        data = np.bincount(self._map, weights=np.asarray(values).ravel(), minlength=self.nnz)
        return sp.csr_matrix((data, self._indices, self._indptr), shape=(self.n, self.n))


class CutThermal:
    """Heat conduction in the metal {phi < 0} of a level-set mesh; see the module docstring."""

    def __init__(self, material, emissivity=0.0, ghost=0.01, lumped=False, newton_tol=1e-8, max_iterations=100,
                 linear_solver="amg", cg_tol=1e-12):
        self.mat, self.eps, self.ghost, self.lumped = material, float(emissivity), float(ghost), lumped
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.linear_solver, self.cg_tol = linear_solver, cg_tol
        self.T = None

    # --- geometry --------------------------------------------------------------------------------------------
    def setup(self, mesh, phi):
        """Build the active cells and nodes, the metal volumes Vin (m3), weights w, masses Me, stiffnesses Ke and the
        patches of the cut cells for the level set phi (m, negative inside) on the mesh."""
        self.mesh, self.phi = mesh, np.asarray(phi, float)
        n = len(mesh.points)
        tets = mesh.tets
        nneg = (self.phi[tets] < 0).sum(axis=1)
        self.cells = np.flatnonzero(nneg > 0)                    # active cells
        full = nneg[self.cells] == 4
        _, Jinv, _ = _geom(mesh)
        V = mesh.volumes[self.cells]
        grad_ref = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        self.G = np.einsum("eji,aj->eai", Jinv[self.cells], grad_ref)          # (E, 4 nodes, 3) grad N_a
        Vin = np.where(full, V, 0.0)
        w = np.full((len(self.cells), 4), 0.25)
        Me = (np.ones((4, 4)) + np.eye(4))[None] * (V / 20.0)[:, None, None]
        cut = np.flatnonzero(~full)
        if len(cut):
            S, own = inside_subtets(mesh, self.phi, self.cells[cut])
            vs = subtet_volumes(S)
            Me[cut] = 0.0
            Iw = np.zeros((len(cut), 4))
            Mc = np.zeros((len(cut), 4, 4))
            for q in _QP:
                xq = np.einsum("j,sjk->sk", q, S)
                Nq = barycentric(mesh, self.cells[cut][own], xq)            # (S, 4)
                np.add.at(Iw, own, (vs / 4)[:, None] * Nq)
                np.add.at(Mc, own, (vs / 4)[:, None, None] * Nq[:, :, None] * Nq[:, None, :])
            vin = np.bincount(own, vs, minlength=len(cut))
            Vin[cut] = vin
            w[cut] = Iw / np.maximum(vin, 1e-300)[:, None]
            Me[cut] = Mc
        self.Vin, self.w, self.Me = Vin, w, Me
        self.Ke = np.einsum("eai,ebi->eab", self.G, self.G) * Vin[:, None, None]
        T = tets[self.cells]
        self.active = np.zeros(n, bool)
        self.active[T.ravel()] = True
        self.pattern = _Pattern(np.repeat(T, 4, axis=1).ravel(), np.tile(T, (1, 4)).ravel(), n)
        self._ghost_faces(~full)
        self.P = cut_patches(mesh, self.phi, self.cells[~full])
        Pc = self.P.cell
        self.Pw = barycentric(mesh, Pc, self.P.centroid)
        self.Pnodes = tets[Pc]
        self.Ppattern = _Pattern(np.repeat(self.Pnodes, 4, axis=1).ravel(), np.tile(self.Pnodes, (1, 4)).ravel(), n)

    def _ghost_faces(self, is_cut):
        """Interior faces between two active cells, at least one of them cut: node slots, normal-gradient jump
        vectors, scale h |F| and the two cells of each face. Faces are matched by sorting each cell's vertex triples."""
        T = self.mesh.tets[self.cells]
        E = len(T)
        keys = np.sort(T[:, _FACES], axis=2).reshape(-1, 3)                   # (4E, 3)
        owner = np.repeat(np.arange(E), 4)
        order = np.lexsort(keys.T[::-1])
        ks = keys[order]
        same = np.all(ks[1:] == ks[:-1], axis=1)
        i1, i2 = order[:-1][same], order[1:][same]
        e1, e2 = owner[i1], owner[i2]
        keep = is_cut[e1] | is_cut[e2]
        e1, e2, i1 = e1[keep], e2[keep], i1[keep]
        f = keys[i1]
        p = self.mesh.points
        nrm = np.cross(p[f[:, 1]] - p[f[:, 0]], p[f[:, 2]] - p[f[:, 0]])
        area = 0.5 * np.linalg.norm(nrm, axis=1)
        nrm = nrm / np.maximum(2 * area, 1e-300)[:, None]
        g1 = np.einsum("fai,fi->fa", self.G[e1], nrm)                       # dN_a/dn in cell 1, its 4 nodes
        g2 = np.einsum("fai,fi->fa", self.G[e2], nrm)
        nodes = np.concatenate([T[e1], T[e2]], axis=1)                        # (F, 8), shared nodes repeat
        jump = np.concatenate([g1, -g2], axis=1)
        h = self.mesh.size[self.cells][np.stack([e1, e2], 1)].max(axis=1)
        self.gf_nodes, self.gf_jump, self.gf_scale, self.gf_cells = nodes, jump, h * area, np.stack([e1, e2], 1)
        self.gpattern = _Pattern(np.repeat(nodes, 8, axis=1).ravel(), np.tile(nodes, (1, 8)).ravel(),
                                 len(self.mesh.points))

    # --- state --------------------------------------------------------------------------------------------------
    def set_temperature(self, T):
        n = len(self.mesh.points)
        self.T = np.full(n, float(T)) if np.ndim(T) == 0 else np.array(T, float)

    def cell_T(self, T=None):
        T = self.T if T is None else T
        return np.einsum("ea,ea->e", self.w, T[self.mesh.tets[self.cells]])

    def energy(self, T=None):
        """Stored enthalpy above the material's reference [J]."""
        return float((self.mat.rho * self.Vin * self.mat.enthalpy(self.cell_T(T))).sum())

    def metal_volume(self):
        return float(self.Vin.sum())

    def patch_T(self, T=None):
        T = self.T if T is None else T
        return np.einsum("pa,pa->p", self.Pw, T[self.Pnodes])

    # --- step ---------------------------------------------------------------------------------------------------
    def operators(self, T, T_old):
        """Stiffness K (W/K) at the cell temperatures of T and the secant-capacity mass M (J/K) between T_old and T,
        both (n, n) sparse over all nodes. K includes the ghost-penalty term."""
        Te, Te_old = self.cell_T(T), self.cell_T(T_old)
        k = self.mat.k(Te)
        c = self.mat.cp_eff(Te)
        dT = Te - Te_old
        moved = np.abs(dT) > 1e-9
        c = np.where(moved, (self.mat.enthalpy(Te) - self.mat.enthalpy(Te_old)) / np.where(moved, dT, 1.0), c)
        K = self.pattern.assemble(self.Ke * k[:, None, None])
        if self.lumped:
            rowsum = self.Me.sum(axis=2) * (self.mat.rho * c)[:, None]
            M = sp.csr_matrix((np.bincount(self.mesh.tets[self.cells].ravel(), rowsum.ravel(),
                                            minlength=len(self.T)), (np.arange(len(self.T)), np.arange(len(self.T)))),
                              shape=(len(self.T), len(self.T)))
        else:
            M = self.pattern.assemble(self.Me * (self.mat.rho * c)[:, None, None])
        kf = k[self.gf_cells].max(axis=1)
        J = np.einsum("fa,fb->fab", self.gf_jump, self.gf_jump) * (self.ghost * kf * self.gf_scale)[:, None, None]
        return K + self.gpattern.assemble(J), M

    def step(self, dt, q0, q1=None, rad=None, T_amb=0.0):
        """One backward-Euler step; q0, q1, rad per patch (heat into the metal). Returns heat in over the step [J]."""
        T_old = self.T
        q0 = np.broadcast_to(np.asarray(q0, float), (len(self.P.area),))
        q1 = np.zeros_like(q0) if q1 is None else np.broadcast_to(np.asarray(q1, float), q0.shape)
        rad = np.zeros_like(q0) if rad is None else np.broadcast_to(np.asarray(rad, float), q0.shape)
        idx = np.flatnonzero(self.active)
        Tk = T_old.copy()
        omega, last = 1.0, None
        for it in range(1, self.max_iterations + 1):
            K, M = self.operators(Tk, T_old)
            Tp = self.patch_T(Tk)
            qp = q0 + q1 * Tp - rad * self.eps * SIGMA_SB * (Tp ** 4 - T_amb ** 4)
            dq = q1 - rad * 4 * self.eps * SIGMA_SB * Tp ** 3                 # d q_p / d T_p
            F = np.bincount(self.Pnodes.ravel(), (self.P.area[:, None] * self.Pw * qp[:, None]).ravel(),
                            minlength=len(Tk))
            B = self.Ppattern.assemble(-(self.P.area * dq)[:, None, None] * self.Pw[:, :, None] * self.Pw[:, None, :])
            A = (M / dt + K + B).tocsr()[idx][:, idx]
            b = (M @ T_old / dt + F + B @ Tk)[idx]
            x = self._solve(A, b, Tk[idx])
            change = np.linalg.norm(x - Tk[idx]) / np.linalg.norm(x)
            if change <= self.newton_tol:
                Tk = Tk.copy(); Tk[idx] = x
                break
            if last is not None and change > 0.9 * last:                    # oscillating across a kink of c_eff
                omega = max(0.5 * omega, 1.0 / 16)
            last = change
            Tk = Tk.copy(); Tk[idx] = Tk[idx] + omega * (x - Tk[idx])
        else:
            raise RuntimeError("Newton did not converge ({:.2e})".format(change))
        self.T = Tk
        Tp = self.patch_T()
        qp = q0 + q1 * Tp - rad * self.eps * SIGMA_SB * (Tp ** 4 - T_amb ** 4)
        self.iterations = it
        return float((self.P.area * qp).sum() * dt)

    def _solve(self, A, b, x0):
        if self.linear_solver == "direct":
            return spla.spsolve(A.tocsc(), b)
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
        x, info = spla.cg(A, b, x0=x0, rtol=self.cg_tol, maxiter=1000, M=ml.aspreconditioner(cycle="V"))
        if info != 0:
            raise RuntimeError("CG did not converge (info {})".format(info))
        return x
