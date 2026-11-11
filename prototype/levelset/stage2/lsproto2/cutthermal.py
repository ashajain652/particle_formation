"""The thermal solve on cut cells (design §7, stage 2 of §20.11): P1 temperature on the adaptive tetrahedral mesh,
integrated over the metal {phi < 0} only, with the boundary conditions applied on the zero contour's patches.

Capacity. Each active cell's metal part is split into sub-tetrahedra (one, or three for a wedge); its consistent mass
is integrated exactly on them. The cell's temperature for its coefficients is T* = sum_a w_a T_a, w_a the mean of the
shape function a over the metal part, so that 1^T M_e dT = rho c V_in dT* exactly; with the secant capacity
c = [h(T*) - h(T*_old)] / (T* - T*_old) the books close to the Newton tolerance whatever h(T) is (latent heat included),
and the stored enthalpy is sum_e rho V_in,e h(T*_e).

Conduction. P1 gradients are constant per cell, so the stiffness is k(T*) grad N_a . grad N_b V_in.

Slivers. A ghost penalty (Burman 2010) on every interior face of an active cut cell: gamma k h |F| [dN/dn][dN/dn]^T,
the jump of the normal gradient across the face. It vanishes on constants, so the energy books are untouched.

Boundary. Per patch p (area A_p, weights w_pa = the shape functions at its centroid, exact for a flat patch):
    q_p = q0_p + q1_p T_p - rad_p eps sigma (T_p^4 - T_amb^4),     T_p = sum_a w_pa T_a,
heat into the metal; q1 <= 0 is the skins' base law (implicit), rad_p 0 or 1. Newton on the radiation term.

Moving surface. `move(phi_new)` keeps the nodal temperatures, gives newly active nodes values from their active
neighbours, recomputes the geometry, and returns the enthalpy the move took out (positive) or brought in.
`remesh(new_mesh, phi_new)` projects T conservatively and closes the stored enthalpy with one uniform offset,
reporting the residual before it."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from reentry_model.levelset.cutcell import patches as cut_patches
from reentry_model.levelset.transfer import _geom, barycentric, project

SIGMA_SB = 5.670374419e-8

_A, _B = 0.5854101966249685, 0.1381966011250105
_QP = np.array([[_A, _B, _B, _B], [_B, _A, _B, _B], [_B, _B, _A, _B], [_B, _B, _B, _A]])
_FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])     # face i is opposite vertex i


def _tvol(P):
    return np.abs(np.einsum("ij,ij->i", P[:, 1] - P[:, 0], np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6.0


def inside_subtets(mesh, phi, cells):
    """The metal part {phi < 0} of each given cut cell as sub-tetrahedra: (S, 4, 3) vertices and (S,) index into cells."""
    T = mesh.tets[cells]
    f, P = phi[T], mesh.points[T]
    neg = f < 0
    nneg = neg.sum(axis=1)
    cut = lambda x, y, fx, fy: x + (fx / (fx - fy))[:, None] * (y - x)
    out, owner = [], []
    m = np.flatnonzero(nneg == 1)
    if len(m):
        lone = np.argmax(neg[m], axis=1)
        others = np.array([[j for j in range(4) if j != l] for l in range(4)])[lone]
        r = np.arange(len(m))
        a, fa = P[m][r, lone], f[m][r, lone]
        pts = [cut(a, P[m][r, others[:, i]], fa, f[m][r, others[:, i]]) for i in range(3)]
        out.append(np.stack([a] + pts, axis=1)); owner.append(m)
    m = np.flatnonzero(nneg == 3)
    if len(m):
        lone = np.argmax(~neg[m], axis=1)                    # the single outside vertex
        others = np.array([[j for j in range(4) if j != l] for l in range(4)])[lone]
        r = np.arange(len(m))
        d, fd = P[m][r, lone], f[m][r, lone]
        a, b, c = (P[m][r, others[:, i]] for i in range(3))
        fa, fb, fc = (f[m][r, others[:, i]] for i in range(3))
        pa, pb, pc = cut(a, d, fa, fd), cut(b, d, fb, fd), cut(c, d, fc, fd)
        for tet in ((a, b, c, pa), (b, c, pa, pb), (c, pa, pb, pc)):
            out.append(np.stack(tet, axis=1)); owner.append(m)
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")
        rr = np.arange(len(m))[:, None]
        Po, fo = P[m][rr, order], f[m][rr, order]
        a, b, c, d = (Po[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pad, pbc, pbd = cut(a, c, fa, fc), cut(a, d, fa, fd), cut(b, c, fb, fc), cut(b, d, fb, fd)
        for tet in ((a, pac, pad, b), (pac, pad, b, pbc), (pad, b, pbc, pbd)):
            out.append(np.stack(tet, axis=1)); owner.append(m)
    if not out:
        return np.zeros((0, 4, 3)), np.zeros(0, int)
    return np.concatenate(out), np.concatenate(owner)


class _Pattern:
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
    def __init__(self, material, emissivity=0.0, ghost=0.01, lumped=False, newton_tol=1e-8, max_iterations=100,
                 linear_solver="amg", cg_tol=1e-12):
        self.mat, self.eps, self.ghost, self.lumped = material, float(emissivity), float(ghost), lumped
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.linear_solver, self.cg_tol = linear_solver, cg_tol
        self.T = None

    # --- geometry ---------------------------------------------------------------------------------------------
    def setup(self, mesh, phi):
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
            vs = _tvol(S)
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
        """Interior faces between two active cells, at least one of them cut: the normal-gradient jump vectors."""
        T = self.mesh.tets[self.cells]
        E = len(T)
        keys = np.sort(T[:, _FACES], axis=2).reshape(-1, 3)                   # (4E, 3)
        owner = np.repeat(np.arange(E), 4)
        opp = np.tile(np.arange(4), E)
        order = np.lexsort(keys.T[::-1])
        ks = keys[order]
        same = np.all(ks[1:] == ks[:-1], axis=1)
        i1, i2 = order[:-1][same], order[1:][same]
        e1, e2 = owner[i1], owner[i2]
        keep = is_cut[e1] | is_cut[e2]
        e1, e2, i1, i2 = e1[keep], e2[keep], i1[keep], i2[keep]
        f = keys[i1]
        p = self.mesh.points
        nrm = np.cross(p[f[:, 1]] - p[f[:, 0]], p[f[:, 2]] - p[f[:, 0]])
        area = 0.5 * np.linalg.norm(nrm, axis=1)
        nrm = nrm / np.maximum(2 * area, 1e-300)[:, None]
        g1 = np.einsum("fai,fi->fa", self.G[e1], nrm)                       # dN_a/dn in cell 1, its 4 nodes
        g2 = np.einsum("fai,fi->fa", self.G[e2], nrm)
        nodes = np.concatenate([T[e1], T[e2]], axis=1)                        # (F, 8) with repeats; the jump vector
        jump = np.concatenate([g1, -g2], axis=1)
        h = self.mesh.size[self.cells][np.stack([e1, e2], 1)].max(axis=1)
        self.gf_nodes, self.gf_jump, self.gf_scale, self.gf_cells = nodes, jump, h * area, np.stack([e1, e2], 1)
        n = len(self.mesh.points)
        self.gpattern = _Pattern(np.repeat(nodes, 8, axis=1).ravel(), np.tile(nodes, (1, 8)).ravel(), n)

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
        K = K + self.gpattern.assemble(J)
        return K, M

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

    # --- moving surface -----------------------------------------------------------------------------------------
    def _extend(self, donors, targets, rings=50):
        """Give the target nodes the mean of their valued neighbours (donors first), ring by ring."""
        T = self.T.copy()
        tets = self.mesh.tets
        nc = self.mesh.node_cells()
        valid = donors.copy()
        todo = np.flatnonzero(targets & ~valid)
        for _ in range(rings):
            if not len(todo):
                break
            sub = nc[todo]
            rows = np.repeat(np.arange(len(todo)), np.diff(sub.indptr))
            nb = tets[sub.indices]
            v = valid[nb]
            s = np.bincount(rows, (T[nb] * v).sum(axis=1), minlength=len(todo))
            c = np.bincount(rows, v.sum(axis=1), minlength=len(todo))
            ok = c > 0
            T[todo[ok]] = s[ok] / c[ok]
            valid[todo[ok]] = True
            todo = todo[~ok]
        return T

    def move(self, phi_new):
        """Move the surface on the same mesh; returns the enthalpy the move removed from the metal [J]."""
        E0 = self.energy()
        donors = self.active.copy()
        self.setup(self.mesh, phi_new)
        self.T = self._extend(donors, self.active)
        return E0 - self.energy()

    def remesh(self, new_mesh, phi_new):
        """Carry T to a new mesh conservatively and close the stored enthalpy with a uniform offset.
        Returns the relative residual before the offset."""
        E0 = self.energy()
        old_mesh = self.mesh
        T_full = self._extend(self.active.copy(), ~self.active, rings=3)     # a band outside, for the projection
        T_new = project(old_mesh, new_mesh, T_full)
        self.setup(new_mesh, phi_new)
        self.T = T_new
        r0 = (self.energy() - E0) / max(abs(E0), 1e-300)
        a, b = -50.0, 50.0                                                   # bisection on the offset (monotone)
        for _ in range(200):
            m = 0.5 * (a + b)
            if self.energy(T_new + m) > E0:
                b = m
            else:
                a = m
            if b - a < 1e-12:
                break
        self.T = T_new + 0.5 * (a + b)
        return r0
