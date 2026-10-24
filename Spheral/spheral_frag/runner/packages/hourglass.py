"""Hourglass: a small stiffness against the non-affine part of the displacement (Spheral M2 Task 10b).

**Why.** SolidSPH's tensile instability (Swegle et al. 1995): under tension a particle row's alternating ("zigzag")
displacement is a mode of negative stiffness when the neighbours sit beyond the kernel's inflection point -- with
Spheral's WendlandC4 (extent 1) at nPerh 2.01 every neighbour does (inflection q 0.235 against the nearest
neighbour's 0.45-0.50) -- and SPH's gradient estimate does not see the alternating mode at all (a zero-energy
mode). On the MVP flight's replay the core's 40-50 MPa axial tension grows it exponentially (Task 10: RZ 1.1 mm
x3.2 per interval at 10 crossings, the core's axial stress lost). A growing exponential mode is only slowed, never
stopped, by velocity damping (with x'' = g^2 x - 2 b x' the growth rate is sqrt(b^2 + g^2) - b > 0); a restoring
stiffness larger than g^2 stops it.

**The stiffness.** Every particle i's displacement from its reference position u_i = x_i - x_ref,i (D1: the
replay's bodies are Lagrangian about x_ref) is compared with the linear (affine) least-squares fit of its
neighbours' displacements at its reference position: d_i = u_i - sum_j c_ij u_j, with c_ij the fit's weights over
the neighbours j != i whose reference positions lie within `RADIUS` spacings of i's (the lattice's nearest and
diagonal neighbours: 8 in RZ, 18 in 3D, fewer at the surface), so that d = A u with A = I - C vanishes for every
affine displacement (translation, rotation, uniform strain) exactly. The energy

    E_hg = 1/2 sum_i k_i |d_i|^2,   k_i = kappa G(T_i) V_i / dx^2

(G the table's shear modulus at the particle's temperature, zero in the melt; V = m / rho the physical volume)
gives the forces F = -A^T (k d): they sum to zero (A 1 = 0: momentum conserved, also axially in RZ), vanish on
affine fields, and on a smooth field u of length scale L act as a relative stiffening of order kappa (dx/L)^2 in the
interior and kappa at the surface. On the zigzag the restoring rate^2 is of order kappa G / (rho dx^2), against the
instability's growth rate^2, of order (sigma/E) c^2/dx^2 -- both scale as 1/dx^2, so the kappa needed does not depend
on the spacing but grows with the tension sigma/E. Measured on check 13's core at 40-50 MPa (sigma/E ~ 6e-4; Task
10b, analysis/spheral_m2_task10b.py): kappa 1e-3 and 3e-3 slow the pairing but do not stop it, 1e-2 removes it at 10,
20 and 40 crossings -- the threshold lies between 3e-3 and 1e-2 at this tension, so the default 0.01 (run.py
--hourglass) has a margin of at most about 3, which a body in more tension (M3) will have to recheck. The stiffness alone at 1e-2 left the two rows next to
the RZ axis oscillating axially at 40 crossings (internal KE/E 9.9e-3): a stiff, undamped non-affine mode, which the
viscous part F_v = -A^T (b A v), b_i = BETA sqrt(k_i m_i), removes (9.2e-6); it too vanishes on affine velocities.
The acceleration is F / m with Spheral's (scaled) masses: the stiffness is a physical one, like
the stresses, invariant under mass scaling.

**Limits.** The reference is x_ref (D1: the replay's bodies are Lagrangian about their build positions): accumulated
non-affine plastic flow on the scale of a spacing (M3's tearing, M4's slurry) would be resisted as well, and a moving
reference is then needed. The fit is over reference positions, so a particle's neighbours are its lattice neighbours
whatever the current positions.

**Built per interval** (as every package here): A, G and V at the interval's start, over the whole body (each rank
holds the global A, from the allgathered reference positions in rank order); each evaluation allgathers the current
positions and velocities, forms A u and A v globally and F = -A^T (k A u + b A v) for its own rows. Cost measured in
Task 10b: 1-9 % of a step (RZ 1.1 mm on 4 processes +3-7 %, 3D 2.2 mm on 18 +7 %, 3D 3.0 mm +4-9 %).
`last_energy` is the last evaluation's E_hg (J; recorded over the elastic energy as intervals.csv's
`hourglass_energy_ratio`: 4-15e-4 on the MVP flight's replay)."""
from __future__ import annotations

import numpy as np

from spheral_frag import material
from spheral_frag.runner.packages import LocalCount, add_acceleration, comm, get_vector, state_velocity

RADIUS = 1.5            # spacings: the neighbours of the fit (lattice nearest and diagonal)
BETA = 1.0              # the viscous part's ratio: b_i = BETA sqrt(k_i m_i) (critical-like on each particle's own mode)

_CLASSES = {}


def _fit_rows(x, dx, radius, rows):
    """(i, j, c, good) of the affine fit for the particles `rows`, neighbours within radius dx (module docstring)."""
    from scipy.spatial import cKDTree
    n, dim = x.shape
    tree = cKDTree(x)
    nb = tree.query_ball_point(x[rows], radius*dx*(1.0 + 1e-9))
    i = np.concatenate([np.full(len(l), r) for r, l in zip(rows, nb)]).astype(np.int64) if len(rows) else \
        np.zeros(0, np.int64)
    j = np.concatenate([np.asarray(l, dtype=np.int64) for l in nb]) if len(rows) else np.zeros(0, np.int64)
    keep = i != j
    i, j = i[keep], j[keep]
    p = np.concatenate([np.ones((len(i), 1)), (x[j] - x[i])/dx], axis=1)        # (nnz, dim + 1)
    M = np.zeros((n, dim + 1, dim + 1))
    np.add.at(M, i, p[:, :, None]*p[:, None, :])
    M = M[rows]
    cond = np.linalg.cond(M) if len(rows) else np.zeros(0)
    good = np.isfinite(cond) & (cond < 1e6)
    e0 = np.zeros(dim + 1)
    e0[0] = 1.0
    q = np.zeros((n, dim + 1))
    rg = np.asarray(rows)[good]
    if len(rg):
        q[rg] = np.linalg.solve(M[good], np.broadcast_to(e0, (len(rg), dim + 1))[..., None])[..., 0]
    c = np.einsum("nk,nk->n", q[i], p)
    return i, j, c, dict(zip(np.asarray(rows).tolist(), good.tolist()))


def fit_operator(x_ref, dx, radius=RADIUS):
    """The sparse matrix A = I - C (scipy CSR, n x n) of the affine fit: (A u)_i = u_i - sum_j c_ij u_j, with c_ij
    the least-squares weights (equal weights) that reproduce u_i from its neighbours j != i within radius dx of
    x_ref,i for every affine u. A particle whose neighbours there do not determine an affine fit (a corner of the
    lattice's staircase) takes its neighbours within 1.5 times the radius, then twice it. Returns (A, the number
    of particles that needed the wider neighbourhood)."""
    from scipy.sparse import csr_matrix, identity
    x = np.asarray(x_ref, dtype=np.float64)
    n = len(x)
    rows = np.arange(n)
    I, J, Cv = [], [], []
    widened = 0
    for f in (1.0, 1.5, 2.0):
        i, j, c, good = _fit_rows(x, dx, radius*f, rows)
        ok = np.array([good[r] for r in i.tolist()], dtype=bool) if len(i) else np.zeros(0, bool)
        I.append(i[ok]); J.append(j[ok]); Cv.append(c[ok])
        rows = np.array([r for r, g in good.items() if not g], dtype=np.int64)
        if f > 1.0:
            widened += len(good) - len(rows)
        if not len(rows):
            break
    if len(rows):
        raise RuntimeError(f"Hourglass: {len(rows)} particles have no affine fit within {2*radius} spacings")
    C = csr_matrix((np.concatenate(Cv), (np.concatenate(I), np.concatenate(J))), shape=(n, n))
    return (identity(n, format="csr") - C).tocsr(), widened


def hourglass_class(S):
    if S in _CLASSES:
        return _CLASSES[S]

    class Hourglass(S.Physics):
        """F = -A^T (k A u) (module docstring). Hourglass(body, kappa, table): built from the body's state at the
        interval's start."""

        def __init__(self, body, kappa, table, beta=BETA):
            S.Physics.__init__(self)
            from spheral_frag.runner.eos import get_internal
            self.body = body
            self.dim = body.dim
            self.kappa = float(kappa)
            self.count = LocalCount(body.nodes)
            nodes = body.nodes
            x_ref = get_vector(body.x_ref, self.dim)
            m = np.array(get_internal(nodes.mass()))
            rho = np.array(get_internal(nodes.massDensity()))
            u = np.array(get_internal(nodes.specificThermalEnergy()))
            G = material.thermal_shear_modulus(table, table.temperature(u))
            k_local = self.kappa*G*(m/rho)/body.dx**2              # m/rho: scaled/scaled = physical volume
            b_local = float(beta)*np.sqrt(k_local*m)                # m: Spheral's (scaled) mass
            self.beta = float(beta)
            parts = comm().allgather((x_ref, k_local, b_local))
            self.b_all = np.concatenate([p[2] for p in parts])
            self.counts = [len(p[0]) for p in parts]
            self.offset = int(sum(self.counts[:comm().rank]))
            xr_all = np.concatenate([p[0] for p in parts], axis=0)
            self.k_all = np.concatenate([p[1] for p in parts])
            self.x_ref_all = xr_all
            self.A, self.n_widened = fit_operator(xr_all, body.dx)
            n_loc = len(x_ref)
            self.AT_local = self.A.T.tocsr()[self.offset:self.offset + n_loc]
            self.inv_m = 1.0/m
            self.last_energy = 0.0
            self.nnz_per_row = self.A.nnz/max(len(xr_all), 1)

        def acceleration(self, x_local, v_local):
            parts = comm().allgather((x_local, v_local))
            u = np.concatenate([p[0] for p in parts], axis=0) - self.x_ref_all
            v = np.concatenate([p[1] for p in parts], axis=0)
            d = self.A @ u
            dv = self.A @ v
            self.last_energy = 0.5*float(np.sum(self.k_all*np.einsum("ij,ij->i", d, d)))
            F = -(self.AT_local @ (self.k_all[:, None]*d + self.b_all[:, None]*dv))
            return F*self.inv_m[:, None]

        def evaluateDerivatives(self, t, dt, db, state, derivs):
            self.count.check("Hourglass")
            fields = list(state.vectorFields(S.HydroFieldNames.position))
            x = get_vector(fields[0], self.dim)
            add_acceleration(S, derivs, self.dim, self.acceleration(x, state_velocity(S, state, self.dim)))

        def dt(self, db, state, derivs, t):
            return (1.0e100, "Hourglass: no limit")

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "Hourglass"

    _CLASSES[S] = Hourglass
    return Hourglass


def make_hourglass(S, body, kappa, table, beta=BETA):
    return hourglass_class(S)(body, kappa, table, beta)
