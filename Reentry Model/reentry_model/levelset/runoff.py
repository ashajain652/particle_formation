"""Runoff between the level set's markers (design section 20.6; level-set stage 3 Task 6): the skyscrapers' flow sub-steps
(their skin.py `push`, `thickness_step`, `transport`, `landing_conduct`, `flow_cycles` and `flow_substeps`, matched by
amendment, never imported) on the markers' faces instead of the corner columns' dual faces.

The sheet. Every skin is a column on one marker; the markers' tessellation (`Markers.tessellate`, Ruling 2) gives the
faces between them, each a length and a unit tangent normal from i toward j; its face vector is f = length * normal [m].
`Sheet.of(markers)` puts the markers' normals, positions and faces in the skins' row order. On the sheet:
- the surface gradient of a per-marker quantity q at marker i is (1/A_i) sum_faces (q_j - q_i)/2 f_ij, projected onto
  i's tangent plane (Green-Gauss with the face value the mean, less q_i so that a constant has none however well the
  cell closes);
- the two-point weight of a face is w = length / |x_j - x_i|, and the curvature's change with the heights is
  L = -div grad, (L q)_i = -(1/A_i) sum_j w_ij (q_j - q_i).

The pushes (`pushes`; spec 6.4 of the skyscrapers, their `push`), a force per unit volume [N/m3] in each marker's
tangent plane, the same on every layer of its skin: the pressure term -dp_w/ds along the meridian away from the
stagnation point (dp_ds the shared `SurfaceFlowResult.dp_ds`), the deceleration rho_l (g - (g.n) n), the thickness term
rho_l (g.n) grad_s H and the capillary term -sigma grad_s kappa, rho_l and sigma the shared Material's liquid's. kappa is
phi's mean curvature kappa_1 + kappa_2 (design 20.6: div(grad phi / |grad phi|), from stage 1's
`criteria.principal_curvatures`), carried with the skins' heights within a heat sub-step as kappa_phi + L (H - H_ref),
H_ref the heights phi was made with -- phi itself does not move under the flow sub-steps, so the curvature's response to
the heights is the film's own. It acts only where the top flows (`flowing_tops`), and in its gradient a marker whose
top does not flow takes its flowing neighbours' mean (`_walled`).

The meridian. The flight direction v_hat is given; the stagnation point lies where the outward normal is v_hat, and the
meridian away from it at a marker is the tangent part of -v_hat normalised -- on a sphere the tangent part of
x - x_stag, as the skyscrapers' `along_meridian`.

Each flow sub-step (`flow_substeps`), in the skyscrapers' order: regrid (from the second on), the pushes, the momentum
through the depth with its response (`Skins.momentum_step(..., response=True)`), the thickness and capillary terms
implicit in the heights at the end of the transport (`thickness_step`), the transport (`transport`), and the landing
conduction of the markers metal landed on (`landing_conduct`). The count of flow sub-steps comes from the last heat
sub-step's largest Courant number and Rayleigh-Taylor growth (`flow_cycles`).

Every number here is a numerical setting of the skyscrapers' rules (named with its source); the physics is the shared
package's."""
import math
import warnings
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

from .skin import Skins

SLIVER = 1.0e-12          # a piece of a donor cell shorter than this share of it is a shared boundary's round-off, and a
                          # landing estimate below this share of the largest stays where it lands (skyscrapers' SLIVER)
LANDED_PASSES = 200       # passes at most of the landed metal's velocity estimate (skyscrapers' LANDED_PASSES)
LANDED_TOL = 1.0e-6       # it stops once no estimate moves by more than this share of the fastest (skyscrapers')
THICKNESS_PASSES = 8      # linearisations of the transport's fluxes in the thickness step at most (skyscrapers', plan)
THICKNESS_TOL = 1.0e-3    # they stop once the flux the last one missed would move no marker by more than this share of
                          # the largest height change (skyscrapers' THICKNESS_TOL, plan)
FLOW_CYCLES = 10          # flow sub-steps a heat sub-step where the skins remember no flow (skyscrapers' FLOW_CYCLES)
FLOW_CYCLES_MIN, FLOW_CYCLES_MAX = 1, 40   # the adaptive count's floor and cap (skyscrapers', plan)
FLOW_COURANT = 6.0        # the adaptive count's target largest Courant number of a flow sub-step (skyscrapers', plan)
FLOW_GROWTH = 0.02        # and its target largest Rayleigh-Taylor growth dt s of one (skyscrapers', plan)
FLOW_RELAX = 2.0          # a flow sub-step at most this many times the last heat sub-step's (skyscrapers', plan)
LANDING_DEPTH = 0.4e-3    # m below a landed skin's top that the landing conduction takes (skyscrapers', plan; the fine
                          # zone's depth)
LANDING_CELLS = 2         # cells the landing conduction takes at least (skyscrapers', plan)
TRANSPORT_TOL = 1.0e-13   # the transport's iterative solve stops once no row's residual exceeds this share of the
                          # largest |b| + |A| |x| of its right-hand side's column (skyscrapers' TRANSPORT_TOL)
TRANSPORT_CORRECTIONS = 200   # corrections of the iterative solve at most (a numerical setting of this module; the
                          # skyscrapers' 40, raised: here every cycle is iterated)
TRANSPORT_FALLBACK = False    # a solve that has not converged after TRANSPORT_CORRECTIONS raises; with this set it is
                          # finished by a direct LU instead, with a warning, and booked as not converged (stage 3 Task 9:
                          # the LU hung on 133 000 cells, so it is never taken silently)


# -- the sheet ---------------------------------------------------------------------------------------------------------
def tangent_frames(normals):
    """Two orthonormal vectors (e1, e2), each (n, 3), spanning every normal's tangent plane, e1 x e2 along the normal."""
    ref = np.where((np.abs(normals[:, 0]) < 0.9)[:, None], np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    e1 = ref - np.einsum("ij,ij->i", ref, normals)[:, None] * normals
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    return e1, np.cross(normals, e1)


def _tangent(v, normals):
    return v - np.einsum("ij,ij->i", v, normals)[:, None] * normals


@dataclass
class Sheet:
    """The markers' geometry in the skins' row order: unit normals `normal` (n, 3), positions `x` (n, 3), the faces'
    ends `i`, `j` (E,) (rows, i < j as markers), their vectors `f` (E, 3) [m] from i toward j, their two-point weights
    `w` (E,), the skins' top areas `area` (n,) [m2] and each row's marker `marker` (n,)."""
    normal: np.ndarray
    x: np.ndarray
    i: np.ndarray
    j: np.ndarray
    f: np.ndarray
    w: np.ndarray
    area: np.ndarray
    marker: np.ndarray

    @classmethod
    def of(cls, markers, skins=None):
        """The sheet of `markers` (tessellated: `markers.faces` set) carrying `skins` (default markers.skins)."""
        skins = markers.skins if skins is None else skins
        if markers.faces is None:
            raise ValueError("the markers have no faces: tessellate them first")
        n, rows = skins.n, markers.skin_index
        if len(np.unique(rows)) != len(rows) or len(rows) != n:
            raise ValueError("every skin row must belong to exactly one marker")
        marker = np.empty(n, dtype=np.int64)
        marker[rows] = np.arange(len(rows))
        normal, x = markers.normal[marker], markers.x[marker]
        fi, fj, length, unit = markers.faces
        i, j = rows[fi], rows[fj]
        f = length[:, None] * unit
        d = np.linalg.norm(x[j] - x[i], axis=1)
        return cls(normal, x, i, j, f, length / d, skins.area, marker)

    @property
    def n(self):
        return len(self.normal)

    def rows(self, per_marker):
        """A per-marker array in the skins' row order."""
        return np.asarray(per_marker)[self.marker]

    def directed(self):
        """Every face each way: senders (2E,), receivers (2E,) and the face vectors from sender to receiver (2E, 3) --
        the first E as (i, j), the rest as (j, i)."""
        return np.r_[self.i, self.j], np.r_[self.j, self.i], np.r_[self.f, -self.f]

    def gradient(self, q):
        """The surface gradient (n, 3) of q (n,) (module docstring)."""
        q = np.asarray(q, dtype=float)
        part = 0.5 * (q[self.j] - q[self.i])[:, None] * self.f
        g = np.stack([np.bincount(self.i, part[:, k], self.n) + np.bincount(self.j, part[:, k], self.n)
                      for k in range(3)], axis=1) / self.area[:, None]
        return _tangent(g, self.normal)

    def laplacian(self):
        """L = -div grad as a sparse (n, n) matrix (module docstring)."""
        W = sp.coo_matrix((np.r_[self.w, self.w], (np.r_[self.i, self.j], np.r_[self.j, self.i])),
                          shape=(self.n, self.n)).tocsr()
        return (-sp.diags(1.0 / self.area) @ (W - sp.diags(np.asarray(W.sum(axis=1)).ravel()))).tocsr()

    def incidence(self):
        """(E, n): +1 at each face's i, -1 at its j."""
        E = len(self.i)
        return sp.csr_matrix((np.r_[np.ones(E), -np.ones(E)], (np.r_[np.arange(E), np.arange(E)], np.r_[self.i, self.j])),
                             shape=(E, self.n))


def _sheet(markers):
    return markers if isinstance(markers, Sheet) else Sheet.of(markers)


# -- the forcing and the pushes ----------------------------------------------------------------------------------------
def along_meridian(normals, v_hat, magnitude):
    """Vectors of the given magnitudes (n,) in each tangent plane along the meridian away from the stagnation point (the
    point whose outward normal is the flight direction v_hat): the tangent part of -v_hat, normalised; zero where it
    vanishes (the stagnation and rear points)."""
    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
    away = -_tangent(np.broadcast_to(v, normals.shape), normals)
    size = np.linalg.norm(away, axis=1)
    unit = np.divide(away, size[:, None], out=np.zeros_like(away), where=size[:, None] > 1e-12)
    return np.asarray(magnitude, dtype=float)[:, None] * unit


@dataclass
class Forcing:
    """What drives the flow over a heat sub-step, per marker (in the markers' order): the shared surface flow's
    `dp_ds` [Pa/m] and wall shear `tau` [Pa] (its magnitude, along the meridian), the acceleration the metal feels in the
    body's frame `g` [m/s2] (3,) (the deceleration along the flight direction), the flight direction `v_hat`, the
    liquid's surface tension `sigma` [N/m], phi's mean curvature `kappa` [1/m] and, optionally, the flow sub-steps a
    heat sub-step `cycles` (None: chosen by `flow_cycles`)."""
    dp_ds: np.ndarray
    tau: np.ndarray
    g: np.ndarray
    v_hat: np.ndarray
    sigma: float
    kappa: np.ndarray
    cycles: int = None

    @classmethod
    def from_surface_flow(cls, flow, decel_vec, v_hat, sigma, phi_curvature, cycles=None):
        """From a `SurfaceFlowResult` (its dp_ds and tau per marker), the deceleration vector [m/s2], the flight
        direction, surface tension [N/m] and phi's mean curvature per marker [1/m]."""
        return cls(np.asarray(flow.dp_ds, float), np.asarray(flow.tau, float), np.asarray(decel_vec, float),
                   np.asarray(v_hat, float), float(sigma), np.asarray(phi_curvature, float), cycles)


def flowing_tops(skins, T=None):
    """Whether each skin's top cell flows (n,): where surface tension acts, not on an anchored rigid top nor a crust."""
    T = skins.T() if T is None else T
    return (skins.count > 0) & skins.flowing(T)[:, 0]


def _walled(sheet, kappa, wet):
    """kappa where a top flows; at a top that does not, its flowing neighbours' mean (the skyscrapers' `_walled`)."""
    if wet is None or wet.all():
        return kappa
    i, j, n = sheet.i, sheet.j, sheet.n
    total = np.bincount(i, kappa[j] * wet[j], n) + np.bincount(j, kappa[i] * wet[i], n)
    count = np.bincount(i, wet[j], n) + np.bincount(j, wet[i], n)
    return np.where(wet, kappa, np.where(count > 0.0, total / np.maximum(count, 1.0), kappa))


def _capillary(sheet, kappa, wet, sigma):
    G = -sigma * sheet.gradient(_walled(sheet, kappa, wet))
    return G if wet is None else np.where(wet[:, None], G, 0.0)


def _push(sheet, skins, dp_ds, g, v_hat, sigma, kappa, wet, H=None):
    """The push and its terms on the rows (`pushes`), every input per row."""
    H = skins.thickness() if H is None else H
    rho_l = skins.rho_l
    g = np.broadcast_to(np.asarray(g, dtype=float), (sheet.n, 3))
    g_n = np.einsum("ij,ij->i", g, sheet.normal)
    terms = {"pressure": -along_meridian(sheet.normal, v_hat, dp_ds),
             "deceleration": rho_l * (g - g_n[:, None] * sheet.normal),
             "thickness": rho_l * g_n[:, None] * sheet.gradient(H),
             "capillary": _capillary(sheet, kappa, wet, sigma)}
    return sum(terms.values()), terms


def pushes(markers, skins, forcing, kappa=None, wet=None):
    """The push [N/m3] (n, 3) on every layer of every skin (rows' order) and its four terms (module docstring), at the
    skins' present heights; `kappa` (n,) per row the curvature to use (default the forcing's phi curvature) and `wet`
    (n,) where the top flows (default `flowing_tops`)."""
    sheet = _sheet(markers)
    kappa = sheet.rows(forcing.kappa) if kappa is None else kappa
    wet = flowing_tops(skins) if wet is None else wet
    return _push(sheet, skins, sheet.rows(forcing.dp_ds), forcing.g, forcing.v_hat, forcing.sigma, kappa, wet)


# -- the thickness terms, implicit in the heights ----------------------------------------------------------------------
@dataclass
class FacePush:
    """What `transport` adds to the cells' velocities across each face: every cell's velocity per unit push `response`
    (n, K) and, per directed face (`Sheet.directed`'s order), the face's own two-point push less the sender's push
    dotted with the face vector, `across` (2E,) [Pa]."""
    response: np.ndarray
    across: np.ndarray


@dataclass
class ThicknessBook:
    push: np.ndarray         # N s (n, 3): the impulse of the thickness terms' implicit change on the moving cells
    wall: np.ndarray         # N s (n, 3): the impulse the anchored metal's drag gives on that change (from the wall)
    heights: np.ndarray      # m (n,): the heights predicted at the end of the transport
    faces: FacePush          # for `transport`
    growth: float            # the largest Rayleigh-Taylor growth dt s over the skins
    passes: int              # the linearisations made


def thickness_step(skins, sheet, dt, mbook, applied, g, sigma, kappa, J):
    """The push's thickness terms, rho_l (g.n) grad_s H and -sigma grad_s kappa, implicit in the heights at the end of
    the sub-step's transport (the skyscrapers' `thickness_step`, matched by amendment): called after
    `Skins.momentum_step(..., response=True)` (its book `mbook`) and before `transport`; `applied` (n, 3) [N/m3] is the
    part of the push those two terms made, `kappa` (n,) the curvature they used and J (n, n) its change per change of
    the heights (the sheet's L). The terms are the gradient of P = sigma kappa - rho_l (g.n) H; across each face P
    pushes by its two-point difference phi = w (P_i - P_j) [Pa] (surface tension only between two flowing tops, the
    deceleration with the faces' mean g.n), linear in the heights. Each sub-step: Newton on the transport's piecewise
    linear flux predicts the heights H + dH at the transport's end (at most THICKNESS_PASSES linearisations, stopping at
    THICKNESS_TOL); every moving cell's velocity is corrected by its response times the change of its column's push; and
    each face's own push goes to the transport (the column form of Rhie-Chow interpolation). Corrects U; returns a
    ThicknessBook."""
    n, K = skins.m.shape
    r = mbook.response
    T = skins.T()
    valid = skins.valid
    dz = skins.dz(T)
    A = skins.area
    H = dz.sum(axis=1)
    moving = valid & ~skins._anchored(T)
    wet = flowing_tops(skins, T)
    rho_l = skins.rho_l
    g = np.broadcast_to(np.asarray(g, dtype=float), (n, 3))
    g_n = np.einsum("ij,ij->i", g, sheet.normal)
    i, j, w = sheet.i, sheet.j, sheet.w
    g_e, both = 0.5 * (g_n[i] + g_n[j]), (wet[i] & wet[j]).astype(float)
    phi0 = w * (sigma * both * (kappa[i] - kappa[j]) - rho_l * g_e * (H[i] - H[j]))
    Bt = sheet.incidence()
    dphi = (sp.diags(w) @ (sp.diags(sigma * both) @ (Bt @ J) - sp.diags(rho_l * g_e) @ Bt)).tocsr()
    a_all, b_all, f_all = sheet.directed()
    keep = np.flatnonzero(moving.any(axis=1)[a_all])
    a, b, f = a_all[keep], b_all[keep], f_all[keep]
    phi_d = np.r_[phi0, -phi0][keep]
    dphi_d = sp.vstack([dphi, -dphi]).tocsr()[keep]
    un = np.where(moving[a], np.einsum("dkx,dx->dk", (skins.U - r[..., None] * applied[:, None, :])[a], f), 0.0)
    ra, da = r[a], np.where(moving[a], dz[a], 0.0)

    def flux(phi):
        return (np.maximum(un + ra * phi[:, None], 0.0) * da).sum(axis=1)                  # m3/s

    thick = np.where(moving, dz, 0.0).sum(axis=1)
    D = len(a)
    C = sp.csr_matrix((np.r_[np.ones(D), -np.ones(D)], (np.r_[a, b], np.r_[np.arange(D), np.arange(D)])), shape=(n, D))
    active = np.zeros(n, dtype=bool)
    active[a], active[b] = True, True
    act = np.flatnonzero(active)
    dH, phi, passes = np.zeros(n), phi_d, 0
    F = flux(phi)
    for passes in range(1, (THICKNESS_PASSES if D else 0) + 1):
        M = (ra * da * (un + ra * phi[:, None] > 0.0)).sum(axis=1)                          # m5/(N s)
        Pa = sp.csr_matrix((F / np.where(thick[a] > 0.0, thick[a], 1.0), (np.arange(D), a)), shape=(D, n))
        N = (sp.diags(A) + dt * C @ (Pa + sp.diags(M) @ dphi_d)).tocsr()[act][:, act]
        dH = np.zeros(n)
        dH[act] = spla.splu(N.tocsc()).solve(-dt * (C @ (F + M * (phi_d - phi)))[act])
        new = phi_d + dphi_d @ dH
        F_new = flux(new)
        missed = dt * np.abs(C @ (F_new - F - M * (new - phi))) / A
        phi, F = new, F_new
        if missed.max() <= THICKNESS_TOL * np.abs(dH).max():
            break
    G = rho_l * g_n[:, None] * sheet.gradient(H + dH) + _capillary(sheet, kappa + J @ dH, wet, sigma)
    dG = G - applied
    skins.U = np.where(moving[..., None], skins.U + r[..., None] * dG[:, None, :], skins.U)
    across = np.zeros(len(a_all))
    across[keep] = phi - np.einsum("dx,dx->d", G[a], f)
    V = np.where(moving, skins.m / np.where(valid, skins.density(T), 1.0), 0.0).sum(axis=1)   # as momentum_step's
    s_tension = sigma * wet
    k2 = np.pi ** 2 / A
    fastest = np.maximum(rho_l * g_n, 0.0) / (2.0 * np.maximum(s_tension, 1e-300))
    k2 = np.where(s_tension > 0.0, np.minimum(k2, fastest), k2)
    s = (r * np.where(moving, dz, 0.0)).sum(axis=1) * k2 * (rho_l * g_n - s_tension * k2)
    return ThicknessBook(push=dt * V[:, None] * dG, wall=-dt * mbook.wall_response[:, None] * dG, heights=H + dH,
                         faces=FacePush(r, across), growth=float(max(dt * s.max(), 0.0)) if n else 0.0, passes=passes)


# -- the transport -----------------------------------------------------------------------------------------------------
@dataclass
class TransportBook:
    """What a transport sub-step carried between the skins, per row: the gain from the neighbours less what was given,
    of mass [kg] and enthalpy [J] (each summing to zero to round-off), and of momentum [kg m/s] (n, 3) as it left the
    donors; `wall` [N s] (n, 3) the normal parts the projections removed (an impulse from the wall), so that each skin's
    momentum changes by momentum + wall; `courant` the largest outflow rate of any cell times the sub-step; `landed`
    (n,) where metal landed on top as a new top cell; `corrections` and `converged` the iterative solves' (`_solve`)."""
    mass: np.ndarray
    enthalpy: np.ndarray
    momentum: np.ndarray
    wall: np.ndarray
    courant: float
    landed: np.ndarray
    corrections: int = 0      # the iterative solves' corrections, the larger of the masses' and the momenta's
    converged: bool = True    # False only where TRANSPORT_FALLBACK finished a solve by the direct LU


def _fractions(depth):
    """Each row's cell boundaries (n, K+1) as fractions of the depth they hold (zero for a cell that does not count),
    0 at the top, held at 1 after the last that counts; 1 throughout where none counts."""
    z = np.concatenate([np.zeros((depth.shape[0], 1)), np.cumsum(depth, axis=1)], axis=1)
    total = z[:, -1:]
    return np.where(total > 0.0, z / np.where(total > 0.0, total, 1.0), 1.0)


def _pieces(sig, tau):
    """The pieces two partitions of [0, 1] cut each other into, row by row (sig, tau (D, K+1) non-decreasing): for every
    piece of positive length its row, its cell of sig, its cell of tau and its length. The rows are searched together,
    each offset by 4 (boundaries lie in [0, 1])."""
    Dn, B = sig.shape
    s = np.sort(np.concatenate([sig, tau], axis=1), axis=1)
    length = np.diff(s, axis=1)
    r, p = np.nonzero(length > 0.0)
    mid = 0.5 * (s[r, p] + s[r, p + 1])
    off = 4.0 * np.arange(Dn)
    ks = np.searchsorted((sig + off[:, None]).ravel(), mid + off[r], side="right") - r * B - 1
    kt = np.searchsorted((tau + off[:, None]).ravel(), mid + off[r], side="right") - r * B - 1
    return r, np.clip(ks, 0, B - 2), np.clip(kt, 0, B - 2), length[r, p]


def _upwind_order(src, dst, N):
    """An order of N unknowns in which every donor (src) comes before its receivers (dst) except within a cycle of
    transfers (a strong component), whose members stay together in their index order: the components in a topological
    order of their condensation (Kahn's algorithm, level by level)."""
    if not len(src):
        return np.arange(N)
    _, label = connected_components(sp.csr_matrix((np.ones(len(src)), (src, dst)), shape=(N, N)), directed=True,
                                    connection="strong")
    nc = int(label.max()) + 1
    ls, ld = label[src], label[dst]
    cross = ls != ld
    edges = np.unique(ls[cross] * np.int64(nc) + ld[cross])
    es, ed = edges // nc, edges % nc
    indeg = np.bincount(ed, minlength=nc)
    out = sp.csr_matrix((np.ones(len(es)), (es, ed)), shape=(nc, nc))
    level = np.full(nc, -1, dtype=np.int64)
    ready = np.flatnonzero(indeg == 0)
    k = 0
    while len(ready):
        level[ready] = k
        hit = out[ready].tocoo().col
        indeg -= np.bincount(hit, minlength=nc)
        ready = np.flatnonzero((indeg == 0) & (level < 0))
        k += 1
    return np.lexsort((np.arange(N), label, level[label]))


def _solve(A, b, order):
    """x with A x = b, b (N,) or (N, p), A an M-matrix (or one by 2x2 blocks) of the transport, `order` an upwind order
    of its unknowns. In that order A is lower triangular but for the back-couplings U within cycles where flows meet;
    M = A - U is factored with no fill (natural order, no pivoting) and x_(k+1) = M^-1 (b - U x_k) from x_0 = M^-1 b
    until no row's residual exceeds TRANSPORT_TOL of the largest |b| + |A| |x| in its column, at most
    TRANSPORT_CORRECTIONS corrections (the skyscrapers' scheme, every cycle iterated). For the masses the splitting is
    regular, so the iterates rise to the solution from below and the installed masses stay positive. Returns (x, the
    corrections made, converged). A solve still short of the tolerance raises RuntimeError, or with TRANSPORT_FALLBACK
    is finished by a direct LU with a RuntimeWarning and returned as not converged."""
    Ap = A.tocsr()[order][:, order].tocoo()
    up = Ap.col > Ap.row
    M = sp.csc_matrix((Ap.data[~up], (Ap.row[~up], Ap.col[~up])), shape=A.shape)
    U = sp.csr_matrix((Ap.data[up], (Ap.row[up], Ap.col[up])), shape=A.shape)
    bp = b[order]
    lu = spla.splu(M, permc_spec="NATURAL", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
    x = lu.solve(bp)
    count, converged = 0, True
    if U.nnz:
        absA = abs(Ap.tocsr())
        while True:
            Ux = U @ x
            r = bp - M @ x - Ux
            limit = TRANSPORT_TOL * (np.abs(bp) + absA @ np.abs(x)).max(axis=0)
            if np.all(np.abs(r).max(axis=0) <= limit):
                break
            if count >= TRANSPORT_CORRECTIONS:
                worst = float((np.abs(r).max(axis=0) / np.maximum(limit, 1e-300)).max())   # in tolerances
                message = ("the transport's iterative solve did not converge in {} corrections (residual {:.3e} of its "
                           "scale against a tolerance of {:.0e}, {:.3e} times it)".format(
                               count, worst * TRANSPORT_TOL, TRANSPORT_TOL, worst))
                if not TRANSPORT_FALLBACK:
                    raise RuntimeError(message)
                warnings.warn(message + "; finished by a direct LU", RuntimeWarning)
                x = spla.splu(Ap.tocsc()).solve(bp)
                converged = False
                break
            x = lu.solve(bp - Ux)
            count += 1
    out = np.empty_like(x)
    out[order] = x
    return out, count, converged


def _landed_velocity(sheet, run, dt, base_m, base_p, stay):
    """The velocity (n, 3) at which the metal landing on each rigid top flows on within the transport's sub-step (the
    skyscrapers' `_landed_velocity`): the mass-weighted mean of what lands there, each arrival projected onto the
    receiver's plane; zero where nothing lands, where the top flows (`run` > 0) or where crust metal landed (`stay`).
    What lands from the cells (base_m [kg], base_p [kg m/s]) is known; what landed metal carries on onto further rigid
    tops is estimated pass by pass, a marker sending dt c m / (1 + dt sum c) of what landed, c = (u.f)+ / A, until no
    velocity moves by more than LANDED_TOL of the fastest. Estimates for the rates alone; the solve carries the books."""
    n = sheet.n
    src, dst, f = sheet.directed()
    can = (run == 0) & ~stay
    onto = run[dst] == 0
    u = np.zeros((n, 3))
    m_est, p_est = base_m, base_p
    for _ in range(LANDED_PASSES):
        on = can & (m_est > SLIVER * m_est.max())
        u_new = np.where(on[:, None], p_est / np.where(on, m_est, 1.0)[:, None], 0.0)
        change = np.abs(u_new - u).max()
        u = u_new
        if change <= LANDED_TOL * np.abs(u).max():
            break
        d = np.flatnonzero(on[src])
        rate = np.maximum(np.einsum("dx,dx->d", u[src[d]], f[d]), 0.0) / sheet.area[src[d]]
        out = np.bincount(src[d], rate, n)
        sent = dt * rate * m_est[src[d]] / (1.0 + dt * out[src[d]])
        d, sent = d[onto[d]], sent[onto[d]]
        v = _tangent(u[src[d]], sheet.normal[dst[d]])
        m_est = base_m + np.bincount(dst[d], sent, n)
        p_est = base_p + np.stack([np.bincount(dst[d], sent * v[:, x], n) for x in range(3)], axis=1)
    return u


def transport(skins, sheet, dt, faces=None):
    """One linearly implicit upwind sub-step dt [s] of the metal moving between neighbouring markers at the cells'
    velocities (the skyscrapers' `transport`, matched by amendment). Across each face a moving cell k of the donor a
    sends c m_k' toward the receiver b, c = (u_k . f)+ / A_a [1/s] and m_k' its mass at the sub-step's end, with its
    specific enthalpy (and trial sensitivity) and velocity: (I + dt C) m' = m, an M-matrix whose columns sum to one, so
    the masses stay positive whatever dt. With `faces` (a FacePush from `thickness_step`) each cell's velocity across a
    face adds its response times the face's own push.

    The arrival rule: a flowing donor cell's span of fractions of the donor's flowing depth joins the receiver's
    flowing layer (its top run of flowing cells) between the same fractions, each in proportion to the overlap; where
    the receiver's top is rigid (anchored or a crust) everything lands on top as a new top cell, and so does crust
    metal wherever it goes; anchored metal neither leaves nor receives. Flowing metal landed on a rigid top flows on in
    the same solve at the velocity it arrived with (`_landed_velocity`), by fraction into a flowing receiver's whole
    flowing layer and on top of a rigid one. Arriving momentum is projected onto the receiver's tangent plane, two
    components per cell in its marker's frame; the normal parts removed are booked to the wall. Every cell is installed
    from the solve's fluxes (what it had, less what it sent, plus what arrived), so the books close exactly whatever the
    solve's own error. Updates the cells; returns a TransportBook. No regrid (the flow sub-steps regrid)."""
    n, K = skins.m.shape
    T = skins.T()
    valid = skins.valid
    flowing = skins.flowing(T)
    free = valid & ~skins._anchored(T)
    moving = free & (skins.U != 0.0).any(axis=2)
    if faces is not None:
        moving |= free & (faces.response != 0.0)
    book = TransportBook(np.zeros(n), np.zeros(n), np.zeros((n, 3)), np.zeros((n, 3)), 0.0, np.zeros(n, dtype=bool))
    if not moving.any():
        return book
    A_all, B_all, F_all = sheet.directed()
    keep = np.flatnonzero(moving.any(axis=1)[A_all])
    a, b, f = A_all[keep], B_all[keep], F_all[keep]
    across = np.einsum("dkx,dx->dk", skins.U[a], f)
    if faces is not None:
        across = across + faces.response[a] * faces.across[keep][:, None]
    c = np.where(moving[a], np.maximum(across, 0.0), 0.0) / sheet.area[a][:, None]
    sel = (c > 0.0).any(axis=1)
    a, b, c = a[sel], b[sel], c[sel]
    if not len(a):
        return book
    # where it joins: by fraction of the flowing depth, or on top
    dz = skins.dz(T)
    run = np.cumprod(flowing, axis=1).sum(axis=1)
    sig = _fractions(np.where(flowing, dz, 0.0))
    tau = _fractions(np.where(np.arange(K)[None, :] < run[:, None], dz, 0.0))
    joins = flowing[a] & (run[b] > 0)[:, None] & (c > 0.0)
    rows = np.flatnonzero(joins.any(axis=1))
    r, kd, kr, length = _pieces(sig[a[rows]], tau[b[rows]])
    ok = joins[rows[r], kd]
    edge, kd, kr, length = rows[r[ok]], kd[ok], kr[ok], length[ok]
    width = sig[a[edge], kd + 1] - sig[a[edge], kd]
    ok = length > SLIVER * width
    edge, kd, kr, length = edge[ok], kd[ok], kr[ok], length[ok]
    group = edge * K + kd
    share = length / np.bincount(group, length, len(a) * K)[group]          # each donor cell's pieces sum to one
    te, tk = np.nonzero((c > 0.0) & ~joins)
    # every arrival: its sending and receiving rows, receiving unknown (n K + b: b's new top cell), sending unknown, dt c
    e_src, e_dst = np.r_[a[edge], a[te]], np.r_[b[edge], b[te]]
    e_row = np.r_[b[edge] * K + kr, n * K + b[te]]
    e_col = np.r_[a[edge] * K + kd, a[te] * K + tk]
    e_val = dt * np.r_[c[edge, kd] * share, c[te, tk]]
    # landed metal flows on in the same solve at the velocity it arrived with
    lands = flowing[a[te], tk] & (run[b[te]] == 0)
    stay = np.zeros(n, dtype=bool)
    stay[b[te[~flowing[a[te], tk]]]] = True
    out = np.bincount(e_col, e_val, n * K)[a[te] * K + tk]
    w_land = (dt * c[te, tk] * skins.m[a[te], tk] / (1.0 + out))[lands]
    rec = b[te[lands]]
    v = _tangent(skins.U[a[te], tk][lands], sheet.normal[rec])
    u_on = _landed_velocity(sheet, run, dt, np.bincount(rec, w_land, n),
                            np.stack([np.bincount(rec, w_land * v[:, x], n) for x in range(3)], axis=1), stay)
    on = np.flatnonzero((u_on != 0.0).any(axis=1)[A_all])
    c_on = np.maximum(np.einsum("dx,dx->d", u_on[A_all[on]], F_all[on]), 0.0) / sheet.area[A_all[on]]
    on, c_on = on[c_on > 0.0], c_on[c_on > 0.0]
    reach = np.zeros(n, dtype=bool)                          # the new top cells metal reaches: only those send
    reach[rec] = True
    while len(on):
        more = reach[A_all[on]] & (run[B_all[on]] == 0) & ~reach[B_all[on]]
        if not more.any():
            break
        reach[B_all[on][more]] = True
    on, c_on = on[reach[A_all[on]]], c_on[reach[A_all[on]]]
    if len(on):
        s_, r_ = A_all[on], B_all[on]
        into = np.flatnonzero(run[r_] > 0)                   # by fraction of the receiver's whole flowing layer
        piece = np.diff(tau[r_[into]], axis=1) * (np.arange(K)[None, :] < run[r_[into]][:, None])
        q, kq = np.nonzero(piece > 0.0)
        onto = np.flatnonzero(run[r_] == 0)                  # on top of a rigid receiver, its new top cell
        e_src = np.r_[e_src, s_[into[q]], s_[onto]]
        e_dst = np.r_[e_dst, r_[into[q]], r_[onto]]
        e_row = np.r_[e_row, r_[into[q]] * K + kq, n * K + r_[onto]]
        e_col = np.r_[e_col, n * K + s_[into[q]], n * K + s_[onto]]
        e_val = np.r_[e_val, dt * c_on[into[q]] * piece[q, kq], dt * c_on[onto]]
    # the unknowns that send or receive
    active = np.zeros(n * K + n, dtype=bool)
    active[e_row], active[e_col] = True, True
    act = np.flatnonzero(active)
    N = len(act)
    local = np.full(n * K + n, -1, dtype=np.int64)
    local[act] = np.arange(N)
    lr, lc = local[e_row], local[e_col]
    diag = 1.0 + np.bincount(lc, e_val, N)
    mat = sp.csc_matrix((np.r_[diag, -e_val], (np.r_[np.arange(N), lr], np.r_[np.arange(N), lc])), shape=(N, N))
    order = _upwind_order(lc, lr, N)

    def cell(v):
        return np.concatenate([v.reshape((n * K,) + v.shape[2:]), np.zeros((n,) + v.shape[2:])])[act]

    trial = skins.dH is not None
    b0 = np.column_stack([cell(skins.m), cell(skins.H)] + ([cell(skins.dH)] if trial else []))
    x, count_m, ok_m = _solve(mat, b0, order)
    moved = e_val[:, None] * x[lc]
    x = b0 + np.stack([np.bincount(lr, moved[:, k], N) - np.bincount(lc, moved[:, k], N) for k in range(x.shape[1])],
                      axis=1)
    # the momentum, each arrival projected onto its receiver's plane, two components per cell in its own frame
    owner = np.where(act < n * K, act // K, act - n * K)
    e1, e2 = tangent_frames(sheet.normal)
    P = cell(skins.m[..., None] * skins.U)
    p0 = np.column_stack([np.einsum("ij,ij->i", P, e1[owner]), np.einsum("ij,ij->i", P, e2[owner])])
    held = P - p0[:, :1] * e1[owner] - p0[:, 1:] * e2[owner]                 # a normal part a cell held: the wall's
    Re = np.stack([np.stack([np.einsum("ij,ij->i", er[e_dst], es[e_src]) for es in (e1, e2)], axis=1)
                   for er in (e1, e2)], axis=1)                              # [arrival, receiver's comp, donor's comp]
    rows2 = np.concatenate([2 * np.arange(N), 2 * np.arange(N) + 1] + [2 * lr + qq for qq in (0, 0, 1, 1)])
    cols2 = np.concatenate([2 * np.arange(N), 2 * np.arange(N) + 1] + [2 * lc + ss for ss in (0, 1, 0, 1)])
    vals2 = np.concatenate([diag, diag] + [-e_val * Re[:, qq, ss] for qq, ss in ((0, 0), (0, 1), (1, 0), (1, 1))])
    order2 = np.column_stack([2 * order, 2 * order + 1]).ravel()
    y, count_p, ok_p = _solve(sp.csc_matrix((vals2, (rows2, cols2)), shape=(2 * N, 2 * N)), p0.reshape(-1), order2)
    book.corrections, book.converged = max(count_m, count_p), bool(ok_m and ok_p)
    sent_p = e_val[:, None] * y.reshape(N, 2)[lc]
    got = np.einsum("eqs,es->eq", Re, sent_p)
    p = p0 + np.stack([np.bincount(lr, got[:, k], N) - np.bincount(lc, sent_p[:, k], N) for k in range(2)], axis=1)
    P_new = p[:, :1] * e1[owner] + p[:, 1:] * e2[owner]
    # the books, arrival by arrival
    carried = sent_p[:, :1] * e1[e_src] + sent_p[:, 1:] * e2[e_src]
    lost = np.einsum("ij,ij->i", carried, sheet.normal[e_dst])[:, None] * sheet.normal[e_dst]

    def net(v):
        return np.bincount(e_dst, v, n) - np.bincount(e_src, v, n)

    book.mass, book.enthalpy = net(moved[:, 0]), net(moved[:, 1])
    book.momentum = np.stack([net(carried[:, k]) for k in range(3)], axis=1)
    book.wall = -np.stack([np.bincount(e_dst, lost[:, k], n) + np.bincount(owner, held[:, k], n) for k in range(3)],
                          axis=1)
    book.courant = float((diag - 1.0).max())
    # install: the cells solved for, and a new top cell where metal landed on top
    old = act < n * K
    m, H, U = skins.m.copy(), skins.H.copy(), skins.U.copy()
    dH = skins.dH.copy() if trial else None
    m.ravel()[act[old]], H.ravel()[act[old]] = x[old, 0], x[old, 1]
    U.reshape(-1, 3)[act[old]] = P_new[old] / x[old, :1]
    if trial:
        dH.ravel()[act[old]] = x[old, 2]
    count = skins.count.copy()
    tops = act[~old] - n * K
    if len(tops):
        new = np.zeros(n, dtype=bool)
        new[tops] = True
        src = np.clip(np.arange(K + 1)[None, :] - new[:, None], 0, K)

        def down(v):
            padded = np.concatenate([v, np.zeros((n, 1) + v.shape[2:])], axis=1)
            return np.take_along_axis(padded, src.reshape(src.shape + (1,) * (v.ndim - 2)), axis=1)

        m, H, U = down(m), down(H), down(U)
        m[tops, 0], H[tops, 0], U[tops, 0] = x[~old, 0], x[~old, 1], P_new[~old] / x[~old, :1]
        if trial:
            dH = down(dH)
            dH[tops, 0] = x[~old, 2]
        count = count + new
        book.landed = new
    skins.m, skins.H, skins.U, skins.count = m, H, U, count
    if trial:
        skins.dH = dH
    skins._trim()
    return book


# -- landing conduction ------------------------------------------------------------------------------------------------
def landing_conduct(skins, landed, dt, depth=None):
    """Metal that landed on a rigid top loses its heat to the metal beneath within the flow sub-steps (the skyscrapers'
    `landing_conduct`): the cells of each skin in `landed` (n,) whose tops lie within `depth` (LANDING_DEPTH) of its top,
    at least LANDING_CELLS, are conducted for dt [s] implicitly (`Skins.conduct` on those cells alone), with no heat
    across the block's top (its heating and radiation are the heat sub-step's) nor its bottom; what the Newton iteration
    leaves of the block's energy is kept by its bottom cell, so the books are exact. An open trial pass's sensitivity is
    carried with the heat. Updates the cells' enthalpies."""
    rows = np.flatnonzero(landed)
    if not rows.size:
        return
    depth = LANDING_DEPTH if depth is None else depth
    T = skins.T()
    dz = skins.dz(T)[rows]
    valid = skins.valid[rows]
    count = skins.count[rows]
    k = np.maximum((valid & (np.cumsum(dz, axis=1) - dz < depth)).sum(axis=1), np.minimum(LANDING_CELLS, count))
    width = max(int(k.max()), 1)
    inside = np.arange(width)[None, :] < k[:, None]
    block = Skins.__new__(Skins)                                    # the blocks alone, sharing the skins' tables
    block.material, block.settings, block.props = skins.material, skins.settings, skins.props
    block.rho_s, block.rho_l = skins.rho_s, skins.rho_l
    block.area, block.kappa = skins.area[rows].copy(), skins.kappa[rows].copy()
    block.momentum_out = np.zeros((len(rows), 3))
    block.dH = None
    block.install(np.where(inside, skins.m[rows, :width], 0.0), np.where(inside, skins.H[rows, :width], 0.0), k,
                  np.where(inside[..., None], skins.U[rows, :width], 0.0))
    trial = skins.dH is not None
    if trial:
        block.dH = np.where(inside, skins.dH[rows, :width], 0.0)
        block._trial_E, block._trial_dE = np.zeros(len(rows)), np.zeros(len(rows))
    zero = np.zeros(len(rows))
    heat = block.conduct(dt, zero, zero, zero, None, trial=trial)
    last = np.maximum(k - 1, 0)
    r = np.arange(len(rows))
    block.H[r, last] += heat.E_base                                # the iteration's error, kept by the block
    cells = np.nonzero(inside)
    skins.H[rows[cells[0]], cells[1]] = block.H[cells]
    if trial:
        skins.dH[rows[cells[0]], cells[1]] = block.dH[cells]


# -- the flow sub-steps ------------------------------------------------------------------------------------------------
def flow_cycles(memory, dt):
    """The number of flow sub-steps in a heat sub-step dt [s] from the skins' flow memory (the last heat sub-step's
    largest Courant number C, Rayleigh-Taylor growth g and flow sub-step h_0): the flow sub-step is the shortest of
    FLOW_RELAX h_0, h_0 FLOW_COURANT / C and h_0 (FLOW_GROWTH / g)^(1/2), the count dt over it rounded up within
    FLOW_CYCLES_MIN and FLOW_CYCLES_MAX; FLOW_CYCLES where the skins remember no flow (the skyscrapers' `flow_cycles`)."""
    courant, growth, length = (float(v) for v in memory)
    if length <= 0.0:
        return FLOW_CYCLES
    h = FLOW_RELAX * length
    if courant > 0.0:
        h = min(h, length * FLOW_COURANT / courant)
    if growth > 0.0:
        h = min(h, length * math.sqrt(FLOW_GROWTH / growth))
    return int(min(max(math.ceil(dt / h - 1e-9), FLOW_CYCLES_MIN), FLOW_CYCLES_MAX))


@dataclass
class FlowBook:
    """What the flow sub-steps did, per row: impulses [N s] and momentum [kg m/s] (n, 3) adding up to the skins'
    momentum change -- the push's (its implicit thickness change included), the gas shear's, the wall's (an impulse
    from the wall: the anchored metal's drag, the momentum cells held on joining it, the normal parts the transport's
    projections removed) and what the transport carried in less what it sent (zero over the sheet), less the momentum
    the release took (`released_p`) -- with what the release took and the sub-steps' worst numbers."""
    push: np.ndarray
    shear: np.ndarray
    wall: np.ndarray
    carried: np.ndarray
    courant: float = 0.0
    growth: float = 0.0
    passes: int = 0
    sub_steps: int = 0
    cycles: int = 0
    released: np.ndarray = None      # kg (n,): what the release took from the skins' tops (zeros without a release)
    released_H: np.ndarray = None    # J (n,): its enthalpy
    released_p: np.ndarray = None    # kg m/s (n, 3): its momentum, which leaves with the droplets
    corrections: int = 0             # the transport's most corrections in one solve
    unconverged: int = 0             # transport solves finished by TRANSPORT_FALLBACK's direct LU

    def __post_init__(self):
        n = len(self.push)
        if self.released is None:
            self.released = np.zeros(n)
        if self.released_H is None:
            self.released_H = np.zeros(n)
        if self.released_p is None:
            self.released_p = np.zeros((n, 3))

    @classmethod
    def zero(cls, n):
        return cls(*(np.zeros((n, 3)) for _ in range(4)))

    def add(self, other):
        """Add another book's impulses and releases and keep the worst numbers of both."""
        for key in ("push", "shear", "wall", "carried", "released", "released_H", "released_p"):
            setattr(self, key, getattr(self, key) + getattr(other, key))
        self.courant, self.growth = max(self.courant, other.courant), max(self.growth, other.growth)
        self.passes, self.sub_steps = max(self.passes, other.passes), self.sub_steps + other.sub_steps
        self.cycles = max(self.cycles, other.cycles)
        self.corrections = max(self.corrections, other.corrections)
        self.unconverged += other.unconverged


def flow_substeps(skins, markers, forcing, dt, cycles=None, H_ref=None, release=None, after=None):
    """The flow sub-steps of one heat sub-step dt [s] (module docstring): `cycles` of them (default forcing.cycles, or
    where that is None the count `flow_cycles` chooses from the skins' `flow_memory`, which they then update), each the
    regrid (from the second on), `pushes` with phi's curvature carried with the heights (kappa_phi + L (H - H_ref), H_ref
    the heights phi was made with, default those at the start), `Skins.momentum_step` with its response,
    `thickness_step`, `transport`, the `landing_conduct` of the markers metal landed on this heat sub-step and, with
    `release` (a callable `release(skins, h)` returning a book with per-row `mass`, `enthalpy` and `momentum`, such as
    `spraying.SprayPlan.release`), the release from the tops (stage 3 Task 9: every flow sub-step, after the landing
    conduction, as the skyscrapers' flow sub-steps release), and `after(skins)`, if given, right after each transport,
    before the landing conduction reads the new depths (the body's per-column checks, stage 3 Task 9: a transport
    deepened a column past its curvature's room within one sub-step). Where nothing flows nothing moves and no sub-step is taken:
    any velocity left is booked to the wall and cleared, and the skins forget the flow. The heat sub-step's own regrid
    follows (the caller's). Returns a FlowBook."""
    sheet = _sheet(markers)
    n = skins.n
    book = FlowBook.zero(n)
    memory = skins.flow_memory
    adaptive = forcing.cycles is None and cycles is None
    if not skins.flowing().any():
        book.wall -= skins.totals()[2]
        skins.U = np.zeros_like(skins.U)
        if adaptive:
            skins.flow_memory = np.zeros(3)
        return book
    if cycles is None:
        cycles = flow_cycles(memory, dt) if forcing.cycles is None else int(forcing.cycles)
    h = dt / cycles
    dp_ds, kappa_phi = sheet.rows(forcing.dp_ds), sheet.rows(forcing.kappa)
    shear = along_meridian(sheet.normal, forcing.v_hat, sheet.rows(forcing.tau))
    L = sheet.laplacian()
    H_ref = skins.thickness() if H_ref is None else np.asarray(H_ref, dtype=float)
    landed = np.zeros(n, dtype=bool)
    for k in range(cycles):
        if k:
            skins.regrid()
        H = skins.thickness()
        wet = flowing_tops(skins)
        kappa = kappa_phi + L @ (H - H_ref)
        G, terms = _push(sheet, skins, dp_ds, forcing.g, forcing.v_hat, forcing.sigma, kappa, wet, H)
        mb = skins.momentum_step(h, shear, G, response=True)
        tb = thickness_step(skins, sheet, h, mb, terms["thickness"] + terms["capillary"], forcing.g, forcing.sigma,
                            kappa, L)
        tr = transport(skins, sheet, h, faces=tb.faces)
        if after is not None:                       # the body's per-column checks before anything reads the new depths
            after(skins)
        landed |= tr.landed
        if landed.any():
            landing_conduct(skins, landed, h)
        step = FlowBook(mb.push + tb.push, mb.shear, -mb.wall + tb.wall + tr.wall, tr.momentum, tr.courant,
                        tb.growth, tb.passes, 1, corrections=tr.corrections, unconverged=int(not tr.converged))
        if release is not None:
            rb = release(skins, h)
            step.released, step.released_H = np.asarray(rb.mass, float), np.asarray(rb.enthalpy, float)
            step.released_p = np.asarray(rb.momentum, float)
        book.add(step)
    book.cycles = cycles
    if adaptive:
        skins.flow_memory = np.array([book.courant, book.growth, h])
    return book
