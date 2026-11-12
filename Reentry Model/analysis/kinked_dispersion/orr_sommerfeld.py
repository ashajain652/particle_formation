"""The viscous check of the kinked-profile dispersion (spec stage 6, Task 10; analysis only): the Orr-Sommerfeld
eigenvalue problem of a layered melt column, solved by Chebyshev collocation in each layer.

The column is a stack of cells, top first: faces y (depths in m, measured down from the top), the base speed U at each
cell's centre (m/s), and per cell a viscosity mu (Pa s) and a density rho (kg/m3). The bottom is a no-slip wall; the top
is a free surface with surface tension sigma (N/m, `top="free"`) or a moving no-slip wall (`top="wall"`, for the Couette
validations). The gas is neglected, as in Girin's theory and the broken-line solver, and so are gravity and the
deceleration. The disturbance sees each cell's viscosity as Newtonian, at the cell's base shear rate.

The base profile (`base_profile`). Inside each cell the base speed is straight in each half-cell, through the speeds at
its top face, its centre and its bottom face. By default the speed at a face between two cells is the skins' own:
stress continuity across the two half-cells (`skin.series_drag`), u_f = (w_a u_a + w_b u_b) / (w_a + w_b) with
w = 2 mu / dz; the surface speed is extrapolated linearly from the top two centres (`melt_profile`'s fit does the same);
the bottom wall is at rest. The shear stress is then continuous at every face, and the kink at a cell's centre, where
the half-cells' slopes differ, is where the skins' finite volumes lump the push (pressure gradient and deceleration)
that bends the profile. `u_faces` replaces all of this with given face speeds (the validations use exact profiles).

The disturbance. With y pointing up, the stream function phi(y) exp[ik(x - ct)] and its vorticity
omega = (D^2 - k^2) phi obey, in each region of constant mu and rho and straight U (U'' = 0),

    mu (D^2 - k^2) omega - i k rho U omega = -i k rho c omega,

the Orr-Sommerfeld equation split into two second-order ones; the growth rate of the wave is k Im c (1/s). Where two
regions meet:

- a face where mu or rho changes is a material interface, displaced by eta = -phi / (U - c). Normal velocity,
  tangential velocity, tangential stress and normal stress are continuous there:
  [phi] = 0;  (U - c)[phi'] = [U'] phi;  [mu (phi'' + k^2 phi)] = 0;  [i k rho (U - c) phi' - i k rho U' phi
  - mu phi''' + 3 k^2 mu phi'] = 0. The base shear stress mu U' is continuous at such a face by construction, so the
  base terms the displacement brings in cancel;
- a kink inside uniform fluid (a cell's centre, or a face between equal cells, such as Girin's kink) is the Orr-Sommerfeld
  equation itself with U'' = [U'] delta(y - y_kink): phi, phi' and phi'' continuous and mu [phi'''] = -i k rho [U'] phi;
- the free surface: zero tangential stress mu (phi'' + k^2 phi) = 0, the kinematic condition (U - c) eta = -phi, and the
  normal stress i k rho (U - c) phi' - i k rho U' phi - mu phi''' + 3 k^2 mu phi' + i sigma k^3 eta = 0 (in the
  inviscid limit the broken-line solver's rho [-(U - c) phi' + U' phi] = sigma k^2 eta); the gas carries the base shear
  stress, so the tilt of the surface adds nothing;
- a no-slip wall: phi = phi' = 0.
In these conditions phi'' = omega + k^2 phi and phi''' = omega' + k^2 phi'.

The method. Each region is mapped onto Chebyshev-Gauss-Lobatto nodes, with phi and omega both unknowns at every node;
both equations are collocated at the interior nodes, and the end nodes carry the conditions (Dongarra, Straughan &
Walker's 1996 "D^2" form, which needs no fourth derivative: built with D^4, the same collocation lost up to 0.4 % of the
growth to round-off on Girin's profile at Re 1e5, from the cancellation of D^4's large entries). The pencil (A, B) is
solved by shift and invert: with R the rows linear in c and E the identity's columns at R, the nonzero eigenvalues of
(A - s B)^-1 B are those of B_R (A - s B)^-1 E, so one LU and a standard eigenproblem of the size of R give
c = s + 1 / mu; the shift s lies above the unstable region, at the middle of the speed range plus SHIFT_UP speed scales
times i, where the eigenvalues that matter are computed accurately however badly B is conditioned. Above DENSE_MAX
unknowns (a flight column of 70 cells has about 2500) ARPACK finds the eigenvalues nearest the shift instead, with a
check that it reached far enough (`_solve`). Both agree with scipy's QZ on the same pencil (tested).

Each region's degree follows its local scales -- its depth against the wavelength, the viscous diffusion length
(nu / (k V))^(1/2) and the critical-layer length (nu / (k |U'|))^(1/3), V the column's speed scale -- and a region that
needs more than N_MAX is split. Spurious eigenvalues are filtered twice: any growth faster than GROWTH_BOUND times the
column's speed scale is discarded (the semicircle theorem bounds Im c by half the speed range in the inviscid problem),
and the chosen eigenvalue's eigenfunction must be resolved (the last two Chebyshev coefficients of phi in every region
below TAIL_TOL of the largest). `resolution=2` doubles every region's degree: the resolution check.

Validated (`tests/test_kinked_dispersion_orr_sommerfeld.py`): Lamb's damped capillary wave, Charru & Hinch's long-wave
sign with Yih's wave speed, Hooper & Boyd's short-wave growth and Girin's profile at Re 1e4 and 1e5 against the
broken-line solver. This module changes no model code."""
import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import NamedTuple

import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

N_MIN = 6               # the lowest Chebyshev degree of a region
N_PER = 2.0             # degree per unit of a region's depth over its local scale
N_MAX = 64              # a region needing more is split into equal parts
KINK_TOL = 1e-9         # a slope jump below this times the column's slope scale is no kink
RHO_TOL = 0.01          # a density jump below this fraction is ignored in the interface conditions (no gravity)
GROWTH_BOUND = 2.0      # an eigenvalue with Im c above this times the speed scale is spurious
TAIL_TOL = 1e-3         # an eigenfunction is resolved where the last two Chebyshev coefficients of phi in every region
                        # are below this times its largest coefficient
CANDIDATES = 6          # eigenvalues tried, from the top, for a resolved eigenfunction
SHIFT_UP = 1.0          # the shift lies this many speed scales above the real axis
SURFACE_SHEAR = "stress"  # how the gas's base shear acts on the tilted surface: "stress", the gas's stress tensor
                          # carries the base shear, so the tilt adds nothing (the default, as at every internal face);
                          # "traction", the gas pulls along the tangent only, adding 2 k^2 mu U' eta to the normal stress
                          # row (a sensitivity: the two differ by a term of relative order 1 / Reynolds number)
DENSE_MAX = 1500        # unknowns up to which the eigenproblem is solved densely; ARPACK above
ARPACK_NEV = 120        # eigenvalues ARPACK starts with
ARPACK_REACH = 2.0      # ARPACK's farthest eigenvalue must lie at least this many speed scales from the shift


@lru_cache(maxsize=None)
def _cheb(N):
    """Chebyshev-Gauss-Lobatto nodes x_j = cos(pi j / N), j = 0..N (x_0 = 1), the differentiation matrix (Trefethen
    2000, cheb.m) and the matrix from node values to Chebyshev coefficients."""
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    c = np.where((j == 0) | (j == N), 2.0, 1.0) * (-1.0) ** j
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1.0 / c) / (X - X.T + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    T = np.cos(np.outer(np.arccos(np.clip(x, -1.0, 1.0)), j))      # T_n(x_j)
    return x, D, np.linalg.inv(T)


@dataclass
class Profile:
    """A column's base profile: faces y (L + 1) m downward from 0, the speeds at the faces u_faces and the centres U
    (m/s), mu and rho per cell, `kinked`, the cells whose two half-cells' slopes differ (split at the centre), and the
    top ("free" or "wall")."""
    y: np.ndarray
    u_faces: np.ndarray
    U: np.ndarray
    mu: np.ndarray
    rho: np.ndarray
    kinked: np.ndarray
    top: str


def base_profile(y, U, mu, rho, top="free", u_faces=None):
    """The base profile of a column (module docstring); `u_faces` (L + 1) m/s overrides the default face speeds."""
    y = np.asarray(y, dtype=float)
    U, mu, rho = (np.asarray(a, dtype=float).reshape(-1) for a in (U, mu, rho))
    L = U.size
    if y.size != L + 1 or mu.size != L or rho.size != L:
        raise ValueError("y needs one more entry than U, mu and rho")
    if (np.diff(y) <= 0.0).any():
        raise ValueError("the faces must deepen downward")
    if not (np.isfinite(U).all() and np.isfinite(mu).all() and np.isfinite(rho).all()) or (mu <= 0).any() \
            or (rho <= 0).any():
        raise ValueError("U, mu and rho must be finite, and mu and rho positive")
    if top not in ("free", "wall"):
        raise ValueError("top is 'free' or 'wall'")
    y = y - y[0]
    dz = np.diff(y)
    if u_faces is None:
        uf = np.empty(L + 1)
        w = 2.0 * mu / dz
        uf[1:-1] = (w[:-1] * U[:-1] + w[1:] * U[1:]) / (w[:-1] + w[1:])
        uf[-1] = 0.0
        if L >= 2:
            zc = y[:-1] + 0.5 * dz
            uf[0] = U[0] - zc[0] * (U[1] - U[0]) / (zc[1] - zc[0])
        else:
            uf[0] = 2.0 * U[0] - uf[-1]
    else:
        uf = np.asarray(u_faces, dtype=float).copy()
        if uf.size != L + 1:
            raise ValueError("u_faces needs one entry per face")
    s_up, s_lo = _half_slopes(y, uf, U)
    scale = max(np.max(np.abs(np.concatenate([s_up, s_lo]))), (uf.max() - uf.min()) / y[-1], 1e-300)
    kinked = np.abs(s_up - s_lo) > KINK_TOL * scale
    return Profile(y=y, u_faces=uf, U=U, mu=mu, rho=rho, kinked=kinked, top=top)


def _half_slopes(y, uf, U):
    """dU/dy (y up, 1/s) of each cell's upper and lower half."""
    half = 0.5 * np.diff(y)
    return (uf[:-1] - U) / half, (U - uf[1:]) / half


@dataclass
class _Region:
    top: float          # y (up) of the region's top, m
    bot: float          # y of its bottom
    u: float            # speed at its top, m/s
    s: float            # slope dU/dy, 1/s
    mu: float
    rho: float
    N: int = N_MIN
    below: str = "wall"  # how it meets the region below: 'kink' (uniform fluid), 'material', or 'wall'


def _regions(prof):
    """The column's regions of straight speed and uniform fluid, top first, before the degrees are set."""
    y, uf, U = prof.y, prof.u_faces, prof.U
    s_up, s_lo = _half_slopes(y, uf, U)
    out = []
    L = U.size
    for i in range(L):
        d0, d1 = y[i], y[i + 1]
        if prof.kinked[i]:
            dm = 0.5 * (d0 + d1)
            out.append(_Region(-d0, -dm, uf[i], s_up[i], prof.mu[i], prof.rho[i], below="kink"))
            out.append(_Region(-dm, -d1, U[i], s_lo[i], prof.mu[i], prof.rho[i]))
        else:
            out.append(_Region(-d0, -d1, uf[i], (uf[i] - uf[i + 1]) / (d1 - d0), prof.mu[i], prof.rho[i]))
        if i < L - 1:
            same = prof.mu[i] == prof.mu[i + 1] and prof.rho[i] == prof.rho[i + 1]
            out[-1].below = "kink" if same else "material"
    return out


def _speed_scale(prof, k, sigma):
    """The column's speed scale (m/s): its speed range, or the capillary wave's speed (sigma k / rho)^(1/2) at a free
    surface if that is larger."""
    speeds = np.concatenate([prof.u_faces, prof.U])
    V = speeds.max() - speeds.min()
    if prof.top == "free" and sigma > 0.0:
        V = max(V, np.sqrt(sigma * k / prof.rho[0]))
    return max(V, 1e-300)


def _degrees(regions, k, V, resolution):
    """Each region's Chebyshev degree from its local scales, splitting those that need more than N_MAX."""
    out = []
    for r in regions:
        h = r.top - r.bot
        nu = r.mu / r.rho
        need = k * h + 2.0 * np.sqrt(h / np.sqrt(nu / (k * V)))
        if r.s != 0.0:
            need += h / (nu / (k * abs(r.s))) ** (1.0 / 3.0)
        total = (N_MIN + N_PER * need) * resolution
        parts = max(1, int(np.ceil(total / N_MAX)))
        N = max(4, int(np.ceil(total / parts)))
        edges = np.linspace(r.top, r.bot, parts + 1)
        for p in range(parts):
            out.append(_Region(edges[p], edges[p + 1], r.u + r.s * (edges[p] - r.top), r.s, r.mu, r.rho, N=N,
                               below=r.below if p == parts - 1 else "kink"))
    return out


def _assemble(regions, k, sigma, top):
    """The pencil (A, B) on the unknowns phi and omega at every region's nodes (region by region, top first, each region
    its N + 1 phi then its N + 1 omega, node 0 at its top) and, at a free surface, eta last; with the regions' offsets."""
    offsets = np.cumsum([0] + [2 * (r.N + 1) for r in regions])
    n = offsets[-1] + (1 if top == "free" else 0)
    A = np.zeros((n, n), dtype=complex)
    B = np.zeros((n, n), dtype=complex)
    ik = 1j * k
    k2 = k * k
    row = 0
    D1s = []
    for i, r in enumerate(regions):
        x, D, _ = _cheb(r.N)
        h = r.top - r.bot
        yn = r.bot + 0.5 * (x + 1.0) * h
        D1 = D * (2.0 / h)
        D2 = D1 @ D1
        D1s.append(D1)
        m = r.N + 1
        o = offsets[i]
        P, W = slice(o, o + m), slice(o + m, o + 2 * m)
        inner = np.arange(1, r.N)
        nr = inner.size
        I = np.eye(m)
        # omega - (D^2 - k^2) phi = 0
        A[row:row + nr, W] = I[inner]
        A[row:row + nr, P] = -(D2 - k2 * I)[inner]
        row += nr
        # mu (D^2 - k^2) omega - i k rho U omega = -i k rho c omega
        Un = r.u + r.s * (yn - r.top)
        A[row:row + nr, W] = r.mu * (D2 - k2 * I)[inner] - (ik * r.rho * Un[:, None] * I)[inner]
        B[row:row + nr, W] = -ik * r.rho * I[inner]
        row += nr

    def term(i, end, what, coef):
        """(column indices, values) of coef times a quantity at region i's top (end 0) or bottom (end -1): what is
        'phi', 'phi1' (phi'), 'phi2' (phi''), 'phi3' (phi'''), 'om' or 'om1' (omega')."""
        m = regions[i].N + 1
        o = offsets[i]
        j = 0 if end == 0 else m - 1
        P = np.arange(o, o + m)
        Wc = np.arange(o + m, o + 2 * m)
        d = D1s[i][j]
        if what == "phi":
            return [(np.array([o + j]), np.array([coef]))]
        if what == "phi1":
            return [(P, coef * d)]
        if what == "om":
            return [(np.array([o + m + j]), np.array([coef]))]
        if what == "om1":
            return [(Wc, coef * d)]
        if what == "phi2":                                  # omega + k^2 phi
            return term(i, end, "om", coef) + term(i, end, "phi", coef * k2)
        if what == "phi3":                                  # omega' + k^2 phi'
            return term(i, end, "om1", coef) + term(i, end, "phi1", coef * k2)
        raise ValueError(what)

    def put(M, rr, parts):
        for cols, vals in parts:
            np.add.at(M[rr], cols, vals)

    if top == "free":
        r0 = regions[0]
        e = n - 1
        put(A, row, term(0, 0, "phi2", r0.mu) + term(0, 0, "phi", r0.mu * k2))          # tangential stress
        row += 1
        A[row, e] = r0.u                                                                 # kinematic
        put(A, row, term(0, 0, "phi", 1.0))
        B[row, e] = 1.0
        row += 1
        put(A, row, term(0, 0, "phi1", ik * r0.rho * r0.u + 3.0 * k2 * r0.mu) + term(0, 0, "phi", -ik * r0.rho * r0.s)
            + term(0, 0, "phi3", -r0.mu))                                                # normal stress
        A[row, e] = 1j * sigma * k ** 3 + (2.0 * k2 * r0.mu * r0.s if SURFACE_SHEAR == "traction" else 0.0)
        put(B, row, term(0, 0, "phi1", ik * r0.rho))
        row += 1
    else:
        put(A, row, term(0, 0, "phi", 1.0))
        put(A, row + 1, term(0, 0, "phi1", 1.0))
        row += 2
    s_scale = max(max(abs(r.s) for r in regions), 1e-300)
    for i in range(len(regions) - 1):
        a, b = regions[i], regions[i + 1]
        uf = b.u
        ds = a.s - b.s
        kink = abs(ds) > KINK_TOL * s_scale
        if a.below == "kink":
            for what in ("phi", "phi1", "phi2"):
                put(A, row, term(i, -1, what, 1.0) + term(i + 1, 0, what, -1.0))
                row += 1
            put(A, row, term(i, -1, "phi3", a.mu) + term(i + 1, 0, "phi3", -a.mu) + term(i, -1, "phi", ik * a.rho * ds))
            row += 1
            continue
        put(A, row, term(i, -1, "phi", 1.0) + term(i + 1, 0, "phi", -1.0))                 # [phi] = 0
        row += 1
        if kink:                                                                          # tangential velocity
            put(A, row, term(i, -1, "phi1", uf) + term(i + 1, 0, "phi1", -uf) + term(i, -1, "phi", -ds))
            put(B, row, term(i, -1, "phi1", 1.0) + term(i + 1, 0, "phi1", -1.0))
        else:
            put(A, row, term(i, -1, "phi1", 1.0) + term(i + 1, 0, "phi1", -1.0))
        row += 1
        put(A, row, term(i, -1, "phi2", a.mu) + term(i, -1, "phi", a.mu * k2)             # tangential stress
            + term(i + 1, 0, "phi2", -b.mu) + term(i + 1, 0, "phi", -b.mu * k2))
        row += 1
        dr = a.rho - b.rho                         # normal stress, less i k rho_b times the tangential velocity's row
        if abs(dr) <= RHO_TOL * max(a.rho, b.rho):
            dr = 0.0
        put(A, row, term(i, -1, "phi1", ik * dr * uf + 3.0 * k2 * a.mu) + term(i, -1, "phi", -ik * dr * a.s)
            + term(i, -1, "phi3", -a.mu) + term(i + 1, 0, "phi3", b.mu) + term(i + 1, 0, "phi1", -3.0 * k2 * b.mu))
        if dr != 0.0:
            put(B, row, term(i, -1, "phi1", ik * dr))
        row += 1
    last = len(regions) - 1
    put(A, row, term(last, -1, "phi", 1.0))
    put(A, row + 1, term(last, -1, "phi1", 1.0))
    row += 2
    assert row == n, (row, n)
    scale = np.maximum(np.abs(A).max(axis=1), np.abs(B).max(axis=1))
    A /= scale[:, None]
    B /= scale[:, None]
    return A, B, offsets


class Spectrum(NamedTuple):
    """The finite eigenvalues c (m/s) of a column at one wavenumber, sorted by falling Im c, spurious ones removed;
    `tail`, the resolution measure of each of the first CANDIDATES (NaN beyond); `size`, the number of unknowns."""
    c: np.ndarray
    tail: np.ndarray
    size: int


def _solve(A, B, shift, V, method="auto"):
    """(c, vectors): the finite eigenvalues of A x = c B x and a function giving the eigenvector of the j-th.

    "dense": shift and invert with a dense eigenproblem of the size of R (module docstring), every eigenvalue.
    "arpack": the same operator (A - s B)^-1 B by ARPACK (sparse LU, Arnoldi), the nev eigenvalues of largest |mu|,
    i.e. nearest the shift, nev = ARPACK_NEV, kept only if the farthest found lies at least ARPACK_REACH speed scales V
    from the shift, so that every eigenvalue within that distance -- every one with Im c above about -0.1 V and Re c
    within 1.5 V of the middle of the speed range -- is among them; otherwise (many weakly damped modes, as in a deep
    sheared layer at high Reynolds number) the dense solve takes over. "auto": dense up to DENSE_MAX unknowns, arpack above. "qz": scipy's QZ on the whole pencil
    (slow; the tests' cross-check)."""
    n = A.shape[0]
    if method == "auto":
        method = "dense" if n <= DENSE_MAX else "arpack"
    if method == "qz":
        (alpha, beta), vr = scipy.linalg.eig(A, B, right=True, homogeneous_eigvals=True)
        finite = np.abs(beta) > 1e3 * np.finfo(float).eps * np.abs(alpha)
        c = np.where(finite, alpha / np.where(finite, beta, 1.0), np.nan)
        return c, lambda j: vr[:, j]
    lu = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(A - shift * B))
    R = np.flatnonzero(np.abs(B).max(axis=1) > 0.0)
    if method == "arpack":
        Bs = scipy.sparse.csr_matrix(B)
        op = scipy.sparse.linalg.LinearOperator((n, n), matvec=lambda x: lu.solve(Bs @ x), dtype=complex)
        if 3 * ARPACK_NEV <= R.size:
            mu, vecs = scipy.sparse.linalg.eigs(op, k=ARPACK_NEV, which="LM", v0=np.ones(n, dtype=complex))
            c = shift + 1.0 / mu
            if np.abs(c - shift).max() >= ARPACK_REACH * V:
                return c, lambda j: vecs[:, j]
    E = np.zeros((n, R.size), dtype=complex)
    E[R, np.arange(R.size)] = 1.0
    Wm = lu.solve(E)
    K = B[R] @ Wm
    mu = np.linalg.eigvals(K)
    finite = np.abs(mu) > 1e3 * np.finfo(float).eps * np.abs(mu).max()
    c = np.where(finite, shift + 1.0 / np.where(finite, mu, 1.0), np.nan)

    def vector(j):
        """The eigenvector of K at mu_j by inverse iteration, mapped back: x = W v."""
        lk = scipy.linalg.lu_factor(K - mu[j] * (1.0 + 1e-10) * np.eye(R.size), check_finite=False)
        v = np.ones(R.size, dtype=complex)
        for _ in range(3):
            v = scipy.linalg.lu_solve(lk, v, check_finite=False)
            v /= np.linalg.norm(v)
        return Wm @ v
    return c, vector


def _tail(x, regions, offsets):
    """The eigenfunction's resolution: the largest, over the regions, of the last two Chebyshev coefficients' moduli of
    phi, over the largest coefficient of phi in any region."""
    tails, peak = [], 0.0
    for i, r in enumerate(regions):
        a = _cheb(r.N)[2] @ x[offsets[i]:offsets[i] + r.N + 1]
        peak = max(peak, np.abs(a).max())
        tails.append(np.abs(a[-2:]).sum())
    return max(tails) / peak if peak > 0.0 else np.inf


def _setup(k, y, U, mu, rho, sigma, top, u_faces, resolution):
    prof = base_profile(y, U, mu, rho, top, u_faces)
    V = _speed_scale(prof, k, sigma)
    regions = _degrees(_regions(prof), k, V, resolution)
    speeds = np.concatenate([prof.u_faces, prof.U])
    shift = 0.5 * (speeds.max() + speeds.min()) + 1j * SHIFT_UP * V
    return prof, V, regions, shift


def spectrum(k, y, U, mu, rho, sigma, top="free", u_faces=None, resolution=1, method="auto"):
    """The column's eigenvalues at the wavenumber k (1/m): `Spectrum`. Arguments as for `growth_rate`."""
    prof, V, regions, shift = _setup(k, y, U, mu, rho, sigma, top, u_faces, resolution)
    A, B, offsets = _assemble(regions, k, sigma, top)
    c, vector = _solve(A, B, shift, V, method)
    keep = np.flatnonzero(np.isfinite(c) & (np.where(np.isfinite(c), c.imag, np.inf) <= GROWTH_BOUND * V))
    keep = keep[np.argsort(-c[keep].imag, kind="stable")]
    tail = np.full(keep.size, np.nan)
    for j in range(min(CANDIDATES, keep.size)):
        tail[j] = _tail(vector(keep[j]), regions, offsets)
        if tail[j] <= TAIL_TOL:
            break
    return Spectrum(c=c[keep], tail=tail, size=A.shape[0])


def growth_rate(k, y, U, mu, rho, sigma, top="free", u_faces=None, resolution=1, method="auto"):
    """The most unstable Orr-Sommerfeld eigenvalue c (complex, m/s) of a column at the wavenumber k (1/m); its growth
    rate is k Im c (1/s).

    y: the cell faces, m, downward (L + 1); U: the base speed at the cell centres, m/s (L); mu (Pa s) and rho (kg/m3)
    per cell; sigma: the surface tension, N/m (used at a free top). top: "free" (a free surface with tension) or "wall"
    (a no-slip wall moving at the top face's speed). The bottom is a no-slip wall moving at the bottom face's speed.
    u_faces: the face speeds (L + 1), m/s; by default the skins' stress continuity at the faces, the surface speed
    extrapolated from the top two centres and a wall at rest (`base_profile`). resolution: a factor on every region's
    Chebyshev degree; 2 is the doubling check. method: "auto" (default), "dense", "arpack" or "qz" (`_solve`).

    The eigenvalue taken is the one with the largest Im c whose eigenfunction is resolved (module docstring); if none of
    the first CANDIDATES is, the first is returned with a warning."""
    spec = spectrum(k, y, U, mu, rho, sigma, top, u_faces, resolution, method)
    if spec.c.size == 0:
        return complex(np.nan, np.nan)
    ok = np.flatnonzero(spec.tail <= TAIL_TOL)
    if ok.size:
        return complex(spec.c[ok[0]])
    warnings.warn("no resolved eigenfunction among the first {} eigenvalues at k = {:.6g} (tails {}); raise the "
                  "resolution".format(CANDIDATES, k, np.array2string(spec.tail[:CANDIDATES], precision=2)))
    return complex(spec.c[0])


class Fastest(NamedTuple):
    """The fastest-growing wave over a wavenumber range: k (1/m), growth k Im c (1/s), c (m/s), and `edge`, True where
    the fastest lies at an end of the range (the growth may rise beyond it)."""
    k: float
    growth: float
    c: complex
    edge: bool


def fastest(y, U, mu, rho, sigma, k_range, top="free", u_faces=None, resolution=1, n_grid=13, n_refine=10,
            method="auto"):
    """The fastest-growing wave over k_range = (k_lo, k_hi) 1/m: n_grid log-spaced wavenumbers, then n_refine
    golden-section steps in log k between the best point's neighbours. Other arguments as for `growth_rate`."""
    lo, hi = np.log(k_range[0]), np.log(k_range[1])
    xs = np.linspace(lo, hi, n_grid)
    cache = {}

    def g(x):
        if x not in cache:
            c = growth_rate(np.exp(x), y, U, mu, rho, sigma, top, u_faces, resolution, method)
            cache[x] = (np.exp(x) * c.imag if np.isfinite(c) else -np.inf, c)
        return cache[x][0]

    vals = np.array([g(x) for x in xs])
    b = int(np.argmax(vals))
    a_, c_ = xs[max(b - 1, 0)], xs[min(b + 1, n_grid - 1)]
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = c_ - gr * (c_ - a_), a_ + gr * (c_ - a_)
    for _ in range(n_refine):
        if g(x1) >= g(x2):
            c_, x2 = x2, x1
            x1 = c_ - gr * (c_ - a_)
        else:
            a_, x1 = x1, x2
            x2 = a_ + gr * (c_ - a_)
    best = max(cache, key=lambda x: cache[x][0])
    return Fastest(k=float(np.exp(best)), growth=float(cache[best][0]), c=complex(cache[best][1]),
                   edge=bool(best <= lo + 1e-12 or best >= hi - 1e-12))


def column_cells(cols, i):
    """(faces y [m, downward], U [m/s], mu [Pa s], rho [kg/m3]) of column i of a `melt_profile.LayerColumns`: its
    contiguous flowing cells from the top. The bottom face of the last flowing cell is the no-slip wall, as in the skins'
    momentum step (a flowing bottom cell drags on the rigid metal over its own half-thickness)."""
    f = np.asarray(cols.flowing[i], dtype=bool)
    n = int(np.argmin(f)) if not f.all() else f.size
    dz = np.asarray(cols.dz[i, :n], dtype=float)
    y = np.concatenate([[0.0], np.cumsum(dz)])
    return y, np.asarray(cols.u[i, :n], dtype=float), np.asarray(cols.mu[i, :n], dtype=float), \
        np.asarray(cols.rho[i, :n], dtype=float)
