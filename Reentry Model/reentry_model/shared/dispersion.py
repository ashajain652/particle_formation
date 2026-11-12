"""Girin's (2017, A&A 606, A63, Eq. 1) dispersion relation of the gradient instability of a sheared liquid layer,
solved numerically for the fastest-growing disturbance (spec Step 3 section 9).

    (Omega - Delta) [(Omega - Delta)(Omega + K Delta) + (1 - K) Delta] = Delta^3 We_s^-1 (Omega - K Delta),
    K = [1 - exp(-2 Delta)] / (2 Delta),   Omega = omega delta_m / V_s,   Delta = 2 pi delta_m / lambda,
    We_s = rho_m V_s^2 delta_m / Sigma.

For every We_s the cubic in Omega is solved on a grid of Delta (companion-matrix eigenvalues, vectorised) and the
Delta with the largest Im Omega is the fastest mode: Delta_f (wavelength 2 pi delta_m / Delta_f), Im Omega_f (growth
rate Im Omega_f V_s / delta_m) and Re Omega_f. Measured 2026-09-20: instability appears between We_s 3.00 and 3.08
(Girin: 3.08); Delta_f -> 1.226, Im Omega_f -> 0.247 at We_s 1e4 (Girin's Fig. 5 asymptotes 1.225, 0.24);
Re/Im = 1.68 at We_s 20 and 1.56 at 1e4 (Girin: 1.5); at We_s 4.62: 0.58 / 0.049; 10: 1.34 / 0.18; 20: 1.29 / 0.22.
The table is cached in the package data (`girin_dispersion.json`) and interpolated in log We_s."""
import json
import os
from dataclasses import dataclass, fields

import numpy as np
import scipy.linalg

from .. import DATA_DIR

TABLE_PATH = os.path.join(DATA_DIR, "girin_dispersion.json")
WE_CRITICAL_THEORY = 3.08          # Girin's critical surface Weber number; the model's default (fact 161)
WE_CRITICAL_PRACTICAL = 4.62       # the value he recommends in practice (GI >= 0.6) for mm meteoroids; girin_case only
DELTA_GRID = np.linspace(0.02, 5.0, 2500)
WE_GRID = np.logspace(np.log10(3.0), 4.0, 200)
UNSTABLE_TOL = 1e-6                # Im Omega above which a root counts as unstable


def cubic_coefficients(delta, we):
    """Coefficients (a2, a1, a0) of Omega^3 + a2 Omega^2 + a1 Omega + a0 = 0, Girin's Eq. (1) expanded."""
    D = np.asarray(delta, dtype=float)
    K = (1.0 - np.exp(-2.0 * D)) / (2.0 * D)
    a2 = (K - 2.0) * D
    a1 = (1.0 - 2.0 * K) * D ** 2 + (1.0 - K) * D - D ** 3 / we
    a0 = K * D ** 3 - (1.0 - K) * D ** 2 + K * D ** 4 / we
    return a2, a1, a0


def growth_rates(we, delta=DELTA_GRID):
    """(Im Omega, Re Omega) of the most unstable root of Eq. (1) at each Delta for one We_s."""
    a2, a1, a0 = cubic_coefficients(delta, we)
    C = np.zeros((len(delta), 3, 3), dtype=complex)
    C[:, 0, 1] = 1.0
    C[:, 1, 2] = 1.0
    C[:, 2, 0], C[:, 2, 1], C[:, 2, 2] = -a0, -a1, -a2
    roots = np.linalg.eigvals(C)
    k = np.argmax(roots.imag, axis=1)
    best = roots[np.arange(len(delta)), k]
    return best.imag, best.real


def fastest_mode(we, delta=DELTA_GRID):
    """(Delta_f, Im Omega_f, Re Omega_f) of the fastest-growing disturbance; (nan, 0, 0) when nothing is unstable."""
    im, re = growth_rates(we, delta)
    j = int(np.argmax(im))
    if im[j] <= UNSTABLE_TOL:
        return float("nan"), 0.0, 0.0
    return float(delta[j]), float(im[j]), float(re[j])


def build_table(we_grid=WE_GRID, delta=DELTA_GRID):
    rows = np.array([fastest_mode(we, delta) for we in we_grid])
    return {"we": list(map(float, we_grid)), "delta_f": list(map(float, rows[:, 0])),
            "im_omega_f": list(map(float, rows[:, 1])), "re_omega_f": list(map(float, rows[:, 2])),
            "delta_grid": [float(delta[0]), float(delta[-1]), int(len(delta))]}


def write_table(path=TABLE_PATH):
    table = build_table()
    with open(path, "w") as fh:
        json.dump(table, fh)
    return table


class DispersionTable:
    """Delta_f(We_s), Im Omega_f(We_s), Re Omega_f(We_s) interpolated in log We_s from the cached table."""

    def __init__(self, path=TABLE_PATH):
        if not os.path.isfile(path):
            write_table(path)
        with open(path) as fh:
            t = json.load(fh)
        self.we = np.array(t["we"])
        self.delta_f, self.im_omega_f, self.re_omega_f = (np.array(t[k]) for k in ("delta_f", "im_omega_f", "re_omega_f"))
        unstable = self.im_omega_f > 0.0
        self.we_onset = float(self.we[unstable][0])              # first tabulated We_s with an unstable mode
        self._log_we = np.log(self.we[unstable])
        self._delta, self._im, self._re = self.delta_f[unstable], self.im_omega_f[unstable], self.re_omega_f[unstable]

    def __call__(self, we):
        """Arrays (Delta_f, Im Omega_f, Re Omega_f) for We_s (scalar or array); nan/0 below the onset."""
        we = np.asarray(we, dtype=float)
        x = np.log(np.maximum(we, self.we_onset))
        d, im, re = (np.interp(x, self._log_we, v) for v in (self._delta, self._im, self._re))
        stable = we < self.we_onset
        return np.where(stable, np.nan, d), np.where(stable, 0.0, im), np.where(stable, 0.0, re)


# ---------------------------------------------------------------------------------------------------------------------
# The broken-line solver (kinked-profile dispersion design, docs/superpowers/specs/2026-11-11-kinked-profile-dispersion-
# design.md, sections 3 and 6.2): Girin's gradient instability solved on any piecewise-linear speed profile with a free
# surface on top and a wall (or an unbounded uniform bottom) below. The table above stays the `table` path's.

DISPERSION_NAMES = ("table", "profile")
DEFAULT_DISPERSION = "profile"     # the command line's default on spraying skin runs (the spec's decision 9)
UNSTABLE_REL = 1e-6                # an eigenvalue is unstable where Im c exceeds this times the profile's largest |U|
K_SPAN = (0.02, 5.0)               # the fastest-wave search runs from K_SPAN[0] / D to K_SPAN[1] / s_min
K_PER_DECADE = 24                  # search-grid points per decade of wavenumber
REFINE_ROUNDS = 3                  # refinement rounds about the best grid point
REFINE_POINTS = 9                  # log-spaced points per round, between the best point's two neighbours
PENCIL_COND_MAX = 1e8              # B's condition number above which a pencil goes to scipy's QZ, one at a time
RELEASE_SAMPLES = 129              # samples per segment in each of the release depth's three sub-grids
_KH_FLAT = 40.0                    # k times a thickness beyond which exp(-k h), at most 4e-18, is dropped against 1
_INFINITE_BETA = 1e3 * np.finfo(float).eps   # |beta| / |alpha| at or below which a QZ eigenvalue counts as infinite


@dataclass
class KinkedProfiles:
    """n piecewise-linear speed profiles, each with a free surface on top and a wall, or an unbounded uniform bottom.

    Depths are in m, measured down from the surface; speeds in m/s; densities in kg/m3. Column i has count[i]
    breakpoints at the increasing depths depth[i, :count[i]], where the speed is u_break[i, :count[i]]; its surface
    moves at u_top[i]. Segment j, for j from 0 to count[i], runs from breakpoint j - 1 to breakpoint j (from the surface
    for j = 0, and from the last breakpoint to the wall for j = count[i]) and has density rho[i, j]. The wall lies at
    depth wall[i] and moves at u_wall[i]; wall = inf means an unbounded bottom, which must be uniform (u_wall equal to
    the last breakpoint's speed). Entries beyond a column's count are NaN padding.

    The dataclass accepts any arrays, including unresolved rows (count 0 with NaN depths, wall and speeds), which a
    profile fit produces; the solver functions are meant for valid rows, which `one` checks."""
    count: np.ndarray       # (n,) int: breakpoints per column
    depth: np.ndarray       # (n, B) m
    u_top: np.ndarray       # (n,) m/s
    u_break: np.ndarray     # (n, B) m/s
    wall: np.ndarray        # (n,) m, inf for an unbounded bottom
    u_wall: np.ndarray      # (n,) m/s
    rho: np.ndarray         # (n, B + 1) kg/m3, one per segment

    def __post_init__(self):
        self.count = np.asarray(self.count, dtype=int)
        for name in ("depth", "u_top", "u_break", "wall", "u_wall", "rho"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))

    def rows(self, idx):
        """The profiles idx (an index, an index array or a boolean mask), as a new KinkedProfiles."""
        if np.ndim(idx) == 0:
            idx = [idx]
        return KinkedProfiles(**{f.name: getattr(self, f.name)[idx] for f in fields(self)})

    def segments(self):
        """(top depth, bottom depth, top speed, bottom speed) of every segment, each (n, B + 1); the last segment's
        bottom is the wall (inf for an unbounded bottom); NaN beyond each column's count."""
        n, B = self.depth.shape
        j = np.arange(B + 1)[None, :]
        cnt = self.count[:, None]
        top = np.concatenate([np.zeros((n, 1)), self.depth], axis=1)
        u_hi = np.concatenate([self.u_top[:, None], self.u_break], axis=1)
        bot = np.where(j == cnt, self.wall[:, None], np.concatenate([self.depth, np.full((n, 1), np.nan)], axis=1))
        u_lo = np.where(j == cnt, self.u_wall[:, None], np.concatenate([self.u_break, np.full((n, 1), np.nan)], axis=1))
        return tuple(np.where(j > cnt, np.nan, a) for a in (top, bot, u_hi, u_lo))

    def slopes(self):
        """dU/dy per segment, (n, B + 1) in 1/s, with y pointing up, so a speed falling with depth has a positive
        slope; 0 in an unbounded bottom and NaN beyond each column's count."""
        top, bot, u_hi, u_lo = self.segments()
        return (u_hi - u_lo) / (bot - top)

    @classmethod
    def one(cls, depths, speeds, wall=np.inf, rho=1.0):
        """One profile: breakpoint depths (m), speeds = [u_top, *u_break, u_wall] (m/s), the wall's depth (m, inf for
        an unbounded uniform bottom) and the density per segment (kg/m3, a scalar or one per segment). Raises
        ValueError if the depths do not increase from the surface, if the wall is not below the last breakpoint, if an
        unbounded bottom is not uniform, or if a profile without breakpoints has no wall."""
        depths = np.asarray(depths, dtype=float).reshape(-1)
        speeds = np.asarray(speeds, dtype=float).reshape(-1)
        n_break = depths.size
        if speeds.size != n_break + 2:
            raise ValueError("speeds must be [u_top, *u_break, u_wall]: %d values for %d breakpoints"
                             % (speeds.size, n_break))
        rho = np.array(np.broadcast_to(np.asarray(rho, dtype=float), (n_break + 1,)))
        wall = float(wall)
        if not (np.isfinite(depths).all() and np.isfinite(speeds).all() and np.isfinite(rho).all()):
            raise ValueError("depths, speeds and densities must be finite")
        if (rho <= 0.0).any():
            raise ValueError("densities must be positive")
        if n_break and (depths[0] <= 0.0 or (np.diff(depths) <= 0.0).any()):
            raise ValueError("the breakpoint depths must increase from the surface")
        if np.isnan(wall) or wall <= (depths[-1] if n_break else 0.0):
            raise ValueError("the wall must lie below the deepest breakpoint")
        if np.isinf(wall) and n_break == 0:
            raise ValueError("a profile without breakpoints needs a wall")
        if np.isinf(wall) and speeds[-1] != speeds[-2]:
            raise ValueError("an unbounded bottom must be uniform: u_wall must equal the last breakpoint's speed")
        return cls(count=np.array([n_break]), depth=depths[None, :], u_top=speeds[:1], u_break=speeds[None, 1:-1],
                   wall=np.array([wall]), u_wall=speeds[-1:], rho=rho[None, :])


@dataclass
class Waves:
    """The fastest-growing wave of each profile (`fastest_waves`)."""
    k: np.ndarray               # (n,) wavenumber, 1/m; NaN where stable
    growth: np.ndarray          # (n,) growth rate Im omega_f = k Im c, 1/s; 0 where stable
    c: np.ndarray               # (n,) complex wave speed, m/s; NaN where stable
    release_depth: np.ndarray   # (n,) m; NaN where stable


def _speed_scale(p):
    """Each profile's largest |U| over its surface, breakpoints and wall (m/s), or 1 m/s for a profile at rest."""
    n, B = p.depth.shape
    valid = np.arange(B)[None, :] < p.count[:, None]
    speeds = np.concatenate([p.u_top[:, None], np.where(valid, p.u_break, 0.0), p.u_wall[:, None]], axis=1)
    V = np.max(np.abs(speeds), axis=1)
    return np.where(V > 0.0, V, 1.0)


def _layout(p):
    """The per-segment arrays the pencil needs, each (n, B + 1): thickness (m, inf for an unbounded bottom), the speed
    at each segment's top (m/s), the slope (1/s) and the density (kg/m3); with the speed scale (n,)."""
    top, bot, u_hi, u_lo = p.segments()
    return dict(top=top, bot=bot, h=bot - top, U=u_hi, S=(u_hi - u_lo) / (bot - top), rho=p.rho, V=_speed_scale(p))


def _pencil(lay, rows, N, k, sigma):
    """The pencil (A, B), each (r, m, N + 2, N + 2), of the profiles `rows` (all with N breakpoints) at wavenumbers
    k (r, m) and surface tension sigma (r,), in scaled variables: lengths times k, speeds over the profile's speed
    scale V, densities over the top segment's.

    In segment j the stream function is P_j exp(-k s) + Q_j exp(-k (h_j - s)), s the depth below the segment's top
    and h_j its thickness: each exponential is 1 at its own end of the segment and decays across it, so no entry grows
    with k h. Continuity of phi at each breakpoint and phi = 0 at the wall (Q_N = -exp(-k h_N) P_N; Q_N = 0 for an
    unbounded bottom, where exp(-k h_N) = 0) are free of c and are eliminated, leaving the unknowns x = (eta, P_0, ...,
    P_N) and the c-linear conditions: kinematic and normal stress at the surface, and pressure at each breakpoint.
    Returns (A, B, aux) with aux holding e = exp(-k h) (r, m, N + 1) and Phi (r, m, N + 2, N + 1), the node values of
    phi (surface, breakpoints, wall) as linear functions of P."""
    r, m = k.shape
    h, U, S, rho = (lay[name][rows, :N + 1] for name in ("h", "U", "S", "rho"))
    V = lay["V"][rows]
    kk = k[:, :, None]
    kh = kk * h[:, None, :]
    khc = np.minimum(kh, _KH_FLAT)
    e = np.where(kh < _KH_FLAT, np.exp(-khc), 0.0)
    om = -np.expm1(-2.0 * khc)                              # 1 - e^2, accurate for thin segments
    Ut = U / V[:, None]
    St = S[:, None, :] / (V[:, None, None] * kk)
    rt = rho / rho[:, :1]
    sg = sigma[:, None] * k / (rho[:, :1] * V[:, None] ** 2)
    g = np.sqrt(1.0 + sg)                                   # balances eta against phi when capillarity dominates
    # Node values: phi_j = om_j P_j + e_j phi_{j+1}, from phi_{N+1} = 0 at the wall (or far below) upward.
    Phi = np.zeros((r, m, N + 2, N + 1))
    for j in range(N, -1, -1):
        Phi[:, :, j, :] = e[:, :, j, None] * Phi[:, :, j + 1, :]
        Phi[:, :, j, j] += om[:, :, j]
    eye = np.eye(N + 1)
    Dt = (1.0 + e[..., None] ** 2) * eye - e[..., None] * Phi[:, :, 1:, :]   # phi'/k at each segment's top
    Db = 2.0 * e[..., None] * eye - Phi[:, :, 1:, :]                         # phi'/k at each segment's bottom
    A = np.zeros((r, m, N + 2, N + 2))
    B = np.zeros_like(A)
    # Kinematic: (U_0 - c) eta = -phi(0), with eta scaled by g.
    A[:, :, 0, 0] = Ut[:, None, 0]
    A[:, :, 0, 1:] = g[..., None] * Phi[:, :, 0, :]
    B[:, :, 0, 0] = 1.0
    # Normal stress: rho_0 [-(U_0 - c) phi'(0) + S_0 phi(0)] = sigma k^2 eta.
    A[:, :, 1, 0] = -sg / g
    A[:, :, 1, 1:] = -Ut[:, None, 0, None] * Dt[:, :, 0, :] + St[:, :, 0, None] * Phi[:, :, 0, :]
    B[:, :, 1, 1:] = -Dt[:, :, 0, :]
    # Pressure at breakpoint j: (U_j - c)(rho_above phi'_above - rho_below phi'_below) = (rho_above S_above -
    # rho_below S_below) phi_j.
    for j in range(1, N + 1):
        J = rt[:, None, j - 1, None] * Db[:, :, j - 1, :] - rt[:, None, j, None] * Dt[:, :, j, :]
        dS = rt[:, None, j - 1] * St[:, :, j - 1] - rt[:, None, j] * St[:, :, j]
        A[:, :, 1 + j, 1:] = Ut[:, None, j, None] * J - dS[..., None] * Phi[:, :, j, :]
        B[:, :, 1 + j, 1:] = J
    return A, B, dict(e=e, Phi=Phi)


def _qz(A, B):
    """The finite eigenvalues of A x = c B x by scipy's QZ (NaN where infinite, or if QZ fails)."""
    out = np.full(len(A), complex(np.nan, np.nan))
    sa, sb = np.linalg.norm(A), np.linalg.norm(B)
    try:
        alpha, beta = scipy.linalg.eig(A / sa, B / sb, right=False, homogeneous_eigvals=True)
    except (np.linalg.LinAlgError, ValueError):
        return out
    finite = np.abs(beta) > _INFINITE_BETA * np.abs(alpha)
    out[finite] = alpha[finite] / beta[finite] * (sa / sb)
    return out


def _pencil_eigenvalues(A, B):
    """The finite eigenvalues of a stack of pencils (s, q, q): batched through B^-1 A where B's 1-norm condition
    number is below PENCIL_COND_MAX, by QZ one at a time otherwise; NaN where an eigenvalue is infinite or a pencil is
    not finite."""
    s, q, _ = A.shape
    out = np.full((s, q), complex(np.nan, np.nan))
    ok = np.flatnonzero(np.isfinite(A).all(axis=(1, 2)) & np.isfinite(B).all(axis=(1, 2)))
    if ok.size == 0:
        return out
    slow = ok
    try:
        B_inv = np.linalg.inv(B[ok])
    except np.linalg.LinAlgError:
        B_inv = None
    if B_inv is not None:
        norm_b, norm_inv = (np.abs(X).sum(axis=-2).max(axis=-1) for X in (B[ok], B_inv))
        well = np.isfinite(norm_inv) & (norm_inv < PENCIL_COND_MAX / norm_b)
        slow = ok[~well]
        if well.any():
            try:
                out[ok[well]] = np.linalg.eigvals(B_inv[well] @ A[ok[well]])
            except np.linalg.LinAlgError:
                slow = ok
    for i in slow:
        out[i] = _qz(A[i], B[i])
    return out


def _as_wavenumbers(k, n):
    k = np.asarray(k, dtype=float)
    if k.ndim == 0:
        k = k[None]
    return np.broadcast_to(k, (n, k.shape[-1])) if k.ndim == 1 else k


def broken_line_eigenvalues(profiles, k, sigma):
    """Every finite complex wave speed c (m/s) of the inviscid broken-line problem, (n, m, B + 2), sorted by falling
    Im c and padded with NaN (a column with N breakpoints has N + 2). k is the wavenumber in 1/m, (m,) or (n, m) (a
    scalar is one wavenumber); sigma is the surface tension in N/m, a scalar or (n,). The growth rate of a wave is
    k Im c.

    The problem (spec section 3). y points up, the free surface is at y = 0 and the wall at depth H. A disturbance
    stream function phi(y) exp[ik(x - ct)] obeys the Rayleigh equation, and since U'' = 0 inside each straight segment,
    phi'' = k^2 phi there. The conditions:
      - wall: phi(-H) = 0; an unbounded bottom keeps only the exponential that decays downward;
      - each breakpoint y_j: phi continuous, and the pressure continuous,
        rho_above [(U_j - c) phi'_above - S_above phi_j] = rho_below [(U_j - c) phi'_below - S_below phi_j],
        with S the segment's slope dU/dy, which for equal densities is Girin's (U - c)[phi'] = [U'] phi;
      - free surface, displacement eta, gas pressure neglected: kinematic, (U_0 - c) eta = -phi(0); normal stress
        with surface tension sigma, rho_0 [-(U_0 - c) phi'(0) + S_0 phi(0)] = sigma k^2 eta.
    These form a generalised linear eigenvalue problem (A - c B) x = 0 with N + 2 finite eigenvalues. The model is
    inviscid and neglects the gas, deceleration and any tension at the breakpoints.

    Sources: the side session's derivation (research_notes/Kinked profile film dispersion/fulltext_girin_papers.md,
    sections 1b and 4; its report `reports/Kinked profile film dispersion 2.md`), which found the extension to slurry
    and a wall agrees with Ghosh et al. (2017, Meccanica) to between 1e-16 and 1e-12. With one breakpoint and no wall
    (Girin's profile) the problem reproduces Girin's (2017, A&A 606, A63) Eq. 1 to 1e-10 (tested), and with slurry over
    a wall the notes' closed form to 1e-10 (tested); Ghosh's peak growth rates are reproduced to 0.003 (tested).
    The derivation is not published: derivation ours, checked against Girin Eq. 1 and Ghosh 2017.

    The method (spec section 8, the builder's): `_pencil` writes each segment's exponentials in a basis normalised
    inside that segment and eliminates the c-free conditions, leaving an (N + 2)-square pencil, solved batched by
    numpy where B is well conditioned and by scipy's QZ otherwise; infinite eigenvalues are discarded."""
    p = profiles
    n, B = p.depth.shape
    k = _as_wavenumbers(k, n)
    m = k.shape[1]
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), (n,))
    out = np.full((n, m, B + 2), complex(np.nan, np.nan))
    lay = _layout(p)
    for N in np.unique(p.count):
        rows = np.flatnonzero(p.count == N)
        A, Bm, _ = _pencil(lay, rows, N, k[rows], sigma[rows])
        c = _pencil_eigenvalues(A.reshape(-1, N + 2, N + 2), Bm.reshape(-1, N + 2, N + 2))
        c = c.reshape(len(rows), m, N + 2) * lay["V"][rows, None, None]
        order = np.argsort(np.where(np.isnan(c), np.inf, -c.imag), axis=-1, kind="stable")
        out[rows, :, :N + 2] = np.take_along_axis(c, order, axis=-1)
    return out


def broken_line_c(profiles, k, sigma):
    """The eigenvalue with the largest imaginary part, (n, m) complex m/s (NaN where none is finite); arguments as for
    `broken_line_eigenvalues`."""
    return broken_line_eigenvalues(profiles, k, sigma)[..., 0]


def _growth(p, k, sigma):
    """(growth k Im c in 1/s, c in m/s) of the most unstable eigenvalue at wavenumbers k (n, m); growth -inf where no
    eigenvalue is finite."""
    c = broken_line_c(p, k, sigma)
    g = k * c.imag
    return np.where(np.isfinite(g), g, -np.inf), c


def wave_structure(profiles, k, c, y):
    """The depth structure of the wave (k, c) of each profile: (phi complex (n, L), displacement (n, L)) at depths y in
    m, (L,) or (n, L). phi, the stream function (the vertical velocity is -ik phi), is normalised so phi(0) = 1; the
    displacement is the vertical displacement |phi / (U - c)| over its surface value. NaN where k or c is not finite
    and at depths outside the fluid. k is in 1/m and c in m/s, each (n,).

    Given c, the pressure conditions at the N breakpoints fix the N + 1 segment amplitudes up to a factor (the surface
    conditions are what make c an eigenvalue), so phi needs no surface tension."""
    p = profiles
    n, B = p.depth.shape
    k = np.broadcast_to(np.asarray(k, dtype=float), (n,))
    c = np.broadcast_to(np.asarray(c, dtype=complex), (n,))
    y = np.asarray(y, dtype=float)
    y = np.broadcast_to(y, (n, y.shape[-1]))
    L = y.shape[1]
    phi = np.full((n, L), complex(np.nan, np.nan))
    disp = np.full((n, L), np.nan)
    good = np.isfinite(k) & np.isfinite(c)
    good[good] = k[good] > 0.0
    lay = _layout(p)
    for N in np.unique(p.count[good]):
        rows = np.flatnonzero(good & (p.count == N))
        r = len(rows)
        A, Bm, aux = _pencil(lay, rows, N, k[rows, None], np.zeros(r))
        e, Phi = aux["e"][:, 0], aux["Phi"][:, 0]
        ct = c[rows] / lay["V"][rows]
        if N:
            M = A[:, 0, 2:, 1:] - ct[:, None, None] * Bm[:, 0, 2:, 1:]
            P = np.linalg.svd(M)[2][:, -1, :].conj()          # (r, N + 1): the null vector of the pressure rows
        else:
            P = np.ones((r, 1), dtype=complex)
        nodes = np.einsum("rjm,rm->rj", Phi, P)              # phi at the surface, each breakpoint and the wall
        Q = nodes[:, 1:] - e * P
        scale = nodes[:, 0]
        ok = np.abs(scale) > 0.0
        P[ok] /= scale[ok, None]
        Q[ok] /= scale[ok, None]
        P[~ok] = np.nan
        # The segment holding each depth: the number of breakpoints strictly above it.
        yr = y[rows]
        brk = np.where(np.arange(B)[None, :] < N, p.depth[rows], np.inf)
        seg = (brk[:, None, :] < yr[:, :, None]).sum(axis=-1)
        rr = np.arange(r)[:, None]
        top, bot = lay["top"][rows][rr, seg], lay["bot"][rows][rr, seg]
        inside = (yr >= 0.0) & (yr <= bot)
        kr = k[rows, None]
        x1 = np.where(inside, kr * (yr - top), 0.0)
        x2 = np.where(inside, kr * (bot - yr), np.inf)
        f1 = np.where(x1 < _KH_FLAT, np.exp(-np.minimum(x1, _KH_FLAT)), 0.0)
        f2 = np.where(x2 < _KH_FLAT, np.exp(-np.minimum(x2, _KH_FLAT)), 0.0)
        ph = P[rr, seg] * f1 + Q[rr, seg] * f2
        U = lay["U"][rows][rr, seg] - lay["S"][rows][rr, seg] * (yr - top)
        num = np.abs(ph) * np.abs(lay["U"][rows, :1] - c[rows, None])
        den = np.abs(U - c[rows, None])
        dsp = np.full(num.shape, np.inf)
        np.divide(num, den, out=dsp, where=den > 0.0)
        phi[rows] = np.where(inside, ph, np.nan)
        disp[rows] = np.where(inside, dsp, np.nan)
    return phi, disp


def _release_depth(p, k, c):
    """The release depth (m) of each profile's wave (k, c): the first depth below the displacement's peak where the
    displacement falls to 1/e of its surface value. Each segment is sampled on three sub-grids (across it, and over
    10/k below its top and above its bottom), an unbounded bottom down to 10/k below the last breakpoint; the crossing
    is interpolated linearly and capped at the wall (or at the deepest sample without one)."""
    n, B = p.depth.shape
    top, bot, _, _ = p.segments()
    seg = np.arange(B + 1)[None, :] <= p.count[:, None]
    reach = 10.0 / k[:, None]
    bot = np.where(seg & np.isinf(bot), top + reach, bot)
    last = bot[np.arange(n), p.count]
    top, bot = np.where(seg, top, last[:, None]), np.where(seg, bot, last[:, None])
    h = bot - top
    span = np.minimum(h, reach)
    t = np.linspace(0.0, 1.0, RELEASE_SAMPLES)
    y = np.concatenate([top[..., None] + h[..., None] * t, top[..., None] + span[..., None] * t,
                        bot[..., None] - span[..., None] * t], axis=-1)
    y = np.sort(y, axis=-1).reshape(n, -1)
    _, disp = wave_structure(p, k, c, y)
    thr = np.exp(-1.0)
    peak = np.argmax(disp, axis=1)
    below = (np.arange(y.shape[1])[None, :] >= peak[:, None]) & (disp <= thr)
    found = below.any(axis=1)
    depth = np.where(np.isfinite(p.wall), p.wall, y[:, -1])
    if found.any():
        f = np.flatnonzero(found)
        i = np.argmax(below[f], axis=1)
        y0, y1, d0, d1 = y[f, i - 1], y[f, i], disp[f, i - 1], disp[f, i]
        depth[f] = y0 + (d0 - thr) * (y1 - y0) / (d0 - d1)
    return depth


def fastest_waves(profiles, sigma):
    """The fastest-growing wave of each profile (`Waves`), at surface tension sigma (N/m, a scalar or (n,)).

    The growth k Im c is searched on a grid log-spaced at K_PER_DECADE points per decade from K_SPAN[0] / D to
    K_SPAN[1] / s_min, with D the wall's depth (the deepest breakpoint's without a wall) and s_min the thinnest
    segment; then REFINE_ROUNDS rounds of REFINE_POINTS log-spaced points between the best point's two neighbours;
    then the vertex of the parabola (in log k) through the best point and its neighbours, kept where it grows faster.
    A profile is stable where the fastest wave's Im c is at most UNSTABLE_REL times its largest |U|; then k, c and the
    release depth are NaN and the growth 0."""
    p = profiles
    n, B = p.depth.shape
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), (n,))
    waves = Waves(k=np.full(n, np.nan), growth=np.zeros(n), c=np.full(n, complex(np.nan, np.nan)),
                  release_depth=np.full(n, np.nan))
    top, bot, _, _ = p.segments()
    h = bot - top
    seg = np.arange(B + 1)[None, :] <= p.count[:, None]
    s_min = np.min(np.where(seg & np.isfinite(h), h, np.inf), axis=1)
    deepest = p.depth[np.arange(n), np.maximum(p.count - 1, 0)] if B else np.full(n, np.nan)
    D = np.where(np.isfinite(p.wall), p.wall, np.where(p.count > 0, deepest, np.nan))
    # Rows without a finite depth scale (unresolved rows, which callers normally filter out) stay stable.
    valid = np.isfinite(D) & np.isfinite(s_min)
    valid[valid] = (D[valid] > 0.0) & (s_min[valid] > 0.0)
    v = np.flatnonzero(valid)
    if v.size == 0:
        return waves
    q, sig = p.rows(v), sigma[v]
    lo, hi = np.log(K_SPAN[0] / D[v]), np.log(K_SPAN[1] / s_min[v])
    nv = len(v)
    rr = np.arange(nv)
    step = np.log(10.0) / K_PER_DECADE
    last = np.ceil((hi - lo) / step - 1e-9).astype(int)          # the last grid point, clipped to hi
    x = np.minimum(lo[:, None] + step * np.arange(last.max() + 1)[None, :], hi[:, None])   # ln k
    g, c = _growth(q, np.exp(x), sig)
    b = np.argmax(g, axis=1)
    xb, gb, cb = x[rr, b], g[rr, b], c[rr, b]
    xl, xr = x[rr, np.maximum(b - 1, 0)], x[rr, np.minimum(b + 1, last)]
    t = np.linspace(0.0, 1.0, REFINE_POINTS)
    for _ in range(REFINE_ROUNDS):
        xs = xl[:, None] + (xr - xl)[:, None] * t[None, :]
        gs, cs = _growth(q, np.exp(xs), sig)
        b = np.argmax(gs, axis=1)
        better = gs[rr, b] > gb
        xb, gb, cb = np.where(better, xs[rr, b], xb), np.where(better, gs[rr, b], gb), np.where(better, cs[rr, b], cb)
        xl, xr = xs[rr, np.maximum(b - 1, 0)], xs[rr, np.minimum(b + 1, REFINE_POINTS - 1)]
    # The parabola through the last round's best point and its neighbours, where they bracket a concave peak.
    inner = (b > 0) & (b < REFINE_POINTS - 1)
    gm, g0, gp = gs[rr, np.maximum(b - 1, 0)], gs[rr, b], gs[rr, np.minimum(b + 1, REFINE_POINTS - 1)]
    curv = np.where(inner & np.isfinite(gm) & np.isfinite(gp), gp - 2.0 * g0 + gm, 0.0)
    para = curv < 0.0
    if para.any():
        dx = (xr - xl) / 2.0
        xv = np.where(para, xs[rr, b] - dx * (gp - gm) / (2.0 * np.where(para, curv, -1.0)), xb)
        gv, cv = _growth(q, np.exp(xv[:, None]), sig)
        better = para & (gv[:, 0] > gb)
        xb, gb, cb = np.where(better, xv, xb), np.where(better, gv[:, 0], gb), np.where(better, cv[:, 0], cb)
    finite = np.isfinite(cb)
    unstable = finite & (np.where(finite, cb.imag, 0.0) > UNSTABLE_REL * _speed_scale(q))
    u = np.flatnonzero(unstable)
    if u.size:
        kf = np.exp(xb[u])
        waves.k[v[u]] = kf
        waves.growth[v[u]] = gb[u]
        waves.c[v[u]] = cb[u]
        waves.release_depth[v[u]] = _release_depth(q.rows(u), kf, cb[u])
    return waves
