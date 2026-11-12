"""The level set's skins (design section 20.6, level-set stage 3 Task 2): columns of thin cells standing on the surface
markers, here only as a container of cells -- their storage, their grading with depth, their curved cross-section and
the regrid that keeps every cell near its target thickness. Conduction (Task 3), the layered flow (Task 4), the markers
(Task 5) and runoff (Task 6) build on it.

Storage. n columns side by side in padded arrays: cell 0 is a column's top (at the gas surface phi), cells
0 .. count-1 are real and every padded cell holds zeros. Per cell: mass `m` (n, K) [kg], enthalpy `H` (n, K) [J] above
the shared T_REF (m times the shared Material's specific enthalpy) and velocity `U` (n, K, 3) [m/s]; per column: the top
area `area` (n,) [m2], A0, and phi's principal curvatures `kappa` (n, 2) [1/m], positive where the surface is convex.
Temperature, liquid fraction, density and phase follow from H/m through the material every time they are asked for;
phase labels are never stored (as the skyscrapers' skins).

Cross-section. A column narrows (convex) or widens (concave) with the depth z below its top:
A(z) = A0 (1 - k1 z)(1 - k2 z). A cell of mass m at density rho fills the span [z_top, z_bot] with
int A(z) dz = m / rho; the volume from the top to depth z is A0 (z - (k1 + k2) z^2 / 2 + k1 k2 z^3 / 3), a cubic, solved
for each cell's bottom by a bracketed Newton iteration from the flat value. A flat column (k1 = k2 = 0) takes the flat
thickness m / (rho A0) directly. A column whose volume would reach beyond 1/k of a positive curvature (a skin deeper than
its radius of curvature) has no positive cross-section there: `dz` raises ValueError rather than give negative areas.

Grading and regrid (the skyscrapers' option B' rule, matched by amendment, never imported). Each cell's target
thickness is cell + (growth - 1) max(z_top - fine_zone, 0), z_top the depth of its top: 10 um (production) through the
top 0.4 mm, then the continuous form of "each cell 10 % thicker than the one above". `regrid` merges a cell below
MERGE_BELOW of its target into a neighbour and splits one above SPLIT_ABOVE in two, exactly in mass, enthalpy and
momentum, for at most REGRID_PASSES rounds. Merges pair cells (in a run of small cells every other one takes the cell
below it; a small bottom cell joins the cell above), so a run of thin cells coarsens by halving rather than collapsing
into one cell to be split again.

Every physical property comes from the material passed in (the shared `reentry_model.shared.material.Material` or a
test material with its interface); this module's numbers are numerical settings only."""
from dataclasses import dataclass

import numpy as np

from reentry_model.shared.thermal import SIGMA_SB

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (design 20.6; skyscrapers 6.2)
MERGE_BELOW = 0.5          # a cell thinner than this times its target joins a neighbour (skyscrapers' rule)
SPLIT_ABOVE = 1.5          # a cell thicker than this times its target is split in two (skyscrapers' rule)
REGRID_PASSES = 64         # merge-or-split rounds of one regrid at most (skyscrapers' rule)
DEPTH_NEWTON_MAX = 100     # iterations of the cross-section's depth solve at most (it takes about five)
DEPTH_NEWTON_RTOL = 4.0 * np.finfo(float).eps   # the depth solve stops when a step is below this times the depth
                                                # or the volume's residual below this times the volume (round-off)
NEWTON_TOL = 1.0e-6        # K: a heat sub-step has converged when no cell moves by more than this (skyscrapers' value)
NEWTON_MAX = 40            # Newton iterations of a heat sub-step at most (skyscrapers' value)
LINE_SEARCH_MAX = 12       # halvings of a Newton step before it is taken anyway (skyscrapers' value)
DZ_FLOOR = 1.0e-12         # m: the conduction length of a vanishing cell
DK_STEP = 1.0e-3           # K: half-width of the central difference giving dk/dT for the Newton Jacobian
PROFILE_POINTS = 8         # Gauss-Legendre points of the enthalpy a draw averages along the volume's profile


@dataclass
class SkinSettings:
    cell: float = SKIN_PRESETS["production"]   # m, the fine cell
    fine_zone: float = 0.4e-3      # m below the top held in fine cells (design 20.6)
    growth: float = 1.10           # below the fine zone each cell this much thicker than the one above
    stub: float = 0.8e-3           # m of rigid metal a full-depth skin keeps below psi (`--skin-stub`, design 20.6)
    substep: float = 0.010         # s, the skins' heat sub-step (skyscrapers' spec 7.1)
    swell: bool = True             # each cell at its own phase's mixture density; False holds every cell at the solid's

    def __post_init__(self):
        if not (self.cell > 0.0 and self.fine_zone > 0.0 and self.stub > 0.0 and self.substep > 0.0):
            raise ValueError("skin cell, fine zone, stub and substep must be positive")
        if self.growth < 1.0:
            raise ValueError("skin growth must be at least 1")
        if 2.0 * self.cell > min(self.stub, self.fine_zone) * (1.0 + 1e-12):
            raise ValueError("a skin needs at least two fine cells: cell <= min(stub, fine_zone) / 2")


def _volume(z, s, p):
    """Volume per unit top area from the top to depth z [m], s = k1 + k2, p = k1 k2."""
    return z * (1.0 - z * (0.5 * s - z * p / 3.0))


def _depth_of_volume(v, k1, k2):
    """The depth z [m] holding the volume per unit top area v [m] under the cross-section (1 - k1 z)(1 - k2 z), for
    flat arrays of v, k1, k2. Raises ValueError where v exceeds what the cross-section holds before it closes."""
    s, p = k1 + k2, k1 * k2
    kmax = np.maximum(k1, k2)
    closes = kmax > 0.0
    z_close = np.divide(1.0, kmax, out=np.full_like(v, np.inf), where=closes)
    v_close = np.where(closes, _volume(np.where(closes, z_close, 0.0), s, p), np.inf)
    deep = v > v_close * (1.0 + 1e-12)
    if deep.any():
        raise ValueError("a skin deeper than its radius of curvature: {} cell(s) lie beyond 1/kappa, where the "
                         "cross-section A0 (1 - k1 z)(1 - k2 z) is no longer positive".format(int(deep.sum())))
    # bracket: a column that never closes widens with depth, so the flat value is the deepest it can be (both
    # curvatures <= 0); one that closes lies above 1/kappa_max
    lo, hi = np.zeros_like(v), np.where(closes, z_close, v)
    z = np.minimum(v, hi)
    for _ in range(DEPTH_NEWTON_MAX):
        f = _volume(z, s, p) - v
        lo, hi = np.where(f < 0.0, z, lo), np.where(f > 0.0, z, hi)
        d = 1.0 - s * z + p * z * z
        step = np.divide(f, d, out=np.full_like(z, np.nan), where=d > 0.0)
        z_new = z - step
        bad = ~np.isfinite(z_new) | (z_new < lo) | (z_new > hi)
        z_new = np.where(bad, 0.5 * (lo + hi), z_new)
        z_new = np.where(f == 0.0, z, z_new)
        done = (np.abs(z_new - z) <= DEPTH_NEWTON_RTOL * np.abs(z_new)) | (np.abs(f) <= DEPTH_NEWTON_RTOL * v)
        z = z_new
        if done.all():
            return z
    raise RuntimeError("the skin's cross-section depth solve did not converge in {} iterations".format(DEPTH_NEWTON_MAX))


@dataclass
class HeatBook:
    """What one heat sub-step (`Skins.conduct`) did, per column [J]: the convective heat taken in at the top, the heat
    radiated from it, the heat conducted into the volume through the base (the residual of the books, so they close
    exactly whatever the iteration error), the heat gained sideways (zero: no sideways conduction in the level set's
    skins yet) and, in a trial pass, the base heat's derivative with respect to the base temperature [J/K]."""
    E_top: np.ndarray
    E_rad: np.ndarray
    E_base: np.ndarray
    E_side: np.ndarray
    dE_base_dTb: np.ndarray = None


def _tridiagonal(lower, diag, upper, rhs):
    """Rows of tridiagonal systems (n, K) solved by the Thomas algorithm: row j reads
    lower[:, j] x[:, j-1] + diag[:, j] x[:, j] + upper[:, j] x[:, j+1] = rhs[:, j] (lower[:, 0], upper[:, -1] unused)."""
    n, K = diag.shape
    c, g, x = np.zeros((n, K)), np.zeros((n, K)), np.zeros((n, K))
    den = diag[:, 0]
    c[:, 0], g[:, 0] = upper[:, 0] / den, rhs[:, 0] / den
    for j in range(1, K):
        den = diag[:, j] - lower[:, j] * c[:, j - 1]
        c[:, j] = upper[:, j] / den
        g[:, j] = (rhs[:, j] - lower[:, j] * g[:, j - 1]) / den
    x[:, K - 1] = g[:, K - 1]
    for j in range(K - 2, -1, -1):
        x[:, j] = g[:, j] - c[:, j] * x[:, j + 1]
    return x


def _profile_points(T_b, g, s0, s1, cap):
    """The Gauss-Legendre points' temperatures [K] (n, PROFILE_POINTS) along T(s) = T_b - g s over the depths s0 to s1
    [m] below a base, each capped at `cap`, with the weights (summing to 1) and whether each point lies below the cap."""
    x, w = np.polynomial.legendre.leggauss(PROFILE_POINTS)
    s = 0.5 * (s0 + s1)[:, None] + 0.5 * (s1 - s0)[:, None] * x[None, :]
    T = T_b[:, None] - g[:, None] * s
    return np.minimum(T, cap), 0.5 * w, T < cap


class Skins:
    """n columns of cells (module docstring). `material` gives every property; `area` (n,) the top areas [m2];
    `kappa` (n, 2) phi's principal curvatures [1/m], positive convex. The columns start empty (`fill` or `install`)."""

    def __init__(self, material, area, kappa, settings=None):
        self.material = material
        self.settings = settings if settings is not None else SkinSettings()
        self.area = np.array(area, dtype=float).reshape(-1)
        self.kappa = np.array(kappa, dtype=float).reshape(len(self.area), 2)
        self.rho_s = float(material.rho)
        liquid = getattr(material, "liquid", None)
        self.rho_l = float(liquid.rho) if liquid is not None else self.rho_s   # a material without a liquid never swells
        self.dH = None             # (n, K) J/K, each cell's dH/dT_b while a trial pass is open (`conduct`), else None
        self.momentum_out = np.zeros((self.n, 3))   # kg m/s, what the last `remove_feed` took with the metal
        self.install(np.zeros((self.n, 1)), np.zeros((self.n, 1)), np.zeros(self.n, dtype=np.int64))

    # -- storage ----------------------------------------------------------------------------------------------------
    @property
    def n(self):
        return len(self.area)

    @property
    def valid(self):
        return np.arange(self.m.shape[1])[None, :] < self.count[:, None]

    def install(self, m, H, count, U=None):
        """Set the cells from padded arrays m, H (n, K), count (n,) and velocities U (n, K, 3) (at rest if None);
        cells beyond each count are zeroed."""
        m, H = np.array(m, dtype=float), np.array(H, dtype=float)
        count = np.array(count, dtype=np.int64).reshape(self.n)
        if m.shape != H.shape or m.shape[0] != self.n or np.any(count < 0) or np.any(count > m.shape[1]):
            raise ValueError("cell arrays must be (n, K) with 0 <= count <= K")
        valid = np.arange(m.shape[1])[None, :] < count[:, None]
        if np.any(m[valid] <= 0.0):
            raise ValueError("every real cell must hold a positive mass")
        U = np.zeros(m.shape + (3,)) if U is None else np.array(U, dtype=float).reshape(m.shape + (3,))
        self.m, self.H = np.where(valid, m, 0.0), np.where(valid, H, 0.0)
        self.U, self.count = np.where(valid[..., None], U, 0.0), count
        self._end_trial()
        self._trim()

    def fill(self, mass, T):
        """Every column filled with `mass` [kg] of metal at temperature T [K] (each per column or one for all), at rest,
        in equal cells of about the fine size by the flat thickness."""
        mass = np.broadcast_to(np.asarray(mass, dtype=float), (self.n,))
        T = np.broadcast_to(np.asarray(T, dtype=float), (self.n,))
        flat = mass / (self.density(T) * self.area)
        cells = np.maximum(np.rint(flat / self.settings.cell).astype(np.int64), 1)
        K = int(cells.max()) if self.n else 1
        valid = np.arange(K)[None, :] < cells[:, None]
        m = np.where(valid, (mass / cells)[:, None], 0.0)
        self.install(m, m * self.material.enthalpy(T)[:, None], cells)

    def _width(self, K):
        """Pad the arrays to at least K cells."""
        extra = K - self.m.shape[1]
        if extra > 0:
            self.m = np.pad(self.m, ((0, 0), (0, extra)))
            self.H = np.pad(self.H, ((0, 0), (0, extra)))
            self.U = np.pad(self.U, ((0, 0), (0, extra), (0, 0)))
            if self.dH is not None:
                self.dH = np.pad(self.dH, ((0, 0), (0, extra)))

    def _trim(self):
        """Drop padding beyond the longest column (keeping at least one cell of width)."""
        K = max(int(self.count.max()) if self.n else 0, 1)
        if K < self.m.shape[1]:
            self.m, self.H, self.U = self.m[:, :K], self.H[:, :K], self.U[:, :K]
            if self.dH is not None:
                self.dH = self.dH[:, :K]

    def totals(self):
        """The books per column: (mass [kg] (n,), enthalpy above T_REF [J] (n,), momentum [kg m/s] (n, 3)), summed
        over the real cells."""
        return self.m.sum(axis=1), self.H.sum(axis=1), (self.m[..., None] * self.U).sum(axis=1)

    # -- what follows from the enthalpy -----------------------------------------------------------------------------
    def T(self):
        """Cell temperatures [K] (n, K) from the specific enthalpies through the material; nan beyond each count."""
        valid = self.valid
        T = np.full(self.m.shape, np.nan)
        if valid.any():
            T[valid] = self.material.temperature_from_enthalpy(self.H[valid] / self.m[valid])
        return T

    def density(self, T):
        """Mixture density [kg/m3]: specific volumes add, 1 / ((1 - f_l)/rho_s + f_l/rho_l), f_l the material's liquid
        fraction; rho_s everywhere with `swell=False`."""
        T = np.asarray(T, dtype=float)
        if not self.settings.swell:
            return np.full(T.shape, self.rho_s)
        f = self.material.liquid_fraction(T)
        return 1.0 / ((1.0 - f) / self.rho_s + f / self.rho_l)

    def cross_section(self, z):
        """The cross-section [m2] at depths z [m] below each column's top, z (n, ...)."""
        z = np.asarray(z, dtype=float)
        shape = (self.n,) + (1,) * (z.ndim - 1)
        k1, k2 = self.kappa[:, 0].reshape(shape), self.kappa[:, 1].reshape(shape)
        return self.area.reshape(shape) * (1.0 - k1 * z) * (1.0 - k2 * z)

    def dz(self, T=None):
        """Cell thicknesses [m] (n, K), 0 beyond each count: each cell's span holds its volume m/rho(T) under the
        column's cross-section (module docstring). ValueError for a skin deeper than its radius of curvature."""
        valid = self.valid
        T = self.T() if T is None else T
        flat = np.zeros(self.m.shape)
        flat[valid] = self.m[valid] / (self.density(T[valid]) * np.broadcast_to(self.area[:, None], self.m.shape)[valid])
        curved = np.any(self.kappa != 0.0, axis=1)
        if not curved.any():
            return flat
        dz = flat.copy()
        rows = np.flatnonzero(curved)
        v = np.cumsum(flat[rows], axis=1)                     # volume per unit top area down to each cell's bottom
        sel = valid[rows]
        k1 = np.broadcast_to(self.kappa[rows, 0][:, None], v.shape)[sel]
        k2 = np.broadcast_to(self.kappa[rows, 1][:, None], v.shape)[sel]
        z_bot = np.zeros(v.shape)
        z_bot[sel] = _depth_of_volume(v[sel], k1, k2)
        z_top = np.concatenate([np.zeros((len(rows), 1)), z_bot[:, :-1]], axis=1)
        dz[rows] = np.where(sel, z_bot - z_top, 0.0)
        return dz

    def thickness(self):
        """Each column's thickness [m] (n,)."""
        return self.dz().sum(axis=1)

    def targets(self, dz=None):
        """Target thickness of every cell [m] (n, K) from the depth of its top below the column's top."""
        dz = self.dz() if dz is None else dz
        z_top = np.cumsum(dz, axis=1) - dz
        s = self.settings
        return s.cell + (s.growth - 1.0) * np.maximum(z_top - s.fine_zone, 0.0)

    def depth_above(self, T_level):
        """Per column, the depth [m] from the top to where the unbroken run of cells hotter than T_level ends,
        interpolated linearly in temperature between the centres of the last cell in the run and the first below it;
        0 if the top cell is not above T_level, the whole thickness if every cell is."""
        T, valid = self.T(), self.valid
        run = np.cumprod(valid & (T > T_level), axis=1).astype(bool).sum(axis=1)
        dz = self.dz(T)
        centre = np.cumsum(dz, axis=1) - 0.5 * dz
        r, K = np.arange(self.n), self.m.shape[1]
        a, b = np.clip(run - 1, 0, K - 1), np.clip(run, 0, K - 1)
        Ta, Tb = T[r, a], T[r, b]
        cooler = (run > 0) & (run < self.count)
        frac = np.zeros(self.n)
        frac[cooler] = np.clip((Ta[cooler] - T_level) / (Ta[cooler] - Tb[cooler]), 0.0, 1.0)
        between = centre[r, a] + frac * (centre[r, b] - centre[r, a])
        return np.where(run == 0, 0.0, np.where(run >= self.count, dz.sum(axis=1), between))

    # -- regrid -----------------------------------------------------------------------------------------------------
    def regrid(self):
        """Merge cells below MERGE_BELOW of their target, else split those above SPLIT_ABOVE, until none qualifies."""
        for _ in range(REGRID_PASSES):
            valid = self.valid
            dz = self.dz()
            t = self.targets(dz)
            small = valid & (dz < MERGE_BELOW * t) & (self.count >= 2)[:, None]
            big = valid & (dz > SPLIT_ABOVE * t)
            if small.any():
                self._merge(np.flatnonzero(small.any(axis=1)), small)
            elif big.any():
                self._split(np.flatnonzero(big.any(axis=1)), big)
            else:
                return
        raise RuntimeError("the skin cells did not settle in {} merge-and-split passes".format(REGRID_PASSES))

    def _merge(self, rows, small):
        """In each run of small cells (in `rows`), the first, third, ... cell takes the cell below it; a small bottom
        cell with none below to take joins the cell above. Mass, enthalpy and momentum add."""
        s, valid = small[rows], self.valid[rows]
        nr, K = s.shape
        j = np.arange(K)[None, :]
        last_other = np.maximum.accumulate(np.where(s, -1, j), axis=1)    # the last cell at or above that is not small
        leader = s & ((j - last_other - 1) % 2 == 0)                        # even place in its run of small cells
        join = np.zeros_like(s)
        join[:, 1:] = leader[:, :-1] & valid[:, 1:]
        r = np.arange(nr)
        bottom = self.count[rows] - 1
        join[r, bottom] |= leader[r, bottom] & (bottom >= 1)
        start = valid & ~join
        gid = np.cumsum(start, axis=1) - 1
        flat = (r[:, None] * K + gid)[valid]
        m0 = self.m[rows]
        m, H, P = np.zeros(nr * K), np.zeros(nr * K), np.zeros((nr * K, 3))
        np.add.at(m, flat, m0[valid])
        np.add.at(H, flat, self.H[rows][valid])
        np.add.at(P, flat, (m0[..., None] * self.U[rows])[valid])
        if self.dH is not None:
            dH = np.zeros(nr * K)
            np.add.at(dH, flat, self.dH[rows][valid])
            self.dH[rows] = dH.reshape(nr, K)
        U = np.divide(P, m[:, None], out=np.zeros_like(P), where=m[:, None] > 0.0)
        self.m[rows], self.H[rows] = m.reshape(nr, K), H.reshape(nr, K)
        self.U[rows], self.count[rows] = U.reshape(nr, K, 3), start.sum(axis=1)
        self._trim()

    def _split(self, rows, big):
        """Every big cell in `rows` becomes two halves of its mass and enthalpy, both at its velocity."""
        valid = self.valid[rows]
        reps = np.where(valid, 1 + big[rows], 0)
        nr, K = reps.shape
        count = reps.sum(axis=1)
        self._width(int(count.max()))
        K2 = self.m.shape[1]
        dest = np.cumsum(reps, axis=1) - reps
        rr = np.broadcast_to(np.arange(nr)[:, None], (nr, K))
        share = 1.0 / np.maximum(reps, 1)
        m, H, U = np.zeros((nr, K2)), np.zeros((nr, K2)), np.zeros((nr, K2, 3))
        m0, H0, U0 = self.m[rows, :K] * share, self.H[rows, :K] * share, self.U[rows, :K]
        trial = self.dH is not None
        if trial:
            dH, dH0 = np.zeros((nr, K2)), self.dH[rows, :K] * share
        for half in (0, 1):
            sel = valid & (half < reps)
            m[rr[sel], dest[sel] + half] = m0[sel]
            H[rr[sel], dest[sel] + half] = H0[sel]
            U[rr[sel], dest[sel] + half] = U0[sel]
            if trial:
                dH[rr[sel], dest[sel] + half] = dH0[sel]
        if trial:
            self.dH[rows] = dH
        self.m[rows], self.H[rows], self.U[rows], self.count[rows] = m, H, U, count

    # -- the bottom and the top -------------------------------------------------------------------------------------
    def add_bottom(self, i, m, H):
        """Metal of mass m >= 0 [kg] and enthalpy H [J] (above T_REF) joins the bottom cell of each column i (an index
        array; a repeated index adds each entry; an empty column gets a new cell). It arrives at rest, as rigid metal
        from the volume does: the bottom cell keeps its momentum."""
        i = np.asarray(i, dtype=np.int64).reshape(-1)
        m = np.broadcast_to(np.asarray(m, dtype=float), i.shape)
        H = np.broadcast_to(np.asarray(H, dtype=float), i.shape)
        if np.any(~np.isfinite(m) | (m < 0.0)) or np.any(~np.isfinite(H)):
            raise ValueError("metal added at a skin's bottom must be a finite, non-negative mass with finite enthalpy")
        some = m > 0.0
        i, m, H = i[some], m[some], H[some]
        if not len(i):
            return
        self._width(1)
        cols = np.unique(i)
        b = np.maximum(self.count[cols] - 1, 0)
        P = self.m[cols, b, None] * self.U[cols, b]
        bottom = np.zeros(self.n, dtype=np.int64)
        bottom[cols] = b
        np.add.at(self.m, (i, bottom[i]), m)
        np.add.at(self.H, (i, bottom[i]), H)
        self.U[cols, b] = P / self.m[cols, b, None]
        self.count[cols] = np.maximum(self.count[cols], 1)

    def take_top(self, i, k):
        """Remove the top k[j] cells of each column i[j] (distinct indices, k <= count). Returns per entry the
        (mass [kg], enthalpy [J], momentum [kg m/s] (., 3)) removed; the cells beneath move up."""
        i = np.asarray(i, dtype=np.int64).reshape(-1)
        k = np.broadcast_to(np.asarray(k, dtype=np.int64), i.shape)
        if len(np.unique(i)) != len(i) or np.any(k < 0) or np.any(k > self.count[i]):
            raise ValueError("take_top needs distinct columns and 0 <= k <= count")
        K = self.m.shape[1]
        top = np.arange(K)[None, :] < k[:, None]
        m0 = self.m[i]
        mass, H = (m0 * top).sum(axis=1), (self.H[i] * top).sum(axis=1)
        P = ((m0 * top)[..., None] * self.U[i]).sum(axis=1)
        src = np.minimum(np.arange(K)[None, :] + k[:, None], K - 1)
        count = self.count[i] - k
        keep = np.arange(K)[None, :] < count[:, None]
        self.m[i] = np.where(keep, np.take_along_axis(m0, src, axis=1), 0.0)
        self.H[i] = np.where(keep, np.take_along_axis(self.H[i], src, axis=1), 0.0)
        self.U[i] = np.where(keep[..., None], np.take_along_axis(self.U[i], src[..., None], axis=1), 0.0)
        if self.dH is not None:
            self.dH[i] = np.where(keep, np.take_along_axis(self.dH[i], src, axis=1), 0.0)
        self.count[i] = count
        self._trim()
        return mass, H, P

    # -- trial passes -----------------------------------------------------------------------------------------------
    def state(self):
        """A copy of the cells, for a trial pass to be undone with `set_state`."""
        return self.m.copy(), self.H.copy(), self.U.copy(), self.count.copy()

    def set_state(self, state):
        """The cells as `state` saved them; an open trial pass is closed."""
        m, H, U, count = state
        self.m, self.H, self.U, self.count = m.copy(), H.copy(), U.copy(), count.copy()
        self._end_trial()

    def _end_trial(self):
        self.dH = None
        self._trial_E = self._trial_dE = self._trial_Tb = None

    @property
    def S(self):
        """The trial pass's sensitivity dT/dT_b of every cell (n, K), 0 beyond each count; None outside a trial. Carried
        as each cell's dH/dT_b (`dH`, extensive, so that merges, splits, takes and draws keep it exactly) and read here
        as dH / (m cp_eff(T))."""
        if self.dH is None:
            return None
        valid, T = self.valid, self.T()
        S = np.zeros(self.m.shape)
        S[valid] = self.dH[valid] / (self.m[valid] * self.material.cp_eff(T[valid]))
        return S

    def base_law(self, dt):
        """The linear base law of the trial pass now open, (a [W], b [W/K]) per column: the heat into the volume over
        the pass, E_base / dt = a + b T_b, from the pass's summed base heat and its derivative with respect to the base
        temperature (the skyscrapers' trial pass), the slope clipped at zero so that the volume's matrix keeps its
        diagonal dominance (what the clip misses is the interface mismatch's to book). `dt` [s] is the pass's length;
        T_b the base temperature its sub-steps held. Closes the trial (the caller restores the cells with `set_state`)."""
        if self._trial_E is None:
            raise ValueError("base_law needs a trial pass: conduct(..., trial=True) first")
        b = np.minimum(self._trial_dE / dt, 0.0)
        a = self._trial_E / dt - b * self._trial_Tb
        self._end_trial()
        return a, b

    # -- heat through the depth -------------------------------------------------------------------------------------
    def conduct(self, dt, q_top, T_amb, emissivity, T_b, trial=False):
        """One implicit (backward-Euler) heat sub-step dt [s] of every column through its depth: the convective flux
        q_top [W/m2] over the top area A(0) in at the top cell, radiation emissivity * SIGMA_SB (T0^4 - T_amb^4) over
        the same area out of it (the shared Stefan-Boltzmann constant), conduction between cells through harmonic
        conductances A(z_i) / (dz_i / 2k_i + dz_j / 2k_j) with the cross-section at the interface's depth z_i, and the
        base held at T_b [K] through the bottom half-cell, A(z_bottom) k / (dz / 2) (T_b None: an adiabatic base).
        The cells' thicknesses (and so the interface depths) are those at the sub-step's start. Newton in the
        temperature with a backtracking line search per column on the residual's norm, the Jacobian with the
        conductivity's slope (dk/dT by a central difference of the material's k). An empty column passes its heat
        straight to the base.

        With trial=True the sub-step belongs to a trial pass (opened by the first such call, closed by `base_law` or
        `set_state`): each cell's dH/dT_b is carried from the last sub-step and updated by the converged step's own
        linearisation, J S = dH_old / dt + G_b e_bottom, so the book's dE_base_dTb is the exact derivative of E_base
        for the linearised step (and of the pass, through the carried dH). Returns a HeatBook."""
        n = self.n
        valid, rows = self.valid, np.arange(n)
        has = self.count > 0
        last = np.maximum(self.count - 1, 0)
        mat = self.material
        q = np.broadcast_to(np.asarray(q_top, dtype=float), (n,)) * self.area
        es = np.broadcast_to(np.asarray(emissivity, dtype=float), (n,)) * SIGMA_SB * self.area
        T_amb = np.broadcast_to(np.asarray(T_amb, dtype=float), (n,))
        base = T_b is not None
        T_b = np.broadcast_to(np.asarray(T_b if base else 0.0, dtype=float), (n,))
        if trial and self.dH is None:
            self.dH = np.zeros(self.m.shape)
            self._trial_E, self._trial_dE = np.zeros(n), np.zeros(n)
        elif not trial and self.dH is not None:
            raise ValueError("a trial pass is open: close it with base_law or set_state before a real sub-step")

        T_old = self.T()                                  # once per sub-step: the material's inversion is the cost
        fill = np.where(has, T_old[rows, last], 1.0)
        T_old = np.where(valid, T_old, fill[:, None])     # padded cells sit at their column's bottom temperature
        dz = np.where(valid, np.maximum(self.dz(T_old), DZ_FLOOR), 1.0)
        z_bot = np.cumsum(np.where(valid, dz, 0.0), axis=1)
        link = valid[:, 1:]
        A_i = np.where(link, self.cross_section(z_bot[:, :-1]), 1.0)          # at each interface's depth
        A_b = self.cross_section(z_bot[rows, last][:, None])[:, 0] if base else np.zeros(n)
        bottom = has & base
        m, E_old = self.m, self.H

        def system(T):
            """The residual R, the Jacobian's (lower, diag, upper) and the base conductance at T."""
            h, cp, k = mat.enthalpy(T), mat.cp_eff(T), mat.k(T)
            kp = (mat.k(T + DK_STEP) - mat.k(T - DK_STEP)) / (2.0 * DK_STEP)
            r = 0.5 * dz / k
            G = np.where(link, A_i / (r[:, :-1] + r[:, 1:]), 0.0)
            dT = T[:, :-1] - T[:, 1:]
            F = G * dT                                                          # down through each interface
            dG = G * G / A_i * 0.5                                              # dG/dk_j = dG * dz_j / k_j^2
            dFl = G + dT * dG * dz[:, :-1] * kp[:, :-1] / k[:, :-1] ** 2
            dFr = -G + dT * dG * dz[:, 1:] * kp[:, 1:] / k[:, 1:] ** 2
            Tl, kl, kpl, dzl = T[rows, last], k[rows, last], kp[rows, last], dz[rows, last]
            Gb = np.where(bottom, A_b * kl / (0.5 * dzl), 0.0)
            Fb = Gb * (Tl - T_b)
            dFb = Gb + np.where(bottom, (Tl - T_b) * A_b * kpl / (0.5 * dzl), 0.0)
            R = (m * h - E_old) / dt
            R[:, :-1] += F
            R[:, 1:] -= F
            R[:, 0] -= q - es * (T[:, 0] ** 4 - T_amb ** 4)
            R[rows, last] += Fb
            D = m * cp / dt
            D[:, :-1] += dFl
            D[:, 1:] -= dFr
            D[:, 0] += 4.0 * es * T[:, 0] ** 3
            D[rows, last] += dFb
            lower, upper = np.zeros_like(D), np.zeros_like(D)
            upper[:, :-1] = dFr
            lower[:, 1:] = -dFl
            return (np.where(valid, R, 0.0), np.where(valid, lower, 0.0), np.where(valid, D, 1.0),
                    np.where(valid, upper, 0.0), Gb, cp)

        T = T_old
        R, lo, D, up, Gb, cp = system(T)
        moved = np.zeros(n)
        for _ in range(NEWTON_MAX):
            dT = _tridiagonal(lo, D, up, -R)
            r0 = np.sqrt((R * R).sum(axis=1))
            lam = np.ones(n)
            for _ in range(LINE_SEARCH_MAX):
                T_try = T + lam[:, None] * dT
                R_try, lo_t, D_t, up_t, Gb_t, cp_t = system(T_try)
                moved = np.abs(T_try - T).max(axis=1)
                worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
                if not worse.any():
                    break
                lam = np.where(worse, 0.5 * lam, lam)
            T, R, lo, D, up, Gb, cp = T_try, R_try, lo_t, D_t, up_t, Gb_t, cp_t
            if moved.max() <= NEWTON_TOL:
                break
        else:
            raise RuntimeError("a skin heat sub-step did not converge in {} iterations (last change {:.2e} K)".format(
                NEWTON_MAX, moved.max()))
        E_new = np.where(valid, m * mat.enthalpy(T), 0.0)
        E_top = q * dt
        E_rad = np.where(has, es * (T[:, 0] ** 4 - T_amb ** 4) * dt, 0.0)
        E_side = np.zeros(n)
        E_base = E_top - E_rad + E_side - (E_new.sum(axis=1) - E_old.sum(axis=1))
        book = HeatBook(E_top, E_rad, E_base, E_side)
        if trial:
            rhs = self.dH / dt
            rhs[rows, last] += Gb                                  # -dR/dT_b: the base's pull on the bottom cell
            S = _tridiagonal(lo, D, up, np.where(valid, rhs, 0.0))   # the converged step's own Jacobian
            dH_new = np.where(valid, m * cp * S, 0.0)
            dE = np.where(has, -4.0 * es * T[:, 0] ** 3 * S[:, 0] * dt, 0.0) - (dH_new - self.dH).sum(axis=1)
            book.dE_base_dTb = dE
            self.dH = dH_new
            self._trial_E += E_base
            self._trial_dE += dE
            self._trial_Tb = np.array(T_b, dtype=float)
        self.H = E_new
        return book

    # -- removal, draw and hand-down --------------------------------------------------------------------------------
    def remove_feed(self):
        """The flow-off removal (the plates' case): the contiguous top cells at or above the material's T_feed leave,
        with their momentum (kept in `momentum_out`). Returns (mass [kg], enthalpy [J]) per column."""
        hot = self.valid & (self.T() >= self.material.T_feed)
        k = np.cumprod(hot, axis=1).sum(axis=1)
        mass, H = np.zeros(self.n), np.zeros(self.n)
        self.momentum_out = np.zeros((self.n, 3))
        i = np.flatnonzero(k > 0)
        if len(i):
            mass[i], H[i], self.momentum_out[i] = self.take_top(i, k[i])
        return mass, H

    def base_area(self):
        """The cross-section [m2] (n,) at each column's bottom, A(z_bottom); A0 for an empty column."""
        return self.cross_section(self.thickness()[:, None])[:, 0]

    def minimum_mass(self):
        """Each column's minimum mass [kg] (n,): rho_s A(z_bottom) stub, a stub of solid metal over the base's
        cross-section (the skyscrapers' minimum as flown, without Girin's delta_m term)."""
        return self.rho_s * self.base_area() * self.settings.stub

    def draw(self, T_b, g_base, m_min, s0=0.0):
        """A column below m_min [kg] takes the shortfall at its bottom from the volume: rigid metal whose specific
        enthalpy is the average, by PROFILE_POINTS-point Gauss-Legendre, of h along the volume's profile
        T(s) = T_b - g_base s [K, K/m] capped at the material's T_rigid, over the depths s0 to s0 + dd below the base,
        dd = shortfall / (rho_s A(z_bottom)) (the layer it takes, as solid). Returns (mass [kg], enthalpy [J]) per
        column; in a trial pass the drawn metal carries its enthalpy's derivative with respect to T_b."""
        n = self.n
        T_b = np.broadcast_to(np.asarray(T_b, dtype=float), (n,)).copy()
        g = np.broadcast_to(np.asarray(0.0 if g_base is None else g_base, dtype=float), (n,)).copy()
        m_min = np.broadcast_to(np.asarray(m_min, dtype=float), (n,))
        s0 = np.broadcast_to(np.asarray(s0, dtype=float), (n,)).copy()
        short = np.maximum(m_min - self.m.sum(axis=1), 0.0)
        H = np.zeros(n)
        i = np.flatnonzero(short > 0.0)
        if not len(i):
            return short, H
        dd = short[i] / (self.rho_s * self.base_area()[i])
        T, w, below = _profile_points(T_b[i], g[i], s0[i], s0[i] + dd, self.material.T_rigid)
        h = (self.material.enthalpy(T) * w).sum(axis=1)
        H[i] = short[i] * h
        self.add_bottom(i, short[i], H[i])
        if self.dH is not None:
            self._width(1)
            dh = (self.material.cp_eff(T) * below * w).sum(axis=1)
            self.dH[i, self.count[i] - 1] += short[i] * dh
        return short, H

    def anchored(self):
        """The rigid cells joined to the floor (n, K): each column's unbroken run of cells at or below the material's
        T_rigid counted from its bottom cell up. A column whose bottom cell is not rigid has none."""
        return self._anchored(self.T())

    def _anchored(self, T):
        valid, K = self.valid, self.m.shape[1]
        j = np.arange(K)[None, :]
        up = np.clip(self.count[:, None] - 1 - j, 0, K - 1)                       # the j-th cell from the bottom
        rigid_up = np.take_along_axis(valid & (T <= self.material.T_rigid), up, axis=1) & (j < self.count[:, None])
        run = np.cumprod(rigid_up, axis=1).sum(axis=1)
        return valid & (j >= (self.count - run)[:, None])

    def crust(self):
        """Rigid cells (at or below T_rigid) not joined to the floor: non-rigid metal lies between them and it."""
        T = self.T()
        return self.valid & (T <= self.material.T_rigid) & ~self._anchored(T)

    def hand_down(self, stub):
        """The anchored cells lying wholly below the shallowest `stub` [m] of their anchored run leave the skin for the
        volume, whole and at their own enthalpy (design 20.6, tightened 2026-11-10); a crust never goes, and a cell
        reaching above that depth stays. Returns (mass [kg], enthalpy [J], volume [m3]) per column, the volume each
        cell's m / rho(T), the span it filled under the column's cross-section."""
        T = self.T()
        valid = self.valid
        anch = self._anchored(T)
        dz = self.dz(T)
        z_top = np.cumsum(dz, axis=1) - dz
        run_top = np.where(anch, z_top, np.inf).min(axis=1)
        go = anch & (z_top >= (run_top + stub)[:, None] * (1.0 - 1e-12))
        mass = np.where(go, self.m, 0.0).sum(axis=1)
        H = np.where(go, self.H, 0.0).sum(axis=1)
        vol = np.where(go, self.m / np.where(valid, self.density(T), 1.0), 0.0).sum(axis=1)
        k = go.sum(axis=1)
        if k.any():                                            # the leaving cells are the bottom k of each column
            keep = np.arange(self.m.shape[1])[None, :] < (self.count - k)[:, None]
            self.m, self.H = np.where(keep, self.m, 0.0), np.where(keep, self.H, 0.0)
            self.U = np.where(keep[..., None], self.U, 0.0)
            if self.dH is not None:
                self.dH = np.where(keep, self.dH, 0.0)
            self.count = self.count - k
            self._trim()
        return mass, H, vol
