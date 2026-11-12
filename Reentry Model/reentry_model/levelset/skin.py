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

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (design 20.6; skyscrapers 6.2)
MERGE_BELOW = 0.5          # a cell thinner than this times its target joins a neighbour (skyscrapers' rule)
SPLIT_ABOVE = 1.5          # a cell thicker than this times its target is split in two (skyscrapers' rule)
REGRID_PASSES = 64         # merge-or-split rounds of one regrid at most (skyscrapers' rule)
DEPTH_NEWTON_MAX = 100     # iterations of the cross-section's depth solve at most (it takes about five)
DEPTH_NEWTON_RTOL = 4.0 * np.finfo(float).eps   # the depth solve stops when a step is below this times the depth
                                                # or the volume's residual below this times the volume (round-off)


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

    def _trim(self):
        """Drop padding beyond the longest column (keeping at least one cell of width)."""
        K = max(int(self.count.max()) if self.n else 0, 1)
        if K < self.m.shape[1]:
            self.m, self.H, self.U = self.m[:, :K], self.H[:, :K], self.U[:, :K]

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
        for half in (0, 1):
            sel = valid & (half < reps)
            m[rr[sel], dest[sel] + half] = m0[sel]
            H[rr[sel], dest[sel] + half] = H0[sel]
            U[rr[sel], dest[sel] + half] = U0[sel]
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
        self.count[i] = count
        self._trim()
        return mass, H, P
