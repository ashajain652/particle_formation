"""The layered skins on the corner columns (layered-skin design, docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-
design.md sections 6 and 7; stage 2, sub-plan 20: heat only, liquid removed as it forms).

Every corner column carries a skin from time zero: a stack of cells along the column's direction holding all metal
outward of the rigid front, at least the column's minimum mass. A cell holds its mass [kg] and enthalpy [J above
material.T_REF]; its temperature, liquid fraction and density follow from the material's enthalpy law, and its thickness
is its mass over the column's area and its own density (liquid is about 15 % less dense than solid). All columns are
stored side by side in padded arrays (n, K): cell 0 is the top, cells 0 .. count - 1 are real, the rest hold nothing.

Cells move with the metal: nothing re-divides a whole column. After every sub-step's exchanges a cell more than
SPLIT_ABOVE times its target thickness is split in two and one less than MERGE_BELOW times it is merged into the cell
beneath (the bottom cell into the one above), exactly in mass and enthalpy. The target is the fine cell within the fine
zone (`fine_zone` below the top) and grows by (growth - 1) per unit depth below it -- the continuous form of each cell
being `growth` times thicker than the one above.

The minimum and the fine zone are set apart (Asha, 2026-10-16, option B'): a skin holds at least 0.8 mm of metal, so
that where the non-rigid layer is thinner -- a cold body melting fast -- its bottom is rigid metal and the skin-3D
interface lies below the rigid front, off the steepest part of the temperature profile, while the finest cells stay in
the top 0.4 mm. With 0.4 mm for both the level-2 cold plate's interface sat on the rigid front (its non-rigid layer is
0.39-0.49 mm), the 3D base node rose 1.6 K above T_rigid within each part before absorption reset it, and the plate
melted 3.6 % slow from 1.0 to 2.0 s; with B' 1.7 % (measured 2026-10-16).

Stage 2 has no flow and no spraying: after each sub-step the contiguous cells at the top that are at or above T_feed
leave the body with their enthalpy (`removed`), and a column below its minimum draws rigid metal from the 3D model
(`drawn`) at the base temperature, capped at T_rigid: the drawn metal's enthalpy is h(min(T_b, T_rigid)), or, given the
3D model's temperature gradient below the base, h(min(T, T_rigid)) averaged over the depth each draw takes along the
3D model's linear profile (`profile_enthalpy`; Asha, 2026-10-12). The body
(`skinbody.SkinBody`) turns the drawn, absorbed and handed-down metal into the corners' displacements and reconciles
each column with what the moving mesh actually swept."""
import math
from dataclasses import dataclass

import numpy as np

from .thermal import SIGMA_SB

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (spec section 6.2)
SPLIT_ABOVE = 1.5                 # a cell thicker than this times its target is split in two (spec 6.2)
MERGE_BELOW = 0.5                 # a cell thinner than this times its target is merged into its neighbour (spec 6.2)
REGRID_PASSES = 64                # merge-and-split rounds after an exchange: one unless a cell grew many-fold, each
                                  # round halving it (2.6 mm added at once needs 9)
DELTA_M_MARGIN = 50.0e-6          # m of liquid kept below Girin's conjugate depth in the minimum mass (spec 6.1)
INTERFACE_TOLERANCE = 1.0e-3      # |mismatch| / (|E_top| + |E_base|): a macro step repeats above it (spec 7.2)
JUMP_TOLERANCE = 5.0              # K between the metal crossing the base and the skin's bottom cell: a fault (spec 6.6,
                                  # provisional)
NEWTON_TOL = 1.0e-6               # K: a sub-step has converged when no cell moves by more than this
NEWTON_MAX = 40
LINE_SEARCH_MAX = 12              # halvings of a Newton step before it is taken anyway
DZ_FLOOR = 1.0e-12                # m: the conduction length of a vanishing cell
SIDEWAYS_CFL = 0.9                # explicit sideways conduction: sub-cycles keep dt k sum(w) / (rho c_p A) below this
TABLE_DT = 0.01                   # K, spacing of the skins' property tables (PropertyTable)
TABLE_RANGE = (150.0, 3000.0)     # K, their extent; beyond it the end intervals are extended linearly
PROFILE_POINTS = 8                # Gauss-Legendre points of an enthalpy averaged over a depth of the 3D profile


# -- settings and cells (Task 2) ---------------------------------------------------------------------------------------
@dataclass
class SkinSettings:
    min_thickness: float = 0.8e-3          # m of solid metal: the minimum mass (spec 6.1); was 0.4 mm before option B'
    cell: float = SKIN_PRESETS["production"]   # m, the fine cell
    fine_zone: float = 0.4e-3              # m below the top held in fine cells (spec 6.2); the minimum's depth before B'
    growth: float = 1.10                   # each cell below the fine zone this much thicker than the one above
    substep: float = 0.010                 # s, the skins' sub-step (spec 7.1)
    swell: bool = True                     # each cell at its own phase's density, so melting swells the skin (spec 5.2);
                                           # False holds every cell at the solid's density, the setting levels 2 and 2b
                                           # verify with (Asha, 2026-10-14: their references hold one density, and the
                                           # swelling alone melts the exact hot plate 0.094 mm sooner by 1.0 s)

    def __post_init__(self):
        if not (self.min_thickness > 0.0 and self.fine_zone > 0.0 and self.cell > 0.0 and self.substep > 0.0):
            raise ValueError("skin minimum thickness, fine zone, cell and substep must be positive")
        if self.growth < 1.0:
            raise ValueError("skin growth must be at least 1")
        if 2.0 * self.cell > min(self.min_thickness, self.fine_zone) * (1.0 + 1e-12):
            raise ValueError("a skin needs at least two fine cells: cell <= min(min_thickness, fine_zone) / 2")


def thomas(sub, diag, sup, rhs):
    """Row-wise tridiagonal solve: sub and sup (n, k-1), diag and rhs (n, k). Returns x (n, k)."""
    n, k = diag.shape
    c = np.zeros((n, max(k - 1, 1)))
    d = np.zeros((n, k))
    if k > 1:
        c[:, 0] = sup[:, 0] / diag[:, 0]
    d[:, 0] = rhs[:, 0] / diag[:, 0]
    for i in range(1, k):
        denom = diag[:, i] - sub[:, i - 1] * c[:, i - 1]
        if i < k - 1:
            c[:, i] = sup[:, i] / denom
        d[:, i] = (rhs[:, i] - sub[:, i - 1] * d[:, i - 1]) / denom
    x = np.zeros((n, k))
    x[:, -1] = d[:, -1]
    for i in range(k - 2, -1, -1):
        x[:, i] = d[:, i] - c[:, i] * x[:, i + 1]
    return x


class PropertyTable:
    """The material's h(T), c_p,eff(T), k(T) and f_l(T) tabulated on a uniform TABLE_DT grid and interpolated
    linearly: a lookup is index arithmetic instead of the material's binary search, which was half of a skin pass's
    time (measured while planning: 11 s of a 4 mm sphere's macro step, 5.5 s of it in the property calls). The
    interpolation moves h by at most about 60 J/kg (0.05 K) across a Scheil node, and the books never see it: every pass
    starts from each cell's stored enthalpy."""

    def __init__(self, material, dT=TABLE_DT, T_range=TABLE_RANGE):
        self.T0, self.dT = float(T_range[0]), float(dT)
        T = np.arange(T_range[0], T_range[1] + dT, dT)
        self.n = len(T)
        self._h, self._cp, self._k, self._f = (material.enthalpy(T), material.cp_eff(T), material.k(T),
                                               material.liquid_fraction(T))

    def _at(self, table, T):
        x = (np.asarray(T, dtype=float) - self.T0) / self.dT
        i = np.clip(x.astype(np.int64), 0, self.n - 2)
        return table[i] + (x - i) * (table[i + 1] - table[i])

    def conduction(self, T):
        """(h, dh/dT, k) at T with one index computation, for the conduction's Newton iteration. The heat capacity is
        the slope of the tabulated h itself, so the Newton tangent is exact: interpolating c_p,eff instead left it off
        the slope near the Scheil curve's 1 K nodes and took five iterations a sub-step where three do (measured while
        planning)."""
        x = (np.asarray(T, dtype=float) - self.T0) / self.dT
        i = np.clip(x.astype(np.int64), 0, self.n - 2)
        dh = self._h[i + 1] - self._h[i]
        return self._h[i] + (x - i) * dh, dh / self.dT, self._k[i] + (x - i) * (self._k[i + 1] - self._k[i])

    def enthalpy(self, T):
        return self._at(self._h, T)

    def cp_eff(self, T):
        return self._at(self._cp, T)

    def k(self, T):
        return self._at(self._k, T)

    def liquid_fraction(self, T):
        return np.clip(self._at(self._f, T), 0.0, 1.0)


class Skins:
    """The skins of every column, as padded arrays (module docstring). `area` [m2] is each column's area, the dual
    area of its base corner; the body updates it when the corners move (the cells keep their mass). `S`, when not None,
    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b)."""

    def __init__(self, material, area, settings=None):
        self.material = material
        self.settings = settings or SkinSettings()
        self.area = np.array(area, dtype=float)
        n = len(self.area)
        self.m, self.H, self.T = np.zeros((n, 1)), np.zeros((n, 1)), np.full((n, 1), 300.0)
        self.count = np.zeros(n, dtype=np.int64)
        self.S = None
        self.rho_l = material.liquid.rho if material.liquid is not None and self.settings.swell else material.rho
        self.props = PropertyTable(material)

    # -- books ------------------------------------------------------------------------------------------------------
    @property
    def n(self):
        return len(self.area)

    @property
    def valid(self):
        return np.arange(self.m.shape[1])[None, :] < self.count[:, None]

    def mass(self):
        return float(self.m.sum())

    def energy(self):
        return float(self.H.sum())

    def column_mass(self):
        return self.m.sum(axis=1)

    def state(self):
        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count))

    def set_state(self, state):
        self.area, self.m, self.H, self.T, self.count = (a.copy() for a in state)
        self.S = None

    # -- geometry of the cells --------------------------------------------------------------------------------------
    def density(self, T):
        """Mixture density [kg/m3]: specific volumes add, (1 - f_l)/rho_s + f_l/rho_l."""
        f = self.props.liquid_fraction(T)
        return 1.0 / ((1.0 - f) / self.material.rho + f / self.rho_l)

    def dz(self):
        """Cell thicknesses [m] (0 beyond each column's count)."""
        return np.where(self.valid, self.m / (self.density(self.T) * self.area[:, None]), 0.0)

    def thickness(self):
        return self.dz().sum(axis=1)

    def targets(self, dz=None):
        """Target thickness of every cell [m] from the depth of its top below the column's top (module docstring)."""
        dz = self.dz() if dz is None else dz
        z = np.cumsum(dz, axis=1) - dz
        s = self.settings
        return s.cell + (s.growth - 1.0) * np.maximum(z - s.fine_zone, 0.0)

    # -- installing cells -------------------------------------------------------------------------------------------
    def _install(self, m, H, count, S=None, T=None):
        """Install padded cell arrays. Temperatures follow from the enthalpies, except where `T` is given for cells
        whose specific enthalpy the operation did not change (a shifted, split or partly taken cell); nan in T marks
        the cells to recompute."""
        valid = np.arange(m.shape[1])[None, :] < count[:, None]
        self.m, self.H, self.count = np.where(valid, m, 0.0), np.where(valid, H, 0.0), np.asarray(count, dtype=np.int64)
        redo = valid if T is None else valid & np.isnan(T)
        self.T = np.full(self.m.shape, 300.0) if T is None else np.where(valid & ~redo, T, 300.0)
        if redo.any():
            self.T[redo] = self.material.temperature_from_enthalpy(self.H[redo] / self.m[redo])
        if self.S is not None:
            self.S = np.where(valid, S, 0.0)

    def fill(self, mass, T):
        """Every column filled with `mass` [kg] (per column) of metal at temperature T [K] (per column), in equal cells of
        about the fine size (the initial skins, spec 6.1)."""
        mass, T = np.broadcast_to(np.asarray(mass, dtype=float), (self.n,)), np.broadcast_to(np.asarray(T, dtype=float), (self.n,))
        thick = mass / (self.density(T) * self.area)
        nc = np.maximum(np.rint(thick / self.settings.cell).astype(np.int64), 1)
        K = int(nc.max())
        valid = np.arange(K)[None, :] < nc[:, None]
        m = np.where(valid, (mass / nc)[:, None], 0.0)
        self.S = None
        self._install(m, m * self.material.enthalpy(T)[:, None], nc)

    def _width(self, K):
        """Pad the arrays to at least K cells."""
        extra = K - self.m.shape[1]
        if extra > 0:
            pad = lambda a, v: np.hstack([a, np.full((self.n, extra), v)])
            self.m, self.H, self.T = pad(self.m, 0.0), pad(self.H, 0.0), pad(self.T, 300.0)
            if self.S is not None:
                self.S = pad(self.S, 0.0)

    # -- split and merge --------------------------------------------------------------------------------------------
    def regrid(self):
        """Merge cells below MERGE_BELOW of their target, split those above SPLIT_ABOVE, until none qualifies. Only the
        columns that need it are rebuilt (a split every sub-step on every melting column is the common case)."""
        for _ in range(REGRID_PASSES):
            valid = self.valid
            dz = self.dz()
            t = self.targets(dz)
            small = valid & (dz < MERGE_BELOW * t) & (self.count >= 2)[:, None]
            big = valid & (dz > SPLIT_ABOVE * t)
            if not (small.any() or big.any()):
                return
            if small.any():
                self._merge(np.flatnonzero(small.any(axis=1)), small)
            else:
                self._split(np.flatnonzero(big.any(axis=1)), big)
        raise RuntimeError("the skin cells did not settle in {} split-and-merge passes".format(REGRID_PASSES))

    def _put_rows(self, rows, m, H, T, S, count):
        """Write rebuilt cell arrays for `rows` back (temperatures of nan are recomputed from the enthalpies), padding
        or trimming the arrays to the longest column."""
        self._width(m.shape[1])
        K = self.m.shape[1]
        pad = lambda a: np.hstack([a, np.zeros((len(rows), K - a.shape[1]))]) if a.shape[1] < K else a
        valid = np.arange(K)[None, :] < count[:, None]
        m, H, T = pad(m), pad(H), pad(T)
        redo = valid & np.isnan(T)
        if redo.any():
            T[redo] = self.material.temperature_from_enthalpy(H[redo] / m[redo])
        self.m[rows], self.H[rows] = np.where(valid, m, 0.0), np.where(valid, H, 0.0)
        self.T[rows] = np.where(valid, T, 300.0)
        self.count[rows] = count
        if self.S is not None:
            self.S[rows] = np.where(valid, pad(S), 0.0)
        width = int(self.count.max()) if self.n else 1
        if width < K:
            self.m, self.H, self.T = self.m[:, :width], self.H[:, :width], self.T[:, :width]
            if self.S is not None:
                self.S = self.S[:, :width]

    def _merge(self, rows, small):
        """Each run of small cells joins the cell beneath that ends it; a run at the bottom with no cell to end it joins
        the group above (in `rows` only)."""
        small = small[rows]
        valid = self.valid[rows]
        m0, H0, T0 = self.m[rows], self.H[rows], self.T[rows]
        n, K = m0.shape
        start = valid.copy()
        start[:, 1:] &= ~small[:, :-1]
        gid = np.cumsum(start, axis=1) - 1
        groups = start.sum(axis=1)
        r = np.arange(n)
        trailing = small[r, np.maximum(self.count[rows] - 1, 0)] & (groups >= 2)
        gid = np.where(trailing[:, None] & (gid == (groups - 1)[:, None]), gid - 1, gid)
        groups = groups - trailing
        flat = (r[:, None] * K + gid)[valid]
        m, H = np.zeros(n * K), np.zeros(n * K)
        np.add.at(m, flat, m0[valid])
        np.add.at(H, flat, H0[valid])
        size = np.bincount(flat, minlength=n * K)
        T = np.full(n * K, np.nan)
        T[flat] = np.where(size[flat] == 1, T0[valid], np.nan)             # a cell that merged with none keeps its T
        S = None
        if self.S is not None:
            mS = np.zeros(n * K)
            np.add.at(mS, flat, (m0 * self.S[rows])[valid])
            S = np.divide(mS, m, out=np.zeros_like(m), where=m > 0.0).reshape(n, K)
        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups)

    def _split(self, rows, big):
        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T and S)."""
        valid = self.valid[rows]
        reps = np.where(valid, 1 + big[rows], 0)
        n, K = reps.shape
        count = reps.sum(axis=1)
        K2 = int(count.max())
        pos = np.cumsum(reps, axis=1) - reps
        r_all = np.repeat(np.arange(n)[:, None], K, axis=1)
        m, H, T = np.zeros((n, K2)), np.zeros((n, K2)), np.zeros((n, K2))
        S = np.zeros((n, K2)) if self.S is not None else None
        share = 1.0 / reps.clip(1)
        for part in (0, 1):
            sel = valid & (part < reps)
            r, c = r_all[sel], pos[sel] + part
            m[r, c] = (self.m[rows] * share)[sel]
            H[r, c] = (self.H[rows] * share)[sel]
            T[r, c] = self.T[rows][sel]                                    # halves keep the specific enthalpy
            if S is not None:
                S[r, c] = self.S[rows][sel]
        self._put_rows(rows, m, H, T, S, count)

    # -- the top and the bottom -------------------------------------------------------------------------------------
    def take_top(self, k):
        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column."""
        k = np.asarray(k, dtype=np.int64)
        n, K = self.m.shape
        top = np.arange(K)[None, :] < k[:, None]
        mass, H = (self.m * top).sum(axis=1), (self.H * top).sum(axis=1)
        if k.any():
            src = np.minimum(np.arange(K)[None, :] + k[:, None], K - 1)
            shift = lambda a: np.take_along_axis(a, src, axis=1)
            S = shift(self.S) if self.S is not None else None
            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T))
        return mass, H

    def add_bottom(self, dm, h, s=1.0):
        """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
        a new cell); returns the enthalpy added. With S tracked, the metal carries the sensitivity s per column
        (default 1: it is at the base temperature and follows it; 0 for metal whose enthalpy does not depend on the
        base temperature, such as rigid metal drawn at T_rigid's enthalpy from a base above T_rigid)."""
        dm, h = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,)), np.broadcast_to(np.asarray(h, dtype=float), (self.n,))
        s = np.broadcast_to(np.asarray(s, dtype=float), (self.n,))
        some = dm > 0.0
        if not some.any():
            return np.zeros(self.n)
        self._width(int(np.max(np.maximum(self.count, 1))))
        rows = np.flatnonzero(some)
        b = np.maximum(self.count[rows] - 1, 0)
        m_old = self.m[rows, b]
        self.m[rows, b] += dm[rows]
        self.H[rows, b] += dm[rows] * h[rows]
        if self.S is not None:
            self.S[rows, b] = (m_old * self.S[rows, b] + dm[rows] * s[rows]) / self.m[rows, b]
        self.count[rows] = np.maximum(self.count[rows], 1)
        self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])
        return np.where(some, dm * h, 0.0)

    def take_bottom(self, dm):
        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy (a partly
        taken cell keeps its specific enthalpy); returns the enthalpy removed per column."""
        dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
        if np.any(dm > self.column_mass() * (1.0 + 1e-12) + 1e-300):
            raise ValueError("cannot take more metal from a skin's bottom than the skin holds")
        if not np.any(dm > 0.0):
            return np.zeros(self.n)
        n, K = self.m.shape
        valid = self.valid
        below = np.cumsum(self.m[:, ::-1], axis=1)[:, ::-1] - self.m      # mass beneath each cell
        taken = np.clip(dm[:, None] - below, 0.0, self.m)                 # how much of each cell goes
        frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
        H_taken = (self.H * frac).sum(axis=1)
        m, H = self.m - taken, self.H * (1.0 - frac)
        keep = valid & (m > 0.0)
        count = keep.sum(axis=1)                                          # the kept cells are the top ones
        self._install(m, H, count, self.S, self.T)                        # a partly taken cell keeps its T
        return H_taken

    def exchange_bottom(self, dm, dH):
        """Signed reconciliation at the bottom of every column: the skin gains dm [kg] (negative: loses) and its
        enthalpy changes by exactly dH [J]. Gained metal joins the bottom cell; lost metal leaves from the bottom up at
        the cells' own specific enthalpy, and the difference to dH stays in the new bottom cell."""
        dm, dH = np.asarray(dm, dtype=float), np.asarray(dH, dtype=float)
        gain = np.where(dm > 0.0, dm, 0.0)
        h = np.divide(dH, gain, out=np.zeros_like(dH), where=gain > 0.0)
        self.add_bottom(gain, h)
        rest = np.where(dm > 0.0, 0.0, dH)
        H_taken = self.take_bottom(np.where(dm < 0.0, -dm, 0.0))
        rest = rest + H_taken                                              # what is left over after the cells left
        rows = np.flatnonzero(rest != 0.0)
        if rows.size:
            b = self.count[rows] - 1
            if np.any(b < 0):
                raise ValueError("a skin emptied by a reconciliation cannot take the enthalpy difference")
            self.H[rows, b] += rest[rows]
            self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])

    def depth_above(self, T_level):
        """Per column: the depth [m] from the top to where the contiguous cells above T_level end, interpolated between
        cell centres (0 if the top cell is below it, the whole skin if every cell is above)."""
        valid = self.valid
        above = np.cumprod(valid & (self.T > T_level), axis=1).astype(bool)
        k = above.sum(axis=1)
        dz = self.dz()
        z_mid = np.cumsum(dz, axis=1) - 0.5 * dz
        rows = np.arange(self.n)
        Ka = np.clip(k - 1, 0, self.m.shape[1] - 1)
        Kb = np.clip(k, 0, self.m.shape[1] - 1)
        Ta, Tb = self.T[rows, Ka], self.T[rows, Kb]
        frac = np.clip(np.divide(Ta - T_level, Ta - Tb, out=np.zeros_like(Ta), where=Ta > Tb), 0.0, 1.0)
        inside = z_mid[rows, Ka] + frac * (z_mid[rows, Kb] - z_mid[rows, Ka])
        return np.where(k == 0, 0.0, np.where(k >= self.count, self.thickness(), inside))


# -- the minimum mass (Task 3) -----------------------------------------------------------------------------------------
def minimum_mass(material, area, settings, delta_m=None):
    """Each column's minimum mass [kg] (spec 6.1): the larger of `min_thickness` of solid metal over the column's area
    and the liquid filling Girin's conjugate depth delta_m [m] plus DELTA_M_MARGIN, so that the sheared layer always lies
    inside the fine cells. delta_m None (or nan where Girin's closure does not apply) leaves the solid term alone."""
    area = np.asarray(area, dtype=float)
    solid = material.rho * area * settings.min_thickness
    if delta_m is None:
        return solid
    rho_l = material.liquid.rho if material.liquid is not None else material.rho
    liquid = rho_l * area * (np.nan_to_num(np.asarray(delta_m, dtype=float), nan=-DELTA_M_MARGIN) + DELTA_M_MARGIN)
    return np.maximum(solid, liquid)


# -- heat through the depth (Task 4) -----------------------------------------------------------------------------------
def profile_enthalpy(material, T_top, g, s0, s1, cap=None):
    """The specific enthalpy [J/kg] averaged over the depths s0 to s1 [m] below a corner along the 3D model's linear
    temperature profile T(s) = T_top - g s (g [K/m] the temperature gradient along the outward normal), each depth's
    temperature capped at `cap` if given; s1 equal to s0 gives the value there. Metal swept or drawn from a layer whose
    temperature falls with depth carries this, not the surface node's enthalpy: on the steep Scheil curve the two differ
    by about 1e5 J/kg over a 1.5 mm layer falling 22 K/mm (the plate at 0.5 s, measured 2026-10-12)."""
    x, w = np.polynomial.legendre.leggauss(PROFILE_POINTS)
    T_top, g, s0, s1 = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (T_top, g, s0, s1)))
    s = 0.5 * (s0 + s1)[..., None] + 0.5 * (s1 - s0)[..., None] * x
    T = T_top[..., None] - g[..., None] * s
    if cap is not None:
        T = np.minimum(T, cap)
    return (material.enthalpy(T.ravel()).reshape(T.shape) * w).sum(axis=-1) / 2.0


@dataclass
class SkinBook:
    """What a pass did, per column: energies [J] and masses [kg]."""
    E_top: np.ndarray        # convective heat taken in at the top
    E_rad: np.ndarray        # radiated from the top
    E_base: np.ndarray       # conducted down into the 3D model through the base
    E_side: np.ndarray       # gained from the neighbouring columns by sideways conduction
    removed: np.ndarray      # liquid removed as it formed (stage 2)
    removed_H: np.ndarray    # the enthalpy it took
    drawn: np.ndarray        # rigid metal drawn from the 3D model in the sub-steps
    drawn_H: np.ndarray      # the enthalpy it brought


def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, Q_side=None):
    """One implicit (backward-Euler) sub-step for every column: the convective flux q_top [W/m2] over the top area
    A_top in at the top cell, radiation from it over the same area, conduction through the cells (harmonic
    conductances), the explicit sideways heat Q_side [W per cell], and the base held at T_b [K] through the bottom
    half-cell (T_b None: an adiabatic base, for the tests). Newton in the temperature with a backtracking line search
    per column (stage 1's solver, sub-plan 18). Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the
    skins track S, updates it to dT/dT_b and returns the base energy's derivative with respect to T_b as a fifth item.
    E_base is what the books leave, so they close exactly whatever the iteration error."""
    mat = sk.props
    n, K = sk.m.shape
    valid = sk.valid
    rows = np.arange(n)
    last = np.maximum(sk.count - 1, 0)
    A = sk.area
    dz = np.where(valid, np.maximum(sk.dz(), DZ_FLOOR), 1.0)
    m = sk.m
    link = valid[:, 1:]
    T_old = sk.T.copy()
    E_old = sk.H.copy()                              # the stored enthalpy: the books close whatever h(T) is tabulated
    eps_s = np.asarray(eps, dtype=float) * SIGMA_SB * np.asarray(A_top, dtype=float)
    q = np.asarray(q_top, dtype=float) * np.asarray(A_top, dtype=float)
    base = T_b is not None
    T_b = np.asarray(T_b if base else 0.0, dtype=float) * np.ones(n)
    Qs = np.zeros((n, K)) if Q_side is None else np.where(valid, Q_side, 0.0)

    def system(T):
        h, cp, k = mat.conduction(T)
        G = np.where(link, A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:]), 0.0)
        Gb = A * k[rows, last] / (0.5 * dz[rows, last]) if base else np.zeros(n)
        flux = G * (T[:, :-1] - T[:, 1:])
        net = Qs.copy()
        net[:, :-1] -= flux
        net[:, 1:] += flux
        net[:, 0] += q - eps_s * (T[:, 0] ** 4 - T_amb ** 4)
        net[rows, last] -= Gb * (T[rows, last] - T_b)
        diag = m * cp / dt
        diag[:, 0] += 4.0 * eps_s * T[:, 0] ** 3
        diag[:, :-1] += G
        diag[:, 1:] += G
        diag[rows, last] += Gb
        R = np.where(valid, (m * h - E_old) / dt - net, 0.0)
        return R, G, Gb, np.where(valid, diag, 1.0)

    # Newton in the temperature with a backtracking line search per column on the residual's norm: in 10-micrometre
    # cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step, so a cell must move
    # with its neighbours, and the 3D solver's per-node enthalpy-consistent update stalled here (sub-plan 18)
    T = T_old.copy()
    R, G, Gb, diag = system(T)
    moved = np.zeros(n)
    for _ in range(NEWTON_MAX):
        dT = thomas(-G, diag, -G, -R)
        r0 = np.sqrt((R * R).sum(axis=1))
        lam = np.ones(n)
        for _ in range(LINE_SEARCH_MAX):
            T_try = T + lam[:, None] * dT
            R_try, G_try, Gb_try, diag_try = system(T_try)
            moved = np.abs(T_try - T).max(axis=1)
            worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
            if not worse.any():
                break
            lam = np.where(worse, 0.5 * lam, lam)
        T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try
        if moved.max() <= NEWTON_TOL:
            break
    else:
        raise RuntimeError("a skin sub-step did not converge in {} iterations (last change {:.2e} K)".format(
            NEWTON_MAX, moved.max()))
    T = np.where(valid, T, T_old)
    E_new = m * mat.enthalpy(T)
    E_top = q * dt
    E_rad = eps_s * (T[:, 0] ** 4 - T_amb ** 4) * dt
    E_side = Qs.sum(axis=1) * dt
    E_base = E_top - E_rad + E_side - (E_new.sum(axis=1) - E_old.sum(axis=1))
    out = [E_top, E_rad, E_base, E_side]
    if sk.S is not None:
        rhs = m * mat.cp_eff(T_old) / dt * sk.S
        rhs[rows, last] += Gb
        S = thomas(-G, diag, -G, np.where(valid, rhs, 0.0))
        out.append(dt * Gb * (S[rows, last] - 1.0))
        sk.S = np.where(valid, S, 0.0)
    sk.T, sk.H = T, np.where(valid, E_new, 0.0)
    return tuple(out)


def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None):
    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the sideways heat `side(sk, h)` [W per
    cell] if given (Task 5), the implicit conduction, then -- with `remove` and a melting material -- the contiguous top
    cells at or above T_feed leave, a column below its minimum mass m_min [kg] draws back to it, and the cells are split
    and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but never above T_rigid -- a 3D base
    node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it reached 930 K in one 0.5 s step, found
    while planning), and drawing at its temperature would draw liquid that the next sub-step removes, without limit;
    what the 3D model actually loses is reconciled by the body after the base moves. With `g_base` [K/m per column], the
    3D model's temperature gradient along the outward normal below the base, each draw takes the next layer down: its
    enthalpy is h(min(T, T_rigid)) averaged along T(s) = T_b - g_base s over the depth it takes (`profile_enthalpy`).

    T_b1 None: the trial pass, on a copy (the skins are left as they were), base held at T_b0; returns the base law
    (a, b) per column, E_base / dt = a + b T_b [W into the 3D model], from the pass's base energy and its derivative with
    respect to the base temperature, the slope clipped at zero (stage 1's choice: a law rising with T_b would cost the 3D
    matrix its diagonal dominance; what the clip misses is the interface mismatch's to book). Otherwise the real pass,
    base held at T_b1, the 3D model's end-of-step base temperature; returns a SkinBook."""
    mat = sk.material
    N = max(1, int(math.ceil(dt / sk.settings.substep - 1e-9)))
    h = dt / N
    trial = T_b1 is None
    saved = sk.state() if trial else None
    T_b = np.asarray(T_b0 if trial else T_b1, dtype=float) * np.ones(sk.n)
    book = SkinBook(*(np.zeros(sk.n) for _ in range(8)))
    dEb = np.zeros(sk.n)
    if trial:
        sk.S = np.zeros_like(sk.m)
    h_b = mat.enthalpy(np.minimum(T_b, mat.T_rigid))
    s_b = (T_b < mat.T_rigid).astype(float)           # drawn metal follows the base only where the draw is not capped
    drawn_depth = np.zeros(sk.n)                      # m, how far down the 3D profile this pass has drawn
    melting = remove and mat.melts
    for _ in range(N):
        Q_side = side(sk, h) if side is not None else None
        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, Q_side)
        book.E_top += res[0]
        book.E_rad += res[1]
        book.E_base += res[2]
        book.E_side += res[3]
        if trial:
            dEb += res[4]
        if melting:
            k = np.cumprod(sk.valid & (sk.T >= mat.T_feed), axis=1).sum(axis=1)
            gone, gone_H = sk.take_top(k)
            book.removed += gone
            book.removed_H += gone_H
        short = np.maximum(m_min - sk.column_mass(), 0.0)
        if g_base is not None and np.any(short > 0.0):
            rows = np.flatnonzero(short > 0.0)
            dd = short[rows] / (mat.rho * sk.area[rows])
            h_b, s_b = h_b.copy(), s_b.copy()
            h_b[rows] = profile_enthalpy(mat, T_b[rows], g_base[rows], drawn_depth[rows], drawn_depth[rows] + dd,
                                         cap=mat.T_rigid)
            s_b[rows] = (T_b[rows] - g_base[rows] * (drawn_depth[rows] + 0.5 * dd) < mat.T_rigid).astype(float)
            drawn_depth[rows] += dd
        book.drawn += short
        book.drawn_H += sk.add_bottom(short, h_b, s_b)
        sk.regrid()
    if trial:
        sk.set_state(saved)
        b = np.minimum(dEb / dt, 0.0)
        return book.E_base / dt - b * T_b, b
    sk.S = None
    return book


# -- heat sideways (Task 5) --------------------------------------------------------------------------------------------
_SPAN = 1.0                       # m, larger than any column: separates the rows of a flattened row-wise search


def _interp_rows(x, xp, fp):
    """Row-wise linear interpolation of fp (n, q) given at non-decreasing xp (n, q), at x (n, p), held at the ends."""
    n, q = xp.shape
    off = (np.arange(n) * _SPAN)[:, None]
    i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
    i = np.clip(i - (np.arange(n) * q)[:, None], 0, q - 2)
    r = np.arange(n)[:, None]
    x0, x1, f0, f1 = xp[r, i], xp[r, i + 1], fp[r, i], fp[r, i + 1]
    t = np.clip(np.divide(x - x0, x1 - x0, out=np.zeros_like(x, dtype=float), where=x1 > x0), 0.0, 1.0)
    return f0 + t * (f1 - f0)


def depth_grid(settings, depth):
    """The common depth grid of the sideways exchange: breakpoints from the top down to `depth` [m], spaced by the cell
    target at each depth (the fine cell, growing below the fine zone)."""
    z = [0.0]
    while z[-1] < depth:
        z.append(z[-1] + settings.cell + (settings.growth - 1.0) * max(z[-1] - settings.fine_zone, 0.0))
    return np.array(z)


def _sideways_once(sk, columns, T, dz):
    """Sideways heat [W per cell] at temperatures T: each column's cells averaged onto the depth grid, every surface
    edge exchanging k w l (T_j - T_i) in every bin both columns reach (k the mean of the two bins' conductivities, w the
    edge's cotangent weight, l the depth the two columns share in the bin), and each column's bin heat given back to its
    cells in proportion to their overlap with the bin. Antisymmetric edge by edge, so the total is zero to round-off;
    exact for a field linear along the surface where neighbouring columns have the same cells."""
    n = T.shape[0]
    zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)          # cell boundaries, from the top
    H = zb[:, -1]
    zg = depth_grid(sk.settings, float(H.max()))
    zc = np.minimum(zg[None, :], H[:, None])                                        # the grid clipped to each column
    C = np.concatenate([np.zeros((n, 1)), np.cumsum(T * dz, axis=1)], axis=1)      # integral of T from the top
    Cg = _interp_rows(zc, zb, C)
    length = np.diff(zc, axis=1)                                                    # each column's share of each bin
    Tbar = np.divide(np.diff(Cg, axis=1), length, out=np.zeros_like(length), where=length > 0.0)
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    k = sk.props.k(Tbar)
    flux = (0.5 * (k[i] + k[j]) * columns.conduct_weight[:, None] * np.minimum(length[i], length[j])
            * (Tbar[j] - Tbar[i]))
    Qg = np.zeros_like(length)
    np.add.at(Qg, i, flux)
    np.add.at(Qg, j, -flux)
    Qcum = np.concatenate([np.zeros((n, 1)), np.cumsum(Qg, axis=1)], axis=1)
    return np.diff(_interp_rows(zb, zc, Qcum), axis=1)


def sideways_heat(sk, columns, dt):
    """The sideways conduction of one sub-step dt [s] as an explicit heat source [W per cell], the time average over
    as many sub-cycles as keep dt k_max sum(w) / (rho c_p,min A) below SIDEWAYS_CFL (about 25 ms per cycle at 2 mm
    spacing, spec 6.3); the cells' temperatures move between sub-cycles, the source being what `_conduct` then adds."""
    mat = sk.material
    valid = sk.valid
    dz = sk.dz()
    wsum = np.bincount(columns.edges.ravel(), weights=np.repeat(columns.conduct_weight, 2), minlength=sk.n)
    rate = float((mat.k_table.max() * wsum / (mat.rho * mat.cp_table.min() * sk.area)).max())
    cycles = max(1, int(math.ceil(dt * rate / SIDEWAYS_CFL)))
    T, H = sk.T.copy(), sk.H.copy()
    total = np.zeros_like(T)
    for _ in range(cycles):
        Q = np.where(valid, _sideways_once(sk, columns, T, dz), 0.0)
        total += Q / cycles
        if cycles > 1:
            H = H + Q * dt / cycles
            T = np.where(valid, mat.temperature_from_enthalpy(np.divide(H, sk.m, out=np.zeros_like(H), where=sk.m > 0.0)), T)
    return total


# -- the base (Task 7) -------------------------------------------------------------------------------------------------
def hand_down(sk, m_min, T_rigid):
    """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
    oldest first, as far as the column stays at or above its minimum mass m_min (spec 6.6); a cell may go in part.
    Returns (mass, enthalpy) handed down per column."""
    n, K = sk.m.shape
    j = np.arange(K)[None, :]
    idx = np.clip(sk.count[:, None] - 1 - j, 0, K - 1)                 # the cells from the bottom up
    rows = np.arange(n)[:, None]
    rigid = (j < sk.count[:, None]) & (sk.T[rows, idx] <= T_rigid)
    run = np.cumprod(rigid, axis=1).astype(bool)
    available = np.where(run, sk.m[rows, idx], 0.0).sum(axis=1)
    take = np.clip(np.minimum(available, sk.column_mass() - m_min), 0.0, None)
    return take, sk.take_bottom(take)
