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
import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

from . import rheology
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
KRYLOV_TOL = 1.0e-3               # each Newton step's coupled solve (sideways conduction in the implicit sub-step) stops
                                  # at this residual relative to its right-hand side: the Newton iteration corrects the rest
KRYLOV_MAX = 200                  # conjugate-gradient iterations of one coupled solve at most
KRYLOV = {"solves": 0, "iterations": 0}   # the coupled solves made and their iterations, for measuring the cost
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
    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b). `U` (n, K, 3) is every
    cell's velocity [m/s] (sub-plan 21 Task 2): part of the cells' state, carried by every operation with its momentum
    exact -- split halves keep it, merged cells take the mass-weighted mean, metal drawn from the 3D model arrives at
    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top`, `take_top_mass` or `take_bottom` took away
    with the metal.
    `flow_memory` (3,) is what the adaptive count of flow sub-steps remembers of the last heat sub-step's flow (sub-plan
    21 Task 6d, `flow_cycles`): its largest Courant number and Rayleigh-Taylor growth and its flow sub-step's length [s],
    zero before any flow; part of the state, so that a trial pass, a repeated or a shortened part starts from it."""

    def __init__(self, material, area, settings=None):
        self.material = material
        self.settings = settings or SkinSettings()
        self.area = np.array(area, dtype=float)
        n = len(self.area)
        self.m, self.H, self.T = np.zeros((n, 1)), np.zeros((n, 1)), np.full((n, 1), 300.0)
        self.count = np.zeros(n, dtype=np.int64)
        self.S = None
        self.U = np.zeros((n, 1, 3))
        self.momentum_out = np.zeros((n, 3))
        self.flow_memory = np.zeros(3)
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

    def column_momentum(self):
        """Each column's momentum [kg m/s], (n, 3)."""
        return (self.m[..., None] * self.U).sum(axis=1)

    def momentum(self):
        return self.column_momentum().sum(axis=0)

    def state(self):
        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count, self.U, self.flow_memory))

    def rows(self, idx):
        """Skins holding copies of the columns idx alone, with this one's material, settings and property tables: for
        a solve that needs only some columns (sub-plan 21 Task 6d)."""
        sub = object.__new__(Skins)
        sub.material, sub.settings, sub.rho_l, sub.props = self.material, self.settings, self.rho_l, self.props
        sub.area, sub.m, sub.H, sub.T, sub.count, sub.U = (a[idx] for a in (self.area, self.m, self.H, self.T, self.count,
                                                                           self.U))
        sub.S = None if self.S is None else self.S[idx]
        sub.momentum_out, sub.flow_memory = np.zeros((len(sub.area), 3)), self.flow_memory.copy()
        return sub

    def set_state(self, state):
        """Restore a `state()`; one of the six arrays before Task 6d's memory starts the flow with none."""
        self.area, self.m, self.H, self.T, self.count, self.U = (a.copy() for a in state[:6])
        self.flow_memory = np.array(state[6], dtype=float) if len(state) > 6 else np.zeros(3)
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
    def _install(self, m, H, count, S=None, T=None, U=None):
        """Install padded cell arrays. Temperatures follow from the enthalpies, except where `T` is given for cells
        whose specific enthalpy the operation did not change (a shifted, split or partly taken cell); nan in T marks
        the cells to recompute. `U` is the cells' velocity (at rest if not given)."""
        valid = np.arange(m.shape[1])[None, :] < count[:, None]
        self.m, self.H, self.count = np.where(valid, m, 0.0), np.where(valid, H, 0.0), np.asarray(count, dtype=np.int64)
        redo = valid if T is None else valid & np.isnan(T)
        self.T = np.full(self.m.shape, 300.0) if T is None else np.where(valid & ~redo, T, 300.0)
        if redo.any():
            self.T[redo] = self.material.temperature_from_enthalpy(self.H[redo] / self.m[redo])
        if self.S is not None:
            self.S = np.where(valid, S, 0.0)
        self.U = np.zeros(self.m.shape + (3,)) if U is None else np.where(valid[..., None], U, 0.0)

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
            self.U = np.concatenate([self.U, np.zeros((self.n, extra, 3))], axis=1)

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

    def _put_rows(self, rows, m, H, T, S, count, U):
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
        U = np.concatenate([U, np.zeros((len(rows), K - U.shape[1], 3))], axis=1) if U.shape[1] < K else U
        self.U[rows] = np.where(valid[..., None], U, 0.0)
        width = int(self.count.max()) if self.n else 1
        if width < K:
            self.m, self.H, self.T = self.m[:, :width], self.H[:, :width], self.T[:, :width]
            if self.S is not None:
                self.S = self.S[:, :width]
            self.U = self.U[:, :width]

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
        mU = np.zeros((n * K, 3))
        np.add.at(mU, flat, (m0[..., None] * self.U[rows])[valid])
        U = np.divide(mU, m[:, None], out=np.zeros_like(mU), where=m[:, None] > 0.0).reshape(n, K, 3)
        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups, U)

    def _split(self, rows, big):
        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T, S and velocity)."""
        valid = self.valid[rows]
        reps = np.where(valid, 1 + big[rows], 0)
        n, K = reps.shape
        count = reps.sum(axis=1)
        K2 = int(count.max())
        pos = np.cumsum(reps, axis=1) - reps
        r_all = np.repeat(np.arange(n)[:, None], K, axis=1)
        m, H, T = np.zeros((n, K2)), np.zeros((n, K2)), np.zeros((n, K2))
        S = np.zeros((n, K2)) if self.S is not None else None
        U = np.zeros((n, K2, 3))
        share = 1.0 / reps.clip(1)
        for part in (0, 1):
            sel = valid & (part < reps)
            r, c = r_all[sel], pos[sel] + part
            m[r, c] = (self.m[rows] * share)[sel]
            H[r, c] = (self.H[rows] * share)[sel]
            T[r, c] = self.T[rows][sel]                                    # halves keep the specific enthalpy
            if S is not None:
                S[r, c] = self.S[rows][sel]
            U[r, c] = self.U[rows][sel]
        self._put_rows(rows, m, H, T, S, count, U)

    # -- the top and the bottom -------------------------------------------------------------------------------------
    def take_top(self, k):
        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column, and leaves the momentum
        that went with them in `momentum_out`."""
        k = np.asarray(k, dtype=np.int64)
        n, K = self.m.shape
        top = np.arange(K)[None, :] < k[:, None]
        mass, H = (self.m * top).sum(axis=1), (self.H * top).sum(axis=1)
        self.momentum_out = ((self.m * top)[..., None] * self.U).sum(axis=1)
        if k.any():
            src = np.minimum(np.arange(K)[None, :] + k[:, None], K - 1)
            shift = lambda a: np.take_along_axis(a, src, axis=1)
            S = shift(self.S) if self.S is not None else None
            U = np.take_along_axis(self.U, src[..., None], axis=1)
            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T), U)
        return mass, H

    def take_top_mass(self, dm, admit, phase=None):
        """Remove dm >= 0 [kg] per column from its top down through the run of admissible cells from its top (`admit`,
        (n, K) bool), at most what that run holds, each cell's metal at its own specific enthalpy and velocity: a partly
        taken cell keeps both, its temperature and its sensitivity, and the cells beneath are untouched (sub-plan 21 Task
        7: the release of spec 8, which thins the top cell until it merges). Returns (mass, enthalpy, solid) taken per
        column, solid the mass of solid in it by the cells' liquid fractions, and leaves the momentum that went in
        `momentum_out`. With `phase` ((n, K) bool), what was taken is also split by the cells it came from, those in
        `phase` and the others, and left in `phase_out` (sub-plan 22 Task 2: each spray mode's droplets)."""
        dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
        if not np.all(np.isfinite(dm) & (dm >= 0.0)):
            raise ValueError("a release from a skin's top must be a finite, non-negative mass")
        run = np.cumprod(self.valid & admit, axis=1).astype(bool)
        above = np.concatenate([np.zeros((self.n, 1)), np.cumsum(self.m, axis=1)[:, :-1]], axis=1)   # mass above each
        taken = np.where(run, np.clip(dm[:, None] - above, 0.0, self.m), 0.0)
        frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
        mass, H = taken.sum(axis=1), (self.H * frac).sum(axis=1)
        f_s = 1.0 - self.props.liquid_fraction(self.T)
        solid = (taken * f_s).sum(axis=1)
        momentum = (taken[..., None] * self.U).sum(axis=1)
        if phase is not None:
            self.phase_out = self._by_phase(taken, frac, f_s, np.asarray(phase, dtype=bool), (mass, H, solid, momentum))
        whole = np.cumprod(run & (taken >= self.m), axis=1).sum(axis=1)        # the cells wholly taken, from the top
        part = (taken > 0.0) & (np.arange(self.m.shape[1])[None, :] >= whole[:, None])
        self.m, self.H = np.where(part, self.m - taken, self.m), np.where(part, self.H * (1.0 - frac), self.H)
        self.take_top(whole)
        self.momentum_out = momentum
        return mass, H, solid

    def _by_phase(self, taken, frac, f_s, phase, totals):
        """`take_top_mass`'s (mass, enthalpy, solid, momentum) taken, each with a leading axis of two: [0] from the cells
        in `phase`, [1] from the others. A column that took from cells of one kind gives that kind its totals to the last
        bit and the other nothing; one that took from both sums each kind's cells."""
        got = taken > 0.0
        inside, outside = (got & phase).any(axis=1), (got & ~phase).any(axis=1)
        out = [np.zeros((2,) + x.shape) for x in totals]
        for j, only in ((0, inside & ~outside), (1, outside & ~inside)):
            for o, x in zip(out, totals):
                o[j, only] = x[only]
        both = np.flatnonzero(inside & outside)
        for j, on in enumerate((phase[both], ~phase[both])):
            t = taken[both] * on
            parts = (t.sum(axis=1), (self.H[both] * frac[both] * on).sum(axis=1), (t * f_s[both]).sum(axis=1),
                     (t[..., None] * self.U[both]).sum(axis=1))
            for o, x in zip(out, parts):
                o[j, both] = x
        return tuple(out)

    def add_bottom(self, dm, h, s=1.0, u=None):
        """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
        a new cell); returns the enthalpy added. With S tracked, the metal carries the sensitivity s per column
        (default 1: it is at the base temperature and follows it; 0 for metal whose enthalpy does not depend on the
        base temperature, such as rigid metal drawn at T_rigid's enthalpy from a base above T_rigid). It arrives with
        the velocity u (n, 3) [m/s], at rest by default: metal from the 3D model is rigid."""
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
        u = np.zeros((self.n, 3)) if u is None else np.broadcast_to(np.asarray(u, dtype=float), (self.n, 3))
        self.U[rows, b] = (m_old[:, None] * self.U[rows, b] + dm[rows, None] * u[rows]) / self.m[rows, b][:, None]
        self.count[rows] = np.maximum(self.count[rows], 1)
        self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])
        return np.where(some, dm * h, 0.0)

    def take_bottom(self, dm):
        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy and velocity
        (a partly taken cell keeps both); returns the enthalpy removed per column, and leaves the momentum that went in
        `momentum_out`."""
        dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
        if np.any(dm > self.column_mass() * (1.0 + 1e-12) + 1e-300):
            raise ValueError("cannot take more metal from a skin's bottom than the skin holds")
        self.momentum_out = np.zeros((self.n, 3))
        if not np.any(dm > 0.0):
            return np.zeros(self.n)
        n, K = self.m.shape
        valid = self.valid
        below = np.cumsum(self.m[:, ::-1], axis=1)[:, ::-1] - self.m      # mass beneath each cell
        taken = np.clip(dm[:, None] - below, 0.0, self.m)                 # how much of each cell goes
        frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
        H_taken = (self.H * frac).sum(axis=1)
        self.momentum_out = (taken[..., None] * self.U).sum(axis=1)
        m, H = self.m - taken, self.H * (1.0 - frac)
        keep = valid & (m > 0.0)
        count = keep.sum(axis=1)                                          # the kept cells are the top ones
        self._install(m, H, count, self.S, self.T, self.U)                # a partly taken cell keeps its T and velocity
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
    momentum_out: np.ndarray = None   # kg m/s (n, 3): what the removed liquid took with it (sub-plan 21 Task 6c)
    flow: "FlowBook" = None           # the flow sub-steps' books, with the flow on (sub-plan 21 Task 6c)


def _coupled_solve(G, diag, rhs, valid, sideways):
    """One Newton step of the coupled sub-step: J x = rhs with J the columns' tridiagonal Jacobians less the sideways
    map (`Sideways.heat`, linear), symmetric positive definite, by conjugate gradients preconditioned by the columns' own
    tridiagonal solves (`thomas`) and started from them, to KRYLOV_TOL of the right-hand side. Returns x and the sideways
    heat at x, built from the products the iteration made (the map is linear), which the line search reuses."""
    def depth(x):
        y = diag * x
        y[:, :-1] -= G * x[:, 1:]
        y[:, 1:] -= G * x[:, :-1]
        return y

    def solve_columns(r):
        return np.where(valid, thomas(-G, diag, -G, r), 0.0)

    x = solve_columns(rhs)
    Hx = sideways.heat(x)
    r = np.where(valid, rhs - depth(x) + Hx, 0.0)
    target = KRYLOV_TOL * math.sqrt(float((rhs * rhs).sum()))
    z = solve_columns(r)
    p = z.copy()
    rz = float((r * z).sum())
    KRYLOV["solves"] += 1
    for it in range(KRYLOV_MAX + 1):
        if math.sqrt(float((r * r).sum())) <= target:
            KRYLOV["iterations"] += it
            return x, Hx
        Hp = sideways.heat(p)
        Jp = np.where(valid, depth(p) - Hp, 0.0)
        a = rz / float((p * Jp).sum())
        x, Hx, r = x + a * p, Hx + a * Hp, r - a * Jp
        z = solve_columns(r)
        rz, rz_old = float((r * z).sum()), rz
        p = z + (rz / rz_old) * p
    raise RuntimeError("a skin sub-step's coupled sideways solve did not converge in {} iterations".format(KRYLOV_MAX))


def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, sideways=None):
    """One implicit (backward-Euler) sub-step for every column: the convective flux q_top [W/m2] over the top area
    A_top in at the top cell, radiation from it over the same area, conduction through the cells (harmonic
    conductances), the sideways conduction `sideways` (a `Sideways` map made at the sub-step's start, implicit with the
    rest: sub-plan 21, design (a)), and the base held at T_b [K] through the bottom half-cell (T_b None: an adiabatic
    base, for the tests). Newton in the temperature with a backtracking line search per column (stage 1's solver,
    sub-plan 18); with the sideways map the columns are coupled and each Newton step is solved by `_coupled_solve`.
    Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the skins track S, updates it to dT/dT_b and returns
    the base energy's derivative with respect to T_b as a fifth item -- through the depth only: the base law each column
    gives the 3D model is its own, and what the neighbours' share misses is the interface mismatch's to book. E_base is
    what the books leave, so they close exactly whatever the iteration error."""
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
    side = (lambda T: np.zeros((n, K))) if sideways is None else sideways.heat

    def system(T, Q=None):
        """The residual, the conductances and the Jacobian's diagonal at T, and the sideways heat Q at T (given when the
        caller has it already)."""
        Q = side(T) if Q is None else Q
        h, cp, k = mat.conduction(T)
        G = np.where(link, A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:]), 0.0)
        Gb = A * k[rows, last] / (0.5 * dz[rows, last]) if base else np.zeros(n)
        flux = G * (T[:, :-1] - T[:, 1:])
        net = Q.copy()
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
        return R, G, Gb, np.where(valid, diag, 1.0), Q

    # Newton in the temperature with a backtracking line search per column on the residual's norm: in 10-micrometre
    # cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step, so a cell must move
    # with its neighbours, and the 3D solver's per-node enthalpy-consistent update stalled here (sub-plan 18)
    T = T_old.copy()
    R, G, Gb, diag, Q = system(T)
    moved = np.zeros(n)
    for _ in range(NEWTON_MAX):
        if sideways is None:
            dT, Q_dT = thomas(-G, diag, -G, -R), None
        else:
            dT, Q_dT = _coupled_solve(G, diag, -R, valid, sideways)
        r0 = np.sqrt((R * R).sum(axis=1))
        lam = np.ones(n)
        for _ in range(LINE_SEARCH_MAX):
            T_try = T + lam[:, None] * dT
            full = Q_dT is not None and bool(np.all(lam == 1.0))         # the map is linear: Q(T + dT) = Q(T) + Q(dT)
            R_try, G_try, Gb_try, diag_try, Q_try = system(T_try, Q + Q_dT if full else None)
            moved = np.abs(T_try - T).max(axis=1)
            worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
            if not worse.any():
                break
            lam = np.where(worse, 0.5 * lam, lam)
        T, R, G, Gb, diag, Q = T_try, R_try, G_try, Gb_try, diag_try, Q_try
        if moved.max() <= NEWTON_TOL:
            break
    else:
        raise RuntimeError("a skin sub-step did not converge in {} iterations (last change {:.2e} K)".format(
            NEWTON_MAX, moved.max()))
    T = np.where(valid, T, T_old)
    E_new = m * mat.enthalpy(T)
    E_top = q * dt
    E_rad = eps_s * (T[:, 0] ** 4 - T_amb ** 4) * dt
    E_side = Q.sum(axis=1) * dt
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


def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None, flow=None):
    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the implicit conduction, through the depth
    and, with `side`, sideways (`side(sk, previous)` makes the sub-step's `Sideways` map from the last one, Task 5),
    then, with `flow` (a `Flow`), the flow's sub-steps (`flow_substeps`: heat, then flow, sideways exchange and spray,
    spec 7.1; sub-plan 21 Task 6c, the spray a release from the top cells with `flow.release` since Task 7), then --
    with `remove` and a melting material -- the contiguous top cells at or above T_feed leave, taking their momentum
    with them, a column below its minimum mass m_min [kg] draws
    back to it, and the cells are split and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but
    never above T_rigid -- a 3D base node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it
    reached 930 K in one 0.5 s step, found while planning), and drawing at its temperature would draw liquid that the
    next sub-step removes, without limit; what the 3D model actually loses is reconciled by the body after the base
    moves. With `g_base` [K/m per column], the
    3D model's temperature gradient along the outward normal below the base, each draw takes the next layer down: its
    enthalpy is h(min(T, T_rigid)) averaged along T(s) = T_b - g_base s over the depth it takes (`profile_enthalpy`).

    T_b1 None: the trial pass, on a copy (the skins are left as they were), base held at T_b0; returns the base law
    (a, b) per column, E_base / dt = a + b T_b [W into the 3D model], from the pass's base energy and its derivative with
    respect to the base temperature, the slope clipped at zero (stage 1's choice: a law rising with T_b would cost the 3D
    matrix its diagonal dominance; what the clip misses is the interface mismatch's to book). Otherwise the real pass,
    base held at T_b1, the 3D model's end-of-step base temperature; returns a SkinBook. The trial pass runs the flow too
    (spec 7.2's sub-steps are spec 7.1's: heat, flow, sideways exchange and spray), the sensitivity S carried with the
    metal the transport moves. Where the flow's count is chosen (`flow.cycles` None, Task 6d), the trial pass leaves the
    counts it took in `flow.plan` (0 where nothing flowed, so it took none), and the real pass after it takes, in each
    heat sub-step, the larger of the trial's count and the count its own memory chooses, and its own where the trial
    took none (spec 16, 2026-10-26; sub-plan 21 Task 9): so the base law is made by the pass it predicts wherever the
    two flow alike (6d), while a real pass whose base melts more than its trial pass's -- held at the start-of-step
    temperature, a trial can flow slower at a melting onset -- runs at the Courant target rather than at the trial's
    counts (on Task 8's pushed plate one flow sub-step a heat sub-step at a Courant number of 18.8), and does not replay
    counts chosen where nothing flowed (380 flow sub-steps there in which nothing moved). A real pass with no trial
    before it chooses its own."""
    mat = sk.material
    N = max(1, int(math.ceil(dt / sk.settings.substep - 1e-9)))
    h = dt / N
    trial = T_b1 is None
    saved = sk.state() if trial else None
    T_b = np.asarray(T_b0 if trial else T_b1, dtype=float) * np.ones(sk.n)
    book = SkinBook(*(np.zeros(sk.n) for _ in range(8)), momentum_out=np.zeros((sk.n, 3)),
                    flow=FlowBook.zero(sk.n) if flow is not None else None)
    dEb = np.zeros(sk.n)
    if trial:
        sk.S = np.zeros_like(sk.m)
    h_b = mat.enthalpy(np.minimum(T_b, mat.T_rigid))
    s_b = (T_b < mat.T_rigid).astype(float)           # drawn metal follows the base only where the draw is not capped
    drawn_depth = np.zeros(sk.n)                      # m, how far down the 3D profile this pass has drawn
    melting = remove and mat.melts
    sideways = None
    plan = replay = None
    if flow is not None and flow.cycles is None:                  # the trial pass's counts, for its real pass
        replay = not trial and flow.plan is not None and len(flow.plan) == N
        plan = list(flow.plan) if replay else []
    for i in range(N):
        sideways = side(sk, sideways) if side is not None else None
        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, sideways)
        book.E_top += res[0]
        book.E_rad += res[1]
        book.E_base += res[2]
        book.E_side += res[3]
        if trial:
            dEb += res[4]
        if flow is not None:
            cycles = None
            if plan is not None:                                  # chosen: the trial's (or the real's own) count
                cycles = flow_cycles(sk.flow_memory, h)
                if replay and plan[i]:
                    cycles = max(plan[i], cycles)
            fb = flow_substeps(sk, flow, h, cycles=cycles)
            if plan is not None and not replay:
                plan.append(fb.cycles)                           # 0 where nothing flowed: the real pass takes its own
            book.flow.add(fb)
        if melting:
            k = np.cumprod(sk.valid & (sk.T >= mat.T_feed), axis=1).sum(axis=1)
            gone, gone_H = sk.take_top(k)
            book.removed += gone
            book.removed_H += gone_H
            book.momentum_out += sk.momentum_out
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
    if plan is not None:
        flow.plan = plan if trial else None
    if trial:
        sk.set_state(saved)
        b = np.minimum(dEb / dt, 0.0)
        return book.E_base / dt - b * T_b, b
    sk.S = None
    return book


# -- heat sideways (Task 5) --------------------------------------------------------------------------------------------
_SPAN = 1.0                       # m, larger than any column: separates the rows of a flattened row-wise search


def _interp_weights(x, xp):
    """Row-wise linear interpolation at x (n, p) of a function given at non-decreasing xp (n, q), held at the ends, as
    (rows, index, fraction) for `_interp_apply`: they depend on the abscissae only, so a sub-step's sideways map makes
    them once and applies them at every iteration of its solve."""
    n, q = xp.shape
    off = (np.arange(n) * _SPAN)[:, None]
    i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
    i = np.clip(i - (np.arange(n) * q)[:, None], 0, q - 2)
    r = np.arange(n)[:, None]
    x0, x1 = xp[r, i], xp[r, i + 1]
    t = np.clip(np.divide(x - x0, x1 - x0, out=np.zeros_like(x, dtype=float), where=x1 > x0), 0.0, 1.0)
    return r, i, t


def _interp_apply(weights, fp):
    """fp (n, q), given at the abscissae `weights` were made from, interpolated at their points."""
    r, i, t = weights
    f0, f1 = fp[r, i], fp[r, i + 1]
    return f0 + t * (f1 - f0)


def depth_grid(settings, depth):
    """The common depth grid of the sideways exchange: breakpoints from the top down to `depth` [m], spaced by the cell
    target at each depth (the fine cell, growing below the fine zone)."""
    z = [0.0]
    while z[-1] < depth:
        z.append(z[-1] + settings.cell + (settings.growth - 1.0) * max(z[-1] - settings.fine_zone, 0.0))
    return np.array(z)


def _sideways_geometry(sk, columns, dz):
    """What a sub-step's sideways map is built on, from the cells' thicknesses dz alone: the cells' boundaries from the
    top, each column's share of each bin of the depth grid, the depth each edge's two columns share in each bin, and the
    interpolation weights from cells to bins and back."""
    n = dz.shape[0]
    zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)          # cell boundaries, from the top
    H = zb[:, -1]
    zg = depth_grid(sk.settings, float(H.max()))
    zc = np.minimum(zg[None, :], H[:, None])                                        # the grid clipped to each column
    length = np.diff(zc, axis=1)                                                    # each column's share of each bin
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    return {"length": length, "shared": np.minimum(length[i], length[j]),
            "to_bins": _interp_weights(zc, zb), "to_cells": _interp_weights(zb, zc)}


class Sideways:
    """The sideways conduction of one sub-step as a linear map on the cells' temperatures (spec 6.3; implicit since
    sub-plan 21, design (a), Asha 2026-10-18). Each column's cells are averaged onto the common depth grid, every surface
    edge exchanges c (x_j - x_i) in every bin both columns reach, with c = k w l (k the mean of the two bins'
    conductivities at the sub-step's start, w the edge's cotangent weight, l the depth the two columns share in the bin),
    and each column's bin heat goes back to its cells in proportion to their overlap with the bin. `heat(x)` [W per cell]
    is antisymmetric edge by edge, so it sums to zero to round-off, and exact for a field linear along the surface where
    neighbouring columns have the same cells. As a map it is symmetric and negative semi-definite -- the bins' averaging
    and the heat's return to the cells are each other's transposes -- so `_conduct` takes it into its implicit solve
    (`_coupled_solve`), and the sub-step has no sideways stability limit: the explicit exchange it replaces needed three
    sub-cycles a 10 ms sub-step on the 100 mm sphere's smallest columns (facts 147-148)."""

    def __init__(self, sk, columns, previous=None):
        """`previous`: the last sub-step's map, whose geometry is reused while the cells keep their thicknesses (in a
        body that does not melt, all of its sub-steps)."""
        self.valid, self.dz = sk.valid, sk.dz()
        same = previous is not None and previous.dz.shape == self.dz.shape and np.array_equal(previous.dz, self.dz)
        self.geometry = previous.geometry if same else _sideways_geometry(sk, columns, self.dz)
        self.i, self.j = columns.edges[:, 0], columns.edges[:, 1]
        k = sk.props.k(self.bins(sk.T))
        self.c = 0.5 * (k[self.i] + k[self.j]) * columns.conduct_weight[:, None] * self.geometry["shared"]
        self.B = columns.incidence

    def bins(self, x):
        """The mean of x (per cell) over each bin of the depth grid, (n, bins); zero in bins beyond a column."""
        n = x.shape[0]
        g = self.geometry
        C = np.concatenate([np.zeros((n, 1)), np.cumsum(x * self.dz, axis=1)], axis=1)      # integral from the top
        return np.divide(np.diff(_interp_apply(g["to_bins"], C), axis=1), g["length"], out=np.zeros_like(g["length"]),
                         where=g["length"] > 0.0)

    def heat(self, x):
        """The heat [W per cell] the neighbouring columns conduct into each cell at cell temperatures x (linear in x)."""
        xb = self.bins(x)
        Qg = self.B @ (self.c * (xb[self.j] - xb[self.i]))
        Qcum = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(Qg, axis=1)], axis=1)
        return np.where(self.valid, np.diff(_interp_apply(self.geometry["to_cells"], Qcum), axis=1), 0.0)


# -- viscosity (sub-plan 21 Task 3) ------------------------------------------------------------------------------------
def cell_viscosity(sk, gamma):
    """Every cell's viscosity [Pa s], (n, K), sheared at gamma [1/s] (per cell, or one for all): `rheology.viscosity` at
    the cell's temperature and Scheil solid fraction; infinite in a rigid cell (half solid or more) and beyond each
    column's cells, and everywhere for a material that does not melt."""
    mat = sk.material
    if not mat.melts or mat.liquid is None:
        return np.full(sk.m.shape, np.inf)
    f_s = 1.0 - sk.props.liquid_fraction(sk.T)
    f_s_bridge = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    eta = rheology.viscosity(sk.T, f_s, gamma, mat.liquid.mu, f_s_bridge)
    return np.where(sk.valid, eta, np.inf)


def series_drag(sk, eta):
    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1), from the
    cells' viscosities eta [Pa s] (n, K): 1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), the two half-cells in series, as
    the shear stress's continuity across the interface requires, the way `_conduct` combines conductivities. A rigid
    cell resists without limit and adds nothing to the sum, so a flowing cell drags on the rigid cell beneath it over
    its own half-thickness alone: no slip at the rigid front. Zero between two rigid cells and across a missing one."""
    if sk.m.shape[1] < 2:
        return np.zeros((sk.n, 0))
    dz = sk.dz()
    resist = 0.5 * dz[:, :-1] / eta[:, :-1] + 0.5 * dz[:, 1:] / eta[:, 1:]
    link = sk.valid[:, 1:] & (np.isfinite(eta[:, :-1]) | np.isfinite(eta[:, 1:]))
    return np.where(link, 1.0 / np.where(link, resist, 1.0), 0.0)


def interface_drag(sk, gamma):
    """`series_drag` with every cell's viscosity at its own shear rate gamma [1/s] (per cell, or one for all). Each
    cell, not each interface: across a jump from liquid to slurry the liquid shears some thousand times faster than the
    slurry, and a rate taken from the velocity difference between the two centres would give the slurry's half-cell the
    liquid's rate (found while writing Task 4)."""
    return series_drag(sk, cell_viscosity(sk, gamma))


# -- momentum through the depth (sub-plan 21 Task 4) --------------------------------------------------------------------
SHEAR_PASSES = 5                  # passes of the fixed point on the cells' shear rates at most (spec 6.4)
SHEAR_TOL = 1.0e-3                # a pass that moves no flowing cell's viscosity by more than this, relatively, ends it
CRUST_VISCOSITY = 1.0e9           # Pa s: a floating crust's half-cell at its interface with flowing metal (some 1e8 times
                                  # the slurry's, adding nothing to the drag); the crust itself is one unknown (Task 6b)


@dataclass
class MomentumBook:
    """What a momentum sub-step did, per column: impulses [N s], (n, 3), which add up to the momentum change."""
    push: np.ndarray         # of the push on the moving cells
    shear: np.ndarray        # of the gas shear on the top cell, if it moves
    wall: np.ndarray         # of the anchored rigid metal -- the skin's rigid cells joined to the substrate and the 3D
                             # model beneath the base -- on the moving cells, and the momentum a cell loses on joining it
    passes: int              # the fixed point's passes
    response: np.ndarray = None       # with response=True, (n, K) [m/s per N/m3]: every cell's velocity per unit push,
                                      # the last pass's system solved for the push alone (zero where a cell is held)
    wall_response: np.ndarray = None  # with response=True, (n,) [m3]: the anchored metal's drag on that response, so
                                      # that a change dG of the push over the sub-step gives the wall -dt wall_response dG


def anchored_cells(sk, flowing):
    """The rigid cells joined to the substrate, (n, K): each column's run of rigid cells from its bottom cell up, which
    rests on the rigid 3D model. A column whose bottom cell flows has none."""
    n, K = sk.m.shape
    j = np.arange(K)[None, :]
    from_bottom = np.clip(sk.count[:, None] - 1 - j, 0, K - 1)                 # the j-th cell from the bottom
    rigid_up = np.take_along_axis(sk.valid & ~flowing, from_bottom, axis=1) & (j < sk.count[:, None])
    held = np.cumprod(rigid_up, axis=1).sum(axis=1)
    return sk.valid & (j >= (sk.count - held)[:, None])


def momentum_step(sk, dt, push, shear, A_top, response=False):
    """One implicit sub-step dt [s] of every column's momentum through its depth (spec 6.4), each velocity component
    apart:

        m_i (u_i' - u_i) / dt = A [D_(i-1) (u'_(i-1) - u'_i) + D_i (u'_(i+1) - u'_i)] + G V_i + tau A_top [i the top],

    A the column's area, D the drag between neighbouring cells (`series_drag`), G the push per unit volume [N/m3,
    (n, 3)] -- every layer of a column feels the same (spec 6.4) -- V_i the cell's volume and tau the gas shear [Pa,
    (n, 3)] on the top cell over its top area A_top [m2, (n,)]. The rigid metal joined to the substrate (`anchored_cells`)
    is held at rest; a moving cell drags on it, and a flowing bottom cell on the rigid 3D metal beneath the base, with no
    slip over its own half-thickness. A crust -- rigid cells with flowing metal between them and the substrate -- floats:
    it moves as one piece, one unknown of the system (its mass the cells', its drags those at its top and bottom faces),
    dragged by the gas shear on its top, pushed like the rest and dragged by the slurry beneath (Asha, 2026-10-20,
    option B; the column model has no sideways stress in a solid, so a crust over many columns moves as column-sized
    plates, as a warm, barely solid crust would). Task 4 coupled the crust's cells by a viscosity of 1e9 Pa s, whose rows,
    some 1e12 times the cells' inertia, amplified the solver's round-off: on mixed sphere skins the books closed to
    1.2e-5 of the largest impulse, now to 6e-14 (1 ms) and 6e-13 (10 ms), and the velocities moved by up to 1.1e-3 of the
    largest -- the coupled crust's own round-off, which fell as 1.1e-4, 1.1e-5 and 9.6e-7 at 1e8, 1e7 and 1e6 Pa s;
    columns without a crust are unchanged to the last bit (sub-plan 21 Task 6b, measured 2026-10-21). Each flowing
    cell's viscosity is at its own shear rate, the mean over its two faces of the stress there over its viscosity, found
    by fixed-point iteration from the start-of-sub-step velocities: at most SHEAR_PASSES passes, ending when no flowing
    cell's viscosity moves by more than SHEAR_TOL. Updates sk.U; returns a MomentumBook whose impulses add up to the
    columns' momentum change to round-off. With `response`, the book also carries the velocities a unit push alone
    gives through the last pass's system and the wall's drag on them: with the viscosities held, the sub-step is linear
    in the push, so a change of the push changes every velocity by its response times the change (sub-plan 21 Task
    6b).

    Only the columns where something moves are solved (sub-plan 21 Task 6d): in the others every cell is anchored, so it
    ends at rest and the momentum it held is the wall's, as the solve over every column gave -- to the last bit, each
    column's system being its own (6c's step solved all 9 417 columns of the 100 mm sphere where half of them flowed)."""
    n = sk.n
    push = np.broadcast_to(np.asarray(push, dtype=float), (n, 3))
    shear = np.broadcast_to(np.asarray(shear, dtype=float), (n, 3))
    A_top = np.broadcast_to(np.asarray(A_top, dtype=float), (n,))
    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
    live = np.flatnonzero((sk.valid & ~anchored_cells(sk, flowing)).any(axis=1))
    if len(live) == n:
        return _momentum_columns(sk, dt, push, shear, A_top, response)
    book = MomentumBook(np.zeros((n, 3)), np.zeros((n, 3)), -sk.column_momentum(), 0,
                        np.zeros(sk.m.shape) if response else None, np.zeros(n) if response else None)
    U = np.zeros_like(sk.U)
    if len(live):
        sub = sk.rows(live)
        part = _momentum_columns(sub, dt, push[live], shear[live], A_top[live], response)
        U[live] = sub.U
        book.push[live], book.shear[live], book.wall[live], book.passes = part.push, part.shear, part.wall, part.passes
        if response:
            book.response[live], book.wall_response[live] = part.response, part.wall_response
    sk.U = U
    return book


def _momentum_columns(sk, dt, push, shear, A_top, response):
    """`momentum_step` on every column of sk, push (n, 3), shear (n, 3) and A_top (n,) given whole."""
    n, K = sk.m.shape
    valid = sk.valid
    rows = np.arange(n)
    last = np.maximum(sk.count - 1, 0)
    dz, A, m = sk.dz(), sk.area, sk.m
    flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
    anchored = anchored_cells(sk, flowing)
    moving = valid & ~anchored                                            # the flowing cells and the floating crust
    crust = moving & ~flowing
    held = ((m * anchored)[..., None] * sk.U).sum(axis=1)                 # what anchored cells carried, now at rest
    U_old = np.where(moving[..., None], sk.U, 0.0)
    top_moves = moving[:, 0]
    force = np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0)
    force[:, 0] += np.where(top_moves[:, None], shear * A_top[:, None], 0.0)
    bottom_flows = flowing[rows, last] & (sk.count > 0)
    # the unknowns: every cell, but a crust's cells one unknown together (sub-plan 21 Task 6b)
    same = np.zeros((n, K), dtype=bool)
    same[:, 1:] = crust[:, 1:] & crust[:, :-1]
    unit = np.cumsum(~same, axis=1) - 1
    Kr = int(unit.max()) + 1
    lump = lambda x: np.bincount((rows[:, None] * Kr + unit).ravel(), x.ravel(), n * Kr).reshape(n, Kr)
    moving_r, valid_r, last_r = lump(moving.astype(float)) > 0.0, lump(valid.astype(float)) > 0.0, unit[rows, last]
    cross_r, cross_k = np.nonzero(~same[:, 1:])                          # the interfaces between unknowns

    def viscosities(gamma):
        return np.where(crust, CRUST_VISCOSITY, cell_viscosity(sk, gamma))

    def drags(eta):
        Db = np.where(bottom_flows, eta[rows, last] / (0.5 * np.maximum(dz[rows, last], DZ_FLOOR)), 0.0)
        return series_drag(sk, eta), np.where(np.isfinite(Db), Db, 0.0)

    def system(D, Db):
        Dr = np.zeros((n, Kr - 1))
        Dr[cross_r, unit[cross_r, cross_k]] = D[cross_r, cross_k]
        diag = np.where(valid_r, lump(m) / dt, 1.0)
        diag[:, :-1] += A[:, None] * Dr
        diag[:, 1:] += A[:, None] * Dr
        diag[rows, last_r] += A * Db
        off = -A[:, None] * Dr
        sub, sup = np.where(moving_r[:, 1:], off, 0.0), np.where(moving_r[:, :-1], off, 0.0)
        return sub, np.where(moving_r, diag, 1.0), sup

    def unknowns(sub, diag, sup, rhs):
        """Each cell's value from the system on the unknowns, its right-hand side summed over each unknown's cells."""
        return np.where(moving, thomas(sub, diag, sup, lump(np.where(moving, rhs, 0.0)))[rows[:, None], unit], 0.0)

    def solve(eta):
        D, Db = drags(eta)
        sub, diag, sup = system(D, Db)
        rhs = (m / dt)[..., None] * U_old + force
        return np.stack([unknowns(sub, diag, sup, rhs[..., c]) for c in range(3)], axis=-1), D, Db

    def rates(U, D, Db, eta):
        """Every flowing cell's shear rate: the mean over its faces of the stress there over its viscosity."""
        stress = np.zeros((n, K + 1))                                     # at every face, the top face first
        stress[:, 1:K] = D * np.linalg.norm(U[:, :-1] - U[:, 1:], axis=2)
        stress[:, 0] = np.where(top_moves, np.linalg.norm(shear, axis=1), 0.0)
        stress[rows, last + 1] = Db * np.linalg.norm(U[rows, last], axis=1)
        face = 0.5 * (stress[:, :-1] + stress[:, 1:])
        return np.where(flowing, np.divide(face, eta, out=np.zeros_like(face), where=np.isfinite(eta)), 0.0)

    eta = viscosities(0.0)
    eta = viscosities(rates(U_old, *drags(eta), eta))                       # at the start-of-sub-step rates
    for passes in range(1, SHEAR_PASSES + 1):
        U, D, Db = solve(eta)
        eta_new = viscosities(rates(U, D, Db, eta))
        change = np.abs(eta_new[flowing] - eta[flowing]) / eta[flowing]
        if change.size == 0 or change.max() <= SHEAR_TOL:
            break
        eta = eta_new
    wall = np.zeros((n, 3))
    if K > 1:
        down = moving[:, :-1] & anchored[:, 1:]                           # a moving cell on the anchored metal beneath
        wall -= ((A[:, None] * D * down)[..., None] * U[:, :-1]).sum(axis=1)
    wall -= (A * Db)[:, None] * U[rows, last]
    book = MomentumBook(push=dt * np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0).sum(axis=1),
                        shear=dt * np.where(top_moves[:, None], shear * A_top[:, None], 0.0),
                        wall=dt * wall - held, passes=passes)
    if response:                                                          # the last pass's system, for a unit push alone
        r = unknowns(*system(D, Db), dz * A[:, None])
        book.response = r
        book.wall_response = A * Db * r[rows, last]
        if K > 1:
            book.wall_response += (A[:, None] * D * (moving[:, :-1] & anchored[:, 1:]) * r[:, :-1]).sum(axis=1)
    sk.U = U
    return book


# -- the push (sub-plan 21 Task 5) -------------------------------------------------------------------------------------
def along_meridian(normals, v_hat, magnitude):
    """Vectors of the given magnitudes (n,) in each corner's tangent plane, pointing along the meridian away from the
    stagnation point -- the direction of minus the flight direction's tangential part, the way the gas flows and the
    surface flow's s runs; zero where that part vanishes (the stagnation and rear points)."""
    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
    away = -(v[None, :] - (normals @ v)[:, None] * normals)
    size = np.linalg.norm(away, axis=1)
    unit = np.divide(away, size[:, None], out=np.zeros_like(away), where=size[:, None] > 1e-12)
    return np.asarray(magnitude, dtype=float)[:, None] * unit


def push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma, wet=None):
    """The push along the surface on every layer of every column [N/m3], (n, 3), and its four terms (spec 6.4), each in
    the column's own tangent plane -- the base's, where `momentum_step` solves for the velocities:

        pressure      -dp_w/ds along the meridian away from the stagnation point, dp_ds [Pa/m, (n,)] the surface flow's
                      (`SurfaceFlowResult.dp_ds`), differentiated on its own fine bins of the angle. The gradient of the
                      corners' own p_w would carry the scatter of their normals' angles, amplified some 25 times (the
                      radius over the spacing) by the differentiation: on the 2 mm sphere a Newtonian pressure's push
                      came out 1.1 % off on average and 13 % at worst that way, from normals 0.06 deg off radial on
                      average (measured 2026-10-20).
        deceleration  rho_l g_t, g [m/s2, (3,) or (n, 3)] the acceleration the metal feels in the body's frame, minus the
                      body's own -- for a body slowed by drag, the deceleration along the flight direction v_hat -- so
                      that it pulls the metal forward, toward the nose along the surface and out of the front face;
                      the spec's -rho_l a_t, and the element model's -rho_l W sin(theta) on the windward face.
        thickness     rho_l (g . n) grad_s H, H [m, (n,)] the columns' thicknesses: a thicker column under a pull out
                      of the surface has the lower pressure beneath it and draws metal in (the front face's
                      Rayleigh-Taylor mechanism).
        capillary     -sigma grad_s kappa, kappa the thin film's mean curvature (`Columns.film_curvature`): the base's
                      plus minus the cotangent Laplacian of the thicknesses on the top surface over each column's area
                      (sub-plan 21 Task 6b; Asha, 2026-10-21, option A) -- not the top surface's own curvature, whose
                      corners touching an obtuse triangle take their neighbours' mean, cannot see their own height, and
                      make surface tension draw metal onto a pile-up there. It acts only where a column's top flows
                      (`wet`, (n,) bool, `flowing_tops`; None for everywhere): not on an anchored rigid top nor on a
                      floating crust; a rigid neighbour's height enters a flowing column's curvature as a wall's would,
                      and in the gradient a wall takes its flowing neighbours' mean curvature (`_walled`).

    The two thickness terms are evaluated at the given H; `thickness_step` takes them implicitly in the column heights."""
    H = np.asarray(thickness, dtype=float)
    normals = columns.normal
    g = np.broadcast_to(np.asarray(g, dtype=float), (columns.n, 3))
    g_n = np.einsum("ij,ij->i", g, normals)
    terms = {"pressure": -along_meridian(normals, v_hat, dp_ds),
             "deceleration": rho_l * (g - g_n[:, None] * normals),
             "thickness": rho_l * g_n[:, None] * columns.gradient(H),
             "capillary": _capillary(columns, columns.film_curvature(H), wet, sigma)}
    return sum(terms.values()), terms


# -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
SLIVER = 1.0e-12                  # a piece of a donor cell shorter than this share of it is the round-off of a boundary
                                  # the donor's and the receiver's partitions share, and is left out
TRANSPORT_TOL = 1.0e-13           # the iterative solve stops once no row's residual exceeds this share of the largest
                                  # |b| + |A| |x| of its right-hand side's column (a normwise backward error; Task 6d)
TRANSPORT_CORRECTIONS = 40        # corrections of the iterative solve at most, before the direct LU takes over
TRANSPORT_CYCLE = 500             # cells: a cycle of meeting flows larger than this is iterated, not factored (Task 6d;
                                  # one of 1 700 cells filled its triangles 1.1 times, one of 10 000 5-6 times)
TRANSPORT = {"solves": 0, "corrections": 0, "fallbacks": 0}   # the transport's systems solved, the iteration's
                                  # corrections and the solves the direct LU finished, for measuring the cost
LANDED_PASSES = 200               # passes at most of the estimate of the landed metal's velocities over the columns
LANDED_TOL = 1.0e-6               # it stops once no estimated velocity moves by more than this share of the fastest


@dataclass
class TransportBook:
    """What a transport sub-step carried between the columns, per column: each column's gain from its neighbours less
    what it gave them. The mass, the enthalpy and the momentum as it left the donors each sum to zero over the body to
    round-off; every arrival's velocity is projected onto its receiver's plane, and the normal parts removed are the
    wall's impulse, so that each column's momentum changes by momentum + wall."""
    mass: np.ndarray         # kg (n,)
    enthalpy: np.ndarray     # J (n,)
    momentum: np.ndarray     # kg m/s (n, 3), as it left the donors, before the projection
    wall: np.ndarray         # N s (n, 3), the normal parts removed, each along its own column's normal
    courant: float           # the largest outflow rate of any cell [1/s] times the sub-step
    landed: np.ndarray = None   # (n,) bool: the columns metal landed on top of, a new top cell (sub-plan 21 Task 9)


def tangent_frames(normals):
    """Two orthonormal vectors (e1, e2), each (n, 3), spanning every column's tangent plane, e1 x e2 along the normal:
    e1 the global x axis's part in the plane, or y's where the normal lies within 26 degrees of x."""
    ref = np.where((np.abs(normals[:, 0]) < 0.9)[:, None], np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    e1 = ref - np.einsum("ij,ij->i", ref, normals)[:, None] * normals
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    return e1, np.cross(normals, e1)


def _fractions(depth):
    """Per column, its cells' boundaries (n, K+1) as fractions of the depth they hold (`depth` (n, K) each cell's share,
    zero for a cell that does not count), from 0 at the top to 1 at the last that counts and held there; 1 throughout
    where none counts."""
    z = np.concatenate([np.zeros((depth.shape[0], 1)), np.cumsum(depth, axis=1)], axis=1)
    total = z[:, -1:]
    return np.where(total > 0.0, z / np.where(total > 0.0, total, 1.0), 1.0)


def _pieces(sig, tau):
    """The pieces two partitions of [0, 1] cut each other into, row by row: sig and tau (D, K+1) non-decreasing
    boundaries. Returns, for every piece of positive length, its row, its cell of sig, its cell of tau and its length.
    A rank merge: each boundary's place in the merged row is its own index plus the number of the other partition's
    boundaries before it (tau's first where two are equal)."""
    D, B = sig.shape
    rows = np.arange(D)[:, None]
    off = (np.arange(D) * 4.0)[:, None]                      # boundaries lie in [0, 1]: rows 4 apart never mix
    before = lambda xp, x, side: (np.searchsorted((xp + off).ravel(), (x + off).ravel(), side=side).reshape(x.shape)
                                  - rows * B)
    at_s = np.arange(B)[None, :] + before(tau, sig, "right")
    at_t = np.arange(B)[None, :] + before(sig, tau, "left")
    merged = np.empty((D, 2 * B))
    merged[rows, at_s], merged[rows, at_t] = sig, tau
    from_sig = np.zeros((D, 2 * B), dtype=np.int64)
    from_sig[rows, at_s] = 1
    seen = np.cumsum(from_sig, axis=1)[:, :-1]                # sig's boundaries at or before each piece's start
    length = np.diff(merged, axis=1)
    r, p = np.nonzero(length > 0.0)
    return r, seen[r, p] - 1, p - seen[r, p], length[r, p]


def _upwind_order(src, dst, n):
    """An order of n unknowns in which every donor (src) comes before its receivers (dst), except within a cycle of
    transfers, which stays together; whether it holds; and the size of every unknown's cycle (its strong component; 1
    outside any cycle). The order is scipy's strong components, which its (Pearce's) algorithm numbers in reverse
    topological order -- how scipy numbers them, not what it promises, so it is checked."""
    graph = sp.csr_matrix((np.ones(len(src)), (src, dst)), shape=(n, n))
    _, label = connected_components(graph, directed=True, connection="strong")
    return np.argsort(-label, kind="stable"), bool(np.all(label[src] >= label[dst])), np.bincount(label)[label]


def _factor(A, upwind):
    """A sparse LU of A, an M-matrix or one by blocks, in the given order with no pivoting where the order is upwind (the
    matrix is then triangular but for its cycles: no fill, every sign kept), else in SuperLU's own fill-reducing order."""
    if upwind:
        return spla.splu(A, permc_spec="NATURAL", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
    return spla.splu(A)


def _solve(A, rhs, upwind, big=None):
    """x with A x = rhs, rhs (N,) or (N, p), for the transport's systems (sub-plan 21 Task 6d). In upwind order A is
    lower triangular but for the back-couplings of its cycles, where flows meet, and its LU fills only within them:
    a cycle of a few cells costs nothing, but one where meeting flows join thousands of cells filled the triangles 5-6
    times, and the LU then took most of a flow sub-step (6c). `big` (N,) marks the unknowns in such cycles, of more than
    TRANSPORT_CYCLE cells: the back-couplings among them, U, are left out of the factored matrix M = A - U, which fills
    only within the small cycles, and x_(k+1) = M^-1 (b - U x_k) from x_0 = M^-1 b, whose residual b - A x_(k+1) is
    U (x_k - x_(k+1)), until no row's residual exceeds TRANSPORT_TOL of the largest |b| + |A| |x| in its column, at most
    TRANSPORT_CORRECTIONS corrections; then the LU of A itself. Where no cycle is large, or `big` is None, M is A and x
    its solve, as 6a's; where the order is not upwind, the LU in SuperLU's own order (`_factor`).

    For the masses M is an M-matrix and N = -U >= 0 (a regular splitting), so M^-1 >= 0 and from x_0 >= 0 the iterates
    rise monotonically to the solution, every residual b - A x being non-negative: the masses installed from an
    iterate's fluxes, x + (b - A x), are positive (`transport`).

    Measured 2026-10-22 on 6c's states of the 100 mm sphere, one flow sub-step of 1 ms, the mass's system and the
    momentum's: with every column melted (one cycle of 18 379 cells) 6c's LU took 4.65 and 25.3 s, filling its
    triangles 6.0 and 6.3 times, and this solve 0.16 and 0.25 s in 9 corrections each; on a quiescent melt (10 190
    cells) 0.20 and 0.76 s against 0.10 and 0.28 s (3); on the heated melting state (11 109) 0.43 and 2.04 s against
    0.14 and 0.21 s (5); driven at 93 s (1 700 cells, filling 1.1 times) 0.03 and 0.07 s against 0.07 and 0.10 s (4).
    Each within 2e-13 of the LU's answer, the masses non-negative."""
    TRANSPORT["solves"] += 1
    if not upwind:
        return _factor(A, upwind).solve(rhs)
    U = sp.triu(A, k=1, format="coo")
    keep = np.zeros(U.nnz, dtype=bool) if big is None else big[U.row] & big[U.col]
    U = sp.csr_matrix((U.data[keep], (U.row[keep], U.col[keep])), shape=A.shape)
    if U.nnz == 0:
        return _factor(A, upwind).solve(rhs)
    lu = _factor((A - U).tocsc(), upwind)
    x = lu.solve(rhs)
    Ux = U @ x
    limit = TRANSPORT_TOL * (np.abs(rhs) + abs(A) @ np.abs(x)).max(axis=0)
    r = -Ux
    for k in range(TRANSPORT_CORRECTIONS + 1):
        if np.all(np.abs(r).max(axis=0) <= limit):
            TRANSPORT["corrections"] += k
            return x
        if k == TRANSPORT_CORRECTIONS:
            break
        x = lu.solve(rhs - Ux)
        Ux, r = U @ x, Ux
        r = r - Ux
    TRANSPORT["fallbacks"] += 1
    return _factor(A, upwind).solve(rhs)


def _landed_velocity(columns, area, run, dt, base_m, base_p, stay):
    """The velocity [m/s, (n, 3)] at which the metal landing on each column's rigid top flows on within the transport's
    sub-step (sub-plan 21 Task 6d): the mass-weighted mean of what lands there, each arrival's velocity projected onto
    the column's plane -- the landed metal coasts at the velocity it arrived with until the next momentum step drags it
    on the rigid top -- zero where nothing lands, where the top flows (`run` > 0: the metal joins it by fraction) or
    where a crust landed too (`stay`: rigid metal on a rigid top is anchored). What lands from the cells, base_m [kg,
    (n,)] and base_p [kg m/s, (n, 3)] with its velocities projected, is known before the solve; what the landed metal
    carries on onto further rigid tops is estimated column by column -- a column sends dt c m / (1 + dt sum c) of what
    landed, c = (u . f)+ / A at its estimated velocity u -- pass by pass, a pass reaching one column further, until no
    velocity moves by more than LANDED_TOL of the fastest, at most LANDED_PASSES passes; metal estimated at less than
    SLIVER of the largest landing stays where it lands (the far end of the implicit scheme's geometric tail). These are
    estimates for the rates alone: the solve carries the masses, the enthalpy and the momentum exactly."""
    n = len(area)
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    src, dst, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
    normal = columns.normal
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
        rate = np.maximum(np.einsum("dx,dx->d", u[src[d]], f[d]), 0.0) / area[src[d]]
        out = np.bincount(src[d], rate, n)
        sent = dt * rate * m_est[src[d]] / (1.0 + dt * out[src[d]])
        d, sent = d[onto[d]], sent[onto[d]]
        v = u[src[d]] - np.einsum("dx,dx->d", u[src[d]], normal[dst[d]])[:, None] * normal[dst[d]]
        m_est = base_m + np.bincount(dst[d], sent, n)
        p_est = base_p + np.stack([np.bincount(dst[d], sent * v[:, x], n) for x in range(3)], axis=1)
    return u


def transport(sk, columns, dt, regrid=True, faces=None):
    """One linearly implicit upwind sub-step dt [s] of the metal moving sideways between neighbouring columns at the
    cells' velocities sk.U (spec 6.5; Asha, 2026-10-21: its own step, after the momentum step). Updates the cells and
    returns a TransportBook. `faces`, a FacePush from `thickness_step` (sub-plan 21 Task 6b), adds to every moving cell's
    velocity across each face its response times the face's own thickness push, less its column's; None leaves the
    velocities as they are.

    Across each face, cell k of the donor column a sends c m_k' toward the receiver b, c = (u_k . f)+ / A_a [1/s], f the
    face's vector (`Columns.face_vector`) and m_k' the cell's mass at the sub-step's end; it carries the cell's specific
    enthalpy and velocity. Every cell is then m_k' (1 + dt sum c) = m_k + dt (what arrives): (I + dt C) m' = m, the
    film's runoff scheme (`film.Runoff`) on the cells -- an M-matrix whose columns sum to one, so the masses stay
    positive and nothing limits dt. The enthalpy, the carried sensitivity S and the momentum solve the same system, the
    momentum with every arrival projected onto its receiver's tangent plane, in two components per cell
    (`tangent_frames`), the normal part booked to the wall (spec 6.4, section 14 assumption 2). Every cell is installed
    from the solve's fluxes -- what it had, less what it sent, plus what arrived, each arrival taking dt c times its
    donor's solved value -- so the books close exactly whatever the solve's own error (Task 6d).

    Where metal joins the receiver (spec 6.5, section 14 assumption 3): a flowing donor cell's span of fractions of the
    donor's flowing depth (its flowing cells, top down) joins the receiver's cells between the same fractions of its
    flowing layer, the run of flowing cells from its top, each in proportion to the overlap. Where the receiver's top is
    rigid -- anchored metal, or a crust -- everything lands on top of it, as a new top cell; so does a floating crust's
    metal (Task 4, Asha 2026-10-20) wherever it goes, being no part of any flowing depth. Rigid metal joined to the
    substrate (`anchored_cells`) neither leaves nor receives, whatever velocity it holds: nothing passes through rigid
    metal. Only the cells that send or receive are solved for; with every velocity zero nothing changes.

    Flowing metal that lands on a rigid top flows on in the same solve (Task 6d): its new top cell is a flowing donor
    of the same system, sending across its column's faces at the velocity it arrived with (`_landed_velocity`: the
    arrivals' mean, projected onto the column's plane, fixed before the solve so that the step stays linear), by
    fraction into a flowing receiver's whole flowing layer and on top of a rigid one, whose new cell flows on in turn.
    6a's rule held it there until the next sub-step, so a front crossing rigid metal advanced at most one column a
    sub-step: 2 m/s on 2 mm columns at 1 ms, which set 6c's flow sub-step. Its momentum is the solve's, the arrivals'
    projected; the drag of the rigid top beneath reaches it at the next momentum step. A new top cell that crust metal
    reached stays where it lands (rigid metal on a rigid top is anchored).

    The cells are taken in upwind order (`_upwind_order`), in which the system is triangular but for its cycles where
    flows meet, and solved by `_solve`: the small cycles factored, those of more than TRANSPORT_CYCLE cells iterated.
    With `regrid` the cells are then split and merged (spec 6.2); a sub-step that regrids anyway passes False.

    Measured 2026-10-21 on the 100 mm sphere's skins (9 417 columns, 10 um cells), flowing along the meridian away from
    the stagnation point: 0.27 s a sub-step with the front hemisphere's top 0.3 mm flowing (124 000 of 553 000 cells
    moving), the same at Courant numbers of 1 and 96, and 0.95 s with all 584 000 cells moving, against 0.6-0.8 s for a
    coupled conduction sub-step on the sphere's skins at 700 K. With every cell moving, the mass's LU in the cells' own
    order took 19.6 s and in SuperLU's fill-reducing order 0.66 s, against 0.06 s upwind. Projected once at the
    sub-step's end instead of at every arrival, a
    uniform meridian flow kept 90 % of its kinetic energy in one sub-step at Courant number 9.6 against 96 %, and 24 %
    against 49 % at 96: metal crossing many columns in a sub-step would carry its velocity straight on, off the curved
    surface."""
    n, K = sk.m.shape
    valid = sk.valid
    flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
    moving = valid & ~anchored_cells(sk, flowing) & (sk.U != 0.0).any(axis=2)
    if faces is not None:
        moving |= valid & ~anchored_cells(sk, flowing) & (faces.response != 0.0)
    book = TransportBook(np.zeros(n), np.zeros(n), np.zeros((n, 3)), np.zeros((n, 3)), 0.0)
    if not moving.any():
        return book
    # every face a moving column has, each way, and every cell's outflow rate across it [1/s]
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    a, b, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
    keep = np.flatnonzero(moving.any(axis=1)[a])
    a, b, f = a[keep], b[keep], f[keep]
    across = np.einsum("dkx,dx->dk", sk.U[a], f)
    if faces is not None:
        across = across + faces.response[a] * faces.across[keep][:, None]
    c = np.where(moving[a], np.maximum(across, 0.0), 0.0) / sk.area[a][:, None]
    keep = np.flatnonzero((c > 0.0).any(axis=1))
    if not len(keep):
        return book
    a, b, c = a[keep], b[keep], c[keep]
    D = len(a)
    # where it joins: by fraction of the flowing depth, or on top
    dz = sk.dz()
    run = np.cumprod(flowing, axis=1).sum(axis=1)             # the receiver's flowing layer: its top run of flowing cells
    sig = _fractions(np.where(flowing, dz, 0.0))
    tau = _fractions(np.where(np.arange(K)[None, :] < run[:, None], dz, 0.0))
    joins = np.where(flowing[a] & (run[b] > 0)[:, None], c, 0.0)
    rows = np.flatnonzero((joins > 0.0).any(axis=1))
    r, kd, kr, length = _pieces(sig[a[rows]], tau[b[rows]])
    some = joins[rows[r], kd] > 0.0
    edge, kd, kr, length = rows[r[some]], kd[some], kr[some], length[some]
    group = edge * K + kd
    span = np.bincount(group, length, D * K)[group]
    some = length > SLIVER * span                             # not the round-off of a boundary the two share
    edge, kd, kr, length, group = edge[some], kd[some], kr[some], length[some], group[some]
    share = length / np.bincount(group, length, D * K)[group]  # each donor cell's pieces sum to one
    te, tk = np.nonzero((c > 0.0) & ~(flowing[a] & (run[b] > 0)[:, None]))
    # every arrival: its face, the receiving cell (n K + b: a new top cell of column b), the sending cell, dt c share
    e_edge = np.r_[edge, te]
    e_row = np.r_[b[edge] * K + kr, n * K + b[te]]
    e_col = np.r_[a[edge] * K + kd, a[te] * K + tk]
    e_val = dt * np.r_[c[edge, kd] * share, c[te, tk]]
    # metal landing on a rigid top flows on in the same solve (Task 6d): the new top cell sends across its column's
    # faces at the velocity of what lands there, by fraction into a flowing receiver and on top of a rigid one
    lands = flowing[a[te], tk] & (run[b[te]] == 0)
    stay = np.zeros(n, dtype=bool)
    stay[b[te[~flowing[a[te], tk]]]] = True                  # a crust landed: rigid metal on a rigid top stays
    out = np.bincount(e_col, e_val, n * K)[a[te] * K + tk]
    w = (e_val[len(edge):] * sk.m[a[te], tk] / (1.0 + out))[lands]   # what lands from each cell, estimated
    rec = b[te[lands]]
    v = sk.U[a[te], tk][lands]
    v = v - np.einsum("dx,dx->d", v, columns.normal[rec])[:, None] * columns.normal[rec]
    u_on = _landed_velocity(columns, sk.area, run, dt, np.bincount(rec, w, n),
                            np.stack([np.bincount(rec, w * v[:, x], n) for x in range(3)], axis=1), stay)
    s_all, r_all, f_all = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
    on = np.flatnonzero((u_on != 0.0).any(axis=1)[s_all])
    c_on = np.maximum(np.einsum("dx,dx->d", u_on[s_all[on]], f_all[on]), 0.0) / sk.area[s_all[on]]
    on, c_on = on[c_on > 0.0], c_on[c_on > 0.0]
    reach = np.zeros(n, dtype=bool)                           # the new top cells metal reaches: only those send
    reach[rec] = True
    while len(on):
        more = reach[s_all[on]] & (run[r_all[on]] == 0) & ~reach[r_all[on]]
        if not more.any():
            break
        reach[r_all[on][more]] = True
    on, c_on = on[reach[s_all[on]]], c_on[reach[s_all[on]]]
    if len(on):
        s, r = s_all[on], r_all[on]
        into = np.flatnonzero(run[r] > 0)                     # by fraction of the receiver's flowing layer, whole
        piece = np.diff(tau[r[into]], axis=1) * (np.arange(K)[None, :] < run[r[into]][:, None])
        q, kq = np.nonzero(piece > 0.0)
        onto = np.flatnonzero(run[r] == 0)                    # on top of a rigid receiver, its new top cell
        e_edge = np.r_[e_edge, D + into[q], D + onto]
        e_row = np.r_[e_row, r[into[q]] * K + kq, n * K + r[onto]]
        e_col = np.r_[e_col, n * K + s[into[q]], n * K + s[onto]]
        e_val = np.r_[e_val, dt * c_on[into[q]] * piece[q, kq], dt * c_on[onto]]
        a, b = np.r_[a, s], np.r_[b, r]
        D = len(a)
    # the cells that send or receive, in upwind order; each sends what its arrivals take, so the books close exactly
    active = np.zeros(n * K + n, dtype=bool)
    active[e_row], active[e_col] = True, True
    act = np.flatnonzero(active)
    N = len(act)
    local = np.full(n * K + n, -1, dtype=np.int64)
    local[act] = np.arange(N)
    order, upwind, size = _upwind_order(local[e_col], local[e_row], N)
    big = size[order] > TRANSPORT_CYCLE                       # the unknowns in large cycles, in the solve's order
    at = np.empty(N, dtype=np.int64)
    at[order] = np.arange(N)
    to, frm = at[local[e_row]], at[local[e_col]]
    lr, lc = local[e_row], local[e_col]
    diag = 1.0 + np.bincount(lc, e_val, N)
    # mass, enthalpy and the carried sensitivity, installed from the solve's fluxes: each cell keeps what it had, less
    # what it sent, plus what arrived, so the books close exactly whatever the solve's own error (Task 6d)
    cell = lambda v: np.r_[v.reshape((n * K,) + v.shape[2:]), np.zeros((n,) + v.shape[2:])][act]
    b0 = np.column_stack([cell(sk.m), cell(sk.H)] + ([cell(sk.m * sk.S)] if sk.S is not None else []))
    x = np.empty_like(b0)
    x[order] = _solve(sp.csc_matrix((np.r_[diag, -e_val], (np.r_[at, to], np.r_[at, frm])), shape=(N, N)), b0[order],
                      upwind, big)
    moved = e_val[:, None] * x[lc]
    x = b0 + np.stack([np.bincount(lr, moved[:, v], N) - np.bincount(lc, moved[:, v], N) for v in range(x.shape[1])],
                      axis=1)
    # the momentum, each arrival projected onto its receiver's plane: two components per cell in its own plane's frame
    owner = np.where(act < n * K, act // K, act - n * K)
    e1, e2 = tangent_frames(columns.normal)
    P = cell(sk.m[..., None] * sk.U)
    p0 = np.column_stack([np.einsum("ij,ij->i", P, e1[owner]), np.einsum("ij,ij->i", P, e2[owner])])
    held = P - p0[:, :1] * e1[owner] - p0[:, 1:] * e2[owner]  # any normal part a cell held: the wall's
    R = np.stack([np.stack([np.einsum("ij,ij->i", e[b], g[a]) for g in (e1, e2)], axis=1) for e in (e1, e2)], axis=1)
    Re = R[e_edge]                                            # [arrival, receiver's component, donor's component]
    rows2 = [2 * at, 2 * at + 1] + [2 * to + q for q in (0, 0, 1, 1)]
    cols2 = [2 * at, 2 * at + 1] + [2 * frm + q for q in (0, 1, 0, 1)]
    vals2 = [diag, diag] + [-e_val * Re[:, q, s] for q, s in ((0, 0), (0, 1), (1, 0), (1, 1))]
    rhs2 = np.empty(2 * N)
    rhs2[2 * at], rhs2[2 * at + 1] = p0[:, 0], p0[:, 1]
    y = _solve(sp.csc_matrix((np.concatenate(vals2), (np.concatenate(rows2), np.concatenate(cols2))),
                             shape=(2 * N, 2 * N)), rhs2, upwind, np.repeat(big, 2))
    sent_p = e_val[:, None] * np.column_stack([y[2 * at], y[2 * at + 1]])[lc]
    got = np.einsum("eqs,es->eq", Re, sent_p)
    p = p0 + np.stack([np.bincount(lr, got[:, v], N) - np.bincount(lc, sent_p[:, v], N) for v in range(2)], axis=1)
    P = p[:, :1] * e1[owner] + p[:, 1:] * e2[owner]
    # the books, face by face: what every arrival took from its sending cell at the sub-step's end
    sent = np.stack([np.bincount(e_edge, v, D) for v in (moved[:, 0], moved[:, 1], sent_p[:, 0], sent_p[:, 1])], axis=1)
    carried = sent[:, 2:3] * e1[a] + sent[:, 3:4] * e2[a]
    lost = np.einsum("ij,ij->i", carried, columns.normal[b])[:, None] * columns.normal[b]
    net = lambda v: np.bincount(b, v, n) - np.bincount(a, v, n)
    book = TransportBook(net(sent[:, 0]), net(sent[:, 1]), np.stack([net(carried[:, v]) for v in range(3)], axis=1),
                         -np.stack([np.bincount(b, lost[:, v], n) + np.bincount(owner, held[:, v], n)
                                    for v in range(3)], axis=1), float((diag - 1.0).max()))
    # install: new values for the cells solved for, a new top cell where metal landed on top
    m, H, T, U = sk.m.copy(), sk.H.copy(), sk.T.copy(), sk.U.copy()
    S = sk.S.copy() if sk.S is not None else None
    old = act < n * K
    m.ravel()[act[old]], H.ravel()[act[old]], T.ravel()[act[old]] = x[old, 0], x[old, 1], np.nan
    U.reshape(-1, 3)[act[old]] = P[old] / x[old, :1]
    if S is not None:
        S.ravel()[act[old]] = x[old, 2] / x[old, 0]
    count = sk.count.copy()
    tops = act[~old] - n * K
    if len(tops):
        new = np.zeros(n, dtype=bool)
        new[tops] = True
        src = np.clip(np.arange(K + 1)[None, :] - new[:, None], 0, K)
        down = lambda v: np.take_along_axis(np.concatenate([v, np.zeros((n, 1) + v.shape[2:])], axis=1),
                                            src.reshape(src.shape + (1,) * (v.ndim - 2)), axis=1)
        m, H, T, U = down(m), down(H), down(T), down(U)
        m[tops, 0], H[tops, 0], T[tops, 0], U[tops, 0] = x[~old, 0], x[~old, 1], np.nan, P[~old] / x[~old, :1]
        if S is not None:
            S = down(S)
            S[tops, 0] = x[~old, 2] / x[~old, 0]
        count = count + new
        book.landed = new
    sk._install(m, H, count, S, T, U)
    if regrid:
        sk.regrid()
    return book


# -- the thickness terms, implicit (sub-plan 21 Task 6b) ----------------------------------------------------------------
THICKNESS_PASSES = 8              # linearisations of the transport's fluxes at most
THICKNESS_TOL = 1.0e-3            # the linearisations stop once the flux the last one missed would move no column by
                                  # more than this share of the largest height change it predicts


def flowing_tops(sk):
    """Whether each column's top cell flows, (n,) bool: where surface tension acts (Asha, 2026-10-21), not on an
    anchored rigid top nor on a floating crust."""
    return (sk.count > 0) & np.isfinite(cell_viscosity(sk, 0.0)[:, 0])


def _walled(columns, kappa, wet):
    """The curvatures `kappa` (n,) where a column's top flows; at a column whose top does not -- a wall -- the mean of
    its flowing neighbours', so that the gradient beside a wall carries no capillary push of the wall's own."""
    if wet is None or wet.all():
        return kappa
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    total = np.bincount(i, kappa[j] * wet[j], columns.n) + np.bincount(j, kappa[i] * wet[i], columns.n)
    count = np.bincount(i, wet[j], columns.n) + np.bincount(j, wet[i], columns.n)
    return np.where(wet, kappa, np.where(count > 0.0, total / np.maximum(count, 1.0), kappa))


def _capillary(columns, kappa, wet, sigma):
    """-sigma grad_s kappa where a column's top flows, zero elsewhere (`push`'s capillary term)."""
    G = -sigma * columns.gradient(_walled(columns, kappa, wet))
    return G if wet is None else np.where(wet[:, None], G, 0.0)


@dataclass
class FacePush:
    """What `transport` adds to the cells' velocities across the faces."""
    response: np.ndarray     # (n, K) [m/s per N/m3]: every cell's velocity per unit push (MomentumBook.response)
    across: np.ndarray       # (2E,) [Pa]: for each face, from its edge's first column to its second and then the other
                             # way (the order `transport` builds), the face's own thickness push times the face, less
                             # the sending column's push dotted with the face's vector


@dataclass
class ThicknessBook:
    """What the implicit thickness terms did in a sub-step, per column."""
    push: np.ndarray         # N s (n, 3): the impulse of their change on the moving cells
    wall: np.ndarray         # N s (n, 3): the anchored metal's drag on the velocities that change gave
    heights: np.ndarray      # m (n,): the thicknesses predicted at the end of the transport
    faces: FacePush          # for `transport`
    growth: float            # the largest Rayleigh-Taylor growth in a sub-step, dt s, over the columns
    passes: int              # the linearisations made


def thickness_step(sk, columns, dt, book, applied, g, rho_l, sigma):
    """The push's thickness terms, rho_l (g . n) grad_s H and -sigma grad_s kappa, implicit in the column heights at
    the end of the sub-step's transport (spec 6.4's last paragraph; Asha, 2026-10-21: the transport is its own step
    after the momentum step). Called between `momentum_step`, given `push`'s G with `wet=flowing_tops(sk)` and
    response=True (its book is `book`), and `transport`, given the returned book's `faces`; `applied` (n, 3) [N/m3] is the
    part of that push the two terms made, terms['thickness'] + terms['capillary']. Corrects sk.U; returns a
    ThicknessBook.

    The two terms are the gradient of P = sigma kappa - rho_l (g . n) H, kappa the thin film's curvature
    (`Columns.film_curvature`), surface tension only where a column's top flows (`flowing_tops`). Across each dual face
    P pushes by its two-point difference, phi = w (P_a - P_b) [Pa] with w the face's clipped cotangent weight
    (`Columns.conduct_weight`; the capillary part only between two flowing tops, the deceleration's with the mean of the
    two columns' g . n), linear in the heights through `Columns.film_operator`. Each sub-step:

      1. predicts the heights at the end of the transport, H + dH, from the transport linearised on the columns: a column
         sends across a face what its cells carry at their velocity without the thickness push plus their response
         times the face's push -- implicitly in its moving metal, as the transport does -- linearised in the face's push
         through the mobility of the cells that cross, sum r dz. Newton's method on that piecewise linear flux: the
         system is solved again at the new face pushes until the flux the last solve missed would move no column by more
         than THICKNESS_TOL of the largest change (1-3 solves on the test sphere and 3 on the 100 mm sphere's melting
         state, measured 2026-10-21);
      2. corrects every moving cell's velocity by its response times the change of the column's push -- `push`'s terms at
         the predicted heights, less `applied` -- so that the cells' momentum follows the push of Task 5's form, and books
         the impulse of that change and the anchored metal's drag on it (`MomentumBook.wall_response`);
      3. hands `transport` each face's own push (ruling, 2026-10-21; the column form of Rhie-Chow interpolation): a
         cell's velocity across a face becomes its velocity without the thickness push plus its response times that
         face's push. The two-point difference sees a disturbance one column wide; the columns' pushes are centred
         differences, blind to it, and with them alone a bump 1 um high on one column of 1.1 mm of slurry stayed at
         99.8 % of its height while the faces' push took it to 12.5 % in one 10 ms sub-step.

    Measured 2026-10-21. A cosine 24 mm long on 1 mm of slurry at 902 K pulled out of its surface at 100 m/s2 grows at
    the thin-film rate h^3 k^2 (rho g - sigma k^2) / 3 mu within 0.42 % at sub-steps of 1-10 ms; one of 8 mm, below the
    11.5 mm cutoff, decays at it within 1.2 % at 1 ms and within 0.26 % of backward Euler's own rate at 1-10 ms; a mode
    alternating from column to column on 0.3 mm of liquid on 2 mm columns falls to at most 0.12 of itself in one sub-step
    of 10 ms to 1 s. On the 100 mm sphere's melting state (the face where g . n > 0 melted 0.3 mm deep, 141 000 cells
    moving) the step takes 0.30 s and the response 0.03-0.05 s more in the momentum step, against 0.31-0.37 s for the
    transport. The growth reported is the largest dt M k^2 (rho_l g . n - sigma k^2) over the columns, M = sum r dz their
    mobility and k at its fastest up to the columns' pi / sqrt(A): backward Euler follows the Rayleigh-Taylor growth only
    where it is well under one."""
    n, K = sk.m.shape
    r = book.response
    dz, A = sk.dz(), sk.area
    H = dz.sum(axis=1)
    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
    moving = sk.valid & ~anchored_cells(sk, flowing)
    wet = flowing_tops(sk)
    normals = columns.normal
    g = np.broadcast_to(np.asarray(g, dtype=float), (n, 3))
    g_n = np.einsum("ij,ij->i", g, normals)
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    w, g_e, both = columns.conduct_weight, 0.5 * (g_n[i] + g_n[j]), (wet[i] & wet[j]).astype(float)
    kappa, J = columns.film_curvature(H), columns.film_operator(H)
    phi0 = w * (sigma * both * (kappa[i] - kappa[j]) - rho_l * g_e * (H[i] - H[j]))      # each face's, first to second
    Bt = columns.incidence.T.tocsr()
    dphi = sp.diags(w) @ (sp.diags(sigma * both) @ (Bt @ J) - sp.diags(rho_l * g_e) @ Bt)  # its change per change of H
    # every face each way whose sender moves; its cells' velocity across it without the thickness push [m2/s]
    a, b, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
    keep = np.flatnonzero(moving.any(axis=1)[a])
    a, b, f, phi_d = a[keep], b[keep], f[keep], np.r_[phi0, -phi0][keep]
    dphi_d = sp.vstack([dphi, -dphi]).tocsr()[keep]
    un = np.where(moving[a], np.einsum("dkx,dx->dk", (sk.U - r[..., None] * applied[:, None, :])[a], f), 0.0)
    ra, da = r[a], np.where(moving[a], dz[a], 0.0)
    flux = lambda phi: (np.maximum(un + ra * phi[:, None], 0.0) * da).sum(axis=1)       # m3/s
    thick = np.where(moving, dz, 0.0).sum(axis=1)
    D = len(a)
    C = sp.csr_matrix((np.r_[np.ones(D), -np.ones(D)], (np.r_[a, b], np.r_[np.arange(D), np.arange(D)])), shape=(n, D))
    active = np.zeros(n, dtype=bool)
    active[a], active[b] = True, True
    act = np.flatnonzero(active)
    dH, phi, passes = np.zeros(n), phi_d, 0
    F = flux(phi)
    for passes in range(1, (THICKNESS_PASSES if D else 0) + 1):
        M = (ra * da * (un + ra * phi[:, None] > 0.0)).sum(axis=1)                       # m5/(N s): the crossing cells'
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
    # the cells' velocities by the change of the columns' push; the faces' own push for the transport
    G = rho_l * g_n[:, None] * columns.gradient(H + dH) + _capillary(columns, kappa + J @ dH, wet, sigma)
    dG = G - applied
    sk.U = np.where(moving[..., None], sk.U + r[..., None] * dG[:, None, :], sk.U)
    across = np.zeros(2 * len(i))
    across[keep] = phi - np.einsum("dx,dx->d", G[a], f)
    # the books, and the fastest Rayleigh-Taylor growth in a sub-step
    V = np.where(moving, dz, 0.0).sum(axis=1) * A
    s_tension = sigma * wet
    k2 = np.pi ** 2 / A
    fastest = np.maximum(rho_l * g_n, 0.0) / (2.0 * np.maximum(s_tension, 1e-300))         # k^2 of the fastest growth
    k2 = np.where(s_tension > 0.0, np.minimum(k2, fastest), k2)
    s = (r * np.where(moving, dz, 0.0)).sum(axis=1) * k2 * (rho_l * g_n - s_tension * k2)
    return ThicknessBook(push=dt * V[:, None] * dG, wall=-dt * book.wall_response[:, None] * dG, heights=H + dH,
                         faces=FacePush(r, across), growth=float(max(dt * s.max(), 0.0)) if n else 0.0, passes=passes)


# -- the flow in the sub-steps (sub-plan 21 Task 6c) --------------------------------------------------------------------
FLOW_CYCLES = 10                  # flow sub-steps in a heat sub-step where the skins remember no flow (6c's fixed count,
                                  # 1 ms under the 10 ms heat sub-step, chosen on the 100 mm sphere's runoff when a front
                                  # crossing rigid metal advanced one column a flow sub-step; sub-plan 21 Task 6c)
FLOW_CYCLES_MIN, FLOW_CYCLES_MAX = 1, 40   # the adaptive count's floor and cap (Task 6d): 10 ms and 0.25 ms flow
                                  # sub-steps under the 10 ms heat sub-step
FLOW_COURANT = 6.0                # the adaptive count's target: the largest Courant number of a flow sub-step (Task 6d).
                                  # On the 100 mm sphere's runoff at 93 s, with the front limit lifted, 50 ms of flow at
                                  # 1.25 ms (Courant 6.1) put the metal within 2.5 % of the 0.25 ms run's places and
                                  # its kinetic energy within 5 %, the thickness step's linearisations at 4 at most; at
                                  # 1.5 ms (7.3) 3.6 % and 12 %, the linearisations at their cap of 8; at 1 ms 2.2 %
FLOW_GROWTH = 0.02                # and the largest Rayleigh-Taylor growth dt s of one: backward Euler's growth rate is
                                  # then within about dt s / 2 = 1 % (6b's cosine met it within 0.42 % at 0.016)
FLOW_RELAX = 2.0                  # a flow sub-step at most this many times the last heat sub-step's
FLOW_SECONDS = {"push": 0.0, "momentum": 0.0, "thickness": 0.0, "transport": 0.0, "regrid": 0.0, "release": 0.0,
                "landing": 0.0}
                                  # wall time, for measuring the cost (as KRYLOV counts the coupled solves)


def flow_cycles(memory, dt):
    """The number of flow sub-steps in a heat sub-step dt [s] (sub-plan 21 Task 6d), from the skins' `flow_memory` --
    the last heat sub-step's largest Courant number C and Rayleigh-Taylor growth g over its flow sub-steps, and their
    length h_0 -- by the scalings measured on the 100 mm sphere: the Courant number in proportion to the flow sub-step,
    the growth to its square. The flow sub-step is the longest of FLOW_RELAX h_0, h_0 FLOW_COURANT / C and
    h_0 (FLOW_GROWTH / g)^(1/2), so that the next one meets both targets and grows by at most FLOW_RELAX a heat
    sub-step, and the count dt over it, rounded up, within FLOW_CYCLES_MIN and FLOW_CYCLES_MAX; FLOW_CYCLES where the
    skins remember no flow (at the start, or when nothing flowed in the last heat sub-step). Deterministic: the same
    memory and dt give the same count."""
    courant, growth, length = (float(x) for x in memory)
    if length <= 0.0:
        return FLOW_CYCLES
    h = FLOW_RELAX * length
    if courant > 0.0:
        h = min(h, length * FLOW_COURANT / courant)
    if growth > 0.0:
        h = min(h, length * math.sqrt(FLOW_GROWTH / growth))
    return int(min(max(math.ceil(dt / h - 1e-9), FLOW_CYCLES_MIN), FLOW_CYCLES_MAX))


LANDING_DEPTH = 0.4e-3             # m below a landed column's top that the landing conduction takes with it (Task 9)


def landing_conduct(sk, landed, dt, depth=None):
    """Metal that landed on a rigid top loses its heat to the metal beneath it within the flow sub-steps (sub-plan 21
    Task 9; spec 6.5: runoff landing on cold metal freezes from below in the skin's own heat equation). The cells of
    each column in `landed` (n,) whose tops lie within `depth` (LANDING_DEPTH) of its top -- the landed metal and the
    receiver's top cells, at least two cells -- are conducted through the depth for dt [s], implicitly (`_conduct` on
    those cells alone), with no heat across the block's top, whose heating and radiation are the heat sub-step's
    conduction's, nor across its bottom: the block's energy is its own, what the Newton iteration leaves (about 1e-10 J
    a cell) kept by its bottom cell, so the books are exact and the base law the trial pass makes is untouched; the
    trial pass's sensitivity is carried with the heat. The flow sub-steps conducted no heat before Task 9, so that
    landed metal stayed liquid and flowed on until the next heat sub-step's conduction: on Task 8's pushed plate a film
    of 10-20 um at 0.6-1.2 m/s ran 6-12 mm over metal at about 700 K in a 10 ms heat sub-step, where it freezes within
    some 0.005-0.03 ms (a contact estimate). With it, at 10 ms heat sub-steps 3.0 % of that plate's deposit lands 5-10 mm
    past the heat against 9.9 % without (0.9 and 0.7 % at 5 and 2.5 ms, where the two agreed), at 2-3 % of the flow's
    cost (measured 2026-10-26). The block's own cells were conducted with their column in the heat sub-step already, so
    their diffusion within the block is counted again: LANDING_DEPTH, the fine zone's depth, bounds that, about the
    depth heat reaches in a flow sub-step (sqrt(alpha dt), 0.33 mm at 2 ms on 700 K metal); at 0.2 and 0.8 mm 3.3 and
    2.7 % of the deposit land 5-10 mm past the heat. Updates the cells."""
    rows = np.flatnonzero(landed)
    if not rows.size:
        return
    depth = LANDING_DEPTH if depth is None else depth
    sub = sk.rows(rows)
    dz = sub.dz()
    k = np.maximum((sub.valid & (np.cumsum(dz, axis=1) - dz < depth)).sum(axis=1), np.minimum(2, sub.count))
    width = int(k.max())
    valid = np.arange(width)[None, :] < k[:, None]
    block = object.__new__(Skins)                                       # the blocks alone: the cells beneath hold nothing
    block.material, block.settings, block.rho_l, block.props = sk.material, sk.settings, sk.rho_l, sk.props
    block.area, block.count = sub.area, k
    block.m, block.H = np.where(valid, sub.m[:, :width], 0.0), np.where(valid, sub.H[:, :width], 0.0)
    block.T = np.where(valid, sub.T[:, :width], 300.0)
    block.S = None if sub.S is None else np.where(valid, sub.S[:, :width], 0.0)
    block.U, block.momentum_out, block.flow_memory = sub.U[:, :width], np.zeros((len(rows), 3)), sub.flow_memory
    zero = np.zeros(len(rows))
    E_base = _conduct(block, dt, zero, block.area, 0.0, zero, None)[2]  # an adiabatic block: only the iteration's error
    last = np.maximum(k - 1, 0)
    r = np.arange(len(rows))
    block.H[r, last] += E_base                                         # kept by the block, so its energy is its own
    block.T[r, last] = sk.material.temperature_from_enthalpy(block.H[r, last] / block.m[r, last])
    cells = np.nonzero(valid)
    sk.H[rows[cells[0]], cells[1]] = block.H[cells]
    sk.T[rows[cells[0]], cells[1]] = block.T[cells]
    if sk.S is not None:
        sk.S[rows[cells[0]], cells[1]] = block.S[cells]


@dataclass
class Flow:
    """What drives the skins' flow over a part of a macro step, held over it (spec 6.4): `dp_ds` [Pa/m, (n,)] the
    surface flow's wall-pressure gradient along the meridian, `shear` [Pa, (n, 3)] the gas shear on the top cells in
    each column's plane (`along_meridian` of the surface flow's tau), `g` [m/s2, (3,)] the acceleration the metal feels
    in the body's frame, minus the body's own, `v_hat` the flight direction, `rho_l` [kg/m3] and `sigma` [N/m] the
    liquid's density and surface tension, `A_top` [m2, (n,)] the top areas the shear acts over, and `cycles`, the flow
    sub-steps in each heat sub-step: None (the default since Task 6d) chooses them each heat sub-step (`flow_cycles`),
    a number fixes them. `plan` is what a trial pass chose, for its real pass (`advance`).

    `release` is the spray's placeholder (sub-plan 21 Task 7), what sub-plan 22's spraying supplies: None releases
    nothing; otherwise a callable `release(sk) -> (n,)` [kg/(m2 s)], the release rate per unit of each column's top area
    A_top at the skins' present state, evaluated in every flow sub-step after the transport and taken from the top cells
    by `release_top` -- liquid first, the exposed slurry too where `slurry` (spec 8's `--slurry-spray on`, a setting of
    sub-plan 22's). It is called in the trial pass, in every repeat of a part and in parts run again, so it keeps no
    records of its own: what was released is in the pass's FlowBook (`released`, `released_H`, `released_p`,
    `released_solid`). Since sub-plan 22 Task 1c it may give a `Cascade` instead, the liquid's rate and the exposed
    slurry's own per column, which the release takes cell by cell in time order (`cascade_mass`; Girin's,
    `skinbody.SprayPlan.release`), and the book then carries what the slurry's rate took (`released_slurry`). Every
    release's book carries it by the cells it came from too, liquid and slurry (`phase_mass`, ...), and a Cascade with a
    `histogram` has each flow sub-step's droplets binned by size (`hist_n`, `hist_m`; sub-plan 22 Task 2)."""
    columns: object
    dp_ds: np.ndarray
    shear: np.ndarray
    g: np.ndarray
    v_hat: np.ndarray
    rho_l: float
    sigma: float
    A_top: np.ndarray
    cycles: int = None
    plan: list = None
    release: object = None
    slurry: bool = False


@dataclass
class FlowBook:
    """What the flow sub-steps did, per column: impulses [N s] and momentum [kg m/s], (n, 3), which add up to the
    columns' momentum change less what the release took (`released_p`), the sub-steps' worst numbers, and the metal
    released from the top cells (Task 7; None in a single sub-step's book that released nothing)."""
    push: np.ndarray         # the push's impulse on the moving cells, its thickness terms' implicit change included
    shear: np.ndarray        # the gas shear's on the top cells
    wall: np.ndarray         # the anchored metal's drag (the momentum and thickness steps'), the momentum cells held on
                             # joining it, and the normal parts the transport's projections removed
    carried: np.ndarray      # what the transport carried in less what it sent, as it left the donors: zero over the body
    courant: float = 0.0     # the transport's largest outflow rate of any cell times its sub-step
    growth: float = 0.0      # the thickness step's largest Rayleigh-Taylor growth dt s
    passes: int = 0          # the thickness step's most linearisations
    sub_steps: int = 0       # the flow sub-steps taken (none where nothing flows)
    cycles: int = 0          # the most flow sub-steps in one heat sub-step (Task 6d)
    released: np.ndarray = None        # kg (n,): the metal released from the top cells (Task 7)
    released_H: np.ndarray = None      # J (n,): its enthalpy, at its cells' own specific enthalpies
    released_p: np.ndarray = None      # kg m/s (n, 3): its momentum, at its cells' own velocities -- fact 29's sink
    released_solid: np.ndarray = None  # kg (n,): the solid in it, by its cells' liquid fractions
    released_slurry: np.ndarray = None  # kg (n,): what the cascade took at the exposed slurry's own rate (sub-plan 22
                                        # Task 1c; zero for a release given as one rate)
    # the release by the cells it came from (sub-plan 22 Task 2): [0] fully liquid, at or above the liquidus, [1] the
    # others (the exposed slurry, only with slurry spraying) -- each spray mode's droplets
    phase_mass: np.ndarray = None      # kg (2, n)
    phase_H: np.ndarray = None         # J (2, n): its enthalpy, at its cells' own specific enthalpies
    phase_p: np.ndarray = None         # kg m/s (2, n, 3): its momentum, at its cells' own velocities
    phase_solid: np.ndarray = None     # kg (2, n): the solid in it
    hist_n: np.ndarray = None          # (bins,): the droplets released, binned by size each flow sub-step at the radii of
    hist_m: np.ndarray = None          # kg (bins,): the Cascade's `histogram`, and their mass (None without one)

    @classmethod
    def zero(cls, n):
        return cls(*(np.zeros((n, 3)) for _ in range(4)), released=np.zeros(n), released_H=np.zeros(n),
                   released_p=np.zeros((n, 3)), released_solid=np.zeros(n), released_slurry=np.zeros(n),
                   phase_mass=np.zeros((2, n)), phase_H=np.zeros((2, n)), phase_p=np.zeros((2, n, 3)),
                   phase_solid=np.zeros((2, n)))

    def add(self, other):
        """Add another book's impulses and releases and keep the worst numbers of both."""
        for key in ("push", "shear", "wall", "carried"):
            setattr(self, key, getattr(self, key) + getattr(other, key))
        for key in ("released", "released_H", "released_p", "released_solid", "released_slurry", "phase_mass", "phase_H",
                    "phase_p", "phase_solid", "hist_n", "hist_m"):
            mine, theirs = getattr(self, key), getattr(other, key)
            if theirs is not None:
                setattr(self, key, theirs.copy() if mine is None else mine + theirs)
        self.courant, self.growth = max(self.courant, other.courant), max(self.growth, other.growth)
        self.passes, self.sub_steps = max(self.passes, other.passes), self.sub_steps + other.sub_steps
        self.cycles = max(self.cycles, other.cycles)


def flow_substeps(sk, flow, dt, cycles=None):
    """`cycles` flow sub-steps of dt / cycles [s] each -- by default `flow.cycles`, or, where that is None, as many as
    `flow_cycles` chooses from the skins' `flow_memory`, which they then update (sub-plan 21 Task 6d) -- every one in the
    order of 6b's ruling: the push, surface tension only where a column's top flows (`push`, `flowing_tops`); the
    momentum through the depth with its response (`momentum_step`); the thickness terms implicit in the heights at the
    end of the transport (`thickness_step`); the transport with the faces' own push (`transport`, regrid=False); then,
    with `flow.release`, the spray's placeholder (sub-plan 21 Task 7): the release the rate gives at that moment over
    the flow sub-step, rate A_top h, from the top cells (`release_top`) -- spec 7.1's heat, flow, sideways exchange and
    spray, the last three in every flow sub-step. The cells are split and merged between the flow sub-steps, and the
    heat sub-step's own regrid follows the last (spec 6.2: once after each sub-step's exchanges). Updates the cells;
    returns a FlowBook. Where nothing flows nothing moves, nothing can be released and no sub-step is taken: any
    velocity left in the rigid cells is booked to the wall and cleared, as `momentum_step` books an anchored cell's, and
    the skins forget the flow.

    The flow's sub-step was fixed in 6c at 1 ms (FLOW_CYCLES) and is chosen since 6d from the last heat sub-step's
    Courant number and growth, against targets set from the sub-step study on the 100 mm sphere (FLOW_COURANT,
    FLOW_GROWTH). The release runs in every flow sub-step, not once a heat sub-step (Task 7's ruling): the top cells
    cross a column for each unit of their Courant number, which the count keeps near 6 a flow sub-step, so they cross
    some 25 columns of 2 mm in a heat sub-step at the 100 mm sphere's 5 m/s and up to 240 at the count's cap; a release
    once a heat sub-step would spray metal from where the runoff has taken it rather than from where its rate was
    found, and fact 29's split between the runoff and the stripping would be first order in the heat sub-step rather
    than in the flow's. Its cost: 0.05 s of processor time a release on the 100 mm
    sphere's 9 417 columns, against about 2 s for a 1 ms flow sub-step of its driven melting state (measured
    2026-10-24)."""
    n = sk.n
    book = FlowBook.zero(n)
    adaptive = flow.cycles is None
    if not (sk.valid & np.isfinite(cell_viscosity(sk, 0.0))).any():
        book.wall -= sk.column_momentum()
        sk.U = np.zeros_like(sk.U)
        if adaptive:
            sk.flow_memory = np.zeros(3)
        return book
    if cycles is None:
        cycles = flow_cycles(sk.flow_memory, dt) if adaptive else flow.cycles
    h = dt / cycles
    col = flow.columns
    clock = lambda key, started: FLOW_SECONDS.__setitem__(key, FLOW_SECONDS[key] + time.perf_counter() - started)
    landed = np.zeros(n, dtype=bool)                     # the columns metal landed on top of in this heat sub-step
    for k in range(cycles):
        if k:
            started = time.perf_counter()
            sk.regrid()
            clock("regrid", started)
        started = time.perf_counter()
        G, terms = push(col, sk.thickness(), flow.dp_ds, flow.g, flow.v_hat, flow.rho_l, flow.sigma, wet=flowing_tops(sk))
        clock("push", started)
        started = time.perf_counter()
        mb = momentum_step(sk, h, G, flow.shear, flow.A_top, response=True)
        clock("momentum", started)
        started = time.perf_counter()
        tb = thickness_step(sk, col, h, mb, terms["thickness"] + terms["capillary"], flow.g, flow.rho_l, flow.sigma)
        clock("thickness", started)
        started = time.perf_counter()
        tr = transport(sk, col, h, regrid=False, faces=tb.faces)
        clock("transport", started)
        if tr.landed is not None:
            landed |= tr.landed
        if landed.any():                                 # what landed this heat sub-step loses its heat (Task 9)
            started = time.perf_counter()
            landing_conduct(sk, landed, h)
            clock("landing", started)
        step = FlowBook(mb.push + tb.push, mb.shear, mb.wall + tb.wall + tr.wall, tr.momentum, tr.courant, tb.growth,
                        tb.passes, 1)
        if flow.release is not None:
            started = time.perf_counter()
            given = flow.release(sk)
            if isinstance(given, Cascade):               # the liquid's rate and the exposed slurry's (sub-plan 22 Task 1c)
                rate = np.broadcast_to(np.array([given.rate, given.rate_slurry], dtype=float), (2, n))
            else:
                rate = np.broadcast_to(np.asarray(given, dtype=float), (n,))
            if not np.all(np.isfinite(rate) & (rate >= 0.0)):
                raise ValueError("a release rate must be finite and non-negative")
            if isinstance(given, Cascade):
                dm, step.released_slurry = cascade_mass(sk, rate[0], rate[1], flow.A_top, h, flow.slurry, given.T_slurry)
            else:
                dm, step.released_slurry = rate * flow.A_top * h, np.zeros(n)
            rb = release_top(sk, dm, flow.slurry)
            step.released, step.released_H = rb.mass, rb.enthalpy
            step.released_p, step.released_solid = rb.momentum, rb.solid
            ph = rb.by_phase                             # by the cells it came from, and binned by size (Task 2)
            step.phase_mass, step.phase_H, step.phase_p, step.phase_solid = ph.mass, ph.enthalpy, ph.momentum, ph.solid
            if isinstance(given, Cascade) and given.histogram is not None:
                step.hist_n, step.hist_m = given.histogram(ph.mass)
            clock("release", started)
        book.add(step)
    book.cycles = cycles
    if adaptive:
        sk.flow_memory = np.array([book.courant, book.growth, h])
    return book


# -- the droplet momentum sink at the top cells (sub-plan 21 Task 7) ----------------------------------------------------
PI_FLAG = 0.2                     # fact 29: above this the droplets' momentum is no longer small against the gas shear
                                  # that drives the film, and the sink must be solved with the profile (as here it is)


@dataclass
class ReleaseBook:
    """What a release took from the top cells, per column."""
    mass: np.ndarray         # kg (n,)
    enthalpy: np.ndarray     # J (n,), at the cells' own specific enthalpies
    momentum: np.ndarray     # kg m/s (n, 3), at the cells' own velocities: the droplets' momentum, fact 29's sink
    solid: np.ndarray        # kg (n,): the solid in it, by the cells' liquid fractions (the droplets' solid fraction)
    by_phase: object = None  # the same by the cells it came from, a ReleaseBook of (2, ...) arrays: [0] fully liquid,
                             # [1] the others (sub-plan 22 Task 2)


def release_top(sk, dm, slurry=False):
    """Release dm [kg, (n,)] per column from its top cells (spec 8; sub-plan 21 Task 7, the placeholder for stage 4's
    spraying, whose rate sub-plan 22 supplies): liquid first -- the run of fully liquid cells from the top, at or above
    the liquidus -- and, with `slurry` (spec 8's `--slurry-spray on`), the exposed slurry too, every flowing cell down
    to the first rigid one; never rigid metal (spec 6.5: nothing passes through it), so a top that is slurry without
    the flag, a floating crust or anchored metal releases nothing. At most what those cells hold leaves; the order
    within them, and each one's rate, is the cascade rule's (sub-plan 22). The metal leaves at its cells' own specific
    enthalpy and velocity (`Skins.take_top_mass`): the droplets carry momentum away from the top cells at their release
    speed (spec 6.4), the speed of the metal they are made of, the top of the resolved profile -- fact 29's sink
    mdot V_s with V_s the sheared layer's surface speed -- so each column's mass, enthalpy and momentum fall by exactly
    what the book says left. Updates the cells; returns a ReleaseBook, with what came from the fully liquid cells and
    what from the others in `by_phase` (sub-plan 22 Task 2: the liquid's mode's droplets and the slurry's)."""
    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
    liquid = flowing & (sk.T >= sk.material.T_liquidus)
    admit = flowing if slurry else liquid
    mass, H, solid = sk.take_top_mass(dm, admit, phase=liquid)
    pm, pH, ps, pp = sk.phase_out
    return ReleaseBook(mass, H, sk.momentum_out.copy(), solid, ReleaseBook(pm, pH, pp, ps))


@dataclass
class Cascade:
    """A release given phase by phase (`Flow.release`; sub-plan 22 Task 1c): per column the liquid's rate and the exposed
    slurry's own [kg/(m2 s) of top area, (n,)], the slurry's holding for slurry cells at least as warm as `T_slurry` [K,
    (n,)], the cell whose properties gave it -- a colder slurry cell, more solid and more viscous, ends the column's
    cascade -- or for every slurry cell where None. `histogram`, if given, bins each flow sub-step's droplets by size:
    a callable of the release by phase [kg, (2, n): the liquid's mode's, the slurry's] giving (droplets, mass [kg]) per
    size bin (sub-plan 22 Task 2: `skinbody.SprayPlan.histogram`, at the macro step's radii)."""
    rate: np.ndarray
    rate_slurry: np.ndarray
    T_slurry: np.ndarray = None
    histogram: object = None


def cascade_mass(sk, rate, rate_slurry, A_top, dt, slurry=False, T_slurry=None):
    """The mass [kg, (n,)] the release takes from each column's top cells in dt [s] when each cell goes at its own phase's
    rate, and the part of it taken from slurry cells (spec 8's cascade, the 2026-10-07 design's 8.1; sub-plan 22 Task 1c):
    down the run of cells `release_top` may take -- the liquid cells from the top, at or above the liquidus, and with
    `slurry` every flowing cell down to the first rigid one -- in time order, a liquid cell at `rate` and a slurry cell at
    `rate_slurry` [kg/(m2 s) of top area, (n,)] times the top area A_top [m2], until dt is spent: the liquid first where
    it covers the slurry, the slurry once nothing liquid covers it, then whatever lies beneath. A cell whose rate is zero
    ends the column's cascade, as rigid metal does, and so does a slurry cell colder than `T_slurry` [K, (n,)], if given
    (`Cascade`). `release_top(sk, mass, slurry)` then takes exactly these cells."""
    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
    liquid = flowing & (sk.T >= sk.material.T_liquidus)
    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
    rate_s = np.asarray(rate_slurry, dtype=float)[:, None]
    if T_slurry is not None:
        rate_s = np.where(sk.T >= np.asarray(T_slurry, dtype=float)[:, None], rate_s, 0.0)
    r = np.where(liquid, np.asarray(rate, dtype=float)[:, None], rate_s)
    r = r * np.asarray(A_top, dtype=float)[:, None]                           # kg/s from each cell while it is exposed
    with np.errstate(divide="ignore", invalid="ignore"):
        lasts = np.where(run & (r > 0.0), sk.m / np.where(r > 0.0, r, 1.0), np.inf)    # s to strip each cell
    starts = np.concatenate([np.zeros((sk.n, 1)), np.cumsum(lasts, axis=1)[:, :-1]], axis=1)
    spare = np.where(run & (starts < dt), dt - starts, 0.0)
    taken = np.where(run, np.minimum(spare * r, sk.m), 0.0)
    return taken.sum(axis=1), (taken * ~liquid).sum(axis=1)


def prescribed_release(rate):
    """A release (`Flow.release`) at the rate given [kg/(m2 s), per column or one for all] whatever the skins' state:
    the tests' device, and any prescribed spraying's, which sub-plan 22 replaces with Girin's rate per column."""
    rate = np.array(rate, dtype=float)
    return lambda sk: np.broadcast_to(rate, (sk.n,))


def droplet_momentum_ratio(released, shear):
    """Fact 29's droplet-momentum ratio Pi = mdot V_s / tau_gas: the momentum the droplets take from the film against the
    gas shear that drives it, from the books of a time (a macro step in the history): released (n, 3) [kg m/s] the
    momentum the release took from each column's top cells (`FlowBook.released_p`), shear (n, 3) [N s] the gas shear's
    impulse on its top over the same time (`FlowBook.shear`), both integrated over the column's top and the time.
    Returns (Pi, each, flagged): each column's Pi_i = |released_i| / |shear_i|, zero where nothing was released and
    infinite where metal left a top the gas did not shear; Pi = sum |released_i| / sum |shear_i| over the sheared
    columns, the share of the gas shear's momentum that left with droplets -- each column's ratio weighted by its shear,
    so that columns barely sheared do not dominate it -- zero without a release; and flagged, the number of sheared
    columns above PI_FLAG, where fact 29 calls the stripping and the runoff no longer separable. Pi is defined against
    the shear, so a top the gas does not shear has none and is neither in Pi nor flagged: Girin's modes release only on
    the windward face, where the gas shears, but a prescribed release on the leeward runoff (moving at metres a second
    on the 100 mm sphere's rigid rear, under no shear) would otherwise make Pi infinite. On a layer at its steady
    Couette profile released at the top, Pi_i = mdot (H - dz_0 / 2) / mu: fact 29's v_melt delta_m / nu with delta_m the
    top cell's height above the rigid front.

    Measured 2026-10-24 on the 100 mm sphere's driven melting state (the front's top 0.3 mm liquid, driven from rest by
    the flight's forcing at 93 s for 20 ms), released at rho_l v_melt over the windward face for 10 ms more: Pi = 0.145
    at v_melt = 1 mm/s, 1 146 of the 4 705 sheared columns above 0.2 (the largest 0.26), and 0.015 at 0.1 mm/s, none.
    Fact 29 estimated 0.48 and 0.048, taking the film at its steady Couette speed tau delta_m / mu; the 0.3 mm layer,
    whose viscous time h^2 / nu is 0.17 s, has not reached it in 30 ms (its top at about 4.3 m/s under 64 Pa, against
    15 m/s)."""
    p, J = np.linalg.norm(released, axis=1), np.linalg.norm(shear, axis=1)
    sheared = J > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        each = np.where(p > 0.0, p / J, 0.0)
    Pi = float(p[sheared].sum() / J.sum()) if sheared.any() else 0.0
    return Pi, each, int((each[sheared] > PI_FLAG).sum())


# -- the base (Task 7) -------------------------------------------------------------------------------------------------
def _bottom_run(sk, T_rigid, T_limit=None):
    """The cells a hand-down may take (`hand_down`), from each column's bottom up: their temperatures, whether each is in
    the run of rigid cells (no warmer than T_limit, with one) contiguous from the bottom, its mass where it is, and
    whether the limit stops the run below a rigid cell."""
    n, K = sk.m.shape
    j = np.arange(K)[None, :]
    idx = np.clip(sk.count[:, None] - 1 - j, 0, K - 1)                 # the cells from the bottom up
    rows = np.arange(n)[:, None]
    T_up = sk.T[rows, idx]
    rigid = (j < sk.count[:, None]) & (T_up <= T_rigid)
    bound = np.zeros(n, dtype=bool)
    if T_limit is not None:
        free = np.cumprod(rigid, axis=1).sum(axis=1)
        rigid &= T_up <= np.asarray(T_limit, dtype=float)[:, None]
        bound = np.cumprod(rigid, axis=1).sum(axis=1) < free           # the limit stops the run below a rigid cell
    run = np.cumprod(rigid, axis=1).astype(bool)
    return T_up, run, np.where(run, sk.m[rows, idx], 0.0), bound


def handable(sk, m_min, T_rigid, T_limit=None):
    """The mass [kg per column] `hand_down` would take with these arguments, taking nothing. With no limit, less what a
    limited hand-down took, it is the metal the fault tolerance holds in the skin (sub-plan 21 Task 9; the skin body's
    held mass)."""
    m_up = _bottom_run(sk, T_rigid, T_limit)[2]
    return np.clip(np.minimum(m_up.sum(axis=1), sk.column_mass() - m_min), 0.0, None)


def hand_down(sk, m_min, T_rigid, T_limit=None, span=False, cap=None):
    """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
    oldest first, as far as the column stays at or above its minimum mass m_min (spec 6.6); a cell may go in part. With
    `T_limit` [K, per column], only cells no warmer than it go: the base node's temperature plus the fault tolerance
    (spec 6.6, amended 2026-10-26, sub-plan 21 Task 9), so that metal warmer than the node -- a runoff deposit whose heat
    is still conducting down -- waits in the skin until it has cooled, rather than leaving the node a slab 35-67 K/mm
    warmer than itself (Task 8's pushed plate) and its move, 1.1-3.2 move allowances a part, rejected. With `cap` [kg,
    per column], at most that much goes, a cell going in part where the cap ends inside it, and the rest waits: the
    skin body caps a part's hand-down at its move allowance (spec 6.6 and 7.1, amended 2026-10-27), because the limit
    releases a held cell whole once it is near the node's temperature, and on the full pushed plate at 1.5 s those cells
    were up to 0.79 mm of solid metal against a 0.556 mm allowance -- a release that does not shrink with the part, so
    parts were split in vain and the 12.4 g held came down in 0.039 s. A cap that does not bind changes nothing to the
    last bit. Returns (mass, enthalpy) handed down per column and, with `span`, the coldest and warmest temperatures [K]
    of the cells it took from, the one taken in part among them, nan where none went (the fault counter's range for a
    hand-down)."""
    n = sk.n
    T_up, run, m_up, bound = _bottom_run(sk, T_rigid, T_limit)
    take = np.clip(np.minimum(m_up.sum(axis=1), sk.column_mass() - m_min), 0.0, None)
    if cap is not None:                                                # at most the part's allowance (2026-10-27)
        take = np.minimum(take, np.asarray(cap, dtype=float))
    touched = run & (np.cumsum(m_up, axis=1) - m_up < take[:, None])  # the cells the take reaches, wholly or in part
    top = np.where(touched, m_up, 0.0)[np.arange(n), np.maximum(touched.sum(axis=1) - 1, 0)]
    H = sk.take_bottom(take)
    whole = bound & (take > 0.0) & (take >= m_up.sum(axis=1))         # the limit ends the take at a cell's top: what
    if whole.any():                                                    # the take's round-off leaves of that cell goes too
        rest = sk.m[np.arange(n), np.maximum(sk.count - 1, 0)]
        sliver = np.where(whole & (sk.count > 0) & (rest <= SLIVER * top), rest, 0.0)
        if sliver.any():
            momentum = sk.momentum_out
            H = H + sk.take_bottom(sliver)
            take = take + sliver
            sk.momentum_out = momentum + sk.momentum_out                # what both takes took away
    if not span:
        return take, H
    with np.errstate(invalid="ignore"):
        lo = np.where(touched.any(axis=1), np.where(touched, T_up, np.inf).min(axis=1), np.nan)
        hi = np.where(touched.any(axis=1), np.where(touched, T_up, -np.inf).max(axis=1), np.nan)
    return take, H, lo, hi
