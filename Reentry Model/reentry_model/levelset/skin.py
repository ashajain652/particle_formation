"""The level set's skins (design section 20.6, level-set stage 3 Task 2): columns of thin cells standing on the surface
markers, here only as a container of cells -- their storage, their grading with depth, their cross-section and the
regrid that keeps every cell near its target thickness. Conduction (Task 3), the layered flow (Task 4), the markers
(Task 5) and runoff (Task 6) build on it.

Storage. n columns side by side in padded arrays: cell 0 is a column's top (at the gas surface phi), cells
0 .. count-1 are real and every padded cell holds zeros. Per cell: mass `m` (n, K) [kg], enthalpy `H` (n, K) [J] above
the shared T_REF (m times the shared Material's specific enthalpy) and velocity `U` (n, K, 3) [m/s]; per column: the top
area `area` (n,) [m2], A0, and phi's principal curvatures `kappa` (n, 2) [1/m], positive where the surface is convex,
carried with the columns for the record (plan 03b Task 5 retired the curved cross-section that read them; the capillary
push reads the markers' own). Temperature, liquid fraction, density and phase follow from H/m through the material every
time they are asked for; phase labels are never stored (as the skyscrapers' skins).

Cross-section. Without a measured profile a column is flat: A(z) = A0 at every depth, and a cell of mass m at density
rho is m / (rho A0) thick. (Its curved form A0 (1 - k1 z)(1 - k2 z), with a cap on the curvatures, was retired by plan
03b Task 5: design section 21.5 measures every column's cross-section instead.)

The measured cross-section (design section 21.5, plan 03b Task 3). Every skin is one skin depth H deep
(`SkinSettings.depth`, 1.2 mm), a column of the layer between the gas surface and the skin base that the division
(`division.divide`) measures. `set_profile(levels, area)` gives each column that measure, its volume per unit depth at
the division's levels from 0 to H; from then on every place that takes
a cell's or an interface's area takes it from the profile: `cross_section` interpolates it linearly between the levels
and is zero below the last level with a nonzero area and beyond H (a column beside a convex edge pinches before H, where
it meets the column across the edge), so the volume from the top to depth z is piecewise quadratic, and `dz` inverts it
for each cell's bottom; the fill, the conduction's interfaces and base, the base area and the layered flow's interfaces
and floor follow. The top area `area` stays what the body gives: the area of the gas surface the column's marker owns,
the one the top fluxes act on (the gas heat, radiation, the gas shear and the spraying), which differs from the
profile's top by up to about 10 % at edge columns (the controller's ruling T1-g). A column's cells hold exactly its
volume after each fill; between a gain at the top (runoff, swelling on melting) and the next fill they may hold more,
and the cells below the column's end then take the larger of the end's cross-section and the column's mean
cross-section, continued (`cell_section`), so that every cell has a positive area and the excess lies no deeper below
the end than its share of the column's volume times the column's depth. (Plan 03b Task 5: with the end's section alone,
the controller's ruling T3-d, a column beside an edge whose base held 3.4e-6 of its top took 22 % excess 26 m below its
end on the test cube, and its regrid never settled.) Cells within END_SLACK of their column's volume
reach its end exactly (at a pinch the depth's sensitivity to the volume is infinite). The regrid splits a cell of a
profiled column at its mid-depth rather than in halves of its mass, so that where the cross-section closes a split never
leaves a cell small enough to be merged straight back. `exchange_base` draws metal from the volume below a column's
bottom cell or returns its bottom volume, `take_between` takes the metal between two depths, and `combine(...,
by="depth")` gives a column's metal to another at the same depths.

Grading and regrid (the skyscrapers' option B' rule, matched by amendment, never imported). Each cell's target
thickness is cell + (growth - 1) max(z_top - fine_zone, 0), z_top the depth of its top: 10 um (production) through the
top 0.4 mm, then the continuous form of "each cell 10 % thicker than the one above". `regrid` merges a cell below
MERGE_BELOW of its target into a neighbour and splits one above SPLIT_ABOVE in two, exactly in mass, enthalpy and
momentum, for at most REGRID_PASSES rounds. Merges pair cells (in a run of small cells every other one takes the cell
below it; a small bottom cell joins the cell above), so a run of thin cells coarsens by halving rather than collapsing
into one cell to be split again.

Every physical property comes from the material passed in (the shared `reentry_model.shared.material.Material` or a
test material with its interface), read through a `PropertyTable` of it built once per Skins (Task 4: h, cp_eff, k and
the liquid fraction on a 0.01 K grid; the cells' enthalpies and temperatures are the table's, so the books close on its
own h), and the viscosity from the shared rheology; this module's numbers are numerical settings only.

The layered flow (Task 4, the skyscrapers' momentum through the depth, matched by amendment). Phase labels follow from
T every time: a cell flows above the shared T_rigid (solid fraction below the rheology's RIGID_SOLID_FRACTION, a half),
is liquid at or above the liquidus, anchored in the unbroken run of rigid cells from the bottom up, and a crust where it
is rigid with flowing metal between it and the floor. `momentum_step` advances every column's velocities implicitly,
each component apart, with the drag between neighbouring cells through the cross-section at their interface's depth,
each flowing cell's viscosity at its own shear rate, the anchored cells at rest and a crust one unknown.

Sideways conduction (plan 03b Task 4, design 21.5, planning ruling 5 as the controller's ruling T6-b of Task 6 amends
it). Neighbouring columns share boundaries in the layer (`division.Division`'s i, j and length, `SideFaces`): a point on
a shared boundary lies at the same depth below phi on either side, so the cells that meet across it are those at that
depth in each column. Cell k of column i and cell l of column j exchange heat through the conductance
G_kl = Lambda_kl k_h / d_kl, Lambda_kl the shared length integrated over the overlap of the two cells' depth spans (the
length piecewise linear between the levels, and held at its last level's value below H, where only an overfull column's
cells reach, as `cell_section` holds the end's cross-section), k_h the harmonic mean 2 k_k k_l / (k_k + k_l) of their
conductivities and d_kl the distance between the two columns' centroids over that overlap (each column's first moment
integrated over it over its cross-section's, both linear between the levels: `SideFaces.section` and `moment`, the
division's), a two-point flux between the parts of the two columns that meet there. (Task 6: with the markers' distance
|x_j - x_i| instead, the stage 3 box's columns across an edge, whose centroids at 0.5 mm depth are 0.56 of their
markers' distance apart, conducted about half as well as their geometry does, and the probe column beside the edge, 0.60
to 1.84 times its neighbours' along the face, lost 9.3 % of its rise against the 3D solve's.) Where a column has no
cross-section over the overlap (beyond a pinch, where only an overfull column's cells reach), and for faces given with
no centroids, the markers' distance `SideFaces.distance` stands instead. `conduct(..., side=SideFaces)` puts these fluxes in the same backward-Euler step
as the depth conduction: one Newton iteration over every column's cells together, its Jacobian (the columns'
tridiagonal blocks and the sideways couplings, with their dk/dT terms) solved by restarted GMRES preconditioned by the
columns' own tridiagonal solves (`_Sideways.solve`). Why that solver (plan 03b Task 4, design 21.14 leaves it to the
builder; measured once on the stage 3 box's Jacobian, 132 columns, 5676 unknowns): a sparse LU (SuperLU) factorised it
in 0.115 s, its fill 1.9 million entries for 64 thousand, against 0.014 s for GMRES in 18 iterations; with sideways
conductances eight times larger (markers 0.25 mm apart) 0.031 s in 45. The LU's fill grows faster than the cells,
while an iteration costs one Thomas solve and one product, and the depth conduction, the stiff part, is solved exactly
by the preconditioner. The heat each column gains sideways is booked per column (`HeatBook.E_side`), each pair's flux
given by one and taken by the other, so it sums to zero to round-off whatever the solve's error."""
from collections import namedtuple
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from reentry_model.shared import rheology
from reentry_model.shared.thermal import SIGMA_SB

from .properties import PropertyTable

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (design 20.6; skyscrapers 6.2)
MERGE_BELOW = 0.5          # a cell thinner than this times its target joins a neighbour (skyscrapers' rule)
SPLIT_ABOVE = 1.5          # a cell thicker than this times its target is split in two (skyscrapers' rule)
REGRID_PASSES = 64         # merge-or-split rounds of one regrid at most (skyscrapers' rule)
NEWTON_TOL = 1.0e-6        # K: a heat sub-step has converged when no cell moves by more than this (skyscrapers' value)
NEWTON_MAX = 40            # Newton iterations of a heat sub-step at most (skyscrapers' value)
LINE_SEARCH_MAX = 12       # halvings of a Newton step before it is taken anyway (skyscrapers' value)
DZ_FLOOR = 1.0e-12         # m: the conduction (and drag) length of a vanishing cell
SHEAR_PASSES = 5           # passes of the fixed point on the cells' shear rates at most (skyscrapers', spec 6.4)
SHEAR_TOL = 1.0e-3         # a pass moving no flowing viscosity by more than this, relatively, ends it (skyscrapers')
CRUST_VISCOSITY = 1.0e9    # Pa s: a floating crust's half-cell at its faces with flowing metal, a numerical setting
                           # of the skin rule (skyscrapers' CRUST_VISCOSITY, sub-plan 21 Task 6b); the crust is one unknown
CUT_SLACK = 4.0 * np.finfo(float).eps   # fractions cut from a column summing to one within this empty it (`cut`)
FLOWING_SOLID_MAX = np.nextafter(rheology.RIGID_SOLID_FRACTION, 0.0)   # a flowing cell's solid fraction at most: guards
                           # the table's round-off at T_rigid, where the rheology turns infinite
END_SLACK = 1e-12          # cells within this share of their column's volume reach its end; shares of a cell below it
                           # (round-off where two depths nearly coincide) are none, and above one less it the whole
                           # (plan 03b Task 3: the profile's depths, `exchange_base`, `take_between`, by-depth combine)
LEVELS_RTOL = 1e-9         # a profile's deepest level must be the skins' depth to within this share (`set_profile`)
SIDE_RTOL = 1e-12          # the sideways heat's linear solve stops at this share of its right-hand side's norm (plan 03b
                           # Task 4; Newton's NEWTON_TOL still decides the sub-step's convergence)
SIDE_RESTART = 50          # GMRES iterations between restarts in that solve (plan 03b Task 4: on the stage 3 box it took
                           # 18, and 45 with the sideways conductances eight times larger, measured once each)
SIDE_CYCLES = 20           # restarts at most: SIDE_RESTART * SIDE_CYCLES iterations before it raises


@dataclass
class SkinSettings:
    cell: float = SKIN_PRESETS["production"]   # m, the fine cell
    fine_zone: float = 0.4e-3      # m below the top held in fine cells (design 20.6)
    growth: float = 1.10           # below the fine zone each cell this much thicker than the one above
    substep: float = 0.010         # s, the skins' heat sub-step (skyscrapers' spec 7.1)
    swell: bool = True             # each cell at its own phase's mixture density; False holds every cell at the solid's
    depth: float = 1.2e-3          # m, the skin depth H: every skin's column reaches this far below the gas surface
                                   # (design 21.3: today's 0.4 mm fine zone plus 0.8 mm stub; recorded with every run)

    def __post_init__(self):
        if not (self.cell > 0.0 and self.fine_zone > 0.0 and self.substep > 0.0 and self.depth > 0.0):
            raise ValueError("skin cell, fine zone, substep and depth must be positive")
        if self.growth < 1.0:
            raise ValueError("skin growth must be at least 1")
        if self.fine_zone > self.depth * (1.0 + 1e-12):
            raise ValueError("a skin's fine zone must lie within its depth: fine_zone <= depth")
        if 2.0 * self.cell > min(self.fine_zone, self.depth) * (1.0 + 1e-12):
            raise ValueError("a skin needs at least two fine cells: cell <= min(fine_zone, depth) / 2")


@dataclass
class HeatBook:
    """What one heat sub-step (`Skins.conduct`) did, per column [J]: the convective heat taken in at the top, the heat
    radiated from it, the heat conducted into the volume through the base (the residual of the books, so they close
    exactly whatever the iteration error), the heat gained sideways from the neighbouring columns (zero without
    `side`; with it, each pair's flux given by one column and taken by the other, so the columns' sum is zero to
    round-off) and, in a trial pass, the base heat's derivative with respect to the base temperature [J/K] (with
    sideways conduction, for every column's base temperature moving together, `Skins.conduct`)."""
    E_top: np.ndarray
    E_rad: np.ndarray
    E_base: np.ndarray
    E_side: np.ndarray
    dE_base_dTb: np.ndarray = None


@dataclass
class MomentumBook:
    """What one `Skins.momentum_step` did, per column [kg m/s] (n, 3): the impulses of the push and of the gas shear on
    the moving cells, and the momentum given to the wall -- the anchored metal and the volume beneath the base: the
    momentum the anchored cells held and the drag the moving cells lost to them -- so that the momentum change is
    push + shear - wall; and the shear-rate fixed point's passes. With response=True (Task 6, for the implicit thickness
    terms): every cell's velocity per unit push, `response` (n, K) [m/s per N/m3], from the last pass's system with the
    push alone as its right-hand side, and the drag of the anchored metal and the volume beneath on it, `wall_response`
    (n,) [m3] (the drag coefficients times the response), so that a change dG [N/m3] of the push changes the cells'
    momentum by dt (V - wall_response) dG, V the column's moving volume, the rest going to the wall."""
    push: np.ndarray
    shear: np.ndarray
    wall: np.ndarray
    passes: int
    response: np.ndarray = None
    wall_response: np.ndarray = None


@dataclass(frozen=True)
class BaseExchange:
    """What `Skins.exchange_base` moved across the base into each column it was given (plan 03b Task 3), in the order
    of its rows: the mass [kg], enthalpy above T_REF [J] and momentum [kg m/s] (N, 3), positive for a draw from the
    volume and negative for a return to it, so that each column's totals change by exactly these and the volume's by
    minus them; and `cells` (N,) int, the change in each column's count of cells: +1 for a draw (one new cell), minus
    the number of whole cells for a return (a cell it takes in part stays), 0 where nothing moved."""
    mass: np.ndarray
    H: np.ndarray
    P: np.ndarray
    cells: np.ndarray


SideFaces = namedtuple("SideFaces", "i j length levels distance section moment", defaults=(None, None))
SideFaces.__doc__ = """The boundaries neighbouring columns share in the layer, for sideways conduction (module docstring;
plan 03b Task 4): the two columns of each boundary `i`, `j` (E,) (skins' rows), its length at each level `length`
(E, B + 1) [m], the levels' depths `levels` (B + 1,) [m] from 0 to H (a `division.Division`'s, in the skins' row
order: `runoff.Sheet.side_faces`), the distance between the two columns' markers `distance` (E,) [m], |x_j - x_i|,
and (plan 03b Task 6, the controller's ruling T6-b) every column's cross-section `section` (n, B + 1) [m2] and first
moment `moment` (n, B + 1, 3) [m3] at the levels, in the skins' rows (the division's `area` and `moment`), whose
ratio over a depth span is the column's centroid there; without them the markers' distance stands."""


def level_values(levels, values, z):
    """Row by row, `values` (E, B + 1) given at the depths `levels` (B + 1,) [m] and linear between them, at the depths
    z (E, ...) [m]; held at the end values above the first level and below the last."""
    z = np.asarray(z, dtype=float)
    flat = z.reshape(len(values), int(np.prod(z.shape[1:])))
    k = np.clip(np.searchsorted(levels, flat, side="right") - 1, 0, len(levels) - 2)
    s = np.clip((flat - levels[k]) / (levels[k + 1] - levels[k]), 0.0, 1.0)
    out = (1.0 - s) * np.take_along_axis(values, k, axis=1) + s * np.take_along_axis(values, k + 1, axis=1)
    return out.reshape(z.shape)


def level_integral(levels, values, z):
    """Row by row, the integral of `level_values` from 0 to the depths z (E, ...) [m] (zero above 0): piecewise quadratic
    to the last level and linear below it, at the last value."""
    z = np.asarray(z, dtype=float)
    flat = z.reshape(len(values), int(np.prod(z.shape[1:])))
    h = np.diff(levels)
    W = np.zeros(values.shape)
    W[:, 1:] = np.cumsum(0.5 * (values[:, :-1] + values[:, 1:]) * h, axis=1)
    zc = np.clip(flat, levels[0], levels[-1])
    k = np.clip(np.searchsorted(levels, zc, side="right") - 1, 0, len(levels) - 2)
    t = zc - levels[k]
    v0, v1 = np.take_along_axis(values, k, axis=1), np.take_along_axis(values, k + 1, axis=1)
    inside = np.take_along_axis(W, k, axis=1) + t * (v0 + 0.5 * (v1 - v0) / h[k] * t)
    return (inside + values[:, -1:] * np.maximum(flat - levels[-1], 0.0)).reshape(z.shape)


class _Sideways:
    """The sideways conduction of one heat sub-step (module docstring; planning ruling 5, amended by ruling T6-b), between
    the cells of `side`'s columns as they lie at the sub-step's start: `z` (n, K + 1) each column's cell boundaries [m]
    below its top (held at its bottom beyond its count). For every pair of cells that overlap in depth across a shared
    boundary, the flat indices `a`, `b` of the cells (in i's and j's columns), their columns `ci`, `cj` and `lam` [m],
    the shared length integrated over the overlap over the distance between the columns' centroids there
    (`_distances`); with the Jacobian's sparsity pattern of the whole sub-step."""

    def __init__(self, side, z, count):
        n, K1 = z.shape
        K = K1 - 1
        i, j = np.asarray(side.i, dtype=np.int64), np.asarray(side.j, dtype=np.int64)
        length = np.asarray(side.length, dtype=float)
        levels = np.asarray(side.levels, dtype=float).reshape(-1)
        distance = np.asarray(side.distance, dtype=float).reshape(-1)
        E = len(i)
        if (j.shape != (E,) or length.shape != (E, len(levels)) or distance.shape != (E,)
                or (E and (min(i.min(), j.min()) < 0 or max(i.max(), j.max()) >= n))):
            raise ValueError("side faces need columns i, j (E,) of the skins, lengths (E, B + 1) at the levels and "
                             "distances (E,)")
        if not (np.all(np.diff(levels) > 0.0) and np.all(np.isfinite(length) & (length >= 0.0))
                and np.all(np.isfinite(distance) & (distance > 0.0))):
            raise ValueError("side faces need rising levels, finite non-negative lengths and positive distances")
        if side.section is not None and (np.shape(side.section) != (n, len(levels))
                                         or np.shape(side.moment) != (n, len(levels), 3)):
            raise ValueError("side faces' centroids need every column's cross-section (n, B + 1) and first moment "
                             "(n, B + 1, 3) at the levels")
        # the pieces both columns' cell boundaries cut the depth into: a piece starting at the p-th boundary in the
        # merged order lies in the cell of i whose top is the last of i's boundaries up to it (likewise for j)
        s = np.concatenate([z[i], z[j]], axis=1)
        order = np.argsort(s, axis=1, kind="stable")
        s = np.take_along_axis(s, order, axis=1)
        from_i = order < K1
        ki, kj = np.cumsum(from_i, axis=1)[:, :-1] - 1, np.cumsum(~from_i, axis=1)[:, :-1] - 1
        W = level_integral(levels, length, s)
        lam = np.diff(W, axis=1) / self._distances(side, i, j, levels, s, distance)
        keep = ((s[:, 1:] > s[:, :-1]) & (lam > 0.0) & (ki >= 0) & (kj >= 0)
                & (ki < count[i][:, None]) & (kj < count[j][:, None]))
        e, p = np.nonzero(keep)
        self.ci, self.cj, self.lam = i[e], j[e], lam[e, p]
        self.a, self.b = self.ci * K + ki[e, p], self.cj * K + kj[e, p]
        self.n, self.K = n, K
        idx = np.arange(n * K).reshape(n, K)
        self._rows = np.concatenate([idx.ravel(), idx[:, 1:].ravel(), idx[:, :-1].ravel(),
                                     self.a, self.a, self.b, self.b])
        self._cols = np.concatenate([idx.ravel(), idx[:, :-1].ravel(), idx[:, 1:].ravel(),
                                     self.a, self.b, self.a, self.b])

    @staticmethod
    def _distances(side, i, j, levels, s, distance):
        """Each overlap's distance [m] (E, P), the pieces between consecutive depths s (E, P + 1) of the boundaries
        between columns i and j (E,): the distance between the two columns' centroids over the piece (module docstring),
        or the markers' `distance` (E,) where either column has no cross-section over it, or `side` gives none."""
        out = np.broadcast_to(distance[:, None], (len(i), s.shape[1] - 1)).copy()
        if side.section is None or side.moment is None or not len(i):
            return out
        section = np.asarray(side.section, dtype=float)
        moment = np.asarray(side.moment, dtype=float)
        centre = []
        for c in (i, j):
            A = np.diff(level_integral(levels, section[c], s), axis=1)
            M = np.stack([np.diff(level_integral(levels, moment[c][..., x], s), axis=1) for x in range(3)], axis=-1)
            with np.errstate(divide="ignore", invalid="ignore"):
                centre.append(np.where((A > 0.0)[..., None], M / A[..., None], np.nan))
        d = np.linalg.norm(centre[1] - centre[0], axis=-1)
        ok = np.isfinite(d) & (d > 0.0)
        out[ok] = d[ok]
        return out

    def fluxes(self, T, k, kp):
        """At the cells' temperatures T, conductivities k and their slopes kp (n, K): each pair's flux F [W] from cell
        a to cell b, and its derivatives with respect to T_a and T_b [W/K]."""
        Ta, Tb, ka, kb = T.ravel()[self.a], T.ravel()[self.b], k.ravel()[self.a], k.ravel()[self.b]
        s = ka + kb
        G = self.lam * 2.0 * ka * kb / s
        dT = Ta - Tb
        F = G * dT
        dFa = G + dT * self.lam * 2.0 * kb * kb / (s * s) * kp.ravel()[self.a]
        dFb = -G + dT * self.lam * 2.0 * ka * ka / (s * s) * kp.ravel()[self.b]
        return F, dFa, dFb

    def residual(self, F):
        """The fluxes' share of the residual (n, K): out of each pair's a, into its b."""
        N = self.n * self.K
        return (np.bincount(self.a, F, N) - np.bincount(self.b, F, N)).reshape(self.n, self.K)

    def gained(self, F):
        """The heat rate each column gains sideways [W] (n,) from the pairs' quantities F (fluxes, or their changes)."""
        return np.bincount(self.cj, F, self.n) - np.bincount(self.ci, F, self.n)

    def solve(self, lower, diag, upper, jac, rhs):
        """x (n, K) with J x = rhs, J the sub-step's whole Jacobian: the columns' tridiagonal rows (lower, diag, upper,
        as `_tridiagonal` takes them) and the pairs' couplings jac = (dFa, dFb). Restarted GMRES preconditioned by the
        columns' own tridiagonal solves (block Jacobi: the depth conduction exactly, the sideways couplings left to
        the iteration), to SIDE_RTOL of the right-hand side; RuntimeError if it has not converged within SIDE_RESTART
        times SIDE_CYCLES iterations (module docstring)."""
        dFa, dFb = jac
        n, K = self.n, self.K
        N = n * K
        b = rhs.ravel()
        if not np.any(b):
            return np.zeros((n, K))
        data = np.concatenate([diag.ravel(), lower[:, 1:].ravel(), upper[:, :-1].ravel(), dFa, dFb, -dFa, -dFb])
        J = sp.csr_matrix((data, (self._rows, self._cols)), shape=(N, N))
        M = spla.LinearOperator((N, N), matvec=lambda v: _tridiagonal(lower, diag, upper, v.reshape(n, K)).ravel(),
                                dtype=float)
        x, info = spla.gmres(J, b, M=M, rtol=SIDE_RTOL, atol=0.0, restart=SIDE_RESTART, maxiter=SIDE_CYCLES)
        if info != 0:
            raise RuntimeError("the skins' sideways heat solve did not converge in {} GMRES iterations (residual {:.2e} "
                               "of the right-hand side, against {:.0e})".format(
                                   SIDE_RESTART * SIDE_CYCLES, np.linalg.norm(b - J @ x) / np.linalg.norm(b),
                                   SIDE_RTOL))
        return x.reshape(n, K)


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


class Skins:
    """n columns of cells (module docstring). `material` gives every property; `area` (n,) the top areas [m2], the
    area the top fluxes act on; `kappa` (n, 2) phi's principal curvatures [1/m], positive convex, carried for the
    record (no cross-section reads them, module docstring). Each column is flat until a measured profile is set
    (`set_profile`). The columns start empty (`fill` or `install`)."""

    profile_levels = None      # (B + 1,) m, the measured profile's levels from 0 to H, None until `set_profile`; class
    profile = None             # (n, B + 1) m2 its cross-sections; defaults, so that a Skins built by __new__ (runoff's
                               # landing blocks) or pickled before plan 03b Task 3 reads as having none

    def __init__(self, material, area, kappa, settings=None):
        self.material = material
        self.settings = settings if settings is not None else SkinSettings()
        self.area = np.array(area, dtype=float).reshape(-1)
        self.kappa = np.array(kappa, dtype=float).reshape(len(self.area), 2)
        self.rho_s = float(material.rho)
        liquid = getattr(material, "liquid", None)
        self.rho_l = float(liquid.rho) if liquid is not None else self.rho_s   # a material without a liquid never swells
        self.props = PropertyTable(material)      # every h, cp_eff, k and f_l the skins read (module docstring)
        self.dH = None             # (n, K) J/K, each cell's dH/dT_b while a trial pass is open (`conduct`), else None
        self.momentum_out = np.zeros((self.n, 3))   # kg m/s, what the last `remove_feed` took with the metal
        self.flow_memory = np.zeros(3)  # the last heat sub-step's largest Courant number, Rayleigh-Taylor growth and flow
                                        # sub-step [s] (Task 6, `runoff.flow_cycles`); zeros: no flow remembered
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
        in equal cells of about the fine size by the flat thickness; with a measured profile, in cells of equal depth
        about the fine size, each holding the metal its span holds (a span of no cross-section holds no cell)."""
        mass = np.broadcast_to(np.asarray(mass, dtype=float), (self.n,))
        T = np.broadcast_to(np.asarray(T, dtype=float), (self.n,))
        if self.profile is not None:
            self._fill_profile(mass, T)
            return
        flat = mass / (self.density(T) * self.area)
        cells = np.maximum(np.rint(flat / self.settings.cell).astype(np.int64), 1)
        K = int(cells.max()) if self.n else 1
        valid = np.arange(K)[None, :] < cells[:, None]
        m = np.where(valid, (mass / cells)[:, None], 0.0)
        self.install(m, m * self.props.enthalpy(T)[:, None], cells)

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
        """Cell temperatures [K] (n, K) from the specific enthalpies through the property table; nan beyond each
        count."""
        valid = self.valid
        T = np.full(self.m.shape, np.nan)
        if valid.any():
            T[valid] = self.props.temperature_from_enthalpy(self.H[valid] / self.m[valid])
        return T

    def density(self, T):
        """Mixture density [kg/m3]: specific volumes add, 1 / ((1 - f_l)/rho_s + f_l/rho_l), f_l the material's liquid
        fraction; rho_s everywhere with `swell=False`."""
        T = np.asarray(T, dtype=float)
        if not self.settings.swell:
            return np.full(T.shape, self.rho_s)
        f = self.props.liquid_fraction(T)
        return 1.0 / ((1.0 - f) / self.rho_s + f / self.rho_l)

    def cross_section(self, z):
        """The cross-section [m2] at depths z [m] below each column's top, z (n, ...): the top area at every depth
        (flat), or the measured profile once set, linear between its levels, zero below its last nonzero level and
        beyond H."""
        z = np.asarray(z, dtype=float)
        if self.profile is not None:
            return self._section(z)
        shape = (self.n,) + (1,) * (z.ndim - 1)
        return np.broadcast_to(self.area.reshape(shape), z.shape).copy()

    def dz(self, T=None):
        """Cell thicknesses [m] (n, K), 0 beyond each count: each cell's span holds its volume m/rho(T) under the
        column's cross-section (module docstring): the flat thickness m / (rho A0), or with a measured profile the
        inverse of its piecewise-quadratic volume (`cell_section` below the column's end), and ValueError for metal in
        a column with neither a profile volume nor a top area (`_depth_of`)."""
        valid = self.valid
        T = self.T() if T is None else T
        if self.profile is not None:
            return self._profile_dz(T)
        flat = np.zeros(self.m.shape)
        flat[valid] = self.m[valid] / (self.density(T[valid]) * np.broadcast_to(self.area[:, None], self.m.shape)[valid])
        return flat

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

    # -- the measured cross-section (plan 03b Task 3) ---------------------------------------------------------------
    def set_profile(self, levels, area):
        """Every column's cross-section from now on the measured profile `area` (n, B + 1) [m2], each column's volume
        per unit depth at the depths `levels` (B + 1,) [m], rising from 0 to the skins' depth H (the division's
        `levels` and `area`, plan 03b Task 1), linear between them (module docstring). The top area `area` the top fluxes act on is left as it is (the controller's ruling T1-g). The rows follow
        the columns through `divide`, `cut`, `combine` and `remove`."""
        levels = np.array(levels, dtype=float).reshape(-1)
        area = np.array(area, dtype=float)
        if len(levels) < 2 or levels[0] != 0.0 or np.any(np.diff(levels) <= 0.0):
            raise ValueError("a profile's levels must rise strictly from 0")
        if abs(levels[-1] - self.settings.depth) > LEVELS_RTOL * self.settings.depth:
            raise ValueError("a profile's deepest level ({:.6g} m) must be the skins' depth ({:.6g} m)".format(
                levels[-1], self.settings.depth))
        if area.shape != (self.n, len(levels)) or not np.all(np.isfinite(area) & (area >= 0.0)):
            raise ValueError("a profile needs one finite, non-negative cross-section per column and level, (n, B + 1)")
        self.profile_levels, self.profile = levels, area

    def column_volume(self):
        """Each column's volume [m3] (n,): the trapezoid integral of its measured profile (the division's own); without
        one, its top area times H (flat)."""
        if self.profile is not None:
            return self._ends()[2]
        return self.area * self.settings.depth

    def cells_volume(self, T=None):
        """The volume [m3] (n,) each column's cells fill, the sum of m / rho(T) over its cells (T (n, K) their
        temperatures if already known)."""
        valid = self.valid
        T = self.T() if T is None else T
        vol = np.zeros(self.m.shape)
        vol[valid] = self.m[valid] / self.density(T[valid])
        return vol.sum(axis=1)

    def cell_section(self, z):
        """The cross-section [m2] the cells take at depths z [m] (n, ...): `cross_section` within each column, and below
        the column's end -- where only an overfull column's cells reach, between a gain at its top and the next fill --
        the larger of the end's cross-section and the column's mean, continued (a pinch's is the mean), or its top area
        where it has no profile at all (module docstring); for the cells' thicknesses and interfaces, not the
        base (`base_area`). Without a profile, `cross_section` itself."""
        z = np.asarray(z, dtype=float)
        return self._section(z, extend=True) if self.profile is not None else self.cross_section(z)

    def _rows(self, rows):
        return np.arange(self.n) if rows is None else np.asarray(rows, dtype=np.int64).reshape(-1)

    def _ends(self, rows=None):
        """For the columns `rows` (all by default), from the profile: the volumes W (r, B + 1) [m3] from the top to each
        level, each column's end z_end [m] (H, or the level after its last nonzero cross-section; 0 for none), its
        volume V [m3], its cross-section at the end A_end and the one below it, A_ext [m2] (`cell_section`): the larger
        of A_end and the column's mean cross-section V / z_end (its top area for a column with no profile volume)."""
        rows = self._rows(rows)
        P, lv = self.profile[rows], self.profile_levels
        B = len(lv) - 1
        W = np.zeros(P.shape)
        W[:, 1:] = np.cumsum(0.5 * (P[:, :-1] + P[:, 1:]) * np.diff(lv), axis=1)
        nonzero = P > 0.0
        last = np.where(nonzero.any(axis=1), B - np.argmax(nonzero[:, ::-1], axis=1), -1)
        end = np.clip(last + 1, 0, B)
        z_end, V = lv[end], W[:, -1]
        A_end = P[np.arange(len(rows)), end]
        mean = np.divide(V, z_end, out=np.zeros_like(V), where=z_end > 0.0)
        A_ext = np.where(z_end > 0.0, np.maximum(A_end, mean), self.area[rows])
        return W, z_end, V, A_end, A_ext

    def _section(self, z, rows=None, extend=False):
        """The profile's cross-section [m2] at depths z [m] (r, ...) of the columns `rows`; with `extend`, the cells'
        (`cell_section`)."""
        rows = self._rows(rows)
        flat = z.reshape(len(rows), int(np.prod(z.shape[1:])))
        P, lv = self.profile[rows], self.profile_levels
        k = np.clip(np.searchsorted(lv, flat, side="right") - 1, 0, len(lv) - 2)
        s = np.clip((flat - lv[k]) / (lv[k + 1] - lv[k]), 0.0, 1.0)
        out = (1.0 - s) * np.take_along_axis(P, k, axis=1) + s * np.take_along_axis(P, k + 1, axis=1)
        out = np.where(flat > lv[-1], 0.0, out)
        if extend:
            _, z_end, _, _, A_ext = self._ends(rows)
            out = np.where(flat > z_end[:, None], A_ext[:, None], out)
        return out.reshape(z.shape)

    def _volume_at(self, z, rows=None):
        """The volume [m3] from the top of each of the columns `rows` down to the depths z [m] (r, ...), under the
        cells' cross-section (`cell_section`; without a profile flat, the top area at every depth); 0 above the
        top."""
        rows = self._rows(rows)
        z = np.asarray(z, dtype=float)
        flat = z.reshape(len(rows), int(np.prod(z.shape[1:])))
        if self.profile is None:
            return (self.area[rows][:, None] * np.maximum(flat, 0.0)).reshape(z.shape)
        W, z_end, V, _, A_ext = self._ends(rows)
        P, lv = self.profile[rows], self.profile_levels
        k = np.clip(np.searchsorted(lv, flat, side="right") - 1, 0, len(lv) - 2)
        h = np.diff(lv)[k]
        t = np.clip(flat - lv[k], 0.0, h)
        Ak, Ak1 = np.take_along_axis(P, k, axis=1), np.take_along_axis(P, k + 1, axis=1)
        inside = np.take_along_axis(W, k, axis=1) + t * (Ak + 0.5 * (Ak1 - Ak) / h * t)
        beyond = V[:, None] + A_ext[:, None] * (flat - z_end[:, None])
        return np.where(flat > z_end[:, None], beyond, inside).reshape(z.shape)

    def _depth_of(self, v, rows=None):
        """The depths [m] (r, ...) holding the volumes v [m3] (r, ...) from the top of each of the columns `rows` under
        the profile: in the band of levels holding v, the root of the quadratic W_k + A_k t + (A_k+1 - A_k) t^2 / 2h, in
        the form 2 r / (A_k + sqrt(A_k^2 + 2 s r)) that stays exact where the cross-section closes; at the column's end
        for v within END_SLACK of its volume, and below it at `cell_section`'s cross-section beyond. ValueError for a
        volume to place in a column with neither a profile volume nor a top area (no cross-section anywhere: Task 3 fix
        round 1, where it gave NaN thicknesses)."""
        rows = self._rows(rows)
        flat = np.asarray(v, dtype=float).reshape(len(rows), int(np.prod(np.shape(v)[1:])))
        W, z_end, V, _, A_ext = self._ends(rows)
        nowhere = np.any(flat > V[:, None], axis=1) & (A_ext <= 0.0)
        if nowhere.any():
            raise ValueError("{} column(s) hold metal but have neither a profile volume nor a top area to hold it "
                             "(rows {})".format(int(nowhere.sum()), rows[nowhere][:8].tolist()))
        P, lv = self.profile[rows], self.profile_levels
        k = np.zeros(flat.shape, dtype=np.int64)
        for level in range(1, len(lv) - 1):                    # the band: the levels whose volume v reaches
            k += flat >= W[:, level, None]
        h = np.diff(lv)[k]
        Ak, Ak1 = np.take_along_axis(P, k, axis=1), np.take_along_axis(P, k + 1, axis=1)
        r = np.maximum(flat - np.take_along_axis(W, k, axis=1), 0.0)
        den = Ak + np.sqrt(np.maximum(Ak * Ak + 2.0 * (Ak1 - Ak) / h * r, 0.0))
        z = lv[k] + np.minimum(np.divide(2.0 * r, den, out=np.zeros_like(r), where=den > 0.0), h)
        excess = np.maximum(flat - V[:, None], 0.0)
        A = np.broadcast_to(A_ext[:, None], flat.shape)
        below = np.divide(excess, A, out=np.zeros_like(excess), where=A > 0.0)
        full = flat >= V[:, None] * (1.0 - END_SLACK)
        return np.where(full, z_end[:, None] + below, z).reshape(np.shape(v))

    def _profile_dz(self, T, rows=None):
        """`dz` under the profile for the columns `rows` (all by default), at the temperatures T (n, K): each cell's
        bottom the depth holding the volume m / rho(T) of it and all above (`_depth_of`), less its top."""
        rows = self._rows(rows)
        valid = self.valid[rows]
        vol = np.zeros(valid.shape)
        vol[valid] = self.m[rows][valid] / self.density(T[rows][valid])
        z_bot = self._depth_of(np.cumsum(vol, axis=1), rows)
        z_top = np.concatenate([np.zeros((len(rows), 1)), z_bot[:, :-1]], axis=1)
        return np.where(valid, z_bot - z_top, 0.0)

    def _base_section(self, T, rows=None):
        """The cross-section [m2] at the bottom of each of the columns `rows` under the profile, from its cells' volume
        at the temperatures T (n, K): the profile's at the column's end, A(H), where they hold its volume to END_SLACK
        or more (none for a column that pinches before H: it has no base, design 21.5, even while overfull); the
        profile's at their depth where they hold less (A(0) for an empty column). Judged by volume, not by the summed
        thicknesses, which round-off leaves either side of H."""
        rows = self._rows(rows)
        valid = self.valid[rows]
        Tr = T[rows]
        vol = np.zeros(valid.shape)
        vol[valid] = self.m[rows][valid] / self.density(Tr[valid])
        v = vol.sum(axis=1)
        _, _, V, A_end, _ = self._ends(rows)
        A = self._section(self._depth_of(v[:, None], rows), rows)[:, 0]
        return np.where(v >= V * (1.0 - END_SLACK), A_end, A)

    def _fill_profile(self, mass, T):
        """`fill` under the profile: the depth holding mass / rho(T) cut into equal spans of about the fine size, each
        cell the mass its span holds (the shares of the volume, so they sum to `mass`); cells of no mass are dropped."""
        depth = self._depth_of((mass / self.density(T))[:, None])[:, 0]
        cells = np.maximum(np.rint(depth / self.settings.cell).astype(np.int64), 1)
        K = int(cells.max()) if self.n else 1
        edges = depth[:, None] * np.minimum(np.arange(K + 1)[None, :], cells[:, None]) / cells[:, None]
        V = self._volume_at(edges)
        share = np.divide(np.maximum(np.diff(V, axis=1), 0.0), V[:, -1:], out=np.zeros((self.n, K)),
                          where=V[:, -1:] > 0.0)
        m = mass[:, None] * share
        order = np.argsort(m <= 0.0, axis=1, kind="stable")          # the cells with metal first, in their order
        m = np.take_along_axis(m, order, axis=1)
        self.install(m, m * self.props.enthalpy(T)[:, None], (m > 0.0).sum(axis=1))

    # -- the base exchange and the metal between two depths (plan 03b Task 3) ----------------------------------------
    def exchange_base(self, rows, dV, T_draw=None):
        """Metal crosses the base of each column rows[j] (distinct), between it and the volume beneath (design 21.3:
        every skin keeps its column full by drawing at its base and returning there). dV[j] > 0 [m3] draws that volume
        of the volume's metal, at the solid's density rho_s (the volume's one density, design 21.6) and at T_draw[j]
        [K], as a new cell at rest below the bottom cell. dV[j] < 0 returns the bottom -dV[j] of the column's volume
        (each cell's m / rho(T)) whatever its phase: the whole cells it covers, from the bottom up, and the share left
        of the one cell above them, which keeps the rest at its own specific enthalpy and velocity. T_draw (one or per
        row) is needed only where a column draws. Returns a BaseExchange: what crossed into each column, signed
        (negative for a return), and the change in its count of cells. The caller regrids."""
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        dV = np.broadcast_to(np.asarray(dV, dtype=float), rows.shape)
        if len(np.unique(rows)) != len(rows) or not np.all(np.isfinite(dV)):
            raise ValueError("a base exchange needs distinct columns and finite volumes")
        self._no_trial("exchange metal at their bases")
        draw, back = np.flatnonzero(dV > 0.0), np.flatnonzero(dV < 0.0)
        if len(draw):
            if T_draw is None:
                raise ValueError("a draw at a skin's base needs the drawn metal's temperature T_draw")
            T_d = np.broadcast_to(np.asarray(T_draw, dtype=float), rows.shape)[draw]
            if not np.all(np.isfinite(T_d)):
                raise ValueError("a draw at a skin's base needs a finite T_draw")
        N = len(rows)
        mass, H, P, cells = np.zeros(N), np.zeros(N), np.zeros((N, 3)), np.zeros(N, dtype=np.int64)
        if len(back):
            r, want = rows[back], -dV[back]
            T, valid = self.T()[r], self.valid[r]
            vol = np.zeros(valid.shape)
            vol[valid] = self.m[r][valid] / self.density(T[valid])
            if np.any(want > vol.sum(axis=1) * (1.0 + END_SLACK)):
                raise ValueError("a return at a skin's base can take at most the volume its cells hold")
            below = np.cumsum(vol[:, ::-1], axis=1)[:, ::-1]                  # each cell's volume and all beneath it
            whole = valid & (below <= want[:, None] * (1.0 + END_SLACK))
            n_whole = whole.sum(axis=1)
            take = whole.astype(float)
            j = self.count[r] - 1 - n_whole                                    # the cell above the whole ones
            part = np.flatnonzero(j >= 0)
            rest = want[part] - (vol[part] * whole[part]).sum(axis=1)
            f = np.clip(rest / vol[part, j[part]], 0.0, 1.0)
            take[part, j[part]] = np.where(f > END_SLACK, f, 0.0)
            m_out = self.m[r] * take
            mass[back], H[back] = -m_out.sum(axis=1), -(self.H[r] * take).sum(axis=1)
            P[back] = -(m_out[..., None] * self.U[r]).sum(axis=1)
            kept = 1.0 - take[part, j[part]]
            self.m[r[part], j[part]] *= kept
            self.H[r[part], j[part]] *= kept
            self._compact(r, whole)
            cells[back] = -n_whole
        if len(draw):
            m_d = self.rho_s * dV[draw]
            H_d = m_d * self.props.enthalpy(T_d)
            self._append_bottom(rows[draw], m_d, H_d)
            mass[draw], H[draw], cells[draw] = m_d, H_d, 1
        self._trim()
        return BaseExchange(mass, H, P, cells)

    def _append_bottom(self, rows, m, H, U=None):
        """One new cell below the bottom of each column rows[j] (distinct), of mass m[j] > 0 [kg], enthalpy H[j] [J]
        and velocity U[j] (3,) [m/s] (at rest by default); the caller regrids. (A draw's cell, `exchange_base`, and the
        metal a column with no base takes from its neighbour, the body's fill.)"""
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        m = np.broadcast_to(np.asarray(m, dtype=float), rows.shape)
        H = np.broadcast_to(np.asarray(H, dtype=float), rows.shape)
        if len(np.unique(rows)) != len(rows) or np.any(~np.isfinite(m) | (m <= 0.0)) or np.any(~np.isfinite(H)):
            raise ValueError("a new bottom cell needs distinct columns and a finite, positive mass with finite enthalpy")
        U = np.zeros((len(rows), 3)) if U is None else np.broadcast_to(np.asarray(U, dtype=float), (len(rows), 3))
        self._width(int((self.count[rows] + 1).max(initial=1)))
        c = self.count[rows]
        self.m[rows, c], self.H[rows, c], self.U[rows, c] = m, H, U
        self.count[rows] = c + 1

    def take_between(self, rows, z_top, z_bottom):
        """Remove the metal lying between the depths z_top and z_bottom [m] (one or per row, z_top <= z_bottom) below
        the tops of the columns rows (distinct): each cell's share of its volume between them, by the volume its span
        holds under the column's cross-section (`cell_section`). Wholly covered cells leave, the cells the depths cut
        keep the rest at their own specific enthalpy and velocity, and the cells beneath move up. Returns per row the
        (mass [kg], enthalpy [J], momentum [kg m/s] (., 3)) taken; the caller regrids."""
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        z_top = np.broadcast_to(np.asarray(z_top, dtype=float), rows.shape)
        z_bottom = np.broadcast_to(np.asarray(z_bottom, dtype=float), rows.shape)
        if len(np.unique(rows)) != len(rows) or not np.all(np.isfinite(z_top) & (z_bottom >= z_top)):
            raise ValueError("take_between needs distinct columns and finite depths z_top <= z_bottom")
        self._no_trial("take metal between depths")
        T, valid = self.T()[rows], self.valid[rows]
        vol = np.zeros(valid.shape)
        vol[valid] = self.m[rows][valid] / self.density(T[valid])
        v1 = np.cumsum(vol, axis=1)
        Va = self._volume_at(np.maximum(z_top, 0.0)[:, None], rows)
        Vb = self._volume_at(np.maximum(z_bottom, 0.0)[:, None], rows)
        overlap = np.clip(np.minimum(v1, Vb) - np.maximum(v1 - vol, Va), 0.0, vol)
        f = np.divide(overlap, vol, out=np.zeros_like(vol), where=vol > 0.0)
        f = np.where(f >= 1.0 - END_SLACK, 1.0, np.where(f > END_SLACK, f, 0.0))
        m_out = self.m[rows] * f
        out = m_out.sum(axis=1), (self.H[rows] * f).sum(axis=1), (m_out[..., None] * self.U[rows]).sum(axis=1)
        self.m[rows] = self.m[rows] * (1.0 - f)
        self.H[rows] = self.H[rows] * (1.0 - f)
        self._compact(rows, valid & (f == 1.0))
        self._trim()
        return out

    def _compact(self, rows, drop):
        """Remove the cells `drop` (len(rows), K) from the columns `rows` (distinct); the cells beneath move up, in
        their order. No trial pass may be open."""
        if not drop.any():
            return
        order = np.argsort(drop, axis=1, kind="stable")                    # the cells kept first, in their order
        count = self.count[rows] - drop.sum(axis=1)
        keep = np.arange(drop.shape[1])[None, :] < count[:, None]
        self.m[rows] = np.where(keep, np.take_along_axis(self.m[rows], order, axis=1), 0.0)
        self.H[rows] = np.where(keep, np.take_along_axis(self.H[rows], order, axis=1), 0.0)
        self.U[rows] = np.where(keep[..., None], np.take_along_axis(self.U[rows], order[..., None], axis=1), 0.0)
        self.count[rows] = count

    # -- regrid -----------------------------------------------------------------------------------------------------
    def regrid(self):
        """Merge cells below MERGE_BELOW of their target, else split those above SPLIT_ABOVE, until none qualifies.
        Under a measured profile a cell splits at its mid-depth (module docstring), unless one of its halves would be
        empty or small enough to merge straight back (`_splits_settle`): then it is left whole."""
        for _ in range(REGRID_PASSES):
            valid = self.valid
            dz = self.dz()
            t = self.targets(dz)
            small = valid & (dz < MERGE_BELOW * t) & (self.count >= 2)[:, None]
            big = valid & (dz > SPLIT_ABOVE * t)
            if small.any():
                self._merge(np.flatnonzero(small.any(axis=1)), small)
                continue
            rows = np.flatnonzero(big.any(axis=1))
            upper = None
            if len(rows) and self.profile is not None:
                upper = self._upper_halves(rows, dz)
                big[rows] &= self._splits_settle(rows, dz, upper)
                some = big[rows].any(axis=1)
                rows, upper = rows[some], upper[some]
            if not len(rows):
                return
            self._split(rows, big, upper)
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

    def _splits_settle(self, rows, dz, upper):
        """Which cells (len(rows), K) of a profiled column may split at their mid-depth: those whose two halves, by
        the depths their volumes would hold (`_depth_of`), are each at least MERGE_BELOW of their targets. Beside an edge
        a column's measured cross-section is jagged and may be zero over a band (plan 03b Task 5, on the test cube); a
        cell over such a band split at its mid-depth leaves one half within the band, small, merged straight back, and
        the regrid would never settle; or an empty half (`upper` (len(rows), K), each cell's share above its
        mid-depth, at 0 or 1). A mid-depth inside such a band is never a split's (the depth a volume reaches there is
        either end of the band, by round-off)."""
        d = dz[rows]
        z0 = np.cumsum(d, axis=1) - d
        V0, V1 = self._volume_at(z0, rows), self._volume_at(z0 + d, rows)
        zs = self._depth_of(V0 + upper * (V1 - V0), rows)
        s = self.settings
        target = lambda z: s.cell + (s.growth - 1.0) * np.maximum(z - s.fine_zone, 0.0)
        inside = self._section(z0 + 0.5 * d, rows, extend=True) > 0.0      # a mid-depth within a band of no
        return (inside & (upper > END_SLACK) & (upper < 1.0 - END_SLACK)        # cross-section is no place to split
                & (zs - z0 >= MERGE_BELOW * target(z0)) & (z0 + d - zs >= MERGE_BELOW * target(zs)))

    def _upper_halves(self, rows, dz):
        """Each cell's share (len(rows), K) of its volume above its mid-depth, under the profile, from the thicknesses
        dz (n, K); one half for a cell of no volume."""
        d = dz[rows]
        z0 = np.cumsum(d, axis=1) - d
        V0, Vm, V1 = (self._volume_at(z, rows) for z in (z0, z0 + 0.5 * d, z0 + d))
        return np.divide(Vm - V0, V1 - V0, out=np.full(d.shape, 0.5), where=V1 > V0)

    def _split(self, rows, big, upper=None):
        """Every big cell in `rows` becomes two, both at its velocity: halves of its mass and enthalpy, or with `upper`
        (len(rows), K) that share of each for the upper one and the rest for the lower (a profile's mid-depth split)."""
        valid = self.valid[rows]
        reps = np.where(valid, 1 + big[rows], 0)
        nr, K = reps.shape
        count = reps.sum(axis=1)
        self._width(int(count.max()))
        K2 = self.m.shape[1]
        dest = np.cumsum(reps, axis=1) - reps
        rr = np.broadcast_to(np.arange(nr)[:, None], (nr, K))
        share = 1.0 / np.maximum(reps, 1)
        if upper is None:
            shares = (share, share)
        else:
            shares = (np.where(reps == 2, upper, 1.0), np.where(reps == 2, 1.0 - upper, 1.0))
        m, H, U = np.zeros((nr, K2)), np.zeros((nr, K2)), np.zeros((nr, K2, 3))
        U0 = self.U[rows, :K]
        trial = self.dH is not None
        if trial:
            dH = np.zeros((nr, K2))
        for half in (0, 1):
            sel = valid & (half < reps)
            f = shares[half]
            m[rr[sel], dest[sel] + half] = (self.m[rows, :K] * f)[sel]
            H[rr[sel], dest[sel] + half] = (self.H[rows, :K] * f)[sel]
            U[rr[sel], dest[sel] + half] = U0[sel]
            if trial:
                dH[rr[sel], dest[sel] + half] = (self.dH[rows, :K] * f)[sel]
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

    def take_top_mass(self, dm, admit, phase=None):
        """Remove dm >= 0 [kg] (n,) per column from its top down through the run of admissible cells from the top
        (`admit`, (n, K) bool), at most what that run holds (stage 3 Task 8, the release; the skyscrapers' rule, matched
        by amendment). Each cell's metal leaves at its own specific enthalpy and velocity: a partly taken cell keeps
        both, and its trial sensitivity per kilogram, and the cells beneath are untouched; the cells wholly taken leave
        by `take_top`. Returns (mass [kg], enthalpy [J], solid [kg], momentum [kg m/s] (n, 3), by_phase): solid the
        mass of solid in what left, by the cells' liquid fractions; by_phase None, or with `phase` ((n, K) bool) the
        same four split by the cells they came from, each with a leading axis of two, [0] the cells in `phase` and [1]
        the others. `momentum_out` is left holding the momentum."""
        dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
        if not np.all(np.isfinite(dm) & (dm >= 0.0)):
            raise ValueError("a release from a skin's top must be a finite, non-negative mass")
        valid, T = self.valid, self.T()
        run = np.cumprod(valid & np.asarray(admit, dtype=bool), axis=1).astype(bool)
        above = np.concatenate([np.zeros((self.n, 1)), np.cumsum(self.m, axis=1)[:, :-1]], axis=1)
        taken = np.where(run, np.clip(dm[:, None] - above, 0.0, self.m), 0.0)
        frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
        f_s = np.zeros(self.m.shape)
        f_s[valid] = 1.0 - self.props.liquid_fraction(T[valid])

        def parts(on):
            t = taken * on
            return t.sum(axis=1), (self.H * frac * on).sum(axis=1), (t * f_s).sum(axis=1), (t[..., None] * self.U).sum(axis=1)

        mass, H, solid, momentum = parts(1.0)
        by_phase = None
        if phase is not None:
            ph = np.asarray(phase, dtype=bool)
            by_phase = tuple(np.stack(pair) for pair in zip(parts(ph), parts(~ph)))
        whole = np.cumprod(run & (taken >= self.m), axis=1).sum(axis=1)       # the cells wholly taken, from the top
        part = (taken > 0.0) & (np.arange(self.m.shape[1])[None, :] >= whole[:, None])
        self.m, self.H = np.where(part, self.m - taken, self.m), np.where(part, self.H * (1.0 - frac), self.H)
        if self.dH is not None:
            self.dH = np.where(part, self.dH * (1.0 - frac), self.dH)
        rows = np.flatnonzero(whole > 0)
        if len(rows):
            self.take_top(rows, whole[rows])
        self.momentum_out = momentum
        return mass, H, solid, momentum, by_phase

    # -- columns appear, combine and leave with their markers (Task 5) -----------------------------------------------
    def _no_trial(self, what):
        if self.dH is not None:
            raise RuntimeError("columns cannot {} while a trial pass is open".format(what))

    def divide(self, rows):
        """Each column in `rows` (distinct) becomes two with the same profile: every cell's mass and enthalpy halved
        (exactly, a power of two), the velocities kept, the top area halved, the curvatures kept, the measured
        cross-sections halved (ruling T3-a). The new columns are appended; returns their indices, in the order of
        `rows`."""
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        if len(np.unique(rows)) != len(rows):
            raise ValueError("divide needs distinct columns")
        self._no_trial("divide")
        self.m[rows] *= 0.5
        self.H[rows] *= 0.5
        self.area[rows] *= 0.5
        if self.profile is not None:
            self.profile[rows] *= 0.5
            self.profile = np.concatenate([self.profile, self.profile[rows]])
        new = self.n + np.arange(len(rows))
        self.m = np.concatenate([self.m, self.m[rows]])
        self.H = np.concatenate([self.H, self.H[rows]])
        self.U = np.concatenate([self.U, self.U[rows]])
        self.count = np.concatenate([self.count, self.count[rows]])
        self.area = np.concatenate([self.area, self.area[rows]])
        self.kappa = np.concatenate([self.kappa, self.kappa[rows]])
        self.momentum_out = np.concatenate([self.momentum_out, np.zeros((len(rows), 3))])
        return new

    def cut(self, rows, fractions):
        """Fractions of columns (stage 3 Task 9b, the step end's re-partition): for each k a new column holding
        fractions[k] > 0 of column rows[k] with its profile -- every cell's mass and enthalpy times the fraction, its
        velocity kept -- and that fraction of its top area and of its measured cross-sections (ruling T3-a), its
        curvatures kept. A row may be cut several times; each source keeps the rest, one less the sum of its fractions,
        of every cell, of its area and of its cross-sections. A source whose fractions sum to one (to within CUT_SLACK)
        is emptied and its fractions are scaled to sum to one. Mass, enthalpy and momentum add to round-off. The new
        columns are appended; returns their indices, in the order of `rows`."""
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        f = np.asarray(fractions, dtype=float).reshape(-1)
        if len(f) != len(rows) or not np.all(np.isfinite(f) & (f > 0.0)):
            raise ValueError("cut needs one positive, finite fraction per row")
        self._no_trial("cut")
        total = np.bincount(rows, f, minlength=self.n)
        if np.any(total > 1.0 + CUT_SLACK):
            raise ValueError("the fractions cut from a column must sum to at most one")
        whole = (total > 0.0) & (total >= 1.0 - CUT_SLACK)
        f = f / np.where(whole, total, 1.0)[rows]
        rest = np.where(whole, 0.0, 1.0 - total)
        new = self.n + np.arange(len(rows))
        m, H, U = self.m[rows] * f[:, None], self.H[rows] * f[:, None], self.U[rows]
        count, area, kappa = self.count[rows], self.area[rows] * f, self.kappa[rows]
        profiled = self.profile is not None
        if profiled:
            profile = self.profile[rows] * f[:, None]
        src = np.flatnonzero(total > 0.0)
        self.m[src] *= rest[src, None]
        self.H[src] *= rest[src, None]
        self.area[src] *= rest[src]
        if profiled:
            self.profile[src] *= rest[src, None]
            self.profile = np.concatenate([self.profile, profile])
        empty = src[whole[src]]
        self.U[empty], self.count[empty] = 0.0, 0
        self.m = np.concatenate([self.m, m])
        self.H = np.concatenate([self.H, H])
        self.U = np.concatenate([self.U, U])
        self.count = np.concatenate([self.count, count])
        self.area = np.concatenate([self.area, area])
        self.kappa = np.concatenate([self.kappa, kappa])
        self.momentum_out = np.concatenate([self.momentum_out, np.zeros((len(rows), 3))])
        return new

    def combine(self, keep, gone, by="mass"):
        """The hand-over rule (design 9.2): column gone[k] is mixed into column keep[k] layer by layer at the same
        relative depth, the depth measured by the mass above as a fraction of the column's mass. The merged column's
        cell boundaries are the union of both columns' boundaries in that fraction; each merged cell takes from each
        column the share of the cell spanning it, mass, enthalpy and momentum in proportion (that cell's specific
        enthalpy and velocity), so the books add to round-off. Cells left below 1e-9 of the column's largest by nearly
        coinciding boundaries join a neighbour (`_merge`), exactly. The top areas add, the curvatures are averaged by
        area, the measured cross-sections add (ruling T3-a), and column gone[k] is left empty with no area (remove it
        with `remove`). Pairs must be distinct and disjoint; the caller regrids.

        by="depth" (plan 03b Task 3, design 21.5: metal stays where it is, only its owner changes) is a conservative
        remap on depth instead: each of the keeper's cells takes the gone column's metal lying at its own depths (the
        share of each gone cell's volume between the keeper cell's top and bottom, from the gone column's own
        cross-section), and the gone metal deeper than the keeper's bottom joins it as new cells below, the gone cells'
        parts in their order. The keeper's cell boundaries are kept: with both columns exactly full and their profiles
        summed, every cell spans the depths it spanned. Shares below END_SLACK of a cell are dropped and the rest scaled
        to sum to one, so mass, enthalpy and momentum add to round-off; a gone cell whose every share is so dropped (a
        sliver at the round-off of the volumes the overlaps are reckoned from: plan 03b Task 6 found one of 2e-27 m3 in a
        column of 7e-10 m3, whose shares were 0/0) joins whole the keeper's cell at its depth, or the metal below the
        keeper's bottom."""
        keep = np.asarray(keep, dtype=np.int64).reshape(-1)
        gone = np.asarray(gone, dtype=np.int64).reshape(-1)
        both = np.concatenate([keep, gone])
        if len(keep) != len(gone) or len(np.unique(both)) != len(both):
            raise ValueError("combine needs distinct, disjoint pairs of columns")
        if by not in ("mass", "depth"):
            raise ValueError("combine is by 'mass' or by 'depth', not {!r}".format(by))
        self._no_trial("combine")
        merged = self._by_depth(keep, gone) if by == "depth" else self._by_mass(keep, gone)
        self._width(max([len(t[0]) for t in merged] + [1]))
        for a, b, (m, H, P) in zip(keep, gone, merged):
            k = len(m)
            self.m[a], self.H[a], self.U[a] = 0.0, 0.0, 0.0
            self.m[a, :k], self.H[a, :k] = m, H
            self.U[a, :k] = np.divide(P, m[:, None], out=np.zeros_like(P), where=m[:, None] > 0.0)
            self.count[a] = k
            self.m[b], self.H[b], self.U[b], self.count[b] = 0.0, 0.0, 0.0, 0
            A = self.area[a] + self.area[b]
            if A > 0.0:
                self.kappa[a] = (self.area[a] * self.kappa[a] + self.area[b] * self.kappa[b]) / A
            self.area[a], self.area[b] = A, 0.0
            if self.profile is not None:
                self.profile[a] += self.profile[b]
                self.profile[b] = 0.0
            self.momentum_out[a] += self.momentum_out[b]
            self.momentum_out[b] = 0.0
        tiny = np.zeros(self.m.shape, bool)
        sub = self.m[keep]
        tiny[keep] = (self.valid[keep] & (sub < 1e-9 * sub.max(axis=1, keepdims=True))
                      & (self.count[keep] >= 2)[:, None])
        if tiny.any():
            self._merge(np.flatnonzero(tiny.any(axis=1)), tiny)
        self._trim()

    def _by_mass(self, keep, gone):
        """`combine`'s merged cells by relative depth in mass, (m, H, P) per pair, before any column changes."""
        merged = []
        for a, b in zip(keep, gone):
            cols = [(self.m[c, :self.count[c]], self.H[c, :self.count[c]], self.U[c, :self.count[c]])
                    for c in (a, b) if self.count[c] > 0]
            F = [np.concatenate([[0.0], np.cumsum(m) / m.sum()]) for m, _, _ in cols]
            for f in F:
                f[-1] = 1.0
            s = np.unique(np.concatenate(F)) if F else np.zeros(1)
            lo, hi = s[:-1], s[1:]
            mid = 0.5 * (lo + hi)
            k = len(lo)
            m, H, P = np.zeros(k), np.zeros(k), np.zeros((k, 3))
            for (mc, Hc, Uc), f in zip(cols, F):
                c = np.clip(np.searchsorted(f, mid, side="right") - 1, 0, len(mc) - 1)
                w = (hi - lo) / (f[c + 1] - f[c])
                m += mc[c] * w
                H += Hc[c] * w
                P += (mc[c] * w)[:, None] * Uc[c]
            merged.append((m, H, P))
        return merged

    def _by_depth(self, keep, gone):
        """`combine(..., by="depth")`'s merged cells, (m, H, P) per pair, before any column changes. A gone column whose
        own cross-section holds no volume (an all-zero profile row, whatever its top area, or no top area without a
        profile: a marker left with no area, ruling 6) has no depths to remap, and its pair takes the by-mass rule
        (`_by_mass`), which keeps its top metal at the keeper's top instead of putting it below the keeper's bottom
        (Task 3 fix round 1). Only the keepers' depths are taken, so such a column's own depths are never asked for."""
        T = self.T()
        dz = self._profile_dz(T, keep) if self.profile is not None else self.dz(T)[keep]     # the keepers' alone
        no_volume = self.column_volume()[gone] <= 0.0
        valid = self.valid
        vol = np.zeros(self.m.shape)
        vol[valid] = self.m[valid] / self.density(T[valid])
        merged = []
        for i, (a, b) in enumerate(zip(keep, gone)):
            ka, kb = self.count[a], self.count[b]
            m, H = self.m[a, :ka].copy(), self.H[a, :ka].copy()
            P = m[:, None] * self.U[a, :ka]
            if kb == 0:
                merged.append((m, H, P))
                continue
            if no_volume[i]:
                merged.append(self._by_mass([a], [b])[0])
                continue
            za = np.concatenate([[0.0], np.cumsum(dz[i, :ka])])           # the keeper's boundaries
            Va = self._volume_at(za[None, :], [b])[0]                      # as volumes from the gone column's top
            vb = np.concatenate([[0.0], np.cumsum(vol[b, :kb])])
            lo = np.maximum(vb[:-1, None], Va[None, :-1])
            hi = np.minimum(vb[1:, None], Va[None, 1:])
            share = np.column_stack([np.maximum(hi - lo, 0.0), np.maximum(vb[1:] - np.maximum(vb[:-1], Va[-1]), 0.0)])
            share /= vol[b, :kb, None]
            share[share < END_SLACK] = 0.0
            sliver = np.flatnonzero(share.sum(axis=1) <= 0.0)             # wholly at the volumes' round-off (docstring)
            if len(sliver):
                at = np.clip(np.searchsorted(Va, 0.5 * (vb[sliver] + vb[sliver + 1]), side="right") - 1, 0, ka)
                share[sliver, at] = 1.0
            share /= share.sum(axis=1, keepdims=True)                      # gone cell c's share in each keeper cell,
            mg, Hg, Ug = self.m[b, :kb, None] * share, self.H[b, :kb, None] * share, self.U[b, :kb]    # and beyond
            m += mg[:, :ka].sum(axis=0)
            H += Hg[:, :ka].sum(axis=0)
            P += np.einsum("cj,ck->jk", mg[:, :ka], Ug)
            deeper = np.flatnonzero(mg[:, ka] > 0.0)                       # the gone metal below the keeper's bottom
            merged.append((np.concatenate([m, mg[deeper, ka]]), np.concatenate([H, Hg[deeper, ka]]),
                           np.concatenate([P, mg[deeper, ka, None] * Ug[deeper]])))
        return merged

    def remove(self, rows):
        """The columns `rows` (distinct) leave; the others keep their order. Returns per removed column, in the order
        of `rows`, its (mass [kg], enthalpy [J], momentum [kg m/s] (., 3)) from `totals`."""
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        if len(np.unique(rows)) != len(rows):
            raise ValueError("remove needs distinct columns")
        self._no_trial("remove")
        M, E, P = self.totals()
        out = M[rows], E[rows], P[rows]
        keep = np.ones(self.n, bool)
        keep[rows] = False
        self.m, self.H, self.U, self.count = self.m[keep], self.H[keep], self.U[keep], self.count[keep]
        self.area, self.kappa, self.momentum_out = self.area[keep], self.kappa[keep], self.momentum_out[keep]
        if self.profile is not None:
            self.profile = self.profile[keep]
        self._trim()
        return out

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
        as dH / (m dh/dT), the slope of the table's h."""
        if self.dH is None:
            return None
        valid, T = self.valid, self.T()
        S = np.zeros(self.m.shape)
        S[valid] = self.dH[valid] / (self.m[valid] * self.props.dh_dT(T[valid]))
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
    def conduct(self, dt, q_top, T_amb, emissivity, T_b, trial=False, side=None):
        """One implicit (backward-Euler) heat sub-step dt [s] of every column through its depth: the convective flux
        q_top [W/m2] over the top area A(0) in at the top cell, radiation emissivity * SIGMA_SB (T0^4 - T_amb^4) over
        the same area out of it (the shared Stefan-Boltzmann constant), conduction between cells through harmonic
        conductances A(z_i) / (dz_i / 2k_i + dz_j / 2k_j) with the cross-section at the interface's depth z_i, and the
        base held at T_b [K] through the bottom half-cell, A(z_bottom) k / (dz / 2) (T_b None: every base adiabatic;
        a column whose T_b is not finite has an insulated base alone, plan 03b Task 5: a column with no base).
        Under a measured profile A(z_i) is the cells' cross-section (`cell_section`) and A(z_bottom) the base's
        (`base_area`: none for a column that pinches, so its bottom is insulated); the top fluxes stay on A0 = `area`.
        The cells' thicknesses (and so the interface depths) are those at the sub-step's start. Newton in the
        temperature with a backtracking line search per column on the residual's norm, every property from the table
        (`PropertyTable.conduction`), the Jacobian with the slopes of its h and k across each cell's grid interval, so
        it is the exact tangent of the residual solved. An empty column passes its heat straight to the base.

        With `side` (a SideFaces, plan 03b Task 4) the columns also conduct sideways across their shared boundaries,
        in the same backward-Euler step (module docstring; planning ruling 5): each pair of cells overlapping in depth
        across a boundary through Lambda k_h / |x_j - x_i| at the sub-step's start's depths, the conductivities' harmonic
        mean at the iterate's temperatures. The Newton iteration is then one over every column's cells together, its
        Jacobian (the columns' tridiagonal blocks, the sideways couplings and their dk/dT terms) solved by GMRES
        preconditioned by the columns' tridiagonal solves (`_Sideways.solve`), the line search still per column; the
        heat each column gains sideways is the book's E_side, and it enters E_base's residual. Without `side` every
        column is solved apart by the Thomas algorithm, as before.

        With trial=True the sub-step belongs to a trial pass (opened by the first such call, closed by `base_law` or
        `set_state`): each cell's dH/dT_b is carried from the last sub-step and updated by the converged step's own
        linearisation, J S = dH_old / dt + G_b e_bottom, so the book's dE_base_dTb is the exact derivative of E_base
        for the linearised step (and of the pass, through the carried dH). With `side` the right-hand side holds every
        column's G_b at once, so S and dE_base_dTb (its sideways part included) are the derivatives for every column's
        T_b moving together (an exact derivative with respect to a column's own T_b alone would need one solve per
        column). Returns a HeatBook."""
        n = self.n
        valid, rows = self.valid, np.arange(n)
        has = self.count > 0
        last = np.maximum(self.count - 1, 0)
        props = self.props
        q = np.broadcast_to(np.asarray(q_top, dtype=float), (n,)) * self.area
        es = np.broadcast_to(np.asarray(emissivity, dtype=float), (n,)) * SIGMA_SB * self.area
        T_amb = np.broadcast_to(np.asarray(T_amb, dtype=float), (n,))
        base = T_b is not None
        T_b = np.broadcast_to(np.asarray(T_b if base else 0.0, dtype=float), (n,))
        insulated = ~np.isfinite(T_b)                     # a column with no base: no conduction through its bottom
        T_b = np.where(insulated, 0.0, T_b)
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
        A_i = np.where(link, self.cell_section(z_bot[:, :-1]), 1.0)           # at each interface's depth
        if not base:
            A_b = np.zeros(n)
        elif self.profile is not None:
            A_b = self._base_section(T_old)
        else:
            A_b = self.cross_section(z_bot[rows, last][:, None])[:, 0]
        bottom = has & base & ~insulated
        m, E_old = self.m, self.H
        sideways = None
        if side is not None:
            z = np.concatenate([np.zeros((n, 1)), z_bot], axis=1)            # each column's cell boundaries
            sideways = _Sideways(side, z, self.count)
            solve = sideways.solve
        else:
            solve = lambda lower, diag, upper, jac, rhs: _tridiagonal(lower, diag, upper, rhs)

        def system(T):
            """The residual R, the Jacobian's (lower, diag, upper), the base conductance, cp and, with sideways
            conduction, the pairs' fluxes and their derivatives (F, dFa, dFb) at T."""
            h, cp, k, kp = props.conduction(T)
            r = 0.5 * dz / k
            G = np.where(link, A_i / (r[:, :-1] + r[:, 1:]), 0.0)
            dT = T[:, :-1] - T[:, 1:]
            F = G * dT                                                          # down through each interface
            dG = np.divide(G * G, A_i, out=np.zeros_like(G), where=A_i > 0.0) * 0.5   # dG/dk_j = dG * dz_j / k_j^2
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
            pairs = None
            if sideways is not None:
                pairs = sideways.fluxes(T, k, kp)
                R += sideways.residual(pairs[0])
            return (np.where(valid, R, 0.0), np.where(valid, lower, 0.0), np.where(valid, D, 1.0),
                    np.where(valid, upper, 0.0), Gb, cp, pairs)

        T = T_old
        R, lo, D, up, Gb, cp, pairs = system(T)
        moved = np.zeros(n)
        for _ in range(NEWTON_MAX):
            dT = solve(lo, D, up, None if pairs is None else pairs[1:], -R)
            r0 = np.sqrt((R * R).sum(axis=1))
            lam = np.ones(n)
            for _ in range(LINE_SEARCH_MAX):
                T_try = T + lam[:, None] * dT
                trying = system(T_try)
                R_try = trying[0]
                moved = np.abs(T_try - T).max(axis=1)
                worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
                if not worse.any():
                    break
                lam = np.where(worse, 0.5 * lam, lam)
            T, (R, lo, D, up, Gb, cp, pairs) = T_try, trying
            if moved.max() <= NEWTON_TOL:
                break
        else:
            raise RuntimeError("a skin heat sub-step did not converge in {} iterations (last change {:.2e} K)".format(
                NEWTON_MAX, moved.max()))
        E_new = np.where(valid, m * props.enthalpy(T), 0.0)
        E_top = q * dt
        E_rad = np.where(has, es * (T[:, 0] ** 4 - T_amb ** 4) * dt, 0.0)
        E_side = np.zeros(n) if pairs is None else sideways.gained(pairs[0]) * dt
        E_base = E_top - E_rad + E_side - (E_new.sum(axis=1) - E_old.sum(axis=1))
        book = HeatBook(E_top, E_rad, E_base, E_side)
        if trial:
            rhs = self.dH / dt
            rhs[rows, last] += Gb                                  # -dR/dT_b: the base's pull on the bottom cell
            S = solve(lo, D, up, None if pairs is None else pairs[1:],
                      np.where(valid, rhs, 0.0))                   # the converged step's own Jacobian
            dH_new = np.where(valid, m * cp * S, 0.0)
            dE = np.where(has, -4.0 * es * T[:, 0] ** 3 * S[:, 0] * dt, 0.0) - (dH_new - self.dH).sum(axis=1)
            if pairs is not None:                                  # the sideways heat's own change
                _, dFa, dFb = pairs
                dE = dE + sideways.gained(dFa * S.ravel()[sideways.a] + dFb * S.ravel()[sideways.b]) * dt
            book.dE_base_dTb = dE
            self.dH = dH_new
            self._trial_E += E_base
            self._trial_dE += dE
            self._trial_Tb = np.array(T_b, dtype=float)
        self.H = E_new
        return book

    # -- removal and the base ------------------------------------------------------------------------------------
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
        """The cross-section [m2] (n,) at each column's bottom, A(z_bottom); A0 for an empty column. Under a measured
        profile, from the cells' volume (`_base_section`): the profile's base for a full column, none where it
        pinches."""
        if self.profile is not None:
            return self._base_section(self.T())
        return self.cross_section(self.thickness()[:, None])[:, 0]

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

    # -- the layered flow through the depth (Task 4) ----------------------------------------------------------------
    def flowing(self, T=None):
        """The cells that flow (n, K): above the material's T_rigid, where the solid fraction is below the shared
        rheology's RIGID_SOLID_FRACTION (a half). `T` (n, K) the cells' temperatures if already known."""
        T = self.T() if T is None else T
        return self.valid & (T > self.material.T_rigid)

    def liquid(self, T=None):
        """The liquid cells (n, K): at or above the material's liquidus."""
        T = self.T() if T is None else T
        return self.valid & (T >= self.material.T_liquidus)

    def viscosity(self, gamma, T=None):
        """Every cell's viscosity [Pa s] (n, K) sheared at gamma [1/s] (per cell, or one for all): the shared
        `rheology.viscosity` at the cell's temperature and solid fraction (1 - the table's f_l) with the material's
        liquid viscosity and its solid fraction at the rheology's BRIDGE_LOW; infinite in a cell that does not flow
        (rigid, padded, or of a material without a liquid)."""
        T = self.T() if T is None else T
        return self._viscosity(T, self.flowing(T), gamma)

    def _viscosity(self, T, flowing, gamma):
        """`viscosity` for cells at temperatures T whose flowing cells are `flowing` (arrays of one shape)."""
        eta = np.full(T.shape, np.inf)
        liquid = getattr(self.material, "liquid", None)
        if liquid is None or not flowing.any():
            return eta
        gamma = np.broadcast_to(np.asarray(gamma, dtype=float), T.shape)
        f_s = np.minimum(1.0 - self.props.liquid_fraction(T[flowing]), FLOWING_SOLID_MAX)
        f_s_bridge = 1.0 - float(self.props.liquid_fraction(rheology.BRIDGE_LOW))
        eta[flowing] = rheology.viscosity(T[flowing], f_s, gamma[flowing], liquid.mu, f_s_bridge)
        return eta

    def momentum_step(self, dt, tau, push, response=False):
        """One implicit (backward-Euler) step dt [s] of every column's velocities through its depth, each component
        apart (the skyscrapers' momentum step, matched by amendment):

            m_i (u_i' - u_i) / dt = A(z_i) D_i (u'_(i+1) - u'_i) + A(z_(i-1)) D_(i-1) (u'_(i-1) - u'_i)
                                    + push V_i + tau A_top [i the top cell],

        A(z) the column's cross-section at each interface's depth (the top area for a flat column, or the measured
        profile's `cell_section`, the thicknesses at the step's start), D = 1 / (dz_i / 2 eta_i + dz_j / 2 eta_j) the drag per
        unit area between two cells' centres (the two half-cells in series), V_i = m_i / rho_i the cell's volume, `push`
        [N/m3, (n, 3) or (3,)] the force per unit volume, the same on every layer of a column (built by Task 6), and
        `tau` [Pa, (n, 3) or (3,)] the gas shear on the top cell over the top area A_top = A(0) = `area`, a top flux
        even under a profile (ruling T1-g); under a profile the floor's A(z_bottom) is `base_area`'s, none where the
        column pinches. The anchored cells (`anchored`) are held at rest: a moving cell
        drags on them, and a flowing bottom cell on the rigid volume beneath the base, over its own half-thickness
        alone, A(z_bottom) eta / (dz / 2) (no slip at the rigid front). A crust (`crust`) floats: its contiguous cells are
        one unknown with their summed mass, momentum and forces, dragged at its faces through a half-cell of
        CRUST_VISCOSITY. Each flowing cell's viscosity (`viscosity`) is at its own shear rate, the mean over its two
        faces of the stress there over its viscosity, found by fixed point from the step's starting velocities: at
        most SHEAR_PASSES passes, ending when no flowing cell's viscosity moves by more than SHEAR_TOL relatively.

        With `response` the book also carries each cell's velocity per unit push and the wall's share of it (MomentumBook).
        Updates U and returns a MomentumBook: the impulses of the push and the shear on the moving cells and the
        momentum given to the wall (the anchored metal and the volume beneath), so that the columns' momentum change is
        push + shear - wall to round-off. Only the columns where a cell moves are solved; in the others every cell is
        anchored, ends at rest, and gives the wall the momentum it held."""
        n = self.n
        tau = np.broadcast_to(np.asarray(tau, dtype=float), (n, 3))
        push = np.broadcast_to(np.asarray(push, dtype=float), (n, 3))
        T = self.T()
        valid = self.valid
        anchored = self._anchored(T)
        moving = valid & ~anchored
        wall = ((self.m * anchored)[..., None] * self.U).sum(axis=1)       # what the anchored cells held, now at rest
        book = MomentumBook(np.zeros((n, 3)), np.zeros((n, 3)), wall, 0)
        if response:
            book.response, book.wall_response = np.zeros(self.m.shape), np.zeros(n)
        U = np.zeros_like(self.U)
        live = np.flatnonzero(moving.any(axis=1))
        if len(live):
            U[live], push_i, shear_i, wall_i, book.passes, resp = self._momentum_columns(
                live, T, anchored, dt, tau[live], push[live], response)
            book.push[live], book.shear[live] = push_i, shear_i
            book.wall[live] += wall_i
            if response:
                book.response[live], book.wall_response[live] = resp
        self.U = U
        return book

    def _momentum_columns(self, r, T_all, anchored, dt, tau, push, response=False):
        """`momentum_step` on the columns r (each with a moving cell): (U (len r, K, 3), push impulse, shear impulse,
        the wall's share of the moving cells' momentum (both (len r, 3)), the fixed point's passes, and with `response`
        (the cells' velocity per unit push (len r, K), its drag on the wall (len r,)), else None)."""
        nr, K = len(r), self.m.shape[1]
        rows = np.arange(nr)
        valid, T = self.valid[r], T_all[r]
        m, count, U0 = self.m[r], self.count[r], self.U[r]
        last = count - 1
        flowing = valid & (T > self.material.T_rigid)                   # as `flowing`
        anchored = anchored[r]
        moving = valid & ~anchored                                       # the flowing cells and the floating crusts
        crust = moving & ~flowing
        dz = np.where(valid, np.maximum(self.dz(T_all)[r], DZ_FLOOR), 0.0)       # the thicknesses at the start
        z_bot = np.cumsum(dz, axis=1)
        A_top = self.area[r]                                             # the gas shear's area, a top flux (T1-g)
        if self.profile is None:                                         # flat: the top area at every depth
            A_z = np.broadcast_to(A_top[:, None], z_bot.shape)
            A_i, A_b = A_z[:, :-1], A_z[rows, last]                      # at each interface's depth, at the base
        else:                                                            # the measured profile's (plan 03b Task 3)
            A_i, A_b = self._section(z_bot[:, :-1], r, extend=True), self._base_section(T_all, r)
        V = np.where(valid, m / np.where(valid, self.density(T), 1.0), 0.0)
        force = np.where(moving[..., None], push[:, None, :] * V[..., None], 0.0)
        top_moves = moving[:, 0]
        force[:, 0] += np.where(top_moves[:, None], tau * A_top[:, None], 0.0)
        bottom_flows = flowing[rows, last]

        # the unknowns: every cell, but each crust's contiguous cells one unknown together
        same = np.zeros((nr, K), dtype=bool)
        same[:, 1:] = crust[:, 1:] & crust[:, :-1]
        unit = np.cumsum(~same, axis=1) - 1
        Ku = int(unit.max()) + 1
        flat = (rows[:, None] * Ku + unit).ravel()

        def lump(x):
            return np.bincount(flat, x.ravel(), nr * Ku).reshape(nr, Ku)

        moving_u = lump(moving.astype(float)) > 0.0
        last_u = unit[rows, last]
        between = ~same[:, 1:]                                            # interfaces joining two different unknowns
        ir, ik = np.nonzero(between)
        mass_u = lump(np.where(moving, m, 0.0))
        rhs_u = np.stack([lump(np.where(moving, m * U0[..., c] / dt + force[..., c], 0.0)) for c in range(3)])

        def viscosities(gamma):
            return np.where(crust, CRUST_VISCOSITY, self._viscosity(T, flowing, gamma))

        def drags(eta):
            """The drag per unit area at each interface (nr, K-1) and at the base (nr,) [Pa s/m]."""
            resist = 0.5 * dz[:, :-1] / eta[:, :-1] + 0.5 * dz[:, 1:] / eta[:, 1:]
            link = valid[:, 1:] & (np.isfinite(eta[:, :-1]) | np.isfinite(eta[:, 1:]))
            D = np.where(link, 1.0 / np.where(link, resist, 1.0), 0.0)
            eta_b = eta[rows, last]
            Db = np.where(bottom_flows, np.where(np.isfinite(eta_b), eta_b, 0.0) / (0.5 * dz[rows, last]), 0.0)
            return D, Db

        def solve(D, Db, rhs_u=rhs_u):
            """The cells' velocities (nr, K, c) for the right-hand sides rhs_u (c, nr, Ku), one per component."""
            nc = rhs_u.shape[0]
            C = A_i * D                                                   # N s/m: the drag through each interface
            Cu = np.zeros((nr, max(Ku - 1, 0)))
            Cu[ir, unit[ir, ik]] = C[ir, ik]                              # unit a's face with unit a + 1
            diag = mass_u / dt
            diag[:, :-1] += Cu
            diag[:, 1:] += Cu
            diag[rows, last_u] += A_b * Db
            lower, upper = np.zeros((nr, Ku)), np.zeros((nr, Ku))
            lower[:, 1:] = np.where(moving_u[:, :-1], -Cu, 0.0)
            upper[:, :-1] = np.where(moving_u[:, 1:], -Cu, 0.0)
            held = ~moving_u
            diag, lower, upper = np.where(held, 1.0, diag), np.where(held, 0.0, lower), np.where(held, 0.0, upper)
            def three(x):                                                 # the components solved as rows
                return np.concatenate([x] * nc)
            rhs = np.where(three(held), 0.0, rhs_u.reshape(nc * nr, Ku))
            u = _tridiagonal(three(lower), three(diag), three(upper), rhs)
            u = np.moveaxis(u.reshape(nc, nr, Ku), 0, -1)
            return np.where(moving[..., None], u[rows[:, None], unit], 0.0)

        def rates(U, D, Db, eta):
            """Every flowing cell's shear rate [1/s]: the mean over its two faces of the stress there over its
            viscosity."""
            stress = np.zeros((nr, K + 1))                                # at every face, the top face first
            stress[:, 1:K] = D * np.linalg.norm(U[:, :-1] - U[:, 1:], axis=2)
            stress[:, 0] = np.where(top_moves, np.linalg.norm(tau, axis=1), 0.0)
            stress[rows, last + 1] = Db * np.linalg.norm(U[rows, last], axis=1)
            face = 0.5 * (stress[:, :-1] + stress[:, 1:])
            return np.where(flowing, face / np.where(flowing, eta, 1.0), 0.0)

        eta = viscosities(0.0)
        eta = viscosities(rates(U0 * moving[..., None], *drags(eta), eta))   # at the step's starting rates
        for passes in range(1, SHEAR_PASSES + 1):
            D, Db = drags(eta)
            U = solve(D, Db)
            eta_new = viscosities(rates(U, D, Db, eta))
            change = np.abs(eta_new[flowing] - eta[flowing]) / eta[flowing]
            if change.size == 0 or change.max() <= SHEAR_TOL:
                break
            eta = eta_new
        on_anchored = np.zeros((nr, K), dtype=bool)                       # a moving cell on the anchored metal beneath
        on_anchored[:, :-1] = moving[:, :-1] & anchored[:, 1:]
        drag = np.zeros((nr, K))
        drag[:, :-1] = np.where(on_anchored[:, :-1], A_i * D, 0.0)
        drag[rows, last] += np.where(bottom_flows, A_b * Db, 0.0)
        to_wall = dt * (drag[..., None] * U).sum(axis=1)
        push_impulse = dt * np.where(moving[..., None], push[:, None, :] * V[..., None], 0.0).sum(axis=1)
        shear_impulse = dt * np.where(top_moves[:, None], tau * A_top[:, None], 0.0)
        resp = None
        if response:                                                       # the last pass's system, a unit push alone
            rr = solve(D, Db, lump(np.where(moving, V, 0.0))[None])[..., 0]
            resp = rr, (drag * rr).sum(axis=1)
        return U, push_impulse, shear_impulse, to_wall, passes, resp
