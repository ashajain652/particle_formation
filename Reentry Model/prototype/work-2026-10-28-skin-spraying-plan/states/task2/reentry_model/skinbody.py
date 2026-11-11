"""The body of the layered-skin design (docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md sections 4,
6.6 and 7; stage 2, sub-plan 20: heat only, liquid removed as it forms).

A moving 3D mesh holds the rigid metal; every corner of its surface carries a column with a layered skin holding the
metal outward of the rigid front (`skin.Skins`). At time zero each column takes its minimum mass from the 3D model at
the initial temperature: the corners move inward by it and the skins receive exactly what they swept. Per macro step,
after the trajectory and the heating:

1. the trial pass gives every column its base law, E_base / dt = a + b T_b, with the base held at the 3D model's
   start-of-step temperature (`skin.advance`);
2. the 3D step takes that law at the corner nodes, implicitly, with no convective load and no radiation on its
   boundary (the skins' tops radiate), and last step's deferred mismatch as a nodal load;
3. the real pass, base held at the new 3D temperature, books the heat, removes the liquid and draws to the minimum;
   the difference between what it sent through the base and what the 3D step received is deferred to the node. By
   default ("damped", Asha, 2026-10-12) the step repeats, each pass re-linearised at a temperature moved only part of
   the way to the last 3D answer, until the total mismatch is within INTERFACE_TOLERANCE of the throughput; the other
   rules repeat at most once (REPEAT_RULES). Every draw takes the next layer down the 3D model's linear profile below
   the base, at that layer's enthalpy (`skin.profile_enthalpy`), the profile's slope taken at the step's start;
4. the base moves by the metal crossing it (spec 6.6): drawn (in the sub-steps) and absorbed metal move the corner in,
   handed-down metal moves it out; with the volume-neutral smoothing this is one displacement per corner, which
   `MovingMesh.move_with_policy` applies (rebuilding when needed), and each column is then reconciled with what its
   corner actually swept: its bottom gains the swept mass less what it already exchanged, and its enthalpy changes by
   the swept enthalpy less what it already exchanged, so mass and energy are exact. A sweep that takes metal out of
   the 3D model is booked along the same profile as the draws, one that adds metal to it at the enthalpy of what the
   skins handed down there; and only rigid metal no warmer than the base node by more than the fault tolerance is
   handed down, warmer metal waiting in the skin until it has cooled (sub-plan 21 Task 9, Asha 2026-10-26, options A
   and C). Booked at the node's enthalpy, as sub-plan 20 built it, the handed metal's excess over the node stayed in
   the skin's new bottom cell, often a sliver of a cell that the hand-down took in part: where runoff froze onto cold
   metal on Task 8's pushed plate the deposit went down 35-67 K/mm warmer than its node, a sliver of 1.4 % of a cell
   reached 3 709 K and bottom cells 1 253 K, and the deposit's growth, 1.1-3.2 move allowances a 0.0625 s part, had
   every melting part rejected for its moves down to 0.0104 s. A part hands down at most what moves its base by the
   move allowance, a cell going in part (Asha, 2026-10-27, option 1 of Task 9's third hand-back): the limit releases a
   held cell whole once it is near the node's temperature, and the deep cells of the full pushed plate's thick skins
   were 0.2-0.8 mm of solid metal against its 0.556 mm allowance, a release that did not shrink with the part, so from
   1.56 s parts were split six times in vain, eight of the last fifteen were kept over the allowance and the 12.4 g
   held came down in 0.039 s. What the tolerance and the cap keep in the skins is recorded at every base move
   (`held_mass`; unbounded, Asha 2026-10-26), and the bottom cell, read after the reconciliation, is judged against the
   metal that crossed, for a hand-down the handed cells and the cell it uncovered (`jump_faults`).

A macro step is run as equal parts so that no base moves more than MOVE_MAX of the mean surface edge in one part
(Asha, 2026-10-12): a part whose draws would move a base further is run again as shorter parts, and the next macro
step starts with as many parts as this one's moves called for; a hand-down never splits a part, being capped at the
allowance (HAND_DOWN_SHARE of it, 2026-10-27). Absorption takes what is left of a part's
allowance after its draws and hand-downs, and the rest of the non-rigid metal waits for the next part (Asha,
2026-10-14: a cap of a quarter-edge per macro step stranded up to 1.3 mm of non-rigid metal in the 3D model of a plate
started at 880 K, whose exact non-rigid layer is 1.5-1.7 mm thick, and the plate melted 4-8 % slow). While any skin's top is non-rigid (at or above T_rigid at
the macro step's start), no part is longer than PART_MAX (Asha, 2026-10-12: long parts during melting were rejected for
their moves and needed many passes), and also in the macro step in which a top is about to become non-rigid: its rise
over the last macro step, carried on across this one, reaches T_rigid (Asha, 2026-10-16: with the skins' 0.8 mm minimum
the level-2 cold plate's top stood at 884.8 K at 0.5 s, so the step in which it melted ran as one part, and the plate
was 14.2 % ahead at 1.0 s; capped one step early, 3.2 %). Within a part the base heat, the draws and the move all take the part's own step;
the heating is interpolated in time between the last macro step's and this one's, at the part's middle, so a macro
step receives the trapezoid of the two (Asha, 2026-10-12).

The body's surface for the coupled run is the top surface, one 'patch' per corner (`columns.CornerSurface`), so the
heating, the history and the trajectory's loop are unchanged.

The flow (sub-plan 21 Task 6c; `flow`, off by default): every macro step the surface flow is evaluated on the corners
as the heating is, from the step's aero state and the tops' temperatures -- its wall-pressure gradient along the
meridian and its gas shear -- with the deceleration, which the metal feels along the flight direction; each part takes
them interpolated in time at its middle, as it takes the heating, and every sub-step of both passes runs the flow after
its conduction (`skin.flow_substeps`), as many flow sub-steps to each as `skin.flow_cycles` chooses from the last
heat sub-step's Courant number and growth (sub-plan 21 Task 6d; `flow_cycles` fixes them). With the flow on the liquid
stays and flows
(`remove` False in a run; spraying is stage 4's). The momentum books: the skins' momentum changes by the flow's impulses
less what leaves with removed or handed-down metal, and `momentum_residual` reports what they miss.

The droplet momentum sink (sub-plan 21 Task 7; `release`, None by default): the placeholder for stage 4's spraying, a
release rate per column (`skin.Flow.release`, what sub-plan 22 supplies from Girin's model) taken in every flow
sub-step from the top cells, liquid first, the exposed slurry too with `slurry_spray`. The released metal leaves the
body with its cells' enthalpy and momentum, booked as the removed liquid is, and the history reports the released mass
and fact 29's droplet-momentum ratio Pi over each macro step (`skin.droplet_momentum_ratio`).

Girin's spraying from the top cells (sub-plan 22 Task 1c; spec 8; `spray`, a `spray.SprayModel`, with the flow on): at
the start of every macro step Girin's model is evaluated once on the corners (`girin_plan`), from the surface flow the
flow's forcing evaluates there and the skins' depths then -- b and the regime test on the flowing depth (spec 6.7), the
viscous film's floor there (Asha, 2026-10-28), the Rayleigh-Taylor mode on the liquid pool and the wave on the liquid's
molten region, a liquid on slurry held from the shear modes off Girin's closure (Asha, 2026-10-29, option B), the droplet
capped by one crest of its wave on the depth the release may take (2026-10-24) -- and, with `slurry_spray`, the exposed
slurry's own mode by the cascade's properties (spec 8; Asha, 2026-10-29, option i). The plan's two rates per column are
the release of every flow sub-step of the macro step (`SprayPlan.release`, through `skin.cascade_mass`); what the
accepted passes released is the plan's record, per column and mode (`spray_plan`).

The droplet records (sub-plan 22 Task 2; spec 8 and 10.2): what the accepted passes released, by the cells it came from
-- the liquid's mode's from the fully liquid cells, the slurry's from the others -- is combined per column, per mode, per
macro step (`SprayPlan.book`), summing mass, enthalpy, solid and momentum over the sub-steps; a column's radius is the
macro step's, from its plan, so combining changes no radius. At the step's end the records become source-table rows
(`SprayPlan.source_block`, spray.SOURCE_COLUMNS: the element model's columns, the column's top corner and angle where
the plan was evaluated, b its flowing depth, the droplets leaving at the mean velocity of the metal they are made of --
the momentum the books lose with them -- and their solid fraction and enthalpy per kilogram, mass-weighted over the
cells released), kept one block per macro step in `source_rows`; the size histograms (`hist_n`, `hist_m`) accumulate
each flow sub-step's release at the plan's radii (`SprayPlan.histogram`, through `skin.Cascade.histogram`)."""
import copy
import dataclasses
import math
import time
import types

import numpy as np

from . import columns as columns_mod
from . import film as film_mod
from . import material as material_mod
from . import movemesh
from . import rheology
from . import skin as skin_mod
from . import spray as spray_mod
from .body import ThermalBody
from .thermal import SIGMA_SB

FLOW_DEFAULT = False              # the flow is off until stage 4's spraying gives the liquid a sink: with it on and nothing
                                  # sprayed, a flight would only pile the liquid up (Task 6c's ruling, 2026-10-22)

GRADIENT_FLOOR = 1.0              # K/m: the smallest gradient along the normal that is divided by. A corner above T_rigid
                                  # whose gradient is below it (a flat or reversed field) is given this one, so its depth
                                  # is huge and it absorbs up to its part's move allowance: it does not absorb nothing
REPEAT_RULES = ("column", "total", "never", "damped")   # when a step's interface mismatch repeats it: "damped"
                                              # repeats until converged (below); the others at most once -- "column" (the
                                              # spec's first rule) when any column's mismatch exceeds the tolerance,
                                              # "total" when the total does -- or never (sub-plan 20 Task 10)
DEFAULT_REPEAT_RULE = "damped"    # Asha, 2026-10-12: the only rule whose coupling converges at 0.5 s macro steps
# The damped repetition (2026-10-12). The skins' base heat has a cliff a little above the start-of-step base
# temperature once the skins near melting: a base held above the melting band for a whole step melts the skin, which
# then draws metal at T_rigid and melts that too, so the heat sent down falls from +5 to -7 MW/m2 over 30 K (the
# plate's centre at 0.5 s, measured). A base law fitted at one temperature cannot follow it, and the 3D surface, whose
# admittance is about four times the law's slope, overshoots it: undamped repetition alternates between two states.
# Damping moves the linearisation temperature only a fraction omega of the way to the 3D answer, Aitken's method
# choosing omega from the last two passes.
DAMPED_MAX_PASSES = 12            # passes of a macro step at most (the last pass's mismatch is deferred as before)
DAMPED_OMEGA0 = 0.5               # the first pass's relaxation factor
DAMPED_OMEGA_RANGE = (0.02, 1.0)  # Aitken's factor kept within these
DAMPED_TOL_K = 0.05               # K: converged when no base temperature moves more than this between passes
# The move limit (2026-10-12). The metal a base sweeps in one step is valued along a straight 3D profile, and the moving
# mesh refuses moves of more than about half an element (1.3 mm on 2 mm elements, 5.2 mm on 4 mm, measured): a step
# whose draws would move a base further than this is split; a hand-down is capped at it (2026-10-27).
MOVE_MAX = 0.25                   # of the mean surface edge: the most a base moves in one part (draws split a part
                                  # beyond it, hand-downs are capped at it, absorption takes what they leave)
HAND_DOWN_SHARE = 1.0 - 1e-6      # of MOVE_MAX: the most a part's hand-down moves a base (Asha, 2026-10-27): just under
                                  # the whole allowance, so that the round-off of the move ratio, which divides the
                                  # capped mass by the area again, never makes a capped hand-down split a part
SPLIT_DEPTH_MAX = 6               # a part is split at most this many times over
MOVE_HINT_MARGIN = 1.25           # the next macro step's parts: this times the parts this one's largest move called for
# The part length while melting (2026-10-12). On the 100 mm sphere's first 1.5 s of melting, 0.5 s parts cost 1 454 s
# (53 % on parts rejected for their moves, 4.4 passes a kept part) against 341 s for 0.0625 s macro steps (1.8 passes);
# before melting short parts cost about 18 % more, so the cap holds only while some skin's top is non-rigid. (Those
# steps also removed 21.28 g against 16.6 g, but because ABSORB_MAX is per macro step, not because the parts were short.)
PART_MAX = 0.0625                 # s: the longest part while any skin's top is at or above T_rigid, or about to be
SKIN_CG_TOL = 1e-13               # the 3D solver's linear tolerance in a skin run: at its default 1e-10 the skin body's
                                  # books drifted by 1e-8 of the heat absorbed over 17 short parts (about 5e-5 J a solve
                                  # on the 4 mm test sphere); at 1e-13, as with a direct solve, 2e-11, at no cost
                                  # (measured 2026-10-14). The element model's runs keep their solver's own tolerance.


SLURRY_SIGMA = 0.80               # N/m: the exposed slurry's surface tension in the cascade (the 2026-10-07 design 8.2, the
                                  # large-fragment design's labelled assumption for slurry)


@dataclasses.dataclass
class FlowForcing:
    """The gas side's push on the skins' metal at a macro step, per corner: `dp_ds` [Pa/m, (n,)] the wall pressure's
    gradient along the meridian and `tau` [Pa, (n,)] the gas shear (the surface flow's), and `g` [m/s2, (3,)] the
    acceleration the metal feels in the body's frame -- minus the body's own: the deceleration along the flight
    direction. `result`, the surface flow's whole evaluation (a `surface_flow.SurfaceFlowResult` on the corners), is what
    Girin's spraying reads (sub-plan 22 Task 1c); None for a push given without one."""
    dp_ds: np.ndarray
    tau: np.ndarray
    g: np.ndarray
    result: object = None

    def toward(self, other, frac):
        """This forcing moved the fraction `frac` of the way to `other`, linearly."""
        return FlowForcing(*(a + (b - a) * frac for a, b in ((self.dp_ds, other.dp_ds), (self.tau, other.tau),
                                                             (self.g, other.g))))


@dataclasses.dataclass
class SprayPlan:
    """Girin's spraying from the top cells over one macro step, evaluated at its start (`girin_plan`; spec 8, sub-plan 22
    Task 1c), per column: `rate` [kg/(m2 s) of top area] the liquid's mode's and `rate_slurry` the exposed slurry's own
    (the cascade; zero without slurry spraying); `result` the liquid's `spray.SprayResult` (branch, radius, wavelength);
    `slurry` the exposed slurry's per column -- its `branch`, radius `r` [m], `wavelength` [m], viscosity `mu` [Pa s] at
    the shear rate `gamma` [1/s] of its own sheared layer, density `rho` [kg/m3], conjugate depth `delta_m` [m], layer
    `depth` [m] from its top to T_rigid, the viscous film's `gamma_tau` and floor `h_floor` [m], its `we_s` and the
    cell's temperature `T` [K] (nan, -1 where there is none); the depths read
    at the start, `flowing`, `liquid` and the cap's `cap_depth` [m], and `on_slurry`, where a flowing cell that is not
    fully liquid lies beneath the liquid run. The record the macro step's accepted passes fill (`book`), per column and
    mode -- the liquid's and the exposed slurry's separate -- from the release by the cells it came from: `mass` and
    `mass_slurry` [kg], `enthalpy` and `enthalpy_slurry` [J], `solid` and `solid_slurry` [kg] and `momentum` and
    `momentum_slurry` [kg m/s, (n, 3)]; and, at the step's end, `number` and `number_slurry`, the droplets: the mass over
    one droplet's at the mode's radius. Where the plan was evaluated, for the records' rows (`source_block`; sub-plan 22
    Task 2): `theta` [rad], `closure` (the surface flow's), `points` [m, (n, 3)] the top surface's corners, `state` the
    aero state and `props` the liquid's properties."""
    rate: np.ndarray
    rate_slurry: np.ndarray
    result: object
    slurry: dict
    flowing: np.ndarray
    liquid: np.ndarray
    cap_depth: np.ndarray
    on_slurry: np.ndarray
    mass: np.ndarray = None
    mass_slurry: np.ndarray = None
    number: np.ndarray = None
    number_slurry: np.ndarray = None
    theta: np.ndarray = None
    closure: np.ndarray = None
    points: np.ndarray = None
    state: object = None
    props: object = None

    def __post_init__(self):
        T = self.slurry["T"]
        self._cascade = skin_mod.Cascade(self.rate, self.rate_slurry, np.where(np.isfinite(T), T, np.inf),
                                         histogram=self.histogram)
        n = len(self.rate)
        self.mass, self.mass_slurry = np.zeros(n), np.zeros(n)
        self.number, self.number_slurry = np.zeros(n), np.zeros(n)
        self.enthalpy, self.enthalpy_slurry = np.zeros(n), np.zeros(n)
        self.solid, self.solid_slurry = np.zeros(n), np.zeros(n)
        self.momentum, self.momentum_slurry = np.zeros((n, 3)), np.zeros((n, 3))

    def release(self, sk):
        """`skin.Flow.release`: the liquid's rate and the exposed slurry's (`skin.Cascade`), whatever the skins hold within
        the macro step -- Girin's modes and radii are the macro step's, from its start (spec 8) -- the slurry's for slurry
        cells at least as warm as the cell that gave it."""
        return self._cascade

    def droplet_mass(self, rho_liquid):
        """One droplet's mass [kg, (2, n)] at each mode's radius, the liquid's mode's at the density `rho_liquid`, the
        slurry's at its cell's; nan where the mode has no droplet."""
        r = np.vstack([self.result.r, self.slurry["r"]])
        rho = np.vstack([np.broadcast_to(np.asarray(rho_liquid, dtype=float), r[0].shape), self.slurry["rho"]])
        ok = np.isfinite(r) & (r > 0.0)
        return np.where(ok, 4.0 / 3.0 * np.pi * rho * np.where(ok, r, 1.0) ** 3, np.nan)

    def count(self, rho_liquid):
        """The record's droplets: each mode's mass over one droplet's at its radius."""
        one = self.droplet_mass(rho_liquid)
        for j, (mass, key) in enumerate(((self.mass, "number"), (self.mass_slurry, "number_slurry"))):
            ok = (mass > 0.0) & np.isfinite(one[j])
            setattr(self, key, np.where(ok, mass / np.where(ok, one[j], 1.0), 0.0))

    def histogram(self, mass):
        """(droplets, mass [kg]) per size bin (spray.BIN_EDGES) of a release by mode [kg, (2, n): the liquid's mode's and
        the slurry's] at the plan's radii -- each flow sub-step's release (`skin.Cascade.histogram`; sub-plan 22 Task 2),
        which, a column's radius being the macro step's, is the histogram of the step's combined records to round-off."""
        mass = np.asarray(mass, dtype=float)
        one = self.droplet_mass(self.props.rho)
        ok = (mass > 0.0) & np.isfinite(one)
        number = np.where(ok, mass / np.where(ok, one, 1.0), 0.0)
        r = np.vstack([self.result.r, self.slurry["r"]])
        return spray_mod.histogram(np.where(ok, r, np.nan).ravel(), number.ravel(), mass.ravel())

    def book(self, fb):
        """Add an accepted pass's release (`skin.FlowBook`, by the cells it came from: the liquid's mode took the fully
        liquid cells, the slurry's the others) to the record."""
        self.mass += fb.phase_mass[0]
        self.mass_slurry += fb.phase_mass[1]
        self.enthalpy += fb.phase_H[0]
        self.enthalpy_slurry += fb.phase_H[1]
        self.solid += fb.phase_solid[0]
        self.solid_slurry += fb.phase_solid[1]
        self.momentum += fb.phase_p[0]
        self.momentum_slurry += fb.phase_p[1]

    def source_block(self, t):
        """The macro step's records as source-table rows (spray.SOURCE_COLUMNS, by `spray.source_block`), one per column
        and mode that released, the liquid's mode's first: at the end of the step `t` [s] and the plan's aero state, the
        column's top corner `points` and angle `theta` where the plan was evaluated, its `closure`; the mode's branch,
        radius, delta_m, We_s, the viscous film's floor and gamma tau*; b, for the liquid's mode the flowing depth and for
        the slurry's its own layer, the slurry from its top to T_rigid; the droplets (`count` first) and their mass; their
        speed and direction the mean velocity of the metal released, its momentum over its mass, the momentum the books
        lose with the droplets (spec 6.4: they leave at the speed of the cells they are made of; zero direction where it
        is at rest); We_d, Oh and the breakup flag at their radius, by their phase's density, viscosity and surface tension
        (the slurry's mixture density, viscosity at its sheared layer's rate and SLURRY_SIGMA, the 2026-10-07 design 8.4);
        and their solid fraction and enthalpy per kilogram, mass-weighted over the cells released."""
        st, s = self.state, self.slurry
        flow = types.SimpleNamespace(closure=self.closure)
        blocks = []
        for mass, H, solid, p, number, res, b, props in (
                (self.mass, self.enthalpy, self.solid, self.momentum, self.number, self.result, self.flowing, self.props),
                (self.mass_slurry, self.enthalpy_slurry, self.solid_slurry, self.momentum_slurry, self.number_slurry,
                 types.SimpleNamespace(**{k: s[k] for k in ("branch", "r", "delta_m", "we_s", "h_floor", "gamma_tau")}),
                 s["depth"], material_mod.LiquidProperties(s["rho"], s["mu"], SLURRY_SIGMA))):
            speed = np.linalg.norm(p, axis=1)
            on = mass > 0.0
            per = lambda x: np.divide(x, mass, out=np.zeros_like(x), where=on)
            direction = np.divide(p, speed[:, None], out=np.zeros_like(p), where=speed[:, None] > 0.0)
            rows = types.SimpleNamespace(branch=res.branch, r=res.r, delta_m=res.delta_m, we_s=res.we_s, h_floor=res.h_floor,
                                         gamma_tau=res.gamma_tau, dm=mass, dn=number, v_s=per(speed))
            blocks.append(spray_mod.source_block(t, st.h, st.V, self.theta, self.points, direction, flow, rows, b, props, st,
                                                 solid_fraction=per(solid), h_kg=per(H)))
        return np.concatenate(blocks, axis=0)


def region_extent(edges, selected, areas):
    """Lateral extent [m] of each column's own contiguous selected region on the columns' graph (`edges`, (E, 2)), as the
    region's equivalent diameter 2 sqrt(A / pi) with A the selected columns' `areas` summed; zero where not selected: the
    element model's `MeltingBody.region_extent` on the corners, the domain a surface wave must fit inside (2026-09-24)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    sel = np.asarray(selected, dtype=bool)
    n = len(sel)
    out = np.zeros(n)
    if not sel.any():
        return out
    i, j = edges[:, 0], edges[:, 1]
    keep = sel[i] & sel[j]
    i, j = i[keep], j[keep]
    g = coo_matrix((np.ones(2 * i.size), (np.concatenate([i, j]), np.concatenate([j, i]))), shape=(n, n))
    count, label = connected_components(g, directed=False)
    area = np.bincount(label[sel], np.asarray(areas, dtype=float)[sel], count)
    out[sel] = 2.0 * np.sqrt(area[label[sel]] / np.pi)
    return out


def _rows(flow, idx):
    """The surface flow's evaluation on the corners idx alone: every per-corner array taken at idx, the rest kept."""
    n = len(flow.closure)
    out = copy.copy(flow)
    for key, value in vars(flow).items():
        if isinstance(value, np.ndarray) and value.ndim >= 1 and value.shape[0] == n:
            setattr(out, key, value[idx])
    return out


def girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False, points=None):
    """Girin's spraying from each column's top cells for the macro step about to run (spec 8; sub-plan 22 Task 1c), from
    the surface flow's evaluation on the corners (`flow`), the aero state (`state`, for the rarefied branch), the shared
    spray model (`model`) and the skins as they are at its start: a SprayPlan.

    The liquid's mode: b and the regime test read the flowing depth `sk.depth_above(T_rigid)` (spec 6.7), and so does
    the viscous film's floor (Asha, 2026-10-28); delta_m is Girin's Eq. (2) and his surface velocity his closure's
    (`film.lubrication` on b, as the element model's); the front-surface Rayleigh-Taylor mode reads the liquid pool, the
    liquid depth `sk.depth_above(T_liquidus)` (`rt_layer`), and every wave must fit in the contiguous region of liquid
    columns, and where Girin's closure does not apply a liquid lying on slurry takes no shear mode (`on_slurry`) -- the
    element model's reading (amendment of 2026-10-06; Asha, 2026-10-29, option B); the droplet is capped by one crest of
    its wave on the depth the release may take, the liquid's, or with `slurry` the flowing depth (2026-10-24); a quarter
    of `radius` caps it too. Columns on the windward face, theta <= 90 degrees, only; `A_top` [m2] the top areas, `edges`
    the columns' graph, `dt` [s] the macro step (the nominal release in the result). With `slurry`, the exposed slurry's
    own mode too (`slurry_modes`). The plan keeps where it was evaluated for its records' rows: `theta`, the closure,
    the aero state, the liquid's properties and `points` [m, (n, 3)], the columns' top corners (sub-plan 22 Task 2)."""
    mat, liq = sk.material, model.liquid
    n = sk.n
    F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
    flowing = sk.valid & np.isfinite(skin_mod.cell_viscosity(sk, 0.0))
    liquid = flowing & (sk.T >= mat.T_liquidus)
    n_flow = np.cumprod(flowing, axis=1).sum(axis=1)
    n_liq = np.cumprod(liquid, axis=1).sum(axis=1)
    on_slurry = n_flow > n_liq
    theta = np.asarray(theta, dtype=float)
    windward = theta <= 0.5 * np.pi
    w_n = flow.deceleration * np.maximum(np.cos(theta), 0.0)
    delta_m, _ = spray_mod.melt_layer(flow, liq)
    v_s = film_mod.lubrication(flow.tau, flow.G, F, delta_m, liq.mu, b_layer=F)[0]
    extent = region_extent(edges, L > 0.0, A_top)
    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
    cap = F if slurry else L
    res = model.evaluate(flow, state, F, delta_m, v_s, windward, dt, A_top, (sk.m * run).sum(axis=1), radius,
                         b_layer=F, extent=extent, deceleration_n=w_n, regime_layer=F, on_slurry=on_slurry, rt_layer=L,
                         cap_depth=cap)
    if slurry:
        rate_s, modes = slurry_modes(sk, flow, state, model, windward, w_n, A_top, edges, radius, dt, F, n_liq, n_flow)
    else:
        rate_s, modes = np.zeros(n), _no_slurry(n)
    return SprayPlan(res.mdot, rate_s, res, modes, F, L, cap, on_slurry, theta=theta, closure=np.asarray(flow.closure),
                     points=points, state=state, props=liq)


def _no_slurry(n):
    out = {key: np.full(n, np.nan) for key in ("r", "wavelength", "mu", "gamma", "rho", "delta_m", "depth", "gamma_tau",
                                               "T", "we_s", "h_floor")}
    out["branch"] = np.full(n, -1)
    return out


def slurry_modes(sk, flow, state, model, windward, w_n, A_top, edges, radius, dt, F, n_liq, n_flow):
    """The exposed slurry's own mode per column, the cascade of spec 8 (the 2026-10-07 design's 8.1-8.2; Asha, 2026-10-29,
    option i): (rate [kg/(m2 s)], the per-column dict of SprayPlan.slurry). Where a flowing cell that is not fully liquid
    lies beneath the liquid run on the windward face, its cell -- the slurry the release exposes first -- gives the
    properties: Li et al.'s viscosity (`rheology.viscosity`, with the bridge to the liquid) at the cell's solid fraction
    and the shear rate of the slurry's own sheared layer, V_s over delta_m or over the layer, found by fixed-point
    iteration from rest (at most skin.SHEAR_PASSES passes, until no viscosity moves by skin.SHEAR_TOL, the rate held at
    rheology.GAMMA_MAX); the cell's mixture density; SLURRY_SIGMA. The layer runs from the cell's top to T_rigid, and the
    regime test, the viscous film's floor and b are it. Girin's formulas are the shared model's with these properties
    (its k_r, k_t, critical Weber number and dispersion table), the Rayleigh-Taylor mode left to the liquid; the wave
    must fit in the flowing region; the cap is the flowing depth. Evaluated once a macro step, as the liquid's mode: its
    rate holds for slurry cells at least as warm as the one whose properties gave it, and a colder cell, more solid and
    more viscous, waits for the next macro step's evaluation (`skin.Cascade.T_slurry`)."""
    mat = sk.material
    n = sk.n
    out = _no_slurry(n)
    rate = np.zeros(n)
    rows = np.flatnonzero((n_flow > n_liq) & windward)
    if not rows.size:
        return rate, out
    k = n_liq[rows]
    T = sk.T[rows, k]
    dz = sk.dz()
    depth = F[rows] - (np.cumsum(dz, axis=1) - dz)[rows, k]
    f_s = 1.0 - sk.props.liquid_fraction(T)
    rho = sk.density(T)
    f_bridge = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    sub = _rows(flow, rows)
    viscosity = lambda gamma: rheology.viscosity(T, f_s, gamma, mat.liquid.mu, f_bridge)

    def layer(mu):
        d, _ = spray_mod.melt_layer(sub, material_mod.LiquidProperties(rho, mu, SLURRY_SIGMA))
        v, _, shear, _ = film_mod.lubrication(sub.tau, sub.G, depth, d, mu, b_layer=depth)
        return d, v, np.minimum(shear, rheology.GAMMA_MAX)

    gamma = np.zeros(rows.size)
    mu = viscosity(gamma)
    for _ in range(skin_mod.SHEAR_PASSES):
        gamma = layer(mu)[2]
        mu_new = viscosity(gamma)
        done = np.all(np.abs(mu_new - mu) <= skin_mod.SHEAR_TOL * mu)
        mu = mu_new
        if done:
            break
    d_m, v_s, _ = layer(mu)
    extent = region_extent(edges, F >= spray_mod.B_MIN, A_top)
    cells = np.arange(sk.m.shape[1])[None, :]
    m_s = (sk.m * ((cells >= n_liq[:, None]) & (cells < n_flow[:, None]))).sum(axis=1)     # the slurry run, nominally
    for j, i in enumerate(rows):
        props = material_mod.LiquidProperties(float(rho[j]), float(mu[j]), SLURRY_SIGMA)
        m = spray_mod.SprayModel(props, model.k_r, model.k_t, model.we_critical, model.table, rt_spray=False)
        one, cell = slice(j, j + 1), slice(i, i + 1)
        r = m.evaluate(_rows(sub, one), state, depth[one], d_m[one], v_s[one], windward[cell], dt, A_top[cell], m_s[cell],
                       radius, b_layer=depth[one], extent=extent[cell], deceleration_n=w_n[cell], regime_layer=depth[one],
                       cap_depth=F[cell])
        rate[i] = r.mdot[0]
        out["branch"][i], out["r"][i], out["wavelength"][i] = r.branch[0], r.r[0], r.wavelength[0]
        out["gamma_tau"][i], out["we_s"][i], out["h_floor"][i] = r.gamma_tau[0], r.we_s[0], r.h_floor[0]
    out["mu"][rows], out["gamma"][rows], out["rho"][rows], out["T"][rows] = mu, gamma, rho, T
    out["delta_m"][rows], out["depth"][rows] = d_m, depth
    return rate, out


def jump_faults(T_bottom, T_node, T_deep, crossed, tol=None, handed=None):
    """The columns whose skin's bottom cell lies more than `tol` (JUMP_TOLERANCE) outside the temperatures of the metal
    that crossed their base (spec 6.6): where metal was drawn or absorbed (`crossed`), between the base node's and the
    3D profile's at the depth swept, `T_deep` (equal to T_node where nothing was taken from the 3D model); with `handed`,
    (mask, coldest, warmest[, uncovered]), where metal was handed down, between the handed cells' coldest and warmest
    temperatures and that of the cell the hand-down uncovered, as the hand-down left it; a column that did both is a
    fault only outside both. Metal drawn from below the node is colder than the node by the profile's slope times the
    depth, which is not a fault (Asha, 2026-10-12: comparing with the node alone counted 812-1 694 faults on the level-2
    plate at 0.1-0.5 s steps, none at 0.05 s). The body reads the bottom cell after the reconciliation and judges a
    hand-down by the handed cells and the cell it uncovered (Asha, 2026-10-26, sub-plan 21 Task 9: judged by the node, a
    hand-down counted the warmer cell it uncovers where runoff deposits, 199 faults on a reduced pushed plate to 1.0 s,
    rather than the 270 bottom cells the reconciliation then superheated; and judged by the handed cells alone, after the
    build, it still counted that cell where the skin warms steeply above its base, 1 900 of the 1 942 faults on Task 8's
    pushed plate to 1.5 s, none of them heated by the reconciliation)."""
    tol = skin_mod.JUMP_TOLERANCE if tol is None else tol
    lo, hi = np.minimum(T_node, T_deep), np.maximum(T_node, T_deep)
    inside = crossed & (T_bottom >= lo - tol) & (T_bottom <= hi + tol)
    judged = crossed
    if handed is not None:
        mask, coldest, warmest = handed[:3]
        if len(handed) > 3:                          # the cell it uncovered, as the hand-down left it, widens the range
            coldest, warmest = np.fmin(coldest, handed[3]), np.fmax(warmest, handed[3])
        inside = inside | (mask & (T_bottom >= coldest - tol) & (T_bottom <= warmest + tol))
        judged = judged | mask
    return judged & ~inside


class SkinBody(ThermalBody):
    """A ThermalBody whose surface is resolved by corner columns with layered skins on a moving mesh (module docstring).
    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass;
    `flow`: the metal flows (sub-plan 21), `flow_cycles` flow sub-steps in each heat sub-step (None: as many as
    `skin.flow_cycles` chooses each heat sub-step, Task 6d), driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its
    defaults if None) from each macro step's aero state; `release`, with the flow on, a release from the top cells in
    every flow sub-step (`skin.Flow.release`: a callable of the skins giving a rate [kg/(m2 s)] per unit top area per
    column; sub-plan 21 Task 7), the exposed slurry released too where `slurry_spray`; or `spray`, a `spray.SprayModel`
    with the flow on, Girin's spraying from the top cells evaluated at the start of every macro step (`girin_plan`; sub-plan
    22 Task 1c), with `slurry_spray` the exposed slurry's own mode too (the cascade), and its record in `spray_plan`; its
    droplet records, one block of source-table rows per macro step, in `source_rows` where `particles`, and the size
    histograms in `hist_n` and `hist_m` (sub-plan 22 Task 2: the element model's interface to the coupled writer)."""

    def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                 v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                 repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
                 surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True):
        if repeat_rule not in REPEAT_RULES:
            raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
        self.flow_cycles = flow_cycles
        if not (flow_cycles is None or isinstance(flow_cycles, (int, np.integer)) and flow_cycles >= 1):
            raise ValueError("flow_cycles must be None or an integer of at least 1, got {!r}".format(flow_cycles))
        self.flow_on, self.surface_flow = bool(flow), surface_flow
        if release is not None and not self.flow_on:
            raise ValueError("a release from the top cells needs the flow on: without it the heat-only skins remove the "
                             "liquid as it forms (stage 2's stand-in for flow and spray)")
        if spray is not None and not self.flow_on:
            raise ValueError("Girin's spraying from the top cells needs the flow on: without it the heat-only skins "
                             "remove the liquid as it forms")
        if spray is not None and release is not None:
            raise ValueError("a skin body sprays by Girin's model or by a release given, not both")
        self.release, self.slurry_spray = release, bool(slurry_spray)
        self.spray, self.spray_plan = spray, None            # Girin's model, and its plan of the last macro step (Task 1c)
        self.particles = bool(particles)                     # keep the droplet records' rows (Task 2)
        self.source_rows = []                                # one block of source-table rows per macro step (Task 2)
        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)   # droplets, kg per size bin
        self._step_release = release                         # what the parts' Flow releases this macro step
        self.repeat_rule = repeat_rule
        self.part_max = part_max                             # s, or None: no cap on the parts while melting
        self.material, self.solver, self.mass_kg = material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient, self.remove, self.demise_fraction = float(T_ambient), bool(remove), float(demise_fraction)
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.settings = settings or skin_mod.SkinSettings()
        if getattr(solver, "cg_tol", None) is not None:
            solver.cg_tol = min(solver.cg_tol, SKIN_CG_TOL)
        solver.setup(mesh, material, 0.0)                    # the 3D boundary lies under the skins: it does not radiate
        solver.set_temperature(T0)
        self.mass0 = material.rho * mesh.volume()
        self.moving = movemesh.MovingMesh(mesh, solver, material, depth=deform_depth)
        self.columns = columns_mod.Columns(mesh.points[self.moving.corners], self.moving.faces_corner)
        self.skins = skin_mod.Skins(material, self.columns.area, self.settings)
        self.pending = np.zeros(self.columns.n)              # J per corner: the deferred interface mismatch
        self.column_removed = np.zeros(self.columns.n)       # kg per column removed as liquid since time zero
        self.integrated_heat = self.absorbed_heat = self.radiated_heat = 0.0
        self.removed_mass = self.removed_enthalpy = 0.0
        self.drawn_mass = self.absorbed_mass = self.handed_mass = 0.0
        self.held_mass = self.held_depth_max = 0.0           # kg, m: what the tolerance and the cap hold at the last move
        self.jump_faults = self.interface_repeats = self.absorb_capped = 0
        self.last_mismatch = self.last_shortfall = 0.0
        self.last_trigger = {"columns_over": 0, "worst_ratio": 0.0, "total_ratio": 0.0}    # the first pass's mismatch
        self.last_passes, self.last_total = 1, 0.0           # the damped rule's passes and final total mismatch ratio
        self.step_parts, self._hint, self._ratio = 1, 1, 0.0   # parts of the last macro step; of the next; its move ratio
        self._q_prev = None                                  # the last macro step's heating, for the interpolation in time
        self._top_start = None                               # the tops' temperatures [K] at the last macro step's start
        self._dt_last = None                                 # s, and its length: the look-ahead of the part cap
        self.capped_ahead = 0                                # macro steps capped because a top was about to be non-rigid
        self._heating = None                                 # this macro step's (t0, dt, q at its start, q at its end)
        self.wall = {"skins": 0.0, "3d": 0.0, "move": 0.0, "rebuild": 0.0, "flow": 0.0}
        self.last_flow = None                                # the surface flow's last evaluation (SurfaceFlowResult)
        self._forcing = self._forcing_prev = None            # this macro step's (start, end) FlowForcing; the last end
        self.impulse, self.momentum_lost = np.zeros(3), np.zeros(3)   # the momentum books since time zero [N s, kg m/s]
        self.impulse_size = 0.0                              # N s: the push's, shear's and walls' impulses, in size
        self.released_mass = self.released_enthalpy = 0.0    # kg, J: released from the top cells since time zero (Task 7)
        self.droplet_number = 0.0                            # droplets Girin's spraying released since time zero (Task 1c)
        self.droplet_momentum = np.zeros((self.columns.n, 3))   # kg m/s: what the droplets took from each column, and
        self.droplet_shear = np.zeros((self.columns.n, 3))      # N s: the gas shear's impulse on its top, this macro step
        self.flow_courant = self.flow_growth = 0.0           # this macro step's worst flow sub-step
        self.flow_passes = self.flow_cycles_max = 0          # and the most flow sub-steps in a heat sub-step
        self._volume = None                                  # m3, the 3D model's volume since the last move (mass)
        self.iterations, self.last, self.last_move = [], None, None
        self._fill(T0)
        self.energy0 = self.energy()
        self.momentum0 = self.skins.momentum()
        self._top0, self._normal0 = self.surface.points.copy(), self.surface.normals.copy()

    # -- set-up -----------------------------------------------------------------------------------------------------
    def _fill(self, T0):
        """Each column takes its minimum mass from the 3D model (spec 6.1): the corners move in by it, and the skins hold
        exactly what they swept, at the initial temperature."""
        m_min = skin_mod.minimum_mass(self.material, self.columns.area, self.settings)
        u = -self.columns.normal * (m_min / (self.material.rho * self.columns.area))[:, None]
        r = self.moving.move_with_policy(u)
        self.columns.update(self.moving.mesh.points[self.moving.corners])
        self.skins.area = self.columns.area.copy()
        self.skins.fill(self.material.rho * r.swept_volume, T0)
        self.skins.exchange_bottom(np.zeros(self.columns.n), r.swept_enthalpy - self.skins.H.sum(axis=1))
        self._refresh_geometry()                       # the top surface: the bases moved out by the filled skins

    def _refresh_geometry(self):
        """The columns, their areas and the top surface after the corners moved; the 3D model's volume is measured again
        when next asked for (`mass`)."""
        mm = self.moving
        self._volume = None
        self.columns.update(mm.mesh.points[mm.corners])
        self.skins.area = self.columns.area.copy()
        self.surface = columns_mod.CornerSurface.from_points(self.columns.top_points(self.skins.thickness()),
                                                             self.columns.faces)
        self.theta = self.surface.angles_to(self.v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)

    @property
    def mesh(self):
        return self.moving.mesh

    # -- the macro step ---------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None, forcing=None):
        """The macro step from t - dt to t (module docstring) as equal parts, as many as the last macro step's largest
        move called for. With the flow on, the gas side's push is the surface flow's on the corners from the aero state
        `state` (`flow_forcing`), or `forcing`, a FlowForcing given instead (a prescribed push, in the tests)."""
        n, self._ratio, self.step_parts = self._hint, 0.0, 0
        if self.flow_on:
            end = forcing if forcing is not None else self.flow_forcing(state)
            last = self._forcing_prev
            start = last if last is not None and last.dp_ds.shape == end.dp_ds.shape else end
            self._forcing, self._forcing_prev = (start, end), end
            self.flow_courant = self.flow_growth = 0.0
            self.flow_passes = self.flow_cycles_max = 0
            self.droplet_momentum, self.droplet_shear = np.zeros((self.columns.n, 3)), np.zeros((self.columns.n, 3))
            if self.spray is not None:                       # Girin once, at the macro step's start (sub-plan 22 Task 1c)
                if end.result is None or state is None:
                    raise ValueError("Girin's spraying reads the surface flow's evaluation and the aero state (its records' "
                                     "altitude and speed): give the aero state, with a forcing only if it has its result")
                self.spray_plan = girin_plan(self.skins, end.result, state, self.spray, self.theta, self.surface.areas,
                                             self.columns.edges, self.transverse_radius(), dt, self.slurry_spray,
                                             points=self.surface.points.copy())
                self._step_release = self.spray_plan.release
        sk, mat = self.skins, self.material
        top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
        if self.part_max is not None and mat.melts:
            if np.any(sk.valid[:, 0] & (sk.T[:, 0] >= mat.T_rigid)):
                n = max(n, int(math.ceil(dt / self.part_max - 1e-9)))
            elif self._top_start is not None and self._top_start.shape == top.shape:
                with np.errstate(invalid="ignore"):          # a top that rose over the last step, carried on over this one
                    ahead = top + (top - self._top_start) * (dt / self._dt_last)
                if np.any(np.nan_to_num(ahead, nan=-np.inf) >= mat.T_rigid):
                    n = max(n, int(math.ceil(dt / self.part_max - 1e-9)))
                    self.capped_ahead += 1
        self._top_start, self._dt_last = top, dt
        q_end = np.asarray(loads.q_conv, dtype=float)
        q_start = self._q_prev if self._q_prev is not None and self._q_prev.shape == q_end.shape else q_end
        self._heating = (t - dt, dt, q_start, q_end)
        for k in range(n):
            self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt, loads, 0 if n == 1 else 1)
        self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * self._ratio)))
        self._q_prev = q_end.copy()
        if self.spray is not None and self.flow_on:          # the macro step's droplets, from its accepted passes' mass
            self.spray_plan.count(self.spray.liquid.rho)
            self.droplet_number += float(self.spray_plan.number.sum() + self.spray_plan.number_slurry.sum())
            if self.particles:                               # and its records' rows (Task 2)
                block = self.spray_plan.source_block(t)
                if len(block):
                    self.source_rows.append(block)

    def flow_forcing(self, state):
        """The gas side's push on the metal at this macro step (FlowForcing): the surface flow evaluated on the corners
        from the aero state, as the heating is -- the top surface's angles to the flight direction, the nose radius and
        the tops' temperatures at the step's start -- and the deceleration along the flight direction."""
        if state is None:
            raise ValueError("a skin body with the flow on needs the macro step's aero state, or a forcing given")
        if self.surface_flow is None:
            from . import surface_flow
            self.surface_flow = surface_flow.SurfaceFlow()
        mat = self.material
        rho_l = mat.liquid.rho if mat.liquid is not None else mat.rho
        res = self.surface_flow.evaluate(state, self.theta, self.nose_radius(), rho_l, self.surface_temperature())
        self.last_flow = res
        return FlowForcing(np.asarray(res.dp_ds, dtype=float), np.asarray(res.tau, dtype=float),
                           res.deceleration * self.v_hat, res)

    def _advance_part(self, t, dt, dt_macro, loads, depth):
        """One part of a macro step, ending at t: the passes, then -- if its draws would move a base more than MOVE_MAX
        of the mean surface edge (its hand-downs are capped within it) -- the part again as shorter parts, else the
        books and the move."""
        mat, sk, mm = self.material, self.skins, self.moving
        t0, dt_m, q_start, q_end = self._heating
        frac = ((t - 0.5 * dt) - t0) / dt_m
        q = q_start + (q_end - q_start) * frac                              # the heating at the part's middle
        A_top = self.surface.areas
        flow = None
        if self.flow_on:                                                    # the gas side's push at the part's middle
            f = self._forcing[0].toward(self._forcing[1], frac)
            sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
            flow = skin_mod.Flow(self.columns, f.dp_ds, skin_mod.along_meridian(self.columns.normal, self.v_hat, f.tau),
                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles, release=self._step_release,
                                 slurry=self.slurry_spray)
        eps = np.full(sk.n, self.emissivity)
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
        nodes = mm.corners
        T_b0 = self.solver.temperature()[nodes]
        g0 = self._normal_gradient()                         # the profile's slope below every base, for this part
        load = np.zeros(mm.mesh.n_nodes)
        load[nodes] = self.pending / dt
        solver_state, skin_state = self.solver.state(), sk.state()
        T_lin = T_b0
        damped = self.repeat_rule == "damped"
        omega, r_prev = DAMPED_OMEGA0, None
        flow_seconds = sum(skin_mod.FLOW_SECONDS.values())
        for attempt in range(DAMPED_MAX_PASSES if damped else 2):
            started = time.perf_counter()
            a, b = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_lin, m_min, side=side, remove=self.remove,
                                    g_base=g0, flow=flow)
            self.wall["skins"] += time.perf_counter() - started
            started = time.perf_counter()
            res = self.solver.step(dt, np.zeros(len(self.solver.areas)), self.T_ambient, nodal_load=load,
                                   nodal_law=(nodes, a, b))
            self.wall["3d"] += time.perf_counter() - started
            T_b1 = self.solver.temperature()[nodes]
            started = time.perf_counter()
            book = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_b0, m_min, T_b1=T_b1, side=side,
                                    remove=self.remove, g_base=g0, flow=flow)
            self.wall["skins"] += time.perf_counter() - started
            mismatch = book.E_base - (a + b * T_b1) * dt
            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
            ratio = np.abs(mismatch) / scale
            total = float(np.abs(mismatch).sum() / scale.sum())
            if attempt == 0:
                self.last_trigger = {"columns_over": int((ratio > skin_mod.INTERFACE_TOLERANCE).sum()),
                                     "worst_ratio": float(ratio.max()), "total_ratio": total}
            if damped:
                r = T_b1 - T_lin
                self.last_passes, self.last_total = attempt + 1, total
                if (total <= skin_mod.INTERFACE_TOLERANCE or np.abs(r).max() <= DAMPED_TOL_K
                        or attempt == DAMPED_MAX_PASSES - 1):
                    break
                if r_prev is not None:                       # Aitken: omega from the change of the residual
                    dr = r - r_prev
                    denom = float(dr @ dr)
                    if denom > 0.0:
                        omega = float(np.clip(-omega * float(r_prev @ dr) / denom, *DAMPED_OMEGA_RANGE))
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin, r_prev = T_lin + omega * r, r
                self.interface_repeats += 1
                continue
            repeat = {"column": bool((ratio > skin_mod.INTERFACE_TOLERANCE).any()),
                      "total": total > skin_mod.INTERFACE_TOLERANCE, "never": False}[self.repeat_rule]
            if attempt == 0 and repeat:
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin = T_b1
                self.interface_repeats += 1
                continue
            break
        self.wall["flow"] += sum(skin_mod.FLOW_SECONDS.values()) - flow_seconds
        ratio = self._move_ratio(book, m_min)
        if ratio > 1.0 and depth < SPLIT_DEPTH_MAX:          # too far for one part: run it again as shorter parts
            self.solver.set_state(solver_state)
            sk.set_state(skin_state)
            n = int(math.ceil(ratio))
            for k in range(n):
                self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt_macro, loads, depth + 1)
            return
        self._ratio = max(self._ratio, ratio * dt_macro / dt)
        self.step_parts += 1
        self.pending = mismatch
        self.last_mismatch = float(np.abs(mismatch).sum())
        E_top, E_rad = float(book.E_top.sum()), float(book.E_rad.sum())
        self.integrated_heat += E_top
        self.absorbed_heat += E_top - E_rad
        self.radiated_heat += E_rad
        self.removed_mass += float(book.removed.sum())
        self.column_removed += book.removed
        self.removed_enthalpy += float(book.removed_H.sum())
        self.drawn_mass += float(book.drawn.sum())
        self.momentum_lost += book.momentum_out.sum(axis=0)
        if book.flow is not None:                            # the flow's impulses (the transport's carried sum to zero)
            fb = book.flow
            self.impulse += (fb.push + fb.shear + fb.wall + fb.carried).sum(axis=0)
            self.impulse_size += float(sum(np.linalg.norm(x, axis=1).sum() for x in (fb.push, fb.shear, fb.wall)))
            self.flow_courant, self.flow_growth = max(self.flow_courant, fb.courant), max(self.flow_growth, fb.growth)
            self.flow_passes, self.flow_cycles_max = max(self.flow_passes, fb.passes), max(self.flow_cycles_max, fb.cycles)
            self.released_mass += float(fb.released.sum())   # the release from the top cells (Task 7): it leaves the body
            self.released_enthalpy += float(fb.released_H.sum())
            self.momentum_lost += fb.released_p.sum(axis=0)  # the droplets' momentum, at their cells' own velocities
            self.droplet_momentum += fb.released_p
            self.droplet_shear += fb.shear
            if self.spray is not None:                       # the plan's record per column and mode (sub-plan 22 Task 1c),
                self.spray_plan.book(fb)                     # by the cells released, and each flow sub-step's droplets
                if fb.hist_n is not None:                    # binned by size (Task 2)
                    self.hist_n += fb.hist_n
                    self.hist_m += fb.hist_m
        self.iterations.append(res.iterations)
        self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
        self._move_base(book, m_min, g0)
        self.last_shortfall = float(np.maximum(skin_mod.minimum_mass(mat, sk.area, self.settings) - sk.column_mass(),
                                               0.0).sum())

    def _absorb_depth(self, allowance, count=True):
        """Per corner, the depth [m] along the column's direction from the base node to where the 3D temperature falls
        to T_rigid (spec 6.6), from the node's temperature and the 3D gradient along the normal, at most `allowance` [m,
        per corner or one for all]: what the part's draws and hand-downs left of its move. `count`: add the corners the
        allowance held back to absorb_capped. (Against the 3D field traced along the normal, the line method, the
        extrapolation agrees within 0.07 mm on 2 mm elements, measured 2026-10-14.)"""
        mat, mm = self.material, self.moving
        if not mat.melts:
            return np.zeros(self.columns.n)
        T = self.solver.temperature()
        excess = T[mm.corners] - mat.T_rigid
        if not np.any(excess > 0.0):
            return np.zeros(self.columns.n)
        g = self._normal_gradient(T)
        depth = np.where(excess > 0.0, excess / np.maximum(g, GRADIENT_FLOOR), 0.0)
        if count:
            self.absorb_capped += int((depth > allowance).sum())
        return np.minimum(depth, allowance)

    def _normal_gradient(self, T=None):
        """dT/dn [K/m] at every corner along its column's outward direction, from the 3D field's nodal gradient (exact
        for a linear field; positive when the surface is hotter than beneath it)."""
        mm = self.moving
        T = self.solver.temperature() if T is None else T
        mesh = mm.mesh
        touching = np.zeros(mesh.n_nodes, dtype=bool)
        touching[mm.corners] = True
        elements = np.flatnonzero(touching[mesh.tets].any(axis=1))
        return np.einsum("ij,ij->i", movemesh.nodal_gradient(mesh.points, mesh.tets, T, elements)[mm.corners],
                         self.columns.normal)

    def _hand_cap(self):
        """The most a part hands down at each column [kg]: the metal that moves its base by HAND_DOWN_SHARE of the move
        allowance, MOVE_MAX of the mean surface edge (spec 6.6 and 7.1, amended 2026-10-27)."""
        return HAND_DOWN_SHARE * MOVE_MAX * self.columns.edge_length.mean() * self.material.rho * self.columns.area

    def _move_ratio(self, book, m_min):
        """The largest distance a base would move in this part by drawing and handing down, over MOVE_MAX of the mean
        surface edge (absorption takes what is left and the hand-down is capped within it, so neither splits a part).
        The hand-down is found on a copy of the skins."""
        mat, sk = self.material, self.skins
        handed = np.zeros(sk.n)
        if mat.melts:
            saved = sk.state()
            limit = self.solver.temperature()[self.moving.corners] + skin_mod.JUMP_TOLERANCE
            handed, _ = skin_mod.hand_down(sk, m_min, mat.T_rigid, limit, cap=self._hand_cap())
            sk.set_state(saved)
        d = np.abs(book.drawn - handed) / (mat.rho * self.columns.area)
        return float(d.max() / (MOVE_MAX * self.columns.edge_length.mean()))

    def _move_base(self, book, m_min, g0):
        """Draw (already in the skins), absorb and hand down become one corner displacement; the move's actual sweep
        is then reconciled on every column's bottom (module docstring, step 4). A hand-down moves a base by at most
        the allowance, MOVE_MAX of the mean edge (`_hand_cap`, amended 2026-10-27), and absorption takes what the draws
        and hand-downs leave of it. Metal swept out of the 3D model is booked along the profile T(s) = T_node - g0 s
        that the draws used; metal swept into it at the mean enthalpy of what the column handed down, which is only
        metal within the fault tolerance of the node (spec 6.6, amended 2026-10-26); what the tolerance and the cap keep
        back, against a hand-down with neither, is recorded (`held_mass`). The bottom cell is then judged against the
        metal that crossed and, for a hand-down, the cell it uncovered -- the rest of a cell the cap took in part, at its
        own temperature, where the cap ends inside one (`jump_faults`)."""
        mat, sk, mm, col = self.material, self.skins, self.moving, self.columns
        started = time.perf_counter()
        T_node = self.solver.temperature()[mm.corners]
        none = np.full(sk.n, np.nan)
        able = skin_mod.handable(sk, m_min, mat.T_rigid) if mat.melts else np.zeros(sk.n)     # with no limit or cap
        handed, handed_H, coldest, warmest = (
            skin_mod.hand_down(sk, m_min, mat.T_rigid, T_node + skin_mod.JUMP_TOLERANCE, span=True, cap=self._hand_cap())
            if mat.melts else (np.zeros(sk.n), np.zeros(sk.n), none, none))
        uncovered = np.where(handed > 0.0, sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)], np.nan)  # as it left it
        held = np.maximum(able - handed, 0.0)                 # what the tolerance and the cap keep in the skins
        self.held_mass, self.held_depth_max = float(held.sum()), float((held / (mat.rho * col.area)).max())
        if mat.melts:
            self.momentum_lost += sk.momentum_out.sum(axis=0)                     # what the handed-down cells held
        moved = (book.drawn - handed) / (mat.rho * col.area)                      # m inward by drawing and handing down
        depth = self._absorb_depth(np.maximum(MOVE_MAX * col.edge_length.mean() - moved, 0.0))
        inward = book.drawn / mat.rho + depth * col.area - handed / mat.rho          # m3 per corner
        u = -col.normal * (inward / col.area)[:, None] + mm.neutral_smoothing(col.normal, col.area)
        rebuilds = mm.rebuilds
        h_inward = lambda d: skin_mod.profile_enthalpy(mat, T_node, g0, 0.0, d)
        h_outward = np.where(handed > 0.0, handed_H / np.where(handed > 0.0, handed, 1.0), np.nan)
        r = mm.move_with_policy(u, h_inward, h_outward)
        if mm.rebuilds > rebuilds:
            self.wall["rebuild"] += mm.last_rebuild.seconds
        sk.exchange_bottom(mat.rho * r.swept_volume - (book.drawn - handed), r.swept_enthalpy - (book.drawn_H - handed_H))
        self.momentum_lost += sk.momentum_out.sum(axis=0)                         # what a reconciliation took away
        T_bottom = sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)]          # read after the reconciliation
        T_deep = T_node - g0 * (book.drawn / (mat.rho * col.area) + depth)   # the profile at the depth swept
        self.jump_faults += int(jump_faults(T_bottom, T_node, T_deep, (book.drawn > 0.0) | (depth > 0.0),
                                            handed=(handed > 0.0, coldest, warmest, uncovered)).sum())
        sk.regrid()
        self.absorbed_mass += float((mat.rho * depth * col.area).sum())
        self.handed_mass += float(handed.sum())
        self.last_move = r
        self._refresh_geometry()
        self.wall["move"] += time.perf_counter() - started

    # -- what the coupled run reads ---------------------------------------------------------------------------------
    def mass(self, t=None):
        """The 3D model's metal and the skins'. The 3D volume changes only when the corners move, so it is measured once
        per move (`_refresh_geometry`), not at every call: the trajectory asks at every evaluation of its right-hand
        side, and on the 1.2 million elements of the 100 mm sphere's fine band those measurements took 12 % of the time
        of a flight's first six steps (Asha, 2026-10-17)."""
        if self._volume is None:
            self._volume = self.mesh.volume()
        return self.material.rho * self._volume + self.skins.mass()

    def energy(self):
        """Enthalpy held by the 3D model and the skins, and the deferred mismatch in transit to the 3D nodes [J]."""
        return self.solver.energy() + self.skins.energy() + float(self.pending.sum())

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat + removed and released enthalpy) / absorbed heat: zero for exact books."""
        if not self.absorbed_heat:
            return 0.0
        return ((self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy + self.released_enthalpy)
                / self.absorbed_heat)

    def momentum_residual(self):
        """|P - P0 - impulses + momentum lost| over the impulses' total size: zero for exact books (spec 10.2)."""
        if self.impulse_size <= 0.0:
            return 0.0
        miss = self.skins.momentum() - self.momentum0 - self.impulse + self.momentum_lost
        return float(np.linalg.norm(miss) / self.impulse_size)

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / self.mass()))

    def surface_temperature(self):
        """The top cell of every column: what the gas sees (spec 6.3)."""
        return self.skins.T[:, 0].copy()

    def radiated_power(self):
        return float((self.emissivity * SIGMA_SB * (self.skins.T[:, 0] ** 4 - self.T_ambient ** 4) * self.surface.areas).sum())

    def field(self):
        return self.solver.temperature()

    def nose_radius(self):
        return float(self.mesh.params.get("radius_m", 0.05))

    def transverse_radius(self):
        """The top surface's largest distance from the axis through the 3D mesh's mean point along the flight direction
        [m]: the body radius whose quarter caps a droplet (spec 8; the element model's transverse radius)."""
        X = self.surface.points - self.mesh.points.mean(axis=0)
        return float(np.sqrt(np.maximum(np.einsum("ij,ij->i", X, X) - (X @ self.v_hat) ** 2, 0.0)).max())

    def demised(self):
        return self.mass() < self.demise_fraction * self.mass0

    def surface_stats(self):
        T, Tt = self.solver.temperature(), self.surface_temperature()
        centre = self.mesh.points.mean(axis=0)
        i_centre = int(np.argmin(np.linalg.norm(self.mesh.points - centre, axis=1)))
        A = self.surface.areas
        return {"surface_T_max_K": float(Tt.max()), "surface_T_min_K": float(Tt.min()),
                "surface_T_mean_K": float((Tt * A).sum() / A.sum()), "T_stagnation_K": float(Tt[self.i_stag]),
                "T_back_K": float(Tt[self.i_back]), "T_centre_K": float(T[i_centre])}

    def skin_stats(self):
        """The history columns of spec 10.2 that exist at stage 2."""
        sk, mm, mat = self.skins, self.moving, self.material
        th = sk.thickness()
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        flowing = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
        liquid = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
        r = self.last_move
        i, j = self.columns.edges[:, 0], self.columns.edges[:, 1]
        Pi, _, flagged = skin_mod.droplet_momentum_ratio(self.droplet_momentum, self.droplet_shear)
        return {"skin_cells": int(sk.count.sum()), "skin_thickness_max_mm": float(th.max() * 1e3),
                "skin_above_min": int((sk.column_mass() > m_min * (1.0 + 1e-9)).sum()),
                "flowing_depth_mean_mm": float(flowing.mean() * 1e3), "flowing_depth_max_mm": float(flowing.max() * 1e3),
                "liquid_depth_mean_mm": float(liquid.mean() * 1e3), "liquid_depth_max_mm": float(liquid.max() * 1e3),
                "surface_quality_min": mm.surface_quality(), "surface_fold_deg": mm.surface_fold(),
                "zone_ratio_min": float(r.worst_ratio) if r is not None else 1.0, "rebuilds": mm.rebuilds,
                "rebuild_residual_J": mm.rebuild_residual_total,
                "rebuild_residual_last_J": mm.last_rebuild.residual if mm.last_rebuild is not None else 0.0,
                "drawn_mass_kg": self.drawn_mass, "absorbed_mass_kg": self.absorbed_mass, "handed_mass_kg": self.handed_mass,
                "shortfall_kg": self.last_shortfall, "jump_faults": self.jump_faults, "absorb_capped": self.absorb_capped,
                "interface_mismatch_J": self.last_mismatch, "interface_repeats": self.interface_repeats,
                "repeat_columns_over": self.last_trigger["columns_over"],
                "mismatch_worst_ratio": self.last_trigger["worst_ratio"],
                "mismatch_total_ratio": self.last_trigger["total_ratio"],
                "removed_mass_kg": self.removed_mass, "removed_enthalpy_J": self.removed_enthalpy,
                "base_recession_max_mm": float(self._recession().max() * 1e3),
                "wall_skins_s": self.wall["skins"], "wall_3d_s": self.wall["3d"], "wall_move_s": self.wall["move"],
                "wall_rebuild_s": self.wall["rebuild"], "step_parts": self.step_parts,
                # the flow (sub-plan 21 Task 6c): its momentum books, the fastest cell, the top's roughness -- the largest
                # height difference between neighbours over their spacing -- and the flow sub-steps' worst Courant
                # number, Rayleigh-Taylor growth and linearisations in this macro step, and the most of them a heat
                # sub-step took (Task 6d)
                "momentum_residual": self.momentum_residual(),
                "speed_max_ms": float(np.linalg.norm(sk.U, axis=2).max()) if sk.U.size else 0.0,
                "roughness_max": float((np.abs(th[i] - th[j]) / self.columns.edge_length).max()),
                "flow_courant_max": self.flow_courant, "flow_growth_max": self.flow_growth,
                "flow_passes_max": self.flow_passes, "flow_cycles_max": self.flow_cycles_max,
                "wall_flow_s": self.wall["flow"],
                # the droplet momentum sink (Task 7): the mass released from the top cells since time zero, fact 29's
                # droplet-momentum ratio Pi over this macro step and the columns whose own Pi exceeds skin.PI_FLAG
                "released_mass_kg": self.released_mass, "droplet_momentum_ratio": Pi,
                "droplet_momentum_flagged": flagged,
                # what the fault tolerance, and the cap on a part's hand-down, hold in the skins at the last base move
                # (Task 9, spec 6.6, 2026-10-26 and 2026-10-27): unbounded (Asha, 2026-10-26), its mass and its deepest
                # column's depth recorded for sub-plan 22's flow-on flight to report
                "held_mass_kg": self.held_mass, "held_depth_max_mm": self.held_depth_max * 1e3}

    def _recession(self):
        """How far each corner's top has moved inward along its starting direction since time zero [m]."""
        return np.einsum("ij,ij->i", self._top0 - self.surface.points, self._normal0)
