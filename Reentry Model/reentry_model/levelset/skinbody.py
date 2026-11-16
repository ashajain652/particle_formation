"""The level set's skinned body (design sections 10.1, 12.2, 20 and 21; level-set stage 3 Task 9, plan 03b Task 5): one
macro step coupling the skins on the gas-surface markers to the cut-cell 3D solve of stage 2. It wires together the
skins (`skin`), the markers (`markers`), the runoff (`runoff`), the gas loads (`loads`), the spraying (`spraying`), the
thermal solve (`cutthermal`) and the geometry (`boxmesh`, `levelset`, `surfaces`, `division`, `criteria`); the physics
is theirs and the shared package's. It mirrors the skyscrapers' macro step (their `SkinBody.advance` and its parts: the
trial and real passes, the base law, the damped repetition, the deferred mismatch) by amendment, never by import.

The body (design 21.3-21.6). phi, the gas surface, is the evolved level set; beta, the skin base, is derived from it one
skin depth H = `SkinSettings.depth` inside (`surfaces.skin_base`), and the CutThermal solve runs on {beta < 0}; psi, the
molten-rigid surface, is rebuilt from the temperature (`surfaces.molten_rigid_surface`), the volume's on its nodes and
each skin's cell profile between beta and phi. The markers sit on phi, one skin each (a row of `Skins`). The layer
between phi and beta is divided among them (`division.divide`): each skin is the column its marker owns, its
cross-section the division's (`Skins.set_profile`), and every column is exactly full after every acceptance. There are no
skin modes: every skin is H deep, keeps its column full by drawing at its base from the volume and returning there, and
treats its base as a no-slip floor. Mass ownership: the volume owns rho_s V(beta < 0) and the carry of the density
difference (`stepend`), the skins their cells; nothing is counted twice.

The interface (the skyscrapers' per-column law converted per area). Every patch of beta (the solve's boundary pieces)
belongs to the marker that owns its foot on phi: beta's patches are the division's base level, and `_beta_owners` reads
the division's own ownership of it (`Division.base_owner`, by the plain area; the controller's ruling T5-j), so that the
skin above a piece of the base and the group below it are one marker, its own column's base: marker m's group G_m.
GROUP_AREA_RANGE stays as a guard: a group outside it is grown from the patches nearest its base point (`patch_groups`). A column with no base (`stepend.BASELESS_SHARE`) has no group and an insulated bottom (its base
temperature NaN, `Skins.conduct`). A marker's base law a_m + b_m T_b [W, W/K, from its trial pass] reaches each patch p of
its group by the area W_mp of p it holds, q0 = sum_m a_m W_mp / (A(G_m) A_p), q1 likewise with b_m, and its base
temperature T_b is the area mean over its group of the patch temperatures the solve applies its law at -- so the volume
receives exactly sum_m (a_m + b_m T_b,m) dt.

One macro step (`advance`, design 10.1 with 20 and 21): the gas loads per marker from the step's one aero state
(`loads.GasLoads` at the pieces of phi's patches each marker owns, `Markers.owned_patches`, each marker the mean over its
own, `loads.by_patches`; the leeward test, the fitted nose), or a prescribed heating, taken at the step's start and
again on the moved phi before every part, under the ownership that sets the markers' areas there (`_take_loads`; plan 03b
Task 7b, the controller's ruling T7-i); the surface flow's forcing for the runoff (with an aero state or a deceleration),
taken at the step's start and held over it (the controller's ruling T7b-a, `_take_loads`); Girin's plan once
(`spraying.spray_plan`, with a spray model and a state; b is always the skin's flowing depth, the wave's regions on the
division's boundaries at phi; Girin's theory or, with the spray model's `dispersion` "profile", the kinked profile of the
skins' own cell speeds, whose records are 36 columns wide, `record_names`, and whose history ends with the shared
melt_profile.DISPERSION_COLUMNS, `_summary`; plan 03b Task 10); then the step as parts (Ruling 7: none longer than
PART_MAX while any top is at or above T_rigid or about to be, and a part whose move of phi exceeds MOVE_SHARE of the
markers' target spacing -- the part ratio max |dV_i / A_i(0)| / (MOVE_SHARE spacing), A_i(0) the column's cross-section
at phi -- rerun as shorter parts, at most SPLIT_DEPTH_MAX levels and at most PARTS_MAX parts a split, past which the
step stops naming the column; the controller's ruling T7-g). Each part:

1. the trial pass on saved skins (their cells, their flow memory, their profile, the plan's books): every heat sub-step's
   conduction with the trial sensitivity, the columns conducting sideways across the division's shared boundaries
   (`runoff.Sheet.side_faces`, `Skins.conduct(..., side=...)`; plan 03b Task 6), the flow sub-steps with the release
   across the same boundaries (`runoff.Sheet.from_division`, the part's sheet), the feed removal (remove="feed") and
   the regrid, at the linearisation temperature T_lin; its base law (a, b) per marker (none for a column with no base;
   with the sideways heat the derivative for every base temperature moving together, the controller's ruling T4-c,
   the rest of the coupling between neighbours left to the interface mismatch);
2. `CutThermal.step` with each patch's law and the pending load: the heat owed to the volume, delivered at its amount
   over the macro step's time left (at least PENDING_TIME_MIN), spread over the owing markers' patch groups;
3. the real pass from the same saved state with T_b read from the 3D field as above;
4. the mismatch per marker, E_base - (a + b T_b) dt, and the repetition from the saved states by the damped Aitken rule
   (omega0 DAMPED_OMEGA0, omega within DAMPED_OMEGA_RANGE, until the total relative mismatch is at most
   INTERFACE_TOLERANCE, no base moves more than DAMPED_TOL_K, or DAMPED_MAX_PASSES passes);
5. the accepted part's books: what the volume still owes -- each marker's mismatch, and the solve's own residual (the
   heat its law applied less its stored enthalpy's change, the Newton iteration's error), shared by the markers in
   proportion to their base heat -- is the next part's pending load, so the enthalpy closes exactly;
6. the layer (planning ruling 3, `stepend.StepEnd._layer`): each column's volume change dV_i (its cells' change and
   its carry), the move of phi by it, the per-piece correction to the metal's volume, beta derived, `CutThermal.move`,
   the division and the profile, and the fill (the columns with no base from their neighbours, then those with one at
   their bases); the enthalpy the move removed against what the fill took booked as pending at the bases; the
   markers' areas, the top fluxes', given afresh on the moved phi (the markers themselves settle at the step's end),
   and the top fluxes taken again on it before the next part (`_take_loads`). No metal crosses beta within a pass (the
   in-pass draw of stage 3 is retired with its minimum mass).

At the macro step's end (`stepend.StepEnd._end_step`): the last part's carry moves phi, psi, the mesh rebuilt where it
no longer meets its refinement targets (phi carried), the markers settled on phi, the layer divided and filled, and the
thin-metal stop (`stepend.THIN_METAL`, also raised at construction); the carry the step end's own settle and fill leave
waits for the next step's first move (the controller's ruling T5-i); see that module. A step that finds no rigid metal
left anywhere stops with RuntimeError (`stepend.NO_CORE`): at its end, or when a part fails and, the part having
restored the state it started from, psi then finds no rigid metal (the part's error chained as the cause; if psi cannot
be found the part's error stands). No run driver (Ruling 6): tests and scripts drive `advance` at prescribed gas
states.

A step that fails leaves the body in an undefined state (fix round 1 of Task 10): a part that fails in its passes restores
the cells, the profile and the volume's temperatures it started from, but the parts accepted before it and a failed step
end are not undone. A body pickled after a failure is for diagnosis only, not for resuming.

The part count (Ruling 7; fix rounds 1 and 5 of Task 9; `split_part`, `next_hint`): the next step's parts are 1.25 times
the largest move ratio, scaled to the macro step, of the parts that count (the skyscrapers' estimate), but at most 1.25
times the parts that count, at least 1 and at most PARTS_MAX. A part accepted within its allowance counts, whether it
ran as one of the step's parts or was split from one. A part accepted at the deepest split while still over its
allowance does not, and neither does anything the same splitting produced: the part of the step it was split from counts
as the one part it started as, with no move ratio. So the next step may take more than 1.25 times this step's parts only
when parts that met their allowance after splitting ask for it, and a move that splitting does not shrink lets the count
fall back instead of raising it by 1.25 every step to PARTS_MAX, as the skyscrapers' rule did in their stall of
2026-11-08. The look-ahead to T_rigid follows each marker from the nearest marker of the last step's start, so it
survives splits and merges.

Every number here is a numerical setting (named with its source); the physics is the shared package's."""
import math
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

from ..shared import melt_profile, nose
from ..shared import spray as spray_mod
from ..shared.droplets import SKIN_SUMMARY_COLUMNS
from . import loads as loads_mod
from .boxmesh import build
from .cutcell import inside_volumes
from .cutcell import patches as cut_patches
from .cutthermal import CutThermal
from .levelset import pieces
from .markers import Markers
from .runoff import FlowBook, Forcing, Sheet, flow_substeps
from .skin import SkinSettings
from .spraying import spray_plan
from .stepend import (BASELESS_SHARE, FILL_ROUNDS, FILL_RTOL, LEAD, NO_CORE, REGION_SLACK,  # noqa: F401 (re-exported)
                      THIN_METAL, StepEnd)
from .surfaces import skin_base

# the parts of a macro step (Ruling 7, the skyscrapers' skinbody.py rules, matched by amendment)
PART_MAX = 0.0625                 # s: the longest part while any skin's top is at or above T_rigid, or about to be
MOVE_SHARE = 0.25                 # of the markers' target spacing: the most phi moves in one part (Ruling 7, its
                                  # "local spacing" read as the target spacing, fix round 1 of Task 9)
SPLIT_DEPTH_MAX = 6               # a part is rerun shorter at most this many times over (Ruling 7)
MOVE_HINT_MARGIN = 1.25           # the next macro step's parts: this times the parts this one's largest move called for
GROUP_AREA_RANGE = (0.5, 2.0)     # a marker's group of patches of beta has an area within these times its skin's base area
                                  # (the ruling of fix round 2: a whole law on a group 20 times smaller than the base was a
                                  # Robin condition with a Biot number of about 800, and the volume's Newton iteration stalled)
GROUP_CANDIDATES = 64             # patches nearest a base searched first when its group is formed again (then four times more)
PENDING_TIME_MIN = 0.010          # s: the pending heat is delivered at its amount over the macro step's time left, but never
                                  # over less than this (the heat sub-step); what the step's last parts leave stays owed
                                  # into the next step (fix round 3: a whole deferral within a 2 ns part was 1e14 W/m2)
PARTS_MAX = 1000                  # a macro step never starts with more parts (the next step's hint is capped here), and
                                  # a split never asks for more (`split_part` raises instead; plan 03b Task 7b, the
                                  # controller's ruling T7-g: a column's move that no part's length shrinks once asked for
                                  # ceil(4.3e9) parts at every level and the step never ended)
# the interface's damped repetition (the skyscrapers' "damped" rule, their skinbody.py and skin.py)
INTERFACE_TOLERANCE = 1e-3        # of sum |E_top| + |E_base|: the total relative mismatch that ends the repetition
DAMPED_MAX_PASSES = 12            # passes of a part at most (the last pass's mismatch is deferred)
DAMPED_OMEGA0 = 0.5               # the first repetition's relaxation factor
DAMPED_OMEGA_RANGE = (0.02, 1.0)  # Aitken's factor kept within these
DAMPED_TOL_K = 0.05               # K: converged when no base temperature moves more than this between passes
# the markers and the geometry
MARKER_SPACING = 2.0e-3           # m: the markers' target spacing (design 20.13: about 30 000 at 2 mm on the 100 mm sphere)
INITIAL_REFINES = 2               # rebuilds of the first mesh from phi0 alone (psi = phi) at the start
INITIAL_REBUILDS = 4              # then rebuilds of the whole body from phi0 at most, until the step end's test would
                                  # leave the mesh as it is (a cold step's mesh)
INITIAL_SETTLES = 10              # settlings of the new markers at most, until one splits, merges and releases none (on
                                  # the test cube's sharp edges a split and a merge alternate for two more settlings)
BAND_FACTOR = 2.0                 # distance resets are limited to this times (fine_margin + the base cell)
REMOVE_MODES = ("spray", "feed")  # liquid leaves by the spray release only, or as it reaches T_feed (the plates' flow off)

SKIN_SUMMARY_KEYS = tuple(k for _, k in SKIN_SUMMARY_COLUMNS)


@dataclass
class StepBook:
    """One macro step's books [kg, J, kg m/s] and the history's columns (design 12.2). In: `enthalpy_in` the convective
    heat taken at the tops (no mass enters). Out: `mass_sprayed`/`enthalpy_sprayed`/`momentum_sprayed` the release,
    `mass_removed`/`enthalpy_removed` the feed removal, `enthalpy_radiated`. Stored at the step's end: the volume's
    (rho_s (V(beta < 0) + the density carry) and its enthalpy) and the skins'; `pending_J` the heat still owed to the
    volume (deferred to the next part). The residuals: (mass now - mass then + out) / mass then, and (stored + pending
    now - stored + pending then - in + out) / (|stored then| + |in|). `parts` the parts the step took, `parts_next` the
    next step's hint, and `parts_skyscrapers` the hint the skyscrapers' rule would have given (every accepted part's move
    ratio, unlimited), `split_depth_max` the deepest split level a part reached, `parts_over` the parts accepted at the
    deepest split while still over their move allowance. `snap_fallback` the markers the snap's last round
    (`stepend.StepEnd._snap`) brought onto phi by Newton steps from farther than `markers.NEWTON_REACH`, and
    `snap_slide_max` [m] the largest sideways slide (along the surface) of that round's Newton steps.
    `pieces_used_up` the pieces of phi the step's moves and its remesh parted from the body holding no metal (no
    marker on them and none of the volume's metal), lifted out (`stepend.StepEnd._lift_specks`).
    `volume_repartitioned_share` the layer's volume passed between columns at the step end, summed over its
    re-partitions, over the layer's volume after it (an area's share, `area_repartitioned_share`, until plan 03b Task 6
    moved the re-partition to the division's volumes; renamed in Task 7, and a book pickled under the old name loads
    under the new one, `__setstate__`), and
    `stretch_min`, `stretch_max` the least and largest stretch of an old column in them (`repartition`, Task 9b;
    `LevelSetSkinBody.last_settle` keeps each re-partition).
    `group_area_ratio_min`, `group_area_ratio_max` the least and largest area of a marker's group of patches of beta over
    its skin's base area, over the step's parts and the columns with a base (`patch_groups`), and
    `group_area_ratio_trimmed_max` the largest with the group's largest patch left out (a group formed again stops at
    the first patch that takes it past the base area, so only that one patch may carry it past the range; fix round
    5); `groups_guarded_max` the most groups that guard (GROUP_AREA_RANGE) replaced in any part (Task 11's fix, item 1).
    `pending_flux_max` the largest flux [W/m2] the pending heat put on a patch of beta in any part, `surface_flux_max`
    the largest convective flux on a top. `mass_drawn`, `mass_handed` [kg] the fills' exchange at the bases, each way;
    `momentum_to_volume` the momentum the held volume took over the step (the cells the fills returned;
    `LevelSetSkinBody.momentum_to_volume` keeps the total since the start), `momentum_to_wall` what the skins' flow gave
    the wall over the step (every skin's floor is no-slip) and `momentum_removed` what the feed removal took.
    `island_nodes` the volume's nodes the moves made apart from the old metal, which took the skins' temperatures there
    (stage 3 Task 11, item 3). The layer (plan 03b Task 5): `layer_factor_max` the largest departure from one of a
    piece's closure factor over the step's divisions (`division.Division.factor`; design 21.8 item 5, DIVISION_TOL the
    tests' bound), `fill_miss_max` the largest share of its volume by which a column's cells missed it after a fill,
    `move_miss_max` the largest miss of a group's swept volume against what it was asked (`levelset.Moved.miss`; what a
    group cannot sweep the per-piece correction makes up, the controller's ruling T2-e), `n_baseless` the columns with no
    base at the step's end, `density_carry_m3` [m3] the carry left for the next step's first move (planning ruling 4:
    the step end's own settle's and fill's, in the volume's mass, by which phi misses its metal's true volume until that
    move; the controller's ruling T5-i) and `density_carry_share_max` the largest share of its column's volume a row's
    carry is then. `move_corrections` every volume correction's miss before it, in the order made: one entry per move
    of phi and for the remesh, {"kind": "part", "carry" (the step end's move by the
    last part's carry) or "remesh"; "worst": the largest miss over its pieces as a share of the piece's target; "body":
    the body's, the largest piece's}; their largest is `volume_correction_max`. The heat deferred to the volume over the
    step, by source, each the sum of the absolute values of its bookings [J]: `residual_newton_J` the volume solve's own
    residual (the heat its law applied less its stored enthalpy's change), `residual_move_J` the moves' and fills' (the
    volume's move against what the skins took and gave), `residual_remesh_J` the remesh's; the interface mismatch is
    `interface_mismatch_J`. History columns: the markers and their counters (cumulative), `interface_mismatch_J` (the
    last part's sum of |mismatch|), `volume_correction_max`, `n_pieces`, `frontal_area_m2` (phi's area projected on the
    flight direction) and `summary` (SKIN_SUMMARY_COLUMNS' keys and the release angle's)."""
    time: float
    dt: float
    parts: int
    parts_next: int
    parts_skyscrapers: int
    split_depth_max: int
    parts_over: int
    snap_fallback: int
    snap_slide_max: float
    pieces_used_up: int
    volume_repartitioned_share: float
    stretch_min: float
    stretch_max: float
    group_area_ratio_min: float
    group_area_ratio_max: float
    group_area_ratio_trimmed_max: float
    pending_flux_max: float
    surface_flux_max: float
    passes: int
    enthalpy_in: float
    mass_sprayed: float
    enthalpy_sprayed: float
    mass_removed: float
    enthalpy_removed: float
    enthalpy_radiated: float
    mass_volume: float
    mass_skins: float
    enthalpy_volume: float
    enthalpy_skins: float
    pending_J: float
    mass_residual: float
    enthalpy_residual: float
    mass_drawn: float
    mass_handed: float
    momentum_sprayed: np.ndarray
    momentum_removed: np.ndarray
    momentum_to_volume: np.ndarray
    momentum_to_wall: np.ndarray
    residual_newton_J: float
    residual_move_J: float
    residual_remesh_J: float
    n_markers: int
    marker_splits: int
    marker_merges: int
    markers_trapped: int
    interface_mismatch_J: float
    volume_correction_max: float
    n_pieces: int
    frontal_area_m2: float
    layer_factor_max: float
    fill_miss_max: float
    move_miss_max: float
    n_baseless: int
    density_carry_m3: float = 0.0
    density_carry_share_max: float = 0.0
    island_nodes: int = 0
    groups_guarded_max: int = 0
    move_corrections: list = field(default_factory=list)
    summary: dict = field(default_factory=dict)

    def __setstate__(self, state):
        """Unpickling: a book from before plan 03b Task 7 holds the re-partitioned share under its old name
        (`area_repartitioned_share`), taken under the new one."""
        state = dict(state)
        if "area_repartitioned_share" in state:
            state.setdefault("volume_repartitioned_share", state.pop("area_repartitioned_share"))
        self.__dict__.update(state)

    def row(self):
        """The step as one history row (a flat dict; vectors as their three components)."""
        out = {}
        for key, value in self.__dict__.items():
            if key == "summary":
                out.update(value)
            elif isinstance(value, np.ndarray):
                out.update({"{}_{}".format(key, c): float(v) for c, v in zip("xyz", value)})
            else:
                out[key] = value
        return out


class _Tally:
    """The running sums of one macro step."""

    def __init__(self, n):
        self.E_top = self.E_rad = 0.0
        self.sprayed = self.sprayed_H = self.removed = self.removed_H = self.drawn = self.handed = 0.0
        self.sprayed_p, self.wall, self.removed_p, self.to_volume = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)
        self.parts, self.passes, self.mismatch, self.correction = 0, 0, 0.0, 0.0
        self.depth = self.over = 0
        self.snap_fallback, self.snap_slide = 0, 0.0
        self.ratio_all = 0.0
        self.group_min, self.group_max, self.group_trimmed = np.inf, 0.0, 0.0
        self.newton = self.move = self.remesh = 0.0
        self.elapsed, self.pending_flux = 0.0, 0.0
        self.islands = 0                               # nodes the moves' islands took from the skins' temperatures
        self.guarded = 0                               # the most groups the guard replaced in a part
        self.corrections = []                          # every volume correction's miss before it, labelled
        self.layer_factor = 0.0                        # the divisions' largest |factor - 1|
        self.fill_miss = 0.0                           # the fills' largest miss, a share of the column's volume
        self.move_miss = 0.0                           # the moves' largest group miss
        self.baseless = 0                              # the columns with no base at the last fill
        self.pieces_used = 0                           # pieces of phi the moves left with no metal


@dataclass
class _Pass:
    """What one real pass did, per skin row: the heat books [J], the feed removal and the flow's book."""
    E_top: np.ndarray
    E_rad: np.ndarray
    E_base: np.ndarray
    removed: np.ndarray
    removed_H: np.ndarray
    removed_p: np.ndarray
    flow: FlowBook


@dataclass
class _StillGas:
    """The surface flow's two fields the runoff reads, for a deceleration with no gas state: no pressure gradient
    [Pa/m] and no shear [Pa] per marker."""
    dp_ds: np.ndarray
    tau: np.ndarray


@dataclass(frozen=True)
class PartCount:
    """What one part, with the shorter parts it was split into, tells the next step's part count (`next_hint`):
    `counted` the accepted parts that count, `rate` [-] the largest move ratio among them scaled to the macro step
    (ratio x dt_macro / dt, the parts the macro step would need for that move), `failed` whether some part was accepted
    at the deepest split while still over its allowance (then nothing below this part counts), `depth` the deepest
    split level reached and `over` the parts accepted over their allowance."""
    counted: int
    rate: float
    failed: bool
    depth: int
    over: int


def split_part(run, dt, dt_macro, depth, worst=None):
    """Ruling 7's splitting: `run(dt)` runs a part of length dt [s] and returns (its move ratio, accept, undo), the
    largest move of phi over MOVE_SHARE of the target spacing; a ratio over one below SPLIT_DEPTH_MAX levels undoes the
    part and runs it again as ceil(ratio) parts one level down, otherwise the part is accepted. Returns the `PartCount`.
    A part split into shorter parts of which one was accepted at the deepest level still over its allowance counts as
    failed, with nothing below it counted (fix round 5 of Task 9: a move that does not shrink with the part's length is
    not helped by splitting, and the parts and move ratios the splitting produced around it -- sibling parts of a few
    nanoseconds -- must not set the next step's count). A split that would ask for more than PARTS_MAX parts undoes the
    part and raises RuntimeError instead, naming the column that sets its ratio, `worst()` (a description, called before
    the undo; the controller's ruling T7-g)."""
    ratio, accept, undo = run(dt)
    if ratio > 1.0 and depth < SPLIT_DEPTH_MAX:
        k = int(math.ceil(ratio))
        if k > PARTS_MAX:
            named = worst() if worst is not None else "ratio {:.6e}".format(ratio)
            undo()
            raise RuntimeError("a part of {:.6e} s would be split into {} parts, more than PARTS_MAX ({}): its move of "
                               "phi does not fit the allowance; the worst column, {}".format(dt, k, PARTS_MAX, named))
        undo()
        kids = [split_part(run, dt / k, dt_macro, depth + 1, worst) for _ in range(k)]
        deepest, over = max(c.depth for c in kids), sum(c.over for c in kids)
        if any(c.failed for c in kids):
            return PartCount(0, 0.0, True, deepest, over)
        return PartCount(sum(c.counted for c in kids), max(c.rate for c in kids), False, deepest, over)
    accept()
    if ratio <= 1.0:
        return PartCount(1, ratio * dt_macro / dt, False, depth, 0)
    return PartCount(0, 0.0, True, depth, 1)


def next_hint(counts):
    """The next macro step's part count from its parts' `PartCount`s, one per part the step started with (module
    docstring): MOVE_HINT_MARGIN times the largest scaled move ratio of the counted parts (the skyscrapers' estimate),
    but at most MOVE_HINT_MARGIN times the counted parts, at least 1 and at most PARTS_MAX. A part whose splitting failed
    counts as the one part it started as, with no move ratio."""
    counted = sum(1 if c.failed else c.counted for c in counts)
    rate = max((c.rate for c in counts if not c.failed), default=0.0)
    estimate = int(math.ceil(MOVE_HINT_MARGIN * rate))
    used = int(math.ceil(MOVE_HINT_MARGIN * counted))
    return min(max(1, min(estimate, used)), PARTS_MAX)


def patch_groups(centroid, area, base, A_base, owner=None, report=False, insulated=None):
    """Each marker's group of patches of beta and how its base law is spread over them (fix round 2 of Task 9): (S, A_group).
    S (n, P), sparse and row-stochastic, holds marker m's share of each patch p of its group G_m, the area of p it holds
    over A(G_m); A_group (n,) [m2] the groups' areas. `owner` gives each row the patches beneath its part of phi
    (`LevelSetSkinBody._beta_owners`: a sparse (n, P) of the area of each patch each row holds; or (P,) one row per
    whole patch); without it every patch goes to the marker whose base point (`base`, (n, 3)) is nearest its `centroid`
    (P, 3). A group whose area is outside GROUP_AREA_RANGE times its marker's skin base area `A_base` (n,) [m2] (the
    guard) is replaced by the patches nearest that marker's base, whole, taken in order of distance until their area
    first reaches A_base (shared with other groups where they overlap). The rows `insulated` (n,) bool (columns with no
    base, plan 03b Task 5) take no group and are never formed again. With `report`, the number of groups the guard
    replaced is returned third. A patch a shrunk group gave up and no other group holds takes no law. The law stays exact
    whatever the groups: with q0 = S^T a / A_p, q1 = S^T b / A_p, sum q A over the patches is sum a and sum b over the
    rows with groups, and a marker's base temperature S T_p is the area mean over its group. A zero A_base leaves the
    nearest-base groups (and a marker with none the patch nearest its base)."""
    n, npatch = len(base), len(area)
    A_base = np.broadcast_to(np.asarray(A_base, dtype=float), (n,))
    off = np.zeros(n, dtype=bool) if insulated is None else np.asarray(insulated, dtype=bool)
    if owner is None:
        _, owner = cKDTree(base).query(centroid)
    if sp.issparse(owner):
        G = sp.csr_matrix(owner, shape=(n, npatch))
    else:
        owner = np.asarray(owner, dtype=np.int64)
        G = sp.csr_matrix((area, (owner, np.arange(npatch))), shape=(n, npatch))
    A_near = np.asarray(G.sum(axis=1)).ravel()
    lo, hi = GROUP_AREA_RANGE
    ratio = np.divide(A_near, A_base, out=np.full(n, np.inf), where=A_base > 0.0)
    redo = np.flatnonzero((((ratio < lo) | (ratio > hi)) & (A_base > 0.0) | (A_near == 0.0)) & ~off)
    keep = np.where(off, 0.0, 1.0)
    keep[redo] = 0.0
    G = (sp.diags(keep) @ G).tocoo()
    rows, cols, vals = [G.row], [G.col], [G.data]
    if len(redo):
        tree = cKDTree(centroid)
        k = min(GROUP_CANDIDATES, npatch)
        while True:
            _, idx = tree.query(base[redo], k=k)
            idx = idx.reshape(len(redo), k)
            cum = np.cumsum(area[idx], axis=1)
            want = np.where(A_base[redo] > 0.0, A_base[redo], area[idx[:, 0]])
            enough = cum[:, -1] >= want * (1.0 - 1e-12)
            if enough.all() or k == npatch:
                break
            k = min(4 * k, npatch)
        take = np.argmax(cum >= want[:, None] * (1.0 - 1e-12), axis=1)
        take = np.where(enough, take, k - 1)
        sel = np.arange(k)[None, :] <= take[:, None]
        rows.append(np.broadcast_to(redo[:, None], idx.shape)[sel])
        cols.append(idx[sel])
        vals.append(area[idx[sel]])
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    live = vals > 0.0
    rows, cols, vals = rows[live], cols[live], vals[live]
    A_group = np.bincount(rows, vals, minlength=n)
    S = sp.csr_matrix((vals / A_group[rows], (rows, cols)), shape=(n, npatch))
    return (S, A_group, len(redo)) if report else (S, A_group)


_PLAN_BOOKS = ("mass", "mass_slurry", "number", "number_slurry", "enthalpy", "enthalpy_slurry", "solid",
               "solid_slurry", "momentum", "momentum_slurry", "hist_n", "hist_m")


class LevelSetSkinBody(StepEnd):
    """The skinned level-set body (module docstring). `material` the shared Material (or a test material with its
    interface), `base` the mesh's BoxBase, `phi0` a function of points (k, 3) [m] giving the initial gas surface's level
    set (negative inside the metal; its zero contour is the body), `T0` [K] the initial temperature, `settings` the
    SkinSettings (its `depth` the skin depth H). `loads` a `loads.GasLoads` (needed with an aero state), `spray_model`
    the shared SprayModel (Girin's release, with an aero state), `leeward` the leeward test ("shadow" or "centroid"),
    `emissivity` of the tops, `fine_margin`, `fine_size`, `cell_floor` [m] the refinement (design 20.4-20.5), `remove`
    "spray" or "feed". Also: `spacing` [m] the markers' target spacing, `v_hat` the flight direction in the body's frame,
    `T_ambient` [K] for the tops' radiation, `slurry_spray` (`--slurry-spray`), `linear_solver` the volume's ("amg" or
    "direct").

    The mesh is built from phi0 and refined from it (INITIAL_REFINES, then whole rebuilds until the step end's test would
    keep the mesh, INITIAL_REBUILDS); phi is phi0 on its nodes, reset to a distance keeping its contour; beta is derived;
    the markers are seeded one per voxel of the target spacing on phi's patches and settled; the layer is divided and
    every skin filled at T0 exactly to its column. Metal thinner than twice the skin depth raises RuntimeError
    (`stepend.THIN_METAL`)."""

    def __init__(self, material, base, phi0, T0, settings=None, loads=None, spray_model=None, leeward="shadow",
                 emissivity=0.0, fine_margin=4e-3, fine_size=0.7e-3, cell_floor=0.25e-3, remove="spray",
                 spacing=MARKER_SPACING, v_hat=(1.0, 0.0, 0.0), T_ambient=0.0, slurry_spray=False,
                 linear_solver="amg"):
        if remove not in REMOVE_MODES:
            raise ValueError("remove must be one of {}, got {!r}".format(REMOVE_MODES, remove))
        if leeward not in loads_mod.LEEWARD_NAMES:
            raise ValueError("leeward must be one of {}, got {!r}".format(loads_mod.LEEWARD_NAMES, leeward))
        if not callable(phi0):
            raise ValueError("phi0 must be a function of the points (k, 3) giving the initial gas surface's level set")
        self.material, self.base = material, base
        self.settings = settings if settings is not None else SkinSettings()
        self.loads, self.spray_model, self.leeward = loads, spray_model, leeward
        self.emissivity, self.T_ambient = float(emissivity), float(T_ambient)
        self.fine_margin, self.fine_size, self.cell_floor = float(fine_margin), float(fine_size), float(cell_floor)
        self.remove, self.spacing, self.slurry_spray = remove, float(spacing), bool(slurry_spray)
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.rho = float(material.rho)
        self.band = BAND_FACTOR * (self.fine_margin + base.h0)
        self.time = 0.0
        self.pending_x, self.pending_J = np.zeros((0, 3)), np.zeros(0)
        self.momentum_to_volume = np.zeros(3)
        self._hint, self._top_start, self._dt_last = 1, None, None      # _top_start: (marker positions, top T) at the
                                                                        # last step's start
        self.n_released = 0.0
        self.source_rows = []
        self.last_step = None
        self.plan = None
        self.division = None
        self.last_settle = []

        # the mesh, refined from phi0 (all rigid at the start), then the body built on it until it keeps its mesh
        mesh = build(base)
        phi = self._reset(mesh, np.asarray(phi0(mesh.points), float))
        for _ in range(INITIAL_REFINES):
            psi = phi if T0 <= material.T_rigid else np.full(len(phi), np.inf)
            mesh = build(base, self._target(mesh, phi, psi))
            phi = self._reset(mesh, np.asarray(phi0(mesh.points), float))
        self._build_on(mesh, phi, T0, linear_solver)
        for _ in range(INITIAL_REBUILDS):
            target = self._target(self.mesh, self.phi, self.psi)
            if not self._needs_rebuild(target, self._target(self.mesh, self.phi, self.psi, False,
                                                            LEAD - REGION_SLACK * self.fine_size)):
                break
            mesh = build(base, target)
            self._build_on(mesh, self._reset(mesh, np.asarray(phi0(mesh.points), float)), T0, linear_solver)

    def _build_on(self, mesh, phi, T0, linear_solver):
        """The body on `mesh` with the gas surface `phi` (class docstring): beta, the markers seeded and settled, the
        division (the thin-metal stop), the volume at T0 and every skin filled at T0 to its column; psi."""
        material, base = self.material, self.base
        self.mesh, self.phi = mesh, phi
        self.beta = skin_base(mesh, phi, self.settings.depth, band=self.band)
        P = cut_patches(mesh, phi)
        key = np.floor((P.centroid - np.asarray(base.lo)) / self.spacing).astype(np.int64)
        order = np.lexsort((-P.area, key[:, 2], key[:, 1], key[:, 0]))
        k = key[order]
        first = np.r_[True, np.any(k[1:] != k[:-1], axis=1)]
        seeds = np.sort(order[first])
        mk = Markers.create(mesh, phi, seeds, material=material, settings=self.settings, spacing=self.spacing)
        self.markers, self.skins = mk, mk.skins
        mk.tessellate(mesh, phi)
        while mk.merge(self.spacing):
            mk.tessellate(mesh, phi)
        self._sync()
        self.carry = np.zeros(mk.skins.n)
        tally = _Tally(0)
        self._divide(tally)
        self._check_thin()
        self.solver = CutThermal(material, emissivity=0.0, linear_solver=linear_solver)
        self.solver.setup(mesh, self.beta)
        self.solver.set_temperature(T0)
        self._fill_at(T0)
        for _ in range(INITIAL_SETTLES):                  # (the skins are filled afresh after: no re-partition)
            changes = mk.splits + mk.merges + mk.trapped
            self._settle(tally, repartition=False)
            if mk.splits + mk.merges + mk.trapped == changes:
                break
        self._divide(tally)
        self._check_thin()
        self._fill_at(T0)
        self.psi = self._psi(phi)

    def _fill_at(self, T0):
        """Every skin anew: its column's volume of metal at T0 (at T0's density), regridded, and no carry."""
        sk = self.skins
        sk.fill(sk.density(np.full(sk.n, float(T0))) * sk.column_volume(), T0)
        sk.regrid()
        self.carry = np.zeros(sk.n)

    def __setstate__(self, state):
        """Unpickling: a body from before plan 03b Task 5 holds `mode_settings` (the retired modes', dropped) and has no
        division and no carry: `advance` then adopts it (`_adopt`). So too a body whose division predates its base level's
        ownership (Task 5's fix round 1) or its columns' first moments (Task 6's ruling T6-b), which `_adopt` divides
        afresh."""
        state = {k: v for k, v in state.items() if k != "mode_settings"}
        self.__dict__.update(state)
        if "division" not in state or getattr(state["division"], "base_owner", None) is None or \
                getattr(state["division"], "moment", None) is None:
            self.division = None
        if "carry" not in state:
            self.carry = np.zeros(self.skins.n)

    def _adopt(self):
        """A body from before plan 03b Task 5 on the one skin depth: beta derived from its phi and the volume moved to
        it, the layer divided and every column filled at its base (the enthalpy the move removed against what the fill
        took owed to the volume)."""
        tally = _Tally(self.skins.n)
        beta = skin_base(self.mesh, self.phi, self.settings.depth, band=self.band)
        removed = self.solver.move(beta)
        self.beta = beta
        self._divide(tally)
        owed = removed - self._fill(tally)
        self._book_move(self._bases(), self.skins.area, owed)

    def _no_rigid_metal(self):
        """Whether psi, found now from the body's fields, finds no rigid metal (`stepend.StepEnd._check_core`); False if
        psi cannot be found (the caller's own error then stands)."""
        try:
            return not np.any(self._psi(self.phi) < 0.0)
        except Exception:
            return False

    # -- small helpers ------------------------------------------------------------------------------------------------
    def volume_energy(self):
        return self.solver.energy()

    def books(self):
        """The body's stored mass [kg] and enthalpy [J] now: the volume's (rho_s (V(beta) + the carry): the density
        carry is volume metal phi has yet to take, `stepend`), the skins' and the pending heat."""
        sk = self.skins
        return {"mass_volume": self.rho * (float(inside_volumes(self.mesh, self.beta).sum()) + float(self.carry.sum())),
                "mass_skins": float(sk.m.sum()), "enthalpy_volume": self.solver.energy(),
                "enthalpy_skins": float(sk.H.sum()), "pending_J": float(self.pending_J.sum())}

    def add_pending(self, x, J):
        """Book J [J] of heat owed to the volume at the point x [m] (delivered at the next part, at the nearest patch)."""
        self.pending_x = np.vstack([self.pending_x, np.atleast_2d(np.asarray(x, float))])
        self.pending_J = np.r_[self.pending_J, np.atleast_1d(np.asarray(J, float))]

    def _carry_share(self):
        """The largest share of its column's volume a row's carry is (a row with a carry and no column, inf)."""
        V, c = self.skins.column_volume(), np.abs(self.carry)
        share = np.divide(c, V, out=np.where(c > 0.0, np.inf, 0.0), where=V > 0.0)
        return float(share.max(initial=0.0))

    def _top_T(self):
        T = self.skins.T()
        return np.where(self.skins.valid[:, 0], T[:, 0], np.nan)

    def mass_centroid(self):
        """The metal's mass centre [m] (3,): the volume's cells at their centroids and the skins' cells at mid-depth."""
        s = self.solver
        mv = self.rho * s.Vin
        cen = self.mesh.points[self.mesh.tets[s.cells]].mean(axis=1)
        sk = self.skins
        dz = sk.dz()
        depth = np.cumsum(dz, axis=1) - 0.5 * dz
        r = self._row_marker()
        at = self.markers.x[r][:, None, :] - depth[..., None] * self.markers.normal[r][:, None, :]
        M = mv.sum() + sk.m.sum()
        return ((mv[:, None] * cen).sum(axis=0) + (sk.m[..., None] * at).sum(axis=(0, 1))) / M

    # -- the interface ---------------------------------------------------------------------------------------------------
    def _beta_owners(self, P):
        """The area of each patch of beta each skin row owns, (rows, P) csr [m2] (module docstring): the division's
        ownership of its base level (`Division.base_owner`, by the plain area; the controller's ruling T5-j), whose
        patches are beta's (beta keeps phi + H's values in the cells it cuts, `surfaces.skin_base`), each patch of P
        matched by its cell and its owners' areas scaled to P's area of it. The columns sum to the patches' areas to
        round-off. RuntimeError if the division is not of phi as it is (a patch of beta in a cell its base level does
        not cut): the division is rebuilt whenever phi moves, before the next part's interface."""
        mk, sk, d = self.markers, self.skins, self.division
        npatch = len(P.area)
        if not npatch or not mk.n:
            return sp.csr_matrix((sk.n, npatch))
        order = np.argsort(d.base_cell)
        k = np.searchsorted(d.base_cell[order], P.cell)
        found = k < len(order)
        found[found] = d.base_cell[order[k[found]]] == P.cell[found]
        if not found.all():
            raise RuntimeError("{} of {} patches of beta lie in cells the division's base level does not cut: the "
                               "division is not of phi as it is".format(int((~found).sum()), npatch))
        at = order[k]
        W = d.base_owner[:, at].tocoo()
        scale = np.divide(P.area, d.base_patch_area[at], out=np.ones(npatch), where=d.base_patch_area[at] > 0.0)
        return sp.csr_matrix((W.data * scale[W.col], (mk.skin_index[W.row], W.col)), shape=(sk.n, npatch))

    def _interface(self):
        """(S (n, P) sparse, row-stochastic: row m's share A_p / A(G_m) of each patch of its group; A(G) (n,) [m2], zero
        for a column with no base; the base points (n, 3)) on the solve's present patches; the groups the guard replaced
        in `last_guarded`."""
        P = self.solver.P
        base = self._bases()
        S, A_group, self.last_guarded = patch_groups(P.centroid, P.area, base, self.skins.base_area(),
                                                     owner=self._beta_owners(P), report=True,
                                                     insulated=self._baseless())
        return S, A_group, base

    def _pending_per_patch(self, S, base, has):
        """The pending heat [J] per patch of beta: each pending point's heat goes to the row with a group (`has`) whose
        base point is nearest it and is spread over that row's group as its base law is (`patch_groups`, S), so that the
        sum is exact. (Fix round 2: deposited whole on the single patch nearest each point, a few joules made -6e12 W/m2
        on a sliver over a 0.5 ms part and the volume went to -900 K inside the metal.)"""
        out = np.zeros(S.shape[1])
        rows = np.flatnonzero(has)
        if len(self.pending_J) and len(rows):
            _, k = cKDTree(base[rows]).query(self.pending_x)
            out = S.T @ np.bincount(rows[k], self.pending_J, minlength=S.shape[0])
        return out

    # -- one pass of the skins -----------------------------------------------------------------------------------------
    def _sheet(self):
        """The part's sheet: the skins on the division's shared boundaries, depth-resolved, in the skins' rows
        (`runoff.Sheet.from_division`; the division is of phi and the markers as the part finds them)."""
        return Sheet.from_division(self.division, self.markers, self.skins)

    def _pass(self, dt, T_b, ctx, trial, sheet):
        """The skins over one part dt [s] with the base held at T_b [K] (n,) (NaN: an insulated base), on the part's
        `sheet` (`_sheet`): per heat sub-step, the conduction with the columns conducting sideways across its
        boundaries (`Sheet.side_faces`), the flow sub-steps across them with the release, the feed removal and the
        regrid. A trial pass returns its base law (a, b); a real pass its books (_Pass)."""
        sk, s = self.skins, self.settings
        side = sheet.side_faces()
        n = sk.n
        N = max(int(math.ceil(dt / s.substep - 1e-9)), 1)
        h = dt / N
        z = np.zeros(n)
        E_top, E_rad, E_base = z.copy(), z.copy(), z.copy()
        removed, removed_H, removed_p = z.copy(), z.copy(), np.zeros((n, 3))
        fbook = FlowBook.zero(n)
        release = ctx["plan"].release if ctx["plan"] is not None else None
        for _ in range(N):
            hb = sk.conduct(h, ctx["q"], self.T_ambient, self.emissivity, T_b, trial=trial, side=side)
            E_top += hb.E_top
            E_rad += hb.E_rad
            E_base += hb.E_base
            if ctx["forcing"] is not None:
                fbook.add(flow_substeps(sk, sheet, ctx["forcing"], h, release=release))
            if self.remove == "feed":
                mass, H = sk.remove_feed()
                removed += mass
                removed_H += H
                removed_p += sk.momentum_out
            sk.regrid()
        if trial:
            return sk.base_law(dt)
        return _Pass(E_top, E_rad, E_base, removed, removed_H, removed_p, fbook)

    def _save(self, ctx):
        sk, plan = self.skins, ctx["plan"]
        books = None if plan is None else {k: np.array(getattr(plan, k), copy=True) for k in _PLAN_BOOKS}
        profile = None if sk.profile is None else (sk.profile_levels.copy(), sk.profile.copy())
        return self.solver.T.copy(), sk.state(), sk.flow_memory.copy(), profile, books

    def _restore(self, saved, ctx, solver=True):
        """Undo a pass: the cells, the flow memory, the profile and the plan's books, and with `solver` the volume's
        temperatures."""
        T, state, memory, profile, books = saved
        if solver:
            self.solver.T = T.copy()
        self.skins.set_state(state)
        self.skins.flow_memory = memory.copy()
        if profile is not None:
            self.skins.profile_levels, self.skins.profile = profile[0].copy(), profile[1].copy()
        if books is not None:
            for k, v in books.items():
                setattr(ctx["plan"], k, np.array(v, copy=True))

    # -- one part --------------------------------------------------------------------------------------------------------
    def _part(self, dt, dt_macro, depth, ctx, tally):
        """One part of a macro step (module docstring, steps 1-6), a part whose move of phi exceeds MOVE_SHARE of the
        markers' target spacing run again as shorter parts (Ruling 7, `split_part`; past PARTS_MAX of them it stops,
        naming `_worst_column`). Returns its `PartCount`."""
        return split_part(lambda h: self._part_once(h, dt_macro, ctx, tally), dt, dt_macro, depth,
                          worst=self._worst_column)

    def _moves(self):
        """Each row's move of phi [m] (n,) as the part's passes left it: its volume change since its column was last
        filled, its cells' and its carry, over its column's cross-section at phi, A(0) (a column with none, zero)."""
        sk = self.skins
        dC = sk.cells_volume() - sk.column_volume()
        A0 = sk.profile[:, 0] if sk.profile is not None else sk.area
        return np.divide(np.abs(dC + self.carry), A0, out=np.zeros(sk.n), where=A0 > 0.0)

    def _worst_column(self):
        """The column whose move sets a part's ratio (`_moves`): its marker, position, move and ratio, for a stop."""
        move, allowance = self._moves(), MOVE_SHARE * self.spacing
        r = int(np.argmax(move))
        m = int(self._row_marker()[r])
        return ("marker {} at ({:.3f}, {:.3f}, {:.3f}) mm: its move of phi {:.6e} m against the allowance {:.6e} m "
                "(MOVE_SHARE of the spacing), ratio {:.6e}".format(m, *(self.markers.x[m] * 1e3), move[r], allowance,
                                                                   move[r] / allowance))

    def _part_once(self, dt, dt_macro, ctx, tally):
        """Run one part of length dt [s] (steps 1-4) and return (its move ratio, accept, undo): `accept()` books it and
        settles the layer (steps 5-6), `undo()` restores the state it started from. The top fluxes are first taken again
        if phi has moved since they were last taken (`_take_loads`)."""
        sk, solver = self.skins, self.solver
        if ctx["surface"] is not self.phi:
            self._take_loads(ctx)
        S, A_group, base = self._interface()
        has = A_group > 0.0
        tally.guarded = max(tally.guarded, self.last_guarded)
        A_b = sk.base_area()
        ok = has & (A_b > 0.0)
        if ok.any():
            group_ratio = A_group[ok] / A_b[ok]
            tally.group_min = min(tally.group_min, float(group_ratio.min()))
            tally.group_max = max(tally.group_max, float(group_ratio.max()))
            largest = S.max(axis=1).toarray().ravel() * A_group            # each group's largest patch [m2]
            tally.group_trimmed = max(tally.group_trimmed, float(((A_group - largest)[ok] / A_b[ok]).max()))
        saved = self._save(ctx)
        sheet = self._sheet()
        try:
            T_b0 = np.where(has, S @ solver.patch_T(), np.nan)
            A_p = solver.P.area
            share = min(1.0, dt / max(dt_macro - tally.elapsed, PENDING_TIME_MIN))   # of the pending heat, this part's
            q_pend = share * self._pending_per_patch(S, base, has) / (A_p * dt)         # (at its rate over the rest)
            tally.pending_flux = max(tally.pending_flux, float(np.abs(q_pend).max(initial=0.0)))
            E_vol0 = solver.energy()
            T_lin = T_b0
            omega, r_prev = DAMPED_OMEGA0, None
            for attempt in range(DAMPED_MAX_PASSES):
                a, b = self._pass(dt, T_lin, ctx, True, sheet)
                a, b = np.where(has, a, 0.0), np.where(has, b, 0.0)          # a column with no base gives no law
                self._restore(saved, ctx, solver=False)
                q0 = (S.T @ a) / A_p + q_pend
                q1 = (S.T @ b) / A_p
                Q = solver.step(dt, q0, q1)
                T_b1 = np.where(has, S @ solver.patch_T(), np.nan)
                book = self._pass(dt, T_b1, ctx, False, sheet)
                mismatch = book.E_base - np.where(has, (a + b * np.nan_to_num(T_b1)) * dt, 0.0)
                scale = np.abs(book.E_top) + np.abs(book.E_base)
                total = float(np.abs(mismatch).sum() / max(scale.sum(), 1e-300))
                r = np.where(has, np.nan_to_num(T_b1 - T_lin), 0.0)
                if total <= INTERFACE_TOLERANCE or np.abs(r).max(initial=0.0) <= DAMPED_TOL_K or \
                        attempt == DAMPED_MAX_PASSES - 1:
                    break
                if r_prev is not None:                          # Aitken: omega from the change of the residual
                    dr = r - r_prev
                    denom = float(dr @ dr)
                    if denom > 0.0:
                        omega = float(np.clip(-omega * float(r_prev @ dr) / denom, *DAMPED_OMEGA_RANGE))
                self._restore(saved, ctx)
                T_lin, r_prev = T_lin + omega * r, r
            passes = attempt + 1
        except Exception as err:                            # the part failed: back to where it started
            self._restore(saved, ctx)
            err.part_restored = True                        # `advance` may then ask whether rigid metal is left
            raise
        dC = sk.cells_volume() - sk.column_volume()          # the cells' volume change (planning ruling 3, step 1)
        move = self._moves()
        ratio = float(move.max() / (MOVE_SHARE * self.spacing)) if len(move) else 0.0
        return ratio, lambda: self._accept(dt, dt_macro, ratio, passes, book, mismatch, Q, E_vol0, share, base, has,
                                           dC, tally), \
            lambda: self._restore(saved, ctx)

    def _accept(self, dt, dt_macro, ratio, passes, book, mismatch, Q, E_vol0, share, base, has, dC, tally):
        """The accepted part's books (step 5) and the layer settled (step 6, `stepend.StepEnd._layer`): phi moved by the
        columns' volume changes dC [m3] (n,) and their carry, beta derived, the division and every column filled."""
        solver = self.solver
        tally.parts += 1
        tally.elapsed += dt
        tally.passes = max(tally.passes, passes)
        tally.ratio_all = max(tally.ratio_all, ratio * dt_macro / dt)       # the skyscrapers' rule, for the record
        tally.E_top += float(book.E_top.sum())
        tally.E_rad += float(book.E_rad.sum())
        tally.removed += float(book.removed.sum())
        tally.removed_H += float(book.removed_H.sum())
        tally.removed_p += book.removed_p.sum(axis=0)
        fb = book.flow
        tally.sprayed += float(fb.released.sum())
        tally.sprayed_H += float(fb.released_H.sum())
        tally.sprayed_p += fb.released_p.sum(axis=0)
        tally.wall -= fb.wall.sum(axis=0)                    # every skin's floor is no-slip (design 21.3)
        tally.mismatch = float(np.abs(mismatch).sum())
        # what the volume still owes: the mismatch per marker, and the solve's own residual shared by the base heat
        residual = Q - (solver.energy() - E_vol0)
        tally.newton += abs(residual)
        w = np.abs(book.E_base)
        weight = w / w.sum() if w.sum() > 0.0 else np.full(len(w), 1.0 / max(len(w), 1))
        owed = self.pending_J * (1.0 - share)                    # what this part did not deliver stays owed
        rows = np.flatnonzero(has) if has.any() else np.arange(len(base))
        if len(owed):                                            # (merged onto the base with a group nearest each point)
            _, k = cKDTree(base[rows]).query(self.pending_x)
            owed = np.bincount(rows[k], owed, minlength=len(base))
        else:
            owed = np.zeros(len(base))
        self.pending_x, self.pending_J = base.copy(), owed + mismatch + residual * weight
        if self._needs_move(dC + self.carry):
            self._layer(dC, tally, "part")

    # -- the macro step --------------------------------------------------------------------------------------------------
    def _take_loads(self, ctx):
        """The top fluxes [W/m2] (`ctx["q"]`) taken on phi as it is (module docstring; plan 03b Task 7b, the controller's
        ruling T7-i), each marker's under the ownership that sets its area (`Markers.owned_patches`, as
        `stepend.StepEnd._own` gives the areas), so that the heat the skins take, sum q_i A_i, is the heat on the surface
        as it is; and, with an aero state, at the step's start only, the runoff's forcing (`ctx["forcing"]`, held over the
        step: ruling T7b-a, below). With an aero state, `loads.by_patches` at phi's patches as they are (each at its own
        angle and leeward test) from the step's own loads state (`ctx["loads"]`: the state, the nose fit, the mass
        centroid and the top temperatures at the step's start; no new gas state, the state's flux redistributed over the
        moved surface); a prescribed heating (`ctx["prescribed"]`) replaces the flux: a function of the markers is called
        again, an array is held. Called at the step's start and, before a part, whenever phi has moved since
        (`_part_once`): the areas are given afresh at every acceptance that moves phi. (Taken once at the step's start, the
        fluxes stayed while the areas they multiply followed phi: on the stage 3 box the heat in ran 5.5 % over the beam's
        power over a step and up to 25 % in its last part once the face was uneven, and the excess dug a trough inside the
        rim; under one ownership and taken again, +0.5 %, no trough; diagnosis D1.)

        Why the forcing is held (the controller's ruling T7b-a): ruling T7-i's reason is the heat book, flux times area
        equal to the heat on the surface as it is, which the flux alone satisfies; no book or diagnosis asks for the
        forcing's re-take. Taken again at every part as well, measured once on the spraying fixture's first spraying step
        (plan 03b Task 7b), it folded the liquid into a sheet 0.7 mm under the cube's bottom face and the step stopped on
        THIN_METAL, while with the forcing held the step completed. Within a step the markers stand up to 2 mm off phi
        (diagnosis D2), and a re-taken forcing is theirs. `ctx["surface"]` is the phi the fluxes were taken on,
        `ctx["q_max"]` the largest flux so far."""
        mk, sk = self.markers, self.skins
        src = ctx["loads"]
        if src is not None:
            P_phi = cut_patches(self.mesh, self.phi)
            lee_p = loads_mod.leeward_patches(self.mesh, self.phi, P_phi, self.v_hat, self.leeward, src["centroid"])
            ang_p = nose.load_angles(np.arccos(np.clip(P_phi.normal @ self.v_hat, -1.0, 1.0)), ~lee_p)
            q, flow = loads_mod.by_patches(self.loads, src["state"], mk, P_phi, src["T_top"], ang_p, src["fit"],
                                           self.v_hat, W=mk.owned_patches(self.mesh, self.phi, P_phi))
            ctx["q"] = np.zeros(sk.n)
            ctx["q"][mk.skin_index] = q
            if ctx["surface"] is None:                    # the step's start: the forcing, held over the step (T7b-a)
                g = flow.deceleration * self.v_hat if src["decel_vec"] is None else np.asarray(src["decel_vec"], float)
                ctx["forcing"] = Forcing.from_surface_flow(flow, g, self.v_hat, self.material.liquid.sigma,
                                                           mk.kappa.sum(axis=1))
                ctx["taken"] = (P_phi, lee_p, flow, g)
        prescribed = ctx["prescribed"]
        if prescribed is not None:
            q = prescribed(mk) if callable(prescribed) else np.asarray(prescribed, float)
            ctx["q"] = np.zeros(sk.n)
            ctx["q"][mk.skin_index] = np.broadcast_to(q, (mk.n,))
        ctx["surface"] = self.phi
        ctx["q_max"] = max(ctx["q_max"], float(np.max(ctx["q"], initial=0.0)))

    def advance(self, dt, state=None, decel_vec=None, q_prescribed=None):
        """One macro step dt [s] (module docstring). `state` the aero state (the shared loads, the surface flow's forcing
        and, with a spray model, Girin's plan); `decel_vec` [m/s2] (3,) the acceleration the metal feels in the body's frame
        (default the surface flow's deceleration along v_hat; with no state it alone turns the flow sub-steps on, with no
        pressure gradient or shear); `q_prescribed` the convective heating [W/m2] per marker (an array in the markers'
        order, held over the step, or a function of the markers, called again before every part on a moved phi:
        `_take_loads`), instead of the state's. Returns a StepBook."""
        if self.division is None:                            # a body pickled before plan 03b Task 5
            self._adopt()
        mk, sk, mat = self.markers, self.skins, self.material
        start = self.books()
        M0 = start["mass_volume"] + start["mass_skins"]
        E0 = start["enthalpy_volume"] + start["enthalpy_skins"] + start["pending_J"]
        ctx = {"q": np.zeros(sk.n), "forcing": None, "plan": None, "loads": None, "prescribed": q_prescribed,
               "surface": None, "q_max": 0.0}
        if state is not None:
            if self.loads is None:
                raise ValueError("an aero state needs the body's GasLoads (loads=...)")
            centroid = self.mass_centroid()
            T_top = self._top_T()[mk.skin_index]
            ctx["loads"] = {"state": state, "centroid": centroid, "decel_vec": decel_vec,
                            "T_top": np.nan_to_num(T_top, nan=float(np.nanmean(T_top))),
                            "fit": loads_mod.fit_nose(mk, centroid, self.v_hat, "current", fallback=self._radius(centroid))}
        elif decel_vec is not None:
            still = _StillGas(np.zeros(mk.n), np.zeros(mk.n))
            ctx["forcing"] = Forcing.from_surface_flow(still, np.asarray(decel_vec, float), self.v_hat,
                                                       mat.liquid.sigma, mk.kappa.sum(axis=1))
        self._take_loads(ctx)
        if state is not None and self.spray_model is not None:      # Girin's plan once, at the step's start
            P_phi, lee_p, flow, g = ctx["taken"]
            lee = loads_mod.leeward(self.mesh, self.phi, P_phi, mk, self.v_hat, self.leeward, centroid, patch_lee=lee_p)
            ang = nose.load_angles(loads_mod.marker_angles(mk, self.v_hat), ~lee)
            windward = ang <= 0.5 * np.pi
            w_n = float(np.linalg.norm(g)) * np.maximum(np.cos(ang), 0.0)
            ctx["plan"] = spray_plan(sk, mk, flow, state, self.spray_model, windward, w_n, dt, self._radius(centroid),
                                     ctx["loads"]["fit"].used, slurry=self.slurry_spray, theta=ang, centre=centroid,
                                     v_hat=self.v_hat, division=self.division)
        self.plan = ctx["plan"]

        # the parts (Ruling 7)
        n = self._hint
        top = self._top_T()[mk.skin_index]                  # per marker
        if mat.T_rigid < np.inf:
            hot = np.any(np.nan_to_num(top, nan=-np.inf) >= mat.T_rigid)
            if not hot and self._top_start is not None:      # the last step's rise carried on, marker by marker: each
                x0, top0 = self._top_start                   # marker reads the nearest marker of the last step's start,
                _, near = cKDTree(x0).query(mk.x)            # so the look-ahead survives splits and merges
                with np.errstate(invalid="ignore"):
                    ahead = top + (top - top0[near]) * (dt / self._dt_last)
                hot = np.any(np.nan_to_num(ahead, nan=-np.inf) >= mat.T_rigid)
            if hot:
                n = max(n, int(math.ceil(dt / PART_MAX - 1e-9)))
        self._top_start, self._dt_last = (mk.x.copy(), top), dt
        tally = _Tally(sk.n)
        try:
            counts = [self._part(dt / n, dt, 0 if n == 1 else 1, ctx, tally)    # the skyscrapers' part loop starts a
                      for _ in range(n)]                                         # multi-part step one level down; matched
        except (RuntimeError, ValueError) as err:          # a part failed: say so plainly if no rigid metal is left,
            if getattr(err, "part_restored", False) and self._no_rigid_metal():    # judged on the restored state only
                raise RuntimeError(NO_CORE.format(mat.T_rigid)) from err
            raise
        self._hint = next_hint(counts)
        tally.depth = max(c.depth for c in counts)
        tally.over = sum(c.over for c in counts)
        self.time += dt

        # the step's end: the last carry, psi, the mesh, the markers, the layer
        self._end_step(tally, check_core=True)
        self.momentum_to_volume += tally.to_volume
        reps = self.last_settle
        passed = math.fsum(float(r.transfer.sum() - r.transfer.diagonal().sum()) for r in reps)
        covered = float(reps[-1].volume.sum()) if reps else 0.0
        summary = self._summary(state, ctx["plan"])

        end = self.books()
        M1 = end["mass_volume"] + end["mass_skins"]
        E1 = end["enthalpy_volume"] + end["enthalpy_skins"] + end["pending_J"]
        out_m = tally.sprayed + tally.removed
        out_E = tally.sprayed_H + tally.removed_H + tally.E_rad
        P = cut_patches(self.mesh, self.phi)
        label, n_pieces = pieces(self.mesh, self.phi)
        mk = self.markers
        self.last_step = StepBook(
            time=self.time, dt=dt, parts=tally.parts, parts_next=self._hint,
            parts_skyscrapers=int(math.ceil(MOVE_HINT_MARGIN * tally.ratio_all)), split_depth_max=tally.depth,
            parts_over=tally.over, snap_fallback=tally.snap_fallback, snap_slide_max=tally.snap_slide,
            pieces_used_up=tally.pieces_used, volume_repartitioned_share=passed / covered if covered > 0.0 else 0.0,
            stretch_min=min((float(r.stretch.min()) for r in reps), default=1.0),
            stretch_max=max((float(r.stretch.max()) for r in reps), default=1.0),
            group_area_ratio_min=tally.group_min, group_area_ratio_max=tally.group_max,
            group_area_ratio_trimmed_max=tally.group_trimmed, pending_flux_max=tally.pending_flux,
            surface_flux_max=ctx["q_max"], passes=tally.passes, enthalpy_in=tally.E_top,
            mass_sprayed=tally.sprayed, enthalpy_sprayed=tally.sprayed_H, mass_removed=tally.removed,
            enthalpy_removed=tally.removed_H, enthalpy_radiated=tally.E_rad, mass_volume=end["mass_volume"],
            mass_skins=end["mass_skins"], enthalpy_volume=end["enthalpy_volume"], enthalpy_skins=end["enthalpy_skins"],
            pending_J=end["pending_J"], mass_residual=(M1 - M0 + out_m) / M0,
            enthalpy_residual=(E1 - E0 - tally.E_top + out_E) / (abs(E0) + abs(tally.E_top)),
            mass_drawn=tally.drawn, mass_handed=tally.handed, momentum_sprayed=tally.sprayed_p,
            momentum_removed=tally.removed_p, momentum_to_volume=tally.to_volume, momentum_to_wall=tally.wall,
            residual_newton_J=tally.newton, residual_move_J=tally.move, residual_remesh_J=tally.remesh,
            n_markers=mk.n, marker_splits=mk.splits, marker_merges=mk.merges, markers_trapped=mk.trapped,
            interface_mismatch_J=tally.mismatch, volume_correction_max=tally.correction, n_pieces=int(n_pieces),
            frontal_area_m2=float((P.area * np.maximum(P.normal @ self.v_hat, 0.0)).sum()),
            layer_factor_max=tally.layer_factor, fill_miss_max=tally.fill_miss, move_miss_max=tally.move_miss,
            n_baseless=tally.baseless, density_carry_m3=float(self.carry.sum()),
            density_carry_share_max=self._carry_share(), island_nodes=tally.islands,
            groups_guarded_max=tally.guarded, move_corrections=list(tally.corrections), summary=summary)
        return self.last_step

    def _radius(self, centroid):
        """The body's transverse radius [m]: the largest distance of a marker from the flight axis through `centroid`."""
        X = self.markers.x - centroid
        return float(np.sqrt(np.maximum(np.einsum("ij,ij->i", X, X) - (X @ self.v_hat) ** 2, 0.0)).max())

    def _summary(self, state, plan):
        """The history's droplet columns (SKIN_SUMMARY_COLUMNS' keys and the release angle's) of the step; with a spray
        model on the profile path, last, the shared melt_profile.DISPERSION_COLUMNS (plan 03b Task 10, the skyscrapers'
        `dispersion_stats`): the plan's (`SprayPlan.step_summary`), or an empty step's where the step made no plan, so
        that every step of such a run has the same columns; none on the table path or without a spray model."""
        sk, mat = self.skins, self.material
        if plan is not None:
            block = plan.records(self.time, state)
            out = plan.step_summary(sk, self.n_released)
            self.n_released = out["n_released"]
            if len(block):
                self.source_rows.append(block)
            return out
        out = spray_mod.no_step_spray()
        flowing, liquid = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
        out.update({"time_s": self.time, "altitude_km": float(state.h) / 1e3 if state is not None else float("nan"),
                    "velocity_kms": float(state.V) / 1e3 if state is not None else float("nan"),
                    "released_step_mass_kg": 0.0, "n_released": self.n_released,
                    "liquid_depth_mean_mm": float(np.mean(liquid) * 1e3) if len(liquid) else 0.0,
                    "flowing_depth_mean_mm": float(np.mean(flowing) * 1e3) if len(flowing) else 0.0})
        if self._on_profile():
            out.update(melt_profile.summary(None, None, None, None))
        return out

    def _on_profile(self):
        """Whether the body's spray model takes the profile path (its `dispersion` "profile"); a body without one is on
        the table path."""
        return getattr(self.spray_model, "dispersion", "table") == "profile"

    def record_names(self):
        """The droplet records' column names (plan 03b Task 10; the skyscrapers' `SkinBody.record_names`):
        `spray.record_columns(36)` where the spray model is on the profile path, whose rows carry its PROFILE_COLUMNS, and
        the 27 SOURCE_COLUMNS otherwise -- what the shared writer (`droplets.write_particles`) writes, records or none."""
        return spray_mod.record_columns(36 if self._on_profile() else 27)
