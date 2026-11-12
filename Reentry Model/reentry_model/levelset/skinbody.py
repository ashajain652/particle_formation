"""The level set's skinned body (design sections 10.1, 12.2 and 20; level-set stage 3 Task 9): one macro step coupling the
skins on the gas-surface markers to the cut-cell 3D solve of stage 2. It wires together the skins (`skin`), the markers
(`markers`), the runoff (`runoff`), the gas loads (`loads`), the spraying (`spraying`), the thermal solve (`cutthermal`)
and the geometry (`boxmesh`, `levelset`, `surfaces`, `criteria`); the physics is theirs and the shared package's. It
mirrors the skyscrapers' macro step (their `SkinBody.advance` and its parts: the trial and real passes, the base law, the
damped repetition, the deferred mismatch) by amendment, never by import.

The body. Three surfaces (design 20.3): beta, the skin base, the evolved level set and the boundary of the volume's metal
(the CutThermal solve runs on {beta < 0}); phi, the outer surface, rebuilt each macro step as beta moved outward by the
skins' thickness (`surfaces.outer_surface`); psi, the molten-rigid surface, rebuilt from the temperature
(`surfaces.molten_rigid_surface`), the volume's on its nodes and each skin's cell profile between beta and phi. The
markers sit on phi (Ruling 1), one skin each (a row of `Skins`), carrying their skin's base at x - H n. Mass ownership:
the volume owns rho V(beta < 0), the skins their cells; nothing is counted twice.

The interface (the skyscrapers' per-column law converted per area). Every patch of beta (the solve's boundary pieces)
belongs to the marker whose base point is nearest its centroid, and a marker that no patch chose takes the patch nearest
its base: marker m's group G_m. A marker's base law a_m + b_m T_b [W, W/K, from its trial pass] reaches its group's
patches per area, q0 = sum a_m / A(G_m), q1 = sum b_m / A(G_m) over the markers whose group holds the patch, and its base
temperature T_b is the area mean over its group of the patch temperatures the solve applies its law at -- so the volume
receives exactly sum_m (a_m + b_m T_b,m) dt. (With one patch per marker this is the brief's q0 = sum a / A_patch.)

One macro step (`advance`, design 10.1 with 20): the gas loads once per marker at the step's start (`loads.GasLoads`,
the leeward test, the fitted nose), or a prescribed heating; the surface flow's forcing for the runoff (with an aero state
or a deceleration); Girin's plan once (`spraying.spray_plan`, with a spray model and a state); then the step as parts
(Ruling 7: none longer than PART_MAX while any top is at or above T_rigid or about to be, and a part whose base moves
exceed MOVE_SHARE of the markers' target spacing rerun as shorter parts, at most SPLIT_DEPTH_MAX levels). Each part:

1. the trial pass on saved skins (their cells, their flow memory, their curvatures, the plan's books): every heat
   sub-step's conduction
   with the trial sensitivity, the flow sub-steps with the release, the feed removal (remove="feed"), the draw and the
   regrid, at the linearisation temperature T_lin; its base law (a, b) per marker;
2. `CutThermal.step` with each patch's law and the pending load: the heat owed to the volume, delivered at its amount
   over the macro step's time left (at least PENDING_TIME_MIN), spread over the owing markers' patch groups;
3. the real pass from the same saved state with T_b read from the 3D field as above;
4. the mismatch per marker, E_base - (a + b T_b) dt, and the repetition from the saved states by the damped Aitken rule
   (omega0 DAMPED_OMEGA0, omega within DAMPED_OMEGA_RANGE, until the total relative mismatch is at most
   INTERFACE_TOLERANCE, no base moves more than DAMPED_TOL_K, or DAMPED_MAX_PASSES passes);
5. the accepted part's hand-down (`Skins.hand_down` below the stub, keeping at least the minimum mass, so that a hand-down
   by depth at the cells' own density and a draw by the solid's density never pass the same metal back and forth; a
   shortfall below DRAW_SLACK of the minimum is round-off and not drawn) and its books; what the volume still owes -- each
   marker's mismatch, and the solve's own residual (the heat its law applied less its stored enthalpy's change, the
   Newton iteration's error), shared by the markers in proportion to their base heat -- is the next part's pending load,
   so the enthalpy closes exactly;
6. beta moved by the part's draws and hand-downs: each marker's volume change (handed - drawn) / rho as a displacement
   along the normal over its group's area, extended to the nodes from the nearest base point, beta - delta, reset to a
   distance and corrected per piece to the exact volume target (`levelset.correct_volumes`); `CutThermal.move`, the
   enthalpy it removed against what the skins took and gave booked as pending at the markers' bases.

At the macro step's end (`stepend.StepEnd`, this class's mixin): phi and psi rebuilt, the mesh rebuilt where it no longer
meets its refinement targets, and the markers settled on phi; see that module.

Every skin here is full-depth (one mode): Ruling 3 and the two modes are Task 10's; each marker's mode is an attribute
(`mode`, "full" throughout) that Task 10 sets. No run driver (Ruling 6): tests and scripts drive `advance` at prescribed
gas states.

A skin deeper than its radius of curvature: each column's curvatures are capped (`stepend.cap_curvature`) before
every part, every heat sub-step and every flow sub-step (runoff can gather metal into one column) and at every sync;
the trial pass's caps are undone with its cells.

The part count (Ruling 7, fix round 1 of Task 9): the next step's parts are 1.25 times the largest move ratio of the
accepted parts that met their allowance, scaled to the macro step (the skyscrapers' estimate), but at most 1.25 times
the parts this step took and at most PARTS_MAX; a part accepted at the deepest split while still over its allowance
does not set it. (The skyscrapers' rule let such parts multiply the count step after step, their stall of 2026-11-08;
and a part that met its allowance only after deep splits, a burst of draws in a short window, gave an estimate of
1 000 for a step that had taken 33 parts, measured once on the test cube.) The look-ahead to T_rigid follows each marker
from the nearest marker of the last step's start, so it
survives splits and merges.

Every number here is a numerical setting (named with its source); the physics is the shared package's."""
import math
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

from ..shared import nose
from ..shared import spray as spray_mod
from ..shared.droplets import SKIN_SUMMARY_COLUMNS
from . import loads as loads_mod
from .boxmesh import build
from .cutcell import cell_gradients, inside_volumes
from .cutcell import patches as cut_patches
from .cutthermal import CutThermal
from .levelset import pieces, redistance
from .markers import Markers
from .runoff import FlowBook, Forcing, flow_substeps
from .skin import SkinSettings
from .spraying import spray_plan
from .stepend import KAPPA_DEPTH_SHARE, MODES, StepEnd, cap_curvature  # noqa: F401  (re-exported for the tests)
from .surfaces import outer_surface

# the parts of a macro step (Ruling 7, the skyscrapers' skinbody.py rules, matched by amendment)
PART_MAX = 0.0625                 # s: the longest part while any skin's top is at or above T_rigid, or about to be
MOVE_SHARE = 0.25                 # of the markers' target spacing: the most a base moves in one part (Ruling 7, its
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
PARTS_MAX = 1000                  # a macro step never starts with more parts: the next step's hint is capped here
# the interface's damped repetition (the skyscrapers' "damped" rule, their skinbody.py and skin.py)
INTERFACE_TOLERANCE = 1e-3        # of sum |E_top| + |E_base|: the total relative mismatch that ends the repetition
DAMPED_MAX_PASSES = 12            # passes of a part at most (the last pass's mismatch is deferred)
DAMPED_OMEGA0 = 0.5               # the first repetition's relaxation factor
DAMPED_OMEGA_RANGE = (0.02, 1.0)  # Aitken's factor kept within these
DAMPED_TOL_K = 0.05               # K: converged when no base temperature moves more than this between passes
# the markers and the geometry
MARKER_SPACING = 2.0e-3           # m: the markers' target spacing (design 20.13: about 30 000 at 2 mm on the 100 mm sphere)
INITIAL_REFINES = 2               # rebuilds of the first mesh from its own fields at the start
DRAW_SLACK = 1e-12                # of the minimum mass: a shortfall below this is round-off, and is not drawn (a draw of
                                  # 1e-21 kg would move beta)
FILL_DRAWS = 60                   # draws at most that bring a new skin to its minimum mass
FILL_TOL = 1e-15                  # of a skin's mass: the draw left when a new skin is full
INITIAL_REBUILDS = 4              # then step-end rebuilds at most until one leaves the mesh as it was (a cold step's mesh)
BAND_FACTOR = 2.0                 # distance resets are limited to this times (fine_margin + the base cell)
REMOVE_MODES = ("spray", "feed")  # liquid leaves by the spray release only, or as it reaches T_feed (the plates' flow off)

SKIN_SUMMARY_KEYS = tuple(k for _, k in SKIN_SUMMARY_COLUMNS)




@dataclass
class StepBook:
    """One macro step's books [kg, J, kg m/s] and the history's columns (design 12.2). In: `enthalpy_in` the convective
    heat taken at the tops (no mass enters). Out: `mass_sprayed`/`enthalpy_sprayed`/`momentum_sprayed` the release,
    `mass_removed`/`enthalpy_removed` the feed removal, `enthalpy_radiated`. Stored at the step's end: the volume's
    (rho V(beta < 0) and its enthalpy) and the skins'; `pending_J` the heat still owed to the volume (deferred to
    the next part). The residuals: (mass now - mass then + out) / mass then, and (stored + pending now - stored +
    pending then - in + out) / (|stored then| + |in|). `parts` the parts the step took, `parts_next` the next step's
    hint, and `parts_skyscrapers` the hint the skyscrapers' rule would have given (every accepted part's move ratio,
    unlimited). `group_area_ratio_min`, `group_area_ratio_max` the least and largest
    area of a marker's group of patches of beta over its skin's base area, over the step's parts (`patch_groups`). `pending_flux_max` the largest
    flux [W/m2] the pending heat put on a patch of beta in any part, `surface_flux_max` the largest convective flux on
    a top. `mass_drawn`, `mass_handed` the exchange at the bases; `momentum_to_volume` the momentum the held
    volume took over the step (the trapped markers' skins, Ruling 4; `LevelSetSkinBody.momentum_to_volume` keeps the
    total since the start), `momentum_to_wall` what the skins' flow gave the wall over the step and
    `momentum_removed` what the feed removal took. The heat deferred to the volume over the step, by source, each
    the sum of the absolute values of its bookings [J]: `residual_newton_J` the volume solve's own residual (the
    heat its law applied less its stored enthalpy's change), `residual_move_J` the base moves' (their enthalpy
    against what the skins took and gave), `residual_remesh_J` the remesh's; the interface mismatch is
    `interface_mismatch_J`. History columns: the markers and their counters (cumulative), `interface_mismatch_J`
    (the last part's sum of |mismatch|), `volume_correction_max`, `n_pieces`, `frontal_area_m2` (phi's area
    projected on the flight direction) and `summary` (SKIN_SUMMARY_COLUMNS' keys and the release angle's)."""
    time: float
    dt: float
    parts: int
    parts_next: int
    parts_skyscrapers: int
    group_area_ratio_min: float
    group_area_ratio_max: float
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
    summary: dict = field(default_factory=dict)

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
        self.parts, self.passes, self.ratio, self.mismatch, self.correction = 0, 0, 0.0, 0.0, 0.0
        self.ratio_all = 0.0
        self.group_min, self.group_max = np.inf, 0.0
        self.newton = self.move = self.remesh = 0.0
        self.elapsed, self.pending_flux = 0.0, 0.0


@dataclass
class _Pass:
    """What one real pass did, per skin row: the heat books [J], the draw, the feed removal and the flow's book."""
    E_top: np.ndarray
    E_rad: np.ndarray
    E_base: np.ndarray
    drawn: np.ndarray
    drawn_H: np.ndarray
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


def patch_groups(centroid, area, base, A_base):
    """Each marker's group of patches of beta and how its base law is spread over them (fix round 2 of Task 9): (S, A_group).
    S (n, P), sparse and row-stochastic, holds marker m's share A_p / A(G_m) of each patch p of its group G_m; A_group (n,)
    [m2] the groups' areas. Every patch first goes to the marker whose base point (`base`, (n, 3)) is nearest its
    `centroid` (P, 3); a group whose area is outside GROUP_AREA_RANGE times its marker's skin base area `A_base` (n,) [m2]
    is replaced by the patches nearest that marker's base, taken in order of distance until their area first reaches
    A_base (shared with other groups where they overlap). A patch a shrunk group gave up and no other group holds takes no
    law. The law stays exact whatever the groups: with q0 = S^T a / A_p, q1 = S^T b / A_p, sum q A over the patches is
    sum a and sum b, and a marker's base temperature S T_p is the area mean over its group. A zero A_base leaves the
    nearest-base groups (and a marker with none the patch nearest its base)."""
    n, npatch = len(base), len(area)
    A_base = np.broadcast_to(np.asarray(A_base, dtype=float), (n,))
    _, owner = cKDTree(base).query(centroid)
    A_near = np.bincount(owner, area, minlength=n)
    lo, hi = GROUP_AREA_RANGE
    ratio = np.divide(A_near, A_base, out=np.full(n, np.inf), where=A_base > 0.0)
    redo = np.flatnonzero(((ratio < lo) | (ratio > hi)) & (A_base > 0.0) | (A_near == 0.0))
    keep = ~np.isin(owner, redo)
    rows, cols = [owner[keep]], [np.flatnonzero(keep)]
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
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    A_group = np.bincount(rows, area[cols], minlength=n)
    S = sp.csr_matrix((area[cols] / A_group[rows], (rows, cols)), shape=(n, npatch))
    return S, A_group


_PLAN_BOOKS = ("mass", "mass_slurry", "number", "number_slurry", "enthalpy", "enthalpy_slurry", "solid",
               "solid_slurry", "momentum", "momentum_slurry", "hist_n", "hist_m")


class LevelSetSkinBody(StepEnd):
    """The skinned level-set body (module docstring). `material` the shared Material (or a test material with its
    interface), `base` the mesh's BoxBase, `beta0` a function of points (k, 3) [m] giving the initial skin base's level
    set (negative inside the metal), `T0` [K] the initial temperature, `settings` the SkinSettings. `loads` a
    `loads.GasLoads` (needed with an aero state), `spray_model` the shared SprayModel (Girin's release, with an aero
    state), `leeward` the leeward test ("shadow" or "centroid"), `emissivity` of the tops, `fine_margin`, `fine_size`,
    `cell_floor` [m] the refinement (design 20.4-20.5), `remove` "spray" or "feed". Also: `spacing` [m] the markers'
    target spacing, `v_hat` the flight direction in the body's frame, `T_ambient` [K] for the tops' radiation,
    `slurry_spray` (`--slurry-spray`), `linear_solver` the volume's ("amg" or "direct").

    The skins start full-depth with the stub (`settings.stub`) of metal at T0 standing on beta, so phi starts the stub
    outside beta0."""

    def __init__(self, material, base, beta0, T0, settings=None, loads=None, spray_model=None, leeward="shadow",
                 emissivity=0.0, fine_margin=4e-3, fine_size=0.7e-3, cell_floor=0.25e-3, remove="spray",
                 spacing=MARKER_SPACING, v_hat=(1.0, 0.0, 0.0), T_ambient=0.0, slurry_spray=False,
                 linear_solver="amg"):
        if remove not in REMOVE_MODES:
            raise ValueError("remove must be one of {}, got {!r}".format(REMOVE_MODES, remove))
        if leeward not in loads_mod.LEEWARD_NAMES:
            raise ValueError("leeward must be one of {}, got {!r}".format(loads_mod.LEEWARD_NAMES, leeward))
        if not callable(beta0):
            raise ValueError("beta0 must be a function of the points (k, 3) giving the initial skin base's level set")
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
        stub = self.settings.stub

        # the mesh, refined from its own fields: beta0 and phi = beta0 moved out by the stub (all rigid at the start)
        mesh = build(base)
        beta = redistance(mesh, np.asarray(beta0(mesh.points), float), band=self.band)
        for _ in range(INITIAL_REFINES):
            phi = redistance(mesh, beta - stub, band=self.band)
            psi = phi if T0 <= material.T_rigid else np.full(len(phi), np.inf)
            mesh = build(base, self._target(mesh, phi, psi))
            beta = redistance(mesh, np.asarray(beta0(mesh.points), float), band=self.band)
        self.mesh, self.beta = mesh, beta
        self.solver = CutThermal(material, emissivity=0.0, linear_solver=linear_solver)
        self.solver.setup(mesh, beta)
        self.solver.set_temperature(T0)

        # the markers on phi, seeded one per voxel of the target spacing, each with the stub as its skin
        phi = redistance(mesh, beta - stub, band=self.band)
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
        self._fill_stub(T0)
        self.phi = outer_surface(mesh, beta, self._thickness_nodal(), band=self.band)
        self.psi = self._psi(self.phi)
        self._settle(_Tally(0))
        for _ in range(INITIAL_REBUILDS):                # to the mesh the step's end rebuilds from these fields
            self._fill_stub(T0)
            if self._rebuild(_Tally(0)):
                break
        self._fill_stub(T0)

    def _fill_stub(self, T0):
        """Every skin anew: the stub's depth of metal at T0 under its (capped) cross-section, then drawn at T0 up to its
        minimum mass (a concave column's stub holds less than its minimum, rho_s A(z_bottom) stub) until no draw is
        left above round-off, and regridded."""
        mk, sk, stub = self.markers, self.skins, self.settings.stub
        sk.area[mk.skin_index] = mk.area
        sk.kappa[mk.skin_index] = cap_curvature(sk, mk.kappa, np.full(mk.n, stub))
        s, p = sk.kappa.sum(axis=1), sk.kappa.prod(axis=1)
        sk.fill(self.material.rho * sk.area * stub * (1.0 - stub * (0.5 * s - stub * p / 3.0)), T0)
        for _ in range(FILL_DRAWS):
            dm, _ = sk.draw(T0, 0.0, sk.minimum_mass())
            if not np.any(dm > FILL_TOL * sk.m.sum(axis=1)):
                break
        sk.regrid()

    # -- small helpers ------------------------------------------------------------------------------------------------






    def volume_energy(self):
        return self.solver.energy()

    def books(self):
        """The body's stored mass [kg] and enthalpy [J] now: the volume's, the skins' and the pending heat."""
        sk = self.skins
        return {"mass_volume": self.rho * float(inside_volumes(self.mesh, self.beta).sum()),
                "mass_skins": float(sk.m.sum()), "enthalpy_volume": self.solver.energy(),
                "enthalpy_skins": float(sk.H.sum()), "pending_J": float(self.pending_J.sum())}

    def add_pending(self, x, J):
        """Book J [J] of heat owed to the volume at the point x [m] (delivered at the next part, at the nearest patch)."""
        self.pending_x = np.vstack([self.pending_x, np.atleast_2d(np.asarray(x, float))])
        self.pending_J = np.r_[self.pending_J, np.atleast_1d(np.asarray(J, float))]

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

    # -- psi: the volume's temperatures and the skins' profiles --------------------------------------------------------


    # -- the interface ---------------------------------------------------------------------------------------------------
    def _interface(self):
        """(S (n, P) sparse, row-stochastic: marker row m's share A_p / A(G_m) of each patch of its group; A(G) (n,) [m2];
        the base points (n, 3); the base gradient dT/dn [K/m] per row) on the solve's present patches."""
        P = self.solver.P
        base = self._bases()
        S, A_group = patch_groups(P.centroid, P.area, base, self.skins.base_area())
        g_patch = np.einsum("pk,pk->p", cell_gradients(self.mesh, self.solver.T, P.cell), P.normal)
        return S, A_group, base, S @ g_patch

    def _pending_per_patch(self, S, base):
        """The pending heat [J] per patch of beta: each pending point's heat goes to the marker whose base point is
        nearest it and is spread over that marker's group as its base law is (`patch_groups`, S), so that the sum is
        exact. (Fix round 2: deposited whole on the single patch nearest each point, a few joules made -6e12 W/m2 on a
        sliver over a 0.5 ms part and the volume went to -900 K inside the metal.)"""
        out = np.zeros(S.shape[1])
        if len(self.pending_J):
            _, row = cKDTree(base).query(self.pending_x)
            out = S.T @ np.bincount(row, self.pending_J, minlength=S.shape[0])
        return out

    # -- one pass of the skins -----------------------------------------------------------------------------------------
    def _pass(self, dt, T_b, g, m_min, ctx, trial):
        """The skins over one part dt [s] with the base held at T_b [K] (n,): per heat sub-step, the conduction, the flow
        sub-steps with the release, the feed removal, the draw (each sub-step's layer below the last) and the regrid.
        A trial pass returns its base law (a, b); a real pass its books (_Pass)."""
        sk, s = self.skins, self.settings
        n = sk.n
        N = max(int(math.ceil(dt / s.substep - 1e-9)), 1)
        h = dt / N
        z = np.zeros(n)
        E_top, E_rad, E_base = z.copy(), z.copy(), z.copy()
        drawn, drawn_H, removed, removed_H, removed_p = z.copy(), z.copy(), z.copy(), z.copy(), np.zeros((n, 3))
        fbook = FlowBook.zero(n)
        s0 = np.zeros(n)
        release = ctx["plan"].release if ctx["plan"] is not None else None
        for _ in range(N):
            self._cap()                                     # runoff may have deepened a column since the last sub-step
            hb = sk.conduct(h, ctx["q"], self.T_ambient, self.emissivity, T_b, trial=trial)
            E_top += hb.E_top
            E_rad += hb.E_rad
            E_base += hb.E_base
            if ctx["forcing"] is not None:
                fbook.add(flow_substeps(sk, self.markers, ctx["forcing"], h, release=release,
                                        after=lambda skins: self._cap()))
            if self.remove == "feed":
                mass, H = sk.remove_feed()
                removed += mass
                removed_H += H
                removed_p += sk.momentum_out
            A_b = sk.base_area()
            dm, dH = sk.draw(T_b, g, m_min, s0)
            s0 = s0 + dm / (sk.rho_s * A_b)
            drawn += dm
            drawn_H += dH
            sk.regrid()
        if trial:
            return sk.base_law(dt)
        return _Pass(E_top, E_rad, E_base, drawn, drawn_H, removed, removed_H, removed_p, fbook)

    def _save(self, ctx):
        sk, plan = self.skins, ctx["plan"]
        books = None if plan is None else {k: np.array(getattr(plan, k), copy=True) for k in _PLAN_BOOKS}
        return self.solver.T.copy(), sk.state(), sk.flow_memory.copy(), sk.kappa.copy(), books

    def _restore(self, saved, ctx, solver=True):
        """Undo a pass: the cells, the flow memory, the curvatures (which the per-sub-step cap may have lowered) and
        the plan's books, and with `solver` the volume's temperatures."""
        T, state, memory, kappa, books = saved
        if solver:
            self.solver.T = T.copy()
        self.skins.set_state(state)
        self.skins.flow_memory = memory.copy()
        self.skins.kappa = kappa.copy()
        if books is not None:
            for k, v in books.items():
                setattr(ctx["plan"], k, np.array(v, copy=True))

    # -- one part --------------------------------------------------------------------------------------------------------
    def _part(self, dt, dt_macro, depth, ctx, tally):
        """One part of a macro step (module docstring, steps 1-6); a part whose base moves exceed MOVE_SHARE of the local
        marker spacing is run again as shorter parts (Ruling 7)."""
        sk, solver = self.skins, self.solver
        self._cap()
        S, A_group, base, g0 = self._interface()
        group_ratio = A_group / sk.base_area()
        tally.group_min = min(tally.group_min, float(group_ratio.min()))
        tally.group_max = max(tally.group_max, float(group_ratio.max()))
        saved = self._save(ctx)
        T_b0 = S @ solver.patch_T()
        m_min = sk.minimum_mass() * (1.0 - DRAW_SLACK)
        A_p = solver.P.area
        share = min(1.0, dt / max(dt_macro - tally.elapsed, PENDING_TIME_MIN))   # of the pending heat, this part's
        q_pend = share * self._pending_per_patch(S, base) / (A_p * dt)              # (at its rate over the step's rest)
        tally.pending_flux = max(tally.pending_flux, float(np.abs(q_pend).max(initial=0.0)))
        E_vol0 = solver.energy()
        T_lin = T_b0
        omega, r_prev = DAMPED_OMEGA0, None
        for attempt in range(DAMPED_MAX_PASSES):
            a, b = self._pass(dt, T_lin, g0, m_min, ctx, trial=True)
            self._restore(saved, ctx, solver=False)
            q0 = (S.T @ a) / A_p + q_pend
            q1 = (S.T @ b) / A_p
            Q = solver.step(dt, q0, q1)
            T_b1 = S @ solver.patch_T()
            book = self._pass(dt, T_b1, g0, m_min, ctx, trial=False)
            mismatch = book.E_base - (a + b * T_b1) * dt
            scale = np.abs(book.E_top) + np.abs(book.E_base)
            total = float(np.abs(mismatch).sum() / max(scale.sum(), 1e-300))
            r = T_b1 - T_lin
            if total <= INTERFACE_TOLERANCE or np.abs(r).max() <= DAMPED_TOL_K or attempt == DAMPED_MAX_PASSES - 1:
                break
            if r_prev is not None:                          # Aitken: omega from the change of the residual
                dr = r - r_prev
                denom = float(dr @ dr)
                if denom > 0.0:
                    omega = float(np.clip(-omega * float(r_prev @ dr) / denom, *DAMPED_OMEGA_RANGE))
            self._restore(saved, ctx)
            T_lin, r_prev = T_lin + omega * r, r
        passes = attempt + 1
        handed, handed_H, _ = sk.hand_down(self.settings.stub, m_keep=sk.minimum_mass())
        A_base = sk.base_area()
        move = np.abs(book.drawn - handed) / (self.rho * A_base)
        allowance = MOVE_SHARE * self.spacing
        ratio = float((move / allowance).max()) if len(move) else 0.0
        if ratio > 1.0 and depth < SPLIT_DEPTH_MAX:         # too far for one part: run it again as shorter parts
            self._restore(saved, ctx)
            k = int(math.ceil(ratio))
            for _ in range(k):
                self._part(dt / k, dt_macro, depth + 1, ctx, tally)
            return
        # the accepted part's books
        tally.parts += 1
        tally.elapsed += dt
        tally.passes = max(tally.passes, passes)
        tally.ratio_all = max(tally.ratio_all, ratio * dt_macro / dt)       # the skyscrapers' rule, for the record
        if ratio <= 1.0:                                    # a part still over its allowance at the deepest split is
            tally.ratio = max(tally.ratio, ratio * dt_macro / dt)   # accepted but does not set the next step's parts
        tally.E_top += float(book.E_top.sum())
        tally.E_rad += float(book.E_rad.sum())
        tally.removed += float(book.removed.sum())
        tally.removed_H += float(book.removed_H.sum())
        tally.removed_p += book.removed_p.sum(axis=0)
        tally.drawn += float(book.drawn.sum())
        tally.handed += float(handed.sum())
        fb = book.flow
        tally.sprayed += float(fb.released.sum())
        tally.sprayed_H += float(fb.released_H.sum())
        tally.sprayed_p += fb.released_p.sum(axis=0)
        tally.wall -= fb.wall.sum(axis=0)
        tally.mismatch = float(np.abs(mismatch).sum())
        # what the volume still owes: the mismatch per marker, and the solve's own residual shared by the base heat
        residual = Q - (solver.energy() - E_vol0)
        tally.newton += abs(residual)
        w = np.abs(book.E_base)
        weight = w / w.sum() if w.sum() > 0.0 else np.full(len(w), 1.0 / max(len(w), 1))
        owed = self.pending_J * (1.0 - share)                    # what this part did not deliver stays owed
        if len(owed):                                            # (merged onto the base nearest each point)
            _, row = cKDTree(base).query(self.pending_x)
            owed = np.bincount(row, owed, minlength=len(base))
        else:
            owed = np.zeros(len(base))
        self.pending_x, self.pending_J = base.copy(), owed + mismatch + residual * weight
        # beta moved by the draws and the hand-downs, the move's enthalpy against what the skins took and gave
        dV = (handed - book.drawn) / self.rho
        if np.any(dV != 0.0):
            expected = float(handed_H.sum() - book.drawn_H.sum())
            removed = self._move_beta(base, dV, dV / A_base, tally)
            tally.move += abs(expected + removed)
            self._book_move(base, np.abs(dV), expected + removed)




    # -- the macro step --------------------------------------------------------------------------------------------------
    def advance(self, dt, state=None, decel_vec=None, q_prescribed=None):
        """One macro step dt [s] (module docstring). `state` the aero state (the shared loads, the surface flow's forcing
        and, with a spray model, Girin's plan); `decel_vec` [m/s2] (3,) the acceleration the metal feels in the body's frame
        (default the surface flow's deceleration along v_hat; with no state it alone turns the flow sub-steps on, with no
        pressure gradient or shear); `q_prescribed` the convective heating [W/m2] per marker (an array in the markers'
        order, or a function of the markers), instead of the state's. Returns a StepBook."""
        mk, sk, mat = self.markers, self.skins, self.material
        start = self.books()
        M0 = start["mass_volume"] + start["mass_skins"]
        E0 = start["enthalpy_volume"] + start["enthalpy_skins"] + start["pending_J"]
        splits0, merges0 = mk.splits, mk.merges
        ctx = {"q": np.zeros(sk.n), "forcing": None, "plan": None}
        T_top = self._top_T()[mk.skin_index]
        flow = None
        if state is not None:
            if self.loads is None:
                raise ValueError("an aero state needs the body's GasLoads (loads=...)")
            centroid = self.mass_centroid()
            theta = loads_mod.marker_angles(mk, self.v_hat)
            lee = loads_mod.leeward(self.mesh, self.phi, cut_patches(self.mesh, self.phi), mk, self.v_hat, self.leeward,
                                    centroid)
            ang = nose.load_angles(theta, ~lee)
            fit = loads_mod.fit_nose(mk, centroid, self.v_hat, "current", fallback=self._radius(centroid))
            q, flow = self.loads.evaluate(state, mk, np.nan_to_num(T_top, nan=float(np.nanmean(T_top))), ang, fit)
            ctx["q"][mk.skin_index] = q
            g = flow.deceleration * self.v_hat if decel_vec is None else np.asarray(decel_vec, float)
            ctx["forcing"] = Forcing.from_surface_flow(flow, g, self.v_hat, mat.liquid.sigma, mk.kappa.sum(axis=1))
            if self.spray_model is not None:
                windward = ang <= 0.5 * np.pi
                w_n = float(np.linalg.norm(g)) * np.maximum(np.cos(ang), 0.0)
                ctx["plan"] = spray_plan(sk, mk, flow, state, self.spray_model, windward, w_n, dt,
                                         self._radius(centroid), fit.used, slurry=self.slurry_spray, theta=ang,
                                         centre=centroid, v_hat=self.v_hat)
        elif decel_vec is not None:
            still = _StillGas(np.zeros(mk.n), np.zeros(mk.n))
            ctx["forcing"] = Forcing.from_surface_flow(still, np.asarray(decel_vec, float), self.v_hat,
                                                       mat.liquid.sigma, mk.kappa.sum(axis=1))
        if q_prescribed is not None:
            q = q_prescribed(mk) if callable(q_prescribed) else np.asarray(q_prescribed, float)
            ctx["q"] = np.zeros(sk.n)
            ctx["q"][mk.skin_index] = np.broadcast_to(q, (mk.n,))
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
        for _ in range(n):                                  # the skyscrapers' part loop starts a multi-part step one
            self._part(dt / n, dt, 0 if n == 1 else 1, ctx, tally)   # split level down; matched
        estimate = int(math.ceil(MOVE_HINT_MARGIN * tally.ratio))         # the skyscrapers' estimate, limited by the
        used = int(math.ceil(MOVE_HINT_MARGIN * tally.parts))               # parts this step actually took
        self._hint = min(max(1, min(estimate, used)), PARTS_MAX)
        self.time += dt

        # the step's end: phi, psi, the mesh, the markers
        self._rebuild(tally)
        self.momentum_to_volume += tally.to_volume
        summary = self._summary(state, ctx["plan"])

        end = self.books()
        M1 = end["mass_volume"] + end["mass_skins"]
        E1 = end["enthalpy_volume"] + end["enthalpy_skins"] + end["pending_J"]
        out_m = tally.sprayed + tally.removed
        out_E = tally.sprayed_H + tally.removed_H + tally.E_rad
        P = cut_patches(self.mesh, self.phi)
        label, n_pieces = pieces(self.mesh, self.beta)
        mk = self.markers
        self.last_step = StepBook(
            time=self.time, dt=dt, parts=tally.parts, parts_next=self._hint,
            parts_skyscrapers=int(math.ceil(MOVE_HINT_MARGIN * tally.ratio_all)),
            group_area_ratio_min=tally.group_min, group_area_ratio_max=tally.group_max,
            pending_flux_max=tally.pending_flux, surface_flux_max=float(np.max(ctx["q"], initial=0.0)), passes=tally.passes, enthalpy_in=tally.E_top,
            mass_sprayed=tally.sprayed, enthalpy_sprayed=tally.sprayed_H, mass_removed=tally.removed,
            enthalpy_removed=tally.removed_H, enthalpy_radiated=tally.E_rad, mass_volume=end["mass_volume"],
            mass_skins=end["mass_skins"], enthalpy_volume=end["enthalpy_volume"], enthalpy_skins=end["enthalpy_skins"],
            pending_J=end["pending_J"], mass_residual=(M1 - M0 + out_m) / M0,
            enthalpy_residual=(E1 - E0 - tally.E_top + out_E) / (abs(E0) + abs(tally.E_top)),
            mass_drawn=tally.drawn, mass_handed=tally.handed, momentum_sprayed=tally.sprayed_p,
            momentum_removed=tally.removed_p, momentum_to_volume=tally.to_volume, momentum_to_wall=tally.wall,
            residual_newton_J=tally.newton, residual_move_J=tally.move, residual_remesh_J=tally.remesh,
            n_markers=mk.n,
            marker_splits=mk.splits, marker_merges=mk.merges, markers_trapped=mk.trapped,
            interface_mismatch_J=tally.mismatch, volume_correction_max=tally.correction, n_pieces=int(n_pieces),
            frontal_area_m2=float((P.area * np.maximum(P.normal @ self.v_hat, 0.0)).sum()), summary=summary)
        return self.last_step

    def _radius(self, centroid):
        """The body's transverse radius [m]: the largest distance of a marker from the flight axis through `centroid`."""
        X = self.markers.x - centroid
        return float(np.sqrt(np.maximum(np.einsum("ij,ij->i", X, X) - (X @ self.v_hat) ** 2, 0.0)).max())

    def _summary(self, state, plan):
        """The history's droplet columns (SKIN_SUMMARY_COLUMNS' keys and the release angle's) of the step."""
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
        return out

    # -- the step's end ----------------------------------------------------------------------------------------------------





