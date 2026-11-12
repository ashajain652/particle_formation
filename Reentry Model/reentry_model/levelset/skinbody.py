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
exceed MOVE_SHARE of the local marker spacing rerun as shorter parts, at most SPLIT_DEPTH_MAX levels). Each part:

1. the trial pass on saved skins (their cells, their flow memory, the plan's books): every heat sub-step's conduction
   with the trial sensitivity, the flow sub-steps with the release, the feed removal (remove="feed"), the draw and the
   regrid, at the linearisation temperature T_lin; its base law (a, b) per marker;
2. `CutThermal.step` with each patch's law and the pending load (last part's deferred heat, at its points);
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

At the macro step's end: phi and psi rebuilt; the mesh rebuilt from the refinement targets (`criteria.make_target`: the
fine region psi > -fine_margin at fine_size, curvature down to cell_floor) where it no longer meets them
(`_needs_rebuild`), with `CutThermal.remesh` (beta carried to the new
mesh, reset and corrected to the exact volume; the remesh's enthalpy residual booked as pending); phi and psi again on the
new mesh; the markers snapped onto phi (to SNAP_TOL), the trapped ones released (their skins' metal and enthalpy to the
volume, beta moved by their volume, their momentum booked to the volume), tessellated, split and merged, snapped and
tessellated again, the skins regridded and given the markers' areas and curvatures.

Every skin here is full-depth (one mode): Ruling 3 and the two modes are Task 10's; each marker's mode is an attribute
(`mode`, "full" throughout) that Task 10 sets. No run driver (Ruling 6): tests and scripts drive `advance` at prescribed
gas states.

A skin deeper than its radius of curvature (`Skins.dz` raises for all columns if one is): before every part, every heat
sub-step and every flow sub-step (runoff can gather metal into one column) and at every sync, each column's principal
curvatures are capped at KAPPA_DEPTH_SHARE over its flat depth (`cap_curvature`), so the
cross-section A0 (1 - k1 z)(1 - k2 z) stays positive to three times the column's depth; only the columns over the cap
change. A surface bending more tightly than its skin's depth is sub-grid (design 20.6, known limitation 1); the markers
keep phi's curvature for the capillary push.

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
from .criteria import finite_values, kappa_max_nodal, make_target
from .cutcell import cell_gradients, inside_volumes
from .cutcell import patches as cut_patches
from .cutthermal import CutThermal
from .levelset import correct_volumes, pieces, redistance
from .markers import Markers
from .runoff import FlowBook, Forcing, flow_substeps
from .skin import SkinSettings
from .spraying import spray_plan
from .surfaces import molten_rigid_surface, outer_surface
from .transfer import interpolate, locate

# the parts of a macro step (Ruling 7, the skyscrapers' skinbody.py rules, matched by amendment)
PART_MAX = 0.0625                 # s: the longest part while any skin's top is at or above T_rigid, or about to be
MOVE_SHARE = 0.25                 # of the local marker spacing: the most a base moves in one part (Ruling 7)
SPLIT_DEPTH_MAX = 6               # a part is rerun shorter at most this many times over (Ruling 7)
MOVE_HINT_MARGIN = 1.25           # the next macro step's parts: this times the parts this one's largest move called for
PARTS_MAX = 1000                  # a macro step that would start with more parts stops the run (the skyscrapers')
# the interface's damped repetition (the skyscrapers' "damped" rule, their skinbody.py and skin.py)
INTERFACE_TOLERANCE = 1e-3        # of sum |E_top| + |E_base|: the total relative mismatch that ends the repetition
DAMPED_MAX_PASSES = 12            # passes of a part at most (the last pass's mismatch is deferred)
DAMPED_OMEGA0 = 0.5               # the first repetition's relaxation factor
DAMPED_OMEGA_RANGE = (0.02, 1.0)  # Aitken's factor kept within these
DAMPED_TOL_K = 0.05               # K: converged when no base temperature moves more than this between passes
# the markers and the geometry
MARKER_SPACING = 2.0e-3           # m: the markers' target spacing (design 20.13: about 30 000 at 2 mm on the 100 mm sphere)
KAPPA_DEPTH_SHARE = 0.1           # a column's principal curvature at most this over its flat depth (module docstring)
SNAP_TOL = 1.0e-10                # m: the markers are snapped onto phi until none is farther from its zero contour
SNAP_ROUNDS = 20                  # rounds of `Markers.snap` at most
INITIAL_REFINES = 2               # rebuilds of the first mesh from its own fields at the start
DRAW_SLACK = 1e-12                # of the minimum mass: a shortfall below this is round-off, and is not drawn (a draw of
                                  # 1e-21 kg would move beta)
FILL_DRAWS = 60                   # draws at most that bring a new skin to its minimum mass
FILL_TOL = 1e-15                  # of a skin's mass: the draw left when a new skin is full
LEAD = 0.6e-3                     # m: the fine region's lead beyond the margin (criteria.make_target's default)
REGION_SLACK = 0.5                # fine cells: a rebuild is needed once a leaf is too coarse for a lead this much shorter
CURVATURE_SLACK = 2.0             # or once a leaf is more than this times its curvature's target size
COARSEN_SHARE = 0.1               # of the metal's volume: leaves that could be COARSEN_SIZES times coarser before the mesh
COARSEN_SIZES = 4.0               # is rebuilt (two sizes: the closure's grading keeps one-size steps finer than asked)
INITIAL_REBUILDS = 4              # then step-end rebuilds at most until one leaves the mesh as it was (a cold step's mesh)
BAND_FACTOR = 2.0                 # distance resets are limited to this times (fine_margin + the base cell)
MODES = ("full", "top")           # a marker's skin: full-depth (Task 9, everywhere) or top-zone (Task 10)
REMOVE_MODES = ("spray", "feed")  # liquid leaves by the spray release only, or as it reaches T_feed (the plates' flow off)

SKIN_SUMMARY_KEYS = tuple(k for _, k in SKIN_SUMMARY_COLUMNS)


def cap_curvature(skins, kappa, depth=None):
    """Each column's principal curvatures (n, 2) [1/m] with the positive ones capped at KAPPA_DEPTH_SHARE over the
    column's flat depth (module docstring), `depth` [m] (n,) or the skins' flat depth sum m / (rho(T) A0). Columns under
    the cap, and negative (concave) curvatures, are unchanged."""
    kappa = np.asarray(kappa, dtype=float)
    if depth is None:
        T, valid = skins.T(), skins.valid
        vol = np.zeros(skins.m.shape)
        vol[valid] = skins.m[valid] / skins.density(T[valid])
        depth = vol.sum(axis=1) / skins.area
    limit = KAPPA_DEPTH_SHARE / np.maximum(np.asarray(depth, dtype=float), 1e-300)
    return np.minimum(kappa, limit[:, None])


@dataclass
class StepBook:
    """One macro step's books [kg, J, kg m/s] and the history's columns (design 12.2). In: `enthalpy_in` the convective
    heat taken at the tops (no mass enters). Out: `mass_sprayed`/`enthalpy_sprayed`/`momentum_sprayed` the release,
    `mass_removed`/`enthalpy_removed` the feed removal, `enthalpy_radiated`. Stored at the step's end: the volume's
    (rho V(beta < 0) and its enthalpy) and the skins'; `pending_J` the heat still owed to the volume (deferred to the next
    part). The residuals: (mass now - mass then + out) / mass then, and (stored + pending now - stored + pending then - in
    + out) / (|stored then| + |in|). `mass_drawn`, `mass_handed` the exchange at the bases; `momentum_to_volume` the
    momentum the held volume took (the trapped markers' skins, Ruling 4) and `momentum_to_wall` the impulse the skins'
    flow gave the wall. History columns: the markers and their counters (cumulative), `interface_mismatch_J` (the last
    part's sum of |mismatch|), `volume_correction_max`, `n_pieces`, `frontal_area_m2` (phi's area projected on the flight
    direction) and `summary` (SKIN_SUMMARY_COLUMNS' keys and the release angle's)."""
    time: float
    dt: float
    parts: int
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
    momentum_to_volume: np.ndarray
    momentum_to_wall: np.ndarray
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
        self.sprayed_p, self.wall = np.zeros(3), np.zeros(3)
        self.parts, self.passes, self.ratio, self.mismatch, self.correction = 0, 0, 0.0, 0.0, 0.0


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


_PLAN_BOOKS = ("mass", "mass_slurry", "number", "number_slurry", "enthalpy", "enthalpy_slurry", "solid",
               "solid_slurry", "momentum", "momentum_slurry", "hist_n", "hist_m")


class LevelSetSkinBody:
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
        self._hint, self._top_start, self._dt_last = 1, None, None
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
        self._settle()
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
    def _target(self, mesh, phi, psi, curvature=True, lead=LEAD):
        band = 2.0 * self.base.h0
        if curvature:
            kp, ks = kappa_max_nodal(mesh, phi, band=band), kappa_max_nodal(mesh, finite_values(psi), band=band)
        else:
            kp = ks = np.zeros(len(phi))
        return make_target(mesh, phi, psi, kp, ks, fine=self.fine_size, margin=self.fine_margin, lead=lead,
                           floor=self.cell_floor)

    def _row_marker(self):
        """Each skin row's marker (n,)."""
        out = np.empty(self.skins.n, dtype=np.int64)
        out[self.markers.skin_index] = np.arange(self.markers.n)
        return out

    def _bases(self):
        """Each skin row's base point [m] (n, 3), x - H n (Ruling 1)."""
        mk, sk = self.markers, self.skins
        r = self._row_marker()
        return mk.x[r] - sk.thickness()[:, None] * mk.normal[r]

    def _thickness_nodal(self):
        """The skins' thickness [m] at every node, from the marker whose base point is nearest (design 20.3)."""
        _, near = cKDTree(self._bases()).query(self.mesh.points)
        return self.skins.thickness()[near]

    def _sync(self):
        """The skins take their markers' areas and (capped) curvatures."""
        mk, sk = self.markers, self.skins
        sk.area[mk.skin_index] = mk.area
        sk.kappa[mk.skin_index] = mk.kappa
        sk.kappa = cap_curvature(sk, sk.kappa)

    def _cap(self):
        self.skins.kappa = cap_curvature(self.skins, self.skins.kappa)

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
    def _temperature_field(self, phi):
        """The temperature at every node: the volume's inside beta; between beta and phi each skin's cell at the node's
        depth below phi (the nearest marker's), and outside phi its top cell (so that psi sees no crossing there)."""
        T = self.solver.T.copy()
        out = self.beta >= 0.0
        if not out.any():
            return T
        sk = self.skins
        Tc = sk.T()
        dz = sk.dz(Tc)
        z_bot = np.where(sk.valid, np.cumsum(dz, axis=1), np.inf)
        nodes = np.flatnonzero(out)
        _, near = cKDTree(self.markers.x).query(self.mesh.points[nodes])
        r = self.markers.skin_index[near]
        depth = np.maximum(-phi[nodes], 0.0)
        cell = np.minimum((z_bot[r] < depth[:, None]).sum(axis=1), np.maximum(sk.count[r] - 1, 0))
        T[nodes] = np.where(sk.count[r] > 0, Tc[r, cell], T[nodes])
        return T

    def _psi(self, phi):
        psi, _ = molten_rigid_surface(self.mesh, phi, self._temperature_field(phi), self.material.T_rigid,
                                      band=self.band)
        return psi

    # -- the interface ---------------------------------------------------------------------------------------------------
    def _interface(self):
        """(S (n, P) sparse, row-stochastic: marker row m's share A_p / A(G_m) of each patch of its group; A(G) (n,) [m2];
        the base points (n, 3); the base gradient dT/dn [K/m] per row) on the solve's present patches."""
        P = self.solver.P
        base = self._bases()
        n, npatch = self.skins.n, len(P.area)
        _, owner = cKDTree(base).query(P.centroid)
        lone = np.flatnonzero(np.bincount(owner, minlength=n) == 0)
        rows, cols = owner, np.arange(npatch)
        if len(lone):
            _, pp = cKDTree(P.centroid).query(base[lone])
            rows, cols = np.r_[rows, lone], np.r_[cols, pp]
        A_group = np.bincount(rows, P.area[cols], minlength=n)
        S = sp.csr_matrix((P.area[cols] / A_group[rows], (rows, cols)), shape=(n, npatch))
        g_patch = np.einsum("pk,pk->p", cell_gradients(self.mesh, self.solver.T, P.cell), P.normal)
        return S, A_group, base, S @ g_patch

    def _pending_per_patch(self):
        P = self.solver.P
        out = np.zeros(len(P.area))
        if len(self.pending_J):
            _, near = cKDTree(P.centroid).query(self.pending_x)
            out = np.bincount(near, self.pending_J, minlength=len(P.area))
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
        return self.solver.T.copy(), sk.state(), sk.flow_memory.copy(), books

    def _restore(self, saved, ctx, solver=True):
        T, state, memory, books = saved
        if solver:
            self.solver.T = T.copy()
        self.skins.set_state(state)
        self.skins.flow_memory = memory.copy()
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
        saved = self._save(ctx)
        T_b0 = S @ solver.patch_T()
        m_min = sk.minimum_mass() * (1.0 - DRAW_SLACK)
        A_p = solver.P.area
        q_pend = self._pending_per_patch() / (A_p * dt)
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
        allowance = MOVE_SHARE * np.sqrt(sk.area)
        ratio = float((move / allowance).max()) if len(move) else 0.0
        if ratio > 1.0 and depth < SPLIT_DEPTH_MAX:         # too far for one part: run it again as shorter parts
            self._restore(saved, ctx)
            k = int(math.ceil(ratio))
            for _ in range(k):
                self._part(dt / k, dt_macro, depth + 1, ctx, tally)
            return
        # the accepted part's books
        tally.parts += 1
        tally.passes = max(tally.passes, passes)
        tally.ratio = max(tally.ratio, ratio * dt_macro / dt)
        tally.E_top += float(book.E_top.sum())
        tally.E_rad += float(book.E_rad.sum())
        tally.removed += float(book.removed.sum())
        tally.removed_H += float(book.removed_H.sum())
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
        w = np.abs(book.E_base)
        share = w / w.sum() if w.sum() > 0.0 else np.full(len(w), 1.0 / max(len(w), 1))
        self.pending_x, self.pending_J = base.copy(), mismatch + residual * share
        # beta moved by the draws and the hand-downs, the move's enthalpy against what the skins took and gave
        dV = (handed - book.drawn) / self.rho
        if np.any(dV != 0.0):
            expected = float(handed_H.sum() - book.drawn_H.sum())
            removed = self._move_beta(base, dV, dV / A_base, tally)
            self._book_move(base, np.abs(dV), expected + removed)

    def _book_move(self, points, weight, J):
        """Book J [J] owed to the volume after a move, at the points in proportion to `weight`."""
        if J == 0.0 or not len(points):
            return
        w = weight / weight.sum() if weight.sum() > 0.0 else np.full(len(weight), 1.0 / len(weight))
        keep = w > 0.0
        self.add_pending(points[keep], J * w[keep])

    def _move_beta(self, points, dV, disp, tally):
        """Move beta by the volume changes dV [m3] at the points [m] (k, 3), each displaced along its normal by disp [m]
        (outward positive), extended to the nodes within fine_margin of beta from the nearest point; reset, corrected per
        piece to the exact volume target, and `CutThermal.move`. Returns the enthalpy the move removed from the volume
        [J]."""
        mesh, beta = self.mesh, self.beta
        new = beta.copy()
        if np.any(disp != 0.0):
            sel = np.flatnonzero(np.abs(beta) < self.fine_margin)
            _, near = cKDTree(points).query(mesh.points[sel])
            new[sel] = beta[sel] - disp[near]
            new = redistance(mesh, new, band=self.band)
        label, n = pieces(mesh, new)
        targets = self._piece_targets(beta, new, label, n, points, dV)
        new, worst = correct_volumes(mesh, new, label, targets)
        tally.correction = max(tally.correction, worst)
        removed = self.solver.move(new)
        self.beta = new
        return removed

    def _piece_targets(self, old, new, label, n, points, dV):
        """The exact volume of each piece of `new` [m3]: the old metal in its cells (the old metal left in no piece goes
        to the largest) plus the changes dV at the points that lie in it (the largest's where outside every piece)."""
        mesh = self.mesh
        v_old = inside_volumes(mesh, old)
        cell, _ = locate(mesh, points)
        where = np.maximum(label[cell], 0)
        # summed piece by piece with numpy's pairwise sums: bincount's running sum lost 1e-13 of the volume a move
        t = np.array([v_old[label == p].sum() + dV[where == p].sum() for p in range(n)])
        t[0] += v_old[label < 0].sum()
        return t

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
            zero = type("Flow", (), {"dp_ds": np.zeros(mk.n), "tau": np.zeros(mk.n)})()
            ctx["forcing"] = Forcing.from_surface_flow(zero, np.asarray(decel_vec, float), self.v_hat,
                                                       mat.liquid.sigma, mk.kappa.sum(axis=1))
        if q_prescribed is not None:
            q = q_prescribed(mk) if callable(q_prescribed) else np.asarray(q_prescribed, float)
            ctx["q"] = np.zeros(sk.n)
            ctx["q"][mk.skin_index] = np.broadcast_to(q, (mk.n,))
        self.plan = ctx["plan"]

        # the parts (Ruling 7)
        n = self._hint
        top = self._top_T()
        if mat.T_rigid < np.inf:
            hot = np.any(np.nan_to_num(top, nan=-np.inf) >= mat.T_rigid)
            if not hot and self._top_start is not None and self._top_start.shape == top.shape:
                with np.errstate(invalid="ignore"):
                    ahead = top + (top - self._top_start) * (dt / self._dt_last)
                hot = np.any(np.nan_to_num(ahead, nan=-np.inf) >= mat.T_rigid)
            if hot:
                n = max(n, int(math.ceil(dt / PART_MAX - 1e-9)))
        self._top_start, self._dt_last = top, dt
        if n > PARTS_MAX:
            raise RuntimeError("the macro step of {:.6g} s would run as {} parts, more than PARTS_MAX = {}".format(
                dt, n, PARTS_MAX))
        tally = _Tally(sk.n)
        for _ in range(n):
            self._part(dt / n, dt, 0 if n == 1 else 1, ctx, tally)
        self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * tally.ratio)))
        self.time += dt

        # the step's end: phi, psi, the mesh, the markers
        self._rebuild(tally)
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
            time=self.time, dt=dt, parts=tally.parts, passes=tally.passes, enthalpy_in=tally.E_top,
            mass_sprayed=tally.sprayed, enthalpy_sprayed=tally.sprayed_H, mass_removed=tally.removed,
            enthalpy_removed=tally.removed_H, enthalpy_radiated=tally.E_rad, mass_volume=end["mass_volume"],
            mass_skins=end["mass_skins"], enthalpy_volume=end["enthalpy_volume"], enthalpy_skins=end["enthalpy_skins"],
            pending_J=end["pending_J"], mass_residual=(M1 - M0 + out_m) / M0,
            enthalpy_residual=(E1 - E0 - tally.E_top + out_E) / (abs(E0) + abs(tally.E_top)),
            mass_drawn=tally.drawn, mass_handed=tally.handed, momentum_sprayed=tally.sprayed_p,
            momentum_to_volume=self.momentum_to_volume.copy(), momentum_to_wall=tally.wall, n_markers=mk.n,
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
    def _rebuild(self, tally):
        """phi and psi, the mesh rebuilt from the refinement targets with the solve carried to it, then phi and psi on
        it and the markers settled (module docstring). Returns whether the mesh stayed as it was."""
        mesh = self.mesh
        phi = outer_surface(mesh, self.beta, self._thickness_nodal(), band=self.band)
        psi = self._psi_on(phi)
        target = self._target(mesh, phi, psi)
        same = not self._needs_rebuild(target, self._target(mesh, phi, psi, False, LEAD - REGION_SLACK * self.fine_size))
        if not same:
            new = build(self.base, target)
            beta = redistance(new, interpolate(mesh, self.beta, new.points), band=self.band)
            v_old = inside_volumes(mesh, self.beta)
            label, n = pieces(new, beta)
            if n == 0:
                raise RuntimeError("the skin base has no metal on the rebuilt mesh")
            targets = np.zeros(n)
            targets[0] = v_old.sum()
            if n > 1:                                      # the pieces share the old metal by their volumes here
                share = np.bincount(label[label >= 0], inside_volumes(new, beta)[label >= 0], minlength=n)
                targets = share * (v_old.sum() / share.sum())
            beta, worst = correct_volumes(new, beta, label, targets)
            tally.correction = max(tally.correction, worst)
            E0 = self.solver.energy()
            self.solver.remesh(new, beta)
            self.mesh, self.beta = new, beta
            self._book_move(self._bases(), self.skins.area, E0 - self.solver.energy())
            phi = outer_surface(new, self.beta, self._thickness_nodal(), band=self.band)
        self.phi = phi
        self.psi = self._psi_on(phi)
        self._settle()
        return same

    def _needs_rebuild(self, target, region):
        """Whether the mesh must be rebuilt (a rebuild builds every leaf to `target`): some leaf coarser than the fine
        region asks with its lead shortened by REGION_SLACK fine cells (`region`, no curvature), some leaf more than
        CURVATURE_SLACK times its full target, or more than COARSEN_SHARE of the metal (by volume) in leaves that could be
        COARSEN_SIZES times coarser. A rebuild at every step would not settle: the skin base carried to a new mesh is a
        slightly different contour, and phi's recovered curvature a slightly different field, which move leaves across
        the targets' thresholds and so change the next mesh (measured on the test cube: four rebuilds in a row, each
        with some 36 leaves over their curvature target)."""
        mesh = self.mesh
        size, x = mesh.size, mesh.centroids
        want = np.asarray(target(x, size), float)
        if np.any(size > np.asarray(region(x, size), float) * (1.0 + 1e-9)):
            return True
        if np.any(size > CURVATURE_SLACK * want):
            return True
        v = inside_volumes(mesh, self.beta)                 # the metal's leaves; closure grades the rest anyway
        coarser = COARSEN_SIZES * size <= want
        return float(v[coarser].sum() / v.sum()) > COARSEN_SHARE

    def _psi_on(self, phi):
        self.phi = phi
        return self._psi(phi)

    def _snap(self):
        mk = self.markers
        for _ in range(SNAP_ROUNDS):
            if np.abs(interpolate(self.mesh, self.phi, mk.x)).max() <= SNAP_TOL:
                return
            mk.snap(self.mesh, self.phi)
        raise RuntimeError("the markers did not snap onto phi within {} m in {} rounds".format(SNAP_TOL, SNAP_ROUNDS))

    def _release(self, idx):
        """Release the markers idx: their skins' metal and enthalpy return to the volume (beta moved by their volume,
        the move's enthalpy against theirs booked as pending at their points), their momentum to `momentum_to_volume`."""
        mk = self.markers
        x, mass, H, P = mk.release(idx)
        self.momentum_to_volume += P.sum(axis=0)
        tally = _Tally(0)
        removed = self._move_beta(x, mass / self.rho, np.zeros(len(x)), tally)
        self._book_move(x, mass, float(H.sum()) + removed)
        self.phi = outer_surface(self.mesh, self.beta, self._thickness_nodal(), band=self.band)
        self.psi = self._psi(self.phi)

    def _settle(self):
        """The markers onto phi, the trapped (and the patchless) released, tessellated, split, merged, and the skins
        regridded with their markers' areas and curvatures."""
        mk = self.markers
        self._snap()
        cell, _ = locate(self.mesh, mk.x)
        trapped = mk.trapped_markers(self.mesh, self.phi, self.mesh.size[cell])
        if len(trapped):
            self._release(trapped)
            self._snap()
        mk.tessellate(self.mesh, self.phi)
        bare = np.flatnonzero(mk.area <= 0.0)
        if len(bare):
            self._release(bare)
            self._snap()
            mk.tessellate(self.mesh, self.phi)
        self._sync()
        big = mk.split(self.spacing ** 2)
        merged = mk.merge(self.spacing)
        if len(big) or merged:
            self._snap()
            mk.tessellate(self.mesh, self.phi)
        self.skins.regrid()
        self._sync()
        self.mode = np.full(mk.n, MODES[0])          # every skin full-depth in Task 9; Task 10 sets and carries the modes
