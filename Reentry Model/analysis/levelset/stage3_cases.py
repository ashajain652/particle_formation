"""Level-set stage 3 Task 11: the verification cases and the cost of the skins on the level set (design sections 20.11
item 3 and 20.13; plan docs/superpowers/plans/levelset/03-stage3-skins.md, Task 11). Run from `Reentry Model/`:

    "$PY" analysis/levelset/stage3_cases.py plate-cold [--preset dev|production]
    "$PY" analysis/levelset/stage3_cases.py plate-hot  [--preset dev|production]
    "$PY" analysis/levelset/stage3_cases.py box        [--preset dev]
    "$PY" analysis/levelset/stage3_cases.py shoulder   [--preset dev]
    "$PY" analysis/levelset/stage3_cases.py cost --what sphere [--preset dev|production] [--steps 2] [--dt 0.5]
    "$PY" analysis/levelset/stage3_cases.py cost --what rim
    "$PY" analysis/levelset/stage3_cases.py cost --what skip
    "$PY" analysis/levelset/stage3_cases.py record

Each case prints one line per macro step (with STAGE3_STATE_DIR set in the environment it also pickles the body after
every step there, for diagnosis) and writes its summary to output/levelset/stage3_<case>[_<preset>].json
(`reentry_model.OUTPUT_DIR`); `record` gathers those into analysis/levelset/stage3_recorded.json, the committed record.
The tests are tests/levelset/test_levelset_stage3_cases.py (marked `reference`). Every case uses the shared AA7075 Scheil
material and the level set's own modules only (no skyscrapers module, no old alias module).

The plates (Decision A of stage 2: judged at 0.5 s macro steps here, where the skins hold and release the liquid within
their sub-steps). The level-2 plate's centre column: 8 mm of metal, 5.6 mm square (two base cells, stage 2 fact 8), its
walls the mesh box's own faces (on grid lines; the volume solve's natural boundary condition makes them insulated), no
radiation, `swell=False`, no flow, `remove="feed"`; heated on its top by projection (`projected_heat`). The face
carries ONE marker (the body's marker spacing set to the column's width), because the markers' tessellation assumes a
closed surface: on an open surface ending at the mesh box, a marker at the boundary has no neighbour beyond it and its
Voronoi cell runs out to the cap, and the one factor scaling the areas to phi's then shrinks every other marker's.
Measured on this column with nine markers (once, dev cells): groups of patches 0.50 to 1.24 of their skins' areas, base
temperatures 21 K apart after the first melting part, and the face no longer flat after 1.0 s (beta from 1.29 to 2.02 mm
down). A closed column (walls
in the level set, with skins) fails the other way: its convex edges melt first (5.6 mm wide, 0.35 mm-equivalent of metal
removed by 0.5 s and 1.90 by 1.0 s, against 0.79 for the one-dimensional reference). With one marker the tessellation is
exact (its cell is the whole face) and the column is the one-dimensional problem. Judged: the cold plate (700 K,
6 MW/m2) within RECESSION_TOL of REFERENCE_1D and of REFERENCE_1D_LIQUIDUS at 1.0, 1.5 and 2.0 s; the hot plate (880 K,
2 MW/m2) within RECESSION_TOL of plate_hot_1d.json at 1.5 and 2.0 s (1.0 s reported), its non-rigid depth (phi to
psi) within NONRIGID_TOL and the volume's top beneath the skin at most HOT_BASE_TOL_K above T_rigid at all three; every step's books within MASS_TOL and ENERGY_TOL. The recession is the face's
solid top (beta's zero on the axis plus the skin's mass over the solid density and its area), as the skyscrapers' plate
reads it; the recession by the removed mass is reported beside it.

The q6p box (Ruling 5): 8 x 10 x 10 mm at 700 K, 6 MW/m2 by projection on its 10 x 10 face (flight along -x; each
marker the mean projection of the patches of phi nearest it, `projected_heat`, whose docstring says why), 2 s at
0.5 s steps, dev cells, the design's meshes (0.7 mm fine cells over a 4 mm margin, a 0.25 mm curvature floor), the
volume flow held (stage 3) and the skins' own flow on (the runoff driven by the flight's deceleration, PUSH_DECEL along
the flight direction, with no gas pressure gradient or shear), `remove="feed"`. Judged: the face centre's recession
within RECESSION_TOL of REFERENCE_1D at 1.0, 1.5 and 2.0 s (the liquidus reference reported), the books, no marker
trapped; the edges reported (their recession along the bisector of the heated face's edges, the markers turned
top-zone).

The shoulder: a segment of a sphere of radius SHOULDER_R (the 50 mm sphere, as the marker tests' "25 mm sphere") about
its 45 degree station, cut by flat faces (a closed body: SHOULDER_HALF_X along the meridian, SHOULDER_HALF_Y across,
SHOULDER_DEPTH deep), at the skyscrapers' pushed plate's loads (`PUSHES["flight"]` of analysis/skin_plate.py: dp/ds
and shear down the meridian on its curved face, the deceleration along the flight direction), 700 K, 6 MW/m2 projected
on the curved face upstream of the station (the windward cap; the heat ends at the station) and nothing behind, the
liquid kept, dev cells, fine cells to the design's margin but no curvature floor below the fine size (the cut faces'
edges are not the case; the sphere's own curvature asks for none). The deposit is the solid height gained past the
heat's end (swell off: phi's radial rise over the sphere), by marker; its centroid past the heat's end is set against
the skyscrapers' pushed plate's 1.05 mm (SKYSCRAPERS_DEPOSIT_CENTROID_MM). Judged by its test: the books, a deposit
formed past the heat's end, at least DEPOSIT_RIGID_SHARE of it rigid at 2.0 s and at most the pushed plate's
3.0 per cent of it past DEPOSIT_FAR (5 mm) (the controller's ruling on the review's I2); the farthest reach of non-rigid
metal is reported, not judged.

The cost (Step 4): the 50 mm sphere at a melting state (880 K, the continuum gas state of the body tests, the shared
heating, surface flow and Girin's spray), dev and production cells, each macro step's wall time split by the profiler
below (loads and the shadow test, the skins' trial, real and repeated passes, the runoff within them, the cut-cell
solve, phi and psi, the markers, the re-partition, the remesh, the base moves with their redistance, volume correction
and rounds, the mode changes) and the peak memory; the cost of cells at the 0.175 mm floor on a rim of radius 0.5 mm
(design 20.13: a rounded puck, its mesh and one thermal step against the same puck with no curvature floor); and
whether psi and the remesh can run every second step (the cold plate's recession with them skipped at every second
step end; the shadow test does not act on the plate, which is heated by prescription, and is timed on the sphere).

Carried measurements (earlier tasks' parked items), per macro step in every case: the patch groups' area ratios; the
settling (phi rebuilt from the settled skins over phi's area); how far phi moved over the step at the markers' old
positions; the area re-partitioned at the step end and the stretch; phi against a distance function (the share of its
cut cells with |grad phi| below GRAD_SHARE_BOUND, and the angle between phi's normal at each marker and beta's at its
base); phi's fall where skins changed mode against the depth they handed; the owed heat against the heat in, and the mode
changes' owed heat against what they handed; the mass drawn and returned at the bases and beta's pieces.

Every number here is a setting of the cases or a reference value with its source; the physics is the shared package's."""
import argparse
import glob
import json
import math
import os
import pickle
import resource
import sys
import time
from collections import defaultdict
from types import SimpleNamespace

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from reentry_model import OUTPUT_DIR  # noqa: E402
from reentry_model.levelset import loads as loads_mod  # noqa: E402
from reentry_model.levelset import markers as markers_mod  # noqa: E402
from reentry_model.levelset import modes, skinbody, stepend  # noqa: E402
from reentry_model.levelset.boxmesh import BoxBase, build  # noqa: E402
from reentry_model.levelset.criteria import kappa_max_nodal, make_target, nodal_average  # noqa: E402
from reentry_model.levelset.cutcell import cell_gradients, patches  # noqa: E402
from reentry_model.levelset.cutthermal import CutThermal  # noqa: E402
from reentry_model.levelset.levelset import redistance  # noqa: E402
from reentry_model.levelset.skin import SKIN_PRESETS, SkinSettings  # noqa: E402
from reentry_model.levelset.surfaces import outer_surface  # noqa: E402
from reentry_model.levelset.transfer import interpolate, locate  # noqa: E402
from reentry_model.shared import atmosphere, heating, spray, surface_flow  # noqa: E402
from reentry_model.shared import trajectory as tj  # noqa: E402
from reentry_model.shared.material import LiquidProperties, Material  # noqa: E402

OUT = os.path.join(OUTPUT_DIR, "levelset")
RECORD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stage3_recorded.json")

# -- references (copied as data, with their sources) -------------------------------------------------------------------
REFERENCE_1D = {1.0: 0.789, 1.5: 2.712, 2.0: 5.219}     # mm: the cold plate's recession, analysis/skin_plate.py 86-95
                                                         # (analysis/ablation_1d.py at 1 um and 0.5 ms, removal at
                                                         # T_feed, 2026-10-12); skin_plate.py imports the skyscrapers'
                                                         # skins, so it is read, not imported
REFERENCE_1D_LIQUIDUS = {1.0: 0.837, 1.5: 2.772, 2.0: 5.298}   # mm: the same with the removal at the liquidus (2026-10-24)
HOT_REFERENCE = "plate_hot_1d.json"                      # data/reference_values: the hot plate's recession and non-rigid
                                                         # depth (ablation_1d.py --t0 880 --q 2e6, 1 um, 0.5 ms)
PUSH_ALPHA, PUSH_DP_DS, PUSH_TAU, PUSH_DECEL = 45.0, -3.424e5, 80.29, 29.77   # skin_plate.PUSHES["flight"]: the 100 mm
                                                         # sphere's shoulder 45 degrees from the stagnation point at 93 s:
                                                         # dp/ds [Pa/m], gas shear [Pa], deceleration [m/s2]
SKYSCRAPERS_DEPOSIT_CENTROID_MM = 1.05                   # the pushed plate's deposit centroid past the heat's end at
                                                         # 10 ms heat sub-steps with the landing conduction (sub-plan 21,
                                                         # docs/superpowers/plans/melt-spraying-subplans/21-column-flow.md
                                                         # item 21: 1.84 mm without it)
SKYSCRAPERS_DEPOSIT_FAR_SHARE = 0.030                    # its share landing 5-10 mm past the heat's end (9.9 % without)
HOT_BASE_TOL_K = 1.0                                     # the skyscrapers' hot plate: base at most this above T_rigid
                                                         # (skin_plate.HOT_BASE_TOL_K; judged since Task 11's fix, item 4)

# -- tolerances (the brief's) ------------------------------------------------------------------------------------------
RECESSION_TOL = 0.05
NONRIGID_TOL = 0.10
MASS_TOL = 1e-12
ENERGY_TOL = 1e-8
JUDGED_TIMES = (1.0, 1.5, 2.0)

# -- the cases' settings -----------------------------------------------------------------------------------------------
H0 = 2.8e-3                        # m: the base cells (stage 2's plates, the body tests)
DT, T_END = 0.5, 2.0               # s: the macro step and the cases' end
PLATE_CASES = {"cold": (700.0, 6.0e6), "hot": (880.0, 2.0e6)}    # K, W/m2 (the level-2 plate and the hot plate)
PLATE_WIDTH, PLATE_DEPTH, PLATE_GAS = 5.6e-3, 8e-3, 3.2e-3       # m: the column, its metal, the gas above its face
BOX_LENGTH, BOX_HALF, BOX_OFF = 8e-3, 5e-3, 0.13e-3              # m: the q6p box; its offset off the grid
BOX_Q = 6.0e6                      # W/m2
FACE_NORMAL = 0.8                  # a marker whose normal is within acos(0.8) = 37 degrees of -x is on the box's heated face
COLD_T0 = 700.0                    # K
SHOULDER_R = 25e-3                 # m: the sphere's radius
SHOULDER_HALF_X, SHOULDER_HALF_Y, SHOULDER_DEPTH = 7e-3, 3e-3, 6e-3   # m: the segment about the station (local frame:
                                                                    # station at the origin, normal +z, meridian +x)
SHOULDER_Q = 6.0e6                 # W/m2, upstream of the station (x < 0)
CURVED_FACE = 0.5                  # a marker whose normal's dot with the radial direction exceeds this is on the sphere
DEPOSIT_RIGID_SHARE = 0.9          # of the metal deposited past the heat's end at 2.0 s, at least this share rigid
DEPOSIT_FAR = 5e-3                 # m: of the deposit, at most SKYSCRAPERS_DEPOSIT_FAR_SHARE past this from the heat's end
                                   # (the controller's ruling on the review's I2, 2026-11-13)
SPHERE_R = 25e-3                   # m: the cost case's 50 mm sphere
SPHERE_T0 = 880.0                  # K: just below T_rigid (895.1 K), so that the windward cap melts within the first step
CONTINUUM = dict(rho=1e-4, kn=1e-3)   # the body tests' continuum gas state (tests/levelset/test_levelset_skinbody.py)
RIM_R, RIM_HALF_Z, RIM_EDGE = 12.5e-3, 3e-3, 0.5e-3   # m: the puck of the rim cost, its half height, the rim's radius
GRAD_SHARE_BOUND = 0.9             # |grad phi| below this in a cut cell: phi far from a distance function


def material():
    return Material.from_drama_json("AA7075_scheil")


def skin_settings(preset, swell=False):
    return SkinSettings(cell=SKIN_PRESETS[preset], swell=swell)


def gas_state(rho=2.727e-5, T=203.6, V=7500.0, kn=0.0298, ma=26.2):
    """An aero state at 77.5 km and 7.5 km/s (tests/levelset/cases.py's `gas_state`, copied: the analysis does not import
    the tests): the freestream's density, temperature, speed, Knudsen and Mach numbers, the deceleration a thousandth of
    the speed per second."""
    fs = atmosphere.Freestream(rho, T, rho / 28.9644e-3 * 8.314462618 * T, 28.9644 * 1.66053906660e-27, np.zeros(3))
    return tj.AeroState(77.5e3, 0.0, 0.0, fs, np.array([V, 0.0, 0.0]), V, kn, ma, 0.9, np.array([-V * 1e-3, 0, 0]),
                        0.5 * rho * V * V)


def peak_rss_gb():
    """The process's peak resident memory so far [GB] (ru_maxrss is in bytes on macOS, kilobytes on Linux)."""
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == "darwin" else r * 1024 / 1e9


# -- the profiler ------------------------------------------------------------------------------------------------------
class Profile:
    """Wall time by label over the body's calls, inclusive (outermost call of a label only, so a recursion counts once)
    and exclusive (less the labelled calls inside), with call counts; installed by wrapping the methods and module
    functions named in `install`, restored on exit (a context manager). `attempt` counts the passes of the current part,
    so that the skins' first trial and real passes are told from their repeats."""

    def __init__(self):
        self.incl, self.excl, self.calls = defaultdict(float), defaultdict(float), defaultdict(int)
        self.stack, self.depth, self.saved = [], defaultdict(int), []
        self.attempt = 0

    def _wrap(self, owner, name, label):
        orig = getattr(owner, name)
        prof = self

        def wrapper(*a, **k):
            lab = label(*a, **k) if callable(label) else label
            prof.stack.append(0.0)
            prof.depth[lab] += 1
            t0 = time.perf_counter()
            try:
                return orig(*a, **k)
            finally:
                el = time.perf_counter() - t0
                child = prof.stack.pop()
                if prof.depth[lab] == 1:
                    prof.incl[lab] += el
                prof.depth[lab] -= 1
                prof.excl[lab] += el - child
                prof.calls[lab] += 1
                if prof.stack:
                    prof.stack[-1] += el
        self.saved.append((owner, name, orig))
        setattr(owner, name, wrapper)

    def _pass_label(self, body, dt, T_b, g, m_min, ctx, trial=False):
        if trial:
            self.attempt += 1
        kind = "trial" if trial else "real"
        return "skins {}".format(kind) if self.attempt <= 1 else "skins repeats"

    def _within(self, what):
        """A label for a geometry call made both by the base moves and by the remesh: by whichever is running."""
        return lambda *a, **k: "base moves " + what if self.depth["base moves"] else "remesh " + what

    def _solve_label(self, *a, **k):
        return "cut-cell solve" if self.attempt <= 1 else "cut-cell solve repeats"

    def install(self):
        B, S = skinbody.LevelSetSkinBody, stepend.StepEnd
        orig_part = B._part_once
        prof = self

        def part_once(body, *a, **k):
            prof.attempt = 0
            return orig_part(body, *a, **k)
        self.saved.append((B, "_part_once", orig_part))
        B._part_once = part_once
        self._wrap(B, "_pass", self._pass_label)
        self._wrap(CutThermal, "step", self._solve_label)
        for owner, name, label in (
                (B, "_interface", "interface (patch groups)"), (B, "_accept", "part books"),
                (B, "mass_centroid", "loads"), (B, "_summary", "summary"),
                (S, "_move_beta", "base moves"), (S, "_psi", "surfaces psi"), (S, "_settle", "markers"),
                (S, "_repartition", "re-partition"), (S, "_modes", "mode changes"), (S, "_needs_rebuild", "remesh"),
                (S, "_target", "remesh"), (S, "_fill", "fill"), (S, "_release", "markers"),
                (CutThermal, "move", "base moves volume"), (CutThermal, "remesh", "remesh carry"),
                (markers_mod.Markers, "tessellate", "markers tessellate"), (markers_mod.Markers, "snap", "markers snap"),
                (markers_mod.Markers, "split", "markers"), (markers_mod.Markers, "merge", "markers"),
                (markers_mod.Markers, "trapped_markers", "markers"),
                (stepend, "outer_surface", "surfaces phi"), (skinbody, "outer_surface", "surfaces phi"),
                (stepend, "build", "remesh build"), (stepend, "redistance", self._within("redistance")),
                (stepend, "correct_volumes", self._within("volume correction")), (stepend, "pieces", self._within("pieces")),
                (skinbody, "pieces", "summary"), (skinbody, "cut_patches", "summary"),
                (skinbody, "flow_substeps", "runoff"), (skinbody, "spray_plan", "loads spray plan"),
                (loads_mod, "leeward", "loads shadow"), (loads_mod, "leeward_patches", "loads shadow"),
                (loads_mod, "by_patches", "loads"), (loads_mod, "fit_nose", "loads"),
                (loads_mod, "marker_angles", "loads"), (loads_mod.GasLoads, "evaluate", "loads")):
            self._wrap(owner, name, label)
        return self

    def __enter__(self):
        return self.install()

    def __exit__(self, *exc):
        for owner, name, orig in reversed(self.saved):
            setattr(owner, name, orig)
        self.saved = []

    def snapshot(self):
        return {k: (dict(self.incl), dict(self.excl), dict(self.calls))[i] for i, k in enumerate(("incl", "excl", "calls"))}

    @staticmethod
    def since(before, after):
        """The labels' times and calls between two snapshots."""
        out = {}
        for key in ("incl", "excl", "calls"):
            a, b = before[key], after[key]
            out[key] = {k: b[k] - a.get(k, 0) for k in b if b[k] - a.get(k, 0)}
        return out


# the brief's split of a macro step, each the exclusive time of these labels (they sum to the step less "other")
SPLIT = {"loads": ("loads", "loads shadow", "loads spray plan"),
         "skins trial": ("skins trial",), "skins real": ("skins real",), "skins repeats": ("skins repeats",),
         "runoff (within the skins' passes)": ("runoff",),
         "cut-cell solve": ("cut-cell solve",), "cut-cell solve repeats": ("cut-cell solve repeats",),
         "interface and part books": ("interface (patch groups)", "part books"),
         "surfaces phi": ("surfaces phi",), "surfaces psi": ("surfaces psi",),
         "markers": ("markers", "markers tessellate", "markers snap"), "re-partition": ("re-partition",),
         "remesh": ("remesh", "remesh build", "remesh carry", "remesh redistance", "remesh volume correction",
                    "remesh pieces"),
         "base moves": ("base moves", "base moves volume", "base moves redistance", "base moves volume correction",
                        "base moves pieces", "fill"),
         "mode changes (own)": ("mode changes",), "step summary": ("summary",)}


def split_step(delta, wall):
    """A step's wall time [s] by the brief's split (exclusive times, SPLIT) with the rest as `other`, and the inclusive
    times of the base moves, the mode changes, the re-partition, the shadow test and the runoff, with the calls."""
    ex, inc, calls = delta["excl"], delta["incl"], delta["calls"]
    out = {k: round(sum(ex.get(lab, 0.0) for lab in labs), 3) for k, labs in SPLIT.items()}
    out["other"] = round(wall - sum(out.values()), 3)
    out["inclusive"] = {k: round(inc.get(k, 0.0), 3) for k in ("base moves", "mode changes", "re-partition",
                                                               "loads shadow", "runoff", "markers", "base moves volume",
                                                               "base moves redistance", "base moves volume correction")}
    out["calls"] = {k: calls.get(k, 0) for k in ("base moves", "base moves volume", "cut-cell solve",
                                                 "cut-cell solve repeats", "skins trial", "skins repeats", "runoff",
                                                 "re-partition", "remesh build", "mode changes")}
    return out


# -- measures on the body ----------------------------------------------------------------------------------------------
def unit_normals(mesh, f, x):
    """The unit recovered gradient of the field f at the points x (n, 3) (the markers' normal: the nodal average of the
    cells' gradients, interpolated)."""
    cell, b = locate(mesh, x)
    near = mesh.ring(cell, 1)
    G = nodal_average(mesh, cell_gradients(mesh, f, near), near)
    g = np.einsum("ni,nik->nk", b, G[mesh.tets[cell]])
    return g / np.maximum(np.linalg.norm(g, axis=1), 1e-300)[:, None]


def zero_along(mesh, f, start, direction, length, n=4001):
    """The distance [m] from `start` along the unit `direction` to the first point where f falls below zero (linear
    between samples), or nan."""
    s = np.linspace(0.0, length, n)
    pts = np.asarray(start, float)[None, :] + s[:, None] * np.asarray(direction, float)[None, :]
    v = interpolate(mesh, f, pts)
    k = np.flatnonzero(v < 0.0)
    if not len(k) or k[0] == 0:
        return float("nan")
    k = k[0]
    return float(s[k - 1] + (s[k] - s[k - 1]) * v[k - 1] / (v[k - 1] - v[k]))


def solid_thickness(body, rows):
    """The skins' rows' metal as solid [m]: mass over the solid density and the top area."""
    sk = body.skins
    return sk.m[rows].sum(axis=1) / (sk.rho_s * sk.area[rows])


def phi_far_from_distance(body):
    """(share of phi's cut cells with |grad phi| below GRAD_SHARE_BOUND, median, 95th percentile and largest angle [deg]
    between phi's normal at each marker and beta's normal at its base point)."""
    mesh, phi = body.mesh, body.phi
    p = phi[mesh.tets]
    cut = np.flatnonzero((p.min(axis=1) < 0.0) & (p.max(axis=1) > 0.0))
    g = np.linalg.norm(cell_gradients(mesh, phi, cut), axis=1)
    share = float((g < GRAD_SHARE_BOUND).mean()) if len(cut) else 0.0
    mk = body.markers
    base = body._bases()[mk.skin_index]
    n_phi = unit_normals(mesh, phi, mk.x)
    n_beta = unit_normals(mesh, body.beta, base)
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n_phi, n_beta), -1.0, 1.0)))
    return share, float(np.median(ang)), float(np.percentile(ang, 95)), float(ang.max())


def settling_ratio(body):
    """phi rebuilt from the settled skins (each node the thickness of the nearest base, as the step end builds it) over
    the area of the phi they were settled on (Task 9b's settling)."""
    rebuilt = outer_surface(body.mesh, body.beta, body._thickness_nodal(), band=body.band)
    return float(patches(body.mesh, rebuilt).area.sum() / patches(body.mesh, body.phi).area.sum())


class ModeWatch:
    """phi's fall where markers change mode at a step end (the density limitation, Task 10): wraps the body's `_modes`;
    for each step end that changed modes and kept the markers' count, the changed markers' phi after the change at their
    positions before it (phi is a distance near the surface: its value is the fall, positive inward) against the depth
    their skins handed (thickness before less after, positive to the volume)."""

    def __init__(self, body):
        self.events = []
        orig = body._modes

        def watched(tally, update=True, check_core=False):
            mk, sk = body.markers, body.skins
            x, top, n = mk.x.copy(), mk.top.copy(), mk.n
            H = sk.thickness()[mk.skin_index].copy()
            out = orig(tally, update=update, check_core=check_core)
            if body.markers.n == n:
                changed = np.flatnonzero(body.markers.top != top)
                if len(changed):
                    fall = interpolate(body.mesh, body.phi, x[changed])
                    handed = H[changed] - body.skins.thickness()[body.markers.skin_index[changed]]
                    self.events.append({"changed": int(len(changed)), "to_top": int(body.markers.top[changed].sum()),
                                        "fall_mm": fall * 1e3, "handed_mm": handed * 1e3})
            return out
        body._modes = watched

    def take(self):
        out, self.events = self.events, []
        if not out:
            return None
        fall = np.concatenate([e["fall_mm"] for e in out])
        handed = np.concatenate([e["handed_mm"] for e in out])
        ok = np.abs(handed) > 1e-3
        ratio = fall[ok] / handed[ok]
        return {"changed": int(sum(e["changed"] for e in out)), "to_top": int(sum(e["to_top"] for e in out)),
                "fall_mm_median": float(np.median(fall)), "fall_mm_max": float(np.max(fall)),
                "handed_mm_median": float(np.median(handed)),
                "fall_over_handed_median": float(np.median(ratio)) if len(ratio) else None,
                "fall_over_handed_range": [float(ratio.min()), float(ratio.max())] if len(ratio) else None}


def step_row(body, book, wall, x_before, delta, watch):
    """One macro step's record: the books, the parts, the cost split and the carried measurements."""
    mk = body.markers
    share, ang_med, ang_95, ang_max = phi_far_from_distance(body)
    move = np.abs(interpolate(body.mesh, body.phi, x_before))
    pend = float(np.abs(body.pending_J).sum())
    heat = max(book.enthalpy_in, 1e-300)
    row = {"t": round(book.time, 9), "wall_s": round(wall, 2), "peak_rss_gb": round(peak_rss_gb(), 3),
           "parts": book.parts, "passes_max": book.passes, "parts_next": book.parts_next,
           "split_depth_max": book.split_depth_max, "parts_over": book.parts_over,
           "mass_residual": book.mass_residual, "enthalpy_residual": book.enthalpy_residual,
           "enthalpy_in_J": book.enthalpy_in, "enthalpy_radiated_J": book.enthalpy_radiated,
           "mass_removed_kg": book.mass_removed, "mass_sprayed_kg": book.mass_sprayed,
           "mass_volume_kg": book.mass_volume, "mass_skins_kg": book.mass_skins,
           "n_markers": book.n_markers, "markers_trapped": book.markers_trapped, "marker_splits": book.marker_splits,
           "marker_merges": book.marker_merges, "n_top": book.n_top, "modes_to_top": book.modes_to_top,
           "modes_to_full": book.modes_to_full, "cells": int(len(body.mesh.tets)),
           "cells_at_floor": int((body.mesh.size <= body.cell_floor * (1 + 1e-9)).sum()),
           "snap_fallback": book.snap_fallback, "snap_slide_max_m": book.snap_slide_max,
           "volume_correction_max": book.volume_correction_max,
           # carried 1: uneven patch groups
           "group_area_ratio_min": book.group_area_ratio_min, "group_area_ratio_max": book.group_area_ratio_max,
           "group_area_ratio_trimmed_max": book.group_area_ratio_trimmed_max,
           # carried 3: settling, and phi's move at the markers' old positions
           "settling_area_ratio": settling_ratio(body),
           "phi_move_max_m": float(move.max()), "phi_move_median_m": float(np.median(move)),
           # carried 4: sideways mixing
           "area_repartitioned_share": book.area_repartitioned_share, "stretch_min": book.stretch_min,
           "stretch_max": book.stretch_max,
           # carried 5: phi against a distance function
           "grad_phi_below_bound_share": share, "phi_beta_angle_median_deg": ang_med,
           "phi_beta_angle_p95_deg": ang_95, "phi_beta_angle_max_deg": ang_max,
           # carried 6: the density limitation where skins changed mode
           "mode_phi_fall": watch.take() if watch is not None else None,
           # carried 7: owed heat
           "owed_share_of_heat_in": max(pend - book.residual_modes_J, 0.0) / heat,
           "pending_share_of_heat_in": pend / heat, "interface_mismatch_share": book.interface_mismatch_J / heat,
           "residual_newton_J": book.residual_newton_J, "residual_move_J": book.residual_move_J,
           "residual_remesh_J": book.residual_remesh_J,
           "enthalpy_modes_J": book.enthalpy_modes_J, "residual_modes_J": book.residual_modes_J,
           "modes_owed_share": book.residual_modes_J / book.enthalpy_modes_J if book.enthalpy_modes_J else None,
           # carried 8: churn and pieces
           "mass_drawn_kg": book.mass_drawn, "mass_handed_kg": book.mass_handed,
           "mass_modes_handed_kg": book.mass_modes_handed, "mass_modes_absorbed_kg": book.mass_modes_absorbed,
           "n_pieces": book.n_pieces, "pieces_used_up": book.pieces_used_up,
           # Task 11's fix: the islands a base move filled from the handed columns (item 3), the groups the area guard
           # replaced (item 1), the non-rigid metal the full-depth skins absorbed (item 4)
           "island_nodes": getattr(book, "island_nodes", 0), "groups_guarded_max": getattr(book, "groups_guarded_max", 0),
           "mass_absorbed_kg": getattr(book, "mass_absorbed", 0.0),
           "pending_flux_max_Wm2": book.pending_flux_max, "surface_flux_max_Wm2": book.surface_flux_max}
    if delta is not None:
        row["cost_split_s"] = split_step(delta, wall)
    return row


def projected_heat(body, q0, direction, where=None):
    """A flux q0 [W/m2] arriving along -`direction` (a beam from the side `direction` points to) by projection, given to
    the markers by the patches of phi they own, as the gas loads now are (`loads.patch_mean`, stage 3 Task 11's fix, item
    2) and as the skyscrapers' plate gives each column the projection of its own dual cell (`projected_flux`): each marker
    takes q0 times the area-weighted mean, over phi's patches nearest it, of their own projection max(n_p . direction, 0),
    on its own area; `where(centroids)` (bool) heats only some patches. Returns q per marker [W/m2]. (Each marker's own
    normal, q0 max(n . direction, 0), gave the q6p box's markers 73.6 mm2 of projected area against phi's 99.3 mm2,
    measured once after its first step: its edges, rounded to the stub's 0.8 mm, lie between markers 2 mm apart.)"""
    P = patches(body.mesh, body.phi)
    proj = np.maximum(P.normal @ np.asarray(direction, float), 0.0)
    if where is not None:
        proj = np.where(where(P.centroid), proj, 0.0)
    return q0 * loads_mod.patch_mean(body.markers, P, proj)


def heat_coverage(body, direction, q=None, q0=1.0, where=None):
    """The heat the markers take over what phi's projected area would take (both along `direction`): by their own normals
    (q None), or by the fluxes q [W/m2] per marker of a beam q0 [W/m2] (`projected_heat`), over phi's patches `where`."""
    mk = body.markers
    P = patches(body.mesh, body.phi)
    d = np.asarray(direction, float)
    proj = P.area * np.maximum(P.normal @ d, 0.0)
    if where is not None:
        proj = np.where(where(P.centroid), proj, 0.0)
    took = (mk.area * np.maximum(mk.normal @ d, 0.0)).sum() if q is None else (np.asarray(q) * mk.area).sum() / q0
    return float(took / proj.sum())


def books_pass(rows):
    return all(abs(r["mass_residual"]) <= MASS_TOL and abs(r["enthalpy_residual"]) <= ENERGY_TOL for r in rows)


def write(name, out):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "stage3_{}.json".format(name))
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=_json_default)
    print("wrote", path, flush=True)
    return path


def _json_default(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    raise TypeError(type(v))


def save_state(body, path):
    """Pickle the body to `path` (large; never committed), the cases' own wrappers of its methods left out."""
    over = {k: body.__dict__.pop(k) for k in ("_modes", "_psi", "_target", "_needs_rebuild") if k in body.__dict__}
    try:
        with open(path, "wb") as fh:
            pickle.dump(body, fh)
    finally:
        body.__dict__.update(over)
    return path


def run_steps(body, advance, steps, label, extra=None, before_step=None):
    """`steps` macro steps of `advance(body)` under the profiler, a row each (`step_row`, plus extra(body, book) if given);
    a step that raises ends the run and is recorded. With STAGE3_STATE_DIR set, the body is pickled there before each
    step and kept, its path in the failure, when the step raises (Task 11's review, M6: a failing case left nothing to
    diagnose), and after each step that succeeds. Returns (rows, failure or None)."""
    rows, failure = [], None
    tag = label.replace(" ", "_")
    states = os.environ.get("STAGE3_STATE_DIR")
    with Profile() as prof:
        watch = ModeWatch(body)                            # inside the profiler, so that it wraps the profiled _modes
        for k in range(steps):
            if before_step is not None:
                before_step(k)
            x0 = body.markers.x.copy()
            t_start = float(body.time)
            before = save_state(body, os.path.join(states, "{}_before_step.pkl".format(tag))) if states else None
            snap = prof.snapshot()
            t0 = time.perf_counter()
            try:
                book = advance(body)
            except Exception as err:                       # recorded, not hidden: the run stops here
                failure = {"step": k + 1, "t_start": round(t_start, 9), "error": "{}: {}".format(type(err).__name__, err)}
                if before is not None:                     # the body the failing step started from, kept
                    failure["state"] = os.path.join(states, "{}_failed_from_t{:.2f}.pkl".format(tag, t_start))
                    os.replace(before, failure["state"])
                print(label, "FAILED:", failure, flush=True)
                break
            wall = time.perf_counter() - t0
            row = step_row(body, book, wall, x0, Profile.since(snap, prof.snapshot()), watch)
            if extra is not None:
                row.update(extra(body, book))
            rows.append(row)
            os.makedirs(OUT, exist_ok=True)                  # the steps so far, in case the run is stopped
            with open(os.path.join(OUT, "partial_{}.json".format(tag)), "w") as fh:
                json.dump({"steps": rows, "failure": failure}, fh, indent=1, default=_json_default)
            if states:                                       # each step's body kept for diagnosis
                save_state(body, os.path.join(states, "{}_t{:.2f}.pkl".format(tag, row["t"])))
            print("{} t={:.2f}s wall {:.1f}s parts {} mass {:.1e} energy {:.1e} {}".format(
                label, row["t"], wall, row["parts"], row["mass_residual"], row["enthalpy_residual"],
                {k: v for k, v in row.items() if k.endswith("_mm") or k in ("n_top", "n_markers")}), flush=True)
    return rows, failure


# -- the plates --------------------------------------------------------------------------------------------------------
def hot_reference():
    path = os.path.join(REPO_ROOT, "data", "reference_values", HOT_REFERENCE)
    with open(path) as fh:
        doc = json.load(fh)
    return {t: (float(np.interp(t, doc["t_s"], doc["recession_mm"])), float(np.interp(t, doc["t_s"], doc["nonrigid_mm"])))
            for t in JUDGED_TIMES}


def plate_body(case, preset, mat=None):
    """The plate's centre column (module docstring) at the case's starting temperature, its face one marker."""
    mat = material() if mat is None else mat
    T0, _ = PLATE_CASES[case]
    st = skin_settings(preset)
    w = PLATE_WIDTH
    base = BoxBase.around((-w / 2, -w / 2, -PLATE_DEPTH), (w / 2, w / 2, PLATE_GAS), H0)
    body = skinbody.LevelSetSkinBody(mat, base, lambda x: x[:, 2] + st.stub, T0, st, remove="feed", spacing=w)
    if body.markers.n != 1:
        raise RuntimeError("the plate's face should carry one marker, got {}".format(body.markers.n))
    return body


def plate_top(body):
    """The face's solid top [m] on the column's axis: beta's zero there plus the skin's metal as solid."""
    zb = PLATE_GAS - zero_along(body.mesh, body.beta, (0.0, 0.0, PLATE_GAS - 1e-6), (0.0, 0.0, -1.0),
                                PLATE_GAS + PLATE_DEPTH - 2e-6) - 1e-6
    return zb + float(solid_thickness(body, body.markers.skin_index)[0])


class EverySecond:
    """psi and the remesh run only at every second step end (`skip` set by the caller on the others): the body's `_psi`
    gives the last step end's psi while the mesh is the same, and `_needs_rebuild` says no (`_target` is not evaluated).
    The shadow test needs a gas state and does not act on the plate."""

    def __init__(self, body):
        self.skip = False
        cls = type(body)

        def psi(phi):
            if self.skip and body.psi is not None and len(body.psi) == len(body.mesh.points):
                return body.psi
            return cls._psi(body, phi)

        def target(*a, **k):
            return (lambda x, s: np.full(len(x), np.inf)) if self.skip else cls._target(body, *a, **k)

        def needs(target_, region):
            return False if self.skip else cls._needs_rebuild(body, target_, region)
        body._psi, body._target, body._needs_rebuild = psi, target, needs


def plate(case, preset="dev", every_second=False, steps=None, save=True):
    """The plate case (module docstring). Returns its summary (and writes it)."""
    mat = material()
    T0, q0 = PLATE_CASES[case]
    t_build = time.perf_counter()
    body = plate_body(case, preset, mat)
    t_build = time.perf_counter() - t_build
    area = float(body.markers.area.sum())
    top0 = plate_top(body)
    removed = [0.0]
    skipper = EverySecond(body) if every_second else None

    def advance(b):
        return b.advance(DT, q_prescribed=lambda mk: projected_heat(b, q0, (0.0, 0.0, 1.0)))

    def extra(b, book):
        removed[0] += book.mass_removed
        mk = b.markers
        T_patch = b.solver.patch_T()
        return {"recession_mm": (top0 - plate_top(b)) * 1e3,
                "removed_recession_mm": removed[0] / (mat.rho * area) * 1e3,
                "nonrigid_mm": float(modes.melt_depth(b.mesh, b.psi, mk.x)[0] * 1e3),
                "skin_flowing_mm": float(b.skins.depth_above(mat.T_rigid)[mk.skin_index][0] * 1e3),
                "skin_thickness_mm": float(b.skins.thickness()[mk.skin_index][0] * 1e3),
                "base_excess_K": float((T_patch * b.solver.P.area).sum() / b.solver.P.area.sum() - mat.T_rigid),
                "base_excess_max_K": float(T_patch.max() - mat.T_rigid)}

    def before(k):
        if skipper is not None:
            skipper.skip = (k % 2 == 1)                    # the step ends at 1.0 and 2.0 s skip psi and the remesh

    n = int(round(T_END / DT)) if steps is None else steps
    label = "plate-{} {}{}".format(case, preset, " every-second" if every_second else "")
    rows, failure = run_steps(body, advance, n, label, extra, before)
    at = {round(r["t"], 6): r for r in rows}
    compare, verdict = {}, []
    if case == "cold":
        for t in JUDGED_TIMES:
            if t not in at:
                continue
            r = at[t]["recession_mm"]
            rel, rel_l = r / REFERENCE_1D[t] - 1.0, r / REFERENCE_1D_LIQUIDUS[t] - 1.0
            compare[str(t)] = {"recession_mm": r, "removed_recession_mm": at[t]["removed_recession_mm"],
                               "ref_1d_mm": REFERENCE_1D[t], "rel_1d": rel, "ref_liquidus_mm": REFERENCE_1D_LIQUIDUS[t],
                               "rel_liquidus": rel_l, "nonrigid_mm": at[t]["nonrigid_mm"]}
            verdict.append(abs(rel) <= RECESSION_TOL and abs(rel_l) <= RECESSION_TOL)
    else:
        ref = hot_reference()
        for t in JUDGED_TIMES:
            if t not in at:
                continue
            r, nr = at[t]["recession_mm"], at[t]["nonrigid_mm"]
            rel, rel_nr = r / ref[t][0] - 1.0, nr / ref[t][1] - 1.0
            compare[str(t)] = {"recession_mm": r, "removed_recession_mm": at[t]["removed_recession_mm"],
                               "ref_mm": ref[t][0], "rel": rel, "judged_recession": t > 1.25, "nonrigid_mm": nr,
                               "nonrigid_ref_mm": ref[t][1], "nonrigid_rel": rel_nr,
                               "skin_flowing_mm": at[t]["skin_flowing_mm"], "base_excess_K": at[t]["base_excess_K"],
                               "base_within_skyscrapers_tol": at[t]["base_excess_K"] <= HOT_BASE_TOL_K}
            verdict.append((t < 1.25 or abs(rel) <= RECESSION_TOL) and abs(rel_nr) <= NONRIGID_TOL
                           and at[t]["base_excess_K"] <= HOT_BASE_TOL_K)
    out = {"case": "plate-" + case, "preset": preset, "every_second": every_second, "T0_K": T0, "q_Wm2": q0,
           "dt_s": DT, "width_mm": PLATE_WIDTH * 1e3, "depth_mm": PLATE_DEPTH * 1e3, "markers": int(body.markers.n),
           "build_s": round(t_build, 2), "books_pass": books_pass(rows), "compare": compare,
           "passes": bool(verdict) and len(verdict) == len(JUDGED_TIMES) and all(verdict) and books_pass(rows)
           and failure is None, "failure": failure, "steps": rows}
    if save:
        write("plate_{}_{}{}".format(case, preset, "_every_second" if every_second else ""), out)
    return out


# -- the q6p box -------------------------------------------------------------------------------------------------------
def box_beta0(stub):
    def beta0(x):
        y, z = x[:, 1] - BOX_OFF, x[:, 2] - BOX_OFF
        return np.maximum.reduce([-(x[:, 0] - BOX_OFF), x[:, 0] - BOX_OFF - BOX_LENGTH, np.abs(y) - BOX_HALF,
                                  np.abs(z) - BOX_HALF]) + stub
    return beta0


def box_body(preset="dev", mat=None):
    mat = material() if mat is None else mat
    st = skin_settings(preset)
    base = BoxBase.around((-H0, -BOX_HALF - H0, -BOX_HALF - H0), (BOX_LENGTH + H0, BOX_HALF + H0, BOX_HALF + H0), H0)
    return skinbody.LevelSetSkinBody(mat, base, box_beta0(st.stub), COLD_T0, st, remove="feed", v_hat=(-1.0, 0.0, 0.0))


def box_face(body):
    """The heated face: the face markers' solid tops along -x [m] (their bases moved out by their metal as solid), their
    distances from the face's axis [m], and the edge's recession along the bisector of the face's +y edge [m]."""
    mk = body.markers
    face = mk.normal[:, 0] < -FACE_NORMAL
    rows = mk.skin_index[face]
    base = body._bases()[rows]
    top = base[:, 0] - solid_thickness(body, rows)
    d = np.hypot(mk.x[face, 1] - BOX_OFF, mk.x[face, 2] - BOX_OFF)
    corner = np.array([BOX_OFF, BOX_OFF + BOX_HALF, BOX_OFF])
    u = np.array([1.0, -1.0, 0.0]) / math.sqrt(2.0)
    edge = zero_along(body.mesh, body.phi, corner - 2e-3 * u, u, 6e-3) - 2e-3
    return top, d, edge


def box_centre(body):
    """The face centre's solid top [m] (beta's zero on the axis plus the skin of the face marker nearest the axis) and
    that marker's distance from the axis [m]."""
    mk = body.markers
    face = np.flatnonzero(mk.normal[:, 0] < -FACE_NORMAL)
    d = np.hypot(mk.x[face, 1] - BOX_OFF, mk.x[face, 2] - BOX_OFF)
    i = face[np.argmin(d)]
    xb = -H0 + zero_along(body.mesh, body.beta, (-H0, BOX_OFF, BOX_OFF), (1.0, 0.0, 0.0), H0 + BOX_LENGTH)
    return xb - float(solid_thickness(body, mk.skin_index[[i]])[0]), float(d.min())


def box(preset="dev", steps=None, save=True):
    """The q6p box (module docstring)."""
    mat = material()
    t_build = time.perf_counter()
    body = box_body(preset, mat)
    t_build = time.perf_counter() - t_build
    v_hat = body.v_hat
    c0, _ = box_centre(body)
    top0, d0, edge0 = box_face(body)
    face_x0 = float(np.median(top0))
    markers0, coverage0 = int(body.markers.n), heat_coverage(body, v_hat)
    cells0 = {"cells": int(len(body.mesh.tets)),
              "by_size_mm": {str(round(s * 1e3, 4)): int(c) for s, c in zip(*np.unique(body.mesh.size, return_counts=True))}}

    def advance(b):
        return b.advance(DT, decel_vec=PUSH_DECEL * v_hat, q_prescribed=lambda mk: projected_heat(b, BOX_Q, v_hat))

    def extra(b, book):
        c, d = box_centre(b)
        top, dist, edge = box_face(b)
        rec = (top - face_x0) * 1e3
        mk = b.markers
        side = np.abs(mk.normal[:, 0]) < 0.5
        tilted = (mk.normal[:, 0] <= -0.5) & (mk.normal[:, 0] >= -FACE_NORMAL)
        S, A_group, _, _ = b._interface()                   # the face skins' groups now (Task 11's fix, item 1)
        face_rows = b.markers.skin_index[b.markers.normal[:, 0] < -FACE_NORMAL]
        fr = A_group[face_rows] / b.skins.base_area()[face_rows]
        return {"recession_mm": (c - c0) * 1e3, "centre_marker_off_axis_mm": d * 1e3,
                "face_group_ratio": [float(fr.min()), float(np.median(fr)), float(fr.max())] if len(fr) else None,
                "face_flatness_mm": float(np.ptp(rec)) if len(rec) else None,
                "face_recession_mm": {"within_2mm": float(np.median(rec[dist < 2e-3])) if (dist < 2e-3).any() else None,
                                      "2_to_4mm": float(np.median(rec[(dist >= 2e-3) & (dist < 4e-3)]))
                                      if ((dist >= 2e-3) & (dist < 4e-3)).any() else None,
                                      "beyond_4mm": float(np.median(rec[dist >= 4e-3])) if (dist >= 4e-3).any() else None,
                                      "min": float(rec.min()), "max": float(rec.max())},
                "edge_recession_mm": (edge - edge0) * 1e3,
                "face_markers": int((mk.normal[:, 0] < -FACE_NORMAL).sum()),
                "face_tilt_max_deg": float(np.degrees(np.arccos(np.clip(-mk.normal[mk.normal[:, 0] < -0.5, 0], -1, 1))).max())
                if (mk.normal[:, 0] < -0.5).any() else None,
                "heat_coverage_by_normals": heat_coverage(b, v_hat),
                "heat_coverage_projected": heat_coverage(b, v_hat, projected_heat(b, BOX_Q, v_hat), BOX_Q),
                "top_zone": {"face": int(mk.top[mk.normal[:, 0] < -FACE_NORMAL].sum()), "edges": int(mk.top[tilted].sum()),
                             "sides": int(mk.top[side].sum()), "back": int(mk.top[mk.normal[:, 0] > 0.5].sum())},
                "skin_thickness_max_mm": float(b.skins.thickness().max() * 1e3),
                "edge_skin_thickness_max_mm": float(b.skins.thickness()[mk.skin_index[tilted]].max() * 1e3)
                if tilted.any() else None}

    n = int(round(T_END / DT)) if steps is None else steps
    rows, failure = run_steps(body, advance, n, "box " + preset, extra)
    at = {round(r["t"], 6): r for r in rows}
    compare, verdict = {}, []
    for t in JUDGED_TIMES:
        if t not in at:
            continue
        r = at[t]["recession_mm"]
        rel, rel_l = r / REFERENCE_1D[t] - 1.0, r / REFERENCE_1D_LIQUIDUS[t] - 1.0
        compare[str(t)] = {"recession_mm": r, "ref_1d_mm": REFERENCE_1D[t], "rel_1d": rel,
                           "ref_liquidus_mm": REFERENCE_1D_LIQUIDUS[t], "rel_liquidus": rel_l,
                           "edge_recession_mm": at[t]["edge_recession_mm"], "top_zone": at[t]["top_zone"]}
        verdict.append(abs(rel) <= RECESSION_TOL)
    trapped = rows[-1]["markers_trapped"] if rows else None
    out = {"case": "box", "preset": preset, "T0_K": COLD_T0, "q_Wm2": BOX_Q, "dt_s": DT, "decel_ms2": PUSH_DECEL,
           "size_mm": [BOX_LENGTH * 1e3, 2 * BOX_HALF * 1e3, 2 * BOX_HALF * 1e3],
           "mesh": {"fine_mm": body.fine_size * 1e3, "floor_mm": body.cell_floor * 1e3,
                    "margin_mm": body.fine_margin * 1e3, "at_start": cells0},
           "markers_at_start": markers0, "face_markers_at_start": int(len(top0)),
           "heat_coverage_by_normals_at_start": coverage0, "build_s": round(t_build, 2),
           "books_pass": books_pass(rows), "markers_trapped": trapped, "compare": compare,
           "passes": len(verdict) == len(JUDGED_TIMES) and all(verdict) and books_pass(rows) and trapped == 0
           and failure is None, "failure": failure, "steps": rows}
    if save:
        write("box_" + preset, out)
    return out


# -- the shoulder ------------------------------------------------------------------------------------------------------
SHOULDER_CENTRE = np.array([0.0, 0.0, -SHOULDER_R])


def shoulder_v_hat():
    """The flight direction in the local frame: the stagnation point's outward normal, PUSH_ALPHA upstream (toward -x)
    of the station's +z."""
    a = math.radians(PUSH_ALPHA)
    return np.array([-math.sin(a), 0.0, math.cos(a)])


def shoulder_beta0(stub):
    c, R = SHOULDER_CENTRE, SHOULDER_R

    def beta0(x):
        return np.maximum.reduce([np.linalg.norm(x - c, axis=1) - R, np.abs(x[:, 0]) - SHOULDER_HALF_X,
                                  np.abs(x[:, 1]) - SHOULDER_HALF_Y, -x[:, 2] - SHOULDER_DEPTH]) + stub
    return beta0


def curved(mk):
    """Which markers stand on the sphere (normal within acos(CURVED_FACE) of the radial direction), and their arc
    coordinate along the meridian past the station [m] (positive downstream: past the heat's end)."""
    r = mk.x - SHOULDER_CENTRE
    rhat = r / np.linalg.norm(r, axis=1)[:, None]
    on = np.einsum("ij,ij->i", mk.normal, rhat) > CURVED_FACE
    s = SHOULDER_R * np.arctan2(r[:, 0], r[:, 2])
    return on, s


class PushedLoads:
    """The pushed plate's loads on the segment's curved face (`loads.GasLoads`' interface): dp/ds and shear down the
    meridian (PUSH_DP_DS, PUSH_TAU) on the markers standing on the sphere, none on the cut faces, the deceleration
    PUSH_DECEL; the convective heating is prescribed by the case (q_prescribed), so the heating returned is zero."""

    def evaluate(self, state, markers, T_top, theta, nose_fit):
        on, _ = curved(markers)
        flow = SimpleNamespace(dp_ds=np.where(on, PUSH_DP_DS, 0.0), tau=np.where(on, PUSH_TAU, 0.0),
                               deceleration=PUSH_DECEL)
        return np.zeros(markers.n), flow


def shoulder_body(preset="dev", mat=None):
    mat = material() if mat is None else mat
    st = skin_settings(preset)
    lo = (-SHOULDER_HALF_X - H0, -SHOULDER_HALF_Y - H0, -SHOULDER_DEPTH - H0)
    hi = (SHOULDER_HALF_X + H0, SHOULDER_HALF_Y + H0, H0)
    base = BoxBase.around(lo, hi, H0)
    body = skinbody.LevelSetSkinBody(mat, base, shoulder_beta0(st.stub), COLD_T0, st, loads=PushedLoads(),
                                     remove="spray", v_hat=shoulder_v_hat(), cell_floor=0.7e-3)
    return body


def upstream(c):
    """Whether the points c (k, 3) lie upstream of the station along the meridian (before the heat's end)."""
    r = c - SHOULDER_CENTRE
    return np.arctan2(r[:, 0], r[:, 2]) < 0.0


def shoulder_heat(body):
    """SHOULDER_Q by projection on the station's normal (+z) over phi's patches upstream of the station on the curved
    face (`projected_heat`); none past the heat's end."""
    return projected_heat(body, SHOULDER_Q, (0.0, 0.0, 1.0), where=upstream)


def deposit(body):
    """The deposit past the heat's end: per marker on the sphere downstream of the station, the solid height gained
    (phi's radial rise over the sphere, the skins at the solid's density) times its area and the solid density [kg];
    its mass, centroid past the heat's end [mm], farthest reach [mm], share 5-10 mm past and share past DEPOSIT_FAR; its
    rigid share, the deposit less the non-rigid metal it holds (each marker's skin cells above T_rigid, at most its
    deposit) over the deposit (the ruling on the review's I2); the mass of non-rigid cells in the skins past the heat's
    end, and the largest distance past the heat's end at which a skin holds metal above T_rigid [mm] (the instantaneous
    reach of non-rigid metal, reported, not judged)."""
    mk, sk, mat = body.markers, body.skins, body.material
    on, s = curved(mk)
    rise = np.linalg.norm(mk.x - SHOULDER_CENTRE, axis=1) - SHOULDER_R
    dm = np.where(on & (s > 0.0) & (rise > 0.0), rise * mk.area * mat.rho, 0.0)
    T = sk.T()
    hot_mass = np.where(sk.valid & (T > mat.T_rigid), sk.m, 0.0).sum(axis=1)[mk.skin_index]
    hot_past = on & (s > 0.0) & (hot_mass > 0.0)
    M = float(dm.sum())
    up = on & (s < 0.0)
    soft = float(np.minimum(hot_mass, dm).sum())
    return {"mass_kg": M,
            "rigid_share": float(1.0 - soft / M) if M > 0 else None,
            "share_past_5mm": float(dm[s >= DEPOSIT_FAR].sum() / M) if M > 0 else None,
            "centroid_past_heat_mm": float((s * dm).sum() / M * 1e3) if M > 0 else None,
            "reach_mm": float(s[dm > 0].max() * 1e3) if M > 0 else None,
            "share_5_to_10mm": float(dm[(s >= 5e-3) & (s < 10e-3)].sum() / M) if M > 0 else None,
            "nonrigid_mass_past_kg": float(hot_mass[on & (s > 0.0)].sum()),
            "hot_reach_past_mm": float(s[hot_past].max() * 1e3) if hot_past.any() else 0.0,
            "upstream_recession_mm": float(-np.median(rise[up]) * 1e3) if up.any() else None,
            "liquid_on_cap_mm": float(np.median(sk.depth_above(mat.T_liquidus)[mk.skin_index[up]]) * 1e3)
            if up.any() else None,
            "markers_past": int((on & (s > 0.0)).sum())}


def shoulder(preset="dev", steps=None, save=True):
    """Runoff over a shoulder (module docstring)."""
    mat = material()
    t_build = time.perf_counter()
    body = shoulder_body(preset, mat)
    t_build = time.perf_counter() - t_build
    state = gas_state()
    markers0 = int(body.markers.n)

    def advance(b):
        return b.advance(DT, state=state, q_prescribed=lambda mk: shoulder_heat(b))

    def extra(b, book):
        d = deposit(b)
        return {"deposit": d, "heat_coverage_projected": heat_coverage(b, (0.0, 0.0, 1.0), shoulder_heat(b), SHOULDER_Q,
                                                                      upstream), "momentum_to_wall": book.momentum_to_wall.tolist(),
                "momentum_to_volume": book.momentum_to_volume.tolist()}

    n = int(round(T_END / DT)) if steps is None else steps
    rows, failure = run_steps(body, advance, n, "shoulder " + preset, extra)
    last = rows[-1]["deposit"] if rows and abs(rows[-1]["t"] - T_END) < 1e-9 else {}
    formed = bool(last) and last["mass_kg"] > 0.0
    rigid = formed and last["rigid_share"] >= DEPOSIT_RIGID_SHARE
    near = formed and last["share_past_5mm"] <= SKYSCRAPERS_DEPOSIT_FAR_SHARE
    out = {"case": "shoulder", "preset": preset, "T0_K": COLD_T0, "q_Wm2": SHOULDER_Q, "dt_s": DT,
           "sphere_radius_mm": SHOULDER_R * 1e3,
           "segment_mm": [2 * SHOULDER_HALF_X * 1e3, 2 * SHOULDER_HALF_Y * 1e3, SHOULDER_DEPTH * 1e3],
           "push": {"alpha_deg": PUSH_ALPHA, "dp_ds_Pa_m": PUSH_DP_DS, "tau_Pa": PUSH_TAU, "decel_ms2": PUSH_DECEL},
           "markers_at_start": markers0, "build_s": round(t_build, 2),
           "books_pass": books_pass(rows), "deposit_formed": formed, "deposit_rigid": rigid, "deposit_near": near,
           "deposit_rigid_share_min": DEPOSIT_RIGID_SHARE, "deposit_far_mm": DEPOSIT_FAR * 1e3,
           "skyscrapers_pushed_plate": {"centroid_past_heat_mm": SKYSCRAPERS_DEPOSIT_CENTROID_MM,
                                        "share_5_to_10mm": SKYSCRAPERS_DEPOSIT_FAR_SHARE},
           "passes": books_pass(rows) and formed and rigid and near and failure is None, "failure": failure,
           "steps": rows}
    if save:
        write("shoulder_" + preset, out)
    return out


# -- the cost ----------------------------------------------------------------------------------------------------------
def sphere_body(preset, mat=None):
    """The 50 mm sphere at the melting state (module docstring), with the shared loads and Girin's spray model."""
    mat = material() if mat is None else mat
    gl = loads_mod.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), mat.liquid.rho)
    model = spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), spray_gates="2026-11-10")
    st = SkinSettings(cell=SKIN_PRESETS[preset])
    L = SPHERE_R + 2 * H0
    base = BoxBase.around((-L,) * 3, (L,) * 3, H0)
    return skinbody.LevelSetSkinBody(mat, base, lambda x: np.linalg.norm(x, axis=1) - (SPHERE_R - st.stub), SPHERE_T0,
                                     st, loads=gl, spray_model=model, emissivity=0.2, T_ambient=200.0,
                                     v_hat=(1.0, 0.0, 0.0))


def cost_sphere(preset="dev", steps=2, dt=DT, save=True):
    mat = material()
    t0 = time.perf_counter()
    body = sphere_body(preset, mat)
    build_s = time.perf_counter() - t0
    built = {"build_s": round(build_s, 1), "peak_rss_gb": round(peak_rss_gb(), 3), "markers": int(body.markers.n),
             "cells": int(len(body.mesh.tets)), "active_nodes": int(body.solver.active.sum()),
             "by_size_mm": {str(round(s * 1e3, 4)): int(c) for s, c in zip(*np.unique(body.mesh.size, return_counts=True))}}
    print("sphere {} built in {:.1f} s: {}".format(preset, build_s, built), flush=True)
    state = gas_state(**CONTINUUM)
    rows, failure = run_steps(body, lambda b: b.advance(dt, state=state), steps, "sphere " + preset,
                              lambda b, book: {"skins_cells": int(b.skins.count.sum()),
                                               "skins_cells_max": int(b.skins.count.max())})
    out = {"case": "cost-sphere", "preset": preset, "diameter_mm": 2 * SPHERE_R * 1e3, "T0_K": SPHERE_T0, "dt_s": dt,
           "gas_state": CONTINUUM, "built": built, "failure": failure, "steps": rows}
    if save:
        probe = steps != 2 or dt != DT                 # a shorter probe is kept apart from the measurement
        write("cost_sphere_{}{}".format(preset, "_probe_{}x{:g}s".format(steps, dt) if probe else ""), out)
    return out


def puck_sdf(x):
    """A puck of radius RIM_R and half height RIM_HALF_Z whose rims are rounded to RIM_EDGE (a signed distance)."""
    rho = np.hypot(x[:, 0], x[:, 1])
    q = np.column_stack([rho - (RIM_R - RIM_EDGE), np.abs(x[:, 2]) - (RIM_HALF_Z - RIM_EDGE)])
    return np.linalg.norm(np.maximum(q, 0.0), axis=1) + np.minimum(q.max(axis=1), 0.0) - RIM_EDGE


def cost_rim(save=True):
    """The cost of cells at the 0.175 mm floor on a rim of radius 0.5 mm (design 20.13): the puck meshed as the body
    meshes (the targets of `criteria.make_target` from phi and its curvature, rebuilt three times from its own fields, the
    design's 0.7 mm fine cells over a 4 mm margin) with the design's 0.25 mm floor and with none below the fine size;
    their cells, build times, the cut-cell solve's setup and one 0.5 s step (2 MW/m2 by projection along +x), and the
    process's peak memory after each."""
    mat = material()
    L = RIM_R + 2 * H0
    base = BoxBase.around((-L, -L, -RIM_HALF_Z - 2 * H0), (L, L, RIM_HALF_Z + 2 * H0), H0)
    out = {"case": "cost-rim", "puck_mm": {"radius": RIM_R * 1e3, "height": 2 * RIM_HALF_Z * 1e3,
                                           "rim_radius": RIM_EDGE * 1e3},
           "rim_length_mm": 2 * 2 * math.pi * (RIM_R - RIM_EDGE) * 1e3}
    for name, floor in (("no_floor", 0.7e-3), ("floor_0p25", 0.25e-3)):
        t0 = time.perf_counter()
        mesh = build(base)
        phi = redistance(mesh, puck_sdf(mesh.points), band=2 * (4e-3 + H0))
        for _ in range(3):
            km = kappa_max_nodal(mesh, phi, band=2 * H0)
            mesh = build(base, make_target(mesh, phi, phi, km, km, fine=0.7e-3, margin=4e-3, floor=floor))
            phi = redistance(mesh, puck_sdf(mesh.points), band=2 * (4e-3 + H0))
        t_build = time.perf_counter() - t0
        solver = CutThermal(mat, emissivity=0.0, linear_solver="amg")
        t0 = time.perf_counter()
        solver.setup(mesh, phi)
        solver.set_temperature(700.0)
        t_setup = time.perf_counter() - t0
        q = 2e6 * np.maximum(solver.P.normal[:, 0], 0.0)
        t0 = time.perf_counter()
        solver.step(0.5, q)
        t_step = time.perf_counter() - t0
        sizes, counts = np.unique(mesh.size, return_counts=True)
        out[name] = {"floor_mm": floor * 1e3, "cells": int(len(mesh.tets)), "active_nodes": int(solver.active.sum()),
                     "by_size_mm": {str(round(s * 1e3, 4)): int(c) for s, c in zip(sizes, counts)},
                     "build_s": round(t_build, 2), "setup_s": round(t_setup, 2), "step_s": round(t_step, 2),
                     "newton_iterations": int(solver.iterations), "peak_rss_gb": round(peak_rss_gb(), 3)}
        print("rim", name, out[name], flush=True)
        del solver, mesh
    a, b = out["no_floor"], out["floor_0p25"]
    floor_cells = b["by_size_mm"].get("0.175", 0)
    out["extra"] = {"cells": b["cells"] - a["cells"], "cells_at_0p175": floor_cells,
                    "cells_at_0p175_per_mm_of_rim": floor_cells / out["rim_length_mm"],
                    "extra_cells_per_mm_of_rim": (b["cells"] - a["cells"]) / out["rim_length_mm"],
                    "step_ratio": b["step_s"] / max(a["step_s"], 1e-9), "build_ratio": b["build_s"] / max(a["build_s"], 1e-9)}
    if save:
        write("cost_rim", out)
    return out


def cost_skip(preset="dev", save=True):
    """psi and the remesh every second step (module docstring): the cold plate both ways."""
    every = plate("cold", preset, save=False)
    second = plate("cold", preset, every_second=True, save=False)
    out = {"case": "cost-skip", "preset": preset, "compare": {}}
    for a, b in zip(every["steps"], second["steps"]):
        out["compare"][str(a["t"])] = {"recession_every_mm": a["recession_mm"], "recession_second_mm": b["recession_mm"],
                                       "change_mm": b["recession_mm"] - a["recession_mm"],
                                       "change_rel": b["recession_mm"] / a["recession_mm"] - 1.0
                                       if a["recession_mm"] else None,
                                       "wall_every_s": a["wall_s"], "wall_second_s": b["wall_s"],
                                       "psi_every_s": a["cost_split_s"]["surfaces psi"],
                                       "psi_second_s": b["cost_split_s"]["surfaces psi"],
                                       "remesh_every_s": a["cost_split_s"]["remesh"],
                                       "remesh_second_s": b["cost_split_s"]["remesh"],
                                       "cells_every": a["cells"], "cells_second": b["cells"]}
    out["books_pass"] = every["books_pass"] and second["books_pass"]
    out["second_passes_the_plate"] = second["passes"]
    out["shadow"] = "not exercised: the plate is heated by prescription (no gas state); timed on the sphere"
    if save:
        write("cost_skip_" + preset, out)
    return out


# -- the record --------------------------------------------------------------------------------------------------------
def record():
    """Gather output/levelset/stage3_*.json into analysis/levelset/stage3_recorded.json (each case's summary and steps)."""
    out = {}
    for path in sorted(glob.glob(os.path.join(OUT, "stage3_*.json"))):
        with open(path) as fh:
            out[os.path.basename(path)[len("stage3_"):-len(".json")]] = json.load(fh)
    with open(RECORD, "w") as fh:
        json.dump(out, fh, indent=1)
    print("recorded", sorted(out), "in", RECORD)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("case", choices=("plate-cold", "plate-hot", "box", "shoulder", "cost", "record"))
    p.add_argument("--preset", choices=tuple(SKIN_PRESETS), default="dev")
    p.add_argument("--every-second", action="store_true", help="plates: psi and the remesh at every second step end")
    p.add_argument("--what", choices=("sphere", "rim", "skip"), default="sphere", help="cost: which measurement")
    p.add_argument("--steps", type=int, default=None, help="macro steps (default: to 2 s; the sphere 2)")
    p.add_argument("--dt", type=float, default=DT, help="cost sphere: the macro step [s]")
    a = p.parse_args(argv)
    if a.case in ("plate-cold", "plate-hot"):
        out = plate(a.case.split("-")[1], a.preset, a.every_second, a.steps)
    elif a.case == "box":
        out = box(a.preset, a.steps)
    elif a.case == "shoulder":
        out = shoulder(a.preset, a.steps)
    elif a.case == "cost":
        if a.what == "sphere":
            out = cost_sphere(a.preset, 2 if a.steps is None else a.steps, a.dt)
        elif a.what == "rim":
            out = cost_rim()
        else:
            out = cost_skip(a.preset)
    else:
        record()
        return 0
    if "passes" in out:
        print("passes:", out["passes"], json.dumps(out.get("compare", {}), default=_json_default))
    return 0


if __name__ == "__main__":
    sys.exit(main())
