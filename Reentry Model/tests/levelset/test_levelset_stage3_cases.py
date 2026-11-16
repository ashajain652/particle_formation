"""Level-set stage 3's verification cases (Task 11; design section 20.11 item 3, plan Task 11): the plates at 0.5 s macro
steps (Decision A of stage 2), the q6p box (Ruling 5) and runoff over a shoulder, each run by
`analysis/levelset/stage3_cases.py` (whose docstring describes the cases and why the plate's face carries one marker) and
judged here against the brief's tolerances. The cases are marked `reference` (minutes each); the import check runs in
the unit tier. The measured numbers are recorded in `analysis/levelset/stage3_recorded.json`."""
import ast
import importlib.util
import os
import pickle
from types import SimpleNamespace

import numpy as np
import pytest

from helpers import PROJECT_ROOT

SCRIPT = os.path.join(PROJECT_ROOT, "analysis", "levelset", "stage3_cases.py")
ROOT_CONSTANTS = ("OUTPUT_DIR", "DATA_DIR", "PROJECT_ROOT")


def _cases():
    spec = importlib.util.spec_from_file_location("stage3_cases", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _escapes(node):
    """True when an import statement reaches a reentry_model module other than the shared package, the level-set package
    or a package-root constant (the skyscrapers' modules and the old top-level aliases are all outside)."""
    ok = ("reentry_model.shared", "reentry_model.levelset")
    if isinstance(node, ast.Import):
        return any(a.name.split(".")[0] == "reentry_model" and not a.name.startswith(ok) for a in node.names)
    mod = node.module or ""
    if node.level or mod.split(".")[0] != "reentry_model":
        return False
    if mod == "reentry_model":
        return not all(a.name in ROOT_CONSTANTS for a in node.names)
    return not (mod == ok[0] or mod == ok[1] or mod.startswith(ok[0] + ".") or mod.startswith(ok[1] + "."))


def test_the_cases_use_only_the_shared_physics_and_the_level_set():
    """The cases' script and these tests import nothing of the skyscrapers (nor an old alias module): only the shared
    package, the level set and package-root constants."""
    for path in (SCRIPT, os.path.abspath(__file__)):
        with open(path) as fh:
            tree = ast.parse(fh.read())
        bad = [node.lineno for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom)) and _escapes(node)]
        assert bad == [], (path, bad)
    for src in ("from reentry_model import skin", "from reentry_model.skyscrapers import skinbody", "import reentry_model.aero",
                "from reentry_model.material import Material"):
        assert _escapes(ast.parse(src).body[0]), src
    for src in ("from reentry_model.shared.material import Material", "from reentry_model.levelset import division",
                "from reentry_model import OUTPUT_DIR"):
        assert not _escapes(ast.parse(src).body[0]), src


class _Stalled:
    """A stand-in body for the case script's step loop whose step raises (Task 11's review, M6)."""

    def __init__(self):
        self.time = 1.5
        self.markers = SimpleNamespace(x=np.zeros((2, 3)))

    def _rebuild(self, tally, surface=None, check_core=False):
        return None


def _stall(body):
    raise RuntimeError("volume correction of piece 1 did not converge")


def test_a_failing_step_keeps_the_body_it_started_from(monkeypatch, tmp_path):
    """Task 11's review, M6: with STAGE3_STATE_DIR set, the case script's step loop pickles the body before each step and,
    when the step raises, keeps that state and names it in the failure, so that the failure can be diagnosed from the
    body it started from (the q6p box's volume correction failure could not be)."""
    c = _cases()
    monkeypatch.setenv("STAGE3_STATE_DIR", str(tmp_path))
    monkeypatch.setattr(c, "OUT", str(tmp_path))
    rows, failure = c.run_steps(_Stalled(), _stall, 2, "stalled case")
    assert rows == [] and failure["step"] == 1 and "did not converge" in failure["error"]
    with open(failure["state"], "rb") as fh:
        kept = pickle.load(fh)
    assert kept.time == 1.5 and "_rebuild" not in kept.__dict__


def test_a_stopped_step_s_failure_record_carries_the_case_s_measures(monkeypatch, tmp_path):
    """Plan 03b Task 7b (Task 7's carried minor M2): a step that stops leaves in its failure record the case's own
    measures of the body as the stop left it (the box's face grid and edge radii), where the state allows, so that a
    strict expected failure's reason cites recorded numbers only; a measure the state does not allow is recorded as the
    error it raised."""
    c = _cases()
    monkeypatch.setattr(c, "OUT", str(tmp_path))
    monkeypatch.delenv("STAGE3_STATE_DIR", raising=False)
    _, failure = c.run_steps(_Stalled(), _stall, 1, "stalled case", failed=lambda b: {"time_s": b.time})
    assert failure["measures"] == {"time_s": 1.5}

    def unmeasurable(b):
        raise ValueError("no face left to read")
    _, failure = c.run_steps(_Stalled(), _stall, 1, "stalled case", failed=unmeasurable)
    assert failure["measures"] == {"error": "ValueError: no face left to read"}


def test_the_projected_heat_is_meaned_under_the_ownership_that_sets_the_markers_areas(monkeypatch):
    """The controller's ruling T7-i (diagnosis D1): the cases' beam (`projected_heat`) gives each marker the mean
    projection of the patches of phi it owns under `Markers.owned_patches`, the ownership (with the same-side and
    same-piece rules) that sets its area, not under `loads.owned` (the nearest marker alone), so that the heat the
    markers take, sum q_i A_i, is the beam's power on phi. The hot plate's column, its face one marker: the projection
    reads phi's patches through `owned_patches` on phi as it is and never through `loads.owned`, and the heat taken is
    the beam's to 1e-12."""
    from reentry_model.levelset import loads, markers
    from reentry_model.levelset.cutcell import patches
    c = _cases()
    body = c.plate_body("hot", "dev")
    used = []
    owned_patches = markers.Markers.owned_patches

    def spy(self, mesh, phi, P=None, live=None):
        used.append(phi is body.phi)
        return owned_patches(self, mesh, phi, P, live)
    monkeypatch.setattr(markers.Markers, "owned_patches", spy)
    monkeypatch.setattr(loads, "owned", lambda *a, **k: pytest.fail("projected_heat used loads.owned"))
    q0 = c.PLATE_CASES["hot"][1]
    q = c.projected_heat(body, q0, (0.0, 0.0, 1.0))
    P = patches(body.mesh, body.phi)
    beam = q0 * float((P.area * np.maximum(P.normal[:, 2], 0.0)).sum())
    assert used == [True]
    assert float((q * body.markers.area).sum()) == pytest.approx(beam, rel=1e-12)


SHOULDER_PROBE_DEPTH = 3e-3   # m: below the station, where the shoulder's faces are probed (inside the sphere there)
SHOULDER_FINEST = 8           # the shoulder's planes of nodes checked down to its base cells over this (0.35 mm, half its
                              # finest cells: fine size and curvature floor 0.7 mm)


def test_the_shoulder_s_flat_faces_lie_off_its_meshes_planes_of_nodes():
    """The controller's ruling T7-h (diagnosis D2): the shoulder is built SHOULDER_OFF (0.13 mm) off its mesh's planes of
    nodes in each axis, as the box (BOX_OFF) and the test cube are: with its cut faces on those planes, 451 nodes had phi
    = 0 exactly and its division closed 3.7e-3 from one at 24 levels (first order in the level spacing), against 3.7e-5
    shifted. Each of its five flat faces (x = +-SHOULDER_HALF_X, y = +-SHOULDER_HALF_Y and z = -SHOULDER_DEPTH about the
    station, which stands SHOULDER_OFF from the origin along each axis) bounds the metal where it should (phi0 changes
    sign across it, probed SHOULDER_PROBE_DEPTH below the station) and lies at least SHOULDER_OFF from every plane of
    nodes of the case's mesh box, down to cells of its base cells over SHOULDER_FINEST; the sphere's centre lies
    SHOULDER_R below the station."""
    c = _cases()
    station = np.full(3, c.SHOULDER_OFF)
    lo = np.array([-c.SHOULDER_HALF_X, -c.SHOULDER_HALF_Y, -c.SHOULDER_DEPTH]) - c.H0      # the mesh box's corner
    phi0 = c.shoulder_phi0()
    eps = 1e-7
    faces = [(0, s, station[0] + s * c.SHOULDER_HALF_X) for s in (-1, 1)] + \
            [(1, s, station[1] + s * c.SHOULDER_HALF_Y) for s in (-1, 1)] + [(2, -1, station[2] - c.SHOULDER_DEPTH)]
    for axis, s, at in faces:
        for h in c.H0 / 2.0 ** np.arange(int(np.log2(SHOULDER_FINEST)) + 1):
            k = (at - lo[axis]) / h
            assert abs(k - round(k)) * h >= c.SHOULDER_OFF * (1.0 - 1e-9), (axis, at, h)
        p = station - [0.0, 0.0, SHOULDER_PROBE_DEPTH]
        p[axis] = at
        inside, outside = p.copy(), p.copy()
        inside[axis] -= s * eps
        outside[axis] += s * eps
        assert phi0(inside[None])[0] < 0.0 < phi0(outside[None])[0], (axis, s)
    np.testing.assert_allclose(c.SHOULDER_CENTRE, station - [0.0, 0.0, c.SHOULDER_R], rtol=0.0, atol=1e-15)


def test_a_case_s_summary_is_stamped_with_the_code_it_ran_on(monkeypatch, tmp_path):
    """Task 11's review, M5: each case's summary carries the commit it ran on and whether the level-set package or the
    case script differed from it, so that the record cannot mix runs from different code unseen; the stamp is taken when
    the script starts, not when the summary is written (a box run outlives commits made meanwhile); `record` leaves out
    a summary without a stamp."""
    import json
    import re
    import subprocess
    c = _cases()
    monkeypatch.setattr(c, "OUT", str(tmp_path))
    monkeypatch.setattr(c, "RECORD", str(tmp_path / "recorded.json"))
    monkeypatch.setattr(c, "code_stamp", lambda: {"commit": "0" * 40, "dirty": True})    # a commit made since the start
    path = c.write("stamped", {"case": "a"})
    with open(path) as fh:
        doc = json.load(fh)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=c.REPO_ROOT, capture_output=True, text=True).stdout.strip()
    assert re.fullmatch("[0-9a-f]{40}", doc["code"]["commit"]) and doc["code"]["commit"] == head
    assert isinstance(doc["code"]["dirty"], bool)
    with open(tmp_path / "stage3_old.json", "w") as fh:
        json.dump({"case": "b"}, fh)
    assert sorted(c.record()) == ["stamped"]


def test_the_code_stamp_s_difference_check_covers_the_data_the_cases_read(monkeypatch):
    """Task 11's re-review, minor item 4: a case's stamp says whether the code it ran on differed from its commit; the
    check covered only `reentry_model` and `analysis/levelset`, so a changed material table or reference value (the cases
    read `data/materials`, `data/reference_values`, and through the shared package the US76 table and the drag tables in
    `data/`) would not have marked the run. The stamp's `git status` now covers `data` too, and a change there alone
    marks the run as differing."""
    import subprocess
    c = _cases()
    calls = []

    def fake(args, **kw):
        calls.append(list(args))
        if args[:2] == ["git", "rev-parse"]:
            return subprocess.CompletedProcess(args, 0, stdout="0" * 40 + "\n", stderr="")
        changed = " M data/materials/AA7075_scheil.json" if "data" in args else ""
        return subprocess.CompletedProcess(args, 0, stdout=changed, stderr="")
    monkeypatch.setattr(c.subprocess, "run", fake)
    stamp = c.code_stamp()
    status = [a for a in calls if a[:2] == ["git", "status"]][0]
    assert {"reentry_model", "analysis/levelset", "data"} <= set(status[status.index("--") + 1:])
    assert stamp == {"commit": "0" * 40, "dirty": True}


def test_the_shoulder_s_deposit_is_judged_by_its_rigid_share_and_its_far_share():
    """The controller's ruling on the review's I2: the shoulder is judged on the metal deposited past the heat's end at
    2.0 s, at least DEPOSIT_RIGID_SHARE of it rigid and at most the skyscrapers' 3.0 % of it past DEPOSIT_FAR (5 mm);
    the instantaneous reach of non-rigid metal is reported, not judged. Three markers past the heat's end on the
    sphere, 1, 3 and 6 mm along the meridian, each with 0.1 mm of deposit on 4 mm2; the first holds 2 mg above T_rigid
    (more than its deposit's 0.97 mg, so its whole deposit counts as non-rigid), the second 0.5 mg: the rigid share is
    1 - (0.97 + 0.5) / 2.9 of the deposit, a third of it lies past 5 mm, and the reach of non-rigid metal is 3 mm."""
    c = _cases()
    rho, R = 2420.0, c.SHOULDER_R
    s = np.array([1e-3, 3e-3, 6e-3])
    ang = s / R
    rhat = np.column_stack([np.sin(ang), np.zeros(3), np.cos(ang)])
    rise = 0.1e-3
    x = c.SHOULDER_CENTRE + (R + rise) * rhat
    area = np.full(3, 4e-6)
    hot = np.array([2e-6, 0.5e-6, 0.0])
    m = np.column_stack([hot, np.full(3, 1e-5)])
    T = np.column_stack([np.full(3, 950.0), np.full(3, 700.0)])
    mk = SimpleNamespace(x=x, normal=rhat, area=area, skin_index=np.arange(3))
    sk = SimpleNamespace(T=lambda: T, valid=np.ones((3, 2), bool), m=m, depth_above=lambda t: np.zeros(3))
    mat = SimpleNamespace(rho=rho, T_rigid=895.0, T_liquidus=908.0)
    d = c.deposit(SimpleNamespace(markers=mk, skins=sk, material=mat))
    dm = rise * area * rho
    assert d["mass_kg"] == pytest.approx(dm.sum(), rel=1e-12)
    assert d["rigid_share"] == pytest.approx(1.0 - (dm[0] + 0.5e-6) / dm.sum(), rel=1e-12)
    assert d["share_past_5mm"] == pytest.approx(1.0 / 3.0, rel=1e-12)
    assert d["hot_reach_past_mm"] == pytest.approx(3.0, rel=1e-9)


def test_the_solver_watch_keeps_each_step_s_largest_interface_mismatch_on_the_tolerance_s_measure(monkeypatch,
                                                                                                     tmp_path):
    """Plan 03b Task 7's fix round 1 (its review's I1, the controller's carry 5): the solver watch's
    `interface_mismatch_max` is, over a step's accepted parts, the largest of the measure the interface's repetition
    ends on against INTERFACE_TOLERANCE, sum |mismatch| over sum (|E_top| + |E_base|) of the part's real pass
    (`skinbody.LevelSetSkinBody._part_once`), not the last part's mismatch over the step's heat in; and
    `passes_accepted_max` the most passes an accepted part took. The hot plate's first two 0.5 s steps (one part, then
    several), each part's measure and passes taken here from what `_accept` receives."""
    from reentry_model.levelset import skinbody
    c = _cases()
    monkeypatch.setattr(c, "OUT", str(tmp_path))
    monkeypatch.delenv("STAGE3_STATE_DIR", raising=False)
    parts = {}
    accept = skinbody.LevelSetSkinBody._accept

    def measured(body, dt, dt_macro, ratio, passes, book, mismatch, *a, **k):
        scale = float((np.abs(book.E_top) + np.abs(book.E_base)).sum())
        parts.setdefault(round(body.time, 9), []).append((float(np.abs(mismatch).sum()) / scale, passes))
        return accept(body, dt, dt_macro, ratio, passes, book, mismatch, *a, **k)
    monkeypatch.setattr(skinbody.LevelSetSkinBody, "_accept", measured)
    _, q = c.PLATE_CASES["hot"]
    body = c.plate_body("hot", "dev")
    rows, failure = c.run_steps(body, lambda b: b.advance(c.DT, q_prescribed=lambda mk: c.projected_heat(
        b, q, (0.0, 0.0, 1.0))), 2, "mismatch watch")
    assert failure is None and len(rows) == len(parts) == 2 and len(parts[c.DT]) > 1
    for row, step in zip(rows, parts.values()):
        assert row["solver_watch"]["interface_mismatch_max"] == pytest.approx(max(m for m, _ in step), rel=1e-12)
        assert row["solver_watch"]["interface_mismatch_max"] > 0.0
        assert row["solver_watch"]["passes_accepted_max"] == max(n for _, n in step) == row["passes_max"]


def _books(out):
    """Every step's books (mass to MASS_TOL, enthalpy to ENERGY_TOL) and its division's closure (planning ruling 1: a
    piece's closure factor more than DIVISION_TOL from one fails the step's check, raised by the tests; plan 03b Task 7's
    fix round 1, its review's M5)."""
    from reentry_model.levelset.division import DIVISION_TOL
    c = _cases()
    for row in out["steps"]:
        assert abs(row["mass_residual"]) <= c.MASS_TOL, row["t"]
        assert abs(row["enthalpy_residual"]) <= c.ENERGY_TOL, row["t"]
        assert abs(row["layer_factor_max"]) <= DIVISION_TOL, row["t"]


@pytest.mark.reference
@pytest.mark.parametrize("preset", ["dev", "production"])
def test_cold_plate_at_half_second_steps(preset):
    """The level-2 plate's centre column (8 mm deep, 700 K, 6 MW/m2, its walls the mesh box's insulated faces, the
    shared AA7075 Scheil alloy, no radiation, swell off, flow off, the liquid leaving at T_feed within the skins'
    sub-steps) at 0.5 s macro steps: the face's recession within 5 % of the one-dimensional references with the removal
    at T_feed (0.789, 2.712, 5.219 mm) and at the liquidus (0.837, 2.772, 5.298 mm) at 1.0, 1.5 and 2.0 s, and every
    step's mass to 1e-12 and enthalpy to 1e-8. Measured 2026-11-15 on the one skin depth with phi evolved (plan 03b
    Task 7, 721efce; the same to every digit on Task 7b's bfb6344): 0.800, 2.719 and 5.202 mm (dev; production 0.810,
    2.724, 5.205), 1.4 and 0.25 % ahead of the first and 0.33 % behind it (production 2.7, 0.44 % ahead, 0.27 %
    behind), and 4.4, 1.9 and 1.8 % behind the second
    (production 3.2, 1.7, 1.8 %); books at most 5.4e-15 (mass) and 1.5e-14 (enthalpy). (2026-11-12, the old design:
    0.820, 2.762, 5.258 mm.) Re-measured 2026-11-15 on plan 03b Task 11's code (e57945f), under the composite rheology
    and main's merged physics: every number as on Task 7d's code (73b5a9a) to every digit, the books at most 1.7e-15
    (mass) and 1.4e-14 (enthalpy). The plate's face 0.27 mm higher on the same mesh (8.27 mm of metal; the controller's
    ruling T7-m) is recorded, not tested here: plan 03b Task 11's report."""
    c = _cases()
    out = c.plate("cold", preset)
    assert out["failure"] is None and out["markers"] == 1
    _books(out)
    for t in c.JUDGED_TIMES:
        r = out["compare"][str(t)]
        assert abs(r["rel_1d"]) <= c.RECESSION_TOL, (t, r)
        assert abs(r["rel_liquidus"]) <= c.RECESSION_TOL, (t, r)
    assert out["passes"]


@pytest.mark.reference
@pytest.mark.parametrize("preset", ["dev", "production"])
def test_hot_plate_at_half_second_steps(preset):
    """The hot plate (880 K, 2 MW/m2, otherwise the cold plate's column) at 0.5 s macro steps: the recession within 5 %
    of data/reference_values/plate_hot_1d.json at 1.5 and 2.0 s (1.0 s reported: the melting onset), the non-rigid depth
    (phi to psi on the axis) within 10 % at 1.0, 1.5 and 2.0 s, and the books (plan 03b Task 7's judgment). Its
    non-rigid depth passes the skin depth H = 1.2 mm from 1.0 s, so its melting front crosses the skin base (the risk of
    design 21.12); the volume's top beneath the skin is then above T_rigid by construction, so the skyscrapers' criterion
    on it (HOT_BASE_TOL_K, judged here from Task 11's fix, item 4, under design 20.6's full-depth skins, which design
    21.13 retires) is reported, not judged. Measured 2026-11-15 on 721efce, the same to every digit on Task 7b's
    bfb6344 (dev; production within 0.7 % of it): recession 0.760, 1.795, 2.894 mm against 0.789, 1.817, 2.917 (3.7 %
    behind, reported; 1.2 and 0.8 % behind); non-rigid depth 1.489, 1.583, 1.683 mm against 1.487, 1.570, 1.674 (0.1,
    0.8 and 0.5 % thicker); the base 1.75, 2.13 and 2.59 K above T_rigid (reported); books at most 1.6e-15 (mass) and
    1.3e-14 (enthalpy). Re-measured 2026-11-15 on plan 03b Task 11's code (e57945f), under the composite rheology and
    main's merged physics: every number as on Task 7d's code (73b5a9a) to every digit."""
    c = _cases()
    out = c.plate("hot", preset)
    assert out["failure"] is None and out["markers"] == 1
    _books(out)
    for t in c.JUDGED_TIMES:
        r = out["compare"][str(t)]
        if t > 1.25:
            assert abs(r["rel"]) <= c.RECESSION_TOL, (t, r)
        assert abs(r["nonrigid_rel"]) <= c.NONRIGID_TOL, (t, r)
    assert out["passes"]


BOX_FAILS = ("the 2 mm markers do not resolve the sharp heated edge (diagnosis D1, the controller's ruling T7-j; design "
             "21.10's finer markers near edges are the remedy): markers on the side faces own the heated face's outer "
             "strip, which does not recede. Measured 2026-11-15 on plan 03b Task 11's code (e57945f), under the "
             "composite rheology, to 1.0 s (ruling T7d-d): the first step (0.5 s, nothing yet melting) keeps the face "
             "flat and the books exact; the second closes its books (mass 2.9e-15 of the stored mass, against MASS_TOL "
             "1e-12), but the face centre has receded 0.884 mm at 1.0 s against the reference's 0.789 (12.1 % ahead, "
             "past RECESSION_TOL; 0.889 mm and 12.7 % under the bridge law, which this code reproduces to every digit); "
             "the third step, from 1.0 s, is not run on a laptop: on Task 7d's code (73b5a9a) it had not ended after "
             "115 minutes, when its one run was stopped at the two-hour limit, so 1.5 and 2.0 s are not reached")


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=BOX_FAILS)
def test_q6p_box():
    """The q6p box (8 x 10 x 10 mm, 700 K, 6 MW/m2 by projection on its 10 x 10 face, 2 s at 0.5 s steps, dev cells, the
    design's meshes, the volume flow held and the skins' runoff on under the flight's deceleration, the liquid leaving at
    T_feed): the face centre's recession within 5 % of the cold plate's one-dimensional reference (T_feed) at 1.0, 1.5
    and 2.0 s, the books every step, and no marker trapped; the edges are recorded (each edge's rounding radius from
    phi's zero contour at its midpoint, and the recession along the heated face's +y edge's bisector). It fails
    (BOX_FAILS): the brief's tolerance is kept and the failure marked strict."""
    c = _cases()
    out = c.box("dev")
    assert out["failure"] is None
    _books(out)
    assert out["markers_trapped"] == 0
    for t in c.JUDGED_TIMES:
        r = out["compare"][str(t)]
        assert abs(r["rel_1d"]) <= c.RECESSION_TOL, (t, r)
    assert out["passes"]


SHOULDER_FAILS = ("measured 2026-11-15 on plan 03b Task 11's code (e57945f), under the composite rheology, the "
                  "segment built 0.13 mm off the mesh's planes of nodes (ruling T7-h): the first step (to 0.5 s) is "
                  "Task 7d's to every digit; the second accepts its 8 parts and then stops at its end on the "
                  "thin-metal check (stepend.THIN_METAL, design 21.11): two pairs of columns whose normals turn 146 "
                  "and 142 degrees share a boundary in the layer at the heated cap's +y edge, between a marker on the "
                  "cap, 0.37 mm inside the starting surface and facing back across the cap, and two on the +y cut "
                  "face; the metal through the cap's marker along its normal is 1.39 mm thick and the thin test's "
                  "measures at the pairs' edges are 0.84 to 2.06 mm, all under twice the skin depth (2.4 mm), so the "
                  "site is really thin. Under the bridge law the same code reproduces Task 7d's steps to 1.0 s to "
                  "every digit, so the runnier slurry of the composite law makes the thin site; 2.0 s, where the "
                  "deposit is judged, is not reached")


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=SHOULDER_FAILS)
def test_shoulder_runoff_freezes_a_deposit():
    """Runoff over a shoulder: a segment of the 50 mm sphere about its 45 degree station at the skyscrapers' pushed
    plate's loads (dp/ds -3.424e5 Pa/m and 80.29 Pa of shear down the meridian, 29.77 m/s2 of deceleration along the
    flight direction), 700 K, 6 MW/m2 by projection on the windward cap upstream of the station and nothing behind, the
    liquid kept, dev cells: the books close every step, a deposit forms past the heat's end, at least
    DEPOSIT_RIGID_SHARE of it is rigid at 2.0 s and at most the pushed plate's 3.0 % of it lies past DEPOSIT_FAR (the
    controller's ruling on the review's I2; the farthest reach of non-rigid metal is reported, not judged). The
    deposit's centroid past the heat's end is recorded against the skyscrapers' pushed plate's 1.05 mm. It fails
    (SHOULDER_FAILS): measured 2026-11-15 on plan 03b Task 11's code (e57945f), under the composite rheology, it stops at
    the end of its second step, on metal thinner than twice the skin depth at the heated cap's edge (on Task 7d's code,
    73b5a9a, under the bridge law, it stopped in its third)."""
    c = _cases()
    out = c.shoulder("dev")
    assert out["failure"] is None
    _books(out)
    last = out["steps"][-1]["deposit"]
    assert out["deposit_formed"] and last["mass_kg"] > 0.0
    assert last["rigid_share"] >= c.DEPOSIT_RIGID_SHARE                # the controller's ruling on the review's I2
    assert last["share_past_5mm"] <= c.SKYSCRAPERS_DEPOSIT_FAR_SHARE
    assert out["passes"]
