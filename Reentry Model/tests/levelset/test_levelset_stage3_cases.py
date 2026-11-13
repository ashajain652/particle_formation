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
    for src in ("from reentry_model.shared.material import Material", "from reentry_model.levelset import modes",
                "from reentry_model import OUTPUT_DIR"):
        assert not _escapes(ast.parse(src).body[0]), src


class _Stalled:
    """A stand-in body for the case script's step loop whose step raises (Task 11's review, M6)."""

    def __init__(self):
        self.time = 1.5
        self.markers = SimpleNamespace(x=np.zeros((2, 3)))

    def _modes(self, tally, update=True, check_core=False):
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
    assert kept.time == 1.5 and "_modes" not in kept.__dict__


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


def _books(out):
    c = _cases()
    for row in out["steps"]:
        assert abs(row["mass_residual"]) <= c.MASS_TOL, row["t"]
        assert abs(row["enthalpy_residual"]) <= c.ENERGY_TOL, row["t"]


@pytest.mark.reference
@pytest.mark.parametrize("preset", ["dev", "production"])
def test_cold_plate_at_half_second_steps(preset):
    """The level-2 plate's centre column (8 mm deep, 700 K, 6 MW/m2, its walls the mesh box's insulated faces, the
    shared AA7075 Scheil alloy, no radiation, swell off, flow off, the liquid leaving at T_feed within the skins'
    sub-steps) at 0.5 s macro steps: the face's recession within 5 % of the one-dimensional references with the removal
    at T_feed (0.789, 2.712, 5.219 mm) and at the liquidus (0.837, 2.772, 5.298 mm) at 1.0, 1.5 and 2.0 s, and every
    step's mass to 1e-12 and enthalpy to 1e-8. Measured 2026-11-12 (dev and production alike): 0.820, 2.762 and
    5.258 mm (dev; production 0.820, 2.764, 5.252), 3.9, 1.8 and 0.7 % ahead of the first and 2.0, 0.4 and 0.8 % behind
    the second; books at most 4e-16 and 8e-15."""
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
    (phi to psi on the axis) within 10 % at 1.0, 1.5 and 2.0 s, the volume's top beneath the skin at most 1 K above
    T_rigid (the skyscrapers' criterion, HOT_BASE_TOL_K; Task 11's fix, item 4, design 20.6: the full-depth skin now
    absorbs the volume's non-rigid metal and a rigid stub below it), and the books. Measured 2026-11-12 after the fix
    (dev): recession 0.740, 1.769, 2.905 mm against 0.789, 1.817, 2.917 (6.2 % behind, reported; 2.6 and 0.4 % behind);
    non-rigid depth 1.516, 1.603, 1.704 mm against 1.487, 1.570, 1.674 (1.9, 2.1 and 1.8 % thicker); the base 2.15,
    0.48 and 2.57 K below T_rigid (before the fix 4.7 to 5.2 K above it, the skin holding only its 0.8 mm minimum);
    owed heat under 1e-4 of the heat in (2.4 to 3.6 % before)."""
    c = _cases()
    out = c.plate("hot", preset)
    assert out["failure"] is None and out["markers"] == 1
    _books(out)
    for t in c.JUDGED_TIMES:
        r = out["compare"][str(t)]
        if t > 1.25:
            assert abs(r["rel"]) <= c.RECESSION_TOL, (t, r)
        assert abs(r["nonrigid_rel"]) <= c.NONRIGID_TOL, (t, r)
        assert r["base_excess_K"] <= c.HOT_BASE_TOL_K, (t, r)          # Task 11's fix, item 4 (design 20.6)
    assert out["passes"]


BOX_FAILS = ("measured 2026-11-13 after Task 11's fix round 1 (e129c95): the markers take the projected beam over phi "
             "to round-off and the face stays even through 1.0 s (its solid tops within 0.156 mm), but the face centre "
             "recedes 0.667 mm at 1.0 s against 0.789 (15.4 % behind), with every marker full-depth in that step: the "
             "parts' base moves missed the body's volume by up to 0.9 % each, 23 mm3 net, which the volume "
             "correction took from the whole of beta (the review's I3, not settled); at 1.5 s it recedes "
             "3.006 mm against 2.712 (10.8 % ahead), the face's solid tops 2.15 mm apart; at 2.0 s beta is in nine "
             "pieces (115, 83, 15, 14 and 7 mm3 and four specks) and the centre measure reads 0.452 mm from the slab of "
             "metal left on the axis. The books close every step and no marker is trapped")


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=BOX_FAILS)
def test_q6p_box():
    """The q6p box (8 x 10 x 10 mm, 700 K, 6 MW/m2 by projection on its 10 x 10 face, 2 s at 0.5 s steps, dev cells, the
    design's meshes, the volume flow held and the skins' runoff on under the flight's deceleration, the liquid leaving at
    T_feed): the face centre's recession within 5 % of the cold plate's one-dimensional reference (T_feed) at 1.0, 1.5
    and 2.0 s, the books every step, and no marker trapped; the edges are recorded (their recession along the bisector,
    the top-zone markers). It fails (BOX_FAILS): the brief's tolerance is kept and the failure marked strict. Since
    Task 11's fix round 1 the run reaches 2.0 s (the volume correction of a speck no longer stalls)."""
    c = _cases()
    out = c.box("dev")
    assert out["failure"] is None
    _books(out)
    assert out["markers_trapped"] == 0
    for t in c.JUDGED_TIMES:
        r = out["compare"][str(t)]
        assert abs(r["rel_1d"]) <= c.RECESSION_TOL, (t, r)
    assert out["passes"]


SHOULDER_FAILS = ("measured 2026-11-13 after Task 11's fix round 1 (e129c95), judged as the controller's ruling on the "
                  "review's I2 asks: the books close and at 2.0 s a deposit of 14.2 mg lies past the heat's end, none "
                  "of it past 5 mm (the pushed plate: 3.0 %), but only 0.586 of it is rigid against 0.9: 8.4 mg of "
                  "metal above T_rigid lies in the skins within 1.83 mm past the heat's end. At 1.5 s the deposit was "
                  "40.5 mg and 0.972 rigid; in the step to 2.0 s the base moves left the body (246 mm3) holding, net, "
                  "100 mm3 more than their targets before the volume correction (the review's I3, not settled), which "
                  "took it back from all of beta (244 mm2, about 0.4 mm), phi over the never-heated metal 4 to 6 mm past "
                  "the heat's end fell by about 0.1 to 0.6 mm and the deposit read from phi's rise fell to 14.2 mg, so "
                  "the deposit cannot be judged while the base moves lean on the correction")


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=SHOULDER_FAILS)
def test_shoulder_runoff_freezes_a_deposit():
    """Runoff over a shoulder: a segment of the 50 mm sphere about its 45 degree station at the skyscrapers' pushed
    plate's loads (dp/ds -3.424e5 Pa/m and 80.29 Pa of shear down the meridian, 29.77 m/s2 of deceleration along the
    flight direction), 700 K, 6 MW/m2 by projection on the windward cap upstream of the station and nothing behind, the
    liquid kept, dev cells: the books close every step, a deposit forms past the heat's end, at least
    DEPOSIT_RIGID_SHARE of it is rigid at 2.0 s and at most the pushed plate's 3.0 % of it lies past DEPOSIT_FAR (the
    controller's ruling on the review's I2; the farthest reach of non-rigid metal is reported, not judged). The
    deposit's centroid past the heat's end is recorded against the skyscrapers' pushed plate's 1.05 mm. Measured
    2026-11-13 on e129c95: it fails (SHOULDER_FAILS): a deposit of 14.2 mg, its centroid 1.27 mm past the heat's end,
    0.586 of it rigid; books at most 1.0e-14 (mass) and 6.5e-16 (enthalpy)."""
    c = _cases()
    out = c.shoulder("dev")
    assert out["failure"] is None
    _books(out)
    last = out["steps"][-1]["deposit"]
    assert out["deposit_formed"] and last["mass_kg"] > 0.0
    assert last["rigid_share"] >= c.DEPOSIT_RIGID_SHARE                # the controller's ruling on the review's I2
    assert last["share_past_5mm"] <= c.SKYSCRAPERS_DEPOSIT_FAR_SHARE
    assert out["passes"]
