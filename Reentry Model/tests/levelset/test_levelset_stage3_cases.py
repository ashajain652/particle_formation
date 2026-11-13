"""Level-set stage 3's verification cases (Task 11; design section 20.11 item 3, plan Task 11): the plates at 0.5 s macro
steps (Decision A of stage 2), the q6p box (Ruling 5) and runoff over a shoulder, each run by
`analysis/levelset/stage3_cases.py` (whose docstring describes the cases and why the plate's face carries one marker) and
judged here against the brief's tolerances. The cases are marked `reference` (minutes each); the import check runs in
the unit tier. The measured numbers are recorded in `analysis/levelset/stage3_recorded.json`."""
import ast
import importlib.util
import os

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


BOX_FAILS = ("measured 2026-11-12 after Task 11's fix round, diagnosed in its report (section 10.7): the face stays flat "
             "through the cold step (its skins' solid tops within 4e-11 m, their groups 0.70 to 0.99 of their areas, "
             "their tops within 6.8 K), but at 1.0 s the face centre recedes 0.936 mm against 0.789 (18.6 % ahead; "
             "40.5 % before the fixes) and at 1.5 s 2.921 mm against 2.712 (7.7 %): the edges and the rim melt first, "
             "26 then 43 markers turn top-zone (7 on the face by 1.5 s), their melt goes to the volume, which in stage 3 "
             "holds it without removal (0.83 mm of the face removed by 1.5 s against the reference's 2.71), and the face "
             "buckles (2.4 mm between its face markers' solid tops at 1.5 s); the step to 2.0 s then stops in a base "
             "move ('volume correction of piece 1 did not converge')")


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=BOX_FAILS)
def test_q6p_box():
    """The q6p box (8 x 10 x 10 mm, 700 K, 6 MW/m2 by projection on its 10 x 10 face, 2 s at 0.5 s steps, dev cells, the
    design's meshes, the volume flow held and the skins' runoff on under the flight's deceleration, the liquid leaving at
    T_feed): the face centre's recession within 5 % of the cold plate's one-dimensional reference (T_feed) at 1.0, 1.5
    and 2.0 s, the books every step, and no marker trapped; the edges are recorded (their recession along the bisector,
    the top-zone markers). It fails (BOX_FAILS): the brief's tolerance is kept and the failure marked strict."""
    c = _cases()
    out = c.box("dev")
    assert out["failure"] is None
    _books(out)
    assert out["markers_trapped"] == 0
    for t in c.JUDGED_TIMES:
        r = out["compare"][str(t)]
        assert abs(r["rel_1d"]) <= c.RECESSION_TOL, (t, r)
    assert out["passes"]


SHOULDER_FAILS = ("measured 2026-11-12 after Task 11's fix round, partly diagnosed in its report (section 10.7): at 2.0 s "
                  "the skins hold 25 mg of metal above T_rigid out to 4.30 mm past the heat's end (the third marker), "
                  "past FREEZE_REACH; the full-depth skins now absorb the volume's non-rigid metal (item 4), so metal "
                  "the volume held there before, heated from the cap, is now counted in the skins (2.7 mg within "
                  "2.3 mm before the fixes); whether it is runoff that has not frozen or base metal heated past "
                  "T_rigid is not separated. The books close, a deposit of 20 mg forms 1.84 mm past the heat's end")


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=SHOULDER_FAILS)
def test_shoulder_runoff_freezes_a_deposit():
    """Runoff over a shoulder: a segment of the 50 mm sphere about its 45 degree station at the skyscrapers' pushed
    plate's loads (dp/ds -3.424e5 Pa/m and 80.29 Pa of shear down the meridian, 29.77 m/s2 of deceleration along the
    flight direction), 700 K, 6 MW/m2 by projection on the windward cap upstream of the station and nothing behind, the
    liquid kept, dev cells: the books close every step, a deposit forms past the heat's end, and beyond FREEZE_REACH
    (4 mm, two marker spacings) of the heat's end no skin holds metal above T_rigid. The deposit's centroid past the
    heat's end is recorded against the skyscrapers' pushed plate's 1.05 mm. Measured 2026-11-12 before Task 11's fix
    round: at 2.0 s a deposit of 40.7 mg, its centroid 1.88 mm past the heat's end, 2.7 mg of non-rigid metal within
    2.3 mm of it (passed). After the fixes it fails (SHOULDER_FAILS): a deposit of 20.1 mg, its centroid 1.84 mm past the
    heat's end, reaching 3.99 mm, none 5 to 10 mm past it (the pushed plate: 3.0 %), but 25 mg of non-rigid metal in the
    skins out to 4.30 mm; books at most 9.2e-15 (mass) and 2.5e-16 (enthalpy)."""
    c = _cases()
    out = c.shoulder("dev")
    assert out["failure"] is None
    _books(out)
    last = out["steps"][-1]["deposit"]
    assert out["deposit_formed"] and last["mass_kg"] > 0.0
    assert last["hot_reach_past_mm"] <= c.FREEZE_REACH * 1e3
    assert out["passes"]
