"""`python -m spheral_frag analyse` (Spheral M1 plan, Task 9, Steps 2-4) on `fake_run`s of the prepared synthetic
sphere run (tests/spheral_frag_synthetic.py::fake_run, the stand-in for the Spheral runner until M2-M3).

The fake run fills frame 0's body with a 4 mm lattice, removes at frames 1 and 2 the particles whose centres left
the body (frame 2's ten dead nose elements), and from check 1 carves a ring at 100-130 deg from the nose that moves
away (separate but not yet clear at check 1, clear at check 2) and five particles that become isolated dust. analyse
finds the ring as one fragment released at check 1 (ring flag false in 3D, its exact mass and particle count), logs
the five dust particles once each, closes the accounts to 1e-14 of the starting mass at every check, and exits 1 when
one removal's mass is dropped. Thresholds: 1e-14 m0 (plan Task 9, Step 4; Task 8 measured 2e-18), and exact equality
where the rows sum the same stored floats (math.fsum)."""
import csv
import json
import os
import subprocess
import sys

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import naming, record
from spheral_frag.__main__ import main

pytest.importorskip("pyvista")

RUN_DIR = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag", "sphere_run", syn.SPHERE_RUN_NAME)
PREP_NAME = naming.prepare_name(syn.SPHERE_RUN_NAME, 0, 2)
DX = 4e-3
SEED = 1
ACCOUNTS_RTOL = 1e-14            # plan Task 9, Step 4


@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    d = tmp_path_factory.mktemp("analyse_prep")
    syn.synthetic_material_table(str(d / "table" / "material_table.npz"))
    assert main(["prepare", "--fe-run", RUN_DIR, "--material-table", str(d / "table"), "--outdir", str(d),
                 "--no-thickness", "--quiet"]) == 0
    return str(d / PREP_NAME)


@pytest.fixture(scope="module")
def analysed(prepared, tmp_path_factory):
    runs = tmp_path_factory.mktemp("runs")
    run_dir, answers = syn.fake_run(prepared, 0, 2, DX, SEED, str(runs))
    out = tmp_path_factory.mktemp("analyse")
    assert main(["analyse", "--run", run_dir, "--outdir", str(out), "--quiet"]) == 0
    name = os.path.basename(run_dir)
    return run_dir, answers, os.path.join(str(out), name)


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def test_fake_run_is_a_run_in_the_runner_format(analysed, prepared):
    run_dir, answers, _ = analysed
    meta = record.read_meta(os.path.join(run_dir, "meta.json"))
    parsed = naming.parse_run_name(meta["run_name"])
    assert parsed["prepared"] == PREP_NAME and parsed["mode"] == "synthetic" and parsed["frames"] == (0, 2)
    assert parsed["form"] == "3d" and parsed["dx_mm"] == 4.0 and parsed["seed"] == SEED
    assert [n for n, _ in record.list_checks(run_dir)] == [0, 1, 2]
    removed = record.read_removed(os.path.join(run_dir, "removed.npz"))
    assert len(removed) == answers["n_removals"] > 0 and (removed.reason == record.REASON_CENTRE_CROSSING).all()
    assert (removed.check == 2).all()                      # only frame 2 lost elements (fixture README)
    assert answers["n_ring"] >= naming.BRACKET_DEFAULTS["min_particles"]


def test_analyse_finds_the_ring_and_the_debris(analysed):
    run_dir, answers, out = analysed
    rows = read_csv(os.path.join(out, "fragments.csv"))
    assert list(rows[0]) == record.FRAGMENT_COLUMNS and len(rows) == 1
    ring = rows[0]
    assert ring["ring"] == "false" and ring["form"] == "3d"
    assert float(ring["mass_kg"]) == answers["ring_mass_kg"] and int(ring["n_particles"]) == answers["n_ring"]
    assert float(ring["dust_mass_kg"]) == 0.0
    assert float(ring["t_release_s"]) == answers["t_release_s"] == 0.5 and int(ring["frame"]) == 1
    assert (ring["route"], ring["mechanism"]) == ("synthetic", "ring_release")
    # a full ring around the flight axis: its centroid lies behind the centre, on the axis to within one spacing
    assert float(ring["x_m"]) < 0.0 and np.hypot(float(ring["y_m"]), float(ring["z_m"])) < DX
    assert ring["run"] == os.path.basename(run_dir) and ring["window"] == "k00000-00002"
    debris = read_csv(os.path.join(out, "debris.csv"))
    assert list(debris[0]) == record.DEBRIS_COLUMNS and len(debris) == syn.FAKE_N_DEBRIS
    assert {d["kind"] for d in debris} == {"dust"} and {d["t_s"] for d in debris} == {"0.5"}
    m = record.read_check(record.check_path(run_dir, 1)).mass[0]
    assert all(float(d["mass_kg"]) == m and int(d["n_particles"]) == 1 for d in debris)
    doc = json.load(open(os.path.join(out, "analyse.json")))
    assert doc["counts"]["n_fragments"] == 1 and doc["counts"]["n_debris_rows"] == syn.FAKE_N_DEBRIS
    assert doc["released"] == [{"fragment": 0, "release_check": 1, "clear_check": 2}] and doc["unreleased"] == []
    assert doc["prepared"]["name"] == PREP_NAME and doc["run"]["meta"]["run_name"] == os.path.basename(run_dir)


def test_accounts_close_at_every_check(analysed):
    run_dir, answers, out = analysed
    doc = json.load(open(os.path.join(out, "analyse.json")))
    accounts = doc["accounts"]
    assert [a["check"] for a in accounts] == [0, 1, 2]
    assert doc["passed"] is True and doc["max_abs_rel_residual"] <= ACCOUNTS_RTOL
    for a in accounts:
        assert abs(a["rel_residual"]) <= ACCOUNTS_RTOL
    a0, a1, a2 = accounts
    assert a0["mass_fragments_kg"] == 0.0 and a0["mass_film_kg"] == 0.0
    assert a1["mass_fragments_kg"] == answers["ring_mass_kg"] and a1["mass_debris_kg"] > 0.0
    assert a2["mass_film_kg"] == answers["removed_mass_kg"] == a2["mass_film_centre_crossing_kg"]


def test_skip_force_and_a_renaming_floor(analysed, tmp_path, capsys):
    run_dir, answers, out = analysed
    outdir = os.path.dirname(out)
    capsys.readouterr()
    assert main(["analyse", "--run", run_dir, "--outdir", outdir]) == 0
    assert "analysed already" in capsys.readouterr().out
    before = open(os.path.join(out, "fragments.csv"), "rb").read()
    assert main(["analyse", "--run", run_dir, "--outdir", outdir, "--force", "--quiet"]) == 0
    assert open(os.path.join(out, "fragments.csv"), "rb").read() == before
    # a floor above the ring's particle count: the ring is debris, and the analysis is named for that floor
    floor = answers["n_ring"] + 1
    assert main(["analyse", "--run", run_dir, "--outdir", str(tmp_path), "--min-particles", str(floor),
                 "--quiet"]) == 0
    name = os.path.basename(run_dir) + "_minparticles-{}".format(floor)
    assert naming.parse_run_name(name)["brackets"]["min_particles"] == floor
    doc = json.load(open(tmp_path / name / "analyse.json"))
    assert doc["counts"]["n_fragments"] == 0 and doc["passed"] is True
    kinds = [d["kind"] for d in read_csv(tmp_path / name / "debris.csv")]
    assert kinds.count("small_group") == 1 and kinds.count("dust") == syn.FAKE_N_DEBRIS


def test_a_dropped_removal_exits_1(prepared, tmp_path, capsys):
    run_dir, answers = syn.fake_run(prepared, 0, 2, DX, SEED, str(tmp_path / "runs"), drop_removal=True)
    capsys.readouterr()
    assert main(["analyse", "--run", run_dir, "--outdir", str(tmp_path / "out"), "--quiet"]) == 1
    assert "accounts do not close" in capsys.readouterr().err
    doc = json.load(open(tmp_path / "out" / os.path.basename(run_dir) / "analyse.json"))
    m = record.read_check(record.check_path(run_dir, 0)).mass[0]
    assert doc["passed"] is False and doc["max_abs_residual_kg"] == pytest.approx(m, rel=1e-12)


def test_missing_inputs_exit_2(prepared, analysed, tmp_path):
    assert main(["analyse", "--run", str(tmp_path / "nowhere"), "--outdir", str(tmp_path)]) == 2
    run_dir, _, _ = analysed
    assert main(["analyse", "--run", run_dir, "--prepared", str(tmp_path / "noprep"), "--outdir",
                 str(tmp_path)]) == 2
    assert main(["analyse"]) == 2
    r = subprocess.run([sys.executable, "-m", "spheral_frag", "analyse", "--run", str(tmp_path / "nowhere")],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=120)
    assert r.returncode == 2 and "meta.json missing" in r.stderr
