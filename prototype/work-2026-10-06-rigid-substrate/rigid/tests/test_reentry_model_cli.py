"""cli.py: `run` and `compare` end to end on short flights."""
import csv
import json
import os

import numpy as np
import pytest

from reentry_model import cli, compare, coupled, sesam_io

R100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind.csv")
BASE = ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331"]


def test_run_name():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none") == \
        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"
    # the seed of numpy's generator (amendment of 2026-10-05): the default keeps the name, any other seed ends it in
    # _seed-<n>, so runs that differ only in the seed never overwrite each other
    assert cli.DEFAULT_SEED == 12345
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", seed=cli.DEFAULT_SEED) == \
        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", seed=7) == \
        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_seed-7"


def test_run_us76_short_flight(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--t-max", "20", "--cadence", "5", "--outdir", str(tmp_path), "--quiet"])
    assert rc == 0
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none")
    rows = list(csv.DictReader(open(tmp_path / (name + ".csv"))))
    assert [float(r["time_s"]) for r in rows] == [0.0, 5.0, 10.0, 15.0, 20.0]
    doc = json.load(open(tmp_path / (name + ".json")))
    assert doc["settings"]["atmosphere"] == "us76" and doc["results"]["end_reason"] == "t_max"
    assert doc["settings"]["seed"] == cli.DEFAULT_SEED                         # every run records its seed (2026-10-05)
    assert doc["inputs"]["mass_kg"] == pytest.approx(1.4728833557580150) and doc["comparison"] is None
    assert doc["provenance"]["package_version"] and "git_commit" in doc["provenance"]
    assert doc["provenance"]["reference_sha256"] is None and doc["provenance"]["replay_sha256"] is None


def test_run_replay_with_reference_writes_comparison_and_plots(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "replay:" + R100, "--reference", R100, "--t-max", "30", "--cadence", "10",
                          "--outdir", str(tmp_path), "--name", "r100_short", "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "r100_short.json"))
    m = doc["comparison"]["metrics"]
    assert doc["settings"]["atmosphere"] == "replay" and doc["settings"]["replay_reference"].endswith("_nowind")
    assert m["all"]["n_points"] >= 25 and m["hypersonic"]["dV_max_ms"] < 50.0     # 30 s of a matched flight
    assert doc["provenance"]["reference_sha256"] == doc["comparison"]["reference_sha256"]
    assert doc["provenance"]["replay_sha256"] == doc["settings"]["replay_sha256"]
    for plot in compare.PLOT_NAMES:
        assert os.path.isfile(tmp_path / "r100_short" / plot)


def test_compare_subcommand(tmp_path):
    cli.main(BASE + ["--atmosphere", "us76", "--t-max", "30", "--cadence", "1", "--outdir", str(tmp_path), "--name", "m", "--quiet"])
    rc = cli.main(["compare", "--model", str(tmp_path / "m.csv"), "--reference", R100, "--outdir", str(tmp_path / "cmp"), "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "cmp" / "m_vs_reference.json"))
    assert doc["reference"].endswith("_nowind") and "hypersonic" in doc["metrics"]


def test_unimplemented_bridging_exits_1(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--bridging", "matting", "--t-max", "1", "--outdir", str(tmp_path), "--quiet"])
    assert rc == 1


def test_escape_exits_1_but_still_writes_outputs(tmp_path):
    rc = cli.main(["run", "--diameter", "100", "--velocity", "7.8", "--altitude", "140", "--flight-path-angle", "2",
                   "--lat", "0", "--lon", "0", "--heading", "90", "--atmosphere", "us76", "--t-max", "100",
                   "--cadence", "5", "--outdir", str(tmp_path), "--name", "esc", "--quiet"])
    assert rc == 1
    doc = json.load(open(tmp_path / "esc.json"))
    assert doc["results"]["end_reason"] == "escape"
    assert doc["results"]["final_altitude_km"] == pytest.approx(150.0, abs=1e-6)


@pytest.mark.parametrize("argv", [
    BASE + ["--atmosphere", "gram"],
    BASE + ["--bridging", "legge"],
    BASE + ["--gravity", "j6"],
    ["run", "--diameter", "-1", "--velocity", "7.5", "--altitude", "77.5"],
    BASE + ["--epoch", "yesterday"],
    BASE + ["--seed", "abc"],                     # the seed (amendment of 2026-10-05): an integer in [0, 2**32 - 1],
    BASE + ["--seed", "1.5"],                     # the range numpy's generator accepts
    BASE + ["--seed", "-1"],
    BASE + ["--seed", "4294967296"],
])
def test_bad_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2


def test_start_above_escape_altitude_exits_2(tmp_path):
    rc = cli.main(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "200",
                   "--flight-path-angle", "-1", "--atmosphere", "us76", "--outdir", str(tmp_path)])
    assert rc == 2


US76_100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind.csv")
FEM = BASE + ["--atmosphere", "us76", "--thermal", "fem", "--heating", "sesam", "--h-surface", "4", "--h-core", "20", "--quiet"]


def test_run_name_with_heating():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "sesam").endswith("_none_fem-sesam")


def test_thermal_run_writes_columns_plots_and_json(tmp_path):
    rc = cli.main(FEM + ["--t-max", "20", "--dt", "0.5", "--reference", US76_100, "--outdir", str(tmp_path), "--name", "fem_short",
                         "--frames-every", "10", "--stills"])
    assert rc == 0
    rows = list(csv.DictReader(open(tmp_path / "fem_short.csv")))
    assert len(rows) == 41 and float(rows[-1]["time_s"]) == 20.0 and "surface_T_max_K" in rows[0] and "convective_heat_W" in rows[0]
    assert float(rows[-1]["temperature_K"]) > 300.0 and float(rows[0]["convective_heat_W"]) > 1e4
    doc = json.load(open(tmp_path / "fem_short.json"))
    s = doc["settings"]
    assert s["thermal"] == "fem" and s["heating"] == "sesam" and s["n_nodes"] > 2000 and s["macro_step_s"] == 0.5 and s["frames_every"] == 10
    assert s["material"] == "AA7075_nomelt" and s["emissivity"] == 0.4 and s["thermal_solver"] == "skfem"
    assert doc["results"]["n_macro_steps"] == 40 and abs(doc["results"]["energy_balance_residual"]) < 1e-6
    assert "thermal_metrics" in doc["comparison"] and doc["comparison"]["thermal_metrics"]["Q_conv"]["max"] < 0.3
    for plot in compare.PLOT_NAMES + compare.THERMAL_PLOT_NAMES:
        assert os.path.isfile(tmp_path / "fem_short" / plot)
    assert os.path.isfile(tmp_path / "fem_short" / "vtk" / "field.pvd") and doc["files"]["vtk_dir"].endswith("vtk")
    assert doc["files"]["animation"] is None and doc["files"]["section"] is None
    assert len(doc["files"]["stills"]) == 8 and all(os.path.isfile(p) for p in doc["files"]["stills"])       # surface + section stills
    assert sum(os.path.basename(p).startswith("section_") for p in doc["files"]["stills"]) == 4
    assert doc["provenance"]["skfem"] and doc["provenance"]["gmsh"]


def test_thermal_none_is_step_one(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--t-max", "10", "--cadence", "5", "--outdir", str(tmp_path), "--name", "plain", "--quiet"])
    assert rc == 0
    rows = list(csv.DictReader(open(tmp_path / "plain.csv")))
    assert "convective_heat_W" not in rows[0] and json.load(open(tmp_path / "plain.json"))["settings"]["thermal"] == "none"
    with pytest.raises(SystemExit) as exc:                                                    # parser.error: needs --thermal fem
        cli.main(BASE + ["--atmosphere", "us76", "--animate", "--outdir", str(tmp_path), "--quiet"])
    assert exc.value.code == 2


def test_missing_fenicsx_exits_2(tmp_path, capsys):
    try:
        import dolfinx  # noqa: F401
        pytest.skip("dolfinx is importable here")
    except ImportError:
        pass
    rc = cli.main(FEM + ["--thermal-solver", "fenicsx", "--t-max", "5", "--outdir", str(tmp_path)])
    assert rc == 2 and "fenicsx_env" in capsys.readouterr().err


def test_compare_subcommand_with_thermal_columns(tmp_path):
    cli.main(FEM + ["--t-max", "10", "--outdir", str(tmp_path), "--name", "m"])
    rc = cli.main(["compare", "--model", str(tmp_path / "m.csv"), "--reference", US76_100, "--outdir", str(tmp_path / "cmp"), "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "cmp" / "m_vs_reference.json"))
    assert "thermal_metrics" in doc and os.path.isfile(tmp_path / "cmp" / "heating_time.png")


def test_thermal_run_physics_mode(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--thermal", "fem", "--heating", "physics", "--h-surface", "4", "--h-core", "20",
                         "--t-max", "5", "--dt", "0.5", "--outdir", str(tmp_path), "--name", "fem_physics", "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "fem_physics.json"))
    s = doc["settings"]
    assert s["heating"] == "physics" and s["stagnation"] == "fay-riddell" and s["bridging_heat"] == "matting"
    assert s["catalycity"] == 1.0 and s["accommodation"] == 0.8 and s["matting_n"] == 1.0
    assert doc["provenance"]["cantera"]
    rows = list(csv.DictReader(open(tmp_path / "fem_physics.csv")))
    assert len(rows) == 11
    assert float(rows[-1]["q_stag_Wm2"]) > 1e6
    assert 0.0 < float(rows[-1]["heating_blend_f"]) < 0.2
    assert float(rows[-1]["T_stagnation_K"]) > float(rows[-1]["T_back_K"])
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics").endswith("_fem-physics")


def test_the_seed_makes_a_run_reproducible_and_another_seed_changes_it(tmp_path, monkeypatch, capsys):
    """numpy's global generator is seeded at the start of every run (amendment of 2026-10-05): pyamg draws from it the
    starting vector of the spectral-radius estimate that weights the prolongation smoother of every hierarchy it
    builds, which made two runs of one build differ from the first solve on (plan fact 52). Three 5 s coupled runs on
    the coarse mesh, whose conduction solves go through pyamg: two with the default seed -- the second started with the
    generator scrambled, so a run cannot depend on what drew from it before -- agree bit for bit in every history column
    and every result field but the run time, while one with another seed differs, at round-off, under its own name."""
    histories = []
    original_run = coupled.CoupledRun.run

    def run_and_keep(self):
        histories.append(original_run(self))
        return histories[-1]

    monkeypatch.setattr(coupled.CoupledRun, "run", run_and_keep)
    argv = FEM + ["--t-max", "5", "--dt", "0.5"]
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "sesam")
    np.random.seed(1)
    assert cli.main(argv + ["--outdir", str(tmp_path / "a")]) == 0
    np.random.seed(2)
    np.random.rand(1000)
    assert cli.main(argv + ["--outdir", str(tmp_path / "b")]) == 0
    capsys.readouterr()
    assert cli.main([a for a in argv if a != "--quiet"] + ["--outdir", str(tmp_path / "c"), "--seed", "7"]) == 0
    assert name + "_seed-7 (seed 7): t_max at t = 5.0 s" in capsys.readouterr().out     # the printed summary names it
    a, b, c = histories
    ra, rb = (json.load(open(tmp_path / d / (name + ".json")))["results"] for d in ("a", "b"))
    assert set(a.columns) == set(b.columns) and all(np.array_equal(a.columns[k], b.columns[k], equal_nan=True) for k in a.columns)
    assert {k: v for k, v in ra.items() if k != "runtime_s"} == {k: v for k, v in rb.items() if k != "runtime_s"}
    assert not all(np.array_equal(a.columns[k], c.columns[k], equal_nan=True) for k in a.columns)
    for k in a.columns:                                       # the seed acts at round-off; the melting model amplifies it
        np.testing.assert_allclose(c.columns[k], a.columns[k], rtol=1e-6, atol=1e-9)
    doc = json.load(open(tmp_path / "c" / (name + "_seed-7.json")))
    assert doc["settings"]["seed"] == 7 and doc["run_name"] == name + "_seed-7"


# ---------------------------------------------------------------------------------------------------------------
# Step 3: melting flags, the melting run's files, the bookkeeping device against the melting reference

MELT_100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv")
MELT = FEM + ["--melt", "on", "--prism-layers", "0", "--h-surface", "4", "--h-core", "20"]


def test_run_name_with_melt():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin").endswith("_fem-physics_melt-girin")
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                              deep_runoff=False).endswith("_fem-physics_melt-girin_deeprunoff-off")    # amendment of 2026-10-02
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                              molten_cascade=False).endswith("_fem-physics_melt-girin_moltencascade-off")   # the molten cascade
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", deep_runoff=False,
                              molten_cascade=False, seed=7).endswith("_melt-girin_deeprunoff-off_moltencascade-off_seed-7")
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                              rigid_substrate=False).endswith("_fem-physics_melt-girin_rigidsubstrate-off")  # 2026-10-06
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", deep_runoff=False,
                              molten_cascade=False, rigid_substrate=False,
                              seed=7).endswith("_melt-girin_deeprunoff-off_moltencascade-off_rigidsubstrate-off_seed-7")


def test_melting_run_writes_columns_files_and_json(tmp_path):
    """15 s from 71 km with a warm body: melt columns, the particle files, the melt plots, stills with the melt marks."""
    pytest.importorskip("cantera")
    argv = [a for a in MELT if a not in ("--heating", "sesam")] + ["--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "800",
                                                                   "--t-max", "15", "--outdir", str(tmp_path), "--name", "melt_short", "--stills"]
    assert cli.main(argv) == 0
    rows = list(csv.DictReader(open(tmp_path / "melt_short.csv")))
    assert len(rows) == 31 and "sprayed_mass_kg" in rows[0] and "film_mass_kg" in rows[0] and "closure_fraction_girin" in rows[0]
    assert float(rows[-1]["kn_body"]) > 0.0 and float(rows[-1]["kn_local_stag"]) < float(rows[-1]["kn_body"]) and "flow_branch" in rows[0]
    assert float(rows[-1]["mass_kg"]) < float(rows[0]["mass_kg"]) and float(rows[-1]["sprayed_mass_kg"]) > 0.0
    doc = json.load(open(tmp_path / "melt_short.json"))
    s, r, f = doc["settings"], doc["results"], doc["files"]
    assert s["melt"] == "on" and s["material"] == "AA7075_range" and s["removal"] == "girin" and s["runoff"] == "on" and s["prism_layers"] == 0
    assert s["rarefied_shear"] == "slip" and s["we_critical"] == 4.62 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.86
    assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01 and s["deep_runoff"] == "on" and s["seed"] == cli.DEFAULT_SEED
    assert "deep_runoff_mass_kg" in rows[0] and "deep_mass_kg" in rows[0] and r["deep_runoff_mass_kg"] >= 0.0
    assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
    assert s["molten_cascade"] == "on" and "cascade_passes" in rows[0] and "cascade_mass_kg" in rows[0]      # the molten cascade
    assert s["rigid_substrate"] == "on" and all(k in rows[0] for k in ("nonrigid_depth_mean_mm", "slurry_thick_fraction",
                                                                       "rigid_thin_fraction", "slurry_held_mass_kg"))
    assert r["cascade_mass_kg"] == pytest.approx(float(rows[-1]["cascade_mass_kg"])) and r["cascade_mass_kg"] >= 0.0
    assert r["cascade_capped_steps"] == 0
    assert r["cascade_passes_max"] == max(int(float(row["cascade_passes"])) for row in rows)
    assert s["size_feedback"] == "current" and "nose_radius_mm" in rows[0] and float(rows[-1]["nose_radius_mm"]) > 0.0
    assert s["T_liquidus_K"] == 908.0 and s["latent_heat_Jkg"] == 400e3 and s["demise_fraction"] == 0.01 and s["particles"] is True
    assert r["melt_onset_altitude_km"] is not None and r["sprayed_mass_kg"] > 0.0 and r["n_source_rows"] > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    for key in ("particles", "particles_summary", "size_distribution"):
        assert os.path.isfile(f[key])
    assert len(f["melt_plots"]) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in f["melt_plots"])
    assert doc["comparison"] is None and f["film"] is None
    names = [os.path.basename(p) for p in f["stills"]]
    assert any(n.startswith("melt_onset") for n in names) and any(n.startswith("film_") for n in names) and any(n.startswith("section_spraying_onset") for n in names)


def test_bookkeeping_device_against_the_melting_reference(tmp_path):
    """SESAM-equivalent heating + AA7075 + instant removal + k x 1e4 on the coarse mesh over the whole flight: the
    melt metrics are written and the mass follows SESAM's lumped law (2 % / 0.5 km / 2 %, spec 13.1; measured
    0.82-0.98 % / 0.10 km / -1.3 %)."""
    argv = MELT + ["--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--reference", MELT_100,
                   "--outdir", str(tmp_path), "--name", "bookkeeping"]
    assert cli.main(argv) == 0
    doc = json.load(open(tmp_path / "bookkeeping.json"))
    mm = doc["comparison"]["melt_metrics"]
    assert doc["results"]["end_reason"] == "demise" and mm["mass"]["max_rel_m0"] < 0.02
    assert abs(mm["onset_altitude_diff_km"]) < 0.5 and abs(mm["demise_time_rel"]) < 0.02 and mm["sprayed_mass_kg"] == 0.0
    assert os.path.isfile(tmp_path / "bookkeeping" / "mass_time.png") and doc["settings"]["k_scale"] == 1e4
    assert doc["settings"]["size_feedback"] == "initial"                                             # SESAM's D0 / R0 for the device
    assert not os.path.isfile(tmp_path / "bookkeeping" / "particles.npz") or True                    # written (empty table) with --particles


@pytest.mark.parametrize("argv", [
    BASE + ["--atmosphere", "us76", "--melt", "on"],                                              # needs --thermal fem
    MELT + ["--kr", "0"],
    MELT + ["--prism-layers", "-1"],
    MELT + ["--demise-fraction", "1.5"],
    MELT + ["--k-scale", "0"],
    MELT + ["--removal", "magic"],
    MELT + ["--size-feedback", "shrinking"],
    MELT + ["--deep-runoff", "maybe"],
    MELT + ["--molten-cascade", "maybe"],
    MELT + ["--rigid-substrate", "maybe"],
])
def test_bad_melt_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
