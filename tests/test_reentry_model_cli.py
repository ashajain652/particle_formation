"""cli.py: `run` and `compare` end to end on short flights."""
import csv
import json
import os

import pytest

from reentry_model import cli, compare, sesam_io

R100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind.csv")
BASE = ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331"]


def test_run_name():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none") == \
        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"


def test_run_us76_short_flight(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--t-max", "20", "--cadence", "5", "--outdir", str(tmp_path), "--quiet"])
    assert rc == 0
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none")
    rows = list(csv.DictReader(open(tmp_path / (name + ".csv"))))
    assert [float(r["time_s"]) for r in rows] == [0.0, 5.0, 10.0, 15.0, 20.0]
    doc = json.load(open(tmp_path / (name + ".json")))
    assert doc["settings"]["atmosphere"] == "us76" and doc["results"]["end_reason"] == "t_max"
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
])
def test_bad_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2


def test_start_above_escape_altitude_exits_2(tmp_path):
    rc = cli.main(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "200",
                   "--flight-path-angle", "-1", "--atmosphere", "us76", "--outdir", str(tmp_path)])
    assert rc == 2
