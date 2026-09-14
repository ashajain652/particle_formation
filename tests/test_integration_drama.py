"""Integration tests that run the real SESAM through pyDRAMA (marker: drama)."""
import csv
import json
import math
import os
import subprocess
import sys
from datetime import datetime

import pytest

import sphere_reentry as sr
from helpers import REPO_ROOT, PY

pytestmark = pytest.mark.drama

FAP_DAY = sr.read_lines(os.path.join(REPO_ROOT, "data", "fap_day.dat"))
FAP_MON = sr.read_lines(os.path.join(REPO_ROOT, "data", "fap_mon.dat"))

# parent state at 7.5 km/s (spike case T2/T5) and at 0.3 km/s (spike case E3)
ENTRY_STATE = dict(altitude_km=77.5, lat_deg=29.546, lon_deg=-82.134, flight_path_deg=-0.959,
                   heading_deg=347.168, epoch=datetime(2024, 8, 1, 12, 53, 7))
LOW_STATE = dict(altitude_km=33.229, lat_deg=35.99, lon_deg=-83.90, flight_path_deg=-60.2,
                 heading_deg=347.3, epoch=datetime(2024, 8, 1, 12, 56, 17))


def first_csv_row(path):
    with open(path) as fh:
        return next(csv.DictReader(fh))


def test_real_run_50mm_at_entry_speed_demises(tmp_path):
    run = sr.SphereRun(velocity_kms=7.5, temperature_K=300.0, diameter_mm=50.0, **ENTRY_STATE)
    doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, FAP_DAY, FAP_MON)
    assert doc["status"] == "ok", doc["error"]
    r = doc["results"]
    for key in ("max_temperature_K", "final_mass_kg", "time_at_melting_temperature_s",
                "final_velocity_kms", "final_radius_mm"):
        assert r[key] is not None and math.isfinite(r[key]), key
    assert r["max_temperature_K"] == 850.0
    assert r["outcome"] == "demised" and r["end_of_life_reason"] in ("ballooning", "uncritical")
    assert r["time_at_melting_temperature_s"] > 0.0
    assert r["final_radius_mm"] < 1.0
    assert not (tmp_path / "raw" / doc["run_name"]).exists()
    assert doc["provenance"]["sesam_version"] == "2.3.0"
    row = first_csv_row(doc["files"]["csv"])
    assert float(row["time_s"]) == 0.0
    assert float(row["altitude_km"]) == pytest.approx(77.5, abs=1e-3)
    assert float(row["velocity_kms"]) == pytest.approx(7.5, abs=1e-3)
    assert float(row["temperature_K"]) == 300.0
    assert float(row["lat_deg"]) == pytest.approx(29.546, abs=1e-3)
    assert float(row["lon_deg"]) == pytest.approx(-82.134, abs=1e-3)
    assert float(row["flight_path_deg"]) == pytest.approx(-0.959, abs=1e-3)
    assert float(row["heading_deg"]) == pytest.approx(347.168, abs=1e-3)


def test_real_run_5mm_low_altitude_survives_to_ground(tmp_path):
    run = sr.SphereRun(velocity_kms=0.3, temperature_K=300.0, diameter_mm=5.0, **LOW_STATE)
    doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, FAP_DAY, FAP_MON)
    assert doc["status"] == "ok", doc["error"]
    r = doc["results"]
    assert r["outcome"] == "survived" and r["end_of_life_reason"] == "ground impact"
    assert r["final_mass_source"] == "impacting_fragments_xml"
    assert r["final_mass_kg"] == pytest.approx(doc["inputs"]["initial_mass_kg"], abs=1e-9)
    assert r["final_radius_mm"] == pytest.approx(2.5, abs=1e-6)
    assert r["final_altitude_km"] == pytest.approx(0.0, abs=1e-3)
    assert r["max_temperature_K"] < 320.0
    assert r["n_rows"] > 100
    assert any("coarse" in w for w in doc["warnings"])


def test_cli_real_run_and_keep_raw(tmp_path):
    cmd = [PY, os.path.join(REPO_ROOT, "sphere_reentry.py"),
           "--velocity", "0.3", "--altitude", "33.229", "--temperature", "300", "--diameter", "5",
           "--flight-path-angle", "-60.2", "--heading", "347.3", "--lat", "35.99", "--lon", "-83.90",
           "--epoch", "2024-08-01T12:56:17", "--outdir", str(tmp_path / "runs"),
           "--raw-dir", str(tmp_path / "raw"), "--keep-raw"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    name = "sphere_d005.00mm_T0300.0K_v00.30000kms_h033.229km"
    assert (tmp_path / "runs" / (name + ".csv")).is_file()
    doc = json.load(open(tmp_path / "runs" / (name + ".json")))
    assert doc["status"] == "ok" and doc["files"]["raw_dir"] == str(tmp_path / "raw" / name)
    assert (tmp_path / "raw" / name / "run_0" / "reentry" / "sesam.log").is_file()
    assert "survived" in proc.stdout
