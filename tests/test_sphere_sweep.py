"""Unit tests for sphere_sweep.py (no DRAMA needed)."""
import copy
import csv
import hashlib
import json
import os
from datetime import datetime

import pytest

import sphere_reentry as sr
import sphere_sweep as sw
from helpers import FIXTURES

MINI = os.path.join(FIXTURES, "mini_dmf_output.json")


@pytest.fixture
def parent():
    return sw.load_parent(MINI, "Main Body")


class TestLoadParent:
    def test_epoch_object_and_branch(self, parent):
        assert parent.epoch == datetime(2024, 8, 1, 12, 0, 0)
        assert parent.object_name == "Main Body"
        assert len(parent.rows) == 14
        assert parent.branch_start == 3
        assert parent.v_max == 7.905 and parent.v_min == 0.05
        assert parent.v_min_row == 13
        assert [r["velocity"] for r in parent.branch] == [7.905, 7.9, 7.8, 7.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]

    def test_sha256_matches_file(self, parent):
        with open(MINI, "rb") as fh:
            assert parent.sha256 == hashlib.sha256(fh.read()).hexdigest()
        assert parent.path == MINI

    def test_describe_v_min(self, parent):
        assert parent.describe_v_min() == "v_min = 0.050000 km/s at t = 1300.000 s, altitude 0.000 km (parent row 13)"

    def test_other_object_by_name_with_spaces(self):
        other = sw.load_parent(MINI, "Other Part")
        assert len(other.rows) == 3 and other.branch_start == 0 and other.v_max == 7.5

    def test_unknown_object_lists_available_keys(self):
        with pytest.raises(ValueError) as exc:
            sw.load_parent(MINI, "Nope")
        assert "Nope" in str(exc.value) and "Main_Body" in str(exc.value)

    def test_two_satellites_is_an_error(self, tmp_path):
        doc = json.load(open(MINI))
        doc["singleModuleOutputs"]["satellites"].append(copy.deepcopy(doc["singleModuleOutputs"]["satellites"][0]))
        p = tmp_path / "two.json"
        p.write_text(json.dumps(doc))
        with pytest.raises(ValueError, match="exactly one satellite"):
            sw.load_parent(str(p), "Main Body")

    def test_branch_with_one_row_is_an_error(self, tmp_path):
        doc = json.load(open(MINI))
        key = [k for k in doc["singleModuleOutputs"]["satellites"][0]["missionPhases"][0]["epochs"][0]["analysisModules"][0]["results"] if "Main_Body" in k and "Trajectory" in k][0]
        rows = doc["singleModuleOutputs"]["satellites"][0]["missionPhases"][0]["epochs"][0]["analysisModules"][0]["results"][key]
        rows[-1]["velocity"] = 9.0          # velocity increases at the very end -> branch of 1 row
        p = tmp_path / "bad.json"
        p.write_text(json.dumps(doc))
        with pytest.raises(ValueError, match="fewer than 2 rows"):
            sw.load_parent(str(p), "Main Body")


class TestDescendingBranch:
    def test_start_index(self):
        assert sw.descending_branch_start([7.9, 7.8, 7.7]) == 0
        assert sw.descending_branch_start([7.9, 7.8, 7.85, 7.7, 7.6]) == 2
        assert sw.descending_branch_start([7.9, 7.95, 7.8, 7.85, 7.7]) == 3
        assert sw.descending_branch_start([1.0]) == 0

    def test_plateau_is_not_mistaken_for_an_increase(self):
        # A trailing plateau (equal consecutive velocities, as SESAM prints to 3
        # decimals near terminal velocity) must never register as an "increase":
        # with a strict `>` comparison the whole non-increasing list is one branch.
        assert sw.descending_branch_start([7.9, 7.8, 7.7, 0.05, 0.05, 0.05]) == 0
        # A genuine increase earlier still moves branch_start forward as usual,
        # and a plateau right after it must not push branch_start any further.
        assert sw.descending_branch_start([7.9, 7.95, 7.8, 7.7, 7.7, 7.7]) == 1


class TestInterpolation:
    def test_angle_lerp_and_wrap(self):
        assert sw.angle_lerp(350.0, 10.0, 0.5) == pytest.approx(0.0, abs=1e-9)
        assert sw.angle_lerp(90.0, 270.0, 0.0) == pytest.approx(90.0)
        assert sw.angle_lerp(10.0, 20.0, 0.5) == pytest.approx(15.0)
        assert sw.wrap_longitude(190.0) == -170.0
        assert sw.wrap_longitude(-190.0) == 170.0
        assert sw.wrap_longitude(180.0) == -180.0
        assert sw.wrap_longitude(-83.9) == pytest.approx(-83.9)
        assert sw.lerp(10.0, 20.0, 0.25) == 12.5

    def test_midpoint_between_rows(self, parent):
        st = sw.interpolate_state(parent, 7.4)          # between 7.8 (t=500) and 7.0 (t=600)
        assert st["altitude_km"] == pytest.approx(90.0)
        assert st["lat_deg"] == pytest.approx(-5.0)
        assert st["flight_path_deg"] == pytest.approx(-0.75)
        assert st["time_s"] == pytest.approx(550.0)
        assert min(st["heading_deg"], 360.0 - st["heading_deg"]) == pytest.approx(0.0, abs=1e-9)  # 350 -> 10 wraps through 0
        assert st["lon_deg"] == pytest.approx(149.0)
        assert st["epoch"] == datetime(2024, 8, 1, 12, 9, 10)

    def test_longitude_wraps_across_antimeridian(self, parent):
        st = sw.interpolate_state(parent, 6.0)          # lon 178 -> -178
        assert abs(st["lon_deg"]) == pytest.approx(180.0)
        assert -180.0 <= st["lon_deg"] < 180.0

    def test_exact_row_velocities(self, parent):
        assert sw.interpolate_state(parent, 5.0)["altitude_km"] == pytest.approx(60.0)
        assert sw.interpolate_state(parent, 7.905)["altitude_km"] == pytest.approx(128.0)
        assert sw.interpolate_state(parent, 0.05)["altitude_km"] == pytest.approx(0.0)
        assert sw.interpolate_state(parent, 0.05)["time_s"] == pytest.approx(1300.0)

    def test_outside_range_is_none(self, parent):
        assert sw.interpolate_state(parent, 7.906) is None
        assert sw.interpolate_state(parent, 0.049) is None
        assert sw.interpolate_state(parent, 7.905) is not None
        assert sw.interpolate_state(parent, 0.05) is not None

    def test_orbital_phase_is_never_sampled(self, parent):
        # 7.902 km/s occurs in the orbital phase (rows 1-2) AND on the branch (rows 3-4);
        # only the branch is used: altitude between 128 and 120, not between 124 and 126.
        st = sw.interpolate_state(parent, 7.902)
        assert 120.0 < st["altitude_km"] < 128.0 and st["time_s"] > 300.0

    def test_velocity_plateau_maps_to_the_first_highest_altitude_row(self):
        # Mimics the real parent's terminal-velocity plateau: SESAM prints velocity
        # to 3 decimals, so several consecutive rows near impact share one value
        # (here 0.028 km/s, rows 1-3). interpolate_state must take the FIRST
        # bracketing pair, i.e. resolve to the plateau's highest-altitude / earliest
        # row (row 1: altitude 0.153, t=100), not a later row of the same plateau
        # (row 2: altitude 0.080, t=110) or the last one (row 3: altitude 0.012, t=120).
        rows = [
            {"time": 0.0,   "altitude": 1.000, "lat": 10.0, "lon": 20.0, "velocity": 0.030, "path": -80.0, "heading": 40.0},
            {"time": 100.0, "altitude": 0.153, "lat": 10.1, "lon": 20.1, "velocity": 0.028, "path": -85.0, "heading": 41.0},
            {"time": 110.0, "altitude": 0.080, "lat": 10.2, "lon": 20.2, "velocity": 0.028, "path": -86.0, "heading": 42.0},
            {"time": 120.0, "altitude": 0.012, "lat": 10.3, "lon": 20.3, "velocity": 0.028, "path": -87.0, "heading": 43.0},
            {"time": 130.0, "altitude": 0.000, "lat": 10.4, "lon": 20.4, "velocity": 0.012, "path": -90.0, "heading": 44.0},
        ]
        assert sw.descending_branch_start([r["velocity"] for r in rows]) == 0
        plateau_parent = sw.ParentTrajectory(
            path="synthetic", sha256="deadbeef", object_name="Synthetic",
            epoch=datetime(2024, 1, 1, 0, 0, 0), rows=rows, branch_start=0)

        st = sw.interpolate_state(plateau_parent, 0.028)

        assert st["altitude_km"] == pytest.approx(0.153)
        assert st["lat_deg"] == pytest.approx(10.1)
        assert st["lon_deg"] == pytest.approx(20.1)
        assert st["flight_path_deg"] == pytest.approx(-85.0)
        assert st["heading_deg"] == pytest.approx(41.0)
        assert st["time_s"] == pytest.approx(100.0)
        assert st["epoch"] == datetime(2024, 1, 1, 0, 1, 40)


class TestGrid:
    def test_velocity_grid(self):
        grid = sw.velocity_grid(0.05)
        assert len(grid) == 100 and grid[0] == 7.5 and grid[-1] == pytest.approx(0.05)
        assert all(a > b for a, b in zip(grid, grid[1:]))
        assert grid[0] - grid[1] == pytest.approx((7.5 - 0.05) / 99)
        assert all(isinstance(v, float) for v in grid)

    def test_default_diameters_and_temperatures(self):
        assert sw.DIAMETERS_MM == [float(d) for d in range(5, 101, 5)] and len(sw.DIAMETERS_MM) == 20
        assert sw.TEMPERATURES_K[0] == 300.0 and sw.TEMPERATURES_K[-1] == 750.0 and len(sw.TEMPERATURES_K) == 46
        assert sw.V_TOP_KMS == 7.5 and sw.N_VELOCITIES == 100

    def test_parse_float_list(self):
        assert sw.parse_float_list("5, 10,15") == [5.0, 10.0, 15.0]
        with pytest.raises(Exception):
            sw.parse_float_list("")

    def test_format_arg(self):
        assert sw.format_arg(7.5) == "7.500000"
        assert sw.format_arg(0.0283333333) == "0.028333"
        assert sw.format_arg(-82.134) == "-82.134000"


class TestBuildPoints:
    def test_states_skips_and_ordering(self, parent):
        points = sw.build_points(parent, [10.0, 5.0], [310.0, 300.0], [7.4, 9.0, 0.01])
        assert [(p.diameter_mm, p.temperature_K, p.velocity_kms) for p in points] == [
            (5.0, 300.0, 9.0), (5.0, 300.0, 7.4), (5.0, 300.0, 0.01),
            (5.0, 310.0, 9.0), (5.0, 310.0, 7.4), (5.0, 310.0, 0.01),
            (10.0, 300.0, 9.0), (10.0, 300.0, 7.4), (10.0, 300.0, 0.01),
            (10.0, 310.0, 9.0), (10.0, 310.0, 7.4), (10.0, 310.0, 0.01)]
        skipped, ok = points[0], points[1]
        assert skipped.status == "skipped" and skipped.skip_reason == "parent never reaches velocity"
        assert skipped.run_name is None and skipped.altitude_km is None and skipped.epoch_utc is None
        assert ok.status == "pending" and ok.skip_reason is None
        assert ok.altitude_km == pytest.approx(90.0) and ok.flight_path_deg == pytest.approx(-0.75)
        assert ok.lat_deg == pytest.approx(-5.0) and ok.lon_deg == pytest.approx(149.0)
        assert min(ok.heading_deg, 360.0 - ok.heading_deg) == pytest.approx(0.0, abs=1e-9)
        assert ok.epoch_utc == "2024-08-01T12:09:10"
        assert ok.run_name == "sphere_d005.00mm_T0300.0K_v07.40000kms_h090.000km"
        assert points[2].status == "skipped"

    def test_run_name_matches_script1_on_formatted_args(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        run = sw.sphere_run_for(point)
        assert run.velocity_kms == 7.4 and run.altitude_km == 90.0 and run.epoch == datetime(2024, 8, 1, 12, 9, 10)
        assert run.flight_path_deg == float("-0.750000") and run.heading_deg % 360.0 == 0.0
        assert sr.run_name(run) == point.run_name

    def test_limit(self, parent):
        assert len(sw.build_points(parent, [5.0, 10.0], [300.0], [7.4, 5.0], limit=3)) == 3

    def test_full_default_matrix_size(self, parent):
        points = sw.build_points(parent, sw.DIAMETERS_MM, sw.TEMPERATURES_K, sw.velocity_grid(parent.v_min))
        assert len(points) == 20 * 46 * 100
        assert sum(p.status == "skipped" for p in points) == 0
        assert len({p.run_name for p in points}) == len(points)

    def test_colliding_run_names_from_close_velocity_overrides_raise(self, parent):
        # 7.4 and 7.4000001 round to the same run_name under format_arg's 6-decimal / run_name's
        # 5-decimal rounding, but interpolate_state() is called with the raw, un-rounded velocities
        # -- so under the thread-pool executor these two points would run CONCURRENTLY and clobber
        # each other's output files instead of a safe sequential overwrite.
        with pytest.raises(ValueError) as excinfo:
            sw.build_points(parent, [5.0], [300.0], [7.4, 7.4000001])
        assert "sphere_d005.00mm_T0300.0K_v07.40000kms_h090.000km" in str(excinfo.value)


def write_run_json(runs_dir, run_name, status="ok", results=None):
    os.makedirs(runs_dir, exist_ok=True)
    doc = {"schema_version": 1, "run_name": run_name, "status": status, "results": results}
    with open(os.path.join(runs_dir, run_name + ".json"), "w") as fh:
        json.dump(doc, fh)


class TestManifest:
    def test_manifest_document(self, parent):
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 9.0])
        grid = {"diameters_mm": [5.0], "temperatures_K": [300.0], "velocities_kms": [7.4, 9.0]}
        settings = {"timeout_s": 600, "batch_size": None, "sphere_reentry_version": sr.SCRIPT_VERSION}
        doc = sw.manifest_document(parent, grid, settings, points, created_utc="2026-09-14T00:00:00Z")
        assert doc["created_utc"] == "2026-09-14T00:00:00Z" and doc["updated_utc"].endswith("Z")
        assert doc["parent"] == {
            "path": os.path.abspath(MINI), "sha256": parent.sha256, "object": "Main Body",
            "epoch_utc": "2024-08-01T12:00:00", "branch_first_row": 3, "branch_last_row": 13,
            "v_min_kms": 0.05, "v_min_time_s": 1300.0, "v_min_altitude_km": 0.0, "v_max_kms": 7.905}
        assert doc["grid"] == grid and doc["settings"] == settings
        assert [p["status"] for p in doc["points"]] == ["skipped", "pending"]
        assert doc["points"][1]["run_name"] == "sphere_d005.00mm_T0300.0K_v07.40000kms_h090.000km"
        assert set(doc["points"][1]) == set(sw.Point.__dataclass_fields__)

    def test_write_and_load_round_trip(self, parent, tmp_path):
        points = sw.build_points(parent, [5.0], [300.0], [7.4])
        doc = sw.manifest_document(parent, {}, {}, points)
        path = str(tmp_path / "sweep_manifest.json")
        sw.write_manifest(path, doc)
        assert sw.load_manifest(path) == json.loads(json.dumps(doc))
        assert not os.path.exists(path + ".tmp")


class TestResume:
    def test_statuses_from_existing_json(self, parent, tmp_path):
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0, 9.0])   # ordered 9.0, 7.4, 5.0, 2.0
        skipped, ok, failed, fresh = points
        write_run_json(runs, ok.run_name, "ok")
        write_run_json(runs, failed.run_name, "error")
        sw.apply_resume(points, runs)
        assert (ok.status, failed.status, fresh.status, skipped.status) == ("done", "failed", "pending", "skipped")

    def test_force_and_retry_failed(self, parent, tmp_path):
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0])
        write_run_json(runs, points[0].run_name, "ok")
        write_run_json(runs, points[1].run_name, "timeout")
        sw.apply_resume(points, runs, retry_failed=True)
        assert [p.status for p in points] == ["done", "pending", "pending"]
        sw.apply_resume(points, runs, force=True)
        assert [p.status for p in points] == ["pending", "pending", "pending"]

    def test_corrupt_json_counts_as_failed(self, parent, tmp_path):
        runs = tmp_path / "runs"
        runs.mkdir()
        points = sw.build_points(parent, [5.0], [300.0], [7.4])
        (runs / (points[0].run_name + ".json")).write_text("{not json")
        sw.apply_resume(points, str(runs))
        assert points[0].status == "failed"


class TestSummary:
    def test_columns_are_the_spec_list(self):
        assert sw.SUMMARY_COLUMNS == [
            "diameter_mm", "initial_temperature_K", "initial_velocity_kms", "initial_altitude_km",
            "flight_path_angle_deg", "heading_deg", "latitude_deg", "longitude_deg", "epoch_utc",
            "status", "skip_reason", "max_temperature_K", "final_mass_kg", "final_mass_source",
            "mass_loss_fraction", "time_at_melting_temperature_s", "final_velocity_kms",
            "final_radius_mm", "final_altitude_km", "end_of_life_reason", "outcome", "wall_time_s", "run_name"]

    def test_rows_for_each_status(self, parent, tmp_path):
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0, 9.0])   # ordered 9.0, 7.4, 5.0, 2.0
        skipped, done, failed, pending = points
        write_run_json(runs, done.run_name, "ok", results={
            "max_temperature_K": 850.0, "final_mass_kg": 3e-06, "final_mass_source": "sesam_log_event_end",
            "mass_loss_fraction": 0.98, "time_at_melting_temperature_s": 6.774, "final_velocity_kms": 2.998,
            "final_radius_mm": 0.633, "final_altitude_km": 76.548, "end_of_life_reason": "uncritical",
            "outcome": "demised"})
        done.status, done.wall_time_s = "done", 0.31
        failed.status, failed.returncode = "failed", 1
        path = str(tmp_path / "sweep_summary.csv")
        sw.write_summary(path, points, runs)
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            assert reader.fieldnames == sw.SUMMARY_COLUMNS
            rows = list(reader)
        assert len(rows) == 4
        s = rows[0]
        assert s["status"] == "skipped" and s["skip_reason"] == "parent never reaches velocity"
        assert s["run_name"] == "" and s["initial_altitude_km"] == "" and s["initial_velocity_kms"] == "9.0"
        d = rows[1]
        assert d["status"] == "done" and d["run_name"] == done.run_name and d["wall_time_s"] == "0.31"
        assert d["diameter_mm"] == "5.0" and d["initial_velocity_kms"] == "7.4"
        assert float(d["initial_altitude_km"]) == pytest.approx(90.0)
        assert d["epoch_utc"] == "2024-08-01T12:09:10"
        assert d["max_temperature_K"] == "850.0" and d["final_mass_source"] == "sesam_log_event_end"
        assert d["outcome"] == "demised" and d["final_radius_mm"] == "0.633"
        f = rows[2]
        assert f["status"] == "failed" and f["max_temperature_K"] == "" and f["outcome"] == "" and f["run_name"] == failed.run_name
        assert rows[3]["status"] == "pending" and rows[3]["final_mass_kg"] == ""
        assert not os.path.exists(path + ".tmp")

    def test_stale_json_ignored_for_non_done_status(self, parent, tmp_path):
        """Verify status gate prevents stale JSON results from populating stats for non-done points."""
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4])
        point = points[0]
        point.status = "failed"
        # Write a stale JSON file (leftover from previous attempt) with results
        write_run_json(runs, point.run_name, "ok", results={
            "max_temperature_K": 850.0, "final_mass_kg": 3e-06, "final_mass_source": "sesam_log_event_end",
            "mass_loss_fraction": 0.98, "time_at_melting_temperature_s": 6.774, "final_velocity_kms": 2.998,
            "final_radius_mm": 0.633, "final_altitude_km": 76.548, "end_of_life_reason": "uncritical",
            "outcome": "demised"})
        # Summary row should ignore the stale JSON because status is "failed", not "done"
        row = sw.summary_row(point, runs)
        assert row["status"] == "failed"
        assert row["run_name"] == point.run_name  # run_name is always populated
        for col in sw._RESULT_COLUMNS:
            assert row[col] == "", f"Expected {col} to be blank for non-done status, got {row[col]!r}"

    def test_pending_with_stale_json_also_ignores_results(self, parent, tmp_path):
        """Verify the gate also works for pending status."""
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [5.0])
        point = points[0]
        point.status = "pending"
        # Write a stale JSON file with results
        write_run_json(runs, point.run_name, "ok", results={
            "max_temperature_K": 1000.0, "final_mass_kg": 1e-06, "final_mass_source": "sesam",
            "mass_loss_fraction": 0.95, "time_at_melting_temperature_s": 10.0, "final_velocity_kms": 3.0,
            "final_radius_mm": 1.0, "final_altitude_km": 50.0, "end_of_life_reason": "critical",
            "outcome": "fragmented"})
        row = sw.summary_row(point, runs)
        assert row["status"] == "pending"
        for col in sw._RESULT_COLUMNS:
            assert row[col] == "", f"Expected {col} to be blank for pending status, got {row[col]!r}"


import subprocess
import sys


class TestCommand:
    def test_build_command(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        cmd = sw.build_command(point, "/out/runs", "/out/raw", 600)
        assert cmd[0] == sys.executable and cmd[1] == sw.SPHERE_REENTRY
        assert cmd[2:] == [
            "--velocity", "7.400000", "--altitude", "90.000000", "--temperature", "300.000000",
            "--diameter", "5.000000", "--flight-path-angle", "-0.750000", "--heading", sw.format_arg(point.heading_deg),
            "--lat", "-5.000000", "--lon", "149.000000", "--epoch", "2024-08-01T12:09:10",
            "--outdir", "/out/runs", "--raw-dir", "/out/raw", "--timeout", "600", "--quiet"]
        assert sw.build_command(point, "/o", "/r", 5, python="/usr/bin/python3")[0] == "/usr/bin/python3"


class TestRunPoint:
    def fake_subprocess(self, monkeypatch, runs_dir, returncode=0, stderr="", json_status="ok", raise_timeout=False):
        calls = []

        def fake_run(cmd, capture_output, text, timeout):
            calls.append(cmd)
            if raise_timeout:
                raise subprocess.TimeoutExpired(cmd, timeout)
            if json_status is not None:                       # imitate sphere_reentry writing its JSON
                write_run_json(runs_dir, self.current_name, json_status, results={})
            return subprocess.CompletedProcess(cmd, returncode, stdout="", stderr=stderr)

        monkeypatch.setattr(sw.subprocess, "run", fake_run)
        return calls

    def test_success(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        calls = self.fake_subprocess(monkeypatch, str(tmp_path / "runs"))
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "done" and point.returncode == 0 and point.wall_time_s >= 0.0
        assert point.stderr_tail is None
        assert calls[0][:2] == [sys.executable, sw.SPHERE_REENTRY]

    def test_nonzero_return_code_is_failed_with_stderr_tail(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        stderr = "\n".join("line %d" % i for i in range(30))
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), returncode=1, stderr=stderr, json_status="error")
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed" and point.returncode == 1
        assert point.stderr_tail.splitlines() == ["line %d" % i for i in range(10, 30)]

    def test_timeout_is_failed(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), raise_timeout=True, json_status=None)
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed" and point.returncode == -1 and "timed out" in point.stderr_tail

    def test_zero_return_code_but_json_not_ok_is_failed(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), returncode=0, json_status="error")
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed"

    def test_zero_return_code_without_json_is_failed(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), returncode=0, json_status=None)
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed"


class TestBatchesAndPrompts:
    def test_split_batches(self):
        pts = list(range(7))
        assert sw.split_batches(pts, None) == [pts]
        assert sw.split_batches(pts, 3) == [[0, 1, 2], [3, 4, 5], [6]]
        assert sw.split_batches([], None) == [] and sw.split_batches([], 3) == []

    def test_format_duration_and_default_cores(self):
        assert sw.format_duration(3661) == "1:01:01" and sw.format_duration(0.4) == "0:00:00"
        assert 1 <= sw.default_cores() <= (os.cpu_count() or 1)

    def test_confirm_batch(self, parent, caplog):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0])
        answers = iter(["n"])
        assert sw.confirm_batch(1, 3, batch, 0.3, 7, False, ask=lambda prompt: next(answers)) is False
        answers = iter([""])
        assert sw.confirm_batch(1, 3, batch, 0.3, 7, False, ask=lambda prompt: next(answers)) is True
        answers = iter(["maybe", "y"])
        assert sw.confirm_batch(2, 3, batch, 0.3, 7, False, ask=lambda prompt: next(answers)) is True

        def never(prompt):
            raise AssertionError("must not prompt with --yes")
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.confirm_batch(3, 3, batch, 0.3, 7, True, ask=never) is True
        assert "batch 3/3: 2 runs" in caplog.text and batch[0].run_name in caplog.text

    def test_ask_cores(self):
        assert sw.ask_cores(7, preset=3) == 3
        with pytest.raises(ValueError):
            sw.ask_cores(7, preset=0)
        answers = iter(["abc", "99999", "2"])
        assert sw.ask_cores(7, ask=lambda prompt: next(answers)) == 2
        answers = iter([""])
        assert sw.ask_cores(5, ask=lambda prompt: next(answers)) == 5


class FakeProgress:
    def __init__(self):
        self.updates = 0
        self.postfixes = []

    def update(self, n):
        self.updates += n

    def set_postfix_str(self, s, refresh=True):
        self.postfixes.append(s)

    def close(self):
        pass


class TestRunBatch:
    def test_counts_and_progress(self, parent):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0])

        def runner(point, runs_dir, raw_dir, timeout):
            point.status = "failed" if point.velocity_kms == 5.0 else "done"
            point.returncode = 1 if point.status == "failed" else 0
            return point

        progress = FakeProgress()
        n_ok, n_failed, wall = sw.run_batch(batch, 2, "/r", "/w", 600, progress=progress, runner=runner)
        assert (n_ok, n_failed) == (2, 1) and wall >= 0.0
        assert progress.updates == 3 and len(progress.postfixes) == 3
        assert "d=5mm T=300K v=7.40000km/s" in progress.postfixes

    def test_default_runner_is_looked_up_at_call_time(self, parent, monkeypatch):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4])
        seen = []

        def fake_run_point(point, runs_dir, raw_dir, timeout, material_args=()):
            seen.append(point.run_name)
            point.status = "done"
            return point

        monkeypatch.setattr(sw, "run_point", fake_run_point)
        assert sw.run_batch(batch, 1, "/r", "/w", 600)[:2] == (1, 0)
        assert seen == [batch[0].run_name]

    def test_keyboard_interrupt_propagates_after_shutdown(self, parent):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0])

        def runner(point, runs_dir, raw_dir, timeout):
            raise KeyboardInterrupt

        with pytest.raises(KeyboardInterrupt):
            sw.run_batch(batch, 1, "/r", "/w", 600, runner=runner)

    def test_make_progress_is_a_tqdm_bar(self):
        bar = sw.make_progress(total=10, initial=3)
        assert bar.n == 3 and bar.total == 10
        bar.close()


def fake_run_point_factory(fail_velocities=()):
    """A run_point stand-in that writes an ok JSON (or an error JSON) like sphere_reentry would."""
    def fake_run_point(point, runs_dir, raw_dir, timeout, material_args=()):
        if point.velocity_kms in fail_velocities:
            write_run_json(runs_dir, point.run_name, "error")
            point.status, point.returncode, point.stderr_tail = "failed", 1, "ERROR (error): boom"
        else:
            write_run_json(runs_dir, point.run_name, "ok", results={
                "max_temperature_K": 850.0, "final_mass_kg": 0.1, "final_mass_source": "history_file",
                "mass_loss_fraction": 0.5, "time_at_melting_temperature_s": 3.0, "final_velocity_kms": 1.0,
                "final_radius_mm": 20.0, "final_altitude_km": 0.0, "end_of_life_reason": "ground impact",
                "outcome": "survived"})
            point.status, point.returncode = "done", 0
        point.wall_time_s = 0.01
        return point
    return fake_run_point


SUBSET = ["--diameters", "5", "--temperatures", "300", "--velocities", "7.4,5.0,9.0"]


class TestMain:
    def test_dry_run_logs_v_min_and_writes_nothing(self, tmp_path, caplog):
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--dry-run"] + SUBSET)
        assert rc == 0
        assert "v_min = 0.050000 km/s at t = 1300.000 s, altitude 0.000 km (parent row 13)" in caplog.text
        assert "3 points" in caplog.text and "1 skipped" in caplog.text and "2 pending" in caplog.text
        assert "batch 1/1: 2 runs" in caplog.text and "dry run" in caplog.text
        assert not out.exists()

    def test_full_run_with_yes_and_cores(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "2"] + SUBSET)
        assert rc == 0
        manifest = sw.load_manifest(str(out / "sweep_manifest.json"))
        statuses = [p["status"] for p in manifest["points"]]
        assert sorted(statuses) == ["done", "done", "skipped"]
        assert manifest["parent"]["v_min_kms"] == 0.05 and manifest["settings"]["batch_size"] is None
        assert manifest["grid"]["velocities_kms"] == [7.4, 5.0, 9.0]
        with open(out / "sweep_summary.csv", newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 3 and sum(r["outcome"] == "survived" for r in rows) == 2
        log_text = (out / "sweep.log").read_text()
        assert "v_min = 0.050000 km/s" in log_text and "batch 1/1: 2 runs on 2 cores" in log_text
        assert "sweep finished: 2 ok, 0 failed, 0 pending" in caplog.text

    def test_second_invocation_has_nothing_to_do(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        assert "nothing to do" in caplog.text

    def test_failed_runs_are_retried_only_with_retry_failed(self, tmp_path, monkeypatch, caplog):
        out = tmp_path / "out"
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory(fail_velocities=(5.0,)))
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        assert "1 ok, 1 failed" in caplog.text and "FAILED" in caplog.text
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        caplog.clear()
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        assert "nothing to do" in caplog.text
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1", "--retry-failed"] + SUBSET) == 0
        manifest = sw.load_manifest(str(out / "sweep_manifest.json"))
        assert sorted(p["status"] for p in manifest["points"]) == ["done", "done", "skipped"]

    def test_batches_prompt_for_confirmation_and_cores(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        answers = iter(["", "3", "y", ""])          # batch 1: confirm, 3 cores; batch 2: confirm, default cores
        monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--batch-size", "1"] + SUBSET)
        assert rc == 0
        assert "batch 1/2: 1 runs on 3 cores" in caplog.text
        assert "batch 2/2: 1 runs on {} cores".format(sw.default_cores()) in caplog.text
        assert sw.load_manifest(str(out / "sweep_manifest.json"))["settings"]["batch_size"] == 1

    def test_user_declines_first_batch(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        monkeypatch.setattr("builtins.input", lambda prompt="": "n")
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out)] + SUBSET)
        assert rc == 0 and "stopped before batch 1" in caplog.text
        statuses = [p["status"] for p in sw.load_manifest(str(out / "sweep_manifest.json"))["points"]]
        assert sorted(statuses) == ["pending", "pending", "skipped"]

    def test_yes_without_cores_still_asks_for_cores(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        asked = []
        monkeypatch.setattr("builtins.input", lambda prompt="": asked.append(prompt) or "2")
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes"] + SUBSET) == 0
        assert asked == ["Cores for this batch [{}]: ".format(sw.default_cores())]
        assert "on 2 cores" in caplog.text

    def test_keyboard_interrupt_saves_state_and_returns_130(self, tmp_path, monkeypatch):
        def interrupt(*args, **kwargs):
            raise KeyboardInterrupt
        monkeypatch.setattr(sw, "run_batch", interrupt)
        out = tmp_path / "out"
        rc = sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET)
        assert rc == 130
        assert (out / "sweep_manifest.json").is_file() and (out / "sweep_summary.csv").is_file()

    @pytest.mark.parametrize("extra", [["--batch-size", "0"], ["--limit", "0"], ["--cores", "0"], ["--cores", "100000"]])
    def test_bad_arguments_exit_2(self, tmp_path, extra):
        with pytest.raises(SystemExit) as exc:
            sw.main(["--parent", MINI, "--outdir", str(tmp_path / "out")] + extra)
        assert exc.value.code == 2

    def test_missing_parent_returns_2(self, tmp_path, caplog):
        with caplog.at_level("ERROR", logger="sweep"):
            rc = sw.main(["--parent", str(tmp_path / "nope.json"), "--outdir", str(tmp_path / "out"), "--dry-run"])
        assert rc == 2 and "cannot load parent" in caplog.text

    def test_limit(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1", "--limit", "2"] + SUBSET) == 0
        assert len(sw.load_manifest(str(out / "sweep_manifest.json"))["points"]) == 2


# =============================================================================
# Materials (--material / --material-file, one per sweep invocation)
# =============================================================================

TITANIUM = sr.Material("drama-TiAl6v4", 4417.0, 1905.0)
CUSTOM_MATERIAL = {
    "name": "user-AlLi2195", "materialType": "metal", "catalycity": 1.0, "density": 2700.0,
    "specificHeatCapacity": [[293.0, 900.0], [823.0, 1140.0]], "meltingHeat": 390000.0,
    "meltingTemperature": 823.0, "emissivity": [[50.0, 0.3]],
    "heatConductivity": [[293.0, 130.0], [823.0, 129.5]], "oxideActivationTemperature": 0.0,
    "oxideEmissivity": [[50.0, 0.3]], "oxideHeatOfFormation": 0.0, "oxideReactionProbability": 0.0,
}


def write_material_file(tmp_path, data=CUSTOM_MATERIAL):
    path = tmp_path / "material.json"
    path.write_text(json.dumps(data))
    return str(path)


class TestMaterialArgs:
    def test_default_material_forwards_nothing(self):
        assert sw.material_cli_args("drama-AA7075", None) == []

    def test_builtin_material_is_forwarded_by_name(self):
        assert sw.material_cli_args("drama-TiAl6v4", None) == ["--material", "drama-TiAl6v4"]

    def test_material_file_is_forwarded_as_an_absolute_path(self, tmp_path):
        path = write_material_file(tmp_path)
        assert sw.material_cli_args("drama-AA7075", path) == ["--material-file", os.path.abspath(path)]


class TestMaterialInPoints:
    def test_default_run_names_are_unchanged(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        assert point.run_name == "sphere_d005.00mm_T0300.0K_v07.40000kms_h090.000km"

    def test_non_default_material_suffixes_every_run_name(self, parent):
        points = sw.build_points(parent, [5.0, 10.0], [300.0], [7.4, 5.0], material=TITANIUM)
        assert all(p.run_name.endswith("_mdrama-TiAl6v4") for p in points)
        assert sw.sphere_run_for(points[0], material=TITANIUM).material == TITANIUM
        assert sr.run_name(sw.sphere_run_for(points[0], material=TITANIUM)) == points[0].run_name

    def test_skipped_points_still_have_no_run_name(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [9.0], material=TITANIUM)[0]
        assert point.status == "skipped" and point.run_name is None


class TestMaterialInCommand:
    def test_build_command_forwards_material_args(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        cmd = sw.build_command(point, "/out/runs", "/out/raw", 600, material_args=["--material", "drama-TiAl6v4"])
        assert cmd[-3:] == ["--material", "drama-TiAl6v4", "--quiet"] or cmd[-2:] == ["--material", "drama-TiAl6v4"]
        assert "--material" in cmd and cmd[cmd.index("--material") + 1] == "drama-TiAl6v4"
        assert sw.build_command(point, "/out/runs", "/out/raw", 600) == \
            sw.build_command(point, "/out/runs", "/out/raw", 600, material_args=[])

    def test_run_point_forwards_material_args_to_the_subprocess(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4], material=TITANIUM)[0]
        seen = {}

        def fake_run(cmd, capture_output, text, timeout):
            seen["cmd"] = cmd
            write_run_json(str(tmp_path / "runs"), point.run_name, "ok", results={})
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        monkeypatch.setattr(sw.subprocess, "run", fake_run)
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600,
                     material_args=["--material", "drama-TiAl6v4"])
        assert point.status == "done"
        assert seen["cmd"][seen["cmd"].index("--material") + 1] == "drama-TiAl6v4"


class TestMaterialInMain:
    def test_builtin_material_flows_into_names_subprocess_args_and_manifest(self, tmp_path, monkeypatch, caplog):
        received = []

        def fake_run_point(point, runs_dir, raw_dir, timeout, material_args=()):
            received.append(list(material_args))
            write_run_json(runs_dir, point.run_name, "ok", results={"outcome": "survived"})
            point.status, point.returncode, point.wall_time_s = "done", 0, 0.01
            return point

        monkeypatch.setattr(sw, "run_point", fake_run_point)
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1",
                          "--material", "drama-TiAl6v4"] + SUBSET)
        assert rc == 0
        assert received and all(a == ["--material", "drama-TiAl6v4"] for a in received)
        manifest = sw.load_manifest(str(out / "sweep_manifest.json"))
        assert manifest["settings"]["material"] == {
            "name": "drama-TiAl6v4", "source": "builtin", "density_kgm3": 4417.0,
            "melting_temperature_K": 1905.0, "file": None, "file_sha256": None}
        done = [p for p in manifest["points"] if p["status"] == "done"]
        assert done and all(p["run_name"].endswith("_mdrama-TiAl6v4") for p in done)
        assert (out / "runs" / (done[0]["run_name"] + ".json")).is_file()
        assert "material drama-TiAl6v4 (builtin" in caplog.text

    def test_material_file_is_recorded_with_its_hash(self, tmp_path, monkeypatch):
        path = write_material_file(tmp_path)
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1",
                        "--material-file", path] + SUBSET) == 0
        material = sw.load_manifest(str(out / "sweep_manifest.json"))["settings"]["material"]
        assert material["name"] == "user-AlLi2195" and material["source"] == "custom_file"
        assert material["density_kgm3"] == 2700.0 and material["melting_temperature_K"] == 823.0
        assert material["file"] == os.path.abspath(path)
        assert material["file_sha256"] == sw.sha256_of_file(path)

    def test_default_material_is_recorded_too(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        material = sw.load_manifest(str(out / "sweep_manifest.json"))["settings"]["material"]
        assert material["name"] == "drama-AA7075" and material["source"] == "builtin"

    def test_unknown_material_returns_2_before_anything_runs(self, tmp_path, monkeypatch, caplog):
        called = []
        monkeypatch.setattr(sw, "run_point", lambda *a, **k: called.append(1))
        monkeypatch.setattr(sr, "MATERIAL_DB_PATH", str(tmp_path / "nope.xml"))   # no database at all
        out = tmp_path / "out"
        with caplog.at_level("ERROR", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1",
                          "--material", "drama-TiAl6v4"] + SUBSET)
        assert rc == 2 and not called
        assert "material" in caplog.text.lower()
        assert not (out / "sweep_manifest.json").exists()

    def test_material_flags_are_mutually_exclusive(self, tmp_path, capsys):
        with pytest.raises(SystemExit) as exc:
            sw.main(["--parent", MINI, "--outdir", str(tmp_path / "out"), "--material", "drama-TiAl6v4",
                     "--material-file", write_material_file(tmp_path), "--dry-run"])
        assert exc.value.code == 2
        assert "not allowed with" in capsys.readouterr().err

    def test_dry_run_logs_the_material(self, tmp_path, caplog):
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(tmp_path / "out"), "--dry-run",
                            "--material-file", write_material_file(tmp_path)] + SUBSET) == 0
        assert "material user-AlLi2195 (custom_file" in caplog.text
        assert "_muser-AlLi2195" in caplog.text
