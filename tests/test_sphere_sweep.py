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
