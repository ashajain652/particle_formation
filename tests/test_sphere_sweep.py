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
