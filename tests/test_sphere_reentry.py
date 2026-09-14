"""Unit tests for sphere_reentry.py (no DRAMA needed)."""
import json
import math
from datetime import datetime

import pytest

import sphere_reentry as sr
from helpers import fixture_dir, fixture_file


def make_run(**kw):
    base = dict(velocity_kms=7.5, altitude_km=77.5, temperature_K=300.0, diameter_mm=50.0)
    base.update(kw)
    return sr.SphereRun(**base)


class TestConstants:
    def test_material_and_run_constants(self):
        assert sr.RHO_AA7075 == 2813.0
        assert sr.T_MELT_AA7075 == 850.0
        assert sr.MELT_TOLERANCE_K == 0.5
        assert sr.ENERGY_THRESHOLD_J == 1e-9
        assert sr.MATERIAL_NAME == "drama-AA7075"
        assert sr.PARENT_EPOCH == datetime(2024, 8, 1, 12, 0, 0)
        assert sr.DEFAULT_OUTDIR.endswith("sphere_sweep_output/runs")
        assert sr.DEFAULT_RAW_DIR.endswith("sphere_sweep_output/raw")
        assert sr.DEFAULT_FAP_DAY.endswith("data/fap_day.dat")


class TestGeometry:
    def test_mass_of_5mm_sphere(self):
        assert sr.sphere_mass_kg(5.0) == pytest.approx(1.8411e-4, rel=1e-4)

    def test_mass_of_100mm_sphere(self):
        assert sr.sphere_mass_kg(100.0) == pytest.approx(1.4729, rel=1e-4)

    def test_cross_section_of_50mm_sphere(self):
        assert sr.sphere_cross_section_m2(50.0) == pytest.approx(1.9635e-3, rel=1e-4)

    def test_run_properties(self):
        run = make_run()
        assert run.radius_m == 0.025
        assert run.mass_kg == pytest.approx(0.18411, rel=1e-4)
        assert run.cross_section_m2 == pytest.approx(1.9635e-3, rel=1e-4)
        assert run.epoch == datetime(2024, 8, 1, 12, 0, 0)
        assert run.flight_path_deg == 0.0 and run.heading_deg == 0.0
        assert run.lat_deg == 0.0 and run.lon_deg == 0.0


class TestRunName:
    def test_exact_name(self):
        assert sr.run_name(make_run()) == "sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km"

    def test_zero_padding_of_small_values(self):
        run = make_run(diameter_mm=5.0, temperature_K=750.0, velocity_kms=0.028, altitude_km=0.012)
        assert sr.run_name(run) == "sphere_d005.00mm_T0750.0K_v00.02800kms_h000.012km"

    def test_lexicographic_order_matches_parameter_order(self):
        names = [sr.run_name(make_run(diameter_mm=d, velocity_kms=v))
                 for d in (5.0, 10.0, 100.0) for v in (0.5, 7.5)]
        assert names == sorted(names)

    def test_adjacent_grid_velocities_are_distinct(self):
        step = (7.5 - 0.028) / 99
        assert sr.run_name(make_run(velocity_kms=7.5)) != sr.run_name(make_run(velocity_kms=7.5 - step))


def test_read_lines_drops_blank_lines_and_newlines(tmp_path):
    p = tmp_path / "f.dat"
    p.write_text("# header\n\n01/08/2024 170 170 100 8 3 3 3 3 3 3 3 3\n   \n")
    assert sr.read_lines(str(p)) == ["# header", "01/08/2024 170 170 100 8 3 3 3 3 3 3 3 3"]


class TestBuildConfig:
    def test_geodetic_elements_in_order(self):
        run = make_run(flight_path_deg=-0.959, heading_deg=347.168, lat_deg=29.546, lon_deg=-82.134)
        cfg = sr.build_config(run)
        assert cfg["coordinateSystem"] == "geodetic"
        assert (cfg["element1"], cfg["element2"], cfg["element3"]) == (77.5, 29.546, -82.134)
        assert (cfg["element4"], cfg["element5"], cfg["element6"]) == (7.5, -0.959, 347.168)

    def test_heading_is_normalised_to_0_360(self):
        assert sr.build_config(make_run(heading_deg=-10.0))["element6"] == pytest.approx(350.0)
        assert sr.build_config(make_run(heading_deg=370.0))["element6"] == pytest.approx(10.0)

    def test_object_definition(self):
        cfg = sr.build_config(make_run())
        assert isinstance(cfg["objects"], list) and len(cfg["objects"]) == 1
        obj = cfg["objects"][0]
        assert obj["primitive"] == {"sphere": {"radius": 0.025}}
        assert obj["solid"] is True
        assert obj["material"] == "drama-AA7075"
        assert obj["mass"] == pytest.approx(sr.sphere_mass_kg(50.0))
        assert obj["attitude"] == "tumbling"
        assert obj["quantity"] == 1
        assert obj["name"] == sr.OBJECT_NAME and obj["uniqueID"] == sr.OBJECT_UUID
        assert "wallThickness" not in obj
        assert "materialList" not in cfg

    def test_temperatures_and_epoch(self):
        epoch = datetime(2024, 8, 1, 12, 53, 7)
        cfg = sr.build_config(make_run(temperature_K=450.0, epoch=epoch))
        assert cfg["objects"][0]["temperature"] == 450.0
        assert cfg["globalSpacecraftTemperature"] == 450.0
        assert cfg["beginDate"] == epoch and cfg["initialDate"] == epoch
        assert isinstance(cfg["beginDate"], datetime)

    def test_run_control_settings_mirror_parent_except_energy_threshold(self):
        cfg = sr.build_config(make_run())
        assert cfg["runMode"] == "reentry-only"
        assert cfg["monteCarlo"] is False
        assert cfg["energyThreshold"] == 1e-9
        assert cfg["assumedCrossSection"] == pytest.approx(sr.sphere_cross_section_m2(50.0))
        assert cfg["comment1"] == sr.run_name(make_run())
        expected = {"dragCoefficient": 2.2, "reflectivityCoefficient": 1.3, "attitude": "tumbling",
                    "fragmentsAttitudeAfterBreakup": "inherited", "densityScalingFactor": 1.0,
                    "dynamicEnvironment": True, "useWind": True, "solarActivityFromFile": True,
                    "useEnvironmentCSV": False, "ap": 8, "f107a": 170, "voxelatorMode": 1,
                    "plotVisibilityMaps": False, "plotObjectTrajectories": False,
                    "propagationWithOscar": False, "runID": "SPHERE"}
        for key, value in expected.items():
            assert cfg[key] == value, key

    def test_config_is_json_serialisable_with_default_str(self):
        json.dumps(sr.build_config(make_run()), default=str)
