"""Unit tests for sphere_reentry.py (no DRAMA needed)."""
import json
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


class TestReadSaraTable:
    def test_aero_rows_T1(self):
        rows = sr.read_sara_table(fixture_file("T1_demised_50mm", "_AeroThermalHistory.txt"), sr.AERO_COLUMNS)
        assert len(rows) == 373
        assert list(rows[0].keys()) == sr.AERO_COLUMNS
        first, last = rows[0], rows[-1]
        assert (first["time"], first["altitude"], first["temp"], first["mass"], first["thick"]) == (0.0, 101.247, 300.0, 0.184, 25.0)
        assert (last["time"], last["altitude"], last["temp"], last["mass"], last["thick"]) == (370.55, 80.14, 850.0, 0.0, 0.0)

    def test_traj_rows_T1(self):
        rows = sr.read_sara_table(fixture_file("T1_demised_50mm", "_Trajectory.txt"), sr.TRAJ_COLUMNS)
        assert len(rows) == 373
        first = rows[0]
        assert (first["velocity"], first["path"], first["heading"], first["density"]) == (7.907, -0.18985, 347.95995, 4.493e-07)
        assert (first["lat"], first["lon"]) == (-3.749, -74.64)
        assert rows[-1]["velocity"] == 6.467

    def test_comment_and_malformed_lines_are_skipped(self, tmp_path):
        p = tmp_path / "t.txt"
        p.write_text("# header\n#  Time [s]  x  y\n1.0 2.0 3.0\n4.0 5.0\nabc 1.0 2.0\n\n7.0 8.0 9.0\n")
        rows = sr.read_sara_table(str(p), ["t", "x", "y"])
        assert rows == [{"t": 1.0, "x": 2.0, "y": 3.0}, {"t": 7.0, "x": 8.0, "y": 9.0}]

    def test_empty_file(self, tmp_path):
        p = tmp_path / "e.txt"
        p.write_text("")
        assert sr.read_sara_table(str(p), ["t"]) == []


class TestParseSesamVersion:
    def test_version_from_history_header(self):
        assert sr.parse_sesam_version(fixture_file("T1_demised_50mm", "_AeroThermalHistory.txt")) == "2.3.0"

    def test_no_header_or_missing_file(self, tmp_path):
        p = tmp_path / "x.txt"
        p.write_text("1 2 3\n")
        assert sr.parse_sesam_version(str(p)) is None
        assert sr.parse_sesam_version(str(tmp_path / "nope.txt")) is None


class TestParseImpactingFragments:
    def test_survivor_T3(self):
        frag = sr.parse_impacting_fragments(fixture_file("T3_survivor_50mm", "ImpactingFragments.xml"))
        assert frag["mass_kg"] == 0.18411041946975187
        assert frag["velocity_kms"] == pytest.approx(0.058229426237796787)
        assert frag["lat_deg"] == pytest.approx(36.041459697933099)
        assert frag["lon_deg"] == pytest.approx(-83.913412639736748)
        assert frag["epoch"] == "2024-08-01T13:00:12.080"

    def test_no_fragment_T1(self):
        assert sr.parse_impacting_fragments(fixture_file("T1_demised_50mm", "ImpactingFragments.xml")) is None

    def test_missing_or_malformed_file(self, tmp_path):
        assert sr.parse_impacting_fragments(str(tmp_path / "nope.xml")) is None
        assert sr.parse_impacting_fragments(None) is None
        bad = tmp_path / "bad.xml"
        bad.write_text("<fragments><fragment>")
        assert sr.parse_impacting_fragments(str(bad)) is None


class TestParseSesamLog:
    @pytest.mark.parametrize("case, reason, end_mass, start_mass", [
        ("T1_demised_50mm", "uncritical", 0.0, 0.18411),
        ("T3_survivor_50mm", "ground impact", 0.18411, 0.18411),
        ("T5_5mm_750K", "uncritical", 3e-06, 0.000184),
        ("E1_ballooning_5mm", "ballooning", 3e-06, 0.000184),
    ])
    def test_fixture_logs(self, case, reason, end_mass, start_mass):
        info = sr.parse_sesam_log(fixture_file(case, "sesam.log"))
        assert info["end_of_life_reason"] == reason
        assert info["event_end_mass_kg"] == pytest.approx(end_mass, abs=1e-12)
        assert info["event_start_mass_kg"] == pytest.approx(start_mass, abs=1e-12)

    def test_missing_file(self, tmp_path):
        assert sr.parse_sesam_log(str(tmp_path / "nope.log")) == {
            "end_of_life_reason": "unknown", "event_end_mass_kg": None, "event_start_mass_kg": None}
        assert sr.parse_sesam_log(None)["end_of_life_reason"] == "unknown"


class TestFindOutputFiles:
    def test_finds_files_recursively(self, tmp_path):
        import shutil
        dest = tmp_path / "raw" / "run_0" / "reentry"
        shutil.copytree(fixture_dir("T3_survivor_50mm"), dest)
        files = sr.find_output_files(str(tmp_path / "raw"))
        assert files["aero"].endswith("_AeroThermalHistory.txt")
        assert files["traj"].endswith("_Trajectory.txt")
        assert files["fragments"].endswith("ImpactingFragments.xml")
        assert files["log"].endswith("sesam.log")

    def test_missing_files_are_none(self, tmp_path):
        assert sr.find_output_files(str(tmp_path)) == {"aero": None, "traj": None, "fragments": None, "log": None}


def aero_row(t, temp=300.0, mass=0.184, alt=100.0, thick=25.0):
    return {"time": t, "altitude": alt, "temp": temp, "mass": mass, "thick": thick,
            "convectiveHeat": 1.0, "radiativeHeat": 0.0, "oxidationHeat": 0.0,
            "radCooling": -1.0, "integratedHeat": 10.0, "visibilityFactor": 1.0}


def traj_row(t, v=7.5, alt=100.0, lat=1.0, lon=2.0, downrange=3.0, path=-1.0, heading=350.0):
    return {"time": t, "altitude": alt, "lat": lat, "lon": lon, "velocity": v,
            "downrange": downrange, "drag": 2.2, "lift": 0.0, "side": 0.0, "knudsen": 0.1,
            "mach": 20.0, "path": path, "heading": heading, "density": 1e-6,
            "dynamicPressure": 100.0, "loadFactor": 0.5}


def load_fixture_rows(case):
    aero = sr.read_sara_table(fixture_file(case, "_AeroThermalHistory.txt"), sr.AERO_COLUMNS)
    traj = sr.read_sara_table(fixture_file(case, "_Trajectory.txt"), sr.TRAJ_COLUMNS)
    return aero, traj


class TestMergeHistories:
    def test_csv_columns_are_the_spec_list(self):
        assert sr.CSV_COLUMNS == [
            "time_s", "altitude_km", "velocity_kms", "temperature_K", "mass_kg", "thick_mm",
            "lat_deg", "lon_deg", "downrange_km", "flight_path_deg", "heading_deg",
            "drag", "lift", "side", "knudsen", "mach", "density_kgm3", "dynamic_pressure_Pa",
            "load_factor_g", "convective_heat_W", "radiative_heat_W", "oxidation_heat_W",
            "rad_cooling_W", "integrated_heat_J", "visibility_factor"]

    def test_identical_time_grids_zip_T1(self):
        aero, traj = load_fixture_rows("T1_demised_50mm")
        rows, warnings = sr.merge_histories(aero, traj)
        assert warnings == []
        assert len(rows) == 373
        assert set(rows[0]) == set(sr.CSV_COLUMNS)
        first = rows[0]
        assert (first["time_s"], first["altitude_km"], first["velocity_kms"]) == (0.0, 101.247, 7.907)
        assert (first["temperature_K"], first["mass_kg"], first["thick_mm"]) == (300.0, 0.184, 25.0)
        assert (first["flight_path_deg"], first["heading_deg"], first["density_kgm3"]) == (-0.18985, 347.95995, 4.493e-07)
        assert rows[-1]["velocity_kms"] == 6.467 and rows[-1]["mass_kg"] == 0.0

    def test_mismatched_grids_outer_join_with_warning(self):
        aero = [aero_row(0.0), aero_row(1.0, temp=400.0), aero_row(2.0)]
        traj = [traj_row(0.0), traj_row(2.0, v=7.0), traj_row(3.0, v=6.0)]
        rows, warnings = sr.merge_histories(aero, traj)
        assert warnings == ["time grids differ (aero 3 rows, traj 3 rows)"]
        assert [r["time_s"] for r in rows] == [0.0, 1.0, 2.0, 3.0]
        assert rows[1]["temperature_K"] == 400.0 and rows[1]["velocity_kms"] is None
        assert rows[3]["velocity_kms"] == 6.0 and rows[3]["temperature_K"] is None
        assert rows[1]["altitude_km"] == 100.0  # falls back to the aero altitude


class TestWriteCsv:
    def test_header_and_first_row_T5(self, tmp_path):
        aero, traj = load_fixture_rows("T5_5mm_750K")
        rows, _ = sr.merge_histories(aero, traj)
        path = tmp_path / "out.csv"
        sr.write_csv(str(path), rows)
        lines = path.read_text().splitlines()
        assert lines[0] == ",".join(sr.CSV_COLUMNS)
        assert lines[1].startswith("0.0,77.5,7.5,750.0,0.0,2.5,29.546,-82.134,")
        assert ",2.727e-05," in lines[1]
        assert len(lines) == 1 + 16

    def test_none_is_written_as_empty_field(self, tmp_path):
        rows, _ = sr.merge_histories([aero_row(0.0)], [traj_row(1.0)])
        path = tmp_path / "out.csv"
        sr.write_csv(str(path), rows)
        line = path.read_text().splitlines()[1]
        assert line.startswith("0.0,100.0,,300.0,0.184,25.0,,,")


FIXTURE_RUNS = {
    "T1_demised_50mm": dict(velocity_kms=7.907, altitude_km=101.247, temperature_K=300.0, diameter_mm=50.0),
    "T3_survivor_50mm": dict(velocity_kms=0.5, altitude_km=39.94, temperature_K=300.0, diameter_mm=50.0),
    "T5_5mm_750K": dict(velocity_kms=7.5, altitude_km=77.5, temperature_K=750.0, diameter_mm=5.0),
    "E1_ballooning_5mm": dict(velocity_kms=7.5, altitude_km=77.5, temperature_K=750.0, diameter_mm=5.0),
}


def stats_for(case):
    aero, traj = load_fixture_rows(case)
    rows, _ = sr.merge_histories(aero, traj)
    fragments = sr.parse_impacting_fragments(fixture_file(case, "ImpactingFragments.xml"))
    log_info = sr.parse_sesam_log(fixture_file(case, "sesam.log"))
    return sr.compute_stats(make_run(**FIXTURE_RUNS[case]), rows, fragments, log_info)


def synthetic_stats(temps, dt=1.0, mass=0.184, fragments=None, reason="ground impact", event_end=None):
    aero = [aero_row(i * dt, temp=T, mass=mass) for i, T in enumerate(temps)]
    traj = [traj_row(i * dt) for i in range(len(temps))]
    rows, _ = sr.merge_histories(aero, traj)
    log_info = {"end_of_life_reason": reason, "event_end_mass_kg": event_end, "event_start_mass_kg": None}
    return sr.compute_stats(make_run(), rows, fragments, log_info)


class TestMaxTemperature:
    def test_T1_first_occurrence_of_maximum(self):
        r, _ = stats_for("T1_demised_50mm")
        assert (r["max_temperature_K"], r["time_of_max_temperature_s"], r["altitude_of_max_temperature_km"]) == (850.0, 337.55, 84.46)

    def test_T3_survivor(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert (r["max_temperature_K"], r["time_of_max_temperature_s"], r["altitude_of_max_temperature_km"]) == (304.178, 86.08, 15.424)


class TestMeltDuration:
    def test_T1_sums_the_1s_steps_at_850K(self):
        r, _ = stats_for("T1_demised_50mm")
        assert r["time_at_melting_temperature_s"] == pytest.approx(33.0)
        assert r["n_rows_at_melt"] == 34
        assert (r["first_time_at_melt_s"], r["last_time_at_melt_s"]) == (337.55, 370.55)
        assert (r["altitude_first_melt_km"], r["altitude_last_melt_km"]) == (84.46, 80.14)
        assert (r["melting_temperature_K"], r["melt_tolerance_K"]) == (850.0, 0.5)

    def test_T3_never_melts(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert r["time_at_melting_temperature_s"] == 0.0 and r["n_rows_at_melt"] == 0
        assert r["first_time_at_melt_s"] is None and r["altitude_last_melt_km"] is None

    def test_T5_variable_steps(self):
        r, _ = stats_for("T5_5mm_750K")
        assert r["time_at_melting_temperature_s"] == pytest.approx(6.774)
        assert r["n_rows_at_melt"] == 13

    def test_remelt_sums_both_intervals(self):
        r, _ = synthetic_stats([300.0, 850.0, 850.0, 850.0, 700.0, 850.0, 850.0])
        assert r["time_at_melting_temperature_s"] == pytest.approx(3.0)
        assert (r["first_time_at_melt_s"], r["last_time_at_melt_s"]) == (1.0, 6.0)

    def test_tolerance_boundary(self):
        r, _ = synthetic_stats([849.4, 849.6])
        assert r["time_at_melting_temperature_s"] == 0.0 and r["n_rows_at_melt"] == 1
        r, _ = synthetic_stats([849.6, 850.4])
        assert r["time_at_melting_temperature_s"] == pytest.approx(1.0) and r["n_rows_at_melt"] == 2

    def test_single_row_at_melt_has_zero_duration(self):
        r, _ = synthetic_stats([850.0])
        assert r["time_at_melting_temperature_s"] == 0.0 and r["n_rows_at_melt"] == 1


class TestFinalMassAndRadius:
    def test_T3_from_xml(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert r["final_mass_kg"] == 0.18411041946975187
        assert r["final_mass_source"] == "impacting_fragments_xml"
        assert r["final_radius_mm"] == pytest.approx(25.0, abs=1e-6)
        assert r["mass_loss_fraction"] == pytest.approx(0.0, abs=1e-9)
        assert r["final_thick_mm"] == 25.0
        assert r["initial_mass_kg"] == pytest.approx(0.18411041946975187)

    def test_T1_from_log(self):
        r, _ = stats_for("T1_demised_50mm")
        assert r["final_mass_kg"] == 0.0 and r["final_mass_source"] == "sesam_log_event_end"
        assert r["final_radius_mm"] == 0.0 and r["mass_loss_fraction"] == 1.0
        assert r["final_thick_mm"] == 0.0

    def test_T5_residual_from_log(self):
        r, _ = stats_for("T5_5mm_750K")
        assert r["final_mass_kg"] == pytest.approx(3e-6)
        assert r["final_mass_source"] == "sesam_log_event_end"
        assert r["final_radius_mm"] == pytest.approx(0.633, abs=0.005)

    def test_source_priority(self):
        r, _ = synthetic_stats([300.0], mass=0.123, fragments={"mass_kg": 0.1, "velocity_kms": None}, event_end=0.11)
        assert (r["final_mass_kg"], r["final_mass_source"]) == (0.1, "impacting_fragments_xml")
        r, _ = synthetic_stats([300.0], mass=0.123, fragments=None, event_end=0.11)
        assert (r["final_mass_kg"], r["final_mass_source"]) == (0.11, "sesam_log_event_end")
        r, _ = synthetic_stats([300.0], mass=0.123, fragments=None, event_end=None)
        assert (r["final_mass_kg"], r["final_mass_source"]) == (0.123, "history_file")

    def test_radius_from_mass_round_trips(self):
        assert sr.radius_mm_from_mass(sr.sphere_mass_kg(37.0)) == pytest.approx(18.5)
        assert sr.radius_mm_from_mass(0.0) == 0.0


class TestFinalVelocityAndTrajectoryEnd:
    def test_T3_velocity_from_xml(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert r["final_velocity_kms"] == pytest.approx(0.058229426237796787)
        assert r["final_velocity_source"] == "impacting_fragments_xml"
        assert r["final_altitude_km"] == 0.0 and r["final_time_s"] == 261.08

    def test_T1_velocity_from_trajectory(self):
        r, _ = stats_for("T1_demised_50mm")
        assert (r["final_velocity_kms"], r["final_velocity_source"]) == (6.467, "trajectory_file")
        assert (r["final_time_s"], r["final_altitude_km"]) == (370.55, 80.14)
        assert (r["final_latitude_deg"], r["final_longitude_deg"], r["downrange_km"]) == (21.297, -80.098, 2850.943)
        assert r["n_rows"] == 373


class TestOutcome:
    def test_fixture_outcomes(self):
        assert stats_for("T3_survivor_50mm")[0]["outcome"] == "survived"
        assert stats_for("E1_ballooning_5mm")[0]["outcome"] == "demised"
        assert stats_for("T5_5mm_750K")[0]["outcome"] == "demised"      # uncritical, 1.6 % residual
        assert stats_for("T1_demised_50mm")[0]["outcome"] == "demised"
        assert stats_for("E1_ballooning_5mm")[0]["end_of_life_reason"] == "ballooning"

    def test_classify_outcome_rules(self):
        assert sr.classify_outcome("ground impact", 0.0, 0.1) == "survived"
        assert sr.classify_outcome("ballooning", 0.09, 0.1) == "demised"
        assert sr.classify_outcome("uncritical", 0.184, 0.184) == "other"
        assert sr.classify_outcome("uncritical", 0.004, 0.1) == "demised"
        assert sr.classify_outcome("unknown", None, 0.1) == "other"


class TestWarnings:
    def test_coarse_csv_mass_warning_only_for_small_spheres(self):
        assert any("coarse" in w for w in stats_for("T5_5mm_750K")[1])
        assert not any("coarse" in w for w in stats_for("T1_demised_50mm")[1])

    def test_survivor_without_xml_and_unknown_reason_and_single_row(self):
        _, w = synthetic_stats([300.0, 300.0], fragments=None, reason="ground impact")
        assert any("ImpactingFragments" in x for x in w)
        _, w = synthetic_stats([300.0, 300.0], reason="unknown")
        assert any("end-of-life reason" in x for x in w)
        _, w = synthetic_stats([300.0])
        assert any("single row" in x for x in w)
        _, w = stats_for("T1_demised_50mm")
        assert not any("ImpactingFragments" in x for x in w)


import os
import shutil
import sys
import types


def install_fake_drama(monkeypatch, run_impl, version="stub-4.1.4"):
    drama = types.ModuleType("drama")
    drama.__version__ = version
    sara = types.ModuleType("drama.sara")
    sara.run = run_impl
    drama.sara = sara
    monkeypatch.setitem(sys.modules, "drama", drama)
    monkeypatch.setitem(sys.modules, "drama.sara", sara)


def fake_sara_run(case=None, mode="ok", calls=None):
    """A stand-in for drama.sara.run: copies fixture `case` where pyDRAMA would write its output."""
    def run(config, save_output_dirs, keep_output_files, fap_day_content, fap_mon_content,
            parallel, timeout, log_level, spell_check):
        if calls is not None:
            calls.append(dict(config=config, save_output_dirs=save_output_dirs,
                              keep_output_files=keep_output_files, fap_day_content=fap_day_content,
                              fap_mon_content=fap_mon_content, parallel=parallel, timeout=timeout,
                              log_level=log_level, spell_check=spell_check))
        dest = os.path.join(save_output_dirs, "run_0", "reentry")
        if mode == "raise":
            raise RuntimeError("boom")
        if mode == "error":
            os.makedirs(dest)
            return {"config": config, "errors": [{"status": "error in reentry", "reentry_logfile": "line1\nline2\nfatal"}], "results": []}
        if mode == "timeout":
            os.makedirs(dest)
            return {"config": config, "errors": [{"status": "error: timeout in reentry (600s)", "reentry_logfile": ""}], "results": []}
        if mode == "empty":
            os.makedirs(dest)
            return {"config": config, "errors": [], "results": []}
        if mode == "nofiles":
            os.makedirs(dest)
            return {"config": config, "errors": [], "results": [{"status": "success", "reentry_logfile": "no files"}]}
        shutil.copytree(fixture_dir(case), dest)
        if mode == "headeronly":
            # Strip all data rows from the aero/traj history files, leaving only '#' header lines
            # (simulates SESAM writing header-only history files with zero data rows).
            for fname in os.listdir(dest):
                if fname.endswith(("_AeroThermalHistory.txt", "_Trajectory.txt")):
                    fpath = os.path.join(dest, fname)
                    with open(fpath) as fh:
                        header_lines = [ln for ln in fh if ln.startswith("#")]
                    with open(fpath, "w") as fh:
                        fh.writelines(header_lines)
        return {"config": config, "errors": [], "results": [{"status": "success", "reentry_logfile": "", "config": {"output_dir": dest}}]}
    return run


T3_RUN = dict(velocity_kms=0.5, altitude_km=39.94, temperature_K=300.0, diameter_mm=50.0,
              flight_path_deg=-32.79047, heading_deg=347.12153, lat_deg=35.911, lon_deg=-83.878,
              epoch=datetime(2024, 8, 1, 12, 55, 51))


class TestRunSphereOk:
    def test_writes_csv_and_json_and_deletes_raw_tree(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm"))
        run = make_run(**T3_RUN)
        outdir, raw = tmp_path / "runs", tmp_path / "raw"
        doc = sr.run_sphere(run, str(outdir), str(raw), 600, False, ["# fap day"], ["# fap mon"])
        name = sr.run_name(run)
        assert doc["status"] == "ok" and doc["error"] is None and doc["run_name"] == name
        assert doc["schema_version"] == 1
        assert (outdir / (name + ".csv")).is_file() and (outdir / (name + ".json")).is_file()
        assert not (raw / name).exists() and doc["files"]["raw_dir"] is None
        assert doc["files"]["csv"] == os.path.abspath(str(outdir / (name + ".csv")))
        r = doc["results"]
        assert r["outcome"] == "survived" and r["final_mass_source"] == "impacting_fragments_xml"
        assert doc["inputs"]["initial_mass_kg"] == pytest.approx(0.18411, rel=1e-4)
        assert doc["inputs"]["epoch_utc"] == "2024-08-01T12:55:51"
        assert doc["inputs"]["sesam_settings"]["energyThreshold"] == 1e-9
        assert "objects" not in doc["inputs"]["sesam_settings"]
        assert doc["inputs"]["material"] == "drama-AA7075" and doc["inputs"]["melting_temperature_K"] == 850.0
        assert doc["provenance"]["sesam_version"] == "2.3.0"
        assert doc["provenance"]["pydrama_version"] == "stub-4.1.4"
        assert doc["provenance"]["script_version"] == sr.SCRIPT_VERSION
        assert doc["provenance"]["wall_time_s"] >= 0.0 and doc["provenance"]["created_utc"].endswith("Z")
        on_disk = json.load(open(outdir / (name + ".json")))
        assert on_disk["results"]["final_mass_kg"] == r["final_mass_kg"] and on_disk["status"] == "ok"
        with open(outdir / (name + ".csv")) as fh:
            header, first = fh.readline().strip(), fh.readline().strip().split(",")
        assert header == ",".join(sr.CSV_COLUMNS)
        assert first[:4] == ["0.0", "39.94", "0.5", "300.0"]

    def test_keep_raw(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run("T1_demised_50mm"))
        run = make_run(velocity_kms=7.907, altitude_km=101.247)
        raw = tmp_path / "raw"
        (raw / sr.run_name(run)).mkdir(parents=True)
        (raw / sr.run_name(run) / "stale.txt").write_text("old")
        doc = sr.run_sphere(run, str(tmp_path / "runs"), str(raw), 600, True, [], [])
        assert doc["status"] == "ok"
        assert doc["files"]["raw_dir"] == str(raw / sr.run_name(run))
        assert (raw / sr.run_name(run) / "run_0" / "reentry" / "sesam.log").is_file()
        assert not (raw / sr.run_name(run) / "stale.txt").exists()
        assert doc["results"]["outcome"] == "demised"

    def test_pydrama_call_arguments(self, tmp_path, monkeypatch):
        calls = []
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm", calls=calls))
        sr.run_sphere(make_run(**T3_RUN), str(tmp_path / "runs"), str(tmp_path / "raw"), 123, False, ["d1", "d2"], ["m1"])
        call = calls[0]
        assert isinstance(call["config"], list) and len(call["config"]) == 1
        assert call["config"][0]["element1"] == 39.94 and call["config"][0]["objects"][0]["solid"] is True
        assert call["save_output_dirs"] == str(tmp_path / "raw" / sr.run_name(make_run(**T3_RUN)))
        assert call["keep_output_files"] == "all" and call["parallel"] is False
        assert call["fap_day_content"] == ["d1", "d2"] and call["fap_mon_content"] == ["m1"]
        assert call["timeout"] == 123 and call["log_level"] == "ERROR" and call["spell_check"] is False
        assert os.environ.get("DRAMA_INSTALL_PATH")


class TestRunSphereErrors:
    @pytest.mark.parametrize("mode, status, needle", [
        ("error", "error", "error in reentry"),
        ("timeout", "timeout", "timeout in reentry"),
        ("empty", "error", "no result"),
        ("raise", "error", "pyDRAMA raised"),
        ("nofiles", "error", "history files not found"),
    ])
    def test_failure_modes_write_error_json_and_keep_raw(self, tmp_path, monkeypatch, mode, status, needle):
        install_fake_drama(monkeypatch, fake_sara_run(mode=mode))
        run = make_run(**T3_RUN)
        doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
        assert doc["status"] == status and needle in doc["error"]
        assert doc["results"] is None
        assert doc["files"]["raw_dir"] == str(tmp_path / "raw" / sr.run_name(run))
        on_disk = json.load(open(tmp_path / "runs" / (sr.run_name(run) + ".json")))
        assert on_disk["status"] == status
        assert not (tmp_path / "runs" / (sr.run_name(run) + ".csv")).exists()

    def test_error_includes_log_tail(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run(mode="error"))
        doc = sr.run_sphere(make_run(**T3_RUN), str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
        assert "fatal" in doc["error"]

    def test_header_only_history_is_an_error_not_a_silent_ok(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm", mode="headeronly"))
        run = make_run(**T3_RUN)
        doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
        assert doc["status"] == "error"
        assert "no data rows" in doc["error"]
        assert doc["results"] is None
        assert doc["files"]["raw_dir"] == str(tmp_path / "raw" / sr.run_name(run))
        on_disk = json.load(open(tmp_path / "runs" / (sr.run_name(run) + ".json")))
        assert on_disk["status"] == "error"
        assert not (tmp_path / "runs" / (sr.run_name(run) + ".csv")).exists()

    def test_pydrama_not_importable(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "drama", None)
        with pytest.raises(sr.DramaNotAvailable):
            sr.run_sphere(make_run(**T3_RUN), str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
import subprocess
from helpers import REPO_ROOT

BASE_ARGS = ["--velocity", "7.5", "--altitude", "77.5", "--temperature", "300", "--diameter", "50"]


def fap_files(tmp_path):
    day, mon = tmp_path / "fap_day.dat", tmp_path / "fap_mon.dat"
    day.write_text("# fap day\n01/08/2024 170 170 100 8 3 3 3 3 3 3 3 3\n")
    mon.write_text("# fap mon\n")
    return ["--fap-day", str(day), "--fap-mon", str(mon)]


class TestCli:
    def test_dry_run_prints_config_and_creates_nothing(self, tmp_path, capsys):
        rc = sr.main(BASE_ARGS + ["--dry-run", "--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")])
        assert rc == 0
        out = capsys.readouterr().out
        assert '"coordinateSystem": "geodetic"' in out
        assert "run name: sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km" in out
        assert not (tmp_path / "runs").exists() and not (tmp_path / "raw").exists()

    @pytest.mark.parametrize("argv", [
        ["--velocity", "7.5"],                                   # missing required
        BASE_ARGS[:6] + ["--diameter", "-5"],                    # negative diameter
        BASE_ARGS[:6] + ["--diameter", "50", "--epoch", "yesterday"],
        ["--velocity", "0", "--altitude", "77.5", "--temperature", "300", "--diameter", "50"],
    ])
    def test_bad_arguments_exit_2(self, argv):
        with pytest.raises(SystemExit) as exc:
            sr.main(argv)
        assert exc.value.code == 2

    def test_missing_fap_file_exits_2(self, tmp_path):
        with pytest.raises(SystemExit) as exc:
            sr.main(BASE_ARGS + ["--fap-day", str(tmp_path / "nope.dat"), "--outdir", str(tmp_path)])
        assert exc.value.code == 2

    def test_parse_epoch(self):
        assert sr.parse_epoch("2024-08-01T12:53:07Z") == datetime(2024, 8, 1, 12, 53, 7)
        assert sr.parse_epoch("2024-08-01T12:53:07") == datetime(2024, 8, 1, 12, 53, 7)

    def test_defaults_of_optional_state(self):
        args = sr.build_parser().parse_args(BASE_ARGS)
        assert (args.flight_path_angle, args.heading, args.lat, args.lon) == (0.0, 0.0, 0.0, 0.0)
        assert args.epoch == sr.PARENT_EPOCH and args.timeout == 600
        assert args.outdir == sr.DEFAULT_OUTDIR and args.raw_dir == sr.DEFAULT_RAW_DIR

    def test_ok_run_returns_0_and_prints_summary(self, tmp_path, monkeypatch, capsys):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm"))
        rc = sr.main(["--velocity", "0.5", "--altitude", "39.94", "--temperature", "300", "--diameter", "50",
                      "--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 0
        out = capsys.readouterr().out
        assert "survived" in out and "csv  ->" in out
        assert (tmp_path / "runs" / "sphere_d050.00mm_T0300.0K_v00.50000kms_h039.940km.json").is_file()

    def test_quiet_prints_nothing_on_success(self, tmp_path, monkeypatch, capsys):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm"))
        rc = sr.main(["--velocity", "0.5", "--altitude", "39.94", "--temperature", "300", "--diameter", "50", "--quiet",
                      "--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 0 and capsys.readouterr().out == ""

    def test_failed_run_returns_1(self, tmp_path, monkeypatch, capsys):
        install_fake_drama(monkeypatch, fake_sara_run(mode="error"))
        rc = sr.main(BASE_ARGS + ["--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 1 and "error in reentry" in capsys.readouterr().err

    def test_drama_not_importable_returns_2(self, tmp_path, monkeypatch, capsys):
        monkeypatch.setitem(sys.modules, "drama", None)
        rc = sr.main(BASE_ARGS + ["--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 2 and "pyDRAMA not importable" in capsys.readouterr().err


FAKE_DRAMA_SARA = '''
import os, shutil

def run(config, save_output_dirs, **kwargs):
    dest = os.path.join(save_output_dirs, "run_0", "reentry")
    if os.environ.get("FAKE_DRAMA_MODE") == "error":
        os.makedirs(dest, exist_ok=True)
        return {"config": config, "errors": [{"status": "error in reentry", "reentry_logfile": "boom"}], "results": []}
    shutil.copytree(os.environ["FAKE_DRAMA_FIXTURE"], dest)
    return {"config": config, "errors": [], "results": [{"status": "success", "reentry_logfile": ""}]}
'''


class TestCliSubprocess:
    def make_fake_package(self, tmp_path):
        pkg = tmp_path / "fakepkg" / "drama"
        pkg.mkdir(parents=True, exist_ok=True)
        (pkg / "__init__.py").write_text('__version__ = "fake"\n')
        (pkg / "sara.py").write_text(FAKE_DRAMA_SARA)
        return str(tmp_path / "fakepkg")

    def run_cli(self, tmp_path, extra, mode="ok"):
        env = dict(os.environ, PYTHONPATH=self.make_fake_package(tmp_path),
                   FAKE_DRAMA_FIXTURE=fixture_dir("T3_survivor_50mm"), FAKE_DRAMA_MODE=mode)
        cmd = [sys.executable, os.path.join(REPO_ROOT, "sphere_reentry.py")] + extra
        return subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(tmp_path))

    def test_exit_codes(self, tmp_path):
        common = ["--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path)
        ok = self.run_cli(tmp_path, ["--velocity", "0.5", "--altitude", "39.94", "--temperature", "300", "--diameter", "50"] + common)
        assert ok.returncode == 0, ok.stderr
        assert (tmp_path / "runs" / "sphere_d050.00mm_T0300.0K_v00.50000kms_h039.940km.csv").is_file()
        failed = self.run_cli(tmp_path, BASE_ARGS + common, mode="error")
        assert failed.returncode == 1
        usage = self.run_cli(tmp_path, ["--velocity", "7.5"])
        assert usage.returncode == 2


# =============================================================================
# Materials (--material / --material-file)
# =============================================================================

MINI_MATERIAL_DB = """<?xml version="1.0"?>
<materialList>
  <creationDate>2023-06-28</creationDate>
  <metalMaterial>
    <name>drama-AA7075</name>
    <density>2813.0</density>
    <meltingTemperature>850.0</meltingTemperature>
  </metalMaterial>
  <metalMaterial>
    <name>drama-TiAl6v4</name>
    <density>4417.0</density>
    <meltingTemperature>1905</meltingTemperature>
  </metalMaterial>
  <cfrpMaterial>
    <name>drama-CFRP</name>
    <densityFibre>1800.0</densityFibre>
  </cfrpMaterial>
</materialList>
"""


@pytest.fixture
def mini_db(tmp_path):
    path = tmp_path / "material_database.xml"
    path.write_text(MINI_MATERIAL_DB)
    return str(path)


class TestBuiltinMaterials:
    def test_default_material_is_aa7075_from_the_module_constants(self):
        m = sr.DEFAULT_MATERIAL
        assert (m.name, m.density_kgm3, m.melting_temperature_K) == ("drama-AA7075", 2813.0, 850.0)
        assert m.source == "builtin" and m.definition is None

    def test_list_builtin_materials_is_metals_only(self, mini_db):
        names = sorted(sr.list_builtin_materials(mini_db))
        assert names == ["drama-AA7075", "drama-TiAl6v4"]      # the cfrpMaterial is excluded

    def test_load_builtin_material_reads_density_and_melting_point(self, mini_db):
        m = sr.load_builtin_material("drama-TiAl6v4", mini_db)
        assert m == sr.Material("drama-TiAl6v4", 4417.0, 1905.0)
        assert m.source == "builtin" and m.definition is None

    def test_unknown_name_lists_the_valid_ones(self, mini_db):
        with pytest.raises(sr.MaterialError) as exc:
            sr.load_builtin_material("drama-Unobtainium", mini_db)
        assert "drama-Unobtainium" in str(exc.value)
        assert "drama-AA7075" in str(exc.value) and "drama-TiAl6v4" in str(exc.value)

    def test_missing_database_is_a_material_error(self, tmp_path):
        with pytest.raises(sr.MaterialError, match="material database"):
            sr.load_builtin_material("drama-AA7075", str(tmp_path / "nope.xml"))

    @pytest.mark.skipif(not os.path.isfile(getattr(sr, "MATERIAL_DB_PATH", "")), reason="DRAMA install not present")
    def test_real_database_agrees_with_the_hardcoded_default(self):
        assert sr.load_builtin_material("drama-AA7075") == sr.DEFAULT_MATERIAL
        assert len(sr.list_builtin_materials()) == 21


CUSTOM_MATERIAL = {
    "name": "user-AlLi2195",
    "materialType": "metal",
    "catalycity": 1.0,
    "density": 2700.0,
    "specificHeatCapacity": [[293.0, 900.0], [823.0, 1140.0]],
    "meltingHeat": 390000.0,
    "meltingTemperature": 823.0,
    "emissivity": [[50.0, 0.3]],
    "heatConductivity": [[293.0, 130.0], [823.0, 129.5]],
    "oxideActivationTemperature": 0.0,
    "oxideEmissivity": [[50.0, 0.3]],
    "oxideHeatOfFormation": 0.0,
    "oxideReactionProbability": 0.0,
}


def write_material_file(tmp_path, data, name="material.json"):
    path = tmp_path / name
    path.write_text(json.dumps(data))
    return str(path)


class TestMaterialFile:
    def test_loads_a_complete_definition(self, tmp_path):
        m = sr.load_material_file(write_material_file(tmp_path, CUSTOM_MATERIAL))
        assert (m.name, m.density_kgm3, m.melting_temperature_K) == ("user-AlLi2195", 2700.0, 823.0)
        assert m.source == "custom_file"
        assert m.definition == CUSTOM_MATERIAL

    def test_missing_required_fields_are_named(self, tmp_path):
        data = {k: v for k, v in CUSTOM_MATERIAL.items() if k not in ("meltingHeat", "oxideEmissivity")}
        with pytest.raises(sr.MaterialError) as exc:
            sr.load_material_file(write_material_file(tmp_path, data))
        assert "meltingHeat" in str(exc.value) and "oxideEmissivity" in str(exc.value)

    def test_non_metal_type_is_rejected(self, tmp_path):
        data = dict(CUSTOM_MATERIAL, materialType="cfrp")
        with pytest.raises(sr.MaterialError, match="metal"):
            sr.load_material_file(write_material_file(tmp_path, data))

    def test_material_type_defaults_to_metal_when_omitted(self, tmp_path):
        data = {k: v for k, v in CUSTOM_MATERIAL.items() if k != "materialType"}
        m = sr.load_material_file(write_material_file(tmp_path, data))
        assert m.definition["materialType"] == "metal"

    def test_non_positive_density_or_melting_point_is_rejected(self, tmp_path):
        with pytest.raises(sr.MaterialError, match="density"):
            sr.load_material_file(write_material_file(tmp_path, dict(CUSTOM_MATERIAL, density=0)))
        with pytest.raises(sr.MaterialError, match="meltingTemperature"):
            sr.load_material_file(write_material_file(tmp_path, dict(CUSTOM_MATERIAL, meltingTemperature=-5)))

    def test_missing_file_and_bad_json_are_material_errors(self, tmp_path):
        with pytest.raises(sr.MaterialError):
            sr.load_material_file(str(tmp_path / "nope.json"))
        bad = tmp_path / "bad.json"
        bad.write_text("{not json")
        with pytest.raises(sr.MaterialError):
            sr.load_material_file(str(bad))
        top_level_list = tmp_path / "list.json"
        top_level_list.write_text("[1, 2]")
        with pytest.raises(sr.MaterialError):
            sr.load_material_file(str(top_level_list))


TITANIUM = sr.Material("drama-TiAl6v4", 4417.0, 1905.0)


class TestMaterialInRun:
    def test_run_defaults_to_aa7075(self):
        assert make_run().material == sr.DEFAULT_MATERIAL

    def test_mass_uses_the_material_density(self):
        run = make_run(material=TITANIUM)
        assert run.mass_kg == pytest.approx(sr.sphere_mass_kg(50.0, 4417.0))
        assert run.mass_kg == pytest.approx(0.289092, rel=1e-5)   # 4417 kg/m3 * 4/3 pi (0.025 m)^3
        assert make_run().mass_kg == pytest.approx(sr.sphere_mass_kg(50.0))          # unchanged default

    def test_radius_from_mass_uses_the_material_density(self):
        run = make_run(material=TITANIUM)
        assert sr.radius_mm_from_mass(run.mass_kg, run.material) == pytest.approx(25.0)
        assert sr.radius_mm_from_mass(sr.sphere_mass_kg(50.0)) == pytest.approx(25.0)     # default still works

    def test_build_config_names_the_builtin_material_without_a_material_list(self):
        cfg = sr.build_config(make_run(material=TITANIUM))
        assert cfg["objects"][0]["material"] == "drama-TiAl6v4"
        assert cfg["objects"][0]["mass"] == pytest.approx(make_run(material=TITANIUM).mass_kg)
        assert "materialList" not in cfg

    def test_build_config_injects_a_custom_material_list(self):
        custom = sr.Material("user-AlLi2195", 2700.0, 823.0, source="custom_file", definition=CUSTOM_MATERIAL)
        cfg = sr.build_config(make_run(material=custom))
        assert cfg["objects"][0]["material"] == "user-AlLi2195"
        assert cfg["materialList"] == [CUSTOM_MATERIAL]
        assert cfg["materialList"][0] is not CUSTOM_MATERIAL     # a copy, so pyDRAMA cannot mutate ours

    def test_run_name_unchanged_for_the_default_material(self):
        assert sr.run_name(make_run()) == "sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km"

    def test_run_name_gets_a_material_suffix_otherwise(self):
        assert sr.run_name(make_run(material=TITANIUM)) == \
            "sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km_mdrama-TiAl6v4"
        custom = sr.Material("user AlLi 2195/T8", 2700.0, 823.0, source="custom_file", definition=CUSTOM_MATERIAL)
        assert sr.run_name(make_run(material=custom)).endswith("_muser_AlLi_2195_T8")   # filename-safe slug

    def test_melt_duration_uses_the_material_melting_point(self):
        # 1905 K plateau: at melt for titanium, never at melt for the AA7075 default
        aero = [aero_row(i, temp=T) for i, T in enumerate([300.0, 1905.0, 1905.0, 1905.0, 900.0])]
        traj = [traj_row(i) for i in range(5)]
        rows, _ = sr.merge_histories(aero, traj)
        log_info = {"end_of_life_reason": "ground impact", "event_end_mass_kg": None, "event_start_mass_kg": None}
        r_ti, _ = sr.compute_stats(make_run(material=TITANIUM), rows, None, log_info)
        assert r_ti["time_at_melting_temperature_s"] == pytest.approx(2.0)
        assert r_ti["melting_temperature_K"] == 1905.0
        r_al, _ = sr.compute_stats(make_run(), rows, None, log_info)
        assert r_al["time_at_melting_temperature_s"] == 0.0 and r_al["melting_temperature_K"] == 850.0

    def test_final_radius_statistic_uses_the_material_density(self):
        run = make_run(material=TITANIUM)
        aero = [aero_row(0.0, mass=run.mass_kg), aero_row(1.0, mass=run.mass_kg)]
        traj = [traj_row(0.0), traj_row(1.0)]
        rows, _ = sr.merge_histories(aero, traj)
        fragments = {"mass_kg": run.mass_kg, "velocity_kms": 0.05, "lat_deg": 0.0, "lon_deg": 0.0, "epoch": None}
        log_info = {"end_of_life_reason": "ground impact", "event_end_mass_kg": None, "event_start_mass_kg": None}
        r, _ = sr.compute_stats(run, rows, fragments, log_info)
        assert r["final_radius_mm"] == pytest.approx(25.0)
        assert r["initial_mass_kg"] == pytest.approx(run.mass_kg)
