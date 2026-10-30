"""Sweep, sensitivity and tables (Task 9, spec 4.5, 4.6 and 5): the resumable sweep and its summary file, the
sensitivity set, and the supporting tables."""
import csv, filecmp, math, os, shutil
import pytest
from reentry_model.parent_follower_separation_interactions.analysis import sensitivity, sweep, tables

ARGS = ["--diameters", "10", "--altitudes", "60", "--thetas", "90", "180", "--bounds", "low",
        "--wake-flow", "off", "--workers", "1"]


def test_short_sweep_is_reproducible_and_resumable(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    sweep.main(ARGS + ["--out", str(a)]); sweep.main(ARGS + ["--out", str(b)])
    assert filecmp.cmp(a / "summary.csv", b / "summary.csv", shallow=False)
    rows = list(csv.DictReader(open(a / "summary.csv")))
    assert {r["affix_zone"] for r in rows} == {"far_wake", "escape"} and len(rows) == 2
    mtime = os.path.getmtime(next((a / "runs").glob("*.json")))
    sweep.main(ARGS + ["--out", str(a)])                                   # resumes: nothing rerun
    assert os.path.getmtime(next((a / "runs").glob("*.json"))) == mtime


def test_summary_carries_the_validity_flags_and_no_wall_clock(tmp_path):
    sweep.main(ARGS + ["--out", str(tmp_path)])
    with open(tmp_path / "summary.csv") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    needed = {"diameter_mm", "altitude_km", "theta_deg", "bound", "wake_flow", "affix_zone", "residence_time",
              "censored", "exit_altitude_m", "exit_speed", "exit_lag_m", "plume_distance_out_of_range",
              "plume_speed_out_of_range", "follower_regime", "affix_status", "release_kind", "shock_layer_exit",
              "t_shock_layer_exit", "kn_f_release", "kn_f_min_in_layer", "flux_mismatch", "flux_mismatch_above_tol",
              "follower_continuum_in_layer", "total_inside", "t_first_entry"}
    assert needed <= set(reader.fieldnames)
    assert not [c for c in reader.fieldnames if "time_s" in c and c != "t_shock_layer_exit"
                or c in ("runtime_s", "wall_time_s", "git_commit", "git_dirty")]
    # 10 mm at 60 km: at 90 degrees the shoulder gas is thin (transitional follower, released from the surface), at
    # 180 degrees the stagnation gas is dense (continuum, escape)
    by_theta = {float(r["theta_deg"]): r for r in rows}
    assert by_theta[90.0]["follower_regime"] == "transitional" and by_theta[90.0]["release_kind"] == "surface"
    assert by_theta[180.0]["follower_regime"] == "continuum" and by_theta[180.0]["release_kind"] == "escape"
    assert by_theta[90.0]["shock_layer_exit"] and not by_theta[180.0]["shock_layer_exit"]


def test_altitudes_outside_the_twin_are_dropped_with_a_message(tmp_path, capsys):
    sweep.main(["--diameters", "10", "--altitudes", "60", "20", "--thetas", "180", "--bounds", "low",
                "--wake-flow", "off", "--workers", "1", "--out", str(tmp_path)])
    assert "20" in capsys.readouterr().out
    assert {r["altitude_km"] for r in csv.DictReader(open(tmp_path / "summary.csv"))} == {"60.0"}


def test_a_failing_case_is_reported_and_the_rest_run(tmp_path, monkeypatch, capsys):
    real = sweep.sphere_in_spacecraft_plume
    def flaky(d_mm, h_km, theta, settings, twin=None):
        if theta == 90.0:
            raise RuntimeError("deliberate")
        return real(d_mm, h_km, theta, settings, twin=twin)
    monkeypatch.setattr(sweep, "sphere_in_spacecraft_plume", flaky)
    status = sweep.main(ARGS + ["--out", str(tmp_path)])
    assert status == 1 and "deliberate" in capsys.readouterr().out
    assert len(list(csv.DictReader(open(tmp_path / "summary.csv")))) == 1


def test_crossover_bracket_from_lags():
    # lag > 0: the sphere is behind the parent at 1 s; the bracket is [largest positive, smallest negative]
    low, high, monotone = sweep.crossover_bracket({1.0: 5.0, 10.0: 3.0, 100.0: 0.5, 200.0: -2.0, 400.0: -9.0})
    assert (low, high, monotone) == (100.0, 200.0, True)
    assert sweep.crossover_bracket({1.0: 5.0, 10.0: 3.0}) == (10.0, None, True)          # not bracketed from above
    assert sweep.crossover_bracket({1.0: 2.0, 10.0: -1.0, 50.0: 1.0, 100.0: -3.0})[2] is False   # lag changes sign twice


def test_sensitivity_set_is_the_spec_set():
    names = [v.parameter for v in sensitivity.variations()]
    assert names.count("taylor_K") == 2 and names.count("taylor_n") == 2 and names.count("turbulent_law") == 1
    assert names.count("T_axis_hot") == 2 and names.count("feldman_f") == 1
    assert names.count("frontal_area") == 2 and names.count("beta_scale") == 2 and names.count("density") == 2
    assert names.count("wake_flow") == 2
    base = sweep.ScenarioSettings()
    by = {(v.parameter, v.value): v.apply(base) for v in sensitivity.variations()}
    assert by[("taylor_K", 0.66 * 0.5 / 0.7)].law.taylor_K == pytest.approx(0.66 * 0.5 / 0.7)
    assert by[("turbulent_law", "park")].law.turbulent_law == "park"
    assert by[("density", 2700.0)].density == 2700.0
    assert by[("wake_flow", "velocity")].wake_flow == "velocity"
    rare = {(v.parameter, v.value): v.apply(base) for v in sensitivity.rarefied_variations()}
    assert rare[("shock_standoff", "billig")].shock.shock_standoff == "billig"
    assert rare[("pressure_gradient", False)].shock.pressure_gradient is False
    assert rare[("follower_wall_T_K", 1000.0)].shock.follower_wall_T_K == 1000.0
    assert rare[("shock_thickness_mfp", 1.0)].shock.shock_thickness_mfp == 1.0
    assert rare[("shock_thickness_mfp", 10.0)].shock.shock_thickness_mfp == 10.0


def test_short_sensitivity_reports_change_against_the_baseline(tmp_path):
    sensitivity.main(["--diameters", "10", "--altitudes", "60", "--theta", "120", "--parameters", "density",
                      "--no-rarefied", "--bounds", "low", "--workers", "1", "--out", str(tmp_path)])
    with open(tmp_path / "sensitivity.csv") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    assert reader.fieldnames[:8] == ["case", "parameter", "value", "bound", "residence_time", "change_vs_baseline",
                                     "shock_layer_exit", "t_shock_layer_exit"]
    assert [r["parameter"] for r in rows] == ["baseline", "density", "density"]
    base = float(rows[0]["residence_time"])
    assert float(rows[0]["change_vs_baseline"]) == 0.0
    for r in rows[1:]:
        assert float(r["change_vs_baseline"]) == pytest.approx(float(r["residence_time"]) - base, abs=1e-12)
    assert os.path.exists(tmp_path / "sensitivity_ranking.csv")


def test_tables_report_regimes_as_outputs(tmp_path):
    tables.main(["--out", str(tmp_path)])
    rows = list(csv.DictReader(open(tmp_path / "plume_size.csv")))
    assert {r["regime"] for r in rows} <= {"rarefied", "laminar", "transitional", "turbulent", "too_slow"}
    for name in ("twin_vs_drama.csv", "ballistic_coefficient.csv", "affix_regimes.csv", "drag_ratio_check.csv",
                 "wake_gas_sanity.csv"):
        assert (tmp_path / name).exists(), name
    assert "twin_downrange_km" in open(tmp_path / "twin_vs_drama.csv").read()


def test_table_contents_are_physically_consistent(tmp_path):
    tables.main(["--out", str(tmp_path)])
    plume = list(csv.DictReader(open(tmp_path / "plume_size.csv")))
    sizes = {(r["altitude_km"], float(r["s_m"])): r for r in plume}
    for r in plume:                                                        # time since passage is s / U
        assert float(r["t_since_passage_s"]) == pytest.approx(float(r["s_m"]) / float(r["parent_speed_ms"]))
    for alt in {r["altitude_km"] for r in plume if r["regime"] != "rarefied"}:
        lows = [float(sizes[(alt, s)]["r_b_low_m"]) for s in (10.0, 100.0, 1000.0, 10000.0)]
        assert lows == sorted(lows)                                        # the plume widens downstream
        assert float(sizes[(alt, 100.0)]["r_b_high_m"]) >= float(sizes[(alt, 100.0)]["r_b_low_m"])
    drama = [l for l in open(tmp_path / "twin_vs_drama.csv") if not l.startswith("#")]
    assert open(tmp_path / "twin_vs_drama.csv").readline().startswith("#")        # the header comment
    rows = list(csv.DictReader(drama))
    first = rows[0]
    assert float(first["twin_downrange_km"]) == pytest.approx(0.0, abs=1e-6)
    beta = list(csv.DictReader(open(tmp_path / "ballistic_coefficient.csv")))
    assert all(float(r["sphere_beta_kg_m2"]) > 0 and float(r["parent_beta_kg_m2"]) > 0 for r in beta)
    ratio = list(csv.DictReader(open(tmp_path / "drag_ratio_check.csv")))
    assert all(float(r["register_drag_ratio"]) == 0.14 for r in ratio)
    assert all(0.0 < float(r["drag_ratio"]) <= 1.0 for r in ratio if r["drag_ratio"])
    gas = list(csv.DictReader(open(tmp_path / "wake_gas_sanity.csv")))
    assert gas and all(float(r["core_speed_30km_ms"]) >= 0.0 for r in gas if r["core_speed_30km_ms"])


def test_crossover_rows_read_the_lag_one_second_after_release(tmp_path):
    from reentry_model.parent_follower_separation_interactions.scenarios import ScenarioSettings, run_name
    (tmp_path / "runs").mkdir()
    lags = {10.0: 3.0, 100.0: 0.5, 200.0: -2.0}
    rows = []
    for d, lag in lags.items():
        name = run_name("sphere", d, 80.0, ScenarioSettings(), release="far_wake", theta_deg=120.0)
        with open(tmp_path / "runs" / (name + ".csv"), "w") as fh:
            fh.write("t_since_release,lag_m\n0.0,5.6\n0.98,%r\n1.0,%r\n1.5,%r\n" % (lag + 0.1, lag, lag - 1))
        rows.append({"diameter_mm": d, "altitude_km": 80.0, "theta_deg": 120.0, "bound": "low", "wake_flow": "off",
                     "release_kind": "far_wake"})
    rows.append({"diameter_mm": 400.0, "altitude_km": 80.0, "theta_deg": 120.0, "bound": "low", "wake_flow": "off",
                 "release_kind": "surface"})                          # a rarefied follower is not a far-wake run
    (out,) = sweep.crossover_rows(str(tmp_path), rows)
    assert out["largest_diameter_positive_lag_mm"] == 100.0 and out["smallest_diameter_negative_lag_mm"] == 200.0
    assert out["lag_there_m"] == 0.5 and out["n_far_wake_diameters"] == 3 and out["monotone"] is True
