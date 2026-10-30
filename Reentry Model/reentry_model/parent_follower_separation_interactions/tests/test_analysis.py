"""Sweep, sensitivity and tables (Task 9, spec 4.5, 4.6 and 5): the resumable sweep and its summary file, the
sensitivity set, and the supporting tables."""
import concurrent.futures, csv, filecmp, glob, json, math, os
import pytest
from reentry_model.parent_follower_separation_interactions import scenarios
from reentry_model.parent_follower_separation_interactions.analysis import plots, sensitivity, sweep, tables

ARGS = ["--diameters", "10", "--altitudes", "60", "--thetas", "90", "180", "--bounds", "low",
        "--wake-flow", "off", "--workers", "1"]
FAR_WAKE_ARGS = ["--diameters", "25", "--altitudes", "60", "--thetas", "120", "130", "--bounds", "low",
                 "--wake-flow", "off", "--workers", "1"]


@pytest.fixture(autouse=True)
def clean_git(monkeypatch):
    """The cases record the commit and whether the tree is dirty, and a resumed sweep trusts only clean, current ones.
    The working tree is dirty while the code is edited, so the tests run as if on one clean commit."""
    monkeypatch.setattr(scenarios, "_git_commit", lambda: "test-commit")
    monkeypatch.setattr(scenarios, "_git_dirty", lambda: False)


def mtimes(folder):
    """{file: modification time in ns} of every file under folder/runs."""
    return {f: os.stat(f).st_mtime_ns for f in glob.glob(str(folder / "runs" / "*"))}


def test_short_sweep_is_reproducible_and_resumable(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    sweep.main(ARGS + ["--out", str(a)]); sweep.main(ARGS + ["--out", str(b)])
    assert filecmp.cmp(a / "summary.csv", b / "summary.csv", shallow=False)
    rows = list(csv.DictReader(open(a / "summary.csv")))
    assert {r["affix_zone"] for r in rows} == {"far_wake", "escape"} and len(rows) == 2
    before = mtimes(a)
    assert len([f for f in before if f.endswith(".json")]) == 2 and len(before) == 4      # two cases, JSON and history
    sweep.main(ARGS + ["--out", str(a)])                                   # resumes: nothing rerun, no file touched
    assert mtimes(a) == before


def test_summary_carries_the_validity_flags_and_no_wall_clock(tmp_path):
    sweep.main(ARGS + ["--out", str(tmp_path)])
    with open(tmp_path / "summary.csv") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    needed = {"name", "diameter_mm", "altitude_km", "theta_deg", "bound", "wake_flow", "affix_zone", "residence_time",
              "censored", "exit_altitude_m", "exit_speed", "exit_lag_m", "plume_distance_out_of_range",
              "plume_speed_out_of_range", "follower_regime", "affix_status", "release_kind", "shock_layer_exit",
              "t_shock_layer_exit", "kn_f_release", "kn_f_min_in_layer", "flux_mismatch", "flux_mismatch_above_tol",
              "follower_continuum_in_layer", "total_inside", "t_first_entry"}
    case_keys = [set(json.load(open(f))) for f in glob.glob(str(tmp_path / "runs" / "*.json"))]
    all_keys = set.union(*case_keys)
    assert needed <= all_keys and set(sweep.SUMMARY_COLUMNS) <= all_keys     # no column is a stale or misspelt name
    assert needed <= set(reader.fieldnames)
    assert (all_keys - sweep.SUMMARY_EXCLUDED - {"config"}) <= set(reader.fieldnames)    # nothing the cases hold is lost
    assert not {"runtime_s", "wall_time_s", "field_build_time_s", "git_commit", "git_dirty"} & set(reader.fieldnames)
    for r in rows:                                                         # the name is the deterministic run name
        assert r["name"] == scenarios.run_name("sphere", 10.0, 60.0, scenarios.ScenarioSettings(),
                                               release=r["release_kind"], theta_deg=float(r["theta_deg"]))
    # 10 mm at 60 km: at 90 degrees the shoulder gas is thin (transitional follower, released from the surface), at
    # 180 degrees the stagnation gas is dense (continuum, escape)
    by_theta = {float(r["theta_deg"]): r for r in rows}
    assert by_theta[90.0]["follower_regime"] == "transitional" and by_theta[90.0]["release_kind"] == "surface"
    assert by_theta[180.0]["follower_regime"] == "continuum" and by_theta[180.0]["release_kind"] == "escape"
    assert by_theta[90.0]["shock_layer_exit"] and not by_theta[180.0]["shock_layer_exit"]
    assert by_theta[90.0]["kn_f_min_in_layer"] and by_theta[90.0]["flux_mismatch_above_tol"] in ("True", "False")


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


def _case_files(folder):
    return sorted(glob.glob(str(folder / "runs" / "*.json")))


def test_stale_cases_are_rerun_and_force_reruns_everything(tmp_path, capsys, monkeypatch):
    args = ["--diameters", "10", "--altitudes", "60", "--thetas", "180", "--bounds", "low", "--wake-flow", "off",
            "--workers", "1", "--out", str(tmp_path)]
    sweep.main(args)
    (path,) = _case_files(tmp_path)
    before = os.stat(path).st_mtime_ns
    capsys.readouterr()
    for edit in ({"git_commit": "another-commit"}, {"git_dirty": True}, {"git_dirty": None}):
        summary = json.load(open(path))
        summary.update(edit)
        json.dump(summary, open(path, "w"))
        stamp = os.stat(path).st_mtime_ns
        sweep.main(args)
        assert os.stat(path).st_mtime_ns != stamp, edit                   # rerun, so the file was rewritten
        out = capsys.readouterr().out
        assert out.count("stale case") == 1 and "stale case sphere_d0010.0mm" in out, edit
        assert json.load(open(path))["git_commit"] == "test-commit"
    stamp = os.stat(path).st_mtime_ns
    sweep.main(args)                                                       # current again: not rerun, no stale line
    assert os.stat(path).st_mtime_ns == stamp and "stale case" not in capsys.readouterr().out
    sweep.main(args + ["--force"])
    assert os.stat(path).st_mtime_ns != stamp
    # the commit the cases are compared with is the current one
    monkeypatch.setattr(scenarios, "_git_commit", lambda: "next-commit")
    stamp = os.stat(path).st_mtime_ns
    sweep.main(args)
    assert os.stat(path).st_mtime_ns != stamp and "stale case" in capsys.readouterr().out


def _plant_orphan(folder, release, theta=180.0, commit="test-commit"):
    """A JSON and a history CSV for the 10 mm, 60 km, low, off case under another release kind than the real one."""
    name = scenarios.run_name("sphere", 10.0, 60.0, scenarios.ScenarioSettings(), release=release, theta_deg=theta)
    (folder / "runs").mkdir(parents=True, exist_ok=True)
    json.dump({"name": name, "git_commit": commit, "git_dirty": False, "residence_time": 99.0,
               "diameter_mm": 10.0, "altitude_km": 60.0, "theta_deg": theta, "bound": "low", "wake_flow": "off",
               "release_kind": release}, open(folder / "runs" / (name + ".json"), "w"))
    (folder / "runs" / (name + ".csv")).write_text("t_since_release,lag_m\n0.0,1.0\n")
    return name


def test_orphans_of_an_earlier_release_kind_are_deleted_and_never_summarised(tmp_path, capsys):
    args = ["--diameters", "10", "--altitudes", "60", "--thetas", "180", "--bounds", "low", "--wake-flow", "off",
            "--workers", "1", "--out", str(tmp_path)]
    real = scenarios.run_name("sphere", 10.0, 60.0, scenarios.ScenarioSettings(), release="escape", theta_deg=180.0)
    # (1) a stale orphan of an earlier classification and no real file: the case is run and the orphan removed
    stale = _plant_orphan(tmp_path, "far_wake", commit="old-commit")
    sweep.main(args)
    assert f"removed orphaned case {stale}" in capsys.readouterr().out
    left = sorted(os.path.basename(f) for f in glob.glob(str(tmp_path / "runs" / "*")))
    assert left == [real + ".csv", real + ".json"]
    (row,) = csv.DictReader(open(tmp_path / "summary.csv"))
    assert row["name"] == real and row["residence_time"] != "99.0"
    # (2) a current-looking orphan beside the current real file (a partly failed earlier run): it cannot be told from the
    # real result, so the case is rerun, and the orphan removed
    orphan = _plant_orphan(tmp_path, "surface")
    stamp = mtimes(tmp_path)[str(tmp_path / "runs" / (real + ".json"))]
    capsys.readouterr()
    sweep.main(args)
    out = capsys.readouterr().out
    assert "several release kinds" in out and f"removed orphaned case {orphan}" in out
    assert not glob.glob(str(tmp_path / "runs" / (orphan + "*")))
    assert os.stat(tmp_path / "runs" / (real + ".json")).st_mtime_ns != stamp
    text = open(tmp_path / "summary.csv").read()
    assert orphan not in text and "99.0" not in text and real in text
    # (3) a case with files under two release kinds is ambiguous and is left out of summary_rows
    _plant_orphan(tmp_path, "shock_band")
    assert sweep.summary_rows(str(tmp_path), [(10.0, 60.0, 180.0, "low", "off")]) == []


def test_a_dirty_tree_warns_once_and_reruns_everything_but_still_summarises(tmp_path, capsys, monkeypatch):
    args = ["--diameters", "10", "--altitudes", "60", "--thetas", "180", "--bounds", "low", "--wake-flow", "off",
            "--workers", "1", "--out", str(tmp_path)]
    monkeypatch.setattr(scenarios, "_git_dirty", lambda: True)
    outs = []
    for _ in range(2):
        capsys.readouterr()
        sweep.main(args)
        outs.append(capsys.readouterr().out)
    assert [o.count("WARNING: the working tree is dirty") for o in outs] == [1, 1]
    assert "stale case" not in outs[0] and "1 cases run, 0 resumed" in outs[0]
    assert "stale case" in outs[1] and "1 cases run, 0 resumed" in outs[1]     # the second run reran what the first wrote
    (path,) = _case_files(tmp_path)
    assert json.load(open(path))["git_dirty"] is True                          # stamped with the state of the run
    assert len(list(csv.DictReader(open(tmp_path / "summary.csv")))) == 1      # and its own case is summarised


def test_far_wake_angles_share_one_integration_and_resume(tmp_path, monkeypatch):
    twin = sweep.get_twin()
    for theta in (120.0, 130.0):                                           # both classify as far wake for 25 mm at 60 km
        assert scenarios.classify_release(25, 60, theta, twin, scenarios.ScenarioSettings()).zone == "far_wake"
    calls = []
    real = scenarios.dense_solve
    monkeypatch.setattr(scenarios, "dense_solve", lambda *a, **k: calls.append(1) or real(*a, **k))
    scenarios._FAR_WAKE_RUNS.clear()                                       # an earlier test may have left this case
    sweep.main(FAR_WAKE_ARGS + ["--out", str(tmp_path)])
    assert len(calls) == 1                                                 # one integration serves both thetas
    rows = list(csv.DictReader(open(tmp_path / "summary.csv")))
    assert [(r["theta_deg"], r["release_kind"], r["follower_regime"]) for r in rows] == [
        ("120.0", "far_wake", "continuum"), ("130.0", "far_wake", "continuum")]
    assert rows[0]["residence_time"] == rows[1]["residence_time"] and rows[0]["name"] != rows[1]["name"]
    before = mtimes(tmp_path)
    assert len([f for f in before if f.endswith(".json")]) == 2 and len(before) == 4
    sweep.main(FAR_WAKE_ARGS + ["--out", str(tmp_path)])
    assert mtimes(tmp_path) == before and len(calls) == 1                  # every file untouched, nothing integrated
    cases = [(25.0, 60.0, th, "low", "off") for th in (120.0, 130.0)]
    (cross,) = sweep.crossover_rows(str(tmp_path), sweep.summary_rows(str(tmp_path), cases))
    assert cross["n_far_wake_diameters"] == 1 and cross["lag_time_s"] == 1.0
    ends = [cross["largest_diameter_positive_lag_mm"], cross["smallest_diameter_negative_lag_mm"]]
    assert ends.count(25.0) == 1 and ends.count(None) == 1                 # the one diameter is on one side of zero
    lag = cross["lag_there_m"] if ends[0] == 25.0 else cross["lag_there_negative_m"]
    assert lag != 0.0 and (lag > 0.0) == (ends[0] == 25.0)
    out = list(csv.DictReader(open(tmp_path / "crossover.csv")))
    assert len(out) == 1 and out[0]["n_far_wake_diameters"] == "1"


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
    assert names.count("wake_flow") == 3
    base = sweep.ScenarioSettings()
    by = {(v.parameter, v.value): v.apply(base) for v in sensitivity.variations()}
    assert by[("taylor_K", 0.66 * 0.5 / 0.7)].law.taylor_K == pytest.approx(0.66 * 0.5 / 0.7)
    assert by[("turbulent_law", "park")].law.turbulent_law == "park"
    assert by[("density", 2700.0)].density == 2700.0
    assert by[("wake_flow", "velocity")].wake_flow == "velocity" and by[("wake_flow", "both")].wake_flow == "both"
    assert sensitivity.TAYLOR_K == sweep.ScenarioSettings().law.taylor_K
    rare = {(v.parameter, v.value): v.apply(base) for v in sensitivity.rarefied_variations()}
    assert rare[("shock_standoff", "billig")].shock.shock_standoff == "billig"
    assert rare[("pressure_gradient", False)].shock.pressure_gradient is False
    assert rare[("follower_wall_T_K", 1000.0)].shock.follower_wall_T_K == 1000.0
    assert rare[("shock_thickness_mfp", 1.0)].shock.shock_thickness_mfp == 1.0
    assert rare[("shock_thickness_mfp", 10.0)].shock.shock_thickness_mfp == 10.0


def test_short_sensitivity_reports_change_against_the_baseline(tmp_path):
    # an orphan of an earlier classification of the baseline case (10 mm, 60 km, 120 degrees is a far-wake release)
    orphan = scenarios.run_name("sphere", 10.0, 60.0, scenarios.ScenarioSettings(), release="escape", theta_deg=120.0)
    runs = tmp_path / "sensitivity" / "runs"
    runs.mkdir(parents=True)
    json.dump({"name": orphan, "git_commit": "old-commit", "git_dirty": False}, open(runs / (orphan + ".json"), "w"))
    sensitivity.main(["--diameters", "10", "--altitudes", "60", "--theta", "120", "--parameters", "density",
                      "--no-rarefied", "--bounds", "low", "--workers", "1", "--out", str(tmp_path)])
    assert not os.path.exists(runs / (orphan + ".json"))
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


SHORT_SENSITIVITY = ["--diameters", "10", "--altitudes", "60", "--theta", "120", "--parameters", "density",
                     "--bounds", "low"]


def test_dead_sensitivity_worker_is_reported_and_the_rest_written(tmp_path, monkeypatch, capsys):
    class FakePool:                       # a pool whose worker dies on the 2700 kg/m3 case
        def __init__(self, max_workers=None): pass
        def __enter__(self): return self
        def __exit__(self, *exc): return False
        def submit(self, fn, job):
            f = concurrent.futures.Future()
            if job[7] == 2700.0:
                f.set_exception(RuntimeError("worker died"))
            else:
                f.set_result(fn(job))
            return f
    monkeypatch.setattr(sensitivity.concurrent.futures, "ProcessPoolExecutor", FakePool)
    status = sensitivity.main(SHORT_SENSITIVITY + ["--no-rarefied", "--workers", "2", "--out", str(tmp_path)])
    assert status == 1 and "worker died" in capsys.readouterr().out
    rows = list(csv.DictReader(open(tmp_path / "sensitivity.csv")))
    assert [(r["parameter"], r["value"]) for r in rows] == [("baseline", ""), ("density", "2900.0")]


def test_failed_rarefied_baseline_skips_the_rarefied_set(tmp_path, monkeypatch, capsys):
    def broken(*a, **k):
        raise RuntimeError("shock layer failed")
    monkeypatch.setattr(sensitivity, "rarefied_diameter", broken)
    status = sensitivity.main(SHORT_SENSITIVITY + ["--workers", "1", "--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert status == 0 and out.count("skipped: RuntimeError: shock layer failed") == 2     # theta 180 and 120
    assert len(list(csv.DictReader(open(tmp_path / "sensitivity.csv")))) == 3


def test_rarefied_diameter_names_why_ten_mm_was_not_used(monkeypatch, tmp_path):
    def fake(regime, reason):
        def run_or_load(out_dir, d, h, theta, settings, force=False, git=None):
            return {"follower_regime": regime if d >= 10 else "transitional", "reason": reason}, True
        return run_or_load
    monkeypatch.setattr(sensitivity, "run_or_load", fake("continuum", "ok"))
    assert sensitivity.rarefied_diameter(90.0, 180.0, out_dir=str(tmp_path)) == (5.0, "a continuum follower")
    monkeypatch.setattr(sensitivity, "run_or_load", fake(None, "rarefied_no_law"))
    d, why = sensitivity.rarefied_diameter(90.0, 180.0, out_dir=str(tmp_path))
    assert d == 5.0 and "no plume law" in why and "rarefied_no_law" in why and "continuum" not in why
    monkeypatch.setattr(sensitivity, "run_or_load", lambda *a, **k: ({"follower_regime": "transitional"}, True))
    assert sensitivity.rarefied_diameter(90.0, 180.0, out_dir=str(tmp_path)) == (10.0, None)


def test_twin_time_lag_recovers_a_known_shift_below_90_km():
    from datetime import datetime, timedelta
    from types import SimpleNamespace
    t0 = datetime(2026, 1, 1)
    times = [100.0 + 5.0 * i for i in range(60)]                   # DRAMA rows 5 s apart
    def drama_h(T):                                                # km, falling 0.1 km/s from 111 km
        return 111.0 - 0.1 * (T - 100.0)
    lag, offset = 6.3, 20.0                                        # the twin starts 20 s after DRAMA's epoch
    parent = SimpleNamespace(epoch=t0, branch=[{"time": T, "altitude": drama_h(T)} for T in times])
    def altitude_m(t):                                             # twin time t is DRAMA time t + offset; it lags by `lag`
        h = drama_h(t + offset - lag)
        return 1e3 * (h + 5.0 if h > 92.0 else h)                  # distorted above 90 km: the fit must ignore those rows
    twin = SimpleNamespace(epoch=t0 + timedelta(seconds=offset), altitude=altitude_m,
                           track=SimpleNamespace(t0=-offset + times[0], t1=1e6))
    got, resid = tables.twin_time_lag(twin, parent)
    assert got == pytest.approx(lag, abs=1e-4) and resid < 1e-6
    whole, _ = tables.twin_time_lag(twin, parent, h_max_km=1e3)    # fitted over the whole branch it is another number
    assert abs(whole - lag) > 0.1


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
    assert float(first["twin_lag_s"]) == pytest.approx(0.0, abs=1e-6)         # the twin starts at the branch top
    lags = [float(r["twin_lag_s"]) for r in rows if r["twin_lag_s"] and float(r["drama_altitude_km"]) <= 90.0]
    assert lags and all(6.0 < lag < 7.0 for lag in lags)
    header = "".join(l for l in open(tmp_path / "twin_vs_drama.csv") if l.startswith("#"))
    assert "at and below 90 km" in header and "twin_lag_s" in header
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


# ---------------------------------------------------------------------------------------------------------- plots

def test_plots_are_written_from_a_short_sweep(tmp_path, capsys):
    sweep.main(ARGS + ["--out", str(tmp_path)]); tables.main(["--out", str(tmp_path)])
    capsys.readouterr()
    plots.main(["--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert "2 runs found" in out and "shock curve not drawn" not in out and "FAILED" not in out
    # ARGS: 10 mm at 60 km is a rarefied (surface) follower at 90 degrees and a continuum escape release at 180, with no
    # far-wake case, and no case at the map's default 120 degrees (the nearest, 90, is used and said so)
    for name in ("plume_size.png", "residence_map.png", "margin_60km.png", "trajectories_60km.png",
                 "affix_zones.png", "shock_layer_paths_60km.png"):
        assert (tmp_path / "plots" / name).stat().st_size > 0, name
    assert "no case at theta 120 degrees; using the nearest, 90" in out


def _row(**kw):
    base = {"release_kind": "far_wake", "residence_time": 2.0, "first_interval": 2.0, "t_first_entry": 0.0,
            "total_inside": 2.0, "censored": False, "reason": "ok", "diameter_mm": 10.0, "altitude_km": 60.0,
            "theta_deg": 120.0, "bound": "low", "wake_flow": "off"}
    base.update(kw)
    return base


def test_map_cell_shows_residence_for_far_wake_and_total_inside_for_releases_that_start_outside():
    nan = float("nan")
    assert plots.map_cell(_row())["value"] == 2.0
    # an escape or surface release starts outside the plume: residence_time is 0 and the first interval may begin late,
    # so the total time inside is shown
    escape = plots.map_cell(_row(release_kind="escape", residence_time=0.0, first_interval=5.0, t_first_entry=0.4,
                                 total_inside=6.5, reason="starts_outside"))
    assert escape == {"value": 6.5, "kind": "escape", "flags": []}
    surf = plots.map_cell(_row(release_kind="surface", residence_time=0.0, first_interval=1.1, t_first_entry=17.0,
                               total_inside=3.0))
    assert surf["value"] == 3.0 and surf["kind"] == "surface"
    never = plots.map_cell(_row(release_kind="surface", residence_time=0.0, first_interval=nan, t_first_entry=nan,
                                total_inside=nan))
    assert never["value"] == 0.0
    assert plots.map_cell(_row(reason="reached_parent", censored=True))["flags"] == ["c", "r"]
    assert plots.map_cell(_row(reason="hit_parent", release_kind="surface"))["flags"] == ["h"]
    assert plots.map_cell(_row(release_kind="", reason="rarefied_no_law", residence_time=nan))["flags"] == ["n"]
    assert plots.map_cell(_row(release_kind="", reason="shock_band", residence_time=nan))["flags"] == ["b"]
    assert math.isnan(plots.map_cell(_row(release_kind="", residence_time=nan))["value"])      # not integrated: NaN


def test_parent_frame_path_recovers_x_forward_and_y_from_the_axis():
    import numpy as np
    t = np.linspace(0.0, 2.0, 41)
    e = np.array([0.6, 0.0, 0.8])                          # the parent flies along e at 4,000 m/s from a point on the Earth
    P = np.array([6.5e6, 0.0, 0.0]) + 4000.0 * t[:, None] * e
    lateral = np.array([0.0, 1.0, 0.0])                    # perpendicular to e
    x_true, y_true = -3.0 + 2.0 * t, 1.4 + 0.5 * t
    F = P + x_true[:, None] * e + y_true[:, None] * lateral
    hist = {"t": t, "parent_x": P[:, 0], "parent_y": P[:, 1], "parent_z": P[:, 2],
            "follower_x": F[:, 0], "follower_y": F[:, 1], "follower_z": F[:, 2]}
    x, y = plots.parent_frame_path(hist)
    assert np.allclose(x, x_true, atol=1e-6) and np.allclose(y, y_true, atol=1e-6)


def test_an_altitude_without_an_integrated_case_prints_a_note_and_writes_nothing(tmp_path, capsys):
    rows = [_row(release_kind="", altitude_km=70.0, reason="shock_band", residence_time=float("nan"))]
    assert plots.plot_margin(str(tmp_path), rows, 70.0, "low") is None
    assert plots.plot_trajectories(str(tmp_path), rows, 70.0, "low", "off") is None
    assert plots.plot_shock_layer_paths(str(tmp_path), rows, 70.0, "low", "off") is None
    out = capsys.readouterr().out
    assert "margin_70km.png skipped" in out and "trajectories_70km.png skipped" in out
    assert "shock_layer_paths_70km.png skipped" in out
    assert not (tmp_path / "plots").exists()


def test_affix_zones_has_one_panel_unless_the_flags_differ_by_altitude():
    def rows(flag_at_50):
        out = []
        for h in (80.0, 50.0):
            for d in (10.0, 100.0):
                out.append({"altitude_km": str(h), "diameter_mm": str(d), "theta_deg": "120.0", "zone": "far_wake",
                            "theta_bif": "139.0", "flag_not_hypersonic": str(flag_at_50 and h == 50.0),
                            "flag_parent_bulk_density": "63.5"})
        return out
    same = plots.affix_panels(rows(False))
    assert [title for title, _ in same] == ["all altitudes"]
    differ = plots.affix_panels(rows(True))
    assert [title for title, _ in differ] == ["80 km", "50 km"] and all(len(rs) == 2 for _, rs in differ)
    assert plots.affix_flag_names(rows(True)) == ["flag_not_hypersonic"]       # a numeric column is not a flag


SUMMARY_FIELDS = ["name", "diameter_mm", "altitude_km", "theta_deg", "bound", "wake_flow", "release_kind", "reason",
                  "residence_time", "first_interval", "t_first_entry", "total_inside", "censored"]


def _fake_output(folder, affix_rows=True, history=False):
    """A folder with a summary of one far-wake and one not-integrated case, a plume_size.csv and an affix_regimes.csv
    (header only unless affix_rows); the far-wake case has a history file only if `history`."""
    folder.mkdir(parents=True, exist_ok=True)
    rows = [dict(zip(SUMMARY_FIELDS, ["ghost_far", 10.0, 60.0, 120.0, "low", "off", "far_wake", "ok", 2.0, 2.0, 0.0,
                                      2.0, False])),
            dict(zip(SUMMARY_FIELDS, ["ghost_band", 100.0, 60.0, 135.0, "low", "off", "", "shock_band", "", "", "", "",
                                      False]))]
    with open(folder / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, SUMMARY_FIELDS); w.writeheader(); w.writerows(rows)
    with open(folder / "plume_size.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["altitude_km", "regime", "s_m", "t_since_passage_s", "parent_speed_ms", "r_b_low_m", "r_b_high_m"])
        for s_m in (10.0, 100.0, 1000.0):
            w.writerow([60.0, "turbulent", s_m, s_m / 4000.0, 4000.0, 1.4 + s_m / 500.0, 1.5 + s_m / 400.0])
    with open(folder / "affix_regimes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["diameter_mm", "altitude_km", "theta_deg", "zone", "theta_bif", "flag_not_hypersonic"])
        if affix_rows:
            w.writerow([10.0, 60.0, 120.0, "far_wake", 139.0, "False"])
    if history:
        (folder / "runs").mkdir()
        with open(folder / "runs" / "ghost_far.csv", "w") as fh:
            fh.write("t,t_since_release,follower_x,follower_y,follower_z,parent_x,parent_y,parent_z,lag_m,distance_m,r_b,margin,regime\n")
            for i in range(5):
                t = 0.1 * i
                fh.write(f"{t},{t},{6.5e6 + 4000 * t},0,0,{6.5e6 + 4000 * t + 5},0,0,5,{1.4 + t},{1.4 + 2 * t},1,turbulent\n")


def test_a_missing_summary_returns_one(tmp_path, capsys):
    assert plots.main(["--out", str(tmp_path)]) == 1
    assert "no summary.csv" in capsys.readouterr().out and not (tmp_path / "plots").exists()


def test_a_header_only_affix_table_skips_that_plot_only(tmp_path, capsys):
    _fake_output(tmp_path, affix_rows=False)
    assert plots.main(["--out", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "affix_zones.png skipped" in out and "header and no rows" in out
    names = {f.name for f in (tmp_path / "plots").iterdir()}
    assert names == {"plume_size.png", "residence_map.png", "margin_60km.png", "trajectories_60km.png"}


def test_a_missing_run_history_is_named_and_the_other_plots_are_written(tmp_path, capsys):
    _fake_output(tmp_path)
    assert plots.main(["--out", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "runs/ghost_far.csv not found" in out
    assert {"plume_size.png", "residence_map.png", "affix_zones.png", "margin_60km.png"} <= {
        f.name for f in (tmp_path / "plots").iterdir()}


def test_a_failing_plot_is_reported_and_the_others_still_written(tmp_path, capsys, monkeypatch):
    _fake_output(tmp_path, history=True)
    def broken(*a, **k):
        raise RuntimeError("deliberate")
    monkeypatch.setattr(plots, "plot_margin", broken)
    assert plots.main(["--out", str(tmp_path)]) == 1
    out = capsys.readouterr().out
    assert "FAILED margin_60km.png: RuntimeError: deliberate" in out and "1 plot(s) failed" in out
    names = {f.name for f in (tmp_path / "plots").iterdir()}
    assert {"plume_size.png", "trajectories_60km.png", "residence_map.png", "affix_zones.png"} <= names
    assert "margin_60km.png" not in names


def test_shock_curve_rebuild_failure_is_a_printed_note_not_an_exception(tmp_path, capsys):
    R, shape, note = plots.rebuild_shock_curve(str(tmp_path), "no_such_case", 60.0)     # no JSON to read
    assert shape is None and note.startswith("shock curve not drawn") and R == pytest.approx(1.4065, abs=1e-3)
    # a rarefied (surface) case with a history but no run JSON: the plot is still written and the note is printed
    (tmp_path / "runs").mkdir()
    with open(tmp_path / "runs" / "ghost_surface.csv", "w") as fh:
        fh.write("t,t_since_release,follower_x,follower_y,follower_z,parent_x,parent_y,parent_z,lag_m,distance_m,r_b,margin,regime\n")
        for i in range(5):
            t = 0.1 * i
            fh.write(f"{t},{t},{6.5e6 + 4000 * t + 1.4 - 3 * t},{0.1 * t},0,{6.5e6 + 4000 * t},0,0,0,1,1,0,turbulent\n")
    row = _row(name="ghost_surface", release_kind="surface", shock_layer_exit="through_shock", t_shock_layer_exit=0.2,
               x_exit_parent_frame=1.2, y_exit_parent_frame=0.5)
    capsys.readouterr()
    assert plots.plot_shock_layer_paths(str(tmp_path), [row], 60.0, "low", "off") is not None
    assert "note: shock_layer_paths_60km.png: shock curve not drawn" in capsys.readouterr().out
    assert (tmp_path / "plots" / "shock_layer_paths_60km.png").stat().st_size > 0
