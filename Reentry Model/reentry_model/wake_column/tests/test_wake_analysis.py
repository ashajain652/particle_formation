"""The sweep and sensitivity drivers and their plots (plan Task 10): the case lists, the settings set, a tiny sweep
through the plot scripts, and the selection of baseline runs from a folder that also holds other runs."""
import dataclasses
import json

import pytest

from reentry_model.wake_column import run
from reentry_model.wake_column.analysis import results, sensitivity, sweep
from reentry_model.wake_column.settings import Case, Settings, run_name

TINY = ["--workers", "1", "--t-end", "60", "--n-test-particles", "200"]


def test_default_sweep_is_126_cases():
    c = sweep.default_cases()
    assert len(c) == 126 and len({(x.d0_um, x.h0_km, x.time_of_day) for x in c}) == 126
    assert {x.d0_um for x in c} == {1.0, 14.1, 200.0}
    assert sorted({x.h0_km for x in c}) == [float(h) for h in range(50, 91, 2)]
    assert {x.time_of_day for x in c} == {"day", "night"}


def test_sensitivity_set_is_15_settings_over_6_cases():
    s = sensitivity.settings_set()
    assert len(s) == 15 and len({run_name(Case(14.1, 76.0, "day"), x) for _, x in s}) == 15
    assert len(sensitivity.cases()) == 6 and {c.h0_km for c in sensitivity.cases()} == {60.0, 76.0, 90.0}
    assert {c.d0_um for c in sensitivity.cases()} == {14.1}
    labels = dict(s)
    assert labels["droplet-temperature start"].hot_start == "droplet_temperature"
    assert labels["molecular only"].eddy_scale == 0.0
    # the two settings added by the controller's ruling (the 12 h result is dominated by smoke uptake)
    assert labels["smoke sticking x0.1"].smoke_sticking == 0.1
    assert labels["no smoke uptake"].smoke_sticking == 0.0
    # each setting changes exactly the field its label names, and no other
    default = Settings()
    changed = {label: [f.name for f in dataclasses.fields(Settings)
                       if getattr(st, f.name) != getattr(default, f.name)] for label, st in s}
    assert all(len(v) == 1 for v in changed.values())
    assert changed["formation x80"] == ["formation_scale"] and labels["formation x80"].formation_scale == 80.0
    assert labels["inner scale 10 m"].l_inner_m == 10.0 and labels["inner scale 100 m"].l_inner_m == 100.0
    assert labels["hot reactions x0.1"].hot_reaction_scale == 0.1 and labels["hot reactions x3"].hot_reaction_scale == 3.0
    assert labels["sticking x0.1"].sticking_scale == 0.1 and labels["oxidation x0.1"].oxidation_eff == 0.1
    assert labels["eddy x0.1"].eddy_scale == 0.1 and labels["eddy x100"].eddy_scale == 100.0
    assert labels["formation x0.1"].formation_scale == 0.1 and labels["formation x0.01"].formation_scale == 0.01


def test_settings_set_applies_common_overrides():
    s = sensitivity.settings_set(t_end_s=60.0, n_test_particles=200)
    assert len(s) == 15 and all(st.t_end_s == 60.0 and st.n_test_particles == 200 for _, st in s)


def test_baseline_signature_ignores_numerics_but_not_physics():
    base = dataclasses.asdict(Settings())
    shorter = dataclasses.asdict(Settings(t_end_s=60.0, n_test_particles=200, seed=1))
    assert results.physics_signature(base) == results.physics_signature(shorter)
    assert results.is_baseline({"settings": shorter})
    assert not results.is_baseline({"settings": dataclasses.asdict(Settings(formation_scale=0.1))})
    assert not results.is_baseline({"settings": dataclasses.asdict(Settings(hot_start="droplet_temperature"))})
    assert not results.is_baseline({"settings": dataclasses.asdict(Settings(smoke=False))})
    # the signature survives the JSON round trip that summary.json makes (tuples become lists)
    with_tuple = Settings(sticking=(("AlO", 0.5),))
    assert results.physics_signature(json.loads(json.dumps(dataclasses.asdict(with_tuple)))) \
        == results.physics_signature(dataclasses.asdict(with_tuple))


def test_load_summaries_skips_skipped_and_unreadable(tmp_path):
    good = tmp_path / "wake_good"
    good.mkdir()
    (good / "summary.json").write_text(json.dumps({"run_name": "wake_good", "condensed_fraction_end": 0.3,
                                                   "case": {"d0_um": 1.0, "h0_km": 80.0, "time_of_day": "day"},
                                                   "settings": dataclasses.asdict(Settings())}))
    skipped = tmp_path / "wake_skipped"
    skipped.mkdir()
    (skipped / "summary.json").write_text(json.dumps({"run_name": "wake_skipped", "skipped": "no aluminium vaporised"}))
    broken = tmp_path / "wake_broken"
    broken.mkdir()
    (broken / "summary.json").write_text('{"run_name": "wake_br')
    (tmp_path / "plots").mkdir()
    got = results.load_summaries(tmp_path)
    assert [s["run_name"] for s in got] == ["wake_good"] and got[0]["_dir"] == str(good)


def test_plot_scripts_run_on_a_tiny_sweep(tmp_path, parent):
    sweep.main(["--diameters", "1.0", "--altitudes", "80", "--out", str(tmp_path)] + TINY)
    from reentry_model.wake_column.analysis import plot_sweep
    plot_sweep.main(["--out", str(tmp_path)])
    pngs = {p.name for p in (tmp_path / "plots").glob("*.png")}
    assert pngs
    assert "condensed_vs_time_d01.0um.png" in pngs and "provenance_d01.0um.png" in pngs
    assert "size_distribution_d01.0um.png" in pngs and "day_minus_night.png" in pngs
    assert "waccm_vs_msis_afgl.png" not in pngs        # no NRLMSIS + AFGL runs in this folder: nothing to compare
    # a second call finds every run finished and writes nothing new to the sweep (resumable)
    names = sorted(p.name for p in tmp_path.iterdir())
    sweep.main(["--diameters", "1.0", "--altitudes", "80", "--out", str(tmp_path)] + TINY)
    assert sorted(p.name for p in tmp_path.iterdir()) == names


def test_plot_scripts_ignore_a_skipped_run(tmp_path, parent):
    sweep.main(["--diameters", "1.0", "--altitudes", "80", "--out", str(tmp_path)] + TINY)
    d = tmp_path / "wake_d200.00um_h050.0km_day_waccm_tend60_ntp200"
    d.mkdir()
    (d / "summary.json").write_text(json.dumps({
        "run_name": d.name, "skipped": "no aluminium vaporised",
        "case": {"d0_um": 200.0, "h0_km": 50.0, "time_of_day": "day"},
        "settings": dataclasses.asdict(Settings(t_end_s=60.0, n_test_particles=200))}))
    from reentry_model.wake_column.analysis import plot_sweep
    plot_sweep.main(["--out", str(tmp_path)])
    assert not any("d200" in p.name for p in (tmp_path / "plots").glob("*.png"))


def test_sensitivity_runs_and_tornado_plot(tmp_path, parent):
    rc = sensitivity.main(["--out", str(tmp_path), "--altitudes", "76", "--time-of-day", "night",
                           "--labels", "formation x0.1", "droplet-temperature start"] + TINY)
    assert rc == 0
    names = sorted(p.name for p in tmp_path.iterdir())
    assert len(names) == 3                                  # the baseline and two settings, one case
    got = {s["run_name"]: s for s in results.load_summaries(tmp_path)}
    hot = [s for s in got.values() if s["settings"]["hot_start"] == "droplet_temperature"]
    assert len(hot) == 1 and hot[0]["hot_fractions"] and hot[0]["share_above_AlO_bond"] == 0.0
    from reentry_model.wake_column.analysis import plot_sensitivity
    plot_sensitivity.main(["--out", str(tmp_path)])
    pngs = {p.name for p in (tmp_path / "plots").glob("*.png")}
    assert {"tornado_condensed_fraction.png", "tornado_volume_mean_diameter.png"} <= pngs


def test_a_failed_case_does_not_abort_the_sweep(tmp_path, monkeypatch):
    def fake_run(case, settings, out_dir, **kw):
        if case.h0_km == 80.0:
            raise RuntimeError("boom")
        return {"run_name": run_name(case, settings)}
    monkeypatch.setattr(run, "run", fake_run)
    jobs = [(Case(1.0, 80.0, "day"), Settings()), (Case(1.0, 78.0, "day"), Settings())]
    done, failed = sweep.run_jobs(jobs, str(tmp_path), workers=1)
    assert len(done) == 1 and len(failed) == 1 and "boom" in failed[0][1]
