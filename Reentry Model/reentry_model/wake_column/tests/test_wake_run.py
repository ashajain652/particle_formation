"""One-case runner: the outputs it writes, reuse of a finished run, the empty-inventory case and the steady-state
time (plan Task 9; decisions 17 and 18; Review Focus 1)."""
import csv
import json

import numpy as np
import pytest

from reentry_model.wake_column import run
from reentry_model.wake_column.settings import Case, Settings
from reentry_model.wake_column.trail import Trail


def test_short_run_writes_all_outputs(tmp_path, parent):
    st = Settings(t_end_s=60.0, n_test_particles=500, n_emission_samples=4)
    s = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    d = tmp_path / s["run_name"]
    assert (d / "history.csv").exists() and (d / "size_distribution.csv").exists()
    assert (d / "summary.json").exists()
    for k in ("condensed_fraction_end", "t50_s", "provenance_end", "hot_fractions", "r0_m", "steady_state_time_s"):
        assert k in s
    assert 0.0 <= s["condensed_fraction_end"] <= 1.0


def test_short_run_file_contents(tmp_path, parent):
    st = Settings(t_end_s=60.0, n_test_particles=500, n_emission_samples=4)
    s = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    d = tmp_path / s["run_name"]
    with open(d / "history.csv", newline="") as f:
        hist = list(csv.DictReader(f))
    cols = hist[0].keys()
    for c in ("t_s", "r_m", "V_m3", "T_K", "frac_Al", "frac_AlOH2", "frac_particles", "frac_smoke", "frac_settled",
              "frac_condensed", "frac_AlO_formation", "frac_OAlO2_growth", "frac_AlCO3_smoke",
              "c_O_m3", "c_O2_m3", "c_O3_m3", "c_H_m3", "c_H2O_m3"):
        assert c in cols
    assert sum(c.startswith("frac_") for c in cols) == 7 + 4 + 21
    # the aluminium is conserved: gas + particles + smoke + settled = 1 at every output
    gas = ["frac_" + n for n in ("Al", "AlO", "OAlO", "OAlO2", "AlCO3", "AlOH", "AlOH2")]
    for row in hist:
        tot = sum(float(row[c]) for c in gas) + float(row["frac_particles"]) + float(row["frac_smoke"]) \
            + float(row["frac_settled"])
        assert tot == pytest.approx(1.0, abs=1e-5)
        assert float(row["frac_condensed"]) == pytest.approx(
            float(row["frac_particles"]) + float(row["frac_smoke"]) + float(row["frac_settled"]), abs=1e-12)
    # the final condensed fraction in the file matches the summary
    assert float(hist[-1]["frac_condensed"]) == pytest.approx(s["condensed_fraction_end"], abs=1e-12)
    # provenance keys are "source/route", 21 of them, summing to the condensed fraction
    assert len(s["provenance_end"]) == 21 and "AlO/formation" in s["provenance_end"]
    assert sum(s["provenance_end"].values()) == pytest.approx(s["condensed_fraction_end"], abs=1e-9)

    with open(d / "size_distribution.csv", newline="") as f:
        sd = list(csv.DictReader(f))
    times = sorted({float(r["t_s"]) for r in sd})
    assert len(times) == 3 and times[-1] == pytest.approx(60.0)       # 1e-3 s, 1 s, 60 s (later targets > t_end)
    assert "Al_AlO_formation" in sd[0] and "Al_OAlO2_growth" in sd[0]
    # masked bins have empty composition columns but keep N and a diameter
    masked = [r for r in sd if r["O_per_Al"] == ""]
    assert masked and all(r["H_per_Al"] == "" and r["alumina_like"] == "" and r["d_nm"] != "" for r in masked)
    # the summary JSON is valid with no NaN tokens
    json.loads((d / "summary.json").read_text(), parse_constant=lambda c: pytest.fail(f"bad token {c}"))


def test_existing_summary_is_reused(tmp_path, parent):
    st = Settings(t_end_s=1.0, n_test_particles=200, n_emission_samples=2)
    a = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    b = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    assert a["runtime_s"] == b["runtime_s"]


def test_zero_inventory_is_skipped_cleanly(tmp_path):           # Review Focus 1
    t = Trail(1.0, 50.0, 0.0, 0.0, 0.0, 1e-15, 1e-16, 50.0, ())
    s = run.run(Case(1.0, 50.0, "day"), Settings(), out_dir=tmp_path, trail=t)
    assert s["skipped"] == "no aluminium vaporised"
    assert s["m_Zn_kg"] == 1e-15 and s["m_Mg_kg"] == 1e-16
    assert (tmp_path / s["run_name"] / "summary.json").exists()
    assert not (tmp_path / s["run_name"] / "history.csv").exists()


def test_steady_state_time():
    t = np.logspace(-6, 4.6, 107)
    flat = {"a": np.where(t < 10, t / 10, 1.0)}
    assert run.steady_state_time(t, flat) == pytest.approx(10.0, rel=0.3)
    assert run.steady_state_time(t, {"a": np.log10(t) / 5}) is None


def test_steady_state_needs_a_decade_of_evidence():
    t = np.concatenate(([0.0], np.logspace(-6, 0, 61)))          # span of six decades, but flat from the start
    assert run.steady_state_time(t, {"a": np.ones_like(t)}) == pytest.approx(1e-6)
    short = np.array([0.0, 1e-6, 5e-6])                           # less than one decade after the first point
    assert run.steady_state_time(short, {"a": np.ones(3)}) is None


@pytest.mark.parametrize("junk", ['{"run_name": "x", "t50', "", "[1, 2]", '{"no_name": 1}'])
def test_unreadable_summary_is_rerun(tmp_path, junk):
    t = Trail(1.0, 50.0, 0.0, 0.0, 0.0, 1e-15, 1e-16, 50.0, ())
    case, st = Case(1.0, 50.0, "day"), Settings()
    d = tmp_path / run.run_name(case, st)
    d.mkdir()
    (d / "summary.json").write_text(junk)
    s = run.run(case, st, out_dir=tmp_path, trail=t)
    assert s["skipped"] == "no aluminium vaporised" and s["run_name"] == d.name
    assert json.loads((d / "summary.json").read_text())["run_name"] == d.name


def test_truncated_summary_of_a_real_run_is_rerun(tmp_path, parent):
    st = Settings(t_end_s=1.0, n_test_particles=200, n_emission_samples=2)
    case = Case(1.0, 80.0, "night")
    a = run.run(case, st, out_dir=tmp_path)
    path = tmp_path / a["run_name"] / "summary.json"
    path.write_text(path.read_text()[:200])                      # cut off mid-file
    b = run.run(case, st, out_dir=tmp_path)
    assert b["run_name"] == a["run_name"] and 0.0 <= b["condensed_fraction_end"] <= 1.0
    assert json.loads(path.read_text())["runtime_s"] == b["runtime_s"]


def test_t50_and_t90_are_the_first_output_times_reaching_the_level():     # decision 17
    t = np.array([0.0, 1.0, 10.0, 100.0, 1000.0])
    condensed = np.array([0.0, 0.2, 0.5, 0.89, 0.95])
    assert run._first_time_reaching(t, condensed, 0.5) == 10.0       # reaching the level exactly counts
    assert run._first_time_reaching(t, condensed, 0.9) == 1000.0     # 0.89 at t = 100 s has not reached 0.9
    assert run._first_time_reaching(t, condensed, 0.99) is None      # never reached: null
    assert run._first_time_reaching(t, condensed[::-1] * 0.0, 0.5) is None


def test_condensed_is_particles_plus_smoke_plus_settled():                # decision 17
    budget = {"Al": 0.1, "AlO": 0.2, "AlOH": 0.3, "particles": 0.05, "smoke": 0.25, "settled": 0.1}
    assert run._condensed(budget) == pytest.approx(0.4)             # no gas species in it
