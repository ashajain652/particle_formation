"""Spheral M2 Task 10b's measured tables (analysis/spheral_m2_task10b.py, measured 2026-10-12): the remedy for
SolidSPH's tensile pairing on the MVP flight's replay, the hourglass stiffness (packages/hourglass.py, run.py
--hourglass, default 0.01). Thresholds are the measured values rounded up to one significant figure (bracketed
values in the comments); they change only together with the plan's Task 10b section and the run.py docstring."""
import csv
import json
import os

import pytest

from helpers import REPO_ROOT


def table(name):
    with open(os.path.join(REPO_ROOT, "data", name), newline="") as fh:
        return list(csv.DictReader(fh))


def probe(name, n=None):
    rows = [r for r in table("m2_task10b_probe.csv") if r["probe"] == name]
    if n is not None:
        rows = [r for r in rows if int(r["crossings"]) == n]
    assert rows, name
    return rows


def test_the_reproducer_grows_per_step():
    """Check 13's RZ 1.1 mm interval 4 -> 5 at the runner's settings before Task 10b: the pairing grows over the
    interval x3.20, 7.85, 17.2 at 10, 20, 40 crossings, and the core's axial stress is lost at 40."""
    g = {int(r["crossings"]): float(r["growth"]) for r in probe("base")}
    assert g[10] > 3.0 and g[20] > 7.0 and g[40] > 15.0                 # [3.20, 7.85, 17.2]
    r40, = probe("base", 40)
    assert float(r40["core_sxx_Pa"]) < 0.75*float(r40["core_sxx_rz22_Pa"])   # [26.9 against 38.1 MPa]


def test_the_hourglass_stops_the_growth_at_10_20_and_40_crossings():
    """kappa 0.01 with its viscous part: growth x0.18, 0.08, 0.08; the core's axial stress within 0.05 % of RZ
    2.2 mm's; the internal KE ratio as without it (9.2e-6 at 40 crossings; the stiffness alone 9.9e-3)."""
    for r in probe("hg1e-2"):
        assert float(r["growth"]) < 0.2                                  # [0.18, 0.08, 0.08]
        assert abs(float(r["core_sxx_Pa"])/float(r["core_sxx_rz22_Pa"]) - 1.0) < 2e-3     # [5e-4]
    r40, = probe("hg1e-2", 40)
    assert float(r40["ke_internal_ratio"]) < 2e-5                        # [9.24e-6]
    s40, = probe("hg1e-2_stiff", 40)
    assert float(s40["ke_internal_ratio"]) > 5e-3                        # [9.86e-3]: the axis rows' oscillation
    # weaker stiffnesses only slow it
    assert all(float(r["growth"]) > 1.0 for r in probe("hg3e-3"))       # [1.31, 1.45, 1.61]


def test_the_alternatives():
    """nPerh 3.01 is stable (x0.69, 0.46, 0.65) but 25-35 % dearer per step in RZ; WendlandC6 at 2.01 is worse; the
    damping, XSPH and the unlimited Monaghan-Gingold viscosity leave it growing; small epsTensile makes it worse and
    loses the core's tension; epsTensile 0.3 cuts Spheral's step 21 times."""
    assert all(float(r["growth"]) < 1.0 for r in probe("nperh3"))
    assert float(probe("c6", 10)[0]["growth"]) > 5.0                     # [8.91]
    for name in ("zeta0.2", "xsph", "q_mg", "eps0.02", "eps0.05"):
        assert float(probe(name, 10)[0]["growth"]) > 1.5, name           # [3.14, 1.66, 3.45, 3.89, 9.70]
    assert float(probe("eps0.05", 10)[0]["core_sxx_Pa"]) < 1e7           # [7.9 MPa]
    t0, = probe("timing_eps0")
    t3, = probe("timing_eps0.3")
    assert int(t3["steps"]) > 15*int(t0["steps"])                       # [7,128 against 334]
    assert float(t3["ms_per_step"]) > 1.3*float(t0["ms_per_step"])      # [19.5 against 13.1 ms]


def test_the_hourglass_on_the_other_forms():
    """After check 13's switches (k 39 -> 40): 3D 2.2 mm at 20 crossings x1.21 without, x0.76 with, the core within
    0.3 % of the run's; RZ 2.2 mm does not grow there either way."""
    b, = probe("3d22_k39_base", 20)
    h, = probe("3d22_k39_hg1e-2", 20)
    assert float(b["growth"]) > 1.1 and float(h["growth"]) < 0.9         # [1.21, 0.76]
    assert abs(float(h["core_sxx_Pa"])/float(b["core_sxx_Pa"]) - 1.0) < 5e-3     # [2.7e-3]
    assert all(float(r["growth"]) < 0.6 for r in probe("rz22_k39_hg1e-2"))       # [0.33, 0.53]


def test_the_replays_from_k0():
    """--hourglass 0.01 from k 0: RZ 1.1 mm's pairing stays at 3e-4 dx and its core axial stress within 0.8 % of RZ
    2.2 mm's to k 12 (check 13: 0.24-0.34 of 3D's from k 12); RZ 2.2 mm and 3D 2.2 mm differ from check 13's runs by
    at most 4.7 % and 2.0 % of the peak (RMS 1.1 % and 0.23 %; the outer spacing)."""
    rows = table("m2_task10b_replay.csv")
    rz11 = [r for r in rows if r["form"] == "rz" and r["dx_mm"] == "1.1"]
    assert len(rz11) == 12
    assert max(float(r["pairing_dx"]) for r in rz11) < 5e-4                          # [3.1e-4]
    assert max(abs(float(r["core_sxx_Pa"])/float(r["core_sxx_rz22_base_Pa"]) - 1.0) for r in rz11) < 1e-2  # [8.2e-3]
    rz22 = [r for r in rows if r["form"] == "rz" and r["dx_mm"] == "2.2"]
    assert max(float(r["diff_vm_max"]) for r in rz22) < 0.05                          # [0.047]
    assert max(float(r["diff_vm_rms"]) for r in rz22) < 0.02                          # [0.011]
    d3 = [r for r in rows if r["form"] == "3d"]
    assert len(d3) == 3 and max(float(r["diff_vm_max"]) for r in d3) < 0.03         # [0.020]
    assert max(float(r["hourglass_energy_ratio"]) for r in rows) < 2e-3               # [1.5e-3]


def test_checks_4_to_6_are_not_degraded():
    """Checks 4-6 rerun with the hourglass against their committed rows: check 4's errors against the exact solution
    and check 6's against the finite-element reference are equal or smaller (largest increase: check 6's interior
    pressure, +0.0035 of the peak), check 5's RZ 1.1 mm plateau within 4.2e-4 Y of the committed one, check 6's drag
    unchanged (2e-11), the cost 1-10 % per step."""
    rows = table("m2_task10b_checks.csv")

    def pick(check, metric, form=None, dx=None):
        return [r for r in rows if r["check"] == check and r["metric"] == metric and
                (form is None or r["form"] == form) and (dx is None or r["dx_mm"] == dx)]
    for metric in ("interior_err", "shell_err", "interior_rms_particle"):
        assert all(float(r["difference"]) <= 0.0 for r in pick("check4", metric)), metric
    for metric in ("int_vm_max", "int_vm_rms", "core_vm_max", "core_vm_rms", "shell_vm_max", "shell_vm_rms"):
        assert all(float(r["difference"]) <= 0.0 for r in pick("check6", metric)), metric
    assert max(float(r["difference"]) for r in pick("check6", "int_p_max")) < 5e-3            # [3.5e-3]
    assert max(abs(float(r["difference"])) for r in pick("check6", "drag_rel")) < 1e-10      # [2.1e-11]
    assert max(abs(float(r["difference"])) for r in pick("check5", "sxx_err_over_Y", "rz", "1.1")) < 5e-4  # [4.1e-4]
    assert max(float(r["hourglass"]) for r in pick("check5", "lat_max_over_Y", "rz", "1.1")) < 2e-2       # [0.0139]
    ratios = [float(r["hourglass"])/float(r["committed"]) for r in rows if r["metric"] == "wall_per_step_ms"]
    assert max(ratios) < 1.2                                                                   # [1.105]


def test_the_runs_made_before_carry_hourglass_0():
    """Every committed table's runs were made before the default changed and are named _hourglass-0 (the
    migration of Task 10b, as 1d82f9e's for the area shares)."""
    for name, cols in (("m2_check4.csv", ("run",)), ("m2_check5.csv", ("run",)), ("m2_check6.csv", ("run",)),
                       ("m2_check7.csv", ("run",)), ("m2_monitor.csv", ("run",)), ("m2_check13_runs.csv", ("run",)),
                       ("m2_check13.csv", ("run_target", "run_source"))):
        for r in table(name):
            for c in cols:
                assert r[c].endswith("_hourglass-0"), (name, r[c][-80:])
