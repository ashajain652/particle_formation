"""Check 5 (Spheral M2 plan, Task 7; spec §12): a heated bar held at both ends, heated uniformly by 300 K over 10
intervals, Y = 100 MPa, mass-scaled at the default crossings.

The normal-loop tests read the committed table `data/spheral/m2_check5.csv` (written by
`analysis/spheral_m2_check5.py`, whose docstring defines every measure, the expectation and its two
qualifications: the held ends' lateral constraint raises the middle's elastic stress by the finite-element factor
`fe_ratio`, and perfect plasticity leaves the distribution of plastic flow along the bar open, so the middle third's
plastic strain is held to its increment per interval and to the kinematic definition check, with its lag behind
alpha_eff dT - Y/E bounded). The RZ 2.2 mm rod collapses (SolidSPHRZ, five rows across the radius; recorded with
`collapsed` 1 and tested to be flagged, not held to the thresholds). The `spheral` test reruns RZ 1.1 mm over frames
0-3 and requires the committed rows back.

THRESHOLDS: measured 2026-10-10 (the CSV), each the worst case over the frames of its form and spacing, rounded up to
one significant figure; worst case in brackets. Change them only together with the CSV, the README and
docs/spheral_frag_assumptions.md (CLAUDE.md)."""
import csv
import os
import shutil
import subprocess
import sys

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))

CSV = os.path.join(REPO_ROOT, "data", "spheral", "m2_check5.csv")
Y = 100e6
E = 71.7e9
ALPHA_EFF = 2.34e-5/(1.0 + 2.34e-5*7.0)
INCREMENT = ALPHA_EFF*30.0                     # 7.02e-4: the uniaxial plastic strain per interval once yielded
STABLE = [("rz", 1.1), ("3d", 2.2), ("3d", 1.1)]

# --- thresholds (fractions of Y unless stated), over the stable runs; worst case in brackets
ELASTIC_ERR = {("rz", 1.1): 5e-3,      # [-4.5e-3] k 1: (mean sigma_xx - fe_ratio (-E alpha_eff dT)) / Y
               ("3d", 2.2): 2e-2,      # [-1.46e-2]
               ("3d", 1.1): 4e-3}      # [-3.9e-3]
TRANSITION_ERR = 2e-2                  # [+1.01e-2, 3D 1.1 mm] k 2 (dT 60 K; the middle reaches Y at about 56 K)
PLATEAU_ERR = 2e-3                     # [1.35e-3, 3D 2.2 mm] k >= 3: |mean sigma_xx + Y| / Y
SXX_SD = {("rz", 1.1): 5e-3,           # [4.5e-3] k >= 3: the middle third's sigma_xx spread / Y
          ("3d", 2.2): 4e-2,           # [3.08e-2]
          ("3d", 1.1): 6e-3}           # [5.98e-3]
LATERAL_MAX = {("rz", 1.1): 2e-2,      # [1.81e-2] every k: largest |lateral stress| in the middle third / Y
               ("3d", 2.2): 0.2,       # [0.108] the 9-particle section's corners
               ("3d", 1.1): 4e-2}      # [3.09e-2]
LATERAL_MEAN = 2e-3                    # [1.71e-3, 3D 2.2 mm] every k: |volume-mean lateral stress| / Y
KE_INTERNAL = 2e-4                     # [1.18e-4, 3D 1.1 mm, k 1] internal KE / elastic energy
MONITOR = 4e-2                         # [3.78e-2, 3D 2.2 mm, k 9->10, pressure] below the 5 % limit: no switch
INCREMENT_ERR = 0.2                    # [0.136, RZ 1.1 mm] 1.1 mm, k >= 5: |ps increment / 7.02e-4 - 1|
KINEMATIC = {("rz", 1.1): 5e-2,        # [4.7e-2] k >= 4: |ps_mean - ps_kinematic| / ps_expected
             ("3d", 2.2): 0.2,         # [0.114]
             ("3d", 1.1): 8e-2}        # [7.3e-2]
LAG = {("rz", 1.1): 8e-4,              # [7.06e-4] k 10: alpha_eff dT - Y/E - ps_mean (the ends' early flow)
       ("3d", 2.2): 3e-3,              # [2.72e-3]
       ("3d", 1.1): 6e-4}              # [5.28e-4]


def rows():
    with open(CSV, newline="") as fh:
        out = []
        for r in csv.DictReader(fh):
            for c in r:
                if c not in ("form", "run"):
                    r[c] = float(r[c])
            for c in ("k", "collapsed", "n_particles", "n_mid", "n_held", "processes", "monitored", "run_crossings",
                      "steps"):
                r[c] = int(r[c])
            out.append(r)
    return out


def pick(rs, **kw):
    return [r for r in rs if all(r[k] == v for k, v in kw.items())]


def stable(key):
    return sorted(pick(rows(), form=key[0], dx_mm=key[1]), key=lambda r: r["k"])


def test_table_is_complete():
    import spheral_m2_check5 as c5
    rs = rows()
    assert {(r["form"], r["dx_mm"], r["k"]) for r in rs} == {(f, dx, k) for f, dx, _ in c5.FORMS for k in range(1, 11)}
    for r in rs:
        assert r["dT_K"] == 30.0*r["k"]
        assert r["ps_expected"] == pytest.approx(max(ALPHA_EFF*r["dT_K"] - Y/E, 0.0), rel=1e-12, abs=1e-15)
        assert 1.05 < r["fe_ratio"] < 1.07          # the finite-element middle-third factor (fe.json)


def test_holds_are_exact():
    """The Hold package: every held particle at rest and at its reference position, bitwise, in every run."""
    for r in rows():
        assert r["held_max_speed_ms"] == 0.0 and r["held_max_disp_m"] == 0.0, r["run"]
        assert r["n_held"] > 0


def test_rz_22_collapses_and_is_flagged():
    """Five rows across the RZ 2.2 mm rod: the axis rows slide against the outer ones (the driver's docstring)."""
    rz = sorted(pick(rows(), form="rz", dx_mm=2.2), key=lambda r: r["k"])
    assert all(r["collapsed"] == 1 for r in rz)
    assert max(r["ke_internal_ratio"] for r in rz) > 0.05
    assert max(abs(r["sxx_err_over_Y"]) for r in rz if r["k"] >= 3) > 0.2
    assert rz[-1]["run_crossings"] == 40                        # the monitor doubled twice
    assert all(r["collapsed"] == 0 for key in STABLE for r in stable(key))


@pytest.mark.parametrize("key", STABLE)
def test_axial_stress(key):
    rs = stable(key)
    assert abs(rs[0]["sxx_err_over_Y"]) <= ELASTIC_ERR[key]
    assert abs(rs[1]["sxx_err_over_Y"]) <= TRANSITION_ERR
    for r in rs[2:]:
        assert abs(r["sxx_mean_Pa"] + Y)/Y <= PLATEAU_ERR
        assert r["sxx_sd_Pa"]/Y <= SXX_SD[key]
    # the finite-element factor matters: the uniaxial -E alpha dT alone is off by 5-10 % at k 1
    assert abs(rs[0]["sxx_mean_Pa"]/(-E*ALPHA_EFF*30.0) - 1.0) > 0.04


@pytest.mark.parametrize("key", STABLE)
def test_lateral_stresses(key):
    for r in stable(key):
        assert r["lat_max_over_Y"] <= LATERAL_MAX[key]
        assert abs(r["lat_mean_over_Y"]) <= LATERAL_MEAN


@pytest.mark.parametrize("key", STABLE)
def test_plastic_strain(key):
    """Elastic at k 1; at k >= 4 the particles' mean equals the kinematic eps_p (the radial return's equivalent
    plastic strain is the uniaxial plastic shortening); it lags alpha_eff dT - Y/E by the ends' early flow; at
    1.1 mm the increment per interval is the uniaxial 7.02e-4."""
    rs = stable(key)
    assert rs[0]["ps_max"] == 0.0
    for r in rs[3:]:
        assert abs(r["ps_mean"] - r["ps_kinematic"])/r["ps_expected"] <= KINEMATIC[key]
    lag = rs[-1]["ps_expected"] - rs[-1]["ps_mean"]
    assert 0.0 < lag <= LAG[key]
    if key[1] == 1.1:
        for r in rs[4:]:
            assert abs(r["ps_increment"]/INCREMENT - 1.0) <= INCREMENT_ERR


@pytest.mark.parametrize("key", STABLE)
def test_slowness_and_monitor(key):
    rs = stable(key)
    assert max(r["ke_internal_ratio"] for r in rs) <= KE_INTERNAL
    mon = [r for r in rs if r["monitored"]]
    assert [r["k"] for r in mon] == [1, 10]
    assert all(r["crossings_err_limit"] <= MONITOR and r["run_crossings"] == 10 for r in mon)


def test_bar_runs_name_the_geometry():
    from spheral_frag import naming
    import spheral_m2_check5 as c5
    b = naming.canonical_brackets({"strength": "constant:100e6", "sampler": "uniform", "loads": "off",
                                   "geometry": "bar"})
    name = naming.run_name("prep_x_k00000-00010", "replay", (0, 10), "rz", 1.1, b, 1)
    assert name.endswith("_strength-constant-100000000_sampler-uniform_loads-off_geometry-bar")
    assert naming.parse_run_name(name)["brackets"] == b
    assert naming.BRACKET_DEFAULTS["geometry"] == "frame"
    d = naming.run_dir("runs", c5.run_name("3d", 1.1))         # nested under the prepared name (Task 7b)
    assert all(len(p) <= naming.NAME_MAX for p in d.split(os.sep))
    assert naming.run_name_from_dir(d) == c5.run_name("3d", 1.1)


@pytest.mark.spheral
def test_rz_11_rows_reproduce():
    """RZ 1.1 mm over frames 0-3 on 4 processes (the study's decomposition): frames 1-3 equal the committed rows
    (a monitored run's accepted states equal an unmonitored run's, and both monitor interval 0 -> 1), the holds are
    exact, and a bar with --loads on is refused."""
    import spheral_m2_check5 as c5
    from spheral_frag import material, naming
    from spheral_frag.runner import run
    prep = os.path.join(REPO_ROOT, c5.PREPARED)
    body = os.path.join(REPO_ROOT, c5.body_dir("rz", 1.1))
    if not (os.path.isfile(os.path.join(prep, "prepare.json")) and os.path.isfile(os.path.join(body, "body.json"))):
        pytest.skip("the MVP flight's prepared input or the RZ 1.1 mm bar body is not built (git-ignored; "
                    "analysis/spheral_m2_check5.py builds it)")
    out = os.path.join(REPO_ROOT, "spheral_output", "tests", "check5_rz")
    shutil.rmtree(out, ignore_errors=True)
    args = ["spheral_frag/runner/run.py", "replay", "--body", c5.body_dir("rz", 1.1), "--prepared", c5.PREPARED,
            "--frames", "0-3", "--strength", c5.STRENGTH, "--ke-fraction", "1e300"]
    r = subprocess.run([spheral_launcher(), "-n", "4"] + args + ["--loads", "off", "--out",
                                                                 os.path.relpath(out, REPO_ROOT)],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    bad = subprocess.run([spheral_launcher()] + args + ["--loads", "on", "--out", "spheral_output/tests/x"],
                         cwd=REPO_ROOT, capture_output=True, text=True, timeout=600)
    assert bad.returncode == 2 and "--loads off" in bad.stderr
    table = material.MaterialTable.load(prep)
    E_, nu, alpha = c5.constants(table)
    b = naming.canonical_brackets({"strength": c5.STRENGTH, "sampler": "uniform", "loads": "off", "geometry": "bar"})
    d = naming.run_dir(out, naming.run_name(os.path.basename(c5.PREPARED), "replay", (0, 3), "rz", 1.1, b, 1))
    iv = {int(x["k1"]): x for x in run.read_intervals(os.path.join(d, "intervals.csv"))}
    for k in (1, 2, 3):
        m = c5.measure(d, k, table, E_, alpha, 1.1e-3, "rz")
        row, = pick(rows(), form="rz", dx_mm=1.1, k=k)
        assert m["sxx_mean"] == pytest.approx(row["sxx_mean_Pa"], rel=1e-9)
        assert m["ps_mean"] == pytest.approx(row["ps_mean"], rel=1e-9, abs=1e-15)
        assert m["lat_max"]/Y == pytest.approx(row["lat_max_over_Y"], rel=1e-9)
        assert float(iv[k]["held_max_disp_m"]) == 0.0 and float(iv[k]["held_max_speed_ms"]) == 0.0
    shutil.rmtree(out, ignore_errors=True)
