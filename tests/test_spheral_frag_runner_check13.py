"""Check 13 (Spheral M2 plan, Task 10; spec §12 check 13): RZ against 3D on the MVP flight's frames 0-47, the replay
at every default (flight strength, loads and deceleration, patch area shares, SolidSPH, 10 crossings with the
crossings monitor at its defaults).

The normal-loop tests read the committed tables written by `analysis/spheral_m2_check13.py` (whose docstring defines
every measure and records the findings): `data/spheral/m2_check13.csv` (the comparisons at k 12, 24, 36, 47, per
pair and region, over the 3D run's largest von Mises at k), `m2_check13_runs.csv` (one row per run and interval:
the monitor, the cost, the core's axial stress and the pairing amplitude) and `m2_check13_probe.csv` (the interval
reruns at more crossings). The `spheral` test reruns RZ 2.2 mm over frames 0-2 and requires the study's checkpoints
back bitwise.

RZ 1.1 mm is not held to the RZ-against-3D thresholds: from k 5 SolidSPHRZ's tensile pairing grows in its core
(under the core's 40-50 MPa axial tension) and takes the core's axial stress down to a quarter of the 3D's by k 47;
the tests hold that finding instead.

THRESHOLDS: measured 2026-10-11 (the CSVs), each the worst case over k 12, 24, 36, 47 (or over the intervals),
rounded up to one significant figure; worst case in brackets. Change them only together with the CSVs, the README
and docs/spheral_frag_assumptions.md (CLAUDE.md)."""
import csv
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))

DATA = os.path.join(REPO_ROOT, "data", "spheral")
CSV = os.path.join(DATA, "m2_check13.csv")
RUNS_CSV = os.path.join(DATA, "m2_check13_runs.csv")
PROBE_CSV = os.path.join(DATA, "m2_check13_probe.csv")
COMPARE_K = (12, 24, 36, 47)
PAIRS = [("rz2.2", "3d2.2"), ("rz1.1", "3d2.2"), ("rz2.2", "rz1.1")]
REGIONS = ("axis", "axis_interior", "surface", "subsurface", "shell", "interior", "core", "all")
RUNS = [("rz", 2.2), ("rz", 1.1), ("3d", 2.2)]

# --- RZ 2.2 mm against 3D 2.2 mm: per region (vm, p) maxima of |difference| over the 3D peak von Mises at k
RZ_3D_MAX = {"core": (0.2, 0.07),              # [0.128 k 36, 0.0659 k 36]
             "axis_interior": (0.2, 0.3),      # [0.118 k 36, 0.296 k 36: the axis row at the nose's eutectic]
             "axis": (0.2, 0.4),               # [0.118, 0.328 k 36]
             "interior": (0.2, 0.3),           # [0.162 k 12, 0.296 k 36]
             "subsurface": (0.2, 0.4),         # [0.196 k 36, 0.328 k 36]
             "surface": (0.4, 0.3)}            # [0.332 k 24, 0.283 k 36] the outer spacing: both forms' free surface
RZ_3D_RMS = {"core": (0.06, 0.03),             # [0.0590, 0.0252] volume-weighted RMS
             "interior": (0.04, 0.04),         # [0.0390, 0.0329]
             "all": (0.07, 0.06)}              # [0.0607, 0.0521]
CORE_P_MEAN = (-0.02, 0.0)                     # [-0.0080 .. -0.0126] RZ's core slightly more tensile than 3D's
CORE_SXX_RATIO = (1.04, 1.05)                  # [1.0405 .. 1.0465] RZ's core axial stress over 3D's (check 6: RZ
                                               # carries 5-8 % more axial force, 3D 2 % less)
FLOOR = (0.06, 0.2)                            # [0.067 .. 0.113] the body stage's azimuth spread as a stress
# --- RZ 1.1 mm's core (the pairing): from k 9 on
RZ11_PAIRING = 0.05                            # [0.0672 at k 9, 0.104 at k 47] at least this, every k >= 9
RZ11_CORE_SXX_RATIO = 0.4                      # [0.24 .. 0.34] at most this, RZ 1.1's core axial stress over 3D's
STABLE_PAIRING = 0.02                          # [0.0128 RZ 2.2 k 44, 0.0145 3D 2.2 k 42] at most this at every k
# --- runs
KE_INTERNAL = 1e-3                             # [9.99e-4] internal KE / elastic energy, every interval
MONITORED_K0 = (0, 9, 19, 29, 39)
SWITCHES = {("rz", 2.2): [(29, 20), (39, 40)],     # (k0, to crossings)
            ("rz", 1.1): [(9, 20), (29, 40)],
            ("3d", 2.2): [(39, 20)]}
STEPS_REL = 0.02                               # [+1.32 %, RZ 1.1 at 40] steps per interval against the gate's x N_c/10
MS_PER_STEP = {("rz", 2.2): 5.0,               # [4.51] ms/step, every undisturbed interval (gate: 3.0 built-in, 5.4
               ("rz", 1.1): 9.0,               # [8.35]  with the Python classes; RZ 1.1 on 4 processes, gate on 6:
               ("3d", 2.2): 200.0}             # [137.0] 4.2 / 7.5; 3D gate 113 / 135)


def rows(path=CSV):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def pick(target, source, region, k=None):
    return [r for r in rows() if r["target"] == target and r["source"] == source and r["region"] == region
            and (k is None or int(r["k"]) == k)]


def run_rows(form, dx):
    return [r for r in rows(RUNS_CSV) if r["form"] == form and float(r["dx_mm"]) == dx]


def test_tables_are_complete():
    r = rows()
    assert len(r) == len(PAIRS)*len(COMPARE_K)*len(REGIONS)
    for t, s in PAIRS:
        for k in COMPARE_K:
            assert {x["region"] for x in r if (x["target"], x["source"], int(x["k"])) == (t, s, k)} == set(REGIONS)
    for form, dx in RUNS:
        rr = run_rows(form, dx)
        assert [int(x["k1"]) for x in rr] == list(range(1, 48))
        assert all(x["run"].endswith("__replay_k00000-00047_{}_dx{:.2f}mm_seed1_hourglass-0".format(form, dx))
                   for x in rr)
    assert len(rows(PROBE_CSV)) > 0


@pytest.mark.parametrize("region", sorted(RZ_3D_MAX))
def test_rz_agrees_with_3d(region):
    vm, p = RZ_3D_MAX[region]
    for r in pick("rz2.2", "3d2.2", region):
        assert float(r["vm_max"]) <= vm and float(r["p_max"]) <= p, (region, r["k"], r["vm_max"], r["p_max"])
    if region in RZ_3D_RMS:
        vm, p = RZ_3D_RMS[region]
        for r in pick("rz2.2", "3d2.2", region):
            assert float(r["vm_rms"]) <= vm and float(r["p_rms"]) <= p, (region, r["k"])


def test_the_core_agrees_within_the_azimuth_floor():
    """The core's largest pressure difference and both RMS differences lie below the body stage's azimuth spread as
    a stress (E alpha dT_azimuth / (1 - nu) over the peak): the floor of what RZ and 3D can agree to."""
    for r in pick("rz2.2", "3d2.2", "core"):
        floor = float(r["floor_stress"])
        assert FLOOR[0] <= floor <= FLOOR[1]
        assert float(r["p_max"]) < floor and float(r["vm_rms"]) < floor and float(r["p_rms"]) < floor, r["k"]


def test_rz_offset_in_the_core():
    """SolidSPHRZ's core carries 4-5 % more axial stress than 3D's (check 6's section excess, here under thermal
    stress), and its mean pressure is 1 % of the peak lower (more tensile; check 4: RZ +1.4 %, 3D -2 %); the loads'
    own share of it is negligible (check 6's offset scaled to these loads is below 2e-6 of the peak)."""
    for r in pick("rz2.2", "3d2.2", "core"):
        assert CORE_P_MEAN[0] <= float(r["p_mean"]) <= CORE_P_MEAN[1], r["k"]
        assert CORE_SXX_RATIO[0] <= float(r["core_sxx_ratio"]) <= CORE_SXX_RATIO[1], r["k"]
        assert float(r["rz_load_offset"]) < 2e-6


def test_rz_11_core_pairs():
    """RZ 1.1 mm: the tensile pairing (a zigzag of the axial displacement along the rows) grows from k 5 and its core
    loses most of its axial tension; RZ 2.2 mm and 3D 2.2 mm stay below 0.02 spacings."""
    p11 = {int(r["k1"]): float(r["pairing_dx"]) for r in run_rows("rz", 1.1)}
    assert p11[3] < 0.005 and all(p11[k] >= RZ11_PAIRING for k in range(9, 48))
    for key in (("rz", 2.2), ("3d", 2.2)):
        assert max(float(r["pairing_dx"]) for r in run_rows(*key)) <= STABLE_PAIRING, key
    for k in COMPARE_K:
        r, = pick("rz1.1", "3d2.2", "core", k)
        assert float(r["core_sxx_ratio"]) <= RZ11_CORE_SXX_RATIO, k


def test_rz_11_is_off_where_it_pairs():
    """Against 3D and against RZ 2.2 mm, RZ 1.1 mm's core pressure is off by far more than RZ 2.2 mm's (flagged,
    not held to the thresholds)."""
    for k in COMPARE_K:
        good, = pick("rz2.2", "3d2.2", "core", k)
        bad, = pick("rz1.1", "3d2.2", "core", k)
        assert float(bad["p_mean"]) > 10*abs(float(good["p_mean"])), k


@pytest.mark.parametrize("key", RUNS)
def test_runs_are_quiet_and_complete(key):
    rr = run_rows(*key)
    assert max(float(r["ke_internal_ratio"]) for r in rr) <= KE_INTERNAL
    assert all(int(r["repeats"]) == 0 for r in rr)
    assert all(int(r["n_particles"]) == int(rr[0]["n_particles"]) for r in rr)      # nothing removed before k 48


@pytest.mark.parametrize("key", RUNS)
def test_monitor(key):
    """Monitored on the first interval and every 10th (j = 9, 19, ...); a switch exactly where the larger of the von
    Mises and pressure estimates exceeds 5 %, to twice the crossings."""
    rr = run_rows(*key)
    mon = [r for r in rr if int(r["monitored"])]
    assert tuple(int(r["k0"]) for r in mon) == MONITORED_K0
    sw = [(int(r["k0"]), int(r["run_crossings"])) for r in rr if int(r["crossings_switch"])]
    assert sw == SWITCHES[key]
    for r in mon:
        assert (float(r["crossings_err_limit"]) > 0.05) == bool(int(r["crossings_switch"])), r["k0"]
        assert float(r["crossings_err_limit"]) == max(float(r["crossings_err_est"]), float(r["crossings_err_pressure"]))
    crossings = 10
    for r in rr:
        if int(r["crossings_switch"]):
            crossings *= 2
        assert int(r["run_crossings"]) == crossings


@pytest.mark.parametrize("key", RUNS)
def test_cost_against_the_gate(key):
    rr = [r for r in run_rows(*key) if not int(r["disturbed"])]
    at10 = [r for r in rr if int(r["crossings"]) == 10]
    gate = int(at10[0]["gate_steps"])
    assert all(abs(int(r["steps"])/gate - 1.0) <= STEPS_REL for r in at10)
    assert max(float(r["wall_per_step_ms"]) for r in rr) <= MS_PER_STEP[key]
    for r in rr:                                  # steps scale with the crossings
        assert abs(int(r["steps"])/(gate*int(r["crossings"])/10) - 1.0) <= STEPS_REL


def test_probes_converge_with_crossings():
    """The reruns at more crossings converge: at every probed interval the error against the most crossings falls
    with N, and the monitor's 2 |sigma_N - sigma_2N| bounds it within a factor 2."""
    pr = rows(PROBE_CSV)
    for name in sorted({r["probe"] for r in pr if r["probe"].startswith("monitor_")}):
        q = sorted((r for r in pr if r["probe"] == name and r["err_ref_vm"] not in ("", "nan")),
                   key=lambda r: int(r["crossings"]))
        errs = [max(float(r["err_ref_vm"]), float(r["err_ref_p"])) for r in q]
        assert all(b < a for a, b in zip(errs, errs[1:])), (name, errs)
        for r in q:
            if r["limit_2n"] not in ("", "nan"):
                assert float(r["limit_2n"]) <= 2.0*max(float(r["err_ref_vm"]), float(r["err_ref_p"])) + 0.03, name


def test_pairing_grows_with_crossings():
    """RZ 1.1 mm, interval 4 -> 5 from the run's checkpoint: the pairing amplitude grows by more per interval at more
    crossings (it grows per step, not per unit of flight time), so the monitor's switches feed it."""
    q = sorted((r for r in rows(PROBE_CSV) if r["probe"] == "pairing_rz11_k4"), key=lambda r: int(r["crossings"]))
    ends = [float(r["pairing_end_dx"]) for r in q]
    assert all(b > 2.0*a for a, b in zip(ends, ends[1:]))
    assert ends[0] > 2.0*float(q[0]["pairing_start_dx"])


def test_names_and_paths():
    import spheral_m2_check13 as c13
    from spheral_frag import naming
    for form, dx, _ in c13.FORMS:
        name = c13.run_name(form, dx)
        assert naming.parse_run_name(name)["brackets"] == naming.canonical_brackets({"hourglass": c13.HOURGLASS})
        assert name.endswith("_hourglass-0")              # Task 10's runs predate Task 10b's hourglass default
        d = c13.run_dir(form, dx)
        assert all(len(p) <= naming.NAME_MAX for p in d.split(os.sep))
        assert naming.run_name_from_dir(d) == name


def test_cylindrical_components():
    """The 3D stress resolved into (axial, radial, hoop) at each particle's azimuth recovers an axisymmetric field."""
    import spheral_m2_check13 as c13
    rng = np.random.default_rng(3)
    phi = rng.uniform(0.0, 2.0*np.pi, 64)
    rho = rng.uniform(1e-3, 4e-2, 64)
    xr = np.stack([rng.uniform(-0.04, 0.04, 64), rho*np.cos(phi), rho*np.sin(phi)], axis=1)
    A = np.array([[1.0, 0.3, 0.0], [0.3, 2.0, 0.0], [0.0, 0.0, 3.0]])
    Q = np.zeros((64, 3, 3))
    Q[:, 0, 0] = 1.0
    Q[:, 1, 1], Q[:, 1, 2], Q[:, 2, 1], Q[:, 2, 2] = np.cos(phi), np.sin(phi), -np.sin(phi), np.cos(phi)
    sigma = np.einsum("nji,jk,nkl->nil", Q, A, Q)
    x, r, s = c13.cylindrical(xr, sigma)
    assert np.allclose(s, A[None], atol=1e-14) and np.allclose(r, rho, rtol=1e-14)


@pytest.mark.spheral
def test_rz_22_reproduces():
    """RZ 2.2 mm over frames 0-2 on 1 process (the study's decomposition) gives the study's checkpoints 1 and 2
    bitwise, including the monitored interval 0 -> 1."""
    import spheral_m2_check13 as c13
    from spheral_frag.runner import checkpoint
    prep = os.path.join(REPO_ROOT, c13.PREPARED)
    body = os.path.join(REPO_ROOT, c13.body_dir("rz", 2.2))
    study = c13.run_dir("rz", 2.2)
    if not all(os.path.isfile(p) for p in (os.path.join(prep, "prepare.json"), os.path.join(body, "body.json"),
                                            os.path.join(study, "summary.json"))):
        pytest.skip("the MVP flight's prepared input, the RZ 2.2 mm body or the study's run is not there "
                    "(git-ignored; analysis/spheral_m2_check13.py builds them)")
    out = os.path.join(REPO_ROOT, "spheral_output", "tests", "check13_rz")
    shutil.rmtree(out, ignore_errors=True)
    r = subprocess.run([spheral_launcher(), "spheral_frag/runner/run.py", "replay", "--body", c13.body_dir("rz", 2.2),
                        "--prepared", c13.PREPARED, "--frames", "0-2", "--hourglass", "0", "--out",
                        os.path.relpath(out, REPO_ROOT)],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    from spheral_frag import naming
    d = naming.run_dir(out, naming.run_name(os.path.basename(c13.PREPARED), "replay", (0, 2), "rz", 2.2,
                                            naming.canonical_brackets({"hourglass": 0.0}), 1))
    for k in (1, 2):
        a = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
        b = checkpoint.read(os.path.join(study, "checkpoints", checkpoint.filename(k)))
        for key in ("pid", "x", "v", "S", "S_zz", "rho", "u", "plastic_strain"):
            assert np.array_equal(a[key], b[key]), (k, key)
    shutil.rmtree(out, ignore_errors=True)
