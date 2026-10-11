"""Check 6 (Spheral M2 plan, Task 8; spec §12 check 6): frame 182's loads on the elastic 100 mm sphere at 300 K, the
deceleration, and the finite-element reference (analysis/spheral_m2_elastic_fe.py).

The normal-loop tests read the committed table `data/spheral/m2_check6.csv` (written by
`analysis/spheral_m2_check6.py`, whose docstring defines every measure): eight runs -- RZ 2.2 and 1.1 mm, 3D 3.0 and
2.2 mm, each with the default area shares ("volume": V^(2/3), one calibration factor) and with the frame's patch
areas ("patch") -- over frames 0-6 (the table ramped in over interval 0 -> 1, held for five more), one row per run
and frame k = 1 ... 6. VM0 = 3,510.5 Pa, the reference's largest von Mises stress; stress errors are over VM0.
The `spheral` test reruns RZ 2.2 mm over frames 0-2 with both shares and requires the committed rows back.

THRESHOLDS: measured 2026-10-10/11 (the CSV), each the worst case over frames 1-6 of its run (or over the runs where
one value serves several), rounded up to one significant figure; worst case in brackets. Change them only together
with the CSV, the README and docs/spheral_frag_assumptions.md (CLAUDE.md)."""
import csv
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))

CSV = os.path.join(REPO_ROOT, "data", "spheral", "m2_check6.csv")
DRAG_REF = 30.716800986256963          # N: frame 182's table (with the lee extension) over the exact sphere
FE_DECEL = 20.8689                     # m/s^2: the finite-element reference's a_x (1.667 mm mesh)
HISTORY_DECEL = 86.687                 # m/s^2: the history at k 182 (a diagnostic: its 0.572 kg body)
RUNS = [(sh, f, dx) for sh in ("volume", "patch") for f, dx in (("rz", 2.2), ("rz", 1.1), ("3d", 3.0), ("3d", 2.2))]

# --- thresholds, per (shares, form, dx mm); stress errors over VM0
INT_VM_MAX = {("volume", "rz", 2.2): 0.7,    # [0.622] interior (|x_ref| <= 45 mm): max vM(sigma_p - sigma_FE)
              ("volume", "rz", 1.1): 0.4,    # [0.347]
              ("volume", "3d", 3.0): 2.0,    # [1.199] (the interior reaches 1.7 spacings from the surface)
              ("volume", "3d", 2.2): 0.6,    # [0.577]
              ("patch", "rz", 2.2): 0.3,     # [0.284]
              ("patch", "rz", 1.1): 0.2,     # [0.187]
              ("patch", "3d", 3.0): 0.5,     # [0.460]
              ("patch", "3d", 2.2): 0.7}     # [0.646]
INT_VM_RMS = {("volume", "rz", 2.2): 0.2,    # [0.102] interior: volume-weighted RMS
              ("volume", "rz", 1.1): 0.09,   # [0.0883]
              ("volume", "3d", 3.0): 0.2,    # [0.128]
              ("volume", "3d", 2.2): 0.07,   # [0.0678]
              ("patch", "rz", 2.2): 0.06,    # [0.0517]
              ("patch", "rz", 1.1): 0.05,    # [0.0430]
              ("patch", "3d", 3.0): 0.06,    # [0.0556]
              ("patch", "3d", 2.2): 0.04}    # [0.0395]
CORE_VM_MAX = {("volume", "rz", 2.2): 0.3,   # [0.212] core (|x_ref| <= 25 mm): max vM of the difference
               ("volume", "rz", 1.1): 0.3,   # [0.218]
               ("volume", "3d", 3.0): 0.2,   # [0.129]
               ("volume", "3d", 2.2): 0.07,  # [0.0662]
               ("patch", "rz", 2.2): 0.06,   # [0.0597]
               ("patch", "rz", 1.1): 0.08,   # [0.0703]
               ("patch", "3d", 3.0): 0.09,   # [0.0839]
               ("patch", "3d", 2.2): 0.04}   # [0.0400]
CORE_P_MEAN_3D = 4e-3                  # [3.4e-3, patch 3D 2.2] 3D: |mean pressure offset| of the core
CORE_P_MEAN_RZ = (0.01, 0.05)          # [-0.0153 .. -0.0493] RZ: the core more compressed (SolidSPHRZ's section)
SHELL_VM_MAX = {"volume": 2.0,         # [1.70, RZ 1.1] shell (|x_ref| > 45 mm): max vM of the difference
                "patch": 3.0}          # [2.99, 3D 2.2] (zero-area inner corners next to loaded particles)
SHELL_VM_RMS = 0.4                     # [0.314, volume 3D 3.0] shell: RMS
DRAG_REL = {("volume", "rz", 2.2): 3e-2,     # [-0.0295] drag (= deceleration x M) against the table's integral
            ("volume", "rz", 1.1): 4e-2,     # [-0.0302]
            ("volume", "3d", 3.0): 7e-3,     # [-0.0063]
            ("volume", "3d", 2.2): 2e-3,     # [-0.0017]
            ("patch", "rz", 2.2): 4e-3,      # [+0.0035]
            ("patch", "rz", 1.1): 3e-3,      # [+0.0020]
            ("patch", "3d", 3.0): 6e-3,      # [-0.0057]
            ("patch", "3d", 2.2): 8e-4}      # [+0.0008]
RZ_VOLUME_DRAG_DEFICIT = 0.02          # the volume shares' RZ drag deficit is real: below -2 % at both spacings
NORMAL_MAX_DEG = {"rz": 6.0,           # [5.74] surface normals' largest angle from radial
                  "3d": 3.0}           # [2.39]
SECTION_RZ = (1.04, 1.08)              # [1.0500 .. 1.077] RZ mid-plane axial force / the load it must carry
SECTION_3D = (0.97, 0.99)              # [0.974 .. 0.984] 3D
MOMENTUM_3D = 2e-16                    # [1.3e-16] net momentum after an interval / its drag impulse
MOMENTUM_RZ = 2e-3                     # [1.9e-3] (SolidSPHRZ does not conserve axial momentum)
RECENTRE_RZ_DX = 4e-7                  # [3.5e-7] RZ: the axial drift removed at a boundary, in spacings
KE_RAMP = 7e-4                         # [6.5e-4] internal KE / elastic energy at k 1 (the end of the ramp)
KE_HELD = 2e-5                         # [1.5e-5] k >= 2 (held), every run: far below the plan's 1 %
SETTLE = 3e-4                          # [2.0e-4] k 6: max vM(sigma_6 - sigma_5) over VM0
SETTLE_P = 2e-4                        # [1.2e-4] k 6: the same for the pressure
MONITOR_PATCH = 2e-2                   # [0.0150] the monitor's estimate on interval 0 -> 1 (no switch)
MONITOR_VOLUME_RZ = 3e-2               # [0.0231]
MONITOR_VOLUME_3D = (0.05, 0.1)        # [0.065, 0.090] above the 5 % limit: switched to 20 crossings


def rows():
    with open(CSV, newline="") as fh:
        out = []
        for r in csv.DictReader(fh):
            for c in r:
                if c not in ("form", "area_shares", "run"):
                    r[c] = float(r[c])
            for c in ("k", "n_particles", "n_interior", "n_shell", "processes", "n_surface", "n_exposed",
                      "n_zero_area", "monitored", "run_crossings", "steps"):
                r[c] = int(r[c])
            out.append(r)
    return out


def run_rows(key):
    sh, f, dx = key
    return sorted([r for r in rows() if (r["area_shares"], r["form"], r["dx_mm"]) == key], key=lambda r: r["k"])


def test_table_is_complete():
    import spheral_m2_check6 as c6
    rs = rows()
    assert {(r["area_shares"], r["form"], r["dx_mm"], r["k"]) for r in rs} == \
        {(sh, f, dx, k) for sh, f, dx in RUNS for k in range(1, 7)}
    assert c6.FRAMES == (0, 6) and c6.COMPARE_K == 6 and c6.LOADS_FRAME == 182
    for r in rs:
        assert r["drag_ref_N"] == pytest.approx(DRAG_REF, rel=1e-12)
        assert r["fe_decel_ms2"] == pytest.approx(FE_DECEL, abs=1e-4)
        assert r["flight_decel_ms2"] == pytest.approx(HISTORY_DECEL, abs=1e-3)
        assert r["mass_kg"] == pytest.approx(1.47202, rel=1e-5)     # frame 0's finite-element solid mass


@pytest.mark.parametrize("key", RUNS)
def test_settles(key):
    """The internal KE falls below 1 % at once (the ramp's end) and to 1e-5 held; the stress stops changing."""
    rs = run_rows(key)
    assert rs[0]["ke_internal_ratio"] <= KE_RAMP
    assert max(r["ke_internal_ratio"] for r in rs[1:]) <= KE_HELD
    assert rs[-1]["settle_vm"] <= SETTLE and rs[-1]["settle_p"] <= SETTLE_P
    assert rs[-1]["settle_vm"] < rs[1]["settle_vm"]/100.0


@pytest.mark.parametrize("key", RUNS)
def test_body_stays_at_rest(key):
    for r in run_rows(key):
        if key[1] == "3d":
            assert abs(r["momentum_over_impulse"]) <= MOMENTUM_3D
        else:
            assert abs(r["momentum_over_impulse"]) <= MOMENTUM_RZ
            assert abs(r["recentre_shift_dx"]) <= RECENTRE_RZ_DX
        assert r["max_disp_ref_dx"] < 1e-5                     # the elastic displacements are nanometres


@pytest.mark.parametrize("key", RUNS)
def test_deceleration_and_drag(key):
    """sum F / M is the deceleration (exactly), within the shares' drag error of the table's integral over the exact
    sphere; the history's deceleration is a different body (diagnostic only, 4.1-4.3 times larger)."""
    for r in run_rows(key):
        assert r["decel_ms2"]*r["mass_kg"] == pytest.approx(r["drag_N"], rel=1e-12)
        assert r["decel_ref_ms2"] == pytest.approx(DRAG_REF/r["mass_kg"], rel=1e-12)
        assert abs(r["drag_rel"]) <= DRAG_REL[key]
        assert abs(r["decel_ms2"]/FE_DECEL - 1.0) <= DRAG_REL[key] + 1e-3
        assert 4.0 < r["flight_decel_ms2"]/r["decel_ms2"] < 4.4
    if key[0] == "volume" and key[1] == "rz":
        assert all(r["drag_rel"] < -RZ_VOLUME_DRAG_DEFICIT for r in run_rows(key))


@pytest.mark.parametrize("key", RUNS)
def test_surface_detection(key):
    """Every surface particle of the convex sphere is exposed; the shares cover frame 0's faceted area (by
    construction); the normals lie within a few degrees of radial; the normals' error costs at most 0.6 % of drag."""
    sh, f, dx = key
    for r in run_rows(key):
        assert r["n_exposed"] == r["n_surface"]
        assert r["area_cover"] == pytest.approx(0.99967, abs=1e-5)
        assert r["normal_err_max_deg"] <= NORMAL_MAX_DEG[f]
        assert abs(r["drag_surface_N"]/r["drag_radial_normals_N"] - 1.0) <= 6e-3     # [5.7e-3, 3D 3.0 patch]
        if sh == "volume":
            assert r["n_zero_area"] == 0
    n = {("rz", 2.2): 63, ("rz", 1.1): 127, ("3d", 3.0): 3350, ("3d", 2.2): 6194}[(f, dx)]
    assert all(r["n_surface"] == n for r in run_rows(key))


@pytest.mark.parametrize("key", RUNS)
def test_stresses_against_the_finite_element_reference(key):
    sh, f, dx = key
    for r in run_rows(key):
        assert r["int_vm_max"] <= INT_VM_MAX[key]
        assert r["int_vm_rms"] <= INT_VM_RMS[key]
        assert r["core_vm_max"] <= CORE_VM_MAX[key]
        assert r["shell_vm_max"] <= SHELL_VM_MAX[sh]
        assert r["shell_vm_rms"] <= SHELL_VM_RMS
        if f == "3d":
            assert abs(r["core_p_mean"]) <= CORE_P_MEAN_3D
            assert SECTION_3D[0] <= r["section_ratio"] <= SECTION_3D[1]
        else:
            assert CORE_P_MEAN_RZ[0] <= -r["core_p_mean"] <= CORE_P_MEAN_RZ[1]
            assert SECTION_RZ[0] <= r["section_ratio"] <= SECTION_RZ[1]
    last = run_rows(key)[-1]
    # the reference resolves far finer than the comparison: its own mesh difference at the particles
    assert last["fe_conv_int_vm_max"] < 0.01 and last["fe_conv_shell_vm_max"] < 0.05


def test_patch_shares_improve_the_rz_load_and_every_core():
    """The patch shares remove the volume shares' RZ drag deficit and lower the interior RMS error in every run."""
    for f, dx in (("rz", 2.2), ("rz", 1.1), ("3d", 3.0), ("3d", 2.2)):
        v, p = run_rows(("volume", f, dx))[-1], run_rows(("patch", f, dx))[-1]
        assert p["int_vm_rms"] < v["int_vm_rms"]
        if f == "rz":
            assert abs(p["drag_rel"]) < abs(v["drag_rel"])/5.0
            assert p["core_vm_max"] < v["core_vm_max"]/2.0


@pytest.mark.parametrize("key", RUNS)
def test_monitor(key):
    sh, f, dx = key
    rs = run_rows(key)
    mon = [r for r in rs if r["monitored"]]
    assert [r["k"] for r in mon] == [1]
    est = mon[0]["crossings_err_limit"]
    if sh == "patch":
        assert est <= MONITOR_PATCH and rs[-1]["run_crossings"] == 10
    elif f == "rz":
        assert est <= MONITOR_VOLUME_RZ and rs[-1]["run_crossings"] == 10
    else:
        assert MONITOR_VOLUME_3D[0] < est <= MONITOR_VOLUME_3D[1] and rs[-1]["run_crossings"] == 20


def test_patch_shares_assign_each_patch_to_the_nearest_surface_particle():
    from spheral_frag.runner.packages import loads as pl
    xs = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]])
    c = np.array([[0.1, 0.0, 0.9], [0.0, 0.1, -0.8], [0.9, 0.1, 0.0], [0.0, 0.0, 0.8]])
    a = np.array([1.0, 2.0, 4.0, 8.0])
    assert np.array_equal(pl.patch_shares("3d", xs, c, a), [9.0, 2.0, 4.0])
    # RZ: the ring (x, r) collects every patch around the axis nearest it
    xr = np.array([[0.0, 1.0], [1.0, 0.1]])
    c = np.array([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [0.0, 0.6, 0.8], [1.0, 0.0, 0.1]])
    assert np.array_equal(pl.patch_shares("rz", xr, c, a), [7.0, 8.0])
    assert pl.AREA_SHARES == ("volume", "patch")


def test_run_names_carry_the_held_frame_and_the_shares():
    from spheral_frag import naming
    import spheral_m2_check6 as c6
    name = c6.run_name("3d", 2.2, "volume")
    assert name.endswith("_strength-elastic_sampler-uniformdt0_loadsframe-182_areashares-volume")
    assert c6.run_name("rz", 1.1, "patch").endswith("_sampler-uniformdt0_loadsframe-182")     # decision 9's default
    b = naming.parse_run_name(name)["brackets"]
    assert b["loads_frame"] == 182 and b["area_shares"] == "volume"
    assert naming.BRACKET_DEFAULTS["loads_frame"] == -1 and naming.BRACKET_DEFAULTS["area_shares"] == "patch"
    d = c6.run_dir("3d", 2.2, "volume")
    assert all(len(p) <= naming.NAME_MAX for p in d.split(os.sep))
    assert naming.run_name_from_dir(d) == name


@pytest.mark.spheral
def test_rz_22_rows_reproduce():
    """RZ 2.2 mm over frames 0-2 with both shares on 1 process (the study's decomposition): frames 1-2 equal the
    committed rows, and --loads-frame and --area-shares need --loads on."""
    import spheral_m2_check6 as c6
    from spheral_frag import material
    prep = os.path.join(REPO_ROOT, c6.PREPARED)
    body = os.path.join(REPO_ROOT, c6.body_dir("rz", 2.2))
    if not all(os.path.isfile(p) for p in (os.path.join(prep, "prepare.json"), os.path.join(body, "body.json"),
                                            os.path.join(REPO_ROOT, c6.FE_FINE))):
        pytest.skip("the MVP flight's prepared input, the RZ 2.2 mm body or the finite-element solution is not "
                    "built (git-ignored; analysis/spheral_m2_check6.py and spheral_m2_elastic_fe.py build them)")
    out = os.path.join(REPO_ROOT, "spheral_output", "tests", "check6_rz")
    shutil.rmtree(out, ignore_errors=True)
    rel = os.path.relpath(out, REPO_ROOT)
    for sh in ("volume", "patch"):
        r = subprocess.run([spheral_launcher()] + c6.run_command("rz", 2.2, sh, rel, frames=(0, 2)),
                           cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    for bad in (["--loads", "off", "--loads-frame", "182"], ["--loads", "off", "--area-shares", "volume"]):
        b = subprocess.run([spheral_launcher(), "spheral_frag/runner/run.py", "replay", "--body", c6.body_dir("rz", 2.2),
                            "--prepared", c6.PREPARED, "--frames", "0-2", "--strength", "elastic",
                            "--out", "spheral_output/tests/x"] + bad,
                           cwd=REPO_ROOT, capture_output=True, text=True, timeout=600)
        assert b.returncode == 2 and "--loads on" in b.stderr
    table = material.MaterialTable.load(prep)
    fine = c6.reference(c6.FE_FINE)
    for sh in ("volume", "patch"):
        d = c6.run_dir("rz", 2.2, sh, out=out, frames=(0, 2))
        _, prev = c6.state(d, 0, table, "rz")
        for k in (1, 2):
            m = c6.measure(d, k, table, "rz", fine, None, prev)
            prev = m.pop("_state")
            row, = [x for x in run_rows((sh, "rz", 2.2)) if x["k"] == k]
            for c in ("int_vm_max", "int_vm_rms", "core_vm_max", "shell_vm_max", "settle_vm"):
                assert m[c] == pytest.approx(row[c], rel=1e-9), (sh, k, c)
            h = np.load(os.path.join(d, "history", "interval_{:05d}.npz".format(k - 1)))
            assert -float(h["history"][-1, 1]) == pytest.approx(row["drag_N"], rel=1e-12)
    shutil.rmtree(out, ignore_errors=True)
