"""Check 7 (Spheral M2 plan, Task 9; spec §12 check 7, §14.1): tension behaviour and the choice between SolidSPH and
SolidCRKSPH.

The normal-loop tests read the committed tables `data/m2_check7.csv` and `data/m2_check7_damage.csv`
(written by `analysis/spheral_m2_check7.py`, whose docstring defines every measure): checks 4, 5 and 6 at one or two
spacings per form with both schemes (the SPH rows measured from checks 4-6's own runs, the CRKSPH rows from this
study's), plus the free uniform expansion ("fx", 50 K over one interval, the patch test), one row per run and frame;
the damage probe one row per scheme and form. The `spheral` test reruns check 6's RZ 2.2 mm CRKSPH case over frames
0-2, stopped after frame 1 and resumed, and requires the committed rows back.

THRESHOLDS: measured 2026-10-11 (the CSVs), each the worst case over the run's frames (or over the runs one value
serves), rounded up to one significant figure; worst case in brackets. Change them only together with the CSVs,
the README and docs/spheral_frag_assumptions.md (CLAUDE.md)."""
import csv
import math
import os
import shutil
import subprocess
import sys

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))

CSV = os.path.join(REPO_ROOT, "data", "m2_check7.csv")
DAMAGE_CSV = os.path.join(REPO_ROOT, "data", "m2_check7_damage.csv")
CASES = [("fx", "rz", 2.2), ("fx", "3d", 3.0), ("c4", "rz", 2.2), ("c4", "rz", 1.1), ("c4", "3d", 3.0),
         ("c5", "rz", 2.2), ("c5", "rz", 1.1), ("c5", "3d", 2.2), ("c6", "rz", 2.2), ("c6", "rz", 1.1),
         ("c6", "3d", 2.2)]
FRAMES = {"fx": 2, "c4": 3, "c5": 10, "c6": 3}

# --- the free uniform expansion (patch test): largest von Mises over E alpha_eff dT at k 2 (held)
FX_PASS = 6e-4                  # [SPH RZ 5.7e-4, SPH 3D 3.0e-4, CRKSPH 3D 2.8e-4; the held frame k 2]
FX_RZ_CRKSPH_FAIL = 1.0         # [1.07 interior] SolidCRKSPHRZ: stressed by about E alpha dT, not stress free
FX_RZ_CRKSPH_DILATATION = 1.2   # [1.232-1.238] its outer radial strain over alpha_eff dT (SPH 0.9996-0.9998)
# --- check 4 (over sigma_0): 3D CRKSPH against SPH at 3.0 mm
C4_3D = {"sph": (0.04, 0.07),      # [interior 0.0363, shell 0.0624]
         "crksph": (0.02, 0.08)}   # [0.0162, 0.0740]: the interior halves, the shell does not improve
C4_RZ_CRKSPH_INTERIOR = 0.3     # [0.33-0.62] SolidCRKSPHRZ's interior error is above this at both spacings
C4_RZ_SPH_INTERIOR = 0.05       # [0.0448]
# --- check 5 (over Y): 3D 2.2 mm
C5_3D_ELASTIC = {"sph": 0.02, "crksph": 0.01}          # [k 1: 0.0146, 0.0093]
C5_3D_PLATEAU = {"sph": 2e-3, "crksph": 0.09}          # [k >= 3: 0.00136, 0.0812 (k 10, after the monitor's switch)]
C5_3D_LATERAL = {"sph": 0.2, "crksph": 0.4}            # [0.108, 0.353]
C5_RZ_CRKSPH_LATERAL = 0.7      # [0.73-0.99] SolidCRKSPHRZ's lateral stresses exceed this (SPH RZ 1.1 mm 0.018)
# --- check 6 (over VM0), 3D 2.2 mm, patch shares
C6_3D = {"sph": dict(err_interior=0.7, err_interior_rms=0.04, err_core=0.04, err_shell=3.0),    # [0.646, 0.0395, 0.0400, 2.99]
         "crksph": dict(err_interior=1.0, err_interior_rms=0.06, err_core=0.05, err_shell=3.0)}  # [0.949, 0.0584, 0.0421, 2.58]
C6_DRAG = 9e-4                  # [7.8e-4] 3D, both schemes: the drag against the table's integral
C6_RZ_CRKSPH_DIVERGED = 1e3     # [1.1e5] RZ 1.1 mm CRKSPH: the interior error after frame 2 (the run diverged)
# --- clumping (Step 2): no pair under tension closer than half a spacing, in any run
TENSION_MIN_PAIR_DX = 0.9       # [0.913, CRKSPH RZ 1.1 mm check 4/6; SPH 0.937, RZ 1.1 mm check 5 k 10]


def rows():
    with open(CSV, newline="") as fh:
        out = []
        for r in csv.DictReader(fh):
            for c in r:
                if c not in ("check", "scheme", "form", "run"):
                    r[c] = float(r[c])
            r["k"] = int(r["k"])
            out.append(r)
    return out


def run_rows(check, scheme, form, dx):
    return sorted([r for r in rows() if (r["check"], r["scheme"], r["form"], r["dx_mm"]) == (check, scheme, form, dx)],
                  key=lambda r: r["k"])


def worst(rs, c, sel=lambda r: True):
    return max(abs(r[c]) for r in rs if sel(r))


def test_table_is_complete():
    for check, form, dx in CASES:
        for scheme in ("sph", "crksph"):
            rs = run_rows(check, scheme, form, dx)
            assert [r["k"] for r in rs] == list(range(1, FRAMES[check] + 1)), (check, scheme, form, dx)


def test_free_expansion_is_stress_free_except_in_solidcrksphrz():
    for scheme, form, dx in (("sph", "rz", 2.2), ("sph", "3d", 3.0), ("crksph", "3d", 3.0)):
        r = run_rows("fx", scheme, form, dx)[-1]
        assert max(r["err_interior"], r["err_shell"]) < FX_PASS, (scheme, form)
    for r in run_rows("fx", "crksph", "rz", 2.2):
        assert r["err_interior"] > FX_RZ_CRKSPH_FAIL
        assert r["fx_dilatation_ratio"] > FX_RZ_CRKSPH_DILATATION


def test_check4_heated_sphere():
    for scheme, (ti, ts) in C4_3D.items():
        rs = run_rows("c4", scheme, "3d", 3.0)
        assert worst(rs, "err_interior") < ti and worst(rs, "err_shell") < ts, scheme
    sph, crk = run_rows("c4", "sph", "3d", 3.0), run_rows("c4", "crksph", "3d", 3.0)
    assert all(c["err_interior"] < s["err_interior"] for s, c in zip(sph, crk))     # CRKSPH's 3D interior is better
    assert all(c["err_shell"] > s["err_shell"] for s, c in zip(sph, crk))           # ... its shell is not
    for dx in (2.2, 1.1):
        assert min(r["err_interior"] for r in run_rows("c4", "crksph", "rz", dx)) > C4_RZ_CRKSPH_INTERIOR
        assert worst(run_rows("c4", "sph", "rz", dx), "err_interior") < C4_RZ_SPH_INTERIOR


def test_check5_bar():
    for scheme in ("sph", "crksph"):
        rs = run_rows("c5", scheme, "3d", 2.2)
        assert abs(rs[0]["sxx_err_over_Y"]) < C5_3D_ELASTIC[scheme], scheme
        assert worst(rs, "sxx_err_over_Y", lambda r: r["k"] >= 3) < C5_3D_PLATEAU[scheme], scheme
        assert worst(rs, "lat_max_over_Y") < C5_3D_LATERAL[scheme], scheme
    for dx in (2.2, 1.1):
        assert worst(run_rows("c5", "crksph", "rz", dx), "lat_max_over_Y") > C5_RZ_CRKSPH_LATERAL


def test_check6_loads():
    for scheme, lim in C6_3D.items():
        rs = run_rows("c6", scheme, "3d", 2.2)
        for c, t in lim.items():
            assert worst(rs, c) < t, (scheme, c)
        assert worst(rs, "drag_rel") < C6_DRAG
    assert run_rows("c6", "crksph", "rz", 1.1)[-1]["err_interior"] > C6_RZ_CRKSPH_DIVERGED


def test_no_clumping_under_tension():
    rs = rows()
    assert all(r["tension_share_lt_half"] == 0.0 for r in rs if not math.isnan(r["tension_share_lt_half"]))
    assert min(r["tension_min_pair_dx"] for r in rs if not math.isnan(r["tension_min_pair_dx"])) > \
        TENSION_MIN_PAIR_DX
    # tension is present to test: every check-4 and check-6 frame has particles under tension
    assert all(r["n_tension"] > 0 for r in rs if r["check"] in ("c4", "c6"))


def test_damage_runs_with_both_schemes():
    with open(DAMAGE_CSV, newline="") as fh:
        d = {(r["scheme"], r["form"]): r for r in csv.DictReader(fh)}
    assert set(d) == {(s, f) for s in ("sph", "crksph") for f in ("rz", "3d")}
    for key, r in d.items():
        assert r["constructed"] == "1" and r["stepped"] == "1" and r["steps_done"] == "400" and not r["error"], key
        assert "ProbabilisticDamageModel" in r["packages"], key
        assert float(r["damage_max"]) > 2e-3, key        # [0.0028, CRKSPH RZ; 1.0 in the other three]


def test_run_names_carry_the_scheme():
    import spheral_m2_check7 as c7
    from spheral_frag import naming
    for check, form, dx in CASES:
        d = c7.run_dir(check, form, dx, "crksph")
        assert "_scheme-crksph_" in naming.run_name_from_dir(d)
        assert all(len(p) <= naming.NAME_MAX for p in d.split(os.sep))


@pytest.mark.spheral
def test_crksph_rz_rows_reproduce_through_a_restart():
    """Check 6's RZ 2.2 mm CRKSPH case over frames 0-2, stopped after frame 1 and resumed (1 process, the study's
    decomposition): frames 1-2 equal the committed rows (the restart at a frame boundary works under SolidCRKSPH)."""
    import spheral_m2_check7 as c7
    prep = os.path.join(REPO_ROOT, c7.PREPARED)
    body = os.path.join(REPO_ROOT, c7.body_dir("c6", "rz", 2.2))
    if not all(os.path.isfile(p) for p in (os.path.join(prep, "prepare.json"), os.path.join(body, "body.json"))):
        pytest.skip("the MVP flight's prepared input or the RZ 2.2 mm body is not built (git-ignored)")
    out = os.path.join(REPO_ROOT, "spheral_output", "tests", "check7_rz")
    shutil.rmtree(out, ignore_errors=True)
    rel = os.path.relpath(out, REPO_ROOT)
    # the committed rows were made under decision 7's first rule: interval 0 switched once to 20 crossings (11.5 %)
    # and interval 1 ran at 20 unmonitored. Task 10b's rule keeps doubling there (16.1 % at 20), so the rows are
    # reproduced at 20 crossings without the monitor (a switched interval's state is a plain 2 N run's, bitwise)
    cmd = [spheral_launcher()] + c7.command("c6", "rz", 2.2, out=rel, fr=(0, 2)) + ["--crossings", "20",
                                                                                     "--monitor-every", "0"]
    r = subprocess.run(cmd + ["--stop-at", "1"], cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    r = subprocess.run(cmd + ["--resume"], cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    from spheral_frag import naming
    d = naming.run_dir(out, naming.run_name(os.path.basename(c7.PREPARED), "replay", (0, 2), "rz", 2.2,
                                            naming.canonical_brackets(dict(c7.brackets("c6", "crksph"), crossings=20,
                                                                           monitor_every=0)), 1))
    ctx = c7.Context()
    orig = c7.run_dir, c7.frames
    try:
        c7.run_dir = lambda *a: d
        c7.frames = lambda check: (0, 2)
        got = c7.measure_run(ctx, "c6", "rz", 2.2, "crksph")
    finally:
        c7.run_dir, c7.frames = orig
    assert len(got) == 2
    want = run_rows("c6", "crksph", "rz", 2.2)
    for g, w in zip(got, want):
        for c in ("err_interior", "err_interior_rms", "err_core", "err_shell", "drag_rel"):
            assert g[c] == pytest.approx(w[c], rel=1e-9), (g["k"], c)
    shutil.rmtree(out, ignore_errors=True)
