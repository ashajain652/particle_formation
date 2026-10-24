"""Check 4 (Spheral M2 plan, Task 6; spec §12): the heated sphere against Timoshenko and Goodier, with and without mass
scaling, and the crossings per interval and damping studies.

The normal-loop tests read the committed table `data/m2_check4.csv` (written by
`analysis/spheral_m2_check4.py`, whose docstring defines every measure) and hold it to the thresholds below, and
recompute the plan's three answers from it:

  - option A (the mass-scaled replay) is sound: each scaled run and its unscaled twin agree within the interior
    error of either against the exact solution;
  - the default crossings per interval: the smallest of 5, 10, 20 whose interior error is within 1.2 times that of
    20 and whose internal KE ratio stays below 5 %, at every form and spacing and both compared frames;
  - damping (0 and 0.2 against the default 0.05, at 10 crossings) changes the held state (k = 3) by less than the
    interior error.

The `spheral` test reruns one row (RZ 2.2 mm, 10 crossings, damping 0.05) and requires it back.

THRESHOLDS: measured 2026-10-10 (the CSV), each the worst case over the compared frames and the crossings and
damping runs of its form and spacing, rounded up to one significant figure; worst case in brackets. Change them only
together with the CSV, the README and docs/spheral_frag_assumptions.md (CLAUDE.md)."""
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

CSV = os.path.join(REPO_ROOT, "data", "m2_check4.csv")

# --- thresholds: fractions of the exact peak sigma_0 = 50.07 MPa, over the runs at the devices' working values
# (crossings 10 and 20, damping 0.05 and 0.2, scaled and unscaled), both compared frames; worst case in brackets.
# Crossings 5 and damping 0 are measured too and recorded in the CSV, but are not working values (tests below).
INTERIOR_ERR = {("rz", 2.2): 3e-2,     # [2.360e-2, crossings 20, k 1]
                ("rz", 1.1): 3e-2,     # [2.090e-2, crossings 10, k 1]
                ("3d", 3.0): 4e-2,     # [3.675e-2, damping 0.2, k 1]
                ("3d", 2.2): 4e-2}     # [3.177e-2, damping 0.2, k 1]
SHELL_ERR = {("rz", 2.2): 0.1,         # [9.775e-2] the outermost bin's hoop stress: SPH's free surface
             ("rz", 1.1): 0.09,        # [8.551e-2]
             ("3d", 3.0): 0.07,        # [6.242e-2]
             ("3d", 2.2): 0.08}        # [7.974e-2]
UNSCALED_DIFF = 4e-5                   # [3.089e-5, 3D 3.0 mm, crossings 5, k 3]: round-off and step sequence
KE_INTERNAL_DEFAULT = 4e-3             # [3.749e-3, RZ 1.1 mm, k 1] internal KE / elastic energy at the defaults
DAMPING_02_DIFF = 6e-3                 # [4.924e-3, RZ 2.2 mm] held state (k 3), damping 0.2 against 0.05
DEFAULT_CROSSINGS = 10                 # naming.BRACKET_DEFAULTS; the rule's per-case answers: test_default_crossings
DEFAULT_DAMPING = 0.05
CROSSINGS_FACTOR = 1.2                 # plan Task 6 Step 4
KE_FRACTION = 0.05                     # spec §9.1


def rows():
    with open(CSV, newline="") as fh:
        out = []
        for r in csv.DictReader(fh):
            for c in r:
                if c not in ("form", "run"):
                    r[c] = float(r[c])
            for c in ("crossings", "scaled", "k", "n_particles", "processes", "steps_run"):
                r[c] = int(r[c])
            out.append(r)
    return out


def pick(rs, **kw):
    return [r for r in rs if all(r[k] == v for k, v in kw.items())]


def test_table_is_complete():
    import spheral_m2_check4 as c4
    rs = rows()
    want = {(f, dx, n, z, int(s), k) for f, dx, _, n, z, s in c4.design() for k in c4.COMPARE_K}
    have = {(r["form"], r["dx_mm"], r["crossings"], r["damping"], r["scaled"], r["k"]) for r in rs}
    assert have == want
    assert all(abs(r["peak_exact_Pa"] - 50.0748e6) < 1e3 for r in rs)          # 2 alpha E dT / (5 (1 - nu))


def working(r):
    return r["crossings"] >= 10 and r["damping"] > 0.0


def test_errors_against_the_exact_solution():
    rs = rows()
    for r in rs:
        if working(r):
            key = (r["form"], r["dx_mm"])
            assert r["interior_err"] <= INTERIOR_ERR[key], r["run"]
            assert r["shell_err"] <= SHELL_ERR[key], r["run"]
    # the settled interior error (20 crossings, k 3) falls with spacing in both forms, if slowly (measured
    # 2.33 -> 1.69 % RZ, 2.64 -> 2.07 % 3D: a near-uniform isotropic offset of the core, +1.4 % RZ, -2 % 3D)
    held = {(r["form"], r["dx_mm"]): r["interior_err"]
            for r in pick(rs, crossings=20, damping=0.05, scaled=1, k=3)}
    assert held[("rz", 1.1)] < held[("rz", 2.2)] and held[("3d", 2.2)] < held[("3d", 3.0)]


def test_option_a_scaled_and_unscaled_agree():
    """Plan Task 6 Step 4: within the interior error of either against the exact solution (and below the measured
    threshold)."""
    rs = rows()
    twins = [r for r in rs if not r["scaled"]]
    assert len(twins) == 12
    for u in twins:
        s, = pick(rs, form=u["form"], dx_mm=u["dx_mm"], crossings=u["crossings"], damping=u["damping"], scaled=1,
                  k=u["k"])
        assert u["unscaled_diff_interior"] == pytest.approx(s["unscaled_diff_interior"], rel=1e-9)
        assert u["unscaled_diff_interior"] < min(u["interior_err"], s["interior_err"])
        assert u["unscaled_diff_interior"] <= UNSCALED_DIFF


def crossings_rule(rs, form, dx, k):
    """Plan Task 6 Step 4: the smallest of 5, 10, 20 whose interior error is within 1.2 times that of 20 and whose
    internal KE ratio is below 5 %."""
    e20 = pick(rs, form=form, dx_mm=dx, crossings=20, damping=0.05, scaled=1, k=k)[0]["interior_err"]
    for n in (5, 10, 20):
        r = pick(rs, form=form, dx_mm=dx, crossings=n, damping=0.05, scaled=1, k=k)[0]
        if r["interior_err"] <= CROSSINGS_FACTOR*e20 and r["ke_internal_ratio"] < KE_FRACTION:
            return n


def test_default_crossings():
    """Measured: 10 at every form, spacing and frame but one -- RZ 1.1 mm at the end of the ramp (k 1), where 10
    crossings give 2.09 % against 1.2 x 1.48 % (the ramp's quasi-static lag, 0.6 % of sigma_0); at the held state
    (k 3) 10 everywhere. 5 crossings never qualify (3-11 %). Recorded, not decided here (Task 6 report)."""
    rs = rows()
    got = {(f, dx, k): crossings_rule(rs, f, dx, k) for f, dx in INTERIOR_ERR for k in (1, 3)}
    assert {key for key, n in got.items() if n != DEFAULT_CROSSINGS} == {("rz", 1.1, 1)}
    assert got[("rz", 1.1, 1)] == 20
    assert all(n == DEFAULT_CROSSINGS for (f, dx, k), n in got.items() if k == 3)


def test_damping_does_not_move_the_held_state():
    """Damping 0.2 against 0.05 moves the held state (k 3) by less than the interior error; damping 0 does not hold
    a state at all (undamped ringing, internal KE/E 2e-3 to 1e-2 at k 3) and is not a working value."""
    rs = rows()
    for r in pick(rs, crossings=10, scaled=1, k=3):
        if r["damping"] == DEFAULT_DAMPING:
            continue
        ref = pick(rs, form=r["form"], dx_mm=r["dx_mm"], crossings=10, damping=DEFAULT_DAMPING, scaled=1, k=3)[0]
        if r["damping"] == 0.2:
            assert r["damping_diff_interior"] < min(r["interior_err"], ref["interior_err"]), r["run"]
            assert r["damping_diff_interior"] <= DAMPING_02_DIFF
        else:
            assert r["ke_internal_ratio"] > 1e-3 > ref["ke_internal_ratio"]
    for r in pick(rs, crossings=DEFAULT_CROSSINGS, damping=DEFAULT_DAMPING, scaled=1):
        assert r["ke_internal_ratio"] <= KE_INTERNAL_DEFAULT


@pytest.mark.spheral
def test_rz_row_reproduces(tmp_path):
    """RZ 2.2 mm, 10 crossings, damping 0.05, scaled: rerun into a scratch directory under spheral_output (the
    container sees only the repository) and measured as the driver does; the committed row comes back (one process:
    the same step sequence, so to round-off)."""
    import spheral_m2_check4 as c4
    from spheral_frag import material
    prep = os.path.join(REPO_ROOT, c4.PREPARED)
    body = os.path.join(REPO_ROOT, c4.body_dir("rz", 2.2))
    if not (os.path.isfile(os.path.join(prep, "prepare.json")) and os.path.isfile(os.path.join(body, "body.json"))):
        pytest.skip("the MVP flight's prepared input or the radial RZ 2.2 mm body is not built (git-ignored; "
                    "analysis/spheral_m2_check4.py builds the body)")
    out = os.path.join(REPO_ROOT, "spheral_output", "tests", "check4_rz")
    shutil.rmtree(out, ignore_errors=True)
    cmd = [spheral_launcher(), "spheral_frag/runner/run.py", "replay", "--body", c4.body_dir("rz", 2.2),
           "--prepared", c4.PREPARED, "--frames", "0-3", "--crossings", "10", "--max-crossings", "10",
           "--damping", "0.05", "--strength", "elastic", "--loads", "off", "--ke-fraction", "1e300",
           "--hourglass", repr(c4.HOURGLASS), "--out", os.path.relpath(out, REPO_ROOT)]
    r = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    bad = subprocess.run(cmd[:-2] + ["--mass-scaling", "off", "--loads", "on", "--out", "spheral_output/tests/x"],
                         cwd=REPO_ROOT, capture_output=True, text=True, timeout=600)
    assert bad.returncode == 2 and "--loads off" in bad.stderr
    table = material.MaterialTable.load(prep)
    E, nu, alpha = c4.constants(table)
    from spheral_frag import naming
    d = naming.run_dir(out, c4.run_name("rz", 2.2, 10, 0.05, True))
    for k in (1, 3):
        m = c4.measure(d, k, table, E, nu, alpha, 2.2e-3)
        row, = pick(rows(), form="rz", dx_mm=2.2, crossings=10, damping=0.05, scaled=1, k=k)
        assert m["interior"][0] == pytest.approx(row["interior_err"], rel=1e-9, abs=1e-12)
        assert m["shell"][0] == pytest.approx(row["shell_err"], rel=1e-9, abs=1e-12)
    shutil.rmtree(out, ignore_errors=True)


def test_check4_switches_name_the_run():
    """--loads off and --mass-scaling off are run-name brackets (default on) and parse back."""
    from spheral_frag import naming
    b = naming.canonical_brackets({"strength": "elastic", "sampler": "radial", "loads": "off", "mass_scaling": "off"})
    name = naming.run_name("prep_x_k00000-00003", "replay", (0, 3), "rz", 2.2, b, 1)
    assert name.endswith("_strength-elastic_sampler-radial_loads-off_massscaling-off")
    assert naming.parse_run_name(name)["brackets"] == b
    assert naming.BRACKET_DEFAULTS["loads"] == naming.BRACKET_DEFAULTS["mass_scaling"] == "on"
