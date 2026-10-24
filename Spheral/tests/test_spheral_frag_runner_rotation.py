"""spheral_frag/m0/rigid_rotation.py: SolidSPH turns the deviatoric stress with a rigidly rotating body.

Reduced size (12 x 12 particles, 64 interior, 50 steps per quarter turn), same 1 % tolerances as the full check.
At 50 steps the measured error is 0.155 % of |S0| and the von Mises drift 0.010 % (2026-10-05, 116c71f), both
second-order time-integration error; a wrong-sign spin term would be off by 2 |S0| at an eighth of a turn."""
import json
import os
import shutil
import subprocess
import tempfile

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

TOL = 0.01


@pytest.mark.spheral
def test_rigid_rotation_turns_stress_with_body():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")
    os.makedirs(base, exist_ok=True)
    out = tempfile.mkdtemp(prefix="rotation_", dir=base)
    try:
        r = subprocess.run([spheral_launcher(), "spheral_frag/m0/rigid_rotation.py", "--nx", "12", "--steps", "50",
                            "--tol", str(TOL), "--outdir", os.path.relpath(out, REPO_ROOT)],
                           cwd=REPO_ROOT, capture_output=True, text=True, timeout=600)
        runs = os.listdir(out)
        assert len(runs) == 1, (r.returncode, r.stdout[-2000:], r.stderr[-2000:])
        s = json.load(open(os.path.join(out, runs[0], "summary.json")))
    finally:
        shutil.rmtree(out, ignore_errors=True)

    assert s["steps"] == 50 and s["interior_particles"] == 64
    assert abs(s["theta_final_rad"] - 1.5707963267948966) < 0.01        # a quarter turn, measured from positions
    assert s["max_err_over_history"] < TOL                              # S = R S0 R^T within 1 % of |S0|
    assert s["max_von_mises_drift_over_history"] < TOL                  # von Mises within 1 % of its start
    assert s["wrong_sign_max_err_over_history"] > 1.5                   # the test would catch a wrong-sign spin term
    assert s["no_spin_max_err_over_history"] > 1.5                      # ... and a missing one
    assert s["max_radial_strain"] < 1e-4                                # the imposed rotation stayed rigid
    assert r.returncode == 0, r.stdout[-2000:]
