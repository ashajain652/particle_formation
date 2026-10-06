"""LLNL's regression examples on the pinned arm64 image (Spheral M0 plan, Task 3; spec §12 check 1).

The two cheap members of the set; spheral_frag/m0/run_regressions.sh runs the whole set (4-process domain
independence, restarts, Taylor impact in 2d, RZ and 3d on 1 and 8 processes) and records the numbers."""
import os
import shutil
import subprocess
import tempfile

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

ROD = "/opt/spheral/tests/functional/Damage/TensileRod"
TAYLOR = "/opt/spheral/tests/functional/Strength/TaylorImpact/TaylorImpact.py"


@pytest.fixture
def workdir():
    """A scratch directory inside the repository (the container only sees the repository)."""
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")
    os.makedirs(base, exist_ok=True)
    case = tempfile.mkdtemp(prefix="regress_", dir=base)
    yield case
    shutil.rmtree(case, ignore_errors=True)


def run(args, cwd, timeout):
    return subprocess.run([spheral_launcher()] + args, cwd=cwd, capture_output=True, text=True, timeout=timeout)


@pytest.mark.spheral
@pytest.mark.parametrize("model", ["GradyKippTensorDamageOwen", "ProbabilisticDamageModel"])
def test_tensile_rod_matches_reference(workdir, model):
    """ATS t10 / t20: the serial 1-D tensile rod matches LLNL's stored reference to filearraycmp's 1e-4."""
    r = run(["-c", f"import shutil; shutil.copytree('{ROD}/Reference', 'Reference')"], workdir, 300)
    assert r.returncode == 0, r.stderr
    r = run([f"{ROD}/TensileRod-1d.py", "--DamageModelConstructor", model, "--graphics", "False",
             "--clearDirectories", "True", "--domainIndependent", "True",
             "--outputFile", "TensileRod-1d-1proc.gnu", "--checkRef", "True"], workdir, 900)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    assert "Floating point comparison test passed." in r.stdout


@pytest.mark.spheral
def test_taylor_impact_rz_runs(workdir):
    """TaylorImpact in RZ with SPH, 20 steps on one process, runs to completion and writes its snapshot."""
    r = run([TAYLOR, "--geometry", "RZ", "--hydroType", "SPH", "--steps", "20", "--compatibleEnergy", "False",
             "--clearDirectories", "True", "--siloSnapShotFile", "Spheral_sph_rz_state_snapshot_1proc"], workdir, 900)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    snapshots = [f for _, _, files in os.walk(os.path.join(workdir, "dumps-TaylorImpact")) for f in files
                 if f.startswith("Spheral_sph_rz_state_snapshot_1proc")]
    assert snapshots, "no snapshot written"
