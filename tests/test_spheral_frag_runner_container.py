"""spheral_frag/container/spheral: the launcher that runs Spheral inside the pinned image.

Spheral silences print() on ranks other than 0 once it is imported, so the MPI test writes with os.write."""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

PIN = json.load(open(os.path.join(REPO_ROOT, "spheral_frag", "container", "pin.json")))


def run(args, cwd=REPO_ROOT, env=None, timeout=300):
    return subprocess.run([spheral_launcher()] + args, cwd=cwd, capture_output=True, text=True, timeout=timeout,
                          env=env)


@pytest.mark.spheral
def test_imports_pinned_commit():
    """The short hash Spheral was built from (SpheralConfigs.git_hash, also in the banner) is pin.json's commit."""
    r = run(["-c", "import Spheral, SpheralConfigs, os; os.write(1, ('hash=' + SpheralConfigs.git_hash() + '\\n').encode())"])
    assert r.returncode == 0, r.stderr
    short = [line.split("=", 1)[1] for line in r.stdout.splitlines() if line.startswith("hash=")]
    assert len(short) == 1 and len(short[0]) >= 7 and PIN["spheral_commit"].startswith(short[0]), r.stdout
    assert f"risky {short[0]} HEAD" in r.stdout


@pytest.mark.spheral
def test_two_ranks():
    """-n 2: both ranks see size 2 and an allreduce over shared memory sums their ranks."""
    code = ("import Spheral, os\n"
            "from mpi4py import MPI\n"
            "c = MPI.COMM_WORLD\n"
            "os.write(1, f'rank {c.rank} size {c.size} sum {c.allreduce(c.rank)}\\n'.encode())\n")
    r = run(["-n", "2", "-c", code])
    assert r.returncode == 0, r.stderr
    lines = sorted(line for line in r.stdout.splitlines() if line.startswith("rank "))
    assert lines == ["rank 0 size 2 sum 1", "rank 1 size 2 sum 1"], r.stdout


@pytest.mark.spheral
def test_relative_paths_from_subdirectory():
    """Run from a subdirectory of the repository, a relative script path and a relative input file both resolve."""
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")
    os.makedirs(base, exist_ok=True)
    case = tempfile.mkdtemp(prefix="relpath_", dir=base)
    try:
        os.makedirs(os.path.join(case, "data"))
        with open(os.path.join(case, "data", "input.txt"), "w") as f:
            f.write("42\n")
        with open(os.path.join(case, "read_input.py"), "w") as f:
            f.write("import os\nos.write(1, ('value=' + open('data/input.txt').read().strip() + '\\n').encode())\n")
        r = run(["read_input.py"], cwd=case)
        assert r.returncode == 0, r.stderr
        assert "value=42" in r.stdout.splitlines()
    finally:
        shutil.rmtree(case)


@pytest.mark.spheral
def test_missing_image_exits_2():
    r = run(["-c", "pass"], env={**os.environ, "SPHERAL_IMAGE": "particle-formation/spheral:does-not-exist"})
    assert r.returncode == 2
    assert len(r.stderr.strip().splitlines()) == 1 and "missing" in r.stderr


def test_probe_skips_when_launcher_fails():
    """With the launcher replaced by `false`, the conftest probe skips the 'spheral' tests instead of erroring."""
    env = {**os.environ, "SPHERAL": shutil.which("false")}
    r = subprocess.run([sys.executable, "-m", "pytest", __file__, "-m", "spheral", "-q", "-p", "no:cacheprovider"],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=120, env=env)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "4 skipped" in r.stdout and "exited 1" in r.stdout, r.stdout
