"""spheral_frag/m0/bench.py: the timing script of Spheral M0 Task 5, at reduced size.

3D at 10 mm (515 particles, -1.6 % of V/dx^3) and RZ at 5 mm (156 particles, -0.7 % of (pi R^2/2)/dx^2), 5 timed
steps each with damage on, so the strength, the damage model and both geometries are exercised. These are not
measurements (the timed matrix is Task 6), so they pass --force and do not need an idle machine."""
import json
import math
import os
import shutil
import subprocess
import tempfile

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

SCRIPT = "spheral_frag/m0/bench.py"


def bench(out, *args, nproc=1, timeout=600):
    cmd = [spheral_launcher()] + (["-n", str(nproc)] if nproc > 1 else []) + [SCRIPT, *args, "--out",
                                                                              os.path.relpath(out, REPO_ROOT)]
    return subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout)


@pytest.fixture
def outdir():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")
    os.makedirs(base, exist_ok=True)
    out = tempfile.mkdtemp(prefix="bench_", dir=base)
    yield out
    shutil.rmtree(out, ignore_errors=True)


@pytest.mark.spheral
@pytest.mark.parametrize("geometry, dx_mm, expected", [
    ("3d", 10.0, 4.0/3.0*math.pi*50.0**3/10.0**3),
    ("rz", 5.0, 0.5*math.pi*50.0**2/5.0**2),
])
def test_bench_summary_and_skip(outdir, geometry, dx_mm, expected):
    args = ["--geometry", geometry, "--dx-mm", str(dx_mm), "--damage", "on", "--steps", "5", "--warmup", "2",
            "--force"]
    r = bench(outdir, *args)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    runs = os.listdir(outdir)
    assert len(runs) == 1 and runs[0].startswith(f"bench_{geometry}_dx{dx_mm:.3f}mm_R50mm_n1_dmgon_steps5_wu2_")
    s = json.load(open(os.path.join(outdir, runs[0], "summary.json")))

    assert abs(s["particles"]/expected - 1.0) < 0.05                    # the plan's 5 % on the particle count
    assert s["particles_expected"] == pytest.approx(expected)
    assert s["processes"] == 1 and s["steps"] == 5 and s["damage"] == "on" and s["geometry"] == geometry
    for k in ("import_s", "material_s", "body_generation_s", "setup_s", "startup_s", "warmup_s",
              "step_wall_total_s", "processor_s_per_particle_step", "peak_rss_bytes_total",
              "memory_bytes_per_particle"):
        assert s[k] > 0, k
    for k in ("mean", "median", "min", "max"):
        assert s["step_wall_s"][k] > 0 and s["dt_s"][k] > 0, k
    assert s["processor_s_per_particle_step"] == pytest.approx(s["step_wall_s"]["median"]/s["particles"])
    assert s["material"] == "timing only" and s["machine"] == "aarch64"
    assert s["pin"]["spheral_commit"].startswith(s["spheral_commit"])
    # Spheral's own dt is the Courant step of the real material, not one inflated by a wrong sound speed:
    # c_l = sqrt((K + 4/3 G)/rho0) = 6.17 km/s, h = 2.01 dx, so 0.25 h/c_l = 8.1e-7 s per 10 mm.
    assert 0.2 < s["dt_s"]["max"]/(0.25*2.01*1e-3*dx_mm/6171.0) < 1.5
    with open(os.path.join(outdir, runs[0], "steps.csv")) as f:
        assert len(f.read().strip().splitlines()) == 1 + 5

    r2 = bench(outdir, *args)
    assert r2.returncode == 0 and "skipped" in r2.stdout, r2.stdout[-2000:]


@pytest.mark.spheral
def test_bench_refuses_busy_machine(outdir):
    """Above --max-load the run refuses to start (exit 2) and writes nothing; -1 makes any machine 'busy'."""
    r = bench(outdir, "--geometry", "rz", "--dx-mm", "5", "--damage", "off", "--steps", "1", "--max-load", "-1")
    assert r.returncode == 2, r.stdout[-2000:] + r.stderr[-2000:]
    assert "not starting" in r.stdout and os.listdir(outdir) == []


@pytest.mark.spheral
def test_bench_bad_arguments_exit_2(outdir):
    r = bench(outdir, "--geometry", "2d", "--dx-mm", "5", "--damage", "off")
    assert r.returncode == 2
