"""spheral_frag/m0/subclass_check.py: the Python subclasses of Spheral M0 Task 7 against the built-ins, at reduced size.

Method by method on synthetic inputs that reach every branch (both density clamps, melted and cold energies, damage
eigenvalues beyond [0, 1], excluded nodes), and end to end on Task 5's body at reduced size on two processes (so
that ghost points and the distributed boundary are exercised): 3D at 10 mm (515 particles) and RZ at 5 mm (156),
5 timed steps after 2 warm-up steps, each Python class alone against the built-in run. The plan's tolerance is
1e-12 relative; on the pinned image every result is bitwise identical, and the tests also say so where they can.
Not timed: they pass --force."""
import json
import os
import shutil
import subprocess
import tempfile

import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

SCRIPT = "spheral_frag/m0/subclass_check.py"
FIELDS = ("position", "velocity", "deviatoric_stress", "pressure", "damage", "density", "specific_thermal_energy",
          "plastic_strain")
RUNS = {"3d": ("10", ("builtin", "eos", "strength", "damage")),
        "rz": ("5", ("builtin", "eos", "strength", "damage", "linpoly_builtin", "linpoly_python"))}
ALONE = {"3d": ("eos", "strength", "damage"), "rz": ("eos", "strength", "damage", "linpoly_python")}


def check(*args, nproc=1, timeout=900):
    cmd = [spheral_launcher()] + (["-n", str(nproc)] if nproc > 1 else []) + [SCRIPT, *args]
    return subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout)


def rel(path):
    return os.path.relpath(path, REPO_ROOT)


@pytest.fixture(scope="module")
def outdir():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")
    os.makedirs(base, exist_ok=True)
    out = tempfile.mkdtemp(prefix="subclass_", dir=base)
    yield out
    shutil.rmtree(out, ignore_errors=True)


@pytest.fixture(scope="module")
def comparison(outdir):
    """Every reduced-size run, then compare -> the comparison summary."""
    for geometry, (dx, variants) in RUNS.items():
        for v in variants:
            r = check("run", "--geometry", geometry, "--dx-mm", dx, "--variant", v, "--steps", "5", "--warmup", "2",
                      "--force", "--out", rel(outdir), nproc=2)
            assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    r = check("compare", "--out", rel(outdir))
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    s = json.load(open(os.path.join(outdir, "summary.json")))
    return {(c["geometry"], cmp["variant"]): (c, cmp) for c in s["cases"] for cmp in c["comparisons"]}


@pytest.mark.spheral
@pytest.mark.parametrize("geometry", ["3d", "rz"])
def test_methods_bitwise(outdir, geometry):
    path = os.path.join(outdir, f"methods_{geometry}.json")
    r = check("methods", "--geometry", geometry, "--n", "2000", "--json", rel(path))
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    m = json.load(open(path))
    groups = m["methods"]
    assert set(groups) == {"eos_murnaghan", "eos_linpoly", "strength_steinberg_guinan", "damage_rate"}
    assert len(groups["eos_murnaghan"]) == 11 and len(groups["strength_steinberg_guinan"]) == 5
    for grp, d in groups.items():
        for name, diff in d.items():
            assert diff["points"] == 2000
            assert diff["identical"], f"{grp}.{name}: {diff}"


@pytest.mark.spheral
@pytest.mark.parametrize("geometry, variant", [(g, v) for g, vs in ALONE.items() for v in vs])
def test_each_subclass_alone_gives_the_builtin_result(comparison, geometry, variant):
    case, c = comparison[(geometry, variant)]
    assert case["processes"] == 2 and case["steps"] == 5
    assert c["same_particles"] and c["same_decomposition"]
    for f in FIELDS:
        d = c["fields"][f]
        assert d["max_pointwise_relative"] <= 1e-12, (f, d)       # the plan's tolerance
        assert d["identical"], (f, d)                             # what the pinned image gives
    assert c["same_dt_statistics"]
    assert c["median_step_s"][1] > 0 and c["python_override_s_per_step_max_rank"] > 0


@pytest.mark.spheral
def test_run_skipped_when_summary_exists(outdir, comparison):
    r = check("run", "--geometry", "rz", "--dx-mm", "5", "--variant", "builtin", "--steps", "5", "--warmup", "2",
              "--force", "--out", rel(outdir), nproc=2)
    assert r.returncode == 0 and "skipped" in r.stdout, r.stdout[-2000:]


@pytest.mark.spheral
def test_refuses_busy_machine(outdir):
    r = check("run", "--geometry", "rz", "--dx-mm", "5", "--variant", "builtin", "--steps", "1", "--max-load", "-1",
              "--out", rel(os.path.join(outdir, "busy")))
    assert r.returncode == 2 and "not starting" in r.stdout, r.stdout[-2000:] + r.stderr[-2000:]
    assert not os.path.exists(os.path.join(outdir, "busy"))


@pytest.mark.spheral
@pytest.mark.parametrize("args", [["run", "--geometry", "2d", "--dx-mm", "5", "--variant", "builtin"],
                                  ["run", "--geometry", "rz", "--dx-mm", "5", "--variant", "nonsense"],
                                  ["nonsense"]])
def test_bad_arguments_exit_2(args):
    assert check(*args).returncode == 2
