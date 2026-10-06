"""Pytest configuration: repo root on sys.path, 'drama' and 'spheral' marker auto-skip.

The 'spheral' probe runs `$SPHERAL -c "import Spheral"` (default: the container launcher) once per session, and only
when a 'spheral' test was collected. The launcher exits 2 within a second when Docker is stopped or the image is
absent, so the probe costs the normal loop nothing then; the 120 s timeout bounds a Docker daemon that hangs.

Both mesh fixtures pin band = 0: they are the analytic and fast-solver devices and must not be re-meshed by a
change to the default band. A default-band mesh perturbed the Carslaw-Jaeger fixture by five nodes and tipped its
centre-temperature margin, which was only 0.06 percentage points wide (measured 2026-09-27)."""
import importlib.util
import os
import subprocess
import sys

import pytest

from helpers import REPO_ROOT

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

HAVE_DRAMA = importlib.util.find_spec("drama") is not None
SPHERAL_PROBE_TIMEOUT_S = 120
_spheral_probe = {}


def spheral_launcher():
    return os.environ.get("SPHERAL", os.path.join(REPO_ROOT, "spheral_frag", "container", "spheral"))


def probe_spheral(launcher):
    """(ok, message) for `launcher -c "import Spheral"`, cached per launcher for the session."""
    if launcher not in _spheral_probe:
        try:
            r = subprocess.run([launcher, "-c", "import Spheral"], cwd=REPO_ROOT, capture_output=True, text=True,
                               timeout=SPHERAL_PROBE_TIMEOUT_S)
            lines = (r.stderr or r.stdout).strip().splitlines()
            _spheral_probe[launcher] = (r.returncode == 0, f"{launcher} -c 'import Spheral' exited {r.returncode}"
                                        + (f": {lines[-1]}" if lines else ""))
        except (OSError, subprocess.TimeoutExpired) as e:
            _spheral_probe[launcher] = (False, f"{launcher} -c 'import Spheral' failed: {e}")
    return _spheral_probe[launcher]


def pytest_collection_modifyitems(config, items):
    if not HAVE_DRAMA:
        skip = pytest.mark.skip(reason="pyDRAMA (package 'drama') is not importable in this interpreter")
        for item in items:
            if "drama" in item.keywords:
                item.add_marker(skip)
    spheral_items = [item for item in items if item.get_closest_marker("spheral")]
    if spheral_items:
        ok, message = probe_spheral(spheral_launcher())
        if not ok:
            skip = pytest.mark.skip(reason=message)
            for item in spheral_items:
                item.add_marker(skip)


@pytest.fixture(scope="session")
def coarse_sphere_mesh(tmp_path_factory):
    """100 mm sphere, 4 mm surface / 20 mm core (~2.9 k nodes): the mesh for the fast solver tests."""
    from reentry_model import mesh
    return mesh.sphere_mesh(0.05, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), band=0.0)


@pytest.fixture(scope="session")
def uniform_test_mesh(tmp_path_factory):
    """100 mm sphere, uniform 4 mm elements (~6.7 k nodes): resolves the centre for the Carslaw-Jaeger test."""
    from reentry_model import mesh
    return mesh.sphere_mesh(0.05, 4e-3, 4e-3, str(tmp_path_factory.mktemp("meshes")), band=0.0)
