"""Pytest configuration: repo root on sys.path, 'drama' marker auto-skip.

Both mesh fixtures pin band = 0: they are the analytic and fast-solver devices and must not be re-meshed by a
change to the default band. A default-band mesh perturbed the Carslaw-Jaeger fixture by five nodes and tipped its
centre-temperature margin, which was only 0.06 percentage points wide (measured 2026-09-27)."""
import importlib.util
import sys

import pytest

from helpers import REPO_ROOT

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

HAVE_DRAMA = importlib.util.find_spec("drama") is not None


def pytest_collection_modifyitems(config, items):
    if HAVE_DRAMA:
        return
    skip = pytest.mark.skip(reason="pyDRAMA (package 'drama') is not importable in this interpreter")
    for item in items:
        if "drama" in item.keywords:
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
