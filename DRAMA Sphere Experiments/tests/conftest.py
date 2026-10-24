"""Pytest configuration: the project folder ("DRAMA Sphere Experiments/", which holds sphere_reentry.py and
sphere_sweep.py) on sys.path, and the 'drama' marker auto-skip when pyDRAMA is not importable."""
import importlib.util
import sys

import pytest

from helpers import PROJECT_ROOT

if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

HAVE_DRAMA = importlib.util.find_spec("drama") is not None


def pytest_collection_modifyitems(config, items):
    if not HAVE_DRAMA:
        skip = pytest.mark.skip(reason="pyDRAMA (package 'drama') is not importable in this interpreter")
        for item in items:
            if "drama" in item.keywords:
                item.add_marker(skip)
