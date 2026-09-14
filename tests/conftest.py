"""Pytest configuration: repo root on sys.path, 'drama' marker auto-skip."""
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
