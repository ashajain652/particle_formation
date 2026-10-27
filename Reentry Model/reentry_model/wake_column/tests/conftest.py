"""Pytest configuration for the wake-column tests: the project folder ("Reentry Model/") on sys.path, the fixture
folder, and the parent trajectory loaded once per session (same pattern as the droplet tests)."""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # "Reentry Model/"
FIXTURES = os.path.join(HERE, "fixtures")

if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


@pytest.fixture(scope="session")
def parent():
    from reentry_model.droplet_reentry import parent_trajectory
    return parent_trajectory.load_parent()


@pytest.fixture(scope="session")
def fixtures_dir():
    return FIXTURES
