"""Pytest configuration for the wake-residence tests: the project folder ("Reentry Model/") on sys.path, and the
parent trajectory and its twin (spec 4.3) built once per session (same pattern as the wake-column tests)."""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # "Reentry Model/"

if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


@pytest.fixture(scope="session")
def parent():
    from reentry_model.droplet_reentry import parent_trajectory
    return parent_trajectory.load_parent()


@pytest.fixture(scope="session")
def twin(parent):
    """The twin of the DRAMA parent (spec 4.3), integrated once per session from the branch top (about 111 km) down
    to 45 km; it takes some seconds, so every test that needs the parent track shares this one."""
    from reentry_model.parent_follower_separation_interactions.adapters import twin_parent
    return twin_parent(parent)
