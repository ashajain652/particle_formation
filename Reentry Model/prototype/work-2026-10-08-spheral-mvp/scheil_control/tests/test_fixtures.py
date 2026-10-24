import os

import pytest

from helpers import REPO_ROOT, fixture_file

CASES = ["T1_demised_50mm", "T3_survivor_50mm", "T5_5mm_750K", "E1_ballooning_5mm"]
SUFFIXES = ("_AeroThermalHistory.txt", "_Trajectory.txt", "ImpactingFragments.xml", "sesam.log")


@pytest.mark.parametrize("case", CASES)
def test_fixture_has_the_four_sesam_files(case):
    for suffix in SUFFIXES:
        assert os.path.getsize(fixture_file(case, suffix)) > 0


def test_survivor_fixture_has_full_precision_mass():
    with open(fixture_file("T3_survivor_50mm", "ImpactingFragments.xml")) as fh:
        assert '<mass unit="kg">1.8411041946975187e-01</mass>' in fh.read()


@pytest.mark.parametrize("name", ["fap_day.dat", "fap_mon.dat"])
def test_space_weather_files_present(name):
    with open(os.path.join(REPO_ROOT, "data", name)) as fh:
        assert fh.readline().startswith("#")


@pytest.mark.drama
def test_drama_marker_is_registered():
    import drama  # noqa: F401
