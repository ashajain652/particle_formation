"""Paths and helpers shared by the tests (importable because pytest puts tests/ on sys.path)."""
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(REPO_ROOT, "tests", "fixtures")
PY = "/Users/ashajain/miniforge3/envs/drama_env/bin/python"


def fixture_dir(name):
    return os.path.join(FIXTURES, name)


def fixture_file(name, suffix):
    """The single file in tests/fixtures/<name>/ whose name ends with `suffix`."""
    hits = [f for f in os.listdir(fixture_dir(name)) if f.endswith(suffix)]
    assert len(hits) == 1, (name, suffix, hits)
    return os.path.join(fixture_dir(name), hits[0])
