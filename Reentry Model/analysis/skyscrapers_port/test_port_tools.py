"""Tests of the skyscrapers port's tools (path map, checks, resolver).

Not collected by ``pytest tests``; run explicitly from ``Reentry Model``:
``python -m pytest analysis/skyscrapers_port/test_port_tools.py``.  They read refs
(``skyscrapers-build-end``, ``MERGE_BASE``, ``MAIN0``) from the shared repository.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

from path_map import map_path, unmap_path
from resolve_both_sides import three_way

TOOLS = Path(__file__).resolve().parent
MERGE_BASE = "d16018bf6bb19944a88c19cf0987c33951ea3485"
MAIN0 = "8ed75cfeeda2364f14dba2bf73cc48f58af6e72b"
BUILD_END = "skyscrapers-build-end"


def git(*args, cwd=TOOLS, env=None, check=True):
    return subprocess.run(["git", *args], cwd=cwd, env=env, check=check,
                          capture_output=True).stdout


def ls_tree(ref):
    """{path: (mode, blob)} of a tree-ish in the shared repository."""
    out = git("ls-tree", "-r", "-z", "--full-tree", ref)
    result = {}
    for entry in out.split(b"\0"):
        if entry:
            meta, path = entry.split(b"\t", 1)
            mode, _kind, blob = meta.decode().split()
            result[path.decode()] = (mode, blob)
    return result


def run_index_filter(line):
    return subprocess.run([sys.executable, "-I", str(TOOLS / "path_map.py"), "index-filter"],
                          input=line, capture_output=True, text=True, check=True).stdout


def test_build_head_files_map_to_the_designed_places():
    assert map_path("reentry_model/skinbody.py") == "Reentry Model/reentry_model/skinbody.py"
    assert map_path("reentry_model/data/atdb_disc.json") == "Reentry Model/data/atdb_disc.json"
    assert map_path("data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv") == \
        "Reentry Model/data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv"
    assert map_path("data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis.csv") == \
        "DRAMA Sphere Experiments/data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis.csv"
    assert map_path("tests/conftest.py") == "Reentry Model/tests/conftest.py"     # git's rename detection says Spheral/
    assert map_path("tests/test_sphere_sweep.py") == "DRAMA Sphere Experiments/tests/test_sphere_sweep.py"
    assert map_path("analysis/plot_style.py") == "DRAMA Sphere Experiments/analysis/plot_style.py"
    assert map_path(".gitignore") == "Reentry Model/.gitignore"


def test_files_mapped_outside_reentry_model_are_on_main_and_unchanged_by_the_build():
    build, base, main = ls_tree(BUILD_END), ls_tree(MERGE_BASE), ls_tree(MAIN0)
    outside = [p for p in build if not map_path(p).startswith("Reentry Model/")]
    assert outside
    for p in outside:
        q = map_path(p)
        assert q in main, f"main has no {q}"
        assert p in base and base[p][1] == build[p][1], f"the build changed {p}"


def test_the_map_is_one_to_one_and_unmap_inverts_it():
    paths = {p for p in git("log", "--format=", "--name-only", "-z", BUILD_END).decode().split("\0")
             if p.strip("\n")}
    paths = {p.strip("\n") for p in paths}
    assert len(paths) >= 148
    seen = {}
    for p in sorted(paths):
        q = map_path(p)
        assert q not in seen, f"{p} and {seen[q]} both map to {q}"
        seen[q] = p
        assert unmap_path(q) == p
    assert unmap_path("DRAMA Sphere Experiments/tests/test_unrelated.py") is None
    assert unmap_path("Reentry Model/tests/test_sphere_sweep.py") is None


def test_index_filter_maps_each_line():
    line = "100644 " + "a" * 40 + " 0\treentry_model/data/us76_static_environment.csv\n"
    assert run_index_filter(line) == "100644 " + "a" * 40 + " 0\tReentry Model/data/us76_static_environment.csv\n"


def test_three_way_takes_the_build_and_reapplies_main():
    assert three_way(b"a\nB\nc\n", b"a\nb\nc\n", b"a\nb\nc\nd\n") == (b"a\nB\nc\nd\n", True)
    assert three_way(b"a\nX\nc\n", b"a\nb\nc\n", b"a\nY\nc\n")[1] is False


def test_the_rewrite_check_passes_a_faithful_rewrite_and_fails_a_dropped_file(tmp_path):
    n = [0]

    def run(*args, cwd, **extra):
        n[0] += 1
        stamp = f"2026-01-0{n[0] % 9 + 1}T12:00:00 +0000" if "date" not in extra else extra.pop("date")
        env = dict(os.environ, GIT_AUTHOR_NAME="A U", GIT_AUTHOR_EMAIL="a@example.org",
                   GIT_COMMITTER_NAME="C O", GIT_COMMITTER_EMAIL="c@example.org",
                   GIT_AUTHOR_DATE=stamp, GIT_COMMITTER_DATE=stamp, PY=sys.executable, TOOLS=str(TOOLS),
                   FILTER_BRANCH_SQUELCH_WARNING="1", **extra)
        return git(*args, cwd=cwd, env=env)

    orig, clone = tmp_path / "orig", tmp_path / "clone"
    orig.mkdir()
    run("init", "-q", "-b", "skin-spraying", cwd=orig)
    (orig / "reentry_model").mkdir()
    (orig / "data").mkdir()
    (orig / "reentry_model/x.py").write_text("x = 1\n")
    (orig / "data/fap_day.dat").write_text("1 2 3\n")
    run("add", "-A", cwd=orig)
    run("commit", "-q", "-m", "first", cwd=orig)
    (orig / "reentry_model/x.py").write_text("x = 2\n")
    run("commit", "-q", "-a", "-m", "second\n\nbody line", cwd=orig)
    (orig / "reentry_model/y.py").write_text("y = 1\n")
    run("add", "-A", cwd=orig)
    run("commit", "-q", "-m", "third", cwd=orig)

    run("clone", "-q", str(orig), str(clone), cwd=tmp_path)
    run("filter-branch",
        "--index-filter",
        'git ls-files -s | "$PY" -I "$TOOLS/path_map.py" index-filter | GIT_INDEX_FILE="$GIT_INDEX_FILE.new" '
        'git update-index --index-info && mv "$GIT_INDEX_FILE.new" "$GIT_INDEX_FILE"',
        "--msg-filter", 'git interpret-trailers --trailer "Build-repository-commit: $GIT_COMMIT"',
        "--", "skin-spraying", cwd=clone)

    def check():
        return subprocess.run([sys.executable, str(TOOLS / "checks.py"), "rewrite", str(orig), str(clone)],
                              capture_output=True, text=True)

    ok = check()
    assert ok.returncode == 0, ok.stdout + ok.stderr
    assert (clone / "DRAMA Sphere Experiments/data/fap_day.dat").exists()

    # a rewritten tip that lost a file (same count, parents, dates and messages)
    tip = run("rev-parse", "HEAD", cwd=clone).decode().strip()
    run("rm", "-q", "--cached", "--", "DRAMA Sphere Experiments/data/fap_day.dat", cwd=clone)
    run("commit", "-q", "--amend", "--no-edit", cwd=clone)
    bad = check()
    assert bad.returncode != 0
    assert "fap_day.dat" in bad.stdout
    run("reset", "-q", "--hard", tip, cwd=clone)
    assert check().returncode == 0

    # a third commit on the clone only, deleting a file: the histories differ in length
    run("rm", "-q", "--", "Reentry Model/reentry_model/y.py", cwd=clone)
    run("commit", "-q", "-m", "delete y", cwd=clone)
    assert check().returncode != 0
