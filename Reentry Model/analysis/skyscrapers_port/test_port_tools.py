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


# ---------------------------------------------------------------------------------------------------------------
# same_results.py: the comparison of two result folders

import csv
import hashlib
import json

import numpy as np

import same_results


def _write_run(folder, *, wall=1.5, provenance="a", files=("x.vtu",), mesh_dir="m1", mesh_bytes=b"mesh",
               cell="0.1234567890123456", npz=None):
    """A small result folder: history.csv, summary.json (with settings.mesh_file under mesh_dir), fields.npz."""
    folder.mkdir(parents=True)
    mesh = folder.parent.parent / f"meshes_{folder.parent.name}_{mesh_dir}" / "sphere.msh"
    mesh.parent.mkdir(exist_ok=True)
    mesh.write_bytes(mesh_bytes)
    material = mesh.parent / "aa7075.json"
    material.write_bytes(b"{}")
    with open(folder / "history.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t_s", "T_K", "wall_skins_s"])
        writer.writerow(["0.025", cell, repr(wall)])
        writer.writerow(["0.05", "0.5", repr(wall * 2)])
    (folder / "summary.json").write_text(json.dumps({
        "provenance": {"git": provenance}, "files": list(files), "wall_total_s": wall,
        "result": {"peak": 1.0, "wall_inner": wall},
        "settings": {"mesh_file": str(mesh), "material_file": str(material), "dt": 0.025}}))
    np.savez(folder / "fields.npz", **(npz if npz is not None else {"T": np.array([1.0, np.nan, 3.0])}))
    (folder / "run.log").write_text(f"log {wall}")
    return mesh


def test_compare_ignores_timers_paths_and_provenance(tmp_path):
    _write_run(tmp_path / "a" / "run", wall=1.5, provenance="a", files=("x.vtu",), mesh_dir="m1")
    _write_run(tmp_path / "b" / "run", wall=9.75, provenance="b", files=("y.vtu", "z.vtu"), mesh_dir="m2")
    assert same_results.compare(tmp_path / "a", tmp_path / "b") == []


def test_compare_finds_one_changed_history_cell(tmp_path):
    _write_run(tmp_path / "a" / "run")
    _write_run(tmp_path / "b" / "run", cell="0.12345678901234568")      # the 17th significant digit
    found = same_results.compare(tmp_path / "a", tmp_path / "b")
    assert len(found) == 1 and "history.csv" in found[0] and "T_K" in found[0]


def test_compare_finds_a_changed_npz_array(tmp_path):
    _write_run(tmp_path / "a" / "run")
    _write_run(tmp_path / "b" / "run", npz={"T": np.array([1.0, np.nan, np.nextafter(3.0, 4.0)])})
    found = same_results.compare(tmp_path / "a", tmp_path / "b")
    assert len(found) == 1 and "fields.npz" in found[0]


def test_compare_finds_a_different_mesh_file_under_the_same_name(tmp_path):
    _write_run(tmp_path / "a" / "run", mesh_dir="m1", mesh_bytes=b"mesh")
    _write_run(tmp_path / "b" / "run", mesh_dir="m2", mesh_bytes=b"mesh!")
    found = same_results.compare(tmp_path / "a", tmp_path / "b")
    assert len(found) == 1 and "mesh_file" in found[0]


def test_compare_finds_a_file_in_one_folder_only_and_cli_exit_codes(tmp_path):
    _write_run(tmp_path / "a" / "run")
    _write_run(tmp_path / "b" / "run")
    assert same_results.main(["compare", str(tmp_path / "a"), str(tmp_path / "b")]) == 0
    (tmp_path / "b" / "run" / "extra.csv").write_text("a\n1\n")
    assert any("extra.csv" in line for line in same_results.compare(tmp_path / "a", tmp_path / "b"))
    assert same_results.main(["compare", str(tmp_path / "a"), str(tmp_path / "b")]) == 1
    assert same_results.main(["compare", str(tmp_path / "a"), str(tmp_path / "b"), "--report"]) == 0


def test_report_gives_largest_absolute_and_relative_difference(tmp_path, capsys):
    _write_run(tmp_path / "a" / "run", cell="2.0")
    _write_run(tmp_path / "b" / "run", cell="2.5")
    assert same_results.main(["compare", str(tmp_path / "a"), str(tmp_path / "b"), "--report"]) == 0
    out = capsys.readouterr().out
    assert "T_K" in out and "0.5" in out and "0.2" in out      # absolute 0.5, relative 0.5/2.5 = 0.2


def test_run_builds_the_documented_command_lines():
    cases = same_results.case_arguments(Path("/state.pkl"))
    assert set(cases) == {"skins", "element", "resume", "nomelt", "nomelt_elements"}
    assert cases["skins"][:1] == ["run"] and "--skin-preset" in cases["skins"]
    assert cases["element"][-2:] == ["--surface-model", "elements"]
    assert cases["resume"][:2] == ["resume", "/state.pkl"]
    assert "--melt" not in cases["nomelt"]
    assert cases["nomelt_elements"] == cases["nomelt"] + ["--surface-model", "elements"]


def test_compare_ignores_the_run_time(tmp_path):
    for name, runtime in (("a", 3.25), ("b", 7.5)):
        _write_run(tmp_path / name / "run")
        path = tmp_path / name / "run" / "summary.json"
        data = json.loads(path.read_text())
        data["results"] = {"runtime_s": runtime, "peak": 2.0}
        path.write_text(json.dumps(data))
    assert same_results.compare(tmp_path / "a", tmp_path / "b") == []
    path = tmp_path / "b" / "run" / "summary.json"
    data = json.loads(path.read_text())
    data["results"]["peak"] = 2.5
    path.write_text(json.dumps(data))
    assert len(same_results.compare(tmp_path / "a", tmp_path / "b")) == 1


def test_compare_refuses_a_missing_folder(tmp_path, capsys):
    _write_run(tmp_path / "a" / "run")
    for extra in ([], ["--report"]):
        assert same_results.main(["compare", str(tmp_path / "a"), str(tmp_path / "nope")] + extra) != 0
        assert same_results.main(["compare", str(tmp_path / "nope"), str(tmp_path / "nope2")] + extra) != 0
    assert "does not exist" in capsys.readouterr().out


def test_compare_refuses_an_empty_folder(tmp_path, capsys):
    _write_run(tmp_path / "a" / "run")
    (tmp_path / "empty").mkdir()
    (tmp_path / "logs_only").mkdir()
    (tmp_path / "logs_only" / "x.log").write_text("log")        # only ignored kinds: nothing to compare
    for other in ("empty", "logs_only"):
        for extra in ([], ["--report"]):
            assert same_results.main(["compare", str(tmp_path / "a"), str(tmp_path / other)] + extra) != 0
            assert same_results.main(["compare", str(tmp_path / other), str(tmp_path / other)] + extra) != 0
    assert "no file to compare" in capsys.readouterr().out


def test_run_refuses_a_case_folder_that_already_holds_results(tmp_path, capsys):
    held = tmp_path / "out" / "nomelt"
    held.mkdir(parents=True)
    (held / "nomelt.csv").write_text("a\n1\n")
    code = tmp_path / "code"
    code.mkdir()                                  # nothing is run: the refusal comes first
    assert same_results.main(["run", "--code", str(code), "--out", str(tmp_path / "out"), "--cases", "nomelt"]) != 0
    assert "already holds" in capsys.readouterr().out
    assert (held / "nomelt.csv").read_text() == "a\n1\n"       # nothing deleted
    assert not (tmp_path / "out" / "nomelt.log").exists()


def test_compare_treats_the_resumed_state_path_by_file_name(tmp_path):
    for side, folder in (("a", "x"), ("b", "y")):
        _write_run(tmp_path / side / "run", mesh_dir="m1")
        path = tmp_path / side / "run" / "summary.json"
        doc = json.loads(path.read_text())
        doc["resumed_from"] = {"state": f"/{folder}/saved/d100.state.pkl", "step": 3}
        path.write_text(json.dumps(doc))
    assert same_results.compare(tmp_path / "a", tmp_path / "b") == []        # a relocated state file is not a difference
    path = tmp_path / "b" / "run" / "summary.json"
    doc = json.loads(path.read_text())
    doc["resumed_from"]["state"] = "/y/saved/other.state.pkl"                # another file name is
    path.write_text(json.dumps(doc))
    assert len(same_results.compare(tmp_path / "a", tmp_path / "b")) == 1


def test_run_resolves_relative_out_against_the_callers_folder(tmp_path, monkeypatch):
    code = tmp_path / "code"
    code.mkdir()
    calls = []

    class Done:
        returncode = 0

    monkeypatch.setattr(same_results.subprocess, "run",
                        lambda command, **kw: calls.append((command, kw.get("cwd"))) or Done())
    monkeypatch.setattr(same_results, "power_line", lambda: "AC Power")
    monkeypatch.chdir(tmp_path)
    assert same_results.main(["run", "--code", str(code), "--out", "rel_out", "--cases", "nomelt"]) == 0
    command, cwd = calls[0]
    assert Path(cwd) == code
    outdir = Path(command[command.index("--outdir") + 1])
    assert outdir == tmp_path / "rel_out" / "nomelt"                          # not relative to --code
    assert (tmp_path / "rel_out" / "nomelt.log").exists()
