"""Generate sub-plan 21 (docs/superpowers/plans/melt-spraying-subplans/21-column-flow.md) from its prose, part_21.md,
and the tested code in code/, and check that the blocks rebuild the tested files from the base exactly.

The base is prototype/work-2026-10-12-skin-heat-plan/code, sub-plan 20's tested copy, as committed: its tracked files
must be unchanged in the working tree. Every {{NAME}} in part_21.md becomes a fenced block cut from the tested files: a
task's diff of a file runs from the state the tasks before it left (the base for the first) to its own end state, kept
in states/task<N>/ for every file a later task changes again (the build repository's commit of that task) and read
from code/ otherwise. The checks then replay the tasks in order on a copy of the base -- each diff applied with `patch`
-- and compare every file byte for byte with code/. Tasks are added here as they are written at code level. Run:

    "$PY" make_plan21.py            # writes the sub-plan and prints the checks
"""
import difflib
import filecmp
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.join(HERE, "code")
BASE = os.path.join(os.path.dirname(HERE), "work-2026-10-12-skin-heat-plan", "code")
OUT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "docs", "superpowers", "plans", "melt-spraying-subplans",
                   "21-column-flow.md")
STATES = os.path.join(HERE, "states")
SKIP = ("__pycache__", ".pytest_cache", "reentry_model_output", ".git")
# The tasks written at code level so far, in order, and the files each changes ({{key}} in part_21.md, path)
TASKS = (
    (1, (("DIFF_COLUMNS_T1", "reentry_model/columns.py"), ("DIFF_SKIN_T1", "reentry_model/skin.py"),
         ("DIFF_SKINBODY_T1", "reentry_model/skinbody.py"), ("DIFF_TEST_SKIN_T1", "tests/test_reentry_model_skin.py"))),
    (2, (("DIFF_COLUMNS_T2", "reentry_model/columns.py"), ("DIFF_TEST_COLUMNS_T2", "tests/test_reentry_model_columns.py"),
         ("DIFF_SKIN_T2", "reentry_model/skin.py"), ("DIFF_TEST_SKIN_T2", "tests/test_reentry_model_skin.py"))),
    (3, (("RHEOLOGY_T3", "reentry_model/rheology.py"), ("TEST_RHEOLOGY_T3", "tests/test_reentry_model_rheology.py"),
         ("DIFF_SKIN_T3", "reentry_model/skin.py"), ("DIFF_TEST_SKIN_T3", "tests/test_reentry_model_skin.py"))),
    (4, (("DIFF_SKIN_T4", "reentry_model/skin.py"), ("DIFF_TEST_SKIN_T4", "tests/test_reentry_model_skin.py"))),
    (5, (("DIFF_SURFACE_FLOW_T5", "reentry_model/surface_flow.py"),
         ("DIFF_TEST_SURFACE_FLOW_T5", "tests/test_reentry_model_surface_flow.py"),
         ("DIFF_SKIN_T5", "reentry_model/skin.py"), ("DIFF_TEST_SKIN_T5", "tests/test_reentry_model_skin.py"))),
    ("6a", (("DIFF_TEST_COLUMNS_T6A", "tests/test_reentry_model_columns.py"),
            ("DIFF_TEST_SKIN_T6A", "tests/test_reentry_model_skin.py"),
            ("DIFF_COLUMNS_T6A", "reentry_model/columns.py"), ("DIFF_SKIN_T6A", "reentry_model/skin.py"))),
    ("6b", (("DIFF_TEST_COLUMNS_T6B", "tests/test_reentry_model_columns.py"),
            ("DIFF_TEST_SKIN_T6B", "tests/test_reentry_model_skin.py"),
            ("DIFF_COLUMNS_T6B", "reentry_model/columns.py"), ("DIFF_SKIN_T6B", "reentry_model/skin.py"))),
    ("6c", (("DIFF_TEST_SKIN_T6C", "tests/test_reentry_model_skin.py"),
            ("DIFF_TEST_SKINBODY_T6C", "tests/test_reentry_model_skinbody.py"),
            ("DIFF_TEST_CLI_T6C", "tests/test_reentry_model_cli.py"),
            ("DIFF_SKIN_T6C", "reentry_model/skin.py"), ("DIFF_SKINBODY_T6C", "reentry_model/skinbody.py"),
            ("DIFF_COUPLED_T6C", "reentry_model/coupled.py"), ("DIFF_CLI_T6C", "reentry_model/cli.py"))),
    ("6d", (("DIFF_TEST_SKIN_T6D", "tests/test_reentry_model_skin.py"),
            ("DIFF_TEST_SKINBODY_T6D", "tests/test_reentry_model_skinbody.py"),
            ("DIFF_TEST_CLI_T6D", "tests/test_reentry_model_cli.py"),
            ("DIFF_SKIN_T6D", "reentry_model/skin.py"), ("DIFF_SKINBODY_T6D", "reentry_model/skinbody.py"),
            ("DIFF_COUPLED_T6D", "reentry_model/coupled.py"), ("DIFF_CLI_T6D", "reentry_model/cli.py"))),
    (7, (("DIFF_TEST_SKIN_T7", "tests/test_reentry_model_skin.py"),
         ("DIFF_TEST_SKINBODY_T7", "tests/test_reentry_model_skinbody.py"),
         ("DIFF_SKIN_T7", "reentry_model/skin.py"), ("DIFF_SKINBODY_T7", "reentry_model/skinbody.py"),
         ("DIFF_COUPLED_T7", "reentry_model/coupled.py"))),
)


def read(root, rel):
    with open(os.path.join(root, rel)) as fh:
        return fh.read()


def write(root, rel, text):
    path = os.path.join(root, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def fence(text, lang="python"):
    return "```{}\n{}```".format(lang, text if text.endswith("\n") else text + "\n")


def unified(a, b, rel):
    return "".join(difflib.unified_diff(a.splitlines(keepends=True), b.splitlines(keepends=True), "a/" + rel, "b/" + rel))


def files(root):
    out = set()
    for d, dirs, names in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP]
        out.update(os.path.relpath(os.path.join(d, f), root) for f in names)
    return out


def main():
    dirty = subprocess.run(["git", "-C", BASE, "status", "--porcelain", "--untracked-files=no", "--", "."],
                           capture_output=True, text=True).stdout.strip()
    if dirty:
        raise SystemExit("the base must be as committed (dirty: {!r})".format(dirty))
    blocks, steps = {}, []                # steps: (kind, rel, payload) replayed on a copy of the base, in task order

    current = {}                          # each file's text as the tasks so far leave it
    for task, changes in TASKS:
        for key, rel in changes:
            snapshot = os.path.join(STATES, "task{}".format(task), rel)
            after = read(STATES, os.path.join("task{}".format(task), rel)) if os.path.exists(snapshot) else read(NEW, rel)
            if rel not in current and not os.path.exists(os.path.join(BASE, rel)):     # a new file: written whole
                blocks[key] = fence(after)
                steps.append(("write", rel, after))
                current[rel] = after
                continue
            before = current[rel] if rel in current else read(BASE, rel)
            d = unified(before, after, rel)
            if not d:
                raise SystemExit("empty diff for {}".format(key))
            blocks[key] = fence(d, "diff")
            steps.append(("patch", rel, d))
            current[rel] = after
    kept = files(STATES)
    named = {os.path.join("task{}".format(task), rel) for task, changes in TASKS for _, rel in changes}
    if kept - named:
        raise SystemExit("states/ holds files no task changes: {}".format(sorted(kept - named)))

    # -- replay on a copy of the base and compare with code/ ------------------------------------------------------------
    tmp = tempfile.mkdtemp()
    try:
        work = os.path.join(tmp, "code")
        shutil.copytree(BASE, work, ignore=shutil.ignore_patterns(*SKIP))
        for kind, rel, payload in steps:
            if kind == "patch":
                done = subprocess.run(["patch", "-p1", "-s", "-d", work], input=payload, text=True, capture_output=True)
                if done.returncode:
                    raise SystemExit("a diff for {} does not apply: {}".format(rel, done.stderr or done.stdout))
            else:
                write(work, rel, payload)
        differ = sorted(rel for rel in files(NEW) | files(work)
                        if not (os.path.exists(os.path.join(NEW, rel)) and os.path.exists(os.path.join(work, rel))
                                and filecmp.cmp(os.path.join(NEW, rel), os.path.join(work, rel), shallow=False)))
        if differ:
            raise SystemExit("the replayed tasks do not reproduce code/: {}".format(differ))
    finally:
        shutil.rmtree(tmp)

    text = read(HERE, "part_21.md")
    for key, block in blocks.items():
        marker = "{{" + key + "}}"
        if text.count(marker) != 1:
            raise SystemExit("part_21.md must use {} exactly once".format(marker))
        text = text.replace(marker, block)
    if "{{" in text:
        raise SystemExit("unfilled marker in part_21.md")
    with open(OUT, "w") as fh:
        fh.write(text)
    print("checks passed: the {} steps written so far, replayed on sub-plan 20's tested copy, reproduce every file of "
          "code/ byte for byte".format(len(steps)))
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
