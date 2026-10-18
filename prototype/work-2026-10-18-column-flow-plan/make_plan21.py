"""Generate sub-plan 21 (docs/superpowers/plans/melt-spraying-subplans/21-column-flow.md) from its prose, part_21.md,
and the tested code in code/, and check that the blocks rebuild the tested files from the base exactly.

The base is prototype/work-2026-10-12-skin-heat-plan/code, sub-plan 20's tested copy, as committed: its tracked files
must be unchanged in the working tree. Every {{NAME}} in part_21.md becomes a fenced block cut from the tested files. The
checks then replay the tasks in order on a copy of the base -- each diff applied with `patch`, each new file written as
the plan writes it -- and compare every file that differs from the base byte for byte with code/. Tasks are added here as
they are written at code level. Run:

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
SKIP = ("__pycache__", ".pytest_cache", "reentry_model_output", ".git")


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

    def diff_step(key, rel):
        d = unified(read(BASE, rel), read(NEW, rel), rel)
        if not d:
            raise SystemExit("empty diff for {}".format(key))
        blocks[key] = fence(d, "diff")
        steps.append(("patch", rel, d))

    # -- Task 1: the implicit sideways conduction (design (a), Asha 2026-10-18) -----------------------------------------
    for key, rel in (("DIFF_COLUMNS_T1", "reentry_model/columns.py"), ("DIFF_SKIN_T1", "reentry_model/skin.py"),
                     ("DIFF_SKINBODY_T1", "reentry_model/skinbody.py"),
                     ("DIFF_TEST_SKIN_T1", "tests/test_reentry_model_skin.py")):
        diff_step(key, rel)

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
            elif kind == "write":
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
