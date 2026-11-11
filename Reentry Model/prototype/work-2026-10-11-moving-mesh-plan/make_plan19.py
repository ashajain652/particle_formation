"""Generate sub-plan 19 (docs/superpowers/plans/melt-spraying-subplans/19-moving-mesh.md) from its prose, part_19.md,
and the tested code in code/, and check that the blocks reproduce the tested files from the base exactly.

The base is prototype/work-2026-10-08-skin/code at 49d88e0 (its working tree, which must be clean). Every {{NAME}} in
part_19.md becomes a fenced block cut from the tested files; the checks below then rebuild each changed file from the
base and the blocks, in the order the tasks apply them, and compare byte for byte. Run:

    "$PY" make_plan19.py            # writes the sub-plan and prints the checks
"""
import ast
import difflib
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.join(HERE, "code")
BASE = os.path.join(os.path.dirname(HERE), "work-2026-10-08-skin", "code")
OUT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "docs", "superpowers", "plans", "melt-spraying-subplans",
                   "19-moving-mesh.md")
BASE_HEAD = "49d88e0"
EXPECTED_CHANGES = {"reentry_model/movemesh.py", "reentry_model/thermal/skfem_backend.py",
                    "reentry_model/thermal/fenicsx_backend.py", "tests/test_reentry_model_movemesh.py",
                    "tests/test_reentry_model_fenicsx.py", "analysis/movemesh_recession.py"}


def read(root, rel):
    with open(os.path.join(root, rel)) as fh:
        return fh.read()


def fence(text, lang="python"):
    return "```{}\n{}```".format(lang, text if text.endswith("\n") else text + "\n")


def cut(lines, starts):
    """Contiguous slices of `lines` beginning at each index in `starts` (sorted), the last running to the end."""
    starts = list(starts) + [len(lines)]
    return ["".join(lines[a:b]) for a, b in zip(starts[:-1], starts[1:])]


def line_starting(lines, prefix):
    hits = [i for i, line in enumerate(lines) if line.startswith(prefix)]
    if len(hits) != 1:
        raise SystemExit("expected one line starting {!r}, found {}".format(prefix, len(hits)))
    return hits[0]


def main():
    head = subprocess.run(["git", "-C", BASE, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", BASE, "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                           text=True).stdout.strip()
    if head != BASE_HEAD or dirty:
        raise SystemExit("the base must be {} with a clean tree (found {}, dirty: {!r})".format(BASE_HEAD, head, dirty))
    changed = set()
    for root, dirs, files in os.walk(NEW):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", ".pytest_cache", "reentry_model_output")]
        for f in files:
            rel = os.path.relpath(os.path.join(root, f), NEW)
            if not os.path.exists(os.path.join(BASE, rel)) or read(BASE, rel) != read(NEW, rel):
                changed.add(rel)
    if changed != EXPECTED_CHANGES:
        raise SystemExit("files changed against the base: {}".format(sorted(changed ^ EXPECTED_CHANGES)))
    blocks = {}
    # movemesh.py: the slices the tasks append or insert, in file order T1 T3 T4 T5A T6 T5B
    mm = read(NEW, "reentry_model/movemesh.py")
    L = mm.splitlines(keepends=True)
    starts = [0, line_starting(L, "# -- the deforming zone"), line_starting(L, "# -- moving the mesh"),
              line_starting(L, "    def rebuild("), line_starting(L, "    def needs_rebuild("),
              line_starting(L, "# -- rebuilding the mesh")]
    t1, t3, t4, t5a, t6, t5b = cut(L, starts)
    assert t1 + t3 + t4 + t5a + t6 + t5b == mm
    blocks.update(MM_T1=fence(t1), MM_T3=fence(t3), MM_T4=fence(t4), MM_T5A=fence(t5a), MM_T6=fence(t6), MM_T5B=fence(t5b))
    # the new test module: one section per task, marked "# -- ... (Task N)"
    tm = read(NEW, "tests/test_reentry_model_movemesh.py")
    L = tm.splitlines(keepends=True)
    marks = [line_starting(L, "# -- the backends on moved geometry (Task 2)"),
             line_starting(L, "# -- the deforming zone (Task 3)"),
             line_starting(L, "# -- moving the mesh: the remap and its books (Task 4)"),
             line_starting(L, "# -- the rebuild (Task 5)"), line_starting(L, "# -- the rebuild policy (Task 6)"),
             line_starting(L, "# -- a moving front against an exact solution (Task 7)")]
    sections = cut(L, [0] + marks)
    assert "".join(sections) == tm
    for k, text in enumerate(sections, start=1):
        blocks["TEST_T{}".format(k)] = fence(text)
    # the FEniCSx test module: three functions appended to the base
    fx_base, fx_new = read(BASE, "tests/test_reentry_model_fenicsx.py"), read(NEW, "tests/test_reentry_model_fenicsx.py")
    assert fx_new.startswith(fx_base)
    tail = fx_new[len(fx_base):]
    TL = tail.splitlines(keepends=True)
    names = ["test_moved_geometry_steps_like_a_fresh_solver", "test_a_rebuild_keeps_the_energy_on_fenicsx",
             "test_a_moving_front_matches_the_skfem_backend"]
    fx_parts = cut(TL, [0] + [line_starting(TL, "def " + n) for n in names[1:]])
    assert "".join(fx_parts) == tail and all(n in p for n, p in zip(names, fx_parts))
    blocks.update(FX_T2=fence(fx_parts[0]), FX_T5=fence(fx_parts[1]), FX_T7=fence(fx_parts[2]))
    # the backends: unified diffs against the base, checked by applying them with patch
    tmp = tempfile.mkdtemp()
    try:
        for key, rel in (("DIFF_SKFEM", "reentry_model/thermal/skfem_backend.py"),
                         ("DIFF_FENICSX", "reentry_model/thermal/fenicsx_backend.py")):
            a, b = read(BASE, rel), read(NEW, rel)
            diff = "".join(difflib.unified_diff(a.splitlines(keepends=True), b.splitlines(keepends=True),
                                                "a/" + rel, "b/" + rel))
            os.makedirs(os.path.join(tmp, os.path.dirname(rel)), exist_ok=True)
            shutil.copy(os.path.join(BASE, rel), os.path.join(tmp, rel))
            done = subprocess.run(["patch", "-p1", "-s", "-d", tmp], input=diff, text=True, capture_output=True)
            if done.returncode or read(tmp, rel) != b:
                raise SystemExit("the {} diff does not reproduce the tested file: {}".format(rel, done.stderr))
            blocks[key] = fence(diff, "diff")
    finally:
        shutil.rmtree(tmp)
    blocks["HARNESS"] = fence(read(NEW, "analysis/movemesh_recession.py"))
    # assemble
    text = read(HERE, "part_19.md")
    for key, block in blocks.items():
        marker = "{{" + key + "}}"
        if text.count(marker) != 1:
            raise SystemExit("part_19.md must use {} exactly once".format(marker))
        text = text.replace(marker, block)
    if "{{" in text:
        raise SystemExit("unfilled marker in part_19.md")
    with open(OUT, "w") as fh:
        fh.write(text)
    print("checks passed: movemesh.py, the test module and the FEniCSx tests reassemble exactly; both diffs apply; "
          "only the six expected files differ from {}".format(BASE_HEAD))
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
