"""Generate sub-plan 20 (docs/superpowers/plans/melt-spraying-subplans/20-layered-skin-heat.md) from its prose, part_20.md,
and the tested code in code/, and check that the blocks rebuild the tested files from the base exactly.

The base is prototype/work-2026-10-12-moving-mesh/code at 70ecbdd (its working tree, which must be clean). Every {{NAME}}
in part_20.md becomes a fenced block cut from the tested files. The checks then replay the tasks in order on a copy of
the base -- each diff applied with `patch`, each new file or section written as the plan writes it, the renames and the
reference extraction done as the plan does them -- and compare every file that differs from the base byte for byte
with code/. Files that several tasks change (body.py, cli.py, two test modules) pass through the intermediate states
the earlier tasks leave, built here by the replacements those tasks make. Run:

    "$PY" make_plan20.py            # writes the sub-plan and prints the checks
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
BASE = os.path.join(os.path.dirname(HERE), "work-2026-10-12-moving-mesh", "code")
OUT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "docs", "superpowers", "plans", "melt-spraying-subplans",
                   "20-layered-skin-heat.md")
BASE_HEAD = "70ecbdd"
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


def replace(text, pairs):
    for a, b in pairs:
        if text.count(a) < 1:
            raise SystemExit("an intermediate state's replacement does not apply: {!r}".format(a[:70]))
        text = text.replace(a, b)
    return text


def unified(a, b, rel_a, rel_b=None):
    return "".join(difflib.unified_diff(a.splitlines(keepends=True), b.splitlines(keepends=True), "a/" + rel_a,
                                        "b/" + (rel_b or rel_a)))


def sections(text, markers):
    """`text` cut at the lines starting with each of `markers` (in order): len(markers) + 1 pieces that rejoin to it."""
    lines = text.splitlines(keepends=True)
    starts = [0]
    for m in markers:
        hits = [i for i, line in enumerate(lines) if line.startswith(m)]
        if len(hits) != 1:
            raise SystemExit("expected one line starting {!r}, found {}".format(m, len(hits)))
        starts.append(hits[0])
    if starts != sorted(starts):
        raise SystemExit("markers out of order: {}".format(markers))
    starts.append(len(lines))
    parts = ["".join(lines[a:b]) for a, b in zip(starts[:-1], starts[1:])]
    assert "".join(parts) == text
    return parts


def main():
    head = subprocess.run(["git", "-C", BASE, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", BASE, "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                           text=True).stdout.strip()
    if head != BASE_HEAD or dirty:
        raise SystemExit("the base must be {} with a clean tree (found {}, dirty: {!r})".format(BASE_HEAD, head, dirty))
    blocks, steps = {}, []                # steps: (kind, rel, payload) replayed on a copy of the base, in task order

    def diff_step(key, rel, before, after, rel_after=None):
        d = unified(before, after, rel, rel_after)
        if not d:
            raise SystemExit("empty diff for {}".format(key))
        blocks[key] = fence(d, "diff")
        steps.append(("patch", rel_after or rel, d))

    # -- Task 0: the moving mesh's books ------------------------------------------------------------------------------
    for key, rel in (("DIFF_MOVEMESH", "reentry_model/movemesh.py"), ("DIFF_TEST_MOVEMESH", "tests/test_reentry_model_movemesh.py"),
                     ("DIFF_RECESSION", "analysis/movemesh_recession.py")):
        diff_step(key, rel, read(BASE, rel), read(NEW, rel))
    # -- Task 1: the columns ------------------------------------------------------------------------------------------
    for key, rel in (("COLUMNS", "reentry_model/columns.py"), ("TEST_COLUMNS", "tests/test_reentry_model_columns.py")):
        blocks[key] = fence(read(NEW, rel))
        steps.append(("write", rel, read(NEW, rel)))
    # -- Task 2: stage 1 renamed, the new skin's cells ----------------------------------------------------------------
    steps.append(("move", "reentry_model/skin.py", "reentry_model/skin1d.py"))
    steps.append(("move", "tests/test_reentry_model_skin.py", "tests/test_reentry_model_skin1d.py"))
    diff_step("DIFF_SKIN1D", "reentry_model/skin1d.py", read(BASE, "reentry_model/skin.py"), read(NEW, "reentry_model/skin1d.py"))
    diff_step("DIFF_TEST_SKIN1D", "tests/test_reentry_model_skin1d.py", read(BASE, "tests/test_reentry_model_skin.py"),
              read(NEW, "tests/test_reentry_model_skin1d.py"))
    body_a = replace(read(BASE, "reentry_model/body.py"), [
        ("from . import skin as skin_mod\n", "from . import skin1d as skin_mod\n"),
        ("        from .skin import MAX_THICKNESS_FACTOR\n", "        from .skin1d import MAX_THICKNESS_FACTOR\n")])
    cli_a = replace(read(BASE, "reentry_model/cli.py"), [("sesam_io, skin, spray, surface_flow", "sesam_io, skin1d as skin, spray, surface_flow")])
    melt_a = replace(read(BASE, "tests/test_reentry_model_melting_skin.py"), [
        ("from reentry_model import body, heating, material, mesh, skin, thermal\n",
         "from reentry_model import body, heating, material, mesh, thermal\nfrom reentry_model import skin1d as skin\n")])
    fx_base = read(BASE, "tests/test_reentry_model_fenicsx.py")
    fx_a = replace(fx_base, [("    from reentry_model import body, heating, mesh, skin\n",
                              "    from reentry_model import body, heating, mesh\n    from reentry_model import skin1d as skin\n")])
    diff_step("DIFF_BODY_A", "reentry_model/body.py", read(BASE, "reentry_model/body.py"), body_a)
    diff_step("DIFF_CLI_A", "reentry_model/cli.py", read(BASE, "reentry_model/cli.py"), cli_a)
    diff_step("DIFF_MELTSKIN_A", "tests/test_reentry_model_melting_skin.py", read(BASE, "tests/test_reentry_model_melting_skin.py"), melt_a)
    diff_step("DIFF_FENICSX_A", "tests/test_reentry_model_fenicsx.py", fx_base, fx_a)
    sk = sections(read(NEW, "reentry_model/skin.py"), ["# -- the minimum mass (Task 3)", "# -- heat through the depth (Task 4)",
                                                       "# -- heat sideways (Task 5)", "# -- the base (Task 7)"])
    ts = sections(read(NEW, "tests/test_reentry_model_skin.py"), ["# -- the minimum mass (Task 3)",
                                                                   "# -- heat through the depth (Task 4)",
                                                                   "# -- heat sideways (Task 5)"])
    blocks["SKIN_T2"], blocks["TEST_SKIN_T2"] = fence(sk[0]), fence(ts[0])
    steps.append(("write", "reentry_model/skin.py", sk[0]))
    steps.append(("write", "tests/test_reentry_model_skin.py", ts[0]))
    # -- Tasks 3-5: appended to the new skin and its tests ------------------------------------------------------------
    for task, (s_part, t_part) in zip((3, 4, 5), zip(sk[1:4], ts[1:4])):
        blocks["SKIN_T{}".format(task)], blocks["TEST_SKIN_T{}".format(task)] = fence(s_part), fence(t_part)
        steps.append(("append", "reentry_model/skin.py", s_part))
        steps.append(("append", "tests/test_reentry_model_skin.py", t_part))
    # -- Task 6: the nodal base law -----------------------------------------------------------------------------------
    for key, rel in (("DIFF_SKFEM", "reentry_model/thermal/skfem_backend.py"),
                     ("DIFF_FENICSX_BACKEND", "reentry_model/thermal/fenicsx_backend.py"),
                     ("DIFF_THERMAL_INIT", "reentry_model/thermal/__init__.py"),
                     ("DIFF_TEST_THERMAL", "tests/test_reentry_model_thermal.py")):
        diff_step(key, rel, read(BASE, rel), read(NEW, rel))
    fx_new = read(NEW, "tests/test_reentry_model_fenicsx.py")
    i6 = fx_new.index("\n\ndef test_backends_agree_with_a_nodal_law")
    i7 = fx_new.index("\n\ndef test_the_skin_body_matches_the_skfem_backend")
    fx_t6, fx_t7 = fx_new[i6:i7], fx_new[i7:]
    blocks["FX_T6"], blocks["FX_T7"] = fence(fx_t6.lstrip("\n")), fence(fx_t7.lstrip("\n"))
    steps.append(("append", "tests/test_reentry_model_fenicsx.py", fx_t6))
    # -- Task 7: the body -----------------------------------------------------------------------------------------------
    blocks["SKIN_T7"] = fence(sk[4])
    steps.append(("append", "reentry_model/skin.py", sk[4]))
    for key, rel in (("SKINBODY", "reentry_model/skinbody.py"), ("TEST_SKINBODY", "tests/test_reentry_model_skinbody.py")):
        blocks[key] = fence(read(NEW, rel))
        steps.append(("write", rel, read(NEW, rel)))
    steps.append(("append", "tests/test_reentry_model_fenicsx.py", fx_t7))
    # -- Task 8: the command line, the coupled run's outputs -----------------------------------------------------------
    diff_step("DIFF_BODY_B", "reentry_model/body.py", body_a, read(NEW, "reentry_model/body.py"))
    diff_step("DIFF_CLI_B", "reentry_model/cli.py", cli_a, read(NEW, "reentry_model/cli.py"))
    for key, rel in (("DIFF_COUPLED", "reentry_model/coupled.py"), ("DIFF_TEST_CLI", "tests/test_reentry_model_cli.py")):
        diff_step(key, rel, read(BASE, rel), read(NEW, rel))
    diff_step("DIFF_MELTSKIN_B", "tests/test_reentry_model_melting_skin.py", melt_a, read(NEW, "tests/test_reentry_model_melting_skin.py"))
    fx_c = fx_a + fx_t6 + fx_t7
    diff_step("DIFF_FENICSX_D", "tests/test_reentry_model_fenicsx.py", fx_c, fx_new)
    # -- Tasks 9 and 10: harnesses and the reference ----------------------------------------------------------------
    for key, rel in (("ABLATION", "analysis/ablation_1d.py"), ("PLATE", "analysis/skin_plate.py"),
                     ("TIMING", "analysis/skin_timing.py")):
        blocks[key] = fence(read(NEW, rel))
        steps.append(("write", rel, read(NEW, rel)))
    rel = "data/reference_values/plate_hot_1d.json"                       # the hot plate's reference (2026-10-14)
    blocks["PLATE_HOT_REF"] = fence(read(NEW, rel), "json")
    steps.append(("write", rel, read(NEW, rel)))
    blocks["EXTRACT"] = fence(read(HERE, "extract_fine_box.py"))
    steps.append(("extract", "data/reference_values/box_q6p_fine.json", None))

    # -- replay on a copy of the base and compare with code/ ----------------------------------------------------------
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
            elif kind == "append":
                write(work, rel, read(work, rel) + payload)
            elif kind == "move":
                os.rename(os.path.join(work, rel), os.path.join(work, payload))
            elif kind == "extract":
                run = os.path.join(os.path.dirname(HERE), "work-2026-10-08-skin", "analysis", "box", "runs", "q6p__fine.json")
                if not os.path.exists(run):          # the run outputs are git-ignored: in a worktree, use the main checkout's
                    main = subprocess.run(["git", "-C", HERE, "worktree", "list", "--porcelain"], capture_output=True,
                                          text=True).stdout.split("\n")[0].split(" ", 1)[1]
                    run = os.path.join(main, "prototype", "work-2026-10-08-skin", "analysis", "box", "runs", "q6p__fine.json")
                subprocess.run([sys.executable, os.path.join(HERE, "extract_fine_box.py"), run, os.path.join(work, rel)], check=True)
        differ = []
        for root, dirs, files in os.walk(NEW):
            dirs[:] = [d for d in dirs if d not in SKIP]
            for f in files:
                rel = os.path.relpath(os.path.join(root, f), NEW)
                other = os.path.join(work, rel)
                if not os.path.exists(other) or not filecmp.cmp(os.path.join(NEW, rel), other, shallow=False):
                    differ.append(rel)
        for root, dirs, files in os.walk(work):
            dirs[:] = [d for d in dirs if d not in SKIP]
            for f in files:
                rel = os.path.relpath(os.path.join(root, f), work)
                if not os.path.exists(os.path.join(NEW, rel)):
                    differ.append(rel + " (only in the replay)")
        if differ:
            raise SystemExit("the replayed tasks do not reproduce code/: {}".format(sorted(differ)))
    finally:
        shutil.rmtree(tmp)

    text = read(HERE, "part_20.md")
    for key, block in blocks.items():
        marker = "{{" + key + "}}"
        if text.count(marker) != 1:
            raise SystemExit("part_20.md must use {} exactly once".format(marker))
        text = text.replace(marker, block)
    if "{{" in text:
        raise SystemExit("unfilled marker in part_20.md")
    with open(OUT, "w") as fh:
        fh.write(text)
    print("checks passed: the {} steps of Tasks 0-10, replayed on {}, reproduce every file of code/ byte for byte".format(
        len(steps), BASE_HEAD))
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
