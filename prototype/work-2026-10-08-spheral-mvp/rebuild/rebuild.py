#!/usr/bin/env python3
"""Rebuild prototype/proto3 -- the Step 3 prototype -- mechanically from the committed documents.

    python3 prototype/rebuild/rebuild.py [--out prototype/proto3] [--upto STAGE] [--post] [--force]

The prototype and its plan generator were lost on 2026-10-06/07 (not version-controlled, no backup). Commit de3c22e
records the recovery route: every prototype file appears in the committed plan
docs/superpowers/plans/2026-09-20-melt-spraying.md as fenced blocks ("Replace", "Create", "Append to", targeted
"replace ... with ..." edits), and the later amendments live in docs/superpowers/plans/melt-spraying-subplans/ as
unified diffs against the prototype. This script does exactly that, in order, and fails loudly if any block is not
where it expects it (each block is addressed by the line its fence opens on AND by a sentinel in the prose before it):

  stage base      `git archive ee70f6b` of reentry_model/ tests/ data/ analysis/ (the Step 2 state the plan was
                  written against: reentry_model/ and tests/ are unchanged from 07b2856/e376a82 to 187cd51)
  stage task1     Task 1 is superseded by sub-plan 01, implemented in the prototype and ported to main in 187cd51:
                  mesh.py and test_reentry_model_mesh.py from 187cd51 (only those two files: sub-plan 01's plan
                  edits nothing else in the prototype; 187cd51's band = 0 pins in cli.py/conftest.py/the reference
                  test were made for the Step 2 package only)
  stage plan      Tasks 2-15 of the master plan (code blocks verbatim; Task 11's SESAM references and Task 9's DRAMA
                  extraction cannot be run -- see README)
  stage fact37    measured fact 37 (2026-09-27): MeltingBody.last_face_ids / on_current_surface and the VTK writer
                  reading through it (sub-plans 09 and 10), placed where the 2026-10-02 diffs' context puts them
  stage d1002     the nine diff blocks of the amendment of 2026-10-02 (deep runoff, delta_m and deep_thickness frames)
  stage d1003     the eight diff blocks of 2026-10-03 (molten cascade)
  stage d1005     the three diff blocks of 2026-10-05 (--seed)
  stage post      (only with --post) post-reconstruction additions, NOT part of the recovered state:
                  01 sub-plan 02's Scheil and empirical-data materials (post_subplan02), then
                  prototype/rebuild/post/*.diff in name order (02: p_w/tau/closure on every evaluated step, v_hat in the JSON;
                  03: the material in the run name; 04: n_derived in the frames; 05: the runoff's emptying-number bound;
                  06: the gas-side frame fields on faces the step's deaths exposed)

Every operation is logged to <out>/../rebuild/rebuild_log.txt (with patch's own hunk report, so offsets show).
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tarfile
import io

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
PLAN = os.path.join(REPO, "docs", "superpowers", "plans", "2026-09-20-melt-spraying.md")
SUB = os.path.join(REPO, "docs", "superpowers", "plans", "melt-spraying-subplans")
BASE_COMMIT = "ee70f6b"          # "Add the Step 3 implementation plan": Step 2 code, the plan's base
TASK1_COMMIT = "187cd51"         # sub-plan 01 ported to main (mesh.py + its tests)
PY = "/Users/ashajain/miniforge3/envs/drama_env/bin/python"
STAGES = ["base", "task1", "plan", "band13", "fact37", "d1002", "d1003", "d1005", "post"]

LOG = []


def log(msg):
    LOG.append(msg)
    print(msg)


# ------------------------------------------------------------------------------------------------- block parsing
def fenced_blocks(path):
    """{opening line number (1-based): (lang, text)} of every fenced block. A ```markdown block may contain inner
    fences (Task 15), which are tracked by depth so that they do not close it."""
    lines = open(path, encoding="utf-8").read().split("\n")
    out, i = {}, 0
    while i < len(lines):
        m = re.match(r"^```(\S*)\s*$", lines[i])
        if not m:
            i += 1
            continue
        lang, start, depth, body = m.group(1), i + 1, 0, []
        i += 1
        while i < len(lines):
            line = lines[i]
            if lang == "markdown" and re.match(r"^```\S+", line):
                depth += 1
            elif line.rstrip() == "```":
                if depth == 0:
                    break
                depth -= 1
            body.append(line)
            i += 1
        out[start] = (lang, "\n".join(body) + "\n")
        i += 1
    return out, lines


class Doc:
    def __init__(self, path):
        self.path = path
        self.blocks, self.lines = fenced_blocks(path)

    def block(self, line, sentinel=None, lang=None):
        if line not in self.blocks:
            raise SystemExit("{}: no fenced block opens on line {}".format(os.path.basename(self.path), line))
        blang, text = self.blocks[line]
        if lang is not None and blang != lang:
            raise SystemExit("{}:{}: block language {!r}, expected {!r}".format(self.path, line, blang, lang))
        if sentinel is not None:
            prose = "\n".join(self.lines[max(0, line - 12):line - 1])
            if sentinel not in prose:
                raise SystemExit("{}:{}: sentinel {!r} not found in the prose before the block".format(self.path, line, sentinel))
        return text


# ------------------------------------------------------------------------------------------------- file operations
class Tree:
    def __init__(self, root):
        self.root = root

    def p(self, rel):
        return os.path.join(self.root, rel)

    def read(self, rel):
        with open(self.p(rel), encoding="utf-8") as fh:
            return fh.read()

    def write(self, rel, text):
        os.makedirs(os.path.dirname(self.p(rel)), exist_ok=True)
        with open(self.p(rel), "w", encoding="utf-8") as fh:
            fh.write(text)

    def replace_file(self, rel, text, what="replace"):
        if what == "create" and os.path.exists(self.p(rel)):
            raise SystemExit("create {}: already exists".format(rel))
        if what == "replace" and not os.path.exists(self.p(rel)):
            raise SystemExit("replace {}: does not exist".format(rel))
        self.write(rel, text)
        log("  {} {} ({} lines)".format(what, rel, text.count("\n")))

    def append(self, rel, text):
        """'Append': the block goes at the end of the existing file, after the two blank lines PEP 8 puts before a
        top-level block (the line numbers of the 2026-10-02/03 test diffs confirm this separator)."""
        old = self.read(rel)
        new = old.rstrip("\n") + "\n\n\n" + text
        self.write(rel, new)
        log("  append {} (+{} lines, now {})".format(rel, text.count("\n"), new.count("\n")))

    def edit(self, rel, old, new, label=""):
        text = self.read(rel)
        n = text.count(old)
        if n != 1:
            raise SystemExit("edit {} {}: the old text occurs {} times".format(rel, label, n))
        self.write(rel, text.replace(old, new))
        log("  edit {} {}".format(rel, label))

    def insert_after_line(self, rel, anchor_line, new, label=""):
        """Insert `new` (ending in a newline) after the unique line that starts (after indentation) with anchor_line."""
        text = self.read(rel)
        lines = text.split("\n")
        hits = [i for i, l in enumerate(lines) if l.strip().startswith(anchor_line.strip())]
        if len(hits) != 1:
            raise SystemExit("insert {} {}: anchor found {} times".format(rel, label, len(hits)))
        i = hits[0]
        out = "\n".join(lines[:i + 1]) + "\n" + new + "\n".join(lines[i + 1:])
        self.write(rel, out)
        log("  insert after line {} of {} {}".format(i + 1, rel, label))

    def insert_before_line(self, rel, anchor_line, new, label=""):
        text = self.read(rel)
        lines = text.split("\n")
        hits = [i for i, l in enumerate(lines) if l.strip().startswith(anchor_line.strip())]
        if len(hits) != 1:
            raise SystemExit("insert {} {}: anchor found {} times".format(rel, label, len(hits)))
        i = hits[0]
        out = "\n".join(lines[:i]) + "\n" + new + "\n".join(lines[i:])
        self.write(rel, out)
        log("  insert before line {} of {} {}".format(i + 1, rel, label))

    def replace_def(self, rel, header, new, indent=""):
        """Replace a whole top-level (or method, by `indent`) definition starting with `header` up to the next line at
        the same or lower indentation that is not blank."""
        lines = self.read(rel).split("\n")
        hits = [i for i, l in enumerate(lines) if l.startswith(indent + header)]
        if len(hits) != 1:
            raise SystemExit("replace_def {} {!r}: found {} times".format(rel, header, len(hits)))
        i = j = hits[0]
        j += 1
        while j < len(lines):
            l = lines[j]
            if l.strip() and (len(l) - len(l.lstrip())) <= len(indent):
                break
            j += 1
        while j > i + 1 and not lines[j - 1].strip():          # keep the blank lines that separate definitions
            j -= 1
        out = lines[:i] + new.rstrip("\n").split("\n") + lines[j:]
        self.write(rel, "\n".join(out))
        log("  replace definition {!r} in {} (lines {}-{})".format(header, rel, i + 1, j))


def run_bash_block(tree, text, label):
    env = dict(os.environ, PY=PY)
    r = subprocess.run(["bash", "-c", text], cwd=tree.root, env=env, capture_output=True, text=True)
    log("  ran {} -> exit {}; stdout: {}".format(label, r.returncode, r.stdout.strip()[-300:]))
    if r.returncode:
        log(r.stderr[-2000:])
        raise SystemExit("bash block {} failed".format(label))


def apply_diff(tree, text, label):
    """Apply one unified diff block strictly: every hunk's old lines must match the file exactly at the line its header
    states (shifted only by the earlier hunks of the same block); a hunk that matches only elsewhere is applied there
    and its offset reported, one that matches nowhere stops the rebuild. (macOS's patch 2.0 applies at an offset
    silently, which would hide exactly what this reconstruction needs to see.)"""
    files = re.split(r"^--- a/(\S+)\n\+\+\+ b/\S+\n", text, flags=re.M)
    report = []
    for rel, body in zip(files[1::2], files[2::2]):
        lines = tree.read(rel).split("\n")
        delta = 0
        for h in re.finditer(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@.*\n((?:(?:[ +\-\\].*|)\n)*?)(?=^@@|\Z)", body, re.M):
            s_old, n_old, s_new, n_new = int(h.group(1)), int(h.group(2) or 1), int(h.group(3)), int(h.group(4) or 1)
            hl = h.group(5).split("\n")[:-1]
            old = [l[1:] if l else "" for l in hl if l[:1] in (" ", "-", "")]
            new = [l[1:] if l else "" for l in hl if l[:1] in (" ", "+", "")]
            while old and len(old) > n_old and old[-1] == "" and new and new[-1] == "":   # trailing blank padding
                old.pop(); new.pop()
            if len(old) != n_old or len(new) != n_new:
                raise SystemExit("{} {}: hunk @@ -{},{} +{},{}: counts {} / {}".format(label, rel, s_old, n_old, s_new, n_new, len(old), len(new)))
            if s_new != s_old + delta:
                report.append("{} hunk -{} header's new start {} != {}".format(rel, s_old, s_new, s_old + delta))
            pos = s_old - 1 + delta
            if lines[pos:pos + len(old)] != old:
                cands = [i for i in range(len(lines) - len(old) + 1) if lines[i:i + len(old)] == old]
                if len(cands) != 1:
                    raise SystemExit("{} {}: hunk @@ -{} matches {} places (expected at line {})".format(label, rel, s_old, len(cands), pos + 1))
                report.append("{} hunk -{} applied at line {} (offset {:+d})".format(rel, s_old, cands[0] + 1, cands[0] - pos))
                pos = cands[0]
            lines[pos:pos + len(old)] = new
            delta += len(new) - len(old)
        tree.write(rel, "\n".join(lines))
    log("  diff {}: {}".format(label, "; ".join(report) if report else "every hunk exactly at its stated line"))


# ------------------------------------------------------------------------------------------------- stages
def stage_base(tree):
    log("stage base: git archive {} reentry_model tests data analysis".format(BASE_COMMIT))
    data = subprocess.run(["git", "-C", REPO, "archive", "--format=tar", BASE_COMMIT, "reentry_model", "tests", "data", "analysis"],
                          check=True, capture_output=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        tar.extractall(tree.root, filter="data")
    n = sum(len(f) for _, _, f in os.walk(tree.root))
    log("  {} files".format(n))


def stage_task1(tree):
    log("stage task1: mesh.py, its tests and the band = 0 test pins from {} (sub-plan 01 = Task 1 + the dense band, ported from the prototype)".format(TASK1_COMMIT))
    # conftest.py and the Step 2 reference test take 187cd51's band = 0 pins as well: with the default band the
    # Carslaw-Jaeger fixture fails its margin (measured here: 5.005 K against 5 K, exactly the failure 187cd51's message
    # records), and the prototype's unit tier had no such failure after sub-plan 01 (216 passed + the five known)
    for rel in ("reentry_model/mesh.py", "tests/test_reentry_model_mesh.py", "tests/conftest.py",
                "tests/test_reentry_model_reference_thermal.py"):
        text = subprocess.run(["git", "-C", REPO, "show", "{}:{}".format(TASK1_COMMIT, rel)], check=True,
                              capture_output=True, text=True).stdout
        tree.write(rel, text)
        log("  {} ({} lines)".format(rel, text.count("\n")))


def stage_plan(tree):
    d = Doc(PLAN)
    log("stage plan: Tasks 2-15 of {}".format(os.path.relpath(PLAN, REPO)))
    # Task 2 -- material
    log(" Task 2")
    run_bash_block(tree, d.block(796, "Generate them from the packaged no-melt file", "bash"), "plan:796 (material files)")
    tree.append("tests/test_reentry_model_material.py", d.block(822, "Append to `tests/test_reentry_model_material.py`"))
    tree.replace_file("reentry_model/material.py", d.block(894, "Replace `reentry_model/material.py`"))
    # Task 3 -- thermal core
    log(" Task 3")
    tree.edit("tests/test_reentry_model_thermal.py", d.block(1138, "replace the body of `test_operators_match_scikit_fem_assembly`"),
              d.block(1151, "with"), "(plan:1138 -> plan:1151)")
    tree.append("tests/test_reentry_model_thermal.py", d.block(1171, "Then append to the file"))
    tree.edit("tests/test_reentry_model_fenicsx.py", 'thermal.thermal_solver("fenicsx", lumped_mass=True)\n',
              'thermal.thermal_solver("fenicsx", lumped_mass=False)           # the nodal-enthalpy (lumped) capacity only\n',
              "(plan:1316 prose)")
    tree.append("tests/test_reentry_model_fenicsx.py", d.block(1318, "In `tests/test_reentry_model_fenicsx.py` change the constructor"))
    tree.replace_file("reentry_model/thermal/__init__.py", d.block(1370, "Replace `reentry_model/thermal/__init__.py`"))
    tree.replace_file("reentry_model/thermal/skfem_backend.py", d.block(1426, "Replace `reentry_model/thermal/skfem_backend.py`"))
    tree.replace_file("reentry_model/thermal/fenicsx_backend.py", d.block(1767, "Replace `reentry_model/thermal/fenicsx_backend.py`"))
    tree.edit("reentry_model/cli.py", '    th.add_argument("--lumped-mass", action="store_true")\n',
              d.block(2027, 'replace `th.add_argument("--lumped-mass", action="store_true")` with'), "(plan:2027)")
    tree.edit("reentry_model/cli.py", "lumped_mass=args.lumped_mass)", "lumped_mass=not args.consistent_mass)", "(plan:2031 a)")
    tree.edit("reentry_model/cli.py", '"lumped_mass": args.lumped_mass}', '"lumped_mass": not args.consistent_mass}', "(plan:2031 b)")
    # Task 4 -- dispersion
    log(" Task 4")
    tree.replace_file("tests/test_reentry_model_dispersion.py", d.block(2068, "Create `tests/test_reentry_model_dispersion.py`"), "create")
    tree.replace_file("reentry_model/dispersion.py", d.block(2124, "Create `reentry_model/dispersion.py`"), "create")
    gen = d.block(2226, "Generate the packaged table and run the tests", "bash").split("\n")[0]
    run_bash_block(tree, gen, "plan:2226 line 1 (girin_dispersion.json)")
    # Task 5 -- gas + surface flow
    log(" Task 5")
    tree.replace_file("tests/test_reentry_model_surface_flow.py", d.block(2261, "Create `tests/test_reentry_model_surface_flow.py`"), "create")
    tree.insert_after_line("reentry_model/gas.py", 'AIR = "N2:0.79, O2:0.21"', "AVOGADRO = 6.02214076e23\n", "(plan:2400 a)")
    tree.insert_after_line("reentry_model/gas.py", "X: dict           # mole fractions > 1e-6",
                           d.block(2402, "In `GasState`, after `X: dict"), "(plan:2402 b)")
    gas = tree.read("reentry_model/gas.py")
    m = re.search(r"\n        return GasState\(.*?\)\n", gas, re.S)
    tree.edit("reentry_model/gas.py", m.group(0)[1:], d.block(2410, "In `_state()` replace the `return GasState(...)` statement"), "(plan:2410 c)")
    tree.insert_before_line("reentry_model/gas.py", "def wall(self, T_w, p):", d.block(2418, "Before `def wall(self, T_w, p):` insert"), "(plan:2418 d)")
    tree.replace_file("reentry_model/surface_flow.py", d.block(2434, "Create `reentry_model/surface_flow.py`"), "create")
    # Task 6 -- film
    log(" Task 6")
    tree.replace_file("tests/test_reentry_model_film.py", d.block(2804, "Create `tests/test_reentry_model_film.py`"), "create")
    tree.replace_file("reentry_model/film.py", d.block(2890, "Create `reentry_model/film.py`"), "create")
    # Task 7 -- spray
    log(" Task 7")
    tree.replace_file("data/reference_values/girin2017_table1.json", d.block(3036, "`data/reference_values/girin2017_table1.json`"), "create")
    tree.replace_file("data/reference_values/girin1994_tables.json", d.block(3123, "`data/reference_values/girin1994_tables.json`"), "create")
    tree.replace_file("tests/test_reentry_model_spray.py", d.block(3214, "Create `tests/test_reentry_model_spray.py`"), "create")
    tree.replace_file("reentry_model/spray.py", d.block(3539, "Create `reentry_model/spray.py`"), "create")
    # Task 8 -- Girin's cases
    log(" Task 8")
    tree.replace_file("tests/test_reentry_model_girin.py", d.block(3976, "Create `tests/test_reentry_model_girin.py`"), "create")
    tree.replace_file("reentry_model/girin_case.py", d.block(4038, "Create `reentry_model/girin_case.py`"), "create")
    tree.replace_file("analysis/girin_reference.py", d.block(4177, "Create `analysis/girin_reference.py`"), "create")
    # Task 9 -- the melting body
    log(" Task 9")
    tree.replace_file("tests/test_reentry_model_melting.py", d.block(4379, "Create `tests/test_reentry_model_melting.py`"), "create")
    tree.replace_file("reentry_model/body.py", d.block(4667, "Replace `reentry_model/body.py`"))
    write_atdb_disc(tree, d.block(5596, "Extract the flat-disc endpoint of the drag family", "bash"))
    tree.insert_after_line("reentry_model/aero.py", 'DEFAULT_ATDB = os.path.join(DATA_DIR, "atdb_sphere.json")',
                           d.block(5625, "add, next to `DEFAULT_ATDB`"), "(plan:5625 DISC_ATDB)")
    tree.replace_def("reentry_model/aero.py", "def drag_coefficient(", d.block(5631, "give `drag_coefficient` the shape factor"))
    tree.edit("reentry_model/trajectory.py", d.block(5650, "replace the three lines"), d.block(5658, "with"), "(plan:5650 -> 5658)")
    tree.insert_after_line("reentry_model/aero.py", 'DEFAULT_ATDB = os.path.join(DATA_DIR, "atdb_sphere.json")',
                           d.block(5670, "the factor-2 step in `cd_continuum` stalls"), "(plan:5670 MACH_SWITCH)")
    tree.replace_def("reentry_model/aero.py", "def cd_continuum(", d.block(5676, "replace `SphereDragTables.cd_continuum` with"), indent="    ")
    # the class docstring's last sentence (plan:5688, prose); wrapped at the module's 120 columns
    tree.edit("reentry_model/aero.py", '    simply clamped (Kn is negligible wherever Ma < 5)."""\n',
              "    simply clamped (Kn is negligible wherever Ma < 5). The factor-2 step is applied as a smooth (cubic) ramp over\n"
              "    Ma 0.98-1.02: a discontinuous C_D stalls the adaptive integrator when a light body hovers at its transonic\n"
              "    terminal velocity (a melting remnant, Step 3: 1.3e5 RHS evaluations in one macro step, measured 2026-09-21); the\n"
              "    reference spheres cross Ma 1 in a fraction of a second, where the ramp changes nothing measurable (the drag test\n"
              "    excludes |Ma - 1| <= 0.02 rows for SESAM's 3-decimal Mach column).\"\"\"\n", "(plan:5688 docstring)")
    tree.edit("tests/test_reentry_model_aero.py",
              "        assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.0) == 0.898818\n",
              d.block(5690, "In `tests/test_reentry_model_aero.py` replace the assertion"), "(plan:5690)")
    # Task 10 -- the coupled loop
    log(" Task 10")
    tree.append("tests/test_reentry_model_coupled.py", d.block(5730, "Append to `tests/test_reentry_model_coupled.py`"))
    tree.replace_file("reentry_model/coupled.py", d.block(5796, "Replace `reentry_model/coupled.py`"))
    # Task 11 -- the melting references (not generated: no DRAMA on this machine; the tests are kept)
    log(" Task 11 (references NOT generated: DRAMA 4.1.4 is not installed; tests transcribed)")
    tree.append("tests/test_reentry_model_data.py", d.block(6077, "Append to `tests/test_reentry_model_data.py`"))
    tree.edit("tests/test_reentry_model_aero.py", d.block(6115, "replace the last assertion"), d.block(6121, "with"), "(plan:6115 -> 6121)")
    # the appended block uses pytest.approx, which test_reentry_model_data.py never imported (Step 2 did not need it);
    # the prototype's test passed (the five known failures do not include it), so it had the import: added here
    tree.edit("tests/test_reentry_model_data.py", "import os\n\nimport reentry_model\n", "import os\n\nimport pytest\n\nimport reentry_model\n",
              "(reconstructed: import pytest)")
    # Task 12 -- metrics, plots, viz
    log(" Task 12")
    tree.append("tests/test_reentry_model_compare.py", d.block(6157, "Append to `tests/test_reentry_model_compare.py`"))
    tree.append("tests/test_reentry_model_viz.py", d.block(6222, "Append to `tests/test_reentry_model_viz.py`"))
    tree.insert_after_line("reentry_model/sesam_io.py", "integrated_heat: np.ndarray = None      # J, integral of the convective heat",
                           d.block(6284, "In the `Reference` dataclass, after `integrated_heat"), "(plan:6284)")
    tree.insert_after_line("reentry_model/sesam_io.py", 'convective_heat=optional("convective_heat_W"), rad_cooling=optional("rad_cooling_W"), integrated_heat=optional("integrated_heat_J"),',
                           d.block(6291, "in `load_reference`, after the `convective_heat=..."), "(plan:6291)")
    tree.edit("reentry_model/compare.py",
              'THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")\nCONTINUUM_KN = 0.01\n',
              d.block(6299, "Replace the two constant lines"), "(plan:6299)")
    tree.append("reentry_model/compare.py", d.block(6309, "and append to the file"))
    tree.replace_file("reentry_model/viz.py", d.block(6469, "Replace `reentry_model/viz.py`"))
    # Task 13 -- command line
    log(" Task 13")
    tree.append("tests/test_reentry_model_cli.py", d.block(6678, "Append to `tests/test_reentry_model_cli.py`"))
    tree.replace_file("reentry_model/cli.py", d.block(6755, "Replace `reentry_model/cli.py`"))
    # Task 14 -- verification drivers
    log(" Task 14")
    tree.replace_file("tests/test_reentry_model_reference_melt.py", d.block(7218, "Create `tests/test_reentry_model_reference_melt.py`"), "create")
    tree.replace_file("analysis/melt_verification.py", d.block(7270, "Create `analysis/melt_verification.py`"), "create")
    tree.replace_file("analysis/melt_sensitivity.py", d.block(7380, "Create `analysis/melt_sensitivity.py`"), "create")
    # Task 15 -- only the package docstring is code; README/assumptions/spec/facts-note edits are repo documents
    log(" Task 15 (package docstring only)")
    init = tree.read("reentry_model/__init__.py")
    first, second = init.split("\n")[0], init.split("\n")[1]
    tree.edit("reentry_model/__init__.py", first + "\n" + second + "\n",
              '"""First-principles re-entry model of a solid sphere, Steps 1-3: trajectory, coupled 3D heat transfer, melting and\n'
              'melt spraying (specs under docs/superpowers/specs/)."""\n', "(plan:8107)")


def write_atdb_disc(tree, script):
    """Task 9 step 4 extracts reentry_model/data/atdb_disc.json from DRAMA's ATDB_CYLINDER.nc. DRAMA 4.1.4 is not
    installed (memory note of 2026-10-06), so the file is written from the values the committed plan itself states,
    and nothing else: cd_continuum at all six Mach numbers (plan line 5621), cd_free_molecular at Ma 10 = 2.202072 and
    heat_flux_factor_continuum = 0.329561 at every Mach number (the Task 11 test, plan lines 6104-6107). The values the
    documents do not state are null; the two tests that read them fail until the file is re-extracted on a machine
    with DRAMA (the script below, kept verbatim in the provenance)."""
    import json
    if os.path.exists("/Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_CYLINDER.nc"):
        run_bash_block(tree, script, "plan:5596 (atdb_disc.json from DRAMA)")
        return
    doc = {
        "_provenance": ("RECONSTRUCTED 2026-10-07 without DRAMA (prototype/rebuild/rebuild.py): only the values the committed "
                        "plan states are present, the rest are null. The original was the values of "
                        "/Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_CYLINDER.nc (ESA DRAMA 4.1.4 aerothermal database for "
                        "the primitive CYLINDER, Hyperschall Technologie Goettingen GmbH) at angle of attack 0 (the flat face normal "
                        "to the flow) and the thinnest tabulated aspect ratio; re-extract it with the script of plan Task 9 step 4."),
        "mach": [5.0, 10.0, 15.0, 20.0, 25.0, 30.0],
        "cd_free_molecular": [None, 2.202072, None, None, None, None],
        "cd_continuum": [1.801341, 1.824148, 1.828404, 1.829896, 1.830587, 1.830962],
        "heat_flux_factor_free_molecular": [None] * 6,
        "heat_flux_factor_continuum": [0.329561] * 6,
    }
    with open(tree.p("reentry_model/data/atdb_disc.json"), "w") as fh:
        json.dump(doc, fh, indent=2)
    log("  wrote reentry_model/data/atdb_disc.json PARTIALLY (no DRAMA: the plan's stated values only, the rest null)")


def stage_band13(tree):
    """Sub-plan 13's amendment of 2026-09-27 put "the dense-band line" into the prototype's cli.py (its 2026-10-02 diff
    says it is against "the body below plus the dense-band line"), but no document carries the line itself. What is
    known: the 2026-10-02/03/05 cli.py diffs apply exactly up to line 157 and exactly 3 lines late from line 207 on,
    so it added 3 lines in build_thermal's region; and the melting runs of facts 46-68 mesh the 100 mm sphere with
    177 363 tetrahedra and 18 830 patches, i.e. the 15 mm default band. Reconstructed as the band for melting runs and
    band 0 (the Step 2 field) otherwise, mirroring 187cd51's pin of the Step 2 package -- 3 lines."""
    log("stage band13: sub-plan 13's dense-band line (reconstructed, 3 lines)")
    tree.edit("reentry_model/cli.py",
              "    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, layers=layers, layer_thickness=args.layer_thickness * 1e-3)\n",
              "    # the dense band (sub-plan 01, amendment of 2026-09-27) is for melting runs; a run that does not melt keeps the\n"
              "    # Step 2 field (band 0), on which the committed Step 2 verification numbers were measured\n"
              "    band = mesh.DEFAULT_BAND if melting else 0.0\n"
              "    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, layers=layers, layer_thickness=args.layer_thickness * 1e-3, band=band)\n",
              "(reconstructed)")


def stage_fact37(tree):
    """Fact 37 (2026-09-27): regenerated into plan3/plan.md but never copied over the committed master plan; carried as
    snippets in sub-plans 09 and 10. The exact placement and the surrounding lines are fixed by the context lines of the
    2026-10-02 diffs (body.py hunks at 232/400/730, coupled.py hunk at 149/183, test_coupled hunk at 115/146); the
    test lines of fact 37 that no document carries are reconstructed from its prose (README)."""
    s9, s10 = Doc(os.path.join(SUB, "09-melting-body.md")), Doc(os.path.join(SUB, "10-coupled-loop.md"))
    log("stage fact37")
    snippet = s9.block(33, "`MeltingBody` gains the attribute `last_face_ids` and the method `on_current_surface`", "python")
    init_line, step3_line = snippet.split("\n")[0], snippet.split("\n")[1]
    method = "\n".join(snippet.split("\n")[3:]).rstrip("\n") + "\n\n"     # the snippet's own trailing blank line
    assert init_line.strip().startswith("self.last_face_ids = None") and "face_ids.copy()" in step3_line
    # (1) in __init__, beside last_flow/last_spray -- the comment as the 2026-10-02 diff's context line has it
    tree.edit("reentry_model/body.py", "        self.last_flow = self.last_spray = None\n        self.last_melt = {",
              "        self.last_flow = self.last_spray = None\n"
              "        self.last_face_ids = None                                          # face ids of the surface they were evaluated on\n"
              "        self.last_melt = {", "(fact 37: __init__)")
    # (2) at step (iii): where the spray stores last_flow/last_spray, before the deaths of step (v)
    tree.insert_after_line("reentry_model/body.py", "self.last_flow, self.last_spray, self.last_b = flow, res, b",
                           "        self.last_face_ids = self.surface.face_ids.copy()                  # the surface they belong to: the deaths come later\n",
                           "(fact 37: step iii)")
    # (3) the method, first under the reporting header (the 2026-10-02 diff's context at its line 486-487)
    tree.insert_after_line("reentry_model/body.py", "# -- reporting ------", method, "(fact 37: on_current_surface)")
    # the writer: sub-plan 10's snippet replaces the index test, and the docstring is the one the 2026-10-02 diff removes
    carry = s10.block(28, "`MeltingBody.on_current_surface`, which carries them across the step's element deaths", "python")
    old = ("        n = surface.n_patches\n"
           "        res, flow = body.last_spray, body.last_flow\n"
           "        same = res is not None and res.r.size == n\n"
           '        poly.cell_data["we_s"] = res.we_s if same else np.zeros(n)\n'
           '        poly.cell_data["closure"] = flow.closure.astype(float) if same else np.ones(n)\n'
           '        poly.cell_data["kn_local"] = flow.kn_local if same else np.full(n, np.nan)\n'
           '        poly.cell_data["p_w"] = flow.p_w if same else np.zeros(n)\n'
           '        poly.cell_data["tau"] = flow.tau if same else np.zeros(n)\n'
           '        poly.cell_data["r_droplet"] = np.where(res.dm > 0.0, res.r, np.nan) if same else np.full(n, np.nan)\n'
           '        poly.cell_data["release_rate"] = res.dm / surface.areas if same else np.zeros(n)\n')
    tree.edit("reentry_model/coupled.py", old, carry, "(fact 37: the writer)")
    d1002 = Doc(os.path.join(SUB, "10-coupled-loop.md")).block(81, None, "diff")
    hunk = d1002.split("@@ -149,9 +152,10 @@\n")[1].split("@@")[0].split("\n")
    ctx = [l[1:] for l in hunk if l.startswith(" ") or l.startswith("-")]
    docstring = "\n".join(ctx[1:6]) + "\n"                      # the five docstring lines of the fact-37 writer
    coupled = tree.read("reentry_model/coupled.py")
    m = re.search(r'(def write_vtk_frame\(output_dir, k, body, loads\):\n)(    """.*?"""\n)', coupled, re.S)
    tree.edit("reentry_model/coupled.py", m.group(2), docstring, "(fact 37: write_vtk_frame docstring from the 2026-10-02 context)")
    # the test: the two comment lines are the 2026-10-02 diff's context; the twelve lines after them are reconstructed
    tree.edit("tests/test_reentry_model_coupled.py",
              "        assert key in poly.cell_data\n    files = coupled.write_particles(run_dir, b, hist)\n",
              "        assert key in poly.cell_data\n" + FACT37_TEST + "    files = coupled.write_particles(run_dir, b, hist)\n",
              "(fact 37: test lines, 2 from the 2026-10-02 context + 12 reconstructed from fact 37's prose)")


FACT37_TEST = '''\
    # the flow fields are the frame's own step's, carried across that step's element deaths: the wall pressure is positive
    # on the windward patches and nowhere above the stagnation value the history recorded from the same evaluation
    # (the twelve lines below are reconstructed 2026-10-07 from plan fact 37's prose: see prototype/README.md)
    from scipy.stats import spearmanr
    assert c["n_dead_elements"][20] > c["n_dead_elements"][19]                     # the frame's own step had a death
    tri = np.asarray(poly.faces).reshape(-1, 4)[:, 1:]
    pts = np.asarray(poly.points)[tri]
    normal = np.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    angle = np.degrees(np.arccos(np.clip(normal[:, 0], -1.0, 1.0)))              # from +x, the flight direction
    p_w = np.asarray(poly.cell_data["p_w"])
    wet = (angle < 90.0) & (p_w > 0.0)                                             # faces the deaths exposed carry 0
    assert wet.sum() > 0.5 * (angle < 90.0).sum() and p_w.max() <= c["p_w_stag_Pa"][20] * (1.0 + 1e-9)
    assert spearmanr(angle[wet], p_w[wet])[0] < -0.9                               # falls with the angle: not mis-mapped
'''


AMENDMENTS = {
    # date: [(sub-plan file, opening line of the diff block, sentinel), ...] in the order the documents give them
    "d1002": [("06-melt-film.md", 49, "Code (the tested change; apply to `reentry_model/film.py`"),
              ("06-melt-film.md", 90, "Test (the tested code, appended to `tests/test_reentry_model_film.py`)"),
              ("09-melting-body.md", 174, "Code (the tested change to `reentry_model/body.py`)"),
              ("09-melting-body.md", 583, "Tests (the tested code, appended to `tests/test_reentry_model_melting.py`)"),
              ("10-coupled-loop.md", 81, "Code (the tested change to `reentry_model/coupled.py`)"),
              ("10-coupled-loop.md", 132, "Tests (the tested change to `tests/test_reentry_model_coupled.py`)"),
              ("13-cli-wiring.md", 56, "Code (the tested change to `reentry_model/cli.py`"),
              ("13-cli-wiring.md", 121, "Tests (the tested change to `tests/test_reentry_model_cli.py`)"),
              ("14-verification-runs.md", 73, "A sensitivity row, `nodeep`")],
    "d1003": [("09-melting-body.md", 890, "Code (the tested change to `reentry_model/body.py`)"),
              ("09-melting-body.md", 1074, "Tests (the tested code, appended to `tests/test_reentry_model_melting.py`"),
              ("09-melting-body.md", 1268, None),
              ("10-coupled-loop.md", 216, "Code (the tested change to `reentry_model/coupled.py`)"),
              ("10-coupled-loop.md", 241, "Test (the tested change to `tests/test_reentry_model_coupled.py`)"),
              ("13-cli-wiring.md", 187, "Code (the tested change to `reentry_model/cli.py`"),
              ("13-cli-wiring.md", 254, "Tests (the tested change to `tests/test_reentry_model_cli.py`)"),
              ("14-verification-runs.md", 141, "A sensitivity row, `nocascade`")],
    "d1005": [("13-cli-wiring.md", 340, "Code (the tested change to `reentry_model/cli.py`"),
              ("13-cli-wiring.md", 445, "Tests (the tested change to `tests/test_reentry_model_cli.py`)"),
              ("14-verification-runs.md", 224, "A sensitivity row, `seed1`")],
}


def stage_diffs(tree, key):
    log("stage {}: {} diff blocks".format(key, len(AMENDMENTS[key])))
    for name, line, sentinel in AMENDMENTS[key]:
        text = Doc(os.path.join(SUB, name)).block(line, sentinel, "diff")
        apply_diff(tree, text, "{}:{}".format(name, line))


def post_subplan02(tree):
    """Post-reconstruction addition 01 (2026-10-07, requested for Spheral M1, which needs AA7075_scheil): sub-plan 02's
    material amendments of 2026-09-27/28 -- the AA7075_scheil variant and AA7075-empiricaldata -- which the documents say
    were tested in throwaway copies and never merged into proto3 (facts 46 and 54 list the copies' contents). Its body is
    the master plan's Task 2 re-tested: the material-file generator (writing AA7075.json and AA7075_range.json
    byte-identical, checked), material.py replaced, and two test blocks appended after Task 2's (identical to plan:822)."""
    d = Doc(os.path.join(SUB, "02-material-properties.md"))
    log(" post 01: sub-plan 02 (AA7075_scheil, AA7075-empiricaldata)")
    mats = tree.p("reentry_model/data/materials")
    before = {n: open(os.path.join(mats, n), "rb").read() for n in ("AA7075.json", "AA7075_range.json")}
    run_bash_block(tree, d.block(134, "Generate them from the packaged no-melt file", "bash"), "02-material-properties.md:134 (four material files)")
    for n, b in before.items():
        if open(os.path.join(mats, n), "rb").read() != b:
            raise SystemExit("sub-plan 02's generator changed " + n)
    log("  AA7075.json and AA7075_range.json byte-identical after sub-plan 02's generator")
    if d.block(172, "Append to `tests/test_reentry_model_material.py`") != Doc(PLAN).block(822):
        raise SystemExit("sub-plan 02's first test block is not the master plan's")
    tree.append("tests/test_reentry_model_material.py", d.block(238, "Then append the Scheil variant's tests"))
    tree.append("tests/test_reentry_model_material.py", d.block(306, "Then append the empirical-data file's test"))
    tree.replace_file("reentry_model/material.py", d.block(334, "Replace `reentry_model/material.py`"))


def stage_post(tree, materials_only=False, skip=()):
    post = os.path.join(HERE, "post")
    names = sorted(n for n in os.listdir(post) if n.endswith(".diff")) if os.path.isdir(post) else []
    log("stage post: {} post-reconstruction additions (NOT part of the recovered state)".format(len(names) + 1))
    post_subplan02(tree)
    for n in ([] if materials_only else names):
        if n[:2] in skip:
            log("  post/{}: skipped (--post-skip)".format(n))
            continue
        apply_diff(tree, open(os.path.join(post, n)).read(), "post/" + n)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(REPO, "prototype", "proto3"))
    ap.add_argument("--upto", choices=STAGES[:-1], default="d1005", help="last recovered stage to apply")
    ap.add_argument("--post", action="store_true", help="also apply the post-reconstruction additions")
    ap.add_argument("--post-materials-only", action="store_true",
                    help="with --post: apply only post 01 (sub-plan 02's materials), not the export diffs -- the control tree for the "
                         "write-only check of the export change")
    ap.add_argument("--post-skip", default="", help="with --post: comma-separated numbers of post/*.diff to leave out, e.g. "
                    "04,06 -- the control tree for the write-only check of the export diffs (2026-10-08)")
    ap.add_argument("--force", action="store_true", help="delete an existing output tree first")
    args = ap.parse_args(argv)
    if os.path.exists(args.out):
        if not args.force:
            raise SystemExit("{} exists (use --force to rebuild it)".format(args.out))
        shutil.rmtree(args.out)
    os.makedirs(args.out)
    tree = Tree(args.out)
    head = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    log("rebuild of {} from the documents at {} (HEAD {})".format(args.out, REPO, head))
    funcs = {"base": stage_base, "task1": stage_task1, "plan": stage_plan, "band13": stage_band13, "fact37": stage_fact37,
             "d1002": lambda t: stage_diffs(t, "d1002"), "d1003": lambda t: stage_diffs(t, "d1003"),
             "d1005": lambda t: stage_diffs(t, "d1005")}
    for st in STAGES[:STAGES.index(args.upto) + 1]:
        funcs[st](tree)
    if args.post:
        stage_post(tree, materials_only=args.post_materials_only, skip=tuple(x.strip() for x in args.post_skip.split(",") if x.strip()))
    for root, dirs, files in os.walk(args.out):                 # no caches in the tree
        for dname in list(dirs):
            if dname == "__pycache__":
                shutil.rmtree(os.path.join(root, dname))
                dirs.remove(dname)
    n = sum(len(f) for _, _, f in os.walk(args.out))
    log("done: {} files".format(n))
    with open(os.path.join(HERE, "rebuild_log.txt" if os.path.abspath(args.out) == os.path.join(REPO, "prototype", "proto3")
                           else "rebuild_log_{}.txt".format(os.path.basename(args.out))), "w") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
