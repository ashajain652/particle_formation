# Stage 0 + S — the combined refactor: the shared physics package and the architecture flag Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the architecture-independent modules of `reentry_model` — the Steps 1–2 modules and Step 3's ported
physics — into a new package `reentry_model/shared/`, keep every old module path working as an alias of the same
module object, and add `run --architecture skyscrapers|levelset` (levelset exits 2 until built), with every compared
result identical before and after.

**Architecture:** Each moved module is `git mv`-ed whole into `reentry_model/shared/` and its old path becomes a
four-line alias that replaces itself in `sys.modules` with the shared module object, so `reentry_model.spray is
reentry_model.shared.spray` and pickled states that name the old path still load. Four pieces are split out of
skyscrapers modules: the `Body` seam (`body.py`), the thermal backend protocol (`thermal/__init__.py`), the droplet
output writer (`coupled.py`) and the common run-command parts (`cli.py`); each skyscrapers module imports them from
`shared` and re-exports the names it used to define. The skins, the moving mesh, the coupled loop, the element
model's melting body and the thermal backends stay where they are.

**Tech Stack:** Python 3.12 (conda `drama_env`), pytest with pytest-xdist, argparse; no new dependency.

**Spec:** `Reentry Model/docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4.2 (the shared
package, both tables), §4.3 (the flag), §4.4 (coordination and acceptance), §14 item "0+S". The handoff note
`Reentry Model/docs/superpowers/plans/2026-11-11-skyscrapers-refactor-handoff.md` gives the verification tool, the
baseline and the known defects. This plan supersedes `00-stage0-isolation-and-shared-package.md` and reuses its tested
pieces from branch `levelset-stage0` (head `da4c549`; read with `git show da4c549:<path>`):
`reentry_model/shared/__init__.py`, the alias files (e.g. `reentry_model/atmosphere.py`),
`reentry_model/shared/{body,thermal,cli_common}.py`, `tests/test_reentry_model_shared.py`, and the `--architecture`
tests appended to `tests/test_reentry_model_cli.py`. Those files were written for the pre-port layout (repository-root
paths, `reentry_model_output/`); adapt them, do not copy them blind.

## Global Constraints

- Paths: every path below is relative to `Reentry Model/` of the worktree unless it starts with the repository top.
  Interpreter: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. Run every command from `Reentry Model/`.
- **Behaviour-preserving.** No physics, default, option, help text of an existing flag, run name, history column, JSON
  key, output file name or exit code changes. The proof is `analysis/skyscrapers_port/same_results.py compare` against
  Task 1's baseline: `identical`, exit 0. If anything differs, stop and report; never "fix" a result.
- **Work on branch `shared-refactor` from `main` at e912057 (the port, merged 2026-11-11), in the worktree
  `.claude/worktrees/shared-refactor`.** Never switch another checkout's branch. Do not push.
- `reentry_model/shared/` imports nothing from `reentry_model` outside itself except these package-root constants:
  `__version__`, `PROJECT_ROOT`, `REPO_ROOT`, `DATA_DIR`, `OUTPUT_DIR`, `DRAMA_DIR`, `DRAMA_DATA_DIR` (spec §4.2).
- The skins are never shared (spec §4.2): `skin.py`, `skinbody.py`, `columns.py`, `movemesh.py`, `mesh.py`,
  `spreading.py`, `film.py`, `restart.py`, `coupled.py` (except the droplet writer of Task 4), `thermal/*_backend.py`,
  `body.py`'s `ThermalBody`/`MeltingBody` and `cli.py`'s run driver stay in place.
- Do not touch `DRAMA Sphere Experiments/`, `Spheral/`, `reentry_model/droplet_reentry/`, `reentry_model/wake_column/`
  (handoff note). Spheral loads `reentry_model.material` and `reentry_model.cli.make_atmosphere` by name
  (`Spheral/spheral_frag/fe.py`); the aliases must keep both working.
- Known defects 1–5 of the handoff note are left as they are; do not fix them here.
- Do not edit the package while a measurement (Task 1, Task 7) is running. Run measurements on mains power
  (`pmset -g batt | head -1` says `AC Power`).
- Unit tier: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -n 4 -p no:cacheprovider -q`
  (about 24 minutes). Run it once serially first in a fresh worktree or after clearing `output/meshes/`, because of the
  mesh-cache race (defect 2); a single failure that passes on rerun with a warm cache is that race.
- Commit messages: one sentence saying what moved or changed, ending with a blank line and
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

- **A saved state pickled before the refactor** (it names `reentry_model.spray.SprayModel`, `reentry_model.material.
  Material`, …) must load and resume to the same result — pinned by `test_pickles_naming_old_paths_load` (Task 3) and
  the `resume` case of `same_results.py` (Task 7).
- **Spheral's by-name use** of `reentry_model.material` (its private node tables `_T_h`, `_h_nodes`, `_latent_slope`)
  and `reentry_model.cli.make_atmosphere` must keep working — pinned by `test_spheral_entry_points_resolve` (Task 6)
  and the Spheral tests run once in Task 7.
- **Import order:** importing `reentry_model.shared.X` before `reentry_model.X` must still give one module object —
  pinned by `test_import_order_does_not_matter` (Task 2), parametrised over every moved module.
- **`--architecture=levelset` in the equals form and `--architecture` given to `resume` or `compare`** — pinned by
  `test_architecture_levelset_exits_2_until_built` (parametrised over both forms) and
  `test_architecture_only_on_run` (Task 6).
- **Every `reentry_model` import in tests, analysis scripts, prototype-retirement READMEs' commands and Spheral** must
  resolve — pinned by `test_every_reentry_model_import_in_tests_and_analysis_resolves` (Task 7), which walks
  `tests/`, `analysis/` and `../Spheral/spheral_frag/fe.py`.

---

## File structure

| Path | Change | Responsibility |
|---|---|---|
| `reentry_model/shared/__init__.py` | create | the independence rule |
| `reentry_model/shared/{trajectory,earth,atmosphere,fap,aero,gas,constants,heating,material,sesam_io,compare,viz}.py` | `git mv` (Task 2) | Steps 1–2 modules |
| `reentry_model/shared/{rheology,dispersion,surface_flow,spray,girin_case}.py` | `git mv` (Task 3) | Step 3 physics (groups S1–S5) |
| `reentry_model/shared/droplets.py` | create (Task 4) | `write_particles`, `SUMMARY_COLUMNS`, `SKIN_SUMMARY_COLUMNS` (group S6) |
| `reentry_model/shared/body.py` | create (Task 5) | `Body` protocol, `sphere_mass`, `ConstantBody` |
| `reentry_model/shared/thermal.py` | create (Task 5) | `SIGMA_SB`, `MissingBackend`, `StepResult`, `ThermalSolver` |
| `reentry_model/shared/cli_common.py` | create (Task 6) | flight defaults and arguments, `parse_epoch`, `make_atmosphere`, `git_commit`, `provenance`, `ARCHITECTURES`, `architecture_of` |
| `reentry_model/{the 17 moved names}.py` | create (alias) | old path → the shared module object |
| `reentry_model/{body,coupled,cli}.py`, `reentry_model/thermal/__init__.py` | modify | import from `shared`, re-export the moved names |
| `tests/test_reentry_model_shared.py` | create, grown per task | aliases, re-exports, independence, pickles, imports |
| `tests/test_reentry_model_cli.py` | modify (append) | the `--architecture` flag |
| `analysis/skyscrapers_port/same_results.py` + its test | modify (Task 1) | add case `step1` |
| `README.md` (repository top), `Reentry Model/CLAUDE.md`, `docs/model_assumptions.md`, design §14 | modify (Task 7) | the two architectures and the shared package |

---

### Task 1: Worktree, branch, the plan, and the baseline

**Files:**
- Create: `docs/superpowers/plans/levelset/00s-combined-refactor-shared-package.md` (this plan)
- Modify: `analysis/skyscrapers_port/same_results.py`, its test file, its `README.md` (case `step1`)

**Interfaces:**
- Produces: the baseline folder `output/refactor_check/before/` (git-ignored) holding cases `step1`, `element`,
  `skins`, `resume`, `nomelt_elements`; the unit-tier count before the change. Every later task compares against it.

- [ ] **Step 1:** Check that `origin/main` holds the port (`git cat-file -e origin/main:"Reentry Model/docs/superpowers/plans/2026-11-11-skyscrapers-refactor-handoff.md"`), then from the repository top:
  `git worktree add .claude/worktrees/shared-refactor -b shared-refactor origin/main`. Record `git rev-parse --short HEAD` as BASE in the ledger.
- [ ] **Step 2:** Copy this plan into the path above, commit it alone.
- [ ] **Step 3: Write the failing test** `test_step1_case_is_a_whole_step1_flight` in
  `analysis/skyscrapers_port/test_port_tools.py`: `same_results.case_arguments()["step1"] == ["run", "--diameter", "100",
  "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331"]` (Step 1 only: no `--thermal`, the
  default NRLMSISE atmosphere, the whole flight; the tool adds `--outdir` and `--name`). Run it; expect FAIL (KeyError).
- [ ] **Step 4:** Add the case to the tool and one line to its README. Run the tool's tests; expect PASS. Commit.
- [ ] **Step 5: Record the baseline** (BASE code, mains power, nothing else editing):
  `"$PY" analysis/skyscrapers_port/same_results.py run --code "$PWD" --out output/refactor_check/before --cases step1,element,skins,resume,nomelt_elements --state "<repository top>/.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl"`.
  Expected: exit 0, about 6–7 minutes. Never write to the state file.
- [ ] **Step 6: Confirm the baseline is main = port and deterministic.**
  (a) Copy `element`, `skins`, `resume` into a scratch folder and `compare` it with
  `<repository top>/.claude/worktrees/skyscrapers-port/Reentry Model/output/port_check/port` → `identical`, exit 0.
  (b) Run `step1,nomelt_elements` a second time into `output/refactor_check/before2` and `compare` with the matching
  part of `before` → `identical`. If (b) differs (the old stage 0 found Step 2 JSON summaries not bit-reproducible
  run to run), record which fields differ and stop for the controller; do not widen the tool's skip list.
- [ ] **Step 7:** Unit tier, serial once to warm the mesh cache, then parallel. Record `N passed, M skipped` (port:
  702 passed, 1 skipped). Record the reference tier too:
  `"$PY" -m pytest tests/test_reentry_model_reference.py tests/test_reentry_model_reference_thermal.py -m reference -q -p no:cacheprovider`
  (port: 13 passed in about 9 minutes) and, separately, the Step 3 reference tests
  `"$PY" -m pytest tests/test_reentry_model_reference_melt.py tests/test_reentry_model_skinbody.py -m reference -q -p no:cacheprovider`
  (collected but never run on the port, run time unknown: run it in the background and record its BASE outcome,
  whatever it is; Task 7 must reproduce that outcome test for test).

---

### Task 2: The shared package and the twelve Steps 1–2 modules

**Files:**
- Create: `reentry_model/shared/__init__.py`; `tests/test_reentry_model_shared.py`
- Move: `reentry_model/{trajectory,earth,atmosphere,fap,aero,gas,constants,heating,material,sesam_io,compare,viz}.py` → `reentry_model/shared/`
- Create (aliases): the twelve old paths

**Interfaces:**
- Produces: `reentry_model.shared.<name>` for the twelve; `reentry_model.<name> is reentry_model.shared.<name>`.
  The alias file, exactly (one per module, name substituted):

```python
"""Moved to reentry_model.shared.atmosphere (stage 0 + S of the level-set design). This old path is kept as an alias:
importing it returns the shared module object itself, so attributes, private names, pickles and identity all match."""
import sys

from .shared import atmosphere as _shared

sys.modules[__name__] = _shared
```

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_shared.py` (adapt from
  `git show da4c549:tests/test_reentry_model_shared.py`): `MOVED = [the twelve]` (Tasks 3 grows it);
  `test_old_path_is_the_shared_module(name)` asserts `importlib.import_module("reentry_model." + name) is
  importlib.import_module("reentry_model.shared." + name)`; `test_import_order_does_not_matter(name)` does the same in a
  subprocess importing the shared path first; `test_data_paths_resolve` asserts `material`'s default material file,
  `aero`'s data tables, `fap`'s file and `sesam_io.MELTING_REFERENCE_DIR` exist. Run; expect FAIL (`reentry_model.shared`
  missing).
- [ ] **Step 2:** `git mv` the twelve, add `shared/__init__.py` (from da4c549, with the constant list of Global
  Constraints), change each moved module's `from . import <ROOT CONSTANT>` to `from .. import …` (only the
  package root uses `__file__`, so no other path changes), write the twelve aliases. Inside `reentry_model/` (outside
  `shared/`), switch imports of the moved modules to `from .shared import …` / `from .shared.X import …` only where the
  module is touched anyway in later tasks; the aliases make the rest work, and a later task's diff stays small.
- [ ] **Step 3:** Run the new tests; expect PASS. Run `same_results.py run` for `step1,element,nomelt_elements` into
  `output/refactor_check/task2/` and `compare` each against `before` → `identical`. Commit.

---

### Task 3: Step 3's physics (groups S1–S5): rheology, dispersion, surface flow, spray, Girin's cases

**Files:**
- Move: `reentry_model/{rheology,dispersion,surface_flow,spray,girin_case}.py` → `reentry_model/shared/`
- Create (aliases): those five old paths
- Modify: `tests/test_reentry_model_shared.py`

**Interfaces:**
- Consumes: Task 2's package and alias form.
- Produces: `reentry_model.shared.{rheology,dispersion,surface_flow,spray,girin_case}`. S1's material half and S3's
  "additions to gas" are already in `shared/material.py` and `shared/gas.py` from Task 2 (the port put Scheil's curve
  in `material.py` and the surface-flow gas helpers in `gas.py`; verify by reading them and say so in the report).

**The spray split (decided here, from the ported code):** `spray.py` is already only physics and records —
`SprayModel.evaluate` takes per-patch arrays (film depth, mass, `cap_depth`), and `source_block`/`source_rows`/
`SourceSpool`/`histogram`/`register_angle`/`enclosed_centroid` take arrays and paths, not a body. The skyscrapers
bookkeeping that decides those arrays is already outside it: `body.MeltingBody`'s spray step (element model) and
`skinbody.SprayPlan` (layered skins). So the split is: `spray.py` moves whole; nothing of `body.py` or `skinbody.py`
moves. Before moving, the implementer greps `spray.py` for any import of or attribute access on a skyscrapers object
(`skin`, `columns`, `movemesh`, `mesh`, `body`, `skinbody`); if any is found, stop and report rather than move it.
`analysis/girin_reference.py` stays in `analysis/` (already where the spec puts it).

- [ ] **Step 1: Write the failing tests:** extend `MOVED` with the five; add `test_pickles_naming_old_paths_load`:
  build a `spray.SprayModel(material.Material.<default AA7075 constructor used by cli>.liquid)` (read `cli.spray_model_of`
  for the call), pickle it with `pickle.dumps`, assert the bytes contain `b"reentry_model.shared.spray"`; then replace
  that substring with `b"reentry_model.spray"` (same pickle protocol; adjust the SHORT_BINUNICODE length byte, or
  use protocol 0 so the module name is plain text) and assert `pickle.loads` gives an object whose type is
  `reentry_model.shared.spray.SprayModel` and whose `__dict__` equals the original's. Run; expect FAIL.
- [ ] **Step 2:** `git mv` the five, fix their `from . import` lines (`from .. import DATA_DIR`; `from . import
  dispersion` stays relative inside `shared`), write the five aliases.
- [ ] **Step 3:** Run `tests/test_reentry_model_shared.py` and the Step 3 module tests
  (`tests/test_reentry_model_{spray,dispersion,surface_flow,rheology,girin}*.py`, whichever exist); expect PASS. Run
  `same_results.py` `element,skins,resume` into `output/refactor_check/task3/` and `compare` with `before` →
  `identical` (the `resume` case is the proof that a pre-refactor pickle loads and runs). Commit.

---

### Task 4: The droplet outputs (group S6)

**Files:**
- Create: `reentry_model/shared/droplets.py`
- Modify: `reentry_model/coupled.py` (remove `SUMMARY_COLUMNS`, `SKIN_SUMMARY_COLUMNS`, `write_particles`; re-export them)
- Modify: `tests/test_reentry_model_shared.py`

**Interfaces:**
- Produces: `reentry_model.shared.droplets.write_particles(run_dir, body, history, window=10.0) -> dict[str, str]`,
  `SUMMARY_COLUMNS: list[str]`, `SKIN_SUMMARY_COLUMNS: list[tuple[str, str]]` — moved verbatim (body, docstring,
  values); `body` is duck-typed (only `body.source_rows` is read). `coupled.write_particles is
  shared.droplets.write_particles`, and likewise the two lists.

- [ ] **Step 1: Write the failing test** `test_coupled_reexports_the_droplet_writer`: the three identities above, and
  `inspect.getsource(shared.droplets)` contains no `skin`/`column` module import. Run; expect FAIL.
- [ ] **Step 2:** Move the three names into `shared/droplets.py` (it imports `spray` from `.spray`); in `coupled.py`,
  `from .shared.droplets import SKIN_SUMMARY_COLUMNS, SUMMARY_COLUMNS, write_particles  # noqa: F401  re-exported`.
  Keep every other use in `coupled.py` working unchanged.
- [ ] **Step 3:** Run the test, `tests/test_reentry_model_coupled*.py` and any test that calls `write_particles`
  (`grep -rl write_particles tests`); expect PASS. `same_results.py` `element,skins` into `output/refactor_check/task4/`,
  `compare` → `identical`. Commit.

---

### Task 5: The body seam and the thermal protocol

**Files:**
- Create: `reentry_model/shared/body.py`, `reentry_model/shared/thermal.py`
- Modify: `reentry_model/body.py`, `reentry_model/thermal/__init__.py`, `tests/test_reentry_model_shared.py`

**Interfaces:**
- Produces: `shared.body.{Body, sphere_mass, ConstantBody}` (moved verbatim from `body.py` lines 13–62 as ported);
  `shared.thermal.{SIGMA_SB, MissingBackend, StepResult, ThermalSolver}` (moved verbatim from
  `thermal/__init__.py`, including every `StepResult` field the port added: `Q_extra`, `Q_dropped`, `Q_interface`, …).
  `body.py` keeps `ThermalBody`, `MeltingBody` and everything after; `thermal/__init__.py` keeps `SOLVER_NAMES` and
  `thermal_solver`. Both re-export the moved names, so `body.ConstantBody is shared.body.ConstantBody`,
  `thermal.StepResult is shared.thermal.StepResult`, and `from . import SIGMA_SB, MissingBackend, StepResult` in the
  backends keeps working. Adapt from `git show da4c549:reentry_model/shared/body.py` and `…/shared/thermal.py`, but
  take the bodies from the ported files, not from da4c549.

- [ ] **Step 1: Write the failing test** `test_body_and_thermal_reexport_the_shared_names`: the identities above for
  every moved name, and `pickle.loads(pickle.dumps(StepResult(...)))` round-trips. Run; expect FAIL.
- [ ] **Step 2:** Move the names; re-export them with `# noqa: F401` imports at the top of `body.py` and
  `thermal/__init__.py`. Keep `thermal/__init__.py`'s lazy backend imports exactly as they are (the FEniCSx import must
  still raise `MissingBackend` only when chosen).
- [ ] **Step 3:** Run the test, `tests/test_reentry_model_body*.py`, `tests/test_reentry_model_thermal*.py`; expect
  PASS. `same_results.py` `step1,element,skins,nomelt_elements` into `output/refactor_check/task5/`, `compare` →
  `identical`. Commit.

---

### Task 6: The common run-command parts and `--architecture`

**Files:**
- Create: `reentry_model/shared/cli_common.py`
- Modify: `reentry_model/cli.py`, `tests/test_reentry_model_cli.py` (append), `tests/test_reentry_model_shared.py`

**Interfaces:**
- Produces, in `shared/cli_common.py` (moved from `cli.py` with values unchanged): `DEFAULT_OUTDIR` (`= OUTPUT_DIR`),
  `DEFAULT_HEADING_DEG`, `DEFAULT_LAT_DEG`, `DEFAULT_LON_DEG`, `DEFAULT_EPOCH`, `DEFAULT_MATERIAL_DENSITY`,
  `ATMOSPHERES`, `WINDS`, `parse_epoch(text) -> datetime`, `make_atmosphere(spec, epoch, wind_name)`,
  `git_commit() -> str | None`, `provenance() -> dict`, and two new ones:
  `add_flight_arguments(parser: argparse.ArgumentParser) -> None` — adds exactly the arguments `cli.build_parser`'s
  `run` parser adds from `--diameter` through `--name`, in the same order, with the same types, defaults and help
  strings (moved, not rewritten); `ARCHITECTURES = ("skyscrapers", "levelset")`, `DEFAULT_ARCHITECTURE = "skyscrapers"`,
  `architecture_of(run_argv: list[str]) -> str` — reads `--architecture X` or `--architecture=X` from the arguments
  after `run`, default `DEFAULT_ARCHITECTURE`, raises `SystemExit(2)` with the `run` usage on an unknown value or a
  missing one (adapt from `git show da4c549:reentry_model/shared/cli_common.py`).
  `cli.py` re-exports every moved name (`cli.make_atmosphere is cli_common.make_atmosphere`; Spheral needs it).
  `model_run_name` stays in `cli.py` (the port made it skyscrapers-specific: melt, skin and spray suffixes).
- Behaviour: `run` gains `--architecture {skyscrapers,levelset}` (default `skyscrapers`, added right after `--name` by
  `cli.py`, help `"model architecture (default %(default)s); levelset is not built yet"`). With `skyscrapers` the run is
  exactly today's. With `levelset`, `main` prints
  `ERROR: the level-set architecture is not built yet (design: docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md, stage 1 onwards)`
  to stderr and returns 2 before any model work. `resume` and `compare` do not take the flag (argparse rejects it,
  exit 2). Run names are unchanged (the `_arch-levelset` suffix belongs to the level-set driver, later).

- [ ] **Step 1: Write the failing tests** (adapt from da4c549's appended tests):
  `test_architecture_levelset_exits_2_until_built(form)` parametrised over `["--architecture", "levelset"]` and
  `["--architecture=levelset"]`: `cli.main(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.5", *form])
  == 2` and stderr contains `not built yet`; `test_architecture_skyscrapers_is_the_default`: the parsed `run` namespace
  has `architecture == "skyscrapers"`; `test_architecture_only_on_run`: `cli.main(["resume", "x", "--architecture",
  "skyscrapers"])` and the same for `compare` raise `SystemExit` with code 2; in the shared test file,
  `test_cli_common_reexports_and_paths` (the identities, `DEFAULT_OUTDIR == reentry_model.OUTPUT_DIR`) and
  `test_spheral_entry_points_resolve` (`reentry_model.cli.make_atmosphere` is callable;
  `reentry_model.material.Material` has `from_drama_json`; a material built from `data/materials/` has the private
  tables `_T_h`, `_h_nodes`, `_latent_slope`). Also `test_run_help_lists_the_same_flags`: the list of option strings of
  the `run` parser equals BASE's list plus `--architecture` — record BASE's list from `git show BASE:Reentry\ Model/reentry_model/cli.py`
  by running it once before the change and pasting the list into the test as a literal. Run; expect FAIL.
- [ ] **Step 2:** Create `shared/cli_common.py`, move the names, have `build_parser` call `add_flight_arguments(r)` where
  the moved block was, add `--architecture`, route it in `main`.
- [ ] **Step 3:** Run the new tests and `tests/test_reentry_model_cli.py` whole; expect PASS. `same_results.py`
  `step1,element,nomelt_elements` into `output/refactor_check/task6/`, `compare` → `identical`. Commit.

---

### Task 7: Independence, documentation and acceptance

**Files:**
- Modify: `tests/test_reentry_model_shared.py`; repository-top `README.md`; `Reentry Model/CLAUDE.md`;
  `docs/model_assumptions.md`; the design's §14 item "0+S" (one sentence: built, plan path, branch)

**Interfaces:**
- Consumes: everything above.

- [ ] **Step 1: Write the tests** (adapt from da4c549): `test_shared_package_imports_neither_architecture` (in a
  subprocess, import every `reentry_model.shared.*` module and assert none of `reentry_model.{skin,skinbody,columns,
  movemesh,mesh,spreading,film,restart,coupled,cli,body}` and no `reentry_model.thermal.*_backend` is in
  `sys.modules`); `test_shared_modules_import_nothing_outside_shared_statically` (parse each `shared/*.py` with `ast`;
  every relative import is `.`-level inside `shared` or `..` naming only the Global Constraints' root constants);
  `test_every_reentry_model_import_in_tests_and_analysis_resolves` (walk `tests/`, `analysis/` and
  `../Spheral/spheral_frag/fe.py` with `ast`, collect every `reentry_model[.x]` import and `from reentry_model[.x] import
  y`, assert each resolves and each `y` is an attribute). Run; they should PASS at once — if one fails, it is a real
  gap from Tasks 2–6: fix it there, not in the test.
- [ ] **Step 2: Documentation.** `Reentry Model/CLAUDE.md` and the repository `README.md`: a short section "Two model
  architectures and the shared package" — what is in `reentry_model/shared/`, the alias rule, that new code imports
  from `shared` directly, that the skins are never shared, the `--architecture` flag (adapt da4c549's README text to
  the ported layout). `docs/model_assumptions.md`: one paragraph stating the refactor changed no model behaviour and
  how that was shown. No other document changes.
- [ ] **Step 3: Acceptance runs** (mains power, nothing editing):
  (a) `same_results.py run` all five cases into `output/refactor_check/after/`, `compare` with `before` →
  `identical: N files compared`, exit 0 (N recorded in the report);
  (b) unit tier → BASE's count plus the new tests, all passing (1 skipped, FEniCSx);
  (c) both reference commands of Task 1 Step 7 → BASE's outcome, test for test;
  (d) droplet and wake tests: `"$PY" -m pytest reentry_model/droplet_reentry/tests reentry_model/wake_column/tests -m "not reference" -q -p no:cacheprovider`
  → 183 passed (the port's count);
  (e) Spheral: from the repository top, `cd Spheral && "$PY" -m pytest -m "not drama and not reference" -q -p no:cacheprovider`
  → 356 passed, 95 skipped, 1 failed (`test_regenerated_sphere_run_equals_committed_bitwise`, defect 4), as on the port.
  Any other difference: stop and report.
- [ ] **Step 4:** Commit the tests and documentation. Do not merge or push: the controller reports to Asha.
