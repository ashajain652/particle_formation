# Skyscrapers Subpackage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the code that only the skyscrapers architecture uses into `reentry_model/skyscrapers/`, next to
`reentry_model/shared/` and `reentry_model/levelset/`, without changing any result.

**Architecture:** The same mechanism as stage 0 + S. Each moved module is `git mv`-ed verbatim; the only edits are
import lines. Every old path stays importable as an alias file that replaces itself in `sys.modules` with the moved
module, so tests, analysis scripts, Spheral and the saved flight states (which pickle `reentry_model.skinbody.*` etc.)
keep working unchanged.

**Tech Stack:** Python 3.12 (conda `drama_env`), pytest, the `same_results.py` harness.

**Spec:** `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4.1-4.4. This plan **amends**
§4.2's sentence "The skyscrapers-specific code stays where it is (...)": Asha, 2026-11-11, asked for a folder for
skyscrapers-only code. Task 4 rewrites that sentence. Follow-on to
`docs/superpowers/plans/levelset/00s-combined-refactor-shared-package.md` (merged to main, c9d90df).

## Decisions (Asha, 2026-11-11)

- **D1 Folder name:** `reentry_model/skyscrapers/`.
- **D2 What moves to `skyscrapers/` (9 files):** `body.py`, `coupled.py`, `columns.py`, `mesh.py`, `movemesh.py`,
  `restart.py`, `skin.py`, `skinbody.py`, `spreading.py`. **The `thermal/` package stays where it is** (Asha). None of
  them is imported by `levelset/`, `shared/`, `droplet_reentry/` or `wake_column/` (AST scan, 2026-11-11).
- **D2b `film.py` is split** (Asha): the lubrication formulas are shared physics, the film's transport belongs to the
  element model.
  - `lubrication()` (the thin and thick lubrication branches: surface speed, flux, shear rate, branch mask) and
    `deep_flux()` (the pressure-driven flux of the liquid beneath the film) go to **`shared/lubrication.py`**, with the
    parts of the module docstring that derive them. Both the element model and the skins call `lubrication`.
  - `Runoff` (the implicit upwind transport of film mass on the patch graph, used only by the element model,
    `body.MeltingBody`, `--surface-model elements`) goes to **`skyscrapers/elements/film.py`**, with the docstring's
    runoff paragraph. It re-exports `lubrication` and `deep_flux` from `shared.lubrication`, so `reentry_model.film`
    (its alias) still offers all three names and `body.MeltingBody`, which calls them through its film module, is
    unchanged. `skinbody` switches to `from ..shared import lubrication`.
  - **`elements/` sits inside `skyscrapers/`** (Asha, 2026-11-12): `reentry_model/skyscrapers/elements/`, the element
    model's own code (the surface resolved as a grid of patches, run as `--surface-model elements`). Only `film.py`
    moves into it now; `MeltingBody` stays in `skyscrapers/body.py`.
- **D3 `cli.py` and `__main__.py` stay at the top** as the one command line of both architectures; the skyscrapers
  run driver moves later, when the level-set run driver exists.
- **D4 Old paths stay as aliases**; tests and analysis scripts keep their old imports (clean-up tabled).

## Global Constraints

- Paths are relative to `Reentry Model/` of the worktree. Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`.
- **Behaviour-preserving.** No physics, default, option, help text, run name, history column, JSON key, output file
  name or exit code changes. Proof: `analysis/skyscrapers_port/same_results.py compare` against Task 1's baseline
  gives `identical`, exit 0. If anything differs, stop and report; never "fix" a result.
- Branch `skyscrapers-subpackage` from `main` at its head when execution starts (ae20dd3 on 2026-11-11), worktree
  `.claude/worktrees/skyscrapers-subpackage`. Do not push or merge.
- Moved files change only in import lines. Rule inside `skyscrapers/`: a moved sibling stays `from . import X` /
  `from .X import Y`; a shared module becomes `from ..shared import X` / `from ..shared.X import Y` (not via its
  alias); `thermal` becomes `from ..thermal import …`; a package-root constant becomes `from .. import OUTPUT_DIR`.
- `shared/` and `levelset/` must not import `skyscrapers/`.
- Do not touch `DRAMA Sphere Experiments/`, `Spheral/`, `reentry_model/droplet_reentry/`, `reentry_model/wake_column/`,
  `reentry_model/levelset/`, `reentry_model/thermal/`, `prototype/`. Never write to the resume state file below or to
  the sibling worktrees.
- Resume state (read only): `<repository top>/.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl`.
  It was pickled by pre-refactor code, so the `resume` case is the real old-state test.
- Known defects 1-5 of `docs/superpowers/plans/2026-11-11-skyscrapers-refactor-handoff.md` stay unfixed.
- Unit tier: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -n 4 -p no:cacheprovider -q`;
  run on mains power; do not edit the package while a measurement runs.
- Commit messages: one sentence on what moved, blank line, `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A saved state pickled before the move (old module paths) must resume to identical output — the `resume` case.
2. Tests that monkeypatch an old path (`skinbody.SkinBody`, `coupled.CoupledRun.run`, `film.lubrication`) must reach
   the code that runs (for the film, Task 3 Step 2) — guaranteed by `sys.modules` aliases (not re-exports); Task 2's identity tests pin it.
3. Spheral's `fe.py` imports `reentry_model.mesh` by name — the import-resolution test of
   `tests/test_reentry_model_shared.py` already covers Spheral; it must still pass.
4. The no-skyscrapers-modules test in `tests/test_reentry_model_shared.py` bans old names only; after the move it would
   pass vacuously for the new names — Task 2 extends it.
5. Importing a moved module by its new path first, then its old path (or the reverse), must give one module object —
   Task 2's subprocess identity test.

---

### Task 1: Worktree, plan commit and baseline

- [ ] **Step 1:** Create the worktree and branch from `main`; commit this plan file.
- [ ] **Step 2:** `"$PY" analysis/skyscrapers_port/same_results.py run --code "$PWD" --out output/skyscrapers_move/before --cases step1,element,skins,resume,nomelt_elements --state "<resume state>"`.
  Expected: exit 0, about 5 minutes.

### Task 2: The nine modules

**Files:**
- Create: `reentry_model/skyscrapers/__init__.py` (docstring only: "Code only the skyscrapers architecture uses; see
  the README section on the two architectures."); `tests/test_reentry_model_skyscrapers.py`
- Move: `reentry_model/{body,coupled,columns,mesh,movemesh,restart,skin,skinbody,spreading}.py` → `reentry_model/skyscrapers/`
- Create (aliases): the nine old paths, same form as the shared aliases, docstring "Moved to
  reentry_model.skyscrapers.<name> (skyscrapers subpackage, 2026-11-11)…", body
  `from .skyscrapers import <name> as _moved` / `sys.modules[__name__] = _moved`.
- Modify: `reentry_model/cli.py` import lines only (`from .skyscrapers import body, coupled, mesh, movemesh, restart, skin, skinbody`);
  the banned set in `tests/test_reentry_model_shared.py`.

- [ ] **Step 1: Failing tests** in `tests/test_reentry_model_skyscrapers.py`: `MOVED = [the nine]` (the film's alias points at `skyscrapers.elements.film`; Task 3 tests it);
  `test_old_path_is_the_moved_module(name)` asserts
  `import_module("reentry_model." + name) is import_module("reentry_model.skyscrapers." + name)`; the same in a
  subprocess importing the new path first; `test_old_pickle_path_loads_the_moved_class`: pickle a small instance of a
  moved class that needs no mesh (implementer picks one and names it) at protocol 0, replace
  `reentry_model.skyscrapers.` with `reentry_model.` in the bytes, unpickle, assert `type(...) is` the moved class;
  `test_shared_and_levelset_do_not_import_skyscrapers`: AST scan of `shared/` and `levelset/` finds no import resolving
  to `reentry_model.skyscrapers`. Run: FAIL.
- [ ] **Step 2:** `git mv` the nine, apply the import rule (the film is not moved until Task 3, so `body.py` and
  `skinbody.py` import it as `from .. import film as film_mod` for now), write the aliases, switch `cli.py`'s imports. Add
  `reentry_model.skyscrapers` and each `reentry_model.skyscrapers.<name>` to the banned set of
  `tests/test_reentry_model_shared.py`.
- [ ] **Step 3:** New and shared tests PASS; `same_results` for all five cases into `output/skyscrapers_move/task2`,
  `compare` each against `before` → `identical`. Commit.

### Task 3: Splitting the film

**Files:** Create `reentry_model/shared/lubrication.py` (`lubrication`, `deep_flux`, verbatim, with their docstring
paragraphs); move `reentry_model/film.py` → `reentry_model/skyscrapers/elements/film.py` (`git mv`, then cut the two
functions out, keep `Runoff` and the runoff paragraph, add `from ...shared.lubrication import deep_flux, lubrication  # noqa: F401`);
create `reentry_model/skyscrapers/elements/__init__.py` (docstring: "The element model's own code (`--surface-model elements`): the
surface resolved as a grid of patches."); old `reentry_model/film.py` becomes the alias of `skyscrapers.elements.film`. Modify
`skyscrapers/skinbody.py` (`film_mod.lubrication` → `lubrication.lubrication` via `from ..shared import lubrication`)
and `skyscrapers/body.py` (its film import becomes `from .elements import film as film_mod`). Add `lubrication` to `shared/`'s module list and
to the static scan's expected count in `tests/test_reentry_model_shared.py` (22 files becomes 23).

- [ ] **Step 1: Failing tests:** `reentry_model.film is reentry_model.skyscrapers.elements.film`;
  `reentry_model.film.lubrication is reentry_model.shared.lubrication.lubrication` (same for `deep_flux`);
  `reentry_model.film.Runoff.__module__ == "reentry_model.skyscrapers.elements.film"`; `shared.lubrication` imports
  only numpy. Run: FAIL.
- [ ] **Step 2:** Split as above. Confirm the monkeypatch in `tests/test_reentry_model_fenicsx.py` (around line 228:
  `film.lubrication`, `film.Runoff.edge_coefficients`) still reaches the element model's calls; it tests the element
  model only, which calls through its film module.
- [ ] **Step 3:** Tests PASS; `same_results` all five cases into `output/skyscrapers_move/task3`, `identical`. Commit.

### Task 4: Documentation and acceptance

- [ ] **Step 1:** README and `CLAUDE.md` "Two model architectures" sections: add `skyscrapers/` with its module list
  and the alias note (old paths work; states saved after this change pickle `reentry_model.skyscrapers.*` and need a
  checkout from after it); in `CLAUDE.md`'s module list replace `film` under "Physics shared with the other (level-set)
  architecture" by `lubrication`, and add `elements/film` (the element model's runoff) to the skyscrapers list. `docs/model_assumptions.md` §7a: one sentence. Design §4.2: replace "The
  skyscrapers-specific code stays where it is (...)" with the new location and "(moved 2026-11-11, Asha; `cli.py` and
  `thermal/` stay at the top)".
- [ ] **Step 2: Acceptance**, each once, expected numbers from main at execution start (re-measure on the base commit
  if a count differs): `same_results` five cases into `output/skyscrapers_move/after` → `identical`; unit tier (754
  passed, 1 skipped on c9d90df, plus this plan's new tests, plus main's 60 level-set tests if collected); Steps 1-2
  reference 13 passed; Step 3 reference 4 passed; droplet + wake 188; Spheral 356 passed, 95 skipped, 1 known failure.
  Commit.
