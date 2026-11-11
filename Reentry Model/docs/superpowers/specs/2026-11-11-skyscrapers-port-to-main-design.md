# The skyscrapers code onto main, with its history: design

**Status: draft for Asha's review, 2026-11-11, with Asha's answers to the builder check's two questions (section 2,
items 6 and 7).** Decisions marked *Asha* were Asha's in conversation on 2026-11-10/11; those marked *controller's call*
are the writer's, for Asha to check.

## 1. Why

The skyscrapers architecture (the layered-skin moving-mesh model of a melting AA7075 sphere, with the element model
beside it) has been built in prototype copies of the package on branch `step3-rigid-substrate-scheil`, never in the
package on `main`. Since the branch split from `main` (merge base `d16018bf`), `main` was restructured on 2026-10-24
(`d260576`) into three project folders (`DRAMA Sphere Experiments/`, `Reentry Model/`, `Spheral/`), and the two lines
have drifted: the branch is 75 commits ahead of `main` and 90 behind. Asha wants the code refactored before more model
changes (*Asha*, 2026-11-11), and the refactor is to be done by another agent on `main`. So the code must first land on
`main`, in `main`'s layout, without changing what it computes and without losing its history.

## 2. What Asha decided

1. **Port the skyscrapers code onto `main`; another agent refactors it there.** The refactor is the one planned on
   2026-10-12 (*Asha*; recorded in `docs/superpowers/plans/shared/README.md` and the level-set design's §14): the code
   lands at the prototype's module paths (`reentry_model/<module>.py`), and the refactor then moves the
   architecture-independent modules into `reentry_model/shared/`, keeps every old import path working, adds the
   `--architecture skyscrapers|levelset` switch, and checks that every flight's history is unchanged. The refactor is
   not part of this work.
2. **No model change in the port.** Physics changes wait for the two 100 mm test flights now running on the Mac Studio
   (runbook `docs/superpowers/plans/2026-11-10-skin-fixes-test-flights-macstudio.md`, results to branch
   `skin-fixes-flights-macstudio`). The 50 mm flight with the fixes also waits.
3. **After the refactor, skyscrapers work happens on short branches off `main`, with no prototype folder** (*Asha*,
   2026-11-11, replacing her earlier idea of a prototype that takes a copy of the refactored code). One branch per task or
   sub-plan, merged when tested and reviewed; flights run from pushed commits; new behaviour behind options whose old
   setting reproduces earlier runs, as now.
4. **What comes to `main`**: the code, and the designs, plans and specs (`docs/superpowers/`), with the commit history
   kept: this branch is **merged** into `main` (not squashed into one commit), and the build repository's history is
   **merged into `Reentry Model/`** as the port itself (*Asha*, option (c), 2026-11-11), so that the package's files on
   `main` carry their task-by-task history.
5. **The unfinished work of sub-plan 22 is done after the refactor**, from the refactored code on branches: the base
   flights (100 mm to 120 s, 50 mm whole), the convergence repeats (macro step halved, 10 against 20 micrometre cells, 1
   against 2 mm surface spacing, over about 10 s at peak heating), the acceptance test and the comparison with the
   element model, and the documentation (README, `docs/model_assumptions.md`, the facts). The unbuilt ring-and-shedding
   plan (`docs/superpowers/plans/2026-09-27-deformation-ring-shedding.md`) is unaffected.
6. **The build repository's history is rewritten into `main`'s paths before it is merged** (*Asha*, 2026-11-11, option
   (i), chosen over keeping the original commits unchanged). A test on scratch repositories (2026-11-11) showed why:
   merging the original history with git's subtree strategy places the files correctly but leaves `git blame` and
   `git log` on a ported file stopping at the merge commit, because no original commit ever had a file at that path;
   rewriting each commit so its files already sit at the paths they will have on `main`, and then merging, gives plain
   `git log` and `git blame` that reach each task's commit. The cost is that every rewritten commit gets a new
   identifier, while the ledgers, plans and commit messages written so far cite the old ones. So each rewritten commit's
   message ends with a line `Build-repository-commit: <original identifier>`, its author and dates are kept, the rewrite
   is done in a fresh clone (the build repository itself is never rewritten), and the original history stays reachable
   on `main` under the tag `skyscrapers-build-end` (step 2), so an old identifier can be looked up either way.
7. **The five design files written beside the build** (the slip-band spraying plan, sub-plan 23's notes, the
   re-frozen-metal draft, the level-set design's §17 item 8 and the level-set handoff note) are committed on this branch
   so they come to `main` with merge 1 (*Asha*, 2026-11-11; done as `397afae`).
8. **The droplet, wake-column and Spheral tests** (*Asha*, 2026-11-11): "no need to test droplet, wake or spheral tests
   so long as all the changes stay within" `Reentry Model/reentry_model/`, the package itself. The port also changes
   files beside the package (section 4: `Reentry Model/tests/`, `analysis/`, new files in `data/`, `docs/` and
   `prototype/`, `pytest.ini`, `CLAUDE.md` and the READMEs, and the root `.gitignore`), and inside the package it changes
   two modules Spheral loads at run time (section 3). So **these tests are run once, on the merged result, before `main`
   moves** (*Asha*, 2026-11-11, chosen over relying on the reading in section 3). Nothing under
   `DRAMA Sphere Experiments/`, `Spheral/`, `droplet_reentry/` or `wake_column/` is edited, and `DRAMA Sphere
   Experiments/`, which imports nothing from the package (checked 2026-11-11), is not tested.

## 3. The two histories

- **This branch** (`step3-rigid-substrate-scheil`, pushed to `5b68cbe`; local head `397afae` after Task 5k and the five
  design files, 2026-11-11): 78 commits ahead of `main`, 90 behind. Outside
  `prototype/` and `docs/` it changes only `.gitignore`, two tests (`tests/test_reentry_model_aero.py`,
  `tests/test_reentry_model_data.py`) and adds the two melting SESAM references (four files in
  `data/reference_runs/`). Everything else is designs, plans and specs (`docs/superpowers/`) and prototype folders
  (`prototype/work-2026-10-08-skin`, `-10-11-moving-mesh-plan`, `-10-12-moving-mesh`, `-10-12-skin-heat-plan`,
  `-10-18-column-flow-plan`, `-10-28-skin-spraying-plan`, `prototype/levelset/`): about 1 200 files, 43.6 MB. The code
  in these folders is a full copy of the package per sub-plan in the old flat layout; it is the record, not the code of
  record.
- **The build repository** (`prototype/work-2026-10-28-skin-spraying/code/`, its own git, git-ignored, on this laptop
  only; branch `skin-spraying`, head `75519a5` after Task 5k; 9.4 MB): 96 commits from 2026-10-08 (`e7922eb`, "base:
  work-2026-10-07-dt-continuum copy", a whole copy of the package as it stood on 2026-10-07 with Step 3 applied) to
  2026-11-11, one commit per task and review fix: the first skins, the moving mesh (sub-plan 19), the layered skins'
  heat (20), the column flow (21), spraying and the flights (22). Laid out as the package's own repository root:
  `reentry_model/`, `tests/`, `analysis/`, `data/`, `.gitignore`. Its content equals the planning copy's `code/` file
  for file (the build protocol checks it with `diff` at every task). Its "remote" is the column-flow build repository,
  another folder on this laptop.

A dry-run merge of this branch into `main` (`git merge-tree`, 2026-11-11, before Task 5k) gave 425 conflicts: 416 file
locations (files this branch added in folders `main` moved; git proposes the moved folder), 5 rename/rename (the
retired sub-plans 06, 09, 16, 17, 18 moved to `deprecated/` here and to `Reentry Model/` there), 2 directory-rename
splits (inside `prototype/work-2026-10-08-skin/verify/scheil_plus_amendment/reentry_model/data`), and 2 content
conflicts (`.gitignore`, `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`).

`main`'s own changes since the merge base, in the paths the build repository also has, come almost all from the
restructure (`d260576`), and are about where things are found:
- package modules: `__init__.py` (9 lines: `DATA_DIR`, `DRAMA_DATA_DIR`, `OUTPUT_DIR` and friends), `fap.py`,
  `sesam_io.py`, `mesh.py`, `cli.py` (5 lines each), `compare.py` (7);
- tests: `conftest.py` rewritten (64 % similar: `PROJECT_ROOT` from a new `helpers.py`, session mesh fixtures pinned to
  `band=0`), `helpers.py` new, small edits to `test_reentry_model_aero.py`, `_data.py`, `_material.py`,
  `_reference_thermal.py`; the SESAM fixtures and the wrapper scripts' tests moved to `DRAMA Sphere Experiments/tests/`;
- analysis: three scripts edited, the plot scripts and `plot_style.py` moved to `DRAMA Sphere Experiments/analysis/`;
- data: the fap files, the non-melting SESAM references and most material files moved to
  `DRAMA Sphere Experiments/data/`; `AA7075_nomelt.json` to `Reentry Model/data/materials/`; `atdb_sphere.json` and a
  US76 table to `Reentry Model/data/`.
"`main`'s own changes" below means all of these. `main` has also added `reentry_model/droplet_reentry/` and
`reentry_model/wake_column/`, which the build repository does not have.

What the droplet, wake-column and Spheral code use of the shared modules, and whether the port changes it (read from
`main` and the build repository's head, 2026-11-11; not run):
- **The droplet model** subclasses `trajectory.Simulator` but overrides its right-hand side (`DropletSimulator.rhs`)
  with its own drag law (`particle_trajectory.drag_coefficient`), passes no drag tables, and calls `aero.mach`,
  `aero.knudsen`, `earth`, `heating` and `atmosphere`. The build changes `trajectory.py` only in the parent's own drag
  path (it reads the body's reference length, area and shape factor), which the droplet never calls, and gives
  `Simulator.advance` an optional `snap` argument whose default, 0, keeps the old step; it changes `aero.py` only in
  `drag_coefficient` (a shape factor, default 1) and the continuum table's step at Mach 1 (smoothed over Mach 0.98 to
  1.02), neither of which the droplet uses. `aero.mach`, `aero.knudsen`, `earth.py`, `heating.py` and `atmosphere.py`
  are the same on both sides.
- **The wake column** uses only the package's folder constants and the droplet model.
- **Spheral** imports `atmosphere` (the same on both sides) and `mesh`, whose build version differs from `main`'s only by
  `main`'s own path change (the mesh cache under `OUTPUT_DIR`) and a comment, which the port reapplies. It also loads
  two modules by name at run time (`Spheral/spheral_frag/fe.py`, the preparation step only): `material`
  (`Material.from_drama_json`, `T_REF`, and the private node tables `_T_h`, `_h_nodes`, `_latent_slope`, `_fl_T`,
  `_fl`) and `cli.make_atmosphere`. The port changes both modules. This is the change Spheral was written to expect:
  `main`'s `material.py` today lacks `_latent_slope`, `_fl_T` and `_fl`, so Spheral is pointed at a prototype copy of
  the package ("until Step 3 lands on main", in `fe.py`'s docstring); after the port, its default package, `main`'s,
  supplies the melting material. Six Spheral test files load the package; those that need a finite-element run skip
  unless `SPHERAL_FRAG_FE_RUN` and `SPHERAL_FRAG_FE_PACKAGE` are set. (An earlier reading of 2026-11-11 searched import
  statements only and missed these two.)

## 4. The sequence

1. **Task 5k has landed** (2026-11-11, reviewed: outer `fd9cb6f` and `45a421b`, build `8e15153` and `75519a5`; unit tier
   709 passed, 1 skipped). The conflict counts in section 3 are from before it, so the plan re-runs the dry-run merge
   at the start, and first checks with `diff` that the planning copy's `code/` and the build repository's last commit
   are identical. If either line takes another commit before the port starts, the port takes the newer head.
2. **Tags**: `skyscrapers-prototype-end` on this branch's last commit, `skyscrapers-build-end` on the build repository's
   last commit as it is, not rewritten (fetched into the outer repository, so the original identifiers stay reachable).
   After the port, neither line takes new commits (*controller's call*).
3. **A short branch off `main`, `skyscrapers-port`, in its own worktree** (*controller's call*: the main checkout is on
   another session's branch, `wake-reaction-network-figure`, with uncommitted edits; nothing is done there).
4. **Merge 1, this branch into `skyscrapers-port`.** Conflict rules:
   - every file this branch added under `prototype/` goes to `Reentry Model/prototype/` at the same relative path, and
     every file under `docs/` to `Reentry Model/docs/` at the same relative path, whatever git proposes (the prototype
     folders keep their old flat layout inside, as `main` already keeps its older prototype folders);
   - the four melting references go to `Reentry Model/data/reference_runs/` (not to `DRAMA Sphere Experiments/`, where
     `main` keeps the non-melting ones and where git would put them, because nothing there may change: section 2 item
     8), and the package reads them from there through its folder constants;
   - the retired sub-plans go to `Reentry Model/docs/superpowers/plans/melt-spraying-subplans/deprecated/`;
   - `.gitignore`: `main`'s root file, plus this branch's new rules with their paths moved under
     `Reentry Model/prototype/`;
   - `00-shared-context.md`: both sides' text combined (read both diffs from the merge base; neither side's edits lost).
5. **Merge 2, the build repository into `skyscrapers-port`** (unrelated histories). This merge is the port. It is done
   in two stages (*Asha*, section 2 item 6):
   - **5a. The rewrite.** A fresh clone of the build repository (the original stays as it is) has every commit's files
     moved to the paths they take on `main`, by the path map below, and every commit message given the closing line
     `Build-repository-commit: <original identifier>`; authors, dates and the order of commits are kept. Checks: the
     rewritten history has the same number of commits (96); for every commit, the rewritten tree, with the path map
     undone, equals the original tree file for file; and the rewritten head's files, read back through the path map,
     equal the build repository's head.
   - **5b. The merge.** The rewritten history is fetched into the outer repository and merged into `skyscrapers-port`
     with `--allow-unrelated-histories`; each conflict is resolved by the rules below.
   - **The path map** is written out file by file and checked against `main`'s tree, not taken from git's rename
     detection: run on the restructure (`git show -M d260576`, 2026-11-11), that detection pairs files with identical
     copies in other folders (it sends the build's `tests/conftest.py` to `Spheral/tests/conftest.py`, and the SESAM
     fixtures to an old prototype's copies). `reentry_model/`, `tests/`, `analysis/` and `data/` map to the same paths
     under `Reentry Model/`, with two exceptions. The package's own data files (`reentry_model/data/`: the sphere drag
     table, the US76 table, `materials/AA7075_nomelt.json`) follow `main`'s move to `Reentry Model/data/`, where `main`
     holds identical copies. Files whose content `main` already holds outside `Reentry Model/` and that the build never
     changed from the merge base (checked 2026-11-11: the two fap files, the six non-melting SESAM references as twelve
     files, the two user materials, the SESAM fixtures, the wrapper scripts' tests `test_fixtures.py`,
     `test_integration_drama.py`, `test_sphere_reentry.py` and `test_sphere_sweep.py`, and eight DRAMA analysis scripts)
     are not ported; `main`'s copies stay as they are. The merged tree leaves every file outside
     `Reentry Model/` byte-identical to `main`'s apart from the root `.gitignore` (section 2 item 8); where to put those
     files in the rewritten history is the builder's.
   - Where both sides have a file in `Reentry Model/`, the build repository's version is taken, then `main`'s own
     changes since the merge base (section 3) are reapplied to it. The wrapper scripts themselves are never touched.
   - Files `main` has and the build repository does not (`droplet_reentry/`, `wake_column/`, their tests,
     `pytest.ini`) stay as they are.
   - `tests/conftest.py` and `tests/helpers.py` are merged by hand so that both test sets run: `main`'s
     `PROJECT_ROOT` and pinned mesh fixtures kept, the build's `drama`-marker skip hook kept, the build's `REPO_ROOT`
     uses repointed to `main`'s equivalents (*controller's call*). `pytest.ini` gains any marker the build's tests use
     that it lacks; that is configuration, not computation.
   - Data paths in the build's modules, tests and analysis scripts (including modules `main` never had, such as
     `skinbody.py` and `spray.py`) are repointed through `main`'s own mechanism, the folder constants in
     `reentry_model/__init__.py`. These edits, and `main`'s own changes reapplied, are the only differences allowed
     between the build repository's last commit and the port (section 5).
   - The build repository's `.gitignore` folds into `main`'s.
6. **Fix-up commit(s)** on `skyscrapers-port`: whatever the merged code needs to run in `main`'s layout (data paths,
   output folder `Reentry Model/output/`, test paths), and the documentation of the new state: `Reentry Model/CLAUDE.md`
   (it says nothing in the package melts yet, and describes the prototype workflow), the `Reentry Model` README,
   `Reentry Model/prototype/README.md` (the skyscrapers prototype folders are retired as of the port and kept as the
   record), and a note in each retired prototype folder's own README pointing to the package. The model's documentation
   proper (`docs/model_assumptions.md`, facts) is sub-plan 22 Task 6's, after the refactor.
7. **Verification** (section 5) on `skyscrapers-port`.
8. **`main` fast-forwarded to `skyscrapers-port`, locally.** Pushing `main` and the tags waits for Asha.
9. **A handoff note for the refactor agent** (*controller's call*): what landed, where, the verification harness and
   its baselines (so the refactor can show it changes nothing), and the 2026-10-12 decision it implements.

## 5. Verification: what "the port changes nothing" means

- **Same code**: every file of the build repository's last commit, mapped as in step 5, equals the file on
  `skyscrapers-port` apart from the reapplied data-path changes; a script lists every difference and each must be one of
  those.
- **Same results as the prototype, bit for bit**, against the planning copy at the end of Task 5k, on three runs: the
  element model's short command-line melting case, a short skins flight from the command line (the 65 km smoke case at
  coarse settings), and the resume of the 100 mm state saved at 22.0 s to 22.6 s with all the fixes and the slip band on
  (the runbook's smoke check; its values are in the runbook). Every history column and every result is compared, except
  wall-clock timers, file paths and `git_dirty`-style provenance fields.
- **The package's tests pass**: `Reentry Model/tests/` (now including the prototype's tests; baseline 709 passed and 1
  skipped in the prototype after Task 5k), run with `pytest tests`, not the `pytest.ini` default, which would also
  collect the droplet and wake-column tests.
- **The droplet, wake-column and Spheral tests, once, on the merged result** (section 2 item 8), with the commands in
  `Reentry Model/CLAUDE.md`: from `Reentry Model/`, `pytest -m "not reference"` on
  `reentry_model/droplet_reentry/tests` and `reentry_model/wake_column/tests` (the droplet model's reproduction of the
  2026-10-22 sweep to 1e-9 among them); from `Spheral/`, `pytest -m "not drama and not reference"` (tests needing the
  Spheral container or a finite-element run skip, as on `main`). If a test fails, it is run on `main` to see whether it
  already failed there; a failure `main` already has is reported and does not stop the port (*controller's call*). A
  failure the port causes stops the port before `main` moves: no file of those three is edited; the smallest change
  inside `Reentry Model/reentry_model/` that leaves the prototype's results identical (for example a default on a new
  argument) is allowed and recorded; otherwise, and always for a Spheral test whose expected values rest on `main`'s
  non-melting material, the controller brings it to Asha.
- **The merged tree is checked against section 3's reading** as well: `droplet_reentry/`,
  `wake_column/`, `Spheral/` and `DRAMA Sphere Experiments/` byte-identical to `main`'s, and `aero.mach`,
  `aero.knudsen`, `earth.py`, `heating.py`, `atmosphere.py` and `mesh.py` identical to `main`'s. Long tiers on mains
  power only.
- **Runs without melting**: `main`'s current package and the ported one are compared on one non-melting case. The
  prototype rewrote the thermal solver (sub-plan 03: nodal enthalpy, lumped capacity), so these may differ; they are
  accepted if `main`'s pinned Step 2 reference and verification tests pass, and the measured difference is reported
  (*controller's call*).
- **History reachable**: plain `git log` and `git blame` on `Reentry Model/reentry_model/skinbody.py` (no `--follow`)
  reach the rewritten task commits, each naming its original identifier; the tag `skyscrapers-build-end` reaches the
  original history; this branch's commits are ancestors of `main`.

## 6. Constraints

- Paths contain a space (`University Documents /MIT`): quote everything. Interpreter
  `/Users/ashajain/miniforge3/envs/drama_env/bin/python`.
- Never touch: the DRAMA databases, the wrapper scripts `sphere_reentry.py` and `sphere_sweep.py`, `prototype/proto3`,
  `plan3`, `.claude/worktrees/we-critical-308`, saved states and flight outputs. Never commit outputs. Never force-push,
  never bare `git stash`, never `git add -A`. Commit only what the plan names; push nothing without Asha.
- Other sessions' work is not changed: the main checkout and its branch, the wake branches, `levelset-two-interface`,
  `segmented-trail`, `janaf-vapour-pressure`, the Mac Studio's results branch. Their later merges into `main` are theirs.
- The planning copy and the build repository stay where they are, unchanged, after the port; the rewrite works on a
  clone.
- No new tools are installed: `git filter-repo` is not on this laptop, and git's own `filter-branch --index-filter`
  (git 2.50.1) did the rewrite in the scratch test (*controller's call*).
- The jump-fault probe's frozen copy and outputs under `prototype/work-2026-10-28-skin-spraying/measure/` (git-ignored)
  are not touched.
- Laptop: long runs and the unit tiers on mains power only (`pmset -g batt`).

## 7. Goals, non-goals, anti-goals

- **Goals**: the skyscrapers code in `main`'s package with its task history; the designs, plans and specs on `main`;
  `main`'s other projects unaffected; proof that the port computes exactly what the prototype computes.
- **Non-goals**: the refactor; any physics change; the remaining sub-plan 22 work; the level-set code; merging other
  branches.
- **Anti-goals** (a failure even if it works): a port that changes any skyscrapers or element-model result; losing
  either history; breaking the droplet, wake or Spheral code; leaving `main` with two live copies of the skyscrapers code
  that future work could edit by mistake (the prototype folders are marked retired).

## 8. Left to the builder

The exact git commands for the rewrite and the merge, within step 5 (for example an index filter that moves each path
by the map and a message filter that appends the original identifier, then `git merge --allow-unrelated-histories`);
how the path map is written down and checked; the order of the fix-up commits; the comparison scripts; the wording of the README and CLAUDE.md updates
within what section 4 says.
