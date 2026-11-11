# Skyscrapers Port to Main Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The skyscrapers code (the layered-skin moving-mesh model and the element model) lands on `main` in `main`'s
layout, with the prototype line's and the build repository's histories, computing exactly what the prototype computes,
ready for another agent's refactor.

**Architecture:** All work is on branch `skyscrapers-port`, made from `main` in its own worktree. Merge 1 brings
`step3-rigid-substrate-scheil` (designs, plans, specs, prototype folders) with every path under `Reentry Model/`. The
build repository's 96 commits are rewritten into `main`'s paths by a tested path map and merged as unrelated history
(merge 2); a file both sides hold is resolved by a three-way merge against the old merge base, which is the build's own
layout. Port tools in `Reentry Model/analysis/skyscrapers_port/` check every stage and stay on `main` for the refactor.

**Tech Stack:** git 2.50.1 (`merge-tree`, `filter-branch`, `merge-file`, `interpret-trailers`), Python 3.12
(`drama_env`), pytest with pytest-xdist.

**Spec:** `docs/superpowers/specs/2026-11-11-skyscrapers-port-to-main-design.md` (approved by Asha 2026-11-11; with
merge 1 it moves to `Reentry Model/docs/superpowers/specs/`).

## Global Constraints

- Every task uses these shell variables. Paths contain a space: quote every one.
  ```bash
  ROOT="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution"
  W20="$ROOT/.claude/worktrees/subplan20"                            # step3-rigid-substrate-scheil
  WP="$ROOT/.claude/worktrees/skyscrapers-port"                      # skyscrapers-port (Task 1)
  BUILD="$W20/prototype/work-2026-10-28-skin-spraying/code"          # the build repository, branch skin-spraying
  PLAN_CODE="$W20/prototype/work-2026-10-28-skin-spraying-plan/code" # the planning copy (equals BUILD's head)
  REWRITE="$W20/prototype/work-2026-10-28-skin-spraying/rewritten"   # Task 4's clone (git-ignored folder)
  MEASURE="$W20/prototype/work-2026-10-28-skin-spraying/measure/port" # prototype-side outputs (git-ignored)
  S100="$PLAN_CODE/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl"   # read only
  TOOLS="$WP/Reentry Model/analysis/skyscrapers_port"
  PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
  MERGE_BASE=d16018bf6bb19944a88c19cf0987c33951ea3485     # step3/main merge base = the build repository's layout
  MAIN0=8ed75cfeeda2364f14dba2bf73cc48f58af6e72b          # main when this plan was written
  ```
- One `cd` per Bash call, never inside a compound command; otherwise `git -C`.
- Never touch: the DRAMA databases; `sphere_reentry.py`, `sphere_sweep.py`; `proto3`, `plan3` (either location);
  `.claude/worktrees/we-critical-308`; saved states and flight outputs; the main checkout `$ROOT` (another session's
  branch, uncommitted edits).
- Never edit anything under `DRAMA Sphere Experiments/`, `Spheral/`, `Reentry Model/reentry_model/droplet_reentry/` or
  `Reentry Model/reentry_model/wake_column/`.
- The build repository and the planning copy are never changed (their only new files: the mesh cache the CLI writes and
  `$MEASURE`). `step3-rigid-substrate-scheil` takes no commit after Task 1's tag; amendments to this plan from then on
  are committed on `skyscrapers-port`, at the plan's path under `Reentry Model/` once merge 1 has landed.
- No model change: no physics, default, option or numerical change anywhere. The only code edits allowed are `main`'s
  own changes reapplied and data-path repointing through the package's folder constants (`PROJECT_ROOT`, `DATA_DIR`,
  `OUTPUT_DIR`, `DRAMA_DATA_DIR` in `reentry_model/__init__.py`).
- Git: stage named paths only (`git add -- <paths>`); never `git add -A` or `.`, never bare `git stash`, never
  force-push, push nothing. Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Outputs are never committed. Port-side outputs go under `$WP/Reentry Model/output/`, which `main` ignores.
- Every test run: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest ... -p no:cacheprovider`.
- Long runs (the unit tier, the reference tier, Task 7's runs) only on mains power: `pmset -g batt | head -1` says
  `AC Power` before each, or wait.
- No new tools installed (`git filter-repo` is absent; git's own `filter-branch` is used).

## Review Focus

1. **A path from the old layout left where no test reaches it** (an analysis script's default `--outdir`, a usage line,
   a test helper's `REPO_ROOT`): none may remain. Task 5 Step 9 greps for them and runs each build analysis script's
   `--help`.
2. **An output folder created in the old layout during the runs** (`reentry_model_output/` anywhere outside
   `Reentry Model/prototype/`, or a mesh cache outside `Reentry Model/output/`): Task 7 Step 5 looks for one.
3. **The prototype's saved state failing to load in the ported package** (it is pickled with the prototype's module
   paths): Task 7's resume run is that test, and must load and run.
4. **A file left where git proposes rather than where the design puts it** (in `DRAMA Sphere Experiments/`, in
   `Spheral/`, or under the old flat `prototype/` or `docs/`): `checks.py outside` and `checks.py merge1` after each merge
   (Tasks 3, 5, 6, 7).
5. **`main` moving while the port is in progress**: the fast-forward in Task 9 refuses, and the port stops for Asha
   rather than merging new `main` commits unverified.

---

### Task 1: The starting point — checks, tags and the worktree

No code; the deliverable is two tags, a worktree and a report.

- [ ] **Step 1: The two lines are where the design says.** Run and record:
  `git -C "$BUILD" status --short` (empty), `git -C "$BUILD" rev-list --count HEAD` (96),
  `diff -rq -x .git -x reentry_model_output -x __pycache__ -x .pytest_cache "$PLAN_CODE" "$BUILD"` (no output),
  `git -C "$W20" status --short` (empty), `git -C "$W20" rev-parse main` (equals `MAIN0`).
  If `main` is not `MAIN0`, or either line moved: stop with NEEDS_CONTEXT, listing the new commits.
- [ ] **Step 2: Re-run the dry-run merge.** `mkdir -p "$MEASURE"`, then
  `git -C "$W20" merge-tree --write-tree --name-only --messages main step3-rigid-substrate-scheil > "$MEASURE/merge1_dryrun.txt"`
  (it exits 1 on conflicts). Count the `CONFLICT` lines by kind. Expected about 425: about 416 file location,
  5 rename/rename (sub-plans 06, 09, 16, 17, 18), 2 directory-rename splits, 2 content (`.gitignore`,
  `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`). Any other content conflict goes in the report
  with both sides' diffs from `MERGE_BASE`.
- [ ] **Step 3: Tags.**
  ```bash
  git -C "$W20" tag -a skyscrapers-prototype-end -m "Last commit of the skyscrapers prototype line (step3-rigid-substrate-scheil) before the port to main, 2026-11-11" step3-rigid-substrate-scheil
  git -C "$W20" fetch "$BUILD" skin-spraying
  git -C "$W20" tag -a skyscrapers-build-end -m "Last commit of the skyscrapers build repository (prototype/work-2026-10-28-skin-spraying/code, branch skin-spraying), original identifiers, before the port to main, 2026-11-11" FETCH_HEAD
  ```
  Check: `git -C "$W20" rev-parse 'skyscrapers-build-end^{commit}'` equals `git -C "$BUILD" rev-parse HEAD`, and
  `git -C "$W20" rev-list --count skyscrapers-build-end` prints 96.
- [ ] **Step 4: The worktree.** `git -C "$W20" worktree add -b skyscrapers-port "$WP" main`. Check
  `git -C "$WP" rev-parse HEAD` equals `MAIN0`.
- [ ] **Step 5: Report** the five checks, the conflict counts and the tags' commits. Nothing to commit.

### Task 2: The port tools — path map, checks, resolver

**Files** (all new, in `$TOOLS/`): `path_map.py`, `checks.py`, `resolve_both_sides.py`, `README.md`,
`test_port_tools.py`. Not collected by `pytest tests`; run explicitly.

**Interfaces — Produces:**
- `path_map.map_path(build_path: str) -> str`: a path in any commit of the build repository to its path on `main`.
- `path_map.unmap_path(main_path: str) -> str | None`: the inverse; `None` for a path no build path maps to.
- `path_map.prototype_line_path(step3_path: str) -> str`: a path of `step3-rigid-substrate-scheil` to its place after
  merge 1.
- `python path_map.py index-filter`: reads `git ls-files -s` lines on stdin and writes them with each path mapped by
  `map_path` (for `git update-index --index-info`). Imports nothing but the standard library (it runs under `-I`).
- `resolve_both_sides.three_way(ours: bytes, base: bytes, theirs: bytes) -> tuple[bytes, bool]`: `git merge-file -p`
  on the three; returns the merged bytes and whether it was clean.
- `python resolve_both_sides.py`: inside a merge in progress, resolves each add/add conflict (see Task 5 Step 3).
- `python checks.py <command>`, each printing a report and exiting 0 only on a pass:
  `merge1 --ref <tree-ish>`, `rewrite <original repo> <rewritten repo>`, `outside <tree-ish> [--allow PATH]...`,
  `same-code <build-ref> <port tree-ish> <main-ref>`, `history <ref>`.

**`map_path` rules**, first match wins (design §4 step 5; the sets are constants with a comment that they were checked
against `MAIN0`'s tree on 2026-11-11):
1. `reentry_model/data/<rest>` → `Reentry Model/data/<rest>`.
2. `data/fap_day.dat`, `data/fap_mon.dat`, `data/user_materials/<rest>` → `DRAMA Sphere Experiments/<same path>`.
3. `data/reference_runs/<name>`: with `_mAA7075_nomelt` in the name → `DRAMA Sphere Experiments/data/reference_runs/<name>`;
   otherwise (the four melting files) → `Reentry Model/data/reference_runs/<name>`.
4. `tests/fixtures/<rest>`, `tests/test_fixtures.py`, `tests/test_integration_drama.py`, `tests/test_sphere_reentry.py`,
   `tests/test_sphere_sweep.py` → `DRAMA Sphere Experiments/<same path>`.
5. `analysis/<s>` for `s` in `plot_knudsen_number.py`, `plot_max_temperature_by_diameter.py`, `plot_ohnesorge_number.py`,
   `plot_outcome_by_diameter.py`, `plot_outcome_vs_diameter_velocity.py`, `plot_style.py`, `plot_weber_number.py`,
   `summary_from_runs.py` → `DRAMA Sphere Experiments/analysis/<s>`.
6. Anything else → `Reentry Model/<path>` (this includes `.gitignore` → `Reentry Model/.gitignore`, removed in Task 5).

**`prototype_line_path` rules:** `.gitignore` → `.gitignore`; `data/reference_runs/<name>` →
`Reentry Model/data/reference_runs/<name>`; anything else → `Reentry Model/<path>` (the retired sub-plans are already
under `docs/superpowers/plans/melt-spraying-subplans/deprecated/` on this branch, so they map into
`Reentry Model/docs/.../deprecated/`).

**What each check does:**
- `merge1`: for each path `p` that `git diff --name-status -M MERGE_BASE skyscrapers-prototype-end` adds, modifies or
  renames to: `prototype_line_path(p)` exists in the tree; if `main` (`MAIN0`) did not change that file since
  `MERGE_BASE`, its content equals `skyscrapers-prototype-end:p`; otherwise it is listed as *combined*. Each path the
  branch deleted is absent at its mapped place. Fails on a missing or differing file.
- `rewrite`: the two histories have the same number of commits, paired in `git rev-list --topo-order --reverse` order;
  for every pair, the rewritten parents are the paired originals' parents, author and committer names, emails and dates
  are equal, the trailer `Build-repository-commit` (read with `%(trailers:key=Build-repository-commit,valueonly)`) is
  the original identifier, the message with that trailer line removed equals the original (trailing whitespace
  stripped), and the rewritten `git ls-tree -r` equals the original's with every path passed through `map_path` (same
  mode and blob).
- `outside`: every path outside `Reentry Model/` that differs between `MAIN0` and the tree, and every differing path
  under `droplet_reentry/` or `wake_column/`; fails unless each is in `--allow`.
- `same-code`: for each file `p` of `<build-ref>`, with `q = map_path(p)`: *same* if the port's `q` equals the build's
  `p`; *main's change reapplied* if it equals `three_way(main:q, MERGE_BASE:p, build:p)`; *main's copy kept* if `q` is
  outside `Reentry Model/` and equals `main:q`; *removed* for `Reentry Model/.gitignore`; otherwise *fix-up*, printed with
  its diff against the build's file. Fails only on a file missing from the port.
- `history`: `git blame --porcelain` over `Reentry Model/reentry_model/skinbody.py` at the ref: every line's commit
  carries a `Build-repository-commit` trailer; prints how many distinct build commits it reaches;
  `skyscrapers-prototype-end` is an ancestor of the ref; tag `skyscrapers-build-end` exists.

- [ ] **Step 1: Write the failing tests** in `test_port_tools.py` (they read refs from the shared repository):
  ```python
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
      # for every file of skyscrapers-build-end whose map_path is outside "Reentry Model/":
      #   MAIN0 has a file at map_path(p), and the build's blob equals MERGE_BASE's blob at p

  def test_the_map_is_one_to_one_and_unmap_inverts_it():
      # over every path in every commit of skyscrapers-build-end (148 distinct paths on 2026-11-11):
      #   unmap_path(map_path(p)) == p, and no two paths share a map_path

  def test_index_filter_maps_each_line():      # run_index_filter: the line piped through `path_map.py index-filter`
      line = "100644 " + "a" * 40 + " 0\treentry_model/data/us76_static_environment.csv\n"
      assert run_index_filter(line) == "100644 " + "a" * 40 + " 0\tReentry Model/data/us76_static_environment.csv\n"

  def test_three_way_takes_the_build_and_reapplies_main():
      assert three_way(b"a\nB\nc\n", b"a\nb\nc\n", b"a\nb\nc\nd\n") == (b"a\nB\nc\nd\n", True)
      assert three_way(b"a\nX\nc\n", b"a\nb\nc\n", b"a\nY\nc\n")[1] is False

  def test_the_rewrite_check_passes_a_faithful_rewrite_and_fails_a_dropped_file(tmp_path):
      # a two-commit repository in tmp_path with reentry_model/x.py and data/fap_day.dat; a clone rewritten with the
      # Task 4 filter-branch command; checks.py rewrite passes; after a third commit on the clone only that deletes a
      # file, it fails
  ```
- [ ] **Step 2:** `cd "$WP/Reentry Model"`, then
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest analysis/skyscrapers_port/test_port_tools.py -q -p no:cacheprovider`.
  Expected: all fail (modules missing).
- [ ] **Step 3: Implement** `path_map.py`, `resolve_both_sides.py` and `checks.py` to the interfaces above, and
  `README.md` (one paragraph per tool: what it checks and its command line).
- [ ] **Step 4:** the same command; expected: all pass.
- [ ] **Step 5: Commit** on `skyscrapers-port`: `git -C "$WP" add -- "Reentry Model/analysis/skyscrapers_port"`, message
  "Skyscrapers port: the port's tools -- the path map from the build repository's layout to main's, the checks of each
  stage, the three-way resolver; tested".

### Task 3: Merge 1 — the prototype line into `skyscrapers-port`

**Files:** the merge; resolution by hand of `.gitignore` and
`Reentry Model/docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`.

**Interfaces — Consumes:** `checks.py merge1`, `checks.py outside` (Task 2).

- [ ] **Step 1:** `git -C "$WP" merge --no-ff --no-commit skyscrapers-prototype-end` (conflicts expected).
- [ ] **Step 2: Place every file by `prototype_line_path`** (design §4 step 4), whatever git proposed: each file the
  branch added or changed under `prototype/` or `docs/` at `Reentry Model/<same path>` with the branch's content (the
  retired sub-plans in `.../melt-spraying-subplans/deprecated/`, and `main`'s copies of them at
  `.../melt-spraying-subplans/<file>` removed); the four melting reference files at
  `Reentry Model/data/reference_runs/` (git proposes `DRAMA Sphere Experiments/`); nothing left under a top-level
  `prototype/`, `docs/`, `data/` or `tests/`. `Reentry Model/tests/test_reentry_model_aero.py` and `_data.py`: keep git's
  merge of both sides' edits.
- [ ] **Step 3: The two content conflicts.** `.gitignore`: `main`'s file plus the branch's seven new rules with
  `Reentry Model/` in front: `prototype/work-2026-10-08-skin/probe_freeze/run_energy_dt05/`,
  `prototype/work-2026-10-12-moving-mesh/code/`, `prototype/work-2026-10-12-skin-heat/code/`,
  `prototype/work-2026-10-12-skin-heat/measure/`, `prototype/work-2026-10-12-skin-heat/extract_fine_box.py`,
  `prototype/work-2026-10-18-column-flow/`, `prototype/work-2026-10-28-skin-spraying/`. `00-shared-context.md`: both
  sides' changes since `MERGE_BASE` kept (read both diffs; each hunk keeps both texts, in date order where they
  overlap).
- [ ] **Step 4: Check before committing.** `T=$(git -C "$WP" write-tree)` (fails while anything is unmerged), then
  `"$PY" "$TOOLS/checks.py" merge1 --ref "$T"`: passes, *combined* = `.gitignore`, `00-shared-context.md`, the aero and
  data tests (plus any content conflict Task 1 Step 2 listed). `"$PY" "$TOOLS/checks.py" outside "$T" --allow .gitignore`:
  passes. `git -C "$WP" grep -nE '^(<<<<<<<|>>>>>>>) ' "$T" -- "Reentry Model/docs" .gitignore`: no output.
- [ ] **Step 5: Commit the merge**: "Merge the skyscrapers prototype line (step3-rigid-substrate-scheil, tag
  skyscrapers-prototype-end) into main's layout -- designs, plans, specs and prototype folders under Reentry Model/, the
  melting SESAM references in Reentry Model/data/reference_runs/, .gitignore and 00-shared-context.md combined".

### Task 4: The rewrite of the build repository's history

**Files:** none in either worktree; `$REWRITE` is a git-ignored scratch clone.

**Interfaces — Consumes:** `path_map.py index-filter`, `checks.py rewrite` (Task 2).

- [ ] **Step 1: A fresh clone.** `git clone --no-local --branch skin-spraying "$BUILD" "$REWRITE"`, then
  `git -C "$REWRITE" remote remove origin` (nothing can reach the build repository from it). Check its head equals
  `git -C "$W20" rev-parse 'skyscrapers-build-end^{commit}'`.
- [ ] **Step 2: Rewrite** (`cd "$REWRITE"` alone first; `export PY TOOLS`):
  ```bash
  FILTER_BRANCH_SQUELCH_WARNING=1 git filter-branch \
    --index-filter 'git ls-files -s | "$PY" -I "$TOOLS/path_map.py" index-filter | GIT_INDEX_FILE="$GIT_INDEX_FILE.new" git update-index --index-info && mv "$GIT_INDEX_FILE.new" "$GIT_INDEX_FILE"' \
    --msg-filter 'git interpret-trailers --trailer "Build-repository-commit: $GIT_COMMIT"' \
    -- skin-spraying
  ```
- [ ] **Step 3: Check.** `"$PY" "$TOOLS/checks.py" rewrite "$BUILD" "$REWRITE"`: passes (96 commits, every tree, parent,
  identity, date and message). `git -C "$BUILD" rev-parse HEAD` is unchanged.
- [ ] **Step 4: Report** the check's output and the rewritten head. Nothing to commit.

### Task 5: Merge 2 — the port

**Files:** the merge; then fix-ups in `Reentry Model/reentry_model/sesam_io.py`, the build's tests in
`Reentry Model/tests/` and the build's scripts in `Reentry Model/analysis/` (exact list from Step 9's grep).

**Interfaces — Consumes:** `resolve_both_sides.py`, `checks.py outside`, `checks.py same-code` (Task 2); `$REWRITE`
(Task 4). **Produces:** `sesam_io.MELTING_REFERENCE_DIR = os.path.join(DATA_DIR, "reference_runs")`, the melting SESAM
references (Step 3); the non-melting ones stay in `sesam_io.REFERENCE_DIR` (`DRAMA_DATA_DIR/reference_runs`).

- [ ] **Step 1:** `git -C "$WP" fetch "$REWRITE" skin-spraying:refs/heads/build-rewritten`.
- [ ] **Step 2:** `git -C "$WP" merge --allow-unrelated-histories --no-ff --no-commit build-rewritten` (add/add
  conflicts expected wherever both sides hold a file with different content).
- [ ] **Step 3: Resolve.** `cd "$WP"`, then `"$PY" "$TOOLS/resolve_both_sides.py"`: for each add/add conflict at path
  `q`, `three_way(:2:q, MERGE_BASE:unmap_path(q), :3:q)`; clean results written and staged, the rest left with
  conflict markers and listed. By hand, for each listed file: `Reentry Model/tests/conftest.py` and `helpers.py` take
  `MAIN0`'s (the build's `drama` skip hook guards only tests that moved to `DRAMA Sphere Experiments/`; the session mesh
  fixtures are the same apart from docstrings); any package module (`cli.py`, `compare.py`, `mesh.py`, `sesam_io.py`)
  takes the build's code with `main`'s path lines. A file under `DRAMA Sphere Experiments/` must come out as `MAIN0`'s.
  `pytest.ini` needs no new marker: the build's tests use `drama` only in files that moved to
  `DRAMA Sphere Experiments/`, and `reference`, which `main` declares.
- [ ] **Step 4:** `git -C "$WP" rm -q -- "Reentry Model/.gitignore"` (its rules: `__pycache__/` and `.pytest_cache/` are
  in the root `.gitignore`; `reentry_model_output/` is the old layout's output folder).
- [ ] **Step 5: Check before committing.** `T=$(git -C "$WP" write-tree)`;
  `"$PY" "$TOOLS/checks.py" outside "$T" --allow .gitignore`: passes;
  `"$PY" "$TOOLS/checks.py" same-code skyscrapers-build-end "$T" "$MAIN0"`: no file missing, *fix-up* only for Step 3's
  hand resolutions.
- [ ] **Step 6: Commit the merge**: "Port the skyscrapers code onto main: merge the build repository's history
  (prototype/work-2026-10-28-skin-spraying/code, 96 commits, rewritten into main's paths, each naming its original in a
  Build-repository-commit trailer; originals under tag skyscrapers-build-end) ...". Then
  `git -C "$WP" branch -d build-rewritten`.
- [ ] **Step 7: Write the failing check of the old layout's paths.**
  `git -C "$WP" grep -nE 'reentry_model_output|"reentry_model", "data"|import REPO_ROOT' -- "Reentry Model/reentry_model" "Reentry Model/tests" "Reentry Model/analysis" ':!Reentry Model/analysis/skyscrapers_port'`
  prints the lines to repoint; and from `Reentry Model/`,
  `"$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -x --co` shows the collection errors
  (`from helpers import REPO_ROOT`).
- [ ] **Step 8: Repoint**, with no other change:
  - the build's tests: `from helpers import REPO_ROOT` → `PROJECT_ROOT` (main's `helpers.py`); `os.path.join(…,
    "reentry_model", "data", f)` → `os.path.join(PROJECT_ROOT, "data", f)`; `reentry_model_output` → `OUTPUT_DIR` from
    `reentry_model` (`test_reentry_model_reference_melt.py`'s `OUTDIR`);
  - `reentry_model/sesam_io.py`: add `MELTING_REFERENCE_DIR` (Produces, above) with a one-line comment;
    the readers of the four melting files use it: `test_reentry_model_data.py`'s
    `test_melting_reference_runs_are_committed_and_consistent`, `test_reentry_model_cli.py`'s `MELT_100`,
    `test_reentry_model_reference_melt.py`, `analysis/melt_verification.py`, and `test_reentry_model_aero.py`'s loop over
    the reference CSVs, which lists both folders so it covers the same eight CSVs as in the prototype;
  - the build's analysis scripts: default output folders and usage lines from `reentry_model_output/...` to
    `OUTPUT_DIR` (`output/...` in usage text). Their `REPO_ROOT, "data", ...` already resolves to `Reentry Model/data/`
    and stays.
- [ ] **Step 9: Check.** Step 7's grep: no output. From `Reentry Model/`, `"$PY" -B analysis/<s> --help` exits 0 for
  `ablation_1d.py`, `girin_reference.py`, `melt_sensitivity.py`, `melt_verification.py`, `movemesh_recession.py`,
  `reentry_model_verification.py`, `reentry_model_thermal_verification.py`, `skin_plate.py`;
  `"$PY" -B -c "import importlib.util as u; s=u.spec_from_file_location('m','analysis/flow_wavelengths.py'); s.loader.exec_module(u.module_from_spec(s))"`
  exits 0 (`skin_timing.py` has neither a parser nor a main guard: the grep covers it).
- [ ] **Step 10: The unit tier** (mains power). Collect both sides:
  from `$PLAN_CODE`, `"$PY" -m pytest --collect-only -q -m "not drama and not reference" -p no:cacheprovider --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py > "$MEASURE/collect_prototype.txt"`;
  from `$WP/Reentry Model`, the same with `tests` and without the ignores, to `"$MEASURE/collect_port.txt"`.
  `comm -3` of the sorted node-id lines: only `tests/test_fixtures.py::…` lines, on the prototype side (`k` of them).
  Then `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`.
  Expected: `709 − k passed, 1 skipped` (the skip: FEniCSx). A failure is fixed only by a path repointing; anything
  else stops the task (NEEDS_CONTEXT).
- [ ] **Step 11:** the tools' tests (Task 2 Step 2's command) still pass.
- [ ] **Step 12: Commit** the named files: "Skyscrapers port: data paths in main's layout -- the build's tests read
  PROJECT_ROOT and data/, outputs go to OUTPUT_DIR, the melting SESAM references from sesam_io.MELTING_REFERENCE_DIR
  (Reentry Model/data/reference_runs; the non-melting ones stay in DRAMA Sphere Experiments); unit tier <n> passed, 1
  skipped".

### Task 6: The documentation of the new state

**Files:** Modify `Reentry Model/CLAUDE.md`, `README.md` (repository root), `Reentry Model/prototype/README.md`.
Create `README.md` in each retired folder: `Reentry Model/prototype/work-2026-10-08-skin/`,
`work-2026-10-11-moving-mesh-plan/`, `work-2026-10-12-moving-mesh/`, `work-2026-10-12-skin-heat-plan/`,
`work-2026-10-18-column-flow-plan/`, `work-2026-10-28-skin-spraying-plan/` (all under `Reentry Model/prototype/`).

- [ ] **Step 1: `Reentry Model/CLAUDE.md`.** Replace "Step 3 (melting, melt film, Girin melt spraying) is planned but
  not implemented — nothing in `reentry_model/` melts yet" and the paragraph on the generated Step 3 plan and the
  prototype workflow with the state after the port: Step 3 and the skyscrapers architecture are in the package (date,
  merge commits); new skyscrapers work goes on short branches off `main`, no prototype folder (Asha, 2026-11-11); the
  refactor of 2026-10-12 (`docs/superpowers/plans/shared/README.md`) is next and is another agent's. The module list
  gains one line per new module, taken from its docstring. The test section: `pytest tests` runs the package's tests,
  the default adds the droplet and wake tests; the unit tier's count from Task 5. A pointer to
  `analysis/skyscrapers_port/README.md`.
- [ ] **Step 2: The root `README.md`.** The `Reentry Model/` row's "Physics model (Steps 1 and 2)" and any sentence
  saying Step 3 is not implemented, brought up to date in the same terms.
- [ ] **Step 3: `Reentry Model/prototype/README.md`** gains a dated section: the six skyscrapers folders are retired as
  of the port, kept as the record; the package is the code of record; the history is reachable from the tags
  `skyscrapers-prototype-end` and `skyscrapers-build-end`.
- [ ] **Step 4: Each retired folder's `README.md`**: three to five lines: what the folder was (from its sub-plan),
  retired 2026-11-11, the code of record is `Reentry Model/reentry_model/`, do not edit.
- [ ] **Step 5: Check.** `"$PY" "$TOOLS/checks.py" outside HEAD --allow .gitignore --allow README.md` after staging;
  `docs/model_assumptions.md` and the facts unchanged (`git -C "$WP" diff --stat MAIN0 -- "Reentry Model/docs/model_assumptions.md"`
  shows nothing; they are sub-plan 22 Task 6's, after the refactor).
- [ ] **Step 6: Commit** the named files: "Skyscrapers port: the documentation of the new state ...".

### Task 7: Verification — the port computes what the prototype computes

**Files:** Create `$TOOLS/same_results.py`; add its tests to `test_port_tools.py`; a section in `$TOOLS/README.md`.

**Interfaces — Produces** (for the refactor agent, who reruns it to show the refactor changes nothing):
- `python same_results.py run --code <folder holding reentry_model> --out <dir> [--state <file>] --cases <list>`: each
  case runs `"$PY" -B -m reentry_model <args> --outdir <out>/<case> --name <case>` with the working directory `--code`
  and `PYTHONDONTWRITEBYTECODE=1`, logging to `<out>/<case>.log` the first line of `pmset -g batt` at the start and the
  end, the exit code and the wall time. Cases, exactly:
  ```
  COMMON  = --diameter 100 --velocity 7.24 --altitude 65 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem
            --heating physics --h-surface 4 --h-core 20 --melt on --temperature 893 --dt 0.025 --t-max 0.125 --frames-every 0
  skins   = run COMMON --skin-preset dev --fine-size 2            (sub-plan 22's 65 km smoke flight, measure/task4b/run_smoke.sh)
  element = run COMMON --surface-model elements
  resume  = resume <state> --hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10
            --leeward-test centroid --sideways-heat height --slip-band on --t-max 22.6 --frames-every 0   (the runbook's smoke check)
  nomelt  = run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76
            --thermal fem --heating physics --h-surface 4 --h-core 20 --t-max 20
  ```
  If the command line refuses a flag in the element case, that flag is dropped from the element case on both sides and
  the report says so.
- `python same_results.py compare <dirA> <dirB> [--report]`: compares every file of either folder by relative path,
  except `*.log`, `*.pkl`, `*.progress`, `*.vtu`, `*.vtp`, `*.pvd`. `.csv`: same header and rows, the same text in
  every cell except columns whose names start with `wall_`. `.json`: equal throughout except the keys `provenance` and
  `files`, any key starting with `wall_`, and `settings.mesh_file` and `settings.material_file`, compared by file name
  and by the SHA-256 of the files they name. `.npz`: the same arrays, equal element for element (NaN equal to NaN).
  Exit 0 when identical. With `--report`: for every differing numeric column or field, the largest absolute and
  relative difference, and exit 0.

- [ ] **Step 1: Failing tests**, then implement `same_results.py`, then pass:
  `test_compare_ignores_timers_paths_and_provenance` (two folders differing only in `wall_skins_s`, `provenance`,
  `files` and the folder of `settings.mesh_file`, with the same mesh file bytes: identical),
  `test_compare_finds_one_changed_history_cell` (one cell of a `.csv` changed in the 17th significant digit: differs),
  `test_compare_finds_a_changed_npz_array` and `test_compare_finds_a_different_mesh_file_under_the_same_name`. Commit
  the tool, its tests and the README section.
- [ ] **Step 2: The prototype's runs** (mains power):
  `"$PY" "$TOOLS/same_results.py" run --code "$PLAN_CODE" --out "$MEASURE/prototype" --state "$S100" --cases element,skins,resume`.
  Each log ends `exit 0` (the resume took about 4 minutes and 4 GB on the laptop on 2026-11-10).
- [ ] **Step 3: The port's runs:** the same with `--code "$WP/Reentry Model" --out "$WP/Reentry Model/output/port_check/port"`.
- [ ] **Step 4: Bit for bit.** `"$PY" "$TOOLS/same_results.py" compare "$MEASURE/prototype" "$WP/Reentry Model/output/port_check/port"`:
  identical. Any difference: stop, report it (NEEDS_CONTEXT); no fix without the controller's ruling.
- [ ] **Step 5: Review Focus 2.** `find "$WP" -name reentry_model_output -not -path "*/prototype/*"` prints nothing, and
  the meshes the port's runs used are under `$WP/Reentry Model/output/meshes/`.
- [ ] **Step 6: Without melting.** `mkdir -p "$MEASURE/main_snapshot"`, then
  `git -C "$W20" archive "$MAIN0" "Reentry Model" "DRAMA Sphere Experiments/data" "DRAMA Sphere Experiments/analysis" Spheral | tar -x -C "$MEASURE/main_snapshot"`.
  Run `nomelt` with `--code "$MEASURE/main_snapshot/Reentry Model" --out "$MEASURE/main_nomelt"` and with
  `--code "$WP/Reentry Model" --out "$WP/Reentry Model/output/port_check/port_nomelt"`; then
  `compare --report` of the two. Record the differences (they are expected: the prototype rewrote the thermal solver).
- [ ] **Step 7: `main`'s pinned Step 2 tests** (mains power, about 50 minutes): from `$WP/Reentry Model`,
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests/test_reentry_model_reference.py tests/test_reentry_model_reference_thermal.py -m reference -q -p no:cacheprovider`.
  All pass; a failure stops the task.
- [ ] **Step 8: The droplet, wake-column and Spheral tests, once** (design §5). From `$WP/Reentry Model`:
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest reentry_model/droplet_reentry/tests reentry_model/wake_column/tests -m "not reference" -q -p no:cacheprovider`.
  From `$WP/Spheral`: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest -m "not drama and not reference" -q -p no:cacheprovider`.
  A failing test is re-run in the matching folder of `$MEASURE/main_snapshot`: if it fails there too, it is reported
  and does not stop the port; if not, the task stops (NEEDS_CONTEXT) with the cause, and the controller applies design
  §5's rule (the smallest change inside `Reentry Model/reentry_model/` that leaves the prototype's results identical, or
  Asha). No file of those three is edited.
- [ ] **Step 9: The checks.** `"$PY" "$TOOLS/checks.py" same-code skyscrapers-build-end HEAD "$MAIN0"` (every *fix-up* is
  one of Task 5's hand resolutions or path repointings); `"$PY" "$TOOLS/checks.py" outside HEAD --allow .gitignore --allow README.md`;
  `"$PY" "$TOOLS/checks.py" history HEAD`; `git -C "$WP" diff --quiet "$MAIN0" HEAD -- "Reentry Model/reentry_model/earth.py" "Reentry Model/reentry_model/heating.py" "Reentry Model/reentry_model/atmosphere.py"`;
  `mesh.py` the same as `MAIN0`'s apart from comments (equal `ast.dump`); `aero.mach` and `aero.knudsen` the same source
  as `MAIN0`'s (`inspect.getsource`, port against `$MEASURE/main_snapshot`).
- [ ] **Step 10: Report** every number: the three comparisons, the non-melting differences, every test count, the checks'
  outputs, and the wall times.

### Task 8: The handoff note for the refactor agent

**Files:** Create `Reentry Model/docs/superpowers/plans/2026-11-11-skyscrapers-refactor-handoff.md`.

- [ ] **Step 1: Write it**, for an agent that has not seen this work: what landed and where (the merge and fix-up
  commits, the tags); the 2026-10-12 refactor it is to carry out (`docs/superpowers/plans/shared/README.md`, the
  level-set design's §14); how to show the refactor changes nothing (`analysis/skyscrapers_port/same_results.py`, its
  four cases; the port's outputs in `Reentry Model/output/port_check/`, git-ignored, on this laptop; rerun them first
  if absent); Task 7's measured results; the known defects left unfixed (a melting element-model flight cannot save its
  state, "cannot pickle module"; a non-melting element-model resume raises `AttributeError` on `source_rows`; both found
  in sub-plan 22 Task 5k); that the flight launcher `prototype/work-2026-10-28-skin-spraying-plan/flights/` still
  points at the retired planning copy, so a launcher for the package is later work; and that the branches made before
  the port (`levelset-two-interface`, `levelset-stage0`, `wake-residence`, `segmented-trail`, `janaf-vapour-pressure`,
  `we-critical-308`, `skin-fixes-flights-macstudio`) merge later on their own, meeting the moved prototype paths.
- [ ] **Step 2: Commit** it.

### Task 9: `main` (after the final review)

Done by the controller after the whole-branch review passes, not by an implementer.

- [ ] **Step 1:** `git -C "$WP" fetch . skyscrapers-port:main` (moves `main` only by fast-forward; `main` is checked out
  in no worktree). If it refuses because `main` moved past `MAIN0`: stop and bring it to Asha.
- [ ] **Step 2: Check:** `git -C "$WP" rev-parse main` equals `skyscrapers-port`;
  `git -C "$WP" merge-base --is-ancestor skyscrapers-prototype-end main` succeeds. Push nothing (`main`, the tags):
  that waits for Asha.
