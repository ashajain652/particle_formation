# Deposit spreading and the resume command: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`)
> syntax for tracking. Asha chose subagent execution and approved the design on 2026-11-06.

**Goal:** A skin flight can resume from a saved state from the command line, and a pass limits the deposit's fold by
moving metal being handed down sideways to lower neighbours, where it lands hot.

**Architecture:** Task 1 is sub-plan 22's Task 5b: a `resume` command over `CoupledRun.run(resume=...)`. Task 2 is sub-plan
22's Task 5c: a pure pass in `reentry_model/spreading.py`, called from `SkinBody._move_base` between the hand-down and the
corners' move. Metal that lands on a colder corner crosses the base at that node's enthalpy, and the rest of its heat
goes to that node through `SkinBody.pending`, the existing deferred per-corner heat, so it arrives at the corner rather
than being shared over the whole moving zone by `MovingMesh.move`'s uniform correction. Task 3 records it.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, scikit-fem, gmsh, pytest with xdist.

**Spec:** `docs/superpowers/specs/2026-11-06-deposit-spreading-design.md`, and
`docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` sections 5.3, 6.6, 7.2, 10 and 12. Measurements behind
the design: `prototype/work-2026-10-28-skin-spraying/measure/task5/snapshots/findings.md`.

## Global Constraints

- Paths contain a space (`University Documents /MIT`): quote every path. Interpreter:
  `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`.
- Outer repository: this worktree, branch `step3-rigid-substrate-scheil`. Commit locally there after each task and **never
  push**. Stage named paths only, never `git add -A`. Never commit `docs/superpowers/plans/melt-spraying-subplans/23-notes.md`
  or anything under `reentry_model_output/`. Never use bare `git stash`, never force anything.
- Never touch: the DRAMA databases, the wrapper scripts `sphere_reentry.py` and `sphere_sweep.py`, the package at the
  repository root, `prototype/proto3/`, `plan3`, `.claude/worktrees/we-critical-308`, the saved states and frames under
  `code/reentry_model_output/task5_snapshots`, or the generated `22-skin-spraying-flight.md` by hand.
- Code is written and tested in the planning copy `prototype/work-2026-10-28-skin-spraying-plan/code/`. Every task ends by
  following **Finishing a task** below.
- Unit tier, from the planning copy's `code/`, on mains power:
  `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py -p no:cacheprovider -n 4`.
  Baseline 541 passed, 1 skipped, about 11 minutes. No failure is acceptable. Run it once per task, before committing.
- Books: in every skin test and run the mass is closed to 1e-12, the energy balance residual is below 1e-8, and the
  momentum books close. Bad settings exit 2 (`parser.error`); a model failure raises `RuntimeError` (exit 1).
- Without the new options every output of every run is as before to the last bit, column set included.
- Repository style: comments and docstrings say what was decided, by whom and when, with measured numbers; test names are
  full sentences; commit messages are one long descriptive line in the existing pattern, ending with the trailer
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

### Finishing a task

1. **Before your first edit** to a file another task already changed, keep its state: for each file you will change, find
   the last task in `TASKS` in `make_plan22.py` that changes it, and if `states/task<that>/<rel>` does not exist, copy the
   file there unchanged. `make_plan22.py` refuses a replay that does not reproduce `code/`.
2. Add your task to `TASKS`, and a task section to `part_22.md` immediately after the previous task's, before
   `## The outline`, in the pattern of Task 5a's: Files, Interfaces, what was decided, Measured, Rulings, Steps, with
   `{{KEY}}` markers for your diffs. Update the document's opening note ("Written at code level") and the outline's heading.
3. Run `"$PY" make_plan22.py` from the planning copy; it must print its check as passed. It regenerates the sub-plan.
4. Bring the build repository `prototype/work-2026-10-28-skin-spraying/code/` (git-ignored, branch `skin-spraying`, head
   dfe3b66 before these tasks) level with the planning copy's `code/` for your files, confirm `diff` finds no difference
   in them, and commit there with your task's message.
5. Commit in the outer repository: your `code/` files, `states/`, `part_22.md`, `make_plan22.py`, the regenerated
   `docs/superpowers/plans/melt-spraying-subplans/22-skin-spraying-flight.md`, and `flights/fly.py` if changed.

## Review Focus

1. **A state saved before these tasks.** The two real states were pickled by commit 055c313, so a `SkinBody` from them has
   no spreading attributes, and their saved history row has no spreading columns. They must resume, with and without a
   limit. Class-level defaults on `SkinBody` give the attributes; Task 2's test deletes them from an instance to pin this.
2. **Column sets across a resume.** `CoupledRun.run` takes its columns from the first row, the saved one. Resuming with a
   limit must add the spreading columns to it, and resuming with `--spread-limit off` from a state that had one must drop
   them, or the later rows are lost or raise `KeyError`. Pinned in Task 2.
3. **A donor that gives away all its hand-down.** `h_outward` is nan there, and the reconciliation and fault count must
   still hold. Pinned in Task 2.
4. **Heat in transit.** `pending` is overwritten by `self.pending = mismatch` in `_advance_part` and read as a load in the
   next part. The landing heat is added after that line, in `_move_base`, and a part that is split or repeated returns or
   restores before it, so the heat is applied once. The books count `pending`, so the last part of a flight stays closed.
   Pinned by the hot-landing test's energy residual.
5. **A rebuild inside the move.** `move_with_policy` may rebuild the 3D mesh; the corners and so `pending` are unaffected,
   but `T_node` and the node enthalpies must be read before the move. Pinned by reading them with the existing `T_node`.

---

### Task 1: The resume command (sub-plan 22 Task 5b)

**Files:**
- Modify: `reentry_model/cli.py` (the `resume` subcommand, `cmd_resume`; `cmd_run`'s writing of its files shared, not copied)
- Modify: `reentry_model/coupled.py` only if `run(resume=...)` needs it
- Modify: `flights/fly.py` (beside `code/`): `fly.py <outdir> <name> resume <state> [options]` runs the resume command
- Test: `tests/test_reentry_model_restart.py`
- Create (git-ignored): `prototype/work-2026-10-28-skin-spraying/measure/task5b/`, the production check

**Interfaces:**
- Consumes: `restart.load`, `restart.FORMAT`, `CoupledRun.run(resume=state)`, `spray.SourceSpool`, the state dict of Task 5a
  (`format, run, step, row, time_s, rng, meta{name, args}`).
- Produces: `python -m reentry_model resume STATE --outdir DIR --name NAME [--t-max S] [--dt-continuum S] [--seed N]
  [--frames-every N] [--slurry-spray on|off]`, writing the same files as `run`: `<outdir>/<name>.csv`, `.json`, and the
  droplet files under `<outdir>/<name>/`. The JSON gains `"resumed_from": {"state": path, "step": int, "time_s": float,
  "run_name": str, "changed": {option: value}}`. `resume_arguments(parser)`-style wiring so Task 2 can add one option.

**Behaviour:**
- Random state: the saved one, unless `--seed` is given, which seeds numpy's generator instead.
- `run.settings`: `save_at` and `save_path` cleared, `output_dir` set to `<outdir>/<name>/vtk`, `frames` emptied so the new
  series lists only its own frames, `frames_every` set when given (with `run._frame_dt` and `run._next_frame` made
  consistent). `--t-max` sets `run.sim.settings.t_max`.
- `--dt-continuum`: sets `settings.dt_continuum`, and `settings.dt` too if the switch has latched (`run.switched_at`).
- Droplet records go to a new `SourceSpool` in `<outdir>/<name>/source_rows`. Only records after the resume are written.
- `--slurry-spray` sets `run.body.slurry_spray`. The 20 um cell repeat is not part of this task (no flight needs it yet).
- Exit 2 for: a missing state file, a format other than `restart.FORMAT`, `--t-max` not after the saved time,
  `--dt-continuum` not in (0, the run's `dt`]. A non-finite value is refused too.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_restart.py`, using its `FLIGHT` and helpers:
  `test_the_resume_command_repeats_the_flight_to_the_last_bit_and_writes_its_files` (history rows equal the unbroken flight's
  from the saved step, droplet records after the saved time equal, `resumed_from` correct, `changed` empty);
  `test_the_resume_command_changes_what_it_is_told` (`--t-max` flies further, a halved `--dt-continuum` takes twice the
  steps, `--seed` differs from the saved random state and is recorded in `changed`);
  `test_the_resume_command_refuses_what_it_cannot_resume` (parametrized on the four exit-2 cases above).
- [ ] **Step 2: Run them to see them fail** with the subcommand unknown.
- [ ] **Step 3: Implement** `cmd_resume`, moving `cmd_run`'s writing of the csv, json and droplet files into a function both
  call, so `run`'s outputs do not change. Rebuild what the writer needs from `state["meta"]["args"]` with the overrides.
- [ ] **Step 4: `flights/fly.py`** routes `resume`, its rows and progress kept as for `run`.
- [ ] **Step 5: Run the tests to see them pass**, then the unit tier: 544 passed, 1 skipped expected at least.
- [ ] **Step 6: The production check** in `measure/task5b/`, on mains power, through `flights/launch.sh`: resume
  `code/reentry_model_output/task5_snapshots/d050__physics__skin.state.pkl` to 172.0 s and `d100__physics__skin.state.pkl`
  to 23.0 s, into `code/reentry_model_output/task5b_check`, and compare with the snapshot flights' `rows.csv` at full
  precision: no column but the `wall_` timers may differ. Record what was compared and the wall time in `measure/task5b/check.md`.
- [ ] **Step 7: Finish** (Global Constraints, Finishing a task), with diffs keyed `*_T5B`.

### Task 2: Deposit spreading, landing hot (sub-plan 22 Task 5c)

**Files:**
- Create: `reentry_model/spreading.py`
- Modify: `reentry_model/skinbody.py` (`SkinBody.__init__`, class defaults, `_move_base`, `skin_stats`)
- Modify: `reentry_model/cli.py` (`--spread-limit` on `run` and `resume`, the name suffix, the settings' record)
- Modify: `reentry_model/coupled.py` (the results; the column list's note)
- Test: `tests/test_reentry_model_spreading.py` (new), `tests/test_reentry_model_skinbody.py`,
  `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: `Columns` (`edges`, `_tri`, `_opp`, `points`, `normal`, `area`, `faces`), `SkinBody._hand_cap()`, `Task 1`'s resume wiring.
- Produces:
  ```python
  SPREAD_PASSES = 4       # passes a part at most: metal travels at most a ring or two of 2 mm neighbours in one part
  SPREAD_SHARE = 0.25     # of the fold's excess over a donor's four pairs

  @dataclass
  class Spread:
      handed: np.ndarray    # kg per column, the hand-down after spreading; the argument itself if the limit does not bind
      enthalpy: np.ndarray  # J per column, likewise
      moved: float          # kg moved between columns
      unplaced: float       # kg still wanting to move after the last pass
      over: int             # edges over the limit on the first estimate

  def fold_angles(columns, points) -> (angle_deg, convex)   # per edge; convex: a ridge along the edge, else a valley
  def spread_hand_down(columns, points, normal, area, rho, drawn, handed, enthalpy, cap, limit_deg,
                       passes=SPREAD_PASSES) -> Spread
  ```
  `SkinBody(..., spread_limit=None)`; `--spread-limit DEG|off` on `run` and `resume`, default off, DEG in (0, 180), exit 2
  otherwise or without a melting skin run; name suffix `_spread-<DEG>`; history columns present only with a limit:
  `spread_mass_kg` (since time zero), `spread_unplaced_kg` and `spread_edges_over` (at the last move),
  `spread_jump_max_K` (largest landed metal over node temperature since the last row); results `spread_limit_deg`,
  `spread_mass_kg`.

**The pass** (an algorithm the signature does not determine):
```python
h, H = handed.copy(), enthalpy.copy()
floor = np.maximum(drawn - cap, 0.0)                 # a donor keeps what balances its draw, so the part's move ratio holds
for each pass:
    q = points - normal * ((drawn - h) / (rho * area))[:, None]     # the base as it would stand; absorption and the
    angle, convex = fold_angles(columns, q)                         # smoothing are left out of this estimate
    over = np.flatnonzero(angle > limit_deg)
    if not over.size: break
    # four pairs per folded edge (a, b) with opposite corners (c, e), skipping an open boundary's -1:
    # convex: a and b give to c and e; valley: c and e give to a and b
    dz = einsum((q[donor] - q[receiver]), normal[donor]); keep pairs with dz > 0
    want = SPREAD_SHARE * (angle[pair] - limit_deg) / angle[pair] * rho * area[donor] * dz
    want scaled so no donor falls below floor and no receiver above cap (np.bincount over donors, then receivers)
    move want with the donor's mean enthalpy H[d] / h[d]
```
After the last pass `unplaced` is the sum of `want`, scaled by the donors' limit only, on one more evaluation. Sums of `h`
and `H` are kept to round-off. If no edge is over the limit at the first evaluation, return the arguments themselves.

**In `_move_base`**, when `spread_limit` is set, after `hand_down` and the `held` lines:
```python
s = spreading.spread_hand_down(col, col.points, col.normal, col.area, mat.rho, book.drawn, handed, handed_H,
                               self._hand_cap(), self.spread_limit)
handed_s, landed_H = s.handed, s.enthalpy
received = handed_s > handed                         # metal landed here from a neighbour: at this node's enthalpy,
excess = np.where(received, landed_H - handed_s * mat.enthalpy(T_node), 0.0)      # the rest of its heat to this node
landed_H = landed_H - excess
```
`moved`, `inward`, `h_outward` and the reconciliation read `handed_s` and `landed_H`; `jump_faults` keeps the original
`handed`, which is what left each skin; `handed_mass` keeps the original total. After the reconciliation:
`self.pending = self.pending + excess`. Without a limit this code path is not entered, so `handed_s is handed`.

- [ ] **Step 1: Write the failing tests.** In `tests/test_reentry_model_spreading.py`, on a small triangulated plate:
  `test_a_limit_that_does_not_bind_returns_the_hand_down_itself` (`s.handed is handed`, `moved == 0`);
  `test_a_ridge_hands_its_metal_to_its_lower_neighbours_and_keeps_the_sums` (sums to 1e-12 of the total, donors fall, only
  lower corners rise); `test_a_valley_fold_gives_from_the_rim_to_the_floor`;
  `test_repeated_parts_hold_a_growing_ridge_near_the_limit` (a crest handing down every part stays within 5 degrees of a
  60 degree limit over 200 parts, and exceeds 90 degrees with the pass off);
  `test_a_receiving_corner_stays_within_its_move_allowance`; `test_a_donor_keeps_what_balances_its_draw`;
  `test_what_cannot_be_placed_stays_and_is_counted`.
  In `tests/test_reentry_model_skinbody.py`, on the melting test body the spraying tests use:
  `test_a_limit_that_does_not_bind_changes_nothing_to_the_last_bit` (`spread_limit=179.9` against none, every shared column);
  `test_spreading_keeps_the_books` (a limit low enough to bind: mass 1e-12, energy 1e-8, momentum, `spread_mass_kg > 0`);
  `test_metal_landing_hot_takes_its_heat_to_the_receiving_corner_and_counts_no_fault` (a receiver colder than the metal by
  over 5 K: `pending` rises there by the excess, `jump_faults` unchanged, `spread_jump_max_K > 5`, energy closed);
  `test_a_donor_that_gives_away_all_its_hand_down_keeps_the_books`;
  `test_a_body_saved_without_spreading_attributes_flies_on` (delete the instance attributes, then advance).
  In `tests/test_reentry_model_cli.py`: the option's parsing, name suffix, record and refusals (parametrized: 0, 180, nan,
  without melting). In `tests/test_reentry_model_restart.py`:
  `test_resuming_with_a_limit_adds_the_spreading_columns_and_off_drops_them`.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement `spreading.py`**, then the body, then the command line. Comment the hot landing with Asha's decision
  of 2026-11-06 and why `pending` rather than `h_outward`.
- [ ] **Step 4: Run them to see them pass**, then the unit tier.
- [ ] **Step 5: The production check** in `measure/task5c/`: resume the 50 mm state to 172.0 s with `--spread-limit 60` and
  compare with the snapshot rows: every shared column but the timers the same to the last bit (the limit cannot bind at a
  37 degree fold). Then one short run on the 50 mm state to 177.0 s with `--spread-limit 60`, recording the fold, `spread_mass_kg`,
  `spread_unplaced_kg`, books and wall time in `measure/task5c/check.md`. A fold above the limit there is a finding to report, not
  a reason to tune the pass without saying so.
- [ ] **Step 6: Finish**, with diffs keyed `*_T5C`.

### Task 3: The record

**Files:**
- Modify: `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` (§6.6 the pass and the hot landing; §10.1 the
  options; §10.2 the columns; §12 the fold's risk, with where the fold sits; §16 the decisions of 2026-11-06)
- Modify: `docs/superpowers/specs/2026-11-06-deposit-spreading-design.md` (status: built, what the checks showed)
- Modify: `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md` (an amendment section after fact 190, in its
  pattern: fact 191, the snapshot flights, from `measure/task5/snapshots/findings.md`; fact 192, the resume and spreading as built)

- [ ] **Step 1: Write them**, each amendment dated, saying who decided what, every number from the measurement files.
- [ ] **Step 2: Commit** in the outer repository, named paths only.
