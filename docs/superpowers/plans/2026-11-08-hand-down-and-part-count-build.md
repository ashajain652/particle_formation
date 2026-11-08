# The hand-down switched off, and the part count set by melting: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`)
> syntax for tracking. Asha approved the description and chose subagent execution on 2026-11-08.

**Goal:** Two options a skin flight can be compared on: `--hand-down off`, which keeps re-frozen metal in the skins, and
`--part-count-from melting`, which sets the next macro step's part count from the draws alone; and a guard that stops a
run whose part count has run away.

**Architecture:** Both are switches on `SkinBody` (`part_count_from`, `hand_down`), with class-level defaults so states
pickled before them resume, wired through `run` and `resume` in `reentry_model/cli.py` as `--spread-limit` is. Task 1 is
sub-plan 22's Task 5d (the part count and the guard), Task 2 its Task 5e (the hand-down off), Task 3 the short production
checks and the record. After the final review the controller launches the four comparison flights.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, scikit-fem, gmsh, pytest with xdist.

**Spec:** `docs/superpowers/specs/2026-11-08-hand-down-and-part-count-design.md` (the decision and why), and
`docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` sections 6.6 (the base: draw, absorb, hand down,
spread) and 7.1 (the parts and the hint).

## Global Constraints

- Paths contain a space (`University Documents /MIT`): quote every path. Interpreter:
  `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. Do not `cd` inside compound shell commands; use absolute
  paths or `git -C`.
- Outer repository: the worktree `.claude/worktrees/subplan20`, branch `step3-rigid-substrate-scheil`. Commit locally there
  after each task and **never push**. Stage named paths only, never `git add -A`. Never commit
  `docs/superpowers/plans/melt-spraying-subplans/23-notes.md`, `docs/superpowers/plans/2026-11-06-slip-band-spraying.md`
  (another session's files) or anything under `reentry_model_output/`. Never use bare `git stash`, never force anything.
- Never touch: the DRAMA databases, the wrapper scripts `sphere_reentry.py` and `sphere_sweep.py`, the package at the
  repository root, `prototype/proto3/`, `plan3`, `.claude/worktrees/we-critical-308`, the saved states, rows and frames
  under `code/reentry_model_output/` (read them, never write them), or the generated `22-skin-spraying-flight.md` by hand.
- Code is written and tested in the planning copy `prototype/work-2026-10-28-skin-spraying-plan/code/`. Every task that
  changes code ends by following **Finishing a task** below.
- Unit tier, from the planning copy's `code/`, on mains power only (check `pmset -g batt` first; if it says
  `Battery Power`, do not run it: report `DONE_WITH_CONCERNS` saying so):
  `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py -p no:cacheprovider -n 4`.
  Baseline 592 passed, 1 skipped, about 12 minutes. No failure is acceptable. Run it once per task, before committing.
- Books: in every skin test and run the mass is closed to 1e-12, the energy balance residual is below 1e-8, and the
  momentum books close. Bad settings exit 2 (`parser.error`); a model failure raises `RuntimeError` (exit 1).
- Without the new options every output of every run is as before to the last bit, column set included. The one
  exception is the guard of Task 1, which can only stop a run that would have run away.
- No new history columns. JSON records gain a key only when an option is not its default, as `spread_limit_deg` does.
- Repository style: comments and docstrings say what was decided, by whom and when, with measured numbers; test names are
  full sentences; commit messages are one long descriptive line in the existing pattern, ending with the trailer
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

### Finishing a task

1. **Before your first edit** to a file another task already changed, keep its state: for each file you will change, find
   the last task in `TASKS` in `make_plan22.py` that changes it, and if `states/task<that>/<rel>` does not exist, copy the
   file there unchanged. `make_plan22.py` refuses a replay that does not reproduce `code/`. (Task 1: the last task to change
   `skinbody.py`, `cli.py` and the test files is `5c`, and `states/task5c/` does not exist yet.)
2. Add your task to `TASKS`, and a task section to `part_22.md` immediately after the previous task's, before
   `## The outline`, in the pattern of Task 5c's: Files, Interfaces, what was decided, Measured, Rulings, Steps, with
   `{{KEY}}` markers for your diffs. Update the document's opening note ("Written at code level") and the outline's heading.
3. Run `"$PY" make_plan22.py` from the planning copy; it must print its check as passed. It regenerates the sub-plan.
4. Bring the build repository `prototype/work-2026-10-28-skin-spraying/code/` (git-ignored, its own git, branch
   `skin-spraying`, head cba0891 before these tasks) level with the planning copy's `code/` for your files, confirm `diff`
   finds no difference in them, and commit there with your task's message.
5. Commit in the outer repository: your `code/` files, `states/`, `part_22.md`, `make_plan22.py`, and the regenerated
   `docs/superpowers/plans/melt-spraying-subplans/22-skin-spraying-flight.md`.

## Review Focus

1. **States pickled before these tasks.** The snapshot states (commit 055c313) and the spreading states (509b835) have no
   `part_count_from` or `hand_down` on their `SkinBody`. Class-level defaults give them; each task's test deletes the
   instance attributes and advances.
2. **A spread limit carried by a resumed state.** `--hand-down off` on `resume` from a state saved with a limit (the
   spreading states carry 60) must exit 2 unless `--spread-limit off` is also given, and must work with it.
3. **The split test unchanged.** With `melting` a draw over the allowance must still split a part exactly as before;
   only the hint changes. Pinned in Task 1.
4. **The hint's reset.** `_ratio` is reset each macro step and maxed over kept parts only (a split part returns before
   it is maxed); the draw ratio for the hint must follow the same path, so a rejected part's draw does not count twice.
5. **Hand-down off and the `held` record.** With the hand-down off every handable kilogram is held; `held_mass_kg` and
   `held_depth_max_mm` then record all of it, and `jump_faults` sees no hand-down. The books must close with skins that
   only grow.

---

### Task 1: The part count set by melting, and the stall guard (sub-plan 22 Task 5d)

**Files:**
- Modify: `reentry_model/skinbody.py` (constant `PARTS_MAX`, `PART_COUNT_RULES`, `SkinBody(..., part_count_from="moves")`
  with a class default, `_move_ratio` or a sibling giving the draw ratio, the hint line in `_advance_part`, the guard in
  `advance`)
- Modify: `reentry_model/cli.py` (`--part-count-from moves|melting` on `run` and `resume`; name suffix; records; refusals)
- Test: `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: `SkinBody._move_ratio(book, m_min)`, `MOVE_MAX`, `MOVE_HINT_MARGIN`, `HAND_DOWN_SHARE`, `_hand_cap()`; the
  `--spread-limit` wiring in `cli.py` (`parse_spread_limit`, the name suffix at about line 225, the `run` record at about
  line 543, the `resume` handling at about lines 889-962) as the pattern to follow.
- Produces:
  ```python
  PART_COUNT_RULES = ("moves", "melting")
  PARTS_MAX = 1000        # parts: a macro step that would start with more stops the run (the controller, 2026-11-08)
  SkinBody(..., part_count_from="moves")     # class default part_count_from = "moves"
  ```
  `--part-count-from moves|melting` on `run` (default `moves`) and `resume` (default: the state's); exit 2 with `melting`
  unless `--surface-model skin --melt on` (run) or a saved flight of melting skins (resume). Name suffix on `run`:
  `_parts-melting`. Records: `run` JSON `settings.skin.part_count_from = "melting"` only when given; `resume` JSON
  `settings.part_count_from` only when the body's rule is `melting`, and `changed["--part-count-from"]` when given.

**What is decided:**
- With `melting`, the hint reads, for each kept part, `r_draw = max(book.drawn / (rho * area)) / (MOVE_MAX * mean edge)`,
  and `self._ratio = max(self._ratio, r_draw * dt_macro / dt)`; with `moves`, today's line unchanged. The split test
  (`ratio > 1.0 and depth < SPLIT_DEPTH_MAX`, `ratio` from `_move_ratio`) is the same under both rules.
- The guard, in `advance` after `n` is final (the hint and the `part_max` cap) and before the first part: if
  `n > PARTS_MAX`, raise `RuntimeError` saying the macro step from `t - dt` to `t` would run as `n` parts, more than
  `PARTS_MAX`, that the part count has run away, and that with `--part-count-from moves` a capped hand-down feeds it
  (sub-plan 22 Task 5d). It applies under both rules.
- Comment the change with the measured cause from the design (the sequence 837, 1 047, ... 9 765 on the 50 mm spreading
  flight, each the previous times 1.25 rounded up) and Asha's decision of 2026-11-08 to compare both approaches.

- [ ] **Step 1: Write the failing tests.** In `tests/test_reentry_model_skinbody.py`, on the melting test body
  (`melting_body`) or the hand-down cap's test set-up (`test_hand_down_takes_at_most_its_cap_a_cell_going_in_part`), with a
  skin given a held-back deposit deep enough that every macro step's hand-down is capped:
  `test_under_the_moves_rule_a_capped_hand_down_grows_the_next_steps_parts_by_a_quarter` (three macro steps: each
  `step_parts` is the previous times 1.25 rounded up; this pins today's defect);
  `test_under_the_melting_rule_a_capped_hand_down_leaves_the_part_count_alone` (the same set-up: `step_parts` stays at what
  the draws and `part_max` call for, and the books close);
  `test_a_draw_over_the_allowance_splits_a_part_under_either_rule` (equal `step_parts` and equal history rows under both);
  `test_the_moves_rule_by_default_changes_nothing_to_the_last_bit` (`part_count_from="moves"` against the argument
  omitted, every column);
  `test_a_macro_step_that_would_need_more_than_PARTS_MAX_parts_stops_the_run` (set `_hint = PARTS_MAX + 1`; `RuntimeError`
  naming the count; `PARTS_MAX` itself still runs);
  `test_a_body_saved_without_the_part_count_rule_flies_on` (delete the instance attribute, advance).
  In `tests/test_reentry_model_cli.py`: parsing, the suffix, the record, and the refusals (an unknown rule; `melting`
  without melting skins). In `tests/test_reentry_model_restart.py`:
  `test_resuming_with_the_melting_rule_sets_it_on_the_body_and_records_it`.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement** in `skinbody.py`, then `cli.py`.
- [ ] **Step 4: Run them to see them pass**, then the unit tier: at least 600 passed, 1 skipped expected.
- [ ] **Step 5: Finish** (Global Constraints, Finishing a task), with diffs keyed `*_T5D`.

### Task 2: The hand-down switched off (sub-plan 22 Task 5e)

**Files:**
- Modify: `reentry_model/skinbody.py` (`SkinBody(..., hand_down=True)` with a class default, `_move_base`, `_move_ratio`)
- Modify: `reentry_model/cli.py` (`--hand-down on|off` on `run` and `resume`; suffix; records; refusals)
- Test: `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: Task 1's wiring pattern for an option on `run` and `resume`; `skin_mod.hand_down`, `skin_mod.handable`.
- Produces: `SkinBody(..., hand_down=True)`, class default `hand_down = True`; `ValueError` if `hand_down` is False and
  `spread_limit` is not None. `--hand-down on|off` on `run` (default `on`) and `resume` (default: the state's); exit 2 with
  `off` unless melting skins, and with `off` and a spread limit, whether given or carried by the resumed state, with a
  message saying spreading moves only metal being handed down (Asha, 2026-11-08) and to give `--spread-limit off`. Name
  suffix on `run`: `_handdown-off`. Records: `run` JSON `settings.skin.hand_down = "off"` only when off; `resume` JSON
  `settings.hand_down = "off"` only when the body's is off, and `changed["--hand-down"]` when given.

**What is decided:**
- With `hand_down` False, `_move_base` does not call `skin_mod.hand_down`: `handed` and `handed_H` are zeros, `coldest`
  and `warmest` nan, `uncovered` nan; `held = able`, so `held_mass` and `held_depth_max` record every handable kilogram;
  the move, the reconciliation and the fault count read the zero hand-down. `_move_ratio` uses a zero hand-down too.
  Draws and absorption are unchanged, so the base moves only inward (and by the existing neutral smoothing).
- Comment it with Asha's question and decision of 2026-11-08: re-frozen metal can heat again and flow on or be sprayed,
  so it may stay in the skins; and spreading is off with it because it moves only handed metal.

- [ ] **Step 1: Write the failing tests.** In `tests/test_reentry_model_skinbody.py`, on a body that hands metal down
  with the default (the set-up of Task 1, or `melting_body` run until `handed_mass` grows):
  `test_with_the_hand_down_off_no_metal_crosses_the_base_outward` (`handed_mass` stays 0 where the default's grows; no
  corner moves outward along its normal beyond the neutral smoothing's share, compared with a run whose smoothing is the
  same);
  `test_with_the_hand_down_off_the_rigid_metal_stays_in_its_skin_and_is_counted_as_held` (the skins' mass above their
  minimum is at least the default's, and `held_mass` equals what `skin_mod.handable` finds);
  `test_with_the_hand_down_off_the_books_close` (mass 1e-12, energy 1e-8, momentum);
  `test_the_hand_down_on_by_default_changes_nothing_to_the_last_bit`;
  `test_the_hand_down_off_refuses_a_spread_limit` (`ValueError`);
  `test_a_body_saved_without_the_hand_down_switch_flies_on`.
  In `tests/test_reentry_model_cli.py`: parsing, the suffix, the record, the refusals (`off` without melting skins; `off`
  with `--spread-limit 60`). In `tests/test_reentry_model_restart.py`:
  `test_resuming_with_the_hand_down_off_from_a_state_with_a_spread_limit_needs_the_limit_off` (exit 2 without
  `--spread-limit off`, flies with it, recorded).
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement** in `skinbody.py`, then `cli.py`.
- [ ] **Step 4: Run them to see them pass**, then the unit tier.
- [ ] **Step 5: Finish**, with diffs keyed `*_T5E`.

### Task 3: The production checks and the record

**Files:**
- Create (git-ignored): `prototype/work-2026-10-28-skin-spraying/measure/task5d/check.md`, `.../measure/task5e/check.md`
- Modify: `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` (§6.6 the hand-down switch; §7.1 the two
  part-count rules and the guard; §10.1 the options; §16 Asha's decisions of 2026-11-08)
- Modify: `docs/superpowers/specs/2026-11-08-hand-down-and-part-count-design.md` (status: built, what the checks showed)
- Modify: `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md` (an amendment section after fact 192, in its
  pattern: fact 193, the spreading flights and their stall; fact 194, the two options and the guard as built and checked)

**The checks**, on mains power, through `prototype/work-2026-10-28-skin-spraying-plan/flights/launch.sh` from the
planning copy's code, run side by side, output into `code/reentry_model_output/task5d_check` and `task5e_check`:
- **5d:** resume `code/reentry_model_output/task5_spreading/d050__spread60.state.pkl` (178.5 s, limit 60 carried) with
  `--part-count-from melting --t-max 179.7`. Compare `step_parts`, `held_mass_kg`, `surface_fold_deg`, `zone_ratio_min` and
  `rebuilds` row by row with `task5_spreading/d050__spread60.rows.csv`, where the count grew from 24 at 179.3 s to 9 765 at
  179.66 s. The fix works if the count stays near what the draws call for; if it still grows, say so and what the rows show,
  without tuning anything.
- **5e:** resume `code/reentry_model_output/task5_snapshots/d050__physics__skin.state.pkl` (170.5 s) with
  `--hand-down off --t-max 177.0`, the window of Task 5c's check (`measure/task5c/check.md`: the fold 116.75 degrees without
  spreading, 64.86 with it). Record the base fold, `handed_mass_kg` (unchanged from the saved row), the deepest skin
  (`skin_thickness_max_mm`, `flowing_depth_max_mm`), `held_mass_kg`, rebuilds, the lowest radius ratio, `step_parts`.
- Each `check.md` gives the command, the wall time, what was compared and the books from the run's JSON (mass, energy
  residual against 1e-8, momentum).

- [ ] **Step 1: Run the checks** and write both `check.md`.
- [ ] **Step 2: Write the record**, each amendment dated, saying who decided what (Asha: test both and compare, spreading
  on for the melting rule and off with the hand-down off, built with subagents; the controller: the defaults stay today's,
  the guard at 1 000 parts), every number from the measurement files or the design.
- [ ] **Step 3: Commit** in the outer repository, named paths only.

## After the final review (the controller, not a task)

Launch the four comparison flights of the design on mains power, detached, side by side, from the planning copy's code
through `flights/launch.sh`, output into `code/reentry_model_output/task5_compare`:
`d050__handdown-off` and `d100__handdown-off` (`--hand-down off`), `d050__parts-melting` and `d100__parts-melting`
(`--spread-limit 60 --part-count-from melting`), from the snapshot states at 170.5 s and 22.0 s, `--t-max 190` and `45`,
`--save-state-at 178.5` and `32.5`, `--frames-every 1`. A watcher on their process ids.
