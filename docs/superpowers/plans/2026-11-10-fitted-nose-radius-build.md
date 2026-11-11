# The windward loads from a sphere fitted to the skins' tops: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`)
> syntax for tracking. Asha decided the fix and its three choices on 2026-11-10 ("let's stop the flights and implement
> fixes"); the build follows the subagent pattern she chose for Tasks 5d-5f on 2026-11-08.

**Goal:** The skins' windward loads (stagnation heating and surface flow) read the radius of a sphere fitted to the tops of
the skin columns facing the flow, reported every macro step, on by default.

**Architecture:** One task, sub-plan 22's Task 5g. `SkinBody` fits the sphere in `_refresh_geometry` and returns it from
`nose_radius()` when its `size_feedback` is `"current"`; `cli.py` wires `--size-feedback` for skin runs on `run` and
`resume`; the skins' history gains five columns.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, scikit-fem, gmsh, pytest with xdist.

**Spec:** `docs/superpowers/specs/2026-11-10-fitted-nose-radius-design.md` (the decision, the controller's calls and the
measured radii). Pattern for an option on `run` and `resume`: Tasks 5d-5f in
`docs/superpowers/plans/2026-11-08-hand-down-and-part-count-build.md`.

## Global Constraints

- Paths contain a space (`University Documents /MIT`): quote every path. Interpreter:
  `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. Do not `cd` inside compound shell commands; use absolute
  paths or `git -C`.
- Outer repository: the worktree `.claude/worktrees/subplan20`, branch `step3-rigid-substrate-scheil`. Commit locally there
  after the task and **never push**. Stage named paths only, never `git add -A`. Never commit
  `docs/superpowers/plans/melt-spraying-subplans/23-notes.md`, `docs/superpowers/plans/2026-11-06-slip-band-spraying.md`,
  `docs/superpowers/specs/2026-11-08-refrozen-metal-in-skins-design.md` (another session's files) or anything under
  `reentry_model_output/` or `prototype/work-2026-10-28-skin-spraying/measure/`. Never use bare `git stash`, never force
  anything.
- Never touch: the DRAMA databases, the wrapper scripts `sphere_reentry.py` and `sphere_sweep.py`, the package at the
  repository root, `prototype/proto3/`, `plan3`, `.claude/worktrees/we-critical-308`, the saved states, rows and frames
  under `code/reentry_model_output/` (read them, never write them), or the generated `22-skin-spraying-flight.md` by hand.
- **A flight is running from the planning copy's `code/`** (the 100 mm comparison flight, to 45 s). Edit files there as
  usual (its modules are already loaded), but never move, rename or delete anything under `code/`, and never write under
  `code/reentry_model_output/`.
- Code is written and tested in the planning copy `prototype/work-2026-10-28-skin-spraying-plan/code/`, and the task ends
  by following **Finishing a task** below.
- Unit tier, from the planning copy's `code/`, on mains power only (check `pmset -g batt` first; if it says
  `Battery Power`, do not run it: report `DONE_WITH_CONCERNS` saying so):
  `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py -p no:cacheprovider -n 4`.
  Baseline 639 passed, 1 skipped, about 15 minutes (longer now: two long runs share the laptop). No failure is acceptable.
  Run it once, before committing.
- Books: in every skin test and run the mass is closed to 1e-12, the energy balance residual is below 1e-8, and the
  momentum books close. Bad settings exit 2 (`parser.error`); a model failure raises `RuntimeError` (exit 1).
- With `--size-feedback initial` every output of every run is as before to the last bit, apart from the five new columns.
- Repository style: comments and docstrings say what was decided, by whom and when, with measured numbers; test names are
  full sentences; commit messages are one long descriptive line in the existing pattern, ending with the trailer
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

### Finishing a task

1. **Before your first edit** to a file another task already changed, keep its state: for each file you will change, find
   the last task in `TASKS` in `make_plan22.py` that changes it, and if `states/task<that>/<rel>` does not exist, copy the
   file there unchanged. `make_plan22.py` refuses a replay that does not reproduce `code/`.
2. Add your task to `TASKS`, and a task section to `part_22.md` immediately after Task 5f's, before `## The outline`, in
   the pattern of Task 5f's: Files, Interfaces, what was decided, Measured, Rulings, Steps, with `{{KEY}}` markers for your
   diffs. Update the document's opening note ("Written at code level") and the outline's heading.
3. Run `"$PY" make_plan22.py` from the planning copy; it must print its check as passed. It regenerates the sub-plan.
4. Bring the build repository `prototype/work-2026-10-28-skin-spraying/code/` (git-ignored, its own git, branch
   `skin-spraying`, head c33b4a2 before this task) level with the planning copy's `code/` for your files, confirm `diff`
   finds no difference in them, and commit there with your task's message.
5. Commit in the outer repository: your `code/` files, `states/`, `part_22.md`, `make_plan22.py`, the regenerated
   `docs/superpowers/plans/melt-spraying-subplans/22-skin-spraying-flight.md`, this plan and its spec
   `docs/superpowers/specs/2026-11-10-fitted-nose-radius-design.md`.

## Review Focus

1. **States pickled before this task** (the snapshot states of 055c313, the spreading states of 509b835, the comparison
   flights' states of 2026-11-08) have no `size_feedback` and no fit on their `SkinBody`. A class-level default
   `size_feedback = "initial"` must resume them exactly as flown; `--size-feedback current` on `resume` must fit at the
   first refresh.
2. **The loads actually read it**: both the heating (`CoupledRun.loads_at`) and the surface flow (`SkinBody` at about line
   790) call `nose_radius()`; a test must show each receiving the fitted radius, not only that `nose_radius()` returns it.
3. **The bound applies on the flat side** (Asha's option B, 2026-11-10, replacing "reported, never applied"): a convex
   cap flatter than the bound, a concave cap and a cap flat to rounding must give the bound as the used radius; a convex
   cap no flatter than the bound gives its fitted radius; there is no lower bound.
4. **Degenerate caps**: fewer than four cap corners, and a concave cap, must give the documented radius and flags, not a
   crash or a negative radius.
5. **The leeward pool**: the windward transverse radius must ignore corners more than 90 degrees from the direction of
   motion, so a pool standing off the back does not move the bound.

---

### Task 1: The windward loads from a sphere fitted to the skins' tops (sub-plan 22 Task 5g)

**Files:**
- Modify: `reentry_model/skinbody.py` (`SkinBody(..., size_feedback="current")` with class default `"initial"`; the fit in
  `_refresh_geometry`; `nose_radius()`; the five columns in `skin_stats()`)
- Modify: `reentry_model/coupled.py` (`SKINBODY_COLUMNS` gains the five columns)
- Modify: `reentry_model/cli.py` (`--size-feedback` accepted with `--surface-model skin` on `run`, default `current`; on
  `resume`, default the state's; records)
- Test: `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`,
  and `tests/test_reentry_model_coupled.py` if the column list is pinned there

**Interfaces:**
- Consumes: `body.fit_sphere(points) -> (centre, radius)`, `body.NOSE_CAP_ANGLE` (30.0), `body.NOSE_CAP_FACTOR` (1.67),
  `body.SIZE_FEEDBACK_NAMES`; `SkinBody.surface` (a `columns.CornerSurface` of the column tops, `.points`, `.normals`),
  `SkinBody.theta` (angle of each top corner's normal to the direction of motion), `SkinBody.v_hat`, `self.mesh.points`;
  the `--hand-down` wiring in `cli.py` (`hand_down_of`, the `resume` handling) as the pattern.
- Produces:
  ```python
  SkinBody(..., size_feedback="current")      # class default size_feedback = "initial" (states saved before Task 5g)
  SkinBody.nose_radius() -> float             # the fitted radius when "current" (see What is decided), else radius_m
  SkinBody.nose_fit -> NoseFit                # fitted [m] (nan if none), bound [m], points (int), concave (bool), used [m]
  ```
  History columns, after the existing skin columns: `nose_radius_mm`, `fitted_nose_radius_mm`, `nose_bound_mm`,
  `nose_cap_points`, `nose_cap_concave`, written with either setting. `run`: `--size-feedback initial|current` with
  `--surface-model skin`, default `current`; JSON `settings.skin.size_feedback` always (the default changed, so record it).
  `resume`: `--size-feedback`, default the saved body's (`initial` for a state without it); JSON `settings.size_feedback`
  when `current`, and `changed["--size-feedback"]` when given. The refusal of `--size-feedback` with
  `--surface-model skin` at about line 746 goes; the other element-model refusals stay.

**What is decided** (the spec's decision and calls; the definition the signature does not determine):
- The cap: top corners with `theta <= radians(NOSE_CAP_ANGLE)`. Fewer than 4: no fit (fitted nan, used = `radius_m`).
- The fit: `body.fit_sphere(self.surface.points[cap])`. Not finite or not positive: no fit, as above.
- Concave: the fitted centre lies on the flow side of the cap, `(centre - cap_points.mean(axis=0)) @ v_hat > 0`; then
  `concave = True` and used = bound.
- The bound: `NOSE_CAP_FACTOR * R_tw`, `R_tw` the largest distance from the axis through `self.mesh.points.mean(axis=0)`
  along `v_hat` of the top corners with `theta <= pi/2`. Reported only.
- Used radius when `"current"` and a convex fit exists: the fitted radius, with no bound applied. When `"initial"`:
  `radius_m` always (the fit is still computed and reported).
- Comment the change with the spec's measured radii and Asha's decisions of 2026-11-10 (fit the gas-molten interface;
  the bound a comparison, not a limit; on by default).

- [ ] **Step 1: Write the failing tests.** In `tests/test_reentry_model_skinbody.py`, on the existing skin test bodies
  (move the column tops, or give the skins a thickness, to shape the cap):
  `test_the_fitted_nose_radius_of_an_untouched_sphere_is_its_own_radius` (within 2 % of `radius_m`);
  `test_a_flattened_windward_cap_gives_the_radius_of_the_flattened_cap` (tops on a sphere of twice the radius over the
  cap: fitted within 5 % of it, and `nose_radius()` returns it under `"current"`);
  `test_the_bound_is_reported_and_never_applied` (a cap flatter than the bound: used radius above `nose_bound`);
  `test_the_heating_and_the_surface_flow_both_receive_the_fitted_radius` (record the radius each is called with);
  `test_a_cap_with_fewer_than_four_corners_falls_back_to_the_original_radius`;
  `test_a_concave_cap_uses_the_flat_face_bound_and_says_so`;
  `test_a_pool_behind_the_body_does_not_move_the_bound` (thicken leeward columns only: `nose_bound` unchanged);
  `test_size_feedback_initial_changes_nothing_to_the_last_bit` (every pre-existing column equal to a body built without the
  argument and with `"initial"`; the five new columns present);
  `test_a_body_saved_without_size_feedback_flies_on_as_initial` (delete the instance attributes, advance).
  In `tests/test_reentry_model_cli.py`: `--size-feedback` accepted with skins, default `current`, recorded; still refused
  where it was for other reasons, if any. In `tests/test_reentry_model_restart.py`:
  `test_resuming_with_size_feedback_current_fits_the_cap_and_records_it`.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement** in `skinbody.py`, then `coupled.py`, then `cli.py`.
- [ ] **Step 4: Run them to see them pass**, then the unit tier: at least 648 passed, 1 skipped expected.
- [ ] **Step 5: Finish** (Global Constraints, Finishing a task), with diffs keyed `*_T5G`.

### Task 2: The new defaults, and the Rayleigh-Taylor droplet without the body-radius cap (sub-plan 22 Task 5h)

Decided by Asha on 2026-11-10 (spec, "Further decisions", items 1 and 2; the controller's call on options).

**Files:**
- Modify: `reentry_model/cli.py` (defaults of `--hand-down` and `--part-count-from` on `run` for melting skins; the
  spread-limit refusal; `--spray-gates 2026-11-10|before` on `run` and `resume`; name suffixes; records)
- Modify: `reentry_model/spray.py` (the Rayleigh-Taylor branch's radius on the crest-limit path), `reentry_model/skinbody.py`
  (carry the spray-gates setting to the spray plan, class default `"before"` for saved states)
- Test: `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`, `tests/test_reentry_model_spray.py`,
  `tests/test_reentry_model_skinbody.py`

**What is decided:**
- `run` with `--surface-model skin` and melting: `--hand-down` defaults to `off` and `--part-count-from` to `melting`;
  with skins that do not melt the defaults stay `on` and `moves` (the options are refused there with `off`/`melting`, as
  today). `--spread-limit` on such a run without `--hand-down on` exits 2 saying spreading moves only metal being handed
  down (Asha, 2026-11-08) and to give `--hand-down on`. Name suffixes appear only where a value differs from the run's
  default (`_handdown-on`, `_parts-moves` for melting skins); JSON `settings.skin.hand_down` and `.part_count_from` are
  always recorded. `resume` keeps the saved run's values unless given (unchanged).
- `--spray-gates 2026-11-10|before`, default `2026-11-10` on `run`, the saved run's on `resume` (`before` for states
  without it). Under `2026-11-10`, on the crest-limit path (`cap_depth` given: the skins) a Rayleigh-Taylor droplet's
  radius is `min(lambda_RT, crest)`, with no `0.25 * radius` term; every other branch and the element model's path
  (`cap_depth` None) are unchanged. Under `before`, today's line to the last bit. Task 3 adds the thin branch's gate under
  the same setting.
- Comment with Asha's decision (the Rayleigh-Taylor droplet is its wavelength, Girin and Kopyt's front-surface r_d; the
  quarter-radius cap is not from the mode's physics and set the 15.8 mm droplets of 2026-11-10 through a transverse
  radius inflated to 57-64 mm by misplaced column tops; the crest limit, a supply limit, stays).

- [ ] **Step 1: Write the failing tests**: the defaults for melting and non-melting skins; the spread-limit refusal and
  its `--hand-down on` escape; the suffixes and records; `test_a_rayleigh_taylor_droplet_is_not_capped_by_the_body_radius`
  (a deep pool, a small `radius`: r equals the crest-limited wavelength, above `0.25 * radius`);
  `test_the_crest_limit_still_binds_a_rayleigh_taylor_droplet_on_a_shallow_pool`;
  `test_spray_gates_before_reproduces_the_old_radius_to_the_last_bit`; restart: a state without the setting resumes as
  `before`, and `--spray-gates 2026-11-10` on resume is recorded.
- [ ] **Step 2: Run them to see them fail.** **Step 3: Implement.** **Step 4: Pass, then the unit tier.**
- [ ] **Step 5: Finish**, with diffs keyed `*_T5H`.

### Task 3: The thin branch's Knudsen gate, and spraying in the slip band (sub-plan 22 Task 5i)

Decided by Asha on 2026-11-10 (spec, "Further decisions", item 3, as amended the same day): **follow the slip-band plan of
2026-11-06**, `docs/superpowers/plans/2026-11-06-slip-band-spraying.md` (written by a side session; read it in full, its
Tasks 1-3 and Known limits; it is that session's file, never edit or commit it), with Asha's six decisions of 2026-11-06
recorded at its top, and with three changes of 2026-11-10:

1. **The band needs both Knudsen numbers below 0.1**: a patch is in the band where the body is on the merged branch with
   a body Knudsen number below `KN_BODY_SLIP = 0.1` (the 2026-11-06 plan's `in_slip_band`) **and** the patch's own wall
   Knudsen number (`flow.kn_local`) is below 0.1. Outside it, no band depth.
2. **The thin branch's continuum gate**: under `--spray-gates 2026-11-10` (Task 2's setting), the thin branch, and the
   rarefied branch with it, release only where a conjugate depth exists for spraying: Girin's closure below Knudsen 0.01
   (`melt_layer`'s finite `delta_m`), or the band's momentum-balance depth where `--slip-band on` puts a patch in the band
   (`spray_layer`'s finite depth). Elsewhere a wet patch takes no shear mode. The Rayleigh-Taylor mode is not gated (a
   melt-side mode, driven by the deceleration). Under `before`, as today. This applies wherever `SprayModel.evaluate`
   runs, the element model's spraying included.
3. **The slurry in the band**: as the 2026-11-06 plan's decision 4, `--slurry-spray on` lets the slurry spray in the
   band through `slurry_modes`; with it off the release stays capped at the liquid.

What the 2026-11-06 plan says that is out of date, and is replaced here: it calls itself "Task 5c" (taken by spreading;
this is Task 5i, diffs keyed `*_T5I`); its baselines (541 passed) are now 639+ (see Global Constraints); its "Commit
only when Asha asks" is replaced by this build's local commits; its Task 4 measurement (the 74 km short case) is replaced
by the test flights after the final review. Its option name `--slip-band on|off` (default `off`; on `resume` the saved
run's, `off` for states without it) and `KN_BODY_SLIP` are kept. The momentum-balance depth uses `nose_radius()`, so it
reads Task 1's fitted radius.

**Files** (the 2026-11-06 plan's Tasks 1-3, plus the gate and the wall-Knudsen condition):
- Modify: `reentry_model/surface_flow.py`, `reentry_model/spray.py`, `reentry_model/skinbody.py`, `reentry_model/body.py`,
  `reentry_model/cli.py`
- Test: `tests/test_reentry_model_surface_flow.py`, `tests/test_reentry_model_spray.py`,
  `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_melting.py`, `tests/test_reentry_model_cli.py`,
  `tests/test_reentry_model_restart.py`

- [ ] **Step 1: Write the failing tests**: the 2026-11-06 plan's tests of its Tasks 1-3 (adapted to the code as it now
  is), plus `test_a_patch_whose_wall_knudsen_number_is_above_a_tenth_is_outside_the_band` (body Kn 0.04, wall Kn 0.12:
  no band depth, no shear spray); `test_under_the_new_gates_the_thin_branch_does_not_spray_outside_both_depths` (merged
  branch, band off: thin releases nothing, Rayleigh-Taylor still can); `test_with_the_band_on_the_thin_and_thick_branches_spray_in_it`;
  `test_spray_gates_before_and_band_off_change_nothing_to_the_last_bit`; cli and restart wiring and records.
- [ ] **Step 2: Run them to see them fail.** **Step 3: Implement.** **Step 4: Pass, then the unit tier.**
- [ ] **Step 5: Finish**, with diffs keyed `*_T5I`.

### Task 4: The leeward test by the body's centroid (sub-plan 22 Task 5j)

Decided by Asha on 2026-11-10 (spec, "Further decisions", item 4, and her answers to questions 2 and 3; the controller's
calls on the fit region and on which consumers read the test).

**Files:**
- Modify: `reentry_model/skinbody.py` (`SkinBody.mass_centroid()`; the effective angle; `nose_radius` fit points;
  `--leeward-test` setting with class default `"facing"`), `reentry_model/coupled.py` (the heating reads the effective
  angle), `reentry_model/cli.py` (`--leeward-test facing|centroid` on `run`, default `centroid`, and `resume`, default the
  saved run's)
- Test: `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_coupled.py`, `tests/test_reentry_model_cli.py`,
  `tests/test_reentry_model_restart.py`

**What is decided:**
- `SkinBody.mass_centroid()`: the mass-weighted centre of the 3D model's metal (its elements' masses at their centroids)
  and of every skin cell (its mass at the cell's centre along its column, base point plus normal times the depth from the
  base to the cell's mid-height), refreshed with the geometry.
- With `centroid`: a top corner is windward only if `(p - c) @ v_hat > 0` (ahead of the plane through the centroid at
  right angles to the flight direction) and `theta <= pi/2`. The effective angle is `theta` where windward and
  `max(theta, pi/2 + 1e-6)` where not; the heating (`CoupledRun.loads_at`), the surface flow and the spray plan
  (`girin_plan`, so the shear modes and, through `windward` and `w_n`, the Rayleigh-Taylor mode) all read it; `theta`
  itself, the frames and the release angle are unchanged. With `facing`, every consumer reads `theta` as today, bit for
  bit.
- The windward sphere fit (Task 1) takes the corners windward by this test with `theta <= NOSE_CAP_ANGLE`; the windward
  transverse radius for the bound takes the windward corners. Under `facing`, Task 1's selection unchanged.
- History: `centroid_offset_mm` (the centroid's position along the flight direction relative to the 3D mesh's mean point,
  positive forward) and `leeward_corners` (corners facing the flow but behind the plane), with either setting.
- Comment with the measured case: releases from corners up to 82 mm behind the original centre facing 12-20 degrees off
  the flow (median) on the 50 mm sphere at 186.5-187.3 s; the centroid 30 mm behind the original centre by 187.0 s.

- [ ] **Step 1: Write the failing tests**: `test_the_mass_centroid_of_an_untouched_sphere_is_its_centre`;
  `test_a_pool_at_the_back_moves_the_centroid_backward_by_its_mass_moment`;
  `test_a_corner_behind_the_centroid_that_faces_the_flow_gets_no_shear_no_spray_and_no_rayleigh_taylor`;
  `test_a_corner_ahead_of_the_centroid_that_faces_away_gets_no_shear_or_spray`;
  `test_the_heating_reads_the_effective_angle`; `test_the_leeward_test_facing_changes_nothing_to_the_last_bit`;
  `test_the_fit_takes_only_windward_corners_within_the_cap`; cli and restart wiring and records.
- [ ] **Step 2: Run them to see them fail.** **Step 3: Implement.** **Step 4: Pass, then the unit tier.**
- [ ] **Step 5: Finish**, with diffs keyed `*_T5J`.

## After the final review (the controller, not a task)

Re-fly the 100 mm sphere from its snapshot at 22.0 s to 45 s with all the fixes, through `flights/launch.sh` on mains
power: `--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 --leeward-test centroid
--sideways-heat height`, once with `--slip-band off` and once with `--slip-band on` (spec, the controller's call on
question 6), frames every macro step, a state saved at 32.5 s.
