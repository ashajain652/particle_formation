# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1a, 1b, 1c, 2, 3 and 3b (2026-10-28/29),
> 4 and 4b (2026-10-31), 5a (2026-11-01), 5b and 5c (2026-11-06), and 5d, 5e and 5f (2026-11-08).**
> Task 1 is split in three, in the order the amendment of 2026-10-28 sets: 1a brings the element model's two spraying
> amendments of 2026-10-22 across, 1b adds the droplet's cap per crest to the shared spray physics, and 1c sprays from
> each column's top cells, by Asha's decisions of 2026-10-29 on how the skins' depths are read where liquid lies on slurry
> and on the slurry cascade. Task 2 makes the droplet records from the skins. Task 3 retires stage 1's one-dimensional
> skins, makes the layered skins the default surface model and puts the spraying skin flight on the command line; Task
> 3b, by Asha's decisions on the two questions it left, makes the flow the default of a melting skin run and latches
> the skins' continuum step at their first liquid. Task 4 writes a spraying skin flight's outputs -- the summary, the
> rest of spec 10.2 and Register's release angle on both models' records (Asha, 2026-10-29) -- ends the skins' loop
> without a vanishing last step, writes their records to disk as they come and settles the FEniCSx tier's failure;
> 4b, by Asha's decision on its ruling 6, makes the resolved sheared depth the straight-line equivalent, with no cut-off.
> Task 5a, for Task 5's convergence repeats, saves a flight's whole state on its way and resumes from it to the last bit; 5b puts the resume on the command line.
> Task 5c spreads the metal being handed down to lower neighbouring corners where the base would fold past a limit, off by default until the flights have tested it.
> Tasks 5–6 are the outline's, below, and the inputs it listed are carried at the end.

**Goal:** Girin's spraying acts on each column's top cells, with the film thickness b the whole flowing depth, by the
rules of record (the element model's amendments of 2026-10-22) and a droplet no larger than one crest of its own wave;
the old film, deep account and their rules are removed; and the model flies the 100 mm and 50 mm cases to the spec's
acceptance.

**Architecture:** Task 1a changes only the shared spray physics and the element model's callers, exactly as branch
`we-critical-308` changed them; Task 1b gives the shared spray physics a second size cap, used only when its caller
passes the depth the release can take, so the element model is unchanged. Task 1c evaluates Girin per column once at
the start of each macro step, from the surface flow the skins' flow already evaluates there and the skins' depths then,
and hands the rate per column to sub-plan 21 Task 7's release hook, which takes it from the top cells in every flow
sub-step. Task 2 books what the release took by the cells it came from, combines it per column, per mode and per macro
step into the plan's record, and writes it as source-table rows through the shared physics' own row writer, so that one
coupled writer serves both surface models. Task 3 removes the element model's dead branches for stage 1's skins, so
that the element model is unchanged in everything it still writes, and changes only the command line otherwise: the
skins become the default surface model, and with their flow on the command line gives them the shared spray model and
writes their droplets with the element model's writer.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, scikit-fem for the 3D model, gmsh for its rebuilds.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` — §8 (spraying; amended 2026-10-24 and
2026-10-28: the rules of record and the cap per crest), §6.4 (the droplets' momentum), §6.7 (the depths read from the
skin), §9 (what is retired), §10.2 (outputs), §11 (build stage 4; level 3: convergence and the acceptance criterion; the
comparison with the element model), §12 (the deposit's folding, a watched risk for the flight), §14 (the assumptions
the thesis states), §16 (the decisions of 2026-10-24 and 2026-10-28).

**Builds on:** sub-plan 21 — the flowing layer, the droplet momentum sink at the top cells (Task 7), the resolved
sheared depth, the base where runoff lands (Task 9). Build stage 3 is complete: levels 1 and 2 pass with the flow on
(facts 151–160).

**How this plan's code is checked.** Every code and diff block below was written and tested in a planning copy,
`prototype/work-2026-10-28-skin-spraying-plan/code/`, made from sub-plan 21's tested copy,
`prototype/work-2026-10-18-column-flow-plan/code/`, as committed at 0f0425e (identical to the build's 602b262). This
document is generated from it by `prototype/work-2026-10-28-skin-spraying-plan/make_plan22.py`, which replays the tasks
written so far on a copy of the base — every diff through `patch`, each from the state the tasks before it left (kept,
for files a later task changes again, in `states/task<N>/`) — and checks that the result is the tested copy byte for
byte. It replaces the outline of the same name, whose header, amendments, tasks and inputs it carries. Never edit this
document by hand: change `part_22.md` or the tested copy, and regenerate.

## Global Constraints

- Work in the build repository `prototype/work-2026-10-28-skin-spraying/code/`: a git repository of its own, ignored by
  the outer one, branch `skin-spraying`, cloned from the column-flow build at 602b262 (sub-plan 21 Task 9). Commit there
  after every task. Never edit `prototype/proto3/`, `plan3`, the package `reentry_model/` at the repo root, the wrapper
  scripts (`sphere_reentry.py`, `sphere_sweep.py`), the DRAMA databases, the specs or the generated plans. Step 3 code
  stays in the prototype until the skyscrapers architecture is verified (Asha, 2026-10-12).
- The worktree `.claude/worktrees/we-critical-308` (branch `we-critical-308`, head a54520c) is read-only: its
  amendments are brought across from its `prototype/work-2026-10-22-we-critical/amendment/` and
  `prototype/work-2026-10-22-thin-gate/amendment/` diffs, and its tested copies are copied, never run in place.
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the work copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py -p no:cacheprovider`.
  Baseline at 602b262: 517 passed, 1 skipped (about 18–34 minutes). Record the counts after every task; no failure is
  acceptable.
- FEniCSx tier: `FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`. Baseline 19 passed.
- `--surface-model elements` is, from Task 1a on, the element model of branch `we-critical-308` (spec §16,
  2026-10-28): Task 3's bit-for-bit check of the element model is against that branch's flights, not `main` at 612b916.
  Task 1b leaves it unchanged to the last bit.
- Books: in every skin test and run, the mass is closed to 1e-12, the energy balance residual is below 1e-8 and the
  momentum books close as sub-plan 21's do.
- Exit codes: a Newton or coupled solve that does not converge raises `RuntimeError` (exit code 1 at the command line);
  bad settings exit 2.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).
- Quote every path: the repository's path contains a space (`University Documents /MIT`).

## Review Focus

1. **The element model under the cap** — without a cap depth `SprayModel.evaluate` takes the film-mass cap by the same
   arithmetic as before, so `--surface-model elements` is unchanged to the last bit by Task 1b: on the CLI test's short
   melting case every one of the 115 history columns, the 25 266 source rows and the results (fact 170). Pinned by
   `test_without_a_cap_depth_the_film_mass_caps_the_droplet_as_before`.
2. **Two bases** — the amendments were made on `main`'s element model (612b916), this line's differs outside them:
   `main`'s frame export of 2026-10-08 (`p_w_stag_step_Pa`, `v_hat_body`) is not here, sub-plan 18's one-dimensional
   skins and resolved depths and sub-plan 21's `dp_ds` are. Three hunks were taken by hand (Task 1a), and on the short
   melting case this line's element model sprays as that branch's does to the last bit in every column the two share
   (fact 169). The rarefied and regime-2 branches are reached only by the unit tests.
3. **A depth of nothing** — with a cap depth of zero one crest gives nothing: no droplet forms and nothing leaves,
   though the patch keeps its branch (`test_a_cap_depth_of_nothing_releases_nothing`). Task 1c must pass the depth at
   the macro step's start, so a column whose liquid forms within the step waits for the next.
4. **Rayleigh–Taylor droplets set by the crest** — their radius is their wavelength (centimetres), so the cap sets them
   on any layer: 3–7 mm on 0.3–3.3 mm of liquid at 103 m/s², against about 2 mm under the film cap on the element
   model's flight; where the Rayleigh–Taylor criterion reads the cap's own depth h they stay under 2.1 h, inside the
   10 mm top bin while h is under 4.7 mm (fact 170). Task 2's histograms and Task 5's flights report them.
5. **The viscous film's floor on the skins** — the flowing depth (Asha, 2026-10-28), which counts slurry as part of the
   film with the liquid's viscosity; fact 165's open assumption (the effective floor may lie inside the slurry, the
   damping larger) carries over (Task 1c). Pinned, with the rest of the depths' reading (option B, 2026-10-29), by
   `test_girin_gives_each_prescribed_column_its_mode_from_the_skin_s_depths_once_for_the_macro_step`.
6. **Girin once a macro step** (Task 1c's design) — the rate, the mode and the radius per column are fixed for the
   macro step from its start, as the cap is; a column that becomes wet or changes regime within it waits for the next,
   12.5 ms at the flow-on runs' continuum step and up to 0.5 s at the provisional default. Task 5's halving of the macro
   step measures it. On the skin body's test the flag-off step sprayed only the 144 columns liquid at its start (Task
   1c). Pinned by `test_a_skin_body_sprays_by_girin_from_its_top_cells_once_a_macro_step_and_keeps_its_books`.
7. **The cascade near the liquidus** (Task 1c, option i) — the exposed slurry sprays only within about 0.6–0.8 K of the
   liquidus, where the bridge's viscosity (no data) is within about ten times the liquid's, but on a heated surface that
   is where the slurry under the liquid, and a top just entering the bridge, are: with the flag on the skin body's test
   sprayed 4.98 g of slurry against 0.79 g of liquid (1.54 g of liquid with it off). Off by default. Pinned by
   `test_the_slurry_cascade_sprays_exposed_slurry_by_its_own_properties_only_with_the_flag`.
8. **The cascade's floor** (Task 1c, ruling 3) — once a macro step, the slurry's rate holds for slurry at least as warm
   as the cell that gave it: within a macro step 0.36 g on the heated state where a per-cell evaluation each flow
   sub-step would take up to 1.56 g. Pinned by
   `test_the_cascade_takes_each_cell_at_its_own_phase_s_rate_in_time_order_liquid_first_slurry_only_with_the_flag`.
9. **The release velocity of a skin record** (Task 2, ruling 5) -- the mean velocity of the metal released, the momentum
   the books lose with the droplets (spec 6.4), not Girin's closure velocity, which stays the gates' (Task 1c): on the
   test sphere 0.14-0.24 of the closure's, 1.04 against 5.6 m/s by mass. The rows' momentum is the droplets' to 4.7e-17.
   Confirmed by Asha, 2026-10-29 (option A over the closure velocity). Pinned by
   `test_a_spraying_skin_body_s_records_close_on_its_books_and_the_coupled_writer_writes_them`.
10. **The histograms against the records** (Task 2, ruling 2) -- the same bins, to numpy's weighted-histogram round-off,
    which is the whole step's (1.8e-17 of its droplets), not the bin's (up to 3.4e-11 of a small one): "exact" in spec 8
    is this. Pinned by the same test and
    `test_a_plan_s_records_are_one_per_column_and_mode_with_the_released_cells_solid_fraction_enthalpy_and_velocity`.
11. **The element model's two new columns** (Task 2, ruling 8) -- 0 and its film's liquid enthalpy, the sum of dm h_kg
    its removed enthalpy to the last bit; its 115 history columns, its other 24 source columns, its results and its
    summary and size distribution unchanged to the last bit on the short melting case. Pinned by
    `test_film_spraying_death_and_balances`.
12. **The skin flight's files** (Task 2, rulings 9 and 10) -- `particles.npz` and `size_distribution.csv` written by the
    element model's writer from one block of rows a macro step; the summary waits for Task 4's history columns, and the
    CLI's call for Task 3. The rows' volume at the continuum step is an input to settle (below).
13. **The element model without stage 1's skins** (Task 3, rulings 1 and 2) -- every branch removed was dead without
    skins, so on the CLI test's short melting case the 93 history columns kept, the source table, the summary, the size
    distribution, the settings and the 61 results kept are Task 2's to the last bit; the 22 history columns and eight
    results removed held 0 or nan in every row (fact 177). Pinned by `test_the_one_dimensional_skins_are_retired` and
    the element model's own tests.
14. **The default flips** (Task 3) -- a thermal run with no flags now runs the skins: its name loses `_surface-skin` and
    an element run's gains `_surface-elements` (fact 178), and every script or test of the element model names it. A
    Step 2 run with no flags now costs the skins' price (fact 147), and on the CLI tests' coarse mesh exits 1 at the
    first move unless the fine band is coarsened (`--fine-size 2`, as the skin tests give). Pinned by
    `test_the_run_names_mark_the_element_model_and_the_skins_non_default_settings`.
15. **The removal follows the flow on the skins** (Task 3, ruling 3) -- instant removal with the flow off, Girin's with
    it on, so a melting run with no flags runs the heat-only skins while the flow's default is off (an input to settle;
    since Task 3b the flow is on by default, item 18). Pinned by
    `test_the_skins_are_the_default_and_girin_s_removal_follows_their_flow`.
16. **Contradictions exit 2** (Task 3, ruling 4) -- Girin's removal without the flow, `--slurry-spray on` without it,
    and the element model's settings given to the skins, which would otherwise run the skins silently. Pinned by
    `test_contradictory_skin_settings_exit_2`.
17. **The spraying skin flight from the command line** (Task 3) -- the shared spray model with the command line's
    constants, the droplet files and a summary recorded as not written; the books close (fact 179). Pinned by
    `test_a_skin_run_with_flow_sprays_from_its_top_cells_and_writes_the_droplet_files` and
    `test_the_command_line_builds_the_skins_with_the_shared_spray_model`.
18. **The flow on by default** (Task 3b, Asha 2026-10-29) -- a melting skin run with no flags flows and sprays; an
    explicit `--removal instant` selects the heat-only skins (`_skinflow-off`) as an explicit `--skin-flow off` does
    (ruling 1), and `skinbody.FLOW_DEFAULT` stays the class's heat-only default (ruling 2). Pinned by
    `test_the_skins_are_the_default_and_girin_s_removal_follows_their_flow` and the name tests.
19. **The continuum step from the skins' first liquid** (Task 3b) -- the first macro step at whose end a top cell is
    above the liquidus latches it if the Knudsen switch has not; on the element model the body is never asked and its
    outputs are Task 3's to the last bit. At 71 km the step latches at the first liquid and the skins then spray
    nothing while liquid lies on slurry off Girin's closure (fact 181): the steps the trigger adds before the Knudsen
    switch, about 1 950 on each flight (fact 182), buy little spray where option B's hold applies. Pinned by
    `test_the_continuum_switch_latches_at_the_first_liquid_only_where_asked` and
    `test_the_skins_continuum_step_latches_at_their_first_liquid_before_the_knudsen_switch`.
20. **Register's release angle** (Task 4, Asha 2026-10-29) -- from the vector p - c, never the normal's `theta_deg`
    (they differ by -65 to +88 degrees on the element model's staircase); c the centroid of the volume the outer surface
    encloses, from the surface and at the instant of the rows' positions (ruling 1: the mass centre would move the angles
    by at most 0.088 degrees on the test sphere, fact 183); the history's mean, least and largest the step's rows', nan
    where nothing was released. Pinned by
    `test_register_s_release_angle_is_180_at_the_stagnation_point_90_at_the_equator_and_0_behind_from_the_actual_vector`,
    `test_the_release_angle_rides_in_the_source_table_and_a_step_s_statistics_are_by_mass` and the two bodies' tests.
21. **The element model's outputs under Task 4** -- three history columns and one source column appended, every other
    column, field, file and result its own to the last bit on the CLI test's short melting case, its vanishing last step
    kept (ruling 4, fact 185). Pinned by `test_film_spraying_death_and_balances` and the melting run's CLI test.
22. **A skin flight's summary** (Task 4, rulings 2 and 3) -- the element model's droplet columns under their names and
    meanings, the step's released mass from the new `released_step_mass_kg` (the skin history's `released_mass_kg` stays
    since time zero), the skins' liquid and flowing depths where the element model has its film. Pinned by
    `test_the_writer_writes_a_skin_flight_s_summary_from_its_droplet_columns_and_a_spool_as_it_writes_rows`.
23. **The skins' last step** (Task 4, ruling 4) -- with `snap_to_t_max` a step ending within 1e-6 of a step short of t_max
    ends there; without it `Simulator.advance` is as it was. On the 65 km smoke flight the vanishing step had written
    191 of its 1 196 rows (3.4e-17 kg) and read Pi 5.50 (fact 185). Pinned by
    `test_the_skins_loop_ends_at_t_max_without_a_vanishing_step_and_the_element_model_s_is_unchanged`.
24. **The records on disk** (Task 4, ruling 9) -- a spraying skin run's blocks appended to one file per column and read
    back column by column: the same files to the last bit, the writer's allocations about 41 bytes a row against 460
    from memory (fact 186). Pinned by `test_a_spool_on_disk_gives_back_the_blocks_it_was_given_to_the_last_bit` and the
    skin body's records test.
25. **The FEniCSx tier's thick-film test** (Task 4, ruling 5) -- both solvers' gradients at 1e-13, its 1e-8 bound kept:
    the failure since Task 1a was scikit-fem's linear tolerance, not the books (fact 187).
## Amendment of 2026-10-29 (after Task 3) — the flow on by default, and the continuum step from the skins' first liquid (Asha's decisions)

> Decided by Asha on 2026-10-29 on the two questions Task 3 left (spec 10.1 and 16 as amended that day). This adds Task
> 3b; the element model is unchanged.

1. **The flow is on by default for melting skin runs** (option B of Task 3's report). With `--thermal fem --melt on`
   and the skins, `--skin-flow` not given means on: the run flows, sprays by Girin's model and writes the droplet files.
   The heat-only skins are `--skin-flow off`, named `_skinflow-off`, with instant removal as before; a flow run's name
   loses `_skinflow-on`; a run that does not melt is as it was.
2. **The continuum step also starts at melting onset, on the skins only.** On a skin run with a continuum step, the
   macro step latches to it at the first macro step at whose end the skins hold liquid above the liquidus, or at the
   Knudsen switch, whichever comes first; the results record which, with its flight time, Knudsen number and altitude.
   The element model's switch stays exactly as it is.

## Amendment of 2026-10-29 — how the skins' depths are read where liquid lies on slurry, and the slurry cascade (Asha's decisions)

> Decided by Asha on 2026-10-29, on the two choices Task 1c's first writing stopped on (spec §8 and §16 as amended that
> day). This changes Task 1c and adds an optional argument to the shared spray physics; the element model is unchanged.

1. **Decision 1, option B: the element model's reading.** Where liquid lies on slurry, the skins' front-surface
   Rayleigh–Taylor mode reads the liquid depth and the molten region a wave must fit in is the liquid's, and off Girin's
   closure a liquid on slurry is held from the shear modes, as the rigid-substrate amendment of 2026-10-06 has them in
   the element model; b, the regime test and the viscous film's floor stay on the flowing depth. Through a small
   optional argument in the shared spray code (`rt_layer`) that the element model never passes.
2. **Decision 2, option i: the slurry cascade as spec §8 has it**, the 2026-10-07 design's cascade with its slurry
   properties, behind `--slurry-spray on`, off by default, with its tests, and where it releases measured and recorded
   in spec §8 as a measured fact, not a change of design. (Measured in Task 1c: it releases only exposed slurry within
   about 0.6–0.8 K of the liquidus, which on a heated surface is where the slurry under the liquid lies -- not inert.)
3. **The two choices made within the spec stand:** Girin evaluated once a macro step, at its start, only what the
   admitted run holds re-checked in each flow sub-step; Girin's surface velocity his closure's.

## Amendment of 2026-10-28 — the element model's spraying amendments of 2026-10-22 come into this line in Task 1 (Asha's decision)

> Decided by Asha on 2026-10-28, as this sub-plan's code-level writing began (planning copy
> `prototype/work-2026-10-28-skin-spraying-plan/`, from sub-plan 21's tested code as committed at 0f0425e). This
> changes Task 1 and, through the shared spray physics, the element model on this line.

**What comes across.** The two amendments committed on branch `we-critical-308` (cut from `main` at 612b916, not
merged), with their tests, before anything else in Task 1:

- **fa02a4c** (its facts 108–110): Girin's theoretical critical surface Weber number, 3.08, is the default
  (`--we-critical`), against 4.62 before. On the element model's 100 mm flight: sprayed mass +0.13 %, medians under 1 %,
  droplet count +8 %. Its code: `prototype/work-2026-10-22-we-critical/amendment/*.diff` on that branch.
- **a54520c** (its facts 111–115): Girin & Kopyt's thin, rarefied and regime-2 branches no longer borrow Girin's
  critical Weber number, which their Eq. 5 does not have; they are gated by the wave fit and a viscous film on the rigid
  floor (`spray.viscous_film`: Lamb's damping plus Hunt's Stokes layer on the non-rigid depth, Γτ* ≥ 1 blocks; a depth
  under sqrt(2ν/ω) is rejected and recorded in the history and results). On the element model's 100 mm flight: sprayed
  mass +1.4 %, the thin branch 0.37 rising to 0.87 kg and the thick falling from 0.73 to 0.25 kg, the droplet count
  83 rising to 158 million and the median radius by mass falling from 188 to 111 µm. Its code:
  `prototype/work-2026-10-22-thin-gate/amendment/*.diff` on that branch.

Then the per-crest cap (the amendment of 2026-10-24 below) and the spraying from the skins' top cells are built on
them. On the skins the viscous film's floor depth is the skin's flowing (non-rigid) depth.

**Why.** Task 5's flights should spray by the rules of record. Left for the merge, the thin gate alone would move the
flights' droplet sizes by about 40 % and every flight would have to be flown again.

**What it changes elsewhere.** `--surface-model elements` on this line becomes the element model of `we-critical-308`,
so Task 3's bit-for-bit check of the element model is against that branch's flights, not `main` at 612b916. The facts
on that branch (108–115) collide with this line's numbering (101–160) and are renumbered when brought across. Its open
points (its fact 115) travel with it: 22 g of thin spray within 10 degrees of the stagnation point; the mass the
viscous film holds back (Γτ* ≥ 1) is not recorded; the ASCEND abstract's numbers come from the gated flight; subsonic
patches lie outside Girin & Kopyt's derivation; the 50 mm flight was not rerun.

## Amendment of 2026-10-24 — the droplet-size cap is per crest, on the depth the release can take (Asha's decision)

> Decided by Asha on 2026-10-24, after the element model's flights on branch `we-critical-308` (its facts 108–115) showed
> how often the old cap set droplet sizes. This changes Task 1 and the shared spray physics; it is not yet code.

**The old cap.** `SprayModel.evaluate` limits every branch's droplet radius by the film on the patch,
r ≤ (3 m_f / (4π ρ_l))^⅓, and by a quarter of the body radius. In skyscrapers there is no film account (spec §6.1), and
the patch's film depends on the mesh facet's area, not on the flow. Measured on the element model's 100 mm flight with
the thin-gate amendment: it bound on 67 % of the thick branch's release rows (24 % of its mass), 16 % of the thin
branch's (4.4 %), and set essentially every Rayleigh–Taylor droplet.

**The new cap (decided).**

- **Per crest.** A droplet holds no more liquid than one crest of its own wave can give: (4/3)π r³ ≤ λ² h, so
  r ≤ (3 λ² h / (4π))^⅓, where λ is the wavelength of the mode that sprays the column (the thick branch's
  2π δ_m / Δ_f, Girin & Kopyt's λ* on the thin, rarefied and regime-2 branches, the Rayleigh–Taylor λ*) and h is the
  depth below. λ² because a crest repeats every wavelength in both surface directions. The column's area, a mesh
  corner's dual area, no longer enters.
- **On the depth the release can take.** h is the depth of the run of cells the release is allowed to remove, the same
  `admit` run `Skins.take_top_mass` is given (spec §8):
  - `--slurry-spray off` (the default): the liquid run from the top, the skin's **liquid depth** of spec §6.7 (the
    fully liquid cells contiguous from the top, interpolated within the boundary cell), which is shallower than the
    non-rigid layer wherever slurry lies beneath the liquid;
  - `--slurry-spray on`: the whole **flowing (non-rigid) depth**, from the top to the first point at T_rigid.
- **Once per macro step.** The cap is evaluated at the start of the macro step, from that step's gas state, δ_m and the
  column's depth at the start, and holds for every sub-step of it, so spec §8's statement that a column's radius is
  fixed within a macro step (and that combining the sub-steps' records changes no radius) still holds. The release in
  each sub-step still takes at most what the admitted run then holds.
- **Unchanged:** the quarter-body-radius cap.

**Where it goes.** The cap leaves `m_f`: the shared spray physics takes the cap's depth as its own argument (for
example `cap_depth`), and with it absent keeps the film-mass cap, so `--surface-model elements` keeps the current model
bit for bit (Task 3). The skyscrapers call passes the admitted run's depth from the start of the macro step.

**To check when Task 1 is written at code level:** the Rayleigh–Taylor droplets, whose radius is their wavelength
(centimetres) and which the film cap held to about 2 mm on the element model's flight, will be larger under a
per-crest cap on a millimetre-deep layer (about 4.6 mm at λ = 2 cm, h = 1 mm); and the tests of Task 1 gain the cap's
two depths, its independence of the column's area, and its constancy over a macro step's sub-steps.

**Related, not part of this decision:** the element model's two amendments of 2026-10-22 on `we-critical-308` (Girin's
critical surface Weber number 3.08 by default; Girin & Kopyt's branches gated by the wave fit and a viscous film on the
rigid floor, whose floor depth maps onto the skin's flowing depth) change the same spray module and reach this line
only when `main` is merged into it.

**Written at code level (2026-10-28).** The related amendments came across first, in Task 1a, by the amendment of
2026-10-28 above, not with the merge. The cap is Task 1b's `SprayModel.evaluate(..., cap_depth=...)`, the element model
unchanged to the last bit. The check: 4.57 mm at λ = 2 cm and h = 1 mm, against 1.8–2.6 mm under the film cap on the
element model's thin-gate flight (fact 170). The tests: the area in Task 1b, the two depths and the constancy over a
macro step's sub-steps in Task 1c.

## Amendment of 2026-10-12 — port first, move into `shared/` afterwards (Asha's decision; supersedes the paths of the 2026-10-11 amendment below)

> Design of record: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4.2, §4.4 and §14 as revised on 2026-10-12. **When this sub-plan's code is ported from
> the prototype into the package, use the module paths the prototype uses (`reentry_model/<module>.py`), exactly as the
> body and the earlier amendments say, and ignore the `reentry_model/shared/...` paths of the 2026-10-11 amendment
> below.** There is no `reentry_model/shared/` package on `main`. After the skyscrapers port has landed on `main`, one
> combined refactor (the level-set design's stage 0 + S, with its own plan) moves the architecture-independent modules
> into `reentry_model/shared/`, keeps every old path working as an alias of the same module object, and checks that
> every flight's history is unchanged. The 2026-10-11 amendment stays as the record of **which** of this sub-plan's
> modules will be shared; only **when** and **where they are created** is superseded. Helpful, not required: where it
> costs nothing, keep the shared physics free of skyscrapers-only imports, so that the later move stays mechanical.

**Shared later from this sub-plan:** the spray physics and the droplet outputs it calls; the bookkeeping of the top cells is skyscrapers' own.

## Amendment of 2026-10-11 — the shared physics package (Asha's instruction)

> Design of record: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4 and §14 (stages 0 and S). The repository now has two model architectures,
> **skyscrapers** (this plan's) and **level set**, chosen by `reentry_model run --architecture skyscrapers|levelset`
> (default `skyscrapers`, whose runs, run names and results are unchanged). Code that both architectures use lives in
> `reentry_model/shared/`; each moved module's old path stays importable as a one-line re-export, so nothing that
> imports the old path breaks. **This amendment changes where files live, not what they do.** It takes precedence over
> the file paths in the body below and in earlier amendments. Test files keep their names and places and import the
> module from `reentry_model.shared`. The skins are never shared between the architectures. Read the matching amendment
> at the end of `00-shared-context.md` first.

**What changes in this sub-plan.**

- Girin's spraying from the top cells uses the shared physics in **`reentry_model/shared/spray.py`** (S4) and writes
  droplets through **`reentry_model/shared/droplets.py`** (S6); the bookkeeping of the top cells and the removal of the
  old film and deep account are skyscrapers' own.

---

### Task 1a: The element model's spraying amendments of 2026-10-22 brought across

**Files:**
- Modify: `reentry_model/spray.py` (`SprayModel`'s default critical Weber number; `viscous_film`; the gates on Girin &
  Kopyt's branches; `SprayResult.h_floor`, `.gamma_tau`, `.outside`; `SOURCE_COLUMNS` and `source_rows`, 24 columns),
  `reentry_model/dispersion.py` (the constants' comments), `reentry_model/cli.py` (`--we-critical`'s default and help),
  `reentry_model/body.py` (`MeltingBody.thin_invalid_patch_steps`, `.thin_invalid_mass_max`, the step's record and the
  history row), `reentry_model/coupled.py` (two history columns, two results), `analysis/melt_sensitivity.py` (the
  variant `we462` in place of `we308`)
- Test: `tests/test_reentry_model_spray.py` (four tests added, one replaced, five changed),
  `tests/test_reentry_model_melting.py` (one added, two changed), `tests/test_reentry_model_cli.py` (the settings'
  3.08)

**Interfaces:**
- Consumes: the element model's `SprayModel.evaluate(flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius,
  b_layer, extent, deceleration_n, regime_layer, on_slurry)` and `thin_film_mode`, unchanged in their arguments.
- Produces: `spray.SprayModel(liquid, ..., we_critical=dispersion.WE_CRITICAL_THEORY)` (3.08) and `--we-critical`
  defaulting to it, now the thick branch's gate only; `spray.viscous_film(lam, tau, h, liquid) -> (gamma_tau, outside)`;
  `SprayResult.h_floor` [m], `.gamma_tau` and `.outside` (bool), nan / False off Girin & Kopyt's branches;
  `SOURCE_COLUMNS` + `h_floor_m`, `gamma_tau`; `MeltingBody.thin_invalid_patch_steps` and `.thin_invalid_mass_max`
  [kg]; history `thin_invalid_patches`, `thin_invalid_mass_kg`; results `thin_invalid_patch_steps`,
  `thin_invalid_mass_max_kg`. `dispersion.WE_CRITICAL_PRACTICAL` (4.62) stays for `girin_case`.

What comes across is branch `we-critical-308`'s two commits (spec §16, 2026-10-28; facts 161–168), with their tests:
fa02a4c makes Girin's theoretical critical surface Weber number, 3.08, the default against his practical 4.62, since the
wave fit and the growth period already represent what 4.62 stands in for (fact 161); a54520c takes that threshold off
Girin & Kopyt's thin, rarefied and regime-2 branches, whose Eq. 5 has none (fact 164), and gates them by the wave
fitting in the molten region and by `spray.viscous_film` — the dominant wave's damping on a film over a rigid floor,
Γ = 2νk² + k√(νω/2)/sinh(2kh), Lamb's interior term plus Hunt's Stokes layer, the patch releasing nothing where
Γτ* ≥ 1 and a depth under the floor's Stokes depth √(2ν/ω) rejected as outside the inviscid mode and recorded (fact
165). The floor depth is the non-rigid depth the regime test reads; on the skins (Task 1c) the flowing depth.

**How it was carried** (ruling 1). The twelve diffs of the two `amendment/` folders apply with `patch` to this line's
code but for three hunks, where this line's code differs from that branch's base for reasons outside the amendments:
`MeltingBody._film_and_spray`'s step record here ends in sub-plan 18's resolved depths (`**depths`), and the two keys
go before them; `coupled.MELT_COLUMNS` here ends at `slurry_held_mass_kg`, because `main`'s frame export of 2026-10-08
(that branch's facts 101–107: `p_w_stag_step_Pa`, `v_hat_body`) never reached this line, and the two columns follow
it; and the CLI test's settings line here has no `v_hat_body` check, so only its 4.62 becomes 3.08. `spray.py`,
`dispersion.py`, `analysis/melt_sensitivity.py` and the spray and melting tests are then that branch's tested copies
to the byte, but for one comment, `dispersion.WE_CRITICAL_THEORY`'s, which cites this line's fact 161 for that
branch's 108 (ruling 2). Its sub-plan amendments (04, 07, 09, 10, 13, 14, 15) describe `main`'s element-model plans
and come with the merge (ruling 3).

**Measured** (2026-10-28, `prototype/work-2026-10-28-skin-spraying/measure/task1/cmp/`; fact 169). The CLI test's
short melting case (`run_case.sh`: 15 s from 71 km at 7.24 km/s, 800 K, physics heating, the continuum step at 0.1 s;
76 history rows, about 30 s a run), on this line with the amendments and on that branch's tested copy
(`prototype/work-2026-10-22-thin-gate/code/` at a54520c, copied): every one of the 92 history columns they share, all
24 fields of the 25 266 source rows and every result they share are the same to the last bit (`compare_runs.py`). The
columns only this line has are sub-plan 18's 23 skin and resolved-depth columns; the one only that branch has is
`p_w_stag_step_Pa`. The case reaches every branch the amendments touch that liquid aluminium can reach: the thick
branch (12 899 rows, 0.190 kg; 242 of them, 5.0 g, between 3.08 and 4.62, which the old default held back), the thin
branch (11 722 rows, 0.113 kg; 6 985, 8.5 g, below 3.08, which the borrowed gate held back), the viscous film (Γτ* up
to 1.0, median 0.19, the floor's median 132 µm), its validity reject (28 590 patch-steps, at most 0.069 g of film in a
step) and the Rayleigh–Taylor mode (645 rows, 0.213 kg). The rarefied branch needs melting under free-molecular flow and
regime 2 never sprays liquid aluminium (fact 165); both are reached by the unit tests. Against this line before the
amendments the case sprays 0.51611 kg against 0.51587 kg, its thin branch 11 722 rows against 4 586 (0.113 against
0.103 kg) and 8.15 million droplets against 7.63 million.

**Rulings** (2026-10-28; each what, why, cost if wrong):
1. *The three hunks by hand, carried in meaning.* Why: the bases differ outside the amendments, and the amendments'
   intent — two keys in the step record, two history columns at the end of the melt columns, the default in the test —
   does not depend on the context they were written in. Cost if wrong: the history's column order differs from that
   branch's; columns compared by name are identical.
2. *The fact cited in `dispersion.py`'s comment renumbered.* Why: this line's fact 108 is another fact. Cost: none.
3. *Only the code, tests and facts come across, not that branch's sub-plan amendments.* Why: those amend `main`'s
   element-model plans (04–15), which this line's plans do not regenerate; the merge brings them. Cost if wrong: a
   reader of sub-plan 07 on this line finds the old gate described there until the merge.

- [ ] **Step 1: The amendments' spray tests** in `tests/test_reentry_model_spray.py`:

```diff
--- a/tests/test_reentry_model_spray.py
+++ b/tests/test_reentry_model_spray.py
@@ -74,7 +74,7 @@
     res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True, True]), 0.5, areas, m_f)
     assert res.branch.tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN, -1]
     we_s = LIQ.rho * 25.0 * delta_m[0] / LIQ.sigma
-    assert res.we_s[0] == pytest.approx(we_s) and res.unstable[0] == (we_s > 4.62)
+    assert res.we_s[0] == pytest.approx(we_s) and res.unstable[0] == (we_s > 3.08)
     d_f, im_f, _ = dispersion.DispersionTable()(we_s)
     lam = 2 * np.pi * delta_m[0] / d_f
     assert res.r[0] == pytest.approx(min(0.17 * lam, (3 * m_f[0] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
@@ -92,8 +92,11 @@
     res = model.evaluate(flow_c, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
     assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]          # a thick film cannot take Girin's branch here
     # free-molecular branch: no edge state exists, the freestream drives the mode
+    # (a floor a millimetre down: the rarefied lambda* of 7.5 mm has a Stokes layer of 0.14 mm on its floor, deeper
+    # than a 20 um film, which would put that film outside the inviscid mode's validity -- spray.viscous_film)
     flow_fm = FakeFlow(2, closure=sf.CLOSURE_COUETTE, branch=sf.BRANCH_FREE_MOLECULAR)
-    res = model.evaluate(flow_fm, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
+    res = model.evaluate(flow_fm, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9,
+                         b_layer=np.full(2, 1e-3), regime_layer=np.full(2, 1e-3))
     assert res.branch.tolist() == [spray.BRANCH_RAREFIED, spray.BRANCH_RAREFIED]
     lam_r, _ = spray.thin_film_mode(np.array([24.0]), np.array([8e-5 * 7200.0 ** 2]), LIQ)
     assert res.r[0] == pytest.approx(min(lam_r[0] / 4.0, (3 * (m_f[0] + 1e-9) / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
@@ -122,25 +125,95 @@
     assert spray.BIN_EDGES[0] == 1e-6 and spray.BIN_EDGES[-1] == 1e-2
 
 
-def test_girins_critical_weber_number_gates_the_thin_branch_too():
-    """Girin & Kopyt's inviscid side mode is unstable at every wavelength and carries no threshold of its own, so the
-    thin branch used to strip film wherever a wavelength existed -- at any angle, however small the shear. Girin's
-    (2017) criterion We_s > We_cr now gates it as well, measured over the liquid the shear can actually reach,
-    min(delta_m, layer): on a film shallower than the conjugate layer that is the film's own depth."""
+def test_the_thin_branch_is_not_gated_by_girins_critical_weber_number():
+    """Girin & Kopyt's Eq. (5), solved, has no critical Weber number: the growth rate is positive at every We, and the
+    threshold it has once the wavelength is bounded is the wave fitting in the molten region (2026-10-22). So Girin's
+    (2017) We_s > We_cr, which gated the thin branch from 2026-09-24, no longer does; We_s is still reported."""
     flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)                # no delta_m: the thin branch
     model = spray.SprayModel(LIQ)
     delta_m, _ = spray.melt_layer(flow, LIQ)
     b = np.array([2e-5, 2e-5])
     areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
-    v_s = np.array([0.5, 20.0])                                   # subcritical, supercritical
+    v_s = np.array([0.5, 20.0])                                   # below and above Girin's 3.08
     res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True]), 0.5, areas, m_f)
     we = LIQ.rho * v_s ** 2 * b / LIQ.sigma                       # layer = b here, and delta_m is nan
-    assert res.we_s == pytest.approx(we) and we[0] < 4.62 < we[1]
-    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]          # classified either way
-    assert not res.unstable[0] and res.dm[0] == 0.0 and np.isnan(res.r[0])        # but stable: nothing leaves
-    assert res.unstable[1] and res.dm[1] > 0.0                                    # and the supercritical one sprays
-    lam, tau = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
-    assert np.isfinite(lam[0]) and np.isfinite(tau[0])            # a wavelength existed for both -- only We decided
+    assert res.we_s == pytest.approx(we) and we[0] < 3.08 < we[1]
+    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]
+    assert res.unstable.all() and (res.dm > 0.0).all()                            # both spray: We_s does not decide
+    assert res.r[0] == pytest.approx(res.r[1])                    # the mode is the gas's: same lambda*, same droplet
+    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True]), 0.5, areas, m_f,
+                         extent=np.array([1e-4, 1e-2]))           # lambda* = 0.46 mm: fits only in the 10 mm region
+    assert not res.unstable[0] and res.dm[0] == 0.0 and res.unstable[1] and res.dm[1] > 0.0
+
+
+def lamb(lam, tau, liquid):
+    return 2.0 * liquid.mu / liquid.rho * (2.0 * np.pi / lam) ** 2 * tau
+
+
+def capillary_tau(lam, liquid):
+    return spray.CAPILLARY_TAU * lam ** 1.5 * np.sqrt(liquid.rho / liquid.sigma)
+
+
+def test_on_deep_liquid_the_film_damping_is_lambs():
+    """Far above the floor (k h >> 1) Hunt's bottom term vanishes and gamma tau* is Lamb's 2 nu k^2 tau* alone, which
+    binds only below a wavelength of about 3.3 um for this liquid."""
+    lam = np.array([1e-6, 3.0e-6, 3.6e-6, 1e-4])
+    tau = capillary_tau(lam, LIQ)
+    gt, out = spray.viscous_film(lam, tau, np.full(4, 1.0), LIQ)
+    assert gt == pytest.approx(lamb(lam, tau, LIQ), rel=1e-12) and not out.any()
+    assert (gt >= 1.0).tolist() == [True, True, False, False]
+
+
+def test_a_thin_layer_is_damped_more_than_lamb_predicts():
+    """On a layer of a few wavelengths' depth or less the floor's Stokes layer adds k sqrt(nu omega/2)/sinh(2kh), and
+    omega itself falls with tanh(kh); the layer is still deeper than that Stokes layer, so the mode applies."""
+    lam = np.full(3, 1e-4)
+    tau = capillary_tau(lam, LIQ)
+    h = np.array([1.0, 3e-5, 1e-5])                            # deep, k h = 1.9, k h = 0.63
+    gt, out = spray.viscous_film(lam, tau, h, LIQ)
+    assert not out.any()
+    assert gt[0] == pytest.approx(lamb(lam[0], tau[0], LIQ), rel=1e-12) and gt[0] < gt[1] < gt[2]
+    k = 2 * np.pi / lam[2]
+    omega = np.sqrt(LIQ.sigma * k ** 3 / LIQ.rho * np.tanh(k * h[2]))
+    nu = LIQ.mu / LIQ.rho
+    assert gt[2] == pytest.approx((2 * nu * k ** 2 + k * np.sqrt(nu * omega / 2) / np.sinh(2 * k * h[2])) * tau[2], rel=1e-12)
+
+
+def test_a_layer_thinner_than_the_stokes_depth_is_rejected_and_flagged():
+    """Where h < sqrt(2 nu / omega) the floor's boundary layer fills the liquid and Girin & Kopyt's inviscid mode does
+    not apply: nothing leaves by it, and the patch is flagged so the body can count it and the film it holds."""
+    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)                # no delta_m: the thin branch, lambda* 0.46 mm
+    model = spray.SprayModel(LIQ)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.array([2e-6, 5e-5])                                    # 2 um: inside the Stokes layer; 50 um: outside it
+    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
+    res = model.evaluate(flow, FakeState(), b, delta_m, np.full(2, 20.0), np.array([True, True]), 0.5, areas, m_f)
+    lam, tau = spray.thin_film_mode(flow.mach_e, flow.rho_e * flow.u_eff ** 2, LIQ)
+    k = 2 * np.pi / lam[0]
+    omega = np.sqrt(LIQ.sigma * k ** 3 / LIQ.rho * np.tanh(k * b))
+    stokes = np.sqrt(2 * LIQ.mu / LIQ.rho / omega)
+    assert b[0] < stokes[0] and b[1] > stokes[1]
+    assert res.branch.tolist() == [spray.BRANCH_THIN] * 2 and res.outside.tolist() == [True, False]
+    assert res.dm[0] == 0.0 and res.dm[1] > 0.0
+    assert res.h_floor == pytest.approx(b) and np.isfinite(res.gamma_tau).all()
+
+
+def test_the_damping_uses_the_non_rigid_depth_not_the_film():
+    """The floor is where the liquid stops: the non-rigid depth the regime test reads, not the film thickness b."""
+    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
+    model = spray.SprayModel(LIQ)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.full(2, 5e-6)                                          # the same film on both
+    regime = np.array([5e-6, 1e-3])                               # floor right under it, floor a millimetre down
+    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
+    res = model.evaluate(flow, FakeState(), b, delta_m, np.full(2, 20.0), np.array([True, True]), 0.5, areas, m_f,
+                         b_layer=regime, regime_layer=regime)
+    lam, tau = spray.thin_film_mode(flow.mach_e, flow.rho_e * flow.u_eff ** 2, LIQ)
+    gt, out = spray.viscous_film(lam, tau, regime, LIQ)
+    assert res.branch.tolist() == [spray.BRANCH_THIN] * 2         # no delta_m: thin however deep the slurry
+    assert res.h_floor == pytest.approx(regime) and res.gamma_tau == pytest.approx(gt, rel=1e-12)
+    assert res.outside.tolist() == [True, False] and res.dm[0] == 0.0 and res.dm[1] > 0.0
+    assert gt[1] < gt[0]                                          # the deeper floor damps less
 
 
 def test_a_wave_longer_than_its_melt_region_cannot_form():
@@ -188,19 +261,22 @@
     Girin's critical angle gating the shear branches the stagnation cap otherwise has no mechanism at all. The mode is
     applied where it passes *both* criteria -- the depth threshold W cos(phi) h^2 rho_l > 3 Sigma, which is what confines
     it to the nose, and an admissible wave in the molten region -- and where its growth time beats the shear mode's,
-    since two modes of one interface cannot both break it up."""
-    flow = FakeFlow(3, closure=sf.CLOSURE_COUETTE)            # no delta_m: the shear mode here is the weak thin one
+    since two modes of one interface cannot both break it up. Near the stagnation point the edge velocity vanishes and
+    the thin mode's lambda* = 1.5 M Sigma / (rho V^2) outgrows the molten region (2026-10-22: no longer Girin's We_s)."""
+    flow = FakeFlow(3, closure=sf.CLOSURE_COUETTE)            # no delta_m: the shear mode here is the thin one
+    flow.u_eff = flow.u_e = np.array([5.0, 5.0, 2000.0])     # at the nose, at the nose, on the side
     model = spray.SprayModel(LIQ)
     delta_m, _ = spray.melt_layer(flow, LIQ)
     b = np.full(3, 1e-4)
     layer = np.array([5e-3, 1e-3, 5e-3])                      # deep, too shallow for the threshold, deep
-    v_s = np.array([0.1, 0.1, 40.0])                          # stable shear, stable shear, strong shear
+    v_s = np.array([0.1, 0.1, 40.0])
     w_n = np.full(3, 95.0)                                    # at the stagnation point cos(phi) = 1
     areas, m_f = np.full(3, 1e-5), b * LIQ.rho * 1e-5
+    extent = np.full(3, 0.05)                                 # a 5 cm molten cap: the thin lambda* is 7.4 cm at 5 m/s
     deep, _, _ = spray.rayleigh_taylor(w_n, layer, LIQ)
     assert deep.tolist() == [True, False, True]               # h > lambda*/2 pi only for the deep pair
     res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f,
-                         b_layer=layer, deceleration_n=w_n)
+                         b_layer=layer, extent=extent, deceleration_n=w_n)
     assert res.branch.tolist() == [spray.BRANCH_RT, spray.BRANCH_THIN, spray.BRANCH_THIN]
     assert res.dm[0] > 0.0 and res.dm[1] == 0.0 and res.dm[2] > 0.0
     _, lam_rt, tau_rt = spray.rayleigh_taylor_bounded(w_n, layer, np.full(3, np.inf), LIQ)
@@ -211,7 +287,8 @@
     assert lam_rt[0] > 1e-2 and res.r[0] == pytest.approx(cap)   # lambda* is centimetres; the film-mass cap binds
     # with the mode reported but not applied -- every run before 2026-09-24 -- the deep stable patch releases nothing
     off = spray.SprayModel(LIQ, rt_spray=False).evaluate(
-        flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f, b_layer=layer, deceleration_n=w_n)
+        flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f, b_layer=layer, extent=extent,
+        deceleration_n=w_n)
     assert off.branch.tolist() == [spray.BRANCH_THIN] * 3 and off.dm[0] == 0.0 and off.dm[2] > 0.0
     assert off.dm[2] == pytest.approx(res.dm[2])              # the shear-dominated patch is untouched either way
 
@@ -219,7 +296,7 @@
 def test_the_applied_rayleigh_taylor_mode_needs_the_wave_to_fit_its_melt_region():
     """The lateral bound applies to this mode as it does to the shear ones: a pool narrower than half the marginal
     wavelength admits nothing at all, however deep it is."""
-    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
+    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE, u=1.0)     # at the nose: the thin lambda* is 1.8 m and fits nowhere
     model = spray.SprayModel(LIQ)
     delta_m, _ = spray.melt_layer(flow, LIQ)
     b, layer, v_s = np.full(2, 1e-4), np.full(2, 5e-3), np.full(2, 0.1)
@@ -247,27 +324,19 @@
     flow = FakeFlow(2)
     model = spray.SprayModel(MELT_VISCOUS)
     delta_m = np.full(2, 1e-4)                                 # a conjugate depth the film can exceed
-    b, v_s = np.full(2, 1e-3), np.full(2, 10.0)                # deep film; v_s large enough to clear We_cr
+    b, v_s = np.full(2, 1e-3), np.full(2, 10.0)                # deep film
     areas, m_f = np.full(2, 1e-5), np.full(2, 1e-3 * 2400.0 * 1e-5)
     res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f)
     assert res.branch.tolist() == [spray.BRANCH_REGIME2] * 2    # not BRANCH_THICK, and the RT mode does not take it
-    assert res.we_s[0] > model.we_critical                     # Girin's gate still applies, as on every branch
-    # the wavelength and growth time are Girin & Kopyt's Eqs. (11) and (12) -- the thin branch's pair, same mode
+    # its wave is Girin & Kopyt's Eqs. (11) and (12), the thin branch's pair -- and a melt more viscous than the gas
+    # damps it: Lamb's term alone, 2 nu k^2 tau* = 63 nu sqrt(rho/Sigma) lambda^-1/2, is above 1 for any nu above the
+    # smallest edge nu_gas measured (0.018 m2/s) at any wavelength below kilometres, and the floor only adds to it. So
+    # since 2026-10-22 the branch is classified but releases nothing, for any melt that can reach it.
     lam = 1.5 * flow.mach_e[0] * MELT_VISCOUS.sigma / (flow.rho_e[0] * flow.u_eff[0] ** 2)
     tau = spray.CAPILLARY_TAU * lam ** 1.5 * np.sqrt(MELT_VISCOUS.rho / MELT_VISCOUS.sigma)
-    assert res.growth[0] == pytest.approx(tau, rel=1e-12)
-    assert res.r[0] == pytest.approx(spray.K_R * lam, rel=1e-12)            # r = k_r lambda*, below both caps here
-    # the rate is Girin 2017's torus shedding on that wavelength, with no min(b, lambda*/8) factor: one torus of
-    # cross-section r = k_r lambda* per wavelength per growth time, mdot = rho_l pi r^2 / (lambda* tau*)
-    mdot = MELT_VISCOUS.rho * np.pi * (spray.K_R * lam) ** 2 / (lam * tau)
-    assert res.mdot[0] == pytest.approx(mdot, rel=1e-12)
-    # "not depth-limited" means the rate does not reference the film depth at all, as on the thick branch: doubling the
-    # film leaves it unchanged, where the thin branch's min(b, lambda*/8) factor follows a film shallower than the wave
-    deeper = model.evaluate(flow, FakeState(), 2.0 * b, delta_m, v_s, np.full(2, True), 0.5, areas, 2.0 * m_f)
-    assert deeper.branch.tolist() == [spray.BRANCH_REGIME2] * 2
-    assert deeper.mdot[0] == pytest.approx(res.mdot[0], rel=1e-12)
-    # mdot is the instability's demand; the film is what limits the release, exactly as on every other branch
-    assert mdot * areas[0] * 0.5 > m_f[0] and res.dm[0] == pytest.approx(m_f[0])
+    gt, out = spray.viscous_film(np.array([lam]), np.array([tau]), b[:1], MELT_VISCOUS)    # h = the film: no slurry here
+    assert gt[0] >= 1.0 and (res.h_floor == pytest.approx(b))
+    assert not res.unstable.any() and (res.dm == 0.0).all() and np.isnan(res.growth).all()
 
 
 def test_a_kelvin_helmholtz_melt_fails_girins_first_stage_on_this_body():
@@ -292,9 +361,10 @@
     delta_m, _ = spray.melt_layer(flow, liq)
     assert delta_m[0] / flow.delta_a[0] == pytest.approx((liq.rho / flow.rho_e[0]) ** (1.0 / 3.0), rel=1e-12)
     b = np.full(1, 1e-2)                                         # 10 mm of melt: a tenth of the sphere's radius
-    res = spray.SprayModel(liq).evaluate(flow, FakeState(), b, delta_m, np.full(1, 10.0), np.full(1, True), 0.5,
-                                        np.full(1, 1e-5), b * liq.rho * 1e-5)
+    res = spray.SprayModel(liq, rt_spray=False).evaluate(flow, FakeState(), b, delta_m, np.full(1, 10.0),
+                                                        np.full(1, True), 0.5, np.full(1, 1e-5), b * liq.rho * 1e-5)
     assert res.branch.tolist() == [spray.BRANCH_THIN]            # regime 1, the stabilised-core case, as Girin has it
+    assert res.dm[0] == 0.0                   # and a melt this viscous damps the thin mode's wave (2026-10-22)
 
 
 def test_regime_2_never_fires_for_liquid_aluminium():
@@ -365,6 +435,7 @@
         b = np.full(2, 2e-5)
         areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
         res = model.evaluate(flow, FakeState(), b, delta_m, np.full(2, 20.0), np.ones(2, bool), 0.5, areas, m_f,
-                             on_slurry=np.array([True, False]))
+                             on_slurry=np.array([True, False]), b_layer=np.full(2, 1e-3), regime_layer=np.full(2, 1e-3))
+        # (the floor a millimetre down keeps the rarefied 7.5 mm wave inside the inviscid mode's validity, viscous_film)
         assert res.branch.tolist() == [-1, code]
         assert res.dm[0] == 0.0 and not res.unstable[0] and res.dm[1] > 0.0
```

- [ ] **Step 2: The melting body's tests** in `tests/test_reentry_model_melting.py`:

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -9,6 +9,7 @@
 
 from reentry_model import body, heating, material, mesh, thermal
 from reentry_model.body import PHI_DEATH as PHI_DEATH_FOR_TEST
+from reentry_model.spray import BRANCH_THIN
 from test_reentry_model_coupled import MASS_100MM, simulator
 
 pytest.importorskip("cantera")
@@ -94,7 +95,10 @@
     assert b.material.h_liquid < h_out <= b.material.enthalpy(b.surface_temperature().max())       # own temperature, superheat and all
     assert stats["film_T_max_K"] == pytest.approx(b.surface_temperature().max()) and stats["film_T_max_K"] > b.material.T_liquidus
     rows = np.array(b.source_rows)
-    assert rows.shape[1] == 22 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
+    assert rows.shape[1] == 24 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
+    # 2026-10-22: viscous_film's floor depth and gamma tau*, finite on Girin & Kopyt's branches and nan elsewhere
+    gk = np.isin(rows[:, 8], [1, 2, 4])
+    assert np.isfinite(rows[gk, 22:24]).all() and np.isnan(rows[~gk, 22:24]).all() and (rows[gk, 23] < 1.0).all()
     # the histogram bins 1 um to 10 mm and drops the rest by design. With the rigid substrate on (amendment of 2026-10-06)
     # this body, 880 K throughout, is slurry to the core: off Girin's closure its film is held from the shear modes,
     # piles up, and the front-surface Rayleigh-Taylor mode sheds some of it above 10 mm (measured in two unseeded runs:
@@ -668,12 +672,33 @@
     b.solver.set_temperature(T_wall - (T_wall - b.material.T_rigid) * d / depth)
 
 
+def test_films_shallower_than_the_floors_stokes_layer_are_counted_not_sprayed(layered_mesh):
+    """At 70 km under Girin's closure, a 2 um film on the rigid core takes the thin branch, and almost everywhere the
+    liquid is shallower than the floor's Stokes layer (about 16 um at lambda* 0.5 mm): Girin & Kopyt's inviscid mode does
+    not apply there. Nothing leaves those patches by it, and the step records how many there were and the film they
+    held -- a validity limit, not a finding that nothing leaves (2026-10-22)."""
+    b = melting_body(layered_mesh, rigid_substrate=True)
+    sim, a = girin_state(b)
+    flow, delta_m = flow_and_delta_m(b, a)
+    b.m_f = np.where(b.windward, 2e-6 * b.liquid.rho * b.surface.areas, 0.0)
+    h_p = deep_stage_args(b, delta_m)[4]
+    b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
+    res, lm = b.last_spray, b.last_melt
+    wet = b.windward & (b.last_b >= 1e-6)
+    assert (res.branch[wet] == BRANCH_THIN).all()
+    flagged = res.outside & b.windward
+    assert flagged.sum() > 0.9 * wet.sum() and (res.dm[flagged] == 0.0).all()
+    assert lm["thin_invalid_patches"] == flagged.sum() == b.thin_invalid_patch_steps
+    assert lm["thin_invalid_mass_kg"] == pytest.approx(b.m_f[flagged].sum()) and lm["thin_invalid_mass_kg"] > 0.0
+    assert b.thin_invalid_mass_max == lm["thin_invalid_mass_kg"]
+
+
 def test_a_film_on_slurry_deeper_than_the_conjugate_depth_is_thick_in_the_melt_step(layered_mesh):
     """At 70 km, under Girin's closure on every windward patch (conjugate depth 0.17-0.30 mm), a 20 um film over a
     millimetre of slurry: with the rigid substrate on, every wet windward patch is thick, and only because of the
     slurry; with it off the liquid layer is the film alone and every such patch is thin."""
     assert body.MeltSettings().rigid_substrate
-    seen = {}
+    seen, conjugate = {}, {}
     for on in (True, False):
         b = melting_body(layered_mesh, rigid_substrate=on)
         sim, a = girin_state(b)
@@ -684,8 +709,17 @@
         h_p = deep_stage_args(b, delta_m)[4]
         b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
         seen[on] = b
+        conjugate[on] = delta_m
     on, off = seen[True], seen[False]
-    assert on.last_melt["thick_branch_fraction"] == 1.0 and off.last_melt["thick_branch_fraction"] == 0.0
+    assert on.last_melt["thick_branch_fraction"] == 1.0
+    # off: the liquid layer is the film alone, so every wet windward patch whose film is shallower than its conjugate
+    # depth takes the thin branch -- all but the few where the runoff piled the film deeper (57 of 2,296, measured
+    # 2026-10-22). The fraction above is taken over the patches still wet after the spray, and since 2026-10-22, when
+    # Girin's We_s stopped holding the thin branch back, the film sprays away within the step and leaves it nothing to
+    # average; the spray result carries the classification instead.
+    wet = off.windward & (off.last_spray.branch >= 0)
+    shallow = wet & (off.last_b < conjugate[False])
+    assert shallow.sum() > 0.9 * wet.sum() and (off.last_spray.branch[shallow] == BRANCH_THIN).all()
     assert on.last_melt["slurry_thick_fraction"] == 1.0 and on.last_melt["rigid_thin_fraction"] == 0.0
     assert 0.9 < on.last_melt["nonrigid_depth_mean_mm"] < 1.1                 # the millimetre of slurry, read back
     assert on.last_melt["slurry_held_mass_kg"] == 0.0                         # Girin's closure everywhere: nothing held
```

- [ ] **Step 3: The CLI's default** in `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -273,7 +273,7 @@
     s, r, f = doc["settings"], doc["results"], doc["files"]
     # AA7075_scheil is the material of --melt on (amendment of 2026-10-07): Scheil's curve, 390 kJ/kg, 0.80 N/m
     assert s["melt"] == "on" and s["material"] == "AA7075_scheil" and s["removal"] == "girin" and s["runoff"] == "on" and s["prism_layers"] == 0
-    assert s["rarefied_shear"] == "slip" and s["we_critical"] == 4.62 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.80
+    assert s["rarefied_shear"] == "slip" and s["we_critical"] == 3.08 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.80
     assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01 and s["deep_runoff"] == "on" and s["seed"] == cli.DEFAULT_SEED
     assert "deep_runoff_mass_kg" in rows[0] and "deep_mass_kg" in rows[0] and r["deep_runoff_mass_kg"] >= 0.0
     assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
```

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py -q`, then
  `"$PY" -m pytest tests/test_reentry_model_melting.py tests/test_reentry_model_cli.py -q -k
  "film_spraying_death_and_balances or films_shallower_than_the_floors or melting_run_writes_columns_files_and_json"`.
  Expected: 7 failed, 12 passed in the spray tests --
  - `AttributeError: module 'reentry_model.spray' has no attribute 'viscous_film'` (5: the four damping tests' and the
    regime-2 step-profile test);
  - `AttributeError: 'SprayResult' object has no attribute 'outside'` (the Stokes-depth reject);
  - `assert (np.False_)` (the thin branch below 3.08 still gated) and `assert np.float64(0.00024000000000000003) ==
    0.0` (the Kelvin–Helmholtz melt's thin mode not damped);

  and 3 failed in the others: `assert (22 == 24)` (the source table's width), `AttributeError: 'SprayResult' object
  has no attribute 'outside'` (the body's count of the reject), and the CLI settings' `we_critical` (4.62, not 3.08).

- [ ] **Step 5: The constants** in `reentry_model/dispersion.py`:

```diff
--- a/reentry_model/dispersion.py
+++ b/reentry_model/dispersion.py
@@ -19,8 +19,8 @@
 from . import DATA_DIR
 
 TABLE_PATH = os.path.join(DATA_DIR, "girin_dispersion.json")
-WE_CRITICAL_THEORY = 3.08          # Girin's critical surface Weber number
-WE_CRITICAL_PRACTICAL = 4.62       # the value he recommends in practice (GI >= 0.6)
+WE_CRITICAL_THEORY = 3.08          # Girin's critical surface Weber number; the model's default (fact 161)
+WE_CRITICAL_PRACTICAL = 4.62       # the value he recommends in practice (GI >= 0.6) for mm meteoroids; girin_case only
 DELTA_GRID = np.linspace(0.02, 5.0, 2500)
 WE_GRID = np.logspace(np.log10(3.0), 4.0, 200)
 UNSTABLE_TOL = 1e-6                # Im Omega above which a root counts as unstable
```

- [ ] **Step 6: The gates and the viscous film** in `reentry_model/spray.py`:

```diff
--- a/reentry_model/spray.py
+++ b/reentry_model/spray.py
@@ -59,14 +59,20 @@
       not distinguishable there; no history column is added; and `thick_branch_fraction` is `film.lubrication`'s own
       thick flag -- isfinite(delta_m) & (layer > delta_m), the whole deep set -- so it would count a regime-2 patch as
       thick. That column measures how deep the liquid is, not which instability took the patch.
-Two gates apply to every branch (2026-09-24):
+Two gates (2026-09-24; the first on the thick branch only since 2026-10-22):
   * Girin's (2017) stability criterion We_s > We_cr, with We_s = rho_l V_s^2 min(delta_m, layer) / Sigma -- his
     closed-form critical angle phi_cr evaluated against the model's own local flow, which holds on an eroded body where
     a single angle does not. On the thick branch this is his We_s = rho_l V_s^2 delta_m / Sigma unchanged; on a film
     shallower than the conjugate layer it is measured over the film's own depth. The thin branch had no threshold
     before: Girin & Kopyt's inviscid side mode has none of its own (it is unstable at every wavelength, their
     "unlimiting potential intensity of high-frequency disturbances"), so it stripped film at any angle however small
-    the shear -- 2.2 g inside 10 deg on the 100 mm flight, from a wave 6 x wider than its own facet.
+    the shear -- 2.2 g inside 10 deg on the 100 mm flight, from a wave 6 x wider than its own facet. Since 2026-10-22
+    the criterion gates the thick branch only: Girin & Kopyt's Eq. (5), solved, has no critical Weber number, and the
+    threshold it has once the wavelength is bounded is the wave fitting in the molten region, which `wave_fits` tests.
+    Their three branches (thin, rarefied, regime 2) take `viscous_film` instead: the dominant wave's viscous damping
+    on a film of the non-rigid depth over a rigid floor -- Lamb's interior term plus the floor's Stokes layer -- must
+    not outpace its growth 1/tau*, and a film shallower than that Stokes layer is outside the inviscid mode's validity
+    (recorded, not sprayed).
   * `wave_fits`: one whole wavelength must fit inside the *contiguous molten region* the patch belongs to, not on the
     patch (a mesh facet is not a physical boundary and the melt is continuous across it).
 KNOWN LIMITATION, deferred: Girin & Kopyt derive the side mode for the side surface (their sin Theta -> 0,
@@ -103,8 +109,10 @@
 BIN_EDGES = np.logspace(np.log10(R_MIN), np.log10(R_MAX), N_BINS + 1)
 WE_BREAKUP = 12.0                  # Pilch-Erdman: secondary breakup expected above this droplet Weber number
 CAPILLARY_TAU = 4.0 * np.pi / (2.0 * np.pi) ** 1.5      # 0.798: tau* = 2 x 2 pi / omega_cap(lambda*)
+# the last two: viscous_film's floor depth and gamma tau* (nan off Girin & Kopyt's branches; 2026-10-22)
 SOURCE_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "theta_deg", "x_m", "y_m", "z_m", "closure", "branch", "b_m",
-                  "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup"]
+                  "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup",
+                  "h_floor_m", "gamma_tau"]
 BRANCH_THICK, BRANCH_THIN, BRANCH_RAREFIED, BRANCH_RT = 0, 1, 2, 3      # Girin thick / thin-film with the edge state /
 #                                                                        thin-film with the freestream / front-surface Rayleigh-Taylor
 BRANCH_REGIME2 = 4                 # deep film whose melt is the more viscous medium: Girin & Kopyt's step profile again
@@ -193,6 +201,40 @@
     return lam, tau
 
 
+def viscous_film(lam, tau, h, liquid):
+    """The thin mode's wave on a viscous film of depth h over a rigid floor (2026-10-22): (gamma_tau, outside).
+
+    gamma is the amplitude damping rate of a capillary wave of wavenumber k = 2 pi / lambda* on a layer of depth h,
+        gamma = 2 nu k^2 + k sqrt(nu omega / 2) / sinh(2 k h),   omega^2 = (Sigma k^3 / rho_l) tanh(k h),
+    Lamb's interior damping plus the bottom (Stokes) boundary layer's of Hunt (1952), as Mei gives them in The Applied
+    Dynamics of Ocean Surface Waves; nu = mu_l / rho_l. The patch releases nothing where gamma_tau = gamma tau* >= 1:
+    viscosity takes the wave down at least as fast as it grows. On deep liquid the floor term vanishes and this is
+    Lamb's 2 nu k^2 tau* alone, which binds only below a wavelength of about 3.5 um for liquid aluminium; Girin & Kopyt's
+    own dissipation limit, their Eq. (9), carries no liquid viscosity and is a third of lambda* everywhere, so it never
+    binds.
+
+    `outside` marks the patches where h < sqrt(2 nu / omega): the floor's boundary layer reaches through the whole
+    liquid and Girin & Kopyt's inviscid mode does not apply. That is a validity limit, not a finding that no droplets
+    leave -- a sheared viscous film has long-wave instabilities of its own -- so the caller records the patches and
+    the film they hold rather than dropping them silently. The factor 2 under the root is a convention (the Stokes
+    layer's e-folding depth is sqrt(2 nu / omega)).
+
+    h is the depth to the rigid floor, which is where the liquid stops: the non-rigid depth the regime test uses, not
+    the film thickness b, since the thin branch exists because the rigid core lies within delta_m of the surface.
+    OPEN ASSUMPTION: mushy material more than half liquid counts as non-rigid there but is far more viscous than the
+    liquid, so the effective floor may lie inside that layer, and the damping would then be larger than computed."""
+    nu = liquid.mu / liquid.rho
+    lam, tau, h = (np.asarray(x, dtype=float) for x in (lam, tau, h))
+    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
+        k = 2.0 * np.pi / lam
+        kh = k * h
+        omega = np.sqrt(liquid.sigma * k ** 3 / liquid.rho * np.tanh(kh))
+        inv_sinh = np.where(kh > 0.0, 2.0 * np.exp(-2.0 * kh) / -np.expm1(-4.0 * kh), np.inf)   # 1 / sinh(2 k h)
+        gamma = 2.0 * nu * k ** 2 + k * np.sqrt(nu * omega / 2.0) * inv_sinh
+        outside = ~(h >= np.sqrt(2.0 * nu / omega))
+    return gamma * tau, outside
+
+
 def wave_fits(lam, extent):
     """True where one whole wavelength fits inside the contiguous molten region the patch belongs to.
 
@@ -222,10 +264,13 @@
     v_s: np.ndarray
     growth: np.ndarray = None   # s, the selected mode's growth time per patch (nan where no mode): what decides
     #                             whether the shear instability outruns the Rayleigh-Taylor mode on the same liquid
+    h_floor: np.ndarray = None  # m, the depth to the rigid floor viscous_film used (nan off Girin & Kopyt's branches)
+    gamma_tau: np.ndarray = None    # viscous_film's gamma tau* (nan off Girin & Kopyt's branches)
+    outside: np.ndarray = None  # bool, Girin & Kopyt's branches whose film is shallower than the floor's Stokes layer
 
 
 class SprayModel:
-    def __init__(self, liquid, k_r=K_R, k_t=K_T, we_critical=dispersion.WE_CRITICAL_PRACTICAL, table=None,
+    def __init__(self, liquid, k_r=K_R, k_t=K_T, we_critical=dispersion.WE_CRITICAL_THEORY, table=None,
                  rt_spray=True):
         if k_r <= 0.0 or k_t <= 0.0 or we_critical <= 0.0:
             raise ValueError("k_r, k_t and the critical Weber number must be > 0")
@@ -258,6 +303,8 @@
         regime = layer if regime_layer is None else np.maximum(np.asarray(regime_layer, dtype=float), b)
         branch = np.full(n, -1)
         growth = np.full(n, np.nan)          # the mode's growth time per patch [s], for comparison with other modes
+        h_floor, gamma_tau = np.full(n, np.nan), np.full(n, np.nan)
+        outside = np.zeros(n, dtype=bool)
         r = np.full(n, np.nan)
         mdot = np.zeros(n)
         has_film = windward & (b >= B_MIN)
@@ -285,13 +332,18 @@
         held = np.zeros(n, dtype=bool) if on_slurry is None else has_film & ~np.isfinite(delta_m) & np.asarray(on_slurry, dtype=bool)
         rarefied = has_film & ~deep & free_molecular & ~held              # no edge state exists: the freestream drives the mode
         thin = has_film & ~deep & ~rarefied & ~held
-        # Girin's (2017) stability criterion, now on *every* branch. His phi_cr comes from We_s > We_cr with the sphere's
-        # own edge solution substituted in; evaluated locally instead of through that closed form it holds on an eroded
+        # Girin's (2017) stability criterion, on the thick branch only (2026-10-22). His phi_cr comes from We_s > We_cr
+        # with the sphere's own edge solution substituted in; evaluated locally instead of through that closed form it holds on an eroded
         # body too, which a single critical angle does not. The shear acts over the liquid it can reach, min(delta_m,
         # layer), so the thick branch is his We_s = rho_l V_s^2 delta_m / Sigma unchanged while a film shallower than the
         # conjugate layer is measured on its own depth. Below We_cr the surface is stable and nothing leaves. The thin
         # branch had no stability threshold before this (Girin & Kopyt's inviscid side mode has none of its own, being
         # unstable at every wavelength), so it stripped film at any angle however small the shear (added 2026-09-24).
+        # Removed from Girin & Kopyt's branches (thin, rarefied, regime 2) on 2026-10-22: their Eq. (5) has no critical
+        # Weber number (the growth rate is positive at every We), and the threshold it does have once the wavelength is
+        # bounded is a wave that fits inside the molten region -- wave_fits on lambda*, which also stops the stagnation
+        # point, where lambda* grows without bound as the shear vanishes. They are gated instead by wave_fits and by
+        # viscous_film. We_s is still computed on them and written to the source table.
         # A patch the regime test makes deep is measured on delta_m whatever its liquid depth -- a film on slurry deeper
         # than delta_m included (2026-10-06); everywhere else the expression is the one above, bit for bit.
         L = None if extent is None else np.asarray(extent, dtype=float)
@@ -325,7 +377,9 @@
         if regime2.any():
             m = regime2
             lam, tau = thin_film_mode(flow.mach_e[m], (flow.rho_e * flow.u_eff ** 2)[m], liq)
-            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & supercritical[m]
+            gt, out = viscous_film(lam, tau, regime[m], liq)
+            h_floor[m], gamma_tau[m], outside[m] = regime[m], gt, out
+            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & (gt < 1.0) & ~out
                   & wave_fits(lam, None if L is None else L[m]))
             rr = np.where(ok, self.k_r * lam, np.nan)
             with np.errstate(divide="ignore", invalid="ignore"):
@@ -338,7 +392,9 @@
             if not mask.any():
                 continue
             lam, tau = thin_film_mode(mach[mask], flux[mask], liq)
-            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & supercritical[mask]
+            gt, out = viscous_film(lam, tau, regime[mask], liq)
+            h_floor[mask], gamma_tau[mask], outside[mask] = regime[mask], gt, out
+            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & (gt < 1.0) & ~out
                   & wave_fits(lam, None if L is None else L[mask]))
             rr = np.where(ok, lam / 4.0, np.nan)
             rate = np.where(ok, liq.rho * np.minimum(b[mask], lam / 8.0) / np.where(ok, tau, 1.0), 0.0)
@@ -389,7 +445,10 @@
         with np.errstate(divide="ignore", invalid="ignore"):
             dn = np.where(dm > 0.0, dm / (4.0 / 3.0 * np.pi * liq.rho * np.where(dm > 0.0, r, 1.0) ** 3), 0.0)
         rt_active, _, _ = rayleigh_taylor(w_n, layer, liq)     # the unbounded criterion, still reported for continuity
-        return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s, growth)
+        gk = np.isin(branch, (BRANCH_THIN, BRANCH_RAREFIED, BRANCH_REGIME2))    # the Rayleigh-Taylor mode may have taken some
+        h_floor, gamma_tau = np.where(gk, h_floor, np.nan), np.where(gk, gamma_tau, np.nan)
+        return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s, growth,
+                           h_floor, gamma_tau, outside & gk)
 
 
 def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state):
@@ -403,7 +462,7 @@
     rows = np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
                             centroids[k, 0], centroids[k, 1], centroids[k, 2], flow.closure[k], res.branch[k], b[k],
                             res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1], t_hat[k, 2],
-                            we_d, oh, (we_d > WE_BREAKUP).astype(float)])
+                            we_d, oh, (we_d > WE_BREAKUP).astype(float), res.h_floor[k], res.gamma_tau[k]])
     return rows.tolist()
 
 
```

- [ ] **Step 7: The reject recorded** in `reentry_model/body.py`:

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -269,6 +269,8 @@
         self.deep_surfaced_mass = 0.0                                       # kg of deep liquid that became film from the top
         self.cascade_mass = 0.0                                             # kg the molten cascade fed (2026-10-03)
         self.cascade_capped_steps = 0                                       # steps whose cascade stopped at MAX_CASCADE_PASSES
+        self.thin_invalid_patch_steps = 0       # Girin & Kopyt patches outside the inviscid mode's validity (2026-10-22)
+        self.thin_invalid_mass_max = 0.0        # kg, the most film held on such patches in one step
         self.n_released = 0.0
         self.source_rows = []
         self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
@@ -959,6 +961,12 @@
         res = self.spray.evaluate(flow, state, b, delta_m, v_s, self.windward, dt, areas, self.m_f,
                                   self.transverse_radius, b_layer=layer, extent=extent, deceleration_n=w_n,
                                   regime_layer=None if nonrigid is None else regime_layer, on_slurry=on_slurry)
+        # Girin & Kopyt's patches whose film is shallower than the floor's Stokes layer: outside the inviscid mode's
+        # validity, so nothing leaves them by it -- recorded with the film they hold, not dropped silently (2026-10-22)
+        invalid = res.outside & self.windward
+        invalid_mass = float(self.m_f[invalid].sum())
+        self.thin_invalid_patch_steps += int(invalid.sum())
+        self.thin_invalid_mass_max = max(self.thin_invalid_mass_max, invalid_mass)
         released = float(res.dm.sum())
         if released > 0.0:
             if self.spray_onset is None:
@@ -1039,7 +1047,8 @@
             "thick_branch_fraction": float(thick[self.windward & (self.m_f > 0.0)].mean())
             if (self.windward & (self.m_f > 0.0)).any() else float("nan"),
             "nonrigid_depth_mean_mm": nonrigid_mean, "slurry_thick_fraction": slurry_thick,
-            "rigid_thin_fraction": rigid_thin, "slurry_held_mass_kg": held_mass, **depths})
+            "rigid_thin_fraction": rigid_thin, "slurry_held_mass_kg": held_mass,
+            "thin_invalid_patches": int(invalid.sum()), "thin_invalid_mass_kg": invalid_mass, **depths})
         return released
 
     def _kill(self, dead):
@@ -1620,6 +1629,8 @@
                 "slurry_thick_fraction": lm.get("slurry_thick_fraction", float("nan")),
                 "rigid_thin_fraction": lm.get("rigid_thin_fraction", float("nan")),
                 "slurry_held_mass_kg": lm.get("slurry_held_mass_kg", float("nan")),
+                "thin_invalid_patches": float(lm.get("thin_invalid_patches", 0)),
+                "thin_invalid_mass_kg": lm.get("thin_invalid_mass_kg", 0.0),
                 "deep_liquid_kg": lm.get("deep_liquid", 0.0), "deep_mass_kg": float(self.m_d.sum()),
                 "deep_runoff_mass_kg": self.deep_runoff_mass, "deep_surfaced_mass_kg": self.deep_surfaced_mass,
                 "deep_blob_fraction": self.deep_blob_fraction(),
```

- [ ] **Step 8: Its history columns and results** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -29,7 +29,8 @@
                 "unapplied_load_J", "film_blob_fraction", "rt_mass_fraction", "rt_wavelength_over_nose", "rt_growth_ms", "spray_growth_ms", "rt_region_mm", "rt_bounded_fraction", "rt_bounded_growth_ms", "molten_depth_max_mm", "molten_depth_mean_mm",
                 "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements", "deep_liquid_kg", "deep_mass_kg",
                 "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction", "cascade_passes", "cascade_mass_kg",
-                "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg"]
+                "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg",
+                "thin_invalid_patches", "thin_invalid_mass_kg"]
 SKIN_COLUMNS = ["skin_count", "skin_mass_kg", "skin_thickness_min_mm", "skin_thickness_median_mm", "skin_drawn_mass_kg",
                 "skin_handover_mass_kg", "skin_short_events", "interface_mismatch_J", "interface_repeats",
                 "liquid_depth_mean_um", "liquid_depth_max_um", "nonrigid_depth_skin_mean_um", "skin_wall_s",
@@ -169,6 +170,7 @@
                 "deep_mass_kg": float(body.m_d.sum()),
                 "cascade_mass_kg": body.cascade_mass, "cascade_passes_max": int(c["cascade_passes"].max()),
                 "cascade_capped_steps": body.cascade_capped_steps,
+                "thin_invalid_patch_steps": body.thin_invalid_patch_steps, "thin_invalid_mass_max_kg": body.thin_invalid_mass_max,
                 "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                 "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                 "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
```

- [ ] **Step 9: The flag's default and help** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -337,7 +337,9 @@
     me.add_argument("--kn-body-shock", type=float, default=surface_flow.KN_BODY_SHOCK,
                     help="body Knudsen number below which a distinct bow shock is assumed and the shock-layer construction is used "
                          "(default %(default)s); above it the flow is treated as merged and every melt closure is flagged")
-    me.add_argument("--we-critical", type=float, default=dispersion.WE_CRITICAL_PRACTICAL, help="critical surface Weber number (default %(default)s)")
+    me.add_argument("--we-critical", type=float, default=dispersion.WE_CRITICAL_THEORY,
+                    help="critical surface Weber number of Girin's thick branch (default %(default)s); Girin & Kopyt's thin, "
+                         "rarefied and regime-2 branches are gated by the wave fit and the viscous film instead")
     me.add_argument("--kr", type=float, default=spray.K_R, help="droplet radius / wavelength (default %(default)s)")
     me.add_argument("--kt", type=float, default=spray.K_T, help="release period / growth time (default %(default)s)")
     me.add_argument("--prism-layers", type=int, default=None, help="prism layers under the surface (default 4 with --melt on, 0 otherwise)")
```

- [ ] **Step 10: The sensitivity variant** in `analysis/melt_sensitivity.py`:

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -2,7 +2,7 @@
 """Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).
 
     "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
-        [--variants base,seed1,layers2,layers6,dt025,dtcoff,bridged,norunoff,nodeep,nocascade,norigid,we308,kr-30,kr+30,
+        [--variants base,seed1,layers2,layers6,dt025,dtcoff,bridged,norunoff,nodeep,nocascade,norigid,we462,kr-30,kr+30,
                     kt-30,kt+30,AA7075,range,sizeinitial,gammapm14,knbody003,fenicsx]
 
 Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
@@ -36,7 +36,7 @@
 VARIANTS = {
     "base": [], "seed1": ["--seed", "1"], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"], "dtcoff": ["--dt-continuum", "off"],
     "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "nodeep": ["--deep-runoff", "off"],
-    "nocascade": ["--molten-cascade", "off"], "norigid": ["--rigid-substrate", "off"], "we308": ["--we-critical", "3.08"],
+    "nocascade": ["--molten-cascade", "off"], "norigid": ["--rigid-substrate", "off"], "we462": ["--we-critical", "4.62"],
     "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
     "AA7075": ["--material", "AA7075"], "range": ["--material", "AA7075_range"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
     "knbody003": ["--kn-body-shock", "0.003"], "fenicsx": ["--thermal-solver", "fenicsx"],
```

- [ ] **Step 11: Run them to see them pass.** Expected: `tests/test_reentry_model_spray.py` 19 passed (15 before),
  `tests/test_reentry_model_melting.py` 26 (25 before), `tests/test_reentry_model_cli.py` 66,
  `tests/test_reentry_model_coupled.py` 9 (101 for the last three, 263 s); the unit tier 522 passed, 1 skipped.

- [ ] **Step 12: The short case against the branch.** Copy that branch's tested copy
  (`.claude/worktrees/we-critical-308/prototype/work-2026-10-22-thin-gate/code/`, without its `reentry_model_output/`)
  into `measure/task1/wc308_code/` and run the case on both and compare: `./run_case.sh <code> runs line_short`,
  `./run_case.sh .../wc308_code runs wc308_short`, `"$PY" compare_runs.py runs line_short wc308_short` and
  `"$PY" branches_short.py runs line_short` (in `measure/task1/cmp/`). Expected: no shared history column, source field
  or result differing; the branch counts above.

- [ ] **Step 13: Commit.**

```bash
git add reentry_model/spray.py reentry_model/dispersion.py reentry_model/cli.py reentry_model/body.py reentry_model/coupled.py analysis/melt_sensitivity.py tests/test_reentry_model_spray.py tests/test_reentry_model_melting.py tests/test_reentry_model_cli.py
git commit -m "Spray: the element model's two amendments of 2026-10-22 brought across from we-critical-308 (fa02a4c, a54520c) -- Girin's critical surface Weber number 3.08 by default and on the thick branch only, Girin & Kopyt's thin, rarefied and regime-2 branches gated by the wave fit and a viscous film on the rigid floor, the validity reject recorded; on the CLI test's short melting case the element model sprays as that branch's does to the last bit (sub-plan 22 Task 1a)"
```

---

### Task 1b: The droplet's size cap per crest

**Files:**
- Modify: `reentry_model/spray.py` (`SprayModel.evaluate(..., cap_depth=None)`; `SprayResult.wavelength`; the module
  docstring's size caps)
- Test: `tests/test_reentry_model_spray.py` (four tests)

**Interfaces:**
- Consumes: Task 1a's `SprayModel.evaluate` and its modes.
- Produces: `SprayModel.evaluate(..., cap_depth=None)`, `cap_depth` (n,) [m] the depth the release can take from each
  patch, or None for the film-mass cap; `SprayResult.wavelength` (n,) [m], the wavelength of the mode that sprays each
  patch (nan where none). What Task 1c passes: the admitted run's depth at the macro step's start, the liquid depth
  with `--slurry-spray off` and the flowing depth with it on.

Asha decided on 2026-10-24 (the amendment of that date above; spec §8 as amended) that a droplet holds no more liquid
than one crest of its own wave can give, (4/3)π r³ ≤ λ² h, the quarter-body-radius cap unchanged, and that the shared
spray physics takes the cap's depth as its own argument and keeps the film-mass cap without it. `evaluate` therefore
records each branch's wavelength as it selects the mode — the thick branch's λ_f = 2π δ_m / Δ_f, Girin & Kopyt's λ* on
the thin, rarefied and regime-2 branches, the Rayleigh–Taylor wavelength where that mode takes the patch — and, given
`cap_depth`, caps the radius at min((3 λ² h / (4π))^⅓, R/4); a patch with no depth releases nothing. Without it the
radius is capped as before, by the same expression, so the element model is unchanged.

**Where it binds.** The thick branch's droplet is k_r λ_f (k_r = 0.17) and the crest's radius is under it only where h
is under 4π k_r³ λ / 3, about 0.021 λ (18 µm at λ_f = 0.89 mm, the test case's); the thin branch's is λ*/4, bound under
about 0.065 λ* (30 µm at 0.46 mm); the Rayleigh–Taylor droplet is its wavelength, bound under 4π λ / 3, about 4.2 λ,
so always. The outline's check: 4.57 mm at λ = 2 cm and h = 1 mm (the test), against the element model's thin-gate
flight, where the film cap held its Rayleigh–Taylor droplets to 1.8–2.6 mm (10th to 90th percentile of 620 rows; 2.26 mm
the median by mass, 3.15 mm the largest) on films piled 12.7 mm deep at the nose. On the 100 mm sphere's skins with
3.4 mm flowing under 103 m/s² the crest gives 3.05, 4.56 and 6.72 mm on 0.3, 1.0 and 3.3 mm of liquid (fact 171).
For the size bins (1 µm to 10 mm): where the Rayleigh–Taylor depth criterion reads the cap's own depth h it admits only
λ* < 2π h, so r < (3π)^⅓ h, about 2.1 h — inside the top bin while h is under 4.7 mm, and below the 100 mm sphere's
quarter-radius cap of 12.5 mm; a 4.6 mm droplet holds about 1 g, so these move the count by almost nothing and the mass
in the top bins by what the mode releases (fact 170).

**Measured** (2026-10-28, `measure/task1/cmp/`): the short case of Task 1a, rerun with this task's code
(`linecap_short`), is the same as without it to the last bit in all 115 history columns, the source table and the
results.

**Rulings** (2026-10-28):
1. *The argument is `cap_depth`, a depth in metres per patch*, the outline's name. Why: it is what the caller knows —
   the admitted run's depth — and keeps the mass, which needs an area, out of the cap. Cost: none.
2. *The Rayleigh–Taylor λ is the wavelength the mode admits*, `rayleigh_taylor_bounded`'s 2π / max(π / L, k_max): λ*
   where the pool admits the fastest mode, twice the pool's extent where it does not, the wavelength the droplet's radius
   already is. Why: the cap is on the wave that sprays. Cost if wrong: only in pools narrower than λ*/2, where the cap
   would use λ* and be larger.
3. *A depth of nothing releases nothing*, the radius nan. Why: one crest gives nothing on no depth, and a radius of zero
   would make the droplet count infinite. Cost: none.
4. *The result carries the wavelength.* Why: Task 2's records and the tests read it; the element model writes nothing
   new. Cost: none.

- [ ] **Step 1: The cap's tests** — append to `tests/test_reentry_model_spray.py`:

```diff
--- a/tests/test_reentry_model_spray.py
+++ b/tests/test_reentry_model_spray.py
@@ -439,3 +439,117 @@
         # (the floor a millimetre down keeps the rarefied 7.5 mm wave inside the inviscid mode's validity, viscous_film)
         assert res.branch.tolist() == [-1, code]
         assert res.dm[0] == 0.0 and not res.unstable[0] and res.dm[1] > 0.0
+
+
+# -- the droplet-size cap per crest (2026-10-24, Asha; sub-plan 22 Task 1a) ------------------------------------------
+def crest(lam, h):
+    return (3.0 * lam ** 2 * h / (4.0 * np.pi)) ** (1.0 / 3.0)
+
+
+def spraying_cases(cap_depth=None, areas=None):
+    """One patch of each branch the cap reaches, each evaluated as its own body: (result, the mode's own radius per
+    wavelength) for the thick branch (Girin's closure, 1 mm of liquid), the thin branch (20 um on a floor 20 um down),
+    the rarefied branch (20 um on a floor 1 mm down, where its 7.5 mm wave stays inside the inviscid mode's validity) and
+    the front-surface Rayleigh-Taylor mode (5 mm of liquid at the nose, its lambda* 2 cm at W = 106 m/s2)."""
+    model = spray.SprayModel(LIQ)
+    areas = np.full(1, 1e-5) if areas is None else np.asarray(areas, dtype=float)
+    out = {}
+    h = None if cap_depth is None else np.full(1, cap_depth)
+    flow = FakeFlow(1)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.full(1, 1e-3)
+    out["thick"] = (model.evaluate(flow, FakeState(), b, delta_m, np.full(1, 5.0), np.ones(1, bool), 0.5, areas,
+                                   b * LIQ.rho * areas, cap_depth=h), 0.17)
+    b = np.full(1, 2e-5)
+    out["thin"] = (model.evaluate(flow, FakeState(), b, delta_m, np.full(1, 20.0), np.ones(1, bool), 0.5, areas,
+                                  b * LIQ.rho * areas, cap_depth=h), 0.25)
+    flow_fm = FakeFlow(1, closure=sf.CLOSURE_COUETTE, branch=sf.BRANCH_FREE_MOLECULAR)
+    dm_fm, _ = spray.melt_layer(flow_fm, LIQ)
+    out["rarefied"] = (model.evaluate(flow_fm, FakeState(), b, dm_fm, np.full(1, 20.0), np.ones(1, bool), 0.5, areas,
+                                      b * LIQ.rho * areas, b_layer=np.full(1, 1e-3), regime_layer=np.full(1, 1e-3),
+                                      cap_depth=h), 0.25)
+    flow_rt = FakeFlow(1, closure=sf.CLOSURE_COUETTE, u=5.0)       # at the nose: the thin lambda* is metres, no fit
+    dm_rt, _ = spray.melt_layer(flow_rt, LIQ)
+    w_n = np.full(1, 3.0 * LIQ.sigma * (2.0 * np.pi / 0.02) ** 2 / LIQ.rho)     # lambda* = 2 pi (3 Sigma / (W rho))^1/2
+    b = np.full(1, 5e-3)
+    out["Rayleigh-Taylor"] = (model.evaluate(flow_rt, FakeState(), b, dm_rt, np.full(1, 0.1), np.ones(1, bool), 0.5,
+                                             areas, b * LIQ.rho * areas, extent=np.full(1, 0.05), deceleration_n=w_n,
+                                             cap_depth=h), 1.0)
+    return out
+
+
+def test_a_droplet_holds_no_more_than_one_crest_of_its_own_wave_on_the_depth_given():
+    """With a cap depth h (the layered skins, 2026-10-24) every branch's droplet is capped by one crest of the wave that
+    sprays it, (4/3) pi r^3 <= lambda^2 h, and by a quarter of the body radius, instead of by the film on the patch:
+    lambda is the thick branch's 2 pi delta_m / Delta_f, Girin & Kopyt's lambda* on the thin and rarefied branches and the
+    Rayleigh-Taylor lambda*. The outline's check: at lambda = 2 cm and h = 1 mm the Rayleigh-Taylor droplet is 4.6 mm."""
+    codes = {"thick": spray.BRANCH_THICK, "thin": spray.BRANCH_THIN, "rarefied": spray.BRANCH_RAREFIED,
+             "Rayleigh-Taylor": spray.BRANCH_RT}
+    for h in (1e-3, 5e-6):
+        for name, (res, k) in spraying_cases(cap_depth=h).items():
+            lam = res.wavelength[0]
+            assert res.branch[0] == codes[name] and res.dm[0] > 0.0 and np.isfinite(lam)
+            assert res.r[0] == pytest.approx(min(k * lam, crest(lam, h), 0.0125), rel=1e-12)
+    res, _ = spraying_cases(cap_depth=1e-3)["thick"]
+    flow = FakeFlow(1)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    d_f, _, _ = dispersion.DispersionTable()(res.we_s)
+    assert res.wavelength[0] == pytest.approx(2.0 * np.pi * delta_m[0] / d_f[0], rel=1e-12)     # the thick lambda
+    lam_t, _ = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
+    assert spraying_cases(cap_depth=1e-3)["thin"][0].wavelength[0] == pytest.approx(lam_t[0], rel=1e-12)
+    res, _ = spraying_cases(cap_depth=1e-3)["Rayleigh-Taylor"]
+    assert res.wavelength[0] == pytest.approx(0.02, rel=1e-12)
+    assert res.r[0] == pytest.approx(4.57e-3, rel=1e-3)                       # the outline's 4.6 mm
+    film = spraying_cases()["Rayleigh-Taylor"][0]
+    assert film.r[0] == pytest.approx((3.0 * 5e-3 * 1e-5 / (4.0 * np.pi)) ** (1.0 / 3.0), rel=1e-12)   # 2.3 mm on 5 mm
+    # where it binds: on the thick and thin branches only on layers shallower than a few per cent of the wavelength,
+    # on the Rayleigh-Taylor branch on any layer shallower than about four wavelengths
+    for name, (res, k) in spraying_cases(cap_depth=1e-3).items():
+        assert (res.r[0] < k * res.wavelength[0]) == (name == "Rayleigh-Taylor")
+    for name, (res, k) in spraying_cases(cap_depth=5e-6).items():
+        assert res.r[0] < k * res.wavelength[0]
+
+
+def test_the_crest_cap_does_not_see_the_patch_s_area():
+    """lambda^2 h is a crest's volume whatever the mesh facet: four times the area changes no capped radius and no rate,
+    where the film-mass cap grows with the area as its cube root."""
+    for name in ("thick", "thin", "rarefied", "Rayleigh-Taylor"):
+        small, _ = spraying_cases(cap_depth=2e-5, areas=[1e-5])[name]
+        large, _ = spraying_cases(cap_depth=2e-5, areas=[4e-5])[name]
+        assert large.r[0] == small.r[0] and large.mdot[0] == small.mdot[0] and large.dm[0] == pytest.approx(4.0 * small.dm[0])
+    small, _ = spraying_cases(areas=[1e-5])["Rayleigh-Taylor"]
+    large, _ = spraying_cases(areas=[4e-5])["Rayleigh-Taylor"]
+    assert large.r[0] == pytest.approx(4.0 ** (1.0 / 3.0) * small.r[0], rel=1e-12)
+
+
+def test_a_cap_depth_of_nothing_releases_nothing():
+    """One crest gives nothing on no depth: where the release can take nothing, no droplet forms and nothing leaves,
+    though the patch keeps its branch."""
+    flow = FakeFlow(2)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.full(2, 1e-3)
+    areas = np.full(2, 1e-5)
+    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.full(2, 5.0), np.ones(2, bool), 0.5, areas,
+                                         b * LIQ.rho * areas, cap_depth=np.array([0.0, 1e-4]))
+    assert res.branch.tolist() == [spray.BRANCH_THICK] * 2 and np.isfinite(res.wavelength).all()
+    assert res.mdot[0] == 0.0 and res.dm[0] == 0.0 and res.dn[0] == 0.0 and not res.unstable[0] and np.isnan(res.r[0])
+    assert res.mdot[1] > 0.0 and res.dm[1] > 0.0 and np.isfinite(res.r[1])
+
+
+def test_without_a_cap_depth_the_film_mass_caps_the_droplet_as_before():
+    """The element model keeps its cap (2026-10-24): with no cap depth every branch's radius is the mode's own, capped by
+    the film on the patch and a quarter of the body radius, to the last bit, and the rates and releases are untouched by
+    the wavelength the result now carries."""
+    for name, (res, k) in spraying_cases().items():
+        lam = res.wavelength[0]
+        r_mode = lam if name == "Rayleigh-Taylor" else (0.17 * lam if name == "thick" else lam / 4.0)
+        m_f = {"thick": 1e-3, "thin": 2e-5, "rarefied": 2e-5, "Rayleigh-Taylor": 5e-3}[name] * LIQ.rho * 1e-5
+        film = (3.0 * m_f / (4.0 * np.pi * LIQ.rho)) ** (1.0 / 3.0)
+        assert res.r[0] == np.minimum(r_mode, np.minimum(film, 0.25 * 0.05))          # exactly, the same arithmetic
+        assert res.dm[0] == min(res.mdot[0] * 1e-5 * 0.5, m_f)
+    flow = FakeFlow(3)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.array([1e-3, 2e-5, 0.0])
+    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.array([5.0, 20.0, 0.0]), np.ones(3, bool),
+                                         0.5, np.full(3, 1e-5), b * LIQ.rho * 1e-5)
+    assert np.isfinite(res.wavelength[:2]).all() and np.isnan(res.wavelength[2])        # no mode on the dry patch
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py -q -k "crest or
  cap_depth"`. Expected: 4 failed, `TypeError: SprayModel.evaluate() got an unexpected keyword argument 'cap_depth'`.

- [ ] **Step 3: The cap per crest** in `reentry_model/spray.py`:

```diff
--- a/reentry_model/spray.py
+++ b/reentry_model/spray.py
@@ -95,7 +95,21 @@
 nose by 1/cos(phi)).
 
 Release per patch and step: dm = min(mdot A dt, m_f) (the film after feed and runoff), dn = dm / (4/3 pi rho_l r^3);
-one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX]."""
+one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX].
+
+The droplet's size caps (2026-10-24, Asha). Every branch's radius is capped by a quarter of the body radius and by one
+of two volumes. Without `cap_depth` (the element model) it is the film on the patch, r <= (3 m_f / (4 pi rho_l))^1/3,
+which depends on the mesh facet's area: measured on the element model's 100 mm flight with the thin gate it bound on
+67 % of the thick branch's release rows (24 % of its mass), 16 % of the thin branch's (4.4 %) and set essentially every
+Rayleigh-Taylor droplet (about 2 mm). With `cap_depth` h (the layered skins, where there is no film account) a droplet
+holds no more liquid than one crest of its own wave can give, (4/3) pi r^3 <= lambda^2 h, lambda the wavelength of the
+mode that sprays the patch -- the thick branch's 2 pi delta_m / Delta_f, Girin & Kopyt's lambda* on the thin, rarefied
+and regime-2 branches, the admitted Rayleigh-Taylor wavelength -- and h the depth the release can take: lambda^2
+because a crest repeats every wavelength in both surface directions, so the patch's area does not enter. A depth of
+nothing gives no droplet, and nothing leaves. On the thick and thin branches it binds only on layers shallower than
+about 2 % (thick: k_r = 0.17) and 6.5 % (thin: lambda*/4) of the wavelength; on the Rayleigh-Taylor branch, whose
+droplet is its wavelength (centimetres), on any layer shallower than about four wavelengths, so always: at lambda = 2 cm
+and h = 1 mm it is 4.6 mm, against about 2 mm under the film cap."""
 from dataclasses import dataclass, field
 
 import numpy as np
@@ -267,6 +281,8 @@
     h_floor: np.ndarray = None  # m, the depth to the rigid floor viscous_film used (nan off Girin & Kopyt's branches)
     gamma_tau: np.ndarray = None    # viscous_film's gamma tau* (nan off Girin & Kopyt's branches)
     outside: np.ndarray = None  # bool, Girin & Kopyt's branches whose film is shallower than the floor's Stokes layer
+    wavelength: np.ndarray = None   # m, the wavelength of the mode that sprays each patch (nan where none): the
+    #                                 per-crest cap's lambda (2026-10-24)
 
 
 class SprayModel:
@@ -279,7 +295,7 @@
         self.table = table or dispersion.DispersionTable()
 
     def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, b_layer=None,
-                 extent=None, deceleration_n=None, regime_layer=None, on_slurry=None):
+                 extent=None, deceleration_n=None, regime_layer=None, on_slurry=None, cap_depth=None):
         """Per-patch instability and release for the film thickness b, melt-layer thickness delta_m and film surface
         velocity v_s (from film.lubrication), over the step dt; m_f is the film mass available, radius the body's.
 
@@ -296,7 +312,12 @@
         a film on slurry deeper than delta_m is deep melt and takes the thick branch. `on_slurry`, if given, flags the
         patches whose film rests on slurry; where there is no delta_m (off Girin's closure) such a film takes no shear
         mode at all, since the thin mode needs a rigid wall and no deep-melt mode exists there. The Rayleigh-Taylor
-        criteria, a liquid pool's, keep `b_layer` (amendment of 2026-10-06)."""
+        criteria, a liquid pool's, keep `b_layer` (amendment of 2026-10-06).
+
+        `cap_depth`, if given, is the depth [m] the release can take from each patch, and the droplet's radius is capped
+        by one crest of its own wave on it, (3 lambda^2 h / (4 pi))^1/3, instead of by the film mass m_f (the module
+        docstring; 2026-10-24); a patch with no such depth releases nothing. Without it the film-mass cap is kept, bit
+        for bit."""
         liq = self.liquid
         n = b.size
         layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
@@ -304,6 +325,7 @@
         branch = np.full(n, -1)
         growth = np.full(n, np.nan)          # the mode's growth time per patch [s], for comparison with other modes
         h_floor, gamma_tau = np.full(n, np.nan), np.full(n, np.nan)
+        wavelength = np.full(n, np.nan)      # the spraying mode's wavelength per patch [m], for the per-crest cap
         outside = np.zeros(n, dtype=bool)
         r = np.full(n, np.nan)
         mdot = np.zeros(n)
@@ -366,6 +388,7 @@
                 rate = np.where(unstable, liq.rho * np.pi * rr * rr / (lam * t_per), 0.0)
             r[thick], mdot[thick], branch[thick] = rr, rate, BRANCH_THICK
             growth[thick] = t_per
+            wavelength[thick] = lam
         # regime 2: the other half of the deep-film set -- Girin & Kopyt's step profile, in their deep-film limit.
         # Wavelength and growth time are their Eqs. (11) and (12), the same pair the thin branch uses, because the mode is
         # the same one; what differs is the supply. The thin branch multiplies by min(b, lambda*/8) because a film
@@ -386,6 +409,7 @@
                 rate = np.where(ok, liq.rho * np.pi * self.k_r ** 2 * lam / np.where(ok, tau, 1.0), 0.0)
             r[m], mdot[m], branch[m] = rr, rate, BRANCH_REGIME2
             growth[m] = np.where(ok, tau, np.nan)
+            wavelength[m] = np.where(ok, lam, np.nan)
         # thin: Girin & Kopyt 1994 with the edge state
         for mask, code, mach, flux in ((thin, BRANCH_THIN, flow.mach_e, flow.rho_e * flow.u_eff ** 2),
                                        (rarefied, BRANCH_RAREFIED, np.full(n, state.ma), np.full(n, state.freestream.rho * state.V ** 2))):
@@ -400,6 +424,7 @@
             rate = np.where(ok, liq.rho * np.minimum(b[mask], lam / 8.0) / np.where(ok, tau, 1.0), 0.0)
             r[mask], mdot[mask], branch[mask] = rr, rate, code
             growth[mask] = np.where(ok, tau, np.nan)
+            wavelength[mask] = np.where(ok, lam, np.nan)
         # (iv) the front-surface Rayleigh-Taylor mode, now *applied* rather than only reported.
         #
         # Girin & Kopyt assign the front surface -- where the flow influence is negligible, V_0 -> 0 and sin Theta -> 1 --
@@ -436,19 +461,28 @@
             mdot = np.where(faster, rate_rt, mdot)
             branch = np.where(faster, BRANCH_RT, branch)
             growth = np.where(faster, tau_rt, growth)
+            wavelength = np.where(faster, lam_rt, wavelength)
+        if cap_depth is not None:                # one crest gives nothing on no depth: nothing leaves (2026-10-24)
+            h_cap = np.asarray(cap_depth, dtype=float)
+            mdot = np.where(h_cap > 0.0, mdot, 0.0)
         unstable = mdot > 0.0
         dm = np.where(unstable, np.minimum(mdot * areas * dt, m_f), 0.0)
-        # a droplet is never larger than the film on its patch nor than a quarter of the body radius
-        # (the long near-critical waves Girin excludes; the edge dynamic pressure vanishes toward the equator)
+        # a droplet is never larger than a quarter of the body radius (the long near-critical waves Girin excludes; the
+        # edge dynamic pressure vanishes toward the equator) nor than the film on its patch -- or, given cap_depth, than
+        # one crest of its own wave on the depth the release can take, (4/3) pi r^3 <= lambda^2 h (2026-10-24)
         with np.errstate(invalid="ignore"):
-            r = np.minimum(r, np.minimum((3.0 * m_f / (4.0 * np.pi * liq.rho)) ** (1.0 / 3.0), 0.25 * radius))
+            if cap_depth is None:
+                r = np.minimum(r, np.minimum((3.0 * m_f / (4.0 * np.pi * liq.rho)) ** (1.0 / 3.0), 0.25 * radius))
+            else:
+                crest = (3.0 * wavelength ** 2 * h_cap / (4.0 * np.pi)) ** (1.0 / 3.0)
+                r = np.where(h_cap > 0.0, np.minimum(r, np.minimum(crest, 0.25 * radius)), np.nan)
         with np.errstate(divide="ignore", invalid="ignore"):
             dn = np.where(dm > 0.0, dm / (4.0 / 3.0 * np.pi * liq.rho * np.where(dm > 0.0, r, 1.0) ** 3), 0.0)
         rt_active, _, _ = rayleigh_taylor(w_n, layer, liq)     # the unbounded criterion, still reported for continuity
         gk = np.isin(branch, (BRANCH_THIN, BRANCH_RAREFIED, BRANCH_REGIME2))    # the Rayleigh-Taylor mode may have taken some
         h_floor, gamma_tau = np.where(gk, h_floor, np.nan), np.where(gk, gamma_tau, np.nan)
         return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s, growth,
-                           h_floor, gamma_tau, outside & gk)
+                           h_floor, gamma_tau, outside & gk, wavelength)
 
 
 def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state):
```

- [ ] **Step 4: Run them to see them pass.** Expected: `tests/test_reentry_model_spray.py` 23 passed; the unit tier
  526 passed, 1 skipped. Then rerun the short case (`./run_case.sh <code> runs linecap_short`, `"$PY" compare_runs.py
  runs line_short linecap_short`): nothing differs.

- [ ] **Step 5: Commit.**

```bash
git add reentry_model/spray.py tests/test_reentry_model_spray.py
git commit -m "Spray: the droplet's size cap per crest, (4/3) pi r^3 <= lambda^2 h with lambda the spraying mode's wavelength and h the depth the release can take, given as cap_depth; without it the film-mass cap, so the element model is unchanged to the last bit (spec 8, Asha 2026-10-24; sub-plan 22 Task 1b)"
```

---

### Task 1c: Spraying from each column's top cells

**Files:**
- Modify: `reentry_model/spray.py` (`SprayModel.evaluate(..., rt_layer=None)`), `reentry_model/skin.py` (`Cascade`,
  `cascade_mass`; `Flow.release` may give a `Cascade`; `FlowBook.released_slurry`; `flow_substeps`),
  `reentry_model/skinbody.py` (`SLURRY_SIGMA`; `FlowForcing.result`; `SprayPlan`, `region_extent`, `girin_plan`,
  `slurry_modes`; `SkinBody(..., spray=None)`, `.spray_plan`, `.droplet_number`, `.transverse_radius()`; the module
  docstring)
- Test: `tests/test_reentry_model_spray.py` (one test), `tests/test_reentry_model_skin.py` (two),
  `tests/test_reentry_model_skinbody.py` (four)

**Interfaces:**
- Consumes: Task 1a's spray model and its gates; Task 1b's `cap_depth` and `SprayResult.wavelength`; sub-plan 21 Task
  7's `release_top`, `Flow.release`, `FlowBook` and the skin body's release books; 6c's `flow_forcing`;
  `rheology.viscosity`; `film.lubrication`; `spray.melt_layer`.
- Produces: `SprayModel.evaluate(..., rt_layer=None)`, `rt_layer` (n,) [m] the liquid pool the Rayleigh–Taylor mode
  reads; `skin.Cascade(rate, rate_slurry, T_slurry=None)` [kg/(m2 s), kg/(m2 s), K per column];
  `skin.cascade_mass(sk, rate, rate_slurry, A_top, dt, slurry=False, T_slurry=None) -> (mass, from_slurry)` [kg per
  column]; `FlowBook.released_slurry` (n,) [kg]; `skinbody.FlowForcing(dp_ds, tau, g, result=None)`;
  `skinbody.girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False) -> SprayPlan`;
  `SprayPlan` (`rate`, `rate_slurry` [kg/(m2 s)], `result` the liquid's `SprayResult`, `slurry` the exposed slurry's
  per-column dict, `flowing`, `liquid`, `cap_depth` [m], `on_slurry`; the record `mass`, `mass_slurry` [kg], `number`,
  `number_slurry`; `.release(sk) -> Cascade`, `.count(rho)`); `skinbody.region_extent(edges, selected, areas)` [m];
  `skinbody.slurry_modes(...)`; `skinbody.SLURRY_SIGMA` (0.80 N/m); `SkinBody(..., spray=None)` with `.spray_plan`,
  `.droplet_number` and `.transverse_radius()` [m]. No history column, result or command-line setting changes: those
  are Tasks 2–4's.

Asha decided on 2026-10-29 (spec §8 and §16 as amended that day): where liquid lies on slurry the skins read their depths
as the element model does (decision 1, option B), the slurry cascade is built as spec §8 has it (decision 2, option i),
and the two choices Task 1c's outline made within the spec stand — Girin evaluated once a macro step at its start, only
what the admitted run holds re-checked in each flow sub-step; Girin's surface velocity his closure's.

**Girin once a macro step.** At the start of every macro step `SkinBody.advance` evaluates Girin on the corners
(`girin_plan`), from the surface flow the flow's forcing evaluates there (`flow_forcing`, which now hands its whole
result on as `FlowForcing.result`) and the skins as they are then, and gives the plan's release to every part's `Flow`.
`SprayPlan.release` returns the same `Cascade` whatever the skins hold until the next macro step: a column's mode,
wavelength, radius and rates are the macro step's (spec §8), so the records combine per column, per mode, per macro step
without changing a radius. What the accepted passes released is the plan's record, and at the step's end its droplets,
the mass over one droplet's at the mode's radius (`SprayPlan.count`).

**The depths** (option B). b and the regime test read the flowing depth `sk.depth_above(T_rigid)` (spec §6.7), and so
does the viscous film's floor (2026-10-28). The front-surface Rayleigh–Taylor mode reads the liquid depth
`sk.depth_above(T_liquidus)` — its depth criterion, its bounded growth and the depth it sheds within a growth time —
through the shared model's new `rt_layer`, which the element model never passes (ruling 1). The molten region a wave must
fit in is the contiguous region of columns with liquid (`region_extent` on the columns' graph and top areas); and where
Girin's closure does not apply, a column whose liquid run lies on a flowing cell that is not fully liquid is held from
the shear modes (`on_slurry`, judged on the cells, ruling 2). The cap reads the liquid depth with slurry spraying off and
the flowing depth with it on (2026-10-24); the quarter-body cap reads the top surface's transverse radius about the 3D
mesh's mean point (ruling 6). Girin's surface velocity is his closure's, `film.lubrication` on b, as the element model's:
it enters a gate only on the thick branch, where it is τ δ_m / μ_l. The columns on the windward face, θ ≤ 90 degrees,
alone spray.

**The cascade** (option i; spec §8, the 2026-10-07 design's 8.1–8.2). With `slurry_spray`, the exposed slurry has a mode
of its own per column (`slurry_modes`): its first cell beneath the liquid run gives the properties — Li et al.'s
viscosity with the bridge to the liquid (`rheology.viscosity`) at the cell's solid fraction and at the shear rate of the
slurry's own sheared layer (V_s over δ_m or over the layer, film.lubrication's rate), found by fixed point from rest in
at most five passes with the rate held at 367 1/s; the cell's mixture density; 0.80 N/m — and Girin's formulas are the
shared model's with them (its k_r, k_t, critical Weber number and table; the Rayleigh–Taylor mode left to the liquid),
the slurry's layer from its top to T_rigid being its b, regime layer and floor, its wave fitting in the flowing region,
its cap the flowing depth. In every flow sub-step the release works down the top cells in time order, each at its own
phase's rate (`cascade_mass`): the liquid at the liquid's mode's rate, then, once nothing liquid covers it, the slurry at
its own, until the flow sub-step is spent, a cell may not be taken or its rate is zero. Evaluated once a macro step, the
slurry's rate holds for slurry cells at least as warm as the cell that gave it, and a colder cell, more solid and more
viscous, ends the cascade until the next macro step's evaluation (`Cascade.T_slurry`, ruling 3).

**Where the cascade releases** (measured; spec §8 as amended). The viscous film blocks the slurry's thin mode wherever its
viscosity is more than about ten times the liquid's. Swept on the 100 mm sphere's 9 417 columns with the slurry exposed
(`measure/task1c/cascade_where.py`): with the exposed cell at or below 907.0 K (μ 0.033 Pa·s and up, Γτ* 3.2–470) no
column sprays at 93 s or at 120 s under 103 m/s² (Γτ* 3.2–610); from 907.2–907.4 K (μ 9–17 mPa·s, in the upper part
of the bridge from 906.4 to 908 K) the thin mode opens, and from 907.5 K (6.6 mPa·s) nearly every windward column
sprays (4 618–4 636 of 4 705), by the thick branch from 907.7 K. So the cascade releases only slurry within about 0.6–0.8 K of the liquidus, where Scheil's curve leaves a few
per cent of solid. On a heated surface that is the slurry right under the liquid: on sub-plan 21's heated state at 93 s
the first slurry cell beneath every windward column's liquid is at 907.78–908.00 K (median 907.98), so with the flag on
a flow sub-step that strips the 6.6 g of liquid also takes 0.35 g of that slurry (5.1 %; one cell, about 10 µm). On the
skin body's test, the 4 mm sphere heated at 3 MW/m2 from 890 K, the second 0.25 s step sprayed 1.54 g of liquid with the
flag off; with it on, 0.79 g of liquid and 4.98 g of slurry from 615 columns, 471 of them with no liquid at the step's
start but a top in the bridge (907.19–908.00 K): metal heated through the last kelvin before the liquidus is stripped
there. Spec §8 records this as a measured fact; the bridge itself has no data.

**Measured** (2026-10-29, `prototype/work-2026-10-28-skin-spraying/measure/task1c/`):
- *The books.* The skin body's test (`body_run.py`): the 4 mm test sphere heated at 3 MW/m2 from 890 K under the 100 mm
  flight's gas state at 93 s, the flow on (one flow sub-step a heat sub-step), two 0.25 s steps. Its tops melt in the
  first; the second's plan finds liquid on 156 columns and a mode on 144, all thick, and they release 1.54 g as 71 000
  droplets (radius median 126 µm on a liquid depth of 16.5 µm at the start, median); the mass to 1.5e-16, the energy to
  2.9e-11 of the heat absorbed, the momentum to 1.1e-17; Pi 0.008. With slurry spraying on, 5.77 g (above): mass 1.5e-16,
  energy 1.4e-11, momentum 3.8e-17; 89 000 droplets, the slurry's radius median 187 µm.
- *Nothing else moves.* `--surface-model elements` on the CLI test's short melting case (`measure/task1/cmp/`,
  `line1c_short`) is Task 1b's to the last bit in all 115 history columns, the source table and the results; and a skin
  body without a spray model -- flow off with removal, flow on, and flow on with a prescribed release, two 0.25 s steps
  under 3 MW/m2 from 880 K -- is Task 1b's in all 147 arrays (`bitid.py`: the 3D field, the mesh, the skins' state and
  every history column but the wall times).
- *The cost.* On the 100 mm sphere's heated state: the plan 0.09 s of processor time without the cascade and 0.51 s
  with it (4 705 columns' slurry modes, one shared-model evaluation each), against about 2 s a flow sub-step and about
  an hour a melting macro step.

**Rulings** (2026-10-29; each what, why, cost if wrong):
1. *With `rt_layer` the Rayleigh–Taylor mode's rate sheds the liquid pool, ρ_l min(h_pool, λ/2) / τ*, not b.* Why:
   option B's reason, the mode is a liquid pool's, and the element model's own form, its rate on the mobile liquid (the
   film), not on its layer. Cost if wrong: where liquid lies on slurry the mode's demand is smaller by the liquid's share
   of the flowing depth; its release is the liquid's anyway. Confirmed by Asha, 2026-10-29.
2. *A column is on slurry where a flowing cell that is not fully liquid lies beneath its liquid run*, judged on the cells,
   not by comparing the two interpolated depths, which differ by part of a cell at a sharp front. Cost: none.
3. *The cascade's slurry rate holds for slurry cells at least as warm as the first slurry cell it was evaluated on.*
   Why: Asha's confirmed rule (Girin once a macro step, only the admitted run re-checked each flow sub-step) and the
   cascade's per-cell properties both hold: no cell is stripped faster than its own properties allow, to first order.
   Cost if wrong: within one macro step the cascade takes less than the 2026-10-07 design's per-cell evaluation in every
   flow sub-step would -- on the heated state 0.36 g (the first cell, about 10 µm) against up to 1.56 g (the slurry down
   to 907.3 K, about 40 µm; `bridge_mass.py`) -- and catches up over the next macro steps; at the flow-on runs' 12.5 ms
   macro step a heated front advances about one cell a macro step. Confirmed by Asha, 2026-10-29, over the per-cell evaluation.
4. *The slurry's viscosity by fixed point from rest*, the shear rate film.lubrication's for the slurry's layer. Why: the
   2026-10-07 design's rule, started from the most viscous end. Cost: not measured; a start elsewhere could end on
   another root where Li's fit has two (spec §14 assumption 5).
5. *The slurry's droplets at the cell's mixture density, 0.80 N/m*, the design's. Cost: none.
6. *The transverse radius about the 3D mesh's mean point*, as `surface_stats` takes the centre. Why: the skin body has
   no mass centre of its own; on the sphere the two agree. Cost: the quarter-body cap, which binds only above 12.5 mm on
   the 100 mm sphere.
7. *A spray model and a prescribed release are exclusive; a spray model needs the flow and the surface flow's result.*
   Cost: none.
8. *The record is the plan's, per column and mode, not the history's.* Why: Task 2 combines the records into the source
   table and the histograms, Task 4 the history and results (the line drawn in the outline). Cost: none.

- [ ] **Step 1: The Rayleigh–Taylor mode on the liquid pool** — append to `tests/test_reentry_model_spray.py`:

```diff
--- a/tests/test_reentry_model_spray.py
+++ b/tests/test_reentry_model_spray.py
@@ -553,3 +553,31 @@
     res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.array([5.0, 20.0, 0.0]), np.ones(3, bool),
                                          0.5, np.full(3, 1e-5), b * LIQ.rho * 1e-5)
     assert np.isfinite(res.wavelength[:2]).all() and np.isnan(res.wavelength[2])        # no mode on the dry patch
+
+
+def test_the_rayleigh_taylor_mode_can_read_a_liquid_pool_shallower_than_the_film():
+    """Asha, 2026-10-29 (sub-plan 22 Task 1c, decision 1, option B): where liquid lies on slurry the layered skins'
+    front-surface Rayleigh-Taylor mode reads the liquid pool, not the flowing film b -- its depth criterion, its bounded
+    growth and the depth it sheds within a growth time -- as the element model's mode reads its liquid layer (the
+    amendment of 2026-10-06). `rt_layer` is that depth; without it the mode reads the liquid layer as before, and given
+    that same depth it is the same to the last bit. 5 mm of flowing metal at the nose under a lambda* of 2 cm: as one
+    pool the mode takes it; with 0.3 mm of liquid on top it does not (h < lambda*/2 pi, 3.2 mm here), and nothing else
+    can (the thin lambda* at 5 m/s is metres); with 4 mm of liquid it does, at that pool's growth time and rate."""
+    flow = FakeFlow(1, closure=sf.CLOSURE_COUETTE, u=5.0)
+    model = spray.SprayModel(LIQ)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b, areas = np.full(1, 5e-3), np.full(1, 1e-5)
+    w_n = np.full(1, 3.0 * LIQ.sigma * (2.0 * np.pi / 0.02) ** 2 / LIQ.rho)
+    run = lambda **k: model.evaluate(flow, FakeState(), b, delta_m, np.full(1, 0.1), np.ones(1, bool), 0.5, areas,
+                                     b * LIQ.rho * areas, extent=np.full(1, 0.05), deceleration_n=w_n, **k)
+    whole, same = run(), run(rt_layer=b)
+    assert whole.branch[0] == spray.BRANCH_RT and whole.rt_active and whole.dm[0] > 0.0
+    for key in ("branch", "unstable", "r", "mdot", "dm", "dn", "growth", "wavelength"):
+        assert np.array_equal(getattr(whole, key), getattr(same, key), equal_nan=True), key
+    shallow = run(rt_layer=np.full(1, 3e-4))
+    assert shallow.branch[0] == spray.BRANCH_THIN and shallow.mdot[0] == 0.0 and shallow.dm[0] == 0.0
+    assert not shallow.rt_active
+    pool = run(rt_layer=np.full(1, 4e-3))
+    _, lam, tau = spray.rayleigh_taylor_bounded(w_n, np.full(1, 4e-3), np.full(1, 0.05), LIQ)
+    assert pool.branch[0] == spray.BRANCH_RT and pool.growth[0] == tau[0] > whole.growth[0]
+    assert pool.wavelength[0] == lam[0] and pool.mdot[0] == pytest.approx(LIQ.rho * 4e-3 / tau[0], rel=1e-12)
```

- [ ] **Step 2: The cascade's tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -2113,3 +2113,74 @@
     assert max(c for _, c in out["ones"][1:]) < 1.5 * skin.FLOW_COURANT
     assert real_pass("none flowed", [0] * 10) == alone
     assert real_pass("more", [20] * 10) == [20] * 10
+
+
+# -- the cascade: Girin's release from the top cells (sub-plan 22 Task 1c) --------------------------------------------
+def test_the_cascade_takes_each_cell_at_its_own_phase_s_rate_in_time_order_liquid_first_slurry_only_with_the_flag():
+    """Spec 8 and the 2026-10-07 design's 8.1: within a flow sub-step the release works down the top cells in time
+    order, each at its own phase's rate per unit top area -- a liquid cell at the liquid's mode's, a slurry cell at the
+    exposed slurry's own -- until the sub-step is spent or a cell it may not take is reached; slurry only with the flag,
+    never rigid metal, and a cell whose rate is zero ends it, as does a slurry cell colder than the cell whose properties
+    gave the slurry's rate (`T_slurry`). Six columns of 10 um cells: liquid over slurry (902 K) over rigid metal, slurry
+    exposed over liquid, liquid alone over rigid metal, rigid metal, liquid over a warm slurry cell (907 K) over the
+    slurry, and rigid metal. The liquid's rate takes three of a column's cells in the sub-step, the slurry's one and a
+    half: the first column takes its two liquid cells in two thirds of it and half a slurry cell in the last third.
+    `release_top` then takes exactly that mass from the same cells."""
+    col, _ = strip_grid(3, 2)
+    temps = [[LIQUID] * 2 + [SLURRY] * 3 + [RIGID] * 3, [SLURRY] * 2 + [LIQUID] * 2 + [RIGID] * 3,
+             [LIQUID] * 2 + [RIGID] * 3, [RIGID] * 4, [LIQUID] * 2 + [907.0] + [SLURRY] * 2 + [RIGID] * 2, [RIGID] * 4]
+    dt = 0.01
+    cases = ((False, (3.0, 1.5), None, [(2, 0), (0, 0), (2, 0), (0, 0), (2, 0), (0, 0)]),
+             (True, (3.0, 1.5), None, [(2.5, 0.5), (1.5, 1.5), (2, 0), (0, 0), (2.5, 0.5), (0, 0)]),
+             (True, (3.0, 0.0), None, [(2, 0), (0, 0), (2, 0), (0, 0), (2, 0), (0, 0)]),
+             (True, (0.0, 1.5), None, [(0, 0), (1.5, 1.5), (0, 0), (0, 0), (0, 0), (0, 0)]),
+             (True, (3.0, 4.5), 905.0, [(2, 0), (0, 0), (2, 0), (0, 0), (3, 1), (0, 0)]))
+    for slurry, (liquid_cells, slurry_cells), floor, want in cases:
+        sk = stacked(col, temps)
+        m = sk.m[:, 0].copy()                    # every cell of a column holds 10 um of solid metal's mass
+        rate = lambda cells: cells * m / (col.area * dt)
+        T_slurry = None if floor is None else np.full(col.n, floor)
+        mass, from_slurry = skin.cascade_mass(sk, rate(liquid_cells), rate(slurry_cells), col.area, dt, slurry, T_slurry)
+        np.testing.assert_allclose(mass, [w[0] for w in want] * m, rtol=1e-12, atol=0.0)
+        np.testing.assert_allclose(from_slurry, [w[1] for w in want] * m, rtol=1e-12, atol=0.0)
+        M0 = sk.column_mass()
+        book = skin.release_top(sk, mass, slurry=slurry)
+        np.testing.assert_allclose(book.mass, mass, rtol=1e-14, atol=0.0)
+        np.testing.assert_allclose(sk.column_mass(), M0 - mass, rtol=1e-14, atol=1e-30)
+
+
+def test_a_flow_sub_step_with_the_liquid_s_rate_and_the_slurry_s_releases_the_cascade_and_books_its_slurry(
+        coarse_sphere_mesh):
+    """A release that gives two rates per column (sub-plan 22 Task 1c), a `skin.Cascade` of the liquid's mode's and the
+    exposed slurry's [kg/(m2 s)], is taken in each flow sub-step by the cascade (`cascade_mass`) rather than as one rate, after
+    the transport and the landing conduction as before, and the flow's book carries what the slurry's rate took
+    (`released_slurry`); a release of one rate books none. On the test sphere's mixed skins at one and three cycles,
+    without and with the flag: the pieces by hand to the last bit, and the mass and energy fall by what was released."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    rng = np.random.default_rng(10)
+    rates = np.vstack([rng.uniform(0.0, 8.0, col.n), rng.uniform(0.0, 4.0, col.n)])
+    for cycles, slurry in ((1, False), (3, True)):
+        first, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        second, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        flow = sphere_flow(col, cycles, A_top=col.area * 1.01)
+        flow.release, flow.slurry = (lambda sk: skin.Cascade(rates[0], rates[1])), slurry
+        M0, E0 = first.mass(), first.energy()
+        book = skin.flow_substeps(first, flow, 0.01)
+        landed, taken, from_slurry = np.zeros(col.n, dtype=bool), np.zeros(col.n), np.zeros(col.n)
+        for k in range(cycles):
+            if k:
+                second.regrid()
+            by_hand(second, flow, 0.01 / cycles, landed)
+            dm, ds = skin.cascade_mass(second, rates[0], rates[1], flow.A_top, 0.01 / cycles, slurry)
+            taken += skin.release_top(second, dm, slurry=slurry).mass
+            from_slurry += ds
+        assert all(np.array_equal(x, y) for x, y in zip(first.state(), second.state()))
+        np.testing.assert_allclose(book.released, taken, rtol=1e-14, atol=0.0)
+        np.testing.assert_array_equal(book.released_slurry, from_slurry)
+        assert (book.released_slurry.sum() > 0.0) == slurry and book.released.sum() > 0.0
+        assert first.mass() == pytest.approx(M0 - book.released.sum(), rel=3e-16)
+        assert first.energy() == pytest.approx(E0 - book.released_H.sum(), rel=8e-16)
+    flow.release = skin.prescribed_release(rates[0])
+    one = skin.flow_substeps(mixed_sphere_skins(coarse_sphere_mesh)[0], flow, 0.01)
+    assert one.released.sum() > 0.0 and not one.released_slurry.any()
```

- [ ] **Step 3: Girin on the skins' tests** — append to `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -925,3 +925,192 @@
     body, rows, _, failure = sp_.run(dt=0.1, t_end=0.2, T_start=880.0, q0=6e6, **SMALL_PLATE)
     assert failure is None and body.handed_mass > 0.0 and held and not any(held)
     assert body.skin_stats()["held_mass_kg"] == 0.0 and body.skin_stats()["held_depth_max_mm"] == 0.0
+
+
+# -- Girin's spraying from the top cells (sub-plan 22 Task 1c) ----------------------------------------------------------
+def girin_flow(n, closure, u, tau, decel=106.0, branch=None):
+    """The fields Girin's evaluation reads from the surface flow, per column: the spray tests' FakeFlow (edge state,
+    closure, deceleration) with the gas shear and push the lubrication closure reads."""
+    from reentry_model import surface_flow
+    from test_reentry_model_spray import FakeFlow
+    f = FakeFlow(n, branch=surface_flow.BRANCH_SHOCK_LAYER if branch is None else branch, decel=decel)
+    f.closure, f.u_eff = np.asarray(closure), np.asarray(u, dtype=float)
+    f.u_e = f.u_eff
+    f.tau, f.G = np.asarray(tau, dtype=float), np.zeros(n)
+    return f
+
+
+def plan_on(sk, col, flow, slurry=False, theta=None):
+    from reentry_model import spray
+    from test_reentry_model_spray import FakeState
+    return skinbody.girin_plan(sk, flow, FakeState(), spray.SprayModel(sk.material.liquid), np.zeros(sk.n) if theta is None
+                               else theta, col.area, col.edges, 0.05, 0.0125, slurry=slurry)
+
+
+def crest(lam, h):
+    return (3.0 * lam ** 2 * h / (4.0 * np.pi)) ** (1.0 / 3.0)
+
+
+def test_girin_gives_each_prescribed_column_its_mode_from_the_skin_s_depths_once_for_the_macro_step():
+    """Spec 8 and 6.7 with Asha's decision of 2026-10-29 (option B): b and the regime test read the flowing depth, the
+    viscous film's floor too; the Rayleigh-Taylor mode reads the liquid pool and the molten region the wave must fit in
+    is the liquid's; off Girin's closure a liquid on slurry is held from the shear modes; the cap reads the liquid with
+    slurry spraying off. Six columns 20 mm apart at the nose (one molten region about 30 mm across): 1 mm of liquid over
+    rigid metal under Girin's closure (thick); 50 um of liquid on the rigid floor (thin); 5 mm of liquid with no shear
+    mode at 5 m/s (Rayleigh-Taylor, lambda* 2 cm at 106 m/s2); 0.3 mm of liquid over 4 mm of slurry the same (held: read
+    as one flowing film the Rayleigh-Taylor mode would take it); 50 um of liquid over slurry under the Couette closure
+    (held); exposed slurry (nothing: no liquid to take). The plan's rates are the shared spray model's on those depths,
+    and its release gives them whatever the skins hold until the next macro step."""
+    from reentry_model import spray, surface_flow
+    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
+    from test_reentry_model_spray import FakeState
+    col, _ = strip_grid(3, 2, h=20e-3)
+    temps = [[LIQUID] * 86 + [RIGID] * 10, [LIQUID] * 4 + [RIGID] * 10, [LIQUID] * 428 + [RIGID] * 10,
+             [LIQUID] * 26 + [SLURRY] * 340 + [RIGID] * 10, [LIQUID] * 4 + [SLURRY] * 50 + [RIGID] * 10,
+             [SLURRY] * 50 + [RIGID] * 10]
+    sk = stacked(col, temps)
+    G, C = surface_flow.CLOSURE_GIRIN, surface_flow.CLOSURE_COUETTE
+    flow = girin_flow(6, [G, G, C, C, C, G], [2000.0, 2000.0, 5.0, 5.0, 2000.0, 2000.0], [50.0, 50.0, 0.01, 0.01, 50.0, 50.0])
+    plan = plan_on(sk, col, flow)
+    mat, liq = sk.material, sk.material.liquid
+    F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
+    assert plan.result.branch.tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN, spray.BRANCH_RT, -1, -1,
+                                           spray.BRANCH_THICK]
+    assert (plan.rate > 0.0).tolist() == [True, True, True, False, False, False] and not plan.rate_slurry.any()
+    np.testing.assert_array_equal(plan.flowing, F)
+    np.testing.assert_array_equal(plan.liquid, L)
+    np.testing.assert_array_equal(plan.cap_depth, L)
+    assert plan.on_slurry.tolist() == [False, False, False, True, True, True]
+    # the shared model on the skin's depths, by hand
+    delta_m, _ = spray.melt_layer(flow, liq)
+    from reentry_model import film
+    v_s = film.lubrication(flow.tau, flow.G, F, delta_m, liq.mu, b_layer=F)[0]
+    extent = np.where(L > 0.0, 2.0 * np.sqrt(col.area[L > 0.0].sum() / np.pi), 0.0)     # one liquid region
+    run = np.cumprod(sk.valid & (sk.T >= mat.T_liquidus), axis=1).astype(bool)
+    res = spray.SprayModel(liq).evaluate(flow, FakeState(), F, delta_m, v_s, np.ones(6, bool), 0.0125, col.area,
+                                         (sk.m * run).sum(axis=1), 0.05, b_layer=F, extent=extent,
+                                         deceleration_n=np.full(6, 106.0), regime_layer=F, on_slurry=plan.on_slurry,
+                                         rt_layer=L, cap_depth=L)
+    np.testing.assert_allclose(plan.rate, res.mdot, rtol=1e-12, atol=0.0)
+    np.testing.assert_allclose(plan.result.r, res.r, rtol=1e-12, atol=0.0)
+    assert plan.rate[2] == pytest.approx(liq.rho * L[2] / res.growth[2], rel=1e-12)     # the pool shed in a growth time
+    # read as one flowing film the fourth column would be the Rayleigh-Taylor mode's (option A, not taken)
+    one = spray.SprayModel(liq).evaluate(flow, FakeState(), F, delta_m, v_s, np.ones(6, bool), 0.0125, col.area,
+                                         (sk.m * run).sum(axis=1), 0.05, b_layer=F, extent=np.full(6, extent.max()),
+                                         deceleration_n=np.full(6, 106.0), regime_layer=F, cap_depth=L)
+    assert one.branch[3] == spray.BRANCH_RT and one.mdot[3] > 0.0
+    # held for the macro step: the release gives the plan's rates whatever the skins then hold
+    before = plan.release(sk)
+    skin.release_top(sk, sk.column_mass())
+    assert plan.release(sk) is before and np.array_equal(before.rate, plan.rate)
+    assert np.array_equal(before.rate_slurry, plan.rate_slurry) and np.isinf(before.T_slurry).all()
+
+
+def test_the_cap_reads_the_liquid_with_slurry_spraying_off_the_flowing_depth_with_it_on_and_never_the_area():
+    """The decision of 2026-10-24: a droplet holds at most one crest of its own wave on the depth the release may take --
+    the liquid with slurry spraying off, the whole flowing depth with it on -- and the column's area does not enter. Six
+    columns of one liquid cell over 0.5 mm of slurry under Girin's closure, sheared near the thick branch's threshold so
+    that its wavelength is millimetres: the liquid's 15 um binds, the flowing 0.6 mm does not; the corner columns of the
+    strip have smaller areas than the others and the same droplet."""
+    from reentry_model import surface_flow
+    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
+    col, _ = strip_grid(3, 2, h=20e-3)
+    sk = stacked(col, [[LIQUID] + [SLURRY] * 50 + [RIGID] * 10] * 6)
+    flow = girin_flow(6, [surface_flow.CLOSURE_GIRIN] * 6, [2000.0] * 6, [20.0] * 6)
+    off, on = plan_on(sk, col, flow), plan_on(sk, col, flow, slurry=True)
+    mat = sk.material
+    L, F = sk.depth_above(mat.T_liquidus), sk.depth_above(mat.T_rigid)
+    lam = off.result.wavelength
+    assert np.isfinite(lam).all() and (off.rate > 0.0).all() and np.ptp(col.area) > 0.1 * col.area.max()
+    np.testing.assert_allclose(off.result.r, np.minimum(0.17 * lam, crest(lam, L)), rtol=1e-12)
+    np.testing.assert_allclose(on.result.r, np.minimum(0.17 * lam, crest(lam, F)), rtol=1e-12)
+    assert (off.result.r < 0.17 * lam).all() and (on.result.r == 0.17 * lam).all()
+    np.testing.assert_array_equal(on.cap_depth, F)
+    assert np.ptp(off.result.r) <= 1e-12 * off.result.r.max() and np.ptp(on.result.r) <= 1e-12 * on.result.r.max()
+
+
+def test_the_slurry_cascade_sprays_exposed_slurry_by_its_own_properties_only_with_the_flag():
+    """Spec 8's cascade (the 2026-10-07 design's 8.1-8.2, Asha 2026-10-29, option i): with `--slurry-spray on` the slurry
+    the release exposes sprays by Girin's model with its own properties -- Li et al.'s viscosity at its solid fraction and
+    the shear rate of its own sheared layer, found by fixed point (five passes at most, the rate held at 367 1/s), the
+    cell's mixture density and 0.80 N/m -- its layer the slurry from its top down to T_rigid; the Rayleigh-Taylor mode is
+    the liquid's. Five columns of slurry under 15 um of liquid, the slurry's top cell at 902 K (30 % solid), 906.0 K, 907.0
+    K, 907.5 K and 907.9 K (the last three in the bridge to the liquid). Off, nothing; on, the slurry's viscous film blocks
+    its thin mode wherever the slurry is more than a few times the liquid's viscosity, and it sprays only near the
+    liquidus."""
+    from reentry_model import rheology, surface_flow
+    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
+    col, _ = strip_grid(5, 2, h=20e-3)
+    tops = [902.0, 906.0, 907.0, 907.5, 907.9]
+    temps = [[LIQUID] + [T] + [SLURRY] * 48 + [RIGID] * 10 for T in tops] * 2
+    sk = stacked(col, temps)
+    flow = girin_flow(col.n, [surface_flow.CLOSURE_GIRIN] * col.n, [2000.0] * col.n, [50.0] * col.n)
+    off, on = plan_on(sk, col, flow), plan_on(sk, col, flow, slurry=True)
+    assert not off.rate_slurry.any() and np.isnan(off.slurry["mu"]).all()
+    mat = sk.material
+    s = on.slurry
+    f_s = 1.0 - sk.props.liquid_fraction(sk.T[:, 1])                    # the cells' own, as their viscosity's
+    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    np.testing.assert_allclose(s["mu"], rheology.viscosity(sk.T[:, 1], f_s, s["gamma"], mat.liquid.mu, f_b), rtol=1e-12)
+    np.testing.assert_allclose(s["rho"], sk.density(sk.T[:, 1]), rtol=1e-12)
+    assert (s["gamma"] <= rheology.GAMMA_MAX).all() and (s["mu"] > mat.liquid.mu).all()
+    dz = sk.dz()
+    np.testing.assert_allclose(s["depth"], on.flowing - dz[:, 0], rtol=1e-12)          # the slurry, top to T_rigid
+    blocked = s["gamma_tau"] >= 1.0
+    assert (on.rate_slurry[blocked] == 0.0).all() and blocked[[0, 5]].all()             # 30 % solid: blocked
+    assert (on.rate_slurry[~blocked] > 0.0).all() and (~blocked[[4, 9]]).all()          # 907.9 K: sprays
+    assert (np.diff(s["mu"][:5]) < 0.0).all()                                           # runnier toward the liquidus
+    np.testing.assert_array_equal(s["T"], sk.T[:, 1])
+    np.testing.assert_array_equal(on.release(sk).T_slurry, sk.T[:, 1])     # its rate for slurry at least as warm
+    np.testing.assert_array_equal(on.rate, off.rate)            # the liquid's own mode: the flag moves only its cap
+
+
+def test_a_skin_body_sprays_by_girin_from_its_top_cells_once_a_macro_step_and_keeps_its_books(coarse_sphere_mesh,
+                                                                                              monkeypatch):
+    """The skin body with the flow on and a spray model (sub-plan 22 Task 1c): at the start of every macro step Girin is
+    evaluated once on the corners, from the surface flow the flow's forcing evaluates there and the skins' depths then;
+    every part's Flow releases the plan's two rates, the same in every call, through the cascade in every flow sub-step.
+    The 4 mm test sphere heated at 3 MW/m2 from 890 K under the 100 mm flight's gas state at 93 s, two 0.25 s steps, one
+    flow sub-step a heat sub-step: its tops melt in the first, and in the second the liquid is sprayed; the released
+    metal leaves the body, the mass left and released is the starting mass to round-off, the energy books close with the
+    released enthalpy and the momentum balance with the droplets' momentum. The plan's record counts what each column's
+    mode released in the step's accepted passes, and the droplets: their mass over one droplet's at the plan's radius.
+    A spray needs the flow and is not a prescribed release."""
+    pytest.importorskip("cantera")
+    from test_reentry_model_surface_flow import aero_state_at
+    from reentry_model import spray
+    mat = material.Material.from_drama_json("AA7075_scheil")
+    model = spray.SprayModel(mat.liquid)
+    with pytest.raises(ValueError):
+        make(coarse_sphere_mesh, spray=model)
+    with pytest.raises(ValueError):
+        make(coarse_sphere_mesh, flow=True, remove=False, spray=model, release=skin.prescribed_release(1.0))
+    b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=model)
+    plans, flows = [], []
+    girin, advance = skinbody.girin_plan, skin.advance
+
+    def planned(*a, **k):
+        plans.append(girin(*a, **k))
+        return plans[-1]
+
+    def recorded(sk, dt, *a, **k):
+        flows.append(k["flow"])
+        return advance(sk, dt, *a, **k)
+    monkeypatch.setattr(skinbody, "girin_plan", planned)
+    monkeypatch.setattr(skinbody.skin_mod, "advance", recorded)
+    M0 = b.mass()
+    state = aero_state_at(93.0)
+    for k in range(2):
+        flows.clear()
+        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
+        assert len(plans) == k + 1 and b.spray_plan is plans[-1]
+        assert all(f.release == plans[-1].release for f in flows) and len(flows) >= 2
+    p = b.spray_plan
+    assert b.released_mass > 0.0 and p.mass.sum() > 0.0 and not p.mass_slurry.any()
+    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-13)
+    assert abs(b.energy_balance_residual()) < 1e-8 and b.skin_stats()["momentum_residual"] < 1e-12
+    ok = p.mass > 0.0
+    assert (p.rate[ok] > 0.0).all() and np.isfinite(p.result.r[ok]).all()
+    np.testing.assert_allclose(p.number[ok], p.mass[ok] / (4.0 / 3.0 * np.pi * mat.liquid.rho * p.result.r[ok] ** 3),
+                               rtol=1e-12)
+    assert b.droplet_number >= p.number.sum() > 0.0
```

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py
  tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py -q -k "liquid_pool or the_cascade_takes or
  releases_the_cascade or girin_gives or cap_reads or slurry_cascade or sprays_by_girin"`. Expected: 7 failed --
  - `TypeError: SprayModel.evaluate() got an unexpected keyword argument 'rt_layer'`;
  - `AttributeError: module 'reentry_model.skin' has no attribute 'cascade_mass'`, and `... 'Cascade'`;
  - `AttributeError: module 'reentry_model.skinbody' has no attribute 'girin_plan'` (3);
  - `TypeError: SkinBody.__init__() got an unexpected keyword argument 'spray'`.

- [ ] **Step 5: The Rayleigh–Taylor mode's pool** in `reentry_model/spray.py`:

```diff
--- a/reentry_model/spray.py
+++ b/reentry_model/spray.py
@@ -295,7 +295,7 @@
         self.table = table or dispersion.DispersionTable()
 
     def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, b_layer=None,
-                 extent=None, deceleration_n=None, regime_layer=None, on_slurry=None, cap_depth=None):
+                 extent=None, deceleration_n=None, regime_layer=None, on_slurry=None, cap_depth=None, rt_layer=None):
         """Per-patch instability and release for the film thickness b, melt-layer thickness delta_m and film surface
         velocity v_s (from film.lubrication), over the step dt; m_f is the film mass available, radius the body's.
 
@@ -317,7 +317,13 @@
         `cap_depth`, if given, is the depth [m] the release can take from each patch, and the droplet's radius is capped
         by one crest of its own wave on it, (3 lambda^2 h / (4 pi))^1/3, instead of by the film mass m_f (the module
         docstring; 2026-10-24); a patch with no such depth releases nothing. Without it the film-mass cap is kept, bit
-        for bit."""
+        for bit.
+
+        `rt_layer`, if given, is the depth [m] of the liquid pool the front-surface Rayleigh-Taylor mode reads -- its depth
+        criterion, its bounded growth and the depth it sheds within a growth time -- in place of the liquid layer (for the
+        criteria) and the film b (for the rate): the layered skins' liquid depth, where b is the flowing depth and liquid
+        lies on slurry (Asha, 2026-10-29, decision 1 of sub-plan 22 Task 1c, option B: the mode is a liquid pool's, as the
+        amendment of 2026-10-06 has it in the element model). Without it the mode reads them as before, bit for bit."""
         liq = self.liquid
         n = b.size
         layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
@@ -450,12 +456,14 @@
         # shorter growth time takes the patch and the other does not act on the same film. That is also what keeps the
         # shear from being spent twice (the double-counting checklist of the plan's fact 29).
         w_n = flow.deceleration if deceleration_n is None else deceleration_n   # normal component, W cos(phi)
+        pool = layer if rt_layer is None else np.asarray(rt_layer, dtype=float)  # the liquid the mode reads (2026-10-29)
         if self.rt_spray:
-            deep, _, _ = rayleigh_taylor(w_n, layer, liq)               # the threshold: h > lambda*/2 pi
-            rt_on, lam_rt, tau_rt = rayleigh_taylor_bounded(w_n, layer, np.full(n, np.inf) if L is None else L, liq)
+            deep, _, _ = rayleigh_taylor(w_n, pool, liq)                # the threshold: h > lambda*/2 pi
+            rt_on, lam_rt, tau_rt = rayleigh_taylor_bounded(w_n, pool, np.full(n, np.inf) if L is None else L, liq)
             rt_on = np.asarray(rt_on) & np.asarray(deep) & has_film & np.isfinite(tau_rt)
+            shed = b if rt_layer is None else pool                       # what one growth time shakes off
             with np.errstate(divide="ignore", invalid="ignore"):
-                rate_rt = np.where(rt_on, liq.rho * np.minimum(b, lam_rt / 2.0) / np.where(rt_on, tau_rt, 1.0), 0.0)
+                rate_rt = np.where(rt_on, liq.rho * np.minimum(shed, lam_rt / 2.0) / np.where(rt_on, tau_rt, 1.0), 0.0)
             faster = rt_on & (rate_rt > 0.0) & (~np.isfinite(growth) | (tau_rt < growth))
             r = np.where(faster, lam_rt, r)
             mdot = np.where(faster, rate_rt, mdot)
@@ -478,7 +486,7 @@
                 r = np.where(h_cap > 0.0, np.minimum(r, np.minimum(crest, 0.25 * radius)), np.nan)
         with np.errstate(divide="ignore", invalid="ignore"):
             dn = np.where(dm > 0.0, dm / (4.0 / 3.0 * np.pi * liq.rho * np.where(dm > 0.0, r, 1.0) ** 3), 0.0)
-        rt_active, _, _ = rayleigh_taylor(w_n, layer, liq)     # the unbounded criterion, still reported for continuity
+        rt_active, _, _ = rayleigh_taylor(w_n, pool, liq)      # the unbounded criterion, still reported for continuity
         gk = np.isin(branch, (BRANCH_THIN, BRANCH_RAREFIED, BRANCH_REGIME2))    # the Rayleigh-Taylor mode may have taken some
         h_floor, gamma_tau = np.where(gk, h_floor, np.nan), np.where(gk, gamma_tau, np.nan)
         return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s, growth,
```

- [ ] **Step 6: The cascade in the flow sub-steps** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -1748,7 +1748,9 @@
     by `release_top` -- liquid first, the exposed slurry too where `slurry` (spec 8's `--slurry-spray on`, a setting of
     sub-plan 22's). It is called in the trial pass, in every repeat of a part and in parts run again, so it keeps no
     records of its own: what was released is in the pass's FlowBook (`released`, `released_H`, `released_p`,
-    `released_solid`)."""
+    `released_solid`). Since sub-plan 22 Task 1c it may give a `Cascade` instead, the liquid's rate and the exposed
+    slurry's own per column, which the release takes cell by cell in time order (`cascade_mass`; Girin's,
+    `skinbody.SprayPlan.release`), and the book then carries what the slurry's rate took (`released_slurry`)."""
     columns: object
     dp_ds: np.ndarray
     shear: np.ndarray
@@ -1782,17 +1784,19 @@
     released_H: np.ndarray = None      # J (n,): its enthalpy, at its cells' own specific enthalpies
     released_p: np.ndarray = None      # kg m/s (n, 3): its momentum, at its cells' own velocities -- fact 29's sink
     released_solid: np.ndarray = None  # kg (n,): the solid in it, by its cells' liquid fractions
+    released_slurry: np.ndarray = None  # kg (n,): what the cascade took at the exposed slurry's own rate (sub-plan 22
+                                        # Task 1c; zero for a release given as one rate)
 
     @classmethod
     def zero(cls, n):
         return cls(*(np.zeros((n, 3)) for _ in range(4)), released=np.zeros(n), released_H=np.zeros(n),
-                   released_p=np.zeros((n, 3)), released_solid=np.zeros(n))
+                   released_p=np.zeros((n, 3)), released_solid=np.zeros(n), released_slurry=np.zeros(n))
 
     def add(self, other):
         """Add another book's impulses and releases and keep the worst numbers of both."""
         for key in ("push", "shear", "wall", "carried"):
             setattr(self, key, getattr(self, key) + getattr(other, key))
-        for key in ("released", "released_H", "released_p", "released_solid"):
+        for key in ("released", "released_H", "released_p", "released_solid", "released_slurry"):
             mine, theirs = getattr(self, key), getattr(other, key)
             if theirs is not None:
                 setattr(self, key, theirs.copy() if mine is None else mine + theirs)
@@ -1867,10 +1871,18 @@
                         tb.passes, 1)
         if flow.release is not None:
             started = time.perf_counter()
-            rate = np.broadcast_to(np.asarray(flow.release(sk), dtype=float), (n,))
+            given = flow.release(sk)
+            if isinstance(given, Cascade):               # the liquid's rate and the exposed slurry's (sub-plan 22 Task 1c)
+                rate = np.broadcast_to(np.array([given.rate, given.rate_slurry], dtype=float), (2, n))
+            else:
+                rate = np.broadcast_to(np.asarray(given, dtype=float), (n,))
             if not np.all(np.isfinite(rate) & (rate >= 0.0)):
                 raise ValueError("a release rate must be finite and non-negative")
-            rb = release_top(sk, rate * flow.A_top * h, flow.slurry)
+            if isinstance(given, Cascade):
+                dm, step.released_slurry = cascade_mass(sk, rate[0], rate[1], flow.A_top, h, flow.slurry, given.T_slurry)
+            else:
+                dm, step.released_slurry = rate * flow.A_top * h, np.zeros(n)
+            rb = release_top(sk, dm, flow.slurry)
             step.released, step.released_H = rb.mass, rb.enthalpy
             step.released_p, step.released_solid = rb.momentum, rb.solid
             clock("release", started)
@@ -1910,6 +1922,42 @@
     admit = flowing if slurry else flowing & (sk.T >= sk.material.T_liquidus)
     mass, H, solid = sk.take_top_mass(dm, admit)
     return ReleaseBook(mass, H, sk.momentum_out.copy(), solid)
+
+
+@dataclass
+class Cascade:
+    """A release given phase by phase (`Flow.release`; sub-plan 22 Task 1c): per column the liquid's rate and the exposed
+    slurry's own [kg/(m2 s) of top area, (n,)], the slurry's holding for slurry cells at least as warm as `T_slurry` [K,
+    (n,)], the cell whose properties gave it -- a colder slurry cell, more solid and more viscous, ends the column's
+    cascade -- or for every slurry cell where None."""
+    rate: np.ndarray
+    rate_slurry: np.ndarray
+    T_slurry: np.ndarray = None
+
+
+def cascade_mass(sk, rate, rate_slurry, A_top, dt, slurry=False, T_slurry=None):
+    """The mass [kg, (n,)] the release takes from each column's top cells in dt [s] when each cell goes at its own phase's
+    rate, and the part of it taken from slurry cells (spec 8's cascade, the 2026-10-07 design's 8.1; sub-plan 22 Task 1c):
+    down the run of cells `release_top` may take -- the liquid cells from the top, at or above the liquidus, and with
+    `slurry` every flowing cell down to the first rigid one -- in time order, a liquid cell at `rate` and a slurry cell at
+    `rate_slurry` [kg/(m2 s) of top area, (n,)] times the top area A_top [m2], until dt is spent: the liquid first where
+    it covers the slurry, the slurry once nothing liquid covers it, then whatever lies beneath. A cell whose rate is zero
+    ends the column's cascade, as rigid metal does, and so does a slurry cell colder than `T_slurry` [K, (n,)], if given
+    (`Cascade`). `release_top(sk, mass, slurry)` then takes exactly these cells."""
+    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
+    liquid = flowing & (sk.T >= sk.material.T_liquidus)
+    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
+    rate_s = np.asarray(rate_slurry, dtype=float)[:, None]
+    if T_slurry is not None:
+        rate_s = np.where(sk.T >= np.asarray(T_slurry, dtype=float)[:, None], rate_s, 0.0)
+    r = np.where(liquid, np.asarray(rate, dtype=float)[:, None], rate_s)
+    r = r * np.asarray(A_top, dtype=float)[:, None]                           # kg/s from each cell while it is exposed
+    with np.errstate(divide="ignore", invalid="ignore"):
+        lasts = np.where(run & (r > 0.0), sk.m / np.where(r > 0.0, r, 1.0), np.inf)    # s to strip each cell
+    starts = np.concatenate([np.zeros((sk.n, 1)), np.cumsum(lasts, axis=1)[:, :-1]], axis=1)
+    spare = np.where(run & (starts < dt), dt - starts, 0.0)
+    taken = np.where(run, np.minimum(spare * r, sk.m), 0.0)
+    return taken.sum(axis=1), (taken * ~liquid).sum(axis=1)
 
 
 def prescribed_release(rate):
```

- [ ] **Step 7: Girin once a macro step, the slurry's mode and the record** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -69,7 +69,18 @@
 release rate per column (`skin.Flow.release`, what sub-plan 22 supplies from Girin's model) taken in every flow
 sub-step from the top cells, liquid first, the exposed slurry too with `slurry_spray`. The released metal leaves the
 body with its cells' enthalpy and momentum, booked as the removed liquid is, and the history reports the released mass
-and fact 29's droplet-momentum ratio Pi over each macro step (`skin.droplet_momentum_ratio`)."""
+and fact 29's droplet-momentum ratio Pi over each macro step (`skin.droplet_momentum_ratio`).
+
+Girin's spraying from the top cells (sub-plan 22 Task 1c; spec 8; `spray`, a `spray.SprayModel`, with the flow on): at
+the start of every macro step Girin's model is evaluated once on the corners (`girin_plan`), from the surface flow the
+flow's forcing evaluates there and the skins' depths then -- b and the regime test on the flowing depth (spec 6.7), the
+viscous film's floor there (Asha, 2026-10-28), the Rayleigh-Taylor mode on the liquid pool and the wave on the liquid's
+molten region, a liquid on slurry held from the shear modes off Girin's closure (Asha, 2026-10-29, option B), the droplet
+capped by one crest of its wave on the depth the release may take (2026-10-24) -- and, with `slurry_spray`, the exposed
+slurry's own mode by the cascade's properties (spec 8; Asha, 2026-10-29, option i). The plan's two rates per column are
+the release of every flow sub-step of the macro step (`SprayPlan.release`, through `skin.cascade_mass`); what the
+accepted passes released is the plan's record, per column and mode (`spray_plan`)."""
+import copy
 import dataclasses
 import math
 import time
@@ -77,8 +88,12 @@
 import numpy as np
 
 from . import columns as columns_mod
+from . import film as film_mod
+from . import material as material_mod
 from . import movemesh
+from . import rheology
 from . import skin as skin_mod
+from . import spray as spray_mod
 from .body import ThermalBody
 from .thermal import SIGMA_SB
 
@@ -125,20 +140,216 @@
                                   # (measured 2026-10-14). The element model's runs keep their solver's own tolerance.
 
 
+SLURRY_SIGMA = 0.80               # N/m: the exposed slurry's surface tension in the cascade (the 2026-10-07 design 8.2, the
+                                  # large-fragment design's labelled assumption for slurry)
+
+
 @dataclasses.dataclass
 class FlowForcing:
     """The gas side's push on the skins' metal at a macro step, per corner: `dp_ds` [Pa/m, (n,)] the wall pressure's
     gradient along the meridian and `tau` [Pa, (n,)] the gas shear (the surface flow's), and `g` [m/s2, (3,)] the
     acceleration the metal feels in the body's frame -- minus the body's own: the deceleration along the flight
-    direction."""
+    direction. `result`, the surface flow's whole evaluation (a `surface_flow.SurfaceFlowResult` on the corners), is what
+    Girin's spraying reads (sub-plan 22 Task 1c); None for a push given without one."""
     dp_ds: np.ndarray
     tau: np.ndarray
     g: np.ndarray
+    result: object = None
 
     def toward(self, other, frac):
         """This forcing moved the fraction `frac` of the way to `other`, linearly."""
         return FlowForcing(*(a + (b - a) * frac for a, b in ((self.dp_ds, other.dp_ds), (self.tau, other.tau),
                                                              (self.g, other.g))))
+
+
+@dataclasses.dataclass
+class SprayPlan:
+    """Girin's spraying from the top cells over one macro step, evaluated at its start (`girin_plan`; spec 8, sub-plan 22
+    Task 1c), per column: `rate` [kg/(m2 s) of top area] the liquid's mode's and `rate_slurry` the exposed slurry's own
+    (the cascade; zero without slurry spraying); `result` the liquid's `spray.SprayResult` (branch, radius, wavelength);
+    `slurry` the exposed slurry's per column -- its `branch`, radius `r` [m], `wavelength` [m], viscosity `mu` [Pa s] at
+    the shear rate `gamma` [1/s] of its own sheared layer, density `rho` [kg/m3], conjugate depth `delta_m` [m], layer
+    `depth` [m] from its top to T_rigid, the viscous film's `gamma_tau` and the cell's temperature `T` [K] (nan, -1 where
+    there is none); the depths read
+    at the start, `flowing`, `liquid` and the cap's `cap_depth` [m], and `on_slurry`, where a flowing cell that is not
+    fully liquid lies beneath the liquid run. The record the macro step's accepted passes fill: `mass` and
+    `mass_slurry` [kg] released by each mode, and, at the step's end, `number` and `number_slurry`, the droplets: the mass
+    over one droplet's at the mode's radius."""
+    rate: np.ndarray
+    rate_slurry: np.ndarray
+    result: object
+    slurry: dict
+    flowing: np.ndarray
+    liquid: np.ndarray
+    cap_depth: np.ndarray
+    on_slurry: np.ndarray
+    mass: np.ndarray = None
+    mass_slurry: np.ndarray = None
+    number: np.ndarray = None
+    number_slurry: np.ndarray = None
+
+    def __post_init__(self):
+        T = self.slurry["T"]
+        self._cascade = skin_mod.Cascade(self.rate, self.rate_slurry, np.where(np.isfinite(T), T, np.inf))
+        n = len(self.rate)
+        self.mass, self.mass_slurry = np.zeros(n), np.zeros(n)
+        self.number, self.number_slurry = np.zeros(n), np.zeros(n)
+
+    def release(self, sk):
+        """`skin.Flow.release`: the liquid's rate and the exposed slurry's (`skin.Cascade`), whatever the skins hold within
+        the macro step -- Girin's modes and radii are the macro step's, from its start (spec 8) -- the slurry's for slurry
+        cells at least as warm as the cell that gave it."""
+        return self._cascade
+
+    def count(self, rho_liquid):
+        """The record's droplets: each mode's mass over one droplet's at its radius."""
+        for mass, r, rho, key in ((self.mass, self.result.r, rho_liquid, "number"),
+                                  (self.mass_slurry, self.slurry["r"], self.slurry["rho"], "number_slurry")):
+            ok = (mass > 0.0) & np.isfinite(r) & (r > 0.0)
+            with np.errstate(invalid="ignore", divide="ignore"):
+                setattr(self, key, np.where(ok, mass / (4.0 / 3.0 * np.pi * rho * np.where(ok, r, 1.0) ** 3), 0.0))
+
+
+def region_extent(edges, selected, areas):
+    """Lateral extent [m] of each column's own contiguous selected region on the columns' graph (`edges`, (E, 2)), as the
+    region's equivalent diameter 2 sqrt(A / pi) with A the selected columns' `areas` summed; zero where not selected: the
+    element model's `MeltingBody.region_extent` on the corners, the domain a surface wave must fit inside (2026-09-24)."""
+    from scipy.sparse import coo_matrix
+    from scipy.sparse.csgraph import connected_components
+    sel = np.asarray(selected, dtype=bool)
+    n = len(sel)
+    out = np.zeros(n)
+    if not sel.any():
+        return out
+    i, j = edges[:, 0], edges[:, 1]
+    keep = sel[i] & sel[j]
+    i, j = i[keep], j[keep]
+    g = coo_matrix((np.ones(2 * i.size), (np.concatenate([i, j]), np.concatenate([j, i]))), shape=(n, n))
+    count, label = connected_components(g, directed=False)
+    area = np.bincount(label[sel], np.asarray(areas, dtype=float)[sel], count)
+    out[sel] = 2.0 * np.sqrt(area[label[sel]] / np.pi)
+    return out
+
+
+def _rows(flow, idx):
+    """The surface flow's evaluation on the corners idx alone: every per-corner array taken at idx, the rest kept."""
+    n = len(flow.closure)
+    out = copy.copy(flow)
+    for key, value in vars(flow).items():
+        if isinstance(value, np.ndarray) and value.ndim >= 1 and value.shape[0] == n:
+            setattr(out, key, value[idx])
+    return out
+
+
+def girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False):
+    """Girin's spraying from each column's top cells for the macro step about to run (spec 8; sub-plan 22 Task 1c), from
+    the surface flow's evaluation on the corners (`flow`), the aero state (`state`, for the rarefied branch), the shared
+    spray model (`model`) and the skins as they are at its start: a SprayPlan.
+
+    The liquid's mode: b and the regime test read the flowing depth `sk.depth_above(T_rigid)` (spec 6.7), and so does
+    the viscous film's floor (Asha, 2026-10-28); delta_m is Girin's Eq. (2) and his surface velocity his closure's
+    (`film.lubrication` on b, as the element model's); the front-surface Rayleigh-Taylor mode reads the liquid pool, the
+    liquid depth `sk.depth_above(T_liquidus)` (`rt_layer`), and every wave must fit in the contiguous region of liquid
+    columns, and where Girin's closure does not apply a liquid lying on slurry takes no shear mode (`on_slurry`) -- the
+    element model's reading (amendment of 2026-10-06; Asha, 2026-10-29, option B); the droplet is capped by one crest of
+    its wave on the depth the release may take, the liquid's, or with `slurry` the flowing depth (2026-10-24); a quarter
+    of `radius` caps it too. Columns on the windward face, theta <= 90 degrees, only; `A_top` [m2] the top areas, `edges`
+    the columns' graph, `dt` [s] the macro step (the nominal release in the result). With `slurry`, the exposed slurry's
+    own mode too (`slurry_modes`)."""
+    mat, liq = sk.material, model.liquid
+    n = sk.n
+    F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
+    flowing = sk.valid & np.isfinite(skin_mod.cell_viscosity(sk, 0.0))
+    liquid = flowing & (sk.T >= mat.T_liquidus)
+    n_flow = np.cumprod(flowing, axis=1).sum(axis=1)
+    n_liq = np.cumprod(liquid, axis=1).sum(axis=1)
+    on_slurry = n_flow > n_liq
+    theta = np.asarray(theta, dtype=float)
+    windward = theta <= 0.5 * np.pi
+    w_n = flow.deceleration * np.maximum(np.cos(theta), 0.0)
+    delta_m, _ = spray_mod.melt_layer(flow, liq)
+    v_s = film_mod.lubrication(flow.tau, flow.G, F, delta_m, liq.mu, b_layer=F)[0]
+    extent = region_extent(edges, L > 0.0, A_top)
+    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
+    cap = F if slurry else L
+    res = model.evaluate(flow, state, F, delta_m, v_s, windward, dt, A_top, (sk.m * run).sum(axis=1), radius,
+                         b_layer=F, extent=extent, deceleration_n=w_n, regime_layer=F, on_slurry=on_slurry, rt_layer=L,
+                         cap_depth=cap)
+    if slurry:
+        rate_s, modes = slurry_modes(sk, flow, state, model, windward, w_n, A_top, edges, radius, dt, F, n_liq, n_flow)
+    else:
+        rate_s, modes = np.zeros(n), _no_slurry(n)
+    return SprayPlan(res.mdot, rate_s, res, modes, F, L, cap, on_slurry)
+
+
+def _no_slurry(n):
+    out = {key: np.full(n, np.nan) for key in ("r", "wavelength", "mu", "gamma", "rho", "delta_m", "depth", "gamma_tau",
+                                               "T")}
+    out["branch"] = np.full(n, -1)
+    return out
+
+
+def slurry_modes(sk, flow, state, model, windward, w_n, A_top, edges, radius, dt, F, n_liq, n_flow):
+    """The exposed slurry's own mode per column, the cascade of spec 8 (the 2026-10-07 design's 8.1-8.2; Asha, 2026-10-29,
+    option i): (rate [kg/(m2 s)], the per-column dict of SprayPlan.slurry). Where a flowing cell that is not fully liquid
+    lies beneath the liquid run on the windward face, its cell -- the slurry the release exposes first -- gives the
+    properties: Li et al.'s viscosity (`rheology.viscosity`, with the bridge to the liquid) at the cell's solid fraction
+    and the shear rate of the slurry's own sheared layer, V_s over delta_m or over the layer, found by fixed-point
+    iteration from rest (at most skin.SHEAR_PASSES passes, until no viscosity moves by skin.SHEAR_TOL, the rate held at
+    rheology.GAMMA_MAX); the cell's mixture density; SLURRY_SIGMA. The layer runs from the cell's top to T_rigid, and the
+    regime test, the viscous film's floor and b are it. Girin's formulas are the shared model's with these properties
+    (its k_r, k_t, critical Weber number and dispersion table), the Rayleigh-Taylor mode left to the liquid; the wave
+    must fit in the flowing region; the cap is the flowing depth. Evaluated once a macro step, as the liquid's mode: its
+    rate holds for slurry cells at least as warm as the one whose properties gave it, and a colder cell, more solid and
+    more viscous, waits for the next macro step's evaluation (`skin.Cascade.T_slurry`)."""
+    mat = sk.material
+    n = sk.n
+    out = _no_slurry(n)
+    rate = np.zeros(n)
+    rows = np.flatnonzero((n_flow > n_liq) & windward)
+    if not rows.size:
+        return rate, out
+    k = n_liq[rows]
+    T = sk.T[rows, k]
+    dz = sk.dz()
+    depth = F[rows] - (np.cumsum(dz, axis=1) - dz)[rows, k]
+    f_s = 1.0 - sk.props.liquid_fraction(T)
+    rho = sk.density(T)
+    f_bridge = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    sub = _rows(flow, rows)
+    viscosity = lambda gamma: rheology.viscosity(T, f_s, gamma, mat.liquid.mu, f_bridge)
+
+    def layer(mu):
+        d, _ = spray_mod.melt_layer(sub, material_mod.LiquidProperties(rho, mu, SLURRY_SIGMA))
+        v, _, shear, _ = film_mod.lubrication(sub.tau, sub.G, depth, d, mu, b_layer=depth)
+        return d, v, np.minimum(shear, rheology.GAMMA_MAX)
+
+    gamma = np.zeros(rows.size)
+    mu = viscosity(gamma)
+    for _ in range(skin_mod.SHEAR_PASSES):
+        gamma = layer(mu)[2]
+        mu_new = viscosity(gamma)
+        done = np.all(np.abs(mu_new - mu) <= skin_mod.SHEAR_TOL * mu)
+        mu = mu_new
+        if done:
+            break
+    d_m, v_s, _ = layer(mu)
+    extent = region_extent(edges, F >= spray_mod.B_MIN, A_top)
+    cells = np.arange(sk.m.shape[1])[None, :]
+    m_s = (sk.m * ((cells >= n_liq[:, None]) & (cells < n_flow[:, None]))).sum(axis=1)     # the slurry run, nominally
+    for j, i in enumerate(rows):
+        props = material_mod.LiquidProperties(float(rho[j]), float(mu[j]), SLURRY_SIGMA)
+        m = spray_mod.SprayModel(props, model.k_r, model.k_t, model.we_critical, model.table, rt_spray=False)
+        one, cell = slice(j, j + 1), slice(i, i + 1)
+        r = m.evaluate(_rows(sub, one), state, depth[one], d_m[one], v_s[one], windward[cell], dt, A_top[cell], m_s[cell],
+                       radius, b_layer=depth[one], extent=extent[cell], deceleration_n=w_n[cell], regime_layer=depth[one],
+                       cap_depth=F[cell])
+        rate[i] = r.mdot[0]
+        out["branch"][i], out["r"][i], out["wavelength"][i] = r.branch[0], r.r[0], r.wavelength[0]
+        out["gamma_tau"][i] = r.gamma_tau[0]
+    out["mu"][rows], out["gamma"][rows], out["rho"][rows], out["T"][rows] = mu, gamma, rho, T
+    out["delta_m"][rows], out["depth"][rows] = d_m, depth
+    return rate, out
 
 
 def jump_faults(T_bottom, T_node, T_deep, crossed, tol=None, handed=None):
@@ -175,12 +386,14 @@
     `skin.flow_cycles` chooses each heat sub-step, Task 6d), driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its
     defaults if None) from each macro step's aero state; `release`, with the flow on, a release from the top cells in
     every flow sub-step (`skin.Flow.release`: a callable of the skins giving a rate [kg/(m2 s)] per unit top area per
-    column; sub-plan 21 Task 7), the exposed slurry released too where `slurry_spray`."""
+    column; sub-plan 21 Task 7), the exposed slurry released too where `slurry_spray`; or `spray`, a `spray.SprayModel`
+    with the flow on, Girin's spraying from the top cells evaluated at the start of every macro step (`girin_plan`; sub-plan
+    22 Task 1c), with `slurry_spray` the exposed slurry's own mode too (the cascade), and its record in `spray_plan`."""
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
-                 surface_flow=None, release=None, slurry_spray=False):
+                 surface_flow=None, release=None, slurry_spray=False, spray=None):
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
         self.flow_cycles = flow_cycles
@@ -190,7 +403,14 @@
         if release is not None and not self.flow_on:
             raise ValueError("a release from the top cells needs the flow on: without it the heat-only skins remove the "
                              "liquid as it forms (stage 2's stand-in for flow and spray)")
+        if spray is not None and not self.flow_on:
+            raise ValueError("Girin's spraying from the top cells needs the flow on: without it the heat-only skins "
+                             "remove the liquid as it forms")
+        if spray is not None and release is not None:
+            raise ValueError("a skin body sprays by Girin's model or by a release given, not both")
         self.release, self.slurry_spray = release, bool(slurry_spray)
+        self.spray, self.spray_plan = spray, None            # Girin's model, and its plan of the last macro step (Task 1c)
+        self._step_release = release                         # what the parts' Flow releases this macro step
         self.repeat_rule = repeat_rule
         self.part_max = part_max                             # s, or None: no cap on the parts while melting
         self.material, self.solver, self.mass_kg = material, solver, mass_kg
@@ -228,6 +448,7 @@
         self.impulse, self.momentum_lost = np.zeros(3), np.zeros(3)   # the momentum books since time zero [N s, kg m/s]
         self.impulse_size = 0.0                              # N s: the push's, shear's and walls' impulses, in size
         self.released_mass = self.released_enthalpy = 0.0    # kg, J: released from the top cells since time zero (Task 7)
+        self.droplet_number = 0.0                            # droplets Girin's spraying released since time zero (Task 1c)
         self.droplet_momentum = np.zeros((self.columns.n, 3))   # kg m/s: what the droplets took from each column, and
         self.droplet_shear = np.zeros((self.columns.n, 3))      # N s: the gas shear's impulse on its top, this macro step
         self.flow_courant = self.flow_growth = 0.0           # this macro step's worst flow sub-step
@@ -282,6 +503,13 @@
             self.flow_courant = self.flow_growth = 0.0
             self.flow_passes = self.flow_cycles_max = 0
             self.droplet_momentum, self.droplet_shear = np.zeros((self.columns.n, 3)), np.zeros((self.columns.n, 3))
+            if self.spray is not None:                       # Girin once, at the macro step's start (sub-plan 22 Task 1c)
+                if end.result is None:
+                    raise ValueError("Girin's spraying reads the surface flow's evaluation: give the aero state, or a "
+                                     "forcing with its result")
+                self.spray_plan = girin_plan(self.skins, end.result, state, self.spray, self.theta, self.surface.areas,
+                                             self.columns.edges, self.transverse_radius(), dt, self.slurry_spray)
+                self._step_release = self.spray_plan.release
         sk, mat = self.skins, self.material
         top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
         if self.part_max is not None and mat.melts:
@@ -301,6 +529,9 @@
             self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt, loads, 0 if n == 1 else 1)
         self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * self._ratio)))
         self._q_prev = q_end.copy()
+        if self.spray is not None and self.flow_on:          # the macro step's droplets, from its accepted passes' mass
+            self.spray_plan.count(self.spray.liquid.rho)
+            self.droplet_number += float(self.spray_plan.number.sum() + self.spray_plan.number_slurry.sum())
 
     def flow_forcing(self, state):
         """The gas side's push on the metal at this macro step (FlowForcing): the surface flow evaluated on the corners
@@ -316,7 +547,7 @@
         res = self.surface_flow.evaluate(state, self.theta, self.nose_radius(), rho_l, self.surface_temperature())
         self.last_flow = res
         return FlowForcing(np.asarray(res.dp_ds, dtype=float), np.asarray(res.tau, dtype=float),
-                           res.deceleration * self.v_hat)
+                           res.deceleration * self.v_hat, res)
 
     def _advance_part(self, t, dt, dt_macro, loads, depth):
         """One part of a macro step, ending at t: the passes, then -- if its draws would move a base more than MOVE_MAX
@@ -332,7 +563,7 @@
             f = self._forcing[0].toward(self._forcing[1], frac)
             sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
             flow = skin_mod.Flow(self.columns, f.dp_ds, skin_mod.along_meridian(self.columns.normal, self.v_hat, f.tau),
-                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles, release=self.release,
+                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles, release=self._step_release,
                                  slurry=self.slurry_spray)
         eps = np.full(sk.n, self.emissivity)
         m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
@@ -426,6 +657,9 @@
             self.momentum_lost += fb.released_p.sum(axis=0)  # the droplets' momentum, at their cells' own velocities
             self.droplet_momentum += fb.released_p
             self.droplet_shear += fb.shear
+            if self.spray is not None:                       # the plan's record: each mode's mass (sub-plan 22 Task 1c)
+                self.spray_plan.mass += np.maximum(fb.released - fb.released_slurry, 0.0)
+                self.spray_plan.mass_slurry += fb.released_slurry
         self.iterations.append(res.iterations)
         self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
         self._move_base(book, m_min, g0)
@@ -571,6 +805,12 @@
 
     def nose_radius(self):
         return float(self.mesh.params.get("radius_m", 0.05))
+
+    def transverse_radius(self):
+        """The top surface's largest distance from the axis through the 3D mesh's mean point along the flight direction
+        [m]: the body radius whose quarter caps a droplet (spec 8; the element model's transverse radius)."""
+        X = self.surface.points - self.mesh.points.mean(axis=0)
+        return float(np.sqrt(np.maximum(np.einsum("ij,ij->i", X, X) - (X @ self.v_hat) ** 2, 0.0)).max())
 
     def demised(self):
         return self.mass() < self.demise_fraction * self.mass0
```

- [ ] **Step 8: Run them to see them pass.** Expected: `tests/test_reentry_model_spray.py` 24 passed,
  `tests/test_reentry_model_skin.py` 92, `tests/test_reentry_model_skinbody.py` 32 (2 deselected, the reference tier's;
  the body's spraying test about 100 s); the unit tier 533 passed, 1 skipped.

- [ ] **Step 9: The measurements** (in `measure/task1c/`, each on mains power): `"$PY" body_run.py` and
  `"$PY" body_run.py slurry` (the books and the record; about 2 minutes each); `"$PY" cascade_where.py heated`,
  `"$PY" cascade_where.py sweep 93` and `"$PY" cascade_where.py sweep 120 103` (where the cascade releases);
  `"$PY" bridge_mass.py` (ruling 3); `"$PY" bitid.py run <Task 1b's code> bit_1b.npz`, `"$PY" bitid.py run <code>
  bit_1c.npz` and `"$PY" bitid.py compare bit_1b.npz bit_1c.npz`; and the element model's short case,
  `./run_case.sh <code> runs line1c_short` and `"$PY" compare_runs.py runs linecap_short line1c_short` in
  `measure/task1/cmp/`.

- [ ] **Step 10: Commit.**

```bash
git add reentry_model/spray.py reentry_model/skin.py reentry_model/skinbody.py tests/test_reentry_model_spray.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py
git commit -m "Skins: Girin's spraying from each column's top cells (spec 8) -- evaluated once at the start of every macro step from the surface flow and the skins' depths, b, the regime test and the viscous film's floor on the flowing depth, the Rayleigh-Taylor mode and the molten region on the liquid and a liquid on slurry held off Girin's closure (Asha 2026-10-29, option B), the droplet capped per crest on the depth the release may take; the release in every flow sub-step through the cascade, liquid first and the exposed slurry by its own properties only with slurry spraying on (option i), which sprays only within about a kelvin of the liquidus; a per-column record per mode per macro step; the element model and a skin body without a spray model unchanged to the last bit (sub-plan 22 Task 1c)"
```

---

### Task 2: The droplet records from the skins

**Files:**
- Modify: `reentry_model/spray.py` (`SOURCE_COLUMNS` + `solid_fraction`, `h_kg`; `source_block`; `source_rows(...,
  solid_fraction=None, h_kg=None)`; `source_table`; the module docstring), `reentry_model/skin.py`
  (`Skins.take_top_mass(..., phase=None)`, `._by_phase`; `ReleaseBook.by_phase`; `release_top`; `FlowBook.phase_mass`,
  `.phase_H`, `.phase_p`, `.phase_solid`, `.hist_n`, `.hist_m`; `Cascade.histogram`; `flow_substeps`),
  `reentry_model/skinbody.py` (`SprayPlan.theta`, `.closure`, `.points`, `.state`, `.props`, the record's `enthalpy`,
  `solid` and `momentum` per mode, `.droplet_mass`, `.histogram`, `.book`, `.source_block`; `girin_plan(..., points=None)`;
  the slurry's `we_s` and `h_floor`; `SkinBody(..., particles=True)`, `.source_rows`, `.hist_n`, `.hist_m`; the module
  docstring), `reentry_model/body.py` (the element model's rows' two new columns), `reentry_model/coupled.py`
  (`write_particles`)
- Test: `tests/test_reentry_model_spray.py` (one test), `tests/test_reentry_model_skin.py` (two),
  `tests/test_reentry_model_skinbody.py` (two), `tests/test_reentry_model_melting.py` (one changed),
  `tests/test_reentry_model_coupled.py` (one changed)

**Interfaces:**
- Consumes: Task 1c's `girin_plan`, `SprayPlan` and its record, `skin.Cascade`, `cascade_mass`, `release_top` and
  `FlowBook`; the element model's `spray.source_rows` and `MeltingBody.film_enthalpy`; `coupled.write_particles`.
- Produces: `spray.SOURCE_COLUMNS` (26: + `solid_fraction`, `h_kg` [J/kg]); `spray.source_block(t, h, V, theta,
  centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None) -> (k, 26)`, `liquid`'s `rho`, `mu`,
  `sigma` one value or one per patch; `spray.source_rows(..., solid_fraction=None, h_kg=None)` (its rows as lists);
  `spray.source_table(rows) -> (N, 26)` from a list of rows or of blocks; `Skins.take_top_mass(dm, admit, phase=None)`
  leaving `phase_out`; `ReleaseBook.by_phase` (a ReleaseBook of (2, ...) arrays, [0] the fully liquid cells, [1] the
  others); `FlowBook.phase_mass`, `phase_H`, `phase_solid` [kg, J, kg, (2, n)], `phase_p` [kg m/s, (2, n, 3)],
  `hist_n`, `hist_m` [(N_BINS,), or None]; `skin.Cascade(..., histogram=None)`; `SprayPlan.book(fb)`,
  `.histogram(mass) -> (dn, dM)`, `.droplet_mass(rho_liquid) -> (2, n)` [kg], `.source_block(t) -> (k, 26)`,
  `.enthalpy`, `.enthalpy_slurry` [J], `.solid`, `.solid_slurry` [kg], `.momentum`, `.momentum_slurry` [kg m/s, (n,
  3)]; `girin_plan(..., points=None)`; `SkinBody(..., particles=True)` with `.source_rows` (one (k, 26) array per macro
  step), `.hist_n` and `.hist_m` [kg] per size bin; `write_particles` returning only the files it wrote. No history
  column, frame field, result or command-line setting changes: those are Tasks 3 and 4's.

Spec 8 and 10.2 (the 2026-10-07 design's 7.5): the droplet records are combined per column, per mode, per macro step,
summing number and mass over the sub-steps; a column's radius is fixed within the macro step by its plan, so combining
changes no radius; the size histograms accumulate each sub-step; two source columns are added, the droplets' solid
fraction and the enthalpy per kilogram they carry.

**The records.** Task 1c's cascade takes each top cell at its own phase's rate, so what a flow sub-step released from a
column's fully liquid cells is the liquid's mode's and what it released from the others the exposed slurry's. The release
now books that split from its own cells (`Skins.take_top_mass(..., phase=liquid)`, `ReleaseBook.by_phase`): mass,
enthalpy at each cell's own specific enthalpy, solid by each cell's liquid fraction and momentum at each cell's own
velocity, a column that took cells of one kind giving that kind its totals to the last bit (ruling 1). The flow's book
carries it over the flow sub-steps (`FlowBook.phase_*`), and the skin body adds an accepted pass's book to the plan's
record (`SprayPlan.book`), so the record holds, per column and mode, the macro step's mass, enthalpy, solid and momentum,
and at the step's end its droplets (`count`, as Task 1c's). The liquid's mode and the slurry's are separate records, each
at its own radius from the plan, which holds for the macro step: combining changes no radius.

**The histograms.** The plan's `Cascade` carries the plan's `histogram`, which `flow_substeps` calls on every flow
sub-step's release by mode (`spray.histogram` at the plan's radii and densities), the book summing the bins; the skin body
adds an accepted pass's bins to its `hist_n` and `hist_m`, so the sub-steps of trial, repeated and rejected passes, whose
releases are undone, are not counted (ruling 2). With the radius the macro step's, the bins are the records' own
histogram to round-off.

**The rows.** At the step's end the record becomes rows in the source table's format, written by the shared
`spray.source_block` that the element model's `source_rows` now is (`SprayPlan.source_block`), one per column and mode
that released, the liquid's mode's first. On a skin record (rulings 3-8):
- `time_s`, `altitude_km`, `velocity_kms`: the macro step's end and its aero state, as the element model's;
- `theta_deg`, `x_m`, `y_m`, `z_m`: the column's top corner on the top surface and its angle to the flight direction where
  the plan was evaluated, at the macro step's start (`girin_plan(..., points=)`), where its mode and radius come from;
- `closure`, `branch`, `delta_m_m`, `we_s`, `h_floor_m`, `gamma_tau`, `r_m`: the plan's, for the liquid's mode from the
  shared model's result and for the slurry's from its own evaluation (whose We_s and floor the slurry's record now keeps);
- `b_m`: the depth the mode read as Girin's b: for the liquid's mode the flowing depth (spec 6.7), for the slurry's its
  own layer, from its first cell's top to T_rigid;
- `dn`, `dm_kg`: the record's droplets and mass;
- `v_s_ms` and `tx`, `ty`, `tz`: the release velocity relative to the body, as the source table defines it (the 2026-09-20
  design 9: body velocity + V_s t), the mean velocity of the metal released, its momentum over its mass: the momentum the
  books lose with the droplets, which leave at the speed of the cells they are made of (spec 6.4 and 8; sub-plan 21 Task
  7); a zero direction where it is at rest. Girin's closure velocity remains the gates' (Task 1c), in `we_s`;
- `we_d`, `oh`, `breakup`: at the record's radius, by its phase's density, viscosity and surface tension, the liquid's or
  the slurry's (its mixture density, its viscosity at its sheared layer's rate, 0.80 N/m: the 2026-10-07 design 8.4);
- `solid_fraction`, `h_kg`: the record's solid over its mass and its enthalpy over its mass, mass-weighted over the cells
  released.

The element model's rows carry 0 and its film's own liquid enthalpy `h_p`, which its removed enthalpy books for the
sprayed film (`MeltingBody.film_enthalpy`: the film is liquid by construction and carries its latent heat): its other 24
columns and its history are unchanged to the last bit (ruling 8). The skin body keeps one block of rows per macro step
(`source_rows`, ruling 9) where `particles` (on by default), and `coupled.write_particles` reads either form through
`spray.source_table`: a skin flight writes the same `particles.npz` and `size_distribution.csv`; the per-step summary,
whose columns are the element model's history columns, waits for Task 4 (ruling 10).

**Measured** (2026-10-29, `prototype/work-2026-10-28-skin-spraying/measure/task2/`; facts 174–176):
- *The records against the books* (`body_records.py`, `log_body_records.txt`). Task 1c's spraying run of the 4 mm test
  sphere (3 MW/m2 from 890 K, the gas state at 93 s, two 0.25 s steps, one flow sub-step a heat sub-step). With slurry
  spraying on, the second step's block holds 759 rows: 144 of the liquid's mode (the thick branch) and 615 of the
  slurry's, the 144 liquid columns giving both. Every row's radius is its plan's; the rows' mass is the books' released
  mass, 5.76888 g, to the last bit printed; their enthalpy, the sum of dm h_kg, is the books' released enthalpy to 1.5e-16
  and their momentum, dm v_s t summed, the droplets' momentum to 4.7e-17; the histogram accumulated each flow sub-step has
  the records' 14 bins, the counts within 1.6e-12 droplets of the records' (3.4e-11 of a bin holding 0.065, 1.8e-17 of the
  89 173 in all: numpy's weighted histogram takes a bin as a difference of running sums over every droplet, so its error is
  the whole step's, not the bin's) and the mass within 2.4e-14. The liquid's records: solid fraction 0 to 2.3e-11 (the
  property table's interpolation across the liquidus), h_kg 1.043-1.044 MJ/kg (h(T_liquidus) 1.04302 MJ/kg), radius
  120/130/2 734 um (minimum, median, maximum), We_d 3.7/4.0/84, Oh 4.0e-4/1.8e-3/1.9e-3. The slurry's: solid fraction
  8.5e-6/0.0061/0.041, h_kg 1.026/1.041/1.043 MJ/kg, radius 125/187/3 532 um, We_d 3.8/5.7/108 (36 rows above the
  breakup threshold of 12), Oh 5.2e-4/2.5e-3/9.2e-3. With slurry spraying off: 144 rows, the books to 0, 1.4e-16 and
  2.6e-16, the histograms in 7 bins to 3.4e-13.
- *The release velocity* (the same runs). The liquid's records leave at 0.45/1.28/1.62 m/s, 1.04 m/s weighted by mass,
  within 0.4/3.8/20 degrees of the meridian; Girin's closure velocity on the same columns is 2.7/7.2/9.7 m/s, 5.6 m/s by
  mass, so the resolved speed is 0.14-0.24 of it (0.14-0.23 with the flag off). The slurry's leave at 0.88 m/s by mass. On
  sub-plan 21's heated state at 93 s (Task 1's `speed_steady.py`) the top cells ran at 0.48 of the closure's speed with
  0.15 mm of liquid and fell to 0.06 of it as the liquid ran off.
- *Nothing else moves.* The skin body's state at the end of both runs -- its 3D field, mesh, skins, every history column
  but the wall times, its released mass, enthalpy, momentum and droplets -- is Task 1c's to the last bit; with the flag on
  the plan's record differs from Task 1c's in round-off (1.8e-15 at most, in 97 columns: the split by the released cells
  against the cascade's slurry share subtracted), no record appearing or vanishing. A skin body without a spray model is
  Task 1c's in all 147 arrays of Task 1c's three runs (`bitid.py`, `run_bitid.sh`: none of the 147 differs). `--surface-model elements` on the CLI test's
  short melting case (`runs/line2_short`) is Task 1c's in all 115 history columns (its history file byte for byte), all
  24 old fields of its 25 266 source rows, its results, its summary and its size distribution; its rows carry a solid
  fraction of 0 and h_kg 0.94/1.04/1.35 MJ/kg, whose sum of dm h_kg is its removed enthalpy, 540 843 J, to the last bit.
- *The cost* (`cost_records.py`). On the 100 mm sphere's heated state at 93 s with slurry spraying on, a 1 ms flow
  sub-step releases 6.9 g, the liquid from 4 618 columns and the slurry from 4 613, all of them both: the release costs
  0.079 s of processor time without the split and 0.113 s with it, the histogram 1.1 ms and the macro step's 9 231 rows
  1.4 ms, against about 2 s a flow sub-step; with slurry spraying off no column takes cells of both kinds and the split
  adds two masks. Each the mean of 20 repeats, measured once.

**Rulings** (2026-10-29; each what, why, cost if wrong):
1. *The records split by the cells the release took*, not by the cascade's arithmetic. Why: the cascade gives each cell
   its phase's rate, so the cell's phase is the mode that took it; booked from the release's own cells a record is exact
   and a column of one kind keeps its totals to the last bit, where Task 1c's `released - released_slurry` subtracts two
   sums made differently and could leave a round-off residue as a liquid record with no liquid. Cost: none (the two
   agree to 1.8e-15 on the test).
2. *The histograms accumulate in the flow sub-steps through the Cascade, an accepted pass's added to the body.* Why: spec
   8's "each sub-step", without counting the sub-steps of passes that are undone. `spray.histogram` is unchanged, so the
   element model's files are. Cost: none; the agreement with the records is numpy's round-off (above).
3. *Position and angle: the column's top corner and its angle at the plan's evaluation*, the macro step's start. Why:
   the droplets leave the top, and the mode and radius are from there; the element model's centroid is on its surface.
   Cost if the base corner was meant: the skin's thickness, about a millimetre on the 100 mm sphere.
4. *b_m the depth the mode read*: the flowing depth for the liquid's mode, the slurry's layer for the slurry's. Cost:
   none.
5. *v_s_ms and t: the mean velocity of the metal released.* Why: the source table's column is the release velocity, and
   spec 6.4 and 8 have the droplets leave the top cells at the speed of the metal they are made of, the books' sink, so
   the table's momentum is the momentum the body lost (to 4.7e-17); Girin's closure velocity is a gate's, Task 1c's
   ruling, and is in We_s. Cost if wrong: the table's release speed is 0.14-0.24 of the closure's on the test sphere
   (0.06-0.48 on sub-plan 21's heated state), a factor of 4-7 slower than the element model's convention; the
   alternative is one line (`res.v_s` and the meridian). Confirmed by Asha, 2026-10-29 (option A).
6. *A slurry droplet's We_d and Oh by its slurry's properties*, its mixture density, its viscosity at its sheared
   layer's rate (the mode's) and 0.80 N/m (the 2026-10-07 design 8.4). Cost: Oh with another viscosity scales with it.
7. *solid_fraction and h_kg by the books' cell values*: the solid by the property table's liquid fraction, as
   `released_solid`, so a liquid record may carry up to 2.3e-11 of solid. Cost: none.
8. *The element model's rows: 0 and its film's liquid enthalpy.* Why: its film is liquid by construction, booked at
   `h_p` when it sprays, so the column is what its droplets carry and its books' sprayed enthalpy; nan would hide it.
   Cost: none (its other columns unchanged to the last bit).
9. *A skin body keeps one array of rows per macro step*, read by `spray.source_table`. Why: a list of floats costs about
   0.9 kB a row against 0.2 kB, and a skin flight makes a row per spraying column per macro step (about 4 600 at 93 s,
   so 370 000 a second of flight at the 12.5 ms continuum step). Cost: `len(source_rows)` counts macro steps on a skin
   body (Task 4's results count rows).
10. *No summary for a skin body until Task 4*: `particles_summary.csv`'s columns are the element model's per-step
    history columns, which a skin history lacks, and its `released_mass_kg` is cumulative (sub-plan 21 Task 7) where the
    summary's is per step; `write_particles` writes the summary only for a history that has all its columns and returns
    the paths it wrote. Cost: none until Task 4, which gives the skin history its per-step droplet columns.
11. *A spraying skin body needs the aero state*, now refused without it: the rows' altitude and speed. Cost: none (every
    caller gives it).

**What Tasks 3 and 4 add for a skin flight's files.** Task 3: the CLI builds the skin body with `spray=SprayModel(...)`,
`slurry_spray` from `--slurry-spray` and `particles` from `--particles`, and calls `write_particles` for a skin run (it
does for the element model's only, `melting` excluding the skins). Task 4: the skin history's per-step droplet columns
for the summary (released mass and droplets per step, the radius median and maximum, the spraying area) and the results'
row count, sprayed mass and droplets.

- [ ] **Step 1: The two columns and the per-patch properties** -- append to `tests/test_reentry_model_spray.py`:

```diff
--- a/tests/test_reentry_model_spray.py
+++ b/tests/test_reentry_model_spray.py
@@ -581,3 +581,37 @@
     _, lam, tau = spray.rayleigh_taylor_bounded(w_n, np.full(1, 4e-3), np.full(1, 0.05), LIQ)
     assert pool.branch[0] == spray.BRANCH_RT and pool.growth[0] == tau[0] > whole.growth[0]
     assert pool.wavelength[0] == lam[0] and pool.mdot[0] == pytest.approx(LIQ.rho * 4e-3 / tau[0], rel=1e-12)
+
+
+def test_the_source_table_carries_the_droplets_solid_fraction_and_enthalpy_and_takes_their_properties_per_patch():
+    """Spec 8 and 10.2 (sub-plan 22 Task 2): two columns end the source table, `solid_fraction` and `h_kg`, the droplets'
+    solid share and the enthalpy per kilogram they carry, 0 and nan where the caller gives none; the droplets' density,
+    viscosity and surface tension may be one per patch, a slurry droplet's Weber and Ohnesorge numbers being its own
+    phase's (the 2026-10-07 design 8.4); `source_block` is `source_rows` as an array, and `spray.source_table` reads a
+    body's rows or blocks of rows alike. The element model's form is unchanged in every other column."""
+    assert spray.SOURCE_COLUMNS[-2:] == ["solid_fraction", "h_kg"] and len(spray.SOURCE_COLUMNS) == 26
+    flow = FakeFlow(3)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.full(3, 1e-3)
+    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.full(3, 5.0), np.array([True, False, True]), 0.5,
+                                         np.full(3, 1e-5), b * LIQ.rho * 1e-5)
+    where = (10.0, 71e3, 7200.0, np.radians([30.0, 45.0, 60.0]), np.zeros((3, 3)), np.tile([0.0, 1.0, 0.0], (3, 1)), flow,
+             res, b)
+    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
+    plain = np.array(spray.source_rows(*where, LIQ, FakeState()))
+    assert plain.shape == (2, 26) and (plain[:, C["solid_fraction"]] == 0.0).all() and np.isnan(plain[:, C["h_kg"]]).all()
+    slurry = material.LiquidProperties(np.array([2400.0, 2500.0, 2550.0]), np.array([1.3e-3, 0.02, 0.05]),
+                                       np.array([0.86, 0.80, 0.80]))
+    block = spray.source_block(*where, slurry, FakeState(), solid_fraction=np.array([0.0, 0.1, 0.03]),
+                               h_kg=np.array([1.1e6, 1.0e6, 0.9e6]))
+    assert block.shape == (2, 26) and block[:, C["solid_fraction"]].tolist() == [0.0, 0.03]
+    assert block[:, C["h_kg"]].tolist() == [1.1e6, 0.9e6]
+    np.testing.assert_array_equal(block[:, :C["we_d"]], plain[:, :C["we_d"]])
+    np.testing.assert_array_equal(block[0, :C["solid_fraction"]], plain[0, :C["solid_fraction"]])   # the liquid's patch
+    r = block[1, C["r_m"]]
+    assert block[1, C["we_d"]] == pytest.approx(8e-5 * 7200.0 ** 2 * 2.0 * r / 0.80, rel=1e-14)
+    assert block[1, C["oh"]] == pytest.approx(0.05 / np.sqrt(2550.0 * 0.80 * 2.0 * r), rel=1e-14)
+    np.testing.assert_array_equal(spray.source_block(*where, LIQ, FakeState()), plain)        # rows as an array
+    table = spray.source_table([block, block[:1]])
+    assert table.shape == (3, 26) and np.array_equal(table[:2], block, equal_nan=True)
+    assert np.array_equal(spray.source_table(plain.tolist()), plain, equal_nan=True) and spray.source_table([]).shape == (0, 26)
```

- [ ] **Step 2: The release by the cells it came from, and the histograms** -- append to
  `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -2184,3 +2184,82 @@
     flow.release = skin.prescribed_release(rates[0])
     one = skin.flow_substeps(mixed_sphere_skins(coarse_sphere_mesh)[0], flow, 0.01)
     assert one.released.sum() > 0.0 and not one.released_slurry.any()
+
+
+# -- the droplet records: the release by the cells it came from (sub-plan 22 Task 2) ----------------------------------
+def test_a_release_is_split_by_the_cells_it_came_from_and_a_column_of_one_kind_keeps_its_totals_to_the_last_bit():
+    """Spec 8's separate records for the liquid's mode and the slurry's (sub-plan 22 Task 2): `release_top` books what it
+    took by the cells it came from -- [0] the fully liquid ones, [1] the others -- mass, enthalpy, solid and momentum at
+    each cell's own. Five columns of 10 um cells, each kind of cell moving its own way: liquid over slurry, exposed
+    slurry, liquid alone, rigid metal, slurry over liquid. With slurry spraying on, asked for two and a half cells: the
+    first column gives its two liquid cells to [0] and half a slurry cell to [1], the fifth two slurry cells and half a
+    liquid one; a column that took one kind gives that kind its totals to the last bit and the other nothing; the liquid
+    carries no solid. With it off, everything is liquid's."""
+    col, _ = strip_grid(5, 2)
+    temps = [[LIQUID] * 2 + [SLURRY] * 3 + [RIGID] * 3, [SLURRY] * 3 + [RIGID] * 3, [LIQUID] * 4 + [RIGID] * 3,
+             [RIGID] * 4, [SLURRY] * 2 + [LIQUID] * 2 + [RIGID] * 3] * 2
+    for slurry in (True, False):
+        sk = stacked(col, temps)
+        liquid = sk.valid & (sk.T >= sk.material.T_liquidus)
+        sk.U = np.where(liquid[..., None], [2.0, 0.0, 0.0], [0.0, 0.5, 0.0]) * sk.valid[..., None]
+        m, h, U = sk.m[:, 0].copy(), sk.H / np.where(sk.m > 0.0, sk.m, 1.0), sk.U.copy()
+        f_s = 1.0 - sk.props.liquid_fraction(sk.T)
+        book = skin.release_top(sk, 2.5 * m, slurry=slurry)
+        ph = book.by_phase
+        for total, part in ((book.mass, ph.mass), (book.enthalpy, ph.enthalpy), (book.solid, ph.solid),
+                            (book.momentum, ph.momentum)):
+            assert part.shape == (2,) + total.shape
+            np.testing.assert_allclose(part.sum(axis=0), total, rtol=1e-15, atol=1e-30)
+        assert not ph.solid[0].any()
+        one = [1, 2, 3, 6, 7, 8] if slurry else list(range(10))                # one kind taken, or none
+        kind = (ph.mass[1] > 0.0).astype(int)
+        for i in one:
+            for total, part in ((book.mass, ph.mass), (book.enthalpy, ph.enthalpy), (book.momentum, ph.momentum)):
+                assert np.array_equal(part[kind[i], i], total[i]) and not part[1 - kind[i], i].any()
+        if not slurry:
+            assert not ph.mass[1].any() and np.array_equal(ph.mass[0], book.mass)
+            continue
+        for i, (n_liq, n_sl, first) in ((0, (2.0, 0.5, 0)), (4, (0.5, 2.0, 2))):     # the mixed columns, by hand
+            cells = [(0, 1.0), (1, 1.0), (2, 0.5)]
+            liq = [(k, w) for k, w in cells if liquid[i, k]]
+            sl = [(k, w) for k, w in cells if not liquid[i, k]]
+            for j, run in ((0, liq), (1, sl)):
+                np.testing.assert_allclose(ph.mass[j, i], sum(w for _, w in run) * m[i], rtol=1e-14)
+                np.testing.assert_allclose(ph.enthalpy[j, i], sum(w * m[i] * h[i, k] for k, w in run), rtol=1e-14)
+                np.testing.assert_allclose(ph.solid[j, i], sum(w * m[i] * f_s[i, k] for k, w in run), rtol=1e-14,
+                                           atol=1e-30)
+                np.testing.assert_allclose(ph.momentum[j, i], sum(w * m[i] * U[i, k] for k, w in run), rtol=1e-14)
+            assert ph.mass[0, i] == pytest.approx(n_liq * m[i], rel=1e-14) and ph.mass[1, i] == pytest.approx(n_sl * m[i])
+
+
+def test_a_cascade_with_a_histogram_bins_each_flow_sub_step_s_release_and_the_book_keeps_it_by_phase(coarse_sphere_mesh):
+    """Sub-plan 22 Task 2: the flow's book carries every flow sub-step's release by the cells it came from
+    (`phase_mass`, `phase_H`, `phase_p`, `phase_solid`), which add up to the release, and a `skin.Cascade` with a
+    `histogram` has each flow sub-step's release binned by size, the book summing the bins over the sub-steps. On the
+    test sphere's mixed skins at three cycles with slurry spraying on: the histogram is called once a flow sub-step, on
+    that sub-step's release by phase, and the slurry's part is what the cascade took at the slurry's rate (to round-off:
+    the release's own cells against the cascade's arithmetic)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    rng = np.random.default_rng(10)
+    rates = np.vstack([rng.uniform(0.0, 8.0, col.n), rng.uniform(0.0, 4.0, col.n)])
+    seen = []
+
+    def histogram(mass):
+        seen.append(mass.copy())
+        return np.array([mass[0].sum(), mass[1].sum()]), np.array([1.0, 2.0])
+    sk, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    flow = sphere_flow(col, 3, A_top=col.area * 1.01)
+    flow.release, flow.slurry = (lambda s: skin.Cascade(rates[0], rates[1], histogram=histogram)), True
+    book = skin.flow_substeps(sk, flow, 0.01)
+    assert len(seen) == 3 and book.hist_m.tolist() == [3.0, 6.0]
+    np.testing.assert_allclose(book.hist_n, [sum(s[0].sum() for s in seen), sum(s[1].sum() for s in seen)], rtol=1e-15)
+    np.testing.assert_allclose(book.phase_mass, sum(seen), rtol=1e-15)
+    for total, part in ((book.released, book.phase_mass), (book.released_H, book.phase_H),
+                        (book.released_solid, book.phase_solid), (book.released_p, book.phase_p)):
+        np.testing.assert_allclose(part.sum(axis=0), total, rtol=1e-14, atol=1e-25)
+    assert book.phase_mass[1].sum() > 0.0 and book.phase_mass[0].sum() > 0.0 and not book.phase_solid[0].any()
+    np.testing.assert_allclose(book.phase_mass[1], book.released_slurry, rtol=1e-12, atol=1e-22)
+    flow.release = skin.prescribed_release(rates[0])                          # a release of one rate: no histogram
+    one = skin.flow_substeps(mixed_sphere_skins(coarse_sphere_mesh)[0], flow, 0.01)
+    assert one.hist_n is None and one.phase_mass.sum() > 0.0
```

- [ ] **Step 3: The plan's records and the body's** -- append to `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1114,3 +1114,168 @@
     np.testing.assert_allclose(p.number[ok], p.mass[ok] / (4.0 / 3.0 * np.pi * mat.liquid.rho * p.result.r[ok] ** 3),
                                rtol=1e-12)
     assert b.droplet_number >= p.number.sum() > 0.0
+
+
+# -- the droplet records (sub-plan 22 Task 2) ---------------------------------------------------------------------------
+def test_a_plan_s_records_are_one_per_column_and_mode_with_the_released_cells_solid_fraction_enthalpy_and_velocity():
+    """Spec 8 and 10.2: the droplet records are combined per column, per mode, per macro step, summing mass and number
+    over the sub-steps; a column's radius is its plan's, so combining changes none; the liquid's mode and the slurry's
+    are separate records; `solid_fraction` and `h_kg` are mass-weighted over the cells each record released, and the
+    droplets leave at the mean velocity of the metal they are made of. Three columns of 10 um cells under Girin's
+    closure, slurry spraying on: liquid alone (the thin branch), exposed slurry three cells of it in the bridge at 907.9 K
+    (the slurry's thick mode), and one liquid cell over that slurry -- the liquid moving at 2 m/s along x, the slurry at
+    0.5 m/s along y. Three flow sub-steps of 12 us: the third column gives its liquid cell in the first two and slurry
+    after it. Four rows, the liquid's mode's first; the size histograms accumulated each sub-step are the records'."""
+    from reentry_model import spray, surface_flow
+    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
+    from test_reentry_model_spray import FakeState
+    col, pts = strip_grid(3, 2, h=20e-3)
+    WARM = 907.9
+    sk = stacked(col, [[LIQUID] * 5 + [RIGID] * 10, [WARM] * 3 + [SLURRY] * 47 + [RIGID] * 10,
+                       [LIQUID] + [WARM] * 3 + [SLURRY] * 46 + [RIGID] * 10, [RIGID] * 10, [RIGID] * 10, [RIGID] * 10])
+    mat = sk.material
+    liquid = sk.valid & (sk.T >= mat.T_liquidus)
+    sk.U = np.where(liquid[..., None], [2.0, 0.0, 0.0], [0.0, 0.5, 0.0]) * sk.valid[..., None]
+    flow = girin_flow(6, [surface_flow.CLOSURE_GIRIN] * 6, [2000.0] * 6, [50.0] * 6)
+
+    class State(FakeState):
+        h = 71e3
+    theta = np.radians([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
+    plan = skinbody.girin_plan(sk, flow, State(), spray.SprayModel(mat.liquid), theta, col.area, col.edges, 0.05, 0.0125,
+                               slurry=True, points=pts)
+    assert plan.result.branch[:3].tolist() == [spray.BRANCH_THIN, spray.BRANCH_THICK, spray.BRANCH_THICK]
+    assert (plan.rate > 0.0).tolist() == [True, False, True] + [False] * 3
+    assert (plan.rate_slurry > 0.0).tolist() == [False, True, True] + [False] * 3
+    hist = np.zeros((2, spray.N_BINS))
+    taken, m0 = np.zeros((2, 6)), sk.m[:, 0].copy()
+    for _ in range(3):
+        cascade = plan.release(sk)
+        dm, _ = skin.cascade_mass(sk, cascade.rate, cascade.rate_slurry, col.area, 12e-6, True, cascade.T_slurry)
+        ph = skin.release_top(sk, dm, slurry=True).by_phase
+        fb = skin.FlowBook.zero(6)
+        fb.phase_mass, fb.phase_H, fb.phase_p, fb.phase_solid = ph.mass, ph.enthalpy, ph.momentum, ph.solid
+        plan.book(fb)
+        hist += cascade.histogram(ph.mass)
+        taken += ph.mass
+    assert taken[0, 2] == pytest.approx(m0[2], rel=1e-12) and 0.0 < taken[1, 2] < m0[2]        # all its liquid
+    plan.count(mat.liquid.rho)
+    rows = plan.source_block(0.0125)
+    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
+    i, mode = np.array([0, 2, 1, 2]), np.array([0, 0, 1, 1])
+    assert rows.shape == (4, 26)
+    np.testing.assert_array_equal(rows[:, C["x_m"]:C["z_m"] + 1], pts[i])
+    np.testing.assert_array_equal(rows[:, C["theta_deg"]], np.degrees(theta[i]))
+    assert (rows[:, C["time_s"]] == 0.0125).all() and (rows[:, C["altitude_km"]] == 71.0).all()
+    s = plan.slurry
+    r = np.where(mode == 0, plan.result.r[i], s["r"][i])
+    np.testing.assert_array_equal(rows[:, C["r_m"]], r)                                  # combining changes no radius
+    np.testing.assert_array_equal(rows[:, C["dm_kg"]], taken[mode, i])
+    rho = np.where(mode == 0, mat.liquid.rho, s["rho"][i])
+    np.testing.assert_allclose(rows[:, C["dn"]], taken[mode, i] / (4.0 / 3.0 * np.pi * rho * r ** 3), rtol=1e-14)
+    for key, liq, sl in (("branch", plan.result.branch, s["branch"]), ("delta_m_m", plan.result.delta_m, s["delta_m"]),
+                         ("we_s", plan.result.we_s, s["we_s"]), ("b_m", plan.flowing, s["depth"]),
+                         ("h_floor_m", plan.result.h_floor, s["h_floor"]), ("gamma_tau", plan.result.gamma_tau,
+                                                                             s["gamma_tau"])):
+        np.testing.assert_array_equal(rows[:, C[key]], np.where(mode == 0, liq[i], sl[i]))
+    # the two new columns: liquid at 950 K carries no solid and its enthalpy; the slurry its cells' at 907.9 K
+    f_s = 1.0 - float(sk.props.liquid_fraction(np.array([WARM]))[0])
+    assert 0.0 < f_s < 0.1
+    np.testing.assert_allclose(rows[:, C["solid_fraction"]], np.where(mode == 0, 0.0, f_s), rtol=1e-12, atol=0.0)
+    np.testing.assert_allclose(rows[:, C["h_kg"]], mat.enthalpy(np.where(mode == 0, LIQUID, WARM)), rtol=1e-12)
+    # leaving at the metal's own velocity, and each phase's Weber and Ohnesorge numbers at its radius
+    np.testing.assert_allclose(rows[:, C["v_s_ms"]], np.where(mode == 0, 2.0, 0.5), rtol=1e-14)
+    np.testing.assert_allclose(rows[:, C["tx"]:C["tz"] + 1], np.where(mode[:, None] == 0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
+                               atol=1e-15)
+    sigma = np.where(mode == 0, mat.liquid.sigma, skinbody.SLURRY_SIGMA)
+    mu = np.where(mode == 0, mat.liquid.mu, s["mu"][i])
+    np.testing.assert_allclose(rows[:, C["we_d"]], State.freestream.rho * State.V ** 2 * 2.0 * r / sigma, rtol=1e-14)
+    np.testing.assert_allclose(rows[:, C["oh"]], mu / np.sqrt(rho * sigma * 2.0 * r), rtol=1e-14)
+    assert s["mu"][2] > 1.2 * mat.liquid.mu and s["rho"][2] > mat.liquid.rho
+    # the histograms accumulated each sub-step are the records' (in the same bins, to round-off)
+    n_rec, m_rec = spray.histogram(rows[:, C["r_m"]], rows[:, C["dn"]], rows[:, C["dm_kg"]])
+    np.testing.assert_array_equal(hist[0] > 0.0, n_rec > 0.0)
+    np.testing.assert_allclose(hist[0], n_rec, rtol=1e-13, atol=0.0)
+    np.testing.assert_allclose(hist[1], m_rec, rtol=1e-13, atol=0.0)
+
+
+def test_a_spraying_skin_body_s_records_close_on_its_books_and_the_coupled_writer_writes_them(coarse_sphere_mesh,
+                                                                                            tmp_path, monkeypatch):
+    """The skin body's droplet records (sub-plan 22 Task 2), on Task 1c's spraying test with slurry spraying on: the 4 mm
+    test sphere heated at 3 MW/m2 from 890 K under the 100 mm flight's gas state at 93 s, two 0.25 s steps, one flow
+    sub-step a heat sub-step. The second step sprays the liquid and the slurry under it, and the body keeps one block of
+    rows for it: one row per column and mode, the liquid's mode's first, each at its plan's radius; the rows' mass is
+    what the accepted passes released from each column (the slurry's part what the cascade took at the slurry's rate)
+    and the books' released mass, their enthalpy the books' released enthalpy, and their momentum, mass times speed
+    along their direction, the momentum the droplets took from each column, all to round-off; the size histograms
+    accumulated each flow sub-step are the rows'. The coupled writer writes the rows and their size distribution and
+    reads them back to the last bit; the summary, whose columns are the element model's history, waits for Task 4."""
+    pytest.importorskip("cantera")
+    from test_reentry_model_surface_flow import aero_state_at
+    from reentry_model import coupled, spray
+    from reentry_model import trajectory as tj
+    mat = material.Material.from_drama_json("AA7075_scheil")
+    b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid),
+             slurry_spray=True)
+    plans, books = [], []
+    girin, book = skinbody.girin_plan, skinbody.SprayPlan.book
+
+    def planned(*a, **k):
+        plans.append(girin(*a, **k))
+        return plans[-1]
+
+    def booked(plan, fb):
+        books.append((plan, fb.released.copy(), fb.released_slurry.copy()))
+        return book(plan, fb)
+    monkeypatch.setattr(skinbody, "girin_plan", planned)
+    monkeypatch.setattr(skinbody.SprayPlan, "book", booked)
+    state = aero_state_at(93.0)
+    for k in range(2):
+        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
+    assert len(b.source_rows) == 1 and b.particles                      # the first step melts the tops
+    rows = spray.source_table(b.source_rows)
+    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
+    p = plans[-1]
+    index = {tuple(x): i for i, x in enumerate(p.points)}
+    i = np.array([index[tuple(x)] for x in rows[:, C["x_m"]:C["z_m"] + 1]])
+    n_liq, n_sl = int((p.mass > 0.0).sum()), int((p.mass_slurry > 0.0).sum())
+    mode = (np.arange(len(rows)) >= n_liq).astype(int)
+    assert len(rows) == n_liq + n_sl and n_liq > 100 and n_sl > 100 and (rows[:, C["time_s"]] == 0.5).all()
+    assert len(set(i[mode == 0])) == n_liq and len(set(i[mode == 1])) == n_sl
+    assert ((p.mass > 0.0) & (p.mass_slurry > 0.0)).sum() > 10                    # columns with both records
+    np.testing.assert_array_equal(rows[:, C["r_m"]], np.where(mode == 0, p.result.r[i], p.slurry["r"][i]))
+    # the books: the mass each column released and what the cascade took at the slurry's rate, the enthalpy, the momentum
+    mine = [x for x in books if x[0] is p]
+    released, from_slurry = sum(x[1] for x in mine), sum(x[2] for x in mine)
+    by_column = np.bincount(i, rows[:, C["dm_kg"]], b.columns.n)
+    np.testing.assert_allclose(by_column, released, rtol=1e-12, atol=0.0)
+    np.testing.assert_allclose(np.bincount(i[mode == 1], rows[mode == 1, C["dm_kg"]], b.columns.n), from_slurry,
+                               rtol=1e-12, atol=1e-21)
+    assert rows[:, C["dm_kg"]].sum() == pytest.approx(b.released_mass, rel=1e-12)
+    assert (rows[:, C["dm_kg"]] * rows[:, C["h_kg"]]).sum() == pytest.approx(b.released_enthalpy, rel=1e-12)
+    velocity = rows[:, C["v_s_ms"], None] * rows[:, C["tx"]:C["tz"] + 1]
+    momentum = np.stack([np.bincount(i, rows[:, C["dm_kg"]] * velocity[:, k], b.columns.n) for k in range(3)], axis=1)
+    np.testing.assert_allclose(momentum, b.droplet_momentum, rtol=0.0, atol=1e-12 * np.abs(b.droplet_momentum).max())
+    # the two new columns: the liquid carries no solid (to the property table's interpolation at the liquidus), the
+    # slurry less than half its mass, each at its cells' own enthalpy
+    s_f, h = rows[:, C["solid_fraction"]], rows[:, C["h_kg"]]
+    h_rigid, h_liq = mat.enthalpy(np.array([mat.T_rigid, mat.T_liquidus]))
+    assert (s_f[mode == 0] < 1e-9).all() and (s_f[mode == 1] > 0.0).all() and (s_f[mode == 1] < 0.5).all()
+    assert (h[mode == 0] > h_liq * (1.0 - 1e-6)).all() and (h_rigid < h[mode == 1]).all() and (h[mode == 1] < h_liq).all()
+    # the size histograms, accumulated each flow sub-step, are the rows': in the same bins, to the round-off of numpy's
+    # weighted histogram, which takes each bin as a difference of running sums over all the droplets (measured: 1.6e-12
+    # droplets in a bin of 0.065, 2e-17 of the 71 000 in all)
+    n_rec, m_rec = spray.histogram(rows[:, C["r_m"]], rows[:, C["dn"]], rows[:, C["dm_kg"]])
+    np.testing.assert_array_equal(b.hist_n > 0.0, n_rec > 0.0)
+    np.testing.assert_allclose(b.hist_n, n_rec, rtol=0.0, atol=1e-15 * n_rec.sum())
+    np.testing.assert_allclose(b.hist_m, m_rec, rtol=0.0, atol=1e-15 * m_rec.sum())
+    # the coupled writer: the same files as the element model's, read back to the last bit; no summary until Task 4
+    history = tj.History({"time_s": np.array([0.0, 0.25, 0.5])}, np.zeros((3, 6)), None)
+    files = coupled.write_particles(str(tmp_path), b, history)
+    assert set(files) == {"particles", "size_distribution"}
+    npz = np.load(files["particles"])
+    assert npz.files == spray.SOURCE_COLUMNS
+    assert all(np.array_equal(npz[name], rows[:, k], equal_nan=True) for k, name in enumerate(spray.SOURCE_COLUMNS))
+    with open(files["size_distribution"]) as fh:
+        total = [line.split(",") for line in fh.read().splitlines()[1:]][-spray.N_BINS:]
+    np.testing.assert_allclose([float(x[4]) for x in total], n_rec, rtol=1e-8)
+    np.testing.assert_allclose([float(x[5]) for x in total], m_rec, rtol=1e-8)
```

- [ ] **Step 4: The element model's two columns** in `tests/test_reentry_model_melting.py` and
  `tests/test_reentry_model_coupled.py`:

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -95,7 +95,10 @@
     assert b.material.h_liquid < h_out <= b.material.enthalpy(b.surface_temperature().max())       # own temperature, superheat and all
     assert stats["film_T_max_K"] == pytest.approx(b.surface_temperature().max()) and stats["film_T_max_K"] > b.material.T_liquidus
     rows = np.array(b.source_rows)
-    assert rows.shape[1] == 24 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
+    assert rows.shape[1] == 26 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
+    # sub-plan 22 Task 2: the film is liquid and leaves at its own enthalpy, the removed enthalpy's
+    assert (rows[:, 24] == 0.0).all() and np.isfinite(rows[:, 25]).all()
+    assert (rows[:, 14] * rows[:, 25]).sum() == pytest.approx(b.removed_enthalpy, rel=1e-12)
     # 2026-10-22: viscous_film's floor depth and gamma tau*, finite on Girin & Kopyt's branches and nan elsewhere
     gk = np.isin(rows[:, 8], [1, 2, 4])
     assert np.isfinite(rows[gk, 22:24]).all() and np.isnan(rows[~gk, 22:24]).all() and (rows[gk, 23] < 1.0).all()
```

```diff
--- a/tests/test_reentry_model_coupled.py
+++ b/tests/test_reentry_model_coupled.py
@@ -141,6 +141,8 @@
     from reentry_model import spray
     npz = np.load(files["particles"])
     assert set(npz.files) == set(spray.SOURCE_COLUMNS) and npz["dm_kg"].sum() == pytest.approx(b.sprayed_mass, rel=1e-9)
+    assert np.array_equal(np.column_stack([npz[k] for k in spray.SOURCE_COLUMNS]), np.array(b.source_rows), equal_nan=True)
+    assert not npz["solid_fraction"].any() and np.isfinite(npz["h_kg"]).all()        # liquid film (sub-plan 22 Task 2)
     with open(files["size_distribution"]) as fh:
         lines = fh.read().splitlines()
     assert lines[0] == "window_start_s,window_end_s,r_lo_m,r_hi_m,dn,dm_kg" and (len(lines) - 1) % spray.N_BINS == 0
```

- [ ] **Step 5: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py
  tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_melting.py
  tests/test_reentry_model_coupled.py -q -k "solid_fraction_and_enthalpy or split_by_the_cells or bins_each_flow or
  records_are_one_per_column or close_on_its_books_and_the_coupled_writer or film_spraying_death_and_balances or
  coupled_melting_run_and_writers"`. Expected: 7 failed --
  - `AssertionError: assert (['h_floor_m', 'gamma_tau'] == ['solid_fraction', 'h_kg']` (the source table);
  - `AttributeError: 'ReleaseBook' object has no attribute 'by_phase'` and `TypeError: Cascade.__init__() got an
    unexpected keyword argument 'histogram'` (the release and the flow's book);
  - `TypeError: girin_plan() got an unexpected keyword argument 'points'` and `AttributeError: type object 'SprayPlan'
    has no attribute 'book'` (the plan's records and the body's);
  - `assert (24 == 26)` (the melting body's rows) and `KeyError: 'solid_fraction is not a file in the archive'` (the
    coupled writer's table).

- [ ] **Step 6: The source table's two columns, its blocks and the per-patch properties** in `reentry_model/spray.py`:

```diff
--- a/reentry_model/spray.py
+++ b/reentry_model/spray.py
@@ -95,7 +95,10 @@
 nose by 1/cos(phi)).
 
 Release per patch and step: dm = min(mdot A dt, m_f) (the film after feed and runoff), dn = dm / (4/3 pi rho_l r^3);
-one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX].
+one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX]. The layered
+skins' records are one row per column, per mode, per macro step (`skinbody.SprayPlan.source_block`, sub-plan 22 Task 2):
+the same columns, written by the same `source_block`, the two last columns the droplets' solid fraction and the enthalpy
+per kilogram they carry (spec 8 and 10.2; the 2026-10-07 design 7.5).
 
 The droplet's size caps (2026-10-24, Asha). Every branch's radius is capped by a quarter of the body radius and by one
 of two volumes. Without `cap_depth` (the element model) it is the film on the patch, r <= (3 m_f / (4 pi rho_l))^1/3,
@@ -123,10 +126,13 @@
 BIN_EDGES = np.logspace(np.log10(R_MIN), np.log10(R_MAX), N_BINS + 1)
 WE_BREAKUP = 12.0                  # Pilch-Erdman: secondary breakup expected above this droplet Weber number
 CAPILLARY_TAU = 4.0 * np.pi / (2.0 * np.pi) ** 1.5      # 0.798: tau* = 2 x 2 pi / omega_cap(lambda*)
-# the last two: viscous_film's floor depth and gamma tau* (nan off Girin & Kopyt's branches; 2026-10-22)
+# h_floor_m and gamma_tau: viscous_film's floor depth and gamma tau* (nan off Girin & Kopyt's branches; 2026-10-22);
+# solid_fraction and h_kg: the droplets' solid share of their mass and the enthalpy per kilogram they carry [J/kg, on
+# the material's enthalpy scale] (spec 8 and 10.2, sub-plan 22 Task 2): the element model's film is liquid by
+# construction and booked at its own liquid enthalpy, so its rows carry 0 and that enthalpy (`MeltingBody`)
 SOURCE_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "theta_deg", "x_m", "y_m", "z_m", "closure", "branch", "b_m",
                   "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup",
-                  "h_floor_m", "gamma_tau"]
+                  "h_floor_m", "gamma_tau", "solid_fraction", "h_kg"]
 BRANCH_THICK, BRANCH_THIN, BRANCH_RAREFIED, BRANCH_RT = 0, 1, 2, 3      # Girin thick / thin-film with the edge state /
 #                                                                        thin-film with the freestream / front-surface Rayleigh-Taylor
 BRANCH_REGIME2 = 4                 # deep film whose melt is the more viscous medium: Girin & Kopyt's step profile again
@@ -493,19 +499,45 @@
                            h_floor, gamma_tau, outside & gk, wavelength)
 
 
-def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state):
-    """Source-table rows (list of lists in SOURCE_COLUMNS order) for the patches that released mass this step."""
+def source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None):
+    """Source-table rows, an array (k, len(SOURCE_COLUMNS)), for the patches that released mass (`res.dm` > 0): `h` [m]
+    and `V` [m/s] the altitude and speed, `theta` [rad], `centroids` (n, 3) [m] and `t_hat` (n, 3) per patch, `flow`'s
+    `closure`, `res` a SprayResult (or anything with its per-patch `branch`, `delta_m`, `we_s`, `r`, `dn`, `dm`, `v_s`,
+    `h_floor`, `gamma_tau`), `b` [m] the film thickness, `liquid` the droplets' properties -- each of its `rho`, `mu` and
+    `sigma` one value or one per patch (a slurry droplet's are its slurry's: sub-plan 22 Task 2) -- for their Weber and
+    Ohnesorge numbers at the freestream `state`, and `solid_fraction` and `h_kg` [J/kg] per patch (0 and nan if not given)."""
     k = np.flatnonzero(res.dm > 0.0)
     if k.size == 0:
-        return []
+        return np.zeros((0, len(SOURCE_COLUMNS)))
     r = res.r[k]
-    we_d = state.freestream.rho * state.V ** 2 * 2.0 * r / liquid.sigma
-    oh = liquid.mu / np.sqrt(liquid.rho * liquid.sigma * 2.0 * r)
-    rows = np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
+    n = len(res.dm)
+    rho, mu, sigma = (np.broadcast_to(np.asarray(x, dtype=float), (n,))[k] for x in (liquid.rho, liquid.mu, liquid.sigma))
+    we_d = state.freestream.rho * state.V ** 2 * 2.0 * r / sigma
+    oh = mu / np.sqrt(rho * sigma * 2.0 * r)
+    s_f = np.zeros(k.size) if solid_fraction is None else np.asarray(solid_fraction, dtype=float)[k]
+    h_d = np.full(k.size, np.nan) if h_kg is None else np.asarray(h_kg, dtype=float)[k]
+    return np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
                             centroids[k, 0], centroids[k, 1], centroids[k, 2], flow.closure[k], res.branch[k], b[k],
-                            res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1], t_hat[k, 2],
-                            we_d, oh, (we_d > WE_BREAKUP).astype(float), res.h_floor[k], res.gamma_tau[k]])
-    return rows.tolist()
+                            res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1],
+                            t_hat[k, 2], we_d, oh, (we_d > WE_BREAKUP).astype(float), res.h_floor[k], res.gamma_tau[k],
+                            s_f, h_d])
+
+
+def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None):
+    """Source-table rows (list of lists in SOURCE_COLUMNS order) for the patches that released mass this step
+    (`source_block`'s, as lists: the element model's form, one row per patch-step)."""
+    return source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction, h_kg).tolist()
+
+
+def source_table(rows):
+    """A body's `source_rows` as one array (N, len(SOURCE_COLUMNS)): a list of rows (the element model's, one per
+    patch-step) or of blocks of rows (the layered skins', one array per macro step: a Python list of floats costs about
+    0.9 kB a row against 0.2 kB in an array, and a flight at the continuum step makes of order ten million records)."""
+    if not len(rows):
+        return np.zeros((0, len(SOURCE_COLUMNS)))
+    if isinstance(rows[0], np.ndarray) and rows[0].ndim == 2:
+        return np.concatenate(rows, axis=0)
+    return np.array(rows, dtype=float)
 
 
 def histogram(r, dn, dm):
```

- [ ] **Step 7: The element model's rows** in `reentry_model/body.py`:

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -971,8 +971,9 @@
         if released > 0.0:
             if self.spray_onset is None:
                 self.spray_onset = t
-            if s.particles:
-                self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state)
+            if s.particles:                  # the film is liquid and leaves at its own enthalpy (sub-plan 22 Task 2)
+                self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state,
+                                                          solid_fraction=np.zeros(b.size), h_kg=h_p)
             hn, hm = spray_mod.histogram(res.r, res.dn, res.dm)
             self.hist_n += hn
             self.hist_m += hm
```

- [ ] **Step 8: The release by phase, and the histograms in the flow sub-steps** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -389,13 +389,14 @@
             self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T), U)
         return mass, H
 
-    def take_top_mass(self, dm, admit):
+    def take_top_mass(self, dm, admit, phase=None):
         """Remove dm >= 0 [kg] per column from its top down through the run of admissible cells from its top (`admit`,
         (n, K) bool), at most what that run holds, each cell's metal at its own specific enthalpy and velocity: a partly
         taken cell keeps both, its temperature and its sensitivity, and the cells beneath are untouched (sub-plan 21 Task
         7: the release of spec 8, which thins the top cell until it merges). Returns (mass, enthalpy, solid) taken per
         column, solid the mass of solid in it by the cells' liquid fractions, and leaves the momentum that went in
-        `momentum_out`."""
+        `momentum_out`. With `phase` ((n, K) bool), what was taken is also split by the cells it came from, those in
+        `phase` and the others, and left in `phase_out` (sub-plan 22 Task 2: each spray mode's droplets)."""
         dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
         if not np.all(np.isfinite(dm) & (dm >= 0.0)):
             raise ValueError("a release from a skin's top must be a finite, non-negative mass")
@@ -404,14 +405,36 @@
         taken = np.where(run, np.clip(dm[:, None] - above, 0.0, self.m), 0.0)
         frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
         mass, H = taken.sum(axis=1), (self.H * frac).sum(axis=1)
-        solid = (taken * (1.0 - self.props.liquid_fraction(self.T))).sum(axis=1)
+        f_s = 1.0 - self.props.liquid_fraction(self.T)
+        solid = (taken * f_s).sum(axis=1)
         momentum = (taken[..., None] * self.U).sum(axis=1)
+        if phase is not None:
+            self.phase_out = self._by_phase(taken, frac, f_s, np.asarray(phase, dtype=bool), (mass, H, solid, momentum))
         whole = np.cumprod(run & (taken >= self.m), axis=1).sum(axis=1)        # the cells wholly taken, from the top
         part = (taken > 0.0) & (np.arange(self.m.shape[1])[None, :] >= whole[:, None])
         self.m, self.H = np.where(part, self.m - taken, self.m), np.where(part, self.H * (1.0 - frac), self.H)
         self.take_top(whole)
         self.momentum_out = momentum
         return mass, H, solid
+
+    def _by_phase(self, taken, frac, f_s, phase, totals):
+        """`take_top_mass`'s (mass, enthalpy, solid, momentum) taken, each with a leading axis of two: [0] from the cells
+        in `phase`, [1] from the others. A column that took from cells of one kind gives that kind its totals to the last
+        bit and the other nothing; one that took from both sums each kind's cells."""
+        got = taken > 0.0
+        inside, outside = (got & phase).any(axis=1), (got & ~phase).any(axis=1)
+        out = [np.zeros((2,) + x.shape) for x in totals]
+        for j, only in ((0, inside & ~outside), (1, outside & ~inside)):
+            for o, x in zip(out, totals):
+                o[j, only] = x[only]
+        both = np.flatnonzero(inside & outside)
+        for j, on in enumerate((phase[both], ~phase[both])):
+            t = taken[both] * on
+            parts = (t.sum(axis=1), (self.H[both] * frac[both] * on).sum(axis=1), (t * f_s[both]).sum(axis=1),
+                     (t[..., None] * self.U[both]).sum(axis=1))
+            for o, x in zip(out, parts):
+                o[j, both] = x
+        return tuple(out)
 
     def add_bottom(self, dm, h, s=1.0, u=None):
         """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
@@ -1750,7 +1773,9 @@
     records of its own: what was released is in the pass's FlowBook (`released`, `released_H`, `released_p`,
     `released_solid`). Since sub-plan 22 Task 1c it may give a `Cascade` instead, the liquid's rate and the exposed
     slurry's own per column, which the release takes cell by cell in time order (`cascade_mass`; Girin's,
-    `skinbody.SprayPlan.release`), and the book then carries what the slurry's rate took (`released_slurry`)."""
+    `skinbody.SprayPlan.release`), and the book then carries what the slurry's rate took (`released_slurry`). Every
+    release's book carries it by the cells it came from too, liquid and slurry (`phase_mass`, ...), and a Cascade with a
+    `histogram` has each flow sub-step's droplets binned by size (`hist_n`, `hist_m`; sub-plan 22 Task 2)."""
     columns: object
     dp_ds: np.ndarray
     shear: np.ndarray
@@ -1786,17 +1811,28 @@
     released_solid: np.ndarray = None  # kg (n,): the solid in it, by its cells' liquid fractions
     released_slurry: np.ndarray = None  # kg (n,): what the cascade took at the exposed slurry's own rate (sub-plan 22
                                         # Task 1c; zero for a release given as one rate)
+    # the release by the cells it came from (sub-plan 22 Task 2): [0] fully liquid, at or above the liquidus, [1] the
+    # others (the exposed slurry, only with slurry spraying) -- each spray mode's droplets
+    phase_mass: np.ndarray = None      # kg (2, n)
+    phase_H: np.ndarray = None         # J (2, n): its enthalpy, at its cells' own specific enthalpies
+    phase_p: np.ndarray = None         # kg m/s (2, n, 3): its momentum, at its cells' own velocities
+    phase_solid: np.ndarray = None     # kg (2, n): the solid in it
+    hist_n: np.ndarray = None          # (bins,): the droplets released, binned by size each flow sub-step at the radii of
+    hist_m: np.ndarray = None          # kg (bins,): the Cascade's `histogram`, and their mass (None without one)
 
     @classmethod
     def zero(cls, n):
         return cls(*(np.zeros((n, 3)) for _ in range(4)), released=np.zeros(n), released_H=np.zeros(n),
-                   released_p=np.zeros((n, 3)), released_solid=np.zeros(n), released_slurry=np.zeros(n))
+                   released_p=np.zeros((n, 3)), released_solid=np.zeros(n), released_slurry=np.zeros(n),
+                   phase_mass=np.zeros((2, n)), phase_H=np.zeros((2, n)), phase_p=np.zeros((2, n, 3)),
+                   phase_solid=np.zeros((2, n)))
 
     def add(self, other):
         """Add another book's impulses and releases and keep the worst numbers of both."""
         for key in ("push", "shear", "wall", "carried"):
             setattr(self, key, getattr(self, key) + getattr(other, key))
-        for key in ("released", "released_H", "released_p", "released_solid", "released_slurry"):
+        for key in ("released", "released_H", "released_p", "released_solid", "released_slurry", "phase_mass", "phase_H",
+                    "phase_p", "phase_solid", "hist_n", "hist_m"):
             mine, theirs = getattr(self, key), getattr(other, key)
             if theirs is not None:
                 setattr(self, key, theirs.copy() if mine is None else mine + theirs)
@@ -1885,6 +1921,10 @@
             rb = release_top(sk, dm, flow.slurry)
             step.released, step.released_H = rb.mass, rb.enthalpy
             step.released_p, step.released_solid = rb.momentum, rb.solid
+            ph = rb.by_phase                             # by the cells it came from, and binned by size (Task 2)
+            step.phase_mass, step.phase_H, step.phase_p, step.phase_solid = ph.mass, ph.enthalpy, ph.momentum, ph.solid
+            if isinstance(given, Cascade) and given.histogram is not None:
+                step.hist_n, step.hist_m = given.histogram(ph.mass)
             clock("release", started)
         book.add(step)
     book.cycles = cycles
@@ -1905,6 +1945,8 @@
     enthalpy: np.ndarray     # J (n,), at the cells' own specific enthalpies
     momentum: np.ndarray     # kg m/s (n, 3), at the cells' own velocities: the droplets' momentum, fact 29's sink
     solid: np.ndarray        # kg (n,): the solid in it, by the cells' liquid fractions (the droplets' solid fraction)
+    by_phase: object = None  # the same by the cells it came from, a ReleaseBook of (2, ...) arrays: [0] fully liquid,
+                             # [1] the others (sub-plan 22 Task 2)
 
 
 def release_top(sk, dm, slurry=False):
@@ -1917,11 +1959,14 @@
     enthalpy and velocity (`Skins.take_top_mass`): the droplets carry momentum away from the top cells at their release
     speed (spec 6.4), the speed of the metal they are made of, the top of the resolved profile -- fact 29's sink
     mdot V_s with V_s the sheared layer's surface speed -- so each column's mass, enthalpy and momentum fall by exactly
-    what the book says left. Updates the cells; returns a ReleaseBook."""
+    what the book says left. Updates the cells; returns a ReleaseBook, with what came from the fully liquid cells and
+    what from the others in `by_phase` (sub-plan 22 Task 2: the liquid's mode's droplets and the slurry's)."""
     flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
-    admit = flowing if slurry else flowing & (sk.T >= sk.material.T_liquidus)
-    mass, H, solid = sk.take_top_mass(dm, admit)
-    return ReleaseBook(mass, H, sk.momentum_out.copy(), solid)
+    liquid = flowing & (sk.T >= sk.material.T_liquidus)
+    admit = flowing if slurry else liquid
+    mass, H, solid = sk.take_top_mass(dm, admit, phase=liquid)
+    pm, pH, ps, pp = sk.phase_out
+    return ReleaseBook(mass, H, sk.momentum_out.copy(), solid, ReleaseBook(pm, pH, pp, ps))
 
 
 @dataclass
@@ -1929,10 +1974,13 @@
     """A release given phase by phase (`Flow.release`; sub-plan 22 Task 1c): per column the liquid's rate and the exposed
     slurry's own [kg/(m2 s) of top area, (n,)], the slurry's holding for slurry cells at least as warm as `T_slurry` [K,
     (n,)], the cell whose properties gave it -- a colder slurry cell, more solid and more viscous, ends the column's
-    cascade -- or for every slurry cell where None."""
+    cascade -- or for every slurry cell where None. `histogram`, if given, bins each flow sub-step's droplets by size:
+    a callable of the release by phase [kg, (2, n): the liquid's mode's, the slurry's] giving (droplets, mass [kg]) per
+    size bin (sub-plan 22 Task 2: `skinbody.SprayPlan.histogram`, at the macro step's radii)."""
     rate: np.ndarray
     rate_slurry: np.ndarray
     T_slurry: np.ndarray = None
+    histogram: object = None
 
 
 def cascade_mass(sk, rate, rate_slurry, A_top, dt, slurry=False, T_slurry=None):
```

- [ ] **Step 9: The plan's records, rows and histogram, and the body's** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -79,11 +79,22 @@
 capped by one crest of its wave on the depth the release may take (2026-10-24) -- and, with `slurry_spray`, the exposed
 slurry's own mode by the cascade's properties (spec 8; Asha, 2026-10-29, option i). The plan's two rates per column are
 the release of every flow sub-step of the macro step (`SprayPlan.release`, through `skin.cascade_mass`); what the
-accepted passes released is the plan's record, per column and mode (`spray_plan`)."""
+accepted passes released is the plan's record, per column and mode (`spray_plan`).
+
+The droplet records (sub-plan 22 Task 2; spec 8 and 10.2): what the accepted passes released, by the cells it came from
+-- the liquid's mode's from the fully liquid cells, the slurry's from the others -- is combined per column, per mode, per
+macro step (`SprayPlan.book`), summing mass, enthalpy, solid and momentum over the sub-steps; a column's radius is the
+macro step's, from its plan, so combining changes no radius. At the step's end the records become source-table rows
+(`SprayPlan.source_block`, spray.SOURCE_COLUMNS: the element model's columns, the column's top corner and angle where
+the plan was evaluated, b its flowing depth, the droplets leaving at the mean velocity of the metal they are made of --
+the momentum the books lose with them -- and their solid fraction and enthalpy per kilogram, mass-weighted over the
+cells released), kept one block per macro step in `source_rows`; the size histograms (`hist_n`, `hist_m`) accumulate
+each flow sub-step's release at the plan's radii (`SprayPlan.histogram`, through `skin.Cascade.histogram`)."""
 import copy
 import dataclasses
 import math
 import time
+import types
 
 import numpy as np
 
@@ -169,12 +180,16 @@
     (the cascade; zero without slurry spraying); `result` the liquid's `spray.SprayResult` (branch, radius, wavelength);
     `slurry` the exposed slurry's per column -- its `branch`, radius `r` [m], `wavelength` [m], viscosity `mu` [Pa s] at
     the shear rate `gamma` [1/s] of its own sheared layer, density `rho` [kg/m3], conjugate depth `delta_m` [m], layer
-    `depth` [m] from its top to T_rigid, the viscous film's `gamma_tau` and the cell's temperature `T` [K] (nan, -1 where
-    there is none); the depths read
+    `depth` [m] from its top to T_rigid, the viscous film's `gamma_tau` and floor `h_floor` [m], its `we_s` and the
+    cell's temperature `T` [K] (nan, -1 where there is none); the depths read
     at the start, `flowing`, `liquid` and the cap's `cap_depth` [m], and `on_slurry`, where a flowing cell that is not
-    fully liquid lies beneath the liquid run. The record the macro step's accepted passes fill: `mass` and
-    `mass_slurry` [kg] released by each mode, and, at the step's end, `number` and `number_slurry`, the droplets: the mass
-    over one droplet's at the mode's radius."""
+    fully liquid lies beneath the liquid run. The record the macro step's accepted passes fill (`book`), per column and
+    mode -- the liquid's and the exposed slurry's separate -- from the release by the cells it came from: `mass` and
+    `mass_slurry` [kg], `enthalpy` and `enthalpy_slurry` [J], `solid` and `solid_slurry` [kg] and `momentum` and
+    `momentum_slurry` [kg m/s, (n, 3)]; and, at the step's end, `number` and `number_slurry`, the droplets: the mass over
+    one droplet's at the mode's radius. Where the plan was evaluated, for the records' rows (`source_block`; sub-plan 22
+    Task 2): `theta` [rad], `closure` (the surface flow's), `points` [m, (n, 3)] the top surface's corners, `state` the
+    aero state and `props` the liquid's properties."""
     rate: np.ndarray
     rate_slurry: np.ndarray
     result: object
@@ -187,13 +202,22 @@
     mass_slurry: np.ndarray = None
     number: np.ndarray = None
     number_slurry: np.ndarray = None
+    theta: np.ndarray = None
+    closure: np.ndarray = None
+    points: np.ndarray = None
+    state: object = None
+    props: object = None
 
     def __post_init__(self):
         T = self.slurry["T"]
-        self._cascade = skin_mod.Cascade(self.rate, self.rate_slurry, np.where(np.isfinite(T), T, np.inf))
+        self._cascade = skin_mod.Cascade(self.rate, self.rate_slurry, np.where(np.isfinite(T), T, np.inf),
+                                         histogram=self.histogram)
         n = len(self.rate)
         self.mass, self.mass_slurry = np.zeros(n), np.zeros(n)
         self.number, self.number_slurry = np.zeros(n), np.zeros(n)
+        self.enthalpy, self.enthalpy_slurry = np.zeros(n), np.zeros(n)
+        self.solid, self.solid_slurry = np.zeros(n), np.zeros(n)
+        self.momentum, self.momentum_slurry = np.zeros((n, 3)), np.zeros((n, 3))
 
     def release(self, sk):
         """`skin.Flow.release`: the liquid's rate and the exposed slurry's (`skin.Cascade`), whatever the skins hold within
@@ -201,13 +225,72 @@
         cells at least as warm as the cell that gave it."""
         return self._cascade
 
+    def droplet_mass(self, rho_liquid):
+        """One droplet's mass [kg, (2, n)] at each mode's radius, the liquid's mode's at the density `rho_liquid`, the
+        slurry's at its cell's; nan where the mode has no droplet."""
+        r = np.vstack([self.result.r, self.slurry["r"]])
+        rho = np.vstack([np.broadcast_to(np.asarray(rho_liquid, dtype=float), r[0].shape), self.slurry["rho"]])
+        ok = np.isfinite(r) & (r > 0.0)
+        return np.where(ok, 4.0 / 3.0 * np.pi * rho * np.where(ok, r, 1.0) ** 3, np.nan)
+
     def count(self, rho_liquid):
         """The record's droplets: each mode's mass over one droplet's at its radius."""
-        for mass, r, rho, key in ((self.mass, self.result.r, rho_liquid, "number"),
-                                  (self.mass_slurry, self.slurry["r"], self.slurry["rho"], "number_slurry")):
-            ok = (mass > 0.0) & np.isfinite(r) & (r > 0.0)
-            with np.errstate(invalid="ignore", divide="ignore"):
-                setattr(self, key, np.where(ok, mass / (4.0 / 3.0 * np.pi * rho * np.where(ok, r, 1.0) ** 3), 0.0))
+        one = self.droplet_mass(rho_liquid)
+        for j, (mass, key) in enumerate(((self.mass, "number"), (self.mass_slurry, "number_slurry"))):
+            ok = (mass > 0.0) & np.isfinite(one[j])
+            setattr(self, key, np.where(ok, mass / np.where(ok, one[j], 1.0), 0.0))
+
+    def histogram(self, mass):
+        """(droplets, mass [kg]) per size bin (spray.BIN_EDGES) of a release by mode [kg, (2, n): the liquid's mode's and
+        the slurry's] at the plan's radii -- each flow sub-step's release (`skin.Cascade.histogram`; sub-plan 22 Task 2),
+        which, a column's radius being the macro step's, is the histogram of the step's combined records to round-off."""
+        mass = np.asarray(mass, dtype=float)
+        one = self.droplet_mass(self.props.rho)
+        ok = (mass > 0.0) & np.isfinite(one)
+        number = np.where(ok, mass / np.where(ok, one, 1.0), 0.0)
+        r = np.vstack([self.result.r, self.slurry["r"]])
+        return spray_mod.histogram(np.where(ok, r, np.nan).ravel(), number.ravel(), mass.ravel())
+
+    def book(self, fb):
+        """Add an accepted pass's release (`skin.FlowBook`, by the cells it came from: the liquid's mode took the fully
+        liquid cells, the slurry's the others) to the record."""
+        self.mass += fb.phase_mass[0]
+        self.mass_slurry += fb.phase_mass[1]
+        self.enthalpy += fb.phase_H[0]
+        self.enthalpy_slurry += fb.phase_H[1]
+        self.solid += fb.phase_solid[0]
+        self.solid_slurry += fb.phase_solid[1]
+        self.momentum += fb.phase_p[0]
+        self.momentum_slurry += fb.phase_p[1]
+
+    def source_block(self, t):
+        """The macro step's records as source-table rows (spray.SOURCE_COLUMNS, by `spray.source_block`), one per column
+        and mode that released, the liquid's mode's first: at the end of the step `t` [s] and the plan's aero state, the
+        column's top corner `points` and angle `theta` where the plan was evaluated, its `closure`; the mode's branch,
+        radius, delta_m, We_s, the viscous film's floor and gamma tau*; b, for the liquid's mode the flowing depth and for
+        the slurry's its own layer, the slurry from its top to T_rigid; the droplets (`count` first) and their mass; their
+        speed and direction the mean velocity of the metal released, its momentum over its mass, the momentum the books
+        lose with the droplets (spec 6.4: they leave at the speed of the cells they are made of; zero direction where it
+        is at rest); We_d, Oh and the breakup flag at their radius, by their phase's density, viscosity and surface tension
+        (the slurry's mixture density, viscosity at its sheared layer's rate and SLURRY_SIGMA, the 2026-10-07 design 8.4);
+        and their solid fraction and enthalpy per kilogram, mass-weighted over the cells released."""
+        st, s = self.state, self.slurry
+        flow = types.SimpleNamespace(closure=self.closure)
+        blocks = []
+        for mass, H, solid, p, number, res, b, props in (
+                (self.mass, self.enthalpy, self.solid, self.momentum, self.number, self.result, self.flowing, self.props),
+                (self.mass_slurry, self.enthalpy_slurry, self.solid_slurry, self.momentum_slurry, self.number_slurry,
+                 types.SimpleNamespace(**{k: s[k] for k in ("branch", "r", "delta_m", "we_s", "h_floor", "gamma_tau")}),
+                 s["depth"], material_mod.LiquidProperties(s["rho"], s["mu"], SLURRY_SIGMA))):
+            speed = np.linalg.norm(p, axis=1)
+            on = mass > 0.0
+            per = lambda x: np.divide(x, mass, out=np.zeros_like(x), where=on)
+            direction = np.divide(p, speed[:, None], out=np.zeros_like(p), where=speed[:, None] > 0.0)
+            rows = types.SimpleNamespace(branch=res.branch, r=res.r, delta_m=res.delta_m, we_s=res.we_s, h_floor=res.h_floor,
+                                         gamma_tau=res.gamma_tau, dm=mass, dn=number, v_s=per(speed))
+            blocks.append(spray_mod.source_block(t, st.h, st.V, self.theta, self.points, direction, flow, rows, b, props, st,
+                                                 solid_fraction=per(solid), h_kg=per(H)))
+        return np.concatenate(blocks, axis=0)
 
 
 def region_extent(edges, selected, areas):
@@ -241,7 +324,7 @@
     return out
 
 
-def girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False):
+def girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False, points=None):
     """Girin's spraying from each column's top cells for the macro step about to run (spec 8; sub-plan 22 Task 1c), from
     the surface flow's evaluation on the corners (`flow`), the aero state (`state`, for the rarefied branch), the shared
     spray model (`model`) and the skins as they are at its start: a SprayPlan.
@@ -255,7 +338,8 @@
     its wave on the depth the release may take, the liquid's, or with `slurry` the flowing depth (2026-10-24); a quarter
     of `radius` caps it too. Columns on the windward face, theta <= 90 degrees, only; `A_top` [m2] the top areas, `edges`
     the columns' graph, `dt` [s] the macro step (the nominal release in the result). With `slurry`, the exposed slurry's
-    own mode too (`slurry_modes`)."""
+    own mode too (`slurry_modes`). The plan keeps where it was evaluated for its records' rows: `theta`, the closure,
+    the aero state, the liquid's properties and `points` [m, (n, 3)], the columns' top corners (sub-plan 22 Task 2)."""
     mat, liq = sk.material, model.liquid
     n = sk.n
     F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
@@ -279,12 +363,13 @@
         rate_s, modes = slurry_modes(sk, flow, state, model, windward, w_n, A_top, edges, radius, dt, F, n_liq, n_flow)
     else:
         rate_s, modes = np.zeros(n), _no_slurry(n)
-    return SprayPlan(res.mdot, rate_s, res, modes, F, L, cap, on_slurry)
+    return SprayPlan(res.mdot, rate_s, res, modes, F, L, cap, on_slurry, theta=theta, closure=np.asarray(flow.closure),
+                     points=points, state=state, props=liq)
 
 
 def _no_slurry(n):
     out = {key: np.full(n, np.nan) for key in ("r", "wavelength", "mu", "gamma", "rho", "delta_m", "depth", "gamma_tau",
-                                               "T")}
+                                               "T", "we_s", "h_floor")}
     out["branch"] = np.full(n, -1)
     return out
 
@@ -346,7 +431,7 @@
                        cap_depth=F[cell])
         rate[i] = r.mdot[0]
         out["branch"][i], out["r"][i], out["wavelength"][i] = r.branch[0], r.r[0], r.wavelength[0]
-        out["gamma_tau"][i] = r.gamma_tau[0]
+        out["gamma_tau"][i], out["we_s"][i], out["h_floor"][i] = r.gamma_tau[0], r.we_s[0], r.h_floor[0]
     out["mu"][rows], out["gamma"][rows], out["rho"][rows], out["T"][rows] = mu, gamma, rho, T
     out["delta_m"][rows], out["depth"][rows] = d_m, depth
     return rate, out
@@ -388,12 +473,14 @@
     every flow sub-step (`skin.Flow.release`: a callable of the skins giving a rate [kg/(m2 s)] per unit top area per
     column; sub-plan 21 Task 7), the exposed slurry released too where `slurry_spray`; or `spray`, a `spray.SprayModel`
     with the flow on, Girin's spraying from the top cells evaluated at the start of every macro step (`girin_plan`; sub-plan
-    22 Task 1c), with `slurry_spray` the exposed slurry's own mode too (the cascade), and its record in `spray_plan`."""
+    22 Task 1c), with `slurry_spray` the exposed slurry's own mode too (the cascade), and its record in `spray_plan`; its
+    droplet records, one block of source-table rows per macro step, in `source_rows` where `particles`, and the size
+    histograms in `hist_n` and `hist_m` (sub-plan 22 Task 2: the element model's interface to the coupled writer)."""
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
-                 surface_flow=None, release=None, slurry_spray=False, spray=None):
+                 surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True):
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
         self.flow_cycles = flow_cycles
@@ -410,6 +497,9 @@
             raise ValueError("a skin body sprays by Girin's model or by a release given, not both")
         self.release, self.slurry_spray = release, bool(slurry_spray)
         self.spray, self.spray_plan = spray, None            # Girin's model, and its plan of the last macro step (Task 1c)
+        self.particles = bool(particles)                     # keep the droplet records' rows (Task 2)
+        self.source_rows = []                                # one block of source-table rows per macro step (Task 2)
+        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)   # droplets, kg per size bin
         self._step_release = release                         # what the parts' Flow releases this macro step
         self.repeat_rule = repeat_rule
         self.part_max = part_max                             # s, or None: no cap on the parts while melting
@@ -504,11 +594,12 @@
             self.flow_passes = self.flow_cycles_max = 0
             self.droplet_momentum, self.droplet_shear = np.zeros((self.columns.n, 3)), np.zeros((self.columns.n, 3))
             if self.spray is not None:                       # Girin once, at the macro step's start (sub-plan 22 Task 1c)
-                if end.result is None:
-                    raise ValueError("Girin's spraying reads the surface flow's evaluation: give the aero state, or a "
-                                     "forcing with its result")
+                if end.result is None or state is None:
+                    raise ValueError("Girin's spraying reads the surface flow's evaluation and the aero state (its records' "
+                                     "altitude and speed): give the aero state, with a forcing only if it has its result")
                 self.spray_plan = girin_plan(self.skins, end.result, state, self.spray, self.theta, self.surface.areas,
-                                             self.columns.edges, self.transverse_radius(), dt, self.slurry_spray)
+                                             self.columns.edges, self.transverse_radius(), dt, self.slurry_spray,
+                                             points=self.surface.points.copy())
                 self._step_release = self.spray_plan.release
         sk, mat = self.skins, self.material
         top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
@@ -532,6 +623,10 @@
         if self.spray is not None and self.flow_on:          # the macro step's droplets, from its accepted passes' mass
             self.spray_plan.count(self.spray.liquid.rho)
             self.droplet_number += float(self.spray_plan.number.sum() + self.spray_plan.number_slurry.sum())
+            if self.particles:                               # and its records' rows (Task 2)
+                block = self.spray_plan.source_block(t)
+                if len(block):
+                    self.source_rows.append(block)
 
     def flow_forcing(self, state):
         """The gas side's push on the metal at this macro step (FlowForcing): the surface flow evaluated on the corners
@@ -657,9 +752,11 @@
             self.momentum_lost += fb.released_p.sum(axis=0)  # the droplets' momentum, at their cells' own velocities
             self.droplet_momentum += fb.released_p
             self.droplet_shear += fb.shear
-            if self.spray is not None:                       # the plan's record: each mode's mass (sub-plan 22 Task 1c)
-                self.spray_plan.mass += np.maximum(fb.released - fb.released_slurry, 0.0)
-                self.spray_plan.mass_slurry += fb.released_slurry
+            if self.spray is not None:                       # the plan's record per column and mode (sub-plan 22 Task 1c),
+                self.spray_plan.book(fb)                     # by the cells released, and each flow sub-step's droplets
+                if fb.hist_n is not None:                    # binned by size (Task 2)
+                    self.hist_n += fb.hist_n
+                    self.hist_m += fb.hist_m
         self.iterations.append(res.iterations)
         self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
         self._move_base(book, m_min, g0)
```

- [ ] **Step 10: The writer** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -335,22 +335,30 @@
 
 def write_particles(run_dir, body, history, window=10.0):
     """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
-    size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths."""
+    size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths.
+
+    The table is the body's `source_rows`, the element model's rows or the layered skins' blocks of them
+    (`spray.source_table`; sub-plan 22 Task 2). The summary's columns are the element model's history columns; a history
+    without them -- the skin body's, until sub-plan 22 Task 4 gives it its droplet columns -- writes no summary, and the
+    paths returned then have none."""
     import csv
     from . import spray
     os.makedirs(run_dir, exist_ok=True)
-    rows = np.array(body.source_rows, dtype=float) if body.source_rows else np.zeros((0, len(spray.SOURCE_COLUMNS)))
+    rows = spray.source_table(body.source_rows)
     npz = os.path.join(run_dir, "particles.npz")
     np.savez_compressed(npz, **{name: rows[:, i] for i, name in enumerate(spray.SOURCE_COLUMNS)})
     c = history.columns
     summary = os.path.join(run_dir, "particles_summary.csv")
-    with open(summary, "w", newline="") as fh:
-        w = csv.writer(fh)
-        keys = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg",
-                "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
-        w.writerow(keys)
-        for i in range(len(history)):
-            w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
+    keys = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg",
+            "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
+    if not set(keys) <= set(c):
+        summary = None
+    else:
+        with open(summary, "w", newline="") as fh:
+            w = csv.writer(fh)
+            w.writerow(keys)
+            for i in range(len(history)):
+                w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
     dist = os.path.join(run_dir, "size_distribution.csv")
     edges = spray.BIN_EDGES
     t = rows[:, 0] if len(rows) else np.zeros(0)
@@ -370,4 +378,5 @@
             for j in range(spray.N_BINS):
                 w.writerow(["{:.9g}".format(lo), "{:.9g}".format(hi), "{:.9g}".format(edges[j]), "{:.9g}".format(edges[j + 1]),
                             "{:.9g}".format(n_hist[j]), "{:.9g}".format(m_hist[j])])
-    return {"particles": npz, "particles_summary": summary, "size_distribution": dist}
+    out = {"particles": npz, "particles_summary": summary, "size_distribution": dist}
+    return {k: v for k, v in out.items() if v is not None}
```

- [ ] **Step 11: Run them to see them pass.** Expected: `tests/test_reentry_model_spray.py` 25 passed,
  `tests/test_reentry_model_skin.py` 94, `tests/test_reentry_model_skinbody.py` 34 (2 deselected, the reference tier's;
  the body's writer test about 130 s); the unit tier 538 passed, 1 skipped.

- [ ] **Step 12: The measurements** (in `measure/task2/`, on mains power): `./run_body_records.sh` (Task 1c's code and
  this one's, slurry spraying off and on, and their comparison; about 8 minutes); `./run_bitid.sh` (Task 1c's
  `bitid.py` on this code against `bit_1c.npz`); `"$PY" -B cost_records.py`; and the element model's short case,
  `../task1/cmp/run_case.sh <code> runs line2_short` and `"$PY" ../task1/cmp/compare_runs.py runs line1c_short
  line2_short`, with `line1c_short` copied into `runs/`.

- [ ] **Step 13: Commit.**

```bash
git add reentry_model/spray.py reentry_model/skin.py reentry_model/skinbody.py reentry_model/body.py reentry_model/coupled.py tests/test_reentry_model_spray.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_melting.py tests/test_reentry_model_coupled.py
git commit -m "Skins: the droplet records (spec 8 and 10.2) -- what the accepted passes released, split by the cells it came from, combined per column, per mode (the liquid's and the exposed slurry's), per macro step at the plan's radius, and written by the shared source_block as one block of source-table rows a macro step, the droplets at the mean velocity of the metal released and two new columns, solid_fraction and h_kg, mass-weighted over the cells released; the size histograms accumulated each flow sub-step; coupled.write_particles reads either form; the element model's rows carry 0 and its film's liquid enthalpy, its other columns and history unchanged to the last bit (sub-plan 22 Task 2)"
```

---

### Task 3: Stage 1's skins retired, the skins the default surface model, and the spraying skin flight from the command line

**Files:**
- Delete: `reentry_model/skin1d.py`, `tests/test_reentry_model_skin1d.py`, `tests/test_reentry_model_melting_skin.py`
- Modify: `reentry_model/body.py` (`SURFACE_MODEL_NAMES`, `DEATH_TOLERANCE`, `MeltSettings.surface_model` and `.skin`
  removed; `MeltingBody`'s skins: their counters, `_advance_with_skins` (`advance` is the element step), the branches in
  `melt_step`, `_defer_to_patches`, `_film_nodal`, `_film_and_spray` and `_kill`, the nine skin methods, the
  `surface_temperature` and `radiated_power` overrides, `film_enthalpy`, `film_frozen_fraction`, `mass`, `energy`,
  `demised` and `_skin_stats`), `reentry_model/coupled.py` (`SKIN_COLUMNS` removed, `MELT_COLUMNS` ending in
  `nonrigid_depth_wet_median_um`; eight results and four frame fields removed), `reentry_model/cli.py`
  (`SURFACE_MODELS`, `DEFAULT_SURFACE_MODEL`; `surface_model_of`, `removal_of`, `spray_model_of`, `run_name_of`;
  `model_run_name(..., surface_model="skin")`, `skin_name_suffix(..., slurry_spray=False)`; `--surface-model` and
  `--removal` defaulting to the run's own, `--slurry-spray`; `build_thermal`'s spraying skins; `cmd_run`'s checks and
  the skin run's droplet files), `analysis/reentry_model_thermal_verification.py`, `analysis/melt_verification.py`,
  `analysis/melt_sensitivity.py` (the element model named)
- Test: `tests/test_reentry_model_cli.py` (five tests added, one replaced, five changed),
  `tests/test_reentry_model_melting.py` (one added), `tests/test_reentry_model_thermal_verification.py` (one changed),
  `tests/test_reentry_model_fenicsx.py` (`test_skins_match_the_skfem_backend` removed: the FEniCSx tier is 18)

**Interfaces:**
- Consumes: Task 2's `SkinBody(..., spray=, slurry_spray=, particles=)`, `.source_rows`, `.hist_n`, `.hist_m`,
  `.released_mass`, `.droplet_number` and `coupled.write_particles` (no summary without the element model's history
  columns); `spray.SprayModel`; sub-plan 21's `--skin-flow`, `--skin-flow-cycles` and the continuum step.
- Produces: `cli.SURFACE_MODELS` (`skin`, `elements`), `cli.DEFAULT_SURFACE_MODEL` (`skin`);
  `cli.surface_model_of(args)` (the skins for a thermal run, None without one), `cli.removal_of(args)` (Girin's, but
  instant removal for the heat-only skins); `cli.spray_model_of(args, mat)`; `cli.run_name_of(args, diameter_m,
  velocity_ms, altitude_m, atmosphere_name, dt_continuum)`; `model_run_name(..., surface_model="skin")` ending an
  element-model thermal run in `_surface-elements`; `skin_name_suffix(..., slurry_spray=False)` adding
  `_slurryspray-on`; `--slurry-spray on|off` (off); the skin run's settings `skin.slurry_spray` and, spraying, `rt_spray`,
  `rarefied_shear`, `we_critical`, `k_r`, `k_t`, `gamma_pm`, `kn_body_shock`, `particles`, `liquid`; its files
  `particles` and `size_distribution`, and `particles_summary` null. Removed: `skin1d`, `body.SURFACE_MODEL_NAMES`,
  `body.DEATH_TOLERANCE`, `MeltSettings.surface_model`, `.skin`, `coupled.SKIN_COLUMNS`, and the element model's 22
  history columns and eight results of sub-plan 18's skins. No skin-body history column, frame field or result changes:
  those are Task 4's.

Spec 9 retires the element model's mechanisms from the skin path and keeps `--surface-model elements` unchanged for
comparison; 10.1 makes `skin` the default surface model, `_surface-elements` the suffix when not, and adds
`--slurry-spray on|off` (off, `_slurryspray-on`); 11.2's stage 4 removes "the old film, deep account and their rules".
Since sub-plan 20 the skin path is `skinbody.SkinBody`, which never had the element model's film, deep account, fill
fractions, deaths, feed rule, freeze-back or the stage-1 return device, so what remains to retire is stage 1's
one-dimensional skins of sub-plan 18.

**What used them** (found before anything was removed): `reentry_model/skin1d.py` (503 lines); in `body.py`,
`SURFACE_MODEL_NAMES`, `DEATH_TOLERANCE` and `MeltSettings`'s `surface_model` and `skin`, and in `MeltingBody` the skins
and their thirteen counters, the macro step with skins (`_advance_with_skins`), skin branches in `melt_step` (the step's
own flow evaluation, the draws, the return device, the feed and freeze-back masks, three refreshes of the film's
enthalpy, the death level and the hand-over to skins), `_defer_to_patches`, `_film_nodal`, `_film_and_spray` (the
composed depths, a refresh of the enthalpy, three skin depths), `_kill`, nine methods (`_refresh_skin_index`,
`_skin_available`, `_create_skins`, `_compose_depth`, `_defer_to_facets`, `_return_from_skins`, `_death_to_skins`,
`_hand_over_targets`, `_hand_over_skins`), the skin-aware `surface_temperature`, `radiated_power`, `film_enthalpy`,
`film_frozen_fraction`, `mass`, `energy` and `demised`, and `_skin_stats`; in `coupled.py`, `SKIN_COLUMNS` (23 history
columns), eight results in `melt_results` and four frame fields in `write_vtk_frame`; a comment in `cli.py`; the tests
`test_reentry_model_skin1d.py` (21), `test_reentry_model_melting_skin.py` (14) and the FEniCSx tier's
`test_skins_match_the_skfem_backend`. No analysis script used them. Each branch was dead without skins, so removing it
removes arithmetic the element model never ran: `body.py` falls from 1 671 to 1 251 lines. Kept (ruling 1): the thermal
backends' facet interface law (`step(..., interface=)`, `StepResult.Q_interface`), which only sub-plan 18 called and
which their own tests check, and `nonrigid_depth_wet_median_um`, fact 94's statistic of the regime test's depth over
the wet windward patches, which every element run fills and which now ends `MELT_COLUMNS`. The other 22 of
`SKIN_COLUMNS` and the eight results held 0 or nan in every element run, and go with the skins (ruling 2).

**The skins the default** (spec 10.1). `--surface-model` not given is the skins for a thermal run (`surface_model_of`);
a run without a thermal model has no surface model, and `--surface-model skin` given to such a run is refused as
before. The
run names then follow the rule of 10.1 (`model_run_name`, built from the command line by `run_name_of`):

| Run | Before (Task 2) | After |
|---|---|---|
| `--thermal none` | no surface suffix | no surface suffix |
| thermal, the skins (no melting, or heat-only melting) | `_surface-skin` + the skins' settings | the skins' settings only |
| melting, the skins with the flow | `_melt-girin_surface-skin…_skinflow-on` | `_melt-girin…_skinflow-on` (+ `_flowcycles<n>`, `_slurryspray-on`) |
| thermal, the element model | no surface suffix | `_surface-elements` |
| melting, the element model | `_melt-girin` (and its own suffixes) | `_melt-girin`, its own suffixes, `_surface-elements` |

So a skin run's name loses `_surface-skin` and an element run's gains `_surface-elements`, after the removal's and the
element model's own suffixes and before the seed's. Pinned by
`test_the_run_names_mark_the_element_model_and_the_skins_non_default_settings` and the sub-plan 20 and 21 name tests.
What assumed the old names or the old default in the planning copy, and is changed: three name assertions of
`tests/test_reentry_model_cli.py` (and its element-model tests, which now name `--surface-model elements` in `FEM` and
the physics-mode thermal run); `tests/test_reentry_model_thermal_verification.py`'s command line; and the element
model's verification and sensitivity scripts, `analysis/reentry_model_thermal_verification.py`,
`analysis/melt_verification.py` and `analysis/melt_sensitivity.py`, which name their runs with `--name` but, without the
flag, would have run the skins: each now names its surface model. No document inside the copy names a surface model.
Outside it, unchanged here: the worktree's README, CLAUDE.md and `docs/model_assumptions.md` (no surface model yet:
the package at the repository root has no skins, Task 6), sub-plan 20's level-3 flight plan (it gives `--surface-model
skin`, which still holds), and the wrapper scripts, which drive Step 1 alone.

**The removal follows the flow on the skins** (ruling 3). `--removal` not given is the run's own (`removal_of`):
Girin's for the element model, and for the skins with the flow on, where spec 10.1 has `--skin-flow on` need `--removal
girin` and only spraying removes the liquid; instant removal for the heat-only skins, whose liquid leaves as it forms
(stage 2's stand-in for flow and spray, 11.2). So for the skins Girin's removal is required with the flow, and is the
removal of a flow run, not of every skin melting run: with the flow off by default (10.1, `skinbody.FLOW_DEFAULT`), a
melting run with no flags runs the heat-only skins (`_melt-instant`). A contradiction exits 2: `--removal girin`
without the flow, `--removal instant` with it, `--slurry-spray on` without the flow or with the element model, and the
element model's own settings given to the skins (`--runoff off`, `--deep-runoff off`, `--molten-cascade off`,
`--rigid-substrate off`, `--size-feedback`), which the skins have no use for and would otherwise ignore silently now
that they are the default (ruling 4). The continuum step follows Girin's removal as before (`continuum_step`), so a
spraying skin run takes 0.0125 s once the body's Knudsen number is below 0.01.

**The spraying skin flight.** With `--skin-flow on`, `build_thermal` gives the skin body the shared spray physics with
the command line's constants (`spray_model_of`, the element model's `--kr`, `--kt`, `--we-critical`, `--rt-spray`),
`--slurry-spray` and `--particles`, and records them in the run's settings; and `cmd_run` writes its droplets with the
element model's writer, `coupled.write_particles`: `particles.npz` and `size_distribution.csv`. Its per-step summary,
whose columns are the element model's history columns, waits for Task 4: the run's files record `particles_summary`
as null, and the printed summary says so with the mass, droplets and rows released (ruling 5). The heat-only skins
write no droplet files, having none.

**Measured** (2026-10-29, `prototype/work-2026-10-28-skin-spraying/measure/task3/`; facts 177–179):
- *The element model unchanged* (`element_files.py`; `../task1/cmp/compare_runs.py`; `runs/line3_short` against Task 2's
  `line2_short`). The CLI test's short melting case (fact 169) with `--surface-model elements`: the text of the 93
  history columns kept, in their order, is Task 2's in all 76 rows; `particles.npz`'s 26 fields over its 25 266 rows,
  `particles_summary.csv` and `size_distribution.csv` byte for byte; the settings, and the 61 results kept but the run
  time. The 22 history columns removed held 0 (15) or nan (7) in every row, the eight results 0.
- *The spraying skin flight from the command line* (`run_smoke.sh`, `smoke_check.py`; fact 179). The 100 mm sphere on
  the CLI tests' coarse mesh (2 367 columns, the dev preset's 20 um cells), 125 ms from 65 km at 7.24 km/s (body
  Knudsen number 0.005; q_stag 4.57 MW/m2), the body at 893 K, the flow's sub-steps chosen each skin sub-step, a 25 ms
  first step, then Girin's continuum step of 12.5 ms. The nose's liquid appears in the third step and the last six
  spray. Slurry spraying off: 1.956 g as 76 300 droplets in 1 196 rows, all thick, radius 108/133/1 808 um; the rows'
  mass the released mass exactly, mass 6e-17, energy 4.3e-11, momentum 4.0e-17; `size_distribution.csv` the rows'
  histogram to its nine printed digits; `particles_summary` null and no such file. On: 2.983 g, 915 liquid rows
  (1.246 g) and 1 308 slurry rows (1.737 g, solid fraction up to 0.063, radius median 289 um and largest 7.97 mm); books
  1.3e-16, 4.2e-11, 5.6e-18. The loop's last step is 1.4e-17 s (the steps' sum falls short of t_max), where Pi reads
  5.5 and 8.6: Task 4's. At the CLI tests' 71 km (body Knudsen number 0.011) two 0.25 s steps from 880 K made 0.23 mm
  of liquid over 1.3 mm of flowing metal at the nose and sprayed nothing: off Girin's closure a liquid on slurry is held
  (option B), as in the element model (ruling 7).
- *The cost* (the same runs; `run_cost.sh`). Run alone on mains power, the smoke flight with slurry spraying off took
  61 s (84 s beside the unit tests): 1.5–6.2 s of skins a 12.5 ms macro step once melting (three quarters in the flow;
  12.7 s for the 25 ms first step), 0.09–0.13 s in the 3D solve and 0.46–0.61 s in the moves; peak memory about 740 MB. The same 125 ms at
  62.5 ms macro steps with the continuum step off -- a 0.5 s macro step while melting runs as eight parts of that
  length -- took 116 s: 45.6 and 67.5 s of skins, its passes repeated three times and its flow sub-steps up to ten a
  heat sub-step (the largest Courant number 9.0 against the target 6); it released 1.947 g against 1.956 g (0.5 % less,
  from one plan at 62.5 ms against six), 74 600 droplets against 76 300, in 210 rows against 1 196. So on this case the
  continuum step cost about half as much a second of flight and wrote 5.7 times the rows; measured once, on a coarse
  case, for Task 5's study.
- *The FEniCSx tier* (`fenicsx_tier.txt`, `log_fenicsx_thick_*.txt`): 17 passed and 1 failed,
  `test_the_thick_film_flux_matches_the_skfem_backend`, an element-model energy balance of 1.0496e-8 against the test's
  bound of 1e-8. The same test gives the same residual to the last digit on Tasks 1a, 1b and 2's code and passes on
  sub-plan 21's (602b262): Task 1a's spraying amendments moved it over the bound, and Task 3 leaves it unchanged. Not
  fixed here (the bound is the element model's test's); for the coordinator.
- *A thermal run with no flags on the coarse mesh* (`spike/run_cli.sh`): the CLI tests' physics-mode Step 2 run without
  `--surface-model elements`, now the skins with the default 0.7 mm fine band under 4 mm surface elements, exits 1 at
  the skins' first move ("a corner move of up to 0.0008 m degrades the mesh even after a rebuild"); the skin tests give
  `--fine-size 2`, and the production surface (2 mm) passed level 3 (fact 149).

**Rulings** (2026-10-29; each what, why, cost if wrong):
1. *The thermal backends' facet interface law stays* (`step(..., interface=)`, `StepResult.Q_interface`), though only
   sub-plan 18 called it. Why: it is the backends' own capability, checked by their own tests in both backends, not the
   skins'; the retirement is the brief's three items. Cost: about twenty unused lines in each backend.
2. *The element model's 22 history columns and eight results of sub-plan 18's skins go*, rather than staying as
   constants. Why: they were the retired skins' diagnostics and held 0 or nan in every element run; kept, every element
   history would carry a retired model's names. `nonrigid_depth_wet_median_um`, which every run fills, stays. Cost: an
   element history is 22 columns narrower than Task 2's, so it is the same to the last bit column by column, not as a
   file; nothing in the planning copy reads the columns removed.
3. *On the skins `--removal` not given follows the flow*: instant removal with it off, Girin's with it on. Why: spec 10.1
   pairs the flow with Girin's removal and leaves the flow off by default, so with Girin's removal always the default a
   melting run with no flags would exit 2; this way the default run runs, and if the flow's default changes the removal
   follows it. Cost: a default skin melting run is named `_melt-instant` and runs no droplets while the flow is off.
4. *The element model's settings given to the skins exit 2* (`--runoff off`, `--deep-runoff off`, `--molten-cascade off`,
   `--rigid-substrate off`, any `--size-feedback`). Why: with the skins the default, such a setting given without
   `--surface-model elements` would run the skins and ignore it silently; spec 10.1 refuses contradictory settings. The
   parser cannot tell a default given from none, so only values other than the default are refused. Cost: none; of the
   skins' settings given to the element model only `--skin-flow on` and `--slurry-spray on` are refused, the rest
   ignored as before.
5. *A skin run's summary is recorded as not written* (`particles_summary` null in the run's files, and said in the
   printed summary). Why: Task 2's ruling 10, the summary's columns are the element model's history; Task 4 gives the
   skin history its droplet columns. Cost: none.
6. *The run name is built by `run_name_of`*, taken out of `cmd_run`, so that the rule can be tested from command lines.
   Cost: none (the names of every CLI test are unchanged but for the rule's own change).
7. *The smoke flight's case*: the CLI tests' coarse mesh and the dev preset, 65 km, the body at 893 K, a 25 ms first
   step. Why: the smallest case that melts and sprays within minutes; at the CLI tests' 71 km the body Knudsen number is
   0.011, Girin's closure does not hold, and a liquid on slurry is held from the shear modes (option B), so two 0.25 s
   steps from 880 K sprayed nothing. Cost: not the flight's conditions, which Task 5 flies.

**Open for Asha** (Inputs to settle, below): the flow's default, and the macro step of a spraying skin flight. (Asha
decided on 2026-10-29 that the flow is on by default and that the skins' continuum step also latches at their first
liquid -- the amendment of that day, above, built in Task 3b; the macro step itself waits for Task 5.)

- [ ] **Step 1: The run names, the defaults, the checks and the spraying flight's files** -- in
  `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -99,7 +99,10 @@
 
 
 US76_100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind.csv")
-FEM = BASE + ["--atmosphere", "us76", "--thermal", "fem", "--heating", "sesam", "--h-surface", "4", "--h-core", "20", "--quiet"]
+# a thermal run, with the layered skins, the default surface model since sub-plan 22 Task 3; FEM is the element model, named,
+# and the skins' tests below that start from it give --surface-model skin after it, which overrides it
+THERMAL = BASE + ["--atmosphere", "us76", "--thermal", "fem", "--heating", "sesam", "--h-surface", "4", "--h-core", "20", "--quiet"]
+FEM = THERMAL + ["--surface-model", "elements"]
 
 
 def test_run_name_with_heating():
@@ -158,7 +161,8 @@
 
 def test_thermal_run_physics_mode(tmp_path):
     rc = cli.main(BASE + ["--atmosphere", "us76", "--thermal", "fem", "--heating", "physics", "--h-surface", "4", "--h-core", "20",
-                         "--t-max", "5", "--dt", "0.5", "--outdir", str(tmp_path), "--name", "fem_physics", "--quiet"])
+                         "--t-max", "5", "--dt", "0.5", "--outdir", str(tmp_path), "--name", "fem_physics", "--quiet",
+                         "--surface-model", "elements"])
     assert rc == 0
     doc = json.load(open(tmp_path / "fem_physics.json"))
     s = doc["settings"]
@@ -189,7 +193,7 @@
 
     monkeypatch.setattr(coupled.CoupledRun, "run", run_and_keep)
     argv = FEM + ["--t-max", "5", "--dt", "0.5"]
-    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "sesam")
+    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "sesam", surface_model="elements")
     np.random.seed(1)
     assert cli.main(argv + ["--outdir", str(tmp_path / "a")]) == 0
     np.random.seed(2)
@@ -338,10 +342,10 @@
 def test_skin_settings_reach_the_run_name():
     name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics")
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics",
-                              surface_model="skin") == name + "_surface-skin"
+                              surface_model="skin") == name                    # the default since sub-plan 22 Task 3
     melt = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "instant")
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "instant", surface_model="skin",
-                              skin_suffix="_skin-dev", seed=7) == melt + "_surface-skin_skin-dev_seed-7"
+                              skin_suffix="_skin-dev", seed=7) == melt + "_skin-dev_seed-7"
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", None,
                               surface_model="skin") == cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table",
                                                                           "none", None)        # no thermal model, no skins
@@ -374,12 +378,12 @@
     assert cli.main(argv) == 0
     (path,) = glob.glob(str(tmp_path / "*.json"))
     doc = json.load(open(path))
-    assert doc["run_name"].endswith("_melt-instant_surface-skin_skin-dev_fine4x2mm")
+    assert doc["run_name"].endswith("_fem-sesam_melt-instant_skin-dev_fine4x2mm")
     s, r = doc["settings"], doc["results"]
     assert s["surface_model"] == "skin" and s["prism_layers"] == 0
     assert s["skin"] == {"min_thickness_mm": 0.8, "fine_zone_mm": 0.4, "cell_um": 20.0, "growth": 1.1, "substep_ms": 10.0, "preset": "dev",
                          "deform_depth_mm": 15.0, "n_columns": r["skin_columns"], "fine_band_mm": 4.0, "fine_size_mm": 2.0,
-                         "swell": True, "flow": "off", "flow_cycles": "adaptive"}
+                         "swell": True, "flow": "off", "flow_cycles": "adaptive", "slurry_spray": "off"}
     assert r["skin_columns"] == s["n_patches"] and abs(r["energy_balance_residual"]) < 1e-8
     assert r["held_mass_max_kg"] >= 0.0 and r["held_depth_max_mm"] >= 0.0       # what the fault tolerance held (Task 9)
     rows = list(csv.DictReader(open(path[:-5] + ".csv")))
@@ -400,7 +404,7 @@
     MELT + ["--removal", "instant", "--surface-model", "skin", "--deform-depth", "0"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--fine-band", "-1"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--fine-size", "0"],
-    MELT + ["--surface-model", "skin"],                                     # Girin's removal needs the flow (and spraying)
+    MELT + ["--removal", "girin", "--surface-model", "skin"],               # Girin's spraying needs the flow
     MELT + ["--removal", "instant", "--surface-model", "skin", "--animate"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--prism-layers", "4"],     # the rebuild cannot make prisms
     ["run", "--surface-model", "skin"],                                     # no thermal model
@@ -441,27 +445,51 @@
     assert args.skin_flow == "off" and args.skin_flow_cycles is None
 
 
-def test_a_skin_run_with_flow_keeps_its_liquid_and_writes_the_flow_s_columns_and_frames(tmp_path):
-    """0.5 s from 71 km with a warm body, the layered skins (dev preset) and the flow on at two flow sub-steps a skin
-    sub-step: the liquid stays and flows -- no removal, the spraying that will take it is stage 4's -- the run name ends
-    in _skinflow-on_flowcycles2 after the removal's _melt-girin, the settings record the flow and its sub-cycles, and the
-    history's flow columns and the frames' top speed are written."""
+def test_a_skin_run_with_flow_sprays_from_its_top_cells_and_writes_the_droplet_files(tmp_path):
+    """The spraying skin flight from the command line (sub-plan 22 Task 3). 62.5 ms from 65 km, where the body Knudsen
+    number is 0.005 and Girin's closure holds, of a body at 893 K (2 K below T_rigid) under physics heating, with the
+    layered skins (dev preset) and the flow on at two flow sub-steps a skin sub-step: Girin's removal is the default with
+    the flow, and so is its continuum step, so a 25 ms first step is followed by three of 12.5 ms. The tops melt by the
+    last step's start and its plan sprays them. The run name ends in _melt-girin and the skins' own suffixes, with no
+    surface suffix (the skins are the default); the settings record the flow and the spray model's constants; mass,
+    energy and momentum close; the droplets are written by the element model's writer -- particles.npz, whose rows carry
+    the mass released, and size_distribution.csv, their histogram -- and the per-step summary, which waits for Task 4's
+    history columns, is recorded as not written; the history's flow columns and the frames' top speed are written."""
+    pytest.importorskip("cantera")
     import glob
     import pyvista as pv
-    argv = MELT + ["--altitude", "71", "--velocity", "7.24", "--temperature", "880", "--t-max", "0.5",
-                   "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "1", "--fine-size", "2",
-                   "--skin-flow", "on", "--skin-flow-cycles", "2", "--outdir", str(tmp_path)]
+    from reentry_model import spray
+    argv = MELT + ["--heating", "physics", "--altitude", "65", "--velocity", "7.24", "--temperature", "893", "--dt", "0.025",
+                   "--t-max", "0.0625", "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "1",
+                   "--fine-size", "2", "--skin-flow", "on", "--skin-flow-cycles", "2", "--outdir", str(tmp_path)]
     assert cli.main(argv) == 0
     (path,) = glob.glob(str(tmp_path / "*.json"))
     doc = json.load(open(path))
-    assert doc["run_name"].endswith("_melt-girin_surface-skin_skin-dev_fine4x2mm_skinflow-on_flowcycles2")
-    s, r = doc["settings"], doc["results"]
-    assert s["skin"]["flow"] == "on" and s["skin"]["flow_cycles"] == 2
-    assert r["removed_mass_kg"] == 0.0 and abs(r["energy_balance_residual"]) < 1e-8 and r["momentum_residual"] < 1e-12
+    assert doc["run_name"].endswith("_fem-physics_melt-girin_skin-dev_fine4x2mm_skinflow-on_flowcycles2")
+    s, r, f = doc["settings"], doc["results"], doc["files"]
+    assert s["surface_model"] == "skin" and s["removal"] == "girin" and s["skin"]["flow"] == "on"
+    assert s["skin"]["flow_cycles"] == 2 and s["skin"]["slurry_spray"] == "off" and s["particles"] is True
+    assert s["we_critical"] == 3.08 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["rt_spray"] == "on"
+    assert s["dt_continuum_s"] == 0.0125 and r["dt_switch_time_s"] == pytest.approx(0.025)
     rows = list(csv.DictReader(open(path[:-5] + ".csv")))
-    assert all(k in rows[0] for k in coupled.SKINBODY_COLUMNS)
+    assert len(rows) == 5 and all(k in rows[0] for k in coupled.SKINBODY_COLUMNS)
+    released = r["released_mass_kg"]
+    assert released > 0.0 and float(rows[-2]["released_mass_kg"]) == 0.0 and r["removed_mass_kg"] == 0.0
+    # the books: mass, energy and momentum
+    assert r["initial_mass_kg"] - r["final_mass_kg"] == pytest.approx(released, rel=1e-9)
+    assert abs(r["energy_balance_residual"]) < 1e-8 and r["momentum_residual"] < 1e-12
+    # the droplet files, the element model's, and no summary until Task 4
+    assert f["particles_summary"] is None and not os.path.exists(os.path.join(os.path.dirname(f["particles"]),
+                                                                                  "particles_summary.csv"))
+    table = np.load(f["particles"])
+    assert table.files == spray.SOURCE_COLUMNS and len(table["dm_kg"]) > 0 and (table["time_s"] == 0.0625).all()
+    assert table["dm_kg"].sum() == pytest.approx(released, rel=1e-12)
+    dist = [row for row in csv.DictReader(open(f["size_distribution"])) if float(row["window_end_s"]) == 0.0625]
+    n_hist, m_hist = spray.histogram(table["r_m"], table["dn"], table["dm_kg"])
+    np.testing.assert_allclose([float(row["dn"]) for row in dist], n_hist, rtol=1e-8)
+    np.testing.assert_allclose([float(row["dm_kg"]) for row in dist], m_hist, rtol=1e-8, atol=1e-30)
     surf = sorted(glob.glob(str(tmp_path / "**" / "surface_*.vtp"), recursive=True))
-    assert len(surf) == 2 and "top_speed" in pv.read(surf[-1]).point_data
+    assert len(surf) >= 2 and "top_speed" in pv.read(surf[-1]).point_data
 
 
 @pytest.mark.parametrize("argv", [
@@ -477,3 +505,106 @@
     with pytest.raises(SystemExit) as exc:
         cli.main(argv + ["--outdir", str(tmp_path)])
     assert exc.value.code == 2
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Sub-plan 22 Task 3: the layered skins the default surface model, the run names, the spraying flight's settings
+
+def test_the_skins_are_the_default_and_girin_s_removal_follows_their_flow():
+    """Spec 10.1: --surface-model defaults to the skins for a thermal run (a run without a thermal model has no surface
+    model). --removal not given is the run's own: Girin's for the element model and for the skins with the flow on (spec
+    10.1, --skin-flow on needs --removal girin), instant removal for the heat-only skins (the flow off, its default), whose
+    liquid leaves as it forms; the continuum step (0.0125 s) follows Girin's removal, so a spraying skin run takes it."""
+    p = cli.build_parser()
+    melt = THERMAL + ["--melt", "on"]
+    resolved = lambda argv: (lambda a: (cli.surface_model_of(a), cli.removal_of(a), cli.continuum_step(a)))(p.parse_args(argv))
+    assert cli.DEFAULT_SURFACE_MODEL == "skin" and p.parse_args(THERMAL).surface_model is None
+    assert resolved(BASE) == (None, "girin", None)                               # no thermal model, no surface model
+    assert resolved(THERMAL) == ("skin", "instant", None)                        # nothing melts: the removal is unused
+    assert resolved(melt) == ("skin", "instant", None)                           # the heat-only skins
+    assert resolved(melt + ["--skin-flow", "on"]) == ("skin", "girin", 0.0125)   # the spraying skins
+    assert resolved(melt + ["--surface-model", "elements"]) == ("elements", "girin", 0.0125)
+    assert resolved(melt + ["--surface-model", "elements", "--removal", "instant"]) == ("elements", "instant", None)
+
+
+def test_the_run_names_mark_the_element_model_and_the_skins_non_default_settings():
+    """The run-name rule after the skins became the default (spec 10.1): a thermal run with the element model ends in
+    _surface-elements, after its removal's and its own settings' suffixes and before the seed's; a skin run carries only
+    its non-default settings, the flow's _skinflow-on followed by _flowcycles<n> and _slurryspray-on; a run without a
+    thermal model has no surface suffix. Built from the command line as cmd_run builds it (`run_name_of`)."""
+    p = cli.build_parser()
+    base = "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"
+    name = lambda argv: cli.run_name_of(p.parse_args(argv), 0.1, 7500.0, 77500.133, "us76", cli.continuum_step(p.parse_args(argv)))
+    melt = THERMAL + ["--melt", "on"]
+    assert name(BASE) == base
+    assert name(THERMAL) == base + "_fem-sesam"
+    assert name(FEM) == base + "_fem-sesam_surface-elements"
+    assert name(melt) == base + "_fem-sesam_melt-instant"
+    assert name(melt + ["--skin-preset", "dev"]) == base + "_fem-sesam_melt-instant_skin-dev"
+    assert name(melt + ["--skin-flow", "on"]) == base + "_fem-sesam_melt-girin_skinflow-on"
+    assert name(melt + ["--skin-flow", "on", "--slurry-spray", "on"]) == base + "_fem-sesam_melt-girin_skinflow-on_slurryspray-on"
+    assert (name(melt + ["--skin-flow", "on", "--skin-flow-cycles", "2", "--slurry-spray", "on", "--seed", "7"])
+            == base + "_fem-sesam_melt-girin_skinflow-on_flowcycles2_slurryspray-on_seed-7")
+    assert name(melt + ["--skin-flow", "on", "--dt-continuum", "off"]) == base + "_fem-sesam_melt-girin_dtcontinuum-off_skinflow-on"
+    elements = melt + ["--surface-model", "elements"]
+    assert name(elements) == base + "_fem-sesam_melt-girin_surface-elements"
+    assert name(elements + ["--removal", "instant"]) == base + "_fem-sesam_melt-instant_surface-elements"
+    assert (name(elements + ["--deep-runoff", "off", "--rigid-substrate", "off", "--dt-continuum", "off", "--seed", "7"])
+            == base + "_fem-sesam_melt-girin_deeprunoff-off_rigidsubstrate-off_dtcontinuum-off_surface-elements_seed-7")
+    assert name(elements + ["--skin-preset", "dev"]) == name(elements)          # the skins' settings name skin runs only
+    assert name(melt + ["--name", "mine"]) == "mine"
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, slurry_spray=True) == \
+        "_skinflow-on_slurryspray-on"
+
+
+@pytest.mark.parametrize("argv", [
+    THERMAL + ["--melt", "on", "--removal", "girin"],                                  # Girin's spraying needs the flow
+    THERMAL + ["--melt", "on", "--slurry-spray", "on"],                                # the slurry sprays by the flow's spraying
+    THERMAL + ["--melt", "on", "--skin-flow", "on", "--slurry-spray", "maybe"],
+    FEM + ["--melt", "on", "--slurry-spray", "on"],                                    # the element model has no cascade
+    THERMAL + ["--slurry-spray", "on"],                                                # nothing melts
+    THERMAL + ["--melt", "on", "--runoff", "off"],                                     # the element model's settings ...
+    THERMAL + ["--melt", "on", "--skin-flow", "on", "--deep-runoff", "off"],
+    THERMAL + ["--melt", "on", "--molten-cascade", "off"],
+    THERMAL + ["--melt", "on", "--skin-flow", "on", "--rigid-substrate", "off"],
+    THERMAL + ["--melt", "on", "--size-feedback", "initial"],                         # ... refused with the skins
+    THERMAL + ["--stills"],                                                            # the element model's frames
+])
+def test_contradictory_skin_settings_exit_2(argv, tmp_path):
+    with pytest.raises(SystemExit) as exc:                    # 1 s of flight, were a setting accepted
+        cli.main(argv + ["--t-max", "1", "--outdir", str(tmp_path)])
+    assert exc.value.code == 2
+
+
+def test_an_element_model_setting_given_to_the_default_skins_names_the_surface_model(tmp_path, capsys):
+    """With the skins the default, an element-model setting given without --surface-model elements would otherwise be
+    ignored silently: it is refused, and the message says which surface model it belongs to."""
+    with pytest.raises(SystemExit) as exc:
+        cli.main(THERMAL + ["--melt", "on", "--deep-runoff", "off", "--t-max", "1", "--outdir", str(tmp_path)])
+    assert exc.value.code == 2 and "--surface-model elements" in capsys.readouterr().err
+
+
+class _Built(Exception):
+    pass
+
+
+def test_the_command_line_builds_the_skins_with_the_shared_spray_model(tmp_path, monkeypatch):
+    """--skin-flow on gives the skin body the shared spray physics with the command line's constants (the element model's
+    --kr, --kt, --we-critical, --rt-spray), --slurry-spray and --particles; the heat-only skins get none of them and
+    remove their liquid as it forms. The body is not built: its arguments are."""
+    from reentry_model import skinbody
+
+    def built(*args, **kwargs):
+        raise _Built(kwargs)
+    monkeypatch.setattr(skinbody, "SkinBody", built)
+    melt = THERMAL + ["--melt", "on", "--fine-band", "0", "--t-max", "1", "--outdir", str(tmp_path)]
+    with pytest.raises(_Built) as got:
+        cli.main(melt + ["--skin-flow", "on", "--slurry-spray", "on", "--kr", "0.2", "--kt", "0.9", "--we-critical", "4.62",
+                         "--rt-spray", "off", "--no-particles"])
+    k = got.value.args[0]
+    assert k["flow"] and not k["remove"] and k["surface_flow"] is not None and k["slurry_spray"] and not k["particles"]
+    assert (k["spray"].k_r, k["spray"].k_t, k["spray"].we_critical, k["spray"].rt_spray) == (0.2, 0.9, 4.62, False)
+    with pytest.raises(_Built) as got:
+        cli.main(melt)
+    k = got.value.args[0]
+    assert not k["flow"] and k["remove"] and k["spray"] is None and k["surface_flow"] is None and not k["slurry_spray"]
```

- [ ] **Step 2: Nothing imports the retired skins** -- append to `tests/test_reentry_model_melting.py`; the
  verification script's command line names its surface model, in `tests/test_reentry_model_thermal_verification.py`:

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -792,3 +792,23 @@
     assert b.m_f[i] == 0.0 and b.m_d[i] == 0.0
     assert frozen == pytest.approx(before, rel=1e-15)
     assert (b.phi[owner] - 0.5) * b.element_mass[owner] == pytest.approx(before, rel=1e-9)
+
+
+def test_the_one_dimensional_skins_are_retired():
+    """Sub-plan 22 Task 3 (spec 9): stage 1's one-dimensional skins are gone -- the module, the element model's setting
+    that chose them and its skin columns -- and nothing in the package, its tests or its analysis scripts imports them;
+    the regime test's depth over the wet windward patches (fact 94), which every element run fills, stays."""
+    import dataclasses
+    import importlib.util
+    import pathlib
+    import re
+    from reentry_model import coupled
+    assert importlib.util.find_spec("reentry_model.skin1d") is None and not hasattr(body, "SURFACE_MODEL_NAMES")
+    assert {f.name for f in dataclasses.fields(body.MeltSettings)}.isdisjoint({"surface_model", "skin"})
+    assert not hasattr(coupled, "SKIN_COLUMNS") and not any(k.startswith("skin_") for k in coupled.MELT_COLUMNS)
+    assert coupled.MELT_COLUMNS[-1] == "nonrigid_depth_wet_median_um"
+    imports = re.compile(r"^\s*(from\s+\S*skin1d\S*\s+import|from\s+\S+\s+import\s+[^\n]*\bskin1d\b|import\s+\S*skin1d)", re.M)
+    root = pathlib.Path(__file__).resolve().parent.parent
+    users = [str(p.relative_to(root)) for d in ("reentry_model", "tests", "analysis") for p in sorted((root / d).rglob("*.py"))
+             if imports.search(p.read_text())]
+    assert users == []
```

```diff
--- a/tests/test_reentry_model_thermal_verification.py
+++ b/tests/test_reentry_model_thermal_verification.py
@@ -31,7 +31,8 @@
     assert elements[elements.index("--name") + 1] == "d100__sesam"
     assert skins[skins.index("--name") + 1] == "d100__sesam__skin"
     drop = lambda a: a[:a.index("--name")] + a[a.index("--name") + 2:]
-    assert drop(skins) == drop(elements) + ["--surface-model", "skin"]
+    # each names its surface model, the skins being the model's default since sub-plan 22 Task 3
+    assert elements[-2:] == ["--surface-model", "elements"] and drop(skins) == drop(elements)[:-1] + ["skin"]
     smoke = tv.argv_of("d050", "physics", "out", surface="skin", t_max=2.0)
     assert smoke[smoke.index("--name") + 1] == "d050__physics__skin__t2" and smoke[-2:] == ["--t-max", "2"]
 
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_cli.py
  tests/test_reentry_model_melting.py tests/test_reentry_model_thermal_verification.py -q -k "default_and_girin or
  run_names_mark or contradictory or element_model_setting or builds_the_skins or sprays_from_its_top or
  skin_settings_reach or writes_its_columns_frames or seed_makes or retired or element_run_with_skins"`. Expected: 14
  failed, 7 passed (about 2.5 minutes) --
  - the names: `AssertionError: assert 'model_d100.0..._surface-skin' == 'model_d100.0...e_fem-physics'`, and the two
    skin runs' names ending in `_surface-skin_skin-dev_fine4x2mm...`;
  - `AttributeError: module 'reentry_model.cli' has no attribute 'DEFAULT_SURFACE_MODEL'`, `... 'run_name_of'`;
  - `Failed: DID NOT RAISE SystemExit` (five of `test_contradictory_skin_settings_exit_2`, and the element-model setting
    given to the default: the element model was the default), and `SystemExit: 2` where the skin body was to be built;
  - `AssertionError: assert (ModuleSpec(name='reentry_model.skin1d', ...` (the retirement);
  - `AssertionError: assert (['d100__sesam', '--quiet'] == ['--surface-m...', 'elements']` (the verification script).

- [ ] **Step 4: Remove stage 1's skins.**

```bash
git rm reentry_model/skin1d.py tests/test_reentry_model_skin1d.py tests/test_reentry_model_melting_skin.py
```

  and their test in the FEniCSx tier, `tests/test_reentry_model_fenicsx.py`:

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -308,35 +308,6 @@
     assert np.allclose(Ta, Tb, rtol=0.0, atol=1e-6) and ra.Q_interface == pytest.approx(rb.Q_interface, rel=1e-9)
 
 
-def test_skins_match_the_skfem_backend(coarse_sphere_mesh):
-    """The melt-layer skins (sub-plan 18) in both backends: the body at 800 K under a uniform 1.5 MW/m2 for six steps
-    of 0.5 s, every patch skinned at the first step. Skins, mass, film and temperatures agree and both energy balances
-    are exact."""
-    pytest.importorskip("cantera")
-    from reentry_model import body, heating, mesh
-    from reentry_model import skin1d as skin
-    from test_reentry_model_coupled import MASS_100MM
-    out = {}
-    for name in ("skfem", "fenicsx"):
-        np.random.seed(12345)
-        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
-        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_scheil"), thermal.thermal_solver(name), MASS_100MM,
-                             settings=body.MeltSettings(surface_model="skin1d",
-                                                        skin=skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], substep=0.05)))
-        b.solver.set_temperature(800.0)
-        b.energy0 = b.energy()
-        loads = heating.HeatingResult(np.full(b.surface.n_patches, 1.5e6), 1.5e6, 0.0, 0.0, 0.0)
-        for k in range(6):
-            b.advance(0.5 * (k + 1), 0.5, loads)
-        out[name] = b
-    s, f = out["skfem"], out["fenicsx"]
-    assert s.skins.n == f.skins.n == s.surface.n_patches and np.array_equal(s.mesh.active, f.mesh.active)
-    assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-8) and f.m_f.sum() == pytest.approx(s.m_f.sum(), rel=1e-4)
-    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 1e-3
-    assert np.abs(f.skins.T - s.skins.T).max() < 1e-3
-    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
-
-
 def test_moved_geometry_steps_like_a_fresh_solver(coarse_sphere_mesh):
     """Moving mesh (sub-plan 19): after update_geometry, a step equals a fresh solver's on the moved mesh."""
     from reentry_model import mesh as mesh_mod
```

- [ ] **Step 5: The element model without them** in `reentry_model/body.py`:

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -6,8 +6,6 @@
 from typing import Protocol
 
 import numpy as np
-
-from .thermal import SIGMA_SB
 
 
 class Body(Protocol):
@@ -146,13 +144,12 @@
 PHI_DEATH = 0.05                 # a patch owner below this fraction dies (its remainder goes to the film): keeps the surface
                                  # nodes' thermal mass above 5 % of an element's, which the Newton iteration needs (measured
                                  # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
-DEATH_TOLERANCE = 1.0e-9           # skins drain an element to PHI_DEATH only to rounding: within this it dies (sub-plan 18)
 MAX_CASCADE_PASSES = 32          # safety cap on the molten cascade's passes in one macro step; what it leaves waits for the
                                  # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
 NONRIGID_MAX_CROSSINGS = 64      # elements the non-rigid depth's line may cross before it stops (amendment of 2026-10-06)
 REMOVAL_NAMES = ("girin", "instant")
-SURFACE_MODEL_NAMES = ("elements", "skin1d")  # skin1d: a one-dimensional skin under every melting patch (sub-plan 18, failed
-                                              # its box test; retired in sub-plan 22). The layered skins are skinbody.SkinBody
+# Sub-plan 18's one-dimensional skins under the melting patches (MeltSettings(surface_model="skin1d"), skin1d.py) failed
+# their box test and are retired (sub-plan 22 Task 3, spec 9); the layered skins are skinbody.SkinBody.
 
 
 @dataclass
@@ -171,13 +168,8 @@
     rigid_substrate: bool = True      # the thin branch needs a rigid substrate: the regime test reads the film plus the slurry
                                       # (above 50 % liquid) beneath it, and off Girin's closure a film on slurry does not
                                       # spray (amendment of 2026-10-06)
-    surface_model: str = "elements"   # skin1d: a one-dimensional skin of fine cells under every melting patch resolves the
-                                      # melt layer and feeds the film continuously (sub-plan 18); elements: as before
-    skin: object = None               # skin1d.SkinSettings for surface_model "skin1d"; None: its defaults
 
     def __post_init__(self):
-        if self.surface_model not in SURFACE_MODEL_NAMES:
-            raise ValueError("surface_model must be one of {}, got {!r}".format(SURFACE_MODEL_NAMES, self.surface_model))
         if self.removal not in REMOVAL_NAMES:
             raise ValueError("removal must be one of {}, got {!r}".format(REMOVAL_NAMES, self.removal))
         if self.size_feedback not in SIZE_FEEDBACK_NAMES:
@@ -186,8 +178,6 @@
             raise ValueError("feed_depth must be one of {}, got {!r}".format(FEED_DEPTH_NAMES, self.feed_depth))
         if not 0.0 < self.demise_fraction < 1.0:
             raise ValueError("demise_fraction must be within (0, 1)")
-        if self.surface_model == "skin1d" and self.removal == "instant":
-            raise ValueError("surface_model 'skin1d' needs removal 'girin'")
 
 
 def fit_sphere(points):
@@ -283,22 +273,6 @@
         self.consumed = False
         self._refresh_geometry()
         self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation
-        self.skins = None                                          # sub-plan 18: the skins, with surface_model "skin1d"
-        if self.settings.surface_model == "skin1d" and self.settings.removal != "instant":
-            from . import skin1d as skin_mod
-            self.skins = skin_mod.SkinField(material, self.settings.skin or skin_mod.SkinSettings())
-        self._draw_m, self._draw_H = np.zeros(mesh.n_elements), np.zeros(mesh.n_elements)   # metal the skins drew this step
-        self.skins_created = self.skin_short_events = self.interface_repeats = 0
-        self.skin_drawn_mass = self.skin_handover_mass = self.last_mismatch = 0.0
-        self.skin_frozen_mass = 0.0                                # kg of liquid frozen into skins (spec 5.4, freeze-back)
-        self.skin_returned_mass = 0.0                              # kg of metal skins above 2s returned to their elements
-        self.skin_full_events = 0                                  # skin-steps above 2s whose element had too little room
-        self.skin_thick_feed_events = 0                            # skin-steps that began above 2s and handed liquid on (cumulative)
-        self.last_thick_feeding = (0, float("nan"))                # the step's (count, start-of-step thickness max [m])
-        self.skin_wall = 0.0                                       # s of wall time in the skins' passes (trial and real)
-        self.last_skin_shortfall = 0.0                             # kg of metal the skins lack below their target
-        self._step_flow = None
-        self._refresh_skin_index()
 
     # -- geometry -----------------------------------------------------------------------------------------------
     def _refresh_geometry(self):
@@ -315,7 +289,6 @@
         self.patch_of_face = np.full(len(self.mesh._face_nodes), -1, dtype=np.int64)     # face id -> patch index (-1: not a patch)
         self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
         self._fit_nose()
-        self._refresh_skin_index()
 
     def _fit_nose(self):
         """Mass centre, transverse radius and the windward-cap radius of the current surface (module docstring)."""
@@ -381,8 +354,7 @@
         return 2.0 * self.equivalent_radius()
 
     def mass(self, t):
-        m = float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
-        return m + self.skins.mass() if getattr(self, "skins", None) is not None else m
+        return float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
 
     def reference_area(self):
         """Drag reference area: the hull's projected area, paired with the C_D of drag_shape_factor (the raw surface's
@@ -395,8 +367,7 @@
         return (3.0 * self.mass(0.0) / (4.0 * math.pi * self.material.rho)) ** (1.0 / 3.0)
 
     def energy(self):
-        E = self.solver.energy() + self.film_energy()
-        return E + self.skins.energy() if getattr(self, "skins", None) is not None else E
+        return self.solver.energy() + self.film_energy()
 
     def mean_temperature(self):
         m = self.mass(0.0)
@@ -413,109 +384,6 @@
 
     # -- the step --------------------------------------------------------------------------------------------------
     def advance(self, t, dt, loads, state=None):
-        if self.skins is None:
-            return self._advance_elements(t, dt, loads, state)
-        return self._advance_with_skins(t, dt, loads, state)
-
-    def _advance_with_skins(self, t, dt, loads, state):
-        """The macro step with skins (sub-plan 18; spec sections 4 and 6): create skins; the trial pass gives each skin's
-        base law; the 3D step takes that law on the skins' facets (no heating, no radiation there) with the deferred
-        loads; the real pass, base held at the new facet temperature, feeds, draws and books; the mismatch goes to the 3D
-        nodes, and a step whose mismatch exceeds the tolerance repeats once, re-linearised about that temperature. Then
-        the rest of the melt step. Until the first skin exists this is the element model's step exactly."""
-        import dataclasses
-        import time
-        from . import skin1d as skin_mod
-        from . import spray as spray_mod
-        mat, liq, sk = self.material, self.liquid, self.skins
-        self._create_skins(self.solver.temperature())
-        if not sk.n:
-            self._step_flow = None
-            self.last_thick_feeding = (0, float("nan"))
-            return self._advance_elements(t, dt, loads, state)
-        # the skins above 2s at the start of the step (skin_thick_count's margin): which of them feed this step
-        th0 = sk.thickness().copy()
-        thick0 = th0 > skin_mod.MAX_THICKNESS_FACTOR * sk.settings.thickness * (1.0 + 1e-9)
-        on = self.skin_of_patch >= 0
-        p = self.patch_of_skin
-        flow = delta_m = None
-        if state is not None:                     # evaluated before the passes: the feed's reach needs delta_m
-            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
-            delta_m, _ = spray_mod.melt_layer(flow, liq)
-        self._step_flow = (flow, delta_m)
-        dm_p = delta_m[p] if delta_m is not None else np.full(sk.n, np.nan)
-        room = np.where(np.isfinite(dm_p), mat.rho * sk.area * np.nan_to_num(dm_p), np.inf)
-        avail = self._skin_available(p)
-        L = (self.m_f + self.m_d)[p]
-        q_top = loads.q_conv[p]
-        eps = np.full(sk.n, self.emissivity)
-        T_b0 = self.surface.facet_mean(self.solver.temperature())[p]
-        q3 = np.where(on, 0.0, loads.q_conv)
-        capacity = LOAD_DT_MAX * self.solver.nodal_capacity()
-        applied = np.clip(self.pending_load, -capacity, capacity)
-        solver_state, skin_state = self.solver.state(), sk.state()
-        T_lin = T_b0
-        for attempt in (0, 1):
-            started = time.perf_counter()
-            a, b = sk.advance(dt, q_top, self.T_ambient, eps, T_lin, L, room, avail)
-            self.skin_wall += time.perf_counter() - started
-            A_f, B_f = np.zeros(self.surface.n_patches), np.zeros(self.surface.n_patches)
-            A_f[p], B_f[p] = a, b
-            res = self.solver.step(dt, q3, self.T_ambient, nodal_load=applied / dt, interface=(on, A_f, B_f))
-            T_b1 = self.surface.facet_mean(self.solver.temperature())[p]
-            started = time.perf_counter()
-            book = sk.advance(dt, q_top, self.T_ambient, eps, T_b0, L, room, avail, T_b1=T_b1)
-            self.skin_wall += time.perf_counter() - started
-            mismatch = book.E_base - (a + b * T_b1) * sk.area * dt
-            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
-            if attempt == 0 and (np.abs(mismatch) > skin_mod.INTERFACE_TOLERANCE * scale).any():
-                self.solver.set_state(solver_state)
-                sk.set_state(skin_state)
-                T_lin = T_b1
-                self.interface_repeats += 1
-                continue
-            break
-        # a skin that began the step above 2s and handed liquid to the film or the deep account: the thick skins that melt
-        fed = (book.to_film + book.to_deep) > 0.0
-        sel = thick0 & fed
-        self.last_thick_feeding = (int(sel.sum()), float(th0[sel].max()) if sel.any() else float("nan"))
-        self.skin_thick_feed_events += int(sel.sum())
-        self.last_skin_shortfall = float(np.maximum(mat.rho * sk.area * sk.settings.thickness - sk.m.sum(axis=1), 0.0).sum())
-        self.pending_load = self.pending_load - applied
-        full = np.zeros(self.surface.n_patches)
-        full[p] = mismatch
-        self._defer_to_facets(full)
-        self.last_mismatch = float(np.abs(mismatch).sum())
-        E_top, E_rad = float(book.E_top.sum()), float(book.E_rad.sum())
-        self.integrated_heat += res.Q_conv * dt + E_top
-        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt + E_top - E_rad
-        self.radiated_heat += res.Q_rad * dt + E_rad
-        self.iterations.append(res.iterations)
-        self.last = dataclasses.replace(res, Q_conv=res.Q_conv + E_top / dt, Q_rad=res.Q_rad + E_rad / dt)
-        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
-        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
-        self.m_f[p] = self.m_f[p] + book.to_film
-        self.m_d[p] = self.m_d[p] + book.to_deep
-        # the liquid the sub-steps froze into the skins leaves the patch's accounts (spec 5.4): the deep liquid, which
-        # lies against the solid, first, then the film; both ride the same top node, so the split books nothing
-        from_deep = np.minimum(self.m_d[p], book.frozen)
-        from_film = np.minimum(self.m_f[p], book.frozen - from_deep)
-        self.m_d[p] = np.maximum(self.m_d[p] - from_deep, 0.0)
-        self.m_f[p] = np.maximum(self.m_f[p] - from_film, 0.0)
-        frozen = float(book.frozen.sum())
-        self.skin_frozen_mass += frozen
-        self.frozen_mass += frozen                                 # the model's one re-solidified counter, as with elements
-        e = self.surface.owner[p]
-        np.add.at(self._draw_m, e, book.drawn)
-        np.add.at(self._draw_H, e, book.drawn_H)
-        self.skin_drawn_mass += float(book.drawn.sum())
-        # skin-steps below s/2: a skin thin for ten steps counts ten times
-        self.skin_short_events += int((sk.thickness() < 0.5 * sk.settings.thickness).sum())
-        if book.to_film.sum() + book.to_deep.sum() > 0.0 and self.melt_onset is None:
-            self.melt_onset = t
-        self.melt_step(t, dt, state)
-
-    def _advance_elements(self, t, dt, loads, state=None):
         # A deferred melt load is energy the transferred mass delivered to the nodes it arrived at; a node that has
         # since melted away has nothing to heat with it, so no node is asked to move more than LOAD_DT_MAX in one
         # step and the remainder waits (it is applied as soon as the node has the capacity, dropped into Q_dropped
@@ -558,12 +426,7 @@
         # receded to it. Feeding from any depth let the interior drain through an intact skin -- by 82 s of the 100 mm
         # flight the core was being consumed faster than the outermost shell (measured 2026-09-22).
         flow = delta_m = None
-        if self.skins is not None and self._step_flow is not None:      # with skins: the step's own evaluation
-            flow, delta_m = self._step_flow
-            self._step_flow = None
-            if flow is not None:
-                self.last_face_ids = self.surface.face_ids.copy()
-        elif state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
+        if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
             flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
             delta_m, _ = spray_mod.melt_layer(flow, liq)
             self.last_face_ids = self.surface.face_ids.copy()      # the surface delta_m belongs to (the deaths come later)
@@ -572,12 +435,6 @@
         h_e = h_node[tets].mean(axis=1)                          # as solver.energy() weighs the solid
         h_liq = mat.enthalpy_liquid(T)
         h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
-        if self.skins is not None and self._draw_m.any():        # the metal the skins drew leaves its elements
-            self.phi = self.phi - self._draw_m / self.element_mass
-            self._defer_to_elements(self._draw_m * h_e - self._draw_H)     # the element keeps what the drawn metal left
-            self._draw_m[:], self._draw_H[:] = 0.0, 0.0
-        if self.skins is not None and self.skins.n:              # skins above 2s give their bottom metal back
-            self._return_from_skins(h_e)
         # (i) feed and (iv) re-solidification, netted. The feed hands the molten fraction of what each element still
         # holds to the film (f_feed of phi_e, so an element that is a quarter molten hands over a quarter of its
         # remainder); the mirror rule hands back the fraction of the film that has fallen below the ramp. What stops
@@ -593,17 +450,11 @@
         reach = self._shear_reaches(delta_m) if s.feed_depth == "conjugate" and s.removal != "instant" else None
         if reach is not None:
             f = f * reach
-        if self.skins is not None and self.skins.n:              # a patch's skin feeds its film: its element does not
-            skin_owner = np.zeros(self.mesh.n_elements, dtype=bool)
-            skin_owner[self.surface.owner[self.skin_of_patch >= 0]] = True
-            f = np.where(skin_owner, 0.0, f)
         h_hot = (fn * h_node[tets]).mean(axis=1) * self.mesh.active     # what the molten part carries, per kg of element
         cap = np.where(self.owner_area > 0.0, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
         gross = np.minimum(f * self.phi, cap) * self.element_mass
         # the liquid on each patch below the ramp re-solidifies: the film and the deep liquid beneath it alike
         solid = (1.0 - mat.feed_fraction(self.film_temperature())) * (self.m_f + self.m_d) if s.removal != "instant" else np.zeros(0)
-        if self.skins is not None and self.skins.n and solid.size:      # on a skin it freezes in the skins' sub-steps
-            solid = np.where(self.skin_of_patch >= 0, 0.0, solid)          # (spec 5.4), never into the element
         want = np.bincount(self.surface.owner, solid, self.mesh.n_elements) if solid.size else np.zeros(self.mesh.n_elements)
         net = gross - want
         fed, want = np.maximum(net, 0.0), np.maximum(-net, 0.0)
@@ -629,11 +480,7 @@
             delta, carried = self._add_to_film(fed, fed_h)
             self._defer_to_patches(carried - delta * h_p)        # ... and the melt arrives at its patch's temperature
             if s.deep_runoff and flow is not None:
-                if self.skins is not None:                       # a skin's top moves when energy is booked on it
-                    h_p = self.film_enthalpy(h_liq)
                 deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
-            if self.skins is not None:
-                h_p = self.film_enthalpy(h_liq)
             released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
         frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
         # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed --
@@ -642,10 +489,9 @@
         n_dead, passes, cascade, capped = 0, 0, 0.0, False
         # fully molten as `molten_depth` has it -- the element's mean nodal temperature at or above T_feed -- so the cascade
         # removes exactly the contiguous molten layer the branch test reads, wherever a death exposes it
-        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" and (self.skins is None or not self.skins.n) else None
-        death_level = PHI_DEATH if self.skins is None or not self.skins.n else PHI_DEATH * (1.0 + DEATH_TOLERANCE)
+        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" else None
         while not self.consumed:
-            dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= death_level))
+            dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
             if dead.size == 0:
                 break
             rest = self.phi[dead] * self.element_mass[dead]
@@ -654,15 +500,11 @@
                 self.removed_enthalpy += float(rest.sum()) * mat.h_liquid
                 self._defer_to_elements(rest * (h_e[dead] - mat.h_liquid), dead)
             else:
-                # with skins, an element under a skin passes its remainder to that skin's base (spec section 9)
-                went = self._death_to_skins(dead, rest, h_e) if self.skins is not None and self.skins.n else np.zeros(dead.size, dtype=bool)
                 extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
-                extra[dead[~went]], extra_h[dead[~went]] = rest[~went], rest[~went] * h_e[dead[~went]]
+                extra[dead], extra_h[dead] = rest, rest * h_e[dead]
                 delta, carried = self._add_to_film(extra, extra_h)
                 self._defer_to_patches(carried - delta * h_p)
             self.phi[dead] = 0.0
-            if self.skins is not None:
-                h_p = self.film_enthalpy(h_liq)
             before = float(((self.m_f + self.m_d) * h_p).sum())
             was_owner = self.owner_area > 0.0
             self._kill(dead)                                     # the liquid of a dead patch moves to another patch...
@@ -769,17 +611,8 @@
             self._defer_to_patches(energy * weights / total)
 
     def _defer_to_patches(self, energy):
-        """Book a per-patch energy [J]: on the nodes of a patch without a skin, to be applied over the next step; on a
-        skin's top node at once (sub-plan 18), shared with the liquid that node carries."""
+        """Book a per-patch energy [J] on the nodes of the patches' facets, to be applied over the next step."""
         e = np.asarray(energy, dtype=float)
-        if not e.any():
-            return
-        if getattr(self, "skins", None) is not None and self.skins.n:
-            on = (self.skin_of_patch >= 0) & (e != 0.0)
-            if on.any():
-                p = np.flatnonzero(on)
-                self.skins.book_top(self.skin_of_patch[p], e[p], (self.m_f + self.m_d)[p])
-                e = np.where(self.skin_of_patch >= 0, 0.0, e)
         if e.any():
             np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))
 
@@ -790,8 +623,6 @@
         if not self.m_f.size:
             return np.zeros(self.mesh.n_nodes)
         lm = self.m_f + self.m_d
-        if getattr(self, "skins", None) is not None:              # liquid on a skin rides the skin's top, not the nodes
-            lm = np.where(self.skin_of_patch >= 0, 0.0, lm)
         return np.bincount(self.surface.faces.ravel(), np.repeat(lm / 3.0, 3), self.mesh.n_nodes)
 
     def _deep_runoff(self, dt, flow, delta_m, reach, T, h_node, h_e, h_p):
@@ -905,13 +736,11 @@
         return self.m_f - before, got
 
     def _film_and_spray(self, t, dt, state, h_p, flow=None, delta_m=None):
-        liq, s, mat = self.liquid, self.settings, self.material
+        liq, s = self.liquid, self.settings
         from . import spray as spray_mod
         if state is None or self.m_f.sum() <= 0.0 or flow is None:
             self.last_flow = self.last_spray = self.last_nonrigid = None
-            nan = float("nan")
-            self.last_melt.update({"liquid_depth_mean_um": nan, "liquid_depth_max_um": nan,
-                                   "nonrigid_depth_skin_mean_um": nan, "nonrigid_depth_wet_median_um": nan})
+            self.last_melt.update({"nonrigid_depth_wet_median_um": float("nan")})
             return 0.0
         areas = self.surface.areas
         # the depth of liquid under each patch: its film plus the contiguous molten material beneath it. The branch
@@ -924,10 +753,6 @@
         # instead of the fully molten material alone. The Rayleigh-Taylor criteria and the molten region the wave-fits
         # test measures against keep the liquid.
         nonrigid = self.nonrigid_depth() if s.rigid_substrate else None
-        if self.skins is not None and self.skins.n:              # the skins resolve both depths (spec section 5.5)
-            molten = self._compose_depth(mat.T_feed, molten)
-            if nonrigid is not None:
-                nonrigid = self._compose_depth(mat.T_rigid, nonrigid)
         under = molten if nonrigid is None else nonrigid
         # (ii) lubrication and runoff
         n_sub, moved = 0, 0.0
@@ -941,8 +766,6 @@
             # over the patches that gained film, which is exact in total and second order in the attribution)
             d = self.m_f - before
             self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
-            if self.skins is not None and self.skins.n:          # the booking moved skin tops: the droplets leave at theirs
-                h_p = self.film_enthalpy()
         b = self.m_f / (liq.rho * areas)
         layer = b + molten + deep
         regime_layer = layer if nonrigid is None else b + nonrigid + deep
@@ -998,16 +821,8 @@
             slurry_thick = float((thick & ~liquid_thick)[girin_wet].mean())
             rigid_thin = float((liquid_thick & ~thick)[girin_wet].mean())
         nonrigid_mean = float(nonrigid[wet_w].mean() * 1e3) if nonrigid is not None and wet_w.any() else float("nan")
-        # the resolved depths, at evaluation time on the surface they were evaluated on (final review ruling): over the
-        # skin patches the composed skin + 3D depths [um]; for every run the regime test's depth over the wet windward
-        # patches (fact 94's statistic)
-        nan = float("nan")
-        skinned = self.skin_of_patch >= 0 if self.skins is not None and self.skins.n else np.zeros(molten.shape, dtype=bool)
-        depths = {"liquid_depth_mean_um": float(molten[skinned].mean() * 1e6) if skinned.any() else nan,
-                  "liquid_depth_max_um": float(molten[skinned].max() * 1e6) if skinned.any() else nan,
-                  "nonrigid_depth_skin_mean_um": float(nonrigid[skinned].mean() * 1e6)
-                  if nonrigid is not None and skinned.any() else nan,
-                  "nonrigid_depth_wet_median_um": float(np.median(under[wet_w]) * 1e6) if wet_w.any() else nan}
+        # the regime test's depth over the wet windward patches (fact 94's statistic), at evaluation time
+        depths = {"nonrigid_depth_wet_median_um": float(np.median(under[wet_w]) * 1e6) if wet_w.any() else float("nan")}
         # The Rayleigh-Taylor criterion is reported, never applied, so report it usefully: a body-level "any patch"
         # boolean says nothing about how much melt is involved or whether the unstable wave even fits on the nose.
         # `rt_mass_fraction` is the share of the film on unstable patches; `rt_wavelength_over_nose` is the shortest
@@ -1059,20 +874,13 @@
             self.consumed = True                                 # the liquid still on it leaves with the body
             self.removed_mass += float(self.m_f.sum() + self.m_d.sum())
             self.removed_enthalpy += float(((self.m_f + self.m_d) * self.film_enthalpy()).sum())
-            if self.skins is not None and self.skins.n:          # ... and so does the skins' metal
-                self.removed_mass += self.skins.mass()
-                self.removed_enthalpy += self.skins.energy()
-                self.skins.keep(np.zeros(0, dtype=np.int64))
             self.m_f, self.m_d = np.zeros(0), np.zeros(0)
-            self.skin_of_patch, self.patch_of_skin = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
             return
         self._refresh_geometry()
         keep = self.patch_of_face[old_surface.face_ids]
         kept = keep >= 0
         self.m_f = self._hand_over(old_surface, old_m_f, keep, kept)
         self.m_d = self._hand_over(old_surface, old_m_d, keep, kept)   # the deep liquid beneath it goes the same way
-        if self.skins is not None and self.skins.n:
-            self._hand_over_skins(old_surface)
 
     def _hand_over(self, old_surface, old, keep, kept):
         """A per-patch liquid mass of the surface before the deaths, carried onto the current surface: kept patches keep
@@ -1113,177 +921,6 @@
                 np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old[orphan])
         return m_f
 
-    # -- skins (sub-plan 18) -----------------------------------------------------------------------------------------
-    def _refresh_skin_index(self):
-        """skin_of_patch: the skin row under each patch of the current surface (-1: none); patch_of_skin the reverse
-        (-1 for a skin whose patch has vanished, until the death loop hands it over)."""
-        self.skin_of_patch = np.full(self.surface.n_patches, -1, dtype=np.int64)
-        sk = getattr(self, "skins", None)
-        if sk is None or not sk.n:
-            self.patch_of_skin = np.zeros(0, dtype=np.int64)
-            return
-        self.patch_of_skin = self.patch_of_face[sk.face_id]
-        ok = self.patch_of_skin >= 0
-        self.skin_of_patch[self.patch_of_skin[ok]] = np.flatnonzero(ok)
-
-    def _skin_available(self, patches):
-        """Metal [kg] a skin on each of `patches` may take from its element: its share, by area among the element's
-        patches, of what the element holds above PHI_DEATH."""
-        e = self.surface.owner[patches]
-        share = self.surface.areas[patches] / self.owner_area[e]
-        return np.maximum(self.phi[e] - PHI_DEATH, 0.0) * self.element_mass[e] * share
-
-    def _create_skins(self, T):
-        """Patches whose surface has reached the solidus, or that hold liquid, get a skin (spec 5.2): a slab of the
-        target thickness -- less where the element cannot give that above PHI_DEATH -- taken from the patch's element,
-        its cells filled from the P1 field along the inward normal. The element's mass leaves at its own mean nodal
-        enthalpy and the difference to what the cells hold is booked on its nodes; liquid already on the patch moves
-        onto the skin's top with its energy. Returns how many skins were made."""
-        sk, s, mat = self.skins, self.surface, self.material
-        if sk is None or not s.n_patches:
-            return 0
-        Tf = s.facet_mean(T)
-        idx = np.flatnonzero((self.skin_of_patch < 0) & ((Tf >= mat.T_solidus) | (self.m_f + self.m_d > 0.0)))
-        if not idx.size:
-            return 0
-        thick = np.minimum(sk.settings.thickness, self._skin_available(idx) / (mat.rho * s.areas[idx]))
-        idx, thick = idx[thick > 0.0], thick[thick > 0.0]
-        if not idx.size:
-            return 0
-        e = s.owner[idx]
-        z = (np.arange(sk.nc) + 0.5)[None, :] * (thick / sk.nc)[:, None]
-        x = s.centroids[idx][:, None, :] - s.normals[idx][:, None, :] * z[:, :, None]
-        P = self.mesh.points[self.mesh.tets[e]]
-        inv = np.linalg.inv(np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]], axis=2))
-        lam = np.einsum("kij,kcj->kci", inv, x - P[:, None, 0])
-        lam = np.concatenate([1.0 - lam.sum(axis=2, keepdims=True), lam], axis=2)
-        lam = np.clip(lam, 0.0, None)
-        lam = lam / lam.sum(axis=2, keepdims=True)            # a point past the element takes its nearest face's value
-        Tp = np.einsum("kci,ki->kc", lam, T[self.mesh.tets[e]])
-        L = (self.m_f + self.m_d)[idx]
-        h_nodes = self.film_enthalpy()[idx]                   # the liquid's enthalpy on the 3D nodes, before the skins
-        mass, H = sk.add(s.face_ids[idx], s.areas[idx], Tp, thick)
-        h_e = mat.enthalpy(T)[self.mesh.tets[e]].mean(axis=1)
-        np.subtract.at(self.phi, e, mass / self.element_mass[e])
-        self._defer_to_elements(np.bincount(e, mass * h_e - H, self.mesh.n_elements))
-        self._refresh_skin_index()
-        if (L > 0.0).any():
-            rows = self.skin_of_patch[idx]
-            sk.book_top(rows, L * (h_nodes - mat.enthalpy_liquid(sk.T[rows, 0])), L)
-        self.solver.set_fractions(self.phi)
-        self.solver.set_film_mass(self._film_nodal())
-        self.skins_created += int(idx.size)
-        return int(idx.size)
-
-    def _compose_depth(self, level, depth3d):
-        """Per patch: the depth of the layer above `level` read from the skin where there is one, continued into the 3D
-        model's own depth (measured from the skin's base) when the whole skin is above it (spec 5.5)."""
-        out = np.array(depth3d, dtype=float, copy=True)
-        if self.skins is None or not self.skins.n:
-            return out
-        d, through = self.skins.depth_above(level)
-        on = self.skin_of_patch >= 0
-        r = self.skin_of_patch[on]
-        out[on] = d[r] + np.where(through[r], out[on], 0.0)
-        return out
-
-    def _defer_to_facets(self, energy):
-        """Book a per-patch energy [J] on the 3D nodes of the patches' facets, skins or not (the interface mismatch)."""
-        e = np.asarray(energy, dtype=float)
-        if e.any():
-            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))
-
-    def _return_from_skins(self, h_e):
-        """The resolution device (spec 5.4, "Thickness bounds"): a skin above MAX_THICKNESS_FACTOR times its target
-        thickness has cells too coarse for the melt layer, so it returns its bottom metal in excess of that to its
-        element where the element has room (phi_e < 1), the room shared among the element's skin patches by area. The
-        metal leaves at the enthalpy of the bottom cells and arrives in the element at its mean nodal enthalpy h_e [J/kg];
-        the difference is booked to the element's nodes, as with the draw. Where the element is full the skin keeps the
-        metal and the event is counted. It moves metal between two accounts of the same solid and changes neither how
-        much freezes nor how much melts."""
-        from .skin1d import MAX_THICKNESS_FACTOR
-        sk, s = self.skins, self.surface
-        p = np.flatnonzero(self.skin_of_patch >= 0)
-        r, e = self.skin_of_patch[p], s.owner[p]
-        bound = self.material.rho * sk.area[r] * MAX_THICKNESS_FACTOR * sk.settings.thickness
-        excess = sk.m[r].sum(axis=1) - bound
-        # a skin returned to the bound sits at it only to rounding (up to about 4e-16 of it above): an excess below 1e-9
-        # of the bound, skin_thick_count's margin, is none -- else the next step returns 1e-21 kg, or counts a full
-        # event on a full element. What is given is still measured from the bound itself.
-        excess = np.where(excess > 1e-9 * bound, excess, 0.0)
-        if not excess.any():
-            return
-        skin_area = np.bincount(e, s.areas[p], self.mesh.n_elements)
-        room = np.maximum(1.0 - self.phi[e], 0.0) * self.element_mass[e] * s.areas[p] / skin_area[e]
-        give = np.minimum(excess, room)
-        self.skin_full_events += int(((excess > 0.0) & (give < excess)).sum())
-        some = give > 0.0
-        if not some.any():
-            return
-        p, r, e, give = p[some], r[some], e[some], give[some]
-        H = sk.take_bottom(r, give, (self.m_f + self.m_d)[p])
-        # a new array, as with the draw (the solver may hold the old one until set_fractions); the shares of an
-        # element's room add up to it only to rounding, so phi is held at 1
-        self.phi = np.minimum(self.phi + np.bincount(e, give / self.element_mass[e], self.mesh.n_elements), 1.0)
-        self._defer_to_elements(np.bincount(e, H - give * h_e[e], self.mesh.n_elements))
-        self.skin_returned_mass += float(give.sum())
-
-    def _death_to_skins(self, dead, rest, h_e):
-        """Dying elements under skins pass their remainder `rest` [kg] into those skins' bottoms, shared by area among
-        the element's skin patches, at the element's own mean nodal enthalpy (spec 9). Returns which of `dead` did."""
-        on = self.skin_of_patch >= 0
-        p = np.flatnonzero(on & np.isin(self.surface.owner, dead))
-        went = np.isin(dead, self.surface.owner[p])
-        if p.size:
-            e = self.surface.owner[p]
-            skin_area = np.bincount(e, self.surface.areas[p], self.mesh.n_elements)
-            full = np.zeros(self.mesh.n_elements)
-            full[dead] = rest
-            mass = full[e] * self.surface.areas[p] / skin_area[e]
-            self.skins.add_bottom(self.skin_of_patch[p], mass, mass * h_e[e], (self.m_f + self.m_d)[p])
-        return went
-
-    def _hand_over_targets(self, old_surface, idx):
-        """For vanished patches `idx` of `old_surface`: (targets (k, 4), patch indices of the current surface or -1;
-        weights (k, 4)) -- the faces their own element exposed, weighted by area projected on the vanished patch's
-        normal (raw area where none faces that way), or the nearest patch where the element exposed nothing. The rule
-        of `_hand_over`, kept apart so that the liquid's hand-over keeps its exact arithmetic."""
-        targets = self.patch_of_face[self.mesh.element_faces(old_surface.owner[idx])]
-        ok = targets >= 0
-        safe = np.where(ok, targets, 0)
-        cosine = np.clip(np.einsum("nij,nj->ni", self.surface.normals[safe], old_surface.normals[idx]), 0.0, None)
-        area = np.where(ok, self.surface.areas[safe] * cosine, 0.0)
-        total = area.sum(axis=1)
-        flat = ok & (total[:, None] <= 0.0)
-        if flat.any():
-            area = np.where(flat, self.surface.areas[safe], area)
-            total = area.sum(axis=1)
-        has = total > 0.0
-        weights = np.where(has[:, None], area / np.where(has, total, 1.0)[:, None], 0.0)
-        targets = np.where(ok & has[:, None], targets, -1)
-        if (~has).any():
-            orphan = np.flatnonzero(~has)
-            targets[orphan, 0] = self._patch_tree.query(old_surface.centroids[idx[orphan]])[1]
-            weights[orphan, 0] = 1.0
-        return targets, weights
-
-    def _hand_over_skins(self, old_surface):
-        """Skins whose patch vanished pass to the faces the deaths exposed, layer by layer (spec 9); empty skins go."""
-        sk = self.skins
-        lost = np.flatnonzero(self.patch_of_face[sk.face_id] < 0)
-        if lost.size:
-            old_patch = np.full(len(self.patch_of_face), -1, dtype=np.int64)
-            old_patch[old_surface.face_ids] = np.arange(old_surface.n_patches)
-            targets, weights = self._hand_over_targets(old_surface, old_patch[sk.face_id[lost]])
-            tgt_faces = np.where(targets >= 0, self.surface.face_ids[np.maximum(targets, 0)], -1)
-            _, handed = sk.hand_over(lost, tgt_faces, weights, lambda f: self.surface.areas[self.patch_of_face[f]])
-            self.skin_handover_mass += handed
-        keep = np.ones(sk.n, dtype=bool)
-        keep[lost] = False
-        keep &= sk.m.sum(axis=1) > 0.0
-        sk.keep(np.flatnonzero(keep))
-        self._refresh_skin_index()
-
     # -- reporting -----------------------------------------------------------------------------------------------
     def on_current_surface(self, values, fill):
         """A per-patch array from the last flow and spray evaluation (`last_flow`, `last_spray`), carried onto the
@@ -1300,27 +937,6 @@
         j = index[self.surface.face_ids]
         return np.where(j >= 0, np.asarray(values, dtype=float)[np.maximum(j, 0)], fill)
 
-    def surface_temperature(self):
-        """Wall temperature per patch: the skin's top cell where a patch has a skin (sub-plan 18), else the facet's."""
-        Tf = self.surface.facet_mean(self.solver.temperature())
-        if getattr(self, "skins", None) is not None and self.skins.n:
-            on = self.skin_of_patch >= 0
-            Tf[on] = self.skins.T[self.skin_of_patch[on], 0]
-        return Tf
-
-    def radiated_power(self):
-        """Radiated power [W]. With skins (sub-plan 18) the facets under skins do not radiate -- their skins' tops do,
-        at patch area -- read from the current skin index rather than from the solver's mask of its last step, which
-        the surface refresh at the end of every melt step clears."""
-        sk = getattr(self, "skins", None)
-        if sk is None or not sk.n:
-            return self.solver.radiated_power(self.T_ambient)
-        Tf = self.surface.facet_mean(self.solver.temperature())
-        bare = self.skin_of_patch < 0
-        P = float((self.emissivity * SIGMA_SB * (Tf[bare] ** 4 - self.T_ambient ** 4) * self.surface.areas[bare]).sum())
-        T0 = sk.T[:, 0]
-        return P + float((self.emissivity * SIGMA_SB * (T0 ** 4 - self.T_ambient ** 4) * sk.area).sum())
-
     def film_temperature(self):
         """The film's temperature per patch: the surface's own (the film is thermally thin, class docstring)."""
         return self.surface_temperature()
@@ -1335,11 +951,7 @@
             h_node = self.material.enthalpy_liquid(self.solver.temperature())
         if not self.m_f.size:
             return np.zeros(0)
-        h = h_node[self.surface.faces].mean(axis=1)
-        if getattr(self, "skins", None) is not None and self.skins.n:   # liquid on a skin shares its top cell's temperature
-            on = self.skin_of_patch >= 0
-            h[on] = self.material.enthalpy_liquid(self.skins.T[self.skin_of_patch[on], 0])
-        return h
+        return h_node[self.surface.faces].mean(axis=1)
 
     def film_energy(self):
         """Enthalpy held by the film [J]. Identical to the nodal sum the solver carries (`set_film_mass` gives each
@@ -1350,13 +962,8 @@
     def film_frozen_fraction(self):
         """Fraction of the film sitting on patches below the feed ramp: mass the enthalpy calls solid but the model
         still treats as liquid, because the mesh cannot grow a crust outside itself and its own element may be gone
-        (section 9.5). It is the honest measure of that limitation, so it is recorded every step. Skin patches are left
-        out of both sums (sub-plan 18): their liquid freezes in the skins' sub-steps as far as the top node's energy
-        allows, and what is left waits at T_low -- where the feed fraction is zero -- to freeze in the next ones, so it
-        is not film the model cannot freeze."""
+        (section 9.5). It is the honest measure of that limitation, so it is recorded every step."""
         m_f = self.m_f
-        if getattr(self, "skins", None) is not None and self.skins.n:
-            m_f = np.where(self.skin_of_patch >= 0, 0.0, m_f)
         total = float(m_f.sum())
         if total <= 0.0:
             return 0.0
@@ -1577,7 +1184,7 @@
 
     def demised(self):
         """True once the body's material (film excluded) is below the demise fraction of the initial mass."""
-        metal = float((self.phi * self.element_mass).sum()) + (self.skins.mass() if self.skins is not None else 0.0)
+        metal = float((self.phi * self.element_mass).sum())
         return self.consumed or metal < self.settings.demise_fraction * self.mass0
 
     def energy_balance_residual(self):
@@ -1640,32 +1247,5 @@
                 "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                 "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
                 "fitted_nose_radius_mm": self.fitted_nose_radius * 1e3,
-                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active), **self._skin_stats()}
-
-    def _skin_stats(self):
-        """The skins' history columns (sub-plan 18); zero or nan without skins. The resolved depths are the step's, taken
-        in `_film_and_spray` on the surface they were evaluated on."""
-        from .skin1d import MAX_THICKNESS_FACTOR
-        sk, lm, nan = self.skins, self.last_melt, float("nan")
-        th = sk.thickness() if sk is not None and sk.n else np.zeros(0)
-        return {"skin_count": float(sk.n) if sk is not None else 0.0, "skin_mass_kg": sk.mass() if sk is not None else 0.0,
-                "skin_thickness_min_mm": float(th.min() * 1e3) if th.size else float("nan"),
-                "skin_thickness_median_mm": float(np.median(th) * 1e3) if th.size else float("nan"),
-                "skin_drawn_mass_kg": self.skin_drawn_mass, "skin_handover_mass_kg": self.skin_handover_mass,
-                "skin_short_events": float(self.skin_short_events), "interface_mismatch_J": self.last_mismatch,
-                "interface_repeats": float(self.interface_repeats),
-                "liquid_depth_mean_um": lm.get("liquid_depth_mean_um", nan),
-                "liquid_depth_max_um": lm.get("liquid_depth_max_um", nan),
-                "nonrigid_depth_skin_mean_um": lm.get("nonrigid_depth_skin_mean_um", nan),
-                "skin_wall_s": self.skin_wall,
-                "skin_thickness_max_mm": float(th.max() * 1e3) if th.size else nan,
-                # beyond rounding: a skin that returned its metal down to 2s sits at it to an ulp or two either side
-                "skin_thick_count": float((th > MAX_THICKNESS_FACTOR * sk.settings.thickness * (1.0 + 1e-9)).sum())
-                if th.size else 0.0,
-                "skin_shortfall_kg": self.last_skin_shortfall,
-                "nonrigid_depth_wet_median_um": lm.get("nonrigid_depth_wet_median_um", nan),
-                "skin_frozen_mass_kg": self.skin_frozen_mass, "skin_returned_mass_kg": self.skin_returned_mass,
-                "skin_full_events": float(self.skin_full_events),
-                "skin_thick_feeding_count": float(self.last_thick_feeding[0]),
-                "skin_thick_feeding_max_mm": self.last_thick_feeding[1] * 1e3,
-                "skin_thick_feed_events": float(self.skin_thick_feed_events)}
+                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active),
+                "nonrigid_depth_wet_median_um": lm.get("nonrigid_depth_wet_median_um", float("nan"))}
```

- [ ] **Step 6: Its history, results and frames without them** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -30,14 +30,10 @@
                 "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements", "deep_liquid_kg", "deep_mass_kg",
                 "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction", "cascade_passes", "cascade_mass_kg",
                 "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg",
-                "thin_invalid_patches", "thin_invalid_mass_kg"]
-SKIN_COLUMNS = ["skin_count", "skin_mass_kg", "skin_thickness_min_mm", "skin_thickness_median_mm", "skin_drawn_mass_kg",
-                "skin_handover_mass_kg", "skin_short_events", "interface_mismatch_J", "interface_repeats",
-                "liquid_depth_mean_um", "liquid_depth_max_um", "nonrigid_depth_skin_mean_um", "skin_wall_s",
-                "skin_thickness_max_mm", "skin_thick_count", "skin_shortfall_kg", "nonrigid_depth_wet_median_um",
-                "skin_frozen_mass_kg", "skin_returned_mass_kg", "skin_full_events", "skin_thick_feeding_count",
-                "skin_thick_feeding_max_mm", "skin_thick_feed_events"]   # sub-plan 18
-MELT_COLUMNS = MELT_COLUMNS + SKIN_COLUMNS
+                "thin_invalid_patches", "thin_invalid_mass_kg",
+                # the regime test's depth over the wet windward patches (fact 94); sub-plan 18's skin columns, which held
+                # 0 or nan in every run of the element model, went with its one-dimensional skins (sub-plan 22 Task 3)
+                "nonrigid_depth_wet_median_um"]
 # the layered skins on the corner columns (sub-plan 20): skinbody.SkinBody.skin_stats(), every macro step
 SKINBODY_COLUMNS = ["skin_cells", "skin_thickness_max_mm", "skin_above_min", "flowing_depth_mean_mm", "flowing_depth_max_mm",
                     "liquid_depth_mean_mm", "liquid_depth_max_mm", "surface_quality_min", "surface_fold_deg", "zone_ratio_min",
@@ -175,10 +171,6 @@
                 "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                 "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
                 "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
-                "skins_created": body.skins_created, "skin_drawn_mass_kg": body.skin_drawn_mass,
-                "skin_handover_mass_kg": body.skin_handover_mass, "interface_repeats": body.interface_repeats,
-                "skin_wall_s": body.skin_wall, "skin_frozen_mass_kg": body.skin_frozen_mass,
-                "skin_returned_mass_kg": body.skin_returned_mass, "skin_thick_feed_events": body.skin_thick_feed_events,
                 "melt_energy_balance_residual": body.energy_balance_residual()}
 
     def skin_results(self, history):
@@ -293,19 +285,6 @@
         # the depth the regime test read under each patch: down to the first rigid point, slurry included (amendment of
         # 2026-10-06); the spray step's own, carried like delta_m, nan where it was not evaluated
         poly.cell_data["nonrigid_depth"] = carry(body.last_nonrigid, np.nan)
-        # the skins (sub-plan 18): top temperature, resolved liquid and non-rigid depths, thickness; nan without a skin
-        sk = getattr(body, "skins", None)
-        if sk is not None:
-            on = body.skin_of_patch >= 0
-            r = body.skin_of_patch[on]
-            d_liq, _ = sk.depth_above(body.material.T_feed)
-            d_nr, _ = sk.depth_above(body.material.T_rigid)
-            for key, v in (("skin_T_top", sk.T[:, 0] if sk.n else None), ("skin_liquid_depth", d_liq),
-                           ("skin_nonrigid_depth", d_nr), ("skin_thickness", sk.thickness())):
-                out = np.full(surface.n_patches, np.nan)
-                if sk.n:
-                    out[on] = v[r]
-                poly.cell_data[key] = out
     poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
 
 
```

- [ ] **Step 7: The skins the default, the run names, and the spraying skin flight** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -4,6 +4,9 @@
         [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
         [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
         [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
+        [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off ...]
+                                     (the layered skins, the default with --thermal fem, their flow and spraying with
+                                      --skin-flow on; the element model with --surface-model elements)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
 Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
@@ -47,10 +50,12 @@
 # seeds, the droplet count within 2.7 % of the 0.00625 s step's (plan facts 92-94) -- while the default 0.5 s step makes
 # six times fewer droplets. Before the switch the default step serves.
 DEFAULT_DT_CONTINUUM = 0.0125
-# The surface models a run can choose (sub-plan 20): the element model, or the layered skins on the corners of a moving
-# mesh (skinbody.SkinBody). The default stays the element model until sub-plan 22 gives the skins their flow and spraying;
-# sub-plan 18's one-dimensional skins are reachable only from Python (body.MeltSettings(surface_model="skin1d")).
-SURFACE_MODELS = ("elements", "skin")
+# The surface models a run can choose (sub-plan 20): the layered skins on the corners of a moving mesh
+# (skinbody.SkinBody), the default of a thermal run since sub-plan 22 Task 3 gave them Girin's spraying (spec 10.1), or
+# the element model, kept unchanged for comparison (spec 9 and 11.1 level 3), whose runs end in `_surface-elements`.
+# Sub-plan 18's one-dimensional skins are retired (sub-plan 22 Task 3).
+SURFACE_MODELS = ("skin", "elements")
+DEFAULT_SURFACE_MODEL = "skin"
 # The fine band under a skin run's surface (Asha, 2026-10-14, option 1): the 3D model's elements are finer to this depth
 # while the surface keeps its triangulation, so that the 3D model follows the temperature just below the rigid front
 # (on 2 mm elements throughout the plates melted 2-6 % late), and every rebuild makes the band again
@@ -97,19 +102,37 @@
     (no film; its SESAM thresholds were measured at the default step) or a run that does not melt."""
     if args.dt_continuum is not None:
         return None if args.dt_continuum == "off" else args.dt_continuum
-    return DEFAULT_DT_CONTINUUM if (args.thermal == "fem" and args.melt == "on" and args.removal == "girin") else None
+    return DEFAULT_DT_CONTINUUM if (args.thermal == "fem" and args.melt == "on" and removal_of(args) == "girin") else None
+
+
+def surface_model_of(args):
+    """The run's surface model: --surface-model if given, else the layered skins for a thermal run (spec 10.1; sub-plan
+    22 Task 3); None without a thermal model, which has no surface to model."""
+    if args.surface_model is not None:
+        return args.surface_model
+    return DEFAULT_SURFACE_MODEL if args.thermal == "fem" else None
+
+
+def removal_of(args):
+    """--removal if given, else the run's own (sub-plan 22 Task 3): Girin's, for the element model and for the skins with
+    the flow on, whose liquid only spraying removes (spec 10.1: --skin-flow on needs --removal girin); instant removal for
+    the heat-only skins, whose liquid leaves as it forms (stage 2's stand-in for flow and spray, spec 11.2)."""
+    if args.removal is not None:
+        return args.removal
+    return "instant" if surface_model_of(args) == "skin" and args.skin_flow == "off" else "girin"
 
 
 def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
                    deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True, dt_continuum="default",
-                   surface_model="elements", skin_suffix=""):
+                   surface_model=DEFAULT_SURFACE_MODEL, skin_suffix=""):
     """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
     run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
     (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
     `_rigidsubstrate-off`, one whose continuum step is not its removal mode's default (amendment of 2026-10-07) in
-    `_dtcontinuum-off` or `_dtcontinuum-<s>`, a thermal run with the layered skins (sub-plan 20) in `_surface-skin`
-    followed by `skin_suffix` (the skins' non-default settings, `skin_name_suffix`), and any run with a seed other than
-    DEFAULT_SEED (amendment of 2026-10-05) in `_seed-<n>`; the defaults keep the old name."""
+    `_dtcontinuum-off` or `_dtcontinuum-<s>`; a thermal run with the element model in `_surface-elements`, the layered
+    skins being the default (spec 10.1; sub-plan 22 Task 3), and one with the skins in `skin_suffix`, their non-default
+    settings (`skin_name_suffix`, empty for the defaults); and any run with a seed other than DEFAULT_SEED (amendment of
+    2026-10-05) in `_seed-<n>`. Without a thermal model there is no surface model and no suffix for it."""
     dtc_default = DEFAULT_DT_CONTINUUM if melt == "girin" else None
     dtc = dtc_default if dt_continuum == "default" else dt_continuum
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
@@ -118,19 +141,21 @@
             + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else "")
             + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
             + ("" if not melt or dtc == dtc_default else "_dtcontinuum-off" if dtc is None else "_dtcontinuum-{:g}".format(dtc))
-            + ("_surface-skin" + skin_suffix if heating_name and surface_model == "skin" else "")
+            + (("_surface-elements" if surface_model == "elements" else skin_suffix) if heating_name else "")
             + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))
 
 
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
-                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=False, flow_cycles=None):
+                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=False, flow_cycles=None,
+                     slurry_spray=False):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
     preset's, `_growth<ratio>`, `_sub<ms>ms` and `_deform<mm>mm` for another growth, sub-step or deforming depth,
     `_noswell` for skins held at the solid's density, `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none), and,
     with the flow on (sub-plan 21 Task 6c), `_skinflow-on` followed by `_flowcycles<n>` where the flow sub-steps in each
-    heat sub-step are fixed at n rather than chosen each heat sub-step (Task 6d); empty for the defaults."""
+    heat sub-step are fixed at n rather than chosen each heat sub-step (Task 6d) and by `_slurryspray-on` where the
+    exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3); empty for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -154,6 +179,8 @@
         out += "_skinflow-on"
         if flow_cycles is not None:
             out += "_flowcycles{}".format(flow_cycles)
+    if slurry_spray:
+        out += "_slurryspray-on"
     return out
 
 
@@ -261,9 +288,11 @@
     th.add_argument("--stills", action="store_true", help="only the stills (start, peak heating, peak surface T, end; surface and section)")
     me = r.add_argument_group("melting and spraying (Step 3)")
     me.add_argument("--melt", choices=MELT_MODES, default="off", help="off: Step 2 behaviour (default); on: melting, melt film, spraying")
-    me.add_argument("--removal", choices=body.REMOVAL_NAMES, default="girin",
-                    help="girin: film + runoff + Girin spraying (default); instant: liquid removed as it forms -- the lumped-melting "
-                         "verification device, not physical")
+    me.add_argument("--removal", choices=body.REMOVAL_NAMES, default=None,
+                    help="girin: Girin's spraying -- the element model's film, runoff and spraying, or the skins' flow and "
+                         "spraying from their top cells (default, except for the heat-only skins); instant: liquid removed "
+                         "as it forms -- the lumped-melting verification device, not physical (the default of the "
+                         "heat-only skins, --skin-flow off)")
     me.add_argument("--runoff", choices=("on", "off"), default="on", help="film runoff transport along the surface (default on)")
     me.add_argument("--deep-runoff", choices=("on", "off"), default="on",
                     help="move the contiguous liquid below Girin's conjugate depth along the surface by the pressure-gradient and "
@@ -282,11 +311,11 @@
                          "with the film plus everything more than half liquid beneath it (slurry included), and where there "
                          "is no conjugate depth a film on slurry does not spray (default on). off: the regime test reads the "
                          "film plus the fully molten material, as in every run before the amendment of 2026-10-06")
-    me.add_argument("--surface-model", choices=SURFACE_MODELS, default="elements",
+    me.add_argument("--surface-model", choices=SURFACE_MODELS, default=None,
                     help="skin: a layered skin on every corner of the surface holds the outer metal in fine cells above a "
-                         "moving 3D mesh (layered-skin design, sub-plan 20: heat only, so with --melt on it needs --removal "
-                         "instant, liquid leaving as it forms). elements: the surface elements melt and feed the film, as "
-                         "in every run before it (default)")
+                         "moving 3D mesh (the layered-skin design; default with --thermal fem since sub-plan 22). elements: "
+                         "the surface elements melt and feed the film, as in every run before it, kept for comparison; "
+                         "its runs end in _surface-elements")
     me.add_argument("--skin-min-thickness", type=float, default=0.8,
                     help="the skins' minimum thickness [mm of solid]: where the non-rigid layer is thinner, a skin's bottom "
                          "is rigid metal (default %(default)s)")
@@ -309,9 +338,13 @@
                          "one density")
     me.add_argument("--skin-flow", choices=("on", "off"), default="off",
                     help="with --surface-model skin and --melt on: on, the skins' liquid and slurry flow over the surface "
-                         "(sub-plan 21) and stay -- the spraying that will take them is sub-plan 22's, so give --removal "
+                         "(sub-plan 21) and Girin's model sprays them from the top cells (sub-plan 22), with --removal "
                          "girin, whose continuum step (--dt-continuum) follows; off: the heat-only skins, liquid removed "
                          "as it forms with --removal instant (default)")
+    me.add_argument("--slurry-spray", choices=("on", "off"), default="off",
+                    help="with --skin-flow on: on, the exposed slurry sprays too, once no liquid covers it, by its own mode "
+                         "(the cascade of the layered-skin design, section 8); it sprays only within about a kelvin of the "
+                         "liquidus, on a viscosity law with no data there (default off)")
     me.add_argument("--skin-flow-cycles", type=int, default=None,
                     help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
                          "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
@@ -348,7 +381,8 @@
     me.add_argument("--size-feedback", choices=body.SIZE_FEEDBACK_NAMES, default=None,
                     help="current: Kn on the equivalent diameter of the remaining mass and the stagnation radius fitted to the windward cap "
                          "(default with --removal girin); initial: D0 and R0 throughout, SESAM's convention (default with --removal instant)")
-    me.add_argument("--particles", dest="particles", action="store_true", default=True, help="write the particle source table (default)")
+    me.add_argument("--particles", dest="particles", action="store_true", default=True,
+                    help="write the particle source table and the size distribution of a spraying run (default)")
     me.add_argument("--no-particles", dest="particles", action="store_false")
 
     c = sub.add_parser("compare", help="metrics and plots for an existing model history")
@@ -382,17 +416,20 @@
         mat.k_table = mat.k_table * args.k_scale                  # verification device (near-isothermal body), not physical
     solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=not args.consistent_mass)
     if skins:
-        flow = args.skin_flow == "on"                       # the liquid stays and flows; instant removal is instead of it
+        # with the flow on the liquid stays, flows and is sprayed from the top cells by Girin's model (sub-plan 22), the
+        # exposed slurry too with --slurry-spray on; instant removal is instead of both
+        flow = args.skin_flow == "on"
         flow_model = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm,
                                               kn_body_shock=args.kn_body_shock) if flow else None
         the_body = skinbody.SkinBody(the_mesh, mat, solver, mass, skin_settings(args), T0=args.temperature,
                                      emissivity=args.emissivity, T_ambient=args.t_ambient, deform_depth=args.deform_depth * 1e-3,
                                      remove=melting and not flow, demise_fraction=args.demise_fraction, flow=flow,
-                                     flow_cycles=args.skin_flow_cycles, surface_flow=flow_model)
+                                     flow_cycles=args.skin_flow_cycles, surface_flow=flow_model,
+                                     spray=spray_model_of(args, mat) if flow else None,
+                                     slurry_spray=args.slurry_spray == "on", particles=args.particles)
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
-        spray_model = spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical,
-                                       rt_spray=args.rt_spray == "on")
+        spray_model = spray_model_of(args, mat)
         size_feedback = args.size_feedback or ("current" if args.removal == "girin" else "initial")
         melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                           particles=args.particles, size_feedback=size_feedback,
@@ -420,13 +457,19 @@
                         "deform_depth_mm": args.deform_depth, "n_columns": the_body.columns.n,
                         "fine_band_mm": args.fine_band, "fine_size_mm": min(args.fine_size, h_surface * 1e3),
                         "swell": args.skin_swell == "on", "flow": args.skin_flow,
-                        "flow_cycles": "adaptive" if the_body.flow_cycles is None else the_body.flow_cycles}
+                        "flow_cycles": "adaptive" if the_body.flow_cycles is None else the_body.flow_cycles,
+                        "slurry_spray": args.slurry_spray}
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
     if melting and skins:
         info.update({"removal": args.removal, "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus,
                      "latent_heat_Jkg": mat.latent_heat})
+        if the_body.spray is not None:                  # Girin's spraying from the top cells (sub-plan 22)
+            info.update({"rt_spray": args.rt_spray, "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
+                         "k_r": args.kr, "k_t": args.kt, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
+                         "particles": args.particles,
+                         "liquid": {"rho": mat.liquid.rho, "mu": mat.liquid.mu, "sigma": mat.liquid.sigma}})
     elif melting:
         info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff,
                      "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
@@ -438,6 +481,11 @@
     return the_body, heating_model, info
 
 
+def spray_model_of(args, mat):
+    """The shared spray physics with the command line's constants, for either surface model."""
+    return spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical, rt_spray=args.rt_spray == "on")
+
+
 def skin_settings(args):
     """The layered skins' settings from the command line (sub-plan 20)."""
     cell = args.skin_cell * 1e-6 if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset]
@@ -446,11 +494,32 @@
                              fine_zone=args.skin_fine_zone * 1e-3)
 
 
+def run_name_of(args, diameter_m, velocity_ms, altitude_m, atmosphere_name, dt_continuum):
+    """The run's name (`model_run_name`) from its command line, its surface model and removal resolved (`surface_model_of`,
+    `removal_of`) and its continuum step (`continuum_step`); --name overrides it."""
+    return args.name or model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, args.bridging, args.wind,
+                                       args.heating if args.thermal == "fem" else None,
+                                       removal_of(args) if args.melt == "on" else None,
+                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on",
+                                       seed=args.seed, rigid_substrate=args.rigid_substrate == "on",
+                                       dt_continuum=dt_continuum, surface_model=surface_model_of(args),
+                                       skin_suffix=skin_name_suffix(args.skin_min_thickness, args.skin_cell, args.skin_growth,
+                                                                    args.skin_substep, args.skin_preset, args.deform_depth,
+                                                                    args.skin_swell == "on", args.fine_band, args.fine_size,
+                                                                    args.skin_fine_zone, args.skin_flow == "on",
+                                                                    args.skin_flow_cycles, args.slurry_spray == "on"))
+
+
 def cmd_run(args, parser):
     # first, before anything that could draw from it (amendment of 2026-10-05; plan facts 52 and 62): the only draws in a
     # run are pyamg's, at every hierarchy the skfem backend builds, and with the generator seeded here two runs of one
     # build agree bit for bit
     np.random.seed(args.seed)
+    if args.surface_model == "skin" and args.thermal != "fem":
+        parser.error("--surface-model skin needs --thermal fem")
+    # the surface model and the removal not given take the run's own (sub-plan 22 Task 3): the skins for a thermal run,
+    # and Girin's removal except for the heat-only skins
+    args.surface_model, args.removal = surface_model_of(args), removal_of(args)
     for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                          ("--cadence", args.cadence), ("--t-max", args.t_max)):
         if value <= 0.0:
@@ -496,8 +565,6 @@
     cell_um = args.skin_cell if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset] * 1e6
     if 2.0 * cell_um * 1e-3 > min(args.skin_min_thickness, args.skin_fine_zone) * (1.0 + 1e-12):
         parser.error("--skin-cell must be at most half of --skin-min-thickness and of --skin-fine-zone")
-    if args.surface_model == "skin" and args.thermal != "fem":
-        parser.error("--surface-model skin needs --thermal fem")
     if args.skin_flow == "on" and not (args.surface_model == "skin" and args.melt == "on"):
         parser.error("--skin-flow on needs --surface-model skin and --melt on: it moves the skins' molten metal")
     if args.skin_flow_cycles is not None and args.skin_flow != "on":
@@ -508,8 +575,20 @@
         parser.error("--skin-flow on keeps the liquid, which only spraying removes: it needs --removal girin, not instant "
                      "removal, which stands in for flow and spraying in the heat-only skins (spec 11.2)")
     if args.surface_model == "skin" and args.melt == "on" and args.skin_flow == "off" and args.removal != "instant":
-        parser.error("--surface-model skin with --melt on needs --removal instant, or --skin-flow on with Girin's removal "
-                     "(sub-plan 21: the skins' spraying is sub-plan 22's)")
+        parser.error("Girin's spraying from the skins needs --skin-flow on: the heat-only skins (--skin-flow off, the "
+                     "default) remove their liquid as it forms, --removal instant")
+    if args.slurry_spray == "on" and args.skin_flow != "on":
+        parser.error("--slurry-spray on needs --skin-flow on (with --surface-model skin and --melt on): the exposed slurry "
+                     "sprays by the skins' spraying from their top cells")
+    if args.surface_model == "skin":
+        for flag, given in (("--runoff off", args.runoff == "off"), ("--deep-runoff off", args.deep_runoff == "off"),
+                            ("--molten-cascade off", args.molten_cascade == "off"),
+                            ("--rigid-substrate off", args.rigid_substrate == "off"),
+                            ("--size-feedback", args.size_feedback is not None)):
+            if given:
+                parser.error("{} is the element model's setting: give --surface-model elements (the layered skins, the "
+                             "default since sub-plan 22, have no film, runoff, deep account or element deaths, and their "
+                             "regime test always reads their own depths)".format(flag))
     if args.surface_model == "skin" and (args.animate or args.stills):
         parser.error("--animate/--stills read the element model's frames; with --surface-model skin use --frames-every")
     if args.surface_model == "skin" and args.prism_layers is not None and args.prism_layers > 0:
@@ -524,16 +603,7 @@
     mass = body.sphere_mass(settings.diameter, args.material_density)
     atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
     reference = sesam_io.load_reference(args.reference) if args.reference else None
-    name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
-                                       args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None,
-                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on", seed=args.seed,
-                                       rigid_substrate=args.rigid_substrate == "on", dt_continuum=dt_continuum,
-                                       surface_model=args.surface_model,
-                                       skin_suffix=skin_name_suffix(args.skin_min_thickness, args.skin_cell, args.skin_growth,
-                                                                    args.skin_substep, args.skin_preset, args.deform_depth,
-                                                                    args.skin_swell == "on", args.fine_band, args.fine_size,
-                                                                    args.skin_fine_zone, args.skin_flow == "on",
-                                                                    args.skin_flow_cycles))
+    name = run_name_of(args, settings.diameter, initial.velocity, initial.altitude, atm_name, dt_continuum)
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
@@ -575,10 +645,16 @@
         "files": {"csv": os.path.abspath(csv_path)},
     }
     melting = args.thermal == "fem" and args.melt == "on" and args.surface_model != "skin"
+    # the skins' spraying flight (sub-plan 22 Task 3): the element model's droplet files, from the same writer; the
+    # per-step summary waits for Task 4, which gives the skin history its droplet columns, and is recorded as not written
+    spraying_skins = args.thermal == "fem" and args.melt == "on" and args.surface_model == "skin" and args.skin_flow == "on"
     if melting:
         doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()} if args.particles
                             else {})
         doc["files"]["melt_plots"] = [os.path.abspath(p) for p in compare.plot_melt(history, run_dir, name, reference if compare.has_melt(history, reference) else None)]
+    elif spraying_skins and args.particles:
+        doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()})
+        doc["files"]["particles_summary"] = None
     if reference is not None:
         plots = compare.plot_all(history, reference, run_dir, name)
         doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
@@ -617,6 +693,11 @@
                   "wall time {:.0f} s in the skins, {:.0f} s in the 3D solve".format(
                       res["skin_columns"], res["removed_mass_kg"], res["drawn_mass_kg"], res["rebuilds"],
                       res["interface_repeats"], res["jump_faults"], res["wall_skins_s"], res["wall_3d_s"]))
+        if spraying_skins:
+            print("  spraying: released {:.4g} kg as {:.3g} droplets in {} source rows{}; no particles_summary.csv until "
+                  "sub-plan 22 Task 4 gives the skin history its droplet columns".format(
+                      the_body.released_mass, the_body.droplet_number, sum(len(b) for b in the_body.source_rows),
+                      " (particles.npz, size_distribution.csv)" if args.particles else ""))
         if melting:
             print("  melt: onset {} km, spraying onset {} km, demise {} km at t = {} s; sprayed {:.4f} kg of {:.4f}, film left {:.4f} kg, "
                   "{:.3g} droplets (median r {} um), {} source rows, melt balance {:.1e}".format(
```

- [ ] **Step 8: The element model's scripts name it** in `analysis/reentry_model_thermal_verification.py`,
  `analysis/melt_verification.py` and `analysis/melt_sensitivity.py`:

```diff
--- a/analysis/reentry_model_thermal_verification.py
+++ b/analysis/reentry_model_thermal_verification.py
@@ -54,9 +54,7 @@
     ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
     argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                              "--heating", mode, "--reference", ref, "--outdir", outdir, "--name", label_of(key, mode, surface, t_max),
-                             "--quiet"]
-    if surface == "skin":
-        argv += ["--surface-model", "skin"]
+                             "--quiet", "--surface-model", surface]       # named: the skins are the model's default
     if t_max:
         argv += ["--t-max", "{:g}".format(t_max)]
     if animate:
```

```diff
--- a/analysis/melt_verification.py
+++ b/analysis/melt_verification.py
@@ -45,7 +45,8 @@
     ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
     label = "{}__{}".format(key, mode)
     argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
-                             "--melt", "on"] + MODES[mode] + ["--reference", ref, "--outdir", outdir, "--name", label, "--quiet"]
+                             "--melt", "on", "--surface-model", "elements"] + MODES[mode] + ["--reference", ref, "--outdir",
+                                                                                            outdir, "--name", label, "--quiet"]
     if animate:
         argv.append("--animate")
     rc = cli.main(argv)
```

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -47,7 +47,8 @@
 
 def argv_for(key, variant, outdir):
     return ["run"] + CASES[key] + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
-                                   "--melt", "on", "--heating", "physics", "--material", "AA7075_scheil", "--no-particles",
+                                   "--melt", "on", "--surface-model", "elements", "--heating", "physics", "--material",
+                                   "AA7075_scheil", "--no-particles",
                                    "--outdir", outdir, "--name", "{}__{}".format(key, variant), "--quiet"] + VARIANTS[variant]
 
 
```

- [ ] **Step 9: Run them to see them pass.** Expected: `tests/test_reentry_model_cli.py` 81 passed (the spraying run about 20–30 s),
  `tests/test_reentry_model_melting.py` 27, `tests/test_reentry_model_thermal_verification.py` 6; the unit tier 519
  passed, 1 skipped (Task 2's 538 less the 35 tests of stage 1's skins, plus 16); the FEniCSx tier 17 passed and 1 failed
  (`test_the_thick_film_flux_matches_the_skfem_backend`, below).

- [ ] **Step 10: The measurements** (in `measure/task3/`, on mains power): the element model's short case,
  `../task1/cmp/run_case.sh <code> runs line3_short --surface-model elements`, against Task 2's `line2_short` copied into
  `runs/`: `"$PY" -B ../task1/cmp/compare_runs.py runs line2_short line3_short` and `"$PY" -B element_files.py runs
  line2_short line3_short`; the smoke flights, `./run_smoke.sh` (two runs, about a minute and a half each), and
  `"$PY" -B smoke_check.py smoke/<run>.json smoke/<log>` for each; the macro step's cost, `./run_cost.sh` (about three
  minutes, alone) and `smoke_check.py` on its two runs; and the FEniCSx tier.

- [ ] **Step 11: Commit.**

```bash
git add -A reentry_model/body.py reentry_model/coupled.py reentry_model/cli.py reentry_model/skin1d.py analysis/reentry_model_thermal_verification.py analysis/melt_verification.py analysis/melt_sensitivity.py tests/test_reentry_model_cli.py tests/test_reentry_model_melting.py tests/test_reentry_model_thermal_verification.py tests/test_reentry_model_fenicsx.py tests/test_reentry_model_skin1d.py tests/test_reentry_model_melting_skin.py
git commit -m "Skins: stage 1's one-dimensional skins retired and the layered skins the default surface model (spec 9 and 10.1) -- skin1d.py, its tests and the element model's branches, history columns and results for it removed, the element model unchanged to the last bit in everything it still writes; a thermal run's surface model the skins unless --surface-model elements, whose runs end in _surface-elements, a skin run's name its non-default settings only; --removal not given Girin's but for the heat-only skins; contradictory settings, the element model's own given to the skins among them, exit 2; with --skin-flow on the command line gives the skins the shared spray model, --slurry-spray (_slurryspray-on) and --particles and writes their particles.npz and size distribution, the summary waiting for Task 4 (sub-plan 22 Task 3)"
```

---

### Task 3b: The flow on by default, and the continuum step from the skins' first liquid

**Files:**
- Modify: `reentry_model/cli.py` (`skin_flow_of`; `removal_of` reading it; `--skin-flow` defaulting to the run's own;
  `skin_name_suffix(..., flow=True)`, `_skinflow-off`; `run_name_of`; `cmd_run`'s resolution and its
  `CoupledSettings(..., switch_on_melt=)`; the `--dt-continuum` help), `reentry_model/coupled.py`
  (`CoupledSettings.switch_on_melt`; `CoupledRun.switch_trigger`; `_switch`; the results' `dt_switch_trigger`),
  `reentry_model/skinbody.py` (`SkinBody.liquid_present`; `FLOW_DEFAULT`'s comment)
- Test: `tests/test_reentry_model_cli.py` (one test added, nine changed), `tests/test_reentry_model_coupled.py` (one
  added), `tests/test_reentry_model_skinbody.py` (one added)

**Interfaces:**
- Consumes: Task 3's `surface_model_of`, `removal_of`, `run_name_of`, `skin_name_suffix`, `continuum_step`; the
  coupled loop's continuum switch (the amendment of 2026-10-07); `Skins.T`, `.valid`.
- Produces: `cli.skin_flow_of(args)` (on for the skins' melting runs unless `--removal instant` is given, off where the
  flow does not apply); `--skin-flow` not given resolved, its value in the run's settings; `skin_name_suffix(...,
  flow=True)`, `flow` False giving `_skinflow-off` and None (no melting) nothing; `coupled.CoupledSettings(...,
  switch_on_melt=False)`; `CoupledRun.switch_trigger` (`knudsen`, `melting`, `both` or None); the results'
  `dt_switch_trigger` where `switch_on_melt`; `SkinBody.liquid_present() -> bool`. The element model's settings,
  history and results are unchanged; no skin history column or frame field changes.

Asha decided Task 3's two questions on 2026-10-29 (spec 10.1 and 16 as amended that day; the amendment above): the flow
is on by default for the skins' melting runs (option B), and on the skins the continuum step also latches at the first
macro step at whose end they hold liquid, or at the Knudsen switch if that comes first; the element model's switch is
unchanged.

**The flow on by default.** `--skin-flow` not given is the run's own (`skin_flow_of`, resolved in `cmd_run` beside the
surface model and the removal, and kept in the run's settings): on for a melting run with the skins, so it flows,
sprays by Girin's model and writes the droplet files with no flags; off where `--removal instant` is given, the
heat-only skins' removal (ruling 1); off for a run that does not melt or does not use the skins, which the flow does not
concern, as before. `removal_of` follows it as in Task 3, so a melting skin run with no flags takes Girin's removal and
its continuum step. A flow run's name loses `_skinflow-on`; the heat-only skins of a melting run are named
`_skinflow-off` (spec 10.1), and a run that does not melt carries nothing for the flow. The contradictions of Task 3
stand, each now reached only by explicit settings that disagree (`--skin-flow off` with `--removal girin`, `--skin-flow
on` with `--removal instant`, `--slurry-spray on` without the flow). `skinbody.FLOW_DEFAULT` stays off: it is the
class's default, which the heat-only skins' tests build on (ruling 2).

**The continuum step from the skins' first liquid.** `CoupledRun._switch` latches the continuum step at the first macro
step whose body Knudsen number is below `kn_switch`, as before, or, with `CoupledSettings.switch_on_melt`, at the first
at whose end the body holds liquid (`body.liquid_present()`), whichever comes first. The command line sets it for the
skins (`switch_on_melt=args.surface_model == "skin"`), so it acts where a skin run has a continuum step: Girin's removal,
that is the flow. Onset is liquid above the liquidus (ruling 3): `SkinBody.liquid_present` is true once any column's
top cell is above T_liquidus, which is where the liquid depth -- the depth option B reads for the Rayleigh-Taylor rules
and the hold -- is not zero. The results record the trigger beside the `dt_switch_*` fields, which hold the latching
step's flight time, Knudsen number and altitude: `dt_switch_trigger` is `melting`, `knudsen`, or `both` where the two
come at the same step (ruling 4), and only where `switch_on_melt` (ruling 5). Without it the body is never asked, so
the element model's switch, its outputs and its results are unchanged. No opt-out flag is added (ruling 6).

**Measured** (2026-10-29, `prototype/work-2026-10-28-skin-spraying/measure/task3b/`; facts 180–182):
- *The element model unchanged* (`runs/line3b_short` against Task 3's `line3_short`; `compare_runs.py`,
  `element_files.py`). The CLI test's short melting case with `--surface-model elements`: the history file byte for
  byte (93 columns, 76 rows), all 26 fields of the 25 266 source rows, the summary and the size distribution byte for
  byte, the settings and all 61 results but the run time. Its switch is the Knudsen number's at 9.5 s, after melting
  at 1.5 s, and its results carry no trigger.
- *The melting trigger* (`evidence/`, `onset_check.py`; `run_evidence.sh`). The CLI tests' 71 km case as Task 3's spike
  ran it -- the 100 mm sphere on the coarse mesh (2 367 columns, 20 um cells) from 880 K, physics heating, a 0.25 s first
  step -- now with no flow flag, to 0.5 s: the body Knudsen number is 0.0112-0.0113 throughout, so the Knudsen switch
  never fires; the nose holds 0.12 mm of liquid at 0.25 s and the step latches there, `melting`, at 70.97 km, the next
  20 steps 12.5 ms. The liquid grows to 0.20 mm over a flowing layer of 0.73 rising to 1.29 mm, and nothing is
  released: off Girin's closure a liquid on slurry is held from the shear modes (option B). 49 s, 0.64 GB. The 65 km
  smoke case (Task 3's, Knudsen number 0.005) with no flow flag: `knudsen` at the first step's end (25 ms), the first
  liquid two steps later (50 ms); it sprays as Task 3's smoke flight did, 1.95572 g as 76 300 droplets in 1 196 rows.
- *Onset read elsewhere* (the same runs; the element model's flights). At 71 km a top at T_rigid or at the solidus
  comes at the same macro step as the liquid (the first); at 65 km the Knudsen switch comes first whichever is read
  (T_rigid and the solidus at the first step, the liquid at the third). On the element model's 100 mm flight (the
  thin-gate flight of record, 0.5 s history) the hottest surface reaches the solidus at 13.5 s, T_rigid at 23.0 s and
  the liquidus at 24.5 s: read at T_rigid the trigger would come three 0.5 s steps earlier, at the solidus 22.
- *The cost, estimated* (not measured on a skin flight). The skins' first liquid is taken as the element model's
  hottest surface reaching the liquidus: 24.5 s on the 100 mm flight (74.1 km, Knudsen number 0.018), whose Knudsen
  switch is at 49.5 s (fact 69), so 25.0 s of flight at 12.5 ms, 2 000 macro steps where 0.5 s steps made 50: 1 950
  more; on the 50 mm flight (the element model's of 2026-10-07 from 115 km, with the linear melting range: no Scheil
  50 mm flight exists) 174.5 s (78.3 km, Knudsen number 0.068), and the Knudsen number never falls below 0.01 before the
  demise at 199.5 s, so that flight had no continuum step at all and now has 25.0 s of it, again 1 950 more. The skins'
  tops ran 2-3 K above the element model's surface in level 3 (fact 149), so their onset may come a little earlier. On
  Task 3's coarse case the continuum step cost about half as much a second of flight as 62.5 ms steps (fact 179); what
  each extra step adds regardless is its fixed work -- the surface flow (0.06 s), Girin's plan (0.09-0.51 s) and a 3D
  solve each pass (1.17 s on the 100 mm sphere's fine band, fact 140) -- about 3 s, so about 1.6 h a flight, beside a
  row of the history and the source rows of every spraying column.
- *The tests* fail first on Task 3's code (9 failed, 27 passed, `log_fail_first.txt`) and pass (`log_slow_tests.txt`).

**Rulings** (2026-10-29; each what, why, cost if wrong):
1. *An explicit `--removal instant` turns the skins' flow off*, as an explicit `--skin-flow off` makes the removal
   instant: each setting not given takes the value the other one given needs, and only two given that disagree exit 2.
   Why: Asha's "the way `removal_of` resolves the removal"; a command line of Task 3's heat-only skins
   (`--removal instant`) keeps its meaning instead of exiting 2. Cost if wrong: such a run is the heat-only skins and its
   name says so (`_melt-instant_skinflow-off`); the alternative refuses it and asks for `--skin-flow off`.
2. *`skinbody.FLOW_DEFAULT` stays off*, its comment saying why. Why: it is the Python class's default, on which the
   heat-only skins' tests (most of the skin body's) build; what a run takes is the command line's. Cost: a caller
   building `SkinBody` directly gets the heat-only body unless it passes `flow=True`.
3. *Onset is a top cell above the liquidus*, read at the end of each macro step until the step latches. Why: Asha's
   choice, the liquid depth of option B; a positive liquid depth is exactly a top cell above the liquidus. Cost: a
   comparison per column a macro step. Read at T_rigid or at the solidus it would come at the same macro step on both
   test cases, but 1.5 s (three 0.5 s steps) and 11 s (22 steps) earlier on the element model's 100 mm flight (fact
   181): for Asha.
4. *`both` where the Knudsen number and the liquid come at the same step.* Why: it says both; the step is the same.
   Cost: none.
5. *The trigger is recorded only where the melting trigger is armed (the skins).* Why: the element model's results stay
   the same key for key; its trigger is always the Knudsen number. Cost: an element run's results have no
   `dt_switch_trigger`; a reader takes its absence as `knudsen`.
6. *No opt-out flag.* Why: `--dt-continuum off` already turns the continuum step off, and a skin run that latches on
   the Knudsen number alone is needed only to compare with the element model's rule, which Task 5's study of the macro
   step may want; it can add one then, named by spec 10.1 when not the default. Cost: until then that comparison needs
   code.

- [ ] **Step 1: The flow's default, its names and the melting trigger from the command line** -- in
  `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -283,6 +283,8 @@
     assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
     assert s["molten_cascade"] == "on" and "cascade_passes" in rows[0] and "cascade_mass_kg" in rows[0]      # the molten cascade
     assert s["dt_continuum_s"] == 0.1 and r["dt_continuum_s"] == 0.1 and r["dt_switch_kn"] < s["kn_body_shock"]
+    # the element model's switch is the Knudsen number's alone (sub-plan 22 Task 3b): it melts long before it latches
+    assert r["melt_onset_time_s"] < r["dt_switch_time_s"] - 1.0 and "dt_switch_trigger" not in r
     assert s["rigid_substrate"] == "on" and all(k in rows[0] for k in ("nonrigid_depth_mean_mm", "slurry_thick_fraction",
                                                                        "rigid_thin_fraction", "slurry_held_mass_kg"))
     assert r["cascade_mass_kg"] == pytest.approx(float(rows[-1]["cascade_mass_kg"])) and r["cascade_mass_kg"] >= 0.0
@@ -378,7 +380,7 @@
     assert cli.main(argv) == 0
     (path,) = glob.glob(str(tmp_path / "*.json"))
     doc = json.load(open(path))
-    assert doc["run_name"].endswith("_fem-sesam_melt-instant_skin-dev_fine4x2mm")
+    assert doc["run_name"].endswith("_fem-sesam_melt-instant_skin-dev_fine4x2mm_skinflow-off")
     s, r = doc["settings"], doc["results"]
     assert s["surface_model"] == "skin" and s["prism_layers"] == 0
     assert s["skin"] == {"min_thickness_mm": 0.8, "fine_zone_mm": 0.4, "cell_um": 20.0, "growth": 1.1, "substep_ms": 10.0, "preset": "dev",
@@ -404,7 +406,7 @@
     MELT + ["--removal", "instant", "--surface-model", "skin", "--deform-depth", "0"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--fine-band", "-1"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--fine-size", "0"],
-    MELT + ["--removal", "girin", "--surface-model", "skin"],               # Girin's spraying needs the flow
+    MELT + ["--removal", "girin", "--surface-model", "skin", "--skin-flow", "off"],   # Girin's spraying needs the flow
     MELT + ["--removal", "instant", "--surface-model", "skin", "--animate"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--prism-layers", "4"],     # the rebuild cannot make prisms
     ["run", "--surface-model", "skin"],                                     # no thermal model
@@ -433,25 +435,29 @@
 # Sub-plan 21 Task 6c: the flow in the skins, its setting and its sub-cycles
 
 def test_the_flow_settings_reach_the_run_name_only_when_not_default():
-    """Task 6d: the flow sub-steps are chosen each skin sub-step by default, so any fixed count is named."""
+    """Task 6d: the flow sub-steps are chosen each skin sub-step by default, so any fixed count is named. Sub-plan 22 Task
+    3b: the flow is on by default with melting (Asha, 2026-10-29), so the heat-only skins are named _skinflow-off, and a
+    run that does not melt (flow None) carries nothing for it."""
     from reentry_model import skin
     suffix = lambda **k: cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, **k)
-    assert suffix() == suffix(flow=False) == ""
-    assert suffix(flow=True) == suffix(flow=True, flow_cycles=None) == "_skinflow-on"
-    assert suffix(flow=True, flow_cycles=skin.FLOW_CYCLES) == "_skinflow-on_flowcycles{}".format(skin.FLOW_CYCLES)
-    assert suffix(flow=True, flow_cycles=2) == "_skinflow-on_flowcycles2"
-    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "dev", 15.0, swell=False, flow=True) == "_skin-dev_noswell_skinflow-on"
+    assert suffix() == suffix(flow=True) == suffix(flow=None) == ""
+    assert suffix(flow=False) == "_skinflow-off"
+    assert suffix(flow_cycles=skin.FLOW_CYCLES) == "_flowcycles{}".format(skin.FLOW_CYCLES)
+    assert suffix(flow_cycles=2) == "_flowcycles2"
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "dev", 15.0, swell=False, flow=False) == "_skin-dev_noswell_skinflow-off"
     args = cli.build_parser().parse_args(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.5"])
-    assert args.skin_flow == "off" and args.skin_flow_cycles is None
+    assert args.skin_flow is None and args.skin_flow_cycles is None
 
 
 def test_a_skin_run_with_flow_sprays_from_its_top_cells_and_writes_the_droplet_files(tmp_path):
     """The spraying skin flight from the command line (sub-plan 22 Task 3). 62.5 ms from 65 km, where the body Knudsen
     number is 0.005 and Girin's closure holds, of a body at 893 K (2 K below T_rigid) under physics heating, with the
-    layered skins (dev preset) and the flow on at two flow sub-steps a skin sub-step: Girin's removal is the default with
-    the flow, and so is its continuum step, so a 25 ms first step is followed by three of 12.5 ms. The tops melt by the
-    last step's start and its plan sprays them. The run name ends in _melt-girin and the skins' own suffixes, with no
-    surface suffix (the skins are the default); the settings record the flow and the spray model's constants; mass,
+    layered skins (dev preset) and the flow at two flow sub-steps a skin sub-step: the flow is the default of a melting
+    skin run (Task 3b), Girin's removal the default with the flow, and so is its continuum step, so a 25 ms first step is
+    followed by three of 12.5 ms, latched by the Knudsen number at the first step's end, before the skins hold liquid.
+    The tops melt by the last step's start and its plan sprays them. The run name ends in _melt-girin and the skins' own
+    suffixes, with no surface or flow suffix (both defaults); the settings record the flow and the spray model's
+    constants, the results the switch's trigger; mass,
     energy and momentum close; the droplets are written by the element model's writer -- particles.npz, whose rows carry
     the mass released, and size_distribution.csv, their histogram -- and the per-step summary, which waits for Task 4's
     history columns, is recorded as not written; the history's flow columns and the frames' top speed are written."""
@@ -461,18 +467,20 @@
     from reentry_model import spray
     argv = MELT + ["--heating", "physics", "--altitude", "65", "--velocity", "7.24", "--temperature", "893", "--dt", "0.025",
                    "--t-max", "0.0625", "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "1",
-                   "--fine-size", "2", "--skin-flow", "on", "--skin-flow-cycles", "2", "--outdir", str(tmp_path)]
+                   "--fine-size", "2", "--skin-flow-cycles", "2", "--outdir", str(tmp_path)]
     assert cli.main(argv) == 0
     (path,) = glob.glob(str(tmp_path / "*.json"))
     doc = json.load(open(path))
-    assert doc["run_name"].endswith("_fem-physics_melt-girin_skin-dev_fine4x2mm_skinflow-on_flowcycles2")
+    assert doc["run_name"].endswith("_fem-physics_melt-girin_skin-dev_fine4x2mm_flowcycles2")
     s, r, f = doc["settings"], doc["results"], doc["files"]
     assert s["surface_model"] == "skin" and s["removal"] == "girin" and s["skin"]["flow"] == "on"
     assert s["skin"]["flow_cycles"] == 2 and s["skin"]["slurry_spray"] == "off" and s["particles"] is True
     assert s["we_critical"] == 3.08 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["rt_spray"] == "on"
     assert s["dt_continuum_s"] == 0.0125 and r["dt_switch_time_s"] == pytest.approx(0.025)
+    assert r["dt_switch_trigger"] == "knudsen" and r["dt_switch_kn"] < 0.01
     rows = list(csv.DictReader(open(path[:-5] + ".csv")))
     assert len(rows) == 5 and all(k in rows[0] for k in coupled.SKINBODY_COLUMNS)
+    assert float(rows[1]["liquid_depth_max_mm"]) == 0.0                     # no liquid yet where the step latched
     released = r["released_mass_kg"]
     assert released > 0.0 and float(rows[-2]["released_mass_kg"]) == 0.0 and r["removed_mass_kg"] == 0.0
     # the books: mass, energy and momentum
@@ -512,19 +520,24 @@
 
 def test_the_skins_are_the_default_and_girin_s_removal_follows_their_flow():
     """Spec 10.1: --surface-model defaults to the skins for a thermal run (a run without a thermal model has no surface
-    model). --removal not given is the run's own: Girin's for the element model and for the skins with the flow on (spec
-    10.1, --skin-flow on needs --removal girin), instant removal for the heat-only skins (the flow off, its default), whose
-    liquid leaves as it forms; the continuum step (0.0125 s) follows Girin's removal, so a spraying skin run takes it."""
+    model). --skin-flow not given is on for the skins' melting runs (sub-plan 22 Task 3b; Asha, 2026-10-29) unless instant
+    removal is given, and off where the flow does not apply. --removal not given is the run's own: Girin's for the
+    element model and for the skins with the flow on (spec 10.1, --skin-flow on needs --removal girin), instant removal
+    for the heat-only skins, whose liquid leaves as it forms; the continuum step (0.0125 s) follows Girin's removal, so a
+    skin melting run with no flags flows, sprays and takes it."""
     p = cli.build_parser()
     melt = THERMAL + ["--melt", "on"]
-    resolved = lambda argv: (lambda a: (cli.surface_model_of(a), cli.removal_of(a), cli.continuum_step(a)))(p.parse_args(argv))
+    resolved = lambda argv: (lambda a: (cli.surface_model_of(a), cli.skin_flow_of(a), cli.removal_of(a),
+                                        cli.continuum_step(a)))(p.parse_args(argv))
     assert cli.DEFAULT_SURFACE_MODEL == "skin" and p.parse_args(THERMAL).surface_model is None
-    assert resolved(BASE) == (None, "girin", None)                               # no thermal model, no surface model
-    assert resolved(THERMAL) == ("skin", "instant", None)                        # nothing melts: the removal is unused
-    assert resolved(melt) == ("skin", "instant", None)                           # the heat-only skins
-    assert resolved(melt + ["--skin-flow", "on"]) == ("skin", "girin", 0.0125)   # the spraying skins
-    assert resolved(melt + ["--surface-model", "elements"]) == ("elements", "girin", 0.0125)
-    assert resolved(melt + ["--surface-model", "elements", "--removal", "instant"]) == ("elements", "instant", None)
+    assert resolved(BASE) == (None, "off", "girin", None)                        # no thermal model, no surface model
+    assert resolved(THERMAL) == ("skin", "off", "instant", None)                 # nothing melts: flow and removal unused
+    assert resolved(melt) == ("skin", "on", "girin", 0.0125)                     # the spraying skins, the default
+    assert resolved(melt + ["--skin-flow", "on"]) == resolved(melt)
+    assert resolved(melt + ["--skin-flow", "off"]) == ("skin", "off", "instant", None)    # the heat-only skins
+    assert resolved(melt + ["--removal", "instant"]) == ("skin", "off", "instant", None)  # ... chosen by their removal
+    assert resolved(melt + ["--surface-model", "elements"]) == ("elements", "off", "girin", 0.0125)
+    assert resolved(melt + ["--surface-model", "elements", "--removal", "instant"]) == ("elements", "off", "instant", None)
 
 
 def test_the_run_names_mark_the_element_model_and_the_skins_non_default_settings():
@@ -539,13 +552,15 @@
     assert name(BASE) == base
     assert name(THERMAL) == base + "_fem-sesam"
     assert name(FEM) == base + "_fem-sesam_surface-elements"
-    assert name(melt) == base + "_fem-sesam_melt-instant"
-    assert name(melt + ["--skin-preset", "dev"]) == base + "_fem-sesam_melt-instant_skin-dev"
-    assert name(melt + ["--skin-flow", "on"]) == base + "_fem-sesam_melt-girin_skinflow-on"
-    assert name(melt + ["--skin-flow", "on", "--slurry-spray", "on"]) == base + "_fem-sesam_melt-girin_skinflow-on_slurryspray-on"
-    assert (name(melt + ["--skin-flow", "on", "--skin-flow-cycles", "2", "--slurry-spray", "on", "--seed", "7"])
-            == base + "_fem-sesam_melt-girin_skinflow-on_flowcycles2_slurryspray-on_seed-7")
-    assert name(melt + ["--skin-flow", "on", "--dt-continuum", "off"]) == base + "_fem-sesam_melt-girin_dtcontinuum-off_skinflow-on"
+    assert name(melt) == name(melt + ["--skin-flow", "on"]) == base + "_fem-sesam_melt-girin"    # the flow, Task 3b's default
+    assert name(melt + ["--skin-preset", "dev"]) == base + "_fem-sesam_melt-girin_skin-dev"
+    assert name(melt + ["--slurry-spray", "on"]) == base + "_fem-sesam_melt-girin_slurryspray-on"
+    assert (name(melt + ["--skin-flow-cycles", "2", "--slurry-spray", "on", "--seed", "7"])
+            == base + "_fem-sesam_melt-girin_flowcycles2_slurryspray-on_seed-7")
+    assert name(melt + ["--dt-continuum", "off"]) == base + "_fem-sesam_melt-girin_dtcontinuum-off"
+    assert (name(melt + ["--skin-flow", "off"]) == name(melt + ["--removal", "instant"])
+            == base + "_fem-sesam_melt-instant_skinflow-off")                     # the heat-only skins
+    assert name(THERMAL + ["--skin-flow", "off"]) == name(THERMAL)               # no melting: the flow does not apply
     elements = melt + ["--surface-model", "elements"]
     assert name(elements) == base + "_fem-sesam_melt-girin_surface-elements"
     assert name(elements + ["--removal", "instant"]) == base + "_fem-sesam_melt-instant_surface-elements"
@@ -553,13 +568,13 @@
             == base + "_fem-sesam_melt-girin_deeprunoff-off_rigidsubstrate-off_dtcontinuum-off_surface-elements_seed-7")
     assert name(elements + ["--skin-preset", "dev"]) == name(elements)          # the skins' settings name skin runs only
     assert name(melt + ["--name", "mine"]) == "mine"
-    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, slurry_spray=True) == \
-        "_skinflow-on_slurryspray-on"
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, slurry_spray=True) == "_slurryspray-on"
 
 
 @pytest.mark.parametrize("argv", [
-    THERMAL + ["--melt", "on", "--removal", "girin"],                                  # Girin's spraying needs the flow
-    THERMAL + ["--melt", "on", "--slurry-spray", "on"],                                # the slurry sprays by the flow's spraying
+    THERMAL + ["--melt", "on", "--skin-flow", "off", "--removal", "girin"],            # Girin's spraying needs the flow
+    THERMAL + ["--melt", "on", "--removal", "instant", "--slurry-spray", "on"],        # the slurry sprays by the flow's spraying
+    THERMAL + ["--melt", "on", "--skin-flow", "off", "--slurry-spray", "on"],
     THERMAL + ["--melt", "on", "--skin-flow", "on", "--slurry-spray", "maybe"],
     FEM + ["--melt", "on", "--slurry-spray", "on"],                                    # the element model has no cascade
     THERMAL + ["--slurry-spray", "on"],                                                # nothing melts
@@ -590,8 +605,9 @@
 
 def test_the_command_line_builds_the_skins_with_the_shared_spray_model(tmp_path, monkeypatch):
     """--skin-flow on gives the skin body the shared spray physics with the command line's constants (the element model's
-    --kr, --kt, --we-critical, --rt-spray), --slurry-spray and --particles; the heat-only skins get none of them and
-    remove their liquid as it forms. The body is not built: its arguments are."""
+    --kr, --kt, --we-critical, --rt-spray), --slurry-spray and --particles, and so does a melting run with no flags (the
+    flow on by default, Task 3b); the heat-only skins (--skin-flow off) get none of them and remove their liquid as it
+    forms. The body is not built: its arguments are."""
     from reentry_model import skinbody
 
     def built(*args, **kwargs):
@@ -607,4 +623,32 @@
     with pytest.raises(_Built) as got:
         cli.main(melt)
     k = got.value.args[0]
+    assert k["flow"] and not k["remove"] and k["spray"] is not None and k["spray"].k_r == 0.17 and k["particles"]
+    with pytest.raises(_Built) as got:
+        cli.main(melt + ["--skin-flow", "off"])
+    k = got.value.args[0]
     assert not k["flow"] and k["remove"] and k["spray"] is None and k["surface_flow"] is None and not k["slurry_spray"]
+
+
+def test_the_skins_continuum_step_latches_at_their_first_liquid_before_the_knudsen_switch(tmp_path):
+    """Sub-plan 22 Task 3b (Asha, 2026-10-29): on the skins the continuum step also latches at the first macro step at
+    whose end they hold liquid above the liquidus, where that comes before the Knudsen switch. 112.5 ms from 71 km, where
+    the body Knudsen number is 0.011 and the Knudsen switch does not fire, of a body at 893 K under physics heating, with
+    the skins' default flow (dev preset, one flow sub-step a skin sub-step) and a 0.1 s first step: the nose holds liquid
+    at its end, so the step latches there to 12.5 ms, and the results say it was the melting, at the step's flight time,
+    Knudsen number and altitude."""
+    pytest.importorskip("cantera")
+    import glob
+    argv = MELT + ["--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "893", "--dt", "0.1",
+                   "--t-max", "0.1125", "--surface-model", "skin", "--skin-preset", "dev", "--fine-size", "2",
+                   "--skin-flow-cycles", "1", "--outdir", str(tmp_path)]
+    assert cli.main(argv) == 0
+    (path,) = glob.glob(str(tmp_path / "*.json"))
+    r = json.load(open(path))["results"]
+    rows = list(csv.DictReader(open(path[:-5] + ".csv")))
+    t = [float(row["time_s"]) for row in rows]
+    assert t == pytest.approx([0.0, 0.1, 0.1125]) and float(rows[1]["liquid_depth_max_mm"]) > 0.0
+    assert r["dt_switch_trigger"] == "melting" and r["dt_switch_time_s"] == pytest.approx(0.1)
+    # the history prints nine significant digits
+    assert r["dt_switch_kn"] == pytest.approx(float(rows[1]["knudsen"]), rel=1e-8) and r["dt_switch_kn"] > 0.01
+    assert r["dt_switch_altitude_km"] == pytest.approx(float(rows[1]["altitude_km"]), rel=1e-8)
```

- [ ] **Step 2: The switch's two triggers, and the skins' liquid** -- append to `tests/test_reentry_model_coupled.py`
  and `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_coupled.py
+++ b/tests/test_reentry_model_coupled.py
@@ -249,3 +249,39 @@
         coupled.CoupledSettings(dt=0.5, dt_continuum=1.0)
     with pytest.raises(ValueError):
         coupled.CoupledSettings(dt=0.5, dt_continuum=0.0)
+
+
+def test_the_continuum_switch_latches_at_the_first_liquid_only_where_asked():
+    """Sub-plan 22 Task 3b: with `switch_on_melt` (the layered skins) the continuum step latches at the first macro step
+    at whose end the body holds liquid or whose body Knudsen number is below kn_switch, whichever comes first, and
+    records which (`knudsen`, `melting`, `both`); without it (the element model) the body is never asked and the step
+    latches on the Knudsen number alone, as before; without a continuum step nothing latches."""
+    import types
+
+    class Body:
+        liquid, asked = False, 0
+
+        def liquid_present(self):
+            self.asked += 1
+            return self.liquid
+
+    rarefied, continuum = types.SimpleNamespace(kn=0.02, h=71e3), types.SimpleNamespace(kn=0.005, h=65e3)
+
+    def steps(on_melt, states, dt_continuum=0.0125):
+        run = coupled.CoupledRun(None, Body(), None, coupled.CoupledSettings(dt=0.5, dt_continuum=dt_continuum,
+                                                                             switch_on_melt=on_melt))
+        for k, (liquid, a) in enumerate(states):
+            run.body.liquid = liquid
+            run._switch(0.5 * (k + 1), a)
+        return run
+
+    states = [(False, rarefied), (True, rarefied), (True, continuum)]
+    run = steps(True, states)
+    assert (run.switched_at, run.switch_trigger, run.settings.dt) == ((1.0, 0.02, 71.0), "melting", 0.0125)
+    assert run.body.asked == 2                                   # not asked again once latched
+    run = steps(False, states)                                   # the element model: the Knudsen number alone
+    assert (run.switched_at, run.switch_trigger, run.body.asked) == ((1.5, 0.005, 65.0), "knudsen", 0)
+    assert steps(True, [(False, rarefied), (False, continuum)]).switch_trigger == "knudsen"
+    assert steps(True, [(True, continuum)]).switch_trigger == "both"
+    run = steps(True, states, dt_continuum=None)
+    assert run.switched_at is None and run.settings.dt == 0.5
```

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1279,3 +1279,20 @@
         total = [line.split(",") for line in fh.read().splitlines()[1:]][-spray.N_BINS:]
     np.testing.assert_allclose([float(x[4]) for x in total], n_rec, rtol=1e-8)
     np.testing.assert_allclose([float(x[5]) for x in total], m_rec, rtol=1e-8)
+
+
+def test_the_skins_hold_liquid_once_a_top_cell_is_above_the_liquidus(coarse_sphere_mesh):
+    """Sub-plan 22 Task 3b: `liquid_present`, the melting onset at which a run's continuum step latches, is true once any
+    column's top cell is above the liquidus -- where the liquid depth, the depth option B reads, is not zero -- and not
+    for slurry however warm, nor for a material that does not melt."""
+    b = make(coarse_sphere_mesh, T0=890.0)
+    mat = b.material
+    assert not b.liquid_present() and not (b.skins.depth_above(mat.T_liquidus) > 0.0).any()
+    b.skins.T[3, 0] = mat.T_liquidus - 1e-6                      # slurry a micro-kelvin under the liquidus
+    assert not b.liquid_present()
+    b.skins.T[3, 0] = mat.T_liquidus + 1e-6
+    assert b.liquid_present() and (b.skins.depth_above(mat.T_liquidus) > 0.0).sum() == 1
+    copy = mesh.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
+    cold = skinbody.SkinBody(copy, material.Material.from_drama_json("AA7075_nomelt"), thermal_solver("skfem"), 1.47, DEV,
+                             T0=1200.0, v_hat=(-1.0, 0.0, 0.0))
+    assert not cold.liquid_present()                             # a material that does not melt holds no liquid
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_cli.py
  tests/test_reentry_model_coupled.py tests/test_reentry_model_skinbody.py -q -k "flow_settings_reach or default_and_girin
  or run_names_mark or contradictory or builds_the_skins or bad_skin or sprays_from_its_top or
  latches_at_their_first_liquid or writes_its_columns_frames or melting_run_writes or latches_at_the_first_liquid_only or
  hold_liquid_once"`. Expected: 9 failed, 27 passed (about a minute) --
  - the names: `assert '' == '_skinflow-on'`, `assert 'model_d100.0..._melt-instant' == 'model_d100.0...n_skinflow-on'`
    and the heat-only skin run's name without `_skinflow-off`;
  - `AttributeError: module 'reentry_model.cli' has no attribute 'skin_flow_of'`;
  - `SystemExit: 2` (a melting skin run with no `--skin-flow` given and `--skin-flow-cycles`: the flow was off) in the
    two skin runs, and the heat-only body built where the flow is now the default (`assert (False)`);
  - `TypeError: CoupledSettings.__init__() got an unexpected keyword argument 'switch_on_melt'`;
  - `AttributeError: 'SkinBody' object has no attribute 'liquid_present'`.

- [ ] **Step 4: The switch at the first liquid** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -61,6 +61,9 @@
     dt_continuum: float = None       # s, the macro step from the first step whose body Knudsen number is below kn_switch,
                                      # latched (None: no switch; amendment of 2026-10-07, the continuum step)
     kn_switch: float = KN_BODY_SHOCK  # the body Knudsen number of the switch: the surface flow's own continuum boundary
+    switch_on_melt: bool = False      # the switch also latches at the first macro step at whose end the body holds liquid
+                                      # (`body.liquid_present()`), if that comes first: the layered skins' (sub-plan 22
+                                      # Task 3b, Asha 2026-10-29); the element model's switch is the Knudsen number's alone
 
     def __post_init__(self):
         if self.dt_continuum is not None and not 0.0 < self.dt_continuum <= self.dt:
@@ -75,6 +78,7 @@
         self.settings = settings or CoupledSettings()
         if self.settings.frames_every and not self.settings.output_dir:
             raise ValueError("frames_every > 0 needs an output_dir")
+        self.switched_at = self.switch_trigger = None              # the continuum switch: (t, Kn, altitude km), and why
 
     def loads_at(self, t, y):
         a = self.sim.aero_state(t, y[:3], y[3:])
@@ -106,7 +110,7 @@
     def run(self):
         started = time.perf_counter()
         sim, body, s = self.sim, self.body, self.settings
-        self.switched_at = None                                    # (t, Kn, altitude km) of the continuum switch
+        self.switched_at = self.switch_trigger = None              # (t, Kn, altitude km) of the continuum switch, and why
         self._t0, self._frame_dt, self._next_frame = sim.t, s.frames_every * s.dt, sim.t
         a, loads = self.loads_at(sim.t, sim.y)
         rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
@@ -137,6 +141,8 @@
             t_sw, kn_sw, h_sw = self.switched_at or (None, None, None)
             history.results.update({"dt_continuum_s": s.dt_continuum, "dt_switch_time_s": t_sw, "dt_switch_kn": kn_sw,
                                     "dt_switch_altitude_km": h_sw})
+            if s.switch_on_melt:                                   # what latched it: knudsen, melting or both at once
+                history.results["dt_switch_trigger"] = self.switch_trigger
         return history
 
     def _switch(self, t, a):
@@ -144,12 +150,19 @@
         just given, before any re-evaluation on a surface the step's deaths changed -- has a body Knudsen number below
         kn_switch, the macro step is dt_continuum, latched. Girin's closure begins there and the droplet population needs
         the fine step from there on (plan facts 92-96); before it the default step serves, as the 2026-10-05 harness
-        that this replaces bit for bit had it."""
+        that this replaces bit for bit had it. With `switch_on_melt` (the layered skins, sub-plan 22 Task 3b) it also
+        latches at the first step at whose end the body holds liquid, if that comes first, since the skins spray from
+        their first liquid: the trigger is recorded, `knudsen`, `melting` or `both`."""
         s = self.settings
-        if s.dt_continuum is None or self.switched_at is not None or not (np.isfinite(a.kn) and a.kn < s.kn_switch):
+        if s.dt_continuum is None or self.switched_at is not None:
+            return
+        knudsen = bool(np.isfinite(a.kn) and a.kn < s.kn_switch)
+        melting = bool(s.switch_on_melt and self.body.liquid_present())
+        if not (knudsen or melting):
             return
         s.dt = s.dt_continuum
         self.switched_at = (float(t), float(a.kn), float(a.h) / 1e3)
+        self.switch_trigger = "both" if knudsen and melting else "knudsen" if knudsen else "melting"
 
     def melt_results(self, history):
         c, body = history.columns, self.body
```

- [ ] **Step 5: The skins' liquid** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -108,8 +108,9 @@
 from .body import ThermalBody
 from .thermal import SIGMA_SB
 
-FLOW_DEFAULT = False              # the flow is off until stage 4's spraying gives the liquid a sink: with it on and nothing
-                                  # sprayed, a flight would only pile the liquid up (Task 6c's ruling, 2026-10-22)
+FLOW_DEFAULT = False              # the class's default, the heat-only skins its tests build on: the flow was off until
+                                  # stage 4's spraying gave the liquid a sink (Task 6c's ruling, 2026-10-22); a melting run
+                                  # from the command line flows by default since then (sub-plan 22 Task 3b, Asha 2026-10-29)
 
 GRADIENT_FLOOR = 1.0              # K/m: the smallest gradient along the normal that is divided by. A corner above T_rigid
                                   # whose gradient is below it (a flat or reversed field) is given this one, so its depth
@@ -912,6 +913,13 @@
     def demised(self):
         return self.mass() < self.demise_fraction * self.mass0
 
+    def liquid_present(self):
+        """True once any column's top cell is above the liquidus, that is where the liquid depth (spec 6.7, the depth
+        option B reads for the Rayleigh-Taylor rules and the hold) is not zero: the melting onset at which a run's
+        continuum step latches (`coupled.CoupledSettings.switch_on_melt`; sub-plan 22 Task 3b, Asha 2026-10-29)."""
+        sk = self.skins
+        return bool(self.material.melts and np.any(sk.valid[:, 0] & (sk.T[:, 0] > self.material.T_liquidus)))
+
     def surface_stats(self):
         T, Tt = self.solver.temperature(), self.surface_temperature()
         centre = self.mesh.points.mean(axis=0)
```

- [ ] **Step 6: The flow on by default, and the skins' switch** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -5,8 +5,8 @@
         [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
         [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
         [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off ...]
-                                     (the layered skins, the default with --thermal fem, their flow and spraying with
-                                      --skin-flow on; the element model with --surface-model elements)
+                                     (the layered skins, the default with --thermal fem, their flow and spraying the
+                                      default with --melt on; the element model with --surface-model elements)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
 Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
@@ -113,13 +113,23 @@
     return DEFAULT_SURFACE_MODEL if args.thermal == "fem" else None
 
 
+def skin_flow_of(args):
+    """--skin-flow if given, else the run's own (sub-plan 22 Task 3b; Asha, 2026-10-29): on for the skins' melting runs,
+    which flow and spray by Girin's model, but off where --removal instant is given, the heat-only skins' removal; off
+    for a run that does not melt or does not use the skins, where the flow does not apply."""
+    if args.skin_flow is not None:
+        return args.skin_flow
+    melting_skins = surface_model_of(args) == "skin" and args.thermal == "fem" and args.melt == "on"
+    return "on" if melting_skins and args.removal != "instant" else "off"
+
+
 def removal_of(args):
     """--removal if given, else the run's own (sub-plan 22 Task 3): Girin's, for the element model and for the skins with
     the flow on, whose liquid only spraying removes (spec 10.1: --skin-flow on needs --removal girin); instant removal for
     the heat-only skins, whose liquid leaves as it forms (stage 2's stand-in for flow and spray, spec 11.2)."""
     if args.removal is not None:
         return args.removal
-    return "instant" if surface_model_of(args) == "skin" and args.skin_flow == "off" else "girin"
+    return "instant" if surface_model_of(args) == "skin" and skin_flow_of(args) == "off" else "girin"
 
 
 def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
@@ -146,16 +156,18 @@
 
 
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
-                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=False, flow_cycles=None,
+                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=True, flow_cycles=None,
                      slurry_spray=False):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
     preset's, `_growth<ratio>`, `_sub<ms>ms` and `_deform<mm>mm` for another growth, sub-step or deforming depth,
-    `_noswell` for skins held at the solid's density, `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none), and,
-    with the flow on (sub-plan 21 Task 6c), `_skinflow-on` followed by `_flowcycles<n>` where the flow sub-steps in each
-    heat sub-step are fixed at n rather than chosen each heat sub-step (Task 6d) and by `_slurryspray-on` where the
-    exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3); empty for the defaults."""
+    `_noswell` for skins held at the solid's density, `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none);
+    `_skinflow-off` for the heat-only skins of a melting run (`flow` False: the flow, sub-plan 21 Task 6c, is on by default
+    with melting since sub-plan 22 Task 3b; None for a run that does not melt, which the flow does not concern), and with
+    the flow on `_flowcycles<n>` where the flow sub-steps in each heat sub-step are fixed at n rather than chosen each heat
+    sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3); empty
+    for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -175,10 +187,10 @@
         out += "_fine0"
     elif (fine_band_mm, fine_size_mm) != (FINE_BAND_MM, FINE_SIZE_MM):
         out += "_fine{:g}x{:g}mm".format(fine_band_mm, fine_size_mm)
-    if flow:
-        out += "_skinflow-on"
-        if flow_cycles is not None:
-            out += "_flowcycles{}".format(flow_cycles)
+    if flow is False:
+        out += "_skinflow-off"
+    if flow_cycles is not None:
+        out += "_flowcycles{}".format(flow_cycles)
     if slurry_spray:
         out += "_slurryspray-on"
     return out
@@ -336,11 +348,12 @@
                     help="on: each skin cell at its own phase's density, so melting swells the skin (default). off: every "
                          "cell at the solid's density, the verification setting of levels 2 and 2b, whose references hold "
                          "one density")
-    me.add_argument("--skin-flow", choices=("on", "off"), default="off",
+    me.add_argument("--skin-flow", choices=("on", "off"), default=None,
                     help="with --surface-model skin and --melt on: on, the skins' liquid and slurry flow over the surface "
                          "(sub-plan 21) and Girin's model sprays them from the top cells (sub-plan 22), with --removal "
-                         "girin, whose continuum step (--dt-continuum) follows; off: the heat-only skins, liquid removed "
-                         "as it forms with --removal instant (default)")
+                         "girin, whose continuum step (--dt-continuum) follows (default, unless --removal instant is "
+                         "given); off: the heat-only skins, liquid removed as it forms with --removal instant, named "
+                         "_skinflow-off")
     me.add_argument("--slurry-spray", choices=("on", "off"), default="off",
                     help="with --skin-flow on: on, the exposed slurry sprays too, once no liquid covers it, by its own mode "
                          "(the cascade of the layered-skin design, section 8); it sprays only within about a kelvin of the "
@@ -364,7 +377,8 @@
                     help="effective ratio of specific heats of the Prandtl-Meyer expansion that sets the wall pressure beyond the sonic "
                          "point (default %(default)s; 1.4 frozen air, ~1.15 dissociated -- a factor ~2 on p_w at 90 deg)")
     me.add_argument("--dt-continuum", type=parse_dt_continuum, default=None,
-                    help="macro step [s] from the first step whose body Knudsen number is below --kn-body-shock, latched "
+                    help="macro step [s] from the first step whose body Knudsen number is below --kn-body-shock, or, "
+                         "with the skins, from the first step at whose end they hold liquid if that comes first, latched "
                          "(default {} with --removal girin, where the droplet population converges; off with instant "
                          "removal). off: the --dt step throughout, as every run before 2026-10-07".format(DEFAULT_DT_CONTINUUM))
     me.add_argument("--kn-body-shock", type=float, default=surface_flow.KN_BODY_SHOCK,
@@ -506,7 +520,8 @@
                                        skin_suffix=skin_name_suffix(args.skin_min_thickness, args.skin_cell, args.skin_growth,
                                                                     args.skin_substep, args.skin_preset, args.deform_depth,
                                                                     args.skin_swell == "on", args.fine_band, args.fine_size,
-                                                                    args.skin_fine_zone, args.skin_flow == "on",
+                                                                    args.skin_fine_zone,
+                                                                    skin_flow_of(args) == "on" if args.melt == "on" else None,
                                                                     args.skin_flow_cycles, args.slurry_spray == "on"))
 
 
@@ -517,9 +532,10 @@
     np.random.seed(args.seed)
     if args.surface_model == "skin" and args.thermal != "fem":
         parser.error("--surface-model skin needs --thermal fem")
-    # the surface model and the removal not given take the run's own (sub-plan 22 Task 3): the skins for a thermal run,
-    # and Girin's removal except for the heat-only skins
-    args.surface_model, args.removal = surface_model_of(args), removal_of(args)
+    # the surface model, the flow and the removal not given take the run's own (sub-plan 22 Tasks 3 and 3b): the skins for
+    # a thermal run, their flow on when they melt unless instant removal is given, and Girin's removal except for the
+    # heat-only skins
+    args.surface_model, args.skin_flow, args.removal = surface_model_of(args), skin_flow_of(args), removal_of(args)
     for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                          ("--cadence", args.cadence), ("--t-max", args.t_max)):
         if value <= 0.0:
@@ -575,8 +591,8 @@
         parser.error("--skin-flow on keeps the liquid, which only spraying removes: it needs --removal girin, not instant "
                      "removal, which stands in for flow and spraying in the heat-only skins (spec 11.2)")
     if args.surface_model == "skin" and args.melt == "on" and args.skin_flow == "off" and args.removal != "instant":
-        parser.error("Girin's spraying from the skins needs --skin-flow on: the heat-only skins (--skin-flow off, the "
-                     "default) remove their liquid as it forms, --removal instant")
+        parser.error("Girin's spraying from the skins needs their flow: the heat-only skins (--skin-flow off) remove their "
+                     "liquid as it forms, --removal instant")
     if args.slurry_spray == "on" and args.skin_flow != "on":
         parser.error("--slurry-spray on needs --skin-flow on (with --surface-model skin and --melt on): the exposed slurry "
                      "sprays by the skins' spraying from their top cells")
@@ -614,7 +630,8 @@
         frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
         run = coupled.CoupledRun(sim, the_body, heating_model,
                                  coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk"),
-                                                         dt_continuum=dt_continuum, kn_switch=args.kn_body_shock))
+                                                         dt_continuum=dt_continuum, kn_switch=args.kn_body_shock,
+                                                         switch_on_melt=args.surface_model == "skin"))
         history = run.run()
         thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
         thermal_info["dt_continuum_s"] = dt_continuum
```

- [ ] **Step 7: Run them to see them pass.** Expected: `tests/test_reentry_model_cli.py` 83 passed (the melting
  trigger's run and the spraying run about 10-20 s each), `tests/test_reentry_model_coupled.py` 10,
  `tests/test_reentry_model_skinbody.py` 35 (2 deselected); the unit tier 523 passed, 1 skipped (Task 3's 519 plus
  four: the melting trigger's two, the skins' liquid, and a third heat-only contradiction).

- [ ] **Step 8: The measurements** (in `measure/task3b/`, on mains power, one heavy run at a time): the element model's
  short case, `../task1/cmp/run_case.sh <code> runs line3b_short --surface-model elements` against Task 3's
  `line3_short` copied into `runs/`, `"$PY" -B ../task1/cmp/compare_runs.py runs line3_short line3b_short` and `"$PY" -B
  ../task3/element_files.py runs line3_short line3b_short`; the two triggers, `./run_evidence.sh` (about a minute and a
  half), then `"$PY" -B onset_check.py evidence/<run>.json` and `"$PY" -B smoke_check.py evidence/<run>.json` for each.

- [ ] **Step 9: Commit.**

```bash
git add reentry_model/cli.py reentry_model/coupled.py reentry_model/skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_coupled.py tests/test_reentry_model_skinbody.py
git commit -m "Skins: the flow on by default for a melting skin run and the continuum step from the skins' first liquid (Asha 2026-10-29; spec 10.1) -- --skin-flow not given on for the skins' melting runs unless instant removal is given, the heat-only skins named _skinflow-off and a flow run's name losing _skinflow-on; on the skins the continuum step also latches at the first macro step at whose end a top cell is above the liquidus, or at the Knudsen switch if that comes first, the results recording the trigger (knudsen, melting or both); the element model's switch, outputs and results unchanged to the last bit (sub-plan 22 Task 3b)"
```

---

### Task 4: The outputs of a spraying skin flight, Register's release angle, and the records on disk

**Files:**
- Modify: `reentry_model/spray.py` (`SOURCE_COLUMNS` + `release_angle_register_deg`; `source_block(...,
  release_angle=None)`, `source_rows`; `source_table` and `source_columns` reading a spool; `SourceSpool`;
  `register_angle`; `enclosed_centroid`; `RELEASE_ANGLE_COLUMNS`, `release_angle_stats`; the module docstring),
  `reentry_model/body.py` (the element model's rows' release angle and its history's three columns),
  `reentry_model/skinbody.py` (`SHEARED_FRACTION`, `sheared_depth`; `SprayPlan.area`, `.centre`, `.v_hat`, `.spraying`,
  `.release_angle`; `girin_plan(..., centre=None, v_hat=None)`; `SkinBody(..., spool=None)`, `.sprayed_mass`,
  `.n_source_rows`, `.step_spray`, `._step_columns`, `._depth_ratios`; the history's new columns; the module
  docstring), `reentry_model/coupled.py` (`SKINBODY_COLUMNS`, `SUMMARY_COLUMNS`, `SKIN_SUMMARY_COLUMNS`, `END_SNAP`;
  `CoupledSettings.snap_to_t_max`; the loop's last step; `skin_results`; `write_particles`),
  `reentry_model/trajectory.py` (`Simulator.advance(dt, snap=0.0)`), `reentry_model/movemesh.py`
  (`RebuildResult.mass_residual`, `MovingMesh.rebuild_mass_residual_total`), `reentry_model/cli.py` (`SOURCE_SPOOL`;
  `build_thermal(..., spool=None)`; `cmd_run`'s spool, last step and printed summary)
- Test: `tests/test_reentry_model_spray.py` (four added, one changed), `tests/test_reentry_model_trajectory.py` (one
  added), `tests/test_reentry_model_coupled.py` (two added), `tests/test_reentry_model_movemesh.py` (one added),
  `tests/test_reentry_model_skinbody.py` (one added, three changed), `tests/test_reentry_model_melting.py` (one changed),
  `tests/test_reentry_model_cli.py` (two changed), `tests/test_reentry_model_fenicsx.py` (one changed)

**Interfaces:**
- Consumes: Task 2's `SprayPlan` records and `source_block`, `SkinBody.source_rows`, `coupled.write_particles`; Task
  3's command line and `particles_summary` recorded as null; Task 3b's `CoupledSettings`; sub-plan 21 Task 7's
  `released_mass_kg` and Pi; the element model's `MeltingBody.source_rows`, `melt_stats` and its surface.
- Produces: `spray.SOURCE_COLUMNS` (27: + `release_angle_register_deg` [deg]); `spray.register_angle(points, centre,
  v_hat) -> (n,)` [deg]; `spray.enclosed_centroid(points, faces, normals=None) -> (3,)` [m];
  `spray.RELEASE_ANGLE_COLUMNS`, `spray.release_angle_stats(dm, angle) -> (mean, min, max)`; `spray.source_block(...,
  release_angle=None)`, `spray.source_rows(..., release_angle=None)`; `spray.SourceSpool(directory)` with `.append(block)`,
  `.columns()`, `.n_rows`, `len` (blocks), `.remove()`; `spray.source_columns(rows) -> {name: array}`;
  `skinbody.SHEARED_FRACTION` (0.1), `skinbody.sheared_depth(sk, fraction) -> (n,)` [m]; `SprayPlan.spraying()`,
  `.release_angle()`; `SkinBody(..., spool=None)`, `.sprayed_mass` [kg], `.n_source_rows`; the skin history's
  `released_step_mass_kg`, `n_released`, `r_median_um`, `r_max_um`, `theta_cr_deg`, `spraying_area_m2`, the three
  `release_angle_register_*_deg`, `b_over_delta_m_median`, `sheared_depth_over_delta_m_median`,
  `rebuild_mass_residual_kg`, `rebuild_mass_residual_last_kg`; the element history's three release-angle columns; the
  skin results' `sprayed_mass_kg`, `n_released`, `n_source_rows`, `wall_flow_s`, `rebuild_mass_residual_kg`;
  `coupled.SUMMARY_COLUMNS`, `coupled.SKIN_SUMMARY_COLUMNS`, `coupled.END_SNAP` (1e-6), `CoupledSettings(...,
  snap_to_t_max=False)`; `trajectory.Simulator.advance(dt, snap=0.0)`; `movemesh.RebuildResult.mass_residual` [kg],
  `MovingMesh.rebuild_mass_residual_total`; `cli.SOURCE_SPOOL`, `build_thermal(..., spool=None)`. Every existing column,
  row, field, result and file of the element model is unchanged to the last bit.

The outline's item 4 (outputs) with the 23-notes' first item, which Asha made this task's on 2026-10-29, and the three
loose ends of Tasks 1a-3: the summary Task 2 left (its ruling 10), the vanishing last macro step of fact 179, and the
FEniCSx tier's failure since Task 1a.

**The release angle** (Asha, 2026-10-29; `docs/superpowers/plans/melt-spraying-subplans/23-notes.md`, item 1; Register,
P. J., et al. 2020, Icarus 337, 113468, section 4.2). The angle between the larger body's trajectory axis and the line
joining its centroid to the released fragment's, from downstream: theta_R = 180 degrees - angle(v_hat, p - c), 180 at
the stagnation point, 90 at the equator, 0 directly behind (`spray.register_angle`, from atan2 of the vector's parts
across and along v_hat). It is computed from the vector p - c, never from the source table's `theta_deg`, which is the
surface normal's angle: on the element model's short case the two differ by -65 to +88 degrees on the staircase its
element deaths leave, on the skins' smoothed corner normals by -7.3 to +5.0. p is the row's own position; c is the
centroid of the volume the body's outer surface encloses (`spray.enclosed_centroid`: signed tetrahedra from the points'
mean, the triangles turned by their outward normals where given), from the surface the positions are on and at their
instant -- the element model's surface at its spray step, the skins' top surface where the plan is evaluated, the macro
step's start (`girin_plan(..., centre=, v_hat=)`); v_hat is the body's direction through the air, its fixed attitude
(ruling 1). Each model's rows carry it through the shared `spray.source_block(..., release_angle=)`, the last column of
the source table, and each history its mean by mass, least and largest over the step's rows
(`spray.release_angle_stats`, `spray.RELEASE_ANGLE_COLUMNS`), nan where the step released nothing: the element model's
from the patches that released (`MeltingBody._film_and_spray`), the skins' from the step's block. Register's wake zone
(AFFIx) is not an output (Asha: post-processing).

**The summary** (Task 2's ruling 10). `particles_summary.csv`'s columns are the element model's per-step history
columns (`coupled.SUMMARY_COLUMNS`), and the skin history now has them under the same names and meanings: `n_released`,
the droplets since time zero, as the element model's; `r_median_um` and `r_max_um` over the step's records; and
`theta_cr_deg` and `spraying_area_m2`, the smallest angle and the top area of the columns whose plan sprays
(`SprayPlan.spraying`: a mode with a release rate, as the element model's `unstable` is its rate's). Its
`released_mass_kg` is since time zero (sub-plan 21 Task 7) where the summary's is the step's, so the skin history gains
`released_step_mass_kg`, what the accepted parts released in the macro step, and keeps its own column as it was; the
writer reads the summary's released mass from it (`coupled.SKIN_SUMMARY_COLUMNS`), and where the element model has its
film, which the skins have not, the skins' mean liquid and flowing depths (ruling 2). Neither model's existing columns
change meaning. The results add the element model's `sprayed_mass_kg` (the records' mass, the books' released mass to
round-off), `n_released` and `n_source_rows` (the rows written).

**The rest of spec 10.2.** Audited against `SKINBODY_COLUMNS`, `write_skin_frame` and `skin_results`: every history
item but three was written, every frame field was, and the results lacked the flow's wall time. Built: b over delta_m
and the resolved sheared depth over delta_m, medians over the wet windward columns -- windward, with flowing metal,
under Girin's closure, delta_m from the macro step's surface flow -- at the step's end
(`SkinBody._depth_ratios`); the resolved sheared depth (spec 6.7) the depth at which the computed speed falls to a tenth
of the top cell's, between the cells' centres (`skinbody.sheared_depth`, `SHEARED_FRACTION`; ruling 6); each rebuild's
mass residual beside its energy residual (`RebuildResult.mass_residual`, the density times the old volume less the new,
uncorrected; `rebuild_mass_residual_kg` summed and `_last_kg`); and the results' `wall_flow_s` and
`rebuild_mass_residual_kg`. One item differed from the spec's word: the spec listed "shortfall mass (cumulative)",
but `shortfall_kg` is the current shortfall below the minimum at the step's end, not a sum; Asha corrected the spec's
word on 2026-10-31 and the column is unchanged (ruling 10).

**The vanishing last macro step** (fact 179). The loop asks `Simulator.advance` for `min(t + dt, t_max)` again and again;
the floating-point sum falls short of t_max -- ten 0.0125 s steps reach 0.12499999999999999 s, ten of 0.1 s
0.9999999999999999 s -- so one more step of 1.4e-17 s or 1.1e-16 s follows, a whole macro step of the body (its parts,
its passes, a 3D solve), over which the droplet-momentum ratio, a ratio of two quantities each nearly nothing, read 5.5
and 8.6. The element model has the same step: on its short case the last two rows are both at 15 s (printed), the last
releasing 1.5e-13 kg. The loop is fixed for the skins: with `CoupledSettings.snap_to_t_max`, which the command line
sets for them, a step that would end within `END_SNAP` (1e-6) of a macro step short of t_max ends at t_max
(`Simulator.advance(dt, snap=)`), so there is no vanishing step; without it, the default and the element model's, the
arithmetic is what it was (`snap` 0 leaves `t_target` as it was), its outputs unchanged to the last bit (ruling 4).

**The records on disk** (Task 2's decision 2). With the continuum step from the first liquid (Task 3b), the 100 mm
flight runs 12.5 ms steps from about 24.5 s to 120 s: 7 640 macro steps, and at Task 2's 4 600-9 200 rows a step (slurry
spraying off and on) 35-70 million rows of 27 columns, 216 bytes a row, so 7.6-15.2 GB as arrays and twice that while
`source_table` joins them at the end -- more than the laptop's 24 GB. So the command line gives a spraying skin body a
spool, `<run>/source_rows/` (`cli.SOURCE_SPOOL`): each macro step's block is appended to one file per column
(`spray.SourceSpool`), and the coupled writer reads the columns back from disk (`spray.source_columns`), each once and
contiguous, writes `particles.npz`, the summary and the size distribution with the same arithmetic as from memory, and
the command line then removes the spool. Nothing is combined over more than a macro step, so spec 8 needs no new
sentence (ruling 9). The element model's rows stay in memory as before, and its files are byte for byte.

**The FEniCSx tier's failure** (`test_the_thick_film_flux_matches_the_skfem_backend`, failing since Task 1a). The
element model's books are one `MeltingBody`'s on both backends; printed term by term after every step
(`measure/task4/fenicsx_books.py`) they close to 8.9e-12 of the heat absorbed on the FEniCSx backend and wander from
2.2e-9 to 1.17e-8 on scikit-fem's, ending at 1.05e-8 (5.1e-9 on sub-plan 21's code, which already reached 1.07e-8 at its
sixth step), with the spraying amendments of 2026-10-22, which spray 2.2 % more there, taking another path. With
scikit-fem's conjugate gradients at 1e-13 instead of the default 1e-10, its Newton iterations the same in every step, the
residual is 4.5e-12; a tighter Newton tolerance alone leaves 5.8e-9. So it is the scikit-fem backend's linear solves'
tolerance, not an error in the books, and the test now runs both backends' gradients at 1e-13 (as the skin runs do,
`skinbody.SKIN_CG_TOL`), its 1e-8 bound kept (ruling 5).

**Measured** (2026-10-31, `prototype/work-2026-10-28-skin-spraying/measure/task4/`; facts 183-187):
- *The element model unchanged* (`runs/line4_short` against Task 3b's `line3b_short`; `../task1/cmp/compare_runs.py`,
  `element_files.py`). The CLI test's short melting case with `--surface-model elements`: its 93 history columns first,
  in their order, the same text in all 76 rows, the three release-angle columns after them; the 26 old fields of its
  25 266 source rows identical, the 27th the release angle; the summary and the size distribution byte for byte; the
  settings and the 61 results but the run time. Its rows' angles 92.0/116.1/177.8 degrees (minimum/median/maximum);
  each step's history mean, least and largest the rows' within 4.7e-9 (nine printed digits), nan on the 4 steps that
  released nothing; 180 - theta_deg differs from the angle by -64.8 to +87.6 degrees (median 7.1), the element deaths'
  staircase. Its last two rows print as 15 s, the last releasing 1.5e-13 kg: its vanishing step, kept.
- *The skin body unchanged* (`bitid.py`, `run_bitid.sh`). Task 2's spraying test body with slurry spraying on, pyamg's
  draws seeded: its 3D field, mesh, skins (mass, enthalpy, temperature, velocity, counts), books, histograms and the 26
  old columns of its 759 rows are Task 3b's in all 12 arrays to the last bit, the history read after every step or not.
- *The skin body's outputs* (`body_outputs.py`; the same body). The second step's 759 rows leave at 111.5/135.4/169.9
  degrees (134.9 by mass), within 7.3 degrees of 180 - theta_deg on the corners' smoothed normals; the columns within
  10.9 degrees of the nose do not spray (theta_cr_deg). The history: 5.769 g released in the step, 89 173 droplets, radius
  median 174 um and largest 3.53 mm, 9.6e-3 m2 spraying; b over delta_m 5.5 at the first step's end and 10.7 at the
  second's, the resolved sheared depth over delta_m 0.46 and 1.91 (medians over the wet windward columns). At the step's
  end the body's mass centre (the 3D model and the skins' cells) lay 0.083 mm from the geometric centroid, which would
  move the angles by at most 0.088 degrees; the centroid itself moved back 0.33 mm over the step.
- *The smoke flight* (`run_smoke.sh`, `smoke_check.py`; Task 3's 65 km case, the default flow, slurry spraying off).
  9 macro steps, 25 ms and eight of 12.5 ms, the last ending at 0.125 s exactly, against Task 3b's 10 whose last was
  1.4e-17 s; every old history column of the first 9 rows the same text as Task 3b's, and the last step its same step's
  but the interface mismatch (1.3e-9 relative) and the momentum residual (2.2e-17 against 4.0e-17, round-off both); Pi
  at most 0.083 where Task 3b's vanishing step read 5.50. 1.95571981656 g released in both, in 1 005 rows against 1 196:
  the vanishing step had written 191 rows holding 3.4e-17 kg in all. The summary: 10 rows under the skins' header; each
  step's released mass, droplets since time zero, radius median and largest and release angle's mean, least and largest
  the rows' within 3.4e-9 (nine printed digits); its released mass the step's; the size distribution the rows'
  histogram within 5.4e-10 of the totals; the spool removed. Angles 123.5/145.5/168.7 degrees (148.9 by mass),
  theta_cr_deg 11.5; b over delta_m 2.8 rising to 7.6, the resolved sheared depth over delta_m 2.5 falling to 0.43; no
  rebuild. With slurry spraying on: 2.983 g in 1 811 rows against 2 223 (412 vanishing rows, 8.3e-17 kg), Pi at most
  0.079 against 8.59, the summary and distribution within 3.9e-9 and 2.6e-10. About 18 s each, 1.0 GB at most.
- *The records on disk* (`spool_cost.py`, 3.0 million rows, blocks of 9 231). Appending a block costs 1.4-2.9 ms (a
  melting 12.5 ms macro step's skins take 1.5-6.2 s, fact 179); the writer from the spool and from blocks in memory give
  the same `particles.npz` fields and size distribution byte for byte, from the spool with at most 123 MB of its own
  allocations (41 bytes a row), from memory with the rows' 616 MB held and 739 MB more (453 bytes a row). For the 100 mm
  flight's 35-70 million rows: 7.6-15.2 GB on disk and about 1.4-2.9 GB at the writer, against 16-32 GB from memory.
- *The FEniCSx tier's failure* (`fenicsx_books.py`, `run_fenicsx_books.sh`): the energy residual over the absorbed heat
  after ten 0.05 s steps --

  | code | scikit-fem, default | FEniCSx, default | scikit-fem, CG 1e-13 | FEniCSx, CG 1e-13 |
  |---|---|---|---|---|
  | sub-plan 21 (602b262) | 5.1e-9 (1.07e-8 at step 6) | 1.1e-11 | 1.6e-12 (Newton 1e-10) | -- |
  | Task 3b | 1.05e-8 (1.17e-8 at step 9) | 8.9e-12 | 4.5e-12 (Newton 1e-6), 6.3e-13 (1e-10) | 2.0e-14 |

  and scikit-fem's at a Newton tolerance of 1e-8 with its default gradients 5.8e-9: the residual is the linear solves'.
  The amended code sprays 2.2 % more there (removed enthalpy 1.375e5 J against 1.345e5 J). With both solvers at 1e-13 the
  test passes in 22 s.
- *The rebuild's mass residual*: zero to the last bit on each of two rebuilds of the coarse test sphere (10 191 and
  9 955 elements, their volumes equal in every bit).

**Rulings** (2026-10-31; each what, why, cost if wrong; Asha accepted rulings 1-5 as built on 2026-10-31):
1. *The release angle's centroid is the geometric one*, of the volume the outer surface encloses, from the surface and
   at the instant of the rows' positions. Why: Register's bodies are of one density, whose centroid is the centre of
   mass; it is defined alike for both models from the surface the rows' positions are on, where the element model's
   mass centre (`mass_centre`) and a skin body's (the 3D model and the skins' cells) are computed differently. Cost:
   on the 4 mm test sphere's spraying step the mass centre lay 0.083 mm from the geometric one (the hot liquid at the
   nose is lighter), moving the angles by at most 0.088 degrees.
2. *A skin flight's summary: the element model's droplet columns, its released mass the step's, the skins' liquid and
   flowing depths where the element model has its film.* Why: the same name must mean the same thing in both models'
   files, the skins have no film, and their mean liquid and flowing depths are the history's columns of the reservoir
   the release takes. Cost: a reader comparing the two summaries finds two columns named differently; the options
   are the element model's film names left empty (nan) or given the skins' liquid.
3. *The step's released mass a new column, `released_step_mass_kg`*, the skins' cumulative `released_mass_kg` kept.
   Why: neither model's columns change meaning. Cost: the skin history has both.
4. *The loop's fix for the skins only* (`snap_to_t_max`). Why: on the element model the same fix moves its last step's
   end by 2e-14 s and removes a row, which changes its outputs; on the skins it removes a wasted macro step. Cost: the
   element model keeps its vanishing step (and its own Pi-free columns); the fix can be given to both at a later
   comparison's cost.
5. *The FEniCSx test runs both solvers' gradients at 1e-13*, its bound kept, rather than a wider bound. Why: it then
   checks the books to about 1e-11 on both backends, a tighter test than the bound's 1e-8, and the evidence shows the
   residual is the solver's. Cost: the test compares the backends at a tolerance neither flight uses; the alternative
   is a bound of 3e-8 for scikit-fem's at its default.
6. *The resolved sheared depth at a tenth of the top cell's speed.* Why: spec 6.7 leaves the fraction to be set; on
   Girin's linear profile, which reaches nothing at delta_m, a tenth gives 0.9 delta_m, so the ratio reads near 1 where
   the resolved profile is Girin's. Cost: another fraction scales the column (a half gives 0.5 delta_m). Replaced by
   Task 4b (Asha, 2026-10-31): the straight-line equivalent, with no fraction.
7. *"Wet windward" for the skins: theta <= 90 degrees, a flowing depth above zero and Girin's delta_m*, at the step's
   end. Why: b is the flowing depth (spec 6.7), and delta_m exists only under his closure. Cost: none.
8. *theta_cr_deg and spraying_area_m2 over the columns whose plan sprays, r_median_um and r_max_um over the step's
   rows.* Why: the element model's are over its patches with a rate and its rows. Cost: none.
9. *The spool always for a spraying skin run from the command line*, removed once written; a skin body built in Python
   keeps its rows in memory unless given a directory. Why: one path for every flight; the tests read either. Cost: a
   failed run leaves its partial records on disk.
10. *`shortfall_kg` unchanged* though spec 10.2 calls it cumulative. Why: it is a state (how far the columns are below
    their minimum), whose sum over steps means nothing, and changing it changes an existing column. Asha, 2026-10-31:
    the spec's word is corrected, the column unchanged. Cost: none.
11. *The element model's history gains three columns at its end, its rows one at their end.* Why: Asha's requirement
    for both models; appended, every existing column and field keeps its position and value. Cost: the history file is
    three columns wider, the same column by column.

- [ ] **Step 1: The release angle, its statistics and the spool** -- in `tests/test_reentry_model_spray.py`:

```diff
--- a/tests/test_reentry_model_spray.py
+++ b/tests/test_reentry_model_spray.py
@@ -588,8 +588,9 @@
     solid share and the enthalpy per kilogram they carry, 0 and nan where the caller gives none; the droplets' density,
     viscosity and surface tension may be one per patch, a slurry droplet's Weber and Ohnesorge numbers being its own
     phase's (the 2026-10-07 design 8.4); `source_block` is `source_rows` as an array, and `spray.source_table` reads a
-    body's rows or blocks of rows alike. The element model's form is unchanged in every other column."""
-    assert spray.SOURCE_COLUMNS[-2:] == ["solid_fraction", "h_kg"] and len(spray.SOURCE_COLUMNS) == 26
+    body's rows or blocks of rows alike. The element model's form is unchanged in every other column. (Sub-plan 22 Task
+    4 adds a 27th column after them, the release angle.)"""
+    assert spray.SOURCE_COLUMNS[24:26] == ["solid_fraction", "h_kg"] and len(spray.SOURCE_COLUMNS) == 27
     flow = FakeFlow(3)
     delta_m, _ = spray.melt_layer(flow, LIQ)
     b = np.full(3, 1e-3)
@@ -599,12 +600,12 @@
              res, b)
     C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
     plain = np.array(spray.source_rows(*where, LIQ, FakeState()))
-    assert plain.shape == (2, 26) and (plain[:, C["solid_fraction"]] == 0.0).all() and np.isnan(plain[:, C["h_kg"]]).all()
+    assert plain.shape == (2, 27) and (plain[:, C["solid_fraction"]] == 0.0).all() and np.isnan(plain[:, C["h_kg"]]).all()
     slurry = material.LiquidProperties(np.array([2400.0, 2500.0, 2550.0]), np.array([1.3e-3, 0.02, 0.05]),
                                        np.array([0.86, 0.80, 0.80]))
     block = spray.source_block(*where, slurry, FakeState(), solid_fraction=np.array([0.0, 0.1, 0.03]),
                                h_kg=np.array([1.1e6, 1.0e6, 0.9e6]))
-    assert block.shape == (2, 26) and block[:, C["solid_fraction"]].tolist() == [0.0, 0.03]
+    assert block.shape == (2, 27) and block[:, C["solid_fraction"]].tolist() == [0.0, 0.03]
     assert block[:, C["h_kg"]].tolist() == [1.1e6, 0.9e6]
     np.testing.assert_array_equal(block[:, :C["we_d"]], plain[:, :C["we_d"]])
     np.testing.assert_array_equal(block[0, :C["solid_fraction"]], plain[0, :C["solid_fraction"]])   # the liquid's patch
@@ -613,5 +614,100 @@
     assert block[1, C["oh"]] == pytest.approx(0.05 / np.sqrt(2550.0 * 0.80 * 2.0 * r), rel=1e-14)
     np.testing.assert_array_equal(spray.source_block(*where, LIQ, FakeState()), plain)        # rows as an array
     table = spray.source_table([block, block[:1]])
-    assert table.shape == (3, 26) and np.array_equal(table[:2], block, equal_nan=True)
-    assert np.array_equal(spray.source_table(plain.tolist()), plain, equal_nan=True) and spray.source_table([]).shape == (0, 26)
+    assert table.shape == (3, 27) and np.array_equal(table[:2], block, equal_nan=True)
+    assert np.array_equal(spray.source_table(plain.tolist()), plain, equal_nan=True) and spray.source_table([]).shape == (0, 27)
+
+
+# -- the outputs of a spraying flight (sub-plan 22 Task 4) -------------------------------------------------------------
+def test_register_s_release_angle_is_180_at_the_stagnation_point_90_at_the_equator_and_0_behind_from_the_actual_vector():
+    """Register et al. (2020) section 4.2 (Asha, 2026-10-29; the 23-notes' first two tests): the angle between the larger
+    body's trajectory axis and the line from its centroid to the release point, from downstream: 180 degrees at the
+    stagnation point, 90 at the equator, 0 directly behind, wherever the body is and whichever way it flies. It is the
+    position's angle, not the surface normal's: on the flat face of a flattened body the normal points along the flight
+    (theta_deg 0) but a point 30 mm off the axis 20 mm ahead of the centroid is at 180 - atan(30/20); and a body whose
+    nose has receded, its centroid moved back, sees the same point at a larger angle."""
+    v = np.array([0.0, 0.6, 0.8])
+    c = np.array([0.01, -0.02, 0.03])
+    side = np.cross(v, [1.0, 0.0, 0.0])
+    side /= np.linalg.norm(side)
+    points = c + 0.05 * np.array([v, side, -v, (v + side) / np.sqrt(2.0)])
+    np.testing.assert_allclose(spray.register_angle(points, c, v), [180.0, 90.0, 0.0, 135.0], atol=1e-12)
+    np.testing.assert_allclose(spray.register_angle(points, c, 3.0 * v), [180.0, 90.0, 0.0, 135.0], atol=1e-12)
+    flat = np.array([[0.02, 0.03, 0.0]])                        # on a flat front face x = 20 mm, 30 mm off the axis
+    assert spray.register_angle(flat, np.zeros(3), [1.0, 0.0, 0.0])[0] == pytest.approx(180.0 - np.degrees(np.arctan2(0.03, 0.02)))
+    receded = spray.register_angle(flat, [-0.01, 0.0, 0.0], [1.0, 0.0, 0.0])[0]
+    assert receded == pytest.approx(180.0 - np.degrees(np.arctan2(0.03, 0.03))) and receded > 180.0 - 56.4
+
+
+def test_the_enclosed_centroid_is_the_volume_s_whichever_way_its_triangles_turn(coarse_sphere_mesh):
+    """`spray.enclosed_centroid`: the centroid of the volume a closed triangulated surface encloses, by the divergence
+    theorem -- the 100 mm test sphere's moved anywhere is its tetrahedra's volume-weighted centroid to round-off; turning
+    every triangle gives the same point, and with one outward normal per triangle a surface whose triangles are mixed
+    gives it too (the element model's surface)."""
+    m = coarse_sphere_mesh
+    s = m.surface()
+    shift = np.array([0.3, -0.1, 0.25])
+    P = m.points + shift
+    w = m.element_volumes()
+    exact = (w[:, None] * P[m.tets].mean(axis=1)).sum(axis=0) / w.sum()
+    got = spray.enclosed_centroid(P, s.faces)
+    np.testing.assert_allclose(got, exact, rtol=0.0, atol=1e-15)
+    np.testing.assert_allclose(spray.enclosed_centroid(P, s.faces[:, [0, 2, 1]]), got, rtol=0.0, atol=1e-15)
+    mixed = s.faces.copy()
+    mixed[::3] = mixed[::3][:, [0, 2, 1]]
+    np.testing.assert_allclose(spray.enclosed_centroid(P, mixed, s.normals), got, rtol=0.0, atol=1e-15)
+
+
+def test_the_release_angle_rides_in_the_source_table_and_a_step_s_statistics_are_by_mass():
+    """The source table's last column is the release angle, nan where the caller gives none; a step's history columns
+    are its mass-weighted mean, least and largest over the records that released, and nan where none did: the history's
+    mean equals the mean recomputed from the step's rows (the 23-notes' third test)."""
+    assert spray.SOURCE_COLUMNS[-1] == "release_angle_register_deg" and spray.RELEASE_ANGLE_COLUMNS == (
+        "release_angle_register_mean_deg", "release_angle_register_min_deg", "release_angle_register_max_deg")
+    flow = FakeFlow(3)
+    delta_m, _ = spray.melt_layer(flow, LIQ)
+    b = np.full(3, 1e-3)
+    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.full(3, 5.0), np.array([True, False, True]), 0.5,
+                                         np.full(3, 1e-5), b * LIQ.rho * 1e-5)
+    where = (10.0, 71e3, 7200.0, np.radians([30.0, 45.0, 60.0]), np.zeros((3, 3)), np.tile([0.0, 1.0, 0.0], (3, 1)), flow,
+             res, b, LIQ, FakeState())
+    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
+    angle = np.array([150.0, 135.0, 120.0])
+    rows = spray.source_block(*where, release_angle=angle)
+    assert rows[:, C["release_angle_register_deg"]].tolist() == [150.0, 120.0]
+    assert np.isnan(spray.source_block(*where)[:, C["release_angle_register_deg"]]).all()
+    np.testing.assert_array_equal(spray.source_block(*where)[:, :-1], rows[:, :-1])
+    assert np.array_equal(spray.source_rows(*where, release_angle=angle), rows, equal_nan=True)
+    mean, lo, hi = spray.release_angle_stats(rows[:, C["dm_kg"]], rows[:, C["release_angle_register_deg"]])
+    dm = res.dm[[0, 2]]
+    assert mean == (dm * [150.0, 120.0]).sum() / dm.sum() and (lo, hi) == (120.0, 150.0)
+    assert mean == spray.release_angle_stats(res.dm, angle)[0]                     # the rows' own mean, to the last bit
+    assert all(np.isnan(x) for x in spray.release_angle_stats(np.zeros(3), angle))
+    assert all(np.isnan(x) for x in spray.release_angle_stats([], []))
+
+
+def test_a_spool_on_disk_gives_back_the_blocks_it_was_given_to_the_last_bit(tmp_path):
+    """`spray.SourceSpool` (sub-plan 22 Task 4): the skins' blocks of rows written to disk as they come, one file per
+    column, read back column by column -- the source table's columns to the last bit, nan included, in the order the
+    blocks came; `len` counts the blocks given rows, as the in-memory list does, `n_rows` the rows; `source_table` and
+    `source_columns` read a spool as they read a list; an empty spool gives empty columns; `remove` deletes it."""
+    rng = np.random.default_rng(1)
+    blocks = [rng.normal(size=(k, len(spray.SOURCE_COLUMNS))) for k in (5, 1, 7)]
+    blocks[1][0, 26] = np.nan
+    spool = spray.SourceSpool(str(tmp_path / "spool"))
+    assert len(spool) == 0 and all(len(c) == 0 for c in spool.columns().values()) and spray.source_table(spool).shape == (0, 27)
+    for blk in blocks[:2] + [np.zeros((0, 27))] + blocks[2:]:
+        spool.append(blk)
+    assert len(spool) == 3 and spool.n_rows == 13
+    table = spray.source_table(blocks)
+    cols = spray.source_columns(spool)
+    assert list(cols) == spray.SOURCE_COLUMNS
+    assert all(np.array_equal(cols[name], table[:, i], equal_nan=True) for i, name in enumerate(spray.SOURCE_COLUMNS))
+    assert np.array_equal(spray.source_table(spool), table, equal_nan=True)
+    listed = spray.source_columns(blocks)
+    assert all(np.array_equal(listed[name], cols[name], equal_nan=True) for name in spray.SOURCE_COLUMNS)
+    with pytest.raises(ValueError):
+        spool.append(np.zeros((2, 26)))
+    del cols
+    spool.remove()
+    assert not (tmp_path / "spool").exists()
```

- [ ] **Step 2: The loop's last step, the writer and the rebuild's mass** -- append to
  `tests/test_reentry_model_trajectory.py`, `tests/test_reentry_model_coupled.py` and
  `tests/test_reentry_model_movemesh.py`:

```diff
--- a/tests/test_reentry_model_trajectory.py
+++ b/tests/test_reentry_model_trajectory.py
@@ -169,3 +169,26 @@
     assert header[:len(tj.CSV_COLUMNS)] == tj.CSV_COLUMNS and header[-1] == "convective_heat_W"
     back = tj.read_history_csv(tmp_path / "h.csv")
     assert back.columns["convective_heat_W"].tolist() == [1.0, 2.0, 3.0] and len(back) == 3
+
+
+def test_a_step_that_would_end_within_the_snap_of_t_max_ends_at_it():
+    """Sub-plan 22 Task 4 (fact 179): t + dt added again and again falls short of t_max by round-off -- ten steps of
+    0.1 s reach 0.9999999999999999 s -- and the stepper then takes an eleventh step of 1.1e-16 s. With `snap` a step
+    that would end within it short of t_max ends at t_max: ten steps, the last ending at 1.0 s exactly, with the state of
+    the first nine unchanged to the last bit; with `snap` 0, the default, every step is as it was."""
+    atm = atmosphere.US76TableAtmosphere()
+    runs = {}
+    for snap in (None, 0.0, 1e-7):
+        sim = make(R100, atm, t_max=1.0)
+        steps, states = [], []
+        while sim.end_reason is None:
+            steps.append(sim.advance(0.1) if snap is None else sim.advance(0.1, snap=snap))
+            states.append(sim.y.copy())
+        runs[snap] = (steps, np.array(states), sim.t)
+    steps, states, t = runs[None]
+    assert len(steps) == 11 and 0.0 < steps[-1] < 1e-15 and t == 1.0
+    assert runs[0.0][0] == steps and np.array_equal(runs[0.0][1], states)
+    snapped, snapped_states, t_snap = runs[1e-7]
+    assert len(snapped) == 10 and t_snap == 1.0 and snapped[:9] == steps[:9]
+    assert np.array_equal(snapped_states[:9], states[:9]) and snapped[-1] == pytest.approx(0.1, rel=1e-14)
+    assert np.abs(snapped_states[-1] - states[-1]).max() < 1e-6                   # the eleventh step moved it by nothing
```

```diff
--- a/tests/test_reentry_model_coupled.py
+++ b/tests/test_reentry_model_coupled.py
@@ -285,3 +285,73 @@
     assert steps(True, [(True, continuum)]).switch_trigger == "both"
     run = steps(True, states, dt_continuum=None)
     assert run.switched_at is None and run.settings.dt == 0.5
+
+
+# -- the outputs of a spraying skin flight (sub-plan 22 Task 4) ----------------------------------------------------------
+def test_the_skins_loop_ends_at_t_max_without_a_vanishing_step_and_the_element_model_s_is_unchanged(coarse_sphere_mesh):
+    """Fact 179's vanishing last macro step: the steps' sum falls short of t_max by round-off -- ten 0.1 s steps reach
+    0.9999999999999999 s -- so the loop took an eleventh step of 1.1e-16 s, a whole macro step of the body over which
+    the skins' droplet-momentum ratio read 5.5 and 8.6. With `snap_to_t_max` (the layered skins') the tenth step ends at
+    t_max and there is no eleventh; without it (the element model's loop, the default) the run is as it was, row for
+    row and bit for bit. The first nine steps are the same either way (pyamg's draws seeded alike)."""
+    out = {}
+    for snap in (False, True):
+        np.random.seed(12345)
+        b = thermal_body(coarse_sphere_mesh)
+        sim = simulator(b, t_max=1.0)
+        out[snap] = coupled.CoupledRun(sim, b, heating.SesamEquivalentHeating(),
+                                       coupled.CoupledSettings(dt=0.1, snap_to_t_max=snap)).run()
+    old, new = out[False].columns, out[True].columns
+    assert coupled.CoupledSettings().snap_to_t_max is False and coupled.END_SNAP == 1e-6
+    assert len(old["time_s"]) == 12 and 0.0 < old["time_s"][-1] - old["time_s"][-2] < 1e-15 and old["time_s"][-1] == 1.0
+    assert len(new["time_s"]) == 11 and new["time_s"][-1] == 1.0 and out[True].end_reason == "t_max"
+    for k in old:
+        np.testing.assert_array_equal(new[k][:10], old[k][:10], err_msg=k)
+    assert new["temperature_K"][-1] == pytest.approx(old["temperature_K"][-1], rel=1e-12)
+
+
+def test_the_writer_writes_a_skin_flight_s_summary_from_its_droplet_columns_and_a_spool_as_it_writes_rows(tmp_path):
+    """`coupled.write_particles` (sub-plan 22 Task 4): a skin flight's history has the element model's per-step droplet
+    columns, so the summary is written -- its released mass the step's (the skin history's released_step_mass_kg; its
+    released_mass_kg is since time zero) and, where the element model has its film, the skins' liquid and flowing depths;
+    the element model's history writes the element model's summary as before. A skin body's rows written from a spool on
+    disk give particles.npz, the summary and the size distribution exactly as from blocks in memory."""
+    import types
+    from reentry_model import spray
+    rng = np.random.default_rng(7)
+    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
+
+    def block(t, k):
+        rows = rng.uniform(0.1, 1.0, size=(k, len(spray.SOURCE_COLUMNS)))
+        rows[:, C["time_s"]], rows[:, C["r_m"]] = t, rng.uniform(5e-5, 2e-3, size=k)
+        return rows
+    blocks = [block(0.25, 40), block(0.5, 70)]
+    t = np.array([0.0, 0.25, 0.5])
+    step = np.array([0.0, blocks[0][:, C["dm_kg"]].sum(), blocks[1][:, C["dm_kg"]].sum()])
+    skin = {"time_s": t, "altitude_km": [65.0, 64.9, 64.8], "velocity_kms": [7.24, 7.23, 7.22], "released_mass_kg": np.cumsum(step),
+            "released_step_mass_kg": step, "n_released": [0.0, 1e3, 3e3], "r_median_um": [np.nan, 300.0, 310.0],
+            "r_max_um": [np.nan, 1900.0, 1950.0], "theta_cr_deg": [np.nan, 0.5, 0.4], "spraying_area_m2": [0.0, 1e-3, 2e-3],
+            "liquid_depth_mean_mm": [0.0, 0.01, 0.02], "flowing_depth_mean_mm": [0.0, 0.1, 0.2]}
+    history = tj.History({k: np.asarray(v, dtype=float) for k, v in skin.items()}, np.zeros((3, 6)), "t_max")
+    spool = spray.SourceSpool(str(tmp_path / "spool"))
+    for blk in blocks:
+        spool.append(blk)
+    written = {}
+    for name, rows in (("memory", blocks), ("disk", spool)):
+        written[name] = coupled.write_particles(str(tmp_path / name), types.SimpleNamespace(source_rows=rows), history)
+    a, b = written["memory"], written["disk"]
+    assert set(a) == set(b) == {"particles", "particles_summary", "size_distribution"}
+    pa, pb = np.load(a["particles"]), np.load(b["particles"])
+    assert pa.files == pb.files == spray.SOURCE_COLUMNS
+    assert all(np.array_equal(pa[k], pb[k], equal_nan=True) for k in pa.files)
+    for key in ("particles_summary", "size_distribution"):
+        assert open(a[key]).read() == open(b[key]).read()
+    lines = open(a["particles_summary"]).read().splitlines()
+    assert lines[0].split(",") == [h for h, _ in coupled.SKIN_SUMMARY_COLUMNS] and lines[0].split(",")[9:] == [
+        "liquid_depth_mean_mm", "flowing_depth_mean_mm"]
+    released = [float(line.split(",")[3]) for line in lines[1:]]
+    np.testing.assert_allclose(released, step, rtol=1e-8)                           # the step's, not the running total
+    element = dict(skin, film_mass_kg=[0.0, 1e-3, 2e-3], film_thickness_mean_mm=[0.0, 0.1, 0.2])
+    hist_e = tj.History({k: np.asarray(v, dtype=float) for k, v in element.items()}, np.zeros((3, 6)), "t_max")
+    files = coupled.write_particles(str(tmp_path / "element"), types.SimpleNamespace(source_rows=blocks), hist_e)
+    assert open(files["particles_summary"]).readline().strip().split(",") == coupled.SUMMARY_COLUMNS
```

```diff
--- a/tests/test_reentry_model_movemesh.py
+++ b/tests/test_reentry_model_movemesh.py
@@ -624,3 +624,23 @@
     elem, bary = movemesh.locate(old_points, old_tets, qp)
     assert len(qp) > 30000
     assert bary.min(axis=1).min() > -1e-9
+
+
+def test_a_rebuild_reports_the_mass_its_new_volume_misses_and_sums_it(coarse_sphere_mesh):
+    """Spec 10.2 (sub-plan 22 Task 4): every rebuild reports its mass residual beside its energy residual -- the
+    material's density times the old volume less the new, uncorrected -- and the moving mesh sums them. The surface is
+    kept to the last bit, so it is round-off: under 1e-12 of the body's mass on each of two rebuilds."""
+    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
+    mat = constant_material()
+    solver = solver_on(m, mat)
+    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3, sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
+    assert moving.rebuild_mass_residual_total == 0.0 and movemesh.RebuildResult(0.0, 0.0, 1, 1, 0.0, "x").mass_residual == 0.0
+    residuals = []
+    for _ in range(2):
+        assert not moving.move(inward(moving, 5e-4)).refused
+        V = moving.mesh.volume()
+        result = moving.rebuild()
+        assert result.mass_residual == mat.rho * (V - moving.mesh.volume())
+        assert abs(result.mass_residual) < 1e-12 * mat.rho * V
+        residuals.append(result.mass_residual)
+    assert moving.rebuild_mass_residual_total == residuals[0] + residuals[1]
```

- [ ] **Step 3: The skin body's outputs** -- in `tests/test_reentry_model_skinbody.py` (and sub-plan 21 Task 6c's
  flow-off test comparing the two bodies' history columns, now nan-aware: the new columns are empty in both):

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -661,7 +661,8 @@
     assert make(coarse_sphere_mesh, T0=600.0).flow_on is False
     assert same_body(off, on) and np.all(on.skins.T[on.skins.valid] < on.material.T_rigid)
     a, c = off.skin_stats(), on.skin_stats()
-    assert all(a[k] == c[k] for k in a if not k.startswith("wall_"))
+    same = lambda x, y: x == y or (np.isnan(x) and np.isnan(y))     # a column empty in both (sub-plan 22 Task 4's)
+    assert all(same(a[k], c[k]) for k in a if not k.startswith("wall_"))
     assert c["speed_max_ms"] == 0.0 and c["momentum_residual"] == 0.0 and c["flow_courant_max"] == 0.0
 
 
@@ -1162,7 +1163,7 @@
     rows = plan.source_block(0.0125)
     C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
     i, mode = np.array([0, 2, 1, 2]), np.array([0, 0, 1, 1])
-    assert rows.shape == (4, 26)
+    assert rows.shape == (4, 27) and np.isnan(rows[:, C["release_angle_register_deg"]]).all()    # no centroid given
     np.testing.assert_array_equal(rows[:, C["x_m"]:C["z_m"] + 1], pts[i])
     np.testing.assert_array_equal(rows[:, C["theta_deg"]], np.degrees(theta[i]))
     assert (rows[:, C["time_s"]] == 0.0125).all() and (rows[:, C["altitude_km"]] == 71.0).all()
@@ -1208,16 +1209,25 @@
     and the books' released mass, their enthalpy the books' released enthalpy, and their momentum, mass times speed
     along their direction, the momentum the droplets took from each column, all to round-off; the size histograms
     accumulated each flow sub-step are the rows'. The coupled writer writes the rows and their size distribution and
-    reads them back to the last bit; the summary, whose columns are the element model's history, waits for Task 4."""
+    reads them back to the last bit.
+
+    Sub-plan 22 Task 4: the rows go to a spool on disk as they come, which gives back the block it was given to the
+    last bit; each carries Register et al.'s release angle of its column's top corner from the centroid of the volume
+    the top surface enclosed where the plan was evaluated; the history reports the step's released mass, the droplets
+    since time zero, the radius median and largest over the step's rows, the smallest angle and the area of the
+    columns whose plan sprays, the release angle's mean by mass -- the mean recomputed from the step's rows -- least and
+    largest, and b and the resolved sheared depth over delta_m on the wet windward columns; and from those columns the
+    coupled writer writes the summary, whose released mass is the step's rows'."""
     pytest.importorskip("cantera")
     from test_reentry_model_surface_flow import aero_state_at
     from reentry_model import coupled, spray
     from reentry_model import trajectory as tj
     mat = material.Material.from_drama_json("AA7075_scheil")
+    np.random.seed(12345)              # pyamg's draws, as a run seeds them: the histograms' round-off is the path's
     b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid),
-             slurry_spray=True)
-    plans, books = [], []
-    girin, book = skinbody.girin_plan, skinbody.SprayPlan.book
+             slurry_spray=True, spool=str(tmp_path / "spool"))
+    plans, books, spooled = [], [], []
+    girin, book, append = skinbody.girin_plan, skinbody.SprayPlan.book, spray.SourceSpool.append
 
     def planned(*a, **k):
         plans.append(girin(*a, **k))
@@ -1226,13 +1236,20 @@
     def booked(plan, fb):
         books.append((plan, fb.released.copy(), fb.released_slurry.copy()))
         return book(plan, fb)
+    def kept(spool, block):
+        spooled.append(block.copy())
+        return append(spool, block)
     monkeypatch.setattr(skinbody, "girin_plan", planned)
     monkeypatch.setattr(skinbody.SprayPlan, "book", booked)
+    monkeypatch.setattr(spray.SourceSpool, "append", kept)
     state = aero_state_at(93.0)
+    stats = [b.skin_stats()]
     for k in range(2):
         b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
-    assert len(b.source_rows) == 1 and b.particles                      # the first step melts the tops
+        stats.append(b.skin_stats())
+    assert len(b.source_rows) == 1 and b.particles and isinstance(b.source_rows, spray.SourceSpool)   # the first step melts
     rows = spray.source_table(b.source_rows)
+    assert len(spooled) == 1 and np.array_equal(rows, spooled[0], equal_nan=True)     # the spool's rows: the block's
     C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
     p = plans[-1]
     index = {tuple(x): i for i, x in enumerate(p.points)}
@@ -1263,15 +1280,54 @@
     assert (h[mode == 0] > h_liq * (1.0 - 1e-6)).all() and (h_rigid < h[mode == 1]).all() and (h[mode == 1] < h_liq).all()
     # the size histograms, accumulated each flow sub-step, are the rows': in the same bins, to the round-off of numpy's
     # weighted histogram, which takes each bin as a difference of running sums over all the droplets (measured: 1.6e-12
-    # droplets in a bin of 0.065, 2e-17 of the 71 000 in all)
+    # droplets in a bin of 0.065, 2e-17 of the 71 000 in all; with the draws seeded 5.3e-16 of the 89 173 in all, and
+    # unseeded, after other tests, up to 1.04e-15, sub-plan 22 Task 4)
     n_rec, m_rec = spray.histogram(rows[:, C["r_m"]], rows[:, C["dn"]], rows[:, C["dm_kg"]])
     np.testing.assert_array_equal(b.hist_n > 0.0, n_rec > 0.0)
     np.testing.assert_allclose(b.hist_n, n_rec, rtol=0.0, atol=1e-15 * n_rec.sum())
     np.testing.assert_allclose(b.hist_m, m_rec, rtol=0.0, atol=1e-15 * m_rec.sum())
-    # the coupled writer: the same files as the element model's, read back to the last bit; no summary until Task 4
-    history = tj.History({"time_s": np.array([0.0, 0.25, 0.5])}, np.zeros((3, 6)), None)
+    # Task 4: the release angle of each row's column top from the plan's centroid, the top surface's where it was
+    # evaluated: 111.5-169.9 degrees, the columns within 10.9 degrees of the nose not spraying, none behind the equator,
+    # and within 7.3 degrees of 180 less the corner normal's angle theta_deg (measured)
+    angle = rows[:, C["release_angle_register_deg"]]
+    np.testing.assert_array_equal(p.centre, spray.enclosed_centroid(p.points, b.surface.faces))
+    np.testing.assert_array_equal(angle, spray.register_angle(p.points, p.centre, b.v_hat)[i])
+    assert 160.0 < angle.max() < 180.0 and angle.min() > 90.0
+    assert (np.abs(angle - (180.0 - rows[:, C["theta_deg"]])) < 10.0).all()
+    # the history's droplet columns: nothing in the first step, the second's from its rows and its plan
+    s0, s1, s2 = stats
+    assert s1["released_step_mass_kg"] == 0.0 and np.isnan(s1["r_median_um"]) and np.isnan(s1["release_angle_register_mean_deg"])
+    assert np.isnan(s0["b_over_delta_m_median"]) and s0["n_released"] == 0.0 and s0["spraying_area_m2"] == 0.0
+    dm = rows[:, C["dm_kg"]]
+    assert s2["released_step_mass_kg"] == pytest.approx(dm.sum(), rel=1e-12) and s2["n_released"] == b.droplet_number
+    assert s2["r_median_um"] == np.median(rows[:, C["r_m"]]) * 1e6 and s2["r_max_um"] == rows[:, C["r_m"]].max() * 1e6
+    on = p.spraying()
+    assert s2["theta_cr_deg"] == np.degrees(p.theta[on].min()) and s2["spraying_area_m2"] == p.area[on].sum() > 0.0
+    assert s2["release_angle_register_mean_deg"] == pytest.approx((dm * angle).sum() / dm.sum(), rel=1e-14)
+    assert (s2["release_angle_register_min_deg"], s2["release_angle_register_max_deg"]) == (angle.min(), angle.max())
+    assert b.sprayed_mass == pytest.approx(b.released_mass, rel=1e-12) and b.n_source_rows == len(rows)
+    # b and the resolved sheared depth over delta_m on the wet windward columns, the step's end
+    F = b.skins.depth_above(mat.T_rigid)
+    delta_m, _ = spray.melt_layer(b.last_flow, mat.liquid)
+    wet = (b.theta <= 0.5 * np.pi) & (F > 0.0) & np.isfinite(delta_m) & (delta_m > 0.0)
+    assert wet.sum() > 100 and s2["b_over_delta_m_median"] == np.median(F[wet] / delta_m[wet])
+    sheared = skinbody.sheared_depth(b.skins)
+    moving = wet & np.isfinite(sheared)
+    assert s2["sheared_depth_over_delta_m_median"] == np.median(sheared[moving] / delta_m[moving]) > 0.0
+    # the coupled writer: the same files as the element model's, read back to the last bit, and the summary from the
+    # history's droplet columns
+    keys = {k for _, k in coupled.SKIN_SUMMARY_COLUMNS} - {"time_s", "altitude_km", "velocity_kms"}
+    columns = {k: np.array([x[k] for x in stats], dtype=float) for k in keys}
+    columns.update({"time_s": np.array([0.0, 0.25, 0.5]), "altitude_km": np.full(3, state.h / 1e3),
+                    "velocity_kms": np.full(3, state.V / 1e3)})
+    history = tj.History(columns, np.zeros((3, 6)), None)
     files = coupled.write_particles(str(tmp_path), b, history)
-    assert set(files) == {"particles", "size_distribution"}
+    assert set(files) == {"particles", "particles_summary", "size_distribution"}
+    with open(files["particles_summary"]) as fh:
+        summary = [line.split(",") for line in fh.read().splitlines()]
+    assert summary[0] == [h for h, _ in coupled.SKIN_SUMMARY_COLUMNS] and len(summary) == 4
+    assert float(summary[3][3]) == pytest.approx(dm.sum(), rel=1e-8) and float(summary[2][3]) == 0.0
+    assert float(summary[3][4]) == pytest.approx(rows[:, C["dn"]].sum(), rel=1e-8)
     npz = np.load(files["particles"])
     assert npz.files == spray.SOURCE_COLUMNS
     assert all(np.array_equal(npz[name], rows[:, k], equal_nan=True) for k, name in enumerate(spray.SOURCE_COLUMNS))
@@ -1296,3 +1352,28 @@
     cold = skinbody.SkinBody(copy, material.Material.from_drama_json("AA7075_nomelt"), thermal_solver("skfem"), 1.47, DEV,
                              T0=1200.0, v_hat=(-1.0, 0.0, 0.0))
     assert not cold.liquid_present()                             # a material that does not melt holds no liquid
+
+
+# -- the outputs of a spraying skin flight (sub-plan 22 Task 4) ----------------------------------------------------------
+def test_the_resolved_sheared_depth_is_where_the_speed_falls_to_a_tenth_of_the_top_s():
+    """Spec 6.7's resolved sheared depth (`skinbody.sheared_depth`, SHEARED_FRACTION 0.1): from the top to where the
+    computed speed falls to a tenth of the top cell's, between the cells' centres. Columns of 10 um cells: a
+    profile falling linearly to rest 0.2 mm down, so 0.9 of that depth and a tenth of the top cell centre's, 5.9 um down
+    (on Girin's linear profile over delta_m, 0.9 delta_m); a skin moving as one piece, its whole thickness; a top at
+    rest, nan."""
+    from test_reentry_model_skin import LIQUID, RIGID, stacked, strip_grid
+    col, _ = strip_grid(3, 2, h=20e-3)
+    sk = stacked(col, [[LIQUID] * 30 + [RIGID] * 10, [LIQUID] * 40] + [[LIQUID] * 40] * 4)
+    dz = sk.dz()
+    z = np.cumsum(dz, axis=1) - 0.5 * dz
+    D = 2e-4
+    U = np.zeros(sk.U.shape)
+    U[0, :, 0] = 3.0 * np.clip(1.0 - z[0] / D, 0.0, None)
+    U[1, :, 1] = 0.5
+    sk.U = U * sk.valid[..., None]
+    assert skinbody.SHEARED_FRACTION == 0.1
+    got = skinbody.sheared_depth(sk)
+    assert got[0] == pytest.approx(0.9 * D + 0.1 * z[0, 0], rel=1e-12) and 5e-6 < z[0, 0] < 7e-6
+    assert got[1] == pytest.approx(sk.thickness()[1], rel=1e-14)
+    assert np.isnan(got[2:]).all()
+    assert skinbody.sheared_depth(sk, 0.5)[0] == pytest.approx(0.5 * D + 0.5 * z[0, 0], rel=1e-12)
```

- [ ] **Step 4: The element model's and the command line's** -- in `tests/test_reentry_model_melting.py`,
  `tests/test_reentry_model_cli.py` and `tests/test_reentry_model_fenicsx.py`:

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -95,7 +95,14 @@
     assert b.material.h_liquid < h_out <= b.material.enthalpy(b.surface_temperature().max())       # own temperature, superheat and all
     assert stats["film_T_max_K"] == pytest.approx(b.surface_temperature().max()) and stats["film_T_max_K"] > b.material.T_liquidus
     rows = np.array(b.source_rows)
-    assert rows.shape[1] == 26 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
+    assert rows.shape[1] == 27 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
+    # sub-plan 22 Task 4: Register et al.'s release angle of every row, from the centroid of the volume the surface
+    # enclosed at the spray step, windward; the history's mean by mass, least and largest are the last step's rows'
+    from reentry_model import spray as spray_mod
+    angle, last = rows[:, 26], rows[:, 0] == rows[-1, 0]
+    assert np.isfinite(angle).all() and (angle > 45.0).all() and (angle <= 180.0).all()
+    assert tuple(stats[k] for k in spray_mod.RELEASE_ANGLE_COLUMNS) == spray_mod.release_angle_stats(rows[last, 14],
+                                                                                                     angle[last])
     # sub-plan 22 Task 2: the film is liquid and leaves at its own enthalpy, the removed enthalpy's
     assert (rows[:, 24] == 0.0).all() and np.isfinite(rows[:, 25]).all()
     assert (rows[:, 14] * rows[:, 25]).sum() == pytest.approx(b.removed_enthalpy, rel=1e-12)
```

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -295,6 +295,11 @@
     assert r["melt_onset_altitude_km"] is not None and r["sprayed_mass_kg"] > 0.0 and r["n_source_rows"] > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
     for key in ("particles", "particles_summary", "size_distribution"):
         assert os.path.isfile(f[key])
+    # sub-plan 22 Task 4: the release angle in the history and the rows; the summary the element model's as before
+    from reentry_model import spray
+    assert all(k in rows[0] for k in spray.RELEASE_ANGLE_COLUMNS)
+    assert "release_angle_register_deg" in np.load(f["particles"]).files
+    assert open(f["particles_summary"]).readline().strip().split(",") == coupled.SUMMARY_COLUMNS
     assert len(f["melt_plots"]) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in f["melt_plots"])
     assert doc["comparison"] is None and f["film"] is None
     names = [os.path.basename(p) for p in f["stills"]]
@@ -459,8 +464,10 @@
     suffixes, with no surface or flow suffix (both defaults); the settings record the flow and the spray model's
     constants, the results the switch's trigger; mass,
     energy and momentum close; the droplets are written by the element model's writer -- particles.npz, whose rows carry
-    the mass released, and size_distribution.csv, their histogram -- and the per-step summary, which waits for Task 4's
-    history columns, is recorded as not written; the history's flow columns and the frames' top speed are written."""
+    the mass released, and size_distribution.csv, their histogram -- and, since Task 4, particles_summary.csv from the
+    skin history's droplet columns, each step's released mass its rows', the rows having gone through a spool on disk
+    that the run then removes; every row and step carries Register's release angle; the results count the sprayed
+    mass, the droplets and the rows; the history's flow columns and the frames' top speed are written."""
     pytest.importorskip("cantera")
     import glob
     import pyvista as pv
@@ -486,12 +493,25 @@
     # the books: mass, energy and momentum
     assert r["initial_mass_kg"] - r["final_mass_kg"] == pytest.approx(released, rel=1e-9)
     assert abs(r["energy_balance_residual"]) < 1e-8 and r["momentum_residual"] < 1e-12
-    # the droplet files, the element model's, and no summary until Task 4
-    assert f["particles_summary"] is None and not os.path.exists(os.path.join(os.path.dirname(f["particles"]),
-                                                                                  "particles_summary.csv"))
+    # the droplet files, the element model's, and the summary from the skin history (Task 4); the spool removed
     table = np.load(f["particles"])
     assert table.files == spray.SOURCE_COLUMNS and len(table["dm_kg"]) > 0 and (table["time_s"] == 0.0625).all()
     assert table["dm_kg"].sum() == pytest.approx(released, rel=1e-12)
+    assert not os.path.exists(os.path.join(os.path.dirname(f["particles"]), cli.SOURCE_SPOOL))
+    assert r["sprayed_mass_kg"] == pytest.approx(released, rel=1e-12) and r["n_source_rows"] == len(table["dm_kg"])
+    assert r["n_released"] == pytest.approx(table["dn"].sum(), rel=1e-12) and r["wall_flow_s"] > 0.0
+    with open(f["particles_summary"]) as fh:
+        summary = list(csv.DictReader(fh))
+    assert list(summary[0]) == [h for h, _ in coupled.SKIN_SUMMARY_COLUMNS] and len(summary) == len(rows)
+    assert [float(x["released_mass_kg"]) for x in summary[:-1]] == [0.0] * 4
+    assert float(summary[-1]["released_mass_kg"]) == pytest.approx(released, rel=1e-8)
+    assert float(summary[-1]["n_released"]) == pytest.approx(r["n_released"], rel=1e-8)
+    angle = table["release_angle_register_deg"]
+    assert (angle > 90.0).all() and (angle <= 180.0).all() and np.isnan(float(rows[-2]["release_angle_register_mean_deg"]))
+    mean = (table["dm_kg"] * angle).sum() / table["dm_kg"].sum()
+    assert float(rows[-1]["release_angle_register_mean_deg"]) == pytest.approx(mean, rel=1e-8)
+    assert float(rows[-1]["b_over_delta_m_median"]) > 0.0 and float(rows[-1]["released_step_mass_kg"]) == pytest.approx(
+        released, rel=1e-8)
     dist = [row for row in csv.DictReader(open(f["size_distribution"])) if float(row["window_end_s"]) == 0.0625]
     n_hist, m_hist = spray.histogram(table["r_m"], table["dn"], table["dm_kg"])
     np.testing.assert_allclose([float(row["dn"]) for row in dist], n_hist, rtol=1e-8)
```

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -197,6 +197,11 @@
     (measured 2026-10-05); with it 2.7e7, set by the deepest piles' G b^3 / (3 mu) term, not by a vanishing film.
     Measured: the same 7 983 active elements, mass 5.8e-10, sprayed mass 5.9e-9, runoff 7.5e-10, film 1.7e-8, deep
     account 2.8e-10, temperatures 2.3e-4 K, balances 5.1e-9 and 1.1e-11.
+    Both solvers' conjugate gradients run to 1e-13 here (sub-plan 22 Task 4). At their default 1e-10 the scikit-fem
+    backend's balance is its linear solves' tolerance, not the books': the books are the one MeltingBody's, closing to
+    8.9e-12 on the FEniCSx backend, and on scikit-fem's they wandered 2.2e-9 to 1.17e-8 over the ten steps (5.1e-9 at
+    the last on sub-plan 21's code, 1.05e-8 with the spraying amendments of 2026-10-22, which spray 2.2 % more), but
+    with its gradients at 1e-13 and its Newton iterations unchanged 4.5e-12 (measured 2026-10-31).
     The rigid substrate (amendment of 2026-10-06) is off here: with the body's interior at 850 K, above AA7075_range's
     829 K half-liquid point, the whole body is slurry, every wet windward patch takes the thick branch, and films of
     about a micron that straddle the spray floor B_MIN are released whole by one backend and not by the other from the
@@ -226,8 +231,8 @@
     for name in ("skfem", "fenicsx"):
         np.random.seed(12345)
         m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
-        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM,
-                             settings=body.MeltSettings(rigid_substrate=False))
+        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name, cg_tol=1e-13),
+                             MASS_100MM, settings=body.MeltSettings(rigid_substrate=False))
         sim = simulator(b, t_max=90.0)
         sim.advance(50.0)
         p = b.mesh.points
@@ -242,6 +247,7 @@
                       state=a)
         out[name] = b
     s, f = out["skfem"], out["fenicsx"]
+    assert s.solver.cg_tol == f.solver.cg_tol == 1e-13
     assert seen["thin_on_thick"] > 0 and seen["c_max"] < 1e9 and f.runoff_mass > 0.0 and f.sprayed_mass > 0.0
     assert np.array_equal(s.mesh.active, f.mesh.active)
     assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6) and f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5)
```

- [ ] **Step 5: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py tests/test_reentry_model_trajectory.py
  tests/test_reentry_model_coupled.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_skinbody.py
  tests/test_reentry_model_melting.py tests/test_reentry_model_cli.py -q -k "solid_fraction_and_enthalpy or register_s_release
  or enclosed_centroid or rides_in_the_source or spool_on_disk or within_the_snap or vanishing_step or
  summary_from_its_droplet or mass_its_new_volume or sheared_depth_is or plan_s_records or
  close_on_its_books_and_the_coupled_writer or film_spraying_death_and_balances or sprays_from_its_top_cells or
  melting_run_writes_columns"`. Expected: 15 failed (about a minute) --
  - the source table's 26 columns where 27 are asked: `assert (26 == 27)`, `assert ((4, 26) == (4, 27)`, `IndexError:
    index 26 is out of bounds for axis 1 with size 26`, `assert ('h_kg' == 'release_angle_register_deg'`;
  - `AttributeError: module 'reentry_model.spray' has no attribute 'register_angle'` (and `'enclosed_centroid'`,
    `'SourceSpool'`, `'RELEASE_ANGLE_COLUMNS'`);
  - `TypeError: Simulator.advance() got an unexpected keyword argument 'snap'` and `TypeError: CoupledSettings.__init__()
    got an unexpected keyword argument 'snap_to_t_max'`;
  - `AttributeError: 'MovingMesh' object has no attribute 'rebuild_mass_residual_total'`;
  - `TypeError: SkinBody.__init__() got an unexpected keyword argument 'spool'` and `AttributeError: module
    'reentry_model.skinbody' has no attribute 'SHEARED_FRACTION'`;
  - `AttributeError: module 'reentry_model.cli' has no attribute 'SOURCE_SPOOL'`.

  The FEniCSx test (`FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q -k thick_film_flux`)
  fails on Task 3b's code as it stands, 1.0496e-8 against 1e-8, and passes with this step's version: its change is the
  test's own solvers.

- [ ] **Step 6: The release angle, the centroid and the spool** in `reentry_model/spray.py`:

```diff
--- a/reentry_model/spray.py
+++ b/reentry_model/spray.py
@@ -97,8 +97,21 @@
 Release per patch and step: dm = min(mdot A dt, m_f) (the film after feed and runoff), dn = dm / (4/3 pi rho_l r^3);
 one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX]. The layered
 skins' records are one row per column, per mode, per macro step (`skinbody.SprayPlan.source_block`, sub-plan 22 Task 2):
-the same columns, written by the same `source_block`, the two last columns the droplets' solid fraction and the enthalpy
-per kilogram they carry (spec 8 and 10.2; the 2026-10-07 design 7.5).
+the same columns, written by the same `source_block`, the droplets' solid fraction and the enthalpy per kilogram they
+carry (spec 8 and 10.2; the 2026-10-07 design 7.5) and, last, their release angle as Register et al. (2020) define it
+(`register_angle`, sub-plan 22 Task 4). A skin flight at the continuum step makes tens of millions of records, so the
+command line has the skin body write each macro step's block to disk as it comes (`SourceSpool`), and the coupled writer
+reads either form column by column (`source_columns`), the numbers the same to the last bit.
+
+The release angle (Asha, 2026-10-29; Register, P. J., et al. 2020, "Interactions between asteroid fragments during
+atmospheric entry", Icarus 337, 113468, section 4.2): the angle between the larger body's trajectory axis and the line
+joining its centroid to the released fragment's, measured from downstream, theta_R = 180 degrees - angle(v_hat, p - c):
+180 at the stagnation point, 90 at the equator, 0 directly behind. It is computed from the actual vector p - c, so a
+flattened or receded body gets its own angle (the source table's `theta_deg` is the surface normal's angle to the flight
+direction, a different thing off a sphere). p is the row's position; c is the centroid of the volume the body's outer
+surface encloses (`enclosed_centroid`), from the same surface and at the same instant as the positions the rows carry --
+the element model's surface at its spray step, the skins' top surface where their plan was evaluated, the macro step's
+start -- and v_hat the body's direction of motion through the air, its fixed attitude.
 
 The droplet's size caps (2026-10-24, Asha). Every branch's radius is capped by a quarter of the body radius and by one
 of two volumes. Without `cap_depth` (the element model) it is the film on the patch, r <= (3 m_f / (4 pi rho_l))^1/3,
@@ -113,6 +126,7 @@
 about 2 % (thick: k_r = 0.17) and 6.5 % (thin: lambda*/4) of the wavelength; on the Rayleigh-Taylor branch, whose
 droplet is its wavelength (centimetres), on any layer shallower than about four wavelengths, so always: at lambda = 2 cm
 and h = 1 mm it is 4.6 mm, against about 2 mm under the film cap."""
+import os
 from dataclasses import dataclass, field
 
 import numpy as np
@@ -129,10 +143,12 @@
 # h_floor_m and gamma_tau: viscous_film's floor depth and gamma tau* (nan off Girin & Kopyt's branches; 2026-10-22);
 # solid_fraction and h_kg: the droplets' solid share of their mass and the enthalpy per kilogram they carry [J/kg, on
 # the material's enthalpy scale] (spec 8 and 10.2, sub-plan 22 Task 2): the element model's film is liquid by
-# construction and booked at its own liquid enthalpy, so its rows carry 0 and that enthalpy (`MeltingBody`)
+# construction and booked at its own liquid enthalpy, so its rows carry 0 and that enthalpy (`MeltingBody`);
+# release_angle_register_deg: Register et al.'s (2020) release angle, 180 degrees at the stagnation point, 90 at the
+# equator, 0 directly behind (`register_angle`; Asha 2026-10-29, sub-plan 22 Task 4)
 SOURCE_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "theta_deg", "x_m", "y_m", "z_m", "closure", "branch", "b_m",
                   "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup",
-                  "h_floor_m", "gamma_tau", "solid_fraction", "h_kg"]
+                  "h_floor_m", "gamma_tau", "solid_fraction", "h_kg", "release_angle_register_deg"]
 BRANCH_THICK, BRANCH_THIN, BRANCH_RAREFIED, BRANCH_RT = 0, 1, 2, 3      # Girin thick / thin-film with the edge state /
 #                                                                        thin-film with the freestream / front-surface Rayleigh-Taylor
 BRANCH_REGIME2 = 4                 # deep film whose melt is the more viscous medium: Girin & Kopyt's step profile again
@@ -499,13 +515,15 @@
                            h_floor, gamma_tau, outside & gk, wavelength)
 
 
-def source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None):
+def source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None,
+                 release_angle=None):
     """Source-table rows, an array (k, len(SOURCE_COLUMNS)), for the patches that released mass (`res.dm` > 0): `h` [m]
     and `V` [m/s] the altitude and speed, `theta` [rad], `centroids` (n, 3) [m] and `t_hat` (n, 3) per patch, `flow`'s
     `closure`, `res` a SprayResult (or anything with its per-patch `branch`, `delta_m`, `we_s`, `r`, `dn`, `dm`, `v_s`,
     `h_floor`, `gamma_tau`), `b` [m] the film thickness, `liquid` the droplets' properties -- each of its `rho`, `mu` and
     `sigma` one value or one per patch (a slurry droplet's are its slurry's: sub-plan 22 Task 2) -- for their Weber and
-    Ohnesorge numbers at the freestream `state`, and `solid_fraction` and `h_kg` [J/kg] per patch (0 and nan if not given)."""
+    Ohnesorge numbers at the freestream `state`, `solid_fraction` and `h_kg` [J/kg] per patch (0 and nan if not given),
+    and `release_angle` [deg] per patch, Register et al.'s (`register_angle`; nan if not given; sub-plan 22 Task 4)."""
     k = np.flatnonzero(res.dm > 0.0)
     if k.size == 0:
         return np.zeros((0, len(SOURCE_COLUMNS)))
@@ -516,28 +534,141 @@
     oh = mu / np.sqrt(rho * sigma * 2.0 * r)
     s_f = np.zeros(k.size) if solid_fraction is None else np.asarray(solid_fraction, dtype=float)[k]
     h_d = np.full(k.size, np.nan) if h_kg is None else np.asarray(h_kg, dtype=float)[k]
+    angle = np.full(k.size, np.nan) if release_angle is None else np.asarray(release_angle, dtype=float)[k]
     return np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
                             centroids[k, 0], centroids[k, 1], centroids[k, 2], flow.closure[k], res.branch[k], b[k],
                             res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1],
                             t_hat[k, 2], we_d, oh, (we_d > WE_BREAKUP).astype(float), res.h_floor[k], res.gamma_tau[k],
-                            s_f, h_d])
-
-
-def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None):
+                            s_f, h_d, angle])
+
+
+def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None,
+                release_angle=None):
     """Source-table rows (list of lists in SOURCE_COLUMNS order) for the patches that released mass this step
     (`source_block`'s, as lists: the element model's form, one row per patch-step)."""
-    return source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction, h_kg).tolist()
+    return source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction, h_kg,
+                        release_angle).tolist()
 
 
 def source_table(rows):
     """A body's `source_rows` as one array (N, len(SOURCE_COLUMNS)): a list of rows (the element model's, one per
     patch-step) or of blocks of rows (the layered skins', one array per macro step: a Python list of floats costs about
-    0.9 kB a row against 0.2 kB in an array, and a flight at the continuum step makes of order ten million records)."""
+    0.9 kB a row against 0.2 kB in an array, and a flight at the continuum step makes of order ten million records), or a
+    `SourceSpool` on disk, read into memory."""
+    if isinstance(rows, SourceSpool):
+        return np.column_stack([np.asarray(c) for c in rows.columns().values()]) if rows.n_rows \
+            else np.zeros((0, len(SOURCE_COLUMNS)))
     if not len(rows):
         return np.zeros((0, len(SOURCE_COLUMNS)))
     if isinstance(rows[0], np.ndarray) and rows[0].ndim == 2:
         return np.concatenate(rows, axis=0)
     return np.array(rows, dtype=float)
+
+
+def source_columns(rows):
+    """A body's `source_rows` as {name: one-dimensional array} in SOURCE_COLUMNS order, what the coupled writer writes: a
+    list of rows or of blocks, the columns of `source_table`'s array; a `SourceSpool`, its files mapped from disk, so that
+    no flight's records need be in memory at once. The values are the same to the last bit either way."""
+    if isinstance(rows, SourceSpool):
+        return rows.columns()
+    table = source_table(rows)
+    return {name: table[:, i] for i, name in enumerate(SOURCE_COLUMNS)}
+
+
+class SourceSpool:
+    """The layered skins' droplet records on disk as they come (sub-plan 22 Task 4): one file per source column,
+    `<directory>/<name>.f64` (float64, little-endian), each macro step's block appended to every one. Measured for the 100
+    mm flight at the continuum step from the first liquid (24.5 s) to 120 s: 7 640 macro steps of 4 600-9 200 rows, 35-70
+    million rows of 27 columns, 7.6-15.2 GB as arrays and twice that while `source_table` joins them -- more than the
+    laptop's 24 GB -- so the command line gives the skin body a spool in the run's directory and the coupled writer maps
+    the columns back (`columns`), each read once, contiguous. `len` counts the blocks (macro steps that released), as the
+    in-memory list does; `n_rows` the rows."""
+
+    def __init__(self, directory):
+        os.makedirs(directory, exist_ok=True)
+        self.directory, self.n_rows, self.n_blocks = directory, 0, 0
+        for name in SOURCE_COLUMNS:
+            open(self._path(name), "wb").close()
+
+    def _path(self, name):
+        return os.path.join(self.directory, name + ".f64")
+
+    def __len__(self):
+        return self.n_blocks
+
+    def append(self, block):
+        """Add a macro step's block of rows (k, len(SOURCE_COLUMNS)) to the end of every column's file."""
+        block = np.asarray(block, dtype=np.float64)
+        if block.ndim != 2 or block.shape[1] != len(SOURCE_COLUMNS):
+            raise ValueError("a block of source rows has {} columns, not {}".format(block.shape[-1], len(SOURCE_COLUMNS)))
+        if not len(block):
+            return
+        for i, name in enumerate(SOURCE_COLUMNS):
+            with open(self._path(name), "ab") as fh:
+                fh.write(np.ascontiguousarray(block[:, i], dtype="<f8").tobytes())
+        self.n_rows += len(block)
+        self.n_blocks += 1
+
+    def columns(self):
+        """{name: the column, read-only, mapped from its file} in SOURCE_COLUMNS order (empty arrays with no rows)."""
+        if not self.n_rows:
+            return {name: np.zeros(0) for name in SOURCE_COLUMNS}
+        return {name: np.memmap(self._path(name), dtype="<f8", mode="r", shape=(self.n_rows,)) for name in SOURCE_COLUMNS}
+
+    def remove(self):
+        """Delete the files and the directory, once what they hold has been written elsewhere."""
+        for name in SOURCE_COLUMNS:
+            if os.path.exists(self._path(name)):
+                os.remove(self._path(name))
+        if os.path.isdir(self.directory) and not os.listdir(self.directory):
+            os.rmdir(self.directory)
+
+
+def register_angle(points, centre, v_hat):
+    """Register et al.'s (2020, section 4.2) release angle [deg] of each point (n, 3) [m] on a body whose centroid is
+    `centre` (3,) [m] and which moves through the air along `v_hat`: 180 degrees minus the angle between v_hat and the
+    line from the centroid to the point, so 180 at the stagnation point, 90 at the equator and 0 directly behind,
+    measured from the actual vector so that a flattened or receded body gets its own angle (Asha, 2026-10-29)."""
+    d = np.asarray(points, dtype=float) - np.asarray(centre, dtype=float)
+    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
+    along = d @ v
+    across = np.linalg.norm(d - along[:, None] * v, axis=1)
+    return 180.0 - np.degrees(np.arctan2(across, along))
+
+
+def enclosed_centroid(points, faces, normals=None):
+    """The centroid [m, (3,)] of the volume a closed triangulated surface encloses (the divergence theorem, by signed
+    tetrahedra from the points' mean, so either consistent orientation gives the same point): `faces` (nt, 3) into
+    `points`; with `normals` (nt, 3), one outward normal per triangle, each triangle is first turned to agree with its
+    own, for a surface whose triangles are not consistently ordered (the element model's). Register's angle is
+    measured from this point (`register_angle`): the body's geometric centroid, which for a body of one density is its
+    centre of mass."""
+    points, faces = np.asarray(points, dtype=float), np.asarray(faces, dtype=np.int64)
+    o = points[np.unique(faces)].mean(axis=0)
+    a, b, c = (points[faces[:, i]] - o for i in range(3))
+    n2 = np.cross(b - a, c - a)
+    v = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0
+    if normals is not None:
+        v = np.where(np.einsum("ij,ij->i", n2, np.asarray(normals, dtype=float)) < 0.0, -v, v)
+    return o + (v[:, None] * (a + b + c)).sum(axis=0) / (4.0 * v.sum())
+
+
+# the history's release angle over a step's released mass, both surface models' (sub-plan 22 Task 4): its mean by mass,
+# least and largest [deg], nan where the step released nothing
+RELEASE_ANGLE_COLUMNS = ("release_angle_register_mean_deg", "release_angle_register_min_deg",
+                         "release_angle_register_max_deg")
+
+
+def release_angle_stats(dm, angle):
+    """(mass-weighted mean, minimum, maximum) [deg] of the release angles `angle` of a step's records weighted by their
+    masses `dm` [kg], over the records that released (dm > 0); nan, nan, nan where none did (the history's
+    release_angle_register_*_deg, sub-plan 22 Task 4)."""
+    dm, angle = np.asarray(dm, dtype=float), np.asarray(angle, dtype=float)
+    on = dm > 0.0
+    if not on.any():
+        return float("nan"), float("nan"), float("nan")
+    w, a = dm[on], angle[on]
+    return float((w * a).sum() / w.sum()), float(a.min()), float(a.max())
 
 
 def histogram(r, dn, dm):
```

- [ ] **Step 7: The element model's rows and history** in `reentry_model/body.py`:

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -6,6 +6,8 @@
 from typing import Protocol
 
 import numpy as np
+
+from .spray import RELEASE_ANGLE_COLUMNS
 
 
 class Body(Protocol):
@@ -148,6 +150,13 @@
                                  # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
 NONRIGID_MAX_CROSSINGS = 64      # elements the non-rigid depth's line may cross before it stops (amendment of 2026-10-06)
 REMOVAL_NAMES = ("girin", "instant")
+
+
+def _no_release_angle():
+    """The release angle's history columns (spray.RELEASE_ANGLE_COLUMNS) of a step that released nothing: empty (nan)."""
+    return {k: float("nan") for k in RELEASE_ANGLE_COLUMNS}
+
+
 # Sub-plan 18's one-dimensional skins under the melting patches (MeltSettings(surface_model="skin1d"), skin1d.py) failed
 # their box test and are retired (sub-plan 22 Task 3, spec 9); the layered skins are skinbody.SkinBody.
 
@@ -740,7 +749,7 @@
         from . import spray as spray_mod
         if state is None or self.m_f.sum() <= 0.0 or flow is None:
             self.last_flow = self.last_spray = self.last_nonrigid = None
-            self.last_melt.update({"nonrigid_depth_wet_median_um": float("nan")})
+            self.last_melt.update({"nonrigid_depth_wet_median_um": float("nan"), **_no_release_angle()})
             return 0.0
         areas = self.surface.areas
         # the depth of liquid under each patch: its film plus the contiguous molten material beneath it. The branch
@@ -791,12 +800,18 @@
         self.thin_invalid_patch_steps += int(invalid.sum())
         self.thin_invalid_mass_max = max(self.thin_invalid_mass_max, invalid_mass)
         released = float(res.dm.sum())
+        angles = _no_release_angle()
         if released > 0.0:
             if self.spray_onset is None:
                 self.spray_onset = t
+            # Register et al.'s release angle of each patch, from the centroid of the volume this surface encloses, the
+            # surface whose centroids the rows carry (sub-plan 22 Task 4)
+            angle = spray_mod.register_angle(self.surface.centroids, spray_mod.enclosed_centroid(
+                self.mesh.points, self.surface.faces, self.surface.normals), self.v_hat)
+            angles = dict(zip(RELEASE_ANGLE_COLUMNS, spray_mod.release_angle_stats(res.dm, angle)))
             if s.particles:                  # the film is liquid and leaves at its own enthalpy (sub-plan 22 Task 2)
                 self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state,
-                                                          solid_fraction=np.zeros(b.size), h_kg=h_p)
+                                                          solid_fraction=np.zeros(b.size), h_kg=h_p, release_angle=angle)
             hn, hm = spray_mod.histogram(res.r, res.dn, res.dm)
             self.hist_n += hn
             self.hist_m += hm
@@ -864,7 +879,7 @@
             if (self.windward & (self.m_f > 0.0)).any() else float("nan"),
             "nonrigid_depth_mean_mm": nonrigid_mean, "slurry_thick_fraction": slurry_thick,
             "rigid_thin_fraction": rigid_thin, "slurry_held_mass_kg": held_mass,
-            "thin_invalid_patches": int(invalid.sum()), "thin_invalid_mass_kg": invalid_mass, **depths})
+            "thin_invalid_patches": int(invalid.sum()), "thin_invalid_mass_kg": invalid_mass, **depths, **angles})
         return released
 
     def _kill(self, dead):
@@ -1248,4 +1263,7 @@
                 "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
                 "fitted_nose_radius_mm": self.fitted_nose_radius * 1e3,
                 "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active),
-                "nonrigid_depth_wet_median_um": lm.get("nonrigid_depth_wet_median_um", float("nan"))}
+                "nonrigid_depth_wet_median_um": lm.get("nonrigid_depth_wet_median_um", float("nan")),
+                # Register et al.'s release angle over the step's released mass: its mean by mass, least and largest, nan
+                # where nothing was released (sub-plan 22 Task 4)
+                **{k: lm.get(k, float("nan")) for k in RELEASE_ANGLE_COLUMNS}}
```

- [ ] **Step 8: The last step** in `reentry_model/trajectory.py`, and the summary, the skin results and the loop in
  `reentry_model/coupled.py`:

```diff
--- a/reentry_model/trajectory.py
+++ b/reentry_model/trajectory.py
@@ -137,14 +137,18 @@
         escape.terminal, escape.direction = True, -1
         return [ground, escape]
 
-    def advance(self, dt):
+    def advance(self, dt, snap=0.0):
         """One macro step from the stored state: DOP853 over [t, t + dt] at the Step 1 tolerances, a fresh solve_ivp
         per step (spec section 7). A ground/escape event or t_max truncates the step and sets `end_reason`.
-        Returns the time actually advanced."""
+        Returns the time actually advanced. A step that would end within `snap` [s] short of t_max ends at t_max: the
+        steps' sum, t + dt added again and again, falls short of t_max by round-off, and the loop then took a last step
+        of about 1e-17 s (sub-plan 22 Task 4, fact 179). With `snap` 0, the default, the step is as it always was."""
         if self.end_reason is not None:
             raise RuntimeError("the flight already ended ({})".format(self.end_reason))
         s = self.settings
         t_target = min(self.t + dt, s.t_max)
+        if s.t_max - t_target <= snap:
+            t_target = s.t_max
         sol = solve_ivp(self.rhs, (self.t, t_target), self.y, method="DOP853", rtol=s.rtol,
                         atol=[s.atol_position] * 3 + [s.atol_velocity] * 3, events=self._events())
         if not sol.success:
```

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -14,6 +14,7 @@
 import numpy as np
 
 from . import trajectory as tj
+from .spray import RELEASE_ANGLE_COLUMNS
 from .surface_flow import KN_BODY_SHOCK
 
 THERMAL_COLUMNS = ["convective_heat_W", "rad_cooling_W", "integrated_heat_J", "absorbed_heat_J",
@@ -34,6 +35,8 @@
                 # the regime test's depth over the wet windward patches (fact 94); sub-plan 18's skin columns, which held
                 # 0 or nan in every run of the element model, went with its one-dimensional skins (sub-plan 22 Task 3)
                 "nonrigid_depth_wet_median_um"]
+# both surface models' histories also carry Register et al.'s release angle over each step's released mass (its mean by
+# mass, least and largest), spray.RELEASE_ANGLE_COLUMNS (Asha 2026-10-29, sub-plan 22 Task 4), after their own columns
 # the layered skins on the corner columns (sub-plan 20): skinbody.SkinBody.skin_stats(), every macro step
 SKINBODY_COLUMNS = ["skin_cells", "skin_thickness_max_mm", "skin_above_min", "flowing_depth_mean_mm", "flowing_depth_max_mm",
                     "liquid_depth_mean_mm", "liquid_depth_max_mm", "surface_quality_min", "surface_fold_deg", "zone_ratio_min",
@@ -48,7 +51,30 @@
                     # the droplet momentum sink (sub-plan 21 Task 7, spec 10.2): fact 29's Pi and its flag
                     "released_mass_kg", "droplet_momentum_ratio", "droplet_momentum_flagged",
                     # what the fault tolerance holds in the skins at the last base move (sub-plan 21 Task 9, spec 6.6)
-                    "held_mass_kg", "held_depth_max_mm"]
+                    "held_mass_kg", "held_depth_max_mm",
+                    # Girin's spraying from the top cells (sub-plan 22 Task 4, spec 10.2): the step's released mass -- the
+                    # released_mass_kg above is since time zero -- the droplets since time zero, the radius median and
+                    # largest over the step's records, the smallest angle and the top area of the columns whose plan
+                    # sprays, the element model's per-step columns of the same names, which particles_summary.csv reads
+                    "released_step_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg", "spraying_area_m2",
+                    # Register et al.'s release angle over the step's released mass (Asha 2026-10-29)
+                    *RELEASE_ANGLE_COLUMNS,
+                    # spec 10.2: b and the resolved sheared depth over Girin's delta_m on the wet windward columns, and
+                    # each rebuild's mass residual beside its energy residual
+                    "b_over_delta_m_median", "sheared_depth_over_delta_m_median",
+                    "rebuild_mass_residual_kg", "rebuild_mass_residual_last_kg"]
+# particles_summary.csv, one row a macro step: the element model's history columns
+SUMMARY_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um",
+                   "theta_cr_deg", "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
+# a skin flight's (sub-plan 22 Task 4): (the summary's column, the skin history's column it is read from) -- the same
+# droplet columns, the step's released mass from released_step_mass_kg, since the skin history's released_mass_kg is
+# cumulative, and where the element model has its film, which the skins do not, the skins' liquid and flowing depths
+SKIN_SUMMARY_COLUMNS = [("time_s", "time_s"), ("altitude_km", "altitude_km"), ("velocity_kms", "velocity_kms"),
+                        ("released_mass_kg", "released_step_mass_kg"), ("n_released", "n_released"),
+                        ("r_median_um", "r_median_um"), ("r_max_um", "r_max_um"), ("theta_cr_deg", "theta_cr_deg"),
+                        ("spraying_area_m2", "spraying_area_m2"), ("liquid_depth_mean_mm", "liquid_depth_mean_mm"),
+                        ("flowing_depth_mean_mm", "flowing_depth_mean_mm")]
+END_SNAP = 1e-6                   # of the macro step: a step that would end this close short of t_max ends at it (below)
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -64,6 +90,11 @@
     switch_on_melt: bool = False      # the switch also latches at the first macro step at whose end the body holds liquid
                                       # (`body.liquid_present()`), if that comes first: the layered skins' (sub-plan 22
                                       # Task 3b, Asha 2026-10-29); the element model's switch is the Knudsen number's alone
+    snap_to_t_max: bool = False       # a macro step that would end within END_SNAP of a step short of t_max ends at it:
+                                      # the steps' sum falls short of t_max by round-off, and the loop then took a last
+                                      # step of about 1e-17 s, a whole macro step of the skins over which the droplet-
+                                      # momentum ratio read 5.5 and 8.6 (fact 179); the layered skins' (sub-plan 22 Task
+                                      # 4), the element model keeping its loop, and that step, to the last bit
 
     def __post_init__(self):
         if self.dt_continuum is not None and not 0.0 < self.dt_continuum <= self.dt:
@@ -116,7 +147,7 @@
         rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
         self._frame(step, sim.t, loads)
         while sim.end_reason is None:
-            dt = sim.advance(s.dt)
+            dt = sim.advance(s.dt, snap=END_SNAP * s.dt if s.snap_to_t_max else 0.0)
             a, loads = self.loads_at(sim.t, sim.y)
             body.advance(sim.t, dt, loads, state=a)
             step += 1
@@ -200,10 +231,14 @@
                "droplet_momentum_ratio_max": float(np.max(last["droplet_momentum_ratio"])),
                "droplet_momentum_flagged_steps": int((last["droplet_momentum_flagged"] > 0).sum()),
                "held_mass_max_kg": float(np.max(last["held_mass_kg"])),          # what the fault tolerance held (Task 9)
-               "held_depth_max_mm": float(np.max(last["held_depth_max_mm"]))}
+               "held_depth_max_mm": float(np.max(last["held_depth_max_mm"])),
+               # Girin's spraying (sub-plan 22 Task 4): the droplet records' mass, the droplets and the rows written --
+               # the element model's results of the same names
+               "sprayed_mass_kg": body.sprayed_mass, "n_released": body.droplet_number,
+               "n_source_rows": body.n_source_rows}
         out.update({k: float(last[k][-1]) for k in ("removed_mass_kg", "removed_enthalpy_J", "drawn_mass_kg", "absorbed_mass_kg",
                                                   "handed_mass_kg", "rebuild_residual_J", "wall_skins_s", "wall_3d_s",
-                                                  "wall_move_s", "wall_rebuild_s")})
+                                                  "wall_move_s", "wall_rebuild_s", "wall_flow_s", "rebuild_mass_residual_kg")})
         out.update({k: int(last[k][-1]) for k in ("rebuilds", "jump_faults", "absorb_capped", "interface_repeats")})
         return out
 
@@ -329,31 +364,34 @@
     """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
     size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths.
 
-    The table is the body's `source_rows`, the element model's rows or the layered skins' blocks of them
-    (`spray.source_table`; sub-plan 22 Task 2). The summary's columns are the element model's history columns; a history
-    without them -- the skin body's, until sub-plan 22 Task 4 gives it its droplet columns -- writes no summary, and the
-    paths returned then have none."""
+    The table is the body's `source_rows`, the element model's rows, the layered skins' blocks of them (sub-plan 22 Task
+    2) or their `spray.SourceSpool` on disk (Task 4), read column by column (`spray.source_columns`), so that a flight's
+    records need not be in memory at once: the numbers written are the same either way. The summary's columns are the
+    element model's history columns (SUMMARY_COLUMNS); a skin flight's are SKIN_SUMMARY_COLUMNS, read from its history's
+    droplet columns (Task 4); a history with neither writes no summary, and the paths returned then have none."""
     import csv
     from . import spray
     os.makedirs(run_dir, exist_ok=True)
-    rows = spray.source_table(body.source_rows)
+    cols = spray.source_columns(body.source_rows)
     npz = os.path.join(run_dir, "particles.npz")
-    np.savez_compressed(npz, **{name: rows[:, i] for i, name in enumerate(spray.SOURCE_COLUMNS)})
+    np.savez_compressed(npz, **cols)
     c = history.columns
     summary = os.path.join(run_dir, "particles_summary.csv")
-    keys = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg",
-            "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
-    if not set(keys) <= set(c):
+    if set(SUMMARY_COLUMNS) <= set(c):
+        header, keys = SUMMARY_COLUMNS, SUMMARY_COLUMNS
+    elif {k for _, k in SKIN_SUMMARY_COLUMNS} <= set(c):
+        header, keys = [h for h, _ in SKIN_SUMMARY_COLUMNS], [k for _, k in SKIN_SUMMARY_COLUMNS]
+    else:
         summary = None
-    else:
+    if summary is not None:
         with open(summary, "w", newline="") as fh:
             w = csv.writer(fh)
-            w.writerow(keys)
+            w.writerow(header)
             for i in range(len(history)):
                 w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
     dist = os.path.join(run_dir, "size_distribution.csv")
     edges = spray.BIN_EDGES
-    t = rows[:, 0] if len(rows) else np.zeros(0)
+    t, r, dn, dm = cols["time_s"], cols["r_m"], cols["dn"], cols["dm_kg"]
     t_end = float(c["time_s"][-1])
     starts = np.arange(0.0, t_end + 1e-9, window)
     with open(dist, "w", newline="") as fh:
@@ -361,12 +399,12 @@
         w.writerow(["window_start_s", "window_end_s", "r_lo_m", "r_hi_m", "dn", "dm_kg"])
         for t0 in list(starts) + [None]:
             if t0 is None:
-                sel = np.ones(len(rows), dtype=bool)
+                sel = np.ones(len(t), dtype=bool)
                 lo, hi = 0.0, t_end
             else:
                 lo, hi = t0, t0 + window
                 sel = (t >= lo) & (t < hi)
-            n_hist, m_hist = spray.histogram(rows[sel, 12], rows[sel, 13], rows[sel, 14]) if sel.any() else (np.zeros(spray.N_BINS), np.zeros(spray.N_BINS))
+            n_hist, m_hist = spray.histogram(r[sel], dn[sel], dm[sel]) if sel.any() else (np.zeros(spray.N_BINS), np.zeros(spray.N_BINS))
             for j in range(spray.N_BINS):
                 w.writerow(["{:.9g}".format(lo), "{:.9g}".format(hi), "{:.9g}".format(edges[j]), "{:.9g}".format(edges[j + 1]),
                             "{:.9g}".format(n_hist[j]), "{:.9g}".format(m_hist[j])])
```

- [ ] **Step 9: The rebuild's mass residual** in `reentry_model/movemesh.py`:

```diff
--- a/reentry_model/movemesh.py
+++ b/reentry_model/movemesh.py
@@ -46,8 +46,10 @@
 by preconditioned conjugate gradients. The projection reproduces a linear field exactly and keeps the energy where it
 was, region by region, to the rule's accuracy; what the quadrature still misses is made exact by a uniform offset over
 the new deforming zone only, where the field varies, as the 2026-09-27 remeshing design does for the whole body. The
-pre-correction residual is reported with every rebuild and summed in `rebuild_residual_total`. The validity gate stops
-on a non-manifold surface, more than one connected body, or a non-positive enclosed volume.
+pre-correction residual is reported with every rebuild and summed in `rebuild_residual_total`, and so is the mass the
+rebuild's new volume misses against the old (`RebuildResult.mass_residual`, `rebuild_mass_residual_total`; spec 10.2,
+sub-plan 22 Task 4): the surface is kept, so it is round-off. The validity gate stops on a non-manifold surface, more
+than one connected body, or a non-positive enclosed volume.
 
 The fine band (Asha, 2026-10-14, option 1): on 2 mm elements the 3D model cannot follow the temperature just below the
 rigid front -- the layer ahead of the front is 2-3 mm thick on the hot plate -- and the plates melt 2-6 % late; 1 mm
@@ -308,6 +310,7 @@
         self.compression = np.zeros(len(self.corners))      # inward normal motion of each corner since the last rebuild
         self.steps_since_rebuild, self.rebuilds, self.last_rebuild = 0, 0, None
         self.rebuild_residual_total = 0.0       # J, summed pre-correction residuals of all rebuilds
+        self.rebuild_mass_residual_total = 0.0  # kg, summed mass residuals of all rebuilds (sub-plan 22 Task 4)
 
     def _set_surface(self, faces_nodes):
         local = np.full(self.mesh.n_nodes, -1, dtype=np.int64)
@@ -425,6 +428,7 @@
         h_new = project_enthalpy(mesh.points, mesh.tets, h_old, points, tets)
         params = dict(mesh.params, rebuilds=self.rebuilds + 1)
         new_mesh = mesh_mod.VolumeMesh(points, tets, params)
+        mass_residual = float(mat.rho * (mesh.volume() - new_mesh.volume()))
         check_rebuilt_boundary(new_mesh, node_of_corner[self.faces_corner])
         solver.setup(new_mesh, mat, solver.emissivity)
         solver.set_temperature(mat.temperature_from_enthalpy(h_new))
@@ -438,11 +442,12 @@
         h_new[in_zone] += residual / m_zone.sum()
         solver.set_temperature(mat.temperature_from_enthalpy(h_new))
         self.rebuild_residual_total += float(residual)
+        self.rebuild_mass_residual_total += mass_residual
         self.compression[:] = 0.0
         self.steps_since_rebuild = 0
         self.rebuilds += 1
         self.last_rebuild = RebuildResult(float(residual), float(residual / E_old) if E_old else 0.0, new_mesh.n_nodes,
-                                          new_mesh.n_elements, time.perf_counter() - t0, reason)
+                                          new_mesh.n_elements, time.perf_counter() - t0, reason, mass_residual)
         return self.last_rebuild
 
     def needs_rebuild(self):
@@ -481,6 +486,8 @@
     n_elements: int
     seconds: float
     reason: str
+    mass_residual: float = 0.0          # kg, the material's density times the old volume less the new (spec 10.2): the
+                                        # surface is kept, so round-off; nothing corrects it
 
 
 def check_surface(points, faces):
```

- [ ] **Step 10: The skin body's outputs** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -89,7 +89,17 @@
 the plan was evaluated, b its flowing depth, the droplets leaving at the mean velocity of the metal they are made of --
 the momentum the books lose with them -- and their solid fraction and enthalpy per kilogram, mass-weighted over the
 cells released), kept one block per macro step in `source_rows`; the size histograms (`hist_n`, `hist_m`) accumulate
-each flow sub-step's release at the plan's radii (`SprayPlan.histogram`, through `skin.Cascade.histogram`)."""
+each flow sub-step's release at the plan's radii (`SprayPlan.histogram`, through `skin.Cascade.histogram`).
+
+The outputs of a spraying skin flight (sub-plan 22 Task 4; spec 10.2). Each record carries Register et al.'s (2020)
+release angle, from the centroid of the volume the top surface encloses where the plan was evaluated
+(`spray.register_angle`, `spray.enclosed_centroid`; Asha, 2026-10-29); every macro step the history reports the step's
+released mass, the droplets since time zero, the radius median and largest over the step's records, the smallest angle
+and the area of the columns whose plan sprays -- the element model's per-step columns, so the coupled writer writes the
+same particles_summary.csv -- the release angle's mean by mass, least and largest, b and the resolved sheared depth over
+Girin's delta_m on the wet windward columns (`sheared_depth`), and each rebuild's mass residual. Where `spool` names a
+directory the records go to disk as they come (`spray.SourceSpool`): a flight at the continuum step makes tens of
+millions of them."""
 import copy
 import dataclasses
 import math
@@ -154,6 +164,9 @@
 
 SLURRY_SIGMA = 0.80               # N/m: the exposed slurry's surface tension in the cascade (the 2026-10-07 design 8.2, the
                                   # large-fragment design's labelled assumption for slurry)
+SHEARED_FRACTION = 0.1            # the resolved sheared depth (spec 6.7): from the top to where the computed speed falls to
+                                  # this fraction of the top cell's (sub-plan 22 Task 4): on Girin's linear profile, which
+                                  # falls to nothing at delta_m, it is 0.9 delta_m
 
 
 @dataclasses.dataclass
@@ -190,7 +203,9 @@
     `momentum_slurry` [kg m/s, (n, 3)]; and, at the step's end, `number` and `number_slurry`, the droplets: the mass over
     one droplet's at the mode's radius. Where the plan was evaluated, for the records' rows (`source_block`; sub-plan 22
     Task 2): `theta` [rad], `closure` (the surface flow's), `points` [m, (n, 3)] the top surface's corners, `state` the
-    aero state and `props` the liquid's properties."""
+    aero state and `props` the liquid's properties; and for the history and Register's release angle (Task 4), `area`
+    [m2] the columns' top areas, `centre` [m, (3,)] the centroid of the volume the top surface enclosed and `v_hat` the
+    direction of flight."""
     rate: np.ndarray
     rate_slurry: np.ndarray
     result: object
@@ -208,6 +223,9 @@
     points: np.ndarray = None
     state: object = None
     props: object = None
+    area: np.ndarray = None
+    centre: np.ndarray = None
+    v_hat: np.ndarray = None
 
     def __post_init__(self):
         T = self.slurry["T"]
@@ -264,6 +282,17 @@
         self.momentum += fb.phase_p[0]
         self.momentum_slurry += fb.phase_p[1]
 
+    def spraying(self):
+        """The columns whose plan sprays [bool, (n,)]: a mode with a release rate, the element model's `unstable`."""
+        return (self.rate > 0.0) | (self.rate_slurry > 0.0)
+
+    def release_angle(self):
+        """Register et al.'s release angle [deg, (n,)] of each column's top corner from the plan's centroid (nan without
+        one; `spray.register_angle`, sub-plan 22 Task 4)."""
+        if self.centre is None or self.points is None:
+            return np.full(len(self.rate), np.nan)
+        return spray_mod.register_angle(self.points, self.centre, self.v_hat)
+
     def source_block(self, t):
         """The macro step's records as source-table rows (spray.SOURCE_COLUMNS, by `spray.source_block`), one per column
         and mode that released, the liquid's mode's first: at the end of the step `t` [s] and the plan's aero state, the
@@ -274,9 +303,11 @@
         lose with the droplets (spec 6.4: they leave at the speed of the cells they are made of; zero direction where it
         is at rest); We_d, Oh and the breakup flag at their radius, by their phase's density, viscosity and surface tension
         (the slurry's mixture density, viscosity at its sheared layer's rate and SLURRY_SIGMA, the 2026-10-07 design 8.4);
-        and their solid fraction and enthalpy per kilogram, mass-weighted over the cells released."""
+        their solid fraction and enthalpy per kilogram, mass-weighted over the cells released; and their column's release
+        angle (Task 4)."""
         st, s = self.state, self.slurry
         flow = types.SimpleNamespace(closure=self.closure)
+        angle = self.release_angle()
         blocks = []
         for mass, H, solid, p, number, res, b, props in (
                 (self.mass, self.enthalpy, self.solid, self.momentum, self.number, self.result, self.flowing, self.props),
@@ -290,7 +321,7 @@
             rows = types.SimpleNamespace(branch=res.branch, r=res.r, delta_m=res.delta_m, we_s=res.we_s, h_floor=res.h_floor,
                                          gamma_tau=res.gamma_tau, dm=mass, dn=number, v_s=per(speed))
             blocks.append(spray_mod.source_block(t, st.h, st.V, self.theta, self.points, direction, flow, rows, b, props, st,
-                                                 solid_fraction=per(solid), h_kg=per(H)))
+                                                 solid_fraction=per(solid), h_kg=per(H), release_angle=angle))
         return np.concatenate(blocks, axis=0)
 
 
@@ -315,6 +346,35 @@
     return out
 
 
+def sheared_depth(sk, fraction=SHEARED_FRACTION):
+    """Per column, the resolved sheared depth [m] (spec 6.7; sub-plan 22 Task 4): the depth from the top at which the
+    computed speed falls to `fraction` of the top cell's, interpolated between the cells' centres as `Skins.depth_above`
+    is; the skin's thickness where it does not fall that far within the skin (the whole skin moves), nan where the top
+    cell does not move."""
+    speed = np.where(sk.valid, np.linalg.norm(sk.U, axis=2), 0.0)
+    top = speed[:, 0]
+    level = fraction * top
+    rows = np.arange(sk.n)
+    below = sk.valid & (speed <= level[:, None])
+    below[:, 0] = False
+    fell = below.any(axis=1)
+    k = np.where(fell, np.argmax(below, axis=1), 1)
+    dz = sk.dz()
+    z_mid = np.cumsum(dz, axis=1) - 0.5 * dz
+    Ka, Kb = np.clip(k - 1, 0, sk.m.shape[1] - 1), np.clip(k, 0, sk.m.shape[1] - 1)
+    sa, sb = speed[rows, Ka], speed[rows, Kb]
+    frac = np.clip(np.divide(sa - level, sa - sb, out=np.zeros_like(sa), where=sa > sb), 0.0, 1.0)
+    inside = z_mid[rows, Ka] + frac * (z_mid[rows, Kb] - z_mid[rows, Ka])
+    return np.where(top > 0.0, np.where(fell, inside, sk.thickness()), np.nan)
+
+
+def _no_step_spray():
+    """The history's droplet columns of a macro step that sprayed nothing (the element model's defaults)."""
+    out = {"r_median_um": float("nan"), "r_max_um": float("nan"), "theta_cr_deg": float("nan"), "spraying_area_m2": 0.0}
+    out.update({k: float("nan") for k in spray_mod.RELEASE_ANGLE_COLUMNS})
+    return out
+
+
 def _rows(flow, idx):
     """The surface flow's evaluation on the corners idx alone: every per-corner array taken at idx, the rest kept."""
     n = len(flow.closure)
@@ -325,7 +385,8 @@
     return out
 
 
-def girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False, points=None):
+def girin_plan(sk, flow, state, model, theta, A_top, edges, radius, dt, slurry=False, points=None, centre=None,
+               v_hat=None):
     """Girin's spraying from each column's top cells for the macro step about to run (spec 8; sub-plan 22 Task 1c), from
     the surface flow's evaluation on the corners (`flow`), the aero state (`state`, for the rarefied branch), the shared
     spray model (`model`) and the skins as they are at its start: a SprayPlan.
@@ -340,7 +401,9 @@
     of `radius` caps it too. Columns on the windward face, theta <= 90 degrees, only; `A_top` [m2] the top areas, `edges`
     the columns' graph, `dt` [s] the macro step (the nominal release in the result). With `slurry`, the exposed slurry's
     own mode too (`slurry_modes`). The plan keeps where it was evaluated for its records' rows: `theta`, the closure,
-    the aero state, the liquid's properties and `points` [m, (n, 3)], the columns' top corners (sub-plan 22 Task 2)."""
+    the aero state, the liquid's properties and `points` [m, (n, 3)], the columns' top corners (sub-plan 22 Task 2), and
+    for the history and the release angle the top areas, `centre` [m, (3,)], the centroid of the volume the top surface
+    encloses, and `v_hat`, the direction of flight (Task 4)."""
     mat, liq = sk.material, model.liquid
     n = sk.n
     F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
@@ -365,7 +428,8 @@
     else:
         rate_s, modes = np.zeros(n), _no_slurry(n)
     return SprayPlan(res.mdot, rate_s, res, modes, F, L, cap, on_slurry, theta=theta, closure=np.asarray(flow.closure),
-                     points=points, state=state, props=liq)
+                     points=points, state=state, props=liq, area=np.asarray(A_top, dtype=float), centre=centre,
+                     v_hat=None if v_hat is None else np.asarray(v_hat, dtype=float))
 
 
 def _no_slurry(n):
@@ -476,12 +540,13 @@
     with the flow on, Girin's spraying from the top cells evaluated at the start of every macro step (`girin_plan`; sub-plan
     22 Task 1c), with `slurry_spray` the exposed slurry's own mode too (the cascade), and its record in `spray_plan`; its
     droplet records, one block of source-table rows per macro step, in `source_rows` where `particles`, and the size
-    histograms in `hist_n` and `hist_m` (sub-plan 22 Task 2: the element model's interface to the coupled writer)."""
+    histograms in `hist_n` and `hist_m` (sub-plan 22 Task 2: the element model's interface to the coupled writer); with
+    `spool`, a directory, the rows go there as they come (`spray.SourceSpool`, Task 4) instead of staying in memory."""
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
-                 surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True):
+                 surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True, spool=None):
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
         self.flow_cycles = flow_cycles
@@ -499,7 +564,10 @@
         self.release, self.slurry_spray = release, bool(slurry_spray)
         self.spray, self.spray_plan = spray, None            # Girin's model, and its plan of the last macro step (Task 1c)
         self.particles = bool(particles)                     # keep the droplet records' rows (Task 2)
-        self.source_rows = []                                # one block of source-table rows per macro step (Task 2)
+        # one block of source-table rows per macro step (Task 2), in memory or, with `spool`, on disk (Task 4)
+        self.source_rows = spray_mod.SourceSpool(spool) if spool is not None and self.particles else []
+        self.sprayed_mass, self.n_source_rows = 0.0, 0         # kg of droplet records, and their rows kept (Task 4)
+        self.step_spray, self._released_step = _no_step_spray(), 0.0     # this macro step's droplet columns, release
         self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)   # droplets, kg per size bin
         self._step_release = release                         # what the parts' Flow releases this macro step
         self.repeat_rule = repeat_rule
@@ -586,6 +654,7 @@
         move called for. With the flow on, the gas side's push is the surface flow's on the corners from the aero state
         `state` (`flow_forcing`), or `forcing`, a FlowForcing given instead (a prescribed push, in the tests)."""
         n, self._ratio, self.step_parts = self._hint, 0.0, 0
+        self.step_spray, self._released_step = _no_step_spray(), 0.0
         if self.flow_on:
             end = forcing if forcing is not None else self.flow_forcing(state)
             last = self._forcing_prev
@@ -598,9 +667,12 @@
                 if end.result is None or state is None:
                     raise ValueError("Girin's spraying reads the surface flow's evaluation and the aero state (its records' "
                                      "altitude and speed): give the aero state, with a forcing only if it has its result")
-                self.spray_plan = girin_plan(self.skins, end.result, state, self.spray, self.theta, self.surface.areas,
+                surface = self.surface               # where the plan is evaluated: its corners and its centroid (Task 4)
+                self.spray_plan = girin_plan(self.skins, end.result, state, self.spray, self.theta, surface.areas,
                                              self.columns.edges, self.transverse_radius(), dt, self.slurry_spray,
-                                             points=self.surface.points.copy())
+                                             points=surface.points.copy(),
+                                             centre=spray_mod.enclosed_centroid(surface.points, surface.faces),
+                                             v_hat=self.v_hat)
                 self._step_release = self.spray_plan.release
         sk, mat = self.skins, self.material
         top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
@@ -622,12 +694,35 @@
         self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * self._ratio)))
         self._q_prev = q_end.copy()
         if self.spray is not None and self.flow_on:          # the macro step's droplets, from its accepted passes' mass
-            self.spray_plan.count(self.spray.liquid.rho)
-            self.droplet_number += float(self.spray_plan.number.sum() + self.spray_plan.number_slurry.sum())
-            if self.particles:                               # and its records' rows (Task 2)
-                block = self.spray_plan.source_block(t)
-                if len(block):
+            p = self.spray_plan
+            p.count(self.spray.liquid.rho)
+            self.droplet_number += float(p.number.sum() + p.number_slurry.sum())
+            block = p.source_block(t)                        # its records' rows (Task 2), and the history's columns
+            self.step_spray = self._step_columns(p, block)   # from them (Task 4)
+            if len(block):
+                self.sprayed_mass += float(block[:, spray_mod.SOURCE_COLUMNS.index("dm_kg")].sum())
+                if self.particles:
                     self.source_rows.append(block)
+                    self.n_source_rows += len(block)
+
+    @staticmethod
+    def _step_columns(plan, block):
+        """The history's droplet columns of a macro step from its plan and its records' rows (sub-plan 22 Task 4): the
+        radius median and largest over the rows [um], the smallest angle [deg] and the top area [m2] of the columns whose
+        plan sprays (the element model's theta_cr_deg and spraying_area_m2, over its unstable patches), and the release
+        angle's mean by mass, least and largest over the rows."""
+        out = _no_step_spray()
+        on = plan.spraying()
+        if on.any():
+            out["theta_cr_deg"] = float(np.degrees(plan.theta[on].min()))
+            out["spraying_area_m2"] = float(plan.area[on].sum())
+        if len(block):
+            C = spray_mod.SOURCE_COLUMNS
+            r = block[:, C.index("r_m")]
+            out["r_median_um"], out["r_max_um"] = float(np.median(r) * 1e6), float(r.max() * 1e6)
+            out.update(zip(spray_mod.RELEASE_ANGLE_COLUMNS, spray_mod.release_angle_stats(
+                block[:, C.index("dm_kg")], block[:, C.index("release_angle_register_deg")])))
+        return out
 
     def flow_forcing(self, state):
         """The gas side's push on the metal at this macro step (FlowForcing): the surface flow evaluated on the corners
@@ -749,6 +844,7 @@
             self.flow_courant, self.flow_growth = max(self.flow_courant, fb.courant), max(self.flow_growth, fb.growth)
             self.flow_passes, self.flow_cycles_max = max(self.flow_passes, fb.passes), max(self.flow_cycles_max, fb.cycles)
             self.released_mass += float(fb.released.sum())   # the release from the top cells (Task 7): it leaves the body
+            self._released_step += float(fb.released.sum())  # and this macro step's (sub-plan 22 Task 4)
             self.released_enthalpy += float(fb.released_H.sum())
             self.momentum_lost += fb.released_p.sum(axis=0)  # the droplets' momentum, at their cells' own velocities
             self.droplet_momentum += fb.released_p
@@ -974,7 +1070,36 @@
                 # what the fault tolerance, and the cap on a part's hand-down, hold in the skins at the last base move
                 # (Task 9, spec 6.6, 2026-10-26 and 2026-10-27): unbounded (Asha, 2026-10-26), its mass and its deepest
                 # column's depth recorded for sub-plan 22's flow-on flight to report
-                "held_mass_kg": self.held_mass, "held_depth_max_mm": self.held_depth_max * 1e3}
+                "held_mass_kg": self.held_mass, "held_depth_max_mm": self.held_depth_max * 1e3,
+                # Girin's spraying (sub-plan 22 Task 4, spec 10.2): the step's released mass, the droplets since time
+                # zero, and the step's droplet columns (`_step_columns`), which particles_summary.csv reads
+                "released_step_mass_kg": self._released_step, "n_released": self.droplet_number, **self.step_spray,
+                # b and the resolved sheared depth over Girin's delta_m on the wet windward columns (medians), and each
+                # rebuild's mass residual, summed and the last
+                **self._depth_ratios(flowing),
+                "rebuild_mass_residual_kg": mm.rebuild_mass_residual_total,
+                "rebuild_mass_residual_last_kg": mm.last_rebuild.mass_residual if mm.last_rebuild is not None else 0.0}
+
+    def _depth_ratios(self, flowing):
+        """Spec 10.2's b over delta_m and resolved sheared depth over delta_m (medians) on the wet windward columns: on
+        the windward face (theta <= 90 degrees), with flowing metal (b, the flowing depth, above zero) and Girin's
+        conjugate depth delta_m from the macro step's surface flow (`last_flow`; nan off his closure); the sheared depth
+        `sheared_depth`, on those whose top moves. nan where there are none, or without the flow."""
+        mat, out = self.material, {"b_over_delta_m_median": float("nan"), "sheared_depth_over_delta_m_median": float("nan")}
+        if not (self.flow_on and mat.melts and mat.liquid is not None) or self.last_flow is None:
+            return out
+        delta_m, _ = spray_mod.melt_layer(self.last_flow, mat.liquid)
+        if len(delta_m) != self.skins.n:
+            return out
+        with np.errstate(invalid="ignore"):
+            wet = (self.theta <= 0.5 * np.pi) & (flowing > 0.0) & np.isfinite(delta_m) & (delta_m > 0.0)
+        if wet.any():
+            out["b_over_delta_m_median"] = float(np.median(flowing[wet] / delta_m[wet]))
+            sheared = sheared_depth(self.skins)
+            moving = wet & np.isfinite(sheared)
+            if moving.any():
+                out["sheared_depth_over_delta_m_median"] = float(np.median(sheared[moving] / delta_m[moving]))
+        return out
 
     def _recession(self):
         """How far each corner's top has moved inward along its starting direction since time zero [m]."""
```

- [ ] **Step 11: The spool and the skins' last step from the command line** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -56,6 +56,10 @@
 # Sub-plan 18's one-dimensional skins are retired (sub-plan 22 Task 3).
 SURFACE_MODELS = ("skin", "elements")
 DEFAULT_SURFACE_MODEL = "skin"
+# A spraying skin run writes its droplet records to this directory inside the run's own as they come (spray.SourceSpool;
+# sub-plan 22 Task 4): at the continuum step the 100 mm flight makes 35-70 million of them, 7.6-15.2 GB in memory. The
+# coupled writer reads them from there into particles.npz, and the directory is then removed.
+SOURCE_SPOOL = "source_rows"
 # The fine band under a skin run's surface (Asha, 2026-10-14, option 1): the 3D model's elements are finer to this depth
 # while the surface keeps its triangulation, so that the 3D model follows the temperature just below the rigid front
 # (on 2 mm elements throughout the plates melted 2-6 % late), and every rebuild makes the band again
@@ -408,8 +412,9 @@
     return p
 
 
-def build_thermal(args, settings, mass):
-    """(ThermalBody or MeltingBody, HeatingModel, settings-provenance dict) for --thermal fem [--melt on]."""
+def build_thermal(args, settings, mass, spool=None):
+    """(ThermalBody or MeltingBody, HeatingModel, settings-provenance dict) for --thermal fem [--melt on]; `spool`, a
+    directory for a spraying skin body's droplet records (SOURCE_SPOOL in the run's directory, sub-plan 22 Task 4)."""
     radius = settings.diameter / 2.0
     melting = args.melt == "on"
     skins = args.surface_model == "skin"
@@ -440,7 +445,8 @@
                                      remove=melting and not flow, demise_fraction=args.demise_fraction, flow=flow,
                                      flow_cycles=args.skin_flow_cycles, surface_flow=flow_model,
                                      spray=spray_model_of(args, mat) if flow else None,
-                                     slurry_spray=args.slurry_spray == "on", particles=args.particles)
+                                     slurry_spray=args.slurry_spray == "on", particles=args.particles,
+                                     spool=spool if flow else None)
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray_model_of(args, mat)
@@ -624,14 +630,15 @@
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
     if args.thermal == "fem":
-        the_body, heating_model, thermal_info = build_thermal(args, settings, mass)
+        the_body, heating_model, thermal_info = build_thermal(args, settings, mass, spool=os.path.join(run_dir, SOURCE_SPOOL))
     sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
     if args.thermal == "fem":
         frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
         run = coupled.CoupledRun(sim, the_body, heating_model,
                                  coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk"),
                                                          dt_continuum=dt_continuum, kn_switch=args.kn_body_shock,
-                                                         switch_on_melt=args.surface_model == "skin"))
+                                                         switch_on_melt=args.surface_model == "skin",
+                                                         snap_to_t_max=args.surface_model == "skin"))
         history = run.run()
         thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
         thermal_info["dt_continuum_s"] = dt_continuum
@@ -662,8 +669,8 @@
         "files": {"csv": os.path.abspath(csv_path)},
     }
     melting = args.thermal == "fem" and args.melt == "on" and args.surface_model != "skin"
-    # the skins' spraying flight (sub-plan 22 Task 3): the element model's droplet files, from the same writer; the
-    # per-step summary waits for Task 4, which gives the skin history its droplet columns, and is recorded as not written
+    # the skins' spraying flight (sub-plan 22 Tasks 3 and 4): the element model's droplet files, from the same writer --
+    # the summary from the skin history's droplet columns -- the records read from their spool, which then goes
     spraying_skins = args.thermal == "fem" and args.melt == "on" and args.surface_model == "skin" and args.skin_flow == "on"
     if melting:
         doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()} if args.particles
@@ -671,7 +678,8 @@
         doc["files"]["melt_plots"] = [os.path.abspath(p) for p in compare.plot_melt(history, run_dir, name, reference if compare.has_melt(history, reference) else None)]
     elif spraying_skins and args.particles:
         doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()})
-        doc["files"]["particles_summary"] = None
+        if isinstance(the_body.source_rows, spray.SourceSpool):
+            the_body.source_rows.remove()
     if reference is not None:
         plots = compare.plot_all(history, reference, run_dir, name)
         doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
@@ -711,10 +719,9 @@
                       res["skin_columns"], res["removed_mass_kg"], res["drawn_mass_kg"], res["rebuilds"],
                       res["interface_repeats"], res["jump_faults"], res["wall_skins_s"], res["wall_3d_s"]))
         if spraying_skins:
-            print("  spraying: released {:.4g} kg as {:.3g} droplets in {} source rows{}; no particles_summary.csv until "
-                  "sub-plan 22 Task 4 gives the skin history its droplet columns".format(
-                      the_body.released_mass, the_body.droplet_number, sum(len(b) for b in the_body.source_rows),
-                      " (particles.npz, size_distribution.csv)" if args.particles else ""))
+            print("  spraying: released {:.4g} kg as {:.3g} droplets in {} source rows{}".format(
+                      res["released_mass_kg"], res["n_released"], res["n_source_rows"],
+                      " (particles.npz, particles_summary.csv, size_distribution.csv)" if args.particles else ""))
         if melting:
             print("  melt: onset {} km, spraying onset {} km, demise {} km at t = {} s; sprayed {:.4f} kg of {:.4f}, film left {:.4f} kg, "
                   "{:.3g} droplets (median r {} um), {} source rows, melt balance {:.1e}".format(
```

- [ ] **Step 12: Run them to see them pass.** Expected: the 15 pass (about two and a half minutes; the skin body's records test about 100 s);
  `tests/test_reentry_model_spray.py` 29 passed, `tests/test_reentry_model_trajectory.py` 13,
  `tests/test_reentry_model_coupled.py` 12, `tests/test_reentry_model_movemesh.py` 40, `tests/test_reentry_model_skinbody.py`
  36 (2 deselected, the reference tier's); the unit tier 532 passed, 1 skipped (Task 3b's 523 and nine added); the FEniCSx tier 18 passed (Task 3's 17 and the thick-film test).

- [ ] **Step 13: The measurements** (in `measure/task4/`, on mains power, one heavy run at a time): the element model's
  short case, `../task1/cmp/run_case.sh <code> "$PWD/runs" line4_short --surface-model elements` against Task 3b's
  `line3b_short` copied into `runs/`, `"$PY" -B ../task1/cmp/compare_runs.py runs line3b_short line4_short` and
  `"$PY" -B element_files.py runs line3b_short line4_short`; `./run_bitid.sh` (the skin body against Task 3b's, about 8
  minutes); `./run_smoke.sh` (two runs, about three minutes) and `"$PY" -B smoke_check.py smoke/<run>.json
  ../task3b/evidence/<run>.json`; `"$PY" -B body_outputs.py <code> body_outputs_slurry.npz on`; `"$PY" -B
  spool_cost.py <code> <scratch dir> 3000000`; `./run_fenicsx_books.sh` (about two minutes, the FEniCSx interpreter).

- [ ] **Step 14: Commit.**

```bash
git add reentry_model/spray.py reentry_model/body.py reentry_model/skinbody.py reentry_model/coupled.py reentry_model/trajectory.py reentry_model/movemesh.py reentry_model/cli.py tests/test_reentry_model_spray.py tests/test_reentry_model_trajectory.py tests/test_reentry_model_coupled.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_melting.py tests/test_reentry_model_cli.py tests/test_reentry_model_fenicsx.py
git commit -m "Skins: the outputs of a spraying skin flight (spec 8 and 10.2) -- every droplet record of both surface models carries Register et al.'s (2020) release angle from the centroid of the volume the outer surface encloses (Asha 2026-10-29), each history its mean by mass, least and largest; the skin history gains the element model's per-step droplet columns and released_step_mass_kg, so the coupled writer writes its particles_summary.csv, and b and the resolved sheared depth over delta_m and each rebuild's mass residual; the skin results the sprayed mass, droplets and rows; the skins' loop ends at t_max without a vanishing last step; a spraying skin run's records go to disk as they come and are read back to the last bit; the FEniCSx tier's thick-film test runs both solvers' gradients at 1e-13, its residual having been scikit-fem's linear tolerance; the element model's existing outputs unchanged to the last bit (sub-plan 22 Task 4)"
```

---

### Task 4b: The resolved sheared depth as the straight-line equivalent

**Files:**
- Modify: `reentry_model/skinbody.py` (`sheared_depth(sk)`; `SHEARED_FRACTION` removed)
- Test: `tests/test_reentry_model_skinbody.py` (one test replaced)

**Interfaces:**
- Consumes: Task 4's `sheared_depth`, `SkinBody._depth_ratios` and the history's `sheared_depth_over_delta_m_median`.
- Produces: `skinbody.sheared_depth(sk) -> (n,)` [m], twice the flow per unit width over the top cell's speed, nan where
  the top cell does not move; `skinbody.SHEARED_FRACTION` no longer exists. The history column keeps its name; no other
  column, field or result changes.

Asha decided Task 4's ruling 6 on 2026-10-31, asking why the depth should not be read where the speed is nearly at rest,
as Girin's profile is at delta_m. Measured on Task 4's test sphere (below): below a fast top the slurry creeps on to the
rigid front, so a speed near rest is first reached at the flowing depth b, and any fraction of the top's speed is a
choice. **The resolved sheared depth** (spec 6.7) is therefore the depth of the straight-line profile that falls from
the top cell's speed to rest and carries the column's computed flow: twice the flow per unit width -- the cells' speeds
times their depths, summed over the skin -- over the top cell's speed, with no cut-off, the fast top and the slow creep
each counted by the flow it carries. On Girin's linear profile it is delta_m, and about half a top cell more, since the
top cell's speed is its centre's (ruling 1); on a skin moving as one piece, twice its thickness. `_depth_ratios` reads it
as before, over delta_m on the wet windward columns whose top moves.

**Measured** (2026-10-31, `prototype/work-2026-10-28-skin-spraying/measure/task4b/`; fact 188):
- *The cut-offs* (`shear_cutoffs.py`, `log_shear_cutoffs.txt`). Task 4's test sphere (`bitid.py`'s body, slurry
  spraying on, seeded), on its 919 and 937 wet windward columns after each of its two 0.25 s steps, medians over delta_m
  (about 160 um; cells about 18 um): at half the top cell's speed 0.19 and 0.19; at a tenth 0.46 and 1.91 (quartiles
  0.35-0.52 and 0.86-3.16); at a hundredth 4.2 and 7.5; at a thousandth 5.3 and 10.1, against b over delta_m 5.5 and
  10.7, the median depth there the flowing depth itself. The straight-line equivalent: 0.71 and 1.47 (quartiles
  0.58-0.79 and 1.17-1.85).
- *The smoke flight* (`run_smoke.sh`, `log_smoke_compare.txt`). Task 4's 65 km case, slurry spraying off: every other
  history column of the 10 rows the same as text as Task 4's (the wall times aside); the column 2.9, 3.0, 2.5, 1.8, then
  1.08 falling to 0.72-0.78 over the last five steps, where a tenth read 2.5 falling to 0.43.

**Ruling** (2026-10-31; what, why, cost if wrong):
1. *The top cell's speed is the surface's*, not extrapolated to the surface. Why: it is the model's surface speed in
   every other use (the release, the cut-off this replaces), and an extrapolation is thrown by a profile that rises below
   the top. Cost: on Girin's profile the column reads 1 and about half a top cell over delta_m, 6 % on 18 um cells over
   160 um.

- [ ] **Step 1: The straight-line equivalent's test** -- in `tests/test_reentry_model_skinbody.py`, replacing Task 4's
  test of the tenth:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1355,25 +1355,29 @@
 
 
 # -- the outputs of a spraying skin flight (sub-plan 22 Task 4) ----------------------------------------------------------
-def test_the_resolved_sheared_depth_is_where_the_speed_falls_to_a_tenth_of_the_top_s():
-    """Spec 6.7's resolved sheared depth (`skinbody.sheared_depth`, SHEARED_FRACTION 0.1): from the top to where the
-    computed speed falls to a tenth of the top cell's, between the cells' centres. Columns of 10 um cells: a
-    profile falling linearly to rest 0.2 mm down, so 0.9 of that depth and a tenth of the top cell centre's, 5.9 um down
-    (on Girin's linear profile over delta_m, 0.9 delta_m); a skin moving as one piece, its whole thickness; a top at
-    rest, nan."""
+def test_the_resolved_sheared_depth_is_the_straight_line_depth_carrying_the_computed_flow():
+    """Spec 6.7's resolved sheared depth (`skinbody.sheared_depth`, Task 4b, Asha 2026-10-31): twice the flow per unit
+    width over the top cell's speed -- the depth of the straight line from the top cell's speed to rest that carries the
+    computed flow -- with no cut-off. Columns of 11.7 um cells: a profile falling linearly to rest 0.2 mm down, sampled at
+    the cells' centres, D / (1 - z0 / D) with z0 the top cell centre's depth, 5.9 um (on Girin's linear profile over
+    delta_m, delta_m and about half a top cell), within the sampling of the cell where it comes to rest; a skin moving as
+    one piece, twice its thickness; a fast top of five cells over a creep at a hundredth of its speed, the creep counted
+    by the flow it carries, 10.7 cells where a tenth of the top's speed read 5.4; a top at rest, nan."""
     from test_reentry_model_skin import LIQUID, RIGID, stacked, strip_grid
     col, _ = strip_grid(3, 2, h=20e-3)
-    sk = stacked(col, [[LIQUID] * 30 + [RIGID] * 10, [LIQUID] * 40] + [[LIQUID] * 40] * 4)
+    sk = stacked(col, [[LIQUID] * 30 + [RIGID] * 10, [LIQUID] * 40, [LIQUID] * 40] + [[LIQUID] * 40] * 3)
     dz = sk.dz()
     z = np.cumsum(dz, axis=1) - 0.5 * dz
     D = 2e-4
     U = np.zeros(sk.U.shape)
     U[0, :, 0] = 3.0 * np.clip(1.0 - z[0] / D, 0.0, None)
     U[1, :, 1] = 0.5
+    U[2, :, 0] = np.where(np.arange(40) < 5, 3.0, 0.03)
     sk.U = U * sk.valid[..., None]
-    assert skinbody.SHEARED_FRACTION == 0.1
+    assert not hasattr(skinbody, "SHEARED_FRACTION")
     got = skinbody.sheared_depth(sk)
-    assert got[0] == pytest.approx(0.9 * D + 0.1 * z[0, 0], rel=1e-12) and 5e-6 < z[0, 0] < 7e-6
-    assert got[1] == pytest.approx(sk.thickness()[1], rel=1e-14)
-    assert np.isnan(got[2:]).all()
-    assert skinbody.sheared_depth(sk, 0.5)[0] == pytest.approx(0.5 * D + 0.5 * z[0, 0], rel=1e-12)
+    assert got[0] == pytest.approx(D / (1.0 - z[0, 0] / D), rel=1e-4) and 5e-6 < z[0, 0] < 7e-6
+    assert got[0] - D == pytest.approx(z[0, 0], rel=0.05)
+    assert got[1] == pytest.approx(2.0 * sk.thickness()[1], rel=1e-14)
+    assert np.allclose(dz[2], dz[2, 0]) and got[2] == pytest.approx(10.7 * dz[2, 0], rel=1e-12)
+    assert np.isnan(got[3:]).all()
```

- [ ] **Step 2: Run it to see it fail.** `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -m "not drama and not
  reference" -k straight_line_depth`. Expected: 1 failed -- `AssertionError: assert not True` (`hasattr(skinbody,
  'SHEARED_FRACTION')`).

- [ ] **Step 3: The straight-line equivalent** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -164,9 +164,6 @@
 
 SLURRY_SIGMA = 0.80               # N/m: the exposed slurry's surface tension in the cascade (the 2026-10-07 design 8.2, the
                                   # large-fragment design's labelled assumption for slurry)
-SHEARED_FRACTION = 0.1            # the resolved sheared depth (spec 6.7): from the top to where the computed speed falls to
-                                  # this fraction of the top cell's (sub-plan 22 Task 4): on Girin's linear profile, which
-                                  # falls to nothing at delta_m, it is 0.9 delta_m
 
 
 @dataclasses.dataclass
@@ -346,26 +343,18 @@
     return out
 
 
-def sheared_depth(sk, fraction=SHEARED_FRACTION):
-    """Per column, the resolved sheared depth [m] (spec 6.7; sub-plan 22 Task 4): the depth from the top at which the
-    computed speed falls to `fraction` of the top cell's, interpolated between the cells' centres as `Skins.depth_above`
-    is; the skin's thickness where it does not fall that far within the skin (the whole skin moves), nan where the top
-    cell does not move."""
+def sheared_depth(sk):
+    """Per column, the resolved sheared depth [m] (spec 6.7; sub-plan 22 Task 4b, Asha 2026-10-31): the depth of the
+    straight-line profile that falls from the top cell's speed to rest and carries the column's computed flow, twice the
+    flow per unit width -- the cells' speeds times their depths, summed -- over the top cell's speed. No cut-off: below
+    a fast top the slurry creeps on to the rigid front, so a speed near rest is first reached at the flowing depth b,
+    and a fraction of the top's is a choice; here the fast top and the slow creep each count by the flow they carry. On
+    Girin's linear profile, which falls to rest at delta_m, it is delta_m and about half a top cell (the top cell's speed
+    is its centre's); on a skin moving as one piece, twice its thickness; nan where the top cell does not move."""
     speed = np.where(sk.valid, np.linalg.norm(sk.U, axis=2), 0.0)
     top = speed[:, 0]
-    level = fraction * top
-    rows = np.arange(sk.n)
-    below = sk.valid & (speed <= level[:, None])
-    below[:, 0] = False
-    fell = below.any(axis=1)
-    k = np.where(fell, np.argmax(below, axis=1), 1)
-    dz = sk.dz()
-    z_mid = np.cumsum(dz, axis=1) - 0.5 * dz
-    Ka, Kb = np.clip(k - 1, 0, sk.m.shape[1] - 1), np.clip(k, 0, sk.m.shape[1] - 1)
-    sa, sb = speed[rows, Ka], speed[rows, Kb]
-    frac = np.clip(np.divide(sa - level, sa - sb, out=np.zeros_like(sa), where=sa > sb), 0.0, 1.0)
-    inside = z_mid[rows, Ka] + frac * (z_mid[rows, Kb] - z_mid[rows, Ka])
-    return np.where(top > 0.0, np.where(fell, inside, sk.thickness()), np.nan)
+    flow = (speed * sk.dz()).sum(axis=1)
+    return np.where(top > 0.0, 2.0 * flow / np.where(top > 0.0, top, 1.0), np.nan)
 
 
 def _no_step_spray():
```

- [ ] **Step 4: Run it to see it pass.** Expected: the test passes; `tests/test_reentry_model_skinbody.py` 36 passed (2
  deselected, the reference tier's); the unit tier 532 passed, 1 skipped (Task 4's, one test replaced).

- [ ] **Step 5: The measurements** (in `measure/task4b/`, on mains power): `"$PY" -I shear_cutoffs.py <code>` (about 90 s,
  1.1 GB) and `./run_smoke.sh` (about 20 s), compared with Task 4's `smoke/` run.

- [ ] **Step 6: Commit.**

```bash
git add reentry_model/skinbody.py tests/test_reentry_model_skinbody.py
git commit -m "Skins: the resolved sheared depth (spec 6.7 and 10.2) the straight-line equivalent, twice the flow per unit width over the top cell's speed, with no cut-off (Asha 2026-10-31) -- below a fast top the slurry creeps on to the rigid front, so a cut-off near rest reads the flowing depth and a fraction of the top's speed is a choice; skinbody.SHEARED_FRACTION removed; the history column keeps its name and no other output changes (sub-plan 22 Task 4b)"
```

---

### Task 5a: A flight's whole state saved on its way, and resumed from it

**Files:**
- Create: `reentry_model/restart.py` (`save(path, state)`, `load(path)`, `FORMAT`, `PART`, `COARSE_SOLVERS`)
- Modify: `reentry_model/coupled.py` (`CoupledSettings.save_at`, `save_path`, `save_meta`; `CoupledRun.run(resume=None)`
  and `CoupledRun._save`; the results' `state_saved_at_s` and `state_save_error`)
- Modify: `reentry_model/cli.py` (`--save-state-at`, `STATE_SUFFIX`; the settings' `save_state_at_s`, the files' `state`)
- Test: `tests/test_reentry_model_restart.py` (new: 9 tests)

**Interfaces:**
- Consumes: `coupled.CoupledRun` and its loop; `spray.SourceSpool` (Task 4); the skfem backend's multigrid hierarchy
  (`solver._ml`, built afresh every `amg_rebuild_every`, 30, solves, drawing from numpy's global random state).
- Produces: `restart.save(path, state)` and `restart.load(path) -> state`, any object graph, pyamg's coarse solvers
  included; a flight's state, `{"format": 1, "run": CoupledRun, "step": int, "row": dict, "time_s": float, "rng": numpy's
  global random state, "meta": {"name": str, "args": dict}}`, written by `--save-state-at T` to
  `<outdir>/<name>.state.pkl` at the end of the first macro step ending at or after T (0 < T < --t-max), the flight
  otherwise unchanged; `CoupledRun.run(resume=state)`, the flight from the end of the saved step, its history starting
  with the saved row. Task 5b's resume command reads the state and changes a setting before resuming.

Asha's decisions of 2026-10-31 and 2026-11-01: the convergence repeats are flown over a window of about 10 s around peak
heating (93.5 s on the 100 mm flight); the model could only start a flight from the beginning, so each repeat would
have cost nearly a whole flight, two to four times dearer for a halved step or a finer surface; she chose to have the
restart scoped, then built as scoped -- the halved step, another seed and 20 um cells resumed from the base flight's
state, the 1 mm surface compared on two windows started from a uniform temperature -- with the 100 mm skin flight held on
the Mac Studio until it can save (runbook commit 1169338). **The save** pickles the coupled run whole, with numpy's
global random state, whose only draws in a run are pyamg's (`cli.cmd_run`), and the command line's arguments and run
name. Everything in it pickles but pyamg's coarse-grid solver, a class made inside `pyamg.multilevel.coarse_grid_solver`;
it is saved as its solver's name and the factor it has cached and comes back as a new solver of that name holding that
factor (ruling 1), so a hierarchy built before the save solves after it as it would have. The save alters nothing the
flight holds, and a save that fails is recorded and the flight goes on (ruling 2). **The resume** goes on from the end
of the saved step with the saved run's switch, frames' clock and step count, its history starting with the saved row
(ruling 3), its droplet records in the spool the caller gives it.

**Measured** (2026-11-01, `prototype/work-2026-10-28-skin-spraying/measure/task5/restart_scope/` and `measure/task5a/`):
- *What a flight's state holds* (`restart_scope/unpicklable.py`, the probe case on the production mesh, after one macro
  step): 1 429 MB in 127 parts, the heat solver's matrices most of it (the mass matrix 620 MB, the mesh 254 MB, the
  sparsity pattern 169 MB, the stiffness 155 MB, the hierarchy's levels 64 MB, the skins 38 MB, the deforming zone 37
  MB); the one part that fails, `solver._ml.coarse_solver`, the moving mesh failing only through the same solver.
- *A save and a resume on the production mesh* (`task5a/resume_check.py`, the probe case to 0.1 s through the flight
  launcher, saved at 0.05 s): 1.38 GB; the flight the same as the unsaved launcher run in every history column at full
  precision (8 rows) and every droplet record; resumed in a fresh process, its 5 rows the same at full precision and its
  2 311 droplet records the same in all 27 columns.
- *The test's flight* (`tests/test_reentry_model_restart.py`): two solves a macro step, so the hierarchy saved at the
  end of step 14 (28 solves) is used after the resume and a new one is built at the 31st solve, in step 16: with the
  saved random state the resumed history and records are the flight's; with the generator re-seeded instead, the
  history differs (checked by hand, 2026-11-01).

**Rulings** (2026-11-01; what, why, cost if wrong):
1. *The coarse solver is saved with its factor*, not dropped for the next solve to rebuild (as the scoping test did).
   Why: dropping it would give the base flight a fresh hierarchy at its next solve, moving it at the linear solver's
   tolerance from the save on; kept, the base flight is the unsaved one to the last bit. Cost if wrong: it reads
   pyamg's closure to refuse a solver made with arguments (which a save cannot carry); a pyamg whose internals differ
   makes the save refuse, recorded by ruling 2, not save wrongly (`test_a_coarse_grid_solver_with_arguments_is_refused`).
2. *A failed save is recorded and the flight goes on.* Why: the state is worth less than the flight. Cost: a lost state
   costs a second flight to the window, about two to three days of a core.
3. *The resumed history starts with the saved row*, not one made again. Why: the row is the record the flight wrote;
   making it again would read the step's accumulators a second time. Cost: none found.
4. *The save is the model's* (`CoupledSettings`, `--save-state-at`), not the flight launcher's, as the scoping note had
   it. Why: the command line's arguments, which a resume needs to name and describe its run, exist only there, and the
   launcher stays a wrapper that changes nothing. Cost: model code, through this plan; without `--save-state-at` every
   run is as before, its outputs unchanged.

- [ ] **Step 1: The tests** -- `tests/test_reentry_model_restart.py`:

```python
"""restart.py and the coupled run's save and resume (sub-plan 22 Task 5a; Asha 2026-11-01): a flight saves its whole
state once, at the end of the first macro step at or after a given flight time, without changing anything it does, and
a saved flight resumes from it, in an object graph of its own, to the last bit."""
import csv
import json
import os

import numpy as np
import pytest

from reentry_model import cli, restart, spray
from reentry_model import trajectory as tj

# the spraying skin flight of the command line's test (sub-plan 22 Task 3): 65 km, a body at 893 K under physics heating,
# the dev skins on the 4 mm mesh; a 25 ms first step, then 12.5 ms ones, the tops spraying from the step ending 0.0625 s
FLIGHT = ["run", "--diameter", "100", "--velocity", "7.24", "--altitude", "65", "--flight-path-angle", "-0.959331",
          "--atmosphere", "us76", "--thermal", "fem", "--heating", "physics", "--h-surface", "4", "--h-core", "20",
          "--melt", "on", "--prism-layers", "0", "--temperature", "893", "--dt", "0.025", "--skin-preset", "dev",
          "--fine-size", "2", "--skin-flow-cycles", "2", "--quiet"]


def rows_of(path):
    """A history CSV's rows as text, the wall-time columns left out (they time the run, not the flight)."""
    with open(path) as fh:
        return [{k: v for k, v in r.items() if not k.startswith("wall_")} for r in csv.DictReader(fh)]


def test_a_coarse_grid_solver_is_saved_with_its_factor(tmp_path):
    """pyamg's coarse-grid solver is a class made inside a function, which pickle cannot name; it is saved as its
    solver's name and the factor it has cached, and comes back as the same solver holding the same factor, so a
    hierarchy built before a save solves after it exactly as before."""
    pyamg = pytest.importorskip("pyamg")
    from pyamg.gallery import poisson
    A = poisson((12, 12), format="csr")
    np.random.seed(3)
    ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
    b = np.random.rand(A.shape[0])
    x = ml.solve(b, tol=1e-12)                     # the coarse solver caches its pseudo-inverse on the way
    assert hasattr(ml.coarse_solver, "P")
    path = str(tmp_path / "s.pkl")
    restart.save(path, {"ml": ml})
    assert not os.path.exists(path + restart.PART)
    got = restart.load(path)["ml"]
    assert got is not ml and repr(got.coarse_solver) == repr(ml.coarse_solver)
    np.testing.assert_array_equal(got.coarse_solver.P, ml.coarse_solver.P)
    np.testing.assert_array_equal(got.solve(b, tol=1e-12), x)


def test_a_coarse_grid_solver_with_arguments_is_refused(tmp_path):
    """A coarse solver given arguments keeps them where a save cannot reach them, so saving one is refused rather than
    bringing back a different solver."""
    pyamg = pytest.importorskip("pyamg")
    from pyamg.multilevel import coarse_grid_solver
    with pytest.raises(TypeError, match="coarse"):
        restart.save(str(tmp_path / "s.pkl"), coarse_grid_solver(("pinv", {"rtol": 1e-6})))
    assert not os.path.exists(str(tmp_path / "s.pkl"))


def test_a_flight_saved_on_its_way_is_unchanged_by_the_save_and_resumes_to_the_last_bit(tmp_path):
    """The flight to 0.2125 s twice, once saving its state at 0.1875 s: the same history and droplet records, so a save
    changes nothing in the flight that makes it. The saved state, loaded as a new object graph with the random state it
    carries and a new folder for its droplet records, resumes from the end of its 14th step: the history it makes
    starts with the saved step's row and repeats the flight's from there in every column, and its droplet records are
    the flight's after the save, to the last bit. The heat solver solves twice a step and builds a new multigrid
    hierarchy every 30 solves, drawing from numpy's random state, so the saved one is built before the save and used
    after it, and the 16th step builds a new one: both the saved hierarchy and the saved random state are tested (with
    the generator re-seeded instead, the history differs)."""
    pytest.importorskip("cantera")
    plain, saved = str(tmp_path / "plain"), str(tmp_path / "saved")
    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", plain, "--name", "f"]) == 0
    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", saved, "--name", "f", "--save-state-at", "0.1875"]) == 0
    flight = rows_of(os.path.join(plain, "f.csv"))
    assert len(flight) == 17 and rows_of(os.path.join(saved, "f.csv")) == flight
    a, b = np.load(os.path.join(plain, "f", "particles.npz")), np.load(os.path.join(saved, "f", "particles.npz"))
    assert a.files == b.files and len(a["dm_kg"]) > 0 and all(np.array_equal(a[k], b[k], equal_nan=True) for k in a.files)
    doc = json.load(open(os.path.join(saved, "f.json")))
    assert doc["settings"]["save_state_at_s"] == 0.1875 and doc["results"]["state_saved_at_s"] == pytest.approx(0.1875)
    assert doc["results"]["state_save_error"] is None and doc["files"]["state"] == os.path.join(saved, "f.state.pkl")
    assert "save_state_at_s" not in json.load(open(os.path.join(plain, "f.json")))["settings"]

    state = restart.load(doc["files"]["state"])
    assert state["format"] == restart.FORMAT and state["step"] == 14 and state["time_s"] == pytest.approx(0.1875)
    assert state["meta"]["name"] == "f" and state["meta"]["args"]["save_state_at"] == 0.1875
    run = state["run"]
    assert run.body.solver._solves == 28 and run.body.solver.amg_rebuild_every == 30
    np.random.set_state(state["rng"])
    run.body.source_rows = spray.SourceSpool(str(tmp_path / "resumed_rows"))
    run.settings.save_at = None
    history = run.run(resume=state)
    assert run.body.solver._solves == 32                          # the 31st solve built a new hierarchy
    tj.write_history_csv(history, str(tmp_path / "resumed.csv"))
    assert rows_of(str(tmp_path / "resumed.csv")) == flight[14:]
    after = a["time_s"] > state["time_s"]          # the saved step's records carry its end: 0.1875 and a round-off
    got = run.body.source_rows.columns()
    assert after.any() and all(np.array_equal(got[k], a[k][after], equal_nan=True) for k in a.files)


def test_a_save_that_fails_leaves_the_flight_flying(tmp_path, monkeypatch):
    """A save that fails (a full disk, say) is recorded and the flight goes on: the state is worth less than the flight."""
    pytest.importorskip("cantera")

    def broken(path, state):
        raise OSError("no space left on device")

    monkeypatch.setattr(restart, "save", broken)
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "f", "--save-state-at", "0.05"]) == 0
    doc = json.load(open(str(tmp_path / "f.json")))
    assert doc["results"]["end_reason"] == "t_max" and doc["results"]["state_saved_at_s"] is None
    assert "no space left on device" in doc["results"]["state_save_error"] and "state" not in doc["files"]
    assert len(rows_of(str(tmp_path / "f.csv"))) == 5


@pytest.mark.parametrize("argv", [
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "0"],
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "-0.05"],
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "0.0875"],     # nothing would be left to resume
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "nan"],
    ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.5", "--t-max", "10", "--save-state-at", "5"],
])
def test_bad_save_times_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_restart.py -q`. Expected: a
  collection error -- `ImportError: cannot import name 'restart' from 'reentry_model'`.

- [ ] **Step 3: The save and the load** -- `reentry_model/restart.py`:

```python
"""A flight's whole state on disk, and back (sub-plan 22 Task 5a; Asha 2026-11-01).

The convergence repeats of Task 5 -- a halved macro step, another seed, 20 instead of 10 um skin cells -- are compared
over a window of about 10 s around peak heating, 93.5 s into the 100 mm flight; flown from the start each would cost
nearly a whole flight. So the base flight saves its state once on its way (`coupled.CoupledSettings.save_at`, the
command line's --save-state-at) and each repeat resumes from it (`coupled.CoupledRun.run(resume=...)`).

The state is the coupled run itself -- the trajectory, the body with its 3D mesh, heat solver, skins, moving mesh and
spray model -- with numpy's global random state, whose only draws in a run are pyamg's (cli.cmd_run), pickled whole.
Measured on the 100 mm sphere's production mesh (2026-11-01): 1.32 GB, written in 0.4 s, and everything in it pickles but
pyamg's coarse-grid solver, a class made inside `pyamg.multilevel.coarse_grid_solver`, which pickle cannot name. It is
saved as the name of its solver and the factor it has cached (`P`, `LU` or `L`), and comes back as a new solver of that
name holding that factor, so a hierarchy built before the save solves after it exactly as it would have: a flight
resumed from its state in a fresh process repeated the unbroken one to the last bit. The objects the flight holds are
never altered by a save."""
import ast
import os
import pickle

FORMAT = 1                      # of a flight's saved state (`coupled.CoupledRun.state`); a resume refuses any other
PART = ".part"                  # the state is written beside its path under this suffix, then renamed into place
COARSE_SOLVERS = ("pinv", "pinv2", "lu", "cholesky")    # pyamg's coarse solvers whose caches are plain arrays


def _is_coarse_solver(obj):
    cls = type(obj)
    return cls.__module__ == "pyamg.multilevel" and cls.__qualname__.endswith("coarse_grid_solver.<locals>.GenericSolver")


def _coarse_name(obj):
    """The solver a pyamg coarse-grid solver was made with, which must be one of COARSE_SOLVERS given no arguments: the
    arguments live in the closure that made it, where a save cannot carry them, so a solver with any is refused."""
    try:
        name = ast.literal_eval(type(obj).name())
        cells = dict(zip(type(obj).__call__.__code__.co_freevars, type(obj).__call__.__closure__ or ()))
        solve = cells["solve"].cell_contents
        kwargs = dict(zip(solve.__code__.co_freevars, solve.__closure__ or ())).get("kwargs")
        kwargs = {} if kwargs is None else kwargs.cell_contents
    except (KeyError, ValueError, SyntaxError, AttributeError, TypeError):
        raise TypeError("a pyamg coarse solver whose make cannot be read cannot be saved: {!r}".format(obj))
    if name not in COARSE_SOLVERS or kwargs:
        raise TypeError("a pyamg coarse solver {!r} with arguments {!r} cannot be saved (saved: {} with none)".format(
            name, kwargs, ", ".join(COARSE_SOLVERS)))
    return name


def _coarse_solver(name, cache):
    """A pyamg coarse-grid solver of `name` holding `cache`, the attributes the saved one had cached."""
    from pyamg.multilevel import coarse_grid_solver
    solver = coarse_grid_solver(name)
    vars(solver).update(cache)
    return solver


class _Pickler(pickle.Pickler):
    def reducer_override(self, obj):
        if _is_coarse_solver(obj):
            return _coarse_solver, (_coarse_name(obj), dict(vars(obj)))
        return NotImplemented


def save(path, state):
    """Write `state` to `path`: pickled (protocol 5), the coarse solvers as above, into path + PART and then renamed, so
    a save cut short leaves no file at `path` (and its part removed when the failure is Python's)."""
    part = path + PART
    try:
        with open(part, "wb") as fh:
            _Pickler(fh, protocol=5).dump(state)
        os.replace(part, path)
    finally:
        if os.path.exists(part):
            os.remove(part)


def load(path):
    """The state `save` wrote at `path`."""
    with open(path, "rb") as fh:
        return pickle.load(fh)
```

- [ ] **Step 4: The save on the way, and the resume**, in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -13,6 +13,7 @@
 
 import numpy as np
 
+from . import restart
 from . import trajectory as tj
 from .spray import RELEASE_ANGLE_COLUMNS
 from .surface_flow import KN_BODY_SHOCK
@@ -95,10 +96,16 @@
                                       # step of about 1e-17 s, a whole macro step of the skins over which the droplet-
                                       # momentum ratio read 5.5 and 8.6 (fact 179); the layered skins' (sub-plan 22 Task
                                       # 4), the element model keeping its loop, and that step, to the last bit
+    save_at: float = None             # s of flight: the run saves its whole state once, at the end of the first macro
+                                      # step ending at or after it, to save_path (`restart.save`; sub-plan 22 Task 5a)
+    save_path: str = None
+    save_meta: dict = None            # kept with the state: what made the run (the command line's arguments and name)
 
     def __post_init__(self):
         if self.dt_continuum is not None and not 0.0 < self.dt_continuum <= self.dt:
             raise ValueError("dt_continuum must be > 0 and no longer than dt, got {!r}".format(self.dt_continuum))
+        if self.save_at is not None and not (self.save_at > 0.0 and self.save_path):
+            raise ValueError("save_at must be > 0 and needs a save_path, got {!r} and {!r}".format(self.save_at, self.save_path))
 
 
 class CoupledRun:
@@ -110,6 +117,7 @@
         if self.settings.frames_every and not self.settings.output_dir:
             raise ValueError("frames_every > 0 needs an output_dir")
         self.switched_at = self.switch_trigger = None              # the continuum switch: (t, Kn, altitude km), and why
+        self.state_saved_at = self.state_save_error = None        # the save's flight time [s], or why it failed (Task 5a)
 
     def loads_at(self, t, y):
         a = self.sim.aero_state(t, y[:3], y[3:])
@@ -138,14 +146,22 @@
     def skins(self):
         return hasattr(self.body, "skin_stats")
 
-    def run(self):
+    def run(self, resume=None):
+        """Fly to the end, from the start or, with `resume` -- a state this run saved (`settings.save_at`), read back by
+        `restart.load` (sub-plan 22 Task 5a) -- from the end of the saved macro step: the history then starts with that
+        step's row, and the continuum switch, the frames' clock and the step count are the saved run's."""
         started = time.perf_counter()
         sim, body, s = self.sim, self.body, self.settings
-        self.switched_at = self.switch_trigger = None              # (t, Kn, altitude km) of the continuum switch, and why
-        self._t0, self._frame_dt, self._next_frame = sim.t, s.frames_every * s.dt, sim.t
-        a, loads = self.loads_at(sim.t, sim.y)
-        rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
-        self._frame(step, sim.t, loads)
+        if resume is None:
+            self.switched_at = self.switch_trigger = None          # (t, Kn, altitude km) of the continuum switch, and why
+            self._t0, self._frame_dt, self._next_frame = sim.t, s.frames_every * s.dt, sim.t
+            a, loads = self.loads_at(sim.t, sim.y)
+            rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
+            self._frame(step, sim.t, loads)
+        else:
+            if resume.get("format") != restart.FORMAT or resume.get("run") is not self:
+                raise ValueError("a run resumes only from a state of its own of format {}".format(restart.FORMAT))
+            rows, states, step = [dict(resume["row"])], [sim.y.copy()], resume["step"]
         while sim.end_reason is None:
             dt = sim.advance(s.dt, snap=END_SNAP * s.dt if s.snap_to_t_max else 0.0)
             a, loads = self.loads_at(sim.t, sim.y)
@@ -159,6 +175,9 @@
             self._frame(step, sim.t, loads)
             if (self.melting or self.skins) and body.demised():
                 sim.end_reason = "demise"
+            if (s.save_at is not None and sim.end_reason is None and self.state_saved_at is None
+                    and self.state_save_error is None and sim.t >= s.save_at - END_SNAP * s.dt):
+                self._save(step, rows[-1])
         self._collection()
         columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
         history = tj.History(columns, np.array(states), sim.end_reason)
@@ -174,7 +193,22 @@
                                     "dt_switch_altitude_km": h_sw})
             if s.switch_on_melt:                                   # what latched it: knudsen, melting or both at once
                 history.results["dt_switch_trigger"] = self.switch_trigger
+        if s.save_at is not None:
+            history.results.update({"state_saved_at_s": self.state_saved_at, "state_save_error": self.state_save_error})
         return history
+
+    def _save(self, step, row):
+        """Save the whole run once, at the end of macro step `step`, whose history row is `row` (sub-plan 22 Task 5a):
+        the run, the step, its row and flight time, numpy's global random state and `settings.save_meta`, by
+        `restart.save`, which alters nothing the run holds. A save that fails is recorded and the flight goes on: its
+        state is worth less than the flight."""
+        s = self.settings
+        try:
+            restart.save(s.save_path, {"format": restart.FORMAT, "run": self, "step": step, "row": row,
+                                       "time_s": float(self.sim.t), "rng": np.random.get_state(), "meta": s.save_meta})
+            self.state_saved_at = float(self.sim.t)
+        except Exception as exc:
+            self.state_save_error = "{}: {}".format(type(exc).__name__, exc)
 
     def _switch(self, t, a):
         """The continuum step (amendment of 2026-10-07): from the first step whose aero state -- the one the body was
```

- [ ] **Step 5: `--save-state-at`** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -60,6 +60,9 @@
 # sub-plan 22 Task 4): at the continuum step the 100 mm flight makes 35-70 million of them, 7.6-15.2 GB in memory. The
 # coupled writer reads them from there into particles.npz, and the directory is then removed.
 SOURCE_SPOOL = "source_rows"
+# --save-state-at saves the whole coupled run once on its way to <outdir>/<name> + this (restart.py; sub-plan 22 Task 5a):
+# about 1.3 GB on the 100 mm sphere's production mesh, from which a convergence repeat resumes.
+STATE_SUFFIX = ".state.pkl"
 # The fine band under a skin run's surface (Asha, 2026-10-14, option 1): the 3D model's elements are finer to this depth
 # while the surface keeps its triangulation, so that the 3D model follows the temperature just below the rigid front
 # (on 2 mm elements throughout the plates melted 2-6 % late), and every rebuild makes the band again
@@ -300,6 +303,9 @@
     th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
     th.add_argument("--dt", type=float, default=0.5, help="macro step [s] (default %(default)s)")
     th.add_argument("--frames-every", type=int, default=0, help="VTK frame every n macro steps (default 0: none; 10 with --animate/--stills)")
+    th.add_argument("--save-state-at", type=float, default=None,
+                    help="[s] of flight: save the whole run once, at the end of the first macro step ending at or after it, "
+                         "as <outdir>/<name>{}, for a run to resume from (sub-plan 22 Task 5a)".format(STATE_SUFFIX))
     th.add_argument("--animate", action="store_true", help="MP4/GIF of the surface temperature and of the meridional cross-section, plus stills")
     th.add_argument("--stills", action="store_true", help="only the stills (start, peak heating, peak surface T, end; surface and section)")
     me = r.add_argument_group("melting and spraying (Step 3)")
@@ -559,6 +565,10 @@
         parser.error("--animate/--stills/--frames-every need --thermal fem")
     if args.melt == "on" and args.thermal != "fem":
         parser.error("--melt on needs --thermal fem")
+    if args.save_state_at is not None and args.thermal != "fem":
+        parser.error("--save-state-at needs --thermal fem: it saves a coupled run")
+    if args.save_state_at is not None and not 0.0 < args.save_state_at < args.t_max:
+        parser.error("--save-state-at must be within (0, --t-max): a state saved at the end leaves nothing to resume")
     if not 1.0 < args.gamma_pm < 2.0:
         parser.error("--gamma-pm must be within (1, 2)")
     for label, value in (("--we-critical", args.we_critical), ("--kr", args.kr), ("--kt", args.kt), ("--layer-thickness", args.layer_thickness),
@@ -634,14 +644,20 @@
     sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
     if args.thermal == "fem":
         frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
+        saving = args.save_state_at is not None
+        state_path = os.path.join(args.outdir, name + STATE_SUFFIX) if saving else None
         run = coupled.CoupledRun(sim, the_body, heating_model,
                                  coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk"),
                                                          dt_continuum=dt_continuum, kn_switch=args.kn_body_shock,
                                                          switch_on_melt=args.surface_model == "skin",
-                                                         snap_to_t_max=args.surface_model == "skin"))
+                                                         snap_to_t_max=args.surface_model == "skin",
+                                                         save_at=args.save_state_at, save_path=state_path,
+                                                         save_meta={"name": name, "args": dict(vars(args))} if saving else None))
         history = run.run()
         thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
         thermal_info["dt_continuum_s"] = dt_continuum
+        if saving:
+            thermal_info["save_state_at_s"] = args.save_state_at
     else:
         history = sim.run(extra_times=reference.time if reference is not None else None)
     csv_path = os.path.join(args.outdir, name + ".csv")
@@ -668,6 +684,8 @@
         "provenance": prov,
         "files": {"csv": os.path.abspath(csv_path)},
     }
+    if history.results.get("state_saved_at_s") is not None:
+        doc["files"]["state"] = os.path.abspath(state_path)
     melting = args.thermal == "fem" and args.melt == "on" and args.surface_model != "skin"
     # the skins' spraying flight (sub-plan 22 Tasks 3 and 4): the element model's droplet files, from the same writer --
     # the summary from the skin history's droplet columns -- the records read from their spool, which then goes
```

- [ ] **Step 6: Run them to see them pass.** Expected: 9 passed in about 2 minutes; the unit tier 541 passed, 1 skipped (Task 4b's 532 and these
  9; 10 minutes across four processes with `-n 4`, measured 2026-11-01).

- [ ] **Step 7: The production check** (in `measure/task5a/`, on mains power): the probe case to 0.1 s through
  `flights/launch.sh` with `--save-state-at 0.05` (about 5 minutes, 4.3 GB), then `"$PY" -B resume_check.py <run dir>
  saved_t01`: no column or record differing, from the unsaved launcher run nor after the resume.

- [ ] **Step 8: Commit.**

```bash
git add reentry_model/restart.py reentry_model/coupled.py reentry_model/cli.py tests/test_reentry_model_restart.py
git commit -m "Coupled runs: a flight saves its whole state once on its way (--save-state-at T, at the end of the first macro step ending at or after T, to <outdir>/<name>.state.pkl) and resumes from it to the last bit (CoupledRun.run(resume=...)), for Task 5's convergence repeats over a window at peak heating (Asha 2026-11-01) -- restart.py pickles the run, numpy's random state and the command line's arguments, pyamg's coarse solver as its name and cached factor; the save alters nothing in the flight, a failed one is recorded and the flight goes on; without the option every output is as before (sub-plan 22 Task 5a)"
```

---

### Task 5b: The resume command

**Files:**
- Modify: `reentry_model/cli.py` (`resume_arguments`, `load_state`, `cmd_resume`, the `resume` subcommand; the writing of a
  run's csv, json and droplet files moved out of `cmd_run` into `finish_run`, which both call)
- Modify: `flights/fly.py` (beside `code/`, not part of the planning copy's code): `fly.py <outdir> <name> resume <state>
  [options]` runs the resume command, its rows and progress kept as for `run`
- Test: `tests/test_reentry_model_restart.py` (3 tests, 8 cases, beside Task 5a's 9)
- Create (git-ignored): `prototype/work-2026-10-28-skin-spraying/measure/task5b/`, the production check

**Interfaces:**
- Consumes: Task 5a's `restart.load`, `restart.FORMAT` and `CoupledRun.run(resume=state)`; `spray.SourceSpool` (Task 4);
  the state's `{format, run, step, row, time_s, rng, meta{name, args}}`.
- Produces: `python -m reentry_model resume STATE --outdir DIR --name NAME [--t-max S] [--dt-continuum S] [--seed N]
  [--frames-every N] [--slurry-spray on|off]`, writing the files `run` writes (`<outdir>/<name>.csv` and `.json`, the droplet
  files under `<outdir>/<name>/`); the JSON gains `"resumed_from": {"state", "step", "time_s", "run_name", "changed"}`.
  The options are in `resume_arguments(parser)` and applied in one block of `cmd_resume`, between loading the state and
  `run.run(resume=state)`, so a later task adds one option in those two places.

Asha's decision of 2026-11-01 was to resume the convergence repeats from the base flight's state; Task 5a made the state
and this task puts the resume on the command line, which until now could only be done from a script
(`measure/task5a/resume_check.py`). **The random state** is the saved one unless `--seed` is given, which seeds numpy's
generator instead. **What the run is given:** `save_at`, `save_path` and `save_meta` are cleared (a resumed flight saves
nothing); `output_dir` is `<outdir>/<name>/vtk` and `frames` is emptied, so the new series lists only its own frames (with
`--frames-every` the frames' step and the clock the series follows after the continuum switch are set together); `--t-max`
sets the simulator's `t_max`; `--dt-continuum` sets the continuum step, and the macro step in use too if the switch has
latched; `--slurry-spray` sets the body's flag; the droplet records go to a new spool in `<outdir>/<name>/source_rows`, so
only those after the resume are written. **The files** are written by the code `run` used, now `finish_run`, from the saved
flight's own arguments (`state["meta"]["args"]`) with the options given here, so `run`'s outputs are unchanged. Bad input
exits 2: a missing file, a state of another format or without the command line's arguments, `--t-max` not finite or not
after the saved time, `--dt-continuum` not finite or outside (0, the saved `--dt`], or given for a flight that has no
continuum switch.

**Measured** (2026-11-06): the two real states
saved by Task 5a's build (commit 055c313) resumed through `flights/launch.sh` on mains power (`measure/task5b/check.md`).
The 50 mm flight's state (0.29 GB, saved at 170.5 s, step 341) resumed to 172.0 s in 190 s (3 macro steps of 0.5 s, peak
memory 1.47 GB), and the 100 mm flight's (1.38 GB, saved at 22.0 s, step 44) to 23.0 s in 330 s (2 steps, 3.95 GB): in both,
84 history columns -- every column but the five `wall_` timers -- equal the snapshot flights' own rows at full precision,
and the column sets are the same, the saved states pickled before this task resuming with the command unchanged. Neither
flight had released a droplet yet (0 source rows), so the droplet records are pinned by the unit test only. The unit
tier is 549 passed, 1 skipped.

**Rulings** (2026-11-06; what, why, cost if wrong):
1. *The resume writes `cmd_run`'s files by sharing its writer*, not by a copy. Why: two writers would drift and `run`'s
   outputs must not change (`finish_run` is `cmd_run`'s tail moved whole, with the JSON's `resumed_from` added only when
   resuming). Cost: none found; the unit tier's `run` tests pass unchanged.
2. *`--outdir` and `--name` are required.* Why: a resume that defaulted to the saved flight's directory and name would
   overwrite the files the state came from. Cost: two more words on the command line.
3. *`--dt-continuum` is bounded by the saved `--dt`, the macro step before the switch, not by the step now in use.* Why: it
   is the rule `run` applies, and a coarser repeat than the base flight's continuum step is a legitimate one. Cost: none.
4. *The JSON's `settings.thermal` block of a resumed run is short* -- the heating, melt, surface model, solvers, skin-flow,
   slurry-spray settings, the macro and continuum steps and `frames_every` -- not the mesh and skin details `run`
   records. Why: those come from `build_thermal`, which makes the 3D mesh (minutes, GBs) and is not in the state. Cost: a
   resumed run's JSON lacks `n_nodes` and the like; the saved flight's JSON has them.
5. *`results.n_source_rows` and the other body totals of a resumed run count from the start of the flight, the droplet
   files only from the resume.* Why: the totals are the body's own cumulative accumulators, in the state. Cost: the printed
   count of source rows is higher than the rows in `particles.npz`; the droplet columns of the history agree with the body.

- [ ] **Step 1: The tests** -- in `tests/test_reentry_model_restart.py`:

```diff
--- a/tests/test_reentry_model_restart.py
+++ b/tests/test_reentry_model_restart.py
@@ -121,3 +121,72 @@
     with pytest.raises(SystemExit) as exc:
         cli.main(argv + ["--outdir", str(tmp_path)])
     assert exc.value.code == 2
+
+
+@pytest.fixture(scope="module")
+def flown(tmp_path_factory):
+    """The test flight to 0.2125 s once, saving its state at 0.1875 s (the end of its 14th step), for the resume command's
+    tests: (the flight's directory, its history rows as text, its droplet records, the saved state's path)."""
+    pytest.importorskip("cantera")
+    out = str(tmp_path_factory.mktemp("flown"))
+    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", out, "--name", "f", "--save-state-at", "0.1875"]) == 0
+    rows = rows_of(os.path.join(out, "f.csv"))
+    assert len(rows) == 17
+    return out, rows, dict(np.load(os.path.join(out, "f", "particles.npz"))), os.path.join(out, "f.state.pkl")
+
+
+def test_the_resume_command_repeats_the_flight_to_the_last_bit_and_writes_its_files(flown, tmp_path):
+    """`resume` on the state saved at 0.1875 s, with nothing changed: it writes the files `run` writes under its own name
+    -- the history, the JSON, the droplet files -- the history starting with the saved step's row and repeating the
+    unbroken flight's from there in every column, the droplet records those of the flight after the save, to the last bit,
+    and the JSON naming the state it came from with nothing changed (sub-plan 22 Task 5b, Asha 2026-11-06)."""
+    out, flight, records, state_path = flown
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r"]) == 0
+    assert rows_of(str(tmp_path / "r.csv")) == flight[14:]
+    after = records["time_s"] > 0.1875 + 1e-9        # the saved step's records carry its end: 0.1875 and a round-off
+    got = np.load(str(tmp_path / "r" / "particles.npz"))
+    assert after.any() and got.files == list(records) and all(np.array_equal(got[k], records[k][after], equal_nan=True) for k in records)
+    assert os.path.exists(str(tmp_path / "r" / "particles_summary.csv")) and not os.path.exists(str(tmp_path / "r" / cli.SOURCE_SPOOL))
+    doc = json.load(open(str(tmp_path / "r.json")))
+    assert doc["run_name"] == "r" and doc["settings"]["t_max_s"] == 0.2125
+    assert doc["resumed_from"] == {"state": os.path.abspath(state_path), "step": 14, "time_s": pytest.approx(0.1875),
+                                   "run_name": "f", "changed": {}}
+    assert "resumed_from" not in json.load(open(os.path.join(out, "f.json")))
+    assert doc["results"]["end_reason"] == "t_max" and "state_saved_at_s" not in doc["results"]
+
+
+def test_the_resume_command_changes_what_it_is_told(flown, tmp_path):
+    """Each option changes the resumed flight and is recorded in the JSON's "changed": --t-max flies further (two more
+    12.5 ms steps than the saved flight's end), a halved --dt-continuum takes twice the steps over the same time (4 from
+    0.1875 s to 0.2125 s instead of 2), and --seed seeds numpy's generator instead of resuming the saved random state, so
+    the multigrid hierarchy the heat solver builds in the 16th step is another and the history differs from the flight's."""
+    out, flight, records, state_path = flown
+    run = lambda name, *options: (cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", name, *options]),
+                                  rows_of(str(tmp_path / (name + ".csv"))), json.load(open(str(tmp_path / (name + ".json")))))
+    code, longer, doc = run("longer", "--t-max", "0.2375")
+    assert code == 0 and len(longer) == 5 and longer[:2] == flight[14:16] and doc["resumed_from"]["changed"] == {"--t-max": 0.2375}
+    assert doc["results"]["final_time_s"] == pytest.approx(0.2375)       # (its 3rd row is mid-flight, the flight's the last: only its momentum_residual differs)
+    code, finer, doc = run("finer", "--dt-continuum", "0.00625")
+    assert code == 0 and len(finer) == 5 and doc["resumed_from"]["changed"] == {"--dt-continuum": 0.00625}
+    assert doc["results"]["final_time_s"] == pytest.approx(0.2125) and doc["settings"]["thermal"] == "fem"
+    code, reseeded, doc = run("reseeded", "--seed", "7")
+    assert code == 0 and len(reseeded) == 3 and reseeded != flight[14:] and doc["resumed_from"]["changed"] == {"--seed": 7}
+    assert doc["settings"]["seed"] == 7
+    code, same, doc = run("same")
+    assert same == flight[14:] and doc["settings"]["seed"] == cli.DEFAULT_SEED
+
+
+@pytest.mark.parametrize("case", ["missing", "format", "t_max", "dt_continuum", "t_max_nan", "dt_continuum_nan"])
+def test_the_resume_command_refuses_what_it_cannot_resume(case, flown, tmp_path):
+    """A missing state file, a state of another format, an end of flight not after the saved time, a continuum step outside
+    (0, the saved flight's --dt of 0.025 s] and a value that is not a number are refused with exit 2 (bad input), and
+    nothing is written."""
+    out, flight, records, state_path = flown
+    bad_format = str(tmp_path / "other.pkl")
+    restart.save(bad_format, {"format": restart.FORMAT + 1})
+    state, options = {"missing": (str(tmp_path / "none.pkl"), []), "format": (bad_format, []),
+                      "t_max": (state_path, ["--t-max", "0.1875"]), "dt_continuum": (state_path, ["--dt-continuum", "0.05"]),
+                      "t_max_nan": (state_path, ["--t-max", "nan"]), "dt_continuum_nan": (state_path, ["--dt-continuum", "nan"])}[case]
+    with pytest.raises(SystemExit) as exc:
+        cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", *options])
+    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_restart.py -q -k resume_command`.
  Expected: failures with `argument command: invalid choice: 'resume' (choose from run, compare)`.

- [ ] **Step 3: The command** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -7,6 +7,8 @@
         [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off ...]
                                      (the layered skins, the default with --thermal fem, their flow and spraying the
                                       default with --melt on; the element model with --surface-model elements)
+    python -m reentry_model resume <state> --outdir ... --name ... [--t-max S --dt-continuum S --seed N --frames-every N
+        --slurry-spray on|off]       (a flight saved with --save-state-at, continued; sub-plan 22 Task 5b)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
 Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
@@ -22,7 +24,7 @@
 
 import numpy as np
 
-from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, movemesh, sesam_io, skin, skinbody, spray, surface_flow, thermal, viz
+from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, movemesh, restart, sesam_io, skin, skinbody, spray, surface_flow, thermal, viz
 from . import trajectory as tj
 from .earth import GRAVITY_MODELS
 from .thermal import MissingBackend
@@ -242,6 +244,22 @@
             versions[key] = None
     return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
             "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0], **versions}
+
+
+def resume_arguments(rs):
+    """The resume command's options, in one place: the new run's files and what it may change from the saved flight. A
+    flight resumed with none of the changes repeats the saved flight from its saved step to the last bit."""
+    rs.add_argument("--outdir", required=True, help="directory of the resumed run's files (never the saved flight's by default)")
+    rs.add_argument("--name", required=True, help="the resumed run's name: <outdir>/<name>.csv and .json, droplet files under <outdir>/<name>/")
+    rs.add_argument("--t-max", type=float, default=None, help="[s] end of flight, after the saved time (default: the saved flight's)")
+    rs.add_argument("--dt-continuum", type=float, default=None,
+                    help="[s] macro step after the continuum switch, in (0, the saved --dt] (default: the saved flight's)")
+    rs.add_argument("--seed", type=parse_seed, default=None,
+                    help="seed numpy's generator with this instead of resuming from the saved random state (default: the saved state)")
+    rs.add_argument("--frames-every", type=int, default=None,
+                    help="VTK frame every n macro steps, the series new (default: the saved flight's setting)")
+    rs.add_argument("--slurry-spray", choices=("on", "off"), default=None,
+                    help="the exposed slurry sprays too (default: the saved flight's)")
 
 
 def build_parser():
@@ -409,6 +427,11 @@
                     help="write the particle source table and the size distribution of a spraying run (default)")
     me.add_argument("--no-particles", dest="particles", action="store_false")
 
+    rs = sub.add_parser("resume", help="continue a flight from a state it saved (--save-state-at), to the last bit, "
+                                       "changing only what is given (sub-plan 22 Task 5b)")
+    rs.add_argument("state", help="the saved state, <outdir>/<name>{}".format(STATE_SUFFIX))
+    resume_arguments(rs)
+
     c = sub.add_parser("compare", help="metrics and plots for an existing model history")
     c.add_argument("--model", required=True, help="model history CSV")
     c.add_argument("--reference", required=True, help="SESAM run CSV")
@@ -660,6 +683,16 @@
             thermal_info["save_state_at_s"] = args.save_state_at
     else:
         history = sim.run(extra_times=reference.time if reference is not None else None)
+    return finish_run(args, name, settings, mass, atm_name, atm_info, reference, thermal_info, the_body, history,
+                      state_path if args.thermal == "fem" else None)
+
+
+def finish_run(args, name, settings, mass, atm_name, atm_info, reference, thermal_info, the_body, history, state_path,
+               resumed=None):
+    """Write a finished run's files and print its summary (cmd_run and cmd_resume share this, sub-plan 22 Task 5b):
+    <outdir>/<name>.csv and .json, the droplet files and plots under <outdir>/<name>/. `resumed`, the resume command's
+    record of its state, goes into the JSON as "resumed_from". Returns the exit code."""
+    run_dir = os.path.join(args.outdir, name)
     csv_path = os.path.join(args.outdir, name + ".csv")
     json_path = os.path.join(args.outdir, name + ".json")
     tj.write_history_csv(history, csv_path)
@@ -684,6 +717,8 @@
         "provenance": prov,
         "files": {"csv": os.path.abspath(csv_path)},
     }
+    if resumed is not None:
+        doc["resumed_from"] = resumed
     if history.results.get("state_saved_at_s") is not None:
         doc["files"]["state"] = os.path.abspath(state_path)
     melting = args.thermal == "fem" and args.melt == "on" and args.surface_model != "skin"
@@ -768,6 +803,96 @@
     return 0
 
 
+def load_state(path, parser):
+    """The saved flight at `path`, or a parser error (exit 2) for a missing file, one that is no state of this model's
+    format, or one saved without the command line's arguments."""
+    if not os.path.isfile(path):
+        parser.error("no saved state at {}".format(path))
+    try:
+        state = restart.load(path)
+    except Exception as exc:
+        parser.error("{} is not a saved state: {}: {}".format(path, type(exc).__name__, exc))
+    if not isinstance(state, dict) or state.get("format") != restart.FORMAT:
+        found = state.get("format") if isinstance(state, dict) else type(state).__name__
+        parser.error("{} is a state of format {!r}; this model resumes format {}".format(path, found, restart.FORMAT))
+    if not state.get("meta") or "args" not in state["meta"]:
+        parser.error("{} was saved without the command line's arguments (settings.save_meta), which a resume needs".format(path))
+    return state
+
+
+def cmd_resume(args, parser):
+    """Resume a saved flight (sub-plan 22 Task 5b; Asha, 2026-11-01): the state `--save-state-at` wrote, loaded as a new
+    object graph and flown on by `CoupledRun.run(resume=state)` from the end of the saved macro step, with the saved
+    random state (numpy's only draws are pyamg's, so it is part of the flight) unless --seed is given. Writes the files
+    `run` writes, from the saved flight's arguments with the options given here; the JSON's "resumed_from" records the
+    state and what was changed. Only the droplet records after the resume are written."""
+    state = load_state(args.state, parser)
+    run, saved = state["run"], argparse.Namespace(**state["meta"]["args"])
+    run_dir = os.path.join(args.outdir, args.name)
+    settings, changed = run.settings, {}
+    if args.t_max is not None:
+        if not (math.isfinite(args.t_max) and args.t_max > state["time_s"]):
+            parser.error("--t-max must be a finite time after the saved state's {:.6g} s, got {!r}".format(state["time_s"], args.t_max))
+        changed["--t-max"] = args.t_max
+    if args.dt_continuum is not None:
+        if settings.dt_continuum is None:
+            parser.error("--dt-continuum: the saved flight has no continuum switch (it did not run with --melt on and Girin removal)")
+        if not (math.isfinite(args.dt_continuum) and 0.0 < args.dt_continuum <= saved.dt):
+            parser.error("--dt-continuum must be a finite step in (0, the saved --dt = {} s], got {!r}".format(saved.dt, args.dt_continuum))
+        changed["--dt-continuum"] = args.dt_continuum
+    if args.frames_every is not None:
+        if args.frames_every < 0:
+            parser.error("--frames-every must be >= 0")
+        changed["--frames-every"] = args.frames_every
+    if args.slurry_spray is not None:
+        if args.slurry_spray == "on" and getattr(run.body, "spray", None) is None:
+            parser.error("--slurry-spray on needs a saved flight with the skins' spraying (--skin-flow on)")
+        changed["--slurry-spray"] = args.slurry_spray
+    if args.seed is not None:
+        changed["--seed"] = args.seed
+    # -- the point between loading the state and flying on, where each option is applied; keep the options above and
+    # their application below together, in the order of resume_arguments --
+    if "--t-max" in changed:
+        run.sim.settings.t_max = args.t_max
+    if "--dt-continuum" in changed:
+        settings.dt_continuum = args.dt_continuum
+        if run.switched_at is not None:            # latched: the step in use is the continuum step
+            settings.dt = args.dt_continuum
+    settings.save_at = settings.save_path = settings.save_meta = None
+    settings.output_dir = os.path.join(run_dir, "vtk")
+    settings.frames = []                           # the new series lists only its own frames
+    if "--frames-every" in changed:
+        settings.frames_every = args.frames_every
+        run._frame_dt, run._next_frame = args.frames_every * saved.dt, run.sim.t
+    if "--slurry-spray" in changed:
+        run.body.slurry_spray = args.slurry_spray == "on"
+    if isinstance(run.body.source_rows, spray.SourceSpool):      # only the records after the resume are written
+        run.body.source_rows = spray.SourceSpool(os.path.join(run_dir, SOURCE_SPOOL))
+    if args.seed is not None:
+        np.random.seed(args.seed)
+    else:
+        np.random.set_state(state["rng"])
+    history = run.run(resume=state)
+    # the arguments of the writer: the saved flight's, with what is changed here
+    for key, value in (("outdir", args.outdir), ("name", args.name), ("save_state_at", None)):
+        setattr(saved, key, value)
+    for key, flag in (("t_max", "--t-max"), ("dt_continuum", "--dt-continuum"), ("seed", "--seed"),
+                      ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray")):
+        if flag in changed:
+            setattr(saved, key, changed[flag])
+    atm, atm_name, atm_info = make_atmosphere(saved.atmosphere, saved.epoch, saved.wind)
+    reference = sesam_io.load_reference(saved.reference) if saved.reference else None
+    mass = body.sphere_mass(run.sim.settings.diameter, saved.material_density)
+    thermal_info = {key: getattr(saved, key) for key in ("heating", "melt", "surface_model", "thermal_solver", "linear_solver",
+                                                         "skin_flow", "slurry_spray")}
+    thermal_info.update({"macro_step_s": saved.dt, "frames_every": settings.frames_every, "dt_continuum_s": settings.dt_continuum})
+    resumed = {"state": os.path.abspath(args.state), "step": state["step"], "time_s": state["time_s"],
+               "run_name": state["meta"]["name"], "changed": changed}
+    os.makedirs(args.outdir, exist_ok=True)
+    return finish_run(saved, args.name, run.sim.settings, mass, atm_name, atm_info, reference, thermal_info, run.body, history,
+                      None, resumed=resumed)
+
+
 def _fmt(x):
     return "n/a" if x is None else "{:.1f}".format(x)
 
@@ -801,6 +926,8 @@
     try:
         if args.command == "run":
             return cmd_run(args, parser)
+        if args.command == "resume":
+            return cmd_resume(args, parser)
         return cmd_compare(args, parser)
     except KeyError as exc:
         # e.g. fap.solar_indices(): the run epoch (or the day before it) has no record in
```

- [ ] **Step 4: `flights/fly.py`** runs `resume` as the first of its arguments as the command line's resume, not `run`:

```python
        argv = args if args[:1] == ["resume"] else ["run"] + args
        return cli.main(argv + ["--outdir", outdir, "--name", name])
```

- [ ] **Step 5: Run the tests to see them pass.** Expected: 17 passed in about 2.5 minutes in the file; the unit tier 549
  passed, 1 skipped (Task 5a's 541 and these 8; 10.6 minutes across four processes).

- [ ] **Step 6: The production check** (in `measure/task5b/`, on mains power, through `flights/launch.sh`, then `compare.py`):
  resume `d050__physics__skin.state.pkl` to 172.0 s and `d100__physics__skin.state.pkl` to 23.0 s into
  `reentry_model_output/task5b_check`, and compare with the snapshot flights' `rows.csv` at full precision: no column
  but the `wall_` timers differs. Recorded in `measure/task5b/check.md`.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/cli.py tests/test_reentry_model_restart.py
git commit -m "Coupled runs: the resume command (python -m reentry_model resume STATE --outdir DIR --name NAME [--t-max --dt-continuum --seed --frames-every --slurry-spray]) continues a flight saved with --save-state-at to the last bit, writing the files run writes through one shared writer, the JSON naming the state and what was changed; exit 2 for a missing or foreign state, an end of flight not after the saved time, or a continuum step outside (0, --dt]; flights/fly.py runs it (sub-plan 22 Task 5b)"
```

---

### Task 5c: The deposit spread to lower corners, landing hot

**Files:**
- Create: `reentry_model/spreading.py`
- Modify: `reentry_model/skinbody.py` (`SPREAD_COLUMNS`, the class's defaults, `SkinBody.__init__`, `advance`, `_move_base`,
  `spread_stats`, `skin_stats`)
- Modify: `reentry_model/cli.py` (`parse_spread_limit`, `spread_limit_of`, `--spread-limit` on `run` and in `resume_arguments`, the
  name suffix, the settings' record, `cmd_run`'s refusal, `cmd_resume`'s option and its marked block)
- Modify: `reentry_model/coupled.py` (`skin_results`; the note after `SKINBODY_COLUMNS`)
- Test: `tests/test_reentry_model_spreading.py` (new, 10 tests), `tests/test_reentry_model_skinbody.py` (7),
  `tests/test_reentry_model_cli.py` (4 tests, 16 cases), `tests/test_reentry_model_restart.py` (2 tests, and 3 more cases of Task 5b's
  refusals)
- Create (git-ignored): `prototype/work-2026-10-28-skin-spraying/measure/task5c/`, the production check

**Interfaces:**
- Consumes: `Columns` (`edges`, `_tri`, `_opp`, `points`, `normal`, `area`, `faces`), `SkinBody._hand_cap()`, Task 5b's
  `resume_arguments` and the block of `cmd_resume` where options are applied.
- Produces: `spreading.fold_angles(columns, points)`, `(angle_deg, convex)` per edge; `spreading.spread_hand_down(columns, points,
  normal, area, rho, drawn, handed, enthalpy, cap, limit_deg, passes=SPREAD_PASSES)`, a `Spread(handed, enthalpy, moved, unplaced,
  over, landed)`, `landed` true at a column where metal arrived in any pass; `SkinBody(..., spread_limit=None)`; `skinbody.SPREAD_COLUMNS` and `SkinBody.spread_stats()`; `--spread-limit DEG|off` on
  `run` and `resume`, off by default, DEG within (0, 180), exit 2 otherwise or without a melting skin run; the name suffix
  `_spread-<deg>`; with a limit only, the history's `spread_mass_kg` (since time zero, each hop counted), `spread_unplaced_kg` and
  `spread_edges_over` (at the last move) and `spread_jump_max_K` (the largest landing over its node's temperature since the last
  row), the results' `spread_limit_deg` and `spread_mass_kg`, and the JSON's `settings.skin.spread_limit_deg`
  (`settings.spread_limit_deg` in a resume, whose block is Task 5b's short one).

Asha's idea, in answer to the two skin flights stopping: spread the re-frozen mass to neighbouring cells to keep the angles
reasonable (`docs/superpowers/specs/2026-11-06-deposit-spreading-design.md`). The snapshot flights put the fold at the toe of a ring of
re-frozen metal travelling outward (`measure/task5/snapshots/findings.md`), a ridge rather than a sawtooth, so she decided on
2026-11-06 to build it, and that **spread metal lands hot**. The lower neighbours of a toe are its cold corners: across folded edges the
3D temperature jumps by a median 62 K and 31 K, and about 90 % of the edges exceed the 5 K fault tolerance, so the tolerance does not
apply at a receiving corner. **The pass**, in `SkinBody._move_base` after `hand_down` and the `held` lines: where the base as it would
stand folds past the limit -- on the draws and the hand-down, the part's absorption and smoothing left out -- across a ridge along an
edge (a, b) the corners a and b give to the corners opposite it, c and e, and across a valley c and e give to a and b, each of the four
pairs wanting `SPREAD_SHARE` of the fold's excess over its angle of the metal in the donor's dual cell above the receiver. No donor
falls below what balances its draw and no receiver rises above the part's allowance, so the part's move ratio holds and a hand-down
still never splits a part. The metal keeps the donor's mean enthalpy, at most `SPREAD_PASSES` passes a part. **What crosses:** the
skins have already given up `handed`, so nothing in them changes; `moved`, `inward`, `h_outward` and the reconciliation read the
spread hand-down, `handed_s` and `landed_H`; `jump_faults`, `held_mass` and `handed_mass` keep what left each skin. **The landing:** at
a corner where metal arrived in any pass (`Spread.landed`), whatever it gave onward, the metal crosses the base at the node's enthalpy
and the rest of its heat, the excess, goes to the node through `pending`, after the reconciliation, because `_advance_part` has set `pending` by then and a part split or repeated has
returned or restored before it. `MovingMesh.move` shares a correction over all the moving nodes, so booked through `h_outward` the heat
would go to the whole zone rather than that corner; `pending` is the deferred per-corner heat the next part's 3D solve takes as a load
and the books count. **A limit that does not bind** returns the very arrays it was given. **A state saved before this** has none of the
attributes on its body, so they are class defaults too; and a run's columns are its first row's, the saved row's in a resume, so
`cmd_resume`'s block drops the spreading columns from it and adds `spread_stats()`, which is empty without a limit.

**Measured** (2026-11-06; `measure/task5c/check.md`, each once, on the 50 mm sphere, one limit, one mesh; on the code just before ruling 11,
which changes where some landing heat goes and was not flown). On the test plate a crest handing down 0.025 of an edge a part for 200 parts,
moved as the body moves its corners, passes 90 degrees by part 31 with the pass off and stays between 59.96 (over the last 100 parts) and
60.03 degrees with a limit of 60. On the production mesh, resumed from the state saved at 170.5 s (pickled by 055c313, before spreading
existed): to 172.0 s, in 193 s, all 84 shared columns equal the snapshot flight's at full precision, the limit binding nowhere at a 37.22
degree fold. To 177.0 s, in 2 014 s, the limit first binds in the row at 174.0 s and the largest fold over the 91 rows is 64.86 degrees, at
176.0 s, against 116.75 over the same flight time in the snapshot flight; 82 rows are above 60 degrees, 4 above 62, none above 65, and the last
is 61.82. At 177.0 s 5 base edges are past 60 degrees against 108, in the same travelling band. 1.704 g moved against 13.63 g handed down (12.5 %,
each hop counted), the largest unplaced mass at a move 0.25 mg, the largest landing recorded 57.4 K above its node. Mass closes to 4.5e-16, the
energy balance residual is -1.76e-11, the momentum residual 3.76e-17, and no fault is counted against the snapshot's 2. **After 176.0 s the two
flights are not like with like:** the continuum switch latched at 176.0 s here and 176.5 s in the snapshot flight, so from there they ran different
macro steps, 0.0125 s against 0.5 s, and different time grids. Up to 176.0 s, the same 11 macro steps and 95 parts, the largest fold is 64.86
against 114.59, the lowest 3D radius ratio, which is what stops a flight, 0.0756 against 0.0678 (the floor is 0.05), rebuilds 6 against 10, the
deepest skin 2.22 mm against 3.04 mm. Later figures carry both differences: to 177.0 s the lowest ratio is 0.0533 against 0.0503, near the floor
either way, and at 177.0 s the deepest skin is 5.99 mm against 4.48 mm and 13.63 g is handed down against 15.69 g, and how much of that is the
spreading is unknown. The earlier switch is a threshold: the first liquid in both is one corner at the nose, 0 degrees from the point of highest
heating and far from the folds, at 908.11 K against 907.98 K at 176.0 s, the liquidus being 908.0 K; why spreading moves the nose by 0.12 K is not
traced. The check stops 2.04 s before the snapshot flight's stop, so whether a flight now gets past it is untested. The unit tier is 587 passed, 1 skipped (Task 5b's 549 and these 38; 15.1 minutes across four processes).

**Rulings** (2026-11-06; what, why, cost if wrong):
1. *`SPREAD_SHARE` stays at the 0.25 the design proposed; no tuning was needed.* From 0.0625 to 1.0 the test's largest fold is 60.03 to
   61.02 degrees, 0.25 the closest, and the passes matter more: one a part overshoots to 70.6, two to 60.4, four to 60.03. Cost: none found
   on the plate; the production check flew one value.
2. *The growing ridge's test moves its corners as the body does*, along their normals with the volume-neutral smoothing, irregular
   triangles, the plate's open rim held still and the crest stopping short of it, at a tenth of the allowance a part. Why: moved along tilted
   normals with no smoothing the triangles overturn with or without the pass, which tests nothing, and on an open plate the smoothing draws
   in the rim, which a closed surface does not have. A straight outward move would have passed too, but that is not how the body moves. Cost:
   the test repeats `MovingMesh.neutral_smoothing`'s few lines, which could drift. At twice the rate the fold holds at 60.7 for 140 parts, then
   a valley grows at the foot of the ridge's end, where only 4 of the 16 corners around it hold handed metal: a limit of the design, recorded
   in the test.
3. *A donor that gives all it has ends at exactly zero* (`SPENT`, 1e-14 of what it held). Why: round-off left there would make its mean
   enthalpy round-off over round-off, and `MovingMesh.move` books any sweep there at it; a sweep of 1e-13 m3 at a nonsense enthalpy of 3e10 J/kg
   would move about 8 J. Cost: at most 1e-14 of that donor's hand-down from the sum, within round-off.
4. *The landing's excess may be negative*: metal colder than its node lands at the node's enthalpy too, and `pending` takes the difference
   off it. Why: the same books and the same place for the heat. Cost: none found; the check's largest landing is hot.
5. *`spread_jump_max_K` is reset at the start of `advance`*, as `step_spray` is, not when read: reading the stats twice then gives the same
   answer. Cost: none, a row is made once a macro step.
6. *`--spread-limit` is accepted on the heat-only skins (`--skin-flow off`) too*: they hand metal down the same way. Cost: that case is
   unflown; only its argument is tested.
7. *The settings record the limit only when there is one*, so every run's JSON without the option is unchanged. *On `resume`, `off` asks
   nothing of the flight, and a limit is refused unless the body is melting skins.* Cost: none.
8. *The 177.0 s check used `--frames-every 1`, not 40.* After the continuum switch frames go by flight time, at n times the saved 0.5 s step,
   so 40 would give one frame in 20 s; 1 gives 13 frames, 58 MB, at the snapshots' cadence. Cost: none.
9. *`spread_mass_kg` sums the passes*, so metal that moves twice in a part counts twice. Why: it is what the passes moved, and simple. Cost:
   the column overstates the net sideways mass; the 12.5 % above is an upper bound on it.
10. *At a corner that received, all its hand-down lands at the node's enthalpy, its own metal too*, not only what came from a neighbour. Why:
    `h_outward` is one number a corner, and its own metal is within the fault tolerance of the node anyway. Cost: at most 5 K of that metal's
    heat goes through `pending` rather than across the base.
11. *Where metal landed is read pass by pass, `Spread.landed`, not from `handed_s > handed`* (Task 5c's review; the coordinator's ruling). Why: a
    corner that takes hot metal and gives more onward ends handing down less than it started with, and its neighbours' heat then crossed through
    `h_outward`, over the whole moving zone, the route Asha's decision excluded, and was left out of `spread_jump_max_K`. The review measured, on a
    jittered plate, 16 to 127 corners missed against 174 to 273 caught, carrying 1.9 to 7.0 % as much heat as went through `pending`, and in 4 of 6
    cases a missed landing hotter than the largest recorded (38.0 K against 18.8 K at the worst); on the 4 mm melting test body the recorded largest
    goes from 1.18 K to 1.42 K with the fix. The books were closed either way. A corner in `landed` that passed on all it took has no metal and no
    enthalpy left, so its excess is zero and it records no landing. Cost: the fix was not flown, so the production numbers above understate the
    largest landing by an unknown amount and put some of the heat in the zone.

- [ ] **Step 1: The tests** -- `tests/test_reentry_model_spreading.py`:

```python
"""spreading.py: the metal being handed down spread to lower neighbouring corners where the base would fold past a limit
(sub-plan 22 Task 5c; Asha 2026-11-06), on a small plate of near-equilateral 1 mm triangles, six neighbours to a corner as
on the production surface."""
import numpy as np
import pytest

from reentry_model import columns, movemesh, spreading

RHO = 2813.0                    # kg/m3
EDGE = 1e-3                     # m: the plate's triangles' edge
SPECIFIC = 5e5                  # J/kg: the enthalpy of the metal handed down, unless a test says otherwise


def plate(nx=25, ny=25, edge=EDGE, jitter=0.0):
    """A flat plate in z = 0: `ny` rows of `nx` corners, every other row shifted half an edge and each corner up to `jitter`
    edges off its place in x and y, as a `columns.Columns` (its triangles outward, +z), with `row(j)`, the corners of row j."""
    x, y = np.meshgrid(np.arange(nx) * edge, np.arange(ny) * edge * np.sqrt(3.0) / 2.0, indexing="ij")
    x = x + 0.5 * edge * (np.arange(ny) % 2)[None, :]
    points = np.c_[x.ravel(), y.ravel(), np.zeros(nx * ny)]
    points[:, :2] += jitter * edge * np.random.default_rng(1).uniform(-1.0, 1.0, (nx * ny, 2))
    at = lambda i, j: i * ny + j
    faces = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            a, b, c, d = at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1)
            faces += [[a, b, d], [a, d, c]] if j % 2 == 0 else [[a, b, c], [b, d, c]]
    col = columns.Columns(points, np.array(faces))
    col.row = lambda j: np.array([at(i, j) for i in range(nx)])
    return col


def raised(col, rows, heights):
    """The plate with each of the rows in `rows` lifted to its height in `heights` [m]."""
    p = col.points.copy()
    for j, z in zip(rows, heights):
        p[col.row(j), 2] += z
    col.update(p)
    return col


def spread(col, handed, limit=60.0, drawn=None, cap=None, enthalpy=None, passes=spreading.SPREAD_PASSES):
    drawn = np.zeros(col.n) if drawn is None else drawn
    allowance = 0.25 * col.edge_length.mean() * RHO * col.area
    H = SPECIFIC * handed if enthalpy is None else enthalpy
    return spreading.spread_hand_down(col, col.points, col.normal, col.area, RHO, drawn, handed, H,
                                      allowance if cap is None else cap, limit, passes=passes)


def test_the_fold_is_the_surface_s_measure_with_a_ridge_convex_and_a_valley_not():
    """`fold_angles` reads the same angle the moving mesh's `surface_fold` reports, and says which way the edge folds; an
    edge on an open boundary has one triangle and folds nowhere."""
    col = raised(plate(), [12], [0.5 * EDGE])
    angle, convex = spreading.fold_angles(col, col.points)
    assert angle.max() == pytest.approx(movemesh.max_dihedral(col.points, col.faces), rel=1e-12)
    crest, toe = np.isin(col.edges, col.row(12)).all(axis=1), np.isin(col.edges, col.row(11)).all(axis=1)
    assert angle[crest][2:-2].min() == pytest.approx(60.0, abs=1e-9) and convex[crest].all()     # two 30 degree faces
    assert angle[toe][2:-2].min() == pytest.approx(30.0, abs=1e-9) and not convex[toe].any()
    assert (angle[col._tri[:, 1] < 0] == 0.0).all() and not convex[col._tri[:, 1] < 0].any()


def test_a_limit_that_does_not_bind_returns_the_hand_down_itself():
    """On a flat plate nothing folds, so the hand-down and its enthalpy come back as the very arrays given and nothing has
    moved: a limit that does not bind changes nothing to the last bit."""
    col = plate()
    handed = np.full(col.n, 1e-7)
    H = SPECIFIC * handed
    s = spread(col, handed, enthalpy=H)
    assert s.handed is handed and s.enthalpy is H
    assert (s.moved, s.unplaced, s.over) == (0.0, 0.0, 0) and not s.landed.any() and s.landed.shape == (col.n,)


def test_a_ridge_hands_its_metal_to_its_lower_neighbours_and_keeps_the_sums():
    """A ridge 1 mm high along a row of corners, folded about 98 degrees against a 60 degree limit, handing down a quarter
    of the move allowance, with its own enthalpy per column: the sums of the mass and the enthalpy are kept to round-off,
    every crest corner hands down less than it would have, the metal lands on corners lower than the crest and on no other,
    and what lands first on a row that handed nothing carries the mean enthalpy of the crest's metal."""
    col = raised(plate(), [12], [1e-3])
    crest = col.row(12)
    handed = np.zeros(col.n)
    handed[crest] = 0.25 * 0.25 * EDGE * RHO * col.area[crest]
    unit = np.full(col.n, SPECIFIC)
    unit[crest] = SPECIFIC * (1.0 + 0.01 * np.arange(len(crest)) / len(crest))     # no two columns' metal alike
    H = unit * handed
    s = spread(col, handed, enthalpy=H)
    assert s.handed is not handed and s.over > 0 and s.moved > 0.0
    assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12) and s.enthalpy.sum() == pytest.approx(H.sum(), rel=1e-12)
    inner = crest[2:-2]                                            # clear of the plate's open ends
    assert (s.handed[inner] < handed[inner]).all()
    rose = s.handed > handed
    assert rose.any() and (col.points[rose, 2] < col.points[crest, 2].max()).all()
    assert (s.handed[col.row(11)[2:-2]] > 0.0).all() and (s.handed[col.row(13)[2:-2]] > 0.0).all()
    assert (s.handed >= 0.0).all()
    got = s.enthalpy[rose] / s.handed[rose]
    assert got.min() >= SPECIFIC * (1.0 - 1e-12) and got.max() <= unit[crest].max() * (1.0 + 1e-12)
    assert s.moved <= spreading.SPREAD_PASSES * handed.sum()


def test_a_valley_fold_gives_from_the_rim_to_the_floor():
    """A trench along a row, its rims 0.8 mm above its floor: the fold is the floor's, a valley of about 95 degrees (the rims
    fold about half that, within the 70 degree limit), so the rims hand their metal to the floor and the floor's corners,
    which handed down nothing, receive it."""
    col = plate()
    p = col.points.copy()
    p[:, 2] = 0.8e-3
    p[col.row(12), 2] = 0.0
    col.update(p)
    angle, convex = spreading.fold_angles(col, col.points)
    assert not convex[angle > 70.0].any() and (angle > 70.0).sum() > 0
    floor, rims = col.row(12), np.r_[col.row(11), col.row(13)]
    handed = np.zeros(col.n)
    handed[rims] = 0.25 * 0.25 * EDGE * RHO * col.area[rims]
    s = spread(col, handed, limit=70.0)
    assert (s.handed[floor[2:-2]] > 0.0).all() and handed[floor].sum() == 0.0
    assert s.handed[rims].sum() < handed[rims].sum()
    assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12)
    rose = np.flatnonzero(s.handed > handed)
    assert set(rose) <= set(floor)


def test_a_corner_that_receives_and_gives_more_onward_is_marked_as_landed():
    """A two-step cascade down a staircase whose 49 degree toe is past a 40 degree limit: the crest gives its hotter metal to the
    flank below it, and the flank, which holds far more, gives more onward than it took, so it ends handing down less than it
    started with and its total alone would not show that anything arrived. It carries the crest's heat all the same, so `landed`
    marks it, and the toe it gave to; nothing arrives at the crest (Finding 1 of Task 5c's review: comparing the totals missed
    these corners, and their neighbours' heat crossed the base through `h_outward`)."""
    col = raised(plate(), [10, 11, 12, 13, 14], [0.0, 1e-3, 2e-3, 1e-3, 0.0])
    crest, flank, toe = col.row(12)[2:-2], col.row(11)[2:-2], col.row(10)[2:-2]
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    handed = np.zeros(col.n)
    handed[col.row(12)], handed[col.row(11)], handed[col.row(13)] = (0.1 * cap[col.row(12)], 0.8 * cap[col.row(11)],
                                                                     0.8 * cap[col.row(13)])
    unit = np.full(col.n, SPECIFIC)
    unit[col.row(12)] = 1.1 * SPECIFIC
    s = spread(col, handed, limit=40.0, cap=cap, enthalpy=unit * handed)
    assert (s.handed[flank] < handed[flank]).all() and (s.handed[flank] > 0.0).all()      # net donors, with metal left
    assert s.landed[flank].all() and s.landed[toe].all() and not s.landed[crest].any()
    assert s.landed.dtype == bool and s.landed.shape == (col.n,)
    mean = s.enthalpy[flank] / s.handed[flank]
    assert (mean > SPECIFIC * (1.0 + 1e-6)).all() and (mean < 1.1 * SPECIFIC).all()       # the crest's heat came with it
    assert (s.landed[s.handed > handed]).all()                                            # whatever rose, landed


def growing_ridge(limit, parts=200, nx=37, ny=37):
    """A jittered plate on which the middle of the middle row, 21 corners stopping short of the rim, hands down a tenth of
    the move allowance every part: the spreading pass (with `limit`, or none), then the body's move, `SkinBody._move_base`'s
    -- every corner out along its normal by what lands on it, with the tangential smoothing made volume-neutral as
    `MovingMesh.neutral_smoothing` makes it -- the plate's open rim held still, as a closed surface has no rim for the
    smoothing to draw in. Returns the largest fold [deg] after every part and the plate."""
    col = plate(nx, ny, jitter=0.1)
    crest = col.row(ny // 2)[8:-8]
    rim = np.unique(col.edges[col._tri[:, 1] < 0].ravel())
    edges = movemesh.surface_edges(col.faces)

    def held(u):
        u[rim] = 0.0
        return u
    folds = []
    for _ in range(parts):
        p, normal, area = col.points, col.normal, col.area
        allowance = 0.25 * col.edge_length.mean() * RHO * area
        handed = np.zeros(col.n)
        handed[crest] = 0.1 * 0.25 * EDGE * RHO * area[crest]
        landing = handed
        if limit is not None:
            s = spreading.spread_hand_down(col, p, normal, area, RHO, np.zeros(col.n), handed, SPECIFIC * handed, allowance,
                                           limit)
            assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12) and (s.handed <= allowance * (1.0 + 1e-12)).all()
            landing = s.handed
        u = held(movemesh.tangential_smoothing(p, col.faces, normal, edges=edges))
        for _ in range(movemesh.NEUTRAL_PASSES):
            u = held(u + normal * (movemesh.swept_volumes(p, p + u, col.faces, col.faces, col.n) / area)[:, None])
        col.update(p + u + held(normal * (landing / (RHO * area))[:, None]))
        folds.append(float(spreading.fold_angles(col, col.points)[0].max()))
    return np.array(folds), col


def test_repeated_parts_hold_a_growing_ridge_near_the_limit():
    """The test of the pass (Asha's idea, 2026-11-06): a crest handing down 0.025 of an edge every part for 200 parts, moved
    as the body moves its corners. With the pass off it passes 90 degrees, past spec 12's 75-95 degrees, by part 31 and
    goes on climbing; with a 60 degree limit the base stays within 5 degrees of the limit the whole way, not flattened below
    it either, the ridge widening instead. Measured 2026-11-06 on this plate: 60.03 degrees at most and 59.96 at least over
    the last 100 parts. SPREAD_SHARE stayed at the 0.25 the design proposed: from 0.0625 to 1.0 the largest fold is 60.03 to
    61.02, 0.25 the closest, so the share is not what holds the fold; a limit of 45, 75 or 90 degrees holds to within 0.5; one
    pass a part overshoots to 70.6, two to 60.4 and four, SPREAD_PASSES, to 60.03. At twice this crest's rate the fold holds
    at 60.7 for 140 parts, then a valley grows at the foot of the ridge's end, where only 4 of the 16 corners around it hold
    any handed metal and the smoothing has squeezed the triangles to 0.2 mm edges: only handed metal moves."""
    off = growing_ridge(None, parts=60)[0]
    assert off[:40].max() > 90.0 and off.max() > 120.0
    folds, col = growing_ridge(60.0)
    assert folds.max() <= 65.0
    assert folds[-100:].min() >= 55.0
    assert (col.points[:, 2] > 0.1e-3).sum() > 3 * 21             # it has widened: more than three rows of the crest's length


def test_a_receiving_corner_stays_within_its_move_allowance():
    """However much the ridge would give, no corner hands down more than a part's allowance: a corner's own hand-down plus
    what lands on it moves its base by at most the allowance, the cap of 2026-10-27, so a part is never split for it."""
    col = raised(plate(), [12], [2e-3])
    crest = col.row(12)
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    handed = np.zeros(col.n)
    handed[crest] = cap[crest]
    handed[col.row(11)] = 0.9 * cap[col.row(11)]                   # next to nothing left of the allowance
    s = spread(col, handed, cap=cap)
    assert s.moved > 0.0
    assert (s.handed <= cap * (1.0 + 1e-12)).all()
    assert s.handed[col.row(11)].max() > 0.9 * cap[col.row(11)].max()


def test_a_donor_keeps_what_balances_its_draw():
    """A corner drawing metal moves in by its draw less its hand-down, and that move too is within the allowance, or the part
    is split for it: a donor gives only down to its draw less the allowance, and a corner already at or under that gives
    nothing."""
    col = raised(plate(), [12], [2e-3])
    crest = col.row(12)
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    drawn, handed = np.zeros(col.n), np.zeros(col.n)
    drawn[crest], handed[crest] = 1.6 * cap[crest], cap[crest]     # moves in by 0.6 allowances, and may reach 1
    idle = crest[10]
    drawn[idle], handed[idle] = 2.0 * cap[idle], 0.5 * cap[idle]   # over the allowance already: it gives nothing
    s = spread(col, handed, drawn=drawn, cap=cap)
    floor = drawn - cap
    assert s.moved > 0.0
    others = crest[crest != idle]
    assert (s.handed[others] >= floor[others] * (1.0 - 1e-12)).all()
    assert (s.handed[others[2:-2]] < handed[others[2:-2]]).any()
    assert s.handed[idle] >= handed[idle]


def test_what_cannot_be_placed_stays_and_is_counted():
    """Lower corners with no room left: with every neighbour of the crest at its allowance nothing can land, so the crest
    keeps its metal and the mass that wanted to move is reported, as `held_mass` reports what the tolerance holds. With a
    little room, that much lands and the rest is reported."""
    col = raised(plate(), [12], [2e-3])
    crest, flanks = col.row(12), np.r_[col.row(11), col.row(13)]
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    handed = np.zeros(col.n)
    handed[crest], handed[flanks] = 0.5 * cap[crest], cap[flanks]
    # the flanks' toes are folded too and may give to the next rows, leaving room: a limit only the crest exceeds stops that
    q = col.points + col.normal * (handed / (RHO * col.area))[:, None]
    limit = float(spreading.fold_angles(col, q)[0].max()) - 5.0
    s = spread(col, handed, cap=cap, limit=limit, passes=1)
    assert s.moved == 0.0 and s.unplaced > 0.0 and s.over > 0
    np.testing.assert_array_equal(s.handed, handed)
    handed[flanks] = 0.9 * cap[flanks]
    room = (cap - handed)[flanks].sum()
    s = spread(col, handed, cap=cap, limit=limit, passes=1)
    assert 0.0 < s.moved <= room * (1.0 + 1e-12) and s.unplaced > 0.0


def test_a_donor_that_gives_all_it_has_hands_down_exactly_nothing():
    """A corner with a sliver of a hand-down on a ridge far past the limit gives it all: its hand-down and its enthalpy are
    exactly zero, not round-off, whose quotient would be the mean enthalpy the moving mesh books its sweep at, and the sums
    still hold."""
    col = raised(plate(), [12], [1e-3])
    crest = col.row(12)
    handed = np.zeros(col.n)
    handed[crest] = 1e-3 * 0.25 * EDGE * RHO * col.area[crest]
    H = SPECIFIC * handed
    s = spread(col, handed, enthalpy=H)
    inner = crest[2:-2]
    assert (s.handed[inner] == 0.0).all() and (s.enthalpy[inner] == 0.0).all() and not s.landed[inner].any()
    assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12) and s.enthalpy.sum() == pytest.approx(H.sum(), rel=1e-12)
```

  in `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1381,3 +1381,215 @@
     assert got[1] == pytest.approx(2.0 * sk.thickness()[1], rel=1e-14)
     assert np.allclose(dz[2], dz[2, 0]) and got[2] == pytest.approx(10.7 * dz[2, 0], rel=1e-12)
     assert np.isnan(got[3:]).all()
+
+
+# -- spreading the deposit (sub-plan 22 Task 5c; Asha 2026-11-06) ----------------------------------------------------------
+def melting_body(mesh, **kw):
+    """The melting test body of the spraying tests above, 890 K under the 100 mm flight's gas state at 93 s, one flow
+    sub-step a heat sub-step, flown for two 0.125 s macro steps at 3 MW/m2, numpy's generator seeded as a run seeds it: it
+    hands down 3.8 g. Returns the body and its starting mass."""
+    pytest.importorskip("cantera")
+    from test_reentry_model_surface_flow import aero_state_at
+    from reentry_model import spray
+    np.random.seed(12345)
+    mat = material.Material.from_drama_json("AA7075_scheil")
+    b = make(mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid), **kw)
+    M0 = b.mass()
+    state = aero_state_at(93.0)
+    for k in range(2):
+        b.advance(0.125 * (k + 1), 0.125, front_loads(b, 3e6), state=state)
+    return b, M0
+
+
+def zero_book(sk):
+    """What no draw left: the base's operations alone, as the deposit tests above run them."""
+    return skin.SkinBook(*(np.zeros(sk.n) for _ in range(8)), momentum_out=np.zeros((sk.n, 3)))
+
+
+def bump(b, peak, share=0.9, T=700.0):
+    """A deposit at `T` on the bottoms of the columns `peak`, `share` of a part's allowance deep."""
+    mat, sk, col = b.material, b.skins, b.columns
+    extra = np.zeros(sk.n)
+    extra[peak] = share * skinbody.MOVE_MAX * col.edge_length.mean() * mat.rho * col.area[peak]
+    sk.add_bottom(extra, float(mat.enthalpy(T)))
+    sk.regrid()
+
+
+def spied(b, monkeypatch):
+    """The arguments of the base move's spreading and of its corner move, as they come: ([(hand-down given, Spread)], [the
+    enthalpy each outward sweep is to be booked at])."""
+    spreads, outward = [], []
+    spread, policy = skinbody.spreading.spread_hand_down, b.moving.move_with_policy
+
+    def recorded(*a, **k):
+        spreads.append((a[6], spread(*a, **k)))
+        return spreads[-1][1]
+
+    def recorded_move(u, h_inward=None, h_outward=None):
+        outward.append(h_outward)
+        return policy(u, h_inward, h_outward)
+    monkeypatch.setattr(skinbody.spreading, "spread_hand_down", recorded)
+    b.moving.move_with_policy = recorded_move
+    return spreads, outward
+
+
+def neighbours(b, corners):
+    """The corners sharing an edge with any of `corners`, those left out."""
+    e = b.columns.edges
+    near = np.unique(np.r_[e[np.isin(e[:, 0], corners), 1], e[np.isin(e[:, 1], corners), 0]])
+    return near[~np.isin(near, corners)]
+
+
+def test_a_limit_that_does_not_bind_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
+    """Nothing folds near 179.9 degrees, so the melting body flown with that limit is the body flown without one in every
+    array and every shared history column, and only the four spreading columns are added, recording nothing moved. Without
+    a limit there are no spreading columns: every run without the option keeps its column set."""
+    off, _ = melting_body(coarse_sphere_mesh)
+    on, _ = melting_body(coarse_sphere_mesh, spread_limit=179.9)
+    assert off.handed_mass > 0.0 and on.spread_limit == 179.9 and off.spread_limit is None
+    assert same_body(off, on)
+    a, c = off.skin_stats(), on.skin_stats()
+    same = lambda x, y: x == y or (np.isnan(x) and np.isnan(y))
+    assert all(same(a[k], c[k]) for k in a if not k.startswith("wall_"))
+    assert set(c) - set(a) == set(skinbody.SPREAD_COLUMNS) and not set(a) & set(skinbody.SPREAD_COLUMNS)
+    assert [c[k] for k in skinbody.SPREAD_COLUMNS] == [0.0, 0.0, 0, 0.0]
+
+
+def test_a_limit_outside_a_fold_s_range_is_refused(coarse_sphere_mesh):
+    for bad in (0.0, 180.0, -1.0, float("nan")):
+        with pytest.raises(ValueError):
+            make(coarse_sphere_mesh, spread_limit=bad)
+
+
+def test_spreading_keeps_the_books(coarse_sphere_mesh):
+    """The melting body with a limit low enough to bind on the 4 mm mesh, whose sphere folds 14 degrees before anything happens:
+    0.6 g of the 3.8 g handed down moves between corners, and the mass closes to 1e-12, the energy balance residual is below
+    1e-8 and the momentum books close, the hot landing's heat in transit counted (measured: mass exact, energy 1.7e-11,
+    momentum 1.4e-17)."""
+    b, M0 = melting_body(coarse_sphere_mesh, spread_limit=6.0)
+    stats = b.skin_stats()
+    assert stats["spread_mass_kg"] == b.spread_mass > 0.0 and stats["spread_edges_over"] > 0
+    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-12)
+    assert abs(b.energy_balance_residual()) < 1e-8
+    assert stats["momentum_residual"] < 1e-12
+
+
+def test_metal_landing_hot_takes_its_heat_to_the_receiving_corner_and_counts_no_fault(coarse_sphere_mesh, monkeypatch):
+    """Asha's decision of 2026-11-06: spread metal lands hot, because the lower neighbours of a toe are the cold corners ahead
+    of it. A proud corner of the 4 mm sphere at 700 K with its neighbours' nodes 40 K colder, past the 5 K fault tolerance:
+    the metal crosses each receiving base at that node's enthalpy, so the moving mesh shares no correction over the zone for
+    it, and the rest of its heat, the excess, waits in `pending` at that corner alone, which the books count and the next
+    part's 3D solve takes as a load there. No fault is counted, the largest landing is recorded, and the energy is exact
+    before and after the 3D model takes the heat."""
+    b = make(coarse_sphere_mesh, T0=700.0, spread_limit=10.0)
+    b.emissivity = 0.0
+    mat, sk, col, nodes = b.material, b.skins, b.columns, b.moving.corners
+    peak = np.argsort(b.theta)[:1]
+    ring = neighbours(b, peak)
+    T = b.solver.temperature().copy()
+    T[nodes[ring]] -= 40.0
+    b.solver.set_temperature(T)
+    bump(b, peak)
+    spreads, outward = spied(b, monkeypatch)
+    E0, T_node = b.energy(), T[nodes]
+    b._move_base(zero_book(sk), skin.minimum_mass(mat, sk.area, DEV), b._normal_gradient())
+    handed, s = spreads[0]
+    rose = s.landed
+    cold = rose & np.isin(np.arange(sk.n), ring)
+    assert cold.any() and s.moved > 0.0
+    excess = s.enthalpy[rose] - s.handed[rose] * mat.enthalpy(T_node[rose])
+    assert excess[np.isin(np.flatnonzero(rose), ring)].min() > 0.0
+    np.testing.assert_allclose(b.pending[rose], excess, rtol=1e-12)
+    assert (b.pending[~rose] == 0.0).all()
+    np.testing.assert_allclose(outward[0][cold], mat.enthalpy(T_node[cold]), rtol=1e-12)    # at the node's own enthalpy
+    assert b.jump_faults == 0
+    stats = b.skin_stats()
+    assert 30.0 < stats["spread_jump_max_K"] <= 40.0 + skin.JUMP_TOLERANCE
+    assert b.energy() == pytest.approx(E0, rel=1e-12)
+    T_before = b.solver.temperature()[nodes].copy()
+    zero = heating.HeatingResult(np.zeros(col.n), 0.0, 0.0, 0.0, 0.0)
+    b.advance(0.01, 0.01, zero)
+    assert b.energy() == pytest.approx(E0, rel=1e-12)
+    assert (b.solver.temperature()[nodes[cold]] > T_before[cold]).all()
+    assert b.skin_stats()["spread_jump_max_K"] < stats["spread_jump_max_K"]               # since the last row
+
+
+def test_a_corner_that_takes_hot_metal_and_gives_more_onward_still_sends_its_landing_heat_to_its_node(coarse_sphere_mesh,
+                                                                                                       monkeypatch):
+    """A two-step cascade on the 4 mm sphere at 700 K (Finding 1 of Task 5c's review): a peak with its node and metal 60 K
+    hotter stands above a ring of corners, themselves above the rest. The ring takes the peak's hot metal and gives more onward than
+    it took, so each of its corners that keeps any ends handing down less than it started with. The heat it took is still its
+    landing, as Asha decided: its metal crosses at the node's enthalpy, the rest goes to that node through `pending`, not through
+    `h_outward` and so over the whole moving zone, and the landing counts in the largest recorded. Energy closes and no fault is
+    counted."""
+    b = make(coarse_sphere_mesh, T0=700.0, spread_limit=15.0)
+    b.emissivity = 0.0
+    mat, sk, col, nodes = b.material, b.skins, b.columns, b.moving.corners
+    peak = np.argsort(b.theta)[:1]
+    ring = neighbours(b, peak)
+    u = np.zeros((col.n, 3))
+    u[peak], u[ring] = 2.4e-3 * col.normal[peak], 1.2e-3 * col.normal[ring]
+    b.moving.move_with_policy(u)
+    b._refresh_geometry()
+    nodes = b.moving.corners
+    T = b.solver.temperature().copy()
+    T[nodes[peak]] += 60.0
+    b.solver.set_temperature(T)
+    bump(b, peak, share=0.9, T=760.0)
+    bump(b, ring, share=0.5)
+    spreads, outward = spied(b, monkeypatch)
+    E0, T_node = b.energy(), b.solver.temperature()[nodes]
+    b._move_base(zero_book(sk), skin.minimum_mass(mat, sk.area, DEV), b._normal_gradient())
+    handed, s = spreads[0]
+    on = ring[(s.handed[ring] < handed[ring]) & (s.handed[ring] > 0.0)]
+    assert on.size > 0
+    jump = mat.temperature_from_enthalpy(s.enthalpy[on] / s.handed[on]) - T_node[on]
+    assert jump.max() > skin.JUMP_TOLERANCE
+    excess = s.enthalpy[on] - s.handed[on] * mat.enthalpy(T_node[on])
+    np.testing.assert_allclose(b.pending[on], excess, rtol=1e-12)
+    np.testing.assert_allclose(outward[0][on], mat.enthalpy(T_node[on]), rtol=1e-12)
+    assert b.skin_stats()["spread_jump_max_K"] >= jump.max() and b.jump_faults == 0
+    assert b.energy() == pytest.approx(E0, rel=1e-12)
+
+
+def test_a_donor_that_gives_away_all_its_hand_down_keeps_the_books(coarse_sphere_mesh, monkeypatch):
+    """A corner 1.5 mm proud of its neighbours with little to hand down -- 0.2 of a part's allowance, some of which goes to the
+    minimum its larger dual cell now needs -- gives all it has and nothing lands on it, being the highest: its hand-down is exactly nothing, so the moving mesh has no enthalpy for it (`h_outward` nan) and books any sweep
+    there at the corner's own, and the reconciliation and the fault count hold: mass and energy exact, nothing not finite, and
+    no fault."""
+    b = make(coarse_sphere_mesh, T0=700.0, spread_limit=15.0)
+    b.emissivity = 0.0
+    mat, sk, col = b.material, b.skins, b.columns
+    peak = np.argsort(b.theta)[:1]
+    u = np.zeros((col.n, 3))
+    u[peak] = 1.5e-3 * col.normal[peak]
+    b.moving.move_with_policy(u)
+    b._refresh_geometry()
+    bump(b, peak, share=0.2)
+    spreads, outward = spied(b, monkeypatch)
+    M0, E0 = b.mass(), b.energy()
+    b._move_base(zero_book(sk), skin.minimum_mass(mat, sk.area, DEV), b._normal_gradient())
+    handed, s = spreads[0]
+    assert handed[peak] > 0.0 and s.handed[peak] == 0.0
+    gone = (handed > 0.0) & (s.handed == 0.0)
+    assert gone.any() and (s.enthalpy[gone] == 0.0).all() and np.isnan(outward[0][gone]).all()
+    assert b.mass() == pytest.approx(M0, rel=1e-12) and b.energy() == pytest.approx(E0, rel=1e-12)
+    assert np.isfinite(b.solver.temperature()).all() and np.isfinite(sk.T[sk.valid]).all() and np.isfinite(b.pending).all()
+    assert b.jump_faults == 0
+
+
+def test_a_body_saved_without_spreading_attributes_flies_on(coarse_sphere_mesh):
+    """The two real states were pickled before spreading existed, so a body from them has none of its attributes on the
+    instance: the class's defaults give them, it flies on as a body without a limit, and it takes a limit a resume gives it."""
+    b = make(coarse_sphere_mesh, T0=700.0)
+    b.emissivity = 0.0
+    for key in [k for k in vars(b) if "spread" in k]:
+        del vars(b)[key]
+    assert not [k for k in vars(b) if "spread" in k]
+    assert b.spread_limit is None
+    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    b.advance(0.5, 0.5, zero)
+    assert not set(b.skin_stats()) & set(skinbody.SPREAD_COLUMNS) and b.spread_stats() == {}
+    b.spread_limit = 60.0
+    b.advance(1.0, 0.5, zero)
+    assert set(b.spread_stats()) == set(skinbody.SPREAD_COLUMNS) and set(skinbody.SPREAD_COLUMNS) <= set(b.skin_stats())
```

  in `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -619,8 +619,62 @@
     assert exc.value.code == 2 and "--surface-model elements" in capsys.readouterr().err
 
 
+@pytest.mark.parametrize("text,value", [("off", "off"), ("60", 60.0), ("72.5", 72.5), ("0.5", 0.5), ("179.5", 179.5)])
+def test_the_spread_limit_reads_off_or_a_fold_angle(text, value):
+    assert cli.parse_spread_limit(text) == value
+
+
+def test_the_spread_limit_marks_the_run_s_name_and_is_off_by_default():
+    """--spread-limit (sub-plan 22 Task 5c) is off by default, so that every existing run keeps its name; a limit ends a skin
+    run's settings in _spread-<deg>, after the spraying's and before the seed's."""
+    p = cli.build_parser()
+    base = "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-sesam_melt-girin"
+    name = lambda argv: cli.run_name_of(p.parse_args(argv), 0.1, 7500.0, 77500.133, "us76", cli.continuum_step(p.parse_args(argv)))
+    melt = THERMAL + ["--melt", "on"]
+    assert p.parse_args(melt).spread_limit == "off" and cli.spread_limit_of(p.parse_args(melt)) is None
+    assert name(melt) == name(melt + ["--spread-limit", "off"]) == base
+    assert name(melt + ["--spread-limit", "60"]) == base + "_spread-60"
+    assert name(melt + ["--spread-limit", "72.5"]) == base + "_spread-72.5"
+    assert (name(melt + ["--slurry-spray", "on", "--spread-limit", "60", "--seed", "7"])
+            == base + "_slurryspray-on_spread-60_seed-7")
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, spread_limit=45.0) == "_spread-45"
+
+
+@pytest.mark.parametrize("argv", [
+    THERMAL + ["--melt", "on", "--spread-limit", "0"],
+    THERMAL + ["--melt", "on", "--spread-limit", "180"],
+    THERMAL + ["--melt", "on", "--spread-limit", "-10"],
+    THERMAL + ["--melt", "on", "--spread-limit", "nan"],
+    THERMAL + ["--melt", "on", "--spread-limit", "inf"],
+    THERMAL + ["--melt", "on", "--spread-limit", "steep"],
+    THERMAL + ["--spread-limit", "60"],                                                # nothing melts, so nothing re-freezes
+    FEM + ["--melt", "on", "--spread-limit", "60"],                                    # the element model has no hand-down
+    BASE + ["--spread-limit", "60"],                                                   # no thermal model at all
+])
+def test_a_bad_spread_limit_or_one_without_a_melting_skin_run_exits_2(argv, tmp_path):
+    with pytest.raises(SystemExit) as exc:
+        cli.main(argv + ["--t-max", "1", "--outdir", str(tmp_path)])
+    assert exc.value.code == 2
+
+
 class _Built(Exception):
     pass
+
+
+def test_the_command_line_gives_the_skins_their_spread_limit(tmp_path, monkeypatch):
+    """--spread-limit reaches the skin body as a number [deg], None when off, with the flow off too: the heat-only skins hand
+    metal down as well. The body is not built: its arguments are."""
+    from reentry_model import skinbody
+
+    def built(*args, **kwargs):
+        raise _Built(kwargs)
+    monkeypatch.setattr(skinbody, "SkinBody", built)
+    melt = THERMAL + ["--melt", "on", "--fine-band", "0", "--t-max", "1", "--outdir", str(tmp_path)]
+    for argv, limit in ((melt, None), (melt + ["--spread-limit", "off"], None), (melt + ["--spread-limit", "60"], 60.0),
+                        (melt + ["--skin-flow", "off", "--spread-limit", "75"], 75.0)):
+        with pytest.raises(_Built) as got:
+            cli.main(argv)
+        assert got.value.args[0]["spread_limit"] == limit
 
 
 def test_the_command_line_builds_the_skins_with_the_shared_spray_model(tmp_path, monkeypatch):
```

  and in `tests/test_reentry_model_restart.py`:

```diff
--- a/tests/test_reentry_model_restart.py
+++ b/tests/test_reentry_model_restart.py
@@ -176,17 +176,114 @@
     assert same == flight[14:] and doc["settings"]["seed"] == cli.DEFAULT_SEED
 
 
-@pytest.mark.parametrize("case", ["missing", "format", "t_max", "dt_continuum", "t_max_nan", "dt_continuum_nan"])
+def test_resuming_with_a_limit_adds_the_spreading_columns_and_off_drops_them(flown, tmp_path):
+    """The rows of a run all take the first row's columns, the saved one's in a resume, so a limit given to a resume adds the
+    spreading columns to the saved row, or the later rows' would be lost, and `off` drops them from a state saved with a
+    limit, or the later rows would raise KeyError (sub-plan 22 Task 5c). The test flight's state saved without a limit, as the
+    two real states were, resumed with one: the columns are in every row, zero in the saved one, the JSON records the limit in
+    its settings, its results and what was changed, and every other column is the flight's, the limit binding nowhere. Then
+    the flight saved with a limit: the resume keeps it unless told, and `off` gives back the flight's column set."""
+    out, flight, records, state_path = flown
+    names = ["spread_mass_kg", "spread_unplaced_kg", "spread_edges_over", "spread_jump_max_K"]
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "on", "--spread-limit", "60"]) == 0
+    rows = rows_of(str(tmp_path / "on.csv"))
+    assert len(rows) == 3 and list(rows[0]) == list(flight[14]) + names and all(list(r) == list(rows[0]) for r in rows)
+    assert all(float(rows[0][k]) == 0.0 for k in names)
+    assert [{k: v for k, v in r.items() if k not in names} for r in rows] == flight[14:]
+    doc = json.load(open(str(tmp_path / "on.json")))
+    assert doc["resumed_from"]["changed"] == {"--spread-limit": 60.0} and doc["settings"]["spread_limit_deg"] == 60.0
+    assert doc["results"]["spread_limit_deg"] == 60.0 and doc["results"]["spread_mass_kg"] == 0.0
+
+    limited = str(tmp_path / "limited")
+    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", limited, "--name", "f", "--save-state-at", "0.1875",
+                              "--spread-limit", "60"]) == 0
+    saved = os.path.join(limited, "f.state.pkl")
+    doc = json.load(open(os.path.join(limited, "f.json")))
+    assert doc["settings"]["skin"]["spread_limit_deg"] == 60.0 and doc["results"]["spread_limit_deg"] == 60.0
+    assert list(rows_of(os.path.join(limited, "f.csv"))[0]) == list(flight[0]) + names
+
+    assert cli.main(["resume", saved, "--outdir", str(tmp_path), "--name", "kept"]) == 0
+    assert list(rows_of(str(tmp_path / "kept.csv"))[0]) == list(flight[0]) + names
+    assert json.load(open(str(tmp_path / "kept.json")))["results"]["spread_limit_deg"] == 60.0
+
+    assert cli.main(["resume", saved, "--outdir", str(tmp_path), "--name", "off", "--spread-limit", "off"]) == 0
+    assert rows_of(str(tmp_path / "off.csv")) == flight[14:]
+    doc = json.load(open(str(tmp_path / "off.json")))
+    assert doc["resumed_from"]["changed"] == {"--spread-limit": "off"}
+    assert "spread_limit_deg" not in doc["results"] and "spread_mass_kg" not in doc["results"]
+    assert "spread_limit_deg" not in doc["settings"]
+
+
+@pytest.mark.parametrize("case", ["missing", "format", "t_max", "dt_continuum", "t_max_nan", "dt_continuum_nan",
+                                  "spread_0", "spread_180", "spread_nan", "save_not_after", "save_at_end", "save_nan"])
 def test_the_resume_command_refuses_what_it_cannot_resume(case, flown, tmp_path):
     """A missing state file, a state of another format, an end of flight not after the saved time, a continuum step outside
-    (0, the saved flight's --dt of 0.025 s] and a value that is not a number are refused with exit 2 (bad input), and
-    nothing is written."""
+    (0, the saved flight's --dt of 0.025 s], a spread limit outside (0, 180) degrees and a value that is not a number are
+    refused with exit 2 (bad input), and nothing is written."""
     out, flight, records, state_path = flown
     bad_format = str(tmp_path / "other.pkl")
     restart.save(bad_format, {"format": restart.FORMAT + 1})
     state, options = {"missing": (str(tmp_path / "none.pkl"), []), "format": (bad_format, []),
                       "t_max": (state_path, ["--t-max", "0.1875"]), "dt_continuum": (state_path, ["--dt-continuum", "0.05"]),
-                      "t_max_nan": (state_path, ["--t-max", "nan"]), "dt_continuum_nan": (state_path, ["--dt-continuum", "nan"])}[case]
+                      "t_max_nan": (state_path, ["--t-max", "nan"]), "dt_continuum_nan": (state_path, ["--dt-continuum", "nan"]),
+                      "spread_0": (state_path, ["--spread-limit", "0"]), "spread_180": (state_path, ["--spread-limit", "180"]),
+                      "spread_nan": (state_path, ["--spread-limit", "nan"]),
+                      "save_not_after": (state_path, ["--save-state-at", "0.1875"]),
+                      "save_at_end": (state_path, ["--save-state-at", "0.2125"]),
+                      "save_nan": (state_path, ["--save-state-at", "nan"])}[case]
     with pytest.raises(SystemExit) as exc:
         cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", *options])
     assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
+
+
+def test_a_spread_limit_is_refused_on_a_flight_that_does_not_melt(tmp_path):
+    """A flight that does not melt has nothing re-frozen to spread, so a limit given to its resume exits 2, as `run` refuses
+    it; `off` asks nothing of the flight, so that goes through."""
+    pytest.importorskip("cantera")
+    k = FLIGHT.index("--melt")
+    argv = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--quiet"]              # without --melt on and --skin-flow-cycles 2
+    assert cli.main(argv + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "dry", "--save-state-at", "0.05"]) == 0
+    state = str(tmp_path / "dry.state.pkl")
+    with pytest.raises(SystemExit) as exc:
+        cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--spread-limit", "60"])
+    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
+    assert cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--spread-limit", "off"]) == 0
+
+
+def test_a_resumed_flight_saves_its_state_on_its_way_and_that_state_resumes(flown, tmp_path):
+    """A resumed flight can save its state as the first did (final review, 2026-11-07): `--save-state-at 0.2` on the resume
+    of the state saved at 0.1875 s writes <outdir>/<name>.state.pkl at the end of the 15th step, though the original
+    flight's own save had set the run's record of one; the JSON records the time and the file, the new state carries the
+    resumed run's name and arguments, and resumed in its turn it repeats the unbroken flight from its step to the last bit.
+    The same flight resumed without the option saves nothing."""
+    out, flight, records, state_path = flown
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r", "--save-state-at", "0.2"]) == 0
+    doc = json.load(open(str(tmp_path / "r.json")))
+    new = str(tmp_path / "r.state.pkl")
+    assert doc["settings"]["save_state_at_s"] == 0.2 and doc["results"]["state_saved_at_s"] == pytest.approx(0.2)
+    assert doc["results"]["state_save_error"] is None and doc["files"]["state"] == new
+    assert doc["resumed_from"]["changed"] == {"--save-state-at": 0.2}
+    state = restart.load(new)
+    assert state["step"] == 15 and state["meta"]["name"] == "r" and state["meta"]["args"]["save_state_at"] == 0.2
+    assert state["meta"]["args"]["outdir"] == str(tmp_path)
+    assert cli.main(["resume", new, "--outdir", str(tmp_path), "--name", "again"]) == 0
+    assert rows_of(str(tmp_path / "again.csv")) == flight[15:]
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "plain"]) == 0
+    assert not os.path.exists(str(tmp_path / "plain.state.pkl"))
+    assert "state" not in json.load(open(str(tmp_path / "plain.json")))["files"]
+
+
+def test_a_resume_says_what_it_is_doing_before_it_flies(flown, tmp_path, capsys):
+    """One line at the start of a resume, before the flight, names the state, its saved time, the end of flight and every
+    option in force, the spread limit or off included, so a flight's log says what it is without waiting for its JSON."""
+    out, flight, records, state_path = flown
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r", "--t-max", "0.225", "--seed", "7",
+                     "--spread-limit", "60", "--save-state-at", "0.2"]) == 0
+    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")]
+    assert len(line) == 1
+    for text in (os.path.abspath(state_path), "saved at 0.1875 s", "step 14", "to 0.225 s", "seed 7", "spread limit 60 deg",
+                 "save state at 0.2 s", "dt_continuum 0.0125", "slurry spray off", os.path.join(str(tmp_path), "r")):
+        assert text in line[0]
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "q"]) == 0
+    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
+    assert "saved random state" in line and "spread limit off" in line and "save state at no" in line
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spreading.py -q`. Expected: `ImportError: cannot
  import name 'spreading'`. The others, with the module in place: failures on `spread_limit`, `parse_spread_limit` and `--spread-limit`,
  except the refusals, which an unknown option also exits 2 on, and which mean something once the option exists.

- [ ] **Step 3: The pass**, `reentry_model/spreading.py`:

```python
"""Spreading the re-frozen deposit (sub-plan 22 Task 5c; Asha's idea and decisions of 2026-11-06,
docs/superpowers/specs/2026-11-06-deposit-spreading-design.md).

Both skin flights of 2026-11-05 stopped because a corner's move could no longer be absorbed by the 3D mesh, at the toe of a
ring of re-frozen metal travelling outward ahead of the flow: at the 100 mm sphere's worst edge the base steps about 1.4 mm
across a triangle about 1.5 mm wide. A 2 mm triangle cannot carry a ridge narrower than itself, so where the base as it
would stand after this part's move folds past a limit angle, part of the metal a corner is handing down goes to its lower
neighbouring corners instead: the angle-of-repose relaxation used for loose material on terrain, applied to the deposit
(`SkinBody._move_base` calls `spread_hand_down` after `skin.hand_down`). It moves solid metal sideways with no physical
mechanism, so the deposit's shape at the 2 mm scale is then set by the limit, and the sensitivity to it has to be measured.

Only handed metal moves: draws and absorption stay where they are, and the skins, which have already given the metal up,
do not change. Each pass finds the edges of the base as it would stand that fold past the limit. Across a ridge along an edge
(a, b), a and b give to the corners opposite it, c and e; across a valley, c and e give to a and b. Each of those four pairs
wants the share SPREAD_SHARE of the fold's excess over its angle, of the metal in the donor's dual cell above the receiver.
No donor falls below what balances its draw, so its move stays within the part's allowance, and no receiver rises above that
allowance: the part's move ratio holds, and a hand-down still never splits a part (spec 6.6, amended 2026-10-27). The metal
keeps the donor's mean enthalpy. SPREAD_PASSES passes a part at most, so metal travels at most a ring or two of 2 mm
neighbours in one part. What is left wanting to move after the last pass, and could leave its donor, is reported, as
`held_mass` reports what the tolerance holds. A limit that does not bind returns the very arrays it was given, so a run under
it is the run without it to the last bit. `SPENT` drops a true remainder of at most 1e-14 of one corner's hand-down per
part: a leak by design, far below the books' 1e-12, not round-off."""
from dataclasses import dataclass

import numpy as np

SPREAD_PASSES = 4       # passes a part at most: metal travels at most a ring or two of 2 mm neighbours in one part
SPREAD_SHARE = 0.25     # of the fold's excess over a donor's four pairs
SPENT = 1e-14           # of a donor's hand-down: round-off of a donor that gave all it had, kept as exactly nothing


@dataclass
class Spread:
    handed: np.ndarray    # kg per column, the hand-down after spreading; the argument itself if the limit does not bind
    enthalpy: np.ndarray  # J per column, likewise
    moved: float          # kg moved between columns
    unplaced: float       # kg still wanting to move after the last pass
    over: int             # edges over the limit on the first estimate
    landed: np.ndarray    # bool per column: metal arrived here in some pass, whatever it gave onward


def fold_angles(columns, points):
    """Per edge of the columns' triangulation with corners at `points` [m, (n, 3)]: the angle [deg] between its two
    triangles' normals, the measure of `movemesh.max_dihedral` and so of `surface_fold`, zero on an open boundary; and
    whether it is convex, a ridge along the edge, rather than a valley."""
    p = np.asarray(points, dtype=float)
    a, b, c = (p[columns.faces[:, v]] for v in range(3))
    unit = np.cross(b - a, c - a)
    unit /= np.linalg.norm(unit, axis=1)[:, None]
    angle, convex = np.zeros(len(columns.edges)), np.zeros(len(columns.edges), dtype=bool)
    two = np.flatnonzero(columns._tri[:, 1] >= 0)
    first, second = unit[columns._tri[two, 0]], unit[columns._tri[two, 1]]
    angle[two] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", first, second), -1.0, 1.0)))
    beyond = p[columns._opp[two, 1]] - p[columns.edges[two, 0]]         # the second triangle's far corner, under the
    convex[two] = np.einsum("ij,ij->i", beyond, first) < 0.0             # first's plane where the edge is a ridge
    return angle, convex


def _pairs(columns, q, normal, area, rho, angle, convex, over, limit_deg):
    """The donors, the receivers and the mass [kg] each pair wants to move, on the base at `q`, for the edges `over` the
    limit (an open boundary's edge has no fold, so no -1 corner comes in): the four pairs an edge makes, those whose donor
    stands above its receiver along the donor's normal."""
    ends, opposite = columns.edges[over].T, columns._opp[over].T
    high, low = np.where(convex[over], ends, opposite), np.where(convex[over], opposite, ends)
    donor = np.concatenate([high[0], high[0], high[1], high[1]])
    receiver = np.concatenate([low[0], low[1], low[0], low[1]])
    dz = np.einsum("ij,ij->i", q[donor] - q[receiver], normal[donor])
    above = dz > 0.0
    fold = np.tile(angle[over], 4)[above]
    want = SPREAD_SHARE * (fold - limit_deg) / fold * rho * area[donor[above]] * dz[above]
    return donor[above], receiver[above], want


def _within(want, who, room, n):
    """`want` [kg per pair] scaled so that what the pairs of each column `who` add up to stays within its `room` [kg]."""
    asked = np.bincount(who, want, n)
    return want * np.minimum(room / np.where(asked > 0.0, asked, 1.0), 1.0)[who]


def spread_hand_down(columns, points, normal, area, rho, drawn, handed, enthalpy, cap, limit_deg, passes=SPREAD_PASSES):
    """Spread the part's hand-down (`handed` [kg] and `enthalpy` [J] per column) over the corners at `points` [m], outward
    `normal`, dual `area` [m2], where the base as it would stand folds past `limit_deg` (module docstring): a Spread. `drawn`
    [kg] is the part's draws and `cap` [kg] its allowance per column (`SkinBody._hand_cap`). The base's estimate leaves out
    the part's absorption and smoothing. The sums of the hand-down and its enthalpy are kept to round-off. `landed` is read
    pass by pass, not from the totals: a corner that takes metal and gives more onward ends handing down less than it started
    with, and carries its neighbours' heat all the same."""
    n, first, moved = len(handed), 0, 0.0
    h, H = handed, enthalpy
    landed = np.zeros(n, dtype=bool)
    floor = np.maximum(drawn - cap, 0.0)                 # a donor keeps what balances its draw: the part's move ratio holds
    for k in range(passes + 1):
        q = points - normal * ((drawn - h) / (rho * area))[:, None]     # the base as it would stand
        angle, convex = fold_angles(columns, q)
        over = np.flatnonzero(angle > limit_deg)
        if k == 0:
            if not over.size:
                return Spread(handed, enthalpy, 0.0, 0.0, 0, landed)
            first = int(over.size)
        if not over.size:
            return Spread(h, H, moved, 0.0, first, landed)
        donor, receiver, want = _pairs(columns, q, normal, area, rho, angle, convex, over, limit_deg)
        want = _within(want, donor, np.maximum(h - floor, 0.0), n)
        if k == passes:
            return Spread(h, H, moved, float(want.sum()), first, landed)
        want = _within(want, receiver, np.maximum(cap - h, 0.0), n)
        going = want > 0.0
        donor, receiver, want = donor[going], receiver[going], want[going]
        unit = H[donor] / h[donor]                       # the metal keeps the donor's mean enthalpy
        gave = np.bincount(donor, want, n)
        rest = 1.0 - gave / np.where(h > 0.0, h, 1.0)
        rest = np.where(rest < SPENT, 0.0, rest)         # all it had: exactly nothing, or the mean enthalpy the moving
        got = np.bincount(receiver, want, n)             # mesh books its sweep at would be round-off over round-off
        h = h * rest + got
        H = H * rest + np.bincount(receiver, want * unit, n)
        landed = landed | (got > 0.0)
        moved += float(want.sum())
```

- [ ] **Step 4: The body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -99,7 +99,18 @@
 same particles_summary.csv -- the release angle's mean by mass, least and largest, b and the resolved sheared depth over
 Girin's delta_m on the wet windward columns (`sheared_depth`), and each rebuild's mass residual. Where `spool` names a
 directory the records go to disk as they come (`spray.SourceSpool`): a flight at the continuum step makes tens of
-millions of them."""
+millions of them.
+
+Spreading the deposit (sub-plan 22 Task 5c; `spread_limit`, None by default; Asha's idea and decisions of 2026-11-06, the
+design of that date; `spreading`): where the base as it would stand after a part's move folds past the limit, part of the
+metal a corner hands down goes to its lower neighbours' corners in step 4, after the hand-down and before the corners move.
+The skins have already given up what they handed down, so nothing in them changes: the corners move, and each column is
+reconciled, by the spread hand-down, whose sum is the hand-down's. The metal lands hot: it crosses a receiving base at the
+node's enthalpy and the rest of its heat goes to that node through `pending` (`_move_base`). `held_mass`, `handed_mass` and
+the fault count keep the hand-down as it left each skin. The history's spreading columns are `SPREAD_COLUMNS`, with a limit
+only (`spread_stats`), and a limit that does not bind changes nothing to the last bit. A body pickled before this has none of
+these attributes on the instance, so they are class defaults too: the saved flights of 2026-11-06 resume with or without a
+limit."""
 import copy
 import dataclasses
 import math
@@ -115,6 +126,7 @@
 from . import rheology
 from . import skin as skin_mod
 from . import spray as spray_mod
+from . import spreading
 from .body import ThermalBody
 from .thermal import SIGMA_SB
 
@@ -164,6 +176,11 @@
 
 SLURRY_SIGMA = 0.80               # N/m: the exposed slurry's surface tension in the cascade (the 2026-10-07 design 8.2, the
                                   # large-fragment design's labelled assumption for slurry)
+
+# the history's spreading columns (sub-plan 22 Task 5c), with a limit only (`SkinBody.spread_stats`): kg moved between columns
+# since time zero, kg still wanting to move and edges over the limit at the last move, and the largest landing [K] of metal over
+# its node's temperature since the last row
+SPREAD_COLUMNS = ("spread_mass_kg", "spread_unplaced_kg", "spread_edges_over", "spread_jump_max_K")
 
 
 @dataclasses.dataclass
@@ -530,12 +547,23 @@
     22 Task 1c), with `slurry_spray` the exposed slurry's own mode too (the cascade), and its record in `spray_plan`; its
     droplet records, one block of source-table rows per macro step, in `source_rows` where `particles`, and the size
     histograms in `hist_n` and `hist_m` (sub-plan 22 Task 2: the element model's interface to the coupled writer); with
-    `spool`, a directory, the rows go there as they come (`spray.SourceSpool`, Task 4) instead of staying in memory."""
+    `spool`, a directory, the rows go there as they come (`spray.SourceSpool`, Task 4) instead of staying in memory.
+    `spread_limit` [deg, within (0, 180)]: the fold past which a part's hand-down is spread to lower neighbouring corners
+    (`spreading`; sub-plan 22 Task 5c), None for no spreading."""
+    # class defaults as well as the instance's: the flights saved on 2026-11-06 were pickled before spreading existed
+    spread_limit = None                  # deg
+    spread_mass = spread_unplaced = 0.0  # kg: moved between columns since time zero; still wanting to move at the last move
+    spread_over = 0                      # edges over the limit on the last move's first estimate
+    spread_jump = 0.0                    # K: the largest landing of metal over its node's temperature this macro step
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
-                 surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True, spool=None):
+                 surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True, spool=None,
+                 spread_limit=None):
+        if spread_limit is not None and not 0.0 < spread_limit < 180.0:
+            raise ValueError("spread_limit must be None or an angle within (0, 180) degrees, got {!r}".format(spread_limit))
+        self.spread_limit = spread_limit
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
         self.flow_cycles = flow_cycles
@@ -644,6 +672,7 @@
         `state` (`flow_forcing`), or `forcing`, a FlowForcing given instead (a prescribed push, in the tests)."""
         n, self._ratio, self.step_parts = self._hint, 0.0, 0
         self.step_spray, self._released_step = _no_step_spray(), 0.0
+        self.spread_jump = 0.0
         if self.flow_on:
             end = forcing if forcing is not None else self.flow_forcing(state)
             last = self._forcing_prev
@@ -908,7 +937,9 @@
         metal within the fault tolerance of the node (spec 6.6, amended 2026-10-26); what the tolerance and the cap keep
         back, against a hand-down with neither, is recorded (`held_mass`). The bottom cell is then judged against the
         metal that crossed and, for a hand-down, the cell it uncovered -- the rest of a cell the cap took in part, at its
-        own temperature, where the cap ends inside one (`jump_faults`)."""
+        own temperature, where the cap ends inside one (`jump_faults`). With `spread_limit` the hand-down is first spread
+        over the corners (module docstring): the move and the reconciliation read the spread hand-down, `handed_s` with
+        its enthalpy `landed_H`, and the fault count and `handed_mass` what left each skin, `handed`."""
         mat, sk, mm, col = self.material, self.skins, self.moving, self.columns
         started = time.perf_counter()
         T_node = self.solver.temperature()[mm.corners]
@@ -920,19 +951,49 @@
         uncovered = np.where(handed > 0.0, sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)], np.nan)  # as it left it
         held = np.maximum(able - handed, 0.0)                 # what the tolerance and the cap keep in the skins
         self.held_mass, self.held_depth_max = float(held.sum()), float((held / (mat.rho * col.area)).max())
+        handed_s, landed_H, excess = handed, handed_H, None                       # what crosses each base; the same arrays
+        if self.spread_limit is not None:                                         # without a limit (sub-plan 22 Task 5c)
+            s = spreading.spread_hand_down(col, col.points, col.normal, col.area, mat.rho, book.drawn, handed, handed_H,
+                                           self._hand_cap(), self.spread_limit)
+            self.spread_mass += s.moved
+            self.spread_unplaced, self.spread_over = s.unplaced, s.over
+            handed_s, landed_H = s.handed, s.enthalpy
+            # Spread metal lands hot (Asha, 2026-11-06). The 5 K fault tolerance does not apply where it lands: the lower
+            # neighbours of a toe are the cold corners ahead of it, across folded edges the 3D temperature jumps by a median
+            # 62 K and 31 K on the two spheres and about 90 % of the edges exceed 5 K, so the rule would refuse most of the
+            # spreading it is for, and hot metal freezing on cold metal is what happens at a real toe. Metal landed here from a
+            # neighbour crosses this base at this node's enthalpy, and the rest of its heat goes to this node through
+            # `pending`, the deferred per-corner heat the next part's 3D solve takes as a load and the energy books count:
+            # booked through `h_outward` instead, `MovingMesh.move` would share the correction over all the moving nodes, the
+            # whole zone rather than this corner. Added below, after the reconciliation: `_advance_part` has by then set
+            # `pending` to this part's mismatch, and a part split or repeated has returned or restored before it, so the heat
+            # is applied once. T_node is read above, before the move, which may rebuild the 3D mesh. Where metal landed is
+            # the pass's `landed`, not handed_s > handed: a corner that took hot metal and gave more onward ends handing down
+            # less than it started with, and its neighbours' heat would cross through `h_outward` (Task 5c's review measured
+            # 1.9-7.0 % as much heat as went through `pending` going that way on a test plate). One that passed on all it
+            # took has no metal here and no enthalpy left to book: its excess is zero.
+            metal = s.landed & (handed_s > 0.0)
+            if metal.any():
+                jump = mat.temperature_from_enthalpy(landed_H[metal] / handed_s[metal]) - T_node[metal]
+                self.spread_jump = max(self.spread_jump, float(jump.max()))
+            excess = np.where(s.landed, landed_H - handed_s * mat.enthalpy(T_node), 0.0)
+            landed_H = landed_H - excess
         if mat.melts:
             self.momentum_lost += sk.momentum_out.sum(axis=0)                     # what the handed-down cells held
-        moved = (book.drawn - handed) / (mat.rho * col.area)                      # m inward by drawing and handing down
+        moved = (book.drawn - handed_s) / (mat.rho * col.area)                    # m inward by drawing and handing down
         depth = self._absorb_depth(np.maximum(MOVE_MAX * col.edge_length.mean() - moved, 0.0))
-        inward = book.drawn / mat.rho + depth * col.area - handed / mat.rho          # m3 per corner
+        inward = book.drawn / mat.rho + depth * col.area - handed_s / mat.rho        # m3 per corner
         u = -col.normal * (inward / col.area)[:, None] + mm.neutral_smoothing(col.normal, col.area)
         rebuilds = mm.rebuilds
         h_inward = lambda d: skin_mod.profile_enthalpy(mat, T_node, g0, 0.0, d)
-        h_outward = np.where(handed > 0.0, handed_H / np.where(handed > 0.0, handed, 1.0), np.nan)
+        # nan where nothing lands, as where a donor gave away all it handed down: any sweep there is at the corner's own enthalpy
+        h_outward = np.where(handed_s > 0.0, landed_H / np.where(handed_s > 0.0, handed_s, 1.0), np.nan)
         r = mm.move_with_policy(u, h_inward, h_outward)
         if mm.rebuilds > rebuilds:
             self.wall["rebuild"] += mm.last_rebuild.seconds
-        sk.exchange_bottom(mat.rho * r.swept_volume - (book.drawn - handed), r.swept_enthalpy - (book.drawn_H - handed_H))
+        sk.exchange_bottom(mat.rho * r.swept_volume - (book.drawn - handed_s), r.swept_enthalpy - (book.drawn_H - landed_H))
+        if excess is not None:
+            self.pending = self.pending + excess
         self.momentum_lost += sk.momentum_out.sum(axis=0)                         # what a reconciliation took away
         T_bottom = sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)]          # read after the reconciliation
         T_deep = T_node - g0 * (book.drawn / (mat.rho * col.area) + depth)   # the profile at the depth swept
@@ -1067,7 +1128,16 @@
                 # rebuild's mass residual, summed and the last
                 **self._depth_ratios(flowing),
                 "rebuild_mass_residual_kg": mm.rebuild_mass_residual_total,
-                "rebuild_mass_residual_last_kg": mm.last_rebuild.mass_residual if mm.last_rebuild is not None else 0.0}
+                "rebuild_mass_residual_last_kg": mm.last_rebuild.mass_residual if mm.last_rebuild is not None else 0.0,
+                **self.spread_stats()}
+
+    def spread_stats(self):
+        """The history's spreading columns (SPREAD_COLUMNS), last in `skin_stats`; none without a limit, so that every run
+        without one keeps its column set, and read here alone, as a resume adds them to its saved row or drops them from
+        it: a run's columns are its first row's (`coupled.CoupledRun.run`)."""
+        if self.spread_limit is None:
+            return {}
+        return dict(zip(SPREAD_COLUMNS, (self.spread_mass, self.spread_unplaced, self.spread_over, self.spread_jump)))
 
     def _depth_ratios(self, flowing):
         """Spec 10.2's b over delta_m and resolved sheared depth over delta_m (medians) on the wet windward columns: on
```

- [ ] **Step 5: The results and the column list's note** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -64,6 +64,8 @@
                     # each rebuild's mass residual beside its energy residual
                     "b_over_delta_m_median", "sheared_depth_over_delta_m_median",
                     "rebuild_mass_residual_kg", "rebuild_mass_residual_last_kg"]
+# with --spread-limit only, skinbody.SPREAD_COLUMNS follow these (sub-plan 22 Task 5c): a run's columns are its first row's, so
+# a resume that gives or drops a limit adds them to its saved row or takes them from it (cli.cmd_resume)
 # particles_summary.csv, one row a macro step: the element model's history columns
 SUMMARY_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um",
                    "theta_cr_deg", "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
@@ -274,6 +276,8 @@
                                                   "handed_mass_kg", "rebuild_residual_J", "wall_skins_s", "wall_3d_s",
                                                   "wall_move_s", "wall_rebuild_s", "wall_flow_s", "rebuild_mass_residual_kg")})
         out.update({k: int(last[k][-1]) for k in ("rebuilds", "jump_faults", "absorb_capped", "interface_repeats")})
+        if body.spread_limit is not None:                                       # the deposit's spreading (sub-plan 22 Task 5c)
+            out.update({"spread_limit_deg": body.spread_limit, "spread_mass_kg": body.spread_mass})
         return out
 
     def thermal_results(self, history):
```

- [ ] **Step 6: The command line** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -4,11 +4,12 @@
         [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
         [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
         [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
-        [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off ...]
+        [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off --spread-limit DEG|off ...]
                                      (the layered skins, the default with --thermal fem, their flow and spraying the
                                       default with --melt on; the element model with --surface-model elements)
     python -m reentry_model resume <state> --outdir ... --name ... [--t-max S --dt-continuum S --seed N --frames-every N
-        --slurry-spray on|off]       (a flight saved with --save-state-at, continued; sub-plan 22 Task 5b)
+        --slurry-spray on|off --spread-limit DEG|off]
+                                     (a flight saved with --save-state-at, continued; sub-plan 22 Task 5b)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
 Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
@@ -105,6 +106,24 @@
     return value
 
 
+def parse_spread_limit(text):
+    """--spread-limit: `off`, or a fold angle [degrees] within (0, 180)."""
+    if text == "off":
+        return text
+    try:
+        value = float(text)
+    except ValueError:
+        value = float("nan")
+    if not 0.0 < value < 180.0:
+        raise argparse.ArgumentTypeError("spread-limit must be off or an angle within (0, 180) degrees, got {!r}".format(text))
+    return value
+
+
+def spread_limit_of(args):
+    """The fold angle [degrees] past which the skins spread what they hand down, or None for --spread-limit off."""
+    return None if args.spread_limit == "off" else args.spread_limit
+
+
 def continuum_step(args):
     """The macro step after the continuum switch [s], or None for no switch: --dt-continuum if given (`off`: none), else
     DEFAULT_DT_CONTINUUM for the model proper -- melting with Girin removal -- and none for the instant-removal device
@@ -166,7 +185,7 @@
 
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
                      fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=True, flow_cycles=None,
-                     slurry_spray=False):
+                     slurry_spray=False, spread_limit=None):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
@@ -175,8 +194,9 @@
     `_skinflow-off` for the heat-only skins of a melting run (`flow` False: the flow, sub-plan 21 Task 6c, is on by default
     with melting since sub-plan 22 Task 3b; None for a run that does not melt, which the flow does not concern), and with
     the flow on `_flowcycles<n>` where the flow sub-steps in each heat sub-step are fixed at n rather than chosen each heat
-    sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3); empty
-    for the defaults."""
+    sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3);
+    `_spread-<deg>` where what the skins hand down is spread past a fold limit (sub-plan 22 Task 5c); empty for the
+    defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -202,6 +222,8 @@
         out += "_flowcycles{}".format(flow_cycles)
     if slurry_spray:
         out += "_slurryspray-on"
+    if spread_limit is not None:
+        out += "_spread-{:g}".format(spread_limit)
     return out
 
 
@@ -260,6 +282,12 @@
                     help="VTK frame every n macro steps, the series new (default: the saved flight's setting)")
     rs.add_argument("--slurry-spray", choices=("on", "off"), default=None,
                     help="the exposed slurry sprays too (default: the saved flight's)")
+    rs.add_argument("--spread-limit", type=parse_spread_limit, default=None,
+                    help="[deg] the fold past which the metal handed down spreads to lower corners, or off (default: the saved "
+                         "flight's, off for one saved before the option existed)")
+    rs.add_argument("--save-state-at", type=float, default=None,
+                    help="[s] of flight: save the resumed flight's whole state once, at the end of the first macro step ending at "
+                         "or after it, as <outdir>/<name>{}; after the saved time and before --t-max (final review, 2026-11-07)".format(STATE_SUFFIX))
 
 
 def build_parser():
@@ -386,6 +414,10 @@
                     help="with --skin-flow on: on, the exposed slurry sprays too, once no liquid covers it, by its own mode "
                          "(the cascade of the layered-skin design, section 8); it sprays only within about a kelvin of the "
                          "liquidus, on a viscosity law with no data there (default off)")
+    me.add_argument("--spread-limit", type=parse_spread_limit, default="off",
+                    help="with --surface-model skin and --melt on: where the base as it would stand folds past this angle "
+                         "[deg, within (0, 180)], metal being handed down goes to lower neighbouring corners instead, and "
+                         "lands hot (sub-plan 22 Task 5c; off until the flights have tested it; default off)")
     me.add_argument("--skin-flow-cycles", type=int, default=None,
                     help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
                          "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
@@ -475,7 +507,7 @@
                                      flow_cycles=args.skin_flow_cycles, surface_flow=flow_model,
                                      spray=spray_model_of(args, mat) if flow else None,
                                      slurry_spray=args.slurry_spray == "on", particles=args.particles,
-                                     spool=spool if flow else None)
+                                     spool=spool if flow else None, spread_limit=spread_limit_of(args))
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray_model_of(args, mat)
@@ -508,6 +540,8 @@
                         "swell": args.skin_swell == "on", "flow": args.skin_flow,
                         "flow_cycles": "adaptive" if the_body.flow_cycles is None else the_body.flow_cycles,
                         "slurry_spray": args.slurry_spray}
+        if spread_limit_of(args) is not None:
+            info["skin"]["spread_limit_deg"] = spread_limit_of(args)
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
@@ -557,7 +591,8 @@
                                                                     args.skin_swell == "on", args.fine_band, args.fine_size,
                                                                     args.skin_fine_zone,
                                                                     skin_flow_of(args) == "on" if args.melt == "on" else None,
-                                                                    args.skin_flow_cycles, args.slurry_spray == "on"))
+                                                                    args.skin_flow_cycles, args.slurry_spray == "on",
+                                                                    spread_limit_of(args)))
 
 
 def cmd_run(args, parser):
@@ -635,6 +670,8 @@
     if args.slurry_spray == "on" and args.skin_flow != "on":
         parser.error("--slurry-spray on needs --skin-flow on (with --surface-model skin and --melt on): the exposed slurry "
                      "sprays by the skins' spraying from their top cells")
+    if spread_limit_of(args) is not None and not (args.surface_model == "skin" and args.melt == "on"):
+        parser.error("--spread-limit needs --surface-model skin and --melt on: it spreads the metal the skins hand down")
     if args.surface_model == "skin":
         for flag, given in (("--runoff off", args.runoff == "off"), ("--deep-runoff off", args.deep_runoff == "off"),
                             ("--molten-cascade off", args.molten_cascade == "off"),
@@ -848,8 +885,19 @@
         if args.slurry_spray == "on" and getattr(run.body, "spray", None) is None:
             parser.error("--slurry-spray on needs a saved flight with the skins' spraying (--skin-flow on)")
         changed["--slurry-spray"] = args.slurry_spray
+    skins = isinstance(run.body, skinbody.SkinBody)
+    if args.spread_limit is not None:
+        if spread_limit_of(args) is not None and not (skins and run.body.material.melts):
+            parser.error("--spread-limit needs a saved flight of melting skins: it spreads the metal they hand down")
+        changed["--spread-limit"] = args.spread_limit
     if args.seed is not None:
         changed["--seed"] = args.seed
+    if args.save_state_at is not None:
+        t_end = args.t_max if args.t_max is not None else run.sim.settings.t_max
+        if not (math.isfinite(args.save_state_at) and state["time_s"] < args.save_state_at < t_end):
+            parser.error("--save-state-at must be a finite time after the saved state's {:.6g} s and before the end of flight "
+                         "{:.6g} s, got {!r}".format(state["time_s"], t_end, args.save_state_at))
+        changed["--save-state-at"] = args.save_state_at
     # -- the point between loading the state and flying on, where each option is applied; keep the options above and
     # their application below together, in the order of resume_arguments --
     if "--t-max" in changed:
@@ -859,6 +907,11 @@
         if run.switched_at is not None:            # latched: the step in use is the continuum step
             settings.dt = args.dt_continuum
     settings.save_at = settings.save_path = settings.save_meta = None
+    state_path = None
+    if "--save-state-at" in changed:               # the original flight's save left these set; run() saves only while they are None
+        state_path = os.path.join(args.outdir, args.name + STATE_SUFFIX)
+        settings.save_at, settings.save_path = args.save_state_at, state_path
+        run.state_saved_at = run.state_save_error = None
     settings.output_dir = os.path.join(run_dir, "vtk")
     settings.frames = []                           # the new series lists only its own frames
     if "--frames-every" in changed:
@@ -866,31 +919,52 @@
         run._frame_dt, run._next_frame = args.frames_every * saved.dt, run.sim.t
     if "--slurry-spray" in changed:
         run.body.slurry_spray = args.slurry_spray == "on"
+    if "--spread-limit" in changed and skins:      # every row takes the saved row's columns, which have to say so
+        run.body.spread_limit = spread_limit_of(args)
+        for key in skinbody.SPREAD_COLUMNS:
+            state["row"].pop(key, None)
+        state["row"].update(run.body.spread_stats())
     if isinstance(run.body.source_rows, spray.SourceSpool):      # only the records after the resume are written
         run.body.source_rows = spray.SourceSpool(os.path.join(run_dir, SOURCE_SPOOL))
     if args.seed is not None:
         np.random.seed(args.seed)
     else:
         np.random.set_state(state["rng"])
-    history = run.run(resume=state)
-    # the arguments of the writer: the saved flight's, with what is changed here
-    for key, value in (("outdir", args.outdir), ("name", args.name), ("save_state_at", None)):
+    # the arguments of the writer, and of a state this run saves: the saved flight's, with what is changed here
+    for key, value in (("outdir", args.outdir), ("name", args.name), ("save_state_at", changed.get("--save-state-at"))):
         setattr(saved, key, value)
     for key, flag in (("t_max", "--t-max"), ("dt_continuum", "--dt-continuum"), ("seed", "--seed"),
-                      ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray")):
+                      ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray"), ("spread_limit", "--spread-limit")):
         if flag in changed:
             setattr(saved, key, changed[flag])
+    if state_path is not None:
+        settings.save_meta = {"name": args.name, "args": dict(vars(saved))}
+    limit = run.body.spread_limit if skins else None
+    print("resuming {} (saved at {:.6g} s, step {}) to {:.6g} s: dt_continuum {}, seed {}, frames every {}, slurry spray {}, "
+          "spread limit {}, save state at {}; files {}".format(
+              os.path.abspath(args.state), state["time_s"], state["step"], run.sim.settings.t_max,
+              settings.dt_continuum if settings.dt_continuum is not None else "none",
+              args.seed if args.seed is not None else "saved random state", settings.frames_every,
+              "on" if getattr(run.body, "slurry_spray", False) else "off",
+              "off" if limit is None else "{:g} deg".format(limit),
+              "{:g} s".format(args.save_state_at) if args.save_state_at is not None else "no",
+              os.path.join(args.outdir, args.name)), flush=True)
+    history = run.run(resume=state)
     atm, atm_name, atm_info = make_atmosphere(saved.atmosphere, saved.epoch, saved.wind)
     reference = sesam_io.load_reference(saved.reference) if saved.reference else None
     mass = body.sphere_mass(run.sim.settings.diameter, saved.material_density)
     thermal_info = {key: getattr(saved, key) for key in ("heating", "melt", "surface_model", "thermal_solver", "linear_solver",
                                                          "skin_flow", "slurry_spray")}
     thermal_info.update({"macro_step_s": saved.dt, "frames_every": settings.frames_every, "dt_continuum_s": settings.dt_continuum})
+    if state_path is not None:
+        thermal_info["save_state_at_s"] = args.save_state_at
+    if skins and run.body.spread_limit is not None:
+        thermal_info["spread_limit_deg"] = run.body.spread_limit
     resumed = {"state": os.path.abspath(args.state), "step": state["step"], "time_s": state["time_s"],
                "run_name": state["meta"]["name"], "changed": changed}
     os.makedirs(args.outdir, exist_ok=True)
     return finish_run(saved, args.name, run.sim.settings, mass, atm_name, atm_info, reference, thermal_info, run.body, history,
-                      None, resumed=resumed)
+                      state_path, resumed=resumed)
 
 
 def _fmt(x):
```

- [ ] **Step 7: Run the tests to see them pass.** Expected: the four files' new tests pass in about 2 minutes across three processes; the unit
  tier 587 passed, 1 skipped (Task 5b's 549 and these 38; 15.1 minutes across four processes).

- [ ] **Step 8: The production check** (in `measure/task5c/`, on mains power, through `flights/launch.sh`, then `compare.py` and
  `fold_place.py`): resume `d050__physics__skin.state.pkl` to 172.0 s with `--spread-limit 60 --frames-every 0` and compare with the snapshot
  flight's `rows.csv`: every shared column but the `wall_` timers the same to the last bit. Then to 177.0 s with `--spread-limit 60
  --frames-every 1`, recording the fold, `spread_mass_kg`, `spread_unplaced_kg`, the books and the wall time in `measure/task5c/check.md`.
  A fold above the limit there is a finding to report, not a reason to tune the pass.

**The final review's amendments (2026-11-07).** The whole-branch review found the code ready for flights and asked for
three things that stop a multi-hour flight being wasted; they fall into this task's diffs of `cli.py`, `spreading.py` and
`tests/test_reentry_model_restart.py`. (1) `resume` takes `--save-state-at S` (finite, after the saved time and before the
run's end of flight, else exit 2): the resumed flight saves its whole state once at the end of the first macro step ending
at or after S, to `<outdir>/<name>.state.pkl`, with the resumed run's name and arguments in the state, and the record of
the original flight's own save is reset so that the save happens; the JSON records `save_state_at_s` and the `state` file as
`run` does. A flight that gets past the old stop can then be looked at without re-flying (`test_a_resumed_flight_saves_its_state_on_its_way_and_that_state_resumes`
resumes from the new state to the unbroken flight's rows to the last bit). (2) `cmd_resume` prints one line to stdout after
validation and before the flight naming the state, its saved time and step, the end of flight, the continuum step, seed,
frames, slurry spray, spread limit (or off), the save time and the files (`test_a_resume_says_what_it_is_doing_before_it_flies`).
(3) `spreading.py`'s docstring says that `SPENT` drops a true remainder of at most 1e-14 of one corner's hand-down per part,
a leak by design far below the books' 1e-12, and fact 192 lists `--spread-limit` and `--save-state-at`. The resume
command's options are now `[--t-max --dt-continuum --seed --frames-every --slurry-spray --spread-limit --save-state-at]`.

- [ ] **Step 9: Commit.**

```bash
git add reentry_model/spreading.py reentry_model/skinbody.py reentry_model/coupled.py reentry_model/cli.py tests/test_reentry_model_spreading.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "Skins: what the skins hand down spreads to lower neighbouring corners where the base as it would stand folds past a limit (--spread-limit DEG|off on run and resume, off by default, reentry_model/spreading.py; Asha's idea and decisions of 2026-11-06), the metal landing hot -- across a ridge its corners give to the corners opposite, across a valley the reverse, no donor below what balances its draw and no receiver above the part's allowance, the metal crossing at the node's enthalpy any base where metal arrived in any pass (Spread.landed) and the rest of its heat going to that node through pending; the history's columns and the results with a limit only, a limit that does not bind changing nothing to the last bit, a resume adding the columns to its saved row or dropping them from it; on the 50 mm state to 177.0 s with a limit of 60 degrees the fold peaks at 64.86 against 116.75, 1.70 g spread against 13.63 g handed down, the books closed and no fault; to 176.0 s, the two flights' last shared macro step, the lowest 3D radius ratio is 0.0756 against 0.0678 and the deepest skin 2.22 mm against 3.04 mm; unit tier 587 passed, 1 skipped (sub-plan 22 Task 5c)"
```

---

### Task 5d: The part count set by melting, and the stall guard

**Files:**
- Modify: `reentry_model/skinbody.py` (`PART_COUNT_RULES`, `PARTS_MAX`, the class default and `SkinBody(..., part_count_from="moves")`,
  `_draw_ratio`, the hint line in `_advance_part`, the guard in `advance`)
- Modify: `reentry_model/cli.py` (`--part-count-from moves|melting` on `run` and `resume`, `part_count_rule_of`, the name suffix, the
  records, the refusals, the resume's one line)
- Test: `tests/test_reentry_model_skinbody.py` (7 tests), `tests/test_reentry_model_cli.py` (3 tests, 7 cases),
  `tests/test_reentry_model_restart.py` (2 tests)

**Interfaces:**
- Consumes: `SkinBody._move_ratio(book, m_min)`, `MOVE_MAX`, `MOVE_HINT_MARGIN`, `HAND_DOWN_SHARE`, `_hand_cap()`; Task 5c's
  `--spread-limit` wiring as the pattern for the option on `run` and `resume`.
- Produces: `PART_COUNT_RULES = ("moves", "melting")`; `PARTS_MAX = 1000`; `SkinBody(..., part_count_from="moves")`, also a class
  default for states pickled before it; `--part-count-from moves|melting` on `run` (default `moves`) and `resume` (default: the
  state's), exit 2 with `melting` unless `--surface-model skin --melt on` (run) or a saved flight of melting skins (resume); the run's
  name suffix `_parts-melting`; the run JSON's `settings.skin.part_count_from = "melting"` only when given, the resume JSON's
  `settings.part_count_from` only when the body's rule is `melting` and `changed["--part-count-from"]` when given.

**What was decided** (design `docs/superpowers/specs/2026-11-08-hand-down-and-part-count-design.md`; Asha, 2026-11-08, to fly both
and compare them). The flights with spreading (`--spread-limit 60`) got past the old stops and then stalled as the part count grew
by a quarter a macro step: on the 50 mm sphere 837, 1 047, ..., 9 765 parts, each the previous times 1.25 rounded up, on the 100
mm sphere from 1 to 1 537. The cause: since 2026-10-27 the hand-down is capped just under the move allowance, so it cannot split
a part, but `_move_ratio` still counts it, and where a capped hand-down fires in any part the ratio is just under 1 a part and
the next macro step's hint, `MOVE_HINT_MARGIN` times it, is a quarter more parts, whose hand-downs are capped again. With
`melting` the hint reads, for each kept part, the draws alone, `max(drawn / (rho area)) / (MOVE_MAX mean edge)` times
`dt_macro / dt`, kept parts only (a split part returns before it is maxed, as before); with `moves` the line is today's. **The split
test is the same under both rules.** The guard, in `advance` once `n` is final and before the first part, raises
`RuntimeError` if `n` is over `PARTS_MAX = 1000` (the controller, 2026-11-08), under both rules: it can only stop a run that would have
run away.

**Rulings** (the controller's, 2026-11-08): the guard's test lowers `PARTS_MAX` to 3 instead of running a thousand parts; the split
test's two rules are compared with no hand-down in the compared steps (the cap given to `skin.hand_down` is zero), where they must
give the same parts and rows; the capped-hand-down tests use a deposit nine move allowances deep and 0.05 s macro steps, under
the 0.0625 s part length, so that the length does not set the count: under `moves` 1, 2 and 3 parts, under `melting` 1, 1 and 1.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1593,3 +1593,153 @@
     b.spread_limit = 60.0
     b.advance(1.0, 0.5, zero)
     assert set(b.spread_stats()) == set(skinbody.SPREAD_COLUMNS) and set(skinbody.SPREAD_COLUMNS) <= set(b.skin_stats())
+
+
+# -- what sets the next macro step's part count (sub-plan 22 Task 5d; Asha 2026-11-08) ---------------------------------------
+def held_back_body(mesh, depth_allowances=9.0, **kw):
+    """The 4 mm sphere at 700 K with thirty front columns holding `depth_allowances` move allowances of metal at the node's
+    temperature at their bottoms and no heat: every part's hand-down is capped at just under one allowance (a deposit that
+    deep is not down in the macro steps flown), and nothing is drawn. Returns the body and the columns' deposit [kg]."""
+    b = make(mesh, T0=700.0, **kw)
+    b.emissivity = 0.0
+    mat, sk, col = b.material, b.skins, b.columns
+    cols = np.argsort(b.theta)[:30]
+    extra = np.zeros(sk.n)
+    extra[cols] = depth_allowances * skinbody.MOVE_MAX * col.edge_length.mean() * mat.rho * col.area[cols]
+    sk.add_bottom(extra, float(mat.enthalpy(700.0)))
+    sk.regrid()
+    return b, extra.sum()
+
+
+def timed_out(stats):
+    """A history row without its wall-time columns, which time the run and not the flight."""
+    return {k: v for k, v in stats.items() if not k.startswith("wall_")}
+
+
+def same_rows(a, b):
+    """Two lists of (parts, row) the same to the last bit, a column of not-a-number (nothing sprayed yet) the same as another."""
+    return len(a) == len(b) and all(
+        pa == pb and ra.keys() == rb.keys() and all(ra[k] == rb[k] or (ra[k] != ra[k] and rb[k] != rb[k]) for k in ra)
+        for (pa, ra), (pb, rb) in zip(a, b))
+
+
+def fly_three_macro_steps(b, dt=0.05):
+    """Three macro steps of `dt` s (not longer than the melting part length) with no heat; the parts each ran as."""
+    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    parts = []
+    for k in range(3):
+        b.advance(dt * (k + 1), dt, zero)
+        parts.append(b.step_parts)
+    return parts
+
+
+def test_under_the_moves_rule_a_capped_hand_down_grows_the_next_steps_parts_by_a_quarter(coarse_sphere_mesh):
+    """The defect of the spreading flights, pinned (sub-plan 22 Task 5d): a hand-down capped just under the move allowance
+    gives a move ratio just under 1 in every part, and the hint is 1.25 times the parts of this step, so the next step
+    runs a quarter more parts whose hand-downs are capped again -- 837, 1 047, ..., 9 765 parts on the 50 mm spreading flight,
+    each the previous times 1.25 rounded up (measured 2026-11-08). Here a deposit nine allowances deep, three 0.05 s
+    macro steps (under the 0.0625 s part length, so that the length does not set the count): 1, 2 and 3 parts, each the
+    previous times 1.25 rounded up, the allowance of HAND_DOWN_SHARE's 1e-6 included."""
+    assert skinbody.PART_MAX >= 0.05
+    b, _ = held_back_body(coarse_sphere_mesh)
+    assert b.part_count_from == "moves"
+    parts = fly_three_macro_steps(b)
+    assert parts[0] == 1
+    for before, after in zip(parts, parts[1:]):
+        assert after == int(np.ceil(1.25 * before * skinbody.HAND_DOWN_SHARE)) and after > before
+    assert b.held_mass > 0.0                                        # the deposit still waits: every hand-down was capped
+
+
+def test_under_the_melting_rule_a_capped_hand_down_leaves_the_part_count_alone(coarse_sphere_mesh):
+    """The same flight under `melting`: the hint reads the draws alone, and nothing is drawn, so every macro step runs as the
+    one part the draws and the part length call for, while the same metal comes down (a capped hand-down an allowance a
+    part, the rest held) and mass and energy close (sub-plan 22 Task 5d)."""
+    b, extra = held_back_body(coarse_sphere_mesh, part_count_from="melting")
+    M0, E0 = b.mass(), b.energy()
+    assert fly_three_macro_steps(b) == [1, 1, 1] and b.part_count_from == "melting"
+    assert b.handed_mass > 0.0 and b.held_mass > 0.0
+    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
+    assert abs(b.energy_balance_residual()) < 1e-8 and b.jump_faults == 0
+
+
+def test_a_draw_over_the_allowance_splits_a_part_under_either_rule(coarse_sphere_mesh, monkeypatch):
+    """The split test is the same under both rules (sub-plan 22 Task 5d): 6 MW/m2 from 880 K for three 0.5 s macro steps on
+    4 mm elements, with no melting cap on the parts, draws more than the allowance in a part, which is run again as shorter
+    parts. The two rules differ only through the hand-down, so it is taken out of the comparison (the cap given to
+    `skin.hand_down` is zero: the skins only grow, the first step of such a flight is the one the unbroken run splits
+    into 8 parts), and then the two give the same parts and the same history rows in every step. numpy's generator is
+    seeded as a run seeds it, since the heat solver's multigrid hierarchy draws from it."""
+    hand_down = skin.hand_down
+    monkeypatch.setattr(skin, "hand_down", lambda *a, cap=None, **k: hand_down(*a, cap=np.zeros_like(cap), **k))
+    runs = {}
+    for rule in skinbody.PART_COUNT_RULES:
+        np.random.seed(12345)
+        b = make(coarse_sphere_mesh, part_count_from=rule, part_max=None)
+        rows = []
+        for k in range(3):
+            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 6e6))
+            rows.append((b.step_parts, timed_out(b.skin_stats())))
+        runs[rule] = (b, rows)
+    moves, melting = runs["moves"], runs["melting"]
+    assert moves[0].handed_mass == 0.0 and moves[0].drawn_mass > 0.0
+    assert moves[1][0][0] > 1                                          # the draws split the first step's part
+    assert same_rows(moves[1], melting[1])
+    assert moves[0].mass() == melting[0].mass() and moves[0].energy() == melting[0].energy()
+
+
+def test_the_moves_rule_by_default_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
+    """`part_count_from="moves"` and the argument left out are one flight, every column of every row and the final mass and
+    energy equal to the last bit, the wall-time columns that time the run left out (sub-plan 22 Task 5d), on a flight that melts, hands down and splits parts (numpy's generator
+    seeded as a run seeds it before each)."""
+    runs = []
+    for kw in ({}, {"part_count_from": "moves"}):
+        np.random.seed(12345)
+        b = make(coarse_sphere_mesh, **kw)
+        rows = []
+        for k in range(3):
+            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 3e6))
+            rows.append((b.step_parts, timed_out(b.skin_stats())))
+        runs.append((b, rows))
+    assert same_rows(runs[0][1], runs[1][1]) and runs[0][1][1][0] > 1
+    assert runs[0][0].mass() == runs[1][0].mass() and runs[0][0].energy() == runs[1][0].energy()
+
+
+def test_a_macro_step_that_would_need_more_than_PARTS_MAX_parts_stops_the_run(coarse_sphere_mesh, monkeypatch):
+    """The stall guard (sub-plan 22 Task 5d; the controller, 2026-11-08): a macro step that would start with more than
+    PARTS_MAX parts raises RuntimeError naming the step, the count and the limit, before any part runs, under either rule;
+    PARTS_MAX parts themselves run. The limit is lowered to 3 so that the test does not run a thousand parts."""
+    monkeypatch.setattr(skinbody, "PARTS_MAX", 3)
+    for rule in skinbody.PART_COUNT_RULES:
+        b = make(coarse_sphere_mesh, T0=600.0, part_count_from=rule)
+        b.emissivity = 0.0
+        zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+        M0, E0 = b.mass(), b.energy()
+        b._hint = 4
+        with pytest.raises(RuntimeError, match=r"from 0\.5 s to 1 s would run as 4 parts, more than PARTS_MAX = 3"):
+            b.advance(1.0, 0.5, zero)
+        assert b.step_parts == 0 and b.mass() == M0 and b.energy() == E0       # nothing ran
+        b._hint = 3
+        b.advance(1.0, 0.5, zero)
+        assert b.step_parts == 3
+    assert skinbody.PARTS_MAX == 3 and skinbody.PART_COUNT_RULES == ("moves", "melting")
+
+
+def test_a_body_saved_without_the_part_count_rule_flies_on(coarse_sphere_mesh):
+    """The states saved before the rule existed have no `part_count_from` on the instance: the class's default gives them
+    `moves`, they fly on as before, and a resume can give them `melting` (sub-plan 22 Task 5d)."""
+    b, _ = held_back_body(coarse_sphere_mesh)
+    del vars(b)["part_count_from"]
+    assert "part_count_from" not in vars(b) and b.part_count_from == "moves"
+    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    b.advance(0.05, 0.05, zero)
+    b.advance(0.1, 0.05, zero)
+    assert b.step_parts == 2                                         # the moves rule: the hand-down fed the hint
+    b.part_count_from = "melting"
+    b.advance(0.15, 0.05, zero)
+    b.advance(0.2, 0.05, zero)
+    assert b.step_parts == 1
+
+
+def test_an_unknown_part_count_rule_is_refused(coarse_sphere_mesh):
+    with pytest.raises(ValueError, match="part_count_from"):
+        make(coarse_sphere_mesh, part_count_from="draws")
```

  in `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -726,3 +726,49 @@
     # the history prints nine significant digits
     assert r["dt_switch_kn"] == pytest.approx(float(rows[1]["knudsen"]), rel=1e-8) and r["dt_switch_kn"] > 0.01
     assert r["dt_switch_altitude_km"] == pytest.approx(float(rows[1]["altitude_km"]), rel=1e-8)
+
+
+def test_the_part_count_rule_reads_moves_or_melting_and_marks_the_run_s_name_only_for_melting():
+    """--part-count-from (sub-plan 22 Task 5d) is `moves` by default, so that every existing run keeps its name and its
+    record; `melting` ends a skin run's settings in _parts-melting, after the spread limit's and before the seed's."""
+    p = cli.build_parser()
+    base = "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-sesam_melt-girin"
+    name = lambda argv: cli.run_name_of(p.parse_args(argv), 0.1, 7500.0, 77500.133, "us76", cli.continuum_step(p.parse_args(argv)))
+    melt = THERMAL + ["--melt", "on"]
+    assert p.parse_args(melt).part_count_from == "moves" and cli.part_count_rule_of(p.parse_args(melt)) == "moves"
+    assert name(melt) == name(melt + ["--part-count-from", "moves"]) == base
+    assert name(melt + ["--part-count-from", "melting"]) == base + "_parts-melting"
+    assert (name(melt + ["--spread-limit", "60", "--part-count-from", "melting", "--seed", "7"])
+            == base + "_spread-60_parts-melting_seed-7")
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, part_count_from="melting") == "_parts-melting"
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True) == ""
+
+
+@pytest.mark.parametrize("argv", [
+    THERMAL + ["--melt", "on", "--part-count-from", "steep"],                          # no such rule
+    THERMAL + ["--melt", "on", "--part-count-from", "off"],
+    THERMAL + ["--part-count-from", "melting"],                                        # nothing melts, so nothing is drawn
+    FEM + ["--melt", "on", "--part-count-from", "melting"],                            # the element model has no parts
+    BASE + ["--part-count-from", "melting"],                                           # no thermal model at all
+])
+def test_a_bad_part_count_rule_or_the_melting_rule_without_a_melting_skin_run_exits_2(argv, tmp_path):
+    with pytest.raises(SystemExit) as exc:
+        cli.main(argv + ["--t-max", "1", "--outdir", str(tmp_path)])
+    assert exc.value.code == 2
+
+
+def test_the_command_line_gives_the_skins_their_part_count_rule_and_records_it_only_for_melting(tmp_path, monkeypatch):
+    """--part-count-from reaches the skin body as `moves` or `melting` (sub-plan 22 Task 5d). The body is not built: its
+    arguments are."""
+    from reentry_model import skinbody
+
+    def built(*args, **kwargs):
+        raise _Built(kwargs)
+    monkeypatch.setattr(skinbody, "SkinBody", built)
+    melt = THERMAL + ["--melt", "on", "--fine-band", "0", "--t-max", "1", "--outdir", str(tmp_path)]
+    for argv, rule in ((melt, "moves"), (melt + ["--part-count-from", "moves"], "moves"),
+                       (melt + ["--part-count-from", "melting"], "melting"),
+                       (melt + ["--skin-flow", "off", "--part-count-from", "melting"], "melting")):
+        with pytest.raises(_Built) as got:
+            cli.main(argv)
+        assert got.value.args[0]["part_count_from"] == rule
```

  and in `tests/test_reentry_model_restart.py`:

```diff
--- a/tests/test_reentry_model_restart.py
+++ b/tests/test_reentry_model_restart.py
@@ -287,3 +287,44 @@
     assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "q"]) == 0
     line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
     assert "saved random state" in line and "spread limit off" in line and "save state at no" in line
+
+
+def test_resuming_with_the_melting_rule_sets_it_on_the_body_and_records_it(flown, tmp_path, capsys):
+    """`--part-count-from melting` on a resume (sub-plan 22 Task 5d) sets the rule on the saved flight's skin body, which the
+    one line before the flight names, records it in the JSON's settings and in what was changed, and carries it into the
+    state the resumed flight saves, so a resume of that state keeps it unless told; `moves` asked of it gives back the
+    saved flight's record (no key), and a flight resumed without the option records none."""
+    out, flight, records, state_path = flown
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "m", "--part-count-from", "melting",
+                     "--save-state-at", "0.2"]) == 0
+    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
+    assert "part count from melting" in line
+    doc = json.load(open(str(tmp_path / "m.json")))
+    assert doc["resumed_from"]["changed"] == {"--part-count-from": "melting", "--save-state-at": 0.2}
+    assert doc["settings"]["part_count_from"] == "melting"
+    assert restart.load(str(tmp_path / "m.state.pkl"))["meta"]["args"]["part_count_from"] == "melting"
+    assert cli.main(["resume", str(tmp_path / "m.state.pkl"), "--outdir", str(tmp_path), "--name", "kept"]) == 0
+    assert "part count from melting" in capsys.readouterr().out
+    assert json.load(open(str(tmp_path / "kept.json")))["settings"]["part_count_from"] == "melting"
+    assert cli.main(["resume", str(tmp_path / "m.state.pkl"), "--outdir", str(tmp_path), "--name", "back",
+                     "--part-count-from", "moves"]) == 0
+    doc = json.load(open(str(tmp_path / "back.json")))
+    assert "part_count_from" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--part-count-from": "moves"}
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "plain"]) == 0
+    assert "part_count_from" not in json.load(open(str(tmp_path / "plain.json")))["settings"]
+    assert cli.main(FLIGHT + ["--t-max", "0.05", "--outdir", str(tmp_path), "--name", "run", "--part-count-from", "melting"]) == 0
+    assert json.load(open(str(tmp_path / "run.json")))["settings"]["skin"]["part_count_from"] == "melting"   # run's record
+
+
+def test_the_melting_rule_is_refused_on_a_resume_of_a_flight_that_does_not_melt(tmp_path):
+    """A flight that does not melt has no draws to set a part count by, so the melting rule asked of its resume exits 2,
+    as `run` refuses it; `moves` asks nothing of the flight, so that goes through (sub-plan 22 Task 5d)."""
+    pytest.importorskip("cantera")
+    k = FLIGHT.index("--melt")
+    argv = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--quiet"]              # without --melt on and --skin-flow-cycles 2
+    assert cli.main(argv + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "dry", "--save-state-at", "0.05"]) == 0
+    state = str(tmp_path / "dry.state.pkl")
+    with pytest.raises(SystemExit) as exc:
+        cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--part-count-from", "melting"])
+    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
+    assert cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--part-count-from", "moves"]) == 0
```

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError` and `TypeError` on `part_count_from` and
  `PART_COUNT_RULES`; the refusals of an unknown rule, which exit 2 without the option too, pass.

- [ ] **Step 3: The body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -163,6 +163,15 @@
                                   # capped mass by the area again, never makes a capped hand-down split a part
 SPLIT_DEPTH_MAX = 6               # a part is split at most this many times over
 MOVE_HINT_MARGIN = 1.25           # the next macro step's parts: this times the parts this one's largest move called for
+# What sets the next macro step's part count (sub-plan 22 Task 5d; Asha, 2026-11-08). `moves`, the rule so far, feeds the
+# hint with the move the draws and the capped hand-down together call for (`_move_ratio`); a capped hand-down is
+# just under the whole allowance, so wherever one fires the ratio is just under 1 per part and the next step runs 1.25
+# times as many parts, whose hand-downs are capped again: on the 50 mm spreading flight 837, 1 047, ..., 9 765 parts,
+# each the previous times 1.25 rounded up (measured 2026-11-08), the metal still held back and the flight stalled.
+# `melting` feeds the hint with the draws alone, so the count follows the melting; the split test is the same under
+# both. Asha decided on 2026-11-08 to fly both and compare them.
+PART_COUNT_RULES = ("moves", "melting")
+PARTS_MAX = 1000                  # parts: a macro step that would start with more stops the run (the controller, 2026-11-08)
 # The part length while melting (2026-10-12). On the 100 mm sphere's first 1.5 s of melting, 0.5 s parts cost 1 454 s
 # (53 % on parts rejected for their moves, 4.4 passes a kept part) against 341 s for 0.0625 s macro steps (1.8 passes);
 # before melting short parts cost about 18 % more, so the cap holds only while some skin's top is non-rigid. (Those
@@ -549,9 +558,12 @@
     histograms in `hist_n` and `hist_m` (sub-plan 22 Task 2: the element model's interface to the coupled writer); with
     `spool`, a directory, the rows go there as they come (`spray.SourceSpool`, Task 4) instead of staying in memory.
     `spread_limit` [deg, within (0, 180)]: the fold past which a part's hand-down is spread to lower neighbouring corners
-    (`spreading`; sub-plan 22 Task 5c), None for no spreading."""
+    (`spreading`; sub-plan 22 Task 5c), None for no spreading.
+    `part_count_from`: `moves` (the default) or `melting`, what the next macro step's part count follows
+    (`PART_COUNT_RULES`; sub-plan 22 Task 5d)."""
     # class defaults as well as the instance's: the flights saved on 2026-11-06 were pickled before spreading existed
     spread_limit = None                  # deg
+    part_count_from = "moves"            # the part count's rule (sub-plan 22 Task 5d): the states saved before it have none
     spread_mass = spread_unplaced = 0.0  # kg: moved between columns since time zero; still wanting to move at the last move
     spread_over = 0                      # edges over the limit on the last move's first estimate
     spread_jump = 0.0                    # K: the largest landing of metal over its node's temperature this macro step
@@ -560,7 +572,10 @@
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
                  surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True, spool=None,
-                 spread_limit=None):
+                 spread_limit=None, part_count_from="moves"):
+        if part_count_from not in PART_COUNT_RULES:
+            raise ValueError("part_count_from must be one of {}, got {!r}".format(PART_COUNT_RULES, part_count_from))
+        self.part_count_from = part_count_from
         if spread_limit is not None and not 0.0 < spread_limit < 180.0:
             raise ValueError("spread_limit must be None or an angle within (0, 180) degrees, got {!r}".format(spread_limit))
         self.spread_limit = spread_limit
@@ -704,6 +719,10 @@
                     n = max(n, int(math.ceil(dt / self.part_max - 1e-9)))
                     self.capped_ahead += 1
         self._top_start, self._dt_last = top, dt
+        if n > PARTS_MAX:                                    # the stall of 2026-11-08 (sub-plan 22 Task 5d): stop, do not crawl
+            raise RuntimeError("the macro step from {:.6g} s to {:.6g} s would run as {} parts, more than PARTS_MAX = {}: "
+                               "the part count has run away (with --part-count-from moves a capped hand-down feeds it; "
+                               "sub-plan 22 Task 5d)".format(t - dt, t, n, PARTS_MAX))
         q_end = np.asarray(loads.q_conv, dtype=float)
         q_start = self._q_prev if self._q_prev is not None and self._q_prev.shape == q_end.shape else q_end
         self._heating = (t - dt, dt, q_start, q_end)
@@ -842,6 +861,8 @@
             for k in range(n):
                 self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt_macro, loads, depth + 1)
             return
+        if self.part_count_from == "melting":                # the hint from the draws alone, kept parts only (Task 5d)
+            ratio = self._draw_ratio(book)
         self._ratio = max(self._ratio, ratio * dt_macro / dt)
         self.step_parts += 1
         self.pending = mismatch
@@ -926,6 +947,12 @@
             handed, _ = skin_mod.hand_down(sk, m_min, mat.T_rigid, limit, cap=self._hand_cap())
             sk.set_state(saved)
         d = np.abs(book.drawn - handed) / (mat.rho * self.columns.area)
+        return float(d.max() / (MOVE_MAX * self.columns.edge_length.mean()))
+
+    def _draw_ratio(self, book):
+        """`_move_ratio` with no hand-down: the largest distance a base would move by drawing alone, over MOVE_MAX of the
+        mean surface edge. It sets the next macro step's part count under `part_count_from="melting"` (Task 5d)."""
+        d = np.abs(book.drawn) / (self.material.rho * self.columns.area)
         return float(d.max() / (MOVE_MAX * self.columns.edge_length.mean()))
 
     def _move_base(self, book, m_min, g0):
```

- [ ] **Step 4: The command line** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -4,11 +4,12 @@
         [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
         [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
         [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
-        [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off --spread-limit DEG|off ...]
+        [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off --spread-limit DEG|off
+        --part-count-from moves|melting ...]
                                      (the layered skins, the default with --thermal fem, their flow and spraying the
                                       default with --melt on; the element model with --surface-model elements)
     python -m reentry_model resume <state> --outdir ... --name ... [--t-max S --dt-continuum S --seed N --frames-every N
-        --slurry-spray on|off --spread-limit DEG|off]
+        --slurry-spray on|off --spread-limit DEG|off --part-count-from moves|melting]
                                      (a flight saved with --save-state-at, continued; sub-plan 22 Task 5b)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
@@ -119,6 +120,11 @@
     return value
 
 
+def part_count_rule_of(args):
+    """The rule that sets a macro step's part count, `moves` unless --part-count-from says `melting` (sub-plan 22 Task 5d)."""
+    return args.part_count_from
+
+
 def spread_limit_of(args):
     """The fold angle [degrees] past which the skins spread what they hand down, or None for --spread-limit off."""
     return None if args.spread_limit == "off" else args.spread_limit
@@ -185,7 +191,7 @@
 
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
                      fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=True, flow_cycles=None,
-                     slurry_spray=False, spread_limit=None):
+                     slurry_spray=False, spread_limit=None, part_count_from="moves"):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
@@ -195,8 +201,8 @@
     with melting since sub-plan 22 Task 3b; None for a run that does not melt, which the flow does not concern), and with
     the flow on `_flowcycles<n>` where the flow sub-steps in each heat sub-step are fixed at n rather than chosen each heat
     sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3);
-    `_spread-<deg>` where what the skins hand down is spread past a fold limit (sub-plan 22 Task 5c); empty for the
-    defaults."""
+    `_spread-<deg>` where what the skins hand down is spread past a fold limit (sub-plan 22 Task 5c);
+    `_parts-melting` where the part count follows the melting, not the moves (Task 5d); empty for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -224,6 +230,8 @@
         out += "_slurryspray-on"
     if spread_limit is not None:
         out += "_spread-{:g}".format(spread_limit)
+    if part_count_from != "moves":
+        out += "_parts-{}".format(part_count_from)
     return out
 
 
@@ -285,6 +293,9 @@
     rs.add_argument("--spread-limit", type=parse_spread_limit, default=None,
                     help="[deg] the fold past which the metal handed down spreads to lower corners, or off (default: the saved "
                          "flight's, off for one saved before the option existed)")
+    rs.add_argument("--part-count-from", choices=skinbody.PART_COUNT_RULES, default=None,
+                    help="what the next macro step's part count follows: the moves (draws and hand-down) or the melting (the "
+                         "draws alone) (default: the saved flight's, moves for one saved before the option existed)")
     rs.add_argument("--save-state-at", type=float, default=None,
                     help="[s] of flight: save the resumed flight's whole state once, at the end of the first macro step ending at "
                          "or after it, as <outdir>/<name>{}; after the saved time and before --t-max (final review, 2026-11-07)".format(STATE_SUFFIX))
@@ -418,6 +429,11 @@
                     help="with --surface-model skin and --melt on: where the base as it would stand folds past this angle "
                          "[deg, within (0, 180)], metal being handed down goes to lower neighbouring corners instead, and "
                          "lands hot (sub-plan 22 Task 5c; off until the flights have tested it; default off)")
+    me.add_argument("--part-count-from", choices=skinbody.PART_COUNT_RULES, default="moves",
+                    help="with --surface-model skin and --melt on: what sets the next macro step's part count: moves, the "
+                         "move the draws and the capped hand-down call for, or melting, the draws alone; a capped "
+                         "hand-down made the count grow by a quarter a step on the spreading flights (sub-plan 22 Task 5d; "
+                         "default moves)")
     me.add_argument("--skin-flow-cycles", type=int, default=None,
                     help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
                          "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
@@ -507,7 +523,8 @@
                                      flow_cycles=args.skin_flow_cycles, surface_flow=flow_model,
                                      spray=spray_model_of(args, mat) if flow else None,
                                      slurry_spray=args.slurry_spray == "on", particles=args.particles,
-                                     spool=spool if flow else None, spread_limit=spread_limit_of(args))
+                                     spool=spool if flow else None, spread_limit=spread_limit_of(args),
+                                     part_count_from=part_count_rule_of(args))
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray_model_of(args, mat)
@@ -542,6 +559,8 @@
                         "slurry_spray": args.slurry_spray}
         if spread_limit_of(args) is not None:
             info["skin"]["spread_limit_deg"] = spread_limit_of(args)
+        if part_count_rule_of(args) != "moves":
+            info["skin"]["part_count_from"] = part_count_rule_of(args)
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
@@ -592,7 +611,7 @@
                                                                     args.skin_fine_zone,
                                                                     skin_flow_of(args) == "on" if args.melt == "on" else None,
                                                                     args.skin_flow_cycles, args.slurry_spray == "on",
-                                                                    spread_limit_of(args)))
+                                                                    spread_limit_of(args), part_count_rule_of(args)))
 
 
 def cmd_run(args, parser):
@@ -672,6 +691,8 @@
                      "sprays by the skins' spraying from their top cells")
     if spread_limit_of(args) is not None and not (args.surface_model == "skin" and args.melt == "on"):
         parser.error("--spread-limit needs --surface-model skin and --melt on: it spreads the metal the skins hand down")
+    if part_count_rule_of(args) != "moves" and not (args.surface_model == "skin" and args.melt == "on"):
+        parser.error("--part-count-from melting needs --surface-model skin and --melt on: it sets the skins' part count")
     if args.surface_model == "skin":
         for flag, given in (("--runoff off", args.runoff == "off"), ("--deep-runoff off", args.deep_runoff == "off"),
                             ("--molten-cascade off", args.molten_cascade == "off"),
@@ -890,6 +911,10 @@
         if spread_limit_of(args) is not None and not (skins and run.body.material.melts):
             parser.error("--spread-limit needs a saved flight of melting skins: it spreads the metal they hand down")
         changed["--spread-limit"] = args.spread_limit
+    if args.part_count_from is not None:
+        if args.part_count_from != "moves" and not (skins and run.body.material.melts):
+            parser.error("--part-count-from melting needs a saved flight of melting skins: it sets their part count")
+        changed["--part-count-from"] = args.part_count_from
     if args.seed is not None:
         changed["--seed"] = args.seed
     if args.save_state_at is not None:
@@ -924,6 +949,8 @@
         for key in skinbody.SPREAD_COLUMNS:
             state["row"].pop(key, None)
         state["row"].update(run.body.spread_stats())
+    if "--part-count-from" in changed and skins:
+        run.body.part_count_from = args.part_count_from
     if isinstance(run.body.source_rows, spray.SourceSpool):      # only the records after the resume are written
         run.body.source_rows = spray.SourceSpool(os.path.join(run_dir, SOURCE_SPOOL))
     if args.seed is not None:
@@ -934,19 +961,21 @@
     for key, value in (("outdir", args.outdir), ("name", args.name), ("save_state_at", changed.get("--save-state-at"))):
         setattr(saved, key, value)
     for key, flag in (("t_max", "--t-max"), ("dt_continuum", "--dt-continuum"), ("seed", "--seed"),
-                      ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray"), ("spread_limit", "--spread-limit")):
+                      ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray"), ("spread_limit", "--spread-limit"),
+                      ("part_count_from", "--part-count-from")):
         if flag in changed:
             setattr(saved, key, changed[flag])
     if state_path is not None:
         settings.save_meta = {"name": args.name, "args": dict(vars(saved))}
     limit = run.body.spread_limit if skins else None
     print("resuming {} (saved at {:.6g} s, step {}) to {:.6g} s: dt_continuum {}, seed {}, frames every {}, slurry spray {}, "
-          "spread limit {}, save state at {}; files {}".format(
+          "spread limit {}, part count from {}, save state at {}; files {}".format(
               os.path.abspath(args.state), state["time_s"], state["step"], run.sim.settings.t_max,
               settings.dt_continuum if settings.dt_continuum is not None else "none",
               args.seed if args.seed is not None else "saved random state", settings.frames_every,
               "on" if getattr(run.body, "slurry_spray", False) else "off",
               "off" if limit is None else "{:g} deg".format(limit),
+              run.body.part_count_from if skins else "moves",
               "{:g} s".format(args.save_state_at) if args.save_state_at is not None else "no",
               os.path.join(args.outdir, args.name)), flush=True)
     history = run.run(resume=state)
@@ -960,6 +989,8 @@
         thermal_info["save_state_at_s"] = args.save_state_at
     if skins and run.body.spread_limit is not None:
         thermal_info["spread_limit_deg"] = run.body.spread_limit
+    if skins and run.body.part_count_from == "melting":
+        thermal_info["part_count_from"] = "melting"
     resumed = {"state": os.path.abspath(args.state), "step": state["step"], "time_s": state["time_s"],
                "run_name": state["meta"]["name"], "changed": changed}
     os.makedirs(args.outdir, exist_ok=True)
```

- [ ] **Step 5: Run the tests to see them pass,** then the unit tier: 608 passed, 1 skipped (Task 5c's 592 and these 16; 22.8 minutes across four
  processes, 2026-11-08).

- [ ] **Step 6: Commit.**

```bash
git add reentry_model/skinbody.py reentry_model/cli.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "Skins: the part count can follow the melting, not the moves (--part-count-from moves|melting on run and resume, moves by default), and a macro step that would start with more than 1 000 parts stops the run; sub-plan 22 Task 5d"
```

---

### Task 5e: The hand-down switched off

**Files:**
- Modify: `reentry_model/skinbody.py` (the class default and `SkinBody(..., hand_down=True)`, `_move_base`, `_move_ratio`)
- Modify: `reentry_model/cli.py` (`--hand-down on|off` on `run` and `resume`, `hand_down_of`, the name suffix, the records, the
  refusals, the resume's one line)
- Test: `tests/test_reentry_model_skinbody.py` (6 tests), `tests/test_reentry_model_cli.py` (3 tests, 7 cases),
  `tests/test_reentry_model_restart.py` (1 test)

**Interfaces:**
- Consumes: Task 5d's wiring of an option on `run` and `resume`; `skin_mod.hand_down`, `skin_mod.handable`.
- Produces: `SkinBody(..., hand_down=True)`, also a class default for states pickled before it; `ValueError` if `hand_down` is False
  and `spread_limit` is not None. `--hand-down on|off` on `run` (default `on`) and `resume` (default: the state's); exit 2 with
  `off` unless melting skins, and with `off` and a spread limit, given or carried by the resumed state (and a limit given to a
  resume of a state saved with `off`), with a message saying spreading moves only metal being handed down (Asha, 2026-11-08) and to
  give `--spread-limit off`. The run's name suffix `_handdown-off`; the run JSON's `settings.skin.hand_down = "off"` only when off,
  the resume JSON's `settings.hand_down = "off"` only when the body's is off and `changed["--hand-down"]` when given.

**What was decided** (design `docs/superpowers/specs/2026-11-08-hand-down-and-part-count-design.md`; Asha, 2026-11-08). Re-frozen
metal can heat again and flow on or be sprayed, so it may stay in the skins instead of passing to the 3D model; spreading is off
with it, because it moves only metal being handed down. With `hand_down` False `_move_base` does not call `skin_mod.hand_down`:
`handed` and `handed_H` are zeros, `coldest` and `warmest` nan, `uncovered` nan, and `held = able`, so `held_mass` and
`held_depth_max` record every handable kilogram; the move, the reconciliation and the fault count read the zero hand-down.
`_move_ratio` uses a zero hand-down too. Draws and absorption are unchanged, so the base moves only inward (and by the neutral
smoothing).

**Rulings** (the controller's, 2026-11-08): the test that no metal crosses the base outward asserts that `handed_mass` stays 0 and the
3D model's enclosed volume never grows from one macro step to the next, in a set-up (Task 5d's held-back deposit) where with the
hand-down on it does grow; and `held_mass` and `held_depth_max` record every handable kilogram with the hand-down off, which is not
a fault and adds no column.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1743,3 +1743,101 @@
 def test_an_unknown_part_count_rule_is_refused(coarse_sphere_mesh):
     with pytest.raises(ValueError, match="part_count_from"):
         make(coarse_sphere_mesh, part_count_from="draws")
+
+
+# -- the hand-down switched off (sub-plan 22 Task 5e; Asha 2026-11-08) ------------------------------------------------------
+def fly_macro_steps(b, steps=2, dt=0.05, loads=None):
+    """`steps` macro steps of `dt` s (not longer than the melting part length), with no heat unless `loads` is given; the 3D
+    model's enclosed volume after each [m3]."""
+    loads = loads if loads is not None else heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    volumes = []
+    for k in range(steps):
+        b.advance(dt * (k + 1), dt, loads)
+        volumes.append(b.mesh.volume())
+    return volumes
+
+
+def test_with_the_hand_down_off_no_metal_crosses_the_base_outward(coarse_sphere_mesh):
+    """Asha's decision of 2026-11-08 (sub-plan 22 Task 5e): with `hand_down=False` nothing is handed down, so no metal is
+    added to the 3D model. On the deposit of `held_back_body` (nine move allowances of metal at the node's temperature,
+    no heat, nothing drawn), two 0.05 s macro steps: the default hands down and the enclosed volume of the 3D model grows
+    in each; with the hand-down off `handed_mass` stays 0 and the volume does not grow from one step to the next."""
+    on, _ = held_back_body(coarse_sphere_mesh)
+    off, _ = held_back_body(coarse_sphere_mesh, hand_down=False)
+    start = off.mesh.volume()
+    assert on.mesh.volume() == start
+    v_on, v_off = fly_macro_steps(on), fly_macro_steps(off)
+    assert on.handed_mass > 0.0 and v_on[0] > start and v_on[1] > v_on[0]
+    assert off.handed_mass == 0.0 and v_off[0] <= start * (1.0 + 1e-12) and v_off[1] <= v_off[0] * (1.0 + 1e-12)
+
+
+def test_with_the_hand_down_off_the_rigid_metal_stays_in_its_skin_and_is_counted_as_held(coarse_sphere_mesh):
+    """The metal the default hands down stays in its skin with the hand-down off, so the skins' mass above their minimum is at
+    least the default's after the same steps, and `held_mass` is every kilogram `skin.handable` finds, with no limit or cap
+    (it is found on the skins as they stand before the first step: nothing is drawn or heated, so they do not change)."""
+    on, extra = held_back_body(coarse_sphere_mesh)
+    off, _ = held_back_body(coarse_sphere_mesh, hand_down=False)
+    m_min = skin.minimum_mass(off.material, off.columns.area, off.settings)
+    able = skin.handable(off.skins, m_min, off.material.T_rigid)
+    area = off.columns.area.copy()               # the move changes the areas, and the held depth reads those before it
+    assert able.sum() > 0.0
+    fly_macro_steps(on, steps=1)
+    fly_macro_steps(off, steps=1)
+    above = lambda b: b.skins.mass() - m_min.sum()
+    assert above(off) >= above(on) and above(off) > above(on) + 0.5 * on.handed_mass
+    assert off.held_mass == pytest.approx(able.sum(), rel=1e-9) and off.held_mass > on.held_mass
+    assert off.held_depth_max == pytest.approx((able / (off.material.rho * area)).max(), rel=1e-9)
+
+
+def test_with_the_hand_down_off_the_books_close(coarse_sphere_mesh):
+    """Skins that only grow still close the books (mass 1e-12, energy 1e-8, momentum) on the held-back deposit flown for two
+    0.05 s macro steps, the second heated at 3 MW/m2 so that draws and the move run too."""
+    b, _ = held_back_body(coarse_sphere_mesh, hand_down=False)
+    M0 = b.mass()
+    b.energy0 = E0 = b.energy()                  # the deposit was added after the body was built: its heat is the start's
+    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    b.advance(0.05, 0.05, zero)
+    b.advance(0.1, 0.05, front_loads(b, 3e6))
+    assert b.handed_mass == 0.0 and b.drawn_mass >= 0.0
+    assert b.mass() == pytest.approx(M0, rel=1e-12) and b.jump_faults == 0
+    assert abs(b.energy_balance_residual()) < 1e-8 and abs(b.momentum_residual()) < 1e-8
+    assert b.energy() == pytest.approx(E0 + b.absorbed_heat - b.removed_enthalpy - b.released_enthalpy, rel=1e-8)
+
+
+def test_the_hand_down_on_by_default_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
+    """`hand_down=True` and the argument left out are one flight: every column of every row and the final mass and energy
+    equal to the last bit on the held-back deposit, which hands down (sub-plan 22 Task 5e)."""
+    runs = []
+    for kw in ({}, {"hand_down": True}):
+        b, _ = held_back_body(coarse_sphere_mesh, **kw)
+        zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+        rows = []
+        for k in range(2):
+            b.advance(0.05 * (k + 1), 0.05, zero)
+            rows.append((b.step_parts, timed_out(b.skin_stats())))
+        runs.append((b, rows))
+    assert same_rows(runs[0][1], runs[1][1]) and runs[0][0].handed_mass > 0.0 and runs[0][0].hand_down is True
+    assert runs[0][0].mass() == runs[1][0].mass() and runs[0][0].energy() == runs[1][0].energy()
+
+
+def test_the_hand_down_off_refuses_a_spread_limit(coarse_sphere_mesh):
+    """Spreading moves only metal being handed down (Asha, 2026-11-08), so a body with neither is refused."""
+    with pytest.raises(ValueError, match="hand_down"):
+        make(coarse_sphere_mesh, hand_down=False, spread_limit=60.0)
+    assert make(coarse_sphere_mesh, hand_down=False).spread_limit is None
+    assert make(coarse_sphere_mesh, spread_limit=60.0).hand_down is True
+
+
+def test_a_body_saved_without_the_hand_down_switch_flies_on(coarse_sphere_mesh):
+    """The states saved before the switch existed have no `hand_down` on the instance: the class's default gives them the
+    hand-down on, they fly on as before, and a resume can switch it off (sub-plan 22 Task 5e)."""
+    b, _ = held_back_body(coarse_sphere_mesh)
+    del vars(b)["hand_down"]
+    assert "hand_down" not in vars(b) and b.hand_down is True
+    fly_macro_steps(b, steps=1)
+    assert b.handed_mass > 0.0
+    handed = b.handed_mass
+    b.hand_down = False
+    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    b.advance(0.1, 0.05, zero)
+    assert b.handed_mass == handed
```

  in `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -772,3 +772,51 @@
         with pytest.raises(_Built) as got:
             cli.main(argv)
         assert got.value.args[0]["part_count_from"] == rule
+
+
+def test_the_hand_down_switch_reads_on_or_off_and_marks_the_run_s_name_only_for_off():
+    """--hand-down (sub-plan 22 Task 5e) is `on` by default, so that every existing run keeps its name and its record; `off`
+    ends a skin run's settings in _handdown-off, after the part count's and before the seed's."""
+    p = cli.build_parser()
+    base = "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-sesam_melt-girin"
+    name = lambda argv: cli.run_name_of(p.parse_args(argv), 0.1, 7500.0, 77500.133, "us76", cli.continuum_step(p.parse_args(argv)))
+    melt = THERMAL + ["--melt", "on"]
+    assert p.parse_args(melt).hand_down == "on" and cli.hand_down_of(p.parse_args(melt)) is True
+    assert cli.hand_down_of(p.parse_args(melt + ["--hand-down", "off"])) is False
+    assert name(melt) == name(melt + ["--hand-down", "on"]) == base
+    assert name(melt + ["--hand-down", "off"]) == base + "_handdown-off"
+    assert (name(melt + ["--part-count-from", "melting", "--hand-down", "off", "--seed", "7"])
+            == base + "_parts-melting_handdown-off_seed-7")
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, hand_down=False) == "_handdown-off"
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True) == ""
+
+
+@pytest.mark.parametrize("argv", [
+    THERMAL + ["--melt", "on", "--hand-down", "maybe"],                                # no such setting
+    THERMAL + ["--hand-down", "off"],                                                  # nothing melts, so nothing is handed down
+    FEM + ["--melt", "on", "--hand-down", "off"],                                      # the element model has no skins
+    BASE + ["--hand-down", "off"],                                                     # no thermal model at all
+    THERMAL + ["--melt", "on", "--hand-down", "off", "--spread-limit", "60"],          # spreading moves only handed-down metal
+])
+def test_a_bad_hand_down_switch_or_off_without_melting_skins_or_with_a_spread_limit_exits_2(argv, tmp_path, capsys):
+    with pytest.raises(SystemExit) as exc:
+        cli.main(argv + ["--t-max", "1", "--outdir", str(tmp_path)])
+    assert exc.value.code == 2
+    if "60" in argv:
+        err = capsys.readouterr().err
+        assert "spreading moves only metal being handed down" in err and "--spread-limit off" in err
+
+
+def test_the_command_line_gives_the_skins_their_hand_down_switch_and_records_it_only_for_off(tmp_path, monkeypatch):
+    """--hand-down reaches the skin body as True or False (sub-plan 22 Task 5e). The body is not built: its arguments are."""
+    from reentry_model import skinbody
+
+    def built(*args, **kwargs):
+        raise _Built(kwargs)
+    monkeypatch.setattr(skinbody, "SkinBody", built)
+    melt = THERMAL + ["--melt", "on", "--fine-band", "0", "--t-max", "1", "--outdir", str(tmp_path)]
+    for argv, on in ((melt, True), (melt + ["--hand-down", "on"], True), (melt + ["--hand-down", "off"], False),
+                     (melt + ["--skin-flow", "off", "--hand-down", "off"], False)):
+        with pytest.raises(_Built) as got:
+            cli.main(argv)
+        assert got.value.args[0]["hand_down"] is on
```

  and in `tests/test_reentry_model_restart.py`:

```diff
--- a/tests/test_reentry_model_restart.py
+++ b/tests/test_reentry_model_restart.py
@@ -328,3 +328,38 @@
         cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--part-count-from", "melting"])
     assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
     assert cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--part-count-from", "moves"]) == 0
+
+
+def test_resuming_with_the_hand_down_off_from_a_state_with_a_spread_limit_needs_the_limit_off(tmp_path, capsys):
+    """Spreading moves only metal being handed down (Asha, 2026-11-08; sub-plan 22 Task 5e), so `--hand-down off` on a resume of
+    a state saved with a spread limit (the spreading states carry 60) exits 2 unless `--spread-limit off` is given too, and
+    flies with it: recorded in the JSON's settings and in what was changed, named in the line before the flight, and kept by
+    a resume of the state it saves. The reverse, a limit given to a state saved with the hand-down off, exits 2 as well."""
+    pytest.importorskip("cantera")
+    limited = str(tmp_path / "limited")
+    assert cli.main(FLIGHT + ["--t-max", "0.075", "--outdir", limited, "--name", "f", "--save-state-at", "0.05",
+                              "--spread-limit", "60"]) == 0
+    state, out = os.path.join(limited, "f.state.pkl"), str(tmp_path / "o")
+    capsys.readouterr()
+    for options in (["--hand-down", "off"], ["--hand-down", "off", "--spread-limit", "60"]):
+        with pytest.raises(SystemExit) as exc:
+            cli.main(["resume", state, "--outdir", out, "--name", "r"] + options)
+        assert exc.value.code == 2 and not os.path.exists(out)
+        err = capsys.readouterr().err
+        assert "spreading moves only metal being handed down" in err and "--spread-limit off" in err
+    assert cli.main(["resume", state, "--outdir", out, "--name", "r", "--hand-down", "off", "--spread-limit", "off",
+                     "--save-state-at", "0.0625"]) == 0
+    assert "hand-down off" in capsys.readouterr().out
+    doc = json.load(open(os.path.join(out, "r.json")))
+    assert doc["resumed_from"]["changed"] == {"--hand-down": "off", "--spread-limit": "off", "--save-state-at": 0.0625}
+    assert doc["settings"]["hand_down"] == "off" and "spread_limit_deg" not in doc["settings"]
+    kept = os.path.join(out, "r.state.pkl")
+    assert restart.load(kept)["meta"]["args"]["hand_down"] == "off"
+    with pytest.raises(SystemExit) as exc:                  # a limit on a flight that hands nothing down
+        cli.main(["resume", kept, "--outdir", out, "--name", "again", "--spread-limit", "60"])
+    assert exc.value.code == 2 and "--hand-down on" in capsys.readouterr().err
+    assert cli.main(["resume", kept, "--outdir", out, "--name", "same"]) == 0
+    assert json.load(open(os.path.join(out, "same.json")))["settings"]["hand_down"] == "off"
+    assert cli.main(["resume", kept, "--outdir", out, "--name", "back", "--hand-down", "on"]) == 0
+    doc = json.load(open(os.path.join(out, "back.json")))
+    assert "hand_down" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--hand-down": "on"}
```

- [ ] **Step 2: Run them to see them fail.** Expected: `TypeError` on `hand_down` for the six body tests, `AttributeError` on
  `hand_down_of` and unrecognised-argument exits for the command line's.

- [ ] **Step 3: The body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -560,10 +560,13 @@
     `spread_limit` [deg, within (0, 180)]: the fold past which a part's hand-down is spread to lower neighbouring corners
     (`spreading`; sub-plan 22 Task 5c), None for no spreading.
     `part_count_from`: `moves` (the default) or `melting`, what the next macro step's part count follows
-    (`PART_COUNT_RULES`; sub-plan 22 Task 5d)."""
+    (`PART_COUNT_RULES`; sub-plan 22 Task 5d).
+    `hand_down` [bool]: whether re-frozen metal at the bottom of a skin passes to the 3D model (the default), or stays in
+    the skins (sub-plan 22 Task 5e); False is refused with a `spread_limit`, which moves only metal being handed down."""
     # class defaults as well as the instance's: the flights saved on 2026-11-06 were pickled before spreading existed
     spread_limit = None                  # deg
     part_count_from = "moves"            # the part count's rule (sub-plan 22 Task 5d): the states saved before it have none
+    hand_down = True                     # re-frozen metal passes to the 3D model (Task 5e): the states saved before it have it on
     spread_mass = spread_unplaced = 0.0  # kg: moved between columns since time zero; still wanting to move at the last move
     spread_over = 0                      # edges over the limit on the last move's first estimate
     spread_jump = 0.0                    # K: the largest landing of metal over its node's temperature this macro step
@@ -572,7 +575,11 @@
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
                  surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True, spool=None,
-                 spread_limit=None, part_count_from="moves"):
+                 spread_limit=None, part_count_from="moves", hand_down=True):
+        if not hand_down and spread_limit is not None:
+            raise ValueError("spreading moves only metal being handed down (Asha, 2026-11-08): hand_down=False needs "
+                             "spread_limit=None, got {!r}".format(spread_limit))
+        self.hand_down = bool(hand_down)
         if part_count_from not in PART_COUNT_RULES:
             raise ValueError("part_count_from must be one of {}, got {!r}".format(PART_COUNT_RULES, part_count_from))
         self.part_count_from = part_count_from
@@ -941,7 +948,7 @@
         The hand-down is found on a copy of the skins."""
         mat, sk = self.material, self.skins
         handed = np.zeros(sk.n)
-        if mat.melts:
+        if mat.melts and self.hand_down:                     # with the hand-down off, a zero one (Task 5e)
             saved = sk.state()
             limit = self.solver.temperature()[self.moving.corners] + skin_mod.JUMP_TOLERANCE
             handed, _ = skin_mod.hand_down(sk, m_min, mat.T_rigid, limit, cap=self._hand_cap())
@@ -972,10 +979,17 @@
         T_node = self.solver.temperature()[mm.corners]
         none = np.full(sk.n, np.nan)
         able = skin_mod.handable(sk, m_min, mat.T_rigid) if mat.melts else np.zeros(sk.n)     # with no limit or cap
+        # The hand-down switched off (Asha's question and decision, 2026-11-08; sub-plan 22 Task 5e): re-frozen metal can heat
+        # again and flow on or be sprayed, so it may stay in the skins instead of passing to the 3D model, and spreading is
+        # off with it, because it moves only metal being handed down (the constructor refuses the two together). Nothing is
+        # handed, so `skin_mod.hand_down` is not called: the hand-down and its enthalpy are zeros, there is no coldest or
+        # warmest cell handed and no cell uncovered, and every handable kilogram is held (`held = able`, below). The draws
+        # and the absorption are as they were, so a base moves only inward, by them and the neutral smoothing.
         handed, handed_H, coldest, warmest = (
             skin_mod.hand_down(sk, m_min, mat.T_rigid, T_node + skin_mod.JUMP_TOLERANCE, span=True, cap=self._hand_cap())
-            if mat.melts else (np.zeros(sk.n), np.zeros(sk.n), none, none))
-        uncovered = np.where(handed > 0.0, sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)], np.nan)  # as it left it
+            if mat.melts and self.hand_down else (np.zeros(sk.n), np.zeros(sk.n), none, none))
+        uncovered = (np.where(handed > 0.0, sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)], np.nan)  # as it left it
+                     if self.hand_down else none)
         held = np.maximum(able - handed, 0.0)                 # what the tolerance and the cap keep in the skins
         self.held_mass, self.held_depth_max = float(held.sum()), float((held / (mat.rho * col.area)).max())
         handed_s, landed_H, excess = handed, handed_H, None                       # what crosses each base; the same arrays
```

- [ ] **Step 4: The command line** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -5,11 +5,12 @@
         [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
         [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
         [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off --spread-limit DEG|off
-        --part-count-from moves|melting ...]
+        --part-count-from moves|melting --hand-down on|off ...]
                                      (the layered skins, the default with --thermal fem, their flow and spraying the
                                       default with --melt on; the element model with --surface-model elements)
     python -m reentry_model resume <state> --outdir ... --name ... [--t-max S --dt-continuum S --seed N --frames-every N
-        --slurry-spray on|off --spread-limit DEG|off --part-count-from moves|melting]
+        --slurry-spray on|off --spread-limit DEG|off --part-count-from moves|melting
+        --hand-down on|off]
                                      (a flight saved with --save-state-at, continued; sub-plan 22 Task 5b)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
@@ -125,6 +126,15 @@
     return args.part_count_from
 
 
+NO_SPREAD_WITHOUT_HAND_DOWN = ("--hand-down off cannot go with a spread limit: spreading moves only metal being handed down "
+                               "(Asha, 2026-11-08); give --spread-limit off")
+
+
+def hand_down_of(args):
+    """Whether the skins pass their re-frozen metal to the 3D model, True unless --hand-down says `off` (sub-plan 22 Task 5e)."""
+    return args.hand_down == "on"
+
+
 def spread_limit_of(args):
     """The fold angle [degrees] past which the skins spread what they hand down, or None for --spread-limit off."""
     return None if args.spread_limit == "off" else args.spread_limit
@@ -191,7 +201,8 @@
 
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
                      fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=True, flow_cycles=None,
-                     slurry_spray=False, spread_limit=None, part_count_from="moves"):
+                     slurry_spray=False, spread_limit=None, part_count_from="moves",
+                     hand_down=True):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
@@ -202,7 +213,8 @@
     the flow on `_flowcycles<n>` where the flow sub-steps in each heat sub-step are fixed at n rather than chosen each heat
     sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3);
     `_spread-<deg>` where what the skins hand down is spread past a fold limit (sub-plan 22 Task 5c);
-    `_parts-melting` where the part count follows the melting, not the moves (Task 5d); empty for the defaults."""
+    `_parts-melting` where the part count follows the melting, not the moves (Task 5d); `_handdown-off` where re-frozen
+    metal stays in the skins (Task 5e); empty for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -232,6 +244,8 @@
         out += "_spread-{:g}".format(spread_limit)
     if part_count_from != "moves":
         out += "_parts-{}".format(part_count_from)
+    if not hand_down:
+        out += "_handdown-off"
     return out
 
 
@@ -296,6 +310,9 @@
     rs.add_argument("--part-count-from", choices=skinbody.PART_COUNT_RULES, default=None,
                     help="what the next macro step's part count follows: the moves (draws and hand-down) or the melting (the "
                          "draws alone) (default: the saved flight's, moves for one saved before the option existed)")
+    rs.add_argument("--hand-down", choices=("on", "off"), default=None,
+                    help="re-frozen metal passes to the 3D model, or stays in the skins; off needs the spread limit off, given "
+                         "or carried by the saved flight (default: the saved flight's, on for one saved before the option existed)")
     rs.add_argument("--save-state-at", type=float, default=None,
                     help="[s] of flight: save the resumed flight's whole state once, at the end of the first macro step ending at "
                          "or after it, as <outdir>/<name>{}; after the saved time and before --t-max (final review, 2026-11-07)".format(STATE_SUFFIX))
@@ -434,6 +451,11 @@
                          "move the draws and the capped hand-down call for, or melting, the draws alone; a capped "
                          "hand-down made the count grow by a quarter a step on the spreading flights (sub-plan 22 Task 5d; "
                          "default moves)")
+    me.add_argument("--hand-down", choices=("on", "off"), default="on",
+                    help="with --surface-model skin and --melt on: re-frozen metal at the bottom of a skin passes to the 3D "
+                         "model (on) or stays in the skins, which can heat it again, let it flow on or spray it (off; needs "
+                         "--spread-limit off, spreading moves only metal being handed down; Asha, 2026-11-08, sub-plan 22 "
+                         "Task 5e; default on)")
     me.add_argument("--skin-flow-cycles", type=int, default=None,
                     help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
                          "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
@@ -524,7 +546,7 @@
                                      spray=spray_model_of(args, mat) if flow else None,
                                      slurry_spray=args.slurry_spray == "on", particles=args.particles,
                                      spool=spool if flow else None, spread_limit=spread_limit_of(args),
-                                     part_count_from=part_count_rule_of(args))
+                                     part_count_from=part_count_rule_of(args), hand_down=hand_down_of(args))
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray_model_of(args, mat)
@@ -561,6 +583,8 @@
             info["skin"]["spread_limit_deg"] = spread_limit_of(args)
         if part_count_rule_of(args) != "moves":
             info["skin"]["part_count_from"] = part_count_rule_of(args)
+        if not hand_down_of(args):
+            info["skin"]["hand_down"] = "off"
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
@@ -611,7 +635,8 @@
                                                                     args.skin_fine_zone,
                                                                     skin_flow_of(args) == "on" if args.melt == "on" else None,
                                                                     args.skin_flow_cycles, args.slurry_spray == "on",
-                                                                    spread_limit_of(args), part_count_rule_of(args)))
+                                                                    spread_limit_of(args), part_count_rule_of(args),
+                                                                    hand_down_of(args)))
 
 
 def cmd_run(args, parser):
@@ -693,6 +718,10 @@
         parser.error("--spread-limit needs --surface-model skin and --melt on: it spreads the metal the skins hand down")
     if part_count_rule_of(args) != "moves" and not (args.surface_model == "skin" and args.melt == "on"):
         parser.error("--part-count-from melting needs --surface-model skin and --melt on: it sets the skins' part count")
+    if not hand_down_of(args) and not (args.surface_model == "skin" and args.melt == "on"):
+        parser.error("--hand-down off needs --surface-model skin and --melt on: it keeps the skins' re-frozen metal in them")
+    if not hand_down_of(args) and spread_limit_of(args) is not None:
+        parser.error(NO_SPREAD_WITHOUT_HAND_DOWN)
     if args.surface_model == "skin":
         for flag, given in (("--runoff off", args.runoff == "off"), ("--deep-runoff off", args.deep_runoff == "off"),
                             ("--molten-cascade off", args.molten_cascade == "off"),
@@ -915,6 +944,16 @@
         if args.part_count_from != "moves" and not (skins and run.body.material.melts):
             parser.error("--part-count-from melting needs a saved flight of melting skins: it sets their part count")
         changed["--part-count-from"] = args.part_count_from
+    if args.hand_down is not None:
+        if args.hand_down == "off" and not (skins and run.body.material.melts):
+            parser.error("--hand-down off needs a saved flight of melting skins: it keeps their re-frozen metal in them")
+        changed["--hand-down"] = args.hand_down
+    # the spread limit and the hand-down as they will stand, each as given here or else the saved flight's
+    if skins:
+        will_spread = spread_limit_of(args) if args.spread_limit is not None else run.body.spread_limit
+        will_hand = args.hand_down == "on" if args.hand_down is not None else run.body.hand_down
+        if will_spread is not None and not will_hand:
+            parser.error(NO_SPREAD_WITHOUT_HAND_DOWN + (" (or --hand-down on)" if args.hand_down is None else ""))
     if args.seed is not None:
         changed["--seed"] = args.seed
     if args.save_state_at is not None:
@@ -951,6 +990,8 @@
         state["row"].update(run.body.spread_stats())
     if "--part-count-from" in changed and skins:
         run.body.part_count_from = args.part_count_from
+    if "--hand-down" in changed and skins:
+        run.body.hand_down = args.hand_down == "on"
     if isinstance(run.body.source_rows, spray.SourceSpool):      # only the records after the resume are written
         run.body.source_rows = spray.SourceSpool(os.path.join(run_dir, SOURCE_SPOOL))
     if args.seed is not None:
@@ -962,20 +1003,21 @@
         setattr(saved, key, value)
     for key, flag in (("t_max", "--t-max"), ("dt_continuum", "--dt-continuum"), ("seed", "--seed"),
                       ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray"), ("spread_limit", "--spread-limit"),
-                      ("part_count_from", "--part-count-from")):
+                      ("part_count_from", "--part-count-from"), ("hand_down", "--hand-down")):
         if flag in changed:
             setattr(saved, key, changed[flag])
     if state_path is not None:
         settings.save_meta = {"name": args.name, "args": dict(vars(saved))}
     limit = run.body.spread_limit if skins else None
     print("resuming {} (saved at {:.6g} s, step {}) to {:.6g} s: dt_continuum {}, seed {}, frames every {}, slurry spray {}, "
-          "spread limit {}, part count from {}, save state at {}; files {}".format(
+          "spread limit {}, part count from {}, hand-down {}, save state at {}; files {}".format(
               os.path.abspath(args.state), state["time_s"], state["step"], run.sim.settings.t_max,
               settings.dt_continuum if settings.dt_continuum is not None else "none",
               args.seed if args.seed is not None else "saved random state", settings.frames_every,
               "on" if getattr(run.body, "slurry_spray", False) else "off",
               "off" if limit is None else "{:g} deg".format(limit),
               run.body.part_count_from if skins else "moves",
+              "off" if skins and not run.body.hand_down else "on",
               "{:g} s".format(args.save_state_at) if args.save_state_at is not None else "no",
               os.path.join(args.outdir, args.name)), flush=True)
     history = run.run(resume=state)
@@ -991,6 +1033,8 @@
         thermal_info["spread_limit_deg"] = run.body.spread_limit
     if skins and run.body.part_count_from == "melting":
         thermal_info["part_count_from"] = "melting"
+    if skins and not run.body.hand_down:
+        thermal_info["hand_down"] = "off"
     resumed = {"state": os.path.abspath(args.state), "step": state["step"], "time_s": state["time_s"],
                "run_name": state["meta"]["name"], "changed": changed}
     os.makedirs(args.outdir, exist_ok=True)
```

- [ ] **Step 5: Run the tests to see them pass,** then the unit tier: 622 passed, 1 skipped (Task 5d's 608 and these 14; 17.3 minutes across four processes, 2026-11-08).

- [ ] **Step 6: Commit.**

```bash
git add reentry_model/skinbody.py reentry_model/cli.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "Skins: re-frozen metal can stay in the skins (--hand-down on|off on run and resume, on by default; off refuses a spread limit); sub-plan 22 Task 5e"
```

---

### Task 5f: Sideways heat between cells that touch

**Files:**
- Modify: `reentry_model/skin.py` (`SIDEWAYS_PAIRINGS`, `_height_geometry` and `SidewaysByHeight`, beside `Sideways`, which is
  unchanged)
- Modify: `reentry_model/skinbody.py` (the class default and `SkinBody(..., sideways_pairing="depth")`, the `side` lambda in
  `_advance_part`)
- Modify: `reentry_model/cli.py` (`--sideways-heat depth|height` on `run` and `resume`, the name suffix, the records, the refusals,
  the resume's one line)
- Test: `tests/test_reentry_model_skin.py` (7 tests), `tests/test_reentry_model_skinbody.py` (3 tests),
  `tests/test_reentry_model_cli.py` (3 tests, 6 cases), `tests/test_reentry_model_restart.py` (1 test)

**Interfaces:**
- Consumes: `skin.Sideways(sk, columns, previous)` and its contract with `_coupled_solve` and `_conduct`: `heat(x)` [W per cell,
  (n, K)] linear in x, antisymmetric edge by edge (it sums to zero to round-off), symmetric and negative semi-definite, zero in
  the cells beyond a column's count; `columns.edges`, `columns.conduct_weight`, `columns.points` (the base corners as they stand)
  and `columns.normal`; Task 5d's wiring of an option on `run` and `resume`.
- Produces: `skin.SIDEWAYS_PAIRINGS = ("depth", "height")`; `skin.SidewaysByHeight(sk, columns, previous=None)`, a map with the same
  contract; `SkinBody(..., sideways_pairing="depth")`, also a class default for states pickled before it, `ValueError` for another
  value. `--sideways-heat depth|height` on `run` (default `depth`) and `resume` (default: the state's); exit 2 with `height` unless
  `--surface-model skin` (run) or a saved flight of skins (resume), melting not required, the pairing applying to any skins. The
  run's name suffix `_side-height`; the run JSON's `settings.skin.sideways_heat = "height"` only when given, the resume JSON's
  `settings.sideways_heat = "height"` only when the body's pairing is `height` and `changed["--sideways-heat"]` when given; the
  resume's one line names the sideways heat.

**What was decided** (design `docs/superpowers/specs/2026-11-08-hand-down-and-part-count-design.md`, item 3; Asha's idea of
2026-11-08, first written in the side session's draft `docs/superpowers/specs/2026-11-08-refrozen-metal-in-skins-design.md`
section 3; Asha then asked for the heat to be paired correctly on the comparison flights). Spec 6.3 pairs each cell with the
neighbour's cells at the same depth below the top. Where the tops differ this pairs cells that do not touch: at the ring's outer
edge the top steps about 1.4 mm across a 1.5 mm triangle, so with the deposit kept in the skins (Task 5e) a cell beside the gas
would be paired with metal 1.4 mm lower. With `height` every cell's range is its height above its column's base,
`[H - z_bottom, H - z_top]` (z from the top, H the column's thickness); on each edge (i, j) the neighbour's ranges are shifted into
i's frame by `s_ij = (p_j - p_i) . n_ij`, p the base corners and n_ij the unit mean of the two columns' directions, `s_ji = -s_ij`;
cells a of i and b of j exchange `c_ab (x_b - x_a)` with `c_ab = 0.5 (k_a + k_b) w_ij l_ab`, k the cells' conductivities at the
sub-step's start, w_ij the edge's conduction weight and l_ab the overlap of the two ranges. `c_ab` is symmetric in a and b, so the
map is symmetric and negative semi-definite by construction and the implicit sub-step takes it as it takes the depth map. **What
it brings out, stated:** the part of a taller column above its neighbour's top is a side wall facing the gas; it exchanges no heat
with the neighbour, and the model does not heat or cool it from the gas. The flow's sideways transport (spec 6.5, matched by
fraction of the flowing depth) is not changed. The geometry is rebuilt whenever the cells' thicknesses or the base corners change
(a copy of the corners is kept, so a base moved in place is seen), and the overlaps are found for every edge at once: each edge's
boundaries of both columns are sorted downward in one pass, and between two consecutive boundaries lies one cell of each column, or
none above a top or below a base. The default stays `depth`, and `Sideways` is not touched, so every run without the option is as
before to the last bit.

**Measured** (2026-11-08, each once, on the laptop on mains power; a finding, not a gate, the production cost being the comparison
flights'). On the 4 mm test sphere with the dev cells (2 254 columns, 6 756 edges, 33-49 cells, the columns one to three times their
minimum deep on bases stepped by up to 0.3 mm, every cell at its own temperature): 401 969 pairs, 59.5 an edge, 13 MB; a map's build
0.019 s against the depth map's 0.0095, a `heat` 0.0025 s against 0.0026, one 10 ms sub-step under 3 MW/m2 on the front 0.148 s against
0.158 in 11 Newton steps and 10 conjugate-gradient iterations each. On a 2 mm surface of the 100 mm sphere with 10 um cells (9 519
columns, 28 551 edges, 59-75 cells): 2.66 million pairs, 93 an edge, 85 MB; a build 0.15 s against 0.09 (its peak 0.37 GB against
0.09), with the geometry reused 0.028 s against 0.034, a `heat` 0.020 s against 0.027; one 10 ms sub-step 3.6-3.8 s against 3.2-3.3 over
three tries, the same 11 Newton steps (26 conjugate-gradient iterations against 29) but 49 evaluations of the residual against 32 in
the line search, which with the build is the whole difference. On level columns of cells on the depth grid the two maps agree to
2.0e-15 of the largest heat; moving one boundary by half a bin makes them differ by 12 %. The books test's flowing body closes its
mass, energy and momentum with the height map as with the depth map; the resumed test flight's rows differ from the depth pairing's
from the first step on. The unit tier: 639 passed, 1 skipped in 15.3 minutes across four processes.

**Rulings** (2026-11-08; what, why, cost if wrong):
1. *The offset and the conductance are the controller's* (`s_ij` along the unit mean of the two directions, `c_ab = 0.5 (k_a + k_b) w l`),
   final in the brief. Why: symmetric by construction, so the implicit solve's contract holds. Cost: another offset changes the heat at a
   base's step only.
2. *The equality with the depth map is tested where it must hold* (the controller's): level tops and bases, cells on the depth grid;
   elsewhere the depth map averages each bin over the cells in it and the two legitimately differ.
3. *A second class beside `Sideways`, not an argument on it*, and the body picks the class. Why: `Sideways` stays the same code, so the
   default is the same flight to the last bit without having to show it branch by branch. Cost: the two classes share no code but the
   contract.
4. *A height map reuses only a height map's geometry* (`type(previous) is type(self)`). In a body a pass's maps are all of one pairing,
   so this only guards a caller that mixes them. Cost: none.
5. *The conductivity is each cell's own at the sub-step's start*, `k(T)` cell by cell, where the depth map takes it at each bin's mean
   temperature; on the depth grid the two are the same. Cost: none found.
6. *The hand-built test columns are a namespace* holding the edges, weights, incidence, corners and directions the maps read, not a
   triangulation. Why: the corners of one triangle lie in one plane, so an offset along the mean direction is zero on any triangle, and
   a two-column case with one edge is computed by hand. Cost: the cases do not exercise `Columns`; the sphere tests do.
7. *The body's books test flies a flowing body with a prescribed push* (two 0.0625 s macro steps from 893 K), not the spraying
   `melting_body`, which needs cantera and two 0.125 s steps. Why: the flow moves metal between columns, so the cells change every
   sub-step and the momentum books carry something; 13 s. Cost: spraying is not flown with the height map in the body tests, only in the
   resume test.
8. *No `sideways_heat_of` helper*: `args.sideways_heat` is read where it is used. Why: Task 5d's review called its `part_count_rule_of` a
   trivial wrapper. Cost: none.
9. *The resume refusal is tested on a saved flight of the element model* (a 0.075 s flight, about a second), the cheapest saved flight
   without skins.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -2263,3 +2263,158 @@
     flow.release = skin.prescribed_release(rates[0])                          # a release of one rate: no histogram
     one = skin.flow_substeps(mixed_sphere_skins(coarse_sphere_mesh)[0], flow, 0.01)
     assert one.hist_n is None and one.phase_mass.sum() > 0.0
+
+
+# -- sideways heat between cells that touch (sub-plan 22 Task 5f; Asha 2026-11-08) -------------------------------------
+def side_by_side(mat, dz, T, offset=0.0, w=0.5, normals=None, edge=(0, 1), area=4e-6, settings=None):
+    """Two hand-built columns 1 mm apart with one edge between them, of conduction weight w (no triangulation: the maps read
+    only the edges, their weights and incidence, the base corners and the columns' directions): cells of thicknesses dz[c]
+    [m], top first, at temperatures T[c] [K], both columns along +z unless `normals` are given, the second's base `offset`
+    [m] above the first's."""
+    import types
+    import scipy.sparse as sp
+    count = np.array([len(d) for d in dz])
+    K = int(count.max())
+    sk = skin.Skins(mat, np.full(2, area), settings or skin.SkinSettings(cell=20e-6))
+    m, temps = np.zeros((2, K)), np.full((2, K), 300.0)
+    for c in range(2):
+        m[c, :count[c]] = np.asarray(dz[c]) * mat.rho * area
+        temps[c, :count[c]] = T[c]
+    sk._install(m, m * mat.enthalpy(temps), count)
+    i, j = edge
+    normal = np.array([[0.0, 0.0, 1.0]] * 2) if normals is None else np.asarray(normals, dtype=float)
+    col = types.SimpleNamespace(n=2, edges=np.array([edge]), conduct_weight=np.array([w]), normal=normal,
+                                points=np.array([[0.0, 0.0, 0.0], [1e-3, 0.0, offset]]),
+                                incidence=sp.csr_matrix(([1.0, -1.0], ([i, j], [0, 0])), shape=(2, 1)))
+    return sk, col
+
+
+def stepped_sphere_skins(coarse_sphere_mesh, seed=4):
+    """`random_sphere_skins` with every base moved in or out along its direction by up to 0.3 mm at random, so that both the
+    tops and the bases step between neighbours."""
+    sk, col, rng = random_sphere_skins(coarse_sphere_mesh, seed)
+    col.update(col.points + rng.uniform(-0.3e-3, 0.3e-3, col.n)[:, None] * col.normal)
+    return sk, col, rng
+
+
+def test_height_pairing_equals_depth_pairing_where_tops_and_bases_are_level_and_cells_lie_on_the_depth_grid():
+    """Where neighbouring columns are equally deep on one base and every cell boundary is a breakpoint of the depth grid,
+    a height above the base is a depth below the top and each bin lies inside one cell of each column, so the two pairings
+    are one map: two 0.2 mm columns of 20 um bins (the dev cell), one in cells of 2, 3 and 5 bins and the other of 1, 4, 1
+    and 4, every cell at its own temperature in the Scheil alloy, whose conductivity varies with it; the heat at a random
+    field agrees to round-off."""
+    mat, settings = scheil(), skin.SkinSettings(cell=20e-6)
+    np.testing.assert_allclose(skin.depth_grid(settings, 0.2e-3)[:11], np.arange(11) * 20e-6, rtol=1e-12, atol=1e-18)
+    sk, col = side_by_side(mat, (20e-6 * np.array([2, 3, 5]), 20e-6 * np.array([1, 4, 1, 4])),
+                           ([520.0, 610.0, 700.0], [560.0, 650.0, 590.0, 730.0]), settings=settings)
+    x = np.where(sk.valid, np.random.default_rng(8).uniform(500.0, 800.0, sk.T.shape), 0.0)
+    height, depth = skin.SidewaysByHeight(sk, col).heat(x), skin.Sideways(sk, col).heat(x)
+    assert np.abs(depth).max() > 0.0
+    np.testing.assert_allclose(height, depth, rtol=1e-11, atol=1e-12 * np.abs(depth).max())
+
+
+def test_a_taller_columns_cells_above_its_neighbours_top_exchange_no_heat():
+    """Where the tops differ the part of the taller column above its neighbour's top is a side wall facing the gas: it
+    exchanges no heat with the neighbour (and the model does not heat or cool it from the gas). A column 1.0 mm deep in one
+    cell beside one 1.4 mm deep in cells of 0.2, 0.3 and 0.9 mm (top first), on one base: the top cell, 1.2-1.4 mm above
+    the base, gets no heat at any temperature and gives none; the next, 0.9-1.2 mm, exchanges over its 0.1 mm below the
+    neighbour's top alone, k w l (x_b - x_a), and the bottom one over its 0.9 mm."""
+    mat = sharp()                                                     # k = 100 W/(m K) at every temperature
+    sk, col = side_by_side(mat, ([1.0e-3], [0.2e-3, 0.3e-3, 0.9e-3]), (700.0, 700.0))
+    S = skin.SidewaysByHeight(sk, col)
+    x = np.array([[700.0, 0.0, 0.0], [1500.0, 800.0, 750.0]])
+    Q = S.heat(x)
+    cold = x.copy()
+    cold[1, 0] = 300.0
+    np.testing.assert_array_equal(S.heat(cold), Q)
+    assert Q[1, 0] == 0.0
+    np.testing.assert_allclose(Q[1, 1:], 100.0 * 0.5 * np.array([0.1e-3, 0.9e-3]) * (700.0 - x[1, 1:]), rtol=1e-12)
+    assert Q[0, 0] == pytest.approx(-Q[1, 1:].sum(), rel=1e-12) and np.all(Q[0, 1:] == 0.0)
+
+
+def test_height_pairing_pairs_cells_at_the_same_height_above_a_common_base():
+    """By hand (sub-plan 22 Task 5f): a column 1.0 mm deep in cells of 0.2, 0.3 and 0.5 mm (top first) beside one 1.4 mm
+    deep in cells of 0.4 and 1.0 mm, on one base. Above the base the first's cells span 0.8-1.0, 0.5-0.8 and 0-0.5 mm and
+    the second's 1.0-1.4 and 0-1.0 mm, so each of the first's touches the second's bottom cell alone, over 0.2, 0.3 and
+    0.5 mm, and the second's top cell touches none: heat k w l (x_b - x_a) a pair. Pairing by depth below the tops instead
+    puts the second's top cell beside the first's top two."""
+    mat = sharp()
+    sk, col = side_by_side(mat, ([0.2e-3, 0.3e-3, 0.5e-3], [0.4e-3, 1.0e-3]), (700.0, 700.0))
+    x = np.array([[600.0, 650.0, 700.0], [900.0, 800.0, 0.0]])
+    Q = skin.SidewaysByHeight(sk, col).heat(x)
+    first = 100.0 * 0.5 * np.array([0.2e-3, 0.3e-3, 0.5e-3]) * (800.0 - x[0])
+    np.testing.assert_allclose(Q[0], first, rtol=1e-12)
+    assert Q[1, 0] == 0.0 and Q[1, 2] == 0.0 and Q[1, 1] == pytest.approx(-first.sum(), rel=1e-12)
+    assert abs(skin.Sideways(sk, col).heat(x)[1, 0]) > 0.1 * np.abs(Q).max()
+
+
+def test_an_offset_between_bases_shifts_the_pairing_by_the_offset():
+    """Two columns of two 0.5 mm cells, the second's base 0.3 mm above the first's along their direction, s_ij = (p_j - p_i)
+    . n_ij with n_ij the unit mean of the two directions: in the first's frame the second's cells span 0.8-1.3 and 0.3-0.8
+    mm, so the first's top cell (0.5-1.0) touches both, over 0.2 and 0.3 mm, and its bottom cell (0-0.5) the second's
+    bottom cell over 0.2 mm. The same with the edge the other way round (s_ji = -s_ij) and with the two directions tilted 20
+    degrees either way about their mean, which leaves the offset along it unchanged; with the bases level, cell by cell."""
+    mat = sharp()
+    cells, x = ([0.5e-3, 0.5e-3], [0.5e-3, 0.5e-3]), np.array([[600.0, 650.0], [900.0, 800.0]])
+    c = 100.0 * 0.5 * np.array([0.2e-3, 0.3e-3, 0.2e-3])               # k w l of the three pairs [W/K]
+    to_first = np.array([c[0] * (900.0 - 600.0) + c[1] * (800.0 - 600.0), c[2] * (800.0 - 650.0)])
+    expected = np.array([to_first, [-c[0] * (900.0 - 600.0), -c[1] * (800.0 - 600.0) - c[2] * (800.0 - 650.0)]])
+    t = np.radians(20.0)
+    tilted = [[np.sin(t), 0.0, np.cos(t)], [-np.sin(t), 0.0, np.cos(t)]]
+    for kw in ({}, {"edge": (1, 0)}, {"normals": tilted}):
+        sk, col = side_by_side(mat, cells, (700.0, 700.0), offset=0.3e-3, **kw)
+        np.testing.assert_allclose(skin.SidewaysByHeight(sk, col).heat(x), expected, rtol=1e-12)
+    sk, col = side_by_side(mat, cells, (700.0, 700.0))
+    np.testing.assert_allclose(skin.SidewaysByHeight(sk, col).heat(x),
+                               100.0 * 0.5 * 0.5e-3 * np.array([[300.0, 150.0], [-300.0, -150.0]]), rtol=1e-12)
+
+
+def test_the_height_map_is_symmetric_and_negative_semi_definite_and_its_heat_sums_to_zero(coarse_sphere_mesh):
+    """The contract of `Sideways` that the implicit sub-step's conjugate gradients rely on, for the height map: on the test
+    sphere's columns of one to three times their minimum depth on bases stepped by up to 0.3 mm, y . S x = x . S y,
+    x . S x < 0, no heat at a uniform temperature, the heat summing to zero and none in a cell beyond a column's count."""
+    sk, col, rng = stepped_sphere_skins(coarse_sphere_mesh)
+    S = skin.SidewaysByHeight(sk, col)
+    x, y = (np.where(sk.valid, rng.normal(size=sk.T.shape), 0.0) for _ in range(2))
+    yx, xy = float((y * S.heat(x)).sum()), float((x * S.heat(y)).sum())
+    assert abs(yx - xy) < 1e-10 * abs(yx)
+    assert float((x * S.heat(x)).sum()) < 0.0
+    assert np.all(S.heat(np.where(sk.valid, 700.0, 0.0)) == 0.0)
+    Q = S.heat(sk.T)
+    assert abs(Q.sum()) < 1e-12 * np.abs(Q).sum() and np.all(Q[~sk.valid] == 0.0)
+
+
+def test_the_coupled_solve_converges_with_the_height_map(coarse_sphere_mesh, monkeypatch):
+    """The implicit sub-step with the height map, on columns of different depths on stepped bases melting under 3 MW/m2 on
+    the front, every cell at its own temperature: each Newton step solved to KRYLOV_TOL gives the temperatures that steps
+    solved to 1e-10 give, to the Newton tolerance, and the sideways books, which carry heat, sum to zero."""
+    out = []
+    for tol in (skin.KRYLOV_TOL, 1e-10):
+        monkeypatch.setattr(skin, "KRYLOV_TOL", tol)
+        sk, col, _ = stepped_sphere_skins(coarse_sphere_mesh)
+        q = 3e6 * np.clip(col.normal[:, 0], 0.0, None)
+        books = skin._conduct(sk, 0.01, q, sk.area, 0.0, 0.15, np.full(sk.n, 850.0), skin.SidewaysByHeight(sk, col))
+        assert abs(books[3].sum()) < 1e-12 * np.abs(books[3]).sum() and np.abs(books[3]).max() > 0.0
+        out.append(sk.T[sk.valid].copy())
+    assert np.abs(out[0] - out[1]).max() < 1e-5
+
+
+def test_a_height_map_reuses_the_last_geometry_only_while_the_cells_and_the_bases_stay(coarse_sphere_mesh):
+    """The height map depends on the base corners as well as the cells' thicknesses, so its geometry is reused only while
+    both stay: a base moved in place along its direction, or a cell grown at the bottom, makes it anew, as a map made from
+    nothing makes it; a depth map is never reused for it."""
+    sk, col, rng = random_sphere_skins(coarse_sphere_mesh)
+    first = skin.SidewaysByHeight(sk, col)
+    assert skin.SidewaysByHeight(sk, col, first).geometry is first.geometry
+    col.points += 0.1e-3 * (rng.uniform(size=col.n) > 0.5)[:, None] * col.normal          # moved in place
+    col.update(col.points)
+    moved = skin.SidewaysByHeight(sk, col, first)
+    assert moved.geometry is not first.geometry
+    np.testing.assert_array_equal(moved.heat(sk.T), skin.SidewaysByHeight(sk, col).heat(sk.T))
+    sk.add_bottom(np.full(sk.n, 1e-7), np.full(sk.n, float(sk.material.enthalpy(700.0))), np.zeros(sk.n))
+    sk.regrid()
+    again = skin.SidewaysByHeight(sk, col, moved)
+    assert again.geometry is not moved.geometry
+    np.testing.assert_array_equal(again.heat(sk.T), skin.SidewaysByHeight(sk, col).heat(sk.T))
+    np.testing.assert_array_equal(skin.SidewaysByHeight(sk, col, skin.Sideways(sk, col)).heat(sk.T), again.heat(sk.T))
+    assert skin.SIDEWAYS_PAIRINGS == ("depth", "height")
```

  in `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -1841,3 +1841,75 @@
     zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
     b.advance(0.1, 0.05, zero)
     assert b.handed_mass == handed
+
+
+# -- sideways heat between cells that touch (sub-plan 22 Task 5f; Asha 2026-11-08) ------------------------------------------
+def maps_made(monkeypatch):
+    """The names of the sideways maps the body makes, as it makes them."""
+    made = []
+
+    def counted(cls):
+        def make_map(*args, **kwargs):
+            made.append(cls.__name__)
+            return cls(*args, **kwargs)
+        return make_map
+    for name in ("Sideways", "SidewaysByHeight"):
+        monkeypatch.setattr(skin, name, counted(getattr(skin, name)))
+    return made
+
+
+def test_height_pairing_keeps_the_books(coarse_sphere_mesh, monkeypatch):
+    """`sideways_pairing="height"` (sub-plan 22 Task 5f) on a melting body that flows: 3 MW/m2 on the front of the 4 mm
+    sphere from 893 K, just under T_rigid, for two 0.0625 s macro steps (the melting part length), the liquid kept and
+    flowing under a Newtonian push, one flow sub-step a heat sub-step, so the cells' thicknesses change every sub-step and
+    metal is drawn and handed down. Every sideways map is the height map; the mass closes to 1e-12, the energy to 1e-8 of
+    the heat absorbed and the momentum to round-off."""
+    made = maps_made(monkeypatch)
+    b = make(coarse_sphere_mesh, T0=893.0, flow=True, remove=False, flow_cycles=1, sideways_pairing="height")
+    M0 = b.mass()
+    for k in range(2):
+        b.advance(0.0625 * (k + 1), 0.0625, front_loads(b, 3e6), forcing=newtonian(b))
+    assert b.sideways_pairing == "height" and made and set(made) == {"SidewaysByHeight"}
+    assert b.drawn_mass > 0.0 and b.handed_mass > 0.0 and b.skin_stats()["flowing_depth_max_mm"] > 0.0
+    assert b.mass() == pytest.approx(M0, rel=1e-12)
+    assert abs(b.energy_balance_residual()) < 1e-8
+    assert b.momentum_residual() < 1e-12 and b.impulse_size > 0.0
+
+
+def test_the_depth_pairing_by_default_changes_nothing_to_the_last_bit(coarse_sphere_mesh, monkeypatch):
+    """`sideways_pairing="depth"` and the argument left out are one flight, every column of every row and the final mass
+    and energy equal to the last bit, every map the depth map (sub-plan 22 Task 5f): 3 MW/m2 on the front from 893 K for two
+    0.0625 s macro steps, the liquid removed as it forms, numpy's generator seeded as a run seeds it before each."""
+    made = maps_made(monkeypatch)
+    runs = []
+    for kw in ({}, {"sideways_pairing": "depth"}):
+        np.random.seed(12345)
+        b = make(coarse_sphere_mesh, T0=893.0, **kw)
+        rows = []
+        for k in range(2):
+            b.advance(0.0625 * (k + 1), 0.0625, front_loads(b, 3e6))
+            rows.append((b.step_parts, timed_out(b.skin_stats())))
+        runs.append((b, rows))
+    assert same_rows(runs[0][1], runs[1][1]) and runs[0][0].removed_mass > 0.0 and runs[0][0].sideways_pairing == "depth"
+    assert runs[0][0].mass() == runs[1][0].mass() and runs[0][0].energy() == runs[1][0].energy()
+    assert made and set(made) == {"Sideways"} and skinbody.SkinBody.sideways_pairing == "depth"
+
+
+def test_a_body_saved_without_the_sideways_pairing_flies_on(coarse_sphere_mesh, monkeypatch):
+    """The states saved before the pairing existed have no `sideways_pairing` on the instance: the class's default gives
+    them the depth pairing, they fly on as before, and a resume can give them the height pairing (sub-plan 22 Task 5f). An
+    unknown pairing is refused."""
+    made = maps_made(monkeypatch)
+    b = make(coarse_sphere_mesh, T0=700.0)
+    b.emissivity = 0.0
+    del vars(b)["sideways_pairing"]
+    assert "sideways_pairing" not in vars(b) and b.sideways_pairing == "depth"
+    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
+    b.advance(0.05, 0.05, zero)
+    assert made and set(made) == {"Sideways"}
+    del made[:]
+    b.sideways_pairing = "height"
+    b.advance(0.1, 0.05, zero)
+    assert made and set(made) == {"SidewaysByHeight"}
+    with pytest.raises(ValueError, match="sideways_pairing"):
+        make(coarse_sphere_mesh, sideways_pairing="width")
```

  in `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -820,3 +820,51 @@
         with pytest.raises(_Built) as got:
             cli.main(argv)
         assert got.value.args[0]["hand_down"] is on
+
+
+def test_the_sideways_heat_reads_depth_or_height_and_marks_the_run_s_name_only_for_height():
+    """--sideways-heat (sub-plan 22 Task 5f) is `depth` by default, so that every existing run keeps its name and its record;
+    `height` ends a skin run's settings in _side-height, melting or not, after the hand-down's and before the seed's."""
+    p = cli.build_parser()
+    base = "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-sesam"
+    name = lambda argv: cli.run_name_of(p.parse_args(argv), 0.1, 7500.0, 77500.133, "us76", cli.continuum_step(p.parse_args(argv)))
+    melt = THERMAL + ["--melt", "on"]
+    assert p.parse_args(THERMAL).sideways_heat == "depth" == p.parse_args(melt).sideways_heat
+    assert p.parse_args(THERMAL + ["--sideways-heat", "height"]).sideways_heat == "height"
+    assert name(THERMAL) == name(THERMAL + ["--sideways-heat", "depth"]) == base
+    assert name(THERMAL + ["--sideways-heat", "height"]) == base + "_side-height"
+    assert (name(melt + ["--hand-down", "off", "--sideways-heat", "height", "--seed", "7"])
+            == base + "_melt-girin_handdown-off_side-height_seed-7")
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True, sideways_pairing="height") == "_side-height"
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, flow=True) == ""
+
+
+@pytest.mark.parametrize("argv", [
+    THERMAL + ["--sideways-heat", "width"],                                            # no such pairing
+    FEM + ["--sideways-heat", "height"],                                               # the element model has no skins
+    FEM + ["--melt", "on", "--sideways-heat", "height"],
+    BASE + ["--sideways-heat", "height"],                                              # no thermal model at all
+])
+def test_a_bad_sideways_heat_or_the_height_pairing_without_skins_exits_2(argv, tmp_path, capsys):
+    with pytest.raises(SystemExit) as exc:
+        cli.main(argv + ["--t-max", "1", "--outdir", str(tmp_path)])
+    assert exc.value.code == 2
+    if "width" not in argv:
+        assert "--sideways-heat height needs --surface-model skin" in capsys.readouterr().err
+
+
+def test_the_command_line_gives_the_skins_their_sideways_pairing(tmp_path, monkeypatch):
+    """--sideways-heat reaches the skin body as `depth` or `height`, whether it melts or not: the pairing concerns the heat
+    of any skins (sub-plan 22 Task 5f). The body is not built: its arguments are."""
+    from reentry_model import skinbody
+
+    def built(*args, **kwargs):
+        raise _Built(kwargs)
+    monkeypatch.setattr(skinbody, "SkinBody", built)
+    dry = THERMAL + ["--fine-band", "0", "--t-max", "1", "--outdir", str(tmp_path)]
+    melt = dry + ["--melt", "on"]
+    for argv, pairing in ((dry, "depth"), (dry + ["--sideways-heat", "height"], "height"), (melt, "depth"),
+                          (melt + ["--sideways-heat", "height"], "height")):
+        with pytest.raises(_Built) as got:
+            cli.main(argv)
+        assert got.value.args[0]["sideways_pairing"] == pairing
```

  and in `tests/test_reentry_model_restart.py`:

```diff
--- a/tests/test_reentry_model_restart.py
+++ b/tests/test_reentry_model_restart.py
@@ -363,3 +363,42 @@
     assert cli.main(["resume", kept, "--outdir", out, "--name", "back", "--hand-down", "on"]) == 0
     doc = json.load(open(os.path.join(out, "back.json")))
     assert "hand_down" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--hand-down": "on"}
+
+
+def test_resuming_with_height_pairing_sets_it_on_the_body_and_records_it(flown, tmp_path, capsys):
+    """`--sideways-heat height` on a resume (sub-plan 22 Task 5f) sets the pairing on the saved flight's skin body, which pairs
+    the heat otherwise from the next step on and which the one line before the flight names; it is recorded in the JSON's
+    settings and in what was changed and carried into the state the resumed flight saves, so a resume of that state keeps
+    it unless told, and `depth` asked of it gives back the saved flight's record (no key). `run` records it in the skins'
+    settings only when given, and the resume of a flight without skins refuses it (exit 2)."""
+    out, flight, records, state_path = flown
+    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "h", "--sideways-heat", "height",
+                     "--save-state-at", "0.2"]) == 0
+    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
+    assert "sideways heat height" in line
+    doc = json.load(open(str(tmp_path / "h.json")))
+    assert doc["resumed_from"]["changed"] == {"--sideways-heat": "height", "--save-state-at": 0.2}
+    assert doc["settings"]["sideways_heat"] == "height"
+    rows = rows_of(str(tmp_path / "h.csv"))
+    assert rows[0] == flight[14] and rows[1:] != flight[15:]
+    saved = restart.load(str(tmp_path / "h.state.pkl"))
+    assert saved["meta"]["args"]["sideways_heat"] == "height" and saved["run"].body.sideways_pairing == "height"
+    assert cli.main(["resume", str(tmp_path / "h.state.pkl"), "--outdir", str(tmp_path), "--name", "kept"]) == 0
+    assert "sideways heat height" in capsys.readouterr().out
+    assert json.load(open(str(tmp_path / "kept.json")))["settings"]["sideways_heat"] == "height"
+    assert cli.main(["resume", str(tmp_path / "h.state.pkl"), "--outdir", str(tmp_path), "--name", "back",
+                     "--sideways-heat", "depth"]) == 0
+    doc = json.load(open(str(tmp_path / "back.json")))
+    assert "sideways_heat" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--sideways-heat": "depth"}
+    assert "sideways_heat" not in json.load(open(os.path.join(out, "f.json")))["settings"]["skin"]
+    assert cli.main(FLIGHT + ["--t-max", "0.05", "--outdir", str(tmp_path), "--name", "run", "--sideways-heat", "height"]) == 0
+    assert json.load(open(str(tmp_path / "run.json")))["settings"]["skin"]["sideways_heat"] == "height"
+    k = FLIGHT.index("--melt")
+    elements = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--surface-model", "elements", "--quiet"]
+    assert cli.main(elements + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "el", "--save-state-at", "0.05"]) == 0
+    capsys.readouterr()
+    with pytest.raises(SystemExit) as exc:
+        cli.main(["resume", str(tmp_path / "el.state.pkl"), "--outdir", str(tmp_path / "o"), "--name", "r",
+                  "--sideways-heat", "height"])
+    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
+    assert "--sideways-heat height needs a saved flight of skins" in capsys.readouterr().err
```

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'SidewaysByHeight'` for the seven map tests and the three body tests (whose spy reads it first); on the command line
  `AttributeError` on `sideways_heat`, `KeyError: 'sideways_pairing'` in the body's arguments, and `unrecognized arguments:
  --sideways-heat height`, which exits 2 and so fails the three refusals on their message and the resume test; the refusal of an
  unknown pairing, which exits 2 without the option too, passes (16 failed, 1 passed).

- [ ] **Step 3: The map** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -878,6 +878,83 @@
         Qg = self.B @ (self.c * (xb[self.j] - xb[self.i]))
         Qcum = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(Qg, axis=1)], axis=1)
         return np.where(self.valid, np.diff(_interp_apply(self.geometry["to_cells"], Qcum), axis=1), 0.0)
+
+
+# -- heat sideways between cells that touch (sub-plan 22 Task 5f) ------------------------------------------------------
+SIDEWAYS_PAIRINGS = ("depth", "height")   # what pairs neighbouring columns' cells in the sideways conduction: their depth
+                                          # below the tops (`Sideways`, spec 6.3) or their height above the bases
+                                          # (`SidewaysByHeight`; Asha, 2026-11-08)
+
+
+def _height_geometry(sk, columns, dz):
+    """What a sub-step's height map is built on, from the cells' thicknesses dz and the base corners: every pair of cells on
+    an edge's two columns whose ranges of height overlap, as flat cell indices (cell a of column i is i K + a) with the
+    overlap [m] and the edge. Each edge's boundaries of both columns, in its first column's frame, are sorted downward in
+    one pass; between two consecutive boundaries lies one cell of each column (or none, above a top or below a base),
+    counted by the boundaries of each above it. Vectorised over the edges: no loop over edges or cells."""
+    n, K = dz.shape
+    zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)          # cell boundaries, from the top
+    y = zb[:, -1:] - zb                                     # their heights above the base: H to 0, and 0 beyond the cells
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    mean = columns.normal[i] + columns.normal[j]
+    s = np.einsum("ij,ij->i", columns.points[j] - columns.points[i], mean / np.linalg.norm(mean, axis=1)[:, None])
+    Y = np.concatenate([y[i], y[j] + s[:, None]], axis=1)  # (E, 2 K + 2): i's boundaries, then j's shifted into i's frame
+    order = np.argsort(-Y, axis=1, kind="stable")
+    t = np.take_along_axis(Y, order, axis=1)
+    above_i = np.cumsum(order <= K, axis=1)[:, :-1]         # i's boundaries above each piece between consecutive ones
+    a, b = above_i - 1, np.arange(2 * K + 1)[None, :] - above_i      # the cells of i and of j it lies in, from the top
+    length = t[:, :-1] - t[:, 1:]
+    e, p = np.nonzero((length > 0.0) & (a >= 0) & (a < sk.count[i][:, None]) & (b >= 0) & (b < sk.count[j][:, None]))
+    return {"a": i[e] * K + a[e, p], "b": j[e] * K + b[e, p], "edge": e, "length": length[e, p]}
+
+
+class SidewaysByHeight:
+    """The sideways conduction of one sub-step with the cells paired by their height above the bases, so that heat passes
+    between cells that touch (`--sideways-heat height`; Asha's idea and decision of 2026-11-08, sub-plan 22 Task 5f).
+    `Sideways` pairs cells at the same depth below their tops (spec 6.3), right while neighbouring tops stand level; with
+    the re-frozen metal kept in the skins (`hand_down=False`, Task 5e) the top steps about 1.4 mm across a 1.5 mm triangle
+    at the ring's edge, where pairing by depth puts a cell beside the gas next to metal 1.4 mm lower.
+
+    Here every cell's range is its height above its column's base, [H - z_bottom, H - z_top] (z from the top, H the
+    column's thickness). On each edge (i, j) the neighbour's ranges are shifted into i's frame by the offset between the
+    bases, s_ij = (p_j - p_i) . n_ij, p the base corners (`columns.points`) and n_ij the unit mean of the two columns'
+    directions (`columns.normal`), so s_ji = -s_ij; cells a of i and b of j exchange c_ab (x_b - x_a), with
+    c_ab = 0.5 (k_a + k_b) w l_ab, k the cells' conductivities at the sub-step's start, w the edge's conduction weight and
+    l_ab the overlap of the two ranges (zero where they do not overlap). c_ab is symmetric in a and b, so `heat` keeps
+    `Sideways`' contract with `_coupled_solve` and `_conduct` by construction: linear in x, antisymmetric pair by pair (it
+    sums to zero to round-off), symmetric and negative semi-definite, zero beyond a column's cells. Where neighbours are
+    equally deep on one base and their cells lie on the depth grid, it is the depth map to round-off.
+
+    What it brings out, stated: the part of a taller column above its neighbour's top is a side wall facing the gas; it
+    exchanges no heat with the neighbour, and the model does not heat or cool it from the gas (the gas reaches the tops
+    alone). Tracking every cell's position in space would give the same answer at 2 mm spacing at more cost and is not
+    chosen (the side session's call). The flow's sideways transport (spec 6.5, matched by fraction of the flowing depth)
+    is unchanged.
+
+    The geometry is rebuilt whenever the cells' thicknesses or the bases change, as the depth map's is with the cells. Its
+    cost, measured 2026-11-08 on a 2 mm surface of the 100 mm sphere (9 519 columns, 28 551 edges, 59-75 cells of 10 um,
+    the columns one to three times their minimum deep on bases stepped by up to 0.3 mm): 2.7 million pairs, 85 MB; a build
+    0.15 s against the depth map's 0.09 (0.37 GB at its peak against 0.09), a `heat` 0.020 s against 0.027, and one 10 ms
+    sub-step with every cell at its own temperature 3.6-3.8 s against 3.2-3.3, its line search evaluating the residual 49
+    times against 32 in the same 11 Newton steps."""
+
+    def __init__(self, sk, columns, previous=None):
+        """`previous`: the last sub-step's map, whose geometry is reused while the cells keep their thicknesses and the base
+        corners stand where they stood (a copy is kept: a base moved in place is a move)."""
+        self.dz, self.points = sk.dz(), np.array(columns.points, dtype=float)
+        same = (type(previous) is type(self) and np.array_equal(previous.dz, self.dz)
+                and np.array_equal(previous.points, self.points))
+        self.geometry = previous.geometry if same else _height_geometry(sk, columns, self.dz)
+        g = self.geometry
+        k = sk.props.k(sk.T).ravel()
+        self.c = 0.5 * (k[g["a"]] + k[g["b"]]) * columns.conduct_weight[g["edge"]] * g["length"]
+
+    def heat(self, x):
+        """The heat [W per cell] the neighbouring columns conduct into each cell at cell temperatures x (linear in x)."""
+        g = self.geometry
+        flat = x.ravel()
+        d = self.c * (flat[g["b"]] - flat[g["a"]])           # into cell a, out of cell b
+        return (np.bincount(g["a"], d, x.size) - np.bincount(g["b"], d, x.size)).reshape(x.shape)
 
 
 # -- viscosity (sub-plan 21 Task 3) ------------------------------------------------------------------------------------
```

- [ ] **Step 4: The body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -562,11 +562,15 @@
     `part_count_from`: `moves` (the default) or `melting`, what the next macro step's part count follows
     (`PART_COUNT_RULES`; sub-plan 22 Task 5d).
     `hand_down` [bool]: whether re-frozen metal at the bottom of a skin passes to the 3D model (the default), or stays in
-    the skins (sub-plan 22 Task 5e); False is refused with a `spread_limit`, which moves only metal being handed down."""
+    the skins (sub-plan 22 Task 5e); False is refused with a `spread_limit`, which moves only metal being handed down.
+    `sideways_pairing`: `depth` (the default) or `height`, what pairs neighbouring columns' cells in the sideways conduction,
+    their depth below the tops (`skin.Sideways`) or their height above the bases (`skin.SidewaysByHeight`, heat between
+    cells that touch; `skin.SIDEWAYS_PAIRINGS`, sub-plan 22 Task 5f)."""
     # class defaults as well as the instance's: the flights saved on 2026-11-06 were pickled before spreading existed
     spread_limit = None                  # deg
     part_count_from = "moves"            # the part count's rule (sub-plan 22 Task 5d): the states saved before it have none
     hand_down = True                     # re-frozen metal passes to the 3D model (Task 5e): the states saved before it have it on
+    sideways_pairing = "depth"           # the sideways conduction's pairing (Task 5f): the states saved before it have none
     spread_mass = spread_unplaced = 0.0  # kg: moved between columns since time zero; still wanting to move at the last move
     spread_over = 0                      # edges over the limit on the last move's first estimate
     spread_jump = 0.0                    # K: the largest landing of metal over its node's temperature this macro step
@@ -575,7 +579,11 @@
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
                  surface_flow=None, release=None, slurry_spray=False, spray=None, particles=True, spool=None,
-                 spread_limit=None, part_count_from="moves", hand_down=True):
+                 spread_limit=None, part_count_from="moves", hand_down=True, sideways_pairing="depth"):
+        if sideways_pairing not in skin_mod.SIDEWAYS_PAIRINGS:
+            raise ValueError("sideways_pairing must be one of {}, got {!r}".format(skin_mod.SIDEWAYS_PAIRINGS,
+                                                                                  sideways_pairing))
+        self.sideways_pairing = sideways_pairing
         if not hand_down and spread_limit is not None:
             raise ValueError("spreading moves only metal being handed down (Asha, 2026-11-08): hand_down=False needs "
                              "spread_limit=None, got {!r}".format(spread_limit))
@@ -802,7 +810,11 @@
                                  slurry=self.slurry_spray)
         eps = np.full(sk.n, self.emissivity)
         m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
-        side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
+        # implicit (sub-plan 21); with the height pairing heat passes between cells that touch (Asha's idea and decision of
+        # 2026-11-08, sub-plan 22 Task 5f: the deposit kept in the skins makes the depth pairing wrong at the ring's edge,
+        # where the top steps about 1.4 mm across a 1.5 mm triangle)
+        sideways_map = skin_mod.SidewaysByHeight if self.sideways_pairing == "height" else skin_mod.Sideways
+        side = lambda s, previous: sideways_map(s, self.columns, previous)
         nodes = mm.corners
         T_b0 = self.solver.temperature()[nodes]
         g0 = self._normal_gradient()                         # the profile's slope below every base, for this part
```

- [ ] **Step 5: The command line** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -5,12 +5,12 @@
         [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
         [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
         [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off --spread-limit DEG|off
-        --part-count-from moves|melting --hand-down on|off ...]
+        --part-count-from moves|melting --hand-down on|off --sideways-heat depth|height ...]
                                      (the layered skins, the default with --thermal fem, their flow and spraying the
                                       default with --melt on; the element model with --surface-model elements)
     python -m reentry_model resume <state> --outdir ... --name ... [--t-max S --dt-continuum S --seed N --frames-every N
         --slurry-spray on|off --spread-limit DEG|off --part-count-from moves|melting
-        --hand-down on|off]
+        --hand-down on|off --sideways-heat depth|height]
                                      (a flight saved with --save-state-at, continued; sub-plan 22 Task 5b)
     python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
 
@@ -202,7 +202,7 @@
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
                      fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=True, flow_cycles=None,
                      slurry_spray=False, spread_limit=None, part_count_from="moves",
-                     hand_down=True):
+                     hand_down=True, sideways_pairing="depth"):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
@@ -214,7 +214,8 @@
     sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3);
     `_spread-<deg>` where what the skins hand down is spread past a fold limit (sub-plan 22 Task 5c);
     `_parts-melting` where the part count follows the melting, not the moves (Task 5d); `_handdown-off` where re-frozen
-    metal stays in the skins (Task 5e); empty for the defaults."""
+    metal stays in the skins (Task 5e); `_side-height` where the sideways conduction pairs cells by height above the bases
+    (Task 5f); empty for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -246,6 +247,8 @@
         out += "_parts-{}".format(part_count_from)
     if not hand_down:
         out += "_handdown-off"
+    if sideways_pairing != "depth":
+        out += "_side-{}".format(sideways_pairing)
     return out
 
 
@@ -313,6 +316,10 @@
     rs.add_argument("--hand-down", choices=("on", "off"), default=None,
                     help="re-frozen metal passes to the 3D model, or stays in the skins; off needs the spread limit off, given "
                          "or carried by the saved flight (default: the saved flight's, on for one saved before the option existed)")
+    rs.add_argument("--sideways-heat", choices=skin.SIDEWAYS_PAIRINGS, default=None,
+                    help="what pairs neighbouring skin columns' cells in the sideways conduction: depth below the tops or "
+                         "height above the bases, heat between cells that touch (default: the saved flight's, depth for one "
+                         "saved before the option existed)")
     rs.add_argument("--save-state-at", type=float, default=None,
                     help="[s] of flight: save the resumed flight's whole state once, at the end of the first macro step ending at "
                          "or after it, as <outdir>/<name>{}; after the saved time and before --t-max (final review, 2026-11-07)".format(STATE_SUFFIX))
@@ -456,6 +463,11 @@
                          "model (on) or stays in the skins, which can heat it again, let it flow on or spray it (off; needs "
                          "--spread-limit off, spreading moves only metal being handed down; Asha, 2026-11-08, sub-plan 22 "
                          "Task 5e; default on)")
+    me.add_argument("--sideways-heat", choices=skin.SIDEWAYS_PAIRINGS, default="depth",
+                    help="with --surface-model skin: what pairs neighbouring columns' cells in the sideways conduction: depth "
+                         "below the tops (spec 6.3), or height above the bases, so that heat passes between cells that touch "
+                         "and a taller column's part above its neighbour's top exchanges none (Asha, 2026-11-08, sub-plan 22 "
+                         "Task 5f; default depth)")
     me.add_argument("--skin-flow-cycles", type=int, default=None,
                     help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
                          "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
@@ -546,7 +558,8 @@
                                      spray=spray_model_of(args, mat) if flow else None,
                                      slurry_spray=args.slurry_spray == "on", particles=args.particles,
                                      spool=spool if flow else None, spread_limit=spread_limit_of(args),
-                                     part_count_from=part_count_rule_of(args), hand_down=hand_down_of(args))
+                                     part_count_from=part_count_rule_of(args), hand_down=hand_down_of(args),
+                                     sideways_pairing=args.sideways_heat)
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray_model_of(args, mat)
@@ -585,6 +598,8 @@
             info["skin"]["part_count_from"] = part_count_rule_of(args)
         if not hand_down_of(args):
             info["skin"]["hand_down"] = "off"
+        if args.sideways_heat != "depth":
+            info["skin"]["sideways_heat"] = args.sideways_heat
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
@@ -636,7 +651,7 @@
                                                                     skin_flow_of(args) == "on" if args.melt == "on" else None,
                                                                     args.skin_flow_cycles, args.slurry_spray == "on",
                                                                     spread_limit_of(args), part_count_rule_of(args),
-                                                                    hand_down_of(args)))
+                                                                    hand_down_of(args), args.sideways_heat))
 
 
 def cmd_run(args, parser):
@@ -722,6 +737,8 @@
         parser.error("--hand-down off needs --surface-model skin and --melt on: it keeps the skins' re-frozen metal in them")
     if not hand_down_of(args) and spread_limit_of(args) is not None:
         parser.error(NO_SPREAD_WITHOUT_HAND_DOWN)
+    if args.sideways_heat != "depth" and args.surface_model != "skin":
+        parser.error("--sideways-heat height needs --surface-model skin: it pairs the skins' cells")
     if args.surface_model == "skin":
         for flag, given in (("--runoff off", args.runoff == "off"), ("--deep-runoff off", args.deep_runoff == "off"),
                             ("--molten-cascade off", args.molten_cascade == "off"),
@@ -948,6 +965,10 @@
         if args.hand_down == "off" and not (skins and run.body.material.melts):
             parser.error("--hand-down off needs a saved flight of melting skins: it keeps their re-frozen metal in them")
         changed["--hand-down"] = args.hand_down
+    if args.sideways_heat is not None:
+        if args.sideways_heat != "depth" and not skins:
+            parser.error("--sideways-heat height needs a saved flight of skins: it pairs their cells")
+        changed["--sideways-heat"] = args.sideways_heat
     # the spread limit and the hand-down as they will stand, each as given here or else the saved flight's
     if skins:
         will_spread = spread_limit_of(args) if args.spread_limit is not None else run.body.spread_limit
@@ -992,6 +1013,8 @@
         run.body.part_count_from = args.part_count_from
     if "--hand-down" in changed and skins:
         run.body.hand_down = args.hand_down == "on"
+    if "--sideways-heat" in changed and skins:
+        run.body.sideways_pairing = args.sideways_heat
     if isinstance(run.body.source_rows, spray.SourceSpool):      # only the records after the resume are written
         run.body.source_rows = spray.SourceSpool(os.path.join(run_dir, SOURCE_SPOOL))
     if args.seed is not None:
@@ -1003,14 +1026,15 @@
         setattr(saved, key, value)
     for key, flag in (("t_max", "--t-max"), ("dt_continuum", "--dt-continuum"), ("seed", "--seed"),
                       ("frames_every", "--frames-every"), ("slurry_spray", "--slurry-spray"), ("spread_limit", "--spread-limit"),
-                      ("part_count_from", "--part-count-from"), ("hand_down", "--hand-down")):
+                      ("part_count_from", "--part-count-from"), ("hand_down", "--hand-down"),
+                      ("sideways_heat", "--sideways-heat")):
         if flag in changed:
             setattr(saved, key, changed[flag])
     if state_path is not None:
         settings.save_meta = {"name": args.name, "args": dict(vars(saved))}
     limit = run.body.spread_limit if skins else None
     print("resuming {} (saved at {:.6g} s, step {}) to {:.6g} s: dt_continuum {}, seed {}, frames every {}, slurry spray {}, "
-          "spread limit {}, part count from {}, hand-down {}, save state at {}; files {}".format(
+          "spread limit {}, part count from {}, hand-down {}, sideways heat {}, save state at {}; files {}".format(
               os.path.abspath(args.state), state["time_s"], state["step"], run.sim.settings.t_max,
               settings.dt_continuum if settings.dt_continuum is not None else "none",
               args.seed if args.seed is not None else "saved random state", settings.frames_every,
@@ -1018,6 +1042,7 @@
               "off" if limit is None else "{:g} deg".format(limit),
               run.body.part_count_from if skins else "moves",
               "off" if skins and not run.body.hand_down else "on",
+              run.body.sideways_pairing if skins else "depth",
               "{:g} s".format(args.save_state_at) if args.save_state_at is not None else "no",
               os.path.join(args.outdir, args.name)), flush=True)
     history = run.run(resume=state)
@@ -1035,6 +1060,8 @@
         thermal_info["part_count_from"] = "melting"
     if skins and not run.body.hand_down:
         thermal_info["hand_down"] = "off"
+    if skins and run.body.sideways_pairing != "depth":
+        thermal_info["sideways_heat"] = run.body.sideways_pairing
     resumed = {"state": os.path.abspath(args.state), "step": state["step"], "time_s": state["time_s"],
                "run_name": state["meta"]["name"], "changed": changed}
     os.makedirs(args.outdir, exist_ok=True)
```

- [ ] **Step 6: Run the tests to see them pass,** then the unit tier: 639 passed, 1 skipped (Task 5e's 622 and these 17; 15.3 minutes across four processes, 2026-11-08).

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/skin.py reentry_model/skinbody.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "Skins: sideways heat can pass between cells that touch (--sideways-heat depth|height on run and resume, depth by default; with height the cells are paired by their height above the bases, the neighbour's shifted by the offset between the bases); sub-plan 22 Task 5f"
```

---

## The outline (Tasks 5–6, and Tasks 1–4, 5a, 5b, 5c, 5d, 5e and 5f written out above: 1 as Tasks 1a, 1b and 1c, 3 with 3b, 4 with 4b)

1. **Spraying per column** -- written out above as Tasks 1a (the amendments of 2026-10-22), 1b (the cap per crest)
   and 1c (the spraying from the top cells, by Asha's decisions of 2026-10-29). The outline's item, as it stood:
   > **Spraying per column**: b the flowing depth, δ_m from Girin's Eq. (2), the rigid-substrate regime test reading the
   > skin's depths, the thin, thick and front-surface Rayleigh–Taylor modes; the release leaves the top cells, liquid
   > first; the droplets' momentum leaves with them. `--slurry-spray` carries over off. Tests: each mode's release from a
   > prescribed column; books with spraying.

2. **Droplet records** -- written out above as Task 2. The outline's item, as it stood:
   > **Droplet records**: combined per column, per mode, per macro step; `solid_fraction` and `h_kg` added to the source
   > table. Tests: combining changes no radius; histograms exact.
3. **Retirement** -- written out above as Task 3, with the skins the default surface model and the spraying skin flight
   from the command line. The outline's item, as it stood:
   > **Retirement**: element fill fractions, `PHI_DEATH`, deaths and hand-overs, the film and deep accounts, the feed
   > rule and gate, the deep runoff, deep surfacing, the molten cascade, freeze-back and the stage-1 return device,
   > removed from the skin path; `--surface-model elements` keeps the current model bit for bit. Since sub-plan 20 the
   > skin path is `skinbody.SkinBody`, which never had those mechanisms, and stage 1's skins are `skin1d.py` behind
   > `MeltSettings(surface_model="skin1d")`: the retirement is deleting `skin1d.py`, its tests and `MeltingBody`'s
   > `skin1d` branches. The default surface model becomes `skin` here, and `--surface-model skin` with `--melt on` gains
   > `--removal girin`.
4. **Outputs** -- written out above as Task 4, with Register's release angle (Asha, 2026-10-29), the records on disk and
   the FEniCSx tier's failure, and Task 4b (the resolved sheared depth). The outline's item, as it stood:
   > **Outputs**: the history columns, frame fields and results of spec §10.2; for a spraying skin flight (Task 3) the
   > history's per-step droplet columns, so that `coupled.write_particles` writes its summary, and the results' sprayed
   > mass, droplets and row count; the droplet-momentum ratio of a vanishing last macro step (fact 179).
5. **Flights**: the 100 mm flight to 120 s and the 50 mm flight; convergence under halving of the macro step, 10 against
   20 µm cells, and a skin surface at 2 against 1 mm spacing; the acceptance criterion (droplet count, thick share,
   medians by number and mass, sprayed mass and mass at 120 s within the scatter between seeds under halving and
   refinement); the comparison with the element model's flight; cost and memory. (For the cost: Girin's evaluation
   once a macro step takes 0.04–0.05 s on the 100 mm sphere's 9 417 columns, fact 171 -- nothing beside a melting
   macro step's hour of processor time. The flights also report the metal the fault tolerance holds in the skins and
   the deposit's fold, spec §12, and measure again the share of the flow at short wavelengths, spec §14 assumption 1.)
   Where they run (Asha, 2026-10-31): every flight expected to take more than about two hours runs on the Mac Studio
   (Apple M5 Max, 6 performance cores, 64 GB; fact 149), from a runbook for the agent there in the pattern of
   `docs/superpowers/plans/2026-10-16-skin-level3-flights.md` -- from a tracked commit, no model code written, the
   results pushed on a branch -- the flights side by side, one core each. A skin flight uses about one core and its 3D
   solve is about 2 % of its time (fact 149: skins 81 %, 3D 1.8 %; the 65 km smoke flight of Task 4b: skins 87 %, of
   which the metal's flow 70 % of the whole, 3D 2.4 %), so the cores are used by running flights together, not by a
   parallel thermal solver (Asha declined it, 2026-10-31); side by side the level-3 flights each ran about 1.5 times
   slower a step than alone on the laptop. Asha's three decisions of 2026-10-31 for this task: (1) a cost probe first --
   about a second of melting flight on the production mesh (2 mm surface, 10 um cells), its time and memory a macro
   step measured before the flights are sized (extrapolated from Task 3's coarse smoke flight, 2-7 s a 12.5 ms step on
   2 367 columns of 20 um cells, a 100 mm flight's 7 640 melting steps take half a day to four days on one core);
   (2) the base flights whole, the convergence repeats -- the macro step halved, 10 against 20 um cells, a 1 against
   2 mm surface -- over a window of about 10 s around peak heating, so the convergence of quantities known only at the
   end (the mass sprayed by 120 s) is not checked; (3) a flight that stops on a fold of the surface keeps and judges
   what it made up to then, and reports where and when it folded, with options. The probe (fact 189): 4 425 s for one
   second of the 65 km case on the production mesh alone on the laptop, a 12.5 ms step rising from 15 to 89 s with
   the depth of the creeping slurry, the flow 79 % of the time and the 3D solve 1.6 %, at most 5.4 GB; so the 100 mm
   flight takes at least about 8 days on one core, extrapolated. Asha (2026-10-31): the flights run as they are -- no
   change to the model for their cost -- on the Mac Studio, from the runbook
   `docs/superpowers/plans/2026-10-31-skin-spraying-flights.md`: the 100 mm and 50 mm skin flights and the element
   model's two beside them, through `flights/fly.py` (beside the planning copy's `code/`, not part of it), which runs the
   command line unchanged and writes each history row and a progress line as it goes, so a flight that stops early
   keeps its history and its spooled droplet records.
6. **Documentation**: the README, `docs/model_assumptions.md` and the facts, as the repository's rules require.

## Inputs to settle

Carried from the outline, with one Task 2 adds and one Task 3 adds; Tasks 1–3 settle none of the outline's, their
spraying evaluated and recorded once a macro step at whatever macro step the run takes (Task 1c); Task 3b settles the
flow's default and Task 4 the records' volume, whatever the macro step.

- Sub-plans 19–21's measurements: the fold limit's resolution, the rebuild transfer, the skins' and the sideways solve's
  costs, which together decide the default macro step and whether the flights fit the laptop or need the cluster.
- Settled: the droplet records' volume (Task 2, ruling 9; Task 4): a skin flight writes a row per spraying column and
  mode per macro step, about 4 600 at 93 s on the 100 mm sphere with slurry spraying off and 9 200 with it on, so 35-70
  million rows on the 100 mm flight from its first liquid to 120 s at the continuum step, 7.6-15.2 GB as arrays; since
  Task 4 a spraying skin run writes them to disk as they come and the writer reads them back (fact 186), nothing
  combined over more than a macro step.
- The macro step with the flow on: sub-plan 21 (Task 6c) pairs `--skin-flow on` with `--removal girin`, whose continuum
  step (0.0125 s once the body's Knudsen number is below 0.01) the flow runs keep for now (Asha, 2026-10-22); this
  sub-plan tests other options, such as the 0.5 s macro step with the spraying taken in every skin sub-step. What the
  code does (Task 3): a skin run with the flow takes Girin's removal and so its continuum step, the `--dt` step (0.5 s)
  before the switch, and Girin's plan, the droplets' radius and the records once a macro step. Task 3b (Asha,
  2026-10-29): on the skins the continuum step also latches at their first liquid, which on the 100 mm flight adds
  about 25 s of 12.5 ms steps before the Knudsen switch (1 950 macro steps) and gives the 50 mm flight, whose Knudsen
  number never falls below 0.01, 25 s of them (estimated, fact 182). Task 5 measures; Asha decides the step itself.
- Settled: the flow's default (Task 3; Asha, 2026-10-29, option B, built in Task 3b: on for a melting skin run).
  `--skin-flow` was off (spec 10.1; `skinbody.FLOW_DEFAULT`, sub-plan 21 Task 6c's ruling that the flow is off until
  stage 4's spraying gives the liquid a sink, which it now does). Off, a melting run with no
  flags runs the heat-only skins, its liquid removed as it forms, and writes no droplets; on, it would run the flow and
  Girin's spraying from the top cells at the continuum step and write the droplet files. On the 100 mm sphere a 0.5 s
  melting macro step costs about an hour of processor time with the flow (fact 157: eight 0.0625 s parts, each run by
  both passes 1.8 times on average, at about 239 s) and, by the same count from Task 6d's timings without it (about
  30 s a part's two passes), about seven minutes without.
