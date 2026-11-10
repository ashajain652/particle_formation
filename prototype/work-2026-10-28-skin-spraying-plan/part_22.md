# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1a, 1b, 1c, 2, 3 and 3b (2026-10-28/29),
> 4 and 4b (2026-10-31), 5a (2026-11-01), 5b and 5c (2026-11-06), and 5d, 5e and 5f (2026-11-08), and 5g (2026-11-10).**
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

{{DIFF_TEST_SPRAY_T1A}}

- [ ] **Step 2: The melting body's tests** in `tests/test_reentry_model_melting.py`:

{{DIFF_TEST_MELTING_T1A}}

- [ ] **Step 3: The CLI's default** in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T1A}}

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

{{DIFF_DISPERSION_T1A}}

- [ ] **Step 6: The gates and the viscous film** in `reentry_model/spray.py`:

{{DIFF_SPRAY_T1A}}

- [ ] **Step 7: The reject recorded** in `reentry_model/body.py`:

{{DIFF_BODY_T1A}}

- [ ] **Step 8: Its history columns and results** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T1A}}

- [ ] **Step 9: The flag's default and help** in `reentry_model/cli.py`:

{{DIFF_CLI_T1A}}

- [ ] **Step 10: The sensitivity variant** in `analysis/melt_sensitivity.py`:

{{DIFF_MELT_SENSITIVITY_T1A}}

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

{{DIFF_TEST_SPRAY_T1B}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py -q -k "crest or
  cap_depth"`. Expected: 4 failed, `TypeError: SprayModel.evaluate() got an unexpected keyword argument 'cap_depth'`.

- [ ] **Step 3: The cap per crest** in `reentry_model/spray.py`:

{{DIFF_SPRAY_T1B}}

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

{{DIFF_TEST_SPRAY_T1C}}

- [ ] **Step 2: The cascade's tests** — append to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T1C}}

- [ ] **Step 3: Girin on the skins' tests** — append to `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T1C}}

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spray.py
  tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py -q -k "liquid_pool or the_cascade_takes or
  releases_the_cascade or girin_gives or cap_reads or slurry_cascade or sprays_by_girin"`. Expected: 7 failed --
  - `TypeError: SprayModel.evaluate() got an unexpected keyword argument 'rt_layer'`;
  - `AttributeError: module 'reentry_model.skin' has no attribute 'cascade_mass'`, and `... 'Cascade'`;
  - `AttributeError: module 'reentry_model.skinbody' has no attribute 'girin_plan'` (3);
  - `TypeError: SkinBody.__init__() got an unexpected keyword argument 'spray'`.

- [ ] **Step 5: The Rayleigh–Taylor mode's pool** in `reentry_model/spray.py`:

{{DIFF_SPRAY_T1C}}

- [ ] **Step 6: The cascade in the flow sub-steps** in `reentry_model/skin.py`:

{{DIFF_SKIN_T1C}}

- [ ] **Step 7: Girin once a macro step, the slurry's mode and the record** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T1C}}

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

{{DIFF_TEST_SPRAY_T2}}

- [ ] **Step 2: The release by the cells it came from, and the histograms** -- append to
  `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T2}}

- [ ] **Step 3: The plan's records and the body's** -- append to `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T2}}

- [ ] **Step 4: The element model's two columns** in `tests/test_reentry_model_melting.py` and
  `tests/test_reentry_model_coupled.py`:

{{DIFF_TEST_MELTING_T2}}

{{DIFF_TEST_COUPLED_T2}}

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

{{DIFF_SPRAY_T2}}

- [ ] **Step 7: The element model's rows** in `reentry_model/body.py`:

{{DIFF_BODY_T2}}

- [ ] **Step 8: The release by phase, and the histograms in the flow sub-steps** in `reentry_model/skin.py`:

{{DIFF_SKIN_T2}}

- [ ] **Step 9: The plan's records, rows and histogram, and the body's** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T2}}

- [ ] **Step 10: The writer** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T2}}

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

{{DIFF_TEST_CLI_T3}}

- [ ] **Step 2: Nothing imports the retired skins** -- append to `tests/test_reentry_model_melting.py`; the
  verification script's command line names its surface model, in `tests/test_reentry_model_thermal_verification.py`:

{{DIFF_TEST_MELTING_T3}}

{{DIFF_TEST_THERMAL_VERIFICATION_T3}}

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

{{DIFF_TEST_FENICSX_T3}}

- [ ] **Step 5: The element model without them** in `reentry_model/body.py`:

{{DIFF_BODY_T3}}

- [ ] **Step 6: Its history, results and frames without them** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T3}}

- [ ] **Step 7: The skins the default, the run names, and the spraying skin flight** in `reentry_model/cli.py`:

{{DIFF_CLI_T3}}

- [ ] **Step 8: The element model's scripts name it** in `analysis/reentry_model_thermal_verification.py`,
  `analysis/melt_verification.py` and `analysis/melt_sensitivity.py`:

{{DIFF_THERMAL_VERIFICATION_T3}}

{{DIFF_MELT_VERIFICATION_T3}}

{{DIFF_MELT_SENSITIVITY_T3}}

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

{{DIFF_TEST_CLI_T3B}}

- [ ] **Step 2: The switch's two triggers, and the skins' liquid** -- append to `tests/test_reentry_model_coupled.py`
  and `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_COUPLED_T3B}}

{{DIFF_TEST_SKINBODY_T3B}}

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

{{DIFF_COUPLED_T3B}}

- [ ] **Step 5: The skins' liquid** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T3B}}

- [ ] **Step 6: The flow on by default, and the skins' switch** in `reentry_model/cli.py`:

{{DIFF_CLI_T3B}}

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

{{DIFF_TEST_SPRAY_T4}}

- [ ] **Step 2: The loop's last step, the writer and the rebuild's mass** -- append to
  `tests/test_reentry_model_trajectory.py`, `tests/test_reentry_model_coupled.py` and
  `tests/test_reentry_model_movemesh.py`:

{{DIFF_TEST_TRAJECTORY_T4}}

{{DIFF_TEST_COUPLED_T4}}

{{DIFF_TEST_MOVEMESH_T4}}

- [ ] **Step 3: The skin body's outputs** -- in `tests/test_reentry_model_skinbody.py` (and sub-plan 21 Task 6c's
  flow-off test comparing the two bodies' history columns, now nan-aware: the new columns are empty in both):

{{DIFF_TEST_SKINBODY_T4}}

- [ ] **Step 4: The element model's and the command line's** -- in `tests/test_reentry_model_melting.py`,
  `tests/test_reentry_model_cli.py` and `tests/test_reentry_model_fenicsx.py`:

{{DIFF_TEST_MELTING_T4}}

{{DIFF_TEST_CLI_T4}}

{{DIFF_TEST_FENICSX_T4}}

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

{{DIFF_SPRAY_T4}}

- [ ] **Step 7: The element model's rows and history** in `reentry_model/body.py`:

{{DIFF_BODY_T4}}

- [ ] **Step 8: The last step** in `reentry_model/trajectory.py`, and the summary, the skin results and the loop in
  `reentry_model/coupled.py`:

{{DIFF_TRAJECTORY_T4}}

{{DIFF_COUPLED_T4}}

- [ ] **Step 9: The rebuild's mass residual** in `reentry_model/movemesh.py`:

{{DIFF_MOVEMESH_T4}}

- [ ] **Step 10: The skin body's outputs** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T4}}

- [ ] **Step 11: The spool and the skins' last step from the command line** in `reentry_model/cli.py`:

{{DIFF_CLI_T4}}

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

{{DIFF_TEST_SKINBODY_T4B}}

- [ ] **Step 2: Run it to see it fail.** `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -m "not drama and not
  reference" -k straight_line_depth`. Expected: 1 failed -- `AssertionError: assert not True` (`hasattr(skinbody,
  'SHEARED_FRACTION')`).

- [ ] **Step 3: The straight-line equivalent** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T4B}}

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

{{FILE_TEST_RESTART_T5A}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_restart.py -q`. Expected: a
  collection error -- `ImportError: cannot import name 'restart' from 'reentry_model'`.

- [ ] **Step 3: The save and the load** -- `reentry_model/restart.py`:

{{FILE_RESTART_T5A}}

- [ ] **Step 4: The save on the way, and the resume**, in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T5A}}

- [ ] **Step 5: `--save-state-at`** in `reentry_model/cli.py`:

{{DIFF_CLI_T5A}}

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

{{DIFF_TEST_RESTART_T5B}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_restart.py -q -k resume_command`.
  Expected: failures with `argument command: invalid choice: 'resume' (choose from run, compare)`.

- [ ] **Step 3: The command** in `reentry_model/cli.py`:

{{DIFF_CLI_T5B}}

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

{{FILE_TEST_SPREADING_T5C}}

  in `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T5C}}

  in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T5C}}

  and in `tests/test_reentry_model_restart.py`:

{{DIFF_TEST_RESTART_T5C}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_spreading.py -q`. Expected: `ImportError: cannot
  import name 'spreading'`. The others, with the module in place: failures on `spread_limit`, `parse_spread_limit` and `--spread-limit`,
  except the refusals, which an unknown option also exits 2 on, and which mean something once the option exists.

- [ ] **Step 3: The pass**, `reentry_model/spreading.py`:

{{FILE_SPREADING_T5C}}

- [ ] **Step 4: The body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T5C}}

- [ ] **Step 5: The results and the column list's note** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T5C}}

- [ ] **Step 6: The command line** in `reentry_model/cli.py`:

{{DIFF_CLI_T5C}}

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

{{DIFF_TEST_SKINBODY_T5D}}

  in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T5D}}

  and in `tests/test_reentry_model_restart.py`:

{{DIFF_TEST_RESTART_T5D}}

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError` and `TypeError` on `part_count_from` and
  `PART_COUNT_RULES`; the refusals of an unknown rule, which exit 2 without the option too, pass.

- [ ] **Step 3: The body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T5D}}

- [ ] **Step 4: The command line** in `reentry_model/cli.py`:

{{DIFF_CLI_T5D}}

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

{{DIFF_TEST_SKINBODY_T5E}}

  in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T5E}}

  and in `tests/test_reentry_model_restart.py`:

{{DIFF_TEST_RESTART_T5E}}

- [ ] **Step 2: Run them to see them fail.** Expected: `TypeError` on `hand_down` for the six body tests, `AttributeError` on
  `hand_down_of` and unrecognised-argument exits for the command line's.

- [ ] **Step 3: The body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T5E}}

- [ ] **Step 4: The command line** in `reentry_model/cli.py`:

{{DIFF_CLI_T5E}}

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

{{DIFF_TEST_SKIN_T5F}}

  in `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T5F}}

  in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T5F}}

  and in `tests/test_reentry_model_restart.py`:

{{DIFF_TEST_RESTART_T5F}}

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'SidewaysByHeight'` for the seven map tests and the three body tests (whose spy reads it first); on the command line
  `AttributeError` on `sideways_heat`, `KeyError: 'sideways_pairing'` in the body's arguments, and `unrecognized arguments:
  --sideways-heat height`, which exits 2 and so fails the three refusals on their message and the resume test; the refusal of an
  unknown pairing, which exits 2 without the option too, passes (16 failed, 1 passed).

- [ ] **Step 3: The map** in `reentry_model/skin.py`:

{{DIFF_SKIN_T5F}}

- [ ] **Step 4: The body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T5F}}

- [ ] **Step 5: The command line** in `reentry_model/cli.py`:

{{DIFF_CLI_T5F}}

- [ ] **Step 6: Run the tests to see them pass,** then the unit tier: 639 passed, 1 skipped (Task 5e's 622 and these 17; 15.3 minutes across four processes, 2026-11-08).

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/skin.py reentry_model/skinbody.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "Skins: sideways heat can pass between cells that touch (--sideways-heat depth|height on run and resume, depth by default; with height the cells are paired by their height above the bases, the neighbour's shifted by the offset between the bases); sub-plan 22 Task 5f"
```

---

### Task 5g: The windward loads from a sphere fitted to the skins' tops

**Files:**
- Modify: `reentry_model/skinbody.py` (`NOSE_COLUMNS`, `NoseFit`, the class defaults and `SkinBody(..., size_feedback="current")`,
  `_fit_nose` called from `_refresh_geometry`, `nose_fit`, `nose_radius()`, `nose_stats()`, the five columns in `skin_stats()`)
- Modify: `reentry_model/coupled.py` (`SKINBODY_COLUMNS` gains the five columns)
- Modify: `reentry_model/cli.py` (`--size-feedback current|initial` accepted with `--surface-model skin` on `run`, default `current`,
  and on `resume`; the records; the one refusal of it with the skins goes)
- Test: `tests/test_reentry_model_skinbody.py` (9 tests), `tests/test_reentry_model_cli.py` (2 tests, one refusal case changed),
  `tests/test_reentry_model_restart.py` (1 test)

**Interfaces:**
- Consumes: `body.fit_sphere(points) -> (centre, radius)`, `body.NOSE_CAP_ANGLE` (30.0), `body.NOSE_CAP_FACTOR` (1.67),
  `body.SIZE_FEEDBACK_NAMES`; `SkinBody.surface` (a `columns.CornerSurface` of the column tops), `SkinBody.theta`, `SkinBody.v_hat`;
  Task 5e's and 5f's wiring of an option on `run` and `resume`. The element model is not changed.
- Produces: `SkinBody(..., size_feedback="current")`, also a class default `"initial"` for states pickled before it (`ValueError` for
  another value); `SkinBody.nose_fit -> NoseFit(fitted, bound, points, concave, used)`; `SkinBody.nose_radius()` returns `nose_fit.used`;
  history columns `nose_radius_mm`, `fitted_nose_radius_mm` (not-a-number where there is no fit), `nose_bound_mm`, `nose_cap_points`,
  `nose_cap_concave` with either setting. `run`: `--size-feedback initial|current` with `--surface-model skin`, default `current`, the
  JSON's `settings.skin.size_feedback` always. `resume`: `--size-feedback`, default the saved body's (`initial` for a state without it),
  exit 2 for a saved flight without skins; JSON `settings.size_feedback = "current"` only when current, `changed["--size-feedback"]` when
  given; the one line before the flight names it; a state saved before the columns has them added to its saved row.

**What was decided** (design `docs/superpowers/specs/2026-11-10-fitted-nose-radius-design.md`; Asha, 2026-11-10: fit the gas-molten
interface, the bound a comparison and not a limit, on by default). `SkinBody.nose_radius()` returned the mesh's original radius whatever
the shape, and it feeds the stagnation heating and the surface flow (the pressure gradient, the boundary layer, the shock Reynolds number).
Now the cap is the top corners whose normals are within 30 degrees of the direction of motion, fitted unweighted by the element model's
`fit_sphere`; fewer than four, or a radius not finite and positive, is no fit and the mesh's radius is used. A cap whose fitted centre lies
on the flow side of its mean point is concave and uses the bound. The bound is 1.67 times the windward transverse radius (the largest
distance from the flight axis through the mesh's mean point of the top corners at most 90 degrees from the flight direction, so a pool
behind the body does not enter it), reported and never applied: a flat cap reads above it. The fit is made in `_refresh_geometry` from the
top surface as built, never modifying the surface or `theta`, so a macro step's loads read the shape at its start as the angles do, and
`initial` reproduces every earlier column to the last bit. The droplet cap's radius (`transverse_radius()`) and the body Knudsen number's
length are not changed. **Measured** on the comparison flights' frames, before the build (2026-11-10): the cap's radius is 49.7 mm at
22.5 s and about 95-107 mm from about 30 s on the 100 mm sphere, and 24.7 mm at 171.0 s, 83 mm at 176.0 s, 45 mm at 181.0 s on the 50 mm
sphere; at twice the radius the stagnation heating is about 30 % lower and the pressure gradient half.

**Rulings** (2026-11-10; what, why, cost if wrong):
1. *Concave means `(centre - cap mean) . v_hat > 0`*, the controller's reading. Why: the stagnation corner's normal points along `v_hat`, so
   a convex cap's centre lies behind it. Cost: a nearly flat cap with noise may read either way; concave uses the bound, which is the
   flat-face value, so the two readings are close in what they give.
2. *A flat cap's radius is not bounded* (Asha's decision, item 2). A fitted radius of about 1 m on a flat dished cap would be read as it
   is under `current`; the column `nose_bound_mm` beside `nose_radius_mm` shows it. Cost: heating falls with the square root of the
   radius, so a nearly flat cap is nearly unheated at the stagnation point.
3. *`used` is read from the setting at each ask, the fit made once a refresh*: a resume can change the setting without a refit, and a state
   saved before the task makes its fit at the first ask (`_nose_raw` class default None).
4. *The run's name carries no suffix for `initial`.* Why: the brief names none; the JSON records the setting always. Cost: two flights
   differing only in it share a name unless `--name` is given.
5. *The tests shape the cap by moving the tops* on the coarse sphere and refitting (`with_tops`), not by flying to it; the untouched
   sphere and the resumed flights are the end-to-end cases.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T5G}}

  in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T5G}}

  and in `tests/test_reentry_model_restart.py`:

{{DIFF_TEST_RESTART_T5G}}

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError` on `nose_fit` and `skinbody.NOSE_COLUMNS`, `TypeError` for the
  argument `size_feedback`, and the command line's `--size-feedback` refused with the skins (exit 2).

- [ ] **Step 3: The body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T5G}}

- [ ] **Step 4: The history's columns** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T5G}}

- [ ] **Step 5: The command line** in `reentry_model/cli.py`:

{{DIFF_CLI_T5G}}

- [ ] **Step 6: Run the tests to see them pass,** then the unit tier: 651 passed, 1 skipped (Task 5f's 639 and these 12; 36.9 minutes across four processes with the two long runs sharing the laptop, 2026-11-10; the first run had one failure, the skin run's pinned settings, which now carry `size_feedback`, and that test was re-run alone and passed).

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/skinbody.py reentry_model/coupled.py reentry_model/cli.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "Skins: the windward loads read the radius of a sphere fitted to the tops of the columns facing the flow (--size-feedback current|initial on run and resume, current by default; five nose columns in the history; the bound reported, never applied); sub-plan 22 Task 5g"
```

---

## The outline (Tasks 5–6, and Tasks 1–4, 5a, 5b, 5c, 5d, 5e, 5f and 5g written out above: 1 as Tasks 1a, 1b and 1c, 3 with 3b, 4 with 4b)

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
