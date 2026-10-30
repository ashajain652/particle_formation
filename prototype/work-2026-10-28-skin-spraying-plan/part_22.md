# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1a, 1b, 1c, 2, 3 and 3b (2026-10-28/29).**
> Task 1 is split in three, in the order the amendment of 2026-10-28 sets: 1a brings the element model's two spraying
> amendments of 2026-10-22 across, 1b adds the droplet's cap per crest to the shared spray physics, and 1c sprays from
> each column's top cells, by Asha's decisions of 2026-10-29 on how the skins' depths are read where liquid lies on slurry
> and on the slurry cascade. Task 2 makes the droplet records from the skins. Task 3 retires stage 1's one-dimensional
> skins, makes the layered skins the default surface model and puts the spraying skin flight on the command line; Task
> 3b, by Asha's decisions on the two questions it left, makes the flow the default of a melting skin run and latches
> the skins' continuum step at their first liquid. Tasks 4–6 are the outline's, below, and the inputs it listed are
> carried at the end.

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

## The outline (Tasks 4–6, and Tasks 1–3 written out above: 1 as Tasks 1a, 1b and 1c, 3 with 3b)

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
4. **Outputs**: the history columns, frame fields and results of spec §10.2; for a spraying skin flight (Task 3) the
   history's per-step droplet columns, so that `coupled.write_particles` writes its summary, and the results' sprayed
   mass, droplets and row count; the droplet-momentum ratio of a vanishing last macro step (fact 179).
5. **Flights**: the 100 mm flight to 120 s and the 50 mm flight; convergence under halving of the macro step, 10 against
   20 µm cells, and a skin surface at 2 against 1 mm spacing; the acceptance criterion (droplet count, thick share,
   medians by number and mass, sprayed mass and mass at 120 s within the scatter between seeds under halving and
   refinement); the comparison with the element model's flight; cost and memory. (For the cost: Girin's evaluation
   once a macro step takes 0.04–0.05 s on the 100 mm sphere's 9 417 columns, fact 171 -- nothing beside a melting
   macro step's hour of processor time. The flights also report the metal the fault tolerance holds in the skins and
   the deposit's fold, spec §12, and measure again the share of the flow at short wavelengths, spec §14 assumption 1.)
6. **Documentation**: the README, `docs/model_assumptions.md` and the facts, as the repository's rules require.

## Inputs to settle

Carried from the outline, with one Task 2 adds and one Task 3 adds; none is settled by Tasks 1–3, whose spraying is
evaluated and recorded once a macro step at whatever macro step the run takes (Task 1c). The records' volume follows the
macro step.

- Sub-plans 19–21's measurements: the fold limit's resolution, the rebuild transfer, the skins' and the sideways solve's
  costs, which together decide the default macro step and whether the flights fit the laptop or need the cluster.
- The droplet records' volume (Task 2, ruling 9): a skin flight writes a row per spraying column and mode per macro
  step, about 4 600 at 93 s on the 100 mm sphere with slurry spraying off and 9 200 with it on, so 370 000 to 740 000
  a second of flight at the 12.5 ms continuum step, 75-150 MB a second in memory; at 0.5 s macro steps, 40 times fewer.
  If the flights keep the continuum step, Task 5 writes the blocks out as it goes or combines the records over a longer
  window than the macro step (which spec 8 would then say).
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
