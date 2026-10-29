# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1a, 1b, 1c and 2 (2026-10-28/29).**
> Task 1 is split in three, in the order the amendment of 2026-10-28 sets: 1a brings the element model's two spraying
> amendments of 2026-10-22 across, 1b adds the droplet's cap per crest to the shared spray physics, and 1c sprays from
> each column's top cells, by Asha's decisions of 2026-10-29 on how the skins' depths are read where liquid lies on slurry
> and on the slurry cascade. Task 2 makes the droplet records from the skins. Tasks 3–6 are the outline's, below, and the
> inputs it listed are carried at the end.

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
coupled writer serves both surface models.

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

## The outline (Tasks 3–6, and Tasks 1 and 2 written out above: 1 as Tasks 1a, 1b and 1c)

1. **Spraying per column** -- written out above as Tasks 1a (the amendments of 2026-10-22), 1b (the cap per crest)
   and 1c (the spraying from the top cells, by Asha's decisions of 2026-10-29). The outline's item, as it stood:
   > **Spraying per column**: b the flowing depth, δ_m from Girin's Eq. (2), the rigid-substrate regime test reading the
   > skin's depths, the thin, thick and front-surface Rayleigh–Taylor modes; the release leaves the top cells, liquid
   > first; the droplets' momentum leaves with them. `--slurry-spray` carries over off. Tests: each mode's release from a
   > prescribed column; books with spraying.

2. **Droplet records** -- written out above as Task 2. The outline's item, as it stood:
   > **Droplet records**: combined per column, per mode, per macro step; `solid_fraction` and `h_kg` added to the source
   > table. Tests: combining changes no radius; histograms exact.
3. **Retirement**: element fill fractions, `PHI_DEATH`, deaths and hand-overs, the film and deep accounts, the feed rule
   and gate, the deep runoff, deep surfacing, the molten cascade, freeze-back and the stage-1 return device, removed from
   the skin path; `--surface-model elements` keeps the current model bit for bit. Since sub-plan 20 the skin path is
   `skinbody.SkinBody`, which never had those mechanisms, and stage 1's skins are `skin1d.py` behind
   `MeltSettings(surface_model="skin1d")`: the retirement is deleting `skin1d.py`, its tests and `MeltingBody`'s
   `skin1d` branches. The default surface model becomes `skin` here, and `--surface-model skin` with `--melt on` gains
   `--removal girin`.
4. **Outputs**: the history columns, frame fields and results of spec §10.2.
5. **Flights**: the 100 mm flight to 120 s and the 50 mm flight; convergence under halving of the macro step, 10 against
   20 µm cells, and a skin surface at 2 against 1 mm spacing; the acceptance criterion (droplet count, thick share,
   medians by number and mass, sprayed mass and mass at 120 s within the scatter between seeds under halving and
   refinement); the comparison with the element model's flight; cost and memory. (For the cost: Girin's evaluation
   once a macro step takes 0.04–0.05 s on the 100 mm sphere's 9 417 columns, fact 171 -- nothing beside a melting
   macro step's hour of processor time. The flights also report the metal the fault tolerance holds in the skins and
   the deposit's fold, spec §12, and measure again the share of the flow at short wavelengths, spec §14 assumption 1.)
6. **Documentation**: the README, `docs/model_assumptions.md` and the facts, as the repository's rules require.

## Inputs to settle

Carried from the outline, with one Task 2 adds; none is settled by Tasks 1 and 2, whose spraying is evaluated and
recorded once a macro step at whatever macro step the run takes (Task 1c). The records' volume follows the macro step.

- Sub-plans 19–21's measurements: the fold limit's resolution, the rebuild transfer, the skins' and the sideways solve's
  costs, which together decide the default macro step and whether the flights fit the laptop or need the cluster.
- The droplet records' volume (Task 2, ruling 9): a skin flight writes a row per spraying column and mode per macro
  step, about 4 600 at 93 s on the 100 mm sphere with slurry spraying off and 9 200 with it on, so 370 000 to 740 000
  a second of flight at the 12.5 ms continuum step, 75-150 MB a second in memory; at 0.5 s macro steps, 40 times fewer.
  If the flights keep the continuum step, Task 5 writes the blocks out as it goes or combines the records over a longer
  window than the macro step (which spec 8 would then say).
- The macro step with the flow on: sub-plan 21 (Task 6c) pairs `--skin-flow on` with `--removal girin`, whose continuum
  step (0.0125 s once the body's Knudsen number is below 0.01) the flow runs keep for now (Asha, 2026-10-22); this
  sub-plan tests other options, such as the 0.5 s macro step with the spraying taken in every skin sub-step.
