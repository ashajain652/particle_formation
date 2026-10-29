# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1a, 1b and 1c (2026-10-28/29).** Task
> 1 is split in three, in the order the amendment of 2026-10-28 sets: 1a brings the element model's two spraying
> amendments of 2026-10-22 across, 1b adds the droplet's cap per crest to the shared spray physics, and 1c sprays from
> each column's top cells, by Asha's decisions of 2026-10-29 on how the skins' depths are read where liquid lies on slurry
> and on the slurry cascade. Tasks 2–6 are the outline's, below, and the inputs it listed are carried at the end.

**Goal:** Girin's spraying acts on each column's top cells, with the film thickness b the whole flowing depth, by the
rules of record (the element model's amendments of 2026-10-22) and a droplet no larger than one crest of its own wave;
the old film, deep account and their rules are removed; and the model flies the 100 mm and 50 mm cases to the spec's
acceptance.

**Architecture:** Task 1a changes only the shared spray physics and the element model's callers, exactly as branch
`we-critical-308` changed them; Task 1b gives the shared spray physics a second size cap, used only when its caller
passes the depth the release can take, so the element model is unchanged. Task 1c evaluates Girin per column once at
the start of each macro step, from the surface flow the skins' flow already evaluates there and the skins' depths then,
and hands the rate per column to sub-plan 21 Task 7's release hook, which takes it from the top cells in every flow
sub-step.

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

## The outline (Tasks 2–6, and Task 1 written out above as Tasks 1a, 1b and 1c)

1. **Spraying per column** -- written out above as Tasks 1a (the amendments of 2026-10-22), 1b (the cap per crest)
   and 1c (the spraying from the top cells, by Asha's decisions of 2026-10-29). The outline's item, as it stood:
   > **Spraying per column**: b the flowing depth, δ_m from Girin's Eq. (2), the rigid-substrate regime test reading the
   > skin's depths, the thin, thick and front-surface Rayleigh–Taylor modes; the release leaves the top cells, liquid
   > first; the droplets' momentum leaves with them. `--slurry-spray` carries over off. Tests: each mode's release from a
   > prescribed column; books with spraying.

2. **Droplet records**: combined per column, per mode, per macro step; `solid_fraction` and `h_kg` added to the source
   table. Tests: combining changes no radius; histograms exact.
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

Carried from the outline; neither is settled by Task 1, whose spraying is evaluated once a macro step at whatever
macro step the run takes (Task 1c).

- Sub-plans 19–21's measurements: the fold limit's resolution, the rebuild transfer, the skins' and the sideways solve's
  costs, which together decide the default macro step and whether the flights fit the laptop or need the cluster.
- The macro step with the flow on: sub-plan 21 (Task 6c) pairs `--skin-flow on` with `--removal girin`, whose continuum
  step (0.0125 s once the body's Knudsen number is below 0.01) the flow runs keep for now (Asha, 2026-10-22); this
  sub-plan tests other options, such as the 0.5 s macro step with the spraying taken in every skin sub-step.
