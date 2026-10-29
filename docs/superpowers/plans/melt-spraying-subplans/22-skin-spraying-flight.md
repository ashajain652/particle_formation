# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1a and 1b (2026-10-28).** Task 1 is
> split in three, in the order the amendment of 2026-10-28 sets: 1a brings the element model's two spraying amendments
> of 2026-10-22 across, 1b adds the droplet's cap per crest to the shared spray physics, and 1c sprays from each
> column's top cells. **Task 1c waits for Asha's decision** on how the skins' depths are read where liquid lies on
> slurry (Task 1c below, decision 1) and on the slurry cascade, which the viscous film makes inert (decision 2). Tasks
> 2–6 are the outline's, below, and the inputs it listed are carried at the end.

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
   damping larger) carries over (Task 1c).
6. **Girin once a macro step** (Task 1c's design) — the rate, the mode and the radius per column are fixed for the
   macro step from its start, as the cap is; a column that becomes wet or changes regime within it waits for the next,
   12.5 ms at the flow-on runs' continuum step and up to 0.5 s at the provisional default. Task 5's halving of the macro
   step measures it.

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

### Task 1c: Spraying from each column's top cells (outline; waits for Asha's decisions 1 and 2)

What it builds (spec §8, §6.4, §6.7; the decisions of 2026-10-11, 2026-10-24 and 2026-10-28), on Tasks 1a and 1b and
sub-plan 21 Task 7's release hook (`SkinBody(release=..., slurry_spray=...)`, `skin.Flow.release`, `skin.release_top`):

- **Girin once a macro step** (decided here, within spec §8). At the start of each macro step, in `SkinBody.advance`,
  from the surface flow the flow's forcing evaluates there (`flow_forcing`, the step's aero state and the tops'
  temperatures at its start) and the skins' depths then: δ_m by Girin's Eq. (2) (`spray.melt_layer`); his surface
  velocity by his closure, as the element model's (`film.lubrication` on b, the thick branch's τ δ_m / μ_l, which alone
  enters a gate); b the flowing depth `sk.depth_above(T_rigid)` and the regime test on it (spec §6.7); the viscous
  film's floor on it (2026-10-28); the windward columns θ ≤ 90°; the molten region on the columns' own graph; the
  quarter-body cap from the top surface's transverse radius; and the cap's depth, the liquid depth
  `sk.depth_above(T_liquidus)` with slurry spraying off and the flowing depth with it on (2026-10-24). That gives each
  column its mode, wavelength, radius and rate per unit top area for the macro step, held through its parts, passes and
  sub-steps. Why: spec §8 fixes a column's radius within a macro step by that step's gas state and δ_m, so the records
  combine per column, per mode, per macro step without changing a radius; the cap is evaluated at the start by Asha's
  decision of 2026-10-24; the gas state changes over seconds (§9); and the element model evaluates its spray once a macro
  step too. Measured on the 100 mm sphere's 9 417 columns: 0.04–0.05 s an evaluation, nearly all of it reading the two
  depths, against 0.06 s for the surface flow and about 2 s a flow sub-step (fact 171). Cost if wrong: a column that
  becomes wet or changes regime within a macro step waits for the next — 12.5 ms at the flow-on runs' continuum step, up
  to 0.5 s at the provisional default — which Task 5's halving of the macro step measures; the alternative, the
  depth-dependent parts again in every flow sub-step with the gas state held, costs about 0.04 s each.
- **Re-evaluated in every flow sub-step:** only what the admitted run holds. The release is the macro step's rate times
  the top area and the flow sub-step, from the top cells, liquid first, the exposed slurry too only with
  `--slurry-spray on`, never more than the run holds (`release_top`); the droplets' momentum and enthalpy leave with the
  cells' own (Task 7's ruling 1).
- **A minimal record per column per macro step** for the tests — released mass, number, mode, radius — and nothing of
  Tasks 2–4: the source table's rows per column, per mode, per macro step with `solid_fraction` and `h_kg` and the
  histograms are Task 2's; the CLI's `--slurry-spray` and the skin run's `--removal girin` default are Task 3's; the
  history columns, frames and results of spec §10.2 are Task 4's.
- **Tests** (the brief's): each mode's release from a prescribed column; mass, energy and momentum with spraying on a
  small body; the cap's two depths, its independence of the column's area and its constancy over a macro step's
  sub-steps; liquid first and slurry only with the flag; the element model unchanged (Task 1b's).

**Decision 1 (Asha): how the skins' depths are read where liquid lies on slurry.** Spec §6.7 and §8 settle b and the
regime test (the flowing depth) and the 2026-10-28 decision the viscous film's floor (the flowing depth); the cap's depth
is the admitted run's (2026-10-24). They do not settle what the Rayleigh–Taylor criteria, the molten region the wave
must fit in, and the element model's hold on a film on slurry off Girin's closure (the rigid-substrate amendment of
2026-10-06) read. The element model reads the liquid there; the shared spray physics as it stands reads at least b, so
with b the flowing depth it would read the flowing depth.
- (A) *One flowing film*: everything reads the flowing depth, no column is "on slurry" (the flowing layer rests on the
  rigid front by construction), and the liquid and slurry differ only in what the release may take and the cap's depth.
  No change to the shared physics. On the 100 mm sphere's skins with 0.3 mm of liquid over 3.1 mm of slurry at
  103 m/s² the Rayleigh–Taylor mode takes the 79 columns within 1–10 degrees of the stagnation point from the thick
  branch (droplets 3.05 mm); it reads as inviscid a slurry layer whose viscous depth at the mode's growth time
  (about 7 ms), the Stokes depth √(2ντ*) of 2.1–5.4 mm for Li et al.'s 0.87–5.6 Pa s, is as deep as the layer itself.
- (B) *The element model's split*: the Rayleigh–Taylor criteria and the molten region read the liquid depth, and a
  liquid on slurry off Girin's closure is held from the shear modes, as the 2026-10-06 amendment has them. The shared
  physics gains an optional depth for the Rayleigh–Taylor criteria (the element model never passes it, so it is
  unchanged). On the same state the mode takes none, and 69 of the 79 columns release nothing; the element flight held
  up to 17.9 g this way, at 48 s (fact 171).
- Where neither matters: on sub-plan 21's heated state at 93 s and on 0.3 mm of liquid over 0.9–3.1 mm of slurry at
  39–53 m/s² every windward column takes the thick branch under either.
- Recommended: (B), because the Rayleigh–Taylor mode is a liquid pool's and the slurry beneath is viscous through its
  depth on the mode's growth time, and because Task 5 compares the skins with the element model, which reads the
  depths so.

**Decision 2 (Asha): the slurry cascade under the viscous film.** Spec §8 carries the 2026-10-07 design's cascade over
with its slurry properties: with `--slurry-spray on`, once the liquid is gone, the exposed slurry is stripped at the rate
Girin's model gives for its own properties. An exposed slurry layer takes Girin & Kopyt's thin mode (its own δ_m is
centimetres), and with Li et al.'s viscosity the viscous film damps it at Γτ* = 24–526 on every windward column of the
100 mm sphere at 93 s, against 0.13 for liquid aluminium: the cascade releases nothing (fact 171). With the flag on, the
only change would then be the liquid droplets' cap on the flowing depth.
- (i) Build the cascade as the spec has it, inert: about a hundred lines (Li's viscosity at the slurry layer's own shear
  rate, its δ_m and regime test) and their tests, for a release of nothing.
- (ii) Record in spec §8 that the viscous film makes the cascade inert, and give `--slurry-spray on` the release of
  sub-plan 21 Task 7 (the flowing cells taken after the liquid at the column's own mode's rate) as a sensitivity bound:
  a spec change, small, and the flag no longer the cascade.
- (iii) Defer the cascade to a later task, the flag refused with the flow on until then.
- Recommended: (iii), the flag being off by default and the cascade inert by the measured gate; (i) if the spec's
  cascade must exist for the thesis's statement of it.

---

## The outline (Tasks 2–6, and Task 1 written out above as Tasks 1a, 1b and 1c)

1. **Spraying per column** -- written out above as Tasks 1a (the amendments of 2026-10-22), 1b (the cap per crest)
   and 1c (the spraying from the top cells, which waits for Asha's decisions 1 and 2). The outline's item, as it
   stood:
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
