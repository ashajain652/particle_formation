# work-2026-10-08-spheral-mvp — the reconstructed Step 3 prototype that writes Spheral's MVP frames

**What this is.** The prototype as reconstructed on this machine on 2026-10-07/08 (when `prototype/` was still
git-ignored and the original `proto3/` was missing here), with the post-reconstruction additions Spheral M1 needs.
It produced the 100 mm Scheil frames flight that Spheral M1 is built on (Asha, 2026-10-08: the frames are an MVP
input for the `spheral_frag` pipeline, at the 0.5 s step; accuracy is secondary). Moved in from the git-ignored
`prototype/` on 2026-10-08, when `origin/main` began tracking the prototype's code; renamed on the way:

| here | was | what it is |
|---|---|---|
| `code/` | `prototype/proto3/` | recovered state + additions 01-06 below: **the package that wrote the MVP frames** |
| `recovered/` | `prototype/proto3_recovered/` | the recovered state only |
| `scheil_control/` | `prototype/proto3_scheil_control/` | recovered + addition 01 (control for 02's write-only check) |
| `export_control/` | `prototype/proto3_export_control/` | `code/` without 04 and 06 (control for their write-only check) |
| `rebuild/` | `prototype/rebuild/` | the reconstruction script, the additions as diffs (`rebuild/post/`), the diagnostics |

The text below keeps the old names; `rebuild/rebuild.py` and `rebuild/run_flight.sh` still write to them (run them
with `--out` pointing here). **Provenance:** SHA-256 over `code/reentry_model/*.py` (spheral_frag.fe's
`package_provenance`) = `9629f39c77d1f76c317b8640ae046ae88a0fd576f29fee671df16009b9c3687f`, the value every
`prepare.json` of the frames records. The frames themselves (5.0 GB, git-ignored) are
`reentry_model_output/model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off/`
at the repository root; prepare them with `--fe-package prototype/work-2026-10-08-spheral-mvp/code`.

**Not the Step 3 line of record.** Against the tracked amendments' current state (`../work-2026-10-07-dt-continuum/code`
and its identical copies), this tree lacks the runoff-flux amendment (sub-plan 06, 2026-10-05: the film carries only
its own share of the conjugate layer's flux), the rigid-substrate rule (a film on slurry deeper than the conjugate
depth sprays as deep melt), the Scheil default and `--dt-continuum`; in their place it has addition 05, a cap on the
runoff's emptying number that the runoff-flux amendment makes unnecessary. Additions 02, 03, 04 and 06 (export and
naming only) are what Step 3 must still adopt for Spheral, with 03 changed to name the material always.

## The reconstruction (README as written on 2026-10-07/08)


`prototype/` is git-ignored. The original `proto3/` (the throwaway copy of `reentry_model` with Step 3 applied) and
the plan generator `plan3/make_plan.py` were lost from disk with no backup. Everything here was rebuilt mechanically
from committed documents by `rebuild/rebuild.py`; nothing was hand-pasted. Rebuild with

    python3 prototype/rebuild/rebuild.py --force --post                       # prototype/proto3 (recovered + additions)
    python3 prototype/rebuild/rebuild.py --force --out prototype/proto3_recovered   # the recovered state only
    python3 prototype/rebuild/rebuild.py --force --post --post-materials-only --out prototype/proto3_scheil_control
    python3 prototype/rebuild/rebuild.py --force --post --post-skip 04,06 --out prototype/proto3_export_control

Each run writes `rebuild/rebuild_log*.txt` (every operation, with the line it found each block on).

| tree | what it is |
|---|---|
| `proto3/` | **the recovered state + the post-reconstruction additions below** (use this one) |
| `proto3_recovered/` | the recovered state only (the reference for the fact checks) |
| `proto3_scheil_control/` | recovered + addition 01 only (control for the write-only check of addition 02) |
| `proto3_export_control/` | `proto3/` without the export-only diffs 04 and 06 (control for their write-only check, 2026-10-08) |

Run the prototype from inside its directory, so that `reentry_model` resolves to the prototype's (verified:
`reentry_model.__file__` is `prototype/proto3/reentry_model/__init__.py`; main has no `spray.py`):

    cd prototype/proto3 && /Users/ashajain/miniforge3/envs/drama_env/bin/python -m pytest -m "not drama and not reference" -q \
        --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py
    rebuild/run_flight.sh <tree> 100|50 <outdir> [extra flags]                # one physics flight (melt_verification's `physics` mode)

## What was reconstructed, from where (the recovered state: 114 files, the count of the 2026-10-03 manifest)

Stages, in order (`rebuild.py`, `STAGES`):

1. **base** — `git archive ee70f6b` of `reentry_model/ tests/ data/ analysis/` (93 files). `ee70f6b` adds the Step 3
   plan; `reentry_model/` and `tests/` are unchanged from 07b2856/e376a82 up to 187cd51, so it is the Step 2 state the
   plan was written against. The wrapper scripts are not in the copy (the prototype's own test command ignores their tests).
2. **task1** — Task 1 is superseded by sub-plan 01, which was implemented in the prototype and ported to main in
   **187cd51**: `reentry_model/mesh.py` and `tests/test_reentry_model_mesh.py` from 187cd51, plus its `band = 0` pins of
   `tests/conftest.py` and `tests/test_reentry_model_reference_thermal.py` (judgement call, see below).
3. **plan** — Tasks 2–15 of `docs/superpowers/plans/2026-09-20-melt-spraying.md`, every code block addressed by the
   line its fence opens on and a sentinel in the prose before it (plan lines 796–8107; full list in `rebuild.py`
   `stage_plan`): replace/create blocks (material, thermal ×3, dispersion, surface_flow, film, spray, girin_case,
   body, coupled, viz, cli, the analysis scripts, tests), "append" blocks (two blank lines before, as the later diffs'
   line numbers confirm), the targeted edits (gas.py's four, aero.py's five, trajectory.py, sesam_io.py, compare.py,
   test_thermal, test_fenicsx, test_aero ×2, cli.py's `--consistent-mass`, the package docstring of Task 15), and the
   two bash blocks (the material JSONs, plan:796; `girin_dispersion.json`, plan:2226). Task 15's README/assumptions/spec
   edits are repository documents, not prototype code, and are not applied. The sub-plans' task bodies are identical to
   the master plan's blocks for Tasks 3–15 (checked block by block).
4. **band13** — "the dense-band line" that sub-plan 13's 2026-10-02 diff says the prototype's cli.py had (3 lines,
   reconstructed; see below).
5. **fact37** — measured fact 37 (2026-09-27, `MeltingBody.last_face_ids`/`on_current_surface` and the VTK writer
   reading through it): sub-plan 09 lines 33–52, sub-plan 10 lines 28–37; the writer's docstring from the context of the
   2026-10-02 coupled.py diff.
6. **d1002** — the nine diff blocks of the amendment of 2026-10-02 (deep runoff; `delta_m`, `deep_thickness` in the
   frames): 06:49, 06:90, 09:174, 09:583, 10:81, 10:132, 13:56, 13:121, 14:73.
7. **d1003** — the eight of 2026-10-03 (molten cascade): 09:890, 09:1074, 09:1268, 10:216, 10:241, 13:187, 13:254, 14:141.
8. **d1005** — the three of 2026-10-05 (`--seed`, default 12345): 13:340, 13:445, 14:224.

**The strongest evidence of fidelity:** all 20 diff blocks are applied by a strict applier (`apply_diff`: every hunk's
old lines must match at the exact line its header states) and **every hunk of all 20 matches at its stated line, with
no offset**, in body.py, coupled.py, cli.py, film.py, melt_sensitivity.py and four test files. That pins the line count
of everything before each hunk, including the fact-37 and band-line placements (before those were placed right, the
applier reported offsets of −1, −2 and −3, which is how their size was found). 25 of the 34 master-plan file blocks
are still verbatim in the latest tree; the other 9 are exactly the files the documented later changes touch
(`rebuild/plan_blocks_vs_latest.txt`). Sub-plan 02's generator writes `AA7075.json` and `AA7075_range.json`
byte-identical to the ones Task 2's block wrote, as its text says it should.

**Not included, by the documents' own account of the prototype** (facts 46 and 54, sub-plan 09): sub-plan 02's
materials (tested only in throwaway copies; added below as a post-reconstruction addition, for Spheral), sub-plan 09's
2026-09-27 changes (`PHI_DEATH = 0.50`, shape consumers on the derived surface: plans only — the prototype keeps
`PHI_DEATH = 0.05`), Tasks 16–17 (plans only), Task 11's melting SESAM references (never in the prototype).

## Judgement calls and what could not be recovered

- **`reentry_model/data/atdb_disc.json`** (Task 9) is extracted from DRAMA's `ATDB_CYLINDER.nc`; DRAMA is not installed.
  Written from the values the plan states (cd_continuum at all six Mach numbers, cd_free_molecular at Ma 10,
  heat_flux_factor_continuum); the other 11 values are `null`. Read only by two tests, which fail on the nulls
  (`test_disc_aerothermal_table`, `test_drag_shape_factor_between_the_two_atdb_endpoints`); no model code reads it.
  `rebuild.py` runs the plan's extraction script instead when DRAMA is present.
- **The dense-band line (cli.py, 3 lines)**: reconstructed as `band = mesh.DEFAULT_BAND if melting else 0.0` (+ two
  comment lines). Known: its size and region (the diffs) and that melting runs get the 15 mm band (177 363 tets,
  18 830 patches, reproduced). Non-melting CLI runs pinned to band 0 mirror 187cd51 — not documented for the prototype.
- **conftest.py / reference_thermal pins (band = 0)**: with the default band the Carslaw–Jaeger test fails (5.005 K
  against its 5 K margin, exactly 187cd51's recorded failure), while the prototype's tier after sub-plan 01 had only the
  five known failures — so the prototype must have pinned them.
- **`import pytest` in tests/test_reentry_model_data.py**: Task 11's appended test uses `pytest.approx` and the file
  never imported pytest; the prototype's test passed, so it had the import. Added.
- **Fact 37's test lines**: 14 lines in `test_coupled_melting_run_and_writers` (size from the diff offsets). The first
  two are the 2026-10-02 diff's context; the twelve after them are **reconstructed from fact 37's prose** (a death in the
  frame's step, p_w positive on most windward patches and never above the history's stagnation value, Spearman < −0.9
  against the angle). The step-(iii) `last_face_ids` line's comment is reconstructed (code is the snippet's).
- **Docstring wording** of aero.py's class (Task 9, step 6) is wrapped at 120 columns; the plan gives one sentence.
- **plan3/make_plan.py**: not reconstructed. Its prose parts (`part_header.py` …) exist only as their concatenation in
  the committed plan, so a "regeneration" would be circular, and the committed plan is not the latest prototype anyway
  (57c1129: the generator was 533 lines ahead on 2026-09-28, and three amendments followed). The block-by-block check
  above is the non-circular part of that proof.

## Verification (2026-10-07, drama_env rebuilt 2026-10-06 without pyDRAMA; no fenicsx_env)

**Unit tier**, recovered state: 232 passed, 1 skipped, 4 failed, 3 errors (132 s). Expected from the documents
(2026-10-05): 234 passed, 1 skipped + the five known failures/errors of the missing melting references. Here: the same five, plus the two tests that read the
partial atdb_disc.json — 232 = 234 − 2, and the collected total (239) matches. `test_full_steps_with_the_molten_cascade_keep_the_books_exact`
is **flaky** at its 1e-12 mass tolerance (relative error 1.0–1.2e-12): it passes with numpy seeded 0 or 1 before the
test, fails with 12345 or unseeded in isolation, and passed in one of three full-tier runs. FEniCSx tests skip
(`fenicsx_env` was not recreated). With the additions (`proto3/`): 237 passed + the same (5 more material tests, 1 new
coupled test).

**Flights** (`reentry_model_output/proto3_verification/`, `AA7075_range`, default step and seed 12345 unless stated).
Recorded = the fact in `melt-spraying-subplans/00-shared-context.md`.

| configuration / quantity | recorded | reproduced |
|---|---|---|
| **50 mm, default (deep on, cascade on)** — facts 60, 64, 67 | | |
| demise | 199.0 s, 70.7 km | 199.0 s, 70.74 km |
| sprayed mass | 0.1770 kg | 0.17696 kg |
| droplets; median r by number / by mass; largest | 1.92e6; 173 / 227 µm; 1.81 mm | 1.917e6; 173.3 / 227.5 µm; 1.811 mm |
| cascade: mass, active steps, max passes, capped | 95.6 g, 45 of 398, 10, 0 | 95.64 g, 45 of 398, 10, 0 |
| enthalpy per kg sprayed; energy balance | 1.089 MJ/kg; 7.4e-11 | 1.0894 MJ/kg; 7.39e-11 |
| macro steps; Newton iterations/step; source rows | 398; 2.236; 21 010 | 398; 2.2362; 21 010 |
| **100 mm, deep runoff off, cascade off, whole flight** — fact 49 | | |
| melt onset | 25.5 s (73.96 km, fact 65) | 25.5 s, 73.964 km |
| landing; mass on landing | 582.5 s; 0.453 kg (30.8 %) | 583.6 s (ground event; 1168 steps); 0.4516 kg (30.68 %) |
| sprayed mass | 1.019 kg | 1.0204 kg |
| spraying ends | 225.0 s | 207.0 s (last nonzero release; 195.5 s above 1 µg) |
| re-solidified | 6.6 g | 5.76 g |
| median droplet radius (by number) | 180.2 µm | 180.15 µm |
| **same run, to 120 s** — facts 49, 50, 57 | | |
| sprayed; mass at 120 s | 1.018 kg; 0.454 kg | 1.0201 kg; 0.4519 kg |
| droplets; median r number / mass | 1.72e7; 180 / 189 µm | 1.720e7; 180.15 / 189.08 µm |
| RT release; re-solidified; largest droplet | 53.0 g; 6.2 g; 4.38 mm | 51.98 g; 5.58 g; 4.42 mm |
| thin-branch share of sprayed; thick-branch median share | 9.4 %; 12.4 % | 9.37 %; 12.39 % |
| deepest molten layer; its median per step | 20.9 mm; 10.55 mm | 20.89 mm; 10.53 mm |
| **100 mm, deep runoff on, cascade off, whole flight** — fact 49 | | |
| landing; mass | 602.7 s; 0.330 kg (22.4 %) | 599.2 s; 0.3311 kg (22.5 %) |
| sprayed; re-solidified | 1.142 kg; 21.0 g | 1.1409 kg; 21.47 g |
| deep: taken / became film / froze back; account on landing | 285.1 / 277.1 / 8.0 g; empty | 284.3 / 275.0 / 9.3 g; empty |
| deep liquid surfaced and sprayed by | 229.5 s | deep account empty from 195 s; last release 222.0 s |
| median r; energy balance | 180.2 µm; −8.8e-10 | 179.6 µm; −5.4e-10 |
| to 120 s: sprayed; mass; deep peak; RT; largest; re-solidified | 1.091; 0.381 kg; 117 g @96.5 s; 72.5 g; 5.66 mm; 11.2 g | 1.0860; 0.3861 kg; 113.4 g @95.5 s; 74.5 g; 6.21 mm; 12.35 g |
| thick-branch median; molten-layer median per step | 15.5 %; 10.46 mm | 15.3 %; 10.44 mm |
| **100 mm, default, to 120 s** — facts 59, 65 (range over seeds 12345, 1, 2, 3) | | |
| seed 12345 | completes | **fails at t = 110.0 s** (exit 1: singular runoff matrix → NaN film mass) |
| seed 2 | completes | **fails** the same way |
| seeds 1 / 3: sprayed | 1.0600–1.0614 kg | 1.0606 / 1.0591 kg |
| mass at 120 s | 0.4107–0.4120 kg | 0.4114 / 0.4129 kg |
| droplets; median r number; by mass | 1.762–1.798e7; 179.6–180.1; 190.0–190.3 µm | 1.788 / 1.786e7; 179.89 / 179.85; 190.10 / 190.02 µm |
| re-solidified; RT release | 6.38–6.46 g; 56.2–60.0 g | 6.28 / 5.98 g; 54.05 / 56.20 g |
| largest droplet; thin share; onset | 5.058 mm; 10.2–10.3 %; 73.96 km | 5.058 / 5.058 mm; 10.24 / 10.28 %; 73.96 km |
| cascade mass; deep mobilised; max passes (fact 59) | 157 g; 91 g; 6 | 157.3 / 156.6 g; 91.1 / 90.7 g; 6 |
| **same, with post/05 (2026-10-08)**: seeds 12345 / 1 / 2 / 3 | all complete | **all four complete** to 120 s (`to120_fixed/`) |
| sprayed; mass at 120 s | 1.0600–1.0614; 0.4107–0.4120 kg | 1.0611 / 1.0610 / 1.0610 / 1.0590; 0.4109 / 0.4111 / 0.4110 / 0.4130 kg |
| droplets; median r number; by mass | 1.762–1.798e7; 179.6–180.1; 190.0–190.3 µm | 1.766 / 1.764 / 1.784 / 1.751e7; 179.98 / 180.15 / 180.01 / 180.19; 190.26 / 190.21 / 189.99 / 190.19 µm |
| re-solidified; RT release | 6.38–6.46 g; 56.2–60.0 g | 6.28 / 6.64 / 6.21 / 6.66 g; 60.7 / 55.4 / 52.7 / 56.3 g |
| largest; thin share; cascade; deep mobilised; energy | 5.058 mm; 10.2–10.3 %; 157 g; 91 g | 5.058 mm all; 10.22 / 10.20 / 10.29 / 10.20 %; 157.6 / 156.7 / 157.4 / 157.4 g; 91.1 / 91.0 / 91.1 / 91.3 g; −7.7 / −8.6 / −9.8 / −8.8e-11 |

**Reading.** The 50 mm flight reproduces every recorded figure to every printed digit, energy balance and Newton
count included (fact 65: this flight barely amplifies round-off). The 100 mm flights reproduce every recorded figure
to within roughly the run-to-run scatter between seeds that facts 52 and 65 measured for the same model versions
(sprayed mass within 0.2–0.5 %, masses within 0.5 %, droplet counts and median radii within 0.2 %, re-solidified mass
and RT release within 10 %), with the threshold-sensitive end-of-flight figures off by more (landing ±1–3 s, last
spray 207 against 225 s). They are **not** bit-identical to the recorded seeded runs, although facts 64/65 say a
seed reproduces the run exactly. Best understanding: the 2026-10-06 rebuild of `drama_env` (fresh conda-forge
binaries of the same numpy 2.5.3 / scipy 1.18.1 / pyamg 5.3.0 / gmsh 4.15.2) changes the last bits of the solves, and
the 100 mm flight's melting thresholds amplify them (fact 62) where the 50 mm flight's do not. A reconstruction error
confined to code that runs only on the 100 mm flight (Girin's closure, the thick branch, the deep runoff, the RT mode)
cannot be excluded by measurement — but that code is verbatim plan code whose context the diffs pin line for line.
Two of four seeds of the default model **crashed at ~110 s** (`spsolve` reported the runoff matrix exactly singular):
the *film's* runoff solve, not the deep stage -- diagnosed and fixed 2026-10-08 by post/05 (below); all four seeds
now complete.

## Post-reconstruction additions (in `proto3/` only; NOT part of the recovered state)

1. **Sub-plan 02's materials** (requested 2026-10-07 for Spheral M1, spec 2026-10-02 §13.1): its generator
   (02-material-properties.md:134; AA7075/AA7075_range unchanged byte for byte), `material.py` (02:334) and its two
   test blocks (02:238, 02:306). `AA7075_scheil` as the file has it: solid density 2813 kg/m³; **liquid density
   2400 kg/m³** (μ 1.3e-3 Pa s, σ 0.80 N/m); **solidus 750 K, liquidus 908 K** in the file (`solidusTemperature`,
   `liquidusTemperature`; `Material.T_solidus` is the eutectic ramp foot, 748 K); Scheil k = 0.4, T_pure 933 K; latent
   heat 390 kJ/kg; T_feed 910 K.
2. **`post/02-frame-loads-every-step-and-v-hat.diff`** (export only):
   - `MeltingBody.last_flow_step` keeps the step's own surface-flow evaluation, which `melt_step` makes on every step
     with an aero state, film or not; the writer takes `closure`, `p_w` and `tau` from it, so frames before melt onset
     and after spraying ends carry the wall loads (before: closure 1, p_w = tau = 0 on every step without film,
     because `_film_and_spray` sets `last_flow = None` there). The defaults (closure 1, p_w and tau 0) now remain only
     where the flow was not evaluated: **frame 0** (written before the first melt step), runs with
     `--removal instant` (never evaluates it), and faces the step's own element deaths exposed (no evaluation on them).
     `kn_local`, `we_s`, `r_droplet`, `release_rate` are unchanged (spray-step fields).
   - History column `p_w_stag_step_Pa`: the stagnation-patch wall pressure of that evaluation, every step (NaN only on
     row 0). The existing `p_w_stag_Pa` is the spray step's (NaN before the first film, stale after the last).
   - Run JSON `settings.v_hat_body` = the direction of motion in the body (mesh) frame, [1, 0, 0], and
     `settings.freestream_velocity_direction_body` = −v_hat (the oncoming flow relative to the body).
   - Tests: one new coupled test (frames and the column on a cold body with no film), one assertion in the CLI test.
   - **Write-only, measured**: the 100 mm AA7075_scheil flight below run with and without this diff
     (`proto3_scheil_control/`, no frames): all 1260 rows of all 86 shared history columns identical, all 172 728
     source-table rows identical, every result field identical but `runtime_s` and `n_frames`. By the code: the stored
     object is the `flow` already computed and passed on; nothing reads `last_flow_step` but the writer and
     `melt_stats` (which only reports it).

3. **`post/03-run-name-encodes-the-material.diff`** (2026-10-08; naming only). `model_run_name` gains `material_name`: a
   finite-element run whose material is not its mode's default (`AA7075_range` melting, `AA7075_nomelt` otherwise)
   gets `_material-<name>` after `_fem-<heating>_melt-<removal>` and before `_deeprunoff-off`/`_moltencascade-off`/
   `_seed-<n>` (a file is named by its basename without `.json`). This follows the existing optional parts (`_<key>-<value>`,
   default omitted, so every recorded `AA7075_range` run keeps its name). The Scheil frames run is now
   `..._melt-girin_material-AA7075_scheil_deeprunoff-off`. Renamed outputs (directory, `.csv`, `.json`, with `run_name`
   and paths inside the JSON rewritten): `reentry_model_output/proto3_verification/RENAMED_2026-10-08.txt`
   (`scheil_control/` and four `diag_energy/` runs; no `AA7075_range` run changes name). Test: 6 assertions in
   `test_run_name_with_melt`.
4. **`post/04-derived-surface-normals-in-the-frames.diff`** (2026-10-08; export only; Spheral M1 decision 9).
   `surface_<k>.vtp` cell field **`n_derived`** (3 components): each patch's outward unit normal on sub-plan 01's derived
   surface, `mesh.surface(derived=True).smoothed_normals()` (spec 2026-09-27 §8 step 4: Taubin, 8 passes, clamp 0.35
   edge). The derived surface's patches are the vtp's patches **one to one, same order** (checked by face id in the
   writer, which raises otherwise). Nothing in the prototype displaces `derived_points` (sub-plan 16 is plans only), so
   the derived facets coincide with the staircase's and `n_derived` is the smoothing alone. The vtp triangles keep the
   face table's node order, so their winding is **inward on 8.2 % of the faces at k = 100, 15.2 % at k = 200, 18.0 % at
   k = 300**; `n_derived` is outward by construction (oriented like `surface.normals`, by the patch's opposite vertex):
   measured, `n_derived` · (winding normal) < 0 on exactly the inward-wound faces and `n_derived` · (outward facet
   normal) < 0 on none. Angle to the outward facet normal: median 0.1° / 0.2° / 0.7°, 90th percentile 36° / 54° / 61°
   (k = 100 / 200 / 300). Effect on Task 6 (`spheral_frag.geometry.layer_depths`, `diag/depth_normals.py`): rays that
   leave the body before f_l falls to 0.5 drop from **42.8 % to 18.7 %** of the 3,436 marched patches at k = 100 and
   from 58.6 % to 21.5 % at k = 200 (no wrong-way normals). Tests: two coupled tests.
5. **`post/05-runoff-emptying-number-bound.diff`** (2026-10-08; **physics**, degenerate rows only by intent; see the
   findings below). `film.EMPTYING_MAX = 1e12`: in `Runoff.transport` a donor patch whose emptying number
   dt_s · Σ c exceeds it has its edge coefficients scaled down to it. Test: two nearly dry patches feeding each other
   (reproduces the exactly singular matrix before, finite/positive/conservative after; an ordinary film caps nothing).
6. **`post/06-frames-flow-record-on-exposed-faces.diff`** (2026-10-08; export only). Faces the step's element deaths
   exposed (not on the surface the step's flow was evaluated on) took the defaults (p_w 0, tau 0, closure 1, NaN) --
   **1,112 of 18,616 faces at k = 100** (10 % of the area; the molten cascade kills 812 elements in that step, 3
   passes), 970 of 8,554 at k = 200. Measured: on the 2026-10-07 frames the zero faces were exactly the faces not in
   the previous frame (1,112 = 1,112, 985 = 985), so the mapping by face id was right and nothing was mis-indexed. Now
   `MeltingBody.flow_for_the_record()` evaluates the same `SurfaceFlow` at the step's aero state (`last_state_step`) on
   the current surface, and the writer takes closure, kn_local, p_w, tau and delta_m from it on those faces only; new
   cell field **`flow_eval`** (1 the step's own evaluation, 2 the record's, 0 none: frame 0 and `--removal instant`).
   Spray fields keep their defaults there (nothing sprayed from a face that did not exist). On the faces both cover,
   the record's p_w equals the step's bit for bit (p_w depends on theta and the state only); its tau differs because
   the nose radius is refit after the deaths (median 1-53 % on the 880 K coarse-mesh test body). Over the rerun's
   1,255 frames: 172,644 faces carry `flow_eval` 2 (198 frames), 0 faces with an evaluation carry p_w = 0; frame 0
   alone has `flow_eval` 0. Tests: the coupled melting test (flow_eval present, p_w > 0 on every face with n_x > 0.2),
   the cold-body test (no deaths: all 1), the conjugate-depth test's expectation.
   - **Write-only, measured for 04 + 06 together**: the rerun below against `proto3_export_control/` (no frames): all
     1,256 rows of all 87 history columns identical, all 172,640 source-table rows identical (NaN-aware), every result
     field identical but `runtime_s` and `n_frames`.

**Unit tier with 03-06** (`rebuild/unit_tier_post_20261008.txt`): 238 passed, 1 skipped, 5 failed, 3 errors -- the
recorded 237 + the new film test; the failures/errors are the same eight (the five missing-melting-reference ones, the
two partial-`atdb_disc.json` ones and the flaky `test_full_steps_with_the_molten_cascade_keep_the_books_exact`).

### Findings of 2026-10-08 (diagnostic drivers in `rebuild/diag/`, outputs in `reentry_model_output/diag_*`)

- **The runoff crash** (`diag/singular_trap.py`, `diag/nan_trap.py`): the singular matrix is the **film's**
  `Runoff.transport` (called from `_film_and_spray`), not the deep stage. On the thick branch the shear-driven flux
  V_s δ_m/2 does not vanish with the film, so the emptying rate q/(A b) of a nearly dry patch grows like 1/b: in the
  failing solve (seed 12345, recovered tree, t ≈ 110 s, 8,620 patches) one patch holding 5.3e-39 kg had
  dt_s · out = **5.1e34**, others 1e17-4e21 (no zero diagonal, no NaN, no isolated patch: diagonal ≥ 1 everywhere);
  1 + dt_s · out rounds to dt_s · out, the columns of nearly dry patches that feed each other sum to exactly zero, and
  SuperLU finds the matrix exactly singular. Such rows are **routine**: with post/05, in the four 120 s flights the
  bound holds 96,646-96,871 donor rows in ~400 of the 760 film solves (largest emptying number 2e49-7e61); the deep
  stage never exceeds 5.6e7 and is never capped (`to120_fixed/capped_seed*.json`). Hence post/05 changes every
  default 100 mm flight from 49.5 s on (seeds 1 and 3: histories bit-identical to the pre-fix runs up to row 98,
  t = 49.0 s; first difference at 49.5 s in `thick_branch_fraction`, `rt_bounded_growth_ms`), and the end figures move
  within the seed-to-seed scatter (table above). A bound applied only when SuperLU complains would keep the
  non-failing runs bit-identical but rests on its detection; the always-on bound keeps every solve well conditioned.
- **The Scheil energy residual 7.0e-9** (`diag/energy_trace.py --decomp`, reproducing the run's 6.977757923139848e-09
  bit for bit): **100 % the CG residual of the conduction solve** (CG rtol 1e-10 of |b|, `skfem_backend._solve`);
  storage (true nodal enthalpy increment vs the one the linear system enforced) -9e-15, radiation linearisation
  -7e-12, melt bookkeeping (feed, cascade, spray, deaths) -4e-15. Not the eutectic ramp, not the cascade (cascade off:
  6.2e-9 with 2 threads, see below), not the export (history identical with and without). It accumulates in the
  cooling phase 250-400 s (+7.4e-9), where each step takes one Newton iteration and its CG stops after ~3 iterations
  just under the tolerance with a systematic sign (91 % of steps positive, 3.9e-5 J per step against a net 135 J per
  step, 3e-7 of it); the range material there needs 2 CG iterations and lands far below the tolerance (-2e-7 J per
  step). Same flags with `AA7075_range`: -1.20e-9 (so the Scheil flight is 6x, not 50x, the comparable range flight;
  the recorded ~1e-10 are the 50 mm and 120 s flights). Newton tolerance 1e-8 instead of 1e-6 (more iterations, so
  the last CG starts closer): Scheil 4.7e-10, range -1.2e-9. In absolute terms 9.4 mJ over the flight, against
  1.35 MJ absorbed and 1.23 MJ removed with the spray: **an expected solver-tolerance effect, not a bug**; no
  tolerance was changed.
- **Thread count changes the bits.** The same Scheil flight with `OMP_NUM_THREADS=2` gives 5.9e-9 (default threads:
  6.98e-9) and the default seed-12345 flight does not crash at 110 s with 2 threads where it does with the default:
  the BLAS-threaded reductions are a second source of last-bit differences, besides the 2026-10-06 rebuild of
  `drama_env`. Every run in this README uses the default thread count.
- **The frame's wall loads against the trajectory's drag** (`diag/frame_drag.py`, coordinator's question; nothing
  changed). They are independent models by design: the trajectory's drag is q · C_D · A_ref with C_D from SESAM's
  sphere tables bridged in Kn and the continuum part scaled by `drag_shape_factor` (the convex hull's modified
  Newtonian C_D over the sphere's), A_ref the hull's projected area; the surface flow's p_w is modified Newtonian +
  Prandtl-Meyer from the equilibrium stagnation pressure, a function of each **staircase** patch's theta, and tau the
  boundary-layer shear. Σ A [(p_w − p_inf) cos θ + τ sin θ] on the rerun's frames: k = 100 (50 s): **14.43 N**
  (pressure 13.94, shear 0.50) vs the history's **19.52 N**; k = 200 (100 s): **27.15 N** vs **42.59 N**. On the
  2026-10-07 frames (exposed faces at 0) the same sum was 10.96 / 20.00 N; filling them by p_w(θ) interpolation gave
  14.43 / 27.13 N, which the record now reproduces. The remaining gap is the staircase: its facets scatter theta by up
  to 45°, and p_w ∝ cos²θ per facet, so the projected-area-weighted pressure falls well below the silhouette's. Reading
  p_w at the derived normal's theta on each facet's exact projected area A cos θ_stair gives 21.00 N (k = 100) and
  36.29 N (k = 200), against 19.52 / 42.59 N; modified-Newtonian sphere reference q (Cp_stag/2) π R_t² = 16.57 / 25.59 N
  (Cp_stag 1.94 / 1.91; the history's shape factor 1.25 / 1.74). Integrating with the derived normals but the
  staircase areas over-counts (31.2 / 51.8 N) because the staircase area exceeds the smooth area.

## The frames for Spheral M1

`reentry_model_output/model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off/`
(+ `.csv`, `.json`), rerun 2026-10-08 from `proto3/` with post/01-06 (log `proto3_verification/S2_d100_scheil_frames.log`):
`--diameter 100 --altitude 77.500133 --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on
--heating physics --material AA7075_scheil --removal girin --deep-runoff off --molten-cascade on --frames-every 1`,
dt 0.5 s, seed 12345, default thread count. **5.0 GB, 1255 frames** (`vtk/field_<k>.vtu`, `vtk/surface_<k>.vtp`,
k = 0-1254, every macro step from 0 to the ground at 626.5 s, with `field.pvd`/`surface.pvd`; `spheral_frag.frames.read_fe_run`
reads it with no contract item missing). The 2026-10-07 run (old name, 1260 frames) was deleted after this one completed.

| flight fact | 2026-10-07 run | 2026-10-08 rerun |
|---|---|---|
| melt onset | 24.0 s, 74.19 km | 24.0 s, 74.19 km |
| last nonzero spray release | 214.5 s | 213.0 s |
| landing; mass on landing | 629.2 s (1259 steps); 0.3044 kg (20.7 %) | 626.5 s (1254 steps); 0.3032 kg (20.6 %) |
| sprayed; droplets | 1.1676 kg; 1.798e7 | 1.1688 kg; 1.811e7 |
| median r by number / by mass | 182.30 / 191.45 µm | 182.30 / 191.51 µm |
| cascade; re-solidified; source rows | 325.0 g; 9.24 g; 172,728 | 326.5 g; 7.37 g; 172,640 |
| energy balance | 6.98e-9 | 5.5e-9 |

The difference between the two is post/05 (the only physics change; the film's runoff bound engages on this flight
as on the range flights); 03, 04 and 06 are naming and export (write-only, measured above).

- `field_<k>.vtu`: active tetrahedra (177 363 at k = 0), all 33 355 mesh points; point `T`, `liquid_fraction`; cell `phi`.
- `surface_<k>.vtp`: the current surface triangles (18 830 at k = 0; winding not oriented -- use `n_derived`), all mesh
  points; point `T`; cell `q_conv`, `q_rad`, `T_patch`, `n_derived` [outward unit normal of the derived surface, 3
  components], `film_thickness` [m], `film_T` [K], `we_s`, `flow_eval` [1 step's flow, 2 record, 0 none], `closure`,
  `kn_local`, `p_w` [Pa], `tau` [Pa],
  `r_droplet` [m], `release_rate` [kg/m² released in the step], `delta_m` [m, NaN off Girin's closure],
  `deep_thickness` [m, zero here: deep runoff off].
- History CSV (87 columns): `time_s`, `altitude_km`, `velocity_kms`, `flight_path_deg`, `heading_deg`, `lat_deg`,
  `lon_deg`, `downrange_km`, `mach`, `knudsen`, `density_kgm3`, `dynamic_pressure_Pa`, `load_factor_g` (the
  deceleration, drag / m g0), `mass_kg`, `p_w_stag_Pa`, `p_w_stag_step_Pa`, and the rest of the Step 2/3 columns.
- Run JSON: `settings.v_hat_body`, `settings.freestream_velocity_direction_body`, `settings.seed`, material, flags.
