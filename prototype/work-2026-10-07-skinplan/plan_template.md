# Sub-plan 18 — the melt-layer skin: stage 0 (the grid gate) and stage 1 (skins with melting) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first (required reading for every Step 3 sub-plan), then the spec.

**Goal:** Measure whether Step 3's converged droplet population depends on the surface mesh (stage 0), then give every
melting surface patch a one-dimensional skin that resolves the melt layer and feeds the film continuously, coupled to
the 3D conduction through an implicit, energy-exact interface (stage 1).

**Architecture:** A new module `reentry_model/skin.py` holds every skin as arrays (cells of metal mass and enthalpy along
the patch normal) and solves them together on short sub-steps. `MeltingBody` creates skins where melting begins, runs a
trial pass to get each skin's base-flux law, solves the 3D step with that law as an affine boundary condition (both
thermal backends gain it), runs the real pass, books the interface mismatch, and hands skins on at element deaths.
`--surface-model elements` (the default in stage 1) keeps every earlier run bit for bit.

**Tech Stack:** Python 3.12 (`drama_env`), NumPy, SciPy, scikit-fem backend, FEniCSx 0.11 backend (`fenicsx_env`), pytest.

**Spec:** `docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md` (approved and revised 2026-10-07). This
sub-plan covers its stages 0 and 1 (§12.1). Stages 2 (the film and spray on the sub-steps) and 3 (the slurry cascade)
get sub-plan 19, written after stage 1's measurements fix the presets and the interface tolerance this plan measures —
the spec gates each stage on the previous one's measurement, and stage 0 may change whether stages 1–3 proceed at all.

**How this plan's code was checked (2026-10-07).** Every code block and diff block of Tasks 2–8 was first applied to a
scratch copy of `prototype/work-2026-10-07-scheil-default/code/` and tested there (unit tier 252 to 286 passed, the same
five known failures), then rebased onto `prototype/work-2026-10-07-dt-continuum/code/` (facts 97–100), the latest
committed state, where only Task 8's CLI diff needed merging by hand. On that base: the skin module's 17 tests, the 2 new
scikit-fem backend tests, the 2 new FEniCSx tests, the 8 body tests, the 7 CLI tests, and the whole unit tier, which
goes {{UNIT_REBASE}}. {{BITID}} The blocks below are generated from the tested files, not retyped, and replaying them in
order on a fresh copy of the dt-continuum copy reproduces the tested copy file for file. Measurements made while checking are quoted
where they set a tolerance or explain a choice; they were each made once, on the cases named.

## Global Constraints

- Work in a throwaway copy of the prototype, never in `prototype/proto3/`, the package `reentry_model/` at the repo
  root, the specs or the monolithic plan. The prototype is git-ignored: there are **no commits during the build**.
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Record its counts before any change and after every task. The only acceptable failures are the five known ones that
  need Task 11's melting SESAM references (the `--dt-continuum` amendment's baseline: 258 passed, 1 skipped, those five;
  FEniCSx tier 12 passed).
- `--surface-model elements` (the `MeltSettings` and CLI default) must reproduce the copy's earlier flights **bit for
  bit** in every history column they share, every source-table row and every result field but the run time. New
  history columns and result fields are added; none is changed.
- A run with skins is **bit-identical to the element model until the first skin is created**.
- Bit-identity is always under the same seed and call sequence: pyamg's AMG set-up draws on numpy's global generator, so
  two solvers stepped from identical input differ by about 2e-11 K unless each step is seeded alike (measured while
  planning). The model's `--seed` (fact 62) does this for a run; tests seed before each compared step.
- Mass books exact (relative 1e-12) and the coupled energy balance (`energy_balance_residual`) below 1e-8 in every test
  that melts; flights below 1e-8.
- Melting default material: `AA7075_scheil` (fact 88). Seed: default 12345 (fact 62). Default macro step after the
  continuum switch: 0.0125 s (`--dt-continuum`, facts 97–100), so a default 100 mm flight to 120 s takes about 2.3 h;
  every skin run states its post-switch step explicitly.
- Exit codes: 0 ok, 1 model failure, 2 bad arguments or a missing optional library. Every new CLI value is validated;
  malformed values exit 2.
- Run names encode the configuration: a skin run appends `_surface-skin` and its non-default settings' suffix (Task 8).
- The CLI default stays `--surface-model elements` in stage 1; it becomes `skin` only when stage 2 passes the spec's
  acceptance (§12.3), so no default run changes before the skin is verified.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).
- Results are recorded only in `docs/superpowers/plans/melt-spraying-subplans/`: facts in the shared context, tested
  code as diffs in the affected sub-plans (Task 9). Nothing outside that folder is edited.

## Review Focus

1. A skin whose element cannot supply its target thickness (an element near `PHI_DEATH`): the skin is made thin, never
   negative, counted, and the element dies at `PHI_DEATH` if it reaches it. Pinned by
   `test_a_skin_its_element_cannot_supply_is_made_thin_and_counted` (Task 7).
2. Liquid on a patch with no skin (it ran there during the previous step): it rides the 3D nodes until the next step
   creates the patch's skin, and its energy moves onto the skin's top exactly. Pinned by
   `test_liquid_on_a_cold_patch_moves_onto_its_new_skin_with_its_energy` (Task 7).
3. A macro step whose interface mismatch exceeds the tolerance (heating that changes fast, or a skin feeding — the base
   law is good to about 1 % per kelvin while a skin feeds, measured while planning): the step repeats once and the books
   stay exact. Pinned by `test_a_large_interface_mismatch_repeats_the_step` (Task 7).
4. A death whose element exposes no face (orphan), or skins left without mass: the hand-over falls back to the nearest
   patch and no empty skin survives. Pinned by `test_hand_over_falls_back_and_drops_empty_skins` (Task 2) and
   `test_deaths_hand_skins_to_the_faces_they_expose` (Task 7).
5. A whole skin turning liquid in one sub-step with nothing left to draw: no NaN, no negative mass, the liquid goes to
   the patch's liquid, and the next sub-step still solves. Pinned by `test_an_exhausted_skin_empties_without_nan`
   (Task 4).

## Where the work is done, and how it reaches the sub-plans

The build happens in `prototype/work-<date>-skin/code/`, a copy of the latest tested copy (Task 1). Paths in the tasks
are relative to that `code/` directory. Tasks 6–8 give their changes as unified diff blocks (`a/` and `b/` prefixes):
apply each with `patch -p1` from `code/`, or by hand where a later amendment has moved the context lines (the changes
are additive). When stage 1's tests pass and its measurements are made (Task 9), the tested changes are written into
the affected sub-plans as dated amendment sections with diff blocks against the copy as it stood before the change, and
the measurements into the shared context as new facts — the repository's rule for every Step 3 change.

## File structure

| File (in the copy) | Responsibility | Task |
|---|---|---|
| `reentry_model/skin.py` (new) | `SkinSettings`, `SkinBook`, `SkinField`: storage, rezoning, the top node, depths, hand-over, the implicit sub-step, feed and draw, the trial and real passes | 2–5 |
| `tests/test_reentry_model_skin.py` (new) | the skin alone, against exact solutions and its own books | 2–5 |
| `reentry_model/thermal/__init__.py` | `StepResult.Q_interface`; the protocol's `interface`, `state`, `set_state` | 6 |
| `reentry_model/thermal/skfem_backend.py`, `fenicsx_backend.py` | the per-facet affine interface flux, no radiation under skins, the state round trip | 6 |
| `tests/test_reentry_model_thermal.py`, `tests/test_reentry_model_fenicsx.py` | backend tests for the interface; the two-backend skin test | 6, 7 |
| `reentry_model/body.py` | `MeltSettings.surface_model`/`skin`; skins in `MeltingBody`; the skin history columns | 7, 8 |
| `tests/test_reentry_model_melting_skin.py` (new) | the melting body with skins | 7 |
| `reentry_model/coupled.py`, `reentry_model/cli.py` | columns, frame fields, results, flags, run names | 8 |
| `tests/test_reentry_model_cli.py` | flags, names, a short skin run | 8 |

---

### Task 0: Stage 0 — the grid gate on the current model

**Files:**
- Read: the 2.0 mm, 0.0125 s runs of facts 91–95 (`prototype/work-2026-10-07-scheil-default/runs/`, both seeds) and
  the switched-step comparison script beside them (`.../harness/compare.py`)
- Create: `prototype/work-<date>-gridgate/` (a frozen measurement copy and its runs)

**Interfaces:**
- Consumes: the model as committed (facts 88–100). Its default `--dt-continuum` of 0.0125 s reproduces the fact-91
  harness's 0.0125 s run bit for bit (fact 99), so the existing 2.0 mm runs are the comparison and no harness is needed.
- Produces: one fact in the shared context — whether the converged droplet population moves with the surface mesh.

- [ ] **Step 1: Freeze a measurement copy**

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
G="$(pwd)/work-$(date +%F)-gridgate"; mkdir -p "$G"
rsync -a --exclude __pycache__ --exclude reentry_model_output work-2026-10-07-dt-continuum/code/ "$G/meas/"
cp work-2026-10-07-scheil-default/harness/compare.py "$G/"
```

Confirm the copy is the committed state before using it: its unit tier gives the baseline counts above.

- [ ] **Step 2: Launch the two 1.4 mm runs (seeds 12345 and 1), side by side**

The 100 mm flight of facts 91–95 with every default (Scheil, the rule on, `--dt-continuum` 0.0125 s) and
`--h-surface 1.4`:

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
for SEED in 12345 1; do
  OUT="$G/runs/hs14_seed$SEED"; mkdir -p "$OUT"
  ( cd "$G/meas" && echo "power at start: $(pmset -g batt | head -1)" > "$OUT.log" && \
    caffeinate -i "$PY" -m reentry_model run --diameter 100 --altitude 77.500133 --velocity 7.5 \
      --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --h-surface 1.4 \
      --seed $SEED --t-max 120 --outdir "$OUT" >> "$OUT.log" 2>&1; echo "exit $?" >> "$OUT.log"; \
    echo "power at end: $(pmset -g batt | head -1)" >> "$OUT.log" ) &
done
```

Expected: both logs end `exit 0`, and each run's JSON records the switch at about 49.5 s. The 2.0 mm default run took
2.3 h (fact 99); the 1.4 mm mesh has roughly 2.9 times the surface elements, so expect something like 7 to 10 h each
with both running at once (an estimate, unmeasured).

- [ ] **Step 3: Compare with the 2.0 mm series**

Run `compare.py` on the two 1.4 mm runs and on the 2.0 mm 0.0125 s runs (both seeds) of facts 91–95, over 49.5–120 s.
Tabulate, for each mesh: droplet count, the thick branch's share of the sprayed mass, median radius by number and by
mass, sprayed mass, mass at 120 s, the runoff measures and the run time. Read every change against the scatter between
the two seeds **on the same mesh**, and against fact 93's scatter at 2.0 mm (0.74 % in count, 0.31 points in thick
share, 0.16 % and 0.07 % in the medians, 0.14 % in sprayed mass, 0.53 % in mass at 120 s).

- [ ] **Step 4: Record the verdict and stop for Asha**

Append a section to `00-shared-context.md`, "Amendment of <date> — the melt-layer skin, stage 0: the grid gate (fact
<n>)", opening the way the 2026-09-27 section does: "Design: `docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md`,
approved 2026-10-07. Plan: sub-plan 18 (stages 0 and 1); stages 2 and 3 follow in sub-plan 19." Then one fact (the next
free number) stating which quantities moved beyond the scatter, by how much, and the verdict: **mesh-dependent** (the
skin is needed for correctness) or **mesh-independent** (the skin is a cost measure). Another agent may be appending to
the shared context; re-read the file immediately before the edit and append only. Report to Asha in her reporting style
and wait for her decision before Task 9's measurements. Tasks 1–8 build code only and may proceed while the runs go.

---

### Task 1: The working copy for stage 1

**Files:**
- Create: `prototype/work-<date>-skin/code/` (copy), `prototype/work-<date>-skin/unit_before.log`

- [ ] **Step 1: Copy the latest tested copy and verify it**

Use the tested copy of the latest committed Step 3 amendment, `prototype/work-2026-10-07-dt-continuum/code/` (facts
97–100); this plan's diff blocks are made against it. If a later amendment has been committed since, use its copy and
apply the blocks by hand where the context moved. Verify it first: re-apply every committed amendment's
diff blocks in date order to a fresh copy of `prototype/proto3/` and compare byte for byte (`diff -rq`, excluding
`__pycache__` and `reentry_model_output`). The amendment authors' block-extraction scripts live beside their copies;
reuse one, extended to find every section with a given date, as the 2026-10-05 (runoff flux) amendment did.

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
W="$(pwd)/work-$(date +%F)-skin"; mkdir -p "$W"
rsync -a --exclude __pycache__ --exclude reentry_model_output work-2026-10-07-dt-continuum/code/ "$W/code/"
```

- [ ] **Step 2: Record the unit tier and the FEniCSx tier before any change**

```bash
cd "$W/code" && "$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py | tail -3 | tee ../unit_before.log
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q | tail -2 | tee -a ../unit_before.log
```

Expected: {{UNIT_BEFORE}}; FEniCSx tier 12 passed.

---

### Task 2: `skin.py` — settings, storage, rezoning, the top node, depths and hand-over

**Files:**
- Create: `reentry_model/skin.py`
- Test: `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: `material.Material` (`rho`, `enthalpy`, `enthalpy_liquid`, `enthalpy_mixed`, `cp_mixed`, `k`,
  `temperature_from_enthalpy`, `temperature_from_enthalpy_mixed`, `T_feed`, `T_rigid`); `thermal.SIGMA_SB`.
- Produces: `SKIN_PRESETS` (`{"production": 10e-6, "dev": 20e-6}` m), `INTERFACE_TOLERANCE` (1e-3), `NEWTON_TOL`
  (1e-6 K), `NEWTON_MAX` (40), `LINE_SEARCH_MAX` (12), `DZ_FLOOR`; `SkinSettings(thickness=0.4e-3, cell=10e-6,
  substep=0.010)` with `.n_cells`; `SkinBook` (per-skin arrays `E_top, E_rad, E_base, to_film, to_deep, drawn,
  drawn_H, L`); `thomas(sub, diag, sup, rhs)`; `cumulative(m, Q, x)`; `SkinField(material, settings=None)` with `.n`,
  `.nc`, `.face_id`, `.area`, `.m`, `.H`, `.T` and `mass()`, `energy()`, `thickness()`, `state()`, `set_state(state)`,
  `add(face_id, area, T_profile, thickness) -> (mass, H)`, `keep(rows)`, `rows_of(face_ids)` (-1 for none),
  `top_energy(rows, L)`, `equilibrate_top(rows, L, E_top) -> T0`, `book_top(rows, dE, L)`, `refreeze_top(rows, taken,
  L_after)`, `_rezone(rows, start, D, hD, S=None, S_drawn=1.0) -> (m, H, S)`, `_set_cells(rows, m, H, L, E_liq_top)`,
  `remap(rows, L)`, `add_bottom(rows, mass, H, L)`, `depth_above(T_level) -> (depth, through)`,
  `hand_over(src, tgt_faces, weights, face_area) -> (receiving rows, mass handed)`. Tasks 3–5 add `conduct`,
  `feed_and_draw` and `advance` to the same class.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_skin.py` with this content (Tasks 3–5 append to it):

```python
{{TEST_SKIN_A}}
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: FAIL with `ImportError: cannot import name 'skin'`.

- [ ] **Step 3: Write the module's first half**

Create `reentry_model/skin.py`:

```python
{{SKIN_A}}
```

- [ ] **Step 4: Run the tests to make sure they pass**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: 7 passed.

---

### Task 3: `skin.py` — the implicit sub-step

**Files:**
- Modify: `reentry_model/skin.py` (append `SkinField.conduct` to the class)
- Test: `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Task 2's `SkinField`, `thomas`, the module constants.
- Produces: `SkinField.conduct(dt, q_top, T_amb, emissivity, T_b, L, T_surface=None, sens=None) -> (E_top, E_rad,
  E_base, sens)`, energies in J per skin; `sens = (S, dEb)`, `S` (n, nc) the derivative of the cell temperatures with
  respect to the base temperature and `dEb` (n,) that of the base energy, carried through the sub-steps of a trial pass.

Why the iteration is plain Newton in the temperature with a line search, not the 3D solver's enthalpy-consistent
update: in 10 µm cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step (a cell's
conductance to its neighbour 38 W/K against a capacity rate of 0.07 W/K, measured while planning), so the cells must move
together. Inverting each cell's own h(T) moved neighbouring cells apart by a tenth of a 7 K step wherever the Scheil
curve's effective heat capacity varies, and the iteration stalled at a 0.4 K change on 888 of 2000 skins. The line search
(halving the step per skin until the residual's norm falls) does across a sharp latent ramp what that update did: the
two-phase Neumann test below, with a ±2 K ramp, needs it.

- [ ] **Step 1: Append the failing tests**

```python
{{TEST_SKIN_B}}
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: the Thomas test passes (Task 2 defined `thomas`); the other three FAIL with
`AttributeError: 'SkinField' object has no attribute 'conduct'`.

- [ ] **Step 3: Append `conduct` to `SkinField`**

```python
{{SKIN_B}}
```

- [ ] **Step 4: Run the tests to make sure they pass**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: 11 passed.

---

### Task 4: `skin.py` — feed and draw

**Files:**
- Modify: `reentry_model/skin.py` (append `SkinField.feed_and_draw`)
- Test: `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Tasks 2–3 (`_rezone`, `_set_cells`).
- Produces: `SkinField.feed_and_draw(L, room, avail, T_b, S=None) -> (to_film, to_deep, drawn, drawn_H, L_new,
  S_new)`, per skin. `room` [kg] is ρ A δ_m of metal where Girin's closure applies and `inf` elsewhere; `avail` [kg]
  is what the element below can still give; `S` is the trial pass's sensitivity, rezoned with the cells (drawn metal,
  held at the base temperature, carries 1).

- [ ] **Step 1: Append the failing tests**

```python
{{TEST_SKIN_C}}
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: the three new tests FAIL with `AttributeError: 'SkinField' object has no attribute 'feed_and_draw'`.

- [ ] **Step 3: Append `feed_and_draw` to `SkinField`**

```python
{{SKIN_C}}
```

- [ ] **Step 4: Run the tests to make sure they pass**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: 14 passed (the steady-ablation test takes about 5 s).

---

### Task 5: `skin.py` — the macro step's trial and real passes, and the base law

**Files:**
- Modify: `reentry_model/skin.py` (append `SkinField.advance`)
- Test: `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Tasks 2–4.
- Produces: `SkinField.advance(dt, q_top, T_amb, emissivity, T_b0, L, room, avail, T_b1=None)`. With `T_b1=None`, the
  trial pass, on a copy, base held at `T_b0`, returning `(a, b)` per skin: the law q_b(T) = a + b T [W/m²] of heat into
  the 3D model through the skin's facet, its slope clipped at zero. Otherwise the real pass, base held at `T_b1` (the 3D
  model's end-of-step facet temperature, spec §6.2 as revised), returning a `SkinBook`.

What the law is good for, measured while planning on Scheil skins: exact at the trial pass's own base temperature
(relative error 1e-16); a kelvin away, good to between 9e-6 and 5e-5 where nothing feeds; while a skin feeds, good to
about 1 % per kelvin (0.93 % at 1 K and 3.4 % at 5 K on one skin; in a 2000-skin stress case up to 1.9 % per kelvin from
a mushy start and 12 % from a start at the liquidus), because the feed comes in whole cells. That error is what the
interface mismatch, its booking and the step's repeat exist for; Task 9 measures how often a flight repeats. A slope
above zero occurred on a skin whose base was held inside the top of the mushy range (+1.65e5 W/(m² K)); it is clipped,
since a law rising with the base temperature would cost the 3D matrix its diagonal dominance.

- [ ] **Step 1: Append the failing tests**

```python
{{TEST_SKIN_D}}
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: the three new tests FAIL with `AttributeError: 'SkinField' object has no attribute 'advance'`.

- [ ] **Step 3: Append `advance` to `SkinField`**

```python
{{SKIN_D}}
```

- [ ] **Step 4: Run the tests, then the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`, then the unit tier.
Expected: 17 passed (about 8 s); the unit tier's counts are the baseline's plus these 17.

---

### Task 6: The thermal backends — the interface flux, no radiation under skins, the state round trip

**Files:**
- Modify: `reentry_model/thermal/__init__.py`, `reentry_model/thermal/skfem_backend.py`,
  `reentry_model/thermal/fenicsx_backend.py`
- Test: `tests/test_reentry_model_thermal.py`, `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Produces: `StepResult.Q_interface: float = 0.0` [W]; `step(dt, q_conv, T_amb, dirichlet=None, nodal_load=None,
  interface=None)` with `interface = (mask, a, b)` per facet — on masked facets the heat flux a + b T_f [W/m²] enters the
  body implicitly (T_f the facet-mean temperature, its Jacobian on the facet's node pairs like the radiation's) and the
  facet does not radiate; `None` or an all-False mask is the old step bit for bit (the scalar emissivity is then used
  exactly as before); `radiated_power` leaves out the masked facets of the last step (`_radiating`, reset whenever the
  surface is rebuilt); `state()` and `set_state(state)` restore the temperatures, the predictor's previous field and
  the solve counter, so that a step can be repeated.

- [ ] **Step 1: Append the failing tests**

To `tests/test_reentry_model_thermal.py`:

```python
{{TEST_THERMAL_T6}}
```

To `tests/test_reentry_model_fenicsx.py`:

```python
{{TEST_FENICSX_T6}}
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q -k "masked or round_trip"`
Expected: FAIL with `TypeError: ... step() got an unexpected keyword argument 'interface'` and
`AttributeError: ... has no attribute 'state'`.

- [ ] **Step 3: Apply the three diff blocks**

```diff
{{DIFF_THERMAL_INIT}}
```

```diff
{{DIFF_SKFEM}}
```

```diff
{{DIFF_FENICSX}}
```

- [ ] **Step 4: Run the tests, the FEniCSx tier and the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`; then
`FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`;
then the unit tier.
Expected: the thermal file 18 passed (16 before); the FEniCSx tier 13 passed (about 95 s); the unit tier's counts the
previous task's plus 2.

---

### Task 7: The skins in the melting body

**Files:**
- Modify: `reentry_model/body.py`
- Create: `tests/test_reentry_model_melting_skin.py`
- Modify: `tests/test_reentry_model_fenicsx.py` (append one test)

**Interfaces:**
- Consumes: Tasks 2–6 (`skin.SkinField`, `skin.SkinSettings`, `skin.INTERFACE_TOLERANCE`, the backends' `interface`,
  `state`, `set_state`, `StepResult.Q_interface`).
- Produces: `body.SURFACE_MODEL_NAMES`; `MeltSettings.surface_model` (`"elements"` default, or `"skin"`) and
  `MeltSettings.skin` (`skin.SkinSettings` or `None` for its defaults); on `MeltingBody`: `skins` (`None` with
  `elements` or `--removal instant`), `skin_of_patch`, `patch_of_skin`, `skins_created`, `skin_drawn_mass`,
  `skin_handover_mass`, `skin_short_events`, `interface_repeats`, `last_mismatch`, `skin_wall`, `last_skin_liquid`,
  `last_skin_nonrigid`; the methods `_advance_with_skins`, `_advance_elements` (the old `advance`, unchanged),
  `_create_skins(T) -> int`, `_skin_available(patches)`, `_compose_depth(level, depth3d)`, `_defer_to_facets(energy)`,
  `_freeze_into_skins(solid)`, `_death_to_skins(dead, rest, h_e)`, `_hand_over_targets(old_surface, idx)`,
  `_hand_over_skins(old_surface)`, `_refresh_skin_index()`; and overrides of `surface_temperature` and
  `radiated_power` that read the skins' tops.

What the diff does, so a reviewer can check it against the spec: (a) `advance` dispatches; with skins it creates them,
evaluates the surface flow once before the passes (the feed's reach needs δ_m; the melt step reuses that evaluation),
runs the trial pass, the 3D step with the law on the skins' facets and no heating or radiation there, the real pass,
books the mismatch on the facets' nodes, and repeats once, re-linearised about the new facet temperature, if any skin's
mismatch exceeds `INTERFACE_TOLERANCE` of its energy throughput. (b) The melt step applies the skins' draws to their
elements (the element keeps what the drawn metal left, as a deferred load), stops the feed of elements that own a skin
patch, freezes liquid on skin patches into the skin's top cell, turns the molten cascade off, and sends a dying
element's remainder into the skins above it. (c) Liquid on a skin patch shares the skin's top cell temperature: its
enthalpy is h_l(T₀), energy booked on such a patch goes to the top node at once, and it is no longer on the 3D nodes.
(d) Because a booking moves a skin's top within the melt step, the step's cached liquid enthalpy `h_p` is refreshed
before the deep runoff, before the film and spray, after the runoff's booking inside `_film_and_spray`, and before each
death pass's `before` sum — each refresh only with skins, so the element path's arithmetic is untouched. (e) Skins are
handed on at deaths by the liquid's own rule (the dying element's exposed faces by projected area, the nearest patch for
an orphan), layer by layer; empty skins are dropped. (f) At consumption the skins' metal leaves with the liquid.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_melting_skin.py`:

```python
{{TEST_MELTING_SKIN}}
```

Append to `tests/test_reentry_model_fenicsx.py`:

```python
{{TEST_FENICSX_T7}}
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_melting_skin.py -q`
Expected: FAIL with `TypeError: MeltSettings.__init__() got an unexpected keyword argument 'surface_model'`.

- [ ] **Step 3: Apply the diff block**

```diff
{{DIFF_BODY_T7}}
```

- [ ] **Step 4: Run the tests, the FEniCSx tier and the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_melting_skin.py -q --durations=10`; then the FEniCSx tier; then the unit
tier.
Expected: 8 passed in about 60 s (the coupled flight steps about 22 s, the drawing test about 19 s); FEniCSx 14 passed;
the unit tier's counts the previous task's plus 8, every earlier test unchanged — every one of them runs
`surface_model="elements"`.

---

### Task 8: CLI flags, run names, history columns, results and frame fields

**Files:**
- Modify: `reentry_model/cli.py`, `reentry_model/coupled.py`, `reentry_model/body.py`
- Test: `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: Task 7.
- Produces: CLI `--surface-model {elements,skin}` (default `elements` in stage 1), `--skin-thickness` [mm] (0.4),
  `--skin-cell` [µm] (default: the preset's), `--skin-substep` [ms] (10), `--skin-preset {production,dev}`;
  `cli.skin_settings(args)`; `cli.skin_name_suffix(thickness_mm, cell_um, substep_ms, preset) -> str`;
  `model_run_name(..., surface_model="elements", skin_suffix="")`; `coupled.SKIN_COLUMNS` (in `MELT_COLUMNS`) filled by
  `MeltingBody._skin_stats()` through `melt_stats()`; result fields `skins_created`, `skin_drawn_mass_kg`,
  `skin_handover_mass_kg`, `interface_repeats`, `skin_wall_s`; frame fields `skin_T_top`, `skin_liquid_depth`,
  `skin_nonrigid_depth`, `skin_thickness` (nan on a patch without a skin); the settings record `surface_model` and
  `skin`.

- [ ] **Step 1: Append the failing tests**

To `tests/test_reentry_model_cli.py`:

```python
{{TEST_CLI_T8}}
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q -k skin`
Expected: the run-name test FAILS (`model_run_name() got an unexpected keyword argument 'surface_model'`) and so does
the short run (argparse rejects `--surface-model`); the five bad-argument cases already pass, because argparse rejects the
unknown flags with exit 2 — once the flags exist they pin the validation instead.

- [ ] **Step 3: Apply the three diff blocks**

```diff
{{DIFF_CLI}}
```

```diff
{{DIFF_COUPLED}}
```

```diff
{{DIFF_BODY_T8}}
```

- [ ] **Step 4: Run the tests and the unit tier; confirm the element model's bits**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q -k skin`, then the unit tier.
Expected: {{CLI_EXPECTED}}; the unit tier {{UNIT_AFTER}}.

Then confirm bit-identity with the copy before the change (`work-2026-10-07-dt-continuum/code`): the 50 mm default
melting flight to the ground (it never switches step), and the 100 mm default flight to 55 s (it switches to 0.0125 s at
49.5 s), each run from both copies with the same arguments and the default seed, agree in every history column they
share, every source-table row (`particles.npz`) and every result field but the run time. The comparison script used
while planning is `prototype/work-2026-10-07-skinplan/bitid/compare.py`; record its output.

---

### Task 9: Stage 1 measurements, and the amendment written into the sub-plans

**Files:**
- Create: `prototype/work-<date>-skin/meas/` (frozen copy), `.../runs/`, `.../analysis/`
- Modify (documents only): `00-shared-context.md`, `03-thermal-core.md`, `09-melting-body.md`, `10-coupled-loop.md`,
  `13-cli-wiring.md`, `14-verification-runs.md`, `15-documentation.md`, and this sub-plan's status line.

- [ ] **Step 1: Freeze the measurement copy and record its manifest**

`rsync` the tested `code/` to `meas/`; write `meas_manifest.md5` (`find meas -name '*.py' | sort | xargs md5`); never
edit `meas/` while a run uses it.

- [ ] **Step 2: The box against the fine box (spec §12.2)**

Write `analysis/box_ablation.py`: an 8 × 2 × 2 mm box (`mesh.box_mesh(8e-3, 2e-3, 2e-3, h)`) heated by a uniform
2 MW/m² on its x = 0 face only (the face its analytic tests heat), no heating on the other faces and radiation off
everywhere (`MeltingBody(..., emissivity=0.0)`), so that the test is one-dimensional; Scheil material, starting at 700 K,
for 2 s. Liquid is removed as soon as it reaches the patch's liquid (the script zeroes `m_f` and `m_d` after each step
and adds their mass and enthalpy to `removed_mass` and `removed_enthalpy`, so the books still close). Compare (i) the
coarse box, `h = 2 mm`, with skins at the dev and production presets and macro steps of 0.05, 0.1, 0.25 and 0.5 s, with
(ii) the fine box, `h = 0.1 mm`, element model, at 0.005 s. Report the recession (mass removed per unit area) and the
top temperature against time, the mismatch and repeats per step, and the energy balances. The skin must recover the
fine box's recession within 5 % at every macro step, with no oscillation in the top temperature.

- [ ] **Step 3: The flights**

All on the 100 mm physics flight to 120 s (Scheil default, seed 12345, `--surface-model skin`). The macro step after
the continuum switch is set explicitly on every run: `--dt-continuum off` keeps 0.5 s throughout, `--dt-continuum 0.25`
and `0.125` shorten it after the switch.
(a) macro step after the switch 0.5 s (`off`), 0.25 s and 0.125 s, production preset, sub-step 10 ms;
(b) sub-step 20, 10 and 5 ms, and cell 20, 10 and 5 µm, with `--dt-continuum off`;
(c) skin thickness 0.4 against 0.8 mm, with `--dt-continuum off`;
(d) `--h-surface` 2.0 against 1.4 mm, with `--dt-continuum off`;
(e) seed 1 beside the default, with `--dt-continuum off`, for the scatter;
(f) skins at the default 0.0125 s step after the switch, against the element model's default runs of stage 0 (the
converged reference of facts 92–94), if (a)'s measured cost makes it affordable — it is the comparison spec §12.3 calls
"against the converged reference".
Plus the 50 mm whole flight with skins (it never switches). Read every comparison against the scatter at its own step. Report the
droplet population and branch split (expected still step-dependent in stage 1: the film and spray still run once per
macro step, which is stage 2's to change), the masses, the interface mismatch and repeat rate, the skin thickness
statistics and short events, the books, and the cost per macro step and peak memory at each preset. Compare with stage
0's verdict and facts 92–94. The cost measured while planning, for scale: 2000 skins take 1.9–2.7 s for a trial pass and
1.6–2.5 s for a real pass at 10 µm, 1.0–1.7 s and 0.9–1.5 s at 20 µm (one macro step of 0.5 s, single runs).

- [ ] **Step 4: Bit-identity of the element model**

Re-run the copy's earlier 50 mm and 100 mm (to 120 s) default flights with `--surface-model elements` from `meas/` and
show they reproduce the pre-change runs in every shared history column, every source row and every result field but the
run time (Task 8's check, now from the frozen copy).

- [ ] **Step 5: Write the amendment**

In each affected sub-plan, a section "Amendment of <date> — the melt-layer skin, stage 1" with the tested changes as
diff blocks against the copy before the change: 03 the backends' interface (Task 6); 09 `skin.py` in full as a new file
and the body's changes (Tasks 2–5, 7, 8's body part); 10 the columns, results and frame fields; 13 the flags and names;
14 the runs and their results; 15 README, assumptions and the spec's measured values. Add the measurements as new facts
in `00-shared-context.md` (next free numbers) and verify that every committed amendment's blocks plus these reproduce
`meas/` byte for byte from `prototype/proto3/`. Update this sub-plan's status line to "stage 1 built and measured;
stages 2–3 in sub-plan 19". Report to Asha in her reporting style, with the decisions stage 1 raises: the presets, the
interface tolerance (and whether the repeat rate makes it too tight), and whether to proceed to stage 2.

---

## Self-review

**Spec coverage.** §5.1–5.3 (the skin, creation, ownership): Tasks 2 and 7. §5.4 (the sub-step, feed, draw, rezoning,
bounds): Tasks 3, 4 and 7. §5.5 (resolved depths): Tasks 2 (`depth_above`) and 7 (`_compose_depth`). §6 (two clocks,
the interface, the repeat, the books): Tasks 5–7. §7 (the film and spray on the sub-steps): **stage 2, sub-plan 19** —
in stage 1 the film and spray run once per macro step, fed by the skins (spec §12.1, stage 1). §8 (the slurry cascade):
**stage 3, sub-plan 19**. §9 (deaths, hand-over, end of flight): Tasks 2 and 7; remeshing waits for sub-plan 17. §10
(retired and kept): Task 7 — feed off for skin owners, cascade off, remainder to the skin, film on the skin's top,
freeze-back into the skin, resolved depths in the branch test; interior feed, deep runoff, the rigid-substrate rule,
`PHI_DEATH` and the once-per-step surface flow kept. §11 (settings, outputs, code): Task 8, with the CLI default left at
`elements` until stage 2's acceptance. §12.1 stage 0: Task 0; stage 1: Tasks 1–9. §12.2 tests: Tasks 2–7 and Task 9
Step 2. §12.3 measurements: Task 9 Step 3 (acceptance proper is stage 2's).

**Corrections made while planning and checking the code, all recorded in the spec (2026-10-07):** the real pass holds
the base at the 3D model's end-of-step facet temperature, not a ramp (§6.2); the conjugate-depth reach is measured from
the skin's top in metal (§5.4); the skin iterates in the temperature with a line search rather than the 3D solver's
enthalpy-consistent update (§5.4); the base law's slope is clipped at zero, and the interface tolerance is measured
against each patch's energy throughput rather than its base energy alone (§6.2). The surface flow is evaluated before
the passes, at the start-of-step wall temperatures, as the spec's step diagram already has it; the cached liquid
enthalpy's refreshes (Task 7 (d)) are an implementation detail with no spec counterpart.

**Placeholder scan.** None: every code step carries its code or its diff, generated from the tested files.

**Type consistency.** `feed_and_draw` returns six arrays in the order `advance` unpacks; `advance` returns `(a, b)` in
the trial pass and a `SkinBook` otherwise, as `_advance_with_skins` uses them; `hand_over` returns `(rows, mass)`, of
which `_hand_over_skins` books the mass; `_skin_available` serves creation and the step alike; `SKIN_PRESETS` keys are
the CLI's `--skin-preset` choices.

**Review Focus.** Each of the five lines has its test in the owning task (Tasks 2, 4 and 7).
