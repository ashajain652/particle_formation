# Sub-plan 20 — skins on every corner, heat only (stage 2 of the layered-skin design) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first (required reading for every Step 3 sub-plan), then the spec.

**Goal:** Every corner of the moving 3D mesh's surface carries a layered skin from time zero that holds its outer metal
in fine cells; heat is conducted through the depth and sideways; liquid leaves as it forms; the base moves the 3D mesh by
the metal crossing it; the moving mesh keeps its energy where it is; and the whole is verified against exact solutions
(level 1), against the one-dimensional reference on a plate (level 2), and measured for cost.

**Architecture:** Four new pieces on top of sub-plan 19's moving mesh. `reentry_model/columns.py` holds the corner
columns' geometry (dual areas, normals, the edge graph with cotangent weights). `reentry_model/skin.py` is rewritten as
the layered skins: per-column cells in padded arrays, split and merge, the minimum mass, the implicit through-depth heat
solve with removal and draw, the explicit sideways conduction and the hand-down; stage 1's one-dimensional skins move to
`skin1d.py` until sub-plan 22 deletes them. `reentry_model/skinbody.py` is the body: a `ThermalBody` whose 3D model lies
under the skins, coupled by a base law at the corner nodes through trial and real passes, whose base moves the mesh.
The thermal backends gain the nodal base law; `movemesh.py` gains own-sweep booking, the volume-neutral smoothing and
an energy-local rebuild transfer. The command line selects the skins with `--surface-model skin`.

**Tech Stack:** Python 3.12 (`drama_env`), NumPy, SciPy, gmsh 4.15.2, scikit-fem backend, FEniCSx 0.11 backend
(`fenicsx_env`), pytest.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md`, amended 2026-10-12: §5.1–5.3 (the
corners, the two surfaces, the moving mesh), §6.1–6.3 and §6.6–6.7 (what a skin holds, cells, heat, the base, depths),
§7 (time-stepping and the interface), §10 (settings, outputs, code), §11 (level 1's heat checks, level 2 with flow off,
build stage 2). Sub-plans 21 (flow) and 22 (spraying and the flight) build on this one.

**How this plan's code was checked.** Every code and diff block below was written and tested in a planning copy,
`prototype/work-2026-10-12-skin-heat-plan/code/`, made from `prototype/work-2026-10-12-moving-mesh/code/` at 70ecbdd
(sub-plan 19 built, sub-plans 11 and 12's reference data and loader). This document is generated from that copy by
`prototype/work-2026-10-12-skin-heat-plan/make_plan20.py`, which replays the tasks below in order on a copy of the base —
every diff through `patch`, every file and appended section as written here, the renames and the reference extraction —
and checks that the result is the tested copy byte for byte. The unit tier goes from 341 to 377 passed (1 skipped, no
failures); the FEniCSx tier from 17 to 19. Measurements quoted were each made once, on the case named, on 2026-10-12.

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

- The slurry viscosity law comes from **`reentry_model/shared/rheology.py`** (stage S's **S1**), `material` from
  `reentry_model.shared`. `columns.py` and the rewritten `skin.py` are skyscrapers' own and are **never shared** with
  the level set.

**How it applies to this sub-plan as written** (2026-10-12). This sub-plan is built in the prototype, where Step 3 code
stays until the skyscrapers architecture is verified (Asha, 2026-10-12), and the prototype has no `reentry_model/shared/`
package yet; so the code below imports `material` from its current path, and stage 2 uses no slurry viscosity (that is
sub-plan 21's). When the verified code is ported to the root package, the imports follow this amendment;
`columns.py`, `skin.py` and `skinbody.py` are skyscrapers' own and never shared. The amendment above is the other
session's text, kept verbatim from this sub-plan's outline.

## Decisions and rulings

Asha's decisions while this plan was written (2026-10-12; the spec records them in §16):

- **Each corner is booked its own sweep** (spec §5.2 amended): a triangle's swept prism is shared among its corners by
  what each corner's own motion swept, A_t (u_k · n_t)/3, not a third to each. Found on the level-2 box: crediting
  non-moving side-wall corners with the receding front's metal made them hand it down, bulge out by 0.33 mm and fold the
  surface within 1.3 s. Task 0.
- **Level 2 is a plate with tapered heating**, not the q6p box, whose 90-degree edges fold the fixed surface
  triangulation (1.4 s on 2 mm elements, 0.9 s on 1 mm); local surface remeshing (spec §13) is the back-up. Task 9.
- **A task measures what triggers the interface repeat.** Task 10.

Rulings made while writing, with their reasons (each is stated again where its code is):

1. **The rebuild carries enthalpy by L2 projection** (spec §5.3, under its own fallback clause): interpolating at the new
   nodes moved 2 768 J — 4.9 % of the heat absorbed — from the plate's cold bulk to its melting front in two rebuilds,
   and the centre melted 9.5 % ahead of the exact solution; the projection leaves 14 J. Task 0.
2. **Stage 1's skins are renamed `skin1d.py`**, selected by `MeltSettings(surface_model="skin1d")` from Python only, so
   that the layered skins take the name `skin.py` (spec §10.3) and stage 1's tests keep running until sub-plan 22
   deletes it. Tasks 2 and 8.
3. **The default surface model stays `elements`** until sub-plan 22: stage 2 has no flow or spraying, so a skin run with
   melting is a verification device and needs `--removal instant`. Task 8.
4. **A draw takes rigid metal**: at the base temperature but never above T_rigid. A 4 mm 3D element under 6 MW/m² reached
   930 K in one 0.5 s step, and drawing at its temperature fed liquid that the next sub-step removed, without limit.
   Task 4.
5. **The skins look their material up in 0.01 K tables** (`PropertyTable`), and every pass starts from each cell's stored
   enthalpy, so the books are exact whatever the table's error (at most about 60 J/kg near a Scheil node). The property
   calls were half of a pass's time. Task 2.
6. **Only the columns that need it are re-split or merged**, and the Newton tangent is the table's own slope: together
   with the tables, a macro step on the 4 mm test sphere went from 11 s to 8 s. Tasks 2 and 4.
7. **The surface gradient and curvature operators and the cells' velocity are left to sub-plan 21**, the first stage to
   use them; nothing here could test them.
8. **`SkinBody` is a module of its own**, `skinbody.py`: `body.py` holds the element model, which must stay bit for bit.
9. **The interface keeps the spec's per-column repeat as the default**, with the other rules selectable for Task 10's
   measurement. A predictor for the linearisation temperature was tried and dropped: it ran away on 4 mm elements.
10. **A skin run's frames carry per-corner data**; `--animate` and `--stills`, which read the element model's frames, are
    refused with skins. Task 8.
11. **Skin runs mesh without prism layers and with the dense band**, which the moving mesh's rebuild reproduces. Task 8.
12. **Level 3's Step 2 runs against SESAM are the build's**, not this plan's planning: at the measured cost each takes
    hours (Task 10).

## Global Constraints

- Work in a fresh work copy, `prototype/work-2026-10-12-skin-heat/code/`, a git clone of
  `prototype/work-2026-10-12-moving-mesh/code/` checked out at 70ecbdd (Task 0, Step 0). It is a git repository of its
  own, ignored by the outer one; commit there after every task. Never edit `prototype/proto3/`, the package
  `reentry_model/` at the repo root, the specs or the monolithic plan. Step 3 code stays in the prototype until the
  skyscrapers architecture is verified (Asha, 2026-10-12).
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the work copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Baseline at 70ecbdd: 341 passed, 1 skipped. Record the counts after every task; no failure is acceptable.
- FEniCSx tier: `FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`. Baseline 17 passed.
- `--surface-model elements` (the default) is unchanged: a melting element-model run from the work copy has the same
  history as the base's, bit for bit (Task 8 Step 6).
- Books: in every skin test and run, the mass is closed to 1e-12 and the energy balance residual is below 1e-8.
- Exit codes: a move no rebuild can absorb raises `RuntimeError` (exit code 1 at the CLI); bad skin settings exit 2.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).

## Review Focus

1. **A coarse 3D base overshooting T_rigid or T_feed within a macro step** — the draw is capped at T_rigid, and the
   interface's repeat keeps the coupling stable (without it, 3 MW/m² on 4 mm elements runs away). Pinned by
   `test_a_draw_from_a_base_above_t_rigid_takes_rigid_metal` (Task 4) and
   `test_every_repeat_rule_keeps_the_books_and_records_what_triggered_it` (Task 7).
2. **Neighbouring corners receding at different speeds** — each is booked its own sweep. Pinned by
   `test_a_corner_that_does_not_move_is_credited_with_nothing_its_neighbours_swept` (Task 0).
3. **A rebuild under a steep field** — the projection keeps the energy where it was and a linear field exact. Pinned by
   `test_the_rebuild_offset_stays_in_the_deforming_zone_and_is_summed` and
   `test_a_rebuild_of_a_linear_field_loses_nothing_to_interpolation` (Task 0).
4. **A sub-step beyond the explicit sideways limit** (a finer surface, a smaller column) — sub-cycled and bounded.
   Pinned by `test_a_sub_step_beyond_the_explicit_limit_is_sub_cycled_and_stays_bounded` (Task 5).
5. **A column drained below its minimum, or above it with rigid metal at its bottom** — it draws back, or hands the
   excess down and the base moves out. Pinned by `test_a_skin_at_its_minimum_above_it_and_drained_returns_to_it_by_drawing`
   (Task 4) and `test_handing_down_moves_the_base_out_and_the_3d_model_takes_the_metal` (Task 7).

---

### Task 0: The moving mesh's books

**Files:**
- Modify: `reentry_model/movemesh.py` (`NEUTRAL_PASSES`, `MovingMesh.neutral_smoothing`, `swept_volumes`, the rebuild's
  transfer, new `project_enthalpy`, `QUADRATURE`, `PROJECT_CHUNK`)
- Modify: `analysis/movemesh_recession.py` (`--neutral`)
- Test: `tests/test_reentry_model_movemesh.py` (three tests changed for the projection, three appended)

**Interfaces:**
- Produces: `swept_volumes(old, new, faces, faces_corner, n_corners)` — same signature, each corner booked its own
  sweep (spec §5.2 as amended); `MovingMesh.neutral_smoothing(normal=None, area=None, passes=NEUTRAL_PASSES) -> (nc, 3)`
  (the tangential smoothing plus the normal corrections that cancel what it sweeps per corner);
  `project_enthalpy(old_points, old_tets, h_old, points, tets) -> (n,)`, the rebuild's transfer.

Why: on the level-2 box the third-of-a-prism booking credited the side walls with metal the front swept (Decisions).
With own sweeps, the plain smoothing books at most 0.46 µm at a corner of the 100 mm sphere (1.03 µm before), and three
normal corrections bring that to 3·10⁻¹¹ m (0.37 µm before), each pass removing about 95 % of what is left. The rebuild:
interpolating at the new nodes and offsetting the total uniformly over the zone moved heat from the cold bulk to the
melting front (2 768 J in two rebuilds of the plate, 4.9 % of its absorbed heat); the L2 projection of the old field —
the new mesh's consistent mass matrix against the old field integrated over the new elements by the 4-point rule — keeps
each region's energy to the rule's accuracy and reproduces a linear field exactly (a lumped projection, also tried,
moved a linear field's surface values by up to 7 K). On the tests' kinked field the residual goes from −4.1 % of the
energy to −1.8 J of 63 219 J.

- [ ] **Step 0: The work copy**

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
git clone work-2026-10-12-moving-mesh/code work-2026-10-12-skin-heat/code
git -C work-2026-10-12-skin-heat/code checkout -b skin-heat 70ecbdd
```

  Append to the outer repository's `.gitignore`, under the line for `prototype/work-2026-10-12-moving-mesh/code/`:
  `prototype/work-2026-10-12-skin-heat/code/`. Record the baseline: unit tier 341 passed, 1 skipped; FEniCSx 17 passed.

- [ ] **Step 1: Apply the test changes** (the projection's expectations in two rebuild tests, the volume-neutral
  smoothing and own-sweep tests appended):

{{DIFF_TEST_MOVEMESH}}

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`
Expected: the three new tests fail (`AttributeError: 'MovingMesh' object has no attribute 'neutral_smoothing'`, and
the own-sweep test's neighbours are credited two thirds of the moving corner's prism), and the two rebuild tests fail on
their residuals (the interpolation's −4.1 %).

- [ ] **Step 3: Apply the change to `reentry_model/movemesh.py`**

{{DIFF_MOVEMESH}}

- [ ] **Step 4: Run them to see them pass.** Expected: 34 passed.

- [ ] **Step 5: The recession measurement's switch** (`analysis/movemesh_recession.py`):

{{DIFF_RECESSION}}

  Run `caffeinate -i "$PY" analysis/movemesh_recession.py --neutral --out ../measure/recession_neutral.json` (about
  75 s). Measured while planning, with the third-of-a-prism booking: the run fails at 29.5 mm of nose recession against
  30.0 mm with the plain smoothing, the worst surface triangle quality 0.711 in both, the rebuilds' residuals 0.07–0.26 %
  of the total enthalpy in both. Record the result with own sweeps and the projection beside those.

- [ ] **Step 6: Unit tier and commit.** Expected 344 passed, 1 skipped.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py analysis/movemesh_recession.py
git commit -m "Moving mesh: each corner booked its own sweep, the smoothing made volume-neutral, and the rebuild's enthalpy carried by L2 projection"
```

### Task 1: The corner columns

**Files:**
- Create: `reentry_model/columns.py`
- Test: `tests/test_reentry_model_columns.py`

**Interfaces:**
- Consumes: `movemesh.corner_frame`.
- Produces: `Columns(points (nc, 3), faces (nt, 3) corner numbering, outward)` with `.n`, `.faces`, `.edges (ne, 2)`,
  `.degree`, and from `update(points)`: `.points`, `.normal (nc, 3)`, `.area (nc,)` (barycentric dual areas),
  `.centroid`, `.edge_length`, `.cot_weight` (½(cot a + cot b) per edge), `.conduct_weight` (the same clipped at zero),
  `.face_length` (the dual cell's two segments across each edge); `.top_points(thickness)`; `.laplacian(f, weights=None)`.
  `CornerSurface(points, faces, normals, areas)` with `from_points(points, faces)`, `.n_patches`, `.area`,
  `.angles_to(v)`, `.patch_toward(d)`, `.projected_area(v)` — what the coupled run reads as the body's surface, one
  patch per corner.

A negative cotangent weight (an edge whose opposite angles are obtuse enough) is clipped to zero for the explicit step's
positivity: 3 of 6 756 edges on the 4 mm test sphere, 4 of 28 245 on the default 2 mm sphere.

- [ ] **Step 1: Write the tests**

{{TEST_COLUMNS}}

- [ ] **Step 2: Run them to see them fail.** Run: `"$PY" -m pytest tests/test_reentry_model_columns.py -q`. Expected:
  `ImportError: cannot import name 'columns' from 'reentry_model'`.

- [ ] **Step 3: Create `reentry_model/columns.py`**

{{COLUMNS}}

- [ ] **Step 4: Run them to see them pass.** Expected: 5 passed.

- [ ] **Step 5: Unit tier and commit.** Expected 349 passed, 1 skipped.

```bash
git add reentry_model/columns.py tests/test_reentry_model_columns.py
git commit -m "Corner columns: dual cells, normals and the edge graph with cotangent conduction weights"
```

### Task 2: Stage 1 renamed; the layered skin's cells

**Files:**
- Rename: `reentry_model/skin.py` → `reentry_model/skin1d.py`, `tests/test_reentry_model_skin.py` →
  `tests/test_reentry_model_skin1d.py` (each with a docstring change)
- Modify: `reentry_model/body.py`, `reentry_model/cli.py`, `tests/test_reentry_model_melting_skin.py`,
  `tests/test_reentry_model_fenicsx.py` (they import `skin1d`)
- Create: `reentry_model/skin.py` (its header, constants, `SkinSettings`, `thomas`, `PropertyTable`, `Skins`)
- Create: `tests/test_reentry_model_skin.py` (its header, helpers and the Task 2 tests)

**Interfaces:**
- Produces: `SKIN_PRESETS`, `SPLIT_ABOVE = 1.5`, `MERGE_BELOW = 0.5`, `REGRID_PASSES`, `DELTA_M_MARGIN = 50e-6`,
  `INTERFACE_TOLERANCE = 1e-3`, `JUMP_TOLERANCE = 5.0`, `NEWTON_TOL`, `NEWTON_MAX`, `LINE_SEARCH_MAX`, `DZ_FLOOR`,
  `SIDEWAYS_CFL`, `TABLE_DT`, `TABLE_RANGE`; `SkinSettings(min_thickness=0.4e-3, cell=10e-6, growth=1.10,
  substep=0.010)`; `thomas(sub, diag, sup, rhs)`; `PropertyTable(material)` with `enthalpy`, `cp_eff`, `k`,
  `liquid_fraction`, `conduction(T) -> (h, dh/dT, k)`; `Skins(material, area, settings)` with `.m`, `.H`, `.T` (n, K;
  cell 0 the top), `.count`, `.area`, `.S` (the trial pass's dT/dT_b, or None), `.props`, `.valid`, `mass()`,
  `energy()`, `column_mass()`, `state()`, `set_state()`, `density(T)`, `dz()`, `thickness()`, `targets(dz=None)`,
  `fill(mass, T)`, `regrid()`, `take_top(k) -> (mass, H)`, `add_bottom(dm, h) -> H`, `take_bottom(dm) -> H`,
  `exchange_bottom(dm, dH)`, `depth_above(T_level)`.

Cells hold mass and enthalpy; temperature, liquid fraction and density follow from the material; a cell's thickness is
its mass over the column's area and its own density. A cell's target is the fine cell within the fine zone (the
minimum thickness below the top) and grows by (growth − 1) per unit depth below it, the continuous form of each cell
being `growth` times the one above. `regrid` merges a run of small cells into the cell that ends it (the bottom run
into the cell above) and halves big cells, only in the columns that need it, until none qualifies; every operation is
exact in mass and enthalpy, and carries `S` when the trial pass tracks it. Velocity is not stored: sub-plan 21 adds it
with its momentum books.

- [ ] **Step 1: Rename stage 1's skin**

```bash
git mv reentry_model/skin.py reentry_model/skin1d.py
git mv tests/test_reentry_model_skin.py tests/test_reentry_model_skin1d.py
```

  and apply:

{{DIFF_SKIN1D}}

{{DIFF_TEST_SKIN1D}}

{{DIFF_BODY_A}}

{{DIFF_CLI_A}}

{{DIFF_MELTSKIN_A}}

{{DIFF_FENICSX_A}}

  Run the unit tier: expected 349 passed, 1 skipped (stage 1's 21 tests now in `test_reentry_model_skin1d.py`).

- [ ] **Step 2: Write the new skin's tests** (`tests/test_reentry_model_skin.py`, which grows by a section per task):

{{TEST_SKIN_T2}}

- [ ] **Step 3: Run them to see them fail.** Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `ImportError: cannot import name 'skin' from 'reentry_model'`.

- [ ] **Step 4: Create `reentry_model/skin.py`**

{{SKIN_T2}}

- [ ] **Step 5: Run them to see them pass.** Expected: 7 passed.

- [ ] **Step 6: Unit tier and commit.** Expected 356 passed, 1 skipped.

```bash
git add -A reentry_model/skin.py reentry_model/skin1d.py reentry_model/body.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skin1d.py tests/test_reentry_model_melting_skin.py tests/test_reentry_model_fenicsx.py
git commit -m "Layered skins: stage 1 moves to skin1d, and the new skin's cells split and merge exactly"
```

### Task 3: The minimum mass

**Files:**
- Modify: `reentry_model/skin.py` (append), `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Produces: `minimum_mass(material, area, settings, delta_m=None) -> (n,)` kg: the larger of `min_thickness` of solid
  over the column's area and the liquid filling Girin's δ_m plus `DELTA_M_MARGIN` (spec §6.1); δ_m None or nan leaves the
  solid term alone (stage 2 evaluates no surface flow, so it passes None).

- [ ] **Step 1: Append the test**

{{TEST_SKIN_T3}}

- [ ] **Step 2: Run it to see it fail.** Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'minimum_mass'`.

- [ ] **Step 3: Append to `reentry_model/skin.py`**

{{SKIN_T3}}

- [ ] **Step 4: Run it to see it pass.** Expected: 8 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 357 passed, 1 skipped.

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Layered skins: the minimum mass, 400 um of solid or Girin's layer and 50 um as liquid"
```

### Task 4: Heat through the depth, with liquid removed as it forms

**Files:**
- Modify: `reentry_model/skin.py` (append), `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Produces: `SkinBook(E_top, E_rad, E_base, E_side, removed, removed_H, drawn, drawn_H)` per column;
  `_conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, Q_side=None)` — one backward-Euler sub-step (T_b None: an adiabatic
  base, for tests); `advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True)` — the
  macro step's sub-steps: with `T_b1` None the trial pass, on a copy, returning the base law `(a, b)` per column
  (E_base/dt = a + b T_b, b clipped at zero); otherwise the real pass, returning a `SkinBook`. `side(sk, h)` is the
  sideways heat per cell (Task 5).

Each sub-step: the sideways heat if given, the implicit conduction (Newton in the temperature with a backtracking line
search per column, stage 1's solver, the tangent the property table's own slope), then — with `remove` and a melting
material — the contiguous top cells at or above T_feed leave with their enthalpy, a column below its minimum draws back
to it, and the cells are re-split. The draw takes rigid metal: at the base temperature, never above T_rigid (ruling 4).
Steady ablation checks it against the exact solution: 2.1 % slow with 20 µm cells and 0.5 ms sub-steps, 1.1 % with
10 µm and 0.25 ms — the first-order error of removing whole cells.

- [ ] **Step 1: Append the tests**

{{TEST_SKIN_T4}}

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: ... no attribute 'advance'` (and `_conduct`).

- [ ] **Step 3: Append to `reentry_model/skin.py`**

{{SKIN_T4}}

- [ ] **Step 4: Run them to see them pass.** Expected: 13 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 362 passed, 1 skipped.

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Layered skins: implicit heat through the depth, liquid removed as it forms, rigid metal drawn back to the minimum"
```

### Task 5: Heat sideways

**Files:**
- Modify: `reentry_model/skin.py` (append), `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Task 1's `Columns` (`.edges`, `.conduct_weight`).
- Produces: `depth_grid(settings, depth)`; `sideways_heat(sk, columns, dt) -> (n, K)` W per cell, the explicit source
  `advance` takes as `side=lambda s, h: sideways_heat(s, columns, h)`.

Each column's cells are averaged onto a common depth grid (spaced by the cell target); every surface edge exchanges
k w l (T_j − T_i) in each bin both columns reach (k the mean of the two bins' conductivities, w the cotangent weight,
l the depth they share); each column's bin heat goes back to its cells by their overlap with the bin. That is
antisymmetric edge by edge — the total is zero to round-off — and exact for a field linear along the surface. The
explicit step is sub-cycled to keep dt k_max Σw / (ρ c_p,min A) below `SIDEWAYS_CFL`; on 2 mm spacing one cycle covers
about 25 ms, so the default 10 ms sub-step needs one.

- [ ] **Step 1: Append the tests**

{{TEST_SKIN_T5}}

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: ... no attribute 'sideways_heat'`.

- [ ] **Step 3: Append to `reentry_model/skin.py`**

{{SKIN_T5}}

- [ ] **Step 4: Run them to see them pass.** Expected: 17 passed in the file (the steady two-dimensional field to
  3·10⁻¹³ K).

- [ ] **Step 5: Unit tier and commit.** Expected 366 passed, 1 skipped.

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Layered skins: explicit sideways conduction between columns, conservative and exact for linear fields"
```

### Task 6: The base law at the corner nodes

**Files:**
- Modify: `reentry_model/thermal/skfem_backend.py`, `reentry_model/thermal/fenicsx_backend.py`,
  `reentry_model/thermal/__init__.py` (`step(..., nodal_law=None)`)
- Test: `tests/test_reentry_model_thermal.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Produces: `solver.step(dt, q_conv, T_amb, ..., nodal_law=(nodes, a, b))` on both backends: heat a_i + b_i T_i [W] into
  each listed node, implicit like the radiation, reported in `StepResult.Q_interface`.

- [ ] **Step 1: Append the tests** to `tests/test_reentry_model_thermal.py`:

{{DIFF_TEST_THERMAL}}

  and to `tests/test_reentry_model_fenicsx.py`:

{{FX_T6}}

- [ ] **Step 2: Run them to see them fail.** Expected: `TypeError: ... got an unexpected keyword argument 'nodal_law'`.

- [ ] **Step 3: Apply the backends' changes**

{{DIFF_SKFEM}}

{{DIFF_FENICSX_BACKEND}}

{{DIFF_THERMAL_INIT}}

- [ ] **Step 4: Run them to see them pass.** Expected: the thermal module 19 passed; FEniCSx `-k nodal_law` 1 passed.

- [ ] **Step 5: Both tiers and commit.** Expected unit tier 367 passed, 1 skipped; FEniCSx 18 passed.

```bash
git add reentry_model/thermal tests/test_reentry_model_thermal.py tests/test_reentry_model_fenicsx.py
git commit -m "Thermal backends: the corner columns' base law as an implicit nodal load"
```

### Task 7: The skin body

**Files:**
- Modify: `reentry_model/skin.py` (append `hand_down`)
- Create: `reentry_model/skinbody.py`
- Test: `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Consumes: Tasks 0–6.
- Produces: `skin.hand_down(sk, m_min, T_rigid) -> (mass, H)`; `skinbody.SkinBody(mesh, material, solver, mass_kg,
  settings=None, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1, 0, 0), deform_depth=DEFORM_DEPTH, remove=True,
  demise_fraction=0.01, repeat_rule="column")` with `.moving`, `.columns`, `.skins`, `.surface` (a `CornerSurface`),
  `.theta`, `.pending`, `.column_removed`, `.removed_mass`, `.removed_enthalpy`, `.drawn_mass`, `.absorbed_mass`,
  `.handed_mass`, `.interface_repeats`, `.jump_faults`, `.absorb_capped`, `.wall`, `.last_trigger`; `advance(t, dt,
  loads, state=None)`, `mass(t=None)`, `energy()`, `energy_balance_residual()`, `mean_temperature()`,
  `surface_temperature()` (the top cells), `radiated_power()`, `field()`, `nose_radius()`, `demised()`,
  `surface_stats()`, `skin_stats()` (the history columns of spec §10.2 that exist at stage 2), `mesh` (property);
  `REPEAT_RULES = ("column", "total", "never")`, `ABSORB_MAX`, `GRADIENT_FLOOR`.

At time zero each column takes its minimum from the 3D model: the corners move in by it and the skins hold exactly what
they swept. The macro step (spec §4, §7.2): the trial pass gives each column its base law; the 3D step takes that law
at the corner nodes with no convective load and no radiation on its boundary, and last step's deferred mismatch as a
nodal load; the real pass books the heat, removes the liquid and draws; the mismatch is deferred to the node, and the
step repeats once, re-linearised, under the repeat rule. Then the base moves: draw (already in the skins), absorb (the
base node's excess over T_rigid along the 3D gradient, capped at `ABSORB_MAX` of the mean edge) and hand-down become one
inward displacement per corner, plus the volume-neutral smoothing; `move_with_policy` applies it, and each column is
reconciled with what its corner actually swept — its bottom gains the swept mass less what it already exchanged, its
enthalpy the swept enthalpy less what it already exchanged — so mass and energy are exact. Heating and radiation use the
top surface (the bases moved out by the skins' thicknesses).

- [ ] **Step 1: Write the tests**

{{TEST_SKINBODY}}

  and append to `tests/test_reentry_model_fenicsx.py`:

{{FX_T7}}

- [ ] **Step 2: Run them to see them fail.** Expected: `ImportError: cannot import name 'skinbody'`.

- [ ] **Step 3: Append `hand_down` to `reentry_model/skin.py`**

{{SKIN_T7}}

- [ ] **Step 4: Create `reentry_model/skinbody.py`**

{{SKINBODY}}

- [ ] **Step 5: Run them to see them pass.** Expected: 6 passed (about 70 s); FEniCSx `-k skin_body` 1 passed (on
  that test's two 0.5 s steps the backends' round-off tipped 2 of 2 254 columns across a split or removal threshold;
  mass, removal and energy agree to 1e-9).

- [ ] **Step 6: Both tiers and commit.** Expected unit tier 373 passed, 1 skipped; FEniCSx 19 passed.

```bash
git add reentry_model/skin.py reentry_model/skinbody.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_fenicsx.py
git commit -m "Skin body: skins on every corner from time zero, the trial and real passes around the 3D step, and the base moving the mesh with exact books"
```

### Task 8: The command line and the run's outputs

**Files:**
- Modify: `reentry_model/cli.py` (`SURFACE_MODELS`, the skin flags, `skin_name_suffix`, `model_run_name`,
  `build_thermal`, `skin_settings`, the checks, the printed summary)
- Modify: `reentry_model/coupled.py` (`SKINBODY_COLUMNS`, `CoupledRun.skins`, `skin_results`, demise, `write_skin_frame`)
- Modify: `reentry_model/body.py` (stage 1's value becomes `"skin1d"`)
- Test: `tests/test_reentry_model_cli.py` (stage 1's three CLI tests replaced), `tests/test_reentry_model_melting_skin.py`
  and `tests/test_reentry_model_fenicsx.py` (`surface_model="skin1d"`)

**Interfaces:**
- Produces: `--surface-model elements|skin` (default `elements`), `--skin-min-thickness` (mm, 0.4), `--skin-cell` (µm),
  `--skin-growth` (1.10), `--skin-substep` (ms, 10), `--skin-preset dev|production`, `--deform-depth` (mm, 15); run-name
  suffix `_surface-skin` plus `skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm)`
  (`_skin-dev`, `_skinmin<mm>mm`, `_cell<um>um`, `_growth<r>`, `_sub<ms>ms`, `_deform<mm>mm`); `--surface-model skin`
  needs `--thermal fem` and, with `--melt on`, `--removal instant`; `--animate`/`--stills` are refused with it.
  `coupled.SKINBODY_COLUMNS` in every skin run's history; results `skin_columns`, the totals, `surface_fold_max_deg`,
  `surface_quality_min`, `zone_ratio_min`; frames `surface_<k>.vtp` with per-corner point data `T_top`, `q_conv`,
  `skin_thickness`, `skin_cells`, `liquid_depth`, `flowing_depth`, `base_position`.

- [ ] **Step 1: Apply the test changes**

{{DIFF_TEST_CLI}}

{{DIFF_MELTSKIN_B}}

{{DIFF_FENICSX_D}}

- [ ] **Step 2: Run them to see them fail.** Expected: the new CLI tests fail (`--skin-min-thickness` unrecognised and
  `skin_name_suffix` taking four arguments); stage 1's tests fail on `surface_model must be one of ('elements', 'skin')`.

- [ ] **Step 3: Apply the changes**

{{DIFF_BODY_B}}

{{DIFF_CLI_B}}

{{DIFF_COUPLED}}

- [ ] **Step 4: Run them to see them pass.** Expected: the CLI module 52 passed (the skin run about 40 s).

- [ ] **Step 5: Both tiers.** Expected unit tier 377 passed, 1 skipped; FEniCSx 19 passed.

- [ ] **Step 6: The element model is unchanged.** Run the same short melting flight from the base and from the work
  copy and compare the histories:

```bash
for d in work-2026-10-12-moving-mesh work-2026-10-12-skin-heat; do (cd ../$d/code && "$PY" -m reentry_model run --diameter 100 --velocity 7.24 --altitude 71 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --heating sesam --melt on --prism-layers 0 --h-surface 4 --h-core 20 --temperature 800 --t-max 3 --dt-continuum off --quiet --outdir /tmp/bitid_$d); done
cmp /tmp/bitid_work-2026-10-12-moving-mesh/*.csv /tmp/bitid_work-2026-10-12-skin-heat/*.csv && echo identical
```

  Expected: `identical` (checked while planning).

- [ ] **Step 7: Commit**

```bash
git add reentry_model/cli.py reentry_model/coupled.py reentry_model/body.py tests/test_reentry_model_cli.py tests/test_reentry_model_melting_skin.py tests/test_reentry_model_fenicsx.py
git commit -m "Command line: --surface-model skin with its settings, run names, history columns, results and frames; stage 1 reachable as skin1d"
```

### Task 9: Level 2 — the plate

**Files:**
- Create: `analysis/ablation_1d.py` (the independent one-dimensional reference), `analysis/skin_plate.py` (the plate),
  `data/reference_values/box_q6p_fine.json` (sub-plan 18's fine box, extracted)

Level 2 (spec §11.1 as amended): a plate 8 mm deep and 110 × 110 mm across, of 2 mm elements at 700 K, Scheil, no
radiation, its front heated at 6 MW/m² within 32 mm of its centre and tapering to zero at 45 mm, for 2 s, liquid removed
as it forms. The centre's recession within 5 % at 1.0, 1.5 and 2.0 s of the 1 µm one-dimensional reference and of the
fine box; mass to 1e-12, energy to 1e-8.

- [ ] **Step 1: The one-dimensional reference.** Create `analysis/ablation_1d.py`:

{{ABLATION}}

  Run `"$PY" analysis/ablation_1d.py --out ../measure/ablation_1d_1um.json` (about 12 s). Expected (measured while
  planning): 0.7890, 2.7120 and 5.2190 mm at 1.0, 1.5 and 2.0 s; at 5 µm and 1 ms 0.7900, 2.7100, 5.2150.

- [ ] **Step 2: The fine box's series.** Save as `../extract_fine_box.py` beside the work copy and run it from the work
  copy's root (the run file is git-ignored, so it is read from the main checkout's
  `prototype/work-2026-10-08-skin/analysis/box/runs/q6p__fine.json`):

{{EXTRACT}}

  The fine box's recession: 0.8045, 2.7290 and 5.2060 mm at 1.0, 1.5 and 2.0 s, within 2 % of the reference.

- [ ] **Step 3: The plate.** Create `analysis/skin_plate.py`:

{{PLATE}}

- [ ] **Step 4: Run it** — `caffeinate -i "$PY" analysis/skin_plate.py --out ../measure/plate_dev.json` (about 3 min),
  and with `--preset production` (about 5 min). Expected, measured while planning:

The comparison below is the criterion: every entry within 5 % of both references, and the production cells agree with
the development cells to 1.0 % or better.

| Cells | 1.0 s | 1.5 s | 2.0 s | Books | Wall time |
|---|---|---|---|---|---|
| 20 µm (dev) | 0.8039 mm (+1.9 %) | 2.7097 mm (−0.1 %) | 5.1817 mm (−0.7 %) | mass 1e-15, energy 2.4e-9 | 193 s |
| 10 µm (production) | 0.8120 mm (+2.9 %) | 2.7085 mm (−0.1 %) | 5.1578 mm (−1.2 %) | mass 1e-15, energy 2.4e-9 | 317 s |

  Percentages against the one-dimensional reference; against the fine box the development run is −0.1, −0.7 and
  −0.5 %, the production run +0.9, −0.8 and −0.9 %. The surface's largest fold stays at 93 degrees, the plate's own unheated edges; 2 rebuilds; no base faults.
  How the plate was reached, for the record: with full flux only to 12 mm the centre lost heat sideways and lagged 9 %
  by 2 s; on the wider plate with the old rebuild interpolation it ran 9.5 % ahead (ruling 1); 1 mm elements put the
  melt onset right where 2 mm elements are 12 % late on the narrow plate.

- [ ] **Step 5: Commit**

```bash
git add analysis/ablation_1d.py analysis/skin_plate.py data/reference_values/box_q6p_fine.json
git commit -m "Level 2: the plate with tapered heating against the 1 um one-dimensional reference and the fine box"
```

### Task 10: What triggers the interface repeat, and what a step costs

**Files:**
- Create: `analysis/skin_timing.py`

This task measures (Asha, 2026-10-12); it changes no default. The results below were measured while planning; the task
repeats them on the built code and records them as facts, and Asha decides the interface rule from them.

- [ ] **Step 1: The repeat rules on the plate.** Run `analysis/skin_plate.py --repeat-rule column|total|never` (the
  three can run side by side, about 3 min each). Measured while planning:

The table compares the three rules on the same plate: the spec's rule and the total rule repeat every step and are the
most accurate; no repeat saves 40 % of the time but puts the centre 6.7 % ahead at 1.0 s, outside the 5 % criterion.

| Rule | 1.0 s | 1.5 s | 2.0 s | Repeats | Wall time |
|---|---|---|---|---|---|
| column (spec §7.2, default) | +1.9 % | −0.1 % | −0.7 % | every step | 193 s |
| total | +1.9 % | −0.1 % | −0.8 % | every step | 194 s |
| never | +6.7 % | +2.0 % | −2.4 % | none | 114 s |

  What triggers it (`repeat_columns_over`, `mismatch_worst_ratio`, `mismatch_total_ratio` in the history): a median of
  6 783 of 6 932 columns exceed 10⁻³ of their throughput every step, most of them columns with almost no heat through
  them, so the per-column rule always repeats; the total mismatch is a median 3.4 % of the total throughput (at most
  6.4 %), so a rule on the total repeats every step too.

- [ ] **Step 2: Stability without the repeat.** On the 4 mm test sphere from 880 K with 0.5 s steps, no repeat runs away
  at 3 MW/m² (the 3D base overshoots T_feed, and every sub-step removes and redraws a whole skin: 0.13 kg in one step);
  the repeating rules do not. A predictor for the linearisation temperature (the last step's rate) cut the plate's
  total mismatch from 3.4 % to 0.8 % but ran away at 2 MW/m² on the test sphere, with or without the repeat, and was
  dropped. Repeat this check with `test_every_repeat_rule_keeps_the_books_and_records_what_triggered_it` and a 3 MW/m²
  variant, and record it.

- [ ] **Step 3: The cost of a step on the 100 mm sphere.** Create `analysis/skin_timing.py`:

{{TIMING}}

  Run `caffeinate -i "$PY" analysis/skin_timing.py ../measure/timing_100mm.json` (about 10 min). Measured while
  planning (9 417 columns, 10 µm cells, 0.5 s macro steps, the per-column rule):

The table gives the wall time of each 0.5 s macro step: nearly all of it is in the skins' passes, and every step
repeated once.

| Phase | Cells | Per macro step | 3D solve | Moving the mesh |
|---|---|---|---|---|
| below melting (300 K, 1 MW/m²) | 374 708 | 55–57 s | 0.45 s | 0.50 s |
| melting (880 K, 2 MW/m²) | 449 911 | 64, 110 and 131 s, rising as more columns melt | 0.77 s | 0.65 s |

- [ ] **Step 4: Level 3, the Step 2 check against SESAM with skins** (spec §11.1). At these costs each flight takes
  hours; run them in the build, not before: the Step 2 verification cases of `analysis/reentry_model_thermal_verification.py`
  with `--surface-model skin` (non-melting, so no `--melt`), against Step 2's thresholds in
  `tests/test_reentry_model_reference_thermal.py`; record the differences from the element model's results.

- [ ] **Step 5: Record and commit.** Add the measurements of Tasks 0, 9 and 10 to `00-shared-context.md` as facts, and
  write the Results section of this sub-plan.

```bash
git add analysis/skin_timing.py
git commit -m "Skins: the interface repeat's triggers and the cost of a step measured"
```

## Notes for sub-plans 21 and 22

- Sub-plan 21 adds the cells' velocity, the surface gradient and curvature operators (ruling 7), and level 2 with flow
  on — on the plate, not the box.
- The interface rule (Task 10) is Asha's decision before sub-plan 21 makes a step more expensive.
- Sub-plan 22's retirement of stage 1 is now the deletion of `skin1d.py`, its tests and `MeltingBody`'s `skin1d` path.
- The moves' remap still closes its books by one uniform offset on the moved nodes: 442 J over the plate's 42 moves,
  0.8 % of its absorbed heat, one-signed. It is reported every step; make it local if a flight shows it matters.
