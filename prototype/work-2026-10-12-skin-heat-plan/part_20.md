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
and checks that the result is the tested copy byte for byte. The unit tier goes from 341 to 403 passed (1 skipped, no
failures; 397 by Task 8, six more in Task 10 Step 4); the FEniCSx tier from 17 to 19. Measurements quoted were each made once, on the case named, on 2026-10-12.

**Built** 2026-10-12, subagent-driven, in `prototype/work-2026-10-12-skin-heat/code/` (branch `skin-heat`,
70ecbdd..f5d6633), then **amended after the build** (a45f204, fb4fe25, d1b5636, and on 2026-10-14 the absorption rule
and the hot plate, then the skins' swelling setting for levels 2 and 2b, on 2026-10-15 the fine band under the
surface, and on 2026-10-16 the skins' 0.8 mm minimum with a 0.4 mm fine zone and the part cap's look-ahead) by Asha's
decisions
on the interface, below: damped repetition is the default, the metal crossing the base is valued along the 3D profile,
a macro step runs as parts so that no base moves more than a quarter of the mean surface edge in one and, while
melting, no part is longer than 0.0625 s, with the heating interpolated in time across the parts; and a base fault is
a bottom cell outside the temperatures of the metal that crossed. The planning copy was then
synced from the code and this document regenerated, so the code below is the code as amended, and the measurements of
Tasks 0, 9 and 10 are the build's or the amendment's, as marked. Results at the end lists the corrections, the
amendment's measurements and what is still open.

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

**Shared later from this sub-plan:** the slurry viscosity law it reads; `columns.py` and `skin.py` are skyscrapers' own and are never shared.

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

Asha's decisions after the build (2026-10-12, from a diagnosis of the interface at 0.5 s macro steps; spec §5.3, §6.6,
§7.1, §7.2 and §16 amended; Results):

- **Damped repetition is the default** (`repeat_rule="damped"`): the step repeats, each pass re-linearised at a
  temperature moved only a fraction of the way to the last 3D answer (Aitken's method), until the mismatch is within
  10⁻³ of the throughput. The one-repeat rules stay selectable. Task 7.
- **The metal crossing the base is valued along the 3D model's linear profile below it** (option a): each draw takes
  the next layer down at its own enthalpy (`skin.profile_enthalpy`, Task 4), and the moving mesh books an inward sweep
  the same way (`MovingMesh.move(u, h_inward)`, Task 0).
- **A macro step runs as parts** (option b) so that no base moves more than `MOVE_MAX`, a quarter of the mean surface
  edge, in one part by drawing and handing down. Task 7.
- **Absorption is limited only by what a part's draws and hand-downs leave of its move allowance** (2026-10-14, option 1
  after investigating absorption): the cap of a quarter-edge per macro step stranded non-rigid metal in the 3D model
  and made the plate started at 880 K melt 4-8 % slow. Coarse surface meshes are not supported under strong melting
  (the 4 mm test sphere is refused at the rim of a melting region at 6 MW/m²; the unit test heats it at 3 MW/m², and
  the strong case runs on the production mesh in the reference tier). Task 7.
- **The hot plate is level 2b** (2026-10-14): the plate started at 880 K under 2 MW/m², against its one-dimensional
  reference's recession and non-rigid depth. Task 9.

- **A base fault is a bottom cell outside the temperatures of the metal that crossed** the base, between the node's and
  the 3D profile's at the depth swept (`skinbody.jump_faults`), not one far from the node alone. Task 7.
- **While melting, no part is longer than `PART_MAX` (0.0625 s)**, under a 0.5 s macro step for the trajectory and the
  heating, and **each part takes the heating interpolated in time** between the last macro step's and this one's at
  its middle (option 1, after a diagnosis of the cost: long parts during melting were rejected for their moves and
  needed many passes; the heating's shortest time scale is about 7 s). Task 7.
- **The cap also covers the macro step in which a top is about to become non-rigid** (2026-10-16, with option B′'s
  0.8 mm minimum): a top's rise over the last macro step, carried on across this one, reaching T_rigid caps the step
  too, counted in `SkinBody.capped_ahead`; without it the level-2 cold plate at 0.5 s steps ran the step in which it
  melted as one part and was 14.2 % ahead at 1.0 s (facts 145–146). Task 7.

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
   Drawn metal capped at T_rigid carries no sensitivity to the base temperature (the build found that carrying 1 zeroed
   the base law's slope above T_rigid; `add_bottom(dm, h, s)`). Task 4.
5. **The skins look their material up in 0.01 K tables** (`PropertyTable`), and every pass starts from each cell's stored
   enthalpy, so the books are exact whatever the table's error (at most about 60 J/kg near a Scheil node). The property
   calls were half of a pass's time. Task 2.
6. **Only the columns that need it are re-split or merged**, and the Newton tangent is the table's own slope: together
   with the tables, a macro step on the 4 mm test sphere went from 11 s to 8 s. Tasks 2 and 4.
7. **The surface gradient and curvature operators and the cells' velocity are left to sub-plan 21**, the first stage to
   use them; nothing here could test them.
8. **`SkinBody` is a module of its own**, `skinbody.py`: `body.py` holds the element model, which must stay bit for bit.
9. **The interface kept the spec's per-column repeat as the default**, with the other rules selectable for Task 10's
   measurement. A predictor for the linearisation temperature was tried and dropped: it ran away on 4 mm elements.
   Superseded after the build: the coupling was wrong at 0.5 s macro steps, and damped repetition is the default
   (above).
10. **A skin run's frames carry per-corner data**; `--animate` and `--stills`, which read the element model's frames, are
    refused with skins. Task 8.
11. **Skin runs mesh without prism layers and with the dense band**, which the moving mesh's rebuild reproduces; an
    explicit `--prism-layers` above 0 is refused with exit code 2 (added in the build). Task 8.
12. **Level 3's Step 2 runs against SESAM are the build's**, not this plan's planning: at the measured cost each takes
    hours (Task 10). Still to run: held while the coupling was wrong at 0.5 s, and now a question of cost (Results).

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

1. **A coarse 3D base overshooting T_rigid or T_feed within a macro step** — the draw is capped at T_rigid, the
   damped repetition converges the interface, and the parts keep every move within the mesh's reach (without them,
   3 MW/m² on 4 mm elements asked for a 4.65 mm corner move that the mesh mover refused). Pinned by
   `test_damped_repetition_converges_where_the_undamped_repeat_parks_the_heat`,
   `test_a_step_that_would_move_a_base_too_far_runs_as_parts_and_matches_the_shorter_steps` (Task 7),
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
  transfer, new `project_enthalpy`, `QUADRATURE`, `PROJECT_CHUNK`; the fine band, added 2026-10-15: `rebuild_volume`'s
  `fine_depth` and `h_fine`, `default_sizes` reading them, new `refined_mesh`)
- Modify: `analysis/movemesh_recession.py` (`--neutral`)
- Test: `tests/test_reentry_model_movemesh.py` (three tests changed for the projection, three appended)

**Interfaces:**
- Produces: `swept_volumes(old, new, faces, faces_corner, n_corners)` — same signature, each corner booked its own
  sweep (spec §5.2 as amended); `MovingMesh.neutral_smoothing(normal=None, area=None, passes=NEUTRAL_PASSES) -> (nc, 3)`
  (the tangential smoothing plus the normal corrections that cancel what it sweeps per corner);
  `project_enthalpy(old_points, old_tets, h_old, points, tets) -> (n,)`, the rebuild's transfer;
  `rebuild_volume(corner_points, faces_corner, h_surface, h_core, band, ramp, fine_depth=0.0, h_fine=None)`;
  `refined_mesh(mesh, fine_depth, h_fine) -> VolumeMesh` (its parameters `fine_depth_m`, `h_fine_m`).

Why: on the level-2 box the third-of-a-prism booking credited the side walls with metal the front swept (Decisions).
With own sweeps, the plain smoothing books at most 0.46 µm at a corner of the 100 mm sphere (1.03 µm before), and three
normal corrections bring that to 3·10⁻¹¹ m (0.37 µm before), each pass removing about 95 % of what is left. The rebuild:
interpolating at the new nodes and offsetting the total uniformly over the zone moved heat from the cold bulk to the
melting front (2 768 J in two rebuilds of the plate, 4.9 % of its absorbed heat); the L2 projection of the old field —
the new mesh's consistent mass matrix against the old field integrated over the new elements by the 4-point rule — keeps
each region's energy to the rule's accuracy and reproduces a linear field exactly (a lumped projection, also tried,
moved a linear field's surface values by up to 7 K). On the tests' kinked field the residual goes from −4.1 % of the
energy to −1.8 J of 63 219 J. The fine band (Asha, 2026-10-14, option 1; added 2026-10-15, facts 138–140): on 2 mm
elements the 3D model cannot follow the temperature just below the rigid front and the plates melt 2–6 % late, so a
mesh may be finer to `fine_depth` below its surface while the surface keeps its triangulation — gmsh fills the volume
under it to that size field, at the start (`refined_mesh`) and at every rebuild, the sizes read from the mesh's
parameters — pinned by `test_a_refined_mesh_keeps_its_surface_and_is_finer_to_the_fine_depth` and
`test_a_rebuild_makes_the_fine_band_again`. Prism layers pushed in from the surface were tried first and could not be
built at the plate's corners (fact 138).

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

- [ ] **Step 4: Run them to see them pass.** Expected: 38 passed. The amendment added `move(u, h_inward)`, which books
  an inward sweep at the enthalpy the caller gives for each corner's swept depth (the skins pass the 3D profile's
  average), pinned by `test_an_inward_sweep_is_booked_at_the_enthalpy_the_caller_gives_for_its_depth`. The build added
  the retry in `locate`
  (`LOCATE_RETRY = 256`): 0.3–0.4 % of a rebuild's quadrature points, all interior, had been placed in a wrong element
  and extrapolated from it (worst barycentric coordinate −0.27 on the 100 mm sphere); they are now retried with 256
  candidates, pinned by `test_locate_finds_the_containing_element_for_every_quadrature_point_of_a_rebuild`.

- [ ] **Step 5: The recession measurement's switch** (`analysis/movemesh_recession.py`):

{{DIFF_RECESSION}}

  Run `caffeinate -i "$PY" analysis/movemesh_recession.py --neutral --out ../measure/recession_neutral.json` (about
  75 s). Measured while planning, with the third-of-a-prism booking: the run fails at 29.5 mm of nose recession against
  30.0 mm with the plain smoothing, the worst surface triangle quality 0.711 in both, the rebuilds' residuals 0.07–0.26 %
  of the total enthalpy in both. Record the result with own sweeps and the projection beside those. Built: it still fails
  at 30.0 mm (29.5 mm completed: the fold of the prescribed motion law, fact 104), and the rebuilds' residuals fall to
  4·10⁻⁶–1.6·10⁻⁵ of the total enthalpy, in 167 s (20 rebuilds of 6.2–7.5 s each).

- [ ] **Step 6: Unit tier and commit.** Expected 348 passed, 1 skipped.

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

- [ ] **Step 5: Unit tier and commit.** Expected 353 passed, 1 skipped.

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
  `SIDEWAYS_CFL`, `TABLE_DT`, `TABLE_RANGE`; `SkinSettings(min_thickness=0.8e-3, cell=10e-6, fine_zone=0.4e-3,
  growth=1.10, substep=0.010, swell=True)` (`swell=False` holds every cell at the solid's density, the setting levels 2
  and 2b verify with, added 2026-10-14; the minimum 0.4 mm and the fine zone the minimum's depth until option B′,
  2026-10-16, set them apart); `thomas(sub, diag, sup, rhs)`; `PropertyTable(material)` with `enthalpy`, `cp_eff`, `k`,
  `liquid_fraction`, `conduction(T) -> (h, dh/dT, k)`; `Skins(material, area, settings)` with `.m`, `.H`, `.T` (n, K;
  cell 0 the top), `.count`, `.area`, `.S` (the trial pass's dT/dT_b, or None), `.props`, `.valid`, `mass()`,
  `energy()`, `column_mass()`, `state()`, `set_state()`, `density(T)`, `dz()`, `thickness()`, `targets(dz=None)`,
  `fill(mass, T)`, `regrid()`, `take_top(k) -> (mass, H)`, `add_bottom(dm, h) -> H`, `take_bottom(dm) -> H`,
  `exchange_bottom(dm, dH)`, `depth_above(T_level)`.

Cells hold mass and enthalpy; temperature, liquid fraction and density follow from the material; a cell's thickness is
its mass over the column's area and its own density. A cell's target is the fine cell within the fine zone (`fine_zone`
below the top) and grows by (growth − 1) per unit depth below it, the continuous form of each cell
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

  Run the unit tier: expected 353 passed, 1 skipped (stage 1's 21 tests now in `test_reentry_model_skin1d.py`).

- [ ] **Step 2: Write the new skin's tests** (`tests/test_reentry_model_skin.py`, which grows by a section per task):

{{TEST_SKIN_T2}}

- [ ] **Step 3: Run them to see them fail.** Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `ImportError: cannot import name 'skin' from 'reentry_model'`.

- [ ] **Step 4: Create `reentry_model/skin.py`**

{{SKIN_T2}}

- [ ] **Step 5: Run them to see them pass.** Expected: 7 passed — once Task 3's `minimum_mass` exists, because these
  tests' helper calls it. Do Tasks 2 and 3 together and commit once, after Task 3 (the build did).

- [ ] **Step 6: Unit tier and commit** (with Task 3, below). Expected 362 passed, 1 skipped, after Task 3.

```bash
git add -A reentry_model/skin.py reentry_model/skin1d.py reentry_model/body.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skin1d.py tests/test_reentry_model_melting_skin.py tests/test_reentry_model_fenicsx.py
git commit -m "Layered skins: stage 1 moves to skin1d, and the new skin's cells split and merge exactly"
```

### Task 3: The minimum mass

Done in the same pass as Task 2, whose tests need `minimum_mass`; one commit for both.

**Files:**
- Modify: `reentry_model/skin.py` (append), `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Produces: `minimum_mass(material, area, settings, delta_m=None) -> (n,)` kg: the larger of `min_thickness` of solid
  over the column's area and the liquid filling Girin's δ_m plus `DELTA_M_MARGIN` (spec §6.1); δ_m None or nan leaves the
  solid term alone (stage 2 evaluates no surface flow, so it passes None).

The minimum is 800 µm of solid (option B′, Asha, 2026-10-16; 400 µm before): where the non-rigid layer is thinner the
skin's bottom is rigid metal and the skin-3D interface lies below the rigid front, off the steepest part of the
temperature profile, while the fine zone keeps the finest cells to the top 400 µm, so a column at its minimum holds 59
cells instead of 80 (facts 141–145). Pinned by `test_the_minimum_is_800_um_of_solid_unless_girins_layer_needs_more`
and `test_a_skin_at_its_minimum_holds_fine_cells_only_in_its_fine_zone`; the cell tests keep a 0.4 mm skin of fine
cells (`CELLS_04`), the geometry they were written for.

- [ ] **Step 1: Append the test**

{{TEST_SKIN_T3}}

- [ ] **Step 2: Run it to see it fail.** Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'minimum_mass'`.

- [ ] **Step 3: Append to `reentry_model/skin.py`**

{{SKIN_T3}}

- [ ] **Step 4: Run it to see it pass.** Expected: 9 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 362 passed, 1 skipped (Tasks 2 and 3).

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Layered skins: the minimum mass, 400 um of solid or Girin's layer and 50 um as liquid"
```

### Task 4: Heat through the depth, with liquid removed as it forms

**Files:**
- Modify: `reentry_model/skin.py` (append), `tests/test_reentry_model_skin.py` (append)
- Create: `data/reference_values/skin_swell_1d.json` (the swelling test's reference, added 2026-10-14: what Task 9's
  `analysis/ablation_1d.py --t0 880 --q 2e6 --length-mm 16 --t-end 1.0 --swell off|all` writes, the two runs' recession
  in one file)

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
to it, and the cells are re-split. The draw takes rigid metal: at the base temperature, never above T_rigid (ruling 4),
and metal drawn at the cap carries no sensitivity to the base (`add_bottom(dm, h, s)` with s 0 where capped), pinned by
`test_the_trial_law_keeps_its_slope_when_the_base_is_above_t_rigid` against a ±0.05 K finite difference. Given the 3D
model's gradient below the base (`g_base`, amended after the build), each draw takes the next layer down the 3D
model's linear profile at that layer's enthalpy, `profile_enthalpy` (an 8-point Gauss average, within 12 J/kg of a
200 001-point one on the Scheil curve), pinned by `test_the_profile_enthalpy_is_the_mean_along_the_3d_profile_and_a_draw_takes_it`.
Steady ablation checks it against the exact solution: 2.1 % slow with 20 µm cells and 0.5 ms sub-steps, 1.1 % with
10 µm and 0.25 ms — the first-order error of removing whole cells. The skins' swelling is checked the same way (level 1,
added 2026-10-14, Asha): one column 8 mm deep started at 880 K under 2 MW/m² melts sooner with each cell at its own
phase's density than with `swell=False`, by the amount the one-dimensional slab does with the same two density laws
(0.113 mm on average over 0.7–1.0 s; the skins' advance 1.3 % above it, each run within 2.4 % of its own reference),
`test_the_skins_swelling_melts_a_column_sooner_as_it_does_the_exact_slab_with_their_density_law`. Its reference values:

{{SKIN_SWELL_REF}}

- [ ] **Step 1: Append the tests**

{{TEST_SKIN_T4}}

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: ... no attribute 'advance'` (and `_conduct`).

- [ ] **Step 3: Append to `reentry_model/skin.py`**

{{SKIN_T4}}

- [ ] **Step 4: Run them to see them pass.** Expected: 17 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 370 passed, 1 skipped.

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py data/reference_values/skin_swell_1d.json
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

- [ ] **Step 4: Run them to see them pass.** Expected: 21 passed in the file (the steady two-dimensional field to
  3·10⁻¹³ K).

- [ ] **Step 5: Unit tier and commit.** Expected 374 passed, 1 skipped.

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
  each listed node, implicit like the radiation, reported in `StepResult.Q_interface`; b is expected to be at most zero.
  A law on a pinned node is dropped in both backends and its power at the old temperatures reported in
  `StepResult.Q_dropped` (added in the build, where a pinned node broke the scikit-fem books), pinned by
  `test_a_nodal_law_on_a_pinned_node_is_dropped_and_reported`.

- [ ] **Step 1: Append the tests** to `tests/test_reentry_model_thermal.py`:

{{DIFF_TEST_THERMAL}}

  and to `tests/test_reentry_model_fenicsx.py`:

{{FX_T6}}

- [ ] **Step 2: Run them to see them fail.** Expected: `TypeError: ... got an unexpected keyword argument 'nodal_law'`.

- [ ] **Step 3: Apply the backends' changes**

{{DIFF_SKFEM}}

{{DIFF_FENICSX_BACKEND}}

{{DIFF_THERMAL_INIT}}

- [ ] **Step 4: Run them to see them pass.** Expected: the thermal module 20 passed; FEniCSx `-k nodal_law` 1 passed.

- [ ] **Step 5: Both tiers and commit.** Expected unit tier 376 passed, 1 skipped; FEniCSx 18 passed.

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
  demise_fraction=0.01, repeat_rule=DEFAULT_REPEAT_RULE)` with `.moving`, `.columns`, `.skins`, `.surface` (a `CornerSurface`),
  `.theta`, `.pending`, `.column_removed`, `.removed_mass`, `.removed_enthalpy`, `.drawn_mass`, `.absorbed_mass`,
  `.handed_mass`, `.interface_repeats`, `.jump_faults`, `.absorb_capped`, `.wall`, `.last_trigger`, `.last_passes`,
  `.last_total`, `.step_parts`; `advance(t, dt,
  loads, state=None)`, `mass(t=None)`, `energy()`, `energy_balance_residual()`, `mean_temperature()`,
  `surface_temperature()` (the top cells), `radiated_power()`, `field()`, `nose_radius()`, `demised()`,
  `surface_stats()`, `skin_stats()` (the history columns of spec §10.2 that exist at stage 2), `mesh` (property);
  `REPEAT_RULES = ("column", "total", "never", "damped")`, `DEFAULT_REPEAT_RULE = "damped"`, `DAMPED_MAX_PASSES`,
  `DAMPED_OMEGA0`, `DAMPED_OMEGA_RANGE`, `DAMPED_TOL_K`, `MOVE_MAX`, `SPLIT_DEPTH_MAX`, `MOVE_HINT_MARGIN`, `PART_MAX`,
  `jump_faults(T_bottom, T_node, T_deep, crossed, tol=None)`; the constructor's `part_max=PART_MAX` (None: no cap),
  `SKIN_CG_TOL` (the 3D solver's linear tolerance in a skin run), `GRADIENT_FLOOR`.

At time zero each column takes its minimum from the 3D model: the corners move in by it and the skins hold exactly what
they swept. The macro step (spec §4, §7.2): the trial pass gives each column its base law; the 3D step takes that law
at the corner nodes with no convective load and no radiation on its boundary, and last step's deferred mismatch as a
nodal load; the real pass books the heat, removes the liquid and draws; the mismatch is deferred to the node, and the
step repeats under the repeat rule: by default (amended after the build) until the mismatch is within the tolerance,
each pass re-linearised at a temperature moved only a fraction ω of the way to the last 3D answer, Aitken's method
choosing ω; the draws take the 3D profile below each base, its slope from the part's start. A macro step runs as
equal parts (amended after the build): a part whose draws, absorption and hand-downs would move a base more than
`MOVE_MAX` of the mean edge is restored and run again as shorter parts, and the next macro step starts with as many
parts as this one's largest move called for; while any skin's top is non-rigid at the macro step's start, no part is
longer than `part_max`. Each part takes the heating interpolated linearly in time between the last macro step's and
this one's, at its middle. Then the base moves: draw (already in the skins), absorb (the
base node's excess over T_rigid along the 3D gradient, at most what the part's draws and hand-downs leave of
`MOVE_MAX` of the mean edge, the rest waiting for the next part) and hand-down become one
inward displacement per corner, plus the volume-neutral smoothing; `move_with_policy` applies it, and each column is
reconciled with what its corner actually swept — its bottom gains the swept mass less what it already exchanged, its
enthalpy the swept enthalpy less what it already exchanged, an inward sweep booked along the profile the draws used —
so mass and energy are exact. Heating and radiation use the
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

- [ ] **Step 5: Run them to see them pass.** Expected: 12 passed and 2 deselected (about 4 minutes: the heating test
  runs its 0.5 s steps as 1, 8 and 8 parts at 3 MW/m²;
  `test_the_absorption_depth_is_the_excess_over_the_gradient_along_the_normal_and_is_capped` was added at the final
  review, the damped, parts, fault-counter and part-cap tests with the amendment, the cap's look-ahead,
  `test_the_step_in_which_a_top_becomes_non_rigid_is_capped_too`, on 2026-10-16; the tests of the repeat rules and of
  the move limit run with `part_max=None`, the conditions they study). The two deselected are the reference tier's
  (2026-10-14): `-m reference` runs `test_strong_melting_on_the_production_mesh` (6 MW/m² on the 2 mm production
  mesh, about 11 minutes) and `test_the_hot_plate_passes_level_2b` (about 3 minutes); FEniCSx `-k skin_body` 1 passed (one 0.5 s step,
  5 parts, no rebuild: the backends agree to round-off, mass 8e-12, removed mass 4e-10; across a rebuild they do not,
  because the mesher answers round-off in the corner positions with a different mesh).

- [ ] **Step 6: Both tiers and commit.** Expected unit tier 388 passed, 1 skipped; FEniCSx 19 passed.

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
- Produces: `--surface-model elements|skin` (default `elements`), `--skin-min-thickness` (mm, 0.8; 0.4 until option B′),
  `--skin-fine-zone` (mm, 0.4; added 2026-10-16), `--skin-cell` (µm),
  `--skin-growth` (1.10), `--skin-substep` (ms, 10), `--skin-preset dev|production`, `--deform-depth` (mm, 15),
  `--skin-swell on|off` (on; added 2026-10-14), `--fine-band` (mm, 4; 0 none) and `--fine-size` (mm, 0.7; added
  2026-10-15: a skin run's mesh is `movemesh.refined_mesh` of the sphere); run-name suffix `_surface-skin` plus
  `skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True, fine_band_mm=4.0,
  fine_size_mm=0.7, fine_zone_mm=0.4)` (`_skin-dev`, `_skinmin<mm>mm` other than 0.8, `_skinfine<mm>mm` other than 0.4,
  `_cell<um>um`, `_growth<r>`, `_sub<ms>ms`, `_deform<mm>mm`, `_noswell`, `_fine<depth>x<size>mm`, `_fine0`); `--surface-model skin`
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

- [ ] **Step 4: Run them to see them pass.** Expected: the CLI module 57 passed (the skin run about 40 s; two of them
  pin the refusal of an explicit `--prism-layers` with skins, added at the final review, two a bad fine band and one a bad fine zone).

- [ ] **Step 5: Both tiers.** Expected unit tier 397 passed, 1 skipped; FEniCSx 19 passed. (The history gains
  `step_parts`, the parts the last macro step took, with the amendment.)

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
  `data/reference_values/box_q6p_fine.json` (sub-plan 18's fine box, extracted),
  `data/reference_values/plate_hot_1d.json` (the hot plate's one-dimensional reference, added 2026-10-14)

Level 2 (spec §11.1 as amended): a plate 8 mm deep and 110 × 110 mm across, of 2 mm elements at 700 K, Scheil, no
radiation, its front heated at 6 MW/m² within 32 mm of its centre and tapering to zero at 45 mm, for 2 s, liquid removed
as it forms. The centre's recession within 5 % at 1.0, 1.5 and 2.0 s of the 1 µm one-dimensional reference and of the
fine box; mass to 1e-12, energy to 1e-8. Level 2b, the hot plate (added 2026-10-14, Asha; `--case hot`): the same
plate started at 880 K under 2 MW/m², where the metal between the melt front and T_rigid is 1.5–1.7 mm thick and the
skins' absorption is exercised; the recession within 5 % of its one-dimensional reference at 1.5 and 2.0 s (the onset
at 1.0 s reported), the skins' flowing depth at the centre within 10 % of the exact non-rigid depth, and the centre's
base within 1 K of T_rigid. The harness judges each case and records whether it passes. Both cases run with the skins
at the solid's density (`--swell off`, the harness's default; Asha, 2026-10-14, option A): the references and the fine
box hold their metal at one density, while the skins' own law swells the liquid by 15 % (spec §5.2), which alone melts
the exact plates sooner and had hidden a numerical delay of the same size (facts 133–135); the one-dimensional
reference can be given the skins' law too (`--swell all|nonrigid`), which is how the swelling's own level-1 check
(Task 4) gets its reference. Both plates also have the skin runs' fine band (`--fine-mm 4 --fine-size-mm 0.7`, the
defaults; Asha, option 1, built 2026-10-15): the 3D model's elements finer to 4 mm below the surface, made again at
every rebuild (Task 0, facts 138–140).

- [ ] **Step 1: The one-dimensional reference.** Create `analysis/ablation_1d.py`:

{{ABLATION}}

  Run `"$PY" analysis/ablation_1d.py --out ../measure/ablation_1d_1um.json` (about 12 s). Expected (measured while
  planning): 0.7890, 2.7120 and 5.2190 mm at 1.0, 1.5 and 2.0 s; at 5 µm and 1 ms 0.7900, 2.7100, 5.2150. The hot
  plate's reference: `"$PY" analysis/ablation_1d.py --t0 880 --q 2e6 --out data/reference_values/plate_hot_1d.json`,
  recession 0.789, 1.817 and 2.917 mm and non-rigid depth 1.487, 1.570 and 1.674 mm at 1.0, 1.5 and 2.0 s:

{{PLATE_HOT_REF}}

- [ ] **Step 2: The fine box's series.** Save as `../extract_fine_box.py` beside the work copy and run it from the work
  copy's root (the run file is git-ignored, so it is read from the main checkout's
  `prototype/work-2026-10-08-skin/analysis/box/runs/q6p__fine.json`):

{{EXTRACT}}

  The fine box's recession: 0.8045, 2.7290 and 5.2060 mm at 1.0, 1.5 and 2.0 s, within 2 % of the reference.

- [ ] **Step 3: The plate.** Create `analysis/skin_plate.py`:

{{PLATE}}

- [ ] **Step 4: Run it** — `caffeinate -i "$PY" analysis/skin_plate.py --out ../measure/plate_dev.json` (about 3 min),
  and with `--preset production` (about 5 min). Measured in the build:

The comparison below is the criterion: every entry within 5 % of both references, and the production cells agree with
the development cells to 1.0 % or better. Both hold; the tightest entry is the development run at 1.0 s, 0.83 points
inside the criterion, in a single run.

| Cells | 1.0 s | 1.5 s | 2.0 s | Books | Wall time |
|---|---|---|---|---|---|
| 20 µm (dev) | 0.8219 mm (+4.2 %) | 2.6766 mm (−1.3 %) | 5.1080 mm (−2.1 %) | mass 1.4e-15, energy 2.4e-9 | 163 s |
| 10 µm (production) | 0.8163 mm (+3.5 %) | 2.6978 mm (−0.5 %) | 5.1358 mm (−1.6 %) | mass 1.4e-15, energy 2.4e-9 | 323 s |

  Percentages against the one-dimensional reference; against the fine box the development run is +2.2, −1.9 and
  −1.9 %, the production run +1.5, −1.2 and −1.4 %; the two cell sizes agree to 0.8 %. The surface's largest fold stays
  at 93 degrees, the plate's own unheated edges; 2 rebuilds; no base faults. Planning measured +1.9, −0.1 and −0.7 %
  (dev) on code whose capped draws still carried a sensitivity of 1 (ruling 4) and whose rebuild transfer had no retry
  (Task 0); together the build's corrections moved the onset at 1.0 s 2.3 points earlier and the later times about
  1.3 points later — which of the two did it was not separated. The same development case run again on the same code
  agreed to every printed digit to 1.5 s and gave −2.04 % at 2.0 s against −2.13 %: the run is not bit-reproducible
  after its second rebuild.
  How the plate was reached, for the record: with full flux only to 12 mm the centre lost heat sideways and lagged 9 %
  by 2 s; on the wider plate with the old rebuild interpolation it ran 9.5 % ahead (ruling 1); 1 mm elements put the
  melt onset right where 2 mm elements are 12 % late on the narrow plate.

  **Under the amended default** (damped repetition, the profile valuation and the parts; measured 2026-10-12, beside
  other runs, so wall times are indicative; the skins' work is counted in sub-step passes, each pass times its
  sub-steps), the same plate:

The table shows that at the plate's own 0.05 s steps the production cells pass against both references, the development
cells miss the fine box at 1.5 s by half a point, and the two cell sizes no longer agree to 1 %; and that at the
flight's 0.5 s step the plate passes, where before the amendment it never melted.

| Macro step, cells | 1.0 s | 1.5 s | 2.0 s | Against the fine box | Skins' work |
|---|---|---|---|---|---|
| 0.05 s, 20 µm (dev) | 0.7916 mm (+0.3 %) | 2.5796 mm (−4.9 %) | 5.0481 mm (−3.3 %) | −1.6, −5.5, −3.0 % | 1 346 sub-step passes |
| 0.05 s, 10 µm (production) | 0.7900 mm (+0.1 %) | 2.6162 mm (−3.5 %) | 4.9790 mm (−4.6 %) | −1.8, −4.1, −4.4 % | 781 (each pass twice the cells) |
| 0.1 s, dev | 0.8522 mm (+8.0 %) | 2.6968 mm (−0.6 %) | 5.0291 mm (−3.6 %) | +5.9, −1.2, −3.4 % | 1 630 |
| 0.25 s, dev | 0.8888 mm (+12.7 %) | 2.7227 mm (+0.4 %) | 5.0814 mm (−2.6 %) | +10.5, −0.2, −2.4 % | 2 151 |
| 0.5 s, dev | 0.8245 mm (+4.5 %) | 2.6994 mm (−0.5 %) | 5.0108 mm (−4.0 %) | +2.5, −1.1, −3.7 % | 2 437 |
| 0.5 s, dev, with the melting cap and the heating in time | 0.8245 mm (+4.5 %) | 2.7079 mm (−0.2 %) | 5.0872 mm (−2.5 %) | +2.5, −0.8, −2.3 % | 256 passes |
| 0.05 s, dev, with the heating in time | 0.7916 mm (+0.3 %) | 2.6132 mm (−3.6 %) | 5.0189 mm (−3.8 %) | −1.6, −4.3, −3.6 % | 275 passes |

  The last two rows are the default since option 1 (d1b5636): at 0.5 s macro steps the parts are 1, 4, 8 and 8 a step,
  and both rows pass against both references (energy to 4.5e-10 and 2.3e-9). The rows at 0.1 and 0.25 s were measured
  before the cap, which now makes their melting parts 0.0625 s long too.

  Every run keeps the books (mass to 1.4e-15, energy to 2.4e-9) and leaves at most 0.44 % of a step's heat deferred.
  The two cell sizes differ by 1.4 % at 1.5 and 2.0 s (0.8 % before the amendment). At 0.1 and 0.25 s the excess at
  1.0 s (0.06–0.10 mm) arises in the steps where melting starts and then stays as an offset; taking the profile's slope
  at the step's end instead of its start made it slightly worse (+14.6 % at 0.25 s), so its cause is not found. The
  base-fault counter, when it compared the skin's bottom cell with the 3D node alone, read 0 at 0.05 s and 812–1 694 at
  0.1–0.5 s: the drawn metal comes from below the node, colder by the profile's slope times the depth drawn (about
  0.35 mm at 22 K/mm is 7.7 K, over the 5 K tolerance), so it flagged the valuation itself. Redefined (fb4fe25) to the
  range of the metal that crossed, it reads 0 at every step of the 0.25 s run.

  **The hot plate and absorption** (2026-10-14). Under the quarter-edge cap per macro step the hot plate failed: at
  0.5 s macro steps −7.7, −3.9 and −5.3 % at 1.0, 1.5 and 2.0 s, the skins 0.54–0.59 mm thick against the exact
  non-rigid depth of 1.49–1.67 mm, and the centre's base 7–8 K above T_rigid. With absorption limited only by the move
  allowance (the default since then), measured once each:

The table shows both plates passing at both macro steps with the absorption of option 1.

| Plate, macro step | 1.0 s | 1.5 s | 2.0 s | Other criteria |
|---|---|---|---|---|
| hot, 0.5 s | +1.4 % | +2.5 % | +1.2 % | flowing depth within 3.6 %, base at most 0.08 K above T_rigid |
| hot, 0.05 s | +7.1 % (onset, not judged) | +0.8 % | +1.3 % | flowing depth within 5.9 %, base at most 0.08 K above |
| cold, 0.5 s | +4.5 % | −0.2 % | −2.8 % | fine box +2.5, −0.8, −2.6 % |
| cold, 0.05 s | +0.3 % | −3.6 % | −4.2 % | fine box −1.6, −4.3, −3.9 % |

  **With the skins at the solid's density** (2026-10-14, the default since option A; facts 133–135). The swelling had
  been ahead of the references by itself: given to the non-rigid metal, it melts the exact cold plate 0.033 mm sooner
  (4.2 % at 1.0 s) and the hot plate 0.094 mm sooner (11.9 %). Without it, measured once each on 2 mm elements:

The table shows the hot plate passing level 2b at both macro steps, the cold plate passing level 2 at 0.5 s and failing
it at 0.05 s by a delay of about 6 %, which is the 3D model's resolution in depth (fact 134).

| Plate, macro step | 1.0 s | 1.5 s | 2.0 s | Other criteria |
|---|---|---|---|---|
| hot, 0.5 s | −10.2 % (onset, not judged) | −2.9 % | −2.0 % | flowing depth within 3.4 %, base at most 0.08 K above T_rigid |
| hot, 0.05 s | −6.8 % (onset, not judged) | −3.6 % | −1.7 % | flowing depth within 8.5 %, base at most 0.08 K above |
| cold, 0.5 s | −1.2 % | −2.2 % | −3.9 % | fine box −3.1, −2.8, −3.7 % |
| cold, 0.05 s | −6.2 % | −6.2 % | −6.4 % | fine box −8.0, −6.8, −6.1 % (fails) |

  **With the fine band** (2026-10-15, the default since option 1; facts 139–140), measured once each:

The table shows both plates passing at both macro steps once the 3D model is finer below the surface, the cold plate at
0.05 s turning from 6.2 % late to 1.6 % ahead at 1.0 s.

| Plate, macro step | 1.0 s | 1.5 s | 2.0 s | Other criteria | Wall time |
|---|---|---|---|---|---|
| hot, 0.5 s | −5.6 % (onset, not judged) | −1.9 % | −0.7 % | flowing depth within 0.7 %, base at most 0.10 K above T_rigid | 539 s (363 s) |
| hot, 0.05 s | −2.5 % (onset, not judged) | −1.6 % | −1.0 % | flowing depth within 0.5 %, base at most 0.08 K above | 624 s (357 s) |
| cold, 0.5 s | +4.9 % | −2.5 % | −3.3 % | fine box +2.8, −3.1, −3.0 % | 2 036 s (1 280 s) |
| cold, 0.05 s | +1.6 % | −2.2 % | −2.7 % | fine box −0.3, −2.8, −2.5 % | 2 022 s (1 116 s) |

  The wall times in brackets are without the band. The plates rebuild 4–5 times instead of 1–2, every 1 mm or so of
  recession, because their 4 mm deforming zone compresses the band's elements directly; each rebuild takes about 40 s.

  **With the skins' 0.8 mm minimum and 0.4 mm fine zone, and the cap's look-ahead** (2026-10-16, the default since
  option B′; facts 141–146), measured once each on battery power:

The table shows both plates passing at both macro steps, the cold plate's late slowness halved (from 2.8 to 0.9 % late
at 2.0 s at 0.05 s steps), and the hot plate's onset at 0.5 s later than before (not judged).

| Plate, macro step | 1.0 s | 1.5 s | 2.0 s | Other criteria |
|---|---|---|---|---|
| cold, 0.05 s | +3.6 % | 0.0 % | −0.9 % | fine box +1.6, −0.6, −0.6 % |
| cold, 0.5 s | +3.2 % | +0.4 % | −0.9 % | fine box +1.2, −0.2, −0.6 % (14.2 % ahead at 1.0 s without the look-ahead) |
| hot, 0.05 s | −3.1 % (onset, not judged) | −1.5 % | −0.8 % | flowing depth within 1.2 %, base at most 0.08 K above T_rigid |
| hot, 0.5 s | −9.0 % (onset, not judged) | −2.8 % | −1.1 % | flowing depth within 0.9 %, base at most 0.09 K above |

- [ ] **Step 5: Commit**

```bash
git add analysis/ablation_1d.py analysis/skin_plate.py data/reference_values/box_q6p_fine.json data/reference_values/plate_hot_1d.json
git commit -m "Level 2: the plate with tapered heating against the 1 um one-dimensional reference and the fine box"
```

### Task 10: What triggers the interface repeat, and what a step costs

**Files:**
- Create: `analysis/skin_timing.py`

This task measures (Asha, 2026-10-12); it changes no default. The results below were measured in the build (the
planning copy's in brackets where they differ), and Asha decides the interface rule from them.

- [ ] **Step 1: The repeat rules on the plate.** Run `analysis/skin_plate.py --repeat-rule column|total|never` (the
  three can run side by side, about 3 min each). Measured in the build:

The table compares the three rules on the same plate: the spec's rule and the total rule repeat every step and agree to
1.5 s (they part by 0.16 points at 2.0 s); no repeat saves 44 % of the time but puts the centre 6.7 % ahead at 1.0 s,
outside the 5 % criterion. Wall times are of runs made alone, except the total rule's, made beside the other two.

| Rule | 1.0 s | 1.5 s | 2.0 s | Repeats | Wall time |
|---|---|---|---|---|---|
| column (spec §7.2, default) | +4.17 % [+1.9] | −1.30 % [−0.1] | −2.04 % [−0.7] | every step | 161 s [193] |
| total | +4.17 % [+1.9] | −1.30 % [−0.1] | −2.20 % [−0.8] | every step | 245 s, side by side [194] |
| never | +6.71 % | +1.96 % | −2.38 % | none | 90 s [114] |
| damped (the default since the amendment; before the profile valuation and the parts) | +4.17 % | −1.48 % | +0.34 % | a median of 4 passes a step | 346 s |

  The damped row was measured after the build, alone, before the profile valuation and the parts were added; with
  them the plate at 0.05 s is the first row of Task 9's amended table.

  What triggers it (`repeat_columns_over`, `mismatch_worst_ratio`, `mismatch_total_ratio` in the history): a median of
  6 783 of 6 932 columns exceed 10⁻³ of their throughput every step (never fewer than 6 248), most of them columns with
  almost no heat through them, so the per-column rule always repeats; the total mismatch is a median 3.3 % of the total
  throughput (at most 6.0 %, at least 0.28 % on the first step), so a rule on the total at 10⁻³ repeats every step too.

- [ ] **Step 2: Stability without the repeat.** On the 4 mm test sphere from 880 K with 0.5 s steps, built as
  `tests/test_reentry_model_skinbody.py`'s `make()` builds it (20 µm cells), three steps under `front_loads`. Planning
  found that no repeat runs away at 3 MW/m² (0.13 kg in one step) and that a predictor for the linearisation temperature
  ran away at 2 MW/m²; the predictor was dropped. Measured in the build (removed mass cumulative over the three steps):

The table shows that at 2 MW/m² both rules complete cleanly, and that at 3 MW/m² the repeat is what lets the steps
complete at all, though not accurately.

| Flux | Rule | Outcome | Removed | Repeats | Base faults | Energy residual |
|---|---|---|---|---|---|---|
| 2 MW/m² | column | completed | 0.0118 kg (0, 0.0011, 0.0107 kg a step) | 3 | 0 | 2.9e-9 |
| 2 MW/m² | never | completed | 0.0126 kg | 0 | 0 | 2.1e-9 |
| 3 MW/m² | column | completed; last step's total mismatch 66 % of the throughput | 0.0194 kg (0, 0.0053, 0.0142 kg a step) | 3 | 75 | 5.5e-9 |
| 3 MW/m² | never | refused in step 1: a 4.65 mm corner move degrades the mesh even after a rebuild | — | 0 | — | — |

  Under the amended default (2026-10-12), the same sphere: at 3 MW/m² the steps complete with 47.4 g removed, 66
  repeats, no base faults and no buried cell above T_feed, energy residual 3.9e-9; at 6 MW/m² with 169.0 g removed
  (near the 190 g the heat could melt from 880 K), 154 repeats, 49 base faults and 289 buried cells above T_feed, the
  hottest 959 K, energy residual 6.9e-9. Each 0.5 s step ran as several parts (5 to 7 at 6 MW/m²).

- [ ] **Step 3: The cost of a step on the 100 mm sphere.** Create `analysis/skin_timing.py`:

{{TIMING}}

  Run `caffeinate -i "$PY" analysis/skin_timing.py ../measure/timing_100mm.json` (about 15 min). Measured in the
  build (9 417 columns, 10 µm cells, 0.5 s macro steps, the per-column rule; planning in brackets):

The table gives the wall time of each 0.5 s macro step: 98–99 % of it is in the skins' passes, and every step repeated
once. No rebuild occurred; the energy residual was 2.5e-11 (solid) and −1.4e-10 (melting); 8.3 g was removed while
melting.

| Phase | Cells | Per macro step | 3D solve | Moving the mesh |
|---|---|---|---|---|
| below melting (300 K, 1 MW/m²) | 374 708 | 62, 70 and 66 s [55–57] | 0.54 s [0.45] | 0.58 s [0.50] |
| melting (880 K, 2 MW/m²) | 456 769 [449 911] | 82, 115 and 133 s, rising as more columns melt [64, 110, 131] | 0.71 s [0.77] | 0.67 s [0.65] |

  Under the amended default (2026-10-12, run alone): below melting 82, 93 and 60 s a macro step (8 passes in all over
  the three steps); melting 123, 307 and 1 041 s, rising
  steeply as more columns melt (48 repeats in all; 16.6 g removed, against 8.3 g under the one-repeat rule; energy
  residual −7.5e-9). The skins took 1 449 s of the 1 470 s and the 3D solves 15 s. At the third step's rate a 600 s
  flight of 1 200 macro steps would take about two weeks, and more as melting spreads. With option 1 (parts of at
  most 0.0625 s while melting, the heating interpolated in time; d1b5636, run alone): below melting 80, 90 and 59 s a
  step (unchanged: the cap is off), melting 122, 121 and 226 s, energy residual −9.9e-9.

- [ ] **Step 4: Level 3, the Step 2 check against SESAM with skins** (spec §11.1). At these costs each flight takes
  hours; run them in the build, not before: the Step 2 verification cases of `analysis/reentry_model_thermal_verification.py`
  with `--surface-model skin` (non-melting, so no `--melt`), against Step 2's thresholds in
  `tests/test_reentry_model_reference_thermal.py`; record the differences from the element model's results.
  Held in the build because the coupling was wrong at 0.5 s; since the amendment the coupling converges and level 2
  passes at 0.5 s. Asha chose to spend the cost on 2026-10-16: the 100 mm flight under SESAM-equivalent heating alone on
  the laptop, the other three on a second machine (`docs/superpowers/plans/2026-10-16-skin-level3-flights.md`).

  **Files:** Modify `analysis/reentry_model_thermal_verification.py`; Create `tests/test_reentry_model_thermal_verification.py`.
  The script gains `--surface-model skin`: each skin run is labelled `<case>__<mode>__skin`, so it shares the output
  directory with the element run of the same case and mode, and `summary_skin.md` puts the two side by side; every run
  is judged against Step 2's thresholds (a test keeps the script's copy equal to the reference test's), reports the
  skins' totals and its wall time a macro step, and appends a line of progress every 20 macro steps to
  `<label>.progress`; `--t-max` ends the flights early as a smoke check, kept out of the summaries, and
  `--summary-only` tabulates the runs gathered in the output directory from several machines.

{{DIFF_THERMAL_VERIFICATION}}

{{TEST_THERMAL_VERIFICATION}}

  Run `"$PY" -m pytest tests/test_reentry_model_thermal_verification.py -q` (6 passed), then the smoke check
  `"$PY" analysis/reentry_model_thermal_verification.py --outdir <scratch> --cases d050 --modes sesam --surface-model skin --t-max 2`.
  Measured on the laptop: 4 macro steps on 2 367 columns in 60 s, 15 s a step at 115 km, the books to 3·10⁻¹², peak
  memory 1.1 GB. Then the flights, each alone or on a machine with cores to spare:
  `caffeinate -i "$PY" analysis/reentry_model_thermal_verification.py --cases d100 --modes sesam --surface-model skin`.
  **Results: not yet measured.** The 100 mm reference reaches the ground at 366 s (732 macro steps), the 50 mm one at
  559 s (1 118).

```bash
git add analysis/reentry_model_thermal_verification.py tests/test_reentry_model_thermal_verification.py
git commit -m "Skins: level 3 flies the Step 2 verification with a skin on every corner, against Step 2's thresholds"
```

- [ ] **Step 5: Record and commit.** Add the measurements of Tasks 0, 9 and 10 to `00-shared-context.md` as facts, and
  write the Results section of this sub-plan. (Done: facts 115–120 for the build, 121–125 for the amendment, and
  Results, below.)

```bash
git add analysis/skin_timing.py
git commit -m "Skins: the interface repeat's triggers and the cost of a step measured"
```

## Notes for sub-plans 21 and 22

- Sub-plan 21 adds the cells' velocity, the surface gradient and curvature operators (ruling 7), and level 2 with flow
  on — on the plate, not the box.
- The interface is damped repetition with parts (Asha, 2026-10-12). Flow will make each pass dearer, and the skins'
  passes are already almost all of a step's time, so the pass count is the first thing to watch in sub-plan 21.
- Sub-plan 22's retirement of stage 1 is now the deletion of `skin1d.py`, its tests and `MeltingBody`'s `skin1d` path.
- The moves' remap still closes its books by one uniform offset on the moved nodes: 442 J over the plate's 42 moves,
  0.8 % of its absorbed heat, one-signed. It is reported every step; make it local if a flight shows it matters.

## Results (built and amended 2026-10-12)

**What was built.** Tasks 0–10, subagent-driven, in `prototype/work-2026-10-12-skin-heat/code/` on branch `skin-heat`
(70ecbdd..f5d6633, 14 commits); every task was reviewed for spec compliance and quality, and the whole branch once more
at the end. Then, the same day, a diagnosis of the interface at 0.5 s macro steps and Asha's three decisions (a45f204):
damped repetition as the default, the metal crossing the base valued along the 3D profile, and the macro step run as
parts; then the fault counter redefined (fb4fe25), the melting part cap with the heating in time (d1b5636), and on
2026-10-14 the absorption rule with the hot plate as level 2b, then the skins' swelling setting: levels 2 and 2b at the
solid's density and the swelling checked at level 1, then the fine band under the surface (option 1), and on 2026-10-16
the skins' 0.8 mm minimum with a 0.4 mm fine zone (option B′) and the part cap's look-ahead, and level 3's option in
the Step 2 verification script. The unit tier went from 341 to 403 passed (1 skipped),
the reference tier gained two tests, the FEniCSx
tier from 17 to 19, and a melting element-model run is bit for bit the base's (Task 8 Step 6). Facts 115–120 record the
build's measurements, 121–146 the amendments'.

**What is verified, and at what step.** Level 1's heat checks (the cells, the conduction, the sideways heat, the base
law, the books, the profile valuation, the parts, the melting cap and its look-ahead, the heating in time, the
swelling, the 0.8 mm minimum with its 0.4 mm fine zone) pass as unit tests. Level 2, the plate (Task 9), with the skins
at the solid's density as the references hold their metal, the fine band under the surface and the skins' 0.8 mm
minimum: the cold plate passes against both references at 0.5 s macro steps (3.2 % ahead, 0.4 % ahead and 0.9 % late
against the one-dimensional reference) and at 0.05 s (3.6 % ahead, 0.0 and 0.9 % late); level 2b, the hot plate,
passes at 0.5 s (2.8 and 1.1 % late at 1.5 and 2.0 s, the skins' flowing depth within 0.9 % of the exact non-rigid
depth) and at 0.05 s (1.5 and 0.8 % late). With the 0.4 mm minimum the cold plate was 2.8 % late at 2.0 s at 0.05 s,
and on 2 mm elements throughout it failed at 0.05 s (6.2, 6.2, 6.4 % late). With the skins swelling, as measured before 2026-10-14, both plates passed at both steps, the
swelling's head start cancelling the delay (fact 133). Earlier measurements: the production cells passed at 0.05 s and
the two cell sizes differed by 1.4 % (1 % asked); the 0.1 and 0.25 s excesses at 1.0 s (+8.0 and +12.7 %) were cured
by the melting part cap (fact 132). Mass closes to 1.4e-15 and energy
to 4e-9 throughout (to 2e-11 since the 3D solver's linear tolerance is tightened in skin runs). Level 3 against SESAM was not run (Task 10 Step 4).

### The final review's findings, and what became of them

1. **At 0.5 s macro steps the interface coupling was wrong — cured by damped repetition.** Found by the final review
   and reproduced: on the plate at 0.5 s, 53 % of the absorbed heat (29.9 of 56.2 kJ) stayed in the deferred mismatch
   by 2 s and nothing melted, while the books closed. The diagnosis (fact 121): once the skins near melting, the heat a
   skin sends down falls off a cliff a little above the start-of-step base temperature — at the plate's centre from
   +5.0 to −7.1 MW/m² between 30 and 60 K above it, as a base held above the melting band melts the skin, draws metal
   at T_rigid and melts that too. The 3D surface's admittance (about 29 000 W/m²K) is four times the law's slope there,
   so each repeat jumped the cliff and back (828 and 927 K, alternately). Damped repetition (Aitken's relaxation)
   converges every step at every macro step tried, with at most 0.44 % of a step's heat deferred.
2. **Large steps then melted too much — cured by the profile valuation and the parts.** With the coupling converged,
   the plate at 0.5 s was 97 % ahead at 1.0 s: the draws and the move valued 1.55 mm of swept metal at the base node's
   enthalpy, about 99 000 J/kg above the 3D profile's average, which predicts 0.72 mm of extra melting against 0.76 mm
   observed (0.46 against 0.45 mm at 0.25 s; fact 122). And single steps asked for moves of 1.3–5.2 mm, which the mesh
   refused. Valued along the profile and run as parts no longer than a quarter-edge move, the plate at 0.5 s is within
   5 % of both references (fact 123).
3. **The base reconciliation buried superheated cells — reduced, still open.** Under 6 MW/m² on 4 mm elements, 289
   buried cells end above T_feed, the hottest at 959 K (2 910 K before); none at 3 MW/m² (fact 124). What the move's
   sweep and the skins' draw or hand-down book differently still goes into one bottom cell (spec §12).
4. **The repeat rule** is decided: damped repetition (Asha).

### Still open after the amendment

- **The cold plate's remaining late slowness** (facts 141–145): 1.7 % slow from 1.0 to 2.0 s with the 0.8 mm minimum,
  about 1 % of it the skins' own cells and sub-steps; the rebuilds' transfer still drops the base 1.3–1.5 K on the fine
  band. The hot plate's onset at 0.5 s steps is 9.0 % late at 1.0 s (not judged): its top crosses T_rigid inside the
  first macro step, before the look-ahead has a rise to read -- a body that starts near its melting range.

- **The fine band's cost** (facts 139–140): it resolved the delay on 2 mm elements (both plates within 3.3 % at 1.5
  and 2.0 s), but a melting step on the 100 mm sphere costs 45 % more (1 010 against 695 s for the strong-melting case),
  the extra in the 3D solves, the moves' remap of the larger mesh and the rebuild; the plates cost 50–80 % more.
  Accepted, with the band in place of prism layers, until level 3 times real flights (Asha, 2026-10-15).
- **Thin elements at the plate's 90-degree edges** (fact 136): the fill misses its request at the box's edges, and on
  0.5 mm elements the corrections exceed the move allowance and grow until a skin runs empty; layers under the faces
  are crushed by the fill unless blended to uniform near the sides. The sphere has no such edges.
- **Coarse surface meshes under strong melting** (fact 131): the 4 mm test sphere is refused at the rim of a melting
  region at 6 MW/m², where a corner absorbs a 0.7 mm lens its neighbours do not have; not supported. The 2 mm
  production mesh follows the same case. A rim limit on absorption made the 4 mm sphere pass and is the fallback.
- **The cost** (facts 125, 126 and 129): the melting 0.5 s steps on the 100 mm sphere took 123, 307 and 1 041 s, 53 %
  of it on parts rejected for their moves; with parts of at most 0.0625 s while melting they take 122, 121 and 226 s
  (82, 115 and 133 s under the one-repeat rule, which left heat parked).
- **Agreement of the backends across a rebuild**: to round-off only up to the first rebuild, because the mesher answers
  round-off in the corner positions with a different mesh.
- **Level 3** against SESAM (Task 10 Step 4).

### Corrections made in the build, with the rulings that made them

- Task 0: `locate` retries the points it misplaced, with 256 candidates; stale comments on the own-sweep booking corrected.
- Tasks 2 and 3 were done together, because Task 2's tests call Task 3's `minimum_mass`.
- Task 4: metal drawn at the T_rigid cap carries no sensitivity to the base, so the base law keeps its slope.
- Task 6: both backends drop a nodal law on a pinned node and report its power in `Q_dropped`.
- Task 7: the FEniCSx skin-body test's tolerances were loosened to mass 1e-9 and removed mass 1e-8 relative, because
  round-off tips two columns across a threshold.
- Final review: an explicit `--prism-layers` is refused with skins (exit code 2); the moving mesh's docstrings describe
  own sweeps and the projection; `_absorb_depth` is pinned by a test that fails on a sign error; comments made true; the
  Task 10 report's two misworded figures were corrected before they became facts.

### Left deferred (each judged harmless at stage 2 by the final review)

- `SkinBody.advance` does not restore the skins if a pass raises (no `try`/`finally`). This is harmless while every
  error ends the run, but it must be fixed before anything catches the error and continues, such as sub-plan 22's
  retries.
- The rebuild's L2 projection is not monotone: measured 0.9 K above and 4.8 K below the old range.
- `--neutral` silently overrides the smoothing flags in the recession harness.
- Mixed-move remainder sharing in `swept_volumes` has no direct test.
- `columns.py` has gaps: the clip and the non-manifold path are untested, there is no guard for a degenerate triangle,
  and `self.points` aliases a float input.
- `skin.py` has gaps:
  - the S branches, the property tables' accuracy claims and a liquid fill are untested;
  - `take_top` has no bound check;
  - `exchange_bottom` mutates the state before it raises.
- The trial-pass test has an assertion that cannot fail.
- The sensitivity's right-hand side uses the interpolated `cp_eff`, not the table's slope.
- The sideways step can overshoot at cell level by 0.03 K where a bin straddles two cells.
- The nodal law: duplicate nodes accumulate, its FEniCSx pinned path has no test of its own, and Dirichlet nodes are
  counted in `Q_int`.
- `absorbed_mass` records the requested mass, not the swept mass (about 2 % apart).
- `Q_dropped` is unchecked.
- `wall["rebuild"]` misses a second rebuild within a step.
- The CLI: the JSON settings for element runs changed shape, skins with `--melt off` or `--dt-continuum` are untested,
  and the skin flags are not validated on element runs.
- `skin_plate.py` prints but does not judge the criterion, and the reference JSON carries no code hash.
