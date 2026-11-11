# Sub-plan 19 — the moving 3D mesh (stage 1 of the layered-skin design) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first (required reading for every Step 3 sub-plan), then the spec.

**Goal:** Make the 3D model's mesh move with its surface — corners moved by a prescribed displacement, a deforming zone
beneath them, the enthalpy remapped with exact books, and the mesh rebuilt inside the same surface before any element
degrades — on both thermal backends, with no skins yet.

**Architecture:** A new module `reentry_model/movemesh.py` holds the geometry helpers, the deforming zone (`ZoneMap`),
and `MovingMesh`, which owns a `VolumeMesh`, its thermal solver and its material: `move(u_corners)` moves the corners and
the zone, carries the nodal enthalpy to the moved nodes by its gradient and closes the books with one uniform offset,
reporting the volume and enthalpy each corner swept out of the 3D model; `rebuild()` regenerates the tetrahedra with
gmsh inside the unchanged surface triangulation and carries the enthalpy over; `move_with_policy` decides when.
Each thermal backend gains `update_geometry(elements=None)`. Nothing in the model's run path uses any of this until
stage 2 (sub-plan 20): stage 1 is verified on its own, with prescribed motion.

**Tech Stack:** Python 3.12 (`drama_env`), NumPy, SciPy, gmsh 4.15.2, scikit-fem backend, FEniCSx 0.11 backend
(`fenicsx_env`), pytest.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` (committed cc15b3f), §5.1 (the corners'
tangential smoothing), §5.3 (the moving 3D mesh) and §11 (level 1's moving-mesh checks, build stage 1). Sub-plans 20
(skins, heat only), 21 (flow) and 22 (spraying and the flight) build on this one and are written in full after it is
built and measured, as the spec gates each stage on the one before.

**How this plan's code was checked, and built.** Every code and diff block below was first tested in a copy of
`prototype/work-2026-10-08-skin/code/` at its head 49d88e0 (2026-10-11), then built subagent-driven in
`prototype/work-2026-10-12-moving-mesh/code/` (2026-10-12; task reviews, a final whole-branch review and their fixes are
in the Results section at the end). After the build the tested copy, `prototype/work-2026-10-11-moving-mesh-plan/code/`,
was made identical to the built one, so the blocks below are the code as built and reviewed. This document is generated
from that copy by `prototype/work-2026-10-11-moving-mesh-plan/make_plan19.py`, which checks that the blocks, applied in
order to the base, reproduce the files exactly. The new module has 31 tests and FEniCSx three more; the unit tier goes
from 303 to 334 passed with 1 skipped and the same five known failures (two failed and three errors, all needing the
melting SESAM references of Task 11 of the original plan); the FEniCSx tier goes from 14 to 17 passed. Measurements
quoted below were each made once, on the cases named; they are recorded in Task 8 as facts.

## Global Constraints

- Work in a fresh work copy, `prototype/work-2026-10-12-moving-mesh/code/`, a git clone of
  `prototype/work-2026-10-08-skin/code/` checked out at 49d88e0 (Task 0). It is a git repository of its own, ignored by
  the outer one; commit there after every task. Never edit `prototype/proto3/`, the package `reentry_model/` at the repo
  root, the specs or the monolithic plan.
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the work copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Baseline at 49d88e0: 303 passed, 1 skipped, and five known failures (`test_reentry_model_compare.py`'s melt plots and
  the melting reference fixtures — two failed, three errors). Record the counts after every task. No other failure is
  acceptable.
- FEniCSx tier: `FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`. Baseline 14 passed.
- Nothing in the existing run path changes behaviour: no CLI flag, no history column, no default. A melting or
  non-melting flight from the work copy must be bit-identical to the base's (Task 8 Step 1 checks one).
- Books: in every test that moves or rebuilds, total enthalpy plus what the corners swept equals the enthalpy before, to
  relative 1e-12; the rebuilt volume equals the old volume to relative 1e-12.
- gmsh runs on one thread in a rebuild: fenicsx_env's gmsh 4.15.2 segfaults meshing a discrete surface's volume with
  HXT on four threads (found while planning, 2026-10-11).
- Exit codes: a move no rebuild can absorb raises `RuntimeError` (exit code 1 when it reaches the CLI in stage 2); a
  surface the validity gate rejects raises `movemesh.RemeshError`, a `RuntimeError`.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).

## Review Focus

1. **Corners moving outward** — rigid metal joining the 3D model, as freezing will in stage 2: the swept volume is
   negative, the books still close and a uniform field stays uniform. Pinned by
   `test_an_outward_move_books_metal_coming_in` (Task 4).
2. **Corners moving along the surface** — the tangential smoothing every move adds: it sweeps no volume and keeps the
   books. Pinned by `test_a_tangential_move_sweeps_nothing_and_keeps_the_books` (Task 4).
3. **A solver set up again on a rebuilt mesh** — the scikit-fem backend kept its AMG hierarchy across `setup` and
   crashed on the first solve that reused it ("matrix and preconditioner have different shapes", found while
   measuring); FEniCSx must keep exact books across a rebuild too. Pinned by
   `test_a_solver_set_up_on_a_new_mesh_drops_its_old_amg_hierarchy` (Task 2) and
   `test_a_rebuild_keeps_the_energy_on_fenicsx` (Task 5).
4. **Corner identity across repeated rebuilds** — stage 2's columns live on the corners in a fixed order, and gmsh
   renumbers nodes. Pinned by `test_corners_keep_their_order_through_repeated_rebuilds` (Task 5).
5. **A move no rebuild can absorb** — a fold in the surface sharper than one triangle (found while measuring, Task 8):
   a clear model failure, not a hang or an inverted element. Pinned by
   `test_a_move_refused_on_a_fresh_mesh_is_a_model_failure` (Task 6).

---

### Task 0: The work copy

**Files:** none in the code; one line in the outer repository's `.gitignore`.

- [ ] **Step 1: Clone the stage-1 skin work copy at its reviewed head**

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
git clone work-2026-10-08-skin/code work-2026-10-12-moving-mesh/code
git -C work-2026-10-12-moving-mesh/code checkout -b moving-mesh 49d88e0
```

- [ ] **Step 2: Ignore it in the outer repository** — append to `.gitignore`, under the existing line for
  `prototype/work-2026-10-08-skin/code/`:

```
prototype/work-2026-10-12-moving-mesh/code/
```

- [ ] **Step 3: Record the baseline.** Run the unit tier and the FEniCSx tier in the work copy (commands in Global
  Constraints). Expected: 303 passed, 1 skipped, two failed and three errors (the five known); FEniCSx 14 passed.

### Task 1: Element and surface geometry

**Files:**
- Create: `reentry_model/movemesh.py` (its header and the geometry section)
- Create: `tests/test_reentry_model_movemesh.py` (its header, helpers and the Task 1 tests)

**Interfaces:**
- Produces: `RemeshError(RuntimeError)`; `signed_volumes(points, tets) -> (ne,)`; `radius_ratio(points, tets) -> (ne,)`
  (3 r_in / R_circ, 1 regular); `triangle_quality(points, faces) -> (nf,)` (2 r_in / R_circ); `max_dihedral(points,
  faces) -> float` degrees; `element_gradients(points, tets) -> ((ne, 4, 3), (ne,))`; `nodal_gradient(points, tets, f,
  elements=None) -> (n, 3)`; `lumped_masses(points, tets, rho) -> (n,)`; `orient_outward(points, faces, normals) ->
  faces`; `corner_frame(points, faces, corners) -> (normals (nc, 3), dual_areas (nc,))`; `surface_edges(faces) -> (ne,
  2)`; `swept_volumes(old, new, faces, faces_corner, n_corners) -> (nc,)` (each surface triangle's prism between its old
  and new positions, a third to each corner, positive inward; added after the Task 4 review); `tangential_smoothing(points, faces, normal, lam=SMOOTH_LAMBDA, clamp=SMOOTH_CLAMP, edges=None) -> (n, 3)`;
  constants `DEFORM_DEPTH = 15e-3`, `RATIO_MIN = 0.05`, `MIN_GAP = 20`, `ZONE_NEIGHBOURS = 12`, `SMOOTH_LAMBDA = 0.5`,
  `SMOOTH_CLAMP = 0.35`.

- [ ] **Step 1: Write the test file's header, its helpers and the Task 1 tests** (`tests/test_reentry_model_movemesh.py`).
  The helper `inward` uses `MovingMesh.frame` from Task 4; it is defined now so that the file grows only by appending.

{{TEST_T1}}

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`
Expected: collection error, `ImportError: cannot import name 'movemesh'`.

- [ ] **Step 3: Create `reentry_model/movemesh.py` with its header and the geometry section**

{{MM_T1}}

- [ ] **Step 4: Run them to see them pass**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`
Expected: 6 passed.

- [ ] **Step 5: Unit tier and commit.** Expected 309 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: element and surface geometry helpers"
```

### Task 2: The backends on moved geometry

**Files:**
- Modify: `reentry_model/thermal/skfem_backend.py` (`setup`, new `update_geometry`)
- Modify: `reentry_model/thermal/fenicsx_backend.py` (`setup`, new `update_geometry`)
- Test: `tests/test_reentry_model_movemesh.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Produces: `solver.update_geometry(elements=None)` on both backends — the mesh's `points` have moved in place;
  scikit-fem recomputes the element matrices of `elements` (all if None) and the surface; FEniCSx copies the coordinates
  into dolfinx's geometry and refreshes its cell volumes and the surface (it recomputes every Jacobian at assembly
  anyway). `setup` on a new mesh now discards the scikit-fem AMG hierarchy. (After the Task 2 review: the scikit-fem
  test updates only the elements touching the moved outer shell, a real subset, and FEniCSx refreshes `cell_volumes`.)

- [ ] **Step 1: Append the Task 2 tests** to `tests/test_reentry_model_movemesh.py`:

{{TEST_T2}}

  and to the end of `tests/test_reentry_model_fenicsx.py`:

{{FX_T2}}

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q -k "fresh_solver or amg"`
Expected: `AttributeError: 'SkfemThermalSolver' object has no attribute 'update_geometry'`, and the AMG test fails with
`ValueError: matrix and preconditioner have different shapes`.

- [ ] **Step 3: Apply the scikit-fem change**

{{DIFF_SKFEM}}

- [ ] **Step 4: Apply the FEniCSx change**

{{DIFF_FENICSX}}

- [ ] **Step 5: Run them to see them pass**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q` — expected 8 passed — and
`FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q -k moved_geometry` — expected 1 passed.

- [ ] **Step 6: Both tiers and commit.** Expected unit tier 311 passed, 1 skipped, the five known; FEniCSx 15 passed.

```bash
git add reentry_model/thermal/skfem_backend.py reentry_model/thermal/fenicsx_backend.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_fenicsx.py
git commit -m "Thermal backends follow moved node coordinates, and a solver set up on a new mesh drops its old AMG hierarchy"
```

### Task 3: The deforming zone

**Files:**
- Modify: `reentry_model/movemesh.py` (append the zone section)
- Test: `tests/test_reentry_model_movemesh.py` (append)

**Interfaces:**
- Consumes: Task 1's `corner_frame`, `orient_outward` (tests).
- Produces: `ZoneMap(points, corners, depth=DEFORM_DEPTH, neighbours=ZONE_NEIGHBOURS)` with `.displacement(u_corners
  (nc, 3)) -> (n, 3)`, `.interior` (zone node ids that are not corners), `.corners`, `.depth`.

- [ ] **Step 1: Append the Task 3 tests**

{{TEST_T3}}

- [ ] **Step 2: Run them to see them fail** — Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q -k zone`.
  Expected: `AttributeError: module 'reentry_model.movemesh' has no attribute 'ZoneMap'`.

- [ ] **Step 3: Append the zone section to `reentry_model/movemesh.py`**

{{MM_T3}}

- [ ] **Step 4: Run them to see them pass** — expected 10 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 313 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: the deforming zone moves the nodes under the corners in proportion to their depth"
```

### Task 4: Moving the mesh — the remap and its books

**Files:**
- Modify: `reentry_model/movemesh.py` (append `MoveResult` and the class `MovingMesh` up to and including `move`)
- Test: `tests/test_reentry_model_movemesh.py` (append)

**Interfaces:**
- Consumes: Tasks 1–3; the solver's `temperature()`, `set_temperature()`, `energy()`, `update_geometry(elements)`,
  `emissivity`; the material's `enthalpy`, `temperature_from_enthalpy`, `rho`.
- Produces: `default_sizes(mesh) -> dict` (gmsh sizes for a rebuild from the generator's parameters; in this section
  because `MovingMesh.__init__` calls it); `MoveResult(refused, swept_volume (nc,), swept_enthalpy (nc,), correction,
  worst_ratio, rebuilt=False)`;
  `MovingMesh(mesh, solver, material, depth=DEFORM_DEPTH, ratio_min=RATIO_MIN, min_gap=MIN_GAP, sizes=None)` with
  `.mesh`, `.solver`, `.material`, `.corners` (current node id of each corner, fixed order), `.faces_corner` (surface
  triangles in corner numbering, outward), `.edges_corner`, `.zone`, `.compression` (inward motion of each corner since
  the last rebuild), `.steps_since_rebuild`, `.rebuilds`, `.last_rebuild`; methods `frame() -> (normals, dual_areas)`,
  `smoothing(normal=None, lam=None, clamp=None) -> (nc, 3)`, `surface_quality() -> float`, `surface_fold() -> float`
  degrees, `move(u_corners) -> MoveResult`. Stage 2 credits each column with `swept_volume` and `swept_enthalpy`.

- [ ] **Step 1: Append the Task 4 tests**

{{TEST_T4}}

- [ ] **Step 2: Run them to see them fail** — Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`.
  Expected: the Task 4 tests fail with `AttributeError: module 'reentry_model.movemesh' has no attribute 'MovingMesh'`.

- [ ] **Step 3: Append `MoveResult` and `MovingMesh` (through `move`) to `reentry_model/movemesh.py`**

{{MM_T4}}

- [ ] **Step 4: Run them to see them pass** — expected 17 passed in the file. The measured sizes behind the tolerances,
  each from one run: on the steep, kinked field of the books test over 4 mm cells the offset that closes the books is
  0.7 % of what the corners swept, every corner books its dual area times the 0.3 mm move to within 1.0 %, and the
  nodes below the zone do not change at all; on a linear field the offset is 0.0035 J against 3877 J swept and moves
  no node by more than 3.3e-6 K; a move along the surface books at most 0.6 % of each corner's dual area times its
  displacement. Each corner's swept volume is booked at its own enthalpy before the move: booking inward slabs at the
  mean of old and remapped enthalpy shrank the offset (from -16.2 J to +0.29 J per 0.3 mm move) but made Task 7's exact
  moving-front test nine times worse (0.0010 to 0.0093 of the rise), so the surface value stays (final review).

- [ ] **Step 5: Unit tier and commit.** Expected 320 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: corners and zone move, the enthalpy is remapped by its gradient and the books close exactly"
```

### Task 5: The rebuild

**Files:**
- Modify: `reentry_model/movemesh.py` (add the method `rebuild` to `MovingMesh`; append the rebuild section)
- Test: `tests/test_reentry_model_movemesh.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Consumes: Tasks 1–4; `mesh.VolumeMesh`; the solver's `setup(mesh, material, emissivity)`.
- Produces: `MovingMesh.rebuild(reason="requested") -> RebuildResult(residual, relative_residual, n_nodes, n_elements,
  seconds, reason)`; `MovingMesh.rebuild_residual_total` (the running sum of the pre-correction residuals); module
  functions `check_surface(points, faces) -> volume` (raises `RemeshError`), `check_rebuilt_boundary(mesh,
  expected_faces)` (raises `RemeshError` unless the rebuilt boundary is exactly the corner triangles), `rebuild_volume(corner_points, faces_corner, h_surface, h_core, band, ramp) -> (points, tets,
  node_of_corner)`, `locate(points, tets, query, candidates=16) -> (elements, barycentric (nq, 4))`. After a rebuild
  `moving.mesh` and `moving.solver.mesh` are the new mesh; `moving.corners` maps each corner to its new node.

- [ ] **Step 1: Append the Task 5 tests**

{{TEST_T5}}

  and to the end of `tests/test_reentry_model_fenicsx.py`:

{{FX_T5}}

- [ ] **Step 2: Run them to see them fail** — expected `AttributeError: 'MovingMesh' object has no attribute 'rebuild'`
  and `module 'reentry_model.movemesh' has no attribute 'check_surface'`.

- [ ] **Step 3: Add the method `rebuild` at the end of the class `MovingMesh`** (after `move`, indented as a method):

{{MM_T5A}}

- [ ] **Step 4: Append the rebuild section at the end of `reentry_model/movemesh.py`**

{{MM_T5B}}

- [ ] **Step 5: Run them to see them pass** — expected 27 passed in the file and the FEniCSx rebuild test passing.
  A rebuild of the coarse test sphere takes well under a second. The offset that makes a rebuild's total exact is spread
  over the new deforming zone only, where the interpolation error arises (final review: the residual is one-signed);
  on the kinked test field over 4–20 mm cells the pre-correction residual is 4.1 % of the total, on a linear field
  6e-8; gmsh errors reach the caller as `RemeshError`.

- [ ] **Step 6: Both tiers and commit.** Expected unit tier 330 passed, 1 skipped, the five known; FEniCSx 16 passed.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_fenicsx.py
git commit -m "Moving mesh: rebuild inside the unchanged surface triangulation, enthalpy carried over with exact totals"
```

### Task 6: The rebuild policy

**Files:**
- Modify: `reentry_model/movemesh.py` (add `needs_rebuild` and `move_with_policy` to `MovingMesh`, after `rebuild`)
- Test: `tests/test_reentry_model_movemesh.py` (append)

**Interfaces:**
- Produces: `MovingMesh.needs_rebuild() -> "compression" | None` (the zone pressed to half its depth, after at least
  `min_gap` moves); `MovingMesh.move_with_policy(u_corners) -> MoveResult` (rebuild first if compression calls for it;
  if the move is refused, rebuild with reason "quality" and move again; refused again raises `RuntimeError` containing
  "even after a rebuild"). Stage 2 calls `move_with_policy` once per macro step.

- [ ] **Step 1: Append the Task 6 tests**

{{TEST_T6}}

- [ ] **Step 2: Run them to see them fail** — expected `AttributeError: 'MovingMesh' object has no attribute
  'move_with_policy'`.

- [ ] **Step 3: Add the two methods to `MovingMesh`, after `rebuild` and before the `# -- rebuilding the mesh`
  section**

{{MM_T6}}

- [ ] **Step 4: Run them to see them pass** — expected 30 passed in the file. The minimum-gap test is built so that compression
  crosses half the zone two moves before the gap is met: the rebuild comes at call 9 with the gap and at call 7 without
  it (Task 6 review).

- [ ] **Step 5: Unit tier and commit.** Expected 333 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: rebuild when the zone is half compressed or a move would degrade an element"
```

### Task 7: A moving front against an exact solution

**Files:**
- Test: `tests/test_reentry_model_movemesh.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Consumes: Tasks 1–6 and each backend's `step(dt, q_conv, T_amb, dirichlet=...)`.
- Produces: `run_moving_front(name)` in the test module, which the FEniCSx test imports.

This is level 1's check of the moving mesh and the conduction together. The face x = 0 of a 20 mm bar of 0.1 mm cells
is held at 600 K and moves into the bar at 10 mm/s for 0.5 s; the exact temperature ahead of a front at a fixed
temperature moving at constant speed into a solid at T0 is (Carslaw and Jaeger) θ = ½ [erfc((ξ − v t)/(2√(α t)))
e^(−v ξ/α) + erfc((ξ + v t)/(2√(α t)))], with θ = (T − T0)/(Ts − T0), ξ the distance ahead of the front and α = 1e-5 m²/s.
It was checked while planning against a 4000-cell finite-difference solution in the front's frame: largest difference
3.6e-6. The bar's sides slide with the zone, so they stay on their planes.

- [ ] **Step 1: Append the Task 7 tests**

{{TEST_T7}}

  and to the end of `tests/test_reentry_model_fenicsx.py`:

{{FX_T7}}

- [ ] **Step 2: Run them** — `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q -k front` (expected 1 passed: the
  largest error was 0.0010 of the temperature rise while planning, against a tolerance of 0.005) and the FEniCSx tier
  (expected 17 passed; the two backends agree to 1e-6 K).

- [ ] **Step 3: Both tiers and commit.** Expected unit tier 334 passed, 1 skipped, the five known; FEniCSx 17 passed.

```bash
git add tests/test_reentry_model_movemesh.py tests/test_reentry_model_fenicsx.py
git commit -m "Moving mesh: a heated front moving into a bar follows the exact solution on both backends"
```

### Task 8: Measurement on the 100 mm sphere, and the record

**Files:**
- Create: `analysis/movemesh_recession.py`
- Modify (outside the work copy, documents only): `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`
  (facts), this sub-plan (a results section)

- [ ] **Step 1: The run path is unchanged.** Nothing in the model's run path calls the new code, so a flight from the
  work copy must be bit-identical to one from the base. Run the 100 mm non-melting flight to 30 s from each, with the
  same seed, and compare:

```bash
P="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
for copy in work-2026-10-08-skin work-2026-10-12-moving-mesh; do
  (cd "$P/$copy/code" && caffeinate -i "$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 \
     --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --heating physics --t-max 30 --seed 12345 \
     --outdir "$P/work-2026-10-12-moving-mesh/runs/bitid_$copy")
done
"$PY" - "$P/work-2026-10-12-moving-mesh/runs/bitid_work-2026-10-08-skin" "$P/work-2026-10-12-moving-mesh/runs/bitid_work-2026-10-12-moving-mesh" <<'PYEOF'
import glob, json, sys
import pandas as pd
a, b = (pd.read_csv(glob.glob(d + "/*.csv")[0]) for d in sys.argv[1:3])
assert list(a.columns) == list(b.columns) and a.equals(b), "history differs"
ra, rb = (json.load(open(glob.glob(d + "/*.json")[0]))["results"] for d in sys.argv[1:3])
differ = [k for k in ra if k != "runtime_s" and ra[k] != rb.get(k)]
assert not differ, differ
print("bit-identical:", len(a), "rows,", len(a.columns), "columns,", len(ra) - 1, "result fields")
PYEOF
```

  Expected: `bit-identical: 61 rows, 31 columns, ...` (checked while planning on the plan's tested copy). Only the
  `results` block is compared, without `runtime_s`; the rest of the JSON holds file paths and the git commit, which
  differ by construction.

- [ ] **Step 2: Create the measurement script**

{{HARNESS}}

- [ ] **Step 3: Run it with and without the tangential smoothing**, one after the other (not side by side: the run time
  is measured), under `caffeinate -i` on mains power:

```bash
"$PY" analysis/movemesh_recession.py --out ../measure/recession_smooth.json > ../measure/recession_smooth.log 2>&1
"$PY" analysis/movemesh_recession.py --no-smoothing --out ../measure/recession_nosmooth.json > ../measure/recession_nosmooth.log 2>&1
```

  While planning (one run each, 2026-10-11): with smoothing the nose receded 29.5 mm before a move failed even after a
  rebuild, 20 rebuilds by the counter (19 recorded; the 20th was tried inside the failed move), the first at 10.5 mm by
  compression and every later one by quality — almost every macro step from 21 mm; without smoothing 25.0 mm and 15 rebuilds. A rebuild took 1.6–2.3 s on one gmsh thread (130 000–177 000
  tetrahedra), against the spec's estimate of about 42 s; its pre-correction residual was 0.07–0.27 % of the total
  enthalpy, and a move's offset at most 1.4 % of what the corners swept. The worst surface triangle stayed at 0.71 with
  smoothing (0.36 without). What ends the run is the shape: the largest angle between neighbouring surface triangles grew
  from 19.5 degrees at 10 mm and 51.9 at 20 mm to 95 degrees, as recession in proportion to how directly a corner faces
  the flow sharpened the shoulder; an element whose four nodes lie on a fold sharper than one triangle is flattened by
  its own corners' motion, and a rebuild inside the same surface cannot undo it. The moment of failure depends on the
  meshes gmsh happens to make: with four gmsh threads the same runs failed at 27.0 and 20.5 mm. A macro step's mesh work
  took a median 0.45 s with smoothing and 0.27 s without, against 0.19–0.20 s for the conduction step.

- [ ] **Step 4: Record the facts** in `00-shared-context.md` (the next numbers after 101): the rebuild cost and residual,
  the policy's rebuild record, the smoothing's effect, the fold at failure and the dependence on gmsh's meshes, the mesh
  work against the conduction step, the AMG fix, and the one-thread gmsh setting with its reason. Mark each measured
  once.

- [ ] **Step 5: Add a results section to this sub-plan** with the counts of every tier after each task, the
  measurements, and anything that differed from the numbers above. The fold limit is recorded, not resolved, here:
  Asha decided on 2026-10-11 to measure it under heating-driven recession in stage 2 (sub-plan 20) before building any
  surface refinement.

## Results (built 2026-10-12)

Built subagent-driven in `prototype/work-2026-10-12-moving-mesh/code/` (branch `moving-mesh`, 15 commits from 49d88e0
to a18b531): one implementer per task on the cheapest model, a reviewer per task, a final whole-branch review on the most
capable model, one fix wave and its re-review. The tested copy beside this plan was then made identical to the build.

**Counts after each task** (movemesh test file / unit tier / FEniCSx tier): baseline – / 303 / 14; Task 1 6 / 309 / 14;
Task 2 8 / 311 / 15; Task 3 10 / 313 / 15; Task 4 17 / 320 / 15; Task 5 27 / 330 / 16; Task 6 30 / 333 / 16; Task 7
31 / 334 / 17; Task 8 and the final fixes 31 / 334 / 17. Every unit tier also had 1 skipped and the five known failures.
The 100 mm flight to 30 s from the build is bit-identical to the base's (61 rows, 31 columns, 28 result fields).

**What review changed against the code first planned.**
- Task 2: the scikit-fem geometry test now updates a real subset of elements; FEniCSx refreshes `cell_volumes` on a move.
- Task 4: each corner's swept volume is the exact prism its faces swept (`swept_volumes`), not the dual area times the
  normal displacement scaled to the total, which blew up wherever the first-order volumes nearly cancelled (a move along
  the surface gave a scale factor of about -7e14). New tests check the per-corner values. A first attempt that spread
  the remainder by dual area times displacement broke Task 7's exact solution and was dropped before review.
- Task 4/5: `default_sizes` sits beside `MovingMesh`, which calls it (the plan had put it in Task 5).
- Task 6: the minimum-gap test now fails if the gap is not enforced.
- Final review: the rebuild's offset is spread over the new deforming zone only and summed in
  `rebuild_residual_total`; `check_rebuilt_boundary` checks the rebuilt boundary is exactly the corner triangles; gmsh
  errors reach the caller as `RemeshError`; the rebuild-residual assertions are two-sided. Booking inward slabs at their
  mean enthalpy was tried and reverted (it shrank the books' offset but made the exact moving-front test nine times
  worse, 0.0010 to 0.0093 of the rise).

**Measurements (Task 8, one run each, both on mains power).** As planned: with the tangential smoothing the nose
receded 29.5 mm before a move failed even after a rebuild; 20 rebuilds by the counter (19 recorded), the first at
10.5 mm by compression, the rest by quality; a rebuild took 1.63–2.06 s, with a pre-correction residual of 0.07–0.27 %
of the total enthalpy; a move's offset at most 1.3 % of what the corners swept; worst element radius ratio 0.054, worst
surface triangle 0.71; the fold 19.5 degrees at 10 mm, 51.9 at 20 mm and 95.4 at most; a macro step's mesh work a median
0.43 s against 0.19 s of conduction. Without the smoothing: 25.0 mm, 15 rebuilds by the counter, worst surface triangle
0.36, fold 22.3 and 73.7 degrees at 10 and 20 mm and 97.8 at most, mesh work 0.26 s. Moving-front error against the
exact solution 0.00103 of the rise; the two backends agree to 4.3e-11 K.

**Left for sub-plan 20** (from the final review; the build's ledger, kept as
`prototype/work-2026-10-11-moving-mesh-plan/sdd-ledger.md`, holds every ruling and every deferred minor finding):
- The tangential smoothing sweeps second-order volume of either sign at each corner (one smoothing-only move on the
  100 mm sphere: up to about 1 µm of equivalent depth either way, 4 856 of 9 417 corners negative), which an empty skin
  cannot supply. Decided 2026-10-12 (Asha): the smoothing is made volume-neutral per corner (spec §5.1), built as
  sub-plan 20's Task 0.
- The enthalpy of metal coming in on an outward move is booked at the corner's enthalpy; stage 2 will supply the skin's.
- The rebuild residual is one-signed (99–186 J per rebuild in the run, 2.3 kJ in all); record it against each step's
  heat input in stage 2's history and set the threshold for an exact (intersection-based) transfer in those terms.
- Integration: `ThermalBody.i_centre` and the centre-node temperature after a rebuild (no centre node is embedded);
  `params["path"]` still names the original mesh file; `MovingMesh` assumes every element is active.
- The mesh work per step (0.43 s) is about twice the conduction step: optimise only if stage 2's run time asks for it.
