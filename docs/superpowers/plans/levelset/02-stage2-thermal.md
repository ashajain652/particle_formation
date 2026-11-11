# Level-set stage 2 — the thermal solve on cut cells — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve heat conduction in the metal {φ < 0} on stage 1's adaptive mesh, with the boundary conditions on the
zero contour, a moving surface and mesh rebuilds, so that ψ (the 895.1 K isotherm) can be built from a computed
temperature — as one tested module `reentry_model/levelset/cutthermal.py`, with prescribed surface motion and no skins
or flow.

**Architecture:** P1 temperature at the nodes of the active cells (cells holding any metal), backward Euler, Newton on
the radiation term with the secant heat capacity of each cell, so the energy books close to the Newton tolerance with
latent heat. Each cut cell's metal part is split into sub-tetrahedra for its exact mass integrals; a ghost penalty on
the faces of cut cells keeps sliver cells conditioned. Heat enters through the patches as q0 + q1·T − rad·εσ(T⁴ − T∞⁴),
which covers today's bare-surface heating and radiation and the skins' base law q_b = a + b·T_b. A moving surface keeps
the nodal temperatures and books the enthalpy of the metal it removes; a mesh rebuild projects the temperature
conservatively and closes the stored enthalpy with one uniform offset, reporting the residual before it.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy (sparse, sparse.linalg), pyamg (already used by
`reentry_model/thermal/skfem_backend.py`); pytest. No new dependency.

**Spec:** `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`: §7 (the thermal solve:
cut cells §7.2, ghost penalty §7.3, boundary conditions §7.5, backend §7.6), §20.11 item 2 (stage 2's verification),
§20.4 (the fine margin d), §20.14 item 1 (stages 1 and 2 before the shared-package refactor).

**Prototype:** every function and test below was written and run in `prototype/levelset/stage2/` on 2026-11-11
(module `lsproto2/cutthermal.py`, test materials `lsproto2/simplemat.py`, tests `tests/test_cutthermal.py`: 12 pass in
about 15 s). Each task ports its part of the prototype to `reentry_model/levelset/`, changing only the import prefix
`lsproto2.` → `reentry_model.levelset.`. The plate comparisons (facts 6–9) used read-only copies of the skyscrapers
build's Scheil material and one-dimensional reference (`prototype/levelset/stage2/ext/`, from build commit 3c16f4a) and
the scripts `scratch_plate.py` and `scratch_margin.py` there; they stay in the prototype (Decision B).

## Measured facts (prototype, 2026-11-11, laptop on mains power, shared with running flights; each measured once)

1. **Books.** With the patch fluxes, radiation and latent heat, the change in stored enthalpy equals the heat that came
   in to 3 × 10⁻¹³ relative with constant properties, and to the Newton tolerance with latent heat: about 4 × 10⁻⁸ on
   the Scheil alloy at a relative tolerance of 10⁻⁸ (6–8 Newton iterations a step), 8 × 10⁻⁸ on a linear melting
   range (up to 67 iterations on the first step, about 11 after).
2. **Accuracy.** A slab under a constant 2 MW/m² for 0.5 s, its sides cut off-grid: the surface temperature rise
   differs from the semi-infinite solution 2q/k·√(αt/π) by 1.38 % at 1 mm cells, 0.37 % at 0.5 mm and 0.04 % at
   0.25 mm — second order. The time step (0.01 against 0.0025 s) and lumped against consistent mass change it by
   under 0.2 %. Where a cell holds both the top and a side wall (a sharp corner) the patch tilts; on 1 mm cells the
   projected heat on a 4 mm wide slab was 3.3 % short. Rounded bodies have no such corners.
3. **Heat by projected area.** Applying a heat flux only to patches tilted less than about 25° from the flow (tried
   first) lost 41 % of a plate column's heat once its receding top crossed cells cut by its walls; applying
   q·max(n·d, 0) to every patch, the patch's share of the flux its projection catches, delivered 98–99 %.
4. **Slivers.** Without the ghost penalty a cell holding 10⁻¹⁰ to 10⁻⁴⁵ of its volume gives a diagonal 10⁻¹⁰ to
   10⁻⁴⁵ times the largest: the sliver node's temperature is undetermined, plain conjugate gradients fail in 5 000
   iterations and AMG converges only because its aggregation hides the node. With the penalty weight γ = 0.01 the
   ratio is about 6 × 10⁻³ and plain CG converges in about 330 iterations (γ = 0.1: 5 × 10⁻², 140).
5. **Jumps.** A surface receding 0.05 mm a step through 1 mm cells under 1 MW/m²: with γ = 0 or 0.01 the temperature
   4 mm below changes smoothly (no step's second difference above three times the median); with γ = 0.1 one jump of
   0.07 K, and with γ = 1 four of up to 0.24 K, each when the surface crossed a layer of nodes and the set of
   penalised faces changed at once. γ = 0.01 is the default.
6. **Plates, uniform 0.7 mm cells, 0.05 s steps** (a column of the level-2 plate's centre with insulated walls, liquid
   removed after every step where T ≥ T_feed = 910 K, against the independent one-dimensional reference at 2 µm
   cells and 1 ms steps): the cold plate (700 K, 6 MW/m²) recedes 0.800, 2.688 and 5.130 mm at 1.0, 1.5 and 2.0 s
   against 0.788, 2.712 and 5.216 mm (1.5 % ahead, 0.9 and 1.6 % behind); the hot plate (880 K, 2 MW/m²) 0.798, 1.806,
   2.904 mm against 0.788, 1.816, 2.916 (1.3 % ahead, 0.6 and 0.4 % behind), its non-rigid depth 1.514, 1.630,
   1.724 mm against 1.488, 1.570, 1.676 (1.7, 3.8 and 2.9 % thicker). Energy books 10⁻⁸. At 0.35 mm cells the cold
   plate's recession moves by at most 2.3 %. 8 s of wall time for 2 s of plate.
7. **Plates at 0.5 s steps fail.** The same columns at 0.5 s: cold 1.8 % ahead, 10.9 and 12.8 % behind; hot 24.6,
   15.1 and 11.9 % behind. The liquid stays on the plate through each long step and absorbs heat until the end-of-step
   removal; the reference removes it every millisecond (Decision A).
8. **Plates on the adaptive mesh with the fine margin** (2.8 mm base cells, 0.7 mm where ψ > −d and above, rebuilt
   after every step; walls on grid lines): with d = 4 mm the cold plate recedes 3.6, 4.0 and 4.2 % behind the
   reference at 1.0, 1.5 and 2.0 s (passes the skyscrapers' 5 %), the hot plate 9.6 % behind at 1.0 s (not judged),
   4.3 and 2.6 % at 1.5 and 2.0 s, its non-rigid depth 3.8–4.3 % thicker (passes 10 %). Each rebuild's stored
   enthalpy before the offset was within 3 × 10⁻⁶ of the old (cold) and 3 × 10⁻⁹ (hot).
9. **Without the margin** (d = 0, 0.7 mm only above ψ): the cold plate 7.4 % behind at 1.0 s (fails), the hot plate's
   non-rigid depth up to 6.6 % thick, and each rebuild's enthalpy residual up to 4.6 × 10⁻⁴ (cold) and 2.4 × 10⁻³
   (hot) — the coarse cells rebuilt just under the front carry the error the offset then spreads over the body. The
   margin earns its place, as §20.4 expected. Two earlier runs of this test were wrong for reasons of the test, not the
   solver: walls through coarse cells, and φ rebuilt as a plane at the patches' mean height, each moved the metal at
   every rebuild (residuals 1–3 %); φ is now carried by interpolation, exact under refinement.
10. **Cost.** The 50 mm sphere (fine 0.7 mm over a 4.6 mm shell, 1 055 790 cells, 659 028 active, 121 191 active
    nodes): setup 8.6 s, one 0.5 s step with radiation 6.6 s (5 Newton iterations, the AMG hierarchy rebuilt each
    iteration), peak 1.79 GB. The 100 mm sphere is *estimated* at about four times that (about 30 s a step), against
    about 330 s of geometry a step (stage 1, fact 12); not measured.

## Decisions (rulings in this plan, for Asha to check)

- **A. The plates are judged at 0.05 s steps in stage 2.** Stage 2 removes the liquid once a step (the skins, which
  take it continuously, are stage 3); at 0.5 s that lags by up to 25 % (fact 7). The 0.5 s judgement moves to stage 3,
  where the skins hold and release the liquid as the skyscrapers' do. Cost if wrong: an in-step removal rule in stage
  2 that stage 3 would replace.
- **B. The plate comparisons stay in the prototype.** The Scheil material and the one-dimensional reference exist only
  in the skyscrapers build (this branch's `reentry_model/material.py` has no melting), and the shared-package
  refactor that brings them (stage 0 + S) waits for the skyscrapers' port (§20.14 item 1). The package's tests use
  test-only materials (constant properties, a linear melting range); the plate numbers are recorded in the design.
  Cost if wrong: copying shared physics into the level-set package, which §4 forbids.
- **C. Heat carried by a flow (§7.4, advection with streamline-upwind stabilisation) moves to stage 4**, where the
  flow solve gives the velocity; §20.11 moved the flow there. Cost if wrong: one task added to this plan.
- **D. My calls:** ghost weight γ = 0.01 (fact 5); Newton damping (halve the update, down to 1/16, whenever an
  iteration changes T by more than 0.9 of the previous) for the kinks of the heat capacity; Newton tolerance 10⁻⁸;
  the module is `reentry_model/levelset/cutthermal.py` (§7.6 named `levelset/thermal/skfem_cut.py`; stage 1's package
  is flat and uses no scikit-fem assembly); the caller supplies each patch's flux, by projected area for a flow
  (fact 3).

## Global Constraints

- Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; tests run from the repository root with
  `"$PY" -m pytest tests/levelset -q`.
- **Isolation (spec §4.1):** new files only, in `reentry_model/levelset/` and `tests/levelset/`; `cutcell.py` gains
  one function (Task 1); no module of `reentry_model/` outside `levelset/` is edited or imported. Task 5 alone edits
  `README.md`, `docs/model_assumptions.md` and the design.
- The material is duck-typed: any object with `rho` (kg/m³) and vectorised `k(T)`, `cp_eff(T)`, `enthalpy(T)` (J/kg);
  the Scheil material of stage 0 + S will satisfy it.
- Sign convention: level sets negative inside; patch normals point out of the metal; heat flux into the metal is
  positive; SI units.
- Defaults: γ = 0.01, Newton tolerance 10⁻⁸ relative, at most 100 iterations, damping floor 1/16, consistent mass
  (lumped optional), linear solver `amg` (CG with a pyamg smoothed-aggregation preconditioner, tolerance 10⁻¹²) or
  `direct`.
- Work on branch `levelset-two-interface`; commit each task only once Asha has approved execution; stage named paths
  only; never push without asking. Commit messages: a sentence saying what changed, ending
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- The laptop is shared with running flights: tests only on mains power (`pmset -g batt`), the task's test file during
  a task, the whole `tests/levelset` tier once at the end of each task.
- Never edit the `subplan20` worktree, `prototype/work-2026-10-28-skin-spraying*`, `code/reentry_model_output/`, the
  DRAMA databases or the wrapper scripts.

## Review Focus

- **A surface exactly through nodes or along faces** (a grid-aligned plate): zero-volume pieces are not active, no
  NaN weights — pinned by `test_mass_rows_and_weights_on_a_grid_aligned_slab` (Task 2).
- **A sliver of metal (10⁻¹⁰ of a cell):** finite, bounded temperatures there —
  `test_ghost_penalty_keeps_slivers_conditioned_and_leaves_constants_alone` (Task 3).
- **A kink in the heat capacity** (the solidus, Scheil's eutectic step): Newton converges —
  `test_energy_books_close_with_radiation_and_latent_heat` (Task 2) starts just below a linear range's solidus.
- **Metal growing back into cells it left** (a deposit, or the surface moving out): new nodes get finite values within
  their neighbours' range — `test_move_books_the_enthalpy_it_removes_and_fills_new_nodes` (Task 4).
- **A rebuild across the melting range:** the residual is reported, not hidden, and the books close —
  `test_remesh_closes_the_enthalpy_and_reports_the_residual` (Task 4).

---

### Task 1: The metal part of a cut cell as sub-tetrahedra

**Files:**
- Modify: `reentry_model/levelset/cutcell.py` (add `inside_subtets` and `subtet_volumes(S) -> (n,)`; the existing
  `_tvol(a, b, c, d)` keeps its signature)
- Test: `tests/levelset/test_levelset_cutcell.py` (add one test)

**Interfaces:**
- Produces: `inside_subtets(mesh, phi, cells) -> (S (n,4,3) float, owner (n,) int)` — the metal part of each given
  cell as sub-tetrahedra, `owner` indexing into `cells`; one for a lone inside vertex, three for the 2–2 prism (the
  same split as `inside_volumes`) and three for the wedge left by a lone outside vertex. Nodes with φ exactly zero
  count as outside, as in `inside_volumes`.

- [ ] **Step 1: Write the failing test** `test_subtets_reproduce_the_cut_volumes`: on a sphere of radius 10 mm in a
  2 mm box mesh, the cut cells contain all three cases (1, 2 and 3 inside vertices) and the summed `subtet_volumes` per cell equal `inside_volumes` to rtol 1e-12.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_cutcell.py -q` — Expected: FAIL (import error).
- [ ] **Step 3: Implement** from the prototype's `inside_subtets`. Wedge (a, b, c inside; d outside; p_x the cut on
  x–d): (a, b, c, p_a), (b, c, p_a, p_b), (c, p_a, p_b, p_c).
- [ ] **Step 4: Run** the same command — Expected: all pass.
- [ ] **Step 5: Commit** the two files.

### Task 2: The cut-cell solver — geometry, operators and the step

**Files:**
- Create: `reentry_model/levelset/cutthermal.py`
- Test: `tests/levelset/test_levelset_cutthermal.py` (with the test-only materials `ConstMaterial(rho=2800, k=150,
  cp=900, T_ref=293)` and `LinearMeltMaterial(L=390e3, T_s=880, T_l=908, ...)` of the prototype's `simplemat.py`,
  defined in the test file)

**Interfaces:**
- Consumes: `transfer._geom`, `transfer.barycentric`, `cutcell.patches`, `cutcell.inside_subtets`.
- Produces: `SIGMA_SB`; `CutThermal(material, emissivity=0.0, ghost=0.01, lumped=False, newton_tol=1e-8,
  max_iterations=100, linear_solver="amg", cg_tol=1e-12)` with
  - `setup(mesh, phi)`: active cells and nodes; per active cell `Vin` (metal volume), `w` (4 weights, the mean of each
    shape function over the metal part, summing to 1), `Me` (the consistent mass of the metal part, unit
    coefficient), `Ke` (∇N_a·∇N_b·Vin); the patches `P` of the cut cells, their centroid weights `Pw` and nodes;
    the ghost faces (Task 3);
  - `set_temperature(T)` (float or nodal array), attribute `T` (all nodes), `active` (bool per node);
  - `cell_T(T=None)` = Σ w_a T_a per active cell; `energy(T=None)` = Σ ρ·Vin·h(cell_T) (J); `metal_volume()`;
    `patch_T(T=None)`;
  - `operators(T, T_old) -> (K, M)`: k at cell_T; the secant capacity [h(cell_T) − h(cell_T_old)] / (cell_T −
    cell_T_old), `cp_eff` where they differ by under 1e-9 K;
  - `step(dt, q0, q1=None, rad=None, T_amb=0.0) -> float` heat in over the step (J): per patch q = q0 + q1·T_p −
    rad·ε·σ·(T_p⁴ − T_amb⁴), T_p = Σ Pw·T; Newton with Jacobian −A_p·(q1 − 4·rad·εσT_p³)·w w^T; damping as in
    Decision D; raises RuntimeError naming the last change if it does not converge; attribute `iterations`.
- [ ] **Step 1: Write the failing tests** (the prototype's, with the analytic expectations; the grid-aligned slab test was added to the prototype at review and passes there):
  - `test_mass_rows_and_weights` — on the 10 mm sphere: weights sum to 1 (1e-12); the mass rows sum to Vin (1e-12);
    `metal_volume()` equals `inside_volumes(...).sum()` (1e-13) and the sphere's volume within 3 %.
  - `test_mass_rows_and_weights_on_a_grid_aligned_slab` — φ = z on a 1 mm mesh (the surface through a node layer):
    no NaN in `w`, `Vin > 0` for every active cell, `metal_volume()` exact to 1e-12.
  - `test_energy_books_close_with_radiation_and_latent_heat[lumped False/True]` — `LinearMeltMaterial`, start 860 K,
    ε 0.3, six 0.05 s steps of 3 MW/m² on the windward half (q·max(−n_x, 0)) with radiation to 250 K: (E − E0 − Q)/Q
    within 1e-6, and some node above 880 K.
  - `test_constant_flux_slab_converges_at_second_order` — a 4 mm wide, 12 mm deep column, walls off-grid by 0.37 of
    a cell, 2 MW/m² on its top (by projected area) for 0.5 s in 0.01 s steps: the area-weighted top-patch rise differs
    from 2q/k·√(αt/π) by under 1.6 % at 1 mm and 0.5 % at 0.5 mm, the ratio above 3.
  - `test_skin_base_law_is_implicit_and_reaches_its_steady_state` — q0 = 3.2e5·900, q1 = −3.2e5 on every patch,
    five 50 s steps from 300 K: every active node at 900 K within 1e-6.
  - `test_isothermal_radiative_cooling_matches_the_lumped_ode` — k = 1e6, ε 0.8, from 1200 K, ten 2 s steps to 0 K:
    the mean active temperature equals the backward-Euler lumped ODE with the solver's own `metal_volume()` and patch
    area, rel 1e-6.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_cutthermal.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement** from the prototype (setup, operators without the ghost term, step, solvers). Mass
  integrals: the 4-point rule exact for quadratics on each sub-tetrahedron, shape functions by `barycentric` on the
  parent cell. Assembly through a fixed CSR pattern summed by `bincount`, as `thermal/skfem_backend.py`'s `_Pattern`
  (copied, not imported — isolation). The unknowns are the active nodes only; the rest keep their values.
- [ ] **Step 4: Run** the same command — Expected: 7 passed.
- [ ] **Step 5: Commit** the two files.

### Task 3: The ghost penalty

**Files:**
- Modify: `reentry_model/levelset/cutthermal.py`
- Test: `tests/levelset/test_levelset_cutthermal.py`

**Interfaces:**
- Produces: in `setup`, the interior faces between two active cells of which at least one is cut, each with its 8
  node slots, its jump vector [∂N/∂n in cell 1, −∂N/∂n in cell 2] (n the face's unit normal), the scale h·|F| (h the
  larger cell's size) and its two cells; `operators` adds γ·max(k of the two cells)·h·|F|·j j^T to K.
- [ ] **Step 1: Write the failing test:**
  - `test_ghost_penalty_keeps_slivers_conditioned_and_leaves_constants_alone` — sphere of radius 10.001 mm on a 1 mm
    mesh (smallest metal fraction about 6e-10): with γ = 0 the ratio of the smallest to the largest diagonal of
    M/1e-3 + K on the active nodes is below 1e-6; with γ = 0.01 above 3e-3; K times a constant is zero on the active
    nodes (to 1e-12 of max |K|).
- [ ] **Step 2: Run** — Expected: FAIL.
- [ ] **Step 3: Implement** from the prototype's `_ghost_faces` (faces matched by sorting each cell's four sorted
  vertex triples) and the term in `operators`.
- [ ] **Step 4: Run** `tests/levelset/test_levelset_cutthermal.py` — Expected: all pass.
- [ ] **Step 5: Commit.**

### Task 4: The moving surface and mesh rebuilds

**Files:**
- Modify: `reentry_model/levelset/cutthermal.py`
- Test: `tests/levelset/test_levelset_cutthermal.py`

**Interfaces:**
- Consumes: `transfer.project`, `AdaptiveMesh.node_cells`.
- Produces: `CutThermal.move(phi_new) -> float` (enthalpy the move removed, J; negative when metal is added): same
  mesh, nodal T kept, nodes newly active filled ring by ring with the mean of their valued neighbours (donors: the
  previously active nodes), geometry recomputed. `CutThermal.remesh(new_mesh, phi_new) -> float` (the relative
  enthalpy residual before the offset): T extended three rings outside the metal, projected with `transfer.project`,
  the geometry rebuilt on the new mesh, then one uniform offset found by bisection (to 1e-12 K) so the stored
  enthalpy equals the old. The caller supplies φ on the new mesh (by interpolation when the new mesh refines the old,
  fact 9).
- [ ] **Step 1: Write the failing tests:**
  - `test_move_books_the_enthalpy_it_removes_and_fills_new_nodes` — sphere 10 mm, T = 300 + 10⁴·(x + 14 mm):
    shrinking to 9.3 mm returns E_before − E_after (rel 1e-14) and a positive number; growing to 10.4 mm gives every
    newly active node a finite value between 300 K and 300 + 10⁴·28e-3 K.
  - `test_remesh_closes_the_enthalpy_and_reports_the_residual` — `LinearMeltMaterial`, φ = z + 0.3 mm on a 1 mm mesh,
    T = 894 + 2·10⁴·(z + 0.3 mm) + 3·10⁷·(z + 0.3 mm)² (curved and crossing the melting range); rebuilt to 0.25 mm
    within 1.5 mm of the surface: the residual is non-zero and under 1e-2, and the energy after equals the energy
    before, rel 1e-12.
  - `test_receding_surface_through_slivers_gives_no_temperature_jumps` — φ = z − top on a 1 mm mesh of
    6.74 × 6.74 × 11 mm, 1 MW/m² by projected area, sixty 0.02 s steps each followed by `move` 0.05 mm down: the
    temperature at (0.2, 0.1, −3.0) mm has no second difference after the third step above three times the median.
- [ ] **Step 2: Run** — Expected: FAIL.
- [ ] **Step 3: Implement** from the prototype's `_extend`, `move`, `remesh`.
- [ ] **Step 4: Run** the task's file, then the whole `tests/levelset` tier — Expected: all pass (43 stage 1 tests
  plus 12).
- [ ] **Step 5: Commit.**

### Task 5: Cost measurement and documentation

**Files:**
- Create: `analysis/levelset/measure_thermal.py` — builds the 50 mm (or `--diameter`) sphere as fact 10 does, times
  setup and three 0.5 s steps with radiation and 2 MW/m² on the windward half, prints cells, active cells and nodes,
  times, Newton iterations and peak memory, writes a JSON beside it in `reentry_model_output/levelset/`.
- Modify: `README.md` (the level-set section: what stage 2 adds, test count, the cost), `docs/model_assumptions.md`
  §9 (the cut-cell solve: metal-only integration, the cell temperature for coefficients, the secant capacity, the
  ghost penalty and its weight, Newton damping, the uniform offset on rebuilds), the design §20.11 item 2 (facts 6–9
  as stage 2's measured results; Decision A moving the 0.5 s plate judgement to stage 3) and §20.13 (fact 10), and
  §7.4 (Decision C).
- [ ] **Step 1: Write and run** `measure_thermal.py` for 50 mm on mains power. Expected: within 50 % of fact 10.
- [ ] **Step 2: Write the documentation** as listed, every number with what was measured and when.
- [ ] **Step 3: Commit** the named files.

## What stage 2 does not do

No skins (stage 3), so no liquid held on the surface and no plate judgement at 0.5 s steps (Decision A); no flow and
no advection (stage 4, Decision C); no coupling to the flight loop. ψ is built from the computed temperature by stage
1's `molten_rigid_surface`; the plate prototype measures where it sits (fact 6, the non-rigid depth). The fine-size
convergence at 0.7, 0.5 and 0.35 mm (Asha, 2026-11-10) is judged in stage 4 on a pool (§20.11 item 4); fact 6's
0.35 mm plate is a first look.
