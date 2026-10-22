# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level so far: Tasks 1–5, 6a, 6b and 6c
> (2026-10-18/22).** Tasks 7 and 8 are still the scoped outline at the end; each is written out with tested code in
> turn, the same way, and this document regenerated.

**Goal:** The flowing layer moves: every skin cell carries a velocity in its column's tangent plane, coupled through the
depth by the viscosity at each interface, pushed by the gas pressure gradient, the deceleration and the two thickness
pressures, dragged by the gas shear at the top and held at the rigid front; mass, momentum and heat cross between
neighbouring columns layer by layer.

**Architecture:** Everything happens inside the skin body's existing sub-steps; the 3D model, the moving mesh and the
coupled run are unchanged. Task 1 first makes the sideways conduction implicit: one coupled solve per Newton iteration
of each sub-step, the columns' tridiagonal depth solves as its preconditioner. The flow tasks then give the cells a
velocity, couple it through the depth by an implicit tridiagonal viscous solve per column, and carry mass, enthalpy and
momentum between neighbouring columns in an implicit upwind transport of their own, after the momentum step (Task 6a;
ruling, 2026-10-21).

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy (sparse matrices; the conjugate gradients are written out),
scikit-fem for the 3D model, gmsh for its rebuilds.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` — §6.3 (heat: the sideways conduction
made implicit, amended 2026-10-17), §6.4 (flow), §6.5 (sideways exchange and the arrival rule), §11 (build stage 3;
level 1's flow and Rayleigh–Taylor checks; level 2 with flow on), §12 (the risks this stage carries: capillary stiffness
at the grid scale, the flowing-layer coordinate, momentum between tangent planes).

**Builds on:** sub-plan 20 — the columns, the packed cells (with their unused velocity), the surface gradient and
curvature operators, and the base operations. Build stage 2 is complete: levels 1, 2 and 3 pass (fact 149).

**How this plan's code is checked.** Every code and diff block below was written and tested in a planning copy,
`prototype/work-2026-10-18-column-flow-plan/code/`, made from sub-plan 20's tested copy,
`prototype/work-2026-10-12-skin-heat-plan/code/`, as committed at 1d01be5 (identical to the skins build's f11e326,
option A included). This document is generated from it by `prototype/work-2026-10-18-column-flow-plan/make_plan21.py`,
which replays the tasks written so far on a copy of the base — every diff through `patch`, each from the state the
tasks before it left (kept, for files a later task changes again, in `states/task<N>/` as the build repository
committed them) — and checks that the result is the tested copy byte for byte. Never edit this document by hand: change `part_21.md` or the tested copy, and
regenerate.

## Global Constraints

- Work in a fresh work copy, `prototype/work-2026-10-18-column-flow/code/`: a git repository of its own, ignored by the
  outer one, cloned from the skins build (`prototype/work-2026-10-12-skin-heat/code/`, branch `skin-heat` at f11e326).
  Commit there after every task. Never edit `prototype/proto3/`, the package `reentry_model/` at the repo root, the
  specs or the generated plans. Step 3 code stays in the prototype until the skyscrapers architecture is verified
  (Asha, 2026-10-12).
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the work copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Baseline at f11e326: 405 passed, 1 skipped. Record the counts after every task; no failure is acceptable.
- FEniCSx tier: `FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`. Baseline 19 passed.
- `--surface-model elements` (the default) is unchanged by this sub-plan.
- Books: in every skin test and run, the mass is closed to 1e-12 and the energy balance residual is below 1e-8.
- Exit codes: a Newton or coupled solve that does not converge raises `RuntimeError` (exit code 1 at the command line);
  bad settings exit 2.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).
- Quote every path: the repository's path contains a space (`University Documents /MIT`).

## Review Focus

1. **A coupled solve that cannot converge** — it stops the run (exit code 1), never a sub-step taken on an unconverged
   Newton step. Pinned by `test_a_coupled_solve_that_cannot_converge_raises` (Task 1).
2. **Cells that change between sub-steps** (draws, removal, a regrid) — a sub-step's sideways map is rebuilt, never the
   last one's geometry reused. Pinned by
   `test_a_sideways_map_reuses_the_last_geometry_only_while_the_cells_keep_their_thicknesses` (Task 1).
3. **Loosely converged coupled solves** — the inexact Newton steps reach the answer tight ones give. Pinned by
   `test_the_inexact_coupled_solves_reach_the_tight_solve_s_answer` (Task 1).
4. **Neighbouring columns of very different depths** — a melting column beside one at its minimum: the exchange still
   sums to zero and stays symmetric. Pinned by `test_the_sideways_map_sums_to_zero_between_columns_of_any_depth` and
   `test_the_sideways_map_is_symmetric_and_negative_semi_definite` (Task 1).
5. **A sub-step far beyond the old explicit limit** (a finer surface, smaller columns) — bounded, the heat kept. Pinned
   by `test_a_sub_step_far_beyond_the_explicit_limit_stays_bounded_and_keeps_its_heat` (Task 1).
6. **Metal crossing many columns in one sub-step** (Courant numbers near 100 in flight) — positive, never overshooting,
   the books exact. Pinned by `test_courant_numbers_far_above_one_stay_positive_bounded_and_conservative` (Task 6a).
7. **Metal arriving at rigid metal** — it lands on top, never inside, and rigid metal joined to the substrate never
   moves. Pinned by `test_metal_arriving_on_a_rigid_top_lands_on_it_and_the_receiver_s_own_cells_are_untouched`,
   `test_a_floating_crust_carried_across_lands_on_top_and_the_slurry_beneath_it_joins_by_fraction` and
   `test_rigid_metal_joined_to_the_substrate_never_leaves_whatever_its_velocity` (Task 6a).
8. **The transport's upwind order** comes from how scipy numbers strong components, checked every call — were it ever
   wrong, the same answer, more slowly. Pinned by `test_the_upwind_order_changes_the_cost_never_the_answer` (Task 6a).
9. **Surface tension at obtuse corners** — the top surface's own curvature ran backwards there; the thin film's must
   rise with each corner's own height and its product with the faces' exchange have no growing mode. Pinned by
   `test_a_pile_up_lowers_its_own_curvature_at_an_obtuse_corner_of_the_top_but_raises_the_film_s` and
   `test_the_faces_exchange_times_the_film_operator_has_no_growing_mode` (Task 6b).
10. **The predicted heights against the transport's** — the thickness terms are implicit in the heights the thickness
    step predicts, which the transport makes to 0.62 % of the change on the test sphere but only to about 4 % at a
    melting state's rim (6.5 µm of 0.15–0.17 mm, at a Courant number of 1.8; not pinned by a test). Measured in 6c on the
    sphere's draining front: 11–16 % of the largest change at 1 ms, at a few columns of the pile at the shoulder
    (90–91°, 1 mm deep: 3.7·10⁻⁵ m predicted, 1–6·10⁻⁶ m made), and no smaller at shorter sub-steps (0.26–0.60 of the
    change over 50 ms at 2 to 0.25 ms), so it is the columns' linearisation, not the sub-step; at 10 ms, with the
    linearisations at their cap, 13 times the change, and the correction then adds kinetic energy (Task 6c). Where it
    matters, a second linearisation at the transport's own heights would close it.
11. **Rayleigh–Taylor growth under-resolved in time** — at dt·s near or above one, backward Euler damps a growing mode
    spuriously instead of following it. `ThicknessBook.growth` reports it every sub-step; no test can pin what 6c's
    choice of sub-step must avoid (Task 6b). At 6c's 1 ms it is at most 0.0067 (3.4 mm of liquid at 48 m/s²); the
    history's `flow_growth_max` reports each macro step's.
12. **Quiescent layers ringing under surface tension** — two-way flows through the depth join most moving cells into
    one cycle of the transport and its LU fills (4 s a sub-step on the test sphere); a cost, not an error. Measured in
    6c on the 100 mm sphere: 10 315 of 141 000 moving cells in one cycle within 12 sub-steps of 1 ms, 4–5 s a transport,
    and with every column melted flows meeting at the rear stagnation point make cycles of 18 000–57 000 cells, 35–651 s
    a transport. The iteration from the upwind solve measured in 6c (2–9 corrections to 10⁻¹³, 80 times faster than the
    LU with every column melted) is the remedy to build if a flight melts that widely (Task 6c).
13. **A runoff front advances at most one column a flow sub-step** — metal landing on a rigid top becomes a new top
    cell after the transport's solve and flows on only in the next flow sub-step (6a's arrival rule), so a front
    crossing rigid metal moves at most the column spacing per flow sub-step: 2 m/s on 2 mm columns at 6c's 1 ms. 6c's
    choice rests on it (the runoff over the sphere's rear converged at 1 ms, held back by 30–75° at 2–10 ms); a faster
    front later in flight would be held back at 1 ms too, and no history column reports it. Letting the landed metal
    flow on within the sub-step — its new top cell an unknown of the same solve — would lift the limit (Task 6c).
14. **The flow's sub-step is fixed** — ten flow sub-steps a heat sub-step, chosen on one state at 93 s; the Courant
    numbers grow with the layer's speed, its depth and the push late in flight. Watch `flow_courant_max`,
    `flow_growth_max` and `flow_passes_max` (the thickness linearisations at their cap of 8 marked the unstable 10 ms
    runs); an adaptive count from the last sub-step's Courant number is the remedy if they grow (Task 6c).
15. **The flow's cost** — with the flow on, a heat sub-step costs ten flow sub-steps, each about as much as its
    conduction (48 s against 4.5 s on the heated melting state, the transport's large cycles three fifths of it); the
    momentum step's shear-rate fixed point ran to its cap of five passes in every melting state measured, slurry being
    where liquid meets rigid metal (Task 6c's cost).

## Amendment of 2026-10-17 — the sideways conduction becomes implicit (Asha's decision)

> Design of record: spec §6.3 as amended on 2026-10-17; facts 147 and 148. **Write and build this task first, before
> the flow**: it is heat only, so the plates verify it with flow off before anything else changes. **Decided
> (Asha, 2026-10-18): design (a)**, measured against (b) on the 100 mm sphere (fact 150); written out as Task 1 below,
> where the steady two-dimensional test is replaced by the map's symmetry and definiteness and the warm spot's
> convergence, which test the same exchange more directly.

**Why.** Sub-plan 20 conducts heat between neighbouring columns explicitly (`skin.sideways_heat`), sub-cycled inside each
sub-step to stay below the explicit stability limit, and the smallest columns set that limit: on the 100 mm sphere's
2 mm surface the columns' dual areas range from 1.26 to 4.59 mm² (median 3.31), so a cycle is 3.8 ms and every 10 ms
sub-step takes three, where the median column alone would need one. On the level-3 flight of the 100 mm sphere (no
melting; 14.6 h, 72 s a macro step) the sideways conduction took three fifths of a step, against an eighth for the
conduction through the depth that the skins exist to resolve (fact 147). Option A (fact 148) removed the work the
sub-cycles repeated without changing a result; what remains is the sub-cycling itself, which an explicit step cannot do
without. Task 6's transport of mass, enthalpy and momentum, layer by layer on the same depth matching, is already
linearly implicit (spec §6.5); the conduction joins it.

**What.** Backward Euler for the sideways conduction, so a sub-step has no sideways stability limit and no sub-cycles.
The code-level plan chooses, with a measurement of each on the 100 mm sphere's first six steps, between:

- (a) one implicit system per sub-step: the columns' depth tridiagonals coupled by the edges' bin exchanges, solved
  inside the existing Newton iteration by a Krylov method preconditioned by the depth solve (Thomas on every column).
  The sub-step stays one backward-Euler step.
- (b) a split sub-step: the depth solve, then an implicit sideways solve on the depth grid's bins (a graph Laplacian over
  the columns in each bin, with the exchange's k w l weights), the bins' heat given back to the cells by the existing
  overlaps. Simpler, but it adds a splitting error of first order in the sub-step.

Either keeps the exchange antisymmetric edge by edge, so the books close to round-off, and exact for a field linear
along the surface; either then carries Task 6's implicit transport.

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

**Shared later from this sub-plan:** the slurry viscosity law it reads; the column flow is skyscrapers' own.

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

- The viscosity in the cells comes from **`reentry_model/shared/rheology.py`** (S1); the column flow is skyscrapers' own.

---

### Task 1: Implicit sideways conduction (design (a))

**Files:**
- Modify: `reentry_model/columns.py` (`Columns.incidence`), `reentry_model/skin.py` (the coupled sub-step and the
  `Sideways` map, which replaces the explicit exchange), `reentry_model/skinbody.py` (the hook that makes the map)
- Test: `tests/test_reentry_model_skin.py` (the sideways tests, Task 5 of sub-plan 20, rewritten)

**Interfaces:**
- Consumes: sub-plan 20's `skin._sideways_geometry`, `skin._interp_weights`, `skin._interp_apply`, `skin.depth_grid`,
  `skin.thomas`, `Columns.edges` and `Columns.conduct_weight`.
- Produces: `skin.Sideways(sk, columns, previous=None)` with `.heat(x) -> (n, K)` [W per cell], linear in the cell
  temperatures x, and `.bins(x) -> (n, bins)`; `skin._conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, sideways=None)`;
  `skin.advance(..., side=...)` with `side(sk, previous) -> Sideways`; `Columns.incidence` (scipy CSR, columns by
  edges); `skin.KRYLOV_TOL`, `skin.KRYLOV_MAX` and the diagnostic counter `skin.KRYLOV`.
- Retires: `skin.sideways_heat`, `skin._sideways_once`, `skin.SIDEWAYS_CFL`, and sub-plan 20's pin of the explicit heat
  to its build (`test_the_sub_cycles_share_their_geometry_and_give_the_heat_as_built_to_the_last_bit`).

The sideways exchange keeps sub-plan 20's construction: each column's cells averaged onto the common depth grid, every
surface edge exchanging k w l (T_j − T_i) in every bin both columns reach, each column's bin heat given back to its
cells by their overlap with the bin. Written as a map on the cells, S = Pᵀ L P (P the averaging onto the bins, L the
edges' conductances, Pᵀ the return to the cells), it is symmetric and negative semi-definite, so it can join the
implicit sub-step: every Newton step solves J dT = −R with J the columns' tridiagonal depth Jacobians less S, by
conjugate gradients preconditioned by the columns' own tridiagonal solves (`_coupled_solve`). Three choices keep it
cheap, each measured on the 100 mm sphere's first six steps (fact 150):

- **Inexact Newton steps.** Each coupled solve stops at `KRYLOV_TOL` = 10⁻³ of its right-hand side and the Newton
  iteration corrects the rest; started from the columns' own solve, it takes about one conjugate-gradient iteration
  (1.13 measured). Solving each step to 10⁻⁶ instead moves no column of the history by more than 6 parts in a billion
  and costs 61 % more.
- **The map is linear**, so the line search's first try takes the sideways heat at T + dT from the solve's own
  products instead of evaluating it again; only a column whose step is cut back costs an evaluation.
- **The geometry is reused** from the last sub-step while the cells keep their thicknesses — in a body that does not
  melt, every sub-step; k is taken afresh at each sub-step's start.

The sensitivity the trial pass tracks for the base law (dT/dT_b) stays through the depth only: the law each column gives
the 3D model is its own, and the neighbours' share is the interface mismatch's to book. Against (b) — the same
exchange split from the depth solve and solved on the bins alone — (a) has no splitting error and costs more: 260 s for
the six steps against 200 s for (b) and 254 s for sub-plan 20's explicit exchange after option A, its processor time
(224 s) a little below the explicit exchange's (235 s).

- [ ] **Step 1: Replace the sideways tests** in `tests/test_reentry_model_skin.py` (everything after `grid_skins`, the
  helper the tests share):

{{DIFF_TEST_SKIN_T1}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q -k "sideways or
  coupled or warm_spot or explicit_limit"`. Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'Sideways'`.

- [ ] **Step 3: The edges' incidence matrix** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T1}}

- [ ] **Step 4: The coupled sub-step and the sideways map** in `reentry_model/skin.py`:

{{DIFF_SKIN_T1}}

- [ ] **Step 5: The skin body's hook** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T1}}

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` and
  `tests/test_reentry_model_columns.py` 31 passed (the warm spot's grid error 1.29 K; its time errors at 40, 20 and
  10 ms 2.79, 1.32 and 0.57 K against a 2.5 ms run, as first order predicts; the checkerboard smoothed to 622–629 K in
  one 0.2 s sub-step with its heat kept to 3.5·10⁻¹³). Then the unit tier: 409 passed, 1 skipped; and the FEniCSx tier:
  19 passed.

- [ ] **Step 7: Level 2 with the implicit exchange.** Run the four plates, flow off, side by side:
  `"$PY" analysis/skin_plate.py --case cold|hot --dt 0.05|0.5 --out <json>`. Each must pass as it did with the
  explicit exchange (spec §11.1: recession within 5 % of the references at 1.5 and 2.0 s, and at 1.0 s for the cold
  plate; for the hot plate the flowing depth within 10 % and the base within 1 K of T_rigid; books to 10⁻¹² and 10⁻⁸).

  Measured (2026-10-18, the four run side by side, each once): all four pass, and the implicit exchange moves no
  judged number by more than 0.2 points from sub-plan 20's explicit exchange (B′ with the look-ahead, facts 145-146).
  The table gives the centre's recession against the one-dimensional reference at 1.0, 1.5 and 2.0 s (negative is
  late), with the explicit exchange's in brackets:

| Plate, macro step | 1.0 s | 1.5 s | 2.0 s | Other criteria |
|---|---|---|---|---|
| Cold, 0.05 s | +3.7 % [+3.6] | +0.2 % [0.0] | −0.9 % [−0.9] | fine box +1.7, −0.5, −0.6 % [+1.6, −0.6, −0.6] |
| Cold, 0.5 s | +3.2 % [+3.2] | +0.5 % [+0.4] | −0.9 % [−0.9] | fine box +1.2, −0.2, −0.7 % [+1.2, −0.2, −0.6] |
| Hot, 0.05 s | −3.1 % (not judged) [−3.1] | −1.5 % [−1.5] | −0.8 % [−0.8] | flowing depth within 0.9 %, base at most 0.08 K above T_rigid |
| Hot, 0.5 s | −9.0 % (not judged) [−9.0] | −3.0 % [−2.8] | −1.1 % [−1.1] | flowing depth within 1.1 %, base at most 0.09 K above T_rigid |

  Mass closed to 4·10⁻¹⁶ and energy to 2.2·10⁻¹⁰ (the cold plate at 0.05 s) or better; no base faults; 84, 76, 12 and
  12 interface repeats. The plates' heated region tapers along a half cosine from 32 to 45 mm, so the implicit exchange
  carries heat at its rim, but the judged centre is one-dimensional: these runs show the change leaves level 2 where
  it was, not that the exchange is more accurate. Their wall times (598–1 546 s) were taken four at a time and are not
  comparable with sub-plan 20's, timed alone.

- [ ] **Step 8: Commit.**

```bash
git add reentry_model/columns.py reentry_model/skin.py reentry_model/skinbody.py tests/test_reentry_model_skin.py
git commit -m "Skins: the sideways conduction is implicit, solved with each sub-step's depth conduction by conjugate gradients preconditioned by the columns' own solves (sub-plan 21 Task 1, design (a))"
```

---

### Task 2: The cells' velocity and the surface operators

**Files:**
- Modify: `reentry_model/skin.py` (`Skins.U`, `Skins.momentum_out`, `Skins.column_momentum`, `Skins.momentum`; every
  cell operation carries the velocity), `reentry_model/columns.py` (`triangle_geometry`, `Columns.gradient`,
  `Columns.curvature`)
- Test: `tests/test_reentry_model_skin.py`, `tests/test_reentry_model_columns.py`

**Interfaces:**
- Consumes: Task 1's `skin.py` and `columns.py`.
- Produces: `Skins.U` (n, K, 3) [m/s]; `Skins.momentum_out` (n, 3) [kg m/s], what the last `take_top` or
  `take_bottom` took away; `Skins.column_momentum() -> (n, 3)`, `Skins.momentum() -> (3,)`;
  `Skins.add_bottom(dm, h, s=1.0, u=None)` (u the arriving metal's velocity, at rest by default);
  `Skins.state()` now six arrays, the velocity last; `columns.triangle_geometry(points, faces) -> (area, unit_normal,
  cot)` with cot (nt, 3) the cotangent at each corner; `Columns.gradient(f, points=None) -> (n, 3)`;
  `Columns.curvature(points=None) -> (n,)` [1/m].

What sub-plan 20 left to this stage (its ruling 7): the velocity the flow will give each cell, and the two surface
operators the push reads. (Its ruling 9's interface rule was settled after that build: damped repetition, Asha,
2026-10-12.) Nothing here moves a cell yet — every velocity stays zero until Task 4 — so this task changes no result: the
100 mm sphere's first six steps give Task 1's history in every column (processor time 221 s against 224 s).

**The velocity** is a vector per cell, stored in three components so that Task 6 can project it between neighbouring
columns' tangent planes. It is part of the cells' state, saved and restored with them, and every cell operation carries
it so that each column's momentum is exact: a split's halves keep their cell's velocity, a merged cell takes the
mass-weighted mean, a partly taken cell keeps its own, and metal drawn from the 3D model arrives at rest (it is rigid)
unless a velocity is given. What leaves at the top or the bottom takes its momentum with it, recorded in
`momentum_out` for the momentum books of Task 8.

**The gradient** of a per-corner field averages the gradients of its linear interpolation on the corner's triangles, by
area, and projects the result onto the corner's tangent plane: exact for a field linear on a flat surface, within
0.09 % of a linear field's tangential part on average on the 4 mm test sphere (0.8 % at its one irregular corner) and
0.02 % on the 2 mm sphere. **The curvature** is the mean curvature, the sum of the principal curvatures, from the
cotangent Laplacian of the position. Sub-plan 20 measured it over the barycentric dual area, right on average but
scattered 7.6 % from corner to corner on the test sphere, which the push's −σ ∇_s κ would turn into flow; over the
mixed Voronoi area (Meyer et al. 2003) the scatter falls to 0.45 % on the 2 mm sphere, and what is left is at the
corners touching an obtuse triangle, where that area is only approximate (33 of the test sphere's 2 254 corners, off by
up to 176 %; 48 of the 2 mm sphere's 9 417, up to 23 %). Every other corner is within 1.1·10⁻⁴ of 2/R, so a corner
touching an obtuse triangle takes the mean of its neighbours that touch none, or of those filled in the pass before:
every corner of both spheres is then within 8.4·10⁻⁵ (test sphere) and 1.3·10⁻⁵ (2 mm) of 2/R (measured 2026-10-18,
`prototype/work-2026-10-18-column-flow/measure/curvature.py`).

- [ ] **Step 1: Append the velocity tests** to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T2}}

- [ ] **Step 2: Append the surface operators' tests** to `tests/test_reentry_model_columns.py`:

{{DIFF_TEST_COLUMNS_T2}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_columns.py -q`. Expected: `AttributeError: 'Skins' object has no attribute 'U'` and
  `AttributeError: module 'reentry_model.columns' has no attribute 'triangle_geometry'`.

- [ ] **Step 4: The velocity in the cells** in `reentry_model/skin.py`:

{{DIFF_SKIN_T2}}

- [ ] **Step 5: The surface operators** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T2}}

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 29 passed and
  `tests/test_reentry_model_columns.py` 9 passed; the unit tier 416 passed, 1 skipped. The first six steps of the
  100 mm level-3 flight (`--heating sesam --surface-model skin --t-max 3`) give Task 1's history in every column.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/skin.py reentry_model/columns.py tests/test_reentry_model_skin.py tests/test_reentry_model_columns.py
git commit -m "Skins: every cell carries a velocity with its momentum exact through split, merge, take and add, and the corner graph gains the surface gradient and the mean curvature over mixed Voronoi areas (sub-plan 21 Task 2)"
```

---

### Task 3: Viscosity per interface

**Files:**
- Create: `reentry_model/rheology.py`, `tests/test_reentry_model_rheology.py`
- Modify: `reentry_model/skin.py` (`cell_viscosity`, `interface_drag`), `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `skin.py`; `material.Material.melts`, `.liquid.mu` and `.liquid_fraction`; `Skins.props`.
- Produces: `rheology.li_viscosity(f_s, gamma)` and `rheology.viscosity(T, f_s, gamma, mu_liquid, f_s_bridge)` [Pa s],
  with `rheology.GAMMA_MAX` (367 1/s), `BRIDGE_LOW` and `BRIDGE_HIGH` (906.4 and 908 K) and `RIGID_SOLID_FRACTION`
  (0.5); `skin.cell_viscosity(sk, gamma) -> (n, K)` [Pa s], infinite in rigid and missing cells;
  `skin.interface_drag(sk, gamma) -> (n, K-1)` [Pa s/m], gamma the interfaces' shear rates [1/s].

The viscosity law is the Step 4 plan's (§4), as spec §2 records it: liquid aluminium's 1.3 mPa·s at and above the
liquidus; a log-linear bridge over 906.4–908 K from Li et al.'s value where Scheil's curve stands at 10 % solid (9.75 %
measured); Li Yageng et al. (2014), η = [0.871 − 0.00849 γ̇^0.74924] exp(3.7311 f_s) Pa·s, from there down to half
solid, the shear rate held at no more than 367 1/s (its bracket would reach zero at 483 1/s); and infinite at and above
half solid, where a cell is rigid and anchored. At rest that is 0.87 Pa·s at f_s = 0 and 5.6 Pa·s at f_s = 0.5, 670 to
4 330 times the liquid. It lives in `rheology.py`, which imports neither architecture's modules, so that the later move
into `reentry_model/shared/` stays mechanical (the 2026-10-12 amendment); the Step 4 plan's Chen branch above half
solid is not needed where rigid cells are anchored, and is left to that plan.

**Ruling (2026-10-19):** the drag between neighbouring cells puts their two half-cells in series,
D = 1 / (dz_i / 2η_i + dz_j / 2η_j) per unit area, each cell's viscosity taken at the interface's shear rate — the
spec's "viscosity at their interface" (§6.4) in the form the stress's continuity across the interface requires, the way
`_conduct` combines conductivities. Across a jump from liquid to slurry the liquid's half dominates and the profile
Task 4 builds on it has the exact kink; a rigid cell resists without limit and adds nothing to the sum, so a flowing
cell drags on the rigid cell beneath over its own half-thickness alone — no slip at the rigid front. The fixed-point
iteration on the shear rate (at most five passes) is Task 4's: the rate comes from the velocities it solves for.
(Corrected in Task 4, 2026-10-20: each cell's viscosity is taken at its own shear rate, not the interface's, which
would give a slurry half-cell under liquid the liquid's rate.)

- [ ] **Step 1: Write the rheology tests**, `tests/test_reentry_model_rheology.py`:

{{TEST_RHEOLOGY_T3}}

- [ ] **Step 2: Append the drag tests** to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T3}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_rheology.py
  tests/test_reentry_model_skin.py -q`. Expected: `ModuleNotFoundError: No module named 'reentry_model.rheology'`.

- [ ] **Step 4: Create `reentry_model/rheology.py`:**

{{RHEOLOGY_T3}}

- [ ] **Step 5: The cells' viscosity and the drag between them** in `reentry_model/skin.py`:

{{DIFF_SKIN_T3}}

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_rheology.py` 5 passed and
  `tests/test_reentry_model_skin.py` 33 passed; the unit tier 425 passed, 1 skipped. No result changes: nothing reads
  the viscosity until Task 4.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/rheology.py tests/test_reentry_model_rheology.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the viscosity of every cell from Li et al.'s slurry law with its bridge to the liquid, and the drag between neighbouring cells with their half-cells in series (sub-plan 21 Task 3)"
```

---

### Task 4: Momentum through the depth

**Files:**
- Modify: `reentry_model/skin.py` (`series_drag`, and `interface_drag` at each cell's own shear rate; `SHEAR_PASSES`,
  `SHEAR_TOL`, `CRUST_VISCOSITY`, `MomentumBook`, `anchored_cells`, `momentum_step`), `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `Skins.U` and `Skins.column_momentum`; Task 3's `cell_viscosity`; `thomas`.
- Produces: `skin.series_drag(sk, eta) -> (n, K-1)` [Pa s/m] from the cells' viscosities; `skin.interface_drag(sk,
  gamma)` now with gamma each cell's own shear rate; `skin.momentum_step(sk, dt, push, shear, A_top) -> MomentumBook`
  with push (n, 3) [N/m3] the same on every layer of a column, shear (n, 3) [Pa] on the top cell and A_top (n,) [m2];
  `MomentumBook(push, shear, wall, passes)`, the impulses [N s] (n, 3); `skin.anchored_cells(sk, flowing) -> (n, K)`
  bool, the rigid cells joined to the substrate; `skin.SHEAR_PASSES` (5), `skin.SHEAR_TOL` (10⁻³) and
  `skin.CRUST_VISCOSITY` (10⁹ Pa·s).

Every column's momentum through its depth in one implicit sub-step (spec §6.4), each velocity component apart: the
cells' mass over the sub-step, the drag between neighbours from Task 3, the push on every flowing cell's volume and the
gas shear on the top cell's top area — a tridiagonal system per column and component, solved by `thomas` as the heat
is. The rigid metal joined to the substrate — each column's run of rigid cells up from its bottom cell, which rests on
the rigid 3D model — is held at rest (its momentum, if it had any, is booked to the wall); a moving cell drags on it,
and a flowing bottom cell on the rigid 3D metal beneath the base, with no slip over its own half-thickness. **A crust
floats** (Asha, 2026-10-20, option B, revising the spec's assumption that a crust is anchored): rigid cells with flowing
metal between them and the substrate move as one piece, given a viscosity of 10⁹ Pa·s, dragged by the gas shear on
their top, pushed like the rest and dragged by the slurry beneath. The columns carry no sideways stress in a solid, so
a crust over many columns moves as column-sized plates — as a warm crust, barely solid near half liquid, would — and a
crust held at its rim is not held; anchoring a crust only where it touches rigid metal at the surface (option C) is
kept as a refinement should Task 8 find crusts over slurry common.
The shear rates come from a fixed point started at the sub-step's velocities: each pass solves the system with the
cells' viscosities, then takes every flowing cell's shear rate as the mean over its two faces of the stress there over
its viscosity, and stops after five passes or when no viscosity moves by more than 10⁻³. **Each cell at its own rate**
(found while writing this task): Task 3 took both half-cells at an interface at that interface's rate — the velocity
difference over the distance between the centres — which at a jump from liquid to slurry mixes the liquid's
770 1/s with the slurry's 0.3 1/s and gives the slurry's half-cell five times too low a viscosity; at its own rate the
two-layer profile is exact. The impulses (push, shear, the wall's drag and held momentum) add up to the momentum
change to round-off, for Task 8's momentum books. Nothing calls `momentum_step` yet: moving the metal is Task 6's, so
no result changes.

Measured against exact solutions (2026-10-19/20): shear alone, Couette's line to 1.1·10⁻¹⁰ in one pass; the push
alone, the half-parabola to 2.8·10⁻⁴ of its top speed at 30 cells, 6.9·10⁻⁵ at 60 and 1.7·10⁻⁵ at 120 (second
order); 1 Pa on liquid over slurry at 900 K (37 % solid), the slurry at 3.46 Pa·s shearing at 0.29 1/s, both
profiles to 1.1·10⁻⁵ in two passes; Stokes' first problem (1 Pa suddenly on 0.9 mm of liquid, 18.5 ms, a diffusion
length of 100 µm), the top cell to 1.3·10⁻³, 6.6·10⁻⁴ and 3.3·10⁻⁴ at 100, 200 and 400 steps (first order, backward
Euler), the momentum the shear's impulse exactly; a crust of three rigid cells on 15 of slurry at 900 K under 1 Pa of
gas shear rides the top of the slurry's Couette line as one piece, to 10⁻⁴.

The layer's response develops over time, not at once (Asha's question, 2026-10-20): the push acts on every layer from
the moment it is applied — the pressure crosses a millimetre of metal in about 10⁻⁷ s and the deceleration acts on all
of it — but the velocity spreads from the anchored front, and the gas shear's from the top, over the viscous time
h²/ν: 0.018 s for 0.1 mm of liquid, 0.46 s for 0.5 mm, 1.85 s for 1 mm, and 0.2 ms for 0.5 mm of slurry at 3.5 Pa·s,
whose high viscosity carries momentum fastest. The momentum sub-step keeps the cells' inertia, so it follows that
development. Started from rest and compared over 1 s with 0.1 ms sub-steps on the same cells, the 10 ms sub-steps were
off by 15 % of the speed at 1 s in the first sub-step for 0.1 mm of liquid under 10 kPa/m, which settles within about
two sub-steps, and by 0.04 % after 0.1 s; by at most 1.0 % for 0.5 mm of liquid under the same push and 1.9 % under
10 Pa of gas shear; and 1.1 % in the first sub-step for 0.5 mm of slurry at 900 K. In flight the forcing changes over
seconds, so only the thinnest liquid layers' first sub-step after a sudden onset lags.

- [ ] **Step 1: Append the momentum tests** to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T4}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `AttributeError: module 'reentry_model.skin' has no attribute 'momentum_step'`.

- [ ] **Step 3: The momentum sub-step, and the drag at each cell's own shear rate,** in `reentry_model/skin.py`:

{{DIFF_SKIN_T4}}

- [ ] **Step 4: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 40 passed (with Task 3's
  drag tests, which hold at one rate for all cells); the unit tier 432 passed, 1 skipped.

- [ ] **Step 5: Commit.**

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: every column's momentum through its depth in one implicit sub-step, each cell's viscosity at its own shear rate by fixed point, the rigid metal joined to the substrate held with no slip and a crust floating as one piece, the books exact (sub-plan 21 Task 4)"
```

---

### Task 5: The push

**Files:**
- Modify: `reentry_model/surface_flow.py` (`SurfaceFlowResult.dp_ds`), `reentry_model/skin.py` (`along_meridian`,
  `push`)
- Test: `tests/test_reentry_model_surface_flow.py`, `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `Columns.gradient`, `Columns.curvature` and `Columns.top_points`; the surface flow's
  `SurfaceFlowResult`.
- Produces: `SurfaceFlowResult.dp_ds` (n,) [Pa/m], the wall pressure's gradient along the meridian on the windward
  face, zero leeward; `skin.along_meridian(normals, v_hat, magnitude) -> (n, 3)`, vectors in each corner's tangent plane
  away from the stagnation point (for the pressure term here and the gas shear, `SurfaceFlowResult.tau`, in Task 6);
  `skin.push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma) -> (G (n, 3) [N/m3], terms)`, the terms a dict of
  `pressure`, `deceleration`, `thickness` and `capillary`.

The push on every layer of a column (spec §6.4), G = −∇_s p_w + ρ_l g_t + ρ_l (g·n) ∇_s H − σ ∇_s κ, each term in the
column's own tangent plane, where Task 4 solves for the velocities. g is the acceleration the metal feels in the body's
frame — minus the body's own; for a body slowed by drag, the deceleration along the flight direction — so it pulls the
metal toward the nose along the surface and out of the front face: the spec's −ρ_l a_t and ρ_l a_n (it writes the
first with the body's acceleration and the second with the deceleration's pull), and the element model's
−ρ_l W sin θ on the windward face. Unlike the element model it acts on the leeward face too, where it pulls the metal
toward the equator. H is each column's thickness and κ the mean curvature of the top surface the thicknesses make.
(Task 6b replaces this κ by its thin-film form, on the top surface's cotangents, and applies surface tension only where a
column's top flows: Task 2's obtuse-corner rule made this term run backwards at those corners, Asha, 2026-10-21.)

**The pressure term takes the surface flow's gradient, not the corners'.** Differentiating each corner's own p_w would
carry the scatter of the corners' normals — 0.06° off radial on average on the 2 mm sphere, 0.87° at worst — amplified
some 25 times, the radius over the spacing: a Newtonian pressure's push came out 1.1 % off on average and 13 % at worst
that way (1.5 % and 10 % on the 4 mm test sphere), against 0.1 % and 1.2 % with each corner's exact angle; normals
taken from the curvature vector, or averaged over the neighbours, did no better (measured 2026-10-20,
`prototype/work-2026-10-18-column-flow/measure/`). The surface flow already differentiates p_w(θ) on its own fine bins
along the meridian for the element model's G; it now also returns that gradient, `dp_ds`, on its own, with every
field it returned before unchanged to the last bit at 20, 80 and 200 s of the 100 mm flight, and the push takes
−dp_ds along the meridian away from the stagnation point.

**Sizes at peak heating** on the 100 mm sphere (estimated, not run: 93 s, 60 km, a dynamic pressure of 6.1 kPa from the
level-3 physics flight, a stagnation pressure of 1.84 times it and a drag coefficient of 0.92): the pressure term peaks
at p_stag / R, about 2.2·10⁵ Pa/m at 45° from the nose; the deceleration, 30 m/s², gives at most ρ_l g = 7.2·10⁴ Pa/m
at 90°, pulling the other way. Their ratio, 1.84 × 4ρ_s / (3 ρ_l C_D) ≈ 3.1, is the same for a solid sphere of any size
at any point of its flight, since both scale with the dynamic pressure. The thickness and capillary terms matter at
short wavelengths: a 0.1 mm bump 4 mm long gives about 3·10⁵ Pa/m of capillary push, which is why Task 6 takes both
implicitly; the sphere's own mesh, with the curvature exact to 10⁻⁴, gives at most 11 Pa/m of spurious capillary push
(0.3 on average).

The two thickness terms are evaluated at the thicknesses given; their implicit treatment, linearised in the column
heights, belongs to Task 6's solve, and so do the Rayleigh–Taylor growth and the surface-tension cutoff's decay, which
need the thickness to evolve (Asha, 2026-10-20). Nothing calls `push` yet: the skin body evaluates the surface flow on
its corner surface and drives `momentum_step` with the push and the gas shear when the flow joins its sub-steps
(Task 6), so no result changes.

- [ ] **Step 1: The surface flow's test** — append to `tests/test_reentry_model_surface_flow.py`:

{{DIFF_TEST_SURFACE_FLOW_T5}}

- [ ] **Step 2: The push's tests** — append to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T5}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_surface_flow.py
  tests/test_reentry_model_skin.py -q`. Expected: `AttributeError: 'SurfaceFlowResult' object has no attribute 'dp_ds'`
  and `AttributeError: module 'reentry_model.skin' has no attribute 'along_meridian'`.

- [ ] **Step 4: The pressure gradient on its own** in `reentry_model/surface_flow.py`:

{{DIFF_SURFACE_FLOW_T5}}

- [ ] **Step 5: The push** in `reentry_model/skin.py`:

{{DIFF_SKIN_T5}}

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_surface_flow.py` 8 passed and
  `tests/test_reentry_model_skin.py` 46 passed; the unit tier 439 passed, 1 skipped.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/surface_flow.py tests/test_reentry_model_surface_flow.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the push along the surface on every layer -- the surface flow's pressure gradient along the meridian, the deceleration, the thickness and the capillary terms -- in each column's tangent plane (sub-plan 21 Task 5)"
```

---

### Task 6a: Sideways transport at given velocities

**Files:**
- Modify: `reentry_model/columns.py` (`Columns.face_vector`), `reentry_model/skin.py` (`SLIVER`, `TransportBook`,
  `tangent_frames`, `_fractions`, `_pieces`, `_upwind_order`, `_factor`, `transport`)
- Test: `tests/test_reentry_model_columns.py`, `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `Skins.U`, Task 3's `cell_viscosity`, Task 4's `anchored_cells`; `Columns.edges`, `.normal`,
  `.area`; `Skins.dz`, `Skins._install` and `Skins.regrid`.
- Produces: `Columns.face_vector` (E, 3) [m], each dual face's normal times its length, from the edge's first corner
  to its second; `skin.transport(sk, columns, dt, regrid=True) -> TransportBook`, which moves the cells at their
  velocities `sk.U` and, with `regrid`, splits and merges them; `TransportBook(mass, enthalpy, momentum, wall, courant)`,
  per column what it gained from its neighbours less what it gave them — mass (n,) [kg], enthalpy (n,) [J], momentum
  (n, 3) [kg m/s] as it left the donors — the wall's impulse (n, 3) [N s], the normal parts removed, and the largest
  cell's outflow rate times dt; `skin.tangent_frames(normals) -> (e1, e2)`, each (n, 3); `skin.SLIVER` (10⁻¹²); and,
  private, `_fractions(depth) -> (n, K+1)`, `_pieces(sig, tau) -> (row, cell of sig, cell of tau, length)`,
  `_upwind_order(src, dst, n) -> (order, upwind)` and `_factor(A, upwind)`.

The metal moves sideways at the velocities the cells already have (spec §6.5); in the tests they are prescribed, and
the skin body drives them with Task 4's `momentum_step` in 6c. The transport is its own step in each sub-step, after
the momentum step (ruling, 2026-10-21); the conduction stays in Task 1's solve. Across each face, cell k of the donor
column a sends c m_k′ toward its neighbour b, with c = (u_k·f)⁺ / A_a [1/s] — the velocity's component across the face
times the face's length, per unit of the column's area — and m_k′ the cell's mass at the sub-step's end; the metal
carries the cell's specific enthalpy and its velocity. Every cell then satisfies m_k′ (1 + dt Σc) = m_k + dt (what
arrives): (I + dt C) m′ = m, the film's runoff scheme (`film.Runoff`) on the cells. The matrix is an M-matrix whose
columns sum to one, so the masses stay positive, the books close to round-off and nothing limits dt; the enthalpy and
the trial pass's sensitivity S solve the same system, and so does the momentum, with every arrival projected onto its
receiver's plane. Only the cells that send or receive are solved for: with every velocity zero, or only rigid metal
joined to the substrate given one, nothing changes at all.

**Where the metal goes** (spec §6.5, §14 assumptions 3 and 4). A flowing donor cell's span of fractions of the donor's
flowing depth joins the receiver's cells between the same fractions of the receiver's flowing layer — the run of
flowing cells from its top — each in proportion to the overlap. Where the receiver's top is rigid (anchored metal, or
a crust over flowing metal) everything arriving lands on top of its top cell, as a new cell; so does a floating
crust's metal, wherever it goes. Rigid metal joined to the substrate neither leaves nor receives, whatever velocity it
holds. After the exchange the cells are split and merged as everywhere else (spec §6.2): a deposit thinner than half
its target merges into the cell beneath it, as any thin cell does.

**The faces.** A velocity crosses each dual face at the face's vector, its two segments' normals times their lengths
(`Columns.face_vector`), not the edge's direction times the face's length: a column's faces close exactly on a plane,
so a uniform flow carries as much into every column as out of it, where the edge's direction gathered or spread a
column's metal at up to 539 1/s under 1 m/s on the columns' tests' jittered plane (2.9·10⁻¹³ with the faces' vectors).
On the 4 mm test sphere the faces close in the tangent plane to 3.7 % of a face length at worst and 0.03 % on average;
the rest is along the normal, the curvature, which a tangential velocity does not see.

**The momentum is projected at every arrival.** Each cell's momentum is held in two components in its own plane
(`tangent_frames`), and every arrival enters its receiver through the 2×2 projection between the two planes, so a
second system of twice the size is solved with the same structure; any normal part a cell held is booked to the wall
too. Projecting once at the sub-step's end instead — the Cartesian components carried by the scalar system, each cell
then projected — costs nothing extra but carries metal that crosses several columns within a sub-step straight on, off
the curved surface: on the 100 mm sphere a uniform meridian flow kept 99.7 % of its kinetic energy in one sub-step at a
Courant number of 1 against 99.8 % projected at arrival, 90 % against 96 % at 9.6, and 24 % against 49 % at 96, where
the two differed by up to 1.8 times the largest cell momentum. Courant numbers near 100 are what the flight will see:
0.5 mm of liquid under the 2.2·10⁵ Pa/m of Task 5's pressure push moves at about 20 m/s at its top (G h²/2μ), 100 columns
of 2 mm in a 10 ms sub-step.

**The solve, measured** (2026-10-21, `prototype/work-2026-10-18-column-flow/measure/task6a/`, on the 100 mm sphere's
skins, 9 417 columns of 10 µm cells, flowing at 0.1 m/s along the meridian, every cell liquid and moving: 584 000
unknowns, 2.3 million arrivals). A sparse direct solve in the cells' own order took 19.6 s to factor; in SuperLU's
fill-reducing order (COLAMD) 0.66 s; minimum-degree 54 s. Taken in upwind order, every donor before its receivers, the
system is triangular but for its cycles where flows meet (here at most 7 cells, near the rear stagnation point), and
SuperLU with no pivoting factors it in 0.06 s with no fill at all — and without pivoting an M-matrix's elimination keeps
every sign, so the masses come out non-negative to the last bit. The order is scipy's strong components, 0.06 s, which
its algorithm (Pearce's) numbers in reverse topological order; that is how scipy numbers them, not what it promises,
so `_upwind_order` checks it every call and `_factor` falls back to COLAMD with pivoting, as exact and slower, should
it fail. Not built (estimated): an upwind sweep written in numpy would need the levels of that order, a hundred or more
relaxations over the 2.3 million arrivals; Gauss–Seidel or a Krylov iteration gains about ν/(1 + ν) a sweep, so far
above a Courant number of one it needs many, and its books close only to its tolerance. The pieces where a donor's and
a receiver's partitions cut each other come from a rank merge (`_pieces`, 0.13 s against 0.22 s for a row-wise sort).

Whole transport sub-steps (without the regrid, which the sub-step does anyway): 0.95 s with all 584 000 cells moving,
of which the momentum's system 0.30 s and the pieces 0.13 s; 0.27 s in a melting state — the front hemisphere's top
0.3 mm liquid over rigid metal, 124 000 of 553 000 cells moving along a half-parabola in depth — the same at Courant
numbers of 1, 9.6 and 96. On the same machine a coupled conduction sub-step on the sphere's skins at 700 K takes
0.6–0.8 s and a momentum sub-step on the melting state 0.54 s; the regrid after the transport 0.09–0.25 s, more at high
Courant numbers, where more cells change size.

**Rulings** (2026-10-21, within spec §6.5 and the rulings of 2026-10-20/21):

- *The face's vector, not the edge's direction times the face's length* — the face is two segments in two triangles,
  and only their vectors close around a column — if wrong, nothing: on symmetric triangles the two are the same.
- *The donor's flowing depth is all its flowing cells, top down; a moving rigid cell is no part of it and lands on top*
  — the spec's "flowing layer"; where flowing metal lies under a floating rigid layer (melt beneath a refrozen deposit),
  it is matched by its share of the flowing metal and the rigid layer lands on the receiver's top — if wrong, such a
  layer arrives above metal that was above it; rare, and transient.
- *A donor cell's pieces are normalised to sum to one, and pieces under 10⁻¹² of the cell are dropped* — the pieces'
  lengths are differences of cumulative fractions, and their rounding left a cell's outflow distributed 5.5·10⁻¹⁴ short
  or over, which grew to 7.6·10⁻¹³ of the mass over ten sub-steps; the slivers are the rounding of a boundary the two
  partitions share — if wrong, a piece of 10⁻¹² of a cell is moved to its neighbour.
- *The diagonal is the sum of the arrivals made* — every cell sends exactly what its arrivals take, so the books close
  whatever the pieces — if wrong, nothing: it equals dt Σc to round-off.
- *The momentum is projected at every arrival* (above) — the spec's text, and, at high Courant numbers, nearer the
  parallel transport it stands for — if wrong, 0.30 s of the worst case's 0.95 s is spent for nothing.
- *The regrid is optional* (`regrid=True` by default) — 6c calls the transport before the sub-step's own regrid and
  passes False, so the cells are split and merged once — if wrong, a sub-step regrids twice.

What the tests measure (2026-10-21): a pulse half again as thick as the liquid around it and 50 K warmer, 3 mm wide, on
1 mm columns at 0.1 m/s for 0.1 s, keeps its mass and heat, its centroids move the 10 mm to 6·10⁻⁸, 10⁻⁵ and
8.3·10⁻⁴ mm at Courant numbers 0.5, 1 and 2 (the scheme's long upwind tail reaching the grid's ends), and its variance
grows by the one-dimensional implicit upwind scheme's h² ν (1 + ν) a sub-step, less 4.2, 3.2 and 2.4 % — the numerical
diffusion u h (1 + ν)/2; depth fractions and both arrival rules exactly; on the 4 mm sphere's mixed skins every book to
4·10⁻¹⁵ of its scale, the wall taking 2.5 % of the largest column's momentum; a 10 s sub-step (Courant numbers to
875) positive, within the old enthalpies and speeds, conservative and regridded; on the sphere a velocity carried into
three neighbours 4.5–4.9° away loses 2.3–7.7 % of its momentum to the wall.

- [ ] **Step 1: The faces' test** — append to `tests/test_reentry_model_columns.py`:

{{DIFF_TEST_COLUMNS_T6A}}

- [ ] **Step 2: The transport's tests** — append to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T6A}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 11 failed — `AttributeError: 'Columns' object has no attribute
  'face_vector'` (5) and `AttributeError: module 'reentry_model.skin' has no attribute 'transport'` (6).

- [ ] **Step 4: The faces' vectors** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T6A}}

- [ ] **Step 5: The transport** in `reentry_model/skin.py`:

{{DIFF_SKIN_T6A}}

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_columns.py` 10 passed and
  `tests/test_reentry_model_skin.py` 56 passed; the unit tier 450 passed, 1 skipped. Nothing calls `transport` yet: the flow
  joins the skin body's sub-steps in 6c, so no result changes.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/columns.py tests/test_reentry_model_columns.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the metal moves sideways between neighbouring columns at the cells' velocities, linearly implicit and upwind, by fraction of the flowing depth and on top of a rigid top, each arrival's velocity projected onto its receiver's plane (sub-plan 21 Task 6a)"
```

---

### Task 6b: The thickness terms implicit in the column heights

**Files:**
- Modify: `reentry_model/columns.py` (`Columns.film_operator`, `Columns.film_curvature`), `reentry_model/skin.py`
  (`MomentumBook.response` and `.wall_response`; `momentum_step(..., response=False)`, each floating crust one unknown,
  the unused `P0` removed; `push(..., wet=None)` with the capillary term in its thin-film form; `transport(...,
  faces=None)`; `THICKNESS_PASSES`, `THICKNESS_TOL`, `flowing_tops`, `_walled`, `_capillary`, `FacePush`,
  `ThicknessBook`, `thickness_step`)
- Test: `tests/test_reentry_model_columns.py`, `tests/test_reentry_model_skin.py` (and Task 5's two capillary tests,
  re-measured)

**Interfaces:**
- Consumes: Task 2's `Columns.gradient`, `.curvature` (on the base alone now), `.top_points`, `.incidence`,
  `.conduct_weight`; Task 4's `momentum_step`, `anchored_cells`, `cell_viscosity`; Task 5's `push`; 6a's `transport` and
  `Columns.face_vector`.
- Produces: `Columns.film_operator(thickness) -> (n, n)` sparse [1/m2], minus the cotangent Laplacian of the top surface
  over each column's area; `Columns.film_curvature(thickness) -> (n,)` [1/m], the base's curvature plus that operator
  times the thicknesses; `skin.flowing_tops(sk) -> (n,)` bool, where a column's top cell flows; `push(columns,
  thickness, dp_ds, g, v_hat, rho_l, sigma, wet=None)`, its capillary term −σ ∇_s κ with κ the film's curvature, zero
  where `wet` is False; `momentum_step(sk, dt, push, shear, A_top, response=False)`, with `response` its book carrying
  `response` (n, K) [m/s per N/m3], each cell's velocity per unit push, and `wall_response` (n,) [m3], the anchored
  metal's drag on it; `transport(sk, columns, dt, regrid=True, faces=None)`, `faces` a `FacePush`;
  `skin.thickness_step(sk, columns, dt, book, applied, g, rho_l, sigma) -> ThicknessBook(push, wall, heights, faces,
  growth, passes)`, push and wall (n, 3) [N s], heights (n,) [m] the thicknesses predicted at the end of the transport,
  faces a `FacePush(response (n, K), across (2E,) [Pa])` for the transport, growth the largest dt·s of the
  Rayleigh–Taylor growth over the columns, passes the linearisations made; `skin.THICKNESS_PASSES` (8),
  `skin.THICKNESS_TOL` (10⁻³).
- **What 6c calls, in this order, every flow sub-step:**
  1. `wet = flowing_tops(sk)`; `G, terms = push(columns, sk.thickness(), dp_ds, g, v_hat, rho_l, sigma, wet=wet)`;
  2. `mbook = momentum_step(sk, dt, G, shear, A_top, response=True)`;
  3. `tbook = thickness_step(sk, columns, dt, mbook, terms['thickness'] + terms['capillary'], g, rho_l, sigma)`;
  4. `trbook = transport(sk, columns, dt, regrid=False, faces=tbook.faces)`;
  5. the sub-step's own regrid.
  The column's momentum changes by `mbook.push + mbook.shear + mbook.wall + tbook.push + tbook.wall + trbook.momentum +
  trbook.wall`. `tbook.growth` must stay well under one: backward Euler follows a Rayleigh–Taylor growth only where the
  sub-step resolves it.

**Decided before this task was written (Asha, 2026-10-21; spec §6.4, §14 assumption 1, §16).** A first attempt stopped
on a finding: Task 5's capillary term made surface tension run backwards. It took κ as Task 2's curvature of the top
surface, and Task 2 gives a corner touching an obtuse triangle its neighbours' mean curvature. Such a corner cannot see
its own height: raising it makes its neighbours concave and lowers its own curvature, so a pile-up lowers its pressure
and draws metal in. Raising each of the test sphere's 33 such corners by 1 µm lowered its own curvature at 32 (by up to
0.083 1/m) and left it unchanged at the other, where every other corner's rose by at least 0.21 1/m.

The faces' exchange times that term's linearisation has eigenvalues down to −1.48·10¹⁰ 1/m⁴ on the test sphere and
−2.40·10¹¹ 1/m⁴ on the 2 mm sphere (48 such corners). That is a growth of about 80 and 1 300 1/s for 0.3 mm of liquid,
unstable at any sub-step and in any time scheme. A uniform 0.3 mm liquid layer with 1 µm of noise on the test sphere,
its velocities reset after each 10 ms sub-step, roughened from 0.48 to 4.0 µm in 30 sub-steps.

Asha's decisions:
- **(option A)** Surface tension's pressure takes the thin-film form on the top surface's geometry. The pressure is
  (1/A_i) dE/dH_i, where E is the thin film's surface energy on the top surface and A_i is the column's area, the one
  that turns the column's volume into its height in the skins.
- Surface tension acts only where a column's top cell flows.
- The depth-uniform pressure stays, stated as a limitation.

Option B, the top's full curvature over each column's own area with no obtuse rule, was also stable but left a uniform
layer rippled by the mesh: 7.7 µm on the test sphere. Keeping Task 5's term, explicitly or with the energy-consistent
increment, stayed unstable.

**The thin film's curvature.** κ_i = κ_base,i + (J H)_i, where:
- κ_base is Task 2's curvature of the base, its obtuse-corner rule acting on the fixed base alone;
- J = −L_top / A is minus the cotangent Laplacian of the top surface over each column's area (`Columns.film_operator`).

A·J is symmetric and non-negative, and so is the faces' exchange B diag(w) Bᵀ. Their product therefore has no growing
mode: on the test sphere with a rough top (0–3.4 mm) every eigenvalue is real and non-negative to round-off. A cosine's
curvature on a plane is the compact second difference of the heights (to 4·10⁻⁷, from its own top's cotangents). A
uniform layer's on the sphere is the base's exactly, so a uniform 0.8 mm skin on the test sphere now gets at most
0.32 Pa/m of capillary push (0.031 on average) from the mesh, where the top's own curvature gave 11 Pa/m (0.28).

A cotangent depends on angles alone, so the top's cotangents differ from the base's only where the thickness varies:
the two operators on a smooth bump at 3.4 mm differ by 0.35 % on average on the test sphere. With the base's area, which
the energy argument requires, a thick layer's bump curvature stays about 2H/R above the top surface's own (14 % at
3.4 mm on the 50 mm sphere): the columns hold their volume over the base's area. Asha kept the prism columns and
stated the error (spec §14 assumption 10, 2026-10-22), over wedge columns whose area grows with height.

**Where surface tension acts.** Surface tension acts only where a column's top cell flows (`flowing_tops`). An anchored
rigid top and a floating crust feel no capillary push, and no capillary push crosses a face with a non-flowing top on
either side. A rigid neighbour's height enters a flowing column's curvature through the top's Laplacian, as a wall's
would. In the column push's gradient, a wall takes its flowing neighbours' mean curvature (`_walled`), so it pushes
nothing of its own. The energy argument holds with the wall's heights as data: the capillary flux runs only between
flowing tops.

**The scheme** (spec §6.4's last paragraph: "inside the transport solve", read since 2026-10-21 as within the sub-step's
transport step, so the thickness terms see the heights at the end of the transport). The momentum step is linear in
the push for fixed viscosities. With `response` it also solves its last pass's system for a unit push alone: one more
tridiagonal solve, giving each cell's response r and the anchored metal's drag on it.

The two terms are the gradient of P = σκ − ρ_l (g·n) H. Across each dual face P pushes by its two-point difference,
φ = w (P_a − P_b), where:
- w is the face's clipped cotangent weight (Task 1's);
- the capillary part acts only between flowing tops;
- the deceleration's part uses the mean of the two columns' g·n.

φ is linear in the heights through J. Each sub-step then does three things.

1. **It predicts the heights at the end of the transport, H + dH**, from the transport linearised on the columns. A
   column sends across a face what its cells carry, at their velocity without the thickness push plus their response
   times the face's push. This is implicit in its moving metal, as the transport is, and linearised in the face's push
   through the mobility Σ r dz of the cells that cross. Newton's method on that piecewise linear flux re-solves at the
   new face pushes until the flux the last solve missed would move no column by more than THICKNESS_TOL of the largest
   change. It took 1–3 solves on the test sphere and 3 on the melting state.
2. **It corrects every moving cell's velocity** by its response times the change of `push`'s two terms at the predicted
   heights. The cells' momentum follows the push of Task 5's form, and the impulse and the wall's share are booked.
3. **It hands `transport` each face's own push** (`FacePush`). A cell's velocity across a face becomes its velocity
   without the thickness push plus its response times that face's push.

Positivity is 6a's: every rate stays clipped at zero, so the masses stay positive at any velocity.

**A floating crust is one unknown.** In `momentum_step` each run of crust cells is now one unknown of the tridiagonal
system, a rigid body: its mass is the cells', and its drags are those at its top and bottom faces.
- Task 4 coupled the crust's cells by 10⁹ Pa s. Those rows, some 10¹² times the cells' inertia, amplified the solver's
  round-off: on mixed sphere skins the books closed to 1.2·10⁻⁵ of the largest impulse. They now close to
  6.6·10⁻¹⁴ (1 ms) and 6.2·10⁻¹³ (10 ms).
- The velocities moved by up to 1.1·10⁻³ of the largest at 10 ms. This is the coupled crust's own round-off: it moved
  by 1.1·10⁻⁴, 1.1·10⁻⁵ and 9.6·10⁻⁷ at 10⁸, 10⁷ and 10⁶ Pa s, and by 4.2·10⁻⁷ and 4.2·10⁻⁶ at 10⁵ and 10⁴, where the
  crust's flexibility shows.
- On Task 4's own crust column the change is 9.8·10⁻⁸ of the top speed, and its books close to 8.9·10⁻¹⁵ instead of
  9.9·10⁻⁸.
- Columns without a crust are unchanged to the last bit.

**Measured** (2026-10-21, `prototype/work-2026-10-18-column-flow/measure/task6b/`; the tests' numbers).

*Rayleigh–Taylor growth* (spec §11.1 level 1). The layer is 1 mm of slurry at 902 K (30 % solid, 2.68 Pa s at rest)
on a strip of 1 mm columns two corners wide, pulled out of its surface at 100 m/s²: a capillary cutoff of 11.5 mm. This
is the lubrication regime the formula assumes, the viscous time ρh²/μ being 0.9 ms against a growth time of 0.63 s.
- A cosine 24 mm long grows at the thin-film rate h³k²(ρg − σk²)/(3μ) = 1.578 1/s to within −0.41, −0.22 and +0.42 % at
  sub-steps of 1, 2.5 and 10 ms. The grid's compact difference accounts for −0.40 %, and backward Euler adds +dt·s/2.
- A cosine 8 mm long, below the cutoff, decays at −19.44 1/s to within 1.2 % at 1 ms. At 10 ms it decays 8.7 % slower,
  backward Euler's own damping, whose rate it meets within 0.10 % (0.26 % at 1 ms).

*The two-column mode.* 0.3 mm of liquid on 2 mm columns, alternating by 1 µm, falls to the following fraction of itself
in one sub-step, and never above 0.009 of itself after that:

| Sub-step | Under 30 m/s² | Without the pull |
|---|---|---|
| 10 ms | 0.112 | 0.117 |
| 100 ms | 5.5·10⁻⁴ | −6.3·10⁻⁵ |
| 1 s | −0.051 | −0.056 |

The table shows the grid-scale mode damped within one sub-step at every sub-step from 10 ms to 1 s, with or without the
pull.

*The faces' own push against the columns' alone.* A bump 1 µm high on one column of 1.1 mm of slurry relaxes to
12.5 % in one 10 ms sub-step and 3.0 % in five. With the columns' pushes alone it stays at 99.8 %: they are centred
differences, zero on the bump's column by symmetry, and the upwind transport takes metal out of a column only at its
own velocity.

*A uniform layer.* 0.3 mm of liquid on the test sphere stays uniform to 26 nm, with speeds of at most 2.2 µm/s, over 20
sub-steps of 10 ms.

*Prediction and Newton passes.* With 1 µm of roughness on 0.3 mm of liquid on the test sphere:
- the transport makes the predicted heights to 0.62 % of the largest change;
- the Newton passes were 1, 3, 3 and 3;
- the roughness falls to 34 % in four sub-steps.

*Books.* On mixed sphere skins, every momentum book closes to 7.8·10⁻¹⁴ (1 ms) and 8.8·10⁻¹³ (10 ms) of the largest
impulse, and mass and enthalpy are exact. With σ = 0 and g = 0 the whole step equals the momentum step followed by the
transport, to the last bit.

*Cost* on the 100 mm sphere's melting state. The melted face is the one where g·n > 0 under the flight direction given
to `push`: its top 0.3 mm is liquid, with 141 000 cells moving at 0.1 m/s along the meridian, under 30 m/s² and
σ = 0.8 N/m. Each figure was measured three times.
- The thickness step takes 0.29–0.30 s at 3 passes.
- The response adds 0.03–0.05 s to the momentum step's 0.24–0.25 s.
- The faces' push adds nothing measurable to the transport's 0.31–0.37 s.
- Together 6b adds about 0.33 s a sub-step.
- At the melt's rim, the transport made the predicted heights to 6.5 µm of a change of 0.15–0.17 mm.

**Rulings** (2026-10-21; each is what, then why, then the cost if wrong):
- *Correction form.* The momentum step is given `push`'s whole G, and the thickness step adds the change. Why: the
  slurry's shear-rate fixed point then sees the whole push, and the curvature's own nonlinearity enters at the start of
  the sub-step. Cost if wrong: none for fixed viscosities, as both forms give the same velocities.
- *The faces' own push in the transport, the stored momentum from the column push* (the coordinator's ruling, the
  column form of Rhie–Chow interpolation). Why: the bump above. Cost if wrong: the mass and the momentum are carried at
  slightly different velocities at the grid scale.
- *The heights predicted from the transport linearised on the columns, by Newton's method on the crossing cells, stopped
  on the flux missed* (the coordinator's ruling). Why: the transport is implicit in the cells' masses, which a column
  model of it can follow but not reproduce. Cost if wrong: where the cells of a column carry very different shares, the
  prediction differs from the transport by a few per cent of the change, and the thickness terms are that much
  explicit.
- *The face's push by Task 1's clipped cotangent weights, the deceleration's term with the two columns' mean g·n.* Why:
  the energy argument needs a non-negative exchange, and the spec's form is ρ_l (g·n) ∇H. Cost if wrong: 4 of the
  28 245 faces on the 2 mm sphere are clipped.
- *A wall takes its flowing neighbours' mean curvature in the column push's gradient.* Why: the face carries no
  capillary push into the wall, and the column's mean of its face velocities is then about half the push from its
  liquid side. Cost if wrong: through the cells' inertia, a velocity toward a wall carries metal onto its top a sub-step
  later.
- *A crust as one unknown* (the coordinator's ruling). Why: section above. Cost if wrong: none, as it is the rigid limit
  of Task 4's coupling.
- *The tests on strips two corners wide* (the coordinator's ruling). Why: with the film's curvature an open edge is a
  mirror for surface tension. Cost if wrong: the strips' ends disturb the cosines, so the tests measure the middle wave.

- [ ] **Step 1: The film curvature's tests** — append to `tests/test_reentry_model_columns.py`:

{{DIFF_TEST_COLUMNS_T6B}}

- [ ] **Step 2: The thickness terms' tests**, and Task 5's two capillary tests re-measured, in
  `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T6B}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 15 failed:
  - `AttributeError: 'Columns' object has no attribute 'film_curvature'` (2);
  - `AttributeError: 'Columns' object has no attribute 'film_operator'` (1);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'flowing_tops'` (9);
  - `TypeError: momentum_step() got an unexpected keyword argument 'response'` (1);
  - the crust's books, `assert ... 1.32e-05 < 1.5e-12` (1);
  - Task 5's sphere test, `assert (... 10.99 < 0.7)` (1).

- [ ] **Step 4: The thin film's curvature** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T6B}}

- [ ] **Step 5: The momentum's response, each crust one unknown, the push's thin-film capillary term, the faces' push in
  the transport and the thickness step** in `reentry_model/skin.py`:

{{DIFF_SKIN_T6B}}

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_columns.py` 13 passed and
  `tests/test_reentry_model_skin.py` 67 passed; the unit tier 464 passed, 1 skipped. Nothing calls `thickness_step`
  yet: the flow joins the skin body's sub-steps in 6c, so no flight's result changes.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/columns.py tests/test_reentry_model_columns.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the push's thickness terms implicit in the column heights at the end of the transport -- surface tension in its thin-film form, only where a column's top flows, each face's own push in the transport and the cells' velocities corrected by their response -- and each floating crust one unknown of the momentum step (sub-plan 21 Task 6b)"
```

---

### Task 6c: The flow in the skin body

**Files:**
- Modify: `reentry_model/skin.py` (`FLOW_CYCLES`, `FLOW_SECONDS`, `Flow`, `FlowBook`, `flow_substeps`;
  `SkinBook.momentum_out` and `.flow`; `advance(..., flow=None)`), `reentry_model/skinbody.py` (`FLOW_DEFAULT`,
  `FlowForcing`; `SkinBody(..., flow, flow_cycles, surface_flow)`, `flow_forcing`, `advance(..., forcing=None)`, the
  momentum books and `momentum_residual`, the flow's history columns), `reentry_model/coupled.py` (`SKINBODY_COLUMNS`,
  `skin_results`, the frames' `top_speed`), `reentry_model/cli.py` (`--skin-flow`, `--skin-flow-cycles`,
  `skin_name_suffix(..., flow, flow_cycles)`, their checks, `build_thermal`)
- Test: `tests/test_reentry_model_skin.py`, `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_cli.py`
  (the skin run's settings gain the flow's two keys)

**Interfaces:**
- Consumes: Task 5's `push` and `along_meridian`, `SurfaceFlowResult.dp_ds` and `.tau`; Task 4's `momentum_step`; 6a's
  `transport`; 6b's `flowing_tops`, `momentum_step(..., response=True)`, `thickness_step` and `FacePush`;
  `Skins.momentum_out`; `SurfaceFlow.evaluate`.
- Produces: `skin.Flow(columns, dp_ds, shear, g, v_hat, rho_l, sigma, A_top, cycles=FLOW_CYCLES)`, what drives the flow
  over a part (dp_ds (n,) [Pa/m], shear (n, 3) [Pa], g (3,) [m/s2], A_top (n,) [m2]); `skin.flow_substeps(sk, flow,
  dt) -> FlowBook`; `FlowBook(push, shear, wall, carried, courant, growth, passes, sub_steps)`, the impulses (n, 3)
  [N s] and the carried momentum [kg m/s], with `.zero(n)` and `.add(other)`; `skin.FLOW_CYCLES` (10) and the
  diagnostic `skin.FLOW_SECONDS`; `skin.advance(..., flow=None)`, its book gaining `momentum_out` (n, 3) [kg m/s], what
  the removed liquid took, and `flow`, the sub-steps' FlowBook; `skinbody.FlowForcing(dp_ds, tau, g)` with
  `.toward(other, frac)`; `SkinBody(..., flow=FLOW_DEFAULT, flow_cycles=None, surface_flow=None)` with `.flow_on`,
  `.flow_cycles`, `.flow_forcing(state)`, `.advance(t, dt, loads, state=None, forcing=None)`, `.momentum_residual()`,
  `.impulse`, `.momentum_lost`, `.impulse_size`, `.last_flow`; `skinbody.FLOW_DEFAULT` (False); history columns
  `momentum_residual`, `speed_max_ms`, `roughness_max`, `flow_courant_max`, `flow_growth_max`, `flow_passes_max`,
  `wall_flow_s`; results `momentum_residual`, `flow_courant_max`; the frames' `top_speed`; `run --skin-flow on|off`
  (run-name suffix `_skinflow-on`) and `--skin-flow-cycles <n>` (`_flowcycles<n>` when not FLOW_CYCLES).

The flow joins the skin body. Every macro step, the surface flow is evaluated on the corners where the heating is —
from the step's aero state, the top surface's angles to the flight direction, the nose radius and the tops'
temperatures at the step's start — for its pressure gradient along the meridian (Task 5's `dp_ds`) and its gas shear
`tau`; with the deceleration, which the metal feels as g, the drag's deceleration along the flight direction (minus the
body's own acceleration, Task 5), that is the step's `FlowForcing`. Each part takes it interpolated linearly in time
between the last macro step's and this one's at the part's middle, as it takes the heating (spec 7.1), the shear laid
along the meridian in each column's plane. Every sub-step of both passes then runs, after its conduction, `cycles` flow
sub-steps of its length over `cycles`, each in 6b's order — `push` with surface tension where a top flows,
`momentum_step` with its response, `thickness_step`, `transport` with the faces' push and no regrid — the cells split
and merged between them; the removal (with `remove`), the draw and the sub-step's own regrid follow the last (spec 7.1:
heat, flow, sideways exchange and spray). The trial pass runs the flow as the real pass does (spec 7.2's sub-steps are
spec 7.1's), the sensitivity carried by the transport. Where no cell flows nothing is computed: a velocity left in rigid
cells is booked to the wall and cleared, as the momentum step would.

**The momentum books.** The flow's impulses — the momentum and thickness steps' push, the gas shear, the walls (the
anchored metal's drag, the momentum cells held on joining it, the transport's projections) and what the transport
carried, which sums to zero over the body — add up in each part's book; what leaves with removed liquid, handed-down
cells or a reconciliation at the bottom is booked as lost; metal drawn arrives at rest. `momentum_residual` is the
skins' momentum less its value at time zero, the impulses and the losses, over the impulses' total size: zero for exact
books.

**The flow's sub-step** (option c, Asha, 2026-10-21, chosen here from measurements;
`prototype/work-2026-10-18-column-flow/measure/task6c/`, 2026-10-22, the laptop under a load average of 5 to 140 from
other sessions, so the times are processor times and indicative). The 100 mm sphere's 9 417 columns of 10 µm cells, its
front hemisphere's top 0.3 mm liquid at 950 K over rigid metal at 850 K (6b's destabilising face), at rest at t = 0 and
driven by the surface flow of the 100 mm flight at 93 s — a stagnation pressure of 11.9 kPa, -dp/ds up to
3.4·10⁵ Pa/m at 45°, gas shear up to 87 Pa at 49°, 29.8 m/s² of deceleration — with surface tension, for 50 ms of
flow alone (no conduction), at flow sub-steps of 10 ms to 0.25 ms:

| Flow sub-step | Courant, largest (median) | KE the transport destroys / KE made | KE at 50 ms / 0.25 ms run | Thickness change, L1 error / 0.25 ms run | Volume past 105° / 0.25 ms run | Linearisations, most |
|---|---|---|---|---|---|---|
| 10 ms | 1 266 (244) | 0.59 | 4.81 | 95 % | 0 | 8 (the cap) |
| 5 ms | 26 (23) | 0.62 | 0.31 | 85 % | 0.11 | 7 |
| 2.5 ms | 12.6 (11) | 0.52 | 0.36 | 59 % | 0.89 | 5 |
| 2 ms | 9.9 (8.4) | 0.43 | 0.48 | 44 % | 0.96 | 3 |
| 1 ms | 4.8 (4.1) | 0.27 | 1.04 | 3.9 % | 1.007 | 4 |
| 0.5 ms | 2.4 (2.0) | 0.28 | 1.02 | 1.4 % | 1.005 | 2 |
| 0.25 ms | 1.2 (1.0) | 0.27 | 1 | — | 1 | 2 |

The table shows the result converging at 1 ms and not above it. In 50 ms the liquid drains off the front at a few metres
a second and crosses the shoulder onto the rigid rear: the volume that crossed 85° is within 1.1 % of the 0.25 ms run's
(4 589 mm³) at 5 ms and below (3 % short at 10 ms), but the runoff over the rear is not. At 1, 0.5 and 0.25 ms it
reaches the rear stagnation point, every band of 15° from 90° to 180° within 7 % (1 ms) and 4 % (0.5 ms) of the finest
run's; at 2, 2.5, 5 and 10 ms it stops at 150°, 135°, 120° and 105°. The cause is 6a's arrival rule, not the scheme's
smearing: metal landing on a rigid top becomes a new top cell after the solve, so it flows on only in the next flow
sub-step, and a front crossing rigid metal advances at most one column a flow sub-step — 0.2 m/s on 2 mm columns at 10
ms, 2 m/s at 1 ms. That predicts the fronts seen: 25, 20, 10 and 5 columns of 2 mm in 50 ms are 57°, 46°, 23° and 11° of
the 50 mm radius beyond the shoulder at 2, 2.5, 5 and 10 ms. Inside the flowing layer the scheme's own smearing at 1 ms
is already at its floor: the transport destroys 27 % of the kinetic energy the push and the shear make, the same at 0.5
and 0.25 ms, so that is the 2 mm grid's upwind diffusion and the projections between planes, not the time step's. At 10
ms the thickness step's linearisations stop at their cap of 8 with the predicted heights 13 times the change off, and
its correction adds kinetic energy: the fastest cell is at 22 m/s at 30 ms against 4.9 at 1 ms. With every column's top
melted (no front; 9 417 columns, 282 000 cells moving), 10 ms flow sub-steps blow up — the fastest cell 72 m/s at the
third sub-step, 2 563 m/s at the fourth, the predicted heights 16 and 91 times the change off — and 5 ms keeps the
linearisations at their cap with the predictions about one change off; at 1 ms it is stable (Courant numbers to 4.7, two
linearisations, the predictions within 0.24 of the change over 10 ms), but its transport took 35-114 s a sub-step, the
flows meeting at the rear stagnation point joining 18 000-23 000 cells into one cycle (below).

*Rayleigh-Taylor growth* (`ThicknessBook.growth`, dt s, one flow sub-step from rest on the sphere's front melted 0.3,
1.0 and 3.4 mm deep, the flowing depths of the element model's flight at about 85, 105 and 120 s, under the forcing at
93 s, 29.8 m/s², and 115 s, 48.4 m/s²):

| Liquid depth, deceleration | 10 ms | 5 ms | 2 ms | 1 ms |
|---|---|---|---|---|
| 0.3 mm, 29.8 m/s² | 0.018 | 0.0048 | 8.3·10⁻⁴ | 2.1·10⁻⁴ |
| 0.3 mm, 48.4 m/s² | 0.047 | 0.013 | 0.0022 | 5.6·10⁻⁴ |
| 1.0 mm, 29.8 m/s² | 0.072 | 0.018 | 0.0030 | 7.6·10⁻⁴ |
| 1.0 mm, 48.4 m/s² | 0.19 | 0.049 | 0.0079 | 0.0020 |
| 3.4 mm, 29.8 m/s² | 0.25 | 0.063 | 0.010 | 0.0025 |
| 3.4 mm, 48.4 m/s² | 0.67 | 0.17 | 0.027 | 0.0067 |

The table shows the growth falling as the square of the sub-step — the layer's response in a sub-step is limited by
its inertia — and well under one at 5 ms and below in every case; at 10 ms the deep layers late in flight reach 0.67,
where backward Euler's growth in a sub-step, 1/(1 - dt s) = 3.0, is half again the exact e^(dt s) = 1.95.

**Ruling: ten flow sub-steps a heat sub-step (`FLOW_CYCLES` = 10, 1 ms under the 10 ms heat sub-step)** — at 1 ms the
measured state is within 3.9 % of the 0.25 ms run in where the metal went and 4 % in its kinetic energy, the runoff
reaches where the finest run's does, the growth is at most 0.0067 even for 3.4 mm of liquid at 48 m/s², the Courant
numbers are at most 4.8 and the kinetic energy the transport destroys is the grid's own 27 %; 2 ms holds the runoff back
by up to 30° of the rear and loses half the kinetic energy by 50 ms, and 10 ms is unstable in a fully melted layer.
Accepted: Courant numbers to about 5, where the implicit upwind scheme's time smearing no longer shows against the
grid's, and the grid's 27 % (the 2 mm surface's; only finer columns reduce it) — if wrong, the flow costs ten times what
one sub-step would: a 10 ms heat sub-step of the heated melting state takes 48 s with ten against 8 s with one and 4.5 s
for its conduction alone (*Cost*, below). The sub-step is fixed, not adapted: the Courant numbers grow with the layer's
speed, which grows with its depth and the push late in flight, and a front faster than 2 m/s on 2 mm columns is held
back at 1 ms too (Review Focus 13 and 14).

**Ruling: the flow is off by default (`FLOW_DEFAULT`, `--skin-flow off`), and with it on the liquid stays** — stage 2's
skins remove the liquid as it forms "instead of flow and spray" (spec 11.2); with the flow on and no spraying until
stage 4, the liquid stays and flows (`remove` False), and a flight would only pile it up — if wrong, a stage-3 flight
must ask for the flow. Flow off is the heat-only skins: the skin body's tests pass unchanged, and the 3 MW/m² removal
run of `test_heating_with_removal_keeps_mass_and_energy_exact` gives Task 6b's code's every array to the last bit
(2026-10-22, numpy's generator seeded as a run seeds it; unseeded, pyamg's starting vectors make any two runs differ in
round-off).

**Ruling: on the command line the flow takes Girin's removal** — `--skin-flow on` needs `--surface-model skin` and
`--melt on`, and `--removal girin` (the default): the liquid stays, and sub-plan 22's spraying will remove it on the
same command line; `--removal instant` is refused with the flow on, being the heat-only skins' stand-in for flow and
spray (spec 11.2), as `--removal girin` is refused with it off. `--skin-flow-cycles` needs the flow on and is an integer
of at least 1; the run name gains `_skinflow-on`, and `_flowcycles<n>` only for other than ten (spec 10.1's convention)
— if wrong, a stage-3 run with the flow on takes the girin removal's continuum step (`--dt-continuum`, 0.0125 s once
the body's Knudsen number is below 0.01) and its name says `_melt-girin` although nothing sprays until sub-plan 22;
`--dt-continuum off` keeps the 0.5 s step.

**Ruling: the forcing interpolated in time as the heating is, the shear in each column's plane** — spec 7.1's
treatment of the heating, at no cost; the push's density is the skins' liquid density (`Skins.rho_l`, the solid's
with `--skin-swell off`, as the cells' masses have it) — if wrong, the push differs by the forcing's change over half a
macro step.

**Ruling: the trial pass runs the flow** (spec 7.1 and 7.2: a sub-step is heat, flow, sideways exchange and spray) — it
changes the base law: on the 100 mm sphere's heated melting state (`heated_state`: the front's skins at the conduction
profile of 3 MW/m², 0.30 mm of liquid over 0.55 mm of slurry, the base held at the profile's 874 K, the flow developed
over 20 ms under the 93 s forcing), one part of 0.0625 s with the flow and without it gives base laws that differ, at
the base temperature, by 6.0 % of the energy through the columns' tops and bases in total; 3 428 of the 9 417 columns
differ by more than the interface tolerance of 10⁻³ of their own (the median of the flowing columns 7·10⁻⁴), and the
slope by up to 46 %. A trial pass without the flow would leave the real pass's base heat that far from its law: repeated
by the damped repetition, and what it leaves deferred to the nodes. The real pass's base heat meets the law made with
the flow to 9·10⁻¹⁵. It costs what the flow costs: the trial pass took 1 216 s of processor time with ten flow sub-steps
a heat sub-step against 17 s without the flow. If wrong, half the flow's cost is spent on the law.

**Ruling: the momentum residual over the impulses' total size**, the push's, the shear's and the walls' summed in size
over every column and part since time zero, as the energy residual is over the heat absorbed — if wrong, only the
scale of a number that should be round-off.

*Cost* (processor time, on `heated_state` with 310 000 cells flowing, the flow developed over 20 ms; each measured twice,
under other sessions' load, so indicative). A 10 ms heat sub-step:

| Flow sub-steps | Whole | Conduction | Push | Momentum | Thickness | Transport | Regrid | Rest (draw, regrid) |
|---|---|---|---|---|---|---|---|---|
| 1 (10 ms) | 7.9-8.5 s | 4.1-4.6 | 0.05 | 1.0 | 1.0 | 1.4 | — | 0.4 |
| 10 (1 ms) | 47.8-48.3 s | 4.5 | 0.5 | 9.4 | 3.2-3.5 | 28-29 | 1.4-1.5 | 0.3 |

The table shows the transport taking three fifths of a heat sub-step at ten flow sub-steps, 2.9 s each, against 1.4 s
for the one 10 ms sub-step, whose Courant number of 490 and eight linearisations (the cap) mark it unresolved; the
momentum step's shear-rate fixed point runs its five passes every sub-step, the slurry where liquid meets rigid metal. A
macro step while melting is eight parts of seven sub-steps, each run by the trial and the real pass, 1.8 times on
average (fact 126): about 200 heat sub-steps, which on this state is about 2.7 h with ten flow sub-steps, 27 minutes
with one and 16 minutes with the flow off — the flow, not the heat, is now the skins' cost. (Fact 126 measured a melting
skin pass at 0.65 s a 10 ms sub-step; these skins are deeper, 1.6 mm, and the laptop was loaded, so the ratios are the
finding, not the seconds.) The remedy for the transport's large cycles below would take its share down; the momentum
step's solve over every column, flowing or not, is the next.

*A quiescent melted region* (Review Focus 12). The sphere's melting state at rest under surface tension alone, no gas
forcing and no deceleration: the melt's rim at the shoulder rounds at 5-17 mm/s (Courant numbers 0.01-0.2), and the
transport's cycles of meeting flows grow from 7 cells to 10 315 of the 141 000 moving within 12 flow sub-steps of 1 ms
(to 7 022 within 6 of 10 ms), its time from 1.1 to 4.0-4.8 s a sub-step. Profiled, 2.3 s of a 3.35 s transport is
SuperLU's factorisation of its two systems, whose lower triangles fill 5.3 and 5.8 times. Under the 93 s forcing the
largest cycle was 1 700 cells after 20 sub-steps of 1 ms and the factorisations took 0.06 and 0.15 s, filling 1.1 times;
with every column melted, flows meeting at the rear stagnation point joined 57 000 cells into one cycle within 9
sub-steps of 5 ms, and a transport took 651 s. *A remedy measured, not built:* solve each system by the iteration x ← x
+ L⁻¹ (b − A x) started from the upwind solve, L its lower triangle in upwind order — the whole system but the cycles'
back-couplings, which an upwind order makes triangular with no fill. On the quiescent state it took 2-3 corrections to
reach 10⁻¹⁰-10⁻¹³ of the right-hand side, 0.05 s against the LU's 0.30 s for the mass and 0.16 s against 1.23 s for the
momentum, the masses positive and the books missing 2·10⁻¹⁷ of the total at 10⁻¹³; under the forcing (small cycles) 4-5
corrections, 0.10 s against the LU's 0.06 s; with every column melted (a cycle of 18 379 cells) 7-9 corrections, 0.19 s
against 7.4 s and 0.41 s against 41.8 s. Installing each cell's mass from the iterate's fluxes, m + dt (arrivals −
departures), would close the books exactly by construction. It changes 6a's solver, so it is reported, not built (Review
Focus 12).

**Measured against the tests** (2026-10-22): one flow sub-step is the push, the momentum step, the thickness step and the
transport called by hand, to the last bit, at one cycle and three, the books to 2.1·10⁻¹² and 9.4·10⁻¹³ of the largest
impulse (6b's on these skins); four flow sub-steps in each 10 ms heat sub-step give sixteen heat sub-steps of 2.5 ms
with one each, to the last bit, with the conduction off; a pass with the conduction on, removal and the flow closes its
momentum to 4.9·10⁻¹⁶; on the 4 mm test sphere under 3 MW/m² from 880 K for two 0.25 s steps with the flow on the mass
is kept to 1.1·10⁻¹⁶, the energy to 9.8·10⁻¹³ and the momentum to 2.8·10⁻¹⁷, the fastest cell 1.6 m/s, Courant numbers
to 6.9.

- [ ] **Step 1: The flow sub-steps' tests** — append to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T6C}}

- [ ] **Step 2: The skin body's tests** — append to `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T6C}}

- [ ] **Step 3: The command line's tests** — in `tests/test_reentry_model_cli.py` (the skin run's settings gain the
  flow's two keys):

{{DIFF_TEST_CLI_T6C}}

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py -q`. Expected: 11 failed, 143 passed —
  - `AttributeError: module 'reentry_model.skin' has no attribute 'Flow'` (4) and `... 'flow_substeps'` (1);
  - `TypeError: SkinBody.__init__() got an unexpected keyword argument 'flow'` (3);
  - `TypeError: skin_name_suffix() got an unexpected keyword argument 'flow'` (1);
  - `SystemExit: 2`, `unrecognized arguments: --skin-flow on --skin-flow-cycles 2` (1);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'FLOW_CYCLES'`, the skin run's settings (1).

  The seven exit-code cases pass already, their flags being unknown to the parser.

- [ ] **Step 5: The flow sub-steps and the pass's books** in `reentry_model/skin.py`:

{{DIFF_SKIN_T6C}}

- [ ] **Step 6: The flow in the skin body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T6C}}

- [ ] **Step 7: The history's flow columns, the results and the frames' top speed** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T6C}}

- [ ] **Step 8: The settings** in `reentry_model/cli.py`:

{{DIFF_CLI_T6C}}

- [ ] **Step 9: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 72 passed,
  `tests/test_reentry_model_skinbody.py` 16 passed (2 deselected, the reference tier's), `tests/test_reentry_model_cli.py`
  66 passed, `tests/test_reentry_model_coupled.py` 9 and `tests/test_reentry_model_columns.py` 13 as before (176 in
  860 s); the unit tier 481 passed, 1 skipped (1 015 s). The flow is off by default, so no flight's result changes:
  the heat-only skins reproduce Task 6b's to the last bit.

- [ ] **Step 10: Commit.**

```bash
git add reentry_model/skin.py reentry_model/skinbody.py reentry_model/coupled.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py
git commit -m "Skins: the flow joins the skin body -- the surface flow's pressure gradient and shear on the corners and the deceleration drive ten flow sub-steps of push, momentum, thickness and transport in every 1 ms of each heat sub-step, in the trial and real passes, the momentum books in the history; off by default until spraying, the liquid kept with it on (sub-plan 21 Task 6c)"
```

---

## Tasks still in outline (to be written in full)

6. **Sideways exchange**, split in three (Asha, 2026-10-21), each with its own tests and commit; 6a, the transport at
   given velocities, 6b, the thickness terms implicit, and 6c, the flow in the skin body, are written out above as
   Tasks 6a, 6b and 6c.
   - **6d. Making the flow affordable** (Asha, 2026-10-22, option a, before Tasks 7 and 8): 6c's flow cost about 48 s
     of processor time a melting 10 ms skin sub-step against 4.5 s for the heat, three fifths of it in the transport's
     solve, which slows sharply where flows meet. An iterative transport solve started from the upwind solve, its books
     exact by construction; metal landing on a rigid top flowing on within the same solve, lifting the one-column-a-flow-
     sub-step front limit, then the flow sub-step measured again; the number of flow sub-steps chosen each skin sub-step
     from the Courant number and the growth; the momentum step on flowing columns only. Then the cost again, the trial
     pass included.
7. **The droplet momentum sink** placeholder at the top cells for stage 4, tested with a prescribed release.
8. **Verification**: level 2's plate (spec §11.1, amended 2026-10-12: the box's 90-degree edges fold the surface; sub-plan
   20's `analysis/skin_plate.py`) with flow on and a push along the face (mass, momentum and energy exact; sensible
   behaviour), the momentum-balance residual and the roughness that 6c's history columns report read on it, and the
   plate's flow sub-step checked against 6c's choice; the cost of the sideways solve; and
   the share of the flight's flow at short wavelengths (kh not small) in flowing layers over about 1 mm, for the
   depth-uniform pressure's limitation (spec §14 assumption 1; Asha, 2026-10-21, option a). If that share is
   significant, a further task replaces the depth-uniform pressure by a non-hydrostatic pressure solved on the cells
   (Asha, 2026-10-22).

## Inputs to settle before the remaining tasks are written

- **The skins' interface** (sub-plan 20's Results; facts 118–125). Amended after the build (Asha, 2026-10-12): damped
  repetition, the profile valuation and the parts make the coupling converge and level 2 pass at 0.5 s, but a melting
  step on the 100 mm sphere costs up to 1 041 s (fact 125), and buried superheated cells remain at 6 MW/m² (fact 124).
  The onset excess at 0.1–0.25 s steps was cured by the melting part cap (fact 132); levels 2 and 2b now verify with the
  skins at the solid's density (`swell=False`), and compared like that the plates melt 2–6 % late on 2 mm elements, the
  3D model's resolution in depth (facts 133–137), now resolved by the fine band, 0.7 mm elements to 4 mm under the
  2 mm surface, which makes a melting step on the 100 mm sphere 45 % dearer (facts 138–140). To settle: whether flow's
  level-2 run uses the same setting, since with it the flowing layer's thickness is that of solid metal, not liquid.
- **The skins' 0.8 mm minimum** (option B′, Asha, 2026-10-16; facts 141–146): where the non-rigid layer is thinner than
  0.8 mm a skin's bottom is rigid metal, so the column flow must hold rigid cells still and count only the cells above
  T_rigid in the flowing depth; the finest cells stay in the top 0.4 mm (`SkinSettings.fine_zone`), and a column at its
  minimum holds 59 cells rather than 40, which the flow's per-cell work will feel. Flow makes each pass dearer,
  so the pass count is the first thing to settle: long parts during melting were rejected and needed many passes
  (fact 126), and parts of at most 0.0625 s while melting now cut a melting step to 122–226 s (fact 129). Absorption
  is limited only by the move allowance and the hot plate (level 2b) checks it; coarse surface meshes are not supported
  under strong melting, so flow's verification should use 2 mm surfaces (facts 130–131).
- Sub-plan 20's measured run time (fact 120 under the one-repeat rule; fact 125 under the amended default: 60–93 s a
  macro step solid, 123–1 041 s melting on the 100 mm sphere, almost all in the skins), which was to decide whether the
  sideways solve can be implicit over every cell or needs sub-cycling. Decided (Asha, 2026-10-17): implicit, the
  amendment above, from the level-3 flight's profile (facts 147 and 148); design (a) (Asha, 2026-10-18; fact 150).
- Whether the implicit thickness pressures are enough to keep grid-scale growth out (level 1's growth test), given
  that the capillary length (about 1.8 mm) is close to the column spacing.
