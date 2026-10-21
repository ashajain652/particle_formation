# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level so far: Tasks 1–5 and 6a (2026-10-18/21).**
> Tasks 6b, 6c, 7 and 8 are still the scoped outline at the end; each is written out with tested code in turn, the same
> way, and this document regenerated.

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

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -351,117 +351,141 @@
     return sk, col, pts, edge
 
 
-def test_sideways_heat_sums_to_zero_between_columns_of_any_depth(coarse_sphere_mesh):
+SPOT_SPACE_ERROR = 1.5            # K: the warm spot's error on 2 mm columns from the grid alone (1.29 K measured, 2026-10-18)
+
+
+def random_sphere_skins(coarse_sphere_mesh, seed=4):
+    """Skins on the 4 mm test sphere's corners, each 1-3 times its minimum deep, every cell at its own temperature."""
     from test_reentry_model_columns import sphere_columns
     mat = scheil()
     col, _ = sphere_columns(coarse_sphere_mesh)
     sk = skin.Skins(mat, col.area, skin.SkinSettings())
-    rng = np.random.default_rng(4)
+    rng = np.random.default_rng(seed)
     sk.fill(skin.minimum_mass(mat, sk.area, sk.settings) * rng.uniform(1.0, 3.0, col.n), rng.uniform(500.0, 900.0, col.n))
     sk.regrid()
     sk.T = np.where(sk.valid, sk.T + rng.normal(0.0, 20.0, sk.T.shape), sk.T)          # every cell different
     sk.H = sk.m * mat.enthalpy(sk.T)
-    Q = skin.sideways_heat(sk, col, 0.01)
+    return sk, col, rng
+
+
+def test_the_sideways_map_sums_to_zero_between_columns_of_any_depth(coarse_sphere_mesh):
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    Q = skin.Sideways(sk, col).heat(sk.T)
     assert abs(Q.sum()) < 1e-12 * np.abs(Q).sum()
     assert np.all(Q[~sk.valid] == 0.0)
 
 
-def test_sideways_heat_is_exact_for_a_field_linear_along_the_surface():
+def test_the_sideways_map_is_exact_for_a_field_linear_along_the_surface():
     mat = sharp()                                                    # constant k: a linear field carries no net flux
     sk, col, pts, edge = grid_skins(mat)
     T = 600.0 + 2e4 * pts[:, 0] - 1e4 * pts[:, 1]                       # 20 K/mm and -10 K/mm, uniform in depth
     sk.T = np.where(sk.valid, T[:, None], sk.T)
     sk.H = sk.m * mat.enthalpy(sk.T)
-    Q = skin.sideways_heat(sk, col, 1e-4)                            # one cycle: the edge columns do not move first
+    Q = skin.Sideways(sk, col).heat(sk.T)
     assert np.abs(Q[~edge]).max() < 1e-9 * np.abs(Q[edge]).max()
 
 
-def test_conduction_between_columns_reaches_the_known_steady_temperature():
-    """Columns on a flat plate, insulated top and bottom, the boundary columns held at a linear field: the steady state
-    is that linear field everywhere (Laplace's equation, with a constant conductivity)."""
-    mat = sharp()
-    sk, col, pts, edge = grid_skins(mat, n=6)
-    T_edge = 600.0 + 3e4 * pts[:, 0] + 1e4 * pts[:, 1]
-    # 4 ms keeps one sideways cycle per step even at the grid's corner columns (a sixth of the interior's area): with
-    # two, the held columns would drift between cycles and the interior settle 1.8 K off (found while planning)
-    for _ in range(200):                                             # 0.8 s: the slowest mode decays e^-22 over it
-        sk.T[edge] = np.where(sk.valid[edge], T_edge[edge][:, None], 300.0)
-        sk.H = sk.m * mat.enthalpy(sk.T)
-        skin._conduct(sk, 0.004, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.sideways_heat(sk, col, 0.004))
-    err = np.abs(sk.T[~edge][sk.valid[~edge]] - np.repeat(T_edge[~edge], sk.count[~edge]))
-    assert err.max() < 1e-9 * (T_edge.max() - T_edge.min())           # 3e-13 K measured: exact
-
-
-def test_a_sub_step_beyond_the_explicit_limit_is_sub_cycled_and_stays_bounded():
+def test_the_sideways_map_is_symmetric_and_negative_semi_definite(coarse_sphere_mesh):
+    """What the conjugate gradients of the implicit sub-step rely on: y . S x = x . S y, x . S x <= 0, and no heat at a
+    temperature uniform over the body."""
+    sk, col, rng = random_sphere_skins(coarse_sphere_mesh)
+    S = skin.Sideways(sk, col)
+    x, y = (np.where(sk.valid, rng.normal(size=sk.T.shape), 0.0) for _ in range(2))
+    yx, xy = float((y * S.heat(x)).sum()), float((x * S.heat(y)).sum())
+    assert abs(yx - xy) < 1e-10 * abs(yx)
+    assert float((x * S.heat(x)).sum()) < 0.0
+    uniform = np.where(sk.valid, 700.0, 0.0)
+    assert np.abs(S.heat(uniform)).max() < 1e-9 * np.abs(S.heat(x)).max() * 700.0
+
+
+def test_a_coupled_sub_step_without_sideways_differences_is_the_uncoupled_one():
+    """Identical columns heated alike carry nothing sideways: the coupled sub-step gives the uncoupled one's
+    temperatures, and its sideways books are zero."""
+    mat = scheil()
+    a, col, _, _ = grid_skins(mat, n=4)
+    b, _, _, _ = grid_skins(mat, n=4)
+    q = np.full(a.n, 2e6)
+    ra = skin._conduct(a, 0.01, q, a.area, 0.0, 0.15, np.full(a.n, 700.0), skin.Sideways(a, col))
+    rb = skin._conduct(b, 0.01, q, b.area, 0.0, 0.15, np.full(b.n, 700.0))
+    assert np.abs(a.T - b.T)[a.valid].max() < 1e-9
+    assert np.abs(ra[3]).max() < 1e-9 * np.abs(ra[0]).max()
+
+
+def test_a_sub_step_far_beyond_the_explicit_limit_stays_bounded_and_keeps_its_heat():
+    """A checkerboard 400 K / 800 K on 1 mm columns over one 0.2 s sub-step, about 30 times the explicit exchange's
+    stability limit there, insulated: implicit, it neither overshoots nor loses heat, and it smooths the board (to
+    622-629 K, measured 2026-10-18). The heat is kept to the Newton iteration's tolerance (3.5e-13 measured)."""
     mat = scheil()
     sk, col, pts, _ = grid_skins(mat, n=6)
-    hot = (np.rint(pts[:, 0] / 1e-3) + np.rint(pts[:, 1] / 1e-3)) % 2 == 0          # a checkerboard
+    hot = (np.rint(pts[:, 0] / 1e-3) + np.rint(pts[:, 1] / 1e-3)) % 2 == 0
     sk.T = np.where(sk.valid, np.where(hot, 800.0, 400.0)[:, None], sk.T)
     sk.H = sk.m * mat.enthalpy(sk.T)
     E0 = sk.energy()
-    Q = skin.sideways_heat(sk, col, 0.2)                             # 0.2 s: about 30 times the limit on 1 mm columns
-    sk.H = sk.H + Q * 0.2
-    T = mat.temperature_from_enthalpy(sk.H[sk.valid] / sk.m[sk.valid])
-    assert sk.energy() == pytest.approx(E0, rel=1e-13)
+    skin._conduct(sk, 0.2, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.Sideways(sk, col))
+    T = sk.T[sk.valid]
+    assert sk.energy() == pytest.approx(E0, rel=1e-11)
     assert T.min() >= 400.0 - 1e-6 and T.max() <= 800.0 + 1e-6
-
-
-def _sideways_as_built(sk, columns, dt):
-    """The sideways heat as sub-plan 20 built it (before 2026-10-17): the depth grid, the interpolation weights and the
-    temperatures made again in every sub-cycle, the last included."""
-    def interp(x, xp, fp):
-        n, q = xp.shape
-        off = (np.arange(n) * skin._SPAN)[:, None]
-        i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
-        i = np.clip(i - (np.arange(n) * q)[:, None], 0, q - 2)
-        r = np.arange(n)[:, None]
-        x0, x1, f0, f1 = xp[r, i], xp[r, i + 1], fp[r, i], fp[r, i + 1]
-        t = np.clip(np.divide(x - x0, x1 - x0, out=np.zeros_like(x, dtype=float), where=x1 > x0), 0.0, 1.0)
-        return f0 + t * (f1 - f0)
-
-    def once(T, dz):
-        n = T.shape[0]
-        zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)
-        H = zb[:, -1]
-        zc = np.minimum(skin.depth_grid(sk.settings, float(H.max()))[None, :], H[:, None])
-        Cg = interp(zc, zb, np.concatenate([np.zeros((n, 1)), np.cumsum(T * dz, axis=1)], axis=1))
-        length = np.diff(zc, axis=1)
-        Tbar = np.divide(np.diff(Cg, axis=1), length, out=np.zeros_like(length), where=length > 0.0)
-        i, j = columns.edges[:, 0], columns.edges[:, 1]
-        k = sk.props.k(Tbar)
-        flux = (0.5 * (k[i] + k[j]) * columns.conduct_weight[:, None] * np.minimum(length[i], length[j])
-                * (Tbar[j] - Tbar[i]))
-        Qg = np.zeros_like(length)
-        np.add.at(Qg, i, flux)
-        np.add.at(Qg, j, -flux)
-        return np.diff(interp(zb, zc, np.concatenate([np.zeros((n, 1)), np.cumsum(Qg, axis=1)], axis=1)), axis=1)
-
-    mat, valid, dz = sk.material, sk.valid, sk.dz()
-    wsum = np.bincount(columns.edges.ravel(), weights=np.repeat(columns.conduct_weight, 2), minlength=sk.n)
-    rate = float((mat.k_table.max() * wsum / (mat.rho * mat.cp_table.min() * sk.area)).max())
-    cycles = max(1, int(np.ceil(dt * rate / skin.SIDEWAYS_CFL)))
-    T, H, total = sk.T.copy(), sk.H.copy(), np.zeros_like(sk.T)
-    for _ in range(cycles):
-        Q = np.where(valid, once(T, dz), 0.0)
-        total += Q / cycles
-        if cycles > 1:
-            H = H + Q * dt / cycles
-            T = np.where(valid, mat.temperature_from_enthalpy(np.divide(H, sk.m, out=np.zeros_like(H), where=sk.m > 0.0)), T)
-    return total, cycles
-
-
-def test_the_sub_cycles_share_their_geometry_and_give_the_heat_as_built_to_the_last_bit():
-    """Asha, 2026-10-17: the sub-cycles of a sideways step share its depth grid and interpolation weights, and the
-    temperatures after the last are not computed -- a third of the sideways heat's time on the 100 mm sphere. Columns of
-    different depths, every cell at its own temperature, over several sub-cycles: the heat is the build's exactly."""
-    mat = scheil()
-    sk, col, pts, _ = grid_skins(mat, n=6)
-    rng = np.random.default_rng(7)
-    sk.fill(skin.minimum_mass(mat, sk.area, sk.settings) * rng.uniform(1.0, 2.5, col.n), rng.uniform(450.0, 850.0, col.n))
-    sk.regrid()
-    sk.T = np.where(sk.valid, sk.T + rng.normal(0.0, 15.0, sk.T.shape), sk.T)
-    sk.H = sk.m * mat.enthalpy(sk.T)
-    for dt in (1e-4, 0.05):
-        reference, cycles = _sideways_as_built(sk, col, dt)
-        np.testing.assert_array_equal(skin.sideways_heat(sk, col, dt), reference)
-    assert cycles >= 3                                               # the second dt is sub-cycled
+    assert T.max() - T.min() < 20.0
+
+
+def test_a_warm_spot_spreads_as_the_heat_equation_says_and_converges_with_the_sub_step():
+    """A Gaussian warm spot (100 K, 4 mm) in a flat array of 2 mm columns, uniform in depth and insulated, spreads over
+    0.2 s as the two-dimensional heat equation's exact solution, to the grid's own spatial error; against the same
+    grid at a 2.5 ms sub-step the time error falls in proportion to the sub-step (backward Euler, first order: the
+    errors at 40, 20 and 10 ms stand as 37.5 : 17.5 : 7.5, their distances from the reference's 2.5 ms; measured 2.79,
+    1.32 and 0.57 K, 2026-10-18)."""
+    mat = sharp()
+    alpha = 100.0 / (2800.0 * 1000.0)
+    t0, T_end = (4e-3) ** 2 / (4.0 * alpha), 0.2
+
+    def run(h):
+        sk, col, pts, _ = grid_skins(mat, n=21, h=2e-3, T=500.0)
+        c = pts[:, :2] - pts[:, :2].mean(axis=0)
+        r2 = (c * c).sum(axis=1)
+        sk.T = np.where(sk.valid, (500.0 + 100.0 * np.exp(-r2 / (4.0 * alpha * t0)))[:, None], sk.T)
+        sk.H = sk.m * mat.enthalpy(sk.T)
+        for _ in range(int(round(T_end / h))):
+            skin._conduct(sk, h, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.Sideways(sk, col))
+        exact = 500.0 + 100.0 * t0 / (t0 + T_end) * np.exp(-r2 / (4.0 * alpha * (t0 + T_end)))
+        return sk.T[:, 0], exact, r2 < (12e-3) ** 2
+
+    ref, exact, inner = run(0.0025)
+    err = {h: np.abs(run(h)[0] - ref)[inner].max() for h in (0.04, 0.02, 0.01)}
+    assert err[0.04] / err[0.02] == pytest.approx(37.5 / 17.5, rel=0.1)
+    assert err[0.02] / err[0.01] == pytest.approx(17.5 / 7.5, rel=0.1)
+    assert np.abs(ref - exact)[inner].max() < SPOT_SPACE_ERROR
+
+
+def test_a_sideways_map_reuses_the_last_geometry_only_while_the_cells_keep_their_thicknesses(coarse_sphere_mesh):
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    first = skin.Sideways(sk, col)
+    assert skin.Sideways(sk, col, first).geometry is first.geometry
+    sk.add_bottom(np.full(sk.n, 1e-7), np.full(sk.n, float(sk.material.enthalpy(700.0))), np.zeros(sk.n))
+    sk.regrid()
+    again, fresh = skin.Sideways(sk, col, first), skin.Sideways(sk, col)
+    assert again.geometry is not first.geometry
+    np.testing.assert_array_equal(again.heat(sk.T), fresh.heat(sk.T))
+
+
+def test_the_inexact_coupled_solves_reach_the_tight_solve_s_answer(coarse_sphere_mesh, monkeypatch):
+    """Columns of different depths melting under 3 MW/m2 on the front, every cell at its own temperature: each Newton
+    step solved to KRYLOV_TOL gives the temperatures that steps solved to 1e-10 give, to the Newton tolerance."""
+    out = []
+    for tol in (skin.KRYLOV_TOL, 1e-10):
+        monkeypatch.setattr(skin, "KRYLOV_TOL", tol)
+        sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+        q = 3e6 * np.clip(col.normal[:, 0], 0.0, None)
+        books = skin._conduct(sk, 0.01, q, sk.area, 0.0, 0.15, np.full(sk.n, 850.0), skin.Sideways(sk, col))
+        assert abs(books[3].sum()) < 1e-12 * np.abs(books[3]).sum()            # the sideways books sum to zero
+        out.append(sk.T[sk.valid].copy())
+    assert np.abs(out[0] - out[1]).max() < 1e-5
+
+
+def test_a_coupled_solve_that_cannot_converge_raises(coarse_sphere_mesh, monkeypatch):
+    """A Newton step whose conjugate gradients run out of iterations stops the run (exit code 1 from the command
+    line), never a sub-step taken on an unconverged step."""
+    monkeypatch.setattr(skin, "KRYLOV_MAX", 0)
+    monkeypatch.setattr(skin, "KRYLOV_TOL", 1e-14)
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    with pytest.raises(RuntimeError, match="coupled sideways solve did not converge"):
+        skin._conduct(sk, 0.01, np.zeros(sk.n), sk.area, 0.0, 0.0, np.full(sk.n, 700.0), skin.Sideways(sk, col))
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q -k "sideways or
  coupled or warm_spot or explicit_limit"`. Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'Sideways'`.

- [ ] **Step 3: The edges' incidence matrix** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -10,12 +10,14 @@
 with a and b the angles opposite the edge in its two triangles: the flux of a piecewise-linear field through the two
 dual segments across the edge is exactly w (T_j - T_i) per unit conductivity and depth, so the discrete operator is
 exact for a field linear along the surface. A weight is clipped at zero where the opposite angles are obtuse enough
-to make it negative: a negative conductance would break the explicit step's positivity, and it is rare (3 of 6 756
-edges on the 4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
+to make it negative: a negative conductance would break the explicit step's positivity (and, since sub-plan 21 makes the
+sideways conduction implicit, the definiteness its conjugate gradients need), and it is rare (3 of 6 756 edges on the
+4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
 The surface gradient and curvature operators the flow needs are added in sub-plan 21, the first stage that uses them."""
 from dataclasses import dataclass
 
 import numpy as np
+import scipy.sparse as sp
 
 from .movemesh import corner_frame
 
@@ -82,6 +84,17 @@
         np.add.at(out, j, -d)
         return out
 
+    @property
+    def incidence(self):
+        """B (n columns x E edges), +1 at each edge's first column and -1 at its second, so B @ f gives every column the
+        sum over its edges of f, with the first column's sign. Built once: the connectivity is fixed."""
+        if getattr(self, "_incidence", None) is None:
+            i, j = self.edges[:, 0], self.edges[:, 1]
+            e = np.arange(len(i))
+            self._incidence = sp.csr_matrix((np.r_[np.ones(len(i)), -np.ones(len(i))], (np.r_[i, j], np.r_[e, e])),
+                                            shape=(self.n, len(i)))
+        return self._incidence
+
 
 @dataclass
 class CornerSurface:
```

- [ ] **Step 4: The coupled sub-step and the sideways map** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -47,7 +47,10 @@
 NEWTON_MAX = 40
 LINE_SEARCH_MAX = 12              # halvings of a Newton step before it is taken anyway
 DZ_FLOOR = 1.0e-12                # m: the conduction length of a vanishing cell
-SIDEWAYS_CFL = 0.9                # explicit sideways conduction: sub-cycles keep dt k sum(w) / (rho c_p A) below this
+KRYLOV_TOL = 1.0e-3               # each Newton step's coupled solve (sideways conduction in the implicit sub-step) stops
+                                  # at this residual relative to its right-hand side: the Newton iteration corrects the rest
+KRYLOV_MAX = 200                  # conjugate-gradient iterations of one coupled solve at most
+KRYLOV = {"solves": 0, "iterations": 0}   # the coupled solves made and their iterations, for measuring the cost
 TABLE_DT = 0.01                   # K, spacing of the skins' property tables (PropertyTable)
 TABLE_RANGE = (150.0, 3000.0)     # K, their extent; beyond it the end intervals are extended linearly
 PROFILE_POINTS = 8                # Gauss-Legendre points of an enthalpy averaged over a depth of the 3D profile
@@ -459,14 +462,53 @@
     drawn_H: np.ndarray      # the enthalpy it brought
 
 
-def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, Q_side=None):
+def _coupled_solve(G, diag, rhs, valid, sideways):
+    """One Newton step of the coupled sub-step: J x = rhs with J the columns' tridiagonal Jacobians less the sideways
+    map (`Sideways.heat`, linear), symmetric positive definite, by conjugate gradients preconditioned by the columns' own
+    tridiagonal solves (`thomas`) and started from them, to KRYLOV_TOL of the right-hand side. Returns x and the sideways
+    heat at x, built from the products the iteration made (the map is linear), which the line search reuses."""
+    def depth(x):
+        y = diag * x
+        y[:, :-1] -= G * x[:, 1:]
+        y[:, 1:] -= G * x[:, :-1]
+        return y
+
+    def solve_columns(r):
+        return np.where(valid, thomas(-G, diag, -G, r), 0.0)
+
+    x = solve_columns(rhs)
+    Hx = sideways.heat(x)
+    r = np.where(valid, rhs - depth(x) + Hx, 0.0)
+    target = KRYLOV_TOL * math.sqrt(float((rhs * rhs).sum()))
+    z = solve_columns(r)
+    p = z.copy()
+    rz = float((r * z).sum())
+    KRYLOV["solves"] += 1
+    for it in range(KRYLOV_MAX + 1):
+        if math.sqrt(float((r * r).sum())) <= target:
+            KRYLOV["iterations"] += it
+            return x, Hx
+        Hp = sideways.heat(p)
+        Jp = np.where(valid, depth(p) - Hp, 0.0)
+        a = rz / float((p * Jp).sum())
+        x, Hx, r = x + a * p, Hx + a * Hp, r - a * Jp
+        z = solve_columns(r)
+        rz, rz_old = float((r * z).sum()), rz
+        p = z + (rz / rz_old) * p
+    raise RuntimeError("a skin sub-step's coupled sideways solve did not converge in {} iterations".format(KRYLOV_MAX))
+
+
+def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, sideways=None):
     """One implicit (backward-Euler) sub-step for every column: the convective flux q_top [W/m2] over the top area
     A_top in at the top cell, radiation from it over the same area, conduction through the cells (harmonic
-    conductances), the explicit sideways heat Q_side [W per cell], and the base held at T_b [K] through the bottom
-    half-cell (T_b None: an adiabatic base, for the tests). Newton in the temperature with a backtracking line search
-    per column (stage 1's solver, sub-plan 18). Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the
-    skins track S, updates it to dT/dT_b and returns the base energy's derivative with respect to T_b as a fifth item.
-    E_base is what the books leave, so they close exactly whatever the iteration error."""
+    conductances), the sideways conduction `sideways` (a `Sideways` map made at the sub-step's start, implicit with the
+    rest: sub-plan 21, design (a)), and the base held at T_b [K] through the bottom half-cell (T_b None: an adiabatic
+    base, for the tests). Newton in the temperature with a backtracking line search per column (stage 1's solver,
+    sub-plan 18); with the sideways map the columns are coupled and each Newton step is solved by `_coupled_solve`.
+    Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the skins track S, updates it to dT/dT_b and returns
+    the base energy's derivative with respect to T_b as a fifth item -- through the depth only: the base law each column
+    gives the 3D model is its own, and what the neighbours' share misses is the interface mismatch's to book. E_base is
+    what the books leave, so they close exactly whatever the iteration error."""
     mat = sk.props
     n, K = sk.m.shape
     valid = sk.valid
@@ -482,14 +524,17 @@
     q = np.asarray(q_top, dtype=float) * np.asarray(A_top, dtype=float)
     base = T_b is not None
     T_b = np.asarray(T_b if base else 0.0, dtype=float) * np.ones(n)
-    Qs = np.zeros((n, K)) if Q_side is None else np.where(valid, Q_side, 0.0)
-
-    def system(T):
+    side = (lambda T: np.zeros((n, K))) if sideways is None else sideways.heat
+
+    def system(T, Q=None):
+        """The residual, the conductances and the Jacobian's diagonal at T, and the sideways heat Q at T (given when the
+        caller has it already)."""
+        Q = side(T) if Q is None else Q
         h, cp, k = mat.conduction(T)
         G = np.where(link, A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:]), 0.0)
         Gb = A * k[rows, last] / (0.5 * dz[rows, last]) if base else np.zeros(n)
         flux = G * (T[:, :-1] - T[:, 1:])
-        net = Qs.copy()
+        net = Q.copy()
         net[:, :-1] -= flux
         net[:, 1:] += flux
         net[:, 0] += q - eps_s * (T[:, 0] ** 4 - T_amb ** 4)
@@ -500,27 +545,31 @@
         diag[:, 1:] += G
         diag[rows, last] += Gb
         R = np.where(valid, (m * h - E_old) / dt - net, 0.0)
-        return R, G, Gb, np.where(valid, diag, 1.0)
+        return R, G, Gb, np.where(valid, diag, 1.0), Q
 
     # Newton in the temperature with a backtracking line search per column on the residual's norm: in 10-micrometre
     # cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step, so a cell must move
     # with its neighbours, and the 3D solver's per-node enthalpy-consistent update stalled here (sub-plan 18)
     T = T_old.copy()
-    R, G, Gb, diag = system(T)
+    R, G, Gb, diag, Q = system(T)
     moved = np.zeros(n)
     for _ in range(NEWTON_MAX):
-        dT = thomas(-G, diag, -G, -R)
+        if sideways is None:
+            dT, Q_dT = thomas(-G, diag, -G, -R), None
+        else:
+            dT, Q_dT = _coupled_solve(G, diag, -R, valid, sideways)
         r0 = np.sqrt((R * R).sum(axis=1))
         lam = np.ones(n)
         for _ in range(LINE_SEARCH_MAX):
             T_try = T + lam[:, None] * dT
-            R_try, G_try, Gb_try, diag_try = system(T_try)
+            full = Q_dT is not None and bool(np.all(lam == 1.0))         # the map is linear: Q(T + dT) = Q(T) + Q(dT)
+            R_try, G_try, Gb_try, diag_try, Q_try = system(T_try, Q + Q_dT if full else None)
             moved = np.abs(T_try - T).max(axis=1)
             worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
             if not worse.any():
                 break
             lam = np.where(worse, 0.5 * lam, lam)
-        T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try
+        T, R, G, Gb, diag, Q = T_try, R_try, G_try, Gb_try, diag_try, Q_try
         if moved.max() <= NEWTON_TOL:
             break
     else:
@@ -530,7 +579,7 @@
     E_new = m * mat.enthalpy(T)
     E_top = q * dt
     E_rad = eps_s * (T[:, 0] ** 4 - T_amb ** 4) * dt
-    E_side = Qs.sum(axis=1) * dt
+    E_side = Q.sum(axis=1) * dt
     E_base = E_top - E_rad + E_side - (E_new.sum(axis=1) - E_old.sum(axis=1))
     out = [E_top, E_rad, E_base, E_side]
     if sk.S is not None:
@@ -544,8 +593,10 @@
 
 
 def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None):
-    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the sideways heat `side(sk, h)` [W per
-    cell] if given (Task 5), the implicit conduction, then -- with `remove` and a melting material -- the contiguous top
+    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the implicit conduction, through the depth
+    and, with `side`, sideways (`side(sk, previous)` makes the sub-step's `Sideways` map from the last one, Task 5),
+    then -- with `remove` and a
+    melting material -- the contiguous top
     cells at or above T_feed leave, a column below its minimum mass m_min [kg] draws back to it, and the cells are split
     and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but never above T_rigid -- a 3D base
     node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it reached 930 K in one 0.5 s step, found
@@ -573,9 +624,10 @@
     s_b = (T_b < mat.T_rigid).astype(float)           # drawn metal follows the base only where the draw is not capped
     drawn_depth = np.zeros(sk.n)                      # m, how far down the 3D profile this pass has drawn
     melting = remove and mat.melts
+    sideways = None
     for _ in range(N):
-        Q_side = side(sk, h) if side is not None else None
-        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, Q_side)
+        sideways = side(sk, sideways) if side is not None else None
+        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, sideways)
         book.E_top += res[0]
         book.E_rad += res[1]
         book.E_base += res[2]
@@ -613,8 +665,8 @@
 
 def _interp_weights(x, xp):
     """Row-wise linear interpolation at x (n, p) of a function given at non-decreasing xp (n, q), held at the ends, as
-    (rows, index, fraction) for `_interp_apply`: they depend on the abscissae only, so the sub-cycles of a sideways step,
-    whose cells keep their thicknesses, share them."""
+    (rows, index, fraction) for `_interp_apply`: they depend on the abscissae only, so a sub-step's sideways map makes
+    them once and applies them at every iteration of its solve."""
     n, q = xp.shape
     off = (np.arange(n) * _SPAN)[:, None]
     i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
@@ -642,9 +694,9 @@
 
 
 def _sideways_geometry(sk, columns, dz):
-    """What every sub-cycle of a sideways step shares, from the cells' thicknesses dz alone: the cells' boundaries from
-    the top, each column's share of each bin of the depth grid, the depth each edge's two columns share in each bin,
-    and the interpolation weights from cells to bins and back."""
+    """What a sub-step's sideways map is built on, from the cells' thicknesses dz alone: the cells' boundaries from the
+    top, each column's share of each bin of the depth grid, the depth each edge's two columns share in each bin, and the
+    interpolation weights from cells to bins and back."""
     n = dz.shape[0]
     zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)          # cell boundaries, from the top
     H = zb[:, -1]
@@ -656,54 +708,43 @@
             "to_bins": _interp_weights(zc, zb), "to_cells": _interp_weights(zb, zc)}
 
 
-def _sideways_once(sk, columns, T, dz, geometry=None):
-    """Sideways heat [W per cell] at temperatures T: each column's cells averaged onto the depth grid, every surface
-    edge exchanging k w l (T_j - T_i) in every bin both columns reach (k the mean of the two bins' conductivities, w the
-    edge's cotangent weight, l the depth the two columns share in the bin), and each column's bin heat given back to its
-    cells in proportion to their overlap with the bin. Antisymmetric edge by edge, so the total is zero to round-off;
-    exact for a field linear along the surface where neighbouring columns have the same cells. `geometry`: the step's
-    `_sideways_geometry`, made here if not given."""
-    n = T.shape[0]
-    g = _sideways_geometry(sk, columns, dz) if geometry is None else geometry
-    C = np.concatenate([np.zeros((n, 1)), np.cumsum(T * dz, axis=1)], axis=1)      # integral of T from the top
-    Cg = _interp_apply(g["to_bins"], C)
-    length = g["length"]
-    Tbar = np.divide(np.diff(Cg, axis=1), length, out=np.zeros_like(length), where=length > 0.0)
-    i, j = columns.edges[:, 0], columns.edges[:, 1]
-    k = sk.props.k(Tbar)
-    flux = 0.5 * (k[i] + k[j]) * columns.conduct_weight[:, None] * g["shared"] * (Tbar[j] - Tbar[i])
-    Qg = np.zeros_like(length)
-    np.add.at(Qg, i, flux)
-    np.add.at(Qg, j, -flux)
-    Qcum = np.concatenate([np.zeros((n, 1)), np.cumsum(Qg, axis=1)], axis=1)
-    return np.diff(_interp_apply(g["to_cells"], Qcum), axis=1)
-
-
-def sideways_heat(sk, columns, dt):
-    """The sideways conduction of one sub-step dt [s] as an explicit heat source [W per cell], the time average over
-    as many sub-cycles as keep dt k_max sum(w) / (rho c_p,min A) below SIDEWAYS_CFL (spec 6.3); the cells' temperatures
-    move between sub-cycles, the source being what `_conduct` then adds. The smallest columns set the count: on the
-    100 mm sphere's 2 mm surface, dual areas of 1.26-4.59 mm2 (median 3.31) give a cycle of 3.8 ms, three to a 10 ms
-    sub-step, where the median column alone would need one (measured 2026-10-17). The sub-cycles share the step's
-    geometry, and the temperatures after the last are not computed (Asha, 2026-10-17: remaking the one in every
-    sub-cycle and computing the other took a third of the sideways heat's time; the result is the same to the last
-    bit)."""
-    mat = sk.material
-    valid = sk.valid
-    dz = sk.dz()
-    wsum = np.bincount(columns.edges.ravel(), weights=np.repeat(columns.conduct_weight, 2), minlength=sk.n)
-    rate = float((mat.k_table.max() * wsum / (mat.rho * mat.cp_table.min() * sk.area)).max())
-    cycles = max(1, int(math.ceil(dt * rate / SIDEWAYS_CFL)))
-    geometry = _sideways_geometry(sk, columns, dz)
-    T, H = sk.T.copy(), sk.H.copy()
-    total = np.zeros_like(T)
-    for c in range(cycles):
-        Q = np.where(valid, _sideways_once(sk, columns, T, dz, geometry), 0.0)
-        total += Q / cycles
-        if c < cycles - 1:
-            H = H + Q * dt / cycles
-            T = np.where(valid, mat.temperature_from_enthalpy(np.divide(H, sk.m, out=np.zeros_like(H), where=sk.m > 0.0)), T)
-    return total
+class Sideways:
+    """The sideways conduction of one sub-step as a linear map on the cells' temperatures (spec 6.3; implicit since
+    sub-plan 21, design (a), Asha 2026-10-18). Each column's cells are averaged onto the common depth grid, every surface
+    edge exchanges c (x_j - x_i) in every bin both columns reach, with c = k w l (k the mean of the two bins'
+    conductivities at the sub-step's start, w the edge's cotangent weight, l the depth the two columns share in the bin),
+    and each column's bin heat goes back to its cells in proportion to their overlap with the bin. `heat(x)` [W per cell]
+    is antisymmetric edge by edge, so it sums to zero to round-off, and exact for a field linear along the surface where
+    neighbouring columns have the same cells. As a map it is symmetric and negative semi-definite -- the bins' averaging
+    and the heat's return to the cells are each other's transposes -- so `_conduct` takes it into its implicit solve
+    (`_coupled_solve`), and the sub-step has no sideways stability limit: the explicit exchange it replaces needed three
+    sub-cycles a 10 ms sub-step on the 100 mm sphere's smallest columns (facts 147-148)."""
+
+    def __init__(self, sk, columns, previous=None):
+        """`previous`: the last sub-step's map, whose geometry is reused while the cells keep their thicknesses (in a
+        body that does not melt, all of its sub-steps)."""
+        self.valid, self.dz = sk.valid, sk.dz()
+        same = previous is not None and previous.dz.shape == self.dz.shape and np.array_equal(previous.dz, self.dz)
+        self.geometry = previous.geometry if same else _sideways_geometry(sk, columns, self.dz)
+        self.i, self.j = columns.edges[:, 0], columns.edges[:, 1]
+        k = sk.props.k(self.bins(sk.T))
+        self.c = 0.5 * (k[self.i] + k[self.j]) * columns.conduct_weight[:, None] * self.geometry["shared"]
+        self.B = columns.incidence
+
+    def bins(self, x):
+        """The mean of x (per cell) over each bin of the depth grid, (n, bins); zero in bins beyond a column."""
+        n = x.shape[0]
+        g = self.geometry
+        C = np.concatenate([np.zeros((n, 1)), np.cumsum(x * self.dz, axis=1)], axis=1)      # integral from the top
+        return np.divide(np.diff(_interp_apply(g["to_bins"], C), axis=1), g["length"], out=np.zeros_like(g["length"]),
+                         where=g["length"] > 0.0)
+
+    def heat(self, x):
+        """The heat [W per cell] the neighbouring columns conduct into each cell at cell temperatures x (linear in x)."""
+        xb = self.bins(x)
+        Qg = self.B @ (self.c * (xb[self.j] - xb[self.i]))
+        Qcum = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(Qg, axis=1)], axis=1)
+        return np.where(self.valid, np.diff(_interp_apply(self.geometry["to_cells"], Qcum), axis=1), 0.0)
 
 
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
```

- [ ] **Step 5: The skin body's hook** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -208,7 +208,7 @@
         A_top = self.surface.areas
         eps = np.full(sk.n, self.emissivity)
         m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
-        side = lambda s, h: skin_mod.sideways_heat(s, self.columns, h)
+        side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
         nodes = mm.corners
         T_b0 = self.solver.temperature()[nodes]
         g0 = self._normal_gradient()                         # the profile's slope below every base, for this part
```

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

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -489,3 +489,58 @@
     sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
     with pytest.raises(RuntimeError, match="coupled sideways solve did not converge"):
         skin._conduct(sk, 0.01, np.zeros(sk.n), sk.area, 0.0, 0.0, np.full(sk.n, 700.0), skin.Sideways(sk, col))
+
+
+# -- the cells' velocity (sub-plan 21 Task 2) ---------------------------------------------------------------------------
+def moving_skins(seed=5):
+    """Five columns filled to two to four times their minimum in fine cells (so the first regrid merges the deep ones),
+    every cell at its own velocity in metres a second."""
+    mat = scheil()
+    s = skin.SkinSettings()
+    sk = skin.Skins(mat, np.full(5, 4e-6), s)
+    rng = np.random.default_rng(seed)
+    sk.fill(skin.minimum_mass(mat, sk.area, s) * rng.uniform(2.0, 4.0, 5), rng.uniform(600.0, 850.0, 5))
+    sk.U = np.where(sk.valid[..., None], rng.normal(0.0, 0.3, sk.m.shape + (3,)), 0.0)
+    return sk, rng
+
+
+def test_a_filled_skin_is_at_rest_and_its_velocity_is_part_of_its_state():
+    mat = scheil()
+    sk = filled(mat)
+    assert np.all(sk.U == 0.0) and sk.U.shape == sk.m.shape + (3,)
+    moving, _ = moving_skins()
+    saved = moving.state()
+    P = moving.column_momentum()
+    moving.U = np.zeros_like(moving.U)
+    moving.set_state(saved)
+    np.testing.assert_array_equal(moving.column_momentum(), P)
+
+
+def test_split_and_merge_keep_every_column_s_momentum():
+    sk, _ = moving_skins()
+    P0, count0 = sk.column_momentum(), sk.count.copy()
+    sk.regrid()                                                       # merges the deep fine cells
+    assert np.all(sk.count < count0)
+    np.testing.assert_allclose(sk.column_momentum(), P0, rtol=1e-12, atol=1e-18)
+    sk.add_bottom(np.full(sk.n, 0.5 * sk.column_mass().min()), np.full(sk.n, float(sk.material.enthalpy(700.0))))
+    P1, count1 = sk.column_momentum(), sk.count.copy()
+    sk.regrid()                                                       # splits the bottom cell the draw thickened
+    assert np.all(sk.count > count1)
+    np.testing.assert_allclose(sk.column_momentum(), P1, rtol=1e-12, atol=1e-18)
+    np.testing.assert_allclose(P1, P0, rtol=1e-12, atol=1e-18)        # the drawn metal arrived at rest
+
+
+def test_what_leaves_at_the_top_or_the_bottom_takes_its_momentum_and_what_arrives_brings_its_own():
+    sk, rng = moving_skins()
+    sk.regrid()
+    P = sk.column_momentum()
+    mass, _ = sk.take_top(np.array([0, 1, 2, 3, 1]))
+    np.testing.assert_allclose(sk.column_momentum() + sk.momentum_out, P, rtol=1e-12, atol=1e-18)
+    assert np.all((np.linalg.norm(sk.momentum_out, axis=1) > 0.0) == (mass > 0.0))
+    P = sk.column_momentum()
+    sk.take_bottom(0.3 * sk.column_mass())                            # whole cells and a part of one
+    np.testing.assert_allclose(sk.column_momentum() + sk.momentum_out, P, rtol=1e-12, atol=1e-18)
+    P = sk.column_momentum()
+    dm, u = np.full(sk.n, 1e-7), rng.normal(0.0, 0.2, (sk.n, 3))
+    sk.add_bottom(dm, np.full(sk.n, float(sk.material.enthalpy(700.0))), u=u)
+    np.testing.assert_allclose(sk.column_momentum(), P + dm[:, None] * u, rtol=1e-12, atol=1e-18)
```

- [ ] **Step 2: Append the surface operators' tests** to `tests/test_reentry_model_columns.py`:

```diff
--- a/tests/test_reentry_model_columns.py
+++ b/tests/test_reentry_model_columns.py
@@ -87,3 +87,46 @@
     v = np.array([-1.0, 0.0, 0.0])
     assert top.angles_to(v)[top.patch_toward(v)] < np.radians(5.0)
     assert top.projected_area(v) == pytest.approx(np.pi * (0.045 + 1e-3) ** 2, rel=0.01)
+
+
+# -- the surface operators (sub-plan 21 Task 2) -------------------------------------------------------------------------
+def obtuse_corners(col, points=None):
+    p = col.points if points is None else points
+    _, _, cot = columns.triangle_geometry(p, col.faces)
+    out = np.zeros(col.n, dtype=bool)
+    for v in range(3):
+        np.logical_or.at(out, col.faces[:, v], (cot < 0.0).any(axis=1))
+    return out
+
+
+def test_the_gradient_is_exact_for_a_field_linear_on_a_flat_surface():
+    pts, faces, inner = flat_grid()
+    col = columns.Columns(pts, faces)
+    g = col.gradient(2.0 + 300.0 * pts[:, 0] - 700.0 * pts[:, 1])
+    assert np.abs(g - np.array([300.0, -700.0, 0.0])).max() < 1e-9          # every corner, the edges' included
+
+
+def test_on_a_sphere_the_gradient_of_a_linear_field_is_its_tangential_part(coarse_sphere_mesh):
+    """On the 4 mm test sphere within 0.09 % of a linear field's tangential part on average and 0.8 % at worst, at the
+    irregular corner gmsh left (0.04 % at worst on the 2 mm sphere; measured 2026-10-18)."""
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    a = np.array([3.0, -1.0, 2.0])
+    exact = a - (col.normal @ a)[:, None] * col.normal
+    err = np.linalg.norm(col.gradient(col.points @ a) - exact, axis=1) / np.linalg.norm(a)
+    assert err.mean() < 0.002 and err.max() < 0.01
+
+
+def test_the_curvature_of_a_sphere_is_two_over_its_radius_at_every_corner(coarse_sphere_mesh):
+    """Every corner within 2e-4 of 2/R, those touching an obtuse triangle included (33 of the 2 254, which the mixed
+    Voronoi area alone left off by up to 176 %; 8.4e-5 measured at worst, 2026-10-18), and the same on a sphere 2 %
+    larger whose corners are given."""
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    assert obtuse_corners(col).sum() > 0                                     # the fallback is exercised
+    R = np.linalg.norm(col.points, axis=1)
+    assert np.abs(col.curvature() * R / 2.0 - 1.0).max() < 2e-4
+    assert np.abs(col.curvature(1.02 * col.points) * 1.02 * R / 2.0 - 1.0).max() < 2e-4
+
+
+def test_the_curvature_of_a_plane_is_zero():
+    pts, faces, inner = flat_grid()
+    assert np.abs(columns.Columns(pts, faces).curvature()[inner]).max() < 1e-9
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_columns.py -q`. Expected: `AttributeError: 'Skins' object has no attribute 'U'` and
  `AttributeError: module 'reentry_model.columns' has no attribute 'triangle_geometry'`.

- [ ] **Step 4: The velocity in the cells** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -143,7 +143,10 @@
 class Skins:
     """The skins of every column, as padded arrays (module docstring). `area` [m2] is each column's area, the dual
     area of its base corner; the body updates it when the corners move (the cells keep their mass). `S`, when not None,
-    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b)."""
+    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b). `U` (n, K, 3) is every
+    cell's velocity [m/s] (sub-plan 21 Task 2): part of the cells' state, carried by every operation with its momentum
+    exact -- split halves keep it, merged cells take the mass-weighted mean, metal drawn from the 3D model arrives at
+    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top` or `take_bottom` took away with the metal."""
 
     def __init__(self, material, area, settings=None):
         self.material = material
@@ -153,6 +156,8 @@
         self.m, self.H, self.T = np.zeros((n, 1)), np.zeros((n, 1)), np.full((n, 1), 300.0)
         self.count = np.zeros(n, dtype=np.int64)
         self.S = None
+        self.U = np.zeros((n, 1, 3))
+        self.momentum_out = np.zeros((n, 3))
         self.rho_l = material.liquid.rho if material.liquid is not None and self.settings.swell else material.rho
         self.props = PropertyTable(material)
 
@@ -174,11 +179,18 @@
     def column_mass(self):
         return self.m.sum(axis=1)
 
+    def column_momentum(self):
+        """Each column's momentum [kg m/s], (n, 3)."""
+        return (self.m[..., None] * self.U).sum(axis=1)
+
+    def momentum(self):
+        return self.column_momentum().sum(axis=0)
+
     def state(self):
-        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count))
+        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count, self.U))
 
     def set_state(self, state):
-        self.area, self.m, self.H, self.T, self.count = (a.copy() for a in state)
+        self.area, self.m, self.H, self.T, self.count, self.U = (a.copy() for a in state)
         self.S = None
 
     # -- geometry of the cells --------------------------------------------------------------------------------------
@@ -202,10 +214,10 @@
         return s.cell + (s.growth - 1.0) * np.maximum(z - s.fine_zone, 0.0)
 
     # -- installing cells -------------------------------------------------------------------------------------------
-    def _install(self, m, H, count, S=None, T=None):
+    def _install(self, m, H, count, S=None, T=None, U=None):
         """Install padded cell arrays. Temperatures follow from the enthalpies, except where `T` is given for cells
         whose specific enthalpy the operation did not change (a shifted, split or partly taken cell); nan in T marks
-        the cells to recompute."""
+        the cells to recompute. `U` is the cells' velocity (at rest if not given)."""
         valid = np.arange(m.shape[1])[None, :] < count[:, None]
         self.m, self.H, self.count = np.where(valid, m, 0.0), np.where(valid, H, 0.0), np.asarray(count, dtype=np.int64)
         redo = valid if T is None else valid & np.isnan(T)
@@ -214,6 +226,7 @@
             self.T[redo] = self.material.temperature_from_enthalpy(self.H[redo] / self.m[redo])
         if self.S is not None:
             self.S = np.where(valid, S, 0.0)
+        self.U = np.zeros(self.m.shape + (3,)) if U is None else np.where(valid[..., None], U, 0.0)
 
     def fill(self, mass, T):
         """Every column filled with `mass` [kg] (per column) of metal at temperature T [K] (per column), in equal cells of
@@ -235,6 +248,7 @@
             self.m, self.H, self.T = pad(self.m, 0.0), pad(self.H, 0.0), pad(self.T, 300.0)
             if self.S is not None:
                 self.S = pad(self.S, 0.0)
+            self.U = np.concatenate([self.U, np.zeros((self.n, extra, 3))], axis=1)
 
     # -- split and merge --------------------------------------------------------------------------------------------
     def regrid(self):
@@ -254,7 +268,7 @@
                 self._split(np.flatnonzero(big.any(axis=1)), big)
         raise RuntimeError("the skin cells did not settle in {} split-and-merge passes".format(REGRID_PASSES))
 
-    def _put_rows(self, rows, m, H, T, S, count):
+    def _put_rows(self, rows, m, H, T, S, count, U):
         """Write rebuilt cell arrays for `rows` back (temperatures of nan are recomputed from the enthalpies), padding
         or trimming the arrays to the longest column."""
         self._width(m.shape[1])
@@ -270,11 +284,14 @@
         self.count[rows] = count
         if self.S is not None:
             self.S[rows] = np.where(valid, pad(S), 0.0)
+        U = np.concatenate([U, np.zeros((len(rows), K - U.shape[1], 3))], axis=1) if U.shape[1] < K else U
+        self.U[rows] = np.where(valid[..., None], U, 0.0)
         width = int(self.count.max()) if self.n else 1
         if width < K:
             self.m, self.H, self.T = self.m[:, :width], self.H[:, :width], self.T[:, :width]
             if self.S is not None:
                 self.S = self.S[:, :width]
+            self.U = self.U[:, :width]
 
     def _merge(self, rows, small):
         """Each run of small cells joins the cell beneath that ends it; a run at the bottom with no cell to end it joins
@@ -303,10 +320,13 @@
             mS = np.zeros(n * K)
             np.add.at(mS, flat, (m0 * self.S[rows])[valid])
             S = np.divide(mS, m, out=np.zeros_like(m), where=m > 0.0).reshape(n, K)
-        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups)
+        mU = np.zeros((n * K, 3))
+        np.add.at(mU, flat, (m0[..., None] * self.U[rows])[valid])
+        U = np.divide(mU, m[:, None], out=np.zeros_like(mU), where=m[:, None] > 0.0).reshape(n, K, 3)
+        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups, U)
 
     def _split(self, rows, big):
-        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T and S)."""
+        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T, S and velocity)."""
         valid = self.valid[rows]
         reps = np.where(valid, 1 + big[rows], 0)
         n, K = reps.shape
@@ -316,6 +336,7 @@
         r_all = np.repeat(np.arange(n)[:, None], K, axis=1)
         m, H, T = np.zeros((n, K2)), np.zeros((n, K2)), np.zeros((n, K2))
         S = np.zeros((n, K2)) if self.S is not None else None
+        U = np.zeros((n, K2, 3))
         share = 1.0 / reps.clip(1)
         for part in (0, 1):
             sel = valid & (part < reps)
@@ -325,27 +346,32 @@
             T[r, c] = self.T[rows][sel]                                    # halves keep the specific enthalpy
             if S is not None:
                 S[r, c] = self.S[rows][sel]
-        self._put_rows(rows, m, H, T, S, count)
+            U[r, c] = self.U[rows][sel]
+        self._put_rows(rows, m, H, T, S, count, U)
 
     # -- the top and the bottom -------------------------------------------------------------------------------------
     def take_top(self, k):
-        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column."""
+        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column, and leaves the momentum
+        that went with them in `momentum_out`."""
         k = np.asarray(k, dtype=np.int64)
         n, K = self.m.shape
         top = np.arange(K)[None, :] < k[:, None]
         mass, H = (self.m * top).sum(axis=1), (self.H * top).sum(axis=1)
+        self.momentum_out = ((self.m * top)[..., None] * self.U).sum(axis=1)
         if k.any():
             src = np.minimum(np.arange(K)[None, :] + k[:, None], K - 1)
             shift = lambda a: np.take_along_axis(a, src, axis=1)
             S = shift(self.S) if self.S is not None else None
-            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T))
+            U = np.take_along_axis(self.U, src[..., None], axis=1)
+            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T), U)
         return mass, H
 
-    def add_bottom(self, dm, h, s=1.0):
+    def add_bottom(self, dm, h, s=1.0, u=None):
         """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
         a new cell); returns the enthalpy added. With S tracked, the metal carries the sensitivity s per column
         (default 1: it is at the base temperature and follows it; 0 for metal whose enthalpy does not depend on the
-        base temperature, such as rigid metal drawn at T_rigid's enthalpy from a base above T_rigid)."""
+        base temperature, such as rigid metal drawn at T_rigid's enthalpy from a base above T_rigid). It arrives with
+        the velocity u (n, 3) [m/s], at rest by default: metal from the 3D model is rigid."""
         dm, h = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,)), np.broadcast_to(np.asarray(h, dtype=float), (self.n,))
         s = np.broadcast_to(np.asarray(s, dtype=float), (self.n,))
         some = dm > 0.0
@@ -359,16 +385,20 @@
         self.H[rows, b] += dm[rows] * h[rows]
         if self.S is not None:
             self.S[rows, b] = (m_old * self.S[rows, b] + dm[rows] * s[rows]) / self.m[rows, b]
+        u = np.zeros((self.n, 3)) if u is None else np.broadcast_to(np.asarray(u, dtype=float), (self.n, 3))
+        self.U[rows, b] = (m_old[:, None] * self.U[rows, b] + dm[rows, None] * u[rows]) / self.m[rows, b][:, None]
         self.count[rows] = np.maximum(self.count[rows], 1)
         self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])
         return np.where(some, dm * h, 0.0)
 
     def take_bottom(self, dm):
-        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy (a partly
-        taken cell keeps its specific enthalpy); returns the enthalpy removed per column."""
+        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy and velocity
+        (a partly taken cell keeps both); returns the enthalpy removed per column, and leaves the momentum that went in
+        `momentum_out`."""
         dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
         if np.any(dm > self.column_mass() * (1.0 + 1e-12) + 1e-300):
             raise ValueError("cannot take more metal from a skin's bottom than the skin holds")
+        self.momentum_out = np.zeros((self.n, 3))
         if not np.any(dm > 0.0):
             return np.zeros(self.n)
         n, K = self.m.shape
@@ -377,10 +407,11 @@
         taken = np.clip(dm[:, None] - below, 0.0, self.m)                 # how much of each cell goes
         frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
         H_taken = (self.H * frac).sum(axis=1)
+        self.momentum_out = (taken[..., None] * self.U).sum(axis=1)
         m, H = self.m - taken, self.H * (1.0 - frac)
         keep = valid & (m > 0.0)
         count = keep.sum(axis=1)                                          # the kept cells are the top ones
-        self._install(m, H, count, self.S, self.T)                        # a partly taken cell keeps its T
+        self._install(m, H, count, self.S, self.T, self.U)                # a partly taken cell keeps its T and velocity
         return H_taken
 
     def exchange_bottom(self, dm, dH):
```

- [ ] **Step 5: The surface operators** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -13,7 +13,17 @@
 to make it negative: a negative conductance would break the explicit step's positivity (and, since sub-plan 21 makes the
 sideways conduction implicit, the definiteness its conjugate gradients need), and it is rare (3 of 6 756 edges on the
 4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
-The surface gradient and curvature operators the flow needs are added in sub-plan 21, the first stage that uses them."""
+
+The flow (sub-plan 21) reads two surface operators on the same graph. `gradient` averages the gradients of a field's
+linear interpolation over each corner's triangles, by area, and projects the result onto the corner's tangent plane.
+`curvature` is the mean curvature (the sum of the principal curvatures) from the cotangent Laplacian of the position
+over each corner's mixed Voronoi area (Meyer et al. 2003). That area is exact only where no touching triangle is
+obtuse, so a corner touching one takes the mean of its neighbours that touch none (33 of the 4 mm test sphere's 2 254
+corners and 48 of the 2 mm sphere's 9 417, which the mixed area alone left off by 15 and 3 % on average and by up to 176
+and 23 %). On these 50 mm-radius spheres the curvature is then within 8.4·10^-5 of 2/R at every corner (1.3·10^-5 on
+the 2 mm sphere), where the barycentric area had scattered 7.6 and 3.9 %, and the gradient of a linear field within
+0.09 % of its tangential part on average (0.02 % on the 2 mm sphere), 0.8 % at the test sphere's one irregular corner
+(measured 2026-10-18)."""
 from dataclasses import dataclass
 
 import numpy as np
@@ -69,6 +79,62 @@
         self.conduct_weight = np.maximum(self.cot_weight, 0.0)
         self.face_length = face
 
+    def gradient(self, f, points=None):
+        """The tangential gradient at every corner of a per-corner field f (n,) [f per m] on the surface through
+        `points` (default: the columns' own corners): the gradients of f's linear interpolation on the touching
+        triangles, averaged by area, projected onto the corner's tangent plane. Exact for a field linear on a flat
+        surface."""
+        p = self.points if points is None else np.asarray(points, dtype=float)
+        normal = self.normal if points is None else corner_frame(p, self.faces, np.arange(self.n))[0]
+        A, unit, _ = triangle_geometry(p, self.faces)
+        a, b, c = (p[self.faces[:, v]] for v in range(3))
+        fa, fb, fc = (np.asarray(f, dtype=float)[self.faces[:, v]] for v in range(3))
+        g = (fa[:, None] * np.cross(unit, c - b) + fb[:, None] * np.cross(unit, a - c)
+             + fc[:, None] * np.cross(unit, b - a)) / (2.0 * A[:, None])
+        G, W = np.zeros((self.n, 3)), np.zeros(self.n)
+        for v in range(3):
+            np.add.at(G, self.faces[:, v], g * A[:, None])
+            np.add.at(W, self.faces[:, v], A)
+        G /= W[:, None]
+        return G - np.einsum("ij,ij->i", G, normal)[:, None] * normal
+
+    def curvature(self, points=None):
+        """The mean curvature, the sum of the principal curvatures [1/m], at every corner of the surface through
+        `points` (default: the columns' own corners), positive where the surface bends away from its outward normal (2/R
+        on a sphere): the cotangent Laplacian of the position along the corner's normal over its mixed Voronoi area. A
+        corner touching an obtuse triangle takes the mean of its neighbours that touch none, or, where every neighbour
+        touches one, of those that took a value in the pass before (module docstring)."""
+        p = self.points if points is None else np.asarray(points, dtype=float)
+        normal = self.normal if points is None else corner_frame(p, self.faces, np.arange(self.n))[0]
+        A, _, cot = triangle_geometry(p, self.faces)
+        lap, area = np.zeros((self.n, 3)), np.zeros(self.n)
+        obtuse = np.zeros(self.n, dtype=bool)
+        blunt = (cot < 0.0).any(axis=1)
+        for v in range(3):
+            i, j, k = self.faces[:, v], self.faces[:, (v + 1) % 3], self.faces[:, (v + 2) % 3]
+            cj, ck = cot[:, (v + 1) % 3], cot[:, (v + 2) % 3]
+            dj, dk = p[j] - p[i], p[k] - p[i]
+            np.add.at(lap, i, 0.5 * (ck[:, None] * dj + cj[:, None] * dk))
+            voronoi = 0.125 * (np.einsum("ij,ij->i", dj, dj) * ck + np.einsum("ij,ij->i", dk, dk) * cj)
+            np.add.at(area, i, np.where(cot[:, v] < 0.0, 0.5 * A, np.where(blunt, 0.25 * A, voronoi)))
+            np.logical_or.at(obtuse, i, blunt)
+        kappa = -np.einsum("ij,ij->i", lap, normal) / area
+        i, j = self.edges[:, 0], self.edges[:, 1]
+        known = ~obtuse
+        while not known.all():                     # outward from the corners whose own value is exact, a ring a pass
+            good = known.astype(float)
+            total, count = np.zeros(self.n), np.zeros(self.n)
+            np.add.at(total, i, kappa[j] * good[j])
+            np.add.at(total, j, kappa[i] * good[i])
+            np.add.at(count, i, good[j])
+            np.add.at(count, j, good[i])
+            fill = ~known & (count > 0.0)
+            if not fill.any():
+                break
+            kappa[fill] = total[fill] / count[fill]
+            known = known | fill
+        return kappa
+
     def top_points(self, thickness):
         """The top surface: each corner moved out along its normal by its column's thickness [m]."""
         return self.points + np.asarray(thickness, dtype=float)[:, None] * self.normal
@@ -96,6 +162,17 @@
         return self._incidence
 
 
+def triangle_geometry(points, faces):
+    """Per triangle: its area [m2], its unit normal (outward for outward-oriented faces) and the cotangent of the angle
+    at each of its three corners, (nt, 3)."""
+    a, b, c = (points[faces[:, v]] for v in range(3))
+    n2 = np.cross(b - a, c - a)
+    twice = np.linalg.norm(n2, axis=1)
+    legs = ((b - a, c - a), (c - b, a - b), (a - c, b - c))
+    cot = np.stack([np.einsum("ij,ij->i", u, w) / np.linalg.norm(np.cross(u, w), axis=1) for u, w in legs], axis=1)
+    return 0.5 * twice, n2 / twice[:, None], cot
+
+
 @dataclass
 class CornerSurface:
     """The gas-facing surface of a corner-column body as the coupled run reads it: one 'patch' per corner (its top dual
```

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

```python
"""The viscosity of liquid and slurry AA7075 (sub-plan 21 Task 3; the Step 4 plan's section 4): Li et al.'s slurry law,
its shear rate held at 367 1/s, the bridge to the liquid at the liquidus, and rigid metal at half solid."""
import numpy as np
import pytest

from reentry_model import material, rheology

MU_L = 1.3e-3


def scheil():
    return material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])


def test_at_rest_li_et_al_runs_from_0_87_to_5_6_pa_s_670_to_4300_times_the_liquid():
    assert rheology.li_viscosity(0.0, 0.0) == pytest.approx(0.871)
    assert rheology.li_viscosity(0.5, 0.0) == pytest.approx(0.871 * np.exp(0.5 * 3.7311))         # 5.63 Pa s
    assert rheology.li_viscosity(0.0, 0.0) / MU_L == pytest.approx(670.0)
    assert rheology.li_viscosity(0.5, 0.0) / MU_L == pytest.approx(4328.0, rel=1e-3)
    # the published form at a measured point: f_s = 0.1, 61 1/s
    assert rheology.li_viscosity(0.1, 61.0) == pytest.approx((0.871 - 0.00849 * 61.0 ** 0.74924) * np.exp(0.37311))


def test_the_shear_rate_is_held_at_367_per_second_where_the_fit_still_holds():
    f_s = np.linspace(0.1, 0.5, 5)
    np.testing.assert_array_equal(rheology.li_viscosity(f_s, 450.0), rheology.li_viscosity(f_s, 367.0))
    assert np.all(rheology.li_viscosity(f_s, 367.0) > 0.0)                    # the bracket's zero is at 483 1/s
    rates = np.array([0.0, 61.0, 150.0, 367.0])
    assert np.all(np.diff(rheology.li_viscosity(0.3, rates)) < 0.0)           # shear-thinning


def test_the_bridge_runs_log_linearly_from_the_slurry_to_the_liquid_at_the_liquidus():
    f_b = 0.0975
    top = rheology.li_viscosity(f_b, 100.0)
    eta = lambda T, f_s: float(rheology.viscosity(T, f_s, 100.0, MU_L, f_b))
    assert eta(906.4, f_b) == pytest.approx(top)
    assert eta(907.2, 0.05) == pytest.approx(np.sqrt(top * MU_L))             # half way: the geometric mean
    assert eta(908.0, 0.0) == MU_L and eta(950.0, 0.0) == MU_L


def test_metal_half_solid_or_more_is_rigid():
    assert np.isinf(rheology.viscosity(895.0, 0.5, 50.0, MU_L, 0.0975))
    assert np.isinf(rheology.viscosity(850.0, 0.9, 50.0, MU_L, 0.0975))
    assert rheology.viscosity(895.2, 0.4999, 50.0, MU_L, 0.0975) == pytest.approx(rheology.li_viscosity(0.4999, 50.0))


def test_on_aa7075_s_scheil_curve_the_viscosity_falls_without_a_break_from_the_rigid_front_to_the_liquid():
    mat = scheil()
    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    assert f_b == pytest.approx(0.0975, abs=1e-4)                              # Scheil's curve at 906.4 K
    T = np.linspace(mat.T_rigid + 1e-3, 910.0, 4001)
    eta = rheology.viscosity(T, 1.0 - mat.liquid_fraction(T), 100.0, mat.liquid.mu, f_b)
    assert np.all(np.isfinite(eta)) and np.all(np.diff(eta) <= 1e-12 * eta[:-1])
    for T_edge in (rheology.BRIDGE_LOW, rheology.BRIDGE_HIGH):                 # no jump where the bridge begins or ends
        Ts = np.array([T_edge - 1e-9, T_edge + 1e-9])
        lo, hi = rheology.viscosity(Ts, 1.0 - mat.liquid_fraction(Ts), 100.0, mat.liquid.mu, f_b)
        assert hi == pytest.approx(lo, rel=1e-6)
```

- [ ] **Step 2: Append the drag tests** to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -544,3 +544,46 @@
     dm, u = np.full(sk.n, 1e-7), rng.normal(0.0, 0.2, (sk.n, 3))
     sk.add_bottom(dm, np.full(sk.n, float(sk.material.enthalpy(700.0))), u=u)
     np.testing.assert_allclose(sk.column_momentum(), P + dm[:, None] * u, rtol=1e-12, atol=1e-18)
+
+
+# -- viscosity (sub-plan 21 Task 3) -------------------------------------------------------------------------------------
+def column_at(T_cells, mat=None):
+    """One column of 10 um cells at the given temperatures, top first."""
+    mat = mat or scheil()
+    T_cells = np.asarray(T_cells, dtype=float)
+    sk = skin.Skins(mat, np.array([4e-6]), skin.SkinSettings(**CELLS_04))
+    m = np.full((1, len(T_cells)), mat.rho * 4e-6 * 1e-5)
+    sk._install(m, m * mat.enthalpy(T_cells)[None, :], np.array([len(T_cells)]))
+    return sk
+
+
+def test_two_liquid_cells_drag_on_each_other_with_the_liquid_s_viscosity_over_their_spacing():
+    sk = column_at([950.0, 940.0, 930.0])
+    D = skin.interface_drag(sk, 200.0)
+    dz = sk.dz()[0]
+    np.testing.assert_allclose(D[0], 1.3e-3 / (0.5 * (dz[:-1] + dz[1:])), rtol=1e-12)
+
+
+def test_across_a_viscosity_jump_the_two_half_cells_act_in_series():
+    """Liquid over slurry: the drag is 1 / (dz_l / 2 mu_l + dz_s / 2 eta_s), which the liquid's half dominates -- the
+    stress is continuous through the jump, and the profile Task 4 builds on it has the exact kink."""
+    sk = column_at([920.0, 900.0])
+    eta = skin.cell_viscosity(sk, 150.0)[0]
+    dz = sk.dz()[0]
+    assert eta[0] == 1.3e-3 and eta[1] > 500 * eta[0]
+    assert skin.interface_drag(sk, 150.0)[0, 0] == pytest.approx(1.0 / (0.5 * dz[0] / eta[0] + 0.5 * dz[1] / eta[1]))
+
+
+def test_a_flowing_cell_drags_on_the_rigid_cell_beneath_over_its_own_half_thickness_alone():
+    sk = column_at([905.0, 890.0, 880.0])                         # slurry, then two rigid cells (T_rigid 895.1 K)
+    eta = skin.cell_viscosity(sk, 50.0)[0]
+    assert np.isfinite(eta[0]) and np.all(np.isinf(eta[1:]))
+    D = skin.interface_drag(sk, 50.0)[0]
+    assert D[0] == pytest.approx(eta[0] / (0.5 * sk.dz()[0, 0]))  # no slip at the rigid front
+    assert D[1] == 0.0                                            # rigid on rigid: nothing moves
+
+
+def test_metal_that_does_not_melt_never_flows():
+    mat = material_mod.Material.from_drama_json("AA7075_nomelt")
+    sk = column_at([1500.0, 1400.0], mat)
+    assert np.all(np.isinf(skin.cell_viscosity(sk, 10.0))) and np.all(skin.interface_drag(sk, 10.0) == 0.0)
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_rheology.py
  tests/test_reentry_model_skin.py -q`. Expected: `ModuleNotFoundError: No module named 'reentry_model.rheology'`.

- [ ] **Step 4: Create `reentry_model/rheology.py`:**

```python
"""The viscosity of liquid and slurry AA7075 (the Step 4 plan's section 4, docs/superpowers/plans/2026-09-27-deformation-
ring-shedding.md; the layered-skin design section 6.4). Shared physics: it imports neither architecture's modules, so
that the later move into reentry_model/shared/ stays mechanical (sub-plan 21's amendment of 2026-10-12).

- Liquid, at and above the liquidus: the material's own viscosity (1.3 mPa s for AA7075).
- Bridge, BRIDGE_LOW to BRIDGE_HIGH (906.4-908 K): log-linear in the temperature from Li et al.'s value at BRIDGE_LOW,
  where Scheil's curve stands at 10 % solid, to the liquid at the liquidus (labelled: no data).
- Slurry, from the bridge down to half solid: Li Yageng et al. (2014, China Foundry 11, 79-84),
  eta = [0.871 - 0.00849 gamma^0.74924] exp(3.7311 f_s) Pa s, measured at f_s = 0.1-0.5 and shear rates gamma of
  61-490 1/s on continuously sheared slurry; gamma is held at no more than GAMMA_MAX = 367 1/s, the highest rate at
  which the fit still matches the data (its bracket would reach zero at 483 1/s). Between zero and 61 1/s it is
  extrapolated, and there it is the runny side of the data (1.4-1.9 times low at 61 1/s; the Step 4 plan's fact 10).
  From 0.87 Pa s (f_s = 0) to 5.6 Pa s (f_s = 0.5) at rest, 670 to 4 300 times the liquid.
- Rigid, at or above half solid: infinite. A rigid cell does not flow; it is anchored (spec 6.4)."""
import numpy as np

LI_A = 0.871                       # Pa s
LI_B = 0.00849                     # Pa s^(1 + LI_N)
LI_N = 0.74924
LI_C = 3.7311
GAMMA_MAX = 367.0                  # 1/s: Li et al.'s law is held at this shear rate
BRIDGE_LOW, BRIDGE_HIGH = 906.4, 908.0     # K: the log-linear bridge from the slurry to the liquid
RIGID_SOLID_FRACTION = 0.5         # at or above it a cell is rigid (material.RIGID_LIQUID_FRACTION's complement)


def li_viscosity(f_s, gamma):
    """Li et al.'s slurry viscosity [Pa s] at solid fraction f_s and shear rate gamma [1/s], the rate held within
    [0, GAMMA_MAX]."""
    g = np.clip(np.asarray(gamma, dtype=float), 0.0, GAMMA_MAX)
    return (LI_A - LI_B * g ** LI_N) * np.exp(LI_C * np.asarray(f_s, dtype=float))


def viscosity(T, f_s, gamma, mu_liquid, f_s_bridge):
    """The viscosity [Pa s] of metal at temperature T [K] and solid fraction f_s, sheared at gamma [1/s] (module
    docstring); `mu_liquid` [Pa s] is the liquid's and `f_s_bridge` the material's solid fraction at BRIDGE_LOW (0.098 on
    AA7075's Scheil curve). Infinite where the metal is rigid."""
    T, f_s = np.asarray(T, dtype=float), np.asarray(f_s, dtype=float)
    top = li_viscosity(f_s_bridge, gamma)
    w = np.clip((T - BRIDGE_LOW) / (BRIDGE_HIGH - BRIDGE_LOW), 0.0, 1.0)
    bridge = np.exp((1.0 - w) * np.log(top) + w * np.log(mu_liquid))
    eta = np.where(T >= BRIDGE_HIGH, mu_liquid, np.where(T > BRIDGE_LOW, bridge, li_viscosity(f_s, gamma)))
    return np.where(f_s >= RIGID_SOLID_FRACTION, np.inf, eta)
```

- [ ] **Step 5: The cells' viscosity and the drag between them** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -32,6 +32,7 @@
 
 import numpy as np
 
+from . import rheology
 from .thermal import SIGMA_SB
 
 SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (spec section 6.2)
@@ -778,6 +779,39 @@
         return np.where(self.valid, np.diff(_interp_apply(self.geometry["to_cells"], Qcum), axis=1), 0.0)
 
 
+# -- viscosity (sub-plan 21 Task 3) ------------------------------------------------------------------------------------
+def cell_viscosity(sk, gamma):
+    """Every cell's viscosity [Pa s], (n, K), sheared at gamma [1/s] (per cell, or one for all): `rheology.viscosity` at
+    the cell's temperature and Scheil solid fraction; infinite in a rigid cell (half solid or more) and beyond each
+    column's cells, and everywhere for a material that does not melt."""
+    mat = sk.material
+    if not mat.melts or mat.liquid is None:
+        return np.full(sk.m.shape, np.inf)
+    f_s = 1.0 - sk.props.liquid_fraction(sk.T)
+    f_s_bridge = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    eta = rheology.viscosity(sk.T, f_s, gamma, mat.liquid.mu, f_s_bridge)
+    return np.where(sk.valid, eta, np.inf)
+
+
+def interface_drag(sk, gamma):
+    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1):
+    1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), both cells' viscosities taken at the interface's shear rate gamma [1/s]
+    (n, K-1) -- the two half-cells in series, as the shear stress's continuity across the interface requires, the way
+    `_conduct` combines conductivities. A rigid cell resists without limit and adds nothing to the sum, so a flowing cell
+    drags on the rigid cell beneath it over its own half-thickness alone: no slip at the rigid front. Zero between two
+    rigid cells and across a missing one."""
+    if sk.m.shape[1] < 2:
+        return np.zeros((sk.n, 0))
+    g = np.asarray(gamma, dtype=float) * np.ones((sk.n, sk.m.shape[1] - 1))
+    eta_up = cell_viscosity(sk, np.pad(g, ((0, 0), (0, 1)), mode="edge"))[:, :-1]
+    eta_dn = cell_viscosity(sk, np.pad(g, ((0, 0), (1, 0)), mode="edge"))[:, 1:]
+    dz = sk.dz()
+    resist = 0.5 * dz[:, :-1] / eta_up + 0.5 * dz[:, 1:] / eta_dn
+    link = sk.valid[:, 1:]
+    flows = np.isfinite(eta_up) | np.isfinite(eta_dn)
+    return np.where(link & flows, 1.0 / np.where(link & flows, resist, 1.0), 0.0)
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

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

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -4,7 +4,7 @@
 import pytest
 
 from reentry_model import material as material_mod
-from reentry_model import skin
+from reentry_model import rheology, skin
 
 
 def scheil():
@@ -587,3 +587,145 @@
     mat = material_mod.Material.from_drama_json("AA7075_nomelt")
     sk = column_at([1500.0, 1400.0], mat)
     assert np.all(np.isinf(skin.cell_viscosity(sk, 10.0))) and np.all(skin.interface_drag(sk, 10.0) == 0.0)
+
+
+# -- momentum through the depth (sub-plan 21 Task 4) --------------------------------------------------------------------
+MU_L = 1.3e-3
+STEADY = 1e9                       # s: a backward-Euler step this long gives the steady profile
+
+
+def heights(sk, k):
+    """The first k cells' centres' height [m] above the rigid front beneath them."""
+    dz = sk.dz()[0][:k]
+    return dz.sum() - (np.cumsum(dz) - 0.5 * dz)
+
+
+def books_close(sk, P0, book):
+    """The impulses add up to the momentum change, to round-off of the largest of them (in a steady step they are
+    large and nearly cancel)."""
+    parts = (book.push, book.shear, book.wall, P0, sk.column_momentum())
+    scale = max(np.abs(x).max() for x in parts) or 1e-30
+    return np.abs(sk.column_momentum() - P0 - book.push - book.shear - book.wall).max() <= 1e-12 * scale
+
+
+def test_shear_alone_gives_a_linear_profile():
+    """Liquid sheared at its top by 1 Pa over rigid metal: u = tau s / mu at every cell centre, s the height above the
+    rigid front (Couette), exactly; the drag to the rigid front over the bottom liquid cell's half-thickness is no slip."""
+    sk = column_at([950.0] * 30 + [880.0] * 10)
+    P0 = sk.column_momentum()
+    book = skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    np.testing.assert_allclose(sk.U[0, :30, 0], heights(sk, 30) / MU_L, rtol=1e-9)
+    assert np.all(sk.U[0, :, 1:] == 0.0) and np.all(sk.U[0, 30:] == 0.0) and book.passes == 1
+    assert books_close(sk, P0, book)
+
+
+def test_a_push_alone_gives_a_half_parabola_converging_at_second_order():
+    """A push of 10 kPa/m on a free liquid layer over rigid metal: u = (G / mu)(h s - s^2 / 2). Measured 2026-10-19:
+    2.8e-4 of the top's speed off at 30 cells, 6.9e-5 at 60."""
+    def error(N):
+        sk = column_at([950.0] * N + [880.0] * 10)
+        skin.momentum_step(sk, STEADY, np.array([1e4, 0.0, 0.0]), np.zeros(3), sk.area)
+        s, h = heights(sk, N), sk.dz()[0][:N].sum()
+        exact = 1e4 / MU_L * (h * s - 0.5 * s * s)
+        return np.abs(sk.U[0, :N, 0] - exact).max() / exact.max()
+    e30, e60 = error(30), error(60)
+    assert e30 < 5e-4 and e30 / e60 == pytest.approx(4.0, rel=0.1)
+
+
+def test_liquid_over_slurry_joins_two_linear_profiles_with_the_kink_at_the_viscosity_jump():
+    """1 Pa of gas shear on 15 liquid cells over 15 of slurry at 900 K (37 % solid) over rigid metal: the stress is 1 Pa
+    all through, so the slurry shears at 1 Pa over its own viscosity at that rate (3.46 Pa s, 0.29 1/s) and the liquid at
+    1 Pa / 1.3 mPa s, the profiles joining at the jump. Each cell at its own shear rate makes it exact; the fixed point
+    takes two passes."""
+    sk = column_at([950.0] * 15 + [900.0] * 15 + [880.0] * 10)
+    mat = sk.material
+    f_s = 1.0 - float(mat.liquid_fraction(np.array([900.0]))[0])
+    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    gamma = 0.0
+    for _ in range(50):
+        gamma = 1.0 / float(rheology.viscosity(900.0, f_s, gamma, MU_L, f_b))
+    book = skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    s = heights(sk, 30)
+    h_s = sk.dz()[0][15:30].sum()
+    np.testing.assert_allclose(sk.U[0, 15:30, 0], gamma * s[15:], rtol=1e-4)
+    np.testing.assert_allclose(sk.U[0, :15, 0], gamma * h_s + (s[:15] - h_s) / MU_L, rtol=1e-6)
+    assert book.passes <= skin.SHEAR_PASSES
+
+
+def test_a_layer_spun_up_from_rest_by_a_sudden_shear_follows_stokes_first_problem():
+    """1 Pa suddenly applied at the top of 0.9 mm of liquid at rest, for 18.5 ms (a diffusion length sqrt(nu t) of
+    100 um): u(y, t) = (2 tau / mu) sqrt(nu t) ierfc(y / 2 sqrt(nu t)), y the depth. The top cell within 1e-3 at 200
+    steps, the error halving with the step (backward Euler; 1.3e-3, 6.6e-4, 3.3e-4 at 100, 200, 400 steps, measured
+    2026-10-19); the momentum is the shear's impulse exactly."""
+    from scipy.special import erfc
+    nu, T_end = MU_L / 2400.0, 0.0185
+    L = np.sqrt(nu * T_end)
+
+    def run(steps):
+        sk = column_at([950.0] * 80 + [880.0] * 5)
+        for _ in range(steps):
+            skin.momentum_step(sk, T_end / steps, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+        dz = sk.dz()[0][:80]
+        y = np.cumsum(dz) - 0.5 * dz
+        exact = 2.0 / MU_L * L * (np.exp(-(y / (2 * L)) ** 2) / np.sqrt(np.pi) - y / (2 * L) * erfc(y / (2 * L)))
+        return sk, abs(sk.U[0, 0, 0] / exact[0] - 1.0), np.abs(sk.U[0, :80, 0] - exact).max() / exact[0]
+    sk, top200, profile = run(200)
+    assert top200 < 1e-3 and profile < 1e-3
+    assert run(100)[1] / top200 == pytest.approx(2.0, rel=0.1)
+    assert sk.column_momentum()[0, 0] == pytest.approx(1.0 * sk.area[0] * T_end, rel=1e-12)
+
+
+def slurry_rate(sk, T, tau):
+    """The steady shear rate [1/s] of slurry at temperature T under the stress tau [Pa]: tau over its own viscosity at
+    that rate, by fixed point."""
+    mat = sk.material
+    f_s = 1.0 - float(mat.liquid_fraction(np.array([T]))[0])
+    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    gamma = 0.0
+    for _ in range(50):
+        gamma = tau / float(rheology.viscosity(T, f_s, gamma, MU_L, f_b))
+    return gamma
+
+
+def test_a_floating_crust_slides_on_the_slurry_beneath_it_as_one_piece():
+    """Asha, 2026-10-20, option B: a crust (three rigid cells at 880 K) on 15 cells of slurry at 900 K over the rigid
+    metal joined to the substrate, 1 Pa of gas shear on its top: the stress is 1 Pa all through, the slurry shears in
+    Couette's line from the anchored front, and the crust rides its top as one piece, at the slurry's speed there."""
+    sk = column_at([880.0] * 3 + [900.0] * 15 + [880.0] * 5)
+    gamma = slurry_rate(sk, 900.0, 1.0)
+    book = skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    s = heights(sk, 18)
+    h_s = sk.dz()[0][3:18].sum()
+    np.testing.assert_allclose(sk.U[0, 3:18, 0], gamma * s[3:], rtol=1e-4)
+    np.testing.assert_allclose(sk.U[0, :3, 0], gamma * h_s, rtol=1e-4)
+    assert np.ptp(sk.U[0, :3, 0]) < 1e-7 * sk.U[0, 0, 0]                    # one piece
+    assert np.all(sk.U[0, 18:] == 0.0) and np.all(book.shear[:, 0] > 0.0)
+
+
+def test_the_rigid_metal_joined_to_the_substrate_is_held_and_the_momentum_books_close():
+    """A crust over slurry and liquid over the anchored rigid metal, every cell moving at the start and pushed along the
+    surface: the anchored cells stop -- their momentum booked to the wall -- while the crust moves and takes the gas
+    shear, and the impulses add up to the momentum change."""
+    sk = column_at([880.0] + [920.0] * 10 + [902.0] * 10 + [880.0] * 5)
+    rng = np.random.default_rng(3)
+    sk.U = rng.normal(0.0, 0.05, sk.U.shape)
+    P0 = sk.column_momentum()
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    assert np.array_equal(skin.anchored_cells(sk, flowing)[0], np.arange(26) >= 21)
+    book = skin.momentum_step(sk, 1e-3, np.array([2e4, -1e4, 0.0]), np.array([3.0, 0.0, 0.0]), sk.area)
+    assert np.all(sk.U[0, 21:] == 0.0) and np.all(sk.U[0, :21] != 0.0)
+    assert book.shear[0, 0] == pytest.approx(3.0 * sk.area[0] * 1e-3)
+    assert books_close(sk, P0, book)
+
+
+def test_a_column_whose_bottom_cell_flows_anchors_nothing():
+    """Slurry down to the base under a crust: no rigid cell is joined to the substrate, so the crust floats and the
+    slurry's bottom cell meets the rigid 3D metal beneath the base with no slip."""
+    sk = column_at([880.0] * 2 + [900.0] * 10)
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    assert not skin.anchored_cells(sk, flowing).any()
+    gamma = slurry_rate(sk, 900.0, 1.0)
+    skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    dz = sk.dz()[0]
+    s = dz[:12].sum() - (np.cumsum(dz[:12]) - 0.5 * dz[:12])                # heights above the base
+    np.testing.assert_allclose(sk.U[0, 2:12, 0], gamma * s[2:], rtol=1e-4)
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `AttributeError: module 'reentry_model.skin' has no attribute 'momentum_step'`.

- [ ] **Step 3: The momentum sub-step, and the drag at each cell's own shear rate,** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -793,23 +793,141 @@
     return np.where(sk.valid, eta, np.inf)
 
 
-def interface_drag(sk, gamma):
-    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1):
-    1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), both cells' viscosities taken at the interface's shear rate gamma [1/s]
-    (n, K-1) -- the two half-cells in series, as the shear stress's continuity across the interface requires, the way
-    `_conduct` combines conductivities. A rigid cell resists without limit and adds nothing to the sum, so a flowing cell
-    drags on the rigid cell beneath it over its own half-thickness alone: no slip at the rigid front. Zero between two
-    rigid cells and across a missing one."""
+def series_drag(sk, eta):
+    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1), from the
+    cells' viscosities eta [Pa s] (n, K): 1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), the two half-cells in series, as
+    the shear stress's continuity across the interface requires, the way `_conduct` combines conductivities. A rigid
+    cell resists without limit and adds nothing to the sum, so a flowing cell drags on the rigid cell beneath it over
+    its own half-thickness alone: no slip at the rigid front. Zero between two rigid cells and across a missing one."""
     if sk.m.shape[1] < 2:
         return np.zeros((sk.n, 0))
-    g = np.asarray(gamma, dtype=float) * np.ones((sk.n, sk.m.shape[1] - 1))
-    eta_up = cell_viscosity(sk, np.pad(g, ((0, 0), (0, 1)), mode="edge"))[:, :-1]
-    eta_dn = cell_viscosity(sk, np.pad(g, ((0, 0), (1, 0)), mode="edge"))[:, 1:]
     dz = sk.dz()
-    resist = 0.5 * dz[:, :-1] / eta_up + 0.5 * dz[:, 1:] / eta_dn
-    link = sk.valid[:, 1:]
-    flows = np.isfinite(eta_up) | np.isfinite(eta_dn)
-    return np.where(link & flows, 1.0 / np.where(link & flows, resist, 1.0), 0.0)
+    resist = 0.5 * dz[:, :-1] / eta[:, :-1] + 0.5 * dz[:, 1:] / eta[:, 1:]
+    link = sk.valid[:, 1:] & (np.isfinite(eta[:, :-1]) | np.isfinite(eta[:, 1:]))
+    return np.where(link, 1.0 / np.where(link, resist, 1.0), 0.0)
+
+
+def interface_drag(sk, gamma):
+    """`series_drag` with every cell's viscosity at its own shear rate gamma [1/s] (per cell, or one for all). Each
+    cell, not each interface: across a jump from liquid to slurry the liquid shears some thousand times faster than the
+    slurry, and a rate taken from the velocity difference between the two centres would give the slurry's half-cell the
+    liquid's rate (found while writing Task 4)."""
+    return series_drag(sk, cell_viscosity(sk, gamma))
+
+
+# -- momentum through the depth (sub-plan 21 Task 4) --------------------------------------------------------------------
+SHEAR_PASSES = 5                  # passes of the fixed point on the cells' shear rates at most (spec 6.4)
+SHEAR_TOL = 1.0e-3                # a pass that moves no flowing cell's viscosity by more than this, relatively, ends it
+CRUST_VISCOSITY = 1.0e9           # Pa s: a floating crust's, so that it moves as one piece (some 1e8 times the slurry's)
+
+
+@dataclass
+class MomentumBook:
+    """What a momentum sub-step did, per column: impulses [N s], (n, 3), which add up to the momentum change."""
+    push: np.ndarray         # of the push on the moving cells
+    shear: np.ndarray        # of the gas shear on the top cell, if it moves
+    wall: np.ndarray         # of the anchored rigid metal -- the skin's rigid cells joined to the substrate and the 3D
+                             # model beneath the base -- on the moving cells, and the momentum a cell loses on joining it
+    passes: int              # the fixed point's passes
+
+
+def anchored_cells(sk, flowing):
+    """The rigid cells joined to the substrate, (n, K): each column's run of rigid cells from its bottom cell up, which
+    rests on the rigid 3D model. A column whose bottom cell flows has none."""
+    n, K = sk.m.shape
+    j = np.arange(K)[None, :]
+    from_bottom = np.clip(sk.count[:, None] - 1 - j, 0, K - 1)                 # the j-th cell from the bottom
+    rigid_up = np.take_along_axis(sk.valid & ~flowing, from_bottom, axis=1) & (j < sk.count[:, None])
+    held = np.cumprod(rigid_up, axis=1).sum(axis=1)
+    return sk.valid & (j >= (sk.count - held)[:, None])
+
+
+def momentum_step(sk, dt, push, shear, A_top):
+    """One implicit sub-step dt [s] of every column's momentum through its depth (spec 6.4), each velocity component
+    apart:
+
+        m_i (u_i' - u_i) / dt = A [D_(i-1) (u'_(i-1) - u'_i) + D_i (u'_(i+1) - u'_i)] + G V_i + tau A_top [i the top],
+
+    A the column's area, D the drag between neighbouring cells (`series_drag`), G the push per unit volume [N/m3,
+    (n, 3)] -- every layer of a column feels the same (spec 6.4) -- V_i the cell's volume and tau the gas shear [Pa,
+    (n, 3)] on the top cell over its top area A_top [m2, (n,)]. The rigid metal joined to the substrate (`anchored_cells`)
+    is held at rest; a moving cell drags on it, and a flowing bottom cell on the rigid 3D metal beneath the base, with no
+    slip over its own half-thickness. A crust -- rigid cells with flowing metal between them and the substrate -- floats:
+    it moves as one piece (CRUST_VISCOSITY), dragged by the gas shear on its top, pushed like the rest and dragged by the
+    slurry beneath (Asha, 2026-10-20, option B; the column model has no sideways stress in a solid, so a crust over many
+    columns moves as column-sized plates, as a warm, barely solid crust would). Each flowing cell's viscosity is at its
+    own shear rate, the mean over its two faces of the stress there over its viscosity, found by fixed-point iteration
+    from the start-of-sub-step velocities: at most SHEAR_PASSES passes, ending when no flowing cell's viscosity moves by
+    more than SHEAR_TOL. Updates sk.U; returns a MomentumBook whose impulses add up to the columns' momentum change to
+    round-off."""
+    n, K = sk.m.shape
+    valid = sk.valid
+    rows = np.arange(n)
+    last = np.maximum(sk.count - 1, 0)
+    dz, A, m = sk.dz(), sk.area, sk.m
+    push = np.broadcast_to(np.asarray(push, dtype=float), (n, 3))
+    shear = np.broadcast_to(np.asarray(shear, dtype=float), (n, 3))
+    A_top = np.broadcast_to(np.asarray(A_top, dtype=float), (n,))
+    flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
+    anchored = anchored_cells(sk, flowing)
+    moving = valid & ~anchored                                            # the flowing cells and the floating crust
+    crust = moving & ~flowing
+    P0 = sk.column_momentum()
+    held = ((m * anchored)[..., None] * sk.U).sum(axis=1)                 # what anchored cells carried, now at rest
+    U_old = np.where(moving[..., None], sk.U, 0.0)
+    top_moves = moving[:, 0]
+    force = np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0)
+    force[:, 0] += np.where(top_moves[:, None], shear * A_top[:, None], 0.0)
+    bottom_flows = flowing[rows, last] & (sk.count > 0)
+
+    def viscosities(gamma):
+        return np.where(crust, CRUST_VISCOSITY, cell_viscosity(sk, gamma))
+
+    def drags(eta):
+        Db = np.where(bottom_flows, eta[rows, last] / (0.5 * np.maximum(dz[rows, last], DZ_FLOOR)), 0.0)
+        return series_drag(sk, eta), np.where(np.isfinite(Db), Db, 0.0)
+
+    def solve(eta):
+        D, Db = drags(eta)
+        diag = np.where(valid, m / dt, 1.0)
+        diag[:, :-1] += A[:, None] * D
+        diag[:, 1:] += A[:, None] * D
+        diag[rows, last] += A * Db
+        off = -A[:, None] * D
+        sub, sup = np.where(moving[:, 1:], off, 0.0), np.where(moving[:, :-1], off, 0.0)
+        diag = np.where(moving, diag, 1.0)
+        rhs = np.where(moving[..., None], (m / dt)[..., None] * U_old + force, 0.0)
+        U = np.stack([thomas(sub, diag, sup, rhs[..., c]) for c in range(3)], axis=-1)
+        return np.where(moving[..., None], U, 0.0), D, Db
+
+    def rates(U, D, Db, eta):
+        """Every flowing cell's shear rate: the mean over its faces of the stress there over its viscosity."""
+        stress = np.zeros((n, K + 1))                                     # at every face, the top face first
+        stress[:, 1:K] = D * np.linalg.norm(U[:, :-1] - U[:, 1:], axis=2)
+        stress[:, 0] = np.where(top_moves, np.linalg.norm(shear, axis=1), 0.0)
+        stress[rows, last + 1] = Db * np.linalg.norm(U[rows, last], axis=1)
+        face = 0.5 * (stress[:, :-1] + stress[:, 1:])
+        return np.where(flowing, np.divide(face, eta, out=np.zeros_like(face), where=np.isfinite(eta)), 0.0)
+
+    eta = viscosities(0.0)
+    eta = viscosities(rates(U_old, *drags(eta), eta))                       # at the start-of-sub-step rates
+    for passes in range(1, SHEAR_PASSES + 1):
+        U, D, Db = solve(eta)
+        eta_new = viscosities(rates(U, D, Db, eta))
+        change = np.abs(eta_new[flowing] - eta[flowing]) / eta[flowing]
+        if change.size == 0 or change.max() <= SHEAR_TOL:
+            break
+        eta = eta_new
+    wall = np.zeros((n, 3))
+    if K > 1:
+        down = moving[:, :-1] & anchored[:, 1:]                           # a moving cell on the anchored metal beneath
+        wall -= ((A[:, None] * D * down)[..., None] * U[:, :-1]).sum(axis=1)
+    wall -= (A * Db)[:, None] * U[rows, last]
+    book = MomentumBook(push=dt * np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0).sum(axis=1),
+                        shear=dt * np.where(top_moves[:, None], shear * A_top[:, None], 0.0),
+                        wall=dt * wall - held, passes=passes)
+    sk.U = U
+    return book
 
 
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
```

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

```diff
--- a/tests/test_reentry_model_surface_flow.py
+++ b/tests/test_reentry_model_surface_flow.py
@@ -126,3 +126,15 @@
     assert r.branch == sf.BRANCH_FREE_MOLECULAR and np.all(np.isnan(r.delta_a)) and np.all(r.u_e == 0.0)
     assert np.isnan(r.p_stag) and np.all(r.p_w == r.p_inf) and np.all(r.closure == sf.CLOSURE_COUETTE)
     assert r.tau[1] == pytest.approx(a.freestream.rho * a.V ** 2 * math.sin(theta[1]) * math.cos(theta[1])) and r.tau[3] == 0.0
+
+
+def test_g_is_the_pressure_gradient_along_the_meridian_and_the_deceleration():
+    """Sub-plan 21 Task 5: the surface flow hands the skins the pressure part of G on its own, dp_ds; on the windward
+    face G = -dp_ds - rho_l W sin(theta), with p_w falling from the stagnation point, and leeward both are zero."""
+    state = aero_state_at(80.0)
+    theta = np.radians(np.linspace(0.0, 150.0, 61))
+    r = sf.SurfaceFlow().evaluate(state, theta, 0.05, 2400.0)
+    windward = theta <= 0.5 * math.pi
+    np.testing.assert_allclose(r.G[windward], -r.dp_ds[windward] - 2400.0 * r.deceleration * np.sin(theta[windward]),
+                               rtol=1e-12, atol=1e-9)
+    assert np.all(r.dp_ds[windward][1:-1] < 0.0) and np.all(r.dp_ds[~windward] == 0.0) and np.all(r.G[~windward] == 0.0)
```

- [ ] **Step 2: The push's tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -729,3 +729,89 @@
     dz = sk.dz()[0]
     s = dz[:12].sum() - (np.cumsum(dz[:12]) - 0.5 * dz[:12])                # heights above the base
     np.testing.assert_allclose(sk.U[0, 2:12, 0], gamma * s[2:], rtol=1e-4)
+
+
+# -- the push (sub-plan 21 Task 5) --------------------------------------------------------------------------------------
+V_HAT = np.array([1.0, 0.0, 0.0])
+RHO_L, SIGMA = 2400.0, 0.8
+
+
+def regular_grid(n_side=81, h=1e-3):
+    """Columns on a regular flat square grid of n_side corners a side, h apart, in z = 0, normals +z."""
+    from reentry_model import columns
+    x, y = np.meshgrid(np.arange(n_side) * h, np.arange(n_side) * h, indexing="ij")
+    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
+    idx = np.arange(n_side * n_side).reshape(n_side, n_side)
+    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
+    return columns.Columns(pts, np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])), pts
+
+
+def test_along_the_meridian_is_in_the_tangent_plane_away_from_the_stagnation_point(coarse_sphere_mesh):
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    s = skin.along_meridian(col.normal, V_HAT, np.ones(col.n))
+    side = np.linalg.norm(V_HAT - (col.normal @ V_HAT)[:, None] * col.normal, axis=1) > 1e-6
+    assert np.abs(np.einsum("ij,ij->i", s, col.normal)).max() < 1e-12
+    np.testing.assert_allclose(np.linalg.norm(s[side], axis=1), 1.0, rtol=1e-12)
+    assert np.all(s[side] @ V_HAT < 0.0)                                       # away from the stagnation point
+    poles = skin.along_meridian(np.array([V_HAT, -V_HAT]), V_HAT, np.ones(2))
+    assert np.all(poles == 0.0)
+
+
+def test_the_deceleration_pulls_the_metal_toward_the_nose_everywhere(coarse_sphere_mesh):
+    """30 m/s2 along the flight direction on the 4 mm sphere: rho_l g sin(theta) toward the nose -- the element model's
+    -rho_l W sin(theta) on the windward face, and on the leeward face too, where it pulls the metal toward the equator."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    G, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), 30.0 * V_HAT, V_HAT, RHO_L, SIGMA)
+    sin = np.linalg.norm(V_HAT - (col.normal @ V_HAT)[:, None] * col.normal, axis=1)
+    away = skin.along_meridian(col.normal, V_HAT, np.ones(col.n))
+    np.testing.assert_allclose(np.einsum("ij,ij->i", terms["deceleration"], away), -RHO_L * 30.0 * sin, rtol=1e-9, atol=1e-9)
+    assert np.all(terms["pressure"] == 0.0) and np.abs(terms["thickness"]).max() < 1e-9          # uniform: round-off
+
+
+def test_the_pressure_term_is_the_surface_flow_s_gradient_away_from_the_stagnation_point(coarse_sphere_mesh):
+    """A Newtonian wall pressure falling from the stagnation point, dp/ds = -2 p0 sin cos / R: its push points away
+    from the stagnation point at the surface flow's magnitude -- the input's own, not re-differentiated."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    cos = np.clip(col.normal @ V_HAT, 0.0, 1.0)
+    dp_ds = -2.0 * 11e3 * np.sqrt(1.0 - cos ** 2) * cos / 0.05
+    _, terms = skin.push(col, np.full(col.n, 8e-4), dp_ds, np.zeros(3), V_HAT, RHO_L, SIGMA)
+    away = skin.along_meridian(col.normal, V_HAT, np.ones(col.n))
+    np.testing.assert_allclose(np.einsum("ij,ij->i", terms["pressure"], away), -dp_ds, rtol=1e-12, atol=1e-9)
+
+
+def test_a_thicker_column_under_a_pull_out_of_the_surface_draws_the_metal_toward_it():
+    """A flat layer thickening by 1 mm a metre under a pull of 100 m/s2 out of its surface: rho_l (g . n) grad H."""
+    col, pts = regular_grid(21)
+    H = 5e-4 + 1e-3 * pts[:, 0]
+    _, terms = skin.push(col, H, np.zeros(col.n), np.array([0.0, 0.0, 100.0]), V_HAT, RHO_L, SIGMA)
+    np.testing.assert_allclose(terms["thickness"], np.tile([RHO_L * 100.0 * 1e-3, 0.0, 0.0], (col.n, 1)), atol=1e-9)
+    assert np.abs(terms["deceleration"]).max() == 0.0                          # the pull lies along the normal
+
+
+def test_a_bump_s_surface_tension_pushes_the_metal_from_its_crest_converging_at_second_order():
+    """A cosine bump of 20 um on a 0.5 mm layer: kappa = -H'' for small slopes, so -sigma grad kappa = sigma H''' =
+    sigma eps k^3 sin(kx). Measured 2026-10-20: 2.45 % of its peak off at 20 spacings a wavelength, 0.62 % at 40."""
+    def error(lam):
+        col, pts = regular_grid()
+        k, eps = 2.0 * np.pi / lam, 2e-5
+        _, terms = skin.push(col, 5e-4 + eps * np.cos(k * pts[:, 0]), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
+        exact = SIGMA * eps * k ** 3 * np.sin(k * pts[:, 0])
+        inner = np.all((pts[:, :2] > 16e-3) & (pts[:, :2] < 64e-3), axis=1)
+        assert np.abs(terms["capillary"][inner, 1:]).max() < 1e-9 * np.abs(exact).max()
+        return np.abs(terms["capillary"][inner, 0] - exact[inner]).max() / np.abs(exact).max()
+    e20, e40 = error(20e-3), error(40e-3)
+    assert e40 < 0.01 and e20 / e40 == pytest.approx(4.0, rel=0.1)
+
+
+def test_a_sphere_s_own_mesh_drives_no_capillary_flow_to_speak_of(coarse_sphere_mesh):
+    """Uniform 0.8 mm skins on the 4 mm sphere: the curvature is uniform, so the capillary push is the mesh's error
+    alone -- 11 Pa/m at worst, 0.28 on average (measured 2026-10-20), against 2e5 Pa/m from the pressure at peak
+    heating (fact 153)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    _, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
+    size = np.linalg.norm(terms["capillary"], axis=1)
+    assert size.max() < 20.0 and size.mean() < 0.5
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_surface_flow.py
  tests/test_reentry_model_skin.py -q`. Expected: `AttributeError: 'SurfaceFlowResult' object has no attribute 'dp_ds'`
  and `AttributeError: module 'reentry_model.skin' has no attribute 'along_meridian'`.

- [ ] **Step 4: The pressure gradient on its own** in `reentry_model/surface_flow.py`:

```diff
--- a/reentry_model/surface_flow.py
+++ b/reentry_model/surface_flow.py
@@ -213,6 +213,8 @@
     tau_continuum: np.ndarray
     tau_fm: np.ndarray
     G: np.ndarray
+    dp_ds: np.ndarray = None      # the wall pressure's gradient along the meridian [Pa/m], from the fine bins (the
+                                  # pressure part of G; the skins' push takes it on its own, sub-plan 21 Task 5)
 
     @property
     def branch_name(self):
@@ -283,7 +285,7 @@
             return SurfaceFlowResult(BRANCH_FREE_MOLECULAR, state.kn, float("nan"), state.ma, float("nan"), p_inf,
                                      float("nan"), decel, self.gamma_pm, np.full(n, p_inf), zero, zero, zero, zero, zero,
                                      zero, nan, nan, np.full(n, np.inf), np.full(n, CLOSURE_COUETTE), np.ones(n, dtype=bool),
-                                     zero, zero, zero, zero)
+                                     zero, zero, zero, zero, zero)
         G_inertial = -rho_liquid * decel * sin
         if branch == BRANCH_FREE_MOLECULAR:
             # Schaaf-Chambre from the freestream directly: no compression model, no edge state, no boundary layer
@@ -292,7 +294,7 @@
                                      self.gamma_pm, np.full(n, p_inf), zero, zero, np.full(n, fs.rho), np.full(n, fs.T),
                                      zero, zero, nan, np.full(n, lam_inf), np.full(n, lam_inf / radius),
                                      np.full(n, CLOSURE_COUETTE), np.ones(n, dtype=bool), tau_fm, zero, tau_fm,
-                                     np.where(windward, G_inertial, 0.0))
+                                     np.where(windward, G_inertial, 0.0), zero)
         # shock-layer construction: the reservoir (also the bridge's continuum endpoint in the merged branch)
         stag = self.air.stagnation(fs.rho, fs.T, V)
         re_shock = float(fs.rho * V * radius / stag.mu)
@@ -328,9 +330,10 @@
         closure = np.where(windward, closure, CLOSURE_COUETTE)
         dp_ds = np.gradient(p_b, s_b)                              # from the actual p_w(theta), not a Newtonian formula
         G = np.where(windward, -at(dp_ds) + G_inertial, 0.0)
+        dp_ds_w = np.where(windward, at(dp_ds), 0.0)
         mach_e = np.where(a_e > 0.0, u_e / np.maximum(a_e, 1e-300), 0.0)
         u_e, rho_e, mu_e = (np.where(np.isfinite(x), x, 0.0) for x in (u_e, rho_e, mu_e))
         return SurfaceFlowResult(branch, float(state.kn), re_shock, float(state.ma), float(stag.p), p_inf, float(phi_star),
                                  decel, self.gamma_pm, p_w, u_e, u_e * slip, rho_e, T_e, mu_e, mach_e,
                                  np.where(windward, delta_a, np.nan), lam_w, np.where(windward, kn_local, np.inf),
-                                 closure, flagged | ~windward, tau, np.where(windward, tau_c, 0.0), tau_fm, G)
+                                 closure, flagged | ~windward, tau, np.where(windward, tau_c, 0.0), tau_fm, G, dp_ds_w)
```

- [ ] **Step 5: The push** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -930,6 +930,50 @@
     return book
 
 
+# -- the push (sub-plan 21 Task 5) -------------------------------------------------------------------------------------
+def along_meridian(normals, v_hat, magnitude):
+    """Vectors of the given magnitudes (n,) in each corner's tangent plane, pointing along the meridian away from the
+    stagnation point -- the direction of minus the flight direction's tangential part, the way the gas flows and the
+    surface flow's s runs; zero where that part vanishes (the stagnation and rear points)."""
+    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
+    away = -(v[None, :] - (normals @ v)[:, None] * normals)
+    size = np.linalg.norm(away, axis=1)
+    unit = np.divide(away, size[:, None], out=np.zeros_like(away), where=size[:, None] > 1e-12)
+    return np.asarray(magnitude, dtype=float)[:, None] * unit
+
+
+def push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma):
+    """The push along the surface on every layer of every column [N/m3], (n, 3), and its four terms (spec 6.4), each in
+    the column's own tangent plane -- the base's, where `momentum_step` solves for the velocities:
+
+        pressure      -dp_w/ds along the meridian away from the stagnation point, dp_ds [Pa/m, (n,)] the surface flow's
+                      (`SurfaceFlowResult.dp_ds`), differentiated on its own fine bins of the angle. The gradient of the
+                      corners' own p_w would carry the scatter of their normals' angles, amplified some 25 times (the
+                      radius over the spacing) by the differentiation: on the 2 mm sphere a Newtonian pressure's push
+                      came out 1.1 % off on average and 13 % at worst that way, from normals 0.06 deg off radial on
+                      average (measured 2026-10-20).
+        deceleration  rho_l g_t, g [m/s2, (3,) or (n, 3)] the acceleration the metal feels in the body's frame, minus the
+                      body's own -- for a body slowed by drag, the deceleration along the flight direction v_hat -- so
+                      that it pulls the metal forward, toward the nose along the surface and out of the front face;
+                      the spec's -rho_l a_t, and the element model's -rho_l W sin(theta) on the windward face.
+        thickness     rho_l (g . n) grad_s H, H [m, (n,)] the columns' thicknesses: a thicker column under a pull out
+                      of the surface has the lower pressure beneath it and draws metal in (the front face's
+                      Rayleigh-Taylor mechanism).
+        capillary     -sigma grad_s kappa, kappa the mean curvature (`Columns.curvature`) of the top surface the
+                      thicknesses make (`Columns.top_points`).
+
+    The two thickness terms are evaluated at the given H; Task 6 takes them implicitly in the column heights."""
+    H = np.asarray(thickness, dtype=float)
+    normals = columns.normal
+    g = np.broadcast_to(np.asarray(g, dtype=float), (columns.n, 3))
+    g_n = np.einsum("ij,ij->i", g, normals)
+    terms = {"pressure": -along_meridian(normals, v_hat, dp_ds),
+             "deceleration": rho_l * (g - g_n[:, None] * normals),
+             "thickness": rho_l * g_n[:, None] * columns.gradient(H),
+             "capillary": -sigma * columns.gradient(columns.curvature(columns.top_points(H)))}
+    return sum(terms.values()), terms
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

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

```diff
--- a/tests/test_reentry_model_columns.py
+++ b/tests/test_reentry_model_columns.py
@@ -130,3 +130,31 @@
 def test_the_curvature_of_a_plane_is_zero():
     pts, faces, inner = flat_grid()
     assert np.abs(columns.Columns(pts, faces).curvature()[inner]).max() < 1e-9
+
+
+# -- the dual faces as vectors (sub-plan 21 Task 6a) --------------------------------------------------------------------
+def test_the_dual_faces_close_every_column_so_a_uniform_flow_neither_gathers_nor_spreads_metal(coarse_sphere_mesh):
+    """`face_vector`: each dual face's two segments, each its normal in its triangle times its length, summed and
+    pointing from the edge's first corner to its second -- what a velocity crossing the face is dotted with. Around
+    every interior corner of a jittered plane the faces close to round-off, so a uniform flow carries as much into each
+    column as out of it; the edge's direction times the face's length does not close where the triangles are not
+    symmetric about their edges, and on this plane a uniform 1 m/s would gather or spread a column's metal at up to
+    539 1/s (measured 2026-10-21; 2.9e-13 with the faces' vectors). On the test sphere the faces close in each column's
+    tangent plane to 3.7 % of the mean face length at worst and 0.03 % on average; what is left is along the normal,
+    the surface's curvature (24 % of a face length on average)."""
+    pts, faces, inner = flat_grid()
+    col = columns.Columns(pts, faces)
+    i, j = col.edges[:, 0], col.edges[:, 1]
+    S = col.face_vector
+    assert np.all(np.einsum("ij,ij->i", S, pts[j] - pts[i]) > 0.0)                    # from the first corner out
+    assert np.all(np.linalg.norm(S, axis=1) <= col.face_length * (1.0 + 1e-12))
+    u = np.array([1.0, 0.3, 0.0])
+    assert np.abs((col.incidence @ (S @ u))[inner] / col.area[inner]).max() < 1e-9
+    e = pts[j] - pts[i]
+    naive = e / np.linalg.norm(e, axis=1)[:, None] * col.face_length[:, None]
+    assert np.abs((col.incidence @ (naive @ u))[inner] / col.area[inner]).max() > 100.0
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    closure = col.incidence @ col.face_vector
+    tangential = closure - np.einsum("ij,ij->i", closure, col.normal)[:, None] * col.normal
+    assert np.linalg.norm(tangential, axis=1).max() < 0.075 * col.face_length.mean()
+    assert np.linalg.norm(tangential, axis=1).mean() < 5e-4 * col.face_length.mean()
```

- [ ] **Step 2: The transport's tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -815,3 +815,296 @@
     _, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
     size = np.linalg.norm(terms["capillary"], axis=1)
     assert size.max() < 20.0 and size.mean() < 0.5
+
+
+# -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
+LIQUID, SLURRY, RIGID = 950.0, 902.0, 880.0       # K: all liquid; 30 % solid, flowing; rigid (T_rigid 895.1 K)
+
+
+def strip_grid(nx, ny, h=1e-3):
+    """Columns on a flat regular grid of nx by ny corners, h apart (regular_grid's triangles), normals +z."""
+    from reentry_model import columns
+    x, y = np.meshgrid(np.arange(nx) * h, np.arange(ny) * h, indexing="ij")
+    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
+    idx = np.arange(nx * ny).reshape(nx, ny)
+    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
+    return columns.Columns(pts, np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])), pts
+
+
+def stacked(col, temps):
+    """Skins on the columns `col`, column i holding a cell of 10 um of solid metal's mass at each temperature of
+    temps[i], top first."""
+    mat = scheil()
+    sk = skin.Skins(mat, col.area, skin.SkinSettings(**CELLS_04))
+    count = np.array([len(t) for t in temps])
+    K = int(count.max())
+    T = np.array([list(t) + [300.0] * (K - len(t)) for t in temps])
+    m = np.where(np.arange(K)[None, :] < count[:, None], mat.rho * col.area[:, None] * 1e-5, 0.0)
+    sk._install(m, m * mat.enthalpy(T), count)
+    return sk
+
+
+def moving_from(sk, columns, donor, u, cells=None):
+    """The velocity u (3,) in the donor column's cells (all, or those listed), every other cell at rest; returns the
+    donor's outflow rate across each edge per unit mass [1/s] (E,) toward each of its neighbours, and the neighbours."""
+    sk.U = np.zeros_like(sk.U)
+    cells = np.arange(sk.count[donor]) if cells is None else np.asarray(cells)
+    sk.U[donor, cells] = u
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    flux = np.where(i == donor, columns.face_vector @ u, np.where(j == donor, -(columns.face_vector @ u), 0.0))
+    other = np.where(i == donor, j, i)
+    out = flux > 0.0
+    return np.maximum(flux[out], 0.0) / sk.area[donor], other[out]
+
+
+def test_with_every_velocity_zero_the_transport_changes_nothing(coarse_sphere_mesh):
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    before = sk.state()
+    book = skin.transport(sk, col, 0.01)
+    assert all(np.array_equal(x, y) for x, y in zip(before, sk.state()))
+    assert not (book.mass.any() or book.enthalpy.any() or book.momentum.any() or book.wall.any()) and book.courant == 0.0
+
+
+def test_rigid_metal_joined_to_the_substrate_never_leaves_whatever_its_velocity():
+    """Anchored cells hold zero velocity after every momentum step (Task 4); given one anyway they do not move."""
+    col, _ = strip_grid(3, 3)
+    sk = stacked(col, [[RIGID] * 6] * col.n)
+    sk.U = np.where(sk.valid[..., None], np.array([0.3, -0.2, 0.0]), 0.0)
+    before = sk.state()
+    skin.transport(sk, col, 0.01)
+    assert all(np.array_equal(x, y) for x, y in zip(before, sk.state()))
+
+
+def test_a_warm_thick_pulse_carried_by_a_uniform_flow_arrives_where_the_flow_takes_it():
+    """A pulse half again as thick as the layer around it and 50 K warmer (a Gaussian 3 mm wide) in liquid on 1 mm
+    columns, every cell moving at 0.1 m/s along x for 0.1 s, against the same grid without the pulse (whose open ends
+    drain and pile up alike): its mass and heat are kept and their centroids move by u t, 10 mm (to 6e-8, 1e-5 and
+    8.3e-4 mm at Courant numbers nu = u dt / h of 0.5, 1 and 2, measured 2026-10-21, as the scheme's long upwind tail
+    reaches the grid's edges), while the implicit upwind scheme spreads them as the one-dimensional scheme does, by
+    h^2 nu (1 + nu) of variance a sub-step -- the numerical diffusion u h (1 + nu) / 2 -- less 4.2, 3.2 and 2.4 %. The
+    books close to round-off and every velocity stays u."""
+    h, u, x0, w = 1e-3, 0.1, 20e-3, 3e-3
+
+    def run(dt, pulse):
+        col, pts = strip_grid(61, 9, h)
+        g = np.exp(-0.5 * ((pts[:, 0] - x0) / w) ** 2) if pulse else np.zeros(col.n)
+        s = skin.SkinSettings(cell=20e-6)
+        sk = skin.Skins(scheil(), col.area, s)
+        sk.fill(skin.minimum_mass(sk.material, sk.area, s) * (1.0 + 0.5 * g), LIQUID + 50.0 * g)
+        sk.regrid()
+        sk.U = np.where(sk.valid[..., None], np.array([u, 0.0, 0.0]), 0.0)
+        M0, E0 = sk.mass(), sk.energy()
+        start = (sk.column_mass(), sk.H.sum(axis=1))
+        for _ in range(int(round(0.1 / dt))):
+            skin.transport(sk, col, dt)
+        assert sk.mass() == pytest.approx(M0, rel=1e-13) and sk.energy() == pytest.approx(E0, rel=1e-13)
+        assert np.abs(sk.U[sk.valid] - np.array([u, 0.0, 0.0])).max() < 1e-12
+        return start, (sk.column_mass(), sk.H.sum(axis=1)), pts
+
+    def moments(f, x):
+        mean = (x * f).sum() / f.sum()
+        return mean, ((x - mean) ** 2 * f).sum() / f.sum(), f.sum()
+
+    for dt in (0.005, 0.01, 0.02):
+        (m0, H0), (m1, H1), pts = run(dt, True)
+        (r0, R0), (r1, R1), _ = run(dt, False)
+        inner = (pts[:, 1] > 0.0) & (pts[:, 1] < pts[:, 1].max())
+        nu = u * dt / h
+        for before, after in (((m0 - r0), (m1 - r1)), ((H0 - R0), (H1 - R1))):
+            a, b = moments(before[inner], pts[inner, 0]), moments(after[inner], pts[inner, 0])
+            assert b[2] == pytest.approx(a[2], rel=1.5e-5)                          # 6.2e-6 at nu = 2
+            assert b[0] - a[0] == pytest.approx(u * 0.1, abs=2e-6)
+            assert b[1] - a[1] == pytest.approx(h * h * nu * (1.0 + nu) * 0.1 / dt, rel=0.09)
+
+
+def test_depth_is_matched_by_fraction_of_the_flowing_layer():
+    """Spec 6.5: metal leaving the donor between fractions s1 and s2 of its flowing depth joins each receiver between
+    the same fractions of the receiver's, which flows from its top down. The donor's ten liquid cells over rigid metal
+    move only in their first and third (fractions 0-0.1 and 0.2-0.3); a receiver of 20 liquid cells takes them in its
+    cells 0-1 and 4-5, a receiver of 4 in its cells 0 and 1 (0.2-0.3 straddles its 0.25), each in proportion to the
+    overlap; nothing reaches rigid metal."""
+    col, _ = strip_grid(3, 3)
+    temps = [[LIQUID] * 20 + [RIGID] * 5 if k % 2 else [LIQUID] * 4 + [RIGID] * 10 for k in range(col.n)]
+    temps[4] = [LIQUID] * 10 + [RIGID] * 5
+    sk = stacked(col, temps)
+    dt, u = 0.01, np.array([0.08, 0.05, 0.0])
+    rate, receivers = moving_from(sk, col, 4, u, cells=[0, 2])
+    before = sk.m.copy()
+    book = skin.transport(sk, col, dt, regrid=False)
+    m_left = before[4, [0, 2]] / (1.0 + dt * rate.sum())               # the donor's cells, implicitly
+    for r, c in zip(receivers, rate):
+        gain = sk.m[r] - before[r]
+        sent = dt * c * m_left                                        # its first and third cells' share
+        if sk.count[r] == 25:
+            np.testing.assert_allclose(gain[[0, 1, 4, 5]], [sent[0] / 2, sent[0] / 2, sent[1] / 2, sent[1] / 2],
+                                       rtol=1e-12)
+            assert np.all(gain[2:4] == 0.0) and np.all(gain[6:] == 0.0)
+        else:
+            np.testing.assert_allclose(gain[:2], [sent[0] + sent[1] / 2, sent[1] / 2], rtol=1e-12)
+            assert np.all(gain[2:] == 0.0)
+        assert book.mass[r] == pytest.approx(sent.sum(), rel=1e-12)
+    assert len(receivers) >= 2 and {25, 14} <= set(sk.count[receivers])
+
+
+def test_metal_arriving_on_a_rigid_top_lands_on_it_and_the_receiver_s_own_cells_are_untouched():
+    """Spec 6.5's arrival rule: where the receiver's top is rigid -- anchored solid, or a crust over flowing metal --
+    the arriving metal lands on top of its top cell, as a new cell holding exactly what arrived; nothing passes
+    through rigid metal."""
+    col, _ = strip_grid(3, 3)
+    temps = [[RIGID] * 12 if k % 2 else [RIGID] * 2 + [LIQUID] * 6 + [RIGID] * 4 for k in range(col.n)]
+    temps[4] = [LIQUID] * 8 + [RIGID] * 4
+    sk = stacked(col, temps)
+    rate, receivers = moving_from(sk, col, 4, np.array([0.1, 0.04, 0.0]), cells=np.arange(8))
+    before = sk.state()
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    h_donor = sk.material.enthalpy(LIQUID)
+    for r in range(col.n):
+        if r == 4:
+            continue
+        if r in receivers:
+            assert sk.count[r] == before[4][r] + 1
+            assert sk.m[r, 0] == pytest.approx(book.mass[r], rel=1e-12)
+            assert sk.H[r, 0] == pytest.approx(book.mass[r] * h_donor, rel=1e-12)
+            np.testing.assert_allclose(sk.U[r, 0], [0.1, 0.04, 0.0], rtol=1e-12)
+            shifted = slice(1, sk.count[r])
+        else:
+            assert book.mass[r] == 0.0
+            shifted = slice(0, sk.count[r])
+        for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3]), (sk.U, before[5])):
+            assert np.array_equal(now[r, shifted], then[r, :before[4][r]])
+    assert any(r % 2 for r in receivers) and any(r % 2 == 0 for r in receivers)          # solid tops and crusts
+
+
+def test_a_floating_crust_carried_across_lands_on_top_and_the_slurry_beneath_it_joins_by_fraction():
+    """Asha, 2026-10-20 (Task 4): a crust over flowing metal floats and moves as one piece. Carried into a neighbour it
+    lands on top of the neighbour's top; the slurry beneath it, ten cells, joins each receiver's 20 liquid cells two by
+    two; the donor's anchored cells, given the same velocity, stay where they are."""
+    col, _ = strip_grid(3, 3)
+    temps = [[LIQUID] * 20 + [RIGID] * 5] * col.n
+    temps[4] = [RIGID] * 2 + [SLURRY] * 10 + [RIGID] * 5
+    sk = stacked(col, temps)
+    dt = 0.01
+    rate, receivers = moving_from(sk, col, 4, np.array([0.1, 0.04, 0.0]))       # every cell, the anchored too
+    before = sk.state()
+    book = skin.transport(sk, col, dt, regrid=False)
+    left = before[1][4, :12] / (1.0 + dt * rate.sum())
+    np.testing.assert_allclose(sk.m[4, :12], left, rtol=1e-12)
+    for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3]), (sk.U, before[5])):
+        assert np.array_equal(now[4, 12:17], then[4, 12:17])                    # the anchored cells
+    h_crust, h_slurry = sk.material.enthalpy(RIGID), sk.material.enthalpy(SLURRY)
+    for r, c in zip(receivers, rate):
+        assert sk.count[r] == 26
+        assert sk.m[r, 0] == pytest.approx(dt * c * left[:2].sum(), rel=1e-12)
+        assert sk.H[r, 0] == pytest.approx(sk.m[r, 0] * h_crust, rel=1e-12)
+        gain = sk.m[r, 1:21] - before[1][r, :20]
+        np.testing.assert_allclose(gain, np.repeat(dt * c * left[2:12] / 2.0, 2), rtol=1e-12)
+        H_gain = sk.H[r, 1:21] - before[2][r, :20]
+        np.testing.assert_allclose(H_gain, gain * h_slurry, rtol=1e-9)
+        assert np.array_equal(sk.m[r, 21:], before[1][r, 20:25])
+        assert book.mass[r] == pytest.approx(dt * c * left.sum(), rel=1e-12)
+
+
+def mixed_sphere_skins(coarse_sphere_mesh, seed=6, speed=0.05):
+    """Skins on the 4 mm test sphere, 1-3 times their minimum deep, every cell at its own temperature between 860 and
+    960 K -- liquid, slurry, crusts and anchored rigid metal in every arrangement -- and moving at its own velocity in
+    its column's tangent plane, about `speed` m/s."""
+    sk, col, rng = random_sphere_skins(coarse_sphere_mesh, seed)
+    sk.T = np.where(sk.valid, rng.uniform(860.0, 960.0, sk.T.shape), sk.T)
+    sk.H = sk.m * sk.material.enthalpy(sk.T)
+    U = rng.normal(0.0, speed, sk.U.shape)
+    U -= np.einsum("ikx,ix->ik", U, col.normal)[..., None] * col.normal[:, None, :]
+    sk.U = np.where(sk.valid[..., None], U, 0.0)
+    return sk, col, rng
+
+
+def test_the_books_close_on_skins_of_every_kind_on_the_sphere(coarse_sphere_mesh):
+    """Mass, enthalpy and the cells' carried sensitivity are conserved to round-off, each column's change is its book;
+    the momentum of each column changes by what was carried in less what was carried out plus the wall's impulse, the
+    normal parts the projection removed, which lie along the column's normal; every velocity ends in its column's
+    plane, and the anchored cells are untouched. Measured 2026-10-21 at Courant numbers up to 0.87: every book to
+    4e-15 of its scale, the wall's impulse 2.5 % of the largest column's momentum."""
+    sk, col, rng = mixed_sphere_skins(coarse_sphere_mesh)
+    sk.S = np.where(sk.valid, rng.uniform(0.0, 1.0, sk.m.shape), 0.0)
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    anchored = skin.anchored_cells(sk, flowing)
+    assert flowing.sum() > 1000 and anchored.sum() > 1000 and (sk.valid & ~flowing & ~anchored).sum() > 100
+    before = sk.state()
+    mS0, P0, M0, E0 = (sk.m * sk.S).sum(), sk.column_momentum(), sk.column_mass(), sk.H.sum(axis=1)
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    assert book.courant > 0.1
+    np.testing.assert_allclose(sk.column_mass() - M0, book.mass, rtol=0, atol=1e-13 * M0.max())
+    np.testing.assert_allclose(sk.H.sum(axis=1) - E0, book.enthalpy, rtol=0, atol=1e-13 * np.abs(E0).max())
+    assert abs(book.mass.sum()) < 1e-13 * M0.sum() and abs(book.enthalpy.sum()) < 1e-13 * np.abs(E0).sum()
+    assert (sk.m * sk.S).sum() == pytest.approx(mS0, rel=1e-13)
+    scale = np.abs(P0).max()
+    np.testing.assert_allclose(sk.column_momentum() - P0, book.momentum + book.wall, rtol=0, atol=1e-13 * scale)
+    assert np.abs(book.momentum.sum(axis=0)).max() < 1e-13 * scale
+    assert np.linalg.norm(np.cross(book.wall, col.normal), axis=1).max() < 1e-13 * scale
+    assert np.abs(book.wall).max() > 1e-6 * scale                                  # the planes differ: it is not idle
+    assert np.abs(np.einsum("ikx,ix->ik", sk.U, col.normal)).max() < 1e-12 * np.abs(sk.U).max()
+    top = (sk.count - before[4])[:, None]                                         # a new top cell shifts the rest
+    k = np.arange(anchored.shape[1])[None, :]
+    rows, cells = np.nonzero(anchored)
+    for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3])):
+        assert np.array_equal(now[rows, cells + top[rows, 0]], then[rows, cells])
+
+
+def test_courant_numbers_far_above_one_stay_positive_bounded_and_conservative(coarse_sphere_mesh):
+    """A 10 s sub-step on the mixed sphere skins (Courant numbers up to 875): every cell keeps a positive mass, every
+    specific enthalpy stays within the old ones' range and every speed below the old fastest (the scheme mixes, it
+    never overshoots), the books close and the cells end within their split and merge bounds."""
+    sk, col, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    v = sk.valid
+    h0, M0, E0 = sk.H[v] / sk.m[v], sk.mass(), sk.energy()
+    fastest = np.linalg.norm(sk.U, axis=2).max()
+    book = skin.transport(sk, col, 10.0)
+    v = sk.valid
+    assert book.courant > 100.0
+    assert np.all(sk.m[v] > 0.0)
+    h = sk.H[v] / sk.m[v]
+    assert h.min() >= h0.min() - 1e-9 * abs(h0.min()) and h.max() <= h0.max() + 1e-9 * abs(h0.max())
+    assert np.linalg.norm(sk.U, axis=2).max() <= fastest * (1.0 + 1e-12)
+    assert sk.mass() == pytest.approx(M0, rel=1e-13) and sk.energy() == pytest.approx(E0, rel=1e-13)
+    dz, t = sk.dz(), sk.targets()
+    assert np.all(dz[v] >= skin.MERGE_BELOW * t[v]) and np.all(dz[v] <= skin.SPLIT_ABOVE * t[v])
+
+
+def test_a_velocity_carried_into_a_column_is_projected_onto_its_plane(coarse_sphere_mesh):
+    """One column on the 4 mm test sphere, liquid throughout, moving at 1 m/s in its plane into its neighbours, which are
+    liquid at rest: each neighbour gains the projection onto its own plane of the momentum that left toward it, and
+    the wall takes the normal part -- the planes of neighbouring 4 mm columns on a 50 mm sphere are 4.6 degrees apart
+    on average, the three receivers' 4.5-4.9 degrees from the donor's, and the normal parts 2.3-7.7 % of the momentum
+    carried (measured 2026-10-21)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    sk = skin.Skins(scheil(), col.area, skin.SkinSettings())
+    sk.fill(skin.minimum_mass(sk.material, sk.area, sk.settings), LIQUID)
+    sk.regrid()
+    donor = int(np.argmin(np.abs(col.normal[:, 2])))                              # on the equator of z
+    u = np.cross(col.normal[donor], [0.0, 0.0, 1.0])
+    rate, receivers = moving_from(sk, col, donor, u / np.linalg.norm(u))
+    P0 = sk.column_momentum()
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    n = col.normal[receivers]
+    carried = book.momentum[receivers]
+    inplane = carried - np.einsum("ij,ij->i", carried, n)[:, None] * n
+    np.testing.assert_allclose(sk.column_momentum()[receivers] - P0[receivers], inplane, rtol=0,
+                               atol=1e-12 * np.abs(carried).max())
+    np.testing.assert_allclose(book.wall[receivers], inplane - carried, rtol=0, atol=1e-12 * np.abs(carried).max())
+    normal_part = np.linalg.norm(carried - inplane, axis=1) / np.linalg.norm(carried, axis=1)
+    assert normal_part.max() > 0.02 and normal_part.max() < 0.16
+
+
+def test_the_upwind_order_changes_the_cost_never_the_answer(coarse_sphere_mesh, monkeypatch):
+    """The solve takes the cells in upwind order -- scipy's strong components, which come out in reverse topological
+    order -- so that its sparse LU has almost no fill; were that order ever not upwind, the LU falls back to its own
+    fill-reducing order and gives the same answer, more slowly."""
+    first, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    second, col, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    a = skin.transport(first, col, 0.05, regrid=False)
+    rng = np.random.default_rng(0)
+    monkeypatch.setattr(skin, "_upwind_order", lambda src, dst, n: (rng.permutation(n), False))
+    b = skin.transport(second, col, 0.05, regrid=False)
+    for x, y in zip(first.state(), second.state()):
+        np.testing.assert_allclose(x, y, rtol=1e-11, atol=1e-11 * np.abs(y).max())
+    np.testing.assert_allclose(a.wall, b.wall, rtol=0, atol=1e-11 * np.abs(a.momentum).max())
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 11 failed — `AttributeError: 'Columns' object has no attribute
  'face_vector'` (5) and `AttributeError: module 'reentry_model.skin' has no attribute 'transport'` (6).

- [ ] **Step 4: The faces' vectors** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -23,7 +23,14 @@
 and 23 %). On these 50 mm-radius spheres the curvature is then within 8.4·10^-5 of 2/R at every corner (1.3·10^-5 on
 the 2 mm sphere), where the barycentric area had scattered 7.6 and 3.9 %, and the gradient of a linear field within
 0.09 % of its tangential part on average (0.02 % on the 2 mm sphere), 0.8 % at the test sphere's one irregular corner
-(measured 2026-10-18)."""
+(measured 2026-10-18).
+
+The flow's sideways transport (sub-plan 21 Task 6a) crosses each dual face at its `face_vector`: the face's two
+segments, each its normal in its own triangle times its length, summed and pointing from the edge's first corner to
+its second, so that a velocity u carries u . face_vector per unit depth across the face. A column's faces close
+exactly on a plane, so a uniform flow carries as much into every column as out of it; the edge's direction times the
+face's length closes only where the triangles are symmetric about their edges, and on the tests' jittered plane of 1 mm
+squares a uniform 1 m/s would gather or spread a column's metal at up to 539 1/s that way (measured 2026-10-21)."""
 from dataclasses import dataclass
 
 import numpy as np
@@ -69,15 +76,21 @@
         mid = 0.5 * (p[i] + p[j])
         cot = np.zeros(len(self.edges))
         face = np.zeros(len(self.edges))
+        vector = np.zeros((len(self.edges), 3))
+        _, unit, _ = triangle_geometry(p, self.faces)
         for s in range(2):
             ok = self._tri[:, s] >= 0
             k = self._opp[ok, s]
             u, v = p[i[ok]] - p[k], p[j[ok]] - p[k]
             cot[ok] += np.einsum("ij,ij->i", u, v) / np.linalg.norm(np.cross(u, v), axis=1)
-            face[ok] += np.linalg.norm(mid[ok] - self.centroid[self._tri[ok, s]], axis=1)
+            segment = self.centroid[self._tri[ok, s]] - mid[ok]
+            face[ok] += np.linalg.norm(segment, axis=1)
+            across = np.cross(segment, unit[self._tri[ok, s]])             # the segment's normal times its length
+            vector[ok] += across * np.sign(np.einsum("ij,ij->i", across, p[j[ok]] - p[i[ok]]))[:, None]
         self.cot_weight = 0.5 * cot
         self.conduct_weight = np.maximum(self.cot_weight, 0.0)
         self.face_length = face
+        self.face_vector = vector
 
     def gradient(self, f, points=None):
         """The tangential gradient at every corner of a per-corner field f (n,) [f per m] on the surface through
```

- [ ] **Step 5: The transport** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -31,6 +31,9 @@
 from dataclasses import dataclass
 
 import numpy as np
+import scipy.sparse as sp
+import scipy.sparse.linalg as spla
+from scipy.sparse.csgraph import connected_components
 
 from . import rheology
 from .thermal import SIGMA_SB
@@ -974,6 +977,227 @@
     return sum(terms.values()), terms
 
 
+# -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
+SLIVER = 1.0e-12                  # a piece of a donor cell shorter than this share of it is the round-off of a boundary
+                                  # the donor's and the receiver's partitions share, and is left out
+
+
+@dataclass
+class TransportBook:
+    """What a transport sub-step carried between the columns, per column: each column's gain from its neighbours less
+    what it gave them. The mass, the enthalpy and the momentum as it left the donors each sum to zero over the body to
+    round-off; every arrival's velocity is projected onto its receiver's plane, and the normal parts removed are the
+    wall's impulse, so that each column's momentum changes by momentum + wall."""
+    mass: np.ndarray         # kg (n,)
+    enthalpy: np.ndarray     # J (n,)
+    momentum: np.ndarray     # kg m/s (n, 3), as it left the donors, before the projection
+    wall: np.ndarray         # N s (n, 3), the normal parts removed, each along its own column's normal
+    courant: float           # the largest outflow rate of any cell [1/s] times the sub-step
+
+
+def tangent_frames(normals):
+    """Two orthonormal vectors (e1, e2), each (n, 3), spanning every column's tangent plane, e1 x e2 along the normal:
+    e1 the global x axis's part in the plane, or y's where the normal lies within 26 degrees of x."""
+    ref = np.where((np.abs(normals[:, 0]) < 0.9)[:, None], np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
+    e1 = ref - np.einsum("ij,ij->i", ref, normals)[:, None] * normals
+    e1 /= np.linalg.norm(e1, axis=1)[:, None]
+    return e1, np.cross(normals, e1)
+
+
+def _fractions(depth):
+    """Per column, its cells' boundaries (n, K+1) as fractions of the depth they hold (`depth` (n, K) each cell's share,
+    zero for a cell that does not count), from 0 at the top to 1 at the last that counts and held there; 1 throughout
+    where none counts."""
+    z = np.concatenate([np.zeros((depth.shape[0], 1)), np.cumsum(depth, axis=1)], axis=1)
+    total = z[:, -1:]
+    return np.where(total > 0.0, z / np.where(total > 0.0, total, 1.0), 1.0)
+
+
+def _pieces(sig, tau):
+    """The pieces two partitions of [0, 1] cut each other into, row by row: sig and tau (D, K+1) non-decreasing
+    boundaries. Returns, for every piece of positive length, its row, its cell of sig, its cell of tau and its length.
+    A rank merge: each boundary's place in the merged row is its own index plus the number of the other partition's
+    boundaries before it (tau's first where two are equal)."""
+    D, B = sig.shape
+    rows = np.arange(D)[:, None]
+    off = (np.arange(D) * 4.0)[:, None]                      # boundaries lie in [0, 1]: rows 4 apart never mix
+    before = lambda xp, x, side: (np.searchsorted((xp + off).ravel(), (x + off).ravel(), side=side).reshape(x.shape)
+                                  - rows * B)
+    at_s = np.arange(B)[None, :] + before(tau, sig, "right")
+    at_t = np.arange(B)[None, :] + before(sig, tau, "left")
+    merged = np.empty((D, 2 * B))
+    merged[rows, at_s], merged[rows, at_t] = sig, tau
+    from_sig = np.zeros((D, 2 * B), dtype=np.int64)
+    from_sig[rows, at_s] = 1
+    seen = np.cumsum(from_sig, axis=1)[:, :-1]                # sig's boundaries at or before each piece's start
+    length = np.diff(merged, axis=1)
+    r, p = np.nonzero(length > 0.0)
+    return r, seen[r, p] - 1, p - seen[r, p], length[r, p]
+
+
+def _upwind_order(src, dst, n):
+    """An order of n unknowns in which every donor (src) comes before its receivers (dst), except within a cycle of
+    transfers, which stays together; and whether it holds. The order is scipy's strong components, which its (Pearce's)
+    algorithm numbers in reverse topological order -- how scipy numbers them, not what it promises, so it is checked."""
+    graph = sp.csr_matrix((np.ones(len(src)), (src, dst)), shape=(n, n))
+    _, label = connected_components(graph, directed=True, connection="strong")
+    return np.argsort(-label, kind="stable"), bool(np.all(label[src] >= label[dst]))
+
+
+def _factor(A, upwind):
+    """A sparse LU of A, an M-matrix or one by blocks, in the given order with no pivoting where the order is upwind (the
+    matrix is then triangular but for its cycles: no fill, every sign kept), else in SuperLU's own fill-reducing order."""
+    if upwind:
+        return spla.splu(A, permc_spec="NATURAL", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
+    return spla.splu(A)
+
+
+def transport(sk, columns, dt, regrid=True):
+    """One linearly implicit upwind sub-step dt [s] of the metal moving sideways between neighbouring columns at the
+    cells' velocities sk.U (spec 6.5; Asha, 2026-10-21: its own step, after the momentum step). Updates the cells and
+    returns a TransportBook.
+
+    Across each face, cell k of the donor column a sends c m_k' toward the receiver b, c = (u_k . f)+ / A_a [1/s], f the
+    face's vector (`Columns.face_vector`) and m_k' the cell's mass at the sub-step's end; it carries the cell's specific
+    enthalpy and velocity. Every cell is then m_k' (1 + dt sum c) = m_k + dt (what arrives): (I + dt C) m' = m, the
+    film's runoff scheme (`film.Runoff`) on the cells -- an M-matrix whose columns sum to one, so the masses stay
+    positive, the books close to round-off and nothing limits dt. The enthalpy, the carried sensitivity S and the
+    momentum solve the same system, the momentum with every arrival projected onto its receiver's tangent plane, in two
+    components per cell (`tangent_frames`), the normal part booked to the wall (spec 6.4, section 14 assumption 2).
+
+    Where metal joins the receiver (spec 6.5, section 14 assumption 3): a flowing donor cell's span of fractions of the
+    donor's flowing depth (its flowing cells, top down) joins the receiver's cells between the same fractions of its
+    flowing layer, the run of flowing cells from its top, each in proportion to the overlap. Where the receiver's top is
+    rigid -- anchored metal, or a crust -- everything lands on top of it, as a new top cell; so does a floating crust's
+    metal (Task 4, Asha 2026-10-20) wherever it goes, being no part of any flowing depth. Rigid metal joined to the
+    substrate (`anchored_cells`) neither leaves nor receives, whatever velocity it holds: nothing passes through rigid
+    metal. Only the cells that send or receive are solved for; with every velocity zero nothing changes.
+
+    The cells are taken in upwind order (`_upwind_order`), in which the system is triangular but for its cycles where
+    flows meet, so its sparse LU has almost no fill. With `regrid` the cells are then split and merged (spec 6.2); a
+    sub-step that regrids anyway passes False.
+
+    Measured 2026-10-21 on the 100 mm sphere's skins (9 417 columns, 10 um cells), flowing along the meridian away from
+    the stagnation point: 0.27 s a sub-step with the front hemisphere's top 0.3 mm flowing (124 000 of 553 000 cells
+    moving), the same at Courant numbers of 1 and 96, and 0.95 s with all 584 000 cells moving, against 0.6-0.8 s for a
+    coupled conduction sub-step on the sphere's skins at 700 K. With every cell moving, the mass's LU in the cells' own
+    order took 19.6 s and in SuperLU's fill-reducing order 0.66 s, against 0.06 s upwind. Projected once at the
+    sub-step's end instead of at every arrival, a
+    uniform meridian flow kept 90 % of its kinetic energy in one sub-step at Courant number 9.6 against 96 %, and 24 %
+    against 49 % at 96: metal crossing many columns in a sub-step would carry its velocity straight on, off the curved
+    surface."""
+    n, K = sk.m.shape
+    valid = sk.valid
+    flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
+    moving = valid & ~anchored_cells(sk, flowing) & (sk.U != 0.0).any(axis=2)
+    book = TransportBook(np.zeros(n), np.zeros(n), np.zeros((n, 3)), np.zeros((n, 3)), 0.0)
+    if not moving.any():
+        return book
+    # every face a moving column has, each way, and every cell's outflow rate across it [1/s]
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    a, b, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
+    keep = np.flatnonzero(moving.any(axis=1)[a])
+    a, b, f = a[keep], b[keep], f[keep]
+    c = np.where(moving[a], np.maximum(np.einsum("dkx,dx->dk", sk.U[a], f), 0.0), 0.0) / sk.area[a][:, None]
+    keep = np.flatnonzero((c > 0.0).any(axis=1))
+    if not len(keep):
+        return book
+    a, b, c = a[keep], b[keep], c[keep]
+    D = len(a)
+    # where it joins: by fraction of the flowing depth, or on top
+    dz = sk.dz()
+    run = np.cumprod(flowing, axis=1).sum(axis=1)             # the receiver's flowing layer: its top run of flowing cells
+    sig = _fractions(np.where(flowing, dz, 0.0))
+    tau = _fractions(np.where(np.arange(K)[None, :] < run[:, None], dz, 0.0))
+    joins = np.where(flowing[a] & (run[b] > 0)[:, None], c, 0.0)
+    rows = np.flatnonzero((joins > 0.0).any(axis=1))
+    r, kd, kr, length = _pieces(sig[a[rows]], tau[b[rows]])
+    some = joins[rows[r], kd] > 0.0
+    edge, kd, kr, length = rows[r[some]], kd[some], kr[some], length[some]
+    group = edge * K + kd
+    span = np.bincount(group, length, D * K)[group]
+    some = length > SLIVER * span                             # not the round-off of a boundary the two share
+    edge, kd, kr, length, group = edge[some], kd[some], kr[some], length[some], group[some]
+    share = length / np.bincount(group, length, D * K)[group]  # each donor cell's pieces sum to one
+    te, tk = np.nonzero((c > 0.0) & ~(flowing[a] & (run[b] > 0)[:, None]))
+    # every arrival: its face, the receiving cell (n K + b: a new top cell of column b), the sending cell, dt c share
+    e_edge = np.r_[edge, te]
+    e_row = np.r_[b[edge] * K + kr, n * K + b[te]]
+    e_col = np.r_[a[edge] * K + kd, a[te] * K + tk]
+    e_val = dt * np.r_[c[edge, kd] * share, c[te, tk]]
+    # the cells that send or receive, in upwind order; each sends what its arrivals take, so the books close exactly
+    active = np.zeros(n * K + n, dtype=bool)
+    active[e_row], active[e_col] = True, True
+    act = np.flatnonzero(active)
+    N = len(act)
+    local = np.full(n * K + n, -1, dtype=np.int64)
+    local[act] = np.arange(N)
+    order, upwind = _upwind_order(local[e_col], local[e_row], N)
+    at = np.empty(N, dtype=np.int64)
+    at[order] = np.arange(N)
+    to, frm = at[local[e_row]], at[local[e_col]]
+    diag = 1.0 + np.bincount(local[e_col], e_val, N)
+    # mass, enthalpy and the carried sensitivity
+    lu = _factor(sp.csc_matrix((np.r_[diag, -e_val], (np.r_[at, to], np.r_[at, frm])), shape=(N, N)), upwind)
+    cell = lambda v: np.r_[v.reshape((n * K,) + v.shape[2:]), np.zeros((n,) + v.shape[2:])][act]
+    rhs = [cell(sk.m), cell(sk.H)] + ([cell(sk.m * sk.S)] if sk.S is not None else [])
+    x = np.empty((N, len(rhs)))
+    x[order] = lu.solve(np.column_stack(rhs)[order])
+    # the momentum, each arrival projected onto its receiver's plane: two components per cell in its own plane's frame
+    owner = np.where(act < n * K, act // K, act - n * K)
+    e1, e2 = tangent_frames(columns.normal)
+    P = cell(sk.m[..., None] * sk.U)
+    p = np.column_stack([np.einsum("ij,ij->i", P, e1[owner]), np.einsum("ij,ij->i", P, e2[owner])])
+    held = P - p[:, :1] * e1[owner] - p[:, 1:] * e2[owner]   # any normal part a cell held: the wall's
+    R = np.stack([np.stack([np.einsum("ij,ij->i", e[b], g[a]) for g in (e1, e2)], axis=1) for e in (e1, e2)], axis=1)
+    Re = R[e_edge]                                            # [arrival, receiver's component, donor's component]
+    rows2 = [2 * at, 2 * at + 1] + [2 * to + q for q in (0, 0, 1, 1)]
+    cols2 = [2 * at, 2 * at + 1] + [2 * frm + q for q in (0, 1, 0, 1)]
+    vals2 = [diag, diag] + [-e_val * Re[:, q, s] for q, s in ((0, 0), (0, 1), (1, 0), (1, 1))]
+    lu2 = _factor(sp.csc_matrix((np.concatenate(vals2), (np.concatenate(rows2), np.concatenate(cols2))),
+                                shape=(2 * N, 2 * N)), upwind)
+    rhs2 = np.empty(2 * N)
+    rhs2[2 * at], rhs2[2 * at + 1] = p[:, 0], p[:, 1]
+    y = lu2.solve(rhs2)
+    p = np.column_stack([y[2 * at], y[2 * at + 1]])
+    P = p[:, :1] * e1[owner] + p[:, 1:] * e2[owner]
+    # the books, face by face: what every arrival took from its sending cell at the sub-step's end
+    moved = e_val[:, None] * np.column_stack([x[:, :2], p])[local[e_col]]
+    sent = np.stack([np.bincount(e_edge, moved[:, v], D) for v in range(4)], axis=1)
+    carried = sent[:, 2:3] * e1[a] + sent[:, 3:4] * e2[a]
+    lost = np.einsum("ij,ij->i", carried, columns.normal[b])[:, None] * columns.normal[b]
+    net = lambda v: np.bincount(b, v, n) - np.bincount(a, v, n)
+    book = TransportBook(net(sent[:, 0]), net(sent[:, 1]), np.stack([net(carried[:, v]) for v in range(3)], axis=1),
+                         -np.stack([np.bincount(b, lost[:, v], n) + np.bincount(owner, held[:, v], n)
+                                    for v in range(3)], axis=1), float((diag - 1.0).max()))
+    # install: new values for the cells solved for, a new top cell where metal landed on top
+    m, H, T, U = sk.m.copy(), sk.H.copy(), sk.T.copy(), sk.U.copy()
+    S = sk.S.copy() if sk.S is not None else None
+    old = act < n * K
+    m.ravel()[act[old]], H.ravel()[act[old]], T.ravel()[act[old]] = x[old, 0], x[old, 1], np.nan
+    U.reshape(-1, 3)[act[old]] = P[old] / x[old, :1]
+    if S is not None:
+        S.ravel()[act[old]] = x[old, 2] / x[old, 0]
+    count = sk.count.copy()
+    tops = act[~old] - n * K
+    if len(tops):
+        new = np.zeros(n, dtype=bool)
+        new[tops] = True
+        src = np.clip(np.arange(K + 1)[None, :] - new[:, None], 0, K)
+        down = lambda v: np.take_along_axis(np.concatenate([v, np.zeros((n, 1) + v.shape[2:])], axis=1),
+                                            src.reshape(src.shape + (1,) * (v.ndim - 2)), axis=1)
+        m, H, T, U = down(m), down(H), down(T), down(U)
+        m[tops, 0], H[tops, 0], T[tops, 0], U[tops, 0] = x[~old, 0], x[~old, 1], np.nan, P[~old] / x[~old, :1]
+        if S is not None:
+            S = down(S)
+            S[tops, 0] = x[~old, 2] / x[~old, 0]
+        count = count + new
+    sk._install(m, H, count, S, T, U)
+    if regrid:
+        sk.regrid()
+    return book
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_columns.py` 10 passed and
  `tests/test_reentry_model_skin.py` 56 passed; the unit tier 450 passed, 1 skipped. Nothing calls `transport` yet: the flow
  joins the skin body's sub-steps in 6c, so no result changes.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/columns.py tests/test_reentry_model_columns.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the metal moves sideways between neighbouring columns at the cells' velocities, linearly implicit and upwind, by fraction of the flowing depth and on top of a rigid top, each arrival's velocity projected onto its receiver's plane (sub-plan 21 Task 6a)"
```

---

## Tasks still in outline (to be written in full)

6. **Sideways exchange**, split in three (Asha, 2026-10-21), each with its own tests and commit; 6a, the transport at
   given velocities, is written out above:
   - **6b. The thickness terms implicit**: ρ_l (g·n) ∇_s H and −σ ∇_s κ (Task 5) linearised in the column heights
     inside 6a's transport (spec §6.4's last paragraph). Tests: a bump on a flat layer under deceleration grows at the
     thin-film Rayleigh–Taylor linear rate; bumps shorter than the surface-tension cutoff decay (from Task 5, Asha,
     2026-10-20).
   - **6c. The flow in the skin body**: the surface flow evaluated on the corner surface (its `dp_ds` and `tau`),
     `push` and `along_meridian(tau)` driving `momentum_step`, then the transport, in every sub-step of `skin.advance`;
     a setting to turn the flow off, which reproduces the heat-only skins; the momentum books in the run's history.
     The flow's sub-step is chosen here (option c, Asha, 2026-10-21): measure the Courant numbers and the momentum lost
     to 6a's smearing on the 100 mm sphere's melting state. Asha expects the flow to need a sub-step shorter than the
     heat's 10 ms (about 100 columns a sub-step are estimated in flight, fact 154), so the momentum step and the
     transport should be able to sub-cycle within a heat sub-step.
7. **The droplet momentum sink** placeholder at the top cells for stage 4, tested with a prescribed release.
8. **Verification**: level 2's plate (spec §11.1, amended 2026-10-12: the box's 90-degree edges fold the surface; sub-plan
   20's `analysis/skin_plate.py`) with flow on and a push along the face (mass, momentum and energy exact; sensible
   behaviour); the momentum-balance residual and the roughness history columns; the cost of the sideways solve.

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
