# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level so far: Task 1 (2026-10-18).** Tasks 2–8
> are still the scoped outline at the end; each is written out with tested code in turn, the same way, and this
> document regenerated.

**Goal:** The flowing layer moves: every skin cell carries a velocity in its column's tangent plane, coupled through the
depth by the viscosity at each interface, pushed by the gas pressure gradient, the deceleration and the two thickness
pressures, dragged by the gas shear at the top and held at the rigid front; mass, momentum and heat cross between
neighbouring columns layer by layer.

**Architecture:** Everything happens inside the skin body's existing sub-steps; the 3D model, the moving mesh and the
coupled run are unchanged. Task 1 first makes the sideways conduction implicit: one coupled solve per Newton iteration
of each sub-step, the columns' tridiagonal depth solves as its preconditioner. The flow tasks then give the cells a
velocity, couple it through the depth by an implicit tridiagonal viscous solve per column, and carry mass, enthalpy and
momentum between neighbouring columns in the same implicit sideways solve.

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
which replays the tasks written so far on a copy of the base — every diff through `patch` — and checks that the result
is the tested copy byte for byte. Never edit this document by hand: change `part_21.md` or the tested copy, and
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

## Tasks still in outline (to be written in full)

2. **What sub-plan 20 left to this stage** (its rulings 7 and 9, 2026-10-12): the cells' velocity, stored and carried
   through `Skins`' split, merge, take and add operations with exact momentum; the surface gradient and curvature
   operators on the corner graph (`columns.Columns`), tested against a sphere's exact curvature and a linear field's
   gradient; and the interface rule Asha chooses from sub-plan 20's Task 10 measurements.
3. **Viscosity per interface**: liquid 1.3 mPa·s; Li et al. (2014) slurry law with the bridge near the liquidus;
   fixed-point iteration on the shear rate (at most five passes, held at 367 1/s); rigid cells anchored. Tests against
   the law's tabulated values.
4. **Momentum through the depth**: implicit tridiagonal viscous coupling per column, no slip at the base, the gas shear
   on the top cell. Tests: shear alone gives a linear profile; the push alone a half-parabola; liquid over slurry the
   joined profiles with the kink at the viscosity jump; Stokes' first problem.
5. **The push**: −∇_s p_w − ρ_l a_t + ρ_l a_n ∇_s H − σ ∇_s κ on every layer, with the two thickness terms implicit
   (linearised in the column heights). Tests: a bump on a flat layer under deceleration grows at the thin-film
   Rayleigh–Taylor linear rate; bumps below the surface-tension cutoff decay.
6. **Sideways exchange** (in Task 1's coupled solve): implicit upwind transport of mass, enthalpy and momentum by
   fraction of the flowing depth; the arrival rule (same fraction where the receiver flows from its top down, on top of
   a rigid top otherwise); velocities projected between tangent planes. Tests: a translating pulse keeps its heat and
   arrives where it should; books to round-off; runoff landing on a solid skin freezes from below and the 3D surface
   moves out by the deposit's volume.
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
