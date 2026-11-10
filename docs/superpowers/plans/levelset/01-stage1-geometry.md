# Level-set stage 1 — geometry: the adaptive box mesh, the level sets, the thin test and the shadow test — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (Asha's choice,
> 2026-11-11) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the geometry layer of the level-set architecture — an adaptively refined tetrahedral box mesh, field
transfer between its refinements, cut cells, the evolved skin-base level set β with its exact per-piece volume
correction, the derived outer surface φ and molten-rigid surface ψ, the refinement criteria, the thin test and the
ray-cast shadow test — as a tested package `reentry_model/levelset/`, with prescribed fields and motions only (no
thermal solve, no skins, no flow).

**Architecture:** The box is a regular grid of cubes, each split into six Kuhn tetrahedra; refinement is Maubach's
tagged newest-vertex bisection, which stays conforming, and coarsening is a rebuild from the base grid with the current
targets. Every level set is a P1 nodal field (negative inside); cut cells give each cell's metal volume and its flat
surface piece (a patch); distances are reset geometrically from the patch triangles. β is the only evolved level set;
φ = β moved out by a skin thickness, ψ = the distance to the core-connected metal below 895.1 K, both rebuilt.

**Tech Stack:** Python 3.12 (`drama_env`), numpy 2.5, scipy 1.18 (cKDTree, sparse, csgraph, cg), scikit-fem 12 (only
for `AdaptiveMesh.to_skfem`); pytest. No new dependency (scikit-fmm and numba are not installed and not used).

**Spec:** `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`, §20 (the amendment of
2026-11-10, approved), with §6 (geometry), §7.2 (cut cells) and §13 item 1 (geometry verification); stage 1 of §20.11.

**Prototype:** every module and test below was written and run in `prototype/levelset/stage1/` on 2026-11-10/11
(package `lsproto/`, tests `tests/`, measurement `measure_sphere.py`): 38 tests pass in about 27 s. Each task ports one
prototype module to `reentry_model/levelset/` with its tests, changing only the import prefix `lsproto.` →
`reentry_model.levelset.`. The prototype files are the reference code for each task; the plan states the decisions
they embody, the interfaces and the facts measured with them.

## Measured facts (prototype, 2026-11-10/11, laptop on mains power, shared with a running flight)

1. **Mesh.** Six Kuhn tetrahedra per cube; tagged bisection stays conforming (every face shared by at most two cells,
   boundary area equal to the box's to 1e-12) and keeps at most three tetrahedron shapes. The achieved size is the
   largest h0/2^k not above the target: with h0 = 2.8 mm the fine size 0.7 mm is exact, the floor 0.25 mm becomes
   0.175 mm; the convergence ladder's 0.5 mm needs h0 = 2.0 mm.
2. **Cut cells.** Sphere volume and area converge at second order (errors 2.6 % and 1.3 % at R/h = 4.4, 0.65 % and
   0.33 % at R/h = 8.8; the P1 volume is always low, the interpolant of a convex distance lying above it).
   A contour through nodes leaves zero-area patches, which are dropped (their centroids would otherwise be the origin:
   the cause of the two rear-pole misses in the first shadow run).
3. **Distances.** The first point-to-triangle routine (Ericson's regions) returned NaN on degenerate patch triangles
   where the contour passes through a node; the edge-plus-foot form used now cannot.
4. **Volume correction.** Exact to 1e-13 relative per piece, including two pieces three cells apart, once each
   piece's shift excludes the other pieces' nodes and cells (it did not before). Pieces closer than about two cells
   share faces of metal-holding cells and count as one piece: pieces join and separate at the cell size.
5. **Refinement reach.** A cell's test region is grown by 0.94 of its size (centroid to farthest vertex of a Kuhn
   tetrahedron, 0.935 h); √3 h (tried first) would have widened the fine band by about 1.2 mm beyond d at 0.7 mm cells.
6. **Shadow rays.** Started a twentieth of a cell outside the patch, with steps of φ but at least a quarter of the
   smallest cell. A whole-cell start offset (tried first) carried leeward rays near a sphere's outline outside the body.
7. **Rebuild cost.** Evaluating targets only on leaves created by the previous pass (a leaf that met its target keeps
   meeting it) cut the 50 mm sphere's rebuild from 150 s to 34 s with the same mesh.
8. **The 50 mm sphere** (box 100.8 mm, h0 2.8 mm, fine 0.7 mm over a 4.6 mm band, no melt): 1 290 192 cells
   (896 784 at 0.7 mm), 222 354 nodes, 109 950 patches; peak memory 1.11 GB; volume error before correction 0.053 %,
   area error −0.030 %; times (serial): build 4.6 s, redistance 12.7 s, patches 0.4 s, pieces 6.2 s, volume
   correction 4.2 s, curvature 12.9 s, rebuild from fields 34.5 s, shadow 21.0 s, thin measures 30.4 s — about two
   minutes of geometry per macro step if all run every step.
9. **The 100 mm sphere** (box 151.2 mm, the design's 150 mm box; same settings): 5 355 438 cells (3 926 544 at
   0.7 mm), 915 685 nodes, 439 788 patches; peak memory 2.80 GB; volume error before correction 0.013 %, area error
   −0.0074 %; share of patches shadowed 48.3 % (half the sphere less the grazing band); times (serial): build 16.7 s,
   redistance 89.0 s, patches 1.4 s, pieces 26.5 s, volume correction 18.3 s, curvature 39.1 s, rebuild from fields
   113.7 s, shadow 57.5 s, thin measures 108.5 s — about 7.8 minutes of geometry per macro step if every operation
   ran every step, so about 52 hours over 400 macro steps of a flight (*estimated* from this one measurement). The
   three largest costs (redistance, rebuild, thin measures) all work on the whole mesh where only the fine band
   needs them; see the cost note at the end of this plan.
10. **The 50 mm sphere with the work limited to the fine band** (Task 9's code; same mesh, 1 289 688 cells after the
    rebuild as before): build 2.3 s, redistance 7.9 s, patches 0.1 s, pieces 3.3 s, volume correction 2.2 s, curvature
    6.7 s, rebuild 13.7 s, shadow 7.9 s, thin measures 6.5 s — about 48 s per macro step (every operation but the
    one-off first build) against about 120 s (fact 8, same convention);
    peak memory 1.16 GB; volume and area errors unchanged. Measured once, on the shared laptop.
11. **The 100 mm sphere with the work limited to the fine band** (same mesh, 5 356 410 cells after the rebuild):
    build 13.2 s, redistance 39.9 s, patches 0.7 s, pieces 19.0 s, volume correction 13.3 s, curvature 27.0 s,
    rebuild 60.7 s, shadow 31.9 s, thin measures 35.4 s — about 3.8 minutes per macro step against about 7.8
    (fact 9), so about 25 hours of geometry over 400 macro steps if every operation ran every step (*estimated*);
    peak memory 2.60 GB; errors unchanged. Measured once, on the shared laptop. The rebuild is now the largest cost;
    whether every operation must run every step (the rebuild only when the front nears the edge of its lead, the
    shadow test only when the surface has moved) is measured in stage 3, not decided here.
12. **After the final review (2026-11-11), with the three surface operations included** (outer surface, molten-rigid
    surface and one advection step, each with the band): 50 mm about 63 s per macro step, peak 1.56 GB; 100 mm about
    330 s (5.5 minutes), peak 3.30 GB, about 37 hours over 400 macro steps (*estimated*); the sum counts about one
    distance reset twice (conservative). Measured once on the shared laptop; 43 tests.

## Global Constraints

- Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; tests run from the repository root with
  `"$PY" -m pytest tests/levelset -q`.
- **Isolation (spec §4.1):** new files only, in `reentry_model/levelset/` and `tests/levelset/`; no existing module
  of `reentry_model/` is edited or imported except `numpy`/`scipy`/`skfem`. Task 9 alone adds a level-set section to
  `README.md` and `docs/model_assumptions.md` (spec §14).
- Sign convention everywhere: a level set is negative inside (metal, or rigid metal for ψ); normals point toward
  positive values; SI units (m, K).
- Defaults from the spec: fine size 0.7 mm, margin d 4 mm, lead 0.6 mm, cell floor 0.25 mm (achieved 0.175 mm with
  h0 2.8 mm), curvature factor 0.25, thin-test c_κ 0.1 and c_s 0.2, hold 3, minimum group 7, T_rigid 895.1 K.
- Work on branch `levelset-two-interface` in `.claude/worktrees/levelset-two-interface`. Commit each task only once
  Asha has approved execution; stage named paths only (never `git add -A`); never push without asking. Commit
  messages: a sentence saying what changed, ending `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- The laptop is shared with running flights: run tests only on mains power (`pmset -g batt`), run only the task's
  test file during a task, and the whole `tests/levelset` tier once at the end of each task.
- Never edit the `subplan20` worktree, `prototype/work-2026-10-28-skin-spraying*`, `code/reentry_model_output/`, the
  DRAMA databases or the wrapper scripts.

## Review Focus

- **A contour through nodes or along cell faces** (common for analytic shapes aligned with the grid): no NaN, no
  zero-area patch, the exact plane area — pinned by `test_a_contour_through_nodes_leaves_no_zero_area_patches`
  (Task 3) and the degenerate-triangle cases of `test_point_triangle_distance_in_every_region` (Task 4).
- **Departure points and rays outside the box:** finite values, rays stop at the box —
  `test_departure_points_outside_the_box_give_finite_values` (Task 4) and the shadow tests (Task 8).
- **Pieces close together:** each piece's correction leaves the others' volumes exact —
  `test_two_pieces_three_cells_apart_are_corrected_independently` (Task 4); closer than two cells they are one piece
  (stated in the docstring of `pieces`).
- **No rigid metal at all (fully molten):** ψ is +inf everywhere and nothing crashes downstream —
  `test_fully_molten_metal_has_no_rigid_region` (Task 5); `make_target` treats ψ = +inf as "all metal fine".
- **A piece whose account reaches zero:** `correct_volumes` refuses it with a message naming the piece rather than
  looping — `test_a_zero_volume_target_is_refused` (Task 4).

---

### Task 1: The box mesh and conforming refinement

**Files:**
- Create: `reentry_model/levelset/__init__.py` (docstring only: the level-set architecture, spec path), `reentry_model/levelset/boxmesh.py` (from `prototype/levelset/stage1/lsproto/boxmesh.py`)
- Test: `tests/levelset/test_levelset_boxmesh.py` (from `prototype/levelset/stage1/tests/test_boxmesh.py`)

**Interfaces:**
- Produces: `BoxBase(lo: tuple, h0: float, n: tuple)` frozen dataclass with `BoxBase.around(lo, hi, h0)` and
  property `hi`; `AdaptiveMesh(base, points (N,3), tets (M,4) int64 Maubach-ordered, tag (M,) int8 in 1..3,
  level (M,) int16)` with properties `size` (h0·2^-(level//3)), `centroids`, `volumes` (cached), method `to_skfem()`,
  methods `node_cells()` (cached CSR nodes × cells) and `ring(cells, rings) -> sorted cells` (cells within `rings`
  rings through shared nodes), and a private `_cache` dict that `transfer` uses; `base_mesh(base) -> AdaptiveMesh`;
  `build(base, target_size=None, max_level=30, max_passes=400) -> AdaptiveMesh`, where
  `target_size(centroids (k,3), sizes (k,)) -> (k,)` largest allowed size, `np.inf` for none.

- [ ] **Step 1: Write the failing tests** — the five tests of the prototype file: six Kuhn tetrahedra per cube and a
  conforming base; refinement toward a point conforming and meeting the target; achieved size the largest halving not
  above the target (0.4 mm asked of h0 1 mm gives 0.25 mm everywhere); at most three shapes with volume/longest-edge³
  above a quarter of the Kuhn value (1/6)/3^1.5; the same targets give the same mesh.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_boxmesh.py -q` — Expected: FAIL (module missing).
- [ ] **Step 3: Implement `boxmesh.py`.** Kuhn split: for each permutation (a, b, c) of the axes, the tetrahedron
  (v, v+e_a, v+e_a+e_b, v+e_a+e_b+e_c) of every cube's corner v, tag 3. Bisection of a leaf (x0, x1, x2, x3) with tag
  k cuts edge x0–x_k at its midpoint z; children (x0..x_{k−1}, z, x_{k+1}..x3) and (x1..x_k, z, x_{k+1}..x3), tag
  k−1, or 3 after 1 (Maubach 1995). Midpoints are shared through sorted int64 edge keys (lo·2³¹ + hi). Closure: a
  leaf with any edge whose midpoint exists is hanging and is bisected in the next pass. Targets are evaluated only on
  leaves created by the previous pass (fact 7); a leaf is marked while level < 3·ceil(log2(h0/target)).
- [ ] **Step 4: Run** the same command — Expected: 5 passed.
- [ ] **Step 5: Commit** `reentry_model/levelset/__init__.py reentry_model/levelset/boxmesh.py tests/levelset/test_levelset_boxmesh.py`.

### Task 2: Point location and field transfer

**Files:**
- Create: `reentry_model/levelset/transfer.py` (from `lsproto/transfer.py`)
- Test: `tests/levelset/test_levelset_transfer.py` (from `tests/test_transfer.py`)

**Interfaces:**
- Consumes: `AdaptiveMesh` (Task 1).
- Produces: `_geom(mesh) -> (cKDTree of centroids, Jinv (M,3,3), x0)` cached in `mesh._cache`;
  `locate(mesh, pts, k=4) -> (cell (n,), bary (n,4))` (k nearest centroids, then 4k, 16k, 64k for the rest;
  points outside the box get the cell with the least negative coordinate); `barycentric(mesh, cells, pts)`;
  `interpolate(mesh, u, pts) -> (n,)`; `mass_matrix(mesh)` (consistent P1, csr); `project(old, new, u_old,
  rtol=1e-13) -> u_new`; `integral(mesh, u) -> float`.

- [ ] **Step 1: Write the failing tests** — locate gives non-negative barycentrics reproducing the points to 1e-15 m;
  a linear field interpolates exactly (rtol 1e-12); projection onto a refinement equals interpolation (atol 1e-10);
  projection onto a coarsening conserves the integral (rel 1e-12).
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_transfer.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `transfer.py`.** `project` assembles the exact right-hand side on the finer mesh of each
  region: new leaves at least as deep as the old leaf holding their centroid are integrated on themselves (old field
  located at their quadrature points); old leaves deeper than the new leaf holding their centroid are integrated on
  themselves (new basis located at their points); 4-point rule exact for quadratics (a = 0.5854101966249685,
  b = 0.1381966011250105); then `scipy.sparse.linalg.cg` on the consistent mass matrix from the interpolated guess,
  `rtol`, `atol=0`, at most 2000 iterations, RuntimeError if not converged.
- [ ] **Step 4: Run** — Expected: 4 passed.
- [ ] **Step 5: Commit** the two files.

### Task 3: Cut cells and patches

**Files:**
- Create: `reentry_model/levelset/cutcell.py` (from `lsproto/cutcell.py`)
- Test: `tests/levelset/test_levelset_cutcell.py` (from `tests/test_cutcell.py`)

**Interfaces:**
- Produces: `Patches` dataclass (`cell`, `area`, `centroid`, `normal`, `tris (Q,3,3)`, `tri_patch (Q,)`);
  `cell_gradients(mesh, phi, cells=None) -> (k,3)` (from the cached inverse edge matrices of `transfer._geom`, a
  multiply instead of a solve; `cutcell` imports `transfer` inside the function); `inside_volumes(mesh, phi, cells=None) -> (k,)` volume of
  {φ < 0}; `patches(mesh, phi, cells=None) -> Patches` (zero-area pieces dropped).

- [ ] **Step 1: Write the failing tests** — inside plus outside equals each cell for random fields (rtol 1e-12); a
  plane gives the exact normal (atol 1e-12) and the same volume and area on a refined and an unrefined mesh
  (rel 1e-12); a sphere's volume and area errors below 1 % and 0.5 % at R/h = 8.8 and falling more than 3.5 times
  per halving; normals point outward (cosine with the radius above 0.9) and centroids lie within h²/R of the sphere;
  a contour through a layer of nodes leaves no zero-area patch and the exact area 4e-6 m².
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_cutcell.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `cutcell.py`.** A node with φ exactly 0 counts as outside. One (three) negative vertices:
  the corner tetrahedron's volume is V·Π t_j, t_j = φ_a/(φ_a − φ_j) over the lone vertex's three edges (inside =
  corner, or V − corner). Two and two: the prism with bottom (a, p_ac, p_ad) and top (b, p_bc, p_bd) as the three
  tetrahedra (a, p_ac, p_ad, b), (p_ac, p_ad, b, p_bc), (p_ad, b, p_bc, p_bd). The patch is the triangle of the three
  cut points, or the quad (p_ac, p_bc, p_bd, p_ad) as two triangles; its normal is ∇φ/|∇φ| of the cell.
- [ ] **Step 4: Run** — Expected: 5 passed.
- [ ] **Step 5: Commit** the two files.

### Task 4: The evolved level set β — distances, motion, pieces and the exact volume correction

**Files:**
- Create: `reentry_model/levelset/levelset.py` (from `lsproto/levelset.py`)
- Test: `tests/levelset/test_levelset_levelset.py` (from `tests/test_levelset.py`)

**Interfaces:**
- Consumes: `patches`, `inside_volumes` (Task 3); `interpolate` (Task 2).
- Produces: `point_triangle_distance(p, a, b, c) -> (n,)`; `redistance(mesh, phi, k=8, band=None) -> phi` (±inf if
  there is no contour; with `band` (m) only nodes within band of a patch-triangle centroid are computed, the rest
  get ±band); `advect(mesh, phi, velocity, dt)` with `velocity(x (n,3)) -> (n,3)` m/s; `pieces(mesh, phi) -> (label
  per cell, −1 without metal; n)`, pieces numbered by decreasing volume; `piece_volumes(mesh, phi, label, n)`;
  `correct_volumes(mesh, phi, label, targets, rings=2, rtol=1e-14, max_iter=60) -> (phi, largest relative
  correction)`; ValueError naming the piece for a non-positive target.

- [ ] **Step 1: Write the failing tests** — the eight tests of the prototype file: point-triangle distance in every
  region (interior 0.5, vertex √2, vertex 1, edge 1, edge √2/2, edge 1); redistancing 7·(r − R) recovers r − R within
  0.2 h; six semi-Lagrangian steps of a 1 mm/s radial recession over 0.05 s each give patch radii within h of 0.8 mm;
  the correction exact to 1e-13 with a first error below 1 %; a dumbbell (two 0.7 mm spheres, a 0.2 mm neck) receding
  0.3 mm gives two pieces corrected to 1e-13; two pieces three cells apart corrected independently to 1e-13; a zero
  target refused; departure points 1 mm outside the box give finite values.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_levelset.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `levelset.py`.** Distance: the nearest of the three edges, or the perpendicular distance
  where the foot falls inside a non-degenerate triangle (fact 3); `redistance` searches the k triangles with the
  nearest centroids (cKDTree) and keeps the old sign. `advect`: φ(x − u(x)Δt) by `interpolate`, then `redistance`.
  `pieces`: cells with metal, joined through shared faces (sorted face triples), scipy `connected_components`.
  `correct_volumes`: per piece, the nodes of its cells dilated by `rings` rings, minus every node of another piece;
  the volume is evaluated over cells touching those nodes that hold no other piece; Newton on the uniform shift c
  (φ − c on those nodes) with the patch area as derivative, safeguarded by bisection between bracketing shifts, until
  |target − V| ≤ rtol·target.
- [ ] **Step 4: Run** — Expected: 8 passed.
- [ ] **Step 5: Commit** the two files.

### Task 5: The derived surfaces φ and ψ

**Files:**
- Create: `reentry_model/levelset/surfaces.py` (from `lsproto/surfaces.py`)
- Test: `tests/levelset/test_levelset_surfaces.py` (the four surface tests of `tests/test_surfaces.py`)

**Interfaces:**
- Consumes: `redistance`, `pieces` (Task 4); `locate` (Task 2).
- Produces: `T_RIGID = 895.1`; `outer_surface(mesh, beta, thickness) -> phi`;
  `molten_rigid_surface(mesh, phi, T, T_rigid=T_RIGID, core_point=None) -> (psi, crust_cells bool per cell)`.

- [ ] **Step 1: Write the failing tests** — the outer surface stands a 0.2 mm thickness outside a 0.9 mm base (patch
  radii within 0.02 mm of 1.1 mm); a radial temperature 850 K + 60 K/mm gives the rigid sphere of radius 0.7517 mm
  within 0.02 mm, no crust, ψ ≥ φ; a cold shell over a molten layer over a cold core is a crust and ψ follows the core
  (within 0.1 mm of 0.6 mm); fully molten metal gives ψ = +inf with no crust.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_surfaces.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `surfaces.py`.** `outer_surface` = `redistance(beta − thickness)`. `molten_rigid_surface`:
  g = `redistance(T − T_rigid)`; q = max(g, φ); pieces of q; keep the piece holding `core_point` (default the node of
  least φ; the largest piece if that point is not rigid); nodes only in other pieces (crusts) get q → −q; ψ =
  max(`redistance(q)`, φ). With no node below T_rigid, ψ = +inf.
- [ ] **Step 4: Run** — Expected: 4 passed.
- [ ] **Step 5: Commit** the two files.

### Task 6: Curvature and the refinement targets

**Files:**
- Create: `reentry_model/levelset/criteria.py` (from `lsproto/criteria.py`)
- Test: `tests/levelset/test_levelset_criteria.py` (the four curvature and refinement tests of `tests/test_surfaces.py`)

**Interfaces:**
- Consumes: `cell_gradients` (Task 3); `locate` (Task 2); `build` (Task 1) in tests.
- Produces: `nodal_average(mesh, cell_values, cells=None)` (NaN at nodes outside the given cells);
  `hessian_field(mesh, phi, cells=None) -> (G (N,3), H (N,3,3))`; `principal_curvatures(mesh, phi, pts, fields=None,
  rings=3) -> (k1, k2)`, k1 ≥ k2, positive where convex — without `fields` the recovery runs only on
  `mesh.ring(cells holding pts, rings)`; `kappa_max_nodal(mesh, phi, band=None, fields=None) -> (N,)` (with `band`,
  only nodes with |φ| < band, others 0); `make_target(mesh, phi, psi, kmax_phi, kmax_psi, fine=0.7e-3,
  margin=4.0e-3, lead=0.6e-3, floor=0.25e-3, curvature_factor=0.25) -> target_size` for `build`.

- [ ] **Step 1: Write the failing tests** — a sphere of 1.1 mm: median k1 and k2 within 5 % of 1/R at h = 0.0625 mm;
  a 0.8 mm cylinder: k1 within 5 % of 1/r, |k2| below 5 % of it; the fine region reaches the margin below the floor
  (every centroid between 0.65 and 1.25 mm from the centre at the fine size; the core inside 0.1 mm at the base size);
  curvature refines a 0.2 mm wire to the floor within 0.05 mm of its surface.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_criteria.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `criteria.py`.** Hessian by gradient recovery twice (cell gradients averaged to nodes by
  volume, then the cell gradients of each component, symmetrised, averaged again). At a point: S = P H P / |∇φ| with
  P = I − nn; drop the eigenvalue whose eigenvector is most nearly n. Targets: a cell is in the fine region when
  φ < r and ψ > −(margin + lead + r), r = 0.94·size (fact 5); within 2r of either surface the target is also at most
  max(floor, curvature_factor / κ_max).
- [ ] **Step 4: Run** — Expected: 4 passed.
- [ ] **Step 5: Commit** the two files.

### Task 7: The thin test with hysteresis

**Files:**
- Create: `reentry_model/levelset/thin.py` (from `lsproto/thin.py`)
- Test: `tests/levelset/test_levelset_thin.py` (the four thin tests of `tests/test_thin_shadow.py`)

**Interfaces:**
- Consumes: `hessian_field`, `nodal_average`, `principal_curvatures` (Task 6); `cell_gradients`, `Patches` (Task 3);
  `locate` (Task 2).
- Produces: `thin_candidate(h, kmax, slope, c_kappa=0.1, c_slope=0.2) -> bool array`; `patch_measures(mesh, phi, psi,
  P) -> (h, kmax, slope)` (h = max(ψ, 0) at the patch; kmax the largest |principal curvature| of φ at the patch and
  of ψ at the foot x − h n; slope = |tangential part of ∇ψ|, the gradient recovered on `mesh.ring(cells, 1)` only); `patch_adjacency(mesh, P) -> csr` (cells sharing a
  node); `ModeState(thin, hold=3, min_group=7)` with `.thin`, `.count`, `update(candidate, adjacency) -> changed`.
  Stage 3 keys `ModeState` by marker instead of patch.

- [ ] **Step 1: Write the failing tests** — the decision uses the largest principal curvature (0.3 mm on a ±1/mm
  saddle is not thin; 0.05 mm on 1/mm is; slope 0.3 fails); a 0.06 mm shell is thin on more than 95 % of patches and
  a 0.5 mm pool on a 1.1 mm sphere thick on more than 95 %; a floor tilted at 0.3 gives slope 0.3/√1.09 within 0.03;
  modes wait three updates and switch only a connected group of at least seven (eight switch, a lone one does not).
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_thin.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `thin.py`** with the signatures above; in `ModeState.update`, a patch whose candidate differs
  counts up (reset to 0 when it agrees); those at `hold` form groups through the adjacency restricted to them, and
  groups of at least `min_group` switch and reset.
- [ ] **Step 4: Run** — Expected: 4 passed.
- [ ] **Step 5: Commit** the two files.

### Task 8: The ray-cast shadow test

**Files:**
- Create: `reentry_model/levelset/shadow.py` (from `lsproto/shadow.py`)
- Test: `tests/levelset/test_levelset_shadow.py` (the two shadow tests of `tests/test_thin_shadow.py`)

**Interfaces:**
- Consumes: `locate` (Task 2); `Patches` (Task 3).
- Produces: `shadowed(mesh, phi, P, direction, max_steps=400) -> bool (P,)`; `direction` is the direction of motion
  (upstream); RuntimeError if rays have not finished in `max_steps`.

- [ ] **Step 1: Write the failing tests** — a sphere moving along +x shadows every patch with n_x < −0.2 and none with
  n_x > 0.2; a stepped stem under a 1 mm cap: the forward-facing ledge annulus (0.4–0.7 mm from the axis at x = 1 mm)
  is shadowed and the cap's front face is not.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_shadow.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement `shadow.py`.** Start at centroid + 0.05·(cell size)·n (fact 6); march along the direction
  by max(φ, a quarter of the smallest cell); φ < 0 marks shadow; leaving the box ends the ray.
- [ ] **Step 4: Run** — Expected: 2 passed.
- [ ] **Step 5: Commit** the two files.

### Task 9: Work limited to the fine band, pinned and measured

Asha's decision of 2026-11-11 (after fact 9): the operations that did whole-mesh work only the fine band needs are
limited to it. The code is already in the modules of Tasks 1–8 (the interfaces above); this task pins that the
limited work gives the same answers, and measures the cost. Profiled first (50 mm sphere, cProfile): the rebuild's time
was point location, the thin measures' was a 3 × 3 solve per cell over every cell, twice per surface, and the
redistance's was the nearest-triangle search, dearest for nodes far from the surface.

**Files:**
- Create: `analysis/levelset/measure_geometry.py` (from `prototype/levelset/stage1/measure_sphere.py`, imports from
  `reentry_model.levelset`; distances reset with `band` = 4.6 mm + 2 h0, curvature with `band` = 1.88 h0)
- Test: `tests/levelset/test_levelset_band.py` (the two band tests of `tests/test_levelset.py` and
  `tests/test_surfaces.py`)

**Interfaces:**
- Consumes: every module above. Usage: `"$PY" analysis/levelset/measure_geometry.py <diameter_mm> [h0_mm] [fine_mm]`,
  one JSON line (cells, nodes, patches, times per operation, peak memory, volume and area errors, shadowed share).

- [ ] **Step 1: Write the tests** — a banded reset (band 0.3 mm) equals the full reset exactly where |φ| < 0.25 mm,
  never exceeds it elsewhere, and reaches exactly 0.3 mm; the local curvature at the patches equals the whole-mesh
  recovery to 1e-12, and `kappa_max_nodal(band=0.3 mm)` is 0 where |φ| ≥ 0.3 mm and positive where |φ| < 0.2 mm.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/levelset/test_levelset_band.py -q` — Expected: 2 passed (the code
  exists from Tasks 1–8; if either fails, the port differs from the prototype).
- [ ] **Step 3: Run the whole tier** `"$PY" -m pytest tests/levelset -q` — Expected: 38 passed.
- [ ] **Step 4: Run the measurement** on mains power for 50 and 100 mm — Expected: cells and memory within 20 % of
  facts 10 and 11 (times vary with the laptop's load).
- [ ] **Step 5: Commit** the two files.

### Task 10: Documentation

**Files:**
- Modify: `README.md` (a "Level-set architecture" section: what exists after stage 1, the test command, the measured
  costs), `docs/model_assumptions.md` (a level-set section: Kuhn mesh and bisection, P1 level sets, the distance
  reset and its band, the per-piece correction, pieces joining below two cells, the thin-test thresholds as starting
  values), `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §20.13 (the measured costs,
  facts 8–10) and §20.8 (the ray start offset of fact 6)

- [ ] **Step 1: Write the documentation** as listed, every number with what was measured and when.
- [ ] **Step 2: Commit** the named files.

## Cost note

Fact 9 put the geometry at about 7.8 minutes per macro step on the 100 mm sphere; Asha chose (2026-11-11) to limit the
whole-mesh work to the fine band before execution (Task 9; facts 10 and 11 are the result).
