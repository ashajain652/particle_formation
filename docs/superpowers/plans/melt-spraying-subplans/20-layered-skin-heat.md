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

```diff
--- a/tests/test_reentry_model_movemesh.py
+++ b/tests/test_reentry_model_movemesh.py
@@ -287,10 +287,9 @@
     assert moving.mesh.volume() == pytest.approx(V, rel=1e-12)
     assert solver.energy() == pytest.approx(E, rel=1e-12)
     assert np.all(movemesh.signed_volumes(moving.mesh.points, moving.mesh.tets) != 0.0)
-    # the interpolation loses -2612 J of 71236 J (-4.1 %) here: a kinked field on 4 to 20 mm cells, whose lumped sum
-    # differs between two meshes (final review: the earlier one-sided check "residual below 1e-3" passed on the sign
-    # alone). The residual is reported, and the offset closes it (the energy assertion above)
-    assert abs(result.relative_residual) < 0.06 and moving.rebuilds == 1 and solver.mesh is moving.mesh
+    # a kinked field on 4 to 20 mm cells: the projection leaves -1.8 J of 63 219 J (3e-5) for the offset to close, where
+    # interpolating at the new nodes left -4.1 % (both measured, sub-plan 20)
+    assert abs(result.relative_residual) < 1e-3 and moving.rebuilds == 1 and solver.mesh is moving.mesh
 
 
 def test_a_rebuild_of_a_linear_field_loses_nothing_to_interpolation(coarse_sphere_mesh):
@@ -302,8 +301,7 @@
     moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3, sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
     assert not moving.move(inward(moving, 1e-3)).refused
     result = moving.rebuild()
-    # linear interpolation and the lumped sum are both exact for a linear field: -6e-8 of the total measured (final
-    # review), two-sided
+    # the projection and the lumped sum are both exact for a linear field (two-sided, to the solve's tolerance)
     assert abs(result.relative_residual) < 1e-6
     # the books' own uniform offsets (the move's, then the rebuild's) leave 2.1e-5 K at most (measured)
     assert np.abs(solver.temperature() - T_of(moving.mesh.points)).max() < 1e-4
@@ -345,22 +343,23 @@
         old_points, old_tets = moving.mesh.points.copy(), moving.mesh.tets.copy()
         h_old = mat.enthalpy(solver.temperature())
         residuals.append(moving.rebuild().residual)
-        # what the interpolation alone gives at the new nodes (the corners keep their own values, and are in the zone)
-        elem, bary = movemesh.locate(old_points, old_tets, moving.mesh.points)
-        h_interp = np.einsum("ij,ij->i", bary, h_old[old_tets[elem]])
+        # what the projection alone gives at the new nodes
+        h_proj = movemesh.project_enthalpy(old_points, old_tets, h_old, moving.mesh.points, moving.mesh.tets)
         h_after = mat.enthalpy(solver.temperature())
         in_zone = np.zeros(moving.mesh.n_nodes, dtype=bool)
         in_zone[moving.zone.corners] = True
         in_zone[moving.zone.interior] = True
         out = ~in_zone
-        # outside the zone (84 and 73 of about 2900 nodes on this coarse mesh) the field is the interpolation's to 3e-11
-        # J per kg; spread over the whole body the offset would have shifted those nodes by 1735 and 1244 J per kg (1.7
-        # and 1.2 K), and in the zone, noncorner nodes carry 2585 and 1804 J per kg (all measured, final review)
+        # outside the zone (about 80 of 2900 nodes on this coarse mesh) the field is the projection's; inside, every node
+        # carries the same offset, the residual over the zone's mass
         assert out.sum() > 50
-        assert np.abs(h_after[out] - h_interp[out]).max() < 1e-6 * np.abs(h_interp).max()      # tolerance 0.08 J per kg
-        zone_inner = in_zone & ~np.isin(np.arange(moving.mesh.n_nodes), moving.corners)
-        assert np.abs(h_after[zone_inner] - h_interp[zone_inner]).min() > 1000.0
-    assert min(abs(x) for x in residuals) > 1.0
+        assert np.abs(h_after[out] - h_proj[out]).max() < 1e-6 * np.abs(h_proj).max()
+        shift = (h_after - h_proj)[in_zone]
+        assert np.ptp(shift) < 1e-6 * np.abs(h_proj).max() + 1e-9 * abs(shift).max()
+        assert np.sign(shift.mean()) == np.sign(residuals[-1]) or residuals[-1] == 0.0
+    # two-signed and small: 0.007 J and 19.7 J here, of 67 000 J (measured, sub-plan 20; interpolating at the new nodes
+    # left one-signed residuals of hundreds to thousands of J)
+    assert max(abs(x) for x in residuals) < 1e-3 * solver.energy()
     assert moving.rebuild_residual_total == pytest.approx(sum(residuals), rel=1e-12, abs=0.0)
 
 
@@ -479,3 +478,51 @@
     theta = (solver.temperature() - 300.0) / 300.0
     exact = exact_front(x - x_front, 0.5, 1e-5, 1e-2)
     assert np.abs(theta - exact).max() < 0.005          # 0.0010 measured while planning: 0.1 % of the rise
+
+
+# -- volume-neutral smoothing (sub-plan 20, Task 0) ----------------------------------------------------------------
+def test_neutral_smoothing_keeps_the_tangential_step_and_cancels_most_of_its_sweep(coarse_sphere_mesh):
+    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
+    mat = constant_material()
+    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat, depth=15e-3)
+    normal, area = moving.frame()
+    points, faces, n = moving.mesh.points[moving.corners], moving.faces_corner, len(moving.corners)
+    plain = moving.smoothing(normal)
+    neutral = moving.neutral_smoothing(normal, area)
+    # the passes add only normal motion: the tangential step is the plain smoothing's
+    extra = neutral - plain
+    assert np.abs(extra - np.einsum("ij,ij->i", extra, normal)[:, None] * normal).max() < 1e-15
+    before = movemesh.swept_volumes(points, points + plain, faces, faces, n) / area
+    after = movemesh.swept_volumes(points, points + neutral, faces, faces, n) / area
+    # equivalent depths of the sweep per corner; on this 4 mm sphere the root mean square falls from 0.78 to 0.27 um
+    # and the worst corner from 11.4 to 3.9 um (measured while planning)
+    assert np.sqrt((after ** 2).mean()) < 0.45 * np.sqrt((before ** 2).mean())
+    assert np.abs(after).max() < 0.5 * np.abs(before).max()
+
+
+def test_a_neutral_smoothing_move_keeps_the_books(coarse_sphere_mesh):
+    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
+    mat = constant_material()
+    solver = solver_on(m, mat)
+    solver.set_temperature(300.0 + 2000.0 * m.points[:, 2])
+    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
+    E0, V0 = solver.energy(), m.volume()
+    r = moving.move(inward(moving, 1e-4) + moving.neutral_smoothing())
+    assert not r.refused
+    assert r.swept_volume.sum() == pytest.approx(V0 - moving.mesh.volume(), rel=1e-12)
+    assert solver.energy() + r.swept_enthalpy.sum() == pytest.approx(E0, rel=1e-12)
+
+
+def test_a_corner_that_does_not_move_is_credited_with_nothing_its_neighbours_swept(coarse_sphere_mesh):
+    """Own sweeps (sub-plan 20, spec 5.2 amended 2026-10-12): one corner moving in sweeps its own triangles' prisms, and
+    every other corner books nothing; a third of each prism to each corner credited its neighbours with two thirds."""
+    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
+    mat = constant_material()
+    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat, depth=15e-3)
+    normal, area = moving.frame()
+    points, faces, n = moving.mesh.points[moving.corners], moving.faces_corner, len(moving.corners)
+    u = np.zeros((n, 3))
+    u[7] = -2e-4 * normal[7]
+    dV = movemesh.swept_volumes(points, points + u, faces, faces, n)
+    assert dV[7] == pytest.approx(area[7] * 2e-4, rel=0.02)
+    assert np.abs(np.delete(dV, 7)).max() == 0.0
```

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`
Expected: the three new tests fail (`AttributeError: 'MovingMesh' object has no attribute 'neutral_smoothing'`, and
the own-sweep test's neighbours are credited two thirds of the moving corner's prism), and the two rebuild tests fail on
their residuals (the interpolation's −4.1 %).

- [ ] **Step 3: Apply the change to `reentry_model/movemesh.py`**

```diff
--- a/reentry_model/movemesh.py
+++ b/reentry_model/movemesh.py
@@ -51,6 +51,9 @@
 ZONE_NEIGHBOURS = 12         # surface nodes averaged for each interior node of the zone
 SMOOTH_LAMBDA = 0.5          # tangential smoothing: share of the way to the neighbours' mean per move
 SMOOTH_CLAMP = 0.35          # of the mean surface edge: the most a corner is nudged per move (2026-09-27 design)
+NEUTRAL_PASSES = 3           # normal corrections that make the smoothing volume-neutral per corner (sub-plan 20 Task 0):
+                             # on the 100 mm sphere each pass removes about 30 % of what is left, the worst corner's
+                             # sweep going from 1.0 to 0.37 um of equivalent depth after three (measured 2026-10-12)
 _REF_GRADS = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
 
 
@@ -185,20 +188,28 @@
 
 def swept_volumes(old, new, faces, faces_corner, n_corners):
     """The volume each corner's surface swept, positive inward: every surface triangle's prism between its old and new
-    positions, split into three tetrahedra, a third of it to each of the triangle's corners. A corner sliding along a
-    flat face sweeps exactly nothing; a corner on a sharp edge books only what the faces it touches swept (Task 4
-    review, 2026-10-12, replacing dual area times normal displacement, which mis-shares at edges and, scaled to the
-    total, blew up wherever the first-order volumes nearly cancelled). `faces` are the outward triangles' node ids,
-    `faces_corner` the same triangles in corner numbering."""
+    positions (exact, split into three tetrahedra), shared among the triangle's corners by what each corner's own motion
+    swept -- to first order A_t (u_k . n_t) / 3 for corner k, n_t the triangle's outward normal -- with the second-order
+    rest of the prism shared in proportion to those. A corner sliding along a flat face sweeps exactly nothing, and a
+    corner that does not move is never credited with what its moving neighbours swept: booking a third of each prism to
+    every corner (sub-plan 19) credited the side walls of the box test with metal its receding front swept, they handed
+    it down, bulged out and folded the surface at 1.3 s (sub-plan 20, 2026-10-12). `faces` are the outward
+    triangles' node ids, `faces_corner` the same triangles in corner numbering."""
     a, b, c = (old[faces[:, i]] for i in range(3))
     A, B, C = (new[faces[:, i]] for i in range(3))
 
     def vol(p, q, r, s):
         return np.einsum("ij,ij->i", q - p, np.cross(r - p, s - p)) / 6.0
     prism = -(vol(a, b, c, A) + vol(b, c, A, B) + vol(c, A, B, C))
+    n2 = np.cross(b - a, c - a)                                   # twice the area times the outward unit normal
+    own = np.stack([-np.einsum("ij,ij->i", Q - P, n2) / 6.0 for P, Q in ((a, A), (b, B), (c, C))], axis=1)
+    weight = np.abs(own)
+    total = weight.sum(axis=1)
+    share = np.where(total[:, None] > 0.0, weight / np.where(total > 0.0, total, 1.0)[:, None], 1.0 / 3.0)
+    booked = own + (prism - own.sum(axis=1))[:, None] * share
     out = np.zeros(n_corners)
     for i in range(3):
-        np.add.at(out, faces_corner[:, i], prism / 3.0)
+        np.add.at(out, faces_corner[:, i], booked[:, i])
     return out
 
 
@@ -292,6 +303,23 @@
         return tangential_smoothing(self.mesh.points[self.corners], self.faces_corner, normal,
                                     SMOOTH_LAMBDA if lam is None else lam, SMOOTH_CLAMP if clamp is None else clamp,
                                     self.edges_corner)
+
+    def neutral_smoothing(self, normal=None, area=None, passes=NEUTRAL_PASSES):
+        """The tangential smoothing made volume-neutral per corner (layered-skin design section 5.1, decided 2026-10-12):
+        each pass moves every corner along its normal by the volume the displacement so far sweeps there, divided by
+        its dual area, so that a column is not credited with metal the smoothing alone swept. A corner's own normal
+        move books only a third of its prisms to itself, so the passes converge slowly for patterns that alternate
+        from corner to corner and never for a pattern summing to zero on every triangle; what is left is booked to the
+        columns like any other sweep (measured 2026-10-12 on the 100 mm sphere: the worst corner 1.03 um of equivalent
+        depth before, 0.37 um after three passes; the root mean square 0.115 um before, 0.038 um after)."""
+        if normal is None or area is None:
+            normal, area = self.frame()
+        points = self.mesh.points[self.corners]
+        u = self.smoothing(normal)
+        for _ in range(passes):
+            dV = swept_volumes(points, points + u, self.faces_corner, self.faces_corner, len(self.corners))
+            u = u + normal * (dV / area)[:, None]
+        return u
 
     def surface_quality(self):
         """The worst surface triangle's quality, 2 r_in / R_circ."""
@@ -357,9 +385,7 @@
         points, tets, node_of_corner = rebuild_volume(mesh.points[self.corners], self.faces_corner, **self.sizes)
         h_old = mat.enthalpy(solver.temperature())
         E_old = solver.energy()
-        elem, bary = locate(mesh.points, mesh.tets, points)
-        h_new = np.einsum("ij,ij->i", bary, h_old[mesh.tets[elem]])
-        h_new[node_of_corner] = h_old[self.corners]          # the surface is the same: corners keep their enthalpy
+        h_new = project_enthalpy(mesh.points, mesh.tets, h_old, points, tets)
         params = dict(mesh.params, rebuilds=self.rebuilds + 1)
         new_mesh = mesh_mod.VolumeMesh(points, tets, params)
         check_rebuilt_boundary(new_mesh, node_of_corner[self.faces_corner])
@@ -368,9 +394,8 @@
         residual = E_old - solver.energy()
         self.mesh, self.corners = new_mesh, node_of_corner
         self._set_surface(node_of_corner[self.faces_corner])             # builds the new deforming zone
-        # the interpolation error arises where the field varies, near the surface (its residuals measured one-signed,
-        # 99 to 186 J per rebuild), so the offset that makes the total exact is spread over the deforming zone only,
-        # by lumped mass, and the deep body keeps the interpolated values untouched
+        # what the quadrature misses arises where the field varies, near the surface, so the offset that makes the
+        # total exact is spread over the deforming zone only, by lumped mass, and the deep body keeps its values
         in_zone = np.concatenate([self.zone.corners, self.zone.interior])
         m_zone = lumped_masses(points, tets, mat.rho)[in_zone]
         h_new[in_zone] += residual / m_zone.sum()
@@ -513,6 +538,45 @@
     return coords, tets, node_of_corner.astype(np.int64)
 
 
+_QUAD_A, _QUAD_B = 0.5854101966249685, 0.1381966011250105      # the 4-point rule exact for quadratics on a tetrahedron
+QUADRATURE = np.full((4, 4), _QUAD_B) + np.eye(4) * (_QUAD_A - _QUAD_B)   # barycentric coordinates, weights 1/4 each
+PROJECT_CHUNK = 200_000                                      # quadrature points located at a time (bounds the memory)
+
+
+def project_enthalpy(old_points, old_tets, h_old, points, tets):
+    """The nodal enthalpy on a new mesh as the L2 projection of the old P1 field: M h = b, M the new mesh's consistent
+    mass matrix and b_i the integral of phi_i times the old field over the new elements by the 4-point rule. The
+    projection reproduces a linear field exactly and keeps the integral of h -- the energy -- to the rule's accuracy,
+    region by region, because the rows of M add up to the lumped masses the solvers' energy uses (sub-plan 20, the spec's
+    fallback once the residuals proved large: interpolating at the new nodes and offsetting the total uniformly over the
+    zone moved 2 768 J -- 4.9 % of the heat absorbed -- from the cold bulk to the melting front of the level-2 plate in two
+    rebuilds, and its centre melted 9.5 % ahead of the exact solution; a lumped projection kept the energy local but
+    moved a linear field's surface values by up to 7 K). An old field point just outside the old mesh -- a quadrature
+    point across the shared surface -- is extrapolated linearly from the nearest element, as `locate` does."""
+    import scipy.sparse as sp
+    import scipy.sparse.linalg as spla
+    x = points[tets]                                         # (ne, 4, 3)
+    vol = np.abs(signed_volumes(points, tets))
+    qp = np.einsum("qa,eai->eqi", QUADRATURE, x).reshape(-1, 3)
+    h_q = np.empty(len(qp))
+    for start in range(0, len(qp), PROJECT_CHUNK):
+        elem, bary = locate(old_points, old_tets, qp[start:start + PROJECT_CHUNK])
+        h_q[start:start + PROJECT_CHUNK] = np.einsum("ij,ij->i", bary, h_old[old_tets[elem]])
+    h_q = h_q.reshape(len(tets), 4)
+    n = len(points)
+    rhs = np.zeros(n)
+    for a in range(4):                                       # integral of phi_a h: V/4 sum_q phi_a(q) h(q)
+        np.add.at(rhs, tets[:, a], vol / 4.0 * (QUADRATURE[:, a][None, :] * h_q).sum(axis=1))
+    local = np.full((4, 4), 1.0 / 20.0) + np.eye(4) / 20.0   # the P1 consistent mass on a tetrahedron, over its volume
+    rows, cols = np.repeat(tets, 4, axis=1).ravel(), np.tile(tets, (1, 4)).ravel()
+    M = sp.csr_matrix(((vol[:, None, None] * local[None]).ravel(), (rows, cols)), shape=(n, n))
+    lumped = np.asarray(M.sum(axis=1)).ravel()
+    h, info = spla.cg(M, rhs, x0=rhs / lumped, rtol=1e-13, maxiter=500, M=sp.diags(1.0 / lumped))
+    if info != 0:
+        raise RuntimeError("the rebuild's projection did not converge ({})".format(info))
+    return h
+
+
 def locate(points, tets, query, candidates=16):
     """For each query point, an element of (points, tets) holding it and its barycentric coordinates there; a point
     just outside every candidate gets the candidate it is least outside of (linear extrapolation)."""
```

- [ ] **Step 4: Run them to see them pass.** Expected: 34 passed.

- [ ] **Step 5: The recession measurement's switch** (`analysis/movemesh_recession.py`):

```diff
--- a/analysis/movemesh_recession.py
+++ b/analysis/movemesh_recession.py
@@ -8,7 +8,9 @@
 material, `AA7075_nomelt`) so the field the remap and the rebuilds carry has the flight's steep near-surface gradient.
 The run stops at NOSE_TOTAL of nose recession.
 
-    python analysis/movemesh_recession.py [--nose-total-mm 30] [--solver skfem|fenicsx] [--out PATH.json]
+    python analysis/movemesh_recession.py [--nose-total-mm 30] [--solver skfem|fenicsx] [--neutral] [--out PATH.json]
+
+`--neutral` uses the volume-neutral smoothing of sub-plan 20 Task 0 (`MovingMesh.neutral_smoothing`).
 """
 import argparse
 import json
@@ -31,6 +33,7 @@
 parser.add_argument("--no-smoothing", action="store_true", help="corners move along their normals only")
 parser.add_argument("--smooth-lambda", type=float, default=None)
 parser.add_argument("--smooth-clamp", type=float, default=None)
+parser.add_argument("--neutral", action="store_true", help="volume-neutral smoothing (sub-plan 20 Task 0)")
 args = parser.parse_args()
 
 mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_nomelt"])
@@ -47,7 +50,9 @@
     normal, _ = moving.frame()
     cos = np.clip(normal @ flight, 0.0, None)
     u = -NOSE_STEP * cos[:, None] * normal
-    if not args.no_smoothing:
+    if args.neutral:
+        u = u + moving.neutral_smoothing(normal)
+    elif not args.no_smoothing:
         u = u + moving.smoothing(normal, args.smooth_lambda, args.smooth_clamp)
     try:
         result = moving.move_with_policy(u)
@@ -82,6 +87,7 @@
            "max_abs_correction_over_swept": max(abs(r["correction_J"]) / r["swept_J"] for r in rows if r["swept_J"] > 0),
            "worst_ratio_min": min(r["worst_ratio"] for r in rows),
            "surface_quality_min": min(r["surface_quality"] for r in rows), "smoothing": not args.no_smoothing,
+           "neutral": args.neutral,
            "move_s_median": float(np.median([r["move_s"] for r in rows])),
            "step_s_median": float(np.median([r["step_s"] for r in rows]))}
 print(json.dumps(summary, indent=1))
```

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

```python
"""The corner columns (layered-skin design section 5.1, sub-plan 20 Task 1): dual cells, normals, the edge graph and its
cotangent conduction weights, and the corner surface the coupled run reads."""
import numpy as np
import pytest

from reentry_model import columns, movemesh


def sphere_columns(m):
    s = m.surface()
    corners = np.unique(s.faces)
    local = np.full(m.n_nodes, -1)
    local[corners] = np.arange(len(corners))
    faces = local[movemesh.orient_outward(m.points, s.faces, s.normals)]
    return columns.Columns(m.points[corners], faces), s


def flat_grid(n=8, h=1e-3, jitter=0.25, seed=2):
    """A jittered flat triangulation of the square [0, n h]^2 in z = 0 (open edges), normals +z."""
    rng = np.random.default_rng(seed)
    x, y = np.meshgrid(np.arange(n + 1) * h, np.arange(n + 1) * h, indexing="ij")
    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
    inner = (x.ravel() > 0) & (x.ravel() < n * h) & (y.ravel() > 0) & (y.ravel() < n * h)
    pts[inner, :2] += rng.uniform(-jitter, jitter, (inner.sum(), 2)) * h
    idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    faces = np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])
    return pts, faces, inner


def test_dual_areas_add_up_to_the_surface_and_normals_point_out(coarse_sphere_mesh):
    col, s = sphere_columns(coarse_sphere_mesh)
    assert col.n == len(np.unique(s.faces))
    assert col.area.sum() == pytest.approx(s.areas.sum(), rel=1e-12)
    radial = col.points / np.linalg.norm(col.points, axis=1)[:, None]
    # the area-weighted normal of a 4 mm triangulation of a 50 mm sphere follows the radius: 0.18 degrees off on
    # average and 4.7 at worst, where gmsh left an irregular corner (measured while planning)
    angle = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", col.normal, radial), -1, 1)))
    assert angle.mean() < 0.5 and angle.max() < 6.0
    assert np.all(col.degree >= 3)
    assert len(col.edges) == 3 * len(col.faces) // 2           # closed: every edge in exactly two triangles


def test_cotangent_weights_are_exact_for_a_linear_field_on_a_flat_surface():
    pts, faces, inner = flat_grid()
    col = columns.Columns(pts, faces)
    f = 2.0 + 300.0 * pts[:, 0] - 700.0 * pts[:, 1]
    assert np.abs(col.laplacian(f)[inner]).max() < 1e-12
    assert np.all(col.normal[:, 2] > 0.999999)


def test_on_an_equilateral_grid_the_weight_is_the_dual_face_over_the_edge():
    h = 1e-3
    rows = 6
    pts = np.array([[i * h + 0.5 * h * (j % 2), j * h * np.sqrt(3) / 2, 0.0] for j in range(rows) for i in range(rows)])
    faces = []
    for j in range(rows - 1):
        for i in range(rows - 1):
            a, b, c, d = j * rows + i, j * rows + i + 1, (j + 1) * rows + i, (j + 1) * rows + i + 1
            faces += [[a, b, c], [b, d, c]] if j % 2 == 0 else [[a, b, d], [a, d, c]]
    col = columns.Columns(pts, np.array(faces))
    interior = np.all(col._tri >= 0, axis=1)
    assert np.allclose(col.edge_length[interior], h)
    assert np.allclose(col.cot_weight[interior], 1.0 / np.sqrt(3.0))
    assert np.allclose(col.face_length[interior] / col.edge_length[interior], col.cot_weight[interior])


def test_the_cotangent_laplacian_of_a_sphere_sees_its_curvature(coarse_sphere_mesh):
    col, _ = sphere_columns(coarse_sphere_mesh)
    # sum_j w_ij (x_j - x_i) = -2 A_i n_i / R on a sphere of radius R; with barycentric areas it is right on average
    # (0.5 % high, with a scatter of 7.6 % from corner to corner, measured while planning)
    lap = col.laplacian(col.points)
    curv = -np.einsum("ij,ij->i", lap, col.normal) / (2.0 * col.area)
    assert np.mean(curv) == pytest.approx(1.0 / 0.05, rel=0.02)


def test_update_follows_moved_corners_and_the_top_surface_sits_above_the_base(coarse_sphere_mesh):
    col, _ = sphere_columns(coarse_sphere_mesh)
    area0 = col.area.sum()
    col.update(col.points * 0.9)                                  # a sphere 10 % smaller
    assert col.area.sum() == pytest.approx(0.81 * area0, rel=1e-12)
    top = columns.CornerSurface.from_points(col.top_points(np.full(col.n, 1e-3)), col.faces)
    r = np.linalg.norm(top.points, axis=1)
    lift = r - np.linalg.norm(col.points, axis=1)                 # 1 mm times the cosine of the normal's tilt
    assert np.all(lift <= 1e-3 + 1e-12) and np.all(lift > 0.99e-3)
    assert top.n_patches == col.n and top.area > col.area.sum()
    v = np.array([-1.0, 0.0, 0.0])
    assert top.angles_to(v)[top.patch_toward(v)] < np.radians(5.0)
    assert top.projected_area(v) == pytest.approx(np.pi * (0.045 + 1e-3) ** 2, rel=0.01)
```

- [ ] **Step 2: Run them to see them fail.** Run: `"$PY" -m pytest tests/test_reentry_model_columns.py -q`. Expected:
  `ImportError: cannot import name 'columns' from 'reentry_model'`.

- [ ] **Step 3: Create `reentry_model/columns.py`**

```python
"""The corner columns of the layered-skin design (docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md
section 5.1): one column on every corner of the 3D mesh's surface triangulation, owning the corner's dual cell.

The triangulation's connectivity is fixed for the whole flight, so everything here that depends only on it (the
edges, the two triangles on each edge, the corner opposite each edge in each triangle) is built once; `update` recomputes
the geometry when the corners move. The dual cell is barycentric: a third of each touching triangle, bounded by the
segments from each touching edge's midpoint to the triangles' centroids, so the areas add up to the surface's exactly.

Sideways conduction between neighbouring columns uses the cotangent weight of their shared edge, w = (cot a + cot b)/2
with a and b the angles opposite the edge in its two triangles: the flux of a piecewise-linear field through the two
dual segments across the edge is exactly w (T_j - T_i) per unit conductivity and depth, so the discrete operator is
exact for a field linear along the surface. A weight is clipped at zero where the opposite angles are obtuse enough
to make it negative: a negative conductance would break the explicit step's positivity, and it is rare (3 of 6 756
edges on the 4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
The surface gradient and curvature operators the flow needs are added in sub-plan 21, the first stage that uses them."""
from dataclasses import dataclass

import numpy as np

from .movemesh import corner_frame


class Columns:
    """Corner columns on a triangulation (points (nc, 3), faces (nt, 3) in corner numbering, oriented outward)."""

    def __init__(self, points, faces):
        self.faces = np.asarray(faces, dtype=np.int64)
        self.n = int(np.asarray(points).shape[0])
        e = np.vstack([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]])
        opposite = np.concatenate([self.faces[:, 2], self.faces[:, 0], self.faces[:, 1]])
        triangle = np.tile(np.arange(len(self.faces)), 3)
        key = np.sort(e, axis=1)
        self.edges, inverse = np.unique(key, axis=0, return_inverse=True)
        inverse = inverse.ravel()
        counts = np.bincount(inverse, minlength=len(self.edges))
        if np.any(counts > 2):
            raise ValueError("an edge is shared by more than two triangles: the surface is not a manifold")
        # each edge's (up to) two triangles and the corner opposite it in each; -1 on an open boundary
        order = np.argsort(inverse, kind="stable")
        first = np.searchsorted(inverse[order], np.arange(len(self.edges)))
        self._tri = np.full((len(self.edges), 2), -1, dtype=np.int64)
        self._opp = np.full((len(self.edges), 2), -1, dtype=np.int64)
        self._tri[:, 0], self._opp[:, 0] = triangle[order][first], opposite[order][first]
        two = counts == 2
        self._tri[two, 1], self._opp[two, 1] = triangle[order][first[two] + 1], opposite[order][first[two] + 1]
        self.degree = np.bincount(self.edges.ravel(), minlength=self.n)
        self.update(points)

    def update(self, points):
        """Recompute the geometry for corners at `points` (the connectivity is unchanged)."""
        p = np.asarray(points, dtype=float)
        self.points = p
        self.normal, self.area = corner_frame(p, self.faces, np.arange(self.n))
        self.centroid = p[self.faces].mean(axis=1)
        i, j = self.edges[:, 0], self.edges[:, 1]
        self.edge_length = np.linalg.norm(p[j] - p[i], axis=1)
        mid = 0.5 * (p[i] + p[j])
        cot = np.zeros(len(self.edges))
        face = np.zeros(len(self.edges))
        for s in range(2):
            ok = self._tri[:, s] >= 0
            k = self._opp[ok, s]
            u, v = p[i[ok]] - p[k], p[j[ok]] - p[k]
            cot[ok] += np.einsum("ij,ij->i", u, v) / np.linalg.norm(np.cross(u, v), axis=1)
            face[ok] += np.linalg.norm(mid[ok] - self.centroid[self._tri[ok, s]], axis=1)
        self.cot_weight = 0.5 * cot
        self.conduct_weight = np.maximum(self.cot_weight, 0.0)
        self.face_length = face

    def top_points(self, thickness):
        """The top surface: each corner moved out along its normal by its column's thickness [m]."""
        return self.points + np.asarray(thickness, dtype=float)[:, None] * self.normal

    def laplacian(self, f, weights=None):
        """sum over neighbours j of w_ij (f_j - f_i) per corner (f per corner, any trailing shape)."""
        w = self.cot_weight if weights is None else weights
        f = np.asarray(f, dtype=float)
        i, j = self.edges[:, 0], self.edges[:, 1]
        d = (f[j] - f[i]) * w.reshape((-1,) + (1,) * (f.ndim - 1))
        out = np.zeros_like(f)
        np.add.at(out, i, d)
        np.add.at(out, j, -d)
        return out


@dataclass
class CornerSurface:
    """The gas-facing surface of a corner-column body as the coupled run reads it: one 'patch' per corner (its top dual
    cell), with outward unit normals and areas, and the triangles for drawing it."""
    points: np.ndarray            # (nc, 3) the top surface's corners
    faces: np.ndarray             # (nt, 3) in corner numbering
    normals: np.ndarray           # (nc, 3)
    areas: np.ndarray             # (nc,) m2

    @classmethod
    def from_points(cls, points, faces):
        points, faces = np.asarray(points, dtype=float), np.asarray(faces, dtype=np.int64)
        normal, area = corner_frame(points, faces, np.arange(len(points)))
        return cls(points, faces, normal, area)

    @property
    def n_patches(self):
        return len(self.areas)

    @property
    def area(self):
        return float(self.areas.sum())

    def angles_to(self, v_hat):
        v = np.asarray(v_hat, dtype=float)
        return np.arccos(np.clip(self.normals @ (v / np.linalg.norm(v)), -1.0, 1.0))

    def patch_toward(self, direction):
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))

    def projected_area(self, v_hat):
        v = np.asarray(v_hat, dtype=float)
        return float((self.areas * np.clip(self.normals @ (v / np.linalg.norm(v)), 0.0, None)).sum())
```

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

```diff
--- a/reentry_model/skin1d.py
+++ b/reentry_model/skin1d.py
@@ -1,5 +1,7 @@
 """One-dimensional skins under the melting surface patches (Step 3 sub-plan 18; design:
-docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md, sections 5-7).
+docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md, sections 5-7). Renamed from `skin.py` in sub-plan 20, which
+gives that name to the layered skins of the corner columns; this stage-1 design failed its box test and is retired in
+sub-plan 22.
 
 A skin belongs to one surface patch and has its area: a slab of `n_cells` equal cells along the patch's inward normal
 that owns the outermost metal of that patch. Its state is each cell's metal mass `m` [kg] and metal enthalpy `H`
```

```diff
--- a/tests/test_reentry_model_skin1d.py
+++ b/tests/test_reentry_model_skin1d.py
@@ -1,4 +1,4 @@
-"""reentry_model.skin: the one-dimensional skins of the 2026-10-07 design (Step 3 sub-plan 18) -- storage and rezoning,
+"""reentry_model.skin1d (stage 1, renamed from skin in sub-plan 20): the one-dimensional skins of the 2026-10-07 design (Step 3 sub-plan 18) -- storage and rezoning,
 the top node, depths, hand-over, the implicit sub-step against exact solutions, feed and draw, and the macro step's
 trial and real passes."""
 import math
@@ -9,7 +9,8 @@
 from scipy.optimize import brentq
 from scipy.special import erf, erfc
 
-from reentry_model import material, skin
+from reentry_model import material
+from reentry_model import skin1d as skin
 
 
 def constant_material(latent=0.0, T_melt=float("inf")):
```

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -282,7 +282,7 @@
         self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation
         self.skins = None                                          # sub-plan 18: the skins, with --surface-model skin
         if self.settings.surface_model == "skin" and self.settings.removal != "instant":
-            from . import skin as skin_mod
+            from . import skin1d as skin_mod
             self.skins = skin_mod.SkinField(material, self.settings.skin or skin_mod.SkinSettings())
         self._draw_m, self._draw_H = np.zeros(mesh.n_elements), np.zeros(mesh.n_elements)   # metal the skins drew this step
         self.skins_created = self.skin_short_events = self.interface_repeats = 0
@@ -422,7 +422,7 @@
         the rest of the melt step. Until the first skin exists this is the element model's step exactly."""
         import dataclasses
         import time
-        from . import skin as skin_mod
+        from . import skin1d as skin_mod
         from . import spray as spray_mod
         mat, liq, sk = self.material, self.liquid, self.skins
         self._create_skins(self.solver.temperature())
@@ -1190,7 +1190,7 @@
         the difference is booked to the element's nodes, as with the draw. Where the element is full the skin keeps the
         metal and the event is counted. It moves metal between two accounts of the same solid and changes neither how
         much freezes nor how much melts."""
-        from .skin import MAX_THICKNESS_FACTOR
+        from .skin1d import MAX_THICKNESS_FACTOR
         sk, s = self.skins, self.surface
         p = np.flatnonzero(self.skin_of_patch >= 0)
         r, e = self.skin_of_patch[p], s.owner[p]
@@ -1632,7 +1632,7 @@
     def _skin_stats(self):
         """The skins' history columns (sub-plan 18); zero or nan without skins. The resolved depths are the step's, taken
         in `_film_and_spray` on the surface they were evaluated on."""
-        from .skin import MAX_THICKNESS_FACTOR
+        from .skin1d import MAX_THICKNESS_FACTOR
         sk, lm, nan = self.skins, self.last_melt, float("nan")
         th = sk.thickness() if sk is not None and sk.n else np.zeros(0)
         return {"skin_count": float(sk.n) if sk is not None else 0.0, "skin_mass_kg": sk.mass() if sk is not None else 0.0,
```

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -19,7 +19,7 @@
 
 import numpy as np
 
-from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, skin, spray, surface_flow, thermal, viz
+from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, skin1d as skin, spray, surface_flow, thermal, viz
 from . import trajectory as tj
 from .earth import GRAVITY_MODELS
 from .thermal import MissingBackend
```

```diff
--- a/tests/test_reentry_model_melting_skin.py
+++ b/tests/test_reentry_model_melting_skin.py
@@ -4,7 +4,8 @@
 import numpy as np
 import pytest
 
-from reentry_model import body, heating, material, mesh, skin, thermal
+from reentry_model import body, heating, material, mesh, thermal
+from reentry_model import skin1d as skin
 from test_reentry_model_coupled import MASS_100MM, simulator
 
 pytest.importorskip("cantera")
```

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -313,7 +313,8 @@
     of 0.5 s, every patch skinned at the first step. Skins, mass, film and temperatures agree and both energy balances
     are exact."""
     pytest.importorskip("cantera")
-    from reentry_model import body, heating, mesh, skin
+    from reentry_model import body, heating, mesh
+    from reentry_model import skin1d as skin
     from test_reentry_model_coupled import MASS_100MM
     out = {}
     for name in ("skfem", "fenicsx"):
```

  Run the unit tier: expected 349 passed, 1 skipped (stage 1's 21 tests now in `test_reentry_model_skin1d.py`).

- [ ] **Step 2: Write the new skin's tests** (`tests/test_reentry_model_skin.py`, which grows by a section per task):

```python
"""The layered skins (layered-skin design section 6, sub-plan 20): cells with variable counts, split and merge, the
minimum mass, heat through the depth with liquid removed as it forms, and heat sideways between columns."""
import numpy as np
import pytest

from reentry_model import material as material_mod
from reentry_model import skin


def scheil():
    return material_mod.Material.from_drama_json(material_mod.MATERIAL_NAMES["AA7075_scheil"])


def sharp(k=100.0, rho=2800.0, cp=1000.0, latent=4.0e5, T_melt=900.0):
    """Constant properties, a melting point (with the material's +-2 K ramp) and liquid as dense as solid: the
    material of the steady-ablation solution."""
    return material_mod.Material("sharp", rho, 0.0, np.array([0.0, 5000.0]), np.array([cp, cp]), np.array([0.0, 5000.0]),
                                 np.array([k, k]), latent, T_melt, T_melt, material_mod.LiquidProperties(rho, 1e-3, 0.8))


def filled(mat, n=3, area=4e-6, T=700.0, settings=None):
    s = settings or skin.SkinSettings()
    sk = skin.Skins(mat, np.full(n, area), s)
    sk.fill(skin.minimum_mass(mat, sk.area, s), T)
    return sk


# -- settings and cells (Task 2) ---------------------------------------------------------------------------------------
def test_settings_reject_impossible_values():
    for bad in ({"min_thickness": 0.0}, {"cell": -1e-6}, {"substep": 0.0}, {"growth": 0.9}, {"cell": 0.3e-3}):
        with pytest.raises(ValueError):
            skin.SkinSettings(**bad)


def test_a_filled_skin_holds_its_mass_in_fine_cells():
    mat = scheil()
    sk = filled(mat)
    assert np.all(sk.count == 40)                                   # 0.4 mm of 10 um cells
    assert sk.thickness() == pytest.approx(np.full(3, 0.4e-3), rel=1e-12)
    assert sk.mass() == pytest.approx(3 * mat.rho * 4e-6 * 0.4e-3, rel=1e-12)
    assert sk.energy() == pytest.approx(sk.mass() * float(mat.enthalpy(700.0)), rel=1e-12)
    assert np.allclose(sk.T[sk.valid], 700.0)


def test_split_and_merge_conserve_mass_and_enthalpy_and_keep_cells_near_their_targets():
    mat = scheil()
    sk = filled(mat, n=4)
    rng = np.random.default_rng(1)
    v = sk.valid
    sk.m = np.where(v, sk.m * rng.uniform(0.1, 3.0, sk.m.shape), 0.0)       # cells of every size
    sk.H = sk.m * mat.enthalpy(rng.uniform(600.0, 890.0, sk.m.shape))
    sk._install(sk.m, sk.H, sk.count)
    M0, E0 = sk.mass(), sk.energy()
    sk.regrid()
    assert sk.mass() == pytest.approx(M0, rel=1e-14) and sk.energy() == pytest.approx(E0, rel=1e-14)
    dz, t = sk.dz(), sk.targets()
    v = sk.valid
    assert np.all(dz[v] >= skin.MERGE_BELOW * t[v]) and np.all(dz[v] <= skin.SPLIT_ABOVE * t[v])


def test_the_top_cells_are_untouched_by_a_split_or_merge_deeper_down():
    mat = scheil()
    sk = filled(mat, n=2)
    top_m, top_H = sk.m[:, :10].copy(), sk.H[:, :10].copy()
    sk.add_bottom(np.array([1e-7, 0.0]), float(mat.enthalpy(700.0)))       # a bottom cell 1.9 x its target: split
    sk.take_bottom(np.array([0.0, 0.7 * sk.m[1, 0]]))                       # a thin bottom cell: merged
    sk.regrid()
    assert np.array_equal(sk.m[:, :10], top_m) and np.array_equal(sk.H[:, :10], top_H)
    assert sk.count[0] > 40 and sk.count[1] == 39


def test_cells_below_the_fine_zone_grow_by_the_growth_factor():
    mat = scheil()
    sk = filled(mat, n=1)
    sk.add_bottom(np.array([mat.rho * 4e-6 * 2.6e-3]), float(mat.enthalpy(700.0)))   # 3 mm in all
    sk.regrid()
    dz = sk.dz()[0, :sk.count[0]]
    assert sk.thickness()[0] == pytest.approx(3.0e-3, rel=1e-12)
    # the old bottom cell took the 2.6 mm and was halved until every piece is within its bounds: 39 fine cells above
    # it, then pieces doubling in steps as the target grows (75 cells in all, measured while planning; a continuous
    # 10 % growth would need 40 + 10 ln(1 + 0.1 x 2600 / 10) = 73)
    assert np.allclose(dz[:39], 10e-6)
    assert 60 <= sk.count[0] <= 90
    assert np.all(np.diff(dz) >= -1e-15)                             # never a thick cell above a thinner one


def test_taking_and_exchanging_at_the_bottom_and_top_are_exact():
    mat = scheil()
    sk = filled(mat)
    M0, E0 = sk.mass(), sk.energy()
    dm = np.array([0.0, 0.25, 1.5]) * sk.m[0, 0]
    H = sk.take_bottom(dm)
    assert sk.mass() == pytest.approx(M0 - dm.sum(), rel=1e-14) and sk.energy() == pytest.approx(E0 - H.sum(), rel=1e-14)
    sk.exchange_bottom(np.array([1e-9, -2e-9, 0.0]), np.array([1e-3, -3e-3, 5e-4]))
    assert sk.mass() == pytest.approx(M0 - dm.sum() - 1e-9, rel=1e-14)
    assert sk.energy() == pytest.approx(E0 - H.sum() - 1.5e-3, rel=1e-14)
    m_top, H_top = sk.take_top(np.array([2, 0, 1]))
    assert np.all(sk.count == np.array([38, 40, 38]))                  # the third had lost a cell and a half
    assert m_top[1] == 0.0 and m_top[0] == pytest.approx(2 * sk.m[0, 0], rel=1e-12)
    with pytest.raises(ValueError):
        sk.take_bottom(np.full(3, 1.0))


def test_depths_read_from_the_cells():
    mat = sharp()                                                    # liquid as dense as solid: the cells keep 10 um
    sk = filled(mat, n=1)
    z_mid = np.arange(40) * 10e-6 + 5e-6
    sk.T[0, :40] = 920.0 - 1e5 * z_mid                               # a straight profile, 10 K per 100 um
    sk.H = sk.m * mat.enthalpy(sk.T)
    assert sk.depth_above(900.0)[0] == pytest.approx(200e-6, abs=1e-9)
    assert sk.depth_above(1000.0)[0] == 0.0 and sk.depth_above(0.0)[0] == pytest.approx(0.4e-3)


```

- [ ] **Step 3: Run them to see them fail.** Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `ImportError: cannot import name 'skin' from 'reentry_model'`.

- [ ] **Step 4: Create `reentry_model/skin.py`**

```python
"""The layered skins on the corner columns (layered-skin design, docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-
design.md sections 6 and 7; stage 2, sub-plan 20: heat only, liquid removed as it forms).

Every corner column carries a skin from time zero: a stack of cells along the column's direction holding all metal
outward of the rigid front, at least the column's minimum mass. A cell holds its mass [kg] and enthalpy [J above
material.T_REF]; its temperature, liquid fraction and density follow from the material's enthalpy law, and its thickness
is its mass over the column's area and its own density (liquid is about 15 % less dense than solid). All columns are
stored side by side in padded arrays (n, K): cell 0 is the top, cells 0 .. count - 1 are real, the rest hold nothing.

Cells move with the metal: nothing re-divides a whole column. After every sub-step's exchanges a cell more than
SPLIT_ABOVE times its target thickness is split in two and one less than MERGE_BELOW times it is merged into the cell
beneath (the bottom cell into the one above), exactly in mass and enthalpy. The target is the fine cell within the fine
zone, which is the minimum thickness below the top, and grows by (growth - 1) per unit depth below it -- the continuous
form of each cell being `growth` times thicker than the one above.

Stage 2 has no flow and no spraying: after each sub-step the contiguous cells at the top that are at or above T_feed
leave the body with their enthalpy (`removed`), and a column below its minimum draws rigid metal from the 3D model at
the base temperature (`drawn`). The body (`skinbody.SkinBody`) turns the drawn, absorbed and handed-down metal into the
corners' displacements and reconciles each column with what the moving mesh actually swept."""
import math
from dataclasses import dataclass

import numpy as np

from .thermal import SIGMA_SB

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (spec section 6.2)
SPLIT_ABOVE = 1.5                 # a cell thicker than this times its target is split in two (spec 6.2)
MERGE_BELOW = 0.5                 # a cell thinner than this times its target is merged into its neighbour (spec 6.2)
REGRID_PASSES = 64                # merge-and-split rounds after an exchange: one unless a cell grew many-fold, each
                                  # round halving it (2.6 mm added at once needs 9)
DELTA_M_MARGIN = 50.0e-6          # m of liquid kept below Girin's conjugate depth in the minimum mass (spec 6.1)
INTERFACE_TOLERANCE = 1.0e-3      # |mismatch| / (|E_top| + |E_base|) above which a macro step repeats (spec 7.2)
JUMP_TOLERANCE = 5.0              # K between the metal crossing the base and the skin's bottom cell: a fault (spec 6.6,
                                  # provisional)
NEWTON_TOL = 1.0e-6               # K: a sub-step has converged when no cell moves by more than this
NEWTON_MAX = 40
LINE_SEARCH_MAX = 12              # halvings of a Newton step before it is taken anyway
DZ_FLOOR = 1.0e-12                # m: the conduction length of a vanishing cell
SIDEWAYS_CFL = 0.9                # explicit sideways conduction: sub-cycles keep dt k sum(w) / (rho c_p A) below this
TABLE_DT = 0.01                   # K, spacing of the skins' property tables (PropertyTable)
TABLE_RANGE = (150.0, 3000.0)     # K, their extent; beyond it the end intervals are extended linearly


# -- settings and cells (Task 2) ---------------------------------------------------------------------------------------
@dataclass
class SkinSettings:
    min_thickness: float = 0.4e-3          # m of solid metal: the minimum mass, and the fine zone's depth (spec 6.1)
    cell: float = SKIN_PRESETS["production"]   # m, the fine cell
    growth: float = 1.10                   # each cell below the fine zone this much thicker than the one above
    substep: float = 0.010                 # s, the skins' sub-step (spec 7.1)

    def __post_init__(self):
        if not (self.min_thickness > 0.0 and self.cell > 0.0 and self.substep > 0.0):
            raise ValueError("skin minimum thickness, cell and substep must be positive")
        if self.growth < 1.0:
            raise ValueError("skin growth must be at least 1")
        if 2.0 * self.cell > self.min_thickness * (1.0 + 1e-12):
            raise ValueError("a skin needs at least two fine cells: cell <= min_thickness / 2")


def thomas(sub, diag, sup, rhs):
    """Row-wise tridiagonal solve: sub and sup (n, k-1), diag and rhs (n, k). Returns x (n, k)."""
    n, k = diag.shape
    c = np.zeros((n, max(k - 1, 1)))
    d = np.zeros((n, k))
    if k > 1:
        c[:, 0] = sup[:, 0] / diag[:, 0]
    d[:, 0] = rhs[:, 0] / diag[:, 0]
    for i in range(1, k):
        denom = diag[:, i] - sub[:, i - 1] * c[:, i - 1]
        if i < k - 1:
            c[:, i] = sup[:, i] / denom
        d[:, i] = (rhs[:, i] - sub[:, i - 1] * d[:, i - 1]) / denom
    x = np.zeros((n, k))
    x[:, -1] = d[:, -1]
    for i in range(k - 2, -1, -1):
        x[:, i] = d[:, i] - c[:, i] * x[:, i + 1]
    return x


class PropertyTable:
    """The material's h(T), c_p,eff(T), k(T) and f_l(T) tabulated on a uniform TABLE_DT grid and interpolated
    linearly: a lookup is index arithmetic instead of the material's binary search, which was half of a skin pass's
    time (measured while planning: 11 s of a 4 mm sphere's macro step, 5.5 s of it in the property calls). The
    interpolation moves h by at most about 60 J/kg (0.05 K) across a Scheil node, and the books never see it: every pass
    starts from each cell's stored enthalpy."""

    def __init__(self, material, dT=TABLE_DT, T_range=TABLE_RANGE):
        self.T0, self.dT = float(T_range[0]), float(dT)
        T = np.arange(T_range[0], T_range[1] + dT, dT)
        self.n = len(T)
        self._h, self._cp, self._k, self._f = (material.enthalpy(T), material.cp_eff(T), material.k(T),
                                               material.liquid_fraction(T))

    def _at(self, table, T):
        x = (np.asarray(T, dtype=float) - self.T0) / self.dT
        i = np.clip(x.astype(np.int64), 0, self.n - 2)
        return table[i] + (x - i) * (table[i + 1] - table[i])

    def conduction(self, T):
        """(h, dh/dT, k) at T with one index computation, for the conduction's Newton iteration. The heat capacity is
        the slope of the tabulated h itself, so the Newton tangent is exact: interpolating c_p,eff instead left it off
        the slope near the Scheil curve's 1 K nodes and took five iterations a sub-step where three do (measured while
        planning)."""
        x = (np.asarray(T, dtype=float) - self.T0) / self.dT
        i = np.clip(x.astype(np.int64), 0, self.n - 2)
        dh = self._h[i + 1] - self._h[i]
        return self._h[i] + (x - i) * dh, dh / self.dT, self._k[i] + (x - i) * (self._k[i + 1] - self._k[i])

    def enthalpy(self, T):
        return self._at(self._h, T)

    def cp_eff(self, T):
        return self._at(self._cp, T)

    def k(self, T):
        return self._at(self._k, T)

    def liquid_fraction(self, T):
        return np.clip(self._at(self._f, T), 0.0, 1.0)


class Skins:
    """The skins of every column, as padded arrays (module docstring). `area` [m2] is each column's area, the dual
    area of its base corner; the body updates it when the corners move (the cells keep their mass). `S`, when not None,
    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b)."""

    def __init__(self, material, area, settings=None):
        self.material = material
        self.settings = settings or SkinSettings()
        self.area = np.array(area, dtype=float)
        n = len(self.area)
        self.m, self.H, self.T = np.zeros((n, 1)), np.zeros((n, 1)), np.full((n, 1), 300.0)
        self.count = np.zeros(n, dtype=np.int64)
        self.S = None
        self.rho_l = material.liquid.rho if material.liquid is not None else material.rho
        self.props = PropertyTable(material)

    # -- books ------------------------------------------------------------------------------------------------------
    @property
    def n(self):
        return len(self.area)

    @property
    def valid(self):
        return np.arange(self.m.shape[1])[None, :] < self.count[:, None]

    def mass(self):
        return float(self.m.sum())

    def energy(self):
        return float(self.H.sum())

    def column_mass(self):
        return self.m.sum(axis=1)

    def state(self):
        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count))

    def set_state(self, state):
        self.area, self.m, self.H, self.T, self.count = (a.copy() for a in state)
        self.S = None

    # -- geometry of the cells --------------------------------------------------------------------------------------
    def density(self, T):
        """Mixture density [kg/m3]: specific volumes add, (1 - f_l)/rho_s + f_l/rho_l."""
        f = self.props.liquid_fraction(T)
        return 1.0 / ((1.0 - f) / self.material.rho + f / self.rho_l)

    def dz(self):
        """Cell thicknesses [m] (0 beyond each column's count)."""
        return np.where(self.valid, self.m / (self.density(self.T) * self.area[:, None]), 0.0)

    def thickness(self):
        return self.dz().sum(axis=1)

    def targets(self, dz=None):
        """Target thickness of every cell [m] from the depth of its top below the column's top (module docstring)."""
        dz = self.dz() if dz is None else dz
        z = np.cumsum(dz, axis=1) - dz
        s = self.settings
        return s.cell + (s.growth - 1.0) * np.maximum(z - s.min_thickness, 0.0)

    # -- installing cells -------------------------------------------------------------------------------------------
    def _install(self, m, H, count, S=None, T=None):
        """Install padded cell arrays. Temperatures follow from the enthalpies, except where `T` is given for cells
        whose specific enthalpy the operation did not change (a shifted, split or partly taken cell); nan in T marks
        the cells to recompute."""
        valid = np.arange(m.shape[1])[None, :] < count[:, None]
        self.m, self.H, self.count = np.where(valid, m, 0.0), np.where(valid, H, 0.0), np.asarray(count, dtype=np.int64)
        redo = valid if T is None else valid & np.isnan(T)
        self.T = np.full(self.m.shape, 300.0) if T is None else np.where(valid & ~redo, T, 300.0)
        if redo.any():
            self.T[redo] = self.material.temperature_from_enthalpy(self.H[redo] / self.m[redo])
        if self.S is not None:
            self.S = np.where(valid, S, 0.0)

    def fill(self, mass, T):
        """Every column filled with `mass` [kg] (per column) of metal at temperature T [K] (per column), in equal cells of
        about the fine size (the initial skins, spec 6.1)."""
        mass, T = np.broadcast_to(np.asarray(mass, dtype=float), (self.n,)), np.broadcast_to(np.asarray(T, dtype=float), (self.n,))
        thick = mass / (self.density(T) * self.area)
        nc = np.maximum(np.rint(thick / self.settings.cell).astype(np.int64), 1)
        K = int(nc.max())
        valid = np.arange(K)[None, :] < nc[:, None]
        m = np.where(valid, (mass / nc)[:, None], 0.0)
        self.S = None
        self._install(m, m * self.material.enthalpy(T)[:, None], nc)

    def _width(self, K):
        """Pad the arrays to at least K cells."""
        extra = K - self.m.shape[1]
        if extra > 0:
            pad = lambda a, v: np.hstack([a, np.full((self.n, extra), v)])
            self.m, self.H, self.T = pad(self.m, 0.0), pad(self.H, 0.0), pad(self.T, 300.0)
            if self.S is not None:
                self.S = pad(self.S, 0.0)

    # -- split and merge --------------------------------------------------------------------------------------------
    def regrid(self):
        """Merge cells below MERGE_BELOW of their target, split those above SPLIT_ABOVE, until none qualifies. Only the
        columns that need it are rebuilt (a split every sub-step on every melting column is the common case)."""
        for _ in range(REGRID_PASSES):
            valid = self.valid
            dz = self.dz()
            t = self.targets(dz)
            small = valid & (dz < MERGE_BELOW * t) & (self.count >= 2)[:, None]
            big = valid & (dz > SPLIT_ABOVE * t)
            if not (small.any() or big.any()):
                return
            if small.any():
                self._merge(np.flatnonzero(small.any(axis=1)), small)
            else:
                self._split(np.flatnonzero(big.any(axis=1)), big)
        raise RuntimeError("the skin cells did not settle in {} split-and-merge passes".format(REGRID_PASSES))

    def _put_rows(self, rows, m, H, T, S, count):
        """Write rebuilt cell arrays for `rows` back (temperatures of nan are recomputed from the enthalpies), padding
        or trimming the arrays to the longest column."""
        self._width(m.shape[1])
        K = self.m.shape[1]
        pad = lambda a: np.hstack([a, np.zeros((len(rows), K - a.shape[1]))]) if a.shape[1] < K else a
        valid = np.arange(K)[None, :] < count[:, None]
        m, H, T = pad(m), pad(H), pad(T)
        redo = valid & np.isnan(T)
        if redo.any():
            T[redo] = self.material.temperature_from_enthalpy(H[redo] / m[redo])
        self.m[rows], self.H[rows] = np.where(valid, m, 0.0), np.where(valid, H, 0.0)
        self.T[rows] = np.where(valid, T, 300.0)
        self.count[rows] = count
        if self.S is not None:
            self.S[rows] = np.where(valid, pad(S), 0.0)
        width = int(self.count.max()) if self.n else 1
        if width < K:
            self.m, self.H, self.T = self.m[:, :width], self.H[:, :width], self.T[:, :width]
            if self.S is not None:
                self.S = self.S[:, :width]

    def _merge(self, rows, small):
        """Each run of small cells joins the cell beneath that ends it; a run at the bottom with no cell to end it joins
        the group above (in `rows` only)."""
        small = small[rows]
        valid = self.valid[rows]
        m0, H0, T0 = self.m[rows], self.H[rows], self.T[rows]
        n, K = m0.shape
        start = valid.copy()
        start[:, 1:] &= ~small[:, :-1]
        gid = np.cumsum(start, axis=1) - 1
        groups = start.sum(axis=1)
        r = np.arange(n)
        trailing = small[r, np.maximum(self.count[rows] - 1, 0)] & (groups >= 2)
        gid = np.where(trailing[:, None] & (gid == (groups - 1)[:, None]), gid - 1, gid)
        groups = groups - trailing
        flat = (r[:, None] * K + gid)[valid]
        m, H = np.zeros(n * K), np.zeros(n * K)
        np.add.at(m, flat, m0[valid])
        np.add.at(H, flat, H0[valid])
        size = np.bincount(flat, minlength=n * K)
        T = np.full(n * K, np.nan)
        T[flat] = np.where(size[flat] == 1, T0[valid], np.nan)             # a cell that merged with none keeps its T
        S = None
        if self.S is not None:
            mS = np.zeros(n * K)
            np.add.at(mS, flat, (m0 * self.S[rows])[valid])
            S = np.divide(mS, m, out=np.zeros_like(m), where=m > 0.0).reshape(n, K)
        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups)

    def _split(self, rows, big):
        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T and S)."""
        valid = self.valid[rows]
        reps = np.where(valid, 1 + big[rows], 0)
        n, K = reps.shape
        count = reps.sum(axis=1)
        K2 = int(count.max())
        pos = np.cumsum(reps, axis=1) - reps
        r_all = np.repeat(np.arange(n)[:, None], K, axis=1)
        m, H, T = np.zeros((n, K2)), np.zeros((n, K2)), np.zeros((n, K2))
        S = np.zeros((n, K2)) if self.S is not None else None
        share = 1.0 / reps.clip(1)
        for part in (0, 1):
            sel = valid & (part < reps)
            r, c = r_all[sel], pos[sel] + part
            m[r, c] = (self.m[rows] * share)[sel]
            H[r, c] = (self.H[rows] * share)[sel]
            T[r, c] = self.T[rows][sel]                                    # halves keep the specific enthalpy
            if S is not None:
                S[r, c] = self.S[rows][sel]
        self._put_rows(rows, m, H, T, S, count)

    # -- the top and the bottom -------------------------------------------------------------------------------------
    def take_top(self, k):
        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column."""
        k = np.asarray(k, dtype=np.int64)
        n, K = self.m.shape
        top = np.arange(K)[None, :] < k[:, None]
        mass, H = (self.m * top).sum(axis=1), (self.H * top).sum(axis=1)
        if k.any():
            src = np.minimum(np.arange(K)[None, :] + k[:, None], K - 1)
            shift = lambda a: np.take_along_axis(a, src, axis=1)
            S = shift(self.S) if self.S is not None else None
            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T))
        return mass, H

    def add_bottom(self, dm, h):
        """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
        a new cell); returns the enthalpy added. With S tracked, the metal carries S = 1 (it is at the base temperature)."""
        dm, h = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,)), np.broadcast_to(np.asarray(h, dtype=float), (self.n,))
        some = dm > 0.0
        if not some.any():
            return np.zeros(self.n)
        self._width(int(np.max(np.maximum(self.count, 1))))
        rows = np.flatnonzero(some)
        b = np.maximum(self.count[rows] - 1, 0)
        m_old = self.m[rows, b]
        self.m[rows, b] += dm[rows]
        self.H[rows, b] += dm[rows] * h[rows]
        if self.S is not None:
            self.S[rows, b] = (m_old * self.S[rows, b] + dm[rows]) / self.m[rows, b]
        self.count[rows] = np.maximum(self.count[rows], 1)
        self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])
        return np.where(some, dm * h, 0.0)

    def take_bottom(self, dm):
        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy (a partly
        taken cell keeps its specific enthalpy); returns the enthalpy removed per column."""
        dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
        if np.any(dm > self.column_mass() * (1.0 + 1e-12) + 1e-300):
            raise ValueError("cannot take more metal from a skin's bottom than the skin holds")
        if not np.any(dm > 0.0):
            return np.zeros(self.n)
        n, K = self.m.shape
        valid = self.valid
        below = np.cumsum(self.m[:, ::-1], axis=1)[:, ::-1] - self.m      # mass beneath each cell
        taken = np.clip(dm[:, None] - below, 0.0, self.m)                 # how much of each cell goes
        frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
        H_taken = (self.H * frac).sum(axis=1)
        m, H = self.m - taken, self.H * (1.0 - frac)
        keep = valid & (m > 0.0)
        count = keep.sum(axis=1)                                          # the kept cells are the top ones
        self._install(m, H, count, self.S, self.T)                        # a partly taken cell keeps its T
        return H_taken

    def exchange_bottom(self, dm, dH):
        """Signed reconciliation at the bottom of every column: the skin gains dm [kg] (negative: loses) and its
        enthalpy changes by exactly dH [J]. Gained metal joins the bottom cell; lost metal leaves from the bottom up at
        the cells' own specific enthalpy, and the difference to dH stays in the new bottom cell."""
        dm, dH = np.asarray(dm, dtype=float), np.asarray(dH, dtype=float)
        gain = np.where(dm > 0.0, dm, 0.0)
        h = np.divide(dH, gain, out=np.zeros_like(dH), where=gain > 0.0)
        self.add_bottom(gain, h)
        rest = np.where(dm > 0.0, 0.0, dH)
        H_taken = self.take_bottom(np.where(dm < 0.0, -dm, 0.0))
        rest = rest + H_taken                                              # what is left over after the cells left
        rows = np.flatnonzero(rest != 0.0)
        if rows.size:
            b = self.count[rows] - 1
            if np.any(b < 0):
                raise ValueError("a skin emptied by a reconciliation cannot take the enthalpy difference")
            self.H[rows, b] += rest[rows]
            self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])

    def depth_above(self, T_level):
        """Per column: the depth [m] from the top to where the contiguous cells above T_level end, interpolated between
        cell centres (0 if the top cell is below it, the whole skin if every cell is above)."""
        valid = self.valid
        above = np.cumprod(valid & (self.T > T_level), axis=1).astype(bool)
        k = above.sum(axis=1)
        dz = self.dz()
        z_mid = np.cumsum(dz, axis=1) - 0.5 * dz
        rows = np.arange(self.n)
        Ka = np.clip(k - 1, 0, self.m.shape[1] - 1)
        Kb = np.clip(k, 0, self.m.shape[1] - 1)
        Ta, Tb = self.T[rows, Ka], self.T[rows, Kb]
        frac = np.clip(np.divide(Ta - T_level, Ta - Tb, out=np.zeros_like(Ta), where=Ta > Tb), 0.0, 1.0)
        inside = z_mid[rows, Ka] + frac * (z_mid[rows, Kb] - z_mid[rows, Ka])
        return np.where(k == 0, 0.0, np.where(k >= self.count, self.thickness(), inside))


```

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

```python
# -- the minimum mass (Task 3) -----------------------------------------------------------------------------------------
def test_the_minimum_is_400_um_of_solid_unless_girins_layer_needs_more():
    mat = scheil()
    s = skin.SkinSettings()
    area = np.full(4, 4e-6)
    solid = mat.rho * 4e-6 * 0.4e-3
    # 400 um of solid is the liquid of 469 um (2813 / 2400 x 400): delta_m + 50 um decides above 419 um
    m = skin.minimum_mass(mat, area, s, np.array([np.nan, 300e-6, 419e-6, 500e-6]))
    assert m[0] == solid and m[1] == solid
    assert m[2] == pytest.approx(solid, rel=2e-3)
    assert m[3] == pytest.approx(mat.liquid.rho * 4e-6 * 550e-6, rel=1e-12) and m[3] > solid
    assert np.all(skin.minimum_mass(mat, area, s) == solid)


```

- [ ] **Step 2: Run it to see it fail.** Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'minimum_mass'`.

- [ ] **Step 3: Append to `reentry_model/skin.py`**

```python
# -- the minimum mass (Task 3) -----------------------------------------------------------------------------------------
def minimum_mass(material, area, settings, delta_m=None):
    """Each column's minimum mass [kg] (spec 6.1): the larger of `min_thickness` of solid metal over the column's area
    and the liquid filling Girin's conjugate depth delta_m [m] plus DELTA_M_MARGIN, so that the sheared layer always lies
    inside the fine cells. delta_m None (or nan where Girin's closure does not apply) leaves the solid term alone."""
    area = np.asarray(area, dtype=float)
    solid = material.rho * area * settings.min_thickness
    if delta_m is None:
        return solid
    rho_l = material.liquid.rho if material.liquid is not None else material.rho
    liquid = rho_l * area * (np.nan_to_num(np.asarray(delta_m, dtype=float), nan=-DELTA_M_MARGIN) + DELTA_M_MARGIN)
    return np.maximum(solid, liquid)


```

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

```python
# -- heat through the depth (Task 4) -----------------------------------------------------------------------------------
def test_a_skin_at_its_minimum_above_it_and_drained_returns_to_it_by_drawing():
    mat = scheil()
    s = skin.SkinSettings(substep=0.01)
    sk = filled(mat, n=3, settings=s)
    m_min = skin.minimum_mass(mat, sk.area, s)
    sk.add_bottom(np.array([0.0, 0.5, 0.0]) * m_min, float(mat.enthalpy(700.0)))   # one above its minimum
    sk.take_bottom(np.array([0.0, 0.0, 0.6]) * m_min)                              # one drained below it
    sk.regrid()
    book = skin.advance(sk, 0.05, np.zeros(3), sk.area, 0.0, 0.0, np.full(3, 700.0), m_min, T_b1=np.full(3, 700.0))
    assert book.drawn == pytest.approx(np.array([0.0, 0.0, 0.6]) * m_min, rel=1e-12)
    assert sk.column_mass() == pytest.approx(np.array([1.0, 1.5, 1.0]) * m_min, rel=1e-12)
    assert np.allclose(sk.T[sk.valid], 700.0)


def test_steady_ablation_recedes_at_the_energy_limited_speed_with_the_exponential_profile():
    """Liquid removed as it forms under a constant flux: in steady state the top recedes at v = q / (rho (h(T_feed) -
    h(T0))) and the temperature beneath it is T0 + (T_feed - T0) exp(-v z / alpha)."""
    mat = sharp()
    q, T0 = 1.0e8, 300.0
    s = skin.SkinSettings(min_thickness=5e-3, cell=20e-6, substep=5e-4)     # 5 decay lengths, base held at T0
    sk = skin.Skins(mat, np.array([1e-6]), s)
    m_min = skin.minimum_mass(mat, sk.area, s)
    sk.fill(m_min, T0)
    v = q / (mat.rho * float(mat.enthalpy(mat.T_feed) - mat.enthalpy(T0)))       # 35.6 mm/s
    alpha = 100.0 / (2800.0 * 1000.0)
    removed = []
    for _ in range(30):
        book = skin.advance(sk, 0.01, np.array([q]), sk.area, 0.0, 0.0, np.array([T0]), m_min, T_b1=np.array([T0]))
        removed.append(book.removed[0])
    speed = np.sum(removed[10:]) / (mat.rho * 1e-6) / 0.2                         # after 0.1 s, 3.5 rise times
    # 2.1 % low with 20 um cells and 0.5 ms sub-steps, 1.1 % with 10 um and 0.25 ms (measured while planning): the
    # first-order error of removing whole cells once they are liquid
    assert speed == pytest.approx(v, rel=0.03)
    n = sk.count[0]
    dz = sk.dz()[0, :n]
    z = np.cumsum(dz) - 0.5 * dz
    exact = T0 + (mat.T_feed - T0) * np.exp(-v * z / alpha)
    near = z < 3e-3
    assert np.abs(sk.T[0, :n][near] - exact[near]).max() < 0.02 * (mat.T_feed - T0)    # 9.2 K measured


def test_the_books_of_a_pass_close_exactly():
    mat = scheil()
    s = skin.SkinSettings(substep=0.01)
    sk = filled(mat, n=3, settings=s, T=880.0)
    m_min = skin.minimum_mass(mat, sk.area, s)
    E0, M0 = sk.energy(), sk.mass()
    q = np.array([5e6, 2e7, 5e7])           # 0.4 mm of skin passes 10 MW/m2 to the base with a 31 K rise: more melts
    book = skin.advance(sk, 0.2, q, sk.area * 1.02, 300.0, 0.3, np.full(3, 870.0), m_min, T_b1=np.full(3, 875.0))
    assert book.removed.sum() > 0.0 and book.drawn.sum() > 0.0
    assert sk.mass() == pytest.approx(M0 - book.removed.sum() + book.drawn.sum(), rel=1e-13)
    gained = book.E_top - book.E_rad - book.E_base + book.E_side - book.removed_H + book.drawn_H
    assert sk.energy() == pytest.approx(E0 + gained.sum(), rel=1e-12)
    assert np.allclose(book.E_top, q * sk.area * 1.02 * 0.2)


def test_the_trial_pass_leaves_the_skins_and_gives_the_real_pass_its_base_energy():
    mat = scheil()
    s = skin.SkinSettings(substep=0.01)
    sk = filled(mat, n=2, settings=s, T=850.0)
    m_min = skin.minimum_mass(mat, sk.area, s)
    before = sk.state()
    q = np.array([1e6, 4e6])
    a, b = skin.advance(sk, 0.1, q, sk.area, 300.0, 0.3, np.full(2, 840.0), m_min)
    assert all(np.array_equal(x, y) for x, y in zip(before, sk.state()))
    assert np.all(b <= 0.0)
    for dT, tol in ((0.0, 1e-9), (2.0, 2e-3)):
        sk.set_state(before)
        book = skin.advance(sk, 0.1, q, sk.area, 300.0, 0.3, np.full(2, 840.0), m_min, T_b1=np.full(2, 840.0 + dT))
        law = (a + b * (840.0 + dT)) * 0.1
        # exact at the linearisation temperature; 2 K away the law misses by a few parts in a thousand of the energy
        # through the top (the mismatch the body books)
        assert np.abs(book.E_base - law).max() <= tol * np.abs(book.E_top).max()


def test_a_draw_from_a_base_above_t_rigid_takes_rigid_metal():
    """The skin draws rigid metal (spec 6.6): from a 3D base above T_rigid (895 K) it takes metal at T_rigid's enthalpy,
    never the base's -- a base overshooting T_feed would otherwise feed liquid that the next sub-step removes, without
    limit. (At a 930 K base the whole 0.4 mm skin passes T_feed within one 10 ms sub-step, is removed and redrawn: the
    runaway the interface's repeat stops, Task 10.)"""
    mat = scheil()
    s = skin.SkinSettings(substep=0.01)
    sk = filled(mat, n=1, settings=s, T=880.0)
    m_min = skin.minimum_mass(mat, sk.area, s)
    sk.take_bottom(0.5 * m_min)
    book = skin.advance(sk, 0.01, np.zeros(1), sk.area, 0.0, 0.0, np.array([900.0]), m_min, T_b1=np.array([900.0]))
    assert book.drawn[0] == pytest.approx(0.5 * m_min[0], rel=1e-12)
    assert book.drawn_H[0] == pytest.approx(book.drawn[0] * float(mat.enthalpy(mat.T_rigid)), rel=1e-12)


```

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: ... no attribute 'advance'` (and `_conduct`).

- [ ] **Step 3: Append to `reentry_model/skin.py`**

```python
# -- heat through the depth (Task 4) -----------------------------------------------------------------------------------
@dataclass
class SkinBook:
    """What a pass did, per column: energies [J] and masses [kg]."""
    E_top: np.ndarray        # convective heat taken in at the top
    E_rad: np.ndarray        # radiated from the top
    E_base: np.ndarray       # conducted down into the 3D model through the base
    E_side: np.ndarray       # gained from the neighbouring columns by sideways conduction
    removed: np.ndarray      # liquid removed as it formed (stage 2)
    removed_H: np.ndarray    # the enthalpy it took
    drawn: np.ndarray        # rigid metal drawn from the 3D model in the sub-steps
    drawn_H: np.ndarray      # the enthalpy it brought


def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, Q_side=None):
    """One implicit (backward-Euler) sub-step for every column: the convective flux q_top [W/m2] over the top area
    A_top in at the top cell, radiation from it over the same area, conduction through the cells (harmonic
    conductances), the explicit sideways heat Q_side [W per cell], and the base held at T_b [K] through the bottom
    half-cell (T_b None: an adiabatic base, for the tests). Newton in the temperature with a backtracking line search
    per column (stage 1's solver, sub-plan 18). Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the
    skins track S, updates it to dT/dT_b and returns the base energy's derivative with respect to T_b as a fifth item.
    E_base is what the books leave, so they close exactly whatever the iteration error."""
    mat = sk.props
    n, K = sk.m.shape
    valid = sk.valid
    rows = np.arange(n)
    last = np.maximum(sk.count - 1, 0)
    A = sk.area
    dz = np.where(valid, np.maximum(sk.dz(), DZ_FLOOR), 1.0)
    m = sk.m
    link = valid[:, 1:]
    T_old = sk.T.copy()
    E_old = sk.H.copy()                              # the stored enthalpy: the books close whatever h(T) is tabulated
    eps_s = np.asarray(eps, dtype=float) * SIGMA_SB * np.asarray(A_top, dtype=float)
    q = np.asarray(q_top, dtype=float) * np.asarray(A_top, dtype=float)
    base = T_b is not None
    T_b = np.asarray(T_b if base else 0.0, dtype=float) * np.ones(n)
    Qs = np.zeros((n, K)) if Q_side is None else np.where(valid, Q_side, 0.0)

    def system(T):
        h, cp, k = mat.conduction(T)
        G = np.where(link, A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:]), 0.0)
        Gb = A * k[rows, last] / (0.5 * dz[rows, last]) if base else np.zeros(n)
        flux = G * (T[:, :-1] - T[:, 1:])
        net = Qs.copy()
        net[:, :-1] -= flux
        net[:, 1:] += flux
        net[:, 0] += q - eps_s * (T[:, 0] ** 4 - T_amb ** 4)
        net[rows, last] -= Gb * (T[rows, last] - T_b)
        diag = m * cp / dt
        diag[:, 0] += 4.0 * eps_s * T[:, 0] ** 3
        diag[:, :-1] += G
        diag[:, 1:] += G
        diag[rows, last] += Gb
        R = np.where(valid, (m * h - E_old) / dt - net, 0.0)
        return R, G, Gb, np.where(valid, diag, 1.0)

    # Newton in the temperature with a backtracking line search per column on the residual's norm: in 10-micrometre
    # cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step, so a cell must move
    # with its neighbours, and the 3D solver's per-node enthalpy-consistent update stalled here (sub-plan 18)
    T = T_old.copy()
    R, G, Gb, diag = system(T)
    moved = np.zeros(n)
    for _ in range(NEWTON_MAX):
        dT = thomas(-G, diag, -G, -R)
        r0 = np.sqrt((R * R).sum(axis=1))
        lam = np.ones(n)
        for _ in range(LINE_SEARCH_MAX):
            T_try = T + lam[:, None] * dT
            R_try, G_try, Gb_try, diag_try = system(T_try)
            moved = np.abs(T_try - T).max(axis=1)
            worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
            if not worse.any():
                break
            lam = np.where(worse, 0.5 * lam, lam)
        T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try
        if moved.max() <= NEWTON_TOL:
            break
    else:
        raise RuntimeError("a skin sub-step did not converge in {} iterations (last change {:.2e} K)".format(
            NEWTON_MAX, moved.max()))
    T = np.where(valid, T, T_old)
    E_new = m * mat.enthalpy(T)
    E_top = q * dt
    E_rad = eps_s * (T[:, 0] ** 4 - T_amb ** 4) * dt
    E_side = Qs.sum(axis=1) * dt
    E_base = E_top - E_rad + E_side - (E_new.sum(axis=1) - E_old.sum(axis=1))
    out = [E_top, E_rad, E_base, E_side]
    if sk.S is not None:
        rhs = m * mat.cp_eff(T_old) / dt * sk.S
        rhs[rows, last] += Gb
        S = thomas(-G, diag, -G, np.where(valid, rhs, 0.0))
        out.append(dt * Gb * (S[rows, last] - 1.0))
        sk.S = np.where(valid, S, 0.0)
    sk.T, sk.H = T, np.where(valid, E_new, 0.0)
    return tuple(out)


def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True):
    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the sideways heat `side(sk, h)` [W per
    cell] if given (Task 5), the implicit conduction, then -- with `remove` and a melting material -- the contiguous top
    cells at or above T_feed leave, a column below its minimum mass m_min [kg] draws back to it, and the cells are split
    and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but never above T_rigid -- a 3D base
    node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it reached 930 K in one 0.5 s step, found
    while planning), and drawing at its temperature would draw liquid that the next sub-step removes, without limit;
    what the 3D model actually loses is reconciled by the body after the base moves.

    T_b1 None: the trial pass, on a copy (the skins are left as they were), base held at T_b0; returns the base law
    (a, b) per column, E_base / dt = a + b T_b [W into the 3D model], from the pass's base energy and its derivative with
    respect to the base temperature, the slope clipped at zero (stage 1's choice: a law rising with T_b would cost the 3D
    matrix its diagonal dominance; what the clip misses is the interface mismatch's to book). Otherwise the real pass,
    base held at T_b1, the 3D model's end-of-step base temperature; returns a SkinBook."""
    mat = sk.material
    N = max(1, int(math.ceil(dt / sk.settings.substep - 1e-9)))
    h = dt / N
    trial = T_b1 is None
    saved = sk.state() if trial else None
    T_b = np.asarray(T_b0 if trial else T_b1, dtype=float) * np.ones(sk.n)
    book = SkinBook(*(np.zeros(sk.n) for _ in range(8)))
    dEb = np.zeros(sk.n)
    if trial:
        sk.S = np.zeros_like(sk.m)
    h_b = mat.enthalpy(np.minimum(T_b, mat.T_rigid))
    melting = remove and mat.melts
    for _ in range(N):
        Q_side = side(sk, h) if side is not None else None
        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, Q_side)
        book.E_top += res[0]
        book.E_rad += res[1]
        book.E_base += res[2]
        book.E_side += res[3]
        if trial:
            dEb += res[4]
        if melting:
            k = np.cumprod(sk.valid & (sk.T >= mat.T_feed), axis=1).sum(axis=1)
            gone, gone_H = sk.take_top(k)
            book.removed += gone
            book.removed_H += gone_H
        short = np.maximum(m_min - sk.column_mass(), 0.0)
        book.drawn += short
        book.drawn_H += sk.add_bottom(short, h_b)
        sk.regrid()
    if trial:
        sk.set_state(saved)
        b = np.minimum(dEb / dt, 0.0)
        return book.E_base / dt - b * T_b, b
    sk.S = None
    return book


```

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

```python
# -- heat sideways (Task 5) --------------------------------------------------------------------------------------------
def grid_skins(mat, n=7, h=1e-3, T=700.0, settings=None):
    """Columns on a flat n x n grid of h squares split in two triangles (open edges), each filled to its minimum."""
    from reentry_model import columns
    x, y = np.meshgrid(np.arange(n) * h, np.arange(n) * h, indexing="ij")
    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
    idx = np.arange(n * n).reshape(n, n)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    col = columns.Columns(pts, np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])]))
    s = settings or skin.SkinSettings(cell=20e-6)
    sk = skin.Skins(mat, col.area, s)
    sk.fill(skin.minimum_mass(mat, sk.area, s), T)
    edge = (x.ravel() == 0) | (y.ravel() == 0) | (x.ravel() == x.max()) | (y.ravel() == y.max())
    return sk, col, pts, edge


def test_sideways_heat_sums_to_zero_between_columns_of_any_depth(coarse_sphere_mesh):
    from test_reentry_model_columns import sphere_columns
    mat = scheil()
    col, _ = sphere_columns(coarse_sphere_mesh)
    sk = skin.Skins(mat, col.area, skin.SkinSettings())
    rng = np.random.default_rng(4)
    sk.fill(skin.minimum_mass(mat, sk.area, sk.settings) * rng.uniform(1.0, 3.0, col.n), rng.uniform(500.0, 900.0, col.n))
    sk.regrid()
    sk.T = np.where(sk.valid, sk.T + rng.normal(0.0, 20.0, sk.T.shape), sk.T)          # every cell different
    sk.H = sk.m * mat.enthalpy(sk.T)
    Q = skin.sideways_heat(sk, col, 0.01)
    assert abs(Q.sum()) < 1e-12 * np.abs(Q).sum()
    assert np.all(Q[~sk.valid] == 0.0)


def test_sideways_heat_is_exact_for_a_field_linear_along_the_surface():
    mat = sharp()                                                    # constant k: a linear field carries no net flux
    sk, col, pts, edge = grid_skins(mat)
    T = 600.0 + 2e4 * pts[:, 0] - 1e4 * pts[:, 1]                       # 20 K/mm and -10 K/mm, uniform in depth
    sk.T = np.where(sk.valid, T[:, None], sk.T)
    sk.H = sk.m * mat.enthalpy(sk.T)
    Q = skin.sideways_heat(sk, col, 1e-4)                            # one cycle: the edge columns do not move first
    assert np.abs(Q[~edge]).max() < 1e-9 * np.abs(Q[edge]).max()


def test_conduction_between_columns_reaches_the_known_steady_temperature():
    """Columns on a flat plate, insulated top and bottom, the boundary columns held at a linear field: the steady state
    is that linear field everywhere (Laplace's equation, with a constant conductivity)."""
    mat = sharp()
    sk, col, pts, edge = grid_skins(mat, n=6)
    T_edge = 600.0 + 3e4 * pts[:, 0] + 1e4 * pts[:, 1]
    # 4 ms keeps one sideways cycle per step even at the grid's corner columns (a sixth of the interior's area): with
    # two, the held columns would drift between cycles and the interior settle 1.8 K off (found while planning)
    for _ in range(200):                                             # 0.8 s: the slowest mode decays e^-22 over it
        sk.T[edge] = np.where(sk.valid[edge], T_edge[edge][:, None], 300.0)
        sk.H = sk.m * mat.enthalpy(sk.T)
        skin._conduct(sk, 0.004, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.sideways_heat(sk, col, 0.004))
    err = np.abs(sk.T[~edge][sk.valid[~edge]] - np.repeat(T_edge[~edge], sk.count[~edge]))
    assert err.max() < 1e-9 * (T_edge.max() - T_edge.min())           # 3e-13 K measured: exact


def test_a_sub_step_beyond_the_explicit_limit_is_sub_cycled_and_stays_bounded():
    mat = scheil()
    sk, col, pts, _ = grid_skins(mat, n=6)
    hot = (np.rint(pts[:, 0] / 1e-3) + np.rint(pts[:, 1] / 1e-3)) % 2 == 0          # a checkerboard
    sk.T = np.where(sk.valid, np.where(hot, 800.0, 400.0)[:, None], sk.T)
    sk.H = sk.m * mat.enthalpy(sk.T)
    E0 = sk.energy()
    Q = skin.sideways_heat(sk, col, 0.2)                             # 0.2 s: about 30 times the limit on 1 mm columns
    sk.H = sk.H + Q * 0.2
    T = mat.temperature_from_enthalpy(sk.H[sk.valid] / sk.m[sk.valid])
    assert sk.energy() == pytest.approx(E0, rel=1e-13)
    assert T.min() >= 400.0 - 1e-6 and T.max() <= 800.0 + 1e-6
```

- [ ] **Step 2: Run them to see them fail.** Expected: `AttributeError: ... no attribute 'sideways_heat'`.

- [ ] **Step 3: Append to `reentry_model/skin.py`**

```python
# -- heat sideways (Task 5) --------------------------------------------------------------------------------------------
_SPAN = 1.0                       # m, larger than any column: separates the rows of a flattened row-wise search


def _interp_rows(x, xp, fp):
    """Row-wise linear interpolation of fp (n, q) given at non-decreasing xp (n, q), at x (n, p), held at the ends."""
    n, q = xp.shape
    off = (np.arange(n) * _SPAN)[:, None]
    i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
    i = np.clip(i - (np.arange(n) * q)[:, None], 0, q - 2)
    r = np.arange(n)[:, None]
    x0, x1, f0, f1 = xp[r, i], xp[r, i + 1], fp[r, i], fp[r, i + 1]
    t = np.clip(np.divide(x - x0, x1 - x0, out=np.zeros_like(x, dtype=float), where=x1 > x0), 0.0, 1.0)
    return f0 + t * (f1 - f0)


def depth_grid(settings, depth):
    """The common depth grid of the sideways exchange: breakpoints from the top down to `depth` [m], spaced by the cell
    target at each depth (the fine cell, growing below the fine zone)."""
    z = [0.0]
    while z[-1] < depth:
        z.append(z[-1] + settings.cell + (settings.growth - 1.0) * max(z[-1] - settings.min_thickness, 0.0))
    return np.array(z)


def _sideways_once(sk, columns, T, dz):
    """Sideways heat [W per cell] at temperatures T: each column's cells averaged onto the depth grid, every surface
    edge exchanging k w l (T_j - T_i) in every bin both columns reach (k the mean of the two bins' conductivities, w the
    edge's cotangent weight, l the depth the two columns share in the bin), and each column's bin heat given back to its
    cells in proportion to their overlap with the bin. Antisymmetric edge by edge, so the total is zero to round-off;
    exact for a field linear along the surface where neighbouring columns have the same cells."""
    n, K = T.shape
    zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)          # cell boundaries, from the top
    H = zb[:, -1]
    zg = depth_grid(sk.settings, float(H.max()))
    zc = np.minimum(zg[None, :], H[:, None])                                        # the grid clipped to each column
    C = np.concatenate([np.zeros((n, 1)), np.cumsum(T * dz, axis=1)], axis=1)      # integral of T from the top
    Cg = _interp_rows(zc, zb, C)
    length = np.diff(zc, axis=1)                                                    # each column's share of each bin
    Tbar = np.divide(np.diff(Cg, axis=1), length, out=np.zeros_like(length), where=length > 0.0)
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    k = sk.props.k(Tbar)
    flux = (0.5 * (k[i] + k[j]) * columns.conduct_weight[:, None] * np.minimum(length[i], length[j])
            * (Tbar[j] - Tbar[i]))
    Qg = np.zeros_like(length)
    np.add.at(Qg, i, flux)
    np.add.at(Qg, j, -flux)
    Qcum = np.concatenate([np.zeros((n, 1)), np.cumsum(Qg, axis=1)], axis=1)
    return np.diff(_interp_rows(zb, zc, Qcum), axis=1)


def sideways_heat(sk, columns, dt):
    """The sideways conduction of one sub-step dt [s] as an explicit heat source [W per cell], the time average over
    as many sub-cycles as keep dt k_max sum(w) / (rho c_p,min A) below SIDEWAYS_CFL (about 25 ms per cycle at 2 mm
    spacing, spec 6.3); the cells' temperatures move between sub-cycles, the source being what `_conduct` then adds."""
    mat = sk.material
    valid = sk.valid
    dz = sk.dz()
    wsum = np.bincount(columns.edges.ravel(), weights=np.repeat(columns.conduct_weight, 2), minlength=sk.n)
    rate = float((mat.k_table.max() * wsum / (mat.rho * mat.cp_table.min() * sk.area)).max())
    cycles = max(1, int(math.ceil(dt * rate / SIDEWAYS_CFL)))
    T, H = sk.T.copy(), sk.H.copy()
    total = np.zeros_like(T)
    for _ in range(cycles):
        Q = np.where(valid, _sideways_once(sk, columns, T, dz), 0.0)
        total += Q / cycles
        if cycles > 1:
            H = H + Q * dt / cycles
            T = np.where(valid, mat.temperature_from_enthalpy(np.divide(H, sk.m, out=np.zeros_like(H), where=sk.m > 0.0)), T)
    return total


```

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

```diff
--- a/tests/test_reentry_model_thermal.py
+++ b/tests/test_reentry_model_thermal.py
@@ -355,3 +355,34 @@
     s.set_state(snap)
     r2 = s.step(0.5, q, 0.0)
     assert np.allclose(r1.T, r2.T, rtol=0.0, atol=1e-6)
+
+
+def test_a_nodal_law_enters_implicitly_and_keeps_the_balance(coarse_sphere_mesh):
+    """The corner columns' base law (sub-plan 20 Task 6): heat a + b T_i [W] into listed nodes, implicit; with b = 0 it
+    is a nodal load exactly; with b < 0 the step's energy change is what the law delivered at the end-of-step
+    temperatures. No radiation (emissivity 0: the boundary lies under the skins)."""
+    mat = material.Material.from_drama_json("AA7075_scheil")
+    corners = np.unique(coarse_sphere_mesh.surface().faces)
+    a = np.full(len(corners), 0.5)
+    out = []
+    for kind in ("load", "law0"):
+        s = thermal.thermal_solver("skfem", linear_solver="direct")
+        s.setup(coarse_sphere_mesh, mat, 0.0)
+        s.set_temperature(700.0)
+        nf = len(s.areas)
+        if kind == "load":
+            load = np.zeros(coarse_sphere_mesh.n_nodes)
+            load[corners] = a
+            out.append(s.step(0.5, np.zeros(nf), 0.0, nodal_load=load))
+        else:
+            out.append(s.step(0.5, np.zeros(nf), 0.0, nodal_law=(corners, a, np.zeros(len(corners)))))
+    assert np.allclose(out[0].T, out[1].T, rtol=0.0, atol=1e-9) and out[1].Q_interface == pytest.approx(a.sum(), rel=1e-12)
+    s = thermal.thermal_solver("skfem", linear_solver="direct")
+    s.setup(coarse_sphere_mesh, mat, 0.0)
+    s.set_temperature(700.0)
+    b = np.full(len(corners), -2e-3)                                 # 2 mW/K per node: the base cools toward 950 K
+    E0 = s.energy()
+    r = s.step(0.5, np.zeros(len(s.areas)), 0.0, nodal_law=(corners, -950.0 * b, b))
+    T = s.temperature()
+    assert r.Q_interface == pytest.approx(float((-950.0 * b + b * T[corners]).sum()), rel=1e-12)
+    assert s.energy() - E0 == pytest.approx(r.Q_interface * 0.5, rel=1e-6) and r.Q_rad == 0.0
```

  and to `tests/test_reentry_model_fenicsx.py`:

```python
def test_backends_agree_with_a_nodal_law(coarse_sphere_mesh):
    """The corner columns' base law (sub-plan 20 Task 6) at the surface nodes, no radiation, in both backends."""
    from reentry_model import material, thermal
    corners = np.unique(coarse_sphere_mesh.surface().faces)
    b = np.full(len(corners), -2e-3)
    out = []
    for name in ("skfem", "fenicsx"):
        s = thermal.thermal_solver(name)
        s.setup(coarse_sphere_mesh, material.Material.from_drama_json("AA7075_scheil"), 0.0)
        s.set_temperature(700.0)
        r = s.step(0.5, np.zeros(len(s.areas)), 0.0, nodal_law=(corners, -950.0 * b, b))
        out.append((r, s.temperature()))
    (ra, Ta), (rb, Tb) = out
    assert np.allclose(Ta, Tb, rtol=0.0, atol=1e-6) and ra.Q_interface == pytest.approx(rb.Q_interface, rel=1e-9)
```

- [ ] **Step 2: Run them to see them fail.** Expected: `TypeError: ... got an unexpected keyword argument 'nodal_law'`.

- [ ] **Step 3: Apply the backends' changes**

```diff
--- a/reentry_model/thermal/skfem_backend.py
+++ b/reentry_model/thermal/skfem_backend.py
@@ -253,11 +253,13 @@
         h = self.material.enthalpy(self.T if T is None else T)
         return self.phi * self.material.rho * self.vol * h[self.tets].mean(axis=1)
 
-    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None):
+    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None, nodal_law=None):
         """One backward-Euler step with the per-facet convective flux q_conv [W/m2] and an optional nodal load vector
         [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped).
         `interface = (mask, a, b)` per facet: the skins' base law, heat a + b T_f [W/m2] into the body through the
-        masked facets, which neither take q_conv's place nor radiate (Step 3 sub-plan 18)."""
+        masked facets, which neither take q_conv's place nor radiate (Step 3 sub-plan 18). `nodal_law = (nodes, a, b)`:
+        the corner columns' base law (layered-skin design section 7.2, sub-plan 20), heat a + b T_i [W] into each listed
+        node, implicit like the radiation."""
         T_old = self.T
         # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
         T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
@@ -278,6 +280,10 @@
             b_int = np.where(mask, np.asarray(interface[2], dtype=float), 0.0)
         eps_f = self.emissivity if mask is None else np.where(mask, 0.0, self.emissivity)
         self._radiating = None if mask is None else ~mask
+        if nodal_law is not None:
+            law_nodes = np.asarray(nodal_law[0], dtype=np.int64)
+            a_law, b_law = np.asarray(nodal_law[1], dtype=float), np.asarray(nodal_law[2], dtype=float)
+            B_law = sp.csr_matrix((-b_law, (law_nodes, law_nodes)), shape=(len(self.points), len(self.points)))
         for iteration in range(1, self.max_iterations + 1):
             K, M, E = self.operators(T_k, T_old)
             Tf = self.facet_temperature(T_k)
@@ -287,6 +293,9 @@
             if mask is not None:
                 B = B + self.facet_pattern.assemble(self.Bf * (-b_int)[:, None, None])
                 F_int = self.facet_load(a_int + b_int * Tf)
+            if nodal_law is not None:
+                B = B + B_law
+                F_int = F_int + np.bincount(law_nodes, weights=a_law + b_law * T_k[law_nodes], minlength=len(self.points))
             # Newton on R(T) = E(T)/dt + K T - F_conv + F_rad(T) - F_extra with the tangent M/dt + K + B:
             # A T_new = A T_k - R(T_k); with a smooth c_p (c_tan = c_sec) this is the Step 2 iteration unchanged
             A = (M / dt + K + B).tocsr()
@@ -323,6 +332,8 @@
         self._T_prev, self.T = T_old, T_new
         self.last_damping = damping
         Q_int = float(((a_int + b_int * self.facet_temperature(T_new)) * self.areas).sum()) if mask is not None else 0.0
+        if nodal_law is not None:
+            Q_int += float((a_law + b_law * T_new[law_nodes]).sum())
         return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()),
                           Q_dropped, Q_int)
 
```

```diff
--- a/reentry_model/thermal/fenicsx_backend.py
+++ b/reentry_model/thermal/fenicsx_backend.py
@@ -205,7 +205,7 @@
         m.sum_duplicates()
         return PETSc.Mat().createAIJ(size=m.shape, csr=(m.indptr.astype(np.int32), m.indices.astype(np.int32), m.data), comm=self.msh.comm)
 
-    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None):
+    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None, nodal_law=None):
         from petsc4py import PETSc
         mat, n = self.material, self.mesh.n_nodes
         T_old = self.T
@@ -229,6 +229,11 @@
             b_int = np.where(mask, np.asarray(interface[2], dtype=float), 0.0)
         eps_f = self.emissivity if mask is None else np.where(mask, 0.0, self.emissivity)
         self._radiating = None if mask is None else ~mask
+        # the corner columns' base law (sub-plan 20): heat a + b T_i [W] into each listed node, implicit
+        if nodal_law is not None:
+            law_nodes = np.asarray(nodal_law[0], dtype=np.int64)
+            a_law, b_law = np.asarray(nodal_law[1], dtype=float), np.asarray(nodal_law[2], dtype=float)
+            B_law = sp.csr_matrix((-b_law, (law_nodes, law_nodes)), shape=(n, n))
         converged, last_relative_change, previous_change, damping, decreases = False, None, None, 1.0, 0
         T_new, iteration = T_k, 0
         for iteration in range(1, self.max_iterations + 1):
@@ -248,6 +253,9 @@
             if mask is not None:
                 B = B + sp.csr_matrix((np.repeat(-b_int * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
                 F_int = self.facet_load(a_int + b_int * Tf)
+            if nodal_law is not None:
+                B = B + B_law
+                F_int = F_int + np.bincount(law_nodes, weights=a_law + b_law * T_k[law_nodes], minlength=n)
             b = M_tan * T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra + F_int
             # operator in dof numbering: K (dolfinx) + diag(M_tan/dt) + B
             A = self.K.copy()
@@ -292,5 +300,7 @@
         self._T_prev, self.T = T_old, T_new
         self.last_damping = damping
         Q_int = float(((a_int + b_int * self.facet_temperature(T_new)) * self.areas).sum()) if mask is not None else 0.0
+        if nodal_law is not None:
+            Q_int += float((a_law + b_law * T_new[law_nodes]).sum())
         return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()),
                           Q_dropped, Q_int)
```

```diff
--- a/reentry_model/thermal/__init__.py
+++ b/reentry_model/thermal/__init__.py
@@ -24,13 +24,14 @@
     iterations: int               # Newton/Picard iterations taken
     Q_extra: float = 0.0          # W, nodal load applied over the step (Step 3 deferred melt energy)
     Q_dropped: float = 0.0        # W, load that fell on pinned (material-free) nodes and was not applied
-    Q_interface: float = 0.0      # W, heat the skins passed in through their facets (Step 3 sub-plan 18)
+    Q_interface: float = 0.0      # W, heat the skins passed in through their facets (Step 3 sub-plan 18) or, with the
+                                  # corner columns, through the base nodes (sub-plan 20)
 
 
 class ThermalSolver(Protocol):
     def setup(self, mesh, material, emissivity) -> None: ...
     def set_temperature(self, T) -> None: ...                       # uniform float or nodal array
-    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None) -> StepResult: ...
+    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None, nodal_law=None) -> StepResult: ...
     def temperature(self) -> np.ndarray: ...
     def energy(self) -> float: ...                                  # stored enthalpy above material.T_REF [J]
     def element_energies(self, T=None) -> np.ndarray: ...           # phi_e rho V_e h(T_e) per element [J]
```

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

```python
"""The skin body (layered-skin design sections 6.6 and 7, sub-plan 20 Task 7): skins on every corner from time zero, the
trial and real passes around the 3D step, the base moving the 3D mesh, and the books."""
import numpy as np
import pytest

from reentry_model import heating, material, mesh, skin, skinbody
from reentry_model.thermal import thermal_solver

DEV = skin.SkinSettings(cell=20e-6)


def make(m, T0=880.0, name="skfem", settings=DEV, **kw):
    copy = mesh.VolumeMesh(m.points.copy(), m.tets.copy(), dict(m.params))
    mat = material.Material.from_drama_json("AA7075_scheil")
    return skinbody.SkinBody(copy, mat, thermal_solver(name), 1.47, settings, T0=T0, v_hat=(-1.0, 0.0, 0.0), **kw)


def front_loads(b, q0):
    return heating.HeatingResult(q0 * np.clip(np.cos(b.theta), 0.0, None), q0, 0.0, 0.0, 0.0)


def test_every_corner_starts_with_its_minimum_taken_from_the_3d_model(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh)
    m_min = skin.minimum_mass(b.material, b.columns.area, DEV)
    V0 = coarse_sphere_mesh.volume()
    assert b.skins.n == b.columns.n == len(np.unique(coarse_sphere_mesh.surface().faces))
    assert b.mass() == pytest.approx(b.material.rho * V0, rel=1e-14)                  # nothing made or lost
    assert b.energy() == pytest.approx(b.material.rho * V0 * float(b.material.enthalpy(880.0)), rel=1e-12)
    # each column holds what its corner swept moving in by its minimum: within 2 % of it (prisms against dual areas
    # on 4 mm triangles), so the first sub-step draws any shortfall
    assert np.abs(b.skins.column_mass() / m_min - 1.0).max() < 0.02
    assert np.allclose(b.skins.T[b.skins.valid], 880.0)
    assert b.surface.n_patches == b.columns.n and b.theta[b.i_stag] < np.radians(5.0)


def test_a_body_with_no_heat_stays_uniform_and_its_mesh_and_skins_keep_their_books(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh, T0=600.0)
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.emissivity = 0.0
    M0, E0 = b.mass(), b.energy()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, zero)
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 600.0).max() < 1e-6
    assert np.abs(b.skins.T[b.skins.valid] - 600.0).max() < 1e-6
    assert b.removed_mass == 0.0 and b.jump_faults == 0


def test_heating_with_removal_keeps_mass_and_energy_exact(coarse_sphere_mesh):
    """6 MW/m2 on the front for 1.5 s from 880 K: liquid forms and leaves, columns draw, the base recedes; the mass
    removed plus the mass left is the mass at start, and the energy books close to 1e-8 of the heat absorbed."""
    b = make(coarse_sphere_mesh)
    M0 = b.mass()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, front_loads(b, 6e6))
    assert b.removed_mass > 0.0 and b.drawn_mass > 0.0
    assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8                  # 1.6e-9 measured while planning
    stats = b.skin_stats()
    assert stats["base_recession_max_mm"] > 0.5                       # 0.88 mm at the nose (measured)
    assert b.surface_temperature().max() <= b.material.T_feed + 1.0   # liquid left as it formed
    assert stats["skin_cells"] == b.skins.count.sum() and b.last.Q_conv > 0.0


def test_hand_down_passes_rigid_bottom_metal_beyond_the_minimum():
    mat = material.Material.from_drama_json("AA7075_scheil")
    sk = skin.Skins(mat, np.full(3, 4e-6), DEV)
    m_min = skin.minimum_mass(mat, sk.area, DEV)
    sk.fill(1.5 * m_min, np.array([700.0, 700.0, 905.0]))            # two cold, one slurry (above T_rigid)
    sk.take_bottom(np.array([0.0, 0.6, 0.0]) * m_min)                 # the second below its minimum
    M0, E0 = sk.column_mass(), sk.H.sum(axis=1)
    handed, H = skin.hand_down(sk, m_min, mat.T_rigid)
    assert handed == pytest.approx(np.array([0.5 * m_min[0], 0.0, 0.0]), rel=1e-12)
    assert sk.column_mass() == pytest.approx(M0 - handed, rel=1e-13)
    assert sk.H.sum(axis=1) == pytest.approx(E0 - H, rel=1e-13)
    assert H[0] == pytest.approx(handed[0] * float(mat.enthalpy(700.0)), rel=1e-9)


def test_handing_down_moves_the_base_out_and_the_3d_model_takes_the_metal(coarse_sphere_mesh):
    """A cold skin above its minimum hands its bottom metal to the 3D model: the corners move out by its volume, the 3D
    model's volume grows by it, and mass and energy are exact."""
    b = make(coarse_sphere_mesh, T0=600.0)
    b.emissivity = 0.0
    extra = 0.5 * skin.minimum_mass(b.material, b.columns.area, DEV)
    b.skins.add_bottom(extra, float(b.material.enthalpy(600.0)))
    b.skins.regrid()
    M0, E0, V0 = b.mass(), b.energy(), b.mesh.volume()
    # every column's excess over its minimum, the added metal and what it held above its minimum from the start
    excess = np.maximum(b.skins.column_mass() - skin.minimum_mass(b.material, b.columns.area, DEV), 0.0).sum()
    b.advance(0.5, 0.5, heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0))
    assert b.handed_mass == pytest.approx(excess, rel=1e-9) and excess > extra.sum()
    assert b.mesh.volume() - V0 == pytest.approx(excess / b.material.rho, rel=0.02)
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 600.0).max() < 1e-6


# -- what triggers the interface repeat (Task 7, measured in Task 10) --------------------------------------------------
def test_every_repeat_rule_keeps_the_books_and_records_what_triggered_it(coarse_sphere_mesh):
    """Each repeat rule (measured in Task 10; the spec's per-column rule is the default) leaves the books exact, and the
    first pass's mismatch is recorded whatever the rule. 2 MW/m2 from 880 K on 4 mm elements: without a repeat, 3 MW/m2
    already runs away -- the 3D base overshoots T_feed and every sub-step removes and redraws a whole skin (measured
    while planning)."""
    results = {}
    for rule in ("column", "total", "never"):
        b = make(coarse_sphere_mesh, repeat_rule=rule)
        M0 = b.mass()
        for k in range(3):
            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 2e6))
        assert b.removed_mass > 0.0
        assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
        assert abs(b.energy_balance_residual()) < 1e-8
        stats = b.skin_stats()
        assert stats["mismatch_total_ratio"] > 0.0 and stats["mismatch_worst_ratio"] >= stats["mismatch_total_ratio"]
        results[rule] = b.interface_repeats
    assert results["never"] == 0 and results["column"] >= results["total"] > 0
    with pytest.raises(ValueError):
        make(coarse_sphere_mesh, repeat_rule="sometimes")
```

  and append to `tests/test_reentry_model_fenicsx.py`:

```python
def test_the_skin_body_matches_the_skfem_backend(coarse_sphere_mesh):
    """The corner columns (sub-plan 20 Task 7) in both backends: 6 MW/m2 on the front for two 0.5 s steps from 880 K,
    with liquid removed as it forms. Masses, removal and energy agree, both books are exact, and the fields agree except
    where the solvers' round-off (about 1e-7 K) tipped a cell across a split, merge or removal threshold: measured while
    planning, 2 of 2254 columns ended with a different cell count, and conduction spread that to 145 of 2934 3D nodes,
    the largest difference 0.65 K."""
    from reentry_model import heating, mesh, skin, skinbody
    out = {}
    for name in ("skfem", "fenicsx"):
        np.random.seed(12345)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
        b = skinbody.SkinBody(m, material.Material.from_drama_json("AA7075_scheil"), thermal.thermal_solver(name), 1.47,
                              skin.SkinSettings(cell=20e-6), T0=880.0, v_hat=(-1.0, 0.0, 0.0))
        for k in range(2):
            q = 6e6 * np.clip(np.cos(b.theta), 0.0, None)
            b.advance(0.5 * (k + 1), 0.5, heating.HeatingResult(q, 6e6, 0.0, 0.0, 0.0))
        out[name] = b
    s, f = out["skfem"], out["fenicsx"]
    assert f.mass() == pytest.approx(s.mass(), rel=1e-12) and f.removed_mass == pytest.approx(s.removed_mass, rel=1e-9)
    assert f.energy() == pytest.approx(s.energy(), rel=1e-8)
    d = np.abs(f.solver.temperature() - s.solver.temperature())
    assert np.mean(d < 1e-3) > 0.9 and d.max() < 2.0
    assert np.mean(f.skins.count == s.skins.count) > 0.99
    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
```

- [ ] **Step 2: Run them to see them fail.** Expected: `ImportError: cannot import name 'skinbody'`.

- [ ] **Step 3: Append `hand_down` to `reentry_model/skin.py`**

```python
# -- the base (Task 7) -------------------------------------------------------------------------------------------------
def hand_down(sk, m_min, T_rigid):
    """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
    oldest first, as far as the column stays at or above its minimum mass m_min (spec 6.6); a cell may go in part.
    Returns (mass, enthalpy) handed down per column."""
    n, K = sk.m.shape
    j = np.arange(K)[None, :]
    idx = np.clip(sk.count[:, None] - 1 - j, 0, K - 1)                 # the cells from the bottom up
    rows = np.arange(n)[:, None]
    rigid = (j < sk.count[:, None]) & (sk.T[rows, idx] <= T_rigid)
    run = np.cumprod(rigid, axis=1).astype(bool)
    available = np.where(run, sk.m[rows, idx], 0.0).sum(axis=1)
    take = np.clip(np.minimum(available, sk.column_mass() - m_min), 0.0, None)
    return take, sk.take_bottom(take)
```

- [ ] **Step 4: Create `reentry_model/skinbody.py`**

```python
"""The body of the layered-skin design (docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md sections 4,
6.6 and 7; stage 2, sub-plan 20: heat only, liquid removed as it forms).

A moving 3D mesh holds the rigid metal; every corner of its surface carries a column with a layered skin holding the
metal outward of the rigid front (`skin.Skins`). At time zero each column takes its minimum mass from the 3D model at
the initial temperature: the corners move inward by it and the skins receive exactly what they swept. Per macro step,
after the trajectory and the heating:

1. the trial pass gives every column its base law, E_base / dt = a + b T_b, with the base held at the 3D model's
   start-of-step temperature (`skin.advance`);
2. the 3D step takes that law at the corner nodes, implicitly, with no convective load and no radiation on its
   boundary (the skins' tops radiate), and last step's deferred mismatch as a nodal load;
3. the real pass, base held at the new 3D temperature, books the heat, removes the liquid and draws to the minimum;
   the difference between what it sent through the base and what the 3D step received is deferred to the node, and a
   step whose mismatch exceeds INTERFACE_TOLERANCE of the column's throughput repeats once, re-linearised;
4. the base moves by the metal crossing it (spec 6.6): drawn (in the sub-steps) and absorbed metal move the corner in,
   handed-down metal moves it out; with the volume-neutral smoothing this is one displacement per corner, which
   `MovingMesh.move_with_policy` applies (rebuilding when needed), and each column is then reconciled with what its
   corner actually swept: its bottom gains the swept mass less what it already exchanged, and its enthalpy changes by
   the swept enthalpy less what it already exchanged, so mass and energy are exact.

The body's surface for the coupled run is the top surface, one 'patch' per corner (`columns.CornerSurface`), so the
heating, the history and the trajectory's loop are unchanged."""
import dataclasses
import time

import numpy as np

from . import columns as columns_mod
from . import movemesh
from . import skin as skin_mod
from .body import ThermalBody
from .thermal import SIGMA_SB

ABSORB_MAX = 0.25                 # of the mean surface edge: the most a base absorbs in one macro step (the rest waits)
GRADIENT_FLOOR = 1.0              # K/m: below this the 3D field gives no direction to absorb along
REPEAT_RULES = ("column", "total", "never")   # when a step's interface mismatch repeats it (spec 7.2: "column"; the others
                                              # are measured in sub-plan 20 Task 10 before any change of default)


class SkinBody(ThermalBody):
    """A ThermalBody whose surface is resolved by corner columns with layered skins on a moving mesh (module docstring).
    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass."""

    def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                 v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                 repeat_rule="column"):
        if repeat_rule not in REPEAT_RULES:
            raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
        self.repeat_rule = repeat_rule
        self.material, self.solver, self.mass_kg = material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient, self.remove, self.demise_fraction = float(T_ambient), bool(remove), float(demise_fraction)
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.settings = settings or skin_mod.SkinSettings()
        solver.setup(mesh, material, 0.0)                    # the 3D boundary lies under the skins: it does not radiate
        solver.set_temperature(T0)
        self.mass0 = material.rho * mesh.volume()
        self.moving = movemesh.MovingMesh(mesh, solver, material, depth=deform_depth)
        self.columns = columns_mod.Columns(mesh.points[self.moving.corners], self.moving.faces_corner)
        self.skins = skin_mod.Skins(material, self.columns.area, self.settings)
        self.pending = np.zeros(self.columns.n)              # J per corner: the deferred interface mismatch
        self.column_removed = np.zeros(self.columns.n)       # kg per column removed as liquid since time zero
        self.integrated_heat = self.absorbed_heat = self.radiated_heat = 0.0
        self.removed_mass = self.removed_enthalpy = 0.0
        self.drawn_mass = self.absorbed_mass = self.handed_mass = 0.0
        self.jump_faults = self.interface_repeats = self.absorb_capped = 0
        self.last_mismatch = self.last_shortfall = 0.0
        self.last_trigger = {"columns_over": 0, "worst_ratio": 0.0, "total_ratio": 0.0}    # the first pass's mismatch
        self.wall = {"skins": 0.0, "3d": 0.0, "move": 0.0, "rebuild": 0.0}
        self.iterations, self.last, self.last_move = [], None, None
        self._fill(T0)
        self.energy0 = self.energy()
        self._top0, self._normal0 = self.surface.points.copy(), self.surface.normals.copy()

    # -- set-up -----------------------------------------------------------------------------------------------------
    def _fill(self, T0):
        """Each column takes its minimum mass from the 3D model (spec 6.1): the corners move in by it, and the skins hold
        exactly what they swept, at the initial temperature."""
        m_min = skin_mod.minimum_mass(self.material, self.columns.area, self.settings)
        u = -self.columns.normal * (m_min / (self.material.rho * self.columns.area))[:, None]
        r = self.moving.move_with_policy(u)
        self.columns.update(self.moving.mesh.points[self.moving.corners])
        self.skins.area = self.columns.area.copy()
        self.skins.fill(self.material.rho * r.swept_volume, T0)
        self.skins.exchange_bottom(np.zeros(self.columns.n), r.swept_enthalpy - self.skins.H.sum(axis=1))
        self._refresh_geometry()                       # the top surface: the bases moved out by the filled skins

    def _refresh_geometry(self):
        """The columns, their areas and the top surface after the corners moved."""
        mm = self.moving
        self.columns.update(mm.mesh.points[mm.corners])
        self.skins.area = self.columns.area.copy()
        self.surface = columns_mod.CornerSurface.from_points(self.columns.top_points(self.skins.thickness()),
                                                             self.columns.faces)
        self.theta = self.surface.angles_to(self.v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)

    @property
    def mesh(self):
        return self.moving.mesh

    # -- the macro step ---------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
        mat, sk, mm = self.material, self.skins, self.moving
        q = np.asarray(loads.q_conv, dtype=float)
        A_top = self.surface.areas
        eps = np.full(sk.n, self.emissivity)
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        side = lambda s, h: skin_mod.sideways_heat(s, self.columns, h)
        nodes = mm.corners
        T_b0 = self.solver.temperature()[nodes]
        load = np.zeros(mm.mesh.n_nodes)
        load[nodes] = self.pending / dt
        solver_state, skin_state = self.solver.state(), sk.state()
        T_lin = T_b0
        for attempt in (0, 1):
            started = time.perf_counter()
            a, b = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_lin, m_min, side=side, remove=self.remove)
            self.wall["skins"] += time.perf_counter() - started
            started = time.perf_counter()
            res = self.solver.step(dt, np.zeros(len(self.solver.areas)), self.T_ambient, nodal_load=load,
                                   nodal_law=(nodes, a, b))
            self.wall["3d"] += time.perf_counter() - started
            T_b1 = self.solver.temperature()[nodes]
            started = time.perf_counter()
            book = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_b0, m_min, T_b1=T_b1, side=side,
                                    remove=self.remove)
            self.wall["skins"] += time.perf_counter() - started
            mismatch = book.E_base - (a + b * T_b1) * dt
            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
            ratio = np.abs(mismatch) / scale
            total = float(np.abs(mismatch).sum() / scale.sum())
            if attempt == 0:
                self.last_trigger = {"columns_over": int((ratio > skin_mod.INTERFACE_TOLERANCE).sum()),
                                     "worst_ratio": float(ratio.max()), "total_ratio": total}
            repeat = {"column": bool((ratio > skin_mod.INTERFACE_TOLERANCE).any()),
                      "total": total > skin_mod.INTERFACE_TOLERANCE, "never": False}[self.repeat_rule]
            if attempt == 0 and repeat:
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin = T_b1
                self.interface_repeats += 1
                continue
            break
        self.pending = mismatch
        self.last_mismatch = float(np.abs(mismatch).sum())
        E_top, E_rad = float(book.E_top.sum()), float(book.E_rad.sum())
        self.integrated_heat += E_top
        self.absorbed_heat += E_top - E_rad
        self.radiated_heat += E_rad
        self.removed_mass += float(book.removed.sum())
        self.column_removed += book.removed
        self.removed_enthalpy += float(book.removed_H.sum())
        self.drawn_mass += float(book.drawn.sum())
        self.iterations.append(res.iterations)
        self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
        self._move_base(book, m_min)
        self.last_shortfall = float(np.maximum(skin_mod.minimum_mass(mat, sk.area, self.settings) - sk.column_mass(),
                                               0.0).sum())

    def _absorb_depth(self):
        """Per corner, the depth [m] along the column's direction from the base node to where the 3D temperature falls
        to T_rigid (spec 6.6), from the node's temperature and the 3D gradient along the normal, capped at ABSORB_MAX
        of the mean surface edge per macro step."""
        mat, mm = self.material, self.moving
        if not mat.melts:
            return np.zeros(self.columns.n)
        T = self.solver.temperature()
        excess = T[mm.corners] - mat.T_rigid
        if not np.any(excess > 0.0):
            return np.zeros(self.columns.n)
        mesh = mm.mesh
        touching = np.zeros(mesh.n_nodes, dtype=bool)
        touching[mm.corners] = True
        elements = np.flatnonzero(touching[mesh.tets].any(axis=1))
        g = np.einsum("ij,ij->i", movemesh.nodal_gradient(mesh.points, mesh.tets, T, elements)[mm.corners],
                      self.columns.normal)
        depth = np.where(excess > 0.0, excess / np.maximum(g, GRADIENT_FLOOR), 0.0)
        cap = ABSORB_MAX * self.columns.edge_length.mean()
        self.absorb_capped += int((depth > cap).sum())
        return np.minimum(depth, cap)

    def _move_base(self, book, m_min):
        """Draw (already in the skins), absorb and hand down become one corner displacement; the move's actual sweep
        is then reconciled on every column's bottom (module docstring, step 4)."""
        mat, sk, mm, col = self.material, self.skins, self.moving, self.columns
        started = time.perf_counter()
        T_node = self.solver.temperature()[mm.corners]
        depth = self._absorb_depth()
        handed, handed_H = skin_mod.hand_down(sk, m_min, mat.T_rigid) if mat.melts else (np.zeros(sk.n), np.zeros(sk.n))
        inward = book.drawn / mat.rho + depth * col.area - handed / mat.rho          # m3 per corner
        u = -col.normal * (inward / col.area)[:, None] + mm.neutral_smoothing(col.normal, col.area)
        rebuilds = mm.rebuilds
        r = mm.move_with_policy(u)
        if mm.rebuilds > rebuilds:
            self.wall["rebuild"] += mm.last_rebuild.seconds
        crossed = (book.drawn > 0.0) | (depth > 0.0) | (handed > 0.0)
        T_bottom = sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)]
        self.jump_faults += int((crossed & (np.abs(T_node - T_bottom) > skin_mod.JUMP_TOLERANCE)).sum())
        sk.exchange_bottom(mat.rho * r.swept_volume - (book.drawn - handed), r.swept_enthalpy - (book.drawn_H - handed_H))
        sk.regrid()
        self.absorbed_mass += float((mat.rho * depth * col.area).sum())
        self.handed_mass += float(handed.sum())
        self.last_move = r
        self._refresh_geometry()
        self.wall["move"] += time.perf_counter() - started

    # -- what the coupled run reads ---------------------------------------------------------------------------------
    def mass(self, t=None):
        return self.material.rho * self.mesh.volume() + self.skins.mass()

    def energy(self):
        """Enthalpy held by the 3D model and the skins, and the deferred mismatch in transit to the 3D nodes [J]."""
        return self.solver.energy() + self.skins.energy() + float(self.pending.sum())

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat + removed enthalpy) / absorbed heat: zero for exact books."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy) / self.absorbed_heat

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / self.mass()))

    def surface_temperature(self):
        """The top cell of every column: what the gas sees (spec 6.3)."""
        return self.skins.T[:, 0].copy()

    def radiated_power(self):
        return float((self.emissivity * SIGMA_SB * (self.skins.T[:, 0] ** 4 - self.T_ambient ** 4) * self.surface.areas).sum())

    def field(self):
        return self.solver.temperature()

    def nose_radius(self):
        return float(self.mesh.params.get("radius_m", 0.05))

    def demised(self):
        return self.mass() < self.demise_fraction * self.mass0

    def surface_stats(self):
        T, Tt = self.solver.temperature(), self.surface_temperature()
        centre = self.mesh.points.mean(axis=0)
        i_centre = int(np.argmin(np.linalg.norm(self.mesh.points - centre, axis=1)))
        A = self.surface.areas
        return {"surface_T_max_K": float(Tt.max()), "surface_T_min_K": float(Tt.min()),
                "surface_T_mean_K": float((Tt * A).sum() / A.sum()), "T_stagnation_K": float(Tt[self.i_stag]),
                "T_back_K": float(Tt[self.i_back]), "T_centre_K": float(T[i_centre])}

    def skin_stats(self):
        """The history columns of spec 10.2 that exist at stage 2."""
        sk, mm, mat = self.skins, self.moving, self.material
        th = sk.thickness()
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        flowing = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
        liquid = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
        r = self.last_move
        return {"skin_cells": int(sk.count.sum()), "skin_thickness_max_mm": float(th.max() * 1e3),
                "skin_above_min": int((sk.column_mass() > m_min * (1.0 + 1e-9)).sum()),
                "flowing_depth_mean_mm": float(flowing.mean() * 1e3), "flowing_depth_max_mm": float(flowing.max() * 1e3),
                "liquid_depth_mean_mm": float(liquid.mean() * 1e3), "liquid_depth_max_mm": float(liquid.max() * 1e3),
                "surface_quality_min": mm.surface_quality(), "surface_fold_deg": mm.surface_fold(),
                "zone_ratio_min": float(r.worst_ratio) if r is not None else 1.0, "rebuilds": mm.rebuilds,
                "rebuild_residual_J": mm.rebuild_residual_total,
                "rebuild_residual_last_J": mm.last_rebuild.residual if mm.last_rebuild is not None else 0.0,
                "drawn_mass_kg": self.drawn_mass, "absorbed_mass_kg": self.absorbed_mass, "handed_mass_kg": self.handed_mass,
                "shortfall_kg": self.last_shortfall, "jump_faults": self.jump_faults, "absorb_capped": self.absorb_capped,
                "interface_mismatch_J": self.last_mismatch, "interface_repeats": self.interface_repeats,
                "repeat_columns_over": self.last_trigger["columns_over"],
                "mismatch_worst_ratio": self.last_trigger["worst_ratio"],
                "mismatch_total_ratio": self.last_trigger["total_ratio"],
                "removed_mass_kg": self.removed_mass, "removed_enthalpy_J": self.removed_enthalpy,
                "base_recession_max_mm": float(self._recession().max() * 1e3),
                "wall_skins_s": self.wall["skins"], "wall_3d_s": self.wall["3d"], "wall_move_s": self.wall["move"],
                "wall_rebuild_s": self.wall["rebuild"]}

    def _recession(self):
        """How far each corner's top has moved inward along its starting direction since time zero [m]."""
        return np.einsum("ij,ij->i", self._top0 - self.surface.points, self._normal0)
```

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

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -333,52 +333,59 @@
 
 
 # ---------------------------------------------------------------------------------------------------------------
-# Step 3 sub-plan 18: the melt-layer skins' flags, run names, columns and frame fields
+# Sub-plan 20: the layered skins' flags, run names, columns and frame fields
 
 def test_skin_settings_reach_the_run_name():
-    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin")
-    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
+    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics")
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics",
                               surface_model="skin") == name + "_surface-skin"
-    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", surface_model="skin",
-                              skin_suffix="_skin-dev", seed=7) == name + "_surface-skin_skin-dev_seed-7"
-    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", None,
+    melt = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "instant")
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "instant", surface_model="skin",
+                              skin_suffix="_skin-dev", seed=7) == melt + "_surface-skin_skin-dev_seed-7"
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", None,
                               surface_model="skin") == cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table",
-                                                                          "none", "physics", None)   # no melt, no skins
-    assert cli.skin_name_suffix(0.4, None, 10.0, "production") == "" and cli.skin_name_suffix(0.4, 10.0, 10.0, "production") == ""
-    assert cli.skin_name_suffix(0.4, None, 10.0, "dev") == "_skin-dev"
-    assert cli.skin_name_suffix(0.8, 5.0, 20.0, "production") == "_skin0.8mm_cell5um_sub20ms"
+                                                                          "none", None)        # no thermal model, no skins
+    assert cli.skin_name_suffix(0.4, None, 1.10, 10.0, "production", 15.0) == ""
+    assert cli.skin_name_suffix(0.4, 10.0, 1.10, 10.0, "production", 15.0) == ""
+    assert cli.skin_name_suffix(0.4, None, 1.10, 10.0, "dev", 15.0) == "_skin-dev"
+    assert (cli.skin_name_suffix(0.8, 5.0, 1.2, 20.0, "production", 10.0)
+            == "_skinmin0.8mm_cell5um_growth1.2_sub20ms_deform10mm")
 
 
 def test_a_skin_run_writes_its_columns_frames_and_settings(tmp_path):
-    """3 s from 71 km with a warm body and skins (dev preset), the 0.5 s step throughout (`--dt-continuum off`, as the
-    short melting run passes its own step): exit 0, the run name, the skin columns, the frame fields."""
-    pytest.importorskip("cantera")
+    """1 s from 71 km with a warm body and the layered skins (dev preset, liquid removed as it forms): exit 0, the run
+    name, the settings, the skin columns and results, the frame fields."""
     import glob
     import pyvista as pv
-    argv = [a for a in MELT if a not in ("--heating", "sesam")] + [
-        "--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "800", "--t-max", "3",
-        "--dt-continuum", "off", "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "2", "--outdir", str(tmp_path)]
+    argv = MELT + ["--removal", "instant", "--altitude", "71", "--velocity", "7.24", "--temperature", "880", "--t-max", "1",
+                   "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "1", "--outdir", str(tmp_path)]
     assert cli.main(argv) == 0
     (path,) = glob.glob(str(tmp_path / "*.json"))
     doc = json.load(open(path))
-    assert doc["run_name"].endswith("_melt-girin_dtcontinuum-off_surface-skin_skin-dev")
+    assert doc["run_name"].endswith("_melt-instant_surface-skin_skin-dev")
     s, r = doc["settings"], doc["results"]
-    assert s["surface_model"] == "skin" and s["skin"] == {"thickness_mm": 0.4, "cell_um": 20.0, "substep_ms": 10.0, "preset": "dev"}
-    assert r["skins_created"] > 0 and r["skin_drawn_mass_kg"] >= 0.0 and abs(r["melt_energy_balance_residual"]) < 1e-6
+    assert s["surface_model"] == "skin" and s["prism_layers"] == 0
+    assert s["skin"] == {"min_thickness_mm": 0.4, "cell_um": 20.0, "growth": 1.1, "substep_ms": 10.0, "preset": "dev",
+                         "deform_depth_mm": 15.0, "n_columns": r["skin_columns"]}
+    assert r["skin_columns"] == s["n_patches"] and abs(r["energy_balance_residual"]) < 1e-8
     rows = list(csv.DictReader(open(path[:-5] + ".csv")))
-    assert all(k in rows[0] for k in coupled.SKIN_COLUMNS) and float(rows[-1]["skin_count"]) > 0.0
+    assert all(k in rows[0] for k in coupled.SKINBODY_COLUMNS) and float(rows[-1]["skin_cells"]) > 0.0
     surf = sorted(glob.glob(str(tmp_path / "**" / "surface_*.vtp"), recursive=True))
-    assert surf and all(k in pv.read(surf[-1]).cell_data for k in ("skin_T_top", "skin_liquid_depth", "skin_nonrigid_depth",
-                                                                    "skin_thickness"))
+    assert len(surf) == 3 and all(k in pv.read(surf[-1]).point_data for k in ("T_top", "q_conv", "skin_thickness", "skin_cells",
+                                                                             "liquid_depth", "flowing_depth", "base_position"))
 
 
 @pytest.mark.parametrize("argv", [
     MELT + ["--surface-model", "magic"],
-    MELT + ["--surface-model", "skin", "--skin-thickness", "0"],
-    MELT + ["--surface-model", "skin", "--skin-substep", "-1"],
-    MELT + ["--surface-model", "skin", "--skin-cell", "300"],          # more than half of the 0.4 mm skin
-    MELT + ["--surface-model", "skin", "--skin-preset", "fine"],
-    MELT + ["--surface-model", "skin", "--removal", "instant"],
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-min-thickness", "0"],
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-substep", "-1"],
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-cell", "300"],      # more than half of 0.4 mm
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-preset", "fine"],
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-growth", "0.9"],
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--deform-depth", "0"],
+    MELT + ["--surface-model", "skin"],                                     # Girin's removal needs flow and spraying
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--animate"],
+    ["run", "--surface-model", "skin"],                                     # no thermal model
 ])
 def test_bad_skin_arguments_exit_2(argv, tmp_path):
     with pytest.raises(SystemExit) as exc:
```

```diff
--- a/tests/test_reentry_model_melting_skin.py
+++ b/tests/test_reentry_model_melting_skin.py
@@ -18,7 +18,7 @@
 EXACT = dict(linear_solver="direct")
 
 
-def make_body(the_mesh, surface_model="skin", name="AA7075_scheil", solver="skfem", solver_options=None, **settings):
+def make_body(the_mesh, surface_model="skin1d", name="AA7075_scheil", solver="skfem", solver_options=None, **settings):
     mat = material.Material.from_drama_json(name)
     m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
     return body.MeltingBody(m, mat, thermal.thermal_solver(solver, **(solver_options or {})), MASS_100MM,
@@ -52,13 +52,13 @@
     with pytest.raises(ValueError):
         body.MeltSettings(surface_model="magic")
     with pytest.raises(ValueError):
-        body.MeltSettings(surface_model="skin", removal="instant")
+        body.MeltSettings(surface_model="skin1d", removal="instant")
     assert body.MeltSettings().surface_model == "elements"
 
 
 def test_identical_to_the_element_model_until_the_first_skin(coarse_sphere_mesh):
     """Seeded alike before every step (the AMG set-up draws on numpy's generator, as the model's --seed handles)."""
-    be, bs = make_body(coarse_sphere_mesh, "elements"), make_body(coarse_sphere_mesh, "skin")
+    be, bs = make_body(coarse_sphere_mesh, "elements"), make_body(coarse_sphere_mesh, "skin1d")
     assert be.skins is None and bs.skins is not None and bs.skins.n == 0
     t, compared = 0.0, 0
     while t < 120.0:
```

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -321,7 +321,7 @@
         np.random.seed(12345)
         m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
         b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_scheil"), thermal.thermal_solver(name), MASS_100MM,
-                             settings=body.MeltSettings(surface_model="skin",
+                             settings=body.MeltSettings(surface_model="skin1d",
                                                         skin=skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], substep=0.05)))
         b.solver.set_temperature(800.0)
         b.energy0 = b.energy()
```

- [ ] **Step 2: Run them to see them fail.** Expected: the new CLI tests fail (`--skin-min-thickness` unrecognised and
  `skin_name_suffix` taking four arguments); stage 1's tests fail on `surface_model must be one of ('elements', 'skin')`.

- [ ] **Step 3: Apply the changes**

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -151,7 +151,8 @@
                                  # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
 NONRIGID_MAX_CROSSINGS = 64      # elements the non-rigid depth's line may cross before it stops (amendment of 2026-10-06)
 REMOVAL_NAMES = ("girin", "instant")
-SURFACE_MODEL_NAMES = ("elements", "skin")    # skin: a one-dimensional skin under every melting patch (sub-plan 18)
+SURFACE_MODEL_NAMES = ("elements", "skin1d")  # skin1d: a one-dimensional skin under every melting patch (sub-plan 18, failed
+                                              # its box test; retired in sub-plan 22). The layered skins are skinbody.SkinBody
 
 
 @dataclass
@@ -170,9 +171,9 @@
     rigid_substrate: bool = True      # the thin branch needs a rigid substrate: the regime test reads the film plus the slurry
                                       # (above 50 % liquid) beneath it, and off Girin's closure a film on slurry does not
                                       # spray (amendment of 2026-10-06)
-    surface_model: str = "elements"   # skin: a one-dimensional skin of fine cells under every melting patch resolves the
+    surface_model: str = "elements"   # skin1d: a one-dimensional skin of fine cells under every melting patch resolves the
                                       # melt layer and feeds the film continuously (sub-plan 18); elements: as before
-    skin: object = None               # skin.SkinSettings for surface_model "skin"; None: its defaults
+    skin: object = None               # skin1d.SkinSettings for surface_model "skin1d"; None: its defaults
 
     def __post_init__(self):
         if self.surface_model not in SURFACE_MODEL_NAMES:
@@ -185,8 +186,8 @@
             raise ValueError("feed_depth must be one of {}, got {!r}".format(FEED_DEPTH_NAMES, self.feed_depth))
         if not 0.0 < self.demise_fraction < 1.0:
             raise ValueError("demise_fraction must be within (0, 1)")
-        if self.surface_model == "skin" and self.removal == "instant":
-            raise ValueError("surface_model 'skin' needs removal 'girin'")
+        if self.surface_model == "skin1d" and self.removal == "instant":
+            raise ValueError("surface_model 'skin1d' needs removal 'girin'")
 
 
 def fit_sphere(points):
@@ -280,8 +281,8 @@
         self.consumed = False
         self._refresh_geometry()
         self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation
-        self.skins = None                                          # sub-plan 18: the skins, with --surface-model skin
-        if self.settings.surface_model == "skin" and self.settings.removal != "instant":
+        self.skins = None                                          # sub-plan 18: the skins, with surface_model "skin1d"
+        if self.settings.surface_model == "skin1d" and self.settings.removal != "instant":
             from . import skin1d as skin_mod
             self.skins = skin_mod.SkinField(material, self.settings.skin or skin_mod.SkinSettings())
         self._draw_m, self._draw_H = np.zeros(mesh.n_elements), np.zeros(mesh.n_elements)   # metal the skins drew this step
```

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -19,7 +19,7 @@
 
 import numpy as np
 
-from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, skin1d as skin, spray, surface_flow, thermal, viz
+from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, skin, skinbody, spray, surface_flow, thermal, viz
 from . import trajectory as tj
 from .earth import GRAVITY_MODELS
 from .thermal import MissingBackend
@@ -47,6 +47,10 @@
 # seeds, the droplet count within 2.7 % of the 0.00625 s step's (plan facts 92-94) -- while the default 0.5 s step makes
 # six times fewer droplets. Before the switch the default step serves.
 DEFAULT_DT_CONTINUUM = 0.0125
+# The surface models a run can choose (sub-plan 20): the element model, or the layered skins on the corners of a moving
+# mesh (skinbody.SkinBody). The default stays the element model until sub-plan 22 gives the skins their flow and spraying;
+# sub-plan 18's one-dimensional skins are reachable only from Python (body.MeltSettings(surface_model="skin1d")).
+SURFACE_MODELS = ("elements", "skin")
 
 
 def parse_epoch(text):
@@ -96,9 +100,9 @@
     run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
     (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
     `_rigidsubstrate-off`, one whose continuum step is not its removal mode's default (amendment of 2026-10-07) in
-    `_dtcontinuum-off` or `_dtcontinuum-<s>`, one with skins (sub-plan 18) in `_surface-skin` followed by `skin_suffix`
-    (the skins' non-default settings, `skin_name_suffix`), and any run with a seed other than DEFAULT_SEED (amendment of
-    2026-10-05) in `_seed-<n>`; the defaults keep the old name."""
+    `_dtcontinuum-off` or `_dtcontinuum-<s>`, a thermal run with the layered skins (sub-plan 20) in `_surface-skin`
+    followed by `skin_suffix` (the skins' non-default settings, `skin_name_suffix`), and any run with a seed other than
+    DEFAULT_SEED (amendment of 2026-10-05) in `_seed-<n>`; the defaults keep the old name."""
     dtc_default = DEFAULT_DT_CONTINUUM if melt == "girin" else None
     dtc = dtc_default if dt_continuum == "default" else dt_continuum
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
@@ -107,21 +111,26 @@
             + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else "")
             + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
             + ("" if not melt or dtc == dtc_default else "_dtcontinuum-off" if dtc is None else "_dtcontinuum-{:g}".format(dtc))
-            + ("_surface-skin" + skin_suffix if melt and surface_model == "skin" else "")
+            + ("_surface-skin" + skin_suffix if heating_name and surface_model == "skin" else "")
             + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))
 
 
-def skin_name_suffix(thickness_mm, cell_um, substep_ms, preset):
-    """The run-name suffix of a skin run's non-default settings (sub-plan 18): `_skin-dev` for the dev preset,
-    `_skin<mm>mm` for another thickness, `_cell<um>um` for an explicit cell other than the preset's, `_sub<ms>ms` for
-    another sub-step; empty for the defaults."""
+def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm):
+    """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
+    dev preset, `_skinmin<mm>mm` for another minimum thickness, `_cell<um>um` for an explicit cell other than the
+    preset's, `_growth<ratio>`, `_sub<ms>ms` and `_deform<mm>mm` for another growth, sub-step or deforming depth; empty
+    for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
-    if thickness_mm != 0.4:
-        out += "_skin{:g}mm".format(thickness_mm)
+    if min_thickness_mm != 0.4:
+        out += "_skinmin{:g}mm".format(min_thickness_mm)
     if cell_um is not None and abs(cell_um * 1e-6 - skin.SKIN_PRESETS[preset]) > 1e-12:
         out += "_cell{:g}um".format(cell_um)
+    if growth != 1.10:
+        out += "_growth{:g}".format(growth)
     if substep_ms != 10.0:
         out += "_sub{:g}ms".format(substep_ms)
+    if deform_mm != 15.0:
+        out += "_deform{:g}mm".format(deform_mm)
     return out
 
 
@@ -250,15 +259,21 @@
                          "with the film plus everything more than half liquid beneath it (slurry included), and where there "
                          "is no conjugate depth a film on slurry does not spray (default on). off: the regime test reads the "
                          "film plus the fully molten material, as in every run before the amendment of 2026-10-06")
-    me.add_argument("--surface-model", choices=body.SURFACE_MODEL_NAMES, default="elements",
-                    help="skin: a one-dimensional skin of fine cells under every melting patch resolves the melt layer and "
-                         "feeds the film continuously, coupled to the 3D conduction through an implicit interface (Step 3 "
-                         "sub-plan 18). elements: the surface elements melt and feed the film, as in every run before it (default)")
-    me.add_argument("--skin-thickness", type=float, default=0.4, help="skin target thickness [mm] (default %(default)s)")
-    me.add_argument("--skin-cell", type=float, default=None, help="skin cell size [um] (default: the preset's)")
+    me.add_argument("--surface-model", choices=SURFACE_MODELS, default="elements",
+                    help="skin: a layered skin on every corner of the surface holds the outer metal in fine cells above a "
+                         "moving 3D mesh (layered-skin design, sub-plan 20: heat only, so with --melt on it needs --removal "
+                         "instant, liquid leaving as it forms). elements: the surface elements melt and feed the film, as "
+                         "in every run before it (default)")
+    me.add_argument("--skin-min-thickness", type=float, default=0.4,
+                    help="the skins' minimum thickness [mm of solid], also the fine zone's depth (default %(default)s)")
+    me.add_argument("--skin-cell", type=float, default=None, help="skin fine cell [um] (default: the preset's)")
+    me.add_argument("--skin-growth", type=float, default=1.10,
+                    help="each skin cell below the fine zone this much thicker than the one above (default %(default)s)")
     me.add_argument("--skin-substep", type=float, default=10.0, help="skin sub-step [ms] (default %(default)s)")
     me.add_argument("--skin-preset", choices=tuple(skin.SKIN_PRESETS), default="production",
                     help="skin cell-size preset: production 10 um, dev 20 um (default %(default)s)")
+    me.add_argument("--deform-depth", type=float, default=15.0,
+                    help="depth of the moving mesh's deforming zone [mm] (default %(default)s, the dense band's)")
     me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                     help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                          "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
@@ -301,10 +316,12 @@
     """(ThermalBody or MeltingBody, HeatingModel, settings-provenance dict) for --thermal fem [--melt on]."""
     radius = settings.diameter / 2.0
     melting = args.melt == "on"
+    skins = args.surface_model == "skin"
     h_surface, h_core = args.h_surface * 1e-3 * args.mesh_size, args.h_core * 1e-3 * args.mesh_size
-    layers = args.prism_layers if args.prism_layers is not None else (mesh.DEFAULT_LAYERS if melting else 0)
-    band = mesh.DEFAULT_BAND if melting else 0.0    # non-melting runs keep the Step 2 field: the committed Step 2
-                                                    # verification numbers were measured on it (fact 12)
+    layers = args.prism_layers if args.prism_layers is not None else (mesh.DEFAULT_LAYERS if melting and not skins else 0)
+    band = mesh.DEFAULT_BAND if melting or skins else 0.0    # non-melting runs keep the Step 2 field: the committed Step 2
+                                                    # verification numbers were measured on it (fact 12); the skins'
+                                                    # moving mesh is rebuilt with the dense band's size field
     the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, layers=layers, layer_thickness=args.layer_thickness * 1e-3,
                                band=band)
     # the melting default is Scheil's curve (decided 2026-09-27, made the default on 2026-10-07); the linear range stays
@@ -314,7 +331,11 @@
     if args.k_scale != 1.0:
         mat.k_table = mat.k_table * args.k_scale                  # verification device (near-isothermal body), not physical
     solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=not args.consistent_mass)
-    if melting:
+    if skins:
+        the_body = skinbody.SkinBody(the_mesh, mat, solver, mass, skin_settings(args), T0=args.temperature,
+                                     emissivity=args.emissivity, T_ambient=args.t_ambient, deform_depth=args.deform_depth * 1e-3,
+                                     remove=melting, demise_fraction=args.demise_fraction)
+    elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical,
                                        rt_spray=args.rt_spray == "on")
@@ -322,8 +343,7 @@
         melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                           particles=args.particles, size_feedback=size_feedback,
                                           deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on",
-                                          rigid_substrate=args.rigid_substrate == "on", surface_model=args.surface_model,
-                                          skin=skin_settings(args))
+                                          rigid_substrate=args.rigid_substrate == "on")
         the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                     emissivity=args.emissivity, T_ambient=args.t_ambient)
     else:
@@ -338,16 +358,20 @@
             "h_surface_mm": h_surface * 1e3, "h_core_mm": h_core * 1e3, "prism_layers": layers, "layer_thickness_mm": args.layer_thickness,
             "n_nodes": the_mesh.n_nodes, "n_elements": the_mesh.n_elements, "n_patches": the_body.surface.n_patches,
             "thermal_solver": args.thermal_solver, "linear_solver": args.linear_solver, "lumped_mass": not args.consistent_mass,
-            "k_scale": args.k_scale, "melt": args.melt}
+            "k_scale": args.k_scale, "melt": args.melt, "surface_model": args.surface_model}
+    if skins:
+        info["skin"] = {"min_thickness_mm": args.skin_min_thickness, "cell_um": skin_settings(args).cell * 1e6,
+                        "growth": args.skin_growth, "substep_ms": args.skin_substep, "preset": args.skin_preset,
+                        "deform_depth_mm": args.deform_depth, "n_columns": the_body.columns.n}
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
-    if melting:
+    if melting and skins:
+        info.update({"removal": args.removal, "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus,
+                     "latent_heat_Jkg": mat.latent_heat})
+    elif melting:
         info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff,
                      "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
-                     "surface_model": args.surface_model,
-                     "skin": {"thickness_mm": args.skin_thickness, "cell_um": skin_settings(args).cell * 1e6,
-                              "substep_ms": args.skin_substep, "preset": args.skin_preset},
                      "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                      "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                      "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
@@ -357,9 +381,10 @@
 
 
 def skin_settings(args):
-    """The skins' settings from the command line (sub-plan 18)."""
+    """The layered skins' settings from the command line (sub-plan 20)."""
     cell = args.skin_cell * 1e-6 if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset]
-    return skin.SkinSettings(thickness=args.skin_thickness * 1e-3, cell=cell, substep=args.skin_substep * 1e-3)
+    return skin.SkinSettings(min_thickness=args.skin_min_thickness * 1e-3, cell=cell, growth=args.skin_growth,
+                             substep=args.skin_substep * 1e-3)
 
 
 def cmd_run(args, parser):
@@ -399,15 +424,22 @@
     dt_continuum = continuum_step(args)
     if dt_continuum is not None and dt_continuum > args.dt:
         parser.error("--dt-continuum must be no longer than --dt")
-    for label, value in (("--skin-thickness", args.skin_thickness), ("--skin-substep", args.skin_substep),
-                         ("--skin-cell", 1.0 if args.skin_cell is None else args.skin_cell)):
+    for label, value in (("--skin-min-thickness", args.skin_min_thickness), ("--skin-substep", args.skin_substep),
+                         ("--skin-cell", 1.0 if args.skin_cell is None else args.skin_cell),
+                         ("--deform-depth", args.deform_depth)):
         if value <= 0.0:
             parser.error("{} must be > 0".format(label))
+    if args.skin_growth < 1.0:
+        parser.error("--skin-growth must be at least 1")
     cell_um = args.skin_cell if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset] * 1e6
-    if 2.0 * cell_um * 1e-3 > args.skin_thickness * (1.0 + 1e-12):
-        parser.error("--skin-cell must be at most half of --skin-thickness")
-    if args.surface_model == "skin" and args.removal != "girin":
-        parser.error("--surface-model skin needs --removal girin")
+    if 2.0 * cell_um * 1e-3 > args.skin_min_thickness * (1.0 + 1e-12):
+        parser.error("--skin-cell must be at most half of --skin-min-thickness")
+    if args.surface_model == "skin" and args.thermal != "fem":
+        parser.error("--surface-model skin needs --thermal fem")
+    if args.surface_model == "skin" and args.melt == "on" and args.removal != "instant":
+        parser.error("--surface-model skin with --melt on needs --removal instant (sub-plan 20: no flow or spraying yet)")
+    if args.surface_model == "skin" and (args.animate or args.stills):
+        parser.error("--animate/--stills read the element model's frames; with --surface-model skin use --frames-every")
 
     initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                               flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
@@ -422,8 +454,8 @@
                                        deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on", seed=args.seed,
                                        rigid_substrate=args.rigid_substrate == "on", dt_continuum=dt_continuum,
                                        surface_model=args.surface_model,
-                                       skin_suffix=skin_name_suffix(args.skin_thickness, args.skin_cell, args.skin_substep,
-                                                                    args.skin_preset))
+                                       skin_suffix=skin_name_suffix(args.skin_min_thickness, args.skin_cell, args.skin_growth,
+                                                                    args.skin_substep, args.skin_preset, args.deform_depth))
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
@@ -464,7 +496,7 @@
         "provenance": prov,
         "files": {"csv": os.path.abspath(csv_path)},
     }
-    melting = args.thermal == "fem" and args.melt == "on"
+    melting = args.thermal == "fem" and args.melt == "on" and args.surface_model != "skin"
     if melting:
         doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()} if args.particles
                             else {})
@@ -502,6 +534,11 @@
                   "energy balance residual {:.1e}, {} macro steps, {:.1f} Newton iterations/step".format(
                       res["peak_surface_T_K"], res["time_of_peak_surface_T_s"], res["peak_mean_T_K"], res["integrated_heat_J"],
                       res["energy_balance_residual"], res["n_macro_steps"], res["mean_newton_iterations"]))
+        if args.thermal == "fem" and args.surface_model == "skin":
+            print("  skins: {} columns, removed {:.4g} kg, drawn {:.4g} kg, {} rebuilds, {} interface repeats, {} base faults; "
+                  "wall time {:.0f} s in the skins, {:.0f} s in the 3D solve".format(
+                      res["skin_columns"], res["removed_mass_kg"], res["drawn_mass_kg"], res["rebuilds"],
+                      res["interface_repeats"], res["jump_faults"], res["wall_skins_s"], res["wall_3d_s"]))
         if melting:
             print("  melt: onset {} km, spraying onset {} km, demise {} km at t = {} s; sprayed {:.4f} kg of {:.4f}, film left {:.4f} kg, "
                   "{:.3g} droplets (median r {} um), {} source rows, melt balance {:.1e}".format(
```

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -37,6 +37,14 @@
                 "skin_frozen_mass_kg", "skin_returned_mass_kg", "skin_full_events", "skin_thick_feeding_count",
                 "skin_thick_feeding_max_mm", "skin_thick_feed_events"]   # sub-plan 18
 MELT_COLUMNS = MELT_COLUMNS + SKIN_COLUMNS
+# the layered skins on the corner columns (sub-plan 20): skinbody.SkinBody.skin_stats(), every macro step
+SKINBODY_COLUMNS = ["skin_cells", "skin_thickness_max_mm", "skin_above_min", "flowing_depth_mean_mm", "flowing_depth_max_mm",
+                    "liquid_depth_mean_mm", "liquid_depth_max_mm", "surface_quality_min", "surface_fold_deg", "zone_ratio_min",
+                    "rebuilds", "rebuild_residual_J", "rebuild_residual_last_J", "drawn_mass_kg", "absorbed_mass_kg",
+                    "handed_mass_kg", "shortfall_kg", "jump_faults", "absorb_capped", "interface_mismatch_J",
+                    "interface_repeats", "repeat_columns_over", "mismatch_worst_ratio", "mismatch_total_ratio",
+                    "removed_mass_kg", "removed_enthalpy_J", "base_recession_max_mm",
+                    "wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s"]
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -79,11 +87,17 @@
         r.update(body.surface_stats())
         if self.melting:
             r.update(body.melt_stats())
+        if self.skins:
+            r.update(body.skin_stats())
         return r
 
     @property
     def melting(self):
         return hasattr(self.body, "melt_stats")
+
+    @property
+    def skins(self):
+        return hasattr(self.body, "skin_stats")
 
     def run(self):
         started = time.perf_counter()
@@ -104,7 +118,7 @@
             rows.append(self.row(sim.t, sim.y, a, loads))
             states.append(sim.y.copy())
             self._frame(step, sim.t, loads)
-            if self.melting and body.demised():
+            if (self.melting or self.skins) and body.demised():
                 sim.end_reason = "demise"
         self._collection()
         columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
@@ -113,6 +127,8 @@
         history.results.update(self.thermal_results(history))
         if self.melting:
             history.results.update(self.melt_results(history))
+        if self.skins:
+            history.results.update(self.skin_results(history))
         if s.dt_continuum is not None:
             t_sw, kn_sw, h_sw = self.switched_at or (None, None, None)
             history.results.update({"dt_continuum_s": s.dt_continuum, "dt_switch_time_s": t_sw, "dt_switch_kn": kn_sw,
@@ -156,6 +172,20 @@
                 "skin_returned_mass_kg": body.skin_returned_mass, "skin_thick_feed_events": body.skin_thick_feed_events,
                 "melt_energy_balance_residual": body.energy_balance_residual()}
 
+    def skin_results(self, history):
+        """The layered skins' totals (sub-plan 20)."""
+        body, last = self.body, history.columns
+        out = {"skin_columns": body.columns.n, "initial_mass_kg": body.mass0, "final_mass_kg": float(last["mass_kg"][-1]),
+               "demise_time_s": float(last["time_s"][-1]) if history.end_reason == "demise" else None,
+               "surface_fold_max_deg": float(np.max(last["surface_fold_deg"])),
+               "surface_quality_min": float(np.min(last["surface_quality_min"])),
+               "zone_ratio_min": float(np.min(last["zone_ratio_min"]))}
+        out.update({k: float(last[k][-1]) for k in ("removed_mass_kg", "removed_enthalpy_J", "drawn_mass_kg", "absorbed_mass_kg",
+                                                  "handed_mass_kg", "rebuild_residual_J", "wall_skins_s", "wall_3d_s",
+                                                  "wall_move_s", "wall_rebuild_s")})
+        out.update({k: int(last[k][-1]) for k in ("rebuilds", "jump_faults", "absorb_capped", "interface_repeats")})
+        return out
+
     def thermal_results(self, history):
         c, body = history.columns, self.body
         i_surf, i_mean, i_q = int(np.argmax(c["surface_T_max_K"])), int(np.argmax(c["temperature_K"])), int(np.argmax(c["convective_heat_W"]))
@@ -210,6 +240,8 @@
     import pyvista as pv
     from .thermal import SIGMA_SB
     os.makedirs(output_dir, exist_ok=True)
+    if hasattr(body, "skin_stats"):
+        return write_skin_frame(output_dir, k, body, loads)
     mesh, surface, T = body.mesh, body.surface, body.field()
     melting = hasattr(body, "melt_stats")
     tets = mesh.tets[mesh.active]
@@ -261,6 +293,28 @@
     poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
 
 
+def write_skin_frame(output_dir, k, body, loads):
+    """A frame of the layered-skin body (sub-plan 20, spec 10.2): field_<k>.vtu, the 3D model's nodal T on its current
+    (moving) mesh; surface_<k>.vtp, the top surface with per-corner point data -- top temperature, convective flux,
+    liquid and flowing depths, skin thickness and cells, and the base's position (the corner on the 3D mesh)."""
+    import pyvista as pv
+    mesh, sk, mat = body.mesh, body.skins, body.material
+    cells = np.hstack([np.full((mesh.n_elements, 1), 4), mesh.tets]).ravel()
+    grid = pv.UnstructuredGrid(cells, np.full(mesh.n_elements, pv.CellType.TETRA), mesh.points)
+    grid.point_data["T"] = body.field()
+    grid.save(os.path.join(output_dir, "field_{}.vtu".format(k)))
+    surface = body.surface
+    poly = pv.PolyData(surface.points, np.hstack([np.full((len(surface.faces), 1), 3), surface.faces]).ravel())
+    poly.point_data["T_top"] = sk.T[:, 0]
+    poly.point_data["q_conv"] = loads.q_conv
+    poly.point_data["skin_thickness"] = sk.thickness()
+    poly.point_data["skin_cells"] = sk.count.astype(float)
+    poly.point_data["liquid_depth"] = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
+    poly.point_data["flowing_depth"] = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
+    poly.point_data["base_position"] = body.columns.points
+    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
+
+
 def write_particles(run_dir, body, history, window=10.0):
     """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
     size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths."""
```

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

```python
"""Sub-plan 20 Task 9: the independent one-dimensional reference of level 2 (layered-skin design section 11.1).

The centre of the level-2 plate (and the middle of sub-plan 18's fine box) is one-dimensional: a slab of `length`
heated at its front by a constant flux, adiabatic behind, no radiation, liquid removed as it forms. This solves it on a
fixed grid of `cell` cells with its own code -- backward Euler in the enthalpy, Newton in the temperature with a banded
solve, and the contiguous front cells at or above T_feed removed after every step -- sharing nothing with
`reentry_model.skin` but the material.

    python analysis/ablation_1d.py [--cell-um 1] [--dt-ms 0.5] [--q 6e6] [--t0 700] [--length-mm 8] [--t-end 2] [--out PATH]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.linalg import solve_banded

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import material  # noqa: E402


def run(mat, q=6.0e6, T0=700.0, length=8.0e-3, cell=1.0e-6, dt=5.0e-4, t_end=2.0, every=0.05):
    """Recession [m] of the slab's front at multiples of `every` [s], as (times, recession) arrays."""
    n = int(round(length / cell))
    T = np.full(n, T0)
    front = 0                                            # index of the first cell still in the slab
    m = mat.rho * cell                                   # kg per m2 per cell
    times, recession = [0.0], [0.0]
    steps = int(round(t_end / dt))
    record = max(1, int(round(every / dt)))
    for step in range(1, steps + 1):
        Tc = T[front:]
        k_old = Tc.size
        if k_old == 0:
            break
        H_old = m * mat.enthalpy(Tc)
        Tk = Tc.copy()
        for _ in range(60):
            k = mat.k(Tk)
            G = 1.0 / (0.5 * cell / k[:-1] + 0.5 * cell / k[1:])          # W/m2/K between neighbours
            flux = G * (Tk[:-1] - Tk[1:])
            R = (m * mat.enthalpy(Tk) - H_old) / dt
            R[0] -= q
            R[:-1] += flux
            R[1:] -= flux
            diag = m * mat.cp_eff(Tk) / dt
            diag[:-1] += G
            diag[1:] += G
            ab = np.zeros((3, k_old))
            ab[0, 1:], ab[1], ab[2, :-1] = -G, diag, -G
            dT = solve_banded((1, 1), ab, -R)
            Tk = mat.temperature_from_enthalpy(mat.enthalpy(Tk) + mat.cp_eff(Tk) * dT)   # enthalpy-consistent update
            if np.abs(dT).max() < 1e-7:
                break
        else:
            raise RuntimeError("the 1D reference did not converge at step {}".format(step))
        T[front:] = Tk
        liquid = np.cumprod(Tk >= mat.T_feed)
        front += int(liquid.sum())
        if step % record == 0:
            times.append(step * dt)
            recession.append(front * cell)
    return np.array(times), np.array(recession)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-um", type=float, default=1.0)
    p.add_argument("--dt-ms", type=float, default=0.5)
    p.add_argument("--q", type=float, default=6.0e6)
    p.add_argument("--t0", type=float, default=700.0)
    p.add_argument("--length-mm", type=float, default=8.0)
    p.add_argument("--t-end", type=float, default=2.0)
    p.add_argument("--out", default="ablation_1d.json")
    a = p.parse_args()
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    t, r = run(mat, a.q, a.t0, a.length_mm * 1e-3, a.cell_um * 1e-6, a.dt_ms * 1e-3, a.t_end)
    doc = {"q_Wm2": a.q, "T0_K": a.t0, "length_mm": a.length_mm, "cell_um": a.cell_um, "dt_ms": a.dt_ms,
           "material": mat.name, "t_s": t.tolist(), "recession_mm": (r * 1e3).tolist()}
    json.dump(doc, open(a.out, "w"), indent=1)
    for t_ in (1.0, 1.5, 2.0):
        i = int(np.argmin(np.abs(t - t_)))
        print("t = {:.2f} s: recession {:.4f} mm".format(t[i], r[i] * 1e3))


if __name__ == "__main__":
    main()
```

  Run `"$PY" analysis/ablation_1d.py --out ../measure/ablation_1d_1um.json` (about 12 s). Expected (measured while
  planning): 0.7890, 2.7120 and 5.2190 mm at 1.0, 1.5 and 2.0 s; at 5 µm and 1 ms 0.7900, 2.7100, 5.2150.

- [ ] **Step 2: The fine box's series.** Save as `../extract_fine_box.py` beside the work copy and run it from the work
  copy's root (the run file is git-ignored, so it is read from the main checkout's
  `prototype/work-2026-10-08-skin/analysis/box/runs/q6p__fine.json`):

```python
"""Write data/reference_values/box_q6p_fine.json from sub-plan 18's fine-box run (sub-plan 20 Task 9). Run from the
work copy's root: python ../extract_fine_box.py (or with the run file's path as the first argument)."""
import json
import os
import sys

RUN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "work-2026-10-08-skin", "analysis", "box", "runs",
                   "q6p__fine.json")
src = sys.argv[1] if len(sys.argv) > 1 else RUN
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join("data", "reference_values", "box_q6p_fine.json")
d = json.load(open(src))
doc = {"description": "The stage-1 fine box of sub-plan 18 (q6p): an 8 x 2 x 2 mm block of 0.1 mm elements at 700 K, AA7075_scheil, "
       "no radiation, 6 MW/m2 on its x = 0 face distributed by projection, liquid removed as it forms (element model, "
       "--removal instant), 0.005 s steps to 2 s. r_total_mm is the total removed mass over the 4 mm2 cross-section. "
       "Measured 2026-10-09 (prototype/work-2026-10-08-skin/analysis/box/runs/q6p__fine.json); a reference for the "
       "level-2 plate of the layered-skin design beside the one-dimensional 1 um reference (sub-plan 20 Task 9).",
       "spec": d["spec"], "t_s": d["series"]["t"], "r_total_mm": d["series"]["r_total_mm"]}
json.dump(doc, open(out, "w"), indent=0)
```

  The fine box's recession: 0.8045, 2.7290 and 5.2060 mm at 1.0, 1.5 and 2.0 s, within 2 % of the reference.

- [ ] **Step 3: The plate.** Create `analysis/skin_plate.py`:

```python
"""Sub-plan 20 Task 9: level 2 of the layered-skin design's verification ladder (section 11.1, amended 2026-10-12), flow
off: the plate.

The q6p box of sub-plan 18 cannot be followed by corner columns: its 90-degree edges fold past the fixed surface
triangulation's limit (1.4 s on 2 mm elements, 0.9 s on 1 mm, found while planning), so the level-2 test is a plate
whose heated region tapers to nothing well inside its edges (decided 2026-10-12, Asha). The plate is 8 mm deep (the
box's length) and 110 x 110 mm across, of 2 mm elements at 700 K, Scheil material, no radiation. Its x = 0 face is
heated at 6 MW/m2 within 32 mm of its centre, the flux falling to zero at 45 mm along a half cosine, every column
receiving what the projection of its top dual cell on the face receives (`projected_flux`); liquid is removed as it
forms. The centre is one-dimensional while sideways conduction from the taper cannot reach it (the diffusion length
2 sqrt(alpha t) is about 19 mm at 2 s), so its recession -- read on the front corners within 3 mm of
the centre, as the top of each column with all its metal solid (`solid_top`) -- is compared at 1.0, 1.5 and 2.0 s with
the one-dimensional 1 um reference (`analysis/ablation_1d.py`) and the stage-1 fine box
(`data/reference_values/box_q6p_fine.json`). Criterion: within 5 % of both; mass closed to 1e-12 and energy to 1e-8.
The crater's wall is where the surface folds under heating-driven recession: the fold, the triangle and element
quality, rebuilds, faults, repeats and the wall time are recorded every step.

    python analysis/skin_plate.py [--dt 0.05] [--preset dev|production] [--h-mm 2] [--t-end 2] [--out PATH.json]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import heating, material, mesh, skin, skinbody  # noqa: E402
from reentry_model.thermal import thermal_solver  # noqa: E402

Q0, T0 = 6.0e6, 700.0
LENGTH, WIDTH = 8e-3, 110e-3                 # m: the plate's depth along the heating and its extent across
R_FLAT, R_ZERO, R_READ = 32e-3, 45e-3, 3e-3  # m: full flux within R_FLAT, none beyond R_ZERO; recession read within R_READ
# R_FLAT - R_READ is half again the diffusion length 2 sqrt(alpha t), about 19 mm at 2 s: with full flux only to 12 mm
# the centre lost heat sideways and fell 9 % behind the exact solution by 2 s (measured while planning)
REFERENCE_1D = {1.0: 0.7890, 1.5: 2.7120, 2.0: 5.2190}       # mm, analysis/ablation_1d.py at 1 um and 0.5 ms (2026-10-12)


def fine_box():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    doc = json.load(open(os.path.join(here, "data", "reference_values", "box_q6p_fine.json")))
    return {t: float(np.interp(t, doc["t_s"], doc["r_total_mm"])) for t in (1.0, 1.5, 2.0)}


def taper(r, r_flat=R_FLAT, r_zero=R_ZERO):
    """1 within r_flat, a half cosine down to 0 at r_zero, 0 beyond."""
    x = np.clip((r - r_flat) / (r_zero - r_flat), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * x))


def projected_flux(body, radius, r_flat=R_FLAT, r_zero=R_ZERO):
    """Each column's flux [W/m2 over its top area]: Q0 times the taper at its corner's starting radius `radius`, times
    the projection of its top dual cell on the face -- a third of each touching top triangle times that triangle's
    -n_x, floored at zero -- over its area. A column on a tilted crater wall receives the heat its projection catches."""
    s = body.surface
    a, b, c = (s.points[s.faces[:, i]] for i in range(3))
    n2 = np.cross(b - a, c - a)                                   # twice the area times the outward normal
    proj = np.zeros(len(s.areas))
    for i in range(3):
        np.add.at(proj, s.faces[:, i], np.clip(-n2[:, 0], 0.0, None) / 6.0)
    return Q0 * taper(radius, r_flat, r_zero) * proj / s.areas


def solid_top(body):
    """Each column's top as if all its metal were solid: the base plus its mass over the solid density and its area,
    along its direction. The references hold their metal at one density, so this is the recession they measure; the
    skins' own top also swells as the mushy layer's density falls (by up to 0.15 mm before melting, measured)."""
    c = body.columns
    return c.points + (body.skins.column_mass() / (body.material.rho * c.area))[:, None] * c.normal


def run(dt=0.05, preset="dev", deform=4e-3, t_end=2.0, solver="skfem", verbose=True, h=2e-3, width=WIDTH, r_flat=R_FLAT,
        r_zero=R_ZERO, repeat_rule="column"):
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    m = mesh.box_mesh(LENGTH, width, width, h)
    m.points[:, 1:] -= 0.5 * width                                 # the face's centre on the x axis
    body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(),
                             skin.SkinSettings(cell=skin.SKIN_PRESETS[preset]), T0=T0, emissivity=0.0,
                             v_hat=(-1.0, 0.0, 0.0), deform_depth=deform, repeat_rule=repeat_rule)
    p = body._top0                                                 # the top surface at time zero: the plate's own
    front = np.abs(p[:, 0]) < 0.05e-3
    radius = np.where(front, np.hypot(p[:, 1], p[:, 2]), np.inf)
    centre = radius <= R_READ
    solid0, area0 = solid_top(body), body.surface.areas.copy()
    M0, rows, t, started, failure = body.mass(), [], 0.0, time.perf_counter(), None
    for k in range(int(round(t_end / dt))):
        try:
            q = projected_flux(body, radius, r_flat, r_zero)
            body.advance(t + dt, dt, heating.HeatingResult(q, Q0, 0.0, 0.0, 0.0))
        except RuntimeError as exc:                        # a move no rebuild can absorb: record where, and stop
            failure = {"t": t + dt, "error": str(exc), "fold_deg": body.moving.surface_fold()}
            print(json.dumps(failure), flush=True)
            break
        t += dt
        st = body.skin_stats()
        rows.append({"t": t, "recession_mm": float(((solid_top(body) - solid0) @ np.array([1.0, 0.0, 0.0]))[centre].mean() * 1e3),
                     "top_recession_mm": float(body._recession()[centre].mean() * 1e3),
                     "removed_recession_mm": float(body.column_removed[centre].sum() / (body.material.rho * area0[centre].sum()) * 1e3),
                     "centre_area_ratio": float(body.surface.areas[centre].sum() / area0[centre].sum()),
                     "T_top_centre": float(body.surface_temperature()[centre].mean()),
                     "mass_error": (body.mass() + body.removed_mass) / M0 - 1.0,
                     "energy_residual": body.energy_balance_residual(), "fold_deg": st["surface_fold_deg"],
                     "quality": st["surface_quality_min"], "zone_ratio": st["zone_ratio_min"], "rebuilds": st["rebuilds"],
                     "faults": st["jump_faults"], "repeats": st["interface_repeats"], "cells": st["skin_cells"],
                     "columns_over": st["repeat_columns_over"], "worst_ratio": st["mismatch_worst_ratio"],
                     "total_ratio": st["mismatch_total_ratio"],
                     "absorbed_kg": st["absorbed_mass_kg"], "handed_kg": st["handed_mass_kg"],
                     "wall_s": time.perf_counter() - started})
        if verbose and (k + 1) % max(1, int(round(0.25 / dt))) == 0:
            print(json.dumps({key: (round(v, 6) if isinstance(v, float) else v) for key, v in rows[-1].items()}), flush=True)
    return body, rows, int(centre.sum()), failure


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--dt", type=float, default=0.05)
    a.add_argument("--preset", default="dev", choices=tuple(skin.SKIN_PRESETS))
    a.add_argument("--deform-mm", type=float, default=4.0)
    a.add_argument("--t-end", type=float, default=2.0)
    a.add_argument("--solver", default="skfem")
    a.add_argument("--h-mm", type=float, default=2.0, help="3D element size [mm]")
    a.add_argument("--width-mm", type=float, default=WIDTH * 1e3)
    a.add_argument("--repeat-rule", default="column", choices=skinbody.REPEAT_RULES)
    a.add_argument("--r-flat-mm", type=float, default=R_FLAT * 1e3)
    a.add_argument("--r-zero-mm", type=float, default=R_ZERO * 1e3)
    a.add_argument("--out", default="skin_plate.json")
    args = a.parse_args()
    body, rows, n_centre, failure = run(args.dt, args.preset, args.deform_mm * 1e-3, args.t_end, args.solver,
                                        h=args.h_mm * 1e-3, width=args.width_mm * 1e-3, r_flat=args.r_flat_mm * 1e-3,
                                        r_zero=args.r_zero_mm * 1e-3, repeat_rule=args.repeat_rule)
    fine, compare = fine_box(), {}
    for t_ in (1.0, 1.5, 2.0):
        if not rows or t_ > rows[-1]["t"] + 1e-9:
            continue
        r = float(np.interp(t_, [x["t"] for x in rows], [x["recession_mm"] for x in rows]))
        g = float(np.interp(t_, [x["t"] for x in rows], [x["removed_recession_mm"] for x in rows]))
        compare[str(t_)] = {"centre_mm": r, "removed_mm": g, "fine_mm": fine[t_], "ref1d_mm": REFERENCE_1D[t_],
                            "rel_fine": r / fine[t_] - 1.0, "rel_1d": r / REFERENCE_1D[t_] - 1.0,
                            "removed_rel_1d": g / REFERENCE_1D[t_] - 1.0}
    summary = {"failure": failure, "completed_s": rows[-1]["t"] if rows else 0.0, "h_mm": args.h_mm, "dt": args.dt,
               "width_mm": args.width_mm, "r_flat_mm": args.r_flat_mm, "r_zero_mm": args.r_zero_mm,
               "repeat_rule": args.repeat_rule,
               "preset": args.preset, "deform_mm": args.deform_mm, "n_columns": body.columns.n, "n_centre": n_centre,
               "compare": compare, "worst_mass_error": max((abs(x["mass_error"]) for x in rows), default=0.0),
               "worst_energy_residual": max((abs(x["energy_residual"]) for x in rows), default=0.0),
               "fold_max_deg": max((x["fold_deg"] for x in rows), default=0.0),
               "rebuilds": rows[-1]["rebuilds"] if rows else 0, "faults": rows[-1]["faults"] if rows else 0,
               "repeats": rows[-1]["repeats"] if rows else 0, "wall_s": rows[-1]["wall_s"] if rows else 0.0,
               "wall_split": body.wall}
    print(json.dumps(summary, indent=1))
    json.dump({"summary": summary, "rows": rows}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
```

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

```python
"""Sub-plan 20 Task 10: the skins' cost per macro step on the default 100 mm sphere (2 mm surface, dense
band), production cells (10 um), 0.5 s macro steps, 10 ms sub-steps: three steps below melting (300 K, 1 MW/m2 at the
nose) and three melting (880 K, 2 MW/m2 at the nose), liquid removed as it forms."""
import json, sys, time
sys.path.insert(0, ".")
import numpy as np
from reentry_model import heating, material, mesh, skin, skinbody
from reentry_model.thermal import thermal_solver
out = {}
for label, T0, q0 in (("solid", 300.0, 1.0e6), ("melting", 880.0, 2.0e6)):
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    m = mesh.sphere_mesh(0.05, band=mesh.DEFAULT_BAND)
    t0 = time.perf_counter()
    b = skinbody.SkinBody(m, mat, thermal_solver("skfem"), mat.rho * m.volume(), skin.SkinSettings(), T0=T0, v_hat=(-1.0, 0.0, 0.0))
    setup = time.perf_counter() - t0
    steps = []
    for k in range(3):
        t1 = time.perf_counter()
        b.advance(0.5 * (k + 1), 0.5, heating.HeatingResult(q0 * np.clip(np.cos(b.theta), 0.0, None), q0, 0.0, 0.0, 0.0))
        steps.append(time.perf_counter() - t1)
    st = b.skin_stats()
    out[label] = {"columns": b.columns.n, "cells": st["skin_cells"], "n_nodes": b.mesh.n_nodes, "setup_s": setup, "step_s": steps,
                  "repeats": b.interface_repeats, "wall": b.wall, "removed_kg": b.removed_mass, "energy_residual": b.energy_balance_residual()}
    print(label, json.dumps(out[label]), flush=True)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "skin_timing.json", "w"), indent=1)
```

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
