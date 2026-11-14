# Level-set stage 3, amended — one skin depth, the gas surface evolved — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild the level-set skins on design §21: every skin 1.2 mm deep with no modes, the gas surface φ the
evolved level set and the skin base β derived 1.2 mm inside it, the columns cut from one nearest-point division of the
layer between them, runoff and sideways heat across that division's shared boundaries, and every column exactly full.
Then merge main and connect the skins to main's composite rheology and kinked-profile spraying.

**Architecture:**
- **Four standalone pieces, each with its own tests:**
  - the division of the layer (`division.py`, with the same-side ownership rule in `markers.py`);
  - deriving β and moving φ by exact volumes (`surfaces.py`, `levelset.py`);
  - skins with measured cross-sections, base exchange and combining by depth (`skin.py`);
  - sideways conduction and depth-resolved runoff faces (`skin.py`, `runoff.py`).
- **Two tasks that switch the body over** (`skinbody.py`, `stepend.py`):
  - the first makes the body move φ, keep every column full and retire the modes;
  - the second puts runoff, heat and re-partition on the division's faces and retires the tangent-plane faces.
- **Then:**
  - the verification cases;
  - the merge of main (gated on two other branches landing first);
  - the rheology law;
  - the kinked-profile spraying;
  - the cases re-judged on main's physics;
  - the documents.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, pyamg, scikit-fem; pytest with pytest-xdist. No new
dependency.

**Spec:** `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §21 (governs), with §20 for
what it does not change. The stage 3 plan this amends is `docs/superpowers/plans/levelset/03-stage3-skins.md`: its Tasks
1–9 and 9b stand as built, Task 10 (modes) is retired by Tasks 5–6 here, and its Tasks 11–12 are redone as Tasks 7, 11
and 12 here. Interface maps of the code at ed719a5 and of origin/main at 4e857e1 are in the plan's SDD workspace (Task
0): `levelset_map.md` and `main_map.md`. Every implementer brief points at the one it needs.

## Global Constraints

- **Where and how to run.** Paths are relative to `Reentry Model/` in the worktree `.claude/worktrees/levelset-stage3`
  (branch `levelset-stage3`). Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`, run from
  `Reentry Model/`. Paths contain spaces: quote them.
- **Shared physics only** (Asha, 2026-11-11). Every physical law, property and constant comes from
  `reentry_model.shared`; `reentry_model/levelset/` defines only numerical settings, each a named module constant or
  settings field with its source. `test_levelset_imports.py` enforces it.
- **Skins are never shared** (§4.2). Nothing in `reentry_model/levelset/` imports or copies the skyscrapers' skin
  modules. Tests may import them to compare.
- **Changes to `shared/` are behaviour-preserving** for the skyscrapers. Task 8's `same_results.py compare` against
  origin/main must say `identical`.
- **Never edit:**
  - the `subplan20` worktree;
  - `prototype/` (including `prototype/work-2026-10-28-skin-spraying-plan/` and `prototype/work-2026-10-28-skin-spraying/`);
  - anything under `code/reentry_model_output/`;
  - the DRAMA databases;
  - `sphere_reentry.py` and `sphere_sweep.py`;
  - the package at the repository root;
  - `prototype/proto3/` and `plan3`;
  - `.claude/worktrees/we-critical-308`;
  - `DRAMA Sphere Experiments/` and `Spheral/`;
  - `reentry_model/droplet_reentry/` and `reentry_model/wake_column/`;
  - the old top-level alias files in `reentry_model/`;
  - anything in `reentry_model/skyscrapers/`, except Task 8's conflict resolution in `skyscrapers/skinbody.py`.
- **Read-only (other sessions' files):** `docs/superpowers/plans/melt-spraying-subplans/23-notes.md`,
  `docs/superpowers/plans/2026-11-06-slip-band-spraying.md`, `docs/superpowers/specs/2026-11-08-refrozen-metal-in-skins-design.md`.
- **Conventions.** SI units. Level sets are negative inside metal and are signed distances after every reset. Normals
  point out of the metal. Enthalpy is above the shared `T_REF` (293 K) per kilogram.
- **Pinned by §21:**
  - the skin depth H = `SkinSettings.depth` = 1.2e-3 m, recorded with every run;
  - the opposed-normals threshold n_i · n_j < −0.5 (`markers.OPPOSED_NORMALS`, kept);
  - the cells as §20.6 (10 µm production, 20 µm development, over the 0.4 mm fine zone, each deeper cell 10 % thicker).
- **Books** (Task 11 of the stage 3 plan, unchanged): mass residual at most 1e-12 and enthalpy residual at most 1e-8
  per step, both relative (`StepBook.mass_residual`, `enthalpy_residual`).
- **Tests.**
  - Tests live in `tests/levelset/`, one file per module (`test_levelset_<module>.py`).
  - Whole-plate and box runs are marked `reference`.
  - Use `tests/levelset/materials.py` (`ConstMaterial`, `LinearMeltMaterial`) for simple properties and the shared
    AA7075 Scheil material where parity needs it.
  - Unit tier: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests/levelset -m "not reference" -q -p no:cacheprovider -n 4`.
- **Power and run length.** Run tests only on mains power (`pmset -g batt | head -1` says `AC Power`). A laptop run
  longer than about two hours goes to the Mac Studio through a runbook, and the task stops for Asha.
- **Commits.** Commit after each task, with named paths only (never `git add -A`, never a bare `git stash`). The
  message is one sentence on the physics or behaviour that changed, a blank line, then
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Do not push.** The merge of main waits for the other sessions' merges to land on main (Task 8's gate).

## Rulings made while planning (each with its cost if wrong)

1. **The division is built on level surfaces.** It uses `DIVISION_LEVELS = 24` uniform levels from 0 to H.
   - **Areas:** each level surface φ = −z_k (`cutcell.patches(mesh, phi + z_k)`) has its triangles mapped to their feet
     on φ (p + z_k ∇φ). The feet are divided among markers by `ownership`'s halfway-plane rule with the same-side
     rule, and each piece's area is taken on the level surface (the affine map keeps area fractions).
   - **Shared lengths:** a shared boundary's length at a level is the length of the dividing segments, mapped back to
     the level surface.
   - **Volumes:** each column's volume is the trapezoid integral of its areas.
   - **The per-piece factor:** the areas are scaled per piece by one factor, so that the columns' volumes sum exactly
     to `inside_volumes(phi) − inside_volumes(beta)` for that piece. That factor's departure from 1 is the closure
     measure §21.8 item 5 asks for. More than `DIVISION_TOL = 1e-3` fails the step's check (recorded, raised by tests,
     not by the body).
   - Cost if wrong: one function (`divide`) rewritten.
2. **φ moves by groups.** Each node within the band takes the marker owning its foot (node − φ∇φ, same-side rule), and
   the group's φ is offset by one distance d_i.
   - d_i starts at ΔV_i / A_i(0) and is corrected by a secant iteration on the group's swept volume (`inside_volumes`
     change over the group's cells), at most `MOVE_ROUNDS = 8` rounds, until the relative miss is at most
     `MOVE_RTOL = 1e-6`.
   - Then the distance is reset (`redistance`, which keeps the zero contour), and the body corrects per piece
     (`correct_volumes`).
   - A move longer than `MOVE_BAND_SHARE = 0.5` of the band is taken in rounds, as `_move_beta` takes moves longer
     than `stepend.MOVE_ROUND_SHARE` (0.5 of the fine margin) today. `levelset.py` defines its own constant, because it
     must not import `stepend`.
   - Cost if wrong: one function.
3. **Order of a part's acceptance** (Task 5):
   1. the cells' volume change ΔV_i;
   2. the φ move;
   3. the per-piece correction to (skins' cells volume + `solver.metal_volume()`);
   4. β = `skin_base`;
   5. `solver.move(beta)`, its enthalpy booked as today's owed heat;
   6. `divide` and `skins.set_profile`;
   7. the fill: columns with no base first, from their neighbours at their deepest shared levels; then every column
      with a base, at its base.

   Cost if wrong: the order inside `_accept`.
4. **The density difference of a molten cell returned to the volume** is carried as a volume change at the top of that
   column into the next part's ΔV_i. The last part's carry is applied by one more move at the step end, before the
   closure check. Cost if wrong: φ settles one part late.
5. **Sideways conductance.**
   - The conductance between cell k of column i and cell l of column j is the shared length integrated over the overlap
     of their depth spans, times the harmonic mean of their conductivities, divided by |x_i − x_j|: a two-point flux,
     like runoff's weight.
   - Runoff keeps the stage 3 arrival rule. The rate of cell k is (u_k · n_ij)⁺ L_ij(z_k) / A_i(z_k), with the length
     and area at the cell's mid-depth.
   - Cost if wrong: two coefficients.
6. **A marker released by the step end** (trapped in a fold, or left with no area) gives its skin to the nearest
   remaining marker by `combine(by="depth")`, not to the volume. The next fill settles both. Cost if wrong: one call
   in `_release`.
7. **Retired modules and code.**
   - `modes.py` and its tests are deleted.
   - `thin.py` and its tests are kept as stage 1's unused geometry: it is on main, and §21.14 leaves this to the plan.
   - The tangent-plane Voronoi faces, `_neighbours`, `patch_paths`, `NEIGHBOUR_REACH`, `NEIGHBOUR_PATH` and
     `Markers.cell_area` are deleted in Task 6.
   - Cost if wrong: nothing that §21 needs.
8. **Task 8's gate and baseline.**
   - The merge runs only when `interim-gas-loads` and `fitted-sphere-loads` are both ancestors of origin/main.
   - The baseline is `same_results.py run` from a detached temporary worktree of origin/main, in the scratchpad.
   - The three level-set parity tests that build skyscrapers skins pass `rheology="bridge"` until Task 9.
   - Cost if wrong: one rerun.
9. **The kinked-profile connection comes after the merge** (it needs main's `melt_profile`), and after the geometry, as
   §21.7 says. Cost: none.
10. **Two verification passes:** Task 7 on the new geometry with today's physics, and Task 11 after the merge with
    main's defaults. So a delayed merge does not hold up judging the geometry. Cost: one extra set of runs (about an
    hour).
11. **Sharpness on a P1 grid is judged against rounding.** An edge counts as sharp when the volume between φ's zero
    contour and the ideal sharp shape is at most a tenth of what rounding to the move's depth d would change,
    (1 − π/4) d² per unit length of edge. Cost: a test tolerance.

## Review Focus

- **A marker squeezed to almost no area by a move** (at an edge, after a split) must give its skin to its neighbour
  with the books exact, not divide by zero. Pinned by
  `test_a_marker_left_with_no_area_gives_its_skin_to_its_neighbour` (Task 5).
- **A part whose spraying or feed removal would take more than a column's depth** must be split (the part ratio),
  never leave a column with negative volume. Pinned by `test_a_part_that_would_empty_a_column_is_split` (Task 5).
- **A φ move longer than the distance-reset band** must be taken in rounds and still sweep its volumes. Pinned by
  `test_a_move_longer_than_the_band_is_taken_in_rounds` (Task 2).
- **A layer point exactly halfway between two markers** must belong to one of them, with the areas still summing to
  the surface's. Pinned by `test_a_point_halfway_between_two_markers_is_owned_by_one` (Task 1).
- **A body in two pieces** must close each piece's layer separately. Pinned by
  `test_two_pieces_close_their_layers_separately` (Task 1).

Carried from the stage 3 plan and kept passing: `test_rigid_motion_leaves_skins_unchanged`,
`test_a_skin_reduced_to_one_cell_keeps_its_books`, `test_a_trapped_marker_returns_its_skin_exactly` (rewritten in
Task 5 for ruling 6), and `test_marker_release_equals_the_shared_model`.

---

## Task 0 (controller only): the workspace

- [ ] Run this skill's `scripts/sdd-workspace` on this plan. Copy into the workspace the two scratchpad maps, as
  `levelset_map.md` and `main_map.md`:
  - `/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/eb7e9b74-487b-431a-863c-b8a3214cb080/scratchpad/plan21/levelset_map.md`
  - `/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/eb7e9b74-487b-431a-863c-b8a3214cb080/scratchpad/plan21/main_map.md`

  Also copy Task 11 fix round 2's nearest-point assignment sketch, `.../scratchpad/t11fix2/i3_assignment/`, as
  prior art for Tasks 1 and 5.
- [ ] Commit the design amendment and this plan:
  `git add "Reentry Model/docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md" "Reentry Model/docs/superpowers/plans/levelset/03b-stage3-one-skin-depth.md"`,
  then `git commit` with the message: "docs: design section 21 (one skin depth, the gas surface evolved) and its
  implementation plan", a blank line and the Co-Authored-By line.
- [ ] Record BASE and run the unit tier once. Expected: 200 passed.

### Task 1: The division of the layer

**Files:**
- Create `reentry_model/levelset/division.py`.
- Modify `reentry_model/levelset/markers.py` (`ownership`, line 205).
- Test `tests/levelset/test_levelset_division.py`, `tests/levelset/test_levelset_markers.py`.

**Interfaces:**
- Consumes:
  - `cutcell.patches(mesh, phi)`, `cutcell.inside_volumes(mesh, phi)`, `cutcell.cell_gradients`;
  - `levelset.pieces(mesh, phi)`, `levelset.redistance`;
  - `markers.ownership(x, tris, tri_patch, n_patch, size)`.
- Produces:
  - `markers.ownership(x, tris, tri_patch, n_patch, size, normals=None, tri_normal=None, boundaries=False)`.
    - With `normals` (n,3) and `tri_normal` (Q,3), marker m may own triangle q only if
      `normals[m] @ tri_normal[q] > OPPOSED_NORMALS`.
    - With `boundaries=True` it returns `(W, (i, j, length))`, the dividing segments' lengths summed per marker pair
      (i < j).
    - The defaults give exactly today's result.
  - `division.DIVISION_LEVELS = 24`, `division.DIVISION_TOL = 1e-3`.
  - `division.Division`, a frozen dataclass:
    - `levels` (B+1,) in m;
    - `area` (n, B+1) in m²;
    - `volume` (n,) in m³;
    - `i`, `j` (E,);
    - `length` (E, B+1) in m;
    - `normal` (E,3), unit, from i toward j;
    - `opposed` (E,) bool;
    - `piece` (n,) int and `factor` (P,);
    - property `base` = `area[:, -1]`;
    - method `cross_section(z)` (n, …), interpolating linearly in z and zero beyond H.
  - `division.divide(mesh, phi, depth, x, normal, levels=DIVISION_LEVELS, size=None) -> Division`, built as planning
    ruling 1 says.

- [ ] **Step 1: Write the failing tests** in `test_levelset_division.py`. Each test's assertions:
  - `test_a_flat_slab_divides_into_prisms`: a flat face with markers on a 2 mm grid.
    - Every interior column's area is 4.0e-6 m² at every level, to 1e-9 relative.
    - Every interior shared length is 2.0e-3 m at every level.
    - `abs(factor - 1) <= 1e-12`.
  - `test_a_spherical_shell_divides_by_the_shell_volume`: a 10 mm sphere, 2 mm spacing, H = 1.2 mm.
    - `abs(factor - 1) <= DIVISION_TOL`.
    - Every column's `area[:, -1] / area[:, 0]` is within 2 % of ((R − H)/R)².
  - `test_columns_beside_a_convex_edge_narrow_and_meet_on_the_halfway_plane`: the stage 3 box's 10 × 10 face edge.
    - A column owning a strip that touches the edge has `area` strictly decreasing with depth.
    - It shares a boundary of positive length at depth with a column on the other face.
    - That pair is not `opposed`.
  - `test_a_concave_corner_column_widens_with_depth`: an L-shaped body. The column on the concave edge has
    `area[:, -1] > area[:, 0]`.
  - `test_ownership_never_crosses_a_gap`: two parallel sheets 1.0 mm apart, with 2 mm marker spacing.
    - No triangle of one sheet is owned by a marker of the other.
    - No shared boundary joins them.
  - `test_metal_thinner_than_twice_the_depth_has_opposed_neighbours`: a slab 2.0 mm thick has `opposed.any()`; a slab
    3.0 mm thick has none.
  - `test_a_point_halfway_between_two_markers_is_owned_by_one`: each patch's owned areas sum to its area to round-off.
  - `test_two_pieces_close_their_layers_separately`: two separate cubes.
    - `len(factor) == 2`, each within `DIVISION_TOL` of 1.
    - Each marker's `piece` is its cube.

  In `test_levelset_markers.py`, add `test_ownership_with_no_normals_is_unchanged` (bitwise equal to the old call on
  the test sphere).
- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest tests/levelset/test_levelset_division.py tests/levelset/test_levelset_markers.py -q`.
  Expected: FAIL (`division` not found).
- [ ] **Step 3: Implement** `ownership`'s two options and `divide` as planning ruling 1 says.
- [ ] **Step 4: Run the tests and the unit tier.** Expected: all pass, and the unit tier still gives 200 plus the new
  tests.
- [ ] **Step 5: Commit** `division.py`, `markers.py` and the two test files. Message: "reentry_model: the layer
  between the gas surface and the skin base divided among the markers by nearest point, with cross-sections and shared
  boundaries at every depth".

### Task 2: The skin base derived from φ, and φ moved by exact volumes

**Files:**
- Modify `reentry_model/levelset/surfaces.py` and `reentry_model/levelset/levelset.py`.
- Test `tests/levelset/test_levelset_surfaces.py`, `tests/levelset/test_levelset_levelset.py`.

**Interfaces:**
- Consumes: `levelset.redistance`, `cutcell.inside_volumes`, and Task 1's `ownership` with `normals`.
- Produces:
  - `surfaces.skin_base(mesh, phi, depth, band=None) -> np.ndarray`, equal to `redistance(mesh, phi + depth, band)`.
    It raises `ValueError` if the band reaches less than `depth` into the metal.
  - `levelset.MOVE_ROUNDS = 8`, `levelset.MOVE_RTOL = 1e-6`, `levelset.MOVE_BAND_SHARE = 0.5`.
  - `levelset.Moved` (namedtuple `phi swept miss rounds`).
  - `levelset.move_by_volumes(mesh, phi, x, normal, dV, band=None, rounds=MOVE_ROUNDS, rtol=MOVE_RTOL) -> Moved`, as
    planning ruling 2 says. dV (n,) in m³ is positive outward, `swept` (n,) is achieved, and `miss` (n,) is relative.

- [ ] **Step 1: Write the failing tests.** The box is the stage 3 box at the cases' mesh: `box_body`'s `BoxBase`,
  `fine_size` 0.7 mm, `cell_floor` 0.25 mm. "Rounding volume" means (1 − π/4) d² × total edge length.
  - `test_the_skin_base_of_a_box_is_the_box_shrunk_by_the_depth_with_sharp_edges`: |V(β) − (8 − 2.4)(10 − 2.4)(10 − 2.4)
    mm³| is at most a tenth of the rounding volume with d = H.
  - `test_an_equal_recession_keeps_a_box_s_edges_sharp`: d = −1.0 mm on every group.
    - Every `miss` ≤ 1e-6.
    - |V − V_sharp| ≤ a tenth of the rounding volume.
  - `test_extra_loss_along_one_edge_rounds_it`: the groups within one spacing of one edge lose twice the others.
    - Every `miss` ≤ 1e-6.
    - φ at that edge's midpoint rises by more than at the face centres.
  - `test_a_deposit_on_a_convex_edge_rounds_it_by_the_depth_added`: d = +0.5 mm on every group.
    - Every `miss` ≤ 1e-6.
    - V is within a tenth of the rounding volume of (V_sharp_offset − rounding volume).
  - `test_a_receding_concave_corner_rounds_by_the_depth_removed`: an L-shaped body, d = −0.5 mm. The removed volume is
    within a tenth of the rounding term of (mitre removal − rounding volume along the concave edge).
  - `test_a_move_longer_than_the_band_is_taken_in_rounds`: `move_by_volumes` with a d of 1.5 times the band.
    - `rounds > 1`.
    - Every `miss` ≤ 1e-6.
- [ ] **Step 2: Run them to see them fail.** Expected: `skin_base` and `move_by_volumes` not found.
- [ ] **Step 3: Implement** `skin_base` and `move_by_volumes`.
- [ ] **Step 4: Run** the two test files. Expected: all pass.
- [ ] **Step 5: Commit** `surfaces.py`, `levelset.py` and the tests. Message: "reentry_model: the skin base derived as
  the surface one skin depth inside the gas surface, and the gas surface moved by each group's exact volume".

### Task 3: Skins with one depth, measured cross-sections, base exchange and combining by depth

**Files:**
- Modify `reentry_model/levelset/skin.py`.
- Test `tests/levelset/test_levelset_skin.py`.

**Interfaces:**
- Consumes: Task 1's `Division.levels` and `area` (as plain arrays; `skin.py` does not import `division`).
- Produces:
  - `SkinSettings.depth: float = 1.2e-3`. `__post_init__` raises unless 2 × cell ≤ min(fine_zone, depth) and
    fine_zone ≤ depth. `stub` stays until Task 5 deletes it.
  - `Skins.set_profile(levels, area)`: once called, `cross_section(z)` interpolates the profile (zero below its last
    nonzero level and beyond H), `area` = `area[:, 0]`, and `dz` inverts the piecewise-quadratic cumulative volume.
    `kappa` is then ignored (the curvature path is deleted in Task 5).
  - `Skins.column_volume() -> (n,)` and `Skins.cells_volume(T=None) -> (n,)`.
  - `Skins.exchange_base(rows, dV, T_draw=None) -> BaseExchange(mass, H, P, cells)`.
    - dV > 0 adds volume metal of that volume at `T_draw`, at the solid's density, below the bottom cell.
    - dV < 0 removes the bottom −dV of volume, splitting at most one cell, whatever the phase, and returns its mass,
      enthalpy and momentum.
  - `Skins.take_between(rows, z_top, z_bottom) -> (m, H, P)` removes the metal between two depths, splitting cells.
  - `Skins.combine(keep, gone, by="mass")`: `by="depth"` adds each gone column's metal to the keeper's cells at the
    same depths (a one-dimensional conservative remap on depth intervals). `by="mass"` is today's rule.

- [ ] **Step 1: Write the failing tests.**
  - `test_a_measured_profile_holds_a_spherical_shell`: levels with A₀((R − z)/R)². The column volume is within 1e-6 of
    the analytic shell's.
  - `test_a_column_that_pinches_fills_exactly_to_its_pinch`: the area falls linearly to 0 at 0.8 mm.
    - The filled column's cells end at 0.8 mm, to 1e-12 m.
    - The volume is exact to 1e-12 relative.
  - `test_conduction_books_close_on_a_measured_profile`.
  - `test_the_base_exchange_is_exact`: draw and return, with mass, enthalpy and momentum to 1e-12; a molten bottom cell
    is returned.
  - `test_combining_by_depth_keeps_each_depth_s_metal`: mass, enthalpy and momentum exact; the hot top stays above the
    cold bottom.
  - `test_taking_metal_between_two_depths_is_exact`.
  - Update `test_settings_defaults_and_presets` (`depth` 1.2e-3).
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `tests/levelset/test_levelset_skin.py` and the unit tier. Expected: all pass. The old curvature
  path is untouched.
- [ ] **Step 5: Commit.** Message: "reentry_model: level-set skins one depth deep, with cross-sections taken from a
  measured profile, an exact exchange at the base and combining by depth".

### Task 4: Sideways conduction and depth-resolved runoff faces

**Files:**
- Modify `reentry_model/levelset/skin.py` (`conduct`, `HeatBook.E_side`) and `reentry_model/levelset/runoff.py`.
- Test `tests/levelset/test_levelset_skin.py`, `tests/levelset/test_levelset_runoff.py`.

**Interfaces:**
- Consumes: Task 1's `Division`; Task 3's profile.
- Produces:
  - `skin.SideFaces` (namedtuple `i j length levels distance`).
  - `Skins.conduct(dt, q_top, T_amb, emissivity, T_b, trial=False, side=None)`: with `side`, the sideways
    conductances of planning ruling 5 enter the same implicit solve. The solver is the builder's (§21.14).
    `HeatBook.E_side` is per column and sums to zero.
  - `runoff.Sheet.from_division(division, markers, skins=None)`. It gives `f = length[:, 0, None] * normal`,
    `w = length[:, 0] / |x_j − x_i|`, plus `length_z` (E, B+1) and `levels`. `transport` uses the cell-depth length
    and area of ruling 5; with constant lengths and areas it equals today's transport exactly.

- [ ] **Step 1: Write the failing tests.**
  - `test_a_lateral_temperature_step_relaxes_at_the_analytic_rate`:
    - set-up: 40 columns in a row at 0.25 mm spacing, `ConstMaterial`, insulated top and bottom, a step of 100 K, run
      until √(αt) = 2 mm;
    - assertion: every cell is within 1 % of the step of the error-function solution.
  - `test_sideways_conduction_closes_its_books`: Σ E_side is zero to 1e-12 of Σ |E_side|.
  - `test_depth_faces_with_constant_lengths_transport_as_before`: on a flat sheet, `from_division` and `of` give the
    same transport, to 1e-14.
  - `test_a_uniform_flow_over_a_flat_face_neither_gathers_nor_spreads`: every interior column's mass change is at most
    1e-12 relative.
  - `test_runoff_crosses_a_box_edge`: metal flowing toward the edge on one face arrives in the other face's columns at
    matching depths, and the books close.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the two test files and the unit tier. Expected: all pass.
- [ ] **Step 5: Commit.** Message: "reentry_model: level-set skins conduct heat sideways, and runoff crosses
  depth-resolved faces between columns".

### Task 5: The body moves φ and keeps every column full

**Files:**
- Modify:
  - `reentry_model/levelset/skinbody.py` and `reentry_model/levelset/stepend.py`;
  - `surfaces.py` (`outer_surface` deleted);
  - `runoff.py` (`Sheet.top`, `mixed`, `_zone_cells` and `thickness_step`'s `level` removed);
  - `spraying.py` (`b_mode`, `depth_phi_psi` and `top` removed; b is always the skin's);
  - `markers.py` (`top` and `hold` removed);
  - `skin.py` (`stub`, the curvature cross-section, `minimum_mass`, `draw`, `hand_down`, `take_bottom` and
    `append_bottom` deleted);
  - `analysis/levelset/stage3_cases.py`, only as far as it must change to import and run: the builders pass `phi0` =
    the shape itself, without `+ stub`; the use of `modes.melt_depth` (line 671), `ModeWatch` and the mode fields in
    `step_row` go. Task 7 does the rest.
  - `tests/levelset/test_levelset_stage3_cases.py`: `test_a_failing_step_keeps_the_body_it_started_from` monkeypatches
    `_rebuild` instead of the deleted `_modes`.
- Delete `reentry_model/levelset/modes.py` and `tests/levelset/test_levelset_modes.py`.
- Test `tests/levelset/test_levelset_skinbody.py`, and the tests the map lists under "two modes", "curvature
  cross-section / cap" and "beta being the evolved surface" (`levelset_map.md` §4), each rewritten to the new
  behaviour or deleted with its reason in the commit message.

**Interfaces:**
- Consumes: Tasks 1–4.
- Produces:
  - `LevelSetSkinBody.__init__(self, material, base, phi0, T0, settings=None, loads=None, spray_model=None,
    leeward="shadow", emissivity=0.0, fine_margin=4e-3, fine_size=0.7e-3, cell_floor=0.25e-3, remove="spray",
    spacing=MARKER_SPACING, v_hat=(1.0, 0.0, 0.0), T_ambient=0.0, slurry_spray=False, linear_solver="amg")`.
    - `phi0` is a callable giving the initial gas surface (negative inside).
    - `mode_settings` is gone, and `__setstate__` drops it from old pickles.
  - Body attributes:
    - `phi`, evolved;
    - `beta`, always `skin_base(mesh, phi, settings.depth, band)`;
    - `division`;
    - `psi`.
  - `skinbody.FILL_RTOL = 1e-12`, `FILL_ROUNDS = 8`, `BASELESS_SHARE = 1e-6` (a column with base area at most this
    share of its top area has no base).
  - `stepend.THIN_METAL` (message prefix): "metal thinner than twice the skin depth".
  - New `StepBook` fields:
    - `layer_factor_max`, the largest |factor − 1|;
    - `fill_miss_max`;
    - `move_miss_max`;
    - `n_baseless`.
  - Deleted `StepBook` fields: `n_top`, `modes_to_top`, `modes_to_full`, `mass_modes_handed`, `mass_modes_absorbed`,
    `enthalpy_modes_J`, `residual_modes_J`, `momentum_modes_normal`, `mass_absorbed`, `enthalpy_absorbed`.
  - The part ratio becomes max |ΔV_i / A_i(0)| / (`MOVE_SHARE` × spacing).

- [ ] **Step 1: Write the failing tests** in `test_levelset_skinbody.py`. The fixtures build bodies from `phi0`
  (the cube itself). The `owed` helper is rewritten without the mode fields.
  - `test_a_cold_box_keeps_its_edges_and_columns`: 2 mm spacing, no heat, one step.
    - Nodal φ within the band changes by at most 1e-12 m.
    - Every column's volume is unchanged to 1e-12 relative.
    - `mass_drawn == mass_handed == 0`.
  - `test_every_step_closes_the_layer_and_the_books`: the heated cube, three steps.
    - `layer_factor_max ≤ DIVISION_TOL`, `fill_miss_max ≤ 1e-12`, and every piece's V(φ) equals its metal's volume to
      1e-12.
    - The books are within the Global Constraints.
  - `test_a_sprayed_column_draws_from_the_volume_at_its_base`: a face column's drawn volume equals its sprayed volume to
    1e-12, and β falls beneath it.
  - `test_a_column_with_no_base_takes_its_shortfall_from_its_neighbours`: a sprayed edge column with `base` 0 never
    exchanges at a base, and the books are exact.
  - `test_a_molten_cell_returned_to_the_volume_lowers_phi_above_its_column`: φ falls by (1 − ρ_l/ρ_s) times the cell's
    volume over the column's area, within 5 %, by the step end.
  - `test_a_remesh_carries_phi_and_rebuilds_beta`: the piece volumes are preserved to 1e-12, and β equals
    `skin_base(phi)`.
  - `test_metal_thinner_than_twice_the_skin_depth_stops_the_run`: a 2.0 mm slab raises `RuntimeError` starting with
    `THIN_METAL` at construction.
  - `test_a_marker_left_with_no_area_gives_its_skin_to_its_neighbour`: mass, enthalpy and momentum exact.
  - `test_a_part_that_would_empty_a_column_is_split`: a feed removal of more than H in one part gives `parts > 1`, and
    no column's volume goes negative.
  - Rewrite `test_a_trapped_marker_returns_its_skin_exactly` (ruling 6) and
    `test_rigid_motion_leaves_skins_unchanged`.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement.**
  - The body is built from `phi0`, with markers on φ and the skins filled at T0 to exactly their columns (this
    replaces `_fill_stub`).
  - Within a part, the passes run as today, without the in-pass draw (`m_min` and the draw cap go). The base groups
    come from the division's base ownership (the owner of each β patch's foot on φ), and columns with no base have an
    insulated bottom.
  - Acceptance follows planning ruling 3, with ruling 4 for the density carry.
  - At the step end: no φ rebuild; ψ from φ; the remesh carries φ (interpolated, distance reset, corrected per piece)
    and re-derives β; then settle (snap, trapped, areas, split, merge, Task 9b's re-partition with
    `combine(by="depth")`); then divide, `set_profile` and fill.
  - The thin-metal stop is raised if `division.opposed.any()`, at construction and at every step end.
- [ ] **Step 4: Run** the unit tier. Expected: all pass.
- [ ] **Step 5: Commit** the named files. Message: "reentry_model: the level-set body moves the gas surface by each
  column's volume and derives the skin base one skin depth inside it, every column exactly full and no skin modes".

### Task 6: The body on the division's faces

**Files:**
- Modify:
  - `reentry_model/levelset/skinbody.py` (`Sheet.from_division`, and `side` passed to `conduct`);
  - `runoff.py` (`Sheet.of` deleted);
  - `spraying.py` (`region_extent` edges from the division's level-0 faces);
  - `repartition.py` (`_fit_columns` across the division's faces);
  - `markers.py` (as planning ruling 7 says; `tessellate` keeps the areas by ownership, the normals and the curvatures);
  - `stepend.py` (`_tessellate` keeps only the bare-marker release);
  - `tests/levelset/cases.py` (`nearest_ratio` and `seeded` on the division's areas).
- Test `tests/levelset/test_levelset_skinbody.py`, `test_levelset_runoff.py`, `test_levelset_markers.py`,
  `test_levelset_repartition.py`, `test_levelset_spraying.py`, `test_levelset_loads.py`.

**Interfaces:**
- Consumes: Tasks 1, 4 and 5.
- Produces: no new names. `Markers.faces`, `Markers.cell_area` and `Sheet.of` are gone, and every face comes from
  `body.division`.

- [ ] **Step 1: Write the failing tests.**
  - `test_heat_passes_a_box_edge_between_the_faces_columns` (marked `reference`):
    - set-up: the box at 700 K, 1 MW/m² on one face for 0.5 s, nothing melting;
    - reference: `CutThermal` alone on the same box with `cell_floor` 0.1 mm near the edge;
    - assertion: the temperature rise at 0.6 mm depth, 0.5 mm from the edge on the unheated face, is within 5 %.
  - `test_runoff_crosses_a_box_edge_through_the_body`: books exact, and metal arrives in the other face's columns.
  - `test_no_metal_runs_across_a_slot`: two sheets 1.0 mm apart exchange nothing.
  - `test_a_repartition_by_depth_moves_no_metal_through_the_base_and_keeps_each_depth`: mass, enthalpy and momentum
    exact, and each level's metal stays at its level, to 1e-12.

  Delete the tessellation tests the map lists (`levelset_map.md` §4), each with its reason.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the unit tier and `test_heat_passes_a_box_edge_between_the_faces_columns`. Expected: all pass.
- [ ] **Step 5: Commit.** Message: "reentry_model: level-set runoff, sideways heat and re-partition cross the shared
  boundaries of one division of the layer, and the tangent-plane faces are gone".

### Task 7: The verification cases on the new geometry

**Files:**
- Modify `analysis/levelset/stage3_cases.py`, `analysis/levelset/stage3_recorded.json` and
  `tests/levelset/test_levelset_stage3_cases.py`.

**Interfaces:**
- Consumes: Task 5's `phi0` and `StepBook` fields.
- Produces: the recorded numbers that Task 12 cites.

- [ ] **Step 1: Bring the script up to date.**
  - `plate_top`, `box_face` and `box_centre` read φ.
  - `step_row` records `layer_factor_max`, `fill_miss_max`, `move_miss_max`, `n_baseless` and a `division` cost line.
  - The box also records each edge's rounding radius, estimated from φ's zero contour at the edge midpoints.
- [ ] **Step 2: Estimate the run time.** The box's first step at 2 mm took 31.4 s on the old design. If the plates,
  box and shoulder together would take over two hours, write a Mac Studio runbook (as
  `docs/superpowers/plans/levelset/2026-11-12-stage3-cost-macstudio.md`) and stop for Asha.
- [ ] **Step 3: Run each case:**
  - the plates, at development and production cells;
  - the box, at development cells for 2 s;
  - the shoulder.
- [ ] **Step 4: Judge each case:**
  - **cold plate:** within 5 % of `REFERENCE_1D` and `REFERENCE_1D_LIQUIDUS` at 1.0, 1.5 and 2.0 s;
  - **hot plate:** recession within 5 % at 1.5 and 2.0 s, and non-rigid depth within 10 %. This plate's non-rigid depth
    passes 1.2 mm, so its melting front crosses the skin base: the risk of §21.12;
  - **box:** the face-centre recession within 5 % of the cold reference at 1.0, 1.5 and 2.0 s; edges reported;
    the cost per step against 31.4, 175.3, 772.1 and 1175.9 s;
  - **shoulder:** the deposit freezes;
  - **all cases:** the books.

  A case that passes loses its strict xfail. A case that fails keeps `xfail(strict=True)`, with the measured reason,
  and is reported to the controller (not hidden).
- [ ] **Step 5: Commit** the script, the recorded numbers and the test. Message: "analysis: stage 3's plates, box and
  shoulder re-judged with one skin depth and the gas surface evolved".

### Task 8: Merge main (controller, gated)

**Files:**
- The merge.
- Conflict resolution in `README.md`, `Reentry Model/CLAUDE.md`, `tests/test_reentry_model_shared.py`,
  `reentry_model/shared/spray.py` and `reentry_model/skyscrapers/skinbody.py` (`main_map.md` §1).
- The three parity tests in `tests/levelset/test_levelset_spraying.py` (lines 137, 170 and 292 at ed719a5).

**Interfaces:**
- Produces:
  - `spray.slurry_modes(..., law="bridge", rate_cap=None)`. `rate_cap=None` means `rheology.rate_cap(law)`.
  - `spray.plan_source_block(..., extra=None, extra_slurry=None)`, passing the rows through to `source_block`.

- [ ] **Step 1: Check the gate.**
  - Run `git fetch origin`.
  - Then `git merge-base --is-ancestor interim-gas-loads origin/main && git merge-base --is-ancestor fitted-sphere-loads origin/main`.
  - If either fails, stop. Report to Asha which branch has not landed, and do not merge.
- [ ] **Step 2: Make the baseline.**
  - `git worktree add --detach "<scratchpad>/main-base" origin/main`.
  - `"$PY" analysis/skyscrapers_port/same_results.py run --code "<scratchpad>/main-base/Reentry Model" --out "<scratchpad>/same/main" --cases skins,element,nomelt_elements`.
- [ ] **Step 3: Merge and resolve.** Run `git merge --no-ff --no-commit origin/main`, then resolve:
  - **The lists:** `README.md`, `CLAUDE.md` and `SHARED_MODULES` each list both `nose` and `melt_profile`.
  - **`spray.py`:** keep both sides.
  - **`skinbody.py`, hunk 1:** the skyscrapers pass main's profile rows through `plan_source_block`'s `extra`.
  - **`skinbody.py`, hunk 2:** keep main's `layer_columns` and the branch's aliases.
  - **`skinbody.py`, hunks 3–4:** the skyscrapers call the shared `slurry_modes` with `law=sk.settings.rheology` and
    their viscosity through `skin_mod.viscosity_of`, so they follow `--rheology` exactly as main does.
  - **The parity tests:** they pass `rheology="bridge"` to the skyscrapers skins they build.
- [ ] **Step 4: Verify.**
  - From `Reentry Model/` in the merged worktree, run `same_results.py run --code . --out "<scratchpad>/same/merged"`
    with the same cases, then `compare` it against the baseline. Expected: `identical`.
  - Then run, all of which must pass:
    - `tests/test_reentry_model_shared.py`;
    - the skyscrapers' `tests/test_reentry_model_skin.py`, `test_reentry_model_skinbody.py` and
      `test_reentry_model_spray.py`;
    - the level-set unit tier.
- [ ] **Step 5: Commit the merge and clean up.** The message names the three main commits and the conflict rule for
  `slurry_modes`. Then run `git worktree remove "<scratchpad>/main-base"`. Do not push.

### Task 9: The rheology law

**Files:**
- Modify `reentry_model/levelset/skin.py` and `reentry_model/levelset/spraying.py`.
- Test `tests/levelset/test_levelset_skin.py` (lines 425–429 and 485–491 at ed719a5) and
  `tests/levelset/test_levelset_spraying.py`.

**Interfaces:**
- Consumes: main's `rheology.LAWS`, `DEFAULT_LAW`, `law_viscosity`, `rate_cap` and `bridge_solid_fraction`.
- Produces:
  - `SkinSettings.rheology: str = rheology.DEFAULT_LAW`, checked against `LAWS`. `__setstate__` gives old pickles
    `"bridge"`.
  - `Skins._viscosity`, through `law_viscosity(self.settings.rheology, …)`. Every rate cap in the level set uses
    `rate_cap(law)`.
  - The `slurry_modes` call, with the skins' law.

- [ ] **Step 1: Write the failing tests.**
  - `test_composite_viscosity_ignores_the_shear_rate`: one column's cells at three shear rates give the same
    viscosities, equal to `law_viscosity("composite", …)` to 1e-15.
  - `test_the_same_column_gives_the_same_viscosities_in_both_architectures`.
  - `test_bridge_viscosity_is_unchanged`: equal to today's to 1e-15.
  - The three parity tests drop `rheology="bridge"`.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the unit tier. Expected: all pass.
- [ ] **Step 5: Commit.** Message: "reentry_model: level-set skins follow main's rheology setting, composite by
  default".

### Task 10: The kinked-profile spraying

**Files:**
- Modify `reentry_model/levelset/spraying.py`, `skinbody.py` (summary and `record_names`) and
  `analysis/levelset/stage3_cases.py` (`SprayModel(..., dispersion="profile")` where the cases spray).
- Test `tests/levelset/test_levelset_spraying.py`.

**Interfaces:**
- Consumes: main's `melt_profile.LayerColumns`, `profile_modes`, `select_run`, `run_properties`, `profile_rows` and
  `summary`, and `spray.SprayModel.evaluate(..., cap_depth=, profile=, rt_extent=)` (`main_map.md` §2a).
- Produces:
  - `spraying.layer_columns(skins, tau) -> melt_profile.LayerColumns`, with `front` = `depth_above(T_rigid)`, and the
    column's own depth where no rigid metal lies within it.
  - The profile branch of `spray_plan`, mirroring the skyscrapers' template (`main_map.md` §2c). It is taken when
    `model.dispersion == "profile"`, and the table path is left byte for byte.
  - `cascade_mass(..., T_run=None)` and `release_top(..., T_run=None)`.
  - `SprayPlan` gains `T_run`, `resolved`, `profile`, `modes` and `record`.
  - The records are 36 wide on the profile path. `LevelSetSkinBody.record_names()` exists, and the step summary gains
    `DISPERSION_COLUMNS`.

- [ ] **Step 1: Write the failing tests.**
  - `test_the_table_path_is_unchanged`: records and released masses are byte-identical to before on the existing
    spraying fixture.
  - `test_the_profile_release_matches_the_skyscrapers_on_one_column`: the same cells in both architectures give mass,
    enthalpy and momentum within 1e-12, and equal records.
  - `test_a_skin_flowing_to_its_base_takes_its_base_as_the_wall`: with no kink, `wall_depth` equals the column depth.
  - `test_profile_records_are_36_wide`.
  - Girin's published cases are unchanged.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the unit tier. Expected: all pass.
- [ ] **Step 5: Commit.** Message: "reentry_model: level-set skins spray by the kinked profile of their own cell
  speeds, the table kept".

### Task 11: The cases re-judged on main's physics

**Files:** `analysis/levelset/stage3_recorded.json` and `tests/levelset/test_levelset_stage3_cases.py`.

- [ ] **Step 1: Rerun Task 7's cases** with the defaults composite and profile, under the same two-hour rule.
- [ ] **Step 2: Judge them by Task 7's criteria.** Record each number beside its Task 7 value, with the change and its
  cause.
- [ ] **Step 3: Commit.** Message: "analysis: stage 3's cases re-judged with the composite rheology and the
  kinked-profile spraying".

### Task 12: Documentation

**Files:**
- Modify:
  - the design (§21.8, each item marked built with its test and measured result; §21.12, the measured costs);
  - `docs/model_assumptions.md` §12 (stage 3 as amended: each rule with its test, or marked assumed);
  - the repository `README.md` (the level-set section);
  - `docs/superpowers/plans/levelset/parity-ledger.md`;
  - `Reentry Model/CLAUDE.md` (the level set through stage 3, one skin depth);
  - `docs/superpowers/plans/levelset/03-stage3-skins.md` (a header note: Tasks 10–12 superseded by this plan).

- [ ] **Step 1: Write the documents.** Every number cites the test or run that measured it, and a measurement made
  once is said to be once.
- [ ] **Step 2: Run** `tests/levelset`, `tests/test_reentry_model_shared.py` and `tests/test_reentry_model_cli.py`.
  Expected: all pass.
- [ ] **Step 3: Commit.** Message: "docs: stage 3 of the level set as amended by design section 21, with its measured
  results".
