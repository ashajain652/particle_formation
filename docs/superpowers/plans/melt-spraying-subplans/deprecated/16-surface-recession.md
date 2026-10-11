# Sub-plan: Task 16 — The derived surface (`recession.py`)

> **New task, added 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md` §8.
> Not present in `2026-09-20-melt-spraying.md`; there is no extracted body to read. Read `00-shared-context.md`
> first, including the 2026-09-27 amendment block (facts 38–45), then Task 1's amendment.

**Depends on:** Task 1's amendment (`derived_points`, `surface(derived=True)`, `smoothed_normals`, `element_depths`)
and Task 9's `φ_e`. Nothing else.
**Produces, for later tasks:** the displaced surface every consumer of the shape reads (Task 9's nose fit, drag and
projected area; Task 5's θ; Task 6's tangents), the surface Task 17 hands to gmsh, and the
`recession_residual_frac` diagnostic.
**Character:** numerics/geometry, pure arrays, no physics of its own.
**Read before implementing:** facts 42 and 43. Fact 43 is the one that shapes the design — Euler forces a closed
triangulation to have **twice as many triangles as vertices** (measured 18 078 against 9 041), so the per-patch
volume constraints outnumber the nodal unknowns 2:1 and **cannot** all be satisfied. Do not try; solve least
squares, force the total exactly, and report the residual.
**Refinement goal for the sub-agent:** turn this into a standalone implementation plan — function signatures, test
plan, acceptance criteria — precise enough to implement without reading the spec.

---

## What this task does

`φ_e` stays the authoritative mass account against the **original** element volumes (spec §5, Design B). This task
computes, each macro step, where the surface *would* be if that mass account were rendered as geometry, writes it
into `mesh.derived_points`, and leaves `mesh.points` bit-identical. Nothing in `thermal/` changes, which is the
whole point: the solver never sees a moved node, so skfem keeps its 0.015 s assembly rather than 0.262 s (fact 45).

The problem being solved is that the shape is currently a **binary** function of which elements are alive, so the
surface recedes in jumps of about 6 mm (fact 39) with staircase steps one cell tall and one cell wide — a 45°
deviation in the patch normal, which is what the Lees distribution and the film tangents read.

## The algorithm

Each macro step, after element death. The displacement is a **state function of `φ_e` and the active set**,
recomputed from scratch, never accumulated, so it cannot drift away from the mass account.

1. **Consumed volume per patch.** For each live surface patch `p`, `ΔV_p` is the volume of the dead elements in its
   column plus the fill deficit `(1 − φ_e) V_e` of the live ones it shadows. This reads `φ_e` directly and is the
   **only** place the consumed depth is derived — there is no second sub-cell correction, and therefore no way to
   double count. Because a half-consumed live element contributes its deficit here, the surface is drawn half way
   through that element, which removes the recession quantum entirely rather than merely reducing it.
2. **Nodal displacements.** Solve in least squares `M s = ΔV`, where `s_i` is node `i`'s displacement from its
   original position along its area-weighted vertex normal and `M[p, i] = (A_p / 3)(n_i · N_p)` for each of patch
   `p`'s three nodes. Use `scipy.sparse.linalg.lsqr`. The system is overdetermined 2:1 and generically inconsistent.
3. **Force the total exactly.** Apply one uniform offset along the vertex normals, found by bisection, so the
   enclosed volume of the derived surface equals `mass / ρ` exactly. Measured convergence: 1.1·10⁻¹⁴ (fact 40).
4. **Sub-step if the motion is large.** If any `s_i` increment would exceed `RECESSION_SUBSTEP_FRAC = 0.35` of the
   local cell size, subdivide the step and record that it happened. Fact 42: safe to about 40 % of a cell, collapse
   past about 70 %. The peak recession rate is **unmeasured** — fact 5's "3.6 layers melt per step" is about the
   rejected surface-only feed rule and about crossing the **liquidus**, not about what leaves, and fact 28's depth
   gate holds melt in place — so this guard is what makes the design safe under either answer.
5. **Write** `mesh.derived_points`, and report `recession_residual_frac`, the median absolute per-patch mismatch.
   **Write boundary nodes only.** This is an invariant, not an incidental: measured 2026-09-27, every one of the
   18 506 boundary faces on the layered fixture has an *interior* node as its tetrahedron's opposite vertex (0 of
   18 506 are themselves boundary nodes), so while only boundary nodes move, `surface()`'s outward-orientation test
   reads the same coordinate whichever array it indexes. Displace an interior node and that stops being true, and
   the orientation of the derived surface begins to depend on which coordinate system that test reads — a
   distinction **no mesh test can catch**, because the two code paths are numerically identical for every
   boundary-only displacement, valid or inverted (measured: 0/4598 patches inward either way at 0.1 mm; 4598/4598
   inward either way at 1.0 mm, where 13 794 elements are inverted).

## Interfaces

- `RECESSION_SUBSTEP_FRAC = 0.35`, `LSQR_ATOL = 1e-12`, `LSQR_ITER_LIM = 4000`
- `patch_consumed_volume(mesh, phi) -> (n_patches,)` — step 1
- `nodal_normal_operator(surface) -> csr_matrix` — the `M` of step 2, built once per mesh and cached
- `solve_recession(M, dV) -> (s, residual_frac)` — steps 2 and 5
- `correct_total(points, faces, s, target_volume) -> s` — step 3
- `update_derived_surface(mesh, phi, rho, mass) -> RecessionResult(residual_frac, substeps, max_displacement_frac)`

## Tests

1. **Uniform recession is exact.** A uniform `ΔV_p = A_p d` on the sphere gives every `s_i = d` to 10⁻¹⁰, and the
   derived volume matches `V₀ − A d` — the one case where the overdetermined system is consistent.
2. **A smooth field fits well.** A √cos θ field gives median per-patch error below 0.5 % (fact 43 measured 0.03 %)
   and total error below 10⁻¹² after correction.
3. **A patchy field still conserves.** A field with 40 % of patches zeroed gives a total error below 10⁻¹² and a
   **non-zero** reported residual. The test asserts the residual is *measured*, not that it is small — fact 43 puts
   it at 22 % for this case and that is the honest answer.
4. **`points` is never modified.** Bit-identical before and after; `mesh.volume()` unchanged.
5. **Sub-stepping fires.** A step demanding 0.8 of a cell produces more than one sub-step and no inverted element in
   the derived surface.
6. **Half a consumed element shows.** One surface element at `φ_e = 0.5` moves its patch half an element inward, not
   zero and not a whole element.

## Acceptance criteria

1. Uniform recession exact to 10⁻¹⁰; total volume exact to 10⁻¹² in every test.
2. `mesh.points` bit-identical after any number of steps.
3. No inverted element in the derived surface under any test, including the sub-stepping one.
4. The solve costs under 0.1 s per step at 15 430 patches — **measure it**; if it does not, the fallback is
   area-weighted averaging of incident patch recessions, one sparse matrix-vector product, whose accuracy has not
   been measured (spec §12 question 3).
5. `"$PY" -m pytest tests/test_reentry_model_recession.py -q` green.

## Not in this task

Remeshing (Task 17). Any change to `thermal/` — Design B means the solver's mesh never moves.
