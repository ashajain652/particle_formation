# Sub-plan: Task 17 — Periodic remeshing (`remesh.py`)

> **New task, added 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md` §7.
> Not present in `2026-09-20-melt-spraying.md`. Read `00-shared-context.md` first, including the 2026-09-27
> amendment block (facts 38–45), then Tasks 1 and 16.

**Depends on:** Task 1's amendment (the band size field, `element_depths`, `edge_lengths`) and Task 16 (the derived
surface this consumes). Task 10 calls it.
**Produces, for later tasks:** a renewed mesh around the receded shape, and the conservation residuals Task 14
reports and Task 15 documents.
**Character:** geometry plus gmsh plumbing plus conservative field transfer. The riskiest task in the amendment.
**Read before implementing:** facts 40 and 41 — they are measurements of this exact pipeline, end to end, on the
real mesh. Do not re-derive it. In particular: the death boundary is **already** closed and manifold, so nothing is
reconstructed; Taubin smoothing is topology-agnostic because it is purely local; and `classifySurfaces` costs 0.1 s
on a *smoothed* surface and would cost far more on a raw staircase, which is why smoothing comes first.
**Refinement goal for the sub-agent:** turn this into a standalone implementation plan with signatures, a test plan
and acceptance criteria.

---

## Triggers

Evaluated once per macro step by Task 10; first to fire wins, subject to `REMESH_MIN_STEPS = 20` so they cannot
thrash.

1. **Band consumed** — the maximum `element_depths()` over elements that currently own a surface patch exceeds
   `band_thickness / 2`. With a 15 mm band and fact 14's 77 mm nose recession this gives about **10 remeshes** over
   the 100 mm flight, at about 42 s each (fact 41: 18 s fixed plus 76 µs per tetrahedron at 310 k), i.e. about
   7 minutes on a 51-minute flight.
2. **Surface coarsened** — the largest live surface patch edge exceeds `1.5 × h_surface`. This is the symptom
   itself rather than a proxy, and catches non-uniform erosion the depth rule misses.
3. **Quality** — the worst radius ratio among elements touching the derived surface falls below 0.05.

## The pipeline

1. **Validity gate, before anything else.** Reject and stop (exit code 1, clear message) if any boundary edge is
   used by other than two triangles (non-manifold — an element hanging by a vertex), if the live set has more than
   one connected component (fragmentation), or if the enclosed volume is non-positive. **Detecting fragmentation is
   in scope; handling it is not.** Guessing would be worse than stopping; if a flight trips this, that is the
   measured trigger to design for it.
2. **Taubin smoothing** of the derived surface: λ = 0.53, μ = −0.55, 8 passes, each displacement clamped to
   `TAUBIN_CLAMP = 0.35` of the mean edge length so a node can never pass its neighbours. Purely local, hence no
   assumption whatever about the shape — any genus, any number of components, concave or convex. Measured effect:
   +0.283 % on the enclosed volume (fact 40).
3. **Exact volume correction**: one uniform offset along the area-weighted vertex normals, by bisection, back to the
   mass account. Measured residual 1.1·10⁻¹⁴.
4. **Generate.** Discrete surface → `classifySurfaces(π, True, True, π)` → `createGeometry()` → `addSurfaceLoop` →
   `addVolume` → `generate(3)` with Task 1's band size field, `Mesh.Algorithm3D = 10`, and **4 threads**
   (`General.NumThreads`, `Mesh.MaxNumThreads3D`; 8 buys nothing, fact 41). Do not bother tuning the algorithm or
   the distance-field `Sampling` — both were measured to make no difference.
5. **Correct gmsh's own loss.** gmsh re-triangulates the surface at its target size and loses **0.091 %** of the
   volume — 1.3 g per remesh on the 1.472 kg sphere, indefensible over ten remeshes. Remove it with the same uniform
   offset applied to the *new* mesh's boundary; here the offset is **outward**, which cannot invert an interior
   element.
6. **Transfer the fields.** Nodal temperature by P1 point location in the old mesh. Because pointwise interpolation
   is not conservative, compute the total energy `Σ M_i h(T_i)` and the total mass on both meshes and apply a
   uniform enthalpy offset so both match exactly. `φ_e` by volume-weighted averaging onto the new elements; film
   mass per patch and `pending_load` per node by nearest-patch and nearest-node assignment, then rescaled to their
   exact totals. **Every transfer reports its pre-correction residual** — that is what replaces the unqualified
   10⁻¹⁰ claim of facts 28 and 31.
7. **Reset the reference.** The new mesh's coordinates become the origin for Task 16's displacements, and `φ_e` on
   the new elements is what the next step reads.

## Interfaces

- `REMESH_MIN_STEPS = 20`, `TAUBIN_PASSES = 8`, `TAUBIN_LAMBDA = 0.53`, `TAUBIN_MU = -0.55`, `TAUBIN_CLAMP = 0.35`,
  `EDGE_TRIGGER_FACTOR = 1.5`, `QUALITY_TRIGGER = 0.05`, `GMSH_THREADS = 4`
- `check_validity(mesh) -> None | raises RemeshRefused`
- `taubin_smooth(points, faces, ...) -> points`
- `correct_volume(points, faces, target) -> points`
- `generate_from_surface(points, faces, h_surface, h_core, band, ramp) -> VolumeMesh`
- `transfer_fields(old, new, T, phi, film, pending) -> (fields, TransferResiduals)`
- `should_remesh(mesh, band, h_surface, steps_since) -> str | None` (the trigger name, for the history)
- `remesh(mesh, state) -> RemeshResult(mesh, fields, volume_residual, energy_residual, trigger)`

## Tests

1. **Validity gate** stops on a hand-built non-manifold set and on a two-component set, with distinguishable
   messages.
2. **Smoothing conserves volume** to 10⁻¹² after correction; the pre-correction figure is recorded, not asserted.
3. **Smoothing reduces roughness** — the mean dihedral-angle deviation across the surface falls by at least half.
4. **A remesh conserves mass and energy** to 10⁻¹² after correction, with pre-correction residuals printed so a
   regression is visible.
5. **A remesh of an unchanged sphere is near-idempotent**: mass, energy and volume within 0.2 %, and the element
   count within 15 %.
6. **Both backends agree across a remesh** — skfem and FEniCSx give the same melting run to fact 12's tolerance
   (10⁻¹⁰ in mass, 0 K in temperature). **Non-negotiable.** Field transfer across a remesh is exactly the kind of
   change that produces plausible wrong answers, and two independent implementations agreeing is the cheapest
   protection available. Marked `@pytest.mark.drama`-style so it is skipped when dolfinx is absent, but it must run
   in the reference tier.

## Acceptance criteria

1. Mass and energy exact to 10⁻¹² after correction; pre-correction residuals reported for every remesh.
2. No inverted elements in any generated mesh.
3. A remesh of the 100 mm sphere's eroded surface completes in under 90 s at 4 threads (fact 41 predicts ≈42 s at
   310 k tets; the margin is for a more convoluted late-flight remnant, whose `classifySurfaces` cost is the first
   thing expected to degrade — spec §12 question 4).
4. Peak memory under 2 GB.
5. The two-backend agreement test passes.

## Not in this task

Moving the solver's mesh between remeshes (Design A, deferred; spec §10.3 defines its trigger). Handling a
fragmented body. Anisotropic layers on a remeshed surface — fact 12's sensitivity says they are not load-bearing,
and an offset self-intersects in a concave dish.
