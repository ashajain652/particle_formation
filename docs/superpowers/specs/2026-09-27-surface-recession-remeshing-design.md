# Resolving the receding surface: derived-geometry node motion and periodic remeshing (Step 3 amendment) — Design

Date: 2026-09-27
Status: design for review, awaiting implementation plan
Amends: `2026-09-20-melt-spraying-design.md` §5 (mesh) and §10 (geometry and accounting), and therefore sub-plans
`01-mesh-active-set.md`, `09-melting-body.md`, `10-coupled-loop.md`, `13-cli-wiring.md`, `14-verification-runs.md`
and `15-documentation.md`. Nothing in Step 1 or Step 2 changes. Every number below was measured in the design
session of 2026-09-27 unless it cites an existing measured fact.

## 1. Purpose

Step 3 as planned represents the body's shape as a **binary** function of which elements are alive: an element
contributes its full-size triangle to the surface until its fill fraction `φ_e` crosses `PHI_DEATH`, whereupon the
whole element vanishes. The fill fraction therefore gives sub-element accuracy in **mass** and none at all in
**shape**. Two consequences follow, and both were already visible in the prototype:

1. **The fine surface mesh is consumed early.** The prism stack is 3.75 mm thick; measured fact 14 has the nose
   receding 77 mm. The skin covers 4.9 % of the recession, and measured fact 28's discussion records that by 82 s
   on the 100 mm flight "the coarse core *is* the surface" on the windward face.
2. **The surface is a staircase whose steps are as tall as they are wide.** A step one cell tall and one cell wide
   deviates from smooth by 45°, and the patch normal angle θ is what the Lees distribution and the film tangents
   read.

This design fixes both without refining the mesh, which measured cost makes unaffordable (§2, row 6). It adds:
a **derived surface** that moves continuously with the recession while the solver's mesh stands still; a **dense
band** defined by a gmsh size field, replacing the prism stack as the mechanism that keeps the surface resolved;
and **periodic remeshing** onto the smoothed receded shape, so the band is renewed as the front eats through it.

Deliverables: a surface geometry accurate below cell size at every macro step; a surface-coarseness diagnostic in
every run's history; and a remeshing path verified to conserve mass and energy exactly after correction.

## 2. Facts this design relies on

| Topic | Fact | Source |
|---|---|---|
| gmsh cannot refine locally | gmsh 4.15.2's only refinement entry point is `gmsh.model.mesh.refine()` — "refine the mesh of the current model by **uniformly** splitting the elements". No argument, no subset, no target size. Size fields drive generation from scratch, not refinement of an existing mesh, and nothing returns a parent-to-child map. Local adaptive refinement must therefore be written by hand or not at all. | measured 2026-09-27 |
| Subdivision cannot give an anisotropic surface layer | Bisection produces children geometrically similar to the parent, so reaching 0.25 mm from an 8 mm core tetrahedron needs five levels, about 32 768 children per parent. Adaptivity pays when the refined band is a small fraction of the domain; here the front sweeps the whole body. | this design |
| How coarse the core is | Measured on `sphere_R50.000mm_hs2.000mm_hc8.000mm.msh` (18 896 nodes, 87 632 tets), equivalent regular-tetrahedron edge by shell: 2.74 mm at 45–50 mm radius, 3.67 at 40–45, 4.44 at 35–40, 5.24 at 30–35, 6.65 at 20–25, 8.26 at 10–15, 8.32 at 0–5. Depth-averaged over the recession the present mesh recedes in jumps of **about 6 mm**, not 0.25 mm. | measured 2026-09-27 |
| The death boundary is already valid | On a deliberately ugly eroded body (8 174 live tets, 3 404 boundary triangles) every one of the 5 106 boundary edges is used by exactly two triangles and the area-weighted normal sum is 4.05·10⁻¹⁶ of the total area. The boundary of a live element set is watertight and manifold; no surface *reconstruction* is needed, only smoothing. | measured 2026-09-27 |
| Taubin smoothing is exactly correctable | Eight λ/μ passes with each displacement clamped to 0.35 of the mean edge length changed the enclosed volume by **+0.283 %**; a single uniform offset along the area-weighted vertex normal of −19 µm restored it to within **1.1·10⁻¹⁴**. Taubin is purely local, so it is topology-agnostic: any genus, any number of components, no shape assumption. | measured 2026-09-27 |
| gmsh meshes the smoothed body | Discrete surface → `classifySurfaces(π, True, True, π)` → `createGeometry()` took **0.1 s** and produced **2** surface entities (on an unsmoothed jagged surface the feature detector would produce thousands). Volume generation gave 27 635 nodes / 138 946 tets, **0 inverted**, volume −0.091 % of the smoothed target. | measured 2026-09-27 |
| Remesh cost | All of it is `generate(3)`; geometry and size-field setup are ≤ 0.1 s. HXT and Delaunay are indistinguishable; distance-field `Sampling` 200 vs 20 makes no difference. Threads: 45.8 s (1), 29.4 s (4), 29.7 s (8) — saturates at 4, a 1.56× speed-up. Scaling at 4 threads: 137 248 tets 29.5 s, 376 183 tets 48.7 s, 748 606 tets 75.1 s, i.e. ~18 s fixed plus ~76 µs per tetrahedron, **markedly sub-linear**. Peak memory 0.5 GB at 749 k tets. | measured 2026-09-27 |
| Mesh size vs flight cost | 100 mm flight, 1 191 macro steps, scaling the 2.0 s step of measured fact 23: 20 mm band at 2.16 mm = 365 k tets, 8 remeshes, **54 min**; **15 mm band at 2.16 mm = 310 k tets, 10 remeshes, 51 min**; 15 mm at 1.5 mm = 910 k, **1 h 46**; 15 mm at 1.0 mm = 3.07 M, **5 h**. Today's baseline is 39 min. Refining cells scales as the cube; the band thickness only sets the remesh count, which is cheap. | estimated 2026-09-27 from measured remesh and step costs |
| Node motion and element quality | Pushing every boundary node of the 2.74 mm mesh inward: 0 inverted to 1.25 mm; worst-element radius ratio 0.289 → 0.232 (0.5 mm) → 0.076 (1.0 mm) → 0.022 (1.25 mm); at 1.5 mm two elements fall below 5 % of their volume; at **2.0 mm, 1 319 elements invert** and 4 908 are slivers. Safe to about **40 % of the local cell**, collapse past about 70 %. Mean quality near the surface degrades far more slowly than the worst element, so the kill rule must be per element. | measured 2026-09-27 |
| Per-patch volume cannot be matched exactly | Euler forces a closed triangulation to have twice as many triangles as vertices: measured **18 078 patches against 9 041 boundary nodes**, and 2·9 041 − 4 = 18 078 exactly. Per-patch volume constraints outnumber nodal unknowns 2:1, so nodal motion generically cannot satisfy them all. | measured 2026-09-27 |
| …but the least-squares fit is excellent when the field is smooth | Solving for nodal displacements that best reproduce each patch's required volume loss: a **smooth** (Lees-like, √cos θ) recession field gives median per-patch error **0.03 %**, 90th percentile 0.56 %, total +0.000 %, leaving `φ_e` a residual of 10⁻⁴ of an element. A **patchy** field (spray gate flipping) gives median 22 % and a residual of 7.3 % of an element; white noise 29 %. The **total** is within 0.017 % before correction in every case and exactly zeroable after. | measured 2026-09-27 |
| Frozen geometry costs about a factor of two in the conduction | 1-D ablating AA7075 bar (ρ 2813, k 170, c_p 1100 from `AA7075_nomelt.json`), backward Euler, P1, lumped capacity, prescribed recession and conducted flux, 2.16 mm against a 0.05 mm reference. Typical (0.40 MW/m², 0.10 mm/step, surface 312 → 600 K): moving nodes err by mean +0.89 K / max 2.43 K, frozen geometry by mean **+1.97 K** / max 6.00 K. Peak flux (2.00 MW/m²): +4.44/12.14 against **+9.87/29.99**. Peak recession (1.00 mm/step): +0.62/4.13 against **+1.59/7.48**. Frozen geometry is consistently ~2.2× the mean and ~2.5× the maximum error, and **both are dominated by the 2.16 mm cell**, not by the choice. | measured 2026-09-27 |
| The death threshold is a free mitigation | Frozen-geometry mean error against `PHI_DEATH`: 0.05 → +2.97 K (typical) / +14.83 K (peak flux); 0.25 → +2.27/+11.34; **0.50 → +1.97/+9.87**; 0.75 → +2.37/+11.84. A shallow optimum at a half; the plan's current 0.05 is the worst of the four. | measured 2026-09-27 |
| Backend behaviour under a moving mesh | skfem's speed comes from rescaling matrices precomputed in `element_matrices`, valid only while the geometry is fixed: 87 632 tets assemble in **0.004 s** fixed but **0.074 s** (17.4×) if the geometry must be recomputed. dolfinx recomputes Jacobians every assembly anyway, so moving **every** node changed its time not at all: **0.031 s** either way. At 310 k tets: skfem 0.015 s fixed, 0.262 s moving; dolfinx 0.11 s always. `mesh.geometry.x` is writable in place; `dolfinx.fem.create_interpolation_data`, `geometry.bb_tree` and `compute_colliding_cells` exist in the installed 0.11. | measured 2026-09-27 |
| Prism layers are not load-bearing | Measured fact 12's sensitivity varied the layer count 2/4/6 — total fine radial depth 0.75 mm to 15.75 mm — and moved the sprayed mass by **< 0.2 %** and the demise altitude by **< 0.3 %**. The stack costs 15 430 × 3 × 4 = **185 160 tets**, 73 % of the present mesh, for 3.75 mm of depth. | measured fact 12; fact 6 |
| The body is nearly isothermal in depth | Thermal diffusivity ≈ 5.5·10⁻⁵ m²/s at 733 K, so heat travels ≈ 5 mm in one 0.5 s macro step — more than the whole prism stack — and ≈ 230 mm over the flight, four times the radius. There is no thin thermal boundary layer inside the solid to resolve; the resolution requirement is geometric, not thermal. | this design, from `AA7075_nomelt.json` |
| Interior drainage cannot be geometric | Measured fact 5 as amended by 28(b) has interior elements feeding liquid to the four nearest patches from up to 2 mm below the surface. That mass leaves a non-boundary element and cannot be represented by moving nodes without opening a cavity, which `PHI_MIN` = 10⁻³ forbids. Measured fact 28 finds surviving elements 94–99 % full, so it is a small term. | facts 5, 28(b), 28 |

## 3. Scope

**In.** A derived surface (moved boundary node positions) computed each macro step from the mass removed, used by
every consumer of the shape; smoothed normals for the heating, film tangents and nose-radius fit; a dense band set
by a gmsh size field; periodic remeshing with Taubin smoothing, exact volume correction and a validity gate;
conservative field transfer onto the new mesh; `PHI_DEATH` raised to 0.5; new diagnostics and CLI flags; the
verification runs and documentation re-measured.

**Out.** Moving the solver's own mesh (Design A of §5; deferred, with its trigger defined in §10.3). Topology
change — a body that fragments is detected and the run stops (§7.4). Anisotropic surface layers after the first
remesh. Lateral refinement: the patch count and therefore the angular resolution of the heating is unchanged.
Tumbling, which measured fact 24 already records as the next iteration's first item and which is a larger
uncertainty on the shape than anything here.

## 4. Architecture

Two new modules and edits to four existing ones. The guiding principle is that **the solver is not touched**:

- `recession.py` (new) — converts the mass removed per patch into nodal displacements (least squares over the
  patch-volume constraints), maintains the derived surface node positions, and provides the smoothed normals.
- `remesh.py` (new) — manifold check, Taubin smoothing with displacement clamp, exact volume correction, gmsh
  volume generation from the discrete surface, and the conservative field transfer.
- `mesh.py` — the dense band size field replaces the prism stack as the default for melting runs; `VolumeMesh`
  gains a `derived_points` array and `surface(derived=True)`; prism layers stay available but default to 0.
- `body.MeltingBody` — `PHI_DEATH` 0.05 → 0.50; the recession solve and the derived surface are updated in the
  melt step; `nose_radius`, `newtonian_drag`, `reference_area` and `drag_shape_factor` read the derived surface.
- `coupled.py` — the remesh trigger is evaluated once per macro step; new history columns.
- `cli.py`, `viz.py` — flags and VTK output of the derived surface.

Per macro step the melt step gains one stage, placed **after** element death and **before** the history row:
recession solve → derived surface update → smoothed normals → remesh trigger check → (remesh if triggered).

## 5. The central decision: `φ_e` stays authoritative, the geometry is derived

Two coherent designs exist. **Design A** lets the geometry carry the surface part of the state: element volumes
change every step, the solver's matrices are rebuilt, and `φ_e` becomes a residual. **Design B** keeps `φ_e`
authoritative against the *original* element volumes — the thermal core is untouched, the capacity matrix stays
`diag(Σ φ_e V_e(0)/4 ρ c_i)`, the mass stays `Σ φ_e ρ V_e(0)`, the energy balance is the one already verified —
and the moved node positions are a derived rendering read only by the surface physics.

**This design adopts Design B**, for four reasons. The thermal core is where every recorded failure in this project
lives (measured facts 3, 4, 25, 26: a Newton iteration cycling across the latent-heat ramp, nodes driven to 5000 K,
a period-three limit cycle on 22 nodes, deferred loads worth 10⁴ K). The exact energy balance would otherwise have
to be re-derived to account for material sweeping across a moving mesh, which is the standard hard part of
moving-mesh formulations. Design B keeps skfem's 0.015 s assembly instead of 0.262 s, so the backend question
(§2, last row) does not arise. And the displacement solver it needs is *the same code* Design A would need, so it
is a stepping stone, not a detour.

**What Design B gets wrong, stated plainly.** The capacity is right — a half-consumed element stores half the heat.
The **conductance** is not: the remaining material is a thinner slab and should conduct about twice as well, but
the solver sees the full original thickness, so the surface sheds heat inward too slowly and runs hot. Measured
fact 4's decision to leave `k` unscaled is correct as far as it goes but does not fix the path length. The
measured size of this is in §2: mean **+1.97 K** in the typical regime against a **±2 K** feed ramp, rising to
+9.87 K at peak flux — where Design A is also past the ramp at +4.44 K, so the *mesh* and not the design is the
limit. Raising `PHI_DEATH` to 0.5 removes about a third of it for free.

**Where `φ_e` and the geometry meet.** Mass = Σ_e ρ `φ_e` V_e(0) is unchanged and remains the only mass account.
The derived surface is fitted to it, and because the fit is overdetermined 2:1 it matches per patch only
approximately (§2: median 0.03 % for a smooth field). The **total** is forced to agree exactly by the uniform
offset correction. The per-patch mismatch is recorded as a diagnostic, not hidden.

## 6. The dense band replaces the prism stack (`mesh.py`)

The prism stack is retired as the default for melting runs. It costs 185 160 tetrahedra — 73 % of the present mesh
— for 3.75 mm of depth, it is consumed in the first 5 % of the recession, and measured fact 12 shows the answer
moves by under 0.3 % when its depth is varied twentyfold. It also cannot survive remeshing: a prism layer is an
*offset*, which self-intersects once the offset exceeds the local concave radius of curvature, and the windward
face becomes a dish.

In its place, the existing `Distance` + `Threshold` size field gains a non-zero `DistMin`:

```
SizeMin  = h_surface      (2.16 mm)
SizeMax  = h_core         (8 mm)
DistMin  = band_thickness (15 mm)      <- was 0.0
DistMax  = band_thickness + ramp
```

which reads: keep cells at `h_surface` everywhere within 15 mm of the surface, then grow them. A distance field is
well defined at every interior point of *any* shape, concave included, because nothing is being offset — which is
why it survives remeshing where prisms do not.

**Defaults:** `DEFAULT_BAND = 15e-3`, `DEFAULT_LAYERS = 0` for melting runs. `--prism-layers` remains for the
pinned verification devices (measured fact 12 already runs the bookkeeping device with `--prism-layers 0`).
Cost: about **310 000 tetrahedra**, against 255 276 today, and a 51-minute flight against 39. Refining the band
cells to 1.0 mm would cost 3.07 M elements and five hours (§2), which this design does not do; the recession
quantum is instead removed by §7 rather than reduced by refinement.

## 7. Remeshing (`remesh.py`)

### 7.1 Triggers

Evaluated once per macro step; the first to fire wins, subject to a minimum gap of `REMESH_MIN_STEPS` (default 20)
so they cannot thrash.

1. **Band consumed.** At each remesh every element is tagged with its distance to the boundary. Fire when the
   maximum tagged depth over elements that currently own a surface patch exceeds `band_thickness / 2`.
2. **Surface coarsened.** Fire when the largest live surface patch edge exceeds `1.5 × h_surface`. This is the
   symptom itself rather than a proxy and catches non-uniform erosion the depth rule misses.
3. **Quality.** Fire when the worst radius ratio among elements touching the derived surface falls below 0.05.

With a 15 mm band and the 77 mm nose recession of measured fact 14, trigger 1 gives about **10 remeshes** over the
100 mm flight at about 42 s each (18 s fixed plus 76 µs × 310 000, from §2), i.e. about 7 minutes.

### 7.2 Smoothing, and why it accommodates an arbitrary shape

The input is the derived surface: the live element set's boundary, with each boundary node already displaced by
the recession (§8). It is closed and manifold by construction (§2, row 4), so nothing is reconstructed.

**Taubin λ/μ smoothing** is applied: a shrinking pass with λ = 0.53 and an inflating pass with μ = −0.55, eight
times. It looks only at a node and its immediate neighbours, so it makes **no assumption whatever about the
shape** — any genus, any number of components, concave or convex. That locality is what delivers the arbitrary
shape requirement; it is not a property that has to be engineered in.

Each node's displacement is clamped to `TAUBIN_CLAMP` = 0.35 of the mean edge length, so a node can never pass its
neighbours and self-intersection is prevented in practice.

**Exact volume correction.** Smoothing moved the volume by +0.283 % in the measured case. Every node is then
displaced along its area-weighted vertex normal by a single distance found by bisection so that the enclosed
volume exactly matches the mass account; measured residual 1.1·10⁻¹⁴.

A **level-set / marching-cubes** reconstruction was considered and rejected. It is more general in the wrong
direction: it can merge, split and fill, so at any practical grid spacing it would weld a thinning lens shut or
erase a fragment — inventing and destroying exactly the features the physics is trying to predict. A
**radial or spherical-harmonic** parametrisation was rejected because it assumes a star-shaped body.

### 7.3 Generation and transfer

The smoothed surface goes to gmsh as a discrete entity, is reparametrised with `classifySurfaces` and
`createGeometry` (0.1 s, 2 entities — smoothing is what keeps this tractable), and the volume is generated with
the §6 size field at 4 threads (`General.NumThreads`, `Mesh.MaxNumThreads3D`; 8 buys nothing).

gmsh re-triangulates the surface at its own target size and so loses **0.091 %** of the volume. That is 1.3 g per
remesh on the 1.472 kg sphere and would accumulate over ten remeshes into something indefensible. It is removed by
the same uniform-normal-offset bisection applied to the *new* mesh's boundary — here the offset is outward, which
cannot invert an interior element.

**Field transfer.** Nodal temperature is interpolated by point location in the old mesh (P1, barycentric).
Because pointwise interpolation is not conservative, the total energy `Σ M_i h(T_i)` and the total mass are
computed on both meshes and a uniform enthalpy offset is applied so that both match exactly. `φ_e` is transferred
by volume-weighted averaging onto the new elements; film mass per patch and `pending_load` per node are
transferred by nearest-patch and nearest-node assignment and then rescaled to their exact totals. **Every transfer
reports its pre-correction residual**, which replaces the current unqualified 10⁻¹⁰ claim (§10.2).

### 7.4 Validity gate

Before smoothing: reject and stop with a clear message if any boundary edge is used by other than two triangles
(non-manifold, e.g. an element hanging by a vertex), if the live set has more than one connected component
(fragmentation), or if the enclosed volume is non-positive. Detecting fragmentation is in scope; *handling* it is
not, and guessing would be worse than stopping. If a flight trips this gate in practice, that is the measured
trigger to design for it.

## 8. The derived surface (`recession.py`)

Each macro step, after element death. The displacement is a **state function of `φ_e` and the active set**,
recomputed from scratch every step rather than accumulated, so it cannot drift away from the mass account.
If the displacement a step would apply exceeds `RECESSION_SUBSTEP_FRAC` = 0.35 of the local cell size (§2, the
safe envelope), the step is sub-stepped and the fact recorded, so an unmeasured peak recession cannot invert
elements (§12, question 2):

1. **Consumed depth per patch.** For each live surface patch, the volume consumed in its column is the volume of
   the dead elements behind it plus the fill deficit `(1 − φ_e) V_e` of the live ones it shadows, giving a target
   volume `ΔV_p` for that patch. This reads `φ_e` directly and is the *only* place the consumed depth is derived,
   so there is no second sub-cell correction to apply and no possibility of double counting.
2. **Nodal displacements.** Solve, in least squares, `M s = ΔV` where `s_i` is node `i`'s displacement from its
   **original** position along its area-weighted vertex normal and `M[p, i] = (A_p / 3)(n_i · N_p)`. The system is
   overdetermined 2:1 (§2) and generically inconsistent; the residual is expected, and is reported rather than
   suppressed. Because a half-consumed live element contributes its deficit in step 1, the surface is drawn half
   way through that element — which is what removes the recession quantum entirely rather than merely reducing it.
3. **Write** into `derived_points`, a separate array. `VolumeMesh.points` is never modified, which is what keeps
   the solver untouched.
4. **Smoothed normals.** A Taubin-smoothed copy of `derived_points` supplies the patch normals, θ and tangents.
   This is one sparse matrix-vector product over 15 430 patches and it is what removes the 45° staircase slopes.
   *Applied to today's mesh it would already be an improvement*, independently of everything else here.

After a remesh the reference resets: the new mesh's own coordinates become the origin for `s`, and `φ_e` on the new
elements (§7.3) is what the next step's step 1 reads.

`PHI_DEATH` rises from 0.05 to **0.50** (§2, row 12), and the residual mass of a dying element joins the film
exactly as measured fact 4 specifies. Note this hands the film half an element at a time rather than a twentieth:
the film cap of measured fact 25 and `film_blob_fraction` must be checked against it (§10.1).

## 9. Diagnostics and CLI

New history columns: `surface_edge_max_mm`, `surface_edge_mean_mm` (the coarseness the run actually had),
`recession_residual_frac` (median absolute per-patch mismatch from §8 step 2), `phi_surface_min`,
`conduction_dT_err_K` (`q L (1 − φ_e) / k` per surface element, the §5 error made visible), `remesh_count`,
`remesh_volume_residual_frac`, `remesh_energy_residual_frac`.

New flags: `--band-thickness` (default 15 mm), `--remesh {auto,off}` (default `auto` for melting runs, `off` for
the pinned verification devices), `--remesh-min-steps`, `--derived-surface {on,off}` (`off` reproduces the current
staircase behaviour, for the A/B comparison of §10.1).

Run names must continue to encode the configuration: append `_band15mm_remesh` (or `_noremesh`) per the repository
convention, so runs with and without this machinery cannot overwrite each other.

## 10. Verification and acceptance

### 10.1 Acceptance tests, in the order they should be written

1. **Unit — the boundary is manifold.** On the eroded fixture, every edge used exactly twice; closure below
   10⁻¹² of the area. (Measured 4.05·10⁻¹⁶.)
2. **Unit — smoothing preserves volume exactly.** Taubin then the offset correction: volume residual below 10⁻¹².
   (Measured 1.1·10⁻¹⁴.)
3. **Unit — the recession solve.** On a smooth synthetic field, median per-patch volume error below 0.5 % and total
   below 10⁻¹² after correction. (Measured 0.03 % median.) On a patchy field, the test asserts only that the total
   is exact and that the reported residual is non-zero — the point is that it is *measured*, not that it is small.
4. **Unit — the validity gate** stops on a hand-built non-manifold set and on a two-component set.
5. **Integration — a remesh conserves.** Mesh, remesh, and assert mass and energy unchanged to 10⁻¹² after
   correction, with the pre-correction residuals recorded in the test's output so a regression is visible.
6. **Integration — both backends agree.** skfem and FEniCSx must give the same melting run across a remesh to the
   tolerance of measured fact 12 (10⁻¹⁰ in mass, 0 K in temperature). Moving-mesh and transfer code is exactly the
   kind of change that produces plausible wrong answers, and two independent implementations agreeing is the
   cheapest protection available. **Non-negotiable.**
7. **Flight — the shape A/B.** The 50 mm physics flight with `--derived-surface off` against `on`, comparing
   sprayed mass, demise time and median droplet radius. This is the measurement that says what the whole design
   bought.
8. **Flight — the film can absorb early death.** With `PHI_DEATH` at 0.5, confirm `film_blob_fraction` and
   `film_frozen_fraction` stay within the ranges measured fact 25 records, and that the stranded-film mass stays
   below the 0.47 % of initial mass measured there.

### 10.2 What the thesis may no longer claim unqualified

Measured facts 28 and 75 quote balance residuals of −1.7·10⁻¹⁰ and +2.2·10⁻¹⁰. After remeshing these become
**exact by construction following an explicit correction**, with a pre-correction residual per remesh that must be
reported. That is the standard position for an ablation code and it is defensible, but it is a weaker statement
than the present one and `README.md`, `docs/model_assumptions.md` and the Task 15 blocks must say so.

### 10.3 The measured trigger to reconsider Design A

Record `conduction_dT_err_K` over the flight. If it stays well under the ±2 K feed ramp for the great majority of
melting steps, Design B is sufficient and the chapter states it. If it sits near ten kelvin for a substantial
fraction, that is the measured reason to move the solver onto the current geometry. The **better** test, which
§10.1.7 does not cover, is the **melt-rate** error under a surface pinned by the latent heat: while melting, a
temperature error becomes a mass-flux error rather than a temperature error, and mass flux is what feeds the
spraying. The 1-D probe of §2 cannot measure it; the real solver can, and it should be the first thing the
implementation reports.

## 11. Assumptions to state in the thesis

1. The surface geometry is fitted to the mass account, not identical to it; per-patch agreement is median 0.03 %
   for a smooth recession field and degrades to tens of per cent if the field is patchy. The total always agrees
   exactly. The realised value is reported per run.
2. The conduction solves on the original element geometry, so heat crosses material that has been removed. Measured
   penalty: mean +1.97 K at the surface in the typical regime, +9.87 K at peak flux, roughly twice what a moving
   mesh would give, and in both cases smaller than the error the 2.16 mm cell size itself contributes.
3. Smoothing rounds genuinely sharp features as well as staircase noise. The displacement clamp bounds this; the
   rim of a thin lens will lose some of its edge.
4. Remeshing conserves mass and energy only after an explicit correction, whose pre-correction size is reported.
5. A body that fragments is not modelled; the run stops.
6. Lateral resolution is unchanged, so the angular resolution of the heating and the film graph is still set by the
   2.16 mm patch size.
7. All of it remains a fixed-attitude calculation, which measured fact 24 already identifies as the dominant
   uncertainty on the shape.

## 12. Open questions, to be answered by the implementation

1. **How smooth is the real recession field?** §2 shows the answer swings between excellent and mediocre on this
   alone. Measurable directly from the prototype and it should be the first number obtained.
2. **What is the peak recession rate?** Partly resolved 2026-09-27. Measured fact 5's "at peak heating 3.6
   layers melt per step" is a statement about the **rejected** surface-only feed rule and about how much material
   crosses the **liquidus**, not about how much *leaves*. Since fact 28(b) gates the feed by depth — and fact 28
   records the stagnant molten inventory held in place rising to about 1 % of the body — melting rate and recession
   rate are different quantities, and the recession is set by what sprays or runs off. An order-of-magnitude
   estimate from the 50 mm flight's 0.173 kg sprayed over roughly 300 melting steps and a windward area of order
   2·10⁻³ m² gives **about 0.1 mm per step**, which is under 5 % of a cell and comfortably inside the safe envelope
   of §2. The **peak** remains unmeasured, so §8 must sub-step the node motion whenever the step's recession would
   exceed `RECESSION_SUBSTEP_FRAC` (default 0.35) of the local cell, and report when it does. That makes the design
   safe under either answer instead of blocking on it.

3. **What does the least-squares solve cost per step?** Three solves completed in seconds in the probe but were
   not timed in isolation. If it is significant, area-weighted averaging of incident patch recessions is a
   one-matvec fallback whose accuracy has not been measured.
4. **Does `classifySurfaces` stay cheap on a late-flight remnant?** It cost 0.1 s on one moderately eroded body.
   A more convoluted shape is the first thing expected to degrade.
5. **Does the 0.091 % gmsh volume loss hold at production size?** Measured once, on one body, at 1 mm cells.

## 13. Future-iteration items (recorded so they are not lost)

- Design A: the solver on the current geometry, with the ALE energy bookkeeping done properly. §10.3 defines when.
- Fragment tracking: a detached piece as its own body with its own trajectory.
- Anisotropic bands after a remesh, via a curvature-limited offset, if the melt front ever needs sub-millimetre
  normal resolution that §2's evidence says it does not.
- FEniCSx as the primary backend if the mesh grows or Design A is adopted: it assembles a moved mesh for free
  (0.031 s either way) and already ships `create_interpolation_data` for the transfer of §7.3.
