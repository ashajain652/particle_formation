# Synthetic finite-element runs for `spheral_frag` (Spheral M1, Task 2)

Small run directories in the frame contract's format (`spheral_frag/contract.py`, `FE_FIELDS`), written by
`tests/spheral_frag_synthetic.py` through `make_fixtures.py`:

```bash
"$PY" tests/fixtures/spheral_frag/make_fixtures.py      # rewrites the three directories; deterministic
```

Each directory is a finite-element output directory holding one run, laid out as the export writes it:
`<run>.csv` (history, 9 significant digits), `<run>.json` (run JSON), `<run>/vtk/field_<k>.vtu`,
`<run>/vtk/surface_<k>.vtp` and `<run>/vtk/{field,surface}.pvd` (times printed `{:.6g}`). The vtu and the vtp hold
the same node array (every mesh node, unreferenced ones included); the vtu's cells are the active tetrahedra; the
vtp's triangles are wound outward. Every node, tet and patch field of the contract is present, optional ones
included. The history carries the contract's columns plus two it does not name (`surface_T_max_K`,
`n_active_elements`), which a reader must ignore. Total size 664 kB (sphere 353 kB, dumbbell 221 kB, slab 90 kB).

The exact answers below are also available in code: the constants and the `*_answers`, `*_patches` and
`frame_*` functions of `tests/spheral_frag_synthetic.py`. Values marked "measured" are properties of the
committed meshes (gmsh 4, drama_env, 2026-10-07), not closed forms.

## Conventions shared by all three

- **Material `synthetic_linear_fl`** (the run JSON's `settings.material`): solid 2,813 kg/m³, liquid 2,400 kg/m³,
  f_l **linear in T from 0 at 750 K to 1 at 908 K**; nodal `liquid_fraction` equals that law exactly. It is not
  AA7075_scheil, so a consistency check against the real material table does not apply to these fixtures.
- **Loads (decision 4 of the M1 review):** `p_w`, `tau` and `closure` are written on every frame whose step evaluated
  the surface flow, film or not. In the patch's own inclination θ (outward facet normal against v_hat = +x) and with
  p_stag = the row's `p_w_stag_Pa` as written in the CSV:
  p_w = p_stag cos²θ for θ < 90°, 0.01 p_stag for 90° ≤ θ < 150°, and **0 beyond 150° ("not evaluated")**;
  tau = 0.02 p_stag sin 2θ for θ < 90°, 0 elsewhere. The history's drag is the frame's own:
  `mass_kg · load_factor_g · g0` = Σ A (p_w cos θ + tau sin θ), g0 = 9.80665 m/s².
- **Girin's closure** (`closure` = 0, `delta_m` finite) on θ < 60° of melting frames, Couette (1, `delta_m` NaN)
  elsewhere: `delta_m` is finite exactly where `closure` is 0.
- **Mass:** `mass_kg` = Σ φ ρ V + Σ film_thickness ρ_l A + Σ deep_thickness ρ_l A of the frame (to the CSV's 9
  digits, i.e. within 5e-9 relative — half a unit of the 9th digit, not 2e-9); `inputs.mass_kg` is the analytic
  solid's ρV, as the finite-element CLI writes it, and `results.initial_mass_kg` the meshed frame 0's.
- **The run JSON** has every contract run item at its dotted name, including `v_hat` = [1, 0, 0] at the top level
  (the contract's name; decision 4 has the export write it), and a `synthetic` block with the device's parameters.

## `sphere_run/` — `synthetic_sphere_d100.00mm_seed12345`

A coarse 100 mm sphere (`reentry_model.mesh.sphere_mesh(0.05, 8e-3, 20e-3, band=0)`): 753 nodes, 2,497 tetrahedra,
1,194 patches. Three frames, one history row each:

| k | t [s] | p_w_stag_Pa | state | T [K] (r = distance from the centre) |
|---|---|---|---|---|
| 0 | 0 | 1552.5 | written before any step: nothing evaluated — p_w = tau = 0, closure 1, delta_m/kn_local/r_droplet NaN, no film, no release | 300 uniform |
| 1 | 0.5 | 1715.32034 | evaluated, film FILM0 cos θ (5e-5 m), release 1e-3 cos θ kg/m² per step on θ < 60° | 300 + 600 (r/R)² |
| 2 | 1.0 | 1895.2166 | as 1, with ten nose elements dead and two partly consumed | 300 + 660 (r/R)² |

Frame 2: the dead elements are the owners of the ten patches nearest θ = 0 (indices into frame 1's tetrahedra:
1404, 678, 1016, 770, 1413, 441, 12, 1182, 1878, 1063), leaving 2,487 tetrahedra, 1,198 patches (a staircase at the
nose) and one node no active tetrahedron uses (752 of 753 referenced). φ = 0.5 and 0.25 on frame 2's tetrahedra
1362 and 429 (owners of the patches nearest θ = 90°). Every frame is closed as an oriented surface, every edge in
exactly two triangles, and its divergence volume equals its tetrahedra's volume to 1e-12.

Answers:

- **Thickness** 2R = 0.1 m. Facet centroids lie between 0.25 % and 0.70 % of R inside the sphere (measured:
  centroid radius 49.651–49.878 mm), so the centroid-to-opposite-facet distance is 2R minus up to about 1 %;
  Task 5's "below 0.5 %" expectation is for a finer sphere. Measured by Task 5: 99.184–99.841 mm, ε_sphere 0.816 %
  (mean 0.369 %), inside the faceting bound 1 − r_in/R = 0.907 % (r_in = 49.546 mm, the nearest facet plane).
- **Drag** (= `D_patch` = `D_hist`): frame 1 6.866822220 N, frame 2 7.311101107 N; row 0 records 6.215014915 N,
  the drag the analytic loads would give, while frame 0 carries none (a frame without loads).
  Continuum: D = π R² p_stag (1/2 + 0.02 − 0.75 · 0.01) = 0.5125 π R² p_stag; the faceted sphere gives 0.509707
  (−0.55 %), its windward pressure part 0.99454 of π R² p_stag / 2 (measured, frames 0–1 geometry).
- **Loads table:** the last bin with p_w > 0 is the one below 150°; bins 150–180° carry nothing.
- **Mass** (CSV): 1.45904958, 1.45998733, 1.45755883 kg; Σ ρV = 1.459049578 kg (frames 0–1); film 0,
  9.37747047e-4, 9.38696106e-4 kg; deep 0; `sprayed_mass_kg` 0, 5.87115569e-6, 1.1702366e-5 kg (each step's
  increment is Σ release_rate · A of that frame). Analytic `inputs.mass_kg` 1.472883356 kg.
- **Girin patches:** 306 (frame 1), 304 (frame 2).
- **Depths of the radial T** (continuum): f_l = 0.5 at 829 K, at depth R(1 − √(529/(T_wall − 300))) = 3.051 mm
  (frame 1) and 5.236 mm (frame 2); liquid (908 K) to 2.010 mm on frame 2, none on frame 1. P1 on 8–20 mm elements
  is far from these; Measured by Task 6
  (`spheral_frag.geometry.layer_depths`, ds 0.05 mm): slurry depth 3.662–9.029 mm on frame 1 (+0.61 to +5.98 mm)
  and 4.211–9.635 mm on frame 2's 1,184 patches on the sphere (−1.02 to +4.40 mm; the 14 staircase patches
  excluded), all within one surface element (8 mm); liquid depth 0 on every patch (the first interior nodes are
  under the liquidus, so P1 holds no liquid below the surface).
- **Centre crossing 1 → 2:** exactly the ten dead elements' interiors leave the body.

## `dumbbell_frame/` — `synthetic_dumbbell_R20.0mm_neck4.0mm`

Two spheres R = 20 mm centred at x = ∓30 mm, fused (gmsh OCC) with a cylinder r_neck = 4 mm from centre to centre,
so the free neck is L_neck = 20 mm long; size 1.5 mm within 3 mm of the neck, 5 mm elsewhere. 1,551 nodes, 6,042
tetrahedra, 2,068 patches; one frame at t = 0, evaluated, no melting (film 0, closure 1, delta_m NaN, T = 600 K).

- **Thickness:** 2 r_neck = 8 mm on the neck's 567 lateral patches (`dumbbell_neck_patches`: every vertex on
  r = 4 mm within 1e-9 m; longest edge 2.04 mm; the free lateral surface runs to the junction circles at
  |x| = 10.404 mm, so it is 20.808 mm long at r = r_neck, not L_neck), 2R = 40 mm on the 1,471 sphere patches whose
  inward ray exits the same sphere (`dumbbell_sphere_patches`: antipode more than 10° outside the neck's cap,
  half-angle asin(r/R), neck patches excluded). Measured by Task 5 (`spheral_frag.geometry.thickness_map`): neck
  7.879–7.992 mm; spheres 35.65–39.98 mm, the short ones under three 13 mm facets at the left sphere's pole
  (−30, 0, 20) mm — every patch within the faceted bounds of `tests/test_spheral_frag_geometry.py`.
- **thin_patches:** at dx = 2.2 mm exactly the neck's 567 patches (8 mm < 4 × 2.2 = 8.8 mm); at dx = 1.5 mm none
  (6 mm < 8 mm). (Before Task 5 this README said 495 neck and 1,484 sphere patches: the old selection kept only
  centroids with |x| ≤ L_neck / 2 and counted 13 neck patches as sphere patches.)
- **Volume:** exact 6.80462e-5 m³ (`dumbbell_volume`); meshed 6.66494e-5 m³ (−2.05 %, faceting).

## `slab_frame/` — `synthetic_slab_L20.0mm_H12.0mm`

A box 20 × 4 × 12 mm (x, y, z) of 1 mm Kuhn-split cubes (`reentry_model.mesh.box_mesh`): 1,365 nodes, 5,760
tetrahedra, 1,472 patches (the whole closed boundary). The exposed surface is the top face z = H (160 patches); the
film, delta_m and deep fields are set there only. Five columns along x, four cells each (the fifth also takes the
x = 20 mm plane); a node's column is that of its x plane, min(i // 4, 4). Nodal
f_l = 0.5 + (D_c − d)/24 mm at depth d = H − z (inside (0.0063, 0.667): never clipped, no liquid), T from the linear
law, so f_l = 0.5 exactly at depth D_c under P1 below every **exact** top patch (both x planes of its cube in one
column: 24, 24, 24, 24 and 32 patches; `slab_exact_patches`).

| column | x [mm] | slurry depth | liquid depth | film | delta_m | deep_thickness | layer (film + slurry) | with deep |
|---|---|---|---|---|---|---|---|---|
| 0 | 0–4 | 0.15 mm | 0 | 0.1 mm | 0.3 mm | 0.2 mm | 0.25 mm | 0.45 mm |
| 1 | 4–8 | 1.0 mm | 0 | 0.1 mm | 0.3 mm | 0.2 mm | 1.1 mm | 1.3 mm |
| 2 | 8–12 | 2.5 mm | 0 | 0.1 mm | 0.3 mm | 0.2 mm | 2.6 mm | 2.8 mm |
| 3 | 12–16 | 4.0 mm | 0 | 0.1 mm | 0.3 mm | 0.2 mm | 4.1 mm | 4.3 mm |
| 4 | 16–20 | 1.0 mm | 0 | 0.1 mm | NaN (closure 1) | 0.2 mm | 1.1 mm | 1.3 mm |

Zones (M1 plan, Task 6): skin = min(film + liquid depth, δ_m) = 0.1 mm in columns 0–3, NaN in column 4; zones 1, 2,
3, 3 at film limit 2 mm and 1, 2, 2, 3 at 3 mm; bulk 0, 0, 0.6, 2.1 mm at 2 mm; column 4 zone 2; under-resolved at
dx = 1.1 mm for the 2.6 and 4.1 mm layers, at dx = 0.5 mm for neither. Both zone-3 readings (decision 6): the
whole layer (2.6 and 4.1 mm) or the part below the limit (`bulk`); under the second, dx = 0.5 mm flags column 2
(0.6 mm < 2 mm) and not column 3. With `include_deep` the layers are 0.45, 1.3, 2.8, 4.3, 1.3 mm, so column 0
(0.45 mm > δ_m) becomes zone 2: zones 2, 2, 3, 3, 2 and bulk 0, 0, 0.8, 2.3 mm at 2 mm. Measured by Task 6: every
exact patch's slurry depth equals its column's to 5e-18 m. Mass: solid 2.70048e-3 kg, film 1.92e-5 kg,
deep 3.84e-5 kg.
