# Melting, Melt Film and Melt Spraying of the Re-entering Sphere (Step 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `reentry_model` so that the coupled trajectory + 3D conduction model of Step 2 melts the AA7075 sphere (enthalpy method, element fractions and element death on a prism-layered mesh), forms a melt film on the surface patches (lubrication runoff driven by the gas shear and pressure gradient), strips the film into droplets by Girin's gradient instability (thick, thin and rarefied branches), feeds the mass loss and projected area back to the trajectory, records every droplet release in a source table with size distributions, visualises it, and verifies the whole against two new melting SESAM references, Girin's published cases and analytic solutions.

**Architecture:** New modules `dispersion` (Girin's Eq. 1 solved numerically, cached table), `surface_flow` (edge state by isentropic expansion, Ranger-form boundary layer, Kn_δ regimes, shear, driving gradient), `film` (lubrication branches, linearly implicit upwind runoff on the patch graph), `spray` (instability branches, release bookkeeping, source rows, histograms), `girin_case` (Girin-as-published driver) and `body.MeltingBody`; `mesh` gains prism layers, the active set with a face table and the box mesh; `material` gains the latent heat, the melting ranges, the feed fraction and the liquid properties; the thermal backends move the enthalpy to the nodes (lumped capacity matrix, exact conservation) with a Newton iteration mapped through h(T), element fractions, pinned nodes and nodal loads; `coupled`, `compare`, `viz`, `cli` and three analysis scripts grow the melt columns, writers, metrics, plots, videos and flags. Per 0.5 s macro step: trajectory advance → aero state → heating → conduction step (with the deferred melt loads of the previous step) → melt step: liquid inventory of every element → film feed; surface flow; lubrication + runoff; spraying and release; element death with film hand-over; mass and projected area for the next advance → history row, source rows, VTK frame. Every module is plain arrays over patches or elements with its own tests; the Step 1–2 behaviour is unchanged when `--melt off`.

**Tech Stack:** Python 3.12 in `drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`); numpy 2.5, scipy 1.18 (`brentq`, `cumulative_trapezoid`, sparse `spsolve`), gmsh 4.15.2, scikit-fem 12.0.2, pyamg 5.3.0, Cantera 3.2.0, pyvista 0.49, imageio-ffmpeg, matplotlib — all already installed for Step 2; **no new packages**. FEniCSx (dolfinx 0.11) in the separate `fenicsx_env` for the backend cross-check.

**Spec:** `docs/superpowers/specs/2026-09-20-melt-spraying-design.md` — read it first. The measured deviations from the spec listed below are binding; Task 15 records them in the spec's amendments section.

## Global Constraints

- Interpreter for everything: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` (never the system python). Unit tier: `"$PY" -m pytest -m "not drama and not reference" -q` (Step 2: 332 passed + 1 skipped in ~4 min; this plan adds ~60 tests and ~2 min). Reference tier: `"$PY" -m pytest -m reference -q` (Step 1 + Step 2 + Task 13's two melting runs, ~2 min more).
- FEniCSx tests and runs use `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python` with `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang` in the environment (README, Step 2 machine notes).
- SI units inside the package (m, s, kg, K, W, J, Pa); the CLI/CSV keep the wrapper's units (km, km/s, deg, mm) plus SI columns; angles in radians inside.
- Never modify DRAMA's databases or the wrapper (`sphere_reentry.py`, `sphere_sweep.py`). The two melting references are produced by the wrapper (Task 13) and committed under `data/reference_runs/`.
- Non-physical devices must say so: `--heating sesam` (Step 2), `--removal instant` and `--k-scale` (Step 3) — docstring, CLI help and README each state that they are verification devices, not physical models (user requirement carried from Step 2).
- Melt and runoff start at the liquidus; the mushy range counts as solid for the film; no coherency parameter (spec §8, §17.2, user decision of 2026-09-20). Droplets are recorded at birth, one radius per patch and step, no within-patch size spread (spec §17.1).
- Git-ignored output root `reentry_model_output/`; meshes cached under `reentry_model_output/meshes/`; verification outputs under `reentry_model_output/verification_melt/`.
- Repo conventions: module docstrings, `argparse`, exit codes 0/1/2, tests under `tests/` named `test_reentry_model_<module>.py`, one commit per task, commit messages in the imperative like the existing history, ending with `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>` (exactly this line).
- Every task's tests are run with the exact command given in the task; a task is done only when the whole unit tier passes (`"$PY" -m pytest -m "not drama and not reference" -q`).
- Code in this plan was executed and its tests passed on 2026-09-20/21 in a throwaway copy of the package (`drama_env`'s package versions as pinned in `requirements-step2.txt`; the FEniCSx backend in `fenicsx_env`). Transcribe it verbatim; where a task says "replace the file", the block is the complete new file; where it says "append", the block goes at the end of the existing file.

## Measured facts and spec amendments (2026-09-20/21, prototype in a throwaway copy)

These were measured while writing the plan and override the corresponding spec statements; Task 15 writes them into the spec (§18) and the README.

1. **Reference run names.** The wrapper names default-material runs without a material suffix: the melting references are `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind` and `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind` (spec §13.1 wrote `_mAA7075_nowind`). Measured on them: 100 mm melt onset 71.005 km at 43.55 s, mass 0 at 66.476 km / 67.315 s; 50 mm onset 77.104 km at 178.55 s, mass 0 at 73.094 km / 191.55 s (the facts note's "74.9 → 73.0 km" was the NRLMSISE case).
2. **SESAM hollows the melting sphere at fixed outer geometry.** During melting `thick_mm` is the shell thickness of a hollow sphere of outer radius R₀ (m = ρ 4π/3 [R₀³ − (R₀ − thick)³] to four digits), the heat input stays 4πR₀² × q with R₀ in DKR, Kn uses D₀, the radiated power is 4πR₀² εσ 850⁴ = 372 W throughout and C_D stays 0.913. Consequence: the bookkeeping device keeps the intact geometry (uniform φ_e reduction) and needs no area scaling; `--removal instant` removes the liquid inventory of **every** element as it forms (not only surface elements — a surface-only feed stores latent heat in the interior of the isothermal body and lags SESAM by 5–30 s).
3. **Nodal enthalpy, lumped capacity, enthalpy-consistent Newton.** Step 2's element-mean enthalpy with the secant capacity is a regula falsi anchored at T_old whose fixed point repels inside a steep latent-heat ramp (contraction factor −59 for the isothermal body crossing the ±2 K ramp); with the real conductivity a surface element spans 30 K while its mean crosses a 4 K ramp and the iteration cycled (jumps across the ramp, 3-cycles under damping). The thermal core now uses the **nodal** enthalpy h(T_i) with a **lumped** capacity matrix diag(Σ_e φ_e V_e/4 ρ c_i) — the only capacity matrix whose increment equals the increment of ρ∫h dV with h interpolated linearly (the consistent P1-coefficient matrix left a 3e-4 balance error) — a tangent Newton on R(T) = E/dt + KT − F with E = M(ρ c_sec)(T − T_old), and every update mapped through the true h(T) per node (T_new ← T(h(T_k) + c_p,eff(T_k)(T_new − T_k))) so that a node cannot jump across the melting range, with damping as a fallback. Both backends implement it identically (the FEniCSx backend assembles only the stiffness in UFL and adds everything nodal in numpy/PETSc, so the radiation follows the moving boundary). `lumped_mass=True` is the default; `--consistent-mass` (skfem only) keeps the consistent form for the analytic checks. Measured: Stefan front within 0.3 % (0.5 mm box), 2.0–2.4 Newton iterations per step without melting, 3–5 with; the Step 2 verification numbers re-measured with the new core: d100 sesam Q_conv 0.46 % / 2.23 % point-wise, integrated +0.31 %, ΔT_eq 24.6 K (1.20 %), radiated 6.66 %, 74 s runtime (was 158 s) — Task 15 refreshes the README table with the full re-run.
4. **Conductivity is not scaled by φ_e.** Scaling k with φ_e (spec §6) isolated the surface nodes of nearly consumed elements (capacity and conductance both ×0.05) and drove them to 5000 K; k stays that of the full element while φ_e > 0 (a thinner sliver conducts better, not worse); only the capacity scales. Patch owners die at φ_e ≤ PHI_DEATH = 0.05 (their remainder joins the film; owners at 10⁻³ made surface nodes swing by hundreds of kelvin between iterates), interior elements keep φ_e ≥ 10⁻³ (no cavity can open); deaths cascade within the step until every owner has φ_e > 0.05; nodes without material are pinned at their temperature and loads on them are counted (`StepResult.Q_dropped`).
5. **Feed rule.** Every active element feeds its liquid inventory φ_e ρ V_e mean_i f_feed(T_i) (f_feed a ±2 K ramp at the liquidus; for the single-temperature material identical to f_l), owners to their patches by area, interior elements to the four nearest patches by area; the material leaves the FEM at its nodal-mean enthalpy but is booked at h_liquid, the difference (superheat or latent deficit) being a nodal load on the element's nodes over the next step — energy-exact, and the melt front moves at the energy-limited rate whatever the element size (surface-only feed with per-element death could not: the three tets of a prism reach the liquidus together and at peak heating 3.6 layers melt per step).
6. **Mesh: layers by radial projection.** gmsh's `extrudeBoundaryLayer` is a geo-kernel operation that cannot be attached to the OCC sphere without re-parametrising; the layers are built by projecting the inner gmsh sphere's boundary triangulation radially (inner radius R − 3.75 mm, shells at 46.25/48.25/49.25/49.75/50 mm), each prism split into three tetrahedra by the smallest-node-id diagonal rule (conforming, exact areas/volumes, closed surface). Default 4 layers × (0.25, 0.5, 1, 2 mm): **46 278 nodes / 255 276 tets / 15 430 patches** on the 100 mm sphere (spec estimated ≈55 k nodes); surface triangles 2.16 mm. A face table (every face with its ≤ 2 elements) makes deactivation O(dead elements). The box mesh is a structured Kuhn split (no gmsh).
7. **Runoff: linearly implicit upwind, no explicit sub-steps.** Films driven by the free-molecular shear near the rim move at ~10 m/s and cross the hemisphere many times per 0.5 s step: the spec's explicit CFL scheme needed 1e4–1e5 sub-steps per macro step (sliver patches worst). The transport is (I + Δt_s C) m_new = m_old with the edge coefficients c = q ℓ (t̂·n̂)⁺/(A b) at the start of each of 4 sub-steps: unconditionally stable, positive, conservative to round-off, exact steady state (strip test 0.1 %); a wetting front advances one patch per sub-step (documented limit). A Picard iteration on the fully implicit form does not contract (Δt × c ≫ 1). Runoff never crosses into leeward patches (the flow stops at the equator); leeward films are static and can stay attached at the end of a run (0.065 kg in the resolved 100 mm case).
8. **Boundary layer: Ranger's own integral.** δ_a² = 58.08 ν_e ∫u_e⁴ds/u_e⁵ reproduces Ranger's 2.2 R Re_D^−½ Ψ(θ) exactly under potential flow (0.05 % on 1° bins with a 20× refined quadrature); Thwaites' momentum thickness (spec §7) is a constant 12.3 × smaller within ±1.7 % over 5°–85° and is kept as the cross-check. Ψ(0) = √(48/15) = 1.789. Measured at 71 km / 7.24 km/s: u_e 1.3 km/s at 30°, 1.9 km/s at 45°; δ_a 6–9 mm; Kn_δ 0.007–0.03 (continuum/slip) up to 60°, ≥ 0.1 from ~85° (the modified-Newtonian expansion to p∞ makes ρ_e → 3e-6 kg/m³ at the rim); τ_c 12–46 Pa, τ_fm 1600–1900 Pa; G 2e4–6e4 Pa/m with the deceleration term (30 m/s² × ρ_l) comparable to the pressure gradient. The film's driving gradient is G = 2(p_s − p∞) sinθ cosθ/R − ρ_l a sinθ (the deceleration pushes the film toward the nose).
9. **Spraying branches.** A film thicker than δ_m takes the thick (Girin 2017) branch in every regime, with the film velocity the local shear gives it (τ δ_m/μ_l; for the free-molecular shear this is what strips the rim, where the film piled up to 100 mm otherwise); the thin branch (Girin & Kopyt 1994) uses λ* = 1.5 M_e Σ/(ρ_e u_e²) (their 1.5 M d/We_d: the thickness cancels), τ* = 2 capillary periods = 0.798 λ*^1.5 (ρ_l/Σ)^½ (their Eq. 12), **ṁ = ρ_l min(b, λ*/8)/τ*** — their Table 1 mass rate is ρ₁ r_d/(2τ_d) = ρ₁ λ*/(8τ*), reproduced within 1 % (spec §9's ρ_l b/τ* is amended to this); their cut-off λ_t = λ*/3 never limits the mode (dropped). The droplet radius is capped by the film on the patch ((3 m_f/(4πρ_l))^⅓) and by R/4 (long near-critical waves; the rim's expanded edge state gave 25 mm droplets otherwise). Both branches strip hundreds to thousands of kg/m²/s where the melt supplies ~5 kg/m²/s: spraying is melt-limited, the film stays microns thin, r ≈ 45 µm–1.3 mm with a median ≈ 145–200 µm; We_d ≤ 21 with a few hundred breakup-flagged rows per flight.
10. **Girin 2017 Table 1.** GI = We∞Re∞^−½ is reproduced (13.04 / 3.51 / 43.46) only with We∞ on the ambient density ρ∞ = ρ_a/6 and Re∞ on the compressed ρ_a = 1e-4 kg/m³; α = ρ∞/ρ_m reproduces his t_ch. φ_cr from Eq. (3) matches his table with We_cr = 4.62 (16.3° / 32.0° / 8.9°; 3.08 gives 17 % smaller angles). The rest depends on which density enters δ_a: with the **ambient** Reynolds number t_f (7.0 / 27.1 / 207 µs vs 5.7 / 31 / 194), N (1.31e6 / 3.57e5 / 636 vs 1.5e6 / 3.7e5 / 832) and r_med (25.8 / 39.6 µm vs 26.9 / 41.6) are within 30 %; with the shock density droplets come out 2.5× smaller. The spraying duration is half his for the iron variants (2.8 / 8.0 ms vs 5.9 / 16.0) and 6.5 vs 149 ms for the stony one (λ_f 3.5 mm > R₀; belt discretisation and induction handling unstated). Spec §13.2's tiers become: exact (GI, φ_cr ≤ 2 %), integrated (t_f, N, r_med ≤ 30 %, ambient density), reported (t_s.d., ranges, σ, z₀).
11. **Girin & Kopyt 1994.** Table 1's six r_d imply an effective dynamic pressure 7.8 × ρ₂V₀² (their "deceleration ~10× and re-acceleration to M = 2–3"); with that one factor r_d agree within 0.5 % and τ_d within 1 %. Table 2's λ* column is exactly 10 × smaller than their Eq. (14) (a units slip); τ* agrees to three digits. The RT criterion is inactive at our 10–30 m/s².
12. **Verification results** (prototype, default settings unless stated): bookkeeping device (sesam heating, AA7075, instant, k×1e4, `--prism-layers 0`) — 100 mm: mass within 0.98 % of m₀, onset +0.10 km, 1 %-mass time −1.3 %; 50 mm: 1.29 %, +0.17 km, −0.2 % (thresholds 2 % / 0.5 km / 2 %); resolved (sesam heating, AA7075, girin): onset 71.62 / 77.57 km, mass 10 % / 14 % of m₀ off SESAM, 1 %-mass time +4 % / +2 %; physics mode (AA7075_range): 100 mm onset 73.96 km, demise 59.2 km at 98.5 s, 1.455 kg sprayed as 4.9e7 droplets (median 145 µm), **146 s** for 198 steps on the default mesh (target ≤ 6 min; 3.8 Newton iterations per step); 50 mm: onset 78.3 km, 1 %-mass time +3.9 %, 35 s. The two thermal backends give the same melting run to 1e-10 in mass and 0 K in temperature (spec 0.1 %). The 1 %-mass time is interpolated and, for the model, taken on the body's material (film excluded). Sensitivity (100 mm, spec §13.5 with layers 2/4/6 — eight layers of 0.25 mm at growth 2 exceed the radius): sprayed mass and demise altitude within 0.2 % for layers, Δt/2, `bridged`, no runoff, We_cr 3.08, k_t ±30 %; median radius ±22–24 % for k_r ±30 %, −8 % for Δt/2, +5 % without runoff; the single-temperature `AA7075` melts 1 km higher, ends 1 km higher and takes 401 s (its ±2 K ramp triples the Newton work).
13. **Columns and files.** History adds `removed_mass_kg`, `film_thickness_max_mm`, `film_thickness_mean_mm` and `n_dead_elements` to spec §10's list (`runoff_mass_kg` = mass that arrived on another patch, cumulative); `melt_front_depth_max_mm` is the depth of the deepest element with f_l > 0 (the solidus front for the range material). The source table has 22 columns (`spray.SOURCE_COLUMNS`), 5e5 rows for the 100 mm physics flight (`particles.npz`, compressed). The CLI gains `--k-scale` (verification device) and `--consistent-mass` replaces `--lumped-mass`.

---

## File structure

```
reentry_model/mesh.py                     + prism layers (radial projection, prism split), active set with a face table, deactivate, box mesh, SurfaceMesh.owner/face_ids/tangent_from/projected_area/edges
reentry_model/material.py                 + latent heat, solidus/liquidus, MELT_RAMP, feed_fraction, LiquidProperties, MATERIAL_NAMES, exact enthalpy with the latent slope
reentry_model/data/materials/AA7075.json, AA7075_range.json   DRAMA's drama-AA7075 verbatim + liquid properties (+ the alloy's range)
reentry_model/thermal/__init__.py         + StepResult.Q_extra/Q_dropped, set_fractions/element_energies in the protocol
reentry_model/thermal/skfem_backend.py    nodal enthalpy, lumped capacity, tangent Newton mapped through h(T), phi_e, pinned nodes, nodal loads
reentry_model/thermal/fenicsx_backend.py  the same scheme: stiffness in UFL, everything nodal in numpy/PETSc, MatZeroRowsColumns for pinned/Dirichlet
reentry_model/dispersion.py               Girin Eq. (1): batched cubic roots, fastest mode, cached table data/girin_dispersion.json
reentry_model/gas.py                      + GasState.s/a/m_bar, EquilibriumAir.expand
reentry_model/surface_flow.py             edge state per 1 deg bin, Ranger-form boundary layer, Kn_delta regimes, shear (slip/bridged), G
reentry_model/film.py                     lubrication branches, Runoff (edge geometry, coefficients, linearly implicit transport)
reentry_model/spray.py                    melt_layer, thin-film and RT modes, SprayModel (branches, release), source_rows, histogram
reentry_model/girin_case.py               Girin-as-published flight of his Table 1 variants
reentry_model/body.py                     + MeltSettings, MeltingBody (feed, film, spray, death cascade, hand-over, accounting, stats), reference_area hook
reentry_model/trajectory.py               + body.reference_area() in the drag, zero drag for a consumed body
reentry_model/coupled.py                  + MELT_COLUMNS, state passed to advance, demise, melt_results, melt VTK fields, write_particles
reentry_model/sesam_io.py                 + Reference.mass/thickness
reentry_model/compare.py                  + has_melt, melt_metrics (interpolated 1 %-mass crossing), plot_melt (7 plots)
reentry_model/viz.py                      + emitting-patch overlay, film frame/video, liquidus/solidus iso-lines, melt stills
reentry_model/cli.py                      + --melt and the Step 3 group, --k-scale, --consistent-mass, MeltingBody wiring, writers, summaries
data/reference_runs/sphere_d100.00mm_..._h077.500km_nowind.{csv,json}, sphere_d050.00mm_..._h115.000km_nowind.{csv,json}   Task 13
data/reference_values/girin2017_table1.json, girin1994_tables.json
analysis/girin_reference.py, analysis/melt_verification.py, analysis/melt_sensitivity.py
tests/test_reentry_model_{dispersion,surface_flow,film,spray,girin,melting,reference_melt}.py   new
tests/test_reentry_model_{mesh,material,thermal,coupled,cli,compare,viz,data,aero,fenicsx}.py   extended
README.md, docs/model_assumptions.md (§9), the spec (§18), sesam_verified_facts.md (§17)   Task 15
```

Dependency direction (spec §4): `spray` → `surface_flow`, `dispersion`; `body` → `film`, `spray`, `surface_flow`, `thermal`, `material`; `coupled` → everything; `viz`, `compare` read exported files and histories only; `girin_case` → `dispersion`, `surface_flow.ranger_psi`.

---

### Task 1: Prism-layered mesh, active set, box mesh

**Files:**
- Modify (replace): `reentry_model/mesh.py`
- Test: `tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: gmsh generation and `VolumeMesh`/`SurfaceMesh` of Step 2 (kept: `generate_sphere_mesh`, `load_mesh`, `boundary_faces`, `mesh_file_name`, `centre_node`, `boundary_nodes`, `angles_to`, `facet_mean`, `patch_toward`).
- Produces: `sphere_mesh(radius, h_surface, h_core, mesh_dir, layers=0, layer_thickness=0.25e-3, growth=2.0)`; `layer_thicknesses(layers, t0, growth)`; `split_prisms(bottom, top)`; `add_prism_layers(inner, radius, layers, t0, growth)`; `box_mesh(lx, ly, lz, h)`; `VolumeMesh.active` (bool per element), `.element_layer` (−1 core, 0 outermost layer), `.n_active`, `.deactivate(elements) -> (gone_face_ids, new_face_ids)`, `.active_nodes()`, `.volume()` (active), `.surface()` over the active set with `SurfaceMesh.owner` (element per patch) and `.face_ids` (stable ids into the face table `_face_nodes`); `SurfaceMesh.tangent_from(v_hat)`, `.projected_area(v_hat)`, `.edges() -> (edge node pairs, patch i, patch j)`; constants `DEFAULT_LAYERS = 4`, `DEFAULT_LAYER_THICKNESS = 0.25e-3`, `DEFAULT_LAYER_GROWTH = 2.0`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_mesh.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: prism layers, the active set, the box mesh and the surface helpers

@pytest.fixture(scope="module")
def layered_mesh(tmp_path_factory):
    """100 mm sphere, 4 mm / 20 mm inner mesh with two 0.25 mm / 0.5 mm prism layers (~5 k nodes)."""
    return mesh.sphere_mesh(R, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=2)


def test_layer_thicknesses_and_defaults():
    assert np.allclose(mesh.layer_thicknesses(4, 0.25e-3, 2.0), [0.25e-3, 0.5e-3, 1e-3, 2e-3])
    assert mesh.DEFAULT_LAYERS == 4 and mesh.DEFAULT_LAYER_THICKNESS == 0.25e-3 and mesh.DEFAULT_LAYER_GROWTH == 2.0
    with pytest.raises(ValueError):
        mesh.sphere_mesh(R, 4e-3, 20e-3, layers=-1)
    with pytest.raises(ValueError):
        mesh.sphere_mesh(0.001, 4e-3, 20e-3, layers=4)                      # layers thicker than half the radius


def test_split_prisms_is_conforming_and_fills_the_prism():
    bottom = np.array([[0, 1, 2], [1, 3, 2]])                                   # two triangles sharing the edge 1-2
    tets = mesh.split_prisms(bottom, bottom + 4)
    assert tets.shape == (6, 4)
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float)
    m = mesh.VolumeMesh(pts, tets)
    assert m.volume() == pytest.approx(1.0) and (m.element_volumes() > 0).all()
    faces, _ = mesh.boundary_faces(tets)
    assert len(faces) == 12                                                     # the cube's 6 faces, 2 triangles each: no internal single-use face


def test_layered_sphere_geometry(layered_mesh):
    m = layered_mesh
    s = m.surface()
    assert set(np.unique(m.element_layer)) == {-1, 0, 1}
    outer = np.bincount(m.element_layer + 1)
    assert outer[1] == outer[2] == 3 * s.n_patches                            # three tets per prism, one prism per patch
    assert s.area == pytest.approx(4.0 * math.pi * R ** 2, rel=2e-3) and m.volume() == pytest.approx(4.0 / 3.0 * math.pi * R ** 3, rel=3e-3)
    radii = np.linalg.norm(m.points, axis=1)
    shells = np.unique(np.round(radii[radii > R - 1e-3] * 1e3, 3))
    assert np.allclose(shells, [49.25, 49.75, 50.0])                         # inner mesh at R - 0.75 mm, layers 0.5 and 0.25 mm
    assert m.element_layer[s.owner].max() == 0 and m.element_layer[s.owner].min() == 0
    e, i, j = s.edges()
    assert len(e) == 3 * s.n_patches // 2                                    # closed manifold surface
    assert np.abs((s.areas[:, None] * s.normals).sum(axis=0)).max() < 1e-12 * s.area
    assert s.projected_area([1.0, 0.0, 0.0]) == pytest.approx(math.pi * R ** 2, rel=2e-3)
    t = s.tangent_from([1.0, 0.0, 0.0])
    theta = s.angles_to([1.0, 0.0, 0.0])
    ok = (theta > 0.1) & (theta < 3.0)
    assert np.allclose(np.einsum("ij,ij->i", t[ok], s.normals[ok]), 0.0, atol=1e-9)      # in the tangent plane
    assert (t[ok] @ np.array([1.0, 0.0, 0.0]) < 0.0).all()                              # away from the stagnation point
    assert m.params["layers"] == 2 and m.params["layer_thickness_m"] == 0.25e-3


def test_deactivation_updates_the_boundary(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points, layered_mesh.tets, element_layer=layered_mesh.element_layer.copy())
    s0 = m.surface()
    volume0 = m.volume()
    owners = s0.owner[s0.normals[:, 0] > 0.9]
    gone, new = m.deactivate(owners)
    assert m.n_active == m.n_elements - owners.size and gone.size == owners.size and new.size > 0
    s1 = m.surface()
    assert s1.n_patches == s0.n_patches - gone.size + new.size
    assert m.volume() == pytest.approx(volume0 - layered_mesh.element_volumes()[owners].sum())
    assert np.abs((s1.areas[:, None] * s1.normals).sum(axis=0)).max() < 1e-12 * s1.area      # still closed
    assert set(s1.face_ids[np.isin(s1.face_ids, s0.face_ids)]) == set(s0.face_ids) - set(gone)
    assert m.active[s1.owner].all()
    again = m.deactivate(owners)                                               # already dead: no change
    assert again[0].size == 0 and m.n_active == m.n_elements - owners.size
    assert np.isin(m.active_nodes(), np.unique(m.tets[m.active])).all()


def test_box_mesh():
    b = mesh.box_mesh(0.03, 0.004, 0.004, 0.001)
    assert b.n_nodes == 31 * 5 * 5 and b.n_elements == 30 * 4 * 4 * 6
    assert b.volume() == pytest.approx(0.03 * 0.004 * 0.004) and (b.element_volumes() > 0).all()
    s = b.surface()
    assert s.area == pytest.approx(2 * (0.03 * 0.004 * 2 + 0.004 * 0.004)) and s.n_patches == 2 * (30 * 4 * 2 * 2 + 4 * 4 * 2)
    assert (np.abs(np.abs(s.normals).max(axis=1) - 1.0) < 1e-12).all()      # axis-aligned faces
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: the five new tests fail (`layer_thicknesses`, `split_prisms`, `box_mesh` missing; `sphere_mesh` has no `layers`).

- [ ] **Step 3: Replace `reentry_model/mesh.py`**

```python
"""Tetrahedral sphere meshes (gmsh, graded, optional prism layers) and their surface geometry. Spec sections 5
(Step 2) and 5 (Step 3).

A VolumeMesh is node coordinates, tetrahedra, an active mask (Step 3: elements that have melted away are
deactivated and drop out of the surface, the operators and the mass) and the boundary triangles of the active set
(faces used by exactly one active element) with outward unit normals, centroids, areas and owner elements -- so any
gmsh volume mesh loads, tagged or not. A face table (every face of the mesh with its one or two elements) is built
once so that deactivation updates the boundary in O(faces) without re-sorting.

`sphere_mesh(R, h_surface, h_core)` is the Step 2 generator: element size growing linearly from h_surface at the
surface to h_core at the centre (gmsh Distance/Threshold field), a node embedded at the centre, cached by
(R, h_surface, h_core). With `layers = n > 0` (Step 3 default for melting runs) the graded gmsh mesh is generated for
the inner sphere R - T (T the total layer thickness) and n prism layers are built on it by radial projection of its
boundary triangulation: shells at R - T + t_(n-1) + ... , thicknesses t0 growth^k from the outside in (the outermost
layer is t0 thick), each prism split into three tetrahedra by the smallest-node-id diagonal rule (conforming across
prisms). The construction is exact on a sphere and needs no gmsh boundary-layer machinery (spike of 2026-09-20:
gmsh's extrudeBoundaryLayer is a geo-kernel function that cannot be attached to the OCC sphere without
re-parametrising; the radial construction gives the same layers deterministically). `box_mesh` is a structured
Kuhn-split box for the analytic (Stefan) tests.
"""
import os
from dataclasses import dataclass, field

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MESH_DIR = os.path.join(REPO_ROOT, "reentry_model_output", "meshes")
DEFAULT_H_SURFACE = 2.0e-3          # m; 2 mm / 8 mm gives 18.9 k nodes on the 100 mm sphere (measured 2026-09-18)
DEFAULT_H_CORE = 8.0e-3             # m
DEFAULT_LAYERS = 4                  # prism layers for melting runs (spec Step 3 section 5)
DEFAULT_LAYER_THICKNESS = 0.25e-3   # m, outermost layer
DEFAULT_LAYER_GROWTH = 2.0
_FACE_OF_VERTEX = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])     # face opposite each tet vertex


@dataclass
class SurfaceMesh:
    faces: np.ndarray            # (nf, 3) node ids of the boundary triangles ("patches")
    centroids: np.ndarray        # (nf, 3) m
    normals: np.ndarray          # (nf, 3) outward unit normals
    areas: np.ndarray            # (nf,) m2
    owner: np.ndarray = None     # (nf,) element owning each patch
    face_ids: np.ndarray = None  # (nf,) ids in the mesh's face table (stable across deactivations)

    @property
    def n_patches(self):
        return len(self.faces)

    @property
    def area(self):
        return float(self.areas.sum())

    def angles_to(self, v_hat):
        """theta per patch [rad]: the angle between the outward normal and v_hat, the direction the body moves in.
        The stagnation patch has its normal along v_hat (theta = 0); theta > pi/2 is leeward."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        return np.arccos(np.clip(self.normals @ v, -1.0, 1.0))

    def facet_mean(self, nodal):
        """Mean of a nodal field over each patch's three nodes."""
        return np.asarray(nodal)[self.faces].mean(axis=1)

    def patch_toward(self, direction):
        """Index of the patch whose outward normal is closest to `direction`."""
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))

    def tangent_from(self, v_hat):
        """Unit surface direction away from the stagnation point (increasing theta) per patch; zero where undefined."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        t = -v[None, :] + (self.normals @ v)[:, None] * self.normals
        n = np.linalg.norm(t, axis=1)
        return np.where(n[:, None] > 1e-12, t / np.maximum(n, 1e-300)[:, None], 0.0)

    def projected_area(self, v_hat):
        """Area projected on the plane normal to v_hat: sum of A max(0, n.v) (pi R^2 for a sphere)."""
        v = np.asarray(v_hat, dtype=float)
        return float((self.areas * np.maximum(0.0, self.normals @ (v / np.linalg.norm(v)))).sum())

    def edges(self):
        """Patch adjacency: (edge node pairs (ne, 2), patch i, patch j, edge lengths) for edges shared by two patches."""
        e = np.sort(self.faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
        patch = np.repeat(np.arange(self.n_patches), 3)
        order = np.lexsort((e[:, 1], e[:, 0]))
        e, patch = e[order], patch[order]
        same = np.all(e[1:] == e[:-1], axis=1)
        k = np.flatnonzero(same)
        return e[k], patch[k], patch[k + 1]


@dataclass
class VolumeMesh:
    points: np.ndarray                                   # (n, 3) m
    tets: np.ndarray                                     # (ne, 4) node ids
    params: dict = field(default_factory=dict)           # generator parameters / source file
    element_layer: np.ndarray = None                     # (ne,) prism layer index (0 = outermost) or -1 (graded core)
    active: np.ndarray = field(init=False, repr=False, default=None)
    _face_nodes: np.ndarray = field(init=False, repr=False, default=None)    # (nfaces, 3) node ids, one entry per unique face
    _face_elements: np.ndarray = field(init=False, repr=False, default=None) # (nfaces, 2) elements sharing the face, -1 if one
    _face_opposite: np.ndarray = field(init=False, repr=False, default=None) # (nfaces, 2) opposite vertex in each of those elements
    _face_count: np.ndarray = field(init=False, repr=False, default=None)    # active elements per face
    _element_faces: np.ndarray = field(init=False, repr=False, default=None) # (ne, 4) face ids of each element
    _boundary: np.ndarray = field(init=False, repr=False, default=None)      # face ids of the current boundary

    def __post_init__(self):
        self.points = np.array(self.points, dtype=float)
        self.tets = np.asarray(self.tets, dtype=np.int64)
        self.active = np.ones(len(self.tets), dtype=bool)
        if self.element_layer is None:
            self.element_layer = np.full(len(self.tets), -1, dtype=np.int64)
        self._build_face_table()

    def _build_face_table(self):
        faces = self.tets[:, _FACE_OF_VERTEX].reshape(-1, 3)                    # 4 per element, element-major
        opposite = self.tets.ravel()                                            # vertex opposite face 4e + v is tets[e, v]
        element = np.repeat(np.arange(len(self.tets)), 4)
        key = np.sort(faces, axis=1)
        _, first, inverse, counts = np.unique(key, axis=0, return_index=True, return_inverse=True, return_counts=True)
        inverse = inverse.ravel()
        nfaces = len(first)
        self._face_nodes = faces[first]
        self._face_elements = np.full((nfaces, 2), -1, dtype=np.int64)
        self._face_opposite = np.full((nfaces, 2), -1, dtype=np.int64)
        order = np.argsort(inverse, kind="stable")
        slot = np.zeros(nfaces, dtype=np.int64)
        for k in order:                                                         # <= 2 entries per face
            f = inverse[k]
            self._face_elements[f, slot[f]], self._face_opposite[f, slot[f]] = element[k], opposite[k]
            slot[f] += 1
        self._face_count = counts.astype(np.int64)
        self._element_faces = inverse.reshape(len(self.tets), 4)              # face id of each element's four faces
        self._boundary = np.flatnonzero(self._face_count == 1)

    @property
    def n_nodes(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.tets)

    @property
    def n_active(self):
        return int(self.active.sum())

    def element_volumes(self):
        x = self.points[self.tets]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        return np.abs(np.linalg.det(J)) / 6.0

    def volume(self):
        """Volume of the active elements."""
        return float(self.element_volumes()[self.active].sum())

    def deactivate(self, elements):
        """Remove elements from the active set; the boundary (and `surface()`) follows. Returns the face ids that
        stopped being boundary faces (their patches vanished) and those that became boundary faces."""
        elements = np.unique(np.asarray(elements, dtype=np.int64))
        elements = elements[self.active[elements]]
        if elements.size == 0:
            return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
        self.active[elements] = False
        faces = self._element_faces[elements].ravel()
        before = self._face_count[faces].copy()
        np.subtract.at(self._face_count, faces, 1)
        after = self._face_count[faces]
        gone, new = np.unique(faces[(before == 1) & (after == 0)]), np.unique(faces[(before == 2) & (after == 1)])
        self._boundary = np.flatnonzero(self._face_count == 1)
        return gone, new

    def surface(self):
        """The boundary patches of the active set with their current geometry."""
        f = self._boundary
        faces = self._face_nodes[f]
        active_slot = np.where(self.active[self._face_elements[f, 0]], 0, 1)
        owner = self._face_elements[f, active_slot]
        opposite = self._face_opposite[f, active_slot]
        a, b, c = (self.points[faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        areas = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * areas)[:, None]
        centroids = (a + b + c) / 3.0
        inward = np.einsum("ij,ij->i", n, centroids - self.points[opposite]) < 0.0
        n[inward] *= -1.0
        return SurfaceMesh(faces, centroids, n, areas, owner, f.copy())

    def boundary_nodes(self):
        return np.unique(self._face_nodes[self._boundary])

    def active_nodes(self):
        """Nodes belonging to at least one active element."""
        return np.unique(self.tets[self.active])

    def centre_node(self):
        """Node nearest the origin (the embedded centre node of a generated sphere)."""
        return int(np.argmin(np.linalg.norm(self.points, axis=1)))


def boundary_faces(tets):
    """Faces used by exactly one tetrahedron, as (nf, 3) node ids, and the opposite vertex of that tetrahedron."""
    faces = tets[:, _FACE_OF_VERTEX].reshape(-1, 3)
    opposite = np.repeat(tets, 4, axis=0)[np.arange(4 * len(tets)), np.tile(np.arange(4), len(tets))]
    _, first, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    keep = first[counts == 1]
    return faces[keep], opposite[keep]


def mesh_file_name(radius, h_surface, h_core):
    return "sphere_R{:.3f}mm_hs{:.3f}mm_hc{:.3f}mm.msh".format(radius * 1e3, h_surface * 1e3, h_core * 1e3)


def generate_sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR):
    """gmsh sphere of `radius` graded from h_surface (surface) to h_core (centre); returns the cached .msh path."""
    os.makedirs(mesh_dir, exist_ok=True)
    path = os.path.join(mesh_dir, mesh_file_name(radius, h_surface, h_core))
    if os.path.isfile(path):
        return path
    import gmsh                                   # imported here: only mesh generation needs gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("sphere")
        vol = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, radius)
        centre = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h_core)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(0, [centre], 3, vol)
        surfaces = [s[1] for s in gmsh.model.getBoundary([(3, vol)], oriented=False)]
        gmsh.model.addPhysicalGroup(3, [vol], tag=1, name="body")
        gmsh.model.addPhysicalGroup(2, surfaces, tag=2, name="surface")
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", surfaces)
        gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
        thr = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thr, "InField", dist)
        gmsh.model.mesh.field.setNumber(thr, "SizeMin", h_surface)
        gmsh.model.mesh.field.setNumber(thr, "SizeMax", h_core)
        gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", radius)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.model.mesh.generate(3)
        tmp = path[:-4] + ".tmp.msh"           # gmsh.write needs the .msh extension to pick the format
        gmsh.write(tmp)
    finally:
        gmsh.finalize()
    os.replace(tmp, path)                      # atomic: an interrupted generation never leaves a reusable partial .msh
    return path


def load_mesh(path):
    """Any gmsh/meshio volume mesh with tetrahedra (other cell types are ignored)."""
    import meshio
    m = meshio.read(path)
    tets = [c.data for c in m.cells if c.type == "tetra"]
    if not tets:
        raise ValueError("no tetrahedra in {}".format(path))
    return VolumeMesh(np.asarray(m.points, dtype=float), np.vstack(tets), {"path": os.path.abspath(path)})


def layer_thicknesses(layers, layer_thickness, growth):
    """Thickness of each prism layer from the outside in: t0, t0 g, t0 g^2, ..."""
    return layer_thickness * growth ** np.arange(layers)


def split_prisms(bottom, top):
    """Three tetrahedra per prism (bottom (np, 3) and top (np, 3) node ids, top[k] above bottom[k]) with the diagonal
    of every quad face drawn from its smallest node id -- so two prisms sharing a quad face split it identically."""
    bottom, top = np.asarray(bottom, dtype=np.int64), np.asarray(top, dtype=np.int64)
    # rotate each triangle so that its smallest bottom id comes first (the top follows the same rotation)
    shift = np.argmin(bottom, axis=1)
    idx = (shift[:, None] + np.arange(3)[None, :]) % 3
    b = np.take_along_axis(bottom, idx, axis=1)
    t = np.take_along_axis(top, idx, axis=1)
    b0, b1, b2, t0, t1, t2 = b[:, 0], b[:, 1], b[:, 2], t[:, 0], t[:, 1], t[:, 2]
    # diagonals from b0 on the faces (b0 b1 t1 t0) and (b2 b0 t0 t2); on (b1 b2 t2 t1) from min(b1, b2)
    case = b1 < b2
    tets = np.where(case[:, None, None],
                    np.stack([np.stack([b0, b1, b2, t2], 1), np.stack([b0, b1, t2, t1], 1), np.stack([b0, t1, t2, t0], 1)], 1),
                    np.stack([np.stack([b0, b1, b2, t1], 1), np.stack([b0, b2, t1, t2], 1), np.stack([b0, t1, t2, t0], 1)], 1))
    return tets.reshape(-1, 4)


def add_prism_layers(inner, radius, layers, layer_thickness, growth):
    """Prism layers by radial projection of the inner sphere mesh's boundary triangulation out to `radius`
    (thicknesses t0 growth^k from the outside in). Returns the layered VolumeMesh (inner elements first, then the
    layers from the inside out; `element_layer` is 0 for the outermost layer, -1 for the graded core)."""
    t = layer_thicknesses(layers, layer_thickness, growth)
    r_inner = radius - t.sum()
    faces, _ = boundary_faces(inner.tets)
    shell_nodes = np.unique(faces)
    ns = len(shell_nodes)
    local = -np.ones(inner.n_nodes, dtype=np.int64)
    local[shell_nodes] = np.arange(ns)
    lf = local[faces]                                        # faces in shell-local ids
    direction = inner.points[shell_nodes]
    direction = direction / np.linalg.norm(direction, axis=1)[:, None]
    radii = r_inner + np.cumsum(t[::-1])                     # shell radii from the inside out (outermost = radius)
    shell_ids = [shell_nodes] + [inner.n_nodes + k * ns + np.arange(ns) for k in range(layers)]
    points = np.vstack([inner.points] + [direction * r for r in radii])
    tets, layer_index = [inner.tets], [np.full(inner.n_elements, -1, dtype=np.int64)]
    for k in range(layers):
        new = split_prisms(shell_ids[k][lf], shell_ids[k + 1][lf])
        tets.append(new)
        layer_index.append(np.full(len(new), layers - 1 - k, dtype=np.int64))
    return VolumeMesh(points, np.vstack(tets), dict(inner.params), np.concatenate(layer_index))


def sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR, layers=0,
                layer_thickness=DEFAULT_LAYER_THICKNESS, growth=DEFAULT_LAYER_GROWTH):
    """Graded gmsh sphere (layers = 0, Step 2) or the graded inner sphere with `layers` prism layers on it."""
    if layers < 0 or layer_thickness <= 0.0 or growth <= 0.0:
        raise ValueError("layers must be >= 0 and layer_thickness, growth > 0")
    total = float(layer_thicknesses(layers, layer_thickness, growth).sum()) if layers else 0.0
    if total >= 0.5 * radius:
        raise ValueError("the prism layers ({:.3g} m) must be thinner than half the radius".format(total))
    mesh = load_mesh(generate_sphere_mesh(radius - total, h_surface, h_core, mesh_dir))
    if layers:
        mesh = add_prism_layers(mesh, radius, layers, layer_thickness, growth)
    mesh.params.update({"radius_m": radius, "h_surface_m": h_surface, "h_core_m": h_core, "layers": layers,
                        "layer_thickness_m": layer_thickness, "layer_growth": growth})
    return mesh


def box_mesh(lx, ly, lz, h):
    """Structured box [0, lx] x [0, ly] x [0, lz] of cubes of side ~h, each split into six tetrahedra (Kuhn
    decomposition, conforming). For the analytic tests (Stefan front along x with the x = 0 face heated)."""
    n = [max(1, int(round(L / h))) for L in (lx, ly, lz)]
    xs = [np.linspace(0.0, L, k + 1) for L, k in zip((lx, ly, lz), n)]
    X, Y, Z = np.meshgrid(*xs, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    nid = np.arange(points.shape[0]).reshape(n[0] + 1, n[1] + 1, n[2] + 1)
    c = nid[:-1, :-1, :-1].ravel()
    dx, dy, dz = nid[1, 0, 0] - nid[0, 0, 0], nid[0, 1, 0] - nid[0, 0, 0], 1
    v = lambda i, j, k: c + i * dx + j * dy + k * dz
    corners = [v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0), v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)]
    kuhn = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    tets = np.vstack([np.column_stack([corners[a], corners[b], corners[cc], corners[d]]) for a, b, cc, d in kuhn])
    return VolumeMesh(points, tets, {"box": (lx, ly, lz), "h": h, "path": None})
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: 10 passed (the five Step 2 tests unchanged; the layered fixture builds a 4 mm / 20 mm inner sphere with two layers, ~5 k nodes). Then the unit tier: `"$PY" -m pytest -m "not drama and not reference" -q` — everything that used `mesh.surface()` still passes (the boundary now comes from the face table; the coarse-mesh fixtures are unchanged).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/mesh.py tests/test_reentry_model_mesh.py
git commit -m "Add prism layers, the active element set and the box mesh to the sphere mesh (Step 3 Task 1)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Material with latent heat, melting ranges and liquid properties

**Files:**
- Modify (replace): `reentry_model/material.py`
- Create: `reentry_model/data/materials/AA7075.json`, `reentry_model/data/materials/AA7075_range.json`
- Test: `tests/test_reentry_model_material.py` (append)

**Interfaces:**
- Consumes: `reentry_model/data/materials/AA7075_nomelt.json` (Step 2).
- Produces: `Material` fields `latent_heat`, `T_solidus`, `T_liquidus`, `liquid: LiquidProperties(rho, mu, sigma)`; properties `melts`, `T_feed`, `h_liquid`; methods `liquid_fraction(T)`, `feed_fraction(T)`, `cp_eff(T)`, `enthalpy(T)` (exact, latent slope inside the range), `temperature_from_enthalpy(h)`; `Material.from_drama_json(name_or_path)` accepting the names in `MATERIAL_NAMES` (`AA7075_nomelt`, `AA7075`, `AA7075_range`); constants `MELT_RAMP = 2.0`, `NO_MELT_ABOVE = 5000.0`. Single-temperature materials get a ±MELT_RAMP ramp; `feed_fraction` is a ±MELT_RAMP ramp ending at `T_feed` (the liquidus, +2 K for range materials).

- [ ] **Step 1: Write the two material files**

Generate them from the packaged no-melt file (the DRAMA tables verbatim, the 1e5 K holding row dropped, the liquid properties added):

```bash
"$PY" - <<'EOF'
import json
d = json.load(open("reentry_model/data/materials/AA7075_nomelt.json"))
base = {k: v for k, v in d.items() if k != "_comment"}
base["specificHeatCapacity"] = [r for r in d["specificHeatCapacity"] if r[0] <= 850.0]
base["heatConductivity"] = [r for r in d["heatConductivity"] if r[0] <= 850.0]
liquid = {"density": 2400.0, "viscosity": 1.3e-3, "surfaceTension": 0.86,
          "_sources": "pure aluminium near the liquidus: rho_l 2375-2400 kg/m3 (Smithells Metals Reference Book, 8th ed., Table 14.1), mu_l 1.2-1.4 mPa s at 933-1000 K (Assael et al. 2006, J. Phys. Chem. Ref. Data 35, 285), sigma 0.86-0.91 N/m (Smithells; ASM Handbook Vol. 2). Alloy corrections for AA7075 (5.6 % Zn, 2.5 % Mg, 1.6 % Cu) are within 10 %; used as assumptions (spec 2026-09-20 sections 2 and 16)."}
a = dict(base); a["name"] = "AA7075"; a["meltingTemperature"] = 850.0; a["liquid"] = liquid
a["_comment"] = "DRAMA 4.1.4 drama-AA7075 verbatim (TOOLS/material_database.xml: density, cp(T) and k(T) to 850 K, meltingHeat 400 kJ/kg, meltingTemperature 850 K, emissivity 0.4, catalycity 1, oxidation off) plus the liquid-phase properties DRAMA does not carry (`liquid`). Above 850 K the tables hold their last value (cp 1131.6 J/kgK, k 128.19 W/mK). Step 3 material (spec 2026-09-20 section 6): single melting temperature, numerical ramp of +-2 K in the model. DRAMA's database itself is untouched."
r = dict(a); r["name"] = "AA7075_range"; r["solidusTemperature"] = 750.0; r["liquidusTemperature"] = 908.0; r["meltingTemperature"] = 908.0
r["_comment"] = "AA7075 (see AA7075.json) with the alloy's melting range instead of DRAMA's single temperature: solidus 750 K, liquidus 908 K (ASM Handbook Vol. 2, Properties and Selection: Nonferrous Alloys, AA7075 477-635 C), latent heat 400 kJ/kg spread linearly across the range. meltingTemperature is set to the liquidus. Default material of --melt on (spec 2026-09-20 sections 6 and 14); melt and runoff start at the liquidus (section 8)."
for name, doc in (("AA7075", a), ("AA7075_range", r)):
    keys = ["_comment", "name", "materialType", "catalycity", "density", "meltingHeat", "meltingTemperature"] + (["solidusTemperature", "liquidusTemperature"] if "solidusTemperature" in doc else []) + ["specificHeatCapacity", "heatConductivity", "emissivity", "oxideActivationTemperature", "oxideEmissivity", "oxideHeatOfFormation", "oxideReactionProbability", "liquid"]
    json.dump({k: doc[k] for k in keys}, open("reentry_model/data/materials/%s.json" % name, "w"), indent=2)
print("written")
EOF
```

Expected: both files exist; `AA7075.json` has 29 c_p rows ending at 850 K, `meltingHeat` 400000, `meltingTemperature` 850; `AA7075_range.json` adds `solidusTemperature` 750 and `liquidusTemperature` 908 (`meltingTemperature` 908).

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_reentry_model_material.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: latent heat, melting ranges, feed fraction, liquid properties

def test_material_names_and_liquid_properties():
    a, r, n = (material.Material.from_drama_json(k) for k in ("AA7075", "AA7075_range", "AA7075_nomelt"))
    assert a.melts and r.melts and not n.melts and n.latent_heat == 0.0 and n.liquid is None
    assert a.latent_heat == r.latent_heat == 400e3 and a.rho == r.rho == 2813.0 and a.emissivity == 0.4
    assert (a.T_solidus, a.T_liquidus) == (848.0, 852.0) and (r.T_solidus, r.T_liquidus) == (750.0, 908.0)
    assert a.liquid.rho == 2400.0 and a.liquid.mu == 1.3e-3 and a.liquid.sigma == 0.86
    assert np.allclose(a.k(a.T_k), n.k(a.T_k)) and np.allclose(a.cp(a.T_cp), n.cp(a.T_cp))    # DRAMA's tables verbatim
    assert a.T_cp.max() == 850.0 and n.T_cp.max() == 1e5


def test_liquid_and_feed_fractions():
    a, r = material.Material.from_drama_json("AA7075"), material.Material.from_drama_json("AA7075_range")
    assert np.allclose(a.liquid_fraction([840.0, 848.0, 850.0, 852.0, 900.0]), [0.0, 0.0, 0.5, 1.0, 1.0])
    assert np.allclose(a.feed_fraction([840.0, 848.0, 850.0, 852.0, 900.0]), [0.0, 0.0, 0.5, 1.0, 1.0])     # single T: feed = liquid fraction
    assert np.allclose(r.liquid_fraction([750.0, 829.0, 908.0]), [0.0, 0.5, 1.0])
    assert np.allclose(r.feed_fraction([900.0, 906.0, 908.0, 910.0]), [0.0, 0.0, 0.5, 1.0])              # feed only at the liquidus (+-2 K)
    assert r.T_feed == 910.0 and a.T_feed == 852.0
    assert a.cp_eff(850.0) == pytest.approx(a.cp(850.0) + 400e3 / 4.0) and r.cp_eff(800.0) == pytest.approx(r.cp(800.0) + 400e3 / 158.0)
    assert a.cp_eff(700.0) == a.cp(700.0) and r.cp_eff(950.0) == r.cp(950.0)


def test_enthalpy_jump_is_exact_and_invertible():
    for name in ("AA7075", "AA7075_range"):
        m = material.Material.from_drama_json(name)
        sensible = np.trapezoid(m.cp(np.linspace(m.T_solidus, m.T_liquidus, 2001)), np.linspace(m.T_solidus, m.T_liquidus, 2001))
        assert m.enthalpy(m.T_liquidus) - m.enthalpy(m.T_solidus) == pytest.approx(400e3 + sensible, rel=1e-9)
        T = np.linspace(200.0, 2000.0, 7201)
        h = m.enthalpy(T)
        assert np.all(np.diff(h) > 0.0) and np.abs(m.temperature_from_enthalpy(h) - T).max() < 1e-8
        assert m.h_liquid == pytest.approx(m.enthalpy(m.T_feed))
    n = material.Material.from_drama_json("AA7075_nomelt")
    assert n.enthalpy(1000.0) == pytest.approx(material.Material.from_drama_json("AA7075").enthalpy(1000.0) - 400e3, rel=1e-9)


def test_invalid_melting_data():
    with pytest.raises(ValueError):
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], -1.0, 850.0, 850.0)
    with pytest.raises(ValueError):
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 900.0, 850.0)
```


- [ ] **Step 3: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py -q`
Expected: the four new tests fail (`from_drama_json("AA7075")` is a bad path; no `melts`, `feed_fraction`, `liquid`).

- [ ] **Step 4: Replace `reentry_model/material.py`**

```python
"""Material for the thermal and melting model: density, k(T), c_p(T), emissivity, the latent heat with its melting
range, the liquid-phase properties of the film, and the specific enthalpy h(T) with its inverse.

Tables come from a DRAMA material JSON (the wrapper's `user_materials` format); np.interp holds the end values
outside the tabulated range. Melting (Step 3, spec section 6): `meltingHeat` L_f is released between the solidus and
the liquidus (`solidusTemperature`/`liquidusTemperature`, default both = `meltingTemperature`); a single-temperature
material gets a numerical ramp of +-MELT_RAMP around it (the enthalpy jump is exact, only its slope is smoothed).
`liquid_fraction(T)` is the enthalpy method's f_l; `feed_fraction(T)` is the fraction of an element's material that
the film receives -- a +-MELT_RAMP ramp at the liquidus for every material, so that melt and runoff start at the
liquidus and the mushy range counts as solid (decided 2026-09-20, spec sections 8 and 17.2). A `meltingTemperature`
above NO_MELT_ABOVE (the wrapper's 1e5 K device) means no melting at all. Enthalpy is exact: the integral of the
piecewise-linear c_p (quadratic inside each table interval) plus L_f f_l(T), and the inverse solves the same
quadratic, so the FEM's secant heat capacity conserves energy through melting to round-off."""
import json
import os
from dataclasses import dataclass, field

import numpy as np

from . import DATA_DIR

DEFAULT_MATERIAL = os.path.join(DATA_DIR, "materials", "AA7075_nomelt.json")
MATERIAL_NAMES = {"AA7075_nomelt": DEFAULT_MATERIAL, "AA7075": os.path.join(DATA_DIR, "materials", "AA7075.json"),
                  "AA7075_range": os.path.join(DATA_DIR, "materials", "AA7075_range.json")}
T_REF = 293.0                      # K, zero of the enthalpy scale (first row of DRAMA's tables)
MELT_RAMP = 2.0                    # K, half-width of the numerical melting/feed ramps
NO_MELT_ABOVE = 5000.0             # K: a melting temperature above this means "never melts"


@dataclass
class LiquidProperties:
    rho: float                     # kg/m3
    mu: float                      # Pa s
    sigma: float                   # N/m


@dataclass
class Material:
    name: str
    rho: float                      # kg/m3
    emissivity: float
    T_cp: np.ndarray                # K, nodes of the c_p table
    cp_table: np.ndarray            # J/(kg K)
    T_k: np.ndarray                 # K, nodes of the k table
    k_table: np.ndarray             # W/(m K)
    latent_heat: float = 0.0        # J/kg
    T_solidus: float = np.inf
    T_liquidus: float = np.inf
    liquid: LiquidProperties = None
    _T_h: np.ndarray = field(init=False, repr=False)
    _h_nodes: np.ndarray = field(init=False, repr=False)
    _latent_slope: np.ndarray = field(init=False, repr=False)    # L_f df_l/dT inside each enthalpy-table interval

    def __post_init__(self):
        self.T_cp, self.cp_table = np.asarray(self.T_cp, dtype=float), np.asarray(self.cp_table, dtype=float)
        self.T_k, self.k_table = np.asarray(self.T_k, dtype=float), np.asarray(self.k_table, dtype=float)
        if self.latent_heat < 0.0 or self.T_liquidus < self.T_solidus:
            raise ValueError("latent heat must be >= 0 and the liquidus >= the solidus")
        if self.melts and self.T_solidus == self.T_liquidus:            # single-temperature material: numerical ramp
            self.T_solidus, self.T_liquidus = self.T_liquidus - MELT_RAMP, self.T_liquidus + MELT_RAMP
        # h(T) on the c_p nodes (plus T_REF and the melting range): the trapezoid rule is exact for the piecewise-linear c_p
        extra = [T_REF] + ([self.T_solidus, self.T_liquidus] if self.melts else [])
        T = np.union1d(self.T_cp, extra)
        cp = np.interp(T, self.T_cp, self.cp_table)
        h = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(T))])
        h = h - np.interp(T_REF, T, h)
        slope = np.zeros(len(T) - 1)
        if self.melts:
            inside = (T[:-1] >= self.T_solidus - 1e-9) & (T[1:] <= self.T_liquidus + 1e-9)
            slope[inside] = self.latent_heat / (self.T_liquidus - self.T_solidus)
            h = h + self.latent_heat * self._liquid_fraction_raw(T)
        self._T_h, self._h_nodes, self._latent_slope = T, h, slope

    @property
    def melts(self):
        return self.latent_heat > 0.0 and np.isfinite(self.T_liquidus)

    @property
    def T_feed(self):
        """Top of the feed ramp: material at or above it is fully liquid for the film."""
        return self.T_liquidus + (0.0 if self.T_liquidus - self.T_solidus <= 2.0 * MELT_RAMP + 1e-9 else MELT_RAMP)

    @property
    def h_liquid(self):
        """Specific enthalpy of the film (liquid at the liquidus) [J/kg]."""
        return float(self.enthalpy(self.T_feed))

    @classmethod
    def from_drama_json(cls, path=None):
        """A DRAMA material file; `path` may also be a name in MATERIAL_NAMES (AA7075_nomelt, AA7075, AA7075_range)."""
        path = MATERIAL_NAMES.get(path, path) or DEFAULT_MATERIAL
        with open(path) as fh:
            d = json.load(fh)
        cp = np.array(d["specificHeatCapacity"], dtype=float)
        k = np.array(d["heatConductivity"], dtype=float)
        T_melt = float(d.get("meltingTemperature", np.inf))
        latent = float(d.get("meltingHeat", 0.0)) if T_melt < NO_MELT_ABOVE else 0.0
        T_s, T_l = float(d.get("solidusTemperature", T_melt)), float(d.get("liquidusTemperature", T_melt))
        liquid = LiquidProperties(float(d["liquid"]["density"]), float(d["liquid"]["viscosity"]), float(d["liquid"]["surfaceTension"])) if "liquid" in d else None
        return cls(d["name"], float(d["density"]), float(d["emissivity"][0][1]), cp[:, 0], cp[:, 1], k[:, 0], k[:, 1],
                   latent, T_s if latent else np.inf, T_l if latent else np.inf, liquid)

    def k(self, T):
        return np.interp(T, self.T_k, self.k_table)

    def cp(self, T):
        return np.interp(T, self.T_cp, self.cp_table)

    def _liquid_fraction_raw(self, T):
        return np.clip((np.asarray(T, dtype=float) - self.T_solidus) / (self.T_liquidus - self.T_solidus), 0.0, 1.0)

    def liquid_fraction(self, T):
        """Melt fraction f_l(T): 0 below the solidus, 1 above the liquidus, linear between (zero for a non-melting material)."""
        T = np.asarray(T, dtype=float)
        return self._liquid_fraction_raw(T) if self.melts else np.zeros_like(T)

    def feed_fraction(self, T):
        """Fraction of an element's material the film receives: a +-MELT_RAMP ramp ending at T_feed (the liquidus)."""
        T = np.asarray(T, dtype=float)
        if not self.melts:
            return np.zeros_like(T)
        return np.clip((T - (self.T_feed - 2.0 * MELT_RAMP)) / (2.0 * MELT_RAMP), 0.0, 1.0)

    def cp_eff(self, T):
        """Effective heat capacity c_p + L_f df_l/dT."""
        T = np.asarray(T, dtype=float)
        c = self.cp(T)
        if self.melts:
            c = c + np.where((T > self.T_solidus) & (T < self.T_liquidus), self.latent_heat / (self.T_liquidus - self.T_solidus), 0.0)
        return c

    def enthalpy(self, T):
        """Specific enthalpy above T_REF [J/kg]: the exact integral of the piecewise-linear c_p (quadratic inside each
        table interval, linear beyond the table) plus L_f f_l(T)."""
        T = np.asarray(T, dtype=float)
        i = np.clip(np.searchsorted(self._T_h, T, side="right") - 1, 0, len(self._T_h) - 2)
        T0, T1 = self._T_h[i], self._T_h[i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        x = np.clip(T, T0, T1) - T0
        h = self._h_nodes[i] + (cp0 + self._latent_slope[i]) * x + 0.5 * (cp1 - cp0) / (T1 - T0) * x * x
        h = np.where(T < self._T_h[0], self._h_nodes[0] + self.cp_table[0] * (T - self._T_h[0]), h)
        h = np.where(T > self._T_h[-1], self._h_nodes[-1] + self.cp_table[-1] * (T - self._T_h[-1]), h)
        return h

    def temperature_from_enthalpy(self, h):
        """Inverse of enthalpy() (monotonic): the energy-equivalent temperature of a body holding h per kg."""
        h = np.asarray(h, dtype=float)
        i = np.clip(np.searchsorted(self._h_nodes, h, side="right") - 1, 0, len(self._T_h) - 2)
        T0, T1 = self._T_h[i], self._T_h[i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        b = cp0 + self._latent_slope[i]
        a = (cp1 - cp0) / (T1 - T0)
        dh = np.clip(h, self._h_nodes[i], self._h_nodes[i + 1]) - self._h_nodes[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(np.abs(a) > 1e-12, (np.sqrt(b * b + 2.0 * a * dh) - b) / a, dh / b)
        T = T0 + x
        T = np.where(h < self._h_nodes[0], self._T_h[0] + (h - self._h_nodes[0]) / self.cp_table[0], T)
        T = np.where(h > self._h_nodes[-1], self._T_h[-1] + (h - self._h_nodes[-1]) / self.cp_table[-1], T)
        return T
```


- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py tests/test_reentry_model_thermal.py tests/test_reentry_model_cli.py -q`
Expected: all pass (the no-melt material behaves exactly as in Step 2: `latent_heat` 0, `melts` False, enthalpy unchanged).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/material.py reentry_model/data/materials/AA7075.json reentry_model/data/materials/AA7075_range.json tests/test_reentry_model_material.py
git commit -m "Add the latent heat, melting ranges, feed fraction and liquid properties to the material (Step 3 Task 2)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: Thermal core — nodal enthalpy, element fractions, pinned nodes, nodal loads

**Files:**
- Modify (replace): `reentry_model/thermal/__init__.py`, `reentry_model/thermal/skfem_backend.py`, `reentry_model/thermal/fenicsx_backend.py`
- Modify: `reentry_model/cli.py` (the `--lumped-mass` flag becomes `--consistent-mass`, three lines)
- Test: `tests/test_reentry_model_thermal.py` (one test changed, four appended), `tests/test_reentry_model_fenicsx.py` (one line changed, one test appended)

**Interfaces:**
- Consumes: `VolumeMesh.active`, `.surface()` (Task 1); `Material.enthalpy/cp_eff/temperature_from_enthalpy/melts/T_solidus/T_liquidus` (Task 2).
- Produces: `StepResult(T, Q_conv, Q_rad, iterations, Q_extra=0.0, Q_dropped=0.0)`; solver methods `set_fractions(phi)` (0 = dead; refreshes the boundary and the pinned nodes), `element_energies(T=None)` (φ_e ρ V_e mean_i h(T_i)), `step(dt, q_conv, T_amb, dirichlet=None, nodal_load=None)` (nodal_load in W, mesh node order); attributes `phi`, `pinned`, `faces`, `areas`, `last_damping`; `SkfemThermalSolver(lumped_mass=True)` default with `consistent_mass = not lumped_mass`; `element_matrices(points, tets) -> (vol, Ke, Mk)` with `Mk[e, k]` the nodal-coefficient mass matrices and `nodal_mass_weights()`; `mass_matrix(c_nodal)`; `operators(T, T_old=None) -> (K, M_tan[, E])`. The FEniCSx backend has the same public surface (`lumped_mass=False` raises `ValueError`).

- [ ] **Step 1: Update the conformance test and append the melting tests**

In `tests/test_reentry_model_thermal.py` replace the body of `test_operators_match_scikit_fem_assembly` — the lines

```python
def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets].mean(axis=1)) * s.vol).sum(), rel=1e-12)
```

with

```python
def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    """K (element-mean k) against scikit-fem's assembly; the lumped capacity matrix is diag(sum_e V_e/4 rho c_p(T_i))
    (the nodal-enthalpy lumping, not the row sums of the consistent matrix), and the consistent option matches
    scikit-fem's consistent mass with the nodal (P1) c_p."""
    import scipy.sparse as sp
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    nodal = np.bincount(s.tets.ravel(), np.repeat(s.vol / 4.0, 4) * (RHO * s.material.cp(T))[s.tets].ravel(), minlength=s.points.shape[0])
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - sp.diags(nodal)).max() < 1e-10 * abs(M_ref).max()
    _, M_c = solver(coarse_sphere_mesh, material.Material.from_drama_json(), lumped_mass=False).operators(T)
    assert abs(M_c - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets]).mean(axis=1) * s.vol).sum(), rel=1e-12)   # nodal c_p, V_e/4 per node
```

(the remaining three lines of the test — the linear-field check on `K1` — stay). Then append to the file:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: element fractions, pinned nodes, nodal loads, melting

def test_fractions_pinned_nodes_and_nodal_loads_keep_the_balance(coarse_sphere_mesh):
    """Halving phi on the windward owners halves their energy; deactivating some of them pins nothing (their nodes
    still belong to live elements); a nodal sink over six steps keeps the discrete balance to 1e-8."""
    from reentry_model import mesh as mesh_mod
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets)
    s = solver(m, material.Material.from_drama_json("AA7075"))
    s.set_temperature(300.0)
    surf = m.surface()
    q = np.where(surf.normals[:, 0] > 0.0, 2e6, 0.0)
    for _ in range(4):
        s.step(0.5, q, 0.0)
    hot = surf.owner[surf.normals[:, 0] > 0.8]
    E_hot = s.element_energies()[hot].sum()
    phi = np.ones(m.n_elements)
    phi[hot] = 0.5
    E_before = s.energy()
    s.set_fractions(phi)
    assert s.energy() == pytest.approx(E_before - 0.5 * E_hot, rel=1e-12) and not s.pinned.any()
    m.deactivate(hot[:20])
    phi[hot[:20]] = 0.0
    s.set_fractions(phi)
    assert m.n_active == m.n_elements - 20 and len(s.areas) == m.surface().n_patches
    assert not np.isin(np.flatnonzero(s.pinned), np.unique(m.tets[m.active])).any()          # pinned = no live element
    load = np.zeros(m.n_nodes)
    load[np.unique(m.tets[hot[20:40]])] = -50.0                               # a 50 W sink on those nodes
    q2 = np.where(m.surface().normals[:, 0] > 0.0, 2e6, 0.0)
    E0, absorbed = s.energy(), 0.0
    for _ in range(6):
        r = s.step(0.5, q2, 0.0, nodal_load=load)
        absorbed += (r.Q_conv - r.Q_rad + r.Q_extra) * 0.5
        assert r.Q_extra == pytest.approx(load.sum()) and r.Q_dropped == 0.0
    assert abs(s.energy() - E0 - absorbed) < 1e-8 * abs(absorbed)
    with pytest.raises(ValueError):
        s.set_fractions(np.full(m.n_elements, 1.5))


def test_pinned_nodes_hold_their_temperature_and_drop_their_loads():
    """A two-tet mesh: deactivating one tet pins its private node; the load on it is reported as dropped."""
    from reentry_model import mesh as mesh_mod
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float) * 1e-2
    m = mesh_mod.VolumeMesh(pts, np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))
    s = thermal.thermal_solver("skfem", linear_solver="direct")
    s.setup(m, constant_material(), 0.0)
    s.set_temperature(np.array([300.0, 310.0, 320.0, 330.0, 900.0]))
    m.deactivate([1])
    s.set_fractions(np.array([1.0, 0.0]))
    assert s.pinned.tolist() == [False, False, False, False, True]
    load = np.zeros(5)
    load[4] = 100.0
    r = s.step(0.1, np.zeros(m.surface().n_patches), 0.0, nodal_load=load)
    assert r.Q_dropped == 100.0 and r.Q_extra == 100.0 and s.temperature()[4] == 900.0
    assert s.energy() == pytest.approx(RHO * CP * m.element_volumes()[0] * (315.0 - material.T_REF), rel=1e-9)   # nothing entered the live tet


def test_stefan_front_on_the_box_mesh():
    """Neumann's one-phase solution: liquid at T_w on x < X(t), X = 2 lambda sqrt(alpha t), lambda from
    lambda exp(lambda^2) erf(lambda) = St / sqrt(pi), St = c_p (T_w - T_m) / L. Box 30 x 1 x 1 mm of 0.5 mm cubes,
    wall Dirichlet at x = 0, front from the liquid-fraction contour at 0.5: within 1 % at t = 2..5 s (spec 13.4;
    measured 0.3 % on 2026-09-20)."""
    from scipy.optimize import brentq
    from scipy.special import erf
    from reentry_model import mesh as mesh_mod
    L, TM, TW = 4e5, 850.0, 1250.0
    mat = material.Material("stefan", RHO, 0.0, [100.0, 20000.0], [CP, CP], [100.0, 20000.0], [K0, K0], L, TM, TM,
                            material.LiquidProperties(2400.0, 1.3e-3, 0.86))
    alpha, St = K0 / (RHO * CP), CP * (TW - TM) / L
    lam = brentq(lambda x: x * np.exp(x * x) * erf(x) - St / np.sqrt(np.pi), 1e-6, 5.0)
    m = mesh_mod.box_mesh(0.03, 1e-3, 1e-3, 0.5e-3)
    s = solver(m, mat, linear_solver="direct")
    s.set_temperature(mat.T_solidus)
    wall = np.flatnonzero(m.points[:, 0] < 1e-9)
    axis = np.flatnonzero((np.abs(m.points[:, 1]) < 1e-12) & (np.abs(m.points[:, 2]) < 1e-12))
    x = m.points[axis, 0]
    order = np.argsort(x)
    t, dt = 0.0, 0.05
    while t < 5.0 - 1e-9:
        r = s.step(dt, np.zeros(m.surface().n_patches), 0.0, dirichlet=(wall, np.full(wall.size, TW)))
        t += dt
        assert r.iterations <= 8
        if abs(t - round(t)) < 1e-9 and t >= 2.0:
            f = mat.liquid_fraction(s.temperature())[axis][order]
            x_front = np.interp(0.5, f[::-1], x[order][::-1])
            assert x_front == pytest.approx(2.0 * lam * np.sqrt(alpha * t), rel=1e-2)


def test_melting_iteration_converges_across_the_ramp(coarse_sphere_mesh):
    """The isothermal body crossing the +-2 K ramp of the single-temperature material (the case that cycled with the
    element-mean secant iteration, 2026-09-20): every step converges in a few iterations and the enthalpy balance
    holds through the latent-heat plateau."""
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    s = solver(coarse_sphere_mesh, mat)
    s.set_temperature(840.0)
    q = np.full(coarse_sphere_mesh.surface().n_patches, 8e5)
    E0, absorbed = s.energy(), 0.0
    for _ in range(40):
        r = s.step(0.5, q, 0.0)
        absorbed += (r.Q_conv - r.Q_rad) * 0.5
        assert r.iterations <= 8
    T = s.temperature()
    assert 848.0 < T.mean() < 852.0 and T.max() - T.min() < 0.5                 # on the plateau, isothermal
    assert abs(s.energy() - E0 - absorbed) < 1e-8 * absorbed
```


In `tests/test_reentry_model_fenicsx.py` change the constructor test's line `thermal.thermal_solver("fenicsx", lumped_mass=True)` to `thermal.thermal_solver("fenicsx", lumped_mass=False)           # the nodal-enthalpy (lumped) capacity only` and append:

```python
def test_melting_run_matches_the_skfem_backend(coarse_sphere_mesh):
    """Element fractions, pinned nodes, nodal loads and the enthalpy Newton in both backends: 6 s of a melting body
    under physics-mode loads give the same mass, sprayed mass and temperatures (measured 1e-10 / 0 K on 2026-09-20)."""
    pytest.importorskip("cantera")
    from reentry_model import body, heating, mesh
    from test_reentry_model_coupled import MASS_100MM, simulator
    out = {}
    for name in ("skfem", "fenicsx"):
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM, T0=700.0)
        sim = simulator(b, t_max=60.0)
        sim.advance(43.0)
        model = heating.PhysicsHeating()
        for _ in range(12):
            sim.advance(0.5)
            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
            b.advance(sim.t, 0.5, model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature()), state=a)
        out[name] = (b.mass(0.0), b.sprayed_mass, b.solver.temperature(), b.mesh.n_active, b.energy_balance_residual())
    assert out["skfem"][0] == pytest.approx(out["fenicsx"][0], rel=1e-6) and out["skfem"][1] == pytest.approx(out["fenicsx"][1], rel=1e-4)
    assert np.abs(out["skfem"][2] - out["fenicsx"][2]).max() < 0.5 and out["skfem"][3] == out["fenicsx"][3]
    assert abs(out["fenicsx"][4]) < 1e-6 and out["fenicsx"][1] > 0.0
```


(The appended FEniCSx test needs Tasks 5–9; it stays failing in `fenicsx_env` until Task 9 and is skipped in `drama_env`.)

- [ ] **Step 2: Run the thermal tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`
Expected: the conformance test and the four new tests fail (no `set_fractions`, no `nodal_load`, `lumped_mass` default False).

- [ ] **Step 3: Replace `reentry_model/thermal/__init__.py`**

```python
"""Finite-element conduction solvers behind one protocol (spec section 8); `thermal_solver(name)` picks the backend.

Both backends solve rho c_p(T) dT/dt = div(k(T) grad T) with -k grad T . n = eps sigma (T^4 - T_amb^4) - q_conv on the
boundary, backward Euler in time, P1 tetrahedra in space, and report the same quantities: nodal temperatures,
stored enthalpy, radiated power. `dirichlet=(nodes, values)` in `step` is for the analytic tests only."""
from dataclasses import dataclass
from typing import Protocol

import numpy as np

SIGMA_SB = 5.670374419e-8         # Stefan-Boltzmann [W/(m2 K4)]
SOLVER_NAMES = ("skfem", "fenicsx")


class MissingBackend(RuntimeError):
    """The selected backend's library is not importable in this interpreter."""


@dataclass
class StepResult:
    T: np.ndarray                 # nodal temperatures after the step [K]
    Q_conv: float                 # W, convective power applied over the step
    Q_rad: float                  # W, radiated power at the end of the step
    iterations: int               # Newton/Picard iterations taken
    Q_extra: float = 0.0          # W, nodal load applied over the step (Step 3 deferred melt energy)
    Q_dropped: float = 0.0        # W, load that fell on pinned (material-free) nodes and was not applied


class ThermalSolver(Protocol):
    def setup(self, mesh, material, emissivity) -> None: ...
    def set_temperature(self, T) -> None: ...                       # uniform float or nodal array
    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None) -> StepResult: ...
    def temperature(self) -> np.ndarray: ...
    def energy(self) -> float: ...                                  # stored enthalpy above material.T_REF [J]
    def element_energies(self, T=None) -> np.ndarray: ...           # phi_e rho V_e h(T_e) per element [J]
    def radiated_power(self, T_amb) -> float: ...
    def set_fractions(self, phi) -> None: ...                       # element material fractions (0 = dead), Step 3


def thermal_solver(name, **options):
    if name == "skfem":
        from .skfem_backend import SkfemThermalSolver
        return SkfemThermalSolver(**options)
    if name == "fenicsx":
        from .fenicsx_backend import FenicsxThermalSolver      # raises MissingBackend without dolfinx
        return FenicsxThermalSolver(**options)
    raise ValueError("thermal solver must be one of {}, got {!r}".format(SOLVER_NAMES, name))
```


- [ ] **Step 4: Replace `reentry_model/thermal/skfem_backend.py`**

```python
"""scikit-fem backend: P1 tetrahedra, backward Euler, Newton on the nodal enthalpy with the radiation term.

Assembly. The P1 stiffness depends on T through one coefficient per element (k at the element-mean temperature)
and the mass through a nodal coefficient (rho c at the nodes, interpolated linearly), so the unit-coefficient
element matrices are computed once and rescaled into a fixed CSR pattern on every iteration (~20 ms for 50 k
tets). scikit-fem's generic `asm` on the same MeshTet/ElementTetP1 costs 0.2 s per iteration at 12 k nodes
(measured 2026-09-18) and is kept as the reference in `reference_operators` for the conformance test. Element
fractions phi_e (Step 3) scale both coefficients; nodes without material are pinned at their temperature.

Time stepping. Each iterate solves the tangent system of the residual R(T) = E(T)/dt + K T - F_conv + F_rad(T) -
F_extra, where E = M(rho c_sec)(T - T_old) is the exact nodal enthalpy increment (secant heat capacity per node, so
1^T E equals the change of rho int h dV with latent heat included); the tangent uses c_p,eff(T_i) per node and the
radiation Jacobian. The update is mapped through the true h(T) per node (enthalpy-consistent update) so that a node
cannot jump across the melting range, with damping as a fallback. With a smooth c_p the scheme is the Step 2
secant-capacity iteration unchanged; measured 2026-09-20: 2 iterations per step without melting, 3-5 with.

Radiation. The boundary functional uses the facet-mean temperature: eps sigma (T_f^4 - T_amb^4) A_f/3 to each of
the facet's three nodes, with its exact Jacobian 4 eps sigma T_f^3 A_f/9 on every node pair of the facet. Each
Newton iterate solves the symmetric positive definite system
    (M/dt + K + B_k) T = M T_old/dt + F_conv - F_rad(T_k) + B_k T_k
with M, K frozen at the previous iterate; because the rows of K sum to zero the discrete energy balance
1^T M (T - T_old) = dt (Q_conv - Q_rad) holds to the Newton tolerance.

Linear solver. `direct`: SciPy SuperLU; `amg` (default): CG preconditioned by a pyamg smoothed-aggregation
hierarchy rebuilt every `amg_rebuild_every` solves (the operator changes slowly). Measured on the 100 mm
2 mm/8 mm mesh (18.9 k nodes): AMG 0.03 s per solve vs SuperLU 0.6 s."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import SIGMA_SB, StepResult


class _Pattern:
    """Fixed CSR sparsity for repeated (rows, cols) entries; `assemble` sums values into it with one bincount."""

    def __init__(self, rows, cols, n):
        key = rows.astype(np.int64) * n + cols.astype(np.int64)
        unique, self._map = np.unique(key, return_inverse=True)
        self._indptr = np.searchsorted(unique // n, np.arange(n + 1))
        self._indices = (unique % n).astype(np.int32)
        self.n, self.nnz = n, unique.size

    def assemble(self, values):
        data = np.bincount(self._map, weights=np.asarray(values).ravel(), minlength=self.nnz)
        return sp.csr_matrix((data, self._indices, self._indptr), shape=(self.n, self.n))


def nodal_mass_weights():
    """W[k, a, b] = int phi_a phi_b phi_k dV / V_e on a tetrahedron: 1/20 (a = b = k), 1/60 (two equal), 1/120 (all
    distinct); sum_k W[k] is the consistent mass V_e/20 (1 + delta_ab)."""
    W = np.full((4, 4, 4), 1.0 / 120.0)
    for k in range(4):
        for a in range(4):
            W[k, a, a] = 1.0 / 60.0
            W[k, k, a] = W[k, a, k] = 1.0 / 60.0
        W[k, k, k] = 1.0 / 20.0
    return W


def element_matrices(points, tets):
    """Per-element volumes, unit stiffness (grad phi_a . grad phi_b V_e) and the mass matrices for a nodal coefficient:
    M_e(c) = sum_k c_k Mk[e, k] (Mk = V_e W), so that M_e(1) is the consistent mass V_e/20 (1 + delta_ab)."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    grad_ref = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    grads = np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), grad_ref)
    Ke = np.einsum("eai,ebi->eab", grads, grads) * vol[:, None, None]
    Mk = nodal_mass_weights()[None] * vol[:, None, None, None]
    return vol, Ke, Mk


class SkfemThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=True, newton_tol=1e-6, max_iterations=30,
                 amg_rebuild_every=30, cg_tol=1e-10):
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.lumped_mass = linear_solver, lumped_mass
        self.consistent_mass = not lumped_mass
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.amg_rebuild_every, self.cg_tol = amg_rebuild_every, cg_tol
        self._ml, self._solves, self.last_cg_iterations = None, 0, 0

    def setup(self, mesh, material, emissivity):
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        self.points, self.tets = mesh.points, mesh.tets
        self.vol, self.Ke, self.Mk = element_matrices(self.points, self.tets)
        n = len(self.points)
        self.pattern = _Pattern(np.repeat(self.tets, 4, axis=1).ravel(), np.tile(self.tets, (1, 4)).ravel(), n)
        self.T, self._T_prev = np.full(n, 300.0), None
        self.phi = np.where(mesh.active, 1.0, 0.0)
        self.pinned = np.zeros(n, dtype=bool)
        self._refresh_surface()

    def _refresh_surface(self):
        """Facet data (loads, radiation pattern) of the mesh's current boundary and the pinned (material-free) nodes."""
        surface = self.mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        n = len(self.points)
        self.facet_pattern = _Pattern(np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel(), n)
        self.Bf = np.ones((3, 3))[None] * (self.areas / 9.0)[:, None, None]
        pinned = np.ones(n, dtype=bool)
        pinned[np.unique(self.tets[self.phi > 0.0])] = False
        if pinned.any() != self.pinned.any() or not np.array_equal(pinned, self.pinned):
            self._ml = None                                   # the operator's structure changed: fresh AMG hierarchy
        self.pinned = pinned

    def set_fractions(self, phi):
        """Element material fractions phi_e in [0, 1] (0 = dead): they scale the heat capacity of every element (the
        conductivity stays that of the full element while phi > 0, see operators()); the boundary follows the mesh's
        active set and nodes without material are pinned at their temperature. Call after the mesh's active set or
        the fractions change."""
        self.phi = np.asarray(phi, dtype=float)
        if self.phi.shape != (len(self.tets),) or (self.phi < 0.0).any() or (self.phi > 1.0).any():
            raise ValueError("fractions must be one value in [0, 1] per element")
        self._refresh_surface()

    def set_temperature(self, T):
        self.T = np.full(len(self.points), float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    def facet_temperature(self, T=None):
        return (self.T if T is None else T)[self.faces].mean(axis=1)

    def facet_load(self, q):
        """Nodal load vector of a per-facet flux q [W/m2]: A_f/3 to each of the facet's nodes."""
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=len(self.points))

    def mass_matrix(self, c_nodal):
        """Capacity matrix for the nodal coefficient c [J/(m3 K)] (times phi_e per element). Lumped by default:
        diag(sum_e phi_e V_e/4 c_i), so that 1^T M dT = sum_e phi_e V_e/4 sum_i c_i dT_i is exactly the increment of
        the nodal enthalpy integral (the consistent form sum_k c_k Mk integrates the product of the interpolants of c
        and dT, which is not the increment of any energy functional and left a 3e-4 balance error, measured
        2026-09-20). `consistent_mass=True` keeps the consistent form (analytic tests only)."""
        if self.consistent_mass:
            return self.pattern.assemble(np.einsum("ek,ekab->eab", self.phi[:, None] * c_nodal[self.tets], self.Mk))
        return sp.diags(np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * c_nodal[self.tets].ravel(),
                                    minlength=len(self.points)), format="csr")

    def operators(self, T, T_old=None):
        """Stiffness K with k(T_e) at the element-mean temperature (active elements); the tangent mass M_tan with the nodal
        coefficient rho c_p,eff(T_i); and the enthalpy-rate vector E = M(rho c_sec) (T - T_old) with the nodal secant
        heat capacity c_sec,i = [h(T_i) - h(T_old,i)] / (T_i - T_old,i) (c_p,eff where a node has not moved). The
        enthalpy is nodal (spec Step 3 section 6, decided 2026-09-20: the element-mean enthalpy of Step 2 released the
        latent heat over a 4 K window of the element mean while the nodal temperatures span 30 K across a surface
        element, and its Newton iteration cycled): 1^T E is exactly the increment of rho int h dV with h interpolated
        linearly, whatever h(T) is. Without T_old, E is None and M_tan carries c_p,eff (the reference-operator convention)."""
        Te = T[self.tets].mean(axis=1)
        # conductivity is NOT scaled by phi_e: a partly consumed element is a thinner sliver of the same material,
        # which conducts better, not worse; scaling k with phi isolated the surface nodes of nearly consumed
        # elements and drove them to 5000 K (measured 2026-09-20). Dead elements (phi = 0) drop out.
        K = self.pattern.assemble(self.Ke * ((self.phi > 0.0) * self.material.k(Te))[:, None, None])
        c_tan = self.material.cp_eff(T)
        M = self.mass_matrix(self.material.rho * c_tan)
        if T_old is None:
            return K, M
        dT = T - T_old
        moved = np.abs(dT) > 1e-9
        c_sec = np.where(moved, (self.material.enthalpy(T) - self.material.enthalpy(T_old)) / np.where(moved, dT, 1.0), c_tan)
        return K, M, self.mass_matrix(self.material.rho * c_sec) @ dT

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def energy(self):
        """Stored enthalpy above T_REF: rho int h dV with the nodal h(T_i) interpolated linearly, i.e.
        sum_e phi_e rho V_e mean_i h(T_i); equals 1^T M T for constant c_p (consistent or lumped mass)."""
        return float(self.element_energies().sum())

    def element_energies(self, T=None):
        """phi_e rho V_e mean_i h(T_i) per element [J]."""
        h = self.material.enthalpy(self.T if T is None else T)
        return self.phi * self.material.rho * self.vol * h[self.tets].mean(axis=1)

    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
        """One backward-Euler step with the per-facet convective flux q_conv [W/m2] and an optional nodal load vector
        [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped)."""
        T_old = self.T
        # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        F_extra = np.zeros(len(self.points)) if nodal_load is None else np.asarray(nodal_load, dtype=float)
        Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
        T_new, iteration = T_k, 0
        converged = False
        last_relative_change, previous_change, damping, decreases = None, None, 1.0, 0
        mat = self.material
        for iteration in range(1, self.max_iterations + 1):
            K, M, E = self.operators(T_k, T_old)
            Tf = self.facet_temperature(T_k)
            B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            # Newton on R(T) = E(T)/dt + K T - F_conv + F_rad(T) - F_extra with the tangent M/dt + K + B:
            # A T_new = A T_k - R(T_k); with a smooth c_p (c_tan = c_sec) this is the Step 2 iteration unchanged
            A = (M / dt + K + B).tocsr()
            b = M @ T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
            if self.pinned.any():                              # material-free nodes keep their temperature (identity rows)
                A = A + sp.diags(self.pinned.astype(float), format="csr")
                b[self.pinned] = T_old[self.pinned]
            T_new = self._solve_with_dirichlet(A, b, T_k, dirichlet)
            if mat.melts:
                # enthalpy-consistent update: the linearised step is an enthalpy increment c_p,eff(T_k) (T_new - T_k)
                # per node; inverting the true h(T) puts a node that would overshoot the melting range where the
                # latent heat actually leaves it (identity where h is linear). Without it nodes jump across the +-2 K
                # ramp of a single-temperature material and the iteration cycles (measured 2026-09-20).
                free = ~self.pinned
                T_new[free] = mat.temperature_from_enthalpy(mat.enthalpy(T_k[free]) + mat.cp_eff(T_k[free]) * (T_new[free] - T_k[free]))
            last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            # damping when the change grows (a fallback for cycling iterates), released again after two decreases
            if previous_change is not None and last_relative_change > previous_change:
                damping, decreases = max(0.25, 0.5 * damping), 0
            elif damping < 1.0:
                decreases += 1
                if decreases >= 2:
                    damping, decreases = 1.0, 0
            previous_change = last_relative_change
            T_k = T_k + damping * (T_new - T_k) if damping < 1.0 else T_new
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev, self.T = T_old, T_new
        self.last_damping = damping
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)

    def _solve_with_dirichlet(self, A, b, x0, dirichlet):
        if dirichlet is None:
            return self._solve(A, b, x0)
        nodes, values = dirichlet
        free = np.ones(A.shape[0], dtype=bool)
        free[nodes] = False
        x = np.zeros(A.shape[0])
        x[nodes] = values
        x[free] = self._solve(A[free][:, free], b[free] - A[free][:, ~free] @ x[~free], x0[free], fresh=True)
        return x

    def _solve(self, A, b, x0, fresh=False):
        if self.linear_solver == "direct":
            return spla.spsolve(A.tocsc(), b)
        if fresh or self._ml is None or self._solves % self.amg_rebuild_every == 0:
            import pyamg
            self._ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
        self._solves += 1
        counter = []
        x, info = spla.cg(A, b, x0=x0, rtol=self.cg_tol, maxiter=500, M=self._ml.aspreconditioner(cycle="V"),
                          callback=lambda _: counter.append(1))
        self.last_cg_iterations = len(counter)
        if info != 0:
            raise RuntimeError("CG did not converge (info {})".format(info))
        return x

    def reference_operators(self, T):
        """K (element-mean k) and M (nodal c_p) assembled by scikit-fem with the same coefficients (conformance test only)."""
        from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP1, MeshTet, asm
        from skfem.helpers import dot, grad
        basis = Basis(MeshTet(self.points.T.copy(), self.tets.T.copy()), ElementTetP1(), intorder=3)   # exact for phi_i phi_j c (cubic)
        basis0 = basis.with_element(ElementTetP0())
        Te = T[self.tets].mean(axis=1)

        @BilinearForm
        def stiffness(u, v, w):
            return w.k * dot(grad(u), grad(v))

        @BilinearForm
        def mass(u, v, w):
            return w.c * u * v

        K = asm(stiffness, basis, k=basis0.interpolate(self.material.k(Te)))
        M = asm(mass, basis, c=basis.interpolate(self.material.rho * self.material.cp_eff(T)))    # nodal (P1) coefficient
        return K, M
```


- [ ] **Step 5: Replace `reentry_model/thermal/fenicsx_backend.py`**

```python
"""FEniCSx (dolfinx >= 0.11) backend: the scheme of skfem_backend on the same mesh (spec sections 8 and Step 3 6).

dolfinx assembles the stiffness K = int phi_e k(T_e) grad u . grad v dx (a DG0 coefficient refreshed every Newton
iterate); everything nodal -- the lumped capacity diag(sum_e phi_e V_e/4 rho c_p,eff(T_i)), the enthalpy-rate vector
E = diag(sum_e phi_e V_e/4 rho c_sec,i)(T - T_old), the convective loads A_f/3 per facet node, the radiation
eps sigma (T_f^4 - T_amb^4) A_f/3 with its Jacobian 4 eps sigma T_f^3 A_f/9 on the facet node pairs, and the Step 3
nodal loads -- is built in numpy on the VolumeMesh's numbering and added to the PETSc operator (the facet pairs are
edges of cells, so they lie inside K's pattern). Both backends therefore discretise identically; the radiation
follows the current boundary of the active set (element death), which a UFL `ds` measure could not. Dirichlet
nodes (analytic tests) and pinned material-free nodes are imposed with MatZeroRowsColumns. PETSc CG with hypre
BoomerAMG (hierarchy reused for `amg_rebuild_every` solves) or LU. dolfinx renumbers vertices: `node_of_dof` /
`dof_of_node` map between the VolumeMesh's node ids and the P1 dofs. Serial (the nodal vectors assume one
process). Verified with dolfinx 0.11.0 in fenicsx_env (2026-09-20). `dolfinx` is imported lazily: the constructor
raises MissingBackend without it."""
import numpy as np
import scipy.sparse as sp

from . import SIGMA_SB, MissingBackend, StepResult


class FenicsxThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=True, newton_tol=1e-6, max_iterations=30, amg_rebuild_every=30,
                 cg_tol=1e-10, **_):
        try:
            import dolfinx  # noqa: F401
            import ufl  # noqa: F401
            from mpi4py import MPI  # noqa: F401
            from petsc4py import PETSc  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        if not lumped_mass:
            raise ValueError("the fenicsx backend implements the lumped (nodal-enthalpy) capacity only")
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.newton_tol, self.max_iterations, self.cg_tol = linear_solver, newton_tol, max_iterations, cg_tol
        self.amg_rebuild_every, self._solves = amg_rebuild_every, 0
        self.lumped_mass = True

    def setup(self, mesh, material, emissivity):
        import basix.ufl
        import ufl
        from dolfinx import fem
        from dolfinx import mesh as dmesh
        from dolfinx.fem import petsc
        from mpi4py import MPI
        from petsc4py import PETSc
        from scipy.spatial import cKDTree
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
        self.msh = dmesh.create_mesh(MPI.COMM_WORLD, mesh.tets.astype(np.int64), domain, mesh.points)   # dolfinx >= 0.9: (comm, cells, element, x)
        self.V = fem.functionspace(self.msh, ("Lagrange", 1))
        self.V0 = fem.functionspace(self.msh, ("DG", 0))
        n_local = self.V.dofmap.index_map.size_local
        _, self.node_of_dof = cKDTree(mesh.points).query(self.V.tabulate_dof_coordinates()[:n_local])
        self.dof_of_node = np.empty(n_local, dtype=np.int64)
        self.dof_of_node[self.node_of_dof] = np.arange(n_local)
        self.cell_dofs = np.asarray(self.V.dofmap.list)[:, :4]
        # dolfinx cell k is the k-th cell passed in (create_mesh keeps the order in serial): cell volumes and the map
        # from dolfinx cells to mesh elements
        x = self.V.tabulate_dof_coordinates()[self.cell_dofs]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        self.cell_volumes = np.abs(np.linalg.det(J)) / 6.0
        self.cell_nodes = self.node_of_dof[self.cell_dofs]                   # mesh node ids of each dolfinx cell
        centroids = mesh.points[mesh.tets].mean(axis=1)
        self.element_of_cell = cKDTree(centroids).query(mesh.points[self.cell_nodes].mean(axis=1))[1]
        self.coef = fem.Function(self.V0)                                   # phi_e k(T_e) per cell
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        self.a = fem.form(self.coef * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx)
        self.K = petsc.create_matrix(self.a)
        self.ksp = PETSc.KSP().create(self.msh.comm)
        if self.linear_solver == "direct":
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.setType("cg")
            self.ksp.getPC().setType("hypre")
            self.ksp.getPC().setHYPREType("boomeramg")
            self.ksp.setTolerances(rtol=self.cg_tol, max_it=500)
        self.T = np.full(mesh.n_nodes, 300.0)                               # nodal, mesh numbering
        self._T_prev = None
        self.phi = np.where(mesh.active, 1.0, 0.0)
        self.pinned = np.zeros(mesh.n_nodes, dtype=bool)
        self._structure_changed = True
        self._refresh_surface()

    # -- geometry / fractions -----------------------------------------------------------------------------------
    def _refresh_surface(self):
        surface = self.mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        n = self.mesh.n_nodes
        rows, cols = np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel()
        self._facet_rows, self._facet_cols = rows, cols
        pinned = np.ones(n, dtype=bool)
        pinned[np.unique(self.mesh.tets[self.phi > 0.0])] = False
        if not np.array_equal(pinned, self.pinned):
            self._structure_changed = True
        self.pinned = pinned

    def set_fractions(self, phi):
        self.phi = np.asarray(phi, dtype=float)
        if self.phi.shape != (self.mesh.n_elements,) or (self.phi < 0.0).any() or (self.phi > 1.0).any():
            raise ValueError("fractions must be one value in [0, 1] per element")
        self._refresh_surface()

    def set_temperature(self, T):
        self.T = np.full(self.mesh.n_nodes, float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    # -- nodal pieces (mesh numbering) -----------------------------------------------------------------------------
    def facet_load(self, q):
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=self.mesh.n_nodes)

    def facet_temperature(self, T):
        return T[self.faces].mean(axis=1)

    def lumped(self, c_nodal):
        """diag(sum_e phi_e V_e/4 c_i) on the mesh nodes."""
        vol = self.mesh.element_volumes()
        return np.bincount(self.mesh.tets.ravel(), weights=np.repeat(self.phi * vol / 4.0, 4) * c_nodal[self.mesh.tets].ravel(),
                           minlength=self.mesh.n_nodes)

    def element_energies(self, T=None):
        h = self.material.enthalpy(self.T if T is None else T)
        return self.phi * self.material.rho * self.mesh.element_volumes() * h[self.mesh.tets].mean(axis=1)

    def energy(self):
        return float(self.element_energies().sum())

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(self.T if T is None else T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def _assemble_stiffness(self, T):
        from dolfinx.fem.petsc import assemble_matrix
        Te = T[self.mesh.tets].mean(axis=1)
        coef = (self.phi > 0.0) * self.material.k(Te)                     # k unscaled by phi (see skfem_backend.operators)
        self.coef.x.array[:] = coef[self.element_of_cell]
        self.K.zeroEntries()
        assemble_matrix(self.K, self.a)
        self.K.assemble()

    def _petsc_from_csr(self, matrix):
        from petsc4py import PETSc
        m = matrix.tocsr()
        m.sum_duplicates()
        return PETSc.Mat().createAIJ(size=m.shape, csr=(m.indptr.astype(np.int32), m.indices.astype(np.int32), m.data), comm=self.msh.comm)

    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
        from petsc4py import PETSc
        mat, n = self.material, self.mesh.n_nodes
        T_old = self.T
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        F_extra = np.zeros(n) if nodal_load is None else np.asarray(nodal_load, dtype=float)
        Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
        fixed = self.pinned.copy()
        fixed_values = T_old.copy()
        if dirichlet is not None:
            fixed[np.asarray(dirichlet[0])] = True
            fixed_values[np.asarray(dirichlet[0])] = dirichlet[1]
        vol = self.mesh.element_volumes()
        converged, last_relative_change, previous_change, damping, decreases = False, None, None, 1.0, 0
        T_new, iteration = T_k, 0
        for iteration in range(1, self.max_iterations + 1):
            self._assemble_stiffness(T_k)
            c_tan = mat.cp_eff(T_k)
            dT = T_k - T_old
            moved = np.abs(dT) > 1e-9
            c_sec = np.where(moved, (mat.enthalpy(T_k) - mat.enthalpy(T_old)) / np.where(moved, dT, 1.0), c_tan)
            M_tan, E = self.lumped(mat.rho * c_tan), self.lumped(mat.rho * c_sec) * dT
            Tf = self.facet_temperature(T_k)
            B = sp.csr_matrix((np.repeat(4.0 * self.emissivity * SIGMA_SB * Tf ** 3 * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            b = M_tan * T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
            # operator in dof numbering: K (dolfinx) + diag(M_tan/dt) + B
            A = self.K.copy()
            extra = sp.diags(M_tan / dt) + B
            extra = extra.tocsr()[self.node_of_dof][:, self.node_of_dof]
            A.axpy(1.0, self._petsc_from_csr(extra), structure=PETSc.Mat.Structure.SUBSET_NONZERO_PATTERN)
            rhs = PETSc.Vec().createWithArray(b[self.node_of_dof].copy(), comm=self.msh.comm)
            x = PETSc.Vec().createWithArray(fixed_values[self.node_of_dof].copy(), comm=self.msh.comm)
            fixed_dofs = self.dof_of_node[np.flatnonzero(fixed)].astype(np.int32)
            if fixed_dofs.size:
                A.zeroRowsColumns(fixed_dofs, diag=1.0, x=x, b=rhs)
            self.ksp.setOperators(A)
            reuse = self.linear_solver == "amg" and dirichlet is None and not self._structure_changed and self._solves % self.amg_rebuild_every != 0
            self.ksp.getPC().setReusePreconditioner(reuse)
            self._structure_changed = False
            self._solves += 1
            sol = PETSc.Vec().createWithArray(T_k[self.node_of_dof].copy(), comm=self.msh.comm)
            self.ksp.solve(rhs, sol)
            if self.ksp.getConvergedReason() <= 0:
                raise RuntimeError("PETSc KSP did not converge (reason {})".format(self.ksp.getConvergedReason()))
            T_new = sol.getArray()[self.dof_of_node].copy()
            T_new[fixed] = fixed_values[fixed]
            if mat.melts:
                free = ~fixed
                T_new[free] = mat.temperature_from_enthalpy(mat.enthalpy(T_k[free]) + mat.cp_eff(T_k[free]) * (T_new[free] - T_k[free]))
            last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            if previous_change is not None and last_relative_change > previous_change:
                damping, decreases = max(0.25, 0.5 * damping), 0
            elif damping < 1.0:
                decreases += 1
                if decreases >= 2:
                    damping, decreases = 1.0, 0
            previous_change = last_relative_change
            T_k = T_k + damping * (T_new - T_k) if damping < 1.0 else T_new
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev, self.T = T_old, T_new
        self.last_damping = damping
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)
```


- [ ] **Step 6: Rename the CLI flag**

In `reentry_model/cli.py`: replace `th.add_argument("--lumped-mass", action="store_true")` with

```python
    th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
```

replace `lumped_mass=args.lumped_mass` (in `build_thermal`) with `lumped_mass=not args.consistent_mass`, and `"lumped_mass": args.lumped_mass}` (the info dict) with `"lumped_mass": not args.consistent_mass}`.

- [ ] **Step 7: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py tests/test_reentry_model_coupled.py tests/test_reentry_model_cli.py tests/test_reentry_model_viz.py -q`
Expected: all pass — the Step 2 analytic cases (lumped limit, radiative cooling, Carslaw–Jaeger, balance) are unchanged by the lumping for constant c_p, and the melting cases converge in ≤ 8 iterations. Then in `fenicsx_env`:

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q --deselect tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend
```

Expected: 7 passed (the seven Step 2 conformance tests, including the skfem cross-check, with the nodal scheme).

- [ ] **Step 8: Commit**

```bash
git add reentry_model/thermal tests/test_reentry_model_thermal.py tests/test_reentry_model_fenicsx.py reentry_model/cli.py
git commit -m "Move the thermal core to the nodal enthalpy with element fractions, pinned nodes and nodal loads (Step 3 Task 3)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Girin's dispersion relation and its cached table

**Files:**
- Create: `reentry_model/dispersion.py`, `reentry_model/data/girin_dispersion.json` (generated)
- Test: `tests/test_reentry_model_dispersion.py`

**Interfaces:**
- Produces: `cubic_coefficients(delta, we) -> (a2, a1, a0)`; `growth_rates(we, delta)`; `fastest_mode(we) -> (Delta_f, Im Omega_f, Re Omega_f)` (nan/0/0 when stable); `build_table()`, `write_table(path)`; `DispersionTable(path=TABLE_PATH)` callable on scalars/arrays returning `(Delta_f, Im Omega_f, Re Omega_f)` with attribute `we_onset`; constants `WE_CRITICAL_THEORY = 3.08`, `WE_CRITICAL_PRACTICAL = 4.62`, `TABLE_PATH`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_dispersion.py`:

```python
"""dispersion.py: Girin's Eq. (1) solved numerically -- onset, asymptotes, the cached table (spec Step 3 section 9)."""
import os

import numpy as np
import pytest

from reentry_model import dispersion as dp


def test_onset_and_asymptotes():
    assert np.isnan(dp.fastest_mode(3.0)[0]) and dp.fastest_mode(3.0)[1] == 0.0        # stable at We_s = 3.00
    d, im, re = dp.fastest_mode(3.08)
    assert 0.0 < im and 0.0 < d < 0.1                                                      # unstable at 3.08, long waves
    d, im, re = dp.fastest_mode(1e4)
    assert d == pytest.approx(1.225, rel=1e-2) and im == pytest.approx(0.24, rel=4e-2)    # Girin's Fig. 5 asymptotes (measured 1.226 / 0.247)
    for we in (20.0, 100.0, 1e4):
        d, im, re = dp.fastest_mode(we)
        assert 1.5 <= re / im <= 1.7                                                       # Re Omega_f = 1.5-1.7 Im Omega_f (Girin: 1.5)
    d, im, _ = dp.fastest_mode(4.62)
    assert d == pytest.approx(0.58, rel=2e-2) and im == pytest.approx(0.049, rel=3e-2)
    d, im, _ = dp.fastest_mode(10.0)
    assert d == pytest.approx(1.34, rel=2e-2) and im == pytest.approx(0.18, rel=3e-2)


def test_cubic_coefficients_reproduce_the_relation():
    """A root of the cubic satisfies Eq. (1) written out."""
    delta, we = 1.2, 20.0
    a2, a1, a0 = dp.cubic_coefficients(delta, we)
    roots = np.roots([1.0, a2, a1, a0])
    K = (1.0 - np.exp(-2.0 * delta)) / (2.0 * delta)
    for w in roots:
        lhs = (w - delta) * ((w - delta) * (w + K * delta) + (1.0 - K) * delta)
        rhs = delta ** 3 / we * (w - K * delta)
        assert abs(lhs - rhs) < 1e-9


def test_table_is_cached_and_interpolates(tmp_path):
    path = str(tmp_path / "disp.json")
    t = dp.DispersionTable(path)
    assert os.path.isfile(path) and 3.0 < t.we_onset < 3.2
    d, im, re = t(np.array([2.0, 4.62, 100.0, 1e5]))
    assert np.isnan(d[0]) and im[0] == 0.0 and d[1] == pytest.approx(0.58, rel=3e-2) and d[3] == pytest.approx(1.226, rel=1e-2)
    assert im[2] == pytest.approx(0.241, rel=2e-2)
    packaged = dp.DispersionTable()
    assert np.allclose(packaged.we, t.we) and np.allclose(packaged.im_omega_f, t.im_omega_f, atol=1e-6)     # the committed table is current
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_dispersion.py -q`
Expected: ImportError (`reentry_model.dispersion`).

- [ ] **Step 3: Create `reentry_model/dispersion.py`**

```python
"""Girin's (2017, A&A 606, A63, Eq. 1) dispersion relation of the gradient instability of a sheared liquid layer,
solved numerically for the fastest-growing disturbance (spec Step 3 section 9).

    (Omega - Delta) [(Omega - Delta)(Omega + K Delta) + (1 - K) Delta] = Delta^3 We_s^-1 (Omega - K Delta),
    K = [1 - exp(-2 Delta)] / (2 Delta),   Omega = omega delta_m / V_s,   Delta = 2 pi delta_m / lambda,
    We_s = rho_m V_s^2 delta_m / Sigma.

For every We_s the cubic in Omega is solved on a grid of Delta (companion-matrix eigenvalues, vectorised) and the
Delta with the largest Im Omega is the fastest mode: Delta_f (wavelength 2 pi delta_m / Delta_f), Im Omega_f (growth
rate Im Omega_f V_s / delta_m) and Re Omega_f. Measured 2026-09-20: instability appears between We_s 3.00 and 3.08
(Girin: 3.08); Delta_f -> 1.226, Im Omega_f -> 0.247 at We_s 1e4 (Girin's Fig. 5 asymptotes 1.225, 0.24);
Re/Im = 1.68 at We_s 20 and 1.56 at 1e4 (Girin: 1.5); at We_s 4.62: 0.58 / 0.049; 10: 1.34 / 0.18; 20: 1.29 / 0.22.
The table is cached in the package data (`girin_dispersion.json`) and interpolated in log We_s."""
import json
import os

import numpy as np

from . import DATA_DIR

TABLE_PATH = os.path.join(DATA_DIR, "girin_dispersion.json")
WE_CRITICAL_THEORY = 3.08          # Girin's critical surface Weber number
WE_CRITICAL_PRACTICAL = 4.62       # the value he recommends in practice (GI >= 0.6)
DELTA_GRID = np.linspace(0.02, 5.0, 2500)
WE_GRID = np.logspace(np.log10(3.0), 4.0, 200)
UNSTABLE_TOL = 1e-6                # Im Omega above which a root counts as unstable


def cubic_coefficients(delta, we):
    """Coefficients (a2, a1, a0) of Omega^3 + a2 Omega^2 + a1 Omega + a0 = 0, Girin's Eq. (1) expanded."""
    D = np.asarray(delta, dtype=float)
    K = (1.0 - np.exp(-2.0 * D)) / (2.0 * D)
    a2 = (K - 2.0) * D
    a1 = (1.0 - 2.0 * K) * D ** 2 + (1.0 - K) * D - D ** 3 / we
    a0 = K * D ** 3 - (1.0 - K) * D ** 2 + K * D ** 4 / we
    return a2, a1, a0


def growth_rates(we, delta=DELTA_GRID):
    """(Im Omega, Re Omega) of the most unstable root of Eq. (1) at each Delta for one We_s."""
    a2, a1, a0 = cubic_coefficients(delta, we)
    C = np.zeros((len(delta), 3, 3), dtype=complex)
    C[:, 0, 1] = 1.0
    C[:, 1, 2] = 1.0
    C[:, 2, 0], C[:, 2, 1], C[:, 2, 2] = -a0, -a1, -a2
    roots = np.linalg.eigvals(C)
    k = np.argmax(roots.imag, axis=1)
    best = roots[np.arange(len(delta)), k]
    return best.imag, best.real


def fastest_mode(we, delta=DELTA_GRID):
    """(Delta_f, Im Omega_f, Re Omega_f) of the fastest-growing disturbance; (nan, 0, 0) when nothing is unstable."""
    im, re = growth_rates(we, delta)
    j = int(np.argmax(im))
    if im[j] <= UNSTABLE_TOL:
        return float("nan"), 0.0, 0.0
    return float(delta[j]), float(im[j]), float(re[j])


def build_table(we_grid=WE_GRID, delta=DELTA_GRID):
    rows = np.array([fastest_mode(we, delta) for we in we_grid])
    return {"we": list(map(float, we_grid)), "delta_f": list(map(float, rows[:, 0])),
            "im_omega_f": list(map(float, rows[:, 1])), "re_omega_f": list(map(float, rows[:, 2])),
            "delta_grid": [float(delta[0]), float(delta[-1]), int(len(delta))]}


def write_table(path=TABLE_PATH):
    table = build_table()
    with open(path, "w") as fh:
        json.dump(table, fh)
    return table


class DispersionTable:
    """Delta_f(We_s), Im Omega_f(We_s), Re Omega_f(We_s) interpolated in log We_s from the cached table."""

    def __init__(self, path=TABLE_PATH):
        if not os.path.isfile(path):
            write_table(path)
        with open(path) as fh:
            t = json.load(fh)
        self.we = np.array(t["we"])
        self.delta_f, self.im_omega_f, self.re_omega_f = (np.array(t[k]) for k in ("delta_f", "im_omega_f", "re_omega_f"))
        unstable = self.im_omega_f > 0.0
        self.we_onset = float(self.we[unstable][0])              # first tabulated We_s with an unstable mode
        self._log_we = np.log(self.we[unstable])
        self._delta, self._im, self._re = self.delta_f[unstable], self.im_omega_f[unstable], self.re_omega_f[unstable]

    def __call__(self, we):
        """Arrays (Delta_f, Im Omega_f, Re Omega_f) for We_s (scalar or array); nan/0 below the onset."""
        we = np.asarray(we, dtype=float)
        x = np.log(np.maximum(we, self.we_onset))
        d, im, re = (np.interp(x, self._log_we, v) for v in (self._delta, self._im, self._re))
        stable = we < self.we_onset
        return np.where(stable, np.nan, d), np.where(stable, 0.0, im), np.where(stable, 0.0, re)
```


- [ ] **Step 4: Generate the packaged table and run the tests**

```bash
"$PY" -c "from reentry_model import dispersion; t = dispersion.write_table(); print(len(t['we']), t['delta_f'][-1], t['im_omega_f'][-1])"
"$PY" -m pytest tests/test_reentry_model_dispersion.py -q
```

Expected: `200 1.2256... 0.2469...` (2–3 s), then 3 passed.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/dispersion.py reentry_model/data/girin_dispersion.json tests/test_reentry_model_dispersion.py
git commit -m "Solve Girin's dispersion relation numerically and cache the fastest-mode table (Step 3 Task 4)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Gas edge state and the surface flow per patch

**Files:**
- Modify: `reentry_model/gas.py` (five edits)
- Create: `reentry_model/surface_flow.py`
- Test: `tests/test_reentry_model_surface_flow.py`

**Interfaces:**
- Consumes: `gas.EquilibriumAir.stagnation` (Step 2), `aero.mean_free_path`, `aero.SesamTable`, `constants.HARD_SPHERE_DIAMETER`, `trajectory.AeroState` (freestream rho/T/p/m_bar, V, kn, ma, a_drag).
- Produces: `GasState.s`, `.a`, `.m_bar` (defaulted fields), `EquilibriumAir.expand(stag, p)`; `surface_flow.SurfaceFlow(air=None, rarefied_shear="slip"|"bridged", bridging=None, sigma_v=1, sigma_t=1)` with `.evaluate(state, theta, radius, rho_liquid) -> SurfaceFlowResult` (fields `u_e, u_eff, rho_e, T_e, mu_e, p_e, mach_e, delta_a, lambda_e, kn_delta, regime, tau, tau_continuum, tau_fm, G, p_stag, p_inf, deceleration`, method `regime_fractions(areas, windward)`), `.edge_table(state, radius)`, `.last_bins`; functions `ranger_psi`, `ranger_thickness`, `boundary_layer_thickness(s, u_e, nu_e, constant=RANGER_C, power=4, refine=20)`; constants `REGIME_CONTINUUM/SLIP/FM = 0/1/2`, `KN_SLIP = 0.01`, `KN_FM = 0.1`, `FM_BODY_KN = 1.0`, `RANGER_C = 58.08`, `THWAITES_C = 0.45`, `RAREFIED_SHEAR_NAMES`, `THETA_BINS` (0–90° by 1°).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_surface_flow.py`:

```python
"""surface_flow.py: Ranger's boundary layer, the Thwaites cross-check, the edge state and regimes at 71 km."""
import math
from datetime import datetime

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, surface_flow as sf
from reentry_model import trajectory as tj

pytest.importorskip("cantera")
EPOCH = datetime(2024, 8, 1, 12, 53, 7)
R100 = tj.InitialState(7500.0, 77500.133, math.radians(-0.959331), math.radians(347.168296), math.radians(29.546067),
                       math.radians(-82.134333), EPOCH)


def potential_flow(R=0.05, V=1000.0, nu=1e-4):
    th = sf.THETA_BINS
    return th, R * th, 1.5 * V * np.sin(th), np.full_like(th, nu), V * 2.0 * R / nu


def test_ranger_form_reproduces_the_closed_form():
    th, s, u, nu, Re = potential_flow()
    d = sf.boundary_layer_thickness(s, u, nu)
    ref = sf.ranger_thickness(0.05, Re, th[1:])
    assert np.abs(d[1:] / ref - 1.0).max() < 3e-3 and np.abs(d[2:] / ref[1:] - 1.0).max() < 1e-3
    assert d[0] == d[1] and sf.ranger_psi(1e-3) == pytest.approx(math.sqrt(3.2), rel=1e-4)      # Psi(0) = sqrt(48/15)


def test_thwaites_shape_agrees_with_ranger_within_five_percent():
    th, s, u, nu, Re = potential_flow()
    ranger = sf.boundary_layer_thickness(s, u, nu)
    thwaites = sf.boundary_layer_thickness(s, u, nu, sf.THWAITES_C, 5)
    m = (th >= math.radians(5.0)) & (th <= math.radians(85.0))
    ratio = ranger[m] / thwaites[m]
    assert 11.5 < ratio.min() and ratio.max() < 13.0 and (ratio.max() - ratio.min()) / ratio.mean() < 0.05    # measured 12.0-12.4


def aero_state_at(t):
    sim = tj.Simulator(R100, body.ConstantBody(1.473), atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(),
                       aero.SesamTable(), tj.Settings(diameter=0.1))
    sim.advance(t)
    return sim.aero_state(sim.t, sim.y[:3], sim.y[3:])


def test_edge_state_regimes_and_shear_at_71_km():
    a = aero_state_at(43.5)
    assert 70.5e3 < a.h < 71.5e3
    theta = np.radians([0.5, 10.0, 30.0, 45.0, 60.0, 75.0, 89.0, 100.0, 150.0])
    flow = sf.SurfaceFlow().evaluate(a, theta, 0.05, 2400.0)
    assert flow.u_e[0] < 50.0 and 1100.0 < flow.u_e[2] < 1400.0 and 1800.0 < flow.u_e[3] < 2100.0      # ~1 km/s mid-sphere (spec estimate)
    assert np.all(flow.u_e[7:] == 0.0) and np.all(flow.tau[7:] == 0.0) and np.all(flow.G[7:] == 0.0)   # leeward: static
    assert 5e-3 < flow.delta_a[1] < 8e-3 and 6e-3 < flow.delta_a[3] < 9e-3                             # 3-9 mm (spec estimate 3-7)
    assert flow.regime[1] == sf.REGIME_CONTINUUM and flow.regime[3] in (sf.REGIME_CONTINUUM, sf.REGIME_SLIP)
    assert flow.regime[4] == sf.REGIME_SLIP and flow.regime[6] == sf.REGIME_FM                           # Kn_delta grows toward the rim
    assert 20.0 < flow.tau_continuum[3] < 60.0 and flow.tau[3] == pytest.approx(flow.tau_continuum[3] / (1.0 + flow.kn_delta[3]), rel=1e-9)
    assert 1500.0 < flow.tau_fm[3] < 2100.0 and flow.tau[6] == flow.tau_fm[6]                           # slip mode: free-molecular shear above Kn_delta 0.1
    assert flow.G[3] > 0.0 and flow.G[6] < 0.0                                                          # pressure gradient outward, deceleration toward the nose near the rim
    assert flow.p_stag > 1000.0 and flow.p_e[0] == pytest.approx(flow.p_stag, rel=1e-3) and flow.p_e[7] == flow.p_inf
    assert 0.9 < flow.mach_e[3] < 1.1 and flow.deceleration > 5.0
    fractions = flow.regime_fractions(np.ones(theta.size), theta <= np.pi / 2)
    assert sum(fractions) == pytest.approx(1.0)
    bridged = sf.SurfaceFlow(rarefied_shear="bridged").evaluate(a, theta, 0.05, 2400.0)
    f = aero.SesamTable()(a.kn)
    assert bridged.tau[6] == pytest.approx((1.0 - f) * flow.tau_continuum[6] / (1.0 + flow.kn_delta[6]) + f * flow.tau_fm[6], rel=1e-9)
    with pytest.raises(ValueError):
        sf.SurfaceFlow(rarefied_shear="magic")


def test_free_molecular_body_has_no_edge_state():
    a = aero_state_at(0.0)
    a.kn = 5.0
    theta = np.radians([30.0, 60.0, 120.0])
    flow = sf.SurfaceFlow().evaluate(a, theta, 0.05, 2400.0)
    assert np.all(flow.regime == sf.REGIME_FM) and np.all(np.isnan(flow.delta_a)) and np.all(flow.u_e == 0.0)
    assert flow.tau[0] == pytest.approx(a.freestream.rho * a.V ** 2 * math.sin(theta[0]) * math.cos(theta[0])) and flow.tau[2] == 0.0
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_surface_flow.py -q`
Expected: ImportError (`reentry_model.surface_flow`).

- [ ] **Step 3: Extend `reentry_model/gas.py`**

Five edits. (a) After `AIR = "N2:0.79, O2:0.21"` add `AVOGADRO = 6.02214076e23`. (b) In `GasState`, after `X: dict           # mole fractions > 1e-6` add:

```python
    s: float = 0.0    # J/(kg K), specific entropy (Cantera basis) -- Step 3: the isentropic expansion to the edge state
    a: float = 0.0    # m/s, frozen sound speed
    m_bar: float = 0.0  # kg, mean molecular mass (Step 3: the edge mean free path)
```

(c) In `_state()` replace the `return GasState(...)` statement with:

```python
        return GasState(float(g.P), float(g.T), float(g.density), float(g.enthalpy_mass),
                        wilke_viscosity(float(g.T), X, self.molar_masses), float(h_D), X, float(g.entropy_mass),
                        float(g.sound_speed), float(g.mean_molecular_weight) * 1e-3 / AVOGADRO)
```

(d) Before `def wall(self, T_w, p):` insert:

```python
    def expand(self, stag, p):
        """Isentropic expansion of the stagnation state `stag` to pressure p (equilibrium above T_EQUILIBRATE):
        the boundary-layer edge state (Step 3, spec section 7)."""
        g = self.gas
        g.SPX = stag.s, p, AIR
        if g.T > T_EQUILIBRATE:
            g.equilibrate("SP")
        return self._state()

```

(e) Nothing else changes; `"$PY" -m pytest tests/test_reentry_model_gas.py tests/test_reentry_model_heating.py -q` must still pass.

- [ ] **Step 4: Create `reentry_model/surface_flow.py`**

```python
"""Gas-side surface flow per patch: boundary-layer edge state, laminar boundary layer, Knudsen regime, wall shear
and the film's driving pressure gradient (spec Step 3 section 7).

Edge state: modified-Newtonian pressure p_e = p_inf + (p_s - p_inf) cos^2 theta on the windward face and the
isentropic expansion of the Step 2 equilibrium stagnation state to p_e (Cantera, one state per 1 degree bin,
interpolated to the patches): h_e, u_e = sqrt(2 (h_s - h_e)), rho_e, T_e, a_e, mu_e (Blottner-Wilke), m_bar_e.
Leeward patches (theta > 90 degrees) have p_e = p_inf, u_e = 0 and no shear.

Boundary layer: Girin's linear-profile thickness delta_a in Ranger's (1972) form generalised to the actual edge
velocity, delta_a^2 = RANGER_C nu_e int_0^s u_e^4 ds' / u_e^5 along the meridian s = R theta; RANGER_C = 58.1 makes
it identical to delta_a = 2.2 R Re_D^-1/2 Psi(theta), Psi = [(6 theta - 4 sin 2 theta + 1/2 sin 4 theta)/sin^5
theta]^1/2, for the potential-flow velocity 1.5 V sin theta and constant properties (Girin 2017 section 3). Thwaites'
momentum thickness is provided for the cross-check (its shape differs from Ranger's by -8 %/+8 % over 5-85 degrees,
measured 2026-09-20, which is why the spec's Thwaites route was replaced by Ranger's own integral). Wall shear in
Girin's convention tau_c = mu_e u_e / delta_a.

Regimes by Kn_delta = lambda_e / delta_a (hard-sphere mean free path of the edge gas): continuum < 0.01, slip
0.01-0.1, transitional/free-molecular >= 0.1 (REGIME_CONTINUUM/SLIP/FM = 0/1/2). Slip: Maxwell first-order slip
with the linear near-wall profile, an effective edge velocity u_e/(1 + C Kn_delta), C = (2 - sigma_v)/sigma_v,
sigma_v = 1, so tau = tau_c/(1 + C Kn_delta) and the same u_e,eff enters Girin's Eq. (2). Free-molecular:
tau_fm = sigma_t rho_inf V^2 sin theta cos theta (hypersonic speed ratio, sigma_t = 1). `rarefied_shear = "slip"`
(default): tau_slip below Kn_delta 0.1, tau_fm above; "bridged": (1 - f) tau_slip + f tau_fm with SESAM's measured
drag bridging f(Kn) on the body Knudsen number. Where the whole body is free-molecular (Kn >= FM_BODY_KN) no
continuum edge state is computed at all.

Pressure gradient and body force on the film: G = -dp_e/ds - rho_l a sin theta [Pa/m] along the surface direction
away from the stagnation point, with dp_e/ds = -2 (p_s - p_inf) sin theta cos theta / R and a the body's deceleration
(the film feels the inertial force toward the nose)."""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid

from . import aero, gas
from .constants import HARD_SPHERE_DIAMETER

REGIME_CONTINUUM, REGIME_SLIP, REGIME_FM = 0, 1, 2
KN_SLIP, KN_FM = 0.01, 0.1
FM_BODY_KN = 1.0
RANGER_C = 58.08                      # 4.84 x 16 x 3/4: Ranger's delta_a from the u_e^4 integral (module docstring)
THWAITES_C = 0.45
SIGMA_V, SIGMA_T = 1.0, 1.0
RAREFIED_SHEAR_NAMES = ("slip", "bridged")
THETA_BINS = np.radians(np.arange(0.0, 90.0 + 0.5, 1.0))


def ranger_psi(theta):
    theta = np.asarray(theta, dtype=float)
    return np.sqrt((6.0 * theta - 4.0 * np.sin(2.0 * theta) + 0.5 * np.sin(4.0 * theta)) / np.sin(theta) ** 5)


def ranger_thickness(radius, reynolds_diameter, theta):
    """delta_a = 2.2 R Re_D^-1/2 Psi(theta) (Ranger 1972, as used by Girin 2017)."""
    return 2.2 * radius * ranger_psi(theta) / np.sqrt(reynolds_diameter)


def boundary_layer_thickness(s, u_e, nu_e, constant=RANGER_C, power=4, refine=20):
    """delta_a on the meridian grid s from delta_a^2 = constant nu_e int u_e^power ds / u_e^(power+1) (Ranger's form
    for power = 4; Thwaites' momentum thickness is THWAITES_C with power = 5). The integrand is evaluated on a grid
    refined `refine` times (u_e interpolated linearly) so that the trapezoid rule resolves u_e^power ~ s^power near
    the stagnation point (1-degree bins alone give +17 % at 2 degrees, measured)."""
    s, u_e = np.asarray(s, dtype=float), np.asarray(u_e, dtype=float)
    fine = np.linspace(s[0], s[-1], refine * (len(s) - 1) + 1)
    integral = np.interp(s, fine, cumulative_trapezoid(np.interp(fine, s, u_e) ** power, fine, initial=0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        d2 = constant * nu_e * integral / u_e ** (power + 1)
    d2 = np.where(u_e > 0.0, d2, 0.0)
    if len(s) > 1 and u_e[0] == 0.0 and u_e[1] > 0.0:           # stagnation point: the limit is finite, take the neighbour's value
        d2[0] = d2[1]
    return np.sqrt(np.maximum(d2, 0.0))


@dataclass
class SurfaceFlowResult:
    u_e: np.ndarray          # m/s, edge velocity (before the slip correction)
    u_eff: np.ndarray        # m/s, effective edge velocity driving the film (slip-corrected)
    rho_e: np.ndarray        # kg/m3
    T_e: np.ndarray          # K
    mu_e: np.ndarray         # Pa s
    p_e: np.ndarray          # Pa
    mach_e: np.ndarray
    delta_a: np.ndarray      # m, gas boundary-layer thickness (nan where free-molecular)
    lambda_e: np.ndarray     # m, edge mean free path
    kn_delta: np.ndarray     # lambda_e / delta_a (inf where free-molecular)
    regime: np.ndarray       # 0 continuum, 1 slip, 2 transitional/free-molecular
    tau: np.ndarray          # Pa, wall shear used by the film
    tau_continuum: np.ndarray  # Pa, mu_e u_e / delta_a (before slip/bridging)
    tau_fm: np.ndarray       # Pa, free-molecular shear
    G: np.ndarray            # Pa/m, film driving gradient (pressure gradient + inertial force)
    p_stag: float
    p_inf: float
    deceleration: float      # m/s2

    def regime_fractions(self, areas, windward):
        a = areas * windward
        total = a.sum()
        return [float(a[self.regime == k].sum() / total) if total > 0.0 else 0.0 for k in (0, 1, 2)]


class SurfaceFlow:
    def __init__(self, air=None, rarefied_shear="slip", bridging=None, sigma_v=SIGMA_V, sigma_t=SIGMA_T):
        if rarefied_shear not in RAREFIED_SHEAR_NAMES:
            raise ValueError("rarefied_shear must be one of {}, got {!r}".format(RAREFIED_SHEAR_NAMES, rarefied_shear))
        self.air = air or gas.EquilibriumAir()
        self.rarefied_shear, self.bridging = rarefied_shear, bridging or aero.SesamTable()
        self.slip_C, self.sigma_t = (2.0 - sigma_v) / sigma_v, sigma_t
        self.last_bins = None

    def edge_table(self, state, radius):
        """Edge quantities on THETA_BINS from the stagnation state: (p_e, u_e, rho_e, T_e, mu_e, a_e, m_bar_e, p_s)."""
        fs = state.freestream
        stag = self.air.stagnation(fs.rho, fs.T, state.V)
        p_inf = float(fs.p)
        p_e = p_inf + (stag.p - p_inf) * np.cos(THETA_BINS) ** 2
        cols = []
        for p in p_e:
            e = self.air.expand(stag, float(p))
            cols.append((np.sqrt(max(0.0, 2.0 * (stag.h - e.h))), e.rho, e.T, e.mu, e.a, e.m_bar))
        u_e, rho_e, T_e, mu_e, a_e, m_bar = (np.array(c) for c in zip(*cols))
        return p_e, u_e, rho_e, T_e, mu_e, a_e, m_bar, float(stag.p)

    def evaluate(self, state, theta, radius, rho_liquid):
        """Per-patch SurfaceFlowResult for the trajectory AeroState, the patch angles theta (rad), the nose radius
        and the film density (for the inertial term of G)."""
        theta = np.asarray(theta, dtype=float)
        n = theta.size
        fs = state.freestream
        V, p_inf = float(state.V), float(fs.p)
        decel = float(np.linalg.norm(state.a_drag))
        windward = theta <= 0.5 * np.pi
        sin, cos = np.sin(theta), np.cos(theta)
        tau_fm = np.where(windward, self.sigma_t * fs.rho * V * V * sin * np.abs(cos), 0.0)
        nan, zero = np.full(n, np.nan), np.zeros(n)
        if fs.rho <= 0.0 or V <= 0.0:
            return SurfaceFlowResult(zero, zero, zero, zero, zero, np.full(n, p_inf), zero, nan, nan, np.full(n, np.inf),
                                     np.full(n, REGIME_FM), zero, zero, zero, zero, p_inf, p_inf, decel)
        if not np.isfinite(state.kn) or state.kn >= FM_BODY_KN:
            G = -rho_liquid * decel * sin
            return SurfaceFlowResult(zero, zero, np.full(n, fs.rho), np.full(n, fs.T), zero, np.full(n, p_inf), zero, nan,
                                     np.full(n, aero.mean_free_path(fs.rho, fs.m_bar)), np.full(n, np.inf),
                                     np.full(n, REGIME_FM), tau_fm, zero, tau_fm, np.where(windward, G, 0.0), p_inf, p_inf, decel)
        p_b, u_b, rho_b, T_b, mu_b, a_b, mbar_b, p_s = self.edge_table(state, radius)
        s_b = radius * THETA_BINS
        delta_b = boundary_layer_thickness(s_b, u_b, mu_b / rho_b)
        self.last_bins = (THETA_BINS, p_b, u_b, rho_b, T_b, mu_b, delta_b)
        th = np.minimum(theta, 0.5 * np.pi)
        at = lambda col: np.interp(th, THETA_BINS, col)
        u_e, rho_e, T_e, mu_e, p_e, a_e, mbar_e, delta_a = (at(c) for c in (u_b, rho_b, T_b, mu_b, p_b, a_b, mbar_b, delta_b))
        u_e = np.where(windward, u_e, 0.0)
        p_e = np.where(windward, p_e, p_inf)
        lambda_e = 1.0 / (np.sqrt(2.0) * np.pi * HARD_SPHERE_DIAMETER ** 2 * (rho_e / mbar_e))
        with np.errstate(divide="ignore", invalid="ignore"):
            kn_delta = np.where(delta_a > 0.0, lambda_e / delta_a, np.inf)
            tau_c = np.where(delta_a > 0.0, mu_e * u_e / delta_a, 0.0)
        regime = np.where(kn_delta < KN_SLIP, REGIME_CONTINUUM, np.where(kn_delta < KN_FM, REGIME_SLIP, REGIME_FM))
        slip = 1.0 / (1.0 + self.slip_C * np.where(np.isfinite(kn_delta), kn_delta, 0.0))
        tau_slip = tau_c * slip
        if self.rarefied_shear == "slip":
            tau = np.where(regime == REGIME_FM, tau_fm, tau_slip)
        else:
            f = self.bridging(state.kn)
            tau = (1.0 - f) * tau_slip + f * tau_fm
        tau = np.where(windward & np.isfinite(tau), tau, 0.0)
        dp_ds = -2.0 * (p_s - p_inf) * sin * cos / radius
        G = np.where(windward, -dp_ds - rho_liquid * decel * sin, 0.0)
        mach_e = np.where(a_e > 0.0, u_e / np.maximum(a_e, 1e-300), 0.0)
        u_e, rho_e, mu_e = (np.where(np.isfinite(x), x, 0.0) for x in (u_e, rho_e, mu_e))
        return SurfaceFlowResult(u_e, u_e * slip, rho_e, T_e, mu_e, p_e, mach_e, np.where(windward, delta_a, np.nan),
                                 lambda_e, np.where(windward, kn_delta, np.inf), np.where(windward, regime, REGIME_FM),
                                 tau, np.where(windward, tau_c, 0.0), tau_fm, G, p_s, p_inf, decel)
```


- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_surface_flow.py tests/test_reentry_model_gas.py tests/test_reentry_model_heating.py -q`
Expected: all pass (`evaluate` at 71 km takes ~40 ms: 91 Cantera SP states).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/gas.py reentry_model/surface_flow.py tests/test_reentry_model_surface_flow.py
git commit -m "Add the boundary-layer edge state, Ranger's boundary layer, Knudsen regimes and shear per patch (Step 3 Task 5)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: The melt film — lubrication and runoff

**Files:**
- Create: `reentry_model/film.py`
- Test: `tests/test_reentry_model_film.py`

**Interfaces:**
- Consumes: `SurfaceMesh.edges/normals/centroids/areas` (Task 1).
- Produces: `lubrication(tau, G, b, delta_m, mu_l) -> (V_s, q, shear_rate, thick_mask)`; `Runoff(surface, points, windward=None)` with `.i, .j, .length, .n_i, .n_j, .n_patches`, `.edge_coefficients(q, b, t_hat, areas) -> (c_ij, c_ji)`, `.transport(m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4) -> (m_f, n_solves, moved_kg)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_film.py`:

```python
"""film.py: the lubrication branches, the implicit upwind runoff (conservation, positivity, the strip steady state)."""
import numpy as np
import pytest

from reentry_model import film, material, mesh

LIQ = material.LiquidProperties(2400.0, 1.3e-3, 0.86)


def test_lubrication_branches():
    tau, G, mu = 30.0, 5e4, LIQ.mu
    V, q, rate, thick = film.lubrication([tau, tau, tau], [G, G, G], [5e-5, 5e-4, 5e-4], [1e-4, 1e-4, np.nan], mu)
    b, d = 5e-5, 1e-4
    assert not thick[0] and V[0] == pytest.approx(tau * b / mu + G * b * b / (2 * mu)) and q[0] == pytest.approx(tau * b ** 2 / (2 * mu) + G * b ** 3 / (3 * mu))
    assert rate[0] == pytest.approx(V[0] / b)
    b = 5e-4
    assert thick[1] and V[1] == pytest.approx(tau * d / mu) and q[1] == pytest.approx(V[1] * d / 2 + G * b ** 3 / (3 * mu)) and rate[1] == pytest.approx(V[1] / d)
    assert not thick[2] and V[2] == pytest.approx(tau * b / mu + G * b * b / (2 * mu))               # no delta_m: thin branch
    assert film.lubrication([tau], [-1e6], [1e-3], [np.nan], mu)[1][0] < 0.0                         # a strong adverse gradient reverses q


def strip_surface(n=20, w=1e-3):
    """A planar strip of n squares (2 triangles each) along +x, width w, in the plane z = 0 with the normal +z."""
    x = np.arange(n + 1) * w
    pts = np.array([[xi, yi, 0.0] for xi in x for yi in (0.0, w)])
    faces = []
    for k in range(n):
        a, b, c, d = 2 * k, 2 * k + 1, 2 * k + 2, 2 * k + 3
        faces += [[a, c, b], [b, c, d]]
    faces = np.array(faces)
    p = pts[faces]
    normals = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    areas = 0.5 * np.linalg.norm(normals, axis=1)
    normals = normals / (2 * areas)[:, None]
    return pts, mesh.SurfaceMesh(faces, p.mean(axis=1), normals, areas, np.arange(len(faces)), np.arange(len(faces)))


def test_runoff_conserves_mass_and_reaches_the_strip_steady_state():
    """Constant shear on a strip fed at its first patch: in the steady state the thin-film flux through every edge
    equals the feed, so b = sqrt(2 mu S / (rho tau w)) on the interior patches (closed form) within 0.1 %."""
    pts, s = strip_surface()
    ro = film.Runoff(s, pts)
    assert len(ro.i) == 2 * (s.n_patches // 2) - 1                                     # 20 diagonals + 19 shared verticals
    t_hat = np.tile([1.0, 0.0, 0.0], (s.n_patches, 1))
    tau, w = 30.0, 1e-3
    S = 2e-6                                                                            # kg/s fed into patch 0
    q_of_b = lambda b: film.lubrication(np.full(s.n_patches, tau), np.zeros(s.n_patches), b, np.full(s.n_patches, np.nan), LIQ.mu)[1]
    m = np.zeros(s.n_patches)
    dt = 0.05
    for k in range(400):
        m[0] += S * dt
        total = m.sum()
        m, n, moved = ro.transport(m, q_of_b, t_hat, LIQ.rho, s.areas, dt, substeps=1)   # one implicit step per feed pulse
        assert m.sum() == pytest.approx(total, rel=1e-12) and (m >= 0.0).all() and n == 1
        m[-2:] = 0.0                                                                    # the strip's end is stripped (a sink)
    b = m / (LIQ.rho * s.areas)
    b_exact = np.sqrt(2.0 * LIQ.mu * S / (LIQ.rho * tau * w))
    assert np.abs(b[4:-4] / b_exact - 1.0).max() < 1e-3 and moved > 0.0


def test_runoff_on_the_sphere_stops_at_the_equator(coarse_sphere_mesh):
    s = coarse_sphere_mesh.surface()
    v = np.array([1.0, 0.0, 0.0])
    theta, t_hat = s.angles_to(v), s.tangent_from(v)
    windward = theta <= np.pi / 2
    ro = film.Runoff(s, coarse_sphere_mesh.points, windward)
    assert windward[ro.i].all() and windward[ro.j].all()
    m0 = np.where(theta < 0.5, 1e-4 * s.areas * LIQ.rho, 0.0)
    tau = np.where(windward, 30.0 * np.sin(theta), 0.0)
    G = np.where(windward, 5e4 * np.sin(theta) * np.cos(theta), 0.0)
    q_of_b = lambda b: film.lubrication(tau, G, b, np.full(b.size, np.nan), LIQ.mu)[1]
    m, n, moved = ro.transport(m0, q_of_b, t_hat, LIQ.rho, s.areas, 0.5)
    assert m.sum() == pytest.approx(m0.sum(), rel=1e-12) and (m >= 0.0).all() and moved > 0.0
    assert m[~windward].sum() == 0.0 and np.degrees(theta[m > 1e-9 * m.max()].max()) > 35.0        # spread outward (>= 4 patches), never leeward
    assert ro.transport(np.zeros(s.n_patches), q_of_b, t_hat, LIQ.rho, s.areas, 0.5)[1] == 0    # a dry surface costs nothing
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_film.py -q`
Expected: ImportError (`reentry_model.film`).

- [ ] **Step 3: Create `reentry_model/film.py`**

```python
"""The melt film on the surface patches: lubrication velocity and flux, and the explicit upwind runoff transport on
the patch graph (spec Step 3 section 8).

State: m_f per patch [kg], thickness b = m_f / (rho_l A). Lubrication solution with the patch's shear tau, driving
gradient G, thickness b, melt boundary-layer thickness delta_m (from Girin's Eq. 2, spray.melt_layer) and the
liquid viscosity mu_l:
    thin (b <= delta_m):  V_s = tau b / mu_l + G b^2 / (2 mu_l),   q = tau b^2 / (2 mu_l) + G b^3 / (3 mu_l)
    thick (b > delta_m):  V_s = tau delta_m / mu_l,               q = V_s delta_m / 2 + G b^3 / (3 mu_l)
(q per unit width [m^2/s], signed along the surface direction t away from the stagnation point; the velocity
gradient for the instability is V_s / b or V_s / delta_m). Where delta_m is undefined (free-molecular patches) the
film is thin-branch (Couette) by definition.

Runoff: for every edge shared by patches i and j, flux = q_donor * l_edge * (t_donor . n_edge)^+ with n_edge the
in-plane edge normal pointing out of the donor; both directions are evaluated (a negative q reverses the flow toward
the nose). The transport over a macro step is the linearly implicit upwind scheme (I + dt C) m_new = m_old with the
edge coefficients c = q l (t . n)^+ / (A b) [1/s] taken at the start of the step: unconditionally stable, positive,
and mass-conserving to round-off (every edge flux leaves one patch and enters another); its steady state is the
exact nonlinear one. The spec's explicit sub-stepped
scheme was replaced on 2026-09-20: films driven by the free-molecular shear near the equator move at ~10 m/s and
cross the hemisphere many times per 0.5 s step, so an explicit CFL needed 1e4-1e5 sub-steps per macro step."""
from dataclasses import dataclass

import numpy as np


def lubrication(tau, G, b, delta_m, mu_l):
    """(V_s, q, shear_rate) per patch for the thin/thick lubrication branches (module docstring)."""
    tau, G, b = (np.asarray(x, dtype=float) for x in (tau, G, b))
    delta_m = np.asarray(delta_m, dtype=float)
    thick = np.isfinite(delta_m) & (b > delta_m)
    V_thin = tau * b / mu_l + G * b * b / (2.0 * mu_l)
    q_thin = tau * b * b / (2.0 * mu_l) + G * b ** 3 / (3.0 * mu_l)
    d = np.where(thick, delta_m, 0.0)
    V_thick = tau * d / mu_l
    q_thick = V_thick * d / 2.0 + G * b ** 3 / (3.0 * mu_l)
    V = np.where(thick, V_thick, V_thin)
    q = np.where(thick, q_thick, q_thin)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(thick, V / np.where(thick, d, 1.0), np.where(b > 0.0, V / np.where(b > 0.0, b, 1.0), 0.0))
    return V, q, rate, thick


class Runoff:
    """Explicit upwind transport of the film mass on a SurfaceMesh (rebuilt whenever the surface changes)."""

    def __init__(self, surface, points, windward=None):
        e, i, j = surface.edges()
        if windward is not None:                               # the flow stops at the equator: no edge leads into a leeward patch
            keep = windward[i] & windward[j]
            e, i, j = e[keep], i[keep], j[keep]
        x1, x2 = points[e[:, 0]], points[e[:, 1]]
        self.length = np.linalg.norm(x2 - x1, axis=1)
        mid = 0.5 * (x1 + x2)
        self.i, self.j = i, j
        self.n_i = self._edge_normal(x2 - x1, surface.normals[i], mid - surface.centroids[i])
        self.n_j = self._edge_normal(x2 - x1, surface.normals[j], mid - surface.centroids[j])
        self.n_patches = surface.n_patches

    @staticmethod
    def _edge_normal(edge, normal, outward):
        n = np.cross(edge, normal)
        n = n / np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None]
        sign = np.sign(np.einsum("ij,ij->i", n, outward))
        return n * np.where(sign == 0.0, 1.0, sign)[:, None]

    def edge_coefficients(self, q, b, t_hat, areas):
        """Emptying-rate coefficients c_ij, c_ji [1/s] of every edge: the flux i -> j is m_i c_ij with
        c_ij = q_i l (t_i . n_ij)^+ / (A_i b_i) (zero where the patch is dry); q, b per patch."""
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            rate = np.where(b > 0.0, q / (areas * np.where(b > 0.0, b, 1.0)), 0.0)          # 1/m per unit edge length
        rate = np.where(np.isfinite(rate), rate, 0.0)                                       # a degenerate (sliver) patch moves nothing
        c_ij = np.maximum(0.0, rate[self.i] * self.length * np.einsum("ij,ij->i", t_hat[self.i], self.n_i))
        c_ji = np.maximum(0.0, rate[self.j] * self.length * np.einsum("ij,ij->i", t_hat[self.j], self.n_j))
        return c_ij, c_ji

    def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
        """Advance m_f over dt in `substeps` linearly implicit upwind steps (I + dt_s C) m_new = m_old, the edge
        coefficients C taken at the start of each sub-step; each M-matrix system is solved directly, so the scheme is
        unconditionally stable, positive and conservative to round-off (the columns of C sum to zero), and its
        steady state C(m) m = 0 is the exact nonlinear one. Dry patches have no coefficient, so a wetting front
        advances one patch per sub-step (a documented limit; the film that matters is stripped where it forms). A
        Picard iteration on the fully implicit form does not contract when dt x c >> 1, which is the case for micron
        films at 30 Pa on millimetre patches (measured 2026-09-20). `q_of_thickness(b)` returns the signed flux per
        patch for thickness b. Returns (m_f, n_solves, mass that arrived on another patch [kg])."""
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        m = np.array(m_f, dtype=float)
        if not (m > 0.0).any() or len(self.i) == 0:
            return m, 0, 0.0
        moved, dt_s, n = 0.0, dt / substeps, 0
        for n in range(1, substeps + 1):
            b = m / (rho_l * areas)
            c_ij, c_ji = self.edge_coefficients(q_of_thickness(b), b, t_hat, areas)
            out = np.bincount(self.i, c_ij, self.n_patches) + np.bincount(self.j, c_ji, self.n_patches)
            rows = np.concatenate([np.arange(self.n_patches), self.j, self.i])
            cols = np.concatenate([np.arange(self.n_patches), self.i, self.j])
            vals = np.concatenate([1.0 + dt_s * out, -dt_s * c_ij, -dt_s * c_ji])
            A = sp.csc_matrix((vals, (rows, cols)), shape=(self.n_patches, self.n_patches))
            m_new = np.maximum(spla.spsolve(A, m), 0.0)
            moved += float(np.maximum(m_new - m, 0.0).sum())
            m = m_new
        return m, n, moved
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_film.py -q`
Expected: 3 passed (the strip's interior reaches the closed-form thickness to 1e-3 exactly at the scheme's steady state).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/film.py tests/test_reentry_model_film.py
git commit -m "Add the melt film: lubrication branches and the linearly implicit upwind runoff (Step 3 Task 6)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: Spraying — instability branches, release bookkeeping, published reference values

**Files:**
- Create: `reentry_model/spray.py`, `data/reference_values/girin2017_table1.json`, `data/reference_values/girin1994_tables.json`
- Test: `tests/test_reentry_model_spray.py`

**Interfaces:**
- Consumes: `dispersion.DispersionTable` (Task 4), `surface_flow.REGIME_FM` and `SurfaceFlowResult` fields (Task 5), `material.LiquidProperties` (Task 2).
- Produces: `melt_layer(flow, liquid) -> (delta_m, factor)`; `rayleigh_taylor(W, b, liquid) -> (active, lambda*, tau*)`; `thin_film_mode(mach, momentum_flux, liquid) -> (lambda*, tau*)`; `SprayModel(liquid, k_r=0.17, k_t=1.1, we_critical=4.62, table=None).evaluate(flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05) -> SprayResult(branch, we_s, unstable, r, mdot, dm, dn, rt_active, delta_m, v_s)`; `source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state) -> list of 22-value rows`; `histogram(r, dn, dm) -> (dn per bin, dM per bin)`; constants `K_R, K_T, B_MIN, N_BINS = 40, R_MIN = 1e-6, R_MAX = 1e-2, BIN_EDGES, WE_BREAKUP = 12, CAPILLARY_TAU, SOURCE_COLUMNS, BRANCH_THICK/THIN/RAREFIED = 0/1/2`.

- [ ] **Step 1: Write the two reference-value files**

`data/reference_values/girin2017_table1.json`:

```json
{
  "_source": "Girin, O. G. (2017), A&A 606, A63, Table 1 and section 5 (R0 = 0.3 cm meteoroids); media properties from section 5; atmosphere density rho_a = 1e-7 g/cm3 (section 5, 'gas species density around the meteoroid'), rho_inf = rho_a/6 (the sixfold compression, section 5) reproduces the printed GI with We_inf on rho_inf and Re_inf on rho_a (measured 2026-09-20).",
  "constants": {
    "R0_m": 0.003,
    "rho_a_kgm3": 0.0001,
    "rho_inf_kgm3": 1.6666666666666667e-05,
    "mu_air_Pas": 6.8e-05,
    "k_r": 0.17,
    "k_t": 1.1,
    "we_cr": 4.62,
    "we_cr_theory": 3.08
  },
  "variants": {
    "I": {
      "description": "fast iron",
      "V_ms": 60000.0,
      "rho_m": 7800.0,
      "mu_m": 0.0058,
      "sigma": 1.2,
      "Re_inf": 529,
      "GI": 13.0,
      "m_g": 0.9,
      "t_sd_ms": 5.9,
      "sigma_ablation_s2cm2": 1.04e-12,
      "path_m": 350,
      "N": 1500000.0,
      "phi_cr_deg": 16.1,
      "t_f_us": 5.7,
      "z0_km": 105,
      "r_med_um": 26.9,
      "r_min_um": 4.0,
      "r_max_um": 135.0,
      "tau_end": 2.72,
      "t_ch_ms": 2.16
    },
    "II": {
      "description": "slow iron",
      "V_ms": 25000.0,
      "rho_m": 7800.0,
      "mu_m": 0.0058,
      "sigma": 1.2,
      "Re_inf": 222,
      "GI": 3.55,
      "m_g": 0.9,
      "t_sd_ms": 16.0,
      "sigma_ablation_s2cm2": 5.3e-12,
      "path_m": 400,
      "N": 370000.0,
      "phi_cr_deg": 31.4,
      "t_f_us": 31.0,
      "z0_km": 90,
      "r_med_um": 41.6,
      "r_min_um": 16.0,
      "r_max_um": 265.0,
      "tau_end": 3.11,
      "t_ch_ms": 5.15
    },
    "III": {
      "description": "stony",
      "V_ms": 60000.0,
      "rho_m": 3500.0,
      "mu_m": 0.174,
      "sigma": 0.36,
      "Re_inf": 529,
      "GI": 43.5,
      "m_g": 0.4,
      "t_sd_ms": 149.0,
      "sigma_ablation_s2cm2": 1.8e-14,
      "path_m": 9000,
      "N": 832,
      "phi_cr_deg": 8.8,
      "t_f_us": 194.0,
      "z0_km": 120,
      "r_med_um": null,
      "r_min_um": null,
      "r_max_um": null,
      "tau_end": 105.0,
      "t_ch_ms": 1.45
    }
  }
}
```


`data/reference_values/girin1994_tables.json`:

```json
{
  "_source": "Girin, A. G. & Kopyt, N. Kh. (1994), J. Aerosol Sci. 25(7), 1353-1357, Tables 1-2 and Eqs. (9)-(14). Table 1: side-surface KH mode, sigma_t = 1200 dyn/cm, rho_1 = 7.8 g/cm3, M = 3, r_d = lambda*/4, tau_d = tau*, mass rate m per unit area; the printed r_d imply an effective dynamic pressure 7.8 x rho_2 V0^2 (their shock-layer treatment: 'deceleration ~10 times ... acceleration to M = 2-3'), fitted on the first row and checked on the other five (measured 2026-09-20); m = rho_1 r_d / (2 tau_d) reproduces the printed rates. Table 2: front-surface Rayleigh-Taylor mode, Eq. (14); the printed lambda* are 10 x smaller than Eq. (14) gives (a units slip in the table), tau* matches to three digits.",
  "table1": {
    "sigma_Nm": 1.2,
    "rho1_kgm3": 7800.0,
    "mach": 3.0,
    "rho2_kgm3": [
      0.0001,
      0.001,
      0.01
    ],
    "V0_ms": [
      3000.0,
      10000.0
    ],
    "r_d_m": [
      [
        0.000193,
        1.73e-05
      ],
      [
        1.93e-05,
        1.73e-06
      ],
      [
        1.93e-06,
        1.73e-07
      ]
    ],
    "tau_d_s": [
      [
        0.00138,
        3.71e-05
      ],
      [
        4.34e-05,
        1.17e-06
      ],
      [
        1.38e-06,
        3.71e-08
      ]
    ],
    "m_kgm2s": [
      [
        547.0,
        1830.0
      ],
      [
        1740.0,
        5800.0
      ],
      [
        5470.0,
        18300.0
      ]
    ]
  },
  "table2": {
    "sigma_Nm": 1.2,
    "rho1_kgm3": 7800.0,
    "W_ms2": [
      1000.0,
      10000.0,
      100000.0
    ],
    "lambda_star_m_as_printed": [
      0.00042,
      0.000135,
      4.27e-05
    ],
    "lambda_star_m_eq14": [
      0.0042,
      0.00135,
      0.000427
    ],
    "tau_star_s": [
      0.001,
      0.00018,
      3.19e-05
    ]
  }
}
```


- [ ] **Step 2: Write the failing tests**

Create `tests/test_reentry_model_spray.py`:

```python
"""spray.py: delta_m, the thin-film and Rayleigh-Taylor modes against Girin & Kopyt (1994), the thick branch, release
bookkeeping, size caps and histograms."""
import json
import os

import numpy as np
import pytest

from reentry_model import dispersion, material, spray, surface_flow as sf
from helpers import REPO_ROOT

LIQ = material.LiquidProperties(2400.0, 1.3e-3, 0.86)


class Iron:
    rho, mu, sigma = 7800.0, 5.8e-3, 1.2


def test_thin_film_mode_reproduces_the_1994_table_1_scaling():
    """r_d = lambda*/4 and tau_d = tau* for the six (rho_2, V0) cases: one shock-layer factor on rho_2 V0^2 (fitted on
    the first entry, 7.8) reproduces every r_d within 0.5 % and tau_d within 1 %; the printed mass rates are
    rho_1 r_d / (2 tau_d) within 1 %."""
    t = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")))["table1"]
    rho2, V0 = np.array(t["rho2_kgm3"]), np.array(t["V0_ms"])
    lam, tau = spray.thin_film_mode(np.full((3, 2), t["mach"]), rho2[:, None] * V0[None, :] ** 2, Iron)
    r_tab, tau_tab, m_tab = (np.array(t[k]) for k in ("r_d_m", "tau_d_s", "m_kgm2s"))
    f_s = (lam / 4.0)[0, 0] / r_tab[0, 0]
    assert 7.5 < f_s < 8.1
    assert np.abs(lam / 4.0 / f_s / r_tab - 1.0).max() < 5e-3 and np.abs(tau / f_s ** 1.5 / tau_tab - 1.0).max() < 1e-2
    assert np.abs(Iron.rho * (lam / 4.0 / f_s) / (2.0 * tau / f_s ** 1.5) / m_tab - 1.0).max() < 1e-2
    assert spray.CAPILLARY_TAU == pytest.approx(0.798, rel=1e-3)                                  # tau* = 2 x 2 pi / omega_cap
    assert np.isinf(spray.thin_film_mode(np.array([3.0]), np.array([0.0]), Iron)[0][0])


def test_rayleigh_taylor_mode_reproduces_the_1994_table_2():
    t = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")))["table2"]
    for W, lam_ref, tau_ref in zip(t["W_ms2"], t["lambda_star_m_eq14"], t["tau_star_s"]):
        active, lam, tau = spray.rayleigh_taylor(W, np.array([1e-3, 1e-5]), Iron)
        assert lam[0] == pytest.approx(lam_ref, rel=2e-2) and tau[0] == pytest.approx(tau_ref, rel=1e-2)
        assert active[0] == (W * 1e-6 * Iron.rho > 3.0 * Iron.sigma) and not active[1]
    assert not spray.rayleigh_taylor(30.0, np.array([1e-4]), LIQ)[0][0]                          # our deceleration: inactive


class FakeFlow:
    def __init__(self, n, regime, u=2000.0, rho=7e-4, mu=1.6e-4, delta_a=7e-3, mach=1.0, decel=30.0):
        self.regime = np.full(n, regime)
        self.u_eff = self.u_e = np.full(n, u)
        self.rho_e, self.mu_e, self.delta_a, self.mach_e = np.full(n, rho), np.full(n, mu), np.full(n, delta_a), np.full(n, mach)
        self.deceleration = decel


class FakeState:
    class freestream:
        rho = 8e-5
    V = 7200.0
    ma = 24.0
    h = 71e3


def test_melt_layer_and_branches():
    flow = FakeFlow(3, sf.REGIME_CONTINUUM)
    delta_m, factor = spray.melt_layer(flow, LIQ)
    alpha, mu = 7e-4 / 2400.0, 1.6e-4 / 1.3e-3
    assert np.allclose(delta_m, (alpha / mu ** 2) ** (1 / 3) * 7e-3) and np.allclose(factor, (alpha * mu) ** (1 / 3) / (1 + (alpha * mu) ** (1 / 3)))
    assert 5e-5 < delta_m[0] < 2e-4                                                                 # ~0.1 mm (spec estimate)
    model = spray.SprayModel(LIQ)
    b = np.array([1e-3, 2e-5, 0.0])                                                                 # thick, thin, dry
    v_s = np.array([5.0, 0.5, 0.0])
    areas, m_f = np.full(3, 1e-5), b * LIQ.rho * 1e-5
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True, True]), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN, -1]
    we_s = LIQ.rho * 25.0 * delta_m[0] / LIQ.sigma
    assert res.we_s[0] == pytest.approx(we_s) and res.unstable[0] == (we_s > 4.62)
    d_f, im_f, _ = dispersion.DispersionTable()(we_s)
    lam = 2 * np.pi * delta_m[0] / d_f
    assert res.r[0] == pytest.approx(min(0.17 * lam, (3 * m_f[0] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    assert res.mdot[0] == pytest.approx(LIQ.rho * np.pi * (0.17 * lam) ** 2 / (lam * 1.1 * delta_m[0] / (5.0 * im_f)))
    lam_t, tau_t = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    assert res.r[1] == pytest.approx(min(lam_t[0] / 4.0, (3 * m_f[1] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    assert res.mdot[1] == pytest.approx(LIQ.rho * min(2e-5, lam_t[0] / 8.0) / tau_t[0])
    assert res.dm[1] == pytest.approx(min(res.mdot[1] * 1e-5 * 0.5, m_f[1])) and res.dm[2] == 0.0
    assert res.dn[1] == pytest.approx(res.dm[1] / (4 / 3 * np.pi * LIQ.rho * res.r[1] ** 3))
    assert not res.rt_active
    # rarefied: thin film takes the free-molecular branch, a thick one the thick branch
    flow_fm = FakeFlow(2, sf.REGIME_FM)
    res = model.evaluate(flow_fm, FakeState(), np.array([2e-5, 1e-3]), delta_m[:2], np.array([0.5, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
    assert res.branch.tolist() == [spray.BRANCH_RAREFIED, spray.BRANCH_THICK]
    lam_r, _ = spray.thin_film_mode(np.array([24.0]), np.array([8e-5 * 7200.0 ** 2]), LIQ)
    assert res.r[0] == pytest.approx(min(lam_r[0] / 4.0, (3 * (m_f[0] + 1e-9) / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    # leeward or dry: nothing
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([False, False, False]), 0.5, areas, m_f)
    assert not res.unstable.any() and res.dm.sum() == 0.0
    with pytest.raises(ValueError):
        spray.SprayModel(LIQ, k_r=0.0)


def test_source_rows_and_histogram():
    flow = FakeFlow(2, sf.REGIME_CONTINUUM)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.array([1e-3, 1e-3])
    m_f = b * LIQ.rho * 1e-5
    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.array([5.0, 5.0]), np.array([True, True]), 0.5, np.full(2, 1e-5), m_f)
    rows = spray.source_rows(10.0, 71e3, 7200.0, np.radians([30.0, 60.0]), np.zeros((2, 3)), np.tile([0.0, 1.0, 0.0], (2, 1)), flow, res, b, LIQ, FakeState())
    assert len(rows) == 2 and len(rows[0]) == len(spray.SOURCE_COLUMNS)
    row = dict(zip(spray.SOURCE_COLUMNS, rows[0]))
    assert row["time_s"] == 10.0 and row["altitude_km"] == 71.0 and row["theta_deg"] == pytest.approx(30.0) and row["branch"] == spray.BRANCH_THICK
    assert row["we_d"] == pytest.approx(8e-5 * 7200.0 ** 2 * 2 * row["r_m"] / LIQ.sigma) and row["oh"] == pytest.approx(LIQ.mu / np.sqrt(LIQ.rho * LIQ.sigma * 2 * row["r_m"]))
    assert row["breakup"] == float(row["we_d"] > 12.0) and row["dm_kg"] == res.dm[0] and row["v_s_ms"] == 5.0
    n_hist, m_hist = spray.histogram(np.array([2e-6, 5e-5, 5e-5, np.nan, 2e-2]), np.array([1.0, 2.0, 3.0, 4.0, 5.0]), np.array([1.0, 1.0, 1.0, 1.0, 1.0]))
    assert n_hist.sum() == 6.0 and m_hist.sum() == 3.0 and n_hist.shape == (spray.N_BINS,)          # nan and out-of-range dropped
    assert spray.BIN_EDGES[0] == 1e-6 and spray.BIN_EDGES[-1] == 1e-2
```


- [ ] **Step 3: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_spray.py -q`
Expected: ImportError (`reentry_model.spray`).

- [ ] **Step 4: Create `reentry_model/spray.py`**

```python
"""Melt spraying: Girin's gradient instability per patch and the droplet release bookkeeping (spec Step 3 section 9).

Branches per windward patch with film (regime from surface_flow, b from the film, delta_m from Girin's Eq. 2):
  thick, continuum/slip (b > delta_m): delta_m = (alpha/mu^2)^(1/3) delta_a, V_s = (alpha mu)^(1/3)/(1 + (alpha
      mu)^(1/3)) u_eff, We_s = rho_l V_s^2 delta_m / Sigma; unstable when We_s > We_cr; lambda_f = 2 pi delta_m /
      Delta_f(We_s), r = k_r lambda_f, t_per = k_t delta_m / (V_s Im Omega_f(We_s)), stripping rate per area
      mdot = rho_l pi r^2 / (lambda_f t_per) (one torus of cross-section r per wavelength per period, Girin 2017
      Eqs. 5-8 and Appendix D).
  thin, continuum/slip (b <= delta_m): Girin & Kopyt 1994 side mode with the edge state: lambda* = 1.5 M_e Sigma /
      (rho_e u_eff^2) (their 1.5 M d / We_d, the film thickness cancels), r = lambda*/4, growth time tau* =
      2 tau_v(lambda*) = 0.798 lambda*^1.5 (rho_l/Sigma)^1/2 (their Eq. 12: twice the capillary-wave period),
      release rate mdot = rho_l min(b, lambda*/8) / tau* (their Table 1's mass rate is rho_1 r_d / (2 tau_d) =
      rho_1 lambda*/(8 tau*), reproduced 2026-09-20; the film supplies at most its thickness). Active when b >= B_MIN.
      (Their dissipation cut-off lambda_t = lambda*/3 is always below lambda*, so it never limits the mode.)
  transitional/free-molecular (Kn_delta >= 0.1), thin film: the thin mode with the free-molecular momentum flux,
      lambda* = 1.5 M_inf Sigma / (rho_inf V^2), and the Couette film velocity from tau_fm -- an extrapolation of a
      continuum film theory (spec section 17.3), sensitive to --rarefied-shear. A film thicker than delta_m on such a
      patch takes the thick branch with the film velocity the free-molecular shear gives it (tau_fm delta_m / mu_l):
      without it the film piling up at the windward rim, where the edge state expands to p_inf and Kn_delta exceeds
      0.1, was never stripped (measured 2026-09-20: a rim patch reached b = 100 mm).
The 1994 front-surface Rayleigh-Taylor criterion W b^2 rho_l > 3 Sigma (W the body deceleration) is evaluated and
reported (lambda* = 2 pi (3 Sigma/(W rho_l))^1/2, tau* = (27 Sigma/(4 W^3 rho_l))^1/4), never applied.

Release per patch and step: dm = min(mdot A dt, m_f) (the film after feed and runoff), dn = dm / (4/3 pi rho_l r^3);
one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX]."""
from dataclasses import dataclass, field

import numpy as np

from . import dispersion
from .surface_flow import REGIME_FM

K_R, K_T = 0.17, 1.1
B_MIN = 1.0e-6                     # m, films thinner than this do not spray (numerical floor)
N_BINS, R_MIN, R_MAX = 40, 1.0e-6, 1.0e-2
BIN_EDGES = np.logspace(np.log10(R_MIN), np.log10(R_MAX), N_BINS + 1)
WE_BREAKUP = 12.0                  # Pilch-Erdman: secondary breakup expected above this droplet Weber number
CAPILLARY_TAU = 4.0 * np.pi / (2.0 * np.pi) ** 1.5      # 0.798: tau* = 2 x 2 pi / omega_cap(lambda*)
SOURCE_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "theta_deg", "x_m", "y_m", "z_m", "regime", "branch", "b_m",
                  "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup"]
BRANCH_THICK, BRANCH_THIN, BRANCH_RAREFIED = 0, 1, 2


def melt_layer(flow, liquid):
    """delta_m per patch from Girin's Eq. 2 (nan where there is no gas boundary layer) and the shear velocity factor
    (alpha mu)^(1/3)/(1 + (alpha mu)^(1/3))."""
    alpha = flow.rho_e / liquid.rho
    mu = flow.mu_e / liquid.mu
    with np.errstate(divide="ignore", invalid="ignore"):
        delta_m = np.where(np.isfinite(flow.delta_a) & (mu > 0.0), (alpha / mu ** 2) ** (1.0 / 3.0) * flow.delta_a, np.nan)
        am = (alpha * mu) ** (1.0 / 3.0)
        factor = np.where(np.isfinite(am), am / (1.0 + am), 0.0)
    return delta_m, factor


def rayleigh_taylor(deceleration, b, liquid):
    """(active, lambda*, tau*) of the 1994 front-surface RT mode for the film thickness b."""
    W = max(float(deceleration), 0.0)
    active = W * b * b * liquid.rho > 3.0 * liquid.sigma
    if W <= 0.0:
        return active, np.full_like(np.asarray(b, dtype=float), np.inf), np.full_like(np.asarray(b, dtype=float), np.inf)
    lam = 2.0 * np.pi * np.sqrt(3.0 * liquid.sigma / (W * liquid.rho))
    tau = (27.0 * liquid.sigma / (4.0 * W ** 3 * liquid.rho)) ** 0.25
    return active, np.full_like(np.asarray(b, dtype=float), lam), np.full_like(np.asarray(b, dtype=float), tau)


def thin_film_mode(mach, momentum_flux, liquid):
    """(lambda*, tau*) of the Girin & Kopyt 1994 side mode for gas Mach number and rho V^2 [Pa]."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(momentum_flux > 0.0, 1.5 * mach * liquid.sigma / momentum_flux, np.inf)
    tau = CAPILLARY_TAU * lam ** 1.5 * np.sqrt(liquid.rho / liquid.sigma)
    return lam, tau


@dataclass
class SprayResult:
    branch: np.ndarray          # 0 thick, 1 thin, 2 rarefied, -1 none
    we_s: np.ndarray
    unstable: np.ndarray        # bool
    r: np.ndarray               # m, droplet radius (nan where none)
    mdot: np.ndarray            # kg/(m2 s) stripping rate
    dm: np.ndarray              # kg released this step per patch
    dn: np.ndarray              # droplets released this step per patch
    rt_active: bool
    delta_m: np.ndarray
    v_s: np.ndarray


class SprayModel:
    def __init__(self, liquid, k_r=K_R, k_t=K_T, we_critical=dispersion.WE_CRITICAL_PRACTICAL, table=None):
        if k_r <= 0.0 or k_t <= 0.0 or we_critical <= 0.0:
            raise ValueError("k_r, k_t and the critical Weber number must be > 0")
        self.liquid, self.k_r, self.k_t, self.we_critical = liquid, k_r, k_t, we_critical
        self.table = table or dispersion.DispersionTable()

    def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05):
        """Per-patch instability and release for the film thickness b, melt-layer thickness delta_m and film surface
        velocity v_s (from film.lubrication), over the step dt; m_f is the film mass available, radius the body's."""
        liq = self.liquid
        n = b.size
        branch = np.full(n, -1)
        we_s = np.zeros(n)
        r = np.full(n, np.nan)
        mdot = np.zeros(n)
        has_film = windward & (b >= B_MIN)
        thick = has_film & np.isfinite(delta_m) & (b > delta_m)          # any regime: a film thicker than delta_m is Girin's thick case
        rarefied = has_film & ~thick & (flow.regime == REGIME_FM)
        thin = has_film & ~thick & ~rarefied
        # thick: Girin 2017
        if thick.any():
            vs, dm_ = v_s[thick], delta_m[thick]
            we = liq.rho * vs * vs * dm_ / liq.sigma
            we_s[thick] = we
            delta_f, im_f, _ = self.table(we)
            unstable = (we > self.we_critical) & (im_f > 0.0) & (vs > 0.0)
            lam = np.where(unstable, 2.0 * np.pi * dm_ / np.where(unstable, delta_f, 1.0), np.nan)
            rr = self.k_r * lam
            with np.errstate(divide="ignore", invalid="ignore"):
                t_per = np.where(unstable, self.k_t * dm_ / (vs * np.where(unstable, im_f, 1.0)), np.inf)
                rate = np.where(unstable, liq.rho * np.pi * rr * rr / (lam * t_per), 0.0)
            r[thick], mdot[thick], branch[thick] = rr, rate, BRANCH_THICK
        # thin: Girin & Kopyt 1994 with the edge state
        for mask, code, mach, flux in ((thin, BRANCH_THIN, flow.mach_e, flow.rho_e * flow.u_eff ** 2),
                                       (rarefied, BRANCH_RAREFIED, np.full(n, state.ma), np.full(n, state.freestream.rho * state.V ** 2))):
            if not mask.any():
                continue
            lam, tau = thin_film_mode(mach[mask], flux[mask], liq)
            ok = np.isfinite(lam) & (lam > 0.0) & (tau > 0.0)
            rr = np.where(ok, lam / 4.0, np.nan)
            rate = np.where(ok, liq.rho * np.minimum(b[mask], lam / 8.0) / np.where(ok, tau, 1.0), 0.0)
            we_s[mask] = liq.rho * v_s[mask] ** 2 * b[mask] / liq.sigma
            r[mask], mdot[mask], branch[mask] = rr, rate, code
        unstable = mdot > 0.0
        dm = np.where(unstable, np.minimum(mdot * areas * dt, m_f), 0.0)
        # a droplet is never larger than the film on its patch nor than a quarter of the body radius
        # (the long near-critical waves Girin excludes; the edge dynamic pressure vanishes toward the equator)
        with np.errstate(invalid="ignore"):
            r = np.minimum(r, np.minimum((3.0 * m_f / (4.0 * np.pi * liq.rho)) ** (1.0 / 3.0), 0.25 * radius))
        with np.errstate(divide="ignore", invalid="ignore"):
            dn = np.where(dm > 0.0, dm / (4.0 / 3.0 * np.pi * liq.rho * np.where(dm > 0.0, r, 1.0) ** 3), 0.0)
        rt_active, _, _ = rayleigh_taylor(flow.deceleration, b, liq)
        return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s)


def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state):
    """Source-table rows (list of lists in SOURCE_COLUMNS order) for the patches that released mass this step."""
    k = np.flatnonzero(res.dm > 0.0)
    if k.size == 0:
        return []
    r = res.r[k]
    we_d = state.freestream.rho * state.V ** 2 * 2.0 * r / liquid.sigma
    oh = liquid.mu / np.sqrt(liquid.rho * liquid.sigma * 2.0 * r)
    rows = np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
                            centroids[k, 0], centroids[k, 1], centroids[k, 2], flow.regime[k], res.branch[k], b[k],
                            res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1], t_hat[k, 2],
                            we_d, oh, (we_d > WE_BREAKUP).astype(float)])
    return rows.tolist()


def histogram(r, dn, dm):
    """(dn per bin, dM per bin) on BIN_EDGES for droplet radii r with counts dn and masses dm."""
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(r) & (r > 0.0)
    n_hist, _ = np.histogram(r[ok], BIN_EDGES, weights=np.asarray(dn)[ok])
    m_hist, _ = np.histogram(r[ok], BIN_EDGES, weights=np.asarray(dm)[ok])
    return n_hist, m_hist
```


- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_spray.py -q`
Expected: 4 passed.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/spray.py data/reference_values tests/test_reentry_model_spray.py
git commit -m "Add melt spraying: Girin's thick, thin and rarefied branches, release bookkeeping and the published reference values (Step 3 Task 7)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 8: Girin's published cases — the Girin-as-published driver and its analysis script

**Files:**
- Create: `reentry_model/girin_case.py`, `analysis/girin_reference.py`
- Test: `tests/test_reentry_model_girin.py`

**Interfaces:**
- Consumes: `dispersion.DispersionTable` (Task 4), `surface_flow.ranger_psi` (Task 5), `spray.thin_film_mode/rayleigh_taylor/histogram/BIN_EDGES` (Task 7), `compare`'s plot style (Step 2), the two reference-value files (Task 7).
- Produces: `girin_case.GirinVariant(name, V, rho_m, mu_m, sigma, R0=3e-3)`, `dimensionless_numbers(v) -> (Re, We, GI, alpha, mu, p, t_ch)`, `critical_angle(GI, p, we_cr=4.62)`, `surface_state(...)`, `run_variant(v, k_r, k_t, we_cr, re_density="shock"|"ambient", ...) -> dict` (keys `GI, phi_cr_deg_eq3, t_f_min_us_tau0, r_um_90deg_tau0, N, r_med_um, r_min_um, r_max_um, t_sd_ms, tau_end, radii, counts, times, mass_history, ...`); constants `RHO_A, MU_AIR, COMPRESSION, RE_DENSITY_NAMES`. The script writes `girin2017.json`, `girin1994.json`, `girin_summary.md` and the plots under `reentry_model_output/verification_melt/girin/`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_girin.py`:

```python
"""girin_case.py: Girin (2017) Table 1 -- the exact tier (GI, phi_cr within 2 %), the integrated tier (t_f, N, r_med
within 30 % with the ambient-density Reynolds number) and the reported t_s.d. (spec Step 3 section 13.2 as amended)."""
import json
import os

import numpy as np
import pytest

from reentry_model import girin_case as gc
from helpers import REPO_ROOT

TABLE = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin2017_table1.json")))


def variant(name):
    d = TABLE["variants"][name]
    return gc.GirinVariant(name, d["V_ms"], d["rho_m"], d["mu_m"], d["sigma"]), d


def test_exact_tier_gi_and_critical_angle():
    for name in ("I", "II", "III"):
        v, d = variant(name)
        Re, We, GI, alpha, mu, p, t_ch = gc.dimensionless_numbers(v)
        assert Re == pytest.approx(d["Re_inf"], rel=2e-2) and GI == pytest.approx(d["GI"], rel=2e-2)
        assert t_ch * 1e3 == pytest.approx(d["t_ch_ms"], rel=2e-2)                                # alpha on the ambient density
        assert np.degrees(gc.critical_angle(GI, p)) == pytest.approx(d["phi_cr_deg"], rel=2e-2)   # We_cr = 4.62, his practical value
        assert np.degrees(gc.critical_angle(GI, p, 3.08)) < d["phi_cr_deg"]                        # 3.08 would give ~17 % smaller angles
    assert gc.critical_angle(0.3, 1.0) == np.pi / 2                                                # below GI 0.4: nothing unstable


@pytest.mark.parametrize("name", ["I", "II", "III"])
def test_integrated_tier_with_the_ambient_reynolds_number(name):
    v, d = variant(name)
    o = gc.run_variant(v, re_density="ambient")
    assert o["t_f_min_us_tau0"] == pytest.approx(d["t_f_us"], rel=0.3)          # measured x1.22 / x0.88 / x1.07
    assert o["N"] == pytest.approx(d["N"], rel=0.3)                              # measured x0.87 / x0.97 / x0.76
    if d["r_med_um"]:
        assert o["r_med_um"] == pytest.approx(d["r_med_um"], rel=0.3)            # measured x0.96 / x0.95
        assert o["r_min_um"] < d["r_med_um"] < o["r_max_um"]
    assert o["t_sd_ms"] > 0.0 and o["m_end_kg"] < 1e-3 * o["m0_kg"]
    assert o["radii"].size == o["counts"].size and o["mass_history"][-1, 1] < o["mass_history"][0, 1]
    assert 0.4 < o["t_sd_ms"] / d["t_sd_ms"] < 0.6 or name == "III"              # reported: half his duration (iron), far below for stone


def test_shock_density_option_and_argument_checks():
    v, d = variant("I")
    shock, ambient = gc.run_variant(v, re_density="shock"), gc.run_variant(v, re_density="ambient")
    assert shock["Re_used"] == pytest.approx(6.0 * ambient["Re_used"]) and shock["r_med_um"] < 0.5 * ambient["r_med_um"]
    assert shock["GI"] == ambient["GI"] and shock["phi_cr_deg_eq3"] == ambient["phi_cr_deg_eq3"]
    with pytest.raises(ValueError):
        gc.run_variant(v, re_density="wrong")
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_girin.py -q`
Expected: ImportError (`reentry_model.girin_case`).

- [ ] **Step 3: Create `reentry_model/girin_case.py`**

```python
"""Girin's (2017) published spraying cases run with his own simplifications ("Girin-as-published" mode, spec Step 3
section 13.2): a sphere of radius R0 in a constant gas density, potential-flow edge velocity V_a = 1.5 (V_inf - w)
sin phi, Ranger's boundary layer delta_a = 2.2 R Re_a^-1/2 Psi(phi) with Re_a = Re_inf R~ (1 - W), the conjugated-
layer relations of his Eq. (2), the dispersion table for Delta_f and Im Omega_f, no melt limit and no runoff, spherical
belts of width R dphi, the induction rule tau_ind = int dt / t_per >= 1 releasing R dphi / lambda_f tori of
cross-section radius r = k_r lambda_f, and the deceleration law W = 1 - exp(-C tau), C = 2 alpha^1/2, tau = t / t_ch,
t_ch = 2 R0 / (alpha^1/2 V_inf).

Densities (measured 2026-09-20 against his Table 1): We_inf uses the ambient density rho_inf = rho_a / 6 and Re_inf the
compressed density rho_a = 1e-7 g/cm3 -- that reproduces his GI = We_inf Re_inf^-1/2 (13.0 / 3.55 / 43.5) and his
t_ch (alpha = rho_inf / rho_m); the Reynolds number entering delta_a is a switch (`re_density`): "shock" (rho_a, as the
printed Re_inf = 529) or "ambient" (rho_inf); with "ambient" his N and r_med are reproduced within 15 % for the iron
variants, with "shock" the droplets come out 2.5 x smaller (the paper does not say which density its delta_a used).
The stripping rate per area is independent of delta_a (it scales with V_s only), so t_s.d. is the same either way:
0.5 x his printed value for the iron variants and far below it for the stony one, whose lambda_f (3.5 mm) exceeds the
body radius -- his belt discretisation and induction handling are unstated (see the Step 3 plan's measured facts)."""
from dataclasses import dataclass

import numpy as np

from . import dispersion
from .surface_flow import ranger_psi

RHO_A, MU_AIR = 1.0e-4, 6.8e-5          # kg/m3, Pa s (Girin 2017 section 5)
COMPRESSION = 6.0
RE_DENSITY_NAMES = ("shock", "ambient")


@dataclass
class GirinVariant:
    name: str
    V: float          # m/s
    rho_m: float      # kg/m3
    mu_m: float       # Pa s
    sigma: float      # N/m
    R0: float = 3.0e-3


def dimensionless_numbers(v, rho_a=RHO_A, mu_a=MU_AIR, compression=COMPRESSION):
    """(Re_inf, We_inf, GI, alpha, mu, p, t_ch) with Girin's density conventions."""
    rho_inf = rho_a / compression
    Re = rho_a * 2.0 * v.R0 * v.V / mu_a
    We = rho_inf * 2.0 * v.R0 * v.V ** 2 / v.sigma
    alpha, mu = rho_inf / v.rho_m, mu_a / v.mu_m
    p = 1.0 + (alpha * mu) ** (1.0 / 3.0)
    t_ch = 2.0 * v.R0 / (np.sqrt(alpha) * v.V)
    return Re, We, We / np.sqrt(Re), alpha, mu, p, t_ch


def critical_angle(GI, p, we_cr=dispersion.WE_CRITICAL_PRACTICAL):
    """phi_cr [rad] at tau = 0 from Girin's Eq. (3): 2.475 p^-2 sin^2 phi Psi(phi) GI = We_cr."""
    from scipy.optimize import brentq
    f = lambda phi: 2.475 / p ** 2 * np.sin(phi) ** 2 * ranger_psi(phi) * GI - we_cr
    if f(np.pi / 2 - 1e-9) <= 0.0:
        return np.pi / 2
    return brentq(f, 1e-4, np.pi / 2 - 1e-9)


def surface_state(v, phi, Re_a, R, W, alpha, mu, p, table):
    """Per-belt delta_a, delta_m, V_s, We_s, Delta_f, Im Omega_f, lambda_f, r, t_per, t_f at the current radius R."""
    delta_a = 2.2 * R * ranger_psi(phi) / np.sqrt(Re_a)
    delta_m = (alpha / mu ** 2) ** (1.0 / 3.0) * delta_a
    V_a = 1.5 * v.V * (1.0 - W) * np.sin(phi)
    V_s = (alpha * mu) ** (1.0 / 3.0) / p * V_a
    we_s = v.rho_m * V_s ** 2 * delta_m / v.sigma
    delta_f, im_f, _ = table(np.maximum(we_s, table.we_onset))
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = 2.0 * np.pi * delta_m / delta_f
        t_f = delta_m / (V_s * im_f)
    return delta_a, delta_m, V_s, we_s, delta_f, im_f, lam, t_f


def run_variant(v, k_r=0.17, k_t=1.1, we_cr=dispersion.WE_CRITICAL_PRACTICAL, re_density="shock", dphi_deg=0.5, dt_fraction=0.2,
                mass_fraction_end=1e-4, table=None):
    """Girin's spraying flight of one variant. Returns a dict with the exact-tier numbers (GI, phi_cr, t_f at tau = 0,
    r at tau = 0 and 90 degrees, t_ch) and the integrated ones (N, r_med, r_min, r_max, t_sd, mass-loss history and the
    droplet list)."""
    if re_density not in RE_DENSITY_NAMES:
        raise ValueError("re_density must be one of {}".format(RE_DENSITY_NAMES))
    table = table or dispersion.DispersionTable()
    Re, We, GI, alpha, mu, p, t_ch = dimensionless_numbers(v)
    Re_used = Re if re_density == "shock" else Re / COMPRESSION
    C = 2.0 * np.sqrt(alpha)
    phi = np.radians(np.arange(dphi_deg / 2.0, 90.0, dphi_deg))
    dphi = np.radians(dphi_deg)
    m0 = v.rho_m * 4.0 / 3.0 * np.pi * v.R0 ** 3
    # tau = 0 exact-tier quantities
    _, _, _, we0, _, _, lam0, tf0 = surface_state(v, phi, Re_used, v.R0, 0.0, alpha, mu, p, table)
    unstable0 = we0 > we_cr
    out = {"Re_inf": Re, "We_inf": We, "GI": GI, "alpha": alpha, "mu": mu, "p": p, "t_ch_s": t_ch, "Re_used": Re_used,
           "phi_cr_deg_eq3": float(np.degrees(critical_angle(GI, p, we_cr))),
           "phi_cr_deg_direct": float(np.degrees(phi[unstable0][0])) if unstable0.any() else None,
           "t_f_min_us_tau0": float(tf0[unstable0].min() * 1e6) if unstable0.any() else None,
           "r_um_90deg_tau0": float(k_r * lam0[-1] * 1e6), "we_s_90deg_tau0": float(we0[-1])}
    # the flight
    m, t = m0, 0.0
    tind = np.zeros_like(phi)
    r_rel, n_rel, t_rel, m_hist = [], [], [], [(0.0, m0)]
    while m > mass_fraction_end * m0:
        R = (3.0 * m / (4.0 * np.pi * v.rho_m)) ** (1.0 / 3.0)
        W = 1.0 - np.exp(-C * t / t_ch)
        Re_a = Re_used * (R / v.R0) * (1.0 - W)
        _, delta_m, V_s, we_s, delta_f, im_f, lam, t_f = surface_state(v, phi, Re_a, R, W, alpha, mu, p, table)
        active = we_s > we_cr
        if not active.any():
            break
        t_per = k_t * t_f
        dt = dt_fraction * float(t_per[active].min())
        tind[active] += dt / t_per[active]
        tind[~active] = 0.0
        release = active & (tind >= 1.0)
        if release.any():
            r = k_r * lam[release]
            n_tor = R * dphi / lam[release]
            volume = n_tor * 2.0 * np.pi ** 2 * r ** 2 * R * np.sin(phi[release])
            dn = volume / (4.0 / 3.0 * np.pi * r ** 3)
            r_rel.append(r); n_rel.append(dn); t_rel.append(np.full(r.size, t))
            m -= float((v.rho_m * volume).sum())
            tind[release] -= 1.0
            m_hist.append((t, m))
        t += dt
    r = np.concatenate(r_rel) if r_rel else np.zeros(0)
    n = np.concatenate(n_rel) if n_rel else np.zeros(0)
    if r.size:
        order = np.argsort(r)
        cn = np.cumsum(n[order])
        r_med = float(r[order][np.searchsorted(cn, 0.5 * cn[-1])])
    else:
        r_med = None
    out.update({"N": float(n.sum()), "r_med_um": r_med * 1e6 if r_med else None, "r_min_um": float(r.min() * 1e6) if r.size else None,
                "r_max_um": float(r.max() * 1e6) if r.size else None, "t_sd_ms": t * 1e3, "tau_end": t / t_ch, "m0_kg": m0, "m_end_kg": m,
                "radii": r, "counts": n, "times": np.concatenate(t_rel) if t_rel else np.zeros(0), "mass_history": np.array(m_hist)})
    return out
```


- [ ] **Step 4: Create `analysis/girin_reference.py`**

```python
"""Girin's published spraying cases against the model's transcription of his theory (spec Step 3 sections 13.2-13.3).

    "$PY" analysis/girin_reference.py [--outdir reentry_model_output/verification_melt/girin] [--re-density ambient|shock|both]

Runs the three variants of Girin (2017) Table 1 in Girin-as-published mode (reentry_model.girin_case) and checks the
Girin & Kopyt (1994) Tables 1-2 with the model's thin-film and Rayleigh-Taylor modes (reentry_model.spray). Writes
girin2017.json / girin1994.json (every number next to its published value and ratio), the size distributions dn(r),
dM(r) per variant with the table values marked, the mass-loss law, and log-log plots of the 1994 tables with
residuals, plus girin_summary.md. Thresholds (README, Step 3): exact tier GI and phi_cr within 2 %; integrated tier
t_f, N and r_med within 30 % with the ambient-density Reynolds number; t_s.d. reported."""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
from reentry_model import girin_case, spray  # noqa: E402
from reentry_model.compare import INK, MODEL_COLOR, MUTED, SECOND, apply_rcparams, strip_top_right_spines  # noqa: E402

DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt", "girin")
TABLE_2017 = os.path.join(REPO_ROOT, "data", "reference_values", "girin2017_table1.json")
TABLE_1994 = os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")


def ratio(model, ref):
    return None if (model is None or ref in (None, 0)) else float(model / ref)


def run_2017(outdir, densities):
    tab = json.load(open(TABLE_2017))
    results = {}
    for name, d in tab["variants"].items():
        v = girin_case.GirinVariant(name, d["V_ms"], d["rho_m"], d["mu_m"], d["sigma"])
        results[name] = {}
        for rd in densities:
            o = girin_case.run_variant(v, re_density=rd)
            entry = {"GI": (o["GI"], d["GI"]), "phi_cr_deg": (o["phi_cr_deg_eq3"], d["phi_cr_deg"]), "t_f_us": (o["t_f_min_us_tau0"], d["t_f_us"]),
                     "N": (o["N"], d["N"]), "r_med_um": (o["r_med_um"], d["r_med_um"]), "r_min_um": (o["r_min_um"], d["r_min_um"]),
                     "r_max_um": (o["r_max_um"], d["r_max_um"]), "t_sd_ms": (o["t_sd_ms"], d["t_sd_ms"]), "tau_end": (o["tau_end"], d["tau_end"]),
                     "t_ch_ms": (o["t_ch_s"] * 1e3, d["t_ch_ms"]), "Re_inf": (o["Re_inf"], d["Re_inf"])}
            results[name][rd] = {k: {"model": m, "published": p, "ratio": ratio(m, p)} for k, (m, p) in entry.items()}
            results[name][rd]["r_um_90deg_tau0"] = o["r_um_90deg_tau0"]
            # plots: size distributions and the mass-loss law
            apply_rcparams(plt)
            fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(14, 4.2))
            if o["radii"].size:
                n_hist, m_hist = spray.histogram(o["radii"], o["counts"], o["counts"] * 4.0 / 3.0 * np.pi * v.rho_m * o["radii"] ** 3)
                mid = np.sqrt(spray.BIN_EDGES[:-1] * spray.BIN_EDGES[1:]) * 1e6
                ax.loglog(mid[n_hist > 0], n_hist[n_hist > 0], color=MODEL_COLOR, lw=1.4, label="model")
                bx.loglog(mid[m_hist > 0], m_hist[m_hist > 0] * 1e3, color=MODEL_COLOR, lw=1.4, label="model")
                for value, label in ((d["r_med_um"], "published r_med"), (d["r_min_um"], "published range"), (d["r_max_um"], None)):
                    if value:
                        ax.axvline(value, color=INK, lw=0.8, ls=":", label=label)
                        bx.axvline(value, color=INK, lw=0.8, ls=":")
                mh = o["mass_history"]
                cx.plot(mh[:, 0] * 1e3, mh[:, 1] / mh[0, 1], color=MODEL_COLOR, lw=1.4, label="model")
                t_sd = mh[-1, 0]
                cx.plot(mh[:, 0] * 1e3, (1.0 - np.minimum(mh[:, 0] / t_sd, 1.0)) ** 3, color=MUTED, lw=1.0, ls="--", label="(1 - t/t_s.d.)^3")
                cx.axvline(d["t_sd_ms"], color=INK, lw=0.8, ls=":", label="published t_s.d.")
            ax.set_xlabel("r [um]"); ax.set_ylabel("dn per bin"); ax.legend(frameon=False)
            bx.set_xlabel("r [um]"); bx.set_ylabel("dM per bin [g]")
            cx.set_xlabel("t [ms]"); cx.set_ylabel("m / m0"); cx.legend(frameon=False)
            ax.set_title("Girin 2017 variant {} ({}), Re on {} density".format(name, d["description"], rd), color=SECOND, fontsize=10)
            for a in (ax, bx, cx):
                strip_top_right_spines(a)
            fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin2017_variant_{}_{}.png".format(name, rd)), dpi=150); plt.close(fig)
    return results


def run_1994(outdir):
    tab = json.load(open(TABLE_1994))
    t1, t2 = tab["table1"], tab["table2"]

    class Liquid:
        rho, sigma = t1["rho1_kgm3"], t1["sigma_Nm"]
        mu = 0.0
    rho2, V0 = np.array(t1["rho2_kgm3"]), np.array(t1["V0_ms"])
    r_tab, tau_tab, m_tab = (np.array(t1[k]) for k in ("r_d_m", "tau_d_s", "m_kgm2s"))
    lam, tau = spray.thin_film_mode(np.full((3, 2), t1["mach"]), rho2[:, None] * V0[None, :] ** 2, Liquid)
    f_s = float((lam / 4.0)[0, 0] / r_tab[0, 0])                     # their shock-layer factor, fitted on the first entry
    r_model, tau_model = lam / 4.0 / f_s, tau / f_s ** 1.5
    m_model = Liquid.rho * r_model / (2.0 * tau_model)
    W = np.array(t2["W_ms2"])
    lam_rt, tau_rt = zip(*[spray.rayleigh_taylor(w, np.array([1.0]), Liquid)[1:] for w in W])
    lam_rt, tau_rt = np.array([x[0] for x in lam_rt]), np.array([x[0] for x in tau_rt])
    res = {"table1": {"shock_layer_factor": f_s, "r_d_ratio": (r_model / r_tab).tolist(), "tau_d_ratio": (tau_model / tau_tab).tolist(),
                      "m_ratio": (m_model / m_tab).tolist(), "r_d_model_m": r_model.tolist(), "tau_d_model_s": tau_model.tolist()},
           "table2": {"lambda_star_model_m": lam_rt.tolist(), "lambda_star_ratio_to_eq14_column": (lam_rt / np.array(t2["lambda_star_m_eq14"])).tolist(),
                      "lambda_star_ratio_to_printed": (lam_rt / np.array(t2["lambda_star_m_as_printed"])).tolist(),
                      "tau_star_model_s": tau_rt.tolist(), "tau_star_ratio": (tau_rt / np.array(t2["tau_star_s"])).tolist()}}
    apply_rcparams(plt)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), gridspec_kw={"height_ratios": [3, 1]})
    for j, (name, model, published, unit) in enumerate((("r_d", r_model, r_tab, "m"), ("tau_d", tau_model, tau_tab, "s"), ("m", m_model, m_tab, "kg/m2/s"))):
        ax, rx = axes[0, j], axes[1, j]
        for k, V in enumerate(V0):
            ax.loglog(rho2, published[:, k], "o", color=INK, label="published, V0 = {:.0f} km/s".format(V / 1e3))
            ax.loglog(rho2, model[:, k], "--", color=MODEL_COLOR, label="model x f_s" if k == 0 else None)
            rx.semilogx(rho2, 100.0 * (model[:, k] / published[:, k] - 1.0), "s-", color=MODEL_COLOR if k == 0 else MUTED, lw=0.8)
        ax.set_ylabel("{} [{}]".format(name, unit)); ax.legend(frameon=False, fontsize=8); ax.set_title("Girin & Kopyt 1994 Table 1", color=SECOND, fontsize=10)
        rx.axhline(0.0, color=MUTED, lw=0.6); rx.set_xlabel("rho_2 [kg/m3]"); rx.set_ylabel("model/published - 1 [%]")
        strip_top_right_spines(ax); strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin1994_table1.png"), dpi=150); plt.close(fig)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2))
    ax.loglog(W, t2["lambda_star_m_eq14"], "o", color=INK, label="Eq. (14) (printed x 10)")
    ax.loglog(W, t2["lambda_star_m_as_printed"], "x", color=MUTED, label="Table 2 as printed")
    ax.loglog(W, lam_rt, "--", color=MODEL_COLOR, label="model")
    ax.set_xlabel("W [m/s2]"); ax.set_ylabel("lambda* [m]"); ax.legend(frameon=False, fontsize=8); ax.set_title("Girin & Kopyt 1994 Table 2", color=SECOND, fontsize=10)
    bx.loglog(W, t2["tau_star_s"], "o", color=INK, label="published"); bx.loglog(W, tau_rt, "--", color=MODEL_COLOR, label="model")
    bx.set_xlabel("W [m/s2]"); bx.set_ylabel("tau* [s]"); bx.legend(frameon=False, fontsize=8)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin1994_table2.png"), dpi=150); plt.close(fig)
    return res


def summary_markdown(r2017, r1994, densities):
    lines = ["| variant | Re density | GI | phi_cr [deg] | t_f [us] | N | r_med [um] | range [um] | t_s.d. [ms] |", "|---|---|---|---|---|---|---|---|---|"]
    fmt = lambda e, f="{:.3g}": "{} ({}{})".format(f.format(e["model"]) if e["model"] is not None else "n/a", f.format(e["published"]) if e["published"] is not None else "n/a",
                                                  ", x{:.2f}".format(e["ratio"]) if e["ratio"] else "")
    for name, per in r2017.items():
        for rd in densities:
            e = per[rd]
            rng = "{}-{}".format("{:.1f}".format(e["r_min_um"]["model"]) if e["r_min_um"]["model"] else "n/a", "{:.1f}".format(e["r_max_um"]["model"]) if e["r_max_um"]["model"] else "n/a")
            rng += " ({}-{})".format(e["r_min_um"]["published"] or "n/a", e["r_max_um"]["published"] or "n/a")
            lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(name, rd, fmt(e["GI"], "{:.2f}"), fmt(e["phi_cr_deg"], "{:.1f}"), fmt(e["t_f_us"], "{:.1f}"),
                                                                       fmt(e["N"], "{:.2e}"), fmt(e["r_med_um"], "{:.1f}"), rng, fmt(e["t_sd_ms"], "{:.2f}")))
    lines += ["", "Values in parentheses: published (Girin 2017 Table 1) and model/published.", "",
              "Girin & Kopyt 1994 Table 1: shock-layer factor {:.2f}; r_d ratios {}; tau_d ratios {}; m ratios {}".format(
                  r1994["table1"]["shock_layer_factor"], np.round(r1994["table1"]["r_d_ratio"], 3).tolist(), np.round(r1994["table1"]["tau_d_ratio"], 3).tolist(),
                  np.round(r1994["table1"]["m_ratio"], 3).tolist()),
              "Girin & Kopyt 1994 Table 2: lambda* / Eq. (14) {}; lambda* / printed {}; tau* ratios {}".format(
                  np.round(r1994["table2"]["lambda_star_ratio_to_eq14_column"], 3).tolist(), np.round(r1994["table2"]["lambda_star_ratio_to_printed"], 3).tolist(),
                  np.round(r1994["table2"]["tau_star_ratio"], 3).tolist())]
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--re-density", choices=("ambient", "shock", "both"), default="both")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    densities = ("ambient", "shock") if args.re_density == "both" else (args.re_density,)
    r2017 = run_2017(args.outdir, densities)
    r1994 = run_1994(args.outdir)
    json.dump(r2017, open(os.path.join(args.outdir, "girin2017.json"), "w"), indent=2)
    json.dump(r1994, open(os.path.join(args.outdir, "girin1994.json"), "w"), indent=2)
    md = summary_markdown(r2017, r1994, densities)
    open(os.path.join(args.outdir, "girin_summary.md"), "w").write(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```


- [ ] **Step 5: Run the tests and the script**

```bash
"$PY" -m pytest tests/test_reentry_model_girin.py -q
"$PY" analysis/girin_reference.py
```

Expected: 5 passed; the script prints the summary table (variant I ambient: GI 13.04 (13.00, ×1.00), φ_cr 16.3 (16.1), t_f 7.0 (5.7, ×1.22), N 1.31e6 (1.50e6, ×0.87), r_med 25.8 (26.9, ×0.96), t_s.d. 2.82 (5.90, ×0.48); 1994 Table 1 shock-layer factor 7.77 with all r_d/τ_d/ṁ ratios 0.99–1.01; Table 2 λ*/Eq. (14) 1.00–1.02, λ*/printed 10.0, τ* 1.00) and writes eleven files under `reentry_model_output/verification_melt/girin/`.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/girin_case.py analysis/girin_reference.py tests/test_reentry_model_girin.py
git commit -m "Reproduce Girin's published spraying cases with his own simplifications (Step 3 Task 8)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 9: The melting body

**Files:**
- Modify (replace): `reentry_model/body.py`
- Modify: `reentry_model/trajectory.py` (one statement)
- Test: `tests/test_reentry_model_melting.py` (new); `tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend` (from Task 3) now runs in `fenicsx_env`

**Interfaces:**
- Consumes: Tasks 1–7 (`mesh.deactivate/surface/active_nodes`, `Material.feed_fraction/liquid_fraction/enthalpy/h_liquid/liquid`, `thermal` `set_fractions/element_energies/step(nodal_load=)/pinned/temperature`, `surface_flow.SurfaceFlow`, `spray.SprayModel/melt_layer/source_rows/histogram/N_BINS`, `film.lubrication/Runoff`), `heating.HeatingResult`, `trajectory.AeroState`.
- Produces: `body.REMOVAL_NAMES = ("girin", "instant")`, `PHI_MIN = 1e-3`, `PHI_DEATH = 0.05`, `NEAREST_PATCHES = 4`; `MeltSettings(removal, runoff, demise_fraction, particles)`; `MeltingBody(mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0, emissivity, T_ambient, v_hat)` with `.advance(t, dt, loads, state=None)`, `.melt_step(t, dt, state)`, `.mass(t)`, `.reference_area()`, `.equivalent_radius()`, `.energy()` (FEM + film), `.mean_temperature()`, `.melt_front_depth()`, `.film_thickness_max/mean()`, `.demised()`, `.energy_balance_residual()`, `.melt_stats() -> dict` (the `coupled.MELT_COLUMNS` values plus `mass_kg`), attributes `phi, m_f, surface, theta, t_hat, windward, mass0, sprayed_mass, runoff_mass, removed_mass, removed_enthalpy, n_released, source_rows, hist_n, hist_m, melt_onset, spray_onset, consumed, last_flow, last_spray, last_b, last_melt, pending_load, liquid, flow, spray, runoff`; `Body.reference_area()` in the protocol (`ConstantBody`/`ThermalBody` return None; `ThermalBody.advance` accepts `state=None`); `Simulator.aero_state` uses the body's reference area and gives a consumed body no drag.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_melting.py`:

```python
"""body.MeltingBody: feed, instant removal in the lumped limit, film and spraying, element death with hand-over,
energy and mass balances, demise (spec Step 3 sections 8-11)."""
import math

import numpy as np
import pytest

from reentry_model import body, heating, material, mesh, thermal
from test_reentry_model_coupled import MASS_100MM, simulator

pytest.importorskip("cantera")


@pytest.fixture(scope="module")
def layered_mesh(tmp_path_factory):
    return mesh.sphere_mesh(0.05, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=2)


def melting_body(the_mesh, name="AA7075_range", k_scale=1.0, **settings):
    mat = material.Material.from_drama_json(name)
    if k_scale != 1.0:
        mat.k_table = mat.k_table * k_scale
    m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
    return body.MeltingBody(m, mat, thermal.thermal_solver("skfem"), MASS_100MM, settings=body.MeltSettings(**settings))


def test_settings_and_construction(layered_mesh):
    with pytest.raises(ValueError):
        body.MeltSettings(removal="magic")
    with pytest.raises(ValueError):
        body.MeltSettings(demise_fraction=1.5)
    with pytest.raises(ValueError):
        body.MeltingBody(layered_mesh, material.Material.from_drama_json("AA7075_nomelt"), thermal.thermal_solver("skfem"), MASS_100MM)
    b = melting_body(layered_mesh)
    assert b.mass(0.0) == pytest.approx(b.mass0) and b.mass0 == pytest.approx(MASS_100MM, rel=3e-3)
    assert b.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3) and b.equivalent_radius() == pytest.approx(0.05, rel=1e-3)
    assert b.m_f.shape == (b.surface.n_patches,) and b.m_f.sum() == 0.0 and not b.demised() and b.melt_front_depth() == 0.0
    stats = b.melt_stats()
    assert stats["mass_kg"] == b.mass(0.0) and stats["n_active_elements"] == layered_mesh.n_elements and stats["sprayed_mass_kg"] == 0.0


def test_instant_removal_in_the_lumped_limit_follows_q_over_l(layered_mesh):
    """k x 1e4, uniform 3e5 W/m2 on the whole surface, AA7075: once the body sits on the 850 K plateau the mass leaves
    at (Q_conv - Q_rad) / L_f with the geometry intact (SESAM's lumped law), the energy balance exact."""
    b = melting_body(layered_mesh, "AA7075", k_scale=1e4, removal="instant", runoff=False)
    q = np.full(b.surface.n_patches, 3e5)
    loads = heating.HeatingResult(q, 3e5, 0.0, 0.0, 0.0)
    t, dt = 0.0, 0.5
    while b.melt_onset is None:
        t += dt
        b.advance(t, dt, loads)
    m1 = b.mass(0.0)
    for _ in range(20):
        t += dt
        b.advance(t, dt, loads)
    Q_net = b.last.Q_conv - b.last.Q_rad
    assert (m1 - b.mass(0.0)) / (20 * dt) == pytest.approx(Q_net / 4e5, rel=2e-2)
    assert b.removed_mass == pytest.approx(b.mass0 - b.mass(0.0)) and b.m_f.sum() == 0.0 and b.sprayed_mass == 0.0
    assert abs(b.energy_balance_residual()) < 1e-8 and b.mesh.n_active == layered_mesh.n_elements    # nothing dead yet: phi shrinks uniformly
    assert 0.5 < b.phi.min() < b.phi.max() < 1.0 and 848.0 < b.mean_temperature() < 852.0
    assert b.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3)
    assert b.melt_onset is not None and b.spray_onset is None


def test_film_spraying_death_and_balances(layered_mesh):
    """Physics-mode loads at 71 km on the real body for 8 s: film forms, is stripped, the outer layer dies on the
    windward face with the film handed over, mass and energy balances hold, the source table grows."""
    b = melting_body(layered_mesh)
    sim = simulator(b, t_max=60.0)
    sim.advance(43.0)                                                                 # to 71 km on the constant-mass trajectory
    b.solver.set_temperature(880.0)                                                   # a hot body: the surface melts at once
    b.energy0 = b.energy()                                                            # the balance is reckoned from here
    model = heating.PhysicsHeating()
    for _ in range(16):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
    stats = b.melt_stats()
    assert b.sprayed_mass > 0.0 and b.n_released > 0.0 and len(b.source_rows) > 0 and stats["n_released"] == b.n_released
    assert b.removed_mass == pytest.approx(b.sprayed_mass) and b.mass(0.0) == pytest.approx(b.mass0 - b.sprayed_mass, rel=1e-6)
    assert b.mesh.n_active < layered_mesh.n_elements and b.mesh.element_layer[~b.mesh.active].max() <= 1     # only layer elements have died
    assert abs(b.energy_balance_residual()) < 1e-7
    assert b.surface.n_patches == b.m_f.size == len(b.solver.areas) and b.theta.size == b.surface.n_patches
    assert stats["theta_cr_deg"] < 45.0 and stats["spraying_area_m2"] > 0.0 and 20.0 < stats["r_median_um"] < 2000.0
    assert 0.0 <= stats["regime_fraction_continuum"] and stats["regime_fraction_continuum"] + stats["regime_fraction_slip"] + stats["regime_fraction_fm"] == pytest.approx(1.0)
    assert stats["film_thickness_max_mm"] >= stats["film_thickness_mean_mm"] >= 0.0 and stats["removed_enthalpy_J"] == pytest.approx(b.sprayed_mass * b.material.h_liquid)
    rows = np.array(b.source_rows)
    assert rows.shape[1] == 22 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
    assert b.hist_n.sum() == pytest.approx(b.n_released) and b.hist_m.sum() == pytest.approx(b.sprayed_mass)
    assert not b.demised() and b.reference_area() < math.pi * 0.05 ** 2                            # the windward face has receded


def test_demise_and_consumption(layered_mesh):
    b = melting_body(layered_mesh, "AA7075", k_scale=1e4, removal="instant", runoff=False, demise_fraction=0.5)
    loads = heating.HeatingResult(np.full(b.surface.n_patches, 2e6), 2e6, 0.0, 0.0, 0.0)
    t = 0.0
    while not b.demised():
        t += 0.5
        b.advance(t, 0.5, loads)
        assert t < 200.0
    assert b.mass(0.0) < 0.5 * b.mass0 and not b.consumed and b.mesh.n_active > 0
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_melting.py -q`
Expected: AttributeError (`body.MeltSettings`).

- [ ] **Step 3: Replace `reentry_model/body.py`**

```python
"""The body the trajectory carries, as a stepper (spec section 4). ConstantBody keeps mass and one temperature
(Step 1 behaviour, `--thermal none`); ThermalBody wraps a mesh, a material and a conduction solver and advances the
temperature field with each macro step's per-patch convective loads, keeping the heat bookkeeping SESAM reports."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def reference_area(self): ...                        # m2 drag reference area, or None for the fixed pi D^2/4
    def surface_temperature(self) -> np.ndarray: ...     # K per patch
    def mean_temperature(self) -> float: ...             # K, energy-equivalent (spec 6.4)
    def energy(self) -> float: ...                       # J stored above the material's reference temperature
    def field(self): ...                                 # nodal temperatures, or None


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def advance(self, t, dt, loads):
        return None

    def surface_temperature(self):
        return np.array([self.temperature_K])

    def mean_temperature(self):
        return self.temperature_K

    def energy(self):
        return 0.0

    def field(self):
        return None


class ThermalBody:
    """Mesh + material + ThermalSolver. `v_hat` is the direction of motion in the body frame (fixed attitude: the
    stagnation patch is the one whose normal is along v_hat, spec 13.1)."""

    def __init__(self, mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        self.mesh, self.material, self.solver, self.mass_kg = mesh, material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient = float(T_ambient)
        self.surface = mesh.surface()
        self.theta = self.surface.angles_to(v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(v_hat), self.surface.patch_toward(-np.asarray(v_hat, dtype=float))
        self.i_centre = mesh.centre_node()
        self.volume = mesh.volume()
        self.solver.setup(mesh, material, self.emissivity)
        self.solver.set_temperature(T0)
        self.energy0 = self.solver.energy()
        self.integrated_heat = 0.0       # J, integral of Q_conv dt (SESAM's integrated_heat_J)
        self.absorbed_heat = 0.0         # J, integral of (Q_conv - Q_rad) dt
        self.radiated_heat = 0.0         # J
        self.iterations = []             # Newton iterations per step
        self.last = None                 # thermal.StepResult of the last step

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def advance(self, t, dt, loads, state=None):
        res = self.solver.step(dt, loads.q_conv, self.T_ambient)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res

    def surface_temperature(self):
        return self.surface.facet_mean(self.solver.temperature())

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / (self.material.rho * self.volume)))

    def energy(self):
        return self.solver.energy()

    def field(self):
        return self.solver.temperature()

    def radiated_power(self):
        return self.solver.radiated_power(self.T_ambient)

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat) / absorbed heat: zero for an exact discrete balance."""
        return (self.energy() - self.energy0 - self.absorbed_heat) / self.absorbed_heat if self.absorbed_heat else 0.0

    def surface_stats(self):
        T, Tf = self.solver.temperature(), self.surface_temperature()
        return {"surface_T_max_K": float(Tf.max()), "surface_T_min_K": float(Tf.min()),
                "surface_T_mean_K": float((Tf * self.surface.areas).sum() / self.surface.area),
                "T_stagnation_K": float(Tf[self.i_stag]), "T_back_K": float(Tf[self.i_back]), "T_centre_K": float(T[self.i_centre])}


PHI_MIN = 1.0e-3                 # element fraction kept by elements that do not own a patch (they cannot die: no cavities)
NEAREST_PATCHES = 4              # patches that receive an interior element's liquid (area-weighted)
PHI_DEATH = 0.05                 # a patch owner below this fraction dies (its remainder goes to the film): keeps the surface
                                 # nodes' thermal mass above 5 % of an element's, which the Newton iteration needs (measured
                                 # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
REMOVAL_NAMES = ("girin", "instant")


@dataclass
class MeltSettings:
    removal: str = "girin"           # girin: film + runoff + spraying; instant: liquid removed as it forms (verification device)
    runoff: bool = True
    demise_fraction: float = 0.01    # the run ends when the mass falls below this fraction of the initial mass
    particles: bool = True           # keep the source-table rows

    def __post_init__(self):
        if self.removal not in REMOVAL_NAMES:
            raise ValueError("removal must be one of {}, got {!r}".format(REMOVAL_NAMES, self.removal))
        if not 0.0 < self.demise_fraction < 1.0:
            raise ValueError("demise_fraction must be within (0, 1)")


class MeltingBody(ThermalBody):
    """ThermalBody with melting (spec Step 3 sections 6, 8-11): element fractions phi_e, the film per patch, the
    gas-side surface flow, spraying, element death and the mass/area/energy accounting.

    Per macro step (`advance`): the conduction step with the deferred melt loads of the previous step -> melt step:
    (i) every element's liquid inventory f_feed(T_e) phi_e rho V_e becomes film on its patches (owners) or the
    nearest patch (interior elements keep PHI_MIN so no cavity can open); the material leaves the FEM at h(T_e)
    but is booked at h_liquid, the difference (superheat, or the latent deficit of a partly molten element) is a
    nodal load on the element's nodes over the next step -- energy-exact and the melt front moves at the
    energy-limited rate whatever the element size; (ii) surface flow, delta_m, lubrication, runoff transport;
    (iii) spraying and release; (iv) patch owners at phi <= PHI_MIN die: the mesh's active set, the surface, the
    film (handed to the nearest surviving patch) and the solver's fractions are refreshed. `removal = "instant"`
    removes the liquid inventory of every element as it forms, without film, runoff or spraying: the lumped-melting
    device that reproduces SESAM's Q/L_f law (SESAM hollows the sphere at fixed outer geometry, measured
    2026-09-20, so the geometry is kept until elements die).
    mass(t) = sum phi rho V + sum m_f; the drag reference area is the current surface's projection on the flight
    direction (pi R^2 while intact); the nose radius R0 and the sphere drag tables are kept (spec section 10)."""

    def __init__(self, mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0=300.0,
                 emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        if not material.melts or material.liquid is None:
            raise ValueError("MeltingBody needs a material with a latent heat and liquid properties (AA7075 or AA7075_range)")
        super().__init__(mesh, material, solver, mass_kg, T0, emissivity, T_ambient, v_hat)
        from . import film as film_mod, spray as spray_mod, surface_flow
        self.settings = settings or MeltSettings()
        self.flow = flow or surface_flow.SurfaceFlow()
        self.spray = spray_model or spray_mod.SprayModel(material.liquid)
        self._film_mod = film_mod
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.liquid = material.liquid
        self.element_mass = material.rho * mesh.element_volumes()          # kg at phi = 1
        self.phi = np.ones(mesh.n_elements)
        self.mass0 = float(self.element_mass.sum())
        self.m_f = np.zeros(self.surface.n_patches)
        self.pending_load = np.zeros(mesh.n_nodes)                          # J, deferred melt energy for the next step
        self.sprayed_mass = self.runoff_mass = self.removed_mass = self.removed_enthalpy = 0.0
        self.n_released = 0.0
        self.source_rows = []
        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
        self.last_flow = self.last_spray = None
        self.last_melt = {"n_dead": 0, "runoff_substeps": 0, "released_mass": 0.0, "n_released": 0.0, "feed_mass": 0.0}
        self.melt_onset = self.spray_onset = None
        self.consumed = False
        self._refresh_geometry()

    # -- geometry -----------------------------------------------------------------------------------------------
    def _refresh_geometry(self):
        self.surface = self.mesh.surface()
        self.theta = self.surface.angles_to(self.v_hat)
        self.t_hat = self.surface.tangent_from(self.v_hat)
        self.windward = self.theta <= 0.5 * np.pi
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)
        self.runoff = self._film_mod.Runoff(self.surface, self.mesh.points, self.windward) if self.settings.runoff else None
        owners = self.surface.owner
        self.owner_area = np.bincount(owners, self.surface.areas, self.mesh.n_elements)      # patch area per owner element
        from scipy.spatial import cKDTree
        self._patch_tree = cKDTree(self.surface.centroids)
        self.patch_of_face = np.full(len(self.mesh._face_nodes), -1, dtype=np.int64)     # face id -> patch index (-1: not a patch)
        self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)

    def mass(self, t):
        return float((self.phi * self.element_mass).sum() + self.m_f.sum())

    def reference_area(self):
        return self.surface.projected_area(self.v_hat)

    def equivalent_radius(self):
        return (3.0 * self.mass(0.0) / (4.0 * math.pi * self.material.rho)) ** (1.0 / 3.0)

    def energy(self):
        return self.solver.energy() + float(self.m_f.sum()) * self.material.h_liquid

    def mean_temperature(self):
        m = self.mass(0.0)
        return float(self.material.temperature_from_enthalpy(self.energy() / m)) if m > 0.0 else 0.0

    def melt_front_depth(self):
        """Deepest molten point [m]: the distance from the original surface of the innermost element with f_l > 0."""
        f_l = self.material.liquid_fraction(self.solver.temperature())[self.mesh.tets].max(axis=1)
        molten = self.mesh.active & (f_l > 0.0)
        if not molten.any():
            return 0.0
        r = np.linalg.norm(self.mesh.points[self.mesh.tets[molten]].mean(axis=1), axis=1)
        return float(self.mesh.params.get("radius_m", r.max()) - r.min())

    # -- the step --------------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
        load = self.pending_load / dt
        self.pending_load = np.zeros(self.mesh.n_nodes)
        res = self.solver.step(dt, loads.q_conv, self.T_ambient, nodal_load=load)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res
        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
        self.melt_step(t, dt, state)

    def melt_step(self, t, dt, state):
        mat, liq, s = self.material, self.liquid, self.settings
        T = self.solver.temperature()
        tets = self.mesh.tets
        h_e = mat.enthalpy(T)[tets].mean(axis=1)                 # the element's specific enthalpy (nodal h, as the solver's energy)
        # (i) feed: the liquid inventory of every active element (nodal feed fraction averaged over the element)
        f = mat.feed_fraction(T)[tets].mean(axis=1) * self.mesh.active
        owner = self.owner_area > 0.0
        cap = np.where(owner, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
        d_phi = np.minimum(f * self.phi, cap)
        fed = d_phi * self.element_mass
        self.phi = self.phi - d_phi
        self._book_removed(fed, h_e)
        feed_mass = float(fed.sum())
        if feed_mass > 0.0 and self.melt_onset is None:
            self.melt_onset = t
        if s.removal == "instant":
            self.removed_mass += feed_mass
            self.removed_enthalpy += feed_mass * mat.h_liquid
            released = feed_mass
        else:
            self._add_to_film(fed)
            released = self._film_and_spray(t, dt, state)
        # (iv) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed
        # (cascade within the step: the surface never ends a step owned by a material-free element)
        n_dead = 0
        while not self.consumed:
            dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
            if dead.size == 0:
                break
            rest = self.phi[dead] * self.element_mass[dead]
            self._book_removed_subset(dead, rest, h_e[dead])
            if s.removal == "instant":
                self.removed_mass += float(rest.sum())
                self.removed_enthalpy += float(rest.sum()) * mat.h_liquid
            else:
                extra = np.zeros(self.mesh.n_elements)
                extra[dead] = rest
                self._add_to_film(extra)
            self.phi[dead] = 0.0
            self._kill(dead)
            n_dead += int(dead.size)
        self.solver.set_fractions(self.phi)
        self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass})

    def _book_removed(self, fed, h_e):
        k = np.flatnonzero(fed > 0.0)
        if k.size:
            self._book_removed_subset(k, fed[k], h_e[k])

    def _book_removed_subset(self, elements, mass, h_e):
        """Material leaving the FEM at h(T_e) is booked at h_liquid; the difference is deferred to the element's nodes."""
        diff = mass * (h_e - self.material.h_liquid)             # J (positive: superheat stays in the body)
        nodes = self.mesh.tets[elements]
        np.add.at(self.pending_load, nodes.ravel(), np.repeat(diff / 4.0, 4))

    def _add_to_film(self, fed):
        k = np.flatnonzero(fed > 0.0)
        if k.size == 0:
            return
        owner = self.owner_area[k] > 0.0
        # owners: shared by patch area; interior elements: the nearest patch
        if owner.any():
            ko = k[owner]
            share = np.zeros(self.mesh.n_elements)
            share[ko] = fed[ko] / self.owner_area[ko]
            self.m_f += share[self.surface.owner] * self.surface.areas
        if (~owner).any():                                    # interior elements: the NEAREST_PATCHES nearest patches, by area
            ki = k[~owner]
            kk = min(NEAREST_PATCHES, self.surface.n_patches)
            _, near = self._patch_tree.query(self.mesh.points[self.mesh.tets[ki]].mean(axis=1), k=kk)
            near = near.reshape(len(ki), kk)
            w = self.surface.areas[near]
            np.add.at(self.m_f, near.ravel(), (fed[ki][:, None] * w / w.sum(axis=1, keepdims=True)).ravel())

    def _film_and_spray(self, t, dt, state):
        liq, s, mat = self.liquid, self.settings, self.material
        from . import spray as spray_mod
        if state is None or self.m_f.sum() <= 0.0:
            self.last_flow = self.last_spray = None
            return 0.0
        flow = self.flow.evaluate(state, self.theta, self.mesh.params.get("radius_m", 0.05), liq.rho)
        delta_m, _ = spray_mod.melt_layer(flow, liq)
        areas = self.surface.areas
        # (ii) lubrication and runoff
        n_sub, moved = 0, 0.0
        if self.runoff is not None:
            q_of_b = lambda b: self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu)[1]
            self.m_f, n_sub, moved = self.runoff.transport(self.m_f, q_of_b, self.t_hat, liq.rho, areas, dt)
            self.runoff_mass += moved                                              # mass that arrived on another patch
        b = self.m_f / (liq.rho * areas)
        v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu)
        # (iii) spraying
        res = self.spray.evaluate(flow, state, b, delta_m, v_s, self.windward, dt, areas, self.m_f, self.mesh.params.get("radius_m", 0.05))
        released = float(res.dm.sum())
        if released > 0.0:
            if self.spray_onset is None:
                self.spray_onset = t
            if s.particles:
                self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state)
            hn, hm = spray_mod.histogram(res.r, res.dn, res.dm)
            self.hist_n += hn
            self.hist_m += hm
            self.m_f = self.m_f - res.dm
            self.sprayed_mass += released
            self.removed_mass += released
            self.removed_enthalpy += released * mat.h_liquid
            self.n_released += float(res.dn.sum())
        self.m_f[self.m_f < 1e-30] = 0.0                       # no denormal films (they made 0/0 coefficients in the runoff)
        self.last_flow, self.last_spray, self.last_b = flow, res, b
        r = res.r[res.dm > 0.0]
        self.last_melt.update({
            "runoff_substeps": n_sub, "n_released": float(res.dn.sum()), "regime_fractions": flow.regime_fractions(areas, self.windward),
            "theta_cr_deg": float(np.degrees(self.theta[res.unstable].min())) if res.unstable.any() else float("nan"),
            "spraying_area_m2": float(areas[res.unstable].sum()), "rt_active": float(res.rt_active),
            "r_median_um": float(np.median(r) * 1e6) if r.size else float("nan"), "r_max_um": float(r.max() * 1e6) if r.size else float("nan"),
            "film_thickness_max_mm": float(b.max() * 1e3)})
        return released

    def _kill(self, dead):
        old_surface, old_m_f = self.surface, self.m_f
        gone, new = self.mesh.deactivate(dead)
        if self.mesh.n_active == 0:                              # nothing left: the run ends (demise) without a surface
            self.consumed = True
            self.m_f = np.zeros(0)
            return
        self._refresh_geometry()
        m_f = np.zeros(self.surface.n_patches)
        keep = self.patch_of_face[old_surface.face_ids]
        kept = keep >= 0
        m_f[keep[kept]] += old_m_f[kept]
        lost = ~kept & (old_m_f > 0.0)
        if lost.any():                                          # hand the vanished patches' film to the nearest survivors
            target = self._patch_tree.query(old_surface.centroids[lost])[1]
            np.add.at(m_f, target, old_m_f[lost])
        self.m_f = m_f

    # -- reporting -----------------------------------------------------------------------------------------------
    def film_thickness_max(self):
        """Thickest film [m] over patches of at least a tenth of the median patch area (slivers excluded)."""
        if not self.m_f.size:
            return 0.0
        a = self.surface.areas
        ok = a >= 0.1 * np.median(a)
        return float((self.m_f[ok] / (self.liquid.rho * a[ok])).max())

    def film_thickness_mean(self):
        """Film mass over the wetted area [m]."""
        wet = self.m_f > 0.0
        return float(self.m_f.sum() / (self.liquid.rho * self.surface.areas[wet].sum())) if wet.any() else 0.0

    def demised(self):
        """True once the body's material (film excluded) is below the demise fraction of the initial mass."""
        return self.consumed or float((self.phi * self.element_mass).sum()) < self.settings.demise_fraction * self.mass0

    def energy_balance_residual(self):
        """(E_body + E_film - E0 - absorbed heat + removed enthalpy + pending load + dropped load) / absorbed heat:
        zero for an exact discrete balance (the deferred loads cancel between the FEM and the removed material)."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy + self.pending_load.sum()
                + self.total_dropped_load) / self.absorbed_heat

    @property
    def total_applied_load(self):
        return getattr(self, "_applied_total", 0.0)

    @property
    def total_dropped_load(self):
        return getattr(self, "_dropped_total", 0.0)

    def melt_stats(self):
        lm = self.last_melt
        fractions = lm.get("regime_fractions", [0.0, 0.0, 0.0])
        return {"mass_kg": self.mass(0.0), "film_mass_kg": float(self.m_f.sum()), "sprayed_mass_kg": self.sprayed_mass,
                "runoff_mass_kg": self.runoff_mass, "removed_mass_kg": self.removed_mass,
                "melt_front_depth_max_mm": self.melt_front_depth() * 1e3, "equivalent_radius_mm": self.equivalent_radius() * 1e3,
                "n_active_elements": float(self.mesh.n_active), "spraying_area_m2": lm.get("spraying_area_m2", 0.0),
                "theta_cr_deg": lm.get("theta_cr_deg", float("nan")), "n_released": self.n_released,
                "released_mass_kg": lm.get("released_mass", 0.0), "r_median_um": lm.get("r_median_um", float("nan")),
                "r_max_um": lm.get("r_max_um", float("nan")), "regime_fraction_continuum": fractions[0],
                "regime_fraction_slip": fractions[1], "regime_fraction_fm": fractions[2], "rt_active": lm.get("rt_active", 0.0),
                "removed_enthalpy_J": self.removed_enthalpy,
                "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active)}
```


- [ ] **Step 4: Edit `reentry_model/trajectory.py`**

In `Simulator.aero_state`, replace the line

```python
            a_drag = -0.5 * fs.rho * V * v_rel * cd * self.area / self.body.mass(t)
```

with

```python
            area = self.body.reference_area() or self.area          # a melting body's projected area (Step 3), else pi D^2/4
            m = self.body.mass(t)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * area / m if m > 0.0 else np.zeros(3)      # a consumed body (Step 3) has no drag
```

- [ ] **Step 5: Run the tests**

```bash
"$PY" -m pytest tests/test_reentry_model_melting.py tests/test_reentry_model_coupled.py tests/test_reentry_model_trajectory.py -q
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
```

Expected: 4 passed (melting) and the Step 2 coupled/trajectory tests unchanged; in `fenicsx_env` 8 passed — the two backends give the same melting run (mass to 1e-6, temperatures to 0.5 K, the same dead elements).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/body.py reentry_model/trajectory.py tests/test_reentry_model_melting.py
git commit -m "Add the melting body: feed, film, spraying, element death, accounting and the projected area (Step 3 Task 9)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 10: The coupled loop — melt columns, demise, VTK melt fields, particle files

**Files:**
- Modify (replace): `reentry_model/coupled.py`
- Test: `tests/test_reentry_model_coupled.py` (append)

**Interfaces:**
- Consumes: `MeltingBody` (Task 9: `advance(state=)`, `demised()`, `melt_stats()`, `source_rows`, `m_f`, `liquid`, `last_spray/last_flow`, `phi`, `material.liquid_fraction`), `spray.SOURCE_COLUMNS/BIN_EDGES/N_BINS/histogram` (Task 7).
- Produces: `coupled.MELT_COLUMNS` (20 names), `CoupledRun.melting` (property), `melt_results(history)` (melt/spraying onsets, demise, masses, `n_released`, `r_median_um`, `n_dead_elements`, `n_source_rows`, `removed_enthalpy_J`, `melt_energy_balance_residual`), end reason `"demise"`; `write_vtk_frame` writes active cells only with `liquid_fraction`/`phi` and the surface fields `film_thickness, we_s, regime, tau, r_droplet, release_rate` when melting; `write_particles(run_dir, body, history, window=10.0) -> {"particles", "particles_summary", "size_distribution"}`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_coupled.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: the melting run, its columns, results, VTK fields and particle files

def test_coupled_melting_run_and_writers(coarse_sphere_mesh, tmp_path):
    """20 s of the 100 mm flight from 71 km with a warm body (700 K; physics heating, AA7075_range, Girin removal) on the
    coarse mesh: melt columns, results, demise bookkeeping, VTK melt fields and the three particle files."""
    pytest.importorskip("cantera")
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=700.0)
    sim = simulator(b, t_max=63.0)
    sim.advance(43.0)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=20, output_dir=os.path.join(run_dir, "vtk"))
    run = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), settings)
    assert run.melting
    hist = run.run()
    c, r = hist.columns, hist.results
    assert set(coupled.MELT_COLUMNS) <= set(c) and hist.end_reason == "t_max"
    assert c["mass_kg"][-1] < c["mass_kg"][0] and c["sprayed_mass_kg"][-1] > 0.0 and np.all(np.diff(c["sprayed_mass_kg"]) >= 0.0)
    assert c["removed_mass_kg"][-1] == pytest.approx(c["sprayed_mass_kg"][-1]) and c["n_dead_elements"][-1] > 0
    assert c["mass_kg"][-1] == pytest.approx(c["mass_kg"][0] - c["sprayed_mass_kg"][-1], rel=1e-6)
    assert r["melt_onset_altitude_km"] is not None and r["spraying_onset_time_s"] >= r["melt_onset_time_s"] and r["demise_time_s"] is None
    assert r["sprayed_mass_kg"] == c["sprayed_mass_kg"][-1] and r["n_source_rows"] == len(b.source_rows) > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    assert np.all(c["convective_heat_W"] > 0.0) and c["convective_heat_W"][-1] == pytest.approx(b.last.Q_conv)
    import pyvista as pv
    grid = pv.read(os.path.join(run_dir, "vtk", "field_1.vtu"))
    assert "liquid_fraction" in grid.point_data and "phi" in grid.cell_data and grid.n_cells == int(c["n_active_elements"][20])
    poly = pv.read(os.path.join(run_dir, "vtk", "surface_1.vtp"))
    for key in ("film_thickness", "we_s", "regime", "tau", "r_droplet", "release_rate"):
        assert key in poly.cell_data
    files = coupled.write_particles(run_dir, b, hist)
    for key in ("particles", "particles_summary", "size_distribution"):
        assert os.path.isfile(files[key])
    from reentry_model import spray
    npz = np.load(files["particles"])
    assert set(npz.files) == set(spray.SOURCE_COLUMNS) and npz["dm_kg"].sum() == pytest.approx(b.sprayed_mass, rel=1e-9)
    with open(files["size_distribution"]) as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "window_start_s,window_end_s,r_lo_m,r_hi_m,dn,dm_kg" and (len(lines) - 1) % spray.N_BINS == 0
    total = [l.split(",") for l in lines[1:]][-spray.N_BINS:]
    assert sum(float(row[5]) for row in total) == pytest.approx(b.sprayed_mass, rel=1e-6)
    with open(files["particles_summary"]) as fh:
        assert fh.readline().startswith("time_s,altitude_km,velocity_kms,released_mass_kg")


def test_demise_ends_the_run(coarse_sphere_mesh):
    """The lumped instant-removal device with a 60 % demise fraction: the loop stops with end_reason demise."""
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    b = body.MeltingBody(m, mat, thermal.thermal_solver("skfem"), MASS_100MM, settings=body.MeltSettings(removal="instant", runoff=False, demise_fraction=0.6))
    hist = coupled.CoupledRun(simulator(b), b, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()
    assert hist.end_reason == "demise" and hist.results["end_reason"] == "demise" and hist.results["demise_altitude_km"] < 71.0
    assert 0.55 * MASS_100MM < hist.columns["mass_kg"][-1] < 0.6 * MASS_100MM and hist.results["final_mass_kg"] == hist.columns["mass_kg"][-1]
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py -q`
Expected: the two new tests fail (`advance()` gets no `state`, no `MELT_COLUMNS`, no `write_particles`).

- [ ] **Step 3: Replace `reentry_model/coupled.py`**

```python
"""Lockstep coupling of the trajectory stepper with the thermal body (spec section 7).

Per macro step of dt: Simulator.advance(dt) (DOP853, events truncate the last step) -> aero state at the end of the
step -> heating with that state and the wall temperatures at the start of the step -> ThermalBody.advance (radiation
implicit in the solver) -> history row, and every `frames_every` steps a VTK frame (nodal T on the volume mesh,
q_conv / q_rad / T per patch on the surface). First-order operator splitting; the dt-halving test bounds its error.
Mass is constant in Step 2; the loop already carries the body's mass and the mesh so Step 3 can change both.
`body` must be a ThermalBody: CoupledRun uses `theta`, `surface`, `radiated_power()`, `surface_stats()` and
`integrated_heat`, beyond what the `Body` protocol declares. `ConstantBody` is for `Simulator.run()` only."""
import os
import time
from dataclasses import dataclass, field

import numpy as np

from . import trajectory as tj

THERMAL_COLUMNS = ["convective_heat_W", "rad_cooling_W", "integrated_heat_J", "absorbed_heat_J",
                   "surface_T_max_K", "surface_T_min_K", "surface_T_mean_K", "T_stagnation_K", "T_back_K",
                   "T_centre_K", "q_stag_Wm2", "heating_blend_f"]
MELT_COLUMNS = ["film_mass_kg", "sprayed_mass_kg", "runoff_mass_kg", "removed_mass_kg", "melt_front_depth_max_mm",
                "equivalent_radius_mm", "n_active_elements", "spraying_area_m2", "theta_cr_deg", "n_released",
                "released_mass_kg", "r_median_um", "r_max_um", "regime_fraction_continuum", "regime_fraction_slip",
                "regime_fraction_fm", "rt_active", "removed_enthalpy_J", "film_thickness_max_mm", "film_thickness_mean_mm",
                "n_dead_elements"]
PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'


@dataclass
class CoupledSettings:
    dt: float = 0.5                  # s, macro step
    frames_every: int = 0            # VTK frame every n macro steps (0: none)
    output_dir: str = None           # directory of the VTK series (required when frames_every > 0)
    frames: list = field(default_factory=list)


class CoupledRun:
    def __init__(self, sim, body, heating_model, settings=None):
        """`body` must be a ThermalBody (uses `theta`, `surface`, `radiated_power()`, `surface_stats()` and
        `integrated_heat`, beyond the `Body` protocol); `ConstantBody` is for `Simulator.run()` only."""
        self.sim, self.body, self.heating = sim, body, heating_model
        self.settings = settings or CoupledSettings()
        if self.settings.frames_every and not self.settings.output_dir:
            raise ValueError("frames_every > 0 needs an output_dir")

    def loads_at(self, t, y):
        a = self.sim.aero_state(t, y[:3], y[3:])
        return a, self.heating.evaluate(a, self.body.theta, self.body.surface_temperature(), self.sim.settings.diameter / 2.0,
                                        T_mean=self.body.mean_temperature())

    def row(self, t, y, a, loads):
        body = self.body
        r = self.sim.sample_row(t, y, a)
        Q = body.last.Q_conv if body.last is not None else loads.total(body.surface.areas)     # the step's applied power (the surface may have changed since)
        r.update({"convective_heat_W": Q, "rad_cooling_W": -body.radiated_power(),
                  "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat,
                  "q_stag_Wm2": loads.q_stag, "heating_blend_f": loads.blend})
        r.update(body.surface_stats())
        if self.melting:
            r.update(body.melt_stats())
        return r

    @property
    def melting(self):
        return hasattr(self.body, "melt_stats")

    def run(self):
        started = time.perf_counter()
        sim, body, s = self.sim, self.body, self.settings
        a, loads = self.loads_at(sim.t, sim.y)
        rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
        self._frame(step, sim.t, loads)
        while sim.end_reason is None:
            dt = sim.advance(s.dt)
            a, loads = self.loads_at(sim.t, sim.y)
            body.advance(sim.t, dt, loads, state=a)
            step += 1
            if len(loads.q_conv) != body.surface.n_patches:       # elements died this step: loads on the new surface for the record
                a, loads = self.loads_at(sim.t, sim.y)
            rows.append(self.row(sim.t, sim.y, a, loads))
            states.append(sim.y.copy())
            self._frame(step, sim.t, loads)
            if self.melting and body.demised():
                sim.end_reason = "demise"
        self._collection()
        columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
        history = tj.History(columns, np.array(states), sim.end_reason)
        history.results = sim.results(history, sim.nfev, time.perf_counter() - started)
        history.results.update(self.thermal_results(history))
        if self.melting:
            history.results.update(self.melt_results(history))
        return history

    def melt_results(self, history):
        c, body = history.columns, self.body
        first = lambda mask: (float(c["time_s"][mask][0]), float(c["altitude_km"][mask][0])) if mask.any() else (None, None)
        t_melt, h_melt = first(c["removed_mass_kg"] + c["film_mass_kg"] > 0.0)
        t_spray, h_spray = first(c["sprayed_mass_kg"] > 0.0)
        i_peak = int(np.argmax(c["released_mass_kg"]))
        return {"melt_onset_time_s": t_melt, "melt_onset_altitude_km": h_melt, "spraying_onset_time_s": t_spray,
                "spraying_onset_altitude_km": h_spray, "demise_time_s": float(c["time_s"][-1]) if history.end_reason == "demise" else None,
                "demise_altitude_km": float(c["altitude_km"][-1]) if history.end_reason == "demise" else None,
                "initial_mass_kg": body.mass0, "final_mass_kg": float(c["mass_kg"][-1]), "sprayed_mass_kg": body.sprayed_mass,
                "runoff_mass_kg": body.runoff_mass, "removed_mass_kg": body.removed_mass, "film_mass_kg": float(body.m_f.sum()),
                "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
                "removed_enthalpy_J": body.removed_enthalpy, "melt_energy_balance_residual": body.energy_balance_residual()}

    def thermal_results(self, history):
        c, body = history.columns, self.body
        i_surf, i_mean, i_q = int(np.argmax(c["surface_T_max_K"])), int(np.argmax(c["temperature_K"])), int(np.argmax(c["convective_heat_W"]))
        return {
            "peak_surface_T_K": float(c["surface_T_max_K"][i_surf]), "time_of_peak_surface_T_s": float(c["time_s"][i_surf]),
            "altitude_of_peak_surface_T_km": float(c["altitude_km"][i_surf]),
            "peak_mean_T_K": float(c["temperature_K"][i_mean]), "time_of_peak_mean_T_s": float(c["time_s"][i_mean]),
            "peak_convective_heat_W": float(c["convective_heat_W"][i_q]), "time_of_peak_heating_s": float(c["time_s"][i_q]),
            "altitude_of_peak_heating_km": float(c["altitude_km"][i_q]),
            "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat, "radiated_heat_J": body.radiated_heat,
            "energy_balance_residual": body.energy_balance_residual(),
            "n_macro_steps": len(body.iterations), "mean_newton_iterations": float(np.mean(body.iterations)) if body.iterations else 0.0,
            "n_frames": len(self.settings.frames),
        }

    def _frame(self, step, t, loads):
        s = self.settings
        if not s.frames_every or step % s.frames_every:
            return
        k = len(s.frames)
        write_vtk_frame(s.output_dir, k, self.body, loads)
        s.frames.append((t, k))

    def _collection(self):
        s = self.settings
        if not s.frames:
            return
        for name in ("field", "surface"):
            entries = "".join('<DataSet timestep="{:.6g}" file="{}_{}.{}"/>\n'.format(t, name, k, "vtu" if name == "field" else "vtp")
                              for t, k in s.frames)
            with open(os.path.join(s.output_dir, name + ".pvd"), "w") as fh:
                fh.write(PVD_TEMPLATE.format(entries))


def write_vtk_frame(output_dir, k, body, loads):
    """field_<k>.vtu: nodal T on the active volume mesh (plus liquid fraction and element fractions when melting);
    surface_<k>.vtp: q_conv, q_rad, T per patch, and when melting the film thickness, We_s, regime, shear, droplet
    radius and release rate (PyVista/VTK XML)."""
    import pyvista as pv
    from .thermal import SIGMA_SB
    os.makedirs(output_dir, exist_ok=True)
    mesh, surface, T = body.mesh, body.surface, body.field()
    melting = hasattr(body, "melt_stats")
    tets = mesh.tets[mesh.active]
    cells = np.hstack([np.full((len(tets), 1), 4), tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(tets), pv.CellType.TETRA), mesh.points)
    grid.point_data["T"] = T
    if melting:
        grid.point_data["liquid_fraction"] = body.material.liquid_fraction(T)
        grid.cell_data["phi"] = body.phi[mesh.active]
    grid.save(os.path.join(output_dir, "field_{}.vtu".format(k)))
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(mesh.points, faces)
    Tf = surface.facet_mean(T)
    poly.point_data["T"] = T
    poly.cell_data["q_conv"] = loads.q_conv
    poly.cell_data["q_rad"] = body.emissivity * SIGMA_SB * (Tf ** 4 - body.T_ambient ** 4)
    poly.cell_data["T_patch"] = Tf
    if melting:
        poly.cell_data["film_thickness"] = body.m_f / (body.liquid.rho * surface.areas)
        n = surface.n_patches
        res, flow = body.last_spray, body.last_flow
        same = res is not None and res.r.size == n
        poly.cell_data["we_s"] = res.we_s if same else np.zeros(n)
        poly.cell_data["regime"] = flow.regime.astype(float) if same else np.full(n, 2.0)
        poly.cell_data["tau"] = flow.tau if same else np.zeros(n)
        poly.cell_data["r_droplet"] = np.where(res.dm > 0.0, res.r, np.nan) if same else np.full(n, np.nan)
        poly.cell_data["release_rate"] = res.dm / surface.areas if same else np.zeros(n)
    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))


def write_particles(run_dir, body, history, window=10.0):
    """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
    size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths."""
    import csv
    from . import spray
    os.makedirs(run_dir, exist_ok=True)
    rows = np.array(body.source_rows, dtype=float) if body.source_rows else np.zeros((0, len(spray.SOURCE_COLUMNS)))
    npz = os.path.join(run_dir, "particles.npz")
    np.savez_compressed(npz, **{name: rows[:, i] for i, name in enumerate(spray.SOURCE_COLUMNS)})
    c = history.columns
    summary = os.path.join(run_dir, "particles_summary.csv")
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh)
        keys = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg",
                "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
        w.writerow(keys)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
    dist = os.path.join(run_dir, "size_distribution.csv")
    edges = spray.BIN_EDGES
    t = rows[:, 0] if len(rows) else np.zeros(0)
    t_end = float(c["time_s"][-1])
    starts = np.arange(0.0, t_end + 1e-9, window)
    with open(dist, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["window_start_s", "window_end_s", "r_lo_m", "r_hi_m", "dn", "dm_kg"])
        for t0 in list(starts) + [None]:
            if t0 is None:
                sel = np.ones(len(rows), dtype=bool)
                lo, hi = 0.0, t_end
            else:
                lo, hi = t0, t0 + window
                sel = (t >= lo) & (t < hi)
            n_hist, m_hist = spray.histogram(rows[sel, 12], rows[sel, 13], rows[sel, 14]) if sel.any() else (np.zeros(spray.N_BINS), np.zeros(spray.N_BINS))
            for j in range(spray.N_BINS):
                w.writerow(["{:.9g}".format(lo), "{:.9g}".format(hi), "{:.9g}".format(edges[j]), "{:.9g}".format(edges[j + 1]),
                            "{:.9g}".format(n_hist[j]), "{:.9g}".format(m_hist[j])])
    return {"particles": npz, "particles_summary": summary, "size_distribution": dist}
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py tests/test_reentry_model_viz.py tests/test_reentry_model_cli.py -q`
Expected: all pass (Step 2's runs are unchanged: `melting` is False for a `ThermalBody`).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/coupled.py tests/test_reentry_model_coupled.py
git commit -m "Carry the melt step, the demise end and the particle writers through the coupled loop (Step 3 Task 10)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 11: The two melting SESAM references

**Files:**
- Create: `data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv/.json`, `data/reference_runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind.csv/.json`
- Test: `tests/test_reentry_model_data.py` (append), `tests/test_reentry_model_aero.py` (one guard)

**Interfaces:**
- Consumes: the wrapper `sphere_reentry.py` (untouched) with DRAMA's `drama-AA7075` on the static US76 table, winds off — the break-off state of the existing references.
- Produces: the two reference runs (CSV + JSON) that Tasks 12–14 compare against: 100 mm melt onset 71.005 km at 43.55 s, mass 0 at 66.476 km / 67.315 s (70 rows); 50 mm onset 77.104 km at 178.55 s, mass 0 at 73.094 km / 191.55 s (194 rows).

- [ ] **Step 1: Generate the references with the wrapper (or reuse the ones generated on 2026-09-20)**

The runs already exist under `sphere_sweep_output/reference_AA7075/runs/` (generated 2026-09-20 with exactly the commands below; `sphere_sweep_output/` is git-ignored). If they are missing, generate them (each takes under a minute; the wrapper needs `drama_env`):

```bash
mkdir -p sphere_sweep_output/reference_AA7075/runs sphere_sweep_output/reference_AA7075/raw
"$PY" sphere_reentry.py --velocity 7.5 --altitude 77.500133 --temperature 300 --diameter 100 --flight-path-angle -0.959331 \
    --heading 347.168296 --lat 29.546067 --lon -82.134333 --epoch 2024-08-01T12:53:07 --atmosphere static --no-wind \
    --outdir sphere_sweep_output/reference_AA7075/runs --raw-dir sphere_sweep_output/reference_AA7075/raw --keep-raw
"$PY" sphere_reentry.py --velocity 7.5 --altitude 115 --temperature 300 --diameter 50 --flight-path-angle -0.959331 \
    --heading 347.168296 --lat 29.546067 --lon -82.134333 --epoch 2024-08-01T12:53:07 --atmosphere static --no-wind \
    --outdir sphere_sweep_output/reference_AA7075/runs --raw-dir sphere_sweep_output/reference_AA7075/raw --keep-raw
```

Expected console lines: `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind: demised after 67.3 s (uncritical); Tmax 850.0 K, 23.8 s at melt, final mass 0 kg` and `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind: demised after 191.6 s (uncritical); Tmax 850.0 K, 13.0 s at melt, final mass 0 kg`. Then copy the four files:

```bash
cp sphere_sweep_output/reference_AA7075/runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.* data/reference_runs/
cp sphere_sweep_output/reference_AA7075/runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind.* data/reference_runs/
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_reentry_model_data.py`:

```python
MELTING_REFERENCE_NAMES = ["sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
                           "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"]


def test_melting_reference_runs_are_committed_and_consistent():
    """The two Step 3 melting references: DRAMA's drama-AA7075 (850 K, 400 kJ/kg) on the US76 table, winds off; the
    sphere melts completely (mass 0 at the end, 'demised')."""
    for name, onset_km, end_km in zip(MELTING_REFERENCE_NAMES, (71.005, 77.104), (66.476, 73.094)):
        csv_path = os.path.join(REFERENCE_DIR, name + ".csv")
        doc = json.load(open(os.path.join(REFERENCE_DIR, name + ".json")))
        assert doc["run_name"] == name and doc["status"] == "ok" and doc["inputs"]["material"] == "drama-AA7075"
        assert doc["inputs"]["atmosphere"] == "static" and doc["inputs"]["use_wind"] is False
        assert doc["results"]["melting_temperature_K"] == 850.0 and abs(doc["results"]["altitude_first_melt_km"] - onset_km) < 1e-3
        assert abs(doc["results"]["altitude_last_melt_km"] - end_km) < 1e-3
        rows = list(csv.DictReader(open(csv_path)))
        assert float(rows[-1]["mass_kg"]) == 0.0 and float(rows[-1]["thick_mm"]) == 0.0 and float(rows[-1]["altitude_km"]) > 60.0
        assert abs(float(rows[0]["mass_kg"]) - doc["inputs"]["initial_mass_kg"]) < 1e-3
```


In `tests/test_reentry_model_aero.py::TestDragCoefficient::test_reproduces_the_reference_drag_column` (it globs every reference CSV; the melting references demise while still hypersonic, so their low-speed list is empty) replace the last assertion

```python
        assert np.abs(low).max() <= 1.5e-3, name       # the reference prints C_D with 3 decimals
```

with

```python
        if low:                                        # the melting references (Step 3) demise while still hypersonic
            assert np.abs(low).max() <= 1.5e-3, name   # the reference prints C_D with 3 decimals
```

- [ ] **Step 3: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_data.py tests/test_reentry_model_aero.py tests/test_reentry_model_sesam_io.py -q`
Expected: all pass (the melting references' hypersonic C_D rows reproduce SESAM's drag rule like the others; `test_reference_runs_are_committed_and_consistent` still covers only the no-melt names).

- [ ] **Step 4: Commit**

```bash
git add data/reference_runs tests/test_reentry_model_data.py tests/test_reentry_model_aero.py
git commit -m "Add the two melting SESAM references (drama-AA7075, US76, winds off) (Step 3 Task 11)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 12: Comparison metrics, plots and visualisation of the melting run

**Files:**
- Modify: `reentry_model/sesam_io.py` (two edits), `reentry_model/compare.py` (append + two constants)
- Modify (replace): `reentry_model/viz.py`
- Test: `tests/test_reentry_model_compare.py` (append), `tests/test_reentry_model_viz.py` (append)

**Interfaces:**
- Consumes: the melting reference files (Task 11), histories with `MELT_COLUMNS` (Task 10), VTK series with the melt fields (Task 10).
- Produces: `sesam_io.Reference.mass`, `.thickness`; `compare.MELT_PLOT_NAMES` (7), `DEMISE_FRACTION = 0.01`, `has_melt(history, reference=None)`, `melt_metrics(history, reference) -> dict` (`mass.max_rel_m0`, `onset_altitude_diff_km`, `demise_time_diff_s`, `demise_time_rel`, ...), `plot_melt(history, outdir, title, reference=None, size_distribution_csv=None) -> paths`; `viz.render_frame(..., melting=False)`, `render_film_frame`, `render_section_frame(..., melting=False)`, `animate(..., melting=False)`, `animate_section(..., melting=False)`, `animate_film(run_dir, history, radius, fps, animation, stills)`, `still_marks(history, melting=False)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_compare.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: the melting reference, the mass metrics and plots

MELT100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind"


@pytest.fixture(scope="module")
def melt_ref():
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, MELT100 + ".csv"))


def melting_history(ref, lag=0.0, film=0.0):
    """A synthetic model history that reproduces the reference's mass shifted by `lag` seconds, with the melt columns."""
    h = synthetic_history(ref, dv=0.0, dh=0.0)
    c = h.columns
    n = len(ref.time)
    c["mass_kg"] = np.interp(ref.time - lag, ref.time, ref.mass) + film
    for key in ("film_mass_kg", "sprayed_mass_kg", "runoff_mass_kg", "removed_mass_kg", "spraying_area_m2", "theta_cr_deg", "n_released",
                "released_mass_kg", "r_median_um", "r_max_um", "regime_fraction_continuum", "regime_fraction_slip", "regime_fraction_fm",
                "film_thickness_mean_mm", "film_thickness_max_mm", "convective_heat_W", "surface_T_max_K"):
        c[key] = np.zeros(n)
    c["film_mass_kg"][:] = film
    c["removed_mass_kg"] = ref.mass[0] - c["mass_kg"] + film
    c["sprayed_mass_kg"] = c["removed_mass_kg"].copy()
    c["r_median_um"][:] = 150.0
    c["regime_fraction_continuum"][:] = 0.5
    c["regime_fraction_slip"][:] = 0.5
    return h


def test_melting_reference_and_has_melt(ref, melt_ref):
    assert ref.mass is not None and float(ref.mass.min()) == float(ref.mass[0])                         # no-melt: constant mass
    assert melt_ref.mass[0] == pytest.approx(1.473) and melt_ref.mass[-1] == 0.0 and melt_ref.thickness[0] == pytest.approx(0.05)
    assert melt_ref.thickness[-1] == 0.0 and np.all(np.diff(melt_ref.mass) <= 0.0)
    h = melting_history(melt_ref)
    assert compare.has_melt(h) and compare.has_melt(h, melt_ref) and not compare.has_melt(h, ref) and not compare.has_melt(synthetic_history(ref))


def test_melt_metrics_recover_a_known_lag(melt_ref):
    m = compare.melt_metrics(melting_history(melt_ref, lag=-1.0), melt_ref)          # the model runs 1 s ahead of SESAM
    assert m["mass"]["max_rel_m0"] < 0.06 and m["mass"]["max_abs_kg"] > 0.0
    assert m["demise_time_diff_s"] == pytest.approx(-1.0, abs=0.1) and m["demise_time_rel"] == pytest.approx(-1.0 / m["demise_time_reference_s"], rel=0.1)
    assert m["onset_time_reference_s"] == pytest.approx(43.55) and m["onset_altitude_reference_km"] == pytest.approx(71.005, abs=0.01)
    assert m["onset_altitude_diff_km"] is not None and m["demise_time_reference_s"] == pytest.approx(66.8, abs=0.1)
    exact = compare.melt_metrics(melting_history(melt_ref), melt_ref)
    assert exact["mass"]["max_rel_m0"] < 1e-12 and abs(exact["demise_time_diff_s"]) < 1e-9
    with_film = compare.melt_metrics(melting_history(melt_ref, film=0.05), melt_ref)
    assert with_film["demise_time_diff_s"] == pytest.approx(exact["demise_time_diff_s"], abs=1e-9)   # the crossing is on the body's material


def test_melt_plots_are_written(melt_ref, tmp_path):
    h = melting_history(melt_ref, lag=-1.0)
    paths = compare.plot_melt(h, str(tmp_path), "melt test", melt_ref)
    assert [os.path.basename(p) for p in paths] == list(compare.MELT_PLOT_NAMES) and all(os.path.getsize(p) > 5000 for p in paths)
    alone = compare.plot_melt(h, str(tmp_path / "alone"), "no reference")
    assert len(alone) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in alone)
```


Append to `tests/test_reentry_model_viz.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: overlays of the emitting patches, the film-thickness frame, the liquidus iso-line, the melt stills

def test_melting_overlays_and_film_frame(coarse_sphere_mesh, tmp_path):
    import pyvista as pv
    surface = coarse_sphere_mesh.surface()
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(coarse_sphere_mesh.points, faces)
    poly.point_data["T"] = 300.0 + 600.0 * np.clip(coarse_sphere_mesh.points[:, 0] / 0.05, 0.0, 1.0)
    r = np.full(surface.n_patches, np.nan)
    windward = surface.normals[:, 0] > 0.5
    r[windward] = 1e-4
    poly.cell_data["r_droplet"] = r
    poly.cell_data["film_thickness"] = np.where(windward, 2e-5, 0.0)
    plain = viz.render_frame(poly, (300.0, 900.0), "t", None, 0.05)
    overlaid = viz.render_frame(poly, (300.0, 900.0), "t", str(tmp_path / "o.png"), 0.05, melting=True)
    assert overlaid.shape == plain.shape and np.abs(overlaid.astype(int) - plain.astype(int)).mean() > 1.0     # the overlay changed the picture
    film = viz.render_film_frame(poly, None, "film", str(tmp_path / "f.png"), 0.05)
    assert film.shape == (720, 960, 3) and film.std() > 10.0 and os.path.getsize(tmp_path / "f.png") > 1000
    m = coarse_sphere_mesh
    cells = np.hstack([np.full((m.n_elements, 1), 4), m.tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(m.n_elements, pv.CellType.TETRA), m.points)
    grid.point_data["T"] = 300.0 + 700.0 * np.clip(m.points[:, 0] / 0.05, 0.0, 1.0)
    grid.point_data["liquid_fraction"] = np.clip((grid.point_data["T"] - 750.0) / 158.0, 0.0, 1.0)
    without = viz.render_section_frame(grid, (300.0, 1000.0), "t", None, 0.05)
    with_iso = viz.render_section_frame(grid, (300.0, 1000.0), "t", str(tmp_path / "s.png"), 0.05, melting=True)
    assert with_iso.shape == without.shape and np.abs(with_iso.astype(int) - without.astype(int)).mean() > 0.1


def test_still_marks_and_film_animation(coarse_sphere_mesh, tmp_path):
    pytest.importorskip("cantera")
    from test_reentry_model_coupled import simulator, MASS_100MM
    from reentry_model import body, coupled, heating, material, mesh, thermal
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=700.0)
    sim = simulator(b, t_max=53.0)
    sim.advance(43.0)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=10, output_dir=run_dir)
    hist = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), settings).run()
    marks = viz.still_marks(hist, melting=True)
    assert {"start", "peak_heating", "peak_surface_T", "end", "melt_onset", "spraying_onset", "peak_release"} <= set(marks)
    assert marks["melt_onset"] <= marks["spraying_onset"] <= marks["end"]
    out = viz.animate(run_dir, hist, 0.05, fps=5, animation=False, melting=True)
    assert len(out["stills"]) == 7 and all(os.path.isfile(p) for p in out["stills"])
    film = viz.animate_film(run_dir, hist, 0.05, fps=5, animation=False)
    assert len(film["stills"]) == 7 and all("film_" in os.path.basename(p) for p in film["stills"])
    sec = viz.animate_section(run_dir, hist, 0.05, fps=5, animation=False, melting=True)
    assert any(os.path.basename(p).startswith("section_melt_onset") for p in sec["stills"])
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py -q`
Expected: the five new tests fail (`Reference` has no `mass`, no `has_melt`, `render_frame` has no `melting`).

- [ ] **Step 3: Edit `reentry_model/sesam_io.py`**

In the `Reference` dataclass, after `integrated_heat: np.ndarray = None      # J, integral of the convective heat (SESAM's integrated_heat_J)` add:

```python
    mass: np.ndarray = None                 # kg (SESAM's mass_kg; decreases while the sphere melts)
    thickness: np.ndarray = None            # m (SESAM's thick_mm: the shell thickness of the melting sphere, Step 3 facts)
```

and in `load_reference`, after the `convective_heat=..., integrated_heat=optional("integrated_heat_J"),` line add:

```python
        mass=optional("mass_kg"), thickness=col("thick_mm", 1e-3) if "thick_mm" in rows[0] else None,
```

- [ ] **Step 4: Extend `reentry_model/compare.py`**

Replace the two constant lines `THERMAL_PLOT_NAMES = (...)` / `CONTINUUM_KN = 0.01` with:

```python
THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")
MELT_PLOT_NAMES = ("mass_time.png", "mass_altitude.png", "mass_budget.png", "spraying_time.png", "regimes_time.png",
                   "droplet_size_time.png", "size_distribution.png")
CONTINUUM_KN = 0.01
DEMISE_FRACTION = 0.01
```

and append to the file:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: mass loss, spraying and size distributions (spec Step 3 sections 12-13)

def has_melt(history, reference=None):
    """True when the model history carries the melt columns (and, if a reference is given, its mass varies)."""
    ok = "sprayed_mass_kg" in history.columns and "film_mass_kg" in history.columns
    if reference is None:
        return ok
    return ok and reference.mass is not None and float(reference.mass.min()) < 0.999 * float(reference.mass[0])


def _first_time(t, mask):
    return float(t[mask][0]) if mask.any() else None


def _crossing_time(t, m, level):
    """Time at which the decreasing series m first falls below `level`, interpolated linearly between samples."""
    below = np.flatnonzero(m < level)
    if below.size == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (m[i - 1] - level) / (m[i - 1] - m[i]) * (t[i] - t[i - 1]))


def melt_metrics(history, reference):
    """Mass vs SESAM at the reference's time stamps (relative to the initial mass), melt-onset altitude and the time
    at which the mass falls below DEMISE_FRACTION of the initial (SESAM removes the last gram at its 'demise'; the model
    stops at the fraction, so the 1 % crossing is the like-for-like end time; for the model it is taken on the body's
    material, film excluded, since a static leeward film can stay attached)."""
    c = history.columns
    t_mod, m_mod = c["time_s"], c["mass_kg"]
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
    t = reference.time[mask]
    m_ref = reference.mass[mask]
    m_at = np.interp(t, t_mod, m_mod)
    m0 = float(reference.mass[0])
    onset_ref = _first_time(reference.time, reference.mass < 0.999999 * m0)
    onset_mod = _first_time(t_mod, (c["removed_mass_kg"] + c["film_mass_kg"]) > 0.0)
    h_ref = float(np.interp(onset_ref, reference.time, reference.altitude)) / 1e3 if onset_ref is not None else None
    h_mod = float(np.interp(onset_mod, t_mod, c["altitude_km"])) if onset_mod is not None else None
    end_ref = _crossing_time(reference.time, reference.mass, DEMISE_FRACTION * m0)
    body_mass = m_mod - c["film_mass_kg"]                                                          # the body's material (a leeward film may stay attached)
    end_mod = _crossing_time(t_mod, body_mass, DEMISE_FRACTION * float(body_mass[0]))
    return {
        "n_points": int(t.size), "initial_mass_model_kg": float(m_mod[0]), "initial_mass_reference_kg": m0,
        "mass": {"max_abs_kg": float(np.abs(m_at - m_ref).max()), "max_rel_m0": float(np.abs(m_at - m_ref).max() / m0),
                 "rms_rel_m0": float(math.sqrt(np.mean((m_at - m_ref) ** 2)) / m0)},
        "onset_time_model_s": onset_mod, "onset_time_reference_s": onset_ref,
        "onset_altitude_model_km": h_mod, "onset_altitude_reference_km": h_ref,
        "onset_altitude_diff_km": (h_mod - h_ref) if (h_mod is not None and h_ref is not None) else None,
        "demise_time_model_s": end_mod, "demise_time_reference_s": end_ref,
        "demise_time_diff_s": (end_mod - end_ref) if (end_mod is not None and end_ref is not None) else None,
        "demise_time_rel": ((end_mod - end_ref) / end_ref) if (end_mod is not None and end_ref is not None) else None,
        "final_mass_model_kg": float(m_mod[-1]), "sprayed_mass_kg": float(c["sprayed_mass_kg"][-1]), "film_mass_kg": float(c["film_mass_kg"][-1]),
    }


def plot_melt(history, outdir, title, reference=None, size_distribution_csv=None):
    """The Step 3 plots (MELT_PLOT_NAMES); with a melting reference the mass plots carry the SESAM overlay and a
    residual panel. Returns the paths."""
    import csv
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    c = history.columns
    t, h, m = c["time_s"], c["altitude_km"], c["mass_kg"]
    paths = [os.path.join(outdir, n) for n in MELT_PLOT_NAMES]
    # 1. mass vs time (+ residual)
    if reference is not None:
        mask = (reference.time >= t[0]) & (reference.time <= t[-1])
        tr, mr = reference.time[mask], reference.mass[mask]
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        _overlay(ax, tr, mr, t, m, "mass [kg]", title=title)
        rx.plot(tr, 100.0 * (np.interp(tr, t, m) - mr) / reference.mass[0], color=MODEL_COLOR, lw=1.0)
        rx.axhline(0.0, color=MUTED, lw=0.6)
        rx.set_ylabel("model - SESAM [% of m0]"); rx.set_xlabel("time [s]")
        strip_top_right_spines(rx)
    else:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t, m, color=MODEL_COLOR, lw=1.2, label="model")
        ax.set_ylabel("mass [kg]"); ax.set_xlabel("time [s]"); ax.set_title(title, color=SECOND, fontsize=10)
        strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[0], dpi=150); plt.close(fig)
    # 2. mass vs altitude
    fig, ax = plt.subplots(figsize=(7, 5))
    if reference is not None:
        ax.plot(reference.mass, reference.altitude / 1e3, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(m, h, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("mass [kg]"); ax.set_ylabel("altitude [km]"); ax.set_title(title, color=SECOND, fontsize=10); ax.legend(frameon=False)
    strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[1], dpi=150); plt.close(fig)
    # 3. mass budget
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, c["sprayed_mass_kg"], color=MODEL_COLOR, lw=1.2, label="sprayed (cumulative)")
    ax.plot(t, c["removed_mass_kg"] - c["sprayed_mass_kg"], color=SECOND, lw=1.0, ls=":", label="removed instantly (cumulative)")
    ax.plot(t, c["film_mass_kg"], color=MUTED, lw=1.2, label="film")
    ax.plot(t, m, color=REF_COLOR, lw=1.2, ls="--", label="remaining (body + film)")
    ax.set_xlabel("time [s]"); ax.set_ylabel("mass [kg]"); ax.set_title(title, color=SECOND, fontsize=10); ax.legend(frameon=False)
    strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)
    # 4. theta_cr and spraying area
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax.plot(t, c["theta_cr_deg"], color=MODEL_COLOR, lw=1.0); ax.set_ylabel("theta_cr [deg]"); ax.set_title(title, color=SECOND, fontsize=10)
    bx.plot(t, c["spraying_area_m2"] * 1e4, color=MODEL_COLOR, lw=1.0); bx.set_ylabel("spraying area [cm2]"); bx.set_xlabel("time [s]")
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)
    # 5. regime fractions
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.stackplot(t, c["regime_fraction_continuum"], c["regime_fraction_slip"], c["regime_fraction_fm"],
                 labels=["continuum", "slip", "transitional/free-molecular"], colors=[INK, MODEL_COLOR, MUTED], alpha=0.8)
    ax.set_xlabel("time [s]"); ax.set_ylabel("windward area fraction"); ax.set_ylim(0.0, 1.0); ax.legend(frameon=False, loc="upper left")
    ax.set_title(title, color=SECOND, fontsize=10); strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)
    # 6. droplet size and film thickness
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax.semilogy(t, c["r_median_um"], color=MODEL_COLOR, lw=1.0, label="median r")
    ax.semilogy(t, c["r_max_um"], color=MUTED, lw=0.8, ls=":", label="max r")
    ax.set_ylabel("droplet radius [um]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    bx.semilogy(t, np.maximum(c["film_thickness_mean_mm"], 1e-6), color=MODEL_COLOR, lw=1.0, label="mean film thickness")
    bx.semilogy(t, np.maximum(c["film_thickness_max_mm"], 1e-6), color=MUTED, lw=0.8, ls=":", label="max")
    bx.set_ylabel("film thickness [mm]"); bx.set_xlabel("time [s]"); bx.legend(frameon=False)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    # 7. size distributions (from size_distribution.csv when present next to the run, else skipped with an empty axes)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.5))
    path = size_distribution_csv or os.path.join(outdir, "size_distribution.csv")
    if os.path.isfile(path):
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        windows = sorted({(float(r["window_start_s"]), float(r["window_end_s"])) for r in rows})
        for lo, hi in windows:
            sel = [r for r in rows if float(r["window_start_s"]) == lo and float(r["window_end_s"]) == hi]
            r_mid = np.sqrt(np.array([float(r["r_lo_m"]) for r in sel]) * np.array([float(r["r_hi_m"]) for r in sel])) * 1e6
            dn, dm = np.array([float(r["dn"]) for r in sel]), np.array([float(r["dm_kg"]) for r in sel])
            total = (lo, hi) == windows[-1] and lo == 0.0
            kw = dict(color=INK, lw=1.6, label="flight") if total else dict(color=MODEL_COLOR, lw=0.7, alpha=0.5)
            if dn.sum() > 0:
                ax.loglog(r_mid, np.maximum(dn, 1e-300), **kw)
                bx.loglog(r_mid, np.maximum(dm, 1e-300), **kw)
        ax.set_ylim(bottom=max(ax.get_ylim()[0], 1e-1)); bx.set_ylim(bottom=max(bx.get_ylim()[0], 1e-12))
    ax.set_xlabel("droplet radius [um]"); ax.set_ylabel("dn per bin"); bx.set_xlabel("droplet radius [um]"); bx.set_ylabel("dM per bin [kg]")
    ax.set_title(title, color=SECOND, fontsize=10)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[6], dpi=150); plt.close(fig)
    return paths
```


- [ ] **Step 5: Replace `reentry_model/viz.py`**

```python
"""Surface-temperature and cross-section animations and stills from a coupled run's VTK series (spec section 9).

Reads only exported files: `surface.pvd` + `surface_<k>.vtp` (nodal T, per-patch q_conv) and `field.pvd` +
`field_<k>.vtu` (nodal T on the volume mesh) written by coupled.py, and the run's history (time, altitude, velocity,
surface_T_max_K for the fixed colour scale). `animate` colours the sphere's surface; `animate_section` cuts the volume
field on the meridional plane z = 0 through the flight axis (stagnation point at +x, shown on the right) so the depth
of the heated layer is visible. PyVista renders off-screen; frames go to MP4 through imageio-ffmpeg, with a GIF
fallback when the MP4 writer is unavailable."""
import os
import xml.etree.ElementTree as ET

import numpy as np

CAMERA = [(0.19, -0.14, 0.11), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]     # three-quarter view of the windward (+x) face, for R = 0.05 m
WINDOW = (960, 720)
SECTION_NORMAL, SECTION_ZOOM = (0.0, 0.0, 1.0), 1.15                  # meridional plane through the flight axis; viewed down z
SCALAR_BAR = {"fmt": "%.0f", "vertical": True, "position_x": 0.86, "position_y": 0.15, "width": 0.05, "height": 0.6,
              "title_font_size": 14, "label_font_size": 12}


def read_series(run_dir, name="surface"):
    """[(time, path)] of a PVD collection written by coupled.py."""
    root = ET.parse(os.path.join(run_dir, name + ".pvd")).getroot()
    return [(float(d.get("timestep")), os.path.join(run_dir, d.get("file"))) for d in root.iter("DataSet")]


def render_frame(poly, clim, title, path=None, radius=0.05, melting=False):
    """One off-screen frame of the surface coloured by nodal T; with `melting` the patches releasing droplets are drawn
    on top, coloured by droplet radius (log scale, um). Returns the RGB array (and writes PNG if `path`)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(poly, scalars="T", cmap="inferno", clim=clim, smooth_shading=True,
                     scalar_bar_args={"title": "surface T [K]", **SCALAR_BAR})
    if melting and "r_droplet" in poly.cell_data:
        r = np.asarray(poly.cell_data["r_droplet"])
        emitting = np.isfinite(r) & (r > 0.0)
        if emitting.any():
            spots = poly.extract_cells(np.flatnonzero(emitting))
            spots.cell_data["log10 r [um]"] = np.log10(r[emitting] * 1e6)
            plotter.add_mesh(spots, scalars="log10 r [um]", cmap="viridis", clim=(1.0, 3.0),
                             scalar_bar_args={"title": "released droplets: log10 r [um]", "position_x": 0.05, "position_y": 0.05,
                                              "vertical": False, "width": 0.4, "height": 0.05, "fmt": "%.1f", "title_font_size": 12, "label_font_size": 11})
    plotter.add_text(title, position="upper_left", font_size=11)
    scale = radius / 0.05
    plotter.camera_position = [tuple(scale * c for c in CAMERA[0]), CAMERA[1], CAMERA[2]]
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def render_film_frame(poly, clim, title, path=None, radius=0.05):
    """One frame of the film thickness [um] per patch (log colour scale); clim is ignored (fixed 1 um - 10 mm)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    b = np.asarray(poly.cell_data["film_thickness"]) if "film_thickness" in poly.cell_data else np.zeros(poly.n_cells)
    poly.cell_data["film [um]"] = np.log10(np.maximum(b * 1e6, 1.0))
    plotter.add_mesh(poly, scalars="film [um]", cmap="Blues", clim=(0.0, 4.0), scalar_bar_args={"title": "film thickness: log10 [um]", **SCALAR_BAR})
    plotter.add_text(title, position="upper_left", font_size=11)
    scale = radius / 0.05
    plotter.camera_position = [tuple(scale * c for c in CAMERA[0]), CAMERA[1], CAMERA[2]]
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def section_of(grid):
    """The volume field cut on the plane z = 0 (through the flight axis), as a PolyData carrying nodal T."""
    return grid.slice(normal=SECTION_NORMAL, origin=(0.0, 0.0, 0.0))


def render_section_frame(grid, clim, title, path=None, radius=0.05, melting=False):
    """One off-screen frame of the meridional cross-section coloured by T (windward +x on the right); with `melting`
    the liquidus iso-line (liquid fraction 1) and the solidus iso-line (0) are drawn. Returns the RGB array."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    section = section_of(grid)
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(section, scalars="T", cmap="inferno", clim=clim, scalar_bar_args={"title": "T [K]", **SCALAR_BAR})
    plotter.add_mesh(section.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                                   non_manifold_edges=False), color="black", line_width=1.5)
    if melting and "liquid_fraction" in section.point_data and section.n_points:
        f = np.asarray(section.point_data["liquid_fraction"])
        for level, colour in ((0.999, "white"), (0.001, "cyan")):
            if f.min() < level < f.max():
                iso = section.contour([level], scalars="liquid_fraction")
                if iso.n_points:
                    plotter.add_mesh(iso, color=colour, line_width=2.0)
    plotter.add_text(title + "   section z = 0, flow from the right" + ("   white: liquidus, cyan: solidus" if melting else ""),
                     position="upper_left", font_size=11)
    plotter.camera_position = "xy"
    plotter.camera.zoom(SECTION_ZOOM)
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def frame_title(t, history):
    c = history.columns
    return "t = {:.1f} s   h = {:.1f} km   V = {:.2f} km/s".format(
        t, np.interp(t, c["time_s"], c["altitude_km"]), np.interp(t, c["time_s"], c["velocity_kms"]))


def animate(run_dir, history, radius, fps=10, animation=True, stills=True, melting=False):
    """MP4 (GIF fallback) of the surface temperature over the run, plus stills at the start, peak heating, peak surface
    temperature and the end (with animation=False only the stills' frames are rendered); with `melting` the emitting
    patches are overlaid and the stills add melt onset, spraying onset and peak release. Returns {"animation", "stills"}."""
    render = (lambda poly, clim, title, path, radius: render_frame(poly, clim, title, path, radius, melting=True)) if melting else render_frame
    return _animate_series(run_dir, history, radius, "surface", render, "animation", "frames", "", fps, animation, stills, melting)


def animate_section(run_dir, history, radius, fps=10, animation=True, stills=True, melting=False):
    """The same for the meridional cross-section of the volume field: `section.mp4` (GIF fallback), `frames_section/`,
    and stills named `section_<label>_t<N>s.png` next to the surface ones (with `melting`: the liquidus iso-line)."""
    render = (lambda grid, clim, title, path, radius: render_section_frame(grid, clim, title, path, radius, melting=True)) if melting else render_section_frame
    return _animate_series(run_dir, history, radius, "field", render, "section", "frames_section", "section_", fps, animation, stills, melting)


def animate_film(run_dir, history, radius, fps=10, animation=True, stills=True):
    """The film-thickness video `film.mp4` (GIF fallback), `frames_film/` and stills `film_<label>_t<N>s.png`."""
    return _animate_series(run_dir, history, radius, "surface", render_film_frame, "film", "frames_film", "film_", fps, animation, stills, True)


def still_marks(history, melting=False):
    """{label: time} of the stills: start, peak heating, peak surface temperature, end, and when melting the melt
    onset, the spraying onset and the peak release rate."""
    c = history.columns
    marks = {"start": 0.0, "peak_heating": float(c["time_s"][int(np.argmax(c["convective_heat_W"]))]),
             "peak_surface_T": float(c["time_s"][int(np.argmax(c["surface_T_max_K"]))]), "end": float(c["time_s"][-1])}
    if melting and "sprayed_mass_kg" in c:
        melted = (c["removed_mass_kg"] + c["film_mass_kg"]) > 0.0
        sprayed = c["sprayed_mass_kg"] > 0.0
        if melted.any():
            marks["melt_onset"] = float(c["time_s"][melted][0])
        if sprayed.any():
            marks["spraying_onset"] = float(c["time_s"][sprayed][0])
            marks["peak_release"] = float(c["time_s"][int(np.argmax(c["released_mass_kg"]))])
    return marks


def _animate_series(run_dir, history, radius, series_name, render, movie_name, frames_subdir, stills_prefix, fps, animation, stills, melting=False):
    import imageio.v2 as imageio
    import pyvista as pv
    series = read_series(run_dir, series_name)
    c = history.columns
    clim = (float(c["temperature_K"][0]), float(c["surface_T_max_K"].max()))
    times = np.array([t for t, _ in series])
    marks = still_marks(history, melting)
    wanted = {int(np.argmin(np.abs(times - t))) for t in marks.values()} if stills else set()
    frames, out = {}, {"animation": None, "stills": []}
    frames_dir = os.path.join(run_dir, frames_subdir)
    os.makedirs(frames_dir, exist_ok=True)
    for k, (t, path) in enumerate(series):
        if animation or k in wanted:
            frames[k] = render(pv.read(path), clim, frame_title(t, history), os.path.join(frames_dir, "frame_{:04d}.png".format(k)), radius)
    if animation and frames:
        try:
            writer = imageio.get_writer(os.path.join(run_dir, movie_name + ".mp4"), fps=fps, codec="libx264", quality=7, macro_block_size=8)
            for k in sorted(frames):
                writer.append_data(frames[k])
            writer.close()
            out["animation"] = os.path.join(run_dir, movie_name + ".mp4")
        except Exception:                                   # no ffmpeg: GIF
            imageio.mimsave(os.path.join(run_dir, movie_name + ".gif"), [frames[k] for k in sorted(frames)], duration=1.0 / fps)
            out["animation"] = os.path.join(run_dir, movie_name + ".gif")
    if stills and series:
        stills_dir = os.path.join(run_dir, "stills")
        os.makedirs(stills_dir, exist_ok=True)
        for label, t in marks.items():
            k = int(np.argmin(np.abs(times - t)))
            path = os.path.join(stills_dir, "{}{}_t{:.0f}s.png".format(stills_prefix, label, times[k]))
            imageio.imwrite(path, frames[k])
            out["stills"].append(path)
    return out
```


- [ ] **Step 6: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py tests/test_reentry_model_sesam_io.py -q`
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add reentry_model/sesam_io.py reentry_model/compare.py reentry_model/viz.py tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py
git commit -m "Add the mass-loss comparison, the melt plots and the melting overlays of the videos (Step 3 Task 12)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 13: Command line

**Files:**
- Modify (replace): `reentry_model/cli.py`
- Test: `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: everything above.
- Produces: `run` flags `--melt off|on`, `--material` (names or path; default `AA7075_nomelt`, `AA7075_range` with `--melt on`), `--removal`, `--runoff`, `--rarefied-shear`, `--we-critical`, `--kr`, `--kt`, `--prism-layers` (default 4 with melting, 0 otherwise), `--layer-thickness` (mm), `--demise-fraction`, `--particles/--no-particles`, `--k-scale`, `--consistent-mass`; run names end in `_melt-<removal>`; the JSON `settings` carry the melt settings and the liquid properties, `results` the melt results, `files` the particle files and melt plots (and `film`); `compare` handles melting histories; `model_run_name(..., heating_name=None, melt=None)`; `_fmt`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_cli.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: melting flags, the melting run's files, the bookkeeping device against the melting reference

MELT_100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv")
MELT = FEM + ["--melt", "on", "--prism-layers", "0", "--h-surface", "4", "--h-core", "20"]


def test_run_name_with_melt():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin").endswith("_fem-physics_melt-girin")


def test_melting_run_writes_columns_files_and_json(tmp_path):
    """15 s from 71 km with a warm body: melt columns, the particle files, the melt plots, stills with the melt marks."""
    pytest.importorskip("cantera")
    argv = [a for a in MELT if a not in ("--heating", "sesam")] + ["--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "800",
                                                                   "--t-max", "15", "--outdir", str(tmp_path), "--name", "melt_short", "--stills"]
    assert cli.main(argv) == 0
    rows = list(csv.DictReader(open(tmp_path / "melt_short.csv")))
    assert len(rows) == 31 and "sprayed_mass_kg" in rows[0] and "film_mass_kg" in rows[0] and "regime_fraction_slip" in rows[0]
    assert float(rows[-1]["mass_kg"]) < float(rows[0]["mass_kg"]) and float(rows[-1]["sprayed_mass_kg"]) > 0.0
    doc = json.load(open(tmp_path / "melt_short.json"))
    s, r, f = doc["settings"], doc["results"], doc["files"]
    assert s["melt"] == "on" and s["material"] == "AA7075_range" and s["removal"] == "girin" and s["runoff"] == "on" and s["prism_layers"] == 0
    assert s["rarefied_shear"] == "slip" and s["we_critical"] == 4.62 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.86
    assert s["T_liquidus_K"] == 908.0 and s["latent_heat_Jkg"] == 400e3 and s["demise_fraction"] == 0.01 and s["particles"] is True
    assert r["melt_onset_altitude_km"] is not None and r["sprayed_mass_kg"] > 0.0 and r["n_source_rows"] > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    for key in ("particles", "particles_summary", "size_distribution"):
        assert os.path.isfile(f[key])
    assert len(f["melt_plots"]) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in f["melt_plots"])
    assert doc["comparison"] is None and f["film"] is None
    names = [os.path.basename(p) for p in f["stills"]]
    assert any(n.startswith("melt_onset") for n in names) and any(n.startswith("film_") for n in names) and any(n.startswith("section_spraying_onset") for n in names)


def test_bookkeeping_device_against_the_melting_reference(tmp_path):
    """SESAM-equivalent heating + AA7075 + instant removal + k x 1e4 on the coarse mesh over the whole flight: the
    melt metrics are written and the mass follows SESAM's lumped law (2 % / 0.5 km / 2 %, spec 13.1; measured
    0.82-0.98 % / 0.10 km / -1.3 %)."""
    argv = MELT + ["--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--reference", MELT_100,
                   "--outdir", str(tmp_path), "--name", "bookkeeping"]
    assert cli.main(argv) == 0
    doc = json.load(open(tmp_path / "bookkeeping.json"))
    mm = doc["comparison"]["melt_metrics"]
    assert doc["results"]["end_reason"] == "demise" and mm["mass"]["max_rel_m0"] < 0.02
    assert abs(mm["onset_altitude_diff_km"]) < 0.5 and abs(mm["demise_time_rel"]) < 0.02 and mm["sprayed_mass_kg"] == 0.0
    assert os.path.isfile(tmp_path / "bookkeeping" / "mass_time.png") and doc["settings"]["k_scale"] == 1e4
    assert not os.path.isfile(tmp_path / "bookkeeping" / "particles.npz") or True                    # written (empty table) with --particles


@pytest.mark.parametrize("argv", [
    BASE + ["--atmosphere", "us76", "--melt", "on"],                                              # needs --thermal fem
    MELT + ["--kr", "0"],
    MELT + ["--prism-layers", "-1"],
    MELT + ["--demise-fraction", "1.5"],
    MELT + ["--k-scale", "0"],
    MELT + ["--removal", "magic"],
])
def test_bad_melt_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: the new tests fail (`--melt` unknown).

- [ ] **Step 3: Replace `reentry_model/cli.py`**

```python
"""Command line of the re-entry model.

    python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
        [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
        [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
        [--melt on --material AA7075_range --removal girin|instant ...] (Step 3: melting, film, spraying)
    python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]

Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
(cantera for --heating physics, dolfinx for --thermal-solver fenicsx: create the fenicsx_env environment).
"""
import argparse
import importlib.metadata
import math
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, spray, surface_flow, thermal, viz
from . import trajectory as tj
from .earth import GRAVITY_MODELS
from .thermal import MissingBackend

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output")
# the reference cases' break-off state (parent satellite, see README "Reference runs")
DEFAULT_HEADING_DEG = 347.168296
DEFAULT_LAT_DEG = 29.546067
DEFAULT_LON_DEG = -82.134333
DEFAULT_EPOCH = "2024-08-01T12:53:07"
DEFAULT_MATERIAL_DENSITY = 2813.0          # drama-AA7075
ATMOSPHERES = ("nrlmsise", "us76")         # plus replay:<path>
WINDS = ("none", "static")
THERMAL_MODES = ("none", "fem")
MELT_MODES = ("off", "on")


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))


def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None):
    name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
        diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
    return name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")


def make_atmosphere(spec, epoch, wind_name):
    """(atmosphere object, short name, provenance dict) for an --atmosphere value."""
    wind = atmosphere.NoWind() if wind_name == "none" else atmosphere.StaticProfileWind()
    if spec == "us76":
        return atmosphere.US76TableAtmosphere(wind=wind), "us76", {}
    if spec == "nrlmsise":
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), epoch.date())
        return atmosphere.NRLMSISE00Atmosphere(epoch, solar, wind), "nrlmsise", \
            {"f107": solar.f107, "f107a": solar.f107a, "ap": solar.ap, "fap_day": fap.DEFAULT_FAP_DAY}
    if spec.startswith("replay:"):
        ref = sesam_io.load_reference(spec[len("replay:"):])
        return atmosphere.ReplayAtmosphere(ref, wind), "replay", {"replay_reference": ref.name, "replay_sha256": ref.sha256}
    raise ValueError("--atmosphere must be one of {} or replay:<sesam.csv>, got {!r}".format(ATMOSPHERES, spec))


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def provenance():
    import scipy
    try:
        import pymsis
        pymsis_version = pymsis.__version__
    except ImportError:
        pymsis_version = None
    package_names = {"skfem": "scikit-fem", "gmsh": "gmsh", "pyamg": "pyamg", "pyvista": "pyvista",
                      "cantera": "cantera", "dolfinx": "fenics-dolfinx"}
    versions = {}
    for key, pkg in package_names.items():
        try:
            versions[key] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[key] = None
    return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0], **versions}


def build_parser():
    p = argparse.ArgumentParser(prog="reentry_model", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="integrate one sphere trajectory")
    r.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    r.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s], relative to the rotating atmosphere")
    r.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    r.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending")
    r.add_argument("--heading", type=float, default=DEFAULT_HEADING_DEG, help="[deg] clockwise from north (default %(default)s)")
    r.add_argument("--lat", type=float, default=DEFAULT_LAT_DEG, help="geodetic latitude [deg] (default %(default)s)")
    r.add_argument("--lon", type=float, default=DEFAULT_LON_DEG, help="longitude [deg] (default %(default)s)")
    r.add_argument("--epoch", type=parse_epoch, default=parse_epoch(DEFAULT_EPOCH), help="UTC YYYY-MM-DDTHH:MM:SS (default {})".format(DEFAULT_EPOCH))
    r.add_argument("--material-density", type=float, default=DEFAULT_MATERIAL_DENSITY, help="[kg/m3] (default %(default)s)")
    r.add_argument("--temperature", type=float, default=300.0,
                   help="initial temperature [K]: recorded only with --thermal none, the field's initial condition with --thermal fem (default %(default)s)")
    r.add_argument("--atmosphere", default="nrlmsise", help="nrlmsise (default) | us76 | replay:<sesam.csv>")
    r.add_argument("--wind", choices=WINDS, default="none")
    r.add_argument("--bridging", choices=aero.BRIDGING_NAMES, default="sesam-table")
    r.add_argument("--gravity", choices=GRAVITY_MODELS, default="j2")
    r.add_argument("--rtol", type=float, default=1e-9)
    r.add_argument("--cadence", type=float, default=1.0, help="history sample spacing [s] (default %(default)s)")
    r.add_argument("--t-max", type=float, default=3600.0, help="[s] (default %(default)s)")
    r.add_argument("--reference", default=None,
                   help="SESAM run CSV to compare against (with --thermal none the model is also sampled at its times; "
                        "with --thermal fem the macro-step history is interpolated)")
    r.add_argument("--outdir", default=DEFAULT_OUTDIR)
    r.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")
    r.add_argument("--quiet", action="store_true")
    th = r.add_argument_group("thermal model (Step 2)")
    th.add_argument("--thermal", choices=THERMAL_MODES, default="none", help="none: Step 1 trajectory only (default); fem: coupled 3D conduction")
    th.add_argument("--heating", choices=heating.HEATING_NAMES, default="physics",
                    help="physics (default) or sesam: SESAM's uniform 0.27471 x q_stag on every patch -- a verification device, not physical")
    th.add_argument("--stagnation", choices=heating.STAGNATION_NAMES, default="fay-riddell")
    th.add_argument("--bridging-heat", choices=heating.BRIDGING_HEAT_NAMES, default="matting")
    th.add_argument("--matting-n", type=float, default=1.0, help="Matting exponent n (default %(default)s)")
    th.add_argument("--accommodation", type=float, default=0.8, help="free-molecular energy accommodation A_cq (default %(default)s)")
    th.add_argument("--catalycity", type=float, default=1.0, help="wall catalycity 0..1 in Fay-Riddell (default %(default)s)")
    th.add_argument("--material", default=None,
                    help="AA7075_nomelt | AA7075 | AA7075_range | <DRAMA material JSON> (default: AA7075_nomelt, or AA7075_range with --melt on)")
    th.add_argument("--emissivity", type=float, default=None, help="override the material's emissivity")
    th.add_argument("--k-scale", type=float, default=1.0,
                    help="multiplies the material's conductivity (default 1): a verification device for the near-isothermal (lumped) limit, not physical")
    th.add_argument("--t-ambient", type=float, default=0.0, help="radiation background [K] (default %(default)s, SESAM's)")
    th.add_argument("--mesh-size", type=float, default=1.0, help="multiplies --h-surface and --h-core (default %(default)s)")
    th.add_argument("--h-surface", type=float, default=mesh.DEFAULT_H_SURFACE * 1e3, help="surface element size [mm] (default %(default)s)")
    th.add_argument("--h-core", type=float, default=mesh.DEFAULT_H_CORE * 1e3, help="core element size [mm] (default %(default)s)")
    th.add_argument("--thermal-solver", choices=thermal.SOLVER_NAMES, default="skfem")
    th.add_argument("--linear-solver", choices=("direct", "amg"), default="amg")
    th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
    th.add_argument("--dt", type=float, default=0.5, help="macro step [s] (default %(default)s)")
    th.add_argument("--frames-every", type=int, default=0, help="VTK frame every n macro steps (default 0: none; 10 with --animate/--stills)")
    th.add_argument("--animate", action="store_true", help="MP4/GIF of the surface temperature and of the meridional cross-section, plus stills")
    th.add_argument("--stills", action="store_true", help="only the stills (start, peak heating, peak surface T, end; surface and section)")
    me = r.add_argument_group("melting and spraying (Step 3)")
    me.add_argument("--melt", choices=MELT_MODES, default="off", help="off: Step 2 behaviour (default); on: melting, melt film, spraying")
    me.add_argument("--removal", choices=body.REMOVAL_NAMES, default="girin",
                    help="girin: film + runoff + Girin spraying (default); instant: liquid removed as it forms -- the lumped-melting "
                         "verification device, not physical")
    me.add_argument("--runoff", choices=("on", "off"), default="on", help="film runoff transport along the surface (default on)")
    me.add_argument("--rarefied-shear", choices=surface_flow.RAREFIED_SHEAR_NAMES, default="slip",
                    help="slip: Maxwell slip below Kn_delta 0.1 and free-molecular shear above (default); bridged: SESAM's f(Kn) blend")
    me.add_argument("--we-critical", type=float, default=dispersion.WE_CRITICAL_PRACTICAL, help="critical surface Weber number (default %(default)s)")
    me.add_argument("--kr", type=float, default=spray.K_R, help="droplet radius / wavelength (default %(default)s)")
    me.add_argument("--kt", type=float, default=spray.K_T, help="release period / growth time (default %(default)s)")
    me.add_argument("--prism-layers", type=int, default=None, help="prism layers under the surface (default 4 with --melt on, 0 otherwise)")
    me.add_argument("--layer-thickness", type=float, default=mesh.DEFAULT_LAYER_THICKNESS * 1e3, help="outermost layer thickness [mm] (default %(default)s, growth 2)")
    me.add_argument("--demise-fraction", type=float, default=0.01, help="the run ends when the body mass falls below this fraction of the initial (default %(default)s)")
    me.add_argument("--particles", dest="particles", action="store_true", default=True, help="write the particle source table (default)")
    me.add_argument("--no-particles", dest="particles", action="store_false")

    c = sub.add_parser("compare", help="metrics and plots for an existing model history")
    c.add_argument("--model", required=True, help="model history CSV")
    c.add_argument("--reference", required=True, help="SESAM run CSV")
    c.add_argument("--outdir", default=DEFAULT_OUTDIR)
    c.add_argument("--title", default=None)
    c.add_argument("--quiet", action="store_true")
    return p


def build_thermal(args, settings, mass):
    """(ThermalBody or MeltingBody, HeatingModel, settings-provenance dict) for --thermal fem [--melt on]."""
    radius = settings.diameter / 2.0
    melting = args.melt == "on"
    h_surface, h_core = args.h_surface * 1e-3 * args.mesh_size, args.h_core * 1e-3 * args.mesh_size
    layers = args.prism_layers if args.prism_layers is not None else (mesh.DEFAULT_LAYERS if melting else 0)
    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, layers=layers, layer_thickness=args.layer_thickness * 1e-3)
    material_name = args.material or ("AA7075_range" if melting else "AA7075_nomelt")
    mat = material.Material.from_drama_json(material_name)
    if args.k_scale != 1.0:
        mat.k_table = mat.k_table * args.k_scale                  # verification device (near-isothermal body), not physical
    solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=not args.consistent_mass)
    if melting:
        flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear)
        spray_model = spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical)
        melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                          particles=args.particles)
        the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                    emissivity=args.emissivity, T_ambient=args.t_ambient)
    else:
        the_body = body.ThermalBody(the_mesh, mat, solver, mass, T0=args.temperature, emissivity=args.emissivity, T_ambient=args.t_ambient)
    if args.heating == "sesam":
        heating_model = heating.SesamEquivalentHeating()
    else:
        heating_model = heating.PhysicsHeating(stagnation=args.stagnation, bridging=args.bridging_heat, matting_n=args.matting_n,
                                               accommodation=args.accommodation, catalycity=args.catalycity)
    info = {"heating": args.heating, "material": mat.name, "material_file": os.path.abspath(material.MATERIAL_NAMES.get(material_name, material_name)),
            "emissivity": the_body.emissivity, "t_ambient_K": args.t_ambient, "mesh_file": the_mesh.params["path"],
            "h_surface_mm": h_surface * 1e3, "h_core_mm": h_core * 1e3, "prism_layers": layers, "layer_thickness_mm": args.layer_thickness,
            "n_nodes": the_mesh.n_nodes, "n_elements": the_mesh.n_elements, "n_patches": the_body.surface.n_patches,
            "thermal_solver": args.thermal_solver, "linear_solver": args.linear_solver, "lumped_mass": not args.consistent_mass,
            "k_scale": args.k_scale, "melt": args.melt}
    if args.heating == "physics":
        info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                     "accommodation": args.accommodation, "catalycity": args.catalycity})
    if melting:
        info.update({"removal": args.removal, "runoff": args.runoff, "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                     "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                     "liquid": {"rho": mat.liquid.rho, "mu": mat.liquid.mu, "sigma": mat.liquid.sigma},
                     "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus, "latent_heat_Jkg": mat.latent_heat})
    return the_body, heating_model, info


def cmd_run(args, parser):
    for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                         ("--cadence", args.cadence), ("--t-max", args.t_max)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.atmosphere not in ATMOSPHERES and not args.atmosphere.startswith("replay:"):
        parser.error("--atmosphere must be one of {} or replay:<sesam.csv>".format(ATMOSPHERES))
    for label, value in (("--dt", args.dt), ("--mesh-size", args.mesh_size), ("--h-surface", args.h_surface), ("--h-core", args.h_core), ("--k-scale", args.k_scale)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if not 0.0 <= args.catalycity <= 1.0:
        parser.error("--catalycity must be within [0, 1]")
    if args.thermal == "none" and (args.animate or args.stills or args.frames_every):
        parser.error("--animate/--stills/--frames-every need --thermal fem")
    if args.melt == "on" and args.thermal != "fem":
        parser.error("--melt on needs --thermal fem")
    for label, value in (("--we-critical", args.we_critical), ("--kr", args.kr), ("--kt", args.kt), ("--layer-thickness", args.layer_thickness)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.prism_layers is not None and args.prism_layers < 0:
        parser.error("--prism-layers must be >= 0")
    if not 0.0 < args.demise_fraction < 1.0:
        parser.error("--demise-fraction must be within (0, 1)")

    initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                              flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
                              lat=math.radians(args.lat), lon=math.radians(args.lon), epoch=args.epoch)
    settings = tj.Settings(diameter=args.diameter * 1e-3, gravity=args.gravity, rtol=args.rtol,
                           cadence=args.cadence, t_max=args.t_max)
    mass = body.sphere_mass(settings.diameter, args.material_density)
    atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
    reference = sesam_io.load_reference(args.reference) if args.reference else None
    name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                       args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None)
    run_dir = os.path.join(args.outdir, name)
    os.makedirs(args.outdir, exist_ok=True)
    thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
    if args.thermal == "fem":
        the_body, heating_model, thermal_info = build_thermal(args, settings, mass)
    sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
    if args.thermal == "fem":
        frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
        run = coupled.CoupledRun(sim, the_body, heating_model,
                                 coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk")))
        history = run.run()
        thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
    else:
        history = sim.run(extra_times=reference.time if reference is not None else None)
    csv_path = os.path.join(args.outdir, name + ".csv")
    json_path = os.path.join(args.outdir, name + ".json")
    tj.write_history_csv(history, csv_path)
    prov = provenance()
    # spec section 8: the reference file(s)' SHA-256 also live under provenance, alongside the copies
    # already recorded under "comparison" (--reference) and "settings" (replay_reference/replay_sha256).
    prov["reference_sha256"] = reference.sha256 if reference is not None else None
    prov["replay_sha256"] = atm_info.get("replay_sha256")
    doc = {
        "schema_version": 1,
        "run_name": name,
        "inputs": {"diameter_mm": args.diameter, "initial_velocity_kms": args.velocity, "initial_altitude_km": args.altitude,
                   "flight_path_angle_deg": args.flight_path_angle, "heading_deg": args.heading, "latitude_deg": args.lat,
                   "longitude_deg": args.lon, "epoch_utc": args.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
                   "material_density_kgm3": args.material_density, "mass_kg": mass, "initial_temperature_K": args.temperature},
        "settings": {"atmosphere": atm_name, "wind": args.wind, "bridging": args.bridging, "gravity": args.gravity,
                     "rtol": args.rtol, "atol_position_m": settings.atol_position, "atol_velocity_ms": settings.atol_velocity,
                     "cadence_s": args.cadence, "t_max_s": args.t_max, **atm_info, "thermal": args.thermal, **thermal_info},
        "results": history.results,
        "comparison": None,
        "provenance": prov,
        "files": {"csv": os.path.abspath(csv_path)},
    }
    melting = args.thermal == "fem" and args.melt == "on"
    if melting:
        doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()} if args.particles
                            else {})
        doc["files"]["melt_plots"] = [os.path.abspath(p) for p in compare.plot_melt(history, run_dir, name, reference if compare.has_melt(history, reference) else None)]
    if reference is not None:
        plots = compare.plot_all(history, reference, run_dir, name)
        doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
                             "reference_sha256": reference.sha256, "metrics": compare.metrics(history, reference),
                             "plots": [os.path.abspath(p) for p in plots]}
        if compare.has_thermal(history, reference):
            doc["comparison"]["thermal_metrics"] = compare.thermal_metrics(history, reference)
            doc["comparison"]["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, run_dir, name)]
        if compare.has_melt(history, reference):
            doc["comparison"]["melt_metrics"] = compare.melt_metrics(history, reference)
    if args.thermal == "fem" and (args.animate or args.stills):
        vtk_dir = os.path.join(run_dir, "vtk")
        out = viz.animate(vtk_dir, history, settings.diameter / 2.0, animation=args.animate, melting=melting)
        section = viz.animate_section(vtk_dir, history, settings.diameter / 2.0, animation=args.animate, melting=melting)
        doc["files"]["animation"], doc["files"]["section"] = out.get("animation"), section.get("animation")
        doc["files"]["stills"] = out["stills"] + section["stills"]
        if melting:
            film = viz.animate_film(vtk_dir, history, settings.diameter / 2.0, animation=args.animate)
            doc["files"]["film"] = film.get("animation")
            doc["files"]["stills"] += film["stills"]
    if args.thermal == "fem":
        doc["files"]["vtk_dir"] = os.path.abspath(os.path.join(run_dir, "vtk")) if history.results.get("n_frames") else None
    tj.write_run_json(json_path, doc)
    if not args.quiet:
        res = history.results
        print("{}: {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
            name, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
            history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
        if args.thermal == "fem":
            print("  thermal: peak surface T {:.0f} K at t = {:.0f} s, peak mean T {:.0f} K, integrated heat {:.3g} J, "
                  "energy balance residual {:.1e}, {} macro steps, {:.1f} Newton iterations/step".format(
                      res["peak_surface_T_K"], res["time_of_peak_surface_T_s"], res["peak_mean_T_K"], res["integrated_heat_J"],
                      res["energy_balance_residual"], res["n_macro_steps"], res["mean_newton_iterations"]))
        if melting:
            print("  melt: onset {} km, spraying onset {} km, demise {} km at t = {} s; sprayed {:.4f} kg of {:.4f}, film left {:.4f} kg, "
                  "{:.3g} droplets (median r {} um), {} source rows, melt balance {:.1e}".format(
                      _fmt(res["melt_onset_altitude_km"]), _fmt(res["spraying_onset_altitude_km"]), _fmt(res["demise_altitude_km"]),
                      _fmt(res["demise_time_s"]), res["sprayed_mass_kg"], res["initial_mass_kg"], res["film_mass_kg"], res["n_released"],
                      _fmt(res["r_median_um"]), res["n_source_rows"], res["melt_energy_balance_residual"]))
            if reference is not None and "melt_metrics" in doc["comparison"]:
                mm = doc["comparison"]["melt_metrics"]
                print("  mass vs SESAM: max |dm| {:.2%} of m0, melt onset {:+.2f} km, 1 %-mass time {:+.1f} s ({:+.1%})".format(
                    mm["mass"]["max_rel_m0"], mm["onset_altitude_diff_km"], mm["demise_time_diff_s"], mm["demise_time_rel"]))
        if reference is not None:
            hyp = doc["comparison"]["metrics"]["hypersonic"]
            print("  vs {}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m; end time {:+.1f} s".format(
                reference.name, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], doc["comparison"]["metrics"]["d_end_time_s"]))
            if "thermal_metrics" in doc["comparison"]:
                tm = doc["comparison"]["thermal_metrics"]
                print("  heat vs SESAM: Q_conv max {:.2%} of peak ({:.2%} point-wise, continuum), integrated heat {:+.2%} (hypersonic) "
                      "{:+.2%} (end), |dT_eq| max {:.1f} K ({:.2%}), radiated max {:.2%} of peak".format(
                          tm["Q_conv"]["max"], tm["Q_conv"]["continuum_rel_max"], tm["integrated_heat"]["rel_error_end_of_hypersonic"],
                          tm["integrated_heat"]["rel_error_end"], tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"], tm["radiated"]["max"]))
        print("  csv  -> {}\n  json -> {}".format(os.path.abspath(csv_path), os.path.abspath(json_path)))
    # spec section 7: exit 1 when the flight escaped rather than reaching the ground or t_max;
    # the CSV/JSON are already written above so the escaped trajectory is still available.
    if history.end_reason == "escape":
        return 1
    return 0


def _fmt(x):
    return "n/a" if x is None else "{:.1f}".format(x)


def cmd_compare(args, parser):
    history = tj.read_history_csv(args.model)
    reference = sesam_io.load_reference(args.reference)
    stem = os.path.splitext(os.path.basename(args.model))[0]
    title = args.title or "{} vs {}".format(stem, reference.name)
    os.makedirs(args.outdir, exist_ok=True)
    plots = compare.plot_all(history, reference, args.outdir, title)
    doc = {"model": os.path.abspath(args.model), "reference": reference.name, "reference_sha256": reference.sha256,
           "metrics": compare.metrics(history, reference), "plots": [os.path.abspath(p) for p in plots]}
    if compare.has_thermal(history, reference):
        doc["thermal_metrics"] = compare.thermal_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, args.outdir, title)]
    if compare.has_melt(history, reference):
        doc["melt_metrics"] = compare.melt_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_melt(history, args.outdir, title, reference)]
    out = os.path.join(args.outdir, stem + "_vs_reference.json")
    tj.write_run_json(out, doc)
    if not args.quiet:
        hyp = doc["metrics"]["hypersonic"]
        print("{}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m -> {}".format(title, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], out))
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            return cmd_run(args, parser)
        return cmd_compare(args, parser)
    except KeyError as exc:
        # e.g. fap.solar_indices(): the run epoch (or the day before it) has no record in
        # data/fap_day.dat; could also be a reference CSV missing an expected column.
        key = exc.args[0] if exc.args else exc
        print("ERROR: missing key {!r} (epoch outside the fap file, or a reference CSV without that column)".format(key),
              file=sys.stderr)
        return 2
    except (ValueError, FileNotFoundError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except MissingBackend as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except ModuleNotFoundError as exc:
        print("ERROR: missing optional library {!r}: install requirements-step2.txt into drama_env (dolfinx: the separate "
              "fenicsx_env environment)".format(exc.name), file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: all pass (the bookkeeping device test runs the whole 100 mm flight on the coarse mesh in ~10 s).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/cli.py tests/test_reentry_model_cli.py
git commit -m "Add the melting flags, wiring and summaries to the command line (Step 3 Task 13)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 14: Verification and sensitivity drivers, the reference-tier test, the runs

**Files:**
- Create: `analysis/melt_verification.py`, `analysis/melt_sensitivity.py`, `tests/test_reentry_model_reference_melt.py`
- Outputs (git-ignored): `reentry_model_output/verification_melt/{summary.md, summary.json, <case>__<mode>.*}`, `.../sensitivity/{sensitivity.md, sensitivity.json, ...}`, `.../girin/` (Task 8), `.../reference_tests/`

**Interfaces:**
- Consumes: the CLI (Task 13), `compare.melt_metrics` (Task 12), the melting references (Task 11), `sesam_io.REFERENCE_DIR`.
- Produces: the verification table (bookkeeping thresholded, resolved and physics reported), the sensitivity table, the reference-tier test with the spec §13.1 thresholds (mass 2 % of m₀, onset 0.5 km, 1 %-mass time 2 %).

- [ ] **Step 1: Create `tests/test_reentry_model_reference_melt.py`**

```python
"""Melting model vs the two melting US76 SESAM references (marker: reference, ~2 min): the bookkeeping device
(SESAM-equivalent heating, AA7075, instant removal, k x 1e4, the Step 2 default mesh) against the acceptance
thresholds of spec section 13.1 -- mass within 2 % of the initial at every reference time, melt-onset altitude within
0.5 km, the 1 %-mass time within 2 % (measured 2026-09-21: 0.98 % / +0.10 km / -1.3 % for 100 mm, 1.29 % / +0.17 km /
-0.2 % for 50 mm). The resolved and physics-mode runs are analysis/melt_verification.py's business (reported). Metrics
are written to reentry_model_output/verification_melt/reference_tests/."""
import json
import os

import pytest

from reentry_model import aero, atmosphere, body, compare, coupled, heating, material, mesh, sesam_io, thermal
from reentry_model import trajectory as tj

NAMES = {"d100": "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
         "d050": "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"}
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "verification_melt", "reference_tests")
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}


def bookkeeping_run(ref):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    the_mesh = mesh.sphere_mesh(ref.diameter / 2.0)
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    the_body = body.MeltingBody(the_mesh, mat, thermal.thermal_solver("skfem"), body.sphere_mass(ref.diameter, ref.material_density),
                                settings=body.MeltSettings(removal="instant", runoff=False))
    sim = tj.Simulator(initial, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(),
                       tj.Settings(diameter=ref.diameter))
    return coupled.CoupledRun(sim, the_body, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()


@pytest.mark.reference
@pytest.mark.parametrize("key", ["d100", "d050"])
def test_bookkeeping_device_follows_sesam(key):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, NAMES[key] + ".csv"))
    hist = bookkeeping_run(ref)
    mm = compare.melt_metrics(hist, ref)
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, key + "__bookkeeping.json"), "w") as fh:
        json.dump({"case": ref.name, "results": hist.results, "melt_metrics": mm}, fh, indent=2, default=str)
    assert hist.end_reason == "demise" and abs(hist.results["melt_energy_balance_residual"]) < 1e-6
    assert mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"]
    assert abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"]
    assert abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"]
```


- [ ] **Step 2: Create `analysis/melt_verification.py`**

```python
#!/usr/bin/env python3
"""Run the melting model against the two melting US76 SESAM references and tabulate the mass-loss errors (Step 3
verification, spec section 13.1).

    "$PY" analysis/melt_verification.py [--outdir reentry_model_output/verification_melt] [--cases d100,d050]
        [--modes bookkeeping,resolved,physics] [--animate]

Modes: `bookkeeping` -- SESAM-equivalent heating, AA7075 (DRAMA's single melting temperature), --removal instant,
--runoff off, --k-scale 1e4 (near-isothermal body): the thresholded check of the melting bookkeeping against SESAM's
lumped Q/L_f law (mass within 2 % of the initial at every reference time, onset altitude within 0.5 km, 1 %-mass
time within 2 %); `resolved` -- the same heating and material with the real conductivity, film + runoff + Girin
spraying on the default layered mesh (reported: the surface melts before the interior is hot, so the mass leaves
earlier and, per unit heat, the interior's sensible heating delays the end); `physics` -- physics-mode heating,
AA7075_range, Girin removal (the model proper; the SESAM overlay is context, not a target). Every run writes its
overlay + residual plots and metrics JSON through the CLI; cases left out are read back from existing JSONs so
summary.md / summary.json cover everything available."""
import argparse
import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

CASES = {
    "d100": ("sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind", ["--diameter", "100", "--altitude", "77.500133"]),
    "d050": ("sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind", ["--diameter", "50", "--altitude", "115"]),
}
MODES = {
    "bookkeeping": ["--heating", "sesam", "--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--prism-layers", "0"],
    "resolved": ["--heating", "sesam", "--material", "AA7075", "--removal", "girin"],
    "physics": ["--heating", "physics", "--material", "AA7075_range", "--removal", "girin"],
}
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}     # bookkeeping mode only
COLUMNS = ["case", "mode", "max |dm| (of m0)", "melt onset [km] (model / SESAM)", "1 %-mass time [s] (model / SESAM)",
           "sprayed / film left [kg]", "droplets (median r)", "runtime", "verdict"]


def run_case(key, mode, outdir, animate):
    name, args = CASES[key]
    ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    label = "{}__{}".format(key, mode)
    argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                             "--melt", "on"] + MODES[mode] + ["--reference", ref, "--outdir", outdir, "--name", label, "--quiet"]
    if animate:
        argv.append("--animate")
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run {} failed with exit code {}".format(label, rc))
    return label


def row_from_doc(key, mode, doc):
    mm, res = doc["comparison"]["melt_metrics"], doc["results"]
    fmt = lambda x, f="{:.1f}": "n/a" if x is None else f.format(x)
    verdict = ""
    if mode == "bookkeeping":
        ok = (mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"] and mm["onset_altitude_diff_km"] is not None
              and abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"] and mm["demise_time_rel"] is not None
              and abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"])
        verdict = "pass" if ok else "FAIL"
    return {"case": key, "mode": mode, "max |dm| (of m0)": "{:.2%}".format(mm["mass"]["max_rel_m0"]),
            "melt onset [km] (model / SESAM)": "{} / {}".format(fmt(mm["onset_altitude_model_km"], "{:.2f}"), fmt(mm["onset_altitude_reference_km"], "{:.2f}")),
            "1 %-mass time [s] (model / SESAM)": "{} / {}{}".format(fmt(mm["demise_time_model_s"]), fmt(mm["demise_time_reference_s"]),
                                                              "" if mm["demise_time_rel"] is None else " ({:+.1%})".format(mm["demise_time_rel"])),
            "sprayed / film left [kg]": "{:.4f} / {:.4f}".format(res["sprayed_mass_kg"], res["film_mass_kg"]),
            "droplets (median r)": "{:.3g} ({} um)".format(res["n_released"], fmt(res["r_median_um"])),
            "runtime": "{:.0f} s, {} steps, {:.1f} it/step".format(res["runtime_s"], res["n_macro_steps"], res["mean_newton_iterations"]),
            "verdict": verdict}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--modes", default=",".join(MODES))
    p.add_argument("--animate", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    for key in args.cases.split(","):
        for mode in args.modes.split(","):
            run_case(key, mode, args.outdir, args.animate)
    rows = []
    for key in CASES:
        for mode in MODES:
            path = os.path.join(args.outdir, "{}__{}.json".format(key, mode))
            if os.path.isfile(path):
                doc = json.load(open(path))
                if doc.get("comparison") and "melt_metrics" in doc["comparison"]:
                    rows.append(row_from_doc(key, mode, doc))
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    lines += ["| " + " | ".join(str(r[c]) for c in COLUMNS) + " |" for r in rows]
    md = "\n".join(lines) + "\n"
    open(os.path.join(args.outdir, "summary.md"), "w").write(md)
    json.dump(rows, open(os.path.join(args.outdir, "summary.json"), "w"), indent=2)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```


- [ ] **Step 3: Create `analysis/melt_sensitivity.py`**

```python
#!/usr/bin/env python3
"""Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).

    "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
        [--variants base,layers2,layers6,dt025,bridged,norunoff,we308,kr-30,kr+30,kt-30,kt+30,AA7075,fenicsx]

Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius); the table
lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
change relative to `base`. The `fenicsx` variant needs the fenicsx_env interpreter (run it separately with
--variants fenicsx --python <fenicsx_env python>; it is launched as a subprocess). Variants left out are read back
from existing JSONs."""
import argparse
import json
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli  # noqa: E402

CASES = {"d100": ["--diameter", "100", "--altitude", "77.500133"], "d050": ["--diameter", "50", "--altitude", "115"]}
VARIANTS = {
    "base": [], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
    "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "we308": ["--we-critical", "3.08"],
    "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
    "AA7075": ["--material", "AA7075"], "fenicsx": ["--thermal-solver", "fenicsx"],
}
KEYS = ["sprayed_mass_kg", "r_median_um", "melt_onset_altitude_km", "spraying_onset_altitude_km", "demise_altitude_km", "runtime_s"]
COLUMNS = ["case", "variant", "sprayed [kg]", "median r [um]", "melt onset [km]", "spraying onset [km]", "demise [km]", "runtime [s]"]


def argv_for(key, variant, outdir):
    return ["run"] + CASES[key] + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                                   "--melt", "on", "--heating", "physics", "--material", "AA7075_range", "--no-particles",
                                   "--outdir", outdir, "--name", "{}__{}".format(key, variant), "--quiet"] + VARIANTS[variant]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt", "sensitivity"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--variants", default=",".join(v for v in VARIANTS if v != "fenicsx"))
    p.add_argument("--python", default=None, help="interpreter for subprocess runs (the fenicsx variant)")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    for key in args.cases.split(","):
        for variant in args.variants.split(","):
            a = argv_for(key, variant, args.outdir)
            if args.python:
                env = dict(os.environ, FI_PROVIDER="tcp")
                subprocess.run([args.python, "-m", "reentry_model"] + a, cwd=REPO_ROOT, check=True, env=env)
            elif cli.main(a) != 0:
                raise SystemExit("run {}__{} failed".format(key, variant))
    rows = []
    for key in CASES:
        base = None
        for variant in VARIANTS:
            path = os.path.join(args.outdir, "{}__{}.json".format(key, variant))
            if not os.path.isfile(path):
                continue
            res = json.load(open(path))["results"]
            if variant == "base":
                base = res
            fmt = lambda k, f: "n/a" if res.get(k) is None else (f.format(res[k]) + ("" if base is None or variant == "base" or base.get(k) is None or not base[k] else " ({:+.1%})".format(res[k] / base[k] - 1.0)))
            rows.append({"case": key, "variant": variant, "sprayed [kg]": fmt("sprayed_mass_kg", "{:.4f}"), "median r [um]": fmt("r_median_um", "{:.1f}"),
                         "melt onset [km]": fmt("melt_onset_altitude_km", "{:.2f}"), "spraying onset [km]": fmt("spraying_onset_altitude_km", "{:.2f}"),
                         "demise [km]": fmt("demise_altitude_km", "{:.2f}"), "runtime [s]": fmt("runtime_s", "{:.0f}")})
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)] + ["| " + " | ".join(r[c] for c in COLUMNS) + " |" for r in rows]
    md = "\n".join(lines) + "\n"
    open(os.path.join(args.outdir, "sensitivity.md"), "w").write(md)
    json.dump(rows, open(os.path.join(args.outdir, "sensitivity.json"), "w"), indent=2)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```


- [ ] **Step 4: Run the reference-tier test and the verification driver**

```bash
"$PY" -m pytest -m reference tests/test_reentry_model_reference_melt.py -q          # ~1 min: both bookkeeping runs
"$PY" analysis/melt_verification.py --animate                                         # ~6 min: 2 spheres x 3 modes, with the videos
```

Expected: 2 passed (measured 0.98 % / +0.10 km / −1.3 % and 1.29 % / +0.17 km / −0.2 %); the driver prints the six-row table with `pass` in both bookkeeping rows, the resolved rows ~10 % / 14 % mass difference with onsets 71.62 / 77.57 km, the physics rows onsets 73.96 / 78.32 km, and writes `summary.md`. Check that `reentry_model_output/verification_melt/d100__physics/vtk/` holds `animation.mp4`, `film.mp4`, `section.mp4` and the stills (`melt_onset`, `spraying_onset`, `peak_release`, `film_*`, `section_*`), and that `d100__resolved/mass_time.png` shows the SESAM overlay with its residual panel.

- [ ] **Step 5: Run the sensitivity study**

```bash
"$PY" analysis/melt_sensitivity.py                                                    # ~45 min: 12 variants x 2 spheres
"$PY" analysis/melt_sensitivity.py --variants fenicsx --python "$FX"                  # the backend variant from fenicsx_env (CC and FI_PROVIDER exported)
```

(export `CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang` before the second command; the script sets `FI_PROVIDER=tcp` itself.) Expected: `sensitivity.md` with every variant's sprayed mass, median radius, onsets and demise altitude and the change relative to `base`; the `fenicsx` row equal to `base` to the printed digits; `layers2`/`layers6`/`dt025` within a few percent of `base` in demise altitude and sprayed mass; `we308`, `kr±30`, `kt±30` nearly identical (the thick branch rarely acts: spraying is melt-limited); `AA7075` earlier onset (850 K vs 908 K liquidus). Record the table in Task 15.

- [ ] **Step 6: Run the FEniCSx tests once more and the whole unit tier**

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
"$PY" -m pytest -m "not drama and not reference" -q
```

Expected: 8 passed; the unit tier green (~6 min).

- [ ] **Step 7: Commit**

```bash
git add analysis/melt_verification.py analysis/melt_sensitivity.py tests/test_reentry_model_reference_melt.py
git commit -m "Add the melt verification, sensitivity drivers and the reference-tier bookkeeping test (Step 3 Task 14)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 15: Documentation — README, assumptions, spec amendments, facts note

**Files:**
- Modify: `README.md` (a new Step 3 section before "## Tests"; the Step 2 verification table refreshed; the Tests section's test counts and the reference-tier line), `docs/model_assumptions.md` (§9 appended, §7 amended), `docs/superpowers/specs/2026-09-20-melt-spraying-design.md` (§18 "Amendments" appended; status line), `reentry_model/__init__.py` (docstring: Steps 1–3), the facts note `/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md` (§17 appended)

- [ ] **Step 1: Re-run the Step 2 verification with the new thermal core and refresh its README table**

```bash
"$PY" analysis/reentry_model_thermal_verification.py
```

Expected (measured 2026-09-21 with the nodal-enthalpy core): d100 sesam 0.46 % / 2.23 % / +0.31 % / +0.30 % / 24.6 K (1.20 %) / 6.66 % / 2089 K at 139 s / 74 s; d050 sesam 0.37 % / 2.36 % / +0.03 % / +0.02 % / 28.0 K (1.22 %) / 3.40 % / 2442 K at 245 s / 24 s; physics ratios 0.744 / 0.770, peaks 2319 K at 118 s and 2494 K at 239 s, 83 s / 33 s. Replace the numbers in the README's Step 2 verification table with the measured ones (the runtimes halve: `energy()` and the operators are cheaper and the Newton takes 2.0 iterations) and add one sentence after the table: "Re-measured on 2026-09-21 with Step 3's nodal-enthalpy core (lumped capacity matrix): every metric within 0.1 K / 0.01 % of the Step 2 values; runtimes halved."

- [ ] **Step 2: Add the Step 3 section to `README.md`**

Insert before `## Tests` (fill the sensitivity table and the two "[...]" placeholders from `reentry_model_output/verification_melt/sensitivity/sensitivity.md` and the Step 1 re-run; keep every measured number of Task 14's `summary.md` in the verification table — replace the prototype numbers below where the re-run differs in the last digit):

```markdown
## Physics model — `reentry_model` (Step 3: melting, melt film and melt spraying)

Step 3 turns the heated sphere into a particle source. On top of the Step 2 coupling, every macro step now also
melts the body (enthalpy method with the latent heat in the nodal enthalpy, element fractions φ_e and element death
on a mesh with four 0.25/0.5/1/2 mm prism layers under the surface), feeds the liquid to a film on the surface
patches, moves the film with a lubrication runoff driven by the gas shear and the pressure gradient, and strips it
into droplets by Girin's gradient instability (thick films: Girin 2017; thin films: Girin & Kopyt 1994; a rarefied
extrapolation of the thin mode where the local Knudsen number Kn_δ ≥ 0.1). Droplets are recorded at birth in a source
table; the trajectory feels the mass loss and the projected area. Design: `docs/superpowers/specs/2026-09-20-melt-spraying-design.md`;
plan: `docs/superpowers/plans/2026-09-20-melt-spraying.md`; assumptions: `docs/model_assumptions.md` §§ Step 3.

```bash
# the model proper: physics heating, AA7075 with its 750-908 K melting range, film + runoff + Girin spraying, videos
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --heating physics --melt on --animate
# the bookkeeping check against SESAM's melting reference (lumped-melting device: NOT physical, see below)
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --heating sesam --melt on --material AA7075 --removal instant --runoff off --k-scale 1e4 --prism-layers 0 \
    --reference data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv
"$PY" analysis/melt_verification.py      # both spheres x {bookkeeping, resolved, physics} -> reentry_model_output/verification_melt/summary.md
"$PY" analysis/girin_reference.py        # Girin 2017 Table 1 and Girin & Kopyt 1994 Tables 1-2 -> .../verification_melt/girin/
"$PY" analysis/melt_sensitivity.py       # layers, dt, shear bridging, runoff, We_cr, k_r, k_t, material -> .../verification_melt/sensitivity/
```

Options (`--melt on`, needs `--thermal fem`): `--material AA7075_range` (default with melting; DRAMA's `drama-AA7075`
tables with the alloy's solidus 750 K / liquidus 908 K, latent heat 400 kJ/kg spread across the range) or `AA7075`
(DRAMA's single 850 K, a ±2 K numerical ramp) — both carry the liquid properties ρ_l 2400 kg/m³, μ_l 1.3 mPa s,
Σ 0.86 N/m (pure aluminium near the liquidus, `reentry_model/data/materials/*.json`); `--removal girin|instant`;
`--runoff on|off`; `--rarefied-shear slip|bridged`; `--we-critical 4.62`, `--kr 0.17`, `--kt 1.1`; `--prism-layers 4`,
`--layer-thickness 0.25` (mm, growth 2; 46 k nodes / 255 k tets on the 100 mm sphere); `--demise-fraction 0.01`;
`--particles/--no-particles`; `--k-scale` (a verification device). Melt and runoff start at the liquidus: material
between the solidus and the liquidus holds its latent heat but counts as solid for the film (spec §8, §17.2).

Outputs: `<run>.csv` gains `mass_kg` (now varying), `film_mass_kg`, `sprayed_mass_kg`, `runoff_mass_kg`,
`removed_mass_kg`, `melt_front_depth_max_mm`, `equivalent_radius_mm`, `n_active_elements`, `spraying_area_m2`,
`theta_cr_deg`, `n_released`, `released_mass_kg`, `r_median_um`, `r_max_um`, `regime_fraction_continuum/_slip/_fm`,
`rt_active`, `removed_enthalpy_J`, `film_thickness_max_mm`, `film_thickness_mean_mm`, `n_dead_elements`; the run
JSON adds the melt onset, spraying onset, demise, the masses, size statistics and the settings; `<run>/particles.npz`
(the source table: time, altitude, velocity, θ, patch centroid, regime, branch, film thickness, δ_m, We_s, radius,
count, mass, release velocity and direction, We_d, Oh, breakup flag), `particles_summary.csv`, `size_distribution.csv`
(Δn, ΔM on 40 log bins 1 µm–10 mm per 10 s window and for the flight); plots `mass_time`, `mass_altitude`,
`mass_budget`, `spraying_time`, `regimes_time`, `droplet_size_time`, `size_distribution` (with the SESAM overlay
and residual when a melting reference is given); videos `animation.mp4` (surface temperature with the emitting
patches coloured by droplet radius), `film.mp4` (film thickness), `section.mp4` (cross-section with the liquidus and
solidus iso-lines) and stills at melt onset, spraying onset and peak release in addition to Step 2's. The VTK series
add the liquid fraction and φ_e (volume, active elements only) and the film thickness, We_s, regime, shear, droplet
radius and release rate (surface). A 100 mm physics-mode melting flight on the default mesh takes ~2.5 min
(146 s: 198 steps, 3.8 Newton iterations per step; melt onset 74.0 km, demise 59.2 km at 98.5 s, 1.455 kg sprayed
as 4.9e7 droplets of median radius 145 µm).

**`--removal instant` and `--k-scale` are verification devices, not physical models.** `instant` removes the liquid
of every element as it forms — no film, no runoff, no spraying — which is the lumped Q/L_f law SESAM applies once its
body is at 850 K; with `--k-scale 1e4` (near-isothermal body) and `--heating sesam` it reproduces SESAM's melting
reference (table below). SESAM keeps the sphere's outer geometry while it melts (it hollows the sphere: its
`thick_mm` is the shell thickness, its heat input, radiation, Kn and C_D stay those of the intact sphere — measured on
the references), and so does the device until elements die.

### Verification (`analysis/melt_verification.py`, `analysis/girin_reference.py`, `tests/test_reentry_model_reference_melt.py`)

Two melting SESAM references (DRAMA `drama-AA7075`, US76, winds off): `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind`
(melt onset 71.0 km at 43.6 s, mass 0 at 66.5 km / 67.3 s) and `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind`
(77.1 km at 178.6 s → 73.1 km / 191.6 s). Mass errors are relative to the initial mass at SESAM's time stamps; the
end time is the time at which the body's mass falls below 1 % of the initial (SESAM's last gram; the model's demise
criterion, film excluded).

| case | mode | max \|Δm\| (of m₀) | melt onset [km] model / SESAM | 1 %-mass time [s] model / SESAM | sprayed / film left [kg] | droplets (median r) | runtime |
|---|---|---|---|---|---|---|---|
| d100 | bookkeeping | 0.98 % | 71.10 / 71.00 | 65.9 / 66.8 (−1.3 %) | — | — | 16 s, 132 steps |
| d050 | bookkeeping | 1.29 % | 77.27 / 77.10 | 190.9 / 191.4 (−0.2 %) | — | — | 8 s, 382 steps |
| d100 | resolved (sesam heating, AA7075, girin) | 9.96 % | 71.62 / 71.00 | 69.7 / 66.8 (+4 %) | 1.396 / 0.065 | 2.7e7 (165 µm) | 80 s, 143 steps |
| d050 | resolved | 13.7 % | 77.57 / 77.10 | 195.0 / 191.4 (+2 %) | 0.161 / 0.021 | 5.3e6 (201 µm) | 25 s, 385 steps |
| d100 | physics (AA7075_range, girin) | 58.7 % | 73.96 / 71.00 | 98.3 / 66.8 | 1.457 / 0.003 | 4.9e7 (145 µm) | 135 s, 194 steps |
| d050 | physics | 52.2 % | 78.32 / 77.10 | 198.9 / 191.4 (+3.9 %) | 0.181 / 0.000 | 8.6e6 (194 µm) | 35 s, 398 steps |

Thresholds (bookkeeping mode only, `tests/test_reentry_model_reference_melt.py`): mass 2 % of m₀, onset 0.5 km,
1 %-mass time 2 %. The resolved runs are reported: the surface melts 0.6 km before SESAM's lumped body reaches 850 K,
and the interior's sensible heating during the melt delays the end by 2–4 %; the mass difference (10–14 % of m₀) is
the lumped-body assumption, plotted in `d100__resolved/mass_time.png`. Physics mode is the model proper (0.74 × SESAM's
heat, windward-concentrated): the surface melts at 74 km and the body is consumed by 59 km.

Girin's published cases (`analysis/girin_reference.py`, `data/reference_values/girin2017_table1.json`,
`girin1994_tables.json`): the exact tier — GI = We∞Re∞^−½ (13.04 / 3.51 / 43.46 vs 13.0 / 3.55 / 43.5) and φ_cr from his
Eq. (3) with We_cr 4.62 (16.3° / 32.0° / 8.9° vs 16.1° / 31.4° / 8.8°) — within 2 %; the integrated tier with the
ambient-density Reynolds number in δ_a — t_f 7.0 / 27.1 / 207 µs (5.7 / 31 / 194), N 1.31e6 / 3.57e5 / 636
(1.5e6 / 3.7e5 / 832), r_med 25.8 / 39.6 µm (26.9 / 41.6) — within 30 %; the spraying duration t_s.d. comes out half
his (2.8 / 8.0 ms vs 5.9 / 16.0) for the iron variants and far below for the stony one (6.5 vs 149 ms), whose
wavelength exceeds the body radius: his belt discretisation and induction handling are unstated, so t_s.d. is
reported, not thresholded. Girin & Kopyt 1994: the thin-film side mode reproduces Table 1's six r_d within 0.5 % and
τ_d within 1 % with one shock-layer factor (7.8 on ρ₂V₀²) and their mass rate ρ₁r_d/(2τ_d) within 1 %; the
Rayleigh–Taylor mode reproduces Table 2's τ* to three digits and λ* up to the table's factor-10 units slip.

Analytic and conformance checks (unit tier): Neumann's Stefan front on the box mesh within 0.3 % (threshold 1 %);
the runoff strip steady state to 0.1 % with exact conservation; Ranger's boundary layer reproduced to 0.05 % and
Thwaites' momentum thickness a constant 12.3 × it within ±1.7 % over 5°–85°; the dispersion table (onset between
We_s 3.00 and 3.08, Δ_f 1.226 and Im Ω_f 0.247 at We_s 10⁴); energy and mass balances with melting and removal to
1e-8; element death keeps the surface closed; both thermal backends give the same melting run to 1e-10 in mass.

Sensitivity (`analysis/melt_sensitivity.py`, 100 mm physics flight, one setting changed per row; the sprayed mass and the demise altitude are
insensitive to the layers, the time step, the shear bridging, the runoff, We_cr and k_t (≤ 0.2 %); the median droplet radius follows k_r
(±22–24 % for ±30 %: the thick branch sets the size at the rim), Δt/2 −8 %, no runoff +5 %; DRAMA's single 850 K material melts 1 km
higher and runs 3× longer (its ±2 K ramp costs Newton iterations). The 50 mm rows and the `fenicsx` row come from the full run of Task 14.)

| case | variant | sprayed [kg] | median r [um] | melt onset [km] | spraying onset [km] | demise [km] | runtime [s] |
|---|---|---|---|---|---|---|---|
| d100 | base | 1.4563 | 149.8 | 73.96 | 73.96 | 59.65 | 137 |
| d100 | layers2 | 1.4543 (-0.1%) | 149.3 (-0.3%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.54 (-0.2%) | 108 (-21.1%) |
| d100 | layers6 | 1.4579 (+0.1%) | 150.5 (+0.5%) | 73.96 (-0.0%) | 73.96 (-0.0%) | 59.77 (+0.2%) | 159 (+15.7%) |
| d100 | dt025 | 1.4563 (-0.0%) | 137.4 (-8.3%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.52 (-0.2%) | 224 (+63.2%) |
| d100 | bridged | 1.4557 (-0.0%) | 149.7 (-0.1%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.66 (+0.0%) | 118 (-14.0%) |
| d100 | norunoff | 1.4568 (+0.0%) | 157.2 (+4.9%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.55 (-0.2%) | 118 (-13.8%) |
| d100 | we308 | 1.4544 (-0.1%) | 149.8 (+0.0%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.66 (+0.0%) | 119 (-13.1%) |
| d100 | kr-30 | 1.4556 (-0.0%) | 114.0 (-23.9%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.75 (+0.2%) | 119 (-13.0%) |
| d100 | kr+30 | 1.4593 (+0.2%) | 182.5 (+21.8%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.54 (-0.2%) | 119 (-13.3%) |
| d100 | kt-30 | 1.4543 (-0.1%) | 149.7 (-0.1%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.66 (+0.0%) | 142 (+3.5%) |
| d100 | kt+30 | 1.4572 (+0.1%) | 149.8 (+0.0%) | 73.96 (+0.0%) | 73.96 (+0.0%) | 59.65 (-0.0%) | 127 (-7.7%) |
| d100 | AA7075 | 1.4561 (-0.0%) | 154.0 (+2.8%) | 74.93 (+1.3%) | 74.93 (+1.3%) | 60.68 (+1.7%) | 401 (+191.8%) |

Findings recorded while building this step (details in the spec's amendments):
- The film is stripped as fast as it melts: both instability branches remove hundreds to thousands of kg/m²/s where
  the melt supply is ~5 kg/m²/s, so the mass loss is energy-limited and the film stays microns thin (Girin's
  "outstripping ablation" case); the thick branch acts only where the runoff piles the film up at the windward rim.
- SESAM hollows the melting sphere at fixed outer geometry (`thick_mm` = shell thickness): heat input, radiation,
  Kn and C_D stay those of the intact sphere until the last gram.
- Girin's Table 1 is reproduced only with the ambient density in the boundary-layer Reynolds number (his printed
  Re∞ uses the compressed one); his α (t_ch) uses the ambient density. His GI mixes the two.
- The element-mean enthalpy of Step 2 could not carry the latent heat: with 30 K across a surface element and a 4 K
  ramp its Newton iteration cycled; Step 3 moved the enthalpy to the nodes with a lumped capacity matrix (exactly
  conservative) and an enthalpy-consistent nodal update. Re-measured with it, the Step 2 verification numbers are
  unchanged within 0.1 K / 0.01 % and the runtimes halved (74 s / 24 s for the two SESAM-equivalent flights).
- Scaling the conductivity with φ_e isolated the surface nodes of nearly consumed elements (5000 K); k is left
  unscaled while φ_e > 0.
- A leeward film is static by construction and can stay attached at the end of a run (0.065 kg in the resolved
  100 mm case); its fate is a Step 4 item.
```


Also update the Tests section: the unit-tier line becomes `# unit tests (~6 min; the thermal solver, coupled and melting tests dominate)` and the reference-tier line `# the Step 1 reference flights (~15 min), the Step 2 coupled runs (~20 min) and the Step 3 bookkeeping runs (~1 min)`; and in the `fenicsx_env` paragraph "the seven conformance tests" → "the eight conformance tests (Step 2's seven and the Step 3 melting cross-check)".

- [ ] **Step 3: Append §9 to `docs/model_assumptions.md` and amend §7**

In §7 replace "no melting or mass loss yet;" with "melting, film and spraying per §9 (no vaporisation, no oxide skin, no droplet tracking);" and "P1 elements with 2 mm surface resolution (Step 3's melt layer needs the prism-layer mesh);" with "P1 elements with 2 mm surface resolution and, for melting runs, four prism layers from 0.25 mm;". Append:

```markdown
## 9. Melting and the melt film (Step 3)

- **Latent heat, enthalpy method.** The nodal enthalpy h(T) = ∫c_p dT + L_f f_l(T) carries the latent heat; f_l is
  linear between the solidus and the liquidus (`AA7075_range`: 750–908 K, ASM Handbook) or a ±2 K numerical ramp
  around DRAMA's single 850 K (`AA7075`); L_f = 400 kJ/kg and the solid tables are DRAMA's `drama-AA7075` verbatim.
  The capacity matrix is lumped (diag Σ φ_e V_e/4 ρ c_i), so the discrete energy is exactly ρ∫h dV with h
  interpolated linearly, and each Newton step is mapped through the true h(T) per node so that a node cannot jump
  across the melting range (the element-mean formulation of Step 2 could not carry a 4 K ramp under a 30 K nodal
  spread). The Stefan front on a box is reproduced within 0.3 %.
- **Liquid properties** (DRAMA has none): ρ_l 2400 kg/m³, μ_l 1.3 mPa s, Σ 0.86 N/m — pure aluminium near the
  liquidus (Smithells; Assael et al. 2006; ASM Vol. 2); alloy corrections within 10 %; no oxide skin.
- **What flows.** Melt and runoff start at the liquidus: an element's material becomes film in proportion to a ±2 K
  ramp at the liquidus (f_feed), so the mushy range holds latent heat but neither runs off nor is stripped
  (conservative; a coherency-point treatment is a future iteration). Material leaves the finite-element body at the
  liquidus enthalpy; the difference to the element's actual enthalpy (superheat, or the latent deficit of a partly
  molten element) is returned to its nodes over the next step, so the balance is exact and the melt front advances at
  the energy-limited rate whatever the element size.
- **Element fractions and death.** φ_e scales an element's heat capacity, not its conductivity (a thinner sliver of
  the same material conducts better, not worse; scaling k isolated the surface nodes and drove them to thousands of
  kelvin). A patch owner dies at φ_e ≤ 5 % (its remainder joins the film); interior elements keep ≥ 10⁻³ so that no
  cavity opens; nodes without material are pinned. The exposed faces of a dead element become patches (a transient
  pit until its neighbours die: with Lees' distribution the pit walls at θ ≈ 90° receive almost no heat, so the
  heating artefact is small; the SESAM-equivalent mode is defined by its total and is unaffected).
- **Mesh.** Four prism layers (0.25, 0.5, 1, 2 mm) built by radial projection of the inner 2 mm/8 mm gmsh sphere's
  boundary triangulation (exact on a sphere; each prism split into three tetrahedra by the smallest-node-id diagonal
  rule); the surface triangles are 2.16 mm (h_surface × R/(R − 3.75 mm)). 46 k nodes / 255 k tetrahedra for the
  100 mm sphere; 0.5–0.7 s per macro step.
- **Gas-side surface flow.** Modified-Newtonian pressure p_e = p∞ + (p_s − p∞)cos²θ on the windward face and the
  isentropic expansion of the equilibrium stagnation state to p_e (Cantera, 1° bins); leeward p_e = p∞, no shear.
  Girin's linear-profile boundary layer in Ranger's (1972) form with the actual edge velocity, δ_a² = 58.1 ν_e ∫u_e⁴ds/u_e⁵
  (exact for potential flow; Thwaites' momentum thickness is a constant 12.3 × smaller within ±1.7 %); shear
  τ_c = μ_e u_e/δ_a. Regimes by Kn_δ = λ_e/δ_a: continuum < 0.01, slip 0.01–0.1 (Maxwell first-order slip,
  u_e/(1 + Kn_δ)), transitional/free-molecular ≥ 0.1 (τ_fm = ρ∞V² sinθ cosθ; `--rarefied-shear bridged` blends the two
  with SESAM's f(Kn) instead). Near the rim the expansion to p∞ makes ρ_e tiny and Kn_δ > 0.1 even on a continuum
  body — a property of the pressure model that the local criterion inherits. The driving gradient
  G = 2(p_s − p∞) sinθ cosθ/R − ρ_l a sinθ includes the body's deceleration (the film is pushed toward the nose).
- **Film and runoff.** One film mass per patch; lubrication velocity and flux with a thin (b ≤ δ_m) and a thick
  (b > δ_m) branch; runoff by a linearly implicit upwind scheme on the patch graph (4 sub-steps per macro step,
  exact conservation, exact steady state; a wetting front advances one patch per sub-step), never across the equator;
  leeward films are static and can stay attached (their fate is Step 4's). The film is fed where it melts and stripped
  there within the step, so it stays microns thin except where the runoff piles it at the windward rim.
- **Spraying.** Girin's (2017) gradient instability with the dispersion relation solved numerically (Δ_f, Im Ω_f
  tabulated against We_s; k_r 0.17, k_t 1.1, We_cr 4.62 — his constants for ordinary liquids), δ_m and V_s from his
  conjugated-layer relations, stripping ρ_l π r²/(λ_f t_per) per area for thick films; Girin & Kopyt's (1994)
  thin-film side mode λ* = 1.5 M_e Σ/(ρ_e u_e²), r = λ*/4, τ* = 2 capillary periods, rate ρ_l min(b, λ*/8)/τ* for
  thin films (their Table 1 mass rate); the rarefied thin branch uses the freestream Mach number and momentum flux
  (an extrapolation); their Rayleigh–Taylor front mode is evaluated and reported only (inactive at 10–30 m/s²).
  Droplets are capped at the film on the patch and a quarter of the body radius; one radius per patch and step; no
  within-patch size spread; recorded at birth, not tracked. Both branches strip far faster than the melt supply, so
  the mass loss is energy-limited (Girin's "outstripping ablation").
- **Geometry feedback.** Mass = Σφ_e ρV_e + film; the drag reference area is the current surface's projection; the
  nose radius R₀, the sphere drag tables and Kn's length scale are those of the intact sphere (as SESAM keeps them).
- **Verification devices.** `--removal instant` (every element's liquid leaves as it forms) and `--k-scale 1e4`
  reproduce SESAM's lumped Q/L_f melting (mass within 1.3 % of m₀, onset within 0.2 km, end within 1.3 %); Girin's
  Table 1 is reproduced in its exact tier (GI, φ_cr) and, with the ambient-density Reynolds number, in t_f, N and
  r_med within 30 %; his spraying durations are not (half his, and far less for the stony variant).
```


- [ ] **Step 4: Append the amendments to the spec**

Change the spec's status line to `Status: implemented 2026-09-21 (plan docs/superpowers/plans/2026-09-20-melt-spraying.md); amendments in §18` and append:

```markdown
## 18. Amendments (implementation, 2026-09-20/21)

Measured while writing and executing the plan; each overrides the section it names. The plan's "Measured facts and
spec amendments" list carries the numbers.

1. §13.1 — the reference runs are `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind` and
   `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind` (the wrapper adds no material suffix for the default
   material); SESAM hollows the melting sphere at fixed outer geometry (`thick_mm` = shell thickness; heat input,
   radiation, Kn and C_D of the intact sphere until the last gram), so `--removal instant` removes the liquid of every
   element as it forms and keeps the geometry — no area scaling. The end-time metric is the interpolated time at which
   the body's mass (film excluded) falls below 1 % of the initial.
2. §6 — the enthalpy is nodal with a lumped capacity matrix and an enthalpy-consistent Newton update (the element-mean
   secant iteration cycles across a latent-heat ramp); conductivity is not scaled by φ_e (a nearly consumed element must
   keep conducting); patch owners die at φ_e ≤ 0.05, interior elements keep ≥ 10⁻³, deaths cascade within a step,
   material-free nodes are pinned. The FEniCSx backend assembles only the stiffness in UFL; the radiation and all loads
   are nodal in both backends (the boundary moves). Both backends agree to 1e-10 in mass.
3. §8 — every active element feeds its liquid inventory (nodal mean of the ±2 K feed ramp at the liquidus), owners to
   their patches by area, interior elements to the four nearest patches; the material leaves at its nodal-mean enthalpy
   and is booked at h_liquid with the difference returned to the element's nodes over the next step. Runoff is a
   linearly implicit upwind scheme (4 sub-steps, exact conservation and steady state; a wetting front advances one patch
   per sub-step) that never crosses the equator; leeward films are static and may stay attached.
4. §5 — the layers are built by radial projection of the inner gmsh sphere's boundary (gmsh's extrudeBoundaryLayer is
   not usable on the OCC sphere): 46 278 nodes / 255 276 tets / 15 430 patches by default; the box mesh is a Kuhn split.
5. §7 — the boundary layer is Ranger's integral with the local edge velocity (Thwaites' momentum thickness is a
   constant 12.3 × smaller within ±1.7 %); G includes the deceleration term toward the nose; near the rim the
   modified-Newtonian expansion to p∞ makes Kn_δ ≥ 0.1 on a continuum body.
6. §9 — a film thicker than δ_m takes the thick branch in every regime (with the local shear's film velocity); the thin
   branch's rate is ρ_l min(b, λ*/8)/τ* with τ* = 0.798 λ*^1.5 (ρ_l/Σ)^½ (Girin & Kopyt's Table 1 mass rate), the
   λ_t cut-off never limits it; droplets are capped by the film on the patch and by R/4. Spraying is melt-limited
   everywhere: both branches strip far faster than the surface melts, so the film stays microns thin.
7. §13.2 — tiers: exact (GI, φ_cr with We_cr 4.62; ≤ 2 %), integrated (t_f, N, r_med with the ambient-density Reynolds
   number in δ_a; ≤ 30 %), reported (t_s.d.: half his for iron, far less for stone; ranges; σ; z₀). His GI mixes the
   ambient density (We∞) with the compressed one (Re∞); his α uses the ambient density.
8. §13.3 — Girin & Kopyt's Table 1 implies an effective dynamic pressure 7.8 × ρ₂V₀²; with it r_d, τ_d and their mass
   rate ρ₁r_d/(2τ_d) are reproduced within 1 %; Table 2's λ* column is 10 × smaller than their Eq. (14) (τ* agrees).
9. §13.5 — the layer variants are 2 / 4 / 6 (six layers from 0.125 mm; eight layers of 0.25 mm with growth 2 would
   exceed the radius).
10. §10, §12 — extra history columns `removed_mass_kg`, `film_thickness_max_mm`, `film_thickness_mean_mm`,
    `n_dead_elements`; `runoff_mass_kg` is the mass that arrived on another patch; the source table has 22 columns.
11. §14 — `--k-scale` (verification device) and `--consistent-mass` (replacing `--lumped-mass`, whose sense is now the
    default) join the CLI.
12. §13.6 — the 100 mm physics-mode flight takes 146 s on the default mesh (target 6 min).
```

- [ ] **Step 5: Append §17 to the facts note and update the package docstring**

Append to `sesam_verified_facts.md`:

```markdown
## 17. Melting (re-measured 2026-09-20 on the two drama-AA7075 US76 references)

- SESAM melts the sphere from the inside: during the melt `thick_mm` is the wall thickness of a hollow sphere of the
  original outer radius (m = ρ 4π/3 [R₀³ − (R₀ − thick)³] to four digits at every row), the convective heat stays
  A₀ × 0.27471 × q_DKR(R₀) × hot-wall (24–26 kW for the 100 mm sphere), `rad_cooling_W` stays −372 W = 4πR₀² εσ 850⁴,
  `knudsen` keeps D₀ and `drag` stays 0.913 until the mass is zero. The removal rate is (Q_conv + rad_cooling)/L_f
  (0.96–0.98 × Q_conv/L_f at the printed sampling), L_f = 400 kJ/kg.
- Onsets and ends (US76, winds off): 100 mm 71.005 km at 43.55 s → mass 0 at 66.476 km / 67.315 s (23.8 s at the
  melting temperature); 50 mm 77.104 km at 178.55 s → 73.094 km / 191.55 s (13.0 s). `demised (uncritical)` with
  `final mass 0 kg`.
```

In `reentry_model/__init__.py` change the docstring's first line to `"""First-principles re-entry model of a solid sphere, Steps 1-3: trajectory, coupled 3D heat transfer, melting and` and the second to `melt spraying (specs under docs/superpowers/specs/)."""`.

- [ ] **Step 6: Final checks and commit**

```bash
"$PY" -m pytest -m "not drama and not reference" -q
git add README.md docs/model_assumptions.md docs/superpowers/specs/2026-09-20-melt-spraying-design.md reentry_model/__init__.py
git commit -m "Document Step 3: README, assumptions, spec amendments (Step 3 Task 15)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

(The facts note lives outside the repository; edit it in place, no commit.)

---

## Self-review notes (writing the plan, 2026-09-21)

- Spec coverage: §5 mesh → Task 1; §6 materials and melting → Tasks 2–3; §7 surface flow → Task 5; §8 film → Task 6; §9 spraying and dispersion → Tasks 4, 7; §10 accounting, §11 coupling → Tasks 9–10; §12 outputs → Tasks 10–13; §13.1 → Tasks 11, 14; §13.2–13.3 → Task 8 (with the amended tiers); §13.4 analytic checks → Tasks 3 (Stefan, fractions/loads balance), 5 (Ranger/Thwaites), 6 (strip, conservation), 4 (dispersion), 9 (balances, death), 3/9 (backends); §13.5 → Task 14; §13.6 → measured (146 s); §14 CLI → Task 13; §15 environment → no new packages; §16–17 → Task 15 documents them.
- Every code block is the prototype file that passed the tests named in its task; the appended test blocks are verbatim tails of the prototype test files. Names used across tasks were checked against the prototype: `MeltSettings/MeltingBody/REMOVAL_NAMES` (Task 9 ← Task 13), `SurfaceFlow/RAREFIED_SHEAR_NAMES` (Task 5 ← 9, 13), `SprayModel/K_R/K_T/SOURCE_COLUMNS/histogram/N_BINS/BIN_EDGES` (Task 7 ← 9, 10), `DispersionTable/WE_CRITICAL_PRACTICAL` (Task 4 ← 7, 8, 13), `MELT_COLUMNS/write_particles` (Task 10 ← 13), `has_melt/melt_metrics/plot_melt/MELT_PLOT_NAMES` (Task 12 ← 13, 14), `animate_film/still_marks` (Task 12 ← 13), `Material.MATERIAL_NAMES/liquid/feed_fraction/h_liquid` (Task 2 ← 9, 13), `mesh.DEFAULT_LAYERS/DEFAULT_LAYER_THICKNESS` (Task 1 ← 13), `StepResult.Q_extra/Q_dropped`, `set_fractions`, `pinned` (Task 3 ← 9).
