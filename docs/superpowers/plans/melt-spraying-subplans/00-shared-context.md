# Melting, Melt Film and Melt Spraying of the Re-entering Sphere (Step 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `reentry_model` so that the coupled trajectory + 3D conduction model of Step 2 melts the AA7075 sphere (enthalpy method, element fractions and element death on a prism-layered mesh), forms a melt film on the surface patches (lubrication runoff driven by the gas shear and pressure gradient), strips the film into droplets by Girin's gradient instability (thick, thin and rarefied branches), feeds the mass loss and projected area back to the trajectory, records every droplet release in a source table with size distributions, visualises it, and verifies the whole against two new melting SESAM references, Girin's published cases and analytic solutions.

**Architecture:** New modules `dispersion` (Girin's Eq. 1 solved numerically, cached table), `surface_flow` (the three-branch flow-regime gate, modified-Newtonian + Prandtl-Meyer wall pressure, edge state by isentropic expansion, Ranger-form boundary layer, wall Knudsen number and the melt closure it selects, shear, driving gradient), `film` (lubrication branches, linearly implicit upwind runoff on the patch graph), `spray` (instability branches, release bookkeeping, source rows, histograms), `girin_case` (Girin-as-published driver) and `body.MeltingBody`; `mesh` gains prism layers, the active set with a face table and the box mesh; `material` gains the latent heat, the melting ranges, the feed fraction and the liquid properties; the thermal backends move the enthalpy to the nodes (lumped capacity matrix, exact conservation) with a Newton iteration mapped through h(T), element fractions, pinned nodes and nodal loads; `coupled`, `compare`, `viz`, `cli` and three analysis scripts grow the melt columns, writers, metrics, plots, videos and flags. Per 0.5 s macro step: trajectory advance → aero state → heating → conduction step (with the deferred melt loads of the previous step) → melt step: liquid inventory of every element → film feed; surface flow; lubrication + runoff; spraying and release; element death with film hand-over; mass and projected area for the next advance → history row, source rows, VTK frame. Every module is plain arrays over patches or elements with its own tests; the Step 1–2 behaviour is unchanged when `--melt off`.

**Tech Stack:** Python 3.12 in `drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`); numpy 2.5, scipy 1.18 (`brentq`, `cumulative_trapezoid`, sparse `spsolve`), gmsh 4.15.2, scikit-fem 12.0.2, pyamg 5.3.0, Cantera 3.2.0, pyvista 0.49, imageio-ffmpeg, matplotlib — all already installed for Step 2; **no new packages**. FEniCSx (dolfinx 0.11) in the separate `fenicsx_env` for the backend cross-check.

**Spec:** `docs/superpowers/specs/2026-09-20-melt-spraying-design.md` — read it first. The measured deviations from the spec listed below are binding; Task 15 records them in the spec's amendments section.

## Global Constraints

- Interpreter for everything: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` (never the system python). Unit tier: `"$PY" -m pytest -m "not drama and not reference" -q` (Step 2: 332 passed + 1 skipped in ~4 min; this plan adds ~60 tests and ~2 min). Reference tier: `"$PY" -m pytest -m reference -q` (Step 1 + Step 2 + Task 13's two melting runs, ~2 min more).
- FEniCSx tests and runs use `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python` with `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang` in the environment (README, Step 2 machine notes).
- SI units inside the package (m, s, kg, K, W, J, Pa); the CLI/CSV keep the wrapper's units (km, km/s, deg, mm) plus SI columns; angles in radians inside.
- Never modify DRAMA's databases or the wrapper (`sphere_reentry.py`, `sphere_sweep.py`). The two melting references are produced by the wrapper (Task 13) and committed under `data/reference_runs/`.
- Non-physical devices must say so: `--heating sesam` (Step 2), `--removal instant` and `--k-scale` (Step 3) — docstring, CLI help and README each state that they are verification devices, not physical models (user requirement carried from Step 2).
- Melt and runoff start at the liquidus; the mushy range counts as solid for the film; no coherency parameter (spec §8, §17.2, user decision of 2026-09-20). Droplets are recorded at birth, one radius per patch and step, no within-patch size spread (spec §17.1).
- Size feedback (user decision of 2026-09-21, replacing spec §10's fixed R₀): in the model proper the body Knudsen number uses the current equivalent diameter of the remaining mass, the stagnation radius of the heating and surface flow is fitted to the current windward cap, and the continuum drag coefficient is the modified-Newtonian drag of the current windward silhouette relative to a sphere's; the SESAM verification devices (`--heating sesam` + `--removal instant`) keep D₀, R₀ and the sphere table, which is what SESAM does. Fixed attitude throughout; tumbling is a later iteration (its bound is recorded, not implemented).
- Flow-regime gate (user decision of 2026-09-22, replacing spec §7's Kn_δ regimes): stage one asks whether a **distinct bow shock exists**, not whether the flow is continuum — a continuum construction must not be allowed to certify itself. The body-scale gate (Kn_body on SESAM's own mean free path) selects one of three branches; within the shock-layer branch a wall-scale Knudsen number selects the **melt closure**, never whether spraying happens: Girin's dispersion relation contains no gas parameters, so the melt is sheared and can spray in every branch. The gate is a **declared conservatism**, not a physical deduction — Kn_body > 0.01 does not imply Kn_local > 0.01 — and every step records `kn_body`, `kn_local_stag`, `re_shock`, `flow_branch` and the closure area fractions so that the cost of the conservatism is measurable rather than invisible.
- Never write that a merged or free-molecular patch "sees no more than freestream density": there is real compression in the merged regime and a cold diffusely reflecting wall raises the number density even in free-molecular flow. The free-molecular branch simply evaluates the surface loads from freestream conditions directly, with no compression model.
- Git-ignored output root `reentry_model_output/`; meshes cached under `reentry_model_output/meshes/`; verification outputs under `reentry_model_output/verification_melt/`.
- Repo conventions: module docstrings, `argparse`, exit codes 0/1/2, tests under `tests/` named `test_reentry_model_<module>.py`, one commit per task, commit messages in the imperative like the existing history, ending with `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>` (exactly this line).
- Every task's tests are run with the exact command given in the task; a task is done only when the whole unit tier passes (`"$PY" -m pytest -m "not drama and not reference" -q`).
- Code in this plan was executed and its tests passed on 2026-09-20/21 in a throwaway copy of the package (`drama_env`'s package versions as pinned in `requirements-step2.txt`; the FEniCSx backend in `fenicsx_env`). Transcribe it verbatim; where a task says "replace the file", the block is the complete new file; where it says "append", the block goes at the end of the existing file.

## Measured facts and spec amendments (2026-09-20 to 2026-09-23, prototype in a throwaway copy)

These were measured while writing the plan and override the corresponding spec statements; Task 15 writes them into the spec (§18) and the README.

1. **Reference run names.** The wrapper names default-material runs without a material suffix: the melting references are `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind` and `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind` (spec §13.1 wrote `_mAA7075_nowind`). Measured on them: 100 mm melt onset 71.005 km at 43.55 s, mass 0 at 66.476 km / 67.315 s; 50 mm onset 77.104 km at 178.55 s, mass 0 at 73.094 km / 191.55 s (the facts note's "74.9 → 73.0 km" was the NRLMSISE case).
2. **SESAM hollows the melting sphere at fixed outer geometry.** During melting `thick_mm` is the shell thickness of a hollow sphere of outer radius R₀ (m = ρ 4π/3 [R₀³ − (R₀ − thick)³] to four digits), the heat input stays 4πR₀² × q with R₀ in DKR, Kn uses D₀, the radiated power is 4πR₀² εσ 850⁴ = 372 W throughout and C_D stays 0.913. Consequence: the bookkeeping device keeps the intact geometry (uniform φ_e reduction) and needs no area scaling; `--removal instant` removes the liquid inventory of **every** element as it forms (not only surface elements — a surface-only feed stores latent heat in the interior of the isothermal body and lags SESAM by 5–30 s).
3. **Nodal enthalpy, lumped capacity, enthalpy-consistent Newton.** Step 2's element-mean enthalpy with the secant capacity is a regula falsi anchored at T_old whose fixed point repels inside a steep latent-heat ramp (contraction factor −59 for the isothermal body crossing the ±2 K ramp); with the real conductivity a surface element spans 30 K while its mean crosses a 4 K ramp and the iteration cycled (jumps across the ramp, 3-cycles under damping). The thermal core now uses the **nodal** enthalpy h(T_i) with a **lumped** capacity matrix diag(Σ_e φ_e V_e/4 ρ c_i) — the only capacity matrix whose increment equals the increment of ρ∫h dV with h interpolated linearly (the consistent P1-coefficient matrix left a 3e-4 balance error) — a tangent Newton on R(T) = E/dt + KT − F with E = M(ρ c_sec)(T − T_old), and every update mapped through the true h(T) per node (T_new ← T(h(T_k) + c_p,eff(T_k)(T_new − T_k))) so that a node cannot jump across the melting range, with damping as a fallback. Where a melt film rides a node the enthalpy inverted is the node's **own mixture** of material and film (`enthalpy_mixed`/`temperature_from_enthalpy_mixed`, weighted by `film_weight`), because the film is liquid and its latent plateau is already spent: inverting the material's h(T) there threw nodes clean across the ramp and the iteration settled into a period-3 limit cycle between 777 K and 864 K on 22 nodes (measured 2026-09-22; with the mixture it converges in 4–9 iterations). Both backends implement it identically (the FEniCSx backend assembles only the stiffness in UFL and adds everything nodal in numpy/PETSc, so the radiation follows the moving boundary). `lumped_mass=True` is the default; `--consistent-mass` (skfem only) keeps the consistent form for the analytic checks. Measured: Stefan front within 0.3 % (0.5 mm box), 2.0–2.4 Newton iterations per step without melting, 3–5 with; the Step 2 verification numbers re-measured with the new core: d100 sesam Q_conv 0.46 % / 2.23 % point-wise, integrated +0.31 %, ΔT_eq 24.6 K (1.20 %), radiated 6.66 %, 74 s runtime (was 158 s) — Task 15 refreshes the README table with the full re-run.
4. **Conductivity is not scaled by φ_e.** Scaling k with φ_e (spec §6) isolated the surface nodes of nearly consumed elements (capacity and conductance both ×0.05) and drove them to 5000 K; k stays that of the full element while φ_e > 0 (a thinner sliver conducts better, not worse); only the capacity scales. Patch owners die at φ_e ≤ PHI_DEATH = 0.05 (their remainder joins the film; owners at 10⁻³ made surface nodes swing by hundreds of kelvin between iterates), interior elements keep φ_e ≥ 10⁻³ (no cavity can open); deaths cascade within the step until every owner has φ_e > 0.05; nodes without material are pinned at their temperature and loads on them are counted (`StepResult.Q_dropped`).
5. **Feed rule** (gated by depth since 2026-09-22, fact 28(b): only material the gas shear can reach leaves its element). Every active element feeds its liquid inventory φ_e ρ V_e mean_i f_feed(T_i) (f_feed a ±2 K ramp at the liquidus; for the single-temperature material identical to f_l), owners to their patches by area, interior elements to the four nearest patches by area. The rate is a *fraction of what the element still holds*, which is what makes the melt front move at the energy-limited rate whatever the element size (surface-only feed with per-element death could not: the three tets of a prism reach the liquidus together and at peak heating 3.6 layers melt per step) — and what limits it is not the rule but the energy, because the mass leaves at the enthalpy the **molten part** of the element carries (the f-weighted nodal mean, not the element mean) and arrives holding `enthalpy_liquid`. With `--removal instant` it leaves at h_liquid instead and the difference is a load on the element's own nodes; that is the lumped device that reproduces SESAM's Q/L_f law, and it is the same rule. Two ways of getting this wrong were measured on 2026-09-22 and are the reason for fact 25: debit only the destination and the element keeps a melt fraction it no longer has and feeds it again next step (the surface melted three times faster than the heat allowed); book the film at the mixture enthalpy h(T) instead of the liquid one and melting is free on the ramp (a 100 mm sphere turned entirely to film on a quarter of its latent heat).
6. **Mesh: layers by radial projection.** gmsh's `extrudeBoundaryLayer` is a geo-kernel operation that cannot be attached to the OCC sphere without re-parametrising; the layers are built by projecting the inner gmsh sphere's boundary triangulation radially (inner radius R − 3.75 mm, shells at 46.25/48.25/49.25/49.75/50 mm), each prism split into three tetrahedra by the smallest-node-id diagonal rule (conforming, exact areas/volumes, closed surface). Default 4 layers × (0.25, 0.5, 1, 2 mm): **46 278 nodes / 255 276 tets / 15 430 patches** on the 100 mm sphere (spec estimated ≈55 k nodes); surface triangles 2.16 mm. A face table (every face with its ≤ 2 elements) makes deactivation O(dead elements). The box mesh is a structured Kuhn split (no gmsh).
7. **Runoff: linearly implicit upwind, no explicit sub-steps.** Films driven by the free-molecular shear near the rim move at ~10 m/s and cross the hemisphere many times per 0.5 s step: the spec's explicit CFL scheme needed 1e4–1e5 sub-steps per macro step (sliver patches worst). The transport is (I + Δt_s C) m_new = m_old with the edge coefficients c = q ℓ (t̂·n̂)⁺/(A b) at the start of each of 4 sub-steps: unconditionally stable, positive, conservative to round-off, exact steady state (strip test 0.1 %); a wetting front advances one patch per sub-step (documented limit). A Picard iteration on the fully implicit form does not contract (Δt × c ≫ 1). Runoff never crosses into leeward patches (the flow stops at the equator); leeward films are static and can stay attached at the end of a run (0.065 kg in the resolved 100 mm case).
8. **Boundary layer: Ranger's own integral.** δ_a² = 58.08 ν_e ∫u_e⁴ds/u_e⁵ reproduces Ranger's 2.2 R Re_D^−½ Ψ(θ) exactly under potential flow (0.05 % on 1° bins with a 20× refined quadrature); Thwaites' momentum thickness (spec §7) is a constant 12.3 × smaller within ±1.7 % over 5°–85° and is kept as the cross-check. Ψ(0) = √(48/15) = 1.789. Measured at 71 km / 7.24 km/s: u_e 1.3 km/s at 30°, 1.9 km/s at 45°; δ_a 6–9 mm; Kn_δ 0.007–0.03 (continuum/slip) up to 60°, ≥ 0.1 from ~85° (the modified-Newtonian expansion to p∞ makes ρ_e → 3e-6 kg/m³ at the rim); τ_c 12–46 Pa, τ_fm 1600–1900 Pa; G 2e4–6e4 Pa/m with the deceleration term (30 m/s² × ρ_l) comparable to the pressure gradient. The film's driving gradient is G = 2(p_s − p∞) sinθ cosθ/R − ρ_l a sinθ (the deceleration pushes the film toward the nose).
9. **Spraying branches.** A film thicker than δ_m takes the thick (Girin 2017) branch in every regime, with the film velocity the local shear gives it (τ δ_m/μ_l; for the free-molecular shear this is what strips the rim, where the film piled up to 100 mm otherwise); the thin branch (Girin & Kopyt 1994) uses λ* = 1.5 M_e Σ/(ρ_e u_e²) (their 1.5 M d/We_d: the thickness cancels), τ* = 2 capillary periods = 0.798 λ*^1.5 (ρ_l/Σ)^½ (their Eq. 12), **ṁ = ρ_l min(b, λ*/8)/τ*** — their Table 1 mass rate is ρ₁ r_d/(2τ_d) = ρ₁ λ*/(8τ*), reproduced within 1 % (spec §9's ρ_l b/τ* is amended to this); their cut-off λ_t = λ*/3 never limits the mode (dropped). The droplet radius is capped by the film on the patch ((3 m_f/(4πρ_l))^⅓) and by R/4 (long near-critical waves; the rim's expanded edge state gave 25 mm droplets otherwise). Both branches strip hundreds to thousands of kg/m²/s where the melt supplies ~5 kg/m²/s: spraying is melt-limited, the film stays microns thin, r ≈ 45 µm–1.3 mm with a median ≈ 145–200 µm; We_d ≤ 21 with a few hundred breakup-flagged rows per flight.
10. **Girin 2017 Table 1.** GI = We∞Re∞^−½ is reproduced (13.04 / 3.51 / 43.46) only with We∞ on the ambient density ρ∞ = ρ_a/6 and Re∞ on the compressed ρ_a = 1e-4 kg/m³; α = ρ∞/ρ_m reproduces his t_ch. φ_cr from Eq. (3) matches his table with We_cr = 4.62 (16.3° / 32.0° / 8.9°; 3.08 gives 17 % smaller angles). The rest depends on which density enters δ_a: with the **ambient** Reynolds number t_f (7.0 / 27.1 / 207 µs vs 5.7 / 31 / 194), N (1.31e6 / 3.57e5 / 636 vs 1.5e6 / 3.7e5 / 832) and r_med (25.8 / 39.6 µm vs 26.9 / 41.6) are within 30 %; with the shock density droplets come out 2.5× smaller. The spraying duration is half his for the iron variants (2.8 / 8.0 ms vs 5.9 / 16.0) and 6.5 vs 149 ms for the stony one (λ_f 3.5 mm > R₀; belt discretisation and induction handling unstated). Spec §13.2's tiers become: exact (GI, φ_cr ≤ 2 %), integrated (t_f, N, r_med ≤ 30 %, ambient density), reported (t_s.d., ranges, σ, z₀).
11. **Girin & Kopyt 1994.** Table 1's six r_d imply an effective dynamic pressure 7.8 × ρ₂V₀² (their "deceleration ~10× and re-acceleration to M = 2–3"); with that one factor r_d agree within 0.5 % and τ_d within 1 %. Table 2's λ* column is exactly 10 × smaller than their Eq. (14) (a units slip); τ* agrees to three digits. The RT criterion is **not** inactive: see fact 30 — the deceleration reaches 95 m/s², the criterion fires in 137 of 217 steps of the 100 mm flight, and the test as written was the wrong one.
23. **Cost, re-measured after the amendments.** A macro step on the default mesh costs ~2.0 s (0.6 s conduction; the rest the melt step's 91 Cantera isentropic expansions, the runoff solve and the hull), up from ~1.5 s: the film-temperature amendment (facts 25-26) adds the liquid-enthalpy evaluations and the mixed inversion, which cost 30-90 % of wall clock across the cases (bookkeeping 39 s / 20 s, resolved 269 s / 54 s, 50 mm physics 324 s). The spec's "100 mm physics flight in <= 6 min" held for a flight that demises (2-5 min) but not for one whose remnant survives: the 100 mm physics case now runs to the ground in 1191 steps / 2357 s. Restate the target per macro step, or bound the flight with `--t-max` when only the spraying phase is wanted.
24. **The physics-mode 100 mm sphere no longer demises.** With the shape feedback and the amended surface flow it sprays 1.032 kg of 1.472 kg as 3.3e7 droplets from 74.0 km and **0.440 kg (29.9 %) reaches the ground**; the 50 mm sphere still demises (203.5 s). Re-measured on 2026-09-22 with the film temperature: the sprayed mass moved by 0.5 % and the median droplet radius from 99 to 122 µm (+23 %, comparable to the +-21 % the sensitivity study spans for a single setting change, and worth a row of its own in Task 14's re-run), so the headline -- that this sphere reaches the ground -- is unchanged by that amendment. This is a fixed-attitude result three times over -- the flattening face is held in its maximum-drag orientation, the fitted nose radius cuts the stagnation flux by a further ~20 %, and the leeward shell is never heated -- and DRAMA's own tumbling-averaged disc C_D (0.60, *below* the sphere's 0.91) shows the spread is an attitude uncertainty, not a drag-law one. It must be reported as such in the README and the spec, with tumbling as the next iteration's first item. The verification devices are pinned to `--size-feedback initial`, so the bookkeeping thresholds are untouched by any of it (0.98 % / 1.29 %, unchanged).
12. **Verification results** (prototype, default settings unless stated): bookkeeping device (sesam heating, AA7075, instant, k×1e4, `--prism-layers 0`) — 100 mm: mass within 0.98 % of m₀, onset +0.10 km, 1 %-mass time −1.3 %; 50 mm: 1.29 %, +0.17 km, −0.2 % (thresholds 2 % / 0.5 km / 2 %; the devices are pinned to `--size-feedback initial`, so no later amendment can move them); resolved (sesam heating, AA7075, girin, D₀/R₀): onset 71.62 / 77.57 km, mass 9.35 % / 7.10 % of m₀ off SESAM, 1 %-mass time +4.8 % / −0.1 % (both moved toward SESAM by the film-temperature amendment: the 50 mm mass error halved from 13.5 %, because melt can no longer become film without paying its latent heat); physics mode (AA7075_range, size feedback, the amended surface flow and the film temperature): 100 mm onset 73.96 km, **no demise** (1.032 kg sprayed, 0.440 kg to the ground, 76.42 % of m0 off SESAM at its own time stamps), 3.3e7 droplets of median radius 122 µm, 2.4 g re-solidified, 1191 steps in 2357 s; 50 mm onset 78.32 km, demise 203.5 s, 0.173 kg sprayed, 66.70 % of m0 off SESAM, 2.5e6 droplets of median radius 236 µm, 407 steps in 324 s. The two thermal backends give the same melting run to 1e-10 in mass and 0 K in temperature (spec 0.1 %). The 1 %-mass time is interpolated and, for the model, taken on the body's material (film excluded). Sensitivity (100 mm, spec §13.5 with layers 2/4/6 — eight layers of 0.25 mm at growth 2 exceed the radius; measured before the film-temperature amendment of facts 25–26, which moved the 50 mm flight by under 1 %, so Task 14's re-run is expected to confirm the ranking rather than change it): sprayed mass within 0.2 % for every variant; demise altitude +8.9 % without the size feedback, +2.8 % for Δt/2, within 0.3 % otherwise; median radius −21 % for k_r −30 % (+1 % for +30 %: the film cap), −21 % for Δt/2, +17 % without runoff, −28 % for the single-temperature `AA7075`, which also melts 1 km higher, takes 753 s (its ±2 K ramp costs Newton iterations) and leaves a 1.8 % leeward remnant.
14. **Size feedback** (measured on the coarse-mesh physics flight): the front erodes fastest — the front-most point recedes from +49 mm to −28 mm while the back stays at −50 mm and the transverse radius at 50 mm until the last 20 % of the mass — so the fitted windward-cap radius grows from 50 mm to 140–200 mm within the first 15 % of mass loss (a flat front) and is capped at 1.67 R_t = 83 mm; the stagnation heating factor (R₀/R_nose)^½ is 0.78–0.82 during most of the melt. A cone from the mass centre selects too few patches on the eroded front (the fit collapsed to 4–16 mm), hence the depth-band cap definition. With the feedback the 100 mm physics flight ends at 114.8 s / 55.7 km (98.5 s without), 230 steps, 1.455 kg sprayed, median 154 µm; the 50 mm flight's 1 %-mass time moves from +3.9 % to +5.4 % of SESAM's. The bookkeeping and resolved runs keep D₀/R₀ (`--size-feedback initial`, the CLI default with `--removal instant`; the verification driver sets it for the resolved mode too).
16. **DRAMA's gamma is constant.** `atmosphereData.xml` and our packaged copy of `StaticEnvironmentData.csv` carry gamma = 1.4000 at every altitude from 0 to 150 km (only the oxygen fraction varies, 0.2317 → 0.1170), so "SARA's own value per altitude" resolves to 1.4 everywhere; the loader reads the column and a test asserts it, so a future varying table would be picked up. The identity Kn = (Ma/Re) sqrt(gamma pi/2) reproduces lambda/L to machine precision; the coefficient is 1.4829 (Maxwell, gamma 1.4) against 1.5105 for the Chapman-Enskog hard sphere — a 2 % difference. The freestream gamma is **not** the Prandtl-Meyer gamma (§20).
17. **SESAM's mean free path, recovered.** From the reference CSVs (lambda = knudsen x D at each row) SESAM's lambda is the hard-sphere value with d = 3.65 A to four digits (ratio 1.0001 over the flight) and within 0.5 % of Maxwell-with-Blottner-viscosity. `aero.mean_free_path` already implements exactly that and is adopted verbatim for Kn_body in the gate and in every SARA comparison.
18. **The three-branch gate and where the reference flights fall.** Thresholds Kn_body < 0.01 (distinct shock: standoff Delta/R 0.08-0.14 against a shock 3-10 mean free paths thick gives Delta > 5 lambda at Kn_D <~ 0.01) and Kn_body >= 10 or Ma <= 1 (free molecular), cross-checked by the shock-layer Reynolds number Re2 = rho_inf V R / mu(T0) (merged below ~100). Measured: **100 mm at 70.0 km sits exactly on both gates** — Kn_body 0.0098, Re2 175 — so it is the marginal case your caveat warned about, and 60 km (Kn_body 0.0026, Re2 599) is the first unambiguous anchor; 55 km gives Re2 1038, not 2000, because V has fallen to 5.7 km/s by then. At the 77.5 km break-off every sweep diameter is merged: Kn_body 0.596 (5 mm), 0.149 (20 mm), 0.0596 (50 mm), 0.0298 (100 mm). On the 100 mm physics melting flight melting begins at 73.9 km in the **merged** branch (Kn_body 0.0172, Re2 102), the gate opens at 69.0 km (Kn_body 0.0090) and the branch stays shock-layer to demise at 58.2 km: **69 % of the melting steps and 84 % of the sprayed mass are Girin-certified**, the first 0.23 kg flagged. Kn_body is not monotone — it falls to 0.0056 at 60 km as the density rises faster than the body shrinks, then rises again as the remnant collapses.
19. **The wall Knudsen number.** Evaluated at the wall state (p_w, T_wall) with a geometric length (the nose radius): the wall gas is ideal at these temperatures, so lambda_w = mu(T_w) sqrt(pi R_s T_w / 2) / p_w exactly and **Kn_local is inversely proportional to the wall pressure** — the rarefied patches are the low-pressure ones near the shoulder, which is the physics the gate is meant to catch. Measured at the 100 mm nose at 71 km: lambda_w/lambda_inf = 1/167 and **Kn_local/Kn_body = 1/84** (the extra factor 2 is D against R), against 1/9 and 1/4.5 had the edge state been used — so the choice of wall over edge changes the conservatism by 20x, and the conservatism is much larger than a factor of 3: even the 5 mm sphere at 77.5 km (Kn_body 0.596, solidly merged) has a nose Kn_local of 0.0066, nominally continuum by two orders of magnitude. Along the 100 mm melting flight Kn_local at the nose stays between 4e-5 and 2e-4.
20. **Modified Newtonian + Prandtl-Meyer.** phi* = 43.38 deg (gamma 1.4) / 40.72 deg (1.15) from p*/p02 = 0.5283 / 0.5744. At the 100 mm sphere at 70.0 km (V 7190 m/s, p02 4165 Pa, q_inf 2141 Pa): p_w(90 deg) = 5.22 Pa Newtonian, **144.3 Pa PM(1.4), 248.8 Pa PM(1.15)**, against 112-219 Pa from measured C_p 0.05-0.1 — a 28-48x correction with a factor ~2 spread from gamma, hence `--gamma-pm` (default 1.15). nu(M) validated against Anderson's Table A.5 (nu(2) = 26.3798 deg, nu(5) = 76.9202 deg) and inverted to 1e-5. **The blend must start at phi*, not straddle it**: the PM slope is singular at the sonic point (nu ~ (M-1)^3/2, so dp/dnu ~ nu^-1/3), and a 5 deg window centred on phi* mixes in the clamped PM value below it and puts back a kink 35x the median curvature; a 10 deg window above phi* leaves 9x and a monotone p_w. Consequence: the "rarefied rim" of the earlier design was an artefact of pure Newtonian (rho_e collapsing to 3e-6 kg/m3 and Kn_delta > 0.1 beyond 85 deg). With PM the rim carries rho_e ~ 1e-4 and Kn_local ~ 2e-3 — continuum — so the *rarefied-rim* pile-up mechanism is gone, and with it the 25 mm droplets it produced. Isolated thick patches are a different matter and are diagnosed in fact 27.
21. **The disc endpoint, and HTG's method.** ATDB_CYLINDER at zero angle of attack is **independent of L/D over six decades** (1.824 at Ma 10 for 7e-5 <= L/D <= 70): face-on, a flat cylinder is a disc to the flow, because hypersonic drag is pressure drag on the frontal area while the side wall is parallel to the flow and the base sits in a near-vacuum wake. Sphere/disc = **0.49897 at every Mach number**, i.e. exactly the Newtonian 1/2 — HTG's continuum database is modified Newtonian and the disc entry is C_p,max(Ma). So integrating modified-Newtonian pressure over our eroded shape is not an approximation beyond what DRAMA already does; it is the same method. `data/atdb_disc.json` is extracted from that file. Free-molecular disc/sphere is 1.03, so only the continuum entry is scaled. Thickness re-enters at angle of attack (edge-on C_D runs 0.025 to 9.8 over the same L/D range) and in the heat factor through the wetted area (0.330 to 0.085).
22. **How C_D moves as the sphere flattens.** For a sphere with a flat front of radius fraction s = r_flat/R the Newtonian integral is C_D = C_p,max (1 + s^4)/2 — 0.92 at s = 0, 1.84 at s = 1 — so the excess over a sphere grows as the **fourth power** of the flattened fraction and nothing happens until the nose is more than half flattened. Measured on the flight (inverting the hull integral): s_eff 0.72 at 90 % mass, 0.97 at 50 %, falling back to 0.81 at 5 % as the rim itself melts and the body becomes a smaller rounded cap. The hull runs ahead of the ring-by-ring profile (s 0.9 there at 50 % mass) because it bridges the central crater and squares off the shoulder — the documented upper bound.
15. **Transonic remnant.** A light remnant (the single-temperature `AA7075` variant leaves 27 g = 1.8 % of m₀ of leeward material that the fixed-attitude, windward-only heating never reaches, so the run continues to the ground) reaches its terminal velocity near Ma 1, where SESAM's factor-2 drag step made DOP853 take 1.3e5 RHS evaluations in one macro step and then stall. The step is now a cubic ramp over Ma 0.98–1.02 (`aero.MACH_SWITCH_LO/HI`); the Step 1 reference tier is unchanged (8 passed) since the intact spheres cross Ma 1 in a fraction of a second. Melting runs may end on the ground with a few percent of leeward remnant; the 1 %-mass time is then n/a.
25. **The film's temperature and its re-solidification.** The film is thermally thin — q b / k_l = 0.22 K across a 10 µm film at 2 MW/m² (22 K even at 1 mm) and b²/α = 0.3 ms against the 0.5 s macro step — so it is given no energy equation. Its **mass** goes to the solver (`set_film_mass`), which carries it on the boundary nodes with the liquid capacity, and it is at the surface temperature by construction: it heats, cools and freezes with the surface, and the droplets leave carrying whatever superheat the surface has (measured on the 50 mm physics flight: the droplets leave at 1.0801 MJ/kg against h_liquid = 1.0553, **+2.4 %**, the film reaching 1013 K against a 908 K liquidus with a mass-weighted mean of 929 K — the earlier design had every droplet leave at h_liquid exactly). Three things had to be true for this to close:
   * **The film holds the liquid enthalpy** h(T) + L_f (1 − f_l), not the mixture h(T) (`Material.enthalpy_liquid`). It is liquid by construction — that is what makes it a film — so it carries its latent heat wherever it sits, and feeding it debits the latent heat the mass has not yet paid. With the mixture enthalpy, mass could be relabelled from solid to film for free wherever the surface sat on the ramp: a 100 mm sphere turned entirely to film on a quarter of its latent heat.
   * **The enthalpy is the patch's nodal mean**, i.e. `mean_i h_liq(T_i)` over the patch's three nodes, matching the nodal sum the solver's mass matrix carries. The mean of the enthalpies and the enthalpy of the mean temperature differ by the whole latent heat across the ramp; using the latter left a **−0.4** residual in the coupled balance.
   * **Every transfer books only the difference it carries, at the destination.** Mass that moves at one temperature books nothing; m kilograms going from h_src to h_dst change the accounted energy by m(h_dst − h_src) and apply m(h_src − h_dst) to the nodes they arrive at. Booking the two halves separately — the solid's m h_src on its element's four nodes, the film's m h_dst on its patch's three — closes the balance just as exactly and wrecks the temperature field: the same mass spread by ¼ on one side and ⅓ on the other leaves ±m h/12 on every surface node, and the body swung to 1618 K and −1202 K in two macro steps.
   Re-solidification is the mirror of the feed: the fraction 1 − f_feed(T_patch) of each patch's film returns to its owner element, and the two directions are **netted per element** (they are one equilibrium seen from opposite sides; run separately they cycled 3 % of the body's mass through the film every step with no net effect, pinned the surface at T_feed and paid Newton iterations for it — netted, the cumulative mass moved is 1.15 × the peak film instead of 3 ×). In flight the mirror rule fires rarely and almost entirely on the leeward side: on the 50 mm physics flight it fires in 12 of 407 steps, all inside a 5 s window before demise, and **99.9 % of the 7.15e-6 kg it returns freezes on leeward patches** (1.0e-8 kg windward). The leeward film runs 10–25 K colder than the windward film (899–914 K against 919–926 K, against a 908 K liquidus and a feed ramp whose foot is 906 K) because a leeward patch gets no convective heat and loses heat by radiation and by conduction into the cold rear, whose mean surface temperature is still 830–860 K at that point. Elsewhere the film sprays away long before it can cool through the ramp, so re-solidification is a leeward and end-of-flight phenomenon — which is where a surviving remnant's melt sits. Two limits are declared, not fixed: φ_e is capped at 1, so film whose owner has no room left — or whose owner has died — stays film for good, because the mesh cannot grow a crust outside itself; `film_frozen_fraction` records exactly how much film the enthalpy calls solid, and it must be read as a mass, not as a fraction of steps: on the 50 mm flight it runs at a median 5 % (90th percentile 17 %) of a film that itself peaks at 5.6 % of the body, i.e. **at most 0.47 % of the initial mass is stranded liquid**, and that stranded film is leeward too (8.68e-4 kg leeward against 4.70e-5 kg windward at their peaks), two orders of magnitude more than the mirror rule manages to freeze — which is the measure of the cap. On the 100 mm flight, which survives to the ground, the fraction sits at 1.0 for most of the 1140 wet steps, because the last gram of film rides a cold body for 400 s with nowhere to go; in mass that is **at most 1.86e-3 kg, 0.126 % of m0**, against the 2.4 g the mirror rule did manage to return. And the droplets' superheat is a flight-dependent few per cent: +2.4 % on h_liquid for the 50 mm case, +0.6 % for the 100 mm one — while a 6 s / 400 s heat-and-cool test with no flow at all (nothing sprays, nothing runs off) ends with all of its remaining 30 g at 641 K and every gram of it stranded. The other limit: a thick crust would conduct, which a lumped nodal capacity does not represent.
26. **Deferred melt loads are bounded.** A deferred load is energy the transferred mass delivered to the nodes it arrived at, and a node whose own mass has since melted away has nothing to heat with it: uncapped, hundreds of drained surface nodes were handed loads worth 1e4–1e5 K of their remaining capacity in a single step (which is how the field reached the excursions in fact 25 before either was fixed). No node is now asked to move more than `LOAD_DT_MAX = 1000 K` in one macro step against its current capacity (material + film, `nodal_capacity`); the remainder waits in `pending_load`, is applied as soon as the node has the capacity, and is dropped into `Q_dropped` if the node dies first. The balance sees both terms, so it stays exact either way, and `unapplied_load_J` reports the queue. Measured under the 100 mm physics loads: capped at 100 K the queue held 17 % of the absorbed heat, at 1000 K it drains to ~2 %, and on the two physics flights the queue never exceeds 618 J against 197 kJ absorbed (50 mm) or 1739 J against 1.37 MJ (100 mm) -- **0.13-0.3 %**, with the temperature field clean throughout; uncapped the books are still exact and the field is wrong.
27. **The single-patch film thickness, diagnosed (2026-09-22).** `film_thickness_max_mm` reached 627 mm on the 100 mm physics flight (38 of 1192 steps above 20 mm) and 2431 mm on the 50 mm one (16 of 408) while the area-mean film stayed at 0.07 mm and 0.00 mm; the same behaviour is in the runs of 2026-09-21, so it is not an effect of the film temperature. Instrumenting every stage of the melt step (the feed, the runoff, spraying, the freeze-back and the death hand-over) on the 50 mm flight settles it:
   * On the **50 mm** flight every sample above 20 mm falls in the last 5.2 % (t = 193–203.5 s of 203.5 s), while the body collapses from 3026 patches to 56 and thousands of elements die per step; the winner sits anywhere up to theta = 86 deg, i.e. in the last windward ring. On the **100 mm** flight (instrumented to t = 108 s) it is not an endgame effect at all: 422 of 1976 samples over t = 27–106 s, the winner at a median **theta of 14.6 deg** — the nose crater the shape feedback digs, where tau and G go as sin(theta) and the film has nothing to drive it out.
   * The winner is always a small *surviving* facet — area 0.10–0.27 x the median on the 50 mm case, 0.10–0.13 x on the 100 mm, i.e. admitted by the 0.1 x median sliver filter — and the melt on it does not fit on it: gathered into a sphere it is 3.9 mm across on a 0.53 mm facet, 6.6 mm on 0.55 mm, 13.9 mm on 0.76 mm (50 mm case), 2.2–6.4 mm on 0.5 mm facets (100 mm case). So `m_f/(rho_l A)` is a volume per area, not a depth, and the lubrication picture the number belongs to has already failed. The mass is real (14 mg–3.4 g) and the mass and energy books stay exact.
   * By stage, the concentrators are the **death hand-over** (67 of the 50 mm case's 226 samples, up to 32 % of the body's whole film onto one 0.6 mm2 facet; 87 of the 100 mm case's 422), the **dead-element feed** to the four nearest patches (66 and 87) and the **runoff** (17 and 58). Only 174 of the 100 mm case's 422 coincide with a death, so the feed and the runoff carry it there.
   * Such a facet cannot drain, by construction: the runoff graph keeps only windward-windward edges (so the last windward ring is a sink — `film.Runoff.__init__`), an outflow needs t_hat . n_ij > 0, and non-finite coefficients are zeroed (`edge_coefficients`: "a degenerate (sliver) patch moves nothing").
   Two things were changed. The hand-over now spreads a vanished patch's film over the `NEAREST_PATCHES` nearest survivors by area, exactly as an interior element's melt is handed out, instead of dumping it on the single nearest: the raw peak ratio falls from 2431 mm to 981 mm and the dominant stage becomes the runoff. That fix is worth more than the diagnostic it was aimed at — film spread over four owners finds room where one owner had none, so in the heat-and-cool test the mass that re-solidifies rose from 0.135 kg to **0.456 kg** and the stranded-film fraction fell from 5 % to **zero**. And `film_thickness_max` now reports only facets where a depth means something (b <= sqrt(A) as well as the sliver filter), with the melt that fails that test surfaced as **`film_blob_fraction`** rather than hidden inside a maximum: the reported maximum falls from 2431 mm to **1.672 mm** on the 50 mm flight and from 626 mm to **4.21 mm** on the 100 mm one, with no step above 20 mm in either, while the blob fraction shows the honest exposure: nonzero in 17 of 408 steps (50 mm, up to 0.999 of the film in the collapse, 4.9 % of m0) and in **110 of 217 steps** (100 mm, median 0.131, up to 0.703, at most 0.21 % of m0) -- so on the 100 mm flight half the melting phase has some melt the film model cannot describe, and that has been invisible until now. At flight level the hand-over change moves nothing much: the 100 mm case's sprayed mass goes 1.0243 -> 1.0227 kg, its median droplet radius 101.9 -> 101.0 µm, its re-solidified mass 1.80 -> 1.86 g. What is *not* fixed, and is declared: the concentration itself, and the fact that `b` still feeds the physics. Tracing the consumers: `film.lubrication`'s thick branch has `q = V d/2 + G b^3/(3 mu)`, **cubic in b**, so a blob facet is handed an enormous runoff flux -- which is self-limiting, because `edge_coefficients` then gives it a rate proportional to b^2/A and the linearly implicit scheme drains it in one sub-step, and it is why the pile-up is transient wherever an outflow edge exists at all; `spray.evaluate` computes `We_s = rho_l v_s^2 b / sigma`, **linear in b**, which inflates the Weber number into Girin's dispersion table and so shifts the fastest mode on that facet, bounded afterwards by the droplet-radius caps (the film mass on the patch, and R/4: r stayed <= 690 µm measured); and `rayleigh_taylor`'s criterion `W b^2 rho_l > 3 sigma` is **quadratic in b**, so `rt_active` flags spuriously -- it is a reported diagnostic only, which is why nothing downstream moves. Clamping the branches to b_eff = min(b, sqrt(A)) would be defensible on exactly the grounds that make the diagnostic honest, but it changes droplet sizes on those facets and is a modelling decision, not a bug fix. A patch graph that merges slivers, or a film model that carries a mass per *element* rather than per facet, is the real answer and belongs with tumbling in the next iteration.
28. **The liquid layer, the branch test and the feed gate (2026-09-22/23).** Asked whether `b` accounts for the elements beneath a patch being molten, the answer was no: `b = m_f/(rho_l A)` is the film account -- the mass the feed rule has mobilised -- and the model's premise was that the film therefore *is* the liquid layer. Measured on the 100 mm physics flight, with depths taken from the **present** wall rather than the original outline (an early attempt measured from the original radius and gave 7-14 mm, meaningless once the nose has receded 77 mm):
   * Girin's conjugate melt-layer depth, how far the gas shear reaches into the liquid: **215-303 um**.
   * The film account: **1.7-3266 um**, thinner than the conjugate depth in 77 % of steps.
   * The **contiguous** molten layer beneath the wall: **350-890 um**, thicker than the conjugate depth in 100 % of steps. Cross-checked independently: the mass-weighted distance melt actually travels from its donating element to the patch it lands on is 2.06 mm on the 50 mm flight, and the contiguous-layer walk gives 2.10 mm on the same flight -- two unrelated measurements agreeing to 2 %.
   So the branch test was comparing the wrong thickness: it put **79-99.8 %** of patches on the thin (Girin & Kopyt 1994) branch where the liquid depth says they are thick (Girin 2017, the outstripping regime). Two changes follow.
   **(a) The contiguous liquid depth decides the branch, and nothing else.** `MeltingBody.molten_depth` marches inward from each patch, hopping to the face-adjacent element furthest along the inward normal, and **stops at the first element that is not fully molten**: molten material separated from the patch by solid is blocked from it and must not be credited to it. Contiguity is not a technicality -- with it enforced the thick branch fires on a median of **26.7 %** of patches (up to 97.9 %), not the 94-100 % an any-nearby-molten measure suggested, because most patches have mushy rather than fully liquid material beneath them and correctly stay thin. What the layer does **not** touch is any Weber number, and the record must be exact here because it was twice written down wrongly while the work was in progress: Girin's surface Weber number in the thick branch is `rho_l V_s^2 delta_m / Sigma`, on the **conjugate depth**, as is everything downstream of it (`lambda_f = 2 pi delta_m / Delta_f`, `t_per = k_t delta_m/(V_s Im Omega_f)`, `r = k_r lambda_f`); the thin branch's Weber number is a reported diagnostic on `b`, which is the whole sheared layer there anyway and which nothing consumes; the droplet Weber number of `spray.source_rows` is on the droplet diameter. The only other quantity moved onto the layer is the Rayleigh-Taylor criterion `W b^2 rho_l > 3 Sigma`, which is reported and never applied.
   **(b) Melt leaves an element only where the shear can reach it** (`MeltSettings.feed_depth`, default `conjugate`, `all` for the old behaviour; `--removal instant` is unaffected, so the SESAM-equivalent device keeps the every-element feed it needs). A wall-owning element always may, being in the sheared layer by definition; a buried element may only if its centroid lies within the conjugate depth of the nearest patch; and where the wall Knudsen gate denies the Girin closure there is no conjugate depth, so only wall-owning elements may feed. The surface flow is now evaluated once per step, at the top of the melt step, because the feed needs delta_m.
   **Measured effect** (flight-integrated medians over all droplets released, not per-step medians -- mixing the two is how a +76 % was briefly mis-reported): 100 mm to 108 s, median droplet radius **130.7 -> 178.2 um (+36 %)**, droplets 32.8 -> 19.8 million, sprayed 1.0227 -> 0.9818 kg (-4.0 %), remaining 0.4491 -> 0.4900 kg, 2.25 s per macro step against 2.0; 50 mm whole flight, median radius **232.0 -> 211.2 um (-9 %)**, droplets 2.56 -> 3.06 million, sprayed 0.17442 -> 0.17925 kg (+2.8 %), demise 203.5 s / 69.32 km -> **207.5 s / 68.08 km**, film left 9.00 -> 2.58 g, runtime 306 -> 172 s (-44 %: far fewer elements pass the gate, and the flow is evaluated once). The two spheres move in opposite directions and that is the gate working: the 100 mm case descends into the branch where Girin's closure is certified, so the branch switch dominates; the 50 mm case never does, so only the feed gate acts and it merely makes the film smaller. The fully molten inventory held *in place* in the mesh rises from about 0.1 % to about **1.0 %** of the body -- the stagnant melt no longer delivered to the surface. Balances stay exact (-1.7e-10, +2.2e-10); tiers 200 unit, 15 reference (10.75 min), 8 FEniCSx.
   **A claim withdrawn.** The shell-retention comparison was read as evidence that the model hollows the body under an intact skin. It is not: the measure that tests hollowness -- how full the *surviving* elements are -- shows them 94-99 % full both before and after, so the body loses whole elements as the surface recedes rather than draining them from inside. The retention crossover at 82 s reflects **where the erosion front is** (by then it has eaten through the thin prism layers on the windward face, so there the coarse core *is* the surface, while the shells survive on the never-heated leeward side), which mixes two places and says nothing about radial drainage. What the gate demonstrably does is stop melt being delivered to the surface from up to two millimetres below it, worth four seconds and 1.2 km of extra life on the 50 mm flight.
29. **Runoff and stripping do not double count the shear (2026-09-23).** Asked how a film can run off and be stripped at once if the same shear drives both, three independent analyses agree that it is not double counting, for a reason worth keeping: **wall shear stress is a flux, not a stock** -- N/m2 is momentum per unit area per unit time, delivered for as long as the body flies, so there is no budget to allocate. The magnitudes: the post-shock gas carries a streamwise momentum flux of order 2e5 Pa of which a 200 Pa wall shear is 0.1 %, and the shear work on the liquid (~8 kW/m2) is ~1/25000 of the gas kinetic energy flux; making droplets is cheaper still (surface creation ~103 W/m2 at 2.4 kg/m2/s and 50 um radii, 1.3 % of the shear work). Physically the shear sets up a mean flow in the sheared sublayer, which carries mass along the wall, and the instability is a perturbation on that flow's free surface, which removes mass from it: the mean and the fluctuation about it, not two claims on one budget. The film momentum balance closes as `tau_gas = tau_wall + h dp/dx + rho h a + mdot_strip V_s + tau_wave`, and the **only** legitimate debit of runoff by stripping is the droplet momentum sink `mdot_strip V_s`, whose size relative to the driving shear is `Pi = v_melt delta_m / nu` = v_melt x 480 s/m for these properties: 5 % at a recession speed of 0.1 mm/s (decoupling safe), 48 % at 1 mm/s (not safe, and the sink must then be solved with the profile). **Recommendation: report Pi every step and flag Pi > 0.2.** What *would* be double counting, as a checklist: applying the sheared layer's velocity to mass that is not in the sheared layer -- the error fact 28(b) removes, and the one that made the question worth asking; computing transport and stripping in two independent passes and adding them without a shared inventory cap (the model applies them in sequence, each capped by the mass present, so it is an operator split with a first-order-in-dt error, not a double count); using one depth limit for both parts of the runoff flux; and charging droplet surface energy to both the shear work and the melting enthalpy. The literature treats the coexistence as routine: annular two-phase flow writes advection, entrainment and deposition as three terms in one film mass balance (Hewitt & Hall-Taylor 1970; entrainment onset, Ishii & Grolmes 1975), as do liquid-film cooling (Gater & L'Ecuyer 1970), melt-layer ablation (Bethe & Adams 1959; Roberts 1959) and -- closest to this work -- Bronshten's *Physics of Meteoric Phenomena* (1983), which partitions meteoroid melt-layer runoff and droplet spray as concurrent channels from one film. **One consequence not implemented, and it is the pressure-response question deferred by the user:** the runoff flux needs *two* depth limits in the same expression -- the shear-driven part over min(delta_m, h), the pressure-gradient and deceleration-driven part over the full contiguous liquid depth h. The thick branch has the first right and uses the film thickness for the second, so the pressure-driven flux is under-integrated. That is the exact form the deferred work should take. (Caveat from the analysis itself: it knows Girin's conjugate construction only from the problem statement given to it and has not read the paper, so whether Girin addresses simultaneous stripping should be checked against the primary source.)
30. **The Rayleigh-Taylor criterion, corrected (2026-09-23), and what chasing it taught about diagnostics** (its zero-film-mass headline is superseded by fact 33: the lateral extent belongs to the molten sheet, not to the flagged subset, and measuring it the second way was circular). Fact 28 moved the criterion onto the liquid layer; asked to check that it had taken, it had -- and it now **fires**, in 137 of 217 steps on the 100 mm flight, from 73.3 down to 56.5 km. Two statements in this plan were therefore wrong. Fact 11's "the RT criterion is inactive at our 10-30 m/s2" is wrong twice: measured from the trajectory, the deceleration runs **6 to 95 m/s2** (median 30), and the criterion is satisfied on patches carrying a median 18 % and up to 81 % of the film. And before fact 28 it fired anyway, in 114 of 217 steps, but off the pile-up artefact of fact 27 -- a true flag from a spurious cause, which is the worst kind.
   **The criterion was also the wrong test.** `W h^2 rho_l > 3 Sigma` rearranges to `h > lambda*/2pi`: it is a *depth* condition, the statement that the pool is deep enough for the fastest mode to see it as deep, not an instability condition -- in an unbounded pool RT is unstable at some wavelength for any depth. The missing condition is lateral, and it is the one that matters. From n^2 = (W k - sigma k^3/rho) tanh(k h): unstable for k < k_c = sqrt(W rho/sigma), i.e. for wavelengths *longer* than lambda_c = 2 pi sqrt(sigma/(W rho)); fastest at k_c/sqrt(3), giving the lambda* the model already reports; and a pool of lateral extent L whose rim pins the interface admits only k >= pi/L. So the pool is **stable** if L < lambda_c/2, unstable but **slower than tau*** if L < lambda*/2, and at full rate only above that. Thresholds: L > **10.9 mm** at 30 m/s2, **6.1 mm** at 95. The tanh(k h) factor the deep-layer formulas omit is not negligible either: 0.52 at h = 2 mm and W = 30 m/s2, so a shallow pool grows at half the deep rate.
   **Measured on the 100 mm flight**: the largest *contiguous* region of deep-enough melt (connected components on the patch graph, as an equivalent diameter) is **4.71 mm across at the median, 7.64 at the 90th percentile, 9.98 at its widest** -- below even the most permissive threshold. So under the bounded test the film mass on genuinely unstable patches is **zero at the median step**, against 35 % under the depth-only test the model had been reporting. Where the pool does widen enough, the bounded growth time is 12.6-83 ms against the shear mode's 0.17-0.21 ms, a ratio near 100.
   **The growth-time race, measured on both spheres**: 50 mm (thin-film mode throughout, the wall Knudsen gate denying Girin's closure) RT 3.1-33.7 ms against shear 0.67-2.31 ms, ratio **2 to 20** (median 11); 100 mm (thick branch on a median 41 % of patches) RT 7.3-39.6 ms against shear **0.17-0.68 ms**, ratio **26 to 137** (median 56). The shear instability wins everywhere, and these ratios are conservative because they use the unbounded RT time.
   So the mode stays **evaluated and reported, never applied**, now for two independent measured reasons: usually no unstable wavelength fits the pool at all, and where one does it is outrun by one to two orders of magnitude. Recorded every step: `rt_mass_fraction`, `rt_region_mm`, `rt_wavelength_over_nose`, `rt_bounded_fraction`, `rt_bounded_growth_ms`, `rt_growth_ms`, `spray_growth_ms`. Comparing lambda* against the *nose* diameter, as the first attempt at this diagnostic did, is wrong twice over -- wrong length, wrong wavelength -- and overstated the mode's availability from zero to 35 % of the film.
31. **A fragile diagnostic, and the reproducibility it appeared to disprove (2026-09-23).** `film_thickness_mean` divided the film mass by the area of **every patch with a nonzero film**, so a patch holding 1e-20 kg contributed its whole area to the denominator. A last-bit difference in the temperature field, enough to nudge one element across the feed ramp, then moved the reported mean by 1-2 % while the film mass was bit-identical to nine figures. Floored at `spray.B_MIN` -- the same depth below which spraying ignores a film -- it is stable to the last bit. The floor needs the same empty-film guard `film_thickness_max` already carries: under `--removal instant` there is no film account at all, and the comparison would otherwise broadcast a (0,) array against the patch areas. The unit tier caught that, no physics run did -- the verification devices exercise paths the physics flights never take, which is a second reason to keep them. The mass-weighted diagnostics (`film_blob_fraction`, `film_frozen_fraction`) never had the problem, because a patch holding nothing contributes nothing to a mass fraction.
   Chasing that 1-2 % wobble cost four wrong explanations in a row -- solver non-determinism, the mesh cache, the epoch, multithreaded arithmetic -- each disproved by measurement before the diagnostic itself turned out to be the whole of it. The lesson for the plan is about which diagnostics to trust: **a ratio whose denominator counts patches by existence rather than by amount will amplify round-off into per-cent scatter.**
   **Reproducibility, measured properly afterwards**: two runs of the 100 mm case agree to **1.4e-9** in the worst column anywhere in an 81-step history and to **1.3e-13** in sprayed mass, 4.3e-13 in median droplet radius; the 50 mm full flight agrees to 1e-13. The model is reproducible and there is **no 1e-4 comparison floor** -- an earlier note to that effect was drawn from the contaminated diagnostic and is withdrawn. One 6e-5 difference between two 100 mm runs remains unexplained by code; the probable cause is a torn source read, because `melt_step` imports `spray` lazily at the first melt step, so editing the package while a run is in flight can mix versions. **Procedural rule: do not edit the package while a measurement is running.**
32. **Girin's critical angle on every branch, the wave-fits test, and the deceleration that drives Rayleigh-Taylor
    (2026-09-24).** Asked why the nose sprays at all when the wall velocity vanishes at the stagnation point, three
    measurements on the 100 mm flight say the flow is right and the gating was not. The edge velocity falls from
    2931 m/s at the shoulder to **46.6 m/s inside two degrees** (a factor of 63), the wall shear from 41.8 to 1.50 Pa
    and the film surface velocity from 5.6 m/s to 0.27 m/s, so everything vanishes toward the stagnation point as it
    must; and Girin's surface Weber number falls to **0.03** there, 150 x below the critical 4.62, crossing it between
    15 and 20 degrees -- his own published phi_cr is 16.1-16.3 deg. The **thick** branch honoured that exactly: zero
    grams released inside 10 deg over the whole flight. Every gram that did leave inside 10 deg -- 2.23 g, 0.23 % of the
    sprayed mass -- came from the **thin** branch, which had no stability threshold at all.
    **Reading Girin & Kopyt (1994) settles why.** Their side-surface analysis is inviscid and unstable at *every*
    wavelength ("unlimiting potential intensity" of high-frequency disturbances), so it carries no critical Weber
    number; the only cut-off they impose is viscous dissipation, their Eq. (9) lambda > 0.5 M d / We, and since their
    dominant wave is lambda* = 1.5 M d / We = 3 x that limit, it can never bind. So the *absence* of a Weber gate on the
    thin branch was faithful to the source. Their wavelength and growth time are faithful too: Eq. (11) is
    lambda* = 1.5 M d / We with We = rho_2 V_0^2 d / Sigma, in which d cancels to give lambda* = 1.5 M Sigma /
    (rho_2 V_0^2) -- exactly what the model computes -- r_d = lambda*/4, and Eq. (12) tau* = 2 tau_v(lambda*) =
    0.798 lambda*^1.5 (rho_l/Sigma)^1/2, matching `CAPILLARY_TAU` = 4 pi/(2 pi)^1.5 = 0.79788 to five figures.
    **What is not faithful is where the mode is applied.** They derive it *for the side surface* (sin Theta -> 0,
    V_0 -> V_inf) and assign the *front* surface, where the influence of the flow is negligible and V_0 -> 0, to a
    different and aperiodic solution: the Rayleigh-Taylor mode of their Eq. (14). The branch selection here has no
    angular dependence -- it turns only on the liquid depth against delta_m -- so the side-surface asymptotics run from
    the stagnation point to the equator, and at 0-2 deg the released wave is **9.36 mm long, 6.1 x wider than the facet
    it sits on**, with the release 100 % supply-limited at every angle (the instability asks for 15 to 460 x the film
    present, so its rate never mattered). Implementing their front/side split is **future work** (fact 33).
    **Three changes made now.** (a) Girin's criterion `We_s > We_cr` gates *every* branch, with
    We_s = rho_l V_s^2 min(delta_m, layer)/Sigma -- his phi_cr evaluated against the model's own local flow rather than
    his closed form, which keeps it valid on an eroded body where one angle does not. (b) `wave_fits`: one whole
    wavelength must fit inside the **contiguous molten region** the patch belongs to, from the new
    `MeltingBody.region_extent`; the region and not the facet, because a facet edge is a bookkeeping boundary and using
    it would make the limit mesh-dependent. (c) The Rayleigh-Taylor driving is the deceleration **normal to the film**,
    their W sin(Theta), which on this body is W cos(phi) -- whole at the stagnation point, zero at the equator where the
    deceleration lies in the surface; the whole deceleration over-drove every patch off the nose by 1/cos(phi).
    **Measured, 100 mm flight, before -> after.** The critical angle now holds: the smallest angle releasing anything
    rises from a median **6.32 to 14.08 deg**, against Girin's 16.1-16.3, and the spraying area halves (61.2 to
    33.5 cm2 median). The unstable share of the wetted area inside 10 deg collapses by about nine times (25.9 to 2.9 %
    at 0-2 deg). But the **mass budget barely moves** -- sprayed 0.98165 to 0.97531 kg, **-0.65 %** -- precisely because
    the release is supply-limited: the film leaves through whichever patches remain unstable. Droplets become fewer and
    larger, count 19.67e6 to 16.47e6 (**-16 %**) and median radius 178.4 to 187.9 um (**+5.4 %**), because the film
    accumulates thicker before it finds an unstable patch. Energy balance stays exact at -8.0e-11, runtime -8 %.
    The near-nose release falls only from 2.23 to 1.98 g: the two gates move *where* the melt leaves, not how much.
    **The 50 mm sphere confirms it**, and it is the harder test because the wall Knudsen gate denies Girin's closure
    there so the flight is thin-branch throughout: demise is unmoved to the step (207.5 s, 68.082 km), sprayed mass
    unmoved to five decimals (0.179251 to 0.179246 kg), droplet count -0.25 %, and the only visible effect is the
    median radius, 211.2 to 224.5 um (**+6.3 %**), with the balance still exact and runtime -28 %. So on both spheres
    the melt *supply* sets the mass loss and the instability rate sets only the droplet size and the release map --
    which is the same supply-limited conclusion from the other side.
33. **The bounded Rayleigh-Taylor criterion, corrected again, and what now decides it (2026-09-24).** Measuring the
    wave against the molten region rather than the flagged subset (fact 32b) **reverses fact 30's headline number** and
    the reason is worth keeping. Fact 30 measured the lateral extent as the largest contiguous region of patches that
    *already satisfied the depth criterion* -- 4.71 mm at the median -- and concluded the bounded criterion held for no
    film mass. That was circular: it used the criterion's own output to define the criterion's domain. The physical
    domain is the connected liquid sheet, whose rim is what pins the interface, and depth enters separately through the
    per-patch tanh(k h) factor and the k < k_c test. Measured on the molten sheet the extent is **125 mm**, essentially
    the whole windward face, and the bounded criterion is then satisfied for a median **69 %** of the film mass instead
    of 0 %.
    **The conclusion is unchanged and much stronger.** A sheet that wide admits long waves, but they are marginal ones:
    the bounded growth time rises from 18.8 ms to **4614 ms** against the shear mode's **0.170 ms**, a ratio of about
    **27000** where fact 30 measured 100. The admitted population is the near-equator patches where W cos(phi) is small,
    k_1 = pi/L sits just below the cutoff k_c, and the growth rate consequently goes to zero. The normal-component fix
    also tightens the depth criterion on its own: it fires in 118 of 217 steps instead of 137, over a median 21.6 % of
    the film mass instead of 34.8 %. So the mode stays **reported, never applied**, now because it is outran by four
    orders of magnitude rather than because no wave fits.
    **One caveat, not implemented.** The equivalent diameter 2 sqrt(A/pi) of a molten sheet wrapped over the windward
    face is 125 mm on a body whose transverse diameter is about 100 mm and falling: a 25 cm wave cannot exist on a 10 cm
    sphere however large the wetted *area*, because the surface curves away. Capping the extent at the body's transverse
    diameter is the obvious correction and is left as a decision, since it changes k_1 only from 25.1 to 31.4 m^-1 and
    would not change the four-orders-of-magnitude conclusion.
34. **Future work, recorded so it is not lost.** (a) **Girin & Kopyt's front/side split for the thin branch**: select
    the side mode where sin Theta = cos(phi) is small and their front-surface Rayleigh-Taylor solution (Eq. 14) where it
    approaches one, instead of applying the side asymptotics everywhere. Their paper gives two limits and no blending,
    so the transition is a modelling choice that must be stated rather than buried. (b) The extent cap of fact 33.
    (c) The pressure-driven runoff term of fact 29, which needs the full contiguous liquid depth where the model uses
    the film. (d) The ALE mesh of fact 28's preamble.
35. **The front-surface Rayleigh-Taylor mode, applied (2026-09-24).** Asked to apply the criterion so it can spray, with
    the expectation that the nose would be the only place it is met. **It is**: the mode releases mass out to 20 deg and
    **exactly zero beyond it**, with 99.1 % of that mass inside 15 deg and 58.3 % inside 10. Switched by
    `--rt-spray on|off`, default on; `off` reproduces every run before this date.
    **The gate needs both criteria, and finding out why was instructive.** The *bounded* criterion has no threshold of
    its own: with a wide molten region k = k_max is always admissible and tanh(k h) > 0 at any depth, so it answers yes
    wherever there is liquid -- depth enters its *rate*, not its yes/no. Gating on it alone made every wet patch
    Rayleigh-Taylor-unstable, which a unit test caught at once. The criterion with a threshold is the depth one,
    `W cos(phi) h^2 rho_l > 3 Sigma`, equivalently h > lambda*/2 pi, and that is exactly what confines the mode to the
    nose: the normal deceleration falls as cos(phi), so the depth required rises as 1/sqrt(cos phi) while the melt gets
    shallower. Both are now required -- the depth criterion for the threshold, the bounded form for the admissible
    wavelength, the lateral fit and the growth rate.
    **Three modelling decisions, none of them Girin & Kopyt's**, because their paper gives no mass-loss rate for the
    front surface (its Table 2 has no mdot column, unlike Table 1 for the side). (a) The rate carries their side-surface
    construction across: for the side they give r_d = lambda*/4 and tau_d = tau*, which this model writes
    mdot = rho_l min(b, r_d/2)/tau*, so for the front, where they give r_d ~ lambda* and tau_d ~ tau*, the same form is
    mdot = rho_l min(b, lambda*/2)/tau*. Since lambda* is centimetre-scale (21 mm at 95 m/s^2) and b millimetre-scale it
    reduces to **rho_l b / tau***: the film is shaken off within one growth time, which is what their text describes, and
    **no new closure constant is introduced**. (b) The droplet radius is lambda*, their front-surface r_d, then capped as
    on every branch by the film mass on the patch and a quarter of the body radius -- and here the film-mass cap is what
    binds, giving 1.2-1.9 mm rather than 21. (c) Two modes of one interface cannot both break it up, so the **shorter
    growth time takes the patch** and the other does not act on the same film; that is also what keeps the shear from
    being spent twice (the double-counting checklist of fact 29).
    **Measured, 100 mm flight.** The mode releases **61.04 g, 6.27 % of the sprayed mass**, and the crossover is sharp
    because the two growth times move in opposite directions: tau_RT is nearly constant at **10-13 ms** across the whole
    cap while tau_shear rises from 0.27 ms at 30-45 deg through 1.41 ms at 10-15 deg to **35.6 ms at 0-2 deg**. They
    cross between 2-4 and 4-6 deg, and beyond 15 deg Girin's thick branch wins outright (58.4 g at 15-20 deg against the
    front mode's 0.55 g). The **total** sprayed mass moves only **-0.16 %** (0.97531 to 0.97375 kg) -- supply-limited
    once again, so what changes is the mechanism and the droplet population, not the mass.
    **The size distribution is where it shows.** The median radius is unmoved, 187.9 to 188.6 um, but the tail is
    transformed: the per-step largest droplet goes from a median **388.6 to 2027.0 um** and the largest anywhere from
    1551 to **4093 um**. So about 6 % of the sprayed mass now leaves as droplets an order of magnitude larger in radius
    and three in volume, which is precisely the inertial "shaking off" Girin & Kopyt describe at the front surface and
    associate with meteoroid flares. **For Step 4 the droplet source is now two populations, not one**: a ~190 um shear
    population from the 20-90 deg annulus and a ~1.7 mm inertial population from the stagnation cap.
    **It is a large-body mechanism here.** On the 50 mm sphere it barely fires: demise unmoved (207.5 s, 68.082 km),
    sprayed mass unmoved to six figures, droplet count to five, median radius -0.10 %, and the only visible effect is the
    largest droplet, 1359 to 2567 um. The threshold h > sqrt(3 Sigma/(W cos(phi) rho_l)) is an **absolute length** --
    3.36 mm at 95 m/s^2 -- so the 100 mm sphere's pool reaches it and the 50 mm sphere's does not. The energy balance
    stays exact on both flights (-1.7e-10 and 2.5e-10 of the absorbed heat).
    **A procedural near-miss worth recording.** The first 50 mm comparison was of the wrong run: the output directory
    already held a run from the previous day, and a wait loop testing only for the file's *existence* found it
    immediately and returned stale numbers that happened to reproduce an older baseline exactly -- which is what gave it
    away. Fact 31's rule extends: **check that an output path is fresh, not merely that the file is there**, and prefer a
    directory that did not exist before the run.
36. **Regime 2: Girin's Kelvin-Helmholtz cell, implemented (2026-09-25).** `spray.BRANCH_REGIME2 = 4` adds the one cell
    of Girin's own regime classification the model did not carry, and his two-stage test is now applied in his order.
    **Stage one, on the depth.** delta_m is a *theoretical* thickness -- the melt velocity boundary layer his Eq. (2)
    predicts -- so what the liquid depth decides is whether that layer can physically form in the melt that is present.
    Girin (2017) Sect. 1: when "the aerodynamic interaction dominates over the rate of heating, the mass loss overtakes
    fusion", the molten layer is thinner than delta_m, "the unstable disturbances are then affected by the stabilizing
    influence of the meteoroid rigid core", and "the results obtained in Girin & Kopyt (1994) are thus valid for the case
    of dominant ablation" -- **regime 1**, this model's thin branch. When melting outstrips the mass loss the layer does
    form, layer > delta_m, and the profile becomes "the completed boundary layer ... consisting of free conjugated
    boundary layers, in air, delta_a, and in the melt, delta_m, independent from the solid core" -- dominant fusion.
    **Stage two, on the kinematic viscosities, only where the layer formed.** Sect. 2: the "classical Kelvin-Helmholtz
    type ... is in action only when the liquid kinematic viscosity is greater than that one of the gas: nu_m > nu_a",
    because then "the liquid boundary layer is much thicker than the gas one" and V_s << V_a, so "the velocity profile is
    close to the discontinuous one (tangential discontinuity) which is the base for the Kelvin-Helmholtz mechanism" --
    **regime 2**. Under the inverse inequality nu_a > nu_m the thicknesses become comparable, delta_m = O(delta_a) and
    V_s = O(V_a), the profile is "inflated", and the mechanism is "completely different from the Kelvin-Helmholtz one",
    being set by his Eq. (1) -- **regime 3**, the gradient instability the thick branch already solves through
    `dispersion`. So the plan's thin and thick branches were Girin's regimes 1 and 3 all along, and 2 was the empty cell.
    **What the branch computes, and the one choice in it that is not Girin's.** Wavelength and growth time are Girin &
    Kopyt's Eqs. (11) and (12), lambda* = 1.5 M_e Sigma/(rho_e u_eff^2) and tau* = 0.798 lambda*^1.5 (rho_l/Sigma)^1/2 --
    the same mode the thin branch uses, because it is the same tangential-discontinuity profile -- with their rigid-wall
    factor cth(Lambda) of Eq. (5) going to 1, the deep-film limit, since the layer has formed and the core no longer
    stabilises it. Neither paper gives a mass-loss rate for this cell, so it takes the torus form Girin (2017) uses on the
    thick branch, mdot = rho_l pi r^2/(lambda* tau*) with r = k_r lambda*, which like that rate and unlike the thin
    branch's does not reference the film depth: "uncapped" means *not depth-limited*, not larger, and it comes out
    **0.73 x** the thin branch's rate at the same wavelength (pi k_r^2 = 0.0908 against 1/8) wherever that branch is not
    itself depth-limited. What bounds the release is the film mass, through dm = min(mdot A dt, m_f), as on every branch.
    That rate transfer is the only modelling choice in the branch, and it introduces no new closure constant.
    **This body is in regime 3, measured.** Liquid aluminium has nu_melt = 5.42e-7 m2/s against an edge nu_gas of
    **0.0178 to 2.04 m2/s** over the 100 mm flight to 110 s and 0.0951 to 3.75 m2/s over the whole 50 mm flight, so
    nu_gas/nu_melt runs **3.3e4 to 6.9e6** and never approaches one: wherever the layer forms the mechanism is his
    gradient instability, which is exactly why that is the right model for this work. Directly: **0 of 533 960** windward
    wet patch-steps on the 100 mm flight and 0 of 42 268 on the 50 mm one take regime 2. (Those flights split 120 682
    regime 3 / 411 949 regime 1 / 1 329 front-surface Rayleigh-Taylor and 40 990 regime 1 / 1 278 Rayleigh-Taylor with no
    regime 3 at all -- the 50 mm sphere being regime 1 throughout is the wall-Knudsen gate denying Girin's closure, as
    fact 32 has it.)
    **A melt viscous enough for regime 2 would fail stage one on this body**, which is his ordering working rather than a
    gap in it. Eq. (2) carries the same ratio, delta_m/delta_a = (alpha/mu^2)^(1/3) = ((nu_melt/nu_gas)^2
    rho_l/rho_e)^(1/3), so at the viscosity threshold itself the predicted layer is already (rho_l/rho_e)^(1/3) = **151**
    air boundary-layer thicknesses, about **1.06 m** at delta_a = 7 mm. Girin records the same effect for his stony
    variant -- "the thick boundary layer, which for a high-viscosity stony melt becomes comparable with the radius of the
    meteoroid remnant" -- and it is why that case ends with lambda_f above R_0 (fact 10). No film a 100 mm sphere can
    hold accommodates a metre-scale delta_m, so such a melt is routed to regime 1, correctly. A unit test pins that
    routing so the branch is not later read as dead code that ought to have fired, and the branch-behaviour test
    therefore supplies delta_m directly.
    **Inert for every run in this plan**, measured against the model's own reproducibility floor rather than asserted: on
    the 50 mm physics flight the worst column anywhere in the 416-row history moves by **3.1e-8** between the build
    without the branch and the build with it, against **1.6e-8** between two runs of the *same* build -- the same order,
    so it is not distinguishable from run-to-run scatter -- while `thick_branch_fraction` is identical in all 67
    comparable steps (maximum absolute difference exactly 0), demise time is identical to the step at 207.5 s, and
    sprayed mass agrees to twelve significant figures. It is not bit-identical and is not claimed to be: the extra array
    operations perturb the last bit, which the adaptive integrator amplifies into 12 more right-hand-side evaluations out
    of 11 702. Where it shows up, checked rather than assumed: the source table's `branch` column (`SOURCE_COLUMNS` index
    8) and nowhere else -- the surface VTK carries surface_flow's `closure`, not the spray branch; no history column is
    added; and `thick_branch_fraction` is `film.lubrication`'s thick flag over the whole formed-layer set, so it would
    count a regime-2 patch as regime 3. That column records whether the layer formed, not which mechanism took the patch.
13. **Columns and files.** History adds `removed_mass_kg`, `film_thickness_max_mm`, `film_thickness_mean_mm`, `nose_radius_mm`, `transverse_radius_mm`, `fitted_nose_radius_mm` and `n_dead_elements` to spec §10's list (`runoff_mass_kg` = mass that arrived on another patch, cumulative); `melt_front_depth_max_mm` is the depth of the deepest element with f_l > 0 (the solidus front for the range material); `film_T_max_K`, `film_T_mean_K`, `film_frozen_fraction` (the share of the film sitting on patches below the feed ramp — mass the enthalpy calls solid that the model still treats as liquid, fact 25), `unapplied_load_J` (the deferred melt energy still queued, fact 26) and `film_blob_fraction` (the share of the film deeper than its patch is wide, fact 27), and `molten_depth_max_mm`, `molten_depth_mean_mm`, `delta_m_mean_um` and `thick_branch_fraction` (the contiguous liquid layer, the conjugate depth, and the share of wet windward patches on Girin's thick branch, fact 28) come with the film's temperature and the layer. The source table has 22 columns (`spray.SOURCE_COLUMNS`), 5e5 rows for the 100 mm physics flight (`particles.npz`, compressed). The CLI gains `--k-scale` (verification device), `--size-feedback current|initial`, `--rt-spray on|off` (fact 35) and `--consistent-mass` replaces `--lumped-mass`.

---

## File structure

```
reentry_model/mesh.py                     + prism layers (radial projection, prism split), active set with a face table, deactivate, box mesh, SurfaceMesh.owner/face_ids/tangent_from/projected_area/edges
reentry_model/material.py                 + latent heat, solidus/liquidus, MELT_RAMP, feed_fraction, LiquidProperties, MATERIAL_NAMES, exact enthalpy with the latent slope, enthalpy_liquid/enthalpy_mixed/cp_mixed/temperature_from_enthalpy_mixed (the film's liquid branch)
reentry_model/data/materials/AA7075.json, AA7075_range.json   DRAMA's drama-AA7075 verbatim + liquid properties (+ the alloy's range)
reentry_model/thermal/__init__.py         + StepResult.Q_extra/Q_dropped, set_fractions/element_energies/set_film_mass/nodal_capacity in the protocol
reentry_model/thermal/skfem_backend.py    nodal enthalpy, lumped capacity, tangent Newton mapped through the node's own material+film enthalpy, phi_e, pinned nodes, nodal loads, film mass on the boundary nodes, CG fallback
reentry_model/thermal/fenicsx_backend.py  the same scheme: stiffness in UFL, everything nodal in numpy/PETSc, MatZeroRowsColumns for pinned/Dirichlet
reentry_model/dispersion.py               Girin Eq. (1): batched cubic roots, fastest mode, cached table data/girin_dispersion.json
reentry_model/gas.py                      + GasState.s/a/m_bar, EquilibriumAir.expand
reentry_model/surface_flow.py             three-branch gate, Newtonian+Prandtl-Meyer wall pressure, edge state per 1 deg bin, Ranger boundary layer, wall Knudsen -> melt closure, shear, G
reentry_model/film.py                     lubrication branches, Runoff (edge geometry, coefficients, linearly implicit transport)
reentry_model/spray.py                    melt_layer (Girin closure only), thin-film and RT modes, SprayModel (branches, release), source_rows, histogram
reentry_model/girin_case.py               Girin-as-published flight of his Table 1 variants
reentry_model/body.py                     + MeltSettings, fit_sphere, MeltingBody (netted feed/freeze, film temperature and energy, spray, death cascade, hand-over, nose-cap fit, destination-booked accounting with bounded deferred loads, stats), reference_area/reference_length/nose_radius hooks
reentry_model/trajectory.py               + body.reference_area() and reference_length() in the drag and Kn, zero drag for a consumed body
reentry_model/aero.py                     + SESAM's Mach-1 drag step smoothed over Ma 0.98-1.02, + shape_factor on the continuum entry
reentry_model/data/atdb_disc.json         ATDB_CYLINDER at zero angle of attack: the flat-disc endpoint of the shape family
reentry_model/coupled.py                  + MELT_COLUMNS, state passed to advance, body.nose_radius() in the heating, demise, melt_results, melt VTK fields, write_particles
reentry_model/sesam_io.py                 + Reference.mass/thickness
reentry_model/compare.py                  + has_melt, melt_metrics (interpolated 1 %-mass crossing), plot_melt (7 plots)
reentry_model/viz.py                      + emitting-patch overlay, film frame/video, liquidus/solidus iso-lines, melt stills
reentry_model/cli.py                      + --melt and the Step 3 group, --k-scale, --consistent-mass, MeltingBody wiring, writers, summaries
data/reference_runs/sphere_d100.00mm_..._h077.500km_nowind.{csv,json}, sphere_d050.00mm_..._h115.000km_nowind.{csv,json}   Task 13
data/reference_values/girin2017_table1.json, girin1994_tables.json
analysis/girin_reference.py, analysis/melt_verification.py, analysis/melt_sensitivity.py
tests/test_reentry_model_{dispersion,surface_flow,film,spray,girin,melting,reference_melt}.py   new
tests/test_reentry_model_{mesh,material,thermal,coupled,cli,compare,viz,data,aero,fenicsx}.py   extended (aero: the smoothed step and the melting references)
README.md, docs/model_assumptions.md (§9), the spec (§18), sesam_verified_facts.md (§17)   Task 15
```

37. **The surface frames carry the step's own flow and spray fields (2026-09-27).** `MeltingBody` stores `last_flow` and
    `last_spray` at step (iii) of the melt step, and the element deaths of step (v) then rebuild the surface, so the
    writer's old test -- the spray arrays have as many entries as the surface has patches -- failed after every step
    with a death, which is every melting step: on the 50 mm physics flight **0 of 104** surface frames carried a
    non-zero `p_w`, `tau` or `we_s`, with `closure` 1 and `kn_local` NaN throughout. The body now keeps the face ids of
    the surface it evaluated (`last_face_ids`), and `on_current_surface(values, fill)` carries each array onto the
    current surface by face id, with the old defaults on faces the deaths exposed. `test_coupled_melting_run_and_writers`
    pins it on a frame whose step had a death: the wall pressure is positive, never above the stagnation value the
    history records from the same evaluation, and falls with the angle from the flight direction (Spearman rank
    correlation below -0.9, which a mis-mapped array fails). Output only, measured rather than asserted: the 50 mm physics
    flight before and after the fix differs by at most **1.0e-8** (relative) in any history column, against the 1.6e-8
    between two runs of one build, and its frames now carry wall pressures to 3.2 kPa and shear of 66-160 Pa over the
    melting phase.

## Amendment of 2026-09-27 — surface recession and remeshing (facts 38–45)

Design: `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`, approved 2026-09-27. It amends
spec §5 and §10 and therefore sub-plans 01, 09, 10, 13, 14, 15, and adds Tasks 16 and 17. Facts 1–37 stand except
where these say otherwise. **Fact 4's `PHI_DEATH = 0.05` is superseded by fact 44 (0.50).** Numbering starts at 38
because fact 37 (`on_current_surface`, 2026-09-27) is immediately above: it was regenerated into
`prototype/plan3/plan.md` but never copied over the committed master plan, so it is carried here instead.

38. **gmsh cannot refine locally, and subdivision is the wrong mechanism anyway.** gmsh 4.15.2's only refinement
    entry point is `gmsh.model.mesh.refine()` — "refine the mesh of the current model by **uniformly** splitting the
    elements". No argument, no subset, no target size, and no parent-to-child map with which to carry temperature or
    `φ_e` across a split. Size fields drive generation from scratch, not refinement. Nor would hand-written local
    refinement help: bisection produces children similar to the parent, so reaching 0.25 mm from an 8 mm core
    tetrahedron needs five levels, about 32 768 children per parent. Adaptivity pays when the refined band is a small
    fraction of the domain; here the front sweeps the whole body (fact 14: the nose recedes 77 mm on a 50 mm radius).

39. **The graded core is much coarser than the skin, and that is what the surface becomes.** Measured on
    `sphere_R50.000mm_hs2.000mm_hc8.000mm.msh` (18 896 nodes, 87 632 tets), equivalent regular-tetrahedron edge by
    shell: **2.74 mm** at 45–50 mm radius, 3.67 at 40–45, 4.44 at 35–40, 5.24 at 30–35, 6.65 at 20–25, 8.26 at 10–15,
    8.32 at 0–5. Depth-averaged over the recession the present mesh therefore recedes in jumps of **about 6 mm**, not
    the 0.25 mm of the outermost prism layer — which covers only 3.75 mm, i.e. **4.9 %**, of the 77 mm the nose
    travels. Fact 28's "the coarse core *is* the surface" by 82 s is this, observed in flight.

40. **The death boundary is already a valid closed surface, and smoothing it is exactly correctable.** On a
    deliberately ugly eroded body (8 174 live tets, 3 404 boundary triangles) every one of the 5 106 boundary edges is
    used by exactly two triangles and the area-weighted normal sum is **4.05·10⁻¹⁶** of the total area. So no surface
    *reconstruction* is needed, only smoothing. Eight Taubin λ/μ passes (λ 0.53, μ −0.55) with each displacement
    clamped to 0.35 of the mean edge length moved the enclosed volume by **+0.283 %**; one uniform offset along the
    area-weighted vertex normal of −19 µm restored it to **1.1·10⁻¹⁴**. Taubin is purely local, hence topology-agnostic:
    any genus, any number of components, no shape assumption — which is what delivers the arbitrary-shape requirement.
    A level-set/marching-cubes reconstruction was rejected as general in the *wrong* direction (it would weld a
    thinning lens shut or erase a fragment); a radial or spherical-harmonic fit was rejected because it assumes a
    star-shaped body.

41. **gmsh will mesh the smoothed body, and the remesh cost is sub-linear.** Discrete surface →
    `classifySurfaces(π, True, True, π)` → `createGeometry()` took **0.1 s** and produced **2** surface entities;
    smoothing first is what keeps this tractable, since the feature detector on a raw staircase produces thousands.
    Volume generation gave 0 inverted elements and a volume **0.091 %** below the smoothed target — gmsh
    re-triangulates the surface at its own size and cuts corners, worth 1.3 g per remesh on the 1.472 kg sphere, so it
    must be corrected by an outward uniform-normal offset. Cost is entirely in `generate(3)`: HXT and Delaunay are
    indistinguishable, `Sampling` 200 vs 20 makes no difference, and threads give 45.8 s (1) → **29.4 s (4)** → 29.7 s
    (8), saturating at four. Scaling at 4 threads: 137 248 tets 29.5 s, 376 183 tets 48.7 s, 748 606 tets 75.1 s,
    i.e. **≈18 s fixed plus ≈76 µs per tetrahedron**; peak memory 0.5 GB at 749 k tets.

42. **Node motion is safe to about 40 % of the local cell and collapses past about 70 %.** Pushing every boundary
    node of the 2.74 mm mesh inward: 0 inverted out to 1.25 mm; worst-element radius ratio 0.289 → 0.232 (0.5 mm) →
    0.076 (1.0 mm) → 0.022 (1.25 mm); at 1.5 mm two elements fall below 5 % of their volume; at **2.0 mm, 1 319
    elements invert** and 4 908 are slivers. Mean quality near the surface degrades far more slowly than the worst
    element (0.761 → 0.711 at 1.0 mm), so the kill rule must be **per element**, and the node motion must sub-step
    when a step would exceed `RECESSION_SUBSTEP_FRAC` = 0.35 of the local cell.

43. **Per-patch volume cannot be matched exactly by nodal motion; the total always can.** Euler forces a closed
    triangulation to have twice as many triangles as vertices — measured **18 078 patches against 9 041 boundary
    nodes**, and 2·9 041 − 4 = 18 078 exactly — so the per-patch volume constraints outnumber the nodal unknowns
    **2:1** and are generically inconsistent. Least-squares fit quality depends entirely on how smooth the recession
    field is: a **smooth** (Lees-like, √cos θ) field gives median per-patch error **0.03 %**, 90th percentile 0.56 %,
    leaving `φ_e` a residual of 10⁻⁴ of an element; a **patchy** field (the spray gate flipping) gives median 22 % and
    a residual of 7.3 % of an element; white noise 29 %. The **total** is within 0.017 % before correction in every
    case and exactly zeroable after. How smooth the real field is has **not** been measured and is the first number
    the implementation should obtain.

44. **Frozen geometry costs about a factor of two in the conduction, and `PHI_DEATH` should rise to 0.50.** 1-D
    ablating AA7075 bar (ρ 2813, k 170, c_p 1100 from `AA7075_nomelt.json`), backward Euler, P1, lumped capacity,
    prescribed recession and conducted flux, 2.16 mm against a 0.05 mm reference. Typical (0.40 MW/m², 0.10 mm/step,
    surface 312 → 600 K): moving nodes err by mean +0.89 K / max 2.43 K, frozen geometry by mean **+1.97 K** / max
    6.00 K. Peak flux (2.00 MW/m²): +4.44/12.14 against **+9.87/29.99**. Peak recession (1.00 mm/step): +0.62/4.13
    against +1.59/7.48. Frozen geometry is consistently ≈2.2× the mean and ≈2.5× the maximum error, and **both are
    dominated by the 2.16 mm cell**, not by the design choice. Against `PHI_DEATH` the frozen-geometry mean error is
    0.05 → +2.97 K (typical) / +14.83 K (peak); 0.25 → +2.27/+11.34; **0.50 → +1.97/+9.87**; 0.75 → +2.37/+11.84 — a
    shallow optimum at a half, and the plan's current 0.05 is the worst of the four. **`PHI_DEATH` becomes 0.50**,
    which also matches what fact 42's quality envelope wants. Consequence to check: the film now receives half an
    element at a time rather than a twentieth (fact 25's cap, `film_blob_fraction`).

45. **The backends behave oppositely under a moving mesh, which is why the solver is left alone.** skfem's speed
    comes from rescaling matrices precomputed in `element_matrices`, valid only while the geometry is fixed: 87 632
    tets assemble in **0.004 s** fixed but **0.074 s** (17.4×) if the geometry must be recomputed. dolfinx recomputes
    Jacobians at every assembly anyway, so moving **every** node changed its time not at all: **0.031 s** either way.
    At 310 k tets: skfem 0.015 s fixed, 0.262 s moving; dolfinx 0.11 s always. `mesh.geometry.x` is writable in place,
    and `dolfinx.fem.create_interpolation_data`, `geometry.bb_tree` and `compute_colliding_cells` exist in the
    installed 0.11, so FEniCSx already has the mesh-to-mesh transfer a remesh needs. Because the adopted design does
    **not** move the solver's mesh, none of this is incurred per step; it is recorded for the deferred Design A and
    for the remesh, where both backends must be re-set-up.


## Amendment of 2026-10-02 — the liquid below the conjugate depth runs off, and the frames carry δ_m (facts 46–53)

Requested by Asha on 2026-10-02 for the large-fragment (Spheral) model that will read this model's exported frames (her
three-zone rule, fact 46). It amends sub-plans 06, 07, 09, 10, 13, 14, 15 and 17; sub-plan 09's amendment holds the
design and the tested code. Facts 1–45 stand except where these say otherwise. **Fact 29's "one consequence not
implemented" and fact 34(c) are implemented by fact 46.** Fact 28(b)'s feed gate is unchanged: the deep runoff moves the
liquid the gate holds back, it does not release it to the spray. **Fact 31's reproducibility statement is revised by
fact 52, and fact 12's Δt/2 sensitivity is extended by fact 50.** Fact 13's column list gains the five columns of
fact 46. Measured in a throwaway copy of `prototype/proto3/` as it stood on 2026-10-02 — Step 3 through the amendments
of 2026-09-25, fact 37, and sub-plan 01's dense band and derived-surface plumbing (177 363 tetrahedra and 18 830 surface
patches on the 100 mm sphere); not sub-plan 02's material amendments, not sub-plan 09's changes of 2026-09-27
(`PHI_DEATH` 0.50, the shape consumers on the derived surface) and not Tasks 16–17, which exist only as plans — with
`AA7075_range`, US76, the physics heating and every other setting at its default.

46. **The liquid below the conjugate depth runs off and is never sprayed (decided 2026-10-02).** Asha's three-zone
    rule: (1) liquid above the liquidus within Girin's conjugate depth δ_m of the surface is the sprayable skin and
    stays in the film; (2) the contiguous liquid below δ_m, down to a film limit of about 2–3 mm, is not sprayed but
    runs off under the pressure gradient along the surface and the deceleration, and stays in this model; (3) material
    thicker than the film limit goes to the large-fragment model. Until now zone 2 never moved: fact 28(b) holds it in
    its elements. **The design** (sub-plan 09): a second liquid account per patch, `m_d`, and one new stage after the
    feed. On windward patches under Girin's closure — the only place a conjugate depth exists — the deep liquid of a
    patch is the fully liquid inventory of the elements of its contiguous molten chain that the feed gate held back
    (shared by patch area among the chains that pass through them) plus its `m_d`; it moves by
    `film.deep_flux` = G ((b + h_D)³ − b³)/(3 μ_l), fact 29's pressure- and deceleration-driven part over the whole
    liquid depth, with h_D the deep liquid's thickness by mass and b the film's, on the film's own linearly implicit
    upwind transport; a patch that loses deep liquid gives it up in proportion from its held elements (never below
    `PHI_MIN`) and its `m_d`; what arrives goes to `m_d`. `m_d` is never offered to the spray: it becomes film only from
    the top, as far as the film is thinner than δ_m — at most one skin's worth per macro step — and all at once where
    there is no conjugate depth. It rides the patch nodes with the film, holds the liquid enthalpy, freezes back first,
    and is handed over at an element death by the film's rule. Every transfer is booked by facts 5 and 25 (sub-plan
    09 sets this out against the two rejected implementations). **Alternatives not taken:** mobilising the whole held
    inventory every step (lifts the pool out of the mesh and puts new melt above a non-sprayable pool),
    element-to-element transfer (the receiving elements of a pool are full and the interior may not empty), a
    non-sprayable flag on the film account (every film consumer would need it), surfacing at the edge of the in-situ
    pool (would spray liquid however deep it lay) and a two-directional depth partition of the film (would change thick
    patches where no liquid lies below δ_m). **What does not change:** `lubrication`, the feed gate and the spray
    module; with no liquid below the conjugate depth a melt step is bit-identical with the deep runoff on and off (unit
    test), and so is the whole 50 mm flight (fact 49). **New outputs:** history columns `deep_liquid_kg` (the deep
    liquid the step saw), `deep_mass_kg` (the deep account), `deep_runoff_mass_kg` (cumulative mass taken out of the
    elements), `deep_surfaced_mass_kg` (cumulative deep liquid that became film) and `deep_blob_fraction` (the share
    of the deep account deeper than its patch is wide); surface-frame fields `delta_m` (fact 51) and `deep_thickness`;
    CLI `--deep-runoff on|off`, default on, `off` reproducing every earlier run (fact 49), run names ending in
    `_deeprunoff-off` when off.

47. **The deep flux is strong, so the deep liquid ends where the flow converges — and on this surface that is a few
    hundred crater facets.** Lubrication gives a layer of liquid aluminium (ν_l = 5.4e-7 m²/s) a mean velocity
    G h²/(3 μ_l): 10 m/s for a 1 mm layer at G = 4e4 Pa/m. The layer's Reynolds number u h/ν_l is then 10³–10⁴, so the
    laminar flux is an overestimate (turbulent wall friction would give about 2 m/s for 1 mm), but either way a
    millimetre layer crosses a 2 mm facet in about a millisecond and the transport reaches its steady state within every
    0.5 s step: measured in a unit-test setting at 69.8 km, a millimetre of deep liquid spread over all 2 296 windward
    patches is gathered within one step onto so few patches that only 100 still hold deep liquid once each has filled
    its skin. On the flight G changes sign: the median G on the patches whose liquid is deeper than δ_m is positive
    until about 78 s — typically +5 to +15 kPa/m, +22 kPa/m when the gate first opens (outward, the pressure gradient
    winning) — and negative after it, typically −10 to −50 kPa/m and down to −100 kPa/m (toward the nose, the
    deceleration winning), so the liquid is pushed to wherever the field converges. On the staircase surface left by
    element death those places are the element-death craters of the eroding front: 88–100 % of the deep account lies
    within 5 mm of the front-most point of the body, about half of it on backward-facing crater walls (facet normals
    beyond 90° from the flight direction, where it arrives by the death hand-over), on 600–930 of the 7 700–17 000
    surface facets. From 60 s on, 96–99.6 % of it sits on facets where it is deeper than 2 mm and 92–99 % where it is
    deeper than the facet is wide (`deep_blob_fraction` has a median of 0.96 over the steps that hold any), with
    equivalent depths up to 0.69 m on facets of 2.5–52 mm². By the three-zone rule that is zone-3 material, but its
    location is the mesh's: these are fact 27's undrainable facets at a larger scale. Two physical limits on piling are
    missing from the model: levelling by surface tension below about 2 cm (on the front face the normal component of the
    deceleration is destabilising — the Rayleigh–Taylor mode — so only capillarity levels), and stripping of a pile's
    exposed surface beyond its facet's area. Task 16's derived surface, which takes θ and the film tangents from
    smoothed normals, is expected to remove most of these sinks; the pile-up must be re-measured after it (Task 14's
    amendment). The piles form only at the default step: at 0.25 and 0.125 s there is too little liquid below the
    conjugate depth to pile (fact 50).

48. **The runoff transport conserves mass only to the conditioning of its direct solve.** Measured on a 3 mm molten pool
    at 69.8 km: the linearly implicit transport's total changed by 1e-14 to 2e-11 of the mass it moved per call, for the
    film transport as well as the deep one, because the coefficient matrix is stiff (dt × c, the fraction of a patch's
    liquid an edge would carry in one sub-step, reaches about 10³ for a millimetre layer and far more for the piles of
    fact 47) and `spsolve`'s error grows with the matrix's condition number. For the film that is at most 1e-13 kg per
    step and invisible; for tens of grams of deep liquid it reached 7.5e-12 kg (5e-12 of the body) in six steps and
    broke the 1e-12 mass test. The deep stage therefore scales its arrivals to its departures — a correction of that
    size — so its books are exact by construction. The film transport is left as it was, because correcting it would
    perturb every existing run at that level.

49. **What the deep runoff does to the flights** (physics mode, everything else default; both copies run with numpy's
    random seed fixed, fact 52, so the differences are the amendment's and nothing else's). **50 mm, whole flight:
    unchanged, bit for bit** — the sphere never has Girin's closure (fact 32), so there is no conjugate depth and no
    deep liquid; all 79 shared history columns in all 415 rows, all 22 source-table columns of 26 922 rows and 48 of the
    49 result fields agree exactly, the 49th being the run time (+1.0 %, 220.9 to 223.1 s: the deep stage's march costs
    that even when it finds nothing). **100 mm to 120 s, flag off: unchanged, bit for bit** (all 79 history columns in
    all 241 rows, all 150 352 source-table rows and every result field but the run time agree exactly).
    **100 mm to 120 s, deep runoff on**: the deep runoff takes **0.285 kg** out of the elements over the 120 s — liquid
    the feed gate held below the conjugate depth, of which the unamended model holds a median of 4.8 g (at most 7.0 g)
    at any one time, in 149 of the 240 steps from 49.5 s, when Girin's closure first applies. Of that, 0.222 kg became
    film from the top and was mostly sprayed, and 0.055 kg is still deep at 120 s (0.117 kg at the peak, 96.5 s).
    **Sprayed mass rises from 1.018 to 1.091 kg (+7.1 %) and the body at 120 s is 16 % lighter (0.454 to 0.381 kg)**,
    against a spread of 0.15 % and 0.33 % between two runs of the unamended model with different random states (fact 52)
    — about fifty times the scatter. It is not a change in how melt is sprayed: the energy carried away per kilogram is
    1.0624 against 1.0647 MJ/kg (−0.2 %), the median radius of all the droplets released (by number) is unchanged (180.2
    to 180.3 µm) and the droplet count rises 3.1 % (1.72e7 to 1.77e7). It is a change in where and when liquid reaches
    the surface: the owner feed falls 15 % (1.009 to 0.856 kg; 27 % over 25.5–80 s), because liquid that used to wait in
    its element until the surface reached it is now drawn from below and the drained elements then die without needing
    to be melted; the absorbed heat rises 1.8 % and the rest of the extra mass loss is heat the lighter body no longer
    stores. The front-surface Rayleigh–Taylor mode releases 37 % more (53.0 to 72.5 g), because the piles deepen the
    layer its criterion sees, and the largest droplet grows from 4.38 to 5.66 mm in radius (+29 %). Re-solidified mass
    rises from 6.2 to 11.2 g, the thick-branch share of wet windward patches from a median 12.4 % to 15.5 %, and the
    deepest contiguous molten layer falls from 20.9 to 17.6 mm (its median per step from 10.55 to 10.46 mm). The energy
    balance stays exact (−1.1e-10 of the absorbed heat, against −8.9e-12), Newton takes 3.01 iterations per step against
    3.00, and the cost is within the machine's noise: 584 s against 623 s seeded and 585 s against 570 s unseeded, while
    two runs of one build differ by 9 %. **The mass reaching the equator** (the Step 4 plan's fact-1 quantity), measured
    on the equatorial ring — the windward patches that border a leeward one *and* lie at least 0.9 of the transverse
    radius from the flight axis, because the adjacency test alone also catches the rims of the craters on the eroded
    front: while the equator is intact (25.5–80 s), 22.9 g of film reaches it from 570 g of melt without the deep runoff
    (4.0 %), and with it 21.7 g of film plus 10.8 g of deep liquid from 577 g of melt delivered (418 g by the owner feed
    and 159 g taken from below the conjugate depth; 5.6 %, i.e. +41 %). Over the 120 s it is 45.9 g against 56.6 g of
    film plus 54.2 g of deep liquid, and the mass sprayed on the ring rises 24 % (72.6 to 90.4 g). Almost none of it is
    still there at the end of any step (at most 0.07 g of film on the whole ring, and no deep liquid) — the equator
    sprays what reaches it (the Step 4 plan's fact 4) — so no rim forms in Step 3 either way. These fractions are five
    to seven times the Step 4 plan's fact 1 (0.8 % of the melt for the 100 mm sphere and 1.5 % for the 50 mm one, over
    the same windows, 25.5–80 s and 174.5–194.5 s, but measured before the dense band replaced the prism layers; the
    50 mm flight here gives 5.3 g of 69.4 g, 7.7 %), so that fact should be re-measured on the current prototype with
    its definition stated. At smaller steps the deep runoff's increase disappears (fact 50). **The whole flight, to the
    ground** (seeded, default step): without the deep runoff the 100 mm sphere lands at 582.5 s with 0.453 kg, 30.8 % of
    its initial mass (fact 24 measured 29.9 % on the prism-layer mesh); with it, at 602.7 s with 0.330 kg, **22.4 % — a
    quarter less** — because the 55 g still deep at 120 s surfaces and sprays by 229.5 s (spraying ends at 225.0 s
    without it), the sprayed mass reaching 1.142 kg against 1.019 kg (+12 %), and 21.0 g re-solidifies instead of 6.6 g.
    Every gram the deep runoff mobilised is accounted for: 285.1 g taken from the elements, 277.1 g became film from the
    top and the remaining 8.0 g froze back (freeze-back being the only other way out of the deep account); the deep
    account is empty on landing and the energy balance exact (−8.8e-10 of the absorbed heat). The median droplet radius
    is again unchanged (180.2 µm) and the run takes 1 897 s against 1 956 s. At the default step the deep runoff
    therefore moves fact 24's headline by a quarter — which fact 50 shows to be a time-step artefact (at 0.125 s it
    moves the sprayed mass by −0.16 %). **Both thermal backends** with the deep runoff active (a 6 mm molten pool at
    69.8 km on the coarse mesh, six coupled steps, each backend seeded): mass agrees to 1.3e-9, sprayed mass to 1.4e-8,
    the deep account and the mass taken from the elements to 9e-7, the film to 1e-8, the number of active elements
    exactly (7 069), the temperatures to 0.021 K, and both energy balances are exact (8.7e-10 and 2.2e-11). That is
    looser than fact 12's 1e-10 and 0 K, most likely because the stiff deep transport carries the backends' last-bit
    differences further (not isolated); the FEniCSx test file, whose melting check starts in the merged branch where
    nothing runs deep, passes with the amended code (9 passed in `fenicsx_env`).

50. **The liquid below the conjugate depth is a time-step artefact on this flight: it vanishes as the macro step
    shrinks, and the deep runoff's effect vanishes with it.** Measured on the 100 mm flight to 120 s at 0.5, 0.25 and
    0.125 s, every other setting unchanged and every run seeded. The liquid the unamended model holds below the
    conjugate depth falls from a median 4.75 g (at most 7.0 g) at any one time at 0.5 s to 0.23 g (at most 1.65 g) at
    0.25 s and 0.004 g (at most 0.08 g) at 0.125 s — integrated over the flight, from 295 to 22.9 to 0.15 g s. The
    reason is the order of the step: the conduction runs first and the melt step only then hands the wall-owning
    element's liquid to the film and lets that element die, so between two melt steps the heat flux superheats the
    molten surface and melts the material beneath it, and the surface can recede through molten material by only one
    element per macro step (an exposed element dies only once it is consumed, and it is fed only at the next step's
    feed); whatever melts below the owner waits there, and a shorter step drains it sooner. That backlog is what the
    deep runoff moves: it takes 285 g, 38 g and 0.53 g from the elements over the flight at the three steps, the deep
    account peaks at 117 g, 3.0 g and 0.014 g, and the piles of fact 47 form only at the default step (a median 96 %,
    7 % and 0 % of the deep account deeper than its facet is wide). Its effect on the results falls accordingly:
    sprayed mass +7.1 %, −0.3 % and −0.16 %, mass at 120 s −16 %, +0.65 % and +0.36 % — at 0.125 s within the spread
    of two runs of the unamended model with different random states (0.15 % and 0.33 %, fact 52). Only small,
    threshold-sensitive quantities still move at 0.125 s: the front-surface Rayleigh–Taylor release (8.5 to 11.4 g,
    +33 %, where two random states differ by 3 % at 0.5 s; it is 0.8 % of the sprayed mass at this step) and the
    largest droplet (+10 %). The equator receives 4.0 %, 4.6 % and 4.7 % of the melt without the deep runoff, and
    5.6 %, 4.7 % and 4.7 % with it.
    **The unamended model is not converged in the step either, for the same reason.** The contiguous molten layer that
    the thick/thin branch test reads (fact 28) is the same backlog: its mean depth over the surface (median over the
    steps) falls from 1.26 to 0.28 to 0.066 mm, and under the patches whose liquid is deeper than δ_m from 1.86 to 0.76
    to 0.59 mm, while δ_m itself stays at 289–295 µm. The share of wet windward patches on Girin's thick branch
    therefore falls from 12.4 % to 8.0 % to 2.6 % and the thin branch's share of the sprayed mass rises from 9 % to
    18 % to 38 %; because that branch makes small droplets (a median radius by number of 70–78 µm, against 176–186 µm
    on the thick branch), the droplet count rises from 1.72e7 to 2.82e7 to 5.85e7 and the median radius by number falls
    from 180 to 105 to 76 µm (by mass only from 189 to 188 to 176 µm). The front-surface Rayleigh–Taylor release falls
    from 53.0 to 17.2 to 8.5 g and the re-solidified mass rises from 6.2 to 16.8 to 48.4 g. The sprayed mass (1.018,
    0.991, 1.018 kg) and the mass at 120 s (0.454, 0.481, 0.454 kg) move by up to 3 % and 6 % without a trend. This
    extends fact 12's Δt/2 sensitivity (−21 % in median radius on the 50 mm flight) and is the larger finding of this
    amendment: the droplet population Step 3 hands to the wake depends on the macro step through the melt step's
    one-element-per-step recession. The cost of a smaller step is proportional to the number of steps — about 2.4 s
    of wall time per step on this machine, so 10, 20 and 40 minutes for these 120 s at 0.5, 0.25 and 0.125 s.
    **Consequences:** at the default step the deep runoff amplifies a splitting artefact into a mesh artefact (fact
    47); with a short enough step there is, on this flight, almost no contiguous liquid below the conjugate depth, so
    zone 2 of the three-zone rule is essentially empty and the frames' `deep_thickness` is close to zero, while their
    `delta_m` (fact 51) does not depend on it. Fact 28's molten layer of 350–890 µm and everything that reads it — the
    branch test, the Rayleigh–Taylor depth criterion, the droplet population — carry the same dependence. Fact 53 (a)
    is the decision this calls for.

51. **The frames carry Girin's conjugate depth per patch.** `surface_<k>.vtp` gains `delta_m` [m], the step's own
    value from `spray.melt_layer`, carried across that step's element deaths by face id like `p_w` and `tau`, and NaN
    wherever the closure is not Girin's (no conjugate depth exists there), on faces the step's deaths exposed and
    before the first evaluation; and `deep_thickness` [m], m_d/(ρ_l A). `delta_m` is written from every step that
    evaluated the flow, unlike `closure`, `p_w` and `tau`, which come from the spray step and so only from steps with
    film. Measured on the 13 frames of the 100 mm flight to 120 s: `delta_m` is finite on exactly the
    25 650 patch-frames with Girin's closure and NaN on all the others (no mismatch either way), 143–421 µm, and NaN
    throughout until the gate opens at 49.5 s; `deep_thickness` reaches 0.69 m (fact 47) — read it as a mass per area.
    The 50 mm flight never has Girin's closure, so its frames carry `delta_m` as NaN on every patch.

52. **The model was not bit-reproducible, and the cause is pyamg's random starting vectors.** Two runs of one build
    differ from the first melting step on — by 2e-10 in φ_e and 4e-11 J in the deferred loads at that step of the 50 mm
    flight, from a temperature field that already differs before anything melts — because pyamg draws the starting
    vector of its spectral-radius estimate (which sets the smoother weights) from numpy's global generator. Ruled out by
    measurement: threaded BLAS (one thread changes nothing) and Python's hash seed. With `np.random.seed` fixed at the
    start of a run, two runs agree in all 98 recorded arrays exactly. How much the difference grows is chance: on the
    50 mm flight a repeat of the unamended run agreed to 5e-13 in sprayed mass (2e-7 in re-solidified mass), while an
    unseeded amended run (in which the deep stage finds nothing) drifted by 0.04 % in sprayed mass, 0.85 % in final mass
    and 11 % in re-solidified mass, all of it in the late collapse; on the 100 mm flight two unamended runs with
    different random states differ by 0.15 % in sprayed mass, 0.33 % in the mass at 120 s, 1.2 % in droplet count, 14 %
    in re-solidified mass and 9 % in run time. This revises fact 31 ("the model is reproducible") and very likely
    explains its unexplained 6e-5 and fact 36's 1.6e-8 floor. Not fixed here (it is outside Step 3's sub-plans' code):
    one line, `np.random.seed(0)` at the start of `cli.cmd_run` or a fixed starting vector passed to pyamg, would make
    every run reproducible; until then compare runs as seeded pairs or quote the spread beside the difference.

53. **Not done, and for Asha to decide.** (a) **The time step (fact 50) — decide this first; it concerns Step 3 as a
    whole, not only the deep runoff.** On the 100 mm flight the liquid below the conjugate depth, and with it the deep
    runoff's effect, vanishes as the step shrinks (285 g, 38 g and 0.5 g mobilised at 0.5, 0.25 and 0.125 s); but the
    molten layer the thick/thin branch test reads shrinks the same way, and the droplet population — count, branch
    split, median radius by number, Rayleigh–Taylor release, re-solidified mass — has not converged even at 0.125 s.
    Options: (1) make `--deep-runoff off` the default (one line in `MeltSettings` and one in the CLI): at the default
    step the deep runoff then no longer turns the backlog into piles, and with a short enough step it makes no
    difference on this flight anyway; (2) run melting flights at a smaller step: the cost is proportional to the number
    of steps (10, 20 and 40 minutes for 120 s of the 100 mm flight at 0.5, 0.25 and 0.125 s), and 0.125 s still does
    not converge the droplet population; (3) let the surface recede through molten material by more than one element
    per macro step — when a death exposes a fully molten element, feed it within the same step, and repeat until the
    exposed element is not fully molten — which removes the backlog at its source and is expected, but not yet shown,
    to make the molten layer, the branch test and the droplet population converge at the default step. Recommendation:
    (3), as a Step 3 amendment of its own with a time-step study of its own, before any droplet population or
    deep-runoff result is quoted; until then (1), keeping the deep runoff available and reading fact 49's default-step
    numbers as an artefact of the step, not as a physical effect. The amendment as built and tested has the default
    `on`, following the repository's pattern for a new mechanism; flipping it is this decision.
    (b) **The piles (fact 47).** They form only at the default step (fact 50). As specified — deep liquid never
    sprayed, moved by a lubrication flux with nothing to level it — the deep runoff gathers the liquid below the
    conjugate depth into the
    patch graph's sinks within each step; on the staircase surface those are crater facets at the eroding front, and
    nearly all of the liquid then exceeds the 2 mm film limit, i.e. it is zone 3 by the rule, in a place the mesh chose.
    Options: (1) keep it, and let the large-fragment model take zone 3 from the frames (`deep_thickness` as a mass per
    area: the mass is exact, the location mesh-dependent); (2) add levelling by surface tension to the deep transport
    (a graph-Laplacian term in the same implicit solve), the physics that limits such piles below about 2 cm; (3) move
    deep liquid beyond the film limit into a separate, frozen zone-3 account, or out of the body as a source for the
    large-fragment model — a decision about where zone 3 lives; (4) re-measure after Task 16 before deciding, since the
    derived surface is expected to remove most of the craters; (5) make `--deep-runoff off` the default until one of
    these is done. Recommendation: settle (a) first — at 0.125 s no piles form on this flight, so this decision is
    needed only if deep liquid survives (a)'s fix, and then (4) before (2).
    (c) **How fast deep liquid becomes skin.** One conjugate depth per macro step is a resolution choice, so the
    exposure rate scales with 1/Δt (at 0.25 and 0.125 s the deep account peaks at 3 g and 0.014 g instead of 117 g,
    though mostly because less liquid is held below δ_m in the first place, fact 50). Physically the shear
    re-establishes over newly exposed liquid in about δ_m²/ν_l ≈ 0.16 s at δ_m = 290 µm, three times faster than one
    skin per 0.5 s, and under a supercritical skin
    Girin's mode strips far faster than either. Options: keep it; tie the top-up to that renewal time (about three skins
    per step); or let deep liquid under a supercritical skin spray at the instability's rate, which is closest to
    Girin's outstripping regime but contradicts the rule's "not sprayed".
    (d) **The Rayleigh–Taylor mode and the deep liquid** (sub-plan 07's amendment): its criterion sees the whole layer,
    its release only the film. Allow it to release the deep liquid where it applies, or not.
    (e) **The Couette closure.** Where Girin's closure is denied no conjugate depth exists and nothing moves: on the
    50 mm flight, which never has his closure, the feed gate holds a median of 3.5 g and at most 15 g of contiguous
    molten liquid below the owner elements at the default step, unmoved and unsprayed as before (by fact 50 most of
    it is presumably the same backlog; not measured at a smaller step). Whether a pressure-driven runoff should
    apply there too, and what the skin is when the whole film is sheared, is open.
    (f) **"The shear passed down from above"** (Asha's description of zone 2) is not included: fact 29's form gives the
    liquid beneath the skin only the pressure- and deceleration-driven part. In a steady lubrication profile the skin's
    base would pass the wall-parallel stress τ on to the liquid beneath, adding a Couette part τ s²/(2 μ_l) for a deep
    layer of thickness s; Girin's conjugate boundary layer says the shear has not reached below δ_m. Decide which.
    (g) **The film/deep split** (sub-plan 06's amendment): the exact half-channel split would give the film
    (G/μ_l) b s (2b + s)/2 more and the deep liquid that much less; the column's total is exact either way.
    (h) **Laminar lubrication at Reynolds numbers of 10³–10⁴** overstates the deep flux (and the film's); a turbulent
    wall-friction closure would cut it about fivefold for a millimetre layer, which within a step changes how fast the
    deep liquid reaches where it collects, not where that is.
    (i) **Seeding** (fact 52): one line in `cli.cmd_run` would make every run reproducible; recommended.
    (j) **"Unchanged wherever h ≤ δ_m"** is implemented as: nothing changes on any patch, step or flight that has no
    liquid below the conjugate depth (bit for bit, tested and measured); patches downstream of liquid that does run off
    receive it, which is the point of the change.

## Amendment of 2026-10-03 — the molten cascade: the surface recedes through molten material within the step (facts 54–61)

Asha's decision of 2026-10-02 on fact 53 (a): option (3), fix the molten backlog at its cause by letting the surface
recede through more than one molten element per macro step. It amends sub-plans 07, 09, 10, 13, 14 and 15; sub-plan 09's
amendment holds the design, the alternatives and the tested code. Facts 1–53 stand except where these say otherwise.
**Fact 50's account of the cause is refined by fact 54, and fact 53 (a) is answered by facts 57–59.** Measured in a
throwaway copy of `prototype/proto3/` with the deep-runoff amendment of 2026-10-02 applied — verified before any change
to be reproduced byte for byte by that amendment's nine diff blocks — so everything fact 49 lists about the copy holds:
`AA7075_range`, US76, physics heating, the dense band (177 363 tetrahedra and 18 830 surface patches on the 100 mm
sphere), every other setting at its default, and **`PHI_DEATH = 0.05`**, because fact 44's 0.50 exists only as a plan
(sub-plan 09's amendment of 2026-09-27). Every run is seeded in the measurement harness (numpy's generator, seed 12345,
fact 52); the model itself is not seeded.

54. **What the backlog is made of: molten elements beside the melt front, not a column of liquid under an intact skin.**
    Measured at the end of the melt step at 50.5, 60.5 and 70.5 s of the 100 mm flight (cascade in its first form, deep
    runoff off), on the chains of `molten_depth`'s march below the wall-owning elements: 1 660–1 920 elements, 10–15
    elements deep at the deepest, holding 4–5 g of liquid. **Not one of them has all four nodes at or above T_feed
    (910 K)**: 43–45 % have one node below it and 50–51 % two, the coldest node at a median 897 K — 93 % liquid by the
    enthalpy, f_l = (897 − 750)/158 — while their mean temperatures run from 910 to 937 K (median 918 K); 97–98 % of
    those cold nodes are shared with a wall-owning element, and 43–55 % lie on the surface itself. Not even a
    wall-owning element is ever fully molten in that strict sense at the end of a step (0 of 14 000–17 000), though
    1 600–1 700 of them start a molten chain. The reason is the feed's own energy debit (fact 5): an element being fed
    gives up the enthalpy of its molten part and keeps the colder rest, which cools its nodes, so every node an element
    shares with an element being fed sits on or just below the feed ramp. Fact 50's "fully molten element" is therefore
    molten in `molten_depth`'s sense — its mean temperature at or above T_feed — and the chain it marches through runs
    along the melt front, through elements that touch the surface at a node or an edge without owning a face, as much as
    down into the body. A first version of the cascade that fed only elements with every node above T_feed fed 32 g in
    120 s (in 129 of 240 steps, one or two passes each) and left the backlog as it was: held liquid a median 4.41 g
    against 4.75 g, molten layer 1.21 mm against 1.26 mm, 16.7 million droplets against 17.2 million, median radius
    180 µm either way (deep runoff off; with it on, 5.9 g fed and no change beyond the run-to-run spread). It was
    replaced before any other measurement by `molten_depth`'s own definition.

55. **The molten cascade.** After each pass of the death loop of step (v), every element the pass's deaths have just
    made a wall owner that is fully molten — mean nodal temperature at or above T_feed, `molten_depth`'s test — and
    would survive the pass (φ > `PHI_DEATH`) is fed whole within the step (`MeltingBody._feed_exposed`): its remainder
    φρV leaves at its mean nodal enthalpy h_e and arrives on the faces it now owns at their liquid enthalpy, the
    difference released there; φ falls to zero, the next pass kills it, and `_kill` hands the liquid down to the faces
    it exposes like any dying patch's film. Where every node is above T_feed this is the owner feed of step (i) exactly;
    where a node or two lies on the ramp it is fact 4's death rule applied to the whole element, the residual latent
    heat of the cold nodes paid by the faces the liquid lands on (facts 25 and 26). The feed gate of fact 28(b) is
    respected by construction — only wall owners are fed, and the cascade stops at the first exposed element that is not
    fully molten — so molten material that cooler material separates from the wall is never fed.
    `MAX_CASCADE_PASSES = 32` caps the passes that feed in one step; `cascade_passes` and `cascade_mass_kg` are new
    history columns, and the run's results carry `cascade_mass_kg`, `cascade_passes_max` and `cascade_capped_steps`.
    `--molten-cascade on|off`, default on, run names ending in `_moltencascade-off` when off. **With the cascade off the
    amended code reproduces the deep-runoff amendment bit for bit** (100 mm to 120 s, seeded: all 241 rows of all 84
    shared history columns, all 167 473 source-table rows and every result field but the run time), and with no fully
    molten element exposed a melt step is bit-identical with it on and off (unit test).

56. **Tests, the two backends, and a hazard the cascade sharpens.** Five new melting tests (a molten column empties
    within one step; nothing changes when the exposed elements are not fully molten; only wall owners are fed and a
    molten core behind a cooler layer is never touched; the deep liquid is still never sprayed; exact books over six
    coupled steps), extended coupled and CLI tests and a FEniCSx comparison: unit tier 229 passed, 1 skipped, and the
    five known failures and errors from the missing melting SESAM references (Task 11); the FEniCSx file 10 passed. With
    the cascade active (a 6 mm pool at 69.8 km, six coupled steps) the backends take the same passes in every step (4,
    3, 16, 4, 5, 4) and keep the same active set; with the deep runoff on they agree to 1.3e-7 in mass, 6.3e-7 in
    sprayed mass, 2.4e-7 in cascade mass, 3.7e-6 in the deep account and 0.11 K, with it off to 1e-10 and 4e-6 K, and
    every energy balance is exact (below 7e-10). **One existing test had to be pinned to the cascade off**:
    `test_film_temperature_freeze_back_and_the_netted_transfer` puts 2.2 MW/m² on every facet — crater walls and lee
    included — of a body 50 K below its melting point, and with the cascade the whole body melts within the 8 s of
    heating (consumed at step 16, the balance exact throughout), leaving no film to freeze back. Tracing it found a
    hazard worth recording: in that device, the cascade fed two full elements (φ = 1) that shared a node with two nearly
    consumed wall owners (φ ≈ 0.11); the node was left with a sliver of heat capacity and almost no conducting
    neighbours, and the prescribed flux heated it to 2 009 K in one solve (it then lost its material and stayed pinned).
    That is fact 4's thin-owner hazard — which is what `PHI_DEATH` exists to limit — made sharper because the cascade
    removes full elements beside thin owners. On the 100 mm flight it shows no excursion of that kind: with the cascade
    the hottest node of the frames reaches 1141 K from 80 s as it does without it, and the 99.9th percentile rises by
    12–15 K (fact 57).

57. **At the default step the cascade removes the backlog — and the droplet population does not move.** The 100 mm
    physics flight to 120 s at 0.5 s, seeded, against the deep-runoff amendment's seeded runs of the same flight (deep
    runoff off unless stated). The cascade acts from 33.5 s to 108 s, in 146 of the 240 steps, feeding a median 1.9 g
    per active step (at most 3.0 g) and 245 g in all, with a median of 3 passes per active step (90th percentile 4, at
    most 6); the cap of 32 passes never binds. **The backlog goes:** the liquid held below the conjugate depth at the
    spray stage falls from a median 4.75 g to 1.23 g (at most 6.98 to 2.08 g; integrated over the flight from 295 to
    80 g s), and what waits at the end of a step for the next step's feed is a median 0.41 g (at most 0.81 g) — against
    3.5 g with the cascade's first, strict form. **The branch test hardly moves:** the film-weighted mean depth of the
    molten layer it reads falls only from 1.26 to 0.97 mm (it is 0.28 mm at 0.25 s and 0.066 mm at 0.125 s without the
    cascade), the deepest layer from 20.9 to 19.1 mm, and the thick-branch share of wet windward patches from 12.4 % to
    11.9 % (8.0 % and 2.6 % at the smaller steps). **So neither does the droplet population:** 18.1 million droplets
    against 17.2 million (28.2 and 58.5 million at 0.25 and 0.125 s), median radius by number 180 µm either way (105 and
    76 µm), by mass 190 against 189 µm (188 and 176 µm), the thin branch's share of the sprayed mass 10.5 % against 9.4
    % (18 % and 38 %), the front-surface Rayleigh–Taylor release 51.9 against 53.0 g (17.2 and 8.5 g), the re-solidified
    mass 6.2 g either way (16.8 and 48.4 g) and the largest droplet 4.7 against 4.4 mm. **The mass budget shifts a
    little:** sprayed 1.055 against 1.018 kg (+3.6 %) and the body at 120 s 0.417 against 0.454 kg (−8 %), with 1.4 %
    less heat absorbed (1.487 against 1.508 MJ with the deep runoff on) and the droplets leaving at 1.058 MJ/kg against
    1.065 — partly because the cascade feeds whole elements whose mean is above T_feed but whose coldest node is in the
    mushy range (fact 54), which the feed of step (i) would have held until it passed the ramp; whether that moves the
    budget toward or away from the small-step answer could not be measured (fact 60). The equator ring (windward patches
    beyond 0.9 of the transverse radius, 25.5–80 s) receives 3.5 % of the melt delivered (owner feed, cascade and deep
    runoff together) against 4.0 %. The peak surface temperature is 973 K either way; the hottest node of the frames
    reaches 1141 K from 80 s as before (1141–1148 K), and their 99.9th percentile rises from 1014–1017 K to 1027–1029 K,
    so the thin-owner hazard of fact 56 shows no flight-level excursion. Energy balance −6.6e-11 of the absorbed heat;
    3.00 Newton iterations per step. **Bottom line: the backlog fact 50 identified is real and the cascade removes it,
    but it was not what makes the droplet population depend on the step.** Fact 58 is.

58. **The next cause: the branch test's liquid depth is counted in whole elements and read at the end of the
    conduction.** What makes a patch "thick" (layer > δ_m) was measured on the wet windward patches under Girin's
    closure over 55–95 s. Without the cascade at 0.5 s, 2 412 of 4 704 such patches are thick, the molten depth under
    them a median 2.14 mm and the film plus deep liquid 0.38 mm; with it, 2 664 of 4 717 — no fewer — with 0.97 mm and
    0.40 mm; at 0.25 s without the cascade 2 182 of 4 953, 0.82 and 0.15 mm; at 0.125 s, 674 of 4 759, 0.62 and 0.06 mm,
    against a conjugate depth of 286–296 µm throughout. The reason the molten depth never falls below about 0.6 mm on a
    thick patch is the mesh: `molten_depth` adds whole elements, and on the production mesh a single molten wall-owning
    element already counts as 0.67–1.16 mm (5th to 95th percentile over the 18 830 patches, median 0.85 mm; 1.5 × the
    owner's centroid depth), above the conjugate depth (0.215–0.421 mm) on **every** patch. The branch test is therefore
    effectively binary on one question — is the wall-owning element's mean temperature at or above T_feed when the
    spray step runs? — and that is a splitting artefact: the conduction heats the surface element for a whole macro
    step before the feed's energy debit pulls it back onto the ramp, so with 0.5 s between feeds it overshoots T_feed
    far more often than with 0.125 s. The film at the spray stage adds a second step-proportional term, because it
    holds one step's melt supply (0.38–0.40 mm on thick patches at 0.5 s against 0.06 mm at 0.125 s). Neither is a
    liquid depth: the owner's liquid has already been fed into the film when the test reads its geometric depth (it is
    counted twice), and the elements are counted whole whatever their φ. A step-independent branch test needs a liquid
    depth by mass — the film, the deep account and the liquid inventory of the elements below the owner, as the deep
    runoff's h_D already is (sub-plan 09's amendment of 2026-10-02) — or a regime test on rates, melting speed against
    stripping speed, which is Girin's own statement of the regimes (fact 36). Both change fact 28(a)'s decision and are
    for Asha (fact 61). Sprayed mass and the mass budget are insensitive to all of this (the release is
    supply-limited, fact 9): what the step moves is the droplet size, the branch split and the Rayleigh–Taylor release.

59. **The deep runoff after the cascade.** With the cascade on, the deep runoff (at 0.5 s) mobilises 91 g from below the
    conjugate depth over the 120 s instead of 285 g, its account peaks at 11.9 g instead of 117 g (two-thirds of it, at
    the median, still deeper than its facet is wide, against 96 %), the liquid it sees at the spray stage is a median
    0.73 g instead of 1.82 g, and the cascade feeds 157 g (245 g with the deep runoff off): the two now share the liquid
    that waits below the wall-owning elements. Its effect on the results — deep runoff on against off, both with the
    cascade — is now small: sprayed mass +0.5 % (1.060 against 1.055 kg; it was +7.1 %), mass at 120 s −1.3 % (0.412
    against 0.417 kg; it was −16 %), droplet count −1.5 %, median radius unchanged, front-surface Rayleigh–Taylor
    release +8 % (56.3 against 51.9 g; it was +37 %), re-solidified mass +3.5 % (6.4 against 6.2 g; it was +80 %), and
    the equator ring's share of the melt 4.1 % against 3.5 % (3.5 g of deep liquid reaching it against 10.8 g). Against
    the run-to-run spread of fact 52 (0.15 % in sprayed mass, 0.33 % in mass, 14 % in re-solidified mass, about 3 % in
    the Rayleigh–Taylor release) the sprayed-mass and mass effects are three to four times the spread and the others
    within or near it. What it still moves is the per-step melt below the wall-owning elements — an artefact of the
    step, which fact 50 shows vanishing at 0.125 s — but at a tenth of the mass, so the reason fact 53 (a) gave for
    switching it off (it turned the backlog into piles) has largely gone.

60. **The 50 mm flight, and what was not measured.** Measured: the 50 mm flight, which never has Girin's closure, so the
    deep runoff is inert there and every patch is on the thin branch, but whose feed gate holds the same kind of backlog
    under its owners (fact 53 (e)): at the default step (seeded) the cascade feeds 95.6 g in 45 of 398 steps — a median
    of 5 passes, at most 10, never capped — the liquid held below the owners falls from a median 3.5 to 1.6 g and the
    end-of-step backlog is 0.18 g. The flight changes more than the 100 mm one: demise 8.0 s earlier and 2.5 km higher
    (199.0 s and 70.7 km against 207.0 s and 68.2 km), sprayed mass +1.0 % (0.1770 against 0.1753 kg), droplets −32 %
    (1.92 against 2.84 million), median radius by number +17 % (173 against 149 µm) and by mass +11 % (227 against
    206 µm), largest droplet 1.81 against 1.47 mm, and the droplets leave with 4.4 % less enthalpy per kilogram (1.089
    against 1.140 MJ/kg); the energy balance is exact (7.4e-11). The larger droplets are consistent with the thin
    branch's cap on the film present on the patch (fact 9), which grows when the cascade's liquid arrives in step-sized
    lumps (not isolated). Whether these moves are toward the converged answer was not measured: the 50 mm flight has not
    been run at a smaller step with the cascade.
    **Not measured, because the machine ran on battery.** From the early morning to past midday on 2026-10-03 this Mac
    ran on battery, at 7 % falling to 4 %, sleeping most of the time (forty minutes of wall time gave under a minute of
    CPU), so the several CPU-hours the rest of the study needs could not be run. Not measured: the 0.25 and 0.125 s
    flights with the cascade (and so the convergence of the cascade's own results — fact 57 compares the cascade at
    0.5 s with the model *without* it at the smaller steps), the whole 100 mm flight to the ground with the cascade, and
    the run time on mains power (the cascade runs' 1 149 and 1 182 s against 946 and 975 s mix in battery throttling;
    the cascade's first form, run on mains power beside a flag-off run, cost nothing measurable, 917–931 s against
    921 s, and each extra pass costs one surface rebuild, 0.085 s on the production mesh, so about three passes in 146
    steps should add about 40 s, 4 %). Expected, not shown: at 0.25 and 0.125 s the liquid held below the conjugate
    depth is already 0.23 and 0.004 g without the cascade, so it has little to feed there and cannot close the gap of
    fact 57 on its own; fact 58's mechanism remains at every step. Sub-plan 14's amendment lists the runs to make.

61. **For Asha to decide.** (a) **The branch test's liquid depth (fact 58) — decide this first; it is now what keeps the
    droplet population from converging in the step.** Options: (1) keep `molten_depth` (whole elements, read at the end
    of the conduction) and run melting flights at a small step — 0.125 s still does not converge it, and costs four
    times the default; (2) measure the layer by mass: the film, the deep account, and the liquid inventory min(f φ, φ −
    PHI_MIN) ρV of the contiguous chain *below* the wall-owning element, divided by ρ_l A — the owner's own liquid is
    already in the film — which removes the element quantisation and the double count; what remains step-dependent is
    the film's one-step supply; (3) decide the regime on rates rather than on a depth — Girin's own statement is that
    regime 1 is the case where the mass loss overtakes fusion (fact 36), and with release capacity hundreds of times the
    melt supply (fact 9) this flight is in it almost everywhere — which would make the thin branch the rule here and the
    thick branch the exception. Recommendation: (2), as an amendment of its own with the time-step study repeated,
    because it is the smallest change that makes the test measure liquid, and (3) as the check on whether the result is
    physical rather than a property of the film account. Until then, no droplet-population number — count, branch split,
    median radius by number, Rayleigh–Taylor release, re-solidified mass — should be quoted from the default step.
    (b) **Keep the molten cascade as the default?** It removes the backlog at its cause (fact 57) and most of what the
    deep runoff was amplifying (fact 59), costs about 4 % of run time, keeps the books exact and agrees across the two
    backends, and changes the droplet population by under 5 %. Its one declared approximation is that it feeds whole
    elements whose mean is above T_feed though their coldest node is in the mushy range (fact 54); it raises the
    sprayed mass by 3.6 % and lowers the mass at 120 s by 8 % at the default step, and whether that is toward the
    small-step answer is the first thing the missing runs (fact 60) should show. Recommendation: keep it on, re-measure
    at 0.25 and 0.125 s, and switch it off only if the cascade's own results do not converge in the step.
    (c) **The deep runoff's default** (fact 59): it now changes the sprayed mass by +0.5 % and the mass at 120 s by
    −1.3 % at the default step, three to four times the run-to-run spread, its piles hold a tenth of the mass they did,
    and at smaller steps the liquid it moves vanishes (fact 50). Recommendation: keep it on, as built — the reason fact
    53 (a) gave for switching it off (the backlog turned into piles in the craters) has largely gone — and quote any
    number that depends on it with the step it was measured at.
    (d) **`PHI_DEATH`** (fact 44's 0.50, not in this copy): with 0.50 a wall-owning element that the feed has taken half
    of dies in the same step, so partly molten owners stop surviving several steps and the cascade would act behind
    them sooner; it also removes most of the thin owners that fact 56's hazard needs. Measure the time-step study again
    with it before choosing between (a)'s options.
    (e) **The cascade and the deep runoff share the liquid below the owners** (fact 59): the deep stage runs first and
    takes what flows; the cascade feeds what is left when the surface reaches it. If deep liquid should only ever be
    what the surface cannot reach within the step, the cascade could run before the deep stage instead — a reordering
    with consequences of its own for the spray step, not measured.
    (f) **Seeding** (fact 52): still recommended as one line in `cli.cmd_run`; every run here was seeded in the harness.

## Amendment of 2026-10-05 — numpy's generator is seeded at the start of every run (facts 62–68)

Asha's decision of 2026-10-05 on facts 53 (i) and 61 (f): a `--seed` option with a fixed default, applied once at the
start of `cli.cmd_run`, before anything that could draw random numbers, and recorded in each run's summary. It amends
sub-plans 13, 14 and 15; sub-plan 13's amendment holds the design and the tested code. No other sub-plan is affected,
and that was checked against the code rather than assumed: the seed is set, checked and recorded in `cli.py` alone — the
run JSON's settings are assembled in `cmd_run`, not in `coupled.py`, so sub-plan 10's history columns and results are
unchanged; no `MeltSettings` field or melt-step behaviour changes, so sub-plan 09 is untouched; and sub-plan 03's
thermal core, where pyamg is called, is left as it is (fact 63 says why). Facts 1–61 stand except where these say
otherwise. **Fact 31's "the model is reproducible" and fact 52's "the model was not bit-reproducible" are both
superseded by fact 64; fact 52's account of what draws is made exact by fact 62, and its spread (0.15 %, 0.33 %, 1.2 %
and 14 %) is replaced as the floor for comparisons by fact 65. Facts 53 (i) and 61 (f) are answered.** Measured in a
throwaway copy of `prototype/proto3/` with the amendments of 2026-10-02 and 2026-10-03 applied — verified before any
change to be reproduced byte for byte by their 17 diff blocks, applied in date order to a fresh copy of
`prototype/proto3/`, itself unchanged since 2026-10-03 (all 114 files match the manifest taken then) — so everything
facts 49 and 54 list about the copy holds: `AA7075_range`, US76, physics heating, the dense band, `PHI_DEATH = 0.05`,
the deep runoff and the molten cascade on, and every other setting at its default. Unlike the runs of facts 46–61, no
run here was seeded by a measurement harness: each is the copy's own command line, which now seeds itself. Every run was
made on mains power under `caffeinate -i` (the power source logged at the start and end of each).

62. **What draws from numpy's generator, measured.** A probe recorded the global generator's state before and after the
    package import, at the entry of `cmd_run` and after the run, attributed every call of numpy's module-level random
    functions to its caller, and replayed the calls on a fresh generator from the entry state. On a 15 s melting run
    (a warm start at 71 km on the coarse mesh) the import and the argument parsing leave the state untouched; the run
    makes 56 draws, every one `np.random.rand` inside pyamg's `approximate_spectral_radius`; and replaying them
    reproduces the run's final state exactly, so nothing draws by any other route. The mechanism, read in pyamg 5.3:
    `smoothed_aggregation_solver`, which the skfem backend builds at its first solve, every 30 linear solves
    (`amg_rebuild_every`) and whenever a CG solve fails and is retried, smooths the tentative prolongator of every level
    but the coarsest by one Jacobi step weighted by 4/3 over the spectral radius of D⁻¹A, and
    `approximate_spectral_radius` estimates that radius by an Arnoldi iteration (15 iterations, 5 restarts, a 1 %
    tolerance) from a random starting vector of the level's size — one draw per level per build, measured. The weight,
    and with it the preconditioner, therefore changes from one build to the next at the level of a 1 % estimate; the
    conjugate gradients still converge to 1e-10, but along another path, so the solution differs in the last bits.
    Fact 52's "smoother weights" are this prolongation smoother's; the pre- and post-smoothers are symmetric
    Gauss–Seidel sweeps and carry no weight. How the difference grows: in a 5 s Step 2 run on the coarse mesh two
    unseeded runs differ in 6 of 31 history columns, by 2e-16 to 6e-16 relative; in a 30 s warm-start melting run on the
    coarse mesh seeds 12345 and 1 differ from the first macro step, by at most 1.6e-9 K in the mean temperature and
    3.7e-10 in the sprayed mass at any time, ending 1.9e-13 apart in sprayed mass and 3.4e-10 in droplet count; on the
    100 mm production flight to 120 s the four seeds of fact 65 part at the first macro step, in the last bits of the
    surface temperatures, and are still within 1e-12 of each other in sprayed mass at 30 s, 9e-7 at 40 s and 4e-7 at
    60 s; the range then opens with the melting, to 5e-4 at 80 s and a peak of 2e-3 at 90 s, and ends at 1.3e-3. It is
    the melting model's thresholds — the feed ramp, element death, the thick/thin branch test — that turn a last-bit
    difference into the spread of fact 65. The FEniCSx backend never calls pyamg (fact 66).

63. **The design** (sub-plan 13's amendment has the code). `--seed` in the main `run` group, an integer in [0, 2³² − 1]
    (`parse_seed`; anything else is argparse's exit 2), default `DEFAULT_SEED = 12345`, applied as the first statement
    of `cmd_run` — before the argument checks and before anything is built, so that no draw, present or future, can come
    before it; on the amended copy fact 62's probe records `np.random.seed(12345)` as the first call after `cmd_run` is
    entered, and then the same 56 draws. The run JSON's `settings` carry `seed` for every run; the printed summary's
    first line names it (`<run name> (seed 12345): ...`); and a seed other than the default ends the run name in
    `_seed-<n>`, in every mode, because a run's name encodes its whole configuration and the seed is applied to every
    run (the default leaves every existing name unchanged). 12345 because it is the harness seed of facts 46–61, so the
    model reproduces those runs by itself (fact 64); fact 52's `np.random.seed(0)` would have left them out of reach.
    Alternatives not taken: a fixed starting vector inside the skfem backend (reproducible wherever the backend is
    built, the tests included, but a change to the thermal core — `smoothed_aggregation_solver` exposes no argument for
    that estimate's starting vector — that would move every run's numbers away from the measured ones); `--seed random`,
    a drawn and recorded seed (a scatter study is reproducible only if its seeds are listed, which `--seed n` already
    allows); no default seed (an unseeded run cannot be repeated and offers nothing a run with another seed does not);
    and the suffix only on `--thermal fem` runs, the only ones that draw (the seed is applied to all). **Where else runs
    are built.** The verification drivers call `cli.main` in-process and the sensitivity driver runs the CLI as
    subprocesses, so every run they make is seeded with the default, whatever ran before it in the same process; the
    sensitivity table gains a `seed1` row (`--seed 1`), the same flight with another seed, so that it carries its own
    floor (sub-plan 14). The reference-tier tests build their runs directly and stay unseeded, by decision (fact 68
    (a)).

64. **Every run now repeats bit for bit, and the default reproduces the harness's runs.** Two default-seed runs of the
    whole 50 mm flight (398 macro steps, demise at 199.0 s) are identical in all 399 rows of all 86 history columns at
    full precision, in the final state (nodal temperatures, element fractions, active set, film and deep accounts), in
    all 21 010 source-table rows of all 22 columns, in every result field but the run time and in every array of the 40
    VTK frame files; two default-seed runs of the 100 mm flight to 120 s likewise (241 rows, 157 663 source-table rows,
    26 frame files). And the model's own default reproduces the harness-seeded runs of 2026-10-03 — the 50 mm flight of
    fact 60 and the 100 mm flight to 120 s with the cascade and the deep runoff on of fact 59 — in every history cell,
    every source-table row, every result field but the run time and every array of every frame (the harness kept the
    history only to the CSV's 9 digits, but the source table and the frames are full precision); and with
    `--deep-runoff off --molten-cascade off` it reproduces the harness-seeded run of 2026-10-02 of the unamended
    prototype — the 100 mm flight to 120 s that facts 49 and 50 start from, made in a copy identical to
    `prototype/proto3/` — in all 241 rows of the 79 history columns that model had, all 150 352 source-table rows, all
    48 result fields it had and every array of every frame, the two surface fields it did not yet write apart. The
    harnesses seeded 12345 before importing the package and the model seeds at the start of `cmd_run`; nothing draws
    between the two, and neither harness's instrumentation changed anything. The seeded runs of facts 49–61 are
    therefore reproducible from the command line alone, with the flags that select the model each was made with —
    verified on these three, which span both harnesses and both the oldest and the newest model. This supersedes fact
    31's "the model is reproducible", which was drawn from runs that happened to agree, and fact 52's "not
    bit-reproducible"; the advice of fact 52 and of sub-plan 14's amendments to compare runs as seeded pairs is now met
    by every run.

65. **The scatter between seeds is the floor for every comparison.** Measured on the 100 mm physics flight to 120 s at
    the default step, every setting at its default (the deep runoff and the cascade on), with the seeds 12345, 1, 2 and
    3: four runs that are equally valid and differ only in the round-off of fact 62. The range across the four, as a
    share of their mean: sprayed mass 1.0600 to 1.0614 kg, **0.13 %**; the mass at 120 s 0.4107 to 0.4120 kg,
    **0.33 %**; the droplet count 1.762e7 to 1.798e7, **2.0 %**; the median radius by number 179.6 to 180.1 µm,
    **0.29 %** (by mass 190.0 to 190.3 µm, 0.16 %); the re-solidified mass 6.38 to 6.46 g, **1.4 %**. Further: the
    front-surface Rayleigh–Taylor release 56.2 to 60.0 g, 6.7 %; the thin branch's share of the sprayed mass 10.2 % to
    10.3 %, a relative 1.6 %; the median depth of the molten layer 0.94 to 0.96 mm, 2.7 %; the cascade's mass 0.4 % and
    the deep runoff's 0.7 %; the heat absorbed 0.04 %, the enthalpy carried away per kilogram 0.01 % and the Newton
    iterations per step 0.4 %. The melt onset (73.96 km), the largest droplet (5.058 mm, within 0.003 %) and the median
    thick-branch share of wet windward patches (11.6 %) do not move. Against fact 52's spread, measured on the model
    before the cascade from two unseeded runs (0.15 % in sprayed mass, 0.33 % in the mass at 120 s, 1.2 % in droplet
    count and 14 % in re-solidified mass), the masses are where they were, the droplet count spreads somewhat more and
    the re-solidified mass ten times less. **What it means for the facts already written:** fact 57's cascade effect on
    the masses (sprayed +3.6 %, mass at 120 s −8 %) is about 25 times this range and its droplet count (+5 %) about 2.5
    times; fact 59's deep-runoff effect with the cascade on the masses (+0.5 % and −1.3 %) is about four times the
    range, as fact 59 said, while its droplet count (−1.5 %) lies inside it, its Rayleigh–Taylor release (+8 %) only
    just outside and its re-solidified mass (+3.5 %) 2.5 times; fact 50's time-step trends (1.7e7 to 5.9e7 droplets, 180
    to 76 µm) are far outside it. The 50 mm whole flight, measured with one pair of seeds (12345 and 1), scatters far
    less: no history column moves by more than 2.5e-6 of its value, the sprayed and final masses by 2e-13, the droplet
    count by 6e-10 and the re-solidified mass by 1.1e-7, with demise at 199.0 s in both — so on this model its late
    collapse did not amplify the round-off as fact 52's 50 mm pair (0.85 % in final mass, on the earlier model) did,
    though one pair is one sample. **How to use it:** quote a difference between two settings against this range, on the
    flight and at the step it was measured at; a smaller difference is not a result. Four seeds give a range, not a
    distribution, and the range of four draws understates the full spread (fact 68 (b)).

66. **The FEniCSx backend is reproducible on its own and does not see the seed.** It solves with PETSc's CG and hypre's
    BoomerAMG and never calls pyamg; the probe of fact 62 finds no draw at all in a FEniCSx run, and the generator's
    state is the same after the run as before it. Measured on the 30 s warm-start melting run on the coarse mesh in
    `fenicsx_env` (60 macro steps, 21 202 source-table rows): two runs with the default seed are identical in every
    history column at full precision, the final state, every source-table row and every result field but the run time,
    and a run with seed 1 is identical to them as well; the skfem backend on the same configuration differs between
    those two seeds (in 27 cells of the history CSV at its nine digits, 10 of the 22 source-table columns and 23 result
    fields, by at most 3.4e-10 of the droplet count and 1.9e-13 of the sprayed mass at the end). Measured in serial
    only: a parallel FEniCSx run partitions the mesh, and whether the partition and hypre's coarsening repeat from run
    to run was not measured. Consequence: the backend comparisons of facts 49 and 56 seeded both backends, but only the
    skfem side ever depended on it, and the `fenicsx` sensitivity row's difference from `base` contains `base`'s own
    scatter (sub-plan 14).

67. **The seed costs nothing.** `np.random.seed(12345)` takes 3.5 µs, once per run, against run times of minutes.
    Measured on the whole 50 mm flight on an otherwise idle machine, alternating the copy without the amendment (no
    seed) and the copy with it: 181.2 s and 181.2 s without the seed, 182.0 s and 180.9 s with it — a mean difference of
    0.25 s (0.14 %), smaller than the 1.1 s between the two seeded runs themselves; 398 macro steps and 2.236 Newton
    iterations per step in all four. The two unseeded runs differed from each other, as fact 52 says they may (by at
    most 1.1e-7 of any history value, 3e-8 in re-solidified mass and 5e-14 in sprayed mass on this flight); the two
    seeded ones did not, and are identical to the default 50 mm run of fact 64, made half an hour earlier beside three
    other runs, so the result does not depend on the machine's load either. The run times of the scatter runs of fact 65
    (597 to 761 s) say nothing about the seed: they ran two to four at a time.

68. **Not done, and for Asha to decide.** (a) **The reference-tier tests.** `tests/test_reentry_model_reference_melt.py`
    and Step 2's `tests/test_reentry_model_reference_thermal.py` build their runs directly rather than through
    `cmd_run`, so they are not seeded and the metrics files they write can differ from one run to the next in the last
    digits. Measured on the bookkeeping device through the CLI (the melting test's configuration on the same
    177 363-element mesh, which the prototype cannot run as a test until Task 11 commits the melting references): seed
    12345 against seed 1 moves no history column by more than 2.4e-9 of its value, the mass by at most 1.4e-11 of the
    initial mass, and the melt onset (71.10 km at 43.0 s) and the demise (66.0 s) not at all — nine to ten orders of
    magnitude below the thresholds of 2 % of the mass, 0.5 km and 2 % of the 1 %-mass time; a 5 s Step 2 run moves by
    2e-16 to 6e-16 (fact 62). Options: (1) leave them unseeded; (2) add `np.random.seed(cli.DEFAULT_SEED)` at the top of
    each run helper, one line per file, so that the files they write repeat to the last bit. Recommendation: (2), in the
    same change in which Task 11 commits the melting references, so that it is tested when it is made; nothing the
    README prints depends on it before then. (b) **How many seeds the floor rests on.** Fact 65's scatter is the range
    of four seeds on one flight at one step, and the 50 mm flight's is one pair. Options: (1) take it as the floor for
    the comparisons already written (facts 57–60) and stop there; (2) measure more seeds — a 100 mm run to 120 s costs
    10 to 13 minutes of one core beside others — to estimate a standard deviation rather than a range; (3) run each
    comparison the thesis quotes at two or three seeds and quote its mean and range, so that every result carries its
    own scatter on its own flight and at its own step. Recommendation: (1) now and (3) for the thesis's numbers, because
    the scatter belongs to the model version, the flight and the step — between the model before the cascade and this
    one the re-solidified mass went from 14 % to 1.4 % — and the sensitivity table's `seed1` row is the cheapest
    standing check. Not done, and recorded so it is not lost: (c) a drawn seed (`--seed random`), which a scatter study
    does not need (fact 63) and which would still have to be recorded in the JSON and the name; (d) parallel FEniCSx
    runs, whose reproducibility was not measured (fact 66).

## Amendment of 2026-10-05 (runoff flux) — the film's flux on a thick patch, and the switched-step study (facts 69–77)

Asha's two decisions of 2026-10-05, taken on the switched-step experiment of fact 69 and the crash of fact 70: (1) fix
the runoff flux on thick patches by option (a), scaling the film's flux on a thick patch by the film actually present,
and verify it by re-running the case that crashed; (2) keep shrinking the time step to make the droplet population
converge, rather than reformulating the spray release, running the switched-step series with the fixed model down to
0.00625 s and saying plainly whether it converges. It amends sub-plans 06 (the change to `film.lubrication` and its
tests), 07 (no code change), 09 (no change to `body.py`; one melting test's tolerance and a two-backend test), 13 (no
code change; an open item), 14 and 15. Facts 1–68 stand except where these say otherwise. **The 2026-10-02 amendment's
statement that the film's flux plus `deep_flux` is fact 29's column flux now holds only where the film fills the
conjugate layer (fact 71); facts 49 and 56's reason for the backends' looser agreement with deep liquid present is
replaced by fact 72's; fact 48's film-transport round-off is re-measured in fact 72; and facts 50 and 57's time-step
findings are extended by facts 69 and 75.** Measured in a throwaway copy of `prototype/proto3/` with the amendments of
2026-10-02, 2026-10-03 and 2026-10-05 (seeding) applied — verified before any change to be reproduced byte for byte by
their 20 diff blocks, applied in date order to a fresh copy of `prototype/proto3/`, itself unchanged since 2026-10-03
(all 114 files match the manifest taken then) — so everything facts 49, 54 and 62 list about the copy holds:
`AA7075_range`, US76, physics heating, the dense band, `PHI_DEATH = 0.05`, the deep runoff and the molten cascade on,
every run seeded by the model itself with the default 12345, and every other setting at its default. The measurement
runs were made from a frozen copy of the amended package, never edited while a run was in flight (fact 31), on mains
power under `caffeinate -i`, the power source logged at the start and end of each.

69. **The switched-step experiment, and what it showed before the fix (2026-10-05).** Two questions: does the melt that
    forms on the windward face run off before it is sprayed, and does the droplet population converge if the macro step
    is shortened where it matters? A harness that wraps the model from outside (`dtswitch.py`; no package file is
    edited, as with the cascade study's measurement harness) runs the command line's step, 0.5 s, until the body Knudsen
    number first falls below 0.01 — the model's own continuum boundary `surface_flow.KN_BODY_SHOCK`, where Girin's
    closure begins to apply — and a fine step from then on, latched; on the 100 mm flight it fires at 49.5 s and
    69.9 km. It writes a frame every 10 s of flight time and records, per melt step, the film available to the spray
    stage, the runoff's arrivals and departures and the release, binned by the polar angle about the flight axis through
    the mass centre, plus the equatorial ring's arrivals and releases and the deep transport's arrivals. Its summary
    (`compare.py`) uses only measures that do not depend on the step: the net runoff across fixed latitude lines (15,
    30, 45, 60, 75, 85 and 90 degrees), the mass-weighted shift between the latitude where film entered and where it was
    sprayed, the ring's share of the sprayed mass and, from the source table, the droplet count, the median radius by
    number and the thick branch's share of the sprayed mass. Per-step net arrivals are not comparable between steps,
    because film that creeps across N patches in N small steps counts N times. Measured on the 100 mm flight to 120 s at
    0.5 s throughout and switched to 0.05 and 0.025 s, all over the continuum window 49.5–120 s unless stated:
    * **Runoff is minor at every step.** At most 6.3 %, 3.1 % and 5.8 % of the film formed inside any latitude cap
      crosses its edge; the mass-weighted shift between where film entered and where it was sprayed is 1.14, 0.55 and
      0.72 degrees of latitude; the equatorial ring sprays 9.5 %, 7.8 % and 8.5 % of the mass.
    * **The droplet population does not converge.** 15.4, 87.2 and 110.9 million droplets; median radius by number 182,
      73 and 71 µm, by mass 191, 158 and 126 µm; the thick branch's share of the sprayed mass 92 %, 47 % and 34 %; the
      mass-weighted mean film depth at release 3.30, 0.34 and 0.26 mm. At 0.025 s the thick share rises through the
      flight, from 29 % over 49.5–60 s to 53 % over 100–120 s.
    * **Two effects vanish at small steps:** the front-surface Rayleigh–Taylor release, 56, 5.2 and 3.3 g, and the deep
      runoff, 91, 0.9 and 0.8 g.
    * **The masses move by 1–3 %:** sprayed 1.060, 1.050 and 1.047 kg, the body at 120 s 0.412, 0.422 and 0.425 kg.
    * **The re-solidification counter is not a result.** It reads 6.4, 146 and 287 g and grows by 0.10 g per fine step
      at both small steps: it counts the same film freezing and re-melting from step to step, not lasting refreezing.
    * **Cost.** The 0.05 and 0.025 s runs took 2 569 and 4 557 s of wall time with three or four runs in parallel; a
      fine step costs about 1.7 s alone. The 87 g sprayed before the switch (49.5 s) is the thin-film mode under the
      Couette closure (94.9 % of the 92.1 g released by then), identical in every run because every run computes that
      part of the flight at 0.5 s; whether it depends on the step there was not measured.

70. **The crash at 0.0125 s, and its cause.** The run switched to 0.0125 s failed at 61.5875 s of flight time, after 967
    fine steps, with this chain in its log: an overflow warning in `film.lubrication`'s shear rate; "Matrix is exactly
    singular" from `spsolve` in `film.Runoff.transport`; a NaN film mass on all 16 966 patches; and the thermal solver's
    input check, "film mass must be one finite non-negative value per node", which `cli.main` reports as exit code 2
    (fact 77). A re-run with a probe that re-implements the transport operation for operation and records every
    sub-step's largest edge emptying rate and smallest positive film depth (`dtswitch_diag.py`) crashed at the same step
    — the run is deterministic — and confirmed the cause. `Runoff.edge_coefficients` uses the emptying rate
    c = q ℓ (t̂·n̂)⁺/(A b). On a thin patch the flux q = τ b²/(2 μ_l) + G b³/(3 μ_l) vanishes with the film b, so c
    stays bounded. On a thick patch — where the liquid layer, film plus molten depth plus deep liquid, is deeper than
    δ_m — `lubrication` gave q = V_s δ_m/2 + G b³/(3 μ_l) with V_s = τ δ_m/μ_l: the flux of a whole conjugate layer,
    which does not vanish with the film, so c grows like 1/b. How often: 95 % of the film transport's sub-steps at
    0.0125 s had a largest rate above 10⁶ per second, 66 % above 10²⁰ and 8.6 % above 10¹⁰⁰; the largest was 4.4·10¹⁹⁷
    per second. The runaway: the implicit sub-step empties such a film by a factor of about c Δt_s, which raises the
    next sub-step's rate by the same factor, so the film's depth roughly squares from one sub-step to the next — in one
    call 5.3·10⁻²⁴, 3.8·10⁻⁴³, 2.1·10⁻⁸² and 5.5·10⁻¹²² m with rates 3.2·10²⁰, 1.2·10⁴², 2.3·10⁸¹ and 8.4·10¹²⁰ per
    second, the rate times the depth staying at 0.46 m/s as the formula says it must. In the failing call, sub-step 2
    (rate 1.7·10⁴³ per second, film 2.3·10⁻⁴⁴ m) returned non-finite values, which the solve spread to every patch. A
    unit test reproduces the failure exactly (fact 72): two such patches that drain into each other form a two-by-two
    block whose second pivot is 1 + c Δt − c Δt, which is 0 in floating point once c Δt exceeds about 10¹⁶; which pair
    failed in the flight was not identified, because the probe kept the state of the last non-finite sub-step, by which
    time every patch was NaN. The 0.05 and 0.025 s runs logged the same overflow warning but finished. The same formula
    was a physics error as well: it moved a film thinner than δ_m on a thick patch as fast as a whole conjugate layer,
    overstating the runoff there. The deep transport was not involved: its largest rate in the probe was 5.3·10⁵ per
    second.

71. **The fix, option (a): the film carries only its own part of the conjugate layer's flux** (sub-plan 06's amendment
    has the code and the alternatives). The thick branch's shear-driven velocity is linear within the conjugate layer,
    V_s at the free surface and zero at depth δ_m, and the film is the top b of that layer, so for b < δ_m its
    shear-driven flux is the integral of that profile over its own depth, V_s b (1 − b/(2 δ_m)); for b ≥ δ_m it stays
    V_s δ_m/2, computed by the same expression as before. The flux is continuous at b = δ_m, tends to V_s b as the film
    vanishes, so that the emptying rate's shear part is at most V_s ℓ/A, and is monotone in b; the pressure and
    deceleration part G b³/(3 μ_l) is unchanged, and so are V_s (which the spray's Weber number reads), the shear rate,
    the branch flag and the whole thin branch. `lubrication` has two callers, both in `MeltingBody._film_and_spray`: the
    film transport, whose flux changes, and the call after it, which uses only V_s and the branch flag. `deep_flux`
    needs no change: its emptying rate per unit edge length, G (b² + b h_D + h_D²/3)/(μ_l A), stays bounded as the deep
    liquid h_D vanishes. One statement of the 2026-10-02 amendment changes: the film's flux plus `deep_flux` is fact
    29's column flux τ d²/(2 μ_l) + G h³/(3 μ_l) where the film fills the conjugate layer, and falls short by
    τ (δ_m − b)²/(2 μ_l) where it does not — the shear flux of the part of the conjugate layer that lies in the
    elements, mostly the wall-owning element that `molten_depth` counts whole after its liquid has been fed into the
    film (fact 58), which no account moves. Where deep liquid lies under a patch the identity holds at the start of
    every film transport, because the deep stage tops the film up to δ_m first. Alternatives rejected: a floor on the
    film depth and a cap on the emptying rate (both treat the symptom, add a constant with no physical meaning and leave
    a micron of film moving like a whole conjugate layer); the layer's flux scaled by the share of it present, V_s b/2
    (moves the skin at half the surface speed the spray credits it with); the thin branch's Couette flux τ b²/(2 μ_l)
    (contradicts the thick branch's own V_s; at b = δ_m/10 nineteen times slower than the profile); moving the lower
    part of the conjugate layer too (fact 53 (f)'s open decision, not a fix); and a switch to keep the old flux (it is a
    defect, not a modelling alternative).

72. **Tests, the two backends, and the film transport's round-off.** Two new film tests and one changed one (sub-plan
    06): a thick patch's flux vanishes with its film and stays below V_s b, matches the profile's integral by quadrature
    to 1e-9, is continuous at δ_m and is bit-identical at b ≥ δ_m and on thin patches; the transport stays finite,
    non-negative and conservative (to 2.5e-14 over eight steps of 0.0125 s) on a pair of nearly dry thick patches that
    drain into each other, for films down to 5e-324 kg, the smallest positive double, with its largest rate 1.15e4 per
    second — where the old flux gave "Matrix is exactly singular" and NaN on every patch for six of seven film masses
    between 1e-25 and 1e-300 kg; and the deep liquid's rate stays bounded as it vanishes. The 2026-10-02 test of the
    column flux now asserts the shortfall of fact 71. All three fail on the unamended module. **One melting test's
    tolerance moved** (sub-plan 09): the cascade's six-step books test held the mass to 1e-12 of the body and failed by
    1.6e-12. The drift is the film transport's direct solve, which conserves only to its conditioning (fact 48); over
    ten states of numpy's generator its largest six-step drift was 6.8e-13 to 2.8e-12 of the body before the fix (one
    state of ten already failing) and 9e-14 to 4.4e-12 after it (six of ten), the worst single transport moving the
    total by 4.0e-12 and 6.4e-12 kg, 1.7e-11 and 2.6e-11 of the mass it moved; the pool drives up to a quarter of a
    kilogram of film per step through pairs of patches draining into each other with c Δt_s near 10⁷. The check passed
    in the seeded unit tier by chance and now holds the mass to 1e-11; the deep-runoff test's 3 mm pool, which keeps its
    1e-12, drifts by at most 4.8e-14 before and after. **The two backends agree better with the fix than without it.** A
    new FEniCSx test (the 6 mm pool at 69.8 km, ten steps of 0.05 s, the change exercised in every film transport) holds
    them to the same 7 983 active elements, 5.8e-10 in mass, 5.9e-9 in sprayed mass, 7.5e-10 in runoff, 1.7e-8 in film,
    2.8e-10 in the deep account and 2.3e-4 K, with both balances exact (5.1e-9 and 1.1e-11), and the largest emptying
    rate at 2.7e7 per second, where it reached 3.9e25 before the fix. In the cascade test's own setting (six steps of
    0.5 s) the passes are the same with and without the fix (4, 3, 16, 4, 5, 4), and the agreement improves from 1.3e-7
    to 7.6e-10 in mass, 6.3e-7 to 3.6e-9 in sprayed mass, 2.4e-7 to 1.0e-10 in cascade mass, 3.7e-6 to 6.9e-11 in the
    deep account and 0.11 K to 2.5e-5 K; the deep transport is unchanged, so the amplifier that facts 49 and 56
    attributed to it was the film's runaway emptying rates. Unit tier: 236 passed, 1 skipped, and the five known
    failures and errors from the missing melting SESAM references (Task 11) — the base's 234 passed and the two new film
    tests; `tests/test_reentry_model_fenicsx.py` in `fenicsx_env`: 11 passed, the ten of before and the new one.

73. **The fix at the default step: the 50 mm flight bit for bit, the 100 mm flight by about the scatter in its masses.**
    The **50 mm whole flight**, which never has Girin's closure, is bit-identical: all 399 rows of all 86 history
    columns at full precision, the final state (nodal temperatures, element fractions, active set, film and deep
    accounts), all 21 010 source-table rows of all 22 columns and every result field but the run time. The **100 mm
    flight to 120 s** is bit-identical until its first step under Girin's closure, at 49.5 s, and then moves; against
    the seeding amendment's default run, and read against the scatter between four seeds of fact 65: sprayed mass 1.0600
    to 1.0585 kg (−0.14 %, the scatter 0.13 %); the body at 120 s 0.4120 to 0.4136 kg (+0.37 %, scatter 0.33 %);
    droplets 1.779e7 to 1.789e7 (+0.6 %, scatter 2.0 %); median radius by number 179.9 to 179.7 µm (−0.1 %, scatter
    0.29 %) and by mass 190.3 to 190.6 µm (+0.2 %, scatter 0.16 %); re-solidified 6.39 to 6.62 g (+3.7 %, scatter
    1.4 %); the front-surface Rayleigh–Taylor release 56.3 to 47.4 g (−16 %, scatter 6.7 %), spread over the whole cap
    inside 15 degrees; the deep runoff 91.1 to 90.9 g (−0.2 %) and the cascade 157.1 to 158.0 g (+0.6 %). What the fix
    changes is which film reaches the spray and where: the source table has 25 % more rows (197 075 against 157 663, a
    median of 1 287 releasing patches per step against 927), the thick branch's share of the sprayed mass rises from
    84.4 % to 85.1 % over the flight (92.0 % to 92.8 % over the continuum window), the mass-weighted mean film depth at
    release falls from 3.30 to 2.85 mm, and the share of wet windward patches on the thick branch (the median over the
    continuum steps) rises from 21.1 % to 27.5 %: thick patches now keep thin films that the old flux swept off within a
    sub-step, so more of them are wet when the spray runs. No film piles up instead: the thickest film where a depth
    means anything is a median 1.97 mm over the continuum steps (2.23 mm before) and at most 3.9 mm (5.6 mm), and the
    share of the film deeper than its facet is wide a median 0.112 (0.121). The runoff measures barely move: at most
    6.1 % of the film formed inside a latitude cap crosses its edge (6.3 % before), the shift between where film entered
    and where it was sprayed is 1.19 degrees (1.14), the ring sprays 9.3 % of the mass (9.5 %). The energy balance stays
    exact (−8.5e-11 of the absorbed heat, against −7.5e-11), Newton takes 3.02 iterations per step (3.03), and the peak
    surface temperature is 981 K (974 K). The run time cannot be compared from these runs: they ran eight at a time
    (fact 76).

74. **The fine steps now run, with bounded emptying rates.** The 0.0125 s flight that crashed at 61.6 s runs to 120 s
    with the fix, and so does the 0.00625 s flight; neither log has an overflow warning or a singular matrix, and the
    energy balance stays at round-off (8.4e-10 and 2.3e-9 of the absorbed heat; it accumulates with the number of steps,
    −8.5e-11 at the default step). Both ran under the same probe as the crash (`dtswitch_diag.py`, unchanged), which
    records every transport sub-step. At 0.0125 s the film transport's largest edge rate has a median of 9.7·10³ per
    second over 22 756 sub-steps, a 99th percentile of 3.1·10⁵ and a maximum of 1.5·10⁶, and 14 sub-steps (0.06 %)
    exceed 10⁶ per second — before the fix a median of 1.0·10²⁶, a maximum of 4.4·10¹⁹⁷ and 95 % above 10⁶. Within one
    call the last sub-step's rate is a median 1.01 times the first's and at most 147 times, against a median of 3.7·10⁴⁰
    before: the runaway is gone. At 0.00625 s, over 45 320 sub-steps, the median is 9.8·10³ per second, the maximum
    1.45·10⁶ (two sub-steps above 10⁶) and the within-call growth a median 1.01, at most 83. Films of 10⁻¹⁹¹ m
    (0.0125 s) and 10⁻¹²⁶ m (0.00625 s) still appear inside the transport — round-off left by the direct solve, cleared
    at the end of each spray stage — and are now harmless, because a vanishing film's rate stays below V_s ℓ/A; the
    rates that remain are bounded by (V_s + |G| b²/(3 μ_l)) ℓ/A, which grows with the film, not as it vanishes. The deep
    transport is unchanged: its largest rate is 5.3·10⁵ per second at 0.0125 s and 4.6·10⁵ at 0.00625 s, as in the
    crashed run's probe (5.3·10⁵).

75. **The switched-step series with the fix: the masses settle within their scatter, the droplet population does not
    converge, even at 0.00625 s.** The 100 mm flight to 120 s at 0.5 s throughout and switched at 49.5 s to 0.05, 0.025,
    0.0125 and 0.00625 s, values in that order and over the continuum window 49.5–120 s unless stated. The scatter
    quoted beside a change is the spread between seeds 12345 and 1, measured for this fact at 0.05 s and at 0.0125 s
    (one pair each), because fact 65's four seeds were run at the default step only.
    * **The masses settle within their scatter, and the scatter grows as the step shrinks.** Sprayed mass 1.0585,
      1.0500, 1.0521, 1.0316 and 1.0353 kg; the body at 120 s 0.4136, 0.4220, 0.4200, 0.4404 and 0.4368 kg. Between the
      two seeds the sprayed mass differs by 0.008 % at 0.05 s but by 1.45 % at 0.0125 s (1.0316 against 1.0467 kg), and
      the mass at 120 s by 3.5 %: at 0.0125 s the two seeds' sprayed masses part by more than 0.1 % from 52.6 s and by
      more than 1 % from 72.6 s, as their fitted nose radii part between 70 and 90 s (a median of 81 against 72 mm) and,
      through the size feedback of fact 14, the stagnation heating with them (the convective heat integrated over
      70–90 s differs by 3.2 %). So the 2 % fall of the sprayed mass at 0.0125 s is round-off amplified by the shape
      feedback, not the step: from 0.05 s down the step-to-step changes (+0.2 %, −1.9 %, +0.4 %) lie within the spread
      at 0.0125 s, and the sprayed mass sits 0.6–2.5 % below the default step's.
    * **Runoff stays minor at every step, and its measures do not resolve between steps.** At most 6.1, 2.8, 3.2, 6.5
      and 7.2 % of the film formed inside any latitude cap leaves it, while the cap inside 15 degrees gains film net
      (2.5–13.5 % of what forms in it, the deceleration pushing the film toward the nose late in the flight); the
      mass-weighted shift between where film entered and where it was sprayed is 1.19, 0.67, 0.62, 0.89 and 0.83 degrees
      of latitude, against a spread of 0.13 and 0.18 degrees between seeds at 0.05 and 0.0125 s; the equatorial ring
      sprays 9.3, 7.9, 8.6, 8.8 and 9.2 % of the mass (spread 0.02 and 0.25 points). Fact 69's runoff question has the
      same answer at every step: the melt is sprayed within about a degree of latitude of where it entered the film, and
      less than a tenth of it is sprayed at the equator.
    * **The droplet population does not converge.** 15.5, 87.7, 106.6, 126.5 and 141.2 million droplets; median radius
      by number 181.6, 72.5, 71.9, 70.4 and 69.9 µm and by mass 190.9, 159.0, 131.1, 105.2 and 97.2 µm; the thick
      branch's share of the sprayed mass 92.8, 47.0, 35.3, 25.1 and 20.5 %; the mass-weighted mean film depth at release
      2.85, 0.33, 0.26, 0.23 and 0.21 mm. Per halving from 0.05 s the droplet count rises by 21.5, 18.6 and 11.6 %, the
      median radius by mass falls by 17.5, 19.8 and 7.6 % and the thick share by 11.7, 10.2 and 4.6 points; the last
      changes are 27, 16 and 29 times the spread between seeds at 0.0125 s (0.43 % in count, 0.47 % in the median by
      mass, 0.16 points in thick share). Only the median radius by number has settled: about 70 µm from 0.05 s down, its
      last change (−0.6 %) inside its 0.96 % spread. The thick share still rises through the flight at every fine step
      (at 0.00625 s from 17.7 % over 49.5–60 s to 35.0 % over 100–120 s).
    * **Within each branch the droplets do not depend on the step; what does is which branch takes the mass.** The thick
      branch's droplets have a median radius by number of 186, 177, 179, 186 and 186 µm (by mass 188–196 µm) and the
      thin branch's of 67, 71, 71, 70 and 70 µm (by mass 85–92 µm) at every step, while the thick branch's mass falls
      from 897 to 450, 339, 236 and 193 g and the thin branch's rises from 22 to 499, 617, 698 and 744 g. The droplet
      count follows the thin branch's mass, which makes about sixteen times as many droplets per gram. The branch split
      is set by fact 58's branch test, which reads a molten depth counted in whole elements at the end of the conduction
      step: the shorter the step, the less often the wall-owning element is above the top of the feed ramp when the
      spray runs, so the fewer patches are thick.
    * **Small quantities.** The front-surface Rayleigh–Taylor release is 47.4, 8.9, 4.1, 5.3 and 5.8 g (spread 17 % at
      0.05 s and 30 % at 0.0125 s, so from 0.025 s down it is 4–6 g within its scatter); the deep runoff 90.9, 0.85,
      0.84, 0.84 and 0.83 g; the re-solidification counter 6.6, 143, 286, 502 and 882 g, still 0.08–0.10 g per fine step
      — a count of freeze-and-re-melt cycles, not a result (fact 69).
    * **Against the pre-fix series** (fact 69), at the steps both have: at 0.05 s the fix moves the droplet count by
      +0.6 %, the median radius by number by −0.6 % and by mass by +0.7 %, and the thick share from 46.8 % to 47.0 %; at
      0.025 s by −3.9 %, +0.9 % and +4.5 %, and from 33.5 % to 35.3 %. The fix neither causes the step dependence nor
      removes it.
    **Verdict: the droplet population is not converged at 0.00625 s, and shrinking the step alone will not converge it
    at an affordable cost.** The last changes are smaller than the ones before (11.6 against 18.6 % in count, 7.6
    against 19.8 % in the median by mass, 4.6 against 10.2 points in thick share), the first sign of slowing. If they
    kept shrinking by those last ratios — 0.62, 0.38 and 0.45 per halving — they would fall below the spread between
    seeds after seven, three and five more halvings: a fine step between about 0.0008 and 0.00005 s, which is 90 000 to
    1.4 million fine steps for the 70.5 s after the switch, from about a day to three weeks of one core at the measured
    1.2 s per fine step. Two ratios do not establish a geometric sequence, so this is an order of magnitude at best.
    Followed to its limit, the same extrapolation points to a population in which the thin-film branch takes most of the
    mass after the switch — about 165 million droplets, a median radius by mass near 92 µm and a thick share near 17 % —
    numbers that show the direction of travel and are not results. Not measured: whether the part of the flight before
    the switch, computed at 0.5 s in every run, depends on the step too.

76. **Run times and power.** Every run of facts 72–75 was made on mains power under `caffeinate -i`; every log records
    "AC Power" at its start and end (the battery at 100 %, charged), from 22:47 on 2026-10-05 to 03:50 on 2026-10-06.
    The runs shared the machine (eight cores, four of them performance cores): eight at a time for the first twenty
    minutes, then four, three and two, so their wall times say as much about the load as about the model. The fix's own
    cost was measured like for like: the unamended and the amended copy run side by side on the 100 mm flight to 120 s
    at the default step, beside two other runs, took 821 and 823 s of model time (+0.2 %) — nothing measurable — and
    each reproduced its earlier run bit for bit, the unamended one the seeding amendment's default run (a check that the
    copy is what the diff blocks say it is). Wall times: the default-step flight 1 142 s eight at a time (761 s for the
    seeding amendment's run four at a time); the switched flights 3 605 s at 0.05 s, 6 169 s at 0.025 s, 10 968 s at
    0.0125 s and 18 139 s (5.0 h) at 0.00625 s; the two runs with seed 1, 2 985 s at 0.05 s and 8 557 s at 0.0125 s; the
    50 mm flight 412 s. A fine step cost a median 1.24 s at 0.00625 s and 1.66 s at 0.0125 s, falling through the flight
    as the body shrinks and the load eased: 2.4 s per step over 50–70 s, eight and four at a time, against 0.9 s over
    110–120 s with the run alone. The probe adds a few array reductions per sub-step and changes no result (it
    re-implements the transport operation for operation).

77. **Not done, and for Asha to decide.** (a) **The time step and the droplet population — decide this first.** Fact 75
    answers decision (2) of 2026-10-05: shrinking the step makes the masses settle within their scatter but not the
    droplet population, and the population's trend at 0.00625 s implies a fine step between about 0.0008 and 0.00005 s —
    a day to three weeks of computing per flight — before it would settle, if it does. Options: (1) keep shrinking: the
    next halving, 0.003125 s, costs six to eight hours alone on this machine and would show whether the slowing seen at
    the last halving continues, but would not reach convergence; (2) fact 61 (a)'s option (2), the branch test on a
    liquid depth by mass — the film, the deep account and the liquid inventory of the contiguous chain below the
    wall-owning element, divided by ρ_l A — which removes the whole-element counting that fact 75 identifies as what the
    step changes; (3) fact 61 (a)'s option (3), Girin's own statement of the regimes on rates, melting speed against
    stripping speed, which would make the thin branch the rule on this flight, the direction the series is travelling in
    anyway; (4) quote the droplet population only with its step, the series as its uncertainty. Recommendation: (2), as
    an amendment of its own with this series repeated (the five runs take about five hours in parallel with the harness
    as it is), and (3) as the check on whether the result is physical; until then (4): the masses may be quoted from
    0.05 s down to within about 2.5 % and the median radius by number as about 70 µm, and nothing else of the droplet
    population.
    (b) **Whether the Knudsen-based step switch should become a model option.** The series was made with a harness
    outside the package. Options: (1) leave it there until a converged step is known, the harness reproducing every run
    of the series meanwhile; (2) add `--dt-continuum <s>`, a second macro step applied from the first step at which the
    body Knudsen number falls below `KN_BODY_SHOCK`, latched, and recorded in the run name and JSON. Recommendation: (1)
    for now, since there is no step to give it; (2) once (a) is settled, because the part of the flight before the
    switch is 41 % of the flight time to 120 s, so the switch saves about 40 % of what the same fine step would cost
    over the whole flight — provided that part needs no fine step, which was not measured: every run of the series
    computes it at 0.5 s.
    (c) **The exit code of a model failure** (sub-plan 13's amendment of this date). The crash reached the user as exit
    code 2, "bad arguments", because the thermal solver's input check raises `ValueError` and `cli.main` maps every
    `ValueError` to 2. Options: leave it; raise `RuntimeError` (exit 1) from the solvers' input checks when the melt
    step calls them; or check the film for finite values at the end of each melt step and raise `RuntimeError` there,
    naming the stage. Recommendation: the last, with a unit test, in a change of its own.
    (d) **The re-solidification counter counts freeze-and-re-melt cycles, not lasting refreezing** (facts 69 and 75: it
    grows by 0.08–0.10 g per fine step whatever the step, from 143 g at 0.05 s to 882 g at 0.00625 s, against 6.6 g at
    0.5 s). The feed and the freeze-back are netted per element within a step (fact 25), but a film that freezes into
    its owner at one step and is fed out again at the next is counted every time. Options: (1) keep the counter and
    document it as gross; (2) count only the frozen-back mass still held in the elements at the end of each step (a
    per-element ledger debited when the element feeds again); (3) report the film's frozen share alone
    (`film_frozen_fraction` already exists). Recommendation: (2), in its own change; until then the column is not a
    result.
    (e) **The 87 g sprayed before the continuum switch** (fact 69). Before 49.5 s the flight is in the merged branch,
    where Girin's closure is not certified and the thin-film mode acts on the Couette closure; it sprays 87.4 g in every
    run of the series, all of which compute that part at 0.5 s. Asha questioned it on 2026-10-05 and it is left as it is
    pending her decision: whether the merged branch should spray at all is the gate's declared conservatism of the
    Global Constraints, not a question this amendment answers.
    (f) **Exact books for the film transport** (fact 72). Options: (1) leave it, its drift at most 3e-11 of the mass it
    moves per call and the cascade's books test held at 1e-11; (2) scale its arrivals to its departures, as fact 48 did
    for the deep stage, making the mass exact by construction, at the cost of changing every run in the last bits, the
    50 mm flight included. Recommendation: (2), in the same change as the next amendment that changes results anyway, so
    that one re-measurement covers both.
    (g) **Not measured, and recorded so it is not lost.** The resolved-mode verification rows (they transport a film
    under Girin's closure and so can move; their comparison needs Task 11's melting references); the whole 100 mm flight
    to the ground with the fix; the 50 mm flight at a fine step; why the front-surface Rayleigh–Taylor release falls by
    16 % at the default step and rises by 73 % at 0.05 s with the fix (a small, threshold-sensitive quantity whose
    scatter between two seeds is 17 % at 0.05 s); and fact 44's `PHI_DEATH = 0.50`, which every time-step study since
    2026-10-02 has left out.
    *Forward pointer (2026-10-07):* (a) is answered for the rule on, on Scheil's curve, by facts 92–94, and (b)'s step
    is proposed in fact 96 (a).
    *Forward pointer (2026-10-07, continuum step):* (b) is answered by fact 97 — `--dt-continuum`, 0.0125 s by default
    with Girin removal.

## Amendment of 2026-10-06 — the thin branch needs a rigid substrate (facts 78–87)

Asha's request of 2026-10-06: Girin's thin branch — his dominant ablation, the Girin & Kopyt (1994) mode, the case where
"the rigid core still stabilises the disturbances" — is valid only where the film rests on a rigid surface. A film on
slurry is deep melt and takes the thick branch; a film on coherent (semi-solid) mush or on solid AA7075 keeps the thin
one. Her decisions while the design was settled, the same day: (1) the boundary between coherent mush and slurry is
50 % liquid, the boundary Step 4 and the large-fragment design already use (Chen et al. 2016's semi-solid strength law
ends there and Li et al. 2014's slurry viscosity data begin there); (2) the regime test's melt layer runs from the top of
the film down to the first point at that boundary, and the thick branch fires where that layer is deeper than δ_m — a
slurry skin shallower than δ_m over rigid material keeps the thin branch (chosen over "any slurry beneath the film
triggers the thick branch"); (3) where there is no δ_m (off Girin's closure) a film whose base is slurry does not spray
by the shear modes at all — the strict reading, chosen over keeping the thin mode there; (4) the design as presented and
approved: lubrication's branch flag and surface velocity read the same layer as the spray's regime test, the
Rayleigh–Taylor criteria and the wave-fits region keep the liquid layer, the front-surface Rayleigh–Taylor mode is not
covered by (3), `--rigid-substrate on|off` defaults to on, and the change is recorded in history columns and a frame
field. It amends sub-plans 02 (`Material.T_rigid`), 06 (no code change), 07 (`spray.evaluate`), 09
(`MeltingBody.nonrigid_depth` and the melt step; one melting test and one two-backend test changed, one two-backend
test added), 10 (four history columns and a surface-frame field) and 13 (the flag). Sub-plans 14 and 15 are not amended
(fact 87 (d)). Facts 1–77 stand except where these say otherwise. **With the rule on, fact 28(a)'s decision that the
contiguous liquid depth decides the branch is replaced by fact 79's non-rigid depth, so fact 58's whole-element count
and double count no longer enter the branch test; decision 77(e) — whether the merged branch should spray — is answered
for films on slurry by (3), and facts 83 and 85 say what that does.** Measured in a throwaway copy
(`prototype/work-2026-10-06-rigid-substrate/`, git-ignored) of `prototype/proto3/` — unchanged since 2026-10-03, its 114
files plus six cached meshes — with the 24 diff blocks of the amendments of 2026-10-02 (9), 2026-10-03 (8), 2026-10-05
seeding (3) and 2026-10-05 runoff flux (4) applied in date order, every one exactly (no fuzz, no offset), and sub-plan
02's material code of 2026-09-28 (the Scheil variant and `AA7075-empiricaldata`) added for the Scheil runs; the
generator in sub-plan 02 rewrote `AA7075.json` and `AA7075_range.json` byte-identical to the copy's (`cmp`). Before any
change the copy reproduced fact 73's flights to the last digit (fact 81), so everything facts 49, 54, 62 and 69 list
about the copy holds: `AA7075_range` unless stated, US76, physics heating, the dense band, `PHI_DEATH = 0.05`, the deep
runoff and the molten cascade on, every run seeded with the default 12345. The six measurement runs were made from the
amended copy, never edited while they ran, all six at once, on mains power under `caffeinate -i`, the power source
logged at the start and end of each. (The baseline 100 mm run's log records battery power at its end — a short
interruption of the mains; its results match fact 73 to the last digit, as a seeded run's must.)

78. **The half-liquid temperature.** `Material.T_rigid` is the temperature at which the material is
    `RIGID_LIQUID_FRACTION = 0.5` liquid, found by bisection on the material's own liquid-fraction law (monotonic):
    829.0 K for `AA7075_range` (750 + 0.5 × 158), 895.10 K for `AA7075_scheil` (Scheil's formula gives 895.11 K; its
    1 K table, 895.10 K), 850.0 K for the single-temperature `AA7075` (the middle of its ±2 K ramp), and infinite for a
    material that does not melt. The slurry band between T_rigid and the liquidus is therefore 79 K wide on the linear
    law and 12.9 K on Scheil's, and that difference is the one the flights depend on most (facts 82 and 84).

79. **The non-rigid depth, and what reads it.** `MeltingBody.nonrigid_depth` follows a line from each patch's centre
    along its inward normal through the tetrahedra it actually crosses, leaving each by the face its barycentric
    coordinates reach first. The P1 temperature is linear along the line inside an element, so the point where it falls
    to T_rigid is found exactly rather than counted in whole elements. Only material continuous with the wall counts: the
    line stops at the first point at or below T_rigid, at the active mesh's boundary, or after
    `NONRIGID_MAX_CROSSINGS = 64` elements. Each element's stretch counts φ_e of its length, because what the element has
    fed to the film is already in the film account added on top — fact 58's double count, removed by an assumption:
    that the fed part came evenly from along the line (fact 87 (f)). Tested exact to 1e-9 relative on a linear field at
    0.4 mm (inside the second prism layer) and 3 mm (several elements down), halved to 1e-12 by φ = 0.5, and stopped at
    the first rigid point when the field is hot again deeper. The regime layer is the film depth plus the deep account's
    plus the non-rigid depth. It replaces `molten_depth` in four places only: `film.lubrication`'s branch flag and its
    surface velocity (both calls, so the runoff's flux and the spray's Weber number see the same branch), the spray's
    regime test, and the shear depth of Girin's Weber number on the patches the test makes deep (δ_m there, as on every
    thick patch). `molten_depth` still feeds the Rayleigh–Taylor criteria, the molten region the wave-fits test measures
    against, the deep runoff's chain, the molten cascade and the `molten_depth_*` columns. Off Girin's closure a film
    whose wall is above T_rigid (non-rigid depth > 0, `on_slurry`) takes neither the thin nor the rarefied mode; the
    front-surface Rayleigh–Taylor mode is not affected (decision (4)). With the switch off every one of these reads the
    liquid layer as before, bit for bit (fact 81).

80. **Tests.** Added: `test_the_rigid_temperature_is_where_half_the_material_is_liquid` (02); three spray tests — a thin
    film on slurry deeper than δ_m takes the thick branch with Girin's Weber number on δ_m while one whose slurry ends
    within δ_m stays thin, the Rayleigh–Taylor criteria keep the liquid layer, and off the closure a film on slurry takes
    no shear mode in the Couette and free-molecular branches alike (07); five melting tests — the three depth tests of
    fact 79, a 20 µm film over a millimetre of slurry at 70 km thick on every wet windward patch with the rule on and
    thin with it off, and at 30 s (no closure) the film held over slurry and sprayed over a rigid wall (09); and a
    two-backend test with the rule on (09). Extended: the CLI's run name, columns, settings and bad-argument tests (13)
    and the coupled run's frame fields (10). Changed, both because their scenarios make the whole body slurry under the
    linear law: (a) `test_film_spraying_death_and_balances` starts the body at 880 K, 82 % liquid, so its held film
    piles up and the front-surface Rayleigh–Taylor mode sheds droplets above the size histogram's 10 mm top edge (fact
    85); its check that the histogram holds every droplet released now compares the histogram with the released droplets
    inside its range and requires the rest to lie above it and below R/4 — the histogram drops out-of-range radii by
    design; (b) the 2026-10-05 two-backend test runs with the rule off, because with its interior at 850 K the two
    backends stop agreeing at the spray floor (fact 86) and it is the record of the runoff-flux agreement. Unit tier
    (`not drama and not reference`): 251 passed, 1 skipped, and the same 2 failures and 3 errors as the copy before the
    change (236 passed) — all five need Task 11's melting SESAM references, which the prototype does not have. The 15
    tests added are these 9, sub-plan 02's five Scheil and empirical-data tests, and one new bad-argument case. FEniCSx
    tier: 12 passed (11 before).

81. **Reproduction, and the switch off bit for bit.** The copy before the change reproduced fact 73: the 100 mm flight
    to 120 s sprays 1.0585 kg, leaves 0.4136 kg, writes 197 075 source-table rows, moves 90.9 g by the deep runoff and
    158.0 g by the cascade, re-solidifies 6.62 g and closes its energy balance to −8.5e-11; the 50 mm flight demises at
    199.0 s and 70.74 km with 399 history rows and 21 010 source rows. With `--rigid-substrate off` the amended code
    reproduces both bit for bit: all 241 and 399 rows of all 86 history columns the runs share, every result field but
    the run time, and all 22 source-table columns.

82. **The 100 mm flight to 120 s (`AA7075_range`), the rule on against off**, each change read against fact 65's
    scatter between four seeds: sprayed mass 1.0585 to 1.0708 kg (+1.2 %, scatter 0.13 %); the body at 120 s 0.4136 to
    0.4012 kg (−3.0 %, scatter 0.33 %); droplets 17.89 to 12.10 million (−32 %, scatter 2.0 %); median radius by number
    179.7 to 186.5 µm (+3.8 %, scatter 0.29 %) and by mass 190.6 to 199.5 µm (+4.6 %, scatter 0.16 %); re-solidified
    6.62 to 4.27 g (−36 %, scatter 1.4 %); largest droplet 5.3 to 5.2 mm, none above 10 mm. By branch, as shares of the
    sprayed mass: thick 85.1 to 88.2 % (901 to 944 g), thin 10.4 to 2.5 % (110 to 27 g), front-surface Rayleigh–Taylor
    4.5 to 9.4 % (47.4 to 100.1 g). The thin branch's droplets fall from 6.61 to 0.92 million, which is 5.69 million of
    the 5.79 million fewer droplets. Before Girin's closure (to 49.5 s) the flight sprays 87.3 g, all by the thin mode,
    with the rule off, and 74.8 g with it on — 16.4 g by the thin mode on walls below T_rigid (the thin mode's sprayed
    mass sits at a mass-weighted 70° from the stagnation point, against 52° with the rule off) and 58.5 g by the
    Rayleigh–Taylor mode; spraying begins at 28.5 s instead of 25.5 s, and the film held from the shear modes peaks at
    16.0 g at 48.5 s, one step before the closure begins (the film on the body peaks at 15.5 g against 6.4 g). Under the
    closure, as medians over its 142 steps: the thick share of the wet windward patches rises from 27.4 to 62.8 %; within
    the run with the rule on, 23.3 % of those patches are thick only because of slurry and 11.6 % are thin by the
    non-rigid depth where the whole-element liquid layer said thick (fact 58's count, removed); and the non-rigid depth
    under them is a median 8.1 mm, against a molten depth of 1.15 mm and a δ_m of 290 µm. Medians of fractions taken
    within one run do not add up to the difference between two runs, whose wet patches differ. The energy balance stays
    exact (−1.2e-10 of the absorbed heat), Newton takes 3.00 iterations per step (3.02), the peak surface temperature is
    975.8 K (981.1 K), the deep runoff moves 91.1 g (90.9 g) and the cascade feeds 156.9 g (158.0 g). The run took
    990 s against 973 s (+1.7 %), but all six runs ran at once, so the cost is not measured cleanly (fact 76).

83. **The 50 mm flight, which never has Girin's closure.** Decision (3) acts on the whole flight. Sprayed mass 0.1770 to
    0.1797 kg (+1.6 %); demise at 199.0 s and 70.74 km against 199.5 s and 70.58 km; droplets 1.917 to 0.555 million
    (−71 %); median radius by number 173.3 to 173.7 µm, but by mass 227 to 1 682 µm, 7.4 times larger; largest droplet
    1.8 to 3.9 mm. The front-surface Rayleigh–Taylor mode releases 0.3 g with the rule off and 119.4 g with it on, 66.4 %
    of the sprayed mass; the thin mode 176.6 g and 60.3 g, the latter from film that reached walls below T_rigid (a
    mass-weighted 56° from the stagnation point, against 47°). The held film peaks at 10.3 g at 191 s; spraying begins
    at 176.5 s instead of 174.5 s; re-solidified 3.7 mg against 12.2 mg; the energy balance 6.0e-11.

84. **The 100 mm flight with `AA7075_scheil`, the rule on against off.** Compare within the material, not across: Scheil
    also differs from the linear range in its latent heat (390 against 400 kJ/kg) and liquid surface tension (0.80
    against 0.86 N/m; sub-plan 02, 2026-09-28). Sprayed mass 1.1724 to 1.1789 kg (+0.6 %); the body at 120 s 0.2996 to
    0.2931 kg (−2.2 %); droplets 18.10 to 13.40 million (−26 %); median radius by number 181.9 to 186.9 µm (+2.7 %) and
    by mass 192.3 to 201.6 µm (+4.8 %); thick 85.7 to 88.4 % of the sprayed mass, thin 9.3 to 2.3 % (109 to 27 g),
    Rayleigh–Taylor 5.0 to 9.3 % (59.2 to 109.4 g); before the closure 93.8 g by the thin mode against 16.0 g thin and
    64.3 g Rayleigh–Taylor; the held film peaks at 17.3 g at 49.0 s. Under the closure (medians over 142 steps) the thick
    share of the wet windward patches rises from 32.2 to 55.9 %, but only 6.2 % are thick only because of slurry
    (23.3 % on the linear law) and 14.4 % are thin by the non-rigid depth where the liquid layer said thick; the
    non-rigid depth is a median 1.45 mm, 5.6 times thinner than the linear law's 8.1 mm, as the 12.9 K band against
    79 K predicts. So the slurry test itself matters four times less on Scheil's curve, while the no-closure rule acts
    the same. Energy balance 1.1e-10, Newton 3.51 iterations per step (3.51), re-solidified 5.23 to 4.55 g.

85. **Where the held film goes: the front-surface Rayleigh–Taylor mode.** Off the closure a film on slurry is held from
    the shear modes. Some of it runs to cooler walls, where the thin mode takes it (facts 82 and 83); the rest collects on
    the cap until its depth passes the Rayleigh–Taylor depth criterion, and the applied front-surface mode sheds it within
    one growth time, as droplets of its wavelength capped by the film mass on the patch. On the flights none exceeds
    5.2 mm. In `test_film_spraying_death_and_balances` (the body at 880 K throughout, 8 s at 71 km, no closure) the mode
    sheds 5 droplets above 10 mm carrying 57.8 g, 5.2 % of the 1.118 kg sprayed (a second, unseeded run of the same
    setting: 6 droplets, 10.1–10.9 mm, 6.1 %), from facets whose film is 18–34 mm deep by mass per area — blob facets
    where the lubrication picture has already failed (fact 27); with the rule off the same setting sprays 45.5 million
    droplets, the largest 3.2 mm, against 9.0 million. The masses barely move because the release is supply-limited
    (fact 9): what the rule changes is which mechanism releases the film, and so the droplet count and sizes.

86. **The two backends, and a threshold the rule sharpens.** In the 2026-10-05 two-backend test's setting — the 6 mm pool
    at 960 K over an interior at 850 K, ten steps of 0.05 s on the coarse mesh — the whole body is slurry under the
    linear law, and the non-rigid line runs 99.95 mm, through the sphere. With the rule on the two backends make the same
    branch decision on every patch for five steps and then differ on 8 to 29 of about 4 000; every difference is a film
    of 0.01–1 µm at the spray floor `B_MIN`, which Girin's thick branch releases whole (its rate does not depend on the
    film's depth) on one side and not the other. After ten steps the runoff differs by 1.6e-4, the sprayed mass by
    4.4e-7 and the temperatures by up to 0.29 K, against about 1e-8 and 2e-4 K with the rule off. With the interior at
    820 K, below T_rigid, the rule on agrees as tightly as ever: identical branch decisions at every step, mass 6.7e-10,
    sprayed mass 1.0e-8, runoff 1.2e-8, film 6.1e-8, temperatures within 1.8e-4 K — the new two-backend test. The
    threshold is the model's, not the backends', and the rule sharpens it because it puts more micron films on thick
    patches. How much it moves a whole flight was not measured: the flights of facts 82–84 ran on one backend, and the
    seed scatter of fact 65 was measured with the rule off.

87. **Not done, and for Asha to decide.** (a) **The front-surface Rayleigh–Taylor mode over slurry off the closure —
    decide this first.** Decision (3) holds the film only from the shear modes, and the Rayleigh–Taylor mode then sheds it
    as millimetre droplets (facts 83 and 85): two-thirds of the 50 mm flight's sprayed mass and a mass median radius 7.4
    times larger. Options: (1) keep it — the mode is a different mechanism, a deceleration-driven instability of a deep
    pool, and film piling up on slurry is such a pool; (2) hold the Rayleigh–Taylor mode as well over slurry off the
    closure — nothing would then spray before 49.5 s on the 100 mm flight and almost nothing on the 50 mm flight, whose
    film (already 10–17 g held) would ride the body until the closure, the deaths hand it on, it freezes, or the body
    demises with it; (3) revisit decision (3) and keep the thin mode there. Recommendation: (1), recorded as an
    assumption, because the mode's depth criterion and growth time are its own and it is the only mechanism this model
    has for a pool too deep for the thin mode; but the droplet sizes it gives on blob facets inherit fact 27's limit, so
    the 50 mm droplet population should be quoted only with that caveat until fact 27's blob issue is addressed.
    (b) **The material for quoted branch splits.** The slurry test's own effect depends on the liquid-fraction law: 23 %
    of wet windward patches thick only because of slurry on the linear law, 6 % on Scheil's (facts 82 and 84). Sub-plan
    02 already records the decision to make `AA7075_scheil` the melting default; until that is done, branch splits
    should be quoted with their material. (c) **The time-step series with the rule on** (fact 77 (a)). The non-rigid
    depth removes the whole-element count fact 75 found the step changing; the film's one-step supply remains in the
    regime layer. Recommendation: repeat the five-run series (about five hours) after (a). (d) **Sub-plans 14 and 15 are
    not amended**: the sensitivity script has no `rigidsubstrate-off` variant, and the README and
    `docs/model_assumptions.md` entries for the rule, `T_rigid` and the four columns are not written. (e) **The 64-element
    cap** is not known to bind on the flights; where the whole body is above T_rigid the line runs through it (99.95 mm,
    fact 86), which leaves the branch unchanged — anything deeper than δ_m is thick — but truncates
    `nonrigid_depth_mean_mm`. Not measured. (f) **The φ weighting assumes the fed part of an element came evenly from
    along the line**; in a wall-owning element it is in fact the hottest part, at the wall. Not measured.
    *Answered on 2026-10-07 (fact 88):* (a) option (1), the mode is not held; (b) `AA7075_scheil` is the melting
    default; (c) the series repeated on Scheil's curve with the rule on (facts 91–95); (d) sub-plans 14 and 15 amended.

## Amendment of 2026-10-07 — the Rayleigh–Taylor mode stays, Scheil's curve is the melting default, a freeze-back round-off fixed, and the time-step series with the rule on (facts 88–96)

Asha's decisions of 2026-10-07 on fact 87: (a) the front-surface Rayleigh–Taylor mode is **not** held back over slurry
off Girin's closure — option (1), so the film held from the shear modes there is still shed by it (facts 83 and 85
stand, and their caveat on the sizes it gives on blob facets with them); (b) `AA7075_scheil` becomes the melting
default; (c) the switched-step series is repeated with the rule on, on Scheil's curve, the new default (chosen over
repeating it on the linear range, which would have isolated the rule's effect from the material's); (d) sub-plans 14 and
15 are brought up to date. The series found a round-off in the deep runoff's freeze-back (fact 89), which is fixed here.
It amends sub-plans 02 (no code change), 09 (the freeze-back), 13 (the default), 14 (the drivers and the series) and 15
(the README, the assumptions and the spec, for this amendment and the one of 2026-10-06). Facts 1–87 stand except where
these say otherwise; **fact 87 (a), (b), (c) and (d) are answered here, and fact 77 (a)'s question — does the droplet
population converge in the step — is answered for the rule on by facts 92–94.** Measured in a throwaway copy
(`prototype/work-2026-10-07-scheil-default/`, git-ignored) made from the copy of fact 78's amendment, which that
amendment's ten diff blocks rebuild exactly from the earlier copy (checked file by file); the measurement copy was
frozen with a checksum manifest before the runs, re-frozen once for the fix of fact 89, and never edited while a run was
in flight; every run was made on mains power under `caffeinate -i`, the power source logged at the start and the end.

88. **What changed, in date order.** (1) The default (sub-plan 13): `--melt on` without `--material` selects
    `AA7075_scheil` — Scheil's curve between 750 and 908 K, 390 kJ/kg, Σ 0.80 N/m — instead of `AA7075_range`; the four
    other packaged materials stay selectable by name, and the help text lists all five. (2) The drivers (sub-plan 14):
    `melt_verification.py`'s physics mode and `melt_sensitivity.py`'s base run use `AA7075_scheil`, and the sensitivity
    table gains `range` (the linear law) and `norigid` (`--rigid-substrate off`). (3) The freeze-back (sub-plan 09, fact
    89). The CLI test of the melting run expects the Scheil default (it failed on the material first, as it should); one
    melting test is added (fact 89). Unit tier 252 passed, 1 skipped, with only Task 11's five known reference failures;
    FEniCSx tier 12 passed.

89. **A negative film from the freeze-back's round-off, and its fix.** The series' 0.00625 s run, first launched on the
    copy with only (1) and (2), stopped at its first fine step, 49.50625 s, with the thermal solver's input check "film
    mass must be one finite non-negative value per node" (exit code 2, fact 77 (c)); the six other runs, then about
    twenty minutes in, were stopped and the seven relaunched after the fix. A probe that re-ran the seeded flight and
    checked the film and deep accounts after every stage of the melt step found the first bad value right after
    `_freeze_back`: a film of −1.29e-25 kg on one windward patch 51° from the stagnation point, whose deep liquid the
    freeze-back had just taken whole; nothing was non-finite, and the non-rigid depth, `film.lubrication`'s outputs and
    the runoff transport were clean. The cause is the 2026-10-02 amendment's arithmetic: freeze-back takes the deep
    liquid first and debits the film by the capped amount less the deep part, m_f − ((m_f + m_d) − m_d), which in
    floating point is negative by a rounding unit for about half of all pairs whose deep liquid dwarfs the film (99 860
    of 200 000 random pairs), and the spray stage's clean-up of films below 1e-30 kg runs before freeze-back. The fix
    takes each account's part directly — from_deep = min(m_d, wanted), from_film = min(m_f, wanted − from_deep) — so no
    subtraction can round below zero; amounts where the cap does not bind are those of before. A unit test with one such
    pair (2.34e-14 kg of film over 1.88e-9 kg of deep liquid) failed before the fix and passes after it. **None of the
    2026-10-06 flights is changed by it:** re-run on the fixed copy, the 100 mm flight to 120 s and the 50 mm flight on
    the linear range with the rule on and off, the 100 mm flight on Scheil's curve with the rule off, and (as the
    series' 0.5 s run, fact 91) with the rule on are bit-identical to the runs of facts 81–84 in all 90 history columns,
    every result field but the run time, and all 22 source-table columns. So facts 81–86 stand as written.

90. **The default changes no explicit run.** The re-runs of fact 89 name `--material AA7075_range` or `AA7075_scheil`
    explicitly and reproduce the runs made under the old default bit for bit; the series' 0.5 s run, which names no
    material, reproduces fact 84's explicit `AA7075_scheil` run bit for bit, so the default selects exactly that
    material. The run name does not carry the material (fact 96 (c)).

91. **The series.** The 2026-10-05 harness, unchanged (`dtswitch.py`, recovered from that session's scratchpad and kept
    with the copy; fact 69): the 100 mm flight to 120 s at 0.5 s throughout, and switched at the body Knudsen number's
    first value below 0.01 — 49.5 s, as before — to 0.05, 0.025, 0.0125 and 0.00625 s, latched; seed 1 beside the
    default 12345 at 0.05 and 0.0125 s, for the scatter between seeds at those steps (fact 75 measured it the same way);
    the same command line as fact 75's runs but without `--material`, so `AA7075_scheil`, the rule on and every other
    setting at its default. Before 49.5 s the runs with the default seed are identical: each sprays 80.3 g, 16.0 g by
    the thin mode where the wall is below T_rigid and 64.3 g by the front-surface Rayleigh–Taylor mode (fact 84), and
    the step at 49.5 s, the first under Girin's closure and still 0.5 s long, 18.1 g (the seed-1 runs 80.1 and 18.2 g).
    Values below are over the continuum window 49.5–120 s unless stated, as fact 75's are, and each is read against the
    scatter between seeds measured at its own step.

92. **The masses settle within their scatter from 0.0125 s down.** At 0.5 s throughout and switched to 0.05, 0.025,
    0.0125 and 0.00625 s, in that order: sprayed mass 1.1789, 1.1458, 1.1602, 1.1654 and 1.1642 kg; the body at 120 s
    0.2931, 0.3262, 0.3118, 0.3066 and 0.3079 kg. Between seeds 12345 and 1 the sprayed mass differs by 0.12 % at 0.05 s
    and 0.14 % at 0.0125 s, the mass at 120 s by 0.42 % and 0.53 %. The last change, from 0.0125 to 0.00625 s, is −0.10
    % in sprayed mass and +0.40 % in the mass at 120 s, both inside the scatter at 0.0125 s; the default step sprays 1.3
    % more and leaves 4.8 % less at 120 s than the finest. Before the switch every run is the same (fact 91).

93. **The droplet population nearly converges, and from 0.0125 s it is within a few per cent.** Over the continuum
    window: droplets 12.97, 51.42, 69.24, 76.49 and 78.57 million; median radius by number 186.9, 69.5, 68.5, 68.5 and
    68.7 µm and by mass 196.8, 190.4, 188.6, 187.3 and 185.8 µm; the thick branch's share of the sprayed mass 94.8,
    80.2, 73.5, 71.2 and 71.4 %; the mass-weighted mean film depth at release 3.59, 0.245, 0.169, 0.134 and 0.102 mm.
    Per halving from 0.05 s the droplet count rises by 34.7, 10.5 and 2.7 %, the thick share changes by −6.7, −2.3 and
    +0.2 points, the median radius by number by −1.4, −0.1 and +0.3 % and by mass by −0.9, −0.7 and −0.8 %. Against the
    scatter between seeds at 0.0125 s — 0.74 % in count, 0.31 points in thick share, 0.16 % and 0.07 % in the medians by
    number and by mass — the last changes are 3.7 times, within, 2 times and 12 times. So the branch split and the
    median by number have converged at 0.0125 s, the count is converging fast (each halving's change a quarter to a
    third of the one before; one more halving would be expected to move it by about 0.7 %, within its scatter — an
    expectation from two ratios, not a measurement), and the median by mass drifts slowly, by under 1 % per halving
    (fact 94 says where from). Against the 2026-10-05 series (linear law, rule off; fact 75), whose last halving still
    moved the count by 11.6 %, the thick share by 4.6 points and the median by mass by 7.6 %, this is the convergence
    fact 77 (a) asked for; but two things changed at once, the rule and the material (decision (c)), so how much of it
    is the rule's was not measured. The default 0.5 s step remains unfit for the droplet population — 6 times fewer
    droplets than the fine steps, because the thin branch is almost absent there (11 g against 304 g) — while its masses
    are within a few per cent (fact 92).

94. **The branch split by mass converges; each branch's droplets nearly do.** The thick branch sprays 1 025, 840, 780,
    760 and 761 g over the window and the thin branch 11, 205, 279, 305 and 304 g: from 0.0125 s the split is fixed to
    within 2 g. The thin branch's droplets have a median radius by number of 72, 65, 66, 66 and 66 µm (by mass 156, 86,
    82, 81 and 81 µm) and the thick branch's of 188, 188, 188, 187 and 184 µm (by mass 194.6, 194.6, 194.4, 193.7 and
    192.0 µm); the thick branch's droplet count rises from 10.1 to 11.4 million at the last halving, which is where the
    last rise in the count and the drift in the median by mass come from — its own droplets getting slightly smaller at
    the finest step, not a shift of mass between branches. Why was not isolated (fact 96 (b)). The thick share still
    varies through the flight at every fine step — 72 % over 49.5–60 s, 60 % over 60–80 s, 78 % over 80–100 s and 97 %
    over 100–120 s at 0.00625 s — and those sub-window shares agree between 0.0125 and 0.00625 s to within a point.
    Under the closure the median share of the wet windward patches on the thick branch is 56, 36, 31, 27 and 24 %, of
    which 6, 33, 28, 25 and 23 % are thick only because of slurry; the non-rigid depth under them is a median 1.47,
    0.41, 0.34, 0.30 and 0.28 mm against a conjugate depth of 288–295 µm. At fine steps the non-rigid depth therefore
    sits on the conjugate depth, so the regime test is decided near its threshold, and the slurry under the film is a
    fifth of what the default step, which lets the surface overshoot the liquidus for half a second between feeds,
    reports.

95. **Runoff, small quantities, cost.** Runoff stays minor at every step: at most 13.0, 10.2, 8.3, 10.0 and 9.7 % of the
    film formed inside any latitude cap leaves it (fact 75: 6.1–7.2 % with the rule off on the linear law — the film on
    a patch made thick by slurry moves as the top of the conjugate layer, sub-plan 06's amendment of 2026-10-06); the
    mass-weighted shift between where film entered and where it was sprayed is 1.92, 0.98, 0.92, 0.97 and 0.96 degrees
    of latitude (0.95 and 1.02 at the second seed's 0.05 and 0.0125 s); the equatorial ring sprays 11.1, 7.7, 8.5, 8.9
    and 9.0 % of the mass. The front-surface Rayleigh–Taylor release over the window is 44.6, 2.9, 2.6, 2.7 and 0.5 g (a
    3 % spread between seeds at 0.05 s and 31 % at 0.0125 s: a few grams, within its scatter from 0.05 s down), on top
    of the 64.3 g it sheds before the switch; the deep runoff moves 116.9 g at 0.5 s and 0.88–0.90 g at every fine step;
    the re-solidification counter reads 4.5, 107, 238, 486 and 941 g, still about 0.08 g per fine step —
    freeze-and-re-melt cycles, not a result (fact 77 (d)). The energy balance is 1.1e-10 at 0.5 s and grows with the
    number of steps to 4.1e-9 at 0.00625 s, as before (fact 74). All seven runs started together at 14:16 on 2026-10-07
    and took 0.37, 1.30, 2.14, 3.74 and 5.92 h of wall time (the seed-1 runs 1.30 and 3.75 h), seven, then six, four,
    three and one at a time on the eight-core machine; every log records mains power at its start and end. The 0.00625 s
    run took 21 298 s against 18 139 s for fact 75's under a different load, so the rule's cost at fine steps was not
    measured cleanly.

96. **Not done, and for Asha to decide.** (a) **The step to quote droplet populations at — decide this first.** From
    0.0125 s the population is within a few per cent of the finest step's, its branch split within its scatter (fact
    93); 0.0125 s costs 3.7 h of wall time for the 100 mm flight to 120 s with up to seven runs sharing the machine.
    Options: (1) quote the population from a 0.0125 s run made with the harness, labelled with its step; (2) make the
    switch a model option — fact 77 (b)'s `--dt-continuum <s>`, a second macro step from the first step at which the
    body Knudsen number falls below `KN_BODY_SHOCK`, latched and recorded in the run name and JSON — now that there is a
    step to give it; (3) keep 0.5 s and quote only masses from it. Recommendation: (2) with 0.0125 s, in its own change
    with a test that the option reproduces the harness's run bit for bit, because every quoted number should come from
    the model's own command line; masses may meanwhile be quoted from the default step within about 5 % (fact 92). (b)
    **The thick branch's own droplets shrink slightly at the finest step** (median by number 187 to 184 µm, fact 94),
    which is the remaining drift in the median by mass; a seed-1 run at 0.00625 s (about 6 h) would show whether it is
    outside the scatter there. (c) **The run name does not carry the material** (fact 90): with the Scheil default a run
    with `--material AA7075_range` and a default run of the same settings have the same name and overwrite each other in
    one output directory — a gap since Step 3 began (it held for `AA7075` against `AA7075_range` too) that the new
    default makes easier to hit, against the repository's rule that run names encode the whole configuration. Options:
    append `_mat-<name>` when the material is not the mode's default, so existing default names stay; or leave it, every
    driver already passing `--name`. Recommendation: the suffix, in its own change. (d) **Not measured.** How much of
    fact 93's convergence is the rule's rather than the material's (a Scheil series with the rule off, five runs, about
    6 h); the 50 mm flight at a fine step; fact 44's `PHI_DEATH = 0.50`; and every physics and sensitivity row of Task
    14 on the new default (sub-plan 14's amendment of this date).
    *Answered on 2026-10-07 (continuum step, fact 97):* (a) option (2), `--dt-continuum`, 0.0125 s by default with
    Girin removal.

## Amendment of 2026-10-07 (continuum step) — the macro step shortens where Girin's closure begins (facts 97–100)

Asha's decision of 2026-10-07 on fact 96 (a): the switched step becomes a model option, `--dt-continuum`, at 0.0125 s,
option (2). It amends sub-plans 10 (the switch in the coupled loop), 13 (the flag), 14 (the drivers) and 15 (the README,
the assumptions and the spec). Facts 1–96 stand except where these say otherwise; **fact 96 (a) is answered, and fact 77
(b)'s question — should the switch become a model option — with it.** Measured in a throwaway copy
(`prototype/work-2026-10-07-dt-continuum/`, git-ignored) made from the copy of facts 88–96, which that amendment's seven
diff blocks rebuild exactly from the earlier copy; the measurement copy was frozen with a checksum manifest before the
runs and never edited while one was in flight; every run was made on mains power under `caffeinate -i`, the power source
logged at the start and the end. The sensitivity driver's `dtcoff` variant (sub-plan 14) was added to the code after the
copy was frozen; no run uses the driver.

97. **The option.** `CoupledSettings.dt_continuum` and `kn_switch`, and `--dt-continuum off|<s>`. After each macro step,
    on the aero state the body was given, the step becomes `dt_continuum` the first time the trajectory's body Knudsen
    number — the mean free path over the body's current reference length, the history's `knudsen` column — is below
    `--kn-body-shock` (0.01, the surface flow's own continuum boundary, where Girin's closure begins), and stays there.
    That is the harness's rule (fact 69). Its default depends on the removal mode, as `--size-feedback`'s does: 0.0125 s
    for `--removal girin`, the model proper, where the series converges (facts 92–94); off for `--removal instant`, the
    bookkeeping device, which has no film and whose SESAM thresholds were measured at the default step. Frames count
    steps until the switch and flight time after it. The results record the value and the switch's time, Knudsen number
    and altitude. Only a non-default value changes the run name (`_dtcontinuum-off`, `_dtcontinuum-<s>`); a value of
    zero or below, one longer than `--dt`, or the flag without `--melt on` exits 2.

98. **Tests.** A coupled test runs a hot body on the coarse mesh from 48 s through the switch twice — with the option
    and with the harness's rule applied by hand — and requires every history column to be identical, the switch recorded
    where the harness switched (48.5 s there, the body's Knudsen number reading its shrinking size), the steps 0.5 s
    before it and 0.05 s after, and frames by step then by flight time. The CLI tests check the default per removal
    mode, `off`, the run names and three bad values. `test_melting_run_writes_columns_files_and_json` crosses the
    continuum boundary at about 9.5 s of its 15 s from 71 km, so with the default it ran 460 macro steps instead of 30;
    it now passes `--dt-continuum 0.1` and checks the switch it records. Unit tier 258 passed, 1 skipped, with only Task
    11's five known reference failures; FEniCSx tier 12 passed. Writing the coupled test showed one thing that is not
    the option's: where fine steps sum to just short of `t_max`, the trajectory stepper ends the run with a sliver step
    (1.1e-13 s in the test; the series' 0.05 and 0.00625 s runs end the same way), which adds one history row and
    changes nothing measurable (fact 100 (b)).

99. **The option reproduces the harness's runs bit for bit.** Run with the defaults — `AA7075_scheil`, the rule on,
    `--dt-continuum` at its default — the 100 mm flight to 120 s switches at 49.5 s, at 69.93 km and a body Knudsen
    number of 0.009905, and is bit-identical to the series' 0.0125 s harness run (fact 91): all 5 740 rows of the 90
    history columns, every result field the harness run has but the run time and the frame count, all 22 source-table
    columns, and its 13 frames at the same flight times. Its energy balance closes to 3.1e-9. With `--dt-continuum off`
    the same flight is bit-identical to the series' 0.5 s run — all 241 rows of the 90 history columns, every result
    field the older run has but the run time and the frame count, all 22 source-table columns, and its 13 frames at the
    same flight times. The 50 mm flight on the linear range never crosses the continuum boundary — it records no switch
    — and is bit-identical to the run of fact 89 with the option at its default. So the harness can retire: a step
    series is a set of ordinary runs. The cost: the default run took 8 207 s (2.3 h) of model time for 5 739 macro
    steps, 5 640 of them fine, about 1.4 s each on average and 1.9 s early in the window, while the body is large
    (measured from its frames: 10 s of flight in 26 min between 50 and 60 s); the run with the switch off took 643 s
    (10.7 min), but it shared the machine with two other runs for most of its time, so the ratio between them, 12.8, is
    approximate.

100. **Not done, and for Asha to decide.** (a) **The cost of a whole flight — decide this first.** A whole 100 mm flight
     now switches at 49.5 s and keeps the fine step to the ground. The committed non-melting references reach the ground
     at 366–368 s, so a whole flight is about 25 000 fine steps; at the 1 to 1.9 s a fine step costs (fact 99; the
     larger figure early in the window, while the body is large) that is about 7 to 13 hours, against about 40 minutes
     at 0.5 s (fact 23), and every 100 mm row of Task 14's tables costs it. The spray has nearly stopped by 120 s — at
     0.0125 s it releases about 20 g/s from 60 to 100 s, 7.5 g/s over 100–110 s and 1.5 g/s over 110–120 s, the hottest
     surface at 910 K, the top of the feed ramp — so most of those fine steps fall where nothing sprays. Options: (1)
     keep it and budget the runs; (2) return to the default step once the melting is over — for example from the first
     step after the switch at which nothing is fed or sprayed and no surface node is above the half-liquid temperature —
     as a change of its own, with a measurement that the remnant's state at the ground does not depend on the step; (3)
     end Task 14's physics runs at 120 s, which loses the survivor's state at the ground. Recommendation: (2), designed
     and measured on its own; until then the step-sensitive droplet quantities, which are settled by 120 s, can be
     quoted from runs to 120 s. (b) **The sliver last step** (fact 98): snap the trajectory stepper's last step to
     `t_max` when what remains is below a small fraction of the step; it would change the last history row of every run
     that ends this way, so it belongs in a change of its own with the runs re-measured. (c) Facts 96 (b) — the thick
     branch's drift at the finest step — and 96 (c) — the run name does not carry the material — stand.

## Amendment of 2026-10-08 (frame export) — the frames carry what the Spheral replay reads (facts 101–107)

The large-fragment model's core (Spheral M1, plan `docs/superpowers/plans/2026-10-07-spheral-m1.md`) was built on frames
written by a reconstruction of the prototype with five write-only export additions (`prototype/work-2026-10-08-spheral-mvp/`,
its `rebuild/post/`). M1's decisions 4, 9 and 13, and Asha's of 2026-10-08 that a run name always carries its material,
ask Step 3 to adopt four of them. This amendment ports them onto the current copy. It amends sub-plans 09 (the body keeps
the step's flow evaluation and evaluates the flow for the record), 10 (the history column and the frame fields), 13 (the
run name and the run JSON) and 15 (the README, the assumptions and the spec), and puts a forward pointer in sub-plan 18.
Facts 1–100 stand except where these say otherwise; **fact 96 (c) and fact 100 (c) — the run name does not carry the
material — are answered.** Made in `prototype/work-2026-10-08-frame-export/` (tracked; run outputs ignored). Its base is
the continuum-step copy of facts 97–100 (`prototype/work-2026-10-07-dt-continuum/code/`, tracked); `code/` is the same
with the additions, and the five diffs in `amendment/` rebuild `code/` from that base exactly. The measurement copy `meas/` was frozen with a checksum manifest before its run, which was
made on mains power under `caffeinate -i`.

101. **The wall loads on every evaluated step** (addition 02).
     - The surface flow is evaluated at the top of every melt step with an aero state, film or not, but only the spray
       step kept it, as `last_flow`, which is None without film.
     - So the frames carried `p_w`, `tau` and `closure` only from the first film, and the history's `p_w_stag_Pa` is NaN
       before it: rows 0–47, 0–23.5 s of the 100 mm Scheil flight.
     - The body now keeps the step's own evaluation and its aero state. The frames write the gas-side fields (`closure`,
       `kn_local`, `p_w`, `tau`, `delta_m`) from it.
     - The history gains `p_w_stag_step_Pa`, every step's own stagnation wall pressure, NaN on row 0 only. The spray
       step's `p_w_stag_Pa` is kept as it was.

102. **The faces a step's deaths exposed** (addition 06).
     - With the molten cascade a step's element deaths expose many faces at once (1,112 of 18,616 at frame 100 of the reconstruction's
       flight, a third of the nose), which the frames carried as "not evaluated".
     - `MeltingBody.flow_for_the_record()` evaluates the same surface flow at the step's aero state on the current
       surface, for the frame only. Its one side effect, `SurfaceFlow.last_bins`, is read by nothing.
     - The exposed faces take the gas-side fields from it, and `flow_eval` says which evaluation a patch carries: 1 the
       step's, 2 the record's, 0 none (frame 0; `--removal instant`).
     - The spray fields keep their defaults there, since nothing sprayed from a face that did not yet exist.

103. **The derived surface's normals and the flight direction** (additions 04 and 02).
     - Every surface frame carries `n_derived`, each patch's outward unit normal on sub-plan 01's derived surface (its
       Taubin-smoothed normal, matched to the body's patches by face id).
     - The Spheral core marches its layer depths along it: along the staircase facets' normals 43 % of the rays left the
       body before f_l fell to 0.5 (M1 decision 9).
     - The triangles keep the face table's node order, whose winding is inward on 8–18 % of the patches; `n_derived` is
       outward by construction.
     - The run JSON's settings record `v_hat_body`, (1, 0, 0), and `freestream_velocity_direction_body`, its negative.

104. **The run name always names the material** (addition 03, made unconditional).
     - A run with a thermal model carries `_material-<name>` after the heating and melt parts, default or not; a material
       file is named by its basename.
     - Runs without a thermal model keep their names. Task 14's drivers pass `--name` and are unaffected.
     - It renames every thermal run, the default ones included, so a default run made before this amendment no longer
       shares a name with one made after it.
     - `resolve_material` chooses the default in one place, for the name and for the body.
     - The reconstruction's addition 03 named only a non-default material; Asha asked for it always (2026-10-08).
     - Addition 05, a bound on the runoff's emptying number, is not ported: the runoff-flux amendment (facts 69–77)
       replaced the flux it bounded.

105. **Tests.**
     - Unit tier: 261 passed, 1 skipped, with only Task 11's five known reference failures. The base copy gives 258
       passed and the same five failures.
     - Two new coupled tests:
       - A cold body for 2 s from 69.8 km: no film, yet the loads appear on every frame after the first, equal to the
         step's evaluation, and `flow_eval` is 1 everywhere.
       - A hot body for 6 s with deaths every step, run with frames at every step and without: every history column is
         identical, and the record evaluation ran. Both runs are seeded as the CLI seeds; without the seed, pyamg's
         random starts alone move the altitude by 1e-12.
     - New CLI test: `test_run_name_always_names_the_material`.
     - The FEniCSx tier was not run: `fenicsx_env` is not installed on this machine, and nothing here touches a thermal
       backend.

106. **The 100 mm Scheil frames flight on this copy** (`launch.sh`).
     - Configuration: the reconstruction's flight flags (`--material AA7075_scheil --removal girin --deep-runoff off
       --molten-cascade on --frames-every 1`) plus `--dt-continuum off`. That is the 0.5 s step Asha asked for the
       pipeline's development, on the current physics: runoff flux, rigid substrate, the freeze-back fix.
     - Name:
       `model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off_dtcontinuum-off`.
     - The run: 1,094 s of model time for 1,262 macro steps; the ground at 630.6 s; 1,263 frames, 5.1 GB.
     - Melting: onset at 24 s (74.2 km), spraying from 25 s. Sprayed 1.1735 kg of 1.4720, a remnant of 0.298 kg,
       1.32e7 droplets (median radius 193 µm), cascade mass 0.323 kg. Energy and melt balances close to 6.3e-9.
     - The frames carry the loads on every frame but frame 0, and `flow_eval` 2 on 201 frames (up to 1,493 faces). The
       largest `p_w` equals `p_w_stag_step_Pa` to a median 6e-10, except where the stagnation patch died that step
       (worst 6.2 %).
     - Release: Σ release_rate · A over a frame's surface is 22 % below the history's sprayed mass over the flight
       (0.913 of 1.174 kg; −40 % to +13 % per step). On the reconstruction's flight it was 12 % (1.025 of 1.169 kg).
     - Spheral M1's `prepare` reads it through its contract with nothing missing, nothing optional
       absent and every gating check passed, on the package of `meas/` (SHA-256 64d949ea…). It processed the 1,263
       frames in 1,190 s and wrote 801.5 MB.
     - Measured by that prepare: mass to 4.3e-9 of the history; open edges 0 once wound outward; f_l bitwise; one
       frame without loads; drag binning 9.3e-4, above the 8e-4 that M1's flight test set on the reconstruction's
       frames; the drag against the history worst at 34 %, which is minor for Asha (M1 decision 15).

107. **Not done, and for Asha to decide.**
     - (a) **The release on removed faces.** `release_rate` is the spray step's, carried onto the frame by face id, so
       what a step released from the faces its own deaths then removed is in the history's sprayed mass but on no
       frame's surface: 22 % of the flight's spray. A Spheral sink that reads `release_rate` (its spec §10, milestone
       MC) inherits that deficit. Options:
       - (1) write the removed faces' release onto the faces their deaths exposed, by the rule the film's hand-over at element
         deaths uses (`MeltingBody._hand_over`);
       - (2) write it as a per-frame scalar (`release_removed_kg`) beside the per-patch field;
       - (3) leave it, since the sink belongs to milestone MC, which is not planned yet.

       Recommendation: (2) now (one number per frame, write-only); (1) if MC needs the distribution.
     - (b) Sub-plan 18's Task 8 was tested on the copy before this one, so its CLI diff merges by hand, and its expected
       run names gain the material (its forward pointer).
     - (c) Spheral M1's thresholds stay on the reconstruction's frames (its decision 13). These frames supersede them as
       Spheral's input once Asha chooses; `prepare`'s contract checks pass on both.

Dependency direction (spec §4): `spray` → `surface_flow`, `dispersion`; `body` → `film`, `spray`, `surface_flow`, `thermal`, `material`; `coupled` → everything; `viz`, `compare` read exported files and histories only; `girin_case` → `dispersion`, `surface_flow.ranger_psi`.

---

