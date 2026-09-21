# Melting, melt film and melt spraying of the re-entering sphere (Step 3 of the physics model) — Design

Date: 2026-09-20
Status: design for review, awaiting implementation plan
Builds on: `2026-09-17-reentry-trajectory-model-design.md` (Step 1, merged at df8aa65) and
`2026-09-18-thermal-fem-design.md` (Step 2, merged at c27da07; amendments in its §14). The assumptions and methods of
Steps 1–2 are summarised in `docs/model_assumptions.md`.

## 1. Purpose

Add to the coupled trajectory + 3D conduction model of Step 2 the three things that turn a heated sphere into a
particle source: **melting** of the AA7075 body on the finite-element mesh, the **melt film** that forms on the surface
(its thickness, its shear-driven runoff along the surface), and **melt spraying** — the quasi-continuous release of
droplets from the film by the hydrodynamic gradient instability of Girin (2017; thin-film case Girin & Kopyt 1994).
Deliverables:

1. the coupled run's history with mass loss (the trajectory now feels it), film and spraying columns;
2. a **particle source table**: every droplet release with its time, place on the sphere, size, number, mass and
   release velocity, plus size distributions Δn(r), ΔM(r) — the input of a later wake step (Step 4);
3. the 3D visualization of temperature, melt front, film thickness and release sites over the flight;
4. verification against SESAM's mass loss (two new melting references), Girin's published cases, and analytic
   solutions, every comparison plotted with its residual and its reference data saved in the repository.

Droplets are recorded at birth and not tracked afterwards (decided 2026-09-20: option (a) of the brainstorm; their
fate — secondary breakup, cooling, the wake — is Step 4).

## 2. Facts this design relies on

| Topic | Fact | Source |
|---|---|---|
| Girin 2017 dispersion relation | (Ω − Δ)[(Ω − Δ)(Ω + KΔ) + (1 − K)Δ] = Δ³ We_s⁻¹ (Ω − KΔ), K = [1 − e^(−2Δ)]/(2Δ), Ω = ωδ_m/V_s, Δ = 2πδ_m/λ, We_s = ρ_m V_s² δ_m/Σ. Unstable above We_s,cr = 3.08 (4.62 recommended in practice); for We_s ≥ 20: Δ_f = 1.225, Im Ω_f = 0.24, Re Ω_f = 1.5 Im Ω_f. **Solved numerically 2026-09-20 (cubic in Ω, scan in Δ): instability appears between We_s 3.00 and 3.08; Δ_f → 1.228 and Im Ω_f → 0.246 at large We_s; Re/Im 1.6; at We_s = 4.62: Δ_f 0.58, Im Ω_f 0.049; at 10: 1.34, 0.18; at 20: 1.29, 0.22.** So Fig. 5 of the paper is a computation, not a transcription. | Girin 2017 Eq. (1), Fig. 5 |
| Conjugated boundary layers | δ_m = (α/μ²)^(1/3) δ_a, V_s = (αμ)^(1/3)/p · V_a, p = 1 + (αμ)^(1/3), α = ρ_a/ρ_m, μ = μ_a/μ_m (Eq. 2, Appendix B); gas boundary layer on a liquid sphere δ_a = 2.2 R Re_a^(−1/2) Ψ(φ), Ψ = [(6φ − 4 sin 2φ + ½ sin 4φ)/sin⁵φ]^(1/2) (Ranger 1972), potential-flow edge velocity V_a = 1.5 V sin φ. | Girin 2017 §3 |
| Droplet size and period | r = k_r λ_f = k_r 2πδ_m/Δ_f, k_r = 0.17; t_per = k_t δ_m/(V_s Im Ω_f), k_t = 1.1 (calibrated on drop-shattering experiments, Girin 2011a); tori of cross-section r, one per wavelength, break into droplets of radius r (Eqs. 5–8). | Girin 2017 §4 |
| Two regimes | Outstripping melting (b > δ_m): the above. Outstripping ablation (b < δ_m): the film is thinner than the melt boundary layer, the rigid core stabilises the disturbances, droplets scale with the film; the thin-film analysis of Girin & Kopyt 1994 applies. | Girin 2017 §1 |
| Thin-film analysis (inviscid film on a rigid base, compressible gas) | Side surface, KH-type: dominant wavelength λ* = 1.5 M d/We_d with We_d = ρ₂ V₀² d/σ (d film thickness, M gas Mach number), dissipation cut-off λ_t = 0.5 M d/We_d, growth time τ* (their Eq. 12); r_d = λ*/4, period τ_d = τ*. Front surface, Rayleigh–Taylor from the deceleration W: λ* = 2π(3σ/(Wρ₁))^(1/2), ω* = i(4W³ρ₁/(27σ))^(1/4), needs W d² ρ₁ > 3σ. Their Table 1 (r_d, τ_d, ṁ for ρ₂ = 10⁻⁷–10⁻⁵ g/cm³ × V₀ = 3, 10 km/s) and Table 2 (RT λ*, τ* for W = 10⁵–10⁷ cm/s²) are the transcription tests; the authors restrict the results to continuum flow. | Girin & Kopyt 1994 Eqs. 9–14, Tables 1–2 |
| Girin 2017 published cases | R₀ = 0.3 cm; variant I iron 60 km/s: GI 13.0, N = 1.5·10⁶, sizes 4–135 µm, r_med 26.9 µm, t_s.d. 5.9 ms, σ = 1.04·10⁻¹² s²/cm², φ_cr 16.1°, t_f 5.7 µs; variant II iron 25 km/s: N 3.7·10⁵, t_s.d. 16 ms, φ_cr 31.4°; variant III stone 60 km/s: N 832, t_s.d. 149 ms. Assumptions: ρ_a = 10⁻⁷ g/cm³ constant, W(τ) = 1 − exp(−Cτ) with C = 2α^(1/2), t_ch = 2R₀/(α^(1/2)V∞), ρ_ir 7.8, ρ_st 3.5 g/cm³, μ∞ 6.8·10⁻⁵, μ_ir 5.8·10⁻³, μ_st 0.174 kg/m/s, Σ_ir 1.2, Σ_st 0.36 N/m. | Girin 2017 Table 1, §5 |
| Our case (100 mm sphere, 60–77 km, 6–7.5 km/s), estimated 2026-09-20 | Post-shock edge state ρ_e ≈ 4·10⁻⁴–2·10⁻³ kg/m³, u_e ≈ 1 km/s mid-sphere, μ_e ≈ 1.7·10⁻⁴ Pa s → Re_e ≈ 300–1200, δ_a ≈ 3–7 mm; α ≈ 2·10⁻⁷, μ ≈ 0.13 → δ_m ≈ 0.07–0.15 mm, V_s ≈ 3–5 m/s (also from τ ≈ μ_e u_e/δ_a ≈ 25–60 Pa and τδ_m/μ_m); We_s ≈ 3–5 — at the threshold; continuum droplet radius k_r 2πδ_m/Δ_f ≈ 90 µm; melt-film equilibration b²/α ≈ 0.3 ms (thermally thin); Kn_δ = λ_e/δ_a ≈ 0.03 at 71 km (slip), < 0.01 below ~65 km; free-molecular shear ρ∞V∞² sin θ cos θ ≈ 1850 Pa at 71 km, 45° (much larger than the continuum 24 Pa — the same ordering as C_D,fm > C_D,c); pressure-gradient runoff term G b²/(2μ) ≈ 25 % of the shear term at b = 0.1 mm, dominant above ~0.3 mm; body deceleration ≈ 30 m/s² → the 1994 RT front mode (needs 10³–10⁵ m/s²) is inactive; droplet freestream Weber number We_d = ρ∞V²·2r/Σ ≈ 1–5 (< 12: no secondary breakup expected). | this design |
| SESAM melting (lumped) | Mass removed at Q_net/L_f once the whole body is at 850 K (drama-AA7075, L_f 400 kJ/kg); 100 mm: onset 71.0 km (Kn 0.011), fully melted 66.5 km; 50 mm: 74.9 → 73.0 km; a fully melted sphere ends "ballooning" with ~3 mg residual; `thick_mm` is the sphere radius in mm and shrinks with the mass. | facts note §§8–9 |
| Step 2 state | Surface temperature at peak heating 2319 K (physics) vs mean 1717 K for the no-melt 100 mm sphere; 2 mm surface elements carry a ~30 K gradient across the first element at 2 MW/m² (0.25 mm layers: ~4 K); step cost 0.10 s at 18.9 k nodes (skfem), 0.31 s (FEniCSx). | Step 2 spec §14, README |
| Liquid AA7075 (assumed, near the liquidus) | ρ_l = 2400 kg/m³, μ_l = 1.3 mPa s, Σ = 0.86 N/m (pure-aluminium values; alloy corrections ≤ 10 %); solidus 750 K, liquidus 908 K, L_f = 400 kJ/kg (DRAMA's value; handbook 380). | Smithells Metals Reference Book; ASM Handbook Vol. 2 (to be cited in the material file) |

## 3. Scope

In: prism-layer surface mesh (optional, default on for melting runs); enthalpy-method melting with the two materials;
per-patch film state with lubrication runoff and transport; gas-side edge state, laminar boundary layer, Kn_δ
regime, shear and pressure gradient per patch; Girin's instability (dispersion relation solved; thick, thin and
rarefied branches); droplet release bookkeeping and the source table; element fractions and element death; mass and
projected area fed back to the trajectory; new history columns, plots, videos; verification against two melting SESAM
references, Girin 2017/1994 tables and analytic solutions; sensitivity study; CLI; README and `docs/model_assumptions.md`.

Out (later steps): droplet tracking, secondary breakup, cooling/solidification and the wake (Step 4); the within-patch
size spread (§17); mushy-zone rheology beyond the coherency threshold (§17); vaporisation and blowing; oxide skin;
tumbling or attitude change; nose-radius and shape feedback on heating and drag tables; non-spherical initial shapes;
MPI runs.

## 4. Architecture

```
reentry_model/dispersion.py    Girin Eq. (1): fastest unstable root -> Δ_f(We_s), Im Ω_f(We_s), Re Ω_f; cached table
reentry_model/surface_flow.py  per patch: edge state, laminar boundary layer (Thwaites), Kn_δ, regime, shear τ, gradient G
reentry_model/film.py          per patch film mass; lubrication velocity and flux (thin/thick branch); upwind runoff transport
reentry_model/spray.py         instability per patch (continuum thick/thin, rarefied thin), release bookkeeping, source rows
reentry_model/body.py          MeltingBody(ThermalBody): element fractions, film feed, element death, mass/area/energy accounting
reentry_model/mesh.py          active-element mask, prism layers, box mesh; surface() from active tets (already generic)
reentry_model/thermal/*        element-fraction scaling of the operators; pattern rebuild on death
reentry_model/material.py      latent heat in enthalpy(); liquid_fraction(T) ramp; liquid properties; AA7075 and AA7075_range
reentry_model/coupled.py       melt step after the thermal step; particle/histogram writers; demise end reason
reentry_model/viz.py, compare.py, cli.py, analysis/   overlays, mass-loss comparison, flags, Girin driver
```

Dependency direction: `spray` → `surface_flow`, `dispersion`, `film`; `body` → `film`, `spray`, `thermal`,
`material`; `coupled` → everything; `viz`, `compare` read exported files only. Every module is plain arrays over
patches or elements with its own tests.

**Data flow per macro step**: trajectory advance → aero state → heating (Step 2) → thermal step → melt step:
(i) liquid inventory of the surface elements from the enthalpy field → film feed; (ii) `surface_flow` from the gas
state and the current geometry; (iii) `film`: velocities/fluxes, runoff transport (sub-stepped); (iv) `spray`:
regime, We_s, release from unstable patches limited by the film; (v) `body`: φ_e update, element death, geometry
refresh, mass/area for the next advance; (vi) history row, source-table rows, VTK frame.

## 5. Mesh (`mesh.py`)

- **Active mask.** `VolumeMesh.active` (bool per tet, all True initially); `surface()`, `boundary_nodes()`,
  `element_volumes()` and the solver's patterns use the active set; deactivating elements exposes the faces of their
  inward neighbours, which become patches. `deactivate(elements)` returns the new patch ids and the map from dead
  patches to the patches that replace them (by nearest centroid), used by the film hand-over.
- **Prism layers.** `sphere_mesh(radius, h_surface, h_core, layers=n, layer_thickness=t0, growth=2.0)`: the sphere
  surface is meshed first, then extruded inward with gmsh's `extrudeBoundaryLayer` (n layers from t0 growing
  geometrically), the prisms split into tetrahedra, and the graded tets meshed inside the new inner surface; if the
  extrusion proves unreliable on this gmsh version, the fallback is a Distance/Threshold size field with a thin
  fine band (h = t0 within n·t0 of the surface) — same node-count scale, isotropic tets instead of prisms; the plan's
  first task is a spike settling which. Default for melting runs n = 4, t0 = 0.25 mm (≈ 55 k nodes on the 100 mm
  sphere), `--prism-layers 0` reverts to the Step 2 mesh; cached by all parameters.
- **Box mesh.** `box_mesh(lx, ly, lz, h)` for the Stefan test (insulated sides, one heated face).
- Melt-relevant node counts are reported in the run JSON; the layer-thickness convergence check is §13.5.

## 6. Materials and melting (`material.py`, `thermal/*`)

- `AA7075` = DRAMA's `drama-AA7075` verbatim (k, c_p tables to 850 K, melting temperature 850 K, L_f 400 kJ/kg,
  ε 0.40) with liquid properties added; `AA7075_range` = the same tables with solidus 750 K, liquidus 908 K, the
  latent heat spread linearly across the range. Both JSON files carry `liquid: {rho, mu, sigma}` and their sources.
- `Material.liquid_fraction(T)`: 0 below the solidus, 1 above the liquidus, linear between; for the single-temperature
  material a ramp of ±2 K around 850 K (a numerical smoothing of the step; the enthalpy jump is exact).
  `enthalpy(T)` = ∫c_p dT + L_f f_l(T) (the Step 2 hook), `temperature_from_enthalpy` its inverse on the augmented
  table, `cp_eff` = c_p + L_f df_l/dT.
- The thermal solvers already use the secant heat capacity [h(T_k) − h(T_old)]/(T_k − T_old), which absorbs the
  latent heat exactly; nothing changes in the Newton loop. Test: the Neumann (Stefan) similarity solution on the box
  mesh, front position x = 2λ√(αt) with λ from the transcendental equation, constant properties, within 1 %.
- Element fractions: the operators are assembled with φ_e ρ c and φ_e k per element (both scale with the material
  present: a thinner remaining sliver has less capacity; scaling k with φ_e keeps the surface node connected to the
  interior without letting an almost-consumed element conduct like a full one — a documented approximation).

## 7. Gas-side surface flow (`surface_flow.py`)

Per patch and step, from the Step 2 stagnation state (p_s, h_s, s_s, ρ_s) and the freestream:

- **Edge state.** p_e(θ) = p∞ + (p_s − p∞) cos²θ for θ ≤ 90° (modified Newtonian); isentropic expansion from the
  stagnation state: h_e = h(s_s, p_e), u_e = √(2(h_s − h_e)), ρ_e, T_e, a_e, M_e, μ_e (Blottner–Wilke), all from
  `gas.EquilibriumAir` (one Cantera state per θ bin of 1°, interpolated to the patches). Leeward (θ > 90°): p_e = p∞,
  u_e = 0, no shear.
- **Boundary layer.** Thwaites' method along the meridian s = Rθ: θ_m²(s) = 0.45 ν_e u_e⁻⁶ ∫₀ˢ u_e⁵ ds′ (momentum
  thickness), converted to Girin's linear-profile boundary-layer thickness δ_a = c·δ_99 with c fixed so that, under the
  potential-flow velocity 1.5 V sin θ and constant properties, δ_a reproduces Ranger's 2.2 R Re^(−1/2) Ψ(φ) (the test:
  5 % over 5°–85°). Wall shear in Girin's convention τ_c = μ_e u_e/δ_a.
- **Knudsen number and regime.** λ_e from ρ_e, T_e with the hard-sphere model of `aero`; Kn_δ = λ_e/δ_a:
  continuum Kn_δ < 0.01; slip 0.01 ≤ Kn_δ < 0.1; transitional/free-molecular Kn_δ ≥ 0.1.
- **Shear.** Continuum τ = τ_c; slip τ = τ_c/(1 + C Kn_δ), C = (2 − σ_v)/σ_v, σ_v = 1 (Maxwell first-order slip with
  the linear near-wall profile, implemented as an effective edge velocity u_e/(1 + C Kn_δ) that also enters Girin's
  Eq. 2); free-molecular τ_fm = σ_t ρ∞ V∞² sin θ cos θ, σ_t = 1 (hypersonic speed ratio, full accommodation).
  `--rarefied-shear slip` (default): τ = τ_slip below Kn_δ 0.1 and τ_fm above; `bridged`: τ = (1 − f) τ_slip + f τ_fm
  with the measured drag bridging f(Kn) of `aero.SesamTable` on the body Knudsen number, continuous across the
  regimes. The difference between the two is reported (§13.5).
- **Pressure gradient and body force.** G = −∂p_e/∂s + ρ_l a_t, with ∂p_e/∂s = −2(p_s − p∞) sin θ cos θ/R and a_t the
  tangential component of the body-frame deceleration (from the trajectory's a_drag).
- Outputs per patch: u_e, ρ_e, μ_e, M_e, δ_a, λ_e, Kn_δ, regime (0/1/2), τ, G — written to the surface VTK series.

## 8. Melt film (`film.py`)

- **State.** m_f per patch (kg), b = m_f/(ρ_l A_patch); carried across steps; hand-over on element death.
- **Feed (solidus–liquidus region, flagged for review §17).** Each surface element's liquid inventory is
  f_l(T̄_e) φ_e ρ V_e. An element feeds its patch(es) only once f_l ≥ f_coh (`--coherency`, default 0.5: the
  coherency point of the dendritic skeleton); the increment since the last step goes into the film. Once an element
  has lost coherency, its remaining solid fraction is treated as carried with the melt as a slurry at the liquid's
  viscosity; below f_coh the element counts as solid. For the single-temperature material f_coh is irrelevant.
  The liquid inventory of an element deeper than the surface element is counted when it becomes the surface.
- **Lubrication solution** with the patch's τ, G, b and δ_m (from §9), μ_l:
  thin branch (b ≤ δ_m): V_s = τb/μ_l + Gb²/(2μ_l), q = τb²/(2μ_l) + Gb³/(3μ_l);
  thick branch (b > δ_m): V_s = τδ_m/μ_l (Girin's Eq. 2 value), q = V_s δ_m/2 + Gb³/(3μ_l).
  The velocity gradient used by the instability is V_s/b (thin) or V_s/δ_m (thick).
- **Runoff transport (`--runoff on|off`, default on).** Explicit upwind finite volume on the patch graph: for each edge
  shared by patches i and j, flux = q_donor · ℓ_edge · (t̂_donor·n̂_edge)⁺, t̂ = unit surface direction of
  −∇p_e + ρ_l a_t (from the nose outward), donor = the upwind patch; sub-steps Δt_sub = 0.5 min(A_patch b/(q ℓ)) inside
  the macro step, vectorised over patches. Mass-conserving to round-off. The flow stops at the equator (G = 0 there);
  leeward films are static. Test: a strip of patches with uniform τ, G against the steady closed form, and global
  conservation on the sphere.

## 9. Spraying (`dispersion.py`, `spray.py`)

- **Dispersion relation.** For We_s on a log grid [3.0, 10⁴] (200 points): for each Δ on [0.02, 5] (2 500 points) the
  cubic in Ω from Eq. (1) is solved (`numpy.roots`), the root with the largest Im Ω is kept, and the Δ maximising
  Im Ω gives Δ_f, Im Ω_f, Re Ω_f; cached in `data/girin_dispersion.json` and interpolated. Tests: no unstable root at
  We_s ≤ 3.00 and one at 3.08; Δ_f(10⁴) = 1.225 ± 1 %, Im Ω_f(10⁴) = 0.24 ± 3 %; Re/Im within 1.5–1.7 for We_s ≥ 10.
- **Per patch, windward, m_f > 0**, after the film update. δ_m is computed from Girin's Eq. 2 for every continuum
  and slip patch (it decides the branch); rarefied patches have no gas boundary layer and are always thin-branch.
  1. *Continuum or slip, thick film (b > δ_m)*: δ_m = (α/μ²)^(1/3) δ_a, V_s = (αμ)^(1/3)/(1 + (αμ)^(1/3)) u_e,eff,
     We_s = ρ_l V_s² δ_m/Σ; unstable if We_s > We_cr (`--we-critical`, default 4.62); λ_f = 2πδ_m/Δ_f(We_s),
     r = k_r λ_f, t_per = k_t δ_m/(V_s Im Ω_f(We_s)); stripping rate per area ṁ = ρ_l π r²/(λ_f t_per) (one torus of
     cross-section r per wavelength per period; Girin's Eq. 7 integrated over the patch).
  2. *Continuum or slip, thin film (b ≤ δ_m)*: Girin & Kopyt 1994 side mode with d = b, M = M_e, We_d = ρ_e u_e² b/Σ:
     λ* = 1.5 M d/We_d, r = λ*/4, τ* per their Eq. (12), release ṁ = ρ_l b/τ* (the film sheds per growth time),
     limited by the film; the mode is active only when λ* ≥ λ_t = 0.5 M d/We_d and when the film thickness exceeds
     the wavelength cut-off's physical floor (d ≥ 1 µm, a numerical guard).
  3. *Transitional/free-molecular (Kn_δ ≥ 0.1)*: branch 2 with the free-molecular shear driving the film (Couette
     V_s = τ_fm b/μ_l) and M = M∞ — an extrapolation of a continuum film theory, flagged (§17), sensitive to
     `--rarefied-shear`.
  4. The 1994 front-surface Rayleigh–Taylor criterion W b² ρ_l > 3Σ (W = body deceleration) is evaluated per step and
     reported (`rt_active` column), expected inactive at ~30 m/s²; if active, its λ*, ω* are recorded, not applied.
- **Release bookkeeping.** Δm = min(ṁ, ṁ_melt + m_f/Δt)·A·Δt; Δn = Δm/(4/3 π ρ_l r³); one source row per emitting
  patch-step: t, h, V, θ, patch centroid (body frame), regime, b, δ_m, We_s, r, Δn, Δm, release velocity
  (body velocity + V_s t̂), We_d = ρ∞V²·2r/Σ and Oh = μ_l/√(ρ_l Σ 2r) with a flag when We_d > 12 (Pilch–Erdman).
  Histograms Δn(r), ΔM(r) on 40 log-spaced bins (1 µm–10 mm) per 10 s window and for the flight.
- **Constants** k_r = 0.17, k_t = 1.1, We_cr = 4.62 (`--kr`, `--kt`, `--we-critical`); Girin's calibration is on
  ordinary liquids, not aluminium (assumption §16).

## 10. Geometry and accounting (`body.MeltingBody`)

- Stripped mass and net runoff leaving a patch reduce φ_e of the patch's owner element (a patch has one owner; an
  element may own several patches — the reduction is shared by patch area). At φ_e < 10⁻³ the element dies (§5);
  the dead patches' film mass is handed to the replacing patches; `surface()` and the solvers' patterns are rebuilt.
- `mass(t)` = Σ_active φ_e ρ V_e + Σ m_f; the drag reference area = Σ A_patch max(0, n̂·v̂) over the current surface
  (πR² for the intact sphere); C_D from the sphere tables and the stagnation nose radius R₀ unchanged (assumptions
  §16, future work §17).
- Energy: E = Σ φ_e ρ V_e h(T̄_e) + Σ m_f h(T_liquidus); enthalpy leaving with stripped/run-off material accumulates in
  `removed_enthalpy_J`; the per-step balance ΔE = (Q_conv − Q_rad)Δt − ΔE_removed holds to the Newton tolerance
  (test 10⁻⁶).
- New history columns: `mass_kg` (varying), `film_mass_kg`, `sprayed_mass_kg` (cumulative), `runoff_mass_kg`
  (cumulative net transport), `melt_front_depth_max_mm`, `equivalent_radius_mm` (from the remaining mass; compared
  with SESAM's `thick_mm`), `n_active_elements`, `spraying_area_m2`, `theta_cr_deg`, `n_released`,
  `released_mass_kg`, `r_median_um`, `r_max_um`, `regime_fraction_continuum`, `_slip`, `_fm`, `rt_active`,
  `removed_enthalpy_J`.
- End of run: ground, escape, or demise when the remaining mass < `--demise-fraction` (default 0.01) of the initial;
  `end_reason` = "demise".

## 11. Coupling and time stepping (`coupled.py`)

- Macro step Δt = 0.5 s (`--dt`); the melt step follows the thermal step (§4 order); the mass and projected area it
  produces enter the next trajectory advance (first-order splitting, as the heating).
- Spraying is a rate over the step: Δm = min(ṁ_strip, ṁ_melt + m_f/Δt) A Δt — stripping-limited (thick branch) or
  melt-limited (thin branch: the film sits at a quasi-steady sub-millimetre thickness); Δt-independent in both limits,
  the Δt-halving check covers the crossover.
- Runoff sub-stepped by its CFL; element death at most once per element per step, after the film update; the
  solver operators use the new φ_e every iterate (patterns rebuilt only on death, ~20 ms).
- Heating on a molten patch unchanged; the hot-wall factor uses the film temperature (≈ liquidus); no vaporisation
  or blowing (aluminium boils at 2743 K; the film stays near 900 K).

## 12. Outputs and visualization

- History CSV: Step 2 columns + §10 columns. Run JSON: melt onset (time/altitude of the first film), spraying onset,
  demise, sprayed/run-off/film/remaining masses, size statistics, regime fractions, sensitivity flags used.
- `<run>/particles.npz` (column arrays of the source table), `particles_summary.csv` (per step),
  `size_distribution.csv` (Δn, ΔM per bin per window and total), written with `--particles` (default on when melting).
- VTK: surface series adds b, We_s, regime, τ, r, release rate; volume series adds the liquid fraction.
- Plots (saved in `<run>/`): mass vs time and vs altitude, cumulative sprayed/run-off/film/remaining, θ_cr and
  spraying area vs time, regime fractions, r_median and We_s range vs time, Δn(r)/ΔM(r) per window and final
  (Girin Fig. 6 style); in verification mode each with a SESAM overlay and residual panel.
- Videos: the surface video with emitting patches overlaid (coloured by r, the θ_cr ring), a film-thickness video, the
  cross-section video with the liquidus iso-line; stills at melt onset, spraying onset, peak release rate,
  demise/impact.

## 13. Verification and acceptance

Every comparison writes an overlay + residual plot and a metrics JSON; reference data live in
`data/reference_runs/` (SESAM) and `data/reference_values/` (published tables, analytic parameters). Thresholds are
expectations in the Step 1/2 sense: measured values are recorded in the README, a threshold changes only with the
measured value and the reason beside it.

1. **SESAM melting references.** Two new wrapper runs with DRAMA's `drama-AA7075` on the US76 table, winds off:
   `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nowind` and
   `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nowind`, committed with their JSONs.
   - *Bookkeeping check (thresholded)*: SESAM-equivalent heating + `AA7075` + `--removal instant` + `--runoff off` +
     k × 10⁴ (near-isothermal body): `mass_kg(t)` within 2 % of the initial mass at every reference time, melt-onset
     and demise altitudes within 0.5 km, demise time within 2 %.
   - *Resolved run (reported)*: the same with the real conductivity: onset earlier than SESAM's (surface melts first),
     demise time expected within ~5 % (energy-limited); the difference is the lumped-body assumption, plotted and
     discussed in the README.
2. **Girin 2017 Table 1** (`analysis/girin_reference.py`, "Girin-as-published" mode: his constant ρ_a, potential-flow
   u_e, W(τ) law, no runoff, iron/stone properties, We_cr 4.62): variants I and II — N, r_med, size range, t_s.d., σ
   within 20 %; variant III reported; `data/reference_values/girin2017_table1.json`; Δn(r), ΔM(r) plotted with his
   table values marked.
3. **Girin & Kopyt 1994 Tables 1–2**: r_d and τ_d for the six (ρ₂, V₀) cases and RT λ*, τ* for W = 10⁵–10⁷ cm/s²
   within 5 % (transcription check; their ṁ is reproduced with the mass-rate definition transcribed from the paper
   and reconciled with §9's ρ_l b/τ* — the reconciliation is recorded, and §9's rate is adjusted if the paper's
   definition differs), `girin1994_tables.json`, log-log plots with residuals.
4. **Analytic and conformance**: Stefan/Neumann front on the box mesh 1 %; steady runoff on a strip 0.1 % and exact
   conservation; Thwaites vs Ranger 5 %; dispersion-relation values (§9); energy and mass balance with melting and
   removal 10⁻⁶; element death keeps the surface closed and the area/volume consistent; both thermal backends give
   the same melting run (0.1 %).
5. **Convergence and sensitivity (reported)**: prism layers 2/4/8 (t0 0.25 mm), Δt 0.5/0.25 s, `--rarefied-shear`,
   `--runoff`, We_cr 3.08/4.62, k_r and k_t ±30 %, coherency 0.3/0.5/0.7, backend skfem/fenicsx — a table of sprayed
   mass, r_median, melt-onset, spraying-onset and demise altitudes for the 100 mm and 50 mm physics-mode flights
   (`analysis/melt_sensitivity.py`).
6. **Cost target**: a 100 mm physics-mode melting flight on the default (4-layer) mesh in ≤ 6 min serial.

## 14. CLI

`run` adds `--melt off|on` (default off: Step 2 behaviour unchanged), `--material AA7075|AA7075_range|<path>`
(default `AA7075_range` with `--melt on`), `--removal girin|instant`, `--runoff on|off`, `--rarefied-shear
slip|bridged`, `--we-critical`, `--kr`, `--kt`, `--coherency`, `--prism-layers`, `--layer-thickness` (mm),
`--demise-fraction`, `--particles/--no-particles`. `compare` accepts the new columns and produces the mass-loss plots
when both files carry `mass_kg` variation. Exit codes as before. `analysis/girin_reference.py`,
`analysis/melt_verification.py`, `analysis/melt_sensitivity.py`.

## 15. Environment

No new packages: numpy, scipy (`roots`, `brentq`, `cumulative_trapezoid`), gmsh (boundary-layer field), Cantera,
scikit-fem/pyamg, PyVista as in Step 2; `fenicsx_env` for the backend cross-check.

## 16. Assumptions to state in the thesis

1. Fixed attitude; axisymmetric heating and shear (the patch model is 3D but the loads are).
2. Girin's gradient-instability theory with its boundary-layer relations and constants (k_r 0.17, k_t 1.1,
   We_cr 4.62) calibrated on ordinary liquids; the thin-film branch from an inviscid analysis.
3. Liquid AA7075 as pure aluminium near the liquidus (ρ_l, μ_l, Σ); no oxide skin (Denis et al. saw it stripped).
4. Solidus–liquidus region: coherency threshold f_coh = 0.5, slurry at the liquid's viscosity (§17).
5. Continuum boundary-layer relations down to Kn_δ 0.1 with first-order slip; the free-molecular branch extrapolates
   the thin-film theory (§17).
6. Modified-Newtonian pressure with isentropic expansion for the edge state; laminar boundary layer (Thwaites);
   leeward films static; runoff stops at the equator.
7. Nose radius R₀ and the sphere drag tables kept for the eroded body; projected area from the current surface.
8. No vaporisation, blowing, oxidation or shock-layer radiation; heating unchanged over the film.
9. One representative droplet radius per patch and step; droplets recorded at birth.
10. First-order operator splitting at 0.5 s; element death in element-size quanta (0.25 mm layers by default).

## 17. Future-iteration review items (recorded here so they are not lost)

1. **Within-patch size spread**: a torus does not break into identical droplets and the unstable band contains a
   range of wavelengths; Girin's 1990 drop-shattering result N(r) ∝ r⁻ˡ, l ≈ 5.5–7.5, could be draped over each
   release, normalised to Δm — as a post-processing option on the source table (each row carries Δn, Δm, r) or a
   switch in `spray`. Left out on 2026-09-20 because it adds an empirical layer with no aluminium data behind it.
2. **Solidus–liquidus rheology**: the coherency convention of §8 ignores the strength of the mushy skeleton and the
   viscosity rising by orders of magnitude toward the solidus; it moves the melt-onset altitude by up to ~1 km.
3. **Rarefied thin-film branch**: no film-instability theory exists for free-molecular driving; the extrapolation and
   the shear bridging are the least certain parts of the spraying onset for small spheres.
4. **Droplet fate** (Step 4): secondary breakup (We_d expected 1–5 here, checked per row), cooling and
   solidification, wake transport.
5. **Shape feedback**: nose radius and drag of the eroded body; local recession finer than the layer quantum (ALE).
6. **Girin's constants for aluminium**: k_r, k_t and We_cr are the sensitivity study's main knobs; an aluminium
   drop-shattering or arc-jet measurement (Kim et al. 2020's AA7075 melt-removal rate under a torch is the closest
   data) would calibrate them.
