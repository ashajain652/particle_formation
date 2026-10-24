## Physics model — `reentry_model` (Step 3: melting, melt film and melt spraying)

Step 3 turns the heated sphere into a particle source. On top of the Step 2 coupling, every macro step now also
melts the body (enthalpy method with the latent heat in the nodal enthalpy, element fractions φ_e and element death
on a mesh with four 0.25/0.5/1/2 mm prism layers under the surface), feeds the liquid to a film on the surface
patches, moves the film with a lubrication runoff driven by the gas shear and the pressure gradient, and strips it
into droplets by whichever instability Girin's own regime classification selects: Girin & Kopyt's (1994) thin-film
mode where the melt is too shallow for his predicted melt boundary layer δ_m to form (regime 1), classical
Kelvin–Helmholtz where it has formed and the melt is the more viscous medium (regime 2), and his (2017) gradient
instability where it has formed and the gas is (regime 3) — with a rarefied variant of regime 1 on freestream
conditions wherever the body-scale gate finds no distinct bow shock, and Girin & Kopyt's front-surface
Rayleigh–Taylor mode on the stagnation cap. Droplets are recorded at birth in a source
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
`--runoff on|off`; `--rarefied-shear slip|bridged`; `--we-critical 4.62`, `--kr 0.17`, `--kt 1.1`;
`--rt-spray on|off` (default on: apply Girin & Kopyt's front-surface Rayleigh–Taylor mode on the patches where it
passes both criteria and outgrows the shear mode; `off` reports it only, as every run before 2026-09-24 did);
`--prism-layers 4`,
`--layer-thickness 0.25` (mm, growth 2; 46 k nodes / 255 k tets on the 100 mm sphere); `--demise-fraction 0.01`;
`--particles/--no-particles`; `--size-feedback current|initial` (`current` with `girin`: the body Knudsen number on
the equivalent diameter of the remaining mass and the stagnation radius fitted to the windward cap, bounded to
1.67 × the transverse radius for a flat front; `initial` with `instant`: D₀ and R₀, SESAM's convention); `--k-scale`
(a verification device). One melt setting has no command-line flag: `MeltSettings.feed_depth`, default `conjugate`,
which lets melt leave an element only where the gas shear can reach it — `all` restores the pre-2026-09-22 behaviour
and is reachable from Python only. Melt and runoff start at the liquidus: material between the solidus and the
liquidus holds its latent heat but counts as solid for the film (spec §8, §17.2). The film has the surface's own temperature (it is
thermally thin), carries the liquid enthalpy wherever it sits, and re-solidifies onto its owner element when the
surface falls back through the ramp — the mirror of the feed rule, netted with it so that no element both melts and
freezes in one step. Film that has nowhere to freeze (a full or dead owner element: the mesh cannot grow a crust
outside itself) stays liquid, and `film_frozen_fraction` reports how much — at most 0.47 % of the initial mass on the
50 mm flight and 0.126 % on the 100 mm one. Re-solidification is a leeward, end-of-flight phenomenon: on the 50 mm
flight the film sprays away before it can cool through the ramp, so only 7.2 µg freezes back — and 99.9 % of that is
leeward, where a patch gets no convective heat and loses heat by radiation and into the cold rear (the leeward film
runs 10–25 K colder than the windward film). The 100 mm remnant, which survives to the ground, freezes 2.4 g back.
Droplets leave at the film's own temperature, so they carry its superheat: +2.4 % on h_liquid for the 50 mm flight
(the film reaching 1013 K against a 908 K liquidus), +0.6 % for the 100 mm one.

Outputs: `<run>.csv` gains `mass_kg` (now varying) and the 48 melt columns of `coupled.MELT_COLUMNS`: the mass
budget (`film_mass_kg`, `sprayed_mass_kg`, `runoff_mass_kg`, `removed_mass_kg`, `released_mass_kg`, `frozen_mass_kg`,
`removed_enthalpy_J`); the shape the trajectory now feels (`melt_front_depth_max_mm`, `equivalent_radius_mm`,
`nose_radius_mm`, `transverse_radius_mm`, `fitted_nose_radius_mm`, `drag_shape_factor`, `n_active_elements`,
`n_dead_elements`); the release (`spraying_area_m2`, `theta_cr_deg`, `n_released`, `r_median_um`, `r_max_um`,
`spray_growth_ms`); the film (`film_thickness_max_mm`, `film_thickness_mean_mm`, `film_T_max_K`, `film_T_mean_K`,
`film_frozen_fraction`, `film_blob_fraction`, `unapplied_load_J`); the liquid layer and the regime it selects
(`molten_depth_max_mm`, `molten_depth_mean_mm`, `delta_m_mean_um`, `thick_branch_fraction`); the flow-regime gate and
the melt closure it certifies (`kn_body`, `kn_local_stag`, `re_shock`, `flow_branch`, `p_w_stag_Pa`, `phi_sonic_deg`,
`closure_fraction_girin`, `closure_fraction_couette_slip`, `closure_fraction_couette_fm`); and the Rayleigh–Taylor
diagnostics (`rt_active`, `rt_mass_fraction`, `rt_region_mm`, `rt_wavelength_over_nose`, `rt_growth_ms`,
`rt_bounded_fraction`, `rt_bounded_growth_ms`). The run
JSON adds the melt onset, spraying onset, demise, the masses, size statistics and the settings; `<run>/particles.npz`
(the source table: time, altitude, velocity, θ, patch centroid, the flow closure, the spray regime, film thickness,
δ_m, We_s, radius, count, mass, release velocity and direction, We_d, Oh, breakup flag), `particles_summary.csv`,
`size_distribution.csv`
(Δn, ΔM on 40 log bins 1 µm–10 mm per 10 s window and for the flight); plots `mass_time`, `mass_altitude`,
`mass_budget`, `spraying_time`, `closures_time`, `droplet_size_time`, `size_distribution` (with the SESAM overlay
and residual when a melting reference is given); videos `animation.mp4` (surface temperature with the emitting
patches coloured by droplet radius), `film.mp4` (film thickness), `section.mp4` (cross-section with the liquidus and
solidus iso-lines) and stills at melt onset, spraying onset and peak release in addition to Step 2's. The VTK series
add the liquid fraction and φ_e (volume, active elements only) and the film thickness, film temperature, We_s, the flow
closure, the wall Knudsen number, the wall pressure, the shear, the droplet radius and the release rate (surface). A
macro step costs about 2 s on the default mesh (0.6 s of it conduction, the rest the melt step's Cantera edge states, the
runoff solve and the contiguous-layer walk); a flight that demises takes 2–5 min, while one whose remnant survives runs
to the ground and takes about 40 min (1191 steps for the 100 mm physics case).

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
| d100 | bookkeeping | 0.98 % | 71.10 / 71.00 | 65.9 / 66.8 (−1.3 %) | — | — | 39 s, 132 steps |
| d050 | bookkeeping | 1.29 % | 77.27 / 77.10 | 190.9 / 191.4 (−0.2 %) | — | — | 20 s, 382 steps |
| d100 | resolved (sesam heating, AA7075, girin, D₀/R₀) | 9.35 % | 71.62 / 71.00 | 70.0 / 66.8 (+4.8 %) | 1.404 / 0.057 | 2.7e7 (160 µm) | 269 s, 141 steps |
| d050 | resolved | 7.10 % | 77.57 / 77.10 | 191.3 / 191.4 (−0.1 %) | 0.161 / 0.021 | 4.5e5 (298 µm) | 54 s, 383 steps |
| d100 | physics (AA7075_range, girin, shape feedback) | 76.4 %† | 73.96 / 71.00† | **no demise** | 0.974 / —† | 1.6e7 (188.6 µm) | 2357 s, 1191 steps† |
| d050 | physics | 66.5 %† | 78.32 / 77.10† | 207.5 (demise) / 191.4 | 0.1792 / 0.0026 | 3.1e6 (224.3 µm) | 130 s, 415 steps |

The two physics rows carry the measurements of 2026-09-25 except where marked †, which are the 2026-09-21/22 values that
the gates of amendments 23–26 have not re-measured; Task 14's `summary.md` replaces the whole table. The d050 physics row
is a fresh run of the committed model, and its end-of-life cell is the **demise time**, not the interpolated 1 %-mass
time, which needs the reference comparison to compute. The four device rows are unaffected by any of it, being pinned to
`--size-feedback initial`.

Thresholds (bookkeeping mode only, `tests/test_reentry_model_reference_melt.py`): mass 2 % of m₀, onset 0.5 km,
1 %-mass time 2 %. The devices are pinned to `--size-feedback initial`, so none of the shape or regime amendments can
move them — which is what makes them a fixed yardstick. The resolved runs are reported: the surface melts 0.6 km before
SESAM's lumped body reaches 850 K, and the interior's sensible heating delays the end by up to 4.8 %; the 7–9 % mass
difference is the lumped-body assumption, plotted in `d100__resolved/mass_time.png`. Giving the film its own latent
heat (2026-09-22) moved both cases toward SESAM — the 50 mm mass error halved, from 13.5 % to 7.1 %, and its 1 %-mass
time from +0.5 % to −0.1 % — because melt can no longer be relabelled as film without paying for itself.

**The physics-mode 100 mm sphere does not demise.** It melts from 74.0 km, sprays 0.974 kg of its 1.472 kg as about
1.6×10⁷ droplets, and the remaining **0.498 kg (33.8 % of the initial mass) reaches the ground**. That is a consequence of the
fixed-attitude assumption acting three times over: the flattening nose is held face-on, which is the maximum-drag
orientation (C_D rises from 0.91 toward 1.8, so the body decelerates high and the heating ∝ ρV³ collapses); the
stagnation radius grows as the face flattens, cutting the stagnation flux by a further ~20 %; and the leeward shell is
never heated at all. SESAM's lumped sphere, which keeps D₀, R₀ and the sphere drag table, demises at 66.5 km. A
tumbling fragment would sit between the two — DRAMA's own tumbling-averaged C_D for a thin disc is 0.60, *below* the
sphere's 0.91 — so the sphere/face-on spread is an attitude uncertainty, not a drag-law one, and tumbling is the first
item of the next iteration (spec §17). The 50 mm sphere still demises, at 207.5 s and 68.08 km against SESAM's 191.6 s and 73.1 km.

**Next iteration** (in order): tumbling, replacing the fixed-attitude assumption; an Arbitrary Lagrangian–Eulerian
mesh in which each patch recedes every step according to its own mass loss, which dissolves the element-removal
question rather than improving it — no threshold, no residual mass to reassign, no hand-over rule, and a shape that
updates continuously instead of in jumps; the pressure-driven part of the runoff flux integrated over the whole liquid
depth while the shear-driven part stays at the conjugate depth, which is how the molten region's response to pressure
differences enters; a patch graph that merges sliver facets, or a film carried per element rather than per facet, which
is the real answer to the single-patch concentration; and reporting the runoff/stripping coupling parameter (melting
speed × conjugate depth / kinematic viscosity) so that the decoupled treatment is checked rather than assumed.

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

Sensitivity (`analysis/melt_sensitivity.py`, 100 mm physics flight, one setting changed per row): fifteen variants —
prism layers 2/4/6, Δt/2, `--rarefied-shear bridged`, no runoff, We_cr 3.08, k_r and k_t ±30 %, the single-temperature
material, `--size-feedback initial`, `--gamma-pm 1.4`, `--kn-body-shock 0.003` and the FEniCSx backend — with the table
written to `reentry_model_output/verification_melt/sensitivity/sensitivity.md`. The two amendments of 2026-09-22 get a
row each so their uncertainty is bounded rather than asserted: γ_PM moves the wall pressure beyond the sonic point by a
factor ~2, and the body gate moves the fraction of the flight that is Girin-certified.

Findings recorded while building this step (details in the spec's amendments):
- The film is stripped as fast as it melts: every instability branch removes hundreds to thousands of kg/m²/s where
  the melt supply is ~5 kg/m²/s, so the mass loss is energy-limited and the film stays microns thin (Girin's
  "outstripping ablation" case).
- Where Girin's closure may be used is a measurable result, not an assumption. On the 100 mm physics flight melting
  starts at 73.9 km inside the *merged* branch (Kn_body 0.0172, Re₂ 102); the gate opens at 69.0 km and stays open to
  the end, so 69 % of the melting steps and 84 % of the sprayed mass carry Girin's Eq. (2) closure and the rest is
  flagged Couette. At the 77.5 km break-off every sweep diameter is merged (Kn_body 0.60 at 5 mm to 0.030 at 100 mm).
- The gate is deliberately conservative. Evaluated at the wall — cold and compressed — the local Knudsen number is
  λ_w/λ∞ ≈ 1/167 of the freestream value, i.e. Kn_local ≈ Kn_body/84, so even the 5 mm sphere at 77.5 km has a nose
  Kn_local of 0.0066 while its body value is 0.60. We decline to certify those patches because the construction that
  would check them is the one we do not trust there.
- Pure modified Newtonian was wrong where it mattered most: it gives p_w = p∞ at 90°, a factor 28–48 below the
  Prandtl-Meyer value (5.2 Pa against 144–249 Pa at 100 mm/70 km) and below the measured sphere C_p band of 112–219 Pa.
  It was also the cause of the "rarefied rim" in the earlier design — with the expansion resolved, the rim carries
  ρ_e ~ 1e-4 kg/m³ and Kn_local ~ 2e-3, so the film no longer piles up there and needs no droplet-size cap.
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
