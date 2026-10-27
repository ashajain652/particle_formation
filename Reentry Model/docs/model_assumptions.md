# `reentry_model` — assumptions and methods (Steps 1–2)

State of the model on `main` as of 2026-09-20 (Step 1: trajectory, spec `superpowers/specs/2026-09-17-reentry-trajectory-model-design.md`;
Step 2: coupled 3D heat transfer, spec `superpowers/specs/2026-09-18-thermal-fem-design.md`). Items marked **(verified)** were
measured against SESAM or an analytic solution; the numbers are in the README's verification tables. The SESAM facts the
model relies on are in `Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md`.

The Spheral-free core of the large-fragment model (`spheral_frag`, milestone M1) has its own assumptions document,
`spheral_frag_assumptions.md`: the frame contract as confirmed on the 2026-10-08 Scheil flight, the material table, the
thickness, depths and zones, the load tables and their measured gap to the history's drag, and the fragment record's
definitions.

The wake column (`reentry_model/wake_column/`, the chemistry and condensation of the aluminium vapour a droplet sheds) is
section 9. There **(verified)** names the pytest test or the measurement that supports the entry, not a comparison with SESAM.

## 1. Trajectory (Step 1)

Assumptions
- Rigid, non-rotating solid sphere of constant mass; point-mass 3-DOF dynamics (no lift, no side force, C_L = 0).
- Drag from the freestream only; drag is attitude-independent (sphere). Attitude matters only for heating (below).
- WGS-84 Earth, rotating frame (Coriolis + centrifugal), gravity with J2 (J2+J4 and point-mass optional).
- Velocity is relative to the rotating atmosphere; winds either none or DRAMA's static profile.
- The flight ends at the ground (0 m) or on escape above 150 km; no fragmentation.

Methods
- State in ECEF Cartesian coordinates, integrated with DOP853 (rtol 1e-9, atol 1e-6 m / 1e-9 m/s); ground and escape
  as terminal events.
- Kn = λ/D with a hard-sphere mean free path (d = 3.65e-10 m) and the local mean molecular mass; Ma with γ = 1.4.
- C_D from DRAMA's ATDB_SPHERE tables (C_D,fm(Ma), C_D,c(Ma), Ma 5–30, clamped; below Ma 5 SESAM's rules: the continuum
  value held to Ma 1, halved below Ma 1), bridged by f(Kn) measured from SESAM's own output (`aero.SesamTable`,
  0.25-decade bins of log10 Kn; erf and sin² alternatives).
- Atmosphere: NRLMSISE-00 via pymsis with the DRAMA fap-file solar indices; DRAMA's US76 table; or a replay of a SESAM
  run's own density/temperature (diagnostic).
- **(verified)** vs SESAM on the same US76 table: 3.9 m/s / 18 m (100 mm) and 3.6 m/s / 4 m (50 mm) over the hypersonic
  phase. NRLMSISE-00 mode differs by ~6 % because SESAM's built-in NRLMSISE-00 has about half the seasonal amplitude of
  the reference implementation (unresolved on SESAM's side).

## 2. Aerothermal heating (Step 2)

Common assumptions
- Fixed attitude: the stagnation point stays at the same material point (+x of the body frame); θ per surface patch is the
  angle between the patch's outward normal and the direction of motion.
- Convective heating only: no shock-layer radiation, no oxidation heat, no ablation blowing (negligible for R ≤ 0.1 m
  below 9 km/s).
- Freestream composition N₂:O₂ = 0.79:0.21 by mole at all altitudes; T∞ clamped to ≥ 200 K (NASA-9 polynomial floor).

Physics mode (`--heating physics`, the model proper)
- Continuum stagnation flux: Fay–Riddell (1958), equilibrium boundary layer, Pr = 0.71, Le = 1.4, Newtonian velocity
  gradient du_e/dx = (1/R) √(2 (p_s − p∞)/ρ_s); wall enthalpy at the patch temperature; the catalycity γ_cat scales the
  Lewis-number term (1 = fully catalytic, SESAM's assumption; 0 leaves the Le = 1 value — the frozen non-catalytic
  reduction 1 − h_D/h_s is not modelled).
- Stagnation state: a normal shock (Rankine–Hugoniot, damped fixed point) followed by an isentropic compression to rest,
  both in **equilibrium air** (Cantera `airNASA9.yaml`, NASA-9 polynomials 200–20 000 K, 11 species including ions)
  whenever the estimated stagnation temperature exceeds 1500 K, frozen composition otherwise; for Ma ≤ 1.1 no shock is
  taken (isentropic compression straight from the freestream: the fixed point stalls at Mach 1 and the entropy jump of
  such a shock is < 0.1 %).
- Viscosity: Blottner curve fits (N₂, O₂, NO, N, O; Blottner, Johnson & Ellis 1971 as tabulated in Gnoffo, Gupta & Shinn,
  NASA TP-2867) with Wilke's mixing rule over the equilibrium mole fractions; ions neglected in the mixture viscosity
  (`airNASA9` carries no transport data).
- Cheap alternatives: Sutton–Graves (1.7415e-4 √(ρ/R) V³) and Detra–Kemp–Riddell, both with the hot-wall factor
  (h_s − h_w)/(h_s − h_w,300 K).
- Free-molecular stagnation flux: q_fm = A_cq ρ V (h_s − h_w) with A_cq = 0.8 (Matting 1971, Eq. 1).
- Bridging: Matting (1971) Eq. 18, q = q_c P(n, [Γ(n+1) q_fm/q_c]^{1/n}) with P the regularised lower incomplete gamma
  function; n = 1 by default (q = q_c [1 − exp(−q_fm/q_c)]); SESAM's drag f(Kn) as an alternative.
- Distribution over the sphere: (1 − w)·Lees (1956) laminar distribution with finite-Mach ε = 1/(γ Ma²), zero leeward,
  + w·cos θ (free-molecular shape), with w = 1 − q_stag/q_c under Matting (f(Kn) under the SESAM table); each patch
  scaled by its own hot-wall factor relative to the stagnation patch.
- **(verified)** Fay–Riddell at the 100 mm start (ρ 2.727e-5 kg/m³, 7.5 km/s): 2.22 MW/m² = 1.29 × Sutton–Graves —
  72 % of the stagnation enthalpy is dissociation at 1.5 kPa, and Sutton–Graves is a Le = 1 fit; the Lees distribution
  integrates to 0.196 (Ma → ∞) / 0.209 (Ma 5) of the stagnation flux over the whole sphere; physics mode delivers
  0.74–0.77 × SESAM's integrated heat.

SESAM-equivalent mode (`--heating sesam`) — a verification device, not a physical model
- SESAM's measured heat input applied **uniformly to every patch, front and back**:
  q = 0.27471 · q_DKR · F_h(Kn) · max(0, 1 − c_p (T_lumped − T∞) / (V²/2)) for Ma ≥ 1 (c_p = 1004.5 J/kg·K), and half
  of 0.27471 · q_DKR · F_h with no hot-wall term below Ma 1. 0.27471 is the ATDB shape factor (§8 below); F_h(Kn) is
  SESAM's heat bridging measured from the two no-melt references in 0.125-decade bins of log10 Kn (`aero.SesamHeatTable`:
  1.005 in the continuum, 0.14 at Kn 1, 0.059 at Kn 40). The hot-wall factor uses the body's energy-equivalent (lumped)
  temperature, as SESAM's lumped model does.
- Findings that shaped it (facts note §15): SESAM's heating carries this hot-wall factor (invisible below 850 K, where
  the Step 1 facts were measured); its transitional heating is not its drag blend; its free-molecular limit is 0.78 × the
  textbook cos θ average of ½ρV³ (the earlier "13× too low q_FM" was the transitional deficit misread at Kn 0.03).
- Purpose: give the finite-element body exactly SESAM's heat input so that the conduction, time stepping, material
  handling and coupling can be compared with SESAM's lumped temperature with no distribution question in between.

## 3. Radiation and material

- Grey-body radiation ε σ (T⁴ − T_amb⁴) from every patch, ε = 0.40 (DRAMA's AA7075), T_amb = 0 K to match SESAM
  (200 K optional); no view factors, no atmospheric absorption or emission.
- Material `AA7075_nomelt` (a copy of `data/user_materials/AA7075_nomelt.json`): DRAMA's drama-AA7075 ρ = 2813 kg/m³,
  k(T) and c_p(T) tables (293–850 K, held at their 850 K values above); the specific enthalpy h(T) is the exact integral
  of the piecewise-linear c_p; no phase change in this step (latent-heat hook present, liquid fraction ≡ 0).
- Uniform initial temperature T₀ (300 K by default).

## 4. Finite-element conduction

- Geometry: gmsh OCC sphere meshed with linear tetrahedra graded from h_surface = 2 mm at the surface to h_core = 8 mm
  at the centre (18.9 k nodes, 18 078 surface patches on the 100 mm sphere); a node embedded at the centre; the boundary
  patches derived from the tetrahedra themselves; node coordinates kept as a mutable array (Step 3 recession).
- Weak form on P1 tetrahedra: ∫ ρ c_p Ṫ v + ∫ k ∇T·∇v = ∫_Γ (q_conv − ε σ (T⁴ − T_amb⁴)) v; consistent mass by default
  (row-sum lumping optional).
- Time integration: backward Euler. Per Newton iterate: k(T) at the element-mean temperature of the previous iterate;
  mass with the **secant heat capacity** [h(T_k) − h(T_old)] / (T_k − T_old), which makes the discrete energy balance
  ΔH = (Q_conv − Q_rad) Δt exact for any h(T) (and is the hook for Step 3's latent heat); the radiation term linearised
  from the facet-mean temperature with its exact Jacobian 4 ε σ T_f³ A_f/9. The systems are symmetric positive
  definite; tolerance 1e-6 relative; the Newton start is extrapolated from the previous step (2.0–2.1 iterations per
  step); non-convergence raises an error rather than returning an unconverged field.
- Linear algebra (scikit-fem backend): unit element stiffness and mass matrices precomputed once and rescaled into a
  fixed CSR pattern every iterate; CG preconditioned by a pyamg smoothed-aggregation hierarchy reused for 30 solves
  (SuperLU optional). FEniCSx backend: the same linearised system in UFL/dolfinx 0.11, in-place assembly, PETSc CG +
  hypre BoomerAMG reused for 30 solves; serial.
- Energy-equivalent mean temperature T_eq from H(T_eq) = (1/V) ∫ H dV — the quantity compared with SESAM's lumped
  temperature; also surface max/min/area-mean, stagnation, back and centre temperatures.
- **(verified)** lumped limit vs the lumped ODE 8e-5; radiative cooling of an isothermal sphere 1.7e-4;
  Carslaw–Jaeger step change: centre 0.5 % / volume mean 0.2 % of the rise (uniform 4 mm test mesh); energy balance
  1e-10 per step with constant c_p and ~1e-8 over full flights with the tabulated c_p; scikit-fem vs FEniCSx
  cross-check 0.1 % and identical coupled-flight results.

## 5. Coupling

- First-order operator splitting with a 0.5 s macro step: the trajectory advances with a fresh DOP853 solve → aero state
  at the end of the step → heating from that state and the wall temperatures at the start of the step → implicit
  thermal step. Mass is constant, so the thermal state does not feed back into the trajectory in this step; the loop
  already carries the body's mass and the mesh for Step 3.
- History every macro step (Step 1 columns plus convective/radiated/integrated/absorbed heat, T_eq, surface and centre
  temperatures, q_stag, the bridging weight); VTK series of the volume and surface fields; surface and cross-section
  animations rendered from the exported files.
- **(verified)** the coupled trajectory reproduces Step 1's `run()` to < 0.01 m/s and < 0.1 m over the whole flight;
  halving h_surface changes the surface-temperature history by ≤ 0.03 % (SESAM-equivalent mode) and ≤ 0.11 % (physics
  mode, stagnation temperature included); halving Δt changes the peak surface temperature by 0.04 %.

## 6. Verification against SESAM — what it shows and what it does not

- SESAM-equivalent mode on the two no-melt US76 references (100 mm from 77.5 km, 50 mm from 115 km): Q_conv within
  0.5 % of peak (2.2–2.4 % point-wise in the continuum), integrated heat +0.31 % / +0.03 %, T_eq within 1.2 %
  (24–28 K), radiated power within 6.7 % / 3.4 % of peak (the resolved surface radiates at its own, hotter temperature).
- Caveat: F_h(Kn) and the hot-wall c_p were measured on those same two references, so the heat rows measure the fit's
  residual (2.6 % / 1.4 % point-wise when the formula is applied to SESAM's own state columns); the independent content
  of the verification is the temperature and the radiated power — i.e. the conduction, time stepping, material handling
  and coupling.
- Physics mode is reported as a ratio to SESAM (0.744 / 0.770), not thresholded: SESAM assumes random tumbling and a
  lumped body; the model assumes a fixed attitude and resolves a ~600 K stagnation-to-mean temperature difference at
  peak heating.

## 7. Simplifications to keep in mind (thesis)

Fixed attitude (no tumbling); laminar, equilibrium boundary layer with catalycity as a parameter; constant freestream
composition; grey emissivity constant with temperature and no oxide layer; no melting or mass loss yet; first-order
splitting at 0.5 s; P1 elements with 2 mm surface resolution (Step 3's melt layer needs the prism-layer mesh);
SESAM's own transitional heating law is unknown and reproduced only as a measured table; the FEniCSx backend is serial
(MPI not exercised).

## 8. The ATDB shape factor

DRAMA's aerothermal database `ATDB_SPHERE.nc` (HTG, 2019) tabulates, for the sphere, two heat factors next to the drag
coefficients: a continuum factor of 0.27471 (constant over Ma 5–30) and a free-molecular factor (0.2707 at Ma 5 →
0.2501 at Ma 30). A heat factor is the ratio of the surface-averaged convective flux to the stagnation-point flux, so
that SESAM's lumped model takes the total heat as Q = A_sphere · factor · q_stag without resolving a distribution.

Attitude does not enter for a sphere: whichever way it faces it presents the same geometry to the flow, so the
instantaneous flux averaged over the whole surface is the same number, tumbling or fixed. (HTG builds its databases
under a random-tumbling convention, which matters for plates, boxes and cylinders, not for spheres.) SESAM with a fixed
attitude would therefore still apply 0.27471 to the sphere — the factor is HTG's surface-average-to-stagnation ratio
for the sphere's flux distribution, not a tumbling correction.

For comparison, Lees' laminar distribution over the windward hemisphere gives 0.196 (Ma → ∞) / 0.209 (Ma 5) of q_stag
over the whole area and nothing on the lee side, and the free-molecular cos θ distribution gives 0.25. HTG's value is
37–40 % above the windward-only laminar integral, so its sphere distribution carries leeward/base heating or a broader
windward distribution than Lees'; the HTG SARA modelling report that would say which is not in hand (measured base
heating on spheres is only a few percent of q_stag, which does not obviously close the gap). That distribution
difference — not attitude — is the main reason the physics mode delivers ~0.74 × SESAM's integrated heat.

## 9. Wake column (`reentry_model/wake_column/`)

Added 2026-10-27 (design: `superpowers/specs/2026-10-27-wake-column-chemistry-design.md`; plan:
`superpowers/plans/2026-10-27-wake-column-chemistry.md`, whose decisions 1-19 are numbered below as "D1" to "D19").
The model follows the aluminium a demising droplet vaporises: Stage A slows the hot vapour to ambient speed (sets the
starting composition and radius), Stage B is a 12-hour box model of a diffusing horizontal cylinder with the gas
chemistry of Plane et al. (2021), sectional particles tallied by source molecule and route, uptake onto background
meteoric smoke, coagulation and settling. Every entry is marked **(verified)** with the test (in
`reentry_model/wake_column/tests/`) or measurement that supports it, or **(assumed)**. The controller's rulings made
during the build (marked "ruling") are binding amendments to the plan and are recorded where they apply. The results
quoted here are the build's diagnostics at default settings; the default sweep (126 runs) and its plots are the
research output and are not summarised here.

### 9.1 Geometry, spreading and ambient air

- **A horizontal cylinder along the trail.** The droplet's flight path angle is about -1 degree, so its vapour is a
  near-horizontal line source. Box length L is the path flown from release to the moment 90 % of the aluminium loss has
  occurred; the inventory is all the aluminium vaporised, deposited at t = 0 (deposition takes 10-260 ms); the altitude
  is the aluminium-mass-weighted mean over the whole loss (D16). **(verified)** `test_trail_length_matches_design_measurements`,
  `test_synthetic_linear_loss_gives_exact_length_altitude_and_samples`, `test_inventory_is_the_whole_aluminium_loss`.
  L is interpolated inside the droplet history's step (ruling): 1 micrometre at 80 km gives 44.6 m and 10 micrometres at
  60 km gives 38.2 m (tests pin them to 2 %), against the design's 51 m and 40 m, which were the path at the first
  history row past 90 % and so biased long by up to one coarse row (about 13 % for the 1 micrometre case; a shorter L
  raises the starting concentrations by that factor). The depositing of the whole inventory at once **(assumed)** is
  valid except against Al + O2 (milliseconds); a segmented trail is future work.
- **Extensive state (D1).** The ODE state is total counts in the plume, so aluminium conservation is a plain sum, and
  entrainment adds ambient air (`c_ambient dV/dt`) while the aluminium-bearing species are diluted only through V.
  **(verified)** `test_aluminium_conserved_in_full_run` (worst departure of gas + particles + smoke + settled from the
  inventory 1.1e-14 by day and 5.2e-15 by night, 100 trillion atoms at 80 km for 1 hour, tolerance 1e-6),
  `test_gas_aluminium_is_conserved_with_full_chemistry`, `test_short_run_file_contents` (1e-5 at every history row),
  `test_inert_tracer_follows_r0_over_r_squared`, `test_ambient_stays_ambient_without_aluminium`.
- **Box volume (D2).** V = pi r^2 L / erf(L / 2r): the centreline concentration of a uniform line segment broadened by
  a Gaussian of mean-square radius r^2 = r0^2 + 4 integral of D dt. It tends to pi r^2 L for r much smaller than L (V
  growing as t) and to pi^1.5 r^3 for r much larger (V as t^1.5). **(verified)** limits, finiteness over L/r from 1e-6
  to 1e6, the analytic derivative against a finite difference and continuity across its series switch:
  `test_volume_limits`, `test_volume_finite_and_monotone_over_extreme_ratios`, `test_dlnV_dr2_matches_a_finite_difference_of_ln_volume`,
  `test_dlnV_dr2_continuous_across_the_series_switch`. **(assumed)** the centreline form rather than a mean-square
  volume, which would be twice as large and halve second-order rates; the spec fixed V = pi r^2 L.
- **Radius growth.** d(r^2)/dt = 4 (D_mol + K_eddy S(r)); D_mol is AlO in air; K_eddy is the Megner et al. (2006)
  reference eddy diffusion profile (Fig. 3c) digitised to `data/wake_column/megner2006_eddy_diffusion.csv`.
  **(verified)** `test_eddy_profile_is_megner_reference` (peak 9.7 m^2/s at 90 km, 0.30 m^2/s at 50 km),
  `test_total_diffusion_adds_the_scaled_switched_eddy_part`. The digitisation is by pixel trace of a published figure,
  **(assumed)** good to about 10 % in K. Particles diffuse with the gas **(assumed)**: Brownian diffusion of particles is
  slower than the gas's, so this is an upper bound on their dilution.
- **Eddy switch (D3).** S(r) is 0 below l_inner/sqrt(2), 1 above sqrt(2) l_inner and a cubic smoothstep in ln r
  between; l_inner = 30 m is **(assumed)** (sensitivity 10 and 100 m) because eddies smaller than the plume cannot
  mix it. **(verified)** `test_switch_limits`, `test_switch_is_smooth_monotone_and_centred_in_ln_r`.
- **Molecular diffusion (D4, D5).** Hard-sphere first-order Chapman-Enskog with the collision integral set to 1 (no
  Lennard-Jones well depths exist for aluminium molecules) and air as a pseudo-species (3.711 angstrom, 28.96 amu).
  Collision diameters are Svehla (1962) for N2, O2, O, H, H2, H2O, CO, CO2 and Al; every one-aluminium molecule is 3.5
  angstrom and O3 3.6 angstrom **(assumed)**. **(verified)** formula on a hand state, `test_molecular_diffusion_formula_on_a_hand_state`;
  AlO at night is 0.0219, 0.270 and 5.90 m^2/s at 50, 70 and 90 km, within the spec's windows
  (`test_molecular_diffusion_matches_spec_magnitudes`). The Svehla values are as cited in the plan, not re-read from the report.
- **Ambient air.** Baseline: the WACCM full-chemistry starting file at 40.7 N, day at 180 E and night at 0 E,
  extracted to `data/wake_column/ambient_waccm_40.7N.csv` on a 1 km grid (45-100 km) with geometric altitude (ruling:
  hypsometric geopotential height converted to geometric; the file has no surface-height field, so terrain is taken as
  zero, a few hundred metres at most). **(verified)** `test_columns_cover_45_to_100_km_and_close`,
  `test_native_column_geometric_altitude_densities_and_n2_remainder`, `test_to_grid_interpolates_temperature_linearly_and_densities_log_linearly`,
  `test_outside_table_raises_with_altitude_in_message` (an altitude outside the table raises, never extrapolates),
  `test_day_night_atomic_oxygen_contrast_below_80_km`. At 71 km WACCM gives O 3.63e15 per cubic metre by day and 4.06e11 by night.
  The plume's ambient comes from WACCM rather than the droplet code's US76 atmosphere, matched by geometric altitude **(assumed)**.
  Held fixed for 12 hours (a day run holds the noon column, a night run the midnight one) **(assumed)**; one date (1
  January 2001, 00 UT) and one latitude only.
- **Cross-check ambient.** NRLMSIS 2.x (pymsis, F10.7 = 150, Ap = 4) for T, n, N2, O2, O, H plus AFGL 1986 midlatitude
  winter for O3, H2O (CO2 is the US Standard profile, 330 ppmv, against WACCM's about 370 ppmv, a 10 % systematic).
  **(verified)** `test_waccm_and_msis_agree_on_total_density_within_factor_2`, `test_no_jump_between_adjacent_levels`
  (no species jumps by 0.7 decades between 1 km levels in either table; the old splice had 9.3). Findings the
  cross-check inherits: MSIS supplies no O below 51 km and no H below 76 km, so those are the WACCM profile scaled to
  MSIS at the lowest valid level (the two sources are not independent there); MSIS night O at 71 km is 1.56e14, about
  380 times WACCM's 4.06e11 (day/night ratio about 60 in MSIS, about 9000 in WACCM), an empirical extrapolation with no
  nocturnal recombination, so the cross-check has far more night oxygen below about 80 km. MSIS closes to 0.991-0.992
  (helium and argon are not network species).

### 9.2 Stage A: the hot start

- **Slowing-down Monte Carlo (spec 4.4, D15).** Hard-sphere test particles (default 10,000 atoms over 20 emission
  samples) leave the droplet with an isotropic speed drawn from the flux-weighted Maxwellian at the droplet temperature
  plus the droplet speed along the track; because the air is isotropic only the droplet's speed matters. Emission
  samples sit at the midpoints of equal-aluminium-mass intervals of the loss history. Scattering is isotropic in the
  centre-of-mass frame with no inelastic loss (this makes the stopping distance slightly long), products keep the
  reactant's velocity, exothermicity and reaction heat are neglected, an atom is thermal below 2 x (3/2) k_B T, and r0
  is the root-mean-square perpendicular stopping distance **(all assumed)**. **(verified)**
  `test_mean_energy_loss_per_collision_hard_sphere` (0.4997 against 2mM/(m+M)^2 = 0.4998 for Al on N2),
  `test_stopping_distance_scales_inversely_with_density` (r0 ratio 80 km to 70 km 4.92 against a density ratio 4.88),
  `test_cold_start_has_no_threshold_reactions` (r0 = 1.12 mean free paths), `test_packet_conservation_and_hot_products`,
  `test_partner_mix_follows_density_times_cross_section` (N2 0.804 and O2 0.196 against 0.807 and 0.193),
  `test_degenerate_inputs_stay_finite` (Review Focus 4). The null-collision acceptance uses the candidate partner's own
  majorant term (ruling); the plan's literal reading undercounted O2 collisions about four times. The flux-weighted
  emission speed (a Gamma(2) draw) is exact by construction but has no separate distribution test.
- **Line-of-centers reactions.** P(E) = P0 (1 - E0/E) above the threshold, calibrated to the thermal rate expressions
  at the top of each validity range, capped at 1, constant for barrierless Al + O2; AlO + M dissociation has P0 = 1 and
  E0 = 498 kJ/mol as an upper bound **(assumed)**. **(verified)** `test_line_of_centers_reproduces_thermal_rate` (0.997-0.999 of
  the analytic rate against 5 % tolerance), `test_calibrated_probabilities` (R1 0.882, AlO + O2 -> OAlO + O 0.0818, D1
  1.0; OAlO + O 0.531 and AlO + O -> Al + O2 0.00996 were checked by hand, not pinned by a test).
- **Result at default settings, 14.1 micrometres, WACCM day (diagnostic).** At 60, 76 and 90 km the fractions Al/AlO/OAlO are
  0.488/0.512/0.0003, 0.365/0.631/0.0034 and 0.370/0.624/0.0062, and r0 is 1.66 mm, 21.6 mm and 181 mm. The design's
  "half or more is AlO" held; its "a few percent OAlO at the highest speeds" did not: OAlO is under one percent, about ten
  times smaller than the design estimated (AlO + O2 -> OAlO + O reacts in only 8 % of above-threshold collisions, and
  about one collision in five is with O2). The OAlO yields rest on 30-60 test particles, so carry about 15 % Monte
  Carlo uncertainty; raise `n_test_particles` to 100,000 if results prove sensitive to it.
- **Ionisation finding (ruling).** The spec's statement that 578 kJ/mol (the aluminium ionisation energy) is "not
  reached by Al + air at 7.6 km/s" holds for the directed beam only: the thermal emission velocity of D15 adds a tail.
  The share of aluminium collisions above 578 kJ/mol is 0, 1.3e-4 and 1.2e-3 at 60, 76 and 90 km (the share of AlO
  collisions above the AlO bond is 0, 3.7e-3 and 9.7e-3), and 0.0166 for an 1800 K droplet at 7.6 km/s (against 0.0
  for a 300 K droplet). The count is of collisions with enough energy, a necessary condition, not of ionisations, so
  ion chemistry stays a non-goal, but about 1 % of top-speed collisions could ionise and are not modelled. The test
  bound is therefore 0.05 at 1800 K plus 1e-3 at 300 K; the approved spec was not edited. **(verified)** the numbers,
  `test_packet_conservation_and_hot_products`.
- **`hot_start = "droplet_temperature"` (Asha's original picture).** The same Stage A with the droplet speed set to zero
  (the reactions stay on; atoms leave with thermal speed only, so threshold reactions are essentially closed and r0 is
  a few mean free paths). An earlier comment in `settings.py` that called this "no hot chemistry" was wrong (ruling,
  corrected). **(verified)** that it runs through `run.run` and gives zero collisions above the AlO bond:
  `test_sensitivity_runs_and_tornado_plot`.

### 9.3 Gas chemistry

- **Network.** The neutral reactions R1-R17 of Plane et al. (2021) plus additions from Gobrecht et al. (2022) Table B.1
  (reactions 26, 33, 36, 38 and 42) and the AlO bond dissociation as a hot-stage-only threshold row. All data are in
  `data/wake_column/reactions.csv` and `species.csv` with a source and a validity range per row. **(verified)**
  `test_plane_table_values_reproduced` (Plane's tabulated rates at 200 and 300 K), `test_species_table_has_every_network_species`.
- **Third bodies (D12).** [M] is the total air number density with Plane's N2-measured coefficients for all air; O2 is
  assumed as efficient a third body as N2 **(assumed)**. **(verified)** `test_AlO_lifetime_against_R6`: the AlO lifetime
  in the box equals 1/(k6 n_total [O2]), which is 0.093, 2.27 and 45.9 s at 70, 80 and 90 km in the night WACCM columns
  (the spec's "about 0.08, 2, 40 s" used round numbers); `test_Al_lifetime_against_R1`.
- **Duplicate removed (D13).** The spec's "AlO + O2 -> AlO2 + O and its reverse": the reverse, OAlO + O -> AlO + O2, is
  Plane's R3, so only R3 is kept. **(verified)** `test_reverse_of_G42_is_not_duplicated`,
  `test_hot_only_dissociation_not_in_box_network`.
- **Validity ranges (D14).** Where the source gives none: Plane's (T/300)^(1/6) entries 100-1000 K, R1 23-295 K (Le
  Picard 1997, Gomez Martin 2017), R6 and R7 100-600 K, Gobrecht/Mangan additions 100-2500 K **(assumed;
  `validity_source = assumed` in the table)**. A rate outside its range is clipped and warned about once per reaction per
  run. **(verified)** `test_clipping_warns_once_per_reaction` (Review Focus 5). The Gobrecht reaction 32 (Al + O2) has
  (T/300)^-0.26; Plane's (T/300)^(1/6) form is used.
- **Species data (D5, D6).** Dipole moments: AlO 4.45 D (Bai & Steimle 2020, as used by Plane), AlOH 0.97 D (Sikorska &
  Skurski 2009); the others are unknown and their pairs use the hard-sphere collision rate **(assumed, flagged)**.
  Multiplicities for odd-electron species are assumed doublet. **(verified)** `test_unknown_dipole_falls_back_to_hard_sphere`,
  `test_zero_dipole_counts_as_unknown_for_the_dipole_term`.
- **Photolysis.** AlOH photolysis (J = 3.3e-3 per second, R17, temperature independent) is on in day runs and off in
  night runs **(assumed)**: the `day_only` flag is in the reaction table and conservation is tested by day and night
  (`test_gas_aluminium_is_conserved_with_full_chemistry`), but no test isolates R17's on-off switch. The day-night
  difference in the AlOH share at 12 hours (gas-only, 76 km: AlOH 0.645 and Al(OH)2 0.344 by day against 0.717 and
  0.283 by night, from the box's gas-phase build run) has not been traced to photolysis rather than the ambient O and
  H contrast.
- **Heat.** Reaction heat and the heat deposited by slowing atoms are neglected: the box is at the ambient temperature
  from t = 0 for 12 hours **(assumed)**, because the aluminium is a trace once r exceeds a few r0.

### 9.4 Particles

- **Bins by aluminium atoms (D8).** Pivots n_k = 2 q^k, from the two-atom cluster (0.476 nm alumina sphere, bin 0) until
  the alumina diameter passes 400 nm; number and atoms conserved by fixed-pivot splitting (Kumar & Ramkrishna 1996);
  above the last pivot an entity goes to the last bin with its number scaled to conserve atoms. Default q is 2^(1/2) (60
  bins), the plan's fallback: q = 2 gives a self-preserving geometric standard deviation of 1.506 by numerical
  diffusion, q = 2^(1/2) gives 1.4755 and q = 2^(1/3) gives 1.4684, against the literature's about 1.46 and the test's
  band of 1.45 +- 0.03 (extrapolating to q = 1 gives about 1.454). **(verified)** `test_coagulation_reaches_self_preserving_width`
  (margin to the 1.48 bound is only 0.0045), `test_split_conserves_number_and_atoms`,
  `test_above_the_last_pivot_the_entity_goes_to_the_last_bin_with_number_scaled`, `test_place_conserves_number_atoms_inside_the_range`,
  `test_bin_count_and_last_pivot_at_ratio_two` (31 pivots, last at 2^31 atoms = 487 nm at q = 2),
  `test_first_bin_is_a_half_nanometre_alumina_dimer`. Every deposit goes through `Bins.place` (ruling), so no atoms are
  lost in the largest bin. Species with two or more aluminium atoms are first-bin clusters **(assumed)**.
- **Atoms delivered per molecule (D9).** (Al, O, H): Al (1,0,0), AlO (1,1,0), OAlO (1,2,0), OAlO2 (1,3,0), AlOH
  (1,1,1), AlOH2 (1,2,2), AlCO3 (1,1,0) (CO2 assumed to leave on sticking). The O released by AlO + AlO -> Al2O + O is
  not tracked, and particle oxidation does not consume gas-phase O, O2 or O3 (the aluminium is a trace) **(all assumed)**.
  **(verified)** `test_seeded_bins_carry_one_oxygen_per_aluminium`.
- **Size and density (D10).** Volume per Al atom = (1 - f) 2.8218e-29 m^3 + f 5.3524e-29 m^3 with f = min(H/Al / 3, 1)
  (bulk alumina 3.0 g/cm^3, aluminium hydroxide 2.42 g/cm^3); the tracked O content is a diagnostic that changes no size
  or mass. **(verified)** `test_hydroxide_is_bulkier_per_aluminium`, `test_diameter_matches_the_sphere_of_the_composition_volume`.
  Bulk densities for nanometre clusters, compact spheres (no aggregates) and no dehydration of hydroxide **(assumed)**.
  The H/Al reporting threshold of 0.5 for "alumina-like" is a reporting choice, not physics.
- **Formation (D7).** One labelled pathway for each of the 28 pairs of aluminium-bearing species (none lumped), at
  max(k_DD, k_HS) x spin x efficiency x `formation_scale`: k_DD = pi Gamma(2/3) v (2 C6 / k_B T)^(1/3) with the
  Keesom-averaged C6, k_HS = pi sigma^2 v, spin = g_P / (g_A g_B). **(verified)** `test_dipole_capture_AlO_AlO_at_298K`
  (1.516e-9 cm^3/s with 4.36 D, 1.558e-9 with the table's 4.45 D), `test_spin_factors`, `test_every_pair_is_its_own_pathway`,
  `test_capture_rate_is_symmetric_and_never_below_hard_sphere`, `test_every_pathway_and_the_smoke_route_are_reported`.
  The 1.5e-9 cm^3/s against Gomez Martin et al. (2017)'s 1.1e-9 is about 40 % higher **(assumed)**; every efficiency
  defaults to 1 and is a sensitivity parameter, so the formation rates are upper bounds with unmeasured efficiency.
- **Growth.** Any aluminium-bearing molecule sticks to a particle of any bin at the free-molecule collision rate
  including its own size times `sticking[species]` x `sticking_scale` (default 1) **(assumed efficiency)**. Smoke uptake
  uses the infinite-mass limit for the smoke particle, slightly overestimating the rate.
- **Oxygen uptake (D11, ruling).** Rate per bin = sum over O, O2 and O3 of collision rate x `oxidation_eff` x g(O/Al),
  one O atom per event; g is 1 below O/Al = 1.485, 0 at 1.5 and above and a cubic smoothstep (3s^2 - 2s^3) between, not
  the plan's linear ramp: the linear ramp's slope jump (-67 to 0) stalled the BDF Newton iteration and let O/Al drift
  to 2-3. The ramp is a numerical cap, not data; O has its own absolute tolerance (1e-18 scaled). Hydrogen is gained only
  through condensing hydroxide molecules (no gas-surface H or H2O reaction) **(assumed)**. **(verified)**
  `test_oxygen_uptake_factor_runs_from_one_to_zero_with_a_continuous_slope`,
  `test_seeded_particles_take_up_oxygen_to_one_and_a_half_per_aluminium` (within 1e-3 in 1 s); with the smoothstep the
  O/Al of every bin holding more than 1e-9 of the inventory agrees with a tight (relative tolerance 1e-10) solution to
  1e-2 in the 76 km day and night, 60 km night and 90 km day cases.
- **Background meteoric smoke.** Number density and size distribution from Megner et al. (2006) Fig. 2, digitised to
  `data/wake_column/megner2006_smoke.csv`. Fig. 2b gives the size distribution at 90 km only and shows only about
  16 % of Fig. 2a's number (the smallest smoke is below the plotted axis), so (ruling) the deficit is put in one bin at
  Megner's smallest radius, 0.2 nm, and other altitudes use the 90 km shape scaled by Fig. 2a **(assumed)**. That bin
  holds 4.0e4 of 4.8e4 per cubic centimetre at 90 km and about 47 % of the smoke surface area, so smoke uptake depends
  on it at about a factor of two. The smoke is a fixed background, not depleted, not coagulated with plume particles,
  and sticks at efficiency `smoke_sticking` (default 1) **(assumed)**. **(verified)** `test_smoke_table_integrates_to_fig2a_number_density`
  (1 % at 60, 75, 90 km), `test_smoke_uptake_equals_the_direct_sum_at_a_table_altitude`, `test_smoke_uptake_outside_the_table_raises_naming_altitude_and_file`.
  Uptake rate of AlO at 200 K is 3.2e-7, 3.7e-6 and 8.8e-6 per second per molecule at 60, 75 and 90 km (a time scale of
  about 32 hours at 90 km).
- **Coagulation.** Free-molecule Brownian kernel between all bins (half for a self-pair), particles compact spheres
  **(assumed)**. **(verified)** `test_coagulation_reaches_self_preserving_width`.
- **Settling.** Bin k leaves at 2 v_s / (pi r), v_s the Epstein (free-molecule) speed, through the lower half of a
  well-mixed cylinder **(assumed)**; settled atoms keep their route and drop their O and H (Al-only settled tally).
  **(verified)** `test_settling_only_decays_at_two_vs_over_pi_r` (worst relative error 1.6e-5 against 1e-4).
- **Irreversible condensation** (Al-O clusters are bound by 400-500 kJ/mol per atom, and the plume is at ambient
  temperature within a fraction of a millisecond) and **particles diluted at the gas's rate** are **(assumed)**.
- **Provenance closure.** The tallies by source and route sum to the condensed aluminium. **(verified)**
  `test_provenance_closes_in_every_bin` (3e-12 atoms in every bin; the 21 provenance shares sum to particles + smoke +
  settled to 9e-19).

### 9.5 Numerics, outputs and definitions

- **Solver (D19).** `solve_ivp` BDF, rtol 1e-6 and atol 1e-12 on scaled variables: aluminium counts divided by the
  inventory, particle number of bin k by inventory / n_k (so the tolerance bounds atoms), each ambient count by
  c_ambient V0, r^2 by r0^2. The Jacobian is analytic for every term except the r^2 column (a one-sided difference).
  **(verified)** `test_analytic_jacobian_matches_finite_differences` (every row outside the ambient gas agrees to 1e-4 of
  the row maximum; the kernels' weak dependence on H/Al is left out, which affects only Newton convergence); condensed
  fractions at 76 km agree to 4-5 significant figures with a solve at relative tolerance 1e-10. A 12-hour run takes
  about 6-25 s depending on the case (measured with seven runs in parallel on an eight-core machine).
- **"Condensed" (D17)** is everything that left the gas phase: plume particles + background-smoke uptake + settled.
  t50 and t90 are the first output times at which it reaches 0.5 and 0.9, or null. **(verified)** `test_short_run_file_contents`.
- **Steady state (D18).** The steady time is the earliest output time t_ss such that, for every t at or after it with
  10 t at most t_end, every aluminium reservoir changes by less than 0.01 of the inventory between t and 10 t; at
  least one decade of evidence is required (null otherwise). Reported, not used to stop. **(verified)**
  `test_steady_state_time`, `test_steady_state_needs_a_decade_of_evidence`.
- **`d_volume_mean_nm` (ruling)** is the volume-weighted mean diameter, the sum of N d^4 over the sum of N d^3 over the
  plume particles at the end of the run (where the condensed volume sits); the alternative, the diameter of the
  mean-volume particle (sum N d^3 / sum N)^(1/3), is not reported. `d_number_mean_nm` is the plain number mean.
  Smoke-taken aluminium is not in either.
- **Output masking (ruling).** In `size_distribution.csv` a bin holding under 1e-9 of the inventory in aluminium has
  empty composition columns (its O/Al and H/Al are set by a handful of atoms) and its alumina diameter; N and the tallies
  are kept. **(verified)** `test_short_run_file_contents`.
- **Skipped runs.** A droplet that vaporises no aluminium gives a summary with `skipped: "no aluminium vaporised"`
  and no outcome keys, and no CSVs; the plot scripts leave such runs out. **(verified)** `test_zero_inventory_is_skipped_cleanly`,
  `test_plot_scripts_ignore_a_skipped_run`. A truncated or unreadable `summary.json` is rerun, not fatal
  (`test_unreadable_summary_is_rerun`).
- **Run names** encode the case, the ambient and every non-default setting, so a sweep skips a finished run
  (`test_every_setting_changes_the_run_name`, `test_existing_summary_is_reused`). The test-only `trail=` argument of
  `run.run` is not part of the name.
- **Test-only switches.** `Settings.particles=False` also turns off smoke uptake (a particle-phase term); the other
  switches (`chemistry`, `coagulation`, `settling`, `smoke`, `freeze_radius`, `reactions`) are in the run name.

### 9.6 Sweeps and the sensitivity set

- **Default sweep.** 21 release altitudes (90 to 50 km in 2 km steps) x diameters 1, 14.1 and 200 micrometres x day and
  night = 126 runs, and the same on the NRLMSIS + AFGL ambient **(assumed extent; the design's call)**.
  **(verified)** `test_default_sweep_is_126_cases`.
- **Sensitivity set (ruling).** 14.1 micrometres at 60, 76 and 90 km by day and night; 15 one-at-a-time settings (the
  design's 13 plus "smoke sticking x0.1" and "no smoke uptake", which test the assumed smoke bin): formation x0.1, x0.01
  and x80 (Plane's polymerisation enhancement), sticking x0.1, oxidation x0.1, hot reactions x0.1 and x3,
  droplet-temperature start, eddy x0.1, x100 and zero (molecular only), inner scale 10 and 100 m, smoke sticking x0.1
  and 0. 15 x 6 = 90 runs; the baselines come from the default sweep. **(verified)**
  `test_sensitivity_set_is_15_settings_over_6_cases` (each setting changes exactly one field). The plot scripts pick
  baseline runs by their physics settings, so shortened test runs and cross-check runs sharing a folder are not mixed
  (`test_baseline_signature_ignores_numerics_but_not_physics`).

### 9.7 Findings from the build, to carry into the research

- **At 76 km almost nothing nucleates in the plume.** A 14.1 micrometre droplet's 2.9e13 aluminium atoms spread over a
  366 m trail with r0 = 21.6 mm: molecular diffusion (0.73 m^2/s) takes the plume to metres within a second, a
  7000-fold dilution in volume, and formation, being second order, stops at once. Plume particles hold about 3.1e-5 of the
  inventory and never grow past the first bin (0.485 nm by day, 0.482 nm by night); the condensed fraction at 12 hours is
  0.168 by day and 0.167 by night, nearly all of it first-order uptake onto background smoke (AlOH 0.114 and Al(OH)2
  0.052 by day), about 82 % is still gas (AlOH 0.54 and Al(OH)2 0.29), and t50, t90 and the steady time are null.
  This rests on the assumed 0.2 nm smoke bin, a smoke sticking of 1, the diffusion coefficient of D4, and one altitude;
  lower altitudes, where r0 is smaller (1.7 mm at 60 km), may nucleate more. **(verified)** as a measurement
  (`run.run` on 14.1 micrometres, 76 km, day, WACCM, 6.6 s); the physical conclusion is **(assumed)** until the
  sensitivity set has been run.
