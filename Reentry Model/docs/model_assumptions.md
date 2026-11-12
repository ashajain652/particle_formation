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

## 7a. The shared-package refactor changed no model behaviour

On 2026-11-11 the modules that both model architectures use moved into `reentry_model/shared/`, with alias files at
the old paths (see `CLAUDE.md`). No physics, default, option or output changed. This was shown by running a Step 1
flight, a Step 2 flight, an element-model melting flight, a layered-skin flight and a resumed skin flight before and
after the move with `analysis/skyscrapers_port/same_results.py` and comparing each case's files, which were identical,
and by the unit and reference test tiers passing with the same results as before. The same was done on 2026-11-11 for
the move of the skyscrapers-only modules into `reentry_model/skyscrapers/` and of the film's lubrication formulas into
`reentry_model/shared/lubrication.py`: the five cases were identical.

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
  **(verified)** `test_aluminium_conserved_in_full_run` (the test pins gas + particles + smoke + settled to the
  inventory within 1e-6 at every output, 100 trillion atoms at 80 km for 1 hour, day and night; the worst departure
  measured once in the Task 8 build was 1.1e-14 by day and 5.2e-15 by night),
  `test_gas_aluminium_is_conserved_with_full_chemistry`, `test_short_run_file_contents` (1e-5 at every history row),
  `test_inert_tracer_follows_r0_over_r_squared`, `test_ambient_stays_ambient_without_aluminium`.
- **Box volume (D2).** V = pi r^2 L / erf(L / 2r): the centreline concentration of a uniform line segment broadened by
  a Gaussian of mean-square radius r^2 = r0^2 + 4 integral of D dt. It tends to pi r^2 L for r much smaller than L (V
  growing as t) and to pi^1.5 r^3 for r much larger (V as t^1.5). **(verified)** limits, finiteness over L/r from 1e-6
  to 1e6, the analytic derivative against a finite difference and continuity across its series switch:
  `test_volume_limits`, `test_volume_finite_and_monotone_over_extreme_ratios`, `test_dlnV_dr2_matches_a_finite_difference_of_ln_volume`,
  `test_dlnV_dr2_continuous_across_the_series_switch`. **(assumed)** the centreline form rather than a mean-square
  volume, which would be twice as large and halve second-order rates; the spec fixed V = pi r^2 L.
- **Radius growth.** d(r^2)/dt = 4 (D_mol + min(K_eddy, K_R(r))); D_mol is AlO in air; K_eddy is the Megner et al.
  (2006) reference eddy diffusion profile (Fig. 3c) digitised to `data/wake_column/megner2006_eddy_diffusion.csv`, which
  acts as the cap on the scale-dependent K_R below. **(verified)** `test_total_diffusion_grows_with_size_up_to_the_capped_eddy_part`,
  `test_box_spreads_with_the_altitude_dissipation_rate`; the table itself is only range-checked by
  `test_eddy_profile_is_megner_reference` (peak between 5 and 15 m^2/s, minimum above 0.05 over 50-97.5 km; the
  values 9.7 m^2/s at 90 km and 0.30 m^2/s at 50 km were read off the table once in the Task 3 build, not pinned).
  The digitisation is by pixel trace of a published figure, **(assumed)** good to about 10 % in K and not
  independently verified. Particles diffuse with the gas **(assumed)**: Brownian diffusion of particles is
  slower than the gas's, so this is an upper bound on their dilution.
- **Scale-dependent eddy diffusion (2026-11-05, replaces the 30 m switch S(r), decision D3).** The review
  `Literature Review/Alumina Nanoparticles/reports/Smallest eddies in reentry wakes.md` found that a fixed 30 m gate was
  too large below about 75 km and that a gate is the wrong shape: eddy mixing grows continuously with plume size. So
  K_R(r) = g^(1/3) 3^(2/3) / 4 epsilon^(1/3) r^(4/3), the K for which d(r^2)/dt = 4 K reproduces Richardson's law
  <l^2> = g epsilon t^3, using r^2 = <l^2> / 3 for a cloud's cross-section. **(verified)** the algebra,
  `test_richardson_K_reproduces_the_t_cubed_law`. **(assumed)** g = 0.52 (`richardson_g`), from DNS, with one experiment
  near 0.1 (K_R goes as g^(1/3), so 0.58 times); this mapping gives K_R = 0.42 epsilon^(1/3) r^(4/3), weaker than
  the review's beta 0.4-0.6 with l = sqrt(2) r (0.63-0.95 on r^(4/3)), which would not reproduce the DNS g.
  **(assumed)** the law is used at all plume sizes: below the Kolmogorov scale it overstates the true (viscous) mixing,
  but there it is already far below D_mol, since turbulence overtakes molecular diffusion above the Kolmogorov scale at
  every altitude. **(assumed)** a sharp min() cap at Megner's K_eddy: the review's mismatch between K_eddy and the measured
  dissipation rates (K_eddy implies epsilon 2-3 times too high at 70 km and 2-10 times too low at 90 km with N = 0.02
  per second) is kept, with Megner's K read as an occurrence-weighted mean and epsilon as the value in turbulence.
  Turbulence is intermittent (for example seen 16 % of the time at 90-100 km in winter) and the law treats it as
  continuous **(assumed)**; `epsilon_scale` x0.1 and x10 are in the sensitivity set.
- **Dissipation rate.** `data/wake_column/mesospheric_dissipation.csv`, central values 1, 2, 1, 10 and 30 mW/kg at 50, 60,
  70, 80 and 90 km, log-linear between rows. **(verified)** the table is read as tabulated,
  `test_dissipation_table_covers_the_ambient_columns`. **(assumed)** applying it at 40.7 N: every measurement is from
  65-69 N, in winter for most of these values, and nothing measured exists at 50-55 km; the 60 km value is one radar
  upper limit.
- **Molecular diffusion (D4, D5).** Hard-sphere first-order Chapman-Enskog with the collision integral set to 1 (no
  Lennard-Jones well depths exist for aluminium molecules) and air as a pseudo-species (3.711 angstrom, 28.96 amu).
  Collision diameters are Svehla (1962) for N2, O2, O, H, H2, H2O, CO, CO2 and Al; every one-aluminium molecule is 3.5
  angstrom and O3 3.6 angstrom **(assumed)**. **(verified)** formula on a hand state, `test_molecular_diffusion_formula_on_a_hand_state`; the spec's
  magnitudes (windows 0.003-0.03 m^2/s at 50 km and 1-10 m^2/s at 90 km, `test_molecular_diffusion_matches_spec_magnitudes`),
  for which AlO at night measured once in the Task 3 build is 0.0219, 0.270 and 5.90 m^2/s at 50, 70 and 90 km. The Svehla values are as cited in the plan, not re-read from the report.
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
  `test_mean_energy_loss_per_collision_hard_sphere` (within 0.01 of 2mM/(m+M)^2 = 0.4998 for Al on N2; measured 0.4997),
  `test_stopping_distance_scales_inversely_with_density` (within 15 % of the density ratio; measured r0 ratio 80 km to
  70 km 4.92 against a density ratio 4.88), `test_cold_start_has_no_threshold_reactions` (no OAlO, no AlO collision
  above the bond, r0 between 0.1 and 10 mean free paths; measured 1.12), `test_packet_conservation_and_hot_products`
  (the fractions are exactly Al, AlO and OAlO and sum to 1),
  `test_partner_mix_follows_density_times_cross_section` (within 3 %; measured N2 0.804 and O2 0.196 against 0.807 and 0.193),
  `test_degenerate_inputs_stay_finite` (Review Focus 4). The null-collision acceptance uses the candidate partner's own
  majorant term (ruling); the plan's literal reading undercounted O2 collisions about four times. The flux-weighted
  emission speed (a Gamma(2) draw) is exact by construction but has no separate distribution test.
- **Line-of-centers reactions.** P(E) = P0 (1 - E0/E) above the threshold, calibrated to the thermal rate expressions
  at the top of each validity range, capped at 1, constant for barrierless Al + O2; AlO + M dissociation has P0 = 1 and
  E0 = 498 kJ/mol as an upper bound **(assumed)**. **(verified)** `test_line_of_centers_reproduces_thermal_rate` (within
  5 % of the analytic rate for R1, G42, R3 and G33; measured 0.997-0.999), `test_calibrated_probabilities` (R1 within
  0.03 of 0.88, AlO + O2 -> OAlO + O within 0.01 of 0.082, D1 exactly 1; the measured values are 0.882 and 0.0818, and
  OAlO + O 0.531 and AlO + O -> Al + O2 0.00996 were checked by hand once, not pinned by a test).
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
  bound is therefore 0.05 at 1800 K plus 1e-3 at 300 K; the approved spec was not edited. **(verified)** the bounds:
  `test_packet_conservation_and_hot_products` asserts the 1800 K share is below 0.05 and the 300 K share below 1e-3.
  The values themselves (0, 1.3e-4, 1.2e-3 and 0.0166) were measured once, in the Task 6 diagnostic, and are not
  pinned by any test.
- **`hot_start = "droplet_temperature"` (Asha's original picture).** The same Stage A with the droplet speed set to zero
  (the reactions stay on; atoms leave with thermal speed only, so the threshold reactions are expected to be nearly
  closed and r0 to be a few mean free paths, as the spec says; Stage A with a droplet at the ambient temperature and
  zero speed gave 1.12 mean free paths in `test_cold_start_has_no_threshold_reactions`, but r0 of this setting itself,
  at the droplet temperature, was not measured). An earlier comment in `settings.py` that called this "no hot chemistry" was wrong (ruling,
  corrected). **(verified)** that the setting runs through `run.run` with a nonempty hot start and a zero share of
  AlO collisions above the bond: `test_sensitivity_runs_and_tornado_plot`.

### 9.3 Gas chemistry

- **Network.** The neutral reactions R1-R17 of Plane et al. (2021) plus additions from Gobrecht et al. (2022) Table B.1
  (reactions 26, 33, 36, 38 and 42) and the AlO bond dissociation as a hot-stage-only threshold row. All data are in
  `data/wake_column/reactions.csv` and `species.csv` with a source and a validity range per row. **(verified)**
  `test_plane_table_values_reproduced` for four rows only (R1, R11, R6 and R17 at 200 and 300 K; the other rows of
  Plane's table and of Gobrecht's are transcribed and not checked by a test),
  `test_species_table_has_every_network_species` (every species of the network is in the table).
- **Third bodies (D12).** [M] is the total air number density with Plane's N2-measured coefficients for all air; O2 is
  assumed as efficient a third body as N2 **(assumed)**. **(verified)** `test_AlO_lifetime_against_R6`: the AlO
  share in the box at the first output time past 1/(k6 n_total [O2]) equals exp(-t/lifetime) to 1e-3 relative at 70, 80 and 90 km in the night WACCM
  columns (the test prints the lifetimes: 0.093, 2.27 and 45.9 s, against the spec's round "about 0.08, 2, 40 s");
  `test_Al_lifetime_against_R1` likewise for Al against R1 with the ambient O2.
- **Duplicate removed (D13).** The spec's "AlO + O2 -> AlO2 + O and its reverse": the reverse, OAlO + O -> AlO + O2, is
  Plane's R3, so only R3 is kept. **(verified)** `test_reverse_of_G42_is_not_duplicated`,
  `test_hot_only_dissociation_not_in_box_network`.
- **Validity ranges (D14).** Where the source gives none: Plane's (T/300)^(1/6) entries 100-1000 K, R1 23-295 K (Le
  Picard 1997, Gomez Martin 2017), R6 and R7 100-600 K, Gobrecht/Mangan additions 100-2500 K **(assumed;
  `validity_source = assumed` in the table)**. A rate outside its range is clipped and warned about once per reaction per
  run. **(verified)** at the unit level `test_clipping_warns_once_per_reaction` (Review Focus 5: 1000 calls of one clipped
  rate leave one message in the collector; the box's collector is `Box.warnings`, whose contents are reported in
  `summary.json` as `rate_warnings`). The Gobrecht reaction 32 (Al + O2) has
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
  bins), the plan's fallback. Measured once in the Task 8 build (the test runs only the default q): q = 2 gives a
  self-preserving geometric standard deviation of 1.506 by numerical diffusion, q = 2^(1/2) gives 1.4755 and q = 2^(1/3)
  gives 1.4684, against the literature's about 1.46 and the test's band of 1.45 +- 0.03 (extrapolating to q = 1 gives
  about 1.454). **(verified)** `test_coagulation_reaches_self_preserving_width`
  (margin to the 1.48 bound is only 0.0045), `test_split_conserves_number_and_atoms`,
  `test_above_the_last_pivot_the_entity_goes_to_the_last_bin_with_number_scaled`, `test_place_conserves_number_atoms_inside_the_range`,
  `test_bin_count_and_last_pivot_at_ratio_two` (31 pivots, last at exactly 2^31 atoms at q = 2, 487 nm as measured),
  `test_first_bin_is_a_half_nanometre_alumina_dimer`. Every deposit goes through `Bins.place` (ruling), so no atoms are
  lost in the largest bin. Species with two or more aluminium atoms are first-bin clusters **(assumed)**.
- **Atoms delivered per molecule (D9).** (Al, O, H): Al (1,0,0), AlO (1,1,0), OAlO (1,2,0), OAlO2 (1,3,0), AlOH
  (1,1,1), AlOH2 (1,2,2), AlCO3 (1,1,0) (CO2 assumed to leave on sticking). The O released by AlO + AlO -> Al2O + O is
  not tracked, and particle oxidation does not consume gas-phase O, O2 or O3 (the aluminium is a trace) **(all assumed)**.
  **(verified)** `test_atoms_delivered_per_molecule_follow_decision_9` (the species table's `deliver_O` and `deliver_H`
  for all seven aluminium species, as listed); `test_seeded_bins_carry_one_oxygen_per_aluminium` checks only that
  seeded bins start with one O per Al and no H. The box's use of the table (a grown or formed particle receives these
  atoms) is covered by the tally and conservation tests of 9.4, not by a separate test of the O and H totals.
- **Size and density (D10).** Volume per Al atom = (1 - f) 2.8218e-29 m^3 + f 5.3524e-29 m^3 with f = min(H/Al / 3, 1)
  (bulk alumina 3.0 g/cm^3, aluminium hydroxide 2.42 g/cm^3); the tracked O content is a diagnostic that changes no size
  or mass. **(verified)** `test_hydroxide_is_bulkier_per_aluminium` (the ratio of the two volumes per Al atom, 5.3524/2.8218; the
  absolute volumes and the two densities are pinned only through that ratio and the 0.4759 nm dimer diameter of
  `test_first_bin_is_a_half_nanometre_alumina_dimer`), `test_diameter_matches_the_sphere_of_the_composition_volume`.
  Bulk densities for nanometre clusters, compact spheres (no aggregates) and no dehydration of hydroxide **(assumed)**.
  The H/Al reporting threshold of 0.5 for "alumina-like" is a reporting choice, not physics.
- **Formation (D7).** One labelled pathway for each of the 28 pairs of aluminium-bearing species (none lumped), at
  max(k_DD, k_HS) x spin x efficiency x `formation_scale`: k_DD = pi Gamma(2/3) v (2 C6 / k_B T)^(1/3) with the
  Keesom-averaged C6, k_HS = pi sigma^2 v, spin = g_P / (g_A g_B). **(verified)** `test_dipole_capture_AlO_AlO_at_298K`
  (1.516e-9 cm^3/s with 4.36 D, within 1e-3; 1.558e-9 with the table's 4.45 D was measured once), `test_spin_factors`, `test_every_pair_is_its_own_pathway`,
  `test_capture_rate_is_symmetric_and_never_below_hard_sphere`, `test_every_pathway_and_the_smoke_route_are_reported`.
  The 1.5e-9 cm^3/s against Gomez Martin et al. (2017)'s 1.1e-9 is about 40 % higher **(assumed)**; every efficiency
  defaults to 1 and is a sensitivity parameter, so the formation rates are upper bounds with unmeasured efficiency.
- **Growth.** Any aluminium-bearing molecule sticks to a particle of any bin at the free-molecule collision rate
  including its own size times `sticking[species]` x `sticking_scale` (default 1) **(assumed efficiency)**. Smoke uptake
  uses the free-molecule kernel with the true reduced mass of the molecule and one smoke particle,
  mu = m_mol m_smoke / (m_mol + m_smoke), the smoke particle of radius r being a sphere of mass (4/3) pi r^3 rho_smoke;
  rho_smoke is `Settings.smoke_density_kgm3`, default 2000 kg/m^3, the density Megner et al. (2006) state for their
  smoke (Sect. 4.5, p. 4421: "the original density of 2 g/cm3"; they vary it between 0.2 and 5 g/cm^3). The earlier
  infinite-mass limit (mu = m_mol) was wrong in sign as well as size: the assumed 0.2 nm smoke particle weighs about
  40 amu at that density, about as much as AlOH (44 amu), so a molecule meets it faster than it would meet an
  immovable wall, by a factor sqrt(1 + m_mol / m_smoke) (1.45 for AlOH and the 0.2 nm bin). Summed over the bins the
  new AlOH uptake rate is 1.30 times the infinite-mass one (the same factor, to three figures, at 60, 76 and 90 km), so
  the infinite-mass limit underestimated the uptake by about 23 % **(measured once)**. The density enters the run name (`smkrho`) when changed, and the sensitivity set does not
  vary it. **(verified)** `test_smoke_of_finite_mass_is_taken_up_faster_than_infinitely_heavy_smoke` (the uptake of a
  single 0.2 nm bin equals a hand computation with the reduced mass to 1e-12, exceeds the infinite-mass value, and
  approaches it for a very heavy smoke), `test_smoke_uptake_equals_the_direct_sum_at_a_table_altitude`.
- **Oxygen uptake (D11, ruling).** Rate per bin = sum over O, O2 and O3 of collision rate x `oxidation_eff` x g(O/Al),
  one O atom per event; g is 1 below O/Al = 1.485, 0 at 1.5 and above and a cubic smoothstep (3s^2 - 2s^3) between, not
  the plan's linear ramp: the linear ramp's slope jump (-67 to 0) stalled the BDF Newton iteration and let O/Al drift
  to 2-3. The ramp is a numerical cap, not data; O has its own absolute tolerance (1e-18 scaled). Hydrogen is gained only
  through condensing hydroxide molecules (no gas-surface H or H2O reaction) **(assumed)**. **(verified)**
  `test_oxygen_uptake_factor_runs_from_one_to_zero_with_a_continuous_slope`,
  `test_seeded_particles_take_up_oxygen_to_one_and_a_half_per_aluminium` (within 1e-3 in 1 s). Measured once in the
  Task 8 build, not by a test: with the smoothstep the O/Al of every bin holding more than 1e-9 of the inventory agrees
  with a tight (relative tolerance 1e-10) solution to 1e-2 in the 76 km day and night, 60 km night and 90 km day cases.
- **Background meteoric smoke.** Number density and size distribution from Megner et al. (2006) Fig. 2, digitised to
  `data/wake_column/megner2006_smoke.csv`. Fig. 2b gives the size distribution at 90 km only and shows only about
  16 % of Fig. 2a's number (the smallest smoke is below the plotted axis), so (ruling) the deficit is put in one bin at
  Megner's smallest radius, 0.2 nm, and other altitudes use the 90 km shape scaled by Fig. 2a **(assumed)**. That bin
  holds 4.0e4 of 4.8e4 per cubic centimetre at 90 km and about 47 % of the smoke surface area. With the
  reduced-mass kernel that bin supplies 71 % of the AlOH uptake (70.7 % at 90 km and 241 K, 70.8 % at 76 km, 70.7 % at
  60 km; 66 % for Al, 72 % for Al(OH)2), so removing it would cut the uptake by a factor of 3.4 **(measured once, not
  by a test)**; the smoke uptake, and with it the 12-hour condensed fraction, therefore rests mainly on this assumed bin. The smoke is a fixed background, not depleted, not coagulated with plume particles,
  and sticks at efficiency `smoke_sticking` (default 1) **(assumed)**. **(verified)** `test_smoke_table_integrates_to_fig2a_number_density`
  (the table's totals match the Fig. 2a values recorded in its CSV header to 1 % at 60, 75 and 90 km: a check of
  the table against the read-off, not of the read-off itself), `test_smoke_uptake_equals_the_direct_sum_at_a_table_altitude`, `test_smoke_uptake_outside_the_table_raises_naming_altitude_and_file`.
  The uptake rate of AlO at 200 K, measured once after the reduced-mass fix (smoke density 2000 kg/m^3), is
  4.1e-7, 4.8e-6 and 1.1e-5 per second per molecule at 60, 75 and 90 km (a time scale of about 24 hours at 90 km;
  with the earlier infinite-mass kernel they were 3.2e-7, 3.7e-6 and 8.8e-6).
- **Coagulation.** Free-molecule Brownian kernel between all bins (half for a self-pair), particles compact spheres
  **(assumed)**. **(verified)** `test_coagulation_reaches_self_preserving_width`.
- **Settling.** Bin k leaves at 2 v_s / (pi r), v_s the Epstein (free-molecule) speed, through the lower half of a
  well-mixed cylinder **(assumed)**; settled atoms keep their route and drop their O and H (Al-only settled tally).
  **(verified)** `test_settling_only_decays_at_two_vs_over_pi_r` (the seeded bin decays as exp(-rate t) within 1e-4
  relative at every output; the worst error measured in the Task 8 build was 1.6e-5).
- **Irreversible condensation** (Al-O clusters are bound by 400-500 kJ/mol per atom, and the plume is at ambient
  temperature within a fraction of a millisecond) and **particles diluted at the gas's rate** are **(assumed)**.
- **Provenance closure.** The tallies by source and route sum to the condensed aluminium. **(verified)**
  `test_provenance_closes_in_every_bin` (the tallies equal N times the bin's atoms within 1e-6 relative or 1e-6 of
  the inventory in atoms, and the 21 provenance shares sum to particles + smoke + settled within 1e-9, at every output;
  the Task 8 build measured 3e-12 atoms and 9e-19).

### 9.5 Numerics, outputs and definitions

- **Solver (D19).** `solve_ivp` BDF, rtol 1e-6 and atol 1e-12 on scaled variables: aluminium counts divided by the
  inventory, particle number of bin k by inventory / n_k (so the tolerance bounds atoms), each ambient count by
  c_ambient V0, r^2 by r0^2. The Jacobian is analytic for every term except the r^2 column (a one-sided difference).
  **(verified)** `test_analytic_jacobian_matches_finite_differences` (every row outside the ambient gas agrees to 1e-4 of
  the row maximum; the kernels' weak dependence on H/Al is left out, which affects only Newton convergence). Measured
  once in the Task 8 build, not by a test: condensed fractions at 76 km agree to 4-5 significant figures with a solve at
  relative tolerance 1e-10. A 12-hour run takes about 6-25 s depending on the case (measured with seven runs in parallel on an eight-core machine).
- **"Condensed" (D17)** is everything that left the gas phase: plume particles + background-smoke uptake + settled.
  t50 and t90 are the first output times at which it reaches 0.5 and 0.9, or null. **(verified)**
  `test_t50_and_t90_are_the_first_output_times_reaching_the_level` (the helper `run._first_time_reaching` on a synthetic
  series: exact crossing counts, 0.89 does not reach 0.9, a level never reached gives null) and
  `test_condensed_is_particles_plus_smoke_plus_settled` (the definition); no test runs a real case through a 0.5 or 0.9
  crossing (the baseline cases never reach them), and `test_short_run_file_contents` checks that the keys exist and that
  the final value matches the history file.
- **Steady state (D18).** The steady time is the earliest output time t_ss such that, for every t at or after it with
  10 t at most t_end, every aluminium reservoir changes by less than 0.01 of the inventory between t and 10 t; at
  least one decade of evidence is required (null otherwise). Reported, not used to stop. **(verified)**
  `test_steady_state_time`, `test_steady_state_needs_a_decade_of_evidence`.
- **`d_volume_mean_nm` (ruling)** is the volume-weighted mean diameter, the sum of N d^4 over the sum of N d^3 over the
  plume particles at the end of the run (where the condensed volume sits); the alternative, the diameter of the
  mean-volume particle (sum N d^3 / sum N)^(1/3), is not reported. `d_number_mean_nm` is the plain number mean.
  Smoke-taken aluminium is not in either. **(assumed)** definition, implemented in `run.summarise` and not pinned by
  a test.
- **Output masking (ruling).** In `size_distribution.csv` a bin holding under 1e-9 of the inventory in aluminium has
  empty composition columns (its O/Al and H/Al are set by a handful of atoms) and its alumina diameter; N and the tallies
  are kept. **(verified)** `test_short_run_file_contents`.
- **Skipped runs.** A droplet that vaporises no aluminium gives a summary with `skipped: "no aluminium vaporised"`
  and no outcome keys, and no CSVs; the plot scripts leave such runs out. **(verified)** `test_zero_inventory_is_skipped_cleanly`,
  `test_plot_scripts_ignore_a_skipped_run`. A truncated or unreadable `summary.json` is rerun, not fatal
  (`test_unreadable_summary_is_rerun`).
- **Run names** encode the case, the ambient and every non-default setting, so a sweep skips a finished run
  (`test_every_setting_changes_the_run_name`, `test_existing_summary_is_reused`). `Settings` rejects an unknown
  `hot_start` or `ambient` value (ruling: a typo would silently change the run name; `test_bad_hot_start_or_ambient_raises`). The test-only `trail=` argument of
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
  droplet-temperature start, eddy x0.1, x100 and zero (molecular only), dissipation x0.1 and x10, smoke sticking x0.1
  and 0. 15 x 6 = 90 runs; `sensitivity.py` also queues the six baselines itself and skips them when finished runs
  with the same name exist (the default sweep's), so the tornado plot works on a folder holding only this set. **(verified)**
  `test_sensitivity_set_is_15_settings_over_6_cases` (each setting changes exactly one field). The plot scripts pick
  baseline runs by their physics settings, so shortened test runs and cross-check runs sharing a folder are not mixed
  (`test_baseline_signature_ignores_numerics_but_not_physics`).

### 9.7 Findings from the build, to carry into the research

- **At 76 km almost nothing nucleates in the plume.** A 14.1 micrometre droplet's 2.9e13 aluminium atoms spread over a
  366 m trail with r0 = 21.6 mm; in the Task 8 diagnostic (3e13 atoms, L = 500 m, r0 = 2 cm) molecular diffusion
  (0.73 m^2/s) took the plume to a radius of 1.7 m within a second, a 7000-fold dilution in volume, and formation, being
  second order, stops at once. Plume particles hold about 3.1e-5 of the
  inventory by day and never grow past the first bin (0.485 nm by day, 0.482 nm by night); the condensed fraction at
  12 hours is 0.217 by day and 0.219 by night, nearly all of it first-order uptake onto background smoke (by day 0.2174
  condensed = 3.1e-5 particles + 0.2174 smoke + 8.9e-7 settled; AlOH 0.144 and Al(OH)2 0.071), about 78 % is still gas
  (AlOH 0.50 and Al(OH)2 0.27), and t50, t90 and the steady time are null. Across altitude (14.1 micrometres, by day,
  same settings) the condensed fraction at 12 hours is 0.0174 at 60 km, 0.217 at 76 km and 0.417 at 90 km, of which the
  plume particles hold 2.0e-4, 3.1e-5 and 4.2e-7 and smoke uptake nearly all the rest (the earlier infinite-mass smoke
  gave 0.0126, 0.168 and 0.340). Re-run 2026-11-05 with the scale-dependent eddy diffusion in place of the 30 m switch
  (all 126 default cases, measured once): the condensed fraction at 12 hours changes by at most 6 % (largest at 200
  micrometres near 62 km; 0.0174, 0.217 and 0.417 above are unchanged to three figures) and plume particles fall by 3 to
  18 % (14.1 micrometres by day: 1.77e-4, 2.78e-5 and 3.73e-7 at 60, 76 and 90 km). Particles from AlO in the first
  milliseconds are unchanged below about 76 km; aluminium reaching particles from Al(OH)2 falls by more than half.
  This rests on the assumed 0.2 nm smoke bin, a smoke sticking of 1, the diffusion coefficient of D4, and one altitude;
  lower altitudes, where r0 is smaller (1.7 mm at 60 km), may nucleate more. Measured once
  (`run.run` on 14.1 micrometres, 76 km, day, WACCM, 5.4 s, from the commit that fixed the smoke mass; no test pins these values); the physical conclusion is **(assumed)** until the
  sensitivity set has been run.

### 9.8 The segmented trail (`segments.py`)

Built 2026-11-05. The default run (`Settings.n_segments = 24`; 0 restores the single column) cuts the path to 99.9 % of
the aluminium loss into equal lengths, each holding the aluminium really released there, with its own hot start,
altitude, ambient air and clock. The last segment also holds the remaining 0.1 %. Each is its own `Box`, and they are
added at the same time since release (`TrailSolution`). Spreading along the trail is applied afterwards (`along_trail`).

- **Verified by tests** (`tests/test_wake_segments.py`): cutting conserves the aluminium, lengths are equal, mass follows
  the release; a quarter of the atoms in a quarter of the volume evolves identically; aluminium is conserved at every
  time since release including the unreleased share; each segment is read at its own age; `along_trail` keeps every
  atom, equals the closed form for a top hat, and matches the box's axis concentration for one segment.
- **Assumed: a segment appears whole at its aluminium-weighted release time.** Laying really takes about 0.03 s near
  the peak and about 0.12 s in the tail, so time series before the droplet finishes (about 1.5 s for 200 um at 70 km)
  are not to be read in detail. In a test without chemistry the particle formation rate at 0.1 s was 19 % of a
  continuous trail's, while the total integrated over time differs only by the variation of release inside a segment.
- **Assumed: spreading along the trail can be separated from the chemistry.** Tested once each at 70 and 60 km
  (200 um, by day, single column, 12 hours) with four times less and four times more aluminium: the ratio of Al(OH)2 to
  AlOH is the same to all printed digits (2.89 at 70 km, 27.62 at 60 km) and smoke uptake per gas atom is unchanged
  (0.01746 at 60 km). Only plume particles respond to concentration. Spreading is assumed equal along and across the
  trail, and neighbouring segments are assumed to hold nearly equal shares by the time they overlap. No coupled run
  has been made; that option is kept in `docs/milestones/2026-11-05-coupled-trail-segments.md`.
- **Measured, 24 against 48 segments** (two cases, by day, 12 hours). Plume particles rise by 0.9 % at 200 um, 70 km and
  by 2.6 % at 14.1 um, 76 km with 48 segments; gas species and smoke uptake change by at most 0.03 %.
- **Measured, 24 segments against the single column** (the same two cases). Plume particles rise by 12.9 % and 9.9 %;
  the condensed share rises by 3 % at 200 um, 70 km and is unchanged at 14.1 um, 76 km; gas species fall by 0.3 % at
  200 um, 70 km. At 200 um, 70 km, the peak atoms per metre of trail at 12 hours is 2.85e13, about half of the 5.6e13
  at 1 second, and its position moves from 389 m to 986 m along the path.
- **Measured, the default sweep with 24 segments against the single column** (all 126 WACCM runs, 12 hours, commit
  eba9a05; aluminium budgets close to 2.4e-13 of the inventory). Plume particles change by the factor
  (24 x the sum of the squared segment shares) x (L at 90 % of the loss / trail length at 99.9 %), which predicts the
  measured ratio to within 4 % in every run: the first factor is 1.39 to 2.22, the second 1/1.34 to 1/1.95. The ratio is
  1.09 to 1.16 from 64 to 90 km and rises to 1.54 at 56 km, because the trail is then only 1.36 times L. Below 56 km the
  particle shares are under 1e-8 of the inventory. The condensed share changes by more than 1 % relative only for 200 um
  at 60 to 66 km (5.64 % to 6.66 % at 62 km, the largest change); gas species and smoke uptake change by at most 1.1 %.
  The undercount from stepping the release over 24 segments, without chemistry, is at most 3.2 % on all 63 trails.
- **Not done:** the sensitivity set and the NRLMSIS + AFGL cross-check have not been re-run with segments.
  `plume_video.py` still shows the single column.
- **Open (found 2026-11-05):** the box applies concentration on the axis to the whole volume, so aluminium-with-
  aluminium rates count c^2 V, while a Gaussian plume has half of that. From Stage A's stopping points at 70 km, with
  molecular diffusion and no chemistry, the time-integrated overlap to 10 ms is 0.52 of the box's. Particle formation
  is therefore over-counted by about twice in both the single column and the segments. Asha chose on 2026-11-05 to leave
  this uncorrected for now, so quoted plume particle shares carry this bias.

## 10. Level-set architecture: geometry (stage 1)

This section describes the geometry layer in `reentry_model/levelset/` (design file, section 20). It is built and
tested but not yet coupled to any physics: every test and measurement uses prescribed fields, with no thermal solve,
skins or flow.

**Mesh.** The box is a grid of cubes of edge h0, each split into six Kuhn tetrahedra along its main diagonal and refined
by Maubach's newest-vertex bisection with tags. The refined mesh stays conforming (no hanging nodes; boundary area equal
to the box's to about 1e-12) and has at most three tetrahedron shapes. Three bisections halve the cube edge, so the
achieved size is the largest h0/2^k that does not exceed the target: with h0 = 2.8 mm the fine size 0.7 mm is met
exactly, and the 0.25 mm cell floor becomes 0.175 mm. Coarsening rebuilds the mesh from the base grid, so a mesh is a
function of the base and the targets only.

**Level sets.** Fields are piecewise linear (P1) on the tetrahedra and negative inside the metal. The zero contour in
each cell is a flat patch; a node exactly on the contour counts as outside, and zero-area patches are dropped. The
distance reset recomputes each node's geometric distance to the patch triangles (not an eikonal solve), searching the
eight nearest triangle centroids; with a band, only nodes within the band are computed and the rest take plus or minus
the band, never more than their true distance. (The search radius is widened by the largest distance from a triangle's centroid to its vertices, so every node within the band of a triangle is computed exactly.)

**Volume correction.** After motion, each connected piece's level set is shifted so that the piece's volume is restored
exactly (relative error about 1e-13 in the tests, including two pieces three cells apart). Each shift excludes the other
pieces' nodes and cells. Pieces closer than about two cells share faces of metal-holding cells and count as one piece,
so pieces join and separate at the scale of the cell, not of the surface.

**Refinement reach.** A cell's test region is grown by 0.94 of its size, the distance from a cell's centroid to its
farthest vertex (0.935 h for a Kuhn tetrahedron). The larger square root of three times h would have widened the fine band
by about 1.2 mm at 0.7 mm cells.

**Shadow test.** A ray starts a twentieth of a cell outside each patch and is marched upstream with steps equal to the
level set but at least a quarter of the smallest cell, for at most as many steps as that minimum step needs to cross the box diagonal (plus ten); the test raises an error only if a ray is still marching after the last step. A whole-cell start offset carried leeward rays near a sphere's
outline outside the body, so they missed it.

**Thin test.** A patch is thin when the molten depth is at most c_kappa times the smallest principal radius of curvature
and the depth's slope along the surface is at most c_slope. The thresholds are starting values only, to be set by
verification: c_kappa 0.1, c_slope 0.2, hysteresis hold of 3 consecutive opposite answers, and a connected group of at least 7 patches
before a mode changes. The molten depth is the nearest distance from the patch to the zero contour of the molten-rigid surface, which is never more than the depth along the normal that the design names; the two agree where the surfaces are nearly parallel, and the approximation is revisited in stage 3. Where the metal is fully molten (no rigid metal, so the molten-rigid surface is plus infinity), the depth is infinite, the patch is thick, and the curvature comes from the outer surface alone.

**Accuracy and cost.** Cut-cell volume and area of a sphere converge at second order: errors of 2.6 per cent (volume)
and 1.3 per cent (area) when the radius spans 4.4 cells, and 0.65 and 0.33 per cent at 8.8 cells; the P1 volume is
always low. Costs are in the README; they were measured once, serially, on the laptop on 2026-11-11: about 63 s of geometry per macro step (every timed operation except the one-off first build, now including the outer-surface rebuild, the molten-rigid-surface rebuild and one motion of the skin base) for the 50 mm sphere and about 330 s, or 5.5 minutes, for the 100 mm sphere, at peaks of 1.56 and 3.30 gigabytes. The cost of cells at the 0.175 mm floor is not yet measured, because the sphere has none.

**Not yet verified.** The thin-test thresholds; behaviour on moving, merging or splitting bodies (the tests move only
prescribed spheres and boxes); correctness of the derived surfaces against real temperature fields; the shadow test on
moving or melting shapes (the overhang test covers a stepped stem under a cap); and whether every operation must run every macro step, which is for stage 3 to measure.

## 11. Level-set architecture: the thermal solve (stage 2)

This section describes the heat solve on cut cells, `CutThermal` in `reentry_model/levelset/cutthermal.py` (design
file, §7 and §20.11 item 2), added on 2026-11-11. It solves heat conduction in the metal of a level-set mesh, with the
boundary conditions applied on the surface pieces (patches) of the cells the surface crosses. It is built and tested but
not coupled to the flight loop, and it has no skins (stage 3), no flow and no heat carried by a flow (stage 4). The
package's tests use test-only materials (constant properties, and a linear melting range); the plate comparisons that
use the real Scheil alloy (AA7075, whose melting follows the Scheil solidification model, with a step of latent heat
at the eutectic) were run outside the package, as described under "Plates" below.

**Metal-only integration.** The temperature is piecewise linear (P1) at the nodes. Each cell holding metal (an active
cell) has its metal part split into sub-tetrahedra, and its heat capacity (the consistent mass matrix) and conduction
are integrated over that part only, so no heat is stored in or conducted through empty space. The unknowns are the
nodes of active cells; other nodes keep their values and leave the system.

**The cell temperature for coefficients.** Each cell evaluates its conductivity and heat capacity at one temperature,
T*, the average of its nodal temperatures weighted by the mean of each shape function over the metal part. With this
choice the heat the cell's mass matrix stores equals density times specific heat times metal volume times the change
of T*, exactly, so the stored enthalpy is the sum over cells of density times metal volume times h(T*), where h is the
specific enthalpy in joules per kilogram.

**The secant capacity.** The heat capacity in a step is the secant [h(T*) − h(T*_old)] / (T* − T*_old) rather than the
tangent c_p(T), so that the change of stored enthalpy equals the heat that came in to the Newton tolerance whatever the
shape of h(T), latent heat included. Measured on the 2026-11-11 prototype (laptop on mains power, shared with running
flights, each case measured once): the books close to 3 × 10⁻¹³ relative with constant properties, and to about
4 × 10⁻⁸ relative on the Scheil alloy at a Newton tolerance of 10⁻⁸ (6 to 8 Newton iterations a step), and 8 × 10⁻⁸ on
a linear melting range (up to 67 iterations on the first step, about 11 after). Closure at the tolerance is the point:
a books error much larger than the tolerance would mean heat created or lost by the scheme.

**Boundary heat.** On each patch the heat into the metal is q0 + q1 T_p − ε σ (T_p⁴ − T_amb⁴), with q0 the prescribed
flux part of the law (watts per square metre into the metal), T_p the temperature at the patch centroid, q1 at most
zero (the skins' base law, implicit), ε the emissivity and σ the
Stefan–Boltzmann constant; the radiation is linearised in the Newton iteration. A flow's heat flux is supplied by the
caller per patch, by projected area: q times max(n·d, 0), with n the patch normal and d the direction the flow comes
from. Applying the flux only to patches tilted less than about 25 degrees from the flow (tried first) lost 41 per cent
of a plate column's heat once its receding top crossed cells cut by its walls; the projected-area rule delivered 98 to
99 per cent (prototype, 2026-11-11, measured once).

**Accuracy.** A slab under a constant 2 megawatts per square metre for 0.5 s, its sides cut off the grid: the rise of
its surface temperature differs from the semi-infinite solution by 1.38 per cent at 1 mm cells, 0.37 per cent at
0.5 mm and 0.04 per cent at 0.25 mm, which is second-order convergence (prototype, 2026-11-11, measured once). The time
step (0.01 against 0.0025 s) and lumped against consistent mass change it by under 0.2 per cent. Where one cell holds
both the top and a side wall of a body (a sharp corner) its patch tilts, and on 1 mm cells the projected heat on a 4 mm
wide slab was 3.3 per cent short; rounded bodies have no such corners.

**The molten-rigid surface tracks a Stefan front.** The design (§20.11 item 2) asks that the molten-rigid surface ψ (the
895.1 K isotherm, built by stage 1's `molten_rigid_surface` from the computed temperature, negative in rigid metal)
follow a melting front to within one fine cell, 0.7 mm. The test
`test_psi_tracks_a_one_phase_stefan_front_within_one_fine_cell`, added on 2026-11-11, heats a 4.2 mm square column,
22 mm deep, on uniform 0.7 mm cells, its top 0.3 mm above a layer of nodes and its sides and bottom insulated. The metal
has constant properties (density 2800 kg/m³, conductivity 150 W/m/K, specific heat 900 J/kg/K) and a latent heat of
390 kJ/kg released linearly between 890 and 900 K; it starts solid at the bottom of the melting range (890 K), so only
the liquid conducts (the one-phase problem). The top is held near 995 K, 100 K above the middle of the melting range, by
the implicit base law with a conductance of 10⁷ W/m²/K, for 2 s in 0.01 s steps, and no metal is removed. The depth of
ψ = 0 on the column's centre line is compared with the exact (Neumann) front depth 2λ√(αt), where α is the thermal
diffusivity, about 6.0 × 10⁻⁵ m²/s, and λ = 0.328 is the root of λ exp(λ²) erf(λ) = St/√π with the Stefan number
St = c_p (T_hot − T_m)/L = 0.231 (the ratio of sensible heat in the liquid to latent heat). The computed front was
0.112 mm deeper than the exact 2.528 mm at 0.25 s, 0.049 mm deeper than 3.575 mm at 0.5 s, 0.006 mm deeper than 5.056 mm
at 1.0 s and 0.010 mm shallower than 7.150 mm at 2.0 s, that is at most 16 per cent of a cell, against the bound of one
cell. The error is largest at the earliest time and within 0.01 mm from 1 s on, so it most likely comes from the first
steps, when the front is within a cell or two of the surface; that cause was not isolated. The test takes 25 to 35 s on
the shared laptop and its errors were recorded from one run. It checks a planar front with a narrow (10 K) melting range
only; the Scheil alloy's wide range is covered through the plates' non-rigid depth, below.

**The ghost penalty and its weight.** A cell holding a sliver of metal would make the system ill-conditioned. A ghost
penalty (Burman 2010) on every interior face between two active cells of which at least one is cut penalises the jump
of the normal temperature gradient, with the ghost weight γ times the larger conductivity, the larger cell size and
the face area; it leaves functions linear across the face, and so constants, untouched. Without it (prototype, 2026-11-11,
measured once) a cell holding 10⁻¹⁰ to 10⁻⁴⁵ of its volume gives a diagonal entry 10⁻¹⁰ to 10⁻⁴⁵ times the largest,
the sliver node's temperature is undetermined, and plain conjugate gradients fail in 5 000 iterations. With γ = 0.01
the ratio is about 6 × 10⁻³ and plain conjugate gradients converge in about 330 iterations. A larger weight conditions
the system better (γ = 0.1 gives 5 × 10⁻² and 140 iterations) but disturbs the solution: with a surface receding
0.05 mm a step through 1 mm cells under 1 megawatt per square metre, the temperature 4 mm below changes smoothly at
γ = 0 or 0.01, while γ = 0.1 gives one jump of 0.07 K and γ = 1 four jumps of up to 0.24 K, each when the surface
crossed a layer of nodes and the set of penalised faces changed at once. The default is therefore γ = 0.01.

**Newton iteration and damping.** Each backward-Euler step iterates to a relative change of the temperature of 10⁻⁸,
at most 100 iterations. The boundary flux is linearised by Newton; the conductivity and the secant capacity are lagged
at the current iterate (a Picard iteration). Because the effective heat capacity has kinks (the solidus, the liquidus),
the update is halved, down to a floor of one sixteenth, whenever an iteration changes the temperature by more than
0.9 of the previous change. Linear systems are solved by conjugate gradients with a smoothed-aggregation algebraic
multigrid (AMG) preconditioner to 10⁻¹², or by a direct solver. A step that does not converge raises an error.

**A moving surface on the same mesh.** `move` changes the level set without changing the mesh. The nodal temperatures
are kept; nodes that become active are filled ring by ring with the mean of their already-valued neighbours, and if any
such node is still without a value after the rings the method raises a RuntimeError rather than keep a stale value.
It returns the enthalpy the move removed with the metal.

**The uniform offset on rebuilds.** `remesh` carries the temperature to a rebuilt mesh: it extends the field three rings
outside the metal, projects it conservatively onto the new mesh, and then adds one uniform temperature offset to every
node, found by bisection to 10⁻¹² K, so that the stored enthalpy equals the old one. The search for the offset starts at
plus or minus 50 K and doubles until it brackets the old enthalpy; if, after that, the enthalpy does not close to 10⁻¹⁰
of the scale below, the method raises a RuntimeError naming the residual. The residual before the offset is returned,
divided by a scale that is the larger of the old stored enthalpy and the energy of a 1 K change of the old body (its
heat capacity, the sum of density times metal volume times effective specific heat over the cells, times 1 K). The
stored enthalpy is measured above the material's reference temperature (293 K for the test materials), so for a body at
or near that temperature it is close to zero and could not serve as the scale alone; for a body hundreds of kelvin above
it, as in every plate, the scale is the stored enthalpy and the residual is relative to it, as before. A move or a
rebuild that raises leaves the solver exactly as it was before the call. The residual is the measure of the transfer's
quality: the offset spreads it over the whole body, so it must be small. On the prototype plates (below) it was within
3 × 10⁻⁶ (cold plate) and 3 × 10⁻⁹ (hot plate) with the fine margin into rigid metal, and up to 4.6 × 10⁻⁴ and
2.4 × 10⁻³ without it, because without the margin coarse cells rebuilt just under the melting front carry the error.

**Plates.** The comparisons with the skyscrapers' plates were run on 2026-11-11 in the prototype folder
`prototype/levelset/stage2/` (commit 7ac46c5 on branch `levelset-two-interface`, not brought to `main`), against read-only copies of the skyscrapers build's Scheil material
and one-dimensional reference (build commit 3c16f4a), because both exist only in that build until the shared-package
refactor (stage 0 + S) brings them (plan Decision B). Each case was measured once on the laptop on mains power, shared
with running flights. The column is the level-2 plate's centre with insulated walls; liquid is removed after every step
wherever T is at or above T_feed = 910 K (the feed temperature); the reference runs at 2 µm cells and 1 ms steps.

- On uniform 0.7 mm cells with 0.05 s steps, the cold plate (700 K, 6 MW/m²) recedes 0.800, 2.688 and 5.130 mm at 1.0,
  1.5 and 2.0 s against the reference's 0.788, 2.712 and 5.216 mm (1.5 per cent ahead, then 0.9 and 1.6 per cent
  behind); the hot plate (880 K, 2 MW/m²) recedes 0.798, 1.806 and 2.904 mm against 0.788, 1.816 and 2.916 mm
  (1.3 per cent ahead, 0.6 and 0.4 per cent behind), and its non-rigid depth (the depth of metal that is not rigid,
  under the molten-rigid surface ψ built by stage 1's `molten_rigid_surface`) is 1.514, 1.630 and 1.724 mm against
  1.488, 1.570 and 1.676 mm (1.7, 3.8 and 2.9 per cent thicker). The energy books close to 10⁻⁸. At 0.35 mm cells the cold plate's
  recession moves by at most 2.3 per cent, a first look at fine-size convergence that stage 4 judges on a pool. These
  runs took 8 s of wall time for 2 s of plate.
- At 0.5 s steps the same columns fail: the cold plate is 1.8 per cent ahead, then 10.9 and 12.8 per cent behind, and
  the hot plate 24.6, 15.1 and 11.9 per cent behind. The liquid stays on the plate through each long step and absorbs
  heat until the removal at the end of the step, where the reference removes it every millisecond. Stage 2 removes
  liquid only once a step; the skins, which take it continuously, are stage 3, so the judgement at 0.5 s steps moves
  to stage 3 (plan Decision A).
- On the adaptive mesh with the fine margin (2.8 mm base cells, 0.7 mm where ψ > −d and above, with d = 4 mm, rebuilt
  after every step, walls on grid lines), the cold plate recedes 3.6, 4.0 and 4.2 per cent behind the reference at
  1.0, 1.5 and 2.0 s, inside the skyscrapers' 5 per cent acceptance; the hot plate is 9.6 per cent behind at 1.0 s (not
  judged) and 4.3 and 2.6 per cent behind at 1.5 and 2.0 s, and its non-rigid depth is 3.8 to 4.3 per cent thicker,
  inside the 10 per cent acceptance.
- Without the margin (d = 0, 0.7 mm only above ψ) the cold plate is 7.4 per cent behind at 1.0 s, which fails, and the
  hot plate's non-rigid depth is up to 6.6 per cent thick. The margin therefore earns its place.

**Cost.** On the 50 mm sphere (0.7 mm fine cells over a 4.6 mm shell; 1 055 790 cells, 659 028 active, 121 191 active
nodes) the prototype measured on 2026-11-11 a setup of 8.6 s, 6.6 s for one 0.5 s step with radiation (5 Newton
iterations, the AMG hierarchy rebuilt every iteration) and a peak memory of 1.79 gigabytes. The package's script
`analysis/levelset/measure_thermal.py`, run once on 2026-11-11 on the same mesh (laptop on mains power, shared with
running flights), measured a setup of 6.3 s, steps of 4.9, 4.7 and 4.7 s (5 Newton iterations each) and a peak of
1.99 gigabytes; the times are about 27 per cent below the prototype's and the memory about 11 per cent above, within
the scatter of a shared laptop. The thermal step is thus under a tenth (about 7.5 per cent) of the 63 s of geometry
per macro step at 50 mm (stage 1), so it does not set the cost of a flight; the two were measured on different meshes
of the 50 mm sphere (1 055 790 cells here, 1 290 192 for the geometry), so the ratio is approximate. The 100 mm sphere is
*estimated* at about four times as much (about 30 s a step), against about 330 s of geometry a step; it was not measured.

**Not yet verified.** The solve with skins on its patches (stage 3), with a flow and advection (stage 4) and in a
flight; the plates at 0.5 s steps (stage 3); fine-size convergence on a pool at 0.7, 0.5 and 0.35 mm (stage 4); the
cost at 100 mm.

## 12. The skins' spraying solved on their own speed profile (`--dispersion`)

This section describes the profile path of the layered skins' spraying, built on 2026-11-12 from
`superpowers/specs/2026-11-11-kinked-profile-dispersion-design.md` (the design, with Asha's nine decisions) and its plan.
Its physics lives in the shared package (`reentry_model/shared/`: the solver and the floor test in `dispersion`, the fit,
the layer's properties, the release-run rule, the mode inputs and the history summary in `melt_profile`, the mode choice
and the record layout in `spray`); the skins only read their cells into those arrays and take cells off the top.

**The setting and its default.** `run --dispersion table|profile` chooses the theory a spraying skin run
(`--surface-model skin --melt on --skin-flow on`) sprays by. The default is `profile` (Asha, 2026-11-11, decision 9).
On it, each windward column whose top cell flows and moves has its cells' speeds fitted as a continuous line of straight
segments (the speed profile), a wall is placed by the floor test below, Girin's gradient mode is solved on the profile
above that wall, and the faster of that mode and Girin & Kopyt's gas-driven mode (with the film depth taken as the depth
to the wall) releases the column's cells from the top down while they are at least as warm as the mode's threshold
temperature, liquid and slurry together (the heat-limited release, decision 5). A wet windward column whose profile is
not resolved (its top cell does not flow, does not move, or has the rigid front at or above its centre) falls back to
Girin's theory. `--dispersion table` switches back to Girin's theory as it stood before, with his
theoretical layer depth and surface speed. The element model has no resolved profile and is always on `table`. Both
shear modes act only where Girin's continuum closure (or the slip band's balance, with `--slip-band on`) gives a layer
depth, as before (decision 6). `--dispersion` given to any other kind of run is refused with exit code 2, and a state
saved before the option existed resumes on `table`.

**Verified (tests in `tests/test_reentry_model_dispersion.py`).** The broken-line solver, the inviscid model of a profile
made of straight segments with a free surface on top and a wall below, is our own derivation (the side session's, not
published). It is checked against Girin's printed equation and against Ghosh et al. (2017): on Girin's profile (one
linear layer over stagnant liquid) it reproduces the roots of Girin's Eq. 1 cubic to 1e-10 at four pairs of wavenumber and
Weber number, and its fastest wave at a surface Weber number of 10 000 matches the cached table of the `table` path
within the table's grid spacing (`test_girin_profile_gives_the_roots_of_his_cubic`,
`test_girin_profile_fastest_wave_matches_the_cached_table`); it reproduces the side session's closed form for a slurry
over a wall to 1e-10 (`test_slurry_over_a_wall_gives_the_closed_form_roots`); and it reproduces Ghosh's limits, a peak
growth rate of 0.229 at a slope ratio of 1.5 and the stability of a single straight film on a wall and of a more viscous
fluid on top (`test_ghosh_limits`).

**The threshold.** Solved exactly, Girin's Eq. 1 becomes unstable above a surface Weber number of 3.00, not the 3.08 Girin
printed (`test_girin_threshold_is_three`: stable at 3.00, unstable at 3.01). On the profile path no critical Weber number
is applied to the gradient mode solved on a column's profile (branch 5): whether it grows is the solver's own answer, and
it sprays only where its growth rate exceeds the viscous damping of its own wave (twice the layer's Couette viscosity
times the wavenumber squared) and one whole wavelength fits in the column's region (the contiguous windward columns
whose layer above the wall is at least `spray.B_MIN` deep), so the threshold there is the solver's 3.00 on Girin's
profile. `--we-critical` (default 3.08) applies on `table`, and on a `profile` run only to the columns whose profile is
not resolved (the wet windward columns whose top cell does not flow, does not move, or has the rigid front at or above
its centre): they keep Girin's theory, their thick branch gated by it as before. A resolved column never falls back to
Girin's theory. One without a finite layer depth (delta_m, Girin's conjugate layer depth), or whose layer above the wall
is shallower than `spray.B_MIN` (1 µm), takes no shear mode at all, neither the gradient mode, the gas-driven mode nor
any of Girin's branches (decision 6; `spray.py` removes every resolved column from those branches and gives the two
modes only to the resolved columns with a finite delta_m and a wall at least `B_MIN` deep, which
`test_no_delta_m_no_shear_mode_on_resolved_columns` pins). Only the Rayleigh–Taylor mode, which has no Weber-number
gate, can still take such a column.

**Assumed.**

- *The floor test is a scaling argument* (the design's section 3), assumed until stage 6, the viscous (Orr-Sommerfeld)
  check that is to calibrate it. A breakpoint of the fit where the shear rate below is at most `KINK_DROP` = 0.2 times the
  shear rate above (a five-fold drop) is a kink. Kinks are tested from the top down, and the first whose floor Reynolds
  number, Re_f = (the speed jump from the surface to the kink) times (the kink's depth) over (1.225 times the kinematic
  viscosity of the cell just below the kink), is below `FLOOR_RE` = 1.0 becomes the wall: the material below it is taken to
  act as a no-slip floor on the wave's time scale, because its viscous diffusion outruns the wave's advection. The test
  assumes a Newtonian slurry; a shear-thinning or yielding slurry would make a weaker floor. Without such a kink the wall
  is the rigid front.
- *The kink rests on the bridge's 906.4 K end point.* In the thick windward columns of the 50 mm flight the steep top
  layer is not liquid: it is the rheology's bridge, the 906.4 to 908 K band in which the viscosity is interpolated, without
  data, from Li et al.'s slurry law down to the liquid's. On the d050 state at 178.5 s (one state), the first kink lies at a
  median 0.936 times the depth of the 906.4 K isotherm on the 28 columns of the committed fixture that have one, and 0.992
  times it on the 558 of the 622 thick columns that have one; the wall lies at a median 1.004 times that depth on the
  fixture (30 columns, 28 of them floor-like) and 1.03 times it on the thick set (609 columns). So where the kink and the
  wall sit, and every input to the floor test, follow the bridge's unmeasured end point. The records carry the
  kink's depth and temperature so that this stays visible.
- *The constants are assumed* until stage 6: `KINK_DROP` = 0.2, `FLOOR_RE` = 1.0, `FIT_TOL` = 0.03 (every fitted cell
  within 3 % of the top cell's speed of the fit) and `N_MAX` = 6 breakpoints at most. Unchanged from `table`: the droplet
  radius over the wavelength `K_R` = 0.17, the release period over the growth time `K_T` = 1.1 and the surface tension of
  0.80 N/m (the liquid's and the exposed slurry's).
- The solver is inviscid and neglects the gas density and any surface tension at the kinks. The d050 thick columns' top
  layer has a Reynolds number of order 1, where an inviscid solver is least trustworthy; stage 6 is the check.

**Records and history.** On `profile` the droplet records gain nine columns after the 27 of `table` (`spray.PROFILE_COLUMNS`:
the depth to the wall, whether it is a floor, the floor Reynolds number, the first kink's depth and temperature, the
gradient mode's release depth and growth rate, and Girin's theoretical layer depth and surface speed beside them), and
the history gains seven per macro step (`melt_profile.DISPERSION_COLUMNS`: the share of resolved columns that are
floor-like, the count that are fluid-like, the medians of the depth to the wall, of the floor Reynolds number and of the
recorded Weber number, and the mass sprayed by the gradient and the gas-driven modes). On `table` neither is written.
Every spraying skin run's JSON records the setting (`settings.skin.dispersion`; a resume's, `settings.dispersion`).

**What `table` keeps (verified).** On 2026-11-12 the three cases of `analysis/skyscrapers_port/same_results.py` were run
with this change and compared with the same cases run before the option existed (the branch's baseline): a 0.125 s
spraying skin flight run with `--dispersion table`, the same flight on the element model, and the 22.6 s resume of the
saved 100 mm flight state (`d100__physics__skin.state.pkl`, saved before the option existed, so it resumes on `table`).
Of the 22 files compared on each side, the only differences were the new `dispersion` key, recorded as `table`, in the
skin flight's JSON (`settings.skin.dispersion`) and in the resume's (`settings.dispersion`). Every history column, every
droplet record and the whole of the element model's output were identical to the last bit (the comparison leaves out
only the wall-clock columns `wall_*` and `runtime_s`). `--dispersion table` is therefore Girin's theory exactly as it
was.

**Cost (measured once, 2026-11-12).** The spray plan is evaluated once per macro step, at its start.
`analysis/kinked_dispersion/plan_cost.py` repeated it five times on each setting on the saved 50 mm flight state
`d050__handdown-off.state.pkl` (178.5 s, 2367 columns, 909 of them resolved on `profile`, 752 of those floor-like): it took
a median 0.025 s on `table` and 0.815 s on `profile`, so the profile path adds about 0.79 s to a macro step. That flight's
macro steps from 177.5 to 179.5 s (four steps of 0.5 s) took a median 3480 s of wall time in the skins, the 3D solve and
the mesh's moves and rebuilds, so the profile path adds about 0.023 % of a macro step, far inside the design's budget of
5 % (about 174 s). The measurement used the state's own last surface-flow evaluation, which was saved before the slip
band's fields existed, so no column had a layer depth and neither path released anything in it; the fit, the floor test
and the solver, which are the profile path's added work, ran on every resolved column all the same. The cost of finding
the release run in every flow sub-step (the cascade's threshold rule) is not in this figure; the comparison flights'
wall times will show it.

**Observed, one state, not yet judged.** On the 100 mm test sphere at 65 km (one 0.25 s macro step, measured 2026-11-12 in
Task 7) the profile path released 11.0 g against the table's 1.54 g, about 7.2 times as much, almost all by the gas-driven
mode on floor-like columns (497 columns, against 1 by the gradient mode, of 919 resolved columns, 78.6 % of them
floor-like). The reviewer found this to be the design's physics, not a defect: the gas-driven mode has no Weber-number
threshold, and the heat-limited release takes every cell that warms to the mode's threshold temperature during the step.
The comparison flights of stage 5 measure it on the d050 state.
