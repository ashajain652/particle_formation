# `reentry_model` — assumptions and methods (Steps 1–2)

State of the model on `main` as of 2026-09-20 (Step 1: trajectory, spec `superpowers/specs/2026-09-17-reentry-trajectory-model-design.md`;
Step 2: coupled 3D heat transfer, spec `superpowers/specs/2026-09-18-thermal-fem-design.md`). Items marked **(verified)** were
measured against SESAM or an analytic solution; the numbers are in the README's verification tables. The SESAM facts the
model relies on are in `Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md`.

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
