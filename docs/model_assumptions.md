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

## 9. Level-set architecture: geometry (stage 1)

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

## 10. Level-set architecture: the thermal solve (stage 2)

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
`prototype/levelset/stage2/` (committed in 7ac46c5), against read-only copies of the skyscrapers build's Scheil material
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
