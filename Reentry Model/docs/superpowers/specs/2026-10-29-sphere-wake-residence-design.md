# Sphere wake residence: how long a released solid sphere stays in the parent's reentry plume — design

Date: 2026-10-29. Status: approved by Asha in conversation ("Looks good. Proceed ..."); this document is the
written form. Package to be built: `reentry_model/wake_residence/`.

## 1. Intent

When a satellite breaks up during reentry, solid fragments are released into the gas trail the parent
leaves behind. Whether the aluminium of a fragment is released into that hot, chemically active trail, or
into ambient air, depends on how long the fragment stays inside the trail before its own trajectory carries it
out. This model answers one question:

> For a solid AA7075 sphere of diameter D released from the Generic Spacecraft at altitude h, how long is the
> sphere inside the parent's reentry plume before its own trajectory takes it out, and how does that depend on
> D and h?

It also produces the size of the plume itself (radius against distance and time behind the parent), using the
NASA equation the Chemical Reaction Networks / Al agent found.

The overlap computation is a **general engine**, not a one-off script (Asha, 2026-10-29). Its inputs are any
parent trajectory with a size, and any follower trajectory; neither has to come from DRAMA, from a sphere run,
or from this task's sweep. The first uses are (a) a solid AA7075 sphere released from the Generic Spacecraft
(this task's sweep) and (b) a **molten AA7075 droplet released from a solid AA7075 sphere**, which is the parent
whose plume it sits in. Use (b) is built as a worked example and a test (§4.0, §6), so that generality is
demonstrated rather than claimed.

## 2. Goals, non-goals, anti-goals

Goals:
1. The plume radius of the Generic Spacecraft as a function of distance behind it and time since it passed,
   at the parent's state at each release altitude (§4.2).
2. For every (diameter, release altitude) pair of a sweep: the **residence time**, the time from release to the
   first moment the sphere is outside the plume, plus the sphere's altitude, speed, lag behind the parent and
   distance to the parent's track at that moment (§4.4).
3. Time histories for each case: parent and sphere positions, the sphere's lag behind the parent, its distance
   from the parent's track, and the plume radius there, so the exit can be seen and checked.
4. The crossover diameter above which a sphere never falls behind the parent (zero residence), per altitude.
5. A statement, per case, of whether the result depends on the plume equation outside its stated range (§4.2).
6. Plots: plume radius against distance and time; parent and sphere trajectories for chosen cases;
   distance-to-track against plume radius for one altitude across diameters; a residence-time map over
   diameter and altitude.

Non-goals (first version):
- Heating, melting, ablation or mass loss of the sphere. The sphere is a solid of constant diameter and mass.
  Real AA7075 spheres of 1 to 50 mm melt during reentry (see the SESAM sweep in the repository root), which
  changes their drag; that coupling is future work and is named in the report as the main limitation.
- The plume's bulk drift: wind, and buoyant rise or sinking of the hot gas. (The wake gas's own forward
  velocity and low density, and their effect on the follower's drag, are in scope as a switchable effect, §4.3a.)
  The plume boundary is the turbulent-core law (or the laminar law, §4.2), not the shock-heated inviscid outer
  wake, although the outer wake's size enters the energy budget of §4.3a.
- Integrating the sideways forces of shock-wave interaction close to the parent. Instead, near-wake trapping and
  shock-riding are classified and flagged with Register et al.'s (2020) fitted model, not integrated (§4.3b).
- Ejection velocity at separation. The sphere is released with the parent's exact position and velocity (a
  deliberate baseline). The code takes the release state as an argument so a later fragment model can pass a
  separation velocity.
- Parent breakup or mass loss. The parent is the intact DRAMA trajectory of "Main Body" from the committed
  snapshot, whatever its own ablation would do.
- Chemistry in the plume (that is `reentry_model/wake_column`).

Anti-goals — the result is a failure, even if it runs, if:
- the separation between sphere and parent is dominated by a mismatch between how DRAMA and this repository's
  model compute drag and density, rather than by the sphere's own physics, and the report does not say so (§4.3);
- a number is reported without saying whether the plume equation was used inside its range;
- the plume equation's coefficient and exponent, which are reconstructed from a poor scan, are hidden in the
  code instead of being named, documented parameters that a sensitivity run can change.

## 3. Constraints and conventions

- Python in the Reentry Model project, run with `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` from
  inside `Reentry Model/`. New package `reentry_model/wake_residence/` with `analysis/` and `tests/`; add the
  tests path to `pytest.ini`, as `droplet_reentry` and `wake_column` do. Outputs go to
  `output/wake_residence_output/` (git-ignored), exposed as `wake_residence.OUTPUT`.
- Reuse, do not rewrite: `reentry_model.trajectory.Simulator` (rotating-Earth 3-DOF integration in Earth-fixed
  Cartesian coordinates with J2 gravity, so orbital-speed curvature is handled correctly — a flat-Earth
  integrator such as `particle_trajectory.integrate` is NOT acceptable here because at 7.7 km/s the
  centripetal term nearly cancels gravity), `body.ConstantBody`, `aero.SphereDragTables` and
  `aero.bridging_by_name`, `atmosphere.US76TableAtmosphere`, `droplet_reentry.parent_trajectory.load_parent`
  and `droplet_reentry.model.parent_release` (the DRAMA parent's state at an altitude as an `InitialState`;
  used only to start the twin at the top of the branch, §4.3), and
  `earth` for geodetic/Earth-fixed conversions. Build the sphere run the way `reentry_model/cli.py` builds a
  `--thermal none` run (line 211 onward), including its default bridging.
- The three working-tree files modified by the oxidation work (`CLAUDE.md`, `droplet_reentry/model.py`,
  `droplet_reentry/tests/test_droplet_reentry_model.py`) and the untracked oxidation scripts and
  `Literature Review/` are not part of this task. Do not edit `model.py`; do not stage or commit those files.
- Material: AA7075, density 2813 kg/m3 (the value in the repository's `droplet_reentry` and DRAMA material
  JSON `data/user_materials`; the implementer confirms it and cites the file).
- Every run name encodes its full configuration (diameter, release altitude, plume law and any non-default
  plume or parent setting, so sensitivity runs never overwrite baseline runs). Summary JSON
  records the git commit, whether the tree was dirty, and the SHA-256 of the parent file.
- Verification rule of this repository: a claim marked "(verified)" in docs cites a test that asserts it.

## 4. The model

### 4.0 The overlap engine and its interface (module `wake_residence/overlap.py`)

The engine is pure geometry plus the plume law. It must not import `trajectory`, `droplet_reentry`,
`parent_trajectory`, DRAMA or anything that knows what a sphere or a droplet is; a test asserts this by
inspecting the module's imports (§6). Everything else in the package is an adapter or a scenario built on it.

Inputs, all in one common Earth-fixed Cartesian frame, metres and seconds:

- `Track`: a time-ordered path. Built from arrays (`t`, `r` of shape N by 3, optional `v` of shape N by 3, optional
  air-relative `speed` of shape N) or from a callable `position(t)` valid on a stated time range. Positions
  are interpolated (cubic Hermite when velocities are given, otherwise a cubic spline) and the engine
  reports an interpolation-error estimate (compare against a track subsampled by two) so a coarse sampling
  cannot silently give a wrong exit time. Parent and follower use the same class.
- The **parent** is a `Track` plus a size and a plume law: `r0` (metres, or a function of the emission time
  tau, because an ablating parent shrinks), per-sample flow-state arrays (ambient density, temperature,
  pressure, viscosity, mean free path; parent speed, diameter and C_D A; §4.2), and a regime-switching
  `PlumeLaw` (rarefied / laminar Feldman / transitional range / turbulent Taylor or Park, chosen per emission
  point from Kn, Re_d and Taylor's pressure parameter with overridable thresholds, §4.2). The regime functions take flow quantities only; they have no altitude
  argument and `overlap.py` contains no altitude constants. The parent's air-relative
  speed along its track is required (it is the path length integrand); if only positions are given, the
  speed defaults to the speed of the track in the Earth-fixed frame.
- The **follower** is a `Track`, plus optional `t_start` and `t_end` for the window to examine.

Outputs: the plume radius function; a `wake_state(position, t, ambient)` query returning, for a point, whether
it is inside, the margin, the regime, and the wake gas velocity, temperature and density (§4.3a; `ambient` is
the local ambient pressure, temperature and density passed in as numbers, so the engine still imports no
atmosphere module); the margin function (§4.4); the list of **all** intervals during which the
follower is inside the plume (entry time, exit time), not only the first exit; the first-interval residence time
measured from the follower's start; total time inside; the state at each exit; interpolation-error estimate;
and a validity record for each interval: whether s' stayed within `s_park_limit`, whether the parent's speed at
emission lay in Park's 5.8 to 7.9 km/s range, and, when the parent track carries a Knudsen number, the
parent's Knudsen number at emission. These are recorded and flagged, never used to refuse a computation.
A follower that starts outside the plume is a valid input (zero residence, with a later entry possible).

Adapters (each a small function returning a `Track`, kept out of the engine, in `wake_residence/adapters.py`):
the DRAMA parent file's descending branch (the `ParentTrajectory` of `parent_trajectory.load_parent`); a
`trajectory.Simulator` `History` (ECEF states, any sphere run); a `DropletSimulator` run (the rows from
`droplet_reentry.model.run_case`, or its simulator's dense output — whichever resolves the exit to the
required accuracy; a small droplet decelerates within a second or two, so one-second rows are too coarse and
the adapter must use the dense output or a finer cadence); a synthetic analytic path (used by tests).

Scenarios (in `wake_residence/scenarios.py`, each a function using the engine and adapters):
1. `sphere_in_spacecraft_plume(d_mm, h_km)`: this task (§4.5); the parent is the twin of §4.3.
2. `droplet_in_sphere_plume(sphere_d_mm, h_km, droplet_d_um)`: a solid AA7075 sphere of diameter
   `sphere_d_mm` released from the Generic Spacecraft at altitude `h_km` is the parent (its `Track` is its
   `trajectory.Simulator` history, its `r0` is half its diameter, its speed is its air-relative speed, and its
   C_D A comes from its own drag tables); a molten droplet of `droplet_d_um` released from that sphere at the
   sphere's state at that altitude is the follower. Because `run_case` builds its own atmosphere, the scenario
   constructs `DropletSimulator(release, d0_um, atmosphere=WakeAtmosphere(...))` directly and repeats
   `run_case`'s Weber-number pre-check (a droplet too large to survive release is reported as a skipped case,
   not an error); `model.py` is not modified. The droplet's own heating and evaporation are the existing
   model's. This scenario is a worked example and a test; its physical
   applicability is limited and is stated in the output (the plume law is a metre-scale turbulent-wake result
   and the parent here is millimetres to centimetres, so the parent Knudsen-number record matters).

### 4.1 Parent

The parent is `Main Body` in `data/generic_satellite_reentry/dmf_output.json` (the 2026-09-13 DRAMA run of
"Generic Satellite Reentry"; the DRAMA project is labelled "3U Cubesat" but its input gives 800 kg total and
740 kg for the main body). Its descending branch starts near 111 km at about 7.9 km/s and a flight-path angle of
-0.22 degrees; it is still at 7.6 km/s at 80 km and 4.3 km/s at 60 km.

**The parent ablates in DRAMA** (the `..._AeroThermalHistory.txt` table in the same file): its mass stays at
740 kg down to about 81 km, then falls to about 555 kg at 75 km and 216 kg at 66 km, and it lands as a 141 kg
body with a cross-section of 6.215 m2 (`PySara.ImpactingFragments.xml`). Its wall thickness falls from 10.7 mm
to 3.1 mm over the same range: DRAMA treats it as a hollow shell that thins while keeping its outer size. That
is why its effective ballistic coefficient (mass over drag coefficient times area) is about 132 to 140 kg/m2
from 100 to 80 km and then drops to about 67 at 70 km, 28 at 60 km and 25 at 50 km (builder check, computed
as dynamic pressure over deceleration from the file's `dynamicPressure` and `loadFactor` columns).

**Decided with Asha (2026-10-29): the DRAMA file calibrates the parent; the geometry comes from the twin
parent (§4.3).** The file's rows are 5 s apart (about 40 km of flight) and its latitude and longitude are
rounded to 0.001 degree (about 55 to 110 m), which is coarser than the plume (about 1.4 m in radius at the
parent and tens of metres at 30 km). Its positions therefore cannot place a sphere relative to the plume. The
file supplies the parent's ballistic coefficient against altitude, the start state of the twin, and the
consistency check; the twin's integrated trajectory is the parent track `P(tau)` used by the engine.

The plume frame is the Earth-fixed (rotating) frame, i.e. the air with no wind.

### 4.2 Plume equation (Park 1971, NASA TM X-62,052, as found by the Chemical Reaction Networks / Al agent)

The source notes are
`Chemical Reaction Networks/Al/research_notes/Reentry trail conditions at 60 km/wake_flowfield_and_cooling.md`
(sections 3 and 4; paper: https://ntrs.nasa.gov/api/citations/19720024212/downloads/19720024212.pdf):

    r_b / r_0 = C * (s' / r_0)^n,   C = 0.7,  n = 1/3

- r_b: plume (turbulent wake) radius at a point on the track, metres.
- r_0: radius of a circle with the parent's frontal area A, r_0 = sqrt(A / pi).
- s': distance along the track from a virtual origin to the point being evaluated, measured from behind the
  parent. With s the path length between the parent's present position and the track point, s' = s + s_0,
  where s_0 = (1/C)^(1/n) * r_0 = 2.915 r_0 is the offset that makes r_b = r_0 at s = 0 (Park's worked case
  has r_b = r_0 at s' = 73 ft with r_0 = 25 ft, which is the same ratio 2.915).
- Validity: Park states the law holds until the wake's mean velocity falls to atmospheric-turbulence level, about
  9 m/s, at s of order 30 km (about 4 s behind a 7.9 km/s vehicle, by the notes' reading of the scan); beyond
  that growth follows the one-half power of time. Park's finite-rate case covered 5.8 to 7.9 km/s at 70 km. The
  notes flag the 0.7 coefficient and the one-third exponent as reconstructed from garbled text. The code must
  define `S_PARK_LIMIT = 30e3` m as a named constant, flag every case whose exit point or whose still-inside
  end state is beyond it, and expose C and n as arguments so the sensitivity run (§4.6) can vary them.
- Parent frontal area A (decided with Asha, 2026-10-29): constant **6.215 m2**, DRAMA's cross-section of the
  main body (`PySara.ImpactingFragments.xml`), so r_0 = 1.407 m. Constant because DRAMA keeps the body's outer
  size while its shell thins. (The 0.767 m2 `crossSectionArea` in the DRAMA input belongs to a different,
  2.26 kg fragment and is not used.) r_b scales with the square root of A; the sensitivity in §4.6 halves and
  doubles A.
- Speed range (decided with Asha, 2026-10-29): releases at 60 and 50 km (parent at about 4.3 and 1.6 km/s) are
  kept and every result carries a flag when the parent's speed at emission is outside the active law's fitted
  range (Taylor: 3.7 to 6.2 km/s; Park: 5.8 to 7.9 km/s), and when the distance behind the parent is outside it
  (Taylor: 1200 (C_D A)^(1/2), about 2.8 km for the Generic Spacecraft; Park: 30 km). The limits are named
  settings and are not rescaled with speed.
- The plume law is pluggable (`PlumeLaw` in §4.0 is an interface). **Decided with Asha (2026-10-29): the
  turbulent baseline is Taylor, Melcher and Washburn (1964)**, AIAA J. 2(10) 1731, equations 4 and 5, read in full
  (`Literature Review/Wake Growth/`; notes in `Particle Wake Evolution/citations/wake_radius_growth_laws.md`, section 0):

      w / (C_D A)^(1/2) = K [x / (C_D A)^(1/2)]^(1/3),   K = 0.66,   plume radius r_b = w / 2

  where w is the full width of the turbulent core, x the distance behind the body, and C_D A the parent's drag
  coefficient times frontal area, taken from the twin as m / beta (about 5.5 m2 for the Generic Spacecraft; for
  a sphere parent, its own drag coefficient times its area). Fit range: 3.7 to 6.2 km/s, 2 <= x/(C_D A)^(1/2)
  <= 1200, root-mean-square scatter 23 percent; beyond 1200 the one-third law is an extrapolation supported by
  Lees (1964, equation 1, "regardless of the body shape or the Reynolds number"). Near the body (x below about
  2 (C_D A)^(1/2)) the radius is not allowed to fall below r_0. Park's law (above) is kept as a second
  implementation, reported alongside as the wider bound. K, the exponent and both laws' constants are
  arguments, for the sensitivity run.
- **Regime selection (Asha's decisions, 2026-10-29, relayed from the side session that studied the laminar
  law).** The law is chosen **per emission point from flow-regime criteria evaluated from the parent's own state
  at that point, never from altitude bands**: "The code should decide at what point the flight changes
  thresholds, rather than hard-coding in altitude bands, because I want the script to work for an arbitrary
  parent-follower pair." In order:
  1. Knudsen number Kn = mean free path / parent diameter. If Kn is above `kn_continuum` (default 0.01; about
     0.01 to 0.1 is the transitional band) the point is **rarefied**: no wake law applies, the interval is
     flagged "rarefied, no law", and no residence is produced from it unless a fallback law is chosen
     explicitly (setting `rarefied_fallback`, default none).
  2. Otherwise two continuum criteria are evaluated: the Reynolds number Re_d = rho_inf V d / mu_inf
     (freestream, parent diameter; threshold `re_crit`, default 5 x 10^4, Lees 1964 p. 418) and Taylor's pressure
     parameter P = p_inf (C_D A)^(1/2) in cm Hg-cm (thresholds `p_lam`, default 0.5, and `p_turb`, default 2.0,
     Taylor 1964 p. 1735). (Corrected 2026-10-29 by the side session: a single Re_d threshold can never select a
     transitional regime.)
     - **laminar** if Re_d < `re_crit` **and** P < `p_lam`: Feldman law, two bounds;
     - **turbulent** if Re_d >= `re_crit` **and** P > `p_turb`: Taylor (baseline), Park (alternative);
     - **transitional** in every other case, including when the two criteria disagree: the result is a range,
       the minimum and maximum over the two Feldman bounds and the turbulent result at that emission point.
  All four thresholds (`kn_continuum`, `re_crit`, `p_lam`, `p_turb`) are named, overridable settings; nothing is
  altitude-based. The validity record of each interval stores Kn, Re_d, P, the regime, and which law or laws
  produced the radius or range; the speed-range and distance-range flags stay as diagnostics.
  For the Generic Spacecraft (equal-area diameter 2.81 m, C_D A about 5.5 m2) the criteria give roughly:
  rarefied above about 95 to 100 km, laminar about 75 to 95 km, transitional about 63 to 75 km, turbulent below
  about 63 km. **These are outputs of the criteria, never inputs.** Small sphere parents become continuum and
  turbulent much lower (a 5 mm sphere has Kn about 0.3 at 80 km).
- **Laminar law: Feldman (1961)**, J. Aerospace Sci. 28:433, equations 6, 8, 11 (p. 442) and 14 (p. 443), assembled
  by the side session (the combined form is not printed in the paper; notes in
  `Particle Wake Evolution/citations/laminar_wake_growth_laws.md`, section 1):

      R^2 = R0^2 + 4 kappa x / (f U),   kappa = mu(T_axis) / (rho_inf Pr)

  R is the 1/e radius of the temperature (enthalpy) excess in Feldman's density-weighted (Howarth) coordinate;
  R0 is about the body radius; f, about 0.75 to 0.8, is the hot core's speed as a fraction of the flight speed U;
  Pr = 0.71; mu from Sutherland's law. The diffusivity is evaluated at two axis temperatures, giving **two
  bounds**: at ambient temperature (lower bound; Feldman p. 441 says the ambient-diffusivity solution brackets
  the answer from below) and at a hot axis temperature `T_axis_hot` (named setting, default 3000 K; upper
  bound). Laminar residence times are reported as this range; **no cap on the laminar radius is applied**
  (Asha has not chosen one) and both bounds are always reported. Feldman's cases reach only about 76 km and
  hemisphere-cylinders of his sizes; use above that is extrapolation and is flagged. Two points the builder
  must document, not decide silently: R is a density-weighted radius and the physical radius is larger where
  the gas is hot (the first version uses R as r_b, a call recorded in §8); and R is a 1/e radius of the
  temperature excess, whereas Taylor's r_b is the edge of the turbulent core, so the "inside" test means
  slightly different things in the two regimes (a named factor `laminar_edge_factor`, default 1, multiplies R).
- **What the parent must carry.** To evaluate these criteria at every emission time tau, the parent `Track`
  carries, per sample: ambient density, temperature, pressure, viscosity and mean free path, the parent's speed,
  diameter and C_D A. Adapters fill these from an atmosphere; the engine only reads them.

### 4.3 Sphere and the DRAMA-versus-repository consistency check

**Twin parent.** The parent is integrated by the same rotating-Earth integrator as the sphere
(`trajectory.Simulator`, J2 gravity, US76 air, no wind), with its drag acceleration set to
-(1/2) rho V v_rel / beta(h), where beta(h) is the parent's effective ballistic coefficient tabulated against
altitude from the DRAMA file (dynamic pressure divided by deceleration, using `dynamicPressure` and
`loadFactor`; the builder decides the smoothing and documents it). This needs a small `Simulator` subclass
that overrides the drag in `aero_state`; it must not change `trajectory.py`. The twin is integrated **once**,
from the DRAMA state at the top of the descending branch (about 111 km), down to the lowest release altitude
plus a margin, with dense output. Its trajectory is the parent track for every sphere.

**Sphere.** Integrated with `trajectory.Simulator`, a `ConstantBody` of mass (2813 * pi * D^3 / 6), the
SphereDragTables with the default bridging, and the same US76 air with no wind, wrapped in `WakeAtmosphere`
(§4.3a) so that with `wake_flow` other than `off` it feels the wake gas. It starts from **the twin's
state at release altitude h** (position and velocity from the twin's dense output, converted to an
`InitialState`), so sphere and parent track coincide exactly at release. Before the release time the twin's
track supplies the plume behind the parent (`tau_pre` in §4.4).

**Consistency check against DRAMA.** The report states how closely the twin reproduces the DRAMA parent: speed
against altitude (target: within 1 percent above 60 km; the builder reports the actual figure), and altitude
and downrange against time, compared with the file's values within their rounding. If the twin cannot meet the
target, that is reported rather than hidden. The tabulated beta(h) is written out next to each sphere's
ballistic coefficient; it tells which sphere sizes can lag the parent at all, since a sphere with a larger
ballistic coefficient than the parent decelerates less and moves ahead. The crossover diameter therefore
changes with altitude (expected near 100 to 200 mm at 80 to 90 km and near 10 to 15 mm at 60 km) and is an
output, not an assumption.

### 4.3a Wake gas state and its effect on the follower (switchable; Asha, 2026-10-29)

The gas inside the plume is not still ambient air. In the frame of the undisturbed air it **moves forward**,
in the parent's direction of travel (the parent's drag pushes it forward; seen from the parent this is the
familiar velocity deficit), and it is **hot and therefore thin**, because the parent's drag work heats it at
roughly ambient pressure. For a follower inside the plume both effects lower its drag: the relative speed is
the follower's velocity minus the wake gas velocity (a tailwind, not a headwind), and the density is lower.
A follower in the plume therefore separates more slowly than in open air, which lengthens residence; the
effect feeds itself. Rough estimates (Claude's, assumption-laden, for orientation only): mean wake gas speed of
about 900 m/s at 100 m behind the Generic Spacecraft, 200 m/s at 1 km and 20 m/s at 30 km (Park gives 9 m/s at
30 km for his vehicle); a hot-wake density several times below ambient within the first few hundred metres.

**Model (first version).** At a track point emitted at time tau and lying at distance s behind the parent:
- Momentum: the parent's drag D(tau) = (1/2) rho_inf U^2 (C_D A) deposits forward momentum D/U per unit length
  of track. Spread over the core of radius r_b with a smooth radial profile, the core-mean forward speed is
  u_bar = D / (U rho_w pi r_b^2), directed along the parent's velocity at emission, and capped at the parent's
  speed U.
- Energy: the drag work D per unit length is deposited as enthalpy at ambient pressure p_inf. Within the
  turbulent core alone this is not consistent close to the parent (the core is too small to hold the energy at
  ambient pressure: the required area is about 0.2 M^2 (C_D A), radius about 16 m at Mach 27); most of the
  energy is in the shock-heated outer wake. The first version therefore spreads the energy over a radius
  R_E = max(r_b, k_E (D / p_inf)^(1/2)), the larger of the core and the line-blast radius (ReVelle), with k_E a
  named parameter (default 1), and takes the core temperature T_w as this well-mixed mean, from
  rho_w c_p (T_w - T_inf) pi R_E^2 = D with rho_w = p_inf / (R T_w), capped at the stagnation temperature
  T_inf + U^2 / (2 c_p). This understates the centreline temperature; it is documented as assumed.
- Profiles: velocity excess and temperature excess follow a smooth radial profile (for example Gaussian)
  scaled so its area integral equals the core-mean value times pi r_b^2 (momentum) or pi R_E^2 (energy). The
  profiles must be smooth in position, because the follower's integrator is adaptive and a step change at the
  plume edge stalls it. The "inside the plume" test of §4.4 is unchanged (r <= r_b).
- Regime: every case is run twice, with `bound="low"` and `bound="high"`. At each emission point the radius
  is the smallest (low) or largest (high) of the laws its regime allows: in a laminar region the two Feldman
  bounds; in a transitional region the two Feldman bounds and the turbulent law; in a turbulent region the
  turbulent law alone (so low equals high). The residence is reported as the range between the two runs. The
  same conservation statements apply to whichever radius is in force. In a rarefied region there is no wake law and
  no wake gas state: the follower sees ambient air there.

**Switch.** `wake_flow` takes the values `off` (ambient air everywhere: the original definition),
`velocity`, `density` (temperature and density only), and `both`. The sweep runs every case with `off` and
`both`; the sensitivity set adds `velocity` and `density`, so each effect's share of the change in residence
time can be read off. The summary records the residence time for each setting run.

**How it reaches the follower.** Through a `WakeAtmosphere` adapter (in `adapters.py`, not in the engine)
that wraps the ordinary atmosphere: its `state(t, h, lat, lon)` returns the ambient `Freestream` outside the
plume and, inside it, the wake's density, temperature, pressure (ambient) and the wake gas velocity added to the
wind (converted to east-north-up). `trajectory.Simulator` and `DropletSimulator` both read the air only through
this call, so **neither `trajectory.py` nor `droplet_reentry/model.py` is modified**. The wrapper holds the
offset between the follower's clock and the parent's. Performance: the nearest-track-point search runs at
every right-hand-side evaluation, so the engine provides a fast lookup (for example a k-d tree on the
densely sampled twin track, with local refinement), and the plan measures the cost of one sphere run with
`wake_flow=both`.

**Effects on the follower other than drag.** The droplet model's heating and evaporation read the same
temperature and density, so a droplet in the plume sees hot gas. That is intended; whether those heating
correlations hold in hot, recirculating wake gas is not verified and is flagged. Sideways forces from
shock-wave interaction within a few parent diameters (reported for pairs of spheres in supersonic flow, for
example Laurence, Parziale and Deiterding 2012, J. Fluid Mech., not yet read) are a non-goal of this
version and are named as a limitation.

**Sanity checks (assumed values, not measurements):** the model's mean wake gas speed at 30 km behind a body
in Park's case should be of order 10 m/s; the mean temperature a few hundred metres behind the Generic
Spacecraft should be of order 1,000 to 2,000 K (the earlier research notes' estimate from Park's well-mixed
model). Large departures are reported, not tuned away.

**The front plane (controller ruling 2026-10-29, accepted by Asha the same day).** The parent is not an obstacle in
the model, and the union of balls reaches r0 ahead of it. With the wake gas on, a follower released at the parent
felt almost no drag, overtook the parent and left through the front. The plume and its wake gas therefore exist
only behind the parent's current base plane: margin = min(union-of-balls margin, -(x - P(t))·v_hat(t)), and the
wake gas excess is ramped from zero at the plane to full at one r0 behind (cubic smoothstep). An exit through
the plane is labelled: `reached_parent` when the follower had fallen behind and caught up (it would hit the
parent or be pushed aside by its shock, neither modelled), `ahead_of_parent` when it never fell behind (a larger
ballistic coefficient than the parent, or a release on or ahead of the plane).

### 4.3b Release angle and the shock-interaction outcome (Register et al. 2020; Asha, 2026-10-29)

Source: Register, Mathias, Wheeler et al. (2020), Icarus 337:113468, the AFFIx model, read in full (notes:
`Particle Wake Evolution/citations/fragment_shock_separation.md`, section 00; supporting Laurence and Deiterding
2011, Laurence, Parziale and Deiterding 2012, Passey and Melosh 1980, section 0). Close to the parent the outcome
of a release is set by shock-interaction forces, which this model does not integrate. AFFIx classifies it from two
numbers: the size ratio r = D2/D1 (follower diameter over parent diameter) and the release angle theta, the angle
between the parent's velocity axis and the line from the parent's centroid to the release point (theta = 0 directly
behind, 180 at the stagnation point).

**Classifier (new pure module `affix.py`).** `affix_regime(d_follower, d_parent, theta_deg, parent_flow,
follower_density, settings)` returns the final zone and the interaction time:
- `trapped` (near wake, z_f = 3): theta < theta_inn, with theta_inn = 122 r + 47 degrees for r <= r* = 0.353
  (for r > r* the paper's other branch; the builder reads it from the notes or flags it unsupported);
- `far_wake` (z_f = 2): theta_inn <= theta < theta_bif - s_w, theta_bif = 100 exp(-1.618 r) + 40 degrees, shock-band
  half-width s_w a named setting (default 2 degrees, range 1 to 3);
- `shock_band` (z_f = 1): |theta - theta_bif| <= s_w;
- `escape` (z_f = 0): theta > theta_bif + s_w.
Interaction time t_fz = D1 (rho_body / q_inf)^(1/2) t_hat, with t_hat = 0 (trapped), A0(r) + 0.5 exp[-0.05(theta -
theta_inn)], A0 = -7.2 r^2 + 3 r + 0.1 for r <= r* (far wake), 0.88 r + 0.32 + 0.01 exp[0.125(theta - theta_bif - 10)]
(escape); rho_body is the follower's density (AFFIx assumes equal densities; the Generic Spacecraft is a hollow
shell, so an `unequal_density` flag is always set for it, with the parent's bulk density m / (pi D1^3 / 6) recorded).
The function takes flow quantities only (no altitude), accepts any theta in 0 to 180, and returns a validity record:
r inside the fitted range 1/16 to 1 (for D1 = 2.81 m that is D2 >= about 176 mm, so smaller spheres are flagged as
extrapolation), both bodies in continuum flow (Kn = mfp / D below `kn_continuum`), hypersonic (parent Mach at least
`affix_mach_min`, default 5), and the equal-density flag.

**Release state by outcome (Asha, 2026-10-29).** The follower's release velocity is the parent's (zero sideways
release speed stays the baseline; Register's dispersion result supports it). The release POSITION depends on the
outcome, because the engine does not integrate the shock-driven motion that delivers the follower to it:
- `far_wake`: on the wake axis, **two parent diameters behind the parent's centroid** (P(t_r) - 2 D1 v_hat), Asha's
  estimate from Register et al. The far-wake start does not depend on theta, so one integration per (diameter,
  altitude, bound, wake_flow) serves every far-wake theta; theta changes only the classification and t_fz.
- `escape`: at the windward surface point, P(t_r) + (D1/2)(-cos(theta) v_hat + sin(theta) n_hat), with n_hat the
  unit vector perpendicular to v_hat in the local vertical plane, pointing up. Its engine residence is expected to
  be near zero (`ahead_of_parent` or `starts_outside`); t_fz is recorded beside it, and a computed residence much
  longer than t_fz is flagged.
- `shock_band`: not integrated; reported with residence NaN and reason `shock_band` ("outside both models").
- `trapped` (unreachable on the windward side for our sizes, kept for general use): not integrated; reported as
  censored with reason `near_wake_trapped`, rather than trusting a drag-only exit.

**Scope for the Generic Spacecraft.** Releases are on the windward side only: theta from 90 to 180 degrees
(Asha). The spacecraft scenario rejects other values; the general function accepts 0 to 180. For 1 to 400 mm
spheres theta_inn is 47 to 64 degrees, so no windward release is trapped; theta_bif is about 139 degrees (10 mm),
134 (100 mm), 129 (200 mm) and 119 (400 mm).

**Cross-check.** AFFIx reports a near-wake force coefficient of about 0.12 against 0.88 in free flow, a ratio of
about 0.14. The tables report our model's drag ratio for a sphere at the far-wake start point (two parent diameters
behind, on the axis), rho_w (V - u)^2 / (rho_inf V^2) with V the parent's speed, at each release altitude and bound,
next to 0.14.

### 4.3c The follower's own flow regime, the shock-layer gas and the bow shock (Asha, 2026-10-29)

Register's and Laurence's outcomes come from continuum forces (the follower's own bow shock interacting with the
parent's). A follower in free-molecular flow has no shock of its own, barely disturbs the parent's flow, and moves
almost ballistically; its separation is set by drag differences, which the engine integrates. So before any
Register branch is used, the follower's own regime is decided (side-session addendum, approved by Asha with the
changes below).

**Order of decisions.**
1. Parent regime first (`laws.classify` on the parent's state at release). Rarefied or too slow: no plume law and no
   AFFIx; the case is reported without integration (reason `rarefied_no_law` / `too_slow_no_law`).
2. Follower regime from Kn_f = mean free path of the LOCAL gas at the release point / follower diameter, with the
   named thresholds `kn_continuum` (default 0.01) and `kn_free_molecular` (default 10, the repository's existing
   free-molecular gate): `continuum` below the first, `free_molecular` at or above the second, `transitional`
   between. Flow quantities only; no altitude.
3. Continuum follower: AFFIx governs (§4.3b) — far wake on the axis two parent diameters behind, escape at the
   windward surface point, shock band not integrated, trapped censored.
   Transitional follower: AFFIx is computed and labelled `indicative`; the case is flagged "transitional follower:
   neither continuum interaction nor ballistic limit"; the integrated engine result is primary.
   Free-molecular follower: AFFIx is labelled `not_applicable` (values recorded for information); the integrated
   engine result is primary and nothing overrides it.

**The local gas at the release point is post-shock (Asha).** Every release point is on the parent's windward
surface (theta 90 to 180 degrees), inside the parent's shock layer, so Kn_f uses the shock-layer gas there, not the
free stream. Model (new adapter, using `reentry_model.gas.EquilibriumAir`, Cantera equilibrium air): the stagnation
state behind the parent's normal shock (`EquilibriumAir.stagnation(rho_inf, T_inf, U)`); the local surface pressure
by modified Newtonian theory, p(phi) = p_inf + (p_stag - p_inf) cos^2(phi), phi = 180 - theta the angle from the
stagnation point; the local state by isentropic (equilibrium) expansion from the stagnation entropy to p(phi); the
mean free path from that local density and mean molecular mass (`aero.mean_free_path`). Requires the parent in
continuum (step 1). Recorded as assumed: Newtonian pressure under-predicts near the shoulder (theta near 90), and
the hollow-shell parent is treated as a sphere of its equal-area diameter. Kn_f at the exit (or end) uses the local
gas there (the `WakeAtmosphere` state: wake core or ambient), so a follower inside the hot core is classified more
rarefied than in ambient air at the same point.

**Rarefied followers (transitional and free-molecular) start at their surface point (Asha).** Release position
P(t_r) + `release_offset(theta, D1, v_hat, up)`, velocity = the parent's. Because the plume lies behind the
parent's base plane, such a follower usually starts outside it and enters as it falls behind, so its result is
reported as: the time of first entry after release, the duration of that first inside interval, the total time
inside, and the exits (the engine's `ResidenceResult.intervals`), with `residence_time` kept as defined in §4.4.

**Does the rarefied follower cross the parent's bow shock? (Asha; superseded the same day by §4.3d.)** The first
version of this section checked only whether the follower's ballistic trajectory left the bow-shock envelope. Asha
then chose to integrate the rarefied follower through a modelled shock-layer flow (§4.3d), so the crossing is now a
result of that trajectory (`shock_layer_exit`), not a separate check.

**Recorded in every summary:** Kn_f at release (shock-layer gas) and at exit or end, the follower regime, the AFFIx
status (`governing` / `indicative` / `not_applicable`) and the branch that governed, the release kind (`far_wake`,
`escape`, `surface`), and for rarefied followers the entry time, first-interval duration, total time inside and the
shock-layer outcome of §4.3d.

**Droplet note (side session).** In free-molecular flow the hot wake lowers drag more than heating, because the wake
gas carries thermal energy; heating may fall less than drag, or even rise, close behind the parent. Droplet results
from the hot core carry this as a stated limitation.

### 4.3d The shock-layer flow and the rarefied follower's motion through it (Asha, 2026-10-29, option 1)

**Why the follower can be integrated through a fixed field.** A transitional or free-molecular follower has no
shock of its own and barely disturbs the parent's flow, so the coupling runs one way: the parent's shock-layer gas
(the gas between the bow shock and the parent's surface) is computed alone, and the follower is integrated through
it the way it is integrated through the wake. Continuum followers keep §4.3b: their own shock interacts with the
parent's, AFFIx governs, and they do not see this field. Transitional followers use it with the §4.3c flag, since
one-way coupling is only approximate for them. The whole section requires the parent in continuum (§4.3c step 1);
there is no shock layer otherwise.

**Frame and quasi-steadiness.** A parent-centred frame: origin at the parent's centroid P(t), x along v_hat(t)
(positive forward), y the distance from the axis; and spherical coordinates (r, phi) about the centroid, with phi
the angle from the stagnation line (phi = 180 - theta). The parent is the equal-area sphere of radius R = D1/2. The
field is built once per release from the parent's free stream at release (rho_inf, T_inf, p_inf, U, Mach M,
mean free path) and then frozen (a named simplification); it moves with P(t) and turns with v_hat(t) along the twin
track. The field gives the gas velocity relative to the parent, w; the Earth-fixed gas velocity handed to the
follower is U(t) v_hat + w, so the gas at the stagnation point moves with the parent and the gas outside the shock
is at rest. The summary records the relative change of the parent's dynamic pressure between release and the
follower's shock-layer exit, flagged above 10 percent.

**Shock shape.** Billig's (1967) sphere correlation for the shape: vertex curvature radius R_c = R * 1.143
exp(0.54/(M - 1)^1.2) and x_s(y) = R + Delta - R_c cot^2(beta) [ (1 + y^2 tan^2(beta)/R_c^2)^(1/2) - 1 ], beta =
arcsin(1/M); a point is inside the envelope when x < x_s(y). The standoff Delta is a setting `shock_standoff`:
`"lighthill"` (default) takes it from the constant-density solution below at the equilibrium density ratio
epsilon = rho_inf / rho_2 across the normal shock; `"billig"` takes Billig's perfect-gas value R * 0.143
exp(3.24/M^2) (sensitivity). Reason for the default (the author's call): Billig's fit is for a perfect gas with a
ratio of specific heats of 1.4, where epsilon is about 1/6, whereas equilibrium air at 7 to 8 km/s has epsilon of
about 0.06 to 0.1, and the standoff scales roughly with epsilon, so Billig's standoff would be about 2 to 3 times
too large. Billig's correlation is fitted for Mach 2 to about 25; use beyond that is flagged.

**State just behind the shock.** At each shock point, the local shock angle sigma (between the shock surface and
the free-stream direction) follows from the slope of x_s(y): 90 degrees on the axis, falling towards beta far out.
The jump is the equilibrium Rankine-Hugoniot relation applied to the normal velocity component U sin(sigma), with
Cantera equilibrium air by the same fixed-point iteration as `gas.EquilibriumAir.stagnation`; the tangential
component U cos(sigma) is unchanged. This gives the post-shock density, temperature, pressure and gas velocity.

**State at the body surface.** Pressure by modified Newtonian theory, p_b(phi) = p_inf + (p_stag - p_inf)
cos^2(phi) (§4.3c); temperature and density by isentropic equilibrium expansion from the stagnation state to p_b;
gas speed tangent to the surface from the total enthalpy, |w| = (2 (h_0 - h))^(1/2), h_0 the stagnation enthalpy.
Inviscid: no boundary layer (named; at 80 km it is roughly a tenth of the layer's thickness, and thicker higher up).

**Between body and shock.** A layer coordinate eta = (r - R)/(r_s(phi) - R), 0 at the body and 1 at the shock,
where r_s(phi) is the shock's distance from the centroid along the ray at phi. Pressure and temperature vary
linearly in eta between the surface and post-shock values; density is the equilibrium density at that pressure and
temperature. The velocity:
- near the nose (phi up to `phi_nose_deg`, default 40), Lighthill's (1957) constant-density solution, as given by
  Hayes and Probstein. The gas behind a spherical shock of radius R_s concentric with the body is taken as
  incompressible at density rho_inf/epsilon, with vorticity omega = C r sin(phi), uniform C (the property of
  axisymmetric inviscid flow at constant density). The Stokes stream function is psi = sin^2(phi) g(r) with g(r) =
  A r^2 + B / r - (C/10) r^4, velocities u_r = 2 g cos(phi) / r^2 and u_phi = -g' sin(phi) / r. The conditions are:
  g(R_s) = -(1/2) epsilon U R_s^2 (normal velocity reduced by epsilon), g'(R_s) = -U R_s (tangential velocity
  continuous), g(R) = 0 (the body is a streamline), and |C| = U (1 - epsilon)^2 / (epsilon R_s^2), the vorticity
  behind a spherical shock from Crocco's theorem. In this sign convention the r^4 coefficient is
  -(C/10) = -U (1 - epsilon)^2 / (10 epsilon R_s^2) (negative; checked numerically by the author on 2026-10-29: the
  resulting standoff is within 2 to 7 percent of Hayes and Probstein's approximation for epsilon 0.05 to 0.2, while
  the opposite sign misses it by 25 to 40 percent). These give A, B and R_s, and hence the standoff R_s - R. For
  small epsilon g(R) = 0 has a second, spurious root further out; take the root closest to the body.
  The solution is evaluated at r_L = R + eta (R_s - R), i.e. stretched to the local layer thickness;
- towards the shoulder (phi from `phi_nose_deg` + `phi_blend_deg`, default 20, onward), linear interpolation in
  eta between the surface velocity and the post-shock velocity;
- between the two, a cubic-smoothstep blend in phi.

Lighthill's solution is the classical one for the nose and is weakest at the shoulder; Newtonian pressure
under-predicts near the shoulder (theta near 90); most windward releases between 90 and 135 degrees lie there.
Both are named as assumptions.

**Edges and smoothing.** The field applies inside the envelope, outside the body and ahead of the base plane
(x > 0). Across the shock, the state passes from the ambient (or wake) gas to the field over a thickness of
`shock_thickness_mfp` (default 3) free-stream mean free paths, measured along the local shock normal, with a cubic
smoothstep. This is physical (a real shock is a few mean free paths thick) and keeps the adaptive integrator from
stalling at a discontinuity. Behind the base plane, the field fades from full at x = 0 to zero at x = -r0 with a
cubic smoothstep, evaluated there as on the shoulder ray (phi = 90) at the same distance from the axis, while the
wake gas of §4.3a ramps in over the same distance (the existing front-plane ramp). Density, temperature, pressure and
velocity are blended linearly with that weight. At trial points inside the body (the integrator may evaluate there)
the state is the surface state (eta clamped to 0), so it is always finite.

**Tabulation.** The field is computed per release on a grid in (phi, eta) and interpolated with smooth (bicubic)
splines, so the Cantera calls happen once per release, not at every step. The pressure gradient is the derivative
of the pressure spline. Grid spacing is the builder's call, with a test that halving it changes the follower's
shock-layer exit time by less than 1 percent. The cost of building one field and of one follower run is measured.

**Forces on the follower.**
- Drag relative to the local gas, with the composed drag law of `droplet_reentry.particle_trajectory`
  (`drag_coefficient`: Schaaf and Chambre free-molecular, a continuum correlation, Cunningham slip, blended on the
  speed ratio), which holds at any speed ratio; the SESAM tables used by `sphere_simulator` are fitted to hypersonic
  conditions and the relative speed in the shock layer is low. The wall-to-gas temperature ratio uses a named
  follower wall temperature `follower_wall_T_K` (default 300 K, assumed: the sphere's temperature is not modelled).
  At zero relative speed the drag force is zero (C_D times V stays finite); the implementation must not produce NaN
  at release, when the relative speed at the stagnation point is near zero. The rarefied follower uses this law for
  its whole flight, inside and outside the shock layer, so its drag comes from one smooth law; the summary records
  C_D from this law and from the SESAM-table law at the release free stream, as a consistency number. The law's
  air viscosity and mean free path are for undissociated air; in the hot shock layer that is approximate (named).
- A pressure-gradient force, a_p = -grad(p) / rho_f (rho_f the follower's density): the net pressure force on a
  body small compared with the distance over which the pressure changes (about R). Switch `pressure_gradient`
  (default on; off in the sensitivity set). The summary records the largest |a_p| and the largest drag acceleration
  inside the shock layer, so its share is visible.
- Gravity, as now.

**How it reaches the follower.** The gas goes through the existing seam: `WakeAtmosphere` takes an optional
shock-layer field and its `state()` returns the blended gas (shock layer, wake or ambient) with the Earth-fixed gas
velocity added to the wind. The pressure-gradient force is not an air property, so the follower is a
`trajectory.Simulator` subclass in `adapters.py` (as the twin is) whose `aero_state` uses the composed law and adds
a_p, read from a new `WakeAtmosphere.pressure_gradient(t, x)`. `trajectory.py` and `droplet_reentry/model.py` are not
modified. Droplets are not given the shock-layer path in this version (the droplet-in-sphere use releases in the
sphere's far wake); that is a named limitation.

**The parent is an obstacle for rarefied followers.** A terminal event fires when the follower's distance from P(t)
falls below R: reason `hit_parent`, with the time, the impact angle theta and the relative speed. The follower starts
on the surface, so the event counts only an inward crossing after release. (Continuum followers keep the front-plane
labels of §4.3a.)

**Outcome recorded.** `shock_layer_exit`: `through_shock` (left the envelope across the shock into the free stream),
`past_base_plane` (crossed the base plane inside the envelope, i.e. swept round towards the wake), `hit_parent`, or
`inside_at_end`; with its time, x and y in the parent frame, and the relative speed. After `through_shock` or
`past_base_plane` the integration continues: the follower may later hit the parent, enter the plume, or neither; the
residence quantities of §4.3c (first entry, first-interval duration, total inside) are reported as before. Also
recorded: epsilon, the standoff used and Billig's for comparison, the Mach number and the validity flags below.

**Validity flags.** Mach outside Billig's fitted range; the parent's Knudsen number above `kn_shock_layer_flag`
(default 1e-3, assumed), where the shock layer is likely out of chemical equilibrium and thickened by viscosity, so
equilibrium air overstates dissociation and the thin-layer picture weakens; the transitional-follower flag of §4.3c;
the dynamic-pressure change above.

### 4.4 In the plume, and residence time

Let P(tau) be the parent's Earth-fixed track position at time tau, and S(t) the sphere's position at time t,
with release at time t_r (both at the same place at t_r). Let s(tau, t) be the parent's path length between tau
and t, and r_b(s' = s + s_0) the plume radius from §4.2. The sphere is **inside the plume** at time t if

    there exists tau in [t_r - tau_pre, t] with  |S(t) - P(tau)| <= r_b( s(tau, t) + s_0 ).

That is, the plume is the union of balls of radius r_b centred on every point of the track, each ball growing
as the parent moves on. The parent track before release (tau < t_r) is included for completeness
(`tau_pre` = 60 s), although with release on the parent's own state the sphere is at first inside the ball at
tau = t_r. The **margin** is the largest value over tau of r_b - |S - P(tau)|; it is positive inside.

**Residence time** is the first t - t_r at which the margin becomes negative, found by bracketing and root
finding on the margin (not by reading the output grid), with a tolerance of 1 ms. Also report the total time
inside and the number of exits and re-entries (expected: one exit), and for the exit moment: the sphere's
altitude, its speed relative to the air, its along-track lag behind the parent, the distance to the nearest
track point, the plume radius there, and s'. If the sphere is still inside when its integration ends (ground,
or the integration's t_max), the case is reported as censored with the time at which it ended, never as an
exit.

A case with a sphere that moves ahead of the parent has residence time equal to the time it takes the margin
to go negative from the start, which can be very short; it must not be reported as "never in the plume" unless
the margin is negative at every sampled time after t_r.

### 4.5 Sweep

Default sweep: diameters 1, 5, 10, 25, 50, 100, 200, 400 mm; release altitudes 100, 90, 80, 70, 60 and
50 km; the twin parent. (The first six diameters span the repository's SESAM sweep, 5 to 100 mm,
extended down to 1 mm; 200 and 400 mm are added so the crossover diameter, expected between 100 and 200 mm at
80 to 90 km and near 10 to 15 mm at 60 km, is bracketed at every altitude.) The
implementer adjusts the altitude list to the parent's branch if any altitude is outside it, and states so. The
spheres are integrated in the continuum-to-free-molecular range the repository's drag tables cover; any case
whose Knudsen number or Mach number leaves the tables' range is flagged.

Release angles (§4.3b): theta = 90, 120, 135, 145, 160 and 180 degrees, which bracket theta_bif for every size.
Every (diameter, altitude, theta) is classified; far-wake cases share one integration per (diameter, altitude,
bound, wake_flow); escape cases are integrated from their surface points; shock-band cases are not integrated.
Rarefied followers (§4.3c, §4.3d) are integrated per theta from their surface points through the shock-layer field.
(The author's estimate, to be confirmed by the sweep: the follower's local Knudsen number at the surface is about a
tenth of its free-stream value near the stagnation point (larger towards the shoulder), so in the default sweep rarefied followers are transitional, at 80 and 90 km for
the smaller diameters, and none is free-molecular; at 100 km the parent itself is rarefied, so nothing is
integrated there.)

Time step: the sphere history is sampled finely enough to resolve the exit to the 1 ms tolerance (use the
integrator's dense output rather than the 1 s cadence).

### 4.6 Sensitivity

Around six reference cases (spheres of 10, 50 and 200 mm, each at 80 km, where the Generic Spacecraft's
wake is expected to be laminar, and at 60 km, where it is expected to be turbulent; the regime is whatever the
criteria give): vary the turbulent law's coefficient by the factors 0.5/0.7, 1 and 0.9/0.7 (Park's original
0.5/0.7/0.9 range, applied to Taylor's K), its exponent (0.30, 1/3, 0.40), the turbulent law itself (Park
instead of Taylor), the Feldman hot-axis temperature (2000, 3000, 5000 K) and core-speed fraction f (0.75,
0.8), the parent frontal area (3.1, 6.2, 12.4 m2), the parent's ballistic-coefficient table scaled by 0.9 and
1.1, the sphere density (2700, 2813, 2900 kg/m3), and `wake_flow` (`velocity` and `density` alone); and, for a
rarefied follower (10 mm at 90 km, theta 180 and 120 degrees), `shock_standoff="billig"`, `pressure_gradient` off,
`follower_wall_T_K` 1000 K and `shock_thickness_mfp` 1 and 10, with the shock-layer outcome reported beside the
residence. Report
each residence time and which parameter moves it most. (Amended 2026-10-29 when Taylor replaced Park as the
baseline: the original set varied Park's C and n at 80 km only, which the regime criteria now place in the
laminar regime.)

## 5. Outputs

Under `output/wake_residence_output/`:
- `plume_size.csv` and plot: r_b against s and against time since passage (s divided by the parent's speed) for
  the parent's state at each release altitude, with the Park-limit distance marked.
- `ballistic_coefficient.csv`: the parent's effective ballistic coefficient against altitude, and each
  sphere's, so the lag direction is visible.
- `runs/<run_name>.csv` and `.json` per case: history and summary (§3, §4.4).
- `twin_vs_drama.csv`: the twin's speed, altitude and downrange against the DRAMA file's, with the differences.
- `summary.csv`: one row per case — diameter, altitude, residence time, censored flag, exit state,
  Park distance-limit flag, Park speed-range flag.
- `sensitivity.csv`.
- plots listed in §2 goal 6, in `plots/`.

Docs: a new section in `docs/model_assumptions.md` (assumptions with verification status, the plume equation's
provenance and its caveat); the commands and layout in `CLAUDE.md` and the repository `README.md`; the
`CLAUDE.md` edit stays uncommitted alongside the existing oxidation edits to that file, and the final report
tells Asha so.

## 6. Verification the plan must include

- Plume equation: unit test that r_b(s=0) = r_0, that r_b scales as the one-third power, and that C and n
  arguments are honoured.
- Geometry: a synthetic straight track with a sphere moving off it at a known rate; the residence time equals
  the hand-computed value; also a curved track, to confirm the union-of-balls margin handles curvature.
- Sphere identical to the twin parent: residence time is the full time the integration runs inside the plume
  (the sphere never leaves), i.e. the censored path is exercised.
- Sphere integration: the repository's `trajectory.Simulator` run reproduces a stored SESAM reference sphere
  run to the tolerance the repository already asserts for that comparison (the implementer finds and cites
  the existing test and tolerance and re-uses its reference), so the sphere integrator is known good in this
  regime before any wake result is trusted.
- Twin parent: reproduces the DRAMA parent's speed against altitude within 1 percent above 60 km (§4.3); if
  it cannot, the test records the achieved figure and the failure is reported rather than hidden. A sphere
  released from the twin starts exactly on the twin's track (margin at release equals r_0 to round-off).
- Reproducibility: a short sweep run twice gives identical `summary.csv`.

Generality tests (the engine is arbitrary-input, so these are required, not optional):
- Import test: `overlap.py` imports nothing from `trajectory`, `droplet_reentry`, `parent_trajectory`, `body`,
  `aero` or `atmosphere`.
- Same code, different inputs: the engine, unchanged, gives the hand-computed answer for (i) a straight parent
  and a straight follower diverging linearly, (ii) a circular-arc parent, (iii) a follower that starts outside
  and enters later, (iv) a follower that exits and re-enters (two intervals), (v) an `r0` that changes with
  emission time.
- Worked example: `droplet_in_sphere_plume(50, 80, 100)` and one more combination run end to end without
  changing the engine, and produce an exit time, a validity record and a non-negative residence; the test
  asserts that the run completes, the record is populated, and the sphere and droplet tracks begin at the
  same point, and does not assert a physical residence value.
- Interpolation: for a track subsampled by two, the exit time changes by less than the stated tolerance, or
  the engine reports an error estimate larger than it.

Regime tests (§4.2):
- Selection is driven only by Kn, Re_d and P: synthetic parent states that cross each threshold at arbitrary,
  deliberately unrealistic altitudes (including a small-body parent that is rarefied low down) give the
  expected regime; the same Kn, Re_d and P placed at different altitudes give the same regime.
- Disagreement: synthetic states where Re_d and P disagree (Re_d above `re_crit` with P below `p_lam`, and the
  reverse) are classified transitional.
- Law switching along one track: a synthetic parent passing rarefied, laminar and turbulent states in one
  flight produces intervals tagged with the right law, and the radius is continuous or its jump is reported at
  each switch.
- Thresholds are overridable: changing any of `kn_continuum`, `re_crit`, `p_lam` or `p_turb` moves the switch
  point accordingly.
- No altitude constants in `overlap.py` (a source check for the absence of altitude parameters in the regime
  functions, plus the test above).
- Feldman bounds: with T_axis equal to ambient the two bounds coincide; R^2 grows linearly in x with slope
  4 kappa/(f U); the hot bound is never narrower than the ambient bound.
- Transitional range: in a transitional interval the record holds the two Feldman results, the turbulent result,
  and their minimum and maximum.

Follower-regime, shock-layer and bow-shock tests (§4.3c):
- With identical r and theta, a follower state with Kn_f < 0.01 gives AFFIx `governing`, one between the
  thresholds `indicative`, one >= 10 `not_applicable`; placed at arbitrary altitudes to show nothing depends on
  altitude.
- The same follower is classified more rarefied inside a `WakeAtmosphere` hot core than in ambient air at the same
  point.
- A free-molecular follower's reported residence equals the engine's integrated value, unchanged by any AFFIx
  output.
- Shock-layer gas: at theta = 180 the local state equals `EquilibriumAir.stagnation`; the local pressure falls
  towards p_inf as theta falls towards 90; the local density exceeds the free-stream density at the stagnation
  point; Kn_f at the release point is smaller than the free-stream Kn_f.
Shock-layer field and rarefied-follower tests (§4.3d):
- Shock: Billig's vertex radius and (with `shock_standoff="billig"`) standoff reproduce the formula; the shock passes
  through (R + Delta, 0); points on the body surface are inside the envelope.
- Jump: the oblique jump at sigma = 90 degrees followed by isentropic compression to rest reproduces
  `EquilibriumAir.stagnation`'s pressure within 0.5 percent; the tangential velocity is unchanged across the jump; a
  shock angle half a degree above the Mach angle gives a density ratio within 5 percent of 1.
- Lighthill: the three shock and body conditions hold to round-off and the stream-function equation holds at
  interior points; the standoff agrees with the approximation Delta/R_s = epsilon / (1 + (8 epsilon/3)^(1/2)) to
  within 10 percent for epsilon from 0.05 to 0.2.
- Field: at the stagnation point the state equals `EquilibriumAir.stagnation` within 1 percent and w = 0; at the
  body the normal velocity is zero and the pressure is the Newtonian value; at eta = 1 the state is the post-shock
  state; the mass flux through the layer at phi = 30 and 60 degrees equals the free-stream flux through the shock
  disc, pi y_s^2 rho_inf U, within 25 percent; beyond the smoothing thickness outside the shock the state is exactly
  the ambient (or wake) state; the state is continuous across the base plane and the shock.
- Pressure gradient: the spline gradient matches a finite difference of the field; on the axis a_p points away from
  the body.
- Follower: with the field and the pressure gradient switched off, the rarefied-follower simulator gives the same
  trajectory as a plain simulator with the composed law (round-off); no NaN at a release with zero relative speed;
  a free-molecular follower released at the stagnation point of a decelerating parent first moves away from the
  surface (forward relative to the parent); a follower given a velocity towards the parent ends with `hit_parent`;
  halving the grid spacing changes the shock-layer exit time by less than 1 percent.

Release classifier tests (§4.3b):
- `affix_regime` reproduces theta_inn and theta_bif at r = 1/16, 1/4 and 1/2 and the zone boundaries around them.
- theta = 90 with a small r gives `far_wake`; theta = 180 gives `escape`; a theta inside the shock band gives
  `shock_band`; theta below theta_inn gives `trapped`.
- t_fz follows the three t_hat formulas; t_hat = 0 when trapped.
- An input outside 90 to 180 is accepted by the general function but rejected by the spacecraft scenario.
- The validity record flags r < 1/16, Kn > `kn_continuum` for either body, Mach below `affix_mach_min`, and unequal
  densities; the classifier has no altitude argument.
- Release positions: the far-wake start is on the axis two parent diameters behind the centroid; the escape start
  is on the parent's surface at theta (distance D1/2 from the centroid, angle theta from -v_hat).

Wake gas state tests (§4.3a):
- Conservation: integrating the velocity profile over a cross-section gives D/U per unit length, and the
  temperature excess gives D per unit length, to quadrature tolerance.
- Caps: the gas speed never exceeds the parent's speed and the temperature never exceeds the stagnation value.
- Smoothness: the returned state is continuous (and has a continuous first derivative) across the plume edge.
- Switch: with `wake_flow="off"`, `WakeAtmosphere.state` returns exactly the wrapped atmosphere's `Freestream`,
  and a sphere run reproduces the plain `trajectory.Simulator` run to round-off.
- Direction: a sphere released inside the wake with `wake_flow="both"` decelerates less than with `off` over
  its first second.
- Neither `trajectory.py` nor `droplet_reentry/model.py` changes (the plan's final check is `git diff` on both
  against the starting commit, ignoring the oxidation edits already in the working tree).

## 7. What is left to the builder

Module layout inside the package; the interpolation scheme for sampled tracks; the smoothing of the
ballistic-coefficient table; the root-finding method; plot styling (use
`DRAMA Sphere Experiments/analysis/plot_style.py` where the repository already does).

## 8. Decisions and calls

Decided with Asha on 2026-10-29, after the builder check: the twin parent carries an altitude-dependent
ballistic coefficient from the DRAMA file (§4.3); the plume uses a constant frontal area of 6.215 m2 (§4.2);
the twin, not the DRAMA file's rounded positions, is the parent track (§4.1); 50 and 60 km releases are kept and
flagged as outside Park's speed range, with the 30 km limit unscaled (§4.2). Later the same day: Taylor et al.
(1964) is the turbulent baseline and Park the alternative (§4.2); the laminar law waits for another agent's
result (§4.2); the wake gas's forward velocity and low density act on the follower as a switchable effect
(§4.3a). From the laminar-law side session (relayed 2026-10-29 with Asha's approval): regime chosen per
emission point from Kn, then Re_d and Taylor's pressure parameter together (laminar only if both say
laminar, turbulent only if both say turbulent, otherwise transitional and reported as a range), never from
altitude bands, with all four thresholds overridable (§4.2; corrected by the side session the same day); Feldman (1961) is the laminar law, run as two bounds
(ambient and hot axis diffusivity) with no cap on the laminar radius (§4.2); rarefied points have no law
unless a fallback is chosen explicitly (§4.2). Later still: the front-plane ruling (plume only behind the
parent's base plane; `reached_parent` / `ahead_of_parent` labels) stands (§4.3a); the Register et al. (2020) release
classifier is added with releases limited to the windward side, theta 90 to 180 degrees, far-wake releases starting
on the wake axis two parent diameters behind the parent (Asha's estimate from Register et al.), escape releases at
the windward surface point, shock-band cases flagged and not integrated (§4.3b). Then: the follower-regime gate
(addendum) approved with Asha's changes — Kn_f from the local post-shock gas at the release point; transitional and
free-molecular followers start at their surface point at theta and are checked for crossing the parent's bow shock
along their ballistic trajectory (§4.3c). Then (option 1 of three offered: engineering model, tabulated continuum
CFD, DSMC): rarefied followers are integrated through an engineering model of the parent's shock-layer flow, which
replaces the ballistic bow-shock check (§4.3d); CFD and DSMC stay as later validation for a few cases.

Calls made by the author (not decided by Asha; the final handover lists them again):
- Union-of-balls definition of "in the plume" (§4.4).
- Plume frame is the Earth-fixed air with no wind (§4.1), and a no-ejection-velocity baseline (§2).
- The 1 percent twin-versus-DRAMA target (§4.3).
- The wake gas model of §4.3a: momentum and energy conservation, the energy spread over the larger of the core
  and the line-blast radius, well-mixed core temperature, smooth profiles, caps, and the four switch settings.
- Integrating the follower through a `WakeAtmosphere` wrapper rather than modifying the integrators (§4.3a).
- Using Feldman's density-weighted 1/e radius R directly as the laminar plume radius (`laminar_edge_factor` = 1)
  (§4.2).
- Running each case twice, on the pointwise smallest and largest admissible radius (`bound="low"`/`"high"`),
  which realises the relayed "minimum and maximum at that emission point" decision (§4.3a).
- The gas-state profiles: Gaussians whose 1/e radii are r_b (momentum density) and R_E (the energy variable
  1 - T_inf/T), so the centre values equal the spec's well-mixed means (§4.3a).
- A parent slower than `u_min_parent` (default 1 m/s) at emission has no wake law there, flagged (§4.2), so a
  sphere parent settling at terminal speed cannot make Feldman's x/(f U) blow up.
- Defaults `T_axis_hot` = 3000 K and f = 0.8 for the Feldman law (§4.2).
- Default diameter and altitude lists (§4.5) and the sensitivity ranges (§4.6).
- Shock-layer surface state (§4.3c): modified Newtonian surface pressure with isentropic equilibrium expansion
  from the stagnation state; the equal-area sphere for the hollow parent; `kn_free_molecular` default 10.
- Shock-layer field (§4.3d): Billig's shape with the standoff from Lighthill's constant-density solution at the
  equilibrium density ratio (Billig's own standoff as the sensitivity); Lighthill velocities near the nose (to 40
  degrees, blended over 20 degrees) and linear interpolation across the layer towards the shoulder; pressure and
  temperature linear across the layer; the field frozen at its release free stream; the shock smoothed over 3 mean
  free paths and the field faded over r0 behind the base plane; the composed drag law for the whole flight of a
  rarefied follower, with a 300 K follower wall; the pressure-gradient force on by default; the parent as an obstacle
  (`hit_parent`) for rarefied followers only; `kn_shock_layer_flag` 1e-3.
- AFFIx details (§4.3b): s_w default 2 degrees, `affix_mach_min` 5, rho_body = the follower's density, n_hat in the
  local vertical plane pointing up, and one shared far-wake integration per case across theta.
- The engine/adapter/scenario split (§4.0), the all-intervals output, and the validity record that flags
  but never refuses (so the droplet-in-sphere-plume use runs even where the plume law is poorly applicable).
- Use of `trajectory.Simulator` rather than the droplet drag law, because the spheres are 1 to 400 mm and
  hypersonic at release, the regime that integrator was verified in against SESAM. Exception: rarefied followers
  released into the shock layer, where the relative speed is low (§4.3d).
