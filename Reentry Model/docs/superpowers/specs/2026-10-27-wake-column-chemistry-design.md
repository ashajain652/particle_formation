# Wake column: aluminium vapour chemistry and condensation behind a demising droplet — design

Date: 2026-10-27. Status: draft for Asha's review. Package to be built: `reentry_model/wake_column/`.

## 1. Intent

A molten AA7075 droplet shed into the wake of a re-entering satellite loses up to about 45 % of its mass
by vaporisation (`reentry_model/droplet_reentry`, 315-case sweep, 1–200 µm, 50–90 km). This model follows
that aluminium vapour after it leaves the droplet: how it cools, how it is oxidised by the ambient air,
and how it condenses into alumina (Al2O3) or aluminium-hydroxide nanoparticles. It answers one research
question:

> Of the aluminium a droplet vaporises, how much ends up in nanoparticles, of what composition and size,
> fed by which gas-phase molecules, and over what time — and how does that change with release altitude
> and between day and night?

The literature behind it is in `Particle Wake Evolution/Literature Review/Alumina Nanoparticles/`
(31 PDFs and `reports/Alumina nanoparticles from aluminum vapor.md`). Its central findings shape this
design: gaseous Al2O3 barely exists, so particles are born by chemistry (two aluminium-bearing molecules
binding), at close to the collision rate, not by classical nucleation of a supersaturated vapour; AlO is
long-lived in cold air because its oxidation by O2 is either uphill or needs a third body; the final
amount and size are set by how fast the plume dilutes compared with how fast aluminium molecules meet.

## 2. Goals, non-goals, anti-goals

Goals (first version):
1. Time histories, over 12 hours, of every gas-phase aluminium species in the plume (Al, AlO, OAlO,
   OAlO2, AlCO3, AlOH, Al(OH)2 — the one-aluminium species of the network in §4.5) and of the condensed
   aluminium. Any species with two or more aluminium atoms is a cluster in the first size bin (§4.6).
2. The condensed fraction of the vaporised aluminium against time, and the times to 50 % and 90 %.
3. A particle size distribution at chosen times, split into alumina-like and hydroxide-like material.
4. Provenance: for the condensed aluminium, which gas molecule delivered it and by which route
   (new-particle formation or growth on an existing particle).
5. The mass that settles out of the plume.
6. All of the above for each release altitude, for a day case and a night case.

Non-goals (first version; written up in §11 as future work):
- Horizontal advection by wind (the planned second phase).
- Ion chemistry (Plane et al. 2021 reactions R18–R35).
- Tier-two (RRKM/master-equation) and tier-three (molecular-dynamics) estimates of reaction and sticking
  efficiencies.
- Aggregate (fractal) particle shape; hydroxide dehydration to alumina; zinc and magnesium chemistry.
- A computational-fluid-dynamics simulation of the wake (rejected: see §4.1).
- Overlapping trails of many droplets: Asha defined the case as one droplet on its independent descent,
  so each trail is modelled in isolation.

Anti-goals — the model is a failure, even if it runs, if:
- an unmeasured rate or efficiency is hidden in the code instead of being a named, documented parameter
  that a sensitivity run can change;
- the composition or provenance of the particles cannot be read back from the output;
- aluminium is not conserved to integration tolerance (gas + particles + settled = vaporised);
- the condensation pathways are collapsed back into a single lumped sink.

## 3. Constraints and conventions

- Python, in the Reentry Model project, run with `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`
  from inside `Reentry Model/`. Tests under `reentry_model/wake_column/tests/`, collected by `pytest.ini`
  (add the path, as for `droplet_reentry`). Analysis scripts under `reentry_model/wake_column/analysis/`,
  run with `-m`, writing to `output/wake_column_output/` (git-ignored).
- Run names encode the full configuration (diameter, release altitude, day/night, ambient source,
  efficiency settings) so sweeps share one output folder and stay resumable (repository convention).
- Every new physical assumption gets an entry in `docs/model_assumptions.md` with its verification status.
- Inputs that the tests depend on (the extracted WACCM columns, the reaction table) are committed under
  `data/wake_column/` or `tests/fixtures/`; outputs are not.

## 4. Decided physics (as Asha decided it, unless marked "my call")

### 4.1 Why a box model, and what it must add

The droplets are in free-molecule or transitional flow during vaporisation (Knudsen number about 0.4 at
50 km to over 100 at 90 km for 200 µm; Reynolds number mostly below 10), so there is no continuum wake and
no vortex shedding; a fluid-dynamics simulation is not needed. Two things a fixed box would miss are
included instead (Asha, 2026-10-27):
- **Diffusive widening.** Molecular diffusion in thin air is fast (diffusion coefficient roughly
  0.01 m²/s at 50 km to 3 m²/s at 90 km). The plume widens by diffusion from the first version.
- **A hot start.** The vapour leaves with the droplet's speed and is slowed by a few collisions with air;
  a separate slowing-down calculation sets the cylinder's starting composition and radius (§4.4).

### 4.2 Geometry: a horizontal cylinder along the vapour trail (Asha accepted, 2026-10-27)

The brainstorm assumed the vapour forms over "a few metres" and spoke of a vertical column. Running the
droplet code (2026-10-27) shows the trail is mostly horizontal and longer than that: the droplet leaves the
parent at 4–7.6 km/s on a flight-path angle of about −1°, and 90 % of its aluminium loss occurs over

| case | path length to 90 % of Al loss | time | altitude drop over it |
|---|---|---|---|
| 1 µm, 80 km | 51 m | 9.9 ms | 0.7 m |
| 10 µm, 80 km | 513 m | 99 ms | 7.3 m |
| 100 µm, 70 km | 1,191 m | 258 ms | 33 m |
| 10 µm, 60 km | 40 m | 16.5 ms | 2.6 m |

(The trail length scales with the droplet diameter because free-molecule deceleration scales with
droplet mass per area.) So the plume is a near-horizontal line source, like a meteor trail. The box is a
cylinder whose axis is the trail:
- **Length** L = the path length containing 90 % of the droplet's aluminium loss (from the droplet history
  rows: integrate `V_ms` over `t_s`; aluminium lost from `m_Al_kg`).
- **Inventory** = all the aluminium the droplet vaporises (so conservation closes against the droplet
  code's total), spread uniformly along L and deposited at t = 0 (deposition takes 10–260 ms, short
  against every chemical time scale except Al + O2; see §11 for a segmented-trail refinement). Every
  fraction in the outputs is a fraction of this total.
- **Altitude** = the mean altitude of the deposited aluminium (mass-weighted).
- **Radius** r(t)² = r0² + 4 ∫ D_total dt, with D_total = D_mol + K_eddy · S(r). D_mol is the binary
  diffusion coefficient of AlO in air from Chapman–Enskog kinetic theory using the ambient n and T.
  K_eddy is the eddy diffusion coefficient profile of Megner et al. (2006, ACP 6, 4415, Fig. 3c,
  "reference profile"; in the papers folder), digitised into a committed CSV (Asha, 2026-10-27). S(r) is a
  smooth switch from 0 to 1 as r passes the turbulence inner scale `l_inner` (my call: default 30 m,
  sensitivity 10 and 100 m), because eddies smaller than the plume do not exist below that scale.
  Megner's own sensitivity range (× 0.1 and × 100) is added to the sensitivity set. A molecular-only run
  (K_eddy = 0) is kept as a comparison case. r0 is the stopping distance from the slowing-down
  calculation (§4.4): millimetres at 60 km to tens of centimetres at 90 km.
- **Axial spreading.** The trail also spreads along its axis, and for short trails (1 µm droplets, about
  50 m) the radius exceeds the length within tens of minutes at 80–90 km. The box volume is therefore the
  effective volume of a line segment of length L broadened by a Gaussian of variance 2 D t in all three
  directions, which tends to cylindrical growth (V ∝ t) while r is much less than L and to spherical growth
  (V ∝ t^1.5) once r exceeds L. The builder derives the expression and tests both limits.
- **Entrainment.** As the cylinder widens, ambient air enters: every gas species i obeys
  dc_i/dt = (chemistry) + (c_i,ambient − c_i) · (1/V) dV/dt, with V = π r² L. Particles share the box
  volume and are diluted at the same rate (my call: particle Brownian diffusion is slower than gas
  diffusion; this is listed as an assumption).
- **Settling.** Particles in a well-mixed horizontal cylinder leave through its lower half at the rate
  2 v_s / (π r), where v_s is the free-molecule (Epstein) settling speed of the bin. Settled mass is counted
  and removed. (Asha asked for gravitational settling out of the column; with a horizontal cylinder it
  leaves through the side, not an end.)

### 4.3 Ambient air (Asha: WACCM baseline; NRLMSIS + AFGL cross-check)

- **Baseline:** the full-chemistry WACCM starting-condition file
  `Space Sustainability/WACCM-Cheyenne/carma-cesm2-inputdata/FWmaCARMAnudged.f19_f19_mg17.carma_trop_strat.2000.001.cam.i.2001-01-01-00000_c220325.nc`
  (64-bit-offset netCDF, readable with `scipy.io.netcdf_file`; 70 levels to about 140 km). Columns at the
  latitude nearest 40.7° N: **day** = longitude 180° (local noon at 00 UT), **night** = longitude 0°
  (local midnight). Variables used: T, O, O2, O3, H, H2, Q (water vapour), CO2, CO, plus PS, hyam/hybm/hyai/hybi/P0
  for pressure, mass mixing ratios converted to number density, geometric altitude from hypsometric
  integration. N2 = remainder of the air. Extract once to a committed CSV in `data/wake_column/`
  (the 1.6 GB file is outside the repository).
- **Cross-check:** NRLMSIS 2.x via `pymsis` (installed in drama_env) for T, n, N2, O2, O, H at the same
  latitude, day of year 1 and local times 12:00 and 00:00; AFGL 1986 midlatitude-winter for O3, H2O, CO2
  (tables from Anderson et al. 1986; the `joseki` package ships them as CSV). The cross-check is a second
  ambient source selectable by flag, run for the same cases.
- **Held fixed for 12 hours** (Asha): a day run holds the noon column, a night run the midnight column.
  AlOH photolysis (J = 3.3 × 10⁻³ s⁻¹, Plane R17) is on in day runs and off in night runs.
- The plume's own ambient (T, n, species) comes from WACCM, not from the droplet code's US76 atmosphere;
  the release altitude is matched by geometric altitude (my call).

### 4.4 The hot start: a slowing-down calculation (Asha, 2026-10-27; details my call — see §10)

Asha's picture was vapour forming at about the droplet temperature (908–1,913 K in the droplet sweep) and
cooling to ambient. The droplet code shows that each atom also leaves with the droplet's speed relative to
the air, 1.5–7.6 km/s, so its first collisions are hyperthermal. In the centre-of-mass frame of the two
colliding molecules, the collision energy is E = ½ μ V², with μ the reduced mass of the pair:

| speed | Al + air (μ = 14.0 amu) | AlO + O2 (μ = 18.3 amu) | equivalent temperature μV²/(3k_B), Al + air |
|---|---|---|---|
| 7.6 km/s | 403 kJ/mol | 530 kJ/mol | about 32,000 K |
| 3.6 km/s | 91 kJ/mol | 119 kJ/mol | about 7,300 K |
| 1.5 km/s | 16 kJ/mol | 21 kJ/mol | about 1,300 K |

Above about 3.5 km/s, AlO + O2 collisions exceed the 106 kJ/mol needed for AlO + O2 → OAlO + O, a channel
that is closed in ambient air. Because aluminium and air molecules have nearly equal masses, a hard-sphere
collision removes on average 2mM/(m+M)² ≈ 0.50 of the atom's energy, so the atom is thermalised in about
5–7 collisions: about 1 µs and 2–3 mm of travel at 60 km, 5 µs and 1 cm at 70 km, 25 µs and 5 cm at 80 km,
110 µs and 20 cm at 90 km (hard-sphere estimates in US76 air). One gas temperature cannot represent this:
it would put the rate expressions at 7,000–32,000 K, far outside their fitted ranges (about 2,000–4,000 K
at most), and clipping them there would remove the threshold reaction that matters. Asha therefore chose
(2026-10-27) to compute the slowing-down separately and to hand its products and stopping distance to the
cylinder:

**Stage A — slowing-down (a test-particle Monte Carlo, run once per trail before the box model).**
- *Emission.* The adapter (§5) supplies the droplet history over the trail: time, altitude, droplet speed
  vector V, droplet temperature T_d and aluminium loss rate. Emission times are sampled in proportion to
  the aluminium loss rate (default 20 time samples, `n_emission_samples`, with `n_test_particles`,
  default 10⁴, atoms in total). Each atom leaves the droplet surface with a half-Maxwellian
  (cosine-law) velocity at T_d in the droplet frame, at a uniformly sampled point on the sphere, plus the
  droplet velocity V.
- *Collisions.* The atom flies freely for an exponentially distributed time set by its total collision
  rate with the ambient species (N2 as the remainder, O2, O and H2O from the same WACCM day or night column as the box, at
  the emission altitude, with Maxwellian partner velocities at the ambient temperature). Collisions are
  hard-sphere, with collision diameters from the same Lennard-Jones parameters as the Chapman–Enskog
  diffusion coefficient (§4.2), and isotropic scattering in the centre-of-mass frame. Inelastic energy loss
  to the partner's vibration is neglected (this makes the stopping distance slightly long; assumption
  entry).
- *Reactions at each collision.* For each reaction open to the current species and partner, the
  probability that a collision reacts is the line-of-centers form P(E) = P0 (1 − E0/E) for E above the
  threshold E0, and zero below it. E0 and P0 are calibrated to the thermal rate expressions already in the
  network, since the line-of-centers thermal rate is exactly k(T) = P0 σ_HS v̄(T) exp(−E0/k_B T): E0 is the
  expression's activation temperature times k_B, and P0 = k(T_ref) exp(E0/k_B T_ref) / (σ_HS v̄(T_ref)) at
  the top of its validity range, capped at 1. Barrierless reactions (Al + O2 → AlO + O) get E0 = 0, that
  is a constant probability per collision. The reactions are: Al + O2 → AlO + O; AlO + O2 → OAlO + O;
  AlO + O → Al + O2; OAlO + O → AlO + O2 (Plane R3); and collisional dissociation AlO + M → Al + O + M with
  E0 = the AlO bond energy (about 5.2 eV) and P0 = 1 as an upper bound (open only near 7.6 km/s with heavy
  partners; my call). After a reaction the product keeps the reactant's velocity (momentum is nearly
  conserved because the departing O is light relative to the pair's total; exothermicity is neglected)
  and continues slowing with its own mass. The probabilities can be scaled together by
  `hot_reaction_scale` (default 1, capped at 1 per collision).
- *End of the stage.* An atom or molecule is thermalised when its kinetic energy in the air's frame falls
  below 2 × (3/2) k_B T_ambient (`thermal_cutoff`, default 2). Its species and its position relative to
  the droplet track at that moment are recorded.
- *Checks it logs.* The share of collisions above the aluminium ionisation energy (578 kJ/mol, not reached
  by Al + air at 7.6 km/s) and above the AlO bond energy, so that leaving ions out (§2) is confirmed per
  trail rather than assumed.

**What Stage A hands to the cylinder.**
- *Initial composition:* the fractions of the inventory thermalised as Al, AlO and OAlO (plus any
  dissociated Al), weighted over the emission history. My rough estimate before building: about half or
  more of the aluminium is already AlO, and a few percent OAlO at the highest speeds; the OAlO share is the
  number that matters, since ambient air turns Al into AlO within milliseconds anyway.
- *Initial radius:* r0 = the root-mean-square distance of the thermalisation points from the droplet track,
  measured perpendicular to it (the displacement along the track only shifts the trail by centimetres
  against a length of tens of metres or more). This replaces the "3 mean free paths" rule.
- *Temperature:* from t = 0 the cylinder is at the ambient temperature of its column, constant for the
  12 hours. Every rate expression in the box is therefore evaluated within its validity range; if one is
  not, the expression at the nearest bound is used and a warning is logged.
- Heat released by the reactions and deposited by the slowing atoms is neglected (the aluminium is a trace
  in the air once r exceeds a few r0); this is an assumption entry.

**Time-scale separation (why the stages can be split).** Slowing down takes 1–110 µs, while laying the trail
down takes 10–260 ms, so each packet of vapour is thermalised long before the next is emitted, and long
before the cylinder's chemistry (whose fastest step, Al + O2 in ambient air, takes milliseconds) acts.

**Sensitivity case (Asha).** `hot_start = "droplet_temperature"` runs the same Stage A with the droplet
velocity set to zero, so that atoms leave at the droplet temperature only. This is Asha's original picture;
it has no threshold reactions and gives an r0 of a few mean free paths. The default is
`hot_start = "slowing_down"`.

### 4.5 Chemistry (gas phase)

Backbone: the neutral reactions R1–R17 of Plane et al. (2021, JGR Space Physics 126, e2020JA028792,
Table 1), exactly as tabulated (in the folder as `Plane_2021_Meteor_ablated_aluminum_MLT.pdf`), with their
bath gas N2 for three-body steps. Added (from Gobrecht et al. 2022 Table B.1 and Mangan
et al. 2020, both used by the research notes):
- AlO + O2 → AlO2 + O, k = 3.54 × 10⁻¹¹ (T/300)^0.19 exp(−13140/T); and its reverse,
  1.90 × 10⁻¹⁰ (T/300)^0.17;
- AlO + O → Al + O2, k = 7.15 × 10⁻¹² exp(−945/T);
- Al + O + M → AlO + M, AlO + O + M → AlO2 + M (from Gobrecht Table B.1);
- AlO + H2 → AlOH + H (exothermic by 45 kJ/mol; Gobrecht Table B.1), since WACCM supplies H2.
AlO + AlO → Al2O + O is not a gas reaction here: its product has two aluminium atoms, so it is a
new-particle pathway (§4.6). The gas network is the species of R1–R17 plus these additions; none of them
needs OH or HO2, which the WACCM file does not carry. In ambient air the endothermic additions are
negligible; they are kept so that the same table, with its thresholds, drives the slowing-down stage (§4.4). Each added reaction carries its source and validity
range in the reaction table (§8).

### 4.6 Condensation pathways and provenance (Asha: track how each molecule feeds each particle)

Particles are sectional (Asha wants sizes): bins of doubling volume from the two-aluminium cluster
(about 0.5 nm) to a few hundred nanometres (about 40 bins; builder's choice of exact count).

**New-particle formation** — a reaction between two aluminium-bearing gas molecules puts one cluster into
the first bin. Each pair is a separately labelled pathway:
- AlO + AlO (+M) → Al2O2 (the (AlO)x route of Gobrecht/Patzer/Lam), and AlO + AlO → Al2O + O;
- AlO + AlOH, AlOH + AlOH, AlO + Al(OH)2, AlOH + Al(OH)2, Al(OH)2 + Al(OH)2 (Plane's polymerisation set,
  R36, split into its pairs);
- AlO2 + AlOH → Al2O3 + H (exothermic by 72.5 kJ/mol);
- OAlO + AlO, OAlO2 + AlO, OAlO2 + AlOH and any other pair of aluminium-bearing species in the network
  (builder: enumerate all pairs; none is dropped silently).
Rate (tier one, Asha 2026-10-27): dipole–dipole capture rate × statistical spin factor × an efficiency
parameter (default 1), one parameter per pathway. Dipole moments and spin multiplicities are a committed
table with sources; where a dipole moment is unknown, the gas-kinetic collision rate is used and flagged.

**Growth** — any aluminium-bearing gas molecule colliding with a particle of any bin sticks (free-molecule
collision rate including the molecule's own size; sticking efficiency parameter per species, default 1)
and moves the particle's volume up.

**Oxygen uptake** — a particle with an oxygen-to-aluminium ratio below 1.5 takes up oxygen from O, O2 and
O3 collisions at the collision rate × an oxidation efficiency (default 1), until it reaches 1.5. Hydrogen
is gained only through the hydroxide molecules that condense (no gas-surface reaction with H or H2O in the
first version; my call).

**Uptake onto background meteoric smoke** (Asha, 2026-10-27; its own route in the provenance) — the
ambient air carries meteoric smoke particles with the number density and size distribution of the Megner
et al. (2006) reference profile (Fig. 2, digitised to a committed CSV by altitude and radius). Any
aluminium-bearing gas molecule hitting a smoke particle sticks at the free-molecule collision rate times
a sticking efficiency (default 1). The smoke is a fixed background inside the box (entrained with the
ambient air, not depleted, not coagulated with plume particles in the first version — my call); the
aluminium it collects is tallied by source molecule under the route "background smoke" and is reported
separately from plume-made particles, since it ends up in particles of meteoric composition.

**Coagulation** — free-molecule Brownian coagulation kernel between bins (the regime holds at every
altitude in the sweep), particles as compact spheres (my call; aggregates are future work).

**Reversibility** — condensation is irreversible (my call: Al–O clusters are bound by 400–500 kJ/mol per
atom and the plume is at ambient temperature within a fraction of a millisecond).

**Composition and provenance per bin.** Each bin carries: particle number; aluminium atoms delivered by
each source molecule, split by route (formation or growth); total oxygen and hydrogen atoms. Coagulation
adds the tallies of the colliding particles. Particle mass follows from the atoms; volume from a
composition-weighted density (alumina 3.0 g/cm³, aluminium hydroxide 2.42 g/cm³; my call, documented).
Alumina-like vs hydroxide-like is read from each bin's H/Al ratio (H/Al below 0.5 counts as alumina-like;
the threshold is a reporting choice, not physics, and the raw ratios are output).

### 4.7 Run matrix and end condition

- 12 hours of simulated time per run (Asha), output on a logarithmic time grid from 1 µs.
- "Steady state" is judged inside that window: the composition is reported as steady once every
  aluminium reservoir changes by less than 1 % per decade of time (reported, not used to stop).
- Default sweep: release altitudes 90 to 50 km in 2 km steps (21, as the droplet sweep), diameters 1, 14.1
  and 200 µm (the droplet sweep's subset), day and night: 126 runs, plus the same with the NRLMSIS + AFGL
  ambient. The full 15-diameter sweep is an option.
- Zinc and magnesium vaporised mass is carried and reported but undergoes no chemistry.
- Sensitivity set (a deliverable, because the tier-two and tier-three decisions depend on it; my call on
  its extent): for the 14.1 µm droplet at 60, 76 and 90 km (76 rather than 75 so the baseline is on the
  sweep's grid), day and night, one-at-a-time changes to (a) all formation efficiencies × 0.1 and × 0.01,
  (b) all sticking efficiencies × 0.1, (c) the oxidation efficiency × 0.1, (d) `hot_reaction_scale` 0.1
  and 3, (e) `hot_start = "droplet_temperature"` (Asha), (f) formation rates × 80 (Plane's
  polymerisation enhancement), (g) K_eddy × 0.1, × 100 (Megner's range) and 0 (molecular only),
  (h) `l_inner` 10 and 100 m. 6 cases × 13 settings = 78 runs, the baselines coming from the default
  sweep, with a tornado plot of the condensed fraction and the mean diameter at 12 hours.

## 5. Inputs from the droplet stage

`reentry_model.droplet_reentry.model.run_case(d0_um, h0_km)` returns `(final_state, rows)`; each row has
`t_s, h_km, lat_deg, lon_deg, V_ms, T_K, m_Al_kg, mdot_Al_kgs, m_Zn_kg, m_Mg_kg, ...`. The wake column takes
from the rows: L (§4.2), the aluminium inventory over L, the mass-weighted mean altitude, and for the
slowing-down stage (§4.4) the history over the loss of time, altitude, the droplet velocity vector (from the
speed and the track between successive rows), droplet temperature `T_K` and aluminium loss rate
`mdot_Al_kgs`. An adapter function does this;
the droplet code is not modified. Running the droplet case inside the wake run is acceptable (seconds per
case) and avoids a stale-cache problem.

## 6. Outputs

Per run, in `output/wake_column_output/<run name>/`: a time-history CSV (every gas species, r, T,
condensed Al by source and route, settled Al, ambient values); a size-distribution CSV at fixed times
(1 ms, 1 s, 1 min, 10 min, 1 h, 12 h) with per-bin composition and provenance; a JSON summary (condensed
fraction at the end, times to 50 % and 90 %, number and volume mean diameters, alumina-like share,
provenance shares, settled share, the steady-state time). Sweep-level plots: condensed fraction against
time for each altitude (day and night panels); provenance stack by altitude; size distributions by
altitude; day-minus-night comparison; WACCM against NRLMSIS + AFGL.

## 7. Verification (tests)

1. Aluminium conservation to integration tolerance in every test run.
2. Inert tracer: with chemistry off, concentration falls as r0²/r(t)² exactly.
3. Ambient stability: with no aluminium, the gas stays at the ambient column.
4. Single-reaction checks: AlO lifetime against R6 in ambient air at 70, 80, 90 km equals
   1/(k6 [N2][O2]) (about 0.08 s, 2 s, 40 s in the night WACCM columns — the builder recomputes for the
   committed columns); Al lifetime against R1.
5. Coagulation only, from a monodisperse start: the distribution approaches the free-molecule
   self-preserving form, geometric standard deviation 1.45 ± 0.03.
6. Settling only: a single bin decays at 2 v_s / (π r).
7. Provenance closure: the source-and-route tallies sum to the condensed aluminium in every bin and time.
8. The reaction table reproduces Plane's tabulated values at 200 K and 300 K.
9. Slowing-down, reactions off, partners at rest: the mean fractional energy loss per collision equals
   2mM/(m+M)² (0.50 for Al on N2) to Monte Carlo tolerance, and the stopping distance scales as 1/n.
10. Line-of-centers calibration: sampling Maxwellian collision energies at a temperature inside each
   reaction's validity range reproduces its tabulated k(T) to within 5 %.
11. Slowing-down with zero droplet speed and the droplet at the ambient temperature: no threshold reaction
   occurs, and r0 is of order one mean free path.
12. Packet conservation: every test particle ends as exactly one of Al, AlO or OAlO, and Stage A's shares
   sum to one.

## 8. The parameter and reaction tables

The reaction network is data, not code: a committed table (CSV or YAML) with, per reaction, reactants,
products, rate form and coefficients, validity range, source, and whether it is Plane R-numbered or added.
The efficiency parameters (one per formation pathway, sticking per species, oxidation, plus
`hot_start`, `hot_reaction_scale`, `thermal_cutoff`, `n_test_particles`, `n_emission_samples`, the random
seed, `l_inner`, the K_eddy scale) live in one settings object with
documented defaults and are part of the run name when changed.

## 9. Left to the builder

Numerical method (stiff ODE solver; scipy `solve_ivp` with BDF or LSODA is expected), exact bin count and
spacing, the data structure for tallies, the adapter's interpolation details, plot styling (reuse
`DRAMA Sphere Experiments/analysis/plot_style.py` as `compare.py` does), parallelism of the sweep.

## 10. Calls I made that Asha has not confirmed

Marked "my call" above: the length/inventory rule and axial spreading of the cylinder (§4.2); the
turbulence inner scale; smoke held fixed and not coagulated with plume particles;
species with two or more aluminium atoms counted as first-bin clusters; the sensitivity set (§4.7); particles diluted with the gas; settling rate
through the side; WACCM ambient rather than US76 for the plume; within the slowing-down stage (§4.4):
hard-sphere isotropic scattering with no inelastic loss, line-of-centers probabilities calibrated at the top
of each validity range, a constant probability for barrierless reactions, AlO dissociation with P0 = 1,
products keeping the reactant's velocity, the thermalisation cutoff of twice the thermal energy, and r0
as the root-mean-square perpendicular stopping distance; the box held at ambient temperature from t = 0;
76 km in the sensitivity set; neglect of reaction heat; no gas–surface H/H2O
reactions; compact spheres; irreversible condensation; particle densities; the H/Al reporting threshold;
the 126-run default sweep.

## 11. Future work (not built now)

- **Horizontal advection** (phase two): a wind carrying species out of the cylinder, as an extra loss
  term with the ambient wind speed and the cylinder's cross-wind size.
- **Tier two — rate theory, if tier one shows sensitivity to the smallest formation steps:** RRKM /
  master-equation (MESMER) falloff for AlO + AlO, AlO + AlOH and AlOH + AlOH at each altitude's pressure
  and temperature, using the energies and frequencies of Gobrecht et al. (2022) and Plane et al. (2021).
- **Tier three — sticking by molecular dynamics, if tier one shows sensitivity to sticking:** reactive or
  machine-learned force fields, molecules fired at alumina and hydroxide clusters by size and temperature.
- **Energy-resolved hot cross-sections, if the sensitivity set shows the slowing-down products matter:**
  quasi-classical trajectories on Al + O2 and AlO + O2 potential surfaces, replacing the line-of-centers
  form above the fitted temperature range.
- **Segmented trail:** split the trail into segments by deposition time, each with its own inventory,
  temperature and altitude, with axial diffusion between them.
- Aggregate particle shape (fractal dimension about 1.6–2) for settling and coagulation.
- Hydroxide dehydration to alumina.
- Ion chemistry (Plane R18–R35) and particle charging.
- Zinc and magnesium chemistry.
- A diurnal cycle instead of fixed day or night columns.
