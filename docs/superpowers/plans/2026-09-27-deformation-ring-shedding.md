# Deforming Body, Melt Ring and Fragment Shedding (Step 4) Plan

> **Status (2026-09-27): plan of record at the design stage — not yet executable.** It fixes the scope, the decisions
> taken so far, the measured facts they rest on, the milestones and the task breakdown of Step 4. Under the
> repository's prototype-then-plan rule, the verbatim code blocks and the measured facts of each task are added from a
> tested prototype once spec sections 6–12 are approved; until then, do not implement from this document. Step 4 is
> the work called "Step 3b" in the design sessions of 2026-09-27.

**Goal.** Let the semi-solid AA7075 body deform; collect the melt that the film carries rearward into a ring or a rear
cap (in flight) or a single tail (in the plasma wind tunnel, under gravity); shed fragments from it by **capillary
pinch-off**, **tearing of semi-solid material** and **forcing by the wake**; and record every fragment's properties
at the moment it separates. The **fate** of the fragments — trajectory, further heating, demise — is a later
iteration, because it depends on when and how fast each fragment leaves; the fragment record carries everything that
iteration will need.

**Builds on.** Step 3 as specified in its sub-plans (`docs/superpowers/plans/melt-spraying-subplans/`), including:
- the surface-recession amendment of 2026-09-27: the derived surface (`recession.py`, Task 16), periodic remeshing
  with conservative field transfer (`remesh.py`, Task 17) and `PHI_DEATH` = 0.50 (shared-context facts 38–45);
- the Scheil material variant `AA7075_scheil` (sub-plan 02, amendment of 2026-09-27, commit e6722ee). Step 4 uses it
  for the heat, the liquid fraction and the viscosity; `AA7075_range` remains Step 3's checked linear material.

**Spec.** To be written as `docs/superpowers/specs/<date>-deformation-ring-shedding-design.md` from the approved
sections recorded below and sections 6–12, which are still to be presented.

**Architecture.** One tetrahedral mesh for the whole body — the temperature mesh itself — moved by a slow viscous
(Stokes) flow whose viscosity is set by the temperature through the Scheil liquid fraction. The melt film stays a mass
per surface triangle, fed at the liquidus, so runny liquid never rides the mesh. Film that the loads carry to the rear
collects in a sub-grid **ring reservoir**; ring material that becomes stiff (its viscosity above the film-to-mesh
handover floor) **accretes** into the mesh, which therefore grows. Necks that thin and cells that tear are cut, and
every disconnected piece becomes a row of the **fragment table**. New modules (names provisional): `lee_loads` (or an
extension of `surface_flow`), `rim`, `accretion`, `rheology`, `stokes`, `surgery`, `shedding`, `fragments`, and
`deforming_body` holding `DeformingBody`, an extension of Step 3's `MeltingBody`; remeshing reuses Step 3's
`remesh.py`. A tunnel configuration (fixed hemisphere, gravity, tunnel conditions) reuses all of it.

**Per 0.5 s macro step** (approved Section 3, extended): Step 3's melt stage — feed, film transport (now over the
whole attached region), Girin spraying → rim collection, rim heat balance and sub-grid shedding → accretion of stiff
rim material into the mesh → deformation sub-steps (Stokes with surface tension and the body force; each sub-step
moves no node more than a quarter of its cell and strains no cell more than 10 %, at most 200 sub-steps) → neck and
tear checks → surgery and fragment rows → remesh when the quality envelope requires it → history row. Girin spraying
stays in the melt stage before the deformation, so each step sees one consistent shape (a one-step lag).

## Global constraints

- Interpreter, test tiers, SI units, repository conventions and the commit trailer exactly as in Step 3's
  `00-shared-context.md` (`Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`).
- DRAMA's databases and the wrapper scripts are never modified; anything non-physical is labelled as a verification
  device in the docstring, the CLI help and the README.
- Fixed attitude (tumbling excluded, labelled). No oxide skin (labelled: the tunnel test ran with very little oxygen,
  so it cannot test the skin, and flight air does contain oxygen). Self-contact of the deforming surface excluded.
- One Scheil material, `AA7075_scheil`, for heat, liquid fraction and viscosity.
- Every transfer — film to ring, ring to mesh, mesh or ring to fragments, old mesh to new mesh — books mass and
  enthalpy exactly; the balances close to round-off as in Step 3.
- Do not edit the package while a measurement is running (Step 3 fact 31).

## Decisions of record

- **Section 1, scope and placement** (approved 2026-09-27): its own step, built on the Step 3 prototype, with its own
  spec, prototype and plan. Film accretion and self-contact were deferred there; **film accretion at the rim was
  brought back into scope on 2026-09-27**, and self-contact stays deferred.
- **Section 2, components** (approved): rheology, stokes, remesh, surgery, shape loads and `DeformingBody`; one mesh
  for the whole body; the film a mass per triangle fed at the liquidus; loads per triangle from its own normal, plus
  a shelter (shadow) test.
- **Section 3, per-step sequence** (approved): as above; the heat capacity tied to each cell's fixed mass rather than
  its changing volume; the skfem backend first.
- **Section 4, flow law** (in progress): below 50 % solid, Li et al. (2014, China Foundry 11, 79–84),
  η = [0.871 − 0.00849 γ̇^0.749] exp(3.731 f_s) Pa·s, i.e. 0.87–6.8 Pa·s, which is film-side material below any
  workable handover floor; a jump at f_s = 0.5 (895.1 K under Scheil); the law above 50 % solid is still to come from
  Asha (a Sellars–Tegart form was proposed; its parameters are not yet known); the film-to-mesh handover floor is
  bracketed at 100 and 1000 Pa·s.
- **Physics scope** (2026-09-27): in flight the collection is axisymmetric (a ring or a rear cap), because gravity is
  not felt in the body frame and the deceleration acts along the axis; in the tunnel, gravity makes one tail. The three
  shedding mechanisms are capillary pinch-off, tearing of semi-solid material and forcing by the wake. Every fragment is
  recorded; its fate is deferred.
- **Milestone order** (approved 2026-09-27): (1) accretion with the loads behind the equator; (2) the deforming mesh
  with tearing; (3) the liquid rim with capillary pinch-off, wake forcing, the tunnel mode and the runs.
- **Material** (2026-09-27): one Scheil representation for the whole body, as its own material variant;
  `AA7075_range` unchanged.
- **Section 5, frame, loads and surface tension** (approved 2026-09-27):
  - Flight stays in the body frame at fixed attitude; the drag deceleration is a body force toward the nose, in the
    film (as in Step 3) and per unit volume in the Stokes flow. The tunnel mode fixes the body, applies one Earth
    gravity in a chosen direction and sets the deceleration to zero.
  - The attached-flow region is decided by each patch's geometric angle from the nose, out to a separation angle
    bracketed at 150° and 180° (no separation); the film moves over the whole attached region. A patch's own normal is
    used only for its local inclination and to shelter facets that face away from the local flow.
  - Behind the equator: the pressure continues Step 3's Prandtl–Meyer expansion down to the base pressure (bracketed
    at 1, 3 and 5 % of the stagnation pressure) and then stays there; the shear continues the boundary-layer
    construction from the local edge state, becoming a weak reverse shear (bracketed) behind a separation line if there
    is one; the heat flux is bracketed at 2, 5 and 10 % of the stagnation value (low confidence; continuum values are
    upper bounds at these Knudsen numbers). Anything protruding is loaded with the local flow at its own height.
  - The wake is a mean load in flight; its fluctuating part stays only in the tunnel mode, scaled by the forced
    response at the actual frequency ratio.
  - Surface tension (0.86 N/m) enters the Stokes solve in a curvature-free weak form; it limits each sub-step to about
    viscosity × cell size / surface tension (about 0.1 s at 100 Pa·s and 1 mm cells).

## Measured facts carried in (2026-09-27)

Measured with the unchanged Step 3 prototype and its linear material; to be re-measured with `AA7075_scheil`.

1. **The ring is fed by about one percent of the melt.** While the equator is intact (50 mm sphere 174.5–194.5 s;
   100 mm sphere 25.5–80 s), 1.0 g of 65.9 g of melt (1.5 %) and 4.4 g of 539 g (0.8 %) reach the last windward row;
   the rest sprays close to where it melts, as Step 3 found. Accumulated, it would form a half-round rim of radius
   1.3 mm and 2.0 mm.
2. **At the equator the film is still being pushed rearward.** It is 0.02–0.31 mm thick (50 mm) and at most 0.09 mm
   (100 mm), against a hold-back thickness 2τ/|G| — the thickness at which the deceleration's forward pull stops a
   sheared film — of 0.7–5.3 mm and 0.5–7.8 mm. The collection therefore forms behind the equator.
3. **The equator is cold.** It receives 0.3–1 % of the stagnation-point heat flux (6–18 kW/m² against 1.9–2.6 MW/m²)
   and sits at 747–890 K (50 mm) and 508–833 K (100 mm, below the solidus until about 65 s), so film arriving there
   becomes mush or solid: the ring starts as an accretion, not a liquid rim.
4. **Step 3 cannot hold that film.** Freeze-back is capped at an element's original mass, so the cold film stays
   booked as liquid and is sprayed at the equator: 0.94 of 1.0 g and 4.3 of 4.4 g.
5. **The late phase is a different regime.** The eroding front reaches the equator at about 195 s (50 mm) and 80 s
   (100 mm); the 100 mm remnant (0.49 kg, a third of the body) then stays at about 820–905 K down to 7.7 km (400 s)
   while decelerating at 1–10 g — semi-solid for about five minutes, and possibly the largest deformation in the model.
6. **The wake's fluctuations do not matter in flight.** Shedding at 0.1–0.2 × speed / diameter is 7–32 kHz, against
   24–760 Hz natural frequencies for liquid pieces 1–10 mm across: at most 1.1 % of the static deflection.
7. **Stripping at the equator is marginal.** The rim Weber number against the full edge flow is 1.6–3.6 (50 mm) and
   4.5–13.6 (100 mm).
8. **Scheil moves Step 3's headline modestly** (50 mm flight, same code): melt onset 0.74 km higher, demise 1 s later,
   sprayed mass +0.8 %, droplets −14 % with a median radius +4.4 % (sub-plan 02's amendment).

**From the literature** (primary sources as cited):
- Separation lies 130–149° from the nose at Mach 2 and Reynolds 300–1000 (Nagata et al. 2020, J. Fluid Mech. 904,
  A36); Hinman's fit θ_sep = 1117.1 Re_µ^−½ + 118.18° gives 126–150° within its range (Hinman 2017, PhD thesis,
  University of Calgary; Hinman & Johansen 2017, AIAA J. 55, 500) and extrapolates to none at the flight's Re_µ of
  20–300; DSMC found none near continuum (Dogra et al. 1994, J. Spacecr. Rockets 31, 713).
- The near wake is steady below an oscillation threshold that the flight sits far under (Fay & Goldburg 1963;
  Goldburg et al. 1965, refitted by Gai et al. 2024).
- Base pressure is about 3 ± 1 % of the pitot pressure for hypersonic laminar cylinder wakes (compiled by Hinman 2017);
  lee-side heating is 2–10 % of stagnation by analogy (Drouet et al. 2021, Acta Astronaut. 181, 446; Moss & Price
  1996, NASA TM-111599), with low confidence for a sphere.
- AA7050 semi-solid (castings): brittle between 90 and 97 % solid, grain coalescence at 94–97 %, strength 0.23–4.92
  MPa, a critical strain of order 0.1 % (Subroto 2014, PhD thesis, TU Delft; Subroto et al. 2021, Metall. Mater.
  Trans. A 52, 871); remelted extruded 7xxx fails brittle with about 2 % liquid (Chen et al. 2014, Mater. Des. 54,
  1); the hot-tearing criterion of Rappaz, Drezet & Gremaud (1999, Metall. Mater. Trans. A 30, 449). The semi-solid
  strength exceeds the aerodynamic loads a hundred- to thousandfold, so tearing needs strain imposed by flow or by
  thermal contraction.
- Wrought 7075 first melts on heating at 769–818 K (Brehm et al. 2022, SAND2022-9908; Gu et al. 2023, Materials 16,
  6145).
- A decelerating rim settles at Bond number one, b = √(σ/ρa): 6.4 mm at 0.9 g, 1.2 mm at 27 g, 6.0 mm at 1 g (Wang et
  al. 2018, Phys. Rev. Lett. 120, 204503); its drops are about 1.5 ligament widths (Wang & Bourouiba 2018, J. Fluid
  Mech. 848, 946); the fastest rim wavelength stays near 9 rim radii (Roisman 2010, J. Fluid Mech. 661, 206);
  wall-attached drops are blown off at a Weber number of 3.45–7.9 measured at their height (White & Schmucker 2008,
  J. Fluids Eng. 130, 061302; 2021, Phys. Rev. Fluids 6, 023601). Prefilming-atomizer correlations were measured at
  gas densities 10⁴–10⁵ times the flight's and transfer only qualitatively.

## File structure (provisional)

```
reentry_model/lee_loads.py (or surface_flow.py)  attached region by geometric angle, theta_sep, lee pressure/shear/heating, shelter
reentry_model/film.py                             Runoff over the attached region instead of windward-windward edges
reentry_model/frame.py                            body force per unit mass: flight (deceleration toward the nose) or tunnel (gravity)
reentry_model/rim.py                              ring reservoir: sectors, collection, along-ring flow, heat balance, Scheil state
reentry_model/accretion.py                        stiff rim material into new cells (derived-surface offset + remesh, exact transfer)
reentry_model/rheology.py                         eta(T, shear rate) from the Scheil solid fraction, handover floor, rigid core
reentry_model/stokes.py                           creeping flow on the tet mesh: body force, surface tension, surface tractions
reentry_model/surgery.py                          neck and tear cutting, connected components, exact accounting of removed pieces
reentry_model/shedding.py                         capillary pinch-off (sub-grid), tearing criterion, wake stripping
reentry_model/fragments.py                        the fragment table and its writers
reentry_model/deforming_body.py                   DeformingBody(MeltingBody): the extended per-step sequence
reentry_model/coupled.py, cli.py, viz.py          history columns, --tunnel and the Step 4 flags, fragment files, overlays
tests/test_reentry_model_{lee_loads,rim,accretion,rheology,stokes,surgery,shedding,fragments,deforming}.py
```

## Tasks

Each task will get, from the prototype: its verbatim code, its measured facts, the exact test commands and a commit.
Here each names its purpose, files, interfaces, tests and acceptance criterion.

### Phase B — the remaining spec sections

- [ ] **B1.** Present sections 6–12 one at a time for approval: 6 the rim and accretion (collection threshold,
  sectors, along-ring flow, heat balance, handover); 7 the shedding criteria (capillary, tearing, wake); 8 the fragment
  table; 9 remeshing and exact transfer (built on Step 3's Task 17); 10 the CFD surface output; 11 verification and
  validation; 12 prototype order and risks.
- [ ] **B2.** Write the spec, self-review it, commit it; Asha reviews it.

### Phase C, milestone 1 — accretion with the loads behind the equator

- [ ] **Task 1: Attached region and lee loads.** Purpose: the film and the loads continue past the equator. Files:
  `lee_loads.py` (or `surface_flow.py`), `film.py`. Produces: an attached mask from the geometric angle of each patch
  on the derived surface (Step 3 Task 16) with `theta_sep` (default 180°, bracket 150°); lee pressure (Prandtl–Meyer
  down to `f_base` × pitot, then constant); lee shear (boundary-layer continuation; reverse shear behind separation);
  lee heat flux (`f_lee` × stagnation); sheltered facets. Tests: the windward loads are unchanged; the pressure
  plateaus at the base value; film mass is conserved across 90°; a pit facet no longer blocks the runoff; the lee heat
  flux integrates to the bracket. Acceptance: both flights run with film reaching the lee side and the balances exact.
- [ ] **Task 2: Frame and body forces.** Files: `frame.py`, `cli.py` (`--tunnel`). Produces: the body force per unit
  mass in the body frame — flight, the deceleration toward the nose; tunnel, gravity in a chosen direction with no
  deceleration and a constant freestream. Tests: direction and magnitude in both modes; the tunnel mode holds the body
  still.
- [ ] **Task 3: Ring reservoir.** Files: `rim.py`. Produces: sectors along the collection line, each with mass,
  enthalpy and cross-section; collection of the film above a thickness threshold of the order of the capillary length
  under the local body force; flow along the ring driven by the tangential body force; a heat balance (gas heating,
  conduction into the body under it, radiation); the Scheil state. Tests: exact mass and enthalpy transfer from film to
  ring; axisymmetric forcing keeps the sectors equal; tunnel gravity drains the ring to its lowest sector.
- [ ] **Task 4: Accretion.** Files: `accretion.py` (with Step 3's `remesh.py`). Produces: ring material above the
  handover floor becomes new cells — the surface offset outward by the accreted volume, remeshed, mass and enthalpy
  booked exactly. Tests: exact transfer; element quality; the new volume equals the accreted mass over its density.
- **Milestone 1 acceptance:** with `AA7075_scheil`, both flights accrete the delivered film into a band or a rear cap
  at the place the lee loads select; the supply agrees with the re-measured probe; the balances close.

### Phase C, milestone 2 — the deforming mesh with tearing

- [ ] **Task 5: Rheology.** Files: `rheology.py`. Produces: η(T, γ̇) from the Scheil solid fraction — Li et al. below
  50 % solid (clamped where its bracket would turn negative, at about 484 1/s), the law above 50 % solid once
  supplied (a labelled placeholder bracket until then), rigid below the solidus, the handover floor. Tests: the jump at
  f_s = 0.5, monotonicity in f_s, the clamp.
- [ ] **Task 6: Stokes solve.** Files: `stokes.py`. Produces: a creeping-flow solve on the body mesh (skfem, MINI or
  P2/P1 elements) with the body force, surface tension in the curvature-free weak form, the aerodynamic tractions from
  Task 1, and the rigid core pinned. Tests against analytic solutions: the relaxation of a slightly deformed viscous
  drop; Papageorgiou's thinning of a viscous thread (the neck radius falls at 0.0709 σ/μ); the sag of a viscous layer
  under a body force; a rigid translation producing no strain.
- [ ] **Task 7: Mesh motion and remeshing.** Files: `deforming_body.py`. Produces: node motion with the Stokes
  velocity in sub-steps (quarter-cell move, 10 % strain, the capillary step limit, at most 200), a quality check, and
  Step 3's remesh with conservative transfer when the envelope is exceeded; the heat capacity on the fixed cell mass.
  Tests: conservation after the correction; the quality envelope; an undeformed translation.
- [ ] **Task 8: Tearing.** Files: `shedding.py`. Produces: the accumulated tensile strain of each cell in the brittle
  window (90–97 % solid, which Scheil puts between 833 K and the eutectic at the solidus), a critical strain
  (placeholder 0.1 %, bracketed), and optionally the Rappaz–Drezet–Gremaud criterion where ring material refreezes.
  Tests: a mushy bar in tension tears at the prescribed strain and at the right place; nothing tears outside the
  window.
- [ ] **Task 9: Surgery.** Files: `surgery.py`. Produces: removal of torn cells and of necks thinner than one cell,
  connected components of what remains, the component holding the core kept as the body, and every other component
  sent to the fragment table, with mass and enthalpy booked exactly. Tests: components found on synthetic meshes; a
  fragment's mass and enthalpy equal its cells'; the body's balance still closes.
- **Milestone 2 acceptance:** the ring deforms and sheds resolved fragments by necking or tearing; the 100 mm
  remnant's deformation over its semi-solid descent is measured, bracketed by the viscosity law; the balances close.

### Phase C, milestone 3 — the liquid rim, capillary pinch-off, wake forcing, runs

- [ ] **Task 10: Capillary pinch-off below the grid.** Files: `shedding.py`. Produces: the rim thickness from the
  Bond-number rule under the local body force; bulges spaced about nine rim radii apart; ligaments stripped when the
  rim's Weber number at its own height exceeds the depinning value (3.5–8, bracketed); drops about 1.5 ligament widths;
  the rate limited by supply and by the capillary or visco-capillary time. Tests: the Rayleigh–Plateau growth rate and
  wavelength of a liquid thread; conservation; no shedding below the threshold.
- [ ] **Task 11: Wake forcing.** Files: `shedding.py`, `lee_loads.py`. Produces: the mean load on protrusions from the
  local flow at their height; in the tunnel mode only, a fluctuating term scaled by the forced-oscillator response at
  the actual frequency ratio (Strouhal number bracketed 0.1–0.2). Tests: the analytic forced response; the flight
  response stays at about one percent or less.
- [ ] **Task 12: The fragment table.** Files: `fragments.py`, `coupled.py`, `cli.py`. Produces: the table below, its
  writers, and history columns (ring mass, fragments shed per mechanism, the largest fragment). Tests: the column list;
  body + film + ring + droplets + fragments = the initial mass and enthalpy, to round-off.
- [ ] **Task 13: Tunnel mode and validation.** Files: `mesh.py` (a hemisphere with a flat back and a sharp edge),
  `frame.py`, `cli.py`. Produces: a fixed hemisphere under gravity across the flow, the tunnel's measured heating and
  its loads (flow regime from the test record) with separation pinned at the edge. Validation against the video `IMG_2888.MOV`: one tail at the
  lowest point, fragments shed from its tip and edge, and fragment sizes relative to the body's height (median about
  2.5 %, 90 % below about 5 %, the largest bright piece about 13 % and one crumpled piece about a third; upper bounds
  because of glare). Needs the test record: sample diameter, heat flux, pressure, speed, gas, duration, frame rate.
- [ ] **Task 14: Flight runs, sensitivity and documentation.** The 50 and 100 mm flights with `AA7075_scheil`, the
  brackets (separation angle, base pressure, lee heating, the viscosity law above 50 % solid, the handover floor, the
  critical strain, the depinning Weber number), the measured facts, and the README, `docs/model_assumptions.md` and
  spec amendments.

### Phase D — the code-level plan

- [ ] Generate the implementation plan from the prototype (verbatim code, measured facts, the three mechanical
  audits), Asha reviews it, then subagent-driven execution in an isolated worktree, merged back to main locally.

## The fragment record

One row per resolved fragment, or per sub-grid drop population per step, parallel to Step 3's droplet source table:
- **When:** the separation time to within a sub-step, the macro step, and the flight state at that moment — altitude,
  flight speed, flight-path angle, freestream density and temperature, deceleration.
- **Where and how fast:** the position on the body (x, y, z and the angles θ and φ) and the velocity relative to the
  body at release.
- **Size and shape:** the mass and the count (1 for a resolved fragment), the volume, the equivalent spherical
  diameter, three principal lengths and the surface area.
- **Thermal state:** the mass-weighted specific enthalpy, the mean, lowest and highest temperature, the mass-weighted
  liquid fraction, and the heat still needed to melt it fully.
- **Origin:** the mechanism (capillary pinch-off, necking, tearing, wake stripping) and whether it was resolved on the
  mesh or handled below the grid.
- **Breakup numbers:** the Weber number on the source table's convention ρ∞V²(2r)/σ, the Ohnesorge number with the
  fragment's own viscosity, and the breakup flag.

## Verification and validation

Analytic tests before any physics run: the Rayleigh–Plateau thread, Papageorgiou thinning, viscous-drop relaxation,
the sagging layer, the mushy bar in tension and the forced oscillator. System checks: mass and energy closed through
film, ring, mesh and fragments; an axisymmetric ring at fixed attitude before it breaks up; one tail at the bottom in
the tunnel mode. Validation: the tunnel video, qualitatively until the test record gives the sample size and the frame
rate. Comparisons: the re-measured probe supply and the Step 3 results under `AA7075_scheil`.

## Open inputs and decisions

- The viscosity law above 50 % solid, with its parameters (Asha).
- The tunnel test record (Asha).
- Spec sections 6–12 (Phase B).
- Whether Step 3's `remesh.py` serves accretion and the deforming mesh as it is, or needs an extension.
- The re-measurement of Step 3's material-dependent facts under `AA7075_scheil` (Step 3's Task 14).

## Risks

- The viscosity law above 50 % solid is open, and the ring's stiffness, neck thinning, tearing and the remnant's
  deformation all hang on it.
- The lee loads are uncertain by a factor of a few until the CFD solver replaces them.
- Hot-tearing limits come from castings, which strain far more slowly than a re-entering body.
- Surgery must stay robust when a ring breaks into many pieces at once.
- The whole-remnant deformation of the 100 mm sphere may dominate the run time and the results.
- Time scales span from milliseconds (liquid pinch-off, below the grid) to minutes (stiff mush): the sub-stepping and
  the handover floor must keep each process on the right side of the grid.

## Out of scope

The fate of the fragments (trajectory, heating, demise) — the next iteration, fed by the fragment record; tumbling;
the oxide skin; self-contact of the deforming surface.
