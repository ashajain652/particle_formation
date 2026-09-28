# Deforming Body, Melt Ring and Fragment Shedding (Step 4) Plan

> **Status (2026-09-28): plan of record at the design stage — not yet executable.** It fixes the scope, the decisions
> taken so far, the measured facts they rest on, the milestones and the task breakdown of Step 4. Under the
> repository's prototype-then-plan rule, the verbatim code blocks and the measured facts of each task are added from a
> tested prototype once spec sections 7–13 are approved; until then, do not implement from this document. Step 4 is
> the work called "Step 3b" in the design sessions of 2026-09-27 and 2026-09-28. Two design sessions contributed; from
> 2026-09-28 this plan is maintained by one of them, and the revision of that date brings in its decisions on the flow
> laws, the film, surface tension and the material values (Sections 1–6 below).

**Goal.** Let the semi-solid AA7075 body change shape; collect the melt that the film carries rearward into a ring or
a rear cap (in flight) or a single tail (in the plasma wind tunnel, under gravity); shed fragments from it by
**capillary pinch-off**, **tearing of semi-solid material** and **forcing by the wake**; and record every fragment's
properties at the moment it separates. The **fate** of the fragments — trajectory, further heating, demise — is a
later iteration, because it depends on when and how fast each fragment leaves; the fragment record carries everything
that iteration will need.

**Builds on.** Step 3 as specified in its sub-plans (`docs/superpowers/plans/melt-spraying-subplans/`), including:
- the surface-recession amendment of 2026-09-27: the derived surface (`recession.py`, Task 16), periodic remeshing
  with conservative field transfer (`remesh.py`, Task 17) and `PHI_DEATH` = 0.50 (shared-context facts 38–45);
- the Scheil material variant `AA7075_scheil` (sub-plan 02, amendment of 2026-09-27, commit e6722ee). Step 4 uses the
  Scheil representation for the heat, the liquid fraction and the viscosity, with two values changed by the decisions
  of 2026-09-28 (latent heat 390 kJ/kg, liquid surface tension 0.80 N/m; see the open decisions for where they live);
  `AA7075_range` remains Step 3's checked linear material.

**Spec.** To be written as `docs/superpowers/specs/<date>-deformation-ring-shedding-design.md` from the approved
sections recorded below and sections 7–13, which are still to be presented.

**Architecture.** One tetrahedral mesh for the whole body — the temperature mesh itself. **In the first version the
mesh is held rigid:** Chen et al.'s (2016) law makes the coherent mush above 50 % solid immobile at the 1–10 kPa
re-entry loads except within about 1–2 K of the handover to the film (Section 4), and a strain check reports every
step whether that holds. The body changes shape in two ways: **recession**, as material leaves the cells for the film
(Step 3's derived surface and remeshing), and **accretion**, as film that freezes where its own cell is already full
grows the surface outward. The melt film stays a mass per surface triangle, handed over from the cells at 50 % solid
(895.1 K on the Scheil curve), and is carried as a **layered film**: a straight-line temperature profile through its
depth, with a liquid skin over a slurry base (Section 6). Film that the loads carry to the rear collects and freezes
into a ring or a rear cap (Section 7, to be presented). Necks that thin and cells that tear are cut, and every
disconnected piece becomes a row of the **fragment table**. A slow viscous (Stokes) deformation of the mesh is built
only as a sensitivity case, and only if the strain check flags a real region. New modules (names provisional):
`lee_loads` (or an extension of `surface_flow`), `frame`, `rheology`, `strain_check`, `accretion`, `rim`, `surgery`,
`shedding`, `fragments`, and `deforming_body` holding `DeformingBody`, an extension of Step 3's `MeltingBody`; `film.py`
is extended; remeshing reuses Step 3's `remesh.py`; `stokes` is deferred. A tunnel configuration (fixed hemisphere,
gravity, tunnel conditions) reuses all of it.

**Per 0.5 s macro step** (Section 3, revised and approved): trajectory, heating on the current shape, conduction;
Step 3's melt stage — hand-over of material at 50 % solid, layered-film transport over the whole attached region with
capillary pressure where the film's top is fully liquid, Girin spraying from the liquid skin, freezing from the film's
base → accretion of frozen material the owner cell cannot absorb (the surface grows; the mesh is rebuilt first if one
step's growth would exceed a quarter of a cell) → refresh of the conduction matrices, surface triangles and runoff
connections → strain check (reports only) → neck and tear checks → surgery and fragment rows → remesh when the quality
envelope requires it → history row. There are no deformation sub-steps in the first version. Girin spraying stays in
the melt stage, so each step's physics sees one consistent shape and responds to shape changes one step late.

## Global constraints

- Interpreter, test tiers, SI units, repository conventions and the commit trailer exactly as in Step 3's
  `00-shared-context.md` (`Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`).
- DRAMA's databases and the wrapper scripts are never modified; anything non-physical is labelled as a verification
  or numerical device in the docstring, the CLI help and the README — in Step 4 at least the eutectic smoothing at the
  solidus, the viscosity bridge just below the liquidus, the film's precursor thickness at wet–dry edges and the
  creep-floor grain size.
- Fixed attitude (tumbling excluded, labelled). No oxide skin as a separate mechanism (labelled: the tunnel test ran
  with very little oxygen, so it cannot test the skin, and flight air does contain oxygen). Self-contact of the surface
  detected and reported, not resolved.
- One Scheil representation for heat, liquid fraction and viscosity.
- Every transfer — cells to film, film to ring, film or ring to new cells, mesh or ring to fragments, old mesh to new
  mesh — books mass and enthalpy exactly; the balances close to round-off as in Step 3.
- Do not edit the package while a measurement is running (Step 3 fact 31).

## Decisions of record

- **Section 1, scope and placement** (approved 2026-09-27, revised and approved 2026-09-28): its own step, built on
  the Step 3 prototype, with its own spec, prototype and plan. In the first version: the Scheil material; the handover
  to the film at 50 % solid; the rigid mesh with the strain check; the film over the whole attached surface; film
  accretion (brought back into scope on 2026-09-27); recession, remeshing, neck surgery and the fragment table; loads
  on any shape by local inclination with the shelter test; the CFD hook (one function supplying per-triangle loads, and
  a watertight surface file written at chosen intervals). Deferred and recorded: the whole-body Stokes flow (sensitivity
  case only); self-contact (detected, reported); the CFD coupling itself; the fate of fragments (the next iteration);
  tumbling; liquid drainage through the mush's skeleton and shear-history (thixotropic) effects; the three other melting
  curves (noted, not run: the linear `AA7075_range`, Chen et al.'s 2016 measured curve, and that curve's shape
  stretched onto 750–908 K); Rokni et al. (2011), still unread.
- **Section 2, components** (approved 2026-09-27, revised and approved 2026-09-28): `material` (the Scheil variant);
  `rheology` (Chen 2016 above 50 % solid, Li et al. below it, the creep floor; a strain rate for a given stress,
  temperature and liquid fraction); `strain_check`; `accretion` (frozen film turned into outward motion of the
  triangle's nodes, the added volume equal to the frozen volume, the adjacent cells' mass and heat content updated
  exactly); `shape_loads`/`lee_loads` (the existing load models per triangle by its own angle, the shelter test, the
  lee loads, the body force — the single entry point a CFD solution will later fill); `film` (the whole attached
  region); `remesh`, `surgery` and the fragment table; `DeformingBody`; `stokes` deferred to the sensitivity case. The
  trajectory still sees only mass, area and a drag factor.
- **Section 3, per-step sequence** (approved 2026-09-27, revised and approved 2026-09-28): as above; each cell's heat
  capacity tied to its fixed mass rather than its changing volume; accretion books the frozen film's heat content into
  the cells it joins, as freeze-back does; a failed rebuild stops the run with exit code 1; the skfem backend first.
- **Section 4, flow law** (approved 2026-09-27 and 2026-09-28):
  - The liquid fraction is Scheil's, partition coefficient 0.4, pure-aluminium melting point 933 K, liquidus 908 K,
    with the 3.6 % eutectic remainder spread over a few kelvin at the solidus (labelled).
  - Below 50 % solid: Li Yageng et al. (2014, China Foundry 11, 79–84), η = [0.871 − 0.00849 γ̇^0.74924]
    exp(3.7311 f_s) Pa·s — 0.87 to 5.6 Pa·s at low shear rates. Measured only at f_s = 0.1–0.5 and γ̇ = 61–490 1/s,
    on continuously sheared slurry; γ̇ held at no more than 367 1/s, the highest rate where the fit still matches the
    data; a log-linear bridge from the f_s = 0.1 value down to liquid aluminium's 1.3 mPa·s over 906.4–908 K
    (labelled; no data). This material is film-side (Section 6).
  - A jump at f_s = 0.5 (895.1 K), where Li et al.'s data end and particle bridging begins. The film-to-mesh handover
    is at this jump; the earlier handover-floor bracket (100 and 1000 Pa·s) is no longer needed.
  - Above 50 % solid: Chen G. et al. (2016, J. Alloys Compd. 674, 26–36, doi 10.1016/j.jallcom.2016.02.254),
    Eq. 20, peak stress of extruded T6 7075 at 350–600 °C and 10⁻³–1 1/s. Semi-solid branch (0 < f_l < 0.5):
    (1 − 2 f_l)^15.81342 ε̇ = 1.28926×10⁸ [sinh(0.03006 σ)]^6.521135 exp(−175131 / RT), σ in MPa. Plastic branch as
    the solid limit: ε̇ = 6.75872×10¹¹ [sinh(0.014677 σ)]^6.85385 exp(−186163 / RT). The zero of strength at
    f_l = 0.5 was chosen by the authors, not measured; their data reach f_l = 0.19.
  - Consequence: evaluated with the model's temperature and liquid fraction, the mush strains less than 1 % in 45 s
    except within 0.4, 1.0 and 1.9 K of the handover at 1, 10 and 40 kPa — a band 0.07–0.25 mm thick at the nose's
    8–15 K/mm, thinner than the outermost 0.25 mm cells and treated as part of the handover. **The mesh is rigid in the
    first version.**
  - Strain check: each step, Chen's law at an estimated local stress (surface shear plus the pressure and body-force
    push times depth) plus a linear creep floor for flow through the liquid films between grains (grain size a
    parameter, 10–30 µm; a rough estimate gives at most 0.3 % strain in 45 s at 10 kPa for 30 µm grains and about 8 %
    for 10 µm); the accumulated strain is reported, and any region above 1 % is flagged. It never stops the run.
- **Physics scope** (2026-09-27): in flight the collection is axisymmetric (a ring or a rear cap), because gravity is
  not felt in the body frame and the deceleration acts along the axis; in the tunnel, gravity makes one tail. The three
  shedding mechanisms are capillary pinch-off, tearing of semi-solid material and forcing by the wake. Every fragment is
  recorded; its fate is deferred.
- **Milestone order** (approved 2026-09-27, milestone 2 revised 2026-09-28): (1) accretion with the loads behind the
  equator; (2) the strain check, tearing and surgery on the rigid mesh, with the Stokes deformation only as a
  sensitivity case if the check flags a region; (3) the liquid rim with capillary pinch-off, wake forcing, the tunnel
  mode and the runs.
- **Material values** (2026-09-28):
  - Heat capacity above 850 K: constant 1131.6 J/(kg K) — the value at which DRAMA's table ends — for solid, mush and
    liquid alike, uncertainty about ±4 % (1085–1180 J/(kg K)). No measured 7075 value above 850 K could be read; the
    compiled liquid value is 1130 (Mills 2002, reprinted in ASM Handbook Vol. 15, Table 4), the mass-weighted element
    sum 1132, a JMatPro curve about 1100, and liquid aluminium 1177 (NIST-JANAF) against 1127 measured (Leitner et al.
    2017). Apparent heat capacities that include the latent heat (peak about 13,400 J/(kg K)) must not be loaded,
    because the Scheil fraction already adds the latent heat.
  - Latent heat: 390 kJ/kg (Asha, 2026-09-28; source: Modulus Metal's AA7075-T6 data sheet, which lists 384–393 kJ/kg
    and cites no primary source). Alternatives noted, not run: 400 kJ/kg (DRAMA and `AA7075_range`), 397 kJ/kg (pure
    aluminium, NIST-JANAF), 358 kJ/kg (Mills 2002 via ASM Vol. 15).
- **Section 5, frame, loads and surface tension** (approved 2026-09-27; surface tension amended 2026-09-28):
  - Flight stays in the body frame at fixed attitude; the drag deceleration is a body force toward the nose, in the
    film (as in Step 3) and, in the Stokes sensitivity case only, per unit volume. The tunnel mode fixes the body,
    applies one Earth gravity in a chosen direction and sets the deceleration to zero.
  - The attached-flow region is decided by each patch's geometric angle from the nose, out to a separation angle
    bracketed at 150° and 180° (no separation); the film moves over the whole attached region. A patch's own normal is
    used only for its local inclination and to shelter facets that face away from the local flow.
  - Behind the equator: the pressure continues Step 3's Prandtl–Meyer expansion down to the base pressure (bracketed
    at 1, 3 and 5 % of the stagnation pressure) and then stays there; the shear continues the boundary-layer
    construction from the local edge state, becoming a weak reverse shear (bracketed) behind a separation line if there
    is one; the heat flux is bracketed at 2, 5 and 10 % of the stagnation value (low confidence; continuum values are
    upper bounds at these Knudsen numbers). Anything protruding is loaded with the local flow at its own height.
  - One per-triangle load interface — pressure, shear vector, convective heat flux and the body-force vector — filled
    by the engineering models now and by a CFD solution later; the film computes the pressure gradient that drives it
    from the triangle pressures itself, so it works on any shape and on CFD output.
  - The wake is a mean load in flight; its fluctuating part stays only in the tunnel mode, scaled by the forced
    response at the actual frequency ratio.
  - Surface tension enters **the film's motion** (there is no Stokes solve in the first version): the pressure inside
    the film is the gas pressure minus the surface tension times the curvature of the film's top surface, which is
    the body's curvature plus the bending of the film thickness. It acts **only where the film's top is fully liquid**
    (above 908 K); a film with a slurry top gets no capillary term (Asha, 2026-09-28; the slurry case is recorded as an
    open assumption, because no measurement of the surface tension of a semi-solid aluminium alloy is known).
  - Value: **0.80 N/m** for fully liquid 7075 (Asha, 2026-09-28), source Bainbridge & Taylor (2013, Metall. Mater.
    Trans. A 44A, 3901–3909, doi 10.1007/s11661-013-1696-9), Table II, commercial 7075 by sessile drop at the liquidus
    plus 50 K: 0.809 ± 0.041 N/m as melted in vacuum, 0.843 ± 0.018 after the oxide skin was broken, 0.777 ± 0.061
    after exposure to dry air, 0.607 ± 0.083 after breaking the skin again in dry air. Constant with temperature. The
    same value is used for Girin spraying and the droplet and fragment Weber numbers in Step 4 (Step 3 keeps
    0.86 N/m); droplets come out up to about 7 % smaller and Weber numbers about 7.5 % higher.
  - Numerics: the capillary term is solved implicitly over the surface triangles each step (a 1 mm liquid rim on 2 mm
    triangles smooths out in about 10⁻⁴ s); a very thin precursor film is kept at wet–dry edges (labelled device); the
    remesher refines the surface where film collects, because a 1–3 mm rim otherwise spans a single 2 mm triangle.
- **Section 6, the film** (approved 2026-09-28):
  - Layered film: each triangle's film has a straight-line temperature profile from its base, at the surface node
    temperature, to a top set by the gas heating conducted through it. Heat crosses 1 mm in 0.02 s and 2 mm in 0.09 s
    (conductivity 128 W/(m K)), far shorter than the 0.5 s step, so the profile is treated as settled within a step;
    the straight line holds while the film is thinner than about 2–3 mm.
  - The viscosity follows the temperature through the depth (Li et al. with the bridge, liquid above the liquidus);
    the flow rate is the depth integral of (shear on top plus the pressure, body-force and capillary push times the
    height above the base) divided by the local viscosity.
  - Girin spraying acts on the liquid skin only, with Step 3's spray model and liquid properties unchanged; the slurry
    part is too viscous to spray (Ohnesorge number 0.6–3.3 for a millimetre slurry film, against 0.009 for a 10 µm
    liquid film). Freezing, and then accretion, happen from the base once it falls below 895.1 K.
  - Energy booking (Asha's rule): each surface node's account holds its own heat content at its temperature plus a
    third of the integrated film energy of every triangle it belongs to; the film energy is the heat content (the
    constant heat capacity above 850 K plus the Scheil latent term) integrated through the profile, not the film mass
    times the heat content at the node temperature. Every transfer — the handover, runoff, spraying, freezing — carries
    the integrated energy of the slice it removes, so the balance stays exact. The solver keeps the film's heat
    capacity in its implicit step (a 2 mm film holds 4.8 kg/m², seven times the outermost 0.25 mm cells), while the
    accounting uses the integrated profile. Radiation stays at the node temperature; a top 50 K hotter would radiate at
    most 0.36 % of the heating more (emissivity 0.4).
  - Rough expected scale (estimates, not model runs): a slurry film 0.7–2.6 mm thick where the viscosity is 1–5.6 Pa·s,
    with 11–40 K across it at 2 MW/m² — more than the 13 K between the handover and the liquidus, hence the liquid
    skin.
- **Section 6 of the first design session, rim and accretion** (presented 2026-09-27, not approved): per-patch stocks
  on the rear half — film, rim (thicker than the Bond thickness or the patch width, moving as a lump by the film's
  lubrication law, pinned below a depinning Weber number of 3.5–8) and crust (frozen, awaiting cells) — with crust
  turned into cells when the thickest crust exceeds 0.35 of a cell. To be re-presented as Section 7, reconciled with
  the layered film (Section 6) and the node-displacement accretion of Sections 1–3.

## Measured facts carried in (2026-09-27)

Measured with the unchanged Step 3 prototype and its linear material; to be re-measured with the Scheil material.

1. **The ring is fed by about one percent of the melt.** While the equator is intact (50 mm sphere 174.5–194.5 s;
   100 mm sphere 25.5–80 s), 1.0 g of 65.9 g of melt (1.5 %) and 4.4 g of 539 g (0.8 %) reach the last windward row;
   the rest sprays close to where it melts, as Step 3 found. Accumulated, it would form a half-round rim of radius
   1.3 mm and 2.0 mm.
2. **At the equator the film is still being pushed rearward.** It is 0.02–0.31 mm thick (50 mm) and at most 0.09 mm
   (100 mm), against a hold-back thickness 2τ/|G| — the thickness at which the deceleration's forward pull stops a
   sheared film — of 0.7–5.3 mm and 0.5–7.8 mm. The collection therefore forms behind the equator. (A Newtonian
   estimate of a stopping line near 81° that ignored the shear on thin films was withdrawn on 2026-09-28.)
3. **The equator is cold.** It receives 0.3–1 % of the stagnation-point heat flux (6–18 kW/m² against 1.9–2.6 MW/m²)
   and sits at 747–890 K (50 mm) and 508–833 K (100 mm, below the solidus until about 65 s), so film arriving there
   becomes mush or solid: the ring starts as an accretion, not a liquid rim.
4. **Step 3 cannot hold that film.** Freeze-back is capped at an element's original mass, so the cold film stays
   booked as liquid and is sprayed at the equator: 0.94 of 1.0 g and 4.3 of 4.4 g.
5. **The late phase is a different regime.** The eroding front reaches the equator at about 195 s (50 mm) and 80 s
   (100 mm); the 100 mm remnant (0.49 kg, a third of the body) then stays at about 820–905 K down to 7.7 km (400 s)
   while decelerating at 1–10 g — semi-solid for about five minutes. Under Chen 2016 it stays rigid at these loads
   unless the strain check says otherwise.
6. **The wake's fluctuations do not matter in flight.** Shedding at 0.1–0.2 × speed / diameter is 7–32 kHz, against
   24–760 Hz natural frequencies for liquid pieces 1–10 mm across: at most 1.1 % of the static deflection.
7. **Stripping at the equator is marginal.** The rim Weber number against the full edge flow is 1.6–3.6 (50 mm) and
   4.5–13.6 (100 mm).
8. **Scheil moves Step 3's headline modestly** (50 mm flight, same code): melt onset 0.74 km higher, demise 1 s later,
   sprayed mass +0.8 %, droplets −14 % with a median radius +4.4 % (sub-plan 02's amendment).

**Evaluated on 2026-09-28** (from the literature and the model's material curves; not model runs):

9. **The coherent mush is stiff at re-entry loads.** Near-solidus compression of extruded 7075 (Gu et al. 2023,
   Materials 16, 6145) converts to 1.1×10⁶–8.7×10⁸ Pa·s; nanoparticle-reinforced 7075 held 4 min and globular still
   needs 2.6–6 MPa at 0.01 1/s at 85–97 % solid (Chen & Yan 2017, J. Alloys Compd. 708, 751); an AI search summary's
   10⁴–8×10⁴ Pa·s near the solidus is not supported and matches only held, globular 7075 at 45–60 % solid tested at
   about 1 1/s and 0.06–0.2 MPa (Zoqui & Torres 2010, Materials Research 13, 305).
10. **Li et al.'s law is the runny side.** At its lowest measured rate (61 1/s) it underestimates the data by 1.4–1.9
    times at f_s = 0.3–0.5, because the exponent is averaged over rates (the 61 1/s fit alone has 4.78, not 3.7311).
11. **Where the latent heat goes.** Half of it is absorbed by 895.1 K on the Scheil curve, 897.8 K on Chen et al.'s
    (2016) measured curve (their Fig. 1, digitised) and 829.0 K on the linear curve; by 850 K the three have absorbed
    54, 32 and 253 kJ/kg. Scheil's peak absorption, 26.5 kJ/(kg K) at 908 K, is a sharp edge the Newton solver must pass.
12. **Latent heat dominates the heat content between 850 and 908 K** (84 % of the heat absorbed there); the 390 versus
    400 kJ/kg choice moves it by about 9 kJ/kg, the heat-capacity uncertainty by about 2.6 kJ/kg.

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
  6145); Chen et al. (2016) put the solidus of their extruded T6 7075 at 752 K with under 1 % liquid below 813 K.
- A decelerating rim settles at Bond number one, b = √(σ/ρa): 6.4 mm at 0.9 g, 1.2 mm at 27 g, 6.0 mm at 1 g (Wang et
  al. 2018, Phys. Rev. Lett. 120, 204503); its drops are about 1.5 ligament widths (Wang & Bourouiba 2018, J. Fluid
  Mech. 848, 946); the fastest rim wavelength stays near 9 rim radii (Roisman 2010, J. Fluid Mech. 661, 206);
  wall-attached drops are blown off at a Weber number of 3.45–7.9 measured at their height (White & Schmucker 2008,
  J. Fluids Eng. 130, 061302; 2021, Phys. Rev. Fluids 6, 023601). Prefilming-atomizer correlations were measured at
  gas densities 10⁴–10⁵ times the flight's and transfer only qualitatively.
- Liquid 7075 surface tension by sessile drop under argon at 923–1073 K rises with temperature from 719–820 to
  917–943 mN/m depending on the substrate, which the authors attribute to oxygen at the surface (Uba & Raush 2025,
  J. Manuf. Mater. Process. 9, 165).

Supporting material outside the repository (Asha's machine): the literature reviews
`/Users/ashajain/reports/Semisolid AA7075 constitutive parameters.md` and `.../AA7075 heat capacity above 850 K.md`,
the digitised curve `.../Chen 2016 Fig 1 liquid fraction.csv` and the figure `.../AA7075 latent heat curves.png`.

## File structure (provisional)

```
reentry_model/lee_loads.py (or surface_flow.py)  attached region by geometric angle, theta_sep, lee pressure/shear/heating, shelter
reentry_model/frame.py                            body force per unit mass: flight (deceleration toward the nose) or tunnel (gravity)
reentry_model/film.py                             layered film over the attached region: profile, depth-integrated flux, capillary term
reentry_model/rheology.py                         Li et al. below 50 % solid (bridge, clamp), Chen 2016 above, creep floor; strain rate for a stress
reentry_model/strain_check.py                     per-cell stress estimate and accumulated strain; flags above 1 %
reentry_model/accretion.py                        frozen film beyond cell capacity into outward surface growth (+ remesh, exact transfer)
reentry_model/rim.py                              ring or cap collection (Section 7, to be reconciled with the layered film)
reentry_model/surgery.py                          neck and tear cutting, connected components, exact accounting of removed pieces
reentry_model/shedding.py                         capillary pinch-off (sub-grid), tearing criterion, wake stripping
reentry_model/fragments.py                        the fragment table and its writers
reentry_model/deforming_body.py                   DeformingBody(MeltingBody): the extended per-step sequence
reentry_model/stokes.py                           deferred: creeping flow for the sensitivity case only
reentry_model/coupled.py, cli.py, viz.py          history columns, --tunnel and the Step 4 flags, fragment files, overlays
tests/test_reentry_model_{lee_loads,film_layered,rheology,strain_check,accretion,rim,surgery,shedding,fragments,deforming}.py
```

## Tasks

Each task will get, from the prototype: its verbatim code, its measured facts, the exact test commands and a commit.
Here each names its purpose, files, interfaces, tests and acceptance criterion. Task numbers are provisional and will
be renumbered in the code-level plan.

### Phase B — the remaining spec sections

- [ ] **B1.** Present sections 7–13 one at a time for approval: 7 the rim and accretion (re-presented, reconciled
  with the layered film and the node-displacement accretion); 8 the shedding criteria (capillary, tearing, wake); 9
  the fragment table; 10 remeshing and exact transfer (built on Step 3's Task 17); 11 the CFD surface output; 12
  verification and validation; 13 prototype order and risks.
- [ ] **B2.** Write the spec, self-review it, commit it; Asha reviews it.

### Phase C, milestone 1 — accretion with the loads behind the equator

- [ ] **Task 1: Attached region and lee loads.** Purpose: the film and the loads continue past the equator. Files:
  `lee_loads.py` (or `surface_flow.py`), `film.py`. Produces: an attached mask from the geometric angle of each patch
  on the derived surface (Step 3 Task 16) with `theta_sep` (default 180°, bracket 150°); lee pressure (Prandtl–Meyer
  down to `f_base` × pitot, then constant); lee shear (boundary-layer continuation; reverse shear behind separation);
  lee heat flux (`f_lee` × stagnation); sheltered facets; the per-triangle load interface. Tests: the windward loads are
  unchanged; the pressure plateaus at the base value; film mass is conserved across 90°; a pit facet no longer blocks
  the runoff; the lee heat flux integrates to the bracket. Acceptance: both flights run with film reaching the lee side
  and the balances exact.
- [ ] **Task 2: Frame and body forces.** Files: `frame.py`, `cli.py` (`--tunnel`). Produces: the body force per unit
  mass in the body frame — flight, the deceleration toward the nose; tunnel, gravity in a chosen direction with no
  deceleration and a constant freestream. Tests: direction and magnitude in both modes; the tunnel mode holds the body
  still.
- [ ] **Task 2a: Layered film, capillary pressure and the energy booking.** Files: `film.py`, `material.py` (the
  Scheil variant with Step 4's values), `deforming_body.py`. Produces: the per-triangle temperature profile, the
  depth-integrated flux with Li et al.'s viscosity and the bridge, the capillary term where the film's top is fully
  liquid (implicit, with the precursor film), spraying from the liquid skin, freezing from the base, and the node
  energy accounts with the integrated film energy. Tests: a film of uniform viscosity reproduces Step 3's lubrication
  flux; the depth integral against quadrature; a sinusoidal thickness perturbation of a liquid film decays at the
  analytic capillary rate; the integrated film energy against quadrature of the profile; no capillary flux under a
  slurry top; the balances close.
- [ ] **Task 3: Ring or cap collection.** Files: `rim.py`. Produces: collection of the film where it thickens and
  freezes behind the equator, per Section 7 once approved (previously: sectors with mass, enthalpy and cross-section,
  collection above a capillary-length threshold, along-ring flow, a heat balance). Tests: exact mass and enthalpy
  transfer from film to ring; axisymmetric forcing keeps the collection axisymmetric; tunnel gravity drains it to the
  lowest point.
- [ ] **Task 4: Accretion.** Files: `accretion.py` (with Step 3's `remesh.py`). Produces: frozen material that its
  owner cell cannot absorb becomes outward surface growth — the triangle's nodes moved so the added volume equals the
  frozen volume, remeshed when the growth or the quality requires it, mass and enthalpy booked exactly. Tests: exact
  transfer; element quality; the new volume equals the accreted mass over its density.
- **Milestone 1 acceptance:** with the Scheil material, both flights accrete the delivered film into a band or a rear
  cap at the place the lee loads select; the supply agrees with the re-measured probe; the balances close.

### Phase C, milestone 2 — strain check, tearing and surgery on the rigid mesh

- [ ] **Task 5: Rheology.** Files: `rheology.py`. Produces: Li et al. below 50 % solid (γ̇ held at no more than
  367 1/s; the bridge to the liquid over 906.4–908 K); the jump at f_s = 0.5; Chen 2016's semi-solid branch above and
  its plastic branch as the solid limit; the linear creep floor; the strain rate for a given stress, temperature and
  liquid fraction. Tests: Chen's peak stresses at his own test conditions (8.2 MPa at 600 °C and 10⁻³ 1/s, 22.2 MPa
  at 1 1/s); the jump; monotonicity in f_s; the clamp and the bridge.
- [ ] **Task 5a: Strain check.** Files: `strain_check.py`, `deforming_body.py`. Produces: the per-cell stress
  estimate, the accumulated strain from Chen's law and the creep floor, the flagged volume above 1 %, and history
  columns. Tests: a cell at a prescribed stress and temperature accumulates the analytic strain; no strain below the
  solidus beyond the creep floor.
- [ ] **Task 6: Stokes solve (sensitivity case only; deferred).** Built only if the strain check flags a real region.
  Files: `stokes.py`. Produces: a creeping-flow solve on the body mesh (skfem, MINI or P2/P1 elements) with the body
  force, surface tension in the curvature-free weak form, the aerodynamic tractions from Task 1, and the rigid core
  pinned. Tests against analytic solutions: the relaxation of a slightly deformed viscous drop; Papageorgiou's thinning
  of a viscous thread (the neck radius falls at 0.0709 σ/μ); the sag of a viscous layer under a body force; a rigid
  translation producing no strain.
- [ ] **Task 7: Mesh motion (sensitivity case only; deferred).** Files: `deforming_body.py`. Produces: node motion
  with the Stokes velocity in sub-steps (quarter-cell move, 10 % strain, the capillary step limit, at most 200), a
  quality check, and Step 3's remesh with conservative transfer; the heat capacity on the fixed cell mass. Tests:
  conservation after the correction; the quality envelope; an undeformed translation.
- [ ] **Task 8: Tearing.** Files: `shedding.py`. Produces: the accumulated tensile strain of each cell in the brittle
  window (90–97 % solid, which Scheil puts between 833 K and the eutectic at the solidus) from imposed strain — thermal
  contraction on refreezing, and flow only in the sensitivity case — a critical strain (placeholder 0.1 %, bracketed),
  and optionally the Rappaz–Drezet–Gremaud criterion where ring material refreezes. Tests: a mushy bar in tension tears
  at the prescribed strain and at the right place; nothing tears outside the window.
- [ ] **Task 9: Surgery.** Files: `surgery.py`. Produces: removal of torn cells and of necks thinner than one cell
  (necks thin by recession at the root of a rim or tail on the rigid mesh), connected components of what remains, the
  component holding the core kept as the body, and every other component sent to the fragment table, with mass and
  enthalpy booked exactly. Tests: components found on synthetic meshes; a fragment's mass and enthalpy equal its
  cells'; the body's balance still closes.
- **Milestone 2 acceptance:** the strain check runs on both flights and reports where the accumulated strain exceeds
  1 %; necking by recession and tearing shed resolved fragments; if the check flags a real region, Tasks 6–7 are built
  and the sensitivity case is run; the balances close.

### Phase C, milestone 3 — the liquid rim, capillary pinch-off, wake forcing, runs

- [ ] **Task 10: Capillary pinch-off below the grid.** Files: `shedding.py`. Produces: the rim thickness from the
  Bond-number rule under the local body force; bulges spaced about nine rim radii apart; ligaments stripped when the
  rim's Weber number at its own height exceeds the depinning value (3.5–8, bracketed); drops about 1.5 ligament widths;
  the rate limited by supply and by the capillary or visco-capillary time; surface tension 0.80 N/m where liquid.
  Tests: the Rayleigh–Plateau growth rate and wavelength of a liquid thread; conservation; no shedding below the
  threshold.
- [ ] **Task 11: Wake forcing.** Files: `shedding.py`, `lee_loads.py`. Produces: the mean load on protrusions from the
  local flow at their height; in the tunnel mode only, a fluctuating term scaled by the forced-oscillator response at
  the actual frequency ratio (Strouhal number bracketed 0.1–0.2). Tests: the analytic forced response; the flight
  response stays at about one percent or less.
- [ ] **Task 12: The fragment table.** Files: `fragments.py`, `coupled.py`, `cli.py`. Produces: the table below, its
  writers, and history columns (ring mass, fragments shed per mechanism, the largest fragment). Tests: the column list;
  body + film + ring + droplets + fragments = the initial mass and enthalpy, to round-off.
- [ ] **Task 13: Tunnel mode and validation.** Files: `mesh.py` (a hemisphere with a flat back and a sharp edge),
  `frame.py`, `cli.py`. Produces: a fixed hemisphere under gravity across the flow, the tunnel's measured heating and
  its loads (flow regime from the test record) with separation pinned at the edge. Validation against the video
  `IMG_2888.MOV`: one tail at the lowest point, fragments shed from its tip and edge, and fragment sizes relative to the
  body's height (median about 2.5 %, 90 % below about 5 %, the largest bright piece about 13 % and one crumpled piece
  about a third; upper bounds because of glare). Needs the test record: sample diameter, heat flux, pressure, speed,
  gas, duration, frame rate.
- [ ] **Task 14: Flight runs, sensitivity and documentation.** The 50 and 100 mm flights with the Scheil material, the
  brackets (separation angle, base pressure, lee heating, the creep-floor grain size, the critical strain, the
  depinning Weber number), the measured facts, and the README, `docs/model_assumptions.md` and spec amendments.

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
- **Breakup numbers:** the Weber number on the source table's convention ρ∞V²(2r)/σ with σ = 0.80 N/m, the Ohnesorge
  number with the fragment's own viscosity, and the breakup flag.

## Verification and validation

Analytic tests before any physics run: the layered film's flux against quadrature and against Step 3's uniform film,
the capillary decay of a sinusoidal film perturbation, Chen's peak stresses at his own conditions, the strain check on
a prescribed cell, the Rayleigh–Plateau thread, the mushy bar in tension and the forced oscillator (and, for the
sensitivity case only, Papageorgiou thinning, viscous-drop relaxation and the sagging layer). System checks: mass and
energy closed through cells, film, ring, accretion and fragments; an axisymmetric ring at fixed attitude before it
breaks up; one tail at the bottom in the tunnel mode. Validation: the tunnel video, qualitatively until the test record
gives the sample size and the frame rate. Comparisons: the re-measured probe supply and the Step 3 results under the
Scheil material.

## Open inputs and decisions

- **Where Step 4's material values live.** The committed `AA7075_scheil` (Step 3 sub-plan 02, e6722ee) carries
  400 kJ/kg and 0.86 N/m; Step 4 needs 390 kJ/kg and 0.80 N/m. Either sub-plan 02 is amended (which changes Step 3's
  Scheil variant as well) or Step 4 gets its own variant (Asha).
- The surface tension of a film with a slurry top: none applied, recorded as an open assumption.
- ESA's measured 7075 heat capacity and heat of fusion through melting (Pagan 2025, University of Stuttgart
  dissertation, doi 10.18419/opus-17557; Bonvoisin et al. 2022, CEAS Space Journal 15, 213–235,
  doi 10.1007/s12567-022-00429-0): not yet read; would supersede the latent heat and heat capacity above.
- The tunnel test record (Asha); the tunnel gas had very little oxygen, so the tunnel mode may need the vacuum surface
  tension (0.81–0.84 N/m in Bainbridge & Taylor's Table II) rather than 0.80 N/m.
- Spec sections 7–13 (Phase B).
- Whether Step 3's `remesh.py` serves accretion as it is, or needs an extension.
- The re-measurement of Step 3's material-dependent facts under the Scheil material (Step 3's Task 14).

## Risks

- Chen 2016 is extrapolated from megapascal tests to kilopascal loads, where slow linear creep through the liquid films
  between grains — grain-size dependent and unmeasured — could make the mush above 50 % solid move; the strain check
  and the creep floor measure that risk rather than remove it.
- The lee loads are uncertain by a factor of a few until the CFD solver replaces them.
- The layered film's straight-line profile fails above about 2–3 mm, which a thick rim can exceed; surface tension
  shapes rims only where the triangles are well under the rim's size.
- Hot-tearing limits come from castings, which strain far more slowly than a re-entering body.
- Surgery must stay robust when a ring breaks into many pieces at once.
- Scheil's sharp latent-heat edge at the liquidus (26.5 kJ/(kg K)) may cost Newton iterations; the run records them.
- Time scales span from milliseconds (liquid pinch-off, below the grid) to minutes (stiff mush): the handover and the
  sub-grid treatments must keep each process on the right side of the grid.

## Out of scope

The fate of the fragments (trajectory, heating, demise) — the next iteration, fed by the fragment record; tumbling;
the oxide skin as a separate mechanism; self-contact of the surface; the whole-body Stokes flow in the first version
(sensitivity case only).
