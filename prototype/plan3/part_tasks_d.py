# ---------------------------------------------------------------------------------------------------------- Task 14
emit(r'''### Task 14: Verification and sensitivity drivers, the reference-tier test, the runs

**Files:**
- Create: `analysis/melt_verification.py`, `analysis/melt_sensitivity.py`, `tests/test_reentry_model_reference_melt.py`
- Outputs (git-ignored): `reentry_model_output/verification_melt/{summary.md, summary.json, <case>__<mode>.*}`, `.../sensitivity/{sensitivity.md, sensitivity.json, ...}`, `.../girin/` (Task 8), `.../reference_tests/`

**Interfaces:**
- Consumes: the CLI (Task 13), `compare.melt_metrics` (Task 12), the melting references (Task 11), `sesam_io.REFERENCE_DIR`.
- Produces: the verification table (bookkeeping thresholded, resolved and physics reported), the sensitivity table, the reference-tier test with the spec §13.1 thresholds (mass 2 % of m₀, onset 0.5 km, 1 %-mass time 2 %).

- [ ] **Step 1: Create `tests/test_reentry_model_reference_melt.py`**

''')
emit(file_block("tests/test_reentry_model_reference_melt.py"))
emit(r'''
- [ ] **Step 2: Create `analysis/melt_verification.py`**

''')
emit(file_block("analysis/melt_verification.py"))
emit(r'''
- [ ] **Step 3: Create `analysis/melt_sensitivity.py`**

''')
emit(file_block("analysis/melt_sensitivity.py"))
emit(r'''
- [ ] **Step 4: Run the reference-tier test and the verification driver**

```bash
"$PY" -m pytest -m reference tests/test_reentry_model_reference_melt.py -q          # ~1 min: both bookkeeping runs
"$PY" analysis/melt_verification.py --animate                                         # ~6 min: 2 spheres x 3 modes, with the videos
```

Expected: 2 passed (measured 0.98 % / +0.10 km / −1.3 % and 1.29 % / +0.17 km / −0.2 %); the driver prints the six-row table with `pass` in both bookkeeping rows and the resolved and physics rows of the README's verification table (Task 15), and writes `summary.md`. The whole reference tier (`"$PY" -m pytest -m reference -q`) is 15 passed in ≈ 24 min. Check that `reentry_model_output/verification_melt/d100__physics/vtk/` holds `animation.mp4`, `film.mp4`, `section.mp4` and the stills (`melt_onset`, `spraying_onset`, `peak_release`, `film_*`, `section_*`), and that `d100__resolved/mass_time.png` shows the SESAM overlay with its residual panel.

- [ ] **Step 5: Run the sensitivity study**

```bash
"$PY" analysis/melt_sensitivity.py                                                    # ~60 min: 13 variants x 2 spheres
"$PY" analysis/melt_sensitivity.py --variants fenicsx --python "$FX"                  # the backend variant from fenicsx_env (CC and FI_PROVIDER exported)
```

(export `CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang` before the second command; the script sets `FI_PROVIDER=tcp` itself.) Expected: `sensitivity.md` with every variant's sprayed mass, median radius, onsets and demise altitude and the change relative to `base`; the `fenicsx` row equal to `base` to the printed digits; `layers2`/`layers6`/`dt025` within a few percent of `base` in demise altitude and sprayed mass; `we308`, `kt±30` nearly identical (spraying is melt-limited), `kr-30` moving the median radius; `AA7075` earlier onset (850 K vs 908 K liquidus) and a leeward remnant that reaches the ground (demise n/a); `sizeinitial` ending higher (the nose stays at R₀, the body keeps the sphere's drag). `gammapm14` (`--gamma-pm 1.4`) and `knbody003` (`--kn-body-shock 0.003`) bound the two new modelling choices: the first halves the wall pressure beyond the sonic point, the second shrinks the Girin-certified band. Record the table in Task 15.

- [ ] **Step 6: Run the FEniCSx tests once more and the whole unit tier**

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
"$PY" -m pytest -m "not drama and not reference" -q
```

Expected: 8 passed; the unit tier green (~6 min).

- [ ] **Step 7: Commit**

```bash
git add analysis/melt_verification.py analysis/melt_sensitivity.py tests/test_reentry_model_reference_melt.py
git commit -m "Add the melt verification, sensitivity drivers and the reference-tier bookkeeping test (Step 3 Task 14)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 15
emit(r'''### Task 15: Documentation — README, assumptions, spec amendments, facts note

**Files:**
- Modify: `README.md` (a new Step 3 section before "## Tests"; the Step 2 verification table refreshed; the Tests section's test counts and the reference-tier line), `docs/model_assumptions.md` (§9 appended, §7 amended), `docs/superpowers/specs/2026-09-20-melt-spraying-design.md` (§18 "Amendments" appended; status line), `reentry_model/__init__.py` (docstring: Steps 1–3), the facts note `/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md` (§17 appended)

- [ ] **Step 1: Re-run the Step 2 verification with the new thermal core and refresh its README table**

```bash
"$PY" analysis/reentry_model_thermal_verification.py
```

Expected (measured 2026-09-21 with the nodal-enthalpy core): d100 sesam 0.46 % / 2.23 % / +0.31 % / +0.30 % / 24.6 K (1.20 %) / 6.66 % / 2089 K at 139 s / 74 s; d050 sesam 0.37 % / 2.36 % / +0.03 % / +0.02 % / 28.0 K (1.22 %) / 3.40 % / 2442 K at 245 s / 24 s; physics ratios 0.744 / 0.770, peaks 2319 K at 118 s and 2494 K at 239 s, 83 s / 33 s. Replace the numbers in the README's Step 2 verification table with the measured ones (the runtimes halve: `energy()` and the operators are cheaper and the Newton takes 2.0 iterations) and add one sentence after the table: "Re-measured on 2026-09-21 with Step 3's nodal-enthalpy core (lumped capacity matrix): every metric within 0.1 K / 0.01 % of the Step 2 values; runtimes halved."

- [ ] **Step 2: Add the Step 3 section to `README.md`**

Insert before `## Tests`. The verification table below carries the prototype's own measurements; replace each row with Task 14's `summary.md` where the re-run differs, and write the sensitivity findings from `reentry_model_output/verification_melt/sensitivity/sensitivity.md` into the paragraph that names the fifteen variants:

''')
emit(block(open(os.path.join(HERE, "readme_step3.md")).read(), "markdown"))
emit(r'''
Also update the Tests section: the unit-tier line becomes `# unit tests (~6 min; the thermal solver, coupled and melting tests dominate)` and the reference-tier line `# the Step 1 reference flights (~15 min), the Step 2 coupled runs (~20 min) and the Step 3 bookkeeping runs (~1 min)`; and in the `fenicsx_env` paragraph "the seven conformance tests" → "the eight conformance tests (Step 2's seven and the Step 3 melting cross-check)".

- [ ] **Step 3: Append §9 to `docs/model_assumptions.md` and amend §7**

In §7 replace "no melting or mass loss yet;" with "melting, film and spraying per §9 (no vaporisation, no oxide skin, no droplet tracking);" and "P1 elements with 2 mm surface resolution (Step 3's melt layer needs the prism-layer mesh);" with "P1 elements with 2 mm surface resolution and, for melting runs, four prism layers from 0.25 mm;". Append:

''')
emit(block(open(os.path.join(HERE, "assumptions_step3.md")).read().lstrip("\n"), "markdown"))
emit(r'''
- [ ] **Step 4: Append the amendments to the spec**

Change the spec's status line to `Status: implemented 2026-09-21 (plan docs/superpowers/plans/2026-09-20-melt-spraying.md); amendments in §18` and append:

```markdown
## 18. Amendments (implementation, 2026-09-20/25)

Measured while writing and executing the plan; each overrides the section it names. The plan's "Measured facts and
spec amendments" list carries the numbers.

0. §7, §16.5 — DRAMA's own gamma is constant: `atmosphereData.xml` and the packaged `StaticEnvironmentData.csv` carry
   gamma = 1.4000 at every altitude from 0 to 150 km, so the Knudsen identity coefficient sqrt(gamma pi/2) = 1.4829 is a
   constant and the freestream gamma never varies along a flight. SESAM's own mean free path, recovered from the reference
   CSVs, is the hard-sphere value with d = 3.65 A (ratio 1.0001), which `aero.mean_free_path` already implements.
1. §13.1 — the reference runs are `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind` and
   `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind` (the wrapper adds no material suffix for the default
   material); SESAM hollows the melting sphere at fixed outer geometry (`thick_mm` = shell thickness; heat input,
   radiation, Kn and C_D of the intact sphere until the last gram), so `--removal instant` removes the liquid of every
   element as it forms and keeps the geometry — no area scaling. The end-time metric is the interpolated time at which
   the body's mass (film excluded) falls below 1 % of the initial.
2. §6 — the enthalpy is nodal with a lumped capacity matrix and an enthalpy-consistent Newton update (the element-mean
   secant iteration cycles across a latent-heat ramp); conductivity is not scaled by φ_e (a nearly consumed element must
   keep conducting); patch owners die at φ_e ≤ 0.05, interior elements keep ≥ 10⁻³, deaths cascade within a step,
   material-free nodes are pinned. The FEniCSx backend assembles only the stiffness in UFL; the radiation and all loads
   are nodal in both backends (the boundary moves). Both backends agree to 1e-10 in mass.
3. §8 — every active element feeds its liquid inventory (nodal mean of the ±2 K feed ramp at the liquidus), owners to
   their patches by area, interior elements to the four nearest patches; what leaves is the element's *molten* part, at
   the enthalpy that part carries (the feed-weighted mean of its nodal enthalpies), with the difference to the element's
   mean returned to its nodes over the next step — `--removal instant` books the mass at h_liquid instead, which is the
   lumped device SESAM's Q/L_f law is compared against. Runoff is a
   linearly implicit upwind scheme (4 sub-steps, exact conservation and steady state; a wetting front advances one patch
   per sub-step) that never crosses the equator; leeward films are static and may stay attached.
4. §5 — the layers are built by radial projection of the inner gmsh sphere's boundary (gmsh's extrudeBoundaryLayer is
   not usable on the OCC sphere): 46 278 nodes / 255 276 tets / 15 430 patches by default; the box mesh is a Kuhn split.
5. §7 — the boundary layer is Ranger's integral with the local edge velocity (Thwaites' momentum thickness is a
   constant 12.3 × smaller within ±1.7 %); G includes the deceleration term toward the nose; near the rim the
   modified-Newtonian expansion to p∞ makes Kn_δ ≥ 0.1 on a continuum body — which amendment 16 then removed, the
   rarefied rim having been an artefact of the pure Newtonian pressure rather than a property of the flow.
6. §9 — a film thicker than δ_m takes the thick branch in every regime (with the local shear's film velocity); the thin
   branch's rate is ρ_l min(b, λ*/8)/τ* with τ* = 0.798 λ*^1.5 (ρ_l/Σ)^½ (Girin & Kopyt's Table 1 mass rate), the
   λ_t cut-off never limits it; droplets are capped by the film on the patch and by R/4. Spraying is melt-limited
   everywhere: both branches strip far faster than the surface melts, so the film stays microns thin.
7. §13.2 — tiers: exact (GI, φ_cr with We_cr 4.62; ≤ 2 %), integrated (t_f, N, r_med with the ambient-density Reynolds
   number in δ_a; ≤ 30 %), reported (t_s.d.: half his for iron, far less for stone; ranges; σ; z₀). His GI mixes the
   ambient density (We∞) with the compressed one (Re∞); his α uses the ambient density.
8. §13.3 — Girin & Kopyt's Table 1 implies an effective dynamic pressure 7.8 × ρ₂V₀²; with it r_d, τ_d and their mass
   rate ρ₁r_d/(2τ_d) are reproduced within 1 %; Table 2's λ* column is 10 × smaller than their Eq. (14) (τ* agrees).
9. §13.5 — the layer variants are 2 / 4 / 6 (six layers from 0.125 mm; eight layers of 0.25 mm with growth 2 would
   exceed the radius).
10. §10, §12 — extra history columns `removed_mass_kg`, `film_thickness_max_mm`, `film_thickness_mean_mm`,
    `film_T_max_K`, `film_T_mean_K`, `film_frozen_fraction`, `unapplied_load_J`, `n_dead_elements`; `runoff_mass_kg` is
    the mass that arrived on another patch; the source table has 22 columns.
11. §14 — `--k-scale` (verification device) and `--consistent-mass` (replacing `--lumped-mass`, whose sense is now the
    default) join the CLI.
12. §13.6 — the cost target is restated per macro step, because a flight's length is no longer bounded: about 2 s on
    the default mesh (0.6 s of it conduction, the rest the melt step's Cantera edge states, the runoff solve and the
    contiguous-layer walk). A flight that demises still fits the spec's 6 minutes — the 50 mm physics case is 130 s in
    415 steps — but one whose remnant survives does not, the 100 mm physics flight running to the ground in 1191 steps
    and 2357 s (about 40 min), so bound it with `--t-max` when only the spraying phase is wanted.
14. §6.4/§6.5 of the Step 1 spec (aero) — SESAM's factor-2 drag step at Ma 1 is applied as a cubic ramp over Ma 0.98–1.02:
    a discontinuous C_D stalled the integrator for a light remnant at its transonic terminal velocity; the Step 1 reference
    flights are unchanged.
15. §7 (decided 2026-09-22, replacing the Kn_delta regimes) — the flow-regime gate has three branches, decided on the body
    Knudsen number with SESAM's own mean free path and the freestream Mach number: distinct shock layer (Kn_body < 0.01 and
    Ma > 1), merged (0.01-10) and free molecular (>= 10 or Ma <= 1). Stage one asks whether a distinct bow shock exists, not
    whether the flow is continuum, because a continuum construction must not certify itself; the 0.01 threshold is the
    standoff criterion (Delta/R 0.08-0.14 against a shock 3-10 mean free paths thick), cross-checked by Re2 = rho_inf V R /
    mu(T0) merged below ~100. The gate selects the MELT CLOSURE, never whether spraying happens: Girin's dispersion relation
    has no gas parameters, and tau_w is defined in every branch (Schaaf-Chambre), so the melt is sheared throughout. Within
    the shock-layer branch a wall Knudsen number Kn_local = lambda_w/R -- geometric length, wall state (p_w, T_wall) --
    selects Girin's Eq. (2) closure below 0.01 and the Couette closure V_s = tau_w b / mu_melt above it. The body gate is a
    DECLARED CONSERVATISM, not a physical deduction: Kn_body > 0.01 does not imply Kn_local > 0.01 (measured
    lambda_w/lambda_inf = 1/167 and Kn_local = Kn_body/84 at the nose), and every step records kn_body, kn_local_stag,
    re_shock, flow_branch and the closure fractions so the cost is measurable. Nothing in the merged or free-molecular
    branch asserts that the patch sees only freestream density; the free-molecular branch evaluates the loads from
    freestream conditions with no compression model.
16. §7 (pressure) — modified Newtonian is replaced by modified Newtonian to the sonic point plus a Prandtl-Meyer expansion
    beyond it, blended over 10 degrees above phi* and floored at p_inf; Prandtl-Meyer supplies only p_w, with the state
    coming from Cantera's isentropic expansion of the post-shock reservoir and u_e from energy conservation. gamma_pm is a
    configuration parameter (`--gamma-pm`, default 1.15) carrying a factor ~2. No oblique-shock or entropy-swallowing
    machinery is built: a mass-flux balance at the shoulder gives y_s/R ~ 0.04-0.1 where the local shock inclination is
    85-88 degrees (M_n = 0.996 M), indistinguishable from normal, because a sphere is all nose and the entropy layer is
    swallowed only many nose radii downstream -- which is also what justifies the normal-shock reservoir.
17. §13.1 (decided 2026-09-22) — the 100 mm sphere at 70 km is the marginal case, not a clean anchor: Kn_body = 0.0098 sits
    exactly on the gate while Re2 = 175 says distinct shock. 60 km (Kn_body 0.0026, Re2 599) is the first unambiguous
    anchor. At the 77.5 km break-off every sweep diameter is merged (5 mm 0.596, 20 mm 0.149, 50 mm 0.0596, 100 mm 0.0298),
    so Girin-certified results exist only for the 100 mm case below ~69 km -- measured as 69 % of its melting steps and 84 %
    of its sprayed mass; everything else is the flagged extension.
18. §8, §10 (decided 2026-09-22) — **the film has a temperature and it re-solidifies.** The film is thermally thin
    (q b/k_l = 0.22 K across a 10 µm film at 2 MW/m², b²/α = 0.3 ms against the 0.5 s macro step), so it is given no
    energy equation: its mass is handed to the solver, which carries it on the boundary nodes with the **liquid** heat
    capacity (tangent c_p, secant of the liquid enthalpy), and it is at the surface temperature by construction. Three
    consequences for §8 and §10. (a) A kilogram of film holds the *liquid* enthalpy h(T) + L_f(1 − f_l(T)), not the
    mixture h(T): a film is liquid by construction, so it carries its latent heat wherever it sits and feeding it
    debits the latent heat its mass has not yet paid. Booked at the mixture enthalpy, mass could be relabelled from
    solid to film for free wherever the surface sat on the melting ramp — measured: a 100 mm sphere turned entirely to
    film on a quarter of its latent heat. (b) The enthalpy is the patch's **nodal mean**, matching the nodal sum the
    capacity matrix carries; the enthalpy of the mean temperature differs from the mean of the enthalpies by the whole
    latent heat across the ramp and left a −0.4 residual in the coupled balance. (c) Every transfer books only the
    enthalpy difference it carries, on the nodes the mass arrives at — mass that moves at one temperature books nothing.
    Booking the two halves separately closes the balance just as exactly and wrecks the field (±m h/12 on every surface
    node; the body swung to 1618 K and −1202 K in two macro steps). Re-solidification is the mirror of the feed rule
    (the fraction 1 − f_feed(T_patch) of a patch's film returns to its owner element), **netted** with the feed per
    element: run as two independent rates they cycled 3 % of the body's mass through the film every step with no net
    effect. Declared limits: φ_e is capped at 1, so film whose owner element is full or dead stays liquid — the mesh
    cannot grow a crust outside itself, and `film_frozen_fraction` reports how much film the enthalpy calls solid; a
    thick crust would conduct, which a lumped nodal capacity does not represent. Droplets leave at the film's own
    temperature, so they carry the surface's superheat (measured: +2.4 % on h_liquid for the 50 mm flight, +0.6 % for the
    100 mm one), not h_liquid exactly.
19. §6, §10 (decided 2026-09-22) — **deferred melt loads are bounded and the Newton update follows the node's own
    mixture.** No node is asked to move more than `LOAD_DT_MAX` = 1000 K in one macro step against its current
    capacity (material + film): a node whose own mass has melted away has nothing to heat with the enthalpy its melt
    delivered, and uncapped such nodes were asked to move 1e4–1e5 K in a step. The remainder waits in `pending_load`
    (reported as `unapplied_load_J`, ~2 % of the absorbed heat; 17 % if the cap is 100 K) and is counted as dropped if
    the node dies first, so the balance is exact either way. The enthalpy-consistent Newton update inverts the node's
    own material-plus-film enthalpy (`enthalpy_mixed`/`temperature_from_enthalpy_mixed`, weighted by the film's share
    of the node's mass), because a film-owned node has a shorter latent plateau: inverting the material's h(T) there
    settled into a period-3 limit cycle between 777 K and 864 K. The linear solve falls back to one fresh AMG hierarchy
    and then a direct solve when CG stalls on the emptied interior (φ ~ 1e-3 with unscaled conduction).
20. §8, §12 (measured 2026-09-22) — **film on one patch can stop being a film, and the diagnostic now says so.** The
    area-mean film is microns to tens of microns, but `film_thickness_max_mm` reached 627 mm on the 100 mm physics
    flight (38 of 1192 steps above 20 mm) and 2431 mm on the 50 mm one. Stage-by-stage instrumentation shows every
    such sample falls in the last 5 % of a flight, while the body collapses (3026 patches to 56 on the 50 mm case),
    always on a small surviving facet (0.10–0.27 x the median area) whose melt would be a 4–14 mm blob on a
    0.5–0.8 mm facet: `m_f/(rho_l A)` is then a volume per area, not a depth. The concentrators are the death
    hand-over, the dead-element feed to the nearest patches, and the runoff into facets the graph cannot drain (only
    windward-windward edges exist, an outflow needs t_hat . n_ij > 0, and non-finite coefficients are zeroed).
    Two changes follow. (a) `_kill` spreads a vanished patch's film over the NEAREST_PATCHES nearest survivors by
    area, as `_add_to_film` already does for interior elements — the raw peak ratio falls 2431 -> 981 mm, and because
    film then reaches four owner elements instead of one, re-solidification in the heat-and-cool test rises from
    0.135 kg to 0.456 kg with the stranded fraction falling from 5 % to zero. (b) `film_thickness_max` reports only
    facets with b <= sqrt(A) (as well as the sliver filter) and the rest is reported as `film_blob_fraction`, a new
    history column: the reported maximum falls to 1.672 mm (50 mm) and 4.21 mm (100 mm) with no step above 20 mm, and
    the blob fraction gives the honest exposure — nonzero in 17 of 408 steps on the 50 mm flight and 110 of 217 on the
    100 mm one (median 0.131 of the film, at most 0.21 % of m0). At flight level nothing much moves: the 100 mm case's
    sprayed mass goes 1.0243 -> 1.0227 kg and its median droplet radius 101.9 -> 101.0 µm. The concentration itself is
    declared, not fixed: `b` still feeds `lubrication` (q cubic in b — self-limiting, since the runoff rate then goes
    as b^2/A), `spray.evaluate` (We_s linear in b, bounded afterwards by the droplet-radius caps, r <= 690 µm measured)
    and `rayleigh_taylor` (quadratic in b; reported only as of this amendment -- applied since amendment 25). Clamping those branches to b_eff = min(b, sqrt(A)) is
    defensible but changes droplet sizes on those facets, so it is left as a modelling decision. Amendment 16's "the
    film pile-up is gone" refers to the rarefied-rim mechanism of pure modified Newtonian, not to these facets.
21. §8, §9 (decided 2026-09-22/23) — **the instability branch is chosen on the depth of liquid, not on the film
    account, and melt only leaves an element where the shear can reach it.** `b = m_f/(rho_l A)` is the mobilised film;
    the liquid the gas shears is the film plus the *contiguous* fully molten material beneath the patch, found by
    marching inward and stopping at the first element that is not fully molten, so melt blocked by solid is never
    credited. Measured on the 100 mm flight: conjugate depth 215-303 um, film 1.7-3266 um (thinner than it in 77 % of
    steps), contiguous layer 350-890 um (thicker in 100 %). The branch test therefore put 79-99.8 % of patches on the
    1994 thin-film mode where Girin's 2017 thick mode applies; with the layer deciding it, the thick branch fires on a
    median 26.7 % of patches (up to 97.9 %) -- lower than a non-contiguous measure suggests, because most patches have
    mushy rather than liquid material beneath them. The layer changes the branch test and the reported-only
    Rayleigh-Taylor criterion and **nothing else**: Girin's surface Weber number stays `rho_l V_s^2 delta_m / Sigma` on
    the conjugate depth, as do his wavelength, period and radius; the thin branch's Weber number is a diagnostic on b;
    the droplet Weber number is on the droplet diameter. Separately, `MeltSettings.feed_depth` (default `conjugate`)
    lets melt leave an element only if a wall-owner or within the conjugate depth of the nearest patch, with
    wall-owners only where the wall Knudsen gate denies the Girin closure; `--removal instant` is unaffected. Effects:
    100 mm median droplet radius 130.7 -> 178.2 um (+36 %) with sprayed mass -4.0 %; 50 mm median 232.0 -> 211.2 um
    (-9 %), demise 203.5 s / 69.32 km -> 207.5 s / 68.08 km, sprayed +2.8 %; molten material held in place in the mesh
    up tenfold to ~1 % of the body. New history columns `molten_depth_max_mm`, `molten_depth_mean_mm`,
    `delta_m_mean_um`, `thick_branch_fraction`.
22. §8 (measured 2026-09-23) — **runoff and stripping are not double counting the shear.** Wall shear stress is a flux,
    not a stock: 200 Pa against a gas momentum flux of order 2e5 Pa, with the shear work ~1/25000 of the gas kinetic
    energy flux and droplet surface creation 1.3 % of the shear work. The shear drives a mean flow in the sheared
    sublayer (runoff) and the instability is a perturbation on that flow's free surface (stripping) -- the mean and the
    fluctuation, not two claims on one budget. The only legitimate coupling is the droplet momentum sink
    `mdot_strip V_s` in `tau_gas = tau_wall + h dp/dx + rho h a + mdot_strip V_s + tau_wave`, whose relative size is
    `Pi = v_melt delta_m / nu` (5 % at 0.1 mm/s of recession, 48 % at 1 mm/s): **report Pi and flag Pi > 0.2**. Not yet
    implemented, and the form the deferred pressure-response work should take: the runoff flux needs two depth limits
    in one expression, the shear-driven part over min(delta_m, h) and the pressure-gradient and deceleration part over
    the full contiguous liquid depth h. The thick branch has the first and uses the film for the second.
23. §9 (decided 2026-09-24) — **Girin's stability criterion gates every branch, and a wave must fit inside its own melt
    region.** The thin branch had no threshold at all: Girin & Kopyt's inviscid side mode is unstable at every wavelength
    (their "unlimiting potential intensity" of high-frequency disturbances), so it stripped film at any angle however
    small the shear — 2.23 g inside 10 deg on the 100 mm flight, 0.23 % of the sprayed mass, from a wave 6.1 x wider than
    its own facet. Girin's We_s > We_cr now gates **every** branch, with We_s = rho_l V_s^2 min(delta_m, layer)/Sigma: his
    closed-form phi_cr evaluated against the model's own local flow, which stays valid on an eroded body where a single
    critical angle does not. Separately `wave_fits` requires one whole wavelength to fit inside the **contiguous molten
    region** the patch belongs to (`MeltingBody.region_extent`) — the region and not the facet, because a facet edge is a
    bookkeeping boundary and using it would make the limit mesh-dependent. Measured on the 100 mm flight: the smallest
    angle releasing anything rises from a median 6.32 to 14.08 deg against Girin's published 16.1-16.3 deg, the spraying
    area halves (61.2 to 33.5 cm2 median) and the unstable share of wetted area inside 10 deg falls about ninefold, while
    the sprayed mass moves by only **-0.65 %** (0.98165 to 0.97531 kg) — because the release is supply-limited, so the film
    leaves through whichever patches remain unstable. Droplets become fewer and larger: count -16 %, median radius +5.4 %
    (178.4 to 187.9 um). The 50 mm sphere, regime 1 throughout, is unmoved in demise (207.5 s) and in sprayed mass to five
    decimals, with only the median radius shifting +6.3 %.
24. §9 (measured 2026-09-23/24) — **the Rayleigh-Taylor criterion, corrected twice, and the deceleration that drives it.**
    Amendment-era fact 11's "the RT criterion is inactive at our 10-30 m/s2" was wrong twice over: the deceleration taken
    from the trajectory runs **6 to 95 m/s2** (median 30), and the criterion fires in 137 of 217 steps of the 100 mm
    flight. `W h^2 rho_l > 3 Sigma` is also the wrong *kind* of test — it rearranges to h > lambda*/2 pi, a statement that
    the pool is deep enough for the fastest mode, not that the pool is unstable. The missing condition is lateral: from
    n^2 = (W k - sigma k^3/rho) tanh(k h) the interface is unstable only for wavelengths longer than
    lambda_c = 2 pi sqrt(sigma/(W rho)), and a pool of lateral extent L whose rim pins the interface admits only
    k >= pi/L, which is what `rayleigh_taylor_bounded` evaluates. Two further corrections followed. The driving is the
    deceleration **normal to the film**, Girin & Kopyt's W sin(Theta) = W cos(phi) on this body — whole at the stagnation
    point, zero at the equator where it lies in the surface — and passing the whole deceleration over-drove every patch off
    the nose by 1/cos(phi); with the normal component the depth criterion fires in 118 of 217 steps over a median 21.6 % of
    the film mass instead of 137 steps and 34.8 %. And the lateral extent belongs to the **connected molten sheet**, whose
    rim is what pins the interface, not to the subset of patches that already satisfy the depth test: measuring it the
    second way used the criterion's own output to define its domain. On the sheet the extent is 125 mm, essentially the
    whole windward face, and the bounded criterion is then satisfied for a median **69 %** of the film mass rather than
    0 %. The conclusion is unchanged and much stronger — the bounded growth time is **4614 ms** against the shear mode's
    **0.170 ms**, a ratio near **27 000**, because the admitted waves are the marginal ones near the equator where
    W cos(phi) is small. Recorded every step: `rt_mass_fraction`, `rt_region_mm`, `rt_wavelength_over_nose`,
    `rt_bounded_fraction`, `rt_bounded_growth_ms`, `rt_growth_ms`, `spray_growth_ms`. One caveat is declared and not
    implemented: an equivalent diameter of 125 mm on a body of about 100 mm transverse diameter cannot host a 25 cm wave,
    because the surface curves away; capping the extent at the transverse diameter would move k_1 only from 25.1 to
    31.4 m^-1 and would not change the conclusion.
25. §9, §14 (decided 2026-09-24) — **the front-surface Rayleigh-Taylor mode is applied, not merely reported**
    (`--rt-spray on|off`, default on; `off` reproduces every run before that date). Girin & Kopyt assign the front surface,
    where V_0 -> 0 and sin Theta -> 1, to their aperiodic Eq. (14) rather than to the side mode, and once amendment 23's
    critical angle gates the shear branches the stagnation cap has no other mechanism at all. The gate needs **both**
    criteria: the bounded form has no threshold of its own — with a wide molten region k = k_max is always admissible and
    tanh(k h) > 0 at any depth, so requiring only that made every wet patch unstable, which a unit test caught at once —
    while the depth criterion `W cos(phi) h^2 rho_l > 3 Sigma` is what confines the mode to the nose, since the normal
    deceleration falls as cos(phi) so the depth required rises as 1/sqrt(cos phi) while the melt gets shallower. Three
    modelling decisions are the model's and not theirs, because their Table 2 carries no mass-rate column: the rate carries
    their side-surface construction across as mdot = rho_l min(b, lambda*/2)/tau*, which reduces to rho_l b/tau* because
    lambda* is centimetre-scale (21 mm at 95 m/s2) against a millimetre-scale b — the film is shaken off within one growth
    time, which is what their text describes, and no new closure constant appears; the droplet radius is lambda*, their
    front-surface r_d, capped as on every branch by the film mass on the patch and a quarter of the body radius, and here
    the film-mass cap binds at 1.2-1.9 mm; and two modes of one interface cannot both break it up, so the **shorter growth
    time takes the patch**. Measured on the 100 mm flight: the mode releases **61.04 g, 6.27 %** of the sprayed mass, out
    to 20 deg and **exactly zero beyond it**, 99.1 % of it inside 15 deg and 58.3 % inside 10. tau_RT is nearly constant at
    10-13 ms across the cap while tau_shear rises from 0.27 ms at 30-45 deg to **35.6 ms at 0-2 deg**, so they cross
    between 2-4 and 4-6 deg. Total sprayed mass moves only **-0.16 %** — supply-limited again — but the size tail is
    transformed: the per-step largest droplet goes from a median 388.6 to **2027.0 um** and the largest anywhere from 1551
    to **4093 um**, which is the inertial "shaking off" Girin & Kopyt associate with meteoroid flares. **For Step 4 the
    droplet source is therefore two populations**, a ~190 um shear population from the 20-90 deg annulus and a ~1.7 mm
    inertial one from the stagnation cap. It is a large-body mechanism: on the 50 mm sphere demise, sprayed mass and
    droplet count are unmoved and only the largest droplet changes (1359 to 2567 um), because the threshold
    h > sqrt(3 Sigma/(W cos(phi) rho_l)) is an **absolute length** — 3.36 mm at 95 m/s2 — that a smaller pool never reaches.
26. §9 (2026-09-25) — **regime 2: Girin's Kelvin-Helmholtz cell** (`spray.BRANCH_REGIME2 = 4`). Girin (2017) classifies the
    melt instability in two stages and the model now applies both, in his order. Stage one asks whether the *theoretical*
    melt boundary layer delta_m of his Eq. (2) can form in the liquid that is present: where it cannot, layer <= delta_m,
    the rigid core still stabilises the disturbances and the case is his dominant ablation — **regime 1**, the
    Girin & Kopyt (1994) thin branch, for which he states their results "are thus valid". Where it has formed the profile
    is the conjugated pair (delta_a in air, delta_m in the melt, independent of the core) and stage two decides the
    mechanism on the kinematic viscosities: nu_melt > nu_gas puts the thicker viscous layer in the melt, with V_s << V_a
    and a profile "close to the discontinuous one (tangential discontinuity)", which is classical Kelvin-Helmholtz —
    **regime 2**; nu_gas > nu_melt makes the thicknesses comparable and the profile "inflated", which is the gradient
    instability of his Eq. (1) — **regime 3**, the existing thick branch. The plan's thin and thick branches were his
    regimes 1 and 3 all along, and 2 was the empty cell. Regime 2 reuses Girin & Kopyt's lambda* and tau* (their Eqs. 11
    and 12), the same mode as regime 1 because it is the same profile, with their rigid-wall factor cth(Lambda) -> 1 in the
    deep-film limit; the release is the thick branch's torus form mdot = rho_l pi r^2/(lambda* tau*) with r = k_r lambda*,
    which does not reference the film depth, and that transfer is the branch's only modelling choice — no new closure
    constant. Measured: liquid aluminium has nu_melt = 5.42e-7 m2/s against an edge nu_gas of 0.0178-2.04 m2/s over the
    100 mm flight to 110 s and 0.0951-3.75 m2/s over the whole 50 mm flight, a ratio of **3.3e4 to 6.9e6**, so this body is
    in regime 3 wherever the layer forms — **0 of 533 960** windward wet patch-steps on the 100 mm flight and 0 of 42 268 on
    the 50 mm one take regime 2 — and a melt viscous enough for regime 2 would fail stage one here in any case, since
    Eq. (2) gives delta_m/delta_a = ((nu_melt/nu_gas)^2 rho_l/rho_e)^(1/3) = **151** at the threshold, about **1.06 m**,
    which is the effect Girin records for his stony variant as the boundary layer becoming "comparable with the radius of
    the meteoroid remnant". The branch is therefore inert on these flights and exercised only by unit test; the 50 mm
    flight is reproduced within the model's own run-to-run floor (worst column 3.1e-8 with the branch against 1.6e-8
    between two runs of a single build, `thick_branch_fraction` identical in all 67 comparable steps, demise unmoved at
    207.5 s).
13. §10, §16.7, §17.5 (decided 2026-09-21) — the body Knudsen number uses the equivalent diameter of the remaining mass and
    the stagnation radius of the heating and the surface flow is fitted to the current windward cap (a least-squares
    sphere through the patch centroids within (1 − cos 30°) R_t of the front-most point, bounded to [0.1, 1.67] × R_t):
    the front erodes fastest and flattens, so the fitted radius reaches the cap within the first 15 % of mass loss and
    the stagnation heating falls to ≈ 0.8 × its sphere value; the sphere drag tables are kept. `--size-feedback initial`
    (D₀, R₀) is SESAM's convention and the verification devices' setting. History columns `nose_radius_mm`,
    `transverse_radius_mm`, `fitted_nose_radius_mm`.
```

- [ ] **Step 5: Append §17 to the facts note and update the package docstring**

Append to `sesam_verified_facts.md`:

```markdown
## 17. Melting (re-measured 2026-09-20 on the two drama-AA7075 US76 references)

- SESAM melts the sphere from the inside: during the melt `thick_mm` is the wall thickness of a hollow sphere of the
  original outer radius (m = ρ 4π/3 [R₀³ − (R₀ − thick)³] to four digits at every row), the convective heat stays
  A₀ × 0.27471 × q_DKR(R₀) × hot-wall (24–26 kW for the 100 mm sphere), `rad_cooling_W` stays −372 W = 4πR₀² εσ 850⁴,
  `knudsen` keeps D₀ and `drag` stays 0.913 until the mass is zero. The removal rate is (Q_conv + rad_cooling)/L_f
  (0.96–0.98 × Q_conv/L_f at the printed sampling), L_f = 400 kJ/kg.
- Onsets and ends (US76, winds off): 100 mm 71.005 km at 43.55 s → mass 0 at 66.476 km / 67.315 s (23.8 s at the
  melting temperature); 50 mm 77.104 km at 178.55 s → 73.094 km / 191.55 s (13.0 s). `demised (uncritical)` with
  `final mass 0 kg`.
```

In `reentry_model/__init__.py` change the docstring's first line to `"""First-principles re-entry model of a solid sphere, Steps 1-3: trajectory, coupled 3D heat transfer, melting and` and the second to `melt spraying (specs under docs/superpowers/specs/)."""`.

- [ ] **Step 6: Final checks and commit**

```bash
"$PY" -m pytest -m "not drama and not reference" -q
git add README.md docs/model_assumptions.md docs/superpowers/specs/2026-09-20-melt-spraying-design.md reentry_model/__init__.py
git commit -m "Document Step 3: README, assumptions, spec amendments (Step 3 Task 15)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

(The facts note lives outside the repository; edit it in place, no commit.)

---

## Self-review notes (writing the plan, 2026-09-21)

- Spec coverage: §5 mesh → Task 1; §6 materials and melting → Tasks 2–3; §7 surface flow → Task 5; §8 film → Task 6; §9 spraying and dispersion → Tasks 4, 7; §10 accounting, §11 coupling → Tasks 9–10; §12 outputs → Tasks 10–13; §13.1 → Tasks 11, 14; §13.2–13.3 → Task 8 (with the amended tiers); §13.4 analytic checks → Tasks 3 (Stefan, fractions/loads balance), 5 (Ranger/Thwaites), 6 (strip, conservation), 4 (dispersion), 9 (balances, death), 3/9 (backends); §13.5 → Task 14; §13.6 → restated per macro step (~2 s; amendment 12); §14 CLI → Task 13; §15 environment → no new packages; §16–17 → Task 15 documents them.
- Every code block is the prototype file that passed the tests named in its task; the appended test blocks are verbatim tails of the prototype test files. Names used across tasks were checked against the prototype: `MeltSettings/MeltingBody/REMOVAL_NAMES` (Task 9 ← Task 13), `SurfaceFlow/RAREFIED_SHEAR_NAMES` (Task 5 ← 9, 13), `SprayModel/K_R/K_T/SOURCE_COLUMNS/histogram/N_BINS/BIN_EDGES` (Task 7 ← 9, 10), `DispersionTable/WE_CRITICAL_PRACTICAL` (Task 4 ← 7, 8, 13), `MELT_COLUMNS/write_particles` (Task 10 ← 13), `has_melt/melt_metrics/plot_melt/MELT_PLOT_NAMES` (Task 12 ← 13, 14), `animate_film/still_marks` (Task 12 ← 13), `Material.MATERIAL_NAMES/liquid/feed_fraction/h_liquid` (Task 2 ← 9, 13), `mesh.DEFAULT_LAYERS/DEFAULT_LAYER_THICKNESS` (Task 1 ← 13), `StepResult.Q_extra/Q_dropped`, `set_fractions`, `pinned` (Task 3 ← 9).
''')
