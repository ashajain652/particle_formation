# Sub-plan: Task 15 — Documentation: README, assumptions, spec amendments, facts note

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 7489–end of file). Read `00-shared-context.md` first. This extract includes everything after the Task 15 heading, since the rest of the source document (physics-model docs, the amendments log, the re-measured melting section, and the self-review notes) is itself the content this task is writing and checking.


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 15.** The documentation blocks are prose, so nothing but this note stops them describing a
model that has moved on — which is exactly what `prototype/plan3/audit_docs.py` exists to catch. Re-run it after
these edits and add the superseded wording to its `SUPERSEDED` list.

1. **The conservation claim is weakened and must say so.** Facts 28 and 31 quote balance residuals of −1.7·10⁻¹⁰ and
   +2.2·10⁻¹⁰. After remeshing these become **exact by construction following an explicit correction**, with a
   reported pre-correction residual per remesh. That is the standard position for an ablation code and defensible,
   but it is a weaker statement, and `README.md`, `docs/model_assumptions.md` and the three Task 15 blocks must all
   carry the change together.
2. **`PHI_DEATH` 0.05 → 0.50** wherever the documentation states it, with fact 44's measurement as the reason.
3. **Prism layers are no longer the default**; the dense band is. State fact 12's sensitivity as the justification
   and fact 39's 6 mm recession quantum as the problem being solved.
4. **The seven assumptions of spec §11** go into `docs/model_assumptions.md` with their verification status: the
   geometry is *fitted* to the mass account, not identical to it; the conduction solves on the original geometry and
   what that costs; smoothing rounds real features as well as noise; remeshing conserves only after correction; a
   fragmenting body is not modelled; lateral resolution is unchanged; and it all remains a fixed-attitude
   calculation, which fact 24 already identifies as the dominant uncertainty on the shape.
5. **Every claim carries its measurement** (CLAUDE.md). The new numbers are facts 38–45.
---

**Depends on:** everything — this is the last task, written after every other plan's results exist.
**Produces:** the README, the model-assumptions section, the spec's amendment log, and an external facts note that lives outside this repository.
**Character:** pure documentation.
**Read before implementing:** this task's own draft amendment currently restates the runtime target as "about 4 minutes on the default mesh (target 6 minutes)," but that number only describes a flight that fully demises — it does not mention the roughly 40-minute non-demising case recorded in measured fact 23 (see the shared context). Resolve this before finalizing the documentation's stated target, so Task 14 and Task 15 don't end up stating two different numbers. Also confirm the external facts-note file's path and write permissions before delegating this task, since it edits a file outside this repository. Carry the fixed-attitude caveat on the 29.9%-mass finding into the write-up, exactly as Task 14 must.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan for the documentation work itself (what to write, where, and how it will be checked for consistency against the other 14 tasks' actual results).

---

### Task 15: Documentation — README, assumptions, spec amendments, facts note

**Files:**
- Modify: `README.md` (a new Step 3 section before "## Tests"; the Step 2 verification table refreshed; the Tests section's test counts and the reference-tier line), `docs/model_assumptions.md` (§9 appended, §7 amended), `docs/superpowers/specs/2026-09-20-melt-spraying-design.md` (§18 "Amendments" appended; status line), `reentry_model/__init__.py` (docstring: Steps 1–3), the facts note `/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md` (§17 appended)

- [ ] **Step 1: Re-run the Step 2 verification with the new thermal core and refresh its README table**

```bash
"$PY" analysis/reentry_model_thermal_verification.py
```

Expected (measured 2026-09-21 with the nodal-enthalpy core): d100 sesam 0.46 % / 2.23 % / +0.31 % / +0.30 % / 24.6 K (1.20 %) / 6.66 % / 2089 K at 139 s / 74 s; d050 sesam 0.37 % / 2.36 % / +0.03 % / +0.02 % / 28.0 K (1.22 %) / 3.40 % / 2442 K at 245 s / 24 s; physics ratios 0.744 / 0.770, peaks 2319 K at 118 s and 2494 K at 239 s, 83 s / 33 s. Replace the numbers in the README's Step 2 verification table with the measured ones (the runtimes halve: `energy()` and the operators are cheaper and the Newton takes 2.0 iterations) and add one sentence after the table: "Re-measured on 2026-09-21 with Step 3's nodal-enthalpy core (lumped capacity matrix): every metric within 0.1 K / 0.01 % of the Step 2 values; runtimes halved."

- [ ] **Step 2: Add the Step 3 section to `README.md`**

Insert before `## Tests`. The verification table below carries the prototype's own measurements; replace each row with Task 14's `summary.md` where the re-run differs, and write the sensitivity findings from `reentry_model_output/verification_melt/sensitivity/sensitivity.md` into the paragraph that names the fifteen variants:

```markdown
## Physics model — `reentry_model` (Step 3: melting, melt film and melt spraying)

Step 3 turns the heated sphere into a particle source. On top of the Step 2 coupling, every macro step now also
melts the body (enthalpy method with the latent heat in the nodal enthalpy, element fractions φ_e and element death
on a mesh with four 0.25/0.5/1/2 mm prism layers under the surface), feeds the liquid to a film on the surface
patches, moves the film with a lubrication runoff driven by the gas shear and the pressure gradient, and strips it
into droplets by Girin's gradient instability (thick films: Girin 2017; thin films: Girin & Kopyt 1994; a rarefied
extrapolation of the thin mode where the local Knudsen number Kn_δ ≥ 0.1). Droplets are recorded at birth in a source
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
`--feed-depth conjugate|all` (default `conjugate`: melt leaves an element only where the gas shear can reach it); `--prism-layers 4`,
`--layer-thickness 0.25` (mm, growth 2; 46 k nodes / 255 k tets on the 100 mm sphere); `--demise-fraction 0.01`;
`--particles/--no-particles`; `--size-feedback current|initial` (`current` with `girin`: the body Knudsen number on
the equivalent diameter of the remaining mass and the stagnation radius fitted to the windward cap, bounded to
1.67 × the transverse radius for a flat front; `initial` with `instant`: D₀ and R₀, SESAM's convention); `--k-scale`
(a verification device). Melt and runoff start at the liquidus: material between the solidus and the liquidus holds
its latent heat but counts as solid for the film (spec §8, §17.2). The film has the surface's own temperature (it is
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

Outputs: `<run>.csv` gains `mass_kg` (now varying), `film_mass_kg`, `sprayed_mass_kg`, `runoff_mass_kg`,
`removed_mass_kg`, `melt_front_depth_max_mm`, `equivalent_radius_mm`, `n_active_elements`, `spraying_area_m2`,
`theta_cr_deg`, `n_released`, `released_mass_kg`, `r_median_um`, `r_max_um`, `regime_fraction_continuum/_slip/_fm`,
`rt_active`, `removed_enthalpy_J`, `film_thickness_max_mm`, `film_thickness_mean_mm`, `nose_radius_mm`,
`transverse_radius_mm`, `fitted_nose_radius_mm`, `n_dead_elements`; the run
JSON adds the melt onset, spraying onset, demise, the masses, size statistics and the settings; `<run>/particles.npz`
(the source table: time, altitude, velocity, θ, patch centroid, regime, branch, film thickness, δ_m, We_s, radius,
count, mass, release velocity and direction, We_d, Oh, breakup flag), `particles_summary.csv`, `size_distribution.csv`
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
- The film is stripped as fast as it melts: both instability branches remove hundreds to thousands of kg/m²/s where
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
```


Also update the Tests section: the unit-tier line becomes `# unit tests (~6 min; the thermal solver, coupled and melting tests dominate)` and the reference-tier line `# the Step 1 reference flights (~15 min), the Step 2 coupled runs (~20 min) and the Step 3 bookkeeping runs (~1 min)`; and in the `fenicsx_env` paragraph "the seven conformance tests" → "the eight conformance tests (Step 2's seven and the Step 3 melting cross-check)".

- [ ] **Step 3: Append §9 to `docs/model_assumptions.md` and amend §7**

In §7 replace "no melting or mass loss yet;" with "melting, film and spraying per §9 (no vaporisation, no oxide skin, no droplet tracking);" and "P1 elements with 2 mm surface resolution (Step 3's melt layer needs the prism-layer mesh);" with "P1 elements with 2 mm surface resolution and, for melting runs, four prism layers from 0.25 mm;". Append:

```markdown
## 9. Melting and the melt film (Step 3)

- **Latent heat, enthalpy method.** The nodal enthalpy h(T) = ∫c_p dT + L_f f_l(T) carries the latent heat; f_l is
  linear between the solidus and the liquidus (`AA7075_range`: 750–908 K, ASM Handbook) or a ±2 K numerical ramp
  around DRAMA's single 850 K (`AA7075`); L_f = 400 kJ/kg and the solid tables are DRAMA's `drama-AA7075` verbatim.
  The capacity matrix is lumped (diag Σ φ_e V_e/4 ρ c_i), so the discrete energy is exactly ρ∫h dV with h
  interpolated linearly, and each Newton step is mapped through the true h(T) per node so that a node cannot jump
  across the melting range (the element-mean formulation of Step 2 could not carry a 4 K ramp under a 30 K nodal
  spread). The Stefan front on a box is reproduced within 0.3 %.
- **Liquid properties** (DRAMA has none): ρ_l 2400 kg/m³, μ_l 1.3 mPa s, Σ 0.86 N/m — pure aluminium near the
  liquidus (Smithells; Assael et al. 2006; ASM Vol. 2); alloy corrections within 10 %; no oxide skin. They are
  constants: the film now has a temperature (below), but its viscosity and surface tension are still evaluated at the
  liquidus, which is where most of it sits. μ_l falls by roughly a third per 200 K of superheat, so the runoff and the
  thin-branch wavelength are the quantities a temperature-dependent μ_l would move.
- **How deep the shear reaches, and what may therefore leave.** The gas shear penetrates the liquid only to Girin's
  conjugate depth, measured at 215 to 303 µm on the 100 mm flight, while the contiguous molten layer beneath the wall
  is 350 to 890 µm. Two rules follow. Melt leaves an element only if the shear can reach it: a wall-owning element
  always, a buried element only within the conjugate depth of the nearest patch, and where the wall-Knudsen gate
  denies Girin's closure, only wall-owning elements. And the thick-versus-thin instability branch is decided on the
  depth of *contiguous* liquid under the patch — the march inward stops at the first element that is not fully molten,
  so melt blocked by solid is never credited — while the film account remains what can actually move and be stripped.
  The layer decides the branch and the reported-only Rayleigh–Taylor criterion and nothing else; Girin's surface Weber
  number, wavelength, period and droplet radius are all on the conjugate depth, as he has them.
- **Runoff and stripping happen at once, and that is not double counting.** Wall shear stress is a flux, not a stock:
  the shear sets up a mean flow in the sheared sublayer, which carries mass along the wall, and the instability is a
  perturbation on that flow's free surface, which removes mass from it. The only coupling is the momentum the departing
  droplets carry, of relative size (melting speed × conjugate depth / kinematic viscosity) — 5 % at a recession speed
  of 0.1 mm/s, 48 % at 1 mm/s. What is *not* yet done, and is deferred with the pressure response of the molten
  region: the pressure-gradient and deceleration-driven part of the runoff flux should integrate over the whole liquid
  depth while the shear-driven part stays limited to the conjugate depth.
- **What flows.** Melt and runoff start at the liquidus: an element's material becomes film in proportion to a ±2 K
  ramp at the liquidus (f_feed), so the mushy range holds latent heat but neither runs off nor is stripped
  (conservative; a coherency-point treatment is a future iteration). What leaves an element is its *molten* part, at
  the enthalpy that part carries (the f_feed-weighted mean of its nodal enthalpies), and what stays behind is the
  colder rest — so the element is debited the difference, which is what makes the melt front advance at the
  energy-limited rate whatever the element size. With `--removal instant` the mass leaves the body at the liquidus
  enthalpy instead, and the difference is returned to the element's nodes.
- **The film's temperature and re-solidification.** The film is thermally thin (q b/k_l = 0.22 K across a 10 µm film
  at 2 MW/m², 22 K even at 1 mm; b²/α = 0.3 ms against the 0.5 s macro step), so it is given no energy equation: its
  mass is carried on the boundary nodes with the liquid heat capacity and it is at the surface temperature by
  construction. It therefore heats, cools and freezes with the surface, and droplets leave carrying the surface's
  superheat (+2.4 % on h_liquid for the 50 mm flight, +0.6 % for the 100 mm one) rather than at the liquidus exactly. Because a film is liquid by
  construction it holds the liquid enthalpy h(T) + L_f(1 − f_l) wherever it sits, so feeding it costs the latent heat
  the mass has not yet paid and no mass can be relabelled from solid to liquid for free. Re-solidification is the
  mirror of the feed rule — the fraction 1 − f_feed(T_patch) of a patch's film returns to its owner element as solid,
  and the two directions are netted into one transfer per element so that the film is never churned. Two limits are
  declared: φ_e is capped at 1, so film whose owner element is full (or has died) stays liquid for good — the mesh
  cannot grow a crust outside itself — and `film_frozen_fraction` reports how much film is in that state; and a thick
  crust would conduct, which a lumped nodal capacity does not represent. A re-solidified crust is therefore still
  available to run off and be stripped, which is conservative for mass loss.
- **How melt energy is booked.** Mass moving at one temperature books nothing; mass that arrives somewhere hotter or
  colder books the enthalpy difference it carries, on the nodes it arrives at. No deferred load may move a node more
  than 1000 K in one macro step (a node whose own mass has melted away has nothing to heat with it); the remainder
  waits and is reported as `unapplied_load_J`, and is counted as dropped if the node dies first. The coupled energy
  balance is exact (1e-9 of the absorbed heat) with the queue and the dropped loads in it.
- **Element fractions and death.** φ_e scales an element's heat capacity, not its conductivity (a thinner sliver of
  the same material conducts better, not worse; scaling k isolated the surface nodes and drove them to thousands of
  kelvin). A patch owner dies at φ_e ≤ 5 % (its remainder joins the film); interior elements keep ≥ 10⁻³ so that no
  cavity opens; nodes without material are pinned. The exposed faces of a dead element become patches (a transient
  pit until its neighbours die: with Lees' distribution the pit walls at θ ≈ 90° receive almost no heat, so the
  heating artefact is small; the SESAM-equivalent mode is defined by its total and is unaffected).
- **Mesh.** Four prism layers (0.25, 0.5, 1, 2 mm) built by radial projection of the inner 2 mm/8 mm gmsh sphere's
  boundary triangulation (exact on a sphere; each prism split into three tetrahedra by the smallest-node-id diagonal
  rule); the surface triangles are 2.16 mm (h_surface × R/(R − 3.75 mm)). 46 k nodes / 255 k tetrahedra for the
  100 mm sphere; 0.5–0.7 s per macro step.
- **Gas-side surface flow.** Modified-Newtonian pressure p_e = p∞ + (p_s − p∞)cos²θ on the windward face and the
  isentropic expansion of the equilibrium stagnation state to p_e (Cantera, 1° bins); leeward p_e = p∞, no shear.
  Girin's linear-profile boundary layer in Ranger's (1972) form with the actual edge velocity, δ_a² = 58.1 ν_e ∫u_e⁴ds/u_e⁵
  (exact for potential flow; Thwaites' momentum thickness is a constant 12.3 × smaller within ±1.7 %); shear
  τ_c = μ_e u_e/δ_a. Regimes by Kn_δ = λ_e/δ_a: continuum < 0.01, slip 0.01–0.1 (Maxwell first-order slip,
  u_e/(1 + Kn_δ)), transitional/free-molecular ≥ 0.1 (τ_fm = ρ∞V² sinθ cosθ; `--rarefied-shear bridged` blends the two
  with SESAM's f(Kn) instead). Near the rim the expansion to p∞ makes ρ_e tiny and Kn_δ > 0.1 even on a continuum
  body — a property of the pressure model that the local criterion inherits. The driving gradient
  G = 2(p_s − p∞) sinθ cosθ/R − ρ_l a sinθ includes the body's deceleration (the film is pushed toward the nose).
- **A film deeper than its patch is wide is not a film.** On a collapsing body the melt of thousands of dying patches
  is concentrated onto a few small survivors — by the death hand-over, by the interior-element feed to the nearest
  patches, and by a runoff graph that has no edge out of the last windward ring — and m_f/(ρ_l A) stops being a depth:
  in the last 5 % of the 50 mm flight a 0.5–0.8 mm facet carries melt that would be a 4–14 mm blob. The reported
  thickness therefore covers only facets with b ≤ √A, and the rest is reported as `film_blob_fraction` (nonzero in 17
  of 408 steps on the 50 mm flight, up to 0.999 of the film in the collapse; in 110 of 217 on the 100 mm one, median
  0.131 — the nose crater, where τ and G ∝ sinθ leave the film nothing to drive it out — and at most 0.21 % of m₀). The mass and the energy stay exactly accounted; what is lost is
  the lubrication picture, and with it the meaning of the runoff and spraying branches on those facets. A patch graph
  that merges slivers, or a film carried per element rather than per facet, is the fix; it is not in this iteration.
- **Film and runoff.** One film mass per patch; lubrication velocity and flux with a thin (b ≤ δ_m) and a thick
  (b > δ_m) branch; runoff by a linearly implicit upwind scheme on the patch graph (4 sub-steps per macro step,
  exact conservation, exact steady state; a wetting front advances one patch per sub-step), never across the equator;
  leeward films are static and can stay attached (their fate is Step 4's). The film is fed where it melts and stripped
  there within the step, so it stays microns thin except where the runoff piles it at the windward rim.
- **Spraying: Girin's three regimes.** His (2017) classification is applied in two stages. δ_m is the melt velocity
  boundary layer his Eq. (2) *predicts*, so the liquid depth decides whether that layer can form at all: where it cannot
  (layer ≤ δ_m) the rigid core still stabilises the disturbances and the case is **regime 1**, Girin & Kopyt's (1994)
  side mode — λ* = 1.5 M_e Σ/(ρ_e u_e²), r = λ*/4, τ* = 2 capillary periods, rate ρ_l min(b, λ*/8)/τ*, their Table 1
  mass rate. Where it has formed, the profile is the conjugated pair and the kinematic viscosities pick the mechanism:
  ν_melt > ν_gas gives a near-discontinuous profile and classical Kelvin–Helmholtz (**regime 2**, Girin & Kopyt's λ* and
  τ* in their deep-film limit with the torus rate; inert for liquid aluminium, which sits 3.3×10⁴ to 6.9×10⁶ from that
  threshold), while ν_gas > ν_melt gives the inflated profile and his **regime 3** gradient instability — the dispersion
  relation solved numerically (Δ_f, Im Ω_f tabulated against We_s; k_r 0.17, k_t 1.1, We_cr 4.62, his constants for
  ordinary liquids), δ_m and V_s from his conjugated-layer relations, stripping ρ_l π r²/(λ_f t_per) per area. The
  rarefied branch uses regime 1's mode on the freestream Mach number and momentum flux (an extrapolation). Two gates
  apply to every regime: Girin's We_s > We_cr, his critical angle evaluated against the model's own local flow, and the
  requirement that one whole wavelength fit inside the contiguous molten region the patch belongs to (the region, not the
  mesh facet). Girin & Kopyt's front-surface Rayleigh–Taylor mode is **applied as well** (`--rt-spray on|off`, default
  on): the body's deceleration is 6–95 m/s², not the 10–30 m/s² first assumed, and the mode needs both its depth
  criterion W cos φ h² ρ_l > 3 Σ, which confines it to the nose, and the laterally bounded form, which supplies the
  admissible wavelength and growth rate. It releases 6.27 % of the sprayed mass on the 100 mm sphere, out to 20° and
  exactly zero beyond, as ~1.7 mm droplets against the shear population's ~190 µm; on the 50 mm sphere it barely fires,
  because its threshold is an absolute depth (3.36 mm at 95 m/s²) a smaller pool never reaches. Where two modes could
  claim one interface the shorter growth time takes the patch. Droplets are capped at the film on the patch and a
  quarter of the body radius; one radius per patch and step; no within-patch size spread; recorded at birth, not
  tracked. Every mode strips far faster than the melt supply, so the mass loss is energy-limited (Girin's "outstripping
  ablation") and the instabilities set the droplet size and the release map rather than the mass.
- **Size feedback** (decided 2026-09-21, replacing the spec's fixed R₀). Mass = Σφ_e ρV_e + film; the drag reference
  area is the current surface's projection; the body Knudsen number uses the equivalent diameter of the remaining
  mass; the stagnation radius for the heating and the surface flow is a least-squares sphere fitted to the current
  windward cap (the patches within (1 − cos 30°) R_t of the front-most point, R_t the transverse radius about the mass
  centre), bounded to [0.1, 1.67] × R_t — the front erodes fastest and flattens, so the fitted radius grows from R₀ to
  the cap and the stagnation heating falls as R^−½ (a flat face heats like a sphere of 1.67 × its radius: its
  stagnation velocity gradient is ≈ 0.6 × a sphere's, Boison & Curtiss 1959). The sphere drag tables are kept. The
  verification devices keep D₀ and R₀ (SESAM's convention).
- **Verification devices.** `--removal instant` (every element's liquid leaves as it forms) and `--k-scale 1e4`
  reproduce SESAM's lumped Q/L_f melting (mass within 1.3 % of m₀, onset within 0.2 km, end within 1.3 %); Girin's
  Table 1 is reproduced in its exact tier (GI, φ_cr) and, with the ambient-density Reynolds number, in t_f, N and
  r_med within 30 %; his spraying durations are not (half his, and far less for the stony variant).
```


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
   modified-Newtonian expansion to p∞ makes Kn_δ ≥ 0.1 on a continuum body.
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
12. §13.6 — the 100 mm physics-mode flight takes ≈ 4 min on the default mesh (target 6 min; 146 s without the size feedback).
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

- Spec coverage: §5 mesh → Task 1; §6 materials and melting → Tasks 2–3; §7 surface flow → Task 5; §8 film → Task 6; §9 spraying and dispersion → Tasks 4, 7; §10 accounting, §11 coupling → Tasks 9–10; §12 outputs → Tasks 10–13; §13.1 → Tasks 11, 14; §13.2–13.3 → Task 8 (with the amended tiers); §13.4 analytic checks → Tasks 3 (Stefan, fractions/loads balance), 5 (Ranger/Thwaites), 6 (strip, conservation), 4 (dispersion), 9 (balances, death), 3/9 (backends); §13.5 → Task 14; §13.6 → measured (146 s); §14 CLI → Task 13; §15 environment → no new packages; §16–17 → Task 15 documents them.
- Every code block is the prototype file that passed the tests named in its task; the appended test blocks are verbatim tails of the prototype test files. Names used across tasks were checked against the prototype: `MeltSettings/MeltingBody/REMOVAL_NAMES` (Task 9 ← Task 13), `SurfaceFlow/RAREFIED_SHEAR_NAMES` (Task 5 ← 9, 13), `SprayModel/K_R/K_T/SOURCE_COLUMNS/histogram/N_BINS/BIN_EDGES` (Task 7 ← 9, 10), `DispersionTable/WE_CRITICAL_PRACTICAL` (Task 4 ← 7, 8, 13), `MELT_COLUMNS/write_particles` (Task 10 ← 13), `has_melt/melt_metrics/plot_melt/MELT_PLOT_NAMES` (Task 12 ← 13, 14), `animate_film/still_marks` (Task 12 ← 13), `Material.MATERIAL_NAMES/liquid/feed_fraction/h_liquid` (Task 2 ← 9, 13), `mesh.DEFAULT_LAYERS/DEFAULT_LAYER_THICKNESS` (Task 1 ← 13), `StepResult.Q_extra/Q_dropped`, `set_fractions`, `pinned` (Task 3 ← 9).
