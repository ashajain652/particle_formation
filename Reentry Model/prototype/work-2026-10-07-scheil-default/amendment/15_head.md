## Amendment of 2026-10-07 — the rigid substrate, the Scheil default and the repeated step series in the README, the assumptions and the spec

> Part of the rigid-substrate amendment of 2026-10-06 (facts 78–87 in `00-shared-context.md`) and of the amendment of
> 2026-10-07 (the Scheil default, the freeze-back fix and the repeated time-step series; facts 88–96). These changes
> stack on the blocks of 2026-09-27, 2026-10-02, 2026-10-03 and both of 2026-10-05 above, which still apply except where
> a change below replaces their wording. The blocks are prose, so re-run `prototype/plan3/audit_docs.py`'s checks
> against them after these edits: `--rigid-substrate` must exist in `cli.py`, and `nonrigid_depth_mean_mm`,
> `slurry_thick_fraction`, `rigid_thin_fraction` and `slurry_held_mass_kg` in `coupled.MELT_COLUMNS`; and add to its
> `SUPERSEDED` list the two wordings change 2 retracts, "`--material AA7075_range` (default with melting" and "both
> carry the liquid properties ρ_l 2400 kg/m³, μ_l 1.3 mPa s, Σ 0.86 N/m".

**README block (`## Physics model — reentry_model (Step 3 ...)`).**

1. *The first example.* Replace the comment "# the model proper: physics heating, AA7075 with its 750-908 K melting
   range, film + runoff + Girin spraying, videos" with "# the model proper: physics heating, AA7075 on Scheil's curve
   between 750 and 908 K (the default with --melt on), film + runoff + Girin spraying, videos".
2. *Options paragraph, the material.* Replace "`--material AA7075_range` (default with melting; DRAMA's `drama-AA7075`
   tables with the alloy's solidus 750 K / liquidus 908 K, latent heat 400 kJ/kg spread across the range) or `AA7075`
   (DRAMA's single 850 K, a ±2 K numerical ramp) — both carry the liquid properties ρ_l 2400 kg/m³, μ_l 1.3 mPa s, Σ
   0.86 N/m (pure aluminium near the liquidus, `reentry_model/data/materials/*.json`);" with: "`--material
   AA7075_scheil` (the default with melting since 2026-10-07: DRAMA's `drama-AA7075` tables with the alloy's solidus 750
   K and liquidus 908 K, the latent heat of 390 kJ/kg released along Scheil's non-equilibrium curve — half of it in the
   13 K below the liquidus — and a liquid surface tension of 0.80 N/m for fully liquid 7075), `AA7075_range` (the same
   range with 400 kJ/kg released linearly and Σ 0.86 N/m), `AA7075` (DRAMA's single 850 K, a ±2 K numerical ramp, Σ 0.86
   N/m) or `AA7075-empiricaldata` (`AA7075` with Σ 0.80 N/m and the heat capacity above 850 K made explicit) — all with
   ρ_l 2400 kg/m³ and μ_l 1.3 mPa s (`reentry_model/data/materials/*.json`). The run name does not carry the material:
   give runs of different materials their own `--name` or `--outdir`;".
3. *Options paragraph, the switch.* After the `--molten-cascade on|off (...)` entry of the 2026-10-03 amendment insert:
   "`--rigid-substrate on|off` (default on: Girin's thin branch needs a rigid substrate, so the regime test compares his
   conjugate depth with the film, the deep liquid and the slurry beneath the wall — everything more than half liquid,
   down to the first point at most half liquid — and off his closure a film on slurry is not sprayed by the shear modes;
   `off` reads the fully molten depth, as every run before 2026-10-06 did);".
4. *Outputs.* Add `nonrigid_depth_mean_mm`, `slurry_thick_fraction`, `rigid_thin_fraction` and `slurry_held_mass_kg` to
   the list of history columns, and `nonrigid_depth` to the surface frame's fields.
5. *Findings.* Append the bullet: "The thin spraying branch needs a rigid substrate (2026-10-06). Girin's thin branch is
   his dominant ablation, where the unmelted core stabilises the disturbances; on an alloy with a 158 K melting range
   the material under the film is mush, and mush more than half liquid is a slurry that flows, not a rigid core. The
   regime test therefore measures the non-rigid depth — the film, the deep liquid, and the material above the
   half-liquid temperature (829 K on the linear law, 895.1 K on Scheil's curve) continuous with the wall, read from the
   temperature field along the inward normal rather than counted in whole elements — and takes the thick branch where it
   exceeds the conjugate depth. Off Girin's closure a film on slurry is held from the shear modes; the front-surface
   Rayleigh–Taylor mode is left to shed it (decided 2026-10-07). On the 100 mm flight to 120 s at the default step the
   rule leaves the masses within a few per cent (sprayed +1.2 % on the linear law, +0.6 % on Scheil's curve) but makes a
   third fewer droplets (−32 % and −26 %), almost all of them the thin branch's; on the 50 mm flight, which never has
   Girin's closure, two-thirds of the sprayed mass leaves by the Rayleigh–Taylor mode, as millimetre droplets from film
   piled on a few facets, and the median radius by mass is 7.4 times larger (plan facts 78–87)."
6. *Findings.* Append the bullet: "The time step, with the rigid substrate on Scheil's curve (2026-10-07). The
   switched-step series of 2026-10-05, repeated on the new default — 0.5 s until the flight enters the continuum regime
   at 49.5 s, then 0.05, 0.025, 0.0125 or 0.00625 s to 120 s — now nearly converges. From 0.0125 s the sprayed mass and
   the mass at 120 s move within the scatter between seeds (1.165 and 1.164 kg; 0.307 and 0.308 kg), the share of the
   mass the thick branch sprays settles at 71 % (71.2 and 71.4 %, against a scatter of 0.3 points) and the median
   droplet radius by number at about 69 µm; the droplet count rises by 2.7 % at the last halving, to 79 million over
   49.5–120 s, and the median radius by mass drifts by under 1 % per halving (187 to 186 µm), from the thick branch's
   own droplets. The default step is not adequate for the droplet population: at 0.5 s the thin branch is almost absent
   and the flight makes 13 million droplets instead of 79 million, although its masses are within about 5 %. Runoff
   stays minor at every step (at most 13 % of the film formed inside any latitude cap leaves it). The rule and the
   material changed together, so how much of the convergence each brings was not measured. Until the fine step is a
   model option, quote the droplet population only from a fine-step run, with its step (plan facts 88–96)."

**Assumptions block (`## 9. Melting and the melt film (Step 3)`).**

7. In the bullet "Latent heat, enthalpy method", replace "f_l is linear between the solidus and the liquidus
   (`AA7075_range`: 750–908 K, ASM Handbook) or a ±2 K numerical ramp around DRAMA's single 850 K (`AA7075`); L_f = 400
   kJ/kg" with "f_l follows Scheil's non-equilibrium curve between the solidus and the liquidus (`AA7075_scheil`, the
   default since 2026-10-07: 750–908 K, ASM Handbook; partition coefficient 0.4, pure-aluminium melting point 933 K;
   half of L_f released in the 13 K below the liquidus; L_f = 390 kJ/kg), is linear across the same range
   (`AA7075_range`, L_f = 400 kJ/kg) or is a ±2 K numerical ramp around DRAMA's single 850 K (`AA7075`, 400 kJ/kg)".
8. In the bullet "Liquid properties", replace "Σ 0.86 N/m — pure aluminium near the liquidus (Smithells; Assael et al.
   2006; ASM Vol. 2)" with "Σ 0.80 N/m for the Scheil default (Bainbridge & Taylor 2013, commercial 7075 by sessile
   drop) and 0.86 N/m for the other variants (pure aluminium near the liquidus; Smithells; Assael et al. 2006; ASM Vol.
   2)".
9. In the bullet "Spraying: Girin's three regimes", after "(layer ≤ δ_m) the rigid core still stabilises the
   disturbances and the case is **regime 1**" insert: "— with the layer measured down to rigid material: since
   2026-10-06 the liquid depth of this test is the non-rigid depth, the film and the deep liquid plus the material above
   the half-liquid temperature continuous with the wall, because slurry is no rigid core (`--rigid-substrate`); off
   Girin's closure a film whose wall is above that temperature takes neither regime 1's mode nor its rarefied
   extrapolation, and only the front-surface Rayleigh–Taylor mode can shed it (2026-10-07)".
10. Replace the bullet "The droplet population and the macro step (measured 2026-10-02, 2026-10-03 and 2026-10-05)" of
   the 2026-10-05 (runoff flux) change 4 with: "- **The droplet population and the macro step (measured 2026-10-02,
   2026-10-03, 2026-10-05 and 2026-10-07).** Measured on the 100 mm physics flight to 120 s. Before the rigid substrate,
   on the linear range, shortening the step — uniformly (2026-10-02/03) or only once the flight enters the continuum
   regime (2026-10-05: 0.5 s until 49.5 s, then 0.05 to 0.00625 s) — never converged the droplet population: over
   49.5–120 s the count rose from 88 to 141 million and the thick branch's share of the sprayed mass fell from 47 % to
   21 % at the fine steps, still moving by 12 % and 4.6 points at the last halving, because the branch test's molten
   depth, counted in whole elements and read at the end of the conduction step, let the step decide which branch took
   the mass. With the rigid substrate on Scheil's curve (2026-10-07) the same series converges from 0.0125 s to within a
   few per cent: the sprayed mass (1.165 and 1.164 kg) and the mass at 120 s (0.307 and 0.308 kg) within their scatter
   between seeds, the thick branch's share of the sprayed mass at 71 % within 0.3 points and each branch's mass fixed to
   2 g, the median radius by number at about 69 µm; the count rises by 2.7 % at the last halving, to 79 million, and the
   median radius by mass drifts by under 1 % per halving, from the thick branch's own droplets (median by number 187 to
   184 µm at the finest step). Each branch makes nearly the same droplets at every fine step (about 187 µm on the thick
   branch, 66 µm on the thin, by number). Which of the two changes, the rule or the material, does the converging was
   not measured. At the default 0.5 s step the thin branch is almost absent — 13 million droplets instead of 79 million
   — while the masses are within about 5 %. Status: converged within a few per cent from 0.0125 s; the droplet count,
   the branch split and the median radius by mass are not to be quoted from the default step, and the re-solidified mass
   not at all (it counts freeze-and-re-melt cycles) (plan facts 58, 61, 69–77 and 88–96)."

**Spec amendments block (`## 18. Amendments`).** Append:

11. "31. §8, §9 (decided 2026-10-06 and 2026-10-07, measured 2026-10-06) — **the thin branch needs a rigid substrate.**
   The regime test's liquid depth is the non-rigid depth: the film, the deep liquid, and the material above the
   half-liquid temperature (`Material.T_rigid`; 829 K linear, 895.1 K Scheil) continuous with the wall, followed along
   the inward normal through the elements the line crosses and read inside each from the P1 field, each element's
   stretch weighed by φ_e; it replaces amendment 21's contiguous molten depth in the lubrication branch, the spray's
   regime test and the thick branch's Weber depth, while the Rayleigh–Taylor criteria keep the liquid layer. Off Girin's
   closure a film on slurry takes no shear mode, and the front-surface Rayleigh–Taylor mode is left to shed it.
   `--rigid-substrate on|off`, default on, `off` reproducing every earlier run bit for bit; history columns
   `nonrigid_depth_mean_mm`, `slurry_thick_fraction`, `rigid_thin_fraction`, `slurry_held_mass_kg`. On the 100 mm flight
   to 120 s at the default step: sprayed +1.2 %, droplets −32 %, the thin branch's share of the sprayed mass 10.4 to 2.5
   %, the thick share of the wet windward patches under the closure 27 to 63 % over a non-rigid depth of a median 8.1 mm
   (linear law); on Scheil's curve the slurry band is 12.9 K instead of 79 K and the depth 1.45 mm. On the 50 mm flight,
   never under the closure, the Rayleigh–Taylor mode releases two-thirds of the mass (plan facts 78–87)."
12. "32. §6 (decided 2026-09-27, applied 2026-10-07) — **Scheil's curve is the melting default**, `AA7075_scheil` with
   Step 4's 390 kJ/kg and 0.80 N/m; the linear range stays selectable. With it and the rigid substrate the switched-step
   series (0.5 s to the continuum onset at 49.5 s, then 0.05 to 0.00625 s) converges from 0.0125 s to within a few per
   cent: the masses within their scatter between seeds, the thick branch's share of the sprayed mass at 71 % within 0.3
   points, the median radius by number at about 69 µm, the droplet count within 2.7 % of the finest step's 79 million
   over 49.5–120 s; the default 0.5 s step makes six times fewer droplets, the thin branch being almost absent there.
   The deep runoff's freeze-back (amendment 27) now debits the film by its own part, not by a difference that rounded
   below zero and stopped a fine-step run; no earlier flight is changed by that (plan facts 88–96)."

