# Viscous-film instabilities for the thin branch (parked 2026-10-22)

This note records why the instabilities of a viscous film sheared by gas came up, what is known about them so far,
and how to investigate them when this is picked up again. Nothing here is implemented. The literature descriptions
below come from memory, not from reading the papers, and must be checked against the sources before they are cited.

## How this came up

1. **The thin branch had no stability threshold of its own.** The thin branch covers films shallower than Girin's
   conjugate depth δ_m (about 0.1 mm on the 100 mm flight). It uses Girin and Kopyt's 1994 side-surface mode
   (J. Aerosol Sci. 25, 1353). Since 2026-09-24 it had borrowed Girin's 2017 critical surface Weber number as its gate
   (3.08 since the critical-Weber amendment of 2026-10-22, facts 108–110, commit `fa02a4c` on branch
   `we-critical-308`).

2. **Girin and Kopyt's equation has no critical Weber number.** Their Eq. 5 was solved numerically following Girin's
   procedure: the largest growth rate over all wavelengths, plotted against the gas Weber number on the film depth.
   The growth rate is positive at every Weber number. In supersonic gas the cause is wave drag, which destabilises
   every wavelength down to nanometres. In subsonic gas it is ordinary Kelvin–Helmholtz instability, where long waves
   always grow. A threshold appears only when the wavelength is bounded, by the molten region's size at the long end
   and their viscous-dissipation rule at the short end. That threshold does not depend on the film depth, and the
   model's existing wave-fit test already applies it. Details are in
   `prototype/work-2026-10-22-thin-threshold/REPORT.md` on branch `we-critical-308`.

3. **The borrowed gate was replaced.** The 3.08 gate is being removed from the thin, rarefied and regime-2 branches.
   The wave-fit test on the dominant wavelength λ* now decides, together with a viscous-damping gate. The work is in
   `prototype/work-2026-10-22-thin-gate/` on branch `we-critical-308`. Its comparison flight against the 3.08 flight
   was not finished when this note was written.

4. **The viscous-damping gate.** The first version used Lamb's damping rate for a wave on deep liquid, 2νk². Here ν is
   the melt's kinematic viscosity (about 5.4 × 10⁻⁷ m²/s for liquid aluminium) and k is the wavenumber. The
   recommended replacement adds the damping in the boundary layer at the rigid floor (Hunt 1952; Mei's *The Applied
   Dynamics of Ocean Surface Waves*):

   Γ = 2νk² + k·√(νω/2) / sinh(2kh),

   where h is the depth of non-rigid material down to the rigid substrate and ω is the capillary frequency on that
   depth, from ω² = (σk³/ρ)·tanh(kh). A patch releases nothing where Γ·τ* ≥ 1, meaning damping outpaces the
   growth time τ*. A patch is also rejected where the liquid is thinner than the floor's boundary layer,
   h < √(2ν/ω).

5. **Why that reject matters.** Below the floor's boundary layer, viscosity reaches through the whole liquid and
   there is no inviscid core, which Girin and Kopyt's derivation needs. Their mode does not apply there. That does
   not mean the film is stable: sheared viscous films have instabilities of their own. The reject therefore marks
   where the model has no valid mechanism. The prompt for the gate asks for the count and mass of rejected patches to
   be recorded, so that what is excluded can be seen.

6. **Measured so far (wecr308_d100, thin releases only).** Using the source table's `b_m`, believed to be the film
   thickness:
   - k·b runs from 0.7 to 5 (1st to 99th percentile), with a median of 2.
   - The film is 10 to 54 times thicker than the floor's boundary layer, with a median of 22.
   - Γ·τ* runs from 0.07 to 0.23, with a maximum of 0.87, against the gate's limit of 1.

   So neither the damping gate nor the reject would have stopped any released mass. However, patches that the old
   gate or the wave-fit test held back are not in the source table, so very thin films near `B_MIN` may reach the
   reject once the borrowed gate is gone.

## The instabilities to investigate

These are recalled from memory and still to be checked against the papers.

- **Yih (1967), "Instability due to viscosity stratification", J. Fluid Mech. 27.** Two fluid layers of different
  viscosity are sheared together. The shear stress is continuous across the interface, so the velocity gradient
  jumps there, and a displaced interface produces a tangential-stress disturbance that feeds the wave. The mechanism
  is purely viscous, needs no inviscid core, is strongest for waves much longer than the layer, and grows slowly at
  low Reynolds number. Surface tension suppresses short waves. Yih (1963, Phys. Fluids) is the related
  gravity-driven film on an incline, not the sheared case.
- **Miesen and Boersma (1995), "Hydrodynamic stability of a sheared liquid film", J. Fluid Mech.** This is the
  configuration closest to ours: a liquid film on a rigid wall, sheared by gas flowing over it. As recalled, they
  found a long-wave mode (Yih's viscosity-jump mechanism) and a short-wave "internal" mode driven by shear inside the
  liquid, which is distinct from Kelvin–Helmholtz.
- **Boomkamp and Miesen (1996).** This paper classifies instabilities in two-phase parallel flow by where the wave's
  energy comes from:
  - the viscosity jump at the interface (Yih);
  - the gas Reynolds stress (Miles' critical layer);
  - the shear inside the liquid (the internal mode);
  - the gas pressure (Kelvin–Helmholtz).

  It would let us state which mechanism each spraying branch assumes.

## What they would mean for the model

- **They are the physics for the rejected films.** Films thinner than the floor's boundary layer, which Girin and
  Kopyt's mode cannot describe, are covered by the Yih and Miesen–Boersma modes.
- **They are precursors, not release rules.** In annular two-phase flow, long-wave film instabilities typically grow
  into large disturbance (roll) waves, and droplets are then torn from those waves' crests, often by a
  Kelvin–Helmholtz-type process. A droplet-release rule would need a second step after the instability.
- **Their wavelengths are long.** They are much longer than the film depth, so the wave-fit test against the molten
  region's size would usually be the binding check, especially near the nose.
- **They may not matter.** On the one flight measured, no thin-film release was in the viscous-film regime. Whether
  these modes need modelling depends on the mass the reject's counter reports.

## Plan for when this is picked up

1. **Find out whether it matters.** Read the reject's count and mass from the comparison flight in
   `prototype/work-2026-10-22-thin-gate/` (and, if needed, the 50 mm flight). If the rejected mass is negligible,
   record that in the facts and in the thesis as a stated limitation, and stop here.
2. **Get the sources.** Check whether Yih 1967, Miesen and Boersma 1995 and Boomkamp and Miesen 1996 are in the
   Literature Review folder. If not, find them. Confirm every claim in the section above before using it.
3. **Map the regime.** For the rejected patches, compute the dimensionless groups those papers use: the viscosity
   ratio (gas to melt, about 0.12 at the edge state), the density ratio (about 10⁻⁷), the film Reynolds number, the
   film depth relative to the gas boundary layer, and a Weber number based on the film. Check whether our parameters
   lie inside the ranges the papers solved; the extreme density ratio may not be covered.
4. **Solve the stability problem if needed.** If the papers do not cover our regime, set up the Orr–Sommerfeld
   problem for a film on a wall under a gas boundary layer (their configuration) and solve it numerically, the way
   `gk_dispersion.py` solves Girin and Kopyt's Eq. 5. Find the growth rate and wavelength of the fastest mode, and
   whether there is a critical condition.
5. **Decide what to do with the result.**
   - (a) Keep the reject and report the excluded mass as a limitation.
   - (b) Add a branch that applies the viscous-film mode as a precursor, with a release rule taken from the
     annular-flow entrainment literature (disturbance waves and crest stripping).
   - (c) Show that the mode grows too slowly to act within the flight, so the reject is physically right.

   Each needs Asha's decision, with the measured mass in hand.
6. **Record it** under the Step 3 sub-plan conventions: a dated amendment with facts in `00-shared-context.md`, the
   code in a work copy with its diffs, and the README and assumptions entries via sub-plan 15.

## Pointers

- Thin-threshold analysis: `prototype/work-2026-10-22-thin-threshold/` on branch `we-critical-308` (`REPORT.md`,
  `gk_dispersion.py`, `scan_we.py`).
- Thin-gate amendment and comparison: `prototype/work-2026-10-22-thin-gate/` on the same branch (`compare_thin.py`,
  `launch.sh`, `meas/`).
- The 3.08 baseline flight: `prototype/work-2026-10-07-dt-continuum/runs/wecr308_d100` in the main checkout (outputs
  git-ignored).
- The spraying code: `reentry_model/spray.py` (`thin_film_mode`, `wave_fits`, the gate function, `SprayModel.evaluate`).
