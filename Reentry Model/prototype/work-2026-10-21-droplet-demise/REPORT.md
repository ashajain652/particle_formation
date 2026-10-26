# Droplet demise model: trajectory, heating and vaporisation of molten AA7075 droplets

> **Moved 2026-10-26.** The code of this folder now lives in `reentry_model/droplet_reentry/` (`droplet_demise.py` became `model.py`; the sweeps, plots and checks are in its `analysis/` and `tests/`), and its outputs are written to `output/droplet_reentry_output/`. This file is kept as the record of the work; script names and `runs/` paths below are as they were.

Session of 2026-10-21. Code and outputs are in `prototype/work-2026-10-21-droplet-demise/`; nothing is committed.
Run outputs are under `runs/demise_sweep/` (git-ignored).

## What the model does

`droplet_demise.py` follows one molten AA7075 droplet from its release with the parent's velocity until it
reaches a steady final state. The steady state is the first moment the droplet's speed is within 1 % of the local
terminal velocity and its temperature is within 1 K of the local equilibrium temperature. The equilibrium
temperature is the one at which heating, radiation and evaporation balance at the current flight conditions. The
final state records temperature, initial mass, mass lost (in total and per element), speed, altitude and time.

- **Trajectory:** the rotating-Earth integrator verified against DRAMA (`reentry_model.trajectory.Simulator`,
  with J2 gravity, Coriolis and centrifugal terms), subclassed to carry the temperature and the element masses.
  Only the drag coefficient is replaced, by the droplet drag law from 2026-10-14, because the validated drag
  table is wrong at low speed.
- **Heating:** the validated physics heating (`reentry_model.heating.PhysicsHeating`) applies in hypersonic,
  not-fully-rarefied flow (Mach number above 6, Knudsen number below 1), and is blended out over Mach 6 to 3 and
  Knudsen 1 to 10. Elsewhere the model uses the exact free-molecule (Schaaf and Chambre) heat exchange, reduced
  towards dense-flow conduction. See "The heating hand-over" below for why.
- **Vaporisation:** aluminium, zinc and magnesium evaporate separately from a well-mixed droplet; copper and
  chromium (1.8 % of the mass) are treated as non-volatile. Each element's partial pressure follows Raoult's law
  (vapour pressure scaled by the element's share of the atoms). The rate is the kinetic (Hertz–Knudsen) rate,
  limited by diffusion through the air only while the summed vapour pressure is below the ambient pressure. Above
  it, the droplet boils and the kinetic rate applies, which covers your point about zinc and magnesium. Each
  element's latent heat is taken out of the droplet's energy.
- **Weber pre-check:** before each run, the release Weber number (air density times speed squared times
  diameter, divided by surface tension) is compared with 12. A droplet at or above 12 is skipped and shaded grey
  in every figure.

The property values are those of the 2026-10-14 sweep:
- density 2,400 kg per cubic metre;
- specific heat 1,180 joules per kilogram per kelvin;
- release temperature 908 K (the alloy's liquidus);
- emissivity 0.1;
- surface tension 0.86 newtons per metre;
- thermal accommodation 0.8, the validated heating's own default (the earlier sweep used 1);
- evaporation coefficient 1;
- composition: aluminium 90.1 %, zinc 5.6 %, magnesium 2.5 %, copper 1.6 %, chromium 0.23 % by mass.

## What was verified

`verify_droplet_demise.py` passes all 12 checks:
- the vapour pressures return one atmosphere at each element's boiling point;
- the evaporation rate in rarefied flow equals the textbook Hertz–Knudsen rate;
- below boiling at sea level, diffusion throttles zinc to 2.6 % of its kinetic rate;
- the Weber check skips a case when the threshold is lowered;
- mass is conserved with evaporation off;
- the element losses add up to the total loss.

The strongest check is the one with evaporation switched off. The new model, with its rotating Earth and reused
heating code, then reproduces the earlier flat-Earth sweep's peak temperature for 10 micrometres at 80 km: 2,818 K
against 2,819 K.

## Results (315 cases: 15 diameters from 1 to 200 micrometres, 21 altitudes from 90 to 50 km)

**No case was skipped.** The largest Weber number anywhere is 1.3, for a 200 micrometre droplet at 62 km, an
order of magnitude below the break-up threshold of 12. **Every case reached a steady state.**

**Vaporisation removes up to 45 % of the droplet.** The loss peaks for releases around 76 km and depends almost
entirely on the release altitude, not the diameter, as the temperature did in the earlier sweep:

| release altitude | 90 km | 82 km | 74 km | 66 km | 58 km | 50 km |
|---|---|---|---|---|---|---|
| share of mass lost | 30 % | 43 % | 45 % | 32 % | 8 % | 6 % |
| peak temperature | 1,690 K | 1,780 K | 1,880 K | 1,910 K | 1,550 K | 908 K |

These are the values at the 1 micrometre end. At 66 km and below, the largest droplets lose a few percentage
points less, because their flow is only partly rarefied.

**Evaporation caps the temperature.** Without vaporisation the earlier sweep reached 3,590 K. With it, the peak
never exceeds about 1,920 K, because evaporation removes heat faster as the droplet gets hotter. That cap also
makes the earlier "unphysical" flags moot: this model is the physical version of those cases.

**Zinc and magnesium are lost first, and almost completely.** Above 58 km, every droplet loses all of its zinc
and magnesium. At 50 km, where the droplet only cools, it still loses 88–89 % of its zinc and about a third of its
magnesium within about a second, because their vapour pressures exceed the ambient pressure at the release
temperature. The final droplet is essentially aluminium with a little copper.

**Final state:**
- **Temperature:** 208 to 272 K, close to the local air temperature.
- **Speed:** the terminal speed, from 0.02 m/s (1 micrometre at 50 km) to 28 m/s. As found before, droplets of
  about 20 micrometres and up released above 80 km meet their 1 % terminal-speed criterion only once it has
  fallen to about 25 metres per second.
- **Time to reach the final state:** from 0.04 seconds to 285 seconds.

## The heating hand-over, and what it costs

At moderate speeds and for hot droplets, the validated heating falls well below the exact free-molecule result.
It is 25 % low for a 2,000 K droplet at Mach 10, and it returns zero instead of cooling once the droplet is hotter
than the air arriving at it. The reason is that it treats rarefied heating as a high-speed limit and subtracts the
heat leaving the surface using equilibrium-air chemistry. Where the flow is only partly rarefied, however, it
captures a real reduction in heating (about 20 % for 200 micrometres at 55 km) that the free-molecule model
cannot. Hence the regime-based hand-over. Where the two are blended, they differ by up to 0.25 of the incoming
energy flux, and that difference is the physics of the transition.

To show what the choice costs, I re-ran three diameters with the validated heating applied everywhere above
Mach 6. The mass lost changes by at most 4.6 percentage points (lower), for 200 micrometres near 75 km.

## Drag law: Singh et al. (2021) is the default from 2026-10-22

The model now uses the general sphere drag correlation of Singh et al. (*AIAA Journal*, 2021, DOI
10.2514/1.J060648; `singh_drag.py`), fitted across continuum, transitional and free-molecule flow at all Mach numbers.
It replaces the earlier composed law, which is kept as `drag_law="blend"` for comparison. `check_singh_drag.py` confirms
that it reproduces Stokes drag, the exact free-molecule drag, and the hypersonic continuum value of 0.9, and is
continuous at Mach 1. The paper defines the post-shock state only above Mach 1; below it the free-stream values are
used, which makes every term continuous at Mach 1.

Re-running both sweeps with it changes the results as follows:
- **Mass lost:** by at most 0.6 percentage points, and by less than 0.1 in all but the largest droplets released at
  66 km and below. During vaporisation the droplets are in fast free-molecule flow, where the two laws are identical.
- **Terminal speeds:** fall by 5 to 11 %, because Singh assumes fully diffuse reflection at low speed (Epstein's
  limit), which gives about 5 % more drag than Millikan's Cunningham fit. The difference grows to 11 % for large
  droplets at 50 km, where the flow is only partly rarefied.
- **The final state** is reached up to 41 % sooner and up to 0.4 km higher.
- **Peak temperatures** rise by at most 44 K.

The numbers in the sections above were produced with the earlier law. The mass-loss and temperature figures in them
stand; the terminal speeds should be read as 5 to 11 % high. The current figures and `final_states.csv` use Singh; the
earlier results are kept with a `_blend` suffix.

## What is assumed and not known

- **Well-mixed droplet.** The model treats the droplet as uniformly mixed. For the larger droplets, liquid
  diffusion (about 2 seconds across 100 micrometres) is slower than the zinc boil-off (milliseconds), so their
  zinc and magnesium losses are upper bounds.
- **Composition at release.** Zinc and magnesium boil at the release temperature. The melt film on the parent
  would therefore already be losing them before any droplet is shed, so the release composition assumed here
  is probably too rich in zinc and magnesium. This is a question for the Step 3 film model.
- **Uncertain inputs.** The evaporation coefficient (1), the vapour pressures of zinc and magnesium far below
  their boiling points (uncertain by tens of percent), and the diffusion coefficient used in dense flow (an
  order-of-magnitude estimate) are all untested. The emissivity was not varied in this sweep.
- **Released into undisturbed air.** As before, each droplet starts in the free stream, not in the parent's wake.
- **Steady state is a snapshot.** It is a quasi-steady state: the droplet keeps settling, so its terminal speed
  and temperature continue to drift slowly as it descends.
