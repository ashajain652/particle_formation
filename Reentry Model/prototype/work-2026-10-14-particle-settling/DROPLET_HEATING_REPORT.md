# Heating of molten AA7075 droplets released into the wake of the descending satellite

Session of 2026-10-14. Everything here lives in `prototype/work-2026-10-14-particle-settling/`; nothing is
committed. The run outputs are under `runs/droplet_sweep/`, which the repository's `.gitignore` already
excludes.

## 1. What was asked, and whether it is done

You asked for a sweep of molten AA7075 droplets over nine diameters (2, 6, 10, 14, 18, 50, 100, 150 and 250
micrometres) and twenty-one release altitudes (90 km down to 50 km in 2 km steps), each released with the
velocity of the parent satellite's Main Body at that altitude, with a lumped thermal model added to the
existing trajectory integrator, and a report of how much each droplet heats while it decelerates. All 189
cases ran without error, and so did three sensitivity sweeps of the same 189 cases (756 runs in total, 341
seconds on all cores). The per-case histories and the summary tables are written. Every verification check
listed in section 3 passes.

Files:

- `droplet_heating_sweep.py` builds the parent state, the thermal model and the sweep. It imports
  `particle_trajectory.py` and `sphere_sweep.py` and calls their functions; it does not modify or copy them.
- `verify_droplet_heating.py` runs the checks in section 3 (11 seconds).
- `summarise_droplet_sweep.py` adds the vaporisation estimates, writes the matrix table and draws the figure.
- `runs/droplet_sweep/matrix_summary_baseline.csv` is the summary table across the matrix, one row per case.
- `runs/droplet_sweep/summary_eps*_alpha*.csv` (and `*_with_vapour.csv`) hold every diagnostic for the
  baseline and each sensitivity sweep.
- `runs/droplet_sweep/histories/droplet_dNNNum_hNNkm.csv` holds the time history of each baseline case.
- `runs/droplet_sweep/droplet_peak_temperature.png` is the figure referred to below.
- `runs/droplet_sweep/parent_states.json` records the release state used at each altitude.

## 2. What was found, in plain language

**Most of this matrix is unphysical, as you expected, and the reason is vaporisation.** Without a
vaporisation model, the droplet temperature rises well past the point where liquid aluminium would boil
away. In 99 of the 189 baseline cases (every case released between 62 and 82 km) the peak temperature
exceeds aluminium's normal boiling point of 2743 K. A further 54 cases (58–60 km and 84–90 km) stay below
2743 K but still exceed the temperature at which aluminium's vapour pressure equals the very low local air
pressure, which is the boiling point that actually applies at these altitudes. Only the 36 cases released
between 50 and 56 km stay below that local boiling point. For the 153 cases above it, the temperature
numbers are not predictions; they say only that the droplet would lose a substantial share of its mass to
vapour, which this model does not represent.

**The droplet diameter makes almost no difference to the temperature.** At every altitude, the peak
temperature of every droplet from 0.01 to 250 micrometres agrees to within 1.4 % (the outlier is the
250 micrometre droplet at 52–56 km; an earlier version of this report said 1 %, which was slightly loose). This is not a
coincidence of the numbers but a property of free-molecule flow (flow so rarefied that air molecules strike
the droplet individually without colliding with each other near it): the heating per kilogram, the drag
deceleration and the radiative cooling per kilogram all scale as one over the diameter, so the temperature
reached at a given speed is the same for every size. What the diameter changes is how *long* and how *far*
the deceleration takes, which grows in proportion to the diameter.

**Peak temperature is highest for release near 70 km, not at the top of the ladder.** A droplet released
at 90 km absorbs more energy per kilogram (it starts faster), but it absorbs it slowly because the air is
thin, and radiation removes most of it as it arrives. At 70 km the droplet is still fast and the air is
about 24 times denser. At the moment of release the incoming heat is about 840 times the radiated power at
70 km, against about 50 times at 90 km, so the 70 km droplet climbs much further before radiation catches up. Below about
60 km the parent itself has already slowed down so much that there is little energy to absorb.

**Every droplet heats and cools within seconds, and has refrozen long before it reaches the ground.** The
heating is complete within 0.4 milliseconds to 23 seconds of release, depending on size and altitude, and
the droplet falls back below its melting temperature within 2 milliseconds to 36 seconds. The fall to the
ground then takes between 4 hours (250 micrometres) and 650 days (2 micrometres).

## 3. What was verified

All of these are in `verify_droplet_heating.py` and all pass.

The free-molecule heat flux on a surface element, which is the building block of the droplet's heating,
reproduces a direct numerical integration over the molecular velocity distribution to a relative error of at
most 7e-16 (machine precision) at four combinations of flow speed and surface angle, including a
surface facing away from the flow. Integrated over the sphere, the heat input approaches the textbook
high-speed limit (accommodation times one half times air density times speed cubed times frontal area) to
within 0.4 % at a speed ratio of 10 and 0.04 % at 30, which is the expected rate of approach. It vanishes
(below 1e-24 W) for a droplet at rest at the air's own temperature, and matches the closed-form result for
pure heat conduction to a resting sphere to every printed digit. Doubling the surface quadrature from 96
to 800 points changes the answer by 1e-15, so the surface integral is converged.

The altitude-keyed parent state agrees with your sanity table at 90 and 80 km. At 70, 60 and 50 km it
differs from your table by 1–4 % in speed (for example, I get 6.835 km/s at 70 km against your 6.91 km/s).
I checked why: your table quotes the parent row just above each altitude (the rows at 70.47, 60.35 and
50.71 km), while my function interpolates to the exact altitude. In each case the interpolated speed and
angle lie between the two bracketing parent rows, and inverting the result through the original
`sphere_sweep.interpolate_state` returns the starting altitude exactly, so I believe the interpolation is
correct and the difference is in the table. I also confirmed that the parent's `velocity` column is the speed
relative to the air: it reproduces the parent's own dynamic pressure column to within 0.8 % at every altitude
in the band, which is what the droplet's release speed relative to the air requires.

The energy books close: the energy absorbed, minus the energy radiated, minus the energy conducted back to
the air, equals the liquid heat capacity times the temperature change, to a relative error of 3e-15. The
hand-off from the coupled thermal integration to your existing `integrate_hybrid` finds the droplet already
on its terminal velocity to a relative mismatch of 5e-7, inside the 1e-6 that function requires.

What was not verified: there is no independent reference (no experiment and no SESAM run) for the
temperature of a micron droplet, so the thermal results are checked for internal correctness only.

## 4. What the model assumes, and where the property values came from

**Molten AA7075 properties.** The density is 2400 kg per cubic metre, the value already chosen for the Step 3
work in `prototype/proto3/reentry_model/data/materials/AA7075.json` (`liquid.density`), which cites the
Smithells Metals Reference Book for pure aluminium near its melting point. The solid value of 2813 is not
used anywhere in this sweep. The liquid specific heat is 1180 joules per kilogram per kelvin, from the
NIST-JANAF value for liquid aluminium (31.75 joules per mole per kelvin, divided by aluminium's molar mass);
this matches your suggested starting value. It is held constant at all temperatures, and the 5.6 % zinc, 2.5 %
magnesium and 1.6 % copper in the alloy are ignored for it, which I expect to be a correction of a few percent.
The release temperature is 908 K, the AA7075 liquidus (the temperature above which the alloy is fully
liquid) from the ASM Handbook, as used in `AA7075_range.json`.

**Emissivity** (the fraction of a perfect radiator's output that the surface emits). I used 0.1 as the
baseline, because clean liquid aluminium measured near 1000–1700 K has a total emissivity of roughly 0.05 to
0.12; the 0.4 in the material file is a solid-surface value. Sensitivity runs used 0.05 and 0.4. An oxidised
droplet (aluminium oxide skin) could be higher still, which 0.4 partly covers. The radiative background is
taken as 250 K; it changes the radiated power by under 1 % at the temperatures that matter.

**Accommodation coefficient** (the fraction of an impinging molecule's excess energy that it leaves with the
droplet). I assumed 1, full accommodation, because that is what the drag law's free-molecule branch already
assumes, and a sensitivity run used 0.8. Air molecules are given two internal (rotational) degrees of
freedom, with vibration frozen.

**Lumped temperature (one temperature for the whole droplet).** I checked rather than assumed this. Taking
the peak heat flux on the front of each droplet, the temperature difference it would drive across the droplet
radius by conduction (liquid aluminium conducts about 90 watts per metre per kelvin) is at most 21 K, in the
250 micrometre case released at 66 km, against a temperature rise of about 2500 K in that case: under 1 %.
For the micron sizes it is below 2 K. The time for heat to diffuse across the radius is at most 0.5
milliseconds, short against the heating times. Internal circulation in the liquid, which is not modelled,
would only make it more uniform.

**Heat input.** I implemented the general Schaaf and Chambre free-molecule heat transfer to a sphere, valid
at any speed ratio and including the energy carried away by molecules re-emitted at the droplet's own
temperature, rather than the high-speed limit. You asked me to check what fraction of the absorbed energy
arrives while the speed ratio (flight speed divided by the air's most-probable molecular speed) is above
3. The answer is at least 99.8 % for every release at 54 km and above, and 93 % for the 50 km releases,
which start at a speed ratio of only 3.8. So by your criterion the high-speed limit would have been nearly
acceptable. However, it would still have been wrong by a large factor, for a different reason: it ignores
the energy carried away by molecules leaving a very hot surface. Integrated over each flight, the high-speed
limit overestimates the absorbed energy by between 14 % and a factor of 5 (median 41 %), the error being
largest when the droplet is hottest relative to the energy of the incoming air. The general expression is
what the results use.

**The droplet's temperature also enters its drag.** The free-molecule drag coefficient depends on the ratio
of the droplet's surface temperature to the air temperature; the coupled integration passes the live
temperature to your `drag_coefficient`, rather than the fixed ratio of 1 the earlier script used.

**Released into undisturbed air.** Each droplet starts at the parent's speed relative to the air, in
undisturbed US76 atmosphere. The real wake behind the parent is hotter, slower and denser than the free
stream for some distance; that is not modelled, and it would reduce the heating of droplets that stay in
the near wake. This is probably the most important physical assumption after the missing vaporisation.

**Not modelled:** vaporisation, mass loss, oxidation, the latent heat released when the droplet refreezes
(400 kilojoules per kilogram, which would extend the cooling time), and any heat-transfer reduction when the
flow is not fully free-molecular (see section 6).

## 5. The evidence

The peak temperatures below are for the 10 micrometre droplet; every other diameter at the same altitude
is within 1.4 %. The local saturation temperature is the temperature at which aluminium's vapour pressure
equals the local air pressure, from the correlation of Alcock, Itkin and Horrigan (1984), which returns
2743 K at one atmosphere. The table shows that the peak temperature exceeds the local saturation temperature
at every altitude from 58 km up, and exceeds even the one-atmosphere boiling point from 62 to 82 km.

| release altitude | release speed | peak temperature (baseline) | local Al saturation temperature | energy absorbed per kilogram |
|---|---|---|---|---|
| 90 km | 7.82 km/s | 2192 K | about 1390 K | 11.8 MJ/kg |
| 80 km | 7.60 km/s | 3035 K | about 1490 K | 10.0 MJ/kg |
| 70 km | 6.84 km/s | 3587 K | about 1590 K | 6.9 MJ/kg |
| 60 km | 4.13 km/s | 2285 K | about 1690 K | 1.9 MJ/kg |
| 50 km | 1.51 km/s | 964 K | about 1800 K | 0.07 MJ/kg |

**Energy absorbed.** Your earlier rough estimate was that a 10 micrometre droplet released at 80 km absorbs
about 9.9 megajoules per kilogram. This model gives 10.0 megajoules per kilogram (interpolating between the
78 km and 82 km cases at 9.56 and 10.4), so the estimate holds. For comparison, heating liquid aluminium from
908 K to its boiling point and then vaporising all of it takes about 13.1 megajoules per kilogram (2.2 to
heat it, 10.9 to vaporise it). No case absorbs that much (the maximum is 11.8 at 90 km), so the energy
budget alone does not say the droplets vaporise completely. The absorbed energy is between 6 % (at 50 km) and
39 % (at 90 km) of the droplet's initial kinetic energy per kilogram; the rest goes into the air.

**How much would vaporise.** As a rough bracket, not a model, I computed two numbers for each case. The lower
one assumes only the heat the droplet holds above the saturation temperature at its peak is turned into
vapour; the upper one assumes everything absorbed beyond heating to saturation is turned into vapour, with
nothing radiated. For a release at 70 km the vaporised share of the mass lies between 22 % and 56 %; at 80 km
between 17 % and 85 %; at 90 km between 9 % and 100 %; at 60 km between 6 % and 9 %; and zero below 58 km.
The width of these brackets is itself the finding: whether the high-altitude droplets survive as liquid at all
cannot be decided without a vaporisation model, because evaporation competes directly with radiation as the
way the droplet sheds heat. These are columns `vap_frac_lower` and `vap_frac_upper` in the summary.

**Alloying elements.** Zinc boils at 1180 K and magnesium at 1363 K at one atmosphere, far below aluminium.
The droplet passes 1363 K in 162 of 189 cases and 1180 K in 171, so even where aluminium survives, the zinc and
magnesium in it are likely to boil out preferentially. The droplet that comes down is probably not AA7075.

**Sensitivity to the uncertain inputs.** Lowering the emissivity to 0.05 raises the peak temperature of a
70 km release from 3587 K to 3990 K; raising it to 0.4 lowers it to 2809 K; reducing the accommodation to 0.8
lowers it to 3295 K. The number of cases above the one-atmosphere boiling point ranges from 36 (emissivity
0.4) to 117 (emissivity 0.05), but the number above the local saturation temperature hardly moves (144 to 153),
because the droplets released above 58 km overshoot it by hundreds of kelvin under every combination. So the
central conclusion, that vaporisation cannot be left out above about 58 km, does not depend on the
uncertain emissivity or accommodation. The left panel of `droplet_peak_temperature.png` shows this.

**Deceleration.** The peak deceleration ranges from about 56 times Earth's gravity (250 micrometres at 90 km) to
about 200,000 times Earth's gravity (2 micrometres at 62 km), and scales as one over the diameter. It peaks
around 62 km, where the product of air density and the square of the release speed is largest. These loads
matter for whether the droplet breaks up, which is not modelled.

**Losing the release speed.** I report two measures. The first is the time and distance until the speed has
fallen to 1/e (about 37 %) of the release speed. The second, which I call the stopping point, is the moment
the horizontal speed has fallen to 1 % of its release value. The time to the stopping point runs from
14 milliseconds (2 micrometres at 50 km) to 91 seconds (250 micrometres at 90 km), and the downrange distance
at that point from 3.5 metres to 157 km; both scale roughly with diameter and rise steeply with altitude. The
downrange distance at which each droplet eventually reaches the ground is only 0.5–4 % larger than at the
stopping point, because after that the droplet falls almost vertically.

**Reaching the terminal velocity** (added after the first report; `terminal_approach.py`, bottom row of
`figures/4_matrix_heatmaps.png`). The terminal velocity is the local settling speed, at which drag balances
weight. It keeps falling as the droplet sinks into denser air, so I define "reached" as the moment after which
the droplet's speed stays within 1 % of the local terminal velocity for the rest of the fall. The crossing is
interpolated between history samples, which are at most 8 % apart in time. Two behaviours appear. Droplets
released low, or small droplets released high, settle onto their terminal velocity right after decelerating,
at their release altitude: 0.05 seconds and 0.04 m/s for 2 micrometres at 50 km, up to 5 seconds and 8 m/s for
250 micrometres at 50 km. Droplets released higher up, in thinner air, decelerate onto their terminal velocity
within about 20 seconds, but then trail slightly above it: the local terminal velocity keeps falling as the
droplet sinks into denser air, and the droplet's speed lags behind. (An earlier version of this paragraph said
these droplets approach their terminal velocity from below; tracing the histories shows they trail it from
above.) For example, a 10 micrometre droplet released at 88 km is within 1.3 % of its terminal velocity 20
seconds after release, but does not get inside the 1 % band until about 113 seconds, when its terminal velocity
has fallen to about 25 m/s. The crossing into the band comes at about 25 to 30 m/s whatever the size, between 64
and 84 km, 50 to 280 seconds after release. That common speed is set by the 1 % tolerance and the atmosphere's density
scale height H (about 6 km), not by the droplet. A droplet can track a terminal velocity that changes over a
fall time H / v_t only if its response time v_t / g is short enough; the lag is about v_t squared divided by
g H, and setting that to 1 % gives v_t of about 0.1 times the square root of g H, roughly 24 m/s. A 10 %
tolerance would move the crossing to about 77 m/s and much earlier times. So the time in this panel measures
when the droplet starts to track its settling speed closely, and it depends on the tolerance chosen.

**Time to the ground** depends almost only on the diameter, not on the release altitude, because it is
dominated by the slow settling at low altitude: about 650 days for 2 micrometres, 98 days for 6, 40 days for 10,
22 days for 14, 14 days for 18, 2.3 days for 50, 17 hours for 100, 9 hours for 150, and 4 hours for 250.

## 5a. Extension to droplets from 0.01 to 2 micrometres (2026-10-15)

At your request I added ten diameters spaced evenly on a logarithmic scale from 0.01 micrometres (10
nanometres) up to, but not including, 2 micrometres: 0.01, 0.017, 0.0289, 0.049, 0.0833, 0.141, 0.24, 0.408,
0.693 and 1.18 micrometres (the exact log-spaced values rounded to three figures). All 210 new cases, and the
same cases in the three sensitivity sweeps, ran without error, and every one reached the ground. The matrix
is now 19 diameters by 21 altitudes, 399 cases. The original 189 cases are unchanged by the extension.

**Temperature.** The new sizes heat exactly like the old ones. Their peak temperatures agree with the 2 micrometre
droplet to within 0.1 K at every altitude, for the free-molecule scaling reason given in section 2. They are deep
in free-molecule flow: the Knudsen number (air mean free path divided by droplet size) while they absorb heat is
at least 129 and reaches four million. So the vaporisation conclusion is unchanged. Of the 399 cases, 209 exceed
the one-atmosphere boiling point, a further 114 exceed the local boiling point, and 76 (all releases at 50–56 km)
stay below it. None of the new cases carries the drag blend-domain warning or the curvature flag.

**Motion.** Everything that took time or distance shrinks in proportion to the diameter. A 10 nanometre droplet
loses 99 % of its horizontal speed within 0.07 milliseconds at 50 km and 17 milliseconds at 90 km, 1.7 centimetres to 12
metres downrange, under a peak deceleration of up to 39 million times Earth's gravity. It reaches its terminal velocity
straight after decelerating, at its release altitude, at 0.0002 to 0.07 metres per second. None of the new sizes shows
the lagging behaviour of the larger droplets released high up, because their terminal velocities are all far below the
24 metres per second threshold explained above. The time to reach the ground grows steeply as size falls: about
4.3 years for 1.18 micrometres, 46 years for 0.24 micrometres, and about 1,560 years for 0.01 micrometres,
again almost independent of release altitude.

**Where the model is not trustworthy for these sizes.** Four things the model leaves out matter more as the
droplet shrinks, and I have not quantified any of them.

- *Settling is not how sub-micron particles come down.* Below a few tenths of a micrometre, the fall speeds
  computed here (micrometres to millimetres per second) are far smaller than the vertical air motions of the
  atmosphere: turbulent mixing, gravity waves, and the large-scale circulation, which in the mesosphere carries
  air downward over the winter pole. Those motions, and Brownian diffusion for the very smallest sizes, will
  bring them down long before settling does. The times to the ground for these sizes, years to centuries,
  should be read only as the settling contribution, not as predictions.
- *Radiation from a particle smaller than the wavelength it emits.* A droplet at 1,000 to 3,000 K radiates
  mostly at wavelengths of 1 to 3 micrometres. A sphere much smaller than that wavelength does not radiate
  like a piece of bulk surface; for a metal its emission per unit area is much weaker. The model uses the bulk
  emissivity for every size, so it over-predicts the cooling of sub-micron droplets, and their true peak
  temperatures would be higher still. This only strengthens the vaporisation conclusion.
- *Nanometre-scale material properties.* A 10 nanometre droplet holds roughly 30,000 atoms. Its vapour pressure
  is raised by its curvature (the Kelvin effect, about 35 % at 10 nanometres and 1,500 K, using a surface
  tension of 0.86 newtons per metre), and its melting point is lowered by tens of kelvin. Both would make it
  evaporate faster than the bulk property values assumed here.
- *Coarser sampling at the smallest sizes.* The crossing into the terminal velocity band is now interpolated
  between history samples up to 11 % apart in time (it was 8 % for the original sizes), so those times carry a
  few percent of uncertainty.

## 6. Where the results are weaker than they look

**Free-molecule heating is an upper bound for the large droplets at the lower altitudes.** Free-molecule
heat transfer is the most a surface can receive; when the air near the droplet becomes collisional (Knudsen
number, the mean free path of air divided by the droplet size, below about 10), a shock layer forms and the
heat reaching the surface falls. For droplets of 18 micrometres and smaller, all the heating happens at
Knudsen number above 10 in every case except one (18 micrometres at 50 km). For 50 to 250 micrometre droplets
released below about 62 to 70 km, essentially all of the heating happens at Knudsen numbers between 0.6 and
10, so their temperatures are overestimates by an unknown amount. Those cases are mostly the ones that heat
least anyway. The column `frac_E_in_Kn_lt10` gives the share for every case.

**The drag blend-domain warning fires in 80 of 189 cases.** That is every 100, 150 and 250 micrometre case, 8
of the 50 micrometre cases, and the 10 to 18 micrometre cases released at 50–54 km. These droplets pass through
the region of Knudsen number 0.1 to 20 together with speed ratio above 0.05, where your composed drag law is
not calibrated. Where the two drag branches were blended, they disagreed by 15 % to 86 %. Because the
heating per kilogram depends on the ratio of heat transfer to drag, the energy absorbed is uncertain to a
similar degree in those cases, and the times and distances more so. The column `blend_warning` marks them.

**Flat Earth.** For 9 cases (100 micrometres and larger, released at 84–90 km) the stopping point is far enough
downrange that the Earth's curvature would put the ground more than 300 metres lower than the flat-Earth
model assumes, up to 1.9 km for 250 micrometres at 90 km. That would put the droplet in somewhat thinner air
and lengthen its deceleration. These are flagged in `flag_curvature_matters`.

**Single parent trajectory.** All release states come from the one DRAMA descent in `dmf_output.json`. The
sweep was run once; the integration itself is deterministic.

## 7. What it means for the research

The temperature numbers in most of the matrix are not usable as predictions, and this should be stated
plainly in the thesis: a molten aluminium droplet shed above about 58 km would, by this model, heat far
past the point where it boils at the local pressure, and its fate is set by vaporisation, which competes
with radiation as the way it sheds heat. Only the releases at 50–56 km give temperatures that can be read
directly, and there the droplets barely warm (about 56 K above release at 50 km, rising to 565 K at 56 km).

Two results are robust and are useful on their own. First, in free-molecule flow the heating history is
independent of droplet size, so the vaporisation question is a function of release altitude only and does
not need to be resolved separately for each diameter. Second, the droplets that survive will have lost much
of their zinc and magnesium, which bears on what material reaches the ground and what is left in the
mesosphere as vapour.

## 8. What you need to decide

**Whether to add vaporisation to the droplet model.** The options are as follows.

1. *Add a Hertz–Knudsen evaporation sink* (mass loss at a rate set by the vapour pressure at the droplet
   temperature, removing the latent heat of vaporisation), with the droplet's mass and diameter shrinking.
   This is perhaps one day of work, uses the same vapour-pressure correlation already in the script, and
   would turn the 153 unphysical cases into predictions of surviving mass and final size. It needs an
   evaporation coefficient, another uncertain number of order 1 for clean liquid metals.
2. *Clamp the temperature at the local saturation temperature* and put every excess joule into vapour. This
   is simpler and gives an upper bound on vaporised mass, but it does not distinguish the zinc and magnesium.
3. *Leave the model as it is and report only the 50–56 km releases* as quantitative, with the rest as a
   statement that vaporisation dominates.

My recommendation is option 1, because the vaporised fraction brackets above are too wide (for example
9 % to 100 % at 90 km) to support any conclusion about what survives, and the evaporation sink is a small
extension to the existing energy equation. A selective-evaporation version for zinc and magnesium could follow
if the aluminium result shows the droplets survive.

**Whether the wake matters.** If droplets are shed into the near wake rather than the free stream, the release
speed relative to the air could be much lower than the parent's speed, and the heating correspondingly
lower. I would want your view on how far behind the parent the melt spray is released before treating the
free-stream assumption as settled.
