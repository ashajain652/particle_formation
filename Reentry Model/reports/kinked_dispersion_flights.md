# Kinked-profile dispersion: the four comparison flights

Written 2026-11-12 for Asha. This is Task 9, Step 2 of the kinked-profile dispersion plan
(`docs/superpowers/plans/2026-11-11-kinked-profile-dispersion.md`). The flights' files are in
`reports/kinked_dispersion_flights_macstudio/`. They came from the Mac Studio on branch
`kinked-dispersion-flights-macstudio`, at commit c3247d3.

## What was asked, and whether it is done

The plan asked for two pairs of flights. Each pair was flown once with Girin's theory as it stood before this branch
(`--dispersion table`) and once with the branch's new default (`--dispersion profile`). The profile path reads each
skin column's own melt speed profile, finds where the shear drops sharply (the kink), and tests whether the slurry
under the kink is stiff enough to act as a rigid floor. The report compares the two settings on four things:

- the share of columns whose wall is floor-like;
- where the walls sit against the 906.4 K isotherm, the low end of the rheology's viscosity bridge;
- the mass each spraying mode releases;
- the droplet sizes.

It also flags anything that argues against making `profile` the default.

All four flights ended normally, at their set end times, with exit status 0, so the comparison is done.

- Flights (a) and (b) are the 50 mm body. Each resumed from the saved state at 178.5 seconds and ran to 179.0 seconds.
- Flights (c) and (d) are the 100 mm sphere, flown from the start of the trajectory (77.5 km altitude, 7.5 km/s) to
  30 seconds. Each saved its state at 25 seconds.

None of the four was restarted. The 25-second saved states (1.39 GB each) stay on the Mac Studio, with their checksums
in `large_files.txt`.

## What was found

On the 50 mm body nothing sprays under either theory. The profile path leaves the body exactly as the table path does.

On the 100 mm sphere the two theories disagree completely:

- **Girin's theory (table)** releases nothing in 30 seconds.
- **The profile path** releases 2.52 grams as about 5,300 droplets. Their median radius is about a third of a millimetre.

Every gram of that comes from one mechanism. The profile path finds a floor-like wall about 66 micrometres below the
surface, inside a flowing layer about 0.85 mm deep. It then sprays the film above that wall with the gas-driven
thin-film mode, which has no critical Weber number. Girin's theory instead treats the whole 0.85 mm as a deep layer,
and on a deep layer the surface is stable at this flight's Weber numbers.

Two things make the 2.5 grams a provisional number rather than a result:

- **The release depends on the near-liquidus viscosity law.** It passes the viscous-damping gate only narrowly. The
  wall that sets the film's depth sits inside the 1.6-kelvin viscosity bridge just below the liquidus, and that bridge
  is the part of the rheology the near-liquidus viscosity branch replaces.
- **The release was still growing when the flight stopped.** It was rising quickly at 30 seconds.

Nothing found here argues against `profile` as the default. Where the profile path sprays, it does so by Girin's own
regime test applied to the depth that can actually move.

## Evidence

### The 50 mm body, flights (a) and (b)

Neither flight released any metal. All 97 history columns the two flights share are identical, apart from the
wall-clock timers.

Profile's row at 179.0 seconds:

- **Floor-like share:** 0.82. Four in five of the columns the profile path resolved have a floor-like wall. The other
  172 columns are fluid-like: they have no floor-like kink above the rigid front.
- **Wall depth against the isotherm:** the median wall lies 0.064 mm below the surface. The median 906.4 K isotherm
  lies at 0.058 mm, and the median temperature at the wall is 906.3 K. So on this body the walls sit at the isotherm,
  within a few hundredths of a millimetre.
- **Floor Reynolds number:** median 0.038. This number measures whether the slurry under the kink can move fast enough
  to take part in a wave. Below 1 the slurry acts as a floor.
- **Surface Weber number:** median 0.00044. Surface tension therefore holds Girin's gradient mode stable everywhere.
- **Cost:** the timing check put the profile path's added cost at 1.2 seconds per macro step. That is 0.035 % of a
  median step of 3,480 seconds of wall time, against a budget of 5 %.

### The 100 mm sphere, flights (c) and (d)

**The two flights are identical until the first droplet leaves.** Metal first starts to flow at 22.5 seconds. Every
shared history column stays identical between the two flights through 24.6 seconds. The profile flight's first droplet
leaves at 24.6125 seconds, at 74.1 km altitude and 7.38 km/s. After that the two bodies differ only through the metal
that was sprayed off.

**How much left, and when.**

- **Total:** the profile flight released 2.524 grams, which is 0.17 % of the sphere's 1.472 kg. It left as 5,337
  droplets, from 3,950 patch-steps (one record per surface patch per step that released metal).
- **Table flight:** released nothing. No patch on it ever released, and its spraying area is zero on every row.
- **Where:** every release came from 5.8 to 13.3 degrees off the stagnation point (the 10th to 90th percentile).

The table below gives the profile flight's release in each one-second window.

| Window (s)         | 24–25 | 25–26 | 26–27  | 27–28 | 28–29 | 29–30 |
|--------------------|-------|-------|--------|-------|-------|-------|
| Released (grams)   | 0.017 | 0.136 | 0.0007 | 0.359 | 0.811 | 1.201 |

The release more than doubled from the 27–28 window to the 28–29 window, then grew by half again in the last second,
which gave almost half the total.
The 2.5 grams is therefore what had left by 30 seconds, not the flight's total. The near-pause between 26 and 27
seconds comes from the heating, not the spraying: the liquid depth drops to almost nothing at 26 seconds on both
flights.

**Why Girin's theory released nothing.** The patches that sprayed on the profile flight hold a flowing layer 0.85 mm
deep (median). That is about nine times Girin's melt boundary-layer thickness, delta_m, the depth over which the gas's
shear is felt in the melt, which is 92 micrometres here. Girin's regime test therefore calls the layer deep, and the
table path sends it to his thick branch.

The thick branch sprays only where the surface Weber number exceeds the critical value of 3.08. The surface Weber
number is the melt's surface inertia measured against surface tension: the liquid density, times Girin's surface speed
squared, times delta_m, divided by the surface tension. On these patches it is 0.13 (median; 0.08 to 0.40 from the
10th to the 90th percentile), which is more than seven times too small. The surface is stable and nothing leaves.

Task 7's one-step check on the 100 mm test sphere does not carry over to this flight. In that check, the table path
released 1.54 grams in a quarter of a second. On this flight the table path never found a patch above the critical
Weber number before 30 seconds.

**Why the profile path released.** The profile path looks inside that 0.85 mm. It finds a kink 66 micrometres below the
surface (62 to 71 micrometres from the 10th to the 90th percentile), where the shear drops at least five-fold. Below the
kink the floor Reynolds number is 0.17 (0.10 to 0.27), far under the threshold of 1, so the column is floor-like and its
wall is at the kink.

A 66-micrometre film over a floor is shallower than delta_m: the wall lies at 0.73 of delta_m. That is Girin's other
case, the thin-film mode of Girin and Kopyt (1994), which is driven by the gas. That mode has no critical Weber number.
Two tests hold it back:

- **The viscous-film gate:** the wave's viscous damping over one growth time must stay below 1.
- **The wave-fit test:** one whole wavelength must fit in the molten region.

Girin's gradient mode never sprayed: its surface Weber number stayed between 0.0001 and 0.0003, so surface tension held
it, as on the 50 mm body. No patch on either flight was taken by the Rayleigh–Taylor mode or by any other branch. Every
release record is the gas-driven mode, and the history's gas-driven column sums to the whole 2.524 grams.

**Floor-like share.**

- At the first spraying step every resolved column was floor-like (share 1.000, no fluid-like columns).
- At 30 seconds the share was 0.979, with 13 fluid-like columns.
- In the release records, 3,949 of 3,950 patch-steps had a floor-like wall. The single fluid-like one released 3.5
  micrograms.

**Where the walls land against the 906.4 K isotherm.** This comes from two sources, which can differ because the two
history medians may be taken over different columns.

- **From the history (medians over all resolved columns):**
  - At the first spraying step, the wall lies 0.056 mm deep and the isotherm 0.088 mm deep, and the wall's temperature
    is 906.97 K.
  - At 30 seconds, the wall lies 0.065 mm deep and the isotherm 0.114 mm deep, at 907.10 K.
- **From the release records (each patch that sprayed):**
  - In 3,879 of the 3,949 floor-like records the wall is the shallowest kink itself.
  - In the other 70 records the shallowest kink was not floor-like, and the wall is a deeper kink, 1.3 to 4.7 times as
    deep. Those 70 released 42 micrograms between them.
  - The kink's temperature is 907.03 K (906.88 to 907.14 K).
  - The one fluid-like record has its wall at 852 micrometres, 15.7 times its kink depth of 54 micrometres, with the
    kink at 907.23 K.

So on the 100 mm sphere the walls sit 0.03 to 0.05 mm shallower than the 906.4 K isotherm, about 0.6 K warmer than it.
That places them inside the bridge between 906.4 K and the 908.0 K liquidus, where the model's viscosity is a log-linear
interpolation between the slurry law and the liquid. On the 50 mm body they sat at the isotherm. The wall is not pinned
to the isotherm: it follows the speed profile, and on the 100 mm sphere the profile kinks inside the bridge.

**The viscous-film gate passes narrowly.** The gas-driven mode sprays only where the damping-to-growth ratio, gamma
tau*, is below 1. On the patches that sprayed, its median is 0.976 (0.93 to 0.996). Patches above 0.9 gave 98 % of the
released mass, and patches above 0.95 gave 77 %.

The most likely reading is a threshold: a patch sprays as soon as warming lowers its film's effective viscosity enough
to take the ratio under 1. The release rate is then set by how fast the top film warms through that threshold. This is
inferred from the distribution, not traced patch by patch. If it holds, the timing and the mass both depend on the
effective viscosity of the top 66 micrometres, and near the liquidus that viscosity comes from the bridge.

**Droplet sizes.**

- **Number-weighted median radius:** 330 micrometres.
- **Mass-weighted median radius:** 409 micrometres.
- **Range:** 239 to 776 micrometres.
- **Composition:** the droplets are a slurry, 2.1 to 3.5 % solid (median 2.6 %).

Their radius is a quarter of Girin and Kopyt's wavelength. That wavelength is set by the gas (the edge Mach number and
momentum flux), not by the profile. So the profile path decides whether a patch sprays and how much leaves, but not how
large the droplets are.

Their Weber number in the freestream is 1.8 to 3.2, well under the breakup threshold of 12, so none is expected to break
up further when it is released.

**What the release did to the body.** At 30 seconds:

- **Mass:** the profile sphere weighs 1.46949 kg against the table's 1.47202 kg.
- **Flowing depth:** its mean flowing depth is 0.053 mm against 0.068 mm, 22 % less, because the sprayed metal has
  gone.
- **Surface fold:** its largest fold is 22.5 degrees against 18.2. Both are far below the folds that stopped earlier
  skin flights.
- **Peak surface temperature:** the same on both, 908.5 K.

### Cost

The 100 mm table flight took 13,453 seconds of wall time (3.7 hours) and the profile flight 13,861 seconds (3.9 hours),
3.0 % longer. Their peak memory footprints were 2.30 and 2.34 GB.

The skins and the 3D solve took the same time on both:

- 9,123 against 9,103 seconds in the skins;
- 1,457 against 1,461 seconds in the 3D solve.

Most of the extra 400 seconds is outside those timers, which fits the profile path's planning cost. All four flights
ran side by side and each ran once, though, so the 3 % includes contention between them.

### Anything unexpected

- No flight failed or was restarted.
- The Mac Studio's software matches the laptop's except for Python, which is 3.12.15 there against 3.12.14 here. The
  numerical packages are the same versions: numpy 2.5.3, scipy 1.18.1, scikit-fem 12.0.2, gmsh 4.15.2, pyamg 5.3.0.
- The 100 mm flights' JSON files record commit 7f9c8e6 as their code. That is the commit holding the 50 mm results,
  made in the same checkout while the 100 mm flights ran. It adds only result files to 45ea185, so the code the
  flights ran is 45ea185's.

## What this means for the research

**The new default changes nothing where nothing should spray.** On the 50 mm body it changes nothing at all, and it
costs a fraction of a percent of a step.

**Where melt flows near the liquidus, it turns a stable surface into a spraying one.** On the 100 mm sphere it
recognises that only the top 66 micrometres of a 0.85 mm flowing layer can move with a wave. That turns Girin's
deep-layer verdict ("stable") into his thin-film verdict ("sprays"). This is the change the design set out to make, and
on this flight it decides whether any metal leaves before 30 seconds.

**The amount that leaves is not yet trustworthy as a number.** It is set by a viscous gate that passes narrowly, and by
a wall that sits inside the interpolated viscosity bridge, so it inherits the bridge's uncertainty directly. The
near-liquidus viscosity law (branch `near-liquidus-viscosity`) replaces exactly that part of the rheology. Its stiffer
or softer near-liquidus slurry could move the wall and the gate either way. These flights therefore test the mechanism
and its cost. They do not give the 100 mm sphere's spray mass.

Two of the open decisions are not tested by these flights at all:

- **The smooth floor factor.** Task 10's smooth floor factor would change only the gradient mode's growth, and the
  gradient mode released nothing on either body. The floor Reynolds numbers here (0.04 to 0.27) are so far below 1 that
  the walls are floor-like under either the switch or the smooth factor.
- **The fluid-like column rule.** It released 3.5 micrograms in all.

## What you need to decide

**1. Merge `kinked-dispersion` with `profile` as the default?** The options are:

- **(a) Merge now.** The near-liquidus viscosity branch, merging second, then makes its one-line change: the skin
  body's layer columns read the stored shear rates.
- **(b) Hold the merge** until the near-liquidus law is built, and re-fly the 100 mm pair on both together first.
- **(c) Merge with `table` as the default.**

I recommend (a). The default behaves as designed, and the 50 mm body shows it is harmless where nothing sprays. Holding
the merge leaves main without the profile path for the length of another build. Option (c) would hide the mechanism this
branch exists to add.

**2. Re-fly the 100 mm pair past 30 seconds, after the near-liquidus law is merged?** This is the number that matters for
spray mass, and it needs both pieces in place. I recommend it, on the Mac Studio, from the start of the trajectory. The
25-second saved states cannot be reused, because metal starts to flow at 22.5 seconds, so the new viscosity law changes
the flight before 25 seconds.

**3. The constants and the smooth factor.** I recommend keeping a floor Reynolds number threshold of 1.0 and a kink
drop of 0.2, and leaving the smooth factor for a later spec amendment. These flights give no reason to change either.
