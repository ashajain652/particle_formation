# Wind shear in the wake column, and mesospheric winds to drive it (milestone, 2026-11-05)

This note records a goal Asha set on 2026-11-05: to replace the wake column's round plume, which grows equally in every
direction, with a sheared plume, and to drive it with real winds for 50 to 90 km. Nothing here is implemented. It
follows the review `Literature Review/Alumina Nanoparticles/reports/Contrail and rocket plume spreading.md` and needs
a design spec and plan before building.

## Why

The model spreads the plume with one diffusivity, r^2 = r0^2 + 4 integral(D dt), D = D_mol + min(K_eddy, K_R(r)), and has
no wind. Contrail work and rocket plume measurements show a line source in stratified, sheared air becomes a thin, wide,
tilted ribbon. Shear alone keeps the cross-section's area; it adds area by stretching what vertical diffusion has mixed.

The dilution error can go either way. At 70 km, against the current model, the review's sheared area is 0.09 to 9 times
at 1 hour and 0.4 to 108 times at 12 hours, for S of 0.005 to 0.01 per second and D_V of 0.1 to 5 m^2/s. D_V and how
long shear keeps one direction decide the sign. These are calculations with assumed inputs, not model results.

## Goals

1. Track the plume's cross-section as a Gaussian covariance, with shear from a wind data set.
2. Keep the single-volume box model and its chemistry and particle terms unchanged.
3. Reduce exactly to the current law at zero shear, so results compare directly.
4. Report how much the 12-hour gas chemistry and condensed fractions change.

Not in this milestone: advection of the plume to other air (phase two of the design, though the same winds would serve
it), a 3D plume model, segregation of reactants inside a contorted plume, or day-night changes in the ambient column.

## The model change

Following Konopka (1995) and Schumann (2012, CoCiP Eqs. 29-32), for y across the trail in the horizontal, z vertical:

- d sigma_yy/dt = 2 D_H + 2 S sigma_yz
- d sigma_yz/dt = S sigma_zz
- d sigma_zz/dt = 2 D_V

A = 2 pi sqrt(sigma_yy sigma_zz - sigma_yz^2) replaces pi r^2, with V = A L, or the existing erf form with
r = sqrt(A / pi) for short trails. Here S(t) is the shear component normal to the trail axis, and sigma is the per-axis
variance (r^2 = 2 sigma in the current convention).

D_V = D_mol + min(K_eddy, K_R(r_z)) and D_H = D_mol + min(f_H K_eddy, K_R(r_y)), as the review proposes. With S = 0 and
f_H = 1 this is the current law. Settling (2 v_s / pi r) needs a new length, the vertical size.

Shear along the trail axis stretches the trail lengthwise instead, which conserves volume. A trail's own altitude span
also matters: a 20 km trail descending 5 degrees spans 1.7 km. Neither has a source; the design must decide them.

## Wind data

The need is S(t) at the trail's altitude and latitude over 12 hours, including tides, with statistics by season.

| Source | Altitudes | Use | Status |
|---|---|---|---|
| WACCM (our CESM2 files) | full 45-105 km | Baseline: same model as the ambient columns | Tested below; time-resolved output not yet found |
| SD-WACCM or WACCM-X output | full | Time-varying winds and tides | From memory; availability to check |
| Reanalyses, ERA5 and MERRA-2 | to about 80 km | Winds at 50-80 km | From memory; no constraint from data near the top, to check |
| HWM14 empirical winds | full | Climatological tides, quick cross-check | From memory; no small-scale shear |
| UA-ICON (Wing et al. 2026) | to 150 km | 20 km, 10 minute winds | One event; output availability unknown |
| SIMONe-Germany + HYPER (Wing et al. 2026) | 80-100 km, 53 N | Measured shear, 1 km, 15 minute | DOI 10.22000/gdp71mh7gm8ufks6 says figure data; contents unchecked |
| Colorado State sodium lidar (Yue et al. 2010) | 80-105 km, 40.6 N | Check at our latitude | Abstract level; large shear defined as above 0.05 per second |

Measured once on 2026-11-05 from the WACCM starting file (2001-01-01 00 UT, 40.7 N, all longitudes), by differencing the
wind vector between adjacent levels: median shear 0.003 to 0.006 per second and 90th percentile 0.007 to 0.016 per
second between 45 and 105 km. That supports the review's judged 0.005 to 0.01 per second. It is one instant, and the
levels are 1.8 to 3.5 km apart, so it misses thinner layers and is a lower limit on true shear. Its file has winds at
one time only, so it cannot show how long shear persists, which is the main unknown.

Proposed route: use a time-resolved WACCM run for S(t) at 50 to 90 km, since it matches the ambient air and its tides,
check it against the lidar at 40.6 N and the radar near 53 N above 80 km, and carry unresolved small-scale shear as a
sensitivity. Fall back to a prescribed tidal shear, S = S0 cos(omega t), if no suitable output exists.

## Inputs still uncertain

- D_V is the largest unknown: the review's cases are Megner's K_eddy and one fifth of it.
- D_H is unmeasured at these altitudes, though the late-time area hardly depends on it.
- The tidal D_eff = D_V S0^2 / (2 omega^2) is our own derivation, checked numerically but not against Young et al. (1982).
- Rocket plumes grew linearly with no link to shear (Beiting 1999, 18-30 km), against contrail findings. Unresolved.

## How we will know it worked

- Zero shear reproduces the current results to solver tolerance.
- Covariance terms match the analytic constant-shear solution, including sigma_yy ~ (2/3) S^2 D_V t^3.
- Aluminium stays conserved, and the sensitivity set gains shear and D_V cases.
- The layer thickness compares sensibly with the Falcon 9 lithium layer (2.3 km at 20 hours near 96 km), a weak check.

## Decisions for Asha

1. Default D_V: Megner capped, or one fifth of it. The review recommends both, Megner as default.
2. Whether time-varying model winds are worth obtaining, or prescribed tidal shear is enough first.
3. Fixed trail orientation against the wind, or a range of orientations.
