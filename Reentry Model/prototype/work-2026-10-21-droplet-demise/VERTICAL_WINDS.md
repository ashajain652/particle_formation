# Vertical winds and the droplets: what they change, how large they are, how to add them

> **Moved 2026-10-26.** The code of this folder now lives in `reentry_model/droplet_reentry/` (`droplet_demise.py` became `model.py`; the sweeps, plots and checks are in its `analysis/` and `tests/`), and its outputs are written to `output/droplet_reentry_output/`. This file is kept as the record of the work; script names and `runs/` paths below are as they were.

2026-10-22. Research note, then implemented as bounding cases (results at the end).

## What an updraft does to a droplet

At the terminal state, drag balances weight, so the droplet's speed relative to the air is fixed by its weight
and size alone. An updraft of speed w does not change that relative speed, and it adds no drag. What it changes
is the droplet's speed over the ground: the droplet sinks at (terminal speed minus w). If w exceeds the terminal
speed, the droplet rises with the air. Updrafts therefore slow the descent and lengthen the time aloft. They do
this by carrying the air, and the droplet with it, upward, not by increasing drag.

During the release and vaporisation phase the droplet moves at 1.5 to 7.8 km/s through the air. A vertical wind
of even 2.5 m/s changes that airspeed by 0.03 to 0.17 %, so mass loss and peak temperature are unaffected. The
hunch that updrafts only matter near the terminal state is right.

The current demise model stops at the terminal state, close to the release altitude. It therefore never reaches
the phase where updrafts matter. To see their effect, the run has to continue through the settling phase, for
example to a fixed lower altitude or to the ground.

## How large the vertical winds are (60-90 km)

- **Mean (residual) circulation:** a slow, persistent flow set up by breaking gravity waves. It rises over the
  summer pole and sinks over the winter pole at a few cm/s; estimates for the southern summer polar mesosphere
  peak at 2-6 cm/s (Vincent 2019). It is weaker at mid and low latitudes. This persistent upwelling keeps the
  ice particles of noctilucent clouds aloft against sedimentation (NASA NTRS 19860018257).
- **Gravity waves:** oscillating vertical winds of about 1 m/s to a few m/s near 70-90 km, with periods of about
  10-40 minutes and roughly zero mean. Measured by the Arecibo incoherent-scatter radar (about 1 m/s accuracy over
  68-88 km, amplitudes of a few m/s) and the EISCAT VHF radar (amplitudes up to about 2.5 m/s).
- **Rare extremes:** localised events reaching tens of m/s (up to 60 m/s reported by MAARSY). These are not
  representative.
- **Below about 50 km:** the stratospheric mean upwelling is much weaker (fractions of a mm/s), and the
  troposphere has convective updrafts of m/s to tens of m/s. Neither matters for the release phase here.

## Measured against the droplets' terminal speeds (Singh-drag sweeps)

| comparison | 1-200 um droplets (315 cases) | 1-500 nm droplets (315 cases) |
|---|---|---|
| terminal speed at the final state | 0.018 to 26 m/s | 0.00002 to 2.8 m/s |
| below a 5 cm/s mean upwelling (would be lofted) | 10 cases (3 %): 1-2 um, high releases | 249 cases (79 %) |
| below 1 m/s (gravity-wave scale) | 98 cases (31 %) | 309 cases (98 %) |
| below 2.5 m/s (gravity-wave peak) | 140 cases (44 %) | 314 cases (100 %) |

The mean upwelling decides whether a droplet can be held up indefinitely. It does so for every nanometre droplet
up to about 6 nm at every altitude, for roughly half of the 100-500 nm droplets, and for micron droplets only at 1-2
micrometres. Gravity-wave winds are larger than most of the terminal speeds but average to zero, so they mainly
spread out the descent times rather than shift their mean. For droplets above about 20 micrometres (terminal
speeds of metres per second up to 26 m/s), vertical winds are a modest correction.

## How the winds could be computed

The atmosphere model already carries a wind vector (east, north, up), and the trajectory code already subtracts it
from the droplet's velocity, so any of these plugs in without changing the dynamics.

1. **Prescribed mean upwelling (recommended first step).** A vertical-wind profile w(altitude), for example a
   summer-polar upwelling of +5 cm/s, a winter downwelling of -5 cm/s, and zero, run as bounding cases. This is
   cheap and deterministic, and it answers which droplets are lofted and how much the descent time changes.
   Uncertainty: only order-of-magnitude values; latitude and season enter only through the choice of case.
2. **Stochastic gravity-wave winds.** A sum of random-phase waves with about 1 m/s amplitude and 10-40 minute
   periods, run as a Monte Carlo over many realisations, gives the spread of descent times. It needs many runs
   per case and a choice of wave spectrum.
3. **A model climatology.** For example the mean vertical velocity from the WACCM whole-atmosphere model (an
   SD-WACCM residual-circulation dataset is on Zenodo, about 300 MB; WACCM-X monthly means are on CEDA). This gives
   realistic latitude and season dependence, but needs a data download and the residual velocity derived from the
   model fields. Empirical wind models such as HWM14 give horizontal winds only.

## Sources

- Vincent (2019), Trends and Variability in Vertical Winds in the Southern Hemisphere Summer Polar Mesosphere and
  Lower Thermosphere, JGR Atmospheres: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JD030735
- Chau et al. (2019), Can VHF radars at polar latitudes measure mean vertical winds in the presence of PMSE?, ACP
  19, 4485: https://acp.copernicus.org/articles/19/4485/2019/
- Connection between mesopause temperature, circulation and noctilucent clouds (1986), NASA NTRS:
  https://ntrs.nasa.gov/citations/19860018257
- Vertical velocities associated with gravity waves measured in the mesosphere (Bath):
  https://researchportal.bath.ac.uk/en/publications/vertical-velocities-associated-with-gravity-waves-measured-in-the/
- Arecibo incoherent-scatter vertical winds (COSPAR 2006 abstract):
  https://meetings.copernicus.org/www.cosis.net/abstracts/COSPAR2006/00185/COSPAR2006-A-00185.pdf
- MAARSY varicose-mode event, ANGEO 42, 29 (2024): https://angeo.copernicus.org/articles/42/29/2024/
- Sato, Kinoshita and Okamoto (2013), residual circulation drivers:
  https://www-aos.eps.s.u-tokyo.ac.jp/~sato-lab/en/publications/publications/PDFs/SatoKinoshitaOkamoto-JAS2013.pdf
- Meteoric smoke and vertical wind (ACP 6, 4415, 2006): https://acp.copernicus.org/articles/6/4415/2006/
- SD-WACCM MLT residual circulations (Zenodo): https://zenodo.org/record/5804538
- WACCM-X monthly means (CEDA): https://catalogue.ceda.ac.uk/uuid/45283390b97c4a27861d74b3d915b0bd/

## Implemented: bounding cases and results (2026-10-22)

`droplet_demise.MeanVerticalWind` carries a mean vertical wind through the trajectory. `settle` then follows each
droplet from its final state to the ground: it sinks at its terminal speed minus the local updraft, at the air
temperature and with its final mass. It also reports the time to the tropopause (15 km). The profile has:
- the mesospheric value from 70 to 90 km, tapered to zero at 100 km;
- a linear change to the stratospheric value at 50 km;
- the stratospheric value from 15 to 50 km;
- no mean wind below 13 km.

Three cases bound the residence time:
- **upward:** +5 cm/s in the mesosphere and +0.3 mm/s in the stratosphere;
- **none;**
- **downward:** -5 cm/s and -0.3 mm/s.

No single place has either the upward or the downward combination; they are limits. Runs are under
`runs/demise_sweep[_nano]_w<meso>cms_s<strat>mms/`, figures 8 and 9. `verify_droplet_demise.py` checks:
- the settling quadrature against a direct integration, agreeing to 5 digits;
- the time to the ground against the 2026-10-14 sweep, 39.88 against 39.81 days;
- a suspended droplet hovering where its terminal speed equals the updraft, to 0.1 %.

**Mass loss and peak temperature do not change at all**, to 0.001 percentage points, as expected.

**Upward case.**
- Every nanometre droplet (1 to 500 nm) is held aloft indefinitely, at 26 to 97 km.
- Micron droplets of 1.46 micrometres and smaller are held in the stratosphere, at 16 to 22.5 km.
- Droplets of 2.13 micrometres and larger still reach the ground, up to 45 % later.

**Downward case.**
- Nothing is held aloft, and times to the ground fall by up to 17 %.
- The large change is in the time to reach the tropopause for small droplets. A 1 nm droplet released at 80 km
  takes 2,300 years with no wind and 3.7 years with the downwelling. A 1 micrometre droplet takes 2.0 years with
  no wind and 1.0 year with the downwelling.

**Above about 10 micrometres** the winds change the times by under 3 %.

**Limits of this result.**
- "Held aloft" is relative to a fixed column of air. The real residual circulation also moves air horizontally,
  from the summer-pole upwelling towards the winter-pole downwelling. So a held droplet is carried sideways into a
  region that brings it down; its residence time is long, but not infinite.
- Gravity-wave oscillations and horizontal transport are not modelled.
- The times to the ground for droplets smaller than a few micrometres are dominated by slow settling through the
  troposphere. In reality rain removes such particles much faster there, so the time to the tropopause is the
  more meaningful result for them.
