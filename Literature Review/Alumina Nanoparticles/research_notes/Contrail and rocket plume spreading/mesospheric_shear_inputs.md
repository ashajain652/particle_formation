# Mesospheric inputs for a contrail-type anisotropic plume model, 50-90 km

Scope and honesty note. About 12 web searches and fetches were done. Most journal papers were behind paywalls or returned only abstracts, so published numbers below are abstract-level. Items marked CALC are my own arithmetic (shown). Items marked MEMORY are from background knowledge and were NOT verified in this session; do not cite them without checking. Calculations use the US Standard Atmosphere 1976 shape (global mean, not 40.7 N seasonal) and the sheared-line formula sigma_s^2 = (2/3) D_v S^2 t^3 from the task statement (standard result for a line in constant shear with vertical diffusion; not re-derived from a source here).

## 1. Measured vertical shear of horizontal wind, 50-90 km

### Takeaway
Direct mid-latitude shear statistics at 50-90 km were not found; the published data sets (chemical tracers, Colorado State sodium lidar) sit mostly at 80-115 km. Shears of 0.05 per second (50 m/s per km) occur and mark the Richardson-number-1/4 edge, but they are the tail, not typical. For model purposes, a typical S of 0.005-0.01 per second with 0.03 as a rare strong case is a defensible but only partly sourced range.

### Cited Findings
- Larsen (2002), a review of more than 400 chemical-tracer wind profiles since 1958 at many latitudes: a wind maximum at 100-110 km at mid and low latitudes, with large shears on its lower flank "often large enough to be either unstable or approaching instability"; features not tied to site or season — [Larsen 2002 via search summaries](https://ntrs.nasa.gov/citations/20030106249) (abstract-level only; DOI 10.1029/2001JA000218).
- Larsen and Fesen (2009, Annales Geophysicae 27, 2277): the Horizontal Wind Model underestimates peak winds and shears in the lower thermosphere (95-115 km) though mean winds agree reasonably; shear statistics there serve as a model test — [Annales Geophysicae](https://angeo.copernicus.org/articles/27/2277/2009/). I could not extract numbers (PDF unreadable). This paper is about 95-115 km, so it is largely outside the 50-90 km range.
- Yue, Jia, She, Liu (2010), Colorado State University (CSU) sodium wind-temperature lidar, mid-latitude (Fort Collins, about 40.6 N, MEMORY), 2002-2005, about 1,600 nocturnal hours, 80-105 km: "large" shear defined as above 50 m/s per km; observed shears fall below the Ri = 1/4 line when plotted against N^2; about 60 percent correlation between large shears and tidal components (24, 12, 8, 6 h) — [UCAR record](https://ckanprod.data-commons.k8s.ucar.edu/dataset/large-wind-shears-and-stabilities-in-the-mesopause-region-observed-by-na-lidar-at-midlatitude).
- Same group: mesopause is most stable in summer (low mesopause, large N, small S); monthly mean N^2/S^2 peaks near 1.06 in mid-summer — [search summary of CSU paper](https://ckanprod.data-commons.k8s.ucar.edu/harvest/object/dfa4e233-976d-41a5-8af1-215f6cfeaa73/html). Inference: monthly mean Ri is about 1 at best, so Ri below 1/4 is reached in the tail of the distribution.
- Andes Lidar Observatory, Chile (about 30 S, 85-100 km), over 2,000 h: dynamic instability probability, with either weak stability or large shear alone able to contribute about 15 percent — [search summary](https://par.nsf.gov/servlets/purl/10633335). This is subtropical southern latitude, not mid-latitude north.
- Tromso sodium lidar (about 70 N): convective instability (N^2 below 0) in 9 percent of cases on average, range 1-24 percent; least instability near 90 km — [UiT record](https://munin.uit.no/handle/10037/28956). High latitude, so only a rough guide.
- Large shears are sustained where static stability is high (lower thermosphere); tides plus gravity waves create thin instability layers — same CSU summaries above.

### Inferences
- CALC: Ri = N^2/S^2 using N about 0.0163-0.022 per second (section 2). For S = 0.005, 0.01, 0.03, 0.05 per second, Ri at 60 km = 0.0163^2/S^2 = 2.66e-4/S^2 gives 10.6, 2.7, 0.30, 0.11. At 90 km (N^2 = 4.84e-4): 19.4, 4.8, 0.54, 0.19. So S = 0.03 per second is already near Ri of 0.3-0.5 and S of about 0.04-0.05 per second is the instability threshold (S = 2N). Shear above roughly 0.03 cannot persist long; this bounds S for a 12-hour model.
- Shear is vertical-scale dependent. Lidar and tracers resolve about 1-2 km; the shear felt by a trail of vertical extent h is U difference over h, so coarse and fine shears differ. Tide with amplitude about 40 m/s and vertical wavelength about 30 km gives peak shear 2 pi x 40 / 30,000 = 0.0084 per second (CALC; amplitude and wavelength are MEMORY-level assumptions, label as illustrative).
- Time scales: tidal shear oscillates with 12 or 24 h periods, gravity-wave shear with minutes to hours (periods between buoyancy period of about 5-6 min and inertial period of about 18 h at 40.7 N: 2 pi /(2 Omega sin lat) = 86164/(2 x 0.652) = 18.4 h, CALC).

### Gaps
- No sourced shear statistics (rms, percentile, vertical scale) specifically at 50-80 km mid-latitude. Candidate sources not retrieved: meteor radar (80-100 km only), MLS (no wind), rocket soundings, MF/MST radar, Rayleigh lidar winds (Doppler Rayleigh lidar to about 60-70 km). Fritts and Alexander (2003) review, Fritts and Williams-type lidar papers, and Liu's Na lidar papers were not opened.
- Fraction of time with Ri below 1/4 at 50-90 km mid-latitude: not found. The only quantified values are high-latitude or subtropical (above).

## 2. Buoyancy frequency N by altitude, 60-90 km

### Takeaway
N is about 0.016-0.022 per second (buoyancy period 4.8-6.4 min), rising with height. These are CALC from US Standard Atmosphere 1976, not measurements; real N varies strongly (summer/winter, waves), with static instability (N^2 below 0) appearing in a few percent of lidar profiles.

### Cited Findings
- Lidar studies derive N^2 from temperature profiles; Tromso static instability fraction about 9 percent (1-24 percent) — [UiT](https://munin.uit.no/handle/10037/28956). No sourced N table was found.

### Inferences
CALC. Formula N^2 = (g/T)(dT/dz + g/cp), with cp = 1004.7 J per kg per K, g(z) = 9.80665 x (6356.766/(6356.766 + z))^2 with z in km. Temperatures from the 1976 standard lapse rates (isothermal 270.65 K to 51 km, minus 2.8 K/km to 71 km, minus 2.0 K/km to 84.85 km, 186.87 K above; geopotential height treated as geometric, an error of order 1 K, label approximate).

| Altitude (km) | g (m/s^2) | T (K) | dT/dz (K/km) | N^2 (1/s^2) | N (1/s) | Period 2 pi/N (min) |
|---|---|---|---|---|---|---|
| 60 | 9.62 | 245.5 | -2.8 | 2.66e-4 | 0.0163 | 6.4 |
| 70 | 9.59 | 217.5 | -2.8 | 2.98e-4 | 0.0173 | 6.1 |
| 75 | 9.58 | 207.5 | -2.0 | 3.48e-4 | 0.0187 | 5.6 |
| 80 | 9.56 | 197.5 | -2.0 | 3.64e-4 | 0.0191 | 5.5 |
| 90 | 9.53 | 186.9 | 0.0 | 4.84e-4 | 0.0220 | 4.8 |

Example at 60 km: g/cp = 9.62/1004.7 = 0.00957 K/m; dT/dz = -0.0028 K/m; sum = 0.00677; times g/T = 9.62/245.5 = 0.03919 gives 2.66e-4. Conclusion from the table: N varies by only a factor of 1.35 over 60-90 km, so a single N of about 0.019 per second is adequate for the model; the 1976 profile is a global mean and the real mesopause near 85-90 km moves with season (about 10 K colder with larger gradient in summer).

### Gaps
- Mid-latitude seasonal N profiles from lidar or NRLMSIS were not retrieved; the existing wake model uses NRLMSIS/WACCM, which can supply N directly at 40.7 N.

## 3. Vertical diffusivity: layer turbulence versus effective transport

### Takeaway
Local turbulent eddy diffusivity at these altitudes is about 0.5-15 m^2/s from the dissipation-rate relation (CALC), consistent with the earlier 2-10 m^2/s from Megner. But eddy-flux lidar measurements and gravity-wave-dissipation transport give effective vertical diffusivities of order 100-300 m^2/s near the mesopause, one to two orders larger, and these are NOT the same quantity as the D_v that controls a sheared trail.

### Cited Findings
- Guo, Liu and Gardner (2017, GRL), Andes Lidar Observatory, 150 h, 85-100 km: eddy heat flux about -0.5 K m/s at 85 km to -0.1 K m/s at 100 km; mean turbulence thermal diffusivity 43 m^2/s and mean energy dissipation rate 37 mW/kg; peak heat flux -0.6 K m/s at 93 km — [ERAU record](https://portfolio.erau.edu/en/publications/first-na-lidar-measurements-of-turbulence-heat-flux-thermal-diffu/).
- Guo and Liu (JGR Atmospheres), over 400 h: annual mean thermal diffusivity from enthalpy flux 303 m^2/s at 85 km, minimum 221 at 90 km, 350 at 99 km — [same search results](https://core.ac.uk/download/217178574.pdf).
- Gardner, Guo and Liu (2019, Earth and Space Science): total dynamical diffusivity from dissipating waves plus turbulence can exceed 300 m^2/s near the mesopause. Gardner and Liu (2014, JGR): lidar needs 20 h or more of averaging for significant eddy flux estimates — [ERAU](https://portfolio.erau.edu/en/publications/vertical-dynamical-transport-of-mesospheric-constituents-by-dissi/).
- Liu and Gardner (2004, JASTP) assumed 200 m^2/s when comparing sodium flux predictions — same ERAU search results.
- Hocking (1999, Earth Planets Space 51, 525) discusses K = c epsilon / omega_B^2 (omega_B is buoyancy frequency), and notes the relation raises many issues; Wilson (2004) review of MST radar diffusivity says estimates show a lack of consensus — [search summaries](https://www.jstage.jst.go.jp/article/eps1998/51/7-8/51_7-8_525/_article/-char/en). The full text was paywalled; I did not obtain the value of c (commonly 0.2-0.8, MEMORY).

### Inferences
- CALC using K = 0.2 epsilon / N^2 (coefficient 0.2 is MEMORY, label as assumed): epsilon = 1 mW/kg = 1e-3 W/kg, N = 0.02: K = 0.2 x 1e-3 / 4e-4 = 0.5 m^2/s. epsilon = 30 mW/kg: K = 0.2 x 0.03 / 4e-4 = 15 m^2/s. So the earlier epsilon range of 1-30 mW/kg gives local K of 0.5-15 m^2/s, matching Megner's 2-10. Note the Guo et al. 37 mW/kg is within the range but their 43 m^2/s uses a larger coefficient and measured flux; the discrepancy is a factor of about 3 and is not resolved here.
- Ozmidov scale L_O = sqrt(epsilon/N^3). CALC N = 0.02, N^3 = 8e-6: epsilon = 1e-3 gives sqrt(125) = 11 m; epsilon = 0.03 gives sqrt(3750) = 61 m. So isotropic turbulence exists only below tens of metres; above that, flow is stratified.
- The D_v in the shear formula is the diffusivity for spreading over the plume's own vertical thickness (tens of metres to kilometres), acting over hours. Reasonable choices: 0.1-1 m^2/s for quiet nighttime layers between turbulent patches, 1-5 for the Megner-type background, and 10-300 only if gravity-wave-driven effective transport applies to the plume scale (unclear; gravity-wave "diffusivity" is a transport by wave dissipation, not spreading of a single tracer parcel). Warning: the 200-300 m^2/s lidar values are at 85-100 km and are heat or constituent transport; the 60-75 km range has no sourced values.
- Thin-layer sheet spacing: not found.

### Gaps
- Lindzen (1981), Garcia, Swenson, Fukao values and thin turbulent layer thickness/spacing (radar typically sees layers hundreds of metres to a few km thick, MEMORY) were not retrieved.
- No measured eddy diffusivity at 50-75 km mid-latitude found.

## 4. Horizontal diffusivity at 100 m to 100 km

### Takeaway
No measured mesospheric horizontal diffusivity was found. A scale-dependent estimate (Richardson-like, K_h about 0.1 epsilon^(1/3) l^(4/3)) gives tens to thousands of m^2/s at 1-10 km, but its applicability in the stratified mesosphere is uncertain, and the shear mechanism almost certainly dominates D_h at times beyond minutes.

### Cited Findings
- Lindborg (1999, J. Fluid Mech. 388, 259), MOZAIC aircraft data (troposphere/lower stratosphere, not mesosphere): third-order structure function positive in 30-300 km, enstrophy flux about 1.8e-13 per s^3; concludes the k^(-5/3) mesoscale range probably not explained by 2-D turbulence — [Cambridge abstract](https://resolve.cambridge.org/core/journals/journal-of-fluid-mechanics/article/can-the-atmospheric-kinetic-energy-spectrum-be-explained-by-twodimensional-turbulence/2835300488AEAA21B17FEF0B372FA3CC).
- Layered strongly stratified regime: layer thickness scales as U/N (Billant and Chomaz 2001 as cited in a conference abstract); Brethouwer, Billant, Lindborg and Chomaz (2007, JFM 585, 343) define two regimes separated by R = Re Fh^2 — [KTH PDF](https://www.mech.kth.se/~geert/papers/jfm585.pdf), [core.ac.uk](https://core.ac.uk/download/15496452.pdf). Not read in full; I did not obtain a diffusivity.

### Inferences
- CALC: the equivalent isotropic diffusivity of the existing model: sigma^2 = 2 D t, so D = (1000 m)^2 / (2 x 43,200 s) = 11.6 m^2/s for a 1 km radius at 12 h. That is the benchmark.
- CALC dimensional: K_h(l) = 0.1 epsilon^(1/3) l^(4/3) (coefficient assumed). For a mesoscale epsilon of 1e-4 W/kg (assumed, MEMORY-order, not mesospheric), epsilon^(1/3) = 0.0464; at l = 1 km, l^(4/3) = 1e4, K_h = 0.1 x 0.0464 x 1e4 = 46 m^2/s; at l = 10 km, l^(4/3) = 2.15e5, K_h = 1000 m^2/s. This is an estimate to show scale only, and depends on an assumed epsilon.
- Estimating from ocean analogy (Okubo diffusion diagram) was not sourced here.

### Gaps
- No mesospheric tracer-trail D_h found (candidates: Lloyd/Zimmerman rocket trail studies, TMA releases, Larsen tracer dispersion). Riley and Lindborg stratified-turbulence diffusivity not retrieved.

## 5. Effective horizontal spreading of a sheared line (calculations)

### Takeaway
Shear dominates for all plausible S and D_v. Along-shear width exceeds the model's 1 km 12-hour radius by a factor of about 12 to about 500, and already at 10 min a 0.01 per second shear gives widths of 40-270 m, comparable to a trail's initial size.

### Cited Findings
- No source is cited for the formula; it is the task statement's relation, sigma_s^2 = (2/3) D_v S^2 t^3. The altitudes 60, 75, 90 km do not enter the formula: results depend only on S, D_v, t. Altitude enters only through which S and D_v are realistic.

### Inferences
CALC. Example at S = 0.01, D_v = 1, t = 3600 s: sigma_s^2 = (2/3)(1)(1e-4)(3600^3 = 4.666e10) = 3.11e6, sigma_s = 1764 m. Widths in metres (1 sigma, along-shear direction):

| S (1/s) | D_v (m^2/s) | 1 min | 10 min | 1 h | 12 h |
|---|---|---|---|---|---|
| 0.005 | 0.1 | 1 | 19 | 279 | 11,592 |
| 0.005 | 1 | 2 | 60 | 882 | 36,656 |
| 0.005 | 5 | 4 | 134 | 1,972 | 81,966 |
| 0.01 | 0.1 | 1 | 38 | 558 | 23,184 |
| 0.01 | 1 | 4 | 120 | 1,764 | 73,313 |
| 0.01 | 5 | 8 | 268 | 3,944 | 163,932 |
| 0.03 | 0.1 | 4 | 114 | 1,673 | 69,551 |
| 0.03 | 1 | 11 | 360 | 5,291 | 219,938 |
| 0.03 | 5 | 25 | 805 | 11,831 | 491,797 |

Conclusion: at 12 h the steady-shear widths are 12 to 490 km, versus the model's 1 km. These 12-hour values assume constant S for 12 h, which is unrealistic (section 7); they are upper bounds. Even the 1-hour widths (0.3-12 km) exceed the model's 1 km. Vertical spread for comparison, sqrt(2 D_v t): D_v = 0.1, 1, 5 gives 11, 35, 77 m at 10 min; 27, 85, 190 m at 1 h; 93, 294, 657 m at 12 h. So the trail stays only hundreds of metres thick vertically but is stretched tens of km horizontally: a thin ribbon, not a round cylinder. Validity: sigma_s grows past 100 km only if S stays coherent across the vertical extent; the Gaussian assumption is addressed in section 7. The formula does not depend on altitude, but N and Ri limit S (S below about 0.04).

### Gaps
- The altitude dependence of realistic S and D_v is the open part (sections 1 and 3).

## 6. Trail length (40 m to 20 km) and orientation against the wind

### Takeaway
Orientation matters in a simple geometric way, and a sloped trail's own vertical extent can matter more than D_v. No sourced treatment of finite near-horizontal trails in shear was found; the following is geometry (CALC), not literature.

### Cited Findings
- None found specific to finite reentry trails in mesospheric shear. Aircraft contrail literature (anisotropic wake models) was not part of this search.

### Inferences
- A shear velocity difference acts between parts of the trail at different heights. If a trail is exactly horizontal, vertical thickness is only the small diffusive thickness (tens to hundreds of metres), and displacement across it is S x thickness: 0.01 x 100 m = 1 m/s, giving about 3.6 km per hour. If the trail descends at flight-path angle gamma, vertical extent = L sin gamma. CALC: L = 20 km, gamma = 5 deg: 20,000 x 0.0872 = 1.74 km; velocity difference across it is 0.01 x 1,740 = 17 m/s, so the ends separate by 17 m/s x 3,600 s = 63 km in one hour. So initial vertical extent dominates over D_v for long, steeply inclined trails; 40 m trails see S x 40 m x t = 0.01 x 40 x 3600 = 1.4 km in an hour.
- If the trail axis is parallel to the shear vector (wind direction at top versus bottom), shear stretches the trail along its own axis; if perpendicular, the trail axis is preserved and the cross-axis width grows as sigma_s above. For a general angle, the cross-section spreads along the shear direction projected onto the plane perpendicular to the axis. The model should therefore carry the angle between trail axis and the shear vector.
- Wind direction changes with height (tidal and wave winds rotate), so the shear vector direction can differ from the mean wind; unknown here.

### Gaps
- No literature on this specific problem was located; the anisotropic contrail-model references (Schumann; Lewellen) were not searched.

## 7. Validity of a Gaussian shear-diffusion model with time-varying shear

### Takeaway
A constant-S t^3 formula overestimates spreading when shear oscillates with tides or waves; the correct treatment is an effective diffusivity growing linearly in time with a reduced coefficient, and the horizontal shape need not stay Gaussian. The classic reference is Young, Rhines and Garrett (1982, JPO, "Shear-flow dispersion, internal waves and horizontal mixing in the ocean", MEMORY) for oscillating shear; this was not retrieved and must be checked.

### Cited Findings
- Young and Jones (1991, Physics of Fluids A, "Shear dispersion") review: the Taylor (1953) 1-D diffusion equation is only an asymptotic large-time result — [search summary, DOI 10.1063/1.858090](https://data.zhenzhi.com/details/1/1003123852). Smith's papers on contaminant dispersion in oscillatory flows (JFM 1982) exist; text not retrieved.
- Estuary literature: oscillatory shear dispersion matters most when the vertical mixing time is comparable to the oscillation period; spring-tide dispersion about 130 m^2/s in the Delaware estuary — [core.ac.uk 222884917](https://core.ac.uk/download/222884917.pdf), [core.ac.uk 4168704](https://core.ac.uk/download/4168704.pdf).

### Inferences
CALC (derived from the shear-diffusion random-walk model, assumption: shear S(t) = S0 cos(omega t), uniform in height, vertical random walk with D_v; not checked against Young et al.): displacement per unit height after one oscillation is (S0/omega) sin(omega t), so the lasting horizontal diffusivity is D_eff of about D_v S0^2 / (2 omega^2), valid when omega is much larger than D_v / h^2. Example, S0 = 0.01 per second, 12 h period, omega = 2 pi / 43,200 = 1.454e-4 per second: S0/omega = 68.8, squared 4,731, half 2,366. D_v = 1 gives D_eff = 2,366 m^2/s; horizontal sigma at 12 h = sqrt(2 x 2,366 x 43,200) = 14.3 km (versus 73 km constant-shear). D_v = 0.1: 237 m^2/s, 4.5 km. D_v = 5: 11,830 m^2/s, 32 km. With S0 = 0.005, divide D_eff by 4. Check of validity: D_v / h^2 = 1 / (2000 m)^2 = 2.5e-7 per second, much less than 1.45e-4, satisfied. Also reversible part: tilt amplitude is (S0/omega) x vertical spread, 68.8 x 294 m = 20 km oscillating with the tide, so the trail leans back and forth by up to about 20 km at 12 h while the irreversible spreading is the D_eff part. Shorter-period gravity-wave shear (period 30-60 min) has far smaller S0/omega and much smaller D_eff; at a buoyancy-limited period of 6 min the factor is (0.01 / 0.0175)^2 / 2 = 0.16, giving D_eff of about 0.16 D_v, a small effect. With vertical wavelengths of order 10 km, the sheet is not coherently sheared across more than half a wavelength, so D_eff is also bounded by wavelength; I did not quantify this.
Practical recommendation (inference): replace sigma_s^2 = (2/3) D_v S^2 t^3 by the time integral of a time-dependent shear, x = z integral of S dt (non-Gaussian when S changes sign), or by D_eff for oscillatory shear; and cap by S at Ri about 1/4 (S at most 2N, about 0.035-0.044 per second). For a first-pass model, use S_rms of 0.005-0.01 per second with the oscillating form (tidal) and report the bounds given above.

### Gaps
- Young, Rhines and Garrett (1982) and Smith (1982) texts not obtained, so the factor 1/2 and applicability conditions are my derivation, not a verified citation.
- Whether horizontal shape remains Gaussian under varying shear with vertical-scale-limited correlation was not located in the literature.

## Summary table of recommended input ranges (inference, not measurement)
S about 0.005-0.01 per second typical, 0.03 strong (Ri about 0.3-0.5 per CALC), 0.04+ unstable; N about 0.016-0.022 per second (CALC, US Standard); D_v 0.1-5 m^2/s local turbulence (CALC from 1-30 mW/kg), with 100-300 m^2/s only for wave-driven transport at 85-100 km (cited, different quantity); D_h not measured, benchmark of the current model 11.6 m^2/s (CALC); at 40.7 N the lidar mid-latitude evidence is Colorado State sodium lidar at about 80-105 km only.
