# Smallest turbulent eddies in the undisturbed mesosphere, 50-90 km: Kolmogorov microscale and measured inner scale

Method note. Everything labelled "[calc]" is my own arithmetic using the 1976 US Standard Atmosphere (temperature, density) and Sutherland's law for dynamic viscosity (mu = 1.458e-6 T^1.5 / (T + 110.4) Pa s), kinematic viscosity nu = mu/rho. Everything with a link is read from the cited source (full text where I could extract it: Strelnikov 2017, Lehmacher 2018, Strelnikov 2009, Hall 2003; abstract only for Lubken 1993, 1997, 2002, Latteck 2005, Vlasov and Kelley 2015). I could not open the full text of Lubken 1992/1993/1997/2002 or the Szewczyk 2015 thesis (12 MB, over the fetch limit), so the summer/winter climatology profiles are known to me only through other papers quoting them. All measured profiles are from high northern latitude (Andoya 69N or Esrange 68N, Poker Flat 65N, Tromso 69N) unless stated; I found no measured low- or mid-latitude epsilon profiles.

## 1. Measured energy dissipation rates (epsilon) from 50 to 90 km

### Takeaway
At high northern latitudes the measured epsilon is about 1 mW/kg near 70 km, about 10 mW/kg near 80 km and about 50 mW/kg near 90 km in winter, and roughly 10 mW/kg as a mean (spread over seven decades) at 80-90 km in summer; below about 70 km there are almost no rocket data and only radar estimates (about 5 mW/kg at 60 km).

### Cited Findings
- Lubken 1993: 17 rocket profiles; for the 12 winter flights at 69N, mean epsilon was 1-2 mW/kg in the lower mesosphere (below 75 km) and 10-20 mW/kg below 100 km, equal to heating of about 0.1 K/day and 1 K/day; this is at least one order of magnitude below most earlier measurements — [Lubken et al. 1993, JGR 98 (abstract)](http://earthref.org/ERR/58965/)
- Lubken 1997: 22 rockets in 6 years, 19 at high latitude (17 from Andoya 69N, 2 from Esrange 68N); winter heating about 0.1 K/day below about 75 km and 1-2 K/day above; summer mean heating 10-20 K/day near the mesopause (about 90 km); turbulence occupies 78-97 km in summer but 60-100 km in winter — [Lubken 1997, JGR 102 (abstract)](http://earthref.org/ERR/52842/)
- Szewczyk (2015) updated winter statistics at Andoya: mean epsilon rises continuously from 1 mW/kg at 70 km to 10 mW/kg at 80 km, up to a maximum of 50 mW/kg at 90 km; the total variability envelope covers almost 4 orders of magnitude (quoted by Lehmacher 2018, I did not see the thesis) — [Lehmacher et al. 2018, Ann. Geophys. 36, 1099](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- MTeX, Poker Flat Alaska (65N), winter, 2015 launches: layers at about 71 km and 75 km with epsilon 0.18-6.4 mW/kg (medians 1.5 and 2.7 mW/kg); downleg maximum 36 mW/kg at 72.8 km (layer median 2.2 mW/kg); layer at 80-81 km up to 1.3 mW/kg; 83-84 km up to 7-12 mW/kg; "very little turbulence above 82 km" on that night; abstract summary: 0.1-10 mW/kg, "typical for the winter mesosphere" — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Summer, Andoya, ECOMA 2007: upper turbulence region (89-96 km) epsilon about 10 to about 2000 mW/kg; lower region (79-86 km) several layers with values typical for 69N summer; the ECOMA in situ range was 1 to 1000 mW/kg; ALWIN radar gave about 16 mW/kg on average (10-40) at 83-89 km — [Strelnikov et al. 2009, Ann. Geophys. 27, 1449](https://www.ann-geophys.net/27/1449/2009/angeo-27-1449-2009.pdf) (read from a local copy of the PDF)
- Summer, Andoya 69N, WADIS-1 (2013): in situ epsilon spans five decades; mean over 80-90 km: upleg 74 mW/kg (Heisenberg model) and 73 mW/kg (Tatarskii model), downleg 492 and 307 mW/kg; MAARSY radar 58 mW/kg, EISCAT 769 mW/kg; thin sub-layers 85.5-85.8 km up to about 7 mW/kg (text as extracted; the quoted downleg mean is far above this, so "up to about 7" probably is a typo for a larger value or a different unit, treat with caution); the statistics of all summer in situ flights (Szewczyk 2015) at 80-90 km have a mean of about 10 mW/kg and a width of one decade (sigma about 1 in log10) spanning seven decades — [Strelnikov et al. 2017, Ann. Geophys. 35, 547](https://angeo.copernicus.org/articles/35/547/2017/)
- Radar-measured epsilon varied by up to 2 orders of magnitude over 24 h (tidal) in WADIS-1 — [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/)
- Lower bound for detectable turbulence: epsilon_min = nu * omega_B^2 (omega_B the Brunt-Vaisala frequency) — [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/)
- ECOMA 2010 (Andoya, November-December, winter-like): peak epsilon about 2 W/kg (2000 mW/kg) in the strong upper mesospheric inversion layer, 86-89 km — [Szewczyk et al. 2013 abstract](https://angeo.copernicus.org/articles/31/775/2013/)
- Radar (3.17 MHz Doppler, Andoya): mean epsilon about 5 mW/kg at 60 km and about 20 mW/kg at 80 km, "in reasonable agreement" with the rocket climatologies (a radar upper-limit type estimate; no season given) — [Latteck, Singer, Hocking 2005, EGU abstract](https://meetings.copernicus.org/www.cosis.net/abstracts/EGU05/09987/EGU05-J-09987.pdf)
- Radar epsilon at 80-90 km in summer: 5-100 mW/kg (MIDAS/MaCWAVE 2002, Engler et al. 2005, quoted via search summary only) and earlier radar medians 0.05-0.1 mW/kg (Hocking 1985, quoted in Strelnikov 2017) — [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/)
- MF-radar fading-time method gives only upper limits of epsilon, because the largest-scale time is assumed to be 1/N and gravity waves leak in — [Hall et al. 2003, Ann. Geophys. 21, 787](https://angeo.copernicus.org/articles/21/787/2003/angeo-21-787-2003.pdf) (read locally)

### Inferences
- Season and altitude contrast. Winter mean epsilon grows from about 1 mW/kg at 70 km to about 50 mW/kg at 90 km; summer has none or little in the lower mesosphere (see section 4) and a mean about 10 mW/kg at 80-90 km with extreme layers of 0.1 to 1000+ mW/kg. For single-valued use, epsilon of 1-10 mW/kg across 70-90 km is the measured "typical" range, with the true mean over a profile dominated by rare strong layers.
- Heating check: 1 mW/kg corresponds to about 0.09 K/day [calc: 1e-3 W/kg / 1004 J/kg/K x 86400 s = 0.086 K/day], consistent with Lubken 1993's "0.1 K/day".

### Gaps
- I could not obtain the Lubken 1997 or Lubken et al. 2002 tabulated profiles; the data are in the Szewczyk 2015 thesis (Rostock, Dissertation_Szewczyk_2015.pdf, 12 MB) and Figure 12 of Strelnikov 2017. Numbers by altitude in 5 km steps for summer are therefore not available to me.
- No measured epsilon from rockets at 50-65 km that I could find except the 17 winter profiles of Lubken 1993/1997 (quoted as "1-2 mW/kg below 75 km"). I did not find MaCWAVE, Becker, Haack or TURBO numbers in text I could read.
- No low- or mid-latitude in situ epsilon found.

## 2. Kolmogorov length and inner scale: definitions, ratio, and values in metres

### Takeaway
Rocket experiments fit the Heisenberg (1948) model, whose inner scale is l0 = 9.9 (nu^3/epsilon)^(1/4), i.e. 9.9 times the Kolmogorov microscale eta; reported mesospheric inner scales are 10-50 m, and the Kolmogorov microscale is therefore about 1-5 m. Because nu grows about e-fold per 8 km (density falls), l0 grows strongly with altitude: about 5-10 m near 60 km to 50-150 m near 90 km for epsilon of 1-50 mW/kg [calc].

### Cited Findings
- Heisenberg-model inner scale: l0 = 9.9 (nu^3/epsilon)^(1/4); epsilon follows from the fitted l0; "Typical inner scales in the mesosphere and lower thermosphere are 10-50 m, which requires measurements at meter-scale resolution" — [Lehmacher et al. 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Worked example at 71.0-71.1 km, winter, Poker Flat: fitted inner scale 21 m (range 17-27 m) corresponds to epsilon of 1.5 mW/kg (range 0.68-3.9 mW/kg), because epsilon goes as the fourth power of l0 — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Spectral fits use Kolmogorov-Heisenberg slopes of -5/3 (inertial) and -7 (viscous); the viscous range lies at 31.6-316 Hz for the payload speed of about 1000 m/s, i.e. scales of about 3 to 30 m (f = v/l) — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Heisenberg and Tatarskii models give very close epsilon in WADIS-1 (74 versus 73 mW/kg upleg), but the Tatarskii model is more accurate when the viscous subrange is resolved over more than 2 decades — [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/); [Strelnikov et al. 2022 (via search summary)](https://essopenarchive.org/users/527614/articles/596672-assessment-of-precision-of-spectral-model-turbulence-analysis-technique-using-dns-data)
- Older radar/rocket estimates: Watkins et al. 1988 microscale 1.5-2.0 m at 80-87 km with l0 about 13 eta (20-26 m); Hall 2000 suggests inner scale about 10 m in moderate turbulence and "a few metres" depending on strength (both from a search-engine summary, not read at source) — [search result](https://nipr.repo.nii.ac.jp/record/6317/files/KJ00000048057.pdf)
- Lubken 1992 (JGR 97, 20385, not GRL) fits the Tatarskii spectral model for inertial-convective plus transition to viscous-diffusive range; epsilon from the model exceeds the older C_n^2 method by factors 1.2 and 2.4 in two flights — [Lubken 1992 abstract](https://earthref.org/ERR/58619/)
- Ratio for the Tatarskii form: l0(Tatarskii) = 7.06 eta is a value I recall from Lubken (1992) and Strelnikov et al. (2003) but I did not find it in a source I could read — treat as unverified.

### Inferences (calculations)
- [calc] Standard atmosphere and Sutherland viscosity (rho in kg/m^3, nu = mu/rho in m^2/s): 

| Altitude | T (K) | rho | nu (m^2/s) |
|---|---|---|---|
| 50 km | 270.6 | 1.03e-3 | 0.017 |
| 60 km | 247.0 | 3.10e-4 | 0.051 |
| 70 km | 219.6 | 8.28e-5 | 0.174 |
| 80 km | 198.6 | 1.85e-5 | 0.716 |
| 90 km | 186.9 | 3.42e-6 | 3.67 |

 Arithmetic example, 80 km: mu = 1.458e-6 x 198.6^1.5 / (198.6 + 110.4) = 1.321e-5 Pa s; nu = 1.321e-5 / 1.846e-5 = 0.716 m^2/s. Then for epsilon = 10 mW/kg: eta = (0.716^3 / 0.01)^(1/4) = (36.7)^(1/4) = 2.46 m; l0(Heisenberg) = 9.9 x 2.46 = 24.4 m; 7.06 x 2.46 = 17.4 m.
- [calc] Kolmogorov microscale eta / Heisenberg l0 = 9.9 eta (all in metres) for fixed epsilon:

| Altitude | eps 0.1 mW/kg | eps 1 mW/kg | eps 10 mW/kg | eps 100 mW/kg |
|---|---|---|---|---|
| 50 km | 0.46 / 4.6 | 0.26 / 2.6 | 0.15 / 1.4 | 0.08 / 0.8 |
| 60 km | 1.08 / 10.6 | 0.60 / 6.0 | 0.34 / 3.4 | 0.19 / 1.9 |
| 70 km | 2.69 / 26.6 | 1.51 / 15.0 | 0.85 / 8.4 | 0.48 / 4.7 |
| 80 km | 7.78 / 77 | 4.38 / 43 | 2.46 / 24.4 | 1.38 / 13.7 |
| 90 km | 26.5 / 262 | 14.9 / 148 | 8.38 / 83 | 4.71 / 46.6 |

 Each cell is eta / l0 in metres. Reading: l0 reaches 30 m at 71.5 km for 0.1 mW/kg, 77 km for 1 mW/kg, 82 km for 10 mW/kg and 87 km for 100 mW/kg [calc]. The values at a given altitude vary by a factor of 10^(1/4) = 1.8 per decade of epsilon.
- [calc] Altitude-dependent function using measured climatology epsilon (winter means 1, 10, 50 mW/kg at 70, 80, 90 km from Szewczyk 2015; 5 mW/kg at 60 km radar from Latteck 2005): eta = 0.40, 1.51, 2.46, 5.6 m and l0 = 4.0, 15, 24, 55 m at 60, 70, 80, 90 km respectively (e.g. 90 km: nu = 3.665, nu^3 = 49.2, /0.05 = 984, ^(1/4) = 5.6 m, x 9.9 = 55 m).
- [calc] Check against the Lehmacher worked example: with 21 m and 1.5 mW/kg, the implied nu is (eta^4 eps)^(1/3) = 0.31 m^2/s, whereas my standard-atmosphere nu at 71 km is about 0.20 m^2/s. With nu = 0.20 an inner scale of 21 m would mean epsilon about 0.4 mW/kg, not 1.5. This discrepancy (factor 1.5 in nu, 3.8 in epsilon) is probably the actual measured density/temperature of that cold-winter profile; I could not resolve it.
- [calc] Detectable-turbulence floor epsilon_min = nu N^2 with N^2 = 4e-4 s^-2 (my assumed value, N about 0.02 s^-1): 0.02 mW/kg at 60 km, 0.07 at 70 km, 0.29 at 80 km, 1.5 mW/kg at 90 km. At 90 km this equals the lower end of the measured winter mean, so quiet turbulence there is near the limit where buoyancy suppresses it (buoyancy scale, not viscosity, would then limit eddies).

### Gaps
- No source I could read states the l0/eta ratio for the Tatarskii form; the 7.06 value is from memory.
- Lubken 1992/1993 originally tabulated l0 profiles; I could not read them.

## 3. Intermittency, patchiness and what a non-turbulent patch does

### Takeaway
Turbulence is intermittent in space and time; in winter 90-100 km it is observed only 16 % of the time, summer turbulence sits in thin, strongly varying layers (sub-layers 300 m thick; 7 of 10 PMSE layers had none), and a non-turbulent volume has zero dissipation, with fossil turbulence and gravity-wave structures left behind.

### Cited Findings
- Winter, 90-100 km: small-scale isotropic turbulence is observed with only 16 % probability (Szewczyk 2015), though epsilon is most likely 10-100 mW/kg when present; strong isotropic turbulence is rarely observed above 95 km — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Lubken 2002: of 10 PMSE layers observed with 8 summer rocket flights at 69N/68N, turbulence was absent in 7; no correlation between PMSE and neutral turbulence — [Lubken, Rapp, Hoffmann 2002, JGR (abstract)](https://zendy.io/title/10.1029/2001jd000915)
- WADIS-1 summer: two main turbulent layers 82-85 km and 85-90 km; sub-layers only 300 m thick (85.5-85.8 km); epsilon distribution log-normal with sigma(log10) about 1 for single in situ points, 0.489 for radar single points and 0.5 for radar vertical means; "turbulence appears mainly in the range 82-95 km"; non-turbulent parts of a radar volume count as zero dissipation, and PMSE is an acceptable proxy for "active or fossil" turbulence — [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/)
- Neutral-density tracer technique cannot see old, already well-mixed (adiabatic) turbulent layers, since persistent mixing washes out the density fluctuations — [Strelnikov 2009](https://www.ann-geophys.net/27/1449/2009/angeo-27-1449-2009.pdf) and [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/)
- Turbulent layers at 71-75 km changed completely over 33 minutes between two Poker Flat flights (a layer at 80-81 km absent on first flight, strong on second); gravity-wave breaking simulations evolve over about five buoyancy periods in 30 min — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Radar estimates may be upper limits because turbulence is "intermittent in both space and time" (Roper 1998, via search summary) — [search result](https://acnpsearch.tweb-dev.unibo.it/singlejournalindex/51823655)

### Inferences
- Non-turbulent patches: no eddy mixing, so the plume would spread only by molecular diffusion (plus wave/shear dispersion, not covered by my sources); molecular diffusivity of a trace vapour is of order nu (Schmidt number order 1 for atoms), i.e. 0.02-4 m^2/s at 50-90 km [calc: equals nu above, order of magnitude] — equal to or larger than Megner's K at 50-70 km, which means Megner's 2-10 m^2/s is not clearly above molecular diffusion below 75 km. This point is my inference; no source states it.
- A time- and volume-averaged diffusivity must be a fraction (occurrence x in-layer K) of the turbulent value, so using a continuous switch-on in the model is consistent with measured intermittency only if K is the average.

### Gaps
- I could not find a percentage of altitude turbulent for summer or for 60-90 km in winter, nor a turbulence-fraction statistic for Lubken 1997/2002 (it should be in the Szewczyk thesis, section "Turbulence occurrence rate").

## 4. Summer versus winter, and the lower mesosphere (50-65 km)

### Takeaway
Summer turbulence is confined to roughly 78-97 km but is much stronger there (mean heating 10-20 K/day at 90 km), while winter turbulence covers 60-100 km but weaker; below about 65 km rocket data are very sparse and only radar estimates (about 5 mW/kg at 60 km) exist.

### Cited Findings
- Summer 78-97 km vs winter 60-100 km; none of the 7 summer flights in Lubken 1997 detected noticeable turbulence in the middle and lower mesosphere — [Lubken 1997](http://earthref.org/ERR/52842/)
- MaCWAVE/MIDAS (Andoya summer 2002): near continuous turbulent layers at 72-90 km; first evidence of summer turbulence below 82 km; linked to gravity-wave breaking and the zonal wind maximum height — [Rapp et al. 2004, GRL 31, via search summary](https://vufind-test.katalog.k.utb.cz/SummonRecord/FETCH-LOGICAL-a4753-9300d5828ad947338892e8f3a0d631f006fc35f272f47abb299ed6121057445b3)
- WADIS-1 summer upleg also saw turbulence layers between 82 and 70 km, "unusual for summer" — [Strelnikov 2017](https://angeo.copernicus.org/articles/35/547/2017/)
- Spring transition (May 2000, Andoya): the turbulent structure changed within about 10 days — [Mullemann et al. 2002, via search summary](https://earthref.org/ERR/77810)
- Turbulent heating in the winter mesosphere is small (0.1 K/day below 75 km); summer polar mesopause heating is largest — [Lubken 1997](http://earthref.org/ERR/52842/)
- Radar 60 km: 5 mW/kg — [Latteck 2005](https://meetings.copernicus.org/www.cosis.net/abstracts/EGU05/09987/EGU05-J-09987.pdf)
- Winter 60-75 km rocket mean 1-2 mW/kg — [Lubken 1993](http://earthref.org/ERR/58965/)

### Inferences
- At 50-65 km the typical epsilon is 1-5 mW/kg (radar, upper-limit-like; rocket values lower) and then, from my calculation, eta is about 0.4-0.6 m and l0 about 4-6 m at 60 km [calc, section 2 table: 1 mW/kg gives eta 0.60 and l0 6.0 m; 5 mW/kg gives eta 0.40 and l0 4.0 m]; at 50 km eta is 0.15-0.26 m and l0 1.4-2.6 m for 1-10 mW/kg.
- Season matters most for altitude range: the choice of one l0 for both seasons hides that summer has a nearly turbulence-free lower mesosphere.

### Gaps
- I found no sub-60 km rocket or radar epsilon with numbers by season; no data from the lidar/ MST radar literature on 50-60 km turbulence; I did not search Becker/Haack/TURBO papers in depth (searches did not return their text).

## 5. Eddy diffusion coefficient K from epsilon versus Megner's 2-10 m^2/s

### Takeaway
Measured values are about 2 to 100 times larger than Megner's profile: Lubken 1997 quotes typical K of 4 m^2/s in the lower mesosphere and 100 m^2/s in the upper mesosphere, whereas Megner's 2 m^2/s at 60 km, 6 at 75 km and 10 at 90 km (as given in the task) are close to the low end only below about 70 km.

### Cited Findings
- Lubken 1997: typical turbulent velocity 0.3 m/s in the lower mesosphere and 1-3 m/s in the upper mesosphere; typical K = 4 m^2/s and 100 m^2/s, with no comparable seasonal variation — [Lubken 1997](http://earthref.org/ERR/52842/)
- Upper limit on K from similarity theory: about 3.2e6 cm^2/s (320 m^2/s) for epsilon at most 2 W/kg; K from rockets at high latitude (Lubken 1997) satisfy the criterion in summer and winter; TOMEX K of 1e7 cm^2/s (1000 m^2/s) does not — [Vlasov and Kelley 2015, Ann. Geophys. 33, 857](https://angeo.copernicus.org/articles/33/857/2015/)
- Eddy diffusion in 1-D models acts as a parameterization of meridional circulation transport more than of small-scale mixing (1988 study, via search summary) — [search result](https://earthref.org/ERR/57253)

### Inferences
- [calc] Using K = 0.81 epsilon / N^2 (the standard Weinstock form I recall; coefficient not verified in a source I read) with N^2 = 4e-4 s^-2 (assumed, N about 0.02 s^-1): epsilon 1 mW/kg gives K = 2.0 m^2/s; 10 mW/kg gives 20 m^2/s; 50 mW/kg gives 101 m^2/s (0.81 x 0.05 / 4e-4 = 101). So winter-mean K is about 2 at 70 km, 20 at 80 km and 100 at 90 km, against Megner's about 4, 8 and 10 (interpolated). Megner's value is about equal at 70 km, 2-3 times too small at 80 km, and about 10 times too small at 90 km relative to the rocket mean. If the plume uses the full measured K the plume grows faster than in the current model at upper altitudes.
- Intermittency and radar values point the other way (turbulent fraction below 1), so the effective mean K could be several times lower than the in-layer values.

### Gaps
- I did not obtain Megner et al. 2006 myself; the 2/6/10 m^2/s figures are from the task statement. Lubken 1997 profile of K by altitude not read.

## 6. Does l0 vary strongly with altitude, and what is a defensible l0(z)?

### Takeaway
Yes: kinematic viscosity rises about e-fold per 8 km, so l0 rises about e-fold per 10 km at fixed epsilon (about 5-10 m at 60 km to 50-150 m at 90 km), and a single 30 m is only right around 77 km for 1 mW/kg; it is too large below about 70 km, and too small above about 85 km.

### Cited Findings
- Kinematic viscosity increases exponentially with the scale height — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)
- Typical measured inner scales: 10-50 m in the mesosphere and lower thermosphere (21 m at 71 km) — [Lehmacher 2018](https://angeo.copernicus.org/articles/36/1099/2018/angeo-36-1099-2018.pdf)

### Inferences (calculation)
- [calc] At fixed epsilon, l0 proportional to nu^(3/4); a fit across 50-90 km gives an e-folding length of 9.98 km for both 1 and 10 mW/kg (l0 at 70 km: 16.8 m for 1 mW/kg and 9.4 m for 10 mW/kg).
- [calc] Defensible function: l0(z) = 9.9 (nu(z)^3 / epsilon(z))^(1/4), with nu from the standard atmosphere and epsilon(z) = 1 mW/kg at 70 km, using these anchors (winter measured means): epsilon about 5 mW/kg at 60 km, 1 mW/kg at 70 km, 10 mW/kg at 80 km, 50 mW/kg at 90 km; giving l0 of about 4, 15, 24 and 55 m. For 50 km, use l0 about 1-3 m (epsilon 1-10 mW/kg; value is an extrapolation with no measurement behind it). A simple interpolation is l0 about 4 m at 60 km, 15 m at 70 km, 25 m at 80 km, 55 m at 90 km, i.e. roughly 15 x exp((z - 70)/14 km) between 60 and 90 km [calc: 15 e^(20/14) = 62 m, close to 55].
- The 30 m value is within the measured range (10-50 m) and equals the l0 expected at about 77-82 km for 1-10 mW/kg; it is too large by factor 2-8 at 50-65 km and too small by factor 2-4 at 90 km. The user's suspicion of "far too large" is right only below about 70 km; the Kolmogorov microscale eta (the actual smallest eddy, about l0/10) is 0.3-0.6 m at 60 km, 1.5 m at 70 km, 2.5-4 m at 80 km and 5-15 m at 90 km [calc].
- For a plume of radius 8 mm growing to 1 km, plume width exceeds l0 (4-55 m depending on altitude) at very different stages: at 60 km after the plume is about 4 m wide, at 90 km at about 55 m; the switch should be altitude-dependent. Also note the plume is far below eta (0.3 m to 15 m) initially, so molecular diffusion alone applies at those scales.

### Gaps
- No source states a recommended l0(z) formula; the function above is my construction from measured epsilon anchors, with the 60 km anchor from a radar and the 50 km value not measured.
- Mean free path at 90 km is about 2 mm (not read from a source) — still far below eta.
