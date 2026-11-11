# Wake radius growth laws below 5.8 km/s

> **Update, 2026-10-29 (later the same day): the full texts of Taylor et al. (1964) and Lees (1964) were read**
> (PDFs saved by Asha from MIT Libraries in `Literature Review/Wake Growth/`). Section 0 below supersedes the
> "abstract only" statements about those two papers in sections 1, 2, 4.1, 4.2 and 7.

## 0. Primary-source check of Taylor et al. (1964) and Lees (1964)

**Taylor, Melcher and Washburn (1964), AIAA Journal 2(10), 1731-1738 (full text read, page images checked).**
Equation (4), printed page 1737, is w/(C_D A)^(1/2) = K [x/(C_D A)^(1/2)]^(1/3). The exponent is one third (the
abstract on printed page 1731 says the same). Fay suggested K = 0.9 from dimensional arguments; Taylor et al.
fitted K = 0.66 to their luminous-wake data for 2 <= x/(C_D A)^(1/2) <= 1200, with a root-mean-square deviation of
23 percent (equation 5). Here w is the full width of the luminous wake (so the radius is w/2), which the authors
conclude coincides with the viscous (turbulent) core for x/(C_D A)^(1/2) greater than about 7 (printed page 1735);
x is distance behind the body; C_D A is drag coefficient times frontal area. The data cover 0.22, 0.55 and 1.0 inch
lexan spheres and quasi-spheres at 12,000 to 20,500 ft/s (3.7 to 6.2 km/s) and freestream pressures of 0.2 to
10 cm Hg. **The law is for the turbulent regime only.** The paper divides regimes with the pressure parameter
p_inf (C_D A)^(1/2) (Goldburg's correlation, transition at a shoulder Reynolds number of 6000): turbulent above
2.0 cm Hg-cm, where the width does not depend on pressure; transitional from 0.5 to 2.0, where the width passes
through a minimum; laminar below 0.5, where the width grows again as pressure falls. For the laminar regime it
gives only a proportionality, equation (1): w/(C_D A)^(1/2) proportional to [1/(C_D A)^(1/2)]^(1/2) times
[nu (x/(C_D A)^(1/2)) / V_inf]^(1/2), with laminar diffusivity nu_L = mu/rho, approximately R T^1.76 / p_inf
(equation 2), and a constant of proportionality fitted to their data at p_inf (C_D A)^(1/2) = 0.275 cm Hg-cm (Fig. 8).
They note the one laminar data set shows a stronger pressure dependence (closer to 1/p) than the predicted
inverse square root.

**Lees (1964), AIAA Journal 2(3), 417-428 (full text read for the growth law and transition; page 2 image
checked).** Equation (1), printed page 418: the density-weighted turbulent-front half-width, scaled by
(C_D A)^(1/2), grows as [K x-bar C_Df/C_D]^(1/3), with x-bar = x/(C_D A)^(1/2); beyond x-bar of about 1000 nearly
all the body's drag has been swallowed and the half-width goes as (K x-bar)^(1/3) "regardless of the body shape or
the Reynolds number". No numerical K is given in the paper (K is inversely proportional to Townsend's universal
Reynolds number). On transition: turbulent wakes appear only above a critical Reynolds number, based on freestream
conditions ahead of the body and body diameter, of about 5 x 10^4 at about 20,000 ft/s; the abstract places this
at about 220,000 ft (67 km) for a 6 ft body, and the text at about 200,000 ft for a 1 ft nose diameter.

**Applied to the Generic Spacecraft** (equal-area diameter 2.81 m, (C_D A)^(1/2) about 2.35 m, US Standard
Atmosphere 1976, parent speed from the DRAMA file; computed 2026-10-29): the freestream Reynolds number is
3.0 x 10^4 at 80 km and 6.0 x 10^4 at 75 km, so by Lees' criterion the wake becomes turbulent near 75 to 78 km and
is laminar above. Taylor's pressure parameter is 0.19 cm Hg-cm at 80 km, 0.42 at 75 km, 0.92 at 70 km, 1.9 at
65 km and 3.9 at 60 km: laminar above about 74 km, transitional between about 74 and 63 km, turbulent below about
63 km. Both criteria were calibrated in ballistic ranges at room temperature and at 3.7 to 6.2 km/s, so the
boundaries are uncertain by several kilometres. Small aluminium sphere parents (5 to 400 mm) are laminar by both
criteria everywhere above 60 km (Reynolds numbers of about 20 to 2 x 10^4 at 3 km/s).

Literature search carried out on 2026-10-29 with web search and online reading only. The one exception to "online only" is Park (1971): the web-fetch tool cached the NASA scan locally, and its pages were rendered to images in the session scratchpad so that the equation could be read off the page itself rather than from the garbled text layer.

## 1. What was asked, and whether it is done

The request was for a published equation for how the radius of the turbulent wake behind a blunt body grows with distance, valid below about 5.8 kilometres per second. The second part of the request was to check Park's coefficient against the original scan.

Both parts are done, with one limit. Park's equation is now confirmed from the page image. I found one drag-based correlation that comes from ballistic-range data at 3.7 to 6.3 km/s, and it is the best candidate for the speed range in question. However, I could reach it only through its abstract, and the exponent is unreadable in the online copy of that abstract. Most of the primary 1960s sources sit on the Defense Technical Information Center (DTIC), the RAND Corporation site, Springer, ScienceDirect or the Canadian Science Publishing site. All of these refused the fetch tool with "403 Forbidden" or a login redirect, so I could not read their full texts. Every entry below says exactly what was read.

## 2. Summary in plain language

Behind a blunt body, the turbulent wake grows roughly as the cube root of the distance behind the body, so that every eightfold increase in distance doubles its width. Every source that states an exponent for a turbulent wake far downstream agrees on this, at hypersonic and at incompressible (low-speed) conditions alike. The sources disagree about the coefficient and about what sets it.

The ballistic-range correlation quoted by Taylor, Melcher and Washburn (1964) includes the body's drag coefficient, so it can be applied at any speed where that drag coefficient is known. It reads

  w / (C_D A)^(1/2) = K [ x / (C_D A)^(1/2) ]^t, with K = 0.66.

The exponent t is unreadable in the online abstract, and a secondary source gives it as one third.

Park's law has no drag coefficient in it. At equal frontal area and distance, it gives a wake radius 1.45 times larger than the Taylor correlation would give for a body with a drag coefficient of 1. The two laws agree only if the body's drag coefficient is about 3, which is unphysically large.

Two departures from the cube-root law matter for the small-sphere parent:

- **Low Reynolds number.** When the wake's local Reynolds number is low, the wake grows as the square root of distance instead of the cube root. This is true both for a viscous "low-Reynolds-number turbulent" state and for a fully laminar wake. At low Reynolds number the laminar wake has an exact textbook solution, read in full here, whose radius is twice the square root of (kinematic viscosity × distance ÷ speed).
- **Rarefied flow.** I found no published wake-growth law for flow that is rarefied, meaning the gas's mean free path is not small compared with the body.

## 3. Park (1971): result of the coefficient check

**Citation.** Park, C. (1971 on the NASA Technical Reports Server record; the scan carries no date on its title page), "Estimates of nitric oxide production for lifting spacecraft reentry" (title as printed on page 2 of the scan), NASA TM X-62,052, NASA Ames Research Center. URL: https://ntrs.nasa.gov/api/citations/19720024212/downloads/19720024212.pdf

**What was read.** I read the full text: rendered page images of PDF pages 14, 18, 19, 21, 25 and 26, together with the scan's text layer for the notation list and the summary.

**The equation is confirmed.** Equation (11) on printed page 17 (PDF page 18) reads

  r_b / r_o = 0.7 (s' / r_o)^(1/3).

The coefficient is 0.7 and the exponent is 1/3, exactly as previously reconstructed. The page says that "its boundary is taken to vary according to the formula (Ref. 2)."

**Symbols.**

- **r_b** is the radius of the wake boundary. The notation list on printed page 4 gives subscript b as "boundary of wake."
- **r_o** is the reference radius, defined on printed page 18 by π r_o² = frontal area of the vehicle. The notation list gives subscript o as "origin of wake."
- **s'** is "the distance measured from a fictitious source (see Fig. 9)" (printed page 18). Both lengths were in feet in the paper; the equation itself has no units because it is a ratio.

**Virtual origin.** Printed page 18 says: "In the present case, the fictitious source is 73 ft upstream of the initiation of the wake." The initial wake diameter is 50 ft (printed page 17).

My own arithmetic (not stated in the paper) shows that these numbers fit together only if r_o is about 25 ft. In that case the wake boundary radius at the start of the wake equals r_o itself, because 0.7 × (73/25)^(1/3) = 1.00. Read this way, the virtual origin sits (1/0.7)³ ≈ 2.92 reference radii upstream of where the wake starts. That is the rule I use below to scale the law to other vehicles.

**Validity statement.** This is on printed page 20, item (3), and is not quite as previously paraphrased:

- The law "assumes that the ambient air is perfectly still and that turbulent diffusion is caused solely by the passage of the spacecraft."
- From the momentum balance and Eq. (11), "at s = 10^5 ft, the mean particle velocity in the wake would decrease to 30 fps which is the same order as the velocity of natural turbulence."
- From that point, dispersion "would proceed as the 1/2-power of time (according to the Einstein's random walk theory) instead of 1/3-power."

So the cut-off is at a distance of 10^5 ft (30.5 km). At that point the wake's mean velocity has fallen to 30 ft/s (9.1 m/s). This cut-off was computed for Park's own vehicle and speed; it is not a universal constant.

**Conditions and source.** Printed page 17 gives the conditions: an altitude of 230,000 ft (70.1 km), flight speeds of 19,000 to 26,000 ft/s (5.8 to 7.9 km/s), and a lifting body at an angle of attack of about 25°. Printed page 19 gives the drag coefficient as 0.356 at 25° angle of attack.

Park assumes that "the wake of the lifting body is the same as that of a sphere" (printed page 17). He takes the growth law from Ref. 2, which is Schapker, R. L. and Camac, M., "NO2 chemiluminescent wake radiation," AIAA Journal 7(12), 2254–2261, December 1969 (reference list, printed page 24). I could not read Schapker and Camac online. Park's 0.7 is therefore second-hand even in Park, and its own range of validity lies in that paper.

Park also writes on printed page 13 that "the diameter of the wake grows approximately as 1/3 power of the distance behind the vehicle (Ref. 2)." Note that the coefficient carries no dependence on drag coefficient or Mach number.

## 4. Equations found for lower speeds

### 4.1 Taylor, Melcher and Washburn (1964): the drag-area correlation (best candidate, abstract only)

**Citation.** Taylor, R. L., Melcher, B. W. and Washburn, W. K. (1964), "Studies of the luminous hypersonic wake," AIAA Journal 2(10), 1731–1738. DOI: https://doi.org/10.2514/3.2656

**What was read.** Abstract only. I read it through the OpenAlex index (https://api.openalex.org/works/doi:10.2514/3.2656), which stores the abstract as optical-character-recognition (OCR) tokens.

**Equation.** The OCR tokens are "w/(CDA)112 = K[x/(CDA)ll2]l t", which I read as

  w / (C_D A)^(1/2) = K [ x / (C_D A)^(1/2) ]^t, with "K = 0.66" (verbatim token).

- **w** is the full width (diameter) of the luminous turbulent wake, in metres. The abstract says this agrees "well with the viscous core width data obtained with shadowgraph and schlieren techniques."
- **C_D** is the body's drag coefficient, a pure number.
- **A** is the body's frontal area, in square metres. So (C_D A)^(1/2) is a "drag length" in metres.
- **x** is the distance behind the body, in metres.
- **K** is a pure number.

**The exponent is uncertain.** The exponent t is garbled ("l t"). A secondary source supports one third. It is a ScienceDirect topic-page excerpt seen only as a search-engine snippet; the page itself refused the fetch. That excerpt attributes "w/√(C_D A) = 0.66 (x/√(C_D A))^(1/3)" to Slattery and Clay, with w the full width and A the base area. **Treat t = 1/3 as probable but unverified.**

If t = 1/3, the wake radius is w/2 = 0.33 (C_D A x)^(1/3).

**Validity stated in the abstract.**

- Bodies: spherical projectiles 0.22, 0.55 and 1.0 inch in diameter.
- Speeds: 12,000 to 20,500 ft/s (3.7 to 6.3 km/s).
- Pressures: free-stream pressures of 0.2 to 10.0 cm of mercury, which is about 270 to 13,000 pascals.
- Distance: wakes observed out to 1,200 body diameters.
- Pressure dependence: the abstract says "the luminous wake growth has a pressure dependence in the transition and laminar flow regimes, but not in the turbulent regime from 2.0 to approximately 76 cm Hg."
- Regime: turbulent, far wake, blunt bodies (spheres).

The abstract calls the result "the simple asymptotic turbulent wake growth correlation" and compares it with the theories of Lees and Hromas and of Lykoudis.

**Why it matters.** This is the only correlation found that does three things at once: it carries the drag coefficient, it was fitted to data below 5.8 km/s, and it states a numerical coefficient.

### 4.2 Lees (1964): the "universal relation" (abstract only, form garbled)

**Citation.** Lees, L. (1964), "Hypersonic wakes and trails," AIAA Journal 2(3), 417–428. DOI: https://doi.org/10.2514/3.2356

**What was read.** Abstract only, through the OpenAlex OCR tokens.

**Equation.** The abstract says "the growth of the turbulent far wake is described by the universal relation" followed by the garbled token "YTf/(Ci)A)^^^'[x/(CDA)^]^". It then says "regardless of body shape or Reynolds number." This is clearly the same drag-length form as section 4.1: the far-wake half-width Y_f divided by (C_D A)^(1/2), as a power of x divided by (C_D A)^(1/2). **Neither the coefficient nor the exponent can be read.**

**Transition criterion in the same abstract.** Laminar-to-turbulent transition in blunt-body wakes, in wind-tunnel and ballistic-range data, is correlated by a constant transition Reynolds number "± 15%" over Mach numbers of about 3.6 to 14.4. This Reynolds number is based on the distance x from the body to the transition point. Its value is garbled ("5 6 X 10<i"). The abstract also gives a minimum critical Reynolds number for a blunt body, based on free-stream conditions and body diameter, of "about 5 X 10" (exponent lost) at Mach 20, corresponding to 220,000 ft. **These numbers cannot be used until the paper is read.**

### 4.3 Lahaye and Doyon (1973): Mach-dependent coefficient (abstract only, exponent missing)

**Citation.** Lahaye, C. and Doyon, P. (1973), "Growth characteristics of turbulent wakes," Canadian Journal of Physics 51(20), 2138–2142. DOI: https://doi.org/10.1139/p73-279. This is work from the Defence Research Establishment Valcartier, formerly the Canadian Armament Research and Development Establishment.

**What was read.** Abstract only, through OpenAlex. The publisher's page refused the fetch.

**Equation.** The equation appears in the abstract as "[Formula: see text]." The surrounding text says three things:

- The wake width D_w and the axial distance X are "normalized by the drag diameter" (also a formula placeholder).
- "The factor C is a function of the Mach number M∞: C = 0.56 + 0.032 M∞."
- "The effect of P∞D (pressure times body diameter) on wake growth is negligible, while the Mach number has a strong effect." Cone and sphere wakes were both studied.

**What the coefficient implies.** The exponent and the exact definition of the drag diameter are not visible. The coefficient alone still says something. C is 0.72 at Mach 5 and 1.01 at Mach 14. So, at the same normalized distance, the wake at Mach 14 would be about 40 percent wider than at Mach 5.

This conflicts with Lykoudis's theoretical statement (section 4.4) that the far wake does not depend on Mach number. The conflict may be because the two refer to different distances: Lykoudis puts the compressibility effects within the first roughly 1,000 diameters. **This is the only source found that gives an explicit Mach-number dependence for the coefficient, and it is worth obtaining in full.**

### 4.4 Lykoudis (1963) and Lees and Hromas (1962): theory (search-engine snippets only)

**Lykoudis citation.** Lykoudis, P. S. (January 1963), "The growth of the hypersonic turbulent wake behind blunt and slender bodies," RAND RM-3270, DTIC AD0295144. URL: https://apps.dtic.mil/sti/tr/pdf/AD0295144.pdf (refused the fetch tool).

**What was read.** Secondary: search-engine summaries of the DTIC text only. According to those summaries:

- Compressibility affects the wake behind a blunt body only within about the first 1,000 diameters.
- Beyond that, the wake "grows with the 1/3 power rule, like a subsonic incompressible wake, with no dependence on the Mach number."
- Very far downstream, only the drag coefficient affects growth.

The equations themselves (numbered 37 and 38 in the snippet) were cut off.

**Lees and Hromas citation.** Lees, L. and Hromas, L. (1962), "Turbulent diffusion in the wake of a blunt-nosed body at hypersonic speeds," Journal of the Aerospace Sciences 29(8), 976–993. DOI: https://doi.org/10.2514/8.9671. The DTIC copy is Space Technology Laboratories Aerodynamics Department Report No. 50, AD0600632.

**What was read.** Abstract only, from search snippets and OpenAlex. The abstract says the predicted growth of the turbulent inner wake agrees well with shadowgraph widths behind spheres in ballistic ranges from 200 to 4,000 body diameters downstream. Its worked example is Mach 22 at 100,000 ft. No equation was visible.

### 4.5 Slattery and Clay (1961, 1962): measured exponent at 2.7 km/s (abstract only)

**Citations.**

- Slattery, R. E. and Clay, W. G. (1961), "Width of the turbulent trail behind a hypervelocity sphere," Physics of Fluids 4(10), 1199–1201. DOI: https://doi.org/10.1063/1.1706195
- Slattery, R. E. and Clay, W. G. (1962), "Measurement of turbulent transition, motion, statistics, and gross radial growth behind hypervelocity objects," Physics of Fluids 5(7), 849–855. DOI: https://doi.org/10.1063/1.1724457

Both are MIT Lincoln Laboratory work.

**What was read.** Abstracts only, through OpenAlex. The fetch tool paraphrased them rather than quoting.

**Findings of the 1961 paper.**

- Conditions: spheres at about 9,000 ft/s (2.7 km/s), in air at 40, 100 and 760 mm of mercury, followed to about 10^4 diameters. This pressure range spans a factor of about 19 in Reynolds number.
- Pressure dependence: growth depends little on pressure. At lower pressures the wake grows slightly more slowly.
- At one atmosphere, the width grows as the 1/2 power of distance.
- At 40 and 100 mm Hg, the width grows as the 1/2 power for about the first 50 diameters and as the 1/3 power after that.

A secondary review snippet (the ScienceDirect topic page again) states the low-pressure case the other way round: 1/3 power over the first 50 calibres. **The ordering needs checking in the paper.**

The 1962 paper covers 0.500-inch spheres at 8,500 ft/s and cones at 5,500 ft/s at reduced pressures. Its abstract gives no formula.

**Why it matters.** This is direct evidence that a sphere's wake at 2.7 km/s settles into cube-root growth after tens of diameters at low pressure. At one atmosphere, the square-root growth may instead reflect a different regime.

### 4.6 Other hypersonic and supersonic ballistic-range measurements (abstracts only, no coefficient)

**Wilson (1967).** Wilson, L. (1967), "Far wake behavior of hypersonic spheres," AIAA Journal 5(7), 1238–1244. DOI: https://doi.org/10.2514/3.4177. The abstract (through OpenAlex) covers ablating spheres at 19,000 to 23,000 ft/s (5.8 to 7.0 km/s). It says growth "approached classical one-third power variation with distance." This overlaps the bottom of Park's range and supports the exponent, not the coefficient.

**Murphy and Dickinson.** The secondary ScienceDirect excerpt (search snippet only) says that the one-third exponent for spheres was "confirmed by Murphy and Dickinson using a steel sphere at speeds of 2650–5300 ft/sec" (0.8 to 1.6 km/s). It also cites data of Dana and Short at 7,400 ft/s (2.3 km/s). These are the lowest-speed confirmations of the exponent found. I could not locate either primary source.

**Levensteins and Krumins (1967).** Levensteins, Z. J. and Krumins, M. V. (1967), "Aerodynamic characteristics of hypersonic wakes," AIAA Journal 5(9), 1596–1602. DOI: https://doi.org/10.2514/3.4256. The abstract describes a ballistic-range programme on how transition and wake structure depend on Mach number, Reynolds number and wall enthalpy. No growth law appears in the abstract.

**Redmond et al. (2026).** Redmond, J. J. et al. (2026), "Experimental analysis of chemically reacting wakes for ballistic spheres at Mach 10 flight conditions," AIAA paper 2026-5068. DOI: https://doi.org/10.2514/6.2026-5068. The abstract (from https://scholar.afit.edu/facpub/3584) covers 25 mm spheres at Mach 10.5 (about 3.4 km/s), at conditions equivalent to 20 km altitude.

- The wake became turbulent 1.16 to 1.22 diameters behind the sphere.
- Large-scale turbulence first appeared at 8.9 diameters.
- Data extend only to 53 diameters, and no width law is given.

This is useful only as a modern statement that, at this Reynolds number, the cube-root far-wake law cannot apply within the first ten or so diameters.

### 4.7 Classical incompressible axisymmetric wake (full text read; no numerical coefficient found)

**Dairay, Obligado and Vassilicos.** Dairay, T., Obligado, M. and Vassilicos, J. C., "Non-equilibrium scaling laws in axisymmetric turbulent wakes." I read the preprint in full, from https://www.imperial.ac.uk/media/imperial-college/research-centres-and-groups/turbulence-mixing-and-flow-control-group/main_production.pdf. It was later published in the Journal of Fluid Mechanics in 2015.

Its equation (1.2), page 1, gives the classical high-Reynolds-number law

  δ(x) / θ ∼ [ (x − x_0) / θ ]^(1/3),

with these symbols:

- δ is the wake width, in metres.
- θ is the momentum thickness, in metres.
- x_0 is a virtual origin, in metres.
- ∼ means "scales as": the paper gives no numerical coefficient.

The momentum thickness θ is fixed by the drag. From momentum conservation, θ² is proportional to C_D A; the numerical factor depends on the convention used. So this is the same physics as the drag-length form of section 4.1.

Its equation (1.4), page 2, gives the low-Reynolds-number regime, in which viscosity matters:

  δ / θ ∼ Re_G^(−1/2) [ (x − x_0) / θ ]^(1/2) (L_b / θ)^(1/2).

Here Re_G is a global Reynolds number set by the inflow, and L_b is the body's size. The paper's equation (1.9) shows that the same square-root growth can also occur at high Reynolds number in a "non-equilibrium" turbulent wake. It reports that this has been seen over large distances behind plates with irregular (multiscale) edges and, in de Stadler et al. (2014), behind a sphere.

**Skidmore, Philip and Monty (2018).** Skidmore, G. M., Philip, J. and Monty, J. P. (2018), "Evolution of the turbulent far wake of a sphere," 21st Australasian Fluid Mechanics Conference, Adelaide. URL: https://people.eng.unimelb.edu.au/imarusic/proceedings/21/Contribution_584_final.pdf. Full text read.

- Conditions: a 30 mm sphere towed through water at a Reynolds number of 13,000 (based on diameter), followed to 340 diameters.
- Results (page 3): the half-width followed the cube-root law only from about 20 to 50 diameters. It then switched to square-root growth from about 75 to 145 diameters, while the flow stayed turbulent.

The paper also quotes three thresholds from other work. I saw these only as quotations, so they are secondary here.

- From Johansson, George and Gourlay (2003): the switch begins when the wake Reynolds number Re_δ* (based on an integral wake width δ* and the centreline velocity defect) falls below about 500.
- From Uberoi and Freymuth (1970): the inertial range of the turbulence vanishes for Re_δ* at or below about 100.
- From Tennekes and Lumley (1972): relaminarisation (the wake becoming laminar again) needs Re_δ* of order one. They estimate where this happens as 2x / c_d ∼ (0.5 u∞ c_d / ν)³, with x in body diameters, c_d the drag coefficient, u∞ the speed and ν the kinematic viscosity. For their own sphere this gives about 4 × 10⁹ diameters.

**What this means for the small sphere.** At the Reynolds numbers of a small sphere in the mesosphere (estimated in section 6), the cube-root law will hold, if at all, over only a short stretch. The square-root law is the better default there.

**Uberoi and Freymuth (1970).** Uberoi, M. S. and Freymuth, P., "Turbulent energy balance and spectra of the axisymmetric wake," Physics of Fluids 13(9), 2205. This is secondary only, from a search snippet. Their sphere wake at a Reynolds number of about 8,600 becomes self-preserving from 50 diameters, with a virtual origin 12 diameters downstream. No coefficient was seen.

### 4.8 Laminar axisymmetric far wake (full text read): for low Reynolds number

**Citation.** Biau, D. (2011), "Exact self-similar solutions for axisymmetric wakes," Comptes Rendus Mécanique 339, 245–249. DOI: https://doi.org/10.1016/j.crme.2011.01.003. URL: https://comptes-rendus.academie-sciences.fr/mecanique/item/10.1016/j.crme.2011.01.003.pdf. Full text read.

**Equation.** The leading (slowest-decaying) mode, which the paper identifies as Schlichting's classical laminar solution (page 247), is

  u_0 ∝ x^(−1) exp(−η), with η = U_0 r² / (4 ν x).

- u_0 is the velocity defect, in metres per second.
- U_0 is the flight speed, in metres per second.
- r is the radius, in metres.
- x is the distance behind the body, in metres.
- ν is the kinematic viscosity of air, in square metres per second.

The radius at which the defect has fallen to 1/e (about 37 percent) of its centreline value is therefore

  r_e = 2 (ν x / U_0)^(1/2).

Because a point at distance x behind the body was passed by the body a time t = x / U_0 earlier, this is the same as r_e = 2 (ν t)^(1/2). In other words, the momentum deficit simply diffuses outward at the molecular rate. For a passive tracer, such as vapour or smoke, the molecular diffusivity replaces ν. The drag sets only the amplitude of the defect, not the width.

**Validity.** The solution assumes a small defect (well downstream, where the boundary-layer equations are linearised), a laminar wake and steady flow. A virtual origin x_0 may be added (page 247). It is also continuum theory, so it needs the wake to be much wider than the mean free path. Park's own late-time statement (square-root-of-time spreading by random walk) has the same form, but uses atmospheric turbulent diffusivity instead of molecular viscosity.

## 5. Comparison with Park at r_o = 1.4 m

These are my own calculations. The frontal area is A = π × 1.4² = 6.16 m².

**How the two laws were set up.**

- **Park.** I used the virtual-origin rule inferred in section 3: s' = x + 2.92 r_o, which is 4.1 m. At the distances below, the virtual origin changes the answer by less than 1 percent.
- **Taylor et al.** I used the radius w/2 = 0.33 (C_D A x)^(1/3), which assumes t = 1/3. The drag coefficient of the spacecraft-sized parent is not fixed here, so three values are shown: 0.356 (Park's lifting body), 1.0 and 1.5 (a blunt capsule).

The table compares predicted wake radii, in metres, at two distances behind the body.

| Distance behind body | Park, wake boundary radius r_b | Taylor et al., C_D = 0.356 | Taylor et al., C_D = 1.0 | Taylor et al., C_D = 1.5 |
|---|---|---|---|---|
| 1 km | 8.8 m | 4.3 m | 6.0 m | 6.9 m |
| 30 km | 27.2 m | 13.3 m | 18.8 m | 21.5 m |

The conclusion to draw from the table is that Park's law gives a wider wake than the ballistic-range correlation for any realistic drag coefficient. The two laws grow at the same cube-root rate, so the ratio between them is the same at every distance.

Written in the same form, Park's law is r_b = 0.478 (A x)^(1/3). Taylor's is 0.33 C_D^(1/3) (A x)^(1/3). They agree only if C_D is about 3.0. For a body with a drag coefficient of 1, Park's wake is 1.45 times wider than Taylor's. For Park's own lifting body (drag coefficient 0.356), Park's law gives a wake about 2.0 times wider than the drag-based correlation would. Over the 3.7 to 5.8 km/s overlap, the gap is therefore a factor of about 1.3 to 2 in radius, which is a factor of about 2 to 4 in wake cross-sectional area and hence in dilution.

**Why the two might differ.** Part of the gap may be definitional. Park's r_b is a "boundary" used for homogeneous mixing. Taylor's w is the luminous or shadowgraph viscous-core width, and a visible edge tends to sit inside the true mixing boundary. I could not check either definition in the primary sources.

**For comparison only: the blast radius.** The line-blast characteristic radius quoted by the earlier agent from ReVelle is R_0 = (E_0 / p_0)^(1/2). Here E_0 is the energy deposited per unit path length, which equals the drag force, and p_0 is the ambient pressure. I did not re-read that source. For this body with C_D = 1, R_0 works out to about 10 m at Mach 5 and 29 m at Mach 14. This is the radius of the shock-heated outer wake, not of the turbulent core. It is the same order as the turbulent radii above at 1 to 30 km, but it does not grow with distance.

## 6. When a turbulent-wake law stops applying (small aluminium sphere parent)

**Reynolds numbers.** These are my estimates, using approximate standard-atmosphere densities and viscosities and not taken from the sources.

The kinematic viscosity of air is about 0.017 m²/s at 50 km, 0.05 at 60 km, 0.2 at 70 km and 0.8 at 80 km. The Reynolds number U d / ν for some cases is therefore:

| Case | Altitude | Speed | Reynolds number |
|---|---|---|---|
| 5 mm sphere | 60 km | 1 km/s | about 100 |
| 5 mm sphere | 80 km | 1 km/s | about 7 |
| 400 mm sphere | 50 km | 1 km/s | about 24,000 |
| 400 mm sphere | 60 km | 300 m/s | about 2,300 |
| Spacecraft parent, 2.8 m across | 60 km | 4.3 km/s | about 2 × 10⁵ |
| Spacecraft parent, 2.8 m across | 50 km | 1.6 km/s | about 3 × 10⁵ |

The spacecraft parent's values are comparable with the ballistic-range sphere tests (for example, a 12.7 mm sphere at 2.7 km/s and 100 mm Hg is about 3 × 10⁵). So the ballistic-range correlations are being used within their own Reynolds-number range for that parent.

**What the literature says.**

1. **Transition to a turbulent wake.** The low-speed sphere wake is steady, and therefore laminar, up to a Reynolds number of a few hundred. The abstract of Bhowmick et al. (2020), *Symmetry* 12(9), 1498 (https://arxiv.org/abs/2008.03074, abstract only) describes the steady regime as covering 25 to 285, with the oblique steady wake above about 225. At hypersonic speed, Lees (1964, abstract) gives a transition criterion based on a distance Reynolds number that holds for Mach numbers of about 3.6 to 14.4, but its value is unreadable online. Slattery and Clay and Taylor et al. both say the wake depends on pressure (that is, on Reynolds number) in the laminar and transitional regimes and not in the turbulent one. **So a 5 to 50 mm sphere at 60 to 90 km will usually have a laminar wake, and the laminar law of section 4.8 applies.**
2. **Late-time departure even when turbulent.** Skidmore et al. found that a sphere wake at a Reynolds number of 13,000 left cube-root growth after only about 50 diameters, for square-root growth. This used Johansson et al.'s threshold of a wake Reynolds number of about 500 (section 4.7). **So for spheres up to 400 mm at these altitudes, the turbulent cube-root law should be expected to hold only over tens of diameters, if at all.**
3. **Rarefied flow.** I found no wake-growth law for transitional or free-molecular flow, meaning flow where the mean free path is comparable to or larger than the body. The nearest online material is direct-simulation Monte Carlo studies of near wakes: computer simulations that track representative gas molecules. One is Moss, Wilmoth and Dogra (1993), "DSMC simulations of Mach 20 nitrogen flows about a 70 degree blunted cone and its wake," NASA TM-107762 (https://ntrs.nasa.gov/api/citations/19940010818, abstract only), which covers free-stream Knudsen numbers (mean free path divided by body size) of 0.001 to 0.03. A search-engine summary of a cylinder study, which I did not read, reports that rarefaction shrinks and eventually removes the recirculating bubble behind the body. For a 5 mm sphere at 90 km, the mean free path is around 30 mm, so the Knudsen number is about 6 and the flow is nearly free-molecular. The body then leaves a momentum and density deficit but no continuum wake structure near it. Farther downstream, once the disturbed region is many mean free paths wide, the laminar diffusion law of section 4.8 (r of order 2 (ν t)^(1/2)) is the physically consistent continuation, but **no source was found that states this for a rarefied body.**

## 7. Recommendation and what remains unknown

**For 1.6 to 5.8 km/s with the spacecraft-sized parent.** Use the drag-length correlation:

  wake radius = 0.33 (C_D A x)^(1/3)

This comes from Taylor, Melcher and Washburn (1964), using K = 0.66 and the probable exponent of one third. It was fitted to sphere wakes at 3.7 to 6.3 km/s, overlaps Park at the top end, and through C_D it carries the body's changing drag as it slows. Wilson (1967) and the Murphy and Dickinson data quoted by a secondary source support the cube-root exponent down to about 0.8 km/s. Because Park's law is 1.45 times wider for C_D = 1, switching laws at 5.8 km/s would make the wake radius jump. There are three ways to deal with that:

- **(a) Keep Park above 5.8 km/s and Taylor below it.** This gives a discontinuity at 5.8 km/s.
- **(b) Use Taylor everywhere.** This is internally consistent, and it is my recommendation, because the 0.7 in Park comes second-hand from Schapker and Camac and has no drag dependence.
- **(c) Scale Taylor's coefficient to match Park at 5.8 km/s.** This keeps continuity with the earlier model but carries Park's unexplained width into lower speeds.

**For a subsonic, transonic or rarefied small-sphere parent.** Use the laminar diffusion law r_e = 2 (ν x / U)^(1/2) (Biau 2011, Schlichting's solution) whenever the Reynolds number is below a few hundred, which covers most of the 5 to 400 mm cases at 60 to 90 km. Here ν is the molecular kinematic viscosity at the local altitude. For the larger spheres at 50 to 60 km, where the Reynolds number is in the thousands, use the turbulent drag-length law only over the first few tens of diameters. After that, switch to square-root growth (section 4.7).

**What remains unknown.**

1. **The exponent in Taylor et al.** It is unreadable in the abstract, and the one third is from a secondary snippet. The definition of w (luminous or viscous core, full width) is also unchecked.
2. **Lahaye and Doyon's formula and their definition of the drag diameter.** These are the only Mach-dependent coefficient found. The coefficient implies about 40 percent wider wakes at Mach 14 than at Mach 5.
3. **The primary sources behind Park's 0.7 and behind the drag-based theory.** These are Schapker and Camac 1969, Lees and Hromas 1962, and Lykoudis 1963. All are paywalled or behind DTIC and RAND servers that refused the fetch tool. A library copy of Lees (1964) and Lykoudis (1966, AIAA Journal 4, 577) would settle points 1 and 2 most cheaply.
4. **Whether any law holds for a rarefied body.** None was found.
5. **Whether wakes at mesospheric Reynolds numbers ever reach cube-root growth.** All the hypersonic correlations come from ballistic ranges at 270 to 100,000 pascals, far above the 0.2 to 80 pascals at 50 to 90 km. They match only in Reynolds number, for the large parent.

Everything in sections 4.1 to 4.6, apart from Park, rests on a single reading of an abstract or a search snippet, never on a full paper.
