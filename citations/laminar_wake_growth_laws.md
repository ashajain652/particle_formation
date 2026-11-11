# Update 2026-10-29: full texts of Feldman (1961), Lykoudis (1966), Pallone et al. (1964)

The three papers were read from the PDFs the user saved in `Literature Review/Wake Growth/`. Every equation quoted below was read from the page image, not from the garbled text layer, and the printed journal page is given. Statements marked "inferred" are my own algebra or judgement and are not on the page.

## Summary

Only Feldman (1961) gives something that can be turned into a closed-form laminar trail-width law. Lykoudis (1966) gives an approximate radial enthalpy profile and a regime-validity condition, plus Lees's transition picture with a single number. Pallone, Erdos and Eckerman (1964) give only numerical (strip-integral) solutions and a graphical transition correlation; they contain no closed-form width law. The usable law (inferred from Feldman's equations) is that the square of the Gaussian 1/e radius of the hot core grows linearly with distance, with an effective diffusivity built from the viscosity at the hot axis temperature and the ambient density. Applied to our parent body, this gives a laminar radius of roughly 1.6 to 2.7 metres at 1 kilometre and 4.7 to 12 metres at 30 kilometres behind the body at 80 kilometres altitude, and roughly 2.3 to 5.4 metres and 10 to 28 metres at 90 kilometres. However, Lees's transition rule, as quoted by Lykoudis, suggests the 80-kilometre wake would already be turbulent within tens of metres, so the laminar law most plausibly applies only at 90 kilometres and above, and even there not certainly.

## 1. Feldman (1961), J. Aerospace Sciences 28:433-449

**What the paper does.** Feldman treats the trail behind an axisymmetric, non-ablating blunt body (a hemisphere-cylinder) in equilibrium air, at 15,000 to 35,000 feet per second and 100,000 to 250,000 feet altitude (about 30 to 76 kilometres; page 433 abstract). He splits the trail into an expansion-controlled region, where the hot gas cools by expanding back to ambient pressure, and a conduction-controlled laminar trail downstream of the point where the axis pressure first reaches ambient pressure (the "p goes to p-infinity" point). Our altitudes, 80 and 90 kilometres, lie above his range.

**The governing equation (page 442, Eq. 6).** Radial heat conduction in a column moving at average axial speed f times the free-stream speed:

(1/r\*) d/dr\* [ (1/Pr)(mu/mu_inf) r\* d(h/h_inf)/dr\* ] = (rho/rho_inf) f d(h/h_inf)/d(x\*/Re)

Here r\* is radius divided by the nose radius r_n; x\* is axial distance divided by r_n, measured from the p-to-ambient point; Re = rho_inf u_inf r_n / mu_inf is the free-stream Reynolds number on the nose radius; h is static enthalpy and h_inf its ambient value; Pr = mu c_p / k is the Prandtl number (taken as about 0.7 for equilibrium air, page 443); mu is viscosity; and f is the axial velocity of the hot column divided by the flight velocity, "in the neighborhood of 0.75" depending weakly on velocity and altitude (pages 442-443). Feldman states (page 442) that this is the energy equation with radial velocity and radial gradient of axial velocity set to zero, i.e. pure conduction in a slug of uniform speed.

**Howarth transform (page 442, Eq. 8).** R^2 = 2 integral from 0 to r\* of (rho/rho_inf) r\* dr\*, and X = x\*/Re. R is a density-weighted radius in units of r_n: in the hot, low-density core the physical radius is larger than R.

**Gaussian solution (page 442, Eqs. 10-11; page 443, Eq. 14).** The enthalpy excess is assumed Gaussian, (h/h_inf) - 1 = alpha(X) exp(-R^2/beta(X)), and energy conservation gives alpha(X) beta(X) = alpha(0) beta(0). Taking the axis viscosity ratio as a slowly varying constant mu-bar/mu_inf (written mu-bar/mu_0 on the page; mu-bar means viscosity on the axis) gives

alpha(X)/alpha(0) = [ 1 + (4/(Pr f beta(0))) (mu-bar(X,0)/mu_inf) X ]^(-1).

Viscosity follows Sutherland's formula (page 443, Eq. 13): mu/mu_0 = (T/273)^(3/2) x 386/(113 + T), with T in kelvin and mu_0 = 1.72e-4 g/(cm s).

**Width law (inferred from Eqs. 8, 11 and 14).** Because beta = alpha(0) beta(0)/alpha, beta(X) = beta(0) + 4 (mu-bar/mu_inf) X / (Pr f). Converting back to dimensional variables (multiply by r_n^2 and substitute X = x mu_inf/(rho_inf u_inf r_n^2)) gives, for the square of the 1/e radius of the enthalpy (temperature) excess in the Howarth coordinate,

R_e^2(x) = R_e^2(0) + 4 kappa_eff x / (f u_inf), with kappa_eff = mu_axis / (rho_inf Pr).

So the effective diffusivity is the viscosity at the hot axis temperature divided by the ambient density and the Prandtl number. The nose radius cancels. This is the width of the temperature (enthalpy) profile, not of the velocity profile. The physical radius exceeds R_e wherever the gas is hotter than ambient; far downstream, once the core has cooled towards ambient, the two coincide. Feldman's initial width at the p-to-ambient point is beta(0) of about 1.0 to 1.4 (R_e of about one nose radius; page 442, Fig. 18) for his 60,000-foot, 17,500-foot-per-second reference case, with peak enthalpy alpha(0) of about 21.9.

**Accuracy (page 441).** The integral (Gaussian) solution with variable properties differed from a fuller solution including convection and viscous terms by "less than 5 percent in temperature at a given axial station". Section 3.3 (page 443, Eqs. 15-19) adds a momentum integral with a Gaussian velocity defect, 1 - u/u_inf = a(X) exp[-R^2/b(X)], where b(X) = D/[a(X)(1 - a(X)/2)] and D is essentially half the drag coefficient based on the maximum cross-section; with constant axial velocity it reduces to Eqs. 11 and 14 with f replaced by (1 - a). That gives a velocity-width law of the same family, but I did not read a closed solution for b(X) beyond Eq. 17.

**Constant-property bounds (page 441, Section 3.1).** Feldman explicitly says that using the axis diffusivity cools the trail faster than using the free-stream diffusivity, and that the two bracket the answer (Fig. 16 curves a and b); the variable-property result lies between them.

**Laminar versus turbulent (page 435).** He notes that a turbulent trail cools faster by several orders of magnitude and gives a shorter observable trail, but gives no transition criterion. He also notes (page 441) that at 250,000 feet and 25,000 feet per second, conduction and expansion must be computed together when the body is smaller than 3 or 30 centimetres; that is, for bodies our size (metres) the separation into expansion then conduction should be valid. No full-scale worked trail widths at 60-100 kilometres were found; his worked case is 60,000 feet.

## 2. Lykoudis (1966), AIAA Journal 4:577-590 (review)

**Closed-form laminar results (page 579, Eqs. 5-9), quoted from Lykoudis's own Ref. 3 (closed-form approximate solutions that "compared well with the numerical computations of Feldman").** At the p-to-ambient point behind a blunt body in equilibrium:

- the axis velocity ratio is f = (1 - M_inf^(-1/3))^(1/2) (Eq. 5), about 0.8 for free-stream Mach number M_inf between 10 and 30;
- the axis enthalpy is h(0)/h_inf of about (1 - f^2)(gamma - 1) M_inf^2/2, about M_inf^2/15 (Eq. 6), with gamma the ratio of specific heats;
- the radial enthalpy profile is h(R)/h(0) of about 1/[1 + 1.2 R^2]^0.8 (Eq. 7), where R is "very nearly the compressible Howarth variable"; the review does not state the normalising length of R here (from Eq. 9 it is presumably the sphere radius; inferred);
- the local Mach number there is about 3, independent of free-stream Mach number (Eq. 8);
- the expansion and conduction regions can be treated separately if M^2/Re is at most 25 (Eq. 9), with Re based on sphere radius and free-stream conditions.

The review does not reproduce a growth law for the downstream conduction region beyond referring to Feldman. It notes (page 579) that these laminar equilibrium models are not consistent with reality at both ends: at low altitude the core is turbulent, and at high altitude, where it may be laminar, it is not in equilibrium.

**Transition (page 582).** Lykoudis quotes Lees's tentative picture: below a minimum critical Reynolds number the wake is laminar; above it, transition first appears about 40 to 50 diameters behind a blunt body and moves forward as ambient pressure increases, "maintaining a constant value of Reynolds number based on the local distance x and conditions prevailing at the turbulent front of 5.6 x 10^4, independent of body diameters". For sharp slender bodies the value is about four times higher. He also quotes that the local wake Reynolds number scales as Re^(-1/4) p_inf d, so stability depends mainly on the product of ambient pressure and body diameter. The minimum critical Reynolds number itself was not given as a number on the page I read. He notes (footnote, page 582) contradictory ballistic-range data with transition beyond 1000 diameters at Mach 17.

## 3. Pallone, Erdos and Eckerman (1964), AIAA Journal 2:855-863

**Method (pages 855-857).** Boundary-layer equations for the laminar far wake (Eqs. 1-3), solved by an N-strip integral method in a Howarth-type transformed coordinate eta (page 856), marched numerically on an IBM 7094. There is no closed-form width law. Variable properties enter through l = rho mu/(rho_1 mu_1) (the Chapman-Rubesin parameter relative to the inviscid edge state, subscript 1), evaluated by Sutherland's law with S = 120 K at low temperature and, at high temperature, by a fit l = [A (h_k/h_r)^(-1/2) - B (h_k/h_r)^(-1)]/[A (h_1/h)^(-1/2) - B (h_1/h)^(-1)] with A = 3.03, B = 2.03 and a reference enthalpy of 4.5e6 square feet per second squared (page 857). So they do not assume rho mu constant.

**Full-scale example (pages 859-860, Figs. 7-9).** A 10-foot, 12-degree half-angle cone at 22,000 feet per second, 200,000 to 100,000 feet, initial enthalpy ratio 0.3, characteristic length L = 20 feet; profiles at 150,000 feet (about 46 kilometres). Read from Figs. 7 and 8: the velocity and temperature profiles extend to y/L of about 0.08 to 0.09 (about 1.6 to 1.8 feet) and this outer edge stays nearly fixed from X/L = 1 to 60 (20 to 1200 feet) while the core profile fills out. My reading (inferred) is that the physical width hardly grows over 1200 feet because the core cools and contracts while the Howarth width grows. This is a slender body, not a blunt one, and an altitude well below ours.

**Transition correlation (page 861, Fig. 13).** Transition data are correlated as a Reynolds number Re_delta = (rho_1 u_1 delta/mu_1)[1 - (u_c/u_1)] against a Mach number M_1 = (u_1/a_1)[1 - (u_c/u_1)], where delta is the wake diameter at transition, u_c the axis velocity, and subscript 1 the local inviscid edge. Transition Re_delta rises from about 10^3 to 10^4 at M_1 near 1-2 to about 10^5 at M_1 near 6-8 (read from the figure; scatter is large and temperature effects were not accounted for). Using it needs the computed axis velocity defect and wake diameter, so it is not a free-stream quick test. Their prediction (Fig. 14, page 861) is for slender cones only.

## 4. Application to our parent body

**Assumptions.** Equal-area diameter 2.81 metres; initial 1/e radius of the hot core R_e(0) taken as 1.4 metres (the body radius; Feldman's initial width is about one nose radius, and our object has no defined nose radius, so this is an assumption); f = 0.8 (Lykoudis Eq. 5); Pr = 0.71; Sutherland viscosity; distance measured from the p-to-ambient point, which I take as close to the body. The axis temperature is the controlling unknown. Lykoudis Eq. 6 gives an initial axis enthalpy of about 48 times ambient at 80 kilometres and 54 times at 90 kilometres (Mach 27 and 28.5), which in equilibrium air corresponds to several thousand kelvin (inferred). Since the core cools downstream, I give three cases for the axis viscosity: ambient temperature (lower bound, which is the Schlichting-type estimate with a Gaussian 1/e radius and the f factor), 2000 K, and 5000 K (upper bound). The radius quoted is sqrt(R_e(0)^2 + 4 kappa_eff x/(f u_inf)), a Howarth (density-weighted) radius of the temperature excess.

At 80 kilometres (rho = 1.85e-5 kilograms per cubic metre, V = 7.6 kilometres per second), the effective diffusivity kappa_eff is 1.0, 4.7 and 7.7 square metres per second for the three viscosity cases. The radius is 1.6, 2.3 and 2.7 metres at 1 kilometre, and 4.7, 9.7 and 12.4 metres at 30 kilometres.

At 90 kilometres (rho = 3.42e-6 kilograms per cubic metre, V = 7.8 kilometres per second), kappa_eff is 5.2, 25 and 42 square metres per second. The radius is 2.3, 4.3 and 5.4 metres at 1 kilometre, and 10, 22 and 28 metres at 30 kilometres.

Comparison: Taylor's turbulent law gives 5.8 metres at 1 kilometre and 18 metres at 30 kilometres. Schlichting r = 2 sqrt(nu x/U) at ambient nu gives 0.6 and 3.4 metres at 80 kilometres and 1.4 and 7.5 metres at 90 kilometres. The Feldman-type law with ambient viscosity is only slightly larger than Schlichting (because of the 1/0.71 Prandtl factor, the 1/0.8 slow-column factor and the initial radius). With hot-axis viscosity it is 3 to 7 times the Schlichting value in diffusivity, about 2 to 3 times in radius; at 90 kilometres with a hot core, the laminar radius at 30 kilometres (22-28 metres) exceeds Taylor's turbulent 18 metres. The conclusion is that at 90 kilometres a laminar wake is not necessarily narrower than Taylor's turbulent estimate.

**Regime check (inferred from Lees's rule, Lykoudis page 582).** Setting rho u x/mu = 5.6e4 with free-stream density and velocity: with ambient viscosity, x is about 5 metres at 80 kilometres and 26 metres at 90 kilometres; with viscosity at 2000 K, about 25 and 130 metres. All are far shorter than 1 kilometre, and the free-stream Reynolds number on diameter is about 30,000 at 80 kilometres and 6,000 at 90 kilometres. This rule therefore predicts turbulence well within the first kilometre at both altitudes, unless the Reynolds number is below Lees's minimum critical value, which the review does not quantify. Lykoudis's M^2/Re condition (Eq. 9, Re on radius 1.4 metres) gives 0.05 at 80 kilometres and 0.27 at 90 kilometres, both well under 25, so the separate expansion-then-conduction treatment would be valid if the wake were laminar.

## 5. What remains unknown

The axis temperature history of a large blunt tumbling body at 80-90 kilometres was not computed; it sets kappa_eff within a factor of about 8. None of the three papers covers altitudes above 76 kilometres (Feldman) or blunt bodies at full scale above 46 kilometres (Pallone), and Lykoudis notes that equilibrium does not hold where the flow is laminar. Lees's minimum critical Reynolds number, which decides whether the 90-kilometre wake ever transitions, was not found on the pages read. Pallone's correlation needs the local velocity defect, which we do not have. Each number above is a single analytical estimate, not checked against data at our conditions.

---

# Laminar wake growth laws for blunt bodies at hypersonic speed

Search date: 2026-10-29. This was a short search. Its main finding is negative: no published laminar hypersonic wake-width law with a fitted numerical constant was found and read in full. The 1960s primary sources that would contain one are technical reports. These are the RAND memoranda by Lykoudis, Feldman's 1961 paper, and the American Rocket Society (ARS) preprints by Fay. The Defense Technical Information Center (DTIC) blocked automated fetching (it returned HTTP 403, access forbidden), so none of these was read.

## 1. Exact citations recovered from the reference lists

These citations were read from the reference lists of the two papers already held locally, which were opened with pymupdf.

The laminar diffusion relation in Taylor, Melcher and Washburn (1964) is attributed to their reference 24: Fay, J. A., "A review of laboratory experiments on hypersonic wakes," ARS Reprint 2660-62 (1962). It was not found online.

Taylor's reference 9 is Fay, J. A. and Goldburg, A., "The unsteady hypersonic wake behind spheres," ARS Preprint 2676-62 (1962). It was not read.

Taylor's reference 3 is Lykoudis, P. S., "The growth of the hypersonic wake behind blunt and slender bodies," RAND RM-3270-PR (January 1963). The RAND catalogue lists this report as the *turbulent* wake (https://www.rand.org/pubs/research_memoranda/RM3270.html). Only the catalogue entry was seen.

Lees (1964), reference 3, is Feldman, S., "On trails of axisymmetric hypersonic blunt bodies flying through the atmosphere," Journal of the Aerospace Sciences 28, 433-448 (1961). This is the most likely source of a laminar "trail" radius law. It was not found online.

Lees (1964), reference 4, is Vaglio-Laurin, R. and Bloom, M. H., "Chemical effects in external hypersonic flows," ARS Preprint 1976-61 (1961), later published in *Hypersonic Flow Research* (Academic Press, 1962). It was not read.

Further leads found by search, none of them read:
- Lykoudis, "Length of the laminar hypersonic wake during the ballistic re-entry of a blunt object," RAND P-2781 (1963), at https://www.rand.org/pubs/papers/P2781.html.
- Lykoudis, "A review of hypersonic wake studies," RAND RM-4493-ARPA (1965), and AIAA Journal 4, 577 (1966). This is the best single review to obtain next.
- DTIC report AD0600632 (https://apps.dtic.mil/sti/pdfs/AD0600632.pdf), whose figure list includes wake-width comparisons. The fetch was blocked.
- Steiger and Bloom's laminar axisymmetric non-equilibrium wake reports from the General Applied Science Laboratories (GASL), 1961.

## 2. The only law with a stated constant: the incompressible laminar far wake

The only law with a stated constant that is available is the classical incompressible axisymmetric far wake. It comes from Schlichting's *Boundary-Layer Theory*, which was read previously in secondary form (Biau 2011), not in this search:

r ≈ 2 (ν x / U)^(1/2)

In this law, r is the wake radius in metres, ν is the kinematic viscosity (the dynamic viscosity divided by density) in square metres per second, x is the distance behind the body in metres, and U is the flight speed in metres per second. The law is valid only in the far wake, where the velocity defect is small compared with U. It assumes constant density and constant viscosity, and it uses no hypersonic or compressibility correction. The constant 2 depends on how the edge of the wake is defined. It also gives no initial radius, so it should be added to an initial width of about one body radius.

Taylor's Eq. (1) has the same form, w proportional to (ν x / V)^(1/2). Its constant was fitted to a single datum and was not published separately; see the existing notes.

## 3. How to evaluate ν: what is known

Taylor's Eq. (2) takes ν = μ/ρ to be proportional to R T^1.76 / p, where μ is dynamic viscosity, ρ is density, R is the gas constant, T is temperature and p is pressure. This makes ν rise about as fast as T^1.76 at constant pressure. Taylor wrote T in this relation as the temperature of the wake gas, not the ambient temperature. That reading comes from the earlier full-text read, and the exact definition of T should be rechecked in that paper.

No source was found that states what reference temperature to use for a laminar hypersonic far wake. The theory papers (Kubota, Lykoudis and Feldman) handle this with a density-weighted (Howarth-Dorodnitsyn) coordinate transform. In that approach the effective ν is roughly ρ_w μ_w / ρ∞², where the subscript w means the wake and ∞ means the ambient freestream. This is my own inference from the standard method and was not read in any of those papers.

## 4. Numbers for the parent spacecraft

The calculation uses r = 2 (ν x / U)^(1/2) with Sutherland viscosity and the US Standard Atmosphere 1976 (US76) values given in the brief. The wake-temperature cases scale the ambient ν by (T_wake / T∞)^1.76 at ambient pressure, following Taylor's Eq. (2). No initial radius is added, so the predicted widths should be increased by about 1.4 m, the body radius.

| Case | ν (m²/s) | r at 1 km (m) | r at 30 km (m) |
|---|---|---|---|
| 80 km, ν at ambient 199 K | 0.71 | 0.6 | 3.4 |
| 80 km, ν at 2000 K | 42 | 4.7 | 26 |
| 80 km, ν at 5000 K | 209 | 10.5 | 57 |
| 90 km, ν at ambient 187 K | 3.7 | 1.4 | 7.5 |
| 90 km, ν at 2000 K | 238 | 11 | 61 |
| 90 km, ν at 5000 K | 1190 | 25 | 135 |

The table shows that the choice of temperature changes the predicted radius by a factor of 10 to 20. With ν taken at ambient temperature, the laminar wake is narrower than Taylor's turbulent radius of 5.8 m at 1 km and 18 m at 30 km. This is the physically expected ordering. With ν taken at a hot-wake temperature of 2000 K or more, the laminar wake comes out *wider* than the turbulent one. That is implausible, because turbulence should spread a wake faster than molecular viscosity, not slower. It suggests that a hot ν applied over the full 30 km overestimates growth, since the wake cools and ν falls with distance. The true answer probably lies between the ambient row and the 2000 K row, and it should probably be capped at the turbulent value.

## 5. Recommendation and unknowns

I recommend using r = r0 + 2 (ν x / U)^(1/2), where r0 is the initial radius, taken as about the body radius. ν should be evaluated at ambient conditions as the lower bound and at a decaying wake temperature as the upper bound. The predicted radius should not be allowed to exceed Taylor's turbulent law. This recommendation rests on an incompressible coefficient; no hypersonic laminar coefficient has been read.

Three things remain unknown:
- A published hypersonic laminar width coefficient. The best sources to obtain next are Feldman (1961) and the Lykoudis 1966 review.
- A justified reference temperature, and how the wake temperature decays with distance.
- How the law applies to the rarefied small-sphere cases. For spheres at a Reynolds number of about 20, continuum wake theory is doubtful, and Direct Simulation Monte Carlo (DSMC) studies were not searched.

Read status: only reference lists from the local full texts were read. Everything else in this file is search snippets or catalogue entries.
