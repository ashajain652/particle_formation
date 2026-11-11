# Scale-dependent turbulent dispersion of a small plume (size at or below the smallest eddies), and how to represent it in the reentry vapour-plume model

Provenance convention used throughout. "Cited" means a claim a source in this session's search results supported (links given). "[Textbook knowledge, not re-verified this session]" marks standard results (Batchelor 1952, Monin and Yaglom, Pope, Csanady) that I state from memory because I could not open the full text; the report writer should treat them as needing a page-level citation before publication. "[Own calculation]" marks arithmetic I did myself (Python, unit-coefficient scalings, N = 0.02 rad/s assumed). Notation: D = molecular diffusivity of the vapour (m^2/s); nu = kinematic viscosity (m^2/s); Sc = nu/D (about 1 here); epsilon = turbulent kinetic energy dissipation rate (W/kg = m^2/s^3; 1 mW/kg = 1e-3); N = buoyancy (Brunt-Vaisala) frequency (rad/s); l = plume size (for a pair of fluid points, their separation); K(l) = relative (scale-dependent) diffusivity (m^2/s); eta = Kolmogorov length = (nu^3/epsilon)^(1/4); tau_eta = Kolmogorov time = (nu/epsilon)^(1/2); L_O = Ozmidov length = (epsilon/N^3)^(1/2); K_eddy = large-scale eddy diffusivity.

## 1. Relative dispersion theory: Richardson, Batchelor, Kraichnan, and how a blob smaller than the Kolmogorov scale behaves

### Takeaway
Relative dispersion has three regimes set by pair separation l: for l below eta the separation grows from the strain rate (about 1/tau_eta, exponential in time in the Batchelor/Kraichnan picture), for eta much less than l much less than the outer scale it follows Richardson-Obukhov (mean-square separation = g epsilon t^3, equivalent to K proportional to epsilon^(1/3) l^(4/3)), and beyond the outer scale it is ordinary Fickian diffusion with K_eddy. For Sc near 1 the Batchelor scale (where molecular diffusion stops strain-driven thinning) equals eta, so there is essentially no separate sub-Kolmogorov scalar range, and below eta molecular diffusion dominates.

### Cited Findings
- Sawford's review (Annu. Rev. Fluid Mech. 33, 289-317, 2001) begins from Kolmogorov similarity scaling for relative dispersion, focuses on stochastic models of pair dispersion, stresses the role of the two-particle acceleration covariance, and uses homogeneous scalar fields to link relative dispersion to turbulent mixing — [Scilit record](https://scilit.com/publications/aad92af5ee329345331a5bfeebbf7e77); [secondary description via search, Salazar and Collins context](https://arxiv.org/pdf/physics/0503169)
- Salazar and Collins (Annu. Rev. Fluid Mech. 41, 405-432, 2009) review experiments, direct numerical simulation (DNS) and theory of two-particle dispersion and complement Sawford's review; a later paper reports they concluded the Richardson t^3 law remains unchallenged although no experiment had clearly confirmed it. I could not access the review text; this is from secondary summaries — [Annual Reviews PDF link (not readable here)](https://annualreviews.org/doi/pdf/10.1146/annurev.fluid.40.111406.102224); [search summary source](https://arxiv.org/pdf/physics/0503169)
- Richardson (1926) first set K(R) = alpha R^(4/3) from atmospheric data on empirical grounds; Obukhov-style dimensional reasoning (Lagrangian velocity differences equal Eulerian ones at the same separation, d(R^2)/dt of order R times the velocity increment at R) gives the same R^(4/3) scaling; the law holds only for R much greater than the dissipation length — [Richardson law and intermittency, arXiv nlin/0107061](https://arxiv.org/pdf/nlin/0107061); [Pair dispersion in synthetic turbulence, arXiv chao-dyn/9807026](https://arxiv.org/pdf/chao-dyn/9807026)
- DNS gives a Richardson-Obukhov constant g of about 0.52 plus or minus 0.05 in <l^2> = g epsilon t^3 after a transient, with a correction of the form g epsilon t^3 (1 + C t0/t) where t0 = (r0^2/epsilon)^(1/3) is the Batchelor time (time over which a pair remembers its initial separation r0) — [Geometry and violent events in turbulent pair dispersion](https://arxiv.org/pdf/1209.2296); [Experimental pair dispersion paper summarising t0](https://arxiv.org/pdf/1205.4140)
- A convective-flow experiment obtained g = 0.10 plus or minus 0.07, much smaller than other flows, and another experimental paper saw no clean Richardson t^3 range at modest ratios of separation to integral scale; a clean t^3 regime needs a wide inertial range — [Experimental pair dispersion, convective flows](https://arxiv.org/pdf/1205.4140); contrast the DNS value above
- In the dissipation range, DNS shows exponential separation only after very long times (beyond about 200 Kolmogorov times), attributed to strain-rate fluctuations; an experiment reported Batchelor's prediction fulfilled over more than two decades in time at high turbulence level and saw no Richardson-Obukhov range up to Taylor-scale Reynolds number 815 — [Pair dispersion in turbulence, arXiv physics/0503169](https://arxiv.org/pdf/physics/0503169); [arXiv 1209.2296](https://arxiv.org/pdf/1209.2296)
- Batchelor's passive-scalar picture for Sc much greater than 1: scalar variance spectrum goes as k^-1 in the viscous-convective range between the viscous cutoff and the scalar diffusive cutoff; Kraichnan argued a steady transfer of scalar variance to high wavenumbers requires E(k) proportional to k^-1; a blob is stretched in one direction and compressed along the most negative strain direction — [Batchelor spectrum of passive scalar turbulence](https://ar5iv.labs.arxiv.org/html/1911.11014); [Passive scalar mixing in the Batchelor regime](https://arxiv.org/pdf/1511.08206)
- Batchelor scale = Kolmogorov scale divided by sqrt(Sc); the experimental cutoff form is E(k) = C k^-1 exp(-k sqrt(kappa/gamma_eff)) with kappa scalar diffusivity and gamma_eff an effective stretching rate; this validity range requires Sc much greater than 1 — [High Schmidt-number turbulent advection](https://arxiv.org/pdf/2112.13115); [arXiv 1911.11014](https://arxiv.org/pdf/1911.11014.pdf)
- The k^-1 law is largely a steady-state result; for decaying scalars the late-time spectrum is of form k^(1/2) — [arXiv 1511.08206](https://arxiv.org/pdf/1511.08206)

### Inferences
Equations (all symbols defined in the header; classical forms are [Textbook knowledge, not re-verified this session]: Batchelor 1952 Proc. Camb. Phil. Soc. 48, 345; Batchelor 1959 J. Fluid Mech. 5, 113; Kraichnan 1974 J. Fluid Mech. 64, 737; Monin and Yaglom vol. 2; Pope 2000 ch. 12 and ch. 6 for eta and the Kolmogorov scales):

- Viscous (Batchelor) range, l less than eta: the velocity difference across l is linear in l, du of order l/tau_eta, so a pair separates as l(t) = l0 exp(lambda t), where lambda is of order 1/tau_eta = (epsilon/nu)^(1/2) times an order-one number (the Lyapunov exponent). Effective relative diffusivity: K(l) of order l^2/tau_eta = (epsilon/nu) l^2.
- Inertial range, eta much less than l much less than the outer scale: <l^2> = g epsilon t^3 with g about 0.5 (cited above); eliminating t gives K(l) = beta epsilon^(1/3) l^(4/3) with beta of order 0.4 to 0.6 depending on the dimensionality convention (my derivation: d<l^2>/dt = 3 g^(1/3) epsilon^(1/3) l^(4/3), divided by 4 for two-dimensional or 6 for three-dimensional diffusion gives about 0.6 or 0.4). [Own calculation]
- Matching: at l = eta the viscous form gives K = eta^2/tau_eta = nu and the inertial form gives epsilon^(1/3) eta^(4/3) = nu, so K(l) is continuous and equals nu at eta. [Own calculation] This is the key point for the model: K(eta) equals the molecular-scale diffusivity, so for Sc near 1 turbulence cannot add much diffusion below eta.
- Competition with molecular diffusion in the viscous range: a blob stretched at rate lambda is thinned to a width where stretching compression balances molecular diffusion, sqrt(D/lambda) = eta/sqrt(Sc) (the Batchelor scale). With Sc = 1, a blob below eta is not further compressed below eta; strain elongates it (increases its length) but its cross-section is set by D. The result: the quantity D_total = D + K(l) with K(l) of order (epsilon/nu) l^2 is smaller than D whenever l is less than eta, so adding it is harmless but also unimportant. [Own reasoning from the cited Batchelor scale]
- Sharp-switch question: the data above say K(l) is a continuous power law, not a switch. A step from 0 to K_eddy at an inner scale is the wrong shape: at 30 m it would jump from approximately D (0.27 m^2/s) to K_eddy (2 to 10 m^2/s), while inertial-range theory gives K(30 m) of 9 to 43 m^2/s for epsilon 1 to 100 mW/kg at unit coefficient (section 2), and K(5 m) below 1 m^2/s. The theory has no "inner scale gate" on K; the inner scale (about 7 to 13 eta, the spectral bend, see section 3) is where the velocity spectrum departs from -5/3, which is inside, not at the edge of, the region where K grows as l^2.

### Gaps
- Could not open Sawford 2001, Salazar and Collins 2009, Batchelor 1952/1959 or Kraichnan 1974 directly; classical equations above are stated from memory and need page-level citation.
- No source found giving a measured numerical Lyapunov-exponent constant (the "order one" prefactor) for lambda tau_eta; DNS literature exists but I did not retrieve it.
- The cited 'g' and the intermittency corrections come from high-Reynolds-number DNS and lab flows, not from stratified mesosphere conditions; their applicability is assumed.

## 2. Effective relative diffusivity versus plume size, and where it reaches K_eddy (outer scale, Ozmidov)

### Takeaway
A smooth replacement for S(r) is K(l) = K_eddy * min[1, (l/l_c)^(4/3)] (inertial law capped at K_eddy), with l_c = (K_eddy^(3/4)) epsilon^(-1/4) / beta^(3/4), and with a viscous-range piece l^2/tau_eta below eta. Because K_eddy is itself of order epsilon/N^2 and K(L_O) = epsilon/N^2 exactly from the inertial law, the cap is reached near the Ozmidov scale, so the plume reaches the full K_eddy only when it is as wide as the largest unstratified eddies (about 10 to 100 m), not at 21 to 42 m by design. For the stated K_eddy and epsilon, though, the cap is reached at about 5 to 30 m at unit coefficient, i.e. the 4/3 law gives MORE diffusion than the current switch at 21 to 42 m for most epsilon, not less.

### Cited Findings
- Eddy diffusion from dissipation: K_ed = b epsilon / N^2 with b = 0.8 is commonly used (Weinstock form); Vlasov and Kelley (2015) show the Weinstock formula follows from similarity theory and implies an outer-scale restriction, with K_ed not exceeding about 3.2e6 cm^2/s = 320 m^2/s for epsilon = 2 W/kg, and eddy scales much less than the scale height; some published K_ed values violate this. Sources disagree on whether the b = 0.8 form is the 1978 or 1981 Weinstock result — [Vlasov and Kelley, Ann. Geophys. 33, 857, 2015](https://angeo.copernicus.org/articles/33/857/2015/); [EarthRef record of Weinstock 1981](https://earthref.org/ERR/45540)
- The 4/3 diffusivity is valid only for scales much larger than the dissipation scale — [arXiv nlin/0107061](https://arxiv.org/pdf/nlin/0107061)
- Mesospheric rocket dissipation rates: roughly 1 to 2 mW/kg below 75 km and 10 to 20 mW/kg in the upper mesosphere and lower thermosphere (Lubken et al. 1993); turbulence covers 60 to 100 km in winter but only 78 to 97 km in summer at high latitude (Lubken 1997); 80 to 90 km rocket/radar values 0.05 to 0.15 W/kg (50 to 150 mW/kg) with occasional peaks near 1 W/kg lasting minutes (Watkins et al. 1988) — [Lubken search results, ACP 15, 2159](https://acp.copernicus.org/articles/15/2159/2015/acp-15-2159-2015.pdf); [Watkins et al. record](https://scilit.com/publications/635bc73d8ae6ebf3fcd70fd62104c305)

### Inferences
[Own calculations] with unit coefficient (beta = 1) in K = epsilon^(1/3) l^(4/3), N = 0.02 rad/s assumed (period about 5 min; mesospheric N is typically 0.02 to 0.03, not sourced here). Scales with beta of 0.4 to 0.6 would make the lengths below larger by a factor 1/beta^(3/4), i.e. about 1.4 to 1.9.

Inertial K(l) at the switch's own scale (m^2/s):

| epsilon (mW/kg) | K at 21 m | K at 30 m | K at 42 m |
|---|---|---|---|
| 1 | 5.8 | 9.3 | 14.6 |
| 10 | 12.5 | 20.1 | 31.5 |
| 100 | 26.9 | 43.3 | 67.8 |

Reading: at the switch scales the 4/3 law gives numbers equal to or larger than the model's K_eddy of 2 to 10 m^2/s, so the stated K_eddy range corresponds to epsilon of at most a few mW/kg if it is interpreted as a plateau at L_O. This is the central consistency check: epsilon/N^2 = 2.5, 25, 250 m^2/s for epsilon = 1, 10, 100 mW/kg. A K_eddy of 2 to 10 m^2/s therefore means epsilon of roughly 1 to 4 mW/kg (or a mixing coefficient b of 0.2 rather than 0.8).

Length at which the inertial law reaches K_eddy, l_c = K_eddy^(3/4) epsilon^(-1/4) (m), versus the Ozmidov scale L_O = (epsilon/N^3)^(1/2):

| epsilon (mW/kg) | l_c for K_eddy = 2, 5, 10 | L_O |
|---|---|---|
| 1 | 9.5, 18.8, 31.6 m | 11 m |
| 10 | 5.3, 10.6, 17.8 m | 35 m |
| 100 | 3.0, 5.9, 10.0 m | 112 m |

Reading: for K_eddy tied to epsilon (K_eddy about epsilon/N^2) the cap sits at L_O, which is 11 to 112 m. The model's 21 to 42 m switch is the right order of magnitude for L_O at epsilon of 3 to 15 mW/kg, but it is the OUTER-scale, not an inner-scale, number. Identity: K(L_O) = epsilon^(1/3) L_O^(4/3) = epsilon/N^2. [Own calculation]

Recommended smooth form (own proposal, not from a source): K(l) = K_eddy * (l^(4/3)/(l^(4/3) + L_c^(4/3))) or the sharper min() version, with L_c about L_O or l_c above, plus the viscous branch l^2/tau_eta below eta (which is below D so can be dropped when Sc is 1). For reference, in this form the effect of the "switch" is nearly identical to using K(l) = beta epsilon^(1/3) l^(4/3) capped at K_eddy; either needs epsilon, which the model must take as an input tied to K_eddy.

Conversion caution: K(l) as written is for pair separation l. The model's r is a plume radius about its centroid. For a Gaussian blob the mean-square pair separation is twice the per-particle variance, so use l about sqrt(2) r (or double the diffusion relation) when mapping; and the model's factor of 4 (r^2 = r0^2 + 4 integral D dt) is the two-dimensional Fickian form. [Own reasoning, [Textbook knowledge]]

### Gaps
- Which coefficient beta to use is unsettled: DNS g = 0.52 versus lab 0.10 (cited), plus the mapping between pair-separation and plume-radius definitions.
- No measured K(l) curve for the mesosphere itself found; the Okubo-style data are oceanic and tropospheric (section 3). The stratified-turbulence effect on K at scales between L_O and the buoyancy scale is not covered here.

## 3. Classical atmospheric and oceanic usage: Gifford, Csanady, Pasquill, Monin and Yaglom, Smith and Hay, Okubo, Thomson; contrails and stratified plumes

### Takeaway
Scale-dependent apparent diffusivity is the standard tool in atmospheric and oceanic dispersion (Richardson 4/3 law, Okubo diagrams, Gifford's Lagrangian theory), but it was built for scales of 100 m to 100 km. I found no source applying it to mesospheric trails or to sub-Kolmogorov plumes; the contrail literature (Duerbeck and Gerz) instead gives effective diffusivities for stratified wake plumes from large-eddy simulation.

### Cited Findings
- Okubo (1971), "Oceanic diffusion diagrams," Deep-Sea Res. 18(8), 789-802: diagrams give the horizontal spread of an instantaneous point release and the apparent diffusivity as a function of patch size, from dye data with time scales of 1 hour to 1 month and length scales of 100 m to 100 km; the earlier 1968 report fits the 4/3 law locally. I could not confirm a commonly quoted l^1.15 exponent from these results (flagged by the search as background knowledge only) — [DTIC 1968 report](https://apps.dtic.mil/sti/html/tr/AD0675269/index.html); [Okubo 1971 record](https://scilit.com/publications/e2d3dee5874a1dd9a8edb6d1913af0e9)
- Gifford (1982), "Horizontal diffusion in the atmosphere: a Lagrangian-dynamical theory," Atmos. Environ. 16, 505-512: a Langevin (random force) theory covering a wide range of diffusion times; reproduces the constant-K Fickian solution plus asymptotic results of both Taylor time-averaged and relative (instantaneous spreading) forms; reproduces accelerating horizontal diffusion; relative diffusion treated with conditioned initial velocities — [Bohrium abstract](https://www.bohrium.com/paper-details/horizontal-diffusion-in-the-atmosphere-a-lagrangian-dynamical-theory/812394205853777921-14619); [LANL report LA-8667-MS](https://sgp.fas.org/othergov/doe/lanl/lib-www/la-pubs/00313710.pdf)
- Gifford's 1959 meandering-plume approach splits total spread into centroid meander plus relative diffusion about the centre of mass; the meandering model fails far downwind where internal (relative) dispersion dominates — [Yee et al. related literature, Springer](https://link.springer.com/doi/10.1007/BF00705506); [AMS conference paper](https://ams.confex.com/ams/pdfpapers/139633.pdf)
- Dürbeck and Gerz (1995 GRL; 1996 JGR) large-eddy simulations of aircraft exhaust plumes in uniform stratification (N 0.006 to 0.03 per second): horizontal effective diffusivity about 11 to 23 m^2/s depending on stratification, vertical about 2.3 m^2/s early falling to 0.15-0.18 m^2/s later and nearly independent of stratification at late times; with shear (Richardson number about 0.1) horizontal rises about 3 times and vertical about 20 times — [Zendy record of 1995 GRL](https://zendy.io/title/10.1029/95gl03021); [DLR elib](https://elib.dlr.de/32304/)
- No result found for Smith and Hay (1961), Pasquill, Csanady (1973) beyond a citation in meandering-plume literature, or Thomson, specifically — see Gaps.

### Inferences
- The accurate historical statement is that Richardson/Okubo/Gifford diffusivities are scale-dependent apparent diffusivities for ensemble clouds: K grows with size because a larger cloud samples larger eddies. This is precisely the S(r) concept, but as a power law not a gate. [Inference from cited findings]
- Aircraft exhaust results show effective horizontal and vertical diffusivities that are very different (factor of order 100 late, 0.15 versus 15 m^2/s), so the model's isotropic K is not defensible in a stratified atmosphere once the plume exceeds the buoyancy scale; see section 5. [Inference from Dürbeck and Gerz]
- Dürbeck and Gerz values (11 to 23 m^2/s horizontal) are for km-scale domains in the stratosphere; they are not a direct anchor for 10 to 100 m scales in the mesosphere.

### Gaps
- Primary text of Okubo 1971, Gifford 1982, Monin and Yaglom, Csanady 1973, Pasquill, Smith and Hay 1961 (Quart. J. R. Met. Soc.), Thomson not retrieved; the viscous-range K proportional to l^2 form is not found stated in any retrieved source (it follows from Batchelor's linear velocity increment).
- No source found applying scale-dependent relative diffusivity to rocket trails or reentry plumes at 60 to 90 km. Rocket chemical-tracer trails (Blamont 1963; Golomb and MacLeod 1966) are reported as the classical way to measure diffusion from trail expansion ([search summary](https://core.ac.uk/works/24876831)) — a worthwhile direct data source for K versus size, not pursued.
- Schumann's contrail-dispersion papers not retrieved (a search hit, [ACP 14, 2713, 2014](https://acp.copernicus.org/articles/14/2713/2014/acp-14-2713-2014.pdf), was not opened).

## 4. Time scale: when does turbulent relative dispersion overtake molecular diffusion? (60, 70, 90 km; epsilon 1, 10, 100 mW/kg)

### Takeaway
All figures here are my own calculations. For these conditions, molecular diffusion and strain-driven spreading are comparable at the Kolmogorov scale by construction; turbulence takes over for plumes larger than about 0.5 to 2 m at 60 km, 1 to 8 m at 70 km and 17 to 53 m at 90 km, reached after about 2 to 20 s, 5 to 48 s and 18 to 175 s respectively. At 90 km turbulence is a weak correction until the plume is tens of metres wide.

### Cited Findings
- Kinematic viscosity inputs (0.05, 0.3 and 4 m^2/s at 60, 70, 90 km) and Sc near 1 (so D about nu) are the assignment's stated inputs; the 70 km molecular diffusivity of 0.27 m^2/s is consistent with them. Dissipation rates 1 to 100 mW/kg bracket the rocket values — [Lubken results above](https://acp.copernicus.org/articles/15/2159/2015/acp-15-2159-2015.pdf); [Watkins et al.](https://scilit.com/publications/635bc73d8ae6ebf3fcd70fd62104c305)
- Measured inner scale at 80 to 90 km was about 20 to 26 m, about 13 times the Kolmogorov microscale, with dissipation 0.05 to 0.15 W/kg — [Watkins et al. 1988](https://scilit.com/publications/635bc73d8ae6ebf3fcd70fd62104c305). Note that 13 eta at 0.1 W/kg and nu about 2 to 4 m^2/s would be greater than 26 m with my formula (about 60 m); the retrieved abstract does not give the viscosity used, so the two are not reconciled here.

### Inferences
[Own calculations]. Definitions: eta = (nu^3/epsilon)^(1/4); tau_eta = (nu/epsilon)^(1/2); L_O = (epsilon/N^3)^(1/2) with N = 0.02 rad/s; inner scale taken as 9.9 eta (the Hill-Clifford/Lubken-type factor, [Textbook knowledge, not re-verified]; the Watkins factor 13 is cited above); crossover time t* from setting the Richardson growth g epsilon t^3 equal to molecular growth 4 D t (two-dimensional, D = nu), giving t* = (4 D/(g epsilon))^(1/2) with g = 0.52 and l* = (4 D t*)^(1/2).

| Altitude | nu (m^2/s) | epsilon (mW/kg) | eta (m) | tau_eta (s) | inner scale 9.9 eta (m) | L_O (m) | t* (s) | l* (m) |
|---|---|---|---|---|---|---|---|---|
| 60 km | 0.05 | 1 | 0.59 | 7.1 | 5.9 | 11 | 19.6 | 2.0 |
| 60 km | 0.05 | 10 | 0.33 | 2.2 | 3.3 | 35 | 6.2 | 1.1 |
| 60 km | 0.05 | 100 | 0.19 | 0.7 | 1.9 | 112 | 2.0 | 0.6 |
| 70 km | 0.3 | 1 | 2.28 | 17.3 | 22.6 | 11 | 48.0 | 7.6 |
| 70 km | 0.3 | 10 | 1.28 | 5.5 | 12.7 | 35 | 15.2 | 4.3 |
| 70 km | 0.3 | 100 | 0.72 | 1.7 | 7.1 | 112 | 4.8 | 2.4 |
| 90 km | 4 | 1 | 15.9 | 63 | 158 | 11 | 175 | 53 |
| 90 km | 4 | 10 | 8.9 | 20 | 89 | 35 | 55.5 | 30 |
| 90 km | 4 | 100 | 5.0 | 6.3 | 50 | 112 | 17.5 | 17 |

Conclusions to draw: (a) at 60 km and 70 km the Kolmogorov scale is below the plume sizes of interest (0.2 to 2.3 m) so a plume of tens of metres is well inside the inertial range and the 4/3 law applies to it; (b) at 90 km eta is 5 to 16 m and the 30 m switch sits inside the viscous-to-inertial transition, so a sharp 30 m inner scale is not wrong in magnitude there, but is wrong in kind; (c) at 90 km and epsilon of 1 mW/kg the Ozmidov scale (11 m) is smaller than eta (16 m), which means no inertial range exists at all (buoyancy-controlled, viscous turbulence); the K_eddy of the model then has no inertial-range meaning there; (d) a 30 m inner scale is a plausible value only at 70 to 90 km with epsilon of order 1 to 10 mW/kg (inner scale 13 to 23 m at 70 km, 89 to 158 m at 90 km) and is far too large for 60 km (3 to 6 m). Altitude-dependent inner scale is therefore needed: one number at 30 m for all of 60 to 90 km is not supported.

The time to overtake molecular spreading is short compared with plume ages of minutes to hours, so unless the plume is observed in its first seconds, turbulent relative dispersion matters as soon as K(l) for the current plume size exceeds D, which per section 1 needs l larger than about eta. Assumption: unit-coefficient laws, N = 0.02 rad/s, D = nu, a stationary, homogeneous epsilon; none of these were measured here.

### Gaps
- The Watkins "13 times" inner scale versus my 9.9 eta values is not reconciled; inner scale definitions (Hill versus Tatarskii) differ by tens of percent.
- N(z) and epsilon(z) for the actual reentry profile are unknown to me; N = 0.02 rad/s is an assumption.
- Temporal intermittency: Watkins reports peaks near 1 W/kg lasting only minutes, so a stationary epsilon over the plume age (minutes to hours) is questionable.

## 5. Line versus point source; stratification and anisotropy; meandering versus true mixing; chemistry rates and segregation

### Takeaway
A line (cylinder, 2D) source differs from a point (3D) source only through the number of dimensions entering the relative-dispersion law (factor 4 versus 6 in the Fickian conversion, and a line is affected only by the cross-section strain), not by a different exponent; stratification makes the vertical K smaller than the horizontal K once plume scales exceed the buoyancy scale, so one isotropic K is not defensible; and a mean-square radius rising with time tells nothing about how fast peak concentration falls, because meandering of the plume centroid increases r^2 (when measured against a fixed frame) without mixing, so mean-field second-order chemistry can be off by a large factor.

### Cited Findings
- Meandering versus relative diffusion: total spread = centroid meander plus relative dispersion about the centroid (Gifford 1959); meandering models ignore the internal vortices and fluctuations and fail far downwind, and in the far field relative dispersion overcomes meander — [AMS conference paper](https://ams.confex.com/ams/pdfpapers/139633.pdf); [Springer abstract](https://link.springer.com/doi/10.1007/BF00705506)
- Concentration fluctuations inside the plume: Yee et al. (1994) added relative in-plume fluctuation intensity; fluctuating-plume models were extended to reacting species, and segregation cannot be neglected for short-term concentrations when reactions occur before turbulence has mixed the pollutants; bivariate gamma distributions are used to model internal mixing — [Yee and related, springer](https://link.springer.com/doi/10.1007/BF00705506); [Bivariate gamma for chemically reacting plume](https://research.uniupo.it/en/publications/application-of-a-bivariate-gamma-distribution-for-a-chemically-re/); [HARMO paper](https://www.harmo.org/Conferences/Proceedings/_Madrid/publishedSections/PPT/H15-159-O.pdf)
- Stratified plume LES (Dürbeck and Gerz): late-time vertical effective diffusivity (0.15 to 0.18 m^2/s) is two orders of magnitude below horizontal (15 to 23 m^2/s), and vertical diffusivity and shear control late plume size and entrainment — [Zendy 1995](https://zendy.io/title/10.1029/95gl03021); [DLR elib](https://elib.dlr.de/32304/)
- Turbulence below the Ozmidov scale in the mesosphere is lower-limited by L_O and small-scale structures there last about a minute versus up to hours for large scales above 75 km (lecture slides, not peer-reviewed) — [Gerding and Avsarkisov slides](https://www.trr-energytransfers.de/fileadmin/user_upload/documents/events/WSP-April-2018/WSPApril2018_Avsarkisov.pdf)

### Inferences
- Line versus point: a line source (a long reentry trail) is an initially one-dimensional object, but once its cross-section is smaller than the along-trail variation scale the relative dispersion in the cross-plane obeys the same K(l) in two dimensions. Cross-section growth: <l^2> grows as the same t^3 law; the Fickian conversion uses 4K rather than 6K (the model's "4" is the 2D one). Along-trail the trail is long compared with eddies so contortions and folds increase the surface area without changing the cross-section law. Not specifically sourced. [Inference; textbook level]
- Stratification: below L_O turbulence is near isotropic; between L_O and the outer scale it becomes layered with vertical K about N^-2 epsilon-limited and horizontal K larger. For plumes of 20 to 100 m at 60 to 90 km, L_O of 11 to 112 m (own calculation, section 4) means a plume reaches the anisotropic regime at about the same time it reaches K_eddy. A defensible representation is K_h(l) with the capped 4/3 form and K_v(l) = min[K(l), epsilon/N^2-based cap or smaller] — e.g. cap K_v at about 0.2 epsilon/N^2 and K_h at several times that. The Dürbeck and Gerz ratio of about 100 is for km-scale plumes and an exhaust wake, and should not be transplanted without caveat. [Inference]
- Mean-square radius and peak concentration: r^2 about the instantaneous centroid, which is what relative dispersion measures, grows as g epsilon t^3. A fixed-frame ensemble r^2 includes meander and is larger. Peak concentration falls with the instantaneous plume cross-section, so the ensemble-mean concentration (which Gaussian models and the model's r^2 produce) is lower in peak and broader than instantaneous values. A model fitted to ensemble-mean r^2 therefore overstates dilution of an individual plume. [Inference from Gifford meander split]
- Chemistry: for an A + B reaction the rate is k times the mean of C_A C_B; with turbulent fluctuations this is k (mean C_A mean C_B + mean c_A' c_B'), and the covariance term is negative for segregated reactants (a vapour plume and ambient oxidant that have not yet interpenetrated). Using mean concentration in a plume whose internal structure is filamentary overestimates reaction rate when scales are above the Batchelor scale, because the stirred-but-unmixed striations are thin sheets whose interiors are shielded. Direction of the error, not the magnitude, is what the cited literature supports. [Inference; cited findings say only that segregation matters]

### Gaps
- No source quantifying the segregation correction (Damkoehler number or reaction-rate reduction factor) for a mesospheric vapour plume; need Damkoehler-number estimate using local second-order rate and the mixing time (Kolmogorov time of 0.7 to 63 s from the section 4 table).
- No source for isotropic versus stratified K in the mesosphere at scales of 10 to 100 m.
- Could not retrieve Pope's Turbulent Flows or Monin and Yaglom for formal definitions of line-source versus point-source relative dispersion, or for the mixing-time-for-a-blob result (Batchelor scale logarithmic mixing time of order tau_eta ln(l0/eta_B)) which should be checked: it implies a plume of size l0 above eta needs only a few tau_eta (a few to tens of seconds) to be thinned to the diffusive scale. [Textbook knowledge, not re-verified]

---
Summary of recommendation for the model (from the above; the decision belongs to the researcher): (1) replace the 21-to-42 m switch S(r) with the continuous form D_total = D + K(l) where K(l) = min[K_eddy, beta epsilon^(1/3) l^(4/3)] using l about sqrt(2) r, (2) tie K_eddy and epsilon together through K_eddy about b epsilon/N^2 with b of 0.2 to 0.8 and take L_O = (epsilon/N^3)^(1/2) as the cap scale, (3) make epsilon altitude dependent (1 to 2 mW/kg below 75 km, 10 to 20 above, per Lubken 1993), (4) note that the suspicion that the switch is too large is not clearly supported: at the same epsilon and unit coefficient the 4/3 law gives more diffusion at 21 to 42 m (6 to 68 m^2/s) than the model's K_eddy of 2 to 10, the real problem is that K_eddy and the 30 m scale are inconsistent with each other and with epsilon; the likely correct fix is to choose epsilon first.
