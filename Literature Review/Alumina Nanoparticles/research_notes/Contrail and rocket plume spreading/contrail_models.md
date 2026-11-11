# Contrail and aircraft plume models: spreading of a thin line source in stratified, sheared air

Scope and provenance. The only source read in full text was Schumann (2012, CoCiP; GMD 5, 543), plus Naiman et al. (2010, ACP 10, 2551) and Unterstrasser and Gierens (2010, ACP 10, 2017). For Schumann et al. (1995), Dürbeck and Gerz (1996), Gerz et al. (1998), Freudenthaler et al. (1995), Pisso et al. (2009), Konopka (1995), Kraabøl et al. (2000) and Kärcher (2018) only abstracts or bibliographic records were retrieved; numbers from them are abstract-level. Schumann (1996) and Schumann et al. (1998, the dilution law) were not retrievable; the 1998 law is quoted through CoCiP and Naiman. Conditions throughout: cruise altitude about 9.4-11.3 km, Brunt-Vaisala frequency N_BV (buoyancy frequency, 1/s) about 0.01-0.02 per second, shear S about 0.001-0.01 per second. Equations are quoted from CoCiP and the text extraction of the PDF is partly garbled; one equation is flagged below.

## 1. Gaussian plume model with horizontal, vertical and shear diffusivities (exact equations)

### Takeaway
The standard model (Konopka 1995, used by Schumann et al. 1995, Dürbeck and Gerz 1996, CoCiP) is an exact analytic solution for a Gaussian line-source cross-section in a constant vertical shear with constant diffusivities. It tracks a 2x2 covariance matrix, not one radius, and the cross-section area grows only through the combination of shear and vertical diffusion. The model is slender (cross-section plane only; along-axis diffusion neglected).

### Cited Findings
- Concentration in the plane normal to the trail: c(x) = (C0/A) exp[-(1/2) x^T sigma^-1 x], with x = (y, z) (y horizontal cross-trail, z vertical, in m), C0 the mass of species per unit trail length (kg/m, conserved), sigma the 2x2 covariance matrix (sigma_yy, sigma_yz, sigma_zz, all in m^2), and effective cross-section A = 2 pi sqrt(det sigma) with det sigma = sigma_yy sigma_zz - sigma_yz^2 (Eqs. 1-7). Width and depth are defined as B^2 = 8 sigma_yy and D^2 = 8 sigma_zz, so B = 2 sqrt(2) sigma_y; effective depth D_eff = A/B. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Constant S (1/s), D_H, D_V, D_S (m^2/s) give, per time step dt (Konopka 1995 solution, quoted as Eqs. 29-31 with segment-stretching factor L(t)/L(t+dt) for trail length L): sigma_yy(t+dt) = [ (2/3) S^2 D_V dt^3 + (S^2 sigma_zz + 2 D_S S) dt^2 + 2 (D_H + S sigma_yz) dt + sigma_yy ] (L(t)/L(t+dt))^2; sigma_zz(t+dt) = 2 D_V dt + sigma_zz; sigma_yz(t+dt) = [ S D_V dt^2 + (2 D_S + S sigma_zz) dt + sigma_yz ] L(t)/L(t+dt). Positive-definiteness requires D_S^2 <= D_V D_H (Eq. 28). Shear is S = d(V_n)/dz, the vertical gradient of the plume-normal horizontal wind. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Analytic cross-section area (Eq. 32): A(t+dt) = 2 pi { (1/3) S^2 D_V^2 dt^4 + (2/3) S^2 D_V sigma_zz dt^3 + (2 S D_V sigma_zz - 2 S D_S sigma_zz + 4 D_H D_V - 4 D_S^2) dt^2 + (2 D_V sigma_zz + 2 D_H sigma_zz - 4 D_S sigma_yz) dt + sigma_yy sigma_zz - sigma_yz^2 }^(1/2). Same as Dürbeck and Gerz (1996) Eq. 9 up to a factor of 2. Shear alone does not increase A (the ellipse gets wider and thinner); with D_V = 0, A = 2 pi [2 D_H sigma_zz dt + sigma_yy sigma_zz - sigma_yz^2]^(1/2) regardless of shear. Shear enhances mixing "considerably" only together with vertical diffusion; D_S matters less as shear increases. Shear diffusivity enlarges A only if S D_S < 0. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Late-time limit read from Eq. 29 (an inference from the equation, see Inferences): with large t the cubic term dominates, sigma_yy ~ (2/3) S^2 D_V t^3.
- Dissipation of Gaussian to elliptical models: Naiman et al.'s subgrid plume model advances an ellipse with radii a, b and angle theta: d theta/dt = (du_s/dz) cos^2 theta; da/dt = a (du_s/dz) sin theta cos theta; db/dt = -b (du_s/dz) sin theta cos theta (area conserved by shear); diffusion: da/dt = D_a/a, db/dt = D_b/b, with D_a = D_v cos theta + D_h sin theta, D_b = D_v sin theta + D_h cos theta (Eqs. 5-9). Mapping to variances: sigma_v^2 = a^2/4 cos^2 theta + b^2/4 sin^2 theta; sigma_h^2 = a^2/4 sin^2 theta + b^2/4 cos^2 theta; sigma_s^2 = (a^2/4 - b^2/4) cos theta sin theta. Compared against Konopka, plume area after 10 h at case-1 shear was 28 times initial (ellipse model) versus 24 times (Konopka). — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html)
- Konopka (1995, J. Non-Equilib. Thermodyn. 20, 78-91) derives Gaussian solutions for linear shear flow with time-dependent diffusion tensor, sedimentation velocity and linear sources and sinks. — [Scilit record](https://scilit.com/publications/0005cf16c1ff4c0c16b8735e61c4936c)
- CoCiP states that any diffusivity that depends on plume scale makes these equations non-linear so the solution is no longer exact and becomes time-step dependent, but stable and positive-definite; time steps of 60 s versus 3600 s differed by about 10 % in dilution. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)

### Inferences
- For the aluminium trail (round, isotropic) the Konopka form reduces to two scalars if S = 0: sigma_yy = sigma_zz = sigma, d sigma/dt = 2D. The user's r^2 = r0^2 + 4 integral D dt corresponds to r^2 = 2 sigma (a 1/e radius of exp(-r^2/(4Dt))), whereas CoCiP's B/2 = sqrt(2 sigma_yy), so the same factor-of-2 convention applies: B/2 is the 1/e half-width. Check this convention before comparing numbers.
- With shear, the area-growth law A ~ S D_V t^2 (long time) implies the mixed-in volume is controlled by the product S sqrt(D_V) times t^(3/2) in width, not by an isotropic D. A single isotropic D cannot reproduce this; two diffusivities plus S are the minimum.

### Gaps
- Schumann (1996) was not located; CoCiP and Naiman credit Konopka (1995) and Schumann et al. (1995) for this model, so the 1996 paper could not be used directly.
- Konopka (1995) full text not read; equations above come through CoCiP.

## 2. Values of D_h, D_v, D_s and S, and their relation to epsilon and N

### Takeaway
At cruise level the measured and LES-derived vertical diffusivity is tiny (0-1 m^2/s, 0-0.6 at plume scales) and the horizontal diffusivity is 5-20 m^2/s for young plumes, growing to about 10^4 m^2/s for plumes tens of km wide. D_S is essentially unconstrained and set to zero. CoCiP's closure ties D_V to N_BV and a fixed turbulent velocity scale, and D_H to depth squared times shear.

### Cited Findings
- Schumann et al. (1995, JGR 100, 14147): ten plumes sampled 5-100 min after emission at 9.4-11.3 km, bulk Richardson number > 10, NO half-widths 500-2000 m; mean vertical plume width about 140 m; vertical diffusivity 0-0.6 m^2/s; horizontal diffusivity 5-20 m^2/s from growth of horizontal scales; turbulence data had practically zero dissipation and weak vertical velocity; with constant diffusivities shear dominates lateral dispersion after about 1 h, even for weak shear 0.002 per second. — [DLR eLib record](https://elib.dlr.de/32042/)
- Dürbeck and Gerz (1996) baseline LES case (case 2): initial a0 = 184 m, b0 = 260 m, D_h = 20.0 m^2/s, D_s = 0.75 m^2/s, D_v = 0.158 m^2/s, shear du/dz = 0.003 per second; seventy-minute runs; their shear cases ran 0.001-0.007 per second. At Richardson number near 0.1, D_h increases about 3 times and D_v about 20 times relative to no shear; shear plus diffusion enlarges cross-section and entrainment; agreement with North Atlantic in situ data. — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html); [Scilit abstract](https://scilit.com/publications/37805b3906274d4af9e7a8e694d9103e)
- CoCiP closure (Eqs. 35-37, 38): D_V = c_V (w'_N^2 / N_BV) + f_T V_T D_eff with c_V = 0.2, f_T = 0.1 (f_T V_T term is ice sedimentation, V_T fall speed); D_H = c_H D^2 S_T with c_H = 0.1 and S_T the total shear sqrt((dU/dz)^2 + (dV/dz)^2); D_S = 0 (or D_S = -C_DS sqrt(D_V D_H) S/|S|); w'_N = 0.1 m/s is the prescribed rms vertical turbulent velocity; N_BV floored at 0.001 per second. The text quotes epsilon = A_S w'_N^2 S (A_S = 0.5; the extracted PDF text of this equation is garbled, so the exact form was not verified) from Hunt et al. (1988) and Schumann and Gerz (1995). — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Arithmetic check with those formulas (my calculation): w'_N = 0.1, N_BV = 0.01-0.03 gives c_V w'^2/N = 0.07-0.2 m^2/s, matching CoCiP's quoted D_V of about 0.07-0.2; D = 100 m and S = 0.001-0.01 gives D_H = 1-10 m^2/s, matching quoted values.
- CoCiP quoted ranges: S = 1e-3 to 1e-2 per second, N_BV = 0.01-0.03 per second, D = 100 m give epsilon about 5e-7 to 0.05 m^2/s^3, buoyancy scale L_B about 3-10 m, Ozmidov scale L_Oz = sqrt(epsilon / N^3) about 0.07-0.7 m, Kolmogorov scale 0.3-6 mm. Dissipation may span 1e-8 to 1e-2 m^2/s^3. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Free-atmosphere D_V typically 0-1 m^2/s (Pavelin et al. 2002, via CoCiP); young-plume D_H 5-20 m^2/s growing to order 10^4 m^2/s for wider, thicker plumes; D_H is about 10-100 times D_V; early-period D_V 0.01-0.5 m^2/s. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Pisso et al. (2009, JGR 114, D19301), free troposphere from ensembles of diffusive backward trajectories for a long-range pollution plume: D_v about 0.3-1 m^2/s and D_h about 1e4 m^2/s. — [IP Paris research portal record](https://researchportal.ip-paris.fr/en/publications/estimation-of-mixing-in-the-troposphere-from-lagrangian-trace-gas/)
- Unterstrasser and Gierens (2010) ambient turbulence in their LES: eddy dissipation rate 3.5e-5 m^2/s^3, rms velocity fluctuation 0.25 m/s, N_BV = 0.01 per second; Dürbeck-Gerz test with initial sigma_hh,0 = 1.4e4 m^2, sigma_vv,0 = 7e3 m^2, N_BV = 0.019 per second, S = 3e-3 per second: sigma_hh grew by 2-3 orders of magnitude while sigma_vv grew only 10-30 %. — [Unterstrasser and Gierens 2010](https://acp.copernicus.org/articles/10/2017/2010/acp-10-2017-2010.html)
- Shear is taken from NWP at coarse resolution and boosted: f_S = (1 + (dz_eff/D)^n)/2, dz_eff = 2000 m (ECMWF effective shear resolution), n = 1/2 (n = 2/3 for inertial-range turbulence, n = 0 for steep stratified spectra). — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)

### Inferences
- The D_V closure is not the classic D_V = epsilon/N^2 form directly, but it is equivalent through the quoted epsilon relation: with epsilon about w'^2 S, D_V = c_V w'^2/N corresponds to a flux-Richardson-like efficiency times epsilon/N^2 up to the S/N factor. A strict Osborn-type D = Gamma epsilon / N^2 (Gamma about 0.2) is the usual textbook equivalent, but this was not found stated in the retrieved sources.
- D_H = c_H D^2 S is a shear-dispersion-style scale dependence (D_H grows as depth squared), different from the user's Richardson law (D_R ~ r^(4/3)). Both increase with plume size; the Schumann values at D = 10 km would be of order 10^3-10^4 m^2/s for S = 0.001-0.01, consistent with the quoted 10^4.

### Gaps
- No source retrieved states D_S measured; only the 0.75 m^2/s LES-fit and the constraint D_S^2 <= D_V D_H.
- Dürbeck and Gerz (1996) full table of D_h, D_v versus Richardson number not retrieved.

## 3. Measured widths and depths against age; which term dominates; shear importance; dilution laws

### Takeaway
Observed contrail width grows roughly linearly in time at tens of m/min in the first hour, reaching 1-10 km after about an hour and up to about 20-22 km after 3 h; depth is limited by the humid layer. Shear dominates lateral growth after about an hour even if weak. The empirical bulk dilution law is N(t) = 7000 (t/1 s)^0.8 with factor-of-3 scatter.

### Cited Findings
- Freudenthaler et al. (1995, GRL 22, 3501): ground lidar, 81 contrails 1-60 min old; horizontal growth 18-140 m/min; vertical growth up to 18 m/min (often capped by the humid layer); cross-section growth 3500 and 25000 m^2/min. Unterstrasser and Gierens note they fitted linear-in-time growth, so the lidar comparison disagrees most at 60 min against nonlinear model growth. — [GRL abstract record](https://earthref.org/ERR/15164/) (abstract numbers via search summary); [Unterstrasser and Gierens 2010](https://acp.copernicus.org/articles/10/2017/2010/acp-10-2017-2010.html)
- Dilution (Schumann et al. 1998, from more than 70 plume observations): N(t) = 7000 (t/t0)^0.8, t in seconds, t0 = 1 s, valid for ages under about 5000-10000 s, individual observations within a factor of 3. N is the ratio of the plume air mass to the emitted exhaust, so excess concentration = EI/N_dil. — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html); [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- CoCiP passive-tracer dilution depends on shear and stratification: stronger for smaller aircraft, stronger shear, weaker stratification; order of magnitude agrees with the empirical law; at one hour a ten times larger D_V gives roughly 3-4 times more dilution; a half and ten times change in D_V alters late dilution by up to a factor 7 in concentration. Reference test cases S = 0.002 per second and N_BV = 0.01 per second; sensitivity to S = 0 and 0.005 per second. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Dürbeck and Gerz LES: area grew by a factor of eight over seventy minutes at the high-shear case (0.007 per second); the area growth is closer to sqrt(t) at low shear and t^2 at high shear in the Konopka solution, but "more linear" in the LES. — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html)
- Gerz, Dürbeck and Konopka (1998, JGR 103, 25905): LES of a large aircraft wake, vortices break up at 1.5-3 min; plume extent 20 km horizontally and 0.3 km vertically at late times; an extra record states dilution to background in 2-12 h for shear 0.002-0.01 per second. — [DLR eLib record](https://elib.dlr.de/8922/) (abstract-level)
- Unterstrasser and Gierens (2010): after about one hour contrail widths typically 1-10 km (via a 2025 preprint citation); substantial visible spreading only if relative humidity over ice is above about 120 %; shear affects optical properties less than expected. — [ACP abstract](https://acp.copernicus.org/articles/10/2017/2010/)
- Satellite spreading rates: about 2.7 km/h (Great Lakes case study, GOES), 7-8 km/h in a post-2001 shutdown study, width about 22 km 3 h after first detection. — [Duda et al. summary via NASA ESPO](https://espo.nasa.gov/glopac/node/15742)
- Schumann (2012): at the time of shear dominance, "width grows stronger than linear" because rising depth makes shear dispersion more efficient; sedimentation mainly increases depth. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)

### Inferences
- Over 12 h and 1 km target width of the aluminium trail, a growth of 1 km in 12 h is about 0.02 m/s, far slower than contrail widening of 0.3-2 m/s, because contrail growth is driven by shear acting on a layer 100 m deep and by ice-supersaturated dynamics. The comparison suggests shear should be included if the mesospheric/thermospheric shear is of order 10^-2 per second.
- The contrail 1 h width of 1-10 km exceeds the aluminium 1 km 12 h target, so the contrail result bounds shear spreading from above; the comparable aluminium case is likely closer to the weak-shear contrail cases.

### Gaps
- Lidar numbers (18-140 m/min) are abstract-level; no table of width versus age with shear values was retrieved.
- Observed widths for individual cases split by shear value not found.
- Schumann et al. (1998) not retrieved.

## 4. Is width controlled by shear rather than turbulence? Line source aligned or across the shear

### Takeaway
Yes, in the contrail literature width is shear-controlled after roughly an hour; turbulence (D_H) controls early growth. Because the model is a slender cross-section model, only the shear component normal to the trail matters (S = d V_n/dz), so a trail aligned with the wind shear vector would see no shear spreading from that component. No paper retrieved tests the aligned versus cross-shear cases explicitly.

### Cited Findings
- With constant diffusivities shear dominates lateral dispersion after about 1 h even for 0.002 per second. — [Schumann et al. 1995](https://elib.dlr.de/32042/)
- Naiman et al.: the difference in vertical variance development between models has little effect on the extent and horizontal spreading, "dominated by the vertical shear". Higher vertical shear gives faster area growth. — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html)
- Plume distortion follows the vertical shear of the plume-normal velocity S = d V_n/dz, whereas turbulence production is driven by total shear S_T (Eq. 33). The slender plume approximation limits shear to the plane perpendicular to the segment. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf); [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html)
- When D_V = 0, shear does not increase A at all. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)

### Inferences
- A trail aligned with the shear vector (S = 0 in the normal plane) widens only by D_H and D_V, and its length is stretched by the along-trail shear rather than spread laterally. For a trail 40 m to 20 km long and nearly vertical (as a descending meteor trail is), the geometry differs from a horizontal contrail: horizontal winds that vary along the vertical trail shear the trail itself, so the relevant "shear" for a vertical line source is a horizontal wind gradient along the line axis, which tilts and stretches the line instead of fattening the cross-section. This vertical-line case is not treated in any retrieved paper.

### Gaps
- No source retrieved addresses aligned versus cross-shear line sources or vertical line sources.

## 5. Early vortex phase and what does not apply to a passive trail

### Takeaway
The wake vortex, jet and downwash phases (first about 2-3 minutes, with sinking of order 100-300 m) set the initial width and depth and do not apply to a droplet's vapour trail; CoCiP's Gaussian begins after that, at t0 of order 30 s plus break-up.

### Cited Findings
- Vortex scales: b0 = (pi/4) s_a (s_a is wingspan, 79.8 m for an A380), t0 = 2 pi b0^2 / Gamma0, w0 = b0/t0; A380 example b0 = 62.7 m, t0 = 30.3 s, w0 = 2.07 m/s, N* = N_BV t0 = 0.363, maximum sinking dz_w = 290 m versus 270 m observed. Parameterization: if N_BV t0 >= 0.8, dz_w = 1.49 w0/N_BV; else dz_w/b0 = 7.68 (1 - 4.07 eps* + 5.67 eps*^2)(0.79 - N_BV t0) + 1.88. Initial depth D1 = C_D0 dz_w with C_D0 = 0.5; initial sinking C_z1 = 0.25 dz_w; initial width from B1 = N_dil(t0) m_F / [(pi/4) rho D1]. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Most exhaust is trapped in the two vortices (primary wake), entrainment rates of primary and secondary wakes differ by orders of magnitude; vortices break up at 1.5-3 min. — [Gerz et al. 1998 record](https://elib.dlr.de/8922/)
- Naiman et al. start comparison at t = 300 s, when the vortex system "can be assumed to have dissipated"; Unterstrasser and Gierens start dispersion at 120 s with break-up near 135 s. — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html); [Unterstrasser and Gierens 2010](https://acp.copernicus.org/articles/10/2017/2010/acp-10-2017-2010.html)
- Ice-related parts: initial ice mass (water emitted plus entrained humidity), crystal survival f_surv, sedimentation (f_T term in D_V), latent heat, radiative heating; Unterstrasser and Gierens state that for contrail vertical expansion "radiation, release of latent heat and sedimentation" are more important than turbulent diffusion. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf); [Unterstrasser and Gierens 2010](https://acp.copernicus.org/articles/10/2017/2010/acp-10-2017-2010.html)

### Inferences
- Do not carry over: vortex sinking and its scaling with N_BV and epsilon, B1 and D1 formulas, N_dil(t0) matching, ice microphysics, sedimentation term in D_V, buoyancy of hot exhaust (Pisso et al. note buoyant plumes locally raise D_V slightly).
- Carry over: the passive-tracer Gaussian cross-section with S, D_H, D_V, D_S after the plume loses its own dynamics; the shear amplification of D_V; the dilution-versus-age benchmark.
- N_BV depends on altitude: the cruise numbers (0.01-0.02 per second) do not apply directly at 50-90 km, where N_BV is of order 0.02 per second in the stratosphere/mesosphere but turbulence, shear and viscosity differ; the closure form D_V = c_V w'^2/N is stratification dependent and needs mesospheric w' or epsilon (the user's Megner 2-10 m^2/s is much larger than the cruise 0-1 m^2/s because of gravity-wave breaking and the low density).

### Gaps
- No retrieved source describes a passive plume starting without an initial vortex phase.

## 6. LES closures for subgrid scale, and scale dependence of D_H with plume size

### Takeaway
LES studies resolve most mixing explicitly (grid 10 m and smaller) with a TKE subgrid model; Eulerian LES has numerical diffusion and cannot identify scale-dependent Lagrangian mixing, so LES-versus-Gaussian-model differences should not be used to tune diffusivities. Observed D_H grows by about three orders of magnitude from small to large plumes.

### Cited Findings
- Unterstrasser and Gierens use EULAG with MPDATA advection (positive definite, minimal dispersion) and a TKE subgrid model; mesh size of order 10 m is much smaller than contrail dimensions so diffusion is mainly by resolved fluctuations, only to lesser extent subgrid turbulence; a Smagorinsky test showed negligible differences; a 2-D model was validated against the 3-D Dürbeck and Gerz passive test. — [Unterstrasser and Gierens 2010](https://acp.copernicus.org/articles/10/2017/2010/acp-10-2017-2010.html)
- Schumann: LES "are prone to numerical diffusion, in particular in the presence of shear and stratification", compute mixing in an Eulerian frame, and "do not identify the Lagrangian and scale-dependent mixing". — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Scale dependence: D_H 5-20 m^2/s for young plumes (hundreds of metres) up to order 10^4 m^2/s for wide plumes; D_H = 0.1 D^2 S_T makes D_H scale with depth squared; replacing D by the width B would make D_H "unrealistically large". Kinetic energy spectra of the region follow k^-3 at large scales and k^-5/3 at smaller scales; overturning only below the Ozmidov scale (0.07-0.7 m for the cruise numbers). Wave motions advect and distort without mixing. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Naiman et al.: in a global model the transition to cirrus is represented by raising the diffusion coefficients passed to the plume model to the plume-scale levels of Schumann et al. (1995). — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html)

### Inferences
- A Richardson 4/3 law K_R = 0.42 epsilon^(1/3) r^(4/3) is in the same spirit as these scale-dependent forms but it is isotropic; the cruise literature uses strongly anisotropic D_H/D_V = 10-100 and caps vertical mixing at small values; the user's `min(K_eddy, K_R)` treats K_eddy as a ceiling on a single D, which has no analogue here.

### Gaps
- No explicit exponent for D_H versus scale (such as 4/3) was retrieved for aircraft plumes; Pisso et al. full text (HTTP 500) was not read.

## 7. Plume concentration for chemistry: effective cross-section and segregation

### Takeaway
The Gaussian models define an effective area A = 2 pi sqrt(det sigma) so that the plume-mean concentration is C0/A and the dilution ratio N_dil = rho A / (mass per length); plume chemistry has been done with expanding layered or box plume models and the plume-averaged results used to correct grid-box emission.

### Cited Findings
- Peak concentration C0/A and total conserved mass C0 (Eqs. 1-4); excess mass-specific concentration = EI/N_dil, N_dil about 7000 (t/1 s)^0.8. — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)
- Naiman et al. define dilution as plume area normalized by its initial area, and then add plume contents to the grid scale once dilution passes a threshold; their ellipse plume is "quasi-uniform" while real plumes lie between uniform and Gaussian. — [Naiman et al. 2010](https://acp.copernicus.org/articles/10/2551/2010/acp-10-2551-2010.html)
- Kraabøl and Stordal (2000): expanding multi-layer (circular or cylindrical) plume model for NOx chemistry tested against two B747 POLINAT plumes (NOx emission indices 12.5 and 26.1 g/kg); NOx left in the plume after 15 h 58-69 % for the 26.1 g/kg case; plume processing reduced grid-mean aircraft NOx by 40 % or less in a regional model (ozone by up to about 30 %), and 25-35 % in a global follow-up (2002). — [Kraabøl and Stordal records](https://elib.dlr.de/53822/)
- Schumann (2012) notes the plume model has been used for air chemistry, sometimes in regional and global models (Schumann and Konopka 1994; Karol et al. 1997; Kraabøl et al. 2000). — [Schumann 2012](https://gmd.copernicus.org/articles/5/543/2012/gmd-5-543-2012.pdf)

### Inferences
- For a reacting vapour (aluminium species) the effective area A fixes the in-plume concentration, but with a Gaussian profile reaction rates averaged over the profile differ from rates at the mean concentration (segregation); the layered Kraabøl approach resolves this at higher cost. No retrieved source quantifies a segregation factor for self-reacting species.

### Gaps
- No quantitative segregation or intensity-of-segregation values found; the Kraabøl full text was not read.

## Summary for the closure design (inferences, not cited findings)

1. Replace the single isotropic D with Konopka's covariance evolution (Eqs. 29-31) driven by D_H, D_V, S, with D_S = 0; this reproduces the well-supported result that cross-section area grows through S times sqrt(D_V) (sigma_yy ~ (2/3) S^2 D_V t^3).
2. Plausible mapping of the user's K_eddy to D_V (and K_R-type scale law to D_H) needs mesospheric S and N_BV; cruise values (S about 0.001-0.01 per second, N_BV about 0.01-0.02) cannot be carried over unchanged. CoCiP's D_V = c_V w'^2/N_BV and D_H = 0.1 D^2 S_T show the functional dependences to test.
3. Calibration targets: dilution law 7000 (t/1 s)^0.8 within factor 3, lidar width growth 18-140 m/min over the first hour, D_V 0-0.6 m^2/s at plume scale, D_H 5-20 m^2/s young and 1e4 m^2/s wide.
