# Sideways separation of a fragment released beside a hypersonic parent body


> **Update 2, 2026-10-29: Register et al. (2020) read in full** (a clean copy that Asha saved as
> `Literature Review/Wake Growth/1-s2.0-S001910351930452X-main.pdf`; the earlier copies had unembedded fonts).
> This section refines section 0.

## 00. Register et al. (2020), Icarus 337:113468 — the AFFIx model

**What they did.**
- Cart3D inviscid force databases for two spheres, with viscous corrections from US3D for the near wake.
- Diameter ratios r = D2/D1 of 1, 1/2, 2/5, 1/3, 1/4, 1/6, 1/8 and 1/16, covering every initial angle θ from
  0 to 180 degrees. θ is measured between the parent's flight axis and the line joining the two centres;
  θ = 0 is directly behind, θ = 180 directly ahead.
- Validated against free-flight tests of two 5 cm spheres in the DLR H2K hypersonic wind tunnel. Streamwise
  separation matched to a fraction of a millimetre in the near-wake test.

**Stated validity (p. 9).**
- The air is a continuum (no rarefied-gas effects).
- The free stream is much faster than the relative motion of the fragments.
- The Mach number is high (hypersonic similarity).
- Both fragments have the same bulk density.

**Outcome: the final zone z_f.** This is the zone the smaller sphere ends in and never leaves (Eq. 5; Table 3):
- **z_f = 3, near wake, "trapped behind F1".** Applies when θ < θ_inn(r). For r ≤ r* = 0.353,
  θ_inn = 122 r + 47 degrees.
- **z_f = 2, far wake.** The small sphere stays in the wake while its streamwise separation grows. Applies when
  θ_inn ≤ θ < θ_bif - s_w, with θ_bif = 100 exp(-1.618 r) + 40 degrees and a shock-band half-width
  s_w of about 1 to 3 degrees.
- **z_f = 1, carried along the shock.** Applies within s_w of θ_bif.
- **z_f = 0, independent flight.** Applies above θ_bif + s_w, if the flight lasts long enough.

**Interaction time (Eqs. 4, 7-9).**
- Dimensional time: t_fz = D1 (rho/q_inf)^(1/2) t_hat. Here rho is the bulk density of both bodies and q_inf
  the dynamic pressure, so (rho/q_inf)^(1/2) has units of seconds per metre.
- Near wake: t_hat = 0.
- Far wake: t_hat = A0(r) + 0.5 exp[-0.05(θ - θ_inn)], with A0 = -7.2 r^2 + 3 r + 0.1 for r ≤ r*.
- Independent flight: t_hat = 0.88 r + 0.32 + 0.01 exp[0.125(θ - θ_bif - 10)].

**Two findings that bear directly on our model.**
- **Drag is much lower in the near wake.** A small sphere trapped there feels a force coefficient of about 0.12,
  against 0.88 in free flow, about 0.14 times as much (p. 18, Fig. 27). The sideways force pulls it toward the
  axis, and it creeps toward the parent (p. 5).
- **Dispersion is slower than Passey and Melosh predict.** Even for identical fragments separating sideways, the
  final speed is 0.34 of the value Passey and Melosh's formula gives with C = 0.5, about three times smaller
  (section 5.1).

**Applied to our spheres** (parent equal-area diameter D1 = 2.81 m):

| Sphere | r = D2/D1 | θ_inn (deg) | θ_bif (deg) |
|---|---|---|---|
| 10 mm | 0.0036 | 47.4 | 139.4 |
| 100 mm | 0.036 | 51.4 | 134.4 |
| 400 mm | 0.142 | 64.3 | 119.5 |

The table means three regions on the parent:
- **Rear and shoulder region** (θ up to about 47 to 64 degrees): a sphere released here stays trapped in the
  near wake.
- **Sides** (up to about 120 to 140 degrees): it ends in the far wake, still inside the wake.
- **Front cap** (beyond θ_bif): it becomes independent.

For the independent case, t_hat is about 0.32 to 0.45. Using rho = 2813 kg/m³ for both bodies, t_fz is then
about 2 s at 80 km (q_inf about 530 Pa) and about 1 s at 60 km (q_inf about 2,900 Pa). This is an
approximation, because our parent is a hollow shell, not an AA7075 sphere.

**Caveats.**
- The 10 mm and 100 mm spheres (r = 0.0036 and 0.036) lie below the smallest database ratio, 1/16 = 0.0625, so
  applying the model to them is an extrapolation.
- Above about 75 to 80 km a millimetre-to-centimetre sphere is not in continuum flow, which the model requires.
- The model assumes both bodies are spheres of equal density.

**Revised recommendation (supersedes section 0, items 1-3).**
1. **Sideways release velocity: keep zero as the baseline.** Register's dispersion result confirms the speed
   would be small.
2. **Release angle θ matters more than sideways speed.** It decides whether the sphere is trapped behind the
   parent, drifts back in the far wake, rides the shock, or escapes. Add θ as a scenario input: the release
   point on the parent's surface, with the parent's velocity. Record the AFFIx final zone and t_fz in each
   case's summary as an independent check.
3. **Use AFFIx's near-wake drag as a check on our wake-gas model.** Our `wake_flow="both"` drag reduction close
   behind the parent should be compared with AFFIx's force-coefficient ratio of about 0.14. Our crude estimate
   of about one-seventh from the low density alone is of the same order.
4. **Add an explicit near-wake capture regime as a later step.** It is the biggest unmodelled effect. A sphere
   released from the parent's aft half may never leave the plume while the parent stays in continuum flow, so
   its residence would equal the time to parent demise or to the end of the continuum regime.


> **Update, 2026-10-29 (later the same day): full texts read.** Asha saved the PDFs in
> `Literature Review/Wake Growth/`. Laurence and Deiterding (2011), Laurence, Parziale and Deiterding (2012) and
> Passey and Melosh (1980) were read in full from their text layers. The Register et al. PDF
> (`1-s2.0-S001910351930452X-main.pdf`) has corrupted embedded fonts: both MuPDF and Ghostscript render every page
> as scrambled symbols, so it is still unread and needs a fresh download. This section supersedes the
> "abstract only" grading of the two Laurence papers and the recommendation in the original report.

## 0. Primary-source results and revised recommendation

**Laurence and Deiterding (2011), J. Fluid Mech. 676:396-431 (full text).**
- They used inviscid simulations and an analytical model at Mach 25 (perfect gas) for two spheres that start
  touching with no relative velocity. A critical radius ratio r2/r1 separates two outcomes. Here r2 is the
  smaller sphere's radius and r1 the larger's.
- **Below the critical ratio, the small sphere is entrained.** It stays inside the region bounded by the
  large sphere's bow shock.
- **Above the critical ratio, it is expelled.** Near the critical ratio it "surfs" the shock and reaches the
  largest sideways speed.
- **The critical ratio depends on where the small sphere starts** (p. 419-420): about 0.19 when it sits at the
  front of the large sphere, 0.48 to 0.51 at the centre, and 0.67 at the rear.
- **An entrained sphere's sideways speed stays small.** For r2/r1 = 0.4 the peak corresponds to a
  Passey-Melosh constant C of about 0.06, against C of about 0.35 at r2/r1 = 0.5 and 0.04 for equal spheres
  (p. 409). Discussion, p. 427: "the smallest fragments are likely to be entrained by larger fragments".
- **The Passey and Melosh scaling does not capture the effect of radius ratio.**
- **Parent spin lowers the critical ratio.** A spinning parent lowers it to about 0.4, and further to about
  0.2 at about 4 revolutions per second (p. 421).
- **Surfing loads the small body hard.** The authors note that surfing brings very high local pressures and
  heat fluxes on the surfing body.

**Laurence, Parziale and Deiterding (2012), J. Fluid Mech. 713:159-182 (full text).**
- These are free-flight experiments at Mach 4 in a Ludwieg tube, a short-duration supersonic wind tunnel,
  checked against simulations. Final separation velocities agree to within 1 to 11 percent.
- **The experiments confirm entrainment for small radius ratios.** Entrainment is more likely the further
  downstream the small sphere starts.
- **The critical ratios they found** (p. 179): slightly below 0.625 at Mach 4 with the centres aligned, against
  about 0.45 for Mach above 10.
- **Inside the shocked region the sideways force reverses** (p. 175). It becomes attractive, pulling the small
  sphere toward the axis, because the Pitot pressure (the stagnation pressure a body would feel) falls toward
  the axis. The sideways velocity of an entrained sphere therefore peaks and then stays bounded.
- **Smallest size ratio tested:** the experiments went down to r2/r1 = 0.25, and that case was entrained
  immediately.

**Applied to our case.** The parent's equal-area radius is 1.41 m, so our spheres of 1 to 400 mm diameter have
r2/r1 = 0.0004 to 0.14. That is below even the smallest critical ratio (0.19, front-aligned, non-spinning,
Mach 25).
- **Expected outcome: entrainment.** The sphere stays within the region bounded by the parent's bow shock,
  which contains the wake. Its sideways velocity peaks at a small value, and once inside it is pulled toward
  the axis rather than pushed out.
- **Passey and Melosh does not apply here.** Its speeds (about 2 to 24 m/s for C = 1) are for expulsion. The
  entrained sideways velocities are smaller. As an order-of-magnitude guess from the 2011 values, C is about
  0.06 or less, which gives a factor of 4 or more below C = 1 in speed.
- **What remains uncertain.**
  - Size ratios below 0.25 were not tested in experiment.
  - Both papers are for inviscid, perfect-gas, continuum flow, whereas our 80 to 100 km cases are laminar or
    rarefied.
  - Only two bodies were considered, with no parent breaking up around them.
  - A fast-spinning parent could still expel the largest (400 mm) spheres.

**Revised recommendation.**
1. **Keep zero sideways release velocity as the baseline.** It is now supported rather than merely assumed.
2. **Do not use Passey and Melosh as the upper bound.** If a sensitivity case is wanted, use a small sideways
   speed equal to Passey and Melosh with C = 0.06. That is 0.25 times the C = 1 values in the table, for
   example about 2.5 m/s for a 10 mm sphere at 80 km.
3. **Note the opposite effect: entrainment holds the sphere near the parent.** The model's drag-only separation
   may therefore *understate* residence close to the parent rather than overstate it. A future refinement could
   add the attractive sideways force inside the parent's shock layer.
4. **Flag the heating.** Surfing and entrainment bring shock-interaction heating that `droplet_reentry` and the
   sphere model do not represent. This matters for when the sphere starts melting.
5. **Re-download Register et al. (2020)**, which is open access, to check the computational separation database
   at conditions closer to ours.

Literature search, 2026-10-29. Web search and online reading only; nothing downloaded on purpose.

## 1. What was asked and whether it is done

The request was to find published laws or data for the sideways (lateral) separation velocity that a
fragment gains from shock-wave interaction when it is released next to a hypersonic parent body, and to
apply the best law to AA7075 spheres of 10, 100 and 400 mm diameter released beside the 2.8 m parent
spacecraft. The search is done, but with an important limitation: **I could not read the full text of any
of the key papers.** The Laurence papers' PDFs either failed to load (certificate error, an unreadable
compressed PDF, a blocked fetch) or the fetch tool refused to process them. Every source below is
therefore read only at abstract level or through secondary citations, and the strict reading level is
stated for each one. No coefficient below was taken from memory; the only numerical coefficient used
(the Passey and Melosh constant range of 0.03 to 2.25) was quoted by a secondary source that was found
in search results, not read in the original.

## 2. Main finding in plain language

There is one closed-form law, from Passey and Melosh (1980), which says the sideways speed of the smaller
body is the flight speed times the square root of a constant, the radius ratio of large to small body,
and the ratio of air density to body density. Applied to our case with the constant equal to one, it gives
sideways speeds of about 2 to 24 metres per second. However, the later and much better experimental and
numerical work by Laurence and co-workers says that this law does **not** describe unequal-sized bodies:
a small body beside a much larger one tends to be **trapped inside the larger body's bow shock and dragged
along ("entrained")**, rather than expelled sideways, especially when it starts at or behind the large
body's centre. Our size ratios (parent diameter divided by sphere diameter of 7, 28 and 280) are far more
unequal than anything tested, so the Passey and Melosh number should be read as an upper bound, not an
estimate.

## 3. Sources

### 3.1 Passey and Melosh (1980)

- Citation: Passey, Q. R. and Melosh, H. J. (1980), "Effects of atmospheric breakup on crater field
  formation", Icarus 42, 211-233. https://doi.org/10.1016/0019-1035(80)90072-X
- Read: **secondary citation only.** The law and the constant range were taken from a search-result
  summary of a University of Southampton preprint (https://eprints.soton.ac.uk/510415/1/main.pdf, the
  fetch itself was refused with HTTP 403, so even this is second-hand) and from the Laurence, Parziale and
  Deiterding (2012) abstract, which states that it tests this law.
- Law: V_T = V × sqrt( C × (r1 / r2) × (rho_a / rho_m) ), where V_T is the final lateral separation
  speed of the smaller body (metres per second), V is the flight speed of the pair at breakup (metres per
  second), r1 is the radius of the larger body and r2 the radius of the smaller body (any common length
  unit; only the ratio matters), rho_a is the air density (kilograms per cubic metre), rho_m is the density
  of the body material (kilograms per cubic metre), and C is a dimensionless constant.
- Constant: C between 0.03 and 2.25, as reported secondhand; later authors (Popova and Podobnaya 2026,
  Icarus, seen only in a search summary) argue C of about 1 is consistent with observed crater fields.
- Range: derived from dimensional reasoning for two touching spheres that separate purely sideways, and
  calibrated against the spread of terrestrial crater fields (so meteoroid sizes of metres, stony or iron
  densities, speeds of about 10 to 20 kilometres per second). It is not a wind-tunnel law, and the size
  ratio range of its calibration is not known to me.
- Outcomes: it only predicts expulsion; it has no capture or drafting branch.

### 3.2 Laurence, Parziale and Deiterding (2012), "Dynamical separation of spherical bodies in supersonic flow"

- Citation: J. Fluid Mech. 713, 159-182. https://doi.org/10.1017/jfm.2012.453 ; open copy at
  https://authors.library.caltech.edu/records/f1zfh-mg238
- Read: **abstract only** (Caltech repository page and search summaries). The full-text PDF could not be
  read.
- Result: two initially touching nylon spheres, hung on weak threads that break, so they are in **free
  flight** in a Mach 4 Ludwieg tube (a short-duration supersonic wind tunnel), filmed with shadowgraphy at
  13 kHz (13,000 frames per second), compared with inviscid simulations. The final lateral speed of the
  smaller sphere depends strongly on both the radius ratio and the initial alignment angle (the angle of
  the line joining the centres relative to the flow). The Passey and Melosh law is found inadequate for
  unequal bodies because the smaller body is not pushed purely sideways but can be entrained.
- Range: Mach 4 only; nylon in air, so the density ratio is that of a wind tunnel, not of flight; radius
  ratios not read (abstract does not give them, but they are modest, of order 1 to a few).
- Outcomes: greater size disparity, or a starting position of the small sphere downstream of the large
  sphere's centre, makes **entrainment** (capture inside the large body's bow shock) more likely. At a
  critical angle or radius ratio the behaviour switches from entrainment to **expulsion**, and right at
  that switch the final lateral speed is close to its maximum, because of **"surfing"** (the small body
  rides along the large body's shock and keeps being pushed outward for longer).

### 3.3 Laurence and Deiterding (2011), "Shock-wave surfing"

- Citation: J. Fluid Mech. 676, 396-431. Cambridge page:
  https://resolve.cambridge.org/core/journals/journal-of-fluid-mechanics/article/shockwave-surfing/5FB3D50E7AEC25DD0B7BED1268524F6B
- Read: **abstract only**, through search summaries.
- Result: numerical and theoretical study. A downstream body that follows ("surfs") the shock of an upstream
  body can build up a much larger lateral speed than it otherwise would. There is a critical radius ratio
  for initially touching fragments that separates entrainment from expulsion and gives the longest surfing.
  The square-root dependence on radius ratio in Passey and Melosh is said not to describe the separation.
  With several fragments, larger fragments repel one another while smaller fragments are entrained in the
  shocks of larger ones. Surfing plus a modest parent rotation explains much of the spread in separation
  speed inferred from crater fields.
- Range: hypersonic, inviscid simulations; Mach and density ratio not read.
- Companion experiment: Parziale et al., "Experimental investigation of shock wave surfing",
  arXiv:1010.3248, Mach 4; PDF fetched but unreadable, so **title only**.

### 3.4 Register et al. (2019/2020), "Interactions between asteroid fragments during atmospheric entry"

- Citation: Register, P. J., Aftosmis, M. J., Stern, E. C., Brock, J. M., Seltner, P. M., Willems, S.,
  Gülhan, A. and Mathias, D. L., Icarus (online October 2019), ScienceDirect S001910351930452X. Summary page
  read: https://karmaka.de/?p=20596
- Read: **abstract and highlights only.**
- Result: computational fluid dynamics of two interacting spheres produced a lookup table of forces;
  trajectories from that table were checked against free-flight hypersonic wind-tunnel tests (the DLR H2K
  tunnel at Cologne, Gülhan's group); the final relative states are categorised and the interaction time
  estimated, and the result is given as semi-analytic equations, non-dimensionalised for general use. The
  highlights say separation speeds are **smaller than previously computed** (that is, smaller than Passey
  and Melosh type estimates). The coefficients were not visible.
- This is the most directly usable source if its full text can be obtained, because it gives explicit
  equations and covers free flight.

### 3.5 Other leads

- Artemieva and Shuvalov (1996), "Interaction of shock waves during the passage of a disrupted meteoroid
  through the atmosphere", Shock Waves 5, 359-367: **citation only** (from a reference list). The 2001
  paper and its coefficient were not found. A 2025 arXiv paper on cometary entry (https://arxiv.org/pdf/2503.22632,
  search summary only) writes the separation speed as v_s = sqrt(C_V rho_atm / rho_m) × v, attributing it
  to Passey and Melosh and Artemieva and Shuvalov, without the radius ratio.
- Whalen et al. (AIAA 2020, Southampton, https://eprints.soton.ac.uk/441565/) describe the Wheeler et al.
  fragment-cloud model as the state of the art and say it is seeded with the Register database: **search
  summary only.**
- Marwege and Gülhan: no paper by Marwege on fragment separation was found. Related DLR H2K free-flight
  multi-body work exists (Seltner, Willems and Gülhan, on cube aerodynamics) and an ESA test case
  (https://indico.esa.int/event/403/page/620-test-case-2022) on a ring crossing a cylinder's shock, all seen
  only in search summaries.
- Park and co-workers on fragment pairs, and "collective wake" / "pancake" drafting: nothing specific was
  found in the time available.

## 4. Applying Passey and Melosh to our case

The inputs are: parent diameter 2.8 m (so r1/r2 equals parent diameter over sphere diameter), sphere
density rho_m of 2813 kilograms per cubic metre, and C equal to 1. At 80 km the speed is 7.6 kilometres
per second and the air density 1.85e-5 kilograms per cubic metre, so rho_a/rho_m is about 6.6e-9. At 60 km
the speed is 4.3 kilometres per second and the air density 3.1e-4 kilograms per cubic metre, so
rho_a/rho_m is about 1.1e-7. The table compares the resulting sideways speed V_T, in metres per second,
for each sphere and altitude.

| Sphere diameter | Size ratio r1/r2 | V_T at 80 km (m/s) | V_T at 60 km (m/s) |
|---|---|---|---|
| 10 mm | 280 | about 10 | about 24 |
| 100 mm | 28 | about 3.3 | about 7.6 |
| 400 mm | 7 | about 1.6 | about 3.8 |

The conclusion to draw is that the predicted sideways speed is between about 1.6 and 24 metres per second,
largest for the smallest sphere and at the lower, denser altitude. Across the reported constant range
(C from 0.03 to 2.25) each value can be multiplied by anything from about 0.17 to 1.5, so the full spread is
about 0.3 to 36 metres per second.

To compare with the plume, which is a few metres in radius near the parent and about 6 m in radius at 1 km
behind: a sideways speed of 10 metres per second moves the sphere 3 m in 0.3 seconds, which is enough to
leave the near-parent plume quickly, whereas 1.6 metres per second takes about 2 seconds to cover 3 m. Over
the time the sphere needs to fall 1 km behind the parent (which depends on its drag difference and is not
computed here), even a few metres per second could carry it outside a 6 m plume. So if the Passey and
Melosh numbers were right, the sideways speed would matter for residence time, especially for the 10 mm
sphere; for the 400 mm sphere at 80 km it would be marginal.

## 5. What it means

The physically better sources say our case is in the regime where the law fails: a sphere 28 to 280 times
smaller than its neighbour is very likely to be entrained in the parent's shock layer if released at or
behind the parent's centre, which pushes the sphere back toward the wake rather than out of it. Expulsion
with surfing happens only near a critical angle and radius ratio that I could not read and which was
measured only at Mach 4 for modest ratios. Register et al. also say separation speeds are smaller than the
older estimates. So the table above is an upper bound for an expulsion case, and the real answer for most
release positions may be close to zero sideways speed plus entrainment.

## 6. What remains unknown, and recommendation

Unknown: the critical angle and radius ratio from the 2011 and 2012 papers; the Register et al. equations
and coefficients; whether any of this holds at size ratios above about 10 or at Mach 5 to 27 with real-gas
effects; and where on the parent a real fragment is released.

Recommendation: add the sideways release velocity as an **optional sensitivity parameter**, not a default.
Use Passey and Melosh with C equal to 1 as an upper bound (about 2 to 24 metres per second here), and run
zero as the baseline, since entrainment is the more likely outcome for our size ratios. The uncertainty is
at least a factor of about six either way from C alone, and the law's form itself is disputed for unequal
bodies. The best next step is to obtain the full texts of Register et al. (2020) and Laurence et al. (2012)
through the MIT library and replace the bound with their equations and the entrainment criterion.
