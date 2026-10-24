
## 9. Melting and the melt film (Step 3)

- **Latent heat, enthalpy method.** The nodal enthalpy h(T) = ∫c_p dT + L_f f_l(T) carries the latent heat; f_l is
  linear between the solidus and the liquidus (`AA7075_range`: 750–908 K, ASM Handbook) or a ±2 K numerical ramp
  around DRAMA's single 850 K (`AA7075`); L_f = 400 kJ/kg and the solid tables are DRAMA's `drama-AA7075` verbatim.
  The capacity matrix is lumped (diag Σ φ_e V_e/4 ρ c_i), so the discrete energy is exactly ρ∫h dV with h
  interpolated linearly, and each Newton step is mapped through the true h(T) per node so that a node cannot jump
  across the melting range (the element-mean formulation of Step 2 could not carry a 4 K ramp under a 30 K nodal
  spread). The Stefan front on a box is reproduced within 0.3 %.
- **Liquid properties** (DRAMA has none): ρ_l 2400 kg/m³, μ_l 1.3 mPa s, Σ 0.86 N/m — pure aluminium near the
  liquidus (Smithells; Assael et al. 2006; ASM Vol. 2); alloy corrections within 10 %; no oxide skin. They are
  constants: the film now has a temperature (below), but its viscosity and surface tension are still evaluated at the
  liquidus, which is where most of it sits. μ_l falls by roughly a third per 200 K of superheat, so the runoff and the
  thin-branch wavelength are the quantities a temperature-dependent μ_l would move.
- **How deep the shear reaches, and what may therefore leave.** The gas shear penetrates the liquid only to Girin's
  conjugate depth, measured at 215 to 303 µm on the 100 mm flight, while the contiguous molten layer beneath the wall
  is 350 to 890 µm. Two rules follow. Melt leaves an element only if the shear can reach it: a wall-owning element
  always, a buried element only within the conjugate depth of the nearest patch, and where the wall-Knudsen gate
  denies Girin's closure, only wall-owning elements. And which of Girin's regimes a patch takes is decided — at his
  first stage — on the depth of *contiguous* liquid under the patch: the march inward stops at the first element that
  is not fully molten, so melt blocked by solid is never credited, while the film account remains what can actually
  move and be stripped. The layer decides that regime choice and the Rayleigh–Taylor depth criterion and nothing else;
  Girin's surface Weber number, wavelength, period and droplet radius are all on the conjugate depth, as he has them.
- **Runoff and stripping happen at once, and that is not double counting.** Wall shear stress is a flux, not a stock:
  the shear sets up a mean flow in the sheared sublayer, which carries mass along the wall, and the instability is a
  perturbation on that flow's free surface, which removes mass from it. The only coupling is the momentum the departing
  droplets carry, of relative size (melting speed × conjugate depth / kinematic viscosity) — 5 % at a recession speed
  of 0.1 mm/s, 48 % at 1 mm/s. What is *not* yet done, and is deferred with the pressure response of the molten
  region: the pressure-gradient and deceleration-driven part of the runoff flux should integrate over the whole liquid
  depth while the shear-driven part stays limited to the conjugate depth.
- **What flows.** Melt and runoff start at the liquidus: an element's material becomes film in proportion to a ±2 K
  ramp at the liquidus (f_feed), so the mushy range holds latent heat but neither runs off nor is stripped
  (conservative; a coherency-point treatment is a future iteration). What leaves an element is its *molten* part, at
  the enthalpy that part carries (the f_feed-weighted mean of its nodal enthalpies), and what stays behind is the
  colder rest — so the element is debited the difference, which is what makes the melt front advance at the
  energy-limited rate whatever the element size. With `--removal instant` the mass leaves the body at the liquidus
  enthalpy instead, and the difference is returned to the element's nodes.
- **The film's temperature and re-solidification.** The film is thermally thin (q b/k_l = 0.22 K across a 10 µm film
  at 2 MW/m², 22 K even at 1 mm; b²/α = 0.3 ms against the 0.5 s macro step), so it is given no energy equation: its
  mass is carried on the boundary nodes with the liquid heat capacity and it is at the surface temperature by
  construction. It therefore heats, cools and freezes with the surface, and droplets leave carrying the surface's
  superheat (+2.4 % on h_liquid for the 50 mm flight, +0.6 % for the 100 mm one) rather than at the liquidus exactly. Because a film is liquid by
  construction it holds the liquid enthalpy h(T) + L_f(1 − f_l) wherever it sits, so feeding it costs the latent heat
  the mass has not yet paid and no mass can be relabelled from solid to liquid for free. Re-solidification is the
  mirror of the feed rule — the fraction 1 − f_feed(T_patch) of a patch's film returns to its owner element as solid,
  and the two directions are netted into one transfer per element so that the film is never churned. Two limits are
  declared: φ_e is capped at 1, so film whose owner element is full (or has died) stays liquid for good — the mesh
  cannot grow a crust outside itself — and `film_frozen_fraction` reports how much film is in that state; and a thick
  crust would conduct, which a lumped nodal capacity does not represent. A re-solidified crust is therefore still
  available to run off and be stripped, which is conservative for mass loss.
- **How melt energy is booked.** Mass moving at one temperature books nothing; mass that arrives somewhere hotter or
  colder books the enthalpy difference it carries, on the nodes it arrives at. No deferred load may move a node more
  than 1000 K in one macro step (a node whose own mass has melted away has nothing to heat with it); the remainder
  waits and is reported as `unapplied_load_J`, and is counted as dropped if the node dies first. The coupled energy
  balance is exact (1e-9 of the absorbed heat) with the queue and the dropped loads in it.
- **Element fractions and death.** φ_e scales an element's heat capacity, not its conductivity (a thinner sliver of
  the same material conducts better, not worse; scaling k isolated the surface nodes and drove them to thousands of
  kelvin). A patch owner dies at φ_e ≤ 5 % (its remainder joins the film); interior elements keep ≥ 10⁻³ so that no
  cavity opens; nodes without material are pinned. The exposed faces of a dead element become patches (a transient
  pit until its neighbours die: with Lees' distribution the pit walls at θ ≈ 90° receive almost no heat, so the
  heating artefact is small; the SESAM-equivalent mode is defined by its total and is unaffected).
- **Mesh.** Four prism layers (0.25, 0.5, 1, 2 mm) built by radial projection of the inner 2 mm/8 mm gmsh sphere's
  boundary triangulation (exact on a sphere; each prism split into three tetrahedra by the smallest-node-id diagonal
  rule); the surface triangles are 2.16 mm (h_surface × R/(R − 3.75 mm)). 46 k nodes / 255 k tetrahedra for the
  100 mm sphere; 0.6 s of conduction per macro step, inside a melt step of about 2 s in all.
- **Gas-side surface flow.** Wall pressure from modified Newtonian p_w = p∞ + (p_s − p∞)cos²θ up to the sonic point
  φ* (40.7° at γ_PM 1.15, 43.4° at 1.4) and a Prandtl–Meyer expansion through the turn beyond it, blended over 10°
  above φ* and floored at p∞: pure modified Newtonian gives p_w = p∞ at 90°, a factor 28–48 below the measured
  sphere pressure, in exactly the band where most of the spraying happens. Prandtl–Meyer supplies only p_w — the edge
  state is the isentropic expansion of the equilibrium stagnation state to it (Cantera, 1° bins), so the perfect-gas
  Mach number never propagates into the state variables. Leeward p_w = p∞, no shear.
  Girin's linear-profile boundary layer in Ranger's (1972) form with the actual edge velocity, δ_a² = 58.1 ν_e ∫u_e⁴ds/u_e⁵
  (exact for potential flow; Thwaites' momentum thickness is a constant 12.3 × smaller within ±1.7 %); shear
  τ_c = μ_e u_e/δ_a. The flow regime is decided at body scale rather than on the boundary layer, because a continuum
  construction must not be allowed to certify itself: a distinct bow shock is assumed only below Kn_body 0.01 (SESAM's
  own mean free path over the body's current diameter), the flow counts as merged from there to Kn_body 10, and
  free-molecular above it or below Ma 1. Inside the shock-layer branch a wall Knudsen number Kn_local = λ_w/R,
  evaluated at the wall state (p_w, T_wall), selects the melt *closure* and never whether spraying happens — Girin's
  conjugate layers below 0.01, a Couette film velocity above — since his dispersion relation contains no gas
  parameters. The shear is the slip-corrected τ_c below Kn_local 0.1 and the free-molecular τ_fm = ρ∞V² sinθ cosθ
  above it, with `--rarefied-shear bridged` blending the two by SESAM's f(Kn) instead; the merged branch always
  bridges. The gate is a declared conservatism whose cost is measured rather than assumed: the wall gas is cold and
  compressed, so Kn_local ≈ Kn_body/84, and every step records `kn_body`, `kn_local_stag`, `re_shock`, `flow_branch`
  and the closure area fractions. With the Prandtl–Meyer expansion resolving the wall pressure the rim carries
  ρ_e ≈ 1e-4 kg/m³ and Kn_local ≈ 2e-3, so the "rarefied rim" of the earlier design — and the film pile-up and
  droplet-size cap it forced — was an artefact of pure modified Newtonian and is gone. The driving gradient
  G = 2(p_s − p∞) sinθ cosθ/R − ρ_l a sinθ includes the body's deceleration (the film is pushed toward the nose).
- **A film deeper than its patch is wide is not a film.** On a collapsing body the melt of thousands of dying patches
  is concentrated onto a few small survivors — by the death hand-over, by the interior-element feed to the nearest
  patches, and by a runoff graph that has no edge out of the last windward ring — and m_f/(ρ_l A) stops being a depth:
  in the last 5 % of the 50 mm flight a 0.5–0.8 mm facet carries melt that would be a 4–14 mm blob. The reported
  thickness therefore covers only facets with b ≤ √A, and the rest is reported as `film_blob_fraction` (nonzero in 17
  of 408 steps on the 50 mm flight, up to 0.999 of the film in the collapse; in 110 of 217 on the 100 mm one, median
  0.131 — the nose crater, where τ and G ∝ sinθ leave the film nothing to drive it out — and at most 0.21 % of m₀). The mass and the energy stay exactly accounted; what is lost is
  the lubrication picture, and with it the meaning of the runoff and spraying branches on those facets. A patch graph
  that merges slivers, or a film carried per element rather than per facet, is the fix; it is not in this iteration.
- **Film and runoff.** One film mass per patch; lubrication velocity and flux with a thin (b ≤ δ_m) and a thick
  (b > δ_m) branch; runoff by a linearly implicit upwind scheme on the patch graph (4 sub-steps per macro step,
  exact conservation, exact steady state; a wetting front advances one patch per sub-step), never across the equator;
  leeward films are static and can stay attached (their fate is Step 4's). The film is fed where it melts and stripped
  there within the step, so it stays microns thin except where the runoff piles it at the windward rim.
- **Spraying: Girin's three regimes.** His (2017) classification is applied in two stages. δ_m is the melt velocity
  boundary layer his Eq. (2) *predicts*, so the liquid depth decides whether that layer can form at all: where it cannot
  (layer ≤ δ_m) the rigid core still stabilises the disturbances and the case is **regime 1**, Girin & Kopyt's (1994)
  side mode — λ* = 1.5 M_e Σ/(ρ_e u_e²), r = λ*/4, τ* = 2 capillary periods, rate ρ_l min(b, λ*/8)/τ*, their Table 1
  mass rate. Where it has formed, the profile is the conjugated pair and the kinematic viscosities pick the mechanism:
  ν_melt > ν_gas gives a near-discontinuous profile and classical Kelvin–Helmholtz (**regime 2**, Girin & Kopyt's λ* and
  τ* in their deep-film limit with the torus rate; inert for liquid aluminium, which sits 3.3×10⁴ to 6.9×10⁶ from that
  threshold), while ν_gas > ν_melt gives the inflated profile and his **regime 3** gradient instability — the dispersion
  relation solved numerically (Δ_f, Im Ω_f tabulated against We_s; k_r 0.17, k_t 1.1, We_cr 4.62, his constants for
  ordinary liquids), δ_m and V_s from his conjugated-layer relations, stripping ρ_l π r²/(λ_f t_per) per area. The
  rarefied branch uses regime 1's mode on the freestream Mach number and momentum flux (an extrapolation). Two gates
  apply to every regime: Girin's We_s > We_cr, his critical angle evaluated against the model's own local flow, and the
  requirement that one whole wavelength fit inside the contiguous molten region the patch belongs to (the region, not the
  mesh facet). Girin & Kopyt's front-surface Rayleigh–Taylor mode is **applied as well** (`--rt-spray on|off`, default
  on): the body's deceleration is 6–95 m/s², not the 10–30 m/s² first assumed, and the mode needs both its depth
  criterion W cos φ h² ρ_l > 3 Σ, which confines it to the nose, and the laterally bounded form, which supplies the
  admissible wavelength and growth rate. It releases 6.27 % of the sprayed mass on the 100 mm sphere, out to 20° and
  exactly zero beyond, as ~1.7 mm droplets against the shear population's ~190 µm; on the 50 mm sphere it barely fires,
  because its threshold is an absolute depth (3.36 mm at 95 m/s²) a smaller pool never reaches. Where two modes could
  claim one interface the shorter growth time takes the patch. Droplets are capped at the film on the patch and a
  quarter of the body radius; one radius per patch and step; no within-patch size spread; recorded at birth, not
  tracked. Every mode strips far faster than the melt supply, so the mass loss is energy-limited (Girin's "outstripping
  ablation") and the instabilities set the droplet size and the release map rather than the mass.
- **Size feedback** (decided 2026-09-21, replacing the spec's fixed R₀). Mass = Σφ_e ρV_e + film; the drag reference
  area is the current surface's projection; the body Knudsen number uses the equivalent diameter of the remaining
  mass; the stagnation radius for the heating and the surface flow is a least-squares sphere fitted to the current
  windward cap (the patches within (1 − cos 30°) R_t of the front-most point, R_t the transverse radius about the mass
  centre), bounded to [0.1, 1.67] × R_t — the front erodes fastest and flattens, so the fitted radius grows from R₀ to
  the cap and the stagnation heating falls as R^−½ (a flat face heats like a sphere of 1.67 × its radius: its
  stagnation velocity gradient is ≈ 0.6 × a sphere's, Boison & Curtiss 1959). The sphere drag tables are kept. The
  verification devices keep D₀ and R₀ (SESAM's convention).
- **Verification devices.** `--removal instant` (every element's liquid leaves as it forms) and `--k-scale 1e4`
  reproduce SESAM's lumped Q/L_f melting (mass within 1.3 % of m₀, onset within 0.2 km, end within 1.3 %); Girin's
  Table 1 is reproduced in its exact tier (GI, φ_cr) and, with the ambient-density Reynolds number, in t_f, N and
  r_med within 30 %; his spraying durations are not (half his, and far less for the stony variant).
