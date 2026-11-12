"""Melt spraying: Girin's gradient instability per patch and the droplet release bookkeeping (spec Step 3 section 9).

Branches per windward patch with film (regime from surface_flow, b from the film, delta_m from Girin's Eq. 2):
  thick, continuum/slip (b > delta_m): delta_m = (alpha/mu^2)^(1/3) delta_a, V_s = (alpha mu)^(1/3)/(1 + (alpha
      mu)^(1/3)) u_eff, We_s = rho_l V_s^2 delta_m / Sigma; unstable when We_s > We_cr; lambda_f = 2 pi delta_m /
      Delta_f(We_s), r = k_r lambda_f, t_per = k_t delta_m / (V_s Im Omega_f(We_s)), stripping rate per area
      mdot = rho_l pi r^2 / (lambda_f t_per) (one torus of cross-section r per wavelength per period, Girin 2017
      Eqs. 5-8 and Appendix D).
  thin, continuum/slip (b <= delta_m): Girin & Kopyt 1994 side mode with the edge state: lambda* = 1.5 M_e Sigma /
      (rho_e u_eff^2) (their 1.5 M d / We_d, the film thickness cancels), r = lambda*/4, growth time tau* =
      2 tau_v(lambda*) = 0.798 lambda*^1.5 (rho_l/Sigma)^1/2 (their Eq. 12: twice the capillary-wave period),
      release rate mdot = rho_l min(b, lambda*/8) / tau* (their Table 1's mass rate is rho_1 r_d / (2 tau_d) =
      rho_1 lambda*/(8 tau*), reproduced 2026-09-20; the film supplies at most its thickness). Active when b >= B_MIN.
      (Their dissipation cut-off lambda_t = lambda*/3 is always below lambda*, so it never limits the mode.)
  regime 2, continuum/slip (liquid layer > delta_m, as for the thick branch, *and* nu_melt > nu_gas): the
      Kelvin-Helmholtz corner of Girin's own two-regime classification -- a deep film whose melt is the more viscous
      medium, where the interfacial profile is a tangential discontinuity rather than a resolved gradient. Wavelength and
      growth time are then Girin & Kopyt's pair, lambda* = 1.5 M_e Sigma/(rho_e u_eff^2) and
      tau* = 0.798 lambda*^1.5 (rho_l/Sigma)^1/2 (their Eqs. 11 and 12), with the rigid-wall factor cth(Lambda) of their
      Eq. (5) -> 1 because the film is deeper than the wave sees, so the unmelted body no longer stiffens the response.
      The release is the torus form Girin (2017) uses on the thick branch, mdot = rho_l pi r^2/(lambda* tau*) with
      r = k_r lambda*, and *not* the thin branch's rho_l min(b, lambda*/8)/tau*: that factor is there because a film
      shallower than the wave can give up no more than its own depth, which is not this film's position. Like the thick
      branch's rate it therefore does not reference b at all -- and it is 0.73 x the thin branch's rate at the same
      wavelength, not a larger one, wherever that branch is not itself depth-limited; what bounds the release here is the
      film mass, through dm = min(mdot A dt, m_f) as on every branch.
      PROVENANCE: both tests are Girin's, and the numbering is his. Girin (2017) Sect. 1 sets out "two regimes of
      meteoroid ablation due to instability". In *dominant ablation* the mass loss overtakes fusion, the molten layer is
      thinner than the melt velocity boundary layer, b < delta_m, the rigid core still stabilises the disturbances, and
      "the results obtained in Girin & Kopyt (1994) are thus valid" -- this module's thin branch. In *dominant fusion*,
      b > delta_m, the profile is the completed boundary layer: free conjugated layers, delta_a in air and delta_m in the
      melt, "independent from the solid core". His Sect. 2 then selects the mechanism inside that conjugated picture on
      the kinematic viscosities. The "classical Kelvin-Helmholtz type ... is in action only when the liquid kinematic
      viscosity is greater than that one of the gas: nu_m > nu_a", because then "the liquid boundary layer is much
      thicker than the gas one" and V_s << V_a, so "the velocity profile is close to the discontinuous one (tangential
      discontinuity) which is the base for the Kelvin-Helmholtz mechanism". For the inverse inequality nu_a > nu_m the
      thicknesses become comparable, delta_m = O(delta_a) and V_s = O(V_a), the profile is "inflated", and the mechanism
      is "completely different from the Kelvin-Helmholtz one", being set by his Eq. (1) -- the gradient instability the
      thick branch solves through `dispersion`. So regime 2 is his Kelvin-Helmholtz cell and the thick branch is his
      gradient cell of the same two-by-two; what this model supplies is only the arithmetic for the cell he does not
      develop, by reusing the 1994 mode that belongs to that profile.
      UNREACHABLE FOR THIS PROBLEM, measured 2026-09-25, and consistent with Girin's own description rather than in
      tension with it. (a) Liquid aluminium is nowhere near the Kelvin-Helmholtz cell: nu_melt = 5.42e-7 m2/s against an
      edge nu_gas measured at 0.0178-2.04 m2/s over the 100 mm flight to 110 s and 0.0951-3.75 m2/s over the whole 50 mm
      flight, so nu_gas/nu_melt runs 3.3e4 to 6.9e6 and never approaches one. This body sits firmly in Girin's
      nu_a > nu_m case, which is exactly why his gradient instability is the right model for it. Directly: 0 of 533 960
      windward wet patch-steps (100 mm) and 0 of 42 268 (50 mm) take this branch. (b) The two tests reinforce each other
      rather than crossing, because Eq. (2) carries the same ratio:
      delta_m/delta_a = (alpha/mu^2)^(1/3) = ((nu_melt/nu_gas)^2 rho_l/rho_e)^(1/3), so at nu_melt = nu_gas the conjugate
      depth is already (rho_l/rho_e)^(1/3) ~ 151 air boundary-layer thicknesses, about a metre at delta_a = 7 mm -- which
      is Girin's "much thicker than the gas one" restated numerically. A melt viscous enough to be Kelvin-Helmholtz
      therefore carries a delta_m no film here can exceed, so the gate puts it at b < delta_m, which is Girin's
      dominant-ablation case and his Girin & Kopyt (1994) mode: the thin branch. That is the correct routing and not a
      gap. Reaching regime 2 needs a viscous melt *and* a melt layer deeper than a metre-scale delta_m, which this
      problem never produces; the branch is kept because it is where that cell belongs if the material or the altitude
      band changes, and it is exercised by test with delta_m supplied directly.
      Where it would show up, checked rather than assumed: the source table's `branch` column (`SOURCE_COLUMNS` index 8)
      and nowhere else. The surface VTK carries surface_flow's `closure`, not the spray branch, so a regime-2 patch is
      not distinguishable there; no history column is added; and `thick_branch_fraction` is `film.lubrication`'s own
      thick flag -- isfinite(delta_m) & (layer > delta_m), the whole deep set -- so it would count a regime-2 patch as
      thick. That column measures how deep the liquid is, not which instability took the patch.
Two gates (2026-09-24; the first on the thick branch only since 2026-10-22):
  * Girin's (2017) stability criterion We_s > We_cr, with We_s = rho_l V_s^2 min(delta_m, layer) / Sigma -- his
    closed-form critical angle phi_cr evaluated against the model's own local flow, which holds on an eroded body where
    a single angle does not. On the thick branch this is his We_s = rho_l V_s^2 delta_m / Sigma unchanged; on a film
    shallower than the conjugate layer it is measured over the film's own depth. The thin branch had no threshold
    before: Girin & Kopyt's inviscid side mode has none of its own (it is unstable at every wavelength, their
    "unlimiting potential intensity of high-frequency disturbances"), so it stripped film at any angle however small
    the shear -- 2.2 g inside 10 deg on the 100 mm flight, from a wave 6 x wider than its own facet. Since 2026-10-22
    the criterion gates the thick branch only: Girin & Kopyt's Eq. (5), solved, has no critical Weber number, and the
    threshold it has once the wavelength is bounded is the wave fitting in the molten region, which `wave_fits` tests.
    Their three branches (thin, rarefied, regime 2) take `viscous_film` instead: the dominant wave's viscous damping
    on a film of the non-rigid depth over a rigid floor -- Lamb's interior term plus the floor's Stokes layer -- must
    not outpace its growth 1/tau*, and a film shallower than that Stokes layer is outside the inviscid mode's validity
    (recorded, not sprayed).
  * `wave_fits`: one whole wavelength must fit inside the *contiguous molten region* the patch belongs to, not on the
    patch (a mesh facet is not a physical boundary and the melt is continuous across it).
KNOWN LIMITATION, deferred: Girin & Kopyt derive the side mode for the side surface (their sin Theta -> 0,
V_0 -> V_inf) and assign the *front* surface (V_0 -> 0, sin Theta -> 1) to a different, aperiodic Rayleigh-Taylor
solution, their Eq. (14). The branch selection here has no angular dependence -- it turns only on layer vs delta_m --
so the side-surface asymptotics are applied from the stagnation point to the equator, including where the paper's own
premise fails. Implementing their front/side split (sin Theta = cos phi on this body) is future work; until then the
near-nose release rests on a solution used outside its stated domain, which is 0.2 % of the sprayed mass.
  free-molecular branch, thin film: the thin mode with the freestream momentum flux, lambda* = 1.5 M_inf Sigma /
      (rho_inf V^2), and the Couette film velocity from tau_fm -- an extrapolation of a continuum film theory
      (spec section 17.3), sensitive to --rarefied-shear.
The branch selection follows the closure gate of surface_flow, not a regime of its own: Girin's dispersion relation
contains no gas parameters (it is a melt-side instability), so spraying is never switched off by rarefaction -- what
the gate decides is which closure supplies the melt velocity scale, and hence whether delta_m exists to define a
"thick" film at all (spec section 18).
The 1994 front-surface Rayleigh-Taylor criterion W b^2 rho_l > 3 Sigma is evaluated and reported (lambda* =
2 pi (3 Sigma/(W rho_l))^1/2, tau* = (27 Sigma/(4 W^3 rho_l))^1/4), never applied. W is the deceleration *normal to the
film* -- their W sin(Theta), which on this body is W cos(phi): the whole deceleration at the stagnation point, zero at
the equator where it lies in the surface (corrected 2026-09-24; the whole deceleration over-drove every patch off the
nose by 1/cos(phi)).

Release per patch and step: dm = min(mdot A dt, m_f) (the film after feed and runoff), dn = dm / (4/3 pi rho_l r^3);
one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX]. The layered
skins' records are one row per column, per mode, per macro step (`skinbody.SprayPlan.source_block`, sub-plan 22 Task 2):
the same columns, written by the same `source_block`, the droplets' solid fraction and the enthalpy per kilogram they
carry (spec 8 and 10.2; the 2026-10-07 design 7.5) and, last, their release angle as Register et al. (2020) define it
(`register_angle`, sub-plan 22 Task 4). A skin flight at the continuum step makes tens of millions of records, so the
command line has the skin body write each macro step's block to disk as it comes (`SourceSpool`), and the coupled writer
reads either form column by column (`source_columns`), the numbers the same to the last bit.

The release angle (Asha, 2026-10-29; Register, P. J., et al. 2020, "Interactions between asteroid fragments during
atmospheric entry", Icarus 337, 113468, section 4.2): the angle between the larger body's trajectory axis and the line
joining its centroid to the released fragment's, measured from downstream, theta_R = 180 degrees - angle(v_hat, p - c):
180 at the stagnation point, 90 at the equator, 0 directly behind. It is computed from the actual vector p - c, so a
flattened or receded body gets its own angle (the source table's `theta_deg` is the surface normal's angle to the flight
direction, a different thing off a sphere). p is the row's position; c is the centroid of the volume the body's outer
surface encloses (`enclosed_centroid`), from the same surface and at the same instant as the positions the rows carry --
the element model's surface at its spray step, the skins' top surface where their plan was evaluated, the macro step's
start -- and v_hat the body's direction of motion through the air, its fixed attitude.

The droplet's size caps (2026-10-24, Asha). Every branch's radius is capped by a quarter of the body radius -- except the
Rayleigh-Taylor droplet on the crest-limit path under `spray_gates` `2026-11-10`, which only the crest limit caps (Asha,
2026-11-10, sub-plan 22 Task 5h) -- and by one of two volumes. Without `cap_depth` (the element model) it is the film on the patch, r <= (3 m_f / (4 pi rho_l))^1/3,
which depends on the mesh facet's area: measured on the element model's 100 mm flight with the thin gate it bound on
67 % of the thick branch's release rows (24 % of its mass), 16 % of the thin branch's (4.4 %) and set essentially every
Rayleigh-Taylor droplet (about 2 mm). With `cap_depth` h (the layered skins, where there is no film account) a droplet
holds no more liquid than one crest of its own wave can give, (4/3) pi r^3 <= lambda^2 h, lambda the wavelength of the
mode that sprays the patch -- the thick branch's 2 pi delta_m / Delta_f, Girin & Kopyt's lambda* on the thin, rarefied
and regime-2 branches, the admitted Rayleigh-Taylor wavelength -- and h the depth the release can take: lambda^2
because a crest repeats every wavelength in both surface directions, so the patch's area does not enter. A depth of
nothing gives no droplet, and nothing leaves. On the thick and thin branches it binds only on layers shallower than
about 2 % (thick: k_r = 0.17) and 6.5 % (thin: lambda*/4) of the wavelength; on the Rayleigh-Taylor branch, whose
droplet is its wavelength (centimetres), on any layer shallower than about four wavelengths, so always: at lambda = 2 cm
and h = 1 mm it is 4.6 mm, against about 2 mm under the film cap.

The thin branch's Knudsen gate and the slip band (sub-plan 22 Task 5i; Asha, 2026-11-06 and 2026-11-10). Girin's Eq. (2)
and Ranger's boundary layer are continuum results, so `melt_layer` gives a conjugate depth only on Girin's closure (body
Knudsen number below 0.01, wall Knudsen number below 0.01). Under `spray_gates` `2026-11-10` the thin and rarefied branches
shear-spray only where a conjugate depth exists, whichever way it is found: elsewhere (the merged branch, the free-molecular
one, a patch off Girin's closure) a wet patch takes no shear mode, because the thin branch sprays only where Girin's continuum
closure applies, or where the band's momentum balance gives a depth, not because the edge state it reads is invalid (the
review of Task 5i, Minor 5, 2026-11-10); the Rayleigh-Taylor mode is not gated. With `slip_band` on, a patch has a depth in the band
too: where the body is on the merged branch with a body Knudsen number below 0.1 and the patch's own wall Knudsen number is
below 0.1, `spray_layer` gives `momentum_layer`'s depth, delta_m^3 = 3 mu^2 I / (rho tau^2 r) from the bridged shear tau and
its integral I from the nose (`surface_flow.shear_moment`), the liquid's properties and the nose radius the integral was taken at (the flow's `nose_radius`, Task 5j: the
radius `nose_radius()` gave the surface flow, the fitted one under --size-feedback current). ASSUMPTIONS the band carries: I assumes a sphere (axisymmetric, r = R sin(theta)),
an approximation on an eroded body; the whole molten depth may spray, with liquid properties, as below 0.01, the
liquid-slurry viscosity contrast being future work; delta_m jumps by about 2 to 2.5 times at body Knudsen number 0.01, where
Eq. (2) hands over to the balance; the thin branch keeps its edge-state wavelength, the bridge's continuum endpoint, which
the flag accepts knowingly; the thick branch has no viscous gate, so a slurry layer at a Reynolds number near 3 can spray by
inviscid theory with `--slurry-spray on`; the droplet constant k_r was calibrated on Ranger's depth, so the band's droplet
sizes carry that; records still label the band's patches as Couette closure. Only the spraying reads the band's depth: the
feed and the film's transport keep `melt_layer`.

The profile path (`SprayModel.dispersion` "profile"; the kinked-profile dispersion design of 2026-11-11, sections 6.5-6.8).
Where the caller has resolved a column's own melt speed profile (`ProfileModes`, built for the layered skins by
`melt_profile`), `evaluate` takes the column out of the thick, regime-2, thin and rarefied branches above and gives it two
shear modes, only where delta_m is finite (decision 6: delta_m decides whether a column may shear-spray, the profile how):
Girin's gradient mode solved on that profile above its wall, branch 5 (`BRANCH_PROFILE`): lambda = 2 pi / k, r = k_r
lambda, growth time k_t / growth and his torus rate, active where the growth beats its own wave's viscous damping
2 nu_eff k^2, with no critical Weber number; and the gas-driven mode, the thin (or rarefied) branch with the depth to the
wall in place of b and nu_eff, the layer's Couette viscosity, in `viscous_film`. The shorter growth time takes the column.
Its droplet's crest cap reads the depth of the run of the mode that took it, and its records carry nine more columns
(`PROFILE_COLUMNS`), while `SOURCE_COLUMNS` stays the element model's 27."""
import os
from dataclasses import dataclass, field

import numpy as np

from . import dispersion
# by name, because SprayModel's `dispersion` setting shadows the module inside the class
from .dispersion import DISPERSION_NAMES, WE_CRITICAL_THEORY, DispersionTable
from .surface_flow import BRANCH_FREE_MOLECULAR, CLOSURE_GIRIN, KN_LOCAL_SLIP

K_R, K_T = 0.17, 1.1
B_MIN = 1.0e-6                     # m, films thinner than this do not spray (numerical floor)
# Size bins, ten to a decade: 1 um to 0.1 m (sub-plan 22 Task 5j, the controller's ruling of 2026-11-10). They stopped at 10 mm
# (40 bins) until then, and the Rayleigh-Taylor droplets, whose wavelength is centimetres on the layered skins, fell above them
# and were not counted; the 40 bins' edges and counts are unchanged, 10 bins were added above them (`size_distribution.csv`
# has ten more rows per window, empty for a run in which no droplet is above 10 mm). One radius moves: a droplet of exactly 10 mm,
# which closed the last bin on the right and now falls in bin 40, not 39; every other radius bins as before. A body saved
# with 40 bins is padded with zeros on resume (`cli.cmd_resume`).
N_BINS, R_MIN, R_MAX = 50, 1.0e-6, 1.0e-1
BIN_EDGES = np.logspace(np.log10(1.0e-6), np.log10(1.0e-2), 41)
BIN_EDGES = np.concatenate([BIN_EDGES, np.logspace(-2.0, np.log10(R_MAX), 11)[1:]])
WE_BREAKUP = 12.0                  # Pilch-Erdman: secondary breakup expected above this droplet Weber number
CAPILLARY_TAU = 4.0 * np.pi / (2.0 * np.pi) ** 1.5      # 0.798: tau* = 2 x 2 pi / omega_cap(lambda*)
# h_floor_m and gamma_tau: viscous_film's floor depth and gamma tau* (nan off Girin & Kopyt's branches; 2026-10-22);
# solid_fraction and h_kg: the droplets' solid share of their mass and the enthalpy per kilogram they carry [J/kg, on
# the material's enthalpy scale] (spec 8 and 10.2, sub-plan 22 Task 2): the element model's film is liquid by
# construction and booked at its own liquid enthalpy, so its rows carry 0 and that enthalpy (`MeltingBody`);
# release_angle_register_deg: Register et al.'s (2020) release angle, 180 degrees at the stagnation point, 90 at the
# equator, 0 directly behind (`register_angle`; Asha 2026-10-29, sub-plan 22 Task 4)
SOURCE_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "theta_deg", "x_m", "y_m", "z_m", "closure", "branch", "b_m",
                  "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup",
                  "h_floor_m", "gamma_tau", "solid_fraction", "h_kg", "release_angle_register_deg"]
BRANCH_THICK, BRANCH_THIN, BRANCH_RAREFIED, BRANCH_RT = 0, 1, 2, 3      # Girin thick / thin-film with the edge state /
#                                                                        thin-film with the freestream / front-surface Rayleigh-Taylor
BRANCH_REGIME2 = 4                 # deep film whose melt is the more viscous medium: Girin & Kopyt's step profile again
BRANCH_PROFILE = 5                 # Girin's gradient mode solved on the column's own resolved profile (the profile path)
# The profile path's record columns, appended after SOURCE_COLUMNS (spec 6.8): the depth to the wall [m], the wall's kind
# (1 floor, 0 the rigid front), the floor Reynolds number at the wall, the first kink's depth [m] and temperature [K], the
# gradient mode's release depth [m], its solved growth rate [1/s], and Girin's theoretical delta_m [m] and V_s [m/s]. No
# name begins with `wall_`, the prefix of the wall-clock timers, which the comparisons of runs skip (the first two were
# renamed on 2026-11-12, in the final review's fix round, before any flight wrote them).
PROFILE_COLUMNS = ("depth_to_wall_m", "floor_like_wall", "floor_re", "kink_depth_m", "kink_T_K", "release_depth_m",
                   "growth_rate_per_s", "girin_delta_m_m", "girin_v_s_ms")


def record_columns(width):
    """The record columns' names for a record `width` columns wide: SOURCE_COLUMNS for 27 (the table path and the element
    model), SOURCE_COLUMNS followed by PROFILE_COLUMNS for 36 (the profile path); any other width is a ValueError."""
    if width == len(SOURCE_COLUMNS):
        return list(SOURCE_COLUMNS)
    if width == len(SOURCE_COLUMNS) + len(PROFILE_COLUMNS):
        return SOURCE_COLUMNS + list(PROFILE_COLUMNS)
    raise ValueError("records are {} or {} columns wide, not {}".format(
        len(SOURCE_COLUMNS), len(SOURCE_COLUMNS) + len(PROFILE_COLUMNS), width))


@dataclass
class ProfileModes:
    """What the profile path needs per column (spec 6.5 and 6.6), all (n,) arrays: `resolved` (bool), the columns whose
    profile the caller resolved; `growth` [1/s], the fastest wave's growth rate on the profile above the wall (0 where it
    is stable or that layer is one straight segment); `k` [1/m], its wavenumber; `wall_depth` [m], the depth to the wall;
    `nu_eff` [m2/s], the layer's Couette viscosity above the wall; `we_s`, the Weber number recorded for the column (spec
    6.7: density times the square of the top segment's velocity jump times its depth, over surface tension); and
    `cap_gradient` and `cap_driven` [m], the depth at the step's start of each mode's release run, for the crest cap."""
    resolved: np.ndarray
    growth: np.ndarray
    k: np.ndarray
    wall_depth: np.ndarray
    nu_eff: np.ndarray
    we_s: np.ndarray
    cap_gradient: np.ndarray
    cap_driven: np.ndarray


def melt_layer(flow, liquid):
    """delta_m per patch from Girin's Eq. (2) and the shear velocity factor (alpha mu)^(1/3)/(1 + (alpha mu)^(1/3)).

    Both are nan/0 wherever the patch's closure is not CLOSURE_GIRIN: Eq. (2) is a conjugate-boundary-layer result and
    Ranger's Psi(phi) for delta_a is a continuum construction, so the wall Knudsen gate of surface_flow decides whether
    they may be used at all (spec section 18). Where they may not, the film takes the Couette closure
    V_s = tau_w b / mu_melt (film.lubrication's thin branch) and the thin-film spray mode."""
    girin = (flow.closure == CLOSURE_GIRIN) & np.isfinite(flow.delta_a)
    alpha = flow.rho_e / liquid.rho
    mu = flow.mu_e / liquid.mu
    with np.errstate(divide="ignore", invalid="ignore"):
        delta_m = np.where(girin & (mu > 0.0), (alpha / np.where(mu > 0.0, mu, 1.0) ** 2) ** (1.0 / 3.0) * flow.delta_a, np.nan)
        am = (alpha * mu) ** (1.0 / 3.0)
        factor = np.where(girin & np.isfinite(am), am / (1.0 + am), 0.0)
    return delta_m, factor


def momentum_layer(tau, moment, radius, liquid):
    """delta_m [m] per patch of the slip band from a momentum balance (sub-plan 22 Task 5i; the slip-band plan of 2026-11-06,
    Asha's decision 2): with I = integral of r tau ds [Pa m2] the gas shear integrated over the meridian from the nose
    (`moment`) and r = `radius` the patch's distance from the axis, the plan's derivation gives delta_m^3 = 3 mu^2 I / (rho
    tau^2 r), the liquid's mu and rho and the patch's shear tau. Then the film's surface speed V_s = tau delta_m / mu makes
    the surface Weber number rho V_s^2 delta_m / sigma equal to 3 I / (r sigma), free of the viscosity. nan where tau or r
    is not above zero (the nose, a patch with no shear) or I is not finite and above zero (a leeward patch), and no warning."""
    tau, moment, radius = (np.asarray(x, dtype=float) for x in (tau, moment, radius))
    with np.errstate(divide="ignore", invalid="ignore"):
        ok = (tau > 0.0) & (radius > 0.0) & (moment > 0.0) & np.isfinite(moment)
        d = np.cbrt(3.0 * liquid.mu ** 2 * moment / (liquid.rho * tau ** 2 * radius))
    return np.where(ok, d, np.nan)


def spray_layer(flow, liquid, band=False, theta=None, radius=None):
    """The depth the spraying reads, delta_m per patch: `melt_layer`'s (Girin's Eq. (2), nan off his closure) except, with
    `band` on, on the patches in the slip band, where it is `momentum_layer`'s on the windward patches at r = `radius` sin
    (theta) from the nose radius the shear integral was taken at (`flow.nose_radius`, which the flow records; `radius` only for
    a flow that does not, a state from before Task 5j) and the angles `theta`. A patch is in the band where the body is on the merged
    branch with a body Knudsen number below 0.1 (`flow.in_slip_band`) and the patch's own wall Knudsen number is below 0.1
    too (Asha, 2026-11-10; the 2026-11-06 plan with the wall condition added). Only the spraying calls it: the film's
    transport and the feed keep `melt_layer`. With `band` off, outside the band, or for a flow with no integrated shear (a
    state from before the band), the result is `melt_layer`'s to the last bit."""
    delta_m = melt_layer(flow, liquid)[0]
    moment = getattr(flow, "shear_integral", None)
    recorded = getattr(flow, "nose_radius", None)                  # the radius the flow's shear integral was taken at
    radius = recorded if recorded is not None else radius
    if not band or moment is None or theta is None or radius is None or not getattr(flow, "in_slip_band", False):
        return delta_m
    theta = np.asarray(theta, dtype=float)
    with np.errstate(invalid="ignore"):
        inside = (np.asarray(flow.kn_local) < KN_LOCAL_SLIP) & (theta <= 0.5 * np.pi)
    return np.where(inside, momentum_layer(flow.tau, moment, radius * np.sin(theta), liquid), delta_m)


def rayleigh_taylor_bounded(deceleration, h, extent, liquid):
    """The front-surface Rayleigh-Taylor mode in a pool of *finite* lateral extent, which is what a melt pool is.

    `deceleration` is the surface-normal component W cos(phi), as in `rayleigh_taylor`; `extent` is the lateral extent
    of the contiguous molten region the patch belongs to, not of the patch (the mesh facet is not a physical boundary).

    For heavy liquid over light gas the dispersion relation with surface tension is n^2 = (W k - sigma k^3/rho) tanh(k h),
    unstable for k below the cutoff k_c = sqrt(W rho/sigma) -- that is, for wavelengths *longer* than
    lambda_c = 2 pi sqrt(sigma/(W rho)) -- and fastest at k_max = k_c/sqrt(3), which is the lambda* the unbounded form
    reports. A pool of lateral extent L with its rim pinning the interface admits only k >= k_1 = pi/L (half a
    wavelength across the pool), so:
      * L < lambda_c/2  : no admissible wavelength is unstable -- the pool is stable however deep it is;
      * L < lambda*/2   : unstable, but the fastest mode does not fit and the growth rate is the one at k_1, not 1/tau*;
      * otherwise       : the fastest mode fits and the unbounded growth rate applies.
    The tanh(k h) factor is the finite-depth correction and is not negligible here: at h = 2 mm and W = 30 m/s^2,
    k h = 0.58 and tanh = 0.52, so a shallow pool grows at about half the deep-layer rate.

    Comparing lambda* against the nose diameter, as an earlier diagnostic did, is wrong twice over: the mode must fit
    in the *unstable region*, not on the nose, and the yes/no question is decided by the marginal wavelength lambda_c,
    not by the fastest-growing one (measured 2026-09-23). Returns (active, wavelength [m], growth time [s])."""
    W = np.maximum(np.asarray(deceleration, dtype=float), 0.0)    # the surface-normal component, W cos(phi)
    h = np.asarray(h, dtype=float)
    L = np.asarray(extent, dtype=float)
    shape = np.broadcast_shapes(W.shape, h.shape, L.shape)
    live = W > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        k_c = np.where(live, np.sqrt(np.where(live, W, 0.0) * liquid.rho / liquid.sigma), 0.0)
        k_max = k_c / np.sqrt(3.0)
        k_1 = np.where(L > 0.0, np.pi / np.where(L > 0.0, L, 1.0), np.inf)
    k = np.maximum(k_1, k_max)                              # the fastest wavelength the pool actually admits
    active = live & (k < k_c) & (h > 0.0) & np.isfinite(k)
    with np.errstate(invalid="ignore", divide="ignore"):
        n2 = (W * k - liquid.sigma * k ** 3 / liquid.rho) * np.tanh(np.clip(k * h, 0.0, 30.0))
        tau = np.where(active & (n2 > 0.0), 1.0 / np.sqrt(np.where(n2 > 0.0, n2, 1.0)), np.inf)
        lam = np.where(active, 2.0 * np.pi / k, np.inf)
    on = np.broadcast_to(active & (n2 > 0.0), shape).copy()
    return on, np.broadcast_to(lam, shape).copy(), np.broadcast_to(tau, shape).copy()


def rayleigh_taylor(deceleration, b, liquid):
    """(active, lambda*, tau*) of the 1994 front-surface RT mode for the film thickness b.

    `deceleration` is the component of the body deceleration **normal to the film**, which is what drives the mode:
    Girin & Kopyt write it W sin(Theta) with Theta the angle between the relative velocity of the media and the
    deceleration vector, and on a sphere that is W cos(phi) -- full at the stagnation point, zero at the equator where
    the deceleration lies in the surface and cannot lift the interface at all. Passing the whole deceleration instead
    over-drives every patch off the nose by 1/cos(phi) (corrected 2026-09-24). A scalar still broadcasts."""
    W = np.maximum(np.asarray(deceleration, dtype=float), 0.0)
    b = np.asarray(b, dtype=float)
    shape = np.broadcast_shapes(W.shape, b.shape)
    active = np.broadcast_to(W * b * b * liquid.rho > 3.0 * liquid.sigma, shape).copy()
    live = W > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(live, 2.0 * np.pi * np.sqrt(3.0 * liquid.sigma / np.where(live, W * liquid.rho, 1.0)), np.inf)
        tau = np.where(live, (27.0 * liquid.sigma / (4.0 * np.where(live, W, 1.0) ** 3 * liquid.rho)) ** 0.25, np.inf)
    return active, np.broadcast_to(lam, shape).copy(), np.broadcast_to(tau, shape).copy()


def thin_film_mode(mach, momentum_flux, liquid):
    """(lambda*, tau*) of the Girin & Kopyt 1994 side mode for gas Mach number and rho V^2 [Pa]."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(momentum_flux > 0.0, 1.5 * mach * liquid.sigma / momentum_flux, np.inf)
    tau = CAPILLARY_TAU * lam ** 1.5 * np.sqrt(liquid.rho / liquid.sigma)
    return lam, tau


def viscous_film(lam, tau, h, liquid, nu=None):
    """The thin mode's wave on a viscous film of depth h over a rigid floor (2026-10-22): (gamma_tau, outside).

    `nu` [m2/s], one value or one per patch, replaces the liquid's mu_l / rho_l where it is given: the profile path's
    nu_eff, the Couette viscosity of the layer above the wall (spec 6.5). Without it the result is unchanged, bit for bit.

    gamma is the amplitude damping rate of a capillary wave of wavenumber k = 2 pi / lambda* on a layer of depth h,
        gamma = 2 nu k^2 + k sqrt(nu omega / 2) / sinh(2 k h),   omega^2 = (Sigma k^3 / rho_l) tanh(k h),
    Lamb's interior damping plus the bottom (Stokes) boundary layer's of Hunt (1952), as Mei gives them in The Applied
    Dynamics of Ocean Surface Waves; nu = mu_l / rho_l. The patch releases nothing where gamma_tau = gamma tau* >= 1:
    viscosity takes the wave down at least as fast as it grows. On deep liquid the floor term vanishes and this is
    Lamb's 2 nu k^2 tau* alone, which binds only below a wavelength of about 3.5 um for liquid aluminium; Girin & Kopyt's
    own dissipation limit, their Eq. (9), carries no liquid viscosity and is a third of lambda* everywhere, so it never
    binds.

    `outside` marks the patches where h < sqrt(2 nu / omega): the floor's boundary layer reaches through the whole
    liquid and Girin & Kopyt's inviscid mode does not apply. That is a validity limit, not a finding that no droplets
    leave -- a sheared viscous film has long-wave instabilities of its own -- so the caller records the patches and
    the film they hold rather than dropping them silently. The factor 2 under the root is a convention (the Stokes
    layer's e-folding depth is sqrt(2 nu / omega)).

    h is the depth to the rigid floor, which is where the liquid stops: the non-rigid depth the regime test uses, not
    the film thickness b, since the thin branch exists because the rigid core lies within delta_m of the surface.
    OPEN ASSUMPTION: mushy material more than half liquid counts as non-rigid there but is far more viscous than the
    liquid, so the effective floor may lie inside that layer, and the damping would then be larger than computed."""
    nu = liquid.mu / liquid.rho if nu is None else np.asarray(nu, dtype=float)
    lam, tau, h = (np.asarray(x, dtype=float) for x in (lam, tau, h))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        k = 2.0 * np.pi / lam
        kh = k * h
        omega = np.sqrt(liquid.sigma * k ** 3 / liquid.rho * np.tanh(kh))
        inv_sinh = np.where(kh > 0.0, 2.0 * np.exp(-2.0 * kh) / -np.expm1(-4.0 * kh), np.inf)   # 1 / sinh(2 k h)
        gamma = 2.0 * nu * k ** 2 + k * np.sqrt(nu * omega / 2.0) * inv_sinh
        outside = ~(h >= np.sqrt(2.0 * nu / omega))
    return gamma * tau, outside


def wave_fits(lam, extent):
    """True where one whole wavelength fits inside the contiguous molten region the patch belongs to.

    A mode cannot form in a puddle narrower than itself. The region, not the patch, is the domain: a mesh facet is a
    bookkeeping boundary, not a physical one, and the melt is continuous across it. `extent` of None disables the test;
    a region of zero extent admits nothing. The Rayleigh-Taylor form of the same argument is in
    `rayleigh_taylor_bounded`, where the rim pins the interface and half a wavelength is enough; a travelling shear
    wave has to complete a period inside the melt, so one whole wavelength is required here (decided 2026-09-24)."""
    if extent is None:
        return np.ones(np.shape(lam), dtype=bool)
    L = np.asarray(extent, dtype=float)
    with np.errstate(invalid="ignore"):
        return (L > 0.0) & np.isfinite(lam) & (lam <= L)


@dataclass
class SprayResult:
    branch: np.ndarray          # 0 thick, 1 thin, 2 rarefied, -1 none
    we_s: np.ndarray
    unstable: np.ndarray        # bool
    r: np.ndarray               # m, droplet radius (nan where none)
    mdot: np.ndarray            # kg/(m2 s) stripping rate
    dm: np.ndarray              # kg released this step per patch
    dn: np.ndarray              # droplets released this step per patch
    rt_active: bool
    delta_m: np.ndarray
    v_s: np.ndarray
    growth: np.ndarray = None   # s, the selected mode's growth time per patch (nan where no mode): what decides
    #                             whether the shear instability outruns the Rayleigh-Taylor mode on the same liquid
    h_floor: np.ndarray = None  # m, the depth to the rigid floor viscous_film used (nan off Girin & Kopyt's branches)
    gamma_tau: np.ndarray = None    # viscous_film's gamma tau* (nan off Girin & Kopyt's branches)
    outside: np.ndarray = None  # bool, Girin & Kopyt's branches whose film is shallower than the floor's Stokes layer
    wavelength: np.ndarray = None   # m, the wavelength of the mode that sprays each patch (nan where none): the
    #                                 per-crest cap's lambda (2026-10-24)
    cap_depth: np.ndarray = None    # m, the crest cap's depth per column on the profile path (spec 6.6); None off it


SPRAY_GATES_NAMES = ("2026-11-10", "before")


class SprayModel:
    """Girin's and Girin and Kopyt's release model for one body's patches (`evaluate`).

    `spray_gates` (sub-plan 22 Task 5h; Asha, 2026-11-10) says which of the release rules decided on 2026-11-10 apply:
    `before` (the class default, so that a model pickled before the setting existed behaves as it did) keeps every rule as
    it stood; `2026-11-10` applies them. So far the one rule: on the crest-limit path (`cap_depth` given, the layered
    skins) a Rayleigh-Taylor droplet is capped by its crest alone, not by a quarter of the body radius as well. Task 5i
    adds the thin branch's Knudsen gate under the same setting (below).

    The thin branch's continuum gate, under `2026-11-10` (Task 5i; Asha, 2026-11-10): the thin branch, and the rarefied
    branch with it, release only on patches with a conjugate depth to spray on, a finite `delta_m` -- Girin's closure below
    a wall Knudsen number of 0.01 (`melt_layer`), or, with `slip_band` on, the slip band's momentum-balance depth
    (`spray_layer`). Elsewhere a wet patch takes no shear mode, because the thin branch sprays only where Girin's continuum
    closure applies (or the band's depth exists), not because the edge state it reads is invalid (review of Task 5i, Minor 5);
    the Rayleigh-Taylor mode, a melt-side mode driven by the
    deceleration, is not gated. `before` keeps the thin branch spraying wherever it did. `slip_band` (class default False,
    so a model pickled before it existed has it off) says whether the callers read the band's depth: the model itself only
    gates on the `delta_m` it is given, which the callers take from `spray_layer` with this flag.

    `dispersion` (spec 6.1; one of `DISPERSION_NAMES`) says which theory the callers spray by: `table`, Girin's tabulated
    dispersion with his delta_m and V_s, or `profile`, his gradient mode solved on the column's own resolved profile. The
    class attribute and the constructor's default are both `table`, so a model pickled before the setting existed, or
    built directly, keeps Girin's theory; the command line passes `profile` for spraying skin runs. Only a setting other
    than the class's is stored on the instance, so a `table` model pickles to the same bytes as before the setting existed
    (a stored "table" would be the very string object that names the `table` attribute, which pickle's memo then shares,
    and a model loaded and pickled again would differ by those bytes). The model only carries the setting: `evaluate` does
    not read it, and takes the profile path wherever it is given a `profile`."""
    spray_gates = "before"
    slip_band = False
    dispersion = "table"

    def __init__(self, liquid, k_r=K_R, k_t=K_T, we_critical=WE_CRITICAL_THEORY, table=None,
                 rt_spray=True, spray_gates="before", slip_band=False, dispersion="table"):
        if spray_gates not in SPRAY_GATES_NAMES:
            raise ValueError("spray_gates must be one of {}, got {!r}".format(SPRAY_GATES_NAMES, spray_gates))
        if dispersion not in DISPERSION_NAMES:
            raise ValueError("dispersion must be one of {}, got {!r}".format(DISPERSION_NAMES, dispersion))
        self.spray_gates, self.slip_band = spray_gates, bool(slip_band)
        if dispersion != SprayModel.dispersion:      # the class default is not stored: a table model pickles as before
            self.dispersion = dispersion
        if k_r <= 0.0 or k_t <= 0.0 or we_critical <= 0.0:
            raise ValueError("k_r, k_t and the critical Weber number must be > 0")
        self.liquid, self.k_r, self.k_t, self.we_critical = liquid, k_r, k_t, we_critical
        self.rt_spray = bool(rt_spray)               # apply the front-surface mode, or only report it
        self.table = table or DispersionTable()

    def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, b_layer=None,
                 extent=None, deceleration_n=None, regime_layer=None, on_slurry=None, cap_depth=None, rt_layer=None,
                 profile=None, rt_extent=None):
        """Per-patch instability and release for the film thickness b, melt-layer thickness delta_m and film surface
        velocity v_s (from film.lubrication), over the step dt; m_f is the film mass available, radius the body's.

        `b_layer`, if given, is the depth of liquid beneath the wall (film plus the contiguous molten material under
        it): the thick/thin test belongs on that, since it asks whether the gas shear reaches the bottom of the liquid,
        while the release rates and the droplet cap stay on the film, which is the mass that can actually leave.
        `extent`, if given, is the lateral extent of the contiguous molten region each patch belongs to, for the
        wave-fits test; `deceleration_n` is the surface-normal deceleration W cos(phi) for the reported
        Rayleigh-Taylor flag, defaulting to the whole deceleration when it is not supplied.

        `regime_layer`, if given, replaces `b_layer` in the regime test only: the depth of non-rigid material under the
        wall, the film plus everything more than half liquid beneath it (`MeltingBody.nonrigid_depth`). The thin branch
        is Girin's dominant ablation, where the rigid core stabilises the disturbances, and slurry is no rigid core, so
        a film on slurry deeper than delta_m is deep melt and takes the thick branch. `on_slurry`, if given, flags the
        patches whose film rests on slurry; where there is no delta_m (off Girin's closure) such a film takes no shear
        mode at all, since the thin mode needs a rigid wall and no deep-melt mode exists there. The Rayleigh-Taylor
        criteria, a liquid pool's, keep `b_layer` (amendment of 2026-10-06).

        `cap_depth`, if given, is the depth [m] the release can take from each patch, and the droplet's radius is capped
        by one crest of its own wave on it, (3 lambda^2 h / (4 pi))^1/3, instead of by the film mass m_f (the module
        docstring; 2026-10-24); a patch with no such depth releases nothing. Without it the film-mass cap is kept, bit
        for bit.

        `rt_layer`, if given, is the depth [m] of the liquid pool the front-surface Rayleigh-Taylor mode reads -- its depth
        criterion, its bounded growth and the depth it sheds within a growth time -- in place of the liquid layer (for the
        criteria) and the film b (for the rate): the layered skins' liquid depth, where b is the flowing depth and liquid
        lies on slurry (Asha, 2026-10-29, decision 1 of sub-plan 22 Task 1c, option B: the mode is a liquid pool's, as the
        amendment of 2026-10-06 has it in the element model). Without it the mode reads them as before, bit for bit.

        `profile`, a `ProfileModes`, puts the columns it resolves on the profile path (the module docstring; spec 6.5 and
        6.6): they leave the thick, regime-2, thin and rarefied branches, whose thick/thin gate is not used there, and
        take the faster of the gradient mode solved on their profile and the gas-driven mode on their wall
        (`_profile_shear`), only where they are windward, delta_m is finite and the layer above the wall is at least B_MIN
        deep, whatever `spray_gates` says. Their crest cap reads `profile.cap_gradient` where the gradient mode took the
        column and `profile.cap_driven` where the gas-driven mode did, and the result records the depth used,
        `cap_depth`. Every other column takes today's code, and the caller's `cap_depth`, which the profile path needs.
        `rt_extent`, if given, is the region the Rayleigh-Taylor mode's wave must fit in, in place of `extent` (the
        skins' liquid pool, where `extent` is the region of the layer above the wall). With neither, nothing changes."""
        if profile is not None and cap_depth is None:
            raise ValueError("the profile path is the crest-limit path: a profile needs cap_depth")
        liq = self.liquid
        n = b.size
        layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
        regime = layer if regime_layer is None else np.maximum(np.asarray(regime_layer, dtype=float), b)
        branch = np.full(n, -1)
        growth = np.full(n, np.nan)          # the mode's growth time per patch [s], for comparison with other modes
        h_floor, gamma_tau = np.full(n, np.nan), np.full(n, np.nan)
        wavelength = np.full(n, np.nan)      # the spraying mode's wavelength per patch [m], for the per-crest cap
        outside = np.zeros(n, dtype=bool)
        r = np.full(n, np.nan)
        mdot = np.zeros(n)
        has_film = windward & (b >= B_MIN)
        # Girin's regime test, in his order, two stages (2026-09-25). Stage one: delta_m is a *theoretical* thickness --
        # the melt velocity boundary layer Eq. (2) predicts -- so the question the liquid depth answers is whether that
        # layer can physically form in the melt that is present. If the liquid cannot accommodate it, layer <= delta_m,
        # the layer did not form, the rigid core still stabilises the disturbances and the case is Girin's dominant
        # ablation: regime 1, his Girin & Kopyt (1994) mode, the thin branch below. If it did form, layer > delta_m, the
        # profile is the conjugated pair (delta_a in air, delta_m in melt) and the case is dominant fusion.
        deep = has_film & np.isfinite(delta_m) & (regime > delta_m)        # delta_m fits in the liquid: it formed
        # Stage two, only where it formed: which mechanism the conjugated profile carries, decided on the kinematic
        # viscosities. nu_melt > nu_gas puts the thicker viscous layer in the melt, V_s << V_a and the profile close to a
        # tangential discontinuity -- classical Kelvin-Helmholtz, regime 2. nu_gas > nu_melt makes the two thicknesses
        # comparable and the profile inflated, which is the gradient instability of his Eq. (1) -- regime 3, the thick
        # branch. Both stages are Girin (2017), Sects. 1 and 2; see the module docstring for the quotations and for the
        # measurement that puts liquid aluminium firmly in regime 3.
        nu_melt = liq.mu / liq.rho
        with np.errstate(divide="ignore", invalid="ignore"):
            nu_gas = np.where(flow.rho_e > 0.0, flow.mu_e / np.where(flow.rho_e > 0.0, flow.rho_e, 1.0), np.inf)
        step_profile = np.asarray(nu_melt > nu_gas)                        # the melt carries the thicker viscous layer
        regime2 = deep & step_profile                                     # Kelvin-Helmholtz
        thick = deep & ~step_profile                                      # regime 3: Girin 2017's gradient instability
        free_molecular = flow.branch == BRANCH_FREE_MOLECULAR
        # a film on slurry with no delta_m: the thin mode needs a rigid wall and there is no deep-melt mode to take it
        held = np.zeros(n, dtype=bool) if on_slurry is None else has_film & ~np.isfinite(delta_m) & np.asarray(on_slurry, dtype=bool)
        rarefied = has_film & ~deep & free_molecular & ~held              # no edge state exists: the freestream drives the mode
        thin = has_film & ~deep & ~rarefied & ~held
        if self.spray_gates != "before":
            # Asha, 2026-11-10 (sub-plan 22 Task 5i): the thin and rarefied branches shear-spray only where a conjugate depth
            # exists, Girin's closure or the slip band's momentum balance; elsewhere the wet patch takes no shear mode
            # (the Rayleigh-Taylor mode below is not gated)
            with_depth = np.isfinite(delta_m)
            rarefied, thin = rarefied & with_depth, thin & with_depth
        if profile is not None:
            # the resolved columns leave today's branches (`held` acts only through thin and rarefied); regime 2 is not
            # evaluated on them, since the thick/thin gate is not used there (spec 6.7, the plan's call 4)
            resolved = np.asarray(profile.resolved, dtype=bool)
            thick, regime2, thin, rarefied = thick & ~resolved, regime2 & ~resolved, thin & ~resolved, rarefied & ~resolved
        # Girin's (2017) stability criterion, on the thick branch only (2026-10-22). His phi_cr comes from We_s > We_cr
        # with the sphere's own edge solution substituted in; evaluated locally instead of through that closed form it holds on an eroded
        # body too, which a single critical angle does not. The shear acts over the liquid it can reach, min(delta_m,
        # layer), so the thick branch is his We_s = rho_l V_s^2 delta_m / Sigma unchanged while a film shallower than the
        # conjugate layer is measured on its own depth. Below We_cr the surface is stable and nothing leaves. The thin
        # branch had no stability threshold before this (Girin & Kopyt's inviscid side mode has none of its own, being
        # unstable at every wavelength), so it stripped film at any angle however small the shear (added 2026-09-24).
        # Removed from Girin & Kopyt's branches (thin, rarefied, regime 2) on 2026-10-22: their Eq. (5) has no critical
        # Weber number (the growth rate is positive at every We), and the threshold it does have once the wavelength is
        # bounded is a wave that fits inside the molten region -- wave_fits on lambda*, which also stops the stagnation
        # point, where lambda* grows without bound as the shear vanishes. They are gated instead by wave_fits and by
        # viscous_film. We_s is still computed on them and written to the source table.
        # A patch the regime test makes deep is measured on delta_m whatever its liquid depth -- a film on slurry deeper
        # than delta_m included (2026-10-06); everywhere else the expression is the one above, bit for bit.
        L = None if extent is None else np.asarray(extent, dtype=float)
        with np.errstate(invalid="ignore"):
            d_shear = np.where(deep, delta_m, np.where(np.isfinite(delta_m), np.minimum(delta_m, layer), layer))
            we_s = np.where(has_film, liq.rho * v_s ** 2 * d_shear / liq.sigma, 0.0)
        supercritical = we_s > self.we_critical
        # thick: Girin 2017
        if thick.any():
            vs, dm_ = v_s[thick], delta_m[thick]
            delta_f, im_f, _ = self.table(we_s[thick])
            with np.errstate(divide="ignore", invalid="ignore"):
                lam_f = np.where(delta_f > 0.0, 2.0 * np.pi * dm_ / np.where(delta_f > 0.0, delta_f, 1.0), np.inf)
            unstable = (supercritical[thick] & (im_f > 0.0) & (vs > 0.0)
                        & wave_fits(lam_f, None if L is None else L[thick]))
            lam = np.where(unstable, lam_f, np.nan)
            rr = self.k_r * lam
            with np.errstate(divide="ignore", invalid="ignore"):
                t_per = np.where(unstable, self.k_t * dm_ / (vs * np.where(unstable, im_f, 1.0)), np.inf)
                rate = np.where(unstable, liq.rho * np.pi * rr * rr / (lam * t_per), 0.0)
            r[thick], mdot[thick], branch[thick] = rr, rate, BRANCH_THICK
            growth[thick] = t_per
            wavelength[thick] = lam
        # regime 2: the other half of the deep-film set -- Girin & Kopyt's step profile, in their deep-film limit.
        # Wavelength and growth time are their Eqs. (11) and (12), the same pair the thin branch uses, because the mode is
        # the same one; what differs is the supply. The thin branch multiplies by min(b, lambda*/8) because a film
        # shallower than the wave can give up no more than its own depth, and here the film is deeper than the sheared
        # layer, so that factor is dropped and the rate is Girin's (2017) torus shedding on this wavelength,
        # mdot = rho_l pi r^2/(lambda tau*) with r = k_r lambda -- one torus of cross-section r per wavelength per growth
        # time. mdot stays the instability's *demand*; what the film can actually supply is imposed once, below, by
        # dm = min(mdot A dt, m_f), exactly as on every other branch.
        if regime2.any():
            m = regime2
            lam, tau = thin_film_mode(flow.mach_e[m], (flow.rho_e * flow.u_eff ** 2)[m], liq)
            gt, out = viscous_film(lam, tau, regime[m], liq)
            h_floor[m], gamma_tau[m], outside[m] = regime[m], gt, out
            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & (gt < 1.0) & ~out
                  & wave_fits(lam, None if L is None else L[m]))
            rr = np.where(ok, self.k_r * lam, np.nan)
            with np.errstate(divide="ignore", invalid="ignore"):
                rate = np.where(ok, liq.rho * np.pi * self.k_r ** 2 * lam / np.where(ok, tau, 1.0), 0.0)
            r[m], mdot[m], branch[m] = rr, rate, BRANCH_REGIME2
            growth[m] = np.where(ok, tau, np.nan)
            wavelength[m] = np.where(ok, lam, np.nan)
        # thin: Girin & Kopyt 1994 with the edge state
        for mask, code, mach, flux in ((thin, BRANCH_THIN, flow.mach_e, flow.rho_e * flow.u_eff ** 2),
                                       (rarefied, BRANCH_RAREFIED, np.full(n, state.ma), np.full(n, state.freestream.rho * state.V ** 2))):
            if not mask.any():
                continue
            lam, tau = thin_film_mode(mach[mask], flux[mask], liq)
            gt, out = viscous_film(lam, tau, regime[mask], liq)
            h_floor[mask], gamma_tau[mask], outside[mask] = regime[mask], gt, out
            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & (gt < 1.0) & ~out
                  & wave_fits(lam, None if L is None else L[mask]))
            rr = np.where(ok, lam / 4.0, np.nan)
            rate = np.where(ok, liq.rho * np.minimum(b[mask], lam / 8.0) / np.where(ok, tau, 1.0), 0.0)
            r[mask], mdot[mask], branch[mask] = rr, rate, code
            growth[mask] = np.where(ok, tau, np.nan)
            wavelength[mask] = np.where(ok, lam, np.nan)
        # the profile path's two shear modes (spec 6.5), on the resolved columns that may shear-spray (decision 6)
        if profile is not None:
            wall = np.asarray(profile.wall_depth, dtype=float)
            on = resolved & np.asarray(windward, dtype=bool) & np.isfinite(delta_m) & (wall >= B_MIN)
            if on.any():
                code, r[on], mdot[on], growth[on], wavelength[on], gamma_tau[on], outside[on] = self._profile_shear(
                    flow, state, profile, on, None if L is None else L[on])
                branch[on], h_floor[on] = code, wall[on]
                we_s[on] = np.where(code == BRANCH_PROFILE, np.asarray(profile.we_s, dtype=float)[on], we_s[on])
        # (iv) the front-surface Rayleigh-Taylor mode, now *applied* rather than only reported.
        #
        # Girin & Kopyt assign the front surface -- where the flow influence is negligible, V_0 -> 0 and sin Theta -> 1 --
        # to this aperiodic mode rather than to their side mode, and with Girin's critical angle now gating the shear
        # branches the stagnation cap has no other mechanism at all. Three modelling decisions, none of them theirs,
        # because their paper gives no mass-loss rate for the front surface (its Table 2 has no mdot column):
        #   * the gate is *both* criteria, because each supplies something the other does not. The depth criterion
        #     W cos(phi) h^2 rho_l > 3 Sigma -- equivalently h > lambda*/2 pi -- is the one with a threshold, and it is
        #     what confines the mode to the nose: the normal deceleration falls as cos(phi) so the depth needed rises as
        #     1/sqrt(cos phi) while the melt gets shallower. The bounded criterion has no threshold of its own (with a
        #     wide molten region k = k_max is always admissible and tanh(k h) > 0 for any depth, so it answers yes
        #     wherever there is liquid); what it supplies is the admissible wavelength and the growth rate, including
        #     the lateral fit and the finite-depth factor. Requiring only the bounded form made every wet patch
        #     Rayleigh-Taylor-unstable, which is how this was caught;
        #   * the rate carries their side-surface construction across. For the side they give r_d = lambda*/4 and
        #     tau_d = tau*, which this model writes mdot = rho_l min(b, r_d/2)/tau*; for the front they give
        #     r_d ~ lambda* and tau_d ~ tau*, so the same form is mdot = rho_l min(b, lambda*/2)/tau*. Since lambda*
        #     is centimetre-scale here and b millimetre-scale, it reduces to rho_l b / tau*: the film is shaken off
        #     within one growth time, which is what their text describes. No new closure constant is introduced.
        #   * the droplet radius is lambda* (their front-surface r_d), then capped as every branch is, by the film mass
        #     on the patch and a quarter of the body radius -- and here the film-mass cap is what actually binds.
        # The two modes are two modes of one interface, so the one that reaches tearing-off first breaks it up: the
        # shorter growth time takes the patch and the other does not act on the same film. That is also what keeps the
        # shear from being spent twice (the double-counting checklist of the plan's fact 29).
        w_n = flow.deceleration if deceleration_n is None else deceleration_n   # normal component, W cos(phi)
        pool = layer if rt_layer is None else np.asarray(rt_layer, dtype=float)  # the liquid the mode reads (2026-10-29)
        if self.rt_spray:
            deep, _, _ = rayleigh_taylor(w_n, pool, liq)                # the threshold: h > lambda*/2 pi
            L_rt = L if rt_extent is None else np.asarray(rt_extent, dtype=float)    # the liquid pool's region (spec 6.6)
            rt_on, lam_rt, tau_rt = rayleigh_taylor_bounded(w_n, pool, np.full(n, np.inf) if L_rt is None else L_rt, liq)
            rt_on = np.asarray(rt_on) & np.asarray(deep) & has_film & np.isfinite(tau_rt)
            shed = b if rt_layer is None else pool                       # what one growth time shakes off
            with np.errstate(divide="ignore", invalid="ignore"):
                rate_rt = np.where(rt_on, liq.rho * np.minimum(shed, lam_rt / 2.0) / np.where(rt_on, tau_rt, 1.0), 0.0)
            faster = rt_on & (rate_rt > 0.0) & (~np.isfinite(growth) | (tau_rt < growth))
            r = np.where(faster, lam_rt, r)
            mdot = np.where(faster, rate_rt, mdot)
            branch = np.where(faster, BRANCH_RT, branch)
            growth = np.where(faster, tau_rt, growth)
            wavelength = np.where(faster, lam_rt, wavelength)
        if cap_depth is not None:                # one crest gives nothing on no depth: nothing leaves (2026-10-24)
            h_cap = np.asarray(cap_depth, dtype=float)
            if profile is not None:              # the depth at the step's start of the run of the mode that took the column
                driven = resolved & np.isin(branch, (BRANCH_THIN, BRANCH_RAREFIED))
                h_cap = np.where(branch == BRANCH_PROFILE, np.asarray(profile.cap_gradient, dtype=float),
                                 np.where(driven, np.asarray(profile.cap_driven, dtype=float), h_cap))
            mdot = np.where(h_cap > 0.0, mdot, 0.0)
        unstable = mdot > 0.0
        dm = np.where(unstable, np.minimum(mdot * areas * dt, m_f), 0.0)
        # a droplet is never larger than a quarter of the body radius (the long near-critical waves Girin excludes; the
        # edge dynamic pressure vanishes toward the equator; not the Rayleigh-Taylor droplet on the crest-limit path under
        # spray_gates 2026-11-10, below) nor than the film on its patch -- or, given cap_depth, than
        # one crest of its own wave on the depth the release can take, (4/3) pi r^3 <= lambda^2 h (2026-10-24); under
        # spray_gates 2026-11-10 the Rayleigh-Taylor droplet on that path has the crest limit alone (2026-11-10)
        with np.errstate(invalid="ignore"):
            if cap_depth is None:
                r = np.minimum(r, np.minimum((3.0 * m_f / (4.0 * np.pi * liq.rho)) ** (1.0 / 3.0), 0.25 * radius))
            else:
                crest = (3.0 * wavelength ** 2 * h_cap / (4.0 * np.pi)) ** (1.0 / 3.0)
                if self.spray_gates == "before":
                    r = np.where(h_cap > 0.0, np.minimum(r, np.minimum(crest, 0.25 * radius)), np.nan)
                else:
                    # Asha, 2026-11-10 (sub-plan 22 Task 5h): the Rayleigh-Taylor droplet is its wavelength, Girin and
                    # Kopyt's front-surface r_d, and the quarter-radius cap is not from the mode's physics: it set the
                    # 15.8 mm droplets of the comparison flights through a transverse radius inflated to 57-64 mm by
                    # misplaced column tops. The crest limit, a supply limit, stays; every other branch keeps both.
                    limit = np.where(branch == BRANCH_RT, crest, np.minimum(crest, 0.25 * radius))
                    r = np.where(h_cap > 0.0, np.minimum(r, limit), np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            dn = np.where(dm > 0.0, dm / (4.0 / 3.0 * np.pi * liq.rho * np.where(dm > 0.0, r, 1.0) ** 3), 0.0)
        rt_active, _, _ = rayleigh_taylor(w_n, pool, liq)      # the unbounded criterion, still reported for continuity
        gk = np.isin(branch, (BRANCH_THIN, BRANCH_RAREFIED, BRANCH_REGIME2))    # the Rayleigh-Taylor mode may have taken some
        h_floor, gamma_tau = np.where(gk, h_floor, np.nan), np.where(gk, gamma_tau, np.nan)
        return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s, growth,
                           h_floor, gamma_tau, outside & gk, wavelength, None if profile is None else h_cap)

    def _profile_shear(self, flow, state, profile, on, L):
        """The profile path's two shear modes on the columns `on` (spec 6.5 and 6.6), whose region extents are `L` (None
        for no wave-fit test): (branch, r [m], mdot [kg/(m2 s)], growth time [s], wavelength [m], gamma_tau, outside), each
        one per column of `on`.

        The gradient mode solved on the profile above the wall: lambda = 2 pi / k, r = k_r lambda, a growth time k_t /
        growth and Girin's torus rate rho pi r^2 / (lambda t), active where its growth beats its own wave's viscous
        damping, growth > 2 nu_eff k^2, and the wave fits in its region, with no critical Weber number. The gas-driven mode:
        Girin & Kopyt's thin branch with the edge state (the rarefied one, with the freestream, on the free-molecular
        branch) on a film as deep as the wall, rho min(wall, lambda*/8) / tau*, with `viscous_film` on that depth and
        nu_eff. The shorter growth time takes the column (the comparison the Rayleigh-Taylor mode already makes, the plan's
        call 3); a tie, and a column neither mode sprays, go to the gas-driven mode, which keeps such a column with no rate
        as today's thin branch keeps its patches, so that its viscous-film record and `outside` flag are kept."""
        liq, n = self.liquid, on.size
        k, rate, nu, wall = (np.asarray(x, dtype=float)[on] for x in (profile.k, profile.growth, profile.nu_eff,
                                                                       profile.wall_depth))
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            lam_g = 2.0 * np.pi / k
            grad = np.isfinite(lam_g) & (lam_g > 0.0) & (rate > 2.0 * nu * k ** 2) & wave_fits(lam_g, L)
            t_g = np.where(grad, self.k_t / np.where(grad, rate, 1.0), np.inf)
            r_g = self.k_r * lam_g
            mdot_g = np.where(grad, liq.rho * np.pi * r_g * r_g / (lam_g * t_g), 0.0)
        free = np.broadcast_to(flow.branch == BRANCH_FREE_MOLECULAR, (n,))[on]
        mach = np.where(free, state.ma, flow.mach_e[on])
        flux = np.where(free, state.freestream.rho * state.V ** 2, (flow.rho_e * flow.u_eff ** 2)[on])
        lam, tau = thin_film_mode(mach, flux, liq)
        gt, out = viscous_film(lam, tau, wall, liq, nu=nu)
        ok = np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & (gt < 1.0) & ~out & wave_fits(lam, L)
        with np.errstate(divide="ignore", invalid="ignore"):
            mdot_d = np.where(ok, liq.rho * np.minimum(wall, lam / 8.0) / np.where(ok, tau, 1.0), 0.0)
        wins = grad & (t_g < np.where(ok, tau, np.inf))
        code = np.where(wins, BRANCH_PROFILE, np.where(free, BRANCH_RAREFIED, BRANCH_THIN))
        return (code, np.where(wins, r_g, np.where(ok, lam / 4.0, np.nan)), np.where(wins, mdot_g, mdot_d),
                np.where(wins, t_g, np.where(ok, tau, np.nan)), np.where(wins, lam_g, np.where(ok, lam, np.nan)), gt, out)


def source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None,
                 release_angle=None, extra=None):
    """Source-table rows, an array (k, len(SOURCE_COLUMNS)), for the patches that released mass (`res.dm` > 0): `h` [m]
    and `V` [m/s] the altitude and speed, `theta` [rad], `centroids` (n, 3) [m] and `t_hat` (n, 3) per patch, `flow`'s
    `closure`, `res` a SprayResult (or anything with its per-patch `branch`, `delta_m`, `we_s`, `r`, `dn`, `dm`, `v_s`,
    `h_floor`, `gamma_tau`), `b` [m] the film thickness, `liquid` the droplets' properties -- each of its `rho`, `mu` and
    `sigma` one value or one per patch (a slurry droplet's are its slurry's: sub-plan 22 Task 2) -- for their Weber and
    Ohnesorge numbers at the freestream `state`, `solid_fraction` and `h_kg` [J/kg] per patch (0 and nan if not given),
    and `release_angle` [deg] per patch, Register et al.'s (`register_angle`; nan if not given; sub-plan 22 Task 4).
    `extra`, an (n, len(PROFILE_COLUMNS)) array, is appended to each released patch's row: the profile path's columns
    (spec 6.8), so that its rows are `record_columns(36)` wide; without it the rows are SOURCE_COLUMNS, as before."""
    k = np.flatnonzero(res.dm > 0.0)
    if k.size == 0:
        return np.zeros((0, len(SOURCE_COLUMNS) + (0 if extra is None else np.shape(extra)[1])))
    r = res.r[k]
    n = len(res.dm)
    rho, mu, sigma = (np.broadcast_to(np.asarray(x, dtype=float), (n,))[k] for x in (liquid.rho, liquid.mu, liquid.sigma))
    we_d = state.freestream.rho * state.V ** 2 * 2.0 * r / sigma
    oh = mu / np.sqrt(rho * sigma * 2.0 * r)
    s_f = np.zeros(k.size) if solid_fraction is None else np.asarray(solid_fraction, dtype=float)[k]
    h_d = np.full(k.size, np.nan) if h_kg is None else np.asarray(h_kg, dtype=float)[k]
    angle = np.full(k.size, np.nan) if release_angle is None else np.asarray(release_angle, dtype=float)[k]
    return np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
                            centroids[k, 0], centroids[k, 1], centroids[k, 2], flow.closure[k], res.branch[k], b[k],
                            res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1],
                            t_hat[k, 2], we_d, oh, (we_d > WE_BREAKUP).astype(float), res.h_floor[k], res.gamma_tau[k],
                            s_f, h_d, angle] + ([] if extra is None else [np.asarray(extra, dtype=float)[k]]))


def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction=None, h_kg=None,
                release_angle=None):
    """Source-table rows (list of lists in SOURCE_COLUMNS order) for the patches that released mass this step
    (`source_block`'s, as lists: the element model's form, one row per patch-step)."""
    return source_block(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state, solid_fraction, h_kg,
                        release_angle).tolist()


def source_table(rows):
    """A body's `source_rows` as one array (N, width): a list of rows (the element model's, one per patch-step) or of
    blocks of rows (the layered skins', one array per macro step: a Python list of floats costs about 0.9 kB a row against
    0.2 kB in an array, and a flight at the continuum step makes of order ten million records), or a `SourceSpool` on disk,
    read into memory. The width is the rows' own, 27 or, on the profile path, 36 (`record_columns`), or the spool's names';
    an empty list has 27 columns."""
    if isinstance(rows, SourceSpool):
        return np.column_stack([np.asarray(c) for c in rows.columns().values()]) if rows.n_rows \
            else np.zeros((0, len(rows.names)))
    if not len(rows):
        return np.zeros((0, len(SOURCE_COLUMNS)))
    if isinstance(rows[0], np.ndarray) and rows[0].ndim == 2:
        return np.concatenate(rows, axis=0)
    return np.array(rows, dtype=float)


def source_columns(rows, names=None):
    """A body's `source_rows` as {name: one-dimensional array} in `record_columns` order for their width, what the coupled
    writer writes: a list of rows or of blocks, the columns of `source_table`'s array; a `SourceSpool`, its files mapped
    from disk under its own names, so that no flight's records need be in memory at once. The values are the same to the
    last bit either way. `names`, if given, are the columns of a list with no rows yet (the body's own record names, a
    profile-path skin body's 36, spec 6.8), which would otherwise read as the 27 SOURCE_COLUMNS."""
    if isinstance(rows, SourceSpool):
        return rows.columns()
    table = source_table(rows)
    if names is not None and not len(table):
        table = np.zeros((0, len(names)))
    return {name: table[:, i] for i, name in enumerate(record_columns(table.shape[1]))}


class SourceSpool:
    """The layered skins' droplet records on disk as they come (sub-plan 22 Task 4): one file per source column,
    `<directory>/<name>.f64` (float64, little-endian), each macro step's block appended to every one. Measured for the 100
    mm flight at the continuum step from the first liquid (24.5 s) to 120 s: 7 640 macro steps of 4 600-9 200 rows, 35-70
    million rows of 27 columns, 7.6-15.2 GB as arrays and twice that while `source_table` joins them -- more than the
    laptop's 24 GB -- so the command line gives the skin body a spool in the run's directory and the coupled writer maps
    the columns back (`columns`), each read once, contiguous. `len` counts the blocks (macro steps that released), as the
    in-memory list does; `n_rows` the rows. `names` are its columns, SOURCE_COLUMNS unless it is given others (the profile
    path's `record_columns(36)`, spec 6.8); a spool pickled before it had names has the class's, the 27 SOURCE_COLUMNS."""
    names = tuple(SOURCE_COLUMNS)

    def __init__(self, directory, names=None):
        os.makedirs(directory, exist_ok=True)
        self.directory, self.n_rows, self.n_blocks = directory, 0, 0
        self.names = tuple(SOURCE_COLUMNS if names is None else names)
        for name in self.names:
            open(self._path(name), "wb").close()

    def _path(self, name):
        return os.path.join(self.directory, name + ".f64")

    def __len__(self):
        return self.n_blocks

    def append(self, block):
        """Add a macro step's block of rows (k, len(names)) to the end of every column's file."""
        block = np.asarray(block, dtype=np.float64)
        if block.ndim != 2 or block.shape[1] != len(self.names):
            raise ValueError("a block of source rows has {} columns, not {}".format(block.shape[-1], len(self.names)))
        if not len(block):
            return
        for i, name in enumerate(self.names):
            with open(self._path(name), "ab") as fh:
                fh.write(np.ascontiguousarray(block[:, i], dtype="<f8").tobytes())
        self.n_rows += len(block)
        self.n_blocks += 1

    def columns(self):
        """{name: the column, read-only, mapped from its file} in `names` order (empty arrays with no rows)."""
        if not self.n_rows:
            return {name: np.zeros(0) for name in self.names}
        return {name: np.memmap(self._path(name), dtype="<f8", mode="r", shape=(self.n_rows,)) for name in self.names}

    def remove(self):
        """Delete the files and the directory, once what they hold has been written elsewhere."""
        for name in self.names:
            if os.path.exists(self._path(name)):
                os.remove(self._path(name))
        if os.path.isdir(self.directory) and not os.listdir(self.directory):
            os.rmdir(self.directory)


def fit_width(rows, width):
    """In-memory records fitted to `width` columns, for a resume that changes `SprayModel.dispersion` (spec 6.8; Review
    Focus 4): each block of a list cut to its first `width` columns (profile to table: the 27 SOURCE_COLUMNS) or padded on
    the right with nan columns (table to profile: the nine PROFILE_COLUMNS unknown), the columns kept bit for bit, as a new
    list. Any other form -- a `SourceSpool`, which a resume replaces, or an empty list -- comes back unchanged."""
    if not isinstance(rows, list) or not rows:
        return rows
    fitted = []
    for block in rows:
        block = np.asarray(block, dtype=float)
        have = block.shape[-1]
        if have > width:
            block = block[..., :width].copy()
        elif have < width:
            block = np.concatenate([block, np.full(block.shape[:-1] + (width - have,), np.nan)], axis=-1)
        fitted.append(block)
    return fitted


def register_angle(points, centre, v_hat):
    """Register et al.'s (2020, section 4.2) release angle [deg] of each point (n, 3) [m] on a body whose centroid is
    `centre` (3,) [m] and which moves through the air along `v_hat`: 180 degrees minus the angle between v_hat and the
    line from the centroid to the point, so 180 at the stagnation point, 90 at the equator and 0 directly behind,
    measured from the actual vector so that a flattened or receded body gets its own angle (Asha, 2026-10-29)."""
    d = np.asarray(points, dtype=float) - np.asarray(centre, dtype=float)
    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
    along = d @ v
    across = np.linalg.norm(d - along[:, None] * v, axis=1)
    return 180.0 - np.degrees(np.arctan2(across, along))


def enclosed_centroid(points, faces, normals=None):
    """The centroid [m, (3,)] of the volume a closed triangulated surface encloses (the divergence theorem, by signed
    tetrahedra from the points' mean, so either consistent orientation gives the same point): `faces` (nt, 3) into
    `points`; with `normals` (nt, 3), one outward normal per triangle, each triangle is first turned to agree with its
    own, for a surface whose triangles are not consistently ordered (the element model's). Register's angle is
    measured from this point (`register_angle`): the body's geometric centroid, which for a body of one density is its
    centre of mass."""
    points, faces = np.asarray(points, dtype=float), np.asarray(faces, dtype=np.int64)
    o = points[np.unique(faces)].mean(axis=0)
    a, b, c = (points[faces[:, i]] - o for i in range(3))
    n2 = np.cross(b - a, c - a)
    v = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0
    if normals is not None:
        v = np.where(np.einsum("ij,ij->i", n2, np.asarray(normals, dtype=float)) < 0.0, -v, v)
    return o + (v[:, None] * (a + b + c)).sum(axis=0) / (4.0 * v.sum())


# the history's release angle over a step's released mass, both surface models' (sub-plan 22 Task 4): its mean by mass,
# least and largest [deg], nan where the step released nothing
RELEASE_ANGLE_COLUMNS = ("release_angle_register_mean_deg", "release_angle_register_min_deg",
                         "release_angle_register_max_deg")


def release_angle_stats(dm, angle):
    """(mass-weighted mean, minimum, maximum) [deg] of the release angles `angle` of a step's records weighted by their
    masses `dm` [kg], over the records that released (dm > 0); nan, nan, nan where none did (the history's
    release_angle_register_*_deg, sub-plan 22 Task 4)."""
    dm, angle = np.asarray(dm, dtype=float), np.asarray(angle, dtype=float)
    on = dm > 0.0
    if not on.any():
        return float("nan"), float("nan"), float("nan")
    w, a = dm[on], angle[on]
    return float((w * a).sum() / w.sum()), float(a.min()), float(a.max())


def histogram(r, dn, dm):
    """(dn per bin, dM per bin) on BIN_EDGES for droplet radii r with counts dn and masses dm."""
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(r) & (r > 0.0)
    n_hist, _ = np.histogram(r[ok], BIN_EDGES, weights=np.asarray(dn)[ok])
    m_hist, _ = np.histogram(r[ok], BIN_EDGES, weights=np.asarray(dm)[ok])
    return n_hist, m_hist
