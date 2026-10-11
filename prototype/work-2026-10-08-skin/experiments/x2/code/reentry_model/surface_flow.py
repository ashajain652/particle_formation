"""Gas-side surface flow per patch: the flow-regime gate, the wall pressure, the boundary-layer edge state, the wall
Knudsen number, the melt closure it selects, the wall shear and the film's driving gradient (spec Step 3 sections 7
and 18, amended 2026-09-22).

REGIME GATE (body scale, three branches). Stage one asks whether a *distinct bow shock* exists, not whether the flow
is "continuum": a continuum construction must not be allowed to certify itself. The gate is on the body Knudsen
number Kn_body = lambda_inf / D (SESAM's own definition, hard-sphere with d = 3.65 A, verified to reproduce its
`knudsen` column to 1e-4) and the freestream Mach number:

  Kn_body < KN_BODY_SHOCK (0.01) and Ma_inf > 1   BRANCH_SHOCK_LAYER: normal shock -> wall pressure -> isentropic
      expansion to the edge state -> boundary layer -> wall Knudsen number -> the melt closure (below).
  KN_BODY_SHOCK <= Kn_body < KN_BODY_FM (10)      BRANCH_MERGED: the shock is not distinct from the shock layer, so no
      post-shock construction is valid. The wall loads are bridged between the free-molecular and continuum limits
      with SESAM's measured f(Kn); the melt closure is Couette; every patch is flagged.
  Kn_body >= KN_BODY_FM, or Ma_inf <= 1           BRANCH_FREE_MOLECULAR: Schaaf-Chambre loads from the freestream
      directly, with no compression model; Couette closure; flagged.

The 0.01 threshold is the standoff criterion: a bow shock is distinct only if it stands off farther than it is thick.
Standoff Delta/R is 0.14 (perfect gas) to 0.08 (real gas) and a strong shock is 3-10 upstream mean free paths thick,
so Delta > 5 lambda gives Kn_D <~ 0.01. The shock-layer Reynolds number Re2 = rho_inf V_inf R / mu(T0) (merged below
~100) agrees: 100 mm at 70 km gives Re2 = 175 and Kn_body = 0.0098 -- exactly on both thresholds -- while 5 mm at
77.5 km gives Re2 = 3.0, Kn_body = 0.60.

The gate is a DECLARED CONSERVATISM, not a physical deduction: Kn_body > 0.01 does not imply Kn_local > 0.01. The
wall gas is compressed and cold, so lambda_w / lambda_inf ~ 1/170 at the nose of the 100 mm sphere at 70 km (the
edge value is ~1/9) and even the 5 mm sphere at 77.5 km has a nose Kn_local of 0.0066 while its Kn_body is 0.60.
We decline to certify a patch as continuum because the construction we would use to check it is itself unreliable
there; `branch` and the closure fractions are recorded every step so the cost of that conservatism is measurable.
Nothing here claims that a merged or free-molecular patch sees only freestream density -- there is real compression
in the merged regime, and a cold diffusely reflecting wall raises the number density even in free-molecular flow;
the free-molecular branch simply evaluates the surface loads from freestream conditions with no compression model.

WALL PRESSURE (modified Newtonian + Prandtl-Meyer). Modified Newtonian alone gives p_w = p_inf at 90 degrees, a
factor ~25 below the measured sphere data (C_p 0.05-0.1), in exactly the band where Girin expects most of the
spraying; that error would propagate into rho_w, u_e, tau_w and We_s. So the windward face uses
  phi <= phi*:  p_w = p_inf + q_inf C_p,max cos^2 phi,    C_p,max = (p02 - p_inf)/q_inf
  phi >  phi*:  Prandtl-Meyer expansion through the turn angle phi - phi* (the surface tangent of a sphere rotates
                one-for-one with phi), p_w/p02 = [1 + (g-1) M^2/2]^(-g/(g-1)) with nu(M) = phi - phi*,
with phi* the sonic point, cos^2 phi* = (p* - p_inf)/(q_inf C_p,max), p*/p02 = (2/(g+1))^(g/(g-1)) -- 43.38 deg at
g = 1.4, 40.72 deg at g = 1.15. The two branches meet at phi* by construction but their slopes do not, so they are
blended over PM_BLEND_DEG degrees above phi* (the kink would otherwise show up in tau_w through the pressure gradient). p_w is
floored at p_inf. Prandtl-Meyer is used ONLY to produce p_w: the state itself comes from Cantera, (rho_e, T_e, h_e)
= isentropic expansion of the post-shock reservoir to p_w and u_e = sqrt(2 (h02 - h_e)), so the perfect-gas M_e never
propagates into the state variables. `gamma_pm` is a configuration parameter because the expansion's effective gamma
varies (1.4 frozen to ~1.15 dissociated): at 90 degrees for the 100 mm sphere at 70.0 km (V 7190 m/s, p02 4165 Pa)
the branch gives 144.3 Pa (g = 1.4) or 248.8 Pa (g = 1.15) against 112-219 Pa from measured C_p and 5.22 Pa from
Newtonian -- a factor ~2 uncertainty on a 28-48x improvement. The branch is not extended past PM_MAX_DEG (real flow separates); our windward
face ends at 90 degrees anyway.

No oblique-shock or entropy-swallowing machinery is built (spec section 18): equating the boundary-layer mass flow at
the shoulder with the mass crossing the shock inside radius y_s gives y_s/R ~ 0.04-0.1, where the local shock
inclination is 85-88 degrees and M_n = 0.996 M -- indistinguishable from normal. A sphere is all nose: the whole
windward hemisphere lies within ~1.6 nose radii of arc and the entropy layer is swallowed only many nose radii
downstream, which is also what justifies the normal-shock reservoir for the Prandtl-Meyer branch.

WALL KNUDSEN NUMBER AND THE MELT CLOSURE. Kn_local = lambda_w / R with lambda_w the Maxwell mean free path of the
gas *touching the surface* (p_w, T_wall) and R the nose radius -- a geometric length, not the boundary-layer
thickness (which only exists in the continuum, so using it would be circular) and not the running length s = R phi
(which vanishes at the stagnation point). Because the wall gas is ideal at these temperatures this is exactly
lambda_w = mu(T_w) sqrt(pi R_s T_w / 2) / p_w, i.e. Kn_local is inversely proportional to the wall pressure. The
identity Kn = (Ma/Re) sqrt(gamma pi/2) holds to machine precision (coefficient 1.4829 at gamma = 1.4; the
hard-sphere Chapman-Enskog form gives 1.5105 -- a 2 % difference).

Kn_local selects the closure that supplies the melt's velocity scale, NOT whether spraying happens: Girin's
dispersion relation contains no gas parameters -- it is a melt-side instability. What needs the continuum is his
Eq. (2) conjugate boundary layers and Ranger's Psi(phi) for delta_a:
  Kn_local < KN_LOCAL_CONTINUUM (0.01)   CLOSURE_GIRIN:   delta_m and V_s from Eq. (2)
  0.01 - KN_LOCAL_SLIP (0.1)             CLOSURE_COUETTE: V_s = tau_w b / mu_melt, slip-corrected tau_w, flagged
  >= 0.1                                 CLOSURE_COUETTE: free-molecular tau_w, flagged
tau_w is well defined in every branch (Schaaf-Chambre, tau_w = sigma_t rho_inf V^2 sin phi cos phi for diffuse
reflection), so the melt is sheared -- and can spray -- in all of them.
"""
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from . import aero, gas
from .constants import GAMMA_AIR, K_BOLTZMANN, M_BAR_AIR

BRANCH_SHOCK_LAYER, BRANCH_MERGED, BRANCH_FREE_MOLECULAR = 0, 1, 2
BRANCH_NAMES = ("shock layer", "merged", "free molecular")
CLOSURE_GIRIN, CLOSURE_COUETTE = 0, 1
KN_BODY_SHOCK, KN_BODY_FM = 0.01, 10.0            # body-scale gate (module docstring)
KN_LOCAL_CONTINUUM, KN_LOCAL_SLIP = 0.01, 0.1     # wall-scale closure gate
RE_SHOCK_MERGED = 100.0                            # Re2 below which the shock layer is merged (diagnostic only)
R_SPECIFIC_AIR = K_BOLTZMANN / M_BAR_AIR          # J/(kg K), 287.06
RANGER_C = 58.08                                   # 4.84 x 12: Ranger's delta_a from the u_e^4 integral (below)
THWAITES_C = 0.45
SIGMA_V, SIGMA_T = 1.0, 1.0                        # tangential momentum accommodation (Maxwell slip, Schaaf-Chambre)
GAMMA_PM = 1.15                                    # effective gamma of the Prandtl-Meyer expansion (config, see docstring)
PM_BLEND_DEG, PM_MAX_DEG = 10.0, 110.0   # blend window above phi* (10 deg: max |d2p| 9x the median, against 35x at 5 deg)
RAREFIED_SHEAR_NAMES = ("slip", "bridged")
THETA_BINS = np.radians(np.arange(0.0, 90.0 + 0.5, 1.0))


def ranger_psi(theta):
    theta = np.asarray(theta, dtype=float)
    return np.sqrt((6.0 * theta - 4.0 * np.sin(2.0 * theta) + 0.5 * np.sin(4.0 * theta)) / np.sin(theta) ** 5)


def ranger_thickness(radius, reynolds_diameter, theta):
    """delta_a = 2.2 R Re_D^-1/2 Psi(theta) (Ranger 1972, as used by Girin 2017)."""
    return 2.2 * radius * ranger_psi(theta) / np.sqrt(reynolds_diameter)


def boundary_layer_thickness(s, u_e, nu_e, constant=RANGER_C, power=4, refine=20):
    """delta_a on the meridian grid s from delta_a^2 = constant nu_e int u_e^power ds / u_e^(power+1) (Ranger's form
    for power = 4; Thwaites' momentum thickness is THWAITES_C with power = 5). The integrand is evaluated on a grid
    refined `refine` times (u_e interpolated linearly) so that the trapezoid rule resolves u_e^power ~ s^power near
    the stagnation point (1-degree bins alone give +17 % at 2 degrees, measured)."""
    s, u_e = np.asarray(s, dtype=float), np.asarray(u_e, dtype=float)
    fine = np.linspace(s[0], s[-1], refine * (len(s) - 1) + 1)
    integral = np.interp(s, fine, cumulative_trapezoid(np.interp(fine, s, u_e) ** power, fine, initial=0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        d2 = constant * nu_e * integral / u_e ** (power + 1)
    d2 = np.where(u_e > 0.0, d2, 0.0)
    if len(s) > 1 and u_e[0] == 0.0 and u_e[1] > 0.0:           # stagnation point: the limit is finite, take the neighbour's value
        d2[0] = d2[1]
    return np.sqrt(np.maximum(d2, 0.0))


def prandtl_meyer(mach, gamma=GAMMA_PM):
    """Prandtl-Meyer function nu(M) [rad] for M >= 1."""
    m = np.asarray(mach, dtype=float)
    b = math.sqrt((gamma + 1.0) / (gamma - 1.0))
    x = np.maximum(m * m - 1.0, 0.0)
    return b * np.arctan(np.sqrt(x) / b) - np.arctan(np.sqrt(x))


def mach_from_turn(turn, gamma=GAMMA_PM):
    """Inverse of prandtl_meyer: the Mach number reached after turning through `turn` [rad] from sonic conditions."""
    turn = np.asarray(turn, dtype=float)
    out = np.ones_like(turn)
    nu_max = prandtl_meyer(60.0, gamma)
    for k, t in enumerate(np.ravel(turn)):
        if t <= 0.0:
            out.ravel()[k] = 1.0
        else:
            out.ravel()[k] = brentq(lambda m: prandtl_meyer(m, gamma) - min(t, 0.999 * nu_max), 1.0 + 1e-9, 60.0, xtol=1e-10)
    return out


def sonic_angle(p_stag, p_inf, q_inf, gamma=GAMMA_PM):
    """The polar angle where modified-Newtonian pressure reaches the sonic value [rad] (43.38 deg at gamma 1.4)."""
    cp_max = (p_stag - p_inf) / q_inf
    p_star = p_stag * (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    c2 = (p_star - p_inf) / (q_inf * cp_max) if cp_max > 0.0 else 0.0
    return math.acos(math.sqrt(min(1.0, max(0.0, c2))))


def wall_pressure(theta, p_stag, p_inf, q_inf, gamma=GAMMA_PM, blend_deg=PM_BLEND_DEG):
    """p_w(theta) on the windward face: modified Newtonian to the sonic point, Prandtl-Meyer beyond, blended over
    `blend_deg` and floored at p_inf. Returns (p_w, phi_sonic [rad])."""
    theta = np.asarray(theta, dtype=float)
    phi_star = sonic_angle(p_stag, p_inf, q_inf, gamma)
    cp_max = (p_stag - p_inf) / q_inf
    newton = p_inf + q_inf * cp_max * np.cos(np.minimum(theta, 0.5 * math.pi)) ** 2
    turn = np.clip(theta - phi_star, 0.0, math.radians(PM_MAX_DEG) - phi_star)
    mach = mach_from_turn(turn, gamma)
    p_star = p_stag * (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    pm = p_stag * (1.0 + 0.5 * (gamma - 1.0) * mach ** 2) ** (-gamma / (gamma - 1.0))
    pm = np.where(turn > 0.0, pm, p_star)
    # the blend runs from phi* upward, never below it: the two branches already meet in value at phi* (both give p*),
    # but the Prandtl-Meyer slope is singular there (nu ~ (M-1)^3/2, so dp/dnu ~ nu^-1/3), so the Newtonian slope is
    # carried through the window by a smoothstep whose own derivative vanishes at both ends. A window straddling phi*
    # instead mixes in the clamped PM value below it and puts back the kink it was meant to remove (measured 2026-09-22).
    w = np.clip((theta - phi_star) / math.radians(blend_deg), 0.0, 1.0)
    w = w * w * (3.0 - 2.0 * w)
    return np.maximum(np.where(theta <= 0.5 * math.pi, (1.0 - w) * newton + w * pm, p_inf), p_inf), phi_star


def mean_free_path_maxwell(mu, rho, T, R_s=R_SPECIFIC_AIR):
    """Maxwell mean free path mu/rho sqrt(pi/(2 R T)); Kn = (Ma/Re) sqrt(gamma pi/2) is the same quantity."""
    return mu / rho * np.sqrt(math.pi / (2.0 * R_s * T))


def wall_knudsen(p_w, T_wall, mu_wall, length, R_s=R_SPECIFIC_AIR):
    """Kn_local = lambda_w / length with the wall gas ideal at (p_w, T_wall): lambda_w = mu sqrt(pi R T/2)/p_w."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = mu_wall * np.sqrt(math.pi * R_s * np.asarray(T_wall, dtype=float) / 2.0) / np.maximum(p_w, 1e-30)
    return lam, lam / length


@dataclass
class SurfaceFlowResult:
    # body-scale diagnostics
    branch: int
    kn_body: float
    re_shock: float               # Re2 = rho_inf V R / mu(T0); nan outside the shock-layer/merged branches
    mach_inf: float
    p_stag: float                 # post-shock stagnation pressure (nan in the free-molecular branch)
    p_inf: float
    phi_sonic: float              # rad (nan where no shock layer is constructed)
    deceleration: float
    gamma_pm: float
    # per patch
    p_w: np.ndarray
    u_e: np.ndarray
    u_eff: np.ndarray             # slip-corrected edge velocity driving the film
    rho_e: np.ndarray
    T_e: np.ndarray
    mu_e: np.ndarray
    mach_e: np.ndarray
    delta_a: np.ndarray           # gas boundary-layer thickness (nan where none is constructed)
    lambda_w: np.ndarray
    kn_local: np.ndarray
    closure: np.ndarray           # CLOSURE_GIRIN / CLOSURE_COUETTE per patch
    flagged: np.ndarray           # True where the closure is outside Girin's validated regime
    tau: np.ndarray
    tau_continuum: np.ndarray
    tau_fm: np.ndarray
    G: np.ndarray

    @property
    def branch_name(self):
        return BRANCH_NAMES[self.branch]

    def closure_fractions(self, areas, windward):
        """Windward area fractions of (Girin closure, Couette slip, Couette free-molecular)."""
        a = np.asarray(areas, dtype=float) * np.asarray(windward, dtype=float)
        total = a.sum()
        if total <= 0.0:
            return [0.0, 0.0, 0.0]
        girin = self.closure == CLOSURE_GIRIN
        fm = (self.closure == CLOSURE_COUETTE) & (self.kn_local >= KN_LOCAL_SLIP)
        slip = (self.closure == CLOSURE_COUETTE) & ~fm
        return [float(a[m].sum() / total) for m in (girin, slip, fm)]


class SurfaceFlow:
    def __init__(self, air=None, rarefied_shear="slip", bridging=None, sigma_v=SIGMA_V, sigma_t=SIGMA_T,
                 gamma_pm=GAMMA_PM, kn_body_shock=KN_BODY_SHOCK, kn_body_fm=KN_BODY_FM):
        if rarefied_shear not in RAREFIED_SHEAR_NAMES:
            raise ValueError("rarefied_shear must be one of {}, got {!r}".format(RAREFIED_SHEAR_NAMES, rarefied_shear))
        if not 1.0 < gamma_pm < 2.0:
            raise ValueError("gamma_pm must be within (1, 2), got {!r}".format(gamma_pm))
        self.air = air or gas.EquilibriumAir()
        self.rarefied_shear, self.bridging = rarefied_shear, bridging or aero.SesamTable()
        self.slip_C, self.sigma_t, self.gamma_pm = (2.0 - sigma_v) / sigma_v, sigma_t, gamma_pm
        self.kn_body_shock, self.kn_body_fm = kn_body_shock, kn_body_fm
        self.last_bins = None

    # -- the gate ---------------------------------------------------------------------------------------------
    def branch_of(self, state):
        """The flow branch of the whole body (module docstring): shock layer, merged, or free molecular."""
        kn = state.kn
        if not np.isfinite(kn) or kn >= self.kn_body_fm or state.ma <= 1.0:
            return BRANCH_FREE_MOLECULAR
        return BRANCH_SHOCK_LAYER if kn < self.kn_body_shock else BRANCH_MERGED

    def edge_table(self, state, stag, radius):
        """Edge quantities on THETA_BINS from the post-shock reservoir: (p_w, u_e, rho_e, T_e, mu_e, a_e, phi_sonic)."""
        fs = state.freestream
        p_inf, q_inf = float(fs.p), 0.5 * fs.rho * state.V ** 2
        p_w, phi_star = wall_pressure(THETA_BINS, stag.p, p_inf, q_inf, self.gamma_pm)
        cols = []
        for p in p_w:
            e = self.air.expand(stag, float(p))
            cols.append((math.sqrt(max(0.0, 2.0 * (stag.h - e.h))), e.rho, e.T, e.mu, e.a))
        u_e, rho_e, T_e, mu_e, a_e = (np.array(c) for c in zip(*cols))
        return p_w, u_e, rho_e, T_e, mu_e, a_e, phi_star

    def evaluate(self, state, theta, radius, rho_liquid, T_wall=None):
        """Per-patch SurfaceFlowResult for the trajectory AeroState, the patch angles theta [rad], the nose radius,
        the film density (for the inertial term of G) and the patch wall temperatures [K]."""
        theta = np.asarray(theta, dtype=float)
        n = theta.size
        fs = state.freestream
        V, p_inf = float(state.V), float(fs.p)
        decel = float(np.linalg.norm(state.a_drag))
        windward = theta <= 0.5 * math.pi
        sin, cos = np.sin(theta), np.cos(theta)
        T_w = np.full(n, 300.0) if T_wall is None else np.asarray(T_wall, dtype=float)
        mu_w = np.array([gas.wilke_viscosity(float(t), {"N2": 0.79, "O2": 0.21}, self.air.molar_masses) for t in np.unique(np.round(T_w, 1))])
        mu_wall = np.interp(T_w, np.unique(np.round(T_w, 1)), mu_w)
        tau_fm = np.where(windward, self.sigma_t * fs.rho * V * V * sin * np.abs(cos), 0.0)
        zero, nan = np.zeros(n), np.full(n, np.nan)
        branch = self.branch_of(state)
        if fs.rho <= 0.0 or V <= 0.0:
            return SurfaceFlowResult(BRANCH_FREE_MOLECULAR, state.kn, float("nan"), state.ma, float("nan"), p_inf,
                                     float("nan"), decel, self.gamma_pm, np.full(n, p_inf), zero, zero, zero, zero, zero,
                                     zero, nan, nan, np.full(n, np.inf), np.full(n, CLOSURE_COUETTE), np.ones(n, dtype=bool),
                                     zero, zero, zero, zero)
        G_inertial = -rho_liquid * decel * sin
        if branch == BRANCH_FREE_MOLECULAR:
            # Schaaf-Chambre from the freestream directly: no compression model, no edge state, no boundary layer
            lam_inf = aero.mean_free_path(fs.rho, fs.m_bar)
            return SurfaceFlowResult(branch, state.kn, float("nan"), state.ma, float("nan"), p_inf, float("nan"), decel,
                                     self.gamma_pm, np.full(n, p_inf), zero, zero, np.full(n, fs.rho), np.full(n, fs.T),
                                     zero, zero, nan, np.full(n, lam_inf), np.full(n, lam_inf / radius),
                                     np.full(n, CLOSURE_COUETTE), np.ones(n, dtype=bool), tau_fm, zero, tau_fm,
                                     np.where(windward, G_inertial, 0.0))
        # shock-layer construction: the reservoir (also the bridge's continuum endpoint in the merged branch)
        stag = self.air.stagnation(fs.rho, fs.T, V)
        re_shock = float(fs.rho * V * radius / stag.mu)
        p_b, u_b, rho_b, T_b, mu_b, a_b, phi_star = self.edge_table(state, stag, radius)
        s_b = radius * THETA_BINS
        delta_b = boundary_layer_thickness(s_b, u_b, mu_b / rho_b)
        self.last_bins = (THETA_BINS, p_b, u_b, rho_b, T_b, mu_b, delta_b)
        th = np.minimum(theta, 0.5 * math.pi)
        at = lambda col: np.interp(th, THETA_BINS, col)
        p_w, u_e, rho_e, T_e, mu_e, a_e, delta_a = (at(c) for c in (p_b, u_b, rho_b, T_b, mu_b, a_b, delta_b))
        p_w = np.where(windward, p_w, p_inf)
        u_e = np.where(windward, u_e, 0.0)
        lam_w, kn_local = wall_knudsen(p_w, T_w, mu_wall, radius)
        with np.errstate(divide="ignore", invalid="ignore"):
            tau_c = np.where(delta_a > 0.0, mu_e * u_e / delta_a, 0.0)
        slip = 1.0 / (1.0 + self.slip_C * kn_local)
        tau_slip = tau_c * slip
        if branch == BRANCH_MERGED:
            # no valid post-shock construction: bridge the continuum and free-molecular limits, Couette closure
            f = self.bridging(state.kn)
            tau = (1.0 - f) * tau_slip + f * tau_fm
            closure = np.full(n, CLOSURE_COUETTE)
            flagged = np.ones(n, dtype=bool)
        else:
            closure = np.where(kn_local < KN_LOCAL_CONTINUUM, CLOSURE_GIRIN, CLOSURE_COUETTE)
            flagged = closure != CLOSURE_GIRIN
            if self.rarefied_shear == "slip":
                tau = np.where(kn_local >= KN_LOCAL_SLIP, tau_fm, tau_slip)
            else:
                f = self.bridging(state.kn)
                tau = (1.0 - f) * tau_slip + f * tau_fm
        tau = np.where(windward & np.isfinite(tau), tau, 0.0)
        closure = np.where(windward, closure, CLOSURE_COUETTE)
        dp_ds = np.gradient(p_b, s_b)                              # from the actual p_w(theta), not a Newtonian formula
        G = np.where(windward, -at(dp_ds) + G_inertial, 0.0)
        mach_e = np.where(a_e > 0.0, u_e / np.maximum(a_e, 1e-300), 0.0)
        u_e, rho_e, mu_e = (np.where(np.isfinite(x), x, 0.0) for x in (u_e, rho_e, mu_e))
        return SurfaceFlowResult(branch, float(state.kn), re_shock, float(state.ma), float(stag.p), p_inf, float(phi_star),
                                 decel, self.gamma_pm, p_w, u_e, u_e * slip, rho_e, T_e, mu_e, mach_e,
                                 np.where(windward, delta_a, np.nan), lam_w, np.where(windward, kn_local, np.inf),
                                 closure, flagged | ~windward, tau, np.where(windward, tau_c, 0.0), tau_fm, G)
