"""The parent's shock layer: the gas between its bow shock and its windward surface, as an engineering field (spec 4.3d).

A transitional or free-molecular follower has no shock of its own and barely disturbs the parent's flow, so the parent's
shock-layer gas is computed alone, once per release, from the parent's free stream at release (`FreeStream`), and then
frozen; the follower is integrated through it (`adapters.WakeAtmosphere`, `adapters.RarefiedSphereSimulator`).

Frame. Parent-centred: x along the parent's velocity (positive forward), y the distance from the axis; spherical
coordinates (r, phi) about the centroid with phi the angle from the stagnation line (phi = 180 - theta, theta the
release angle of `affix`). The parent is the equal-area sphere of radius R. Gas velocities (w_x, w_y) are relative to
the parent: the free stream is (-U, 0) and the gas at the stagnation point is at rest.

The model, piece by piece (each named in the spec as an assumption):
  * Shock shape: Billig (1967), x_s(y) = R + Delta - R_c cot^2(beta) [(1 + y^2 tan^2(beta) / R_c^2)^(1/2) - 1], with
    the vertex radius R_c = R 1.143 exp(0.54 / (M - 1)^1.2) and beta = arcsin(1/M). The standoff Delta is Lighthill's
    (default, at the equilibrium density ratio epsilon across the normal shock) or Billig's R 0.143 exp(3.24 / M^2)
    (`ShockLayerSettings.shock_standoff`). Billig's fit covers Mach 2 to about 25 (`flags`). This is the
    `shock_shape="billig"` sensitivity case and, through the shock angles it implies, the start of the flux-matching
    iteration (each ray's distance balanced alone at Billig's angle, `_balanced_distances`).
  * Flux-matched shock (`shock_shape="flux_matched"`, the default; Asha's option 2): Billig's perfect-gas curvature
    with the equilibrium standoff carried 27 to 73 percent too much gas through the layer, so ahead of the base plane
    the shock distance r_s(phi) on every grid ray is set by mass balance: the layer flux, the integral of
    rho |u_phi| 2 pi r sin(phi) dr from R to r_s, equals the free-stream flux rho_inf U pi (r_s sin(phi))^2 through the
    shock disc of the same lateral radius. The profiles depend on the shock angles, which depend on the shape, so
    the shock distance is found iteratively, started from Billig's shape. The shape is a cubic spline in phi through
    knots every 3 degrees (continuous slope, independent of the field grid), the axis knot held at the chosen
    standoff; its knots are found by a damped (relaxed) Gauss-Newton iteration on the flux residuals of all grid rays,
    with post-shock states from a table of oblique jumps during the iteration and exact jumps for the final field
    (`ShockLayerField._flux_match` says why the plain fixed point of the spec was replaced). The shock angle is
    floored at the Mach angle; behind the base plane the slope relaxes towards the Mach angle (`FluxMatchedShape`).
    `flux_mismatch` (the largest relative mismatch on the grid rays from 10 to 90 degrees; 0.4 percent at 80 km and at
    most 1.5 percent from 48 to 96 km, spec limit 2 percent, `flags["flux_mismatch_above_tol"]` above it) and
    `flux_iterations` record the result.
  * Just behind the shock: equilibrium Rankine-Hugoniot on the normal component U sin(sigma) of the free stream, sigma
    the local shock angle from the slope of x_s, by the damped fixed point of `gas.EquilibriumAir.stagnation`; the
    tangential component U cos(sigma) is unchanged (`oblique_jump`).
  * At the surface: the pressure is the Newtonian-Prandtl-Meyer method (`surface_pressure="newtonian_prandtl_meyer"`,
    the default; Lees and Kubota 1957, Wagner 1959, Kaufman 1963): modified Newtonian, p_N = p_inf + (p_stag - p_inf)
    cos^2(phi), from the stagnation point to the matching point phi_q where its slope equals that of the expansion, and
    a Prandtl-Meyer expansion of the stagnation isentrope from there to the shoulder (`newtonian_prandtl_meyer`).
    Modified Newtonian alone puts the shoulder at p_inf, so its gas is about a tenth of the free-stream density and
    every release near theta = 90 degrees comes out rarefied; the expansion leaves it at tens of p_inf
    (`surface_pressure="modified_newtonian"` is the previous model, a sensitivity variation). The temperature and
    density are then the isentropic equilibrium expansion from the stagnation state to the surface pressure, and the
    gas speed tangent to the surface comes from the total enthalpy, (2 (h_0 - h))^(1/2). Inviscid (no boundary layer).
  * Between: eta = (r - R) / (r_s(phi) - R), 0 at the body and 1 at the shock. p and T linear in eta, density (and
    mean molecular mass) the equilibrium values at that p and T. Velocity: Lighthill's (1957) constant-density solution
    near the nose (phi up to `phi_nose_deg`), evaluated at r_L = R + eta (R_s - R) (stretched to the local layer);
    linear in eta between the surface and post-shock velocities from `phi_nose_deg + phi_blend_deg` on; a cubic
    smoothstep in phi between.
  * Edges: across the shock the field fades to the ambient gas over `shock_thickness_mfp` free-stream mean free paths,
    measured outward along the local shock normal (first order); behind the base plane (x < 0) it fades to zero at
    x = -r0, evaluated there on the shoulder ray (phi = 90 degrees) at the same distance from the axis. There the
    field's phi derivatives are zero, so the pressure-gradient force changes discontinuously at the base plane (the
    weight and the state are continuous). Inside the body the state is the surface state (eta clamped to 0), so it is
    always finite.

Tabulation: the field is computed on an (phi, eta) grid with phi from 0 to 90 degrees (n_phi by n_eta nodes), mirrored
across the axis by three nodes (p, T, rho, m_bar and w_x even in phi, w_y odd) so the splines are symmetric there, and
interpolated with bicubic splines; the Cantera calls happen once per build. The pressure gradient is the analytic
derivative of weight * (p - p_inf) through the pressure spline, the shock-distance spline and the smoothing weights.

Nothing here takes or knows an altitude; it imports `gas` and `aero` only (spec 4.0; a test inspects the signature).
"""
import math
import time
from dataclasses import dataclass, fields

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.optimize import brentq

from reentry_model.shared import aero
from reentry_model.shared.gas import AIR, CP_AIR, T_EQUILIBRATE, T_FLOOR, EquilibriumAir

AVOGADRO_KMOL = 6.02214076e26        # 1/kmol: Cantera's molecular weights are kg/kmol
GAS_CONSTANT = 8314.462618           # J/(kmol K), universal gas constant (Cantera's units)
BILLIG_MACH_RANGE = (2.0, 25.0)      # Mach numbers covered by Billig's sphere correlation
STANDOFF_MODELS = ("lighthill", "billig")
SHAPE_MODELS = ("flux_matched", "billig")
SURFACE_PRESSURE_MODELS = ("newtonian_prandtl_meyer", "modified_newtonian")
ISENTROPE_NODES = 400                # pressures (log-spaced, stagnation to free stream) of the stagnation isentrope table
N_MIRROR = 3                         # grid nodes mirrored across the axis for the splines


@dataclass(frozen=True)
class ShockLayerSettings:
    """Named, overridable choices of the shock-layer model (spec 4.3d)."""
    shock_standoff: str = "lighthill"          # "lighthill" (equilibrium epsilon) or "billig" (perfect gas)
    shock_shape: str = "flux_matched"          # "flux_matched" (mass balance) or "billig" (Billig's shape)
    surface_pressure: str = "newtonian_prandtl_meyer"   # or "modified_newtonian" (the shoulder is then at p_inf)
    phi_nose_deg: float = 40.0                 # Lighthill's velocities up to this angle from the stagnation line
    phi_blend_deg: float = 20.0                # width of the smoothstep blend to the linear-in-eta velocities
    shock_thickness_mfp: float = 3.0           # shock smoothing thickness, free-stream mean free paths
    n_phi: int = 91                            # grid nodes in phi, 0 to 90 degrees
    n_eta: int = 21                            # grid nodes in eta, body to shock
    pressure_gradient: bool = True             # apply the pressure-gradient force to the follower
    follower_wall_T_K: float = 300.0           # follower wall temperature for the drag law (assumed)
    kn_shock_layer_flag: float = 1e-3          # parent Knudsen number above which the shock layer is flagged
    q_change_flag: float = 0.10                # relative change of the parent's dynamic pressure that is flagged

    def __post_init__(self):
        if self.shock_standoff not in STANDOFF_MODELS:
            raise ValueError(f"shock_standoff must be one of {STANDOFF_MODELS}, got {self.shock_standoff!r}")
        if self.shock_shape not in SHAPE_MODELS:
            raise ValueError(f"shock_shape must be one of {SHAPE_MODELS}, got {self.shock_shape!r}")
        if self.surface_pressure not in SURFACE_PRESSURE_MODELS:
            raise ValueError(f"surface_pressure must be one of {SURFACE_PRESSURE_MODELS}, got "
                             f"{self.surface_pressure!r}")
        for name in ("phi_nose_deg", "phi_blend_deg", "shock_thickness_mfp", "follower_wall_T_K",
                     "kn_shock_layer_flag", "q_change_flag"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be a finite positive number, got {v!r}")
        if self.phi_nose_deg + self.phi_blend_deg > 90.0:
            raise ValueError(f"phi_nose_deg + phi_blend_deg must not exceed 90 degrees, got "
                             f"{self.phi_nose_deg} + {self.phi_blend_deg}")
        for name in ("n_phi", "n_eta"):
            v = getattr(self, name)
            if not (isinstance(v, int) and v >= 4):
                raise ValueError(f"{name} must be an integer of at least 4 grid nodes, got {v!r}")

    def non_default(self) -> dict:
        """The fields that differ from the defaults, for run names."""
        default = ShockLayerSettings()
        return {f.name: getattr(self, f.name) for f in fields(self)
                if getattr(self, f.name) != getattr(default, f.name)}


@dataclass(frozen=True)
class FreeStream:
    """The parent's free stream at release (all SI)."""
    rho: float      # kg/m3
    T: float        # K
    p: float        # Pa
    U: float        # parent's air-relative speed, m/s
    mfp: float      # mean free path, m
    m_bar: float    # mean molecular mass, kg


@dataclass(frozen=True)
class PostShock:
    """Equilibrium state just behind an oblique shock; velocities relative to the parent, in the shock's frame."""
    rho: float      # kg/m3
    T: float        # K
    p: float        # Pa
    h: float        # J/kg (Cantera basis)
    s: float        # J/(kg K)
    m_bar: float    # kg
    w_n: float      # normal velocity behind the shock, m/s (into the layer)
    w_t: float      # tangential velocity, m/s (U cos(sigma), unchanged)


def _m_bar(g):
    return float(g.mean_molecular_weight) / AVOGADRO_KMOL


def oblique_jump(air, fs, sigma_rad):
    """Equilibrium Rankine-Hugoniot jump on the normal component fs.U sin(sigma) of the free stream (the damped fixed
    point of `EquilibriumAir.stagnation`, equilibrated only when T_inf + w_n^2 / (2 c_p) exceeds T_EQUILIBRATE, and no
    jump at a normal Mach number of 1.1 or less, as there); the tangential component fs.U cos(sigma) is unchanged."""
    g = air.gas
    T_inf = max(fs.T, T_FLOOR)
    g.TDX = T_inf, fs.rho, AIR
    p_inf, h_inf = float(g.P), float(g.enthalpy_mass)
    v_n, v_t = fs.U * math.sin(sigma_rad), fs.U * math.cos(sigma_rad)
    if v_n / float(g.sound_speed) <= 1.1:
        return PostShock(fs.rho, T_inf, p_inf, h_inf, float(g.entropy_mass), _m_bar(g), v_n, v_t)
    equilibrate = T_inf + 0.5 * v_n * v_n / CP_AIR > T_EQUILIBRATE
    eps, change = 0.1, math.inf
    for _ in range(3000):
        u2 = eps * v_n
        g.HPX = h_inf + 0.5 * (v_n * v_n - u2 * u2), p_inf + fs.rho * v_n * (v_n - u2), AIR
        if equilibrate:
            g.equilibrate("HP")
        eps_new = fs.rho / float(g.density)
        change = abs(eps_new - eps)
        if change < 1e-10:
            break
        eps = 0.5 * (eps + eps_new)
    else:
        raise RuntimeError(f"oblique shock iteration did not converge in 3000 iterations (rho {fs.rho:.3e} kg/m3, "
                           f"T {fs.T:.1f} K, normal velocity {v_n:.0f} m/s, last change {change:.2e})")
    return PostShock(float(g.density), float(g.T), float(g.P), float(g.enthalpy_mass), float(g.entropy_mass),
                     _m_bar(g), u2, v_t)


# ---------------------------------------------------------------- Lighthill's constant-density nose solution

@dataclass(frozen=True)
class Lighthill:
    """Lighthill's (1957) constant-density flow between a sphere of radius R and a concentric shock of radius R_s:
    Stokes stream function psi = sin^2(phi) g(r), g(r) = A r^2 + B / r + K r^4, K = -U (1 - eps)^2 / (10 eps R_s^2)."""
    eps: float
    U: float
    R: float
    R_s: float
    A: float
    B: float
    K: float

    def g(self, r):
        return self.A * r * r + self.B / r + self.K * r ** 4

    def dg(self, r):
        return 2.0 * self.A * r - self.B / (r * r) + 4.0 * self.K * r ** 3

    def velocity(self, r, phi):
        """(u_r, u_phi) relative to the body: u_r = 2 g cos(phi) / r^2, u_phi = -g' sin(phi) / r."""
        return 2.0 * self.g(r) * math.cos(phi) / (r * r), -self.dg(r) * math.sin(phi) / r


def _lighthill_coefficients(eps, U, R_s):
    """(A, B, K) for a shock of radius R_s: g(R_s) = -eps U R_s^2 / 2 and g'(R_s) = -U R_s (2x2 linear in A, B)."""
    K = -U * (1.0 - eps) ** 2 / (10.0 * eps * R_s * R_s)
    c1 = -0.5 * eps * U * R_s * R_s - K * R_s ** 4          # A R_s^2 + B / R_s
    c2 = -U * R_s - 4.0 * K * R_s ** 3                       # 2 A R_s - B / R_s^2
    A = (c1 + c2 * R_s) / (3.0 * R_s * R_s)
    B = R_s * (c1 - A * R_s * R_s)
    return A, B, K


def lighthill(eps, U, R):
    """Lighthill's solution for density ratio eps = rho_inf / rho_2, speed U and body radius R: R_s is the root of
    g(R) = 0 closest to the body (for small eps a second, spurious root lies further out), found by scanning outward
    from R for the first sign change and refining with brentq."""
    if not (0.0 < eps < 1.0 and U > 0.0 and R > 0.0):
        raise ValueError(f"need 0 < eps < 1 and positive U and R, got eps {eps!r}, U {U!r}, R {R!r}")

    def body_g(R_s):
        A, B, K = _lighthill_coefficients(eps, U, R_s)
        return A * R * R + B / R + K * R ** 4

    grid = R * (1.0 + np.linspace(1e-6, 3.0, 30001))
    vals = np.array([body_g(x) for x in grid])
    change = np.nonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if not change.size:
        raise ValueError(f"no shock radius found within 4 R for eps {eps!r}")
    i = int(change[0])
    R_s = brentq(body_g, grid[i], grid[i + 1], xtol=1e-14 * R, rtol=1e-15)
    return Lighthill(float(eps), float(U), float(R), float(R_s), *_lighthill_coefficients(eps, U, R_s))


def lighthill_standoff_approx(eps):
    """Hayes and Probstein's approximation of Lighthill's standoff, as a fraction of R_s: eps / (1 + sqrt(8 eps / 3))."""
    return eps / (1.0 + math.sqrt(8.0 * eps / 3.0))


# ---------------------------------------------------------------- shock shapes

def billig_standoff(R, M):
    """Billig's (1967) sphere standoff, perfect gas with gamma 1.4: R 0.143 exp(3.24 / M^2)."""
    return R * 0.143 * math.exp(3.24 / M ** 2)


def billig_vertex_radius(R, M):
    """Billig's shock vertex radius of curvature for a sphere: R 1.143 exp(0.54 / (M - 1)^1.2)."""
    return R * 1.143 * math.exp(0.54 / (M - 1.0) ** 1.2)


class _ShockCurve:
    """Geometry shared by the shock shapes, from x_s(y), its first two derivatives, R, M and the axis standoff."""

    @property
    def beta(self):
        """Mach angle arcsin(1/M) [rad]."""
        return math.asin(1.0 / self.M)

    def sigma(self, y):
        """Local shock angle [rad] between the shock surface and the free-stream direction: pi/2 on the axis, falling
        towards the Mach angle far out, and never below it."""
        return max(math.atan2(1.0, -self.dx_s(y)), self.beta)

    def r_s(self, phi):
        """Distance [m] from the centroid to the shock along the ray at phi (root of x_s(r sin phi) = r cos phi)."""
        f = lambda r: self.x_s(r * math.sin(phi)) - r * math.cos(phi)
        hi = 2.0 * (self.R + self.standoff)
        while f(hi) > 0.0:
            hi *= 2.0
        return brentq(f, 0.0, hi, xtol=1e-13 * self.R, rtol=1e-15)

    def inside(self, x, y):
        return x < self.x_s(y)

    def normal_distance(self, x, y):
        """Signed distance [m] to the shock, positive inside, first order along the local shock normal:
        (x_s(y) - x) / (1 + x_s'(y)^2)^(1/2)."""
        return (self.x_s(y) - x) / math.sqrt(1.0 + self.dx_s(y) ** 2)

    def normal_distance_gradient(self, x, y):
        """(d/dx, d/dy) of `normal_distance`."""
        d1, d2 = self.dx_s(y), self.d2x_s(y)
        s = 1.0 / math.sqrt(1.0 + d1 * d1)
        return -s, d1 * s - (self.x_s(y) - x) * d1 * d2 * s ** 3


@dataclass(frozen=True)
class ShockShape(_ShockCurve):
    """Billig's hyperbola with vertex radius from Mach M and the given standoff, in the parent frame (x forward)."""
    R: float
    M: float
    standoff: float

    @property
    def R_c(self):
        return billig_vertex_radius(self.R, self.M)

    def _u(self, y):
        return (y * math.tan(self.beta) / self.R_c) ** 2

    def x_s(self, y):
        """Shock position along x at distance y from the axis [m]."""
        u = self._u(y)
        return self.R + self.standoff - (y * y / self.R_c) / (math.sqrt(1.0 + u) + 1.0)

    def dx_s(self, y):
        return -(y / self.R_c) / math.sqrt(1.0 + self._u(y))

    def d2x_s(self, y):
        return -1.0 / (self.R_c * (1.0 + self._u(y)) ** 1.5)


def _mirrored(phi, a, odd=False):
    """(phi, a) extended across the axis by N_MIRROR nodes (a even in phi, or odd)."""
    m = min(N_MIRROR, len(phi) - 1)
    return (np.concatenate([-phi[m:0:-1], phi]),
            np.concatenate([(-1.0 if odd else 1.0) * np.asarray(a)[m:0:-1], np.asarray(a)]))


class FluxMatchedShape(_ShockCurve):
    """The flux-matched shock (spec 4.3d): ahead of the base plane, the curve through the distances r_s at the knots
    phi_k from 0 to 90 degrees (`r_s_spline`, a cubic spline in phi mirrored across the axis, so its slope is
    continuous); behind it, x_s'(y) relaxes exponentially from its value at the 90-degree ray towards -cot(beta) (the
    Mach angle), with the length chosen so that x_s'' is continuous there as well. ValueError if the distance from the
    axis does not grow with phi (the curve x_s(y) would not be single-valued)."""
    N_DENSE = 721                    # samples of the spline in phi from which x_s(y) is interpolated

    def __init__(self, R, M, standoff, phi_nodes, r_s_nodes):
        self.R, self.M, self.standoff = float(R), float(M), float(standoff)
        self.r_s_spline = CubicSpline(*_mirrored(np.asarray(phi_nodes, float), r_s_nodes))
        ph = np.linspace(0.0, 0.5 * math.pi, self.N_DENSE)
        r, dr = self.r_s_spline(ph), self.r_s_spline(ph, 1)
        c, s = np.cos(ph), np.sin(ph)
        y = r * s
        dx, dy = dr * c - r * s, dr * s + r * c
        if not np.all(dy[1:] > 0.0):
            raise ValueError("the shock's distance from the axis must grow with phi (dr_s/dphi too negative)")
        slope_90 = float(dx[-1] / dy[-1])
        self._inner = CubicSpline(y, r * c, bc_type=((1, 0.0), (1, slope_90)))
        self._y90 = float(y[-1])
        cot_b = 1.0 / math.tan(self.beta)
        self._cot_b = cot_b
        self._excess = max(slope_90 + cot_b, 0.0)           # slope above the Mach-angle slope at the 90-degree ray
        curv = float(self._inner(self._y90, 2))
        self._L = self._excess / -curv if curv < 0.0 and self._excess > 0.0 else self._y90

    def _behind(self, y):
        return math.exp(-(y - self._y90) / self._L)

    def x_s(self, y):
        if y <= self._y90:
            return float(self._inner(y))
        q = y - self._y90
        return float(self._inner(self._y90)) - self._cot_b * q + self._excess * self._L * (1.0 - self._behind(y))

    def dx_s(self, y):
        if y <= self._y90:
            return float(self._inner(y, 1))
        return -self._cot_b + self._excess * self._behind(y)

    def d2x_s(self, y):
        if y <= self._y90:
            return float(self._inner(y, 2))
        return -self._excess / self._L * self._behind(y)


# ---------------------------------------------------------------- the surface pressure

@dataclass(frozen=True)
class SurfacePressure:
    """Result of `newtonian_prandtl_meyer`: the surface pressure p_b [Pa] at each requested angle, the matching
    point's angle phi_q [rad] and Mach number mach_q (both nan when `no_match`), and the two flags."""
    p_b: np.ndarray
    phi_q: float
    mach_q: float
    no_match: bool      # no supersonic part, or no root of the slope balance: Newtonian everywhere
    below_p_inf: bool   # the expansion would take the pressure below p_inf somewhere: held at p_inf there


def newtonian_prandtl_meyer(p_stag, p_inf, p, rho, w, a, phi):
    """Surface pressure of a sphere by the Newtonian-Prandtl-Meyer method (spec 4.3d, "State at the body surface").

    `p, rho, w, a` tabulate the stagnation isentrope (pressure strictly decreasing from p_stag to p_inf, density,
    gas speed, equilibrium sound speed); `phi` are angles from the stagnation point [rad], 0 to pi/2. Modified
    Newtonian, p_N = p_inf + (p_stag - p_inf) cos^2(phi), holds up to the matching point phi_q; beyond it the turning
    of the Prandtl-Meyer expansion sets the pressure, nu(p_b) = nu(p_q) + (phi - phi_q), because the body inclination
    delta = pi/2 - phi falls as fast as phi rises. With M = w / a, d nu = -(M^2 - 1)^(1/2) dp / (rho w^2) from the
    sonic point (trapezoidal in p), valid for any gas. At the matching point the pressures agree by construction and
    the slopes dp/d(delta) must too: 2 ((p - p_inf)(p_stag - p))^(1/2) for p_N against rho w^2 / (M^2 - 1)^(1/2) for
    the expansion, which is infinite at the sonic point. Their difference is negative there, so the root sought is
    the first sign change walking down in pressure from the sonic point, refined by brentq on linear interpolants of
    rho, w and a in log p. No supersonic node or no sign change: Newtonian everywhere and `no_match`; where the
    turning exceeds that to p_inf, p_b = p_inf and `below_p_inf`."""
    p, rho, w, a, phi = (np.asarray(x, float) for x in (p, rho, w, a, phi))
    newton = p_inf + (p_stag - p_inf) * np.cos(phi) ** 2
    mach = w / a
    supersonic = np.nonzero(mach > 1.0)[0]
    if p_stag <= p_inf or not supersonic.size:
        return SurfacePressure(newton, math.nan, math.nan, True, False)
    i0 = int(supersonic[0])
    lp = np.log(p)
    lps, rs, ws, as_ = lp[i0:], rho[i0:], w[i0:], a[i0:]
    if i0 > 0:                                     # the sonic point, between nodes i0 - 1 and i0, becomes the first node
        f = (1.0 - mach[i0 - 1]) / (mach[i0] - mach[i0 - 1])
        lerp = lambda x: x[i0 - 1] + f * (x[i0] - x[i0 - 1])
        lps, rs, ws, as_ = (np.concatenate([[v], x]) for v, x in
                            ((lerp(lp), lps), (lerp(rho), rs), (lerp(w), ws), (lerp(w), as_)))
    ps = np.exp(lps)
    m2m1 = np.maximum((ws / as_) ** 2 - 1.0, 0.0)
    with np.errstate(divide="ignore"):
        g_pm = np.where(m2m1 > 0.0, rs * ws * ws / np.sqrt(m2m1), np.inf)
    g_n = 2.0 * np.sqrt(np.maximum((ps - p_inf) * (p_stag - ps), 0.0))
    rate = np.sqrt(m2m1) / (rs * ws * ws)                          # -d nu / dp
    nu = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * (ps[:-1] - ps[1:]))])
    diff = g_n - g_pm
    crossing = np.nonzero((diff[:-1] < 0.0) & (diff[1:] >= 0.0))[0]
    if not crossing.size:
        return SurfacePressure(newton, math.nan, math.nan, True, False)
    j = int(crossing[0])

    def balance(x):
        r_, w_, a_ = (float(np.interp(x, lps[::-1], v[::-1])) for v in (rs, ws, as_))
        p_ = math.exp(x)
        return (2.0 * math.sqrt(max((p_ - p_inf) * (p_stag - p_), 0.0))
                - r_ * w_ * w_ / math.sqrt(max((w_ / a_) ** 2 - 1.0, 1e-300)))

    lo, hi = lps[j + 1], lps[j]                    # log p falls with the index; hi may be the sonic point (infinite)
    lq = lo if diff[j + 1] == 0.0 else brentq(balance, lo, lo + (1.0 - 1e-9) * (hi - lo), xtol=1e-14, rtol=1e-14)
    p_q = math.exp(lq)
    mach_q = float(np.interp(lq, lps[::-1], (ws / as_)[::-1]))
    phi_q = math.acos(math.sqrt((p_q - p_inf) / (p_stag - p_inf)))
    nu_q = float(np.interp(lq, lps[::-1], nu[::-1]))
    beyond = phi > phi_q
    target = nu_q + (phi - phi_q)
    beyond_inf = beyond & (target > nu[-1])
    p_b = newton.copy()
    p_b[beyond] = np.exp(np.interp(target[beyond], nu, lps))
    p_b[beyond_inf] = p_inf
    return SurfacePressure(p_b, phi_q, mach_q, False, bool(beyond_inf.any()))


# ---------------------------------------------------------------- the tabulated field

@dataclass
class FieldState:
    """Field gas at one point: rho [kg/m3], T [K], p [Pa], m_bar [kg], velocity relative to the parent (w_x forward,
    w_y away from the axis) [m/s], and the field's share `weight` in [0, 1]."""
    rho: float
    T: float
    p: float
    m_bar: float
    w_x: float
    w_y: float
    weight: float


def _smoothstep(u):
    """Cubic smoothstep S(u) and its derivative S'(u), clamped to [0, 1]."""
    if u <= 0.0:
        return 0.0, 0.0
    if u >= 1.0:
        return 1.0, 0.0
    return u * u * (3.0 - 2.0 * u), 6.0 * u * (1.0 - u)


_KEYS = ("rho", "T", "p", "m_bar", "w_x", "w_y")    # the FieldState fields, in order
FLUX_TOL = 0.02                      # flux mismatch required on every ray from 10 to 90 degrees (spec 4.3d)
FLUX_TARGET = 0.002                  # the iteration stops below this, or where the knots allow no further improvement
FLUX_MAX_ITER = 30
FLUX_DSIGMA = 1e-4                   # rad, step for the derivative of the flux residuals in the shock angle
FLUX_MAX_STEP = 0.25                 # largest change of a ray's gap in one iteration, as a fraction of the gap
FLUX_KNOT_DEG = 3.0                  # degrees between the knots of the flux-matched shape
FLUX_PHI_MIN = math.radians(10.0)    # rays checked for the mismatch: 10 to 90 degrees
JUMP_TABLE_NODES = 61                # oblique jumps tabulated in sin(sigma) for the iteration
JUMP_TABLE_MN_MIN = 1.2              # lowest normal Mach number in the table (exact jumps below it)
MW_TABLE_T_STEP = 50.0               # K, temperature spacing of the equilibrium molecular-weight table
MW_TABLE_NP = 25                     # pressure nodes (uniform in log p) of that table


class _JumpTable:
    """Oblique jumps at JUMP_TABLE_NODES shock angles, uniform in sin(sigma) from a normal Mach number of
    JUMP_TABLE_MN_MIN to 1, interpolated by cubic splines, for the flux-matching iteration (one exact jump costs
    several milliseconds). Exact jumps below the table."""

    def __init__(self, air, fs):
        self.air, self.fs = air, fs
        air.gas.TDX = max(fs.T, T_FLOOR), fs.rho, AIR
        self.s_min = min(JUMP_TABLE_MN_MIN * float(air.gas.sound_speed) / fs.U, 1.0)
        s = np.linspace(self.s_min, 1.0, JUMP_TABLE_NODES)
        jumps = [oblique_jump(air, fs, math.asin(x)) for x in s]
        self._splines = {k: CubicSpline(s, [getattr(j, k) for j in jumps]) for k in ("rho", "T", "p", "w_n")}

    def __call__(self, sigma):
        s = math.sin(sigma)
        if s < self.s_min:
            return oblique_jump(self.air, self.fs, sigma)
        v = {k: float(sp(s)) for k, sp in self._splines.items()}
        return PostShock(v["rho"], v["T"], v["p"], math.nan, math.nan, math.nan, v["w_n"],
                         self.fs.U * math.cos(sigma))


class _MolecularWeightTable:
    """Equilibrium mean molecular weight [kg/kmol] on a (T, log p) grid covering the layer, bicubic, for the
    flux-matching iteration; constant (frozen air) when the field is not equilibrated."""

    def __init__(self, air, T_max, p_min, p_max, equilibrate):
        g = air.gas
        g.TPX = 300.0, 101325.0, AIR
        self.frozen = float(g.mean_molecular_weight)
        self.equilibrate = equilibrate
        if not equilibrate:
            return
        T = np.arange(T_FLOOR, T_max + 2.0 * MW_TABLE_T_STEP, MW_TABLE_T_STEP)
        lp = np.linspace(math.log(p_min), math.log(p_max), MW_TABLE_NP)
        W = np.empty((T.size, lp.size))
        for i, t in enumerate(T):
            for j, x in enumerate(lp):
                g.TPX = t, math.exp(x), AIR
                if t > T_EQUILIBRATE:
                    g.equilibrate("TP")
                W[i, j] = float(g.mean_molecular_weight)
        self._spline = RectBivariateSpline(T, lp, W, kx=3, ky=3)

    def __call__(self, p, T):
        if not self.equilibrate:
            return np.full(np.shape(p), self.frozen)
        return self._spline.ev(T, np.log(p))


class ShockLayerField:
    """The frozen shock-layer field of one release (module docstring). Build with `ShockLayerField.build`.

    Attributes: fs, R, r0, settings, mach, eps, lighthill (the nose solution), standoff (used on the axis),
    standoff_billig (for comparison), shape (`FluxMatchedShape` or Billig's `ShockShape`), flags, flux_mismatch (the
    largest |layer flux / inflow - 1| on the grid rays from 10 to 90 degrees, for either shape), flux_iterations (0
    for Billig's shape), build_time_s, surface_match_phi_deg and surface_match_mach (the Newtonian-Prandtl-Meyer
    matching point; nan for modified Newtonian or no match), shoulder_pressure_ratio (surface pressure at phi = 90
    degrees over p_inf, either model).
    """

    @classmethod
    def build(cls, fs, R, r0, settings=ShockLayerSettings(), air=None):
        """Field around a sphere of radius R [m] in free stream `fs`, faded over r0 [m] behind the base plane."""
        started = time.perf_counter()
        air = air if air is not None else EquilibriumAir()
        self = cls()
        self.fs, self.R, self.r0, self.settings = fs, float(R), float(r0), settings
        self.mach = fs.U / aero.speed_of_sound(fs.T, fs.m_bar)
        self.eps = fs.rho / oblique_jump(air, fs, 0.5 * math.pi).rho
        self.lighthill = lighthill(self.eps, fs.U, self.R)
        self.standoff_billig = billig_standoff(self.R, self.mach)
        self.standoff = (self.lighthill.R_s - self.R if settings.shock_standoff == "lighthill"
                         else self.standoff_billig)
        self.flags = {
            "mach_outside_billig": not BILLIG_MACH_RANGE[0] <= self.mach <= BILLIG_MACH_RANGE[1],
            "kn_parent_above_flag": fs.mfp / (2.0 * self.R) > settings.kn_shock_layer_flag,
            "npm_no_match": False, "surface_below_p_inf": False,
        }
        self.surface_match_phi_deg = self.surface_match_mach = math.nan
        self._air = air
        stag = air.stagnation(fs.rho, fs.T, fs.U)
        self._s0, self._h0, self._p_stag = float(air.gas.entropy_mass), stag.h, stag.p   # stagnation() leaves the gas there
        self._equilibrate = max(fs.T, T_FLOOR) + 0.5 * fs.U ** 2 / CP_AIR > T_EQUILIBRATE
        self._phi = np.linspace(0.0, 0.5 * math.pi, settings.n_phi)
        self._eta = np.linspace(0.0, 1.0, settings.n_eta)
        self._p_surface = self._surface_pressure()
        self._surface = [self._surface_state(ph, p_b) for ph, p_b in zip(self._phi, self._p_surface)]
        self.shoulder_pressure_ratio = float(self._p_surface[-1] / fs.p)

        billig = ShockShape(self.R, self.mach, self.standoff)
        r_s = np.array([billig.r_s(ph) for ph in self._phi])
        self.shape, self.flux_iterations = billig, 0
        if settings.shock_shape == "flux_matched":
            r_s = self._flux_match(r_s, stag)
        exact = lambda sigma: oblique_jump(air, fs, sigma)
        rays = [self._ray(i, self.shape.sigma(r_s[i] * math.sin(ph)), exact)
                for i, ph in enumerate(self._phi)]
        self.flux_mismatch = self._worst(self._flux_residuals(rays, r_s)[0])
        self.flags["flux_mismatch_above_tol"] = self.flux_mismatch > FLUX_TOL

        self._splines = {}
        for k in _KEYS:
            phi_ext, table = _mirrored(self._phi, np.array([ray[k] for ray in rays]), odd=(k == "w_y"))
            self._splines[k] = RectBivariateSpline(phi_ext, self._eta, table, kx=3, ky=3)
        if isinstance(self.shape, FluxMatchedShape):
            self._r_s = self.shape.r_s_spline
        else:
            self._r_s = CubicSpline(*_mirrored(self._phi, r_s))
        self._r_s_shoulder = float(r_s[-1])
        self.build_time_s = time.perf_counter() - started
        return self

    # ---- construction ----

    def _surface_pressure(self):
        """Surface pressure [Pa] at every node of `self._phi`: modified Newtonian, or the Newtonian-Prandtl-Meyer method
        on the stagnation isentrope (equilibrium air at the stagnation entropy, 400 pressures log-spaced from p_stag
        to p_inf; gas speed from the total enthalpy, sound speed a^2 = dp/drho along it). Sets
        `surface_match_phi_deg`, `surface_match_mach` and the two flags."""
        g, fs, phi = self._air.gas, self.fs, self._phi
        if self.settings.surface_pressure == "modified_newtonian" or self._p_stag <= fs.p:
            return fs.p + (self._p_stag - fs.p) * np.cos(phi) ** 2
        p = np.geomspace(self._p_stag, fs.p, ISENTROPE_NODES)
        rho, h = np.empty(p.size), np.empty(p.size)
        for i, pk in enumerate(p):
            g.SPX = self._s0, pk, AIR
            if self._equilibrate:
                g.equilibrate("SP")
            rho[i], h[i] = float(g.density), float(g.enthalpy_mass)
        w = np.sqrt(np.maximum(2.0 * (self._h0 - h), 0.0))
        r = newtonian_prandtl_meyer(self._p_stag, fs.p, p, rho, w, np.sqrt(np.gradient(p, rho)), phi)
        self.surface_match_phi_deg = math.degrees(r.phi_q)
        self.surface_match_mach = r.mach_q
        self.flags["npm_no_match"], self.flags["surface_below_p_inf"] = r.no_match, r.below_p_inf
        return r.p_b

    def _surface_state(self, ph, p_b):
        """(p_b, T_b, w_b) on the body at phi, given the surface pressure p_b: isentropic equilibrium expansion from
        the stagnation state, tangential speed from the total enthalpy."""
        g = self._air.gas
        c, s = math.cos(ph), math.sin(ph)
        g.SPX = self._s0, p_b, AIR
        if self._equilibrate:
            g.equilibrate("SP")
        speed = math.sqrt(max(2.0 * (self._h0 - float(g.enthalpy_mass)), 0.0))
        return p_b, float(g.T), np.array([-speed * s, speed * c])

    def _flux_match(self, r_s, stag):
        """Ray distances r_s [m] of the flux-matched shock (spec 4.3d), starting from `r_s` (Billig's); sets `shape`
        and `flux_iterations`.

        The shape is a cubic spline in phi through knots every FLUX_KNOT_DEG degrees (mirrored across the axis, the
        axis knot held at the chosen standoff), so it is smooth and independent of the field grid. Its knot values are
        found by a damped Gauss-Newton iteration on the residuals of all grid rays, F_i = layer flux / inflow - 1;
        F_i depends on the ray's distance directly and on the shape through its local shock angle sigma_i (from the
        spline's value and slope). A plain relaxed fixed point (each ray solved for its own distance at the current
        angles, then relaxed and smoothed) is unstable at the shoulder: there the layer flux depends so strongly on
        the shock angle that the balance acts as an equation for the slope, and short-wavelength corrections are
        amplified (measured 2026-10-29: it diverged, or with heavy smoothing stalled at 5 percent, or was stable for
        one setting only). The relaxation here is the step halving until the residuals fall, with each gap changed by
        at most FLUX_MAX_STEP of itself per step. dF/dr is exact at fixed profiles and dF/dsigma a difference of two
        ray evaluations. Post-shock states come from a table of oblique jumps and densities from a table of
        equilibrium molecular weights, so an iteration costs tens of milliseconds; the final field is rebuilt from
        exact jumps and Cantera states."""
        fs, phi = self.fs, self._phi
        knots = np.linspace(0.0, 0.5 * math.pi, int(round(90.0 / FLUX_KNOT_DEG)) + 1)
        m = knots.size
        eye = np.eye(m)
        unit = [CubicSpline(*_mirrored(knots, eye[k])) for k in range(m)]
        B = np.column_stack([u(phi) for u in unit])             # ray values from knot values
        Bd = np.column_stack([u(phi, 1) for u in unit])         # ray slopes from knot values
        jumps = _JumpTable(self._air, fs)
        weights = _MolecularWeightTable(self._air, 1.05 * stag.T, 0.5 * fs.p, 1.05 * stag.p, self._equilibrate)
        q = np.interp(knots, phi, self._balanced_distances(r_s, jumps, weights))    # start: each ray balanced alone
        q[0] = self.R + self.standoff
        beta = math.asin(1.0 / self.mach)
        c, s = np.cos(phi), np.sin(phi)
        off_axis = phi > 0.0

        def angles(q):
            r, dr = B @ q, Bd @ q
            a, b = dr * s + r * c, r * s - dr * c              # shock tangent (dy, -dx) per unit phi
            h = a * a + b * b
            sig = np.arctan2(a, b)
            free = sig > beta                                  # floored at the Mach angle: no dependence there
            return (np.maximum(sig, beta), np.where(free, (b * c - a * s) / h, 0.0),
                    np.where(free, (b * s + a * c) / h, 0.0))

        def residuals(q, dsig=0.0):
            sig = angles(q)[0]
            return self._flux_residuals([self._ray(i, sig[i] + dsig, jumps, weights) for i in range(phi.size)], B @ q)

        def valid(q):
            try:
                FluxMatchedShape(self.R, self.mach, self.standoff, knots, q)
            except ValueError:
                return False
            return bool(np.all(B @ q > self.R))

        F, dF_dr = residuals(q)
        while self._worst(F) >= FLUX_TARGET and self.flux_iterations < FLUX_MAX_ITER:
            _, sig_r, sig_dr = angles(q)
            dF_dsig = (residuals(q, FLUX_DSIGMA)[0] - F) / FLUX_DSIGMA
            J = dF_dr[:, None] * B + dF_dsig[:, None] * (sig_r[:, None] * B + sig_dr[:, None] * Bd)
            step = np.zeros(m)
            step[1:] = np.linalg.lstsq(J[off_axis][:, 1:], -F[off_axis], rcond=None)[0]   # axis knot fixed
            norm = float(np.linalg.norm(F))
            lam = min(1.0, FLUX_MAX_STEP / float(np.max(np.abs(step[1:]) / (q[1:] - self.R))))
            improved = False
            while lam >= 1.0 / 1024.0:
                trial = q + lam * step
                if valid(trial):
                    F_trial, dF_trial = residuals(trial)
                    if np.linalg.norm(F_trial) < norm * (1.0 - 1e-6):
                        improved = True
                        break
                lam *= 0.5
            if not improved:                                   # the least-squares minimum of this knot spacing
                break
            q, F, dF_dr = trial, F_trial, dF_trial
            self.flux_iterations += 1
        self.shape = FluxMatchedShape(self.R, self.mach, self.standoff, knots, q)
        return B @ q

    def _balanced_distances(self, r_s, jump, weights):
        """Ray distances [m] at which each ray's layer flux balances its inflow, with the shock angle the shape `r_s`
        (Billig's) has there: the start of the flux-matching iteration. With the profiles fixed in eta the balance
        2 G (R I0 + G I1) = a (R + G)^2 (see `_flux_residuals`) is a quadratic in the gap G, and its positive root is
        used (Billig's distance where there is none). Billig's own distances carry up to 3 times too much gas at
        high Mach number (0.78 m against 0.27 m at the shoulder at 70 km); started there, the Gauss-Newton iteration
        reached the shoulder only through shapes whose shock angle jumps from 24 to 77 degrees between 85 and 90
        degrees (a flux mismatch of 0.28 at 70, 71, 72, 87 and 94 km, measured 2026-10-31). Requiring the shock
        angle not to rise along phi did not rescue that start (the line search then found no step: mismatch 1.19 at
        70 km) and is not a property of the solutions at low Mach number (3.7 at 48 km), so it is not imposed."""
        R, fs = self.R, self.fs
        billig = self.shape
        out = np.array(r_s, float)
        for i, ph in enumerate(self._phi):
            if ph == 0.0:
                continue
            s = math.sin(ph)
            ray = self._ray(i, billig.sigma(r_s[i] * s), jump, weights)
            f = np.abs(-ray["w_x"] * s + ray["w_y"] * math.cos(ph)) * ray["rho"]
            I0 = CubicSpline(self._eta, f).integrate(0.0, 1.0)
            I1 = CubicSpline(self._eta, f * self._eta).integrate(0.0, 1.0)
            a = fs.rho * fs.U * s
            A, B = 2.0 * I1 - a, 2.0 * R * (I0 - a)
            disc = B * B + 4.0 * A * a * R * R
            roots = [g for g in ((-B + e * math.sqrt(disc)) / (2.0 * A) for e in (1.0, -1.0))
                     if g > 0.0] if disc >= 0.0 and A != 0.0 else []
            if roots:
                out[i] = R + min(roots)
        return out

    def _ray(self, i, sigma, jump, weights=None):
        """Profiles across the layer on ray i, {key: array over eta}, behind a shock of local angle sigma with the
        post-shock state `jump(sigma)`. Densities and molecular masses are Cantera's equilibrium values at each (p, T),
        or from `weights` (a `_MolecularWeightTable`) when given."""
        settings, g, L = self.settings, self._air.gas, self.lighthill
        ph, et = self._phi[i], self._eta
        c, s = math.cos(ph), math.sin(ph)
        p_b, T_b, wb = self._surface[i]
        ps = jump(sigma)
        ws = np.array([-ps.w_t * math.cos(sigma) - ps.w_n * math.sin(sigma),
                       ps.w_t * math.sin(sigma) - ps.w_n * math.cos(sigma)])
        p, T = (1.0 - et) * p_b + et * ps.p, (1.0 - et) * T_b + et * ps.T
        if weights is not None:
            W = weights(p, T)
            rho, m_bar = p * W / (GAS_CONSTANT * T), W / AVOGADRO_KMOL
        else:
            rho, m_bar = np.empty(et.size), np.empty(et.size)
            for j in range(et.size):
                g.TPX = T[j], p[j], AIR
                if self._equilibrate and T[j] > T_EQUILIBRATE:
                    g.equilibrate("TP")
                rho[j], m_bar[j] = float(g.density), _m_bar(g)
        w = np.outer(1.0 - et, wb) + np.outer(et, ws)
        mix = _smoothstep((ph - math.radians(settings.phi_nose_deg)) / math.radians(settings.phi_blend_deg))[0]
        if mix < 1.0:
            r_L = L.R + et * (L.R_s - L.R)
            u_r = 2.0 * L.g(r_L) * c / (r_L * r_L)
            u_phi = -L.dg(r_L) * s / r_L
            w = (1.0 - mix) * np.column_stack([u_r * c - u_phi * s, u_r * s + u_phi * c]) + mix * w
        return {"rho": rho, "T": T, "p": p, "m_bar": m_bar, "w_x": w[:, 0], "w_y": w[:, 1]}

    def _flux_residuals(self, rays, r_s):
        """(F, dF/dr): F_i = layer flux / inflow - 1 on every ray (0 on the axis), and its derivative in r_i with the
        profiles fixed.

        With the profiles fixed in eta, the layer flux on a ray of gap G = r_s - R is 2 pi sin(phi) G (R I0 + G I1),
        I0 and I1 the integrals over eta of rho |u_phi| and rho |u_phi| eta (cubic splines in eta), and the inflow is
        pi rho_inf U sin^2(phi) (R + G)^2."""
        R, fs = self.R, self.fs
        F, dF = np.zeros(len(rays)), np.ones(len(rays))
        for i, ph in enumerate(self._phi):
            if ph == 0.0:
                continue
            s, c = math.sin(ph), math.cos(ph)
            f = np.abs(-rays[i]["w_x"] * s + rays[i]["w_y"] * c) * rays[i]["rho"]
            I0 = CubicSpline(self._eta, f).integrate(0.0, 1.0)
            I1 = CubicSpline(self._eta, f * self._eta).integrate(0.0, 1.0)
            G = r_s[i] - R
            a = fs.rho * fs.U * s
            F[i] = 2.0 * G * (R * I0 + G * I1) / (a * (R + G) ** 2) - 1.0
            dF[i] = 2.0 * ((R * I0 + 2.0 * G * I1) * (R + G) - 2.0 * G * (R * I0 + G * I1)) / (a * (R + G) ** 3)
        return F, dF

    def _worst(self, F):
        """Largest |F| on the rays from 10 to 90 degrees (the flux mismatch)."""
        return float(np.max(np.abs(F[self._phi >= FLUX_PHI_MIN - 1e-12])))

    # ---- coordinates and weights ----

    def _coords(self, x, y):
        """(phi, eta_raw, d eta_raw/dx, d eta_raw/dy, d phi/dx, d phi/dy) of a point. Behind the base plane the
        shoulder ray (phi = 90 degrees) at the same distance from the axis, so the phi derivatives are zero there and
        the field's pressure gradient changes discontinuously at x = 0 (the weight is continuous)."""
        R = self.R
        if x < 0.0:
            gap = self._r_s_shoulder - R
            return 0.5 * math.pi, (y - R) / gap, 0.0, 1.0 / gap, 0.0, 0.0
        r = math.hypot(x, y)
        phi = math.atan2(y, x)
        rs = float(self._r_s(phi))
        gap = rs - R
        eta = (r - R) / gap
        if r == 0.0:
            return phi, eta, 0.0, 0.0, 0.0, 0.0
        de_dr, de_dphi = 1.0 / gap, -(r - R) * float(self._r_s(phi, 1)) / (gap * gap)
        dphi_dx, dphi_dy = -y / (r * r), x / (r * r)
        return (phi, eta, de_dr * x / r + de_dphi * dphi_dx, de_dr * y / r + de_dphi * dphi_dy,
                dphi_dx, dphi_dy)

    def _weight(self, x, y):
        """(weight, d weight/dx, d weight/dy): the shock smoothing times the base-plane fade."""
        d = self.settings.shock_thickness_mfp * self.fs.mfp
        ws, dws = _smoothstep((self.shape.normal_distance(x, y) + d) / d)
        if ws == 0.0:
            return 0.0, 0.0, 0.0
        wb, dwb = _smoothstep((x + self.r0) / self.r0) if x < 0.0 else (1.0, 0.0)
        if wb == 0.0:
            return 0.0, 0.0, 0.0
        if dws:
            nx, ny = self.shape.normal_distance_gradient(x, y)
            gx, gy = dws / d * nx * wb, dws / d * ny * wb
        else:
            gx = gy = 0.0
        return ws * wb, gx + ws * dwb / self.r0, gy

    # ---- queries ----

    def _ev(self, key, phi, eta, dphi=0, deta=0):
        return float(self._splines[key].ev(phi, eta, dx=dphi, dy=deta))

    def state(self, x, y):
        """FieldState at (x, y) in the parent frame; the free stream with weight 0 where the field does not reach."""
        w = self._weight(x, y)[0]
        if w == 0.0:
            fs = self.fs
            return FieldState(fs.rho, fs.T, fs.p, fs.m_bar, 0.0, 0.0, 0.0)
        phi, eta = self._coords(x, y)[:2]
        eta = min(max(eta, 0.0), 1.0)
        return FieldState(*(self._ev(k, phi, eta) for k in _KEYS), w)

    def surface_state(self, theta_deg):
        """FieldState on the body at release angle theta in [0, 180] degrees (eta = 0 at phi = 180 - theta), weight 1;
        ValueError outside that range. Releases with theta below 90 degrees lie behind the base plane, where this is the
        shoulder state's tabulated extension, not a modelled state."""
        theta = float(theta_deg)
        if not 0.0 <= theta <= 180.0:
            raise ValueError(f"theta_deg must lie in [0, 180] degrees, got {theta_deg!r}")
        phi = math.radians(180.0 - theta)
        return FieldState(*(self._ev(k, phi, 0.0) for k in _KEYS), 1.0)

    def excess_pressure_gradient(self, x, y):
        """(d/dx, d/dy) of weight * (p_field - p_inf) [Pa/m]: includes the smoothed jump across the shock and the
        base-plane fade, and is exactly zero where the field does not reach. Behind the base plane the field is the
        shoulder ray's, so its own gradient has no x part there (see `_coords`)."""
        w, dwx, dwy = self._weight(x, y)
        if w == 0.0:
            return 0.0, 0.0
        phi, eta, de_dx, de_dy, dphi_dx, dphi_dy = self._coords(x, y)
        if eta > 1.0:                                 # beyond the shock the post-shock state is held (weight fades)
            eta, de_dx, de_dy = 1.0, 0.0, 0.0
        elif eta < 0.0:                               # inside the body: the gradient at the surface, so the force on
            eta = 0.0                                 # a follower released on the surface does not depend on round-off
        p = self._ev("p", phi, eta)
        p_phi, p_eta = self._ev("p", phi, eta, dphi=1), self._ev("p", phi, eta, deta=1)
        dpx, dpy = p_phi * dphi_dx + p_eta * de_dx, p_phi * dphi_dy + p_eta * de_dy
        excess = p - self.fs.p
        return dwx * excess + w * dpx, dwy * excess + w * dpy
