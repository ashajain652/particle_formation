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
    (`ShockLayerSettings.shock_standoff`). Billig's fit covers Mach 2 to about 25 (`flags`).
  * Just behind the shock: equilibrium Rankine-Hugoniot on the normal component U sin(sigma) of the free stream, sigma
    the local shock angle from the slope of x_s, by the damped fixed point of `gas.EquilibriumAir.stagnation`; the
    tangential component U cos(sigma) is unchanged (`oblique_jump`).
  * At the surface: modified Newtonian pressure p_b = p_inf + (p_stag - p_inf) cos^2(phi), temperature and density by
    isentropic equilibrium expansion from the stagnation state to p_b, and a gas speed tangent to the surface from the
    total enthalpy, (2 (h_0 - h))^(1/2). Inviscid (no boundary layer).
  * Between: eta = (r - R) / (r_s(phi) - R), 0 at the body and 1 at the shock. p and T linear in eta, density (and
    mean molecular mass) the equilibrium values at that p and T. Velocity: Lighthill's (1957) constant-density solution
    near the nose (phi up to `phi_nose_deg`), evaluated at r_L = R + eta (R_s - R) (stretched to the local layer);
    linear in eta between the surface and post-shock velocities from `phi_nose_deg + phi_blend_deg` on; a cubic
    smoothstep in phi between.
  * Edges: across the shock the field fades to the ambient gas over `shock_thickness_mfp` free-stream mean free paths,
    measured outward along the local shock normal (first order); behind the base plane (x < 0) it fades to zero at
    x = -r0, evaluated there on the shoulder ray (phi = 90 degrees) at the same distance from the axis. Inside the body
    the state is the surface state (eta clamped to 0), so it is always finite.

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

from reentry_model import aero
from reentry_model.gas import AIR, CP_AIR, T_EQUILIBRATE, T_FLOOR, EquilibriumAir

AVOGADRO_KMOL = 6.02214076e26        # 1/kmol: Cantera's molecular weights are kg/kmol
BILLIG_MACH_RANGE = (2.0, 25.0)      # Mach numbers covered by Billig's sphere correlation
STANDOFF_MODELS = ("lighthill", "billig")
N_MIRROR = 3                         # grid nodes mirrored across the axis for the splines


@dataclass(frozen=True)
class ShockLayerSettings:
    """Named, overridable choices of the shock-layer model (spec 4.3d)."""
    shock_standoff: str = "lighthill"          # "lighthill" (equilibrium epsilon) or "billig" (perfect gas)
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


# ---------------------------------------------------------------- Billig's shock shape

def billig_standoff(R, M):
    """Billig's (1967) sphere standoff, perfect gas with gamma 1.4: R 0.143 exp(3.24 / M^2)."""
    return R * 0.143 * math.exp(3.24 / M ** 2)


def billig_vertex_radius(R, M):
    """Billig's shock vertex radius of curvature for a sphere: R 1.143 exp(0.54 / (M - 1)^1.2)."""
    return R * 1.143 * math.exp(0.54 / (M - 1.0) ** 1.2)


@dataclass(frozen=True)
class ShockShape:
    """Billig's hyperbola with vertex radius from Mach M and the given standoff, in the parent frame (x forward)."""
    R: float
    M: float
    standoff: float

    @property
    def R_c(self):
        return billig_vertex_radius(self.R, self.M)

    @property
    def beta(self):
        return math.asin(1.0 / self.M)

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

    def sigma(self, y):
        """Local shock angle [rad] between the shock surface and the free-stream direction: pi/2 on the axis, falling
        towards beta far out."""
        return math.atan2(1.0, -self.dx_s(y))

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
        (x_s(y) - x) sin(sigma(y))."""
        return (self.x_s(y) - x) / math.sqrt(1.0 + self.dx_s(y) ** 2)

    def normal_distance_gradient(self, x, y):
        """(d/dx, d/dy) of `normal_distance`."""
        d1, d2 = self.dx_s(y), self.d2x_s(y)
        s = 1.0 / math.sqrt(1.0 + d1 * d1)
        return -s, d1 * s - (self.x_s(y) - x) * d1 * d2 * s ** 3


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


class ShockLayerField:
    """The frozen shock-layer field of one release (module docstring). Build with `ShockLayerField.build`."""

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
        self.shape = ShockShape(self.R, self.mach, self.standoff)
        self.flags = {
            "mach_outside_billig": not BILLIG_MACH_RANGE[0] <= self.mach <= BILLIG_MACH_RANGE[1],
            "kn_parent_above_flag": fs.mfp / (2.0 * self.R) > settings.kn_shock_layer_flag,
        }

        g = air.gas
        stag = air.stagnation(fs.rho, fs.T, fs.U)
        s0, h0 = float(g.entropy_mass), stag.h                   # stagnation() leaves the gas at the stagnation state
        equilibrate = max(fs.T, T_FLOOR) + 0.5 * fs.U ** 2 / CP_AIR > T_EQUILIBRATE
        phi = np.linspace(0.0, 0.5 * math.pi, settings.n_phi)
        eta = np.linspace(0.0, 1.0, settings.n_eta)
        table = {k: np.empty((phi.size, eta.size)) for k in _KEYS}
        r_s = np.empty(phi.size)
        nose, blend = math.radians(settings.phi_nose_deg), math.radians(settings.phi_blend_deg)
        L = self.lighthill
        for i, ph in enumerate(phi):
            c, s = math.cos(ph), math.sin(ph)
            p_b = fs.p + (stag.p - fs.p) * c * c                 # modified Newtonian surface pressure
            g.SPX = s0, p_b, AIR
            if equilibrate:
                g.equilibrate("SP")
            T_b = float(g.T)
            speed_b = math.sqrt(max(2.0 * (h0 - float(g.enthalpy_mass)), 0.0))
            wb = np.array([-speed_b * s, speed_b * c])           # tangent to the surface, away from the stagnation point
            r_s[i] = self.shape.r_s(ph)
            sig = self.shape.sigma(r_s[i] * s)
            ps = oblique_jump(air, fs, sig)
            ws = np.array([-ps.w_t * math.cos(sig) - ps.w_n * math.sin(sig),
                           ps.w_t * math.sin(sig) - ps.w_n * math.cos(sig)])
            mix = _smoothstep((ph - nose) / blend)[0]
            for j, et in enumerate(eta):
                p, T = (1.0 - et) * p_b + et * ps.p, (1.0 - et) * T_b + et * ps.T
                g.TPX = T, p, AIR
                if equilibrate and T > T_EQUILIBRATE:
                    g.equilibrate("TP")
                w = (1.0 - et) * wb + et * ws
                if mix < 1.0:
                    u_r, u_phi = L.velocity(L.R + et * (L.R_s - L.R), ph)
                    w = (1.0 - mix) * np.array([u_r * c - u_phi * s, u_r * s + u_phi * c]) + mix * w
                for k, v in zip(_KEYS, (float(g.density), T, p, _m_bar(g), w[0], w[1])):
                    table[k][i, j] = v

        m = min(N_MIRROR, phi.size - 1)
        phi_ext = np.concatenate([-phi[m:0:-1], phi])
        def mirrored(a, odd=False):
            return np.concatenate([(-1.0 if odd else 1.0) * a[m:0:-1], a])
        self._splines = {k: RectBivariateSpline(phi_ext, eta, mirrored(table[k], odd=(k == "w_y")), kx=3, ky=3)
                         for k in _KEYS}
        self._r_s = CubicSpline(phi_ext, mirrored(r_s))
        self._r_s_shoulder = float(r_s[-1])
        self.build_time_s = time.perf_counter() - started
        return self

    # ---- coordinates and weights ----

    def _coords(self, x, y):
        """(phi, eta_raw, d eta_raw/dx, d eta_raw/dy, d phi/dx, d phi/dy) of a point; behind the base plane the
        shoulder ray at the same distance from the axis."""
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
        """FieldState on the body at release angle theta (eta = 0 at phi = 180 - theta), weight 1."""
        phi = math.radians(180.0 - float(theta_deg))
        return FieldState(*(self._ev(k, phi, 0.0) for k in _KEYS), 1.0)

    def excess_pressure_gradient(self, x, y):
        """(d/dx, d/dy) of weight * (p_field - p_inf) [Pa/m]: includes the smoothed jump across the shock and the
        base-plane fade, and is exactly zero where the field does not reach."""
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
