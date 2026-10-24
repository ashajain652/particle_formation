"""Independent demise of a molten AA7075 droplet shed from the descending parent: trajectory, heating and
multi-component vaporisation, run until the droplet reaches a steady temperature and velocity.

What is reused, and what is new
-------------------------------
* Trajectory: reentry_model.trajectory.Simulator, the rotating-Earth 3-DOF integrator verified against DRAMA/SESAM
  (ECEF state, J2 gravity, Coriolis and centrifugal terms, geodetic conversions, ground/escape events).
  `DropletSimulator` subclasses it and extends the state with the droplet temperature and the masses of the
  volatile elements; only the drag coefficient is replaced, because the validated SESAM drag table is
  hypersonic-only and wrong at low speed (factor 1.94 at Kn 0.1, measured 2026-10-14). The default (drag_law="singh",
  from 2026-10-22) is the general sphere correlation of Singh et al., AIAA Journal 2021 (singh_drag.py), fitted across
  continuum, transitional and free-molecule flow at all Mach numbers. drag_law="blend" keeps the earlier composed law
  of prototype/work-2026-10-14-particle-settling/particle_trajectory.py (Schaaf & Chambre bridged to the continuum
  sphere, blended on the speed ratio into Stokes-Cunningham) for comparison. The two are identical in fast
  free-molecule flow; Singh gives 5-13 % more drag at low speed with the droplet at air temperature (its Epstein
  limit assumes fully diffuse reflection, where Millikan's Cunningham fit reflects real surfaces) and differs by up
  to about 35 % either way in the fast transitional corner.

* Heating: reentry_model.heating.PhysicsHeating (Fay-Riddell continuum + Matting bridging + Lees/cosine
  distribution), the model verified against SESAM, is used where it adds physics the free-molecule formula lacks:
  hypersonic flow (Mach >= 6, blended down to Mach 3) that is not free-molecular (Knudsen number lambda/D <= 1,
  blended up to 10), while the wall is below the recovery temperature. Elsewhere the heat exchange is the general
  Schaaf & Chambre free-molecule sphere at any speed ratio (droplet_heating_sweep.fm_sphere_heat), reduced towards
  the continuum Ranz-Marshall conduction by Sherman's interpolation q = q_fm / (1 + h_fm0 / h_c). Measured
  2026-10-21: in free-molecule flow PhysicsHeating's flux alpha rho V (h_s - h_w) is the high-speed limit of
  Schaaf & Chambre, 3 % low at Mach 20 for a 1000 K wall but 25 % low at Mach 10 for a 2000 K wall (it subtracts
  the equilibrium-air wall enthalpy, which counts dissociation that free-molecule re-emission does not), and it
  returns zero, never a cooling flux, once the wall is above recovery. heating_mode="mach" restores PhysicsHeating
  wherever Mach >= 6 regardless of rarefaction, for sensitivity.

* Vaporisation (new): each volatile element i of the well-mixed droplet (Al, Zn, Mg; Cu and Cr inert) leaves at
      mdot_i = 1 / (1/mdot_HK,i + 1/mdot_diff,i)
  with mdot_HK,i = alpha_e A x_i p_sat,i(T) sqrt(M_i / (2 pi R T)) the Hertz-Knudsen-Langmuir kinetic limit
  (Raoult's law for the partial pressure, no metal vapour in the ambient air), and mdot_diff,i the continuum
  diffusion limit with the Stefan-flow factor ln(1/(1-y))/y, y = sum_i x_i p_sat,i / p_ambient. When the summed
  vapour pressure exceeds the ambient pressure (y >= 1) the droplet boils: the vapour pushes the air away, the
  diffusion limit is removed, and evaporation proceeds at the kinetic rate. In rarefied flow the kinetic limit
  governs anyway. Each element's latent heat is drawn from the droplet's energy.

Final state: the first instant at which the speed is within STEADY_SPEED_TOL of the local terminal velocity and
the temperature within STEADY_T_TOL of the local equilibrium temperature (the temperature at which convective,
radiative and evaporative power balance at the current flight state). Runs also end on the ground, on complete
vaporisation, or at T_MAX.

Before a run, the release Weber number We = rho V^2 d / sigma is checked against WE_CRIT; an unstable droplet is
skipped, not run.
"""
import math
import os
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SETTLING = os.path.join(REPO, "prototype", "work-2026-10-14-particle-settling")
for p in (REPO, SETTLING):
    if p not in sys.path:
        sys.path.insert(0, p)

import particle_trajectory as pt                                   # composed drag law (verified 2026-10-14)
import droplet_heating_sweep as dhs                                # Schaaf & Chambre heat (verified 2026-10-14)
import sphere_sweep
from reentry_model import aero, earth, heating
from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.trajectory import InitialState, Settings, Simulator

import singh_drag                                                 # Singh et al. (2021) general sphere drag

R_GAS = 8.314462618
SIGMA_SB = 5.670374419e-8

# ---------------------------------------------------------------- molten AA7075

RHO_LIQUID = 2400.0      # kg/m3 (proto3 AA7075.json liquid.density)
CP_LIQUID = 1180.0       # J/(kg K) (liquid Al, NIST-JANAF)
T_RELEASE = 908.0        # K, AA7075 liquidus
EPSILON = 0.10           # liquid Al total emissivity (2026-10-14 sweep baseline)
T_BG = 250.0             # K, radiative background
ALPHA_HEAT = 0.8         # thermal accommodation: PhysicsHeating's default, used in both heat branches
ALPHA_EVAP = 1.0         # evaporation coefficient of a clean liquid metal surface
SIGMA_LIQUIDUS = 0.86    # N/m at the liquidus (proto3 AA7075.json liquid.surfaceTension)
DSIGMA_DT = -1.55e-4     # N/(m K), liquid Al; only used to follow We along the flight
WE_CRIT = 12.0           # critical Weber number for low-Ohnesorge droplets (Pilch & Erdman 1987)


@dataclass(frozen=True)
class Element:
    name: str
    M: float            # kg/mol
    mass_fraction: float
    volatile: bool
    T_boil: float       # K, normal boiling point (1 atm)
    H_vap: float        # J/mol at the boiling point


# AA7075 nominal composition (ASM Handbook Vol. 2): Zn 5.6, Mg 2.5, Cu 1.6, Cr 0.23 wt %, balance Al.
ELEMENTS = (
    Element("Al", 0.026982, 0.9007, True, 2743.0, 294.0e3),
    Element("Zn", 0.065380, 0.0560, True, 1180.0, 115.3e3),
    Element("Mg", 0.024305, 0.0250, True, 1363.0, 127.4e3),
    Element("Cu", 0.063546, 0.0160, False, 2835.0, 300.4e3),
    Element("Cr", 0.051996, 0.0023, False, 2944.0, 339.5e3),
)
VOLATILE = tuple(e for e in ELEMENTS if e.volatile)
INERT = tuple(e for e in ELEMENTS if not e.volatile)


def p_sat(el, T):
    """Pure-element vapour pressure [Pa] over the liquid. Al: Alcock, Itkin & Horrigan (1984),
    log10(p/atm) = 5.911 - 16211/T. Zn, Mg: Clausius-Clapeyron through the normal boiling point with the
    enthalpy of vaporisation there (accuracy a few tens of percent far below the boiling point)."""
    if el.name == "Al":
        return 101325.0 * 10.0 ** (5.911 - 16211.0 / T)
    return 101325.0 * math.exp(-el.H_vap / R_GAS * (1.0 / T - 1.0 / el.T_boil))


def surface_tension(T):
    return max(SIGMA_LIQUIDUS + DSIGMA_DT * (T - T_RELEASE), 0.05)


def weber(rho, V, d, T):
    return rho * V * V * d / surface_tension(T)


# ---------------------------------------------------------------- release state from the parent

PARENT_JSON = os.path.join(REPO, "Generic_Satellite Reentry", "output", "dmf_output.json")
PARENT_OBJECT = "Main Body"


def parent_release(parent, h_km):
    """Parent state at altitude h_km on its descending branch (linear in altitude), as an InitialState."""
    b = parent.branch
    for a, c in zip(b, b[1:]):
        if a["altitude"] >= h_km >= c["altitude"]:
            f = (a["altitude"] - h_km) / (a["altitude"] - c["altitude"])
            lerp = sphere_sweep.lerp
            t = lerp(a["time"], c["time"], f)
            return InitialState(
                velocity=lerp(a["velocity"], c["velocity"], f) * 1e3, altitude=h_km * 1e3,
                flight_path=math.radians(lerp(a["path"], c["path"], f)),
                heading=math.radians(sphere_sweep.angle_lerp(a["heading"], c["heading"], f) % 360.0),
                lat=math.radians(lerp(a["lat"], c["lat"], f)),
                lon=math.radians(sphere_sweep.wrap_longitude(sphere_sweep.angle_lerp(a["lon"], c["lon"], f))),
                epoch=parent.epoch + sphere_sweep.timedelta(seconds=t))
    raise ValueError("altitude {} km is outside the parent's descending branch".format(h_km))


# ---------------------------------------------------------------- the droplet

class _MassBody:
    """Simulator's body interface; the droplet's mass lives in the extended state, so this is only used by the
    parent class's own drag path, which DropletSimulator does not call."""

    def __init__(self, m):
        self.m = m

    def mass(self, t):
        return self.m


def smoothstep(x, lo, hi):
    if x <= lo:
        return 0.0
    if x >= hi:
        return 1.0
    u = (x - lo) / (hi - lo)
    return u * u * (3.0 - 2.0 * u)


_MU_W = np.polynomial.legendre.leggauss(24)
MU_WINDWARD = 0.5 * (_MU_W[0] + 1.0)          # cos(theta) on [0, 1]
W_WINDWARD = 0.5 * _MU_W[1]
THETA_WINDWARD = np.arccos(MU_WINDWARD)


class DropletSimulator(Simulator):
    """Simulator with state (r, v, T, m_Al, m_Zn, m_Mg); see the module docstring."""

    STEADY_SPEED_TOL = 0.01
    STEADY_T_TOL = 1.0         # K
    T_MAX = 3600.0             # s
    VAPORISED = 1e-3           # remaining fraction of the initial mass below which the droplet counts as gone
    HYPER_MA = (3.0, 6.0)      # Mach window over which PhysicsHeating takes over
    KN_FM = (1.0, 10.0)        # Knudsen window (lambda / D, PhysicsHeating's convention) over which it hands back

    def __init__(self, initial, d0_um, atmosphere=None, epsilon=EPSILON, alpha_heat=ALPHA_HEAT,
                 alpha_evap=ALPHA_EVAP, physics_heating=None, heating_mode="regime", drag_law="singh"):
        d0 = d0_um * 1e-6
        m0 = RHO_LIQUID * math.pi * d0 ** 3 / 6.0
        self.d0, self.m0 = d0, m0                 # before the parent constructor, which builds the state vector
        self.m_inert = m0 * sum(e.mass_fraction for e in INERT)
        super().__init__(initial, _MassBody(m0), atmosphere or US76TableAtmosphere(), None, None,
                         Settings(diameter=d0, gravity="j2", rotating_frame=True, t_max=self.T_MAX))
        self.eps, self.alpha_heat, self.alpha_evap = epsilon, alpha_heat, alpha_evap
        self.heating = physics_heating or heating.PhysicsHeating(accommodation=alpha_heat)
        if heating_mode not in ("regime", "mach"):
            raise ValueError("heating_mode must be 'regime' or 'mach'")
        self.heating_mode = heating_mode
        if drag_law not in ("blend", "singh"):
            raise ValueError("drag_law must be 'blend' or 'singh'")
        self.drag_law = drag_law

    # -- state layout: 0-2 r, 3-5 v, 6 T, 7.. volatile masses
    def initial_state_vector(self):
        y = super().initial_state_vector()
        return np.concatenate([y, [T_RELEASE], [self.m0 * e.mass_fraction for e in VOLATILE]])

    def unpack(self, y):
        r, v, T = y[:3], y[3:6], float(y[6])
        m_vol = np.maximum(np.asarray(y[7:7 + len(VOLATILE)], dtype=float), 0.0)
        m = float(m_vol.sum()) + self.m_inert
        d = (6.0 * m / (math.pi * RHO_LIQUID)) ** (1.0 / 3.0)
        return r, v, T, m_vol, m, d

    def mole_fractions(self, m_vol):
        n = [mi / e.M for mi, e in zip(m_vol, VOLATILE)] + \
            [self.m0 * e.mass_fraction / e.M for e in INERT]
        tot = sum(n)
        return np.array(n[:len(VOLATILE)]) / tot

    # -- physics at one state
    def flow(self, t, r, v):
        h, lat, lon = earth.ecef_to_geodetic(r)
        fs = self.atmosphere.state(t, min(max(h, 0.0), self.settings.escape_altitude), lat, lon)
        e, n, u = earth.enu_basis(lat, lon)
        wind = fs.wind_enu[0] * e + fs.wind_enu[1] * n + fs.wind_enu[2] * u
        v_rel = v - wind
        return h, lat, lon, fs, v_rel, float(np.linalg.norm(v_rel))

    def q_exchange(self, fs, V, T, d):
        """Free-molecule exchange at any speed ratio (Schaaf & Chambre), reduced towards continuum conduction by
        Sherman's interpolation with the Ranz-Marshall Nusselt number [W, into the droplet]."""
        q_fm = dhs.fm_sphere_heat(fs.rho, fs.T, V, T, d, self.alpha_heat)
        c = math.sqrt(2.0 * pt.R_AIR * fs.T)
        h_fm0 = self.alpha_heat * (2.0 + 0.5 * dhs.ZETA_AIR) * pt.R_AIR * fs.rho * c / (2.0 * math.sqrt(math.pi))
        mu = pt.viscosity(fs.T)
        k_air = mu * 1004.5 / 0.71
        re = fs.rho * V * d / mu
        h_c = (2.0 + 0.6 * math.sqrt(re) * 0.71 ** (1.0 / 3.0)) * k_air / d
        return q_fm / (1.0 + h_fm0 / h_c)

    def q_physics(self, fs, V, T, d):
        """reentry_model.heating.PhysicsHeating integrated over the windward hemisphere at a uniform wall T [W]."""
        ma = aero.mach(V, fs.T, fs.m_bar)
        state = _Aero(fs, V, aero.knudsen(fs.rho, fs.m_bar, d), ma)
        res = self.heating.evaluate(state, THETA_WINDWARD, np.full(THETA_WINDWARD.size, T), 0.5 * d)
        return 2.0 * math.pi * (0.5 * d) ** 2 * float(np.dot(W_WINDWARD, res.q_conv))

    def heat_in(self, fs, V, T, d):
        """Net convective power into the droplet [W] and the PhysicsHeating weight used."""
        q_slow = self.q_exchange(fs, V, T, d)
        w = smoothstep(aero.mach(V, fs.T, fs.m_bar), *self.HYPER_MA)
        if self.heating_mode == "regime":       # free-molecule flow: the exact Schaaf & Chambre exchange instead
            kn = aero.knudsen(fs.rho, fs.m_bar, d)
            w *= 1.0 - smoothstep(math.log10(kn), math.log10(self.KN_FM[0]), math.log10(self.KN_FM[1]))
        if w > 0.0:
            q_ph = self.q_physics(fs, V, T, d)
            if q_ph <= 0.0:          # wall at or above recovery: PhysicsHeating cannot cool, use the exchange model
                w = 0.0
            return w * q_ph + (1.0 - w) * q_slow, w
        return q_slow, 0.0

    def evaporation(self, fs, V, T, m_vol, d):
        """Mass-loss rate of each volatile element [kg/s] (positive = leaving) and the summed vapour fraction."""
        A = math.pi * d * d
        x = self.mole_fractions(m_vol)
        ps = np.array([xi * p_sat(e, T) for xi, e in zip(x, VOLATILE)])
        y = float(ps.sum()) / fs.p
        mdot = np.zeros(len(VOLATILE))
        if y >= 1.0:
            stefan = math.inf                                   # boiling: no diffusion resistance
        else:
            stefan = (-math.log1p(-y) / y) if y > 1e-12 else 1.0
        T_f = 0.5 * (T + fs.T)
        mu = pt.viscosity(fs.T)
        re = fs.rho * V * d / mu
        D = 2.0e-5 * (T_f / 300.0) ** 1.75 * (101325.0 / fs.p)  # metal vapour in air, order-of-magnitude estimate
        sh = 2.0 + 0.6 * math.sqrt(re) * (mu / (fs.rho * D)) ** (1.0 / 3.0)
        for i, (e, p_i) in enumerate(zip(VOLATILE, ps)):
            if m_vol[i] <= 0.0 or p_i <= 0.0:
                continue
            hk = self.alpha_evap * A * p_i * math.sqrt(e.M / (2.0 * math.pi * R_GAS * T))
            if math.isinf(stefan):
                mdot[i] = hk
            else:
                diff = sh * math.pi * d * D * (e.M * p_i / (R_GAS * T_f)) * stefan
                mdot[i] = hk * diff / (hk + diff)
        return mdot, y

    def power_balance(self, fs, V, T, m_vol, d):
        q_conv, w = self.heat_in(fs, V, T, d)
        q_rad = self.eps * SIGMA_SB * math.pi * d * d * (T ** 4 - T_BG ** 4)
        mdot, y = self.evaporation(fs, V, T, m_vol, d)
        q_evap = float(sum(md * e.H_vap / e.M for md, e in zip(mdot, VOLATILE)))
        return q_conv, q_rad, q_evap, mdot, y, w

    def cd(self, V, fs, d, T):
        """Drag coefficient of the selected law, with particle_trajectory's regime diagnostics (Kn = 2 lambda/D,
        Re, Ma, s, blend weight) for the history either way."""
        cd, dg = pt.drag_coefficient(max(V, 1e-30), fs.rho, fs.T, fs.p, d, T / fs.T, "blend")
        if self.drag_law == "singh":
            M = max(V, 1e-30) / math.sqrt(pt.GAMMA * pt.R_AIR * fs.T)
            cd, _ = singh_drag.drag_coefficient(M, fs.rho * max(V, 1e-30) * d / pt.viscosity(fs.T), T / fs.T)
        return cd, dg

    def drag(self, fs, V, T, d, m):
        cd, dg = self.cd(V, fs, d, T)
        return cd, dg, 0.5 * fs.rho * V * cd * (math.pi * d * d / 4.0) / m

    def rhs(self, t, y):
        r, v, T, m_vol, m, d = self.unpack(y)
        h, lat, lon, fs, v_rel, V = self.flow(t, r, v)
        cd, dg, k = self.drag(fs, V, T, d, m)
        a = earth.gravity(r, self.settings.gravity) - 2.0 * np.cross(self.omega, v) \
            - np.cross(self.omega, np.cross(self.omega, r)) - k * v_rel
        q_conv, q_rad, q_evap, mdot, _, _ = self.power_balance(fs, V, T, m_vol, d)
        dT = (q_conv - q_rad - q_evap) / (m * CP_LIQUID)
        dm = np.where(m_vol > 0.0, -mdot, 0.0)
        return np.concatenate([v, a, [dT], dm])

    # -- steady-state diagnostics
    def g_eff(self, r):
        return float(np.linalg.norm(earth.gravity(r, self.settings.gravity)
                                    - np.cross(self.omega, np.cross(self.omega, r))))

    def terminal_velocity(self, fs, r, T, d, m):
        w = m * self.g_eff(r) * (1.0 - fs.rho / RHO_LIQUID)
        area = math.pi * d * d / 4.0

        def res(lv):
            v = math.exp(lv)
            cd, _ = self.cd(v, fs, d, T)
            return 0.5 * fs.rho * v * v * cd * area - w
        return math.exp(brentq(res, math.log(1e-12), math.log(1e4), xtol=1e-12, rtol=1e-12))

    def equilibrium_temperature(self, fs, V, m_vol, d):
        def net(T):
            q_conv, q_rad, q_evap, _, _, _ = self.power_balance(fs, V, T, m_vol, d)
            return q_conv - q_rad - q_evap
        lo, hi = 100.0, 6000.0
        if net(lo) < 0.0:
            return lo
        return brentq(net, lo, hi, xtol=1e-3)

    # -- run
    def diagnostics(self, t, y):
        r, v, T, m_vol, m, d = self.unpack(y)
        h, lat, lon, fs, v_rel, V = self.flow(t, r, v)
        cd, dg, k = self.drag(fs, V, T, d, m)
        q_conv, q_rad, q_evap, mdot, yv, w = self.power_balance(fs, V, T, m_vol, d)
        row = dict(t_s=t, h_km=h / 1e3, lat_deg=math.degrees(lat), lon_deg=math.degrees(lon), V_ms=V,
                   T_K=T, d_um=d * 1e6, mass_kg=m, mass_frac=m / self.m0,
                   q_conv_W=q_conv, q_rad_W=q_rad, q_evap_W=q_evap, w_physics_heating=w,
                   vapour_over_ambient=yv, decel_g=k * V / 9.80665, C_D=cd, Kn=dg["Kn"], Ma=dg["Ma"],
                   speed_ratio=dg["s"], blend_w=dg["w"], We=weber(fs.rho, V, d, T),
                   rho_air=fs.rho, p_air=fs.p, T_air=fs.T)
        for e, mi, md in zip(VOLATILE, m_vol, mdot):
            row["m_{}_kg".format(e.name)] = float(mi)
            row["mdot_{}_kgs".format(e.name)] = float(md)
        return row, (r, v, T, m_vol, m, d, fs, V)

    def run(self):
        """Advance in geometrically growing chunks; after each, test for the steady final state."""
        s = self.settings
        atol = [s.atol_position] * 3 + [s.atol_velocity] * 3 + [1e-6] + [1e-9 * self.m0] * len(VOLATILE)
        t, y = 0.0, self.initial_state_vector()
        dt = 1e-6 * max(self.d0 / 1e-6, 1e-4)   # every time scale shrinks with the diameter in free-molecule flow
        rows = [self.diagnostics(t, y)[0]]
        end = None

        def ground(tt, yy):
            return earth.ecef_to_geodetic(yy[:3])[0] - s.ground_altitude
        ground.terminal, ground.direction = True, -1

        def gone(tt, yy):
            return float(np.maximum(yy[7:], 0.0).sum()) + self.m_inert - self.VAPORISED * self.m0
        gone.terminal, gone.direction = True, -1

        while end is None:
            t_next = min(t + dt, self.T_MAX)
            sol = solve_ivp(self.rhs, (t, t_next), y, method="LSODA", rtol=1e-8, atol=atol,
                            events=(ground, gone))
            if not sol.success:
                raise RuntimeError("integration failed at t = {:.4g} s: {}".format(t, sol.message))
            self.nfev += int(sol.nfev)
            t, y = float(sol.t[-1]), sol.y[:, -1].copy()
            row, (r, v, T, m_vol, m, d, fs, V) = self.diagnostics(t, y)
            if sol.t_events[0].size:
                end = "ground"
            elif sol.t_events[1].size:
                end = "vaporised"
            elif t >= self.T_MAX:
                end = "t_max"
            else:
                vt = self.terminal_velocity(fs, r, T, d, m)
                T_eq = self.equilibrium_temperature(fs, V, m_vol, d)
                row.update(v_terminal_ms=vt, T_eq_K=T_eq)
                if abs(V / vt - 1.0) < self.STEADY_SPEED_TOL and abs(T - T_eq) < self.STEADY_T_TOL:
                    end = "steady"
            rows.append(row)
            dt = min(dt * 1.25, 2.0)
        self.end_reason = end
        return rows


@dataclass
class _Aero:
    """The fields of trajectory.AeroState that PhysicsHeating.evaluate reads."""
    freestream: object
    V: float
    kn: float
    ma: float


def final_state(rows, sim, release, d0_um, h0_km, We0):
    last = rows[-1]
    out = dict(d0_um=d0_um, h0_km=h0_km, V0_ms=release.velocity, gamma0_deg=math.degrees(release.flight_path),
               We_release=We0, skipped=False, end_reason=sim.end_reason,
               t_final_s=last["t_s"], h_final_km=last["h_km"], V_final_ms=last["V_ms"],
               v_terminal_final_ms=last.get("v_terminal_ms", float("nan")),
               T_final_K=last["T_K"], T_eq_final_K=last.get("T_eq_K", float("nan")),
               mass0_kg=sim.m0, mass_final_kg=last["mass_kg"],
               mass_lost_kg=sim.m0 - last["mass_kg"], mass_lost_frac=1.0 - last["mass_kg"] / sim.m0,
               d_final_um=last["d_um"],
               T_peak_K=max(r["T_K"] for r in rows), We_max=max(r["We"] for r in rows),
               peak_decel_g=max(r["decel_g"] for r in rows),
               boiled=any(r["vapour_over_ambient"] >= 1.0 for r in rows), nfev=sim.nfev)
    for e in VOLATILE:
        m0e = sim.m0 * e.mass_fraction
        out["lost_{}_kg".format(e.name)] = m0e - last["m_{}_kg".format(e.name)]
        out["lost_{}_frac_of_element".format(e.name)] = 1.0 - last["m_{}_kg".format(e.name)] / m0e
    return out


# ---------------------------------------------------------------- vertical wind and the settling phase

class MeanVerticalWind:
    """Mean (residual-circulation) vertical wind [m/s, + up] in the middle atmosphere, as bounding cases:

      above 100 km      0
      90-100 km         w_meso tapered linearly to 0
      70-90 km          w_meso (mesospheric upwelling/downwelling, a few cm/s; Vincent 2019, Chau et al. 2019)
      50-70 km          linear from w_strat at 50 km to w_meso at 70 km (the mesospheric flow weakens downward)
      15-50 km          w_strat (Brewer-Dobson circulation, about 0.1-0.5 mm/s; ACP 26, 3743, 2026; CMIP6 study)
      13-15 km          w_strat tapered to 0
      below 13 km       0 (troposphere; particles there are removed mainly by precipitation, not settling)

    The layer shape is an assumption; see citations/vertical_winds_middle_atmosphere.md and VERTICAL_WINDS.md."""

    H_TROPOPAUSE = 15e3

    def __init__(self, w_meso, w_strat=0.0):
        self.w0, self.w_strat = float(w_meso), float(w_strat)

    def w(self, h):
        km = h / 1e3
        if km >= 100.0 or km <= 13.0:
            return 0.0
        if km > 90.0:
            return self.w0 * (100.0 - km) / 10.0
        if km >= 70.0:
            return self.w0
        if km > 50.0:
            return self.w_strat + (self.w0 - self.w_strat) * (km - 50.0) / 20.0
        if km >= 15.0:
            return self.w_strat
        return self.w_strat * (km - 13.0) / 2.0

    def wind_enu(self, h):
        return np.array([0.0, 0.0, self.w(h)])


def settle(sim, rows, wind, dh=50.0, h_top=100e3):
    """From the final state to the ground. Quasi-steady: the droplet sinks at its terminal velocity (Singh or blend
    drag, at the local air temperature, with its final mass) minus the local updraft, dh/dt = -(v_t(h) - w(h)),
    exact while the response time v_t/g is short against the time to fall a scale height. Evaporation after the
    final state is neglected (the droplet is at the air temperature). If v_t falls below w the droplet is held at
    the altitude where they balance ('suspended'). Returns a dict of settling results."""
    last = rows[-1]
    lat, lon = math.radians(last["lat_deg"]), math.radians(last["lon_deg"])
    m, d, h1 = last["mass_kg"], last["d_um"] * 1e-6, last["h_km"] * 1e3
    atm = sim.atmosphere

    def net(h):        # sinking speed over the ground [m/s], positive downward
        fs = atm.state(0.0, h, lat, lon)
        r = earth.geodetic_to_ecef(h, lat, lon)
        return sim.terminal_velocity(fs, r, fs.T, d, m) - wind.w(h)

    hs = np.arange(0.0, h_top + dh, dh)
    hs = np.unique(np.concatenate((hs, [h1])))
    f = np.array([net(h) for h in hs])
    k1 = int(np.searchsorted(hs, h1))
    out = dict(settle_start_km=h1 / 1e3, v_sink_start_ms=float(f[k1]),
               updraft_at_start_ms=wind.w(h1), min_net_sink_ms=float(f[:k1 + 1].min()))
    if f[k1] <= 0.0:                                   # updraft exceeds the terminal velocity: carried upward
        above = np.nonzero(f[k1:] > 0.0)[0]
        if above.size == 0:
            return dict(out, end_settle="lofted above {:.0f} km".format(h_top / 1e3), h_suspended_km=float("nan"),
                        t_ground_s=float("inf"), t_tropopause_s=float("inf"))
        j = k1 + int(above[0])
        h_eq = hs[j - 1] + (hs[j] - hs[j - 1]) * (0.0 - f[j - 1]) / (f[j] - f[j - 1])
        return dict(out, end_settle="suspended", h_suspended_km=h_eq / 1e3, t_ground_s=float("inf"),
                    t_tropopause_s=float("inf"))
    below = np.nonzero(f[:k1] <= 0.0)[0]
    if below.size:                                     # sinks until the updraft balances it, then hovers
        j = int(below[-1])
        h_eq = hs[j] + (hs[j + 1] - hs[j]) * (0.0 - f[j]) / (f[j + 1] - f[j])
        return dict(out, end_settle="suspended", h_suspended_km=h_eq / 1e3, t_ground_s=float("inf"),
                    t_tropopause_s=float("inf") if h_eq > MeanVerticalWind.H_TROPOPAUSE else float("nan"))
    t_fall = float(np.trapezoid(1.0 / f[:k1 + 1], hs[:k1 + 1]))
    kt = int(np.searchsorted(hs, MeanVerticalWind.H_TROPOPAUSE))
    t_trop = last["t_s"] + (float(np.trapezoid(1.0 / f[kt:k1 + 1], hs[kt:k1 + 1])) if k1 > kt else 0.0)
    return dict(out, end_settle="ground", h_suspended_km=float("nan"), t_fall_s=t_fall,
                t_ground_s=last["t_s"] + t_fall, t_tropopause_s=t_trop)


def run_case(d0_um, h0_km, parent=None, w_meso=None, w_strat=0.0, **kw):
    """Weber pre-check, then the flight. With w_meso [m/s] set, the mesosphere carries that mean vertical wind
    (MeanVerticalWind) and the droplet is followed past its final state to the ground (settle). Returns
    (final_state dict, history rows or None)."""
    parent = parent or sphere_sweep.load_parent(PARENT_JSON, PARENT_OBJECT)
    rel = parent_release(parent, h0_km)
    wind = MeanVerticalWind(w_meso, w_strat) if w_meso is not None else None
    atm = US76TableAtmosphere(wind=wind)
    fs = atm.state(0.0, rel.altitude, rel.lat, rel.lon)
    We0 = weber(fs.rho, rel.velocity, d0_um * 1e-6, T_RELEASE)
    if We0 >= WE_CRIT:
        return dict(d0_um=d0_um, h0_km=h0_km, V0_ms=rel.velocity, We_release=We0, skipped=True,
                    end_reason="skipped: release Weber number above {:g}".format(WE_CRIT)), None
    sim = DropletSimulator(rel, d0_um, atm, **kw)
    rows = sim.run()
    out = final_state(rows, sim, rel, d0_um, h0_km, We0)
    if wind is not None:
        out.update(w_meso_ms=wind.w0, w_strat_ms=wind.w_strat, **settle(sim, rows, wind))
    return out, rows
