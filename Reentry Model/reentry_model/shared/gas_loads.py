"""The interim gas-loads provider: pressure, wall shear and similarity heat flux round the whole body (0-180 degrees
from the stagnation point) for the surface flow and the heating, solved once per macro step and shared.

Pressure, by the perfect-gas free-stream Mach number M (`daniels_schwartzentruber.condition(...).mach_inf`): at and
above M = 5 the model's windward `surface_flow.wall_pressure` (modified Newtonian plus Prandtl-Meyer) to 90 degrees; past
90 degrees the same Prandtl-Meyer expansion continued without the free-stream floor, until it falls to Daniels and
Schwartzentruber's (2026) lee pressure, then theirs (`joined_pressure`). That joined law has a step at 90 degrees below
about Mach 3.4, where the Newtonian cosine-squared law and the Prandtl-Meyer expansion no longer agree. At and below
M = 3.5 the pressure is Daniels and Schwartzentruber's own, p_inf + c_p q_inf, rescaled so that its stagnation value is
the unreacted stagnation pressure (their c_p at the stagnation point is only approximately the frozen-gas value, and
the rescaling keeps the stagnation pressure the same quantity at every Mach number); between M = 3.5 and 5 the two
pressures are blended with a weight falling linearly in M from 1 to 0 (a choice, not a fitted law). The joined
expansion starts from the UNREACTED (chemically frozen) stagnation pressure. Shear and heat flux: Lees local similarity (`similarity.march`) on that pressure. Below the shear
switch Kn_body = 0.05 the continuum similarity shear is returned (the caller applies slip and blend); at and above it
Daniels and Schwartzentruber's shear, already rarefied. Ported from analysis/gas_loads/reference/similarity_lee.py
(`joined_pressure`) and similarity_heat.py (`main`). This module must not be imported by surface_flow."""
import math
from dataclasses import dataclass

import numpy as np

from . import daniels_schwartzentruber as ds
from . import surface_flow
from .gas import UnreactedAir
from .similarity import FrozenAirProperties, SimilarityError, SimilaritySolver, march

KN_SHEAR_SWITCH = 0.05
BIN_DEG = 2.0
R_GAS = 287.05
_CACHE_SIZE = 8
# The joined pressure is clamped to P_FLOOR times the freestream pressure before Cantera's isentropic expansion, which
# refuses a pressure at or below zero. Only Daniels and Schwartzentruber's lee pressure p_inf + c_p q_inf can be that low
# (the windward law is floored at p_inf and the Prandtl-Meyer pressure is positive by construction), where their c_p falls
# below -p_inf/q_inf. The floor is more than five orders of magnitude below the lowest lee pressure their model gives
# over 1e-12 to 1.26 kg/m3, 100 to 8000 m/s, 5 to 300 mm, 190 to 1000 K and walls at 300 to 2500 K (7.4e-15 of the
# freestream, at 1e-12 kg/m3), so it changes no result. Cantera cannot expand to the floor either (the isentropic
# temperature falls below zero), so a state that reaches it still fails its step, on Cantera's temperature instead.
P_FLOOR = 1e-20
# Daniels and Schwartzentruber's pressure is used alone at and below DS_PRESSURE_MACH_LOW and not at all at and above
# DS_PRESSURE_MACH_HIGH (perfect-gas free-stream Mach number); linear weight in between.
DS_PRESSURE_MACH_LOW = 3.5
DS_PRESSURE_MACH_HIGH = 5.0


def ds_pressure_weight(mach):
    """Weight of the anchored Daniels and Schwartzentruber pressure: 1 at and below DS_PRESSURE_MACH_LOW, 0 at and
    above DS_PRESSURE_MACH_HIGH, linear in Mach between."""
    if mach <= DS_PRESSURE_MACH_LOW:
        return 1.0
    if mach >= DS_PRESSURE_MACH_HIGH:
        return 0.0
    return (DS_PRESSURE_MACH_HIGH - mach) / (DS_PRESSURE_MACH_HIGH - DS_PRESSURE_MACH_LOW)


def anchored_ds_pressure(cp_ds, p_stag, p_inf, q_inf):
    """Daniels and Schwartzentruber's pressure p_inf + c_p q_inf with c_p rescaled so the stagnation-point value is
    p_stag: p_inf + c_p q_inf (p_stag - p_inf) / (c_p[0] q_inf). Raises ValueError if c_p at the stagnation point
    (cp_ds[0]) is not positive."""
    cp_ds = np.asarray(cp_ds, dtype=float)
    if not cp_ds[0] > 0.0:
        raise ValueError(f"Daniels and Schwartzentruber's stagnation-point pressure coefficient is {cp_ds[0]:.4g}; "
                         "it must be positive to anchor the pressure")
    return p_inf + cp_ds * q_inf * (p_stag - p_inf) / (cp_ds[0] * q_inf)


def joined_pressure(theta_deg, p_stag, p_inf, q_inf, gamma_pm, p_ds):
    """Wall pressure on the bins theta_deg [deg]: the model's wall_pressure to 90 degrees, the uncapped Prandtl-Meyer
    expansion beyond, and from the first angle above 90 degrees where that is at or below p_ds, p_ds itself.
    Returns (p, meet angle in degrees, nan if they never meet)."""
    th_deg = np.asarray(theta_deg, dtype=float)
    th = np.radians(th_deg)
    p_wind, phi_star = surface_flow.wall_pressure(np.minimum(th, 0.5 * math.pi), p_stag, p_inf, q_inf, gamma_pm)
    turn = np.maximum(th - phi_star, 0.0)
    mach = surface_flow.mach_from_turn(turn, gamma_pm)
    p_pm = p_stag * (1.0 + 0.5 * (gamma_pm - 1.0) * mach**2) ** (-gamma_pm / (gamma_pm - 1.0))
    p = np.where(th_deg <= 90, p_wind, p_pm)
    lee = th_deg > 90
    hit = lee & (p_pm <= p_ds)
    meet = int(np.argmax(hit)) if np.any(hit) else len(th_deg)
    p[meet:] = np.asarray(p_ds)[meet:]
    return p, float(th_deg[meet]) if meet < len(th_deg) else float("nan")


@dataclass
class MeridianLoads:
    theta_deg: np.ndarray
    p_w: np.ndarray          # Pa
    dp_ds: np.ndarray        # Pa/m, pressure gradient along the surface
    tau: np.ndarray          # Pa; continuum below the switch, Daniels and Schwartzentruber's at and above it
    q_c: np.ndarray          # W/m2, continuum similarity heat flux
    q_c0: float              # W/m2, stagnation-point similarity flux
    H_e: float               # J/kg, total enthalpy of the unreacted post-shock gas
    h_w: float               # J/kg, enthalpy of the wall-temperature gas
    T_wall: float
    kn_body: float
    shear_source: str        # "similarity" or "daniels-schwartzentruber"
    meet_deg: float
    separation_deg: float
    p_stag: float
    p_inf: float
    q_inf: float
    ds_pressure_weight: float = 0.0   # weight of the anchored Daniels and Schwartzentruber pressure in p_w (0 to 1)

    def at(self, theta_rad, field):
        """Linear interpolation of `field` (an array attribute name) at polar angles theta_rad, clipped to [0, pi]."""
        th = np.degrees(np.clip(np.asarray(theta_rad, dtype=float), 0.0, math.pi))
        return np.interp(th, self.theta_deg, getattr(self, field))


class GasLoads:
    def __init__(self, gamma_pm=surface_flow.GAMMA_PM, kn_switch=KN_SHEAR_SWITCH, bin_deg=BIN_DEG):
        self.gamma_pm, self.kn_switch, self.bin_deg = gamma_pm, kn_switch, bin_deg
        self._cache = {}
        self._air = self._props = self._solver = None

    def __getstate__(self):
        d = dict(self.__dict__)
        d["_cache"] = {}
        d["_air"] = d["_props"] = d["_solver"] = None      # rebuilt on first use (Cantera objects, warm start)
        return d

    def _build(self):
        if self._air is None:
            self._air = UnreactedAir()
            self._props = FrozenAirProperties(self._air)
            self._solver = SimilaritySolver(self._props)

    def wall_enthalpy(self, T):
        """Frozen-air specific enthalpy (J/kg) at wall temperature(s) T (K), the same spline `meridian` uses for h_w;
        vectorised (a scalar in gives a 0-d array)."""
        self._build()
        return np.asarray(self._props._h(np.asarray(T, dtype=float)))

    def meridian(self, state, radius, T_wall):
        fs = state.freestream
        key = (fs.rho, fs.T, state.V, radius, T_wall)
        hit = self._cache.get(key)
        if hit is None:
            hit = self._solve(state, radius, T_wall)
            if len(self._cache) >= _CACHE_SIZE:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = hit
        return hit

    def _solve(self, state, radius, T_wall):
        fs = state.freestream
        rho, T_inf, V, p_inf = fs.rho, fs.T, state.V, fs.p
        th = np.arange(0.0, 180.0 + 0.5 * self.bin_deg, self.bin_deg)
        n = len(th)
        q_inf = 0.5 * rho * V * V
        if not (rho > 0.0 and V > 0.0):
            z = np.zeros(n)
            return MeridianLoads(th, np.full(n, p_inf), z.copy(), z.copy(), z.copy(), 0.0, 0.0, 0.0, T_wall,
                                 state.kn, "similarity", float("nan"), float("nan"), p_inf, p_inf, 0.0, 0.0)
        self._build()
        air, props = self._air, self._props
        stag = air.stagnation(rho, T_inf, V)
        H_e = stag.h
        cond = ds.condition(rho, T_inf, V, 2.0 * radius, T_wall)
        cp_ds, cf_ds = ds.coefficients(np.radians(th), cond)
        p_ds = p_inf + cp_ds * q_inf
        w = ds_pressure_weight(cond.mach_inf)
        if w == 1.0:
            p_w, meet = anchored_ds_pressure(cp_ds, stag.p, p_inf, q_inf), float("nan")
        else:
            joined, meet = joined_pressure(th, stag.p, p_inf, q_inf, self.gamma_pm, p_ds)
            if w == 0.0:
                p_w = joined
            else:
                p_w = w * anchored_ds_pressure(cp_ds, stag.p, p_inf, q_inf) + (1.0 - w) * joined
        p_w = np.maximum(p_w, P_FLOOR * p_inf)            # never a non-positive pressure into Cantera (P_FLOOR)
        u_e = np.array([math.sqrt(max(0.0, 2.0 * (H_e - air.expand(stag, float(p)).h))) for p in p_w])
        u_e[0] = 0.0
        mu_w = float(props.mu_of_T(T_wall))
        rmw = p_w / (R_GAS * T_wall) * mu_w
        try:
            sol = march(th, radius, p_w, u_e, H_e, T_wall, rmw, self._solver)
        except SimilarityError as err:
            raise SimilarityError(f"similarity stagnation-point solution failed at rho={rho:.4g} T={T_inf:.4g} "
                                  f"V={V:.4g} Kn_body={state.kn:.4g}: {err}") from err
        if state.kn < self.kn_switch:
            tau, source = sol.tau_c, "similarity"
        else:
            tau, source = cf_ds * q_inf, "daniels-schwartzentruber"
        dp_ds = np.gradient(p_w, radius * np.radians(th))
        return MeridianLoads(th, p_w, dp_ds, tau, sol.q_c, float(sol.q_c[0]), H_e, props.h_of_T(T_wall), T_wall,
                             state.kn, source, meet, sol.separation_deg, stag.p, p_inf, q_inf, w)
