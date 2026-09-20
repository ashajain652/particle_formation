"""Convective heating on the sphere's surface patches (spec section 6).

Two models behind one call, `evaluate(state, theta, T_wall, radius) -> HeatingResult`, where `state` is the
trajectory's AeroState (freestream rho/T/p, airspeed V, Kn, Ma), `theta` the angle of every patch's outward normal
to the direction of motion (0 at the stagnation point), `T_wall` every patch's temperature and `radius` the nose
radius. `q_conv` is W/m2 into the body per patch.

SesamEquivalentHeating -- VERIFICATION DEVICE, NOT A PHYSICAL MODEL. It reproduces SESAM's lumped-body heat input as
measured on the two no-melt US76 references (2026-09-18, superseding the (1 - f) q_DKR + f q_FM decomposition of
facts note s.4): q = 0.27471 x q_DKR x F_h(Kn) x max(0, 1 - c_p,air (T_lumped - T_inf) / (V^2/2)) for Ma >= 1, with
SESAM's own heat bridging F_h (aero.SesamHeatTable: 1.005 in the continuum, 0.14 at Kn 1, 0.06 at Kn 40) and its
hot-wall factor (c_p 1004.5 J/kg/K: rms 0.008, max 0.015 of the ratio over 189 rows up to 2450 K; the factor is
1/2 with no hot-wall term below Ma 1), the surface-average factor 0.27471 applied uniformly to EVERY patch, front
and back. Real heating is concentrated on the windward face (Lees: 0.196 of the stagnation flux averaged over the
sphere, zero leeward). The uniform distribution exists so that our conduction, time stepping, material curves and
coupling can be checked against SESAM's lumped temperature and heat totals with no distribution question in between.

PhysicsHeating -- the model proper: continuum stagnation flux from Fay-Riddell (equilibrium air, `gas.py`),
Sutton-Graves or DKR (with the hot-wall factor (h_s - h_w)/(h_s - h_w,300K) where the correlation is cold-wall),
free-molecular stagnation flux A_cq rho V (h_s - h_w) (Matting Eq. 1), Matting's bridging (facts note s.14) or
SESAM's f(Kn), and the distribution (1 - w) Lees(theta, Ma) + w cos(theta) with the free-molecular weight
w = 1 - q_stag/q_c (Matting) or f(Kn) (SESAM table); per patch the flux is scaled by its own hot-wall factor
relative to the stagnation patch. Fixed attitude: theta is whatever the caller passes.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.special import gamma, gammainc

from . import aero, gas
from .constants import GAMMA_AIR

SESAM_SHAPE_FACTOR = 0.27471        # ATDB_SPHERE continuum heat factor: surface average / stagnation flux
CP_AIR_SESAM = 1004.5               # J/(kg K): SESAM's hot-wall enthalpy h_w - h_inf = c_p (T_w - T_inf) (best fit 985; 1004.5 within 1.5 %)
T_WALL_COLD = 300.0                 # K, wall temperature the cold-wall correlations are referred to
HEATING_NAMES = ("sesam", "physics")
STAGNATION_NAMES = ("fay-riddell", "sutton-graves", "dkr")
BRIDGING_HEAT_NAMES = ("matting", "sesam-table")


@dataclass
class HeatingResult:
    q_conv: np.ndarray        # W/m2 per patch, positive into the body
    q_stag: float             # bridged stagnation-point flux
    q_stag_c: float           # continuum stagnation value
    q_stag_fm: float          # free-molecular stagnation value
    blend: float              # free-molecular weight of the distribution: w (physics) or f(Kn) (sesam)

    def total(self, areas):
        return float((self.q_conv * areas).sum())


def lees_shape(theta, mach):
    """Lees (1956) laminar sphere distribution q(theta)/q_stag for theta <= pi/2 (zero leeward), finite Mach through
    eps = 1/(gamma Ma^2); Ma is clamped to >= 1. Integrates over the windward hemisphere to 0.196 (Ma -> inf),
    0.200 (Ma 10), 0.209 (Ma 5) of q_stag x sphere area."""
    theta = np.asarray(theta, dtype=float)
    eps = 0.0 if not np.isfinite(mach) else 1.0 / (GAMMA_AIR * max(mach, 1.0) ** 2)
    D = (1.0 - eps) * (theta ** 2 - 0.5 * theta * np.sin(4.0 * theta) + 0.125 * (1.0 - np.cos(4.0 * theta))) \
        + 4.0 * eps * (theta ** 2 - theta * np.sin(2.0 * theta) + 0.5 * (1.0 - np.cos(2.0 * theta)))
    small = theta < 1e-4                                       # limit theta -> 0 is exactly 1
    with np.errstate(divide="ignore", invalid="ignore"):
        q = 2.0 * theta * np.sin(theta) * ((1.0 - eps) * np.cos(theta) ** 2 + eps) / np.sqrt(D)
    q = np.where(small, 1.0, q)
    return np.where(theta <= 0.5 * math.pi, q, 0.0)


def cosine_shape(theta):
    """Free-molecular distribution cos(theta) on the windward hemisphere (integrates to 0.25 of q_stag x area)."""
    theta = np.asarray(theta, dtype=float)
    return np.where(theta <= 0.5 * math.pi, np.cos(theta), 0.0)


def matting_bridge(q_fm, q_c, n=1.0):
    """Matting (1971) Eq. 18: q = q_c P(n, [Gamma(n+1) q_fm/q_c]^(1/n)) with P the regularised lower incomplete gamma
    function; n = 1 gives q_c [1 - exp(-q_fm/q_c)]. Limits: q -> q_fm when q_fm << q_c, q -> q_c when q_fm >> q_c."""
    if q_c <= 0.0:
        return 0.0
    if n == 1.0:
        return q_c * (1.0 - math.exp(-q_fm / q_c))
    return q_c * float(gammainc(n, (gamma(n + 1.0) * q_fm / q_c) ** (1.0 / n)))


class SesamEquivalentHeating:
    """SESAM's heating on every patch: q = 0.27471 x q_DKR x F_h(Kn) x max(0, 1 - c_p (T_lumped - T_inf) / (V^2/2))
    for Ma >= 1, and 0.5 x 0.27471 x q_DKR x F_h(Kn) with no hot-wall term below Ma 1 (both measured, module docstring).

    NOT PHYSICAL -- a verification device: the uniform factor reproduces SESAM's surface-averaged lumped-body input so
    that the conduction/coupling can be compared with SESAM's lumped temperature. No distribution over theta; the
    hot-wall factor uses the body's lumped (energy-equivalent) temperature `T_mean` as SESAM does, falling back to
    the mean of `T_wall`. HeatingResult: q_stag is the effective stagnation flux q_DKR x F_h x hot-wall, q_stag_c the
    DKR value, q_stag_fm the textbook 1/2 rho V^3 for reference only, blend = F_h(Kn)."""

    name = "sesam"

    def __init__(self, heat_bridging=None, shape_factor=SESAM_SHAPE_FACTOR):
        self.heat_bridging = heat_bridging or aero.SesamHeatTable()
        self.shape_factor = shape_factor

    def evaluate(self, state, theta, T_wall, radius, T_mean=None):
        rho, V = state.freestream.rho, state.V
        n_patches = len(np.asarray(theta))
        if rho <= 0.0 or V <= 0.0:
            return HeatingResult(np.zeros(n_patches), 0.0, 0.0, 0.0, 1.0)
        q_c = gas.dkr(rho, V, radius)
        F_h = self.heat_bridging(state.kn) if np.isfinite(state.kn) else float(self.heat_bridging.F[-1])
        T_lumped = float(np.mean(T_wall)) if T_mean is None else float(T_mean)
        if state.ma >= 1.0:
            hot_wall = max(0.0, 1.0 - CP_AIR_SESAM * (T_lumped - state.freestream.T) / (0.5 * V * V))
        else:
            hot_wall = 0.5                       # SESAM below Ma 1: half the continuum value, no hot-wall factor (measured)
        q_stag = q_c * F_h * hot_wall
        # verification only: SESAM's surface-average factor on EVERY patch, windward and leeward alike (not physical)
        return HeatingResult(np.full(n_patches, self.shape_factor * q_stag), q_stag, q_c, 0.5 * rho * V ** 3, F_h)


class PhysicsHeating:
    """Stagnation correlation + Matting/SESAM bridging + Lees/cos(theta) distribution (module docstring)."""

    name = "physics"

    def __init__(self, stagnation="fay-riddell", bridging="matting", matting_n=1.0, accommodation=0.8,
                 catalycity=1.0, air=None):
        if stagnation not in STAGNATION_NAMES:
            raise ValueError("stagnation must be one of {}, got {!r}".format(STAGNATION_NAMES, stagnation))
        if bridging not in BRIDGING_HEAT_NAMES:
            raise ValueError("bridging must be one of {}, got {!r}".format(BRIDGING_HEAT_NAMES, bridging))
        self.stagnation, self.bridging, self.matting_n = stagnation, bridging, matting_n
        self.accommodation, self.catalycity = accommodation, catalycity
        self.air = air or gas.EquilibriumAir()
        self.sesam_table = aero.SesamTable()

    def stagnation_fluxes(self, state, T_wall_stag, radius):
        """(q_c, q_fm, h_s) at the stagnation patch for the freestream `state` and the stagnation wall temperature."""
        fs = state.freestream
        rho, V = fs.rho, state.V
        h_s = float(self.air.air_enthalpy(fs.T)) + 0.5 * V * V
        h_w = float(self.air.air_enthalpy(T_wall_stag))
        hot_wall = (h_s - h_w) / (h_s - float(self.air.air_enthalpy(T_WALL_COLD)))
        if self.stagnation == "fay-riddell":
            stag = self.air.stagnation(rho, fs.T, V)
            wall = self.air.wall(T_wall_stag, stag.p)
            q_c, h_s = gas.fay_riddell(stag, wall, fs.p, radius, self.catalycity), stag.h
        elif self.stagnation == "sutton-graves":
            q_c = gas.sutton_graves(rho, V, radius) * hot_wall
        else:
            q_c = gas.dkr(rho, V, radius) * hot_wall
        q_fm = self.accommodation * rho * V * (h_s - h_w)
        return q_c, q_fm, h_s

    def evaluate(self, state, theta, T_wall, radius, T_mean=None):
        theta, T_wall = np.asarray(theta, dtype=float), np.asarray(T_wall, dtype=float)
        rho, V = state.freestream.rho, state.V
        if rho <= 0.0 or V <= 0.0:
            return HeatingResult(np.zeros(theta.size), 0.0, 0.0, 0.0, 1.0)
        i_stag = int(np.argmin(theta))
        q_c, q_fm, h_s = self.stagnation_fluxes(state, float(T_wall[i_stag]), radius)
        if self.bridging == "matting":
            q_stag = matting_bridge(q_fm, q_c, self.matting_n)
            w = 1.0 - q_stag / q_c if q_c > 0.0 else 1.0
        else:
            w = self.sesam_table(state.kn) if np.isfinite(state.kn) else 1.0
            q_stag = max(0.0, (1.0 - w) * q_c + w * q_fm)
        shape = (1.0 - w) * lees_shape(theta, state.ma) + w * cosine_shape(theta)
        h_w = self.air.air_enthalpy(T_wall)
        per_patch_hot_wall = (h_s - h_w) / (h_s - h_w[i_stag])
        return HeatingResult(q_stag * shape * per_patch_hot_wall, float(q_stag), float(q_c), float(q_fm), float(w))


class TabulatedHeating:
    """Reserved (spec 6.3): q_stag and q(theta)/q_stag interpolated from an offline CFD/DSMC table over (V, rho, T_w).
    Only the interface exists; constructing it raises NotImplementedError."""

    name = "tabulated"

    def __init__(self, table_path):
        raise NotImplementedError("TabulatedHeating is an interface slot; no table format is defined yet")

    def evaluate(self, state, theta, T_wall, radius, T_mean=None):
        raise NotImplementedError


def heating_by_name(name, **options):
    if name == "sesam":
        return SesamEquivalentHeating(**options)
    if name == "physics":
        return PhysicsHeating(**options)
    raise ValueError("heating must be one of {}, got {!r}".format(HEATING_NAMES, name))
