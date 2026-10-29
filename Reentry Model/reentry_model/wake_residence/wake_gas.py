"""Wake gas state: the forward-moving, hot, thin gas inside a parent's plume (design spec section 4.3a).

Seen from the undisturbed air, the gas behind a re-entering parent is not still. The parent's drag D = 1/2 rho U^2 (C_D A)
(newtons) pushes it forward, in the parent's direction of travel, and does work on it that heats it at roughly the
ambient pressure, so it is also thin. A follower inside the plume feels both: the relative speed is its own velocity minus
the wake gas velocity (a tailwind), and the density is lower, so its drag is lower. `gas_state` gives the gas velocity,
temperature, density and pressure at one point of the plume; `Plume.wake_state` (overlap.py) finds that point from the
follower's position and time. Nothing here takes or knows an altitude: the ambient pressure, temperature and density are
numbers passed in, and nothing is imported from the trajectory, droplet, body, aero or atmosphere modules (spec section
4.0; a test inspects the imports).

Model (first version; every choice below is an assumption, not a measurement)
  Per unit length of the parent's track, the drag D deposits forward momentum D / U (kg m/s per m of track, so a mass
  flux of D / U in kg/s through the plume cross-section) and the drag work D (J per m) as enthalpy, at ambient pressure.

  Momentum.  The forward momentum density (momentum per unit time per unit area, kg/m2/s) is a Gaussian about the emission
  point of the track,  m(r) = (D / U) / (pi r_b^2) exp(-(r / r_b)^2),  whose area integral over the plane is exactly
  D / U. Its width is the plume radius r_b, so the profile is the spec's "smooth radial profile scaled so its integral is
  the core-mean value times pi r_b^2". The gas speed is u = m / rho_used, capped at the parent's speed U (gas pushed by a
  body cannot outrun it); rho_used is the hot gas density when temperature is on ("both") and the ambient density when only
  the velocity is on ("velocity"). Where the speed cap does not bind, rho u = m, so the momentum integral is exactly D / U
  in either mode (a test integrates it); where it binds the gas carries less momentum than D / U.

  Energy.  The turbulent core is too small to hold the energy at ambient pressure close behind a fast parent (the required
  area is about 0.2 M^2 C_D A), so the energy is spread over R_E = max(r_b, k_E sqrt(D / p)): the larger of the core and
  the line-blast radius (ReVelle), with k_E a named parameter (default 1) and p the pressure at emission (the radius belongs to the emission point; the
  temperature, density and the energy-density conversion use the local ambient air passed in). The excess
  enthalpy fraction q = 1 - T_inf / T is a Gaussian of width R_E,  q = q_c exp(-(r / R_E)^2),  with
  q_c = D R / (p c_p pi R_E^2). At fixed pressure the enthalpy per unit volume above ambient is
  rho c_p (T - T_inf) = (p c_p / R) q exactly (because rho T = p / R and T - T_inf = T q), so the area integral of q,
  q_c pi R_E^2 = D R / (p c_p), is the statement that the drag work D per unit length is all deposited (a test integrates
  it). T = T_inf / (1 - q) and rho = rho_inf (1 - q) = rho_inf T_inf / T: constant pressure, using the local ambient's
  own gas constant p_inf / (rho_inf T_inf) (the spec's rho_w = p_inf / (R T_w) with that R), so the density is
  continuous with the ambient density at the cutoff even when the table the ambient came from is not exactly
  consistent with R_AIR (a US76 table drifts by a few percent at high altitude). The centreline value q_c is capped at 1 - T_inf / T_st, where
  T_st = T_inf + U^2 / (2 c_p) is the stagnation temperature of the parent's flow, so the gas never exceeds what
  stopping the flow could give (computed from the emission speed U); where the cap binds the energy integral falls short of D R / (p c_p), because the cap
  limits the temperature. This well-mixed mean understates the centreline temperature of a real wake (assumed).

  Smoothness.  Both profiles are Gaussian in the distance from the emission point, with no step at r_b, because the
  follower's adaptive integrator stalls on a step at the plume edge. The "inside the plume" test of spec section 4.4 (the
  margin, r <= r_b) is unchanged and independent of these profiles. Beyond `cutoff` times max(r_b, R_E) the state is
  exactly ambient (at 6 widths a Gaussian is 2e-16 of its peak), so far from the plume nothing depends on the profile.

Modes (`wake_flow`): "off" ambient everywhere (the original definition); "velocity" gas velocity only, with ambient
temperature and density; "density" hot, thin gas but no velocity; "both" both effects. The pressure is the ambient
pressure in every mode.

Units: SI throughout. R_AIR and CP_AIR are for dry air, 287.05 and 1004.5 J/kg/K, and a calorically perfect gas is
assumed (a real hot wake at 5,000 K or more has a larger specific heat, so this overstates the temperature).
"""
from dataclasses import dataclass
import math

import numpy as np

from reentry_model.wake_residence.laws import EmissionState

CP_AIR = 1004.5                                   # J/kg/K: specific heat of dry air at constant pressure
R_AIR = 287.05                                    # J/kg/K: specific gas constant of dry air
WAKE_FLOW_MODES = ("off", "velocity", "density", "both")


@dataclass(frozen=True)
class Ambient:
    """The undisturbed air at the follower: pressure p [Pa], temperature T [K], density rho [kg/m3]."""
    p: float
    T: float
    rho: float


@dataclass(frozen=True)
class GasState:
    """The gas the follower sees: density rho [kg/m3], temperature T [K], pressure p [Pa] and the gas velocity u (3,) [m/s]
    relative to the undisturbed air (zero where the gas is at rest)."""
    rho: float
    T: float
    p: float
    u: np.ndarray


def ambient_state(ambient: Ambient) -> GasState:
    return GasState(ambient.rho, ambient.T, ambient.p, np.zeros(3))


def check_mode(mode: str) -> None:
    if mode not in WAKE_FLOW_MODES:
        raise ValueError(f"wake_flow must be one of {WAKE_FLOW_MODES}, got {mode!r}")


def gas_state(r_axis: float, r_b: float, emission: EmissionState, direction, ambient: Ambient, mode: str,
              k_E: float = 1.0, cutoff: float = 6.0) -> GasState:
    """Wake gas state at distance `r_axis` [m] from the emission point of the track.

    r_b        plume radius there [m] (the width of the momentum profile)
    emission   the parent's flow state at emission (its drag D, speed U and pressure p set the profiles)
    direction  the parent's velocity direction at emission, (3,), normalised here
    ambient    the undisturbed air at the follower (pressure, temperature, density)
    mode       "off", "velocity", "density" or "both" (see the module docstring)
    k_E        scale on the line-blast radius sqrt(D / p) in the energy radius R_E = max(r_b, k_E sqrt(D / p))
    cutoff     the state is exactly ambient beyond cutoff * max(r_b, R_E)
    """
    check_mode(mode)
    if mode == "off":
        return ambient_state(ambient)
    if not (math.isfinite(r_axis) and r_axis >= 0.0):
        raise ValueError(f"r_axis must be a finite, non-negative distance in metres, got {r_axis!r}")
    if not (math.isfinite(r_b) and r_b > 0.0):
        raise ValueError(f"r_b must be a finite, positive plume radius in metres, got {r_b!r}")
    if not (math.isfinite(emission.U) and emission.U > 0.0 and math.isfinite(emission.p) and emission.p > 0.0):
        raise ValueError(f"the emission state needs a positive speed and pressure, got U = {emission.U!r}, "
                         f"p = {emission.p!r}")
    d = np.asarray(direction, dtype=float).reshape(3)
    norm = float(np.linalg.norm(d))
    if not (math.isfinite(norm) and norm > 0.0):
        raise ValueError(f"direction must be a non-zero finite vector, got {direction!r}")
    d = d / norm

    D = 0.5 * emission.rho * emission.U**2 * emission.cda                       # drag at emission [N]
    R_E = max(r_b, k_E * math.sqrt(D / emission.p))                             # energy radius [m]
    if r_axis > cutoff * max(r_b, R_E):
        return ambient_state(ambient)

    T_st = ambient.T + emission.U**2 / (2.0 * CP_AIR)                           # stagnation-temperature cap [K]
    q_c = min(D * R_AIR / (ambient.p * CP_AIR * math.pi * R_E**2), 1.0 - ambient.T / T_st)
    q = q_c * math.exp(-(r_axis / R_E)**2) if mode in ("density", "both") else 0.0
    T = ambient.T / (1.0 - q)
    rho = ambient.rho * (1.0 - q)                                                # constant pressure: rho T = rho_inf T_inf
    if mode == "density":
        return GasState(rho, T, ambient.p, np.zeros(3))

    m = (D / emission.U) / (math.pi * r_b**2) * math.exp(-(r_axis / r_b)**2)   # forward momentum density [kg/m2/s]
    u = min(m / (rho if mode == "both" else ambient.rho), emission.U)
    return GasState(rho, T, ambient.p, u * d)
