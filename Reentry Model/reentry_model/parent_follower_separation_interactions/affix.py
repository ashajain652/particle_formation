"""Release classifier: where does a fragment let go of a hypersonic parent end up, from the size ratio and release angle.

This is the AFFIx fit of Register et al. (2020), "Interactions between asteroid fragments during atmospheric entry",
Icarus 337:113468 (Cart3D inviscid force databases with US3D viscous corrections for the near wake, checked against
free-flight tests of two 5 cm spheres in the DLR H2K tunnel). It answers one question for a small sphere released from
the surface of a larger sphere in hypersonic continuum flow: with r = D2/D1 the ratio of the follower's diameter to the
parent's, and theta the angle at the parent's centre between the direction opposite to the flight velocity and the line
to the release point (theta = 0 directly behind the parent, theta = 180 directly ahead, at the stagnation point), which
of four final zones does the follower end in? (Their Eq. 5 and Table 3; source notes in the main checkout at
citations/fragment_shock_separation.md section 00.)

    trapped    (their z_f = 3, near wake)  theta < theta_inn(r): stays in the near wake behind the parent;
    far_wake   (z_f = 2)                   theta_inn <= theta < theta_bif - s_w: stays in the wake while the streamwise
                                           gap grows;
    shock_band (z_f = 1)                   |theta - theta_bif| <= s_w: carried along on the parent's bow shock;
    escape     (z_f = 0)                   theta > theta_bif + s_w: flies independently.

with theta_inn = 122 r + 47 degrees (fitted only for r <= r* = 0.353) and theta_bif = 100 exp(-1.618 r) + 40 degrees.
s_w, the half-width of the shock band, is "about 1 to 3 degrees" in the paper; we use 2 (`AffixSettings.s_w_deg`).
For r above r* the paper gives no inner threshold: theta_inn is NaN here, `r_above_rstar` is flagged, and the zone is
decided by theta_bif alone (far_wake below the band, no trapped zone, and no time for the far wake, whose A0 is also
defined only for r <= r*).

The interaction time t_fz = D1 sqrt(rho / q_inf) t_hat (their Eq. 4) is the time over which the pair interacts, in
seconds; rho is the bulk density of both bodies, q_inf = 0.5 rho_air U^2 the free-stream dynamic pressure, and t_hat is
a dimensionless time from their Eqs. 7-9: 0 for trapped; A0(r) + 0.5 exp[-0.05 (theta - theta_inn)] with
A0 = -7.2 r^2 + 3 r + 0.1 for the far wake; 0.88 r + 0.32 + 0.01 exp[0.125 (theta - theta_bif - 10)] for escape. In
the shock band no t_hat is given and NaN is returned (the repository does not integrate shock-band cases).

Validity (their p. 9), recorded in `AffixResult.flags` and never used to refuse a computation:
  - the air is a continuum (no rarefied-gas effects): `rarefied_parent` and `rarefied_follower` are set when the mean
    free path divided by that body's diameter exceeds `kn_continuum` (0.01, as in `laws.LawSettings`);
  - the free stream is hypersonic: `not_hypersonic` when the Mach number U / sqrt(gamma R T) is below `mach_min` (5);
  - both fragments have the same bulk density: `unequal_density` when a parent mass is given and its bulk density,
    m / (pi D1^3 / 6), differs from the follower's by more than 1 percent (the value is in `parent_bulk_density`,
    NaN without a mass; our parent is a hollow shell, so for the Generic Spacecraft this is always set);
  - the fitted ratios are r = 1/16 to 1 (`R_FIT_RANGE`): `r_below_fit` below it, so a 10 mm or 100 mm sphere on a
    2.8 m parent (r = 0.0036, 0.036) is an extrapolation; `r_above_rstar` above 0.353.

Nothing here takes or knows an altitude: the flow state comes in as an `EmissionState` and the follower's bulk density
as an argument, and nothing is imported from the trajectory, droplet, body, aero or atmosphere modules (spec section
4.0; a test inspects the imports). `release_offset` and `far_wake_offset` turn a release angle into a position
relative to the parent's centroid for the follower's start.
"""
from dataclasses import dataclass
import math

import numpy as np

R_STAR = 0.353                       # largest size ratio for which the paper gives theta_inn and A0
R_FIT_RANGE = (1.0 / 16.0, 1.0)      # size ratios covered by the Cart3D databases
GAMMA_AIR = 1.4                      # ratio of specific heats, for the Mach number
R_AIR = 287.05                       # specific gas constant of air, J/(kg K)
ZONES = ("trapped", "far_wake", "shock_band", "escape")
DENSITY_TOL = 0.01                   # relative difference of bulk densities above which unequal_density is set


def theta_inn(r):
    """Inner threshold angle in degrees: below it the follower stays in the near wake. NaN for r above R_STAR."""
    return 122.0 * r + 47.0 if r <= R_STAR else math.nan


def theta_bif(r):
    """Bifurcation angle in degrees: the release angle at which the follower stops being wake-bound."""
    return 100.0 * math.exp(-1.618 * r) + 40.0


@dataclass(frozen=True)
class AffixSettings:
    """Named, overridable choices (the thresholds that set the validity flags and the shock-band half-width)."""
    s_w_deg: float = 2.0             # shock-band half-width, degrees (paper: about 1 to 3)
    kn_continuum: float = 0.01       # mean free path / diameter above which that body is rarefied
    mach_min: float = 5.0            # free-stream Mach number below which the flow is not hypersonic


@dataclass(frozen=True)
class AffixResult:
    zone: str                        # one of ZONES
    r: float                         # size ratio D2 / D1
    theta_deg: float                 # release angle as given
    theta_inn: float                 # degrees, NaN for r above R_STAR
    theta_bif: float                 # degrees
    t_hat: float                     # dimensionless interaction time (NaN in the shock band or if undefined)
    t_fz: float                      # interaction time, s (NaN where t_hat is)
    flags: dict                      # validity record, see the module docstring


def _t_hat(zone, r, theta, t_inn, t_bif):
    if zone == "trapped":
        return 0.0
    if zone == "far_wake":
        if r > R_STAR:
            return math.nan
        a0 = -7.2 * r**2 + 3.0 * r + 0.1
        return a0 + 0.5 * math.exp(-0.05 * (theta - t_inn))
    if zone == "escape":
        return 0.88 * r + 0.32 + 0.01 * math.exp(0.125 * (theta - t_bif - 10.0))
    return math.nan


def affix_regime(d_follower, d_parent, theta_deg, emission, follower_density, parent_mass=None,
                 settings=AffixSettings()):
    """Classify a release (Register et al. 2020). `emission` is the parent's `EmissionState` at the release time and
    `follower_density` the follower's bulk density in kg/m3. Raises ValueError for a non-positive diameter or density,
    or for theta outside [0, 180] degrees."""
    if not d_follower > 0.0 or not d_parent > 0.0:
        raise ValueError(f"diameters must be positive, got follower {d_follower}, parent {d_parent}")
    if not follower_density > 0.0:
        raise ValueError(f"follower_density must be positive, got {follower_density}")
    if not 0.0 <= theta_deg <= 180.0:
        raise ValueError(f"theta_deg must lie in [0, 180] degrees, got {theta_deg}")
    r = d_follower / d_parent
    t_inn, t_bif = theta_inn(r), theta_bif(r)
    if abs(theta_deg - t_bif) <= settings.s_w_deg:
        zone = "shock_band"
    elif theta_deg > t_bif:
        zone = "escape"
    elif theta_deg < t_inn:                     # False when t_inn is NaN (r above R_STAR): no trapped zone
        zone = "trapped"
    else:
        zone = "far_wake"
    t_hat = _t_hat(zone, r, theta_deg, t_inn, t_bif)
    q_inf = 0.5 * emission.rho * emission.U**2
    t_fz = d_parent * math.sqrt(follower_density / q_inf) * t_hat if q_inf > 0.0 else math.nan

    rho_parent = math.nan if parent_mass is None else parent_mass / (math.pi * d_parent**3 / 6.0)
    mach = emission.U / math.sqrt(GAMMA_AIR * R_AIR * emission.T)
    flags = {
        "r_below_fit": r < R_FIT_RANGE[0],
        "r_above_rstar": r > R_STAR,
        "rarefied_parent": emission.mfp / d_parent > settings.kn_continuum,
        "rarefied_follower": emission.mfp / d_follower > settings.kn_continuum,
        "not_hypersonic": mach < settings.mach_min,
        "unequal_density": (parent_mass is not None
                            and abs(rho_parent - follower_density) > DENSITY_TOL * follower_density),
        "parent_bulk_density": rho_parent,
    }
    return AffixResult(zone, r, theta_deg, t_inn, t_bif, t_hat, t_fz, flags)


def _unit(v, what):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if not n > 0.0:
        raise ValueError(f"{what} must be a non-zero vector")
    return v / n


def release_offset(theta_deg, d_parent, v_hat, up):
    """Release point on the parent's surface relative to its centroid, for a release angle theta measured from the
    direction opposite to the velocity: (d/2)(-cos(theta) v_hat + sin(theta) n_hat), where n_hat is the part of `up`
    perpendicular to the velocity, normalised. theta = 180 is the stagnation point (d/2) v_hat ahead of the centroid.
    Raises ValueError if `up` is parallel to the velocity (the side is then undefined)."""
    v_hat = _unit(v_hat, "v_hat")
    up = np.asarray(up, dtype=float)
    perp = up - (up @ v_hat) * v_hat
    if not np.linalg.norm(perp) > 1e-12 * max(np.linalg.norm(up), 1.0):
        raise ValueError("up is parallel to the velocity, so the release side is undefined")
    n_hat = perp / np.linalg.norm(perp)
    th = math.radians(theta_deg)
    return 0.5 * d_parent * (-math.cos(th) * v_hat + math.sin(th) * n_hat)


def far_wake_offset(d_parent, v_hat):
    """Start of a far-wake release: on the wake axis two parent diameters behind the centroid, -2 d v_hat."""
    return -2.0 * d_parent * _unit(v_hat, "v_hat")
