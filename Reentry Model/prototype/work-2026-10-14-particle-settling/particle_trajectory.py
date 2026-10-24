"""Trajectory of a single sphere from an initial velocity vector and altitude, with the drag
coefficient re-evaluated at every instant from the local air state.

The sphere is launched into the atmosphere and integrated until it reaches the ground (or a stop
altitude / stop time).  The point of the script is the drag law: over a descent from the upper
mesosphere a micron-scale particle crosses the entire rarefaction map, so no single textbook
expression applies.  Two numbers decide which physics is in force:

    Kn = lambda / d          rarefaction  - how large the air's mean free path is against the sphere
    s  = V / sqrt(2 R T)     speed ratio  - how fast the sphere moves against the air's own thermal motion

The drag coefficient is composed as

    C_D = C_D,cont + (C_D,fm - C_D,cont) * f(Kn)

exactly the structure reentry_model/aero.py uses, but with both branches valid at *any* speed, not
only the hypersonic ones SESAM's tables cover:

  * C_D,fm  is Schaaf & Chambre's free-molecule sphere, diffuse reflection with full thermal
    accommodation.  It is valid for every speed ratio: s -> infinity gives C_D -> 2, and s -> 0
    gives the Epstein creeping-flow drag, so the same formula covers the hypersonic ejection and
    the millimetre-per-second settling months later.

  * C_D,cont is the continuum sphere: Schiller-Naumann at low Mach (which contains Stokes' 24/Re as
    its Re -> 0 limit) blended into DRAMA's ATDB_SPHERE compressible table above Mach 1.

  * f(Kn) is the rarefaction bridge.  The default, "cunningham", is *derived* so that in the
    low-speed limit the composed C_D reproduces the Cunningham-Millikan-Davies slip correction
    exactly -- the standard, Millikan-calibrated result for aerosol settling, and what
    terminal_velocity.py uses.  The alternative, "sesam", is reentry_model.aero.SesamErf, fitted to
    SESAM's own hypersonic drag output; it is offered for comparison and is NOT calibrated at low
    speed.

Air is the US76 table committed in this repository.  Viscosity is Sutherland's law.  The mean free
path is the viscosity-based one that the Cunningham correlation and Schaaf & Chambre both assume,
not the hard-sphere lambda reentry_model.aero uses to reproduce SESAM's Knudsen column; the
hard-sphere value is carried in the output as a cross-reference column only.

Geometry is a flat-Earth vertical plane: x downrange, y up.  This is checked, not assumed -- the
script reports the distance and time over which the sphere loses its initial speed, and warns if
either is large enough for Earth curvature or rotation to matter.
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for cand in (REPO, "/Users/ashajain/Documents/University Documents /MIT Graduate Work/"
                   "Research/Space Sustainability/Particle Wake Evolution"):
    if os.path.isdir(os.path.join(cand, "reentry_model")) and cand not in sys.path:
        sys.path.insert(0, cand)

from reentry_model import aero
from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.constants import G0, K_BOLTZMANN, M_BAR_AIR, WGS84_A

R_AIR = K_BOLTZMANN / M_BAR_AIR            # specific gas constant of air [J/(kg K)]
GAMMA = 1.4
MU_REF, T_REF, S_SUTH = 1.716e-5, 273.15, 110.4     # Sutherland's law for air
LAMBDA_REF, T_LAMBDA_REF, P_LAMBDA_REF = 6.65e-8, 293.15, 101325.0


# ---------------------------------------------------------------- air properties

def viscosity(T):
    """Sutherland's law [Pa s]; within ~2 % over 150-600 K."""
    return MU_REF * (T / T_REF) ** 1.5 * (T_REF + S_SUTH) / (T + S_SUTH)


def mean_free_path(T, p):
    """Willeke (1976) scaling of the reference mean free path, the convention Cunningham's
    coefficients were fitted with (0.0665 um at 293.15 K, 101.325 kPa)."""
    return (LAMBDA_REF * (P_LAMBDA_REF / p) * (T / T_LAMBDA_REF)
            * (1.0 + S_SUTH / T_LAMBDA_REF) / (1.0 + S_SUTH / T))


def gravity(h):
    return G0 * (WGS84_A / (WGS84_A + h)) ** 2


# ---------------------------------------------------------------- drag branches

def cd_free_molecular(s, temperature_ratio=1.0):
    """Schaaf & Chambre (1961) sphere in free-molecule flow, diffuse reflection, full accommodation.

    s is the molecular speed ratio V / sqrt(2 R T_inf); temperature_ratio is T_wall / T_inf.
    Exact for every speed ratio at Kn >> 1:  C_D -> 2 as s -> inf (hypersonic), and
    C_D -> K/s with K = 16/(3 sqrt(pi)) + (2/3) sqrt(pi T_w/T_inf) as s -> 0 (Epstein creeping drag).
    """
    if s <= 0.0:
        return float("inf")
    if s < 1e-3:          # the full expression cancels 1/s^3 terms; use the exact s -> 0 limit
        return (16.0 / (3.0 * math.sqrt(math.pi))
                + 2.0 * math.sqrt(math.pi * temperature_ratio) / 3.0) / s
    return (math.exp(-s * s) * (2.0 * s * s + 1.0) / (math.sqrt(math.pi) * s ** 3)
            + (4.0 * s ** 4 + 4.0 * s * s - 1.0) * math.erf(s) / (2.0 * s ** 4)
            + 2.0 * math.sqrt(math.pi * temperature_ratio) / (3.0 * s))


_ATDB = aero.SphereDragTables.from_json()


def cd_continuum(re, ma):
    """Continuum sphere. Schiller-Naumann below Ma 0.6 (its Re -> 0 limit is Stokes' 24/Re),
    DRAMA's ATDB_SPHERE compressible table above Ma 1, linear in Ma between."""
    if re <= 0.0:
        return float("inf")
    incompressible = (24.0 / re) * (1.0 + 0.15 * re ** 0.687) if re < 1000.0 else 0.44
    if ma <= 0.6:
        return incompressible
    compressible = _ATDB.cd_continuum(max(ma, 1.0))
    if ma >= 1.0:
        return compressible
    w = (ma - 0.6) / 0.4
    return (1.0 - w) * incompressible + w * compressible


def cunningham(kn):
    """Cunningham-Millikan-Davies slip correction, Kn = 2 lambda / d."""
    return 1.0 + kn * (1.257 + 0.400 * math.exp(-1.100 / kn))


_SESAM_ERF = aero.SesamErf()


def cd_slow(re, ma, kn):
    """Low-speed branch: Stokes with the Schiller-Naumann finite-Re factor and the Cunningham slip
    correction. Millikan-calibrated, valid at ANY Knudsen number but only for s << 1. This is
    exactly the law terminal_velocity.py uses."""
    return cd_continuum(re, min(ma, 0.6)) / cunningham(kn)


def cd_fast(re, ma, kn, s, temperature_ratio):
    """High-speed branch: the standard additive rarefaction bridge, the structure
    reentry_model/aero.py uses. Exact at Kn >> 1 (Schaaf & Chambre) and at Kn << 1 (continuum
    sphere); valid at any speed ratio, but the bridge itself is calibrated on hypersonic data."""
    cont = cd_continuum(re, ma)
    fm = cd_free_molecular(s, temperature_ratio)
    f = _SESAM_ERF(kn / 2.0)             # SesamErf is in the repo's Kn = lambda/D convention
    return cont + (fm - cont) * f


# The two branches are blended on the speed ratio, not on Knudsen number. Below S_LO the Cunningham
# law is used; above S_HI the rarefaction bridge. The window sits far below s = 1 because that is
# where the free-molecule drag is still in its creeping-flow limit, so the two branches describe the
# same physics there and the seam costs only the known ~5 % offset between Millikan's coefficients
# and ideal diffuse reflection. The blend is NOT trustworthy at moderate Knudsen number
# (0.1 < Kn < 20) combined with moderate speed ratio; `check_blend_domain` reports whether a
# trajectory ever entered that corner.
S_LO, S_HI = 0.05, 0.5


def speed_blend(s):
    if s <= S_LO:
        return 0.0
    if s >= S_HI:
        return 1.0
    x = (math.log10(s) - math.log10(S_LO)) / (math.log10(S_HI) - math.log10(S_LO))
    return x * x * (3.0 - 2.0 * x)          # smoothstep, zero slope at both ends


def drag_coefficient(V, rho, T, p, d, temperature_ratio=1.0, law="blend"):
    """Composed C_D plus the diagnostics that say which regime produced it."""
    mu = viscosity(T)
    lam = mean_free_path(T, p)
    kn = 2.0 * lam / d                                  # aerosol convention, what Cunningham uses
    re = rho * V * d / mu
    ma = V / math.sqrt(GAMMA * R_AIR * T)
    s = V / math.sqrt(2.0 * R_AIR * T)
    slow = cd_slow(re, ma, kn)
    fast = cd_fast(re, ma, kn, s, temperature_ratio)
    w = speed_blend(s)
    if law == "slow":
        cd = slow
    elif law == "fast":
        cd = fast
    else:
        cd = (1.0 - w) * slow + w * fast
    return cd, dict(Kn=kn, Re=re, Ma=ma, s=s, cd_slow=slow, cd_fast=fast, w=w, mu=mu, lam=lam,
                    cd_fm=cd_free_molecular(s, temperature_ratio), cd_cont=cd_continuum(re, ma))


LAWS = ("blend", "slow", "fast")


# ---------------------------------------------------------------- trajectory

class Particle:
    def __init__(self, diameter, density):
        self.d = diameter
        self.rho_p = density
        self.mass = density * math.pi * diameter ** 3 / 6.0
        self.area = math.pi * diameter ** 2 / 4.0


def integrate(particle, vx0, vy0, h0, atm, law="blend", temperature_ratio=1.0,
              stop_altitude=0.0, max_time=None, rtol=1e-9, atol=1e-12):
    """Integrate in a flat-Earth vertical plane: state = (x, y, vx, vy), y the altitude."""

    def air(y):
        return atm.state(0.0, min(max(y, 0.0), 149.9e3), 0.0, 0.0)

    def rhs(t, u):
        x, y, vx, vy = u
        fs = air(y)
        V = math.hypot(vx, vy)
        g = gravity(max(y, 0.0))
        if V <= 0.0:
            return [vx, vy, 0.0, -g * (1.0 - fs.rho / particle.rho_p)]
        cd, _ = drag_coefficient(V, fs.rho, fs.T, fs.p, particle.d, temperature_ratio, law)
        k = 0.5 * fs.rho * V * cd * particle.area / particle.mass     # drag accel per unit velocity
        return [vx, vy,
                -k * vx,
                -k * vy - g * (1.0 - fs.rho / particle.rho_p)]        # buoyancy included

    def hit(t, u):
        return u[1] - stop_altitude
    hit.terminal, hit.direction = True, -1

    t_end = max_time if max_time is not None else 1e9
    sol = solve_ivp(rhs, (0.0, t_end), [0.0, h0, vx0, vy0], method="LSODA",
                    events=hit, rtol=rtol, atol=atol, dense_output=True, max_step=np.inf)
    return sol


def relaxation_time(particle, h, atm):
    """Velocity relaxation time m/(3 pi mu d / Cc): how long the sphere takes to forget its own
    velocity and adopt the local terminal velocity."""
    fs = atm.state(0.0, h, 0.0, 0.0)
    mu = viscosity(fs.T)
    cc = cunningham(2.0 * mean_free_path(fs.T, fs.p) / particle.d)
    return particle.rho_p * particle.d ** 2 * cc / (18.0 * mu)


def integrate_hybrid(particle, vx0, vy0, h0, atm, law="blend", temperature_ratio=1.0,
                     stop_altitude=0.0, rtol=1e-9, atol=1e-12, converged=1e-6):
    """Ballistic phase as the full 4-state ODE, then the quasi-static descent.

    Once the relaxation time is negligible against the time to fall a scale height the sphere sits
    on its terminal velocity, the momentum equation is slaved to the altitude, and integrating it
    explicitly is both pointless and ruinously stiff -- for a 1 um sphere the descent lasts years
    while the relaxation time is tens of microseconds, and LSODA gives up part way down. The second
    phase therefore integrates the reduced equation dh/dt = -v_t(h), which is exact in that limit.
    The switch is verified, not assumed: the phase-1 solution must already match the local terminal
    velocity to `converged` in relative terms, and phase 1 is extended until it does.
    """
    t_ball = 500.0 * relaxation_time(particle, h0, atm)
    sol1 = mismatch = h_sw = None
    for _ in range(8):
        sol1 = integrate(particle, vx0, vy0, h0, atm, law, temperature_ratio,
                         stop_altitude, t_ball, rtol, atol)
        if sol1.t_events[0].size:                       # already on the ground: nothing to slave
            return sol1, None, dict(switched=False, t_switch=None, mismatch=None, reason="ground")
        x, h_sw, vx, vy = sol1.y[:, -1]
        V = math.hypot(vx, vy)
        vt = terminal_velocity(particle, max(h_sw, 0.0), atm, law, temperature_ratio)
        mismatch = abs(V / vt - 1.0)
        if mismatch < converged and vy < 0.0:
            break
        t_ball *= 4.0
    t_sw = sol1.t[-1]

    def rhs(t, u):
        h = u[0]
        return [-terminal_velocity(particle, max(min(h, 149.9e3), 0.0), atm, law, temperature_ratio)]

    def hit(t, u):
        return u[0] - stop_altitude
    hit.terminal, hit.direction = True, -1

    sol2 = solve_ivp(rhs, (t_sw, t_sw + 1e12), [h_sw], method="LSODA", events=hit,
                     rtol=1e-8, atol=1e-6, dense_output=True)
    return sol1, sol2, dict(switched=True, t_switch=t_sw, h_switch=h_sw, mismatch=mismatch,
                            tau_over_descent=relaxation_time(particle, h_sw, atm)
                            / (7000.0 / max(terminal_velocity(particle, max(h_sw, 0.0), atm,
                                                              law, temperature_ratio), 1e-30)))


class _Stitched:
    """Presents the two phases as one dense solution so `sample` does not care which produced a
    given time. In the quasi-static phase the horizontal velocity is zero and the speed is the
    terminal velocity, which is what the reduced equation integrated."""

    def __init__(self, sol1, sol2, particle, atm, law, temperature_ratio):
        self.sol1, self.sol2 = sol1, sol2
        self.particle, self.atm, self.law, self.tr = particle, atm, law, temperature_ratio
        self.t_switch = sol1.t[-1]
        self.x_switch = sol1.y[0, -1]
        self.t = np.concatenate((sol1.t, sol2.t[1:])) if sol2 is not None else sol1.t
        self.nfev = sol1.nfev + (sol2.nfev if sol2 is not None else 0)
        self.t_events = sol2.t_events if sol2 is not None else sol1.t_events
        self.status = sol2.status if sol2 is not None else sol1.status
        self.success = sol1.success and (sol2 is None or sol2.success)
        self.message = sol2.message if sol2 is not None else sol1.message

    def sol(self, t):
        if self.sol2 is None or t <= self.t_switch:
            return self.sol1.sol(t)
        h = float(self.sol2.sol(t)[0])
        v = terminal_velocity(self.particle, max(h, 0.0), self.atm, self.law, self.tr)
        return np.array([self.x_switch, h, 0.0, -v])


def sample(sol, particle, atm, law, temperature_ratio, n_log=4000):
    """Resample the solution on a log-spaced time grid and attach the regime diagnostics."""
    t_end = sol.t[-1]
    ts = np.unique(np.concatenate((
        [0.0], np.logspace(-6, math.log10(t_end), n_log) if t_end > 1e-6 else [t_end], [t_end])))
    ts = ts[ts <= t_end]
    rows = []
    for t in ts:
        x, y, vx, vy = sol.sol(t)
        fs = atm.state(0.0, min(max(y, 0.0), 149.9e3), 0.0, 0.0)
        V = math.hypot(vx, vy)
        cd, dg = drag_coefficient(max(V, 1e-30), fs.rho, fs.T, fs.p, particle.d,
                                  temperature_ratio, law)
        g = gravity(max(y, 0.0))
        k = 0.5 * fs.rho * V * cd * particle.area / particle.mass
        rows.append(dict(
            t_s=t, x_m=x, h_m=y, vx=vx, vy=vy, V=V, C_D=cd,
            Kn=dg["Kn"], Kn_repo=aero.knudsen(fs.rho, fs.m_bar, particle.d),
            Re=dg["Re"], Ma=dg["Ma"], speed_ratio=dg["s"], blend_w=dg["w"], cd_slow=dg["cd_slow"], cd_fast=dg["cd_fast"],
            cd_fm=dg["cd_fm"], cd_cont=dg["cd_cont"],
            rho_air=fs.rho, T_air=fs.T, decel_g=k * V / 9.80665,
            regime=("free-molecule" if dg["Kn"] > 20.0 else
                    "transitional" if dg["Kn"] > 0.2 else
                    "slip" if dg["Kn"] > 0.02 else "continuum"),
        ))
    return rows


def check_blend_domain(rows):
    """The composed law is trustworthy except at moderate Knudsen number together with moderate
    speed ratio, where neither branch is calibrated. Report whether the trajectory went there, and
    how far apart the two branches were wherever the blend actually mixed them."""
    suspect = [r for r in rows if 0.1 < r["Kn"] < 20.0 and r["speed_ratio"] > S_LO]
    mixed = [r for r in rows if 0.0 < r["blend_w"] < 1.0]
    gap = max((abs(r["cd_fast"] / r["cd_slow"] - 1.0) for r in mixed), default=0.0)
    kn_mixed = [r["Kn"] for r in mixed]
    return dict(suspect=len(suspect), mixed=len(mixed), worst_branch_gap=gap,
                kn_range=(min(kn_mixed), max(kn_mixed)) if kn_mixed else None)


def terminal_velocity(particle, h, atm, law="blend", temperature_ratio=1.0):
    """Steady-state fall speed at altitude h: the v where drag balances weight less buoyancy."""
    fs = atm.state(0.0, h, 0.0, 0.0)
    g = gravity(h)
    w = particle.mass * g * (1.0 - fs.rho / particle.rho_p)

    def residual(log_v):
        v = math.exp(log_v)
        cd, _ = drag_coefficient(v, fs.rho, fs.T, fs.p, particle.d, temperature_ratio, law)
        return 0.5 * fs.rho * v * v * cd * particle.area - w

    return math.exp(brentq(residual, math.log(1e-12), math.log(1e4), xtol=1e-12, rtol=1e-14))


# ---------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--diameter", type=float, default=10.0, help="sphere diameter [um] (default 10)")
    ap.add_argument("--density", type=float, default=2813.0,
                    help="sphere density [kg/m3] (default 2813, AA7075)")
    ap.add_argument("--vx", type=float, default=6.5, help="initial horizontal velocity [km/s]")
    ap.add_argument("--vy", type=float, default=0.5,
                    help="initial vertical velocity [km/s], positive up")
    ap.add_argument("--altitude", type=float, default=80.0, help="initial altitude [km]")
    ap.add_argument("--law", choices=LAWS, default="blend",
                    help="blend (default), or force the slow/fast branch alone")
    ap.add_argument("--wall-temperature-ratio", type=float, default=1.0,
                    help="T_wall/T_air in the free-molecule branch (default 1, isothermal)")
    ap.add_argument("--stop-altitude", type=float, default=0.0, help="stop altitude [km]")
    ap.add_argument("--max-time", type=float, default=None, help="stop time [s]")
    ap.add_argument("--full-ode", action="store_true",
                    help="integrate the momentum equation all the way down instead of slaving it "
                         "to the terminal velocity once relaxation is instantaneous (stiff; fails "
                         "for small spheres, but useful as a check)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args(argv)

    d = a.diameter * 1e-6
    particle = Particle(d, a.density)
    atm = US76TableAtmosphere()
    vx0, vy0, h0 = a.vx * 1e3, a.vy * 1e3, a.altitude * 1e3
    V0 = math.hypot(vx0, vy0)

    print("sphere: d = {:g} um, rho = {:g} kg/m3, mass = {:.4e} kg, frontal area = {:.4e} m2".format(
        a.diameter, a.density, particle.mass, particle.area))
    print("launch: Vx = {:g} km/s, Vy = {:+g} km/s, |V| = {:.1f} m/s at h = {:g} km".format(
        a.vx, a.vy, V0, a.altitude))
    cd0, dg0 = drag_coefficient(V0, *[getattr(atm.state(0.0, h0, 0.0, 0.0), k)
                                      for k in ("rho", "T", "p")], d,
                                a.wall_temperature_ratio, a.law)
    print("at launch: Kn = {:.3e}, speed ratio s = {:.2f}, Ma = {:.2f}, Re = {:.3e}".format(
        dg0["Kn"], dg0["s"], dg0["Ma"], dg0["Re"]))
    print("           C_D = {:.4f}  (free-molecule {:.4f}, continuum {:.4g}, "
          "slow branch {:.4g}, fast branch {:.4f}, speed blend w = {:.3f})".format(
              cd0, dg0["cd_fm"], dg0["cd_cont"], dg0["cd_slow"], dg0["cd_fast"], dg0["w"]))
    print("           law = {}".format(a.law))

    if a.full_ode or a.max_time is not None:
        sol = integrate(particle, vx0, vy0, h0, atm, a.law, a.wall_temperature_ratio,
                        a.stop_altitude * 1e3, a.max_time)
        info = dict(switched=False)
    else:
        sol1, sol2, info = integrate_hybrid(particle, vx0, vy0, h0, atm, a.law,
                                            a.wall_temperature_ratio, a.stop_altitude * 1e3)
        sol = _Stitched(sol1, sol2, particle, atm, a.law, a.wall_temperature_ratio) \
            if sol2 is not None else sol1
    print("\nintegration: {}, {} steps, {} rhs evaluations, t_end = {:.6g} s = {:.4g} days".format(
        "ground reached" if sol.t_events[0].size else sol.status, sol.t.size, sol.nfev,
        sol.t[-1], sol.t[-1] / 86400.0))
    if info.get("switched"):
        print("  ballistic phase integrated to t = {:.4g} s (h = {:.3f} km), where the speed "
              "already matched the local terminal velocity to {:.2e} relative;".format(
                  info["t_switch"], info["h_switch"] / 1e3, info["mismatch"]))
        print("  below that the momentum equation is slaved to it (relaxation time / time to fall "
              "one scale height = {:.2e})".format(info["tau_over_descent"]))
    if not sol.success and not sol.t_events[0].size:
        print("  WARNING: {}".format(sol.message))

    rows = sample(sol, particle, atm, a.law, a.wall_temperature_ratio)
    tag = a.tag or "d{:g}um_vx{:g}_vy{:+g}_h{:g}km_{}".format(a.diameter, a.vx, a.vy, a.altitude, a.law)
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "trajectory_{}.csv".format(tag))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("wrote {} ({} rows)".format(out, len(rows)))

    chk = check_blend_domain(rows)
    print("\nblend-domain check:")
    if chk["kn_range"]:
        print("  the two branches were mixed over Kn = {:.3g} to {:.3g}; worst disagreement "
              "between them there: {:.1f} %".format(chk["kn_range"][0], chk["kn_range"][1],
                                                    100 * chk["worst_branch_gap"]))
    else:
        print("  the blend never mixed the branches (every sample was purely one or the other)")
    if chk["suspect"]:
        print("  WARNING: {} sampled states had moderate Kn (0.1-20) with speed ratio above {:g}; "
              "the composed law is not calibrated there".format(chk["suspect"], S_LO))
    else:
        print("  no sampled state combined moderate Kn (0.1-20) with a speed ratio above "
              "{:g}, so the untrusted corner was never entered".format(S_LO))

    h_final = rows[-1]["h_m"]
    v_term = terminal_velocity(particle, max(h_final, 0.0), atm, a.law, a.wall_temperature_ratio)
    print("  final speed {:.6g} m/s vs terminal velocity {:.6g} m/s at that altitude "
          "(ratio {:.4f})".format(rows[-1]["V"], v_term, rows[-1]["V"] / v_term))
    return sol, rows, particle, atm, a


if __name__ == "__main__":
    main()
