"""Copied 2026-10-26 from prototype/work-2026-10-14-particle-settling/droplet_heating_sweep.py (verified there
2026-10-14); only the imports, the parent file and the default output folder changed. The droplet model uses its
Schaaf & Chambre free-molecule heat exchange.

Molten AA7075 droplets released into the wake of the descending parent satellite: how much each one
heats up while it decelerates, over a diameter x release-altitude matrix.

The momentum physics is particle_trajectory.py's, reused rather than rewritten: its composed drag law
(`drag_coefficient`), gravity, buoyancy and, for the slow settling tail, `integrate_hybrid`.  What this
script adds is

  * the parent state keyed on altitude (`parent_state_at_altitude`), read through
    sphere_sweep.load_parent from the DRAMA descent of the parent's Main Body;

  * a lumped-capacitance energy balance carried alongside the momentum equation,

        m c_p dT/dt = Q_conv(V, air, T) - eps sigma 4 pi r^2 (T^4 - T_bg^4),

    with Q_conv the free-molecule heat transfer to a sphere at *any* speed ratio (Schaaf & Chambre's
    flat-element fluxes integrated over the sphere, `fm_sphere_heat`), and the droplet temperature fed
    back into the free-molecule drag through T_wall / T_air.

Phase 1 integrates (x, y, vx, vy, T) plus running energy tallies until the molecular speed ratio falls
below particle_trajectory.S_LO.  Below that the drag is the Cunningham branch alone, which does not
depend on the wall temperature, so the momentum is handed to `integrate_hybrid` unchanged from that
state; the temperature is continued along the stitched trajectory by itself.

No melting, freezing, vaporisation or mass loss is modelled: the droplet keeps its mass and the
liquid heat capacity at every temperature.  Cases that would vaporise are flagged, not corrected.
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import erf

from reentry_model import OUTPUT_DIR
from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.droplet_reentry import parent_trajectory
from reentry_model.droplet_reentry import particle_trajectory as pt

SIGMA_SB = 5.670374419e-8
PARENT_JSON = parent_trajectory.DEFAULT_PARENT
PARENT_OBJECT = parent_trajectory.DEFAULT_PARENT_OBJECT

# ---------------------------------------------------------------- molten AA7075

RHO_LIQUID = 2400.0     # kg/m3; prototype/proto3 AA7075.json `liquid.density` (Smithells, pure Al near the liquidus)
CP_LIQUID = 1180.0      # J/(kg K); liquid Al 31.75 J/(mol K) (NIST-JANAF) / 0.026982 kg/mol = 1177
K_LIQUID = 90.0         # W/(m K); liquid Al near the liquidus (Mills 2002), only for the Biot check
T_RELEASE = 908.0       # K; AA7075 liquidus (ASM Handbook Vol. 2), prototype/proto3 AA7075_range.json
T_SOLIDUS = 750.0       # K; same source
T_BG = 250.0            # K; radiative background (Earth/atmosphere seen from the mesosphere)
EPS_DEFAULT = 0.10      # liquid Al total hemispherical emissivity is about 0.05-0.12 near 1000-1700 K
ALPHA_DEFAULT = 1.0     # full thermal accommodation, the assumption in the drag law's free-molecule branch
ZETA_AIR = 2.0          # internal (rotational) degrees of freedom of N2/O2, vibration frozen

T_BOIL_AL = 2743.0      # K; normal (1 atm) boiling point of aluminium
T_BOIL_MG = 1363.0      # K; magnesium, 2.5 % of AA7075
T_BOIL_ZN = 1180.0      # K; zinc, 5.6 % of AA7075
H_VAP_AL = 10.9e6       # J/kg; 294 kJ/mol / 0.026982 kg/mol


def al_saturation_temperature(p):
    """Temperature at which liquid Al's vapour pressure equals p [Pa]: Alcock, Itkin & Horrigan (1984),
    log10(p/atm) = 5.911 - 16211/T, which gives 2743 K at 1 atm."""
    return 16211.0 / (5.911 - math.log10(p / 101325.0))


# ---------------------------------------------------------------- parent state keyed on altitude

def parent_state_at_altitude(parent, h_km):
    """Companion of sphere_sweep.interpolate_state keyed on altitude: linear in altitude on the
    descending branch, or None outside it.  Returns the speed, the flight path angle and the planar
    velocity components Vx = V cos(gamma), Vy = V sin(gamma) (gamma negative downward)."""
    branch = parent.branch
    for a, b in zip(branch, branch[1:]):
        if a["altitude"] >= h_km >= b["altitude"]:
            dh = a["altitude"] - b["altitude"]
            f = (a["altitude"] - h_km) / dh if dh else 0.0
            v = parent_trajectory.lerp(a["velocity"], b["velocity"], f)
            gam = parent_trajectory.lerp(a["path"], b["path"], f)
            return dict(altitude_km=h_km, velocity_kms=v, flight_path_deg=gam,
                        time_s=parent_trajectory.lerp(a["time"], b["time"], f),
                        downrange_km=parent_trajectory.lerp(a["downrange"], b["downrange"], f),
                        vx=v * 1e3 * math.cos(math.radians(gam)),
                        vy=v * 1e3 * math.sin(math.radians(gam)))
    return None


# ---------------------------------------------------------------- free-molecule heat transfer

_GL_X, _GL_W = np.polynomial.legendre.leggauss(48)
# mu = cos(angle from the stagnation point), split at mu = 0 where the integrand has its kink at high s
MU = np.concatenate((0.5 * (_GL_X - 1.0), 0.5 * (_GL_X + 1.0)))
MU_W = np.concatenate((0.5 * _GL_W, 0.5 * _GL_W))
SQPI = math.sqrt(math.pi)


def fm_element_fluxes(rho, T, s, sn):
    """Incident fluxes per unit area on a flat element whose inward normal makes speed ratio sn with
    the flow (Schaaf & Chambre 1961; Bird 1994 sec. 7.2).  Returns (mass flux, translational energy
    flux).  beta = 1/sqrt(2RT)."""
    beta = 1.0 / math.sqrt(2.0 * pt.R_AIR * T)
    e = np.exp(-sn * sn)
    g = SQPI * sn * (1.0 + erf(sn))
    mass = rho / (2.0 * SQPI * beta) * (e + g)
    etr = rho / (4.0 * SQPI * beta ** 3) * ((s * s + 2.0) * e + (s * s + 2.5) * g)
    return mass, etr


def fm_sphere_heat(rho, T, V, Tw, d, alpha=ALPHA_DEFAULT):
    """Net free-molecule heat transfer into a sphere [W] at any speed ratio, diffuse re-emission with
    thermal accommodation alpha: q = alpha (E_tr + E_int - (2 + zeta/2) R Tw * mass flux) on every
    element, integrated over the sphere (convex, so no element sees another's re-emission)."""
    s = V / math.sqrt(2.0 * pt.R_AIR * T)
    mass, etr = fm_element_fluxes(rho, T, s, s * MU)
    q = alpha * (etr + mass * pt.R_AIR * (0.5 * ZETA_AIR * T - (2.0 + 0.5 * ZETA_AIR) * Tw))
    r = 0.5 * d
    return 2.0 * math.pi * r * r * float(np.dot(MU_W, q))


def fm_high_speed_limit(rho, V, d, alpha=ALPHA_DEFAULT):
    """alpha * 1/2 rho V^3 * frontal area: the s -> infinity limit of fm_sphere_heat."""
    return alpha * 0.5 * rho * V ** 3 * math.pi * d * d / 4.0


# ---------------------------------------------------------------- one droplet

def _air(atm, y):
    return atm.state(0.0, min(max(y, 0.0), 149.9e3), 0.0, 0.0)


def run_case(d_um, h_km, vx0, vy0, eps=EPS_DEFAULT, alpha=ALPHA_DEFAULT, history=False):
    atm = US76TableAtmosphere()
    d = d_um * 1e-6
    prt = pt.Particle(d, RHO_LIQUID)
    m, A_s = prt.mass, math.pi * d * d
    mcp = m * CP_LIQUID

    def physics(y, vx, vy, T):
        fs = _air(atm, y)
        V = math.hypot(vx, vy)
        cd, dg = pt.drag_coefficient(max(V, 1e-30), fs.rho, fs.T, fs.p, d, T / fs.T, "blend")
        k = 0.5 * fs.rho * V * cd * prt.area / m
        q_conv = fm_sphere_heat(fs.rho, fs.T, V, T, d, alpha)
        q_rad = eps * SIGMA_SB * A_s * (T ** 4 - T_BG ** 4)
        return fs, V, cd, dg, k, q_conv, q_rad

    # state: x, y, vx, vy, T, E_in, E_in(s>3), E_highspeed_limit, E_rad, E_conv_out  (energies per kg)
    def rhs(t, u):
        x, y, vx, vy, T = u[:5]
        fs, V, cd, dg, k, q_conv, q_rad = physics(y, vx, vy, T)
        g = pt.gravity(max(y, 0.0))
        qin = max(q_conv, 0.0) / m
        return [vx, vy, -k * vx, -k * vy - g * (1.0 - fs.rho / prt.rho_p),
                (q_conv - q_rad) / mcp,
                qin, qin if dg["s"] > 3.0 else 0.0,
                fm_high_speed_limit(fs.rho, V, d, alpha) / m,
                q_rad / m, max(-q_conv, 0.0) / m]

    def slow(t, u):
        fs = _air(atm, u[1])
        return math.hypot(u[2], u[3]) / math.sqrt(2.0 * pt.R_AIR * fs.T) - pt.S_LO
    slow.terminal, slow.direction = True, -1

    def ground(t, u):
        return u[1]
    ground.terminal, ground.direction = True, -1

    u0 = [0.0, h_km * 1e3, vx0, vy0, T_RELEASE, 0, 0, 0, 0, 0]
    t_cap = 1e7
    sol1 = solve_ivp(rhs, (0.0, t_cap), u0, method="LSODA", events=(slow, ground),
                     rtol=1e-9, atol=[1e-6, 1e-6, 1e-9, 1e-9, 1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3],
                     dense_output=True)
    if not sol1.success:
        raise RuntimeError("phase 1 failed for d={} um h={} km: {}".format(d_um, h_km, sol1.message))
    t1 = sol1.t[-1]
    x1, y1, vx1, vy1, T1 = sol1.y[:5, -1]
    reached_ground = sol1.t_events[1].size > 0

    # phase 2: momentum by integrate_hybrid from the hand-off state (drag there is independent of T)
    stitched = None
    if not reached_ground:
        s1, s2, info = pt.integrate_hybrid(prt, vx1, vy1, y1, atm)
        stitched = pt._Stitched(s1, s2, prt, atm, "blend", 1.0) if s2 is not None else s1
        # integrate_hybrid stops its settling phase 1e12 s after the switch; a droplet that has not landed by
        # then gets no ground time
        ground_hit = stitched.t_events[0].size > 0
        t_ground = t1 + stitched.t[-1] if ground_hit else float("nan")
        x_ground = x1 + stitched.sol(stitched.t[-1])[0]
        h_last = float(stitched.sol(stitched.t[-1])[1])
        hybrid_mismatch = info.get("mismatch")
    else:
        t_ground, x_ground, hybrid_mismatch = t1, x1, None
        ground_hit, h_last = True, 0.0

    # temperature after the hand-off, along the stitched trajectory, for up to 1e6 s (by then it sits at its
    # equilibrium with the air, and the history holds that last value down to the ground)
    t_liq_cool = t_sol_cool = None
    T_traj = None
    if stitched is not None:
        t2_end = stitched.t[-1]

        def rhs_T(t, u):
            x, y, vx, vy = stitched.sol(min(t, t2_end))
            fs = _air(atm, y)
            V = math.hypot(vx, vy)
            q = fm_sphere_heat(fs.rho, fs.T, V, u[0], d, alpha)
            return [(q - eps * SIGMA_SB * A_s * (u[0] ** 4 - T_BG ** 4)) / mcp]

        def solid(t, u):
            return u[0] - T_SOLIDUS
        solid.terminal, solid.direction = False, -1
        T_traj = solve_ivp(rhs_T, (0.0, min(t2_end, 1e6)), [T1], method="LSODA", events=solid,
                           rtol=1e-8, atol=1e-6, dense_output=True)

    # ------------------------------------------------ diagnostics from a fine sample of phase 1
    ts = np.unique(np.concatenate(([0.0], np.logspace(min(-7.0, math.log10(1e-7 * d_um / 2.0)), math.log10(max(t1, 1e-6)), 1500), [t1])))
    ts = ts[ts <= t1]
    U = sol1.sol(ts)
    rows = []
    for t, u in zip(ts, U.T):
        x, y, vx, vy, T = u[:5]
        fs, V, cd, dg, k, q_conv, q_rad = physics(y, vx, vy, T)
        rows.append(dict(
            t_s=t, x_m=x, h_m=y, vx=vx, vy=vy, V=V, T_K=T, C_D=cd,
            q_conv_W=q_conv, q_rad_W=q_rad, q_conv_over_highspeed=q_conv / max(
                fm_high_speed_limit(fs.rho, V, d, alpha), 1e-300),
            E_in_Jkg=u[5], E_rad_Jkg=u[8], decel_g=k * V / 9.80665,
            Kn=dg["Kn"], Re=dg["Re"], Ma=dg["Ma"], speed_ratio=dg["s"], blend_w=dg["w"],
            cd_slow=dg["cd_slow"], cd_fast=dg["cd_fast"], rho_air=fs.rho, T_air=fs.T, p_air=fs.p,
            T_sat_Al=al_saturation_temperature(fs.p), phase=1))
    if T_traj is not None:
        t2 = np.unique(np.concatenate(([0.0], np.logspace(-6, math.log10(max(t2_end, 2e-6)), 300), [t2_end])))
        t2 = t2[(t2 > 0) & (t2 <= t2_end)]
        for t in t2:
            x, y, vx, vy = stitched.sol(t)
            T = float(T_traj.sol(t)[0]) if t <= T_traj.t[-1] else float(T_traj.y[0, -1])
            fs = _air(atm, y)
            V = math.hypot(vx, vy)
            cd, dg = pt.drag_coefficient(max(V, 1e-30), fs.rho, fs.T, fs.p, d, T / fs.T, "blend")
            q_conv = fm_sphere_heat(fs.rho, fs.T, V, T, d, alpha)
            rows.append(dict(
                t_s=t1 + t, x_m=x1 + x, h_m=y, vx=vx, vy=vy, V=V, T_K=T, C_D=cd,
                q_conv_W=q_conv, q_rad_W=eps * SIGMA_SB * A_s * (T ** 4 - T_BG ** 4),
                q_conv_over_highspeed=q_conv / max(fm_high_speed_limit(fs.rho, V, d, alpha), 1e-300),
                E_in_Jkg=float("nan"), E_rad_Jkg=float("nan"),
                decel_g=0.5 * fs.rho * V * V * cd * prt.area / m / 9.80665,
                Kn=dg["Kn"], Re=dg["Re"], Ma=dg["Ma"], speed_ratio=dg["s"], blend_w=dg["w"],
                cd_slow=dg["cd_slow"], cd_fast=dg["cd_fast"], rho_air=fs.rho, T_air=fs.T, p_air=fs.p,
                T_sat_Al=al_saturation_temperature(fs.p), phase=2))
        if T_traj.t_events[0].size:
            t_sol_cool = t1 + T_traj.t_events[0][0]

    p1 = [r for r in rows if r["phase"] == 1]
    Ts = np.array([r["T_K"] for r in p1])
    i_pk = int(np.argmax(Ts))
    V0 = math.hypot(vx0, vy0)

    def first_time(cond):
        for r in rows:
            if cond(r):
                return r
        return None

    r_e = first_time(lambda r: r["V"] <= V0 / math.e)
    r_99 = first_time(lambda r: abs(r["vx"]) <= 0.01 * abs(vx0))
    # path length to those points (arc length of the sampled phase-1 path)
    xs, hs = np.array([r["x_m"] for r in p1]), np.array([r["h_m"] for r in p1])
    arc = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(xs), np.diff(hs)))))

    def arc_at(r):
        if r is None or r["phase"] != 1:
            return float("nan")
        return float(arc[p1.index(r)])

    # cooling back to the liquidus: first sample after the peak below 908 K
    r_liq = next((r for r in rows[i_pk:] if r["T_K"] < T_RELEASE and r["t_s"] > p1[i_pk]["t_s"]), None)
    E_in, E_hi, E_lim, E_rad, E_cout = sol1.y[5:, -1]
    chk = pt.check_blend_domain(p1)
    over_sat = max(r["T_K"] - r["T_sat_Al"] for r in rows)
    # internal temperature non-uniformity at the peak-flux instant: q_frontal * r / k
    q_peak = max(r["q_conv_W"] for r in p1) / (math.pi * d * d / 4.0)
    dT_internal = q_peak * 0.5 * d / K_LIQUID
    t_cond = (0.5 * d) ** 2 / (K_LIQUID / (RHO_LIQUID * CP_LIQUID))
    # energy-weighted Knudsen number: how rarefied the flow was while the heat actually arrived
    Ein = np.array([r["E_in_Jkg"] for r in p1])
    dE = np.diff(Ein)
    kn_mid = np.array([0.5 * (a["Kn"] + b["Kn"]) for a, b in zip(p1, p1[1:])])
    frac_kn_lt10 = float(dE[kn_mid < 10.0].sum() / max(dE.sum(), 1e-300))
    kn_ew = float(np.exp(np.sum(dE * np.log(kn_mid)) / max(dE.sum(), 1e-300)))
    # time over which the bulk of the heat arrived (5 %-95 % of E_in)
    tt = np.array([r["t_s"] for r in p1])
    t05 = float(np.interp(0.05 * E_in, Ein, tt)) if E_in > 0 else float("nan")
    t95 = float(np.interp(0.95 * E_in, Ein, tt)) if E_in > 0 else float("nan")

    summary = dict(
        d_um=d_um, h0_km=h_km, V0_ms=V0, gamma0_deg=math.degrees(math.atan2(vy0, vx0)),
        eps=eps, alpha=alpha, mass_kg=m,
        T_peak_K=float(Ts[i_pk]), dT_peak_K=float(Ts[i_pk] - T_RELEASE),
        t_T_peak_s=p1[i_pk]["t_s"], h_T_peak_km=p1[i_pk]["h_m"] / 1e3,
        T_handoff_K=float(T1),
        E_in_Jkg=float(E_in), E_rad_Jkg=float(E_rad), E_conv_out_Jkg=float(E_cout),
        E_highspeed_limit_Jkg=float(E_lim), frac_E_in_s_gt3=float(E_hi / E_in) if E_in > 0 else float("nan"),
        E_in_over_vaporise=float(E_in / (CP_LIQUID * (T_BOIL_AL - T_RELEASE) + H_VAP_AL)),
        E_in_over_KE0=float(E_in / (0.5 * V0 * V0)),
        t_heat_05_s=t05, t_heat_95_s=t95,
        peak_decel_g=max(r["decel_g"] for r in p1),
        t_e_s=r_e["t_s"] if r_e else float("nan"), x_e_km=r_e["x_m"] / 1e3 if r_e else float("nan"),
        path_e_km=arc_at(r_e) / 1e3, h_e_km=r_e["h_m"] / 1e3 if r_e else float("nan"),
        t_stop_s=r_99["t_s"] if r_99 else float("nan"),
        x_stop_km=r_99["x_m"] / 1e3 if r_99 else float("nan"),
        path_stop_km=arc_at(r_99) / 1e3, h_stop_km=r_99["h_m"] / 1e3 if r_99 else float("nan"),
        curvature_drop_at_stop_m=((r_99["x_m"] if r_99 else x_ground) ** 2) / (2 * 6.371e6),
        x_ground_km=x_ground / 1e3, t_ground_s=t_ground, t_ground_days=t_ground / 86400.0,
        t_cool_below_liquidus_s=r_liq["t_s"] if r_liq else float("nan"),
        t_cool_below_solidus_s=t_sol_cool if t_sol_cool else float("nan"),
        Kn_min=min(r["Kn"] for r in rows), Kn_max=max(r["Kn"] for r in rows),
        Kn_release=p1[0]["Kn"], Kn_energy_weighted=kn_ew, frac_E_in_Kn_lt10=frac_kn_lt10,
        s_release=p1[0]["speed_ratio"], s_min=min(r["speed_ratio"] for r in rows),
        Ma_release=p1[0]["Ma"], Re_release=p1[0]["Re"],
        blend_suspect_samples=chk["suspect"], blend_warning=chk["suspect"] > 0,
        blend_worst_branch_gap=chk["worst_branch_gap"],
        flag_T_above_Al_boiling_1atm=bool(Ts[i_pk] > T_BOIL_AL),
        flag_T_above_Al_saturation_local=bool(over_sat > 0), T_minus_Tsat_local_max_K=float(over_sat),
        flag_T_above_Mg_boiling=bool(Ts[i_pk] > T_BOIL_MG), flag_T_above_Zn_boiling=bool(Ts[i_pk] > T_BOIL_ZN),
        flag_energy_exceeds_full_vaporisation=bool(E_in > CP_LIQUID * (T_BOIL_AL - T_RELEASE) + H_VAP_AL),
        biot_dT_internal_K=dT_internal, conduction_time_s=t_cond,
        phase1_end_s=t1, phase1_end_h_km=y1 / 1e3, reached_ground_in_phase1=reached_ground, reached_ground=bool(ground_hit),
        h_at_integration_end_km=h_last / 1e3,
        hybrid_switch_mismatch=hybrid_mismatch if hybrid_mismatch is not None else float("nan"),
    )
    return summary, (rows if history else None)


# ---------------------------------------------------------------- sweep

# ten log-spaced points from 0.01 um up to (not including) 2 um, rounded to three figures, then the original nine
DIAMETERS = (0.01, 0.017, 0.0289, 0.049, 0.0833, 0.141, 0.24, 0.408, 0.693, 1.18,
             2, 6, 10, 14, 18, 50, 100, 150, 250)


def history_name(d, h):
    dtag = "{:03d}".format(int(d)) if float(d).is_integer() else "{:g}".format(d)
    return "droplet_d{}um_h{:02d}km.csv".format(dtag, int(h))
ALTITUDES = tuple(range(90, 49, -2))


def _worker(job):
    d, h, vx, vy, eps, alpha, hist_dir = job
    t0 = time.time()
    try:
        s, rows = run_case(d, h, vx, vy, eps, alpha, history=hist_dir is not None)
    except Exception as exc:            # report and keep the sweep going
        return dict(d_um=d, h0_km=h, eps=eps, alpha=alpha, error=repr(exc))
    s["wall_s"] = time.time() - t0
    if rows:
        path = os.path.join(hist_dir, history_name(d, h))
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--outdir", default=os.path.join(OUTPUT_DIR, "droplet_reentry_output", "droplet_heating_sweep"))
    ap.add_argument("--eps", type=float, nargs="+", default=[EPS_DEFAULT, 0.05, 0.4])
    ap.add_argument("--alpha", type=float, nargs="+", default=[ALPHA_DEFAULT, 0.8])
    ap.add_argument("--diameters", type=float, nargs="+", default=DIAMETERS)
    ap.add_argument("--altitudes", type=float, nargs="+", default=ALTITUDES)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)

    parent = parent_trajectory.load_parent(PARENT_JSON, PARENT_OBJECT)
    os.makedirs(a.outdir, exist_ok=True)
    hist_dir = os.path.join(a.outdir, "histories")
    os.makedirs(hist_dir, exist_ok=True)
    states = {}
    for h in a.altitudes:
        st = parent_state_at_altitude(parent, h)
        if st is None:
            raise SystemExit("altitude {} km is outside the parent's descending branch".format(h))
        states[h] = st
    with open(os.path.join(a.outdir, "parent_states.json"), "w") as fh:
        json.dump(dict(parent=parent.path, sha256=parent.sha256, object=PARENT_OBJECT,
                       states=[states[h] for h in a.altitudes]), fh, indent=1)

    # baseline (first eps, first alpha) writes histories; the variants are summary only
    variants = [(a.eps[0], a.alpha[0])] + [(e, a.alpha[0]) for e in a.eps[1:]] \
        + [(a.eps[0], al) for al in a.alpha[1:]]
    jobs = []
    for vi, (eps, alpha) in enumerate(variants):
        for d in a.diameters:
            for h in a.altitudes:
                jobs.append((int(d) if float(d).is_integer() else d, int(h) if float(h).is_integer() else h,
                             states[h]["vx"], states[h]["vy"], eps, alpha, hist_dir if vi == 0 else None))
    t0 = time.time()
    with Pool(a.workers) as pool:
        out = pool.map(_worker, jobs, chunksize=1)
    print("{} cases in {:.0f} s".format(len(out), time.time() - t0))
    errors = [o for o in out if "error" in o]
    for e in errors:
        print("ERROR", e)
    good = [o for o in out if "error" not in o]
    for (eps, alpha) in variants:
        sel = [o for o in good if o["eps"] == eps and o["alpha"] == alpha]
        name = "summary_eps{:g}_alpha{:g}.csv".format(eps, alpha)
        with open(os.path.join(a.outdir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sel[0].keys()))
            w.writeheader(); w.writerows(sel)
        print("wrote", name, len(sel))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
