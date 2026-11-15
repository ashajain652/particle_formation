"""Similarity-based wall shear driven by Holman & Boyd's DSMC PRESSURE distribution instead of the model's.

Everything else is as in similarity_shear.py: the model's own Cantera equilibrium post-shock reservoir, isentropic
equilibrium expansion to the wall pressure for the edge state, Lees local similarity with a cold wall and variable
rho*mu, then the model's slip and bridging. Only the pressure distribution changes:
    p_w(theta) = p_stag_model * p_HB(theta) / p_HB(0),   p_HB = p_inf + Cp_HB q_inf,
i.e. Holman & Boyd's pressure SHAPE on the model's stagnation pressure (their Cp at 0 deg is about 1.87 against the
model's 1.94). The digitized ln(p) is lightly smoothed (Savitzky-Golay, 9 points, cubic, mirrored about 0 deg)
before differentiating, because du_e/ds amplifies digitizing noise.

Run with the model's environment (Cantera): /Users/ashajain/miniforge3/envs/drama_env/bin/python -I similarity_hb_pressure.py
Writes similarity_hb_pressure.csv next to this file.

With the argument `frozen`, the post-shock gas is CHEMICALLY FROZEN instead: the model's own Cantera construction with
its equilibration switched off (gas.T_EQUILIBRATE raised for this process only), so the air keeps its freestream
N2/O2 composition through the shock and the expansion while vibration stays in equilibrium with temperature (NASA-9
enthalpy). Holman & Boyd's DSMC gas is partly reacted, so the equilibrium and frozen results should bracket it if
the chemistry is what separates the model from DSMC. Writes similarity_hb_pressure_frozen.csv.
"""
import csv, math, sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.signal import savgol_filter

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
REPO = str(Path(__file__).resolve().parents[3])  # Reentry Model/
sys.path.insert(1, REPO)
from similarity_shear import Similarity, R, Tw, U, R_AIR  # noqa: E402
from reentry_model.shared import aero, atmosphere, gas, surface_flow, trajectory  # noqa: E402
from reentry_model.shared.constants import M_BAR_AIR  # noqa: E402

T_INF, DIAMETER = 200.0, 2 * R
CASES = (("0.002", 1.007e-4), ("0.01", 2.014e-5), ("0.05", 3.987e-6), ("0.25", 8.057e-7))


def aero_state(rho):
    p = rho * R_AIR * T_INF
    fs = atmosphere.Freestream(rho, T_INF, p, M_BAR_AIR, np.zeros(3))
    v = np.array([U, 0.0, 0.0])
    kn, ma = aero.knudsen(rho, M_BAR_AIR, DIAMETER), aero.mach(U, T_INF, M_BAR_AIR)
    cd = aero.drag_coefficient(kn, ma, aero.SphereDragTables.from_json(), aero.bridging_by_name("sesam-table"))
    a = -0.5 * rho * U * v * cd * math.pi * R**2 / (2813.0 * 4 / 3 * math.pi * R**3)
    return trajectory.AeroState(float("nan"), 0.0, 0.0, fs, v, U, kn, ma, cd, a, 0.5 * rho * U**2)


FROZEN = len(sys.argv) > 1 and sys.argv[1] == "frozen"


def main():
    if FROZEN:
        gas.T_EQUILIBRATE = 1e12            # no equilibrate() anywhere in stagnation/expand: composition stays N2/O2
    hb = {}
    for r in csv.DictReader(open(HERE / "hb_dsmc.csv")):
        if r["Cp"]:                                   # aft cells may be blank; only 0-90 deg is used
            hb.setdefault(r["case_kn"], {})[int(float(r["theta_deg"]))] = float(r["Cp"])
    flow = surface_flow.SurfaceFlow(rarefied_shear="slip", gamma_pm=1.15, kn_body_shock=0.01)
    mu_wall = gas.wilke_viscosity(Tw, {"N2": 0.79, "O2": 0.21}, flow.air.molar_masses)
    out = []
    th = np.arange(0, 91)
    s = R * np.radians(th); r = R * np.sin(np.radians(th))
    for kn, rho in CASES:
        state = aero_state(rho)
        q, p_inf = 0.5 * rho * U**2, rho * R_AIR * T_INF
        stag = flow.air.stagnation(rho, T_INF, U)
        print(f"Kn {kn}: {'frozen' if FROZEN else 'equilibrium'} post-shock T {stag.T:.0f} K, rho {stag.rho:.3e} "
              f"kg/m3 ({stag.rho / rho:.1f} x freestream), p {stag.p:.0f} Pa, mu {stag.mu:.2e} Pa s")
        p_hb = np.array([p_inf + hb[kn][t] * q for t in th])
        lnp = np.log(p_hb / p_hb[0])
        mirrored = np.concatenate([lnp[:0:-1], lnp])                       # symmetric about the stagnation point
        lnp = savgol_filter(mirrored, 9, 3)[len(th) - 1:]
        # near the stagnation point the pressure must fall as theta^2; digitizing noise can make it rise for a
        # degree or two (Kn 0.25), which stalls u_e. Below 6 deg use ln(p/p0) = c2 th^2 + c4 th^4 fitted on 6-20 deg.
        lnp = lnp - lnp[0]
        fit = np.polyfit(np.radians(th[6:21]) ** 2, lnp[6:21], 2)
        lnp[:6] = np.polyval(fit, np.radians(th[:6]) ** 2) - np.polyval(fit, 0.0)
        p_w = stag.p * np.exp(lnp)
        ue, rhoe = np.zeros(len(th)), np.zeros(len(th))
        for i, p in enumerate(p_w):
            e = flow.air.expand(stag, float(p))
            ue[i] = math.sqrt(max(0.0, 2.0 * (stag.h - e.h))); rhoe[i] = e.rho
        ue[0] = 0.0
        rmw = p_w / (R_AIR * Tw) * mu_wall
        fine = np.linspace(0, s[-1], 20 * 90 + 1)
        xi = np.interp(s, fine, cumulative_trapezoid(np.interp(fine, s, rmw) * np.interp(fine, s, ue)
                                                     * (R * np.sin(fine / R))**2, fine, initial=0.0))
        due = np.gradient(ue, s)
        lam_w, kn_local = surface_flow.wall_knudsen(p_w, np.full(len(th), Tw), np.full(len(th), mu_wall), R)
        tau_fm = state.freestream.rho * U * U * np.sin(np.radians(th)) * np.abs(np.cos(np.radians(th)))
        branch = flow.branch_of(state)
        f_br = flow.bridging(state.kn) if branch == surface_flow.BRANCH_MERGED else 0.0
        sim = Similarity()
        for i, t in enumerate(th):
            if t == 0:
                out.append((kn, 0, 0.0, 0.0, 0.0, "", "", "", p_w[0])); continue
            beta = 2 * xi[i] * due[i] / (rmw[i] * ue[i]**2 * r[i]**2)
            gw = rhoe[i] / (p_w[i] / (R_AIR * Tw))
            fpp = sim.fpp0(beta, gw)
            tau_c = rmw[i] * ue[i]**2 * r[i] * fpp / math.sqrt(2 * xi[i])
            tau_slip = tau_c / (1.0 + flow.slip_C * kn_local[i])
            if branch == surface_flow.BRANCH_MERGED:
                tau = (1 - f_br) * tau_slip + f_br * tau_fm[i]
            else:
                tau = tau_fm[i] if kn_local[i] >= surface_flow.KN_LOCAL_SLIP else tau_slip
            out.append((kn, int(t), tau_c, tau_c / q, tau / q, beta, gw, fpp, p_w[i]))
        for t in (2, 10, 20, 30, 45, 60, 75, 90):
            row = out[-91 + t]
            print(f"Kn {kn:5s} {t:2d} deg: beta {row[5]:.3f} g_w {row[6]:.3f} f''(0) {row[7]:.3f} "
                  f"C_tau continuum {row[3]:.4f} pipeline {row[4]:.4f}   (bridge f {f_br:.3f})")
    with open(HERE / ("similarity_hb_pressure_frozen.csv" if FROZEN else "similarity_hb_pressure.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case_kn", "theta_deg", "tau_continuum_Pa", "Ctau_continuum", "Ctau", "beta", "g_w", "fpp0", "p_w_Pa"])
        w.writerows(out)


if __name__ == "__main__":
    main()
