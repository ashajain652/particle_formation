"""Similarity-based shear round the WHOLE sphere (0-180 deg) on the design's decided pressure, at Holman & Boyd's four
conditions, for comparison with their DSMC and with Daniels & Schwartzentruber on the lee.

Pressure (design 2026-11-12-shoulder-melt-loss, A4): the model's windward wall_pressure (modified Newtonian +
Prandtl-Meyer, gamma_pm 1.15, the default) to 90 deg; beyond 90 deg the same Prandtl-Meyer expansion continued (no
p_inf floor, no PM_MAX_DEG cap) until it falls to Daniels & Schwartzentruber's lee pressure, then theirs (ds_model.csv).
Gas: unreacted (chemically frozen) post-shock state, as decided for the shear (A3): the model's Cantera construction
with gas.T_EQUILIBRATE raised in this process only. Shear: Lees local similarity (similarity_shear.Similarity), cold
wall, variable rho*mu; once the rising lee pressure drives f''(0) to zero (or the solution can no longer be continued)
the boundary layer is taken as separated and the shear is zero from there on. Then the model's slip factor
1/(1 + Kn_local); no free-molecular blend (above Kn 0.05 the design uses Daniels & Schwartzentruber instead).

Run: /Users/ashajain/miniforge3/envs/drama_env/bin/python -I similarity_lee.py   -> similarity_lee.csv
"""
import csv, math, sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(Path(__file__).resolve().parents[3]))  # Reentry Model/
from similarity_shear import Similarity, R, Tw, U, R_AIR  # noqa: E402
from reentry_model.shared import gas, surface_flow  # noqa: E402

gas.T_EQUILIBRATE = 1e12                     # unreacted post-shock gas (design A3)
T_INF, GAMMA_PM = 200.0, surface_flow.GAMMA_PM
CASES = (("0.002", 1.007e-4), ("0.01", 2.014e-5), ("0.05", 3.987e-6), ("0.25", 8.057e-7))
FPP_SEPARATED = 0.005


def joined_pressure(th_deg, p_stag, p_inf, q_inf, p_ds):
    th = np.radians(th_deg)
    p_wind, phi_star = surface_flow.wall_pressure(np.minimum(th, 0.5 * math.pi), p_stag, p_inf, q_inf, GAMMA_PM)
    turn = np.maximum(th - phi_star, 0.0)
    mach = surface_flow.mach_from_turn(turn, GAMMA_PM)
    p_pm = p_stag * (1.0 + 0.5 * (GAMMA_PM - 1.0) * mach**2) ** (-GAMMA_PM / (GAMMA_PM - 1.0))
    p = np.where(th_deg <= 90, p_wind, p_pm)
    lee = th_deg > 90
    meet = np.argmax(lee & (p_pm <= p_ds)) if np.any(lee & (p_pm <= p_ds)) else len(th_deg)
    p[meet:] = p_ds[meet:]
    return p, float(th_deg[meet]) if meet < len(th_deg) else float("nan")


def main():
    ds = {}
    for r in csv.DictReader(open(HERE / "ds_model.csv")):
        ds.setdefault(r["case_kn"], {})[int(float(r["theta_deg"]))] = float(r["Cp"])
    flow = surface_flow.SurfaceFlow(rarefied_shear="slip", gamma_pm=GAMMA_PM, kn_body_shock=0.01)
    mu_wall = gas.wilke_viscosity(Tw, {"N2": 0.79, "O2": 0.21}, flow.air.molar_masses)
    th = np.arange(0, 181)
    s = R * np.radians(th); r = R * np.sin(np.radians(th))
    out = []
    for kn, rho in CASES:
        q, p_inf = 0.5 * rho * U**2, rho * R_AIR * T_INF
        stag = flow.air.stagnation(rho, T_INF, U)
        dskey = next(k for k in ds if float(k) == float(kn))
        p_ds = np.array([p_inf + ds[dskey][t] * q for t in th])
        p_w, meet = joined_pressure(th.astype(float), stag.p, p_inf, q, p_ds)
        ue = np.zeros(len(th)); rhoe = np.zeros(len(th))
        for i, p in enumerate(p_w):
            e = flow.air.expand(stag, float(p))
            ue[i] = math.sqrt(max(0.0, 2.0 * (stag.h - e.h))); rhoe[i] = e.rho
        ue[0] = 0.0
        rmw = p_w / (R_AIR * Tw) * mu_wall
        fine = np.linspace(0, s[-1], 20 * 180 + 1)
        xi = np.interp(s, fine, cumulative_trapezoid(np.interp(fine, s, rmw) * np.interp(fine, s, ue)
                                                     * (R * np.sin(fine / R))**2, fine, initial=0.0))
        due = np.gradient(ue, s)
        lam_w, kn_local = surface_flow.wall_knudsen(p_w, np.full(len(th), Tw), np.full(len(th), mu_wall), R)
        sim, separated, sep_at, beta_prev = Similarity(), False, float("nan"), 0.5
        for i, t in enumerate(th):
            if t == 0 or t == 180:
                out.append((kn, int(t), p_w[i] / p_inf, 0.0, "", "", kn_local[i], meet, sep_at)); continue
            beta = 2 * xi[i] * due[i] / (rmw[i] * ue[i]**2 * r[i]**2)
            gw = rhoe[i] / (p_w[i] / (R_AIR * Tw))
            fpp = 0.0
            if not separated:
                try:
                    fpp = sim.fpp0(beta, gw)
                except RuntimeError:
                    fpp = 0.0
                    for b in np.linspace(beta_prev, beta, 8)[1:]:          # continue in small steps toward beta
                        try:
                            fpp = sim.fpp0(b, gw)
                        except RuntimeError:
                            fpp = 0.0; break
                if fpp <= FPP_SEPARATED:
                    separated, sep_at, fpp = True, float(t), 0.0
                beta_prev = beta
            tau_c = rmw[i] * ue[i]**2 * r[i] * fpp / math.sqrt(2 * xi[i]) if fpp > 0 else 0.0
            tau = tau_c / (1.0 + flow.slip_C * kn_local[i])
            out.append((kn, int(t), p_w[i] / p_inf, tau / q, beta, fpp, kn_local[i], meet, sep_at))
        rows = [o for o in out if o[0] == kn]
        print(f"Kn {kn}: PM meets D&S at {meet:.0f} deg; boundary layer separates at {sep_at} deg; "
              f"Kn_local > 0.1 from {next((o[1] for o in rows if o[6] > 0.1), 'never')} deg")
        print("   theta  p/p_inf  C_tau   beta    f''(0)")
        for t in (90, 100, 110, 120, 130, 140, 150, 160):
            o = rows[t]
            print(f"   {t:4d}  {o[2]:7.3f}  {o[3]:.4f}  {o[4] if o[4]=='' else round(o[4],3)!s:>7}  {o[5] if o[5]=='' else round(o[5],3)}")
    with open(HERE / "similarity_lee.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case_kn", "theta_deg", "p_over_pinf", "Ctau", "beta", "fpp0", "kn_local", "meet_deg", "separation_deg"])
        w.writerows(out)


if __name__ == "__main__":
    main()
