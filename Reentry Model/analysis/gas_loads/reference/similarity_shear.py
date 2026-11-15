"""Local-similarity wall shear (Lees 1956 form) on the re-entry model's own windward edge state, at Holman & Boyd's
four conditions. Compressible similar boundary layer in Lees-Dorodnitsyn variables with the wall-based xi:
    xi = int_0^s rho_w mu_w u_e r^2 ds,  beta = (2 xi / u_e) du_e/dxi,
    (C f'')' + f f'' + beta (g - f'^2) = 0,   (C g'/Pr)' + f g' = 0,
    f(0) = f'(0) = 0, f'(inf) = 1;  g(0) = g_w, g(inf) = 1;  C = rho mu / (rho_w mu_w) = (g/g_w)^(omega - 1),
    tau_w = rho_w mu_w u_e^2 r f''(0) / sqrt(2 xi).
g is the density ratio rho_e/rho (perfect-gas link to temperature, used for C); its wall value g_w = rho_e/rho_w is
taken from the model's real-gas edge density and the wall density p_w/(R T_w). Dissipation is neglected in the
energy equation. The continuum result then goes through the model's own slip and bridging, exactly as tau_c does.
"""
import csv, json, math
import numpy as np
from scipy.integrate import solve_bvp, cumulative_trapezoid

from pathlib import Path
SP = str(Path(__file__).resolve().parent) + "/"
R, Tw, U, R_AIR, PR, OMEGA = 0.1524, 1500.0, 7108.0, 287.05, 0.71, 0.7
MU_W = 1.458e-6 * Tw**1.5 / (Tw + 110.4)
ETA = np.linspace(0.0, 10.0, 400)


class Similarity:
    def __init__(self):
        self.sol = None

    def fpp0(self, beta, gw, omega=OMEGA, pr=PR):
        def rhs(eta, y):
            f, fp, fpp, g, gp = y
            g = np.maximum(g, 1e-6)
            C = (g / gw) ** (omega - 1.0)
            dC = (omega - 1.0) * C / g * gp
            fppp = -(f * fpp + beta * (g - fp**2) + dC * fpp) / C
            gpp = -(pr * f * gp / C) - dC / C * gp        # (C g'/Pr)' + f g' = 0
            return np.vstack([fp, fpp, fppp, gp, gpp])
        bc = lambda a, b: np.array([a[0], a[1], b[1] - 1.0, a[3] - gw, b[3] - 1.0])
        if self.sol is None:
            y0 = np.vstack([ETA - 1 + np.exp(-ETA), 1 - np.exp(-ETA), np.exp(-ETA),
                            1 - (1 - gw) * np.exp(-ETA), (1 - gw) * np.exp(-ETA)])
            guess = (ETA, y0)
        else:
            guess = (self.sol.x, self.sol.y)
        s = solve_bvp(rhs, bc, *guess, tol=1e-6, max_nodes=20000)
        if not s.success:
            raise RuntimeError(f"BVP failed beta {beta} gw {gw}: {s.message}")
        self.sol = s
        return float(s.y[2, 0])


def selftest():
    sim = Similarity()
    out = {b: sim.fpp0(b, 1.0, omega=1.0) for b in (0.5, 0.0, 1.0)}   # constant-property Falkner-Skan
    print("Falkner-Skan check f''(0):", {k: round(v, 4) for k, v in out.items()}, "expected 0.9277, 0.4696, 1.2326")


def main():
    selftest()
    diag = list(csv.DictReader(open(SP + "model_loads_diagnostics.csv")))
    summ = json.load(open(SP + "model_loads_summary.json"))
    hb = {(r["case_kn"], r["theta_deg"]): r for r in csv.DictReader(open(SP + "hb_dsmc.csv"))}
    rows_out = []
    for kn, rho_inf in (("0.002", 1.007e-4), ("0.01", 2.014e-5), ("0.05", 3.987e-6), ("0.25", 8.057e-7)):
        q = 0.5 * rho_inf * U**2
        d = {int(float(r["theta_deg"])): r for r in diag if r["case_kn"] == kn}
        th = np.arange(0, 91)
        ue = np.array([float(d[t]["u_e_ms"]) for t in th]); ue[0] = 0.0
        pw = np.array([float(d[t]["p_w_Pa"]) for t in th])
        rhoe = np.array([float(d[t]["rho_e_kgm3"]) for t in th])
        knl = np.array([float(d[t]["kn_local"]) for t in th])
        taufm = np.array([float(d[t]["tau_fm_Pa"]) for t in th])
        s = R * np.radians(th); r = R * np.sin(np.radians(th))
        rmw = pw / (R_AIR * Tw) * MU_W
        # xi on a refined grid (u_e r^2 ~ s^3 near the nose)
        fine = np.linspace(0, s[-1], 20 * 90 + 1)
        integrand = np.interp(fine, s, rmw) * np.interp(fine, s, ue) * (R * np.sin(fine / R))**2
        xi = np.interp(s, fine, cumulative_trapezoid(integrand, fine, initial=0.0))
        due = np.gradient(ue, s)
        sim = Similarity()
        info = summ[kn]; branch = info["body_branch"]; f_br = info["shear_bridging_f_sesam_table"]
        for i, t in enumerate(th):
            if t == 0:
                rows_out.append((kn, 0, 0.0, 0.0, 0.0, "", "", "")); continue
            beta = 2 * xi[i] * due[i] / (rmw[i] * ue[i]**2 * r[i]**2)
            gw = rhoe[i] / (pw[i] / (R_AIR * Tw))
            fpp = sim.fpp0(beta, gw)
            tau_c = rmw[i] * ue[i]**2 * r[i] * fpp / math.sqrt(2 * xi[i])
            tau_slip = tau_c / (1.0 + 1.0 * knl[i])                  # slip_C = 1 (sigma_v = 1), as in the model
            if branch == "merged":
                tau = (1 - f_br) * tau_slip + f_br * taufm[i]
            else:
                tau = taufm[i] if knl[i] >= 0.1 else tau_slip
            rows_out.append((kn, int(t), tau_c, tau_c / q, tau / q, beta, gw, fpp))
        for t in (2, 10, 20, 30, 45, 60, 75, 90):
            row = rows_out[-91 + t]
            print(f"Kn {kn} {t:2d} deg: beta {row[5]:.3f} g_w {row[6]:.3f} f''(0) {row[7]:.3f}  "
                  f"C_tau continuum {row[3]:.4f} pipeline {row[4]:.4f}  DSMC {float(hb[(kn, str(t))]['Ctau']):.4f}")
    with open(SP + "similarity_shear.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["case_kn", "theta_deg", "tau_continuum_Pa", "Ctau_continuum", "Ctau", "beta", "g_w", "fpp0"])
        w.writerows(rows_out)


if __name__ == "__main__":
    main()
