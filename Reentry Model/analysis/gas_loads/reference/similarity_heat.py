"""Similarity-based HEATING (and shear) round the whole sphere, at Holman & Boyd's four conditions, on the design's
decided pressure (model windward pressure; Prandtl-Meyer continued past 90 deg until it meets Daniels &
Schwartzentruber's lee pressure) with the unreacted (chemically frozen) post-shock gas decided for the shear.

Lees local similarity with the energy equation in total enthalpy and real frozen-air properties:
    (C f'')' + f f'' + beta (T/T_e - f'^2) = 0
    (C g'/Pr)' + f g' + (u_e^2/H_e) [ (1 - 1/Pr) C f' f'' ]' = 0,     g = H/H_e
    h = H_e g - u_e^2 f'^2 / 2,  T = T(h) for frozen N2/O2 (Cantera NASA-9, vibration in equilibrium),
    C = rho mu / (rho_w mu_w) = (T_w / T) mu(T)/mu(T_w)  (Wilke-Blottner viscosity, the model's),  Pr = 0.71,
    f(0) = f'(0) = 0, f'(inf) = 1;  g(0) = h_w/H_e, g(inf) = 1;
    tau_w = rho_w mu_w u_e^2 r f''(0) / sqrt(2 xi),   q_w = rho_w mu_w u_e r H_e g'(0) / (Pr sqrt(2 xi)).
For frozen gas rho_e/rho = T/T_e exactly. Beyond separation (f''(0) -> 0 on the rising lee pressure) nothing is
computed. Rarefaction: the model's own Matting bridge (n = 1, accommodation 0.8) applied as the model applies it,
q = q_stag [(1 - w) q_sim(theta)/q_sim(0) + w max(cos theta, 0)], q_stag = q_c (1 - exp(-q_fm/q_c)), w = 1 - q_stag/q_c,
with q_c the similarity stagnation flux in place of Fay-Riddell.

Run: /Users/ashajain/miniforge3/envs/drama_env/bin/python -I similarity_heat.py   -> similarity_heat.csv
"""
import csv, math, sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp, cumulative_trapezoid

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(Path(__file__).resolve().parents[3]))  # Reentry Model/
from similarity_shear import R, Tw, U, R_AIR  # noqa: E402
from similarity_lee import joined_pressure, CASES, T_INF, GAMMA_PM  # noqa: E402  (sets gas.T_EQUILIBRATE: frozen)
from reentry_model.shared import gas, surface_flow  # noqa: E402

PR = 0.71
ETA = np.linspace(0.0, 12.0, 500)
AIR_X = {"N2": 0.79, "O2": 0.21}


class FrozenTables:
    """Smooth (cubic-spline) frozen-air h(T), T(h), mu(T) and their derivatives: the collocation solver needs a
    continuous right-hand side, and piecewise-linear tables make it refine the mesh without end."""

    def __init__(self, air):
        from scipy.interpolate import CubicSpline
        g = air.gas
        self.T = np.arange(150.0, 20000.0 + 1.0, 50.0)
        h = []
        for t in self.T:
            g.TPX = t, 101325.0, "N2:0.79, O2:0.21"
            h.append(g.enthalpy_mass)
        self.h = np.array(h)
        self.mu = np.array([gas.wilke_viscosity(float(t), AIR_X, air.molar_masses) for t in self.T])
        self._h = CubicSpline(self.T, self.h)
        self._T = CubicSpline(self.h, self.T)
        self._mu = CubicSpline(self.T, self.mu)
        self._cp, self._dmu = self._h.derivative(), self._mu.derivative()

    def T_of_h(self, h):
        return self._T(np.clip(h, self.h[0], self.h[-1]))

    def h_of_T(self, T):
        return float(self._h(T))

    def mu_of_T(self, T):
        return self._mu(T)

    def dmu_of_T(self, T):
        return self._dmu(T)

    def cp_of_T(self, T):
        return self._cp(T)


class SimilarityHeat:
    def __init__(self, tab):
        self.tab, self.sol = tab, None

    def solve(self, beta, He, ue, Tw_):
        """f''(0), g'(0), g_w. The first call eases in from a hot wall (homotopy in the wall temperature), because
        a 1500 K wall under a 16,000 K edge is too far from any simple starting guess."""
        if self.sol is None:
            Te = float(self.tab.T_of_h(He - 0.5 * ue * ue))
            for t in np.geomspace(0.8 * Te, Tw_, 12)[:-1]:
                self._solve(beta, He, ue, float(t))
        return self._solve(beta, He, ue, Tw_)

    def _solve(self, beta, He, ue, Tw_):
        tab = self.tab
        hw = tab.h_of_T(Tw_)
        gw = hw / He
        Te = float(tab.T_of_h(He - 0.5 * ue * ue))
        muw = float(tab.mu_of_T(Tw_))
        D = (1.0 - 1.0 / PR) * ue * ue / He

        def props(y):
            f, fp, fpp, g, gp = y
            h = He * g - 0.5 * ue * ue * fp * fp
            T = np.clip(tab.T_of_h(h), 150.0, 20000.0)
            C = (Tw_ / T) * tab.mu_of_T(T) / muw
            dCdT = (-(Tw_ / T**2) * tab.mu_of_T(T) + (Tw_ / T) * tab.dmu_of_T(T)) / muw
            dh = He * gp - ue * ue * fp * fpp
            Cp = dCdT * dh / tab.cp_of_T(T)
            return T, C, Cp

        def rhs(eta, y):
            f, fp, fpp, g, gp = y
            T, C, Cp = props(y)
            fppp = -(Cp * fpp + f * fpp + beta * (T / Te - fp**2)) / C
            gpp = -(PR * (f * gp + D * (Cp * fp * fpp + C * fpp**2 + C * fp * fppp)) + Cp * gp) / C
            return np.vstack([fp, fpp, fppp, gp, gpp])

        bc = lambda a, b: np.array([a[0], a[1], b[1] - 1.0, a[3] - gw, b[3] - 1.0])
        if self.sol is None:
            e = np.exp(-ETA)
            guess = (ETA, np.vstack([ETA - 1 + e, 1 - e, e, 1 - (1 - gw) * e, (1 - gw) * e]))
        else:
            guess = (self.sol.x, self.sol.y)
        s = solve_bvp(rhs, bc, *guess, tol=1e-4, max_nodes=200000)
        if not s.success:
            raise RuntimeError(s.message)
        self.sol = s
        return float(s.y[2, 0]), float(s.y[4, 0]), gw


def main():
    flow = surface_flow.SurfaceFlow(rarefied_shear="slip", gamma_pm=GAMMA_PM, kn_body_shock=0.01)
    tab = FrozenTables(flow.air)
    ds = {}
    for r in csv.DictReader(open(HERE / "ds_model.csv")):
        ds.setdefault(r["case_kn"], {})[int(float(r["theta_deg"]))] = float(r["Cp"])
    th = np.arange(0, 181)
    s = R * np.radians(th); r = R * np.sin(np.radians(th))
    mu_wall = float(tab.mu_of_T(Tw))
    out = []
    for kn, rho in CASES:
        q_inf3, p_inf = 0.5 * rho * U**3, rho * R_AIR * T_INF
        stag = flow.air.stagnation(rho, T_INF, U)
        dskey = next(k for k in ds if float(k) == float(kn))
        p_ds = np.array([p_inf + ds[dskey][t] * 0.5 * rho * U**2 for t in th])
        p_w, meet = joined_pressure(th.astype(float), stag.p, p_inf, 0.5 * rho * U**2, p_ds)
        He = stag.h
        ue = np.array([math.sqrt(max(0.0, 2.0 * (He - flow.air.expand(stag, float(p)).h))) for p in p_w]); ue[0] = 0.0
        rmw = p_w / (R_AIR * Tw) * mu_wall
        fine = np.linspace(0, s[-1], 20 * 180 + 1)
        xi = np.interp(s, fine, cumulative_trapezoid(np.interp(fine, s, rmw) * np.interp(fine, s, ue)
                                                     * (R * np.sin(fine / R))**2, fine, initial=0.0))
        due = np.gradient(ue, s)
        sim, separated = SimilarityHeat(tab), False
        q_sim = np.full(len(th), np.nan); tau_sim = np.full(len(th), np.nan); sep_at = float("nan")
        # stagnation point: beta = 1/2, u_e -> 0, q0 = sqrt(2 rho_w mu_w beta0) H_e g'(0) / Pr
        beta0 = due[1]
        fpp0, gp0, _ = sim.solve(0.5, He, 0.0, Tw)
        q_sim[0] = math.sqrt(2 * rmw[0] * beta0) * He * gp0 / PR
        tau_sim[0] = 0.0
        for i in range(1, 180):
            if separated:
                continue
            beta = 2 * xi[i] * due[i] / (rmw[i] * ue[i]**2 * r[i]**2)
            try:
                fpp, gp, gw = sim.solve(beta, He, ue[i], Tw)
            except RuntimeError:
                fpp, gp = 0.0, 0.0
            if fpp <= 0.005:
                separated, sep_at = True, float(th[i]); continue
            fac = rmw[i] * ue[i] * r[i] / math.sqrt(2 * xi[i])
            tau_sim[i] = fac * ue[i] * fpp
            q_sim[i] = fac * He * gp / PR
        # the model's Matting bridge, with the similarity stagnation flux as q_c
        h_w = tab.h_of_T(Tw)
        q_fm = 0.8 * rho * U * (He - h_w)
        q_c = q_sim[0]
        q_stag = q_c * (1.0 - math.exp(-q_fm / q_c))
        w = 1.0 - q_stag / q_c
        q_bridged = q_stag * ((1 - w) * q_sim / q_c + w * np.maximum(np.cos(np.radians(th)), 0.0))
        print(f"Kn {kn}: q_c(similarity, frozen) {q_c:.3e} W/m2  (DSMC peak "
              f"{ {'0.002': 9.30e5, '0.01': 8.04e5, '0.05': 3.54e5, '0.25': 1.05e5}[kn]:.3e});  bridged {q_stag:.3e}, "
              f"w {w:.3f};  meets D&S at {meet:.0f}; separates at {sep_at}")
        for t in (0, 30, 60, 90, 100, 120, 140, 150):
            print(f"   {t:3d} deg  Cq continuum {q_sim[t] / q_inf3:.5f}  bridged {q_bridged[t] / q_inf3:.5f}")
        for i, t in enumerate(th):
            out.append((kn, int(t), q_sim[i] / q_inf3 if np.isfinite(q_sim[i]) else "",
                        q_bridged[i] / q_inf3 if np.isfinite(q_bridged[i]) else "",
                        tau_sim[i] / (0.5 * rho * U**2) if np.isfinite(tau_sim[i]) else "", meet, sep_at))
    with open(HERE / "similarity_heat.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["case_kn", "theta_deg", "Cq_continuum", "Cq", "Ctau_continuum", "meet_deg", "separation_deg"])
        wr.writerows(out)


if __name__ == "__main__":
    main()
