"""Daniels & Schwartzentruber (2026), arXiv:2609.28658, Appendix A:
analytical surface-pressure (Cp) and skin-friction (Cf) distributions on a sphere.

Equation numbers refer to Appendix A of the preprint.  Grouping of every
equation was checked against 200-400 dpi renders of pages 23-26.

Normalisation (paper eqs 1.1, 1.2):  Cf = tau_w / (0.5 rho U^2),
Cp = (p_w - p_inf) / (0.5 rho U^2)  -- identical to Holman & Boyd (2011)
C_tau and C_p, so no conversion is needed.

Final model (A 26, A 27) is a single composite formula valid for all Kn:
    Cf = Cf_c * psi_tau / (1 + Br^eta) + Cf_fm * Br^eta / (1 + Br^eta)
The slip (beta), transition (phi, W_r^T) and free-molecular (Br) pieces each
tend to unity / zero outside their regime.  regime="piecewise" is an optional
strict-switch variant (by Kn_inf: <0.01 continuum, <0.1 slip, <10 transition,
>=10 free molecular) kept only to measure how much the composite differs.
"""
import math
import numpy as np

SQPI = math.sqrt(math.pi)

# Table 1: incompressible angular positions (A 29a): k0, k1, k2, theta_Stokes
TAB1 = {"max": (49.5, 0.625, 0.50, 90.0),
        "min": (72.0, 0.25, 0.65, 180.0),
        "b":   (91.0, 0.0033, 0.85, 180.0)}
# Table 2: compressible angular positions (A 29b, c): k3..k8, theta_inf (None -> theta* + 5)
TAB2 = {"max": (0.012, 0.1291, 0.0394, 4.480, 2.700, 1.320, None),
        "s":   (0.012, 0.0625, 0.558, 1.200, 1.011, 0.072, 180.0),
        "min": (-0.0186, 0.2338, 0.476, 1.640, 1.660, 0.150, 90.0),
        "b":   (0.0107, 0.0403, 0.019, 1.280, 3.100, 0.700, 180.0)}

_erf = np.vectorize(math.erf)


def kn_from(Re, M, gamma):
    """(A 1): Kn = M/Re * sqrt(gamma*pi/2)  (lambda = mu/rho*sqrt(pi/(2RT)), over diameter)."""
    return M / Re * math.sqrt(gamma * math.pi / 2.0)


def post_shock(Re, M, gamma, omega):
    """(A 2), (A 3).  Bracket in (A 2) is the normal-shock T_s/T_inf."""
    TsT = 1.0 + 2.0 * (gamma - 1.0) / ((gamma + 1.0) ** 2 * M ** 2) * (M ** 2 - 1.0) * (1.0 + gamma * M ** 2)
    Res = Re * TsT ** (-omega)
    Ms = math.sqrt(((gamma - 1.0) * M ** 2 + 2.0) / (2.0 * gamma * M ** 2 - (gamma - 1.0)))
    return Res, Ms, TsT


def theta_ic_s(Re):
    """(A 30) Clift, Grace & Weber separation angle [deg]; no separation (180) for Re <= 20.
    Formula for 400 < Re <= 3e5 is used unchanged above 3e5 (outside paper's range)."""
    if Re <= 20.0:
        return 180.0
    if Re <= 400.0:
        return 180.0 - 42.5 * math.log(Re / 20.0) ** 0.483
    return 77.84 + 275.0 * Re ** -0.37


def angle(name, Re, M, gamma, clip=True):
    """Angular positions theta_max, theta_s, theta_min, theta_b [deg] (A 29a-c, A 30).
    clip=False drops the min(180, theta) of Appendix A (diagnostic only)."""
    if name == "s":
        tic = theta_ic_s(Re)
    else:
        k0, k1, k2, tstk = TAB1[name]
        tic = k0 + (tstk - k0) / (1.0 + k1 * Re) ** k2                      # (A 29a)
    if M < 0.3:
        th = tic
    else:
        k3, k4, k5, k6, k7, k8, tinf = TAB2[name]
        if M <= 1.0:
            th = tic * (1.0 + k3 * M ** 2 + k4 * M ** 4)                     # (A 29b)
        else:
            if tinf is None:                                                 # theta* + 5 deg (2.30)
                tinf = 34.0 + 40.0 * (gamma - 1.0) / (gamma + 1.0) + 5.0
            tM1 = tic * (1.0 + k3 + k4)
            th = tinf * (1.0 + k5 * (M - 1.0) / (1.0 + k6 * (M - 1.0) ** k7)) \
                * (tM1 / tinf) ** (M ** (-k8))                               # (A 29c)
    return min(180.0, th) if clip else th


def parameters(Re, M, Tw_T0, Tw_Tinf, gamma=1.4, omega=0.75, Pr=0.72,
               sig_t=1.0, sig_n=1.0, Kn=None, clip=True, tstar_exp=None):
    """All scalar model parameters.  Re, M, Kn are freestream values (Re on diameter).
    clip, tstar_exp: diagnostic switches only (defaults = paper as printed; tstar_exp
    overrides the omega/2 exponent on T*/T in (A 5))."""
    if Kn is None:
        Kn = kn_from(Re, M, gamma)
    P = dict(Re=Re, M=M, Kn=Kn, Tw_T0=Tw_T0, Tw_Tinf=Tw_Tinf, gamma=gamma, omega=omega, Pr=Pr,
             sig_t=sig_t, sig_n=sig_n)
    # "Re, M, Kn without subscript": freestream if M<1, post-shock if M>=1 (Appendix A preamble)
    if M >= 1.0:
        Res, Ms, TsT = post_shock(Re, M, gamma, omega)
        Rl, Ml = Res, Ms
        Knl = kn_from(Res, Ms, gamma)                                        # (A 1) with post-shock values
        P.update(Re_s=Res, M_s=Ms, Ts_Tinf=TsT, Kn_s=Knl)
    else:
        Rl, Ml, Knl = Re, M, Kn
        P.update(Re_s=float("nan"), M_s=float("nan"), Ts_Tinf=float("nan"), Kn_s=float("nan"))
    rPr = math.sqrt(Pr)
    TstarT = 0.285 - 0.126 * rPr + (0.715 * Tw_T0 + 0.126 * rPr) * (1.0 + (gamma - 1.0) / 2.0 * Ml ** 2)  # (A 6)
    if M < 1.0:
        FM = 1.0 + 0.143 * M ** 2 - 0.243 * M ** 4                           # (A 7a)
    else:
        FM = 0.56 + 0.351 * Ml ** 2 - 0.011 * Ml ** 4                        # (A 7b)
    Cfmax = 6.0 / Rl * (1.0 + Rl / 2.493) ** 0.5 * TstarT ** (omega / 2.0 if tstar_exp is None else tstar_exp) * FM  # (A 5)
    g = gamma
    if M < 1.0:
        Cp0M = 2.0 / (g * M ** 2) * ((1.0 + (g - 1.0) / 2.0 * M ** 2) ** (g / (g - 1.0)) - 1.0)      # (A 10a)
    else:
        Cp0M = 2.0 / (g * M ** 2) * ((2.0 * g * M ** 2 - (g - 1.0)) / (g + 1.0)
                                     * (1.0 + (g - 1.0) / 2.0 * Ml ** 2) ** (g / (g - 1.0)) - 1.0)  # (A 10b)
    Cp0 = Cp0M + 12.0 / (Rl + 0.644 * math.sqrt(Rl))                         # (A 9)
    m2, m4 = min(1.0, M ** 2), min(1.0, M ** 4)
    Cpmin = -6.0 / Re * (1.0 + Re ** 2.5 / 270.0) ** 0.384 \
        * (1.0 + 0.174 * m2 + 0.045 * m4) / max(1.0, M ** 2)                 # (A 11)
    dCpr = (Re / 20.0) ** 1.586 / (1.0 + 168.5 * (Re / 270.0) ** 1.72) \
        * (1.0 + 0.23 * m2 - 0.452 * m4) / max(1.0, M ** 2)                   # (A 12)
    Cpb = Cpmin + dCpr                                                        # (A 13)
    th = {k: angle(k, Re, M, gamma, clip) for k in ("max", "s", "min", "b")}
    # slip (A 14-16): Kn, Re without subscript -> post-shock in supersonic flow
    cs = (2.0 - sig_t) / sig_t
    beta_t = 1.0 / (1.0 + 6.0 * cs * Knl * (1.0 + 0.12 * Rl) ** 0.5)
    beta_p = (1.0 + 12.0 * cs * Knl) / (1.0 + 6.0 * cs * Knl)
    # transition (A 17-21): freestream Kn, M, Re
    phi_t = (1.0 + cs * sig_t * Kn) / (1.0 + Kn)
    phi_p = 1.0 / (1.0 + 3.0 * Kn / (1.0 - 0.1073 * sig_n))
    WrT = M ** (2.0 * omega) / Re * (1.0 + Tw_T0) ** omega
    psi_t = beta_t * phi_t / (1.0 + 1.27 * WrT)
    psi_p = beta_p * phi_p / (1.0 + 1.27 * WrT)
    a = M ** (2.0 * omega - 1.0)
    Br = WrT * (a + 1.0) / a                                                  # (A 28)
    eta = 1.8
    P.update(TstarT=TstarT, F_M=FM, Cf_max=Cfmax, Cp0M=Cp0M, Cp_0=Cp0, Cp_min=Cpmin, dCp_r=dCpr,
             Cp_b=Cpb, theta_max=th["max"], theta_s=th["s"], theta_min=th["min"], theta_b=th["b"],
             c_s=cs, beta_tau=beta_t, beta_p=beta_p, phi_tau=phi_t, phi_p=phi_p, WrT=WrT,
             psi_tau=psi_t, psi_p=psi_p, Br=Br, Br_eta=Br ** eta, w_fm=Br ** eta / (1.0 + Br ** eta))
    return P


def continuum(theta_deg, P):
    """(A 4) and (A 8)."""
    t = np.asarray(theta_deg, float)
    tmax, ts, tmin, tb = P["theta_max"], P["theta_s"], P["theta_min"], P["theta_b"]
    with np.errstate(divide="ignore", invalid="ignore"):
        cf = np.where(t <= tmax, P["Cf_max"] * np.sin(np.pi / 2 * t / tmax),
                      np.where(t <= ts, P["Cf_max"] * np.sin(np.pi / 2 * (t - ts) / (tmax - ts)) ** 2, 0.0))
        cp = np.where(t <= tmin, P["Cp_min"] + (P["Cp_0"] - P["Cp_min"]) * np.cos(np.pi / 2 * t / tmin) ** 2,
                      np.where(t <= tb, P["Cp_min"] + P["dCp_r"] * np.cos(np.pi / 2 * (t - tb) / (tmin - tb)) ** 2,
                               P["Cp_b"]))
    return cf, cp


def free_molecular(theta_deg, P):
    """(A 22)-(A 25), Schaaf & Chambre."""
    t = np.radians(np.asarray(theta_deg, float))
    s = P["M"] * math.sqrt(P["gamma"] / 2.0)
    xi = s * np.cos(t)
    e, E = np.exp(-xi ** 2), 1.0 + _erf(xi)
    st, sn, rT = P["sig_t"], P["sig_n"], math.sqrt(P["Tw_Tinf"])
    cf = st * np.sin(t) / (s * SQPI) * (e + xi * SQPI * E)
    cp = 1.0 / s ** 2 * (((2.0 - sn) / SQPI * xi + sn / 2.0 * rT) * e
                         + ((2.0 - sn) * (xi ** 2 + 0.5) + sn / 2.0 * xi * SQPI * rT) * E - 1.0)
    return cf, cp


def distribution(theta_deg, P, regime="composite"):
    """Cf, Cp at theta (deg from stagnation point).  regime='composite' is the paper's (A 26, A 27)."""
    cfc, cpc = continuum(theta_deg, P)
    cff, cpf = free_molecular(theta_deg, P)
    w = P["w_fm"]
    if regime == "continuum":            # no slip / rarefaction (used for the paper's no-slip CFD figures)
        return cfc, cpc
    if regime == "composite":
        return cfc * P["psi_tau"] / (1 + P["Br_eta"]) + cff * w, cpc * P["psi_p"] / (1 + P["Br_eta"]) + cpf * w
    Kn = P["Kn"]
    if Kn < 0.01:
        return cfc, cpc
    if Kn < 0.1:
        return cfc * P["beta_tau"], cpc * P["beta_p"]
    if Kn < 10.0:
        return cfc * P["psi_tau"] / (1 + P["Br_eta"]) + cff * w, cpc * P["psi_p"] / (1 + P["Br_eta"]) + cpf * w
    return cff, cpf


def drag(P, regime="composite", n=180001):
    """Cd = 2 * int_0^pi (Cp cos + Cf sin) sin dtheta  (reference area pi R^2)."""
    t = np.linspace(0.0, 180.0, n)
    cf, cp = distribution(t, P, regime)
    r = np.radians(t)
    return 2.0 * np.trapezoid((cp * np.cos(r) + cf * np.sin(r)) * np.sin(r), r)


def taw_T0(M, gamma, Pr):
    """(2.28) adiabatic wall temperature over T0, with post-shock Mach for M >= 1."""
    Ml = post_shock(1.0, M, gamma, 0.75)[1] if M >= 1.0 else M
    r = math.sqrt(Pr)
    return (1.0 + r * (gamma - 1.0) / 2.0 * Ml ** 2) / (1.0 + (gamma - 1.0) / 2.0 * Ml ** 2)


# ---------------------------------------------------------------------------
# Holman & Boyd (2011) conditions
# ---------------------------------------------------------------------------
R_AIR, G_AIR, OM_AIR, PR_AIR = 287.05, 1.4, 0.75, 0.72
T_INF, U_INF, D_SPH, T_W = 200.0, 7108.0, 0.3048, 1500.0
HB_CASES = [("0.002", 1.007e-4, 0.959), ("0.01", 2.014e-5, 0.977),
            ("0.05", 3.987e-6, 1.22), ("0.25", 8.057e-7, 1.53)]


def mu_sutherland(T):
    """Air viscosity [Pa s]; the paper gives only the power-law exponent, no reference value."""
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


def hb_parameters(rho, mu=None):
    mu = mu_sutherland(T_INF) if mu is None else mu
    a = math.sqrt(G_AIR * R_AIR * T_INF)
    M = U_INF / a
    Re = rho * U_INF * D_SPH / mu
    T0 = T_INF * (1.0 + (G_AIR - 1.0) / 2.0 * M ** 2)        # perfect-gas stagnation temperature
    return parameters(Re, M, T_W / T0, T_W / T_INF, G_AIR, OM_AIR, PR_AIR, 1.0, 1.0), T0


if __name__ == "__main__":
    import csv, os
    out = os.path.dirname(os.path.abspath(__file__))
    th = np.arange(0, 181, 1.0)
    rows, prow = [], []
    for lab, rho, cd_hb in HB_CASES:
        P, T0 = hb_parameters(rho)
        cf, cp = distribution(th, P)
        for t, a, b in zip(th, cp, cf):
            rows.append([lab, int(t), f"{a:.6g}", f"{b:.6g}", ""])
        P["case_kn"], P["T0"] = lab, T0
        P["Cd_model"], P["Cd_continuum_only"], P["Cd_HB_DSMC"] = drag(P), drag(P, "continuum"), cd_hb
        prow.append(P)
    with open(os.path.join(out, "ds_model.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case_kn", "theta_deg", "Cp", "Ctau", "Cq"])
        w.writerows(rows)
    keys = ["case_kn", "Re", "Kn", "M", "T0", "Tw_T0", "Re_s", "M_s", "Kn_s", "TstarT", "F_M", "Cp0M", "Cp_0",
            "Cp_min", "dCp_r", "Cp_b", "Cf_max", "theta_max", "theta_min", "theta_s", "theta_b",
            "beta_tau", "beta_p", "phi_tau", "phi_p", "WrT", "psi_tau", "psi_p", "Br", "w_fm",
            "Cd_model", "Cd_continuum_only", "Cd_HB_DSMC"]
    with open(os.path.join(out, "ds_params.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for P in prow:
            w.writerow([P[k] if isinstance(P[k], str) else f"{P[k]:.6g}" for k in keys])
    for P in prow:
        print(" ".join(f"{k}={P[k] if isinstance(P[k], str) else format(P[k], '.4g')}" for k in keys))
