"""Daniels & Schwartzentruber (2026), arXiv:2609.28658, Appendix A: analytical surface pressure and skin friction on
a sphere, from the continuum to the free-molecular limit, ported from analysis/gas_loads/reference/ds_model.py.

Cp = (p_w - p_inf) / (0.5 rho U^2) and Cf = tau_w / (0.5 rho U^2) at polar angle theta from the stagnation point
(radians at this module's interface). Equation numbers refer to Appendix A. Perfect gas, Sutherland viscosity at
T_inf, Knudsen number by (A 1), theta_s capped at 180 degrees, (2.29) as printed, and the composite rarefied formula
(A 26)-(A 27) at every Knudsen number (no regime switch).
"""
import math
from dataclasses import dataclass

import numpy as np

SQPI = math.sqrt(math.pi)

# Table 1: incompressible angular positions (A 29a): k0, k1, k2, theta_Stokes [deg]
_TAB1 = {"max": (49.5, 0.625, 0.50, 90.0),
         "min": (72.0, 0.25, 0.65, 180.0),
         "b": (91.0, 0.0033, 0.85, 180.0)}
# Table 2: compressible angular positions (A 29b, c): k3..k8, theta_inf (None -> theta* + 5 deg)
_TAB2 = {"max": (0.012, 0.1291, 0.0394, 4.480, 2.700, 1.320, None),
         "s": (0.012, 0.0625, 0.558, 1.200, 1.011, 0.072, 180.0),
         "min": (-0.0186, 0.2338, 0.476, 1.640, 1.660, 0.150, 90.0),
         "b": (0.0107, 0.0403, 0.019, 1.280, 3.100, 0.700, 180.0)}
_ETA = 1.8
_erf = np.vectorize(math.erf)


def mu_sutherland(T):
    """Air viscosity [Pa s]; the paper gives only the power-law exponent, no reference value."""
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


@dataclass(frozen=True)
class DSCondition:
    """Freestream similarity numbers (Reynolds on the diameter) and gas constants of one condition."""
    re_inf: float
    mach_inf: float
    kn_inf: float
    tw_t0: float      # wall temperature over the perfect-gas stagnation temperature
    tw_tinf: float    # wall temperature over the freestream temperature
    gamma: float = 1.4
    omega: float = 0.75
    pr: float = 0.72
    sigma_t: float = 1.0
    sigma_n: float = 1.0
    t0: float = float("nan")   # perfect-gas stagnation temperature [K]; reported as "T0" by parameters()


def _kn_from(Re, M, gamma):
    """(A 1): Kn = M/Re * sqrt(gamma*pi/2)."""
    return M / Re * math.sqrt(gamma * math.pi / 2.0)


def condition(rho_inf, T_inf, V, diameter, T_wall, gamma=1.4, omega=0.75, R=287.05):
    """DSCondition from freestream density [kg/m3], temperature [K], speed [m/s], body diameter [m] and wall
    temperature [K]; Re and Kn by the paper's definitions, T0 the perfect-gas stagnation temperature."""
    M = V / math.sqrt(gamma * R * T_inf)
    Re = rho_inf * V * diameter / mu_sutherland(T_inf)
    T0 = T_inf * (1.0 + (gamma - 1.0) / 2.0 * M ** 2)
    return DSCondition(re_inf=Re, mach_inf=M, kn_inf=_kn_from(Re, M, gamma), tw_t0=T_wall / T0,
                       tw_tinf=T_wall / T_inf, gamma=gamma, omega=omega, t0=T0)


def _post_shock(Re, M, gamma, omega):
    """(A 2), (A 3): post-shock Reynolds and Mach numbers and T_s/T_inf."""
    TsT = 1.0 + 2.0 * (gamma - 1.0) / ((gamma + 1.0) ** 2 * M ** 2) * (M ** 2 - 1.0) * (1.0 + gamma * M ** 2)
    Res = Re * TsT ** (-omega)
    Ms = math.sqrt(((gamma - 1.0) * M ** 2 + 2.0) / (2.0 * gamma * M ** 2 - (gamma - 1.0)))
    return Res, Ms, TsT


def _theta_ic_s(Re):
    """(A 30) Clift, Grace & Weber separation angle [deg]; no separation (180) for Re <= 20. The formula for
    400 < Re <= 3e5 is used unchanged above 3e5 (outside the paper's range)."""
    if Re <= 20.0:
        return 180.0
    if Re <= 400.0:
        return 180.0 - 42.5 * math.log(Re / 20.0) ** 0.483
    return 77.84 + 275.0 * Re ** -0.37


def _angle(name, Re, M, gamma):
    """Angular positions theta_max, theta_s, theta_min, theta_b [deg] (A 29a-c, A 30), capped at 180."""
    if name == "s":
        tic = _theta_ic_s(Re)
    else:
        k0, k1, k2, tstk = _TAB1[name]
        tic = k0 + (tstk - k0) / (1.0 + k1 * Re) ** k2                       # (A 29a)
    if M < 0.3:
        th = tic
    else:
        k3, k4, k5, k6, k7, k8, tinf = _TAB2[name]
        if M <= 1.0:
            th = tic * (1.0 + k3 * M ** 2 + k4 * M ** 4)                      # (A 29b)
        else:
            if tinf is None:                                                  # theta* + 5 deg (2.30)
                tinf = 34.0 + 40.0 * (gamma - 1.0) / (gamma + 1.0) + 5.0
            tM1 = tic * (1.0 + k3 + k4)
            th = tinf * (1.0 + k5 * (M - 1.0) / (1.0 + k6 * (M - 1.0) ** k7)) \
                * (tM1 / tinf) ** (M ** (-k8))                                # (A 29c)
    return min(180.0, th)


def parameters(cond):
    """All scalar model parameters of `cond`, keyed as the columns of tests/fixtures/gas_loads/ds_params.csv
    (angles in degrees, as in the paper)."""
    Re, M, Kn = cond.re_inf, cond.mach_inf, cond.kn_inf
    Tw_T0, Tw_Tinf = cond.tw_t0, cond.tw_tinf
    gamma, omega, Pr, sig_t, sig_n = cond.gamma, cond.omega, cond.pr, cond.sigma_t, cond.sigma_n
    P = dict(Re=Re, M=M, Kn=Kn, T0=cond.t0, Tw_T0=Tw_T0, Tw_Tinf=Tw_Tinf, gamma=gamma, omega=omega, Pr=Pr,
             sig_t=sig_t, sig_n=sig_n)
    # "Re, M, Kn without subscript": freestream if M < 1, post-shock if M >= 1 (Appendix A preamble)
    if M >= 1.0:
        Res, Ms, TsT = _post_shock(Re, M, gamma, omega)
        Rl, Ml = Res, Ms
        Knl = _kn_from(Res, Ms, gamma)
        P.update(Re_s=Res, M_s=Ms, Ts_Tinf=TsT, Kn_s=Knl)
    else:
        Rl, Ml, Knl = Re, M, Kn
        P.update(Re_s=float("nan"), M_s=float("nan"), Ts_Tinf=float("nan"), Kn_s=float("nan"))
    rPr = math.sqrt(Pr)
    TstarT = 0.285 - 0.126 * rPr + (0.715 * Tw_T0 + 0.126 * rPr) * (1.0 + (gamma - 1.0) / 2.0 * Ml ** 2)  # (A 6)
    if M < 1.0:
        FM = 1.0 + 0.143 * M ** 2 - 0.243 * M ** 4                            # (A 7a)
    else:
        FM = 0.56 + 0.351 * Ml ** 2 - 0.011 * Ml ** 4                         # (A 7b)
    Cfmax = 6.0 / Rl * (1.0 + Rl / 2.493) ** 0.5 * TstarT ** (omega / 2.0) * FM   # (A 5)
    g = gamma
    if M < 1.0:
        Cp0M = 2.0 / (g * M ** 2) * ((1.0 + (g - 1.0) / 2.0 * M ** 2) ** (g / (g - 1.0)) - 1.0)       # (A 10a)
    else:
        Cp0M = 2.0 / (g * M ** 2) * ((2.0 * g * M ** 2 - (g - 1.0)) / (g + 1.0)
                                     * (1.0 + (g - 1.0) / 2.0 * Ml ** 2) ** (g / (g - 1.0)) - 1.0)   # (A 10b)
    Cp0 = Cp0M + 12.0 / (Rl + 0.644 * math.sqrt(Rl))                          # (A 9)
    m2, m4 = min(1.0, M ** 2), min(1.0, M ** 4)
    Cpmin = -6.0 / Re * (1.0 + Re ** 2.5 / 270.0) ** 0.384 \
        * (1.0 + 0.174 * m2 + 0.045 * m4) / max(1.0, M ** 2)                  # (A 11)
    dCpr = (Re / 20.0) ** 1.586 / (1.0 + 168.5 * (Re / 270.0) ** 1.72) \
        * (1.0 + 0.23 * m2 - 0.452 * m4) / max(1.0, M ** 2)                    # (A 12)
    Cpb = Cpmin + dCpr                                                         # (A 13)
    th = {k: _angle(k, Re, M, gamma) for k in ("max", "s", "min", "b")}
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
    Br = WrT * (a + 1.0) / a                                                   # (A 28)
    P.update(TstarT=TstarT, F_M=FM, Cf_max=Cfmax, Cp0M=Cp0M, Cp_0=Cp0, Cp_min=Cpmin, dCp_r=dCpr,
             Cp_b=Cpb, theta_max=th["max"], theta_s=th["s"], theta_min=th["min"], theta_b=th["b"],
             c_s=cs, beta_tau=beta_t, beta_p=beta_p, phi_tau=phi_t, phi_p=phi_p, WrT=WrT,
             psi_tau=psi_t, psi_p=psi_p, Br=Br, Br_eta=Br ** _ETA, w_fm=Br ** _ETA / (1.0 + Br ** _ETA))
    return P


def _continuum(theta_deg, P):
    """(A 4) and (A 8): (Cf, Cp)."""
    t = np.asarray(theta_deg, float)
    tmax, ts, tmin, tb = P["theta_max"], P["theta_s"], P["theta_min"], P["theta_b"]
    with np.errstate(divide="ignore", invalid="ignore"):
        cf = np.where(t <= tmax, P["Cf_max"] * np.sin(np.pi / 2 * t / tmax),
                      np.where(t <= ts, P["Cf_max"] * np.sin(np.pi / 2 * (t - ts) / (tmax - ts)) ** 2, 0.0))
        cp = np.where(t <= tmin, P["Cp_min"] + (P["Cp_0"] - P["Cp_min"]) * np.cos(np.pi / 2 * t / tmin) ** 2,
                      np.where(t <= tb, P["Cp_min"] + P["dCp_r"] * np.cos(np.pi / 2 * (t - tb) / (tmin - tb)) ** 2,
                               P["Cp_b"]))
    return cf, cp


def _free_molecular(theta_deg, P):
    """(A 22)-(A 25), Schaaf & Chambre: (Cf, Cp)."""
    t = np.radians(np.asarray(theta_deg, float))
    s = P["M"] * math.sqrt(P["gamma"] / 2.0)
    xi = s * np.cos(t)
    e, E = np.exp(-xi ** 2), 1.0 + _erf(xi)
    st, sn, rT = P["sig_t"], P["sig_n"], math.sqrt(P["Tw_Tinf"])
    cf = st * np.sin(t) / (s * SQPI) * (e + xi * SQPI * E)
    cp = 1.0 / s ** 2 * (((2.0 - sn) / SQPI * xi + sn / 2.0 * rT) * e
                         + ((2.0 - sn) * (xi ** 2 + 0.5) + sn / 2.0 * xi * SQPI * rT) * E - 1.0)
    return cf, cp


def _composite(theta_deg, P):
    """(Cf, Cp) by the composite formula (A 26, A 27), valid at every Knudsen number."""
    cfc, cpc = _continuum(theta_deg, P)
    cff, cpf = _free_molecular(theta_deg, P)
    w = P["w_fm"]
    return (cfc * P["psi_tau"] / (1 + P["Br_eta"]) + cff * w,
            cpc * P["psi_p"] / (1 + P["Br_eta"]) + cpf * w)


def coefficients(theta, cond):
    """(Cp, Cf) at polar angles `theta` [rad] from the stagnation point: Cp = (p - p_inf)/(0.5 rho U^2) and
    Cf = tau/(0.5 rho U^2)."""
    cf, cp = _composite(np.degrees(np.asarray(theta, float)), parameters(cond))
    return cp, cf


def drag_coefficient(cond, n=180001):
    """Cd = 2 int_0^pi (Cp cos + Cf sin) sin dtheta (reference area pi R^2), composite regime."""
    t = np.linspace(0.0, math.pi, n)
    cp, cf = coefficients(t, cond)
    return float(2.0 * np.trapezoid((cp * np.cos(t) + cf * np.sin(t)) * np.sin(t), t))
