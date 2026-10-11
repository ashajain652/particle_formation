"""Girin's (2017) published spraying cases run with his own simplifications ("Girin-as-published" mode, spec Step 3
section 13.2): a sphere of radius R0 in a constant gas density, potential-flow edge velocity V_a = 1.5 (V_inf - w)
sin phi, Ranger's boundary layer delta_a = 2.2 R Re_a^-1/2 Psi(phi) with Re_a = Re_inf R~ (1 - W), the conjugated-
layer relations of his Eq. (2), the dispersion table for Delta_f and Im Omega_f, no melt limit and no runoff, spherical
belts of width R dphi, the induction rule tau_ind = int dt / t_per >= 1 releasing R dphi / lambda_f tori of
cross-section radius r = k_r lambda_f, and the deceleration law W = 1 - exp(-C tau), C = 2 alpha^1/2, tau = t / t_ch,
t_ch = 2 R0 / (alpha^1/2 V_inf).

Densities (measured 2026-09-20 against his Table 1): We_inf uses the ambient density rho_inf = rho_a / 6 and Re_inf the
compressed density rho_a = 1e-7 g/cm3 -- that reproduces his GI = We_inf Re_inf^-1/2 (13.0 / 3.55 / 43.5) and his
t_ch (alpha = rho_inf / rho_m); the Reynolds number entering delta_a is a switch (`re_density`): "shock" (rho_a, as the
printed Re_inf = 529) or "ambient" (rho_inf); with "ambient" his N and r_med are reproduced within 15 % for the iron
variants, with "shock" the droplets come out 2.5 x smaller (the paper does not say which density its delta_a used).
The stripping rate per area is independent of delta_a (it scales with V_s only), so t_s.d. is the same either way:
0.5 x his printed value for the iron variants and far below it for the stony one, whose lambda_f (3.5 mm) exceeds the
body radius -- his belt discretisation and induction handling are unstated (see the Step 3 plan's measured facts)."""
from dataclasses import dataclass

import numpy as np

from . import dispersion
from .surface_flow import ranger_psi

RHO_A, MU_AIR = 1.0e-4, 6.8e-5          # kg/m3, Pa s (Girin 2017 section 5)
COMPRESSION = 6.0
RE_DENSITY_NAMES = ("shock", "ambient")


@dataclass
class GirinVariant:
    name: str
    V: float          # m/s
    rho_m: float      # kg/m3
    mu_m: float       # Pa s
    sigma: float      # N/m
    R0: float = 3.0e-3


def dimensionless_numbers(v, rho_a=RHO_A, mu_a=MU_AIR, compression=COMPRESSION):
    """(Re_inf, We_inf, GI, alpha, mu, p, t_ch) with Girin's density conventions."""
    rho_inf = rho_a / compression
    Re = rho_a * 2.0 * v.R0 * v.V / mu_a
    We = rho_inf * 2.0 * v.R0 * v.V ** 2 / v.sigma
    alpha, mu = rho_inf / v.rho_m, mu_a / v.mu_m
    p = 1.0 + (alpha * mu) ** (1.0 / 3.0)
    t_ch = 2.0 * v.R0 / (np.sqrt(alpha) * v.V)
    return Re, We, We / np.sqrt(Re), alpha, mu, p, t_ch


def critical_angle(GI, p, we_cr=dispersion.WE_CRITICAL_PRACTICAL):
    """phi_cr [rad] at tau = 0 from Girin's Eq. (3): 2.475 p^-2 sin^2 phi Psi(phi) GI = We_cr."""
    from scipy.optimize import brentq
    f = lambda phi: 2.475 / p ** 2 * np.sin(phi) ** 2 * ranger_psi(phi) * GI - we_cr
    if f(np.pi / 2 - 1e-9) <= 0.0:
        return np.pi / 2
    return brentq(f, 1e-4, np.pi / 2 - 1e-9)


def surface_state(v, phi, Re_a, R, W, alpha, mu, p, table):
    """Per-belt delta_a, delta_m, V_s, We_s, Delta_f, Im Omega_f, lambda_f, r, t_per, t_f at the current radius R."""
    delta_a = 2.2 * R * ranger_psi(phi) / np.sqrt(Re_a)
    delta_m = (alpha / mu ** 2) ** (1.0 / 3.0) * delta_a
    V_a = 1.5 * v.V * (1.0 - W) * np.sin(phi)
    V_s = (alpha * mu) ** (1.0 / 3.0) / p * V_a
    we_s = v.rho_m * V_s ** 2 * delta_m / v.sigma
    delta_f, im_f, _ = table(np.maximum(we_s, table.we_onset))
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = 2.0 * np.pi * delta_m / delta_f
        t_f = delta_m / (V_s * im_f)
    return delta_a, delta_m, V_s, we_s, delta_f, im_f, lam, t_f


def run_variant(v, k_r=0.17, k_t=1.1, we_cr=dispersion.WE_CRITICAL_PRACTICAL, re_density="shock", dphi_deg=0.5, dt_fraction=0.2,
                mass_fraction_end=1e-4, table=None):
    """Girin's spraying flight of one variant. Returns a dict with the exact-tier numbers (GI, phi_cr, t_f at tau = 0,
    r at tau = 0 and 90 degrees, t_ch) and the integrated ones (N, r_med, r_min, r_max, t_sd, mass-loss history and the
    droplet list)."""
    if re_density not in RE_DENSITY_NAMES:
        raise ValueError("re_density must be one of {}".format(RE_DENSITY_NAMES))
    table = table or dispersion.DispersionTable()
    Re, We, GI, alpha, mu, p, t_ch = dimensionless_numbers(v)
    Re_used = Re if re_density == "shock" else Re / COMPRESSION
    C = 2.0 * np.sqrt(alpha)
    phi = np.radians(np.arange(dphi_deg / 2.0, 90.0, dphi_deg))
    dphi = np.radians(dphi_deg)
    m0 = v.rho_m * 4.0 / 3.0 * np.pi * v.R0 ** 3
    # tau = 0 exact-tier quantities
    _, _, _, we0, _, _, lam0, tf0 = surface_state(v, phi, Re_used, v.R0, 0.0, alpha, mu, p, table)
    unstable0 = we0 > we_cr
    out = {"Re_inf": Re, "We_inf": We, "GI": GI, "alpha": alpha, "mu": mu, "p": p, "t_ch_s": t_ch, "Re_used": Re_used,
           "phi_cr_deg_eq3": float(np.degrees(critical_angle(GI, p, we_cr))),
           "phi_cr_deg_direct": float(np.degrees(phi[unstable0][0])) if unstable0.any() else None,
           "t_f_min_us_tau0": float(tf0[unstable0].min() * 1e6) if unstable0.any() else None,
           "r_um_90deg_tau0": float(k_r * lam0[-1] * 1e6), "we_s_90deg_tau0": float(we0[-1])}
    # the flight
    m, t = m0, 0.0
    tind = np.zeros_like(phi)
    r_rel, n_rel, t_rel, m_hist = [], [], [], [(0.0, m0)]
    while m > mass_fraction_end * m0:
        R = (3.0 * m / (4.0 * np.pi * v.rho_m)) ** (1.0 / 3.0)
        W = 1.0 - np.exp(-C * t / t_ch)
        Re_a = Re_used * (R / v.R0) * (1.0 - W)
        _, delta_m, V_s, we_s, delta_f, im_f, lam, t_f = surface_state(v, phi, Re_a, R, W, alpha, mu, p, table)
        active = we_s > we_cr
        if not active.any():
            break
        t_per = k_t * t_f
        dt = dt_fraction * float(t_per[active].min())
        tind[active] += dt / t_per[active]
        tind[~active] = 0.0
        release = active & (tind >= 1.0)
        if release.any():
            r = k_r * lam[release]
            n_tor = R * dphi / lam[release]
            volume = n_tor * 2.0 * np.pi ** 2 * r ** 2 * R * np.sin(phi[release])
            dn = volume / (4.0 / 3.0 * np.pi * r ** 3)
            r_rel.append(r); n_rel.append(dn); t_rel.append(np.full(r.size, t))
            m -= float((v.rho_m * volume).sum())
            tind[release] -= 1.0
            m_hist.append((t, m))
        t += dt
    r = np.concatenate(r_rel) if r_rel else np.zeros(0)
    n = np.concatenate(n_rel) if n_rel else np.zeros(0)
    if r.size:
        order = np.argsort(r)
        cn = np.cumsum(n[order])
        r_med = float(r[order][np.searchsorted(cn, 0.5 * cn[-1])])
    else:
        r_med = None
    out.update({"N": float(n.sum()), "r_med_um": r_med * 1e6 if r_med else None, "r_min_um": float(r.min() * 1e6) if r.size else None,
                "r_max_um": float(r.max() * 1e6) if r.size else None, "t_sd_ms": t * 1e3, "tau_end": t / t_ch, "m0_kg": m0, "m_end_kg": m,
                "radii": r, "counts": n, "times": np.concatenate(t_rel) if t_rel else np.zeros(0), "mass_history": np.array(m_hist)})
    return out
