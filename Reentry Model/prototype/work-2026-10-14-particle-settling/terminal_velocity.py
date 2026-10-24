"""Terminal (gravitational settling) velocity of solid spheres, 0.01 um to 100 um diameter.

Force balance between weight less buoyancy and aerodynamic drag, with the drag law covering the
continuum, slip and free-molecule regimes:

    (pi/6) d^3 (rho_p - rho) g  =  C_D (pi/8) rho v^2 d^2
    C_D = (24/Re) (1 + 0.15 Re^0.687) / C_c          Schiller-Naumann, slip-corrected

which rearranges to the fixed point

    v = (rho_p - rho) d^2 g C_c / (18 mu (1 + 0.15 Re^0.687)),    Re = rho v d / mu

C_c is the Cunningham-Millikan-Davies slip correction, C_c = 1 + Kn (1.257 + 0.400 exp(-1.100/Kn))
with Kn = 2 lambda / d, so the small-particle limit is the Epstein free-molecule drag and the
large-particle limit is Stokes. Air comes from the repository's US76 table; viscosity from
Sutherland's law; the mean free path by Willeke's scaling of the reference value the Davies
coefficients were fitted with (0.0665 um at 293.15 K, 101.325 kPa), which is what the slip
correction as tabulated in Hinds, *Aerosol Technology*, 2nd ed., assumes.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))
REPO = "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution"
sys.path.insert(0, REPO)

from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.constants import G0, K_BOLTZMANN, M_BAR_AIR, WGS84_A

RHO_PARTICLE = 2813.0          # kg/m3, AA7075 (reentry_model/data/materials/AA7075_nomelt.json)

# Sutherland's law for air, valid to a couple of percent over ~150-600 K
MU_REF, T_REF, S_SUTH = 1.716e-5, 273.15, 110.4


def viscosity(T):
    return MU_REF * (T / T_REF) ** 1.5 * (T_REF + S_SUTH) / (T + S_SUTH)


def mean_thermal_speed(T):
    return math.sqrt(8.0 * K_BOLTZMANN * T / (math.pi * M_BAR_AIR))


LAMBDA_REF, T_LAMBDA_REF, P_LAMBDA_REF = 6.65e-8, 293.15, 101325.0


def mean_free_path(T, p):
    """Willeke (1976) scaling of the reference mean free path; the convention Davies' coefficients use."""
    return (LAMBDA_REF * (P_LAMBDA_REF / p) * (T / T_LAMBDA_REF)
            * (1.0 + S_SUTH / T_LAMBDA_REF) / (1.0 + S_SUTH / T))


def kinetic_mean_free_path(rho, T, mu):
    """lambda = 2 mu / (rho c_bar); kept for comparison only, it is ~2 % smaller."""
    return 2.0 * mu / (rho * mean_thermal_speed(T))


def cunningham(Kn):
    return 1.0 + Kn * (1.257 + 0.400 * math.exp(-1.100 / Kn))


def gravity(h):
    return G0 * (WGS84_A / (WGS84_A + h)) ** 2


def terminal_velocity(d, rho, T, p, mu, g, rho_p=RHO_PARTICLE):
    """Settle the fixed point above. Returns (v, Re, Kn, Cc, residual)."""
    lam = mean_free_path(T, p)
    Kn = 2.0 * lam / d
    Cc = cunningham(Kn)
    stokes = (rho_p - rho) * d * d * g * Cc / (18.0 * mu)
    v = stokes
    for it in range(200):
        Re = rho * v * d / mu
        v_new = stokes / (1.0 + 0.15 * Re ** 0.687)
        if abs(v_new - v) <= 1e-12 * max(v_new, 1e-30):
            v = v_new
            break
        v = 0.5 * (v + v_new)          # damped, so the high-Re cases cannot oscillate
    # residual of the force balance itself, normalised by weight: independent check on the iteration
    Re = rho * v * d / mu
    weight = (math.pi / 6.0) * d ** 3 * (rho_p - rho) * g
    C_D = (24.0 / Re) * (1.0 + 0.15 * Re ** 0.687) / Cc
    drag = C_D * (math.pi / 8.0) * rho * v * v * d * d
    return v, Re, Kn, Cc, abs(drag - weight) / weight


def regime(Kn):
    if Kn < 0.01:
        return "continuum"
    if Kn < 0.1:
        return "slip"
    if Kn < 10.0:
        return "transitional"
    return "free-molecule"


ALTITUDES_KM = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
DIAMETERS_UM = [10.0 ** (-2.0 + k / 5.0) for k in range(21)]      # 0.01 to 100 um, 5 per decade

atm = US76TableAtmosphere()
rows = []
residuals = []
for h_km in ALTITUDES_KM:
    h = h_km * 1e3
    fs = atm.state(0.0, h, 0.0, 0.0)
    mu = viscosity(fs.T)
    g = gravity(h)
    c_bar = mean_thermal_speed(fs.T)
    for d_um in DIAMETERS_UM:
        d = d_um * 1e-6
        v, Re, Kn, Cc, resid = terminal_velocity(d, fs.rho, fs.T, fs.p, mu, g)
        residuals.append(resid)
        tau = RHO_PARTICLE * d * d * Cc / (18.0 * mu)             # velocity relaxation time [s]
        rows.append(dict(
            h_km=h_km, d_um=d_um, rho_air=fs.rho, T_air=fs.T, mu=mu,
            lambda_m=mean_free_path(fs.T, fs.p), Kn=Kn, Cc=Cc, regime=regime(Kn),
            v_t=v, v_stokes=(RHO_PARTICLE - fs.rho) * d * d * g / (18.0 * mu),
            Re=Re, tau_s=tau, stop_dist_m=v * tau, speed_ratio=v / c_bar,
            fall_1km_s=1e3 / v,
        ))

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "terminal_velocity.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# sanity checks printed so the numbers can be audited against textbook values
sl = atm.state(0.0, 0.0, 0.0, 0.0)
mu_sl = viscosity(sl.T)
lam_sl = mean_free_path(sl.T, sl.p)
print("sea level: rho = {:.4f} kg/m3, T = {:.2f} K, p = {:.1f} Pa, mu = {:.4e} Pa s".format(
    sl.rho, sl.T, sl.p, mu_sl))
print("  lambda = {:.4e} m (Willeke), {:.4e} m (kinetic)".format(
    lam_sl, kinetic_mean_free_path(sl.rho, sl.T, mu_sl)))
for d_um, hinds in ((10.0, 1.016), (1.0, 1.168), (0.1, 2.928), (0.01, 22.98)):
    print("  Cc({:5.2f} um) = {:8.3f}   Hinds 2nd ed. table 3.1: {:8.3f}".format(
        d_um, cunningham(2 * lam_sl / (d_um * 1e-6)), hinds))
print("max force-balance residual over all {} cases: {:.2e} (relative to weight)".format(
    len(residuals), max(residuals)))
print("wrote {} ({} rows)".format(out, len(rows)))
