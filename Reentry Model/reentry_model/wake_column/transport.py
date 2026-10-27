"""How the plume spreads: eddy and molecular diffusion, the plume volume, and the settling loss rate.

The plume is a horizontal line segment of length L (the trail) broadened into a Gaussian of mean-square radius r^2 =
r0^2 + 4 integral(D dt). Everything here is a pure function of the plume radius r (m) and the ambient state, so `box.py`
can combine them (design decisions 2, 3 and 4 of the plan; spec 4.2).

  EddyProfile   Megner et al. (2006) reference eddy diffusion coefficient K(z), from data/wake_column/megner2006_eddy_diffusion.csv
  eddy_switch   S(r): 0 below l_inner/sqrt(2), 1 above sqrt(2) l_inner, cubic smoothstep in ln r between (decision 3)
  D_molecular   hard-sphere Chapman-Enskog diffusion of a species in air (decision 4)
  D_total       D_mol + eddy_scale K S(r)
  volume        V = pi r^2 L / erf(L / 2r): pi r^2 L for r << L, pi^1.5 r^3 for r >> L (decision 2)
  dlnV_dr2      d ln V / d(r^2), analytic, finite for L/r from 1e-6 to 1e6
  settling_rate 2 v_s / (pi r): the rate at which the settling speed v_s empties a plume of radius r
"""
import csv
import math
import os

import numpy as np

from reentry_model.wake_column import DATA
from reentry_model.wake_column.ambient import _interp_log

K_B = 1.380649e-23              # J/K
AMU = 1.66053906660e-27         # kg
AIR_SIGMA_A = 3.711             # angstrom, air as the pseudo-species (decision 4)
AIR_MASS_AMU = 28.96

EDDY_CSV = os.path.join(DATA, "megner2006_eddy_diffusion.csv")


class EddyProfile:
    """Eddy diffusion coefficient K(z) in m^2/s, log-linear in z between the digitised points."""

    def __init__(self, z_km, K_m2s, path="<memory>"):
        self.z_km = np.asarray(z_km, dtype=float)
        self.K_m2s = np.asarray(K_m2s, dtype=float)
        self.path = path
        if np.any(np.diff(self.z_km) <= 0):
            raise ValueError(f"altitudes in {path} are not increasing")

    @classmethod
    def load(cls, path=None):
        """The Megner et al. (2006) Fig. 3c reference profile (the committed digitised table)."""
        path = path or EDDY_CSV
        with open(path, newline="") as f:
            rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
        return cls([float(r["z_km"]) for r in rows], [float(r["K_m2s"]) for r in rows], path)

    def K(self, z_km):
        """K in m^2/s at altitude z_km; ValueError outside the table (no extrapolation)."""
        z = float(z_km)
        zmin, zmax = float(self.z_km[0]), float(self.z_km[-1])
        if not (zmin <= z <= zmax):         # also false for nan
            raise ValueError(f"altitude {z} km outside the Megner eddy-diffusion table {self.path} ({zmin}-{zmax} km)")
        return _interp_log(z, self.z_km, self.K_m2s)


def eddy_switch(r_m, l_inner_m):
    """Fraction of the eddy diffusion that acts on a plume of radius r_m: 0 below l_inner/sqrt(2), 1 above
    sqrt(2) l_inner, and a cubic smoothstep in ln r between (so 1/2 at r = l_inner)."""
    if r_m <= 0.0:
        return 0.0
    t = (math.log(r_m / l_inner_m) + 0.5 * math.log(2.0)) / math.log(2.0)
    if t <= 0.0:
        return 0.0
    if t >= 1.0:
        return 1.0
    return t * t * (3.0 - 2.0 * t)


def D_molecular(sp, amb):
    """Hard-sphere Chapman-Enskog diffusion coefficient (m^2/s) of species sp in air at ambient state amb:
    D = 3 / (8 n sigma12^2) sqrt(k_B T / (2 pi mu)), sigma12 the mean of the two collision diameters, mu the reduced
    mass, n the air number density; the collision integral is taken as 1 (decision 4)."""
    sigma12 = 0.5 * (sp.sigma_A + AIR_SIGMA_A) * 1e-10
    mu = sp.mass_amu * AIR_MASS_AMU / (sp.mass_amu + AIR_MASS_AMU) * AMU
    return 3.0 / (8.0 * amb.n_total_m3 * sigma12 ** 2) * math.sqrt(K_B * amb.T_K / (2.0 * math.pi * mu))


def D_total(r_m, sp, amb, eddy_K, settings):
    """Total diffusion coefficient (m^2/s) of species sp for a plume of radius r_m:
    D_mol + eddy_scale * K * S(r), with K the eddy coefficient (m^2/s) at the plume altitude."""
    return D_molecular(sp, amb) + settings.eddy_scale * eddy_K * eddy_switch(r_m, settings.l_inner_m)


def volume(r_m, L_m):
    """Plume volume (m^3) at radius r_m for a line segment of length L_m: pi r^2 L / erf(L / 2r)."""
    return math.pi * r_m * r_m * L_m / math.erf(L_m / (2.0 * r_m))


_X_SERIES = 1e-4


def dlnV_dr2(r_m, L_m):
    """d ln V / d(r^2) in 1/m^2. With x = L / 2r the analytic result is
    (1/r^2) [1 + x exp(-x^2) / (sqrt(pi) erf(x))], which tends to 1/r^2 (cylinder) for x large and 1.5/r^2 (sphere)
    for x small. Below x = 1e-4 the series 1.5 - x^2/3 + 4 x^4/45 is used, so nothing is a difference of nearly equal
    numbers."""
    x = L_m / (2.0 * r_m)
    if x < _X_SERIES:
        g = 1.5 - x * x / 3.0 + 4.0 * x ** 4 / 45.0
    else:
        g = 1.0 + x * math.exp(-x * x) / (math.sqrt(math.pi) * math.erf(x))
    return g / (r_m * r_m)


def settling_rate(v_s, r_m):
    """Loss rate (1/s) of the plume content by settling at speed v_s (m/s) out of a plume of radius r_m:
    2 v_s / (pi r) (the fraction of a disc of radius r that lies within v_s dt of its lower edge)."""
    return 2.0 * v_s / (math.pi * r_m)
