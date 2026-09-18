"""Regime numbers and the sphere drag coefficient the way SESAM computes it.

Kn = lambda / D with a hard-sphere mean free path; Ma with gamma = 1.4; C_D = C_D,c(Ma) +
(C_D,fm(Ma) - C_D,c(Ma)) f(Kn) with the ATDB_SPHERE tables and a Knudsen bridging f.
See docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md sections 6.4-6.5.
"""
import json
import math
import os

import numpy as np

from . import DATA_DIR
from .constants import GAMMA_AIR, HARD_SPHERE_DIAMETER, K_BOLTZMANN

DEFAULT_ATDB = os.path.join(DATA_DIR, "atdb_sphere.json")


def mean_free_path(rho, m_bar):
    """Hard-sphere mean free path [m] for mass density rho [kg/m3] and mean molecular mass m_bar [kg]."""
    n = rho / m_bar
    return 1.0 / (math.sqrt(2.0) * math.pi * HARD_SPHERE_DIAMETER ** 2 * n)


def knudsen(rho, m_bar, diameter):
    return mean_free_path(rho, m_bar) / diameter


def speed_of_sound(T, m_bar):
    return math.sqrt(GAMMA_AIR * K_BOLTZMANN * T / m_bar)


def mach(V, T, m_bar):
    return V / speed_of_sound(T, m_bar)


class SphereDragTables:
    """C_D,fm(Ma) and C_D,c(Ma) of DRAMA's sphere aerothermal database (Ma 5..30, linear, clamped).

    Below Ma 5 SESAM holds the continuum value at the Ma-5 entry down to Ma 1 and uses exactly
    half of it below Ma 1 (measured on 466 rows, facts note s.2); the free-molecular value is
    simply clamped (Kn is negligible wherever Ma < 5)."""

    def __init__(self, mach_points, cd_free_molecular, cd_continuum):
        self.mach = np.asarray(mach_points, dtype=float)
        self.cd_fm = np.asarray(cd_free_molecular, dtype=float)
        self.cd_c = np.asarray(cd_continuum, dtype=float)

    @classmethod
    def from_json(cls, path=None):
        with open(path or DEFAULT_ATDB) as fh:
            d = json.load(fh)
        return cls(d["mach"], d["cd_free_molecular"], d["cd_continuum"])

    def cd_free_molecular(self, ma):
        return float(np.interp(ma, self.mach, self.cd_fm))          # np.interp clamps at both ends

    def cd_continuum(self, ma):
        if ma < 1.0:
            return 0.5 * float(self.cd_c[0])
        if ma < self.mach[0]:
            return float(self.cd_c[0])
        return float(np.interp(ma, self.mach, self.cd_c))


class SesamTable:
    """f(Kn) measured from SESAM's own drag output: bin means of (C_D - C_D,c)/(C_D,fm - C_D,c) over
    1544 rows (Ma >= 5, Kn > 1e-3) of the four reference runs and the 7.5 km/s, 300 K sweep runs of
    2026-09-17, in 0.25-decade bins of log10 Kn; linear in log10 Kn between bin centres, 0 below,
    1 above. Reproduces SESAM's C_D to rms < 0.01."""
    LOG10_KN = np.array([-2.125, -1.875, -1.625, -1.375, -1.125, -0.875, -0.625, -0.375, -0.125, 0.125])
    F = np.array([0.0, 0.0021, 0.0287, 0.1021, 0.2517, 0.4567, 0.6850, 0.8897, 0.9760, 1.0])

    def __call__(self, kn):
        if kn <= 0.0:
            return 0.0
        return float(np.interp(math.log10(kn), self.LOG10_KN, self.F, left=0.0, right=1.0))


class SesamErf:
    """Analytic fit of the same data: f = 1/2 [1 + erf((log10 Kn - center)/width)] (rms 0.0065, max error 0.027 near Kn 0.75)."""

    def __init__(self, center=-0.845, width=0.585):
        self.center, self.width = center, width

    def __call__(self, kn):
        if kn <= 0.0:
            return 0.0
        return 0.5 * (1.0 + math.erf((math.log10(kn) - self.center) / self.width))


class Sin2:
    """sin^2 ramp over [kn_lo, kn_hi] in log10 Kn; the textbook form is sin^2[pi (0.5 + 0.25 log10 Kn)] = Sin2(0.01, 1)."""

    def __init__(self, kn_lo=0.01, kn_hi=1.0):
        self.lo, self.hi = math.log10(kn_lo), math.log10(kn_hi)

    def __call__(self, kn):
        if kn <= 0.0:
            return 0.0
        x = (math.log10(kn) - self.lo) / (self.hi - self.lo)
        x = min(1.0, max(0.0, x))
        return math.sin(0.5 * math.pi * x) ** 2


class Matting:
    """Matting (1971) bridging relation: interface reserved; transcribe from J. Spacecraft Rockets 8(1) 35-40."""

    def __call__(self, kn):
        raise NotImplementedError("Matting (1971) bridging is not transcribed yet; use sesam-table, sesam-erf or sin2")


BRIDGING_NAMES = ("sesam-table", "sesam-erf", "sin2", "textbook", "matting")


def bridging_by_name(name):
    if name == "sesam-table":
        return SesamTable()
    if name == "sesam-erf":
        return SesamErf()
    if name in ("sin2", "textbook"):
        return Sin2()
    if name == "matting":
        return Matting()
    raise ValueError("bridging must be one of {}, got {!r}".format(BRIDGING_NAMES, name))


def drag_coefficient(kn, ma, tables, bridging):
    """SESAM's sphere C_D: continuum and free-molecular table values blended by f(Kn)."""
    cd_c = tables.cd_continuum(ma)
    cd_fm = tables.cd_free_molecular(ma)
    f = bridging(kn) if kn > 0.0 else 0.0
    return cd_c + (cd_fm - cd_c) * f
