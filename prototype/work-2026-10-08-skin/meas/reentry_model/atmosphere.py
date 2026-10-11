"""Freestream state along the trajectory: NRLMSISE-00 (pymsis), the US76 table SESAM ships, or a
replay of a SESAM run's own density/temperature (diagnostic). Winds: none or the static profile."""
import csv
import math
import os
from dataclasses import dataclass
from datetime import timedelta

import numpy as np

from . import DATA_DIR
from .constants import ATOMIC_MASS_UNIT, GAMMA_AIR, HARD_SPHERE_DIAMETER, K_BOLTZMANN, M_BAR_AIR

DEFAULT_US76 = os.path.join(DATA_DIR, "us76_static_environment.csv")


@dataclass(frozen=True)
class Freestream:
    rho: float           # kg/m3
    T: float             # K
    p: float             # Pa
    m_bar: float         # kg, mean molecular mass
    wind_enu: np.ndarray  # m/s, (east, north, up)


def _read_us76(path):
    rows = []
    with open(path) as fh:
        for r in csv.reader(fh):
            if r and not r[0].startswith("#"):
                rows.append([float(v) for v in r])
    a = np.array(rows)
    return {"h": a[:, 0], "p": a[:, 1], "rho": a[:, 2], "T": a[:, 3], "wn": a[:, 4], "we": a[:, 5], "wd": a[:, 6]}


class NoWind:
    def wind_enu(self, h):
        return np.zeros(3)


class StaticProfileWind:
    """The wind columns of DRAMA's static environment table (north, east, down -> east, north, up)."""

    def __init__(self, path=None):
        t = _read_us76(path or DEFAULT_US76)
        self._h, self._wn, self._we, self._wd = t["h"], t["wn"], t["we"], t["wd"]

    def wind_enu(self, h):
        return np.array([np.interp(h, self._h, self._we), np.interp(h, self._h, self._wn), -np.interp(h, self._h, self._wd)])


class US76TableAtmosphere:
    """US Standard Atmosphere 1976 as tabulated by DRAMA (0-150 km, 100 m); log-linear rho and p, linear T."""

    def __init__(self, path=None, wind=None):
        t = _read_us76(path or DEFAULT_US76)
        self._h, self._T = t["h"], t["T"]
        self._log_rho, self._log_p = np.log(t["rho"]), np.log(t["p"])
        self.wind = wind or NoWind()

    def state(self, t, h, lat, lon):
        if h < self._h[0] or h > self._h[-1]:
            raise ValueError("altitude {:.1f} m outside the US76 table (0..{:.0f} m)".format(h, self._h[-1]))
        rho = math.exp(np.interp(h, self._h, self._log_rho))
        p = math.exp(np.interp(h, self._h, self._log_p))
        T = float(np.interp(h, self._h, self._T))
        return Freestream(rho, T, p, M_BAR_AIR, self.wind.wind_enu(h))


class NRLMSISE00Atmosphere:
    """NRLMSISE-00 through pymsis (version 0) at the run epoch + t, with the fap-file solar indices."""

    def __init__(self, epoch, solar, wind=None):
        from pymsis import msis        # imported here so the rest of the package works without pymsis
        self._msis = msis
        self.epoch = epoch
        self.solar = solar
        self.wind = wind or NoWind()

    def state(self, t, h, lat, lon):
        when = np.datetime64(self.epoch + timedelta(seconds=float(t)))
        out = self._msis.calculate(np.array([when]), [math.degrees(lon)], [math.degrees(lat)], [h / 1e3],
                                   f107s=[self.solar.f107], f107as=[self.solar.f107a],
                                   aps=[[self.solar.ap] * 7], version=0)
        row = np.asarray(out).reshape(-1, 11)[0]
        rho = float(row[0])
        n = float(np.nansum(row[1:10]))                 # species number densities (NO is NaN for version 0)
        T = float(row[10])
        m_bar = rho / n
        return Freestream(rho, T, n * K_BOLTZMANN * T, m_bar, self.wind.wind_enu(h))


class ReplayAtmosphere:
    """rho(h), m_bar(h), T(h) reconstructed from a SESAM run: rho from its density column; m_bar = air below
    90 km and, above, back-solved from the Knudsen column (lambda = Kn D, n = 1/(sqrt2 pi d^2 lambda));
    T from Ma and V (T = (V/Ma)^2 m_bar / (gamma k_B)) on rows with Ma > 0.3. Interpolated in altitude."""

    def __init__(self, reference, wind=None):
        order = np.argsort(reference.altitude)
        h = reference.altitude[order]
        rho = reference.density[order]
        kn = reference.knudsen[order]
        ma = reference.mach[order]
        V = reference.velocity[order]
        keep = np.concatenate([[True], np.diff(h) > 0])          # strictly increasing altitude nodes
        h, rho, kn, ma, V = h[keep], rho[keep], kn[keep], ma[keep], V[keep]
        m_bar = np.full_like(h, M_BAR_AIR)
        high = (h >= 90e3) & (kn > 0)
        lam = kn[high] * reference.diameter
        n = 1.0 / (math.sqrt(2.0) * math.pi * HARD_SPHERE_DIAMETER ** 2 * lam)
        m_bar[high] = np.clip(rho[high] / n, 16.0 * ATOMIC_MASS_UNIT, 30.0 * ATOMIC_MASS_UNIT)
        ok = ma > 0.3
        T = (V[ok] / ma[ok]) ** 2 * m_bar[ok] / (GAMMA_AIR * K_BOLTZMANN)
        self._h, self._log_rho, self._m_bar = h, np.log(rho), m_bar
        self._hT, self._T = h[ok], T
        self.wind = wind or NoWind()

    def state(self, t, h, lat, lon):
        rho = math.exp(np.interp(h, self._h, self._log_rho))
        m_bar = float(np.interp(h, self._h, self._m_bar))
        T = float(np.interp(h, self._hT, self._T))
        return Freestream(rho, T, rho / m_bar * K_BOLTZMANN * T, m_bar, self.wind.wind_enu(h))


class VacuumAtmosphere:
    """No atmosphere at all (dynamics tests)."""

    def state(self, t, h, lat, lon):
        return Freestream(0.0, 200.0, 0.0, M_BAR_AIR, np.zeros(3))
