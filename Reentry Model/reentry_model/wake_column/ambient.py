"""Ambient air of the wake column: temperature, pressure and the number density of the nine ambient species against
altitude, read from the committed tables built by analysis/extract_waccm_columns.py (WACCM, the baseline) and
analysis/build_msis_afgl_columns.py (NRLMSIS 2.1 + AFGL midlatitude winter, the cross-check).

Each table holds a day block (local noon) and a night block (local midnight) on a 1 km grid from 45 to 100 km at
40.7 N. The model holds the chosen block fixed for the 12-hour run. `AmbientColumn.at(z_km)` interpolates linearly in
T and log-linearly in pressure and every density, and raises ValueError outside the table rather than extrapolating.
"""
import csv
import os
from dataclasses import dataclass

import numpy as np

from reentry_model.wake_column import DATA
from reentry_model.wake_column.species import AMBIENT_SPECIES

SOURCES = ("waccm", "msis_afgl")
TIMES_OF_DAY = ("day", "night")


def table_path(source):
    """Path of the committed ambient table for source ("waccm" or "msis_afgl")."""
    if source not in SOURCES:
        raise ValueError(f"unknown ambient source {source!r}; expected one of {SOURCES}")
    return os.path.join(DATA, f"ambient_{source}_40.7N.csv")


@dataclass(frozen=True)
class AmbientState:
    z_km: float
    T_K: float
    p_Pa: float
    n_total_m3: float
    n: dict                      # species name -> number density in m^-3, every name in AMBIENT_SPECIES


def _interp_log(z, zs, ys):
    """Log-linear interpolation in z (linear wherever either neighbouring value is not positive)."""
    k = int(np.clip(np.searchsorted(zs, z) - 1, 0, len(zs) - 2))
    f = (z - zs[k]) / (zs[k + 1] - zs[k])
    y0, y1 = ys[k], ys[k + 1]
    if y0 > 0 and y1 > 0:
        return float(np.exp((1 - f) * np.log(y0) + f * np.log(y1)))
    return float((1 - f) * y0 + f * y1)


class AmbientColumn:
    """One ambient column (a source and a time of day) on its altitude grid."""

    def __init__(self, z_km, p_Pa, T_K, n_total_m3, n, path="<memory>"):
        self.z_km, self.p_Pa, self.T_K, self.n_total_m3 = (np.asarray(a, dtype=float)
                                                           for a in (z_km, p_Pa, T_K, n_total_m3))
        self.n = {s: np.asarray(n[s], dtype=float) for s in AMBIENT_SPECIES}
        self.path = path
        if np.any(np.diff(self.z_km) <= 0):
            raise ValueError(f"altitudes in {path} are not increasing")

    @classmethod
    def load(cls, source, time_of_day, path=None):
        """Read the day or night block of the ambient table for source ("waccm" or "msis_afgl")."""
        if time_of_day not in TIMES_OF_DAY:
            raise ValueError(f"unknown time of day {time_of_day!r}; expected one of {TIMES_OF_DAY}")
        path = path or table_path(source)
        with open(path, newline="") as f:
            rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))
                    if r["time_of_day"] == time_of_day]
        if not rows:
            raise ValueError(f"no {time_of_day} rows in {path}")
        col = lambda name: [float(r[name]) for r in rows]
        return cls(col("z_km"), col("p_Pa"), col("T_K"), col("n_total_m3"), {s: col(s) for s in AMBIENT_SPECIES}, path)

    def at(self, z_km):
        """The ambient state at altitude z_km; ValueError outside the table."""
        z = float(z_km)
        zmin, zmax = float(self.z_km[0]), float(self.z_km[-1])
        if not (zmin <= z <= zmax):         # also false for nan
            raise ValueError(f"altitude {z} km outside {self.path} ({zmin}-{zmax} km)")
        zs = self.z_km
        k = int(np.clip(np.searchsorted(zs, z) - 1, 0, len(zs) - 2))
        f = (z - zs[k]) / (zs[k + 1] - zs[k])
        T = float((1 - f) * self.T_K[k] + f * self.T_K[k + 1])
        return AmbientState(z_km=z, T_K=T, p_Pa=_interp_log(z, zs, self.p_Pa),
                            n_total_m3=_interp_log(z, zs, self.n_total_m3),
                            n={s: _interp_log(z, zs, self.n[s]) for s in AMBIENT_SPECIES})
