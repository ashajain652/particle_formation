"""Girin's (2017, A&A 606, A63, Eq. 1) dispersion relation of the gradient instability of a sheared liquid layer,
solved numerically for the fastest-growing disturbance (spec Step 3 section 9).

    (Omega - Delta) [(Omega - Delta)(Omega + K Delta) + (1 - K) Delta] = Delta^3 We_s^-1 (Omega - K Delta),
    K = [1 - exp(-2 Delta)] / (2 Delta),   Omega = omega delta_m / V_s,   Delta = 2 pi delta_m / lambda,
    We_s = rho_m V_s^2 delta_m / Sigma.

For every We_s the cubic in Omega is solved on a grid of Delta (companion-matrix eigenvalues, vectorised) and the
Delta with the largest Im Omega is the fastest mode: Delta_f (wavelength 2 pi delta_m / Delta_f), Im Omega_f (growth
rate Im Omega_f V_s / delta_m) and Re Omega_f. Measured 2026-09-20: instability appears between We_s 3.00 and 3.08
(Girin: 3.08); Delta_f -> 1.226, Im Omega_f -> 0.247 at We_s 1e4 (Girin's Fig. 5 asymptotes 1.225, 0.24);
Re/Im = 1.68 at We_s 20 and 1.56 at 1e4 (Girin: 1.5); at We_s 4.62: 0.58 / 0.049; 10: 1.34 / 0.18; 20: 1.29 / 0.22.
The table is cached in the package data (`girin_dispersion.json`) and interpolated in log We_s."""
import json
import os

import numpy as np

from . import DATA_DIR

TABLE_PATH = os.path.join(DATA_DIR, "girin_dispersion.json")
WE_CRITICAL_THEORY = 3.08          # Girin's critical surface Weber number
WE_CRITICAL_PRACTICAL = 4.62       # the value he recommends in practice (GI >= 0.6)
DELTA_GRID = np.linspace(0.02, 5.0, 2500)
WE_GRID = np.logspace(np.log10(3.0), 4.0, 200)
UNSTABLE_TOL = 1e-6                # Im Omega above which a root counts as unstable


def cubic_coefficients(delta, we):
    """Coefficients (a2, a1, a0) of Omega^3 + a2 Omega^2 + a1 Omega + a0 = 0, Girin's Eq. (1) expanded."""
    D = np.asarray(delta, dtype=float)
    K = (1.0 - np.exp(-2.0 * D)) / (2.0 * D)
    a2 = (K - 2.0) * D
    a1 = (1.0 - 2.0 * K) * D ** 2 + (1.0 - K) * D - D ** 3 / we
    a0 = K * D ** 3 - (1.0 - K) * D ** 2 + K * D ** 4 / we
    return a2, a1, a0


def growth_rates(we, delta=DELTA_GRID):
    """(Im Omega, Re Omega) of the most unstable root of Eq. (1) at each Delta for one We_s."""
    a2, a1, a0 = cubic_coefficients(delta, we)
    C = np.zeros((len(delta), 3, 3), dtype=complex)
    C[:, 0, 1] = 1.0
    C[:, 1, 2] = 1.0
    C[:, 2, 0], C[:, 2, 1], C[:, 2, 2] = -a0, -a1, -a2
    roots = np.linalg.eigvals(C)
    k = np.argmax(roots.imag, axis=1)
    best = roots[np.arange(len(delta)), k]
    return best.imag, best.real


def fastest_mode(we, delta=DELTA_GRID):
    """(Delta_f, Im Omega_f, Re Omega_f) of the fastest-growing disturbance; (nan, 0, 0) when nothing is unstable."""
    im, re = growth_rates(we, delta)
    j = int(np.argmax(im))
    if im[j] <= UNSTABLE_TOL:
        return float("nan"), 0.0, 0.0
    return float(delta[j]), float(im[j]), float(re[j])


def build_table(we_grid=WE_GRID, delta=DELTA_GRID):
    rows = np.array([fastest_mode(we, delta) for we in we_grid])
    return {"we": list(map(float, we_grid)), "delta_f": list(map(float, rows[:, 0])),
            "im_omega_f": list(map(float, rows[:, 1])), "re_omega_f": list(map(float, rows[:, 2])),
            "delta_grid": [float(delta[0]), float(delta[-1]), int(len(delta))]}


def write_table(path=TABLE_PATH):
    table = build_table()
    with open(path, "w") as fh:
        json.dump(table, fh)
    return table


class DispersionTable:
    """Delta_f(We_s), Im Omega_f(We_s), Re Omega_f(We_s) interpolated in log We_s from the cached table."""

    def __init__(self, path=TABLE_PATH):
        if not os.path.isfile(path):
            write_table(path)
        with open(path) as fh:
            t = json.load(fh)
        self.we = np.array(t["we"])
        self.delta_f, self.im_omega_f, self.re_omega_f = (np.array(t[k]) for k in ("delta_f", "im_omega_f", "re_omega_f"))
        unstable = self.im_omega_f > 0.0
        self.we_onset = float(self.we[unstable][0])              # first tabulated We_s with an unstable mode
        self._log_we = np.log(self.we[unstable])
        self._delta, self._im, self._re = self.delta_f[unstable], self.im_omega_f[unstable], self.re_omega_f[unstable]

    def __call__(self, we):
        """Arrays (Delta_f, Im Omega_f, Re Omega_f) for We_s (scalar or array); nan/0 below the onset."""
        we = np.asarray(we, dtype=float)
        x = np.log(np.maximum(we, self.we_onset))
        d, im, re = (np.interp(x, self._log_we, v) for v in (self._delta, self._im, self._re))
        stable = we < self.we_onset
        return np.where(stable, np.nan, d), np.where(stable, 0.0, im), np.where(stable, 0.0, re)
