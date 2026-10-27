"""Droplet adapter: turns one droplet run's history rows into the trail that the later stages need.

The droplet model (reentry_model/droplet_reentry) returns rows at integrator steps (irregular spacing, `m_Al_kg`
falling as aluminium evaporates, `mdot_Al_kgs` the positive loss rate). Everything here is built from the cumulative
aluminium loss m_Al(t0) - m_Al(t):

  * the inventory `N_Al`, `m_Al_kg` is the whole loss to the last row (design decision 16);
  * the trail length `L_m` is the path flown, the trapezoid integral of V dt, from release to the moment the
    cumulative loss first reaches `loss_fraction` (0.9) of the total, interpolated within the step;
  * `z_mean_km` is the aluminium-mass-weighted mean altitude over the whole loss (midpoint altitude of each step);
  * the emission samples (decision 15) sit at the midpoints of `n_samples` equal-aluminium-mass intervals of the
    cumulative loss, with time, altitude, speed and droplet temperature interpolated against cumulative loss and
    weight 1/n_samples.

A run that was skipped (rows is None) or that loses no aluminium returns an empty trail: zero inventory, no samples,
and z_mean_km equal to the release altitude.
"""
from dataclasses import dataclass

import numpy as np
from scipy.constants import N_A

M_AL_MOLAR = 26.9815385e-3                 # kg/mol
M_AL_ATOM = M_AL_MOLAR / N_A               # kg per aluminium atom


@dataclass(frozen=True)
class EmissionSample:
    """One equal-mass slice of the aluminium emission: when and where it left, the droplet's speed and temperature
    then, and its share (weight) of the total."""
    t_s: float
    h_km: float
    V_ms: float
    T_d_K: float
    weight: float


@dataclass(frozen=True)
class Trail:
    """The trail one droplet leaves: length [m], aluminium atoms and mass vaporised [kg], zinc and magnesium
    vaporised [kg], aluminium-weighted mean altitude [km] and the emission samples."""
    d0_um: float
    h0_km: float
    L_m: float
    N_Al: float
    m_Al_kg: float
    m_Zn_kg: float
    m_Mg_kg: float
    z_mean_km: float
    samples: tuple


def _loss(rows, key):
    """Non-negative mass lost [kg] by the last row, for rows that carry `key` (0 when they do not)."""
    if key not in rows[0]:
        return 0.0
    return max(float(rows[0][key]) - float(rows[-1][key]), 0.0)


def _at_loss(cum, values, x):
    """`values` interpolated linearly in the step where the cumulative loss `cum` (non-decreasing, cum[0] = 0)
    reaches each x in (0, cum[-1]]. Flat stretches of `cum` are skipped, so the interpolation is never 0/0."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    k = np.clip(np.searchsorted(cum, x, side="left"), 1, len(cum) - 1)
    c0, c1 = cum[k - 1], cum[k]
    f = (x - c0) / (c1 - c0)
    return values[k - 1] + f * (values[k] - values[k - 1])


def trail_from_rows(rows, d0_um, h0_km, n_samples=20, loss_fraction=0.9):
    """The Trail of a droplet run from its history `rows` (None when the run was skipped)."""
    if rows is None or len(rows) < 2:
        return Trail(d0_um, h0_km, 0.0, 0.0, 0.0, 0.0, 0.0, float(h0_km), ())
    m_zn, m_mg = _loss(rows, "m_Zn_kg"), _loss(rows, "m_Mg_kg")
    mAl = np.array([r["m_Al_kg"] for r in rows], dtype=float)
    total = mAl[0] - mAl[-1]
    if not total > 0.0:
        return Trail(d0_um, h0_km, 0.0, 0.0, 0.0, m_zn, m_mg, float(h0_km), ())

    t = np.array([r["t_s"] for r in rows], dtype=float)
    h = np.array([r["h_km"] for r in rows], dtype=float)
    V = np.array([r["V_ms"] for r in rows], dtype=float)
    Td = np.array([r["T_K"] for r in rows], dtype=float)

    cum = np.maximum.accumulate(mAl[0] - mAl)      # cumulative loss, forced non-decreasing against round-off
    path = np.concatenate(([0.0], np.cumsum(0.5 * (V[1:] + V[:-1]) * np.diff(t))))
    L = float(_at_loss(cum, path, loss_fraction * cum[-1])[0])

    dm = np.diff(cum)
    z_mean = float(np.sum(0.5 * (h[1:] + h[:-1]) * dm) / dm.sum())

    x = (np.arange(n_samples) + 0.5) / n_samples * cum[-1]
    ts, hs, Vs, Ts = (_at_loss(cum, a, x) for a in (t, h, V, Td))
    samples = tuple(EmissionSample(float(ts[i]), float(hs[i]), float(Vs[i]), float(Ts[i]), 1.0 / n_samples)
                    for i in range(n_samples))
    return Trail(d0_um, h0_km, L, float(total / M_AL_ATOM), float(total), m_zn, m_mg, z_mean, samples)


def trail_for(d0_um, h0_km, n_samples=20, parent=None):
    """Run the droplet model for a d0_um [micrometre] droplet released at h0_km and return its Trail."""
    from reentry_model.droplet_reentry.model import run_case
    _, rows = run_case(d0_um, h0_km, parent=parent)
    return trail_from_rows(rows, d0_um, h0_km, n_samples=n_samples)
