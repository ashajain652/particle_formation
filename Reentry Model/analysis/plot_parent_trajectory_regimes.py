"""Altitude against flight time of the generic spacecraft's DRAMA re-entry, with its flow-regime bands shaded.

The trajectory is "Main Body" from the committed DRAMA snapshot (data/generic_satellite_reentry/dmf_output.json),
time counted from SESAM's hand-over at 125 km. DRAMA kept the spacecraft intact (every break-up trigger switched
off), so no break-up is marked.

Regimes use the textbook limits on DRAMA's own Knudsen number of the body: free-molecular Kn > 10, transitional
0.01 <= Kn <= 10, continuum Kn < 0.01. A band edge inside the flown altitudes is where DRAMA's Knudsen column crosses
the limit (log Kn interpolated linearly in altitude on the descending branch). The body never reaches Kn 10 (peak
2.79 at 135.7 km), so that edge is extrapolated with Kn proportional to 1 / density on DRAMA's US76 table, anchored
at the highest flown point; the table stops at 150 km, so beyond it the density follows the scale height fitted over
the table's top 10 km (about 16 km).

The continuum band is split by the wake: laminar where the freestream Reynolds number Re = rho V D / mu is below
5 x 10^4, turbulent above it (Lees 1964's minimum critical Reynolds number for blunt-body wakes, on freestream
conditions and body diameter; citations/wake_radius_growth_laws.md). D is the equal-area diameter of DRAMA's
6.215 m^2 cross-section (2.81 m), rho and V are DRAMA's, mu is Sutherland's law at the US76 temperature. Re crosses
the threshold once on the descent and rises all the way to the ground, so the turbulent band runs down to 0 km.

    python analysis/plot_parent_trajectory_regimes.py   ->  output/plots/parent_trajectory_regimes.{png,pdf}
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from reentry_model import DATA_DIR, OUTPUT_DIR
from reentry_model.droplet_reentry.parent_trajectory import load_parent

KN_FREE_MOLECULAR = 10.0
KN_CONTINUUM = 0.01
RE_TURBULENT_WAKE = 5e4
CROSS_SECTION_M2 = 6.215                   # DRAMA "EVENT area" of the main body
US76_TABLE = os.path.join(DATA_DIR, "us76_static_environment.csv")
OUT_STEM = os.path.join(OUTPUT_DIR, "plots", "parent_trajectory_regimes")


def crossing_altitude(alt, kn, limit):
    """Altitude [km] where Kn falls through `limit` on a descending branch (alt decreasing, Kn decreasing)."""
    for i in range(len(alt) - 1):
        if kn[i] >= limit > kn[i + 1]:
            f = (np.log10(limit) - np.log10(kn[i])) / (np.log10(kn[i + 1]) - np.log10(kn[i]))
            return alt[i] + f * (alt[i + 1] - alt[i])
    return None


def reynolds_numbers(rows):
    """Freestream Reynolds number on the equal-area diameter, per trajectory row."""
    table = np.loadtxt(US76_TABLE, delimiter=",", comments="#")
    diameter = np.sqrt(4.0 * CROSS_SECTION_M2 / np.pi)
    h_m = np.array([r["altitude"] for r in rows]) * 1e3
    temperature = np.interp(h_m, table[:, 0], table[:, 3])
    mu = 1.458e-6 * temperature ** 1.5 / (temperature + 110.4)      # Sutherland
    rho = np.array([r["density"] for r in rows])
    v = np.array([r["velocity"] for r in rows]) * 1e3
    return rho * v * diameter / mu


def extrapolated_altitude(alt_anchor_km, kn_anchor, limit):
    """Altitude [km] where Kn = kn_anchor * rho(anchor) / rho(h) reaches `limit`, on DRAMA's US76 table."""
    table = np.loadtxt(US76_TABLE, delimiter=",", comments="#")
    h_km, rho = table[:, 0] / 1e3, table[:, 2]
    rho_anchor = np.exp(np.interp(alt_anchor_km, h_km, np.log(rho)))
    kn = kn_anchor * rho_anchor / rho
    if limit <= kn[-1]:
        return float(np.interp(np.log10(limit), np.log10(kn), h_km))
    # the table stops at 150 km: continue with the density scale height fitted over its top 10 km
    top = h_km >= h_km[-1] - 10.0
    scale_height = -1.0 / np.polyfit(h_km[top], np.log(rho[top]), 1)[0]
    return float(h_km[-1] + scale_height * np.log(limit / kn[-1]))


def main():
    parent = load_parent()
    t = np.array([r["time"] for r in parent.rows])
    h = np.array([r["altitude"] for r in parent.rows])
    kn = np.array([r["knudsen"] for r in parent.rows])

    top = int(np.argmax(h))
    descent_h, descent_kn = h[top:], kn[top:]
    h_continuum = crossing_altitude(descent_h, descent_kn, KN_CONTINUUM)
    h_free_molecular = crossing_altitude(descent_h, descent_kn, KN_FREE_MOLECULAR)
    extrapolated = h_free_molecular is None
    if extrapolated:
        h_free_molecular = extrapolated_altitude(h[top], kn[top], KN_FREE_MOLECULAR)

    # Re rises through the threshold on the descent; crossing_altitude wants a falling quantity, so pass 1/Re
    re = reynolds_numbers(parent.rows)[top:]
    rises = np.flatnonzero((re[:-1] < RE_TURBULENT_WAKE) & (re[1:] >= RE_TURBULENT_WAKE))
    falls = np.flatnonzero((re[:-1] >= RE_TURBULENT_WAKE) & (re[1:] < RE_TURBULENT_WAKE))
    if len(rises) != 1 or len(falls) != 0:
        raise ValueError("expected one upward crossing of Re = {:g} on the descent, found {} up and {} down"
                         .format(RE_TURBULENT_WAKE, len(rises), len(falls)))
    h_turbulent = crossing_altitude(descent_h, 1.0 / re, 1.0 / RE_TURBULENT_WAKE)

    y_max = h_free_molecular + 15.0
    fig, ax = plt.subplots(figsize=(10, 6))
    bands = [
        (h_free_molecular, y_max, "#dbe8f5", "Free-molecular  (Kn > 10)"),
        (h_continuum, h_free_molecular, "#f3e6cf", "Transitional  (0.01 ≤ Kn ≤ 10)"),
        (h_turbulent, h_continuum, "#e2eedb", "Continuum, laminar wake  (Kn < 0.01, Re < 5×10⁴)"),
        (0.0, h_turbulent, "#c9dfbd", "Continuum, turbulent wake  (Re > 5×10⁴)"),
    ]
    for i, (lo, hi, colour, label) in enumerate(bands):
        ax.axhspan(lo, hi, color=colour, zorder=0, lw=0)
        # the continuum bands are labelled on the left, clear of the steep descent on the right
        x, ha = (0.99, "right") if i < 2 else (0.08, "left")
        ax.text(x, (lo + hi) / 2, label, transform=ax.get_yaxis_transform(), ha=ha, va="center",
                fontsize=11, color="#333333")
    for edge in (h_free_molecular, h_continuum):
        ax.axhline(edge, color="#888888", lw=0.8, ls="--", zorder=1)
    ax.axhline(h_turbulent, color="#888888", lw=0.8, ls=":", zorder=1)
    ax.annotate("{:.1f} km".format(h_turbulent), (0.01, h_turbulent), xycoords=("axes fraction", "data"),
                xytext=(0, 3), textcoords="offset points", fontsize=9, color="#555555")
    ax.annotate("{:.1f} km{}".format(h_free_molecular, " (extrapolated)" if extrapolated else ""),
                (0.01, h_free_molecular), xycoords=("axes fraction", "data"), xytext=(0, 3),
                textcoords="offset points", fontsize=9, color="#555555")
    ax.annotate("{:.1f} km".format(h_continuum), (0.01, h_continuum), xycoords=("axes fraction", "data"),
                xytext=(0, 3), textcoords="offset points", fontsize=9, color="#555555")

    ax.plot(t, h, color="#1f3b73", lw=2, zorder=3, label="DRAMA trajectory, Main Body (intact)")
    ax.plot(t[-1], h[-1], "o", ms=8, color="#1f3b73", mec="white", mew=1.5, zorder=4)
    ax.annotate("ground impact, {:.0f} s".format(t[-1]), (t[-1], h[-1]), xytext=(-40, 14),
                textcoords="offset points", ha="right", fontsize=9, color="#333333")

    ax.set_xlim(0, t[-1] * 1.02)
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Flight time since SESAM hand-over at {:.0f} km (s)".format(h[0]))
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Generic spacecraft re-entry (DRAMA 4.1.4, US76 static atmosphere) and flow regimes of the body")
    ax.legend(loc="lower left", frameon=False)
    ax.grid(color="#cccccc", lw=0.5, alpha=0.6, zorder=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()

    os.makedirs(os.path.dirname(OUT_STEM), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_STEM + "." + ext, dpi=200)
    print("turbulent-wake edge (Re = {:g}) {:.2f} km".format(RE_TURBULENT_WAKE, h_turbulent))
    print("free-molecular edge {:.2f} km{}, continuum edge {:.2f} km; peak Kn {:.3f} at {:.2f} km; impact at {:.1f} s"
          .format(h_free_molecular, " (extrapolated)" if extrapolated else "", h_continuum, kn[top], h[top], t[-1]))
    print("wrote", OUT_STEM + ".png")


if __name__ == "__main__":
    main()
