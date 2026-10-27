"""Mass lost to vaporisation against altitude, for each release altitude of the 2026-10-22 demise sweep.

A droplet loses all the mass it will lose while braking from the parent's speed, which takes it only metres to a few
hundred metres below its release altitude; after that it is near the air temperature and settles without losing more.
So on a kilometre scale the ablation profile is a spike at the release altitude, and the left panel gives its size:
the share of the initial mass lost against release altitude. The right panel zooms into the spike: the cumulative
share lost against the depth below release. In free-molecule flow that depth is proportional to the diameter (the
braking length scales as rho_p d / rho_air) while the share lost does not depend on it (plot_size_independence.py), so
the depth is drawn for a 10 um droplet and 1 um and 100 um droplets, rescaled by 10/d, are overlaid to show the
curves coincide. Each case is re-run with dd.run_case, since final_states.csv keeps no trajectories.

  output/droplet_reentry_output/demise_sweep/figures/12_ablation_profile.{png,pdf,svg}
  output/droplet_reentry_output/demise_sweep/ablation_profiles.csv
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

from reentry_model.droplet_reentry import model as dd
from reentry_model.droplet_reentry.analysis import demise_sweep as rs

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
ACCENT = "#2a78d6"
RAMP = LinearSegmentedColormap.from_list("blues", ["#b9d4f5", "#5598e7", "#1c5cab", "#0d366b"])
CHECK_SIZES = (1.0, 100.0)                 # um, overlaid to show size independence
CHECK_ALTITUDES = (90, 78, 66, 54)          # km
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 13, "axes.axisbelow": True})


def profile(d_um, h0_km):
    """Depth below the release altitude (m) and cumulative share of the initial mass lost (%) along one flight."""
    _, rows = dd.run_case(d_um, h0_km)
    h = np.array([r["h_km"] for r in rows])
    m = np.array([r["mass_kg"] for r in rows])
    return 1e3 * (h0_km - h), 100.0 * (1.0 - m / m[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--d", type=float, default=10.0, help="diameter of the main curves, um")
    ap.add_argument("--outdir", default=rs.OUT)
    a = ap.parse_args()
    fig_dir = os.path.join(a.outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    curves = {(a.d, h): profile(a.d, float(h)) for h in rs.ALTITUDES}
    checks = {(d, h): profile(d, float(h)) for d in CHECK_SIZES for h in CHECK_ALTITUDES}
    pd.concat([pd.DataFrame({"d0_um": d, "h0_km": h0, "depth_below_release_m": z, "mass_lost_pct": lost})
               for (d, h0), (z, lost) in {**curves, **checks}.items()]
              ).to_csv(os.path.join(a.outdir, "ablation_profiles.csv"), index=False)

    fig, (L, R) = plt.subplots(1, 2, figsize=(13.333, 7.5), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.3})
    fig.subplots_adjust(left=0.07, right=0.88, top=0.83, bottom=0.1)
    fig.text(0.07, 0.94, "Droplets vaporise within metres of where they are released",
             fontsize=20, weight="bold", color=INK)
    fig.text(0.07, 0.895, "AA7075 droplets released from the parent at its speed; the share lost is the same for "
             "every size in free-molecule flow", fontsize=13, color=INK2)

    # ---- left: how much is lost at each release altitude
    h0s = np.array(rs.ALTITUDES, dtype=float)
    lost = np.array([curves[(a.d, h)][1][-1] for h in rs.ALTITUDES])
    L.barh(h0s, lost, height=1.6, color=ACCENT, edgecolor=SURFACE, linewidth=2)
    for h, v in zip(h0s, lost):
        if h % 4 == 2:
            L.text(v + 0.8, h, "{:.0f} %".format(v), va="center", fontsize=10.5, color=INK2)
    L.set_xlabel("Share of initial mass vaporised (%)")
    L.set_ylabel("Release altitude (km)")
    L.set_xlim(0, 55)
    L.set_ylim(48.5, 91.5)
    L.grid(axis="y", visible=False)
    L.set_title("All of it is lost at the release altitude", loc="left", fontsize=13.5, color=INK)

    # ---- right: the spike resolved in metres
    norm = Normalize(min(rs.ALTITUDES) - 6, max(rs.ALTITUDES))
    for (d, h0), (z, cum) in curves.items():
        if h0 % 4 == 2:                         # every other release altitude, to keep the panel readable
            R.plot(cum, z, color=RAMP(norm(h0)))
    styles = {CHECK_SIZES[0]: (0, (1, 2)), CHECK_SIZES[1]: (0, (5, 3))}
    for (d, h0), (z, cum) in checks.items():
        R.plot(cum, z * a.d / d, color=INK, lw=1.0, ls=styles[d])
    for d in CHECK_SIZES:
        R.plot([], [], color=INK, lw=1.0, ls=styles[d], label="{:g} µm, depth × {:g}".format(d, a.d / d))
    R.legend(title="Same release, other sizes", loc="lower left", frameon=False, fontsize=11, title_fontsize=11)
    R.set_ylim(15, 0)
    R.set_xlim(0, 50)
    R.set_xlabel("Cumulative share of initial mass vaporised (%)")
    R.set_ylabel("Depth below release altitude (m), {:g} µm droplet".format(a.d))
    R.set_title("Resolved in metres (depth scales with diameter)", loc="left", fontsize=13.5,
                color=INK)
    cax = fig.add_axes([0.9, 0.1, 0.015, 0.73])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=RAMP), cax=cax)
    cb.set_label("Release altitude (km)", color=INK2)
    cb.outline.set_visible(False)

    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(fig_dir, "12_ablation_profile." + ext), dpi=200)
    print(os.path.join(fig_dir, "12_ablation_profile.png"))
    for (d, h0), (z, cum) in {**curves, **checks}.items():
        i = np.searchsorted(cum, [0.1 * cum[-1], 0.9 * cum[-1]])
        print("{:5g} um {:3d} km: lost {:5.1f} %, 10-90 % of it between {:.2f} and {:.2f} m below release".format(
            d, h0, cum[-1], *z[i]))


if __name__ == "__main__":
    main()
