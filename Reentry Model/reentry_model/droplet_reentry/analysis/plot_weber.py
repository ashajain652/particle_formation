"""Weber-number stability of every droplet in the demise sweeps (1-500 nm and 1-200 um, Singh drag), into
output/droplet_reentry_output/demise_sweep/figures/10_weber_stability.png.

The Weber number We = rho V^2 d / sigma(T) is followed along the whole flight with the shrinking diameter and a
surface tension that falls with temperature; the break-up threshold is droplet_demise.WE_CRIT (12, Pilch & Erdman
1987 for low Ohnesorge number). Vertical winds do not change it: the maximum is at release, at km/s airspeeds.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

from reentry_model.droplet_reentry import model as dd
from reentry_model.droplet_reentry.analysis import demise_sweep as rs

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
CRITICAL = "#d03b3b"
SEQ = LinearSegmentedColormap.from_list("blue", ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#184f95", "#0d366b"])
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 1.8, "font.size": 14, "axes.titlesize": 16,
    "axes.titleweight": "bold", "legend.frameon": False, "axes.axisbelow": True})

f = pd.concat([pd.read_csv(os.path.join(b, "final_states.csv")) for b in (rs.OUT_NANO, rs.OUT)]).sort_values("d0_um").reset_index(drop=True)
f["skipped"] = f.skipped.astype(str) == "True"
n, n_skip = len(f), int(f.skipped.sum())
we_max = f.We_max.max()
worst = f.loc[f.We_max.idxmax()]

fig, ax = plt.subplots(figsize=(15, 9))
ax.axhspan(dd.WE_CRIT, 1e3, color=CRITICAL, alpha=0.08, lw=0)
ax.axhline(dd.WE_CRIT, color=CRITICAL, lw=1.6, ls="--")
ax.text(1.2, dd.WE_CRIT * 1.25, "break-up threshold, We = {:g} (Pilch & Erdman 1987): droplets above it would be "
        "skipped and shaded grey".format(dd.WE_CRIT), fontsize=13, color=INK2)
norm = Normalize(50, 90)
for h, g in f.groupby("h0_km"):
    g = g[~g.skipped]
    ax.loglog(g.d0_um * 1000.0, g.We_max, color=SEQ(norm(h)))
ax.plot(worst.d0_um * 1000.0, worst.We_max, "o", ms=10, color=INK, mec=SURFACE, mew=2)
ax.annotate("highest anywhere: We = {:.2f}\n({:g} um released at {:g} km), {:.0f} times below the threshold".format(
    we_max, worst.d0_um, worst.h0_km, dd.WE_CRIT / we_max), (worst.d0_um * 1000.0, worst.We_max), xytext=(40, 0.15),
    textcoords="data", fontsize=13, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED))
ax.axvline(1000.0, color=MUTED, lw=0.8, ls=":")
ax.text(1100, 2e-8, "1 um", fontsize=12, color=INK2)
ax.set_xlim(0.8, 3e5); ax.set_ylim(1e-8, 300)
ax.set_xlabel("initial droplet diameter [nanometres]")
ax.set_ylabel("maximum Weber number along the flight")
ax.set_title("Every droplet is stable against aerodynamic break-up: {} of {} cases below the threshold, {} skipped".format(
    n - n_skip - int((f.We_max >= dd.WE_CRIT).sum()), n, n_skip), loc="left")
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=SEQ), ax=ax, pad=0.015)
cb.outline.set_visible(False); cb.set_label("release altitude [km]", color=INK2)
fig.text(0.01, 0.005, "We = air density x speed^2 x diameter / surface tension, followed along each flight with the "
         "shrinking diameter and temperature-dependent surface tension; it peaks at release and only falls afterwards. "
         "One line per release altitude (21), 1 nm to 200 um, Singh et al. drag.", color=INK2, fontsize=11)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(rs.OUT, "figures", "10_weber_stability.png"), dpi=130, bbox_inches="tight")
print("wrote 10_weber_stability.png; max We {:.3f}; cases {}; skipped {}; We_max/We_release max {:.4f}".format(
    we_max, n, n_skip, (f.We_max / f.We_release).max()))
