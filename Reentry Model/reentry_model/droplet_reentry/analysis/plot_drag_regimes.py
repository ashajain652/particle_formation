"""The Singh et al. (2021) sphere drag coefficient (singh_drag.py) across the flow regimes, with the paths of three
droplets from the demise sweep, into output/droplet_reentry_output/demise_sweep/figures/7_drag_law_regimes.png.

(a) C_D over the Knudsen number (lambda/D) and Mach number plane; the two fix the Reynolds number,
    Re = (Ma / Kn) sqrt(pi gamma / 2), so the plane covers every regime. Wall at the air temperature.
(b) C_D against Mach number at fixed Knudsen numbers.  (c) C_D against Knudsen number at fixed Mach numbers.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm

from reentry_model.droplet_reentry.analysis import demise_sweep as rs
from reentry_model.droplet_reentry import singh_drag as sg

FIG = os.path.join(rs.OUT, "figures")
os.makedirs(FIG, exist_ok=True)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
RAMP6 = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#184f95", "#0d366b"]
PATH_COLOURS = ["#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 12, "axes.titlesize": 13.5,
    "axes.titleweight": "bold", "legend.frameon": False, "axes.axisbelow": True})
G = sg.GAMMA
RE_MAX = 1e5          # above this the correlation's data stop and turbulent transition is not modelled


def cd(kn, ma):
    re = (ma / kn) * math.sqrt(math.pi * G / 2.0)
    return sg.drag_coefficient(ma, re, 1.0)[0], re


fig = plt.figure(figsize=(20, 14))
gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.32, wspace=0.22)
ax = fig.add_subplot(gs[:, 0])

# ---------------------------------------------------------------- (a) regime map
KN = np.logspace(-3, 5, 240)
MA = np.logspace(-3, np.log10(30.0), 220)
Z = np.empty((MA.size, KN.size)); RE = np.empty_like(Z)
for i, m in enumerate(MA):
    for j, k in enumerate(KN):
        Z[i, j], RE[i, j] = cd(k, m)
assert not (RE > RE_MAX).any()      # Re = (Ma/Kn) sqrt(pi gamma/2) stays below 1e5 everywhere on this plane
pc = ax.pcolormesh(KN, MA, Z, cmap=SEQ, norm=LogNorm(0.5, Z.max()), shading="auto", rasterized=True)
cs = ax.contour(KN, MA, Z, levels=[1, 2, 10, 100, 1000], colors=INK, linewidths=0.9, alpha=0.6)
ax.clabel(cs, fmt=lambda v: "C_D = {:g}".format(v), fontsize=10.5, inline_spacing=6)
for k, lab in ((0.01, "continuum | slip"), (0.1, "slip | transitional"), (10.0, "transitional | free-molecule")):
    ax.axvline(k, color=INK, lw=1.1, ls="--", alpha=0.7)
    ax.text(k * 1.12, 25, lab, rotation=90, fontsize=10.5, color=INK2, va="top")
for m, lab in ((1.0, "Mach 1"), (5.0, "Mach 5 (hypersonic)")):
    ax.axhline(m, color=INK, lw=1.1, ls=":", alpha=0.7)
    ax.text(1.3e-3, m * 1.08, lab, fontsize=10.5, color=INK2, ha="left")
# droplet paths from the sweep (histories store Kn as 2 lambda / D)
paths = [(1.0, 90), (20.6, 70), (200.0, 50)]
for (d, h), c in zip(paths, PATH_COLOURS):
    p = os.path.join(rs.OUT, "histories", rs.history_name(d, h))
    if not os.path.exists(p):
        continue
    r = pd.read_csv(p)
    kn, ma = r.Kn.values / 2.0, r.Ma.values
    keep = (ma >= MA[0]) & (ma <= MA[-1]) & (kn >= KN[0]) & (kn <= KN[-1])
    ax.plot(kn[keep], ma[keep], color=SURFACE, lw=5.5, solid_capstyle="round")
    ax.plot(kn[keep], ma[keep], color=c, lw=2.8, label="{:g} µm droplet released at {} km".format(d, h))
    ax.plot(kn[0], ma[0], "o", ms=11, color=c, mec=SURFACE, mew=2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(KN[0], KN[-1]); ax.set_ylim(MA[0], MA[-1])
ax.set_xlabel("Knudsen number (air mean free path / droplet diameter)")
ax.set_ylabel("Mach number")
ax.set_title("(a) Drag coefficient across the regimes, with three droplets' paths (dot = release)", loc="left")
ax.legend(loc="upper left", fontsize=11, bbox_to_anchor=(0.0, -0.07), ncol=3, columnspacing=1.2)
ax.grid(False)
cb = fig.colorbar(pc, ax=ax, shrink=0.75, pad=0.015); cb.outline.set_visible(False)
cb.set_label("drag coefficient C_D", color=INK2)

# ---------------------------------------------------------------- (b) against Mach at fixed Knudsen
ax2 = fig.add_subplot(gs[0, 1])
ma = np.logspace(-1, np.log10(30.0), 300)
for k, c in zip((1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0), RAMP6):
    v = np.array([cd(k, m) for m in ma])
    good = v[:, 1] <= RE_MAX
    ax2.plot(ma[good], v[good, 0], color=c, label="{:g}".format(k))
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlim(0.1, 30)
ax2.legend(title="Knudsen number", fontsize=10.5, title_fontsize=10.5, loc="upper right", ncol=2)
ax2.axvline(1.0, color=MUTED, lw=1, ls=":")
ax2.set_xlabel("Mach number"); ax2.set_ylabel("drag coefficient C_D")
ax2.set_title("(b) Against Mach number, at fixed Knudsen number", loc="left")

# ---------------------------------------------------------------- (c) against Knudsen at fixed Mach
ax3 = fig.add_subplot(gs[1, 1])
kn = np.logspace(-3, 5, 300)
for m, c in zip((0.1, 0.5, 1.0, 2.0, 5.0, 25.0), RAMP6):
    v = np.array([cd(k, m) for k in kn])
    good = v[:, 1] <= RE_MAX
    ax3.plot(kn[good], v[good, 0], color=c)
    ax3.annotate("Mach {:g}".format(m), (kn[good][-1], v[good, 0][-1]), xytext=(6, 0), textcoords="offset points",
                 fontsize=10.5, color=INK2, va="center")
for k in (0.01, 0.1, 10.0):
    ax3.axvline(k, color=MUTED, lw=1, ls="--")
ax3.set_xscale("log"); ax3.set_yscale("log")
ax3.set_xlim(1e-3, 1.5e6)
ax3.set_xlabel("Knudsen number (air mean free path / droplet diameter)"); ax3.set_ylabel("drag coefficient C_D")
ax3.set_title("(c) Against Knudsen number, at fixed Mach number", loc="left")

fig.text(0.01, 0.01, "Singh et al. (2021), AIAA Journal, DOI 10.2514/1.J060648, as implemented in singh_drag.py; droplet at "
         "the air temperature, gamma = 1.4. At fixed Mach number, a larger Knudsen number means a smaller Reynolds number, "
         "Re = (Ma/Kn) sqrt(pi gamma / 2), which is why C_D rises at low speed.", color=INK2, fontsize=11)
fig.savefig(os.path.join(FIG, "7_drag_law_regimes.png"), dpi=120, bbox_inches="tight")
print("wrote 7_drag_law_regimes.png")
