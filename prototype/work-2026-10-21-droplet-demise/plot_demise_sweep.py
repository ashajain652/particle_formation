"""Figures for the droplet demise sweep (run_demise_sweep.py) into runs/demise_sweep/figures/.

Every case skipped by the Weber pre-check is shaded grey: a grey cell in the heatmaps, a grey band in the line plots.

  1_final_state.png    final temperature, mass lost, final speed, time to the final state, peak T, Al lost
  2_mass_lost.png      share of the droplet lost to vapour against release altitude, total and by element
  3_histories.png      temperature and remaining mass against time for one diameter at five altitudes
  4_heating_choice.png mass lost with the regime-based heating against PhysicsHeating wherever Mach >= 6
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.patches import Patch, Rectangle

import run_demise_sweep as rs

RUN = sys.argv[1] if len(sys.argv) > 1 else rs.OUT
FIG = os.path.join(RUN, "figures")
os.makedirs(FIG, exist_ok=True)

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SKIP_GREY = "#bdbcb7"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
ORDINAL = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 9.5, "axes.titlesize": 10.5,
    "axes.titleweight": "bold", "legend.frameon": False, "axes.axisbelow": True})

f = pd.read_csv(os.path.join(RUN, "final_states.csv"))
f["skipped"] = f.skipped.astype(str) == "True"
f["mass_lost_pct"] = 100 * f.mass_lost_frac
f["lost_Al_pct"] = 100 * f.lost_Al_frac_of_element
D, H = sorted(f.d0_um.unique()), sorted(f.h0_km.unique())
n_skip = int(f.skipped.sum())
SKIP_LABEL = "skipped: release Weber number at or above {:g} ({} cases)".format(rs.dd.WE_CRIT, n_skip)


def grid(col, frame=f):
    return frame.pivot(index="h0_km", columns="d0_um", values=col).reindex(index=H, columns=D).values


def fmtv(v):
    if not np.isfinite(v):
        return ""
    if v >= 1e4 or (0 < v < 1e-2):
        e = int(np.floor(np.log10(v)))
        return "{:.1f}e{}".format(v / 10 ** e, e)
    return "{:,.0f}".format(v) if v >= 100 else "{:.2g}".format(v)


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 1. final-state heatmaps
panels = [("T_final_K", "Final temperature [K]", False),
          ("mass_lost_pct", "Share of the droplet's mass lost to vapour [%]", False),
          ("V_final_ms", "Final (terminal) speed [m/s]", True),
          ("t_final_s", "Time from release to the final state [s]", True),
          ("T_peak_K", "Peak temperature [K]", False),
          ("lost_Al_pct", "Share of the aluminium lost to vapour [%]", False)]
skip = grid("skipped").astype(bool)
not_steady = (grid("end_reason") != "steady") & ~skip
fig, axs = plt.subplots(3, 2, figsize=(15, 16))
for a, (col, title, log) in zip(axs.flat, panels):
    z = grid(col).astype(float)
    zm = np.ma.masked_where(skip | ~np.isfinite(z), z)
    lo, hi = np.nanmin(zm), np.nanmax(zm)
    norm = LogNorm(lo, hi) if log else Normalize(lo, hi)
    cmap = SEQ.copy(); cmap.set_bad(SKIP_GREY)
    im = a.imshow(zm, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    for j in range(len(D) + 1):
        a.axvline(j - 0.5, color=SURFACE, lw=1.2)
    for i in [H.index(h) for h in (52, 60, 70, 80, 88)]:
        for j in range(len(D)):
            if skip[i, j]:
                continue
            v = z[i, j]
            frac = norm(v)
            a.text(j, i, fmtv(v), ha="center", va="center", fontsize=6.4, color="#ffffff" if frac > 0.55 else INK)
    for i, j in zip(*np.nonzero(not_steady)):
        a.plot(j, i, "x", color=INK, ms=6, mew=1.4)
    a.set_xticks(range(len(D))); a.set_xticklabels(["{:g}".format(d) for d in D], rotation=45, ha="right")
    a.set_yticks(range(0, len(H), 2)); a.set_yticklabels([str(h) for h in H[::2]])
    a.set_xlabel("initial droplet diameter [micrometres]"); a.set_ylabel("release altitude [km]")
    a.grid(False); a.set_title(title, loc="left")
    cb = fig.colorbar(im, ax=a, shrink=0.85, pad=0.02); cb.outline.set_visible(False)
handles = [Patch(fc=SKIP_GREY, label=SKIP_LABEL),
           plt.Line2D([], [], ls="", marker="x", color=INK, mew=1.4,
                      label="run ended before a steady state (ground, full vaporisation or the 1-hour limit; {} cases)".format(
                          int(not_steady.sum())))]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.01, -0.005), fontsize=9)
fig.text(0.01, -0.03, "Final state: speed within 1 % of the local terminal velocity and temperature within 1 K of the local "
         "equilibrium temperature. Numbers printed on the 52, 60, 70, 80 and 88 km rows; every value is in final_states.csv.",
         color=INK2, fontsize=9)
fig.tight_layout(rect=(0, 0.02, 1, 1))
save(fig, "1_final_state.png")

# ---------------------------------------------------------------- 2. mass lost against altitude
sel = [d for d in (1.0, 4.54, 20.6, 64.3, 200.0) if d in D]
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, d in enumerate(sel):
    g = f[f.d0_um == d].sort_values("h0_km")
    for h in g[g.skipped].h0_km:
        ax[0].axvspan(h - 1, h + 1, color=SKIP_GREY, alpha=0.6, lw=0)
    gg = g[~g.skipped]
    ax[0].plot(gg.h0_km, gg.mass_lost_pct, color=ORDINAL[k], label="{:g} micrometres".format(d))
ax[0].set_xlabel("release altitude [km]"); ax[0].set_ylabel("share of initial mass lost to vapour [%]")
ax[0].set_title("Mass lost, by initial diameter", loc="left")
h_, l_ = ax[0].get_legend_handles_labels()
if n_skip:
    h_.append(Patch(fc=SKIP_GREY, label=SKIP_LABEL)); l_.append(SKIP_LABEL)
ax[0].legend(h_, l_, fontsize=8.5, loc="upper left")
g = f[f.d0_um == 20.6].sort_values("h0_km")
for h in g[g.skipped].h0_km:
    ax[1].axvspan(h - 1, h + 1, color=SKIP_GREY, alpha=0.6, lw=0)
gg = g[~g.skipped]
base = np.zeros(len(gg))
for k, el in enumerate(rs.dd.VOLATILE):
    share = 100 * gg["lost_{}_kg".format(el.name)].values / gg.mass0_kg.values
    ax[1].fill_between(gg.h0_km, base, base + share, color=CAT[k], alpha=0.85, lw=0,
                       label="{} ({:.1f} % of the alloy by mass)".format(el.name, 100 * el.mass_fraction))
    base = base + share
ax[1].set_xlabel("release altitude [km]"); ax[1].set_ylabel("share of initial mass lost to vapour [%]")
ax[1].set_title("Mass lost by element, 20.6 micrometre droplet", loc="left")
ax[1].legend(fontsize=8.5, loc="upper left")
save(fig, "2_mass_lost.png")

# ---------------------------------------------------------------- 3. histories
d_hist = 20.6
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, h in enumerate((50, 60, 70, 80, 90)):
    p = os.path.join(RUN, "histories", rs.history_name(d_hist, h))
    if not os.path.exists(p):
        continue
    r = pd.read_csv(p)
    r = r[r.t_s > 0]
    c = CAT[k] if k < 4 else "#e87ba4"
    ax[0].semilogx(r.t_s, r.T_K, color=c, label="released at {} km".format(h))
    ax[1].semilogx(r.t_s, 100 * r.mass_frac, color=c, label="released at {} km".format(h))
ax[0].set_xlabel("time since release [s]"); ax[0].set_ylabel("droplet temperature [K]")
ax[0].set_title("Temperature, {:g} micrometre droplet".format(d_hist), loc="left"); ax[0].legend(fontsize=8.5)
ax[1].set_xlabel("time since release [s]"); ax[1].set_ylabel("remaining mass [% of initial]")
ax[1].set_title("Remaining mass, {:g} micrometre droplet".format(d_hist), loc="left"); ax[1].legend(fontsize=8.5)
save(fig, "3_histories.png")

# ---------------------------------------------------------------- 4. heating choice
pm = os.path.join(RUN, "final_states_mach.csv")
if os.path.exists(pm):
    m = pd.read_csv(pm)
    m["skipped"] = m.skipped.astype(str) == "True"
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 5))
    for k, d in enumerate(sorted(m.d0_um.unique())):
        a_ = f[(f.d0_um == d) & ~f.skipped].sort_values("h0_km")
        b_ = m[(m.d0_um == d) & ~m.skipped].sort_values("h0_km")
        ax.plot(a_.h0_km, 100 * a_.mass_lost_frac, color=CAT[k], label="{:g} um, regime-based heating".format(d))
        ax.plot(b_.h0_km, 100 * b_.mass_lost_frac, color=CAT[k], ls="--", label="{:g} um, validated heating above Mach 6".format(d))
    ax.set_xlabel("release altitude [km]"); ax.set_ylabel("share of initial mass lost to vapour [%]")
    ax.set_title("Effect of the heating hand-over rule", loc="left"); ax.legend(fontsize=8)
    save(fig, "4_heating_choice.png")
