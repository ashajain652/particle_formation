"""Figures for the droplet demise sweep (demise_sweep.py) into output/droplet_reentry_output/demise_sweep/figures/.

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

from reentry_model.droplet_reentry.analysis import demise_sweep as rs

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
# sizes are stored in micrometres; a run whose largest droplet is below 1 um is labelled in nanometres
NANO = max(D) < 1.0
SCALE, UNIT, UNIT_WORD = (1000.0, "nm", "nanometres") if NANO else (1.0, "\u00b5m", "micrometres")


def dlab(d):
    return "{:.3g}".format(d * SCALE)


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
    a.set_xticks(range(len(D))); a.set_xticklabels([dlab(d) for d in D], rotation=45, ha="right")
    a.set_yticks(range(0, len(H), 2)); a.set_yticklabels([str(h) for h in H[::2]])
    a.set_xlabel("initial droplet diameter [{}]".format(UNIT_WORD)); a.set_ylabel("release altitude [km]")
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
sel = [d for d in (1.0, 4.54, 20.6, 64.3, 200.0) if d in D] or [D[int(k)] for k in np.linspace(0, len(D) - 1, 5).round()]
D_MID = 20.6 if 20.6 in D else D[len(D) // 2]
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, d in enumerate(sel):
    g = f[f.d0_um == d].sort_values("h0_km")
    for h in g[g.skipped].h0_km:
        ax[0].axvspan(h - 1, h + 1, color=SKIP_GREY, alpha=0.6, lw=0)
    gg = g[~g.skipped]
    ax[0].plot(gg.h0_km, gg.mass_lost_pct, color=ORDINAL[k], label="{} {}".format(dlab(d), UNIT_WORD))
ax[0].set_xlabel("release altitude [km]"); ax[0].set_ylabel("share of initial mass lost to vapour [%]")
ax[0].set_title("Mass lost, by initial diameter", loc="left")
h_, l_ = ax[0].get_legend_handles_labels()
if n_skip:
    h_.append(Patch(fc=SKIP_GREY, label=SKIP_LABEL)); l_.append(SKIP_LABEL)
ax[0].legend(h_, l_, fontsize=8.5, loc="upper left")
g = f[f.d0_um == D_MID].sort_values("h0_km")
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
ax[1].set_title("Mass lost by element, {} {} droplet".format(dlab(D_MID), UNIT_WORD[:-1]), loc="left")
ax[1].legend(fontsize=8.5, loc="upper left")
save(fig, "2_mass_lost.png")

# ---------------------------------------------------------------- 3. histories
d_hist = D_MID
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
ax[0].set_title("Temperature, {} {} droplet".format(dlab(d_hist), UNIT_WORD[:-1]), loc="left"); ax[0].legend(fontsize=8.5)
ax[1].set_xlabel("time since release [s]"); ax[1].set_ylabel("remaining mass [% of initial]")
ax[1].set_title("Remaining mass, {} {} droplet".format(dlab(d_hist), UNIT_WORD[:-1]), loc="left"); ax[1].legend(fontsize=8.5)
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

# ---------------------------------------------------------------- 5. final terminal velocity, labelled with final mass, size, altitude
# Every second release altitude, every cell labelled in large type; the full 21-row map is figure 1's speed panel.
H5 = [h for h in H if (h - 50) % 4 == 0]
sub = f[f.h0_km.isin(H5)]
def grid5(col):
    return sub.pivot(index="h0_km", columns="d0_um", values=col).reindex(index=H5, columns=D).values
z = grid5("v_terminal_final_ms").astype(float)
sk = grid5("skipped").astype(bool)
zm = np.ma.masked_where(sk | ~np.isfinite(z), z)
norm = LogNorm(np.nanmin(zm), np.nanmax(zm))
cmap = SEQ.copy(); cmap.set_bad(SKIP_GREY)
dfin, mfin, hfin = (grid5(c).astype(float) for c in ("d_final_um", "mass_final_kg", "h_final_km"))
tfin = grid5("t_final_s").astype(float)


def tlab(t):
    """Time to the final state: seconds, in exponent form below 0.01 s."""
    if t < 0.01:
        m, e = "{:.1e}".format(t).split("e")        # exponent of the rounded value, so 0.00999 -> 1.0e-2
        return "{}e{} s".format(m, int(e))
    return "{:.3g} s".format(t)


fig, a = plt.subplots(figsize=(30, 17))
im = a.imshow(zm, origin="lower", aspect="auto", cmap=cmap, norm=norm)
for j in range(len(D) + 1):
    a.axvline(j - 0.5, color=SURFACE, lw=2)
for i in range(len(H5) + 1):
    a.axhline(i - 0.5, color=SURFACE, lw=2)
for i in range(len(H5)):
    for j in range(len(D)):
        if sk[i, j]:
            continue
        dark = norm(z[i, j]) > 0.5
        a.text(j, i, "{} m/s\nat {:.1f} km\n{:.1e} kg\nafter {}".format(fmtv(z[i, j]), hfin[i, j], mfin[i, j],
                                                                       tlab(tfin[i, j])),
               ha="center", va="center", fontsize=15, linespacing=1.2, color="#ffffff" if dark else INK)
a.set_xticks(range(len(D))); a.set_xticklabels([dlab(d) for d in D], fontsize=18)
a.set_yticks(range(len(H5))); a.set_yticklabels([str(h) for h in H5], fontsize=18)
a.set_xlabel("initial droplet diameter [{}]".format(UNIT_WORD), fontsize=20, labelpad=10)
a.set_ylabel("release (starting) altitude [km]", fontsize=20, labelpad=10); a.grid(False)
a.set_title("Final state of each droplet (Singh et al. drag): terminal velocity, altitude reached, final mass, time from release",
            loc="left", fontsize=22, pad=14)
cb = fig.colorbar(im, ax=a, shrink=0.8, pad=0.012); cb.outline.set_visible(False)
cb.set_label("terminal velocity at the final state [m/s]", color=INK2, fontsize=18); cb.ax.tick_params(labelsize=16)
fig.text(0.01, 0.005, "Rows are release altitudes every 4 km; the full 2 km grid is in final_states.csv and in figure 1. "
         "Grey cells would mark cases skipped by the Weber check (none).", color=INK2, fontsize=16)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(FIG, "5_terminal_velocity.png"), dpi=110, bbox_inches="tight")
plt.close(fig)
print("wrote 5_terminal_velocity.png")

# ---------------------------------------------------------------- 6. total mass lost to vaporisation (absolute), same layout as 5
z = grid5("mass_lost_kg").astype(float)
zm = np.ma.masked_where(sk | ~np.isfinite(z), z)
norm = LogNorm(np.nanmin(zm), np.nanmax(zm))
frac, m0 = grid5("mass_lost_frac").astype(float), grid5("mass0_kg").astype(float)
fig, a = plt.subplots(figsize=(30, 17))
im = a.imshow(zm, origin="lower", aspect="auto", cmap=cmap, norm=norm)
for j in range(len(D) + 1):
    a.axvline(j - 0.5, color=SURFACE, lw=2)
for i in range(len(H5) + 1):
    a.axhline(i - 0.5, color=SURFACE, lw=2)
for i in range(len(H5)):
    for j in range(len(D)):
        if sk[i, j]:
            continue
        dark = norm(z[i, j]) > 0.5
        a.text(j, i, "{:.1f} %\nafter {}".format(100 * frac[i, j], tlab(tfin[i, j])),
               ha="center", va="center", fontsize=15, linespacing=1.3, color="#ffffff" if dark else INK)
a.set_xticks(range(len(D))); a.set_xticklabels([dlab(d) for d in D], fontsize=18)
a.set_yticks(range(len(H5))); a.set_yticklabels([str(h) for h in H5], fontsize=18)
a.set_xlabel("initial droplet diameter [{}]".format(UNIT_WORD), fontsize=20, labelpad=10)
a.set_ylabel("release (starting) altitude [km]", fontsize=20, labelpad=10); a.grid(False)
a.set_title("Mass lost to vaporisation (Singh et al. drag): colour in kg; labels give % of initial mass and time to the final state",
            loc="left", fontsize=22, pad=14)
cb = fig.colorbar(im, ax=a, shrink=0.8, pad=0.012); cb.outline.set_visible(False)
cb.set_label("mass lost to vaporisation [kg]", color=INK2, fontsize=18); cb.ax.tick_params(labelsize=16)
fig.text(0.01, 0.005, "Aluminium, zinc and magnesium together; cell labels are the share of the initial mass lost [%]. Rows are release altitudes every 4 km; the full 2 km grid is in "
         "final_states.csv. Grey cells would mark cases skipped by the Weber check (none).", color=INK2, fontsize=16)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(FIG, "6_mass_lost_total.png"), dpi=110, bbox_inches="tight")
plt.close(fig)
print("wrote 6_mass_lost_total.png")
