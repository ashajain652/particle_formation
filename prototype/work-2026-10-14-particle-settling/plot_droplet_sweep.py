"""Figures for the droplet heating sweep (droplet_heating_sweep.py, summarise_droplet_sweep.py).

Reads runs/droplet_sweep/ and writes six PNGs into runs/droplet_sweep/figures/:

  1_temperature.png     peak temperature against release altitude (with sensitivity) and 10 um histories
  2_validity_map.png    which cases are physical, over the whole diameter x altitude matrix
  3_energy.png          energy budget per kilogram and the vaporised-fraction bracket
  4_matrix_heatmaps.png deceleration, stopping time, stopping distance and time to the ground
  5_trajectories.png    speed and path of five diameters released at 80 km
  6_regime_map.png      trajectories in the Knudsen-number / speed-ratio plane against the blend corner
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm
from matplotlib.patches import Patch, Rectangle

import droplet_heating_sweep as dh

RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(dh.HERE, "runs", "droplet_sweep")
FIG = os.path.join(RUN, "figures")
os.makedirs(FIG, exist_ok=True)

# reference palette (dataviz skill, references/palette.md), light mode
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # categorical slots 1-4, fixed order
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORDINAL = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]   # from step 250, light -> dark
GOOD, SERIOUS, CRITICAL = "#0ca30c", "#ec835a", "#d03b3b"
SEQ = LinearSegmentedColormap.from_list("blue", BLUE)

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
    "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "bold", "legend.frameon": False,
    "axes.axisbelow": True,
})

b = pd.read_csv(os.path.join(RUN, "summary_eps0.1_alpha1_with_vapour.csv"))
variants = {f: pd.read_csv(os.path.join(RUN, f)) for f in sorted(os.listdir(RUN))
            if f.startswith("summary_eps") and f.endswith("_with_vapour.csv")}
D = sorted(b.d_um.unique())
H = sorted(b.h0_km.unique())


def hist(d, h):
    return pd.read_csv(os.path.join(RUN, "histories", dh.history_name(d, h)))


def grid(col, frame=b):
    return frame.pivot(index="h0_km", columns="d_um", values=col).loc[H, D].values


def matrix_axes(ax):
    ax.set_xticks(range(len(D)))
    ax.set_xticklabels(["{:g}".format(d) for d in D], rotation=45, ha="right")
    ax.set_yticks(range(0, len(H), 2))
    ax.set_yticklabels([str(h) for h in H[::2]])
    ax.set_xlabel("droplet diameter [micrometres]")
    ax.set_ylabel("release altitude [km]")
    ax.grid(False)
    ax.set_xlim(-0.5, len(D) - 0.5); ax.set_ylim(-0.5, len(H) - 0.5)


def fmtv(v):
    """Compact cell label: thousands separators in the middle range, mantissa-exponent outside it."""
    if v >= 1e4 or v < 1e-3:
        e = int(np.floor(np.log10(v)))
        return "{:.1f}e{}".format(v / 10 ** e, e)
    return "{:,.0f}".format(v) if v >= 100 else "{:.2g}".format(v)


def hatch_cells(ax, mask, hatch="////", color=INK):
    for i, j in zip(*np.nonzero(mask)):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch=hatch, edgecolor=color,
                               linewidth=0, alpha=0.55))


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 1. temperature
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
a = ax[0]
ref = b[b.d_um == 10].sort_values("h0_km")
a.axhspan(dh.T_BOIL_AL, 4300, color=CRITICAL, alpha=0.06, lw=0)
a.axhline(dh.T_BOIL_AL, color=CRITICAL, lw=1.2, ls=":")
a.text(50.3, dh.T_BOIL_AL - 60, "aluminium boiling\npoint at 1 atm,\n2743 K", color=INK2, fontsize=8.5, va="top")
a.plot(ref.h0_km, ref.T_sat_Al_release_K, color=INK2, lw=1.4, ls="--")
a.text(90.8, 1330, "aluminium boils at the\nlocal air pressure", color=INK2, fontsize=8.5, ha="right", va="top")
a.axhline(dh.T_RELEASE, color=MUTED, lw=1, ls="-.")
a.text(90.8, dh.T_RELEASE + 25, "release temperature, 908 K", color=INK2, fontsize=8.5, ha="right")
order = [("summary_eps0.1_alpha1_with_vapour.csv", "baseline: emissivity 0.1, accommodation 1"),
         ("summary_eps0.05_alpha1_with_vapour.csv", "emissivity 0.05"),
         ("summary_eps0.4_alpha1_with_vapour.csv", "emissivity 0.4"),
         ("summary_eps0.1_alpha0.8_with_vapour.csv", "accommodation 0.8")]
for k, (f, lab) in enumerate(order):
    g = variants[f][variants[f].d_um == 10].sort_values("h0_km")
    a.plot(g.h0_km, g.T_peak_K, color=CAT[k], lw=2.6 if k == 0 else 1.6, marker="o", ms=4 if k == 0 else 0,
           label=lab, zorder=3 if k == 0 else 2)
    pk = g.loc[g.T_peak_K.idxmax()]
    a.annotate("{:.0f} K".format(pk.T_peak_K), (pk.h0_km, pk.T_peak_K), xytext=(6, 4),
               textcoords="offset points", fontsize=8, color=INK2)
a.set_xlim(49, 91); a.set_ylim(800, 4300)
a.set_xlabel("release altitude [km]"); a.set_ylabel("peak droplet temperature [K]")
a.set_title("Peak temperature, all diameters (they agree to within 1.4 %)", loc="left")
a.legend(loc="upper left", fontsize=8.5)
a = ax[1]
for k, h in enumerate((50, 60, 70, 80, 90)):
    r = hist(10, h)
    a.semilogx(r.t_s.clip(lower=1e-7), r.T_K, color=CAT[k] if k < 4 else "#e87ba4", label="released at {} km".format(h))
a.axhline(dh.T_BOIL_AL, color=CRITICAL, lw=1.2, ls=":")
a.text(90, dh.T_BOIL_AL + 30, "boiling point at 1 atm", color=INK2, fontsize=8.5, ha="right")
a.axhline(dh.T_RELEASE, color=MUTED, lw=1, ls="-.")
a.text(90, dh.T_RELEASE + 30, "release temperature", color=INK2, fontsize=8.5, ha="right")
a.set_xlim(1e-5, 1e2)
a.set_xlabel("time since release [s]"); a.set_ylabel("droplet temperature [K]")
a.set_title("Temperature history of a 10 micrometre droplet (baseline)", loc="left")
a.legend(loc="upper left", fontsize=8.5)
fig.text(0.01, -0.02, "No vaporisation is modelled: wherever a curve sits above the dashed line (left panel), the temperature is "
         "not a prediction but a sign that the droplet would boil.", color=INK2, fontsize=9)
save(fig, "1_temperature.png")

# ---------------------------------------------------------------- 2. validity map
status = np.where(grid("flag_T_above_Al_boiling_1atm").astype(bool), 2,
                  np.where(grid("flag_T_above_Al_saturation_local").astype(bool), 1, 0))
fig, ax = plt.subplots(figsize=(12, 7.5))
ax.imshow(status, origin="lower", aspect="auto", cmap=ListedColormap([GOOD, SERIOUS, CRITICAL]),
          vmin=-0.5, vmax=2.5, alpha=0.85)
for j in range(len(D) + 1):
    ax.axvline(j - 0.5, color=SURFACE, lw=2)
for i in range(len(H) + 1):
    ax.axhline(i - 0.5, color=SURFACE, lw=2)
hatch_cells(ax, grid("blend_warning").astype(bool))
fm = grid("frac_E_in_Kn_lt10") > 0.5
for i, j in zip(*np.nonzero(fm)):
    ax.plot(j, i, "o", ms=6, mfc=SURFACE, mec=INK, mew=1.2)
for i, j in zip(*np.nonzero(grid("flag_curvature_matters").astype(bool))):
    ax.plot(j, i, "^", ms=7, mfc=INK, mec=SURFACE, mew=1)
matrix_axes(ax)
ax.set_title("Which of the {} baseline cases can be read as predictions".format(len(b)), loc="left")
n = [(status == k).sum() for k in range(3)]
handles = [Patch(fc=GOOD, alpha=0.85, label="below aluminium's local boiling point ({} cases)".format(n[0])),
           Patch(fc=SERIOUS, alpha=0.85, label="above local boiling point, below 2743 K ({} cases)".format(n[1])),
           Patch(fc=CRITICAL, alpha=0.85, label="above 2743 K, the 1-atm boiling point ({} cases)".format(n[2])),
           Patch(fc="none", ec=INK, hatch="////", label="drag blend-domain warning ({} cases)".format(
               int(grid("blend_warning").sum()))),
           plt.Line2D([], [], ls="", marker="o", ms=6, mfc=SURFACE, mec=INK, mew=1.2,
                      label="most heat arrives at Knudsen number below 10: heating overestimated ({})".format(
                          int(fm.sum()))),
           plt.Line2D([], [], ls="", marker="^", ms=7, mfc=INK, mec=SURFACE,
                      label="Earth curvature exceeds 300 m at the stopping point ({})".format(
                          int(grid("flag_curvature_matters").sum())))]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.09), fontsize=8.5, ncol=1)
save(fig, "2_validity_map.png")

# ---------------------------------------------------------------- 3. energy
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
a = ax[0]
g = b[b.d_um == 10].sort_values("h0_km")
ke = 0.5 * g.V0_ms ** 2 / 1e6
a.plot(g.h0_km, ke, color=CAT[0], label="initial kinetic energy")
a.plot(g.h0_km, g.E_in_Jkg / 1e6, color=CAT[1], label="energy absorbed from the air")
a.plot(g.h0_km, g.E_rad_Jkg / 1e6, color=CAT[2], label="energy radiated away")
budget = (dh.CP_LIQUID * (dh.T_BOIL_AL - dh.T_RELEASE) + dh.H_VAP_AL) / 1e6
a.axhline(budget, color=CRITICAL, ls=":", lw=1.2)
a.text(50.3, budget + 0.5, "heat to boiling and vaporise completely, {:.1f} MJ/kg".format(budget), color=INK2, fontsize=8.5)
a.text(89, 6.5, "absorbed and radiated nearly coincide above 70 km:\nradiation removes the heat almost as fast as it arrives",
       color=INK2, fontsize=8.5, ha="right")
a.set_xlim(49, 91)
a.set_xlabel("release altitude [km]"); a.set_ylabel("energy per kilogram of droplet [MJ/kg]")
a.set_title("Energy budget per kilogram (10 micrometres; all sizes alike)", loc="left")
a.legend(loc="upper left", fontsize=8.5, bbox_to_anchor=(0, 0.93))
a = ax[1]
a.fill_between(g.h0_km, 100 * g.vap_frac_lower, 100 * g.vap_frac_upper, color=CAT[1], alpha=0.25, lw=0,
               label="range between the two bounds")
a.plot(g.h0_km, 100 * g.vap_frac_upper, color=CAT[1], lw=1.6,
       label="upper bound: all absorbed heat beyond boiling makes vapour")
a.plot(g.h0_km, 100 * g.vap_frac_lower, color=CAT[1], lw=1.6, ls="--",
       label="lower bound: only the overshoot above boiling makes vapour")
a.set_xlim(49, 91); a.set_ylim(0, 105)
a.set_xlabel("release altitude [km]"); a.set_ylabel("share of droplet mass vaporised [%]")
a.set_title("Rough bracket on the vaporised mass (not a model)", loc="left")
a.legend(loc="upper left", fontsize=8.5)
save(fig, "3_energy.png")

# ---------------------------------------------------------------- 4. matrix heatmaps
panels = [("peak_decel_g", "Peak deceleration [multiples of Earth gravity]", "{:.0f}"),
          ("t_stop_s", "Time to lose 99 % of horizontal speed [s]", "{:.2g}"),
          ("x_stop_km", "Downrange distance at that point [km]", "{:.2g}"),
          ("t_ground_days", "Time to reach the ground [days]", "{:.2g}"),
          ("v_terminal_ms", "Terminal velocity where it is first reached [m/s]", "{:.2g}"),
          ("t_terminal_s", "Time from release to terminal velocity (within 1 %) [s]", "{:.2g}")]
term = pd.read_csv(os.path.join(RUN, "terminal_approach.csv"))
b = b.merge(term[["d_um", "h0_km", "v_terminal_ms", "t_terminal_s", "h_terminal_km"]], on=["d_um", "h0_km"])
fig, axs = plt.subplots(3, 2, figsize=(19, 18))
for a, (col, title, fmt) in zip(axs.flat, panels):
    z = grid(col, b)
    im = a.imshow(z, origin="lower", aspect="auto", cmap=SEQ, norm=LogNorm(z.min(), z.max()))
    for j in range(len(D) + 1):
        a.axvline(j - 0.5, color=SURFACE, lw=1.5)
    hatch_cells(a, grid("blend_warning").astype(bool), color=SURFACE)
    # label every other altitude row so the grid stays readable; the CSV has every value
    for i in [H.index(h) for h in (52, 60, 70, 80, 88)]:
        for j in range(len(D)):
            v = z[i, j]
            dark = np.log(v / z.min()) / np.log(z.max() / z.min()) > 0.55
            a.text(j, i, fmtv(v), ha="center", va="center", fontsize=6.6,
                   color="#ffffff" if dark else INK)
    matrix_axes(a)
    a.set_title(title, loc="left")
    cb = fig.colorbar(im, ax=a, shrink=0.85, pad=0.02)
    cb.outline.set_visible(False)
fig.text(0.01, 0.003, "Hatched cells carry the drag blend-domain warning. Numbers are printed on the 52, 60, 70, 80 and 88 km rows; every value is "
         "in matrix_summary_baseline.csv and terminal_approach.csv.\nTerminal velocity is the local settling speed, which keeps "
         "falling as the droplet sinks; 'reached' means the speed stays within 1 % of it from then on.",
         color=INK2, fontsize=9)
fig.tight_layout(rect=(0, 0.025, 1, 1))
save(fig, "4_matrix_heatmaps.png")

# ---------------------------------------------------------------- 5. trajectories at 80 km
sel = (0.01, 0.141, 2, 18, 250)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, d in enumerate(sel):
    r = hist(d, 80)
    r = r[r.t_s > 0]
    ax[0].loglog(r.t_s, r.V, color=ORDINAL[k], label="{} micrometres".format(d))
    p = r[r.h_m >= 60e3]
    ax[1].semilogx(p.x_m.clip(lower=1) / 1e3, p.h_m / 1e3, color=ORDINAL[k], label="{} micrometres".format(d))
ax[0].set_xlim(1e-6, 1e11); ax[0].set_ylim(1e-7, 1.2e4)
ax[0].set_xlabel("time since release [s]"); ax[0].set_ylabel("speed relative to the air [m/s]")
ax[0].set_title("Speed after release at 80 km", loc="left")
ax[0].legend(loc="lower left", fontsize=8.5)
for d_, lab in ((1, "1 day"), (365, "1 year"), (365e3, "1000 years")):
    ax[0].axvline(d_ * 86400, color=MUTED, lw=0.8, ls=":")
    ax[0].text(d_ * 86400 * 1.15, 5e3, lab, color=INK2, fontsize=8)
ax[1].set_xlim(1e-3, 300)
ax[1].set_xlabel("downrange distance [km]"); ax[1].set_ylabel("altitude [km]")
ax[1].set_title("Path after release at 80 km (down to 60 km)", loc="left")
ax[1].legend(loc="lower left", fontsize=8.5)
save(fig, "5_trajectories.png")

# ---------------------------------------------------------------- 6. regime map
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
for a, h in zip(ax, (90, 60)):
    a.add_patch(Rectangle((0.1, 0.05), 20 - 0.1, 100, color=CRITICAL, alpha=0.08, lw=0))
    a.text(18, 25, "uncalibrated corner of the drag law:\nKnudsen 0.1-20 with speed ratio above 0.05",
           color=INK2, fontsize=8.5)
    a.axhline(3, color=MUTED, lw=0.8, ls=":")
    a.text(1.5e-3, 3.4, "speed ratio 3", color=INK2, fontsize=8, ha="right")
    for k, d in enumerate(sel):
        r = hist(d, h)
        a.loglog(r.Kn, r.speed_ratio, color=ORDINAL[k], label="{} micrometres".format(d))
        a.plot(r.Kn.iloc[0], r.speed_ratio.iloc[0], "o", ms=8, color=ORDINAL[k], mec=SURFACE, mew=1.5)
    a.set_xlim(1e-3, 3e4); a.set_ylim(1e-6, 60)
    a.set_xlabel("Knudsen number (air mean free path / droplet size, x2)")
    a.set_title("Flow regimes crossed after release at {} km (dot = release)".format(h), loc="left")
    a.invert_xaxis()
ax[0].set_ylabel("molecular speed ratio (speed / thermal speed of air)")
ax[0].legend(loc="lower right", fontsize=8.5)
save(fig, "6_regime_map.png")
