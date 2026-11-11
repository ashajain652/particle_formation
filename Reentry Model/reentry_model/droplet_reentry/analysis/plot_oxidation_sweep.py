"""Figures of the oxidation bounding sweep (oxidation_sweep.py), into
output/droplet_reentry_output/oxidation_sweep/figures/:

  1_oxidation_bounds.{png,pdf,svg}  rows: share of the initial mass vaporised, peak temperature, oxide share of the
                                    final particle; columns: share of the oxide left as a surface skin; one line per
                                    oxygen sticking efficiency, against release altitude (10 um droplet)
  2_size_check.{png,pdf,svg}        the vaporised share and oxide share for 1, 10 and 100 um at full sticking
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from reentry_model.droplet_reentry.analysis import oxidation_sweep as ox

FIG = os.path.join(ox.OUT, "figures")
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
STICK_COLOURS = {0.0: MUTED, 0.01: "#86b6ef", 0.1: "#2a78d6", 1.0: "#0d366b"}
SIZE_STYLES = {1.0: (0, (1, 1.5)), 10.0: "-", 100.0: (0, (5, 2))}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 12, "axes.axisbelow": True})
METRICS = [("vaporised_frac", 100.0, "Initial mass\nvaporised (%)"),
           ("T_peak_K", 1.0, "Peak droplet\ntemperature (K)"),
           ("oxide_mass_frac", 100.0, "Oxide share of\nfinal particle (%)")]


def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(FIG, name + "." + ext), dpi=200)


def main():
    os.makedirs(FIG, exist_ok=True)
    f = pd.read_csv(os.path.join(ox.OUT, "final_states.csv"))
    main_d = f[f.d0_um == ox.D_MAIN]
    base = main_d[main_d.sticking == 0.0].sort_values("h0_km")

    fig, axes = plt.subplots(len(METRICS), len(ox.SKIN), figsize=(15, 9.5), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.08, hspace=0.12, wspace=0.08)
    fig.text(0.08, 0.95, "Oxidation bounds: how much the droplet vaporises depends on where the oxide goes",
             fontsize=17, weight="bold")
    fig.text(0.08, 0.915, "{:g} µm AA7075 droplet; Genge-style oxygen capture at sticking efficiency γ, Al → Al₂O₃ "
             "heat into the droplet; a share of the oxide stays on the surface as a skin".format(ox.D_MAIN),
             fontsize=12, color=INK2)
    for j, fs in enumerate(ox.SKIN):
        for i, (col, scale, label) in enumerate(METRICS):
            ax = axes[i, j]
            ax.plot(base.h0_km, base[col] * scale, color=STICK_COLOURS[0.0], lw=1.6, ls=(0, (4, 2)),
                    label="no oxidation")
            for g in ox.STICKING[1:]:
                s = main_d[(main_d.sticking == g) & (main_d.f_skin == fs)].sort_values("h0_km")
                ax.plot(s.h0_km, s[col] * scale, color=STICK_COLOURS[g], label="γ = {:g}".format(g))
            if j == 0:
                ax.set_ylabel(label)
            if i == len(METRICS) - 1:
                ax.set_xlabel("Release altitude (km)")
        title = {0.0: "all entrained (fresh surface)", 1.0: "all on the surface"}.get(fs, "{:g} % on the surface"
                                                                                       .format(100 * fs))
        axes[0, j].set_title("Skin share {:g}: {}".format(fs, title), loc="left", fontsize=12)
    axes[0, 0].legend(frameon=False, fontsize=10.5, loc="upper left")
    save(fig, "1_oxidation_bounds")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.13, wspace=0.22)
    fig.text(0.08, 0.93, "Size check at full sticking (γ = 1)", fontsize=16, weight="bold")
    colours = {0.0: "#2a78d6", 1.0: "#c2410c"}
    for k, (col, scale, label) in enumerate([METRICS[0], METRICS[2]]):
        ax = axes[k]
        for fs in ox.SIZE_CHECK["skin"]:
            for d in (1.0, 10.0, 100.0):
                s = f[(f.d0_um == d) & (f.sticking == 1.0) & (f.f_skin == fs)
                      & f.h0_km.isin(ox.SIZE_CHECK["h"])].sort_values("h0_km")
                ax.plot(s.h0_km, s[col] * scale, color=colours[fs], ls=SIZE_STYLES[d], marker="o", ms=4,
                        label="{:g} µm, skin share {:g}".format(d, fs))
        ax.set_ylabel(label.replace("\n", " "))
        ax.set_xlabel("Release altitude (km)")
    axes[0].legend(frameon=False, fontsize=9.5, ncol=2, loc="upper left", bbox_to_anchor=(0.0, 1.22))
    save(fig, "2_size_check")
    print(FIG)


if __name__ == "__main__":
    main()
