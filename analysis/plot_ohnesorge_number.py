#!/usr/bin/env python3
"""
plot_ohnesorge_number.py -- Ohnesorge number of every sphere at its initial state,
x = initial diameter, y = Oh.

    Oh = mu_liq / sqrt(rho_liq * sigma * d)      (= sqrt(We_liq) / Re_liq)

the standard liquid-drop definition used in breakup regime maps, with the molten
sphere's viscosity mu_liq, density rho_liq and surface tension sigma and the initial
diameter d [m]. It contains no velocity and no air property, so it is one curve in d.
The table written next to the plot still pairs every (diameter, velocity) point's
Weber number (plot_weber_number.py) with its Ohnesorge number.

Usage
-----
    python analysis/plot_ohnesorge_number.py --summary sphere_sweep_output/T850_moltenAA7075/sweep_summary_user-moltenAA7075.csv
    python analysis/plot_ohnesorge_number.py --summary ... --mu 1.2e-3 --rho-liquid 2400 --sigma 0.809

Requires: pandas, numpy, matplotlib.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (BLUE_RAMP, INK, SECOND, MUTED, GRID, DEFAULT_SUMMARY, DEFAULT_PLOT_DIR,
                        apply_rcparams, strip_top_right_spines)
from plot_weber_number import DEFAULT_SIGMA, weber_table

DEFAULT_MU_LIQUID = 1.2e-3        # Pa s, liquid aluminium near its melting point
DEFAULT_RHO_LIQUID = 2400.0       # kg/m3, density of the moltenAA7075 material
VISCOUS_OH = 0.1                  # above this, viscosity starts to raise the critical We


def ohnesorge_number(diameter_mm, mu, rho_liquid, sigma):
    return mu / np.sqrt(rho_liquid * sigma * diameter_mm / 1000.0)


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default=DEFAULT_SUMMARY, help="sweep summary CSV (default %(default)s)")
    p.add_argument("--runs-dir", default=None, help="run histories (default: runs/ next to the summary)")
    p.add_argument("--mu", type=float, default=DEFAULT_MU_LIQUID, help="liquid viscosity [Pa s] (default %(default)s)")
    p.add_argument("--rho-liquid", type=float, default=DEFAULT_RHO_LIQUID, help="liquid density [kg/m3] (default %(default)s)")
    p.add_argument("--sigma", type=float, default=DEFAULT_SIGMA, help="surface tension [N/m] (default %(default)s)")
    p.add_argument("--out", default=None, help="output PNG (default sphere_sweep_output/plots/ohnesorge_number.png)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = args.out or os.path.join(DEFAULT_PLOT_DIR, "ohnesorge_number.png")
    runs_dir = args.runs_dir or os.path.join(os.path.dirname(os.path.abspath(args.summary)), "runs")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    import pandas as pd
    df = pd.read_csv(args.summary)
    table = weber_table(df, runs_dir, args.sigma)
    table["ohnesorge"] = ohnesorge_number(table["diameter_mm"], args.mu, args.rho_liquid, args.sigma)
    table_path = os.path.join(os.path.dirname(out), "ohnesorge_numbers.csv")
    table.to_csv(table_path, index=False, float_format="%.6g")
    oh_min, oh_max = float(table["ohnesorge"].min()), float(table["ohnesorge"].max())
    print("{} points; Oh {:.3g} .. {:.3g} (mu = {} Pa s, rho_liq = {} kg/m3, sigma = {} N/m)".format(
        len(table), oh_min, oh_max, args.mu, args.rho_liquid, args.sigma))

    curve = table.groupby("diameter_mm")["ohnesorge"].first().sort_index()
    d_arr, oh_arr = curve.index.to_numpy(), curve.to_numpy()

    apply_rcparams(plt)
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    ax.plot(d_arr, oh_arr, color=BLUE_RAMP[4], lw=2, marker="o", ms=5, markerfacecolor=BLUE_RAMP[4],
            markeredgewidth=0, zorder=3)
    for d, oh in zip(d_arr, oh_arr):
        ax.annotate("{:.2e}".format(oh), (d, oh), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=7, color=SECOND, rotation=45)
    ax.set_xlabel("initial sphere diameter  [mm]")
    ax.set_ylabel("Ohnesorge number  Oh = mu / sqrt(rho_liq sigma d)")
    ax.set_xticks(d_arr[::2])
    ax.set_xlim(0, d_arr[-1] + 5)
    ax.set_ylim(0, oh_max * 1.18)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    strip_top_right_spines(ax)

    fig.suptitle("Ohnesorge number of each sphere at its initial state", fontsize=13.5, color=INK,
                 x=0.06, ha="left", y=1.04, fontweight="bold")
    ax.set_title("mu = {:g} Pa s, rho_liq = {:g} kg/m3, sigma = {:g} N/m; Oh depends on the diameter only "
                 "(no velocity or air property enters);\nall values are far below Oh = {:g}, where viscosity "
                 "would start to raise the critical Weber number".format(
                     args.mu, args.rho_liquid, args.sigma, VISCOUS_OH),
                 fontsize=9, color=MUTED, loc="left", pad=8)

    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="#fcfcfb")
    print("wrote", out, "and", table_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
