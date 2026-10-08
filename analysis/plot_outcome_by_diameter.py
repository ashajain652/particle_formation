#!/usr/bin/env python3
"""
plot_outcome_by_diameter.py -- small multiples: one heatmap per sphere diameter,
x = initial temperature, y = initial velocity, color = the chosen metric.

Every (temperature, velocity) cell has exactly one run per diameter, so:

  --metric outcome      (default) binary demised(1)/survived(0) outcome of
                         that single run -- the panel is the survive/demise
                         boundary in (temperature, velocity) space
  --metric mass_loss     that run's mass_loss_fraction (continuous 0..1)

Usage
-----
    python analysis/plot_outcome_by_diameter.py
    python analysis/plot_outcome_by_diameter.py --metric mass_loss
    python analysis/plot_outcome_by_diameter.py --summary path/to/sweep_summary.csv --out path/to/plot.png

Requires: pandas, numpy, matplotlib.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (SEQ_BLUE, INK, SECOND, MUTED, DEFAULT_SUMMARY, DEFAULT_PLOT_DIR,
                        apply_rcparams, strip_top_right_spines, style_colorbar, material_from_run_name)

METRIC_LABEL = {"outcome": "demised (1) vs. survived (0)", "mass_loss": "mass-loss fraction"}
METRIC_NOTE = {
    "outcome": "dark = demised, light = survived (one run per cell, binary)",
    "mass_loss": "mass_loss_fraction = 1 - final_mass/initial_mass (one run per cell, continuous)",
}
METRIC_DEFAULT_OUT = {"outcome": "demise_fraction_by_diameter.png",
                      "mass_loss": "mass_loss_fraction_by_diameter.png"}
N_COLS, N_ROWS = 5, 4


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default=DEFAULT_SUMMARY, help="sweep_summary.csv path (default %(default)s)")
    p.add_argument("--metric", choices=["outcome", "mass_loss"], default="outcome",
                   help="'outcome' = binary demised/survived per run; 'mass_loss' = mass_loss_fraction per run")
    p.add_argument("--out", default=None, help="output PNG path (default: sphere_sweep_output/plots/<metric-name>.png)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = args.out or os.path.join(DEFAULT_PLOT_DIR, METRIC_DEFAULT_OUT[args.metric])
    os.makedirs(os.path.dirname(out), exist_ok=True)

    df = pd.read_csv(args.summary)
    df = df[df["status"] == "done"].copy()
    if args.metric == "outcome":
        df["value"] = (df["outcome"] == "demised").astype(float)
    else:
        df["value"] = df["mass_loss_fraction"]
    print("rows:", len(df))

    diameters = sorted(df["diameter_mm"].unique())
    temperatures = sorted(df["initial_temperature_K"].unique())
    velocities = sorted(df["initial_velocity_kms"].unique())
    if len(diameters) != N_COLS * N_ROWS:
        sys.exit("expected {} diameters for a {}x{} grid of panels, found {}".format(
            N_COLS * N_ROWS, N_ROWS, N_COLS, len(diameters)))
    t_index = {t: i for i, t in enumerate(temperatures)}
    v_index = {v: i for i, v in enumerate(velocities)}
    t_arr, v_arr = np.array(temperatures), np.array(velocities)

    grids = {}
    for d, g in df.groupby("diameter_mm"):
        grid = np.full((len(velocities), len(temperatures)), np.nan)
        for r in g.itertuples():
            grid[v_index[r.initial_velocity_kms], t_index[r.initial_temperature_K]] = r.value
        grids[d] = grid

    dt, dv = t_arr[1] - t_arr[0], v_arr[1] - v_arr[0]
    extent = [t_arr[0] - dt / 2, t_arr[-1] + dt / 2, v_arr[0] - dv / 2, v_arr[-1] + dv / 2]

    apply_rcparams(plt, font_size=9.5)
    fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(15.5, 11.5), sharex=True, sharey=True)

    im = None
    for i, d in enumerate(diameters):
        row, col = divmod(i, N_COLS)
        ax = axes[row, col]
        im = ax.imshow(grids[d], origin="lower", aspect="auto", extent=extent,
                       cmap=SEQ_BLUE, vmin=0, vmax=1, interpolation="nearest")
        ax.set_title("d = {:.0f} mm".format(d), fontsize=10, color=INK, loc="left",
                    fontweight="bold", pad=4)
        strip_top_right_spines(ax)
        ax.set_xticks([300, 450, 600, 750])
        ax.set_yticks([0, 2.5, 5, 7.5])
        if row == N_ROWS - 1:
            ax.set_xlabel("initial temperature  [K]", fontsize=9, color=SECOND)
        if col == 0:
            ax.set_ylabel("initial velocity  [km/s]", fontsize=9, color=SECOND)
        ax.tick_params(labelsize=8)

    fig.suptitle("{} sphere fragments: {} by initial temperature and velocity, per sphere size".format(
                 material_from_run_name(df["run_name"].iloc[0]),
                 "demise outcome" if args.metric == "outcome" else "mass-loss fraction"),
                fontsize=15, color=INK, x=0.055, ha="left", y=0.99, fontweight="bold")
    fig.text(0.055, 0.965,
             "{}; {} temperatures ({}-{} K) x {} velocities ({}->{} km/s), {} sphere sizes ({}-{} mm)".format(
                 METRIC_NOTE[args.metric], len(temperatures), int(temperatures[0]), int(temperatures[-1]),
                 len(velocities), velocities[-1], velocities[0], len(diameters), int(diameters[0]), int(diameters[-1])),
             fontsize=10, color=MUTED, ha="left")

    fig.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.06, hspace=0.38, wspace=0.18)
    cax = fig.add_axes([0.925, 0.15, 0.015, 0.65])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(METRIC_LABEL[args.metric], color=SECOND, fontsize=10)
    style_colorbar(cbar)

    fig.savefig(out, dpi=160, facecolor="#fcfcfb")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
