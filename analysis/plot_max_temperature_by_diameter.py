#!/usr/bin/env python3
"""
plot_max_temperature_by_diameter.py -- small multiples: one panel per sphere diameter,
x = initial temperature, y = initial velocity, color = the run's maximum temperature,
shown ONLY where the sphere lost no mass (mass_loss_fraction == 0). Every cell whose
sphere did lose mass is drawn black.

Reads a sweep summary CSV. Because sphere_sweep.py rewrites sweep_summary.csv on every
invocation, point --summary at a per-material file built by analysis/summary_from_runs.py
when several materials share the output directory.

Usage
-----
    python analysis/plot_max_temperature_by_diameter.py --summary sphere_sweep_output/sweep_summary_AA7075.csv
    python analysis/plot_max_temperature_by_diameter.py --summary ... --out path/to/plot.png

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
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (SEQ_BLUE, INK, SECOND, MUTED, DEFAULT_SUMMARY, DEFAULT_PLOT_DIR,
                        apply_rcparams, strip_top_right_spines, style_colorbar)

ZERO_LOSS_TOLERANCE = 1e-9      # |mass_loss_fraction| below this counts as "no mass lost"
MASS_LOSS_COLOR = "#0b0b0b"     # cells whose sphere lost mass
N_COLS, N_ROWS = 5, 4


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default=DEFAULT_SUMMARY, help="summary CSV path (default %(default)s)")
    p.add_argument("--out", default=None,
                   help="output PNG path (default sphere_sweep_output/plots/max_temperature_by_diameter.png)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = args.out or os.path.join(DEFAULT_PLOT_DIR, "max_temperature_by_diameter.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    df = pd.read_csv(args.summary)
    df = df[df["status"] == "done"].copy()
    no_loss = df["mass_loss_fraction"].abs() <= ZERO_LOSS_TOLERANCE
    df["value"] = np.where(no_loss, df["max_temperature_K"], np.nan)   # NaN -> drawn black
    print("rows:", len(df), "| no mass lost:", int(no_loss.sum()), "| lost mass:", int((~no_loss).sum()))

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

    vmin = float(np.nanmin(df["value"]))
    vmax = float(np.nanmax(df["value"]))
    print("max temperature among no-loss runs: {:.1f} .. {:.1f} K".format(vmin, vmax))
    cmap = SEQ_BLUE.copy()
    cmap.set_bad(MASS_LOSS_COLOR)

    dt, dv = t_arr[1] - t_arr[0], v_arr[1] - v_arr[0]
    extent = [t_arr[0] - dt / 2, t_arr[-1] + dt / 2, v_arr[0] - dv / 2, v_arr[-1] + dv / 2]

    apply_rcparams(plt, font_size=9.5)
    fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(15.5, 11.5), sharex=True, sharey=True)

    im = None
    for i, d in enumerate(diameters):
        row, col = divmod(i, N_COLS)
        ax = axes[row, col]
        im = ax.imshow(np.ma.masked_invalid(grids[d]), origin="lower", aspect="auto", extent=extent,
                       cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
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

    material = df["run_name"].iloc[0].split("_m")[1] if "_m" in df["run_name"].iloc[0] else "drama-AA7075"
    fig.suptitle("{} sphere fragments: maximum temperature reached by spheres that lost no mass, per sphere size".format(material),
                fontsize=14.5, color=INK, x=0.055, ha="left", y=0.99, fontweight="bold")
    fig.text(0.055, 0.965,
             "color = max temperature of the run where mass_loss_fraction == 0; black = the sphere lost mass; "
             "{} temperatures ({}-{} K) x {} velocities ({}->{} km/s), {} sphere sizes ({}-{} mm)".format(
                 len(temperatures), int(temperatures[0]), int(temperatures[-1]),
                 len(velocities), velocities[-1], velocities[0], len(diameters), int(diameters[0]), int(diameters[-1])),
             fontsize=10, color=MUTED, ha="left")

    fig.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.06, hspace=0.38, wspace=0.18)
    cax = fig.add_axes([0.925, 0.22, 0.015, 0.58])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("maximum temperature  [K]  (no mass lost)", color=SECOND, fontsize=10)
    style_colorbar(cbar)
    fig.legend(handles=[Patch(facecolor=MASS_LOSS_COLOR, label="mass loss > 0")],
               loc="lower left", bbox_to_anchor=(0.918, 0.12), frameon=False, fontsize=9.5,
               labelcolor=SECOND)

    fig.savefig(out, dpi=160, facecolor="#fcfcfb", bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
