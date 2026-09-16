#!/usr/bin/env python3
"""
plot_outcome_vs_diameter_velocity.py -- demise fraction (or mass-loss fraction)
vs. sphere diameter and initial velocity, averaged over initial temperature.

Reads sweep_summary.csv (written by sphere_sweep.py) and, for each (diameter,
velocity) grid cell, averages the chosen metric over the 46 initial-temperature
runs at that cell:

  --metric outcome      (default) fraction of runs with outcome == "demised"
                         (see spec: ground impact -> survived even with heavy
                         ablation; demised means (near-)total mass loss)
  --metric mass_loss     mean mass_loss_fraction = 1 - final_mass/initial_mass
                         (continuous 0..1, independent of the ground-impact rule)

Produces a heatmap (diameter x velocity) plus the two 1-D marginals.

Usage
-----
    python analysis/plot_outcome_vs_diameter_velocity.py
    python analysis/plot_outcome_vs_diameter_velocity.py --metric mass_loss
    python analysis/plot_outcome_vs_diameter_velocity.py --summary path/to/sweep_summary.csv --out path/to/plot.png

Requires: pandas, numpy, matplotlib (not required by sphere_reentry.py /
sphere_sweep.py themselves -- this is a standalone post-processing script).
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
from plot_style import (BLUE_RAMP, SEQ_BLUE, INK, SECOND, MUTED, GRID,
                        DEFAULT_SUMMARY, DEFAULT_PLOT_DIR, apply_rcparams,
                        strip_top_right_spines, style_colorbar, material_from_run_name)

METRIC_COLUMN = {"outcome": None, "mass_loss": "mass_loss_fraction"}
METRIC_LABEL = {"outcome": "demise fraction", "mass_loss": "mean mass-loss fraction"}
METRIC_TITLE = {"outcome": "Demise fraction", "mass_loss": "Mean mass-loss fraction"}
METRIC_DEFAULT_OUT = {"outcome": "demise_fraction_vs_diameter_velocity.png",
                      "mass_loss": "mass_loss_fraction_vs_diameter_velocity.png"}


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default=DEFAULT_SUMMARY, help="sweep_summary.csv path (default %(default)s)")
    p.add_argument("--metric", choices=["outcome", "mass_loss"], default="outcome",
                   help="'outcome' = fraction demised (binary per run); "
                        "'mass_loss' = mean mass_loss_fraction (continuous per run)")
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
        df["value"] = df[METRIC_COLUMN[args.metric]]
    print("rows:", len(df))

    diameters = sorted(df["diameter_mm"].unique())
    velocities = sorted(df["initial_velocity_kms"].unique())
    d_index = {d: i for i, d in enumerate(diameters)}
    v_index = {v: i for i, v in enumerate(velocities)}

    total = np.zeros((len(diameters), len(velocities)))
    count = np.zeros((len(diameters), len(velocities)))
    for r in df.itertuples():
        d, v = d_index[r.diameter_mm], v_index[r.initial_velocity_kms]
        total[d, v] += r.value
        count[d, v] += 1
    mean_value = np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)
    print("cell counts min/max:", count.min(), count.max())
    print("overall mean:", df["value"].mean())

    frac_by_diameter = np.nanmean(mean_value, axis=1)
    frac_by_velocity = np.nanmean(mean_value, axis=0)

    apply_rcparams(plt)
    fig = plt.figure(figsize=(13.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.42)
    ax_heat = fig.add_subplot(gs[0, 0])
    ax_dia = fig.add_subplot(gs[0, 1])
    ax_vel = fig.add_subplot(gs[0, 2])

    v_arr, d_arr = np.array(velocities), np.array(diameters)
    order = np.argsort(v_arr)
    v_sorted = v_arr[order]
    value_sorted = mean_value[:, order]

    dv, dd = v_sorted[1] - v_sorted[0], d_arr[1] - d_arr[0]
    extent = [v_sorted[0] - dv / 2, v_sorted[-1] + dv / 2, d_arr[0] - dd / 2, d_arr[-1] + dd / 2]

    im = ax_heat.imshow(value_sorted, origin="lower", aspect="auto", extent=extent,
                        cmap=SEQ_BLUE, vmin=0, vmax=1, interpolation="nearest")
    ax_heat.set_xlabel("initial velocity  [km/s]")
    ax_heat.set_ylabel("sphere diameter  [mm]")
    ax_heat.set_title("{} by diameter × velocity".format(METRIC_TITLE[args.metric]), fontsize=11.5,
                      color=INK, loc="left", fontweight="bold", pad=10)
    ax_heat.set_yticks(diameters[::2])
    strip_top_right_spines(ax_heat)
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)
    cbar.set_label("{}  (across 46 initial temperatures)".format(METRIC_LABEL[args.metric]), color=SECOND, fontsize=9.5)
    style_colorbar(cbar)

    ax_dia.plot(d_arr, frac_by_diameter, color=BLUE_RAMP[4], lw=2, marker="o", ms=4,
               markerfacecolor=BLUE_RAMP[4], markeredgewidth=0)
    ax_dia.set_xlabel("sphere diameter  [mm]")
    ax_dia.set_ylabel(METRIC_LABEL[args.metric])
    ax_dia.set_ylim(-0.03, 1.03)
    ax_dia.set_title("vs. diameter", fontsize=11, color=INK, loc="left", fontweight="bold", pad=10)
    ax_dia.grid(True, color=GRID, lw=0.8, zorder=0)
    ax_dia.set_axisbelow(True)
    strip_top_right_spines(ax_dia)

    ax_vel.plot(v_sorted, frac_by_velocity, color=BLUE_RAMP[4], lw=2)
    ax_vel.set_xlabel("initial velocity  [km/s]")
    ax_vel.set_ylabel(METRIC_LABEL[args.metric])
    ax_vel.set_ylim(-0.03, 1.03)
    ax_vel.invert_xaxis()
    ax_vel.set_title("vs. velocity", fontsize=11, color=INK, loc="left", fontweight="bold", pad=10)
    ax_vel.grid(True, color=GRID, lw=0.8, zorder=0)
    ax_vel.set_axisbelow(True)
    strip_top_right_spines(ax_vel)

    fig.suptitle("{} sphere fragments: {} across the sweep".format(
                 material_from_run_name(df["run_name"].iloc[0]), METRIC_LABEL[args.metric]),
                fontsize=13.5, color=INK, x=0.06, ha="left", y=1.03, fontweight="bold")
    metric_note = ("outcome == \"demised\" per run (ground impact always counts as survived)"
                  if args.metric == "outcome" else
                  "mass_loss_fraction = 1 - final_mass/initial_mass per run")
    fig.text(0.06, 0.965,
             "{} diameters x {} temperatures x {} velocities -- {}, averaged over "
             "the 46 temperatures per cell".format(len(diameters), 46, len(velocities), metric_note),
             fontsize=9.5, color=MUTED, ha="left")

    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="#fcfcfb")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
