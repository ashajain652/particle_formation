#!/usr/bin/env python3
"""
plot_weber_number.py -- Weber number of every sphere at its initial state,
diameter (y) x initial velocity (x), color = We on a log scale.

    We = rho_air * v^2 * d / sigma

with rho_air the air density at the sphere's initial altitude, v the initial velocity
[m/s], d the initial diameter [m] and sigma the surface tension of the (molten) sphere
[N/m]. The initial altitude is the parent trajectory's altitude at that velocity (the
sweep's inheritance rule), so rho depends on velocity only, not on initial temperature
or material. rho is DRAMA's own atmosphere model: the `density_kgm3` value at t = 0 of
each run's history CSV under runs/. Only one run per (diameter, velocity) is read --
the first temperature of the summary -- because all temperatures share the state.

Also writes the table of Weber numbers next to the plot (weber_numbers.csv:
diameter_mm, initial_velocity_kms, initial_altitude_km, density_kgm3, weber).

Usage
-----
    python analysis/plot_weber_number.py --summary sphere_sweep_output/T850_moltenAA7075/sweep_summary_user-moltenAA7075.csv
    python analysis/plot_weber_number.py --summary ... --sigma 0.809 --out path/to/plot.png

Requires: pandas, numpy, matplotlib.
"""
import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (SEQ_BLUE, INK, SECOND, MUTED, DEFAULT_SUMMARY, DEFAULT_PLOT_DIR,
                        apply_rcparams, strip_top_right_spines, style_colorbar)

DEFAULT_SIGMA = 0.809          # N/m, surface tension of the molten sphere
CRITICAL_WE = 12.0             # classic bag-breakup threshold, drawn as a reference contour


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default=DEFAULT_SUMMARY, help="sweep summary CSV (default %(default)s)")
    p.add_argument("--runs-dir", default=None, help="run histories (default: runs/ next to the summary)")
    p.add_argument("--sigma", type=float, default=DEFAULT_SIGMA, help="surface tension [N/m] (default %(default)s)")
    p.add_argument("--out", default=None, help="output PNG (default sphere_sweep_output/plots/weber_number.png)")
    return p


def initial_density(runs_dir, run_name):
    """Air density [kg/m3] at t = 0 of one run, from DRAMA's history."""
    h = pd.read_csv(os.path.join(runs_dir, run_name + ".csv"), nrows=1, usecols=["time_s", "density_kgm3"])
    if h.empty or h["time_s"].iloc[0] != 0.0:
        raise ValueError("{}: first history row is not t = 0".format(run_name))
    return float(h["density_kgm3"].iloc[0])


def weber_number(density_kgm3, velocity_kms, diameter_mm, sigma):
    v = velocity_kms * 1000.0          # km/s -> m/s
    d = diameter_mm / 1000.0           # mm -> m
    return density_kgm3 * v ** 2 * d / sigma


def weber_table(df, runs_dir, sigma):
    """One row per (diameter, velocity) with altitude, density and We."""
    first_T = df["initial_temperature_K"].min()
    one = df[(df["status"] == "done") & (df["initial_temperature_K"] == first_T)]
    rows = []
    for r in one.itertuples():
        rho = initial_density(runs_dir, r.run_name)
        rows.append({"diameter_mm": r.diameter_mm, "initial_velocity_kms": r.initial_velocity_kms,
                     "initial_altitude_km": r.initial_altitude_km, "density_kgm3": rho,
                     "weber": weber_number(rho, r.initial_velocity_kms, r.diameter_mm, sigma)})
    return pd.DataFrame(rows).sort_values(["diameter_mm", "initial_velocity_kms"], ascending=[True, False])


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = args.out or os.path.join(DEFAULT_PLOT_DIR, "weber_number.png")
    runs_dir = args.runs_dir or os.path.join(os.path.dirname(os.path.abspath(args.summary)), "runs")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    df = pd.read_csv(args.summary)
    table = weber_table(df, runs_dir, args.sigma)
    table_path = os.path.join(os.path.dirname(out), "weber_numbers.csv")
    table.to_csv(table_path, index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)
    print("{} (diameter, velocity) points; We {:.3g} .. {:.3g}; sigma = {} N/m".format(
        len(table), table["weber"].min(), table["weber"].max(), args.sigma))

    diameters = sorted(table["diameter_mm"].unique())
    velocities = sorted(table["initial_velocity_kms"].unique())
    d_index = {d: i for i, d in enumerate(diameters)}
    v_index = {v: i for i, v in enumerate(velocities)}
    grid = np.full((len(diameters), len(velocities)), np.nan)
    for r in table.itertuples():
        grid[d_index[r.diameter_mm], v_index[r.initial_velocity_kms]] = r.weber
    d_arr, v_arr = np.array(diameters), np.array(velocities)
    dv, dd = v_arr[1] - v_arr[0], d_arr[1] - d_arr[0]
    extent = [v_arr[0] - dv / 2, v_arr[-1] + dv / 2, d_arr[0] - dd / 2, d_arr[-1] + dd / 2]

    apply_rcparams(plt)
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    im = ax.imshow(grid, origin="lower", aspect="auto", extent=extent, cmap=SEQ_BLUE,
                   norm=LogNorm(vmin=table["weber"].min(), vmax=table["weber"].max()), interpolation="nearest")
    cs = ax.contour(v_arr, d_arr, grid, levels=[CRITICAL_WE], colors=["#e34948"], linewidths=1.4)
    ax.clabel(cs, fmt={CRITICAL_WE: "We = {:g}".format(CRITICAL_WE)}, fontsize=9, colors=["#e34948"])
    ax.set_xlabel("initial velocity  [km/s]")
    ax.set_ylabel("initial sphere diameter  [mm]")
    ax.set_yticks(diameters[::2])
    strip_top_right_spines(ax)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Weber number  We = rho v^2 d / sigma  (log scale)", color=SECOND, fontsize=9.5)
    style_colorbar(cbar)

    fig.suptitle("Weber number of each sphere at its initial state", fontsize=13.5, color=INK,
                 x=0.06, ha="left", y=0.995, fontweight="bold")
    ax.set_title("sigma = {} N/m; rho = DRAMA atmosphere at the initial altitude, inherited from the parent "
                 "trajectory at that velocity\n({:.1f} km at {:.1f} km/s down to {:.2f} km at {:.3f} km/s); "
                 "red line: We = {:g}".format(
                     args.sigma, table.loc[table["initial_velocity_kms"].idxmax(), "initial_altitude_km"],
                     v_arr[-1], table.loc[table["initial_velocity_kms"].idxmin(), "initial_altitude_km"], v_arr[0],
                     CRITICAL_WE),
                 fontsize=9, color=MUTED, loc="left", pad=8)

    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="#fcfcfb")
    print("wrote", out, "and", table_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
