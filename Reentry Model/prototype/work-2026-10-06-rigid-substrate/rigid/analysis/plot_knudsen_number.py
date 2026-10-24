#!/usr/bin/env python3
"""
plot_knudsen_number.py -- Knudsen number of every sphere at its initial state,
diameter (y) x initial velocity (x), color = Kn on a log scale banded by flow regime:
blue Kn < 0.01 (continuum), yellow 0.01-0.1 (slip), orange 0.1-10 (transitional),
red >= 10 (free molecular), each band shaded light -> dark with increasing Kn.

    Kn = lambda / d,   lambda = C / rho_air

lambda is the mean free path of the air at the sphere's initial altitude (inherited
from the parent trajectory at that velocity, so it depends on velocity only) and d the
initial diameter [m]. SESAM writes its own Knudsen number (reference length = diameter)
to every history file, but only to 5 decimals, which truncates to 0 low in the
atmosphere for larger spheres. Its values follow lambda = C / rho exactly, so C is
calibrated from the sphere's own runs (every point with SESAM Kn >= 0.01, i.e. at least
three significant digits; the fit reproduces SESAM within 0.2 %) and Kn is then
evaluated at full precision from the t = 0 air density of each run. The calibrated C
and the fit residual are printed and written into the table.

Usage
-----
    python analysis/plot_knudsen_number.py --summary sphere_sweep_output/T850_moltenAA7075/sweep_summary_user-moltenAA7075.csv
    python analysis/plot_knudsen_number.py --summary ... --out path/to/plot.png

Writes knudsen_numbers.csv next to the plot (diameter_mm, initial_velocity_kms,
initial_altitude_km, density_kgm3, mean_free_path_m, knudsen, knudsen_sesam).
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
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (INK, SECOND, MUTED, DEFAULT_SUMMARY, DEFAULT_PLOT_DIR,
                        apply_rcparams, strip_top_right_spines, style_colorbar)
from plot_weber_number import banded_cmap

# flow regimes: (lower Kn bound, name, light shade, dark shade); upper bound = next lower bound
KN_BANDS = [
    (0.0,  "Kn < 0.01: continuum",       "#cde2fb", "#0d366b"),   # blue
    (0.01, "0.01-0.1: slip",             "#fff0b3", "#c99700"),   # yellow
    (0.1,  "0.1-10: transitional",       "#ffd9b0", "#c85400"),   # orange
    (10.0, "Kn >= 10: free molecular",   "#f6b3b0", "#8b1a1a"),   # red
]
MIN_RESOLVED_KN = 0.01         # SESAM prints 5 decimals: >= 3 significant digits above this


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default=DEFAULT_SUMMARY, help="sweep summary CSV (default %(default)s)")
    p.add_argument("--runs-dir", default=None, help="run histories (default: runs/ next to the summary)")
    p.add_argument("--out", default=None, help="output PNG (default sphere_sweep_output/plots/knudsen_number.png)")
    return p


def initial_state(runs_dir, run_name):
    """(air density [kg/m3], SESAM Knudsen number) at t = 0 of one run's history."""
    h = pd.read_csv(os.path.join(runs_dir, run_name + ".csv"), nrows=1,
                    usecols=["time_s", "density_kgm3", "knudsen"])
    if h.empty or h["time_s"].iloc[0] != 0.0:
        raise ValueError("{}: first history row is not t = 0".format(run_name))
    return float(h["density_kgm3"].iloc[0]), float(h["knudsen"].iloc[0])


def knudsen_table(df, runs_dir):
    """One row per (diameter, velocity): altitude, density, calibrated mean free path, Kn."""
    first_T = df["initial_temperature_K"].min()
    one = df[(df["status"] == "done") & (df["initial_temperature_K"] == first_T)]
    rows = []
    for r in one.itertuples():
        rho, kn_sesam = initial_state(runs_dir, r.run_name)
        rows.append({"diameter_mm": r.diameter_mm, "initial_velocity_kms": r.initial_velocity_kms,
                     "initial_altitude_km": r.initial_altitude_km, "density_kgm3": rho, "knudsen_sesam": kn_sesam})
    t = pd.DataFrame(rows)
    resolved = t[t["knudsen_sesam"] >= MIN_RESOLVED_KN]
    if resolved.empty:
        raise ValueError("no run has a SESAM Knudsen number >= {} to calibrate the mean free path".format(MIN_RESOLVED_KN))
    C = float((resolved["knudsen_sesam"] * resolved["diameter_mm"] / 1000.0 * resolved["density_kgm3"]).median())
    t["mean_free_path_m"] = C / t["density_kgm3"]
    t["knudsen"] = t["mean_free_path_m"] / (t["diameter_mm"] / 1000.0)
    fit = t.loc[resolved.index]
    max_err = float((100.0 * (fit["knudsen"] - fit["knudsen_sesam"]) / fit["knudsen_sesam"]).abs().max())
    return t.sort_values(["diameter_mm", "initial_velocity_kms"], ascending=[True, False]), C, max_err, len(resolved)


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = args.out or os.path.join(DEFAULT_PLOT_DIR, "knudsen_number.png")
    runs_dir = args.runs_dir or os.path.join(os.path.dirname(os.path.abspath(args.summary)), "runs")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    df = pd.read_csv(args.summary)
    table, C, max_err, n_fit = knudsen_table(df, runs_dir)
    table_path = os.path.join(os.path.dirname(out), "knudsen_numbers.csv")
    table.to_csv(table_path, index=False, float_format="%.6g")
    kn_min, kn_max = float(table["knudsen"].min()), float(table["knudsen"].max())
    print("{} points; lambda = C / rho with C = {:.4e} kg/m2 (from {} SESAM values, max deviation {:.2f} %); "
          "Kn {:.3g} .. {:.3g}".format(len(table), C, n_fit, max_err, kn_min, kn_max))

    diameters = sorted(table["diameter_mm"].unique())
    velocities = sorted(table["initial_velocity_kms"].unique())
    d_index = {d: i for i, d in enumerate(diameters)}
    v_index = {v: i for i, v in enumerate(velocities)}
    grid = np.full((len(diameters), len(velocities)), np.nan)
    for r in table.itertuples():
        grid[d_index[r.diameter_mm], v_index[r.initial_velocity_kms]] = r.knudsen
    d_arr, v_arr = np.array(diameters), np.array(velocities)
    dv, dd = v_arr[1] - v_arr[0], d_arr[1] - d_arr[0]
    extent = [v_arr[0] - dv / 2, v_arr[-1] + dv / 2, d_arr[0] - dd / 2, d_arr[-1] + dd / 2]

    apply_rcparams(plt)
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    lo, hi = np.log10(kn_min), np.log10(kn_max)
    cmap, edges = banded_cmap(lo, hi, KN_BANDS)
    im = ax.imshow(np.log10(grid), origin="lower", aspect="auto", extent=extent, cmap=cmap,
                   norm=Normalize(vmin=lo, vmax=hi), interpolation="nearest")
    ax.set_xlabel("initial velocity  [km/s]")
    ax.set_ylabel("initial sphere diameter  [mm]")
    ax.set_yticks(diameters[::2])
    strip_top_right_spines(ax)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    decades = [t for t in np.arange(np.ceil(lo), np.floor(hi) + 1)
               if min(abs(t - lo), abs(t - hi)) > 0.3 and all(abs(t - e) > 0.3 for e in edges)]
    ticks = sorted(set([lo, hi] + edges + list(decades)))
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(["{:.2g}".format(10 ** t) for t in ticks])
    for e in edges:
        cbar.ax.axhline(e, color="#fcfcfb", lw=1.2)
    bounds = [lo] + edges + [hi]
    names = [n for (lower, n, _l, _d) in KN_BANDS
             if (np.log10(lower) if lower > 0 else -np.inf) < hi and
                (np.log10(KN_BANDS[[b[0] for b in KN_BANDS].index(lower) + 1][0])
                 if lower != KN_BANDS[-1][0] else np.inf) > lo]
    for name, a, b in zip(names, bounds[:-1], bounds[1:]):
        cbar.ax.text(4.6, (a + b) / 2, name, transform=cbar.ax.get_yaxis_transform(),
                     va="center", ha="left", fontsize=8.5, color=SECOND)
    cbar.set_label("Knudsen number  Kn = lambda / d  (log scale, banded by flow regime)",
                   color=SECOND, fontsize=9.5, labelpad=165)
    style_colorbar(cbar)
    absent = [n.split(":")[0] for (lower, n, _l, _d) in KN_BANDS if lower > 0 and np.log10(lower) >= hi]

    fig.suptitle("Knudsen number of each sphere at its initial state", fontsize=13.5, color=INK,
                 x=0.06, ha="left", y=0.995, fontweight="bold")
    ax.set_title("lambda = mean free path of the air at the initial altitude, inherited from the parent trajectory "
                 "at that velocity\n({:.1f} km at {:.1f} km/s down to {:.2f} km at {:.3f} km/s); "
                 "lambda = {:.3g} / rho_air, calibrated on SESAM's own Knudsen numbers\n(reference length = "
                 "diameter, reproduced within {:.1f} %); colour bands = flow regimes at Kn = 0.01 / 0.1 / 10{}".format(
                     table.loc[table["initial_velocity_kms"].idxmax(), "initial_altitude_km"], v_arr[-1],
                     table.loc[table["initial_velocity_kms"].idxmin(), "initial_altitude_km"], v_arr[0],
                     C, max_err, " ({} not reached)".format(", ".join(absent)) if absent else ""),
                 fontsize=8.5, color=MUTED, loc="left", pad=8)

    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="#fcfcfb")
    print("wrote", out, "and", table_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
