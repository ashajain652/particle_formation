"""Tornado plots of the sensitivity set (spec section 4.7) from the runs in a folder, into <out>/plots/.

    "$PY" -m reentry_model.wake_column.analysis.plot_sensitivity [--out DIR] [--ambient waccm|msis_afgl]

Reads summary.json only. For each of the 6 cases (14.1 micrometres at 60, 76 and 90 km, day and night) it compares
every one-at-a-time setting of sensitivity.settings_set() with the baseline run of the same case:

  tornado_condensed_fraction.png     change in the condensed fraction at the end of the run (an absolute change in
                                     the share of the vaporised aluminium), one bar per setting
  tornado_volume_mean_diameter.png   relative change in the volume-weighted mean diameter of the plume particles
                                     (sum of N d^4 over sum of N d^3) at the end of the run, in percent

Settings are ordered by their largest absolute effect over the cases, so every panel reads in the same order. A bar is
missing where that run is not in the folder or has no particles to measure. Runs are matched on their physics
settings, so shortened test runs and the 12-hour runs are not mixed (the longer wins).
"""
import argparse
import dataclasses
import sys

import matplotlib.pyplot as plt
import numpy as np

from reentry_model.wake_column import OUTPUT
from reentry_model.wake_column.analysis import results, sensitivity
from reentry_model.wake_column.settings import Settings

QUANTITIES = (
    ("condensed_fraction_end", "tornado_condensed_fraction.png", "change in condensed fraction",
     False),
    ("d_volume_mean_nm", "tornado_volume_mean_diameter.png",
     "change in volume-weighted mean diameter (%)", True),
)


def _value(summary, key):
    v = None if summary is None else summary.get(key)
    return None if v is None else float(v)


def effects(summaries, key, relative, ambient="waccm"):
    """({(h0_km, time_of_day): baseline value}, {label: {(h0_km, time_of_day): change}}) for one summary key.
    The change is the setting's value minus the baseline's, or (percent) relative to it when `relative`; NaN where
    either run or the value is missing."""
    base_settings = dataclasses.replace(Settings(), ambient=ambient)
    base_runs = results.pick(summaries, results.physics_signature(dataclasses.asdict(base_settings)), ambient)
    cases = sensitivity.cases()
    keys = [(c.h0_km, c.time_of_day) for c in cases]
    base = {k: _value(base_runs.get((sensitivity.DIAMETER_UM,) + k), key) for k in keys}
    table = {}
    for label, st in sensitivity.settings_set(ambient=ambient):
        runs = results.pick(summaries, results.physics_signature(dataclasses.asdict(st)), ambient)
        row = {}
        for k in keys:
            v, b = _value(runs.get((sensitivity.DIAMETER_UM,) + k), key), base[k]
            if v is None or b is None or (relative and b == 0.0):
                row[k] = np.nan
            else:
                row[k] = 100.0 * (v / b - 1.0) if relative else v - b
        table[label] = row
    return base, table


def plot_tornado(summaries, key, name, xlabel, relative, out, ambient="waccm"):
    base, table = effects(summaries, key, relative, ambient)
    keys = list(base)
    labels = list(table)
    mag = {lab: np.nanmax(np.abs([table[lab][k] for k in keys])) if not np.all(np.isnan(
        [table[lab][k] for k in keys])) else 0.0 for lab in labels}
    order = sorted(labels, key=lambda lab: mag[lab])                    # largest effect at the top of the plot
    heights = sorted({k[0] for k in keys})
    fig, axes = plt.subplots(2, len(heights), figsize=(4.6 * len(heights), 6.6), squeeze=False)
    for row, tod in enumerate(("day", "night")):
        for col, h in enumerate(heights):
            ax = axes[row, col]
            vals = np.array([table[lab][(h, tod)] for lab in order])
            colours = [results.CAT[0] if (np.isnan(v) or v >= 0) else results.CAT[1] for v in vals]
            ax.barh(range(len(order)), np.nan_to_num(vals), color=colours, height=0.7)
            ax.axvline(0.0, color=results.INK, linewidth=0.8)
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels(order)
            ax.tick_params(axis="y", labelleft=(col == 0))
            b = base[(h, tod)]
            what = "no baseline run" if b is None else (f"baseline {b:.3g} nm" if relative else f"baseline {b:.3g}")
            ax.set_title(f"{h:g} km, {tod}\n{what}", fontsize=9.5)
            for i, v in enumerate(vals):
                if np.isnan(v):
                    ax.text(0.0, i, " n/a", va="center", fontsize=7.5, color=results.MUTED)
            if row == 1:
                ax.set_xlabel(xlabel, fontsize=8.5)
            ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return results.save(fig, out, name)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=OUTPUT, help="folder holding the runs (default: the wake-column OUTPUT)")
    ap.add_argument("--ambient", choices=("waccm", "msis_afgl"), default="waccm")
    args = ap.parse_args(argv)
    results.apply_style()
    summaries = results.load_summaries(args.out)
    for key, name, xlabel, relative in QUANTITIES:
        print(plot_tornado(summaries, key, name, xlabel, relative, args.out, args.ambient))
    return 0


if __name__ == "__main__":
    sys.exit(main())
