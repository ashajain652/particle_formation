"""Sweep plots (spec section 6) from the baseline runs in a sweep folder, into <out>/plots/.

    "$PY" -m reentry_model.wake_column.analysis.plot_sweep [--out DIR] [--ambient waccm|msis_afgl]

Reads summary.json, history.csv and size_distribution.csv of the baseline runs (default physics settings, the chosen
ambient; see results.py) and writes, for each diameter present:

  condensed_vs_time_d<D>um.png    condensed fraction of the aluminium against time, one line per release altitude,
                                  day and night panels
  provenance_d<D>um.png           end-of-run provenance stack by altitude, as a share of the inventory: by source
                                  molecule (top row) and by route (bottom row), the remaining gas in grey; day and
                                  night columns
  size_distribution_d<D>um.png    aluminium in plume particles per size bin at the end of the run, as a share of the
                                  inventory, one line per altitude, day and night panels

and across diameters:

  day_minus_night.png             condensed fraction at the end of the run by day minus by night, against altitude
  waccm_vs_msis_afgl.png          condensed fraction at the end of the run on the two ambient atmospheres (written
                                  only when NRLMSIS + AFGL runs are in the folder)

Runs that vaporised no aluminium carry no outcome and are left out. "Condensed" is everything that left the gas phase
(plume particles, background-smoke uptake and settled matter).
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FormatStrFormatter

from reentry_model.wake_column import OUTPUT
from reentry_model.wake_column.analysis import results
from reentry_model.wake_column.box import ROUTES, SOURCES

ROUTE_LABELS = {"formation": "new-particle formation", "growth": "growth of plume particles",
                "smoke": "uptake onto background smoke"}
PANELS = (("day", "Day"), ("night", "Night"))
MIN_SHARE = 1e-15                       # size-distribution points below this share of the inventory are not drawn


def _altitude_colours(runs):
    heights = sorted({h for (_, h, _) in runs})
    norm = Normalize(min(heights), max(heights) if max(heights) > min(heights) else min(heights) + 1.0)
    cmap = ListedColormap(results.ramp_colours(256))
    return norm, cmap


def _altitude_colourbar(fig, axes, norm, cmap):
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes, pad=0.02, fraction=0.04)
    cb.set_label("release altitude (km)")
    cb.outline.set_visible(False)


def plot_condensed_vs_time(runs, d, out):
    sel = {k: v for k, v in runs.items() if k[0] == d}
    norm, cmap = _altitude_colours(sel)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for ax, (tod, title) in zip(axes, PANELS):
        for (_, h, t), s in sorted(sel.items(), key=lambda kv: kv[0][1]):
            if t != tod:
                continue
            hist = results.read_history(s)
            m = hist["t_s"] > 0.0
            ax.plot(hist["t_s"][m], hist["frac_condensed"][m], color=cmap(norm(h)))
        ax.set_xscale("log")
        ax.set_xlabel("time since release (s)")
        ax.set_title(f"{title}, {d:g} µm droplet")
    axes[0].set_ylabel("condensed fraction of the vaporised aluminium")
    axes[0].set_ylim(-0.02, 1.02)
    _altitude_colourbar(fig, list(axes), norm, cmap)
    return results.save(fig, out, f"condensed_vs_time_{results.diameter_tag(d)}.png")


def _stack(ax, heights, layers, labels, colours, gas):
    bottom = np.zeros(len(heights))
    for values, label, colour in zip(layers, labels, colours):
        ax.bar(heights, values, width=1.6, bottom=bottom, color=colour, label=label, linewidth=0)
        bottom += values
    ax.bar(heights, gas, width=1.6, bottom=bottom, color="#dcdbd5", label="still gas", linewidth=0)


def plot_provenance(runs, d, out):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), sharey=True)
    for col, (tod, title) in enumerate(PANELS):
        sel = sorted((k[1], s) for k, s in runs.items() if k[0] == d and k[2] == tod)
        heights = np.array([h for h, _ in sel])
        if heights.size == 0:
            continue
        by_source = {src: np.array([sum(v for k, v in s["provenance_end"].items() if k.split("/")[0] == src)
                                    for _, s in sel]) for src in SOURCES}
        by_route = {r: np.array([sum(v for k, v in s["provenance_end"].items() if k.split("/")[1] == r)
                                 for _, s in sel]) for r in ROUTES + ("smoke",)}
        total = sum(by_source.values())
        gas = np.clip(1.0 - total, 0.0, None)
        _stack(axes[0, col], heights, [by_source[x] for x in SOURCES], list(SOURCES), results.CAT[:len(SOURCES)], gas)
        _stack(axes[1, col], heights, [by_route[r] for r in by_route], [ROUTE_LABELS[r] for r in by_route],
               results.CAT[:3], gas)
        axes[0, col].set_title(f"{title}, {d:g} µm droplet: by source molecule")
        axes[1, col].set_title(f"{title}, {d:g} µm droplet: by route")
        axes[1, col].set_xlabel("release altitude (km)")
    for row in (0, 1):
        axes[row, 0].set_ylabel("share of the vaporised aluminium")
        axes[row, 1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8.5)
    axes[0, 0].set_ylim(0.0, 1.0)
    fig.tight_layout()
    return results.save(fig, out, f"provenance_{results.diameter_tag(d)}.png")


def plot_size_distribution(runs, d, out):
    sel = {k: v for k, v in runs.items() if k[0] == d}
    norm, cmap = _altitude_colours(sel)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    drawn = False
    for ax, (tod, title) in zip(axes, PANELS):
        for (_, h, t), s in sorted(sel.items(), key=lambda kv: kv[0][1]):
            if t != tod:
                continue
            sd = results.read_size_distribution(s)
            last = sd["t_s"] == sd["t_s"].max()
            atoms = sum(v[last] for k, v in sd.items() if k.startswith("Al_"))
            share = atoms / s["N_Al"]
            keep = share > MIN_SHARE
            if keep.any():
                ax.plot(sd["d_nm"][last][keep], share[keep], marker="o", ms=3, color=cmap(norm(h)))
                drawn = True
        ax.set_xscale("log")
        ax.set_yscale("log")
        for fmt in (ax.xaxis.set_major_formatter, ax.xaxis.set_minor_formatter):
            fmt(FormatStrFormatter("%.2g"))
        ax.set_xlabel("particle diameter (nm)")
        ax.set_title(f"{title}, {d:g} µm droplet, end of run")
    axes[0].set_ylabel("aluminium in plume particles per size bin\n(share of the vaporised aluminium)")
    if not drawn:
        axes[0].text(0.5, 0.5, "no plume particles above the plotting floor", transform=axes[0].transAxes,
                     ha="center", color=results.MUTED)
    _altitude_colourbar(fig, list(axes), norm, cmap)
    return results.save(fig, out, f"size_distribution_{results.diameter_tag(d)}.png")


def plot_day_minus_night(runs, out):
    diameters = sorted({k[0] for k in runs})
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    drawn = False
    for d, colour in zip(diameters, results.CAT):
        heights = sorted({k[1] for k in runs if k[0] == d and (k[0], k[1], "day") in runs
                          and (k[0], k[1], "night") in runs})
        if not heights:
            continue
        diff = [runs[(d, h, "day")]["condensed_fraction_end"] - runs[(d, h, "night")]["condensed_fraction_end"]
                for h in heights]
        ax.plot(heights, diff, marker="o", ms=3.5, color=colour, label=f"{d:g} µm")
        drawn = True
    ax.axhline(0.0, color=results.MUTED, linewidth=0.8)
    ax.set_xlabel("release altitude (km)")
    ax.set_ylabel("condensed fraction at the end of the run,\nday minus night")
    ax.set_title("Day minus night")
    if drawn:
        ax.legend(title="droplet diameter")
    fig.tight_layout()
    return results.save(fig, out, "day_minus_night.png")


def plot_ambient_comparison(waccm, msis, out):
    diameters = sorted({k[0] for k in waccm} & {k[0] for k in msis})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for ax, (tod, title) in zip(axes, PANELS):
        for d, colour in zip(diameters, results.CAT):
            for runs, style, label in ((waccm, "-", "WACCM"), (msis, "--", "NRLMSIS + AFGL")):
                pts = sorted((k[1], s["condensed_fraction_end"]) for k, s in runs.items()
                             if k[0] == d and k[2] == tod)
                if pts:
                    ax.plot(*zip(*pts), style, marker="o", ms=3, color=colour, label=f"{d:g} µm, {label}")
        ax.set_xlabel("release altitude (km)")
        ax.set_title(title)
    axes[0].set_ylabel("condensed fraction at the end of the run")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8.5)
    fig.tight_layout()
    return results.save(fig, out, "waccm_vs_msis_afgl.png")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=OUTPUT, help="sweep folder (default: the wake-column OUTPUT)")
    ap.add_argument("--ambient", choices=("waccm", "msis_afgl"), default="waccm", help="ambient of the main plots")
    args = ap.parse_args(argv)
    results.apply_style()
    summaries = results.load_summaries(args.out)
    runs = results.baseline_runs(summaries, args.ambient)
    written = []
    for d in sorted({k[0] for k in runs}):
        written.append(plot_condensed_vs_time(runs, d, args.out))
        written.append(plot_provenance(runs, d, args.out))
        written.append(plot_size_distribution(runs, d, args.out))
    if runs:
        written.append(plot_day_minus_night(runs, args.out))
    other = "msis_afgl" if args.ambient == "waccm" else "waccm"
    others = results.baseline_runs(summaries, other)
    if runs and others:
        waccm, msis = (runs, others) if args.ambient == "waccm" else (others, runs)
        written.append(plot_ambient_comparison(waccm, msis, args.out))
    if not runs:
        print(f"no finished baseline runs with ambient {args.ambient} in {args.out}")
    for p in written:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
