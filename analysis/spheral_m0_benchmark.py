#!/usr/bin/env python3
"""Collect the Spheral M0 benchmark runs (plan 2026-10-05-spheral-m0, Task 6) into the committed table and plot them.

    "$PY" analysis/spheral_m0_benchmark.py [--bench-dir spheral_output/m0/bench] [--csv data/spheral/m0_benchmark.csv]
        [--plot-dir spheral_output/m0]

Reads every <bench-dir>/<run>/summary.json written by spheral_frag/m0/bench.py (through run_benchmark.sh), writes one
row per run to the CSV (the measured fact, committed), prints the table, and writes three plots to --plot-dir:
processor time per particle per step against particle count, parallel efficiency against process count (the two
scaling series), and peak resident memory per particle against particle count.

Parallel efficiency of a scaling series is t(1) / (n t(n)) with t the median step wall time, relative to the same
series' one-process run. Runs on the host (drama_env); nothing here imports Spheral.
"""
import argparse
import csv
import glob
import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))

import plot_style  # noqa: E402

DEFAULT_BENCH_DIR = os.path.join(REPO_ROOT, "spheral_output", "m0", "bench")
DEFAULT_CSV = os.path.join(REPO_ROOT, "data", "spheral", "m0_benchmark.csv")
DEFAULT_PLOT_DIR = os.path.join(REPO_ROOT, "spheral_output", "m0")

# The two scaling series of the plan's matrix: (geometry, dx in mm).
SCALING = (("3d", 2.2), ("rz", 0.55))

COLUMNS = ["run", "geometry", "dx_mm", "radius_mm", "particles", "particles_ratio", "processes",
           "particles_per_process_min", "particles_per_process_max", "damage", "steps", "warmup",
           "startup_s", "step_wall_median_s", "step_wall_mean_s", "step_wall_min_s", "step_wall_max_s",
           "processor_s_per_particle_step", "dt_median_s", "wall_s_per_simulated_us",
           "peak_rss_max_process_MiB", "peak_rss_total_MiB", "memory_bytes_per_particle",
           "vm_load_1m_at_start", "host_load_1m_at_start", "forced", "machine", "spheral_commit"]

# Categorical slots 1 and 2 of the dataviz reference palette (light mode): geometry is the identity. Damage is a
# secondary encoding by line style and marker, so it never rests on colour alone.
GEOMETRY_COLOR = {"3d": "#2a78d6", "rz": "#eb6834"}
GEOMETRY_LABEL = {"3d": "3D", "rz": "RZ"}
DAMAGE_STYLE = {"off": ("-", "o"), "on": ("--", "s")}


def row_from_summary(s):
    host = s.get("host_load_average_at_start")
    return {
        "run": s["run"],
        "geometry": s["geometry"],
        "dx_mm": s["dx_mm"],
        "radius_mm": s["radius_mm"],
        "particles": s["particles"],
        "particles_ratio": "{:.5f}".format(s["particles_ratio"]),
        "processes": s["processes"],
        "particles_per_process_min": s["particles_per_process_min"],
        "particles_per_process_max": s["particles_per_process_max"],
        "damage": s["damage"],
        "steps": s["steps"],
        "warmup": s["warmup"],
        "startup_s": "{:.3f}".format(s["startup_s"]),
        "step_wall_median_s": "{:.6g}".format(s["step_wall_s"]["median"]),
        "step_wall_mean_s": "{:.6g}".format(s["step_wall_s"]["mean"]),
        "step_wall_min_s": "{:.6g}".format(s["step_wall_s"]["min"]),
        "step_wall_max_s": "{:.6g}".format(s["step_wall_s"]["max"]),
        "processor_s_per_particle_step": "{:.4g}".format(s["processor_s_per_particle_step"]),
        "dt_median_s": "{:.5g}".format(s["dt_s"]["median"]),
        # Wall time to advance one microsecond of physical time at the measured median step and dt.
        "wall_s_per_simulated_us": "{:.4g}".format(s["step_wall_s"]["median"]/s["dt_s"]["median"]*1e-6),
        "peak_rss_max_process_MiB": "{:.1f}".format(s["peak_rss_bytes_max_process"]/2**20),
        "peak_rss_total_MiB": "{:.1f}".format(s["peak_rss_bytes_total"]/2**20),
        "memory_bytes_per_particle": "{:.0f}".format(s["memory_bytes_per_particle"]),
        "vm_load_1m_at_start": s["load_average_at_start"][0],
        "host_load_1m_at_start": host[0] if host else "",
        "forced": s["forced"],
        "machine": s["machine"],
        "spheral_commit": s["spheral_commit"],
    }


def collect(bench_dir):
    rows = []
    for path in glob.glob(os.path.join(bench_dir, "*", "summary.json")):
        with open(path) as f:
            rows.append(row_from_summary(json.load(f)))
    rows.sort(key=lambda r: (r["geometry"], r["damage"], -float(r["dx_mm"]), int(r["processes"])))
    return rows


def write_csv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def print_table(rows):
    head = ("geom", "dx mm", "particles", "n", "dmg", "steps", "median s", "proc-s/(particle step)", "dt s",
            "MiB/proc", "KiB/particle")
    fmt = "{:<4} {:>6} {:>9} {:>3} {:>3} {:>5} {:>9} {:>22} {:>9} {:>8} {:>12}"
    print(fmt.format(*head))
    for r in rows:
        print(fmt.format(r["geometry"], "{:.3f}".format(float(r["dx_mm"])), r["particles"], r["processes"],
                         r["damage"], r["steps"], r["step_wall_median_s"], r["processor_s_per_particle_step"],
                         r["dt_median_s"], r["peak_rss_max_process_MiB"],
                         "{:.1f}".format(float(r["memory_bytes_per_particle"])/1024)))


def is_main(r, rows):
    """True for the runs of the main matrix: the most processes run at that geometry and spacing."""
    return int(r["processes"]) == max(int(q["processes"]) for q in rows
                                      if q["geometry"] == r["geometry"] and q["dx_mm"] == r["dx_mm"]
                                      and q["damage"] == r["damage"])


def efficiency(rows, geometry, dx, damage):
    series = sorted((int(r["processes"]), float(r["step_wall_median_s"])) for r in rows
                    if r["geometry"] == geometry and abs(float(r["dx_mm"]) - dx) < 1e-9 and r["damage"] == damage)
    t1 = dict(series).get(1)
    if t1 is None:
        return []
    return [(n, t1/(n*t)) for n, t in series]


def _axes(plt):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    plot_style.strip_top_right_spines(ax)
    ax.grid(True, which="major", color=plot_style.GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def _plot_line(ax, xs, ys, geometry, damage):
    ls, marker = DAMAGE_STYLE[damage]
    ax.plot(xs, ys, linestyle=ls, marker=marker, markersize=6, linewidth=2, color=GEOMETRY_COLOR[geometry],
            markeredgecolor=plot_style.SURFACE, markeredgewidth=1.5,
            label="{}, damage {}".format(GEOMETRY_LABEL[geometry], damage))


def plot(rows, plot_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot_style.apply_rcparams(plt)
    os.makedirs(plot_dir, exist_ok=True)
    written = []

    main = [r for r in rows if is_main(r, rows)]
    for column, ylabel, name, scale in (
            ("processor_s_per_particle_step", "processor time per particle per step (s)",
             "m0_bench_cost_per_particle.png", 1.0),
            ("memory_bytes_per_particle", "peak resident memory per particle (KiB)",
             "m0_bench_memory_per_particle.png", 1/1024)):
        fig, ax = _axes(plt)
        for geometry in ("3d", "rz"):
            for damage in ("off", "on"):
                pts = sorted((int(r["particles"]), float(r[column])*scale, int(r["processes"])) for r in main
                             if r["geometry"] == geometry and r["damage"] == damage)
                if not pts:
                    continue
                _plot_line(ax, [p[0] for p in pts], [p[1] for p in pts], geometry, damage)
                if damage == "off":   # the process count of each point, once per geometry
                    for x, y, n in pts:
                        ax.annotate("n={}".format(n), (x, y), textcoords="offset points", xytext=(4, -13),
                                    fontsize=8, color=plot_style.MUTED)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("particles")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, name), dpi=150)
        plt.close(fig)
        written.append(name)

    fig, ax = _axes(plt)
    for geometry, dx in SCALING:
        for damage in ("off", "on"):
            pts = efficiency(rows, geometry, dx, damage)
            if pts:
                _plot_line(ax, [p[0] for p in pts], [p[1] for p in pts], geometry, damage)
                ax.lines[-1].set_label("{} {} mm, damage {}".format(GEOMETRY_LABEL[geometry], dx, damage))
    ax.axhline(1.0, color=plot_style.AXIS, linewidth=1)
    ax.set_xlabel("MPI processes")
    ax.set_ylabel("parallel efficiency  t(1) / (n t(n))")
    ax.set_ylim(0, 1.15)
    ax.set_xticks([1, 2, 4, 8, 12, 18])
    if ax.lines[1:]:
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "m0_bench_parallel_efficiency.png"), dpi=150)
    plt.close(fig)
    written.append("m0_bench_parallel_efficiency.png")
    return written


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bench-dir", default=DEFAULT_BENCH_DIR)
    p.add_argument("--csv", default=DEFAULT_CSV)
    p.add_argument("--plot-dir", default=DEFAULT_PLOT_DIR)
    a = p.parse_args(argv)
    rows = collect(a.bench_dir)
    if not rows:
        print("spheral_m0_benchmark: no summary.json under {}".format(a.bench_dir), file=sys.stderr)
        return 2
    write_csv(rows, a.csv)
    print_table(rows)
    for geometry, dx in SCALING:
        for damage in ("off", "on"):
            pts = efficiency(rows, geometry, dx, damage)
            if pts:
                print("efficiency {} {} mm damage {}: {}".format(
                    geometry, dx, damage, ", ".join("n={} {:.2f}".format(n, e) for n, e in pts)))
    written = plot(rows, a.plot_dir)
    print("{} runs -> {}; plots in {}: {}".format(len(rows), os.path.relpath(a.csv, REPO_ROOT),
                                                   os.path.relpath(a.plot_dir, REPO_ROOT), ", ".join(written)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
