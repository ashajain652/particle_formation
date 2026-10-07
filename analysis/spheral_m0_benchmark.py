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

    "$PY" analysis/spheral_m0_benchmark.py --gate [--csv data/spheral/m0_benchmark.csv] [--hot-phase-s 199.5]

prints instead the decision gate of spec section 5 (plan Task 8) from the committed CSV alone, without collecting,
writing or plotting: for every damage-on run, the steps and wall time of a mass-scaled replay (spec section 9.1) over
the hot phase of the 100 mm physics flight, with and without the measured cost of the Python material classes
(commit 091c2bc). The arithmetic and its inputs are documented at GATE_* below.
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


# ---------------------------------------------------------------------------------------------------------------
# The decision gate (spec section 5; plan 2026-10-05-spheral-m0, Task 8).
#
# Hot phase of the 100 mm physics flight, melt onset to the end of spraying (Asha, 2026-10-05: the gate is applied
# to the hot phase only). No run directory of that flight exists under reentry_model_output/ (Step 3 is not in
# reentry_model; prototype/ is absent), so both ends are the Step 3 shared context's measured facts
# (docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md): melt onset 25.5 s, the start of the window
# "while the equator is intact (25.5-80 s)", which on the 50 mm flight starts at its melt onset, 174.5 s
# (sub-plan 02, AA7075_range); spraying ends at 225.0 s (the seeded whole flight to the ground at the default step,
# deep runoff off, which lands at 582.5 s). Hot phase 199.5 s.
GATE_HOT_PHASE_S = (25.5, 225.0)
GATE_FLIGHT_S = 582.5            # landing of the same run, for the whole-flight figure
# Mass scaling (spec section 9.1): the density is scaled so that sound crosses the body GATE_CROSSINGS times per
# 0.5 s interval. Sound speed and the explicit step scale alike (as 1/sqrt of the density factor), so the steps per
# crossing are the unscaled ones, D / (c_l dt), with dt Spheral's own measured step (dt_median_s) at that spacing and
# c_l the longitudinal sound speed bench.py's material gives (6,170 m/s, commit 90b38d5: K 70.3 GPa, G0 27.6 GPa,
# 2,813 kg/m^3). The crossing is counted on the 100 mm diameter of the intact sphere at the cold sound speed, so the
# shrinking, softening body makes this an upper bound on the steps (hot metal has a lower c_l and a longer dt).
GATE_DIAMETER_M = 0.100
GATE_CL_M_S = 6170.0
GATE_CROSSINGS = 10
GATE_INTERVAL_S = 0.5
GATE_WEEK_S = 7*86400.0
GATE_PYTHON_LIMIT = 0.25         # spec section 5: "more than about a quarter of a step" sends the classes to C++
# Python material classes, all three together (EOS, strength, damage), against the built-ins, 18 processes, damage
# on, median step against median step (commit 091c2bc; spheral_output/m0/subclass/summary.json):
# geometry -> (dx_mm, particles, processes, built-in median step s, Python median step s).
GATE_PYTHON = {"3d": (3.0, 19381, 18, 0.0475349, 0.0569148),
               "rz": (1.1, 3243, 18, 0.00344987, 0.00619737)}


def python_overhead(geometry):
    """The measured Python overhead: (fraction of the built-in step, seconds per step, seconds per particle-step per
    process). The second and third are the two readings used to carry it to other process counts (unmeasured)."""
    dx, n, procs, t_builtin, t_python = GATE_PYTHON[geometry]
    extra = t_python - t_builtin
    return extra/t_builtin, extra, extra/(n/procs)


def python_overhead_fit():
    """The overhead as a + b * particles per process, fitted through the 3D and the RZ measurement. Two forms, so only
    indicative (3D carries more per particle): returns (a in s per step, b in s per particle-step per process)."""
    pts = []
    for geometry in ("3d", "rz"):
        dx, n, procs, t_builtin, t_python = GATE_PYTHON[geometry]
        pts.append((n/procs, t_python - t_builtin))
    (x1, y1), (x2, y2) = pts
    b = (y1 - y2)/(x1 - x2)
    return y1 - b*x1, b


def steps_per_interval(dt):
    return GATE_CROSSINGS*(GATE_DIAMETER_M/GATE_CL_M_S)/dt


def gate_rows(rows, hot_phase_s):
    intervals = hot_phase_s/GATE_INTERVAL_S
    out = []
    for r in rows:
        if r["damage"] != "on":
            continue
        geometry, n = r["geometry"], int(r["processes"])
        t = float(r["step_wall_median_s"])
        steps = steps_per_interval(float(r["dt_median_s"]))*intervals
        frac, per_step, per_particle = python_overhead(geometry)
        ppp = int(r["particles"])/n
        fit_a, fit_b = python_overhead_fit()
        out.append(dict(geometry=geometry, dx=float(r["dx_mm"]), particles=int(r["particles"]), n=n,
                        step=t, dt=float(r["dt_median_s"]), steps_per_interval=steps/intervals,
                        steps_per_100s=steps/intervals*100.0/GATE_INTERVAL_S, steps=steps, wall=steps*t,
                        # Overhead carried to this run three ways: the measured fraction (exact at the measured run);
                        # a fixed cost per step (if it is per call, as 091c2bc argues); a cost per particle per process.
                        frac_measured=frac, frac_per_step=per_step/t, frac_per_particle=per_particle*ppp/t, frac_fit=(fit_a + fit_b*ppp)/t))
    out.sort(key=lambda g: (g["geometry"], -g["dx"], g["n"]))
    return out


def _h(seconds):
    return "{:.1f}".format(seconds/3600.0)


def print_gate(rows, hot_phase_s):
    g = gate_rows(rows, hot_phase_s)
    print("Decision gate (spec section 5), 100 mm physics flight, hot phase {:.1f} s ({}-{} s), {} crossings of "
          "{:.0f} mm at {:.0f} m/s per {} s interval; damage on".format(
              hot_phase_s, GATE_HOT_PHASE_S[0], GATE_HOT_PHASE_S[1], GATE_CROSSINGS, GATE_DIAMETER_M*1e3,
              GATE_CL_M_S, GATE_INTERVAL_S))
    head = ("geom", "dx mm", "particles", "n", "dt us", "step ms", "steps/0.5s", "steps/100s", "steps", "wall h",
            "+py meas h", "py/step", "py/ptcl", "py/fit")
    fmt = "{:<4} {:>6} {:>9} {:>3} {:>7} {:>8} {:>10} {:>10} {:>10} {:>7} {:>10} {:>7} {:>7} {:>7}"
    print(fmt.format(*head))
    for x in g:
        print(fmt.format(x["geometry"], "{:.3f}".format(x["dx"]), x["particles"], x["n"],
                         "{:.4f}".format(x["dt"]*1e6), "{:.3f}".format(x["step"]*1e3),
                         "{:.0f}".format(x["steps_per_interval"]), "{:.0f}".format(x["steps_per_100s"]),
                         "{:.3g}".format(x["steps"]), _h(x["wall"]), _h(x["wall"]*(1 + x["frac_measured"])),
                         "{:+.0%}".format(x["frac_per_step"]), "{:+.0%}".format(x["frac_per_particle"]),
                         "{:+.0%}".format(x["frac_fit"])))
    print("columns: +py meas = wall time with the measured fraction of the geometry's subclass run; py/step, py/ptcl "
          "and py/fit = the overhead carried to this run as a fixed cost per step, as a cost per particle per process, "
          "or as the fit a + b * particles per process through the 3D and RZ measurements (all three unmeasured "
          "except at the subclass runs themselves)")
    fit_a, fit_b = python_overhead_fit()
    print("fit through 3D and RZ: {:.2f} ms per step + {:.2f} us per particle-step per process".format(
        fit_a*1e3, fit_b*1e6))
    for geometry in ("3d", "rz"):
        frac, per_step, per_particle = python_overhead(geometry)
        dx, n, procs = GATE_PYTHON[geometry][:3]
        print("Python classes {} ({} mm, {} particles, {} processes): {:+.1%} of the built-in step ({:.2f} ms per "
              "step, {:.2f} us per particle-step per process) -> {}".format(
                  geometry, dx, n, procs, frac, per_step*1e3, per_particle*1e6,
                  "within" if frac <= GATE_PYTHON_LIMIT else "over"))
    three = [x for x in g if x["geometry"] == "3d" and abs(x["dx"] - 2.2) < 1e-9 and x["n"] == 18]
    if three:
        x = three[0]
        flight = x["wall"]/hot_phase_s*GATE_FLIGHT_S
        print("Gate, 3D replay at 2.2 mm ({} particles, 18 processes): {} h = {:.2f} days over the hot phase "
              "({} h with the Python classes); a week is {:.0f} h -> {}. Whole flight 0-{} s at the same rate: "
              "{} h.".format(x["particles"], _h(x["wall"]), x["wall"]/86400.0,
                             _h(x["wall"]*(1 + python_overhead("3d")[0])), GATE_WEEK_S/3600.0,
                             "under a week" if x["wall"] <= GATE_WEEK_S else "over a week", GATE_FLIGHT_S,
                             _h(flight)))
    print("Gate, Python classes: 3D {:+.1%} {} {:.0%}; RZ {:+.1%} {} {:.0%} -> split".format(
        python_overhead("3d")[0], "<=" if python_overhead("3d")[0] <= GATE_PYTHON_LIMIT else ">", GATE_PYTHON_LIMIT,
        python_overhead("rz")[0], "<=" if python_overhead("rz")[0] <= GATE_PYTHON_LIMIT else ">", GATE_PYTHON_LIMIT))
    return g


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bench-dir", default=DEFAULT_BENCH_DIR)
    p.add_argument("--csv", default=DEFAULT_CSV)
    p.add_argument("--plot-dir", default=DEFAULT_PLOT_DIR)
    p.add_argument("--gate", action="store_true",
                   help="print the decision gate from --csv only (no collection, CSV writing or plots)")
    p.add_argument("--hot-phase-s", type=float, default=GATE_HOT_PHASE_S[1] - GATE_HOT_PHASE_S[0],
                   help="hot-phase duration for --gate (default %(default)s s, 25.5-225.0 s)")
    a = p.parse_args(argv)
    if a.gate:
        if not os.path.exists(a.csv):
            print("spheral_m0_benchmark: no {}".format(a.csv), file=sys.stderr)
            return 2
        with open(a.csv) as f:
            rows = list(csv.DictReader(f))
        print_gate(rows, a.hot_phase_s)
        return 0
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
