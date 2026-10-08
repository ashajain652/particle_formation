#!/usr/bin/env python3
"""Measure the Spheral M1 core on a prepared real flight (plan 2026-10-07-spheral-m1, Task 10, Step 4).

    "$PY" analysis/spheral_m1_flight.py [--prepared spheral_output/prepare/<prep name>] [--fe-run <outdir>/<run>]
        [--csv data/spheral/m1_frames.csv] [--plot-dir spheral_output/m1] [--fake-runs]

Reads the prepared directory's prepare.json (one entry per frame, written by `spheral_frag prepare`), the
finite-element run's history and, per frame, the prepared frame itself (for the release account and the frame-0
thickness against 2R), and writes:

  --csv            one row per frame (COLUMNS): the frame import, material, loads, thickness, zones, release, size and
                   timing measurements; the committed fact the flight test's thresholds come from
  <plot-dir>/      m1_drag.png (the drag three ways against the history, by flow branch), m1_frame_import.png (mass
                   and volume residuals, non-manifold counts), m1_thickness_zones.png (thinnest patch, layer depths,
                   bulk area at the two film limits), m1_release.png (sum release_rate A against the history's
                   sprayed increment), m1_cost.png (prepared frame size and time per frame)
  <plot-dir>/summary.json  the worst case of every check over the flight, and with --fake-runs the accounts of a
                   synthetic Spheral run (tests/spheral_frag_synthetic.fake_run, 3D, dx 2.2 mm) on two intervals of
                   the flight, analysed by `spheral_frag analyse`

The flow branch is the history's `flow_branch` (the finite-element surface flow model's branch: 0 shock layer, 1
merged, 2 free molecular, which the model also takes at Mach <= 1). Runs on the host (drama_env).
"""
import argparse
import csv
import glob
import json
import math
import os
import subprocess
import sys

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (REPO_ROOT, os.path.join(REPO_ROOT, "analysis"), os.path.join(REPO_ROOT, "tests")):
    if p not in sys.path:
        sys.path.insert(0, p)

import plot_style  # noqa: E402
from spheral_frag import frames  # noqa: E402

DEFAULT_FE_RUN = os.path.join(
    REPO_ROOT, "reentry_model_output",
    "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_"
    "deeprunoff-off")
DEFAULT_CSV = os.path.join(REPO_ROOT, "data", "spheral", "m1_frames.csv")
DEFAULT_PLOT_DIR = os.path.join(REPO_ROOT, "spheral_output", "m1")
FAKE_INTERVALS = ((100, 110), (800, 810))      # one in the hot phase (50-55 s), one late (400-405 s)
FAKE_DX_M = 2.2e-3
BRANCH_NAMES = {-1: "not recorded (before the first film)", 0: "shock layer", 1: "merged",
                2: "free molecular (incl. Ma <= 1)"}
BRANCH_COLOR = {-1: "#898781", 0: "#2a78d6", 1: "#eb6834", 2: "#3f9b5c"}

COLUMNS = ["k", "time_s", "history_row", "altitude_km", "velocity_kms", "mach", "knudsen", "flow_branch",
           "mesh", "n_nodes", "n_tets", "n_faces", "n_flipped", "n_open_directed_edges",
           "n_nonmanifold_edges", "n_nonmanifold_vertices", "volume_rel_diff", "mass_rel_diff", "film_rel_diff",
           "max_abs_dfl", "fl_bitwise", "n_delta_m_mismatch", "n_girin", "n_nan_forbidden",
           "loads_valid", "n_valid", "theta_last_deg", "D_hist_N", "D_patch_N", "D_table_N", "D_smooth_N",
           "rel_table_patch", "rel_patch_hist", "rel_smooth_hist", "lee_share", "smooth_lee_share",
           "thickness_min_m", "thickness_nan", "n_thin_at_2.2mm",
           "depth_normals", "n_marched", "n_wrong_way", "n_left_body", "max_slurry_depth_m", "max_liquid_depth_m",
           "facet_max_slurry_depth_m", "bulk_area_2mm_m2", "bulk_area_3mm_m2", "bulk_area_2mm_with_deep_m2",
           "bulk_area_3mm_with_deep_m2", "max_layer_m", "n_bulk_2mm", "n_runoff_2mm", "n_skin_2mm",
           "release_sum_kg", "sprayed_increment_kg", "film_mass_kg", "p_w_max_Pa", "p_w_stag_step_Pa",
           "rel_p_w_max_stag",
           "size_bytes", "t_read_s", "t_checks_s", "t_loads_s", "t_thickness_s", "t_locator_s", "t_depths_s",
           "t_write_s", "t_total_s"]


def default_prepared(fe_run):
    hits = sorted(glob.glob(os.path.join(REPO_ROOT, "spheral_output", "prepare",
                                         "prep_" + os.path.basename(fe_run) + "_k*")))
    return hits[-1] if hits else None


def frame_row(e, h, prepared):
    r = e["history_row"]
    hv = lambda col: float(h[col][r]) if col in h else float("nan")     # noqa: E731
    sprayed = h["sprayed_mass_kg"]
    row = {"k": e["k"], "time_s": e["time_s"], "history_row": r, "altitude_km": hv("altitude_km"),
           "velocity_kms": hv("velocity_kms"), "mach": hv("mach"), "knudsen": hv("knudsen"),
           "flow_branch": hv("flow_branch"), "film_mass_kg": hv("film_mass_kg"),
           "p_w_stag_step_Pa": hv("p_w_stag_step_Pa"),
           "sprayed_increment_kg": float(sprayed[r] - (sprayed[r - 1] if r else 0.0))}
    if "sha256" not in e:
        return row
    s, m, c, d, t = e["surface"], e["mass"], e["consistency"], e["drag"], e["timing_s"]
    row.update({"mesh": e["mesh"]["index"], "n_nodes": e["n_nodes"], "n_tets": e["n_tets"], "n_faces": e["n_faces"],
                "n_flipped": s["n_flipped"], "n_open_directed_edges": s["n_open_directed_edges"],
                "n_nonmanifold_edges": s["n_nonmanifold_edges"], "n_nonmanifold_vertices": s["n_nonmanifold_vertices"],
                "volume_rel_diff": s["volume_rel_diff"], "mass_rel_diff": e["mass_rel_diff"],
                "film_rel_diff": m.get("film_rel_diff"), "max_abs_dfl": c["max_abs_dfl"],
                "fl_bitwise": int(bool(c["fl_bitwise"])), "n_delta_m_mismatch": c["n_delta_m_mismatch"],
                "n_girin": c["n_girin"], "n_nan_forbidden": sum(c["nan_forbidden"].values()),
                "loads_valid": int(bool(e["loads_valid"])), "n_valid": d["n_valid"],
                "theta_last_deg": d["theta_last_deg"], "D_hist_N": d["history_N"], "D_patch_N": d["patch_N"],
                "D_table_N": d["table_N"], "D_smooth_N": d["smooth_N"], "rel_table_patch": d["rel_table_patch"],
                "rel_patch_hist": d["rel_patch_hist"], "rel_smooth_hist": d["rel_smooth_hist"],
                "lee_share": d["lee_share"], "smooth_lee_share": d["smooth_lee_share"],
                "size_bytes": e["size_bytes"]})
    for key in ("read", "checks", "loads", "thickness", "locator", "depths", "write", "total"):
        row["t_{}_s".format(key)] = t.get(key)
    if "thickness" in e:
        z, dp = e["zones"], e["depths"]
        row.update({"thickness_min_m": e["thickness"]["min_m"], "thickness_nan": e["thickness"]["n_nan"],
                    "n_thin_at_2.2mm": e["thickness"]["n_thin_at_2.2mm"], "depth_normals": dp["normals"],
                    "n_marched": dp["n_marched"], "n_wrong_way": dp["n_wrong_way"], "n_left_body": dp["n_left_body"],
                    "max_slurry_depth_m": dp["max_slurry_depth_m"], "max_liquid_depth_m": dp["max_liquid_depth_m"],
                    "facet_max_slurry_depth_m": dp.get("facet", dp)["max_slurry_depth_m"],
                    "bulk_area_2mm_m2": z["2mm"]["bulk_area_m2"], "bulk_area_3mm_m2": z["3mm"]["bulk_area_m2"],
                    "bulk_area_2mm_with_deep_m2": z["2mm_with_deep"]["bulk_area_m2"],
                    "bulk_area_3mm_with_deep_m2": z["3mm_with_deep"]["bulk_area_m2"],
                    "max_layer_m": z["max_layer_m"], "n_bulk_2mm": z["2mm"]["n_patches"]["bulk"],
                    "n_runoff_2mm": z["2mm"]["n_patches"]["runoff"], "n_skin_2mm": z["2mm"]["n_patches"]["skin"]})
    a = frames.load_prepared_arrays(os.path.join(prepared, e["file"]))
    row["release_sum_kg"] = math.fsum(a["patch_release_rate"] * a["patch_area"])
    row["p_w_max_Pa"] = float(a["patch_p_w"].max()) if len(a["patch_p_w"]) else float("nan")
    if row["p_w_stag_step_Pa"] > 0:     # the stagnation patch's p_w is the frame's largest unless it died that step
        row["rel_p_w_max_stag"] = row["p_w_max_Pa"] / row["p_w_stag_step_Pa"] - 1.0
    return row


def frame0_thickness(prepared, entries, diameter_m):
    e = next(e for e in entries if e["k"] == 0)
    a = frames.load_prepared_arrays(os.path.join(prepared, e["file"]))
    rel = a["patch_thickness"] / diameter_m - 1.0
    return {"diameter_m": diameter_m, "max_abs_rel": float(np.nanmax(np.abs(rel))), "min_rel": float(np.nanmin(rel)),
            "max_rel": float(np.nanmax(rel)), "n_nan": int(np.count_nonzero(~np.isfinite(rel)))}


def write_csv(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else
                            "{:.9g}".format(r[c]) if isinstance(r.get(c), float) else r[c]) for c in COLUMNS})


def col(rows, name):
    return np.array([np.nan if r.get(name) in (None, "") else float(r[name]) for r in rows])


def worst(rows, name, fn=np.nanmax, absolute=True, where=None):
    v = col(rows, name)
    if absolute:
        v = np.abs(v)
    if where is not None:
        v = v[where]
    v = v[np.isfinite(v)]
    return float(fn(v)) if len(v) else None


def summarise(rows, prep, frame0):
    t = col(rows, "time_s")
    valid = col(rows, "loads_valid") == 1
    branch = np.nan_to_num(col(rows, "flow_branch"), nan=-1.0)
    out = {"prepared": prep["name"], "fe_run": prep["fe_run"]["name"], "n_frames": len(rows),
           "fe_package_sha256": prep["fe_run"]["package"]["sha256"],
           "frames_without_loads": [r["k"] for r, v in zip(rows, valid) if not v],
           "frame_import": {"max_open_directed_edges": worst(rows, "n_open_directed_edges"),
                            "max_nonmanifold_edges": worst(rows, "n_nonmanifold_edges"),
                            "frames_with_nonmanifold_edges": int(np.count_nonzero(col(rows, "n_nonmanifold_edges") > 0)),
                            "max_nonmanifold_vertices": worst(rows, "n_nonmanifold_vertices"),
                            "frames_with_nonmanifold_vertices":
                                int(np.count_nonzero(col(rows, "n_nonmanifold_vertices") > 0)),
                            "max_flipped": worst(rows, "n_flipped"),
                            "max_volume_rel_diff": worst(rows, "volume_rel_diff"),
                            "max_mass_rel_diff": worst(rows, "mass_rel_diff"),
                            "max_film_rel_diff": worst(rows, "film_rel_diff"),
                            "n_mesh_versions": int(np.nanmax(col(rows, "mesh")) + 1)},
           "material": {"fl_bitwise_every_frame": bool(np.all(col(rows, "fl_bitwise") == 1)),
                        "max_abs_dfl": worst(rows, "max_abs_dfl"),
                        "table_check": prep["material_table"]},
           "consistency": {"delta_m_mismatches": int(np.nansum(col(rows, "n_delta_m_mismatch"))),
                           "nan_forbidden": int(np.nansum(col(rows, "n_nan_forbidden")))},
           "loads": {"max_abs_rel_table_patch": worst(rows, "rel_table_patch", where=valid),
                     "max_abs_rel_smooth_hist": worst(rows, "rel_smooth_hist", where=valid),
                     "max_abs_rel_patch_hist": worst(rows, "rel_patch_hist", where=valid),
                     "max_smooth_lee_share": worst(rows, "smooth_lee_share", where=valid),
                     "by_branch": {}},
           "thickness": {"frame0_vs_2R": frame0, "thinnest_m": worst(rows, "thickness_min_m", fn=np.nanmin),
                         "k_thinnest": int(rows[int(np.nanargmin(col(rows, "thickness_min_m")))]["k"]),
                         "nan_rays": int(np.nansum(col(rows, "thickness_nan")))},
           "depths": {"wrong_way": int(np.nansum(col(rows, "n_wrong_way"))),
                      "left_body": int(np.nansum(col(rows, "n_left_body"))),
                      "max_slurry_depth_m": worst(rows, "max_slurry_depth_m"),
                      "max_liquid_depth_m": worst(rows, "max_liquid_depth_m")},
           "zones": {"max_bulk_area_2mm_m2": worst(rows, "bulk_area_2mm_m2"),
                     "max_bulk_area_3mm_m2": worst(rows, "bulk_area_3mm_m2"),
                     "max_bulk_area_2mm_with_deep_m2": worst(rows, "bulk_area_2mm_with_deep_m2"),
                     "frames_with_bulk_2mm": int(np.count_nonzero(col(rows, "bulk_area_2mm_m2") > 0)),
                     "max_layer_m": worst(rows, "max_layer_m")},
           "size": {"total_bytes": int(np.nansum(col(rows, "size_bytes"))),
                    "median_frame_bytes": float(np.nanmedian(col(rows, "size_bytes"))),
                    "max_frame_bytes": float(np.nanmax(col(rows, "size_bytes")))},
           "timing": {k: {"median_s": float(np.nanmedian(col(rows, k))), "sum_s": float(np.nansum(col(rows, k)))}
                      for k in ("t_total_s", "t_thickness_s", "t_locator_s", "t_depths_s", "t_checks_s")}}
    for b in sorted(set(branch.astype(int))):
        m = valid & (branch == b)
        if not m.any():
            continue
        rel = col(rows, "rel_smooth_hist")[m]
        out["loads"]["by_branch"][BRANCH_NAMES[b]] = {
            "n_frames": int(m.sum()), "t_range_s": [float(t[m].min()), float(t[m].max())],
            "min_rel_smooth_hist": float(np.nanmin(rel)), "max_rel_smooth_hist": float(np.nanmax(rel)),
            "median_rel_smooth_hist": float(np.nanmedian(rel))}
    rs, di = col(rows, "release_sum_kg"), col(rows, "sprayed_increment_kg")
    spr = di > 0
    out["release"] = {"sum_release_kg": float(np.nansum(rs)), "sum_sprayed_increment_kg": float(np.nansum(di)),
                      "rel_total": float(np.nansum(rs) / np.nansum(di) - 1.0),
                      "rel_per_step_min": float(np.nanmin(rs[spr] / di[spr] - 1.0)),
                      "rel_per_step_max": float(np.nanmax(rs[spr] / di[spr] - 1.0)),
                      "frames_spraying": int(spr.sum()), "released_where_history_zero": int(np.count_nonzero(
                          (rs > 0) & (di <= 0)))}
    out["stagnation"] = {"max_abs_rel_p_w_max_stag": worst(rows, "rel_p_w_max_stag"),
                         "median_abs_rel_p_w_max_stag": worst(rows, "rel_p_w_max_stag", fn=np.nanmedian)}
    out["checks"] = prep["checks"]
    out["timing_prepare_s"] = prep.get("timing_s")
    out["size_prepare_bytes"] = prep.get("size_bytes")
    return out


def fake_runs(prepared, out_root):
    """A synthetic Spheral run on each FAKE_INTERVALS interval (3D, dx 2.2 mm), analysed; the accounts and rows."""
    import spheral_frag_synthetic as syn
    res = []
    runs_dir = os.path.join(out_root, "runs")
    for k0, k1 in FAKE_INTERVALS:
        try:
            run_dir, answers = syn.fake_run(prepared, k0, k1, FAKE_DX_M, seed=1, runs_dir=runs_dir)
            n_debris = syn.FAKE_N_DEBRIS
        except ValueError:      # no shell particle in the debris band on an ablated body: the ring alone
            run_dir, answers = syn.fake_run(prepared, k0, k1, FAKE_DX_M, seed=1, runs_dir=runs_dir, n_debris=0)
            n_debris = 0
        cmd = [sys.executable, "-m", "spheral_frag", "analyse", "--run", run_dir, "--prepared", prepared,
               "--outdir", os.path.join(out_root, "analyse"), "--force"]
        p = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
        entry = {"interval": [k0, k1], "run": os.path.basename(run_dir), "exit": p.returncode,
                 "n_debris_planted": n_debris, "n_particles": answers["n_particles"], "ring_mass_kg": answers["ring_mass_kg"],
                 "ring_count": len(answers["ring_ids"]),
                 "stdout_tail": p.stdout.strip().splitlines()[-3:], "stderr_tail": p.stderr.strip().splitlines()[-3:]}
        if p.returncode == 0:
            adir = os.path.join(out_root, "analyse", os.path.basename(run_dir))
            with open(os.path.join(adir, "analyse.json")) as fh:
                a = json.load(fh)
            entry.update({k: a[k] for k in ("counts", "max_abs_rel_residual", "max_abs_residual_kg", "passed")})
            with open(os.path.join(adir, "fragments.csv")) as fh:
                entry["fragment_rows"] = sum(1 for _ in fh) - 1
            with open(os.path.join(adir, "debris.csv")) as fh:
                entry["debris_rows"] = sum(1 for _ in fh) - 1
        res.append(entry)
    return res


def plot(rows, plot_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot_style.apply_rcparams(plt)
    os.makedirs(plot_dir, exist_ok=True)
    t, valid = col(rows, "time_s"), col(rows, "loads_valid") == 1
    branch = np.nan_to_num(col(rows, "flow_branch"), nan=-1.0)

    def finish(fig, axes, name):
        for ax in np.atleast_1d(axes):
            plot_style.strip_top_right_spines(ax)
            ax.grid(True, which="major", color=plot_style.GRID, linewidth=0.6)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, name), dpi=150)
        plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    ax = axes[0]
    ax.plot(t, col(rows, "D_hist_N"), color=plot_style.INK, lw=1.2, label="history (mass x load factor)")
    ax.plot(t, col(rows, "D_smooth_N"), color="#2a78d6", lw=1.0, label="frame loads on derived normals")
    ax.plot(t, col(rows, "D_patch_N"), color="#eb6834", lw=0.8, ls="--", label="frame loads on facets")
    ax.set_ylabel("drag [N]")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    for b, name in BRANCH_NAMES.items():
        m = valid & (branch == b)
        ax.plot(t[m], 100 * col(rows, "rel_smooth_hist")[m], ".", ms=2.5, color=BRANCH_COLOR[b], label=name)
    ax.plot(t[valid], 100 * col(rows, "rel_table_patch")[valid], color=plot_style.MUTED, lw=0.8,
            label="binning: table / facets - 1")
    for y in (-5, 5):
        ax.axhline(y, color=plot_style.AXIS, lw=0.8, ls=":")
    ax.set_ylabel("relative difference [%]")
    ax.set_xlabel("time [s]")
    ax.legend(frameon=False, fontsize=8, markerscale=3)
    finish(fig, axes, "m1_drag.png")

    fig, axes = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    axes[0].semilogy(t, np.maximum(col(rows, "mass_rel_diff"), 1e-17), ".", ms=2, color="#2a78d6",
                     label="mass vs history")
    axes[0].semilogy(t, np.maximum(col(rows, "volume_rel_diff"), 1e-17), ".", ms=2, color="#eb6834",
                     label="divergence vs tetrahedra volume")
    axes[0].set_ylabel("relative residual")
    axes[0].legend(frameon=False, fontsize=8, markerscale=3)
    axes[1].plot(t, col(rows, "n_nonmanifold_edges"), color="#2a78d6", lw=0.9, label="non-manifold edges")
    axes[1].plot(t, col(rows, "n_nonmanifold_vertices"), color="#eb6834", lw=0.9, label="non-manifold vertices")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].set_ylabel("count")
    axes[2].plot(t, col(rows, "n_flipped"), color=plot_style.SECOND, lw=0.9)
    axes[2].set_ylabel("faces wound inward\nin the export")
    axes[2].set_xlabel("time [s]")
    finish(fig, axes, "m1_frame_import.png")

    fig, axes = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    axes[0].plot(t, 1e3 * col(rows, "thickness_min_m"), color=plot_style.INK, lw=1)
    axes[0].set_ylabel("thinnest patch [mm]")
    axes[1].plot(t, 1e3 * col(rows, "max_slurry_depth_m"), color="#2a78d6", lw=0.9, label="slurry (f_l > 0.5)")
    axes[1].plot(t, 1e3 * col(rows, "max_liquid_depth_m"), color="#eb6834", lw=0.9, label="liquid (f_l = 1)")
    axes[1].plot(t, 1e3 * col(rows, "facet_max_slurry_depth_m"), color=plot_style.MUTED, lw=0.8, ls="--",
                 label="slurry along facet normals")
    axes[1].set_ylabel("deepest layer [mm]")
    axes[1].legend(frameon=False, fontsize=8)
    axes[2].plot(t, 1e4 * col(rows, "bulk_area_2mm_m2"), color="#2a78d6", lw=0.9, label="film limit 2 mm")
    axes[2].plot(t, 1e4 * col(rows, "bulk_area_3mm_m2"), color="#eb6834", lw=0.9, label="film limit 3 mm")
    axes[2].set_ylabel("bulk-slurry area [cm^2]")
    axes[2].set_xlabel("time [s]")
    axes[2].legend(frameon=False, fontsize=8)
    finish(fig, axes, "m1_thickness_zones.png")

    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.plot(t, 1e3 * col(rows, "sprayed_increment_kg"), color=plot_style.INK, lw=1, label="history: sprayed per step")
    ax.plot(t, 1e3 * col(rows, "release_sum_kg"), color="#2a78d6", lw=0.9, label="frame: sum release_rate A")
    ax.set_ylabel("mass per step [g]")
    ax.set_xlabel("time [s]")
    ax.legend(frameon=False, fontsize=8)
    finish(fig, ax, "m1_release.png")

    fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    axes[0].plot(t, col(rows, "size_bytes") / 1e6, color=plot_style.INK, lw=0.9)
    axes[0].set_ylabel("prepared frame [MB]")
    for key, c in (("t_total_s", plot_style.INK), ("t_thickness_s", "#2a78d6"), ("t_depths_s", "#eb6834")):
        axes[1].plot(t, col(rows, key), color=c, lw=0.8, label=key[2:-2])
    axes[1].set_ylabel("time per frame [s]")
    axes[1].set_xlabel("time [s]")
    axes[1].legend(frameon=False, fontsize=8)
    finish(fig, axes, "m1_cost.png")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fe-run", default=DEFAULT_FE_RUN)
    ap.add_argument("--prepared")
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--plot-dir", default=DEFAULT_PLOT_DIR)
    ap.add_argument("--fake-runs", action="store_true")
    args = ap.parse_args(argv)
    prepared = args.prepared or default_prepared(args.fe_run)
    if not prepared or not os.path.isfile(os.path.join(prepared, "prepare.json")):
        sys.exit("no prepared flight with a prepare.json (run spheral_frag prepare first): {}".format(prepared))
    with open(os.path.join(prepared, "prepare.json")) as fh:
        prep = json.load(fh)
    run = frames.read_fe_run(args.fe_run)
    entries = prep["frames"]
    rows = [frame_row(e, run.history, prepared) for e in entries]
    write_csv(args.csv, rows)
    summary = summarise(rows, prep, frame0_thickness(prepared, entries, prep["flight"]["diameter_m"]))
    if args.fake_runs:
        summary["fake_runs"] = fake_runs(prepared, args.plot_dir)
    os.makedirs(args.plot_dir, exist_ok=True)
    with open(os.path.join(args.plot_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
        fh.write("\n")
    plot(rows, args.plot_dir)
    print(json.dumps({k: summary[k] for k in summary if k != "checks"}, indent=1, default=float))
    print("wrote {} ({} rows) and {}".format(args.csv, len(rows), args.plot_dir))


if __name__ == "__main__":
    main()
