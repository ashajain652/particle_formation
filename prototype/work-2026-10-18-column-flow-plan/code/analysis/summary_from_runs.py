#!/usr/bin/env python3
"""
summary_from_runs.py -- rebuild a sweep summary CSV for ONE material from the per-run JSONs.

sphere_sweep.py writes sweep_summary_<material>.csv itself, so this is only needed to
recover a summary that was lost or written by an older version of the sweep (which used
one generic sweep_summary.csv for every material). The per-run JSONs are never overwritten
(a non-default material has its own run-name suffix), so the same 23-column summary can
always be reconstructed from runs/<run_name>.json. Writes to the same file name the sweep
uses, so a rebuild simply replaces it.

Usage
-----
    python analysis/summary_from_runs.py                              # drama-AA7075 -> sweep_summary_AA7075.csv
    python analysis/summary_from_runs.py --material user-moltenAA7075 # -> sweep_summary_user-moltenAA7075.csv
    python analysis/summary_from_runs.py --material drama-TiAl6v4     # -> sweep_summary_TiAl6v4.csv

The default output name is sweep_summary_<material>.csv next to runs/, with DRAMA's
"drama-" prefix dropped for built-in metals.

Runs whose JSON says status != "ok" get status "failed" and blank statistics, like the
sweep's own summary. Requires only the standard library plus sphere_sweep's column lists.
"""
import argparse
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sphere_reentry as sr  # noqa: E402
import sphere_sweep as sw    # noqa: E402

DEFAULT_RUNS_DIR = os.path.join(sw.DEFAULT_OUTDIR, "runs")


def run_json_paths(runs_dir, material_name):
    """The run JSONs of one material, selected by the run-name suffix convention."""
    if material_name == sr.DEFAULT_MATERIAL.name:
        pattern = "sphere_*km.json"                    # no _m<suffix>: the default material
    else:
        pattern = "sphere_*km_m{}.json".format(sr.material_slug(material_name))
    return sorted(glob.glob(os.path.join(runs_dir, pattern)))


def summary_row_from_json(doc):
    inputs, results = doc["inputs"], doc.get("results") or {}
    ok = doc.get("status") == "ok"
    row = {c: "" for c in sw.SUMMARY_COLUMNS}
    row.update({
        "diameter_mm": inputs["diameter_mm"],
        "initial_temperature_K": inputs["initial_temperature_K"],
        "initial_velocity_kms": inputs["initial_velocity_kms"],
        "initial_altitude_km": inputs["initial_altitude_km"],
        "flight_path_angle_deg": inputs["flight_path_angle_deg"],
        "heading_deg": inputs["heading_deg"],
        "latitude_deg": inputs["latitude_deg"],
        "longitude_deg": inputs["longitude_deg"],
        "epoch_utc": inputs["epoch_utc"],
        "status": "done" if ok else "failed",
        "wall_time_s": (doc.get("provenance") or {}).get("wall_time_s", ""),
        "run_name": doc["run_name"],
    })
    if ok:
        for column in sw._RESULT_COLUMNS:
            value = results.get(column)
            row[column] = "" if value is None else value
    return row


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs-dir", default=DEFAULT_RUNS_DIR, help="directory of run JSONs (default %(default)s)")
    p.add_argument("--material", default=sr.DEFAULT_MATERIAL.name,
                   help="material name as recorded in the run JSONs (default %(default)s)")
    p.add_argument("--out", default=None,
                   help="output CSV (default sphere_sweep_output/sweep_summary_<material>.csv)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = args.out or sw.summary_path(os.path.dirname(args.runs_dir.rstrip(os.sep)), args.material)
    paths = run_json_paths(args.runs_dir, args.material)
    if not paths:
        sys.exit("no run JSONs for material {!r} under {}".format(args.material, args.runs_dir))
    rows, mismatched = [], 0
    for path in paths:
        with open(path) as fh:
            doc = json.load(fh)
        if doc["inputs"].get("material") != args.material:
            mismatched += 1
            continue
        rows.append(summary_row_from_json(doc))
    rows.sort(key=lambda r: (r["diameter_mm"], r["initial_temperature_K"], -r["initial_velocity_kms"]))
    tmp = out + ".tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=sw.SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, out)
    print("{} runs of {} ({} done, {} failed{}) -> {}".format(
        len(rows), args.material, sum(r["status"] == "done" for r in rows),
        sum(r["status"] == "failed" for r in rows),
        ", {} skipped: different material".format(mismatched) if mismatched else "", out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
