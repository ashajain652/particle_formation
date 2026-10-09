#!/usr/bin/env python3
"""Check 3 of the Spheral M2 plan (Task 2, Step 5; spec §12): the particle body against the finite-element frame.

    "$PY" analysis/spheral_m2_check3.py [--prepared spheral_output/prepare/<prep name>] [--out spheral_output/body]
        [--csv data/spheral/m2_check3.csv] [--replay-bodies] [--force]

For every frame k of CHECK_FRAMES and every (form, spacing) of CASES, builds the body at k with
`python -m spheral_frag body --frames k-k` (resumable: a body already built is reused with its recorded timings)
and writes one CSV row per case (COLUMNS): the particle count, c, the mass, volume and heat-content relative
differences, the RZ azimuth spread and the timings. With --replay-bodies it also builds the bodies Task 10 needs,
frames 0-47 at REPLAY_CASES, and prints their size and build time.

The heat-content relative difference is undefined at frame 0 (the flight starts at a uniform 300 K, where both heat
contents are zero): its row carries NaN and the absolute difference in J. Runs on the host (drama_env)."""
import argparse
import csv
import json
import math
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from spheral_frag import naming  # noqa: E402
from spheral_frag.__main__ import main  # noqa: E402

FLIGHT = ("model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_"
          "deeprunoff-off")
DEFAULT_PREPARED = os.path.join(REPO_ROOT, "spheral_output", "prepare", "prep_{}_k00000-01254".format(FLIGHT))
DEFAULT_OUT = os.path.join(REPO_ROOT, "spheral_output", "body")
DEFAULT_CSV = os.path.join(REPO_ROOT, "data", "spheral", "m2_check3.csv")
CHECK_FRAMES = (0, 48, 100, 182, 300, 600)       # plan Task 2, Step 5: start, first film, hot, peak q, late
CASES = (("rz", 2.2), ("rz", 1.1), ("rz", 0.55), ("3d", 3.0), ("3d", 2.2), ("3d", 1.5))
REPLAY_FRAMES = (0, 47)                          # plan fact 8 and Task 10: before the first film
REPLAY_CASES = (("rz", 2.2), ("rz", 1.1), ("3d", 2.2))

COLUMNS = ["k", "time_s", "form", "dx_mm", "n_particles", "n_candidates", "c", "mass_rel_diff", "volume_rel_diff",
           "heat_rel_diff", "heat_abs_diff_J", "heat_fe_J", "fe_solid_mass_kg", "azimuth_spread_K",
           "n_outside_other_azimuths", "t_lattice_clip_s", "t_frame_s", "t_total_s", "size_bytes"]


def build(prepared, out, k0, k1, form, dx_mm, force=False):
    args = ["body", "--prepared", prepared, "--frames", "{}-{}".format(k0, k1), "--form", form, "--dx-mm",
            str(dx_mm), "--out", out, "--quiet"] + (["--force"] if force else [])
    rc = main(args)
    if rc != 0:
        raise SystemExit("body {} {} {} exited {}".format(form, dx_mm, (k0, k1), rc))
    name = naming.body_name(os.path.basename(prepared.rstrip(os.sep)), k0, k1, form, dx_mm)
    with open(os.path.join(out, name, "body.json"), encoding="utf-8") as fh:
        return json.load(fh)


def row(doc):
    c3, f = doc["check3"], doc["per_frame"][0]
    t = doc["timing_s"]
    nan = math.nan
    return {"k": f["k"], "time_s": f["time_s"], "form": doc["form"], "dx_mm": doc["dx_mm"],
            "n_particles": c3["n"], "n_candidates": c3["n_candidates"], "c": c3["c"],
            "mass_rel_diff": c3["mass_rel_diff"], "volume_rel_diff": c3["volume_rel_diff"],
            "heat_rel_diff": nan if f["heat_rel_diff"] is None else f["heat_rel_diff"],
            "heat_abs_diff_J": f["heat_abs_diff_J"], "heat_fe_J": f["heat_fe_J"],
            "fe_solid_mass_kg": f["fe_solid_mass_kg"],
            "azimuth_spread_K": f.get("azimuth_spread_K", nan) if doc["form"] == "rz" else nan,
            "n_outside_other_azimuths": f.get("n_outside_other_azimuths", "") if doc["form"] == "rz" else "",
            "t_lattice_clip_s": t["lattice_clip"], "t_frame_s": f["seconds"], "t_total_s": t["total"],
            "size_bytes": sum(v["size_bytes"] for v in doc["files"].values())}


def fmt(v):
    if isinstance(v, float):
        return "nan" if math.isnan(v) else "{:.9g}".format(v)
    return str(v)


def main_(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--prepared", default=DEFAULT_PREPARED)
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--csv", default=DEFAULT_CSV)
    p.add_argument("--replay-bodies", action="store_true")
    p.add_argument("--force", action="store_true", help="rebuild the bodies (fresh timings)")
    a = p.parse_args(argv)
    prepared = os.path.abspath(a.prepared)
    rows = []
    for k in CHECK_FRAMES:
        for form, dx in CASES:
            r = row(build(prepared, a.out, k, k, form, dx, a.force))
            rows.append(r)
            print("k {:4d} {} {:5.2f} mm  n {:7d}  c {:.6f}  mass {:+.1e}  volume {:+.4e}  heat {:+.4e}  "
                  "spread {:7.3f} K  {:.2f} s".format(k, form, dx, r["n_particles"], r["c"], r["mass_rel_diff"],
                                                     r["volume_rel_diff"], r["heat_rel_diff"],
                                                     r["azimuth_spread_K"], r["t_total_s"]), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.csv)), exist_ok=True)
    with open(a.csv, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([fmt(r[c]) for c in COLUMNS])
    print("wrote {}".format(a.csv))
    for form, dx in CASES:
        h = [abs(r["heat_rel_diff"]) for r in rows if (r["form"], r["dx_mm"]) == (form, dx)
             and not math.isnan(r["heat_rel_diff"])]
        print("worst |heat rel diff| {} {:.2f} mm: {:.3e}".format(form, dx, max(h)))
    if a.replay_bodies:
        for form, dx in REPLAY_CASES:
            doc = build(prepared, a.out, *REPLAY_FRAMES, form, dx, a.force)
            size = sum(v["size_bytes"] for v in doc["files"].values())
            print("replay body {} {:.2f} mm k {}-{}: {} particles, {:.1f} MB, {:.1f} s (frames {:.1f} s); max |heat "
                  "rel diff| {:.3e}{}".format(form, dx, *REPLAY_FRAMES, doc["n_particles"], size / 1e6,
                                             doc["timing_s"]["total"], doc["timing_s"]["frames"],
                                             doc["check3"]["max_abs_heat_rel_diff"],
                                             ", max azimuth spread {:.3f} K".format(doc["check3"]["max_azimuth_spread_K"])
                                             if form == "rz" else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main_())
