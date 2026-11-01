#!/usr/bin/env python3
"""Fly one run of the planning copy's command line, its history kept as it goes (sub-plan 22 Task 5).

    "$PY" -B flights/fly.py <outdir> <name> [reentry_model run arguments...]

Runs `python -m reentry_model run <arguments> --outdir <outdir> --name <name>` in this process. A hook on each history
row, which only reads it (as `watch` in analysis/reentry_model_thermal_verification.py does), appends every row to
<outdir>/<name>.rows.csv as it is made and a progress line every PROGRESS_EVERY rows to <outdir>/<name>.progress, so the
results are those of the command line to the last bit. The model writes its history and results only when a flight
ends, and a Task 5 flight takes days with no checkpoint: one that stops early -- on a fold of the surface, a failed
solve or an interrupt -- keeps here its history to its last completed macro step, and its droplet records in
<outdir>/<name>/source_rows (sub-plan 22 Task 4's spool, removed only once the particle files are written). Exit code:
the command line's, or 1 after the traceback on any other exception.
"""
import csv
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(os.path.dirname(HERE), "code")
sys.path.insert(0, CODE)

from reentry_model import cli, coupled  # noqa: E402

PROGRESS_EVERY = 20
PROGRESS = (("time_s", "t {:8.3f} s"), ("altitude_km", "h {:7.3f} km"), ("velocity_kms", "V {:6.4f} km/s"),
            ("mass_kg", "m {:7.4f} kg"), ("temperature_K", "T {:6.1f} K"), ("surface_T_max_K", "T_surface_max {:6.1f} K"),
            ("released_mass_kg", "released {:7.4f} kg"), ("flowing_depth_mean_mm", "flowing depth {:6.3f} mm"))


def keep(outdir, name, every=PROGRESS_EVERY):
    """Hook `coupled.CoupledRun.row` to keep every row in <outdir>/<name>.rows.csv and a progress line every `every`
    rows in <outdir>/<name>.progress. Returns the function that undoes it."""
    original = coupled.CoupledRun.row
    rows_path, progress_path = (os.path.join(outdir, name + ext) for ext in (".rows.csv", ".progress"))
    state = {"n": 0, "fh": None, "writer": None, "started": time.perf_counter()}

    def row(self, t, y, a, loads):
        r = original(self, t, y, a, loads)
        if state["writer"] is None:
            state["fh"] = open(rows_path, "w", newline="")
            state["writer"] = csv.DictWriter(state["fh"], fieldnames=list(r))
            state["writer"].writeheader()
        state["writer"].writerow({k: repr(float(v)) for k, v in r.items()})
        state["fh"].flush()
        if state["n"] % every == 0:
            parts = [fmt.format(float(r[k])) for k, fmt in PROGRESS if k in r]
            with open(progress_path, "a") as fh:
                fh.write("step {:6d}  {}  wall {:9.0f} s\n".format(state["n"], "  ".join(parts),
                                                                   time.perf_counter() - state["started"]))
        state["n"] += 1
        return r

    def undo():
        coupled.CoupledRun.row = original
        if state["fh"] is not None:
            state["fh"].close()

    coupled.CoupledRun.row = row
    return undo


def main(argv):
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    outdir, name, args = os.path.abspath(argv[0]), argv[1], argv[2:]
    os.makedirs(outdir, exist_ok=True)
    undo = keep(outdir, name)
    try:
        return cli.main(["run"] + args + ["--outdir", outdir, "--name", name])
    except BaseException:
        traceback.print_exc()
        return 1
    finally:
        undo()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
