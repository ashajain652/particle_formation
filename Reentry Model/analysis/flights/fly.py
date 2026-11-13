#!/usr/bin/env python3
"""Fly one run of the package's command line, its history kept as it goes (sub-plan 22 Task 5; ported into the package
from the prototype's planning copy, `prototype/work-2026-10-28-skin-spraying-plan/flights/fly.py`, on 2026-11-12).

    "$PY" -B analysis/flights/fly.py <outdir> <name> [reentry_model run arguments...]
    "$PY" -B analysis/flights/fly.py <outdir> <name> resume <state> [reentry_model resume options...]

Runs `python -m reentry_model run <arguments> --outdir <outdir> --name <name>` in this process, with the `reentry_model`
package of the `Reentry Model` folder above this one, or with `resume` first among the arguments (sub-plan 22 Task 5b)
`python -m reentry_model resume <state> <options> --outdir <outdir> --name <name>`, a flight continued from a state it
saved: its rows are those made after the saved step. A hook on each history row
(`reentry_model.skyscrapers.coupled.CoupledRun.row`), which only reads it (as `watch` in
analysis/reentry_model_thermal_verification.py does), appends every row to <outdir>/<name>.rows.csv as it is made and a
progress line every PROGRESS_EVERY rows to <outdir>/<name>.progress, so the results are those of the command line to the
last bit. The model writes its history and results, and saves its state as <outdir>/<name>.end.state.pkl, when a flight
ends, however it ends (since 2026-11-12): one that fails -- on a fold of the surface, a failed solve -- or is stopped by
SIGTERM or Control-C saves its last completed macro step and writes its results from the completed steps (exit 1, 143 or
130). The rows file is the record of a flight the model could not end so: one killed without warning (`kill -9`, a crash
of Python itself), or stopped after its last step or twice within a second or two (README.md). It then holds the rows to
the last one made, and <outdir>/<name>/source_rows the droplet records (sub-plan 22 Task 4's spool, removed only once the
particle files are written). It can also hold a row that <name>.csv does not: that of a step that failed after its row
was made. Exit code: the command line's, or 1 after the traceback on any other exception.
"""
import csv
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(os.path.dirname(HERE))          # the `Reentry Model` folder, which holds `reentry_model`
sys.path.insert(0, PACKAGE_ROOT)

from reentry_model import cli  # noqa: E402
from reentry_model.skyscrapers import coupled  # noqa: E402

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
        argv = args if args[:1] == ["resume"] else ["run"] + args
        return cli.main(argv + ["--outdir", outdir, "--name", name])
    except BaseException:
        traceback.print_exc()
        return 1
    finally:
        undo()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
