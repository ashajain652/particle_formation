"""The default wake-column sweep: 21 release altitudes x 3 droplet diameters x day and night = 126 runs.

    "$PY" -m reentry_model.wake_column.analysis.sweep [--ambient waccm|msis_afgl] [--workers N]
          [--diameters D ...] [--altitudes H ...] [--out DIR] [--t-end S] [--n-test-particles N]

Each run goes to `<out>/<run name>/` (see run.py); a run whose summary.json is finished is not repeated, so the sweep
can be stopped and started again. `--ambient msis_afgl` runs the same cases on the cross-check atmosphere (the run
names differ, so both can share one folder). A case that raises is reported at the end and does not stop the others;
the exit code is then 1. `--out`, `--t-end` and `--n-test-particles` exist for testing and default to OUTPUT and the
Settings defaults (12 hours, 10000 test particles). Workers default to the CPU count less one.
"""
import argparse
import dataclasses
import multiprocessing
import os
import sys
import time
import traceback

from reentry_model.wake_column import OUTPUT, run
from reentry_model.wake_column.settings import Case, Settings, run_name

DIAMETERS_UM = (1.0, 14.1, 200.0)
ALTITUDES_KM = tuple(float(h) for h in range(90, 49, -2))           # 90 down to 50 km in 2 km steps
TIMES_OF_DAY = ("day", "night")


def default_cases(diameters=DIAMETERS_UM, altitudes=ALTITUDES_KM, times=TIMES_OF_DAY):
    """The sweep's cases: every diameter x altitude x time of day (126 with the defaults)."""
    return [Case(float(d), float(h), t) for d in diameters for h in altitudes for t in times]


def _worker(job):
    """Run one (case, settings, out_dir); returns ("ok", run_name, seconds) or ("failed", run_name, traceback text)."""
    case, settings, out = job
    name = run_name(case, settings)
    t0 = time.perf_counter()
    try:
        run.run(case, settings, out_dir=out)
    except Exception:                                              # one bad case must not stop a sweep
        return ("failed", name, traceback.format_exc())
    return ("ok", name, time.perf_counter() - t0)


def run_jobs(jobs, out, workers=1):
    """Run (case, settings) jobs, largest droplets first (they are the slowest, so the pool stays busy to the end).
    Returns (done, failed): done is a list of (run name, seconds); failed a list of (run name, traceback text)."""
    jobs = sorted(jobs, key=lambda j: -j[0].d0_um)
    payload = [(c, s, out) for c, s in jobs]
    done, failed = [], []
    t0 = time.perf_counter()

    def record(result, i):
        status, name, info = result
        if status == "ok":
            done.append((name, info))
            print(f"[{i}/{len(payload)}] {name}  {info:.1f} s", flush=True)
        else:
            failed.append((name, info))
            print(f"[{i}/{len(payload)}] {name}  FAILED\n{info}", flush=True)

    if workers <= 1:
        for i, job in enumerate(payload, 1):
            record(_worker(job), i)
    else:
        with multiprocessing.Pool(workers) as pool:
            for i, result in enumerate(pool.imap_unordered(_worker, payload, chunksize=1), 1):
                record(result, i)
    print(f"{len(done)} finished, {len(failed)} failed, {time.perf_counter() - t0:.0f} s wall", flush=True)
    return done, failed


def common_arguments(ap):
    ap.add_argument("--ambient", choices=("waccm", "msis_afgl"), default="waccm",
                    help="ambient atmosphere: the WACCM baseline or the NRLMSIS + AFGL cross-check")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                    help="worker processes (default: CPU count less one)")
    ap.add_argument("--out", default=OUTPUT, help="output folder (default: the wake-column OUTPUT)")
    ap.add_argument("--t-end", type=float, default=None, help="end time in seconds (testing; default 12 hours)")
    ap.add_argument("--n-test-particles", type=int, default=None,
                    help="Stage A test particles (testing; default 10000)")


def settings_from(args):
    """The baseline Settings of a command line: defaults except the ambient and any testing overrides."""
    over = {"ambient": args.ambient}
    if args.t_end is not None:
        over["t_end_s"] = args.t_end
    if args.n_test_particles is not None:
        over["n_test_particles"] = args.n_test_particles
    return dataclasses.replace(Settings(), **over)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    common_arguments(ap)
    ap.add_argument("--diameters", type=float, nargs="+", default=list(DIAMETERS_UM), help="diameters, micrometres")
    ap.add_argument("--altitudes", type=float, nargs="+", default=list(ALTITUDES_KM), help="release altitudes, km")
    args = ap.parse_args(argv)
    settings = settings_from(args)
    os.makedirs(args.out, exist_ok=True)
    cases = default_cases(args.diameters, args.altitudes)
    print(f"{len(cases)} cases, ambient {args.ambient}, {args.workers} worker(s), output {args.out}", flush=True)
    _, failed = run_jobs([(c, settings) for c in cases], args.out, args.workers)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
