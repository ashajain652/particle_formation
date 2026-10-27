"""The sensitivity set (spec section 4.7): 15 one-at-a-time settings x 6 cases = 90 runs.

    "$PY" -m reentry_model.wake_column.analysis.sensitivity [--workers N] [--out DIR]
          [--ambient waccm|msis_afgl] [--altitudes H ...] [--time-of-day day|night ...] [--labels LABEL ...]
          [--t-end S] [--n-test-particles N]

The cases are the 14.1 micrometre droplet released at 60, 76 and 90 km, by day and by night. The 13 settings of spec
section 4.7 (a) to (h) are joined by two that the controller added after the first full run showed the 12-hour
result is dominated by uptake onto background smoke, which rests on an assumed smallest smoke bin: "smoke sticking
x0.1" and "no smoke uptake". Each setting changes one field of Settings and no other. The baseline of each case is
also run, unless a finished run is already in the folder (the default sweep's), so the tornado plot is self-contained.
`--altitudes`, `--time-of-day` and `--labels` restrict the set (testing); `--out`, `--t-end` and `--n-test-particles`
are as for sweep.py.
"""
import argparse
import dataclasses
import os
import sys

from reentry_model.wake_column.analysis import sweep
from reentry_model.wake_column.settings import Case, Settings

ALTITUDES_KM = (60.0, 76.0, 90.0)
DIAMETER_UM = 14.1


def settings_set(**overrides):
    """[(label, Settings)] for the 15 one-at-a-time settings; `overrides` (for example t_end_s) apply to all."""
    one = {
        "formation x0.1": {"formation_scale": 0.1},
        "formation x0.01": {"formation_scale": 0.01},
        "sticking x0.1": {"sticking_scale": 0.1},
        "oxidation x0.1": {"oxidation_eff": 0.1},
        "hot reactions x0.1": {"hot_reaction_scale": 0.1},
        "hot reactions x3": {"hot_reaction_scale": 3.0},
        "droplet-temperature start": {"hot_start": "droplet_temperature"},
        "formation x80": {"formation_scale": 80.0},
        "eddy x0.1": {"eddy_scale": 0.1},
        "eddy x100": {"eddy_scale": 100.0},
        "molecular only": {"eddy_scale": 0.0},
        "inner scale 10 m": {"l_inner_m": 10.0},
        "inner scale 100 m": {"l_inner_m": 100.0},
        "smoke sticking x0.1": {"smoke_sticking": 0.1},
        "no smoke uptake": {"smoke_sticking": 0.0},
    }
    base = dataclasses.replace(Settings(), **overrides)
    return [(label, dataclasses.replace(base, **change)) for label, change in one.items()]


def cases(altitudes=ALTITUDES_KM, times=sweep.TIMES_OF_DAY):
    """14.1 micrometres at 60, 76 and 90 km by day and by night (6 cases)."""
    return [Case(DIAMETER_UM, float(h), t) for h in altitudes for t in times]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sweep.common_arguments(ap)
    ap.add_argument("--altitudes", type=float, nargs="+", default=list(ALTITUDES_KM), help="release altitudes, km")
    ap.add_argument("--time-of-day", nargs="+", choices=sweep.TIMES_OF_DAY, default=list(sweep.TIMES_OF_DAY))
    ap.add_argument("--labels", nargs="+", default=None, help="only these settings (default: all 15)")
    args = ap.parse_args(argv)
    baseline = sweep.settings_from(args)
    over = {f.name: getattr(baseline, f.name) for f in dataclasses.fields(Settings)
            if getattr(baseline, f.name) != getattr(Settings(), f.name)}
    chosen = settings_set(**over)
    if args.labels is not None:
        unknown = set(args.labels) - {label for label, _ in chosen}
        if unknown:
            ap.error(f"unknown setting label(s): {sorted(unknown)}")
        chosen = [(label, st) for label, st in chosen if label in args.labels]
    os.makedirs(args.out, exist_ok=True)
    case_list = cases(args.altitudes, args.time_of_day)
    jobs = [(c, baseline) for c in case_list] + [(c, st) for c in case_list for _, st in chosen]
    print(f"{len(case_list)} cases x {len(chosen)} settings + {len(case_list)} baselines, {args.workers} worker(s), "
          f"output {args.out}", flush=True)
    _, failed = sweep.run_jobs(jobs, args.out, args.workers)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
