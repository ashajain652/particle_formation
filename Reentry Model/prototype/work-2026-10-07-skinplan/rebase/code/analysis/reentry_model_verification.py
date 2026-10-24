#!/usr/bin/env python3
"""Run reentry_model against the four committed SESAM reference runs, twice each (SESAM density replay =
dynamics only; NRLMSISE-00 = full model), and tabulate the hypersonic-phase errors.

    "$PY" analysis/reentry_model_verification.py [--outdir reentry_model_output/verification]
        [--names <comma-separated substrings>] [--modes replay,nrlmsise]

--names/--modes restrict which cases this invocation actually runs (useful to split the eight runs
across several commands that each stay well under a timeout); cases left out are instead read back from
an existing <outdir>/<name>__<mode>.json if one is present, so summary.md/summary.json always cover every
result available in the output directory, not just the ones this invocation ran.
"""
import argparse
import json
import math
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

NAMES = [
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind",
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis",
]
MODES = ("replay", "nrlmsise")


def _row_from_doc(name, mode, ref, doc):
    m = doc["comparison"]["metrics"]
    return {"case": name, "mode": mode, "winds_in_reference": ref.use_wind,
            "hyp_dV_rel_max": m["hypersonic"]["dV_rel_max"], "hyp_dV_max_ms": m["hypersonic"]["dV_max_ms"],
            "hyp_dh_max_m": m["hypersonic"]["dh_max_m"], "all_dV_max_ms": m["all"]["dV_max_ms"], "all_dh_max_m": m["all"]["dh_max_m"],
            "d_end_time_s": m["d_end_time_s"], "d_end_time_rel": m["d_end_time_rel"],
            "runtime_s": doc["results"]["runtime_s"], "rhs_evaluations": doc["results"]["rhs_evaluations"]}


def run_case(name, mode, outdir):
    """Invoke the CLI for (name, mode), writing <outdir>/<name>__<mode>.{csv,json}, and return its row."""
    csv_path = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    ref = sesam_io.load_reference(csv_path)
    run_name = "{}__{}".format(name, mode)
    argv = ["run", "--diameter", "{:.6g}".format(ref.diameter * 1e3), "--velocity", "{:.6g}".format(ref.initial.velocity / 1e3),
            "--altitude", "{:.9g}".format(ref.initial.altitude / 1e3), "--flight-path-angle", "{:.9g}".format(math.degrees(ref.initial.flight_path)),
            "--heading", "{:.9g}".format(math.degrees(ref.initial.heading)), "--lat", "{:.9g}".format(math.degrees(ref.initial.lat)),
            "--lon", "{:.9g}".format(math.degrees(ref.initial.lon)), "--epoch", ref.initial.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
            "--material-density", "{:.6g}".format(ref.material_density),
            "--atmosphere", "replay:" + csv_path if mode == "replay" else "nrlmsise",
            "--reference", csv_path, "--outdir", outdir, "--name", run_name, "--quiet"]
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run failed for {} ({}): exit {}".format(name, mode, rc))
    with open(os.path.join(outdir, run_name + ".json")) as fh:
        doc = json.load(fh)
    return _row_from_doc(name, mode, ref, doc)


def load_case(name, mode, outdir):
    """Read back an existing <outdir>/<name>__<mode>.json without running anything; None if absent."""
    csv_path = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    run_name = "{}__{}".format(name, mode)
    json_path = os.path.join(outdir, run_name + ".json")
    if not os.path.exists(json_path):
        return None
    ref = sesam_io.load_reference(csv_path)
    with open(json_path) as fh:
        doc = json.load(fh)
    return _row_from_doc(name, mode, ref, doc)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification"))
    p.add_argument("--names", default=None, help="comma-separated substrings; a case runs if any is in its name (default: all)")
    p.add_argument("--modes", default=",".join(MODES), help="comma-separated subset of {} (default: both)".format(MODES))
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    names_filter = [s.strip() for s in args.names.split(",")] if args.names else None
    modes_filter = [s.strip() for s in args.modes.split(",")]
    for mode in modes_filter:
        if mode not in MODES:
            p.error("--modes must be a subset of {}, got {!r}".format(MODES, mode))

    rows = []
    for name in NAMES:
        for mode in MODES:
            should_run = mode in modes_filter and (names_filter is None or any(sub in name for sub in names_filter))
            row = run_case(name, mode, args.outdir) if should_run else load_case(name, mode, args.outdir)
            if row is not None:
                rows.append(row)

    lines = ["| case | mode | winds in ref | hypersonic max dV | hypersonic max dh | whole-flight max dV / dh | end time | runtime |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append("| {} | {} | {} | {:.1f} m/s ({:.3%}) | {:.0f} m | {:.1f} m/s / {:.0f} m | {:+.1f} s ({:+.2%}) | {:.0f} s, {} evals |".format(
            r["case"].replace("sphere_", "").replace("_T0300.0K_v07.50000kms", "").replace("_mAA7075_nomelt_msis", ""), r["mode"],
            "on" if r["winds_in_reference"] else "off", r["hyp_dV_max_ms"], r["hyp_dV_rel_max"], r["hyp_dh_max_m"],
            r["all_dV_max_ms"], r["all_dh_max_m"], r["d_end_time_s"], r["d_end_time_rel"], r["runtime_s"], r["rhs_evaluations"]))
    table = "\n".join(lines)
    with open(os.path.join(args.outdir, "summary.md"), "w") as fh:
        fh.write(table + "\n")
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    print(table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
