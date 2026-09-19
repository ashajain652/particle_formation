#!/usr/bin/env python3
"""Run the coupled thermal model against the two no-melt US76 SESAM references in both heating modes and tabulate
the heat/temperature errors (Step 2 verification, spec section 10).

    "$PY" analysis/reentry_model_thermal_verification.py [--outdir reentry_model_output/verification_thermal]
        [--cases d100,d050] [--modes sesam,physics] [--animate]

Cases left out by --cases/--modes are read back from an existing <outdir>/<case>__<mode>.json, so summary.md /
summary.json always cover every result available in the output directory.
"""
import argparse
import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

CASES = {
    "d100": ("sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind", ["--diameter", "100", "--altitude", "77.500133"]),
    "d050": ("sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_nowind", ["--diameter", "50", "--altitude", "115"]),
}
MODES = ("sesam", "physics")
COLUMNS = ["case", "mode", "Q_conv max (of peak)", "Q_conv max (continuum, point-wise)", "integrated heat (end of hypersonic)",
           "integrated heat (end)", "max |dT_eq|", "radiated max (of peak)", "peak surface T", "runtime"]


def run_case(key, mode, outdir, animate):
    name, args = CASES[key]
    ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    label = "{}__{}".format(key, mode)
    argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                             "--heating", mode, "--reference", ref, "--outdir", outdir, "--name", label, "--quiet"]
    if animate:
        argv.append("--animate")
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run {} failed with exit code {}".format(label, rc))
    return label


def row_from_doc(key, mode, doc):
    tm, res = doc["comparison"]["thermal_metrics"], doc["results"]
    return {"case": key, "mode": mode, "Q_conv max (of peak)": "{:.2%}".format(tm["Q_conv"]["max"]),
            "Q_conv max (continuum, point-wise)": "{:.2%}".format(tm["Q_conv"]["continuum_rel_max"]),
            "integrated heat (end of hypersonic)": "{:+.2%}".format(tm["integrated_heat"]["rel_error_end_of_hypersonic"]),
            "integrated heat (end)": "{:+.2%} (ratio {:.3f})".format(tm["integrated_heat"]["rel_error_end"], tm["integrated_heat"]["ratio_end"]),
            "max |dT_eq|": "{:.1f} K ({:.2%})".format(tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"]),
            "radiated max (of peak)": "{:.2%}".format(tm["radiated"]["max"]),
            "peak surface T": "{:.0f} K at {:.0f} s".format(res["peak_surface_T_K"], res["time_of_peak_surface_T_s"]),
            "runtime": "{:.0f} s, {} steps, {:.1f} it/step".format(res["runtime_s"], res["n_macro_steps"], res["mean_newton_iterations"])}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_thermal"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--modes", default=",".join(MODES))
    p.add_argument("--animate", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    wanted = {(k, m) for k in args.cases.split(",") for m in args.modes.split(",")}
    rows = []
    for key in CASES:
        for mode in MODES:
            label = "{}__{}".format(key, mode)
            if (key, mode) in wanted:
                run_case(key, mode, args.outdir, args.animate)
            path = os.path.join(args.outdir, label + ".json")
            if os.path.isfile(path):
                with open(path) as fh:
                    rows.append(row_from_doc(key, mode, json.load(fh)))
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    lines += ["| " + " | ".join(str(r[c]) for c in COLUMNS) + " |" for r in rows]
    with open(os.path.join(args.outdir, "summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
