#!/usr/bin/env python3
"""Run the melting model against the two melting US76 SESAM references and tabulate the mass-loss errors (Step 3
verification, spec section 13.1).

    "$PY" analysis/melt_verification.py [--outdir reentry_model_output/verification_melt] [--cases d100,d050]
        [--modes bookkeeping,resolved,physics] [--animate]

Modes: `bookkeeping` -- SESAM-equivalent heating, AA7075 (DRAMA's single melting temperature), --removal instant,
--runoff off, --k-scale 1e4 (near-isothermal body): the thresholded check of the melting bookkeeping against SESAM's
lumped Q/L_f law (mass within 2 % of the initial at every reference time, onset altitude within 0.5 km, 1 %-mass
time within 2 %); `resolved` -- the same heating and material with the real conductivity, film + runoff + Girin
spraying on the default layered mesh, D0/R0 kept as SESAM keeps them (reported: the surface melts before the
interior is hot, so the mass leaves earlier and, per unit heat, the interior's sensible heating delays the end);
`physics` -- physics-mode heating, AA7075_scheil (the melting default since 2026-10-07), Girin removal with the size
feedback (the model proper; the SESAM
overlay is context, not a target). Every run writes its
overlay + residual plots and metrics JSON through the CLI; cases left out are read back from existing JSONs so
summary.md / summary.json cover everything available."""
import argparse
import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

CASES = {
    "d100": ("sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind", ["--diameter", "100", "--altitude", "77.500133"]),
    "d050": ("sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind", ["--diameter", "50", "--altitude", "115"]),
}
MODES = {
    "bookkeeping": ["--heating", "sesam", "--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--prism-layers", "0"],
    "resolved": ["--heating", "sesam", "--material", "AA7075", "--removal", "girin", "--size-feedback", "initial"],
    "physics": ["--heating", "physics", "--material", "AA7075_scheil", "--removal", "girin"],
}
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}     # bookkeeping mode only
COLUMNS = ["case", "mode", "max |dm| (of m0)", "melt onset [km] (model / SESAM)", "1 %-mass time [s] (model / SESAM)",
           "sprayed / film left [kg]", "droplets (median r)", "runtime", "verdict"]


def run_case(key, mode, outdir, animate):
    name, args = CASES[key]
    ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    label = "{}__{}".format(key, mode)
    argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                             "--melt", "on"] + MODES[mode] + ["--reference", ref, "--outdir", outdir, "--name", label, "--quiet"]
    if animate:
        argv.append("--animate")
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run {} failed with exit code {}".format(label, rc))
    return label


def row_from_doc(key, mode, doc):
    mm, res = doc["comparison"]["melt_metrics"], doc["results"]
    fmt = lambda x, f="{:.1f}": "n/a" if x is None else f.format(x)
    verdict = ""
    if mode == "bookkeeping":
        ok = (mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"] and mm["onset_altitude_diff_km"] is not None
              and abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"] and mm["demise_time_rel"] is not None
              and abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"])
        verdict = "pass" if ok else "FAIL"
    return {"case": key, "mode": mode, "max |dm| (of m0)": "{:.2%}".format(mm["mass"]["max_rel_m0"]),
            "melt onset [km] (model / SESAM)": "{} / {}".format(fmt(mm["onset_altitude_model_km"], "{:.2f}"), fmt(mm["onset_altitude_reference_km"], "{:.2f}")),
            "1 %-mass time [s] (model / SESAM)": "{} / {}{}".format(fmt(mm["demise_time_model_s"]), fmt(mm["demise_time_reference_s"]),
                                                              "" if mm["demise_time_rel"] is None else " ({:+.1%})".format(mm["demise_time_rel"])),
            "sprayed / film left [kg]": "{:.4f} / {:.4f}".format(res["sprayed_mass_kg"], res["film_mass_kg"]),
            "droplets (median r)": "{:.3g} ({} um)".format(res["n_released"], fmt(res["r_median_um"])),
            "runtime": "{:.0f} s, {} steps, {:.1f} it/step".format(res["runtime_s"], res["n_macro_steps"], res["mean_newton_iterations"]),
            "verdict": verdict}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--modes", default=",".join(MODES))
    p.add_argument("--animate", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    for key in args.cases.split(","):
        for mode in args.modes.split(","):
            run_case(key, mode, args.outdir, args.animate)
    rows = []
    for key in CASES:
        for mode in MODES:
            path = os.path.join(args.outdir, "{}__{}.json".format(key, mode))
            if os.path.isfile(path):
                doc = json.load(open(path))
                if doc.get("comparison") and "melt_metrics" in doc["comparison"]:
                    rows.append(row_from_doc(key, mode, doc))
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    lines += ["| " + " | ".join(str(r[c]) for c in COLUMNS) + " |" for r in rows]
    md = "\n".join(lines) + "\n"
    open(os.path.join(args.outdir, "summary.md"), "w").write(md)
    json.dump(rows, open(os.path.join(args.outdir, "summary.json"), "w"), indent=2)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
