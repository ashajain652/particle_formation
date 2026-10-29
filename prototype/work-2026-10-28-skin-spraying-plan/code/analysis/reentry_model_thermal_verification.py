#!/usr/bin/env python3
"""Run the coupled thermal model against the two no-melt US76 SESAM references in both heating modes and tabulate
the heat/temperature errors (Step 2 verification, spec section 10).

    "$PY" analysis/reentry_model_thermal_verification.py [--outdir reentry_model_output/verification_thermal]
        [--cases d100,d050] [--modes sesam,physics] [--surface-model elements|skin] [--t-max S] [--summary-only]
        [--animate]

Cases left out by --cases/--modes are read back from an existing <outdir>/<case>__<mode>.json, so summary.md /
summary.json always cover every result available in the output directory.

--surface-model skin is level 3 of the layered-skin design (section 11.1, sub-plan 20 Task 10 Step 4): the same flights
with a skin on every corner of a moving 3D mesh, labelled <case>__<mode>__skin, against the same thresholds; its
summary_skin.md / summary_skin.json put each skin run beside the element run of the same case and mode when that is in
the output directory. Every run appends a line to <outdir>/<label>.progress every PROGRESS_EVERY macro steps, so a
flight of many hours can be followed. --t-max ends each flight early (a smoke check of the set-up); such a run is
labelled __t<S> and stays out of the summaries. --summary-only runs nothing and tabulates the runs already in --outdir
(the flights of one check may be run on several machines and their JSONs gathered in one directory).
"""
import argparse
import json
import os
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, coupled, sesam_io  # noqa: E402

CASES = {
    "d100": ("sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind", ["--diameter", "100", "--altitude", "77.500133"]),
    "d050": ("sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_nowind", ["--diameter", "50", "--altitude", "115"]),
}
MODES = ("sesam", "physics")
SURFACE_MODELS = ("elements", "skin")
# Step 2's acceptance thresholds for SESAM-equivalent heating, as tests/test_reentry_model_reference_thermal.py states
# them (a test keeps the two equal); physics mode is reported, bounded only by its integrated heat at the end over SESAM's
THRESHOLDS = {"Q_conv_peak_rel_max": 0.03, "Q_conv_continuum_rel_max": 0.03, "integrated_heat_rel_hypersonic": 0.03,
              "dT_rel_max": 0.02, "radiated_peak_rel_max": 0.08}
PHYSICS_RATIO = (0.4, 1.1)
PROGRESS_EVERY = 20
COLUMNS = ["case", "mode", "Q_conv max (of peak)", "Q_conv max (continuum, point-wise)", "integrated heat (end of hypersonic)",
           "integrated heat (end)", "max |dT_eq|", "radiated max (of peak)", "peak surface T", "runtime", "verdict"]
SKIN_COLUMNS = ["case", "mode", "surface"] + COLUMNS[2:] + ["skins"]


def label_of(key, mode, surface="elements", t_max=None):
    return "{}__{}".format(key, mode) + ("__skin" if surface == "skin" else "") + ("__t{:g}".format(t_max) if t_max else "")


def argv_of(key, mode, outdir, surface="elements", t_max=None, animate=False):
    name, args = CASES[key]
    ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                             "--heating", mode, "--reference", ref, "--outdir", outdir, "--name", label_of(key, mode, surface, t_max),
                             "--quiet"]
    if surface == "skin":
        argv += ["--surface-model", "skin"]
    if t_max:
        argv += ["--t-max", "{:g}".format(t_max)]
    if animate:
        argv.append("--animate")
    return argv


def watch(path, every=PROGRESS_EVERY):
    """Make every coupled run append a line to `path` every `every` macro steps: flight time, altitude, speed, mean
    and hottest surface temperature, and wall time since the run began. It only reads each step's history row, so the
    results are unchanged. Returns the function that undoes it."""
    original = coupled.CoupledRun.row
    seen, started = [0], time.perf_counter()

    def row(self, t, y, a, loads):
        r = original(self, t, y, a, loads)
        if seen[0] % every == 0:
            with open(path, "a") as fh:
                fh.write("step {:6d}  t {:8.2f} s  h {:8.3f} km  V {:7.4f} km/s  T {:7.1f} K  T_surface_max {:7.1f} K  "
                         "wall {:9.0f} s\n".format(seen[0], r["time_s"], r["altitude_km"], r["velocity_kms"],
                                                   r["temperature_K"], r["surface_T_max_K"], time.perf_counter() - started))
        seen[0] += 1
        return r

    coupled.CoupledRun.row = row
    return lambda: setattr(coupled.CoupledRun, "row", original)


def run_case(key, mode, outdir, animate, surface="elements", t_max=None):
    label = label_of(key, mode, surface, t_max)
    undo = watch(os.path.join(outdir, label + ".progress"))
    try:
        rc = cli.main(argv_of(key, mode, outdir, surface, t_max, animate))
    finally:
        undo()
    if rc != 0:
        raise SystemExit("run {} failed with exit code {}".format(label, rc))
    return label


def verdict(mode, tm, results):
    """SESAM-equivalent heating: `pass` or `FAIL (<the metrics outside their threshold>)` against THRESHOLDS, the
    energy books and the ground; physics mode is reported, checked only as the reference test checks it."""
    if results.get("end_reason") != "ground":
        return "FAIL (ended by {})".format(results.get("end_reason"))
    if mode == "physics":
        ratio = tm["integrated_heat"]["ratio_end"]
        return "reported" if PHYSICS_RATIO[0] < ratio < PHYSICS_RATIO[1] else "FAIL (heat ratio {:.3f})".format(ratio)
    measured = {"Q_conv_peak_rel_max": tm["Q_conv"]["max"], "Q_conv_continuum_rel_max": tm["Q_conv"]["continuum_rel_max"],
                "integrated_heat_rel_hypersonic": abs(tm["integrated_heat"]["rel_error_end_of_hypersonic"]),
                "dT_rel_max": tm["temperature"]["dT_rel_max"], "radiated_peak_rel_max": tm["radiated"]["max"]}
    over = [k for k, v in measured.items() if not v <= THRESHOLDS[k]]
    if not abs(results["energy_balance_residual"]) < 1e-6:
        over.append("energy_balance_residual")
    return "pass" if not over else "FAIL ({})".format(", ".join(over))


def row_from_doc(key, mode, doc, surface="elements"):
    tm, res = doc["comparison"]["thermal_metrics"], doc["results"]
    row = {"case": key, "mode": mode, "surface": surface, "Q_conv max (of peak)": "{:.2%}".format(tm["Q_conv"]["max"]),
           "Q_conv max (continuum, point-wise)": "{:.2%}".format(tm["Q_conv"]["continuum_rel_max"]),
           "integrated heat (end of hypersonic)": "{:+.2%}".format(tm["integrated_heat"]["rel_error_end_of_hypersonic"]),
           "integrated heat (end)": "{:+.2%} (ratio {:.3f})".format(tm["integrated_heat"]["rel_error_end"], tm["integrated_heat"]["ratio_end"]),
           "max |dT_eq|": "{:.1f} K ({:.2%})".format(tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"]),
           "radiated max (of peak)": "{:.2%}".format(tm["radiated"]["max"]),
           "peak surface T": "{:.0f} K at {:.0f} s".format(res["peak_surface_T_K"], res["time_of_peak_surface_T_s"]),
           "runtime": "{:.0f} s, {} steps, {:.1f} it/step".format(res["runtime_s"], res["n_macro_steps"], res["mean_newton_iterations"]),
           "verdict": verdict(mode, tm, res), "skins": ""}
    if "skin_columns" in res:
        row["skins"] = ("{} columns, {} rebuilds, {} repeats, {} base faults; {:.0f} s in the skins, {:.0f} s in the 3D solve, "
                        "{:.1f} s a step".format(res["skin_columns"], res["rebuilds"], res["interface_repeats"], res["jump_faults"],
                                                 res["wall_skins_s"], res["wall_3d_s"], res["runtime_s"] / max(res["n_macro_steps"], 1)))
    return row


def read_doc(outdir, label):
    path = os.path.join(outdir, label + ".json")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        doc = json.load(fh)
    if not doc.get("comparison"):
        print("skipping {}: no comparison block (not a run of this script)".format(label))
        return None
    return doc


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_thermal"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--modes", default=",".join(MODES))
    p.add_argument("--surface-model", choices=SURFACE_MODELS, default="elements")
    p.add_argument("--t-max", type=float, default=None, help="end each flight at this time [s]: a smoke check, kept out of "
                                                              "the summaries")
    p.add_argument("--summary-only", action="store_true", help="run nothing: tabulate the runs already in --outdir")
    p.add_argument("--animate", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    wanted = {(k, m) for k in args.cases.split(",") for m in args.modes.split(",")}
    if args.summary_only and args.t_max is not None:
        p.error("--summary-only and --t-max exclude each other")
    unknown = {k for k, _ in wanted if k not in CASES} | {m for _, m in wanted if m not in MODES}
    if unknown:
        p.error("unknown case or mode: {}".format(", ".join(sorted(unknown))))
    if args.t_max is not None:
        if args.t_max <= 0.0:
            p.error("--t-max must be > 0")
        for key, mode in sorted(wanted):
            run_case(key, mode, args.outdir, args.animate, args.surface_model, args.t_max)
        return 0
    skins = args.surface_model == "skin"
    rows = []
    for key in CASES:
        for mode in MODES:
            if (key, mode) in wanted and not args.summary_only:
                run_case(key, mode, args.outdir, args.animate, args.surface_model)
            for surface in (("elements", "skin") if skins else ("elements",)):
                doc = read_doc(args.outdir, label_of(key, mode, surface))
                if doc is not None:
                    rows.append(row_from_doc(key, mode, doc, surface))
    columns = SKIN_COLUMNS if skins else COLUMNS
    stem = "summary_skin" if skins else "summary"
    with open(os.path.join(args.outdir, stem + ".json"), "w") as fh:
        json.dump([{c: r[c] for c in columns} for r in rows], fh, indent=2)
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    lines += ["| " + " | ".join(str(r[c]) for c in columns) + " |" for r in rows]
    if any(r["mode"] == "physics" for r in rows):
        lines.append("")
        lines.append("physics-mode power errors are differences from SESAM's surface-averaged model, not thresholds")
    with open(os.path.join(args.outdir, stem + ".md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
