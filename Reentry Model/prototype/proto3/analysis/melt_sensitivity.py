#!/usr/bin/env python3
"""Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).

    "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
        [--variants base,layers2,layers6,dt025,bridged,norunoff,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]

Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
`knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate); the table
lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
change relative to `base`. Every run is a subprocess of `--python` (default: this interpreter); the `fenicsx` variant
needs the fenicsx_env interpreter (run it separately with --variants fenicsx --python <fenicsx_env python>, with CC
exported). Variants left out are read back from existing JSONs."""
import argparse
import json
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)


CASES = {"d100": ["--diameter", "100", "--altitude", "77.500133"], "d050": ["--diameter", "50", "--altitude", "115"]}
VARIANTS = {
    "base": [], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
    "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "we308": ["--we-critical", "3.08"],
    "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
    "AA7075": ["--material", "AA7075"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
    "knbody003": ["--kn-body-shock", "0.003"], "fenicsx": ["--thermal-solver", "fenicsx"],
}
KEYS = ["sprayed_mass_kg", "r_median_um", "melt_onset_altitude_km", "spraying_onset_altitude_km", "demise_altitude_km", "runtime_s"]
COLUMNS = ["case", "variant", "sprayed [kg]", "median r [um]", "melt onset [km]", "spraying onset [km]", "demise [km]", "runtime [s]"]


def argv_for(key, variant, outdir):
    return ["run"] + CASES[key] + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                                   "--melt", "on", "--heating", "physics", "--material", "AA7075_range", "--no-particles",
                                   "--outdir", outdir, "--name", "{}__{}".format(key, variant), "--quiet"] + VARIANTS[variant]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt", "sensitivity"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--variants", default=",".join(v for v in VARIANTS if v != "fenicsx"))
    p.add_argument("--python", default=sys.executable, help="interpreter for the runs (default: this one; the fenicsx variant needs fenicsx_env's)")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    for key in args.cases.split(","):
        for variant in args.variants.split(","):
            # every run in a fresh interpreter: one flight is ~250 s and a 46 k-node solver; a long in-process series was
            # seen to crawl (measured 2026-09-21)
            env = dict(os.environ, FI_PROVIDER="tcp")
            subprocess.run([args.python, "-m", "reentry_model"] + argv_for(key, variant, args.outdir), cwd=REPO_ROOT, check=True, env=env)
    rows = []
    for key in CASES:
        base = None
        for variant in VARIANTS:
            path = os.path.join(args.outdir, "{}__{}.json".format(key, variant))
            if not os.path.isfile(path):
                continue
            res = json.load(open(path))["results"]
            if variant == "base":
                base = res
            fmt = lambda k, f: "n/a" if res.get(k) is None else (f.format(res[k]) + ("" if base is None or variant == "base" or base.get(k) is None or not base[k] else " ({:+.1%})".format(res[k] / base[k] - 1.0)))
            rows.append({"case": key, "variant": variant, "sprayed [kg]": fmt("sprayed_mass_kg", "{:.4f}"), "median r [um]": fmt("r_median_um", "{:.1f}"),
                         "melt onset [km]": fmt("melt_onset_altitude_km", "{:.2f}"), "spraying onset [km]": fmt("spraying_onset_altitude_km", "{:.2f}"),
                         "demise [km]": fmt("demise_altitude_km", "{:.2f}"), "runtime [s]": fmt("runtime_s", "{:.0f}")})
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)] + ["| " + " | ".join(r[c] for c in COLUMNS) + " |" for r in rows]
    md = "\n".join(lines) + "\n"
    open(os.path.join(args.outdir, "sensitivity.md"), "w").write(md)
    json.dump(rows, open(os.path.join(args.outdir, "sensitivity.json"), "w"), indent=2)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
