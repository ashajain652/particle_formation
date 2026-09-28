# Sub-plan: Task 14 — Verification and sensitivity drivers, the reference-tier test, the runs

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 7194–7488). Read `00-shared-context.md` first.


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 14.** Every number in the verification table was measured with 4 prism layers and
`PHI_DEATH = 0.05`; both change, so the table must be re-run rather than edited.

1. **Re-measure the whole table** with the new defaults. Fact 12's layer sensitivity (sprayed mass within 0.2 %,
   demise altitude within 0.3 % across layers 2/4/6) predicts the move is small, but predicting is not measuring.
2. **The shape A/B**, which is the measurement that says what this whole design bought: the 50 mm physics flight
   with `--derived-surface off` against `on`, comparing sprayed mass, demise time and median droplet radius.
3. **The `PHI_DEATH` A/B**: 0.05 against 0.50 on the same flight, with `film_blob_fraction`,
   `film_frozen_fraction` and the stranded-film mass reported against fact 25's ranges.
4. **The remesh conservation record**: for every remesh in a flight, the pre-correction volume and energy residuals.
   Spec §10.2 — these replace the unqualified 10⁻¹⁰ balance claim of facts 28 and 31 with "exact after an explicit
   correction, pre-correction residual X".
5. **The first number to obtain, before anything else** (fact 43): how smooth the real per-patch recession field is.
   The least-squares fit is excellent at 0.03 % for a smooth field and mediocre at 22 % for a patchy one, and which
   regime the flight sits in is currently unknown.
6. **The melting references KEEP the band — re-measure them and update the README** (decided 2026-09-28).
   `tests/test_reentry_model_reference_melt.py:24` calls `mesh.sphere_mesh(ref.diameter / 2.0)` with no band, so the
   **bookkeeping device** takes `DEFAULT_BAND` and meshes the way production melting runs mesh. That is deliberate
   and it is the decision: a verification device should exercise the configuration actually run, and the band exists
   for melting runs. **Do not pin it with `band=0.0`.** Its Step 2 sibling, `test_reentry_model_reference_thermal.py`,
   *is* pinned, because a no-melt run has no business acquiring a melting-run mesh — the asymmetry is intended.

   Consequences this task must carry out, not merely note:

   * **Re-measure the six committed bookkeeping figures** and update them wherever they appear — measured fact 12 in
     `00-shared-context.md`, the README verification table, and the thresholds in
     `tests/test_reentry_model_reference*.py`. The current values were measured on the old mesh and are superseded:
     100 mm mass within 0.98 % of m₀, onset +0.10 km, 1 %-mass time −1.3 %; 50 mm 1.29 %, +0.17 km, −0.2 %; against
     thresholds 2 % / 0.5 km / 2 %. Fact 12's "no later amendment can move them" no longer holds for this device and
     that sentence must be amended too, with this decision as the reason.
   * **Expect the numbers to move very little, and treat it as a finding if they move a lot.** The device multiplies
     the conductivity by 10⁴ (`mat.k_table * 1e4`), so the body is isothermal and the mesh has almost nothing left to
     influence. If the six figures shift materially, that says the bookkeeping device is more mesh-sensitive than its
     design implies — which is worth knowing before the thesis leans on it.
   * **Record the runtime cost.** The 100 mm reference mesh goes from **18 896 nodes / 87 632 tets** to
     **33 355 / 177 363**, a factor of 2.02 in elements (measured 2026-09-27). The element-proportional part of each
     macro step — principally the conduction solve — roughly doubles for the two melting flights. Time the melting
     reference tier before and after and write both figures into `CLAUDE.md`, whose "~15 min (Step 1) + ~35 min
     (Step 2)" line does not yet cover the melting references at all.

7. **The better Design A trigger** (spec §10.3): the **melt-rate** error under a surface pinned by the latent heat.
   While melting, a temperature error becomes a mass-flux error rather than a temperature error, and mass flux is
   what feeds the spraying. The 1-D probe behind fact 44 cannot measure it; this task can.
---

**Depends on:** Task 13 (CLI) and Task 12 (metrics).
**Produces:** the reference-tier test, the verification/sensitivity driver scripts, and the runs that produce the plan's headline numbers — including the finding that a 100 mm sphere, with size feedback and the amended surface flow, stops demising and has 29.9% of its mass reach the ground.
**Character:** testing/verification/analysis — also the task with the largest real-world runtime cost in the plan.
**Read before implementing:** Measured fact 23 in the shared context: one macro time-step costs about 2.0 seconds of wall-clock time, and a full flight whose remnant survives to the ground can take up to roughly 2,357 seconds (about 40 minutes) — far longer than the spec's original under-6-minute target, which the master plan flags as needing to be restated rather than treated as already fixed. Budget real wall-clock time for whichever agent or CI system executes this task. Measured fact 24 records the 29.9%-mass-surviving finding itself, and is explicit that it is an artifact of the model's fixed body attitude (no tumbling simulated yet) — carry that caveat into every result this task reports; do not state it as a general re-entry finding.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan, and resolve the runtime-target restatement (see the shared context's note on this) before finalizing the plan's acceptance criteria.

---

### Task 14: Verification and sensitivity drivers, the reference-tier test, the runs

**Files:**
- Create: `analysis/melt_verification.py`, `analysis/melt_sensitivity.py`, `tests/test_reentry_model_reference_melt.py`
- Outputs (git-ignored): `reentry_model_output/verification_melt/{summary.md, summary.json, <case>__<mode>.*}`, `.../sensitivity/{sensitivity.md, sensitivity.json, ...}`, `.../girin/` (Task 8), `.../reference_tests/`

**Interfaces:**
- Consumes: the CLI (Task 13), `compare.melt_metrics` (Task 12), the melting references (Task 11), `sesam_io.REFERENCE_DIR`.
- Produces: the verification table (bookkeeping thresholded, resolved and physics reported), the sensitivity table, the reference-tier test with the spec §13.1 thresholds (mass 2 % of m₀, onset 0.5 km, 1 %-mass time 2 %).

- [ ] **Step 1: Create `tests/test_reentry_model_reference_melt.py`**

```python
"""Melting model vs the two melting US76 SESAM references (marker: reference, ~2 min): the bookkeeping device
(SESAM-equivalent heating, AA7075, instant removal, k x 1e4, D0/R0 kept, the Step 2 default mesh) against the acceptance
thresholds of spec section 13.1 -- mass within 2 % of the initial at every reference time, melt-onset altitude within
0.5 km, the 1 %-mass time within 2 % (measured 2026-09-21: 0.98 % / +0.10 km / -1.3 % for 100 mm, 1.29 % / +0.17 km /
-0.2 % for 50 mm). The resolved and physics-mode runs are analysis/melt_verification.py's business (reported). Metrics
are written to reentry_model_output/verification_melt/reference_tests/."""
import json
import os

import pytest

from reentry_model import aero, atmosphere, body, compare, coupled, heating, material, mesh, sesam_io, thermal
from reentry_model import trajectory as tj

NAMES = {"d100": "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
         "d050": "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"}
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "verification_melt", "reference_tests")
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}


def bookkeeping_run(ref):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    the_mesh = mesh.sphere_mesh(ref.diameter / 2.0)
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    the_body = body.MeltingBody(the_mesh, mat, thermal.thermal_solver("skfem"), body.sphere_mass(ref.diameter, ref.material_density),
                                settings=body.MeltSettings(removal="instant", runoff=False, size_feedback="initial"))
    sim = tj.Simulator(initial, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(),
                       tj.Settings(diameter=ref.diameter))
    return coupled.CoupledRun(sim, the_body, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()


@pytest.mark.reference
@pytest.mark.parametrize("key", ["d100", "d050"])
def test_bookkeeping_device_follows_sesam(key):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, NAMES[key] + ".csv"))
    hist = bookkeeping_run(ref)
    mm = compare.melt_metrics(hist, ref)
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, key + "__bookkeeping.json"), "w") as fh:
        json.dump({"case": ref.name, "results": hist.results, "melt_metrics": mm}, fh, indent=2, default=str)
    assert hist.end_reason == "demise" and abs(hist.results["melt_energy_balance_residual"]) < 1e-6
    assert mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"]
    assert abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"]
    assert abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"]
```


- [ ] **Step 2: Create `analysis/melt_verification.py`**

```python
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
`physics` -- physics-mode heating, AA7075_range, Girin removal with the size feedback (the model proper; the SESAM
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
    "physics": ["--heating", "physics", "--material", "AA7075_range", "--removal", "girin"],
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
```


- [ ] **Step 3: Create `analysis/melt_sensitivity.py`**

```python
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
```


- [ ] **Step 4: Run the reference-tier test and the verification driver**

```bash
"$PY" -m pytest -m reference tests/test_reentry_model_reference_melt.py -q          # ~1 min: both bookkeeping runs
"$PY" analysis/melt_verification.py --animate                                         # ~6 min: 2 spheres x 3 modes, with the videos
```

Expected: 2 passed (measured 0.98 % / +0.10 km / −1.3 % and 1.29 % / +0.17 km / −0.2 %); the driver prints the six-row table with `pass` in both bookkeeping rows and the resolved and physics rows of the README's verification table (Task 15), and writes `summary.md`. The whole reference tier (`"$PY" -m pytest -m reference -q`) is 15 passed in ≈ 24 min. Check that `reentry_model_output/verification_melt/d100__physics/vtk/` holds `animation.mp4`, `film.mp4`, `section.mp4` and the stills (`melt_onset`, `spraying_onset`, `peak_release`, `film_*`, `section_*`), and that `d100__resolved/mass_time.png` shows the SESAM overlay with its residual panel.

- [ ] **Step 5: Run the sensitivity study**

```bash
"$PY" analysis/melt_sensitivity.py                                                    # ~60 min: 13 variants x 2 spheres
"$PY" analysis/melt_sensitivity.py --variants fenicsx --python "$FX"                  # the backend variant from fenicsx_env (CC and FI_PROVIDER exported)
```

(export `CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang` before the second command; the script sets `FI_PROVIDER=tcp` itself.) Expected: `sensitivity.md` with every variant's sprayed mass, median radius, onsets and demise altitude and the change relative to `base`; the `fenicsx` row equal to `base` to the printed digits; `layers2`/`layers6`/`dt025` within a few percent of `base` in demise altitude and sprayed mass; `we308`, `kt±30` nearly identical (spraying is melt-limited), `kr-30` moving the median radius; `AA7075` earlier onset (850 K vs 908 K liquidus) and a leeward remnant that reaches the ground (demise n/a); `sizeinitial` ending higher (the nose stays at R₀, the body keeps the sphere's drag). `gammapm14` (`--gamma-pm 1.4`) and `knbody003` (`--kn-body-shock 0.003`) bound the two new modelling choices: the first halves the wall pressure beyond the sonic point, the second shrinks the Girin-certified band. Record the table in Task 15.

- [ ] **Step 6: Run the FEniCSx tests once more and the whole unit tier**

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
"$PY" -m pytest -m "not drama and not reference" -q
```

Expected: 8 passed; the unit tier green (~6 min).

- [ ] **Step 7: Commit**

```bash
git add analysis/melt_verification.py analysis/melt_sensitivity.py tests/test_reentry_model_reference_melt.py
git commit -m "Add the melt verification, sensitivity drivers and the reference-tier bookkeeping test (Step 3 Task 14)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

