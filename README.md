# Sphere fragment re-entry sweep (ESA DRAMA / SARA-SESAM)

Models the atmospheric descent of solid AA7075 spheres that break off a
re-entering satellite, using SESAM (the re-entry module of ESA DRAMA's SARA)
through the pyDRAMA package.

- `sphere_reentry.py` — one sphere from a given initial state → `runs/<run_name>.csv` (full history) + `runs/<run_name>.json` (statistics: maximum temperature, final mass, time at melting temperature, final velocity, final radius, and more).
- `sphere_sweep.py` — runs `sphere_reentry.py` over diameters 5–100 mm, initial temperatures 300–750 K and 100 velocities from 7.5 km/s down to the parent satellite's impact velocity. Altitude, latitude, longitude, flight-path angle, heading and epoch are inherited from the parent run (`Generic_Satellite Reentry/output/dmf_output.json`) at the point where the parent reached that velocity.

Design: `docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md`.

## Environment

```bash
# pyDRAMA lives in the conda env drama_env (Python 3.12); DRAMA 4.1.4 is at /Applications/DRAMA-4.1.4
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
"$PY" -m pip install tqdm pytest          # once
```

`data/fap_day.dat` and `data/fap_mon.dat` are the space-weather files of the parent run.

## One sphere

```bash
"$PY" sphere_reentry.py --velocity 7.5 --altitude 77.5 --temperature 300 --diameter 50 \
    --flight-path-angle -0.96 --heading 347.2 --lat 29.5 --lon -82.1 --epoch 2024-08-01T12:53:07
```

Required: `--velocity` km/s, `--altitude` km, `--temperature` K, `--diameter` mm.
Optional state values default to 0°/0°/0°/0° and the parent epoch. `--dry-run` prints the SESAM
configuration; `--keep-raw` keeps the raw DRAMA tree under `sphere_sweep_output/raw/<run_name>/`.
`--atmosphere static|nrlmsise` selects SESAM's environment (default `static` = DRAMA's US76 table;
`nrlmsise` = NRLMSISE-00 with F10.7/Ap from the fap files and appends `_msis` to the run name).
`--no-wind` runs SESAM without HWM14 winds and appends `_nowind`.

**Environment actually used by SESAM.** pyDRAMA copies boolean settings into `sara.xml` as `True`/`False`,
but SESAM only understands `yes`/`no`, so the parent's `dynamicEnvironment: True` was never honoured: every
run of the sweep used DRAMA's **static environment table** (`TOOLS/SARA/REENTRY/data/StaticEnvironmentData.csv`
= US Standard Atmosphere 1976 with a few m/s of fixed wind; `sesam.log` says "Using static environment from
CSV file") and **no winds** (SESAM's fallback for the unparseable flag; verified 2026-09-17 by identical
histories). The three switches (`dynamicEnvironment`, `useWind`, and `useEnvironmentCSV`, which pyDRAMA
writes as `solarActivity/valuesFromFile`) are now written as `yes`/`no` strings, and the default keeps the
sweep's static, wind-free behaviour explicitly. `sesam.log` reports the choice
("Using dynamic atmosphere model NRLMSISI-00 with dynamic solar activity", "Using dynamic wind model HWM14",
"Simulating without wind").

## Reference runs for the physics model (`sphere_sweep_output/reference_AA7075_nomelt/`)

`data/user_materials/AA7075_nomelt.json` is `drama-AA7075` with the melting temperature raised to 1e5 K
(cp and k tables held at their 850 K values above 850 K), so SESAM never removes mass. Two spheres, both
300 K, 7.5 km/s, γ = −0.959°, heading 347.2°, 29.55° N, 82.13° W, epoch 2024-08-01T12:53:07 (the parent's
break-off state), NRLMSISE-00, winds on and off:

- **100 mm from 77.5 km** — continuum-dominated (Kn 0.03 → 0); ground impact at 368 s.
- **50 mm from 115 km** — starts free-molecular (Kn = 41), Kn < 10 at 106.5 km, < 1 at 93.8 km, < 0.1 at 80.4 km,
  < 0.01 at 65 km; peak heating 9.2 kW at 63 km; ground impact at 561 s.

```bash
"$PY" sphere_reentry.py --velocity 7.5 --altitude 77.500133 --temperature 300 --diameter 100 \
    --flight-path-angle -0.959331 --heading 347.168296 --lat 29.546067 --lon -82.134333 \
    --epoch 2024-08-01T12:53:07 --material-file data/user_materials/AA7075_nomelt.json --atmosphere nrlmsise \
    --outdir sphere_sweep_output/reference_AA7075_nomelt/runs --raw-dir sphere_sweep_output/reference_AA7075_nomelt/raw --keep-raw [--no-wind]
"$PY" sphere_reentry.py --velocity 7.5 --altitude 115 --temperature 300 --diameter 50  ...same options...
```

Run names: `sphere_d100.00mm_..._h077.500km_mAA7075_nomelt_msis[_nowind]` and
`sphere_d050.00mm_..._h115.000km_mAA7075_nomelt_msis[_nowind]`; raw DRAMA trees kept under `raw/`.
The same 100 mm case on the static US76 table (`..._mAA7075_nomelt[_nowind]`) is kept for comparison:
NRLMSISE-00 differs from US76 by −4 % to +4 % in density along that flight (≤ 57 m/s in velocity,
≤ 0.4 km in altitude); HWM14 winds change either sphere's trajectory by ≤ 8 m/s and ≤ 20 m.
Exit codes: 0 ok, 1 the run failed (see the JSON's `error`), 2 bad arguments or pyDRAMA missing.

## The sweep

```bash
"$PY" sphere_sweep.py --dry-run                 # matrix, skipped points, batch plan; runs nothing
"$PY" sphere_sweep.py                           # one batch of everything pending: confirm, choose cores
"$PY" sphere_sweep.py --batch-size 5000         # confirm + choose cores before every batch
"$PY" sphere_sweep.py --yes --cores 7           # unattended
"$PY" sphere_sweep.py --diameters 50 --temperatures 300 --velocities 7.5,0.5 --yes --cores 2   # subset
```

Every point's `runs/<run_name>.json` makes it resumable: re-running skips completed points,
`--retry-failed` re-runs failed ones, `--force` re-runs everything. A tqdm bar shows the current
run, elapsed time, rate and ETA. Ctrl-C finishes the running subprocesses, saves the manifest and
exits with code 130. Outputs in `sphere_sweep_output/`: `runs/`, `raw/` (failed runs only),
`sweep_manifest_<material>.json`, `sweep_summary_<material>.csv` (one row per matrix point) and
`sweep.log`. `<material>` is the sweep's material without DRAMA's `drama-` prefix —
`sweep_summary_AA7075.csv` by default — so sweeps of different materials never overwrite each
other's manifest or summary.

The velocity grid's lower bound is the parent's minimum velocity, read from the parent file and
logged as `v_min = ... km/s at t = ... s, altitude ... km (parent row N)`.

## Materials

Every sphere is `drama-AA7075` unless you say otherwise. Both scripts accept the same
pair of mutually exclusive flags:

```bash
"$PY" sphere_reentry.py ... --material drama-TiAl6v4                       # another DRAMA metal
"$PY" sphere_reentry.py ... --material-file examples/material_al_li_2195.json   # a custom metal
"$PY" sphere_sweep.py --material drama-TiAl6v4 --yes --cores 6              # whole sweep in titanium
```

`--material NAME` must be one of the 21 metals in DRAMA's own database
(`/Applications/DRAMA-4.1.4/TOOLS/material_database.xml`); the sphere's density and melting
temperature are read from there, so mass, final radius and the time-at-melting-temperature
statistic follow the material. A typo lists the valid names: drama-A316, drama-A316_Semi-Empirical,
drama-AA7075, drama-Bat-Li, drama-Bat-NiCd, drama-Beryllium, drama-Carbon-Carbon, drama-Copper,
drama-El-Mat, drama-HC-AA7075, drama-HC-CFRP-4ply, drama-HC-CFRP-8ply, drama-HiperCo,
drama-Inconel718, drama-Inermet, drama-Invar, drama-Iron, drama-SiC, drama-SolarPanel-Mat,
drama-TiAl6v4, drama-Tungsten.

`--material-file PATH` supplies a custom metal as JSON in DRAMA's material format — see
[examples/material_al_li_2195.json](examples/material_al_li_2195.json) (the Al-Li alloy of the
Falcon 9 study) for a complete template. Required fields: `name, density, specificHeatCapacity,
meltingHeat, meltingTemperature, emissivity, heatConductivity, oxideActivationTemperature,
oxideEmissivity, oxideHeatOfFormation, oxideReactionProbability`; curves are `[temperature K, value]`
pairs, a single pair means constant; keys starting with `_` are ignored (use them for notes).
The definition is injected into SESAM as its `materialList` and recorded in the run JSON.

A non-default material appends `_m<name>` to every run name
(`..._h077.500km_mdrama-TiAl6v4`), so sweeps of different materials can share
`sphere_sweep_output/` without overwriting each other and the finished AA7075 sweep stays
resumable. The sweep's material is one setting per invocation (it is not a grid axis) and is
recorded in `sweep_manifest_<material>.json` under `settings.material`.

## Analysis / plotting

```bash
"$PY" -m pip install pandas   # once; not needed by the two main scripts, only by analysis/
"$PY" analysis/plot_outcome_vs_diameter_velocity.py                    # demise fraction, diameter x velocity
"$PY" analysis/plot_outcome_vs_diameter_velocity.py --metric mass_loss # mean mass-loss fraction instead
"$PY" analysis/plot_outcome_vs_diameter_velocity.py --metric mass_loss --vmax auto   # color range = data range, not 0-1
"$PY" analysis/plot_outcome_by_diameter.py                             # one panel per sphere size (temperature x velocity)
"$PY" analysis/plot_outcome_by_diameter.py --metric mass_loss
```

```bash
"$PY" analysis/plot_max_temperature_by_diameter.py                                  # AA7075 by default
"$PY" analysis/plot_outcome_by_diameter.py --summary sphere_sweep_output/sweep_summary_user-moltenAA7075.csv
"$PY" analysis/summary_from_runs.py --material user-moltenAA7075   # rebuild a summary from runs/*.json if ever lost
"$PY" analysis/plot_weber_number.py --summary sphere_sweep_output/T850_moltenAA7075/sweep_summary_user-moltenAA7075.csv --sigma 0.809
"$PY" analysis/plot_ohnesorge_number.py --summary sphere_sweep_output/T850_moltenAA7075/sweep_summary_user-moltenAA7075.csv --mu 1.2e-3 --rho-liquid 2400
"$PY" analysis/plot_knudsen_number.py --summary sphere_sweep_output/T850_moltenAA7075/sweep_summary_user-moltenAA7075.csv
```

The plot scripts read `sphere_sweep_output/sweep_summary_AA7075.csv` by default and write PNGs
to `sphere_sweep_output/plots/` (`--summary`/`--out` override either path — pass another
material's `sweep_summary_<material>.csv` to plot it). `--metric outcome` (default) is the
binary `outcome == "demised"` result per run; `--metric mass_loss` uses the continuous
`mass_loss_fraction` column instead — see the module docstrings for what each represents.
`plot_outcome_vs_diameter_velocity.py` averages over however many initial temperatures the
summary holds (a single-temperature sweep such as `--temperatures 850` is plotted as-is) and
`--vmax auto` stretches its color scale to the data instead of the fixed 0–1 that keeps
materials comparable. To keep a re-run's aggregates separate from an earlier sweep of the same
material, give it its own `--outdir` (e.g. `sphere_sweep_output/T850_moltenAA7075`).
`plot_max_temperature_by_diameter.py` shows, per sphere size, the maximum temperature reached
by spheres that lost no mass, with every mass-losing cell black. `plot_weber_number.py` maps
We = ρ v² d / σ at each sphere's initial state (diameter × velocity, log color banded by the breakup
regimes We < 12 / 12–50 / 50–100 / 100–350 / ≥ 350 in blue / yellow / orange / red / pink), with ρ taken from
SESAM's trajectory density at t = 0 of each run history and σ from `--sigma`; it also writes
`weber_numbers.csv` next to the plot. `plot_ohnesorge_number.py` does the same for the liquid-drop
Ohnesorge number Oh = μ_liq / √(ρ_liq σ d) (`--mu`, `--rho-liquid`, `--sigma`), which depends on the
diameter only, and writes `ohnesorge_numbers.csv` pairing every point's We with its Oh. `plot_knudsen_number.py` maps
Kn = λ / d at the initial state, banded by flow regime (continuum / slip / transitional / free
molecular at Kn = 0.01 / 0.1 / 10); λ = C / ρ_air with C calibrated on SESAM's own Knudsen column
(reference length = diameter), which the history files only print to 5 decimals. Writes `knudsen_numbers.csv`. `summary_from_runs.py`
reconstructs any material's summary from the per-run JSONs (same file name the sweep writes),
for summaries produced by older versions of the sweep or otherwise lost. `analysis/plot_style.py`
holds the shared chart palette (dataviz skill's validated sequential blue ramp).

## Physics model — `reentry_model` (Step 1: trajectory)

A first-principles re-entry model of a solid sphere, built to be verified against SESAM. Step 1 integrates the
trajectory only (constant mass): rotating-Earth 3-DOF in ECEF coordinates with J2 gravity, NRLMSISE-00 (pymsis)
with the fap-file solar activity or the US76 table SESAM ships, SESAM's sphere drag tables blended by a Knudsen
bridging function measured from SESAM's own output, DOP853 integration. Design:
`docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md`; the SESAM facts it relies on:
`Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md`.

```bash
"$PY" -m pip install "pymsis==0.13.0" scipy     # once, in drama_env
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
    --reference data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind.csv
"$PY" -m reentry_model run --diameter 50 --velocity 7.5 --altitude 115 --flight-path-angle -0.959331 \
    --atmosphere replay:data/reference_runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind.csv \
    --reference data/reference_runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind.csv
"$PY" -m reentry_model compare --model reentry_model_output/<run>.csv --reference data/reference_runs/<sesam run>.csv
```

Heading, latitude, longitude and epoch default to the reference cases' break-off state. `--atmosphere` is
`nrlmsise` (default), `us76`, or `replay:<sesam.csv>` (SESAM's own density/temperature, to isolate the dynamics);
`--wind none|static`; `--bridging sesam-table|sesam-erf|sin2|textbook` (matting is reserved and not yet
implemented — selecting it exits 1); `--gravity point|j2|j2j4`. Outputs go to
`reentry_model_output/` (git-ignored): `<run>.csv` with the same columns as the SESAM histories, `<run>.json`
(inputs, settings, results, comparison metrics, provenance), and with `--reference` a folder of six plots
(V(t) and h(t) overlays with residuals, h(V), angles, ground track, Knudsen/C_D). `data/reference_runs/` holds the
four committed SESAM references (100 mm from 77.5 km, 50 mm from 115 km, winds on/off).

## Tests

```bash
"$PY" -m pytest -m "not drama"     # unit tests, no DRAMA needed (real SESAM outputs in tests/fixtures/)
"$PY" -m pytest                    # also the integration tests that run SESAM
```

To refresh a fixture see `tests/fixtures/README.md`.
