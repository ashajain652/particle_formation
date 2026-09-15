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
"$PY" analysis/plot_outcome_by_diameter.py                             # one panel per sphere size (temperature x velocity)
"$PY" analysis/plot_outcome_by_diameter.py --metric mass_loss
```

```bash
"$PY" analysis/plot_max_temperature_by_diameter.py                                  # AA7075 by default
"$PY" analysis/plot_outcome_by_diameter.py --summary sphere_sweep_output/sweep_summary_user-moltenAA7075.csv
"$PY" analysis/summary_from_runs.py --material user-moltenAA7075   # rebuild a summary from runs/*.json if ever lost
```

The plot scripts read `sphere_sweep_output/sweep_summary_AA7075.csv` by default and write PNGs
to `sphere_sweep_output/plots/` (`--summary`/`--out` override either path — pass another
material's `sweep_summary_<material>.csv` to plot it). `--metric outcome` (default) is the
binary `outcome == "demised"` result per run; `--metric mass_loss` uses the continuous
`mass_loss_fraction` column instead — see the module docstrings for what each represents.
`plot_max_temperature_by_diameter.py` shows, per sphere size, the maximum temperature reached
by spheres that lost no mass, with every mass-losing cell black. `summary_from_runs.py`
reconstructs any material's summary from the per-run JSONs (same file name the sweep writes),
for summaries produced by older versions of the sweep or otherwise lost. `analysis/plot_style.py`
holds the shared chart palette (dataviz skill's validated sequential blue ramp).

## Tests

```bash
"$PY" -m pytest -m "not drama"     # unit tests, no DRAMA needed (real SESAM outputs in tests/fixtures/)
"$PY" -m pytest                    # also the integration tests that run SESAM
```

To refresh a fixture see `tests/fixtures/README.md`.
