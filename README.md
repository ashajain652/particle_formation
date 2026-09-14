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
Optional state values default to 0°/0°/0° and the parent epoch. `--dry-run` prints the SESAM
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
`sweep_manifest.json`, `sweep_summary.csv` (one row per matrix point), `sweep.log`.

The velocity grid's lower bound is the parent's minimum velocity, read from the parent file and
logged as `v_min = ... km/s at t = ... s, altitude ... km (parent row N)`.

## Tests

```bash
"$PY" -m pytest -m "not drama"     # unit tests, no DRAMA needed (real SESAM outputs in tests/fixtures/)
"$PY" -m pytest                    # also the integration tests that run SESAM
```

To refresh a fixture see `tests/fixtures/README.md`.
