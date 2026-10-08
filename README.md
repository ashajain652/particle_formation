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
trajectory only (constant mass; Step 2, coupled 3D heat transfer, follows below): rotating-Earth 3-DOF in ECEF coordinates with J2 gravity, NRLMSISE-00 (pymsis)
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
implemented — selecting it exits 1); `--gravity point|j2|j2j4`. Both SESAM bridgings were fitted to
SESAM's own output, which includes these four reference runs, so the C_D agreement in the verification
is in-sample by design; the trajectory comparison is the independent check. Outputs go to
`reentry_model_output/` (git-ignored): `<run>.csv` with the same columns as the SESAM histories, `<run>.json`
(inputs, settings, results, comparison metrics, provenance), and with `--reference` a folder of six plots
(V(t) and h(t) overlays with residuals, h(V), angles, ground track, Knudsen/C_D). `data/reference_runs/` holds the
four committed NRLMSISE-00 SESAM references (100 mm from 77.5 km, 50 mm from 115 km, winds on/off, `_msis`) and the
same two spheres run on SESAM's static US76 table (winds off, no `_msis` in the name). The US76 pairs give the
like-for-like check — SESAM and the model on the identical atmosphere table (`--atmosphere us76`) — and agree over
the hypersonic phase to 3.9 m/s (0.24 %) / 18 m for the 100 mm sphere and 3.6 m/s (0.06 %) / 4 m for the 50 mm sphere,
impact times within 0.2 s and 0.8 s:

```bash
"$PY" -m reentry_model run --diameter 50 --velocity 7.5 --altitude 115 --flight-path-angle -0.959331 --atmosphere us76 \
    --reference data/reference_runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_nowind.csv
```

### Verification (Task 10, `analysis/reentry_model_verification.py`)

Model sampled at SESAM's own time stamps; errors over the hypersonic phase (V > 1 km/s). "replay" feeds SESAM's
density/temperature into the model (dynamics and drag only); "nrlmsise" is the full model. The winds-on
references are compared with the wind-free model, so their rows include SESAM's HWM14 wind effect (≤ 8 m/s).

| case | mode | winds in ref | hypersonic max dV | hypersonic max dh | whole-flight max dV / dh | end time | runtime |
|---|---|---|---|---|---|---|---|
| d100.00mm_h077.500km_nowind | replay | off | 7.3 m/s (0.315%) | 22 m | 7.3 m/s / 23 m | -1.1 s (-0.30%) | 2 s, 24038 evals |
| d100.00mm_h077.500km_nowind | nrlmsise | off | 201.3 m/s (6.802%) | 115 m | 201.3 m/s / 148 m | -0.4 s (-0.11%) | 61 s, 634160 evals |
| d100.00mm_h077.500km | replay | on | 13.4 m/s (0.555%) | 27 m | 13.4 m/s / 28 m | -0.4 s (-0.10%) | 448 s, 6714263 evals |
| d100.00mm_h077.500km | nrlmsise | on | 197.0 m/s (6.553%) | 114 m | 197.0 m/s / 147 m | +0.3 s (+0.08%) | 64 s, 634160 evals |
| d050.00mm_h115.000km_nowind | replay | off | 10.7 m/s (0.428%) | 12 m | 10.7 m/s / 17 m | -0.3 s (-0.06%) | 182 s, 2703485 evals |
| d050.00mm_h115.000km_nowind | nrlmsise | off | 191.9 m/s (6.148%) | 430 m | 191.9 m/s / 430 m | +1.6 s (+0.29%) | 38 s, 386111 evals |
| d050.00mm_h115.000km | replay | on | 16.2 m/s (0.927%) | 22 m | 16.2 m/s / 28 m | -0.8 s (-0.14%) | 2 s, 34247 evals |
| d050.00mm_h115.000km | nrlmsise | on | 190.6 m/s (5.875%) | 420 m | 190.6 m/s / 420 m | +1.1 s (+0.19%) | 39 s, 386111 evals |

Notes:
- Replay mode: the residual builds up in the continuum peak-deceleration phase (Kn < 0.01, V 3-5 km/s), where
  the model's C_D already equals SESAM's printed column to 1e-4; in the bridging region it is <= 1.8 m/s.
  SESAM's own printed `dynamic_pressure_Pa` exceeds 1/2 x density_kgm3 x velocity^2 computed from its own
  printed columns by a factor that grows from 1.000 at 77 km to 1.005 at ~4 km/s (40 km) and 1.010 at
  1.2 km/s, while `load_factor_g` equals `dynamic_pressure_Pa` x C_D x A / (m x g0) to 1e-4 — so the density
  SESAM's drag actually used is 0.1-1 % higher than the density it prints (consistent with a ~0.05 s lag
  between the printed density and the integrated state). A model fed the printed density decelerates
  slightly less: model faster by <= 10.7 m/s winds-off, <= 16.2 m/s winds-on (HWM14 adds ~6 m/s) — the
  observed residual in the observed phase. The replay comparison is limited by the reference file's
  precision, not by the model; `dh_max_m` <= 27 m and `d_end_time_rel` <= 0.30 % confirm the dynamics.
  Threshold set to 1 %. A replay that reconstructs density as `2 x dynamic_pressure / V^2` from the
  reference's own columns would remove this artefact (not implemented).
- NRLMSISE-00 mode: the model's density (pymsis, the NRL reference implementation) differs from the density SESAM's built-in NRLMSISE-00 produced for the same epoch, place and solar inputs by an altitude-structured ratio: pymsis/SESAM = 0.70 at 110 km, 0.63 at 105 km, 0.63 at 100 km, 0.71 at 95 km, 0.85 at 90 km, 0.99 at 85 km, 1.06 at 80 km, 1.05 at 70 km, 1.07 at 60 km, 1.09 at 50 km, 1.07 at 40 km, 1.02 at 20 km (1 Aug 2024 12:53 UT, 29.5° N 82.1° W, F10.7 234 / 194, Ap 19). The solar inputs are not the cause (F10.7 170–246, Ap 8–56, storm mode and pymsis's historical indices move the ratio by < 2 %). A SESAM epoch experiment (same state, 1 Feb vs 1 Aug 2024) shows SESAM's density does respond to the date — Feb/Aug = 1.18–1.26 at 95–110 km, 1.04 at 70–80 km, 0.97–0.98 at 40–50 km — but with roughly half the seasonal amplitude of the reference implementation (pymsis: 1.29–1.71 and 0.91–0.92), and it responds to local time (00 UT/12 UT = 0.93 at 110 km). The gap is therefore in SESAM's NRLMSISE-00 configuration (model variant or variation switches; the SARA modelling report would say which), and it is unresolved on the SESAM side. The `dV_rel_max` threshold for this mode (10 %) is a regression guard set at the measured residual × 1.5, not a verification of the atmosphere; the dynamics are verified by the replay mode.

## Physics model — `reentry_model` (Step 2: coupled 3D heat transfer)

Step 2 solves the trajectory and the temperature field inside the sphere together: every 0.5 s macro step the
trajectory advances (DOP853, Step 1 tolerances), an aerothermal model turns the freestream state into a convective
flux on each of the ~18 000 surface patches (18 078) (100 mm sphere, 2 mm surface elements), and a finite-element conduction step (P1 tetrahedra, backward Euler,
Newton on the ε σ T⁴ radiation term, energy-exact secant heat capacity) advances the field. Design:
`docs/superpowers/specs/2026-09-18-thermal-fem-design.md`; plan: `docs/superpowers/plans/2026-09-18-thermal-fem.md`.

```bash
"$PY" -m pip install -r requirements-step2.txt          # once, in drama_env (scikit-fem, gmsh, pyamg, pyvista, imageio-ffmpeg, cantera)
# verification mode: SESAM's own heat input, uniform over the surface (NOT physical, see below), vs the US76 SESAM reference
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --heating sesam --reference data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind.csv
# physics mode: Fay-Riddell (Cantera equilibrium air) + Matting bridging + Lees distribution, with the surface-temperature animation
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --heating physics --animate
"$PY" analysis/reentry_model_thermal_verification.py     # both spheres x both modes -> reentry_model_output/verification_thermal/summary.md
```

Options (`--thermal fem`): `--heating physics|sesam`; `--stagnation fay-riddell|sutton-graves|dkr`; `--bridging-heat
matting|sesam-table`; `--matting-n` (1); `--accommodation` (0.8); `--catalycity` (1); `--material` (packaged
`AA7075_nomelt`), `--emissivity` (material's 0.40), `--t-ambient` (0 K, SESAM's; 200 K optional); `--h-surface 2`,
`--h-core 8` (mm; 18.9 k nodes on the 100 mm sphere), `--mesh-size` (multiplier); `--thermal-solver skfem|fenicsx`,
`--linear-solver amg|direct`, `--lumped-mass`; `--dt 0.5`; `--frames-every`, `--animate`, `--stills`. Outputs:
`<run>.csv` gains `convective_heat_W, rad_cooling_W, integrated_heat_J, absorbed_heat_J, surface_T_max/min/mean_K,
T_stagnation_K, T_back_K, T_centre_K, q_stag_Wm2, heating_blend_f` (`temperature_K` is the energy-equivalent mean
temperature, the quantity SESAM's lumped model reports; `heating_blend_f`: SESAM-equivalent mode — the measured heat
factor F_h(Kn); physics mode — the free-molecular weight w of the distribution (1 − q_stag/q_c under Matting, f(Kn)
under the SESAM table; w = 1 when the hot-wall clamp zeroes q_c)); `<run>/vtk/` holds `field.pvd` + `field_<k>.vtu` (nodal T) and
`surface.pvd` + `surface_<k>.vtp` (per-patch q_conv, q_rad, T), `<run>/vtk/animation.mp4` (the surface temperature) and
`<run>/vtk/section.mp4` (the meridional cross-section through the flight axis, windward side on the right; GIF fallback
for both), `frames/`, `frames_section/`, `stills/` (surface and `section_*` stills at the start, peak heating, peak surface
temperature and the end); with `--reference`, three more plots (`heating_time`, `temperature_time`, `integrated_heat`) and
`comparison.thermal_metrics` in the JSON. A 100 mm flight takes ~2.5 min (SESAM-equivalent) / ~3 min (182 s with the
animation). The FEniCSx backend (`--thermal-solver fenicsx`) runs from the separate conda environment `fenicsx_env`
(conda-forge `fenics-dolfinx` 0.11 + this repo's pip dependencies); selecting it in `drama_env` exits 2. It is verified:
the seven conformance tests (the analytic cases and a 0.1 % cross-check against the skfem backend) pass, and a 100 s
coupled flight gives the same temperatures, heat totals and energy balance as the skfem backend to every printed
digit. In serial it is ~3x slower (62 s vs 19 s for those 100 s: dolfinx assembles the forms by quadrature every
Newton iterate, where the skfem backend rescales precomputed element matrices); its purpose is MPI scaling for the
larger spheres, which is not exercised yet. Two machine notes for `fenicsx_env` on this Mac: run with
`FI_PROVIDER=tcp` (MPICH's libfabric otherwise aborts at interpreter exit with "OFI poll failed"), and with
`CC=$CONDA_PREFIX/bin/clang` if the environment is not activated (FFCx JIT-compiles the forms; Apple's linker
cannot read the macOS 27.0 SDK on this machine while conda's clang can):

```bash
FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m reentry_model run --diameter 100 --velocity 7.5 \
    --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --thermal-solver fenicsx
```

**`--heating sesam` is a verification device, not a physical model.** It applies SESAM's surface-averaged heat input —
0.27471 × q_DKR × F_h(Kn) × hot-wall factor — uniformly to every patch, front and back, so that the conduction,
time stepping, material curves and coupling can be compared with SESAM's lumped temperature and heat totals with no
distribution question in between. Real heating is concentrated on the windward face (`--heating physics`: Lees'
laminar distribution integrates to 0.196 of the stagnation flux over the sphere, zero leeward).

### Verification (`analysis/reentry_model_thermal_verification.py`, `tests/test_reentry_model_reference_thermal.py`)

Coupled model vs the no-melt US76 SESAM references (winds off), default mesh, Δt 0.5 s. Errors over the hypersonic
phase (V > 1 km/s); power errors are relative to SESAM's peak (its hot-wall factor clamps the heating to zero late in
the flight, where relative errors are unbounded); the point-wise error is over the continuum part (Kn < 0.01, Q > 10 %
of peak); `T_eq` is the energy-equivalent mean temperature vs SESAM's lumped temperature over the whole flight.

| case | heating | Q_conv max (of peak) | Q_conv point-wise (continuum) | integrated heat (end of hypersonic / end) | max ΔT_eq | radiated (of peak) | peak surface T | runtime |
|---|---|---|---|---|---|---|---|---|
| d100.00mm_h077.500km | sesam | 0.46 % | 2.24 % | +0.31 % / +0.30 % | 24.5 K (1.20 %) | 6.65 % | 2089 K at 139 s (SESAM lumped peak 2104 K at 143 s) | 158 s, 732 steps |
| d050.00mm_h115.000km | sesam | 0.37 % | 2.36 % | +0.03 % / +0.02 % | 28.0 K (1.22 %) | 3.40 % | 2442 K at 245 s (SESAM lumped peak 2458 K at 245 s) | 47 s, 1117 steps |
| d100.00mm_h077.500km | physics | — | — | ratio to SESAM 0.744 | — | — | 2319 K stagnation at 118 s, mean 1717 K at 150 s | 182 s incl. animation, 732 steps |
| d050.00mm_h115.000km | physics | — | — | ratio to SESAM 0.770 | — | — | 2493 K stagnation at 239 s, mean 2151 K at 250 s | 54 s, 1117 steps |

Both physics-mode runs now reach the ground (`end_reason == "ground"`, same 732 / 1117 steps as their sesam-mode
counterparts — the trajectory is heating-mode-independent). The Q_conv/point-wise/ΔT_eq/radiated columns are "—"
for physics mode because those thresholds are SESAM-equivalent-mode-specific (spec section 10): physics mode's
raw heat distribution and timing are deliberately different from SESAM's surface average (windward-concentrated
vs. uniform; for a sphere the ATDB factor is the same for any attitude, see docs/model_assumptions.md §8), so a point-wise or peak-power comparison against SESAM isn't meaningful there — only the
integrated-heat ratio is reported, per spec section 10.

Thresholds (`tests/test_reentry_model_reference_thermal.py`): Q_conv 3 % of peak and 3 % point-wise, integrated heat
3 %, T_eq 2 %, radiated power 8 % (= 4 × the temperature margin: the resolved surface radiates at its own, hotter
temperature — measured 6.65 % while T_eq was within 1.20 %). Refinement (100 mm, to 200 s): halving `h_surface`
(1 mm / 8 mm, 76 k nodes) changes the surface-temperature history by ≤ 0.03 % in SESAM-equivalent mode and ≤ 0.11 %
in physics mode (stagnation temperature ≤ 0.11 %); halving Δt changes the peak surface temperature by 0.04 %.
Energy balance closes to 1e-7 or better over every flight.

Note that F_h(Kn) and the hot-wall c_p were measured on these same two references, so the Q_conv and integrated-heat
rows measure the fit's residual (2.6 % / 1.4 % point-wise when the formula is applied to SESAM's own state columns);
the independent content of the verification is the temperature (T_eq within 1.2 %) and the radiated power.

Findings recorded while building this step (details in `sesam_verified_facts.md` §§15–16 and the spec's amendments):
- SESAM's convective heating carries a hot-wall factor max(0, 1 − c_p(T − T∞)/(V²/2)) with c_p ≈ 1004.5 J/kg/K
  (invisible below 850 K, where the earlier facts were measured; it makes SESAM's heating vanish below ~1.8 km/s once
  the no-melt sphere is at 2000 K), and it halves the continuum heating below Mach 1.
- SESAM's transitional heating is not its drag blend: relative to 0.27471 × q_DKR it is 0.94 at Kn 0.04, 0.46 at
  Kn 0.2, 0.14 at Kn 1 and 0.059 at Kn 40 — the last being 0.78 × the textbook free-molecular cos θ average, so SESAM's
  free-molecular limit is the ordinary ½ρV³ with α ≈ 0.8 (the "13× too low q_FM" of the Step 1 facts was the
  transitional deficit misread at Kn 0.03). The measured factor F_h(Kn) (`aero.SesamHeatTable`) is what the
  verification mode uses; a first attempt with the textbook blend gave +12 % / +17 % integrated heat.
- Fay–Riddell with Cantera's equilibrium air is 1.29 × Sutton–Graves at the 100 mm start (72 % of the stagnation
  enthalpy is dissociation at 1.5 kPa; Sutton–Graves is a Le = 1 fit), 1.14 × with the Lewis-number term off
  (independently re-measured this task: q_FR = 2.22e6 W/m², q_SG = 1.72e6 W/m²); `airNASA9.yaml` has no
  transport data, so viscosity comes from Blottner fits with Wilke mixing.
- Physics mode delivers 0.74–0.77 × SESAM's integrated heat (0.196/0.2747 × Fay–Riddell/DKR × hot wall) with a
  600 K stagnation-to-mean temperature difference at peak heating — the reason Step 3 needs the resolved field.
- Both thermal backends raise `RuntimeError` when their Newton iteration fails to converge; the CLI catches this
  (and other model errors) and exits 1 with the message, rather than hanging or producing a silently wrong result.
- The transonic tail is compressed isentropically without a shock because the shock fixed point stalls at Ma 1 —
  found by the first full physics-mode run (facts §16).

## Large fragments — `spheral_frag` core (M1)

The Spheral large-fragment model (spec `docs/superpowers/specs/2026-10-02-spheral-large-fragments-design.md`) replays
the finite-element flight in Spheral's particle model to find what tears, drips or detaches as whole fragments. M1 is
its core without Spheral (plan `docs/superpowers/plans/2026-10-07-spheral-m1.md`):
- frame import through one contract table;
- material tables shared with the finite element;
- thickness, layer depths and the three zones;
- load tables by inclination;
- the fragment record, the debris log and the mass accounts;
- the `prepare` and `analyse` entry points.

Everything runs in drama_env with numpy; nothing here imports Spheral. Assumptions, each marked verified, analytic,
measured or assumed, are in `docs/spheral_frag_assumptions.md`.

`reentry_model` does not melt yet, so the frames come from the reconstructed Step 3 prototype,
`prototype/work-2026-10-08-spheral-mvp/code` (an MVP input, plan decision 13). The flight below is its 100 mm
`AA7075_scheil` physics run: 0.5 s step, seed 12345, 1,255 frames, 5.0 GB, git-ignored.

```bash
FE=reentry_model_output/model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off
PKG=prototype/work-2026-10-08-spheral-mvp/code

# frames -> compact prepared frames, material table, load tables, flight table and prepare.json (26 min, 791 MB)
"$PY" -m spheral_frag prepare --fe-run "$FE" --fe-package "$PKG" [--frames K0:K1] [--every N] [--force]
# a Spheral run's checks -> fragments.csv, debris.csv, analyse.json (accounts must close to 1e-12 m0, else exit 1)
"$PY" -m spheral_frag analyse --run spheral_output/runs/<run> [--prepared <dir>] [--min-particles 30]
# per-frame measurements -> data/spheral/m1_frames.csv (committed) and plots + summary.json in spheral_output/m1/
"$PY" analysis/spheral_m1_flight.py [--fake-runs]
# the thresholds on the real flight (marker fe_flight; skipped unless both variables are set)
SPHERAL_FRAG_FE_RUN="$FE" SPHERAL_FRAG_FE_PACKAGE="$PKG" "$PY" -m pytest tests/test_spheral_frag_flight.py -q
```

**Output layout.** `prepare` writes `spheral_output/prepare/prep_<fe run>_k<K0>-<K1>/`, containing:

| File | Contents |
|---|---|
| `frames/mesh_<i>.npz` | the finite-element nodes and tetrahedra of a mesh version |
| `frames/frame_<k>.npz` (+ `.json` sidecar) | per frame, stored reduced: active-tetrahedron mask, φ ≠ 1, moved nodes, the outward faces, T, f_l, the patch fields, thickness and layer depths; nodes, tetrahedra and patch geometry rebuilt bitwise on loading (`frames.load_prepared_frame`, numpy only) |
| `material_table.npz` | the material table |
| `loads.npz` | the windward load table per frame |
| `flight.npz` | the history's contract columns plus air temperature, deceleration and drag |
| `prepare.json` | provenance (the finite-element run's and package's SHA-256), the contract as found, a per-frame entry, gating checks with a priori thresholds and the measured ones |

`analyse` writes `spheral_output/analyse/<run>/`. Names encode the whole configuration (`spheral_frag/naming.py`).
Exit codes are as for `reentry_model`: 0 ok, 1 a contract, check or accounts failure, 2 bad arguments or missing
inputs.

### Verification (`tests/test_spheral_frag_*.py`, `analysis/spheral_m1_flight.py`)

Synthetic devices (`tests/spheral_frag_synthetic.py`, fixtures in `tests/fixtures/spheral_frag/`):
- a coarse 100 mm sphere whose third frame has ten nose elements dead;
- a dumbbell with a 4 mm neck;
- a slab with five columns of known slurry depth.

The real flight is the one above, prepared whole on 2026-10-08.

| Check | Reference | Threshold (source) | Measured |
|---|---|---|---|
| Contract | every `FE_FIELDS` item on the flight | all present and confirmed (Task 10) | 47 of 47; `p_w` = 0 on frame 0 only |
| Frame import | pyvista's own arrays | bitwise (a priori) | bitwise |
| Open directed edges of the outward-wound surface | closed body | 0 (a priori) | 0 on all 1,255 frames; up to 1,662 faces per frame wound inward in the export, all reoriented |
| Non-manifold edges / vertices | — | 300 / 20 (flight worst, rounded up) | 288 (1,204 frames) / 15 (131 frames) |
| Divergence volume vs Σ tetrahedra | round-off | 1e-12 (a priori) | 2.2e-16 |
| Mass Σ φρV + film + deep vs history `mass_kg` | the CSV's 9 digits | 5e-9 (a priori) | 4.2e-9 |
| Film mass vs `film_mass_kg` | the CSV's 9 digits | 5e-9 | 4.1e-9 |
| Material h, f_l, c_p vs the finite element (100,543 temperatures) | the finite-element material | 1e-12 max\|h\| + 1e-9 J/kg; f_l bitwise (a priori) | h, f_l, c_p bitwise; inverse 1.6e-12 K |
| Nodal `liquid_fraction` vs table f_l(T) | the finite element | bitwise | bitwise on every frame |
| `delta_m` finite exactly where the closure is Girin's | contract | 0 mismatches | 0 |
| Thickness: sphere / dumbbell neck / slab | 2R / 2 r_neck / box extent | faceting bound 0.907 % / within the neck's facets / 1e-12 m | 0.816 % / 7.879–7.992 mm / exact |
| Thickness on the flight's frame 0 | 2R | 6e-4 (flight worst) | 5.8e-4; no escaping ray on any frame; thinnest 0.21 mm (k 166) |
| Slurry depth on the slab | exact columns | 1e-12 m | 5e-18 m |
| Zones on the slab at 2 / 3 mm film limits | the rule by hand | exact | zones 1, 2, 3, 3, 2 / 1, 2, 2, 3, 2 |
| Load table vs Newtonian sphere | p_stag cos²θ | bin averaging; faceting 1 % (plan Task 7) | table within 0.061 %; faceting −0.55 % |
| Binning \|D_table / D_patch − 1\| on the flight | the frame's own loads | 8e-4 (flight worst) | 7.85e-4 |
| Largest `p_w` vs `p_w_stag_step_Pa` | the history | 0.2 (flight worst) | median 6e-10, worst 0.129 where the stagnation patch died (k 64) |
| Σ release_rate · A vs the history's sprayed increment | the history | per step −0.3 to +0.04; flight total 0.2 | −0.295 to +0.038; total −0.123 (1.025 of 1.169 kg; release on removed faces not carried) |
| **Drag from the frame's loads vs the history's drag** | spec §7.3: "a few percent" (5 %) | 0.4 (flight worst; **5 % not met**) | +10 % before the first film, −36 % to +9 % shock layer, −11.0 % free-molecular and subsonic branch |
| Record: lattice ellipsoid 10/6/4 mm at 0.5 mm | exact ellipsoid | principal lengths, equivalent diameter | within 0.09 mm, 0.10 % |
| Record: 10 mm liquid drop at 5 km/s | closed form | We, Oh, threshold | 3,125, 2.967e-4, 12.000 |
| Accounts (synthetic history; flight frames 100–110 and 800–810, 3D, dx 2.2 mm) | starting mass | 1e-14 m₀ (a priori) | 2e-18; 7e-18 and 4e-17 m₀ |

**The drag gap** is M1's answer for the thesis, reported rather than fitted (plan Review focus 6). The history's drag
comes from the trajectory's SESAM-table drag coefficient, and the frame's loads from Step 3's surface-flow model. The
two already differ by 10 % on the intact sphere, before anything melts. It is a minor issue (Asha, 2026-10-08): a
hypersonic CFD code is to supply the loads Spheral receives, retiring both models' drag; until then Spheral takes the
frames' loads.

Measured on the flight:
- **Zones:** the bulk zone exists on 336 frames (24.5–192 s, at most 152 cm² at the 2 mm film limit).
- **Slurry:** it reaches the far side of the body on 250 frames (25.5–150 s).
- **Cost:** a median of 0.45 s per frame, of which the thickness map takes 0.26 s.

These are 0.5 s-step frames, so part of the deep layers is spec §2's molten backlog.

## Tests

```bash
"$PY" -m pytest -m "not drama and not reference" -q   # unit tests (~4 min; the thermal solver and coupled tests dominate)
"$PY" -m pytest -m "not drama"     # unit tests, no DRAMA needed (real SESAM outputs in tests/fixtures/)
"$PY" -m pytest                    # also the integration tests that run SESAM
"$PY" -m pytest -m reference -q   # the Step 1 reference flights (~15 min) and the Step 2 coupled runs (~35 min)
```

`tests/test_reentry_model_fenicsx.py` runs only with an interpreter that can import `dolfinx` (the `fenicsx_env`
environment, see above; the whole unit tier also passes there apart from the sweep wrapper's `tqdm` dependency);
elsewhere it is skipped. To refresh a fixture see `tests/fixtures/README.md`.

The `spheral_frag` tests marked `fe_flight` run only when `SPHERAL_FRAG_FE_RUN` and `SPHERAL_FRAG_FE_PACKAGE` name the
M1 flight and its package (see "Large fragments" above); those marked `spheral` only when the Spheral container
launcher can import Spheral.
