# Sphere fragment re-entry sweep with ESA DRAMA / SARA (SESAM) — Design

Date: 2026-09-13
Status: approved design, awaiting implementation plan

## 1. Purpose

Model the atmospheric descent of solid AA7075 spheres that break off a
re-entering satellite. Two scripts:

1. `sphere_reentry.py` — runs SESAM (SARA's re-entry module, via the pyDRAMA
   package) for one sphere from a given initial state and writes a
   full-resolution CSV history plus a JSON file of statistics.
2. `sphere_sweep.py` — drives script 1 over a grid of sphere diameters,
   initial temperatures and initial velocities. The remaining initial-state
   values (altitude, latitude, longitude, flight-path angle, heading, epoch)
   are inherited from a parent-satellite SESAM run at the point where the
   parent reached the sphere's velocity.

Script 1 must be usable on its own; script 2 calls it as a subprocess.

## 2. Verified facts (spike, 2026-09-13)

Implementers should rely on these; they were measured with DRAMA 4.1.4 /
SESAM 2.3.0 through pyDRAMA in `drama_env`.

| Fact | Detail |
|---|---|
| Geodetic initial state | `coordinateSystem: "geodetic"` with `element1..6` = altitude [km], latitude [deg], longitude [deg], velocity [km/s] (relative to the rotating atmosphere), flight-path angle [deg, negative = descending], heading/azimuth [deg]. SESAM's first trajectory row reproduces all six inputs. Same conventions as the `Trajectory.txt` columns `altitude, lat, lon, velocity, path, heading`. |
| Low starts are fine | Starts at 77.5, 39.9, 12.9 and 1.6 km all ran to ground impact. The 100 km / 140 km floors exist only in pyDRAMA's JSON-schema path (`sara.run(json=...)`); passing `config=[cfg]` bypasses that validation and SESAM has no such floor. |
| `config` must be a list | `sara.run(config=[cfg])`. A bare dict goes through ConfigCreator, which flattens list-valued keys (`objects`) and breaks materials handling. |
| Run time | ~0.15 s per SESAM run including pyDRAMA overhead; Python start + `import drama` ~0.1 s. |
| Parallel safety | pyDRAMA runs each job in its own `tempfile.TemporaryDirectory()`; `save_output_dirs` must not already exist (`shutil.copytree`). |
| Output files | `<save_output_dirs>/run_0/reentry/`: `PySara.<object>.<uuid>_AeroThermalHistory.txt` (11 columns), `PySara.<object>.<uuid>_Trajectory.txt` (16 columns), `PySara.ImpactingFragments.xml`, `sesam.log`, `*.gnu`. File prefix is always `PySara` regardless of `runID`. Both history files use the same time grid. |
| History-file precision | Mass is printed with 3 decimals in kg (a 5 mm sphere, 1.84e-4 kg, prints `0.000`). `thick` has 4 decimals in mm. |
| Precise final mass | `ImpactingFragments.xml` (`<mass unit="kg">1.8411041946975187e-01</mass>`, `<velocity unit="m/s">`) for objects that reach the ground; `sesam.log` line `EVENT end <uuid> <mass>` (6 decimals) for every object. |
| `thick` column | Not the sphere radius: it drops 25.0 → 18.8 mm at the first melt step while mass drops 1.7 %. Record it, do not derive geometry from it. |
| End of life | `sesam.log` line `Fragment <uuid> reached end of life because of ground impact.` / `... because uncritical.` / `Object <uuid> reached end of life because of ballooning.` A fully melted sphere ends by ballooning with ~3 mg residual mass. |
| `energyThreshold` | SESAM stops tracking a fragment once its kinetic energy is below this value ("uncritical"). With the parent's 15 J a 5 mm sphere is dropped at t = 0 below 404 m/s. With `1e-9` J the same sphere is propagated to the ground. It is a risk-bookkeeping cutoff, not physics. |
| Material | `drama-AA7075` exists in the default `materials.xml`: density 2813 kg/m3, melting temperature 850 K, latent heat 400 kJ/kg, emissivity 0.40, oxidation disabled. No `materialList` injection needed. |
| pyDRAMA quirk | `_sara_xml` writes `solarActivity/valuesFromFile` from `useEnvironmentCSV`, not from `solarActivityFromFile`. The parent GUI run went through the same code, so mirroring the parent's keys keeps the two runs consistent. Do not try to fix this. |
| Environment | pyDRAMA lives in conda env `drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`, Python 3.12) with numpy 2.5, matplotlib, jsonschema, dict2xml, python-dateutil. `tqdm` and `pytest` must be pip-installed. `DRAMA_INSTALL_PATH=/Applications/DRAMA-4.1.4`. |

## 3. Architecture and layout

Everything lives in `Particle Wake Evolution/` (this repository).

```
sphere_reentry.py                 script 1 (CLI + importable module)
sphere_sweep.py                   script 2
data/fap_day.dat, data/fap_mon.dat  space-weather files copied once from
                                  "Generic_Satellite Reentry/data/" (DRAMA rewrites that folder)
tests/conftest.py, tests/test_sphere_reentry.py, tests/test_sphere_sweep.py
tests/fixtures/                   small real SESAM outputs + a miniature dmf_output.json
docs/superpowers/specs/, docs/superpowers/plans/
sphere_sweep_output/              git-ignored
  runs/<run_name>.csv, runs/<run_name>.json
  raw/<run_name>/                 raw DRAMA tree, kept only for failed runs or --keep-raw
  sweep_manifest.json, sweep_summary.csv, sweep.log
```

`Generic_Satellite Reentry/` (DRAMA GUI project) is git-ignored; the sweep
records the SHA-256 of the parent `dmf_output.json` it used.

Data flow:

```
dmf_output.json (parent) --> sphere_sweep.py --(subprocess per point)--> sphere_reentry.py --> pyDRAMA/SESAM
                                   |                                          |
                                   v                                          v
                     sweep_manifest.json, sweep_summary.csv          runs/<run_name>.csv + .json
```

## 4. Script 1 — `sphere_reentry.py`

### 4.1 CLI

```
python sphere_reentry.py --velocity V --altitude H --temperature T --diameter D
    [--flight-path-angle 0.0] [--heading 0.0] [--lat 0.0] [--lon 0.0]
    [--epoch 2024-08-01T12:00:00]
    [--outdir sphere_sweep_output/runs] [--raw-dir sphere_sweep_output/raw]
    [--timeout 600] [--keep-raw] [--dry-run] [--quiet]
    [--fap-day data/fap_day.dat] [--fap-mon data/fap_mon.dat]
```

| Argument | Unit | Notes |
|---|---|---|
| `--velocity` | km/s | required, > 0 |
| `--altitude` | km | required, >= 0 |
| `--temperature` | K | required, > 0 (initial bulk temperature) |
| `--diameter` | mm | required, > 0 |
| `--flight-path-angle` | deg | default 0.0; negative = descending |
| `--heading` | deg | default 0.0 |
| `--lat`, `--lon` | deg | default 0.0, 0.0 |
| `--epoch` | ISO-8601 UTC | default `2024-08-01T12:00:00` (parent run's epoch) |
| `--outdir` | path | CSV and JSON go here (created if missing); default relative to the script's directory |
| `--raw-dir` | path | root for raw DRAMA trees (`<raw-dir>/<run_name>/`); default `sphere_sweep_output/raw` relative to the script's directory |
| `--timeout` | s | passed to pyDRAMA; default 600 |
| `--keep-raw` | flag | keep the raw DRAMA tree on success |
| `--dry-run` | flag | print the config as JSON and the run name; create nothing |
| `--quiet` | flag | no per-run console output except errors (the sweep uses it) |
| `--fap-day`, `--fap-mon` | path | DRAMA space-weather files; default `data/` next to the script |

The four required inputs are the ones the sweep varies; the optional state
values default to fixed constants so the script runs standalone, and the sweep
always passes all of them explicitly.

Exit codes: 0 success; 1 run failed (JSON written with the error); 2 usage or
environment error (bad arguments, pyDRAMA not importable) with nothing written.

### 4.2 Module API (importable, pure where possible)

```python
RHO_AA7075 = 2813.0          # kg/m3, drama-AA7075 in the DRAMA material database
T_MELT_AA7075 = 850.0        # K
MELT_TOLERANCE_K = 0.5
ENERGY_THRESHOLD_J = 1e-9
CSV_MASS_RESOLUTION_KG = 1e-3

@dataclass
class SphereRun:              # all inputs of one run
    velocity_kms, altitude_km, temperature_K, diameter_mm,
    flight_path_deg=0.0, heading_deg=0.0, lat_deg=0.0, lon_deg=0.0,
    epoch: datetime = PARENT_EPOCH

def run_name(run: SphereRun) -> str
def sphere_mass_kg(diameter_mm) -> float           # rho * 4/3 pi r^3
def sphere_cross_section_m2(diameter_mm) -> float  # pi r^2
def build_config(run: SphereRun) -> dict           # the pyDRAMA config
def read_sara_table(path, columns) -> list[dict]   # '#'-comment whitespace tables
def parse_impacting_fragments(path) -> dict | None # {'mass_kg', 'velocity_kms', 'lat', 'lon', 'epoch'}
def parse_sesam_log(path) -> dict                  # {'end_of_life_reason', 'event_end_mass_kg', 'event_start_mass_kg'}
def merge_histories(aero_rows, traj_rows) -> (rows, warnings)
def compute_stats(run, merged_rows, fragments, log_info) -> (results: dict, warnings: list)
def write_csv(path, rows) -> None
def run_sphere(run, outdir, raw_dir, timeout, keep_raw, fap_day_lines, fap_mon_lines) -> dict  # the JSON document
```

Only `run_sphere` imports `drama`; everything else is testable without DRAMA.

### 4.3 Run name

```
sphere_d{diameter_mm:06.2f}mm_T{temperature_K:06.1f}K_v{velocity_kms:08.5f}kms_h{altitude_km:07.3f}km
```

Examples: `sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km`,
`sphere_d005.00mm_T0750.0K_v00.02800kms_h000.012km`. Fixed width so
lexicographic order equals parameter order. The other state values are in the
JSON; two runs that differ only in those would overwrite each other
(documented, acceptable — the sweep never does that).

### 4.4 SESAM configuration

Mirrors the parent satellite run (`Generic_Satellite Reentry/input/last_sara.json`)
except where marked.

```python
r = diameter_mm / 2000.0
obj = {
    "name": "sphere", "uniqueID": "5b3e0000-0000-4000-8000-000000000001",
    "primitive": {"sphere": {"radius": r}},
    "mass": RHO_AA7075 * 4/3 * pi * r**3,
    "material": "drama-AA7075",
    "solid": True,
    "relativePosition": {"cartX": 0, "cartY": 0, "cartZ": 0, "yaw": 0, "pitch": 0, "roll": 0},
    "scalingFactors": {"drag": 1, "lift": 1, "sideForce": 1, "averageHeatFlux": 1,
                       "averageHeatFluxCT": 1, "averageHeatFluxTR": 1, "averageHeatFluxFM": 1},
    "attitude": "tumbling", "quantity": 1, "temperature": T0,
}
cfg = {
    "runID": "SPHERE", "comment1": "<run_name>", "comment2": "solid AA7075 sphere fragment",
    "beginDate": epoch, "runMode": "reentry-only", "monteCarlo": False,
    "coordinateSystem": "geodetic", "initialDate": epoch,
    "element1": altitude_km, "element2": lat_deg, "element3": lon_deg,
    "element4": velocity_kms, "element5": flight_path_deg, "element6": heading_deg % 360,
    "assumedCrossSection": pi * r**2, "dragCoefficient": 2.2, "reflectivityCoefficient": 1.3,
    "attitude": "tumbling", "fragmentsAttitudeAfterBreakup": "inherited",
    "globalSpacecraftTemperature": T0,
    "densityScalingFactor": 1.0, "dynamicEnvironment": True, "useWind": True,
    "solarActivityFromFile": True, "useEnvironmentCSV": False, "ap": 8, "f107a": 170,
    "voxelatorMode": 1,
    "energyThreshold": 1e-9,          # DEPARTURE from parent (15 J), see section 2
    "plotVisibilityMaps": False, "plotObjectTrajectories": False,
    "propagationWithOscar": False,
    "objects": [obj],                  # no "materialList": built-in materials.xml is used
}
```

`beginDate`/`initialDate` are `datetime` objects (pyDRAMA calls `strftime`).

### 4.5 Execution

```python
raw_dir = <raw-dir>/<run_name>           # removed first if it exists
results = sara.run(config=[cfg], save_output_dirs=raw_dir, keep_output_files="all",
                   fap_day_content=fap_day_lines, fap_mon_content=fap_mon_lines,
                   parallel=False, timeout=timeout, log_level="ERROR", spell_check=False)
```

`os.environ.setdefault("DRAMA_INSTALL_PATH", "/Applications/DRAMA-4.1.4")` before
importing `drama`. A failed import exits with code 2 and the hint
`pip install /Applications/DRAMA-4.1.4/TOOLS/drama_python_package`.

After the run, locate `*_AeroThermalHistory.txt`, `*_Trajectory.txt`,
`*ImpactingFragments.xml`, `sesam.log` under `raw_dir` (recursive glob).

### 4.6 Parsing and merge

Columns (order as in the files):

- aero: `time, altitude, temp, mass, thick, convectiveHeat, radiativeHeat, oxidationHeat, radCooling, integratedHeat, visibilityFactor`
- traj: `time, altitude, lat, lon, velocity, downrange, drag, lift, side, knudsen, mach, path, heading, density, dynamicPressure, loadFactor`

Lines starting with `#` or with the wrong column count are ignored.

Merge on `time`: if the two time vectors are identical, zip; otherwise outer-join
on exact time values (sorted), leaving missing cells empty, and add the warning
`time grids differ (aero N rows, traj M rows)`.

### 4.7 CSV (`<run_name>.csv`)

Header row, one row per time step:

```
time_s, altitude_km, velocity_kms, temperature_K, mass_kg, thick_mm, lat_deg, lon_deg,
downrange_km, flight_path_deg, heading_deg, drag, lift, side, knudsen, mach,
density_kgm3, dynamic_pressure_Pa, load_factor_g, convective_heat_W, radiative_heat_W,
oxidation_heat_W, rad_cooling_W, integrated_heat_J, visibility_factor
```

Values are the parsed floats written with Python's shortest round-trip `repr`
(`0.184` stays `0.184`, `1.800e-05` becomes `1.8e-05`), so SESAM's printed
precision is preserved and nothing is re-rounded. The `mass_kg` column
inherits SESAM's 0.001 kg resolution (section 2).

### 4.8 Statistics

| Field | Definition |
|---|---|
| `max_temperature_K` | max of `temperature_K`; also `time_of_max_temperature_s`, `altitude_of_max_temperature_km` at its first occurrence |
| `time_at_melting_temperature_s` | sum over consecutive rows (i, i+1) that are both "at melt" of `t[i+1] - t[i]`, where at melt means `abs(T - 850.0) <= 0.5`. Also `n_rows_at_melt`, `first_time_at_melt_s`, `last_time_at_melt_s`, `altitude_first_melt_km`, `altitude_last_melt_km` (null when never at melt; a single row at melt gives duration 0). `melt_tolerance_K = 0.5`, `melting_temperature_K = 850.0` are recorded. |
| `final_mass_kg` | first available of: `ImpactingFragments.xml` mass → `sesam.log` `EVENT end` mass → last history row. `final_mass_source` ∈ {`impacting_fragments_xml`, `sesam_log_event_end`, `history_file`}. `initial_mass_kg` = computed mass. `mass_loss_fraction = 1 - final/initial`. |
| `final_radius_mm` | `1000 * (3 m / (4 pi rho))^(1/3)` from `final_mass_kg` (0 when the mass is 0). `final_thick_mm` = SESAM's last `thick` value, reference only. |
| `final_velocity_kms` | `ImpactingFragments.xml` velocity / 1000 when present, else last trajectory row; `final_velocity_source` recorded. |
| trajectory end | `final_time_s`, `final_altitude_km`, `final_latitude_deg`, `final_longitude_deg`, `downrange_km` from the last row |
| `end_of_life_reason` | from `sesam.log` regex `reached end of life because (?:of )?(.+?)\.` → `ground impact`, `uncritical`, `ballooning`, ...; `unknown` if absent |
| `sesam_version` | provenance, from the history-file header regex `SESAM ([\d.]+)` |
| `outcome` | `survived` if reason is `ground impact`; else `demised` if reason is `ballooning` or `final_mass_kg < 0.05 * initial_mass_kg` (a fully melted 5 mm sphere ends with ~3 mg = 1.6 % residual); otherwise `other` |
| `n_rows` | rows in the merged history |

Warnings (list of strings): coarse CSV mass (`initial_mass_kg < 0.05`), time-grid
mismatch, single-row history, missing `ImpactingFragments.xml` for a survivor,
unknown end-of-life reason.

### 4.9 JSON (`<run_name>.json`)

```json
{
  "schema_version": 1,
  "run_name": "sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km",
  "status": "ok",                       // "ok" | "error" | "timeout"
  "error": null,                        // message + last 40 log lines on failure
  "inputs": {
    "diameter_mm": 50.0, "radius_m": 0.025, "initial_mass_kg": 0.18411,
    "cross_section_m2": 0.0019635,
    "initial_temperature_K": 300.0, "initial_velocity_kms": 7.5, "initial_altitude_km": 77.5,
    "flight_path_angle_deg": -0.959, "heading_deg": 347.168,
    "latitude_deg": 29.546, "longitude_deg": -82.134, "epoch_utc": "2024-08-01T12:53:07",
    "material": "drama-AA7075", "material_density_kgm3": 2813.0,
    "melting_temperature_K": 850.0, "melt_tolerance_K": 0.5,
    "energy_threshold_J": 1e-9,
    "sesam_settings": { "...every non-object key of the config, JSON-serialised..." }
  },
  "results": { "...section 4.8 fields..." },
  "warnings": [],
  "files": { "csv": "runs/sphere_....csv", "raw_dir": null },
  "provenance": {
    "drama_install_path": "/Applications/DRAMA-4.1.4", "sesam_version": "2.3.0",
    "pydrama_version": "...", "script_version": "1.0.0",
    "wall_time_s": 0.31, "created_utc": "2026-09-14T03:10:00Z", "hostname": "..."
  }
}
```

`results` is `null` when `status != "ok"`.

### 4.10 Error handling

Any exception from pyDRAMA, an entry in `results["errors"]`, an empty
`results["results"]`, missing history files, or a timeout → JSON with
`status` `error`/`timeout`, `error` text including the last 40 lines of the
pyDRAMA/SESAM log, raw tree kept, exit code 1. On success the raw tree is
deleted unless `--keep-raw`.

## 5. Script 2 — `sphere_sweep.py`

### 5.1 CLI

```
python sphere_sweep.py
    [--parent "Generic_Satellite Reentry/output/dmf_output.json"] [--parent-object "Main Body"]
    [--outdir sphere_sweep_output] [--batch-size N] [--cores N] [--yes]
    [--diameters 5,10,...] [--temperatures 300,310,...] [--velocities 7.5,3.0,...]
    [--limit N] [--dry-run] [--force] [--retry-failed] [--timeout 600]
```

Defaults: the full grid (section 5.3), one batch of everything pending,
prompts on, `--outdir sphere_sweep_output` relative to the script.

### 5.2 Parent trajectory and state interpolation

1. Load `dmf_output.json`. Navigate
   `singleModuleOutputs.satellites[*].missionPhases[*].epochs[*]`; if there is
   more than one satellite, phase or epoch, abort with a message listing them.
   `epoch` is parsed from `%Y-%m-%dT%H:%M:%S.%fZ`.
2. In `analysisModules[0].results`, pick the key that ends with
   `_Trajectory.txt` and whose second dot-separated segment equals
   `--parent-object` with spaces replaced by underscores (`Main Body` →
   `Main_Body`). None or several matches → abort listing the available keys.
3. **Descending branch**: with `v[i]` the velocity column, let `j` be one past
   the last index where `v[i+1] > v[i]` (`j = 0` if velocity never increases).
   The branch is `rows[j:]`; it must have at least two rows, and velocity is
   strictly decreasing on it.
4. `v_min_parent = branch[-1].velocity`, `v_max_parent = branch[0].velocity`.
   **Log** (console, `sweep.log`, manifest): the parent file, object, epoch,
   the branch's row range, and
   `v_min = <v> km/s at t = <t> s, altitude <h> km (parent row <index>)`.
5. State at grid velocity `v`: find `i` with `v[i] >= v >= v[i+1]`,
   `f = (v[i] - v) / (v[i] - v[i+1])` (0 when the denominator is 0); linear
   interpolation of `altitude`, `lat`, `path`, `time`; angle-aware
   interpolation of `heading` and `lon` (interpolate sin and cos, then
   `atan2`, normalise heading to [0, 360) and longitude to [-180, 180));
   `epoch_run = parent_epoch + timedelta(seconds=time)`.
6. **Skip rule**: `v > v_max_parent` or `v < v_min_parent` → status `skipped`,
   reason `parent never reaches velocity`. With the default grid the lower
   bound is `v_min_parent` itself, so only points above `v_max_parent` can be
   skipped; explicit `--velocities` values may hit either bound.

### 5.3 Grid and matrix

- diameters: 5, 10, ..., 100 mm (20)
- temperatures: 300, 310, ..., 750 K (46)
- velocities: `numpy.linspace(7.5, v_min_parent, 100)` km/s (100, both ends
  included; the upper bound 7.5 km/s is fixed, the lower bound comes from the
  parent as in 5.2)
- `--diameters/--temperatures/--velocities` replace the corresponding list;
  `--limit N` keeps the first N points of the ordered matrix.

Points are ordered diameter → temperature → velocity (descending). Numeric
arguments are passed to the subprocess formatted with `f"{x:.6f}"` and the
epoch as `%Y-%m-%dT%H:%M:%S` (sub-second part dropped; pyDRAMA drops it
anyway). Each point gets its `run_name` by calling
`sphere_reentry.run_name()` on those formatted values parsed back, so the
sweep and script 1 agree on file names. Skipped points have no altitude and
therefore `run_name` and the inherited state fields are `null`.

### 5.4 Manifest and resume

`sweep_manifest.json` is written before any run and rewritten after every batch:

```json
{
  "created_utc": "...", "updated_utc": "...",
  "parent": {"path": "...", "sha256": "...", "object": "Main Body", "epoch_utc": "...",
             "branch_first_row": 469, "branch_last_row": 784,
             "v_min_kms": 0.028, "v_min_time_s": 3870.1, "v_min_altitude_km": 0.012,
             "v_max_kms": 7.918},
  "grid": {"diameters_mm": [...], "temperatures_K": [...], "velocities_kms": [...]},
  "settings": {"timeout_s": 600, "batch_size": null, "sphere_reentry_version": "1.0.0"},
  "points": [
    {"run_name": "...", "diameter_mm": 5.0, "temperature_K": 300.0, "velocity_kms": 7.5,
     "altitude_km": 77.5, "flight_path_deg": -0.959, "heading_deg": 347.168,
     "lat_deg": 29.546, "lon_deg": -82.134, "epoch_utc": "...",
     "status": "pending", "skip_reason": null, "returncode": null,
     "wall_time_s": null, "stderr_tail": null}
  ]
}
```

Statuses: `pending`, `skipped`, `done`, `failed`. On start-up every point whose
`runs/<run_name>.json` exists with `status: ok` becomes `done` (not re-run).
`--force` re-runs every non-skipped point; `--retry-failed` re-runs only
`failed` points. `--dry-run` prints matrix counts, the v_min log line, the
skipped points, the batch plan and the first/last run names, then exits
without running anything or writing the manifest.

### 5.5 Batched interactive execution and progress

`pending` points are split into batches of `--batch-size` (default: one batch
of all pending points). For each batch:

1. Print `Batch k/K: N runs, <first run_name> .. <last run_name>, estimated
   <time>` (rate from runs already completed in this invocation, else 0.3 s/run
   per core).
2. Prompt `Run this batch? [Y/n]` (`n`/`q` exits cleanly, manifest saved),
   then `Cores for this batch [<cpu_count - 1>]:` (integer in
   `1..os.cpu_count()`, re-prompt on invalid input). `--yes --cores N` answers
   both without prompting (`--yes` without `--cores` still prompts for cores).
3. Run the batch with `concurrent.futures.ThreadPoolExecutor(max_workers=cores)`;
   each task calls
   `subprocess.run([sys.executable, sphere_reentry.py, --velocity, ..., --altitude, ...,
   --temperature, ..., --diameter, ..., --flight-path-angle, ..., --heading, ...,
   --lat, ..., --lon, ..., --epoch, ..., --outdir, <outdir>/runs, --raw-dir, <outdir>/raw,
   --timeout, ..., --quiet],
   capture_output=True, timeout=timeout + 60)` and records return code, wall
   time and the last 20 lines of stderr.
4. One `tqdm` bar for the whole sweep: `total` = number of pending points at
   start-up, `initial` = 0, `desc` = "sweep", postfix = the most recently
   started run (`d=25mm T=450K v=3.12335km/s`); tqdm supplies elapsed, ETA and
   rate. The bar is closed/re-opened (with `initial` carried over) around the
   prompts so the prompt text is readable.
5. After the batch: print `ok / failed / wall time`, rewrite the manifest,
   regenerate `sweep_summary.csv`.

Failures never stop a batch. `KeyboardInterrupt` cancels queued tasks, lets
running subprocesses finish (each is ~0.3 s), saves the manifest and summary,
and exits with code 130.

### 5.6 Summary CSV

`sweep_summary.csv`, one row per matrix point including skipped ones:

```
diameter_mm, initial_temperature_K, initial_velocity_kms, initial_altitude_km,
flight_path_angle_deg, heading_deg, latitude_deg, longitude_deg, epoch_utc,
status, skip_reason, max_temperature_K, final_mass_kg, final_mass_source,
mass_loss_fraction, time_at_melting_temperature_s, final_velocity_kms,
final_radius_mm, final_altitude_km, end_of_life_reason, outcome, wall_time_s, run_name
```

Statistic columns are empty for `skipped`, `failed` and `pending` rows. The
summary is regenerated from the JSON files on disk after every batch (and by a
final pass), so it is always consistent with `runs/`.

### 5.7 Logging

`sweep.log` (append) receives: start-up parameters, the parent/v_min lines
from 5.2, per-batch prompts and answers, every failed run with its stderr
tail, and the end-of-sweep totals. Console shows the same at INFO level, above
the progress bar.

## 6. Testing

Framework: `pytest` in `drama_env`. Integration tests are marked
`@pytest.mark.drama` and `conftest.py` skips them automatically when `drama`
is not importable.

### 6.1 Fixtures (`tests/fixtures/`)

Real SESAM outputs captured from the spike (copied verbatim):

- `T1_demised_50mm/` — 50 mm, 300 K, 7.907 km/s from 101 km: history files
  (373 rows), `ImpactingFragments.xml` (empty), `sesam.log` (ends
  `uncritical`, `EVENT end ... 0.000000`)
- `T3_survivor_50mm/` — 50 mm from 39.9 km at 0.5 km/s: reaches the ground,
  XML with full-precision mass and velocity
- `T5_5mm_750K/` — 5 mm at 7.5 km/s: mass column all zeros, ends
  `uncritical`, `EVENT end ... 0.000003`
- `E1_ballooning_5mm/` — same as T5 with `energyThreshold = 1e-6`: ends
  `ballooning`
- `mini_dmf_output.json` — hand-built parent file: one satellite/phase/epoch,
  a `Main_Body` trajectory with an orbital oscillation (velocity rises and
  falls) followed by a strictly decreasing entry down to a non-zero impact
  velocity, a heading crossing 360°, plus a second object to test selection.

### 6.2 Unit tests — no DRAMA (`pytest tests/`)

`tests/test_sphere_reentry.py`

- `run_name`: exact string for known inputs; zero padding; lexicographic
  order matches parameter order; adjacent default-grid velocities give
  distinct names
- `sphere_mass_kg` / `sphere_cross_section_m2`: 5 mm → 1.8411e-4 kg,
  100 mm → 1.4729 kg, 50 mm → 1.9635e-3 m2
- `build_config`: geodetic elements in the right slots; heading normalised;
  `solid` True; material name; `energyThreshold == 1e-9`; both temperatures
  equal T0; epoch is a datetime and equals both `beginDate` and
  `initialDate`; no `materialList`; `objects` is a one-element list; mass and
  cross-section consistent with the diameter
- `read_sara_table`: row counts and column values from the fixtures; comment
  lines and malformed lines skipped; empty file → `[]`
- `parse_impacting_fragments`: T3 → mass `1.8411041946975187e-01`, velocity
  0.0582 km/s; T1 (no fragment) → `None`; missing file → `None`
- `parse_sesam_log`: reasons for T1 (`uncritical`), T3 (`ground impact`), E1
  (`ballooning`); `EVENT end` masses; missing file → `unknown`/`None`
- `merge_histories`: identical grids → one row per step with both column
  sets; deliberately shifted grid → outer join and a warning
- `compute_stats`: max temperature and where it occurs (T1); melt duration
  on T1 equals the sum of the 1 s steps at 850 K (rows both at melt) and
  first/last melt times/altitudes match the fixture; a synthetic re-melt
  series (melt, cool, melt) sums both intervals; tolerance boundary (849.4 K
  excluded, 849.6 K included); final-mass source priority: XML present → XML
  value and source; XML absent → log value; log absent → history value;
  final radius round-trips the initial radius when nothing melted and is 0
  for zero mass; `outcome` for each end-of-life reason: T3 → `survived`,
  E1 (ballooning) → `demised`, T5 (`uncritical`, 1.6 % residual) →
  `demised`, `uncritical` with intact mass → `other`; the coarse-CSV-mass
  warning fires for T5 and not for T1
- `write_csv`: header exactly as in 4.7, one row per merged row, values
  verbatim
- `run_sphere` error path with a stubbed `drama` module (`sys.modules`
  injection) whose `sara.run` returns an `errors` entry → JSON `status`
  `error`, raw dir kept; CLI returns exit code 1 in the same scenario
- CLI `--dry-run` prints the config and run name and creates no files;
  missing required args → exit code 2; negative diameter → exit code 2

`tests/test_sphere_sweep.py`

- parent loader on `mini_dmf_output.json`: epoch parsing; object selection by
  name (`Main Body` → `Main_Body` key); unknown object → error listing keys;
  duplicated satellite → error
- descending-branch detection: branch starts one past the last velocity
  increase; velocity strictly decreasing on it; `v_min`/`v_max` values and the
  logged v_min line contain the right time/altitude/row index
- interpolation: exact at row velocities; midpoint between rows; heading wrap
  (350° and 10° → 0°, not 180°); longitude wrap; epoch offset by the
  interpolated time
- skip rule: above `v_max` and below `v_min` → `skipped` with reason; a value
  equal to a bound is not skipped
- grid: 20 diameters, 46 temperatures, 100 velocities with endpoints 7.5 and
  `v_min_parent`; `--diameters/--temperatures/--velocities` subsets;
  `--limit`
- run names computed by the sweep equal `sphere_reentry.run_name()` applied
  to the parsed subprocess arguments
- batch splitting: default one batch; `--batch-size` remainder handling;
  `--limit` interplay
- resume: pre-existing ok JSON → `done`; failed JSON → `failed`; `--force`
  and `--retry-failed` selection
- manifest: written before execution, round-trips, updated statuses after a
  (fake) batch
- summary CSV from a directory of JSON fixtures (ok + failed + skipped rows;
  statistic columns empty where appropriate; header exactly as in 5.6)
- prompts with stubbed `input`: `n` aborts before submission and leaves the
  manifest with all points pending; invalid core counts re-prompt; `--yes
  --cores` skips both prompts; `--yes` alone still asks for cores
- subprocess command line for a point contains all six state arguments,
  `--outdir`, `--raw-dir`, `--timeout`, `--quiet`, and uses `sys.executable`
- executing a batch with a stubbed `subprocess.run` (fake script that writes
  an ok JSON or returns 1) updates statuses, wall times and stderr tails

### 6.3 Integration tests — `pytest -m drama` (need drama_env + DRAMA)

- script 1 real run: 50 mm, 300 K, 7.5 km/s at the parent's matching state →
  `status ok`, CSV first row equals the six inputs, all five statistics
  present and finite, `outcome demised`, raw tree deleted
- script 1 low-altitude survivor: 5 mm, 300 K, 0.3 km/s at ~33 km with the
  parent's state → `outcome survived`, `final_mass_source
  impacting_fragments_xml`, final mass equals initial to 1e-9, final radius
  2.5 mm
- sweep end-to-end against the real parent file: `--diameters 50
  --temperatures 300 --velocities 7.5,0.5,9.0 --yes --cores 2` → 2 done, 1
  skipped (9.0 > v_max), manifest and summary correct; second invocation
  reports nothing to do; `--dry-run` runs nothing

## 7. Non-goals and known limitations

- No plotting; analysis happens downstream from `sweep_summary.csv` and the
  per-run CSVs.
- The CSV `mass_kg` history is limited to SESAM's 0.001 kg print resolution;
  only the JSON final mass uses the precise sources.
- Runs that differ only in the optional state values share a run name.
- The parent's orbital phase (velocity oscillating above the grid's 7.5 km/s
  ceiling) is never sampled; only the final descending branch is used.
- pyDRAMA's `valuesFromFile` quirk is mirrored, not fixed.

## 8. Decision log

| Decision | Choice |
|---|---|
| Sphere size | diameter, 5–100 mm step 5 mm |
| Mass model | solid sphere, mass from drama-AA7075 density |
| Material | built-in `drama-AA7075` as shipped |
| Environment settings | mirror the parent run; `energyThreshold` lowered to 1e-9 J |
| Initial state | velocity from the grid; altitude, lat, lon, flight-path angle, heading, epoch inherited from the parent at that velocity |
| Velocity grid | `linspace(7.5, v_min_parent, 100)`, `v_min_parent` read from the parent file and logged |
| Skip rule | only velocities outside the parent's descending-branch range |
| Melt duration | sum of intervals with both rows within 0.5 K of 850 K |
| Final radius | from final mass; final mass from XML → log → history |
| Per-run output | one merged CSV + one JSON; raw DRAMA tree deleted on success |
| Sweep | subprocess per run, thread pool of N workers, batches with confirmation and per-batch core count, tqdm progress with ETA, resumable, aggregate summary |
| Location | `Particle Wake Evolution/`, git-initialised; outputs in `sphere_sweep_output/` |
| Tooling | `drama_env` + `pip install tqdm pytest` |
