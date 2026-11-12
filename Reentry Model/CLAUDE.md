# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Research code for a thesis chapter on space sustainability: what happens to the solid AA7075 fragments
that break off a re-entering satellite, how they heat up, melt, and spray molten droplets into the
wake. The repository root holds three projects, each in its own folder with its own `tests/` and
`pytest.ini`, plus the DRAMA GUI project they share:

1. **`DRAMA Sphere Experiments/`** — the SESAM wrapper (`sphere_reentry.py`, `sphere_sweep.py`,
   `analysis/`) drives ESA DRAMA 4.1.4's re-entry module SESAM through the pyDRAMA package. SESAM is
   the established, lumped-parameter industry model; its runs are the *reference data*, not the
   product. Its `data/` holds the four committed SESAM reference runs (`data/reference_runs/`), the
   DRAMA material JSONs (`data/user_materials/`) and the fap space-weather files.
2. **`Reentry Model/`** (this folder) — the first-principles model, the product: the `reentry_model`
   package, a trajectory and 3D finite-element thermal model of the same sphere, built step by step
   and verified against those SESAM references at every step.
3. **`Spheral/`** — the `spheral_frag` package: the fragmentation study run in the Spheral SPH code
   (container under `spheral_frag/container/`), with its own `analysis/`, `data/` and `docs/`.

`Generic_Satellite Reentry/` at the root is the DRAMA GUI project of the parent 3U CubeSat; its
`output/dmf_output.json` is the parent trajectory the sweeps and the droplet work start from. The GUI
rewrites that folder on every run, so the droplet model reads a committed snapshot of the 2026-09-13
run, `Reentry Model/data/generic_satellite_reentry/dmf_output.json`.
`citations/` at the root holds source notes shared by the projects.

Read the root `README.md` for the full command surface and the verification tables, and
`Reentry Model/docs/model_assumptions.md` for the physics assumptions with their verification status.
Those two files are the contract with the thesis; keep them true when you change the model.

### Layout of `Reentry Model/`

- `reentry_model/` — the package (`import reentry_model` needs `Reentry Model/` on `sys.path`; run
  everything from inside this folder). Its `__init__.py` defines the folders the code uses:
  `PROJECT_ROOT` (this folder), `DATA_DIR` (`data/`: the US76 table, the ATDB sphere drag table, the
  material JSON), `OUTPUT_DIR` (`output/`), and `DRAMA_DATA_DIR` (`DRAMA Sphere Experiments/data/`,
  where `sesam_io.REFERENCE_DIR`, the non-melting SESAM reference runs, and `fap.DEFAULT_FAP_DAY` point;
  the two melting SESAM reference runs are in `data/reference_runs/`, which `sesam_io.MELTING_REFERENCE_DIR`
  points at).
- `reentry_model/droplet_reentry/` — the droplet stage, a sub-package with its own `analysis/`
  (sweeps and plots) and `tests/` (see "The droplet stage" below).
- `reentry_model/wake_column/` — the wake column, a sub-package with its own `analysis/` and `tests/`: the
  chemistry and condensation of the aluminium the droplets vaporise (see "The wake column" below). Its committed
  inputs are in `data/wake_column/`.
- `analysis/` — the Step 1 and Step 2 verification tables and the mesh plot; the Step 3 scripts (`melt_verification.py`, `melt_sensitivity.py`, `girin_reference.py`, `flow_wavelengths.py`, `skin_plate.py`, `skin_timing.py`, `ablation_1d.py`, `movemesh_recession.py`); and `analysis/skyscrapers_port/`, the port's own tools (path map, checks, three-way resolver, `same_results.py`) with its README.
- `tests/` — the package's tests (`conftest.py` puts this folder on `sys.path` and holds the session
  mesh fixtures). `pytest.ini` also collects `reentry_model/droplet_reentry/tests/` and
  `reentry_model/wake_column/tests/`.
- `output/` — every run, mesh cache, plot and verification output; git-ignored.
- `docs/` — `model_assumptions.md`, `milestones/`, and `superpowers/{specs,plans}`.
- `prototype/` — the record of the working copies the Step 3 code was developed in; the code of record is now
  `reentry_model/` (see Repository conventions).

`compare.py` borrows the plot styling from `DRAMA Sphere Experiments/analysis/plot_style.py` and falls
back to plain matplotlib without it.

## Interpreters

There is no virtualenv in the repo and no `pip install -e .`; scripts are run with an absolute
interpreter path. Each project's `tests/helpers.py` hardcodes the first one.

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python     # everything, incl. pyDRAMA + SESAM
FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python   # only --thermal-solver fenicsx
```

`drama_env` is Python 3.12 with pyDRAMA, pymsis, scipy, cantera, scikit-fem, gmsh, pyamg, pyvista and
imageio-ffmpeg. `fenicsx_env` adds conda-forge `fenics-dolfinx` 0.11. On this Mac the FEniCSx
environment needs two environment variables or it aborts (run from `Reentry Model/`):

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
```

## Tests

Each project is tested from inside its own folder; the three `tests/` folders each have a `helpers.py`,
so they cannot share one pytest session. In `Reentry Model`, `pytest tests` runs the package's own tests; the
default `pytest` (from `pytest.ini`) also collects the droplet and wake-column tests. After the skyscrapers
port the unit tier, `pytest tests -m "not drama and not reference" -n 4`, gave 702 passed and 1 skipped in about
24 minutes across four processes (the skip is the FEniCSx module, which needs `dolfinx`). The port's own tools
have their tests in `analysis/skyscrapers_port/`, run explicitly (see `analysis/skyscrapers_port/README.md`).

```bash
cd "Reentry Model"
PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4
                                                          # the package's unit tier: 702 passed, 1 skipped, ~24 min on four
                                                          # processes (needs pytest-xdist)
"$PY" -m pytest -m "not reference" -q                     # the default collection, which adds the droplet and wake tests: longer (not timed)
"$PY" -m pytest tests/test_reentry_model_heating.py -q   # one module
"$PY" -m pytest tests/test_reentry_model_mesh.py -q -k prism   # one test
"$PY" -m pytest -m reference -q                          # SESAM reference flights: ~15 min (Step 1) + ~35 min (Step 2); the Step 3 reference tests (below) are
                                                          # collected too but have not yet been run on the port (run time not measured)

cd "DRAMA Sphere Experiments"
"$PY" -m pytest -m "not drama" -q                        # the wrapper and sweep, seconds
"$PY" -m pytest -m drama -q                              # also launches real SESAM runs

cd Spheral
"$PY" -m pytest -m "not drama and not reference" -q      # 'spheral' tests skip without the container
```

Markers, declared in each project's `pytest.ini`: `reference` (Reentry Model: the slow end-to-end
flights against the committed SESAM references in `DRAMA Sphere Experiments/data/reference_runs/`, and the Step 3 melting
references in `Reentry Model/data/reference_runs/`, which `sesam_io.MELTING_REFERENCE_DIR` points at: `tests/test_reentry_model_reference_melt.py` and two tests in `tests/test_reentry_model_skinbody.py`; these Step 3 reference tests were collected but have not yet been run on the port),
`drama` (DRAMA Sphere Experiments: needs pyDRAMA importable; auto-skipped by its `tests/conftest.py`
otherwise), and `spheral` / `fe_flight` (Spheral). The thermal tests dominate the Reentry Model
unit-tier runtime; the session-scoped mesh fixtures in its `conftest.py` exist to keep it that low, so
reuse them rather than meshing again inside a test.

## Running the model

From inside `Reentry Model/`:

```bash
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
    --atmosphere us76 --thermal fem --surface-model elements --heating physics --animate
# a short melting smoke flight (65 km, 0.125 s), not a production one:
"$PY" -m reentry_model run --diameter 100 --velocity 7.24 --altitude 65 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --heating physics --h-surface 4 --h-core 20 --melt on --temperature 893 --dt 0.025 --t-max 0.125 \
    --skin-preset dev --fine-size 2 --frames-every 0
"$PY" -m reentry_model compare --model output/<run>.csv \
    --reference "../DRAMA Sphere Experiments/data/reference_runs/<sesam>.csv"
"$PY" analysis/reentry_model_verification.py          # Step 1 verification table
"$PY" analysis/reentry_model_thermal_verification.py  # Step 2 verification table
```

Exit codes are part of the interface: 0 ok, 1 a model/integration failure (message printed, including
a Newton non-convergence from either thermal backend), 2 bad arguments or a missing optional library.

## Two model architectures and the shared package

The model is built in two architectures that share one command and one package. The **skyscrapers** architecture, the
default, is the mesh model described in this file; the **level set** architecture is designed in
`docs/superpowers/specs/levelset/` and is not yet built. `reentry_model run --architecture skyscrapers|levelset`
chooses between them, and `levelset` exits with code 2 until its run driver exists; `resume` and `compare` do not take
the flag.

Code that both architectures use lives in `reentry_model/shared/`: the Step 1 and Step 2 modules (`trajectory`,
`earth`, `atmosphere`, `fap`, `aero`, `gas`, `constants`, `heating`, `material`, `sesam_io`, `compare`, `viz`), the
Step 3 physics that does not depend on a mesh (`rheology`, `dispersion`, `melt_profile`, `surface_flow`, `spray`,
`girin_case`), the droplet writer (`droplets`), the melt-film lubrication formulas (`lubrication`), the body seam (`body`), the thermal protocol (`thermal`) and the run command's common
arguments (`cli_common`). The old paths of the 17 wholly moved modules (such as `reentry_model.aero`) are alias files that make
the old name the same module object as the shared one. `skyscrapers/body.py`, `skyscrapers/coupled.py`, `thermal/__init__.py` and `cli.py`
keep their skyscrapers-only code and only re-export the moved shared names, so patching a shared name there (for example
`cli.git_commit`) does not reach the shared code; patch it in `reentry_model.shared` instead. Existing imports and
old pickles keep working, but states saved after this change pickle their classes as `reentry_model.shared.*` and
cannot be loaded by a checkout from before it (such as `subplan20`, or a Mac Studio checkout that is behind). New code imports from
`reentry_model.shared` directly. The shared package imports nothing from the rest of `reentry_model` except the
package-root constants (`PROJECT_ROOT`, `DATA_DIR`, and so on), which `tests/test_reentry_model_shared.py` checks both
by reading the source and by importing every shared module in a fresh interpreter. The skyscrapers-only modules live in `reentry_model/skyscrapers/` (`body`, `coupled`, `columns`, `mesh`, `movemesh`,
`restart`, `skin`, `skinbody` and `spreading`, plus `elements/film.py`, the element model's runoff, `--surface-model
elements`); their old paths (such as `reentry_model.skin`) are alias files too, so existing imports keep working, but states
saved after this change pickle `reentry_model.skyscrapers.*` and cannot be loaded by a checkout from before it. `cli.py`,
`__main__.py` and `thermal/` stay at the top level. The skins and the backends are never shared between the architectures. The refactor changed no model behaviour (see
`docs/model_assumptions.md`).

## Architecture of `reentry_model`

The run is a lockstep coupling, not a monolith. `cli.py` builds three things and hands them to
`coupled.CoupledRun`:

- `trajectory.Simulator` — rotating-Earth 3-DOF in ECEF Cartesian coordinates, DOP853, with ground
  and escape as terminal events. It asks `atmosphere` for the freestream state and `aero` for C_D.
- a **body**, which is the extension point. `body.ConstantBody` is Step 1 (fixed mass, one
  temperature, `--thermal none`); `body.ThermalBody` wraps a mesh, a material and a conduction solver
  and advances a temperature field. The trajectory only ever sees mass, area and a temperature, so a
  melting body slots in at the same seam: `body.MeltingBody` (the element model's) and `skinbody` (the layered
  skins') are the melting bodies, in the package since the port of 2026-11-11.
- `heating.*` — turns the freestream state plus every surface patch's angle θ and wall temperature
  into a per-patch convective flux.

Each macro step (`--dt`, 0.5 s): advance the trajectory, evaluate heating at the end-of-step aero
state with the start-of-step wall temperatures, advance the conduction (radiation implicit), write a
history row. That is first-order operator splitting; the Δt-halving test is what bounds its error.

Supporting modules: `mesh` (graded gmsh tetrahedra, cached by geometry in `output/meshes/`; node
coordinates are kept mutable so a receding surface is possible later), `material` (DRAMA material JSON
→ k(T), c_p(T), h(T) with a latent-heat hook), `thermal/` (two interchangeable P1/backward-Euler
backends behind one protocol — `skfem_backend` rescales precomputed element matrices and is ~3× faster
in serial; `fenicsx_backend` exists for MPI scaling and is verified to agree to every printed digit),
`gas` (Cantera equilibrium air + Blottner/Wilke viscosity), `earth`, `fap`, `sesam_io` (a SESAM run as
a `Reference` in SI), `compare` (model sampled at the reference's time stamps → metrics + six plots),
`viz` (PyVista animations read only the exported VTK series, never live objects).

Step 3, melting and melt spraying, joined the package on 2026-11-11 (see Repository conventions). Its modules, each
described from its own docstring. The layered-skin design (the "skyscrapers" architecture):

- `skin` — a layered skin of cells on every corner column of the 3D mesh's surface (stage 2 heat only, then the flow).
- `skinbody` — the body that couples the moving 3D mesh to those skins each macro step.
- `columns` — the corner columns, each owning the dual cell of one corner of the surface triangulation.
- `movemesh` — the moving 3D mesh: surface corners move by the metal crossing the rigid front, nodes in a deforming
  zone beneath follow, and the temperature field is remapped with the energy books closed exactly.
- `spreading` — spreading the re-frozen deposit where the base would fold.
- `restart` — saving a flight's whole state to disk and resuming from it.

- `elements/film` — the element model's `Runoff`, the upwind runoff on the patch graph (`--surface-model elements`); it
  re-exports `lubrication` and `deep_flux`.

These modules sit in `reentry_model/skyscrapers/` since 2026-11-11 (the old paths are aliases).

Physics shared with the other (level-set) architecture:

- `lubrication` — melt-film lubrication flow (`lubrication()` and `deep_flux()`).
- `spray` — Girin's gradient instability per patch and the droplet release bookkeeping.
- `dispersion` — the numerical solution of Girin's dispersion relation.
- `melt_profile` — the melt layer's speed profile fitted as straight segments on plain arrays (`LayerColumns`, which the
  skins fill in `skinbody.layer_columns`), with the kinematic viscosity below each breakpoint and the layer's Couette
  viscosity, for the broken-line solver and the floor test in `dispersion`.
- `surface_flow` — the flow-regime gate, wall pressure, boundary-layer edge state, wall shear and the film's driving gradient.
- `rheology` — the viscosity of liquid and slurry AA7075.
- `girin_case` — Girin's published spraying cases run with his own simplifications, as a validation.

### The droplet stage (`reentry_model/droplet_reentry/`)

The independent demise of the molten droplets sprayed into the wake — trajectory, heating,
aluminium/zinc/magnesium evaporation, the Singh et al. (2021) sphere drag law, mean vertical winds and
settling to the ground. Promoted on 2026-10-26 from `prototype/work-2026-10-21-droplet-demise/`, whose
`REPORT.md`, `DATA_INDEX.md` and `VERTICAL_WINDS.md` still hold its development and results.

- `model.py` — `DropletSimulator` (a `trajectory.Simulator` with temperature and element masses in the
  state) and `run_case(d0_um, h0_km)`, the entry point for one droplet. `run_case(d0_um, release=state)`
  starts from a given `InitialState` instead of the parent's state at `h0_km`; that is the seam for the
  fragment model, which will release each droplet with its velocity at separation (Asha, 2026-10-24).
  Melting (Step 3) joined the package on 2026-11-11; the connection itself is still not built.
- Oxidation (2026-10-27, a bounding model, off by default): `run_case(..., sticking=γ, f_skin=f)` adds
  Genge (2016)-style oxygen capture at sticking efficiency γ, Al → Al₂O₃ heat into the droplet, and two
  oxide masses in the state (skin and entrained). The skin closes the surface (evaporation and further
  oxidation) over a few nanometres. With `sticking=0` the state and the results are the old model
  exactly. `analysis/oxidation_sweep.py` and `plot_oxidation_sweep.py` sweep γ and f.
- Copies of its dependencies, so the package stands alone: `particle_trajectory.py` and
  `droplet_heating_sweep.py` from `prototype/work-2026-10-14-particle-settling/`, and the parent-
  trajectory section of `DRAMA Sphere Experiments/sphere_sweep.py` as `parent_trajectory.py`. Only
  their imports and paths changed; `singh_drag.py` moved in with the model.
- `analysis/` — the demise and terminal-velocity sweeps and their plots, run with
  `"$PY" -m reentry_model.droplet_reentry.analysis.<name>` from this folder; everything they write
  goes to `output/droplet_reentry_output/` (`droplet_reentry.OUTPUT`).
- `tests/` — the prototype's checks as pytest tests, plus a reproduction test that re-runs five cases
  of the 2026-10-22 demise sweep and must match them to 1e-9 (the rows are in `tests/fixtures/`). The
  full re-run of that sweep from the package matched every case.

### The wake column (`reentry_model/wake_column/`)

Follows the aluminium a droplet vaporises (about 45 % of a droplet's mass at most) from hot vapour to alumina
nanoparticles over 12 hours. Built 2026-10-27 from `docs/superpowers/specs/2026-10-27-wake-column-chemistry-design.md`
(the design, with the physics decided with Asha) and its plan in `docs/superpowers/plans/`; the assumptions, each marked
verified (with its test) or assumed, are `docs/model_assumptions.md` section 9. Read that section before changing the
physics or quoting a result.

- Pipeline, one `Case(d0_um, h0_km, time_of_day)` at a time, chained by `run.run`: the droplet run
  (`droplet_reentry`, not modified) gives the trail (`trail.py`: length L is the path containing 90 % of the aluminium
  loss, the inventory, the mean altitude, the emission samples); Stage A (`slowing.py`) slows the hot vapour to
  ambient speed with a test-particle Monte Carlo and returns the starting Al/AlO/OAlO fractions and the radius r0;
  Stage B (`box.py`) is one stiff ODE system (scipy BDF) of extensive counts in a diffusing cylinder:
  entrainment, the Plane et al. (2021) neutral network, sectional particles (bins by aluminium atoms, `particles.py`)
  with tallies by source molecule and route, uptake onto Megner smoke, coagulation and settling. `transport.py`,
  `ambient.py`, `species.py`, `reactions.py` and `settings.py` hold diffusion and geometry, the ambient columns, and
  the tables and settings. By default the trail is 24 equal-length segments (`segments.py`), each with its real
  aluminium, its own Stage A and box, added at a common time since release, with spreading along the trail applied
  afterwards (`along_trail.csv`); `Settings.n_segments = 0` restores one column.
- Data-driven: species, reactions, ambient columns (WACCM baseline, NRLMSIS + AFGL cross-check) and the Megner eddy
  and smoke tables are committed in `data/wake_column/`, each with a source and its assumed values flagged. No
  unmeasured rate or efficiency is hard-coded: each is a named field of `Settings` (formation, sticking, oxidation and
  smoke-sticking efficiencies, `hot_start`, `hot_reaction_scale`, `eddy_scale`, `epsilon_scale`, `richardson_g`, `bin_ratio`, ...).
- Commands, from this folder; everything goes to `output/wake_column_output/` (`wake_column.OUTPUT`, git-ignored):

  ```bash
  "$PY" -m reentry_model.wake_column.analysis.sweep --ambient waccm      # 126 cases (21 altitudes x 1, 14.1, 200 um x day, night); --ambient msis_afgl for the cross-check
  "$PY" -m reentry_model.wake_column.analysis.sensitivity                # 15 settings x 6 cases (14.1 um at 60, 76, 90 km, day, night) = 90 runs
  "$PY" -m reentry_model.wake_column.analysis.plot_sweep                 # plots into <out>/plots/; plot_sensitivity likewise (tornado plots)
  ```

  `sweep` and `sensitivity` take `--workers N` (default CPU count less one) and `--out DIR`; a run whose `summary.json`
  exists is skipped, so a stopped sweep resumes, and a run that raises is reported at the end without stopping the rest.
  The plot scripts choose baseline runs by their physics settings (`analysis/results.py`), so test runs and the two
  ambients can share a folder. Runs that vaporised no aluminium (`skipped` in the summary) have no outcome keys.
- Run names encode the case, the ambient and every non-default setting (`settings.run_name`); keep that when adding a
  setting, and update `test_every_setting_changes_the_run_name`.
- Things that are easy to get wrong: `hot_start="droplet_temperature"` is Stage A with the droplet speed zeroed, not
  "no hot chemistry"; "condensed" includes uptake onto background smoke and settled matter, and at 76 km nearly all of
  the 12-hour condensed fraction is smoke uptake, which rests on an assumed 0.2 nm smoke bin (the sensitivity set has
  "smoke sticking x0.1" and "no smoke uptake" for that); `d_volume_mean_nm` is the volume-weighted mean diameter
  (sum N d^4 over sum N d^3), not the mean-volume diameter; the oxygen-uptake factor is a cubic smoothstep, and the
  default bin ratio is 2^(1/2), both chosen to make BDF converge and the coagulation test pass (assumptions 9.4);
  `Settings.particles=False` also switches off smoke uptake. The WACCM ambient needs `data/wake_column/` only; the
  1.6 GB WACCM file is not in the repository (`analysis/extract_waccm_columns.py` re-extracts the committed table).
- A 12-hour run with 24 segments takes about 1-3 minutes; the default sweep's 126 cases took 65 minutes on three workers.

### Things that are easy to get wrong here

- **`--heating sesam` is a verification device, not physics.** It spreads SESAM's surface-averaged
  heat input uniformly over every patch, front and back, so the conduction and coupling can be
  compared with SESAM's lumped temperature with no distribution question in between. Never quote its
  results as a physical prediction; `--heating physics` (Fay–Riddell + Matting bridging + Lees
  distribution) is the model proper.
- **Some agreement is in-sample by construction.** The Knudsen drag bridging `aero.SesamTable`, the
  heat factor `aero.SesamHeatTable` and the hot-wall c_p were all fitted to SESAM output that
  includes the four committed references. The independent content of the verification is the
  trajectory (replay mode) and the temperature, not C_D or the heat factor.
- **pyDRAMA writes booleans into `sara.xml` as `True`/`False`, which SESAM does not understand** — it
  wants `yes`/`no`. The wrapper now writes the three affected switches as strings; the whole original
  sweep therefore ran on the static US76 table with no wind, whatever its settings said. Check
  `sesam.log`, which states the environment it actually used.
- **`objects.xml` is a SESAM side effect** written to whatever directory SESAM was launched from,
  regardless of `--outdir`; it is git-ignored everywhere, not a source file.

## Repository conventions

- **Work is spec → plan → implementation, under each project's `docs/`.** For this project,
  `docs/superpowers/specs/<date>-<name>-design.md` is the design and
  `docs/superpowers/plans/<date>-<name>.md` the task-by-task implementation plan (DRAMA Sphere
  Experiments and Spheral keep theirs in `docs/specs/` and `docs/plans/`). Step 1 (trajectory)
  and Step 2 (coupled 3D heat transfer) are done and verified. **Step 3 (melting, melt film, Girin melt
  spraying) and the skyscrapers architecture were ported into the package on 2026-11-11** (merge commits 592a77d,
  the prototype line's designs, plans and specs, and d716a24, the build repository's 96 commits rewritten into
  this layout; the data paths were repointed in b92c4e9). The port changed no model behaviour relative to the prototype. Relative to the old `main`, the Step 2 path changed as the prototype had already changed it: the layered skins are the default surface model with `--thermal fem` (non-melting runs need `--surface-model elements`), the capacity is lumped by default, and temperatures differ by about 1 % (stagnation temperature 4.5 K, 1.2 %, at most, on the 20 s non-melting comparison of 2026-11-11; see the handoff note). Whether that default should change is Asha's decision and is pending. How the port was made and
  checked is in `analysis/skyscrapers_port/README.md` and its design
  `docs/superpowers/specs/2026-11-11-skyscrapers-port-to-main-design.md`.
- **Next is the refactor of 2026-10-12, and it is another agent's.** `docs/superpowers/plans/shared/README.md`
  and the level-set design's section 14 record the decision: one combined refactor moves the shared physics into
  `reentry_model/shared/`. After it, skyscrapers work goes on short branches off `main`, with no prototype
  folder (Asha, 2026-11-11). Until then do not restructure the new modules.
- **The Step 3 plan was generated from the prototype; that workflow is retired.** Before the port, `prototype/`
  held `proto3/`, a throwaway copy of the package with Step 3 applied, and `plan3/make_plan.py`, which assembled
  `docs/superpowers/plans/2026-09-20-melt-spraying.md` from it. Those, and the six skyscrapers work copies
  (`work-2026-10-08-skin`, `work-2026-10-11-moving-mesh-plan`, `work-2026-10-12-moving-mesh`,
  `work-2026-10-12-skin-heat-plan`, `work-2026-10-18-column-flow-plan`, `work-2026-10-28-skin-spraying-plan`),
  stay as the record and are not edited; each of the six has a README saying so. Read `prototype/README.md`
  first. `docs/superpowers/plans/melt-spraying-subplans/` holds one file per task; `00-shared-context.md` is
  required reading before any of the others. The history is reachable from the tags
  `skyscrapers-prototype-end` and `skyscrapers-build-end`.
- **Outputs are git-ignored, inputs are committed.** Ignored: `Reentry Model/output/`, the run outputs
  under `Reentry Model/prototype/`, `sphere_sweep_output/` (tens of GB), `spheral_output/` and the
  DRAMA GUI folders. Committed because the verification depends on them: the four SESAM reference
  runs, the material JSONs and fap files in `DRAMA Sphere Experiments/data/`, this project's `data/`,
  `Spheral/data/`, and every project's `tests/fixtures/`. Ignore outputs by folder in `.gitignore`,
  never by file type, so a new input CSV is never silently left untracked.
- **Run names encode the whole configuration** (`sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind`)
  so that runs of different materials, atmospheres and wind settings can share one output directory
  without overwriting each other, and the sweep stays resumable by checking for the JSON. Preserve
  that property when adding a setting.
- **Every claim in `README.md` and `docs/model_assumptions.md` carries its measurement.** When a
  number changes, update the verification table, the assumption entry and the threshold in
  `tests/test_reentry_model_reference*.py` together, and mark what was verified against SESAM or an
  analytic solution versus what was assumed.
- **Commit messages** are a sentence in the imperative or a scope prefix and a sentence —
  `reentry_model: cross-section temperature video ...`, `docs: ...`, `Step 3 plan: ...` — describing
  the physics or behaviour that changed, not the files touched.
