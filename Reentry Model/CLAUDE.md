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
  where `sesam_io.REFERENCE_DIR` and `fap.DEFAULT_FAP_DAY` point).
- `reentry_model/droplet_reentry/` — the droplet stage, a sub-package with its own `analysis/`
  (sweeps and plots) and `tests/` (see "The droplet stage" below).
- `analysis/` — the Step 1 and Step 2 verification tables and the mesh plot.
- `tests/` — the package's tests (`conftest.py` puts this folder on `sys.path` and holds the session
  mesh fixtures). `pytest.ini` also collects `reentry_model/droplet_reentry/tests/`.
- `output/` — every run, mesh cache, plot and verification output; git-ignored.
- `docs/` — `model_assumptions.md`, `milestones/`, and `superpowers/{specs,plans}`.
- `prototype/` — throwaway working copies (see Repository conventions).

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
so they cannot share one pytest session.

```bash
cd "Reentry Model"
"$PY" -m pytest -m "not reference" -q                     # the normal loop, ~3 min (droplet tests ~1 min)
"$PY" -m pytest tests/test_reentry_model_heating.py -q   # one module
"$PY" -m pytest tests/test_reentry_model_mesh.py -q -k prism   # one test
"$PY" -m pytest -m reference -q                          # SESAM reference flights: ~15 min (Step 1) + ~35 min (Step 2)

cd "DRAMA Sphere Experiments"
"$PY" -m pytest -m "not drama" -q                        # the wrapper and sweep, seconds
"$PY" -m pytest -m drama -q                              # also launches real SESAM runs

cd Spheral
"$PY" -m pytest -m "not drama and not reference" -q      # 'spheral' tests skip without the container
```

Markers, declared in each project's `pytest.ini`: `reference` (Reentry Model: the slow end-to-end
flights against the committed SESAM references in `DRAMA Sphere Experiments/data/reference_runs/`),
`drama` (DRAMA Sphere Experiments: needs pyDRAMA importable; auto-skipped by its `tests/conftest.py`
otherwise), and `spheral` / `fe_flight` (Spheral). The thermal tests dominate the Reentry Model
unit-tier runtime; the session-scoped mesh fixtures in its `conftest.py` exist to keep it that low, so
reuse them rather than meshing again inside a test.

## Running the model

From inside `Reentry Model/`:

```bash
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
    --atmosphere us76 --thermal fem --heating physics --animate
"$PY" -m reentry_model compare --model output/<run>.csv \
    --reference "../DRAMA Sphere Experiments/data/reference_runs/<sesam>.csv"
"$PY" analysis/reentry_model_verification.py          # Step 1 verification table
"$PY" analysis/reentry_model_thermal_verification.py  # Step 2 verification table
```

Exit codes are part of the interface: 0 ok, 1 a model/integration failure (message printed, including
a Newton non-convergence from either thermal backend), 2 bad arguments or a missing optional library.

## Architecture of `reentry_model`

The run is a lockstep coupling, not a monolith. `cli.py` builds three things and hands them to
`coupled.CoupledRun`:

- `trajectory.Simulator` — rotating-Earth 3-DOF in ECEF Cartesian coordinates, DOP853, with ground
  and escape as terminal events. It asks `atmosphere` for the freestream state and `aero` for C_D.
- a **body**, which is the extension point. `body.ConstantBody` is Step 1 (fixed mass, one
  temperature, `--thermal none`); `body.ThermalBody` wraps a mesh, a material and a conduction solver
  and advances a temperature field. The trajectory only ever sees mass, area and a temperature, so a
  future melting body slots in at the same seam.
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

### The droplet stage (`reentry_model/droplet_reentry/`)

The independent demise of the molten droplets sprayed into the wake — trajectory, heating,
aluminium/zinc/magnesium evaporation, the Singh et al. (2021) sphere drag law, mean vertical winds and
settling to the ground. Promoted on 2026-10-26 from `prototype/work-2026-10-21-droplet-demise/`, whose
`REPORT.md`, `DATA_INDEX.md` and `VERTICAL_WINDS.md` still hold its development and results.

- `model.py` — `DropletSimulator` (a `trajectory.Simulator` with temperature and element masses in the
  state) and `run_case(d0_um, h0_km)`, the entry point for one droplet. `run_case(d0_um, release=state)`
  starts from a given `InitialState` instead of the parent's state at `h0_km`; that is the seam for the
  fragment model, which will release each droplet with its velocity at separation (Asha, 2026-10-24).
  The connection itself waits until melting (Step 3) joins the package.
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
  and Step 2 (coupled 3D heat transfer) are done and verified; **Step 3 (melting, melt film, Girin
  melt spraying) is planned but not implemented** — nothing in `reentry_model/` melts yet.
- **The Step 3 plan is generated, never hand-edited.** `prototype/` (code tracked since 2026-10-08,
  run outputs git-ignored) holds `proto3/`, a throwaway working copy of the package with Step 3
  applied, and `plan3/make_plan.py`, which assembles `docs/superpowers/plans/2026-09-20-melt-spraying.md`
  (as `plan3/plan.md`) by reading the code blocks straight out of `proto3/`. Regenerating must produce
  a byte-identical file. To change a code block, edit `proto3/`, run its tests there, regenerate, copy
  over. `plan3/audit_docs.py` checks the plan's prose blocks against the prototype's actual flags and
  column names. Read `prototype/README.md` first. `docs/superpowers/plans/melt-spraying-subplans/`
  splits that plan into one file per task for sub-agent refinement; `00-shared-context.md` is
  required reading before any of the others. The prototype copies keep the old flat layout inside
  themselves and were not re-pathed in the 2026-10-24 restructure.
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
