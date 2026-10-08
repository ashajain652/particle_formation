# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Research code for a thesis chapter on space sustainability: what happens to the solid AA7075 fragments
that break off a re-entering satellite, how they heat up, melt, and spray molten droplets into the
wake. It has two halves that exist to be compared with each other.

1. **The SESAM wrapper** (`sphere_reentry.py`, `sphere_sweep.py`, `analysis/`) drives ESA DRAMA 4.1.4's
   re-entry module SESAM through the pyDRAMA package. SESAM is the established, lumped-parameter
   industry model; its runs are the *reference data*, not the product.
2. **The first-principles model** (`reentry_model/`) is the product: a trajectory and 3D
   finite-element thermal model of the same sphere, built step by step and verified against those
   SESAM references at every step.

Read `README.md` for the full command surface and the verification tables, and
`docs/model_assumptions.md` for the physics assumptions with their verification status. Those two
files are the contract with the thesis; keep them true when you change the model.

## Interpreters

There is no virtualenv in the repo and no `pip install -e .`; scripts are run with an absolute
interpreter path. `tests/helpers.py` hardcodes the first one.

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python     # everything, incl. pyDRAMA + SESAM
FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python   # only --thermal-solver fenicsx
```

`drama_env` is Python 3.12 with pyDRAMA, pymsis, scipy, cantera and `requirements-step2.txt`
(scikit-fem, gmsh, pyamg, pyvista, imageio-ffmpeg). `fenicsx_env` adds conda-forge `fenics-dolfinx`
0.11. On this Mac the FEniCSx environment needs two environment variables or it aborts:

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
```

## Tests

```bash
"$PY" -m pytest -m "not drama and not reference" -q      # the normal loop, ~4 min
"$PY" -m pytest tests/test_reentry_model_heating.py -q   # one module
"$PY" -m pytest tests/test_reentry_model_mesh.py -q -k prism   # one test
"$PY" -m pytest -m reference -q                          # SESAM reference flights: ~15 min (Step 1) + ~35 min (Step 2)
"$PY" -m pytest                                          # also the tests that launch real SESAM runs
```

Two markers, both declared in `pytest.ini`: `drama` (needs pyDRAMA importable; auto-skipped by
`tests/conftest.py` otherwise) and `reference` (the slow end-to-end flights against the committed
SESAM references in `data/reference_runs/`). The thermal tests dominate the unit-tier runtime; the
session-scoped mesh fixtures in `conftest.py` exist to keep it that low, so reuse them rather than
meshing again inside a test.

## Running the model

```bash
"$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
    --atmosphere us76 --thermal fem --heating physics --animate
"$PY" -m reentry_model compare --model reentry_model_output/<run>.csv --reference data/reference_runs/<sesam>.csv
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

Supporting modules: `mesh` (graded gmsh tetrahedra, cached by geometry; node coordinates are kept
mutable so a receding surface is possible later), `material` (DRAMA material JSON → k(T), c_p(T),
h(T) with a latent-heat hook), `thermal/` (two interchangeable P1/backward-Euler backends behind one
protocol — `skfem_backend` rescales precomputed element matrices and is ~3× faster in serial;
`fenicsx_backend` exists for MPI scaling and is verified to agree to every printed digit), `gas`
(Cantera equilibrium air + Blottner/Wilke viscosity), `earth`, `fap`, `sesam_io` (a SESAM run as a
`Reference` in SI), `compare` (model sampled at the reference's time stamps → metrics + six plots),
`viz` (PyVista animations read only the exported VTK series, never live objects).

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
- **`objects.xml` at the repo root is a SESAM side effect** written to the working directory
  regardless of `--outdir`; it is git-ignored, not a source file.

## Repository conventions

- **Work is spec → plan → implementation, under `docs/superpowers/`.** `specs/<date>-<name>-design.md`
  is the design; `plans/<date>-<name>.md` is the task-by-task implementation plan. Step 1 (trajectory)
  and Step 2 (coupled 3D heat transfer) are done and verified; **Step 3 (melting, melt film, Girin
  melt spraying) is planned but not implemented** — nothing in `reentry_model/` melts yet.
- **The Step 3 plan is generated, never hand-edited.** `prototype/` (code tracked since 2026-10-08, run outputs
  git-ignored) holds `proto3/`, a
  throwaway working copy of the package with Step 3 applied, and `plan3/make_plan.py`, which
  assembles `docs/superpowers/plans/2026-09-20-melt-spraying.md` by reading the code blocks straight
  out of `proto3/`. Regenerating must produce a byte-identical file. To change a code block, edit
  `proto3/`, run its tests there, regenerate, copy over. `plan3/audit_docs.py` checks the plan's prose
  blocks against the prototype's actual flags and column names. Read `prototype/README.md` first.
  `plans/melt-spraying-subplans/` splits that plan into one file per task for sub-agent refinement;
  `00-shared-context.md` is required reading before any of the others.
- **Outputs are git-ignored, inputs are committed.** `sphere_sweep_output/` (tens of GB),
  `reentry_model_output/`, the run outputs under `prototype/` and the DRAMA GUI folders stay out of git; the four SESAM
  reference runs in `data/reference_runs/`, the material JSONs, the fap space-weather files and the
  test fixtures in `tests/fixtures/` are committed because the verification depends on them.
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
