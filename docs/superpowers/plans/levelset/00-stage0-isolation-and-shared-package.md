# Level-set stage 0 — the shared package and the architecture flag Implementation Plan

> **Superseded on 2026-10-12 (Asha's decision).** The shared-package refactor now waits for the skyscrapers port of the
> prototyped Step 3 code to land on `main`, and is then done as one combined refactor (stage 0 + S) with its own plan,
> written at that time (spec §4.4, §14). This plan was executed on branch `levelset-stage0` (head da4c549; unit 377
> passed, FEniCSx 7, reference 13, bit-identity 3/3), which is **not** merged. Its code — the alias mechanism, the guard
> tests, `prototype/levelset/stage0/bitid.sh` and the `--architecture` flag — is the tested reference the combined plan
> reuses.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the architecture-independent modules of `reentry_model` into a new package `reentry_model/shared/`,
keep every old module path working as an alias, and add `run --architecture skyscrapers|levelset` (levelset exits 2
until built) — with every flight history byte-identical before and after.

**Architecture:** Twelve modules move whole into `reentry_model/shared/` (`git mv`), and each old path becomes a
four-line alias that replaces itself in `sys.modules` with the shared module object, so `reentry_model.aero is
reentry_model.shared.aero`. Three modules are split: the `Body` protocol, `sphere_mass` and `ConstantBody` go to
`shared/body.py`; the Step 2 thermal protocol to `shared/thermal.py`; the run command's flight arguments, defaults,
atmosphere factory, run name, provenance and the `--architecture` pre-parse to `shared/cli_common.py`. The skyscrapers
modules (`body.py`, `thermal/__init__.py`, `cli.py`, `coupled.py`) import from `shared` and re-export the names they
used to define.

**Tech Stack:** Python 3.12 (`drama_env`), pytest, argparse; no new dependency.

**Spec:** `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4 (isolation, the shared
package, the flag, coordination) and §14 (stage 0). The skyscrapers sub-plans carry an "Amendment of 2026-10-11 — the
shared physics package" pointing at this stage; `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`
ends with the matching amendment.

**Prototype:** every code block below was applied and tested in a scratch worktree of `main` (d16018b) on 2026-10-11.
The tested change is `prototype/levelset/stage0/stage0.patch` (against `main`), the bit-identity script
`prototype/levelset/stage0/bitid.sh`. Measured there: the three bit-identity histories identical byte for byte; the unit
tier 352 → 375 tests selected (with the plan's final tests), all passing except the FEniCSx module (skipped in `drama_env`); `tests/test_reentry_model_fenicsx.py`
7 passed in `fenicsx_env`; the reference tier 13 passed, 1 skipped (the FEniCSx module) in 18 min. The history
CSVs print about nine significant digits, and the Step 2 flights are not bit-reproducible run to run in their
full-precision JSON summaries (measured on the unchanged base), so "byte-identical" means identical to the CSV's printed
precision for those flights.

## Global Constraints

- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; FEniCSx tests need
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang /Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`.
- **Behaviour-preserving**: histories byte-identical (Task 1's baseline, Task 5's check); the unit tier
  (`"$PY" -m pytest -m "not drama and not reference" -q`) and the reference tier (`"$PY" -m pytest -m reference -q`) pass.
- **Work on branch `levelset-stage0` from `main`, in its own git worktree** — another session works in the main
  checkout on `step3-rigid-substrate-scheil`; never switch that checkout's branch (spec §4.4).
- Skyscrapers run names, defaults, help texts of existing flags, JSON keys and exit codes do not change.
- `reentry_model/shared/` imports nothing from `reentry_model` outside itself except `DATA_DIR` and `__version__` from
  the package root (spec §4.2).
- The skins are never shared (spec §4.2); no skin, mesh, coupled-loop or thermal-backend code moves.
- Commit messages: a sentence describing the behaviour or structure that changed, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Do not edit the package while a measurement (Task 1, Task 5) is running (Step 3 fact 31).

## Review Focus

- **Analysis scripts and tests that import old paths** (`from reentry_model import cli, sesam_io`, `mesh as meshmod`):
  must keep importing — pinned by `test_every_reentry_model_import_in_tests_and_analysis_resolves` (Task 5).
- **Import order**: importing `reentry_model.shared.X` before `reentry_model.X` (a script that adopts the new path first)
  must still give one module object — pinned by `test_import_order_does_not_matter` (Task 2).
- **`--architecture=levelset` in the equals form** must route the same as `--architecture levelset` — pinned by the
  parametrised `test_architecture_levelset_exits_2_until_built` (Task 4).
- **`--architecture` given to `compare`** must be rejected by argparse (exit 2), not silently routed — pinned by
  `test_architecture_does_not_reach_compare` (Task 4).
- **Data files located relative to a module's own file** (`DATA_DIR`, `REPO_ROOT`, the `analysis/` path `compare.py`
  puts on `sys.path`) must still resolve one directory deeper — pinned by `test_data_and_repository_paths_resolve`
  (Task 2) and `test_cli_common_reexports_and_paths` (Task 4).

---

## File structure

| Path | Change | Responsibility |
|---|---|---|
| `reentry_model/shared/__init__.py` | create | package docstring: the independence rule |
| `reentry_model/shared/{trajectory,earth,atmosphere,fap,aero,gas,constants,heating,material,sesam_io,compare,viz}.py` | `git mv` from `reentry_model/`, path fixes in six | the architecture-independent modules |
| `reentry_model/{same twelve}.py` | create (alias) | old path → the shared module object |
| `reentry_model/shared/body.py` | create | `Body` protocol, `sphere_mass`, `ConstantBody` |
| `reentry_model/body.py` | modify | `ThermalBody` only; re-exports the shared names |
| `reentry_model/shared/thermal.py` | create | `SIGMA_SB`, `MissingBackend`, `StepResult`, `ThermalSolver` |
| `reentry_model/thermal/__init__.py` | modify | `SOLVER_NAMES`, `thermal_solver`; re-exports the shared names |
| `reentry_model/shared/cli_common.py` | create | run defaults and flight arguments, `make_atmosphere`, `model_run_name`, provenance, `ARCHITECTURES`, `architecture_of` |
| `reentry_model/cli.py` | modify | the skyscrapers run driver; `--architecture`; the levelset exit |
| `reentry_model/coupled.py` | modify | imports from `shared` |
| `tests/test_reentry_model_shared.py` | create | aliases, re-exports, independence, paths, imports in tests and analysis |
| `tests/test_reentry_model_cli.py` | modify (append) | the `--architecture` flag |
| `prototype/levelset/stage0/bitid.sh` | create | the three bit-identity flights |
| `README.md`, `CLAUDE.md` | modify | the two architectures and the shared package |

---

### Task 1: Worktree, branch and the bit-identity baseline

**Files:**
- Create: `prototype/levelset/stage0/bitid.sh`

**Interfaces:**
- Produces: `prototype/levelset/stage0/bitid.sh <repo> <outdir>` writing `step1_msis.csv`, `step2_physics.csv`,
  `step2_sesam.csv` (+ JSON) into `<outdir>`; the baseline in `reentry_model_output/stage0_bitid/before/` (git-ignored).

- [ ] **Step 1: Create the worktree on a new branch from `main`**

Use superpowers:using-git-worktrees, or:

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution"
git worktree add .claude/worktrees/levelset-stage0 -b levelset-stage0 main
cd .claude/worktrees/levelset-stage0
git log --oneline -1
```

Expected: the branch starts at `main`'s head; `.claude/worktrees/` is git-ignored. Every later step runs in this worktree.

- [ ] **Step 2: Create `prototype/levelset/stage0/bitid.sh`**

```bash
#!/bin/bash
# Stage 0 bit-identity runs: three flights whose history CSVs must be byte-identical before and after the
# shared-package refactor. Usage: bitid.sh <repo checkout> <output dir>
set -euo pipefail
REPO="$1"; OUT="$2"; PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
mkdir -p "$OUT"; cd "$REPO"
COMMON=(--diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --outdir "$OUT" --quiet)
"$PY" -m reentry_model run "${COMMON[@]}" --name step1_msis                                     # Step 1, whole flight
"$PY" -m reentry_model run "${COMMON[@]}" --atmosphere us76 --thermal fem --heating physics \
    --t-max 20 --dt 0.5 --name step2_physics                                                  # Step 2, physics heating
"$PY" -m reentry_model run "${COMMON[@]}" --atmosphere us76 --thermal fem --heating sesam \
    --t-max 20 --dt 0.5 --name step2_sesam                                                    # Step 2, SESAM-equivalent heating
```

```bash
chmod +x prototype/levelset/stage0/bitid.sh
```

- [ ] **Step 3: Record the baseline twice and confirm it is deterministic**

```bash
prototype/levelset/stage0/bitid.sh "$PWD" "$PWD/reentry_model_output/stage0_bitid/before"
prototype/levelset/stage0/bitid.sh "$PWD" "$PWD/reentry_model_output/stage0_bitid/before2"
for f in step1_msis step2_physics step2_sesam; do cmp reentry_model_output/stage0_bitid/before/$f.csv reentry_model_output/stage0_bitid/before2/$f.csv && echo "$f deterministic"; done
```

Expected: about 50 s per pass; three lines `… deterministic`. (Only the CSVs are compared: the JSONs carry run time and
the git commit.)

- [ ] **Step 4: Record the unit tier's count before the change**

Run: `"$PY" -m pytest -m "not drama and not reference" --collect-only -q | tail -1`
Expected: `352/372 tests collected (20 deselected)`.

- [ ] **Step 5: Commit**

```bash
git add prototype/levelset/stage0/bitid.sh
git commit -m "Level-set stage 0: the bit-identity flights that the shared-package refactor must reproduce byte for byte

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Move the twelve architecture-independent modules into `reentry_model/shared/`

**Files:**
- Create: `reentry_model/shared/__init__.py`, `tests/test_reentry_model_shared.py`
- Move: `reentry_model/{trajectory,earth,atmosphere,fap,aero,gas,constants,heating,material,sesam_io,compare,viz}.py` → `reentry_model/shared/`
- Modify: `reentry_model/shared/{atmosphere,aero,material}.py` (one line each), `reentry_model/shared/{fap,sesam_io,compare}.py` (one line each)
- Create (aliases): `reentry_model/{the same twelve}.py`

**Interfaces:**
- Produces: `reentry_model.shared.<name>` for the twelve names; `reentry_model.<name> is reentry_model.shared.<name>`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_shared.py`:

```python
"""Stage 0 of the level-set design: the shared package. The old module paths are the shared modules themselves,
the split modules re-export their shared names, the shared package imports nothing from either architecture, and the
data and repository paths still resolve after the move one directory down."""
import importlib
import os
import subprocess
import sys

import pytest

from helpers import PY, REPO_ROOT

MOVED = ("trajectory", "earth", "atmosphere", "fap", "aero", "gas", "constants", "heating", "material", "sesam_io",
         "compare", "viz")


@pytest.mark.parametrize("name", MOVED)
def test_old_path_is_the_shared_module(name):
    old = importlib.import_module("reentry_model." + name)
    new = importlib.import_module("reentry_model.shared." + name)
    assert old is new
    assert getattr(importlib.import_module("reentry_model"), name) is new


def test_import_order_does_not_matter():
    code = ("import reentry_model.shared.aero as new\n"
            "import reentry_model.aero as old\n"
            "from reentry_model import aero as attr\n"
            "print(old is new and attr is new)\n")
    out = subprocess.run([PY, "-c", code], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    assert out == "True"


def test_data_and_repository_paths_resolve():
    from reentry_model.shared import aero, atmosphere, compare, fap, material, sesam_io
    for path in (aero.DEFAULT_ATDB, atmosphere.DEFAULT_US76, material.DEFAULT_MATERIAL, fap.DEFAULT_FAP_DAY,
                 sesam_io.REFERENCE_DIR):
        assert os.path.exists(path), path
    for root in (fap.REPO_ROOT, sesam_io.REPO_ROOT, compare.REPO_ROOT):
        assert os.path.samefile(root, REPO_ROOT)
    assert os.path.join(REPO_ROOT, "analysis") in sys.path        # compare.py's import path for the analysis scripts
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py -q`
Expected: FAIL — 14 failed (`ModuleNotFoundError: No module named 'reentry_model.shared'`).

- [ ] **Step 3: Create the package and move the modules**

Create `reentry_model/shared/__init__.py`:

```python
"""Code shared by the model's two architectures, skyscrapers and the level set (design of record:
docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md, section 4).

Modules here depend on nothing in reentry_model outside this package except the package-level constants DATA_DIR and
__version__ (tests/test_reentry_model_shared.py checks it). The old module paths (reentry_model.aero, ...) are aliases
of these modules. A change here must keep both architectures' tests passing."""
```

```bash
for m in trajectory earth atmosphere fap aero gas constants heating material sesam_io compare viz; do
    git mv reentry_model/$m.py reentry_model/shared/$m.py
done
```

The moved modules' relative imports among themselves (`from .constants import …`, `from . import aero, earth`,
`from . import aero, gas`) still resolve, because all their targets moved together.

- [ ] **Step 4: Fix the six paths that are relative to the module's own file**

In `reentry_model/shared/atmosphere.py`, `reentry_model/shared/aero.py` and `reentry_model/shared/material.py`, replace

```python
from . import DATA_DIR
```

with

```python
from .. import DATA_DIR
```

In `reentry_model/shared/fap.py`, `reentry_model/shared/sesam_io.py` and `reentry_model/shared/compare.py`, replace

```python
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
```

with

```python
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
```

- [ ] **Step 5: Write the twelve aliases at the old paths**

```bash
for m in trajectory earth atmosphere fap aero gas constants heating material sesam_io compare viz; do
cat > reentry_model/$m.py <<EOF
"""Moved to reentry_model.shared.$m (stage 0 of the level-set design). This old path is kept as an alias:
importing it returns the shared module object itself, so attributes, private names and identity all match."""
import sys

from .shared import $m as _shared

sys.modules[__name__] = _shared
EOF
done
```

Why an alias and not `from .shared.aero import *` (the spec's example): a star import copies names at import time,
skips private names and makes two module objects, so `monkeypatch.setattr(reentry_model.aero, …)` would not reach the
code. The import system re-reads `sys.modules` after running a module, so the replacement is what every importer gets.

- [ ] **Step 6: Run the new tests and the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py -q`
Expected: PASS — 14 passed.

Run: `"$PY" -m pytest -m "not drama and not reference" -q`
Expected: every test passes except the FEniCSx module's skip (about 70 s).

- [ ] **Step 7: Commit**

```bash
git add -A reentry_model tests/test_reentry_model_shared.py
git commit -m "The twelve architecture-independent modules move to reentry_model/shared, with their old paths kept as aliases of the same module objects

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Split the body seam and the thermal protocol

**Files:**
- Create: `reentry_model/shared/body.py`, `reentry_model/shared/thermal.py`
- Modify: `reentry_model/body.py`, `reentry_model/thermal/__init__.py`, `reentry_model/coupled.py`
- Test: `tests/test_reentry_model_shared.py` (append)

**Interfaces:**
- Consumes: Task 2's `reentry_model.shared`.
- Produces: `reentry_model.shared.body.{Body, sphere_mass, ConstantBody}`;
  `reentry_model.shared.thermal.{SIGMA_SB, MissingBackend, StepResult, ThermalSolver}`; the same names still importable
  from `reentry_model.body` and `reentry_model.thermal` (the backends keep `from . import SIGMA_SB, StepResult`).

- [ ] **Step 1: Write the failing test**

Append to `tests/test_reentry_model_shared.py`:

```python


def test_body_and_thermal_reexport_the_shared_names():
    from reentry_model import body, thermal
    from reentry_model.shared import body as sbody, thermal as sthermal
    assert (body.Body, body.ConstantBody, body.sphere_mass) == (sbody.Body, sbody.ConstantBody, sbody.sphere_mass)
    for name in ("SIGMA_SB", "MissingBackend", "StepResult", "ThermalSolver"):
        assert getattr(thermal, name) is getattr(sthermal, name)
```

- [ ] **Step 2: Run it to verify it fails**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py::test_body_and_thermal_reexport_the_shared_names -q`
Expected: FAIL — `ImportError: cannot import name 'body' from 'reentry_model.shared'`.

- [ ] **Step 3: Create `reentry_model/shared/body.py`**

```python
"""The seam between the trajectory and a body, shared by both architectures (spec section 4 of the trajectory
design): the `Body` protocol every body implements, the solid sphere's mass, and `ConstantBody` (Step 1 behaviour,
`--thermal none`: fixed mass, one temperature), which needs no architecture."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def surface_temperature(self) -> np.ndarray: ...     # K per patch
    def mean_temperature(self) -> float: ...             # K, energy-equivalent (spec 6.4)
    def energy(self) -> float: ...                       # J stored above the material's reference temperature
    def field(self): ...                                 # nodal temperatures, or None


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def advance(self, t, dt, loads):
        return None

    def surface_temperature(self):
        return np.array([self.temperature_K])

    def mean_temperature(self):
        return self.temperature_K

    def energy(self):
        return 0.0

    def field(self):
        return None
```

- [ ] **Step 4: Reduce `reentry_model/body.py` to `ThermalBody`**

Replace everything above `class ThermalBody:` (the old docstring, imports, `Body`, `sphere_mass` and `ConstantBody`)
with:

```python
"""The skyscrapers architecture's thermal body (spec section 4): ThermalBody wraps a mesh, a material and a conduction
solver and advances the temperature field with each macro step's per-patch convective loads, keeping the heat
bookkeeping SESAM reports. The `Body` protocol, `sphere_mass` and `ConstantBody` live in reentry_model.shared.body and
are re-exported here so that their old names keep working."""
import numpy as np

from .shared.body import Body, ConstantBody, sphere_mass  # noqa: F401  (re-exported)


```

`class ThermalBody` and everything after it are unchanged (it uses only `np`).

- [ ] **Step 5: Create `reentry_model/shared/thermal.py`**

```python
"""The Step 2 conduction-solver protocol, shared by both architectures' thermal backends (thermal design section 8):
the Stefan-Boltzmann constant, the error a backend raises when its library is missing, the result of one step and
the protocol itself. The skyscrapers backends and their factory stay in reentry_model.thermal."""
from dataclasses import dataclass
from typing import Protocol

import numpy as np

SIGMA_SB = 5.670374419e-8         # Stefan-Boltzmann [W/(m2 K4)]


class MissingBackend(RuntimeError):
    """The selected backend's library is not importable in this interpreter."""


@dataclass
class StepResult:
    T: np.ndarray                 # nodal temperatures after the step [K]
    Q_conv: float                 # W, convective power applied over the step
    Q_rad: float                  # W, radiated power at the end of the step
    iterations: int               # Newton/Picard iterations taken


class ThermalSolver(Protocol):
    def setup(self, mesh, material, emissivity) -> None: ...
    def set_temperature(self, T) -> None: ...                       # uniform float or nodal array
    def step(self, dt, q_conv, T_amb, dirichlet=None) -> StepResult: ...
    def temperature(self) -> np.ndarray: ...
    def energy(self) -> float: ...                                  # stored enthalpy above material.T_REF [J]
    def radiated_power(self, T_amb) -> float: ...
```

- [ ] **Step 6: Replace `reentry_model/thermal/__init__.py`**

```python
"""Finite-element conduction solvers behind one protocol (spec section 8); `thermal_solver(name)` picks the backend.
The protocol, StepResult, MissingBackend and SIGMA_SB live in reentry_model.shared.thermal and are re-exported here.

Both backends solve rho c_p(T) dT/dt = div(k(T) grad T) with -k grad T . n = eps sigma (T^4 - T_amb^4) - q_conv on the
boundary, backward Euler in time, P1 tetrahedra in space, and report the same quantities: nodal temperatures,
stored enthalpy, radiated power. `dirichlet=(nodes, values)` in `step` is for the analytic tests only."""
from ..shared.thermal import SIGMA_SB, MissingBackend, StepResult, ThermalSolver  # noqa: F401  (re-exported)

SOLVER_NAMES = ("skfem", "fenicsx")


def thermal_solver(name, **options):
    if name == "skfem":
        from .skfem_backend import SkfemThermalSolver
        return SkfemThermalSolver(**options)
    if name == "fenicsx":
        from .fenicsx_backend import FenicsxThermalSolver      # raises MissingBackend without dolfinx
        return FenicsxThermalSolver(**options)
    raise ValueError("thermal solver must be one of {}, got {!r}".format(SOLVER_NAMES, name))
```

- [ ] **Step 7: Point `reentry_model/coupled.py` at `shared`**

Replace `from . import trajectory as tj` with `from .shared import trajectory as tj`, and inside `write_vtk_frame`
replace `    from .thermal import SIGMA_SB` with `    from .shared.thermal import SIGMA_SB`.

- [ ] **Step 8: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py -q`
Expected: PASS — 15 passed.

Run: `"$PY" -m pytest -m "not drama and not reference" -q`
Expected: every test passes except the FEniCSx module's skip.

Run: `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang /Users/ashajain/miniforge3/envs/fenicsx_env/bin/python -m pytest tests/test_reentry_model_fenicsx.py -q`
Expected: `7 passed` (the FEniCSx backend imports `SIGMA_SB, MissingBackend, StepResult` through the re-export).

- [ ] **Step 9: Commit**

```bash
git add reentry_model/shared/body.py reentry_model/shared/thermal.py reentry_model/body.py reentry_model/thermal/__init__.py reentry_model/coupled.py tests/test_reentry_model_shared.py
git commit -m "The body seam (Body, sphere_mass, ConstantBody) and the Step 2 thermal protocol move to reentry_model/shared; ThermalBody and the backends stay with skyscrapers

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: The shared run-command parts and the `--architecture` flag

**Files:**
- Create: `reentry_model/shared/cli_common.py`
- Modify: `reentry_model/cli.py`
- Test: `tests/test_reentry_model_shared.py` (append), `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: Task 2's `shared.{aero, atmosphere, fap, sesam_io, earth}`.
- Produces: `reentry_model.shared.cli_common.{REPO_ROOT, DEFAULT_OUTDIR, DEFAULT_HEADING_DEG, DEFAULT_LAT_DEG,
  DEFAULT_LON_DEG, DEFAULT_EPOCH, DEFAULT_MATERIAL_DENSITY, ATMOSPHERES, WINDS, ARCHITECTURES, DEFAULT_ARCHITECTURE,
  parse_epoch(text), model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name,
  heating_name=None), make_atmosphere(spec, epoch, wind_name), git_commit(), provenance(),
  add_flight_arguments(parser), architecture_of(run_argv) -> str}`; `reentry_model.cli.LEVELSET_NOT_BUILT`. The
  level-set run driver (stage 1) builds its parser with `add_flight_arguments` and is reached from `cli.main` where the
  levelset exit is now.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_shared.py`:

```python


def test_cli_common_reexports_and_paths():
    from reentry_model import cli
    from reentry_model.shared import cli_common
    for name in ("parse_epoch", "model_run_name", "make_atmosphere", "provenance", "git_commit", "DEFAULT_OUTDIR",
                 "REPO_ROOT", "ATMOSPHERES", "WINDS"):
        assert getattr(cli, name) is getattr(cli_common, name)
    assert os.path.samefile(cli_common.REPO_ROOT, REPO_ROOT)
    assert cli_common.DEFAULT_OUTDIR == os.path.join(cli_common.REPO_ROOT, "reentry_model_output")
    assert cli_common.ARCHITECTURES == ("skyscrapers", "levelset") and cli_common.DEFAULT_ARCHITECTURE == "skyscrapers"
```

Append to `tests/test_reentry_model_cli.py`:

```python


# --- the --architecture flag (stage 0 of the level-set design) ---

def test_architecture_skyscrapers_is_the_default_and_changes_nothing(tmp_path):
    plain = BASE + ["--atmosphere", "us76", "--t-max", "20", "--cadence", "5", "--quiet"]
    assert cli.main(plain + ["--outdir", str(tmp_path / "a")]) == 0
    assert cli.main(plain + ["--architecture", "skyscrapers", "--outdir", str(tmp_path / "b")]) == 0
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none")
    assert (tmp_path / "a" / (name + ".csv")).read_bytes() == (tmp_path / "b" / (name + ".csv")).read_bytes()


@pytest.mark.parametrize("flag", [["--architecture", "levelset"], ["--architecture=levelset"]])
def test_architecture_levelset_exits_2_until_built(flag, tmp_path, capsys):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--t-max", "5", "--outdir", str(tmp_path)] + flag)
    assert rc == 2
    assert "levelset is not built yet" in capsys.readouterr().err
    assert not any(tmp_path.iterdir())                       # nothing written


def test_architecture_unknown_exits_2(tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(BASE + ["--architecture", "towers", "--outdir", str(tmp_path)])
    assert exc.value.code == 2


def test_architecture_does_not_reach_compare(tmp_path):
    with pytest.raises(SystemExit) as exc:                   # `compare` has no --architecture: argparse rejects it
        cli.main(["compare", "--model", "m.csv", "--reference", "r.csv", "--architecture", "levelset"])
    assert exc.value.code == 2
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py tests/test_reentry_model_cli.py -q -k "cli_common or architecture"`
Expected: FAIL — `test_cli_common_reexports_and_paths` (no `cli_common`), `test_architecture_skyscrapers_is_the_default_and_changes_nothing`
and both `test_architecture_levelset_exits_2_until_built` (argparse rejects the unknown `--architecture`: `SystemExit: 2`);
`test_architecture_unknown_exits_2` and `test_architecture_does_not_reach_compare` already pass.

- [ ] **Step 3: Create `reentry_model/shared/cli_common.py`**

```python
"""The run command's architecture-independent parts, shared by both architectures' run drivers: the flight's
defaults and arguments, the atmosphere factory, the run name, provenance, and the --architecture pre-parse that routes
`run` to the skyscrapers driver (reentry_model.cli) or the level-set driver (not built yet)."""
import argparse
import importlib.metadata
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from .. import __version__
from . import aero, atmosphere, fap, sesam_io
from .earth import GRAVITY_MODELS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output")
# the reference cases' break-off state (parent satellite, see README "Reference runs")
DEFAULT_HEADING_DEG = 347.168296
DEFAULT_LAT_DEG = 29.546067
DEFAULT_LON_DEG = -82.134333
DEFAULT_EPOCH = "2024-08-01T12:53:07"
DEFAULT_MATERIAL_DENSITY = 2813.0          # drama-AA7075
ATMOSPHERES = ("nrlmsise", "us76")         # plus replay:<path>
WINDS = ("none", "static")
ARCHITECTURES = ("skyscrapers", "levelset")
DEFAULT_ARCHITECTURE = "skyscrapers"


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))


def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None):
    name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
        diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
    return name + ("_fem-" + heating_name if heating_name else "")


def make_atmosphere(spec, epoch, wind_name):
    """(atmosphere object, short name, provenance dict) for an --atmosphere value."""
    wind = atmosphere.NoWind() if wind_name == "none" else atmosphere.StaticProfileWind()
    if spec == "us76":
        return atmosphere.US76TableAtmosphere(wind=wind), "us76", {}
    if spec == "nrlmsise":
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), epoch.date())
        return atmosphere.NRLMSISE00Atmosphere(epoch, solar, wind), "nrlmsise", \
            {"f107": solar.f107, "f107a": solar.f107a, "ap": solar.ap, "fap_day": fap.DEFAULT_FAP_DAY}
    if spec.startswith("replay:"):
        ref = sesam_io.load_reference(spec[len("replay:"):])
        return atmosphere.ReplayAtmosphere(ref, wind), "replay", {"replay_reference": ref.name, "replay_sha256": ref.sha256}
    raise ValueError("--atmosphere must be one of {} or replay:<sesam.csv>, got {!r}".format(ATMOSPHERES, spec))


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def provenance():
    import scipy
    try:
        import pymsis
        pymsis_version = pymsis.__version__
    except ImportError:
        pymsis_version = None
    package_names = {"skfem": "scikit-fem", "gmsh": "gmsh", "pyamg": "pyamg", "pyvista": "pyvista",
                      "cantera": "cantera", "dolfinx": "fenics-dolfinx"}
    versions = {}
    for key, pkg in package_names.items():
        try:
            versions[key] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[key] = None
    return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0], **versions}


def add_flight_arguments(parser):
    """The flight, atmosphere and output arguments every architecture's `run` takes."""
    parser.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    parser.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s], relative to the rotating atmosphere")
    parser.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    parser.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending")
    parser.add_argument("--heading", type=float, default=DEFAULT_HEADING_DEG, help="[deg] clockwise from north (default %(default)s)")
    parser.add_argument("--lat", type=float, default=DEFAULT_LAT_DEG, help="geodetic latitude [deg] (default %(default)s)")
    parser.add_argument("--lon", type=float, default=DEFAULT_LON_DEG, help="longitude [deg] (default %(default)s)")
    parser.add_argument("--epoch", type=parse_epoch, default=parse_epoch(DEFAULT_EPOCH), help="UTC YYYY-MM-DDTHH:MM:SS (default {})".format(DEFAULT_EPOCH))
    parser.add_argument("--material-density", type=float, default=DEFAULT_MATERIAL_DENSITY, help="[kg/m3] (default %(default)s)")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="initial temperature [K]: recorded only with --thermal none, the field's initial condition with --thermal fem (default %(default)s)")
    parser.add_argument("--atmosphere", default="nrlmsise", help="nrlmsise (default) | us76 | replay:<sesam.csv>")
    parser.add_argument("--wind", choices=WINDS, default="none")
    parser.add_argument("--bridging", choices=aero.BRIDGING_NAMES, default="sesam-table")
    parser.add_argument("--gravity", choices=GRAVITY_MODELS, default="j2")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--cadence", type=float, default=1.0, help="history sample spacing [s] (default %(default)s)")
    parser.add_argument("--t-max", type=float, default=3600.0, help="[s] (default %(default)s)")
    parser.add_argument("--reference", default=None,
                        help="SESAM run CSV to compare against (with --thermal none the model is also sampled at its times; "
                             "with --thermal fem the macro-step history is interpolated)")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")
    parser.add_argument("--quiet", action="store_true")


def architecture_of(run_argv):
    """The --architecture value in a `run` command's arguments (unvalidated; the full parser checks the choice)."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--architecture", default=DEFAULT_ARCHITECTURE)
    known, _ = pre.parse_known_args(run_argv)
    return known.architecture
```

- [ ] **Step 4: Rewire `reentry_model/cli.py`**

(a) In the module docstring, after the line
`        [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)` add:

```
        [--architecture skyscrapers|levelset]                          (default skyscrapers; levelset exits 2 until built)
```

(b) Replace everything from `import argparse` down to (not including) `def build_parser():` — the imports, the
constants `REPO_ROOT` … `THERMAL_MODES`, and the functions `parse_epoch`, `model_run_name`, `make_atmosphere`,
`git_commit`, `provenance` — with:

```python
import argparse
import math
import os
import sys

from . import body, coupled, mesh, thermal
from .shared import aero, compare, heating, material, sesam_io, viz
from .shared import trajectory as tj
from .shared.cli_common import (ARCHITECTURES, ATMOSPHERES, DEFAULT_ARCHITECTURE, DEFAULT_EPOCH, DEFAULT_HEADING_DEG,  # noqa: F401
                                DEFAULT_LAT_DEG, DEFAULT_LON_DEG, DEFAULT_MATERIAL_DENSITY, DEFAULT_OUTDIR, REPO_ROOT, WINDS,
                                add_flight_arguments, architecture_of, git_commit, make_atmosphere, model_run_name,
                                parse_epoch, provenance)
from .shared.thermal import MissingBackend

THERMAL_MODES = ("none", "fem")
LEVELSET_NOT_BUILT = ("--architecture levelset is not built yet (stage 1 of "
                      "docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md); use skyscrapers")


```

(c) In `build_parser`, replace the 24 lines from `    r.add_argument("--diameter", …)` through
`    r.add_argument("--quiet", action="store_true")` with:

```python
    add_flight_arguments(r)
    r.add_argument("--architecture", choices=ARCHITECTURES, default=DEFAULT_ARCHITECTURE,
                   help="model architecture: skyscrapers (default; the mesh model of Steps 1-3) | levelset (not built yet)")
```

(d) Replace the first three lines of `main`:

```python
def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
```

with:

```python
def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["run"] and architecture_of(argv[1:]) == "levelset":
        print("ERROR: " + LEVELSET_NOT_BUILT, file=sys.stderr)
        return 2
    parser = build_parser()
    args = parser.parse_args(argv)
```

The rest of `cli.py` (`build_thermal`, `cmd_run`, `cmd_compare`, the error handling in `main`) is unchanged. The
re-exported names (`DEFAULT_*`, `REPO_ROOT`, `WINDS`, `git_commit`, `parse_epoch`) keep `cli.<name>` working for the
tests and scripts that use it.

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py tests/test_reentry_model_cli.py -q`
Expected: PASS — 39 passed (16 + 23).

Run: `"$PY" -m reentry_model run --help | grep -A1 -- --architecture`
Expected: the `--architecture {skyscrapers,levelset}` line and its help.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/shared/cli_common.py reentry_model/cli.py tests/test_reentry_model_shared.py tests/test_reentry_model_cli.py
git commit -m "The run command's flight arguments, defaults and provenance move to reentry_model/shared/cli_common, and run gains --architecture skyscrapers|levelset, levelset exiting 2 until built

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Independence, documentation and the acceptance runs

**Files:**
- Test: `tests/test_reentry_model_shared.py` (append)
- Modify: `README.md`, `CLAUDE.md`

**Interfaces:**
- Consumes: Tasks 2–4; Task 1's baseline and `bitid.sh`.
- Produces: the accepted branch `levelset-stage0`, ready for Asha to merge into `main`.

- [ ] **Step 1: Write the independence and import tests**

Append to `tests/test_reentry_model_shared.py`:

```python


SHARED = MOVED + ("body", "thermal", "cli_common")


def test_shared_package_imports_neither_architecture():
    code = ("import sys\n"
            + "".join("import reentry_model.shared.{}\n".format(n) for n in SHARED)
            + "print('\\n'.join(sorted(m for m in sys.modules if m.startswith('reentry_model'))))\n")
    out = subprocess.run([PY, "-c", code], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.split()
    outside = [m for m in out if m != "reentry_model" and not m.startswith("reentry_model.shared")]
    assert outside == []


def test_every_reentry_model_import_in_tests_and_analysis_resolves():
    import ast
    import glob
    import importlib.util
    files = glob.glob(os.path.join(REPO_ROOT, "tests", "*.py")) + glob.glob(os.path.join(REPO_ROOT, "analysis", "*.py"))
    wanted = set()
    for path in files:
        for node in ast.walk(ast.parse(open(path).read())):
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("reentry_model") and node.level == 0:
                wanted.update((node.module, a.name) for a in node.names)
            elif isinstance(node, ast.Import):
                wanted.update((a.name, None) for a in node.names if a.name.startswith("reentry_model"))
    assert wanted
    for module, name in sorted(wanted, key=str):
        mod = importlib.import_module(module)
        if name is not None:
            assert hasattr(mod, name) or importlib.util.find_spec(module + "." + name), (module, name)
```

- [ ] **Step 2: Run them**

Run: `"$PY" -m pytest tests/test_reentry_model_shared.py -q`
Expected: PASS — 18 passed. (Both pass at once because Tasks 2–4 built what they check; they guard stage 1 onward.
To see the first one bite, add `from .. import mesh` to `reentry_model/shared/gas.py`, rerun, see it fail naming
`reentry_model.mesh`, and revert.)

- [ ] **Step 3: Document the two architectures in `README.md`**

Insert before the line `## Physics model — \`reentry_model\` (Step 1: trajectory)`:

```markdown
## Two model architectures and the shared package

The physics model is built in two architectures that share one command and one package:

- **skyscrapers** (the default): the mesh model of Steps 1–3 below (`reentry_model/`, its specs and plans under
  `docs/superpowers/`).
- **level set**: a level set on a fixed mesh with a flow solve in softened metal and skins on markers, for large
  changes of shape (design: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`; code in
  `reentry_model/levelset/` from its stage 1; plans in `docs/superpowers/plans/levelset/`).

`reentry_model run --architecture skyscrapers|levelset` chooses between them; `skyscrapers` is the default, and
`levelset` exits with code 2 until its run driver exists. Code both use lives in `reentry_model/shared/`: `trajectory`,
`earth`, `atmosphere`, `fap`, `aero`, `gas`, `constants`, `heating`, `material`, `sesam_io`, `compare`, `viz`, the
`Body` protocol and `ConstantBody` (`shared/body.py`), the Step 2 thermal protocol (`shared/thermal.py`) and the run
command's common arguments (`shared/cli_common.py`). The old paths (`reentry_model.aero`, …) are aliases of the shared
modules, so existing imports keep working; new code imports from `reentry_model.shared`. A change to `shared/` must keep
both architectures' tests passing. The skins are never shared between the architectures. Step 3's physics joins
`shared/` through the stage S sub-plans listed in `docs/superpowers/plans/shared/README.md`.
```

- [ ] **Step 4: Update `CLAUDE.md`'s architecture section**

In `CLAUDE.md`, directly under the heading `## Architecture of \`reentry_model\``, insert this paragraph before
`The run is a lockstep coupling, not a monolith.`:

```markdown
Two architectures share this package (since level-set stage 0): **skyscrapers**, described below, and the **level
set** (`docs/superpowers/specs/levelset/`), chosen by `run --architecture`. The supporting modules listed below that
both use — `trajectory`, `earth`, `atmosphere`, `fap`, `aero`, `gas`, `constants`, `heating`, `material`, `sesam_io`,
`compare`, `viz`, the `Body` protocol, the thermal protocol and the common run arguments — live in
`reentry_model/shared/`; their old paths are aliases of the same module objects. Import new code from
`reentry_model.shared`, keep `shared/` free of either architecture's modules, and never share the skins.
```

- [ ] **Step 5: The bit-identity check**

```bash
prototype/levelset/stage0/bitid.sh "$PWD" "$PWD/reentry_model_output/stage0_bitid/after"
for f in step1_msis step2_physics step2_sesam; do cmp reentry_model_output/stage0_bitid/before/$f.csv reentry_model_output/stage0_bitid/after/$f.csv && echo "$f identical"; done
```

Expected: three lines `… identical`. Any difference fails the stage: find it before going on.

- [ ] **Step 6: The full test tiers**

Run: `"$PY" -m pytest -m "not drama and not reference" -q`
Expected: 375 selected, all passing except the FEniCSx module's skip (about 70 s; it was 352 before Task 2: 18 new
tests in `test_reentry_model_shared.py`, 5 in `test_reentry_model_cli.py`).

Run: `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang /Users/ashajain/miniforge3/envs/fenicsx_env/bin/python -m pytest tests/test_reentry_model_fenicsx.py -q`
Expected: `7 passed`.

Run: `"$PY" -m pytest -m reference -q`
Expected: every reference test passes (about 50 min; do not edit the package while it runs).

- [ ] **Step 7: Commit**

```bash
git add tests/test_reentry_model_shared.py README.md CLAUDE.md
git commit -m "Stage 0 accepted: the shared package imports neither architecture, every import in the tests and analysis scripts resolves, and the README and CLAUDE.md describe the two architectures

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 8: Hand back for merging**

Do not merge or push. Report to Asha: the bit-identity result, the three tiers' counts, and that `levelset-stage0` is
ready to merge into `main` before any level-set or stage S code (spec §4.4). After she merges, the skyscrapers branch
(`step3-rigid-substrate-scheil`) needs `main` merged in before its next code task, since its sub-plans now import
from `reentry_model.shared`.
