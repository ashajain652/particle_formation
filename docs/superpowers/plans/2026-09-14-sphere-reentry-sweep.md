# Sphere Fragment Re-entry Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Two scripts — `sphere_reentry.py` runs one solid AA7075 sphere through ESA DRAMA/SARA (SESAM) via pyDRAMA and writes a merged CSV history plus a statistics JSON; `sphere_sweep.py` drives it over a diameter × temperature × velocity grid whose initial states are inherited from a parent-satellite SESAM run, in confirmed batches with a progress bar, resume support and an aggregate summary.

**Architecture:** Script 1 is a single importable module whose pure functions (config builder, parsers, merge, statistics) are unit-tested against real SESAM output fixtures; only `run_sphere()` touches pyDRAMA. Script 2 imports script 1 for naming, maps grid velocities to parent states by interpolation on the parent's descending branch, and launches script 1 as one subprocess per point from a thread pool, tracking everything in `sweep_manifest.json` and `sweep_summary.csv`.

**Tech Stack:** Python 3.12 in conda env `drama_env`, pyDRAMA 4.1.4 (`drama` package) + DRAMA 4.1.4 at `/Applications/DRAMA-4.1.4`, numpy (linspace), tqdm (progress bar), pytest.

**Spec:** `docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md` — read it first; sections are referenced below as "spec §N".

## Global Constraints

- Interpreter for every command: `$PY` = `/Users/ashajain/miniforge3/envs/drama_env/bin/python`. Run every command from the repo root `/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution` (note the space after `Documents` — always quote the path).
- `DRAMA_INSTALL_PATH` is `/Applications/DRAMA-4.1.4`; script 1 sets it with `os.environ.setdefault` before importing `drama`.
- pyDRAMA is called exactly as `sara.run(config=[cfg], save_output_dirs=<dir>, keep_output_files="all", fap_day_content=<lines>, fap_mon_content=<lines>, parallel=False, timeout=<s>, log_level="ERROR", spell_check=False)`. `config` MUST be a one-element list (spec §2).
- Constants (spec §4.2): `RHO_AA7075 = 2813.0`, `T_MELT_AA7075 = 850.0`, `MELT_TOLERANCE_K = 0.5`, `ENERGY_THRESHOLD_J = 1e-9`, `MATERIAL_NAME = "drama-AA7075"`, `CSV_MASS_RESOLUTION_KG = 1e-3`, `COARSE_MASS_LIMIT_KG = 0.05`, `DEMISE_MASS_FRACTION = 0.05`.
- Run name format (spec §4.3): `sphere_d{diameter_mm:06.2f}mm_T{temperature_K:06.1f}K_v{velocity_kms:08.5f}kms_h{altitude_km:07.3f}km`.
- CSV columns (spec §4.7) and summary columns (spec §5.6) exactly as listed there; JSON `schema_version` is `1`.
- Only new dependencies: `tqdm`, `pytest` (installed in Task 1). No pandas/scipy.
- Files: `sphere_reentry.py` and `sphere_sweep.py` at the repo root; tests in `tests/`; fixtures already committed in `tests/fixtures/` (see `tests/fixtures/README.md`). Never commit `sphere_sweep_output/` (git-ignored).
- Tests: unit tests must pass without DRAMA (`$PY -m pytest tests -m "not drama"`); integration tests carry `@pytest.mark.drama`.
- Commit after every task. Commit messages end with the line `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.

---

### Task 1: Environment, data files and test scaffolding

**Files:**
- Create: `pytest.ini`, `tests/helpers.py`, `tests/conftest.py`, `tests/test_fixtures.py`, `data/fap_day.dat`, `data/fap_mon.dat`
- Existing (do not modify): `tests/fixtures/**`

**Interfaces:**
- Produces: `tests/helpers.py` with `REPO_ROOT: str`, `FIXTURES: str`, `PY: str`, `fixture_dir(name) -> str`, `fixture_file(name, suffix) -> str` (the one file in `tests/fixtures/<name>/` ending in `suffix`). Every later test imports these with `from helpers import ...`.
- Produces: pytest marker `drama` (auto-skipped when `drama` is not importable).

- [ ] **Step 1: Install pytest and tqdm into drama_env and set git identity**

Run:
```bash
"$PY" -m pip install pytest tqdm && "$PY" -c "import pytest, tqdm, numpy; print(pytest.__version__, tqdm.__version__, numpy.__version__)"
```
Expected: three version numbers printed, no error.

Run (only sets values if missing):
```bash
git config user.name >/dev/null || git config user.name "Asha Jain"; git config user.email >/dev/null || git config user.email "ashajain652@gmail.com"; git config user.name; git config user.email
```

- [ ] **Step 2: Copy the parent run's space-weather files into data/**

Run:
```bash
mkdir -p data && cp "Generic_Satellite Reentry/data/fap_day.dat" "Generic_Satellite Reentry/data/fap_mon.dat" data/ && head -c 120 data/fap_day.dat && echo && head -c 120 data/fap_mon.dat && echo
```
Expected: both files start with a `#` comment line (`#d/mm/yyyy F10 F3M SSN Ap ...` and `#___Long-Term Activity Forecast ...`).

- [ ] **Step 3: Write pytest.ini**

```ini
[pytest]
testpaths = tests
addopts = -ra
markers =
    drama: needs DRAMA 4.1.4 and pyDRAMA (conda env drama_env); skipped automatically when the drama package is not importable
```

- [ ] **Step 4: Write tests/helpers.py and tests/conftest.py**

`tests/helpers.py`:
```python
"""Paths and helpers shared by the tests (importable because pytest puts tests/ on sys.path)."""
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(REPO_ROOT, "tests", "fixtures")
PY = "/Users/ashajain/miniforge3/envs/drama_env/bin/python"


def fixture_dir(name):
    return os.path.join(FIXTURES, name)


def fixture_file(name, suffix):
    """The single file in tests/fixtures/<name>/ whose name ends with `suffix`."""
    hits = [f for f in os.listdir(fixture_dir(name)) if f.endswith(suffix)]
    assert len(hits) == 1, (name, suffix, hits)
    return os.path.join(fixture_dir(name), hits[0])
```

`tests/conftest.py`:
```python
"""Pytest configuration: repo root on sys.path, 'drama' marker auto-skip."""
import importlib.util
import sys

import pytest

from helpers import REPO_ROOT

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

HAVE_DRAMA = importlib.util.find_spec("drama") is not None


def pytest_collection_modifyitems(config, items):
    if HAVE_DRAMA:
        return
    skip = pytest.mark.skip(reason="pyDRAMA (package 'drama') is not importable in this interpreter")
    for item in items:
        if "drama" in item.keywords:
            item.add_marker(skip)
```

- [ ] **Step 5: Write the fixture sanity tests**

`tests/test_fixtures.py`:
```python
import os

import pytest

from helpers import REPO_ROOT, fixture_file

CASES = ["T1_demised_50mm", "T3_survivor_50mm", "T5_5mm_750K", "E1_ballooning_5mm"]
SUFFIXES = ("_AeroThermalHistory.txt", "_Trajectory.txt", "ImpactingFragments.xml", "sesam.log")


@pytest.mark.parametrize("case", CASES)
def test_fixture_has_the_four_sesam_files(case):
    for suffix in SUFFIXES:
        assert os.path.getsize(fixture_file(case, suffix)) > 0


def test_survivor_fixture_has_full_precision_mass():
    with open(fixture_file("T3_survivor_50mm", "ImpactingFragments.xml")) as fh:
        assert '<mass unit="kg">1.8411041946975187e-01</mass>' in fh.read()


@pytest.mark.parametrize("name", ["fap_day.dat", "fap_mon.dat"])
def test_space_weather_files_present(name):
    with open(os.path.join(REPO_ROOT, "data", name)) as fh:
        assert fh.readline().startswith("#")


@pytest.mark.drama
def test_drama_marker_is_registered():
    import drama  # noqa: F401
```

- [ ] **Step 6: Run the tests**

Run: `"$PY" -m pytest tests/test_fixtures.py -v`
Expected: 8 passed (the `drama` test passes in drama_env; it would be skipped elsewhere).

- [ ] **Step 7: Commit**

```bash
git add pytest.ini tests/helpers.py tests/conftest.py tests/test_fixtures.py data/fap_day.dat data/fap_mon.dat
git commit -m "Scaffold tests, data files and pytest configuration

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: sphere_reentry core — constants, SphereRun, run name, geometry

**Files:**
- Create: `sphere_reentry.py`
- Create: `tests/test_sphere_reentry.py`

**Interfaces:**
- Produces (spec §4.2): module constants `SCRIPT_VERSION`, `SCRIPT_DIR`, `DRAMA_INSTALL_PATH`, `MATERIAL_NAME`, `RHO_AA7075`, `T_MELT_AA7075`, `MELT_TOLERANCE_K`, `ENERGY_THRESHOLD_J`, `CSV_MASS_RESOLUTION_KG`, `COARSE_MASS_LIMIT_KG`, `DEMISE_MASS_FRACTION`, `OBJECT_NAME`, `OBJECT_UUID`, `PARENT_EPOCH`, `DEFAULT_TIMEOUT_S`, `DEFAULT_OUTDIR`, `DEFAULT_RAW_DIR`, `DEFAULT_FAP_DAY`, `DEFAULT_FAP_MON`; `@dataclass(frozen=True) SphereRun(velocity_kms, altitude_km, temperature_K, diameter_mm, flight_path_deg=0.0, heading_deg=0.0, lat_deg=0.0, lon_deg=0.0, epoch=PARENT_EPOCH)` with properties `radius_m`, `mass_kg`, `cross_section_m2`; `sphere_mass_kg(diameter_mm) -> float`; `sphere_cross_section_m2(diameter_mm) -> float`; `run_name(run) -> str`; `read_lines(path) -> list[str]`.

- [ ] **Step 1: Write the failing tests**

`tests/test_sphere_reentry.py`:
```python
"""Unit tests for sphere_reentry.py (no DRAMA needed)."""
import json
import math
from datetime import datetime

import pytest

import sphere_reentry as sr
from helpers import fixture_dir, fixture_file


def make_run(**kw):
    base = dict(velocity_kms=7.5, altitude_km=77.5, temperature_K=300.0, diameter_mm=50.0)
    base.update(kw)
    return sr.SphereRun(**base)


class TestConstants:
    def test_material_and_run_constants(self):
        assert sr.RHO_AA7075 == 2813.0
        assert sr.T_MELT_AA7075 == 850.0
        assert sr.MELT_TOLERANCE_K == 0.5
        assert sr.ENERGY_THRESHOLD_J == 1e-9
        assert sr.MATERIAL_NAME == "drama-AA7075"
        assert sr.PARENT_EPOCH == datetime(2024, 8, 1, 12, 0, 0)
        assert sr.DEFAULT_OUTDIR.endswith("sphere_sweep_output/runs")
        assert sr.DEFAULT_RAW_DIR.endswith("sphere_sweep_output/raw")
        assert sr.DEFAULT_FAP_DAY.endswith("data/fap_day.dat")


class TestGeometry:
    def test_mass_of_5mm_sphere(self):
        assert sr.sphere_mass_kg(5.0) == pytest.approx(1.8411e-4, rel=1e-4)

    def test_mass_of_100mm_sphere(self):
        assert sr.sphere_mass_kg(100.0) == pytest.approx(1.4729, rel=1e-4)

    def test_cross_section_of_50mm_sphere(self):
        assert sr.sphere_cross_section_m2(50.0) == pytest.approx(1.9635e-3, rel=1e-4)

    def test_run_properties(self):
        run = make_run()
        assert run.radius_m == 0.025
        assert run.mass_kg == pytest.approx(0.18411, rel=1e-4)
        assert run.cross_section_m2 == pytest.approx(1.9635e-3, rel=1e-4)
        assert run.epoch == datetime(2024, 8, 1, 12, 0, 0)
        assert run.flight_path_deg == 0.0 and run.heading_deg == 0.0
        assert run.lat_deg == 0.0 and run.lon_deg == 0.0


class TestRunName:
    def test_exact_name(self):
        assert sr.run_name(make_run()) == "sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km"

    def test_zero_padding_of_small_values(self):
        run = make_run(diameter_mm=5.0, temperature_K=750.0, velocity_kms=0.028, altitude_km=0.012)
        assert sr.run_name(run) == "sphere_d005.00mm_T0750.0K_v00.02800kms_h000.012km"

    def test_lexicographic_order_matches_parameter_order(self):
        names = [sr.run_name(make_run(diameter_mm=d, velocity_kms=v))
                 for d in (5.0, 10.0, 100.0) for v in (0.5, 7.5)]
        assert names == sorted(names)

    def test_adjacent_grid_velocities_are_distinct(self):
        step = (7.5 - 0.028) / 99
        assert sr.run_name(make_run(velocity_kms=7.5)) != sr.run_name(make_run(velocity_kms=7.5 - step))


def test_read_lines_drops_blank_lines_and_newlines(tmp_path):
    p = tmp_path / "f.dat"
    p.write_text("# header\n\n01/08/2024 170 170 100 8 3 3 3 3 3 3 3 3\n   \n")
    assert sr.read_lines(str(p)) == ["# header", "01/08/2024 170 170 100 8 3 3 3 3 3 3 3 3"]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'sphere_reentry'`.

- [ ] **Step 3: Write sphere_reentry.py (header, constants, SphereRun, geometry, run name)**

```python
#!/usr/bin/env python3
"""
sphere_reentry.py -- one solid AA7075 sphere through ESA DRAMA / SARA (SESAM)

Runs SESAM through the pyDRAMA package for a single sphere fragment starting
from a geodetic state (altitude, latitude, longitude, velocity, flight-path
angle, heading) and writes

    <outdir>/<run_name>.csv   full-resolution history (trajectory + aerothermal, joined on time)
    <outdir>/<run_name>.json  inputs, statistics, warnings, provenance

Usage
-----
    python sphere_reentry.py --velocity 7.5 --altitude 77.5 --temperature 300 --diameter 50
    python sphere_reentry.py ... --flight-path-angle -0.96 --heading 347.2 --lat 29.5 --lon -82.1 \
                                 --epoch 2024-08-01T12:53:07
    python sphere_reentry.py ... --dry-run          # print the SESAM config and run name, run nothing

Units: velocity km/s, altitude km, temperature K, diameter mm, angles deg.
Requires DRAMA 4.1.4 and pyDRAMA (conda env drama_env); see README.md.
Design: docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import shutil
import socket
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone

SCRIPT_VERSION = "1.0.0"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DRAMA_INSTALL_PATH = "/Applications/DRAMA-4.1.4"

# --- material: drama-AA7075 exactly as shipped in DRAMA 4.1.4 -----------------
MATERIAL_NAME = "drama-AA7075"
RHO_AA7075 = 2813.0            # [kg/m3]
T_MELT_AA7075 = 850.0          # [K]
MELT_TOLERANCE_K = 0.5         # |T - T_melt| <= tol counts as "at melting temperature"

# --- run control ------------------------------------------------------------------
ENERGY_THRESHOLD_J = 1e-9      # SESAM drops fragments below this kinetic energy; parent run used 15 J
CSV_MASS_RESOLUTION_KG = 1e-3  # SESAM prints mass with 3 decimals
COARSE_MASS_LIMIT_KG = 0.05    # warn when the initial mass is below this (CSV mass column too coarse)
DEMISE_MASS_FRACTION = 0.05    # final/initial mass below this -> outcome "demised"
OBJECT_NAME = "sphere"
OBJECT_UUID = "5b3e0000-0000-4000-8000-000000000001"
PARENT_EPOCH = datetime(2024, 8, 1, 12, 0, 0)   # epoch of the parent satellite run
DEFAULT_TIMEOUT_S = 600
DEFAULT_OUTDIR = os.path.join(SCRIPT_DIR, "sphere_sweep_output", "runs")
DEFAULT_RAW_DIR = os.path.join(SCRIPT_DIR, "sphere_sweep_output", "raw")
DEFAULT_FAP_DAY = os.path.join(SCRIPT_DIR, "data", "fap_day.dat")
DEFAULT_FAP_MON = os.path.join(SCRIPT_DIR, "data", "fap_mon.dat")


# =============================================================================
# 1. Inputs
# =============================================================================

def sphere_mass_kg(diameter_mm: float) -> float:
    """Mass of a solid AA7075 sphere: rho * 4/3 * pi * r^3."""
    r = diameter_mm / 2000.0
    return RHO_AA7075 * 4.0 / 3.0 * math.pi * r ** 3


def sphere_cross_section_m2(diameter_mm: float) -> float:
    r = diameter_mm / 2000.0
    return math.pi * r * r


@dataclass(frozen=True)
class SphereRun:
    """All inputs of one sphere run (CLI units)."""
    velocity_kms: float
    altitude_km: float
    temperature_K: float
    diameter_mm: float
    flight_path_deg: float = 0.0
    heading_deg: float = 0.0
    lat_deg: float = 0.0
    lon_deg: float = 0.0
    epoch: datetime = PARENT_EPOCH

    @property
    def radius_m(self) -> float:
        return self.diameter_mm / 2000.0

    @property
    def mass_kg(self) -> float:
        return sphere_mass_kg(self.diameter_mm)

    @property
    def cross_section_m2(self) -> float:
        return sphere_cross_section_m2(self.diameter_mm)


def run_name(run: SphereRun) -> str:
    """Fixed-width name so that lexicographic order equals parameter order (spec 4.3)."""
    return "sphere_d{:06.2f}mm_T{:06.1f}K_v{:08.5f}kms_h{:07.3f}km".format(
        run.diameter_mm, run.temperature_K, run.velocity_kms, run.altitude_km)


def read_lines(path: str) -> list:
    """Non-blank lines of a text file without trailing newlines (DRAMA fap_day/fap_mon content)."""
    with open(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Add sphere_reentry core: constants, SphereRun, run name, sphere geometry

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: build_config — the SESAM configuration

**Files:**
- Modify: `sphere_reentry.py` (append after `read_lines`)
- Modify: `tests/test_sphere_reentry.py` (append)

**Interfaces:**
- Consumes: `SphereRun`, `run_name`, constants from Task 2.
- Produces: `build_config(run: SphereRun) -> dict` — the complete pyDRAMA config of spec §4.4 (`objects` is a one-element list; no `materialList`; `beginDate`/`initialDate` are the `datetime` from the run).

- [ ] **Step 1: Write the failing tests (append to tests/test_sphere_reentry.py)**

```python
class TestBuildConfig:
    def test_geodetic_elements_in_order(self):
        run = make_run(flight_path_deg=-0.959, heading_deg=347.168, lat_deg=29.546, lon_deg=-82.134)
        cfg = sr.build_config(run)
        assert cfg["coordinateSystem"] == "geodetic"
        assert (cfg["element1"], cfg["element2"], cfg["element3"]) == (77.5, 29.546, -82.134)
        assert (cfg["element4"], cfg["element5"], cfg["element6"]) == (7.5, -0.959, 347.168)

    def test_heading_is_normalised_to_0_360(self):
        assert sr.build_config(make_run(heading_deg=-10.0))["element6"] == pytest.approx(350.0)
        assert sr.build_config(make_run(heading_deg=370.0))["element6"] == pytest.approx(10.0)

    def test_object_definition(self):
        cfg = sr.build_config(make_run())
        assert isinstance(cfg["objects"], list) and len(cfg["objects"]) == 1
        obj = cfg["objects"][0]
        assert obj["primitive"] == {"sphere": {"radius": 0.025}}
        assert obj["solid"] is True
        assert obj["material"] == "drama-AA7075"
        assert obj["mass"] == pytest.approx(sr.sphere_mass_kg(50.0))
        assert obj["attitude"] == "tumbling"
        assert obj["quantity"] == 1
        assert obj["name"] == sr.OBJECT_NAME and obj["uniqueID"] == sr.OBJECT_UUID
        assert "wallThickness" not in obj
        assert "materialList" not in cfg

    def test_temperatures_and_epoch(self):
        epoch = datetime(2024, 8, 1, 12, 53, 7)
        cfg = sr.build_config(make_run(temperature_K=450.0, epoch=epoch))
        assert cfg["objects"][0]["temperature"] == 450.0
        assert cfg["globalSpacecraftTemperature"] == 450.0
        assert cfg["beginDate"] == epoch and cfg["initialDate"] == epoch
        assert isinstance(cfg["beginDate"], datetime)

    def test_run_control_settings_mirror_parent_except_energy_threshold(self):
        cfg = sr.build_config(make_run())
        assert cfg["runMode"] == "reentry-only"
        assert cfg["monteCarlo"] is False
        assert cfg["energyThreshold"] == 1e-9
        assert cfg["assumedCrossSection"] == pytest.approx(sr.sphere_cross_section_m2(50.0))
        assert cfg["comment1"] == sr.run_name(make_run())
        expected = {"dragCoefficient": 2.2, "reflectivityCoefficient": 1.3, "attitude": "tumbling",
                    "fragmentsAttitudeAfterBreakup": "inherited", "densityScalingFactor": 1.0,
                    "dynamicEnvironment": True, "useWind": True, "solarActivityFromFile": True,
                    "useEnvironmentCSV": False, "ap": 8, "f107a": 170, "voxelatorMode": 1,
                    "plotVisibilityMaps": False, "plotObjectTrajectories": False,
                    "propagationWithOscar": False, "runID": "SPHERE"}
        for key, value in expected.items():
            assert cfg[key] == value, key

    def test_config_is_json_serialisable_with_default_str(self):
        json.dumps(sr.build_config(make_run()), default=str)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v -k BuildConfig`
Expected: 6 failed with `AttributeError: module 'sphere_reentry' has no attribute 'build_config'`.

- [ ] **Step 3: Implement build_config (append to sphere_reentry.py)**

```python
# =============================================================================
# 2. SESAM configuration (spec 4.4) -- mirrors the parent satellite run except energyThreshold
# =============================================================================

def build_config(run: SphereRun) -> dict:
    """The complete pyDRAMA SARA configuration for one sphere (pass it as config=[cfg])."""
    obj = {
        "name": OBJECT_NAME,
        "uniqueID": OBJECT_UUID,
        "primitive": {"sphere": {"radius": run.radius_m}},
        "mass": run.mass_kg,
        "material": MATERIAL_NAME,
        "solid": True,
        "relativePosition": {"cartX": 0.0, "cartY": 0.0, "cartZ": 0.0,
                             "yaw": 0.0, "pitch": 0.0, "roll": 0.0},
        "scalingFactors": {"drag": 1.0, "lift": 1.0, "sideForce": 1.0,
                           "averageHeatFlux": 1.0, "averageHeatFluxCT": 1.0,
                           "averageHeatFluxTR": 1.0, "averageHeatFluxFM": 1.0},
        "attitude": "tumbling",
        "quantity": 1,
        "temperature": run.temperature_K,
    }
    return {
        # ---- general -------------------------------------------------------
        "runID": "SPHERE",
        "comment1": run_name(run),
        "comment2": "solid AA7075 sphere fragment",
        "beginDate": run.epoch,
        "runMode": "reentry-only",
        "monteCarlo": False,
        # ---- initial state: geodetic = altitude, lat, lon, velocity, flight path, heading
        "coordinateSystem": "geodetic",
        "initialDate": run.epoch,
        "element1": run.altitude_km,
        "element2": run.lat_deg,
        "element3": run.lon_deg,
        "element4": run.velocity_kms,
        "element5": run.flight_path_deg,
        "element6": run.heading_deg % 360.0,
        # ---- spacecraft-level aerodynamics / attitude ----------------------
        "assumedCrossSection": run.cross_section_m2,
        "dragCoefficient": 2.2,
        "reflectivityCoefficient": 1.3,
        "attitude": "tumbling",
        "fragmentsAttitudeAfterBreakup": "inherited",
        "globalSpacecraftTemperature": run.temperature_K,
        # ---- environment (same keys/values as the parent run) ---------------
        "densityScalingFactor": 1.0,
        "dynamicEnvironment": True,
        "useWind": True,
        "solarActivityFromFile": True,
        "useEnvironmentCSV": False,
        "ap": 8,
        "f107a": 170,
        # ---- numerics / output ----------------------------------------------
        "voxelatorMode": 1,
        "energyThreshold": ENERGY_THRESHOLD_J,   # parent used 15 J; see spec section 2
        "plotVisibilityMaps": False,
        "plotObjectTrajectories": False,
        "propagationWithOscar": False,
        # ---- the model (built-in materials.xml is used: no materialList) ----
        "objects": [obj],
    }
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: 16 passed.

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Add SESAM configuration builder for the sphere run

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: SESAM output parsers

**Files:**
- Modify: `sphere_reentry.py` (append)
- Modify: `tests/test_sphere_reentry.py` (append)

**Interfaces:**
- Produces: `AERO_COLUMNS`, `TRAJ_COLUMNS` (lists of SESAM column names); `read_sara_table(path, columns) -> list[dict]`; `parse_sesam_version(path) -> str | None`; `parse_impacting_fragments(path) -> dict | None` with keys `mass_kg, velocity_kms, lat_deg, lon_deg, epoch`; `parse_sesam_log(path) -> dict` with keys `end_of_life_reason` (str, `"unknown"` when absent), `event_end_mass_kg`, `event_start_mass_kg` (float or None); `find_output_files(raw_dir) -> dict` with keys `aero, traj, fragments, log` (path or None).

- [ ] **Step 1: Write the failing tests (append)**

```python
class TestReadSaraTable:
    def test_aero_rows_T1(self):
        rows = sr.read_sara_table(fixture_file("T1_demised_50mm", "_AeroThermalHistory.txt"), sr.AERO_COLUMNS)
        assert len(rows) == 373
        assert list(rows[0].keys()) == sr.AERO_COLUMNS
        first, last = rows[0], rows[-1]
        assert (first["time"], first["altitude"], first["temp"], first["mass"], first["thick"]) == (0.0, 101.247, 300.0, 0.184, 25.0)
        assert (last["time"], last["altitude"], last["temp"], last["mass"], last["thick"]) == (370.55, 80.14, 850.0, 0.0, 0.0)

    def test_traj_rows_T1(self):
        rows = sr.read_sara_table(fixture_file("T1_demised_50mm", "_Trajectory.txt"), sr.TRAJ_COLUMNS)
        assert len(rows) == 373
        first = rows[0]
        assert (first["velocity"], first["path"], first["heading"], first["density"]) == (7.907, -0.18985, 347.95995, 4.493e-07)
        assert (first["lat"], first["lon"]) == (-3.749, -74.64)
        assert rows[-1]["velocity"] == 6.467

    def test_comment_and_malformed_lines_are_skipped(self, tmp_path):
        p = tmp_path / "t.txt"
        p.write_text("# header\n#  Time [s]  x  y\n1.0 2.0 3.0\n4.0 5.0\nabc 1.0 2.0\n\n7.0 8.0 9.0\n")
        rows = sr.read_sara_table(str(p), ["t", "x", "y"])
        assert rows == [{"t": 1.0, "x": 2.0, "y": 3.0}, {"t": 7.0, "x": 8.0, "y": 9.0}]

    def test_empty_file(self, tmp_path):
        p = tmp_path / "e.txt"
        p.write_text("")
        assert sr.read_sara_table(str(p), ["t"]) == []


class TestParseSesamVersion:
    def test_version_from_history_header(self):
        assert sr.parse_sesam_version(fixture_file("T1_demised_50mm", "_AeroThermalHistory.txt")) == "2.3.0"

    def test_no_header_or_missing_file(self, tmp_path):
        p = tmp_path / "x.txt"
        p.write_text("1 2 3\n")
        assert sr.parse_sesam_version(str(p)) is None
        assert sr.parse_sesam_version(str(tmp_path / "nope.txt")) is None


class TestParseImpactingFragments:
    def test_survivor_T3(self):
        frag = sr.parse_impacting_fragments(fixture_file("T3_survivor_50mm", "ImpactingFragments.xml"))
        assert frag["mass_kg"] == 0.18411041946975187
        assert frag["velocity_kms"] == pytest.approx(0.058229426237796787)
        assert frag["lat_deg"] == pytest.approx(36.041459697933099)
        assert frag["lon_deg"] == pytest.approx(-83.913412639736748)
        assert frag["epoch"] == "2024-08-01T13:00:12.080"

    def test_no_fragment_T1(self):
        assert sr.parse_impacting_fragments(fixture_file("T1_demised_50mm", "ImpactingFragments.xml")) is None

    def test_missing_or_malformed_file(self, tmp_path):
        assert sr.parse_impacting_fragments(str(tmp_path / "nope.xml")) is None
        assert sr.parse_impacting_fragments(None) is None
        bad = tmp_path / "bad.xml"
        bad.write_text("<fragments><fragment>")
        assert sr.parse_impacting_fragments(str(bad)) is None


class TestParseSesamLog:
    @pytest.mark.parametrize("case, reason, end_mass, start_mass", [
        ("T1_demised_50mm", "uncritical", 0.0, 0.18411),
        ("T3_survivor_50mm", "ground impact", 0.18411, 0.18411),
        ("T5_5mm_750K", "uncritical", 3e-06, 0.000184),
        ("E1_ballooning_5mm", "ballooning", 3e-06, 0.000184),
    ])
    def test_fixture_logs(self, case, reason, end_mass, start_mass):
        info = sr.parse_sesam_log(fixture_file(case, "sesam.log"))
        assert info["end_of_life_reason"] == reason
        assert info["event_end_mass_kg"] == pytest.approx(end_mass, abs=1e-12)
        assert info["event_start_mass_kg"] == pytest.approx(start_mass, abs=1e-12)

    def test_missing_file(self, tmp_path):
        assert sr.parse_sesam_log(str(tmp_path / "nope.log")) == {
            "end_of_life_reason": "unknown", "event_end_mass_kg": None, "event_start_mass_kg": None}
        assert sr.parse_sesam_log(None)["end_of_life_reason"] == "unknown"


class TestFindOutputFiles:
    def test_finds_files_recursively(self, tmp_path):
        import shutil
        dest = tmp_path / "raw" / "run_0" / "reentry"
        shutil.copytree(fixture_dir("T3_survivor_50mm"), dest)
        files = sr.find_output_files(str(tmp_path / "raw"))
        assert files["aero"].endswith("_AeroThermalHistory.txt")
        assert files["traj"].endswith("_Trajectory.txt")
        assert files["fragments"].endswith("ImpactingFragments.xml")
        assert files["log"].endswith("sesam.log")

    def test_missing_files_are_none(self, tmp_path):
        assert sr.find_output_files(str(tmp_path)) == {"aero": None, "traj": None, "fragments": None, "log": None}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v -k "SaraTable or SesamVersion or ImpactingFragments or SesamLog or FindOutput"`
Expected: failures with `AttributeError: ... has no attribute 'AERO_COLUMNS'` / `'read_sara_table'` etc.

- [ ] **Step 3: Implement the parsers (append to sphere_reentry.py)**

```python
# =============================================================================
# 3. SESAM output parsing (spec 4.6)
# =============================================================================

AERO_COLUMNS = ["time", "altitude", "temp", "mass", "thick", "convectiveHeat",
                "radiativeHeat", "oxidationHeat", "radCooling", "integratedHeat",
                "visibilityFactor"]
TRAJ_COLUMNS = ["time", "altitude", "lat", "lon", "velocity", "downrange", "drag",
                "lift", "side", "knudsen", "mach", "path", "heading", "density",
                "dynamicPressure", "loadFactor"]


def read_sara_table(path: str, columns: list) -> list:
    """Rows of a SESAM whitespace table; '#' lines and rows with the wrong width are skipped."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != len(columns):
                continue
            try:
                rows.append({c: float(p) for c, p in zip(columns, parts)})
            except ValueError:
                continue
    return rows


_SESAM_VERSION_RE = re.compile(r"SESAM\s+([0-9][0-9.]*[0-9])")


def parse_sesam_version(path):
    """SESAM version from the '#  ---- DRAMA ( SESAM 2.3.0 ) ----' header line, or None."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = _SESAM_VERSION_RE.search(line)
            if m:
                return m.group(1)
    return None


def parse_impacting_fragments(path):
    """First <fragment> of PySara.ImpactingFragments.xml (full-precision impact state), or None."""
    if not path or not os.path.isfile(path):
        return None
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    frag = root.find("fragment")
    if frag is None:
        return None

    def number(tag):
        el = frag.find(tag)
        return float(el.text) if el is not None and el.text else None

    velocity_ms = number("velocity")
    epoch = frag.find("epoch")
    return {
        "mass_kg": number("mass"),
        "velocity_kms": velocity_ms / 1000.0 if velocity_ms is not None else None,
        "lat_deg": number("latitude"),
        "lon_deg": number("longitude"),
        "epoch": epoch.text if epoch is not None else None,
    }


_EOL_RE = re.compile(r"reached end of life because (?:of )?(.+?)\.")
_EVENT_END_RE = re.compile(r"EVENT end \S+ ([0-9.eE+-]+)")
_EVENT_START_RE = re.compile(r"EVENT start \S+ ([0-9.eE+-]+)")


def parse_sesam_log(path):
    """End-of-life reason and the 6-decimal EVENT start/end masses from sesam.log."""
    info = {"end_of_life_reason": "unknown", "event_end_mass_kg": None, "event_start_mass_kg": None}
    if not path or not os.path.isfile(path):
        return info
    with open(path, errors="replace") as fh:
        text = fh.read()
    m = _EOL_RE.search(text)
    if m:
        info["end_of_life_reason"] = m.group(1).strip()
    m = _EVENT_END_RE.search(text)
    if m:
        info["event_end_mass_kg"] = float(m.group(1))
    m = _EVENT_START_RE.search(text)
    if m:
        info["event_start_mass_kg"] = float(m.group(1))
    return info


def find_output_files(raw_dir):
    """Locate the four SESAM files anywhere under raw_dir (pyDRAMA writes run_0/reentry/)."""
    def first(pattern):
        hits = sorted(glob.glob(os.path.join(raw_dir, "**", pattern), recursive=True))
        return hits[0] if hits else None
    return {
        "aero": first("*_AeroThermalHistory.txt"),
        "traj": first("*_Trajectory.txt"),
        "fragments": first("*ImpactingFragments.xml"),
        "log": first("sesam.log"),
    }
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: all pass (32 tests).

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Add SESAM output parsers (history tables, impacting fragments, sesam.log)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: merge_histories and write_csv

**Files:**
- Modify: `sphere_reentry.py` (append)
- Modify: `tests/test_sphere_reentry.py` (append)

**Interfaces:**
- Consumes: `AERO_COLUMNS`, `TRAJ_COLUMNS`, `read_sara_table` (Task 4).
- Produces: `CSV_COLUMNS` (the 25 names of spec §4.7, in order); `merge_histories(aero_rows, traj_rows) -> (rows, warnings)` where each row is a dict keyed by `CSV_COLUMNS` (missing cells `None`); `write_csv(path, rows) -> None`.

- [ ] **Step 1: Write the failing tests (append)**

```python
def aero_row(t, temp=300.0, mass=0.184, alt=100.0, thick=25.0):
    return {"time": t, "altitude": alt, "temp": temp, "mass": mass, "thick": thick,
            "convectiveHeat": 1.0, "radiativeHeat": 0.0, "oxidationHeat": 0.0,
            "radCooling": -1.0, "integratedHeat": 10.0, "visibilityFactor": 1.0}


def traj_row(t, v=7.5, alt=100.0, lat=1.0, lon=2.0, downrange=3.0, path=-1.0, heading=350.0):
    return {"time": t, "altitude": alt, "lat": lat, "lon": lon, "velocity": v,
            "downrange": downrange, "drag": 2.2, "lift": 0.0, "side": 0.0, "knudsen": 0.1,
            "mach": 20.0, "path": path, "heading": heading, "density": 1e-6,
            "dynamicPressure": 100.0, "loadFactor": 0.5}


def load_fixture_rows(case):
    aero = sr.read_sara_table(fixture_file(case, "_AeroThermalHistory.txt"), sr.AERO_COLUMNS)
    traj = sr.read_sara_table(fixture_file(case, "_Trajectory.txt"), sr.TRAJ_COLUMNS)
    return aero, traj


class TestMergeHistories:
    def test_csv_columns_are_the_spec_list(self):
        assert sr.CSV_COLUMNS == [
            "time_s", "altitude_km", "velocity_kms", "temperature_K", "mass_kg", "thick_mm",
            "lat_deg", "lon_deg", "downrange_km", "flight_path_deg", "heading_deg",
            "drag", "lift", "side", "knudsen", "mach", "density_kgm3", "dynamic_pressure_Pa",
            "load_factor_g", "convective_heat_W", "radiative_heat_W", "oxidation_heat_W",
            "rad_cooling_W", "integrated_heat_J", "visibility_factor"]

    def test_identical_time_grids_zip_T1(self):
        aero, traj = load_fixture_rows("T1_demised_50mm")
        rows, warnings = sr.merge_histories(aero, traj)
        assert warnings == []
        assert len(rows) == 373
        assert set(rows[0]) == set(sr.CSV_COLUMNS)
        first = rows[0]
        assert (first["time_s"], first["altitude_km"], first["velocity_kms"]) == (0.0, 101.247, 7.907)
        assert (first["temperature_K"], first["mass_kg"], first["thick_mm"]) == (300.0, 0.184, 25.0)
        assert (first["flight_path_deg"], first["heading_deg"], first["density_kgm3"]) == (-0.18985, 347.95995, 4.493e-07)
        assert rows[-1]["velocity_kms"] == 6.467 and rows[-1]["mass_kg"] == 0.0

    def test_mismatched_grids_outer_join_with_warning(self):
        aero = [aero_row(0.0), aero_row(1.0, temp=400.0), aero_row(2.0)]
        traj = [traj_row(0.0), traj_row(2.0, v=7.0), traj_row(3.0, v=6.0)]
        rows, warnings = sr.merge_histories(aero, traj)
        assert warnings == ["time grids differ (aero 3 rows, traj 3 rows)"]
        assert [r["time_s"] for r in rows] == [0.0, 1.0, 2.0, 3.0]
        assert rows[1]["temperature_K"] == 400.0 and rows[1]["velocity_kms"] is None
        assert rows[3]["velocity_kms"] == 6.0 and rows[3]["temperature_K"] is None
        assert rows[1]["altitude_km"] == 100.0  # falls back to the aero altitude


class TestWriteCsv:
    def test_header_and_first_row_T5(self, tmp_path):
        aero, traj = load_fixture_rows("T5_5mm_750K")
        rows, _ = sr.merge_histories(aero, traj)
        path = tmp_path / "out.csv"
        sr.write_csv(str(path), rows)
        lines = path.read_text().splitlines()
        assert lines[0] == ",".join(sr.CSV_COLUMNS)
        assert lines[1].startswith("0.0,77.5,7.5,750.0,0.0,2.5,29.546,-82.134,")
        assert ",2.727e-05," in lines[1]
        assert len(lines) == 1 + 16

    def test_none_is_written_as_empty_field(self, tmp_path):
        rows, _ = sr.merge_histories([aero_row(0.0)], [traj_row(1.0)])
        path = tmp_path / "out.csv"
        sr.write_csv(str(path), rows)
        line = path.read_text().splitlines()[1]
        assert line.startswith("0.0,100.0,,300.0,0.184,25.0,,,")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v -k "Merge or WriteCsv"`
Expected: 5 failed (`AttributeError: ... 'CSV_COLUMNS'` / `'merge_histories'`).

- [ ] **Step 3: Implement merge and CSV writer (append to sphere_reentry.py)**

```python
# =============================================================================
# 4. Merge + CSV (spec 4.6 / 4.7)
# =============================================================================

CSV_COLUMNS = [
    "time_s", "altitude_km", "velocity_kms", "temperature_K", "mass_kg", "thick_mm",
    "lat_deg", "lon_deg", "downrange_km", "flight_path_deg", "heading_deg",
    "drag", "lift", "side", "knudsen", "mach", "density_kgm3", "dynamic_pressure_Pa",
    "load_factor_g", "convective_heat_W", "radiative_heat_W", "oxidation_heat_W",
    "rad_cooling_W", "integrated_heat_J", "visibility_factor",
]
# CSV column -> trajectory-file column ("time_s" / "altitude_km" handled explicitly)
_TRAJ_MAP = {"velocity_kms": "velocity", "lat_deg": "lat", "lon_deg": "lon",
             "downrange_km": "downrange", "flight_path_deg": "path", "heading_deg": "heading",
             "drag": "drag", "lift": "lift", "side": "side", "knudsen": "knudsen", "mach": "mach",
             "density_kgm3": "density", "dynamic_pressure_Pa": "dynamicPressure",
             "load_factor_g": "loadFactor"}
# CSV column -> aerothermal-file column
_AERO_MAP = {"temperature_K": "temp", "mass_kg": "mass", "thick_mm": "thick",
             "convective_heat_W": "convectiveHeat", "radiative_heat_W": "radiativeHeat",
             "oxidation_heat_W": "oxidationHeat", "rad_cooling_W": "radCooling",
             "integrated_heat_J": "integratedHeat", "visibility_factor": "visibilityFactor"}


def _merged_row(t, aero, traj):
    row = {c: None for c in CSV_COLUMNS}
    row["time_s"] = t
    if traj is not None:
        row["altitude_km"] = traj["altitude"]
        for col, src in _TRAJ_MAP.items():
            row[col] = traj[src]
    if aero is not None:
        if row["altitude_km"] is None:
            row["altitude_km"] = aero["altitude"]
        for col, src in _AERO_MAP.items():
            row[col] = aero[src]
    return row


def merge_histories(aero_rows, traj_rows):
    """Join the two SESAM histories on time. Identical grids zip; otherwise outer-join + warning."""
    warnings = []
    a_times = [r["time"] for r in aero_rows]
    t_times = [r["time"] for r in traj_rows]
    if a_times == t_times:
        rows = [_merged_row(t, a, tr) for t, a, tr in zip(a_times, aero_rows, traj_rows)]
    else:
        warnings.append("time grids differ (aero {} rows, traj {} rows)".format(len(aero_rows), len(traj_rows)))
        by_aero = {r["time"]: r for r in aero_rows}
        by_traj = {r["time"]: r for r in traj_rows}
        rows = [_merged_row(t, by_aero.get(t), by_traj.get(t)) for t in sorted(set(by_aero) | set(by_traj))]
    return rows, warnings


def write_csv(path, rows):
    """One row per time step; floats via str() (shortest round-trip repr), None as empty."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row[c]) for c in CSV_COLUMNS})
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: all pass (37 tests).

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Merge SESAM histories on time and write the per-run CSV

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: compute_stats — the five statistics and friends

**Files:**
- Modify: `sphere_reentry.py` (append)
- Modify: `tests/test_sphere_reentry.py` (append)

**Interfaces:**
- Consumes: merged rows (Task 5), `parse_impacting_fragments` / `parse_sesam_log` results (Task 4), `SphereRun` (Task 2).
- Produces: `radius_mm_from_mass(mass_kg) -> float`; `classify_outcome(reason, final_mass_kg, initial_mass_kg) -> str`; `compute_stats(run, rows, fragments, log_info) -> (results: dict, warnings: list[str])` with exactly the result keys of spec §4.8: `max_temperature_K, time_of_max_temperature_s, altitude_of_max_temperature_km, melting_temperature_K, melt_tolerance_K, time_at_melting_temperature_s, n_rows_at_melt, first_time_at_melt_s, last_time_at_melt_s, altitude_first_melt_km, altitude_last_melt_km, initial_mass_kg, final_mass_kg, final_mass_source, mass_loss_fraction, final_radius_mm, final_thick_mm, final_velocity_kms, final_velocity_source, final_time_s, final_altitude_km, final_latitude_deg, final_longitude_deg, downrange_km, end_of_life_reason, outcome, n_rows`.

- [ ] **Step 1: Write the failing tests (append)**

```python
FIXTURE_RUNS = {
    "T1_demised_50mm": dict(velocity_kms=7.907, altitude_km=101.247, temperature_K=300.0, diameter_mm=50.0),
    "T3_survivor_50mm": dict(velocity_kms=0.5, altitude_km=39.94, temperature_K=300.0, diameter_mm=50.0),
    "T5_5mm_750K": dict(velocity_kms=7.5, altitude_km=77.5, temperature_K=750.0, diameter_mm=5.0),
    "E1_ballooning_5mm": dict(velocity_kms=7.5, altitude_km=77.5, temperature_K=750.0, diameter_mm=5.0),
}


def stats_for(case):
    aero, traj = load_fixture_rows(case)
    rows, _ = sr.merge_histories(aero, traj)
    fragments = sr.parse_impacting_fragments(fixture_file(case, "ImpactingFragments.xml"))
    log_info = sr.parse_sesam_log(fixture_file(case, "sesam.log"))
    return sr.compute_stats(make_run(**FIXTURE_RUNS[case]), rows, fragments, log_info)


def synthetic_stats(temps, dt=1.0, mass=0.184, fragments=None, reason="ground impact", event_end=None):
    aero = [aero_row(i * dt, temp=T, mass=mass) for i, T in enumerate(temps)]
    traj = [traj_row(i * dt) for i in range(len(temps))]
    rows, _ = sr.merge_histories(aero, traj)
    log_info = {"end_of_life_reason": reason, "event_end_mass_kg": event_end, "event_start_mass_kg": None}
    return sr.compute_stats(make_run(), rows, fragments, log_info)


class TestMaxTemperature:
    def test_T1_first_occurrence_of_maximum(self):
        r, _ = stats_for("T1_demised_50mm")
        assert (r["max_temperature_K"], r["time_of_max_temperature_s"], r["altitude_of_max_temperature_km"]) == (850.0, 337.55, 84.46)

    def test_T3_survivor(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert (r["max_temperature_K"], r["time_of_max_temperature_s"], r["altitude_of_max_temperature_km"]) == (304.178, 86.08, 15.424)


class TestMeltDuration:
    def test_T1_sums_the_1s_steps_at_850K(self):
        r, _ = stats_for("T1_demised_50mm")
        assert r["time_at_melting_temperature_s"] == pytest.approx(33.0)
        assert r["n_rows_at_melt"] == 34
        assert (r["first_time_at_melt_s"], r["last_time_at_melt_s"]) == (337.55, 370.55)
        assert (r["altitude_first_melt_km"], r["altitude_last_melt_km"]) == (84.46, 80.14)
        assert (r["melting_temperature_K"], r["melt_tolerance_K"]) == (850.0, 0.5)

    def test_T3_never_melts(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert r["time_at_melting_temperature_s"] == 0.0 and r["n_rows_at_melt"] == 0
        assert r["first_time_at_melt_s"] is None and r["altitude_last_melt_km"] is None

    def test_T5_variable_steps(self):
        r, _ = stats_for("T5_5mm_750K")
        assert r["time_at_melting_temperature_s"] == pytest.approx(6.774)
        assert r["n_rows_at_melt"] == 13

    def test_remelt_sums_both_intervals(self):
        r, _ = synthetic_stats([300.0, 850.0, 850.0, 850.0, 700.0, 850.0, 850.0])
        assert r["time_at_melting_temperature_s"] == pytest.approx(3.0)
        assert (r["first_time_at_melt_s"], r["last_time_at_melt_s"]) == (1.0, 6.0)

    def test_tolerance_boundary(self):
        r, _ = synthetic_stats([849.4, 849.6])
        assert r["time_at_melting_temperature_s"] == 0.0 and r["n_rows_at_melt"] == 1
        r, _ = synthetic_stats([849.6, 850.4])
        assert r["time_at_melting_temperature_s"] == pytest.approx(1.0) and r["n_rows_at_melt"] == 2

    def test_single_row_at_melt_has_zero_duration(self):
        r, _ = synthetic_stats([850.0])
        assert r["time_at_melting_temperature_s"] == 0.0 and r["n_rows_at_melt"] == 1


class TestFinalMassAndRadius:
    def test_T3_from_xml(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert r["final_mass_kg"] == 0.18411041946975187
        assert r["final_mass_source"] == "impacting_fragments_xml"
        assert r["final_radius_mm"] == pytest.approx(25.0, abs=1e-6)
        assert r["mass_loss_fraction"] == pytest.approx(0.0, abs=1e-9)
        assert r["final_thick_mm"] == 25.0
        assert r["initial_mass_kg"] == pytest.approx(0.18411041946975187)

    def test_T1_from_log(self):
        r, _ = stats_for("T1_demised_50mm")
        assert r["final_mass_kg"] == 0.0 and r["final_mass_source"] == "sesam_log_event_end"
        assert r["final_radius_mm"] == 0.0 and r["mass_loss_fraction"] == 1.0
        assert r["final_thick_mm"] == 0.0

    def test_T5_residual_from_log(self):
        r, _ = stats_for("T5_5mm_750K")
        assert r["final_mass_kg"] == pytest.approx(3e-6)
        assert r["final_mass_source"] == "sesam_log_event_end"
        assert r["final_radius_mm"] == pytest.approx(0.633, abs=0.005)

    def test_source_priority(self):
        r, _ = synthetic_stats([300.0], mass=0.123, fragments={"mass_kg": 0.1, "velocity_kms": None}, event_end=0.11)
        assert (r["final_mass_kg"], r["final_mass_source"]) == (0.1, "impacting_fragments_xml")
        r, _ = synthetic_stats([300.0], mass=0.123, fragments=None, event_end=0.11)
        assert (r["final_mass_kg"], r["final_mass_source"]) == (0.11, "sesam_log_event_end")
        r, _ = synthetic_stats([300.0], mass=0.123, fragments=None, event_end=None)
        assert (r["final_mass_kg"], r["final_mass_source"]) == (0.123, "history_file")

    def test_radius_from_mass_round_trips(self):
        assert sr.radius_mm_from_mass(sr.sphere_mass_kg(37.0)) == pytest.approx(18.5)
        assert sr.radius_mm_from_mass(0.0) == 0.0


class TestFinalVelocityAndTrajectoryEnd:
    def test_T3_velocity_from_xml(self):
        r, _ = stats_for("T3_survivor_50mm")
        assert r["final_velocity_kms"] == pytest.approx(0.058229426237796787)
        assert r["final_velocity_source"] == "impacting_fragments_xml"
        assert r["final_altitude_km"] == 0.0 and r["final_time_s"] == 261.08

    def test_T1_velocity_from_trajectory(self):
        r, _ = stats_for("T1_demised_50mm")
        assert (r["final_velocity_kms"], r["final_velocity_source"]) == (6.467, "trajectory_file")
        assert (r["final_time_s"], r["final_altitude_km"]) == (370.55, 80.14)
        assert (r["final_latitude_deg"], r["final_longitude_deg"], r["downrange_km"]) == (21.297, -80.098, 2850.943)
        assert r["n_rows"] == 373


class TestOutcome:
    def test_fixture_outcomes(self):
        assert stats_for("T3_survivor_50mm")[0]["outcome"] == "survived"
        assert stats_for("E1_ballooning_5mm")[0]["outcome"] == "demised"
        assert stats_for("T5_5mm_750K")[0]["outcome"] == "demised"      # uncritical, 1.6 % residual
        assert stats_for("T1_demised_50mm")[0]["outcome"] == "demised"
        assert stats_for("E1_ballooning_5mm")[0]["end_of_life_reason"] == "ballooning"

    def test_classify_outcome_rules(self):
        assert sr.classify_outcome("ground impact", 0.0, 0.1) == "survived"
        assert sr.classify_outcome("ballooning", 0.09, 0.1) == "demised"
        assert sr.classify_outcome("uncritical", 0.184, 0.184) == "other"
        assert sr.classify_outcome("uncritical", 0.004, 0.1) == "demised"
        assert sr.classify_outcome("unknown", None, 0.1) == "other"


class TestWarnings:
    def test_coarse_csv_mass_warning_only_for_small_spheres(self):
        assert any("coarse" in w for w in stats_for("T5_5mm_750K")[1])
        assert not any("coarse" in w for w in stats_for("T1_demised_50mm")[1])

    def test_survivor_without_xml_and_unknown_reason_and_single_row(self):
        _, w = synthetic_stats([300.0, 300.0], fragments=None, reason="ground impact")
        assert any("ImpactingFragments" in x for x in w)
        _, w = synthetic_stats([300.0, 300.0], reason="unknown")
        assert any("end-of-life reason" in x for x in w)
        _, w = synthetic_stats([300.0])
        assert any("single row" in x for x in w)
        _, w = stats_for("T1_demised_50mm")
        assert not any("ImpactingFragments" in x for x in w)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v -k "MaxTemperature or MeltDuration or FinalMass or FinalVelocity or Outcome or Warnings"`
Expected: failures with `AttributeError: ... 'compute_stats'`.

- [ ] **Step 3: Implement the statistics (append to sphere_reentry.py)**

```python
# =============================================================================
# 5. Statistics (spec 4.8)
# =============================================================================

def radius_mm_from_mass(mass_kg):
    """Radius of a solid AA7075 sphere of the given mass; 0 for zero mass."""
    if mass_kg is None or mass_kg <= 0.0:
        return 0.0
    return 1000.0 * (3.0 * mass_kg / (4.0 * math.pi * RHO_AA7075)) ** (1.0 / 3.0)


def classify_outcome(reason, final_mass_kg, initial_mass_kg):
    if reason == "ground impact":
        return "survived"
    if reason == "ballooning":
        return "demised"
    if (final_mass_kg is not None and initial_mass_kg > 0.0
            and final_mass_kg < DEMISE_MASS_FRACTION * initial_mass_kg):
        return "demised"
    return "other"


def _at_melt(temp):
    return temp is not None and abs(temp - T_MELT_AA7075) <= MELT_TOLERANCE_K


def _last_with(rows, key):
    for row in reversed(rows):
        if row.get(key) is not None:
            return row
    return None


def compute_stats(run, rows, fragments, log_info):
    """The statistics block of the JSON plus a list of warning strings."""
    results = {}
    warnings = []

    # -- maximum temperature (first occurrence) ------------------------------------
    with_temp = [r for r in rows if r["temperature_K"] is not None]
    if with_temp:
        hottest = max(with_temp, key=lambda r: r["temperature_K"])   # max() keeps the first maximum
        results["max_temperature_K"] = hottest["temperature_K"]
        results["time_of_max_temperature_s"] = hottest["time_s"]
        results["altitude_of_max_temperature_km"] = hottest["altitude_km"]
    else:
        results["max_temperature_K"] = None
        results["time_of_max_temperature_s"] = None
        results["altitude_of_max_temperature_km"] = None

    # -- time at melting temperature -------------------------------------------------
    melt_rows = [r for r in rows if _at_melt(r["temperature_K"])]
    duration = 0.0
    for prev, cur in zip(rows, rows[1:]):
        if _at_melt(prev["temperature_K"]) and _at_melt(cur["temperature_K"]):
            duration += cur["time_s"] - prev["time_s"]
    results["melting_temperature_K"] = T_MELT_AA7075
    results["melt_tolerance_K"] = MELT_TOLERANCE_K
    results["time_at_melting_temperature_s"] = duration
    results["n_rows_at_melt"] = len(melt_rows)
    results["first_time_at_melt_s"] = melt_rows[0]["time_s"] if melt_rows else None
    results["last_time_at_melt_s"] = melt_rows[-1]["time_s"] if melt_rows else None
    results["altitude_first_melt_km"] = melt_rows[0]["altitude_km"] if melt_rows else None
    results["altitude_last_melt_km"] = melt_rows[-1]["altitude_km"] if melt_rows else None

    # -- final mass: XML -> log -> history file ---------------------------------------
    initial_mass = run.mass_kg
    last_mass_row = _last_with(rows, "mass_kg")
    if fragments is not None and fragments.get("mass_kg") is not None:
        final_mass, source = fragments["mass_kg"], "impacting_fragments_xml"
    elif log_info.get("event_end_mass_kg") is not None:
        final_mass, source = log_info["event_end_mass_kg"], "sesam_log_event_end"
    elif last_mass_row is not None:
        final_mass, source = last_mass_row["mass_kg"], "history_file"
    else:
        final_mass, source = None, None
    results["initial_mass_kg"] = initial_mass
    results["final_mass_kg"] = final_mass
    results["final_mass_source"] = source
    results["mass_loss_fraction"] = (1.0 - final_mass / initial_mass) if final_mass is not None else None
    results["final_radius_mm"] = radius_mm_from_mass(final_mass) if final_mass is not None else None
    last_thick_row = _last_with(rows, "thick_mm")
    results["final_thick_mm"] = last_thick_row["thick_mm"] if last_thick_row else None

    # -- final velocity and trajectory end -------------------------------------------
    last_traj = _last_with(rows, "velocity_kms")
    if fragments is not None and fragments.get("velocity_kms") is not None:
        results["final_velocity_kms"] = fragments["velocity_kms"]
        results["final_velocity_source"] = "impacting_fragments_xml"
    elif last_traj is not None:
        results["final_velocity_kms"] = last_traj["velocity_kms"]
        results["final_velocity_source"] = "trajectory_file"
    else:
        results["final_velocity_kms"] = None
        results["final_velocity_source"] = None
    last = rows[-1] if rows else None
    results["final_time_s"] = last["time_s"] if last else None
    results["final_altitude_km"] = last["altitude_km"] if last else None
    results["final_latitude_deg"] = last_traj["lat_deg"] if last_traj else None
    results["final_longitude_deg"] = last_traj["lon_deg"] if last_traj else None
    results["downrange_km"] = last_traj["downrange_km"] if last_traj else None

    # -- end of life --------------------------------------------------------------------
    reason = log_info.get("end_of_life_reason", "unknown")
    results["end_of_life_reason"] = reason
    results["outcome"] = classify_outcome(reason, final_mass, initial_mass)
    results["n_rows"] = len(rows)

    # -- warnings -----------------------------------------------------------------------
    if initial_mass < COARSE_MASS_LIMIT_KG:
        warnings.append("CSV mass column ({:g} kg resolution) is coarse for initial mass {:.6g} kg".format(
            CSV_MASS_RESOLUTION_KG, initial_mass))
    if len(rows) == 1:
        warnings.append("history has a single row")
    if reason == "ground impact" and fragments is None:
        warnings.append("ground impact reported but no ImpactingFragments.xml entry")
    if reason == "unknown":
        warnings.append("end-of-life reason not found in sesam.log")
    return results, warnings
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: all pass (56 tests).

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Compute per-run statistics: max temperature, melt duration, final mass/radius/velocity, outcome

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: run_sphere — pyDRAMA execution and the JSON document

**Files:**
- Modify: `sphere_reentry.py` (append)
- Modify: `tests/test_sphere_reentry.py` (append)

**Interfaces:**
- Consumes: everything above.
- Produces: `class DramaNotAvailable(RuntimeError)`; `run_sphere(run, outdir, raw_dir, timeout, keep_raw, fap_day_lines, fap_mon_lines) -> dict` — writes `<outdir>/<run_name>.csv` and `.json`, returns the JSON document of spec §4.9 (`status` ∈ {`ok`, `error`, `timeout`}); raises `DramaNotAvailable` only when `drama` cannot be imported.
- The tests stub pyDRAMA by inserting fake `drama` / `drama.sara` modules into `sys.modules`; the fake `sara.run` copies a fixture directory to `<save_output_dirs>/run_0/reentry/` exactly like pyDRAMA does.

- [ ] **Step 1: Write the failing tests (append)**

```python
import os
import shutil
import sys
import types


def install_fake_drama(monkeypatch, run_impl, version="stub-4.1.4"):
    drama = types.ModuleType("drama")
    drama.__version__ = version
    sara = types.ModuleType("drama.sara")
    sara.run = run_impl
    drama.sara = sara
    monkeypatch.setitem(sys.modules, "drama", drama)
    monkeypatch.setitem(sys.modules, "drama.sara", sara)


def fake_sara_run(case=None, mode="ok", calls=None):
    """A stand-in for drama.sara.run: copies fixture `case` where pyDRAMA would write its output."""
    def run(config, save_output_dirs, keep_output_files, fap_day_content, fap_mon_content,
            parallel, timeout, log_level, spell_check):
        if calls is not None:
            calls.append(dict(config=config, save_output_dirs=save_output_dirs,
                              keep_output_files=keep_output_files, fap_day_content=fap_day_content,
                              fap_mon_content=fap_mon_content, parallel=parallel, timeout=timeout,
                              log_level=log_level, spell_check=spell_check))
        dest = os.path.join(save_output_dirs, "run_0", "reentry")
        if mode == "raise":
            raise RuntimeError("boom")
        if mode == "error":
            os.makedirs(dest)
            return {"config": config, "errors": [{"status": "error in reentry", "reentry_logfile": "line1\nline2\nfatal"}], "results": []}
        if mode == "timeout":
            os.makedirs(dest)
            return {"config": config, "errors": [{"status": "error: timeout in reentry (600s)", "reentry_logfile": ""}], "results": []}
        if mode == "empty":
            os.makedirs(dest)
            return {"config": config, "errors": [], "results": []}
        if mode == "nofiles":
            os.makedirs(dest)
            return {"config": config, "errors": [], "results": [{"status": "success", "reentry_logfile": "no files"}]}
        shutil.copytree(fixture_dir(case), dest)
        return {"config": config, "errors": [], "results": [{"status": "success", "reentry_logfile": "", "config": {"output_dir": dest}}]}
    return run


T3_RUN = dict(velocity_kms=0.5, altitude_km=39.94, temperature_K=300.0, diameter_mm=50.0,
              flight_path_deg=-32.79047, heading_deg=347.12153, lat_deg=35.911, lon_deg=-83.878,
              epoch=datetime(2024, 8, 1, 12, 55, 51))


class TestRunSphereOk:
    def test_writes_csv_and_json_and_deletes_raw_tree(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm"))
        run = make_run(**T3_RUN)
        outdir, raw = tmp_path / "runs", tmp_path / "raw"
        doc = sr.run_sphere(run, str(outdir), str(raw), 600, False, ["# fap day"], ["# fap mon"])
        name = sr.run_name(run)
        assert doc["status"] == "ok" and doc["error"] is None and doc["run_name"] == name
        assert doc["schema_version"] == 1
        assert (outdir / (name + ".csv")).is_file() and (outdir / (name + ".json")).is_file()
        assert not (raw / name).exists() and doc["files"]["raw_dir"] is None
        assert doc["files"]["csv"] == os.path.abspath(str(outdir / (name + ".csv")))
        r = doc["results"]
        assert r["outcome"] == "survived" and r["final_mass_source"] == "impacting_fragments_xml"
        assert doc["inputs"]["initial_mass_kg"] == pytest.approx(0.18411, rel=1e-4)
        assert doc["inputs"]["epoch_utc"] == "2024-08-01T12:55:51"
        assert doc["inputs"]["sesam_settings"]["energyThreshold"] == 1e-9
        assert "objects" not in doc["inputs"]["sesam_settings"]
        assert doc["inputs"]["material"] == "drama-AA7075" and doc["inputs"]["melting_temperature_K"] == 850.0
        assert doc["provenance"]["sesam_version"] == "2.3.0"
        assert doc["provenance"]["pydrama_version"] == "stub-4.1.4"
        assert doc["provenance"]["script_version"] == sr.SCRIPT_VERSION
        assert doc["provenance"]["wall_time_s"] >= 0.0 and doc["provenance"]["created_utc"].endswith("Z")
        on_disk = json.load(open(outdir / (name + ".json")))
        assert on_disk["results"]["final_mass_kg"] == r["final_mass_kg"] and on_disk["status"] == "ok"
        with open(outdir / (name + ".csv")) as fh:
            header, first = fh.readline().strip(), fh.readline().strip().split(",")
        assert header == ",".join(sr.CSV_COLUMNS)
        assert first[:4] == ["0.0", "39.94", "0.5", "300.0"]

    def test_keep_raw(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run("T1_demised_50mm"))
        run = make_run(velocity_kms=7.907, altitude_km=101.247)
        raw = tmp_path / "raw"
        (raw / sr.run_name(run)).mkdir(parents=True)
        (raw / sr.run_name(run) / "stale.txt").write_text("old")
        doc = sr.run_sphere(run, str(tmp_path / "runs"), str(raw), 600, True, [], [])
        assert doc["status"] == "ok"
        assert doc["files"]["raw_dir"] == str(raw / sr.run_name(run))
        assert (raw / sr.run_name(run) / "run_0" / "reentry" / "sesam.log").is_file()
        assert not (raw / sr.run_name(run) / "stale.txt").exists()
        assert doc["results"]["outcome"] == "demised"

    def test_pydrama_call_arguments(self, tmp_path, monkeypatch):
        calls = []
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm", calls=calls))
        sr.run_sphere(make_run(**T3_RUN), str(tmp_path / "runs"), str(tmp_path / "raw"), 123, False, ["d1", "d2"], ["m1"])
        call = calls[0]
        assert isinstance(call["config"], list) and len(call["config"]) == 1
        assert call["config"][0]["element1"] == 39.94 and call["config"][0]["objects"][0]["solid"] is True
        assert call["save_output_dirs"] == str(tmp_path / "raw" / sr.run_name(make_run(**T3_RUN)))
        assert call["keep_output_files"] == "all" and call["parallel"] is False
        assert call["fap_day_content"] == ["d1", "d2"] and call["fap_mon_content"] == ["m1"]
        assert call["timeout"] == 123 and call["log_level"] == "ERROR" and call["spell_check"] is False
        assert os.environ.get("DRAMA_INSTALL_PATH")


class TestRunSphereErrors:
    @pytest.mark.parametrize("mode, status, needle", [
        ("error", "error", "error in reentry"),
        ("timeout", "timeout", "timeout in reentry"),
        ("empty", "error", "no result"),
        ("raise", "error", "pyDRAMA raised"),
        ("nofiles", "error", "history files not found"),
    ])
    def test_failure_modes_write_error_json_and_keep_raw(self, tmp_path, monkeypatch, mode, status, needle):
        install_fake_drama(monkeypatch, fake_sara_run(mode=mode))
        run = make_run(**T3_RUN)
        doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
        assert doc["status"] == status and needle in doc["error"]
        assert doc["results"] is None
        assert doc["files"]["raw_dir"] == str(tmp_path / "raw" / sr.run_name(run))
        on_disk = json.load(open(tmp_path / "runs" / (sr.run_name(run) + ".json")))
        assert on_disk["status"] == status
        assert not (tmp_path / "runs" / (sr.run_name(run) + ".csv")).exists()

    def test_error_includes_log_tail(self, tmp_path, monkeypatch):
        install_fake_drama(monkeypatch, fake_sara_run(mode="error"))
        doc = sr.run_sphere(make_run(**T3_RUN), str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
        assert "fatal" in doc["error"]

    def test_pydrama_not_importable(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "drama", None)
        with pytest.raises(sr.DramaNotAvailable):
            sr.run_sphere(make_run(**T3_RUN), str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, [], [])
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v -k RunSphere`
Expected: failures with `AttributeError: ... 'run_sphere'` / `'DramaNotAvailable'`.

- [ ] **Step 3: Implement run_sphere (append to sphere_reentry.py)**

```python
# =============================================================================
# 6. Running SESAM through pyDRAMA (spec 4.5 / 4.9 / 4.10)
# =============================================================================

class DramaNotAvailable(RuntimeError):
    """pyDRAMA (package 'drama') cannot be imported."""


def _import_sara():
    os.environ.setdefault("DRAMA_INSTALL_PATH", DRAMA_INSTALL_PATH)
    try:
        import drama
        from drama import sara
    except ImportError as exc:
        raise DramaNotAvailable(
            "pyDRAMA not importable ({}). Install it into this interpreter with:\n"
            "    pip install '{}/TOOLS/drama_python_package'".format(exc, DRAMA_INSTALL_PATH))
    return sara, getattr(drama, "__version__", None)


def _log_tail(text, n=40):
    if text is None:
        return ""
    if isinstance(text, bytes):
        text = text.decode(errors="replace")
    return "\n".join(str(text).splitlines()[-n:])


def _base_document(run, cfg, name, csv_path, run_raw, pydrama_version):
    settings = {k: v for k, v in cfg.items() if k != "objects"}
    return {
        "schema_version": 1,
        "run_name": name,
        "status": None,
        "error": None,
        "inputs": {
            "diameter_mm": run.diameter_mm,
            "radius_m": run.radius_m,
            "initial_mass_kg": run.mass_kg,
            "cross_section_m2": run.cross_section_m2,
            "initial_temperature_K": run.temperature_K,
            "initial_velocity_kms": run.velocity_kms,
            "initial_altitude_km": run.altitude_km,
            "flight_path_angle_deg": run.flight_path_deg,
            "heading_deg": run.heading_deg % 360.0,
            "latitude_deg": run.lat_deg,
            "longitude_deg": run.lon_deg,
            "epoch_utc": run.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
            "material": MATERIAL_NAME,
            "material_density_kgm3": RHO_AA7075,
            "melting_temperature_K": T_MELT_AA7075,
            "melt_tolerance_K": MELT_TOLERANCE_K,
            "energy_threshold_J": ENERGY_THRESHOLD_J,
            "sesam_settings": json.loads(json.dumps(settings, default=str)),
        },
        "results": None,
        "warnings": [],
        "files": {"csv": os.path.abspath(csv_path), "raw_dir": run_raw},
        "provenance": {
            "drama_install_path": os.environ.get("DRAMA_INSTALL_PATH", DRAMA_INSTALL_PATH),
            "sesam_version": None,
            "pydrama_version": pydrama_version,
            "script_version": SCRIPT_VERSION,
            "wall_time_s": None,
            "created_utc": None,
            "hostname": socket.gethostname(),
        },
    }


def _finish(doc, json_path, status, t0, error=None):
    doc["status"] = status
    doc["error"] = error
    doc["provenance"]["wall_time_s"] = round(time.time() - t0, 3)
    doc["provenance"]["created_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(json_path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
    return doc


def run_sphere(run, outdir, raw_dir, timeout, keep_raw, fap_day_lines, fap_mon_lines):
    """Run SESAM for `run`; write <outdir>/<run_name>.csv and .json; return the JSON document."""
    sara, pydrama_version = _import_sara()
    name = run_name(run)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)
    csv_path = os.path.join(outdir, name + ".csv")
    json_path = os.path.join(outdir, name + ".json")
    run_raw = os.path.join(raw_dir, name)
    if os.path.isdir(run_raw):
        shutil.rmtree(run_raw)          # pyDRAMA's copytree refuses an existing destination
    if os.path.isfile(csv_path):
        os.remove(csv_path)
    cfg = build_config(run)
    doc = _base_document(run, cfg, name, csv_path, run_raw, pydrama_version)
    t0 = time.time()
    try:
        results = sara.run(config=[cfg], save_output_dirs=run_raw, keep_output_files="all",
                           fap_day_content=fap_day_lines, fap_mon_content=fap_mon_lines,
                           parallel=False, timeout=timeout, log_level="ERROR", spell_check=False)
    except Exception as exc:  # noqa: BLE001 -- anything pyDRAMA raises is a failed run
        return _finish(doc, json_path, "error", t0, "pyDRAMA raised {!r}".format(exc))

    errors = results.get("errors") or []
    if errors or not results.get("results"):
        status, messages = "error", []
        for err in errors:
            text = str(err.get("status"))
            if "timeout" in text.lower():
                status = "timeout"
            messages.append(text + "\n" + _log_tail(err.get("reentry_logfile") or err.get("logfile")))
        if not errors:
            messages.append("pyDRAMA returned no result")
        return _finish(doc, json_path, status, t0, "\n".join(messages))

    files = find_output_files(run_raw)
    if not files["aero"] or not files["traj"]:
        return _finish(doc, json_path, "error", t0, "history files not found under {}\n{}".format(
            run_raw, _log_tail(results["results"][0].get("reentry_logfile"))))

    aero = read_sara_table(files["aero"], AERO_COLUMNS)
    traj = read_sara_table(files["traj"], TRAJ_COLUMNS)
    rows, warnings = merge_histories(aero, traj)
    fragments = parse_impacting_fragments(files["fragments"])
    log_info = parse_sesam_log(files["log"])
    stats, more_warnings = compute_stats(run, rows, fragments, log_info)
    write_csv(csv_path, rows)
    doc["results"] = stats
    doc["warnings"] = warnings + more_warnings
    doc["provenance"]["sesam_version"] = parse_sesam_version(files["aero"])
    if keep_raw:
        doc["files"]["raw_dir"] = run_raw
    else:
        shutil.rmtree(run_raw, ignore_errors=True)
        doc["files"]["raw_dir"] = None
    return _finish(doc, json_path, "ok", t0)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: all pass (66 tests).

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Run SESAM through pyDRAMA and write the per-run JSON document

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 8: Command-line interface of sphere_reentry.py

**Files:**
- Modify: `sphere_reentry.py` (append)
- Modify: `tests/test_sphere_reentry.py` (append)

**Interfaces:**
- Produces: `parse_epoch(text) -> datetime` (accepts a trailing `Z`); `build_parser() -> argparse.ArgumentParser` with the options of spec §4.1; `main(argv=None) -> int` (0 ok, 1 failed run, 2 usage/environment error via `parser.error` or a printed message); `if __name__ == "__main__": sys.exit(main())`.

- [ ] **Step 1: Write the failing tests (append)**

```python
import subprocess
from helpers import REPO_ROOT, PY

BASE_ARGS = ["--velocity", "7.5", "--altitude", "77.5", "--temperature", "300", "--diameter", "50"]


def fap_files(tmp_path):
    day, mon = tmp_path / "fap_day.dat", tmp_path / "fap_mon.dat"
    day.write_text("# fap day\n01/08/2024 170 170 100 8 3 3 3 3 3 3 3 3\n")
    mon.write_text("# fap mon\n")
    return ["--fap-day", str(day), "--fap-mon", str(mon)]


class TestCli:
    def test_dry_run_prints_config_and_creates_nothing(self, tmp_path, capsys):
        rc = sr.main(BASE_ARGS + ["--dry-run", "--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")])
        assert rc == 0
        out = capsys.readouterr().out
        assert '"coordinateSystem": "geodetic"' in out
        assert "run name: sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km" in out
        assert not (tmp_path / "runs").exists() and not (tmp_path / "raw").exists()

    @pytest.mark.parametrize("argv", [
        ["--velocity", "7.5"],                                   # missing required
        BASE_ARGS[:6] + ["--diameter", "-5"],                    # negative diameter
        BASE_ARGS[:6] + ["--diameter", "50", "--epoch", "yesterday"],
        ["--velocity", "0", "--altitude", "77.5", "--temperature", "300", "--diameter", "50"],
    ])
    def test_bad_arguments_exit_2(self, argv):
        with pytest.raises(SystemExit) as exc:
            sr.main(argv)
        assert exc.value.code == 2

    def test_missing_fap_file_exits_2(self, tmp_path):
        with pytest.raises(SystemExit) as exc:
            sr.main(BASE_ARGS + ["--fap-day", str(tmp_path / "nope.dat"), "--outdir", str(tmp_path)])
        assert exc.value.code == 2

    def test_parse_epoch(self):
        assert sr.parse_epoch("2024-08-01T12:53:07Z") == datetime(2024, 8, 1, 12, 53, 7)
        assert sr.parse_epoch("2024-08-01T12:53:07") == datetime(2024, 8, 1, 12, 53, 7)

    def test_defaults_of_optional_state(self):
        args = sr.build_parser().parse_args(BASE_ARGS)
        assert (args.flight_path_angle, args.heading, args.lat, args.lon) == (0.0, 0.0, 0.0, 0.0)
        assert args.epoch == sr.PARENT_EPOCH and args.timeout == 600
        assert args.outdir == sr.DEFAULT_OUTDIR and args.raw_dir == sr.DEFAULT_RAW_DIR

    def test_ok_run_returns_0_and_prints_summary(self, tmp_path, monkeypatch, capsys):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm"))
        rc = sr.main(["--velocity", "0.5", "--altitude", "39.94", "--temperature", "300", "--diameter", "50",
                      "--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 0
        out = capsys.readouterr().out
        assert "survived" in out and "csv  ->" in out
        assert (tmp_path / "runs" / "sphere_d050.00mm_T0300.0K_v00.50000kms_h039.940km.json").is_file()

    def test_quiet_prints_nothing_on_success(self, tmp_path, monkeypatch, capsys):
        install_fake_drama(monkeypatch, fake_sara_run("T3_survivor_50mm"))
        rc = sr.main(["--velocity", "0.5", "--altitude", "39.94", "--temperature", "300", "--diameter", "50", "--quiet",
                      "--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 0 and capsys.readouterr().out == ""

    def test_failed_run_returns_1(self, tmp_path, monkeypatch, capsys):
        install_fake_drama(monkeypatch, fake_sara_run(mode="error"))
        rc = sr.main(BASE_ARGS + ["--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 1 and "error in reentry" in capsys.readouterr().err

    def test_drama_not_importable_returns_2(self, tmp_path, monkeypatch, capsys):
        monkeypatch.setitem(sys.modules, "drama", None)
        rc = sr.main(BASE_ARGS + ["--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path))
        assert rc == 2 and "pyDRAMA not importable" in capsys.readouterr().err


FAKE_DRAMA_SARA = '''
import os, shutil

def run(config, save_output_dirs, **kwargs):
    dest = os.path.join(save_output_dirs, "run_0", "reentry")
    if os.environ.get("FAKE_DRAMA_MODE") == "error":
        os.makedirs(dest, exist_ok=True)
        return {"config": config, "errors": [{"status": "error in reentry", "reentry_logfile": "boom"}], "results": []}
    shutil.copytree(os.environ["FAKE_DRAMA_FIXTURE"], dest)
    return {"config": config, "errors": [], "results": [{"status": "success", "reentry_logfile": ""}]}
'''


class TestCliSubprocess:
    def make_fake_package(self, tmp_path):
        pkg = tmp_path / "fakepkg" / "drama"
        pkg.mkdir(parents=True, exist_ok=True)
        (pkg / "__init__.py").write_text('__version__ = "fake"\n')
        (pkg / "sara.py").write_text(FAKE_DRAMA_SARA)
        return str(tmp_path / "fakepkg")

    def run_cli(self, tmp_path, extra, mode="ok"):
        env = dict(os.environ, PYTHONPATH=self.make_fake_package(tmp_path),
                   FAKE_DRAMA_FIXTURE=fixture_dir("T3_survivor_50mm"), FAKE_DRAMA_MODE=mode)
        cmd = [sys.executable, os.path.join(REPO_ROOT, "sphere_reentry.py")] + extra
        return subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(tmp_path))

    def test_exit_codes(self, tmp_path):
        common = ["--outdir", str(tmp_path / "runs"), "--raw-dir", str(tmp_path / "raw")] + fap_files(tmp_path)
        ok = self.run_cli(tmp_path, ["--velocity", "0.5", "--altitude", "39.94", "--temperature", "300", "--diameter", "50"] + common)
        assert ok.returncode == 0, ok.stderr
        assert (tmp_path / "runs" / "sphere_d050.00mm_T0300.0K_v00.50000kms_h039.940km.csv").is_file()
        failed = self.run_cli(tmp_path, BASE_ARGS + common, mode="error")
        assert failed.returncode == 1
        usage = self.run_cli(tmp_path, ["--velocity", "7.5"])
        assert usage.returncode == 2
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v -k "Cli"`
Expected: failures with `AttributeError: ... 'main'` / `'parse_epoch'`.

- [ ] **Step 3: Implement the CLI (append to sphere_reentry.py)**

```python
# =============================================================================
# 7. Command line (spec 4.1)
# =============================================================================

def parse_epoch(text):
    """ISO-8601 UTC timestamp (a trailing 'Z' is accepted) -> naive datetime."""
    text = text.strip()
    if text.endswith("Z"):
        text = text[:-1]
    try:
        return datetime.fromisoformat(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("invalid epoch {!r}: {}".format(text, exc))


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s]")
    p.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    p.add_argument("--temperature", type=float, required=True, help="initial bulk temperature [K]")
    p.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    p.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending (default 0)")
    p.add_argument("--heading", type=float, default=0.0, help="[deg] (default 0)")
    p.add_argument("--lat", type=float, default=0.0, help="geodetic latitude [deg] (default 0)")
    p.add_argument("--lon", type=float, default=0.0, help="longitude [deg] (default 0)")
    p.add_argument("--epoch", type=parse_epoch, default=PARENT_EPOCH,
                   help="ISO-8601 UTC (default %(default)s, the parent run's epoch)")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="CSV/JSON directory (default %(default)s)")
    p.add_argument("--raw-dir", default=DEFAULT_RAW_DIR, help="raw DRAMA trees root (default %(default)s)")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_S, help="SESAM timeout [s] (default %(default)s)")
    p.add_argument("--keep-raw", action="store_true", help="keep the raw DRAMA tree on success")
    p.add_argument("--dry-run", action="store_true", help="print the config and run name; run nothing")
    p.add_argument("--quiet", action="store_true", help="no console output except errors")
    p.add_argument("--fap-day", default=DEFAULT_FAP_DAY, help="DRAMA fap_day.dat (default %(default)s)")
    p.add_argument("--fap-mon", default=DEFAULT_FAP_MON, help="DRAMA fap_mon.dat (default %(default)s)")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.velocity <= 0.0:
        parser.error("--velocity must be > 0")
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.temperature <= 0.0:
        parser.error("--temperature must be > 0")
    if args.diameter <= 0.0:
        parser.error("--diameter must be > 0")
    run = SphereRun(velocity_kms=args.velocity, altitude_km=args.altitude,
                    temperature_K=args.temperature, diameter_mm=args.diameter,
                    flight_path_deg=args.flight_path_angle, heading_deg=args.heading,
                    lat_deg=args.lat, lon_deg=args.lon, epoch=args.epoch)
    if args.dry_run:
        print(json.dumps(build_config(run), indent=2, default=str))
        print("run name: " + run_name(run))
        return 0
    for label, path in (("--fap-day", args.fap_day), ("--fap-mon", args.fap_mon)):
        if not os.path.isfile(path):
            parser.error("{} file not found: {}".format(label, path))
    try:
        doc = run_sphere(run, args.outdir, args.raw_dir, args.timeout, args.keep_raw,
                         read_lines(args.fap_day), read_lines(args.fap_mon))
    except DramaNotAvailable as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    if doc["status"] != "ok":
        print("ERROR ({}): {}".format(doc["status"], doc["error"]), file=sys.stderr)
        return 1
    if not args.quiet:
        r = doc["results"]
        print("{}: {} after {:.1f} s ({}); Tmax {:.1f} K, {:.1f} s at melt, final mass {:.6g} kg "
              "(r = {:.3f} mm), final velocity {:.4f} km/s".format(
                  doc["run_name"], r["outcome"], r["final_time_s"], r["end_of_life_reason"],
                  r["max_temperature_K"], r["time_at_melting_temperature_s"], r["final_mass_kg"],
                  r["final_radius_mm"], r["final_velocity_kms"]))
        print("  csv  -> {}".format(doc["files"]["csv"]))
        print("  json -> {}".format(os.path.abspath(os.path.join(args.outdir, doc["run_name"] + ".json"))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_reentry.py -v`
Expected: all pass (79 tests). Also run `"$PY" sphere_reentry.py --velocity 7.5 --altitude 77.5 --temperature 300 --diameter 50 --dry-run | tail -3` and check the last line is `run name: sphere_d050.00mm_T0300.0K_v07.50000kms_h077.500km`.

- [ ] **Step 5: Commit**

```bash
git add sphere_reentry.py tests/test_sphere_reentry.py
git commit -m "Add the sphere_reentry command-line interface

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 9: Integration tests of script 1 against real DRAMA

**Files:**
- Create: `tests/test_integration_drama.py`

**Interfaces:**
- Consumes: `sphere_reentry.run_sphere`, `SphereRun`, `read_lines`, `CSV_COLUMNS`; `data/fap_*.dat`.
- Produces: the `drama`-marked tests that later tasks extend (Task 16 appends the sweep test to this file).

- [ ] **Step 1: Write the integration tests**

`tests/test_integration_drama.py`:
```python
"""Integration tests that run the real SESAM through pyDRAMA (marker: drama)."""
import csv
import json
import math
import os
import subprocess
import sys
from datetime import datetime

import pytest

import sphere_reentry as sr
from helpers import REPO_ROOT, PY

pytestmark = pytest.mark.drama

FAP_DAY = sr.read_lines(os.path.join(REPO_ROOT, "data", "fap_day.dat"))
FAP_MON = sr.read_lines(os.path.join(REPO_ROOT, "data", "fap_mon.dat"))

# parent state at 7.5 km/s (spike case T2/T5) and at 0.3 km/s (spike case E3)
ENTRY_STATE = dict(altitude_km=77.5, lat_deg=29.546, lon_deg=-82.134, flight_path_deg=-0.959,
                   heading_deg=347.168, epoch=datetime(2024, 8, 1, 12, 53, 7))
LOW_STATE = dict(altitude_km=33.229, lat_deg=35.99, lon_deg=-83.90, flight_path_deg=-60.2,
                 heading_deg=347.3, epoch=datetime(2024, 8, 1, 12, 56, 17))


def first_csv_row(path):
    with open(path) as fh:
        return next(csv.DictReader(fh))


def test_real_run_50mm_at_entry_speed_demises(tmp_path):
    run = sr.SphereRun(velocity_kms=7.5, temperature_K=300.0, diameter_mm=50.0, **ENTRY_STATE)
    doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, FAP_DAY, FAP_MON)
    assert doc["status"] == "ok", doc["error"]
    r = doc["results"]
    for key in ("max_temperature_K", "final_mass_kg", "time_at_melting_temperature_s",
                "final_velocity_kms", "final_radius_mm"):
        assert r[key] is not None and math.isfinite(r[key]), key
    assert r["max_temperature_K"] == 850.0
    assert r["outcome"] == "demised" and r["end_of_life_reason"] in ("ballooning", "uncritical")
    assert r["time_at_melting_temperature_s"] > 0.0
    assert r["final_radius_mm"] < 1.0
    assert not (tmp_path / "raw" / doc["run_name"]).exists()
    assert doc["provenance"]["sesam_version"] == "2.3.0"
    row = first_csv_row(doc["files"]["csv"])
    assert float(row["time_s"]) == 0.0
    assert float(row["altitude_km"]) == pytest.approx(77.5, abs=1e-3)
    assert float(row["velocity_kms"]) == pytest.approx(7.5, abs=1e-3)
    assert float(row["temperature_K"]) == 300.0
    assert float(row["lat_deg"]) == pytest.approx(29.546, abs=1e-3)
    assert float(row["lon_deg"]) == pytest.approx(-82.134, abs=1e-3)
    assert float(row["flight_path_deg"]) == pytest.approx(-0.959, abs=1e-3)
    assert float(row["heading_deg"]) == pytest.approx(347.168, abs=1e-3)


def test_real_run_5mm_low_altitude_survives_to_ground(tmp_path):
    run = sr.SphereRun(velocity_kms=0.3, temperature_K=300.0, diameter_mm=5.0, **LOW_STATE)
    doc = sr.run_sphere(run, str(tmp_path / "runs"), str(tmp_path / "raw"), 600, False, FAP_DAY, FAP_MON)
    assert doc["status"] == "ok", doc["error"]
    r = doc["results"]
    assert r["outcome"] == "survived" and r["end_of_life_reason"] == "ground impact"
    assert r["final_mass_source"] == "impacting_fragments_xml"
    assert r["final_mass_kg"] == pytest.approx(doc["inputs"]["initial_mass_kg"], abs=1e-9)
    assert r["final_radius_mm"] == pytest.approx(2.5, abs=1e-6)
    assert r["final_altitude_km"] == pytest.approx(0.0, abs=1e-3)
    assert r["max_temperature_K"] < 320.0
    assert r["n_rows"] > 100
    assert any("coarse" in w for w in doc["warnings"])


def test_cli_real_run_and_keep_raw(tmp_path):
    cmd = [PY, os.path.join(REPO_ROOT, "sphere_reentry.py"),
           "--velocity", "0.3", "--altitude", "33.229", "--temperature", "300", "--diameter", "5",
           "--flight-path-angle", "-60.2", "--heading", "347.3", "--lat", "35.99", "--lon", "-83.90",
           "--epoch", "2024-08-01T12:56:17", "--outdir", str(tmp_path / "runs"),
           "--raw-dir", str(tmp_path / "raw"), "--keep-raw"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    name = "sphere_d005.00mm_T0300.0K_v00.30000kms_h033.229km"
    assert (tmp_path / "runs" / (name + ".csv")).is_file()
    doc = json.load(open(tmp_path / "runs" / (name + ".json")))
    assert doc["status"] == "ok" and doc["files"]["raw_dir"] == str(tmp_path / "raw" / name)
    assert (tmp_path / "raw" / name / "run_0" / "reentry" / "sesam.log").is_file()
    assert "survived" in proc.stdout
```

- [ ] **Step 2: Run the integration tests**

Run: `"$PY" -m pytest tests/test_integration_drama.py -v`
Expected: 3 passed in under a minute (each SESAM run takes well under a second; the CLI test adds interpreter start-up). If SESAM refuses to run, print `doc["error"]` — it contains the last 40 log lines.

- [ ] **Step 3: Run the whole suite**

Run: `"$PY" -m pytest -v`
Expected: all tests pass (unit + integration).

- [ ] **Step 4: Commit**

```bash
git add tests/test_integration_drama.py
git commit -m "Add real-DRAMA integration tests for sphere_reentry

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 10: sphere_sweep — parent trajectory loading and state interpolation

**Files:**
- Create: `sphere_sweep.py`
- Create: `tests/fixtures/mini_dmf_output.json`
- Create: `tests/test_sphere_sweep.py`

**Interfaces:**
- Consumes: `sphere_reentry` (imported as `sr`; used from Task 11 on).
- Plateau rule (spec §5.2): the branch is non-increasing, not strictly decreasing — SESAM prints velocity with 3 decimals, so the real parent's last rows all read 0.028 km/s. `interpolate_state` takes the first bracketing pair, i.e. the highest-altitude row of a plateau (for the real parent, v = 0.028 km/s maps to 0.153 km, while `describe_v_min` reports the last row at 0.012 km). Physically immaterial (141 m at 28 m/s); documented so nobody "fixes" it.
- Produces (spec §5.2): constants `SCRIPT_DIR`, `SPHERE_REENTRY`, `DEFAULT_PARENT`, `DEFAULT_PARENT_OBJECT`, `DEFAULT_OUTDIR`, `DEFAULT_TIMEOUT_S`, `V_TOP_KMS = 7.5`, `N_VELOCITIES = 100`, `DIAMETERS_MM`, `TEMPERATURES_K`, `SECONDS_PER_RUN_ESTIMATE = 0.3`, `EPOCH_ARG_FORMAT = "%Y-%m-%dT%H:%M:%S"`; `sha256_of_file(path) -> str`; `descending_branch_start(velocities) -> int`; `@dataclass ParentTrajectory(path, sha256, object_name, epoch, rows, branch_start)` with properties `branch`, `v_max`, `v_min`, `v_min_row` and method `describe_v_min() -> str`; `load_parent(path, object_name) -> ParentTrajectory` (raises `ValueError` on ambiguity); `lerp(a, b, f)`, `angle_lerp(a_deg, b_deg, f) -> float in (-180, 180]`, `wrap_longitude(lon) -> float in [-180, 180)`; `interpolate_state(parent, v_kms) -> dict | None` with keys `altitude_km, lat_deg, lon_deg, flight_path_deg, heading_deg, time_s, epoch`.

- [ ] **Step 1: Create the miniature parent file**

`tests/fixtures/mini_dmf_output.json` (hand-built; same structure as a DRAMA GUI `dmf_output.json`; the Main_Body velocity rises once during an orbital phase (rows 2→3) and then decreases strictly to a non-zero impact velocity; heading crosses 360° between rows 5 and 6; longitude crosses the antimeridian between rows 6 and 7):
```json
{
  "controlOutputs": {"runStatus": "Finished", "progress": 100, "logs": []},
  "singleModuleOutputs": {
    "satellites": [
      {
        "name": "Mini Sat",
        "id": "sat-1",
        "missionPhases": [
          {
            "name": "Operational",
            "id": "phase-1",
            "epochs": [
              {
                "epoch": "2024-08-01T12:00:00.00Z",
                "analysisModules": [
                  {
                    "analysisModule": "sara",
                    "process_id": ["p-1"],
                    "results": {
                      "PySara.Main_Body.aaaa0000-0000-4000-8000-000000000001_Trajectory.txt": [
                        {"time": 0.0,    "altitude": 125.0, "lat": 10.0,  "lon": 117.0,  "velocity": 7.910, "downrange": 0.0,    "path": -0.06, "heading": 192.0},
                        {"time": 100.0,  "altitude": 124.0, "lat": 5.0,   "lon": 116.0,  "velocity": 7.905, "downrange": 790.0,  "path": 0.02,  "heading": 192.0},
                        {"time": 200.0,  "altitude": 126.0, "lat": 0.0,   "lon": 115.0,  "velocity": 7.900, "downrange": 1580.0, "path": 0.08,  "heading": 193.0},
                        {"time": 300.0,  "altitude": 128.0, "lat": -30.0, "lon": 100.0,  "velocity": 7.905, "downrange": 2370.0, "path": -0.10, "heading": 200.0},
                        {"time": 400.0,  "altitude": 120.0, "lat": -20.0, "lon": 110.0,  "velocity": 7.900, "downrange": 3160.0, "path": -0.20, "heading": 210.0},
                        {"time": 500.0,  "altitude": 100.0, "lat": -10.0, "lon": 120.0,  "velocity": 7.800, "downrange": 3950.0, "path": -0.50, "heading": 350.0},
                        {"time": 600.0,  "altitude": 80.0,  "lat": 0.0,   "lon": 178.0,  "velocity": 7.000, "downrange": 4700.0, "path": -1.00, "heading": 10.0},
                        {"time": 700.0,  "altitude": 60.0,  "lat": 10.0,  "lon": -178.0, "velocity": 5.000, "downrange": 5300.0, "path": -3.00, "heading": 20.0},
                        {"time": 800.0,  "altitude": 50.0,  "lat": 15.0,  "lon": -170.0, "velocity": 2.000, "downrange": 5650.0, "path": -10.0, "heading": 30.0},
                        {"time": 900.0,  "altitude": 40.0,  "lat": 20.0,  "lon": -160.0, "velocity": 1.000, "downrange": 5800.0, "path": -30.0, "heading": 40.0},
                        {"time": 1000.0, "altitude": 30.0,  "lat": 22.0,  "lon": -155.0, "velocity": 0.500, "downrange": 5870.0, "path": -50.0, "heading": 50.0},
                        {"time": 1100.0, "altitude": 20.0,  "lat": 23.0,  "lon": -152.0, "velocity": 0.200, "downrange": 5900.0, "path": -70.0, "heading": 60.0},
                        {"time": 1200.0, "altitude": 10.0,  "lat": 24.0,  "lon": -151.0, "velocity": 0.100, "downrange": 5910.0, "path": -85.0, "heading": 70.0},
                        {"time": 1300.0, "altitude": 0.0,   "lat": 25.0,  "lon": -150.0, "velocity": 0.050, "downrange": 5912.0, "path": -90.0, "heading": 80.0}
                      ],
                      "PySara.Main_Body.aaaa0000-0000-4000-8000-000000000001_AeroThermalHistory.txt": [
                        {"time": 0.0, "altitude": 125.0, "temp": 300.0, "mass": 740.0, "thick": 7.0},
                        {"time": 1300.0, "altitude": 0.0, "temp": 400.0, "mass": 81.0, "thick": 1.0}
                      ],
                      "PySara.Other_Part.bbbb0000-0000-4000-8000-000000000002_Trajectory.txt": [
                        {"time": 0.0,   "altitude": 90.0, "lat": 1.0, "lon": 2.0, "velocity": 7.500, "downrange": 0.0,   "path": -1.0, "heading": 100.0},
                        {"time": 100.0, "altitude": 70.0, "lat": 2.0, "lon": 3.0, "velocity": 6.000, "downrange": 600.0, "path": -2.0, "heading": 100.0},
                        {"time": 200.0, "altitude": 50.0, "lat": 3.0, "lon": 4.0, "velocity": 3.000, "downrange": 900.0, "path": -5.0, "heading": 100.0}
                      ],
                      "PySara.RiskResults.dat": {}
                    }
                  }
                ]
              }
            ]
          }
        ]
      }
    ]
  }
}
```

- [ ] **Step 2: Write the failing tests**

`tests/test_sphere_sweep.py`:
```python
"""Unit tests for sphere_sweep.py (no DRAMA needed)."""
import copy
import csv
import hashlib
import json
import os
from datetime import datetime

import pytest

import sphere_reentry as sr
import sphere_sweep as sw
from helpers import FIXTURES

MINI = os.path.join(FIXTURES, "mini_dmf_output.json")


@pytest.fixture
def parent():
    return sw.load_parent(MINI, "Main Body")


class TestLoadParent:
    def test_epoch_object_and_branch(self, parent):
        assert parent.epoch == datetime(2024, 8, 1, 12, 0, 0)
        assert parent.object_name == "Main Body"
        assert len(parent.rows) == 14
        assert parent.branch_start == 3
        assert parent.v_max == 7.905 and parent.v_min == 0.05
        assert parent.v_min_row == 13
        assert [r["velocity"] for r in parent.branch] == [7.905, 7.9, 7.8, 7.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]

    def test_sha256_matches_file(self, parent):
        with open(MINI, "rb") as fh:
            assert parent.sha256 == hashlib.sha256(fh.read()).hexdigest()
        assert parent.path == MINI

    def test_describe_v_min(self, parent):
        assert parent.describe_v_min() == "v_min = 0.050000 km/s at t = 1300.000 s, altitude 0.000 km (parent row 13)"

    def test_other_object_by_name_with_spaces(self):
        other = sw.load_parent(MINI, "Other Part")
        assert len(other.rows) == 3 and other.branch_start == 0 and other.v_max == 7.5

    def test_unknown_object_lists_available_keys(self):
        with pytest.raises(ValueError) as exc:
            sw.load_parent(MINI, "Nope")
        assert "Nope" in str(exc.value) and "Main_Body" in str(exc.value)

    def test_two_satellites_is_an_error(self, tmp_path):
        doc = json.load(open(MINI))
        doc["singleModuleOutputs"]["satellites"].append(copy.deepcopy(doc["singleModuleOutputs"]["satellites"][0]))
        p = tmp_path / "two.json"
        p.write_text(json.dumps(doc))
        with pytest.raises(ValueError, match="exactly one satellite"):
            sw.load_parent(str(p), "Main Body")

    def test_branch_with_one_row_is_an_error(self, tmp_path):
        doc = json.load(open(MINI))
        key = [k for k in doc["singleModuleOutputs"]["satellites"][0]["missionPhases"][0]["epochs"][0]["analysisModules"][0]["results"] if "Main_Body" in k and "Trajectory" in k][0]
        rows = doc["singleModuleOutputs"]["satellites"][0]["missionPhases"][0]["epochs"][0]["analysisModules"][0]["results"][key]
        rows[-1]["velocity"] = 9.0          # velocity increases at the very end -> branch of 1 row
        p = tmp_path / "bad.json"
        p.write_text(json.dumps(doc))
        with pytest.raises(ValueError, match="fewer than 2 rows"):
            sw.load_parent(str(p), "Main Body")


class TestDescendingBranch:
    def test_start_index(self):
        assert sw.descending_branch_start([7.9, 7.8, 7.7]) == 0
        assert sw.descending_branch_start([7.9, 7.8, 7.85, 7.7, 7.6]) == 2
        assert sw.descending_branch_start([7.9, 7.95, 7.8, 7.85, 7.7]) == 3
        assert sw.descending_branch_start([1.0]) == 0


class TestInterpolation:
    def test_angle_lerp_and_wrap(self):
        assert sw.angle_lerp(350.0, 10.0, 0.5) == pytest.approx(0.0, abs=1e-9)
        assert sw.angle_lerp(90.0, 270.0, 0.0) == pytest.approx(90.0)
        assert sw.angle_lerp(10.0, 20.0, 0.5) == pytest.approx(15.0)
        assert sw.wrap_longitude(190.0) == -170.0
        assert sw.wrap_longitude(-190.0) == 170.0
        assert sw.wrap_longitude(180.0) == -180.0
        assert sw.wrap_longitude(-83.9) == pytest.approx(-83.9)
        assert sw.lerp(10.0, 20.0, 0.25) == 12.5

    def test_midpoint_between_rows(self, parent):
        st = sw.interpolate_state(parent, 7.4)          # between 7.8 (t=500) and 7.0 (t=600)
        assert st["altitude_km"] == pytest.approx(90.0)
        assert st["lat_deg"] == pytest.approx(-5.0)
        assert st["flight_path_deg"] == pytest.approx(-0.75)
        assert st["time_s"] == pytest.approx(550.0)
        assert min(st["heading_deg"], 360.0 - st["heading_deg"]) == pytest.approx(0.0, abs=1e-9)  # 350 -> 10 wraps through 0
        assert st["lon_deg"] == pytest.approx(149.0)
        assert st["epoch"] == datetime(2024, 8, 1, 12, 9, 10)

    def test_longitude_wraps_across_antimeridian(self, parent):
        st = sw.interpolate_state(parent, 6.0)          # lon 178 -> -178
        assert abs(st["lon_deg"]) == pytest.approx(180.0)
        assert -180.0 <= st["lon_deg"] < 180.0

    def test_exact_row_velocities(self, parent):
        assert sw.interpolate_state(parent, 5.0)["altitude_km"] == pytest.approx(60.0)
        assert sw.interpolate_state(parent, 7.905)["altitude_km"] == pytest.approx(128.0)
        assert sw.interpolate_state(parent, 0.05)["altitude_km"] == pytest.approx(0.0)
        assert sw.interpolate_state(parent, 0.05)["time_s"] == pytest.approx(1300.0)

    def test_outside_range_is_none(self, parent):
        assert sw.interpolate_state(parent, 7.906) is None
        assert sw.interpolate_state(parent, 0.049) is None
        assert sw.interpolate_state(parent, 7.905) is not None
        assert sw.interpolate_state(parent, 0.05) is not None

    def test_orbital_phase_is_never_sampled(self, parent):
        # 7.902 km/s occurs in the orbital phase (rows 1-2) AND on the branch (rows 3-4);
        # only the branch is used: altitude between 128 and 120, not between 124 and 126.
        st = sw.interpolate_state(parent, 7.902)
        assert 120.0 < st["altitude_km"] < 128.0 and st["time_s"] > 300.0
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'sphere_sweep'`.

- [ ] **Step 4: Write sphere_sweep.py (header, constants, parent loading, interpolation)**

```python
#!/usr/bin/env python3
"""
sphere_sweep.py -- run sphere_reentry.py over a grid of diameters, temperatures and velocities

The sphere's remaining initial-state values (altitude, latitude, longitude,
flight-path angle, heading, epoch) are inherited from a parent-satellite SESAM
run (a DRAMA GUI dmf_output.json) at the point of its final descent where the
parent reached the sphere's velocity. The velocity grid runs from 7.5 km/s down
to the parent's minimum velocity (its impact velocity), which is read from the
parent file and logged.

Usage
-----
    python sphere_sweep.py --dry-run                       # show the matrix, run nothing
    python sphere_sweep.py                                 # everything pending in one confirmed batch
    python sphere_sweep.py --batch-size 5000               # confirm + choose cores before every batch
    python sphere_sweep.py --diameters 50 --temperatures 300 --velocities 7.5,0.5 --yes --cores 2

Outputs (under --outdir, default sphere_sweep_output/):
    runs/<run_name>.csv|.json   per run (written by sphere_reentry.py)
    raw/<run_name>/             raw DRAMA tree, failed runs only
    sweep_manifest.json         every matrix point with its state and status
    sweep_summary.csv           one row per point with the statistics
    sweep.log                   log of prompts, answers, failures

Requires conda env drama_env (pyDRAMA, numpy, tqdm); see README.md.
Design: docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta, timezone

import numpy as np
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
import sphere_reentry as sr  # noqa: E402

SPHERE_REENTRY = os.path.join(SCRIPT_DIR, "sphere_reentry.py")
DEFAULT_PARENT = os.path.join(SCRIPT_DIR, "Generic_Satellite Reentry", "output", "dmf_output.json")
DEFAULT_PARENT_OBJECT = "Main Body"
DEFAULT_OUTDIR = os.path.join(SCRIPT_DIR, "sphere_sweep_output")
DEFAULT_TIMEOUT_S = 600
V_TOP_KMS = 7.5                    # fixed upper bound of the velocity grid
N_VELOCITIES = 100                 # linspace(V_TOP_KMS, v_min_parent, N_VELOCITIES)
DIAMETERS_MM = [float(d) for d in range(5, 101, 5)]          # 5 .. 100 mm step 5
TEMPERATURES_K = [float(t) for t in range(300, 751, 10)]     # 300 .. 750 K step 10
SECONDS_PER_RUN_ESTIMATE = 0.3     # per run per core, before any run has been timed
STDERR_TAIL_LINES = 20
EPOCH_ARG_FORMAT = "%Y-%m-%dT%H:%M:%S"
PARENT_EPOCH_FORMAT = "%Y-%m-%dT%H:%M:%S.%fZ"

log = logging.getLogger("sweep")


# =============================================================================
# 1. Parent trajectory (spec 5.2)
# =============================================================================

def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def descending_branch_start(velocities):
    """Index one past the last velocity increase (0 when velocity never increases)."""
    start = 0
    for i in range(len(velocities) - 1):
        if velocities[i + 1] > velocities[i]:
            start = i + 1
    return start


@dataclass
class ParentTrajectory:
    path: str
    sha256: str
    object_name: str
    epoch: datetime
    rows: list            # full trajectory rows (dicts straight from dmf_output.json)
    branch_start: int     # first row of the final descending branch

    @property
    def branch(self):
        return self.rows[self.branch_start:]

    @property
    def v_max(self):
        return self.branch[0]["velocity"]

    @property
    def v_min(self):
        return self.branch[-1]["velocity"]

    @property
    def v_min_row(self):
        return len(self.rows) - 1

    def describe_v_min(self):
        last = self.branch[-1]
        return "v_min = {:.6f} km/s at t = {:.3f} s, altitude {:.3f} km (parent row {})".format(
            last["velocity"], last["time"], last["altitude"], self.v_min_row)


def load_parent(path, object_name):
    """Read a DRAMA GUI dmf_output.json and select one object's trajectory."""
    with open(path) as fh:
        doc = json.load(fh)
    satellites = doc["singleModuleOutputs"]["satellites"]
    if len(satellites) != 1:
        raise ValueError("expected exactly one satellite in {}, found {}: {}".format(
            path, len(satellites), [s.get("name") for s in satellites]))
    phases = satellites[0]["missionPhases"]
    if len(phases) != 1:
        raise ValueError("expected exactly one mission phase, found {}: {}".format(
            len(phases), [p.get("name") for p in phases]))
    epochs = phases[0]["epochs"]
    if len(epochs) != 1:
        raise ValueError("expected exactly one epoch, found {}: {}".format(
            len(epochs), [e.get("epoch") for e in epochs]))
    epoch = datetime.strptime(epochs[0]["epoch"], PARENT_EPOCH_FORMAT)
    results = epochs[0]["analysisModules"][0]["results"]
    wanted = object_name.replace(" ", "_")
    trajectory_keys = [k for k in results if k.endswith("_Trajectory.txt")]
    keys = [k for k in trajectory_keys if len(k.split(".")) >= 2 and k.split(".")[1] == wanted]
    if len(keys) != 1:
        raise ValueError("expected one trajectory for object {!r}, found {}; available: {}".format(
            object_name, keys, trajectory_keys))
    rows = results[keys[0]]
    if len(rows) < 2:
        raise ValueError("parent trajectory {} has fewer than 2 rows".format(keys[0]))
    start = descending_branch_start([r["velocity"] for r in rows])
    if len(rows) - start < 2:
        raise ValueError("descending branch of the parent trajectory has fewer than 2 rows")
    return ParentTrajectory(path=path, sha256=sha256_of_file(path), object_name=object_name,
                            epoch=epoch, rows=rows, branch_start=start)


def lerp(a, b, f):
    return a + f * (b - a)


def angle_lerp(a_deg, b_deg, f):
    """Interpolate two angles through their sines/cosines; result in (-180, 180]."""
    ra, rb = math.radians(a_deg), math.radians(b_deg)
    s = lerp(math.sin(ra), math.sin(rb), f)
    c = lerp(math.cos(ra), math.cos(rb), f)
    return math.degrees(math.atan2(s, c))


def wrap_longitude(lon):
    return ((lon + 180.0) % 360.0) - 180.0


def interpolate_state(parent, v_kms):
    """Parent state at velocity v on the descending branch (linear in velocity), or None outside it."""
    if v_kms > parent.v_max or v_kms < parent.v_min:
        return None
    branch = parent.branch
    for a, b in zip(branch, branch[1:]):
        if a["velocity"] >= v_kms >= b["velocity"]:
            dv = a["velocity"] - b["velocity"]
            f = (a["velocity"] - v_kms) / dv if dv else 0.0
            t = lerp(a["time"], b["time"], f)
            return {
                "altitude_km": lerp(a["altitude"], b["altitude"], f),
                "lat_deg": lerp(a["lat"], b["lat"], f),
                "lon_deg": wrap_longitude(angle_lerp(a["lon"], b["lon"], f)),
                "flight_path_deg": lerp(a["path"], b["path"], f),
                "heading_deg": angle_lerp(a["heading"], b["heading"], f) % 360.0,
                "time_s": t,
                "epoch": parent.epoch + timedelta(seconds=t),
            }
    return None
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: 14 passed.

- [ ] **Step 6: Commit**

```bash
git add sphere_sweep.py tests/test_sphere_sweep.py tests/fixtures/mini_dmf_output.json
git commit -m "Add sphere_sweep parent-trajectory loading and state interpolation

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 11: Grid, matrix points and run names

**Files:**
- Modify: `sphere_sweep.py` (append)
- Modify: `tests/test_sphere_sweep.py` (append)

**Interfaces:**
- Consumes: `interpolate_state`, `ParentTrajectory` (Task 10); `sr.SphereRun`, `sr.run_name` (Task 2).
- Produces (spec §5.3): `velocity_grid(v_min_parent) -> list[float]`; `parse_float_list(text) -> list[float]` (argparse type); `format_arg(x) -> str` (`"{:.6f}"`); `@dataclass Point(diameter_mm, temperature_K, velocity_kms, altitude_km=None, flight_path_deg=None, heading_deg=None, lat_deg=None, lon_deg=None, epoch_utc=None, run_name=None, status="pending", skip_reason=None, returncode=None, wall_time_s=None, stderr_tail=None)`; `sphere_run_for(point) -> sr.SphereRun` (built from the formatted CLI strings parsed back); `build_points(parent, diameters, temperatures, velocities, limit=None) -> list[Point]`.

- [ ] **Step 1: Write the failing tests (append)**

```python
class TestGrid:
    def test_velocity_grid(self):
        grid = sw.velocity_grid(0.05)
        assert len(grid) == 100 and grid[0] == 7.5 and grid[-1] == pytest.approx(0.05)
        assert all(a > b for a, b in zip(grid, grid[1:]))
        assert grid[0] - grid[1] == pytest.approx((7.5 - 0.05) / 99)
        assert all(isinstance(v, float) for v in grid)

    def test_default_diameters_and_temperatures(self):
        assert sw.DIAMETERS_MM == [float(d) for d in range(5, 101, 5)] and len(sw.DIAMETERS_MM) == 20
        assert sw.TEMPERATURES_K[0] == 300.0 and sw.TEMPERATURES_K[-1] == 750.0 and len(sw.TEMPERATURES_K) == 46
        assert sw.V_TOP_KMS == 7.5 and sw.N_VELOCITIES == 100

    def test_parse_float_list(self):
        assert sw.parse_float_list("5, 10,15") == [5.0, 10.0, 15.0]
        with pytest.raises(Exception):
            sw.parse_float_list("")

    def test_format_arg(self):
        assert sw.format_arg(7.5) == "7.500000"
        assert sw.format_arg(0.0283333333) == "0.028333"
        assert sw.format_arg(-82.134) == "-82.134000"


class TestBuildPoints:
    def test_states_skips_and_ordering(self, parent):
        points = sw.build_points(parent, [10.0, 5.0], [310.0, 300.0], [7.4, 9.0, 0.01])
        assert [(p.diameter_mm, p.temperature_K, p.velocity_kms) for p in points] == [
            (5.0, 300.0, 9.0), (5.0, 300.0, 7.4), (5.0, 300.0, 0.01),
            (5.0, 310.0, 9.0), (5.0, 310.0, 7.4), (5.0, 310.0, 0.01),
            (10.0, 300.0, 9.0), (10.0, 300.0, 7.4), (10.0, 300.0, 0.01),
            (10.0, 310.0, 9.0), (10.0, 310.0, 7.4), (10.0, 310.0, 0.01)]
        skipped, ok = points[0], points[1]
        assert skipped.status == "skipped" and skipped.skip_reason == "parent never reaches velocity"
        assert skipped.run_name is None and skipped.altitude_km is None and skipped.epoch_utc is None
        assert ok.status == "pending" and ok.skip_reason is None
        assert ok.altitude_km == pytest.approx(90.0) and ok.flight_path_deg == pytest.approx(-0.75)
        assert ok.lat_deg == pytest.approx(-5.0) and ok.lon_deg == pytest.approx(149.0)
        assert min(ok.heading_deg, 360.0 - ok.heading_deg) == pytest.approx(0.0, abs=1e-9)
        assert ok.epoch_utc == "2024-08-01T12:09:10"
        assert ok.run_name == "sphere_d005.00mm_T0300.0K_v07.40000kms_h090.000km"
        assert points[2].status == "skipped"

    def test_run_name_matches_script1_on_formatted_args(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        run = sw.sphere_run_for(point)
        assert run.velocity_kms == 7.4 and run.altitude_km == 90.0 and run.epoch == datetime(2024, 8, 1, 12, 9, 10)
        assert run.flight_path_deg == float("-0.750000") and run.heading_deg % 360.0 == 0.0
        assert sr.run_name(run) == point.run_name

    def test_limit(self, parent):
        assert len(sw.build_points(parent, [5.0, 10.0], [300.0], [7.4, 5.0], limit=3)) == 3

    def test_full_default_matrix_size(self, parent):
        points = sw.build_points(parent, sw.DIAMETERS_MM, sw.TEMPERATURES_K, sw.velocity_grid(parent.v_min))
        assert len(points) == 20 * 46 * 100
        assert sum(p.status == "skipped" for p in points) == 0
        assert len({p.run_name for p in points}) == len(points)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v -k "Grid or BuildPoints"`
Expected: failures with `AttributeError: ... 'velocity_grid'` / `'build_points'`.

- [ ] **Step 3: Implement grid and points (append to sphere_sweep.py)**

```python
# =============================================================================
# 2. Grid and matrix (spec 5.3)
# =============================================================================

def velocity_grid(v_min_parent):
    """100 velocities from 7.5 km/s down to the parent's minimum velocity, both included."""
    return [float(v) for v in np.linspace(V_TOP_KMS, v_min_parent, N_VELOCITIES)]


def parse_float_list(text):
    values = [float(x) for x in text.split(",") if x.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected a comma-separated list of numbers")
    return values


def format_arg(x):
    """How every numeric argument is passed to sphere_reentry.py."""
    return "{:.6f}".format(x)


@dataclass
class Point:
    diameter_mm: float
    temperature_K: float
    velocity_kms: float
    altitude_km: float = None
    flight_path_deg: float = None
    heading_deg: float = None
    lat_deg: float = None
    lon_deg: float = None
    epoch_utc: str = None
    run_name: str = None
    status: str = "pending"        # pending | skipped | done | failed
    skip_reason: str = None
    returncode: int = None
    wall_time_s: float = None
    stderr_tail: str = None


def sphere_run_for(point):
    """The SphereRun script 1 will construct from the formatted CLI strings (so run names agree)."""
    return sr.SphereRun(
        velocity_kms=float(format_arg(point.velocity_kms)),
        altitude_km=float(format_arg(point.altitude_km)),
        temperature_K=float(format_arg(point.temperature_K)),
        diameter_mm=float(format_arg(point.diameter_mm)),
        flight_path_deg=float(format_arg(point.flight_path_deg)),
        heading_deg=float(format_arg(point.heading_deg)),
        lat_deg=float(format_arg(point.lat_deg)),
        lon_deg=float(format_arg(point.lon_deg)),
        epoch=datetime.strptime(point.epoch_utc, EPOCH_ARG_FORMAT),
    )


def build_points(parent, diameters, temperatures, velocities, limit=None):
    """All matrix points ordered diameter -> temperature -> velocity (descending), with parent states."""
    points = []
    for d in sorted(diameters):
        for T in sorted(temperatures):
            for v in sorted(velocities, reverse=True):
                point = Point(diameter_mm=d, temperature_K=T, velocity_kms=v)
                state = interpolate_state(parent, v)
                if state is None:
                    point.status = "skipped"
                    point.skip_reason = "parent never reaches velocity"
                else:
                    point.altitude_km = state["altitude_km"]
                    point.flight_path_deg = state["flight_path_deg"]
                    point.heading_deg = state["heading_deg"]
                    point.lat_deg = state["lat_deg"]
                    point.lon_deg = state["lon_deg"]
                    point.epoch_utc = state["epoch"].strftime(EPOCH_ARG_FORMAT)
                    point.run_name = sr.run_name(sphere_run_for(point))
                points.append(point)
    if limit is not None:
        points = points[:limit]
    return points
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: 22 passed (the full-matrix test builds 92,000 points; it should take a few seconds at most).

- [ ] **Step 5: Commit**

```bash
git add sphere_sweep.py tests/test_sphere_sweep.py
git commit -m "Build the sweep matrix with parent-inherited states and script-1 run names

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 12: Manifest and resume

**Files:**
- Modify: `sphere_sweep.py` (append)
- Modify: `tests/test_sphere_sweep.py` (append)

**Interfaces:**
- Consumes: `Point`, `ParentTrajectory` (Tasks 10–11).
- Produces (spec §5.4): `utc_now_str() -> str`; `manifest_document(parent, grid, settings, points, created_utc=None) -> dict`; `write_manifest(path, doc)` (atomic via `.tmp` + `os.replace`); `load_manifest(path) -> dict`; `apply_resume(points, runs_dir, force=False, retry_failed=False) -> None` (sets `done`/`failed`/`pending` on non-skipped points from the JSON files in `runs_dir`).

Resume semantics (refines spec §5.4 wording): a point whose `runs/<run_name>.json` says `status: ok` is `done`; one whose JSON exists but is not ok (or is unreadable) is `failed`; no JSON → `pending`. `--force` makes every non-skipped point `pending`; `--retry-failed` makes `failed` points `pending` (never-run points are always `pending` and run too). Without `--retry-failed`, `failed` points are left alone so a systematic failure is not retried forever.

- [ ] **Step 1: Write the failing tests (append)**

```python
def write_run_json(runs_dir, run_name, status="ok", results=None):
    os.makedirs(runs_dir, exist_ok=True)
    doc = {"schema_version": 1, "run_name": run_name, "status": status, "results": results}
    with open(os.path.join(runs_dir, run_name + ".json"), "w") as fh:
        json.dump(doc, fh)


class TestManifest:
    def test_manifest_document(self, parent):
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 9.0])
        grid = {"diameters_mm": [5.0], "temperatures_K": [300.0], "velocities_kms": [7.4, 9.0]}
        settings = {"timeout_s": 600, "batch_size": None, "sphere_reentry_version": sr.SCRIPT_VERSION}
        doc = sw.manifest_document(parent, grid, settings, points, created_utc="2026-09-14T00:00:00Z")
        assert doc["created_utc"] == "2026-09-14T00:00:00Z" and doc["updated_utc"].endswith("Z")
        assert doc["parent"] == {
            "path": os.path.abspath(MINI), "sha256": parent.sha256, "object": "Main Body",
            "epoch_utc": "2024-08-01T12:00:00", "branch_first_row": 3, "branch_last_row": 13,
            "v_min_kms": 0.05, "v_min_time_s": 1300.0, "v_min_altitude_km": 0.0, "v_max_kms": 7.905}
        assert doc["grid"] == grid and doc["settings"] == settings
        assert [p["status"] for p in doc["points"]] == ["skipped", "pending"]
        assert doc["points"][1]["run_name"] == "sphere_d005.00mm_T0300.0K_v07.40000kms_h090.000km"
        assert set(doc["points"][1]) == set(sw.Point.__dataclass_fields__)

    def test_write_and_load_round_trip(self, parent, tmp_path):
        points = sw.build_points(parent, [5.0], [300.0], [7.4])
        doc = sw.manifest_document(parent, {}, {}, points)
        path = str(tmp_path / "sweep_manifest.json")
        sw.write_manifest(path, doc)
        assert sw.load_manifest(path) == json.loads(json.dumps(doc))
        assert not os.path.exists(path + ".tmp")


class TestResume:
    def test_statuses_from_existing_json(self, parent, tmp_path):
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0, 9.0])   # ordered 9.0, 7.4, 5.0, 2.0
        skipped, ok, failed, fresh = points
        write_run_json(runs, ok.run_name, "ok")
        write_run_json(runs, failed.run_name, "error")
        sw.apply_resume(points, runs)
        assert (ok.status, failed.status, fresh.status, skipped.status) == ("done", "failed", "pending", "skipped")

    def test_force_and_retry_failed(self, parent, tmp_path):
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0])
        write_run_json(runs, points[0].run_name, "ok")
        write_run_json(runs, points[1].run_name, "timeout")
        sw.apply_resume(points, runs, retry_failed=True)
        assert [p.status for p in points] == ["done", "pending", "pending"]
        sw.apply_resume(points, runs, force=True)
        assert [p.status for p in points] == ["pending", "pending", "pending"]

    def test_corrupt_json_counts_as_failed(self, parent, tmp_path):
        runs = tmp_path / "runs"
        runs.mkdir()
        points = sw.build_points(parent, [5.0], [300.0], [7.4])
        (runs / (points[0].run_name + ".json")).write_text("{not json")
        sw.apply_resume(points, str(runs))
        assert points[0].status == "failed"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v -k "Manifest or Resume"`
Expected: failures with `AttributeError: ... 'manifest_document'` / `'apply_resume'`.

- [ ] **Step 3: Implement manifest and resume (append to sphere_sweep.py)**

```python
# =============================================================================
# 3. Manifest and resume (spec 5.4)
# =============================================================================

def utc_now_str():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def manifest_document(parent, grid, settings, points, created_utc=None):
    now = utc_now_str()
    last = parent.branch[-1]
    return {
        "created_utc": created_utc or now,
        "updated_utc": now,
        "parent": {
            "path": os.path.abspath(parent.path),
            "sha256": parent.sha256,
            "object": parent.object_name,
            "epoch_utc": parent.epoch.strftime(EPOCH_ARG_FORMAT),
            "branch_first_row": parent.branch_start,
            "branch_last_row": len(parent.rows) - 1,
            "v_min_kms": parent.v_min,
            "v_min_time_s": last["time"],
            "v_min_altitude_km": last["altitude"],
            "v_max_kms": parent.v_max,
        },
        "grid": grid,
        "settings": settings,
        "points": [asdict(p) for p in points],
    }


def write_manifest(path, doc):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(doc, fh, indent=2)
    os.replace(tmp, path)


def load_manifest(path):
    with open(path) as fh:
        return json.load(fh)


def _json_status(path):
    """status field of a run JSON, or None when the file is missing/unreadable."""
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh).get("status") or "unreadable"
    except (OSError, ValueError):
        return "unreadable"


def apply_resume(points, runs_dir, force=False, retry_failed=False):
    """Mark points done/failed/pending from the run JSONs already in runs_dir."""
    for point in points:
        if point.status == "skipped":
            continue
        status = _json_status(os.path.join(runs_dir, point.run_name + ".json"))
        if status is None:
            point.status = "pending"
        elif status == "ok":
            point.status = "pending" if force else "done"
        else:
            point.status = "pending" if (force or retry_failed) else "failed"
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: 27 passed.

- [ ] **Step 5: Commit**

```bash
git add sphere_sweep.py tests/test_sphere_sweep.py
git commit -m "Add sweep manifest and resume logic

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 13: Aggregate summary CSV

**Files:**
- Modify: `sphere_sweep.py` (append)
- Modify: `tests/test_sphere_sweep.py` (append)

**Interfaces:**
- Consumes: `Point` (Task 11).
- Produces (spec §5.6): `SUMMARY_COLUMNS` (exact list); `summary_row(point, runs_dir) -> dict`; `write_summary(path, points, runs_dir)` (atomic).

- [ ] **Step 1: Write the failing tests (append)**

```python
class TestSummary:
    def test_columns_are_the_spec_list(self):
        assert sw.SUMMARY_COLUMNS == [
            "diameter_mm", "initial_temperature_K", "initial_velocity_kms", "initial_altitude_km",
            "flight_path_angle_deg", "heading_deg", "latitude_deg", "longitude_deg", "epoch_utc",
            "status", "skip_reason", "max_temperature_K", "final_mass_kg", "final_mass_source",
            "mass_loss_fraction", "time_at_melting_temperature_s", "final_velocity_kms",
            "final_radius_mm", "final_altitude_km", "end_of_life_reason", "outcome", "wall_time_s", "run_name"]

    def test_rows_for_each_status(self, parent, tmp_path):
        runs = str(tmp_path / "runs")
        points = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0, 9.0])   # ordered 9.0, 7.4, 5.0, 2.0
        skipped, done, failed, pending = points
        write_run_json(runs, done.run_name, "ok", results={
            "max_temperature_K": 850.0, "final_mass_kg": 3e-06, "final_mass_source": "sesam_log_event_end",
            "mass_loss_fraction": 0.98, "time_at_melting_temperature_s": 6.774, "final_velocity_kms": 2.998,
            "final_radius_mm": 0.633, "final_altitude_km": 76.548, "end_of_life_reason": "uncritical",
            "outcome": "demised"})
        done.status, done.wall_time_s = "done", 0.31
        failed.status, failed.returncode = "failed", 1
        path = str(tmp_path / "sweep_summary.csv")
        sw.write_summary(path, points, runs)
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            assert reader.fieldnames == sw.SUMMARY_COLUMNS
            rows = list(reader)
        assert len(rows) == 4
        s = rows[0]
        assert s["status"] == "skipped" and s["skip_reason"] == "parent never reaches velocity"
        assert s["run_name"] == "" and s["initial_altitude_km"] == "" and s["initial_velocity_kms"] == "9.0"
        d = rows[1]
        assert d["status"] == "done" and d["run_name"] == done.run_name and d["wall_time_s"] == "0.31"
        assert d["diameter_mm"] == "5.0" and d["initial_velocity_kms"] == "7.4"
        assert float(d["initial_altitude_km"]) == pytest.approx(90.0)
        assert d["epoch_utc"] == "2024-08-01T12:09:10"
        assert d["max_temperature_K"] == "850.0" and d["final_mass_source"] == "sesam_log_event_end"
        assert d["outcome"] == "demised" and d["final_radius_mm"] == "0.633"
        f = rows[2]
        assert f["status"] == "failed" and f["max_temperature_K"] == "" and f["outcome"] == "" and f["run_name"] == failed.run_name
        assert rows[3]["status"] == "pending" and rows[3]["final_mass_kg"] == ""
        assert not os.path.exists(path + ".tmp")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v -k Summary`
Expected: 2 failed (`AttributeError: ... 'SUMMARY_COLUMNS'`).

- [ ] **Step 3: Implement the summary (append to sphere_sweep.py)**

```python
# =============================================================================
# 4. Aggregate summary (spec 5.6)
# =============================================================================

SUMMARY_COLUMNS = [
    "diameter_mm", "initial_temperature_K", "initial_velocity_kms", "initial_altitude_km",
    "flight_path_angle_deg", "heading_deg", "latitude_deg", "longitude_deg", "epoch_utc",
    "status", "skip_reason", "max_temperature_K", "final_mass_kg", "final_mass_source",
    "mass_loss_fraction", "time_at_melting_temperature_s", "final_velocity_kms",
    "final_radius_mm", "final_altitude_km", "end_of_life_reason", "outcome", "wall_time_s", "run_name",
]
_RESULT_COLUMNS = ["max_temperature_K", "final_mass_kg", "final_mass_source", "mass_loss_fraction",
                   "time_at_melting_temperature_s", "final_velocity_kms", "final_radius_mm",
                   "final_altitude_km", "end_of_life_reason", "outcome"]


def _blank(value):
    return "" if value is None else value


def summary_row(point, runs_dir):
    row = {c: "" for c in SUMMARY_COLUMNS}
    row.update({
        "diameter_mm": point.diameter_mm,
        "initial_temperature_K": point.temperature_K,
        "initial_velocity_kms": point.velocity_kms,
        "initial_altitude_km": _blank(point.altitude_km),
        "flight_path_angle_deg": _blank(point.flight_path_deg),
        "heading_deg": _blank(point.heading_deg),
        "latitude_deg": _blank(point.lat_deg),
        "longitude_deg": _blank(point.lon_deg),
        "epoch_utc": _blank(point.epoch_utc),
        "status": point.status,
        "skip_reason": _blank(point.skip_reason),
        "wall_time_s": _blank(point.wall_time_s),
        "run_name": _blank(point.run_name),
    })
    if point.status == "done":
        try:
            with open(os.path.join(runs_dir, point.run_name + ".json")) as fh:
                results = json.load(fh).get("results") or {}
        except (OSError, ValueError):
            results = {}
        for column in _RESULT_COLUMNS:
            row[column] = _blank(results.get(column))
    return row


def write_summary(path, points, runs_dir):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for point in points:
            writer.writerow(summary_row(point, runs_dir))
    os.replace(tmp, path)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: 29 passed.

- [ ] **Step 5: Commit**

```bash
git add sphere_sweep.py tests/test_sphere_sweep.py
git commit -m "Write the sweep summary CSV from run JSON files

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 14: Batched execution — subprocess runner, prompts, thread pool, progress

**Files:**
- Modify: `sphere_sweep.py` (append)
- Modify: `tests/test_sphere_sweep.py` (append)

**Interfaces:**
- Consumes: `Point`, `format_arg`, `SPHERE_REENTRY`, `STDERR_TAIL_LINES`, `log` (Tasks 10–11).
- Produces (spec §5.5): `build_command(point, runs_dir, raw_dir, timeout, python=None) -> list[str]`; `run_point(point, runs_dir, raw_dir, timeout) -> Point` (runs script 1 in a subprocess, sets `returncode`, `wall_time_s`, `stderr_tail`, `status`); `split_batches(points, batch_size) -> list[list[Point]]`; `default_cores() -> int`; `format_duration(seconds) -> "H:MM:SS"`; `confirm_batch(k, n_batches, batch, seconds_per_run, cores_for_estimate, assume_yes, ask=None) -> bool`; `ask_cores(default, preset=None, ask=None) -> int` (both resolve `ask = ask or input` at call time so tests can monkeypatch `builtins.input`); `make_progress(total, initial) -> tqdm`; `run_batch(batch, cores, runs_dir, raw_dir, timeout, progress=None, runner=None) -> (n_ok, n_failed, wall_seconds)`; `class TqdmHandler(logging.Handler)`.
- Prompts read through the injectable `ask` callable so tests never touch stdin; `run_batch` takes an injectable `runner` (defaults to `run_point` looked up at call time so tests can monkeypatch `sphere_sweep.run_point`).

- [ ] **Step 1: Write the failing tests (append)**

```python
import subprocess
import sys


class TestCommand:
    def test_build_command(self, parent):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        cmd = sw.build_command(point, "/out/runs", "/out/raw", 600)
        assert cmd[0] == sys.executable and cmd[1] == sw.SPHERE_REENTRY
        assert cmd[2:] == [
            "--velocity", "7.400000", "--altitude", "90.000000", "--temperature", "300.000000",
            "--diameter", "5.000000", "--flight-path-angle", "-0.750000", "--heading", sw.format_arg(point.heading_deg),
            "--lat", "-5.000000", "--lon", "149.000000", "--epoch", "2024-08-01T12:09:10",
            "--outdir", "/out/runs", "--raw-dir", "/out/raw", "--timeout", "600", "--quiet"]
        assert sw.build_command(point, "/o", "/r", 5, python="/usr/bin/python3")[0] == "/usr/bin/python3"


class TestRunPoint:
    def fake_subprocess(self, monkeypatch, runs_dir, returncode=0, stderr="", json_status="ok", raise_timeout=False):
        calls = []

        def fake_run(cmd, capture_output, text, timeout):
            calls.append(cmd)
            if raise_timeout:
                raise subprocess.TimeoutExpired(cmd, timeout)
            if json_status is not None:                       # imitate sphere_reentry writing its JSON
                write_run_json(runs_dir, self.current_name, json_status, results={})
            return subprocess.CompletedProcess(cmd, returncode, stdout="", stderr=stderr)

        monkeypatch.setattr(sw.subprocess, "run", fake_run)
        return calls

    def test_success(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        calls = self.fake_subprocess(monkeypatch, str(tmp_path / "runs"))
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "done" and point.returncode == 0 and point.wall_time_s >= 0.0
        assert point.stderr_tail is None
        assert calls[0][:2] == [sys.executable, sw.SPHERE_REENTRY]

    def test_nonzero_return_code_is_failed_with_stderr_tail(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        stderr = "\n".join("line %d" % i for i in range(30))
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), returncode=1, stderr=stderr, json_status="error")
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed" and point.returncode == 1
        assert point.stderr_tail.splitlines() == ["line %d" % i for i in range(10, 30)]

    def test_timeout_is_failed(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), raise_timeout=True, json_status=None)
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed" and point.returncode == -1 and "timed out" in point.stderr_tail

    def test_zero_return_code_but_json_not_ok_is_failed(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), returncode=0, json_status="error")
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed"

    def test_zero_return_code_without_json_is_failed(self, parent, tmp_path, monkeypatch):
        point = sw.build_points(parent, [5.0], [300.0], [7.4])[0]
        self.current_name = point.run_name
        self.fake_subprocess(monkeypatch, str(tmp_path / "runs"), returncode=0, json_status=None)
        sw.run_point(point, str(tmp_path / "runs"), str(tmp_path / "raw"), 600)
        assert point.status == "failed"


class TestBatchesAndPrompts:
    def test_split_batches(self):
        pts = list(range(7))
        assert sw.split_batches(pts, None) == [pts]
        assert sw.split_batches(pts, 3) == [[0, 1, 2], [3, 4, 5], [6]]
        assert sw.split_batches([], None) == [] and sw.split_batches([], 3) == []

    def test_format_duration_and_default_cores(self):
        assert sw.format_duration(3661) == "1:01:01" and sw.format_duration(0.4) == "0:00:00"
        assert 1 <= sw.default_cores() <= (os.cpu_count() or 1)

    def test_confirm_batch(self, parent, caplog):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0])
        answers = iter(["n"])
        assert sw.confirm_batch(1, 3, batch, 0.3, 7, False, ask=lambda prompt: next(answers)) is False
        answers = iter([""])
        assert sw.confirm_batch(1, 3, batch, 0.3, 7, False, ask=lambda prompt: next(answers)) is True
        answers = iter(["maybe", "y"])
        assert sw.confirm_batch(2, 3, batch, 0.3, 7, False, ask=lambda prompt: next(answers)) is True

        def never(prompt):
            raise AssertionError("must not prompt with --yes")
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.confirm_batch(3, 3, batch, 0.3, 7, True, ask=never) is True
        assert "batch 3/3: 2 runs" in caplog.text and batch[0].run_name in caplog.text

    def test_ask_cores(self):
        assert sw.ask_cores(7, preset=3) == 3
        with pytest.raises(ValueError):
            sw.ask_cores(7, preset=0)
        answers = iter(["abc", "99999", "2"])
        assert sw.ask_cores(7, ask=lambda prompt: next(answers)) == 2
        answers = iter([""])
        assert sw.ask_cores(5, ask=lambda prompt: next(answers)) == 5


class FakeProgress:
    def __init__(self):
        self.updates = 0
        self.postfixes = []

    def update(self, n):
        self.updates += n

    def set_postfix_str(self, s, refresh=True):
        self.postfixes.append(s)

    def close(self):
        pass


class TestRunBatch:
    def test_counts_and_progress(self, parent):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0, 2.0])

        def runner(point, runs_dir, raw_dir, timeout):
            point.status = "failed" if point.velocity_kms == 5.0 else "done"
            point.returncode = 1 if point.status == "failed" else 0
            return point

        progress = FakeProgress()
        n_ok, n_failed, wall = sw.run_batch(batch, 2, "/r", "/w", 600, progress=progress, runner=runner)
        assert (n_ok, n_failed) == (2, 1) and wall >= 0.0
        assert progress.updates == 3 and len(progress.postfixes) == 3
        assert "d=5mm T=300K v=7.40000km/s" in progress.postfixes

    def test_default_runner_is_looked_up_at_call_time(self, parent, monkeypatch):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4])
        seen = []

        def fake_run_point(point, runs_dir, raw_dir, timeout):
            seen.append(point.run_name)
            point.status = "done"
            return point

        monkeypatch.setattr(sw, "run_point", fake_run_point)
        assert sw.run_batch(batch, 1, "/r", "/w", 600)[:2] == (1, 0)
        assert seen == [batch[0].run_name]

    def test_keyboard_interrupt_propagates_after_shutdown(self, parent):
        batch = sw.build_points(parent, [5.0], [300.0], [7.4, 5.0])

        def runner(point, runs_dir, raw_dir, timeout):
            raise KeyboardInterrupt

        with pytest.raises(KeyboardInterrupt):
            sw.run_batch(batch, 1, "/r", "/w", 600, runner=runner)

    def test_make_progress_is_a_tqdm_bar(self):
        bar = sw.make_progress(total=10, initial=3)
        assert bar.n == 3 and bar.total == 10
        bar.close()
```

Note for `TestRunPoint`: the fake `subprocess.run` writes the run JSON for the point under test (`self.current_name`) to imitate script 1; the command line itself is covered by `TestCommand`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v -k "Command or RunPoint or BatchesAndPrompts or RunBatch"`
Expected: failures with `AttributeError: ... 'build_command'` etc.

- [ ] **Step 3: Implement execution (append to sphere_sweep.py)**

```python
# =============================================================================
# 5. Execution: one subprocess per point, batches, prompts, progress (spec 5.5)
# =============================================================================

def build_command(point, runs_dir, raw_dir, timeout, python=None):
    return [python or sys.executable, SPHERE_REENTRY,
            "--velocity", format_arg(point.velocity_kms),
            "--altitude", format_arg(point.altitude_km),
            "--temperature", format_arg(point.temperature_K),
            "--diameter", format_arg(point.diameter_mm),
            "--flight-path-angle", format_arg(point.flight_path_deg),
            "--heading", format_arg(point.heading_deg),
            "--lat", format_arg(point.lat_deg),
            "--lon", format_arg(point.lon_deg),
            "--epoch", point.epoch_utc,
            "--outdir", runs_dir,
            "--raw-dir", raw_dir,
            "--timeout", str(int(timeout)),
            "--quiet"]


def run_point(point, runs_dir, raw_dir, timeout):
    """Run sphere_reentry.py for one point; mark it done/failed from the exit code and its JSON."""
    cmd = build_command(point, runs_dir, raw_dir, timeout)
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 60)
        point.returncode = proc.returncode
        stderr = proc.stderr or ""
    except subprocess.TimeoutExpired:
        point.returncode = -1
        stderr = "sweep: subprocess timed out after {} s".format(timeout + 60)
    point.wall_time_s = round(time.time() - t0, 3)
    tail = "\n".join(stderr.splitlines()[-STDERR_TAIL_LINES:])
    point.stderr_tail = tail or None
    json_ok = _json_status(os.path.join(runs_dir, point.run_name + ".json")) == "ok"
    point.status = "done" if (point.returncode == 0 and json_ok) else "failed"
    return point


def split_batches(points, batch_size):
    points = list(points)
    if not points:
        return []
    if not batch_size:
        return [points]
    return [points[i:i + batch_size] for i in range(0, len(points), batch_size)]


def default_cores():
    return max(1, (os.cpu_count() or 2) - 1)


def format_duration(seconds):
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return "{:d}:{:02d}:{:02d}".format(hours, minutes, secs)


def confirm_batch(k, n_batches, batch, seconds_per_run, cores_for_estimate, assume_yes, ask=None):
    """Describe the batch, then ask 'Run this batch? [Y/n]' unless assume_yes."""
    ask = ask or input        # resolved at call time (tests monkeypatch builtins.input)
    estimate = len(batch) * seconds_per_run / max(1, cores_for_estimate)
    log.info("batch %d/%d: %d runs, %s .. %s, estimated %s at %d cores",
             k, n_batches, len(batch), batch[0].run_name, batch[-1].run_name,
             format_duration(estimate), cores_for_estimate)
    if assume_yes:
        return True
    while True:
        answer = ask("Run this batch? [Y/n] ").strip().lower()
        if answer in ("", "y", "yes"):
            log.info("user confirmed batch %d", k)
            return True
        if answer in ("n", "no", "q", "quit"):
            log.info("user declined batch %d", k)
            return False


def ask_cores(default, preset=None, ask=None):
    """Cores for the next batch: --cores preset, or a prompt with `default` (1..cpu_count)."""
    ask = ask or input        # resolved at call time (tests monkeypatch builtins.input)
    max_cores = os.cpu_count() or 1
    if preset is not None:
        if 1 <= preset <= max_cores:
            return preset
        raise ValueError("--cores must be between 1 and {}".format(max_cores))
    while True:
        answer = ask("Cores for this batch [{}]: ".format(default)).strip()
        if answer == "":
            return default
        try:
            cores = int(answer)
        except ValueError:
            cores = 0
        if 1 <= cores <= max_cores:
            return cores
        print("please enter an integer between 1 and {}".format(max_cores))


class TqdmHandler(logging.Handler):
    """Console log handler that writes through tqdm so the progress bar is not garbled."""

    def emit(self, record):
        try:
            tqdm.write(self.format(record))
        except Exception:  # noqa: BLE001
            self.handleError(record)


def make_progress(total, initial):
    return tqdm(total=total, initial=initial, desc="sweep", unit="run", dynamic_ncols=True)


def _run_and_mark(point, runs_dir, raw_dir, timeout, progress, runner):
    if progress is not None:
        progress.set_postfix_str("d={:g}mm T={:g}K v={:.5f}km/s".format(
            point.diameter_mm, point.temperature_K, point.velocity_kms), refresh=False)
    return runner(point, runs_dir, raw_dir, timeout)


def run_batch(batch, cores, runs_dir, raw_dir, timeout, progress=None, runner=None):
    """Run one batch on `cores` worker threads. Returns (n_ok, n_failed, wall_seconds).

    KeyboardInterrupt cancels queued points (they stay pending), lets the running
    subprocesses finish, then propagates.
    """
    runner = runner or run_point
    t0 = time.time()
    n_ok = n_failed = 0
    executor = ThreadPoolExecutor(max_workers=cores)
    try:
        futures = [executor.submit(_run_and_mark, p, runs_dir, raw_dir, timeout, progress, runner)
                   for p in batch]
        for future in as_completed(futures):
            point = future.result()
            if point.status == "done":
                n_ok += 1
            else:
                n_failed += 1
                log.warning("FAILED %s (rc=%s): %s", point.run_name, point.returncode, point.stderr_tail)
            if progress is not None:
                progress.update(1)
    except KeyboardInterrupt:
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    executor.shutdown(wait=True)
    return n_ok, n_failed, time.time() - t0
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: 43 passed.

- [ ] **Step 5: Commit**

```bash
git add sphere_sweep.py tests/test_sphere_sweep.py
git commit -m "Add batched subprocess execution with prompts, thread pool and tqdm progress

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 15: sphere_sweep command line — logging, dry run, batch loop, interrupt handling

**Files:**
- Modify: `sphere_sweep.py` (append)
- Modify: `tests/test_sphere_sweep.py` (append)

**Interfaces:**
- Consumes: everything in `sphere_sweep.py` so far.
- Produces (spec §5.1, §5.5, §5.7): `build_parser()`, `setup_logging(outdir, dry_run)`, `main(argv=None) -> int` (0 ok / nothing to do / user declined; 2 usage or parent-file error; 130 after Ctrl-C), `if __name__ == "__main__": sys.exit(main())`.

- [ ] **Step 1: Write the failing tests (append)**

```python
def fake_run_point_factory(fail_velocities=()):
    """A run_point stand-in that writes an ok JSON (or an error JSON) like sphere_reentry would."""
    def fake_run_point(point, runs_dir, raw_dir, timeout):
        if point.velocity_kms in fail_velocities:
            write_run_json(runs_dir, point.run_name, "error")
            point.status, point.returncode, point.stderr_tail = "failed", 1, "ERROR (error): boom"
        else:
            write_run_json(runs_dir, point.run_name, "ok", results={
                "max_temperature_K": 850.0, "final_mass_kg": 0.1, "final_mass_source": "history_file",
                "mass_loss_fraction": 0.5, "time_at_melting_temperature_s": 3.0, "final_velocity_kms": 1.0,
                "final_radius_mm": 20.0, "final_altitude_km": 0.0, "end_of_life_reason": "ground impact",
                "outcome": "survived"})
            point.status, point.returncode = "done", 0
        point.wall_time_s = 0.01
        return point
    return fake_run_point


SUBSET = ["--diameters", "5", "--temperatures", "300", "--velocities", "7.4,5.0,9.0"]


class TestMain:
    def test_dry_run_logs_v_min_and_writes_nothing(self, tmp_path, caplog):
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--dry-run"] + SUBSET)
        assert rc == 0
        assert "v_min = 0.050000 km/s at t = 1300.000 s, altitude 0.000 km (parent row 13)" in caplog.text
        assert "3 points" in caplog.text and "1 skipped" in caplog.text and "2 pending" in caplog.text
        assert "batch 1/1: 2 runs" in caplog.text and "dry run" in caplog.text
        assert not out.exists()

    def test_full_run_with_yes_and_cores(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "2"] + SUBSET)
        assert rc == 0
        manifest = sw.load_manifest(str(out / "sweep_manifest.json"))
        statuses = [p["status"] for p in manifest["points"]]
        assert sorted(statuses) == ["done", "done", "skipped"]
        assert manifest["parent"]["v_min_kms"] == 0.05 and manifest["settings"]["batch_size"] is None
        assert manifest["grid"]["velocities_kms"] == [7.4, 5.0, 9.0]
        with open(out / "sweep_summary.csv", newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 3 and sum(r["outcome"] == "survived" for r in rows) == 2
        log_text = (out / "sweep.log").read_text()
        assert "v_min = 0.050000 km/s" in log_text and "batch 1/1: 2 runs on 2 cores" in log_text
        assert "sweep finished: 2 ok, 0 failed, 0 pending" in caplog.text

    def test_second_invocation_has_nothing_to_do(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        assert "nothing to do" in caplog.text

    def test_failed_runs_are_retried_only_with_retry_failed(self, tmp_path, monkeypatch, caplog):
        out = tmp_path / "out"
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory(fail_velocities=(5.0,)))
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        assert "1 ok, 1 failed" in caplog.text and "FAILED" in caplog.text
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        caplog.clear()
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET) == 0
        assert "nothing to do" in caplog.text
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1", "--retry-failed"] + SUBSET) == 0
        manifest = sw.load_manifest(str(out / "sweep_manifest.json"))
        assert sorted(p["status"] for p in manifest["points"]) == ["done", "done", "skipped"]

    def test_batches_prompt_for_confirmation_and_cores(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        answers = iter(["", "3", "y", ""])          # batch 1: confirm, 3 cores; batch 2: confirm, default cores
        monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out), "--batch-size", "1"] + SUBSET)
        assert rc == 0
        assert "batch 1/2: 1 runs on 3 cores" in caplog.text
        assert "batch 2/2: 1 runs on {} cores".format(sw.default_cores()) in caplog.text
        assert sw.load_manifest(str(out / "sweep_manifest.json"))["settings"]["batch_size"] == 1

    def test_user_declines_first_batch(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        monkeypatch.setattr("builtins.input", lambda prompt="": "n")
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            rc = sw.main(["--parent", MINI, "--outdir", str(out)] + SUBSET)
        assert rc == 0 and "stopped before batch 1" in caplog.text
        statuses = [p["status"] for p in sw.load_manifest(str(out / "sweep_manifest.json"))["points"]]
        assert sorted(statuses) == ["pending", "pending", "skipped"]

    def test_yes_without_cores_still_asks_for_cores(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        asked = []
        monkeypatch.setattr("builtins.input", lambda prompt="": asked.append(prompt) or "2")
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes"] + SUBSET) == 0
        assert asked == ["Cores for this batch [{}]: ".format(sw.default_cores())]
        assert "on 2 cores" in caplog.text

    def test_keyboard_interrupt_saves_state_and_returns_130(self, tmp_path, monkeypatch):
        def interrupt(*args, **kwargs):
            raise KeyboardInterrupt
        monkeypatch.setattr(sw, "run_batch", interrupt)
        out = tmp_path / "out"
        rc = sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1"] + SUBSET)
        assert rc == 130
        assert (out / "sweep_manifest.json").is_file() and (out / "sweep_summary.csv").is_file()

    @pytest.mark.parametrize("extra", [["--batch-size", "0"], ["--limit", "0"], ["--cores", "0"], ["--cores", "100000"]])
    def test_bad_arguments_exit_2(self, tmp_path, extra):
        with pytest.raises(SystemExit) as exc:
            sw.main(["--parent", MINI, "--outdir", str(tmp_path / "out")] + extra)
        assert exc.value.code == 2

    def test_missing_parent_returns_2(self, tmp_path, caplog):
        with caplog.at_level("ERROR", logger="sweep"):
            rc = sw.main(["--parent", str(tmp_path / "nope.json"), "--outdir", str(tmp_path / "out"), "--dry-run"])
        assert rc == 2 and "cannot load parent" in caplog.text

    def test_limit(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(sw, "run_point", fake_run_point_factory())
        out = tmp_path / "out"
        with caplog.at_level("INFO", logger="sweep"):
            assert sw.main(["--parent", MINI, "--outdir", str(out), "--yes", "--cores", "1", "--limit", "2"] + SUBSET) == 0
        assert len(sw.load_manifest(str(out / "sweep_manifest.json"))["points"]) == 2
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v -k TestMain`
Expected: failures with `AttributeError: ... 'main'`.

- [ ] **Step 3: Implement the CLI (append to sphere_sweep.py)**

```python
# =============================================================================
# 6. Command line (spec 5.1, 5.5, 5.7)
# =============================================================================

def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parent", default=DEFAULT_PARENT, help="parent dmf_output.json (default %(default)s)")
    p.add_argument("--parent-object", default=DEFAULT_PARENT_OBJECT, help="object whose trajectory is used")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="sweep output directory (default %(default)s)")
    p.add_argument("--batch-size", type=int, default=None, help="runs per batch (default: all pending runs in one batch)")
    p.add_argument("--cores", type=int, default=None, help="cores for every batch (default: ask before each batch)")
    p.add_argument("--yes", action="store_true", help="do not ask for confirmation before a batch")
    p.add_argument("--diameters", type=parse_float_list, default=None, help="comma-separated diameters [mm]")
    p.add_argument("--temperatures", type=parse_float_list, default=None, help="comma-separated temperatures [K]")
    p.add_argument("--velocities", type=parse_float_list, default=None, help="comma-separated velocities [km/s]")
    p.add_argument("--limit", type=int, default=None, help="only the first N points of the matrix")
    p.add_argument("--dry-run", action="store_true", help="show the matrix and batch plan; run nothing")
    p.add_argument("--force", action="store_true", help="re-run every point, even completed ones")
    p.add_argument("--retry-failed", action="store_true", help="re-run points that failed before")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_S, help="per-run SESAM timeout [s]")
    return p


def setup_logging(outdir, dry_run):
    log.setLevel(logging.INFO)
    for handler in list(log.handlers):
        handler.close()
    log.handlers.clear()
    log.propagate = True
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    console = TqdmHandler()
    console.setFormatter(fmt)
    log.addHandler(console)
    if not dry_run:
        os.makedirs(outdir, exist_ok=True)
        file_handler = logging.FileHandler(os.path.join(outdir, "sweep.log"))
        file_handler.setFormatter(fmt)
        log.addHandler(file_handler)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    max_cores = os.cpu_count() or 1
    if args.batch_size is not None and args.batch_size < 1:
        parser.error("--batch-size must be >= 1")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be >= 1")
    if args.cores is not None and not 1 <= args.cores <= max_cores:
        parser.error("--cores must be between 1 and {}".format(max_cores))

    runs_dir = os.path.join(args.outdir, "runs")
    raw_dir = os.path.join(args.outdir, "raw")
    setup_logging(args.outdir, args.dry_run)
    log.info("sphere_sweep start: parent=%s object=%r outdir=%s batch_size=%s cores=%s yes=%s force=%s retry_failed=%s",
             args.parent, args.parent_object, args.outdir, args.batch_size, args.cores, args.yes,
             args.force, args.retry_failed)

    # -- parent trajectory and the velocity grid's lower bound ------------------------
    try:
        parent = load_parent(args.parent, args.parent_object)
    except (OSError, ValueError, KeyError) as exc:
        log.error("cannot load parent trajectory from %s: %s", args.parent, exc)
        return 2
    log.info("parent %s (sha256 %s) object %r epoch %s; descending branch rows %d..%d, v_max %.6f km/s",
             parent.path, parent.sha256[:12], parent.object_name,
             parent.epoch.strftime(EPOCH_ARG_FORMAT), parent.branch_start, len(parent.rows) - 1, parent.v_max)
    log.info("%s", parent.describe_v_min())

    diameters = args.diameters or DIAMETERS_MM
    temperatures = args.temperatures or TEMPERATURES_K
    velocities = args.velocities or velocity_grid(parent.v_min)
    if args.velocities is None:
        log.info("velocity grid: linspace(%.6f, %.6f, %d) km/s (lower bound = parent v_min)",
                 V_TOP_KMS, parent.v_min, N_VELOCITIES)
    grid = {"diameters_mm": diameters, "temperatures_K": temperatures, "velocities_kms": velocities}

    # -- matrix and resume ------------------------------------------------------------
    points = build_points(parent, diameters, temperatures, velocities, limit=args.limit)
    if not args.dry_run:
        os.makedirs(runs_dir, exist_ok=True)
        os.makedirs(raw_dir, exist_ok=True)
    apply_resume(points, runs_dir, force=args.force, retry_failed=args.retry_failed)
    skipped = [p for p in points if p.status == "skipped"]
    pending = [p for p in points if p.status == "pending"]
    n_done = sum(p.status == "done" for p in points)
    n_failed_before = sum(p.status == "failed" for p in points)
    batches = split_batches(pending, args.batch_size)
    log.info("matrix: %d points = %d diameters x %d temperatures x %d velocities%s; "
             "%d skipped, %d done, %d failed (not retried), %d pending in %d batch(es)",
             len(points), len(diameters), len(temperatures), len(velocities),
             " (limited to {})".format(args.limit) if args.limit else "",
             len(skipped), n_done, n_failed_before, len(pending), len(batches))
    for velocity in sorted({p.velocity_kms for p in skipped}, reverse=True):
        count = sum(p.velocity_kms == velocity for p in skipped)
        log.info("skipped v=%.6f km/s for %d point(s): parent never reaches velocity (range %.6f..%.6f)",
                 velocity, count, parent.v_min, parent.v_max)

    if args.dry_run:
        for k, batch in enumerate(batches, 1):
            log.info("batch %d/%d: %d runs, %s .. %s", k, len(batches), len(batch),
                     batch[0].run_name, batch[-1].run_name)
        log.info("dry run: nothing executed, nothing written")
        return 0

    manifest_path = os.path.join(args.outdir, "sweep_manifest.json")
    summary_path = os.path.join(args.outdir, "sweep_summary.csv")
    settings = {"timeout_s": args.timeout, "batch_size": args.batch_size,
                "sphere_reentry_version": sr.SCRIPT_VERSION}
    created = utc_now_str()

    def save():
        write_manifest(manifest_path, manifest_document(parent, grid, settings, points, created))
        write_summary(summary_path, points, runs_dir)

    save()
    if not pending:
        log.info("nothing to do: %d done, %d failed, %d skipped (use --retry-failed / --force to re-run)",
                 n_done, n_failed_before, len(skipped))
        return 0

    # -- batches ------------------------------------------------------------------------
    completed = 0
    total_ok = total_failed = 0
    total_wall = 0.0
    seconds_per_run = SECONDS_PER_RUN_ESTIMATE      # per run per core
    exit_code = 0
    try:
        for k, batch in enumerate(batches, 1):
            if not confirm_batch(k, len(batches), batch, seconds_per_run,
                                 args.cores or default_cores(), args.yes):
                log.info("stopped before batch %d by user", k)
                break
            cores = ask_cores(default_cores(), preset=args.cores)
            log.info("batch %d/%d: %d runs on %d cores", k, len(batches), len(batch), cores)
            progress = make_progress(total=len(pending), initial=completed)
            try:
                n_ok, n_failed, wall = run_batch(batch, cores, runs_dir, raw_dir, args.timeout, progress)
            finally:
                progress.close()
            completed += n_ok + n_failed
            total_ok += n_ok
            total_failed += n_failed
            total_wall += wall
            if n_ok + n_failed:
                seconds_per_run = wall * cores / (n_ok + n_failed)
            log.info("batch %d/%d done: %d ok, %d failed, %s (%.2f s/run/core)",
                     k, len(batches), n_ok, n_failed, format_duration(wall), seconds_per_run)
            save()
    except KeyboardInterrupt:
        log.warning("interrupted: saving manifest and summary; re-run to resume")
        exit_code = 130
    save()
    still_pending = sum(p.status == "pending" for p in points)
    log.info("sweep finished: %d ok, %d failed, %d pending, total run time %s",
             total_ok, total_failed, still_pending, format_duration(total_wall))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `"$PY" -m pytest tests/test_sphere_sweep.py -v`
Expected: 57 passed. Then check the real parent file end to end without running anything:
`"$PY" sphere_sweep.py --dry-run 2>&1 | tail -8` — expect `v_min = 0.028000 km/s at t = 3870.071 s, altitude 0.012 km (parent row 784)`, `92000 points`, `0 skipped`, one batch of 92000 runs, and `dry run: nothing executed`.

- [ ] **Step 5: Commit**

```bash
git add sphere_sweep.py tests/test_sphere_sweep.py
git commit -m "Add the sphere_sweep command line: logging, dry run, confirmed batches, resume

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 16: Sweep integration test against real DRAMA and the real parent file

**Files:**
- Modify: `tests/test_integration_drama.py` (append)

**Interfaces:**
- Consumes: `sphere_sweep.py` CLI, the real parent file `Generic_Satellite Reentry/output/dmf_output.json` (skipped when absent).

- [ ] **Step 1: Write the integration test (append)**

```python
REAL_PARENT = os.path.join(REPO_ROOT, "Generic_Satellite Reentry", "output", "dmf_output.json")


def sweep(tmp_out, *extra):
    cmd = [PY, os.path.join(REPO_ROOT, "sphere_sweep.py"), "--parent", REAL_PARENT, "--outdir", str(tmp_out),
           "--diameters", "50", "--temperatures", "300", "--velocities", "7.5,0.5,9.0"] + list(extra)
    return subprocess.run(cmd, capture_output=True, text=True)


@pytest.mark.skipif(not os.path.isfile(REAL_PARENT), reason="parent dmf_output.json not present")
def test_sweep_end_to_end_with_real_parent(tmp_path):
    out = tmp_path / "out"
    proc = sweep(out, "--yes", "--cores", "2")
    assert proc.returncode == 0, proc.stderr
    manifest = json.load(open(out / "sweep_manifest.json"))
    assert manifest["parent"]["v_min_kms"] == pytest.approx(0.028, abs=1e-6)
    by_velocity = {p["velocity_kms"]: p for p in manifest["points"]}
    assert by_velocity[9.0]["status"] == "skipped"
    assert by_velocity[7.5]["status"] == "done" and by_velocity[0.5]["status"] == "done"
    assert by_velocity[7.5]["altitude_km"] == pytest.approx(77.5, abs=0.5)
    assert by_velocity[0.5]["altitude_km"] == pytest.approx(39.9, abs=0.5)
    with open(out / "sweep_summary.csv", newline="") as fh:
        rows = {float(r["initial_velocity_kms"]): r for r in csv.DictReader(fh)}
    assert rows[7.5]["outcome"] == "demised" and rows[0.5]["outcome"] == "survived"
    assert rows[9.0]["status"] == "skipped" and rows[9.0]["outcome"] == ""
    assert float(rows[0.5]["final_radius_mm"]) == pytest.approx(25.0, abs=1e-3)
    for p in manifest["points"]:
        if p["status"] == "done":
            assert (out / "runs" / (p["run_name"] + ".csv")).is_file()
            assert (out / "runs" / (p["run_name"] + ".json")).is_file()
            assert not (out / "raw" / p["run_name"]).exists()
    log_text = (out / "sweep.log").read_text()
    assert "v_min = 0.028000 km/s" in log_text and "batch 1/1: 2 runs on 2 cores" in log_text

    again = sweep(out, "--yes", "--cores", "2")
    assert again.returncode == 0 and "nothing to do" in (again.stdout + again.stderr)

    fresh = tmp_path / "fresh"
    dry = sweep(fresh, "--dry-run")
    assert dry.returncode == 0 and not fresh.exists()
    assert "v_min = 0.028000 km/s at t = 3870.071 s" in (dry.stdout + dry.stderr)
```

- [ ] **Step 2: Run the integration tests**

Run: `"$PY" -m pytest tests/test_integration_drama.py -v`
Expected: 4 passed (the sweep test starts three sphere_reentry subprocesses; a few seconds).

- [ ] **Step 3: Commit**

```bash
git add tests/test_integration_drama.py
git commit -m "Add real-DRAMA end-to-end test of sphere_sweep

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 17: README and final verification

**Files:**
- Create: `README.md`

- [ ] **Step 1: Write README.md**

```markdown
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
`sweep_manifest.json`, `sweep_summary.csv` (one row per matrix point), `sweep.log`.

The velocity grid's lower bound is the parent's minimum velocity, read from the parent file and
logged as `v_min = ... km/s at t = ... s, altitude ... km (parent row N)`.

## Tests

```bash
"$PY" -m pytest -m "not drama"     # unit tests, no DRAMA needed (real SESAM outputs in tests/fixtures/)
"$PY" -m pytest                    # also the integration tests that run SESAM
```

To refresh a fixture see `tests/fixtures/README.md`.
```

- [ ] **Step 2: Run the complete test suite and a real dry run**

Run: `"$PY" -m pytest -v 2>&1 | tail -15`
Expected: every test passes (unit + drama).

Run: `"$PY" sphere_sweep.py --dry-run 2>&1 | tail -6`
Expected: the v_min line (`0.028000 km/s ... parent row 784`), `92000 points`, `0 skipped`, `dry run: nothing executed, nothing written`.

Run: `git status --short`
Expected: only `README.md` untracked; no `sphere_sweep_output/` entries (git-ignored).

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "Add README with environment, usage and test instructions

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Notes for executors

- The spike-derived facts in spec §2 are load-bearing: `config=[cfg]` as a list, the geodetic element order, `energyThreshold = 1e-9`, the 3-decimal mass column, `EVENT end` masses, the `ballooning` / `uncritical` / `ground impact` log phrases.
- Never run pytest or the scripts with the system `python3`; only `drama_env` has pyDRAMA, numpy, tqdm and pytest.
- The DRAMA GUI rewrites `Generic_Satellite Reentry/` on every GUI run — never write into it; only read `output/dmf_output.json` and (once, Task 1) copy the two fap files out of `data/`.
- If SESAM output ever gains columns, `read_sara_table` skips rows whose width differs from the column list: a sudden row count of 0 is the symptom.
