# Sphere Re-entry Trajectory Model (Step 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Python package `reentry_model` that integrates the trajectory of a constant-mass sphere from (V₀, h₀, γ₀) through the free-molecular, transitional and continuum regimes the way ESA DRAMA's SESAM does, and verifies it against four SESAM reference runs with overlaid V(t)/h(t) plots and error metrics.

**Architecture:** ECEF-Cartesian 3-DOF integration (rotating frame, J2 gravity, DOP853) with pluggable atmosphere (NRLMSISE-00 via pymsis / US76 table / SESAM replay), SESAM's sphere drag tables blended by a Knudsen bridging function measured from SESAM's own output, a `Body` hook for the later thermal model, and a comparison module that samples the model at the reference's time stamps. Every physics module is a small file with a plain-array interface and its own tests.

**Tech Stack:** Python 3.12 in the conda env `drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`), numpy 2.x, scipy (solve_ivp DOP853), pymsis 0.13.0 (NRLMSISE-00), matplotlib (Agg), pytest.

**Spec:** `docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md` — read it first; the facts it relies on are in `Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md`.

## Global Constraints

- Interpreter for everything: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` (never the system python). Run tests as `"$PY" -m pytest -m "not drama" -q`.
- SI units inside the package (m, s, kg, K, rad); the CLI and CSV use the wrapper's units (km, km/s, deg, mm). Angles in the package are radians; conversion happens only in `cli.py`, `sesam_io.py` and the CSV writer.
- No DRAMA install is needed by the package or its tests: the ATDB tables and the US76 table are copied into `reentry_model/data/`, and the reference runs into `data/reference_runs/`.
- Git-ignored output root: `reentry_model_output/`. Never write elsewhere by default.
- Follow the repo's conventions: module docstrings, `argparse`, exit codes 0/1/2, tests under `tests/`, one commit per task, commit messages in the imperative like the existing history, ending with `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
- Keep the wrapper (`sphere_reentry.py`, `sphere_sweep.py`) untouched; the model is a sibling package.

---

## File structure

```
reentry_model/__init__.py          package marker, __version__
reentry_model/__main__.py          python -m reentry_model -> cli.main()
reentry_model/constants.py         WGS-84, μ, J2, J4, ω, g0, k_B, u, air m̄, γ, hard-sphere d
reentry_model/earth.py             geodetic<->ECEF, ENU basis, (V,γ,ψ)<->velocity vector, gravity, great-circle
reentry_model/fap.py               DRAMA fap_day.dat parsing -> SolarIndices (F10.7, F10.7a, Ap)
reentry_model/aero.py              mean free path, Kn, Ma, sphere drag tables, bridging functions, C_D
reentry_model/sesam_io.py          Reference: a wrapper run CSV+JSON as arrays + initial state
reentry_model/atmosphere.py        Freestream; US76TableAtmosphere, NRLMSISE00Atmosphere, ReplayAtmosphere; winds
reentry_model/body.py              Body protocol, ConstantBody (Step 2 hook)
reentry_model/trajectory.py        InitialState, Settings, Simulator, History, CSV/JSON writers
reentry_model/compare.py           metrics + six plots, model vs reference
reentry_model/cli.py               `run` and `compare` subcommands
reentry_model/data/atdb_sphere.json                 the six-Mach C_D tables (from ATDB_SPHERE.nc)
reentry_model/data/us76_static_environment.csv      copy of DRAMA's StaticEnvironmentData.csv
data/reference_runs/<4 run names>.csv|.json         SESAM reference histories (committed)
tests/test_reentry_model_earth.py ... _fap, _aero, _sesam_io, _atmosphere, _trajectory, _compare, _cli
tests/test_reentry_model_reference.py               integration test, marker `reference`
```

---

### Task 1: Package scaffold, data files, dependencies

**Files:**
- Create: `reentry_model/__init__.py`, `reentry_model/constants.py`, `reentry_model/data/atdb_sphere.json`, `reentry_model/data/us76_static_environment.csv`, `data/reference_runs/` (8 files), `tests/test_reentry_model_data.py`
- Modify: `pytest.ini` (marker), `.gitignore` (output dir)

**Interfaces:**
- Produces: constants `WGS84_A, WGS84_F, WGS84_E2, MU_EARTH, J2, J4, OMEGA_EARTH, G0, K_BOLTZMANN, ATOMIC_MASS_UNIT, M_BAR_AIR, GAMMA_AIR, HARD_SPHERE_DIAMETER` (floats, SI); `reentry_model.DATA_DIR` (str); the two data files; `data/reference_runs/` with the four `..._msis.csv/.json` and `..._msis_nowind.csv/.json` pairs.

- [ ] **Step 1: Install the two new dependencies into drama_env**

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
"$PY" -m pip install "pymsis==0.13.0" scipy
"$PY" -c "import pymsis, scipy; print(pymsis.__version__, scipy.__version__)"
```
Expected: prints `0.13.0` and a scipy version (1.14 or later).

- [ ] **Step 2: Copy the reference runs and the two DRAMA tables into the repo**

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution"
mkdir -p data/reference_runs reentry_model/data
cp sphere_sweep_output/reference_AA7075_nomelt/runs/*_msis.csv sphere_sweep_output/reference_AA7075_nomelt/runs/*_msis.json \
   sphere_sweep_output/reference_AA7075_nomelt/runs/*_msis_nowind.csv sphere_sweep_output/reference_AA7075_nomelt/runs/*_msis_nowind.json data/reference_runs/
ls data/reference_runs            # expect 8 files: 2 spheres x {csv,json} x {winds on, _nowind}
{ echo "# Copy of /Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/StaticEnvironmentData.csv (DRAMA 4.1.4, HTG), verbatim below this line. It is the US Standard Atmosphere 1976 with a fixed wind profile. Columns: alt[m], pressure[Pa], density[kg/m3], temperature[K], wind north[m/s], wind east[m/s], wind down[m/s], gamma, atomic oxygen fraction."; cat /Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/StaticEnvironmentData.csv; } > reentry_model/data/us76_static_environment.csv
```

- [ ] **Step 3: Write the ATDB JSON (values from `ncdump ATDB_SPHERE.nc`, see the facts note §1)**

Create `reentry_model/data/atdb_sphere.json`:
```json
{
  "_provenance": "Values of /Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_SPHERE.nc (ESA DRAMA 3.0 aerothermal database for the primitive SPHERE, Hyperschall Technologie Goettingen GmbH, file date 2019-04-11), read with ncdump on 2026-09-17. cd_free_molecular / cd_continuum are drag coefficients; heat_flux_factor_* are surface-average / stagnation heat-flux ratios (Step 2).",
  "mach": [5.0, 10.0, 15.0, 20.0, 25.0, 30.0],
  "cd_free_molecular": [2.360635, 2.147897, 2.089229, 2.062249, 2.046842, 2.036933],
  "cd_continuum": [0.898818, 0.910198, 0.912322, 0.913067, 0.913411, 0.913599],
  "heat_flux_factor_free_molecular": [0.270686, 0.254836, 0.251899, 0.250872, 0.2504, 0.250149],
  "heat_flux_factor_continuum": [0.27471, 0.27471, 0.27471, 0.27471, 0.27471, 0.27471]
}
```

- [ ] **Step 4: Write the failing data test**

Create `tests/test_reentry_model_data.py`:
```python
"""The package's bundled data and the committed SESAM reference runs."""
import csv
import json
import os

import reentry_model
from reentry_model import constants as c

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(REPO_ROOT, "data", "reference_runs")

REFERENCE_NAMES = [
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis",
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind",
]


def test_atdb_tables_have_the_six_mach_points():
    with open(os.path.join(reentry_model.DATA_DIR, "atdb_sphere.json")) as fh:
        atdb = json.load(fh)
    assert atdb["mach"] == [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    assert atdb["cd_free_molecular"][0] == 2.360635 and atdb["cd_continuum"][-1] == 0.913599
    assert len(atdb["heat_flux_factor_continuum"]) == 6


def test_us76_table_is_the_standard_atmosphere():
    path = os.path.join(reentry_model.DATA_DIR, "us76_static_environment.csv")
    rows = [r for r in csv.reader(open(path)) if r and not r[0].startswith("#")]
    assert len(rows) == 1501                      # 0..150 km every 100 m
    first = [float(v) for v in rows[0]]
    assert first[0] == 0.0 and abs(first[1] - 101325.0) < 1.0 and abs(first[2] - 1.225) < 1e-3
    at70 = [float(v) for v in rows[700]]
    assert at70[0] == 70000.0 and abs(at70[2] - 8.283e-5) / 8.283e-5 < 2e-3   # US76: 8.283e-5 kg/m3


def test_reference_runs_are_committed_and_consistent():
    for name in REFERENCE_NAMES:
        csv_path = os.path.join(REFERENCE_DIR, name + ".csv")
        json_path = os.path.join(REFERENCE_DIR, name + ".json")
        assert os.path.isfile(csv_path) and os.path.isfile(json_path), name
        doc = json.load(open(json_path))
        assert doc["run_name"] == name and doc["status"] == "ok"
        assert doc["inputs"]["material"] == "AA7075_nomelt" and doc["inputs"]["atmosphere"] == "nrlmsise"
        rows = list(csv.DictReader(open(csv_path)))
        assert len(rows) > 100 and float(rows[-1]["altitude_km"]) == 0.0
        assert abs(float(rows[-1]["mass_kg"]) - doc["inputs"]["initial_mass_kg"]) < 1e-3   # no mass loss


def test_constants():
    assert c.WGS84_A == 6378137.0 and abs(c.WGS84_F - 1 / 298.257223563) < 1e-15
    assert c.MU_EARTH == 3.986004418e14 and c.OMEGA_EARTH == 7.2921159e-5
    assert abs(c.M_BAR_AIR / c.ATOMIC_MASS_UNIT - 28.9644) < 1e-6
    assert c.HARD_SPHERE_DIAMETER == 3.65e-10 and c.GAMMA_AIR == 1.4
```

- [ ] **Step 5: Run it to verify it fails**

Run: `"$PY" -m pytest tests/test_reentry_model_data.py -q`
Expected: FAIL / ERROR with `ModuleNotFoundError: No module named 'reentry_model'`

- [ ] **Step 6: Create the package marker and constants**

`reentry_model/__init__.py`:
```python
"""First-principles re-entry model of a solid sphere, Step 1: the trajectory (spec:
docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md)."""
import os

__version__ = "0.1.0"
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
```

`reentry_model/constants.py`:
```python
"""Physical, geodetic and gas constants (SI)."""
# WGS-84 ellipsoid
WGS84_A = 6378137.0                       # semi-major axis [m]
WGS84_F = 1.0 / 298.257223563             # flattening
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)      # first eccentricity squared
# gravity field
MU_EARTH = 3.986004418e14                 # [m^3/s^2]
J2 = 1.08262668e-3
J4 = -1.61962159e-6
OMEGA_EARTH = 7.2921159e-5                # Earth rotation rate [rad/s]
G0 = 9.80665                              # standard gravity [m/s^2], for load factors
# gas
K_BOLTZMANN = 1.380649e-23                # [J/K]
ATOMIC_MASS_UNIT = 1.66053906660e-27      # [kg]
M_BAR_AIR = 28.9644 * ATOMIC_MASS_UNIT    # mean molecular mass of air below ~100 km [kg]
GAMMA_AIR = 1.4
HARD_SPHERE_DIAMETER = 3.65e-10           # [m]; reproduces SESAM's Knudsen column (facts note s.2)
```

- [ ] **Step 7: Register the test marker and ignore the output root**

Append to `pytest.ini` under `markers =`:
```
    reference: runs the full model against the committed SESAM reference runs (seconds, no DRAMA needed)
```
Append to `.gitignore`:
```

# physics-model outputs (histories, comparison plots)
reentry_model_output/
```

- [ ] **Step 8: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_data.py -q`
Expected: 4 passed

- [ ] **Step 9: Commit**

```bash
git add reentry_model/__init__.py reentry_model/constants.py reentry_model/data data/reference_runs tests/test_reentry_model_data.py pytest.ini .gitignore
git commit -m "Scaffold reentry_model: constants, ATDB and US76 tables, committed SESAM reference runs

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Geodesy, local frame and gravity (`earth.py`)

**Files:**
- Create: `reentry_model/earth.py`, `tests/test_reentry_model_earth.py`

**Interfaces:**
- Consumes: `constants`.
- Produces (all SI, radians, numpy arrays of length 3 for vectors):
  - `geodetic_to_ecef(h, lat, lon) -> np.ndarray`
  - `ecef_to_geodetic(r) -> (h, lat, lon)`
  - `enu_basis(lat, lon) -> (e, n, u)` unit vectors in ECEF
  - `velocity_from_flight_angles(V, gamma, heading, lat, lon) -> np.ndarray` (heading clockwise from north)
  - `flight_angles(v, lat, lon) -> (V, gamma, heading)` with heading in [0, 2π)
  - `gravity(r, model="j2") -> np.ndarray` for model in {"point", "j2", "j2j4"}
  - `great_circle_distance(lat1, lon1, lat2, lon2, radius=WGS84_A) -> float`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_earth.py`:
```python
"""earth.py: WGS-84 conversions, local frame, gravity."""
import math

import numpy as np
import pytest

from reentry_model import constants as c
from reentry_model import earth

R100 = dict(h=77500.133, lat=math.radians(29.546067), lon=math.radians(-82.134333))
R50 = dict(h=115000.0, lat=math.radians(29.546067), lon=math.radians(-82.134333))


@pytest.mark.parametrize("state", [R100, R50, dict(h=0.0, lat=0.0, lon=0.0), dict(h=1e5, lat=math.radians(-60.0), lon=math.radians(170.0))])
def test_geodetic_ecef_round_trip(state):
    r = earth.geodetic_to_ecef(state["h"], state["lat"], state["lon"])
    h, lat, lon = earth.ecef_to_geodetic(r)
    assert abs(h - state["h"]) < 1e-6
    assert abs(lat - state["lat"]) < 1e-12 and abs(lon - state["lon"]) < 1e-12


def test_equator_and_pole_radii():
    assert np.linalg.norm(earth.geodetic_to_ecef(0.0, 0.0, 0.0)) == pytest.approx(c.WGS84_A)
    b = c.WGS84_A * (1.0 - c.WGS84_F)
    assert earth.geodetic_to_ecef(0.0, math.pi / 2, 0.0)[2] == pytest.approx(b)


def test_enu_basis_is_orthonormal_and_up_is_radial_at_equator():
    e, n, u = earth.enu_basis(R100["lat"], R100["lon"])
    for a in (e, n, u):
        assert np.linalg.norm(a) == pytest.approx(1.0)
    assert abs(e @ n) < 1e-12 and abs(e @ u) < 1e-12 and abs(n @ u) < 1e-12
    assert np.allclose(np.cross(e, n), u)
    e0, n0, u0 = earth.enu_basis(0.0, 0.0)
    assert np.allclose(u0, [1, 0, 0]) and np.allclose(n0, [0, 0, 1]) and np.allclose(e0, [0, 1, 0])


def test_velocity_from_flight_angles_round_trip():
    V, gamma, heading = 7500.0, math.radians(-0.959331), math.radians(347.168296)
    v = earth.velocity_from_flight_angles(V, gamma, heading, R100["lat"], R100["lon"])
    assert np.linalg.norm(v) == pytest.approx(V)
    V2, gamma2, heading2 = earth.flight_angles(v, R100["lat"], R100["lon"])
    assert V2 == pytest.approx(V) and gamma2 == pytest.approx(gamma, abs=1e-12)
    assert heading2 == pytest.approx(heading, abs=1e-12)


def test_heading_convention_north_and_east():
    e, n, u = earth.enu_basis(R100["lat"], R100["lon"])
    north = earth.velocity_from_flight_angles(100.0, 0.0, 0.0, R100["lat"], R100["lon"])
    east = earth.velocity_from_flight_angles(100.0, 0.0, math.pi / 2, R100["lat"], R100["lon"])
    assert north @ n == pytest.approx(100.0) and abs(north @ e) < 1e-9
    assert east @ e == pytest.approx(100.0) and abs(east @ n) < 1e-9


def test_gravity_point_mass_and_j2_at_equator_and_pole():
    r_eq = earth.geodetic_to_ecef(0.0, 0.0, 0.0)
    g_point = earth.gravity(r_eq, "point")
    assert np.linalg.norm(g_point) == pytest.approx(c.MU_EARTH / c.WGS84_A ** 2)
    g_j2 = earth.gravity(r_eq, "j2")
    assert np.linalg.norm(g_j2) == pytest.approx(c.MU_EARTH / c.WGS84_A ** 2 * (1.0 + 1.5 * c.J2), rel=1e-9)
    assert g_j2 @ r_eq < 0                                  # points inward
    r_pole = earth.geodetic_to_ecef(0.0, math.pi / 2, 0.0)
    g_pole = earth.gravity(r_pole, "j2")
    b = np.linalg.norm(r_pole)
    expected = c.MU_EARTH / b ** 2 * (1.0 - 3.0 * c.J2 * (c.WGS84_A / b) ** 2)
    assert np.linalg.norm(g_pole) == pytest.approx(expected, rel=1e-9)
    assert abs(np.linalg.norm(g_pole) - 9.832) < 0.002       # textbook polar gravity


def test_gravity_j4_is_a_small_correction():
    r = earth.geodetic_to_ecef(77500.0, R100["lat"], R100["lon"])
    g2, g4 = earth.gravity(r, "j2"), earth.gravity(r, "j2j4")
    assert np.linalg.norm(g4 - g2) < 3e-5 * np.linalg.norm(g2)


def test_gravity_rejects_unknown_model():
    with pytest.raises(ValueError):
        earth.gravity(np.array([c.WGS84_A, 0.0, 0.0]), "j6")


def test_great_circle_distance():
    quarter = earth.great_circle_distance(0.0, 0.0, 0.0, math.pi / 2)
    assert quarter == pytest.approx(math.pi / 2 * c.WGS84_A)
    assert earth.great_circle_distance(R100["lat"], R100["lon"], R100["lat"], R100["lon"]) == 0.0
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_earth.py -q`
Expected: ERROR `cannot import name 'earth'`

- [ ] **Step 3: Implement `earth.py`**

```python
"""WGS-84 geodesy, the local East-North-Up frame, gravity with J2/J4, great-circle distance.

Positions are ECEF metres; angles radians; heading is clockwise from north.
"""
import math

import numpy as np

from .constants import J2, J4, MU_EARTH, WGS84_A, WGS84_E2

GRAVITY_MODELS = ("point", "j2", "j2j4")


def geodetic_to_ecef(h, lat, lon):
    """Geodetic altitude [m], latitude and longitude [rad] -> ECEF position [m]."""
    s, c = math.sin(lat), math.cos(lat)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * s * s)
    return np.array([(n + h) * c * math.cos(lon), (n + h) * c * math.sin(lon), (n * (1.0 - WGS84_E2) + h) * s])


def ecef_to_geodetic(r):
    """ECEF position [m] -> (geodetic altitude [m], latitude [rad], longitude [rad]).

    Fixed-point iteration on the latitude (converges to < 1e-9 m in a few iterations for
    |lat| < 89 deg, which covers every trajectory we run)."""
    x, y, z = float(r[0]), float(r[1]), float(r[2])
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    h = 0.0
    for _ in range(10):
        s = math.sin(lat)
        n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * s * s)
        h = p / math.cos(lat) - n
        lat_new = math.atan2(z, p * (1.0 - WGS84_E2 * n / (n + h)))
        if abs(lat_new - lat) < 1e-15:
            lat = lat_new
            break
        lat = lat_new
    return h, lat, lon


def enu_basis(lat, lon):
    """Unit vectors (east, north, up) of the local frame at geodetic (lat, lon), in ECEF."""
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    e = np.array([-so, co, 0.0])
    n = np.array([-sl * co, -sl * so, cl])
    u = np.array([cl * co, cl * so, sl])
    return e, n, u


def velocity_from_flight_angles(V, gamma, heading, lat, lon):
    """Speed [m/s], flight-path angle [rad, + up] and heading [rad, clockwise from north] -> ECEF velocity."""
    e, n, u = enu_basis(lat, lon)
    return V * (math.cos(gamma) * math.sin(heading) * e + math.cos(gamma) * math.cos(heading) * n + math.sin(gamma) * u)


def flight_angles(v, lat, lon):
    """ECEF velocity -> (speed, flight-path angle [rad], heading [rad] in [0, 2 pi))."""
    e, n, u = enu_basis(lat, lon)
    V = float(np.linalg.norm(v))
    if V == 0.0:
        return 0.0, 0.0, 0.0
    ve, vn, vu = float(v @ e), float(v @ n), float(v @ u)
    gamma = math.asin(max(-1.0, min(1.0, vu / V)))
    heading = math.atan2(ve, vn) % (2.0 * math.pi)
    return V, gamma, heading


def gravity(r, model="j2"):
    """Gravitational acceleration [m/s^2] at ECEF position r for a point mass, J2, or J2+J4 field."""
    if model not in GRAVITY_MODELS:
        raise ValueError("gravity model must be one of {}, got {!r}".format(GRAVITY_MODELS, model))
    x, y, z = float(r[0]), float(r[1]), float(r[2])
    r2 = x * x + y * y + z * z
    rr = math.sqrt(r2)
    base = -MU_EARTH / (rr * r2)
    if model == "point":
        return base * np.array([x, y, z])
    z2 = z * z / r2
    k2 = 1.5 * J2 * WGS84_A * WGS84_A / r2
    g = base * np.array([x * (1.0 - k2 * (5.0 * z2 - 1.0)),
                         y * (1.0 - k2 * (5.0 * z2 - 1.0)),
                         z * (1.0 - k2 * (5.0 * z2 - 3.0))])
    if model == "j2j4":
        c4 = 15.0 * MU_EARTH * J4 * WGS84_A ** 4 / (8.0 * rr ** 7)
        g = g + c4 * np.array([x * (1.0 - 14.0 * z2 + 21.0 * z2 * z2),
                               y * (1.0 - 14.0 * z2 + 21.0 * z2 * z2),
                               z * (5.0 - 70.0 * z2 / 3.0 + 21.0 * z2 * z2)])
    return g


def great_circle_distance(lat1, lon1, lat2, lon2, radius=WGS84_A):
    """Haversine distance [m] between two (lat, lon) [rad] points on a sphere of `radius`."""
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2.0 * radius * math.asin(min(1.0, math.sqrt(a)))
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_earth.py -q`
Expected: 11 passed (the parametrized round trip counts as 4)

- [ ] **Step 5: Commit**

```bash
git add reentry_model/earth.py tests/test_reentry_model_earth.py
git commit -m "reentry_model: WGS-84 geodesy, ENU frame, J2/J4 gravity, great-circle distance

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: Space-weather file parsing (`fap.py`)

**Files:**
- Create: `reentry_model/fap.py`, `tests/test_reentry_model_fap.py`

**Interfaces:**
- Consumes: `data/fap_day.dat` (DRAMA format: `dd/mm/yyyy F10 F3M SSN Ap kp...`, `#` comments).
- Produces: `load_fap_day(path) -> dict[datetime.date, FapRecord(f10: float, f3m: float, ssn: int, ap: float)]`; `solar_indices(records, day: datetime.date) -> SolarIndices(f107: float, f107a: float, ap: float)` with f107 = previous day's F10, f107a = that day's F3M (81-day mean), ap = that day's Ap; `DEFAULT_FAP_DAY` (path to `data/fap_day.dat`).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_fap.py`:
```python
"""fap.py: DRAMA fap_day.dat -> NRLMSISE-00 solar/geomagnetic inputs."""
from datetime import date

import pytest

from reentry_model import fap


def test_parses_the_repo_file_and_the_reference_epoch():
    records = fap.load_fap_day(fap.DEFAULT_FAP_DAY)
    assert records[date(2024, 8, 1)].f10 == 234.0
    assert records[date(2024, 8, 1)].f3m == 194.0
    assert records[date(2024, 8, 1)].ap == 19.0
    assert records[date(2024, 7, 31)].f10 == 234.0 and records[date(2024, 8, 2)].f10 == 246.0


def test_solar_indices_use_previous_day_f107_and_same_day_average_and_ap():
    records = fap.load_fap_day(fap.DEFAULT_FAP_DAY)
    s = fap.solar_indices(records, date(2024, 8, 2))
    assert s.f107 == 234.0          # F10.7 of 2024-08-01 (previous day)
    assert s.f107a == 194.0         # 81-day mean column of 2024-08-02
    assert s.ap == 8.0              # daily Ap of 2024-08-02


def test_missing_day_raises(tmp_path):
    p = tmp_path / "fap_day.dat"
    p.write_text("#d/mm/yyyy F10 F3M SSN Ap  3-hr Kp-Indices\n01/08/2024 234 194 259 019 5+4+3+2+3o2+2+2o\n")
    records = fap.load_fap_day(str(p))
    assert list(records) == [date(2024, 8, 1)]
    with pytest.raises(KeyError):
        fap.solar_indices(records, date(2024, 8, 1))      # needs 2024-07-31 for f107
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_fap.py -q`
Expected: ERROR `cannot import name 'fap'`

- [ ] **Step 3: Implement `fap.py`**

```python
"""DRAMA space-weather file (fap_day.dat) -> the inputs NRLMSISE-00 wants.

File format (one line per day): `dd/mm/yyyy F10 F3M SSN Ap kp-string`, `#` comment lines.
F10 = daily F10.7, F3M = its 81-day (3-month) mean, Ap = daily Ap.
"""
import os
from dataclasses import dataclass
from datetime import date, datetime, timedelta

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_FAP_DAY = os.path.join(REPO_ROOT, "data", "fap_day.dat")


@dataclass(frozen=True)
class FapRecord:
    f10: float
    f3m: float
    ssn: int
    ap: float


@dataclass(frozen=True)
class SolarIndices:
    f107: float      # F10.7 of the previous day (NRLMSISE-00 convention)
    f107a: float     # 81-day running mean
    ap: float        # daily Ap


def load_fap_day(path):
    """{date: FapRecord} for every data line of a fap_day.dat file."""
    records = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            day = datetime.strptime(parts[0], "%d/%m/%Y").date()
            records[day] = FapRecord(float(parts[1]), float(parts[2]), int(parts[3]), float(parts[4]))
    return records


def solar_indices(records, day):
    """NRLMSISE-00 inputs for `day`; KeyError if the day or the previous day is missing."""
    previous = records[day - timedelta(days=1)]
    today = records[day]
    return SolarIndices(f107=previous.f10, f107a=today.f3m, ap=today.ap)
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_fap.py -q`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add reentry_model/fap.py tests/test_reentry_model_fap.py
git commit -m "reentry_model: parse DRAMA fap_day.dat into NRLMSISE-00 solar indices

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Regime numbers and sphere drag (`aero.py`)

**Files:**
- Create: `reentry_model/aero.py`, `tests/test_reentry_model_aero.py`

**Interfaces:**
- Consumes: `constants`, `reentry_model/data/atdb_sphere.json`, the reference CSVs (test only).
- Produces:
  - `mean_free_path(rho, m_bar) -> float`, `knudsen(rho, m_bar, diameter) -> float`, `speed_of_sound(T, m_bar) -> float`, `mach(V, T, m_bar) -> float`
  - `SphereDragTables` with `.mach, .cd_fm, .cd_c` (np arrays), `SphereDragTables.from_json(path=None)`, `.cd_free_molecular(ma) -> float`, `.cd_continuum(ma) -> float`
  - bridging callables `SesamTable()`, `SesamErf(center=-0.845, width=0.585)`, `Sin2(kn_lo=0.01, kn_hi=1.0)`, `Matting()`; `bridging_by_name(name) -> callable` for names `sesam-table`, `sesam-erf`, `sin2`, `textbook`, `matting`; `BRIDGING_NAMES`
  - `drag_coefficient(kn, ma, tables, bridging) -> float`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_aero.py`:
```python
"""aero.py: Knudsen/Mach numbers and SESAM's sphere drag (tables + Knudsen bridging)."""
import csv
import math
import os

import numpy as np
import pytest

from reentry_model import aero
from reentry_model import constants as c

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(REPO_ROOT, "data", "reference_runs")
REFERENCE_CSVS = sorted(f for f in os.listdir(REFERENCE_DIR) if f.endswith(".csv"))


def test_knudsen_reproduces_sesam_at_the_break_off_altitude():
    # SESAM: rho = 2.727e-5 kg/m3 at 77.5 km (US76), D = 100 mm -> knudsen column 0.0298
    kn = aero.knudsen(2.727e-5, c.M_BAR_AIR, 0.100)
    assert kn == pytest.approx(0.0298, rel=0.02)
    assert aero.mean_free_path(2.727e-5, c.M_BAR_AIR) == pytest.approx(2.98e-3, rel=0.02)


def test_speed_of_sound_and_mach_at_sea_level():
    assert aero.speed_of_sound(288.15, c.M_BAR_AIR) == pytest.approx(340.3, abs=0.3)
    assert aero.mach(7500.0, 288.15, c.M_BAR_AIR) == pytest.approx(7500.0 / 340.3, rel=1e-3)


class TestTables:
    def test_loads_the_bundled_atdb(self):
        t = aero.SphereDragTables.from_json()
        assert list(t.mach) == [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
        assert t.cd_free_molecular(5.0) == 2.360635 and t.cd_continuum(30.0) == 0.913599

    def test_interpolates_and_clamps_above_mach_30(self):
        t = aero.SphereDragTables.from_json()
        assert t.cd_continuum(27.5) == pytest.approx(0.5 * (0.913411 + 0.913599))
        assert t.cd_free_molecular(21.07) == pytest.approx(2.062249 + (21.07 - 20.0) / 5.0 * (2.046842 - 2.062249))
        assert t.cd_continuum(35.0) == 0.913599 and t.cd_free_molecular(40.0) == 2.036933

    def test_below_mach_5_matches_sesam_clamping(self):
        t = aero.SphereDragTables.from_json()
        assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.0) == 0.898818
        assert t.cd_continuum(0.5) == pytest.approx(0.449409, abs=1e-6)     # exactly half the Ma-5 value
        assert t.cd_free_molecular(2.0) == 2.360635


class TestBridging:
    def test_sesam_table_limits_midpoint_and_monotonicity(self):
        f = aero.SesamTable()
        assert f(1e-3) == 0.0 and f(10.0) == 1.0 and f(41.0) == 1.0
        assert f(0.143) == pytest.approx(0.50, abs=0.03)
        kns = np.logspace(-3, 2, 200)
        values = [f(k) for k in kns]
        assert all(b >= a for a, b in zip(values, values[1:]))

    def test_sesam_erf_values(self):
        f = aero.SesamErf()
        assert f(0.01) == pytest.approx(0.0026, abs=5e-4)
        assert f(0.143) == pytest.approx(0.50, abs=0.01)
        assert f(1.0) == pytest.approx(0.98, abs=0.01)

    def test_sin2_textbook_form(self):
        f = aero.Sin2()
        assert f(0.001) == 0.0 and f(0.01) == 0.0 and f(1.0) == 1.0 and f(5.0) == 1.0
        assert f(0.1) == pytest.approx(0.5)

    def test_matting_is_not_available_yet(self):
        with pytest.raises(NotImplementedError):
            aero.Matting()(0.1)

    def test_factory(self):
        assert isinstance(aero.bridging_by_name("sesam-table"), aero.SesamTable)
        assert isinstance(aero.bridging_by_name("textbook"), aero.Sin2)
        with pytest.raises(ValueError):
            aero.bridging_by_name("legge")


class TestDragCoefficient:
    def test_limits(self):
        t = aero.SphereDragTables.from_json()
        assert aero.drag_coefficient(1e-4, 0.5, t, aero.SesamTable()) == pytest.approx(0.449409, abs=1e-6)
        assert aero.drag_coefficient(41.0, 21.07, t, aero.SesamTable()) == pytest.approx(2.059, abs=2e-3)   # R50 first row
        assert aero.drag_coefficient(0.0596, 26.23, t, aero.SesamTable()) == pytest.approx(1.125, abs=0.02) # 50 mm sweep, t = 0
        assert aero.drag_coefficient(0.0, 26.0, t, aero.SesamTable()) == t.cd_continuum(26.0)

    @pytest.mark.parametrize("name", REFERENCE_CSVS)
    def test_reproduces_the_reference_drag_column(self, name):
        t = aero.SphereDragTables.from_json()
        f = aero.SesamTable()
        rows = list(csv.DictReader(open(os.path.join(REFERENCE_DIR, name))))
        hyper, low = [], []
        for r in rows:
            kn, ma, cd = float(r["knudsen"]), float(r["mach"]), float(r["drag"])
            model = aero.drag_coefficient(kn, ma, t, f)
            if ma >= 5.0:
                hyper.append(model - cd)
            elif abs(ma - 1.0) > 0.02:                 # SESAM's subsonic switch sits at Ma = 1 (3-decimal columns)
                low.append(model - cd)
        hyper = np.array(hyper)
        assert math.sqrt(np.mean(hyper ** 2)) <= 0.010, name
        assert np.abs(hyper).max() <= 0.030, name
        assert np.abs(low).max() <= 1.5e-3, name       # the reference prints C_D with 3 decimals
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_aero.py -q`
Expected: ERROR `cannot import name 'aero'`

- [ ] **Step 3: Implement `aero.py`**

```python
"""Regime numbers and the sphere drag coefficient the way SESAM computes it.

Kn = lambda / D with a hard-sphere mean free path; Ma with gamma = 1.4; C_D = C_D,c(Ma) +
(C_D,fm(Ma) - C_D,c(Ma)) f(Kn) with the ATDB_SPHERE tables and a Knudsen bridging f.
See docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md sections 6.4-6.5.
"""
import json
import math
import os

import numpy as np

from . import DATA_DIR
from .constants import GAMMA_AIR, HARD_SPHERE_DIAMETER, K_BOLTZMANN

DEFAULT_ATDB = os.path.join(DATA_DIR, "atdb_sphere.json")


def mean_free_path(rho, m_bar):
    """Hard-sphere mean free path [m] for mass density rho [kg/m3] and mean molecular mass m_bar [kg]."""
    n = rho / m_bar
    return 1.0 / (math.sqrt(2.0) * math.pi * HARD_SPHERE_DIAMETER ** 2 * n)


def knudsen(rho, m_bar, diameter):
    return mean_free_path(rho, m_bar) / diameter


def speed_of_sound(T, m_bar):
    return math.sqrt(GAMMA_AIR * K_BOLTZMANN * T / m_bar)


def mach(V, T, m_bar):
    return V / speed_of_sound(T, m_bar)


class SphereDragTables:
    """C_D,fm(Ma) and C_D,c(Ma) of DRAMA's sphere aerothermal database (Ma 5..30, linear, clamped).

    Below Ma 5 SESAM holds the continuum value at the Ma-5 entry down to Ma 1 and uses exactly
    half of it below Ma 1 (measured on 466 rows, facts note s.2); the free-molecular value is
    simply clamped (Kn is negligible wherever Ma < 5)."""

    def __init__(self, mach_points, cd_free_molecular, cd_continuum):
        self.mach = np.asarray(mach_points, dtype=float)
        self.cd_fm = np.asarray(cd_free_molecular, dtype=float)
        self.cd_c = np.asarray(cd_continuum, dtype=float)

    @classmethod
    def from_json(cls, path=None):
        with open(path or DEFAULT_ATDB) as fh:
            d = json.load(fh)
        return cls(d["mach"], d["cd_free_molecular"], d["cd_continuum"])

    def cd_free_molecular(self, ma):
        return float(np.interp(ma, self.mach, self.cd_fm))          # np.interp clamps at both ends

    def cd_continuum(self, ma):
        if ma < 1.0:
            return 0.5 * float(self.cd_c[0])
        if ma < self.mach[0]:
            return float(self.cd_c[0])
        return float(np.interp(ma, self.mach, self.cd_c))


class SesamTable:
    """f(Kn) measured from SESAM's own drag output: bin means of (C_D - C_D,c)/(C_D,fm - C_D,c) over
    1544 rows (Ma >= 5, Kn > 1e-3) of the four reference runs and the 7.5 km/s, 300 K sweep runs of
    2026-09-17, in 0.25-decade bins of log10 Kn; linear in log10 Kn between bin centres, 0 below,
    1 above. Reproduces SESAM's C_D to rms < 0.01."""
    LOG10_KN = np.array([-2.125, -1.875, -1.625, -1.375, -1.125, -0.875, -0.625, -0.375, -0.125, 0.125])
    F = np.array([0.0, 0.0021, 0.0287, 0.1021, 0.2517, 0.4567, 0.6850, 0.8897, 0.9760, 1.0])

    def __call__(self, kn):
        if kn <= 0.0:
            return 0.0
        return float(np.interp(math.log10(kn), self.LOG10_KN, self.F, left=0.0, right=1.0))


class SesamErf:
    """Analytic fit of the same data: f = 1/2 [1 + erf((log10 Kn - center)/width)] (rms 0.0065, max error 0.027 near Kn 0.75)."""

    def __init__(self, center=-0.845, width=0.585):
        self.center, self.width = center, width

    def __call__(self, kn):
        if kn <= 0.0:
            return 0.0
        return 0.5 * (1.0 + math.erf((math.log10(kn) - self.center) / self.width))


class Sin2:
    """sin^2 ramp over [kn_lo, kn_hi] in log10 Kn; the textbook form is sin^2[pi (0.5 + 0.25 log10 Kn)] = Sin2(0.01, 1)."""

    def __init__(self, kn_lo=0.01, kn_hi=1.0):
        self.lo, self.hi = math.log10(kn_lo), math.log10(kn_hi)

    def __call__(self, kn):
        if kn <= 0.0:
            return 0.0
        x = (math.log10(kn) - self.lo) / (self.hi - self.lo)
        x = min(1.0, max(0.0, x))
        return math.sin(0.5 * math.pi * x) ** 2


class Matting:
    """Matting (1971) bridging relation: interface reserved; transcribe from J. Spacecraft Rockets 8(1) 35-40."""

    def __call__(self, kn):
        raise NotImplementedError("Matting (1971) bridging is not transcribed yet; use sesam-table, sesam-erf or sin2")


BRIDGING_NAMES = ("sesam-table", "sesam-erf", "sin2", "textbook", "matting")


def bridging_by_name(name):
    if name == "sesam-table":
        return SesamTable()
    if name == "sesam-erf":
        return SesamErf()
    if name in ("sin2", "textbook"):
        return Sin2()
    if name == "matting":
        return Matting()
    raise ValueError("bridging must be one of {}, got {!r}".format(BRIDGING_NAMES, name))


def drag_coefficient(kn, ma, tables, bridging):
    """SESAM's sphere C_D: continuum and free-molecular table values blended by f(Kn)."""
    cd_c = tables.cd_continuum(ma)
    cd_fm = tables.cd_free_molecular(ma)
    f = bridging(kn) if kn > 0.0 else 0.0
    return cd_c + (cd_fm - cd_c) * f
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_aero.py -q`
Expected: 16 passed (4 parametrized reference checks included). If `test_reproduces_the_reference_drag_column` fails on the rms/max bound, print the worst rows (`kn, ma, cd, model`) before changing anything: the table was fitted on exactly these files, so a failure means a transcription error in `SesamTable.F` or `cd_continuum`.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/aero.py tests/test_reentry_model_aero.py
git commit -m "reentry_model: Knudsen/Mach numbers and SESAM sphere drag with measured Knudsen bridging

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Reading SESAM reference runs (`sesam_io.py`)

**Files:**
- Create: `reentry_model/sesam_io.py`, `tests/test_reentry_model_sesam_io.py`

**Interfaces:**
- Consumes: the wrapper's run CSV (columns `time_s, altitude_km, velocity_kms, temperature_K, mass_kg, thick_mm, lat_deg, lon_deg, downrange_km, flight_path_deg, heading_deg, drag, lift, side, knudsen, mach, density_kgm3, dynamic_pressure_Pa, load_factor_g, ...`) and JSON (`inputs.diameter_mm, initial_velocity_kms, initial_altitude_km, flight_path_angle_deg, heading_deg, latitude_deg, longitude_deg, epoch_utc, material_density_kgm3, atmosphere, use_wind`, `run_name`).
- Produces: `Reference` dataclass with SI numpy arrays `time, altitude, velocity, lat, lon, flight_path, heading, downrange, drag, knudsen, mach, density, dynamic_pressure, temperature` plus `name, csv_path, json_path, sha256, diameter, material_density, atmosphere, use_wind, initial` where `initial` is an `InitialStateSpec(velocity, altitude, flight_path, heading, lat, lon, epoch)` in SI/radians; `load_reference(csv_path, json_path=None) -> Reference`; `REFERENCE_DIR`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_sesam_io.py`:
```python
"""sesam_io.py: a wrapper run (CSV + JSON) as a Reference."""
import math
import os
from datetime import datetime

import numpy as np
import pytest

from reentry_model import sesam_io

R50 = "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind"
R100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis"


def test_loads_r50_with_si_units_and_initial_state():
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, R50 + ".csv"))
    assert ref.name == R50 and ref.use_wind is False and ref.atmosphere == "nrlmsise"
    assert len(ref.time) > 500 and ref.time[0] == 0.0 and np.all(np.diff(ref.time) > 0)
    assert ref.altitude[0] == pytest.approx(115000.0, abs=1.0) and ref.altitude[-1] == 0.0
    assert ref.velocity[0] == pytest.approx(7500.0) and ref.diameter == 0.05
    assert ref.material_density == 2813.0
    assert ref.initial.velocity == 7500.0 and ref.initial.altitude == 115000.0
    assert ref.initial.flight_path == pytest.approx(math.radians(-0.959331))
    assert ref.initial.heading == pytest.approx(math.radians(347.168296))
    assert ref.initial.lat == pytest.approx(math.radians(29.546067))
    assert ref.initial.epoch == datetime(2024, 8, 1, 12, 53, 7)
    assert ref.knudsen[0] == pytest.approx(40.95647) and ref.drag[0] == pytest.approx(2.059)
    assert len(ref.sha256) == 64


def test_json_path_defaults_to_the_csv_stem_and_must_exist(tmp_path):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, R100 + ".csv"))
    assert ref.json_path.endswith(R100 + ".json") and ref.use_wind is True
    orphan = tmp_path / "x.csv"
    orphan.write_text("time_s,altitude_km\n0,1\n")
    with pytest.raises(FileNotFoundError):
        sesam_io.load_reference(str(orphan))
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_sesam_io.py -q`
Expected: ERROR `cannot import name 'sesam_io'`

- [ ] **Step 3: Implement `sesam_io.py`**

```python
"""A SESAM run written by sphere_reentry.py (history CSV + run JSON) as a Reference in SI units."""
import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(REPO_ROOT, "data", "reference_runs")


@dataclass(frozen=True)
class InitialStateSpec:
    velocity: float      # m/s, relative to the rotating atmosphere
    altitude: float      # m, geodetic
    flight_path: float   # rad, negative = descending
    heading: float       # rad, clockwise from north
    lat: float           # rad, geodetic
    lon: float           # rad
    epoch: datetime


@dataclass
class Reference:
    name: str
    csv_path: str
    json_path: str
    sha256: str
    diameter: float
    material_density: float
    atmosphere: str
    use_wind: bool
    initial: InitialStateSpec
    time: np.ndarray
    altitude: np.ndarray
    velocity: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    flight_path: np.ndarray
    heading: np.ndarray
    downrange: np.ndarray
    drag: np.ndarray
    knudsen: np.ndarray
    mach: np.ndarray
    density: np.ndarray
    dynamic_pressure: np.ndarray
    temperature: np.ndarray


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_reference(csv_path, json_path=None):
    """Read <run>.csv and <run>.json (default: same stem) into a Reference."""
    if json_path is None:
        json_path = os.path.splitext(csv_path)[0] + ".json"
    if not os.path.isfile(json_path):
        raise FileNotFoundError("reference JSON not found next to {}: {}".format(csv_path, json_path))
    with open(json_path) as fh:
        doc = json.load(fh)
    inp = doc["inputs"]
    initial = InitialStateSpec(
        velocity=float(inp["initial_velocity_kms"]) * 1e3,
        altitude=float(inp["initial_altitude_km"]) * 1e3,
        flight_path=math.radians(float(inp["flight_path_angle_deg"])),
        heading=math.radians(float(inp["heading_deg"])),
        lat=math.radians(float(inp["latitude_deg"])),
        lon=math.radians(float(inp["longitude_deg"])),
        epoch=datetime.strptime(inp["epoch_utc"], "%Y-%m-%dT%H:%M:%S"),
    )
    with open(csv_path) as fh:
        rows = list(csv.DictReader(fh))
    col = lambda key, scale=1.0: np.array([float(r[key]) for r in rows]) * scale
    return Reference(
        name=doc["run_name"], csv_path=os.path.abspath(csv_path), json_path=os.path.abspath(json_path),
        sha256=_sha256(csv_path),
        diameter=float(inp["diameter_mm"]) * 1e-3, material_density=float(inp["material_density_kgm3"]),
        atmosphere=inp.get("atmosphere", "static"), use_wind=bool(inp.get("use_wind", False)),
        initial=initial,
        time=col("time_s"), altitude=col("altitude_km", 1e3), velocity=col("velocity_kms", 1e3),
        lat=np.radians(col("lat_deg")), lon=np.radians(col("lon_deg")),
        flight_path=np.radians(col("flight_path_deg")), heading=np.radians(col("heading_deg")),
        downrange=col("downrange_km", 1e3), drag=col("drag"), knudsen=col("knudsen"), mach=col("mach"),
        density=col("density_kgm3"), dynamic_pressure=col("dynamic_pressure_Pa"), temperature=col("temperature_K"),
    )
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_sesam_io.py -q`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add reentry_model/sesam_io.py tests/test_reentry_model_sesam_io.py
git commit -m "reentry_model: read SESAM reference runs into SI arrays

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: Atmosphere models (`atmosphere.py`)

**Files:**
- Create: `reentry_model/atmosphere.py`, `tests/test_reentry_model_atmosphere.py`

**Interfaces:**
- Consumes: `constants`, `fap.SolarIndices`, `sesam_io.Reference`, `aero.mean_free_path`/`HARD_SPHERE_DIAMETER`, `pymsis`.
- Produces:
  - `Freestream(rho, T, p, m_bar, wind_enu)` frozen dataclass (SI; `wind_enu` np array [east, north, up] m/s)
  - wind models `NoWind()`, `StaticProfileWind(path=None)` with `.wind_enu(h) -> np.ndarray`
  - atmospheres with `.state(t, h, lat, lon) -> Freestream` (t s since epoch, h m, lat/lon rad): `US76TableAtmosphere(path=None, wind=None)`, `NRLMSISE00Atmosphere(epoch, solar, wind=None)`, `ReplayAtmosphere(reference, wind=None)`, `VacuumAtmosphere()` (rho = p = 0, T = 200 K; for dynamics tests)
  - `DEFAULT_US76` path constant

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_atmosphere.py`:
```python
"""atmosphere.py: US76 table, NRLMSISE-00 (pymsis), SESAM replay, winds."""
import math
import os
from datetime import date

import numpy as np
import pytest

from reentry_model import aero, atmosphere, fap, sesam_io
from reentry_model import constants as c

R100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind"
R50 = "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind"
LAT, LON = math.radians(29.546067), math.radians(-82.134333)


def ref(name):
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, name + ".csv"))


class TestUS76:
    def test_sea_level_and_70km(self):
        atm = atmosphere.US76TableAtmosphere()
        s0 = atm.state(0.0, 0.0, LAT, LON)
        assert s0.rho == pytest.approx(1.225, abs=1e-3) and s0.T == pytest.approx(288.15, abs=0.01)
        assert s0.p == pytest.approx(101325.0, abs=1.0) and s0.m_bar == c.M_BAR_AIR
        assert np.all(s0.wind_enu == 0.0)
        s70 = atm.state(0.0, 70000.0, LAT, LON)
        assert s70.rho == pytest.approx(8.283e-5, rel=2e-3) and s70.T == pytest.approx(219.6, abs=0.1)

    def test_matches_sesam_static_density_at_break_off(self):
        atm = atmosphere.US76TableAtmosphere()
        assert atm.state(0.0, 77500.133, LAT, LON).rho == pytest.approx(2.727e-5, rel=0.01)

    def test_log_linear_between_rows_and_top_of_table(self):
        atm = atmosphere.US76TableAtmosphere()
        a, m, b = (atm.state(0.0, h, LAT, LON).rho for h in (70000.0, 70050.0, 70100.0))
        assert m == pytest.approx(math.sqrt(a * b), rel=1e-6)
        with pytest.raises(ValueError):
            atm.state(0.0, 150001.0, LAT, LON)

    def test_static_wind_profile(self):
        wind = atmosphere.StaticProfileWind()
        w = wind.wind_enu(60000.0)
        assert w[0] == pytest.approx(10.90, abs=0.01) and w[1] == pytest.approx(-1.30, abs=0.01) and w[2] == 0.0
        atm = atmosphere.US76TableAtmosphere(wind=wind)
        assert np.allclose(atm.state(0.0, 60000.0, LAT, LON).wind_enu, w)


class TestNRLMSISE00:
    def test_reproduces_sesam_density_at_the_reference_starts(self):
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), date(2024, 8, 1))
        for name, tol in ((R100, 0.03), (R50, 0.05)):
            r = ref(name)
            atm = atmosphere.NRLMSISE00Atmosphere(r.initial.epoch, solar)
            s = atm.state(0.0, r.initial.altitude, r.initial.lat, r.initial.lon)
            assert s.rho == pytest.approx(r.density[0], rel=tol), (name, s.rho, r.density[0])
            assert s.T > 150.0 and s.p > 0.0

    def test_mean_molecular_mass_falls_off_above_100km(self):
        solar = fap.SolarIndices(f107=234.0, f107a=194.0, ap=19.0)
        atm = atmosphere.NRLMSISE00Atmosphere(ref(R50).initial.epoch, solar)
        low = atm.state(0.0, 77500.0, LAT, LON).m_bar / c.ATOMIC_MASS_UNIT
        high = atm.state(0.0, 115000.0, LAT, LON).m_bar / c.ATOMIC_MASS_UNIT
        assert low == pytest.approx(28.96, abs=0.2)
        assert 25.0 < high < 28.0                       # SESAM's Kn column implies ~26.6 u at 115 km


class TestReplay:
    @pytest.mark.parametrize("name", [R100, R50])
    def test_reproduces_density_knudsen_and_mach_columns(self, name):
        r = ref(name)
        atm = atmosphere.ReplayAtmosphere(r)
        for i in range(0, len(r.time), 25):
            s = atm.state(r.time[i], r.altitude[i], r.lat[i], r.lon[i])
            assert s.rho == pytest.approx(r.density[i], rel=1e-6)
            assert aero.knudsen(s.rho, s.m_bar, r.diameter) == pytest.approx(r.knudsen[i], rel=0.02)
            if r.mach[i] > 0.3:
                assert aero.mach(r.velocity[i], s.T, s.m_bar) == pytest.approx(r.mach[i], rel=0.01)


def test_vacuum():
    s = atmosphere.VacuumAtmosphere().state(0.0, 1e5, LAT, LON)
    assert s.rho == 0.0 and s.p == 0.0 and s.T == 200.0 and np.all(s.wind_enu == 0.0)
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_atmosphere.py -q`
Expected: ERROR `cannot import name 'atmosphere'`

- [ ] **Step 3: Implement `atmosphere.py`**

```python
"""Freestream state along the trajectory: NRLMSISE-00 (pymsis), the US76 table SESAM ships, or a
replay of a SESAM run's own density/temperature (diagnostic). Winds: none or the static profile."""
import csv
import math
import os
from dataclasses import dataclass
from datetime import timedelta

import numpy as np

from . import DATA_DIR
from .constants import ATOMIC_MASS_UNIT, GAMMA_AIR, HARD_SPHERE_DIAMETER, K_BOLTZMANN, M_BAR_AIR

DEFAULT_US76 = os.path.join(DATA_DIR, "us76_static_environment.csv")


@dataclass(frozen=True)
class Freestream:
    rho: float           # kg/m3
    T: float             # K
    p: float             # Pa
    m_bar: float         # kg, mean molecular mass
    wind_enu: np.ndarray  # m/s, (east, north, up)


def _read_us76(path):
    rows = []
    with open(path) as fh:
        for r in csv.reader(fh):
            if r and not r[0].startswith("#"):
                rows.append([float(v) for v in r])
    a = np.array(rows)
    return {"h": a[:, 0], "p": a[:, 1], "rho": a[:, 2], "T": a[:, 3], "wn": a[:, 4], "we": a[:, 5], "wd": a[:, 6]}


class NoWind:
    def wind_enu(self, h):
        return np.zeros(3)


class StaticProfileWind:
    """The wind columns of DRAMA's static environment table (north, east, down -> east, north, up)."""

    def __init__(self, path=None):
        t = _read_us76(path or DEFAULT_US76)
        self._h, self._wn, self._we, self._wd = t["h"], t["wn"], t["we"], t["wd"]

    def wind_enu(self, h):
        return np.array([np.interp(h, self._h, self._we), np.interp(h, self._h, self._wn), -np.interp(h, self._h, self._wd)])


class US76TableAtmosphere:
    """US Standard Atmosphere 1976 as tabulated by DRAMA (0-150 km, 100 m); log-linear rho and p, linear T."""

    def __init__(self, path=None, wind=None):
        t = _read_us76(path or DEFAULT_US76)
        self._h, self._T = t["h"], t["T"]
        self._log_rho, self._log_p = np.log(t["rho"]), np.log(t["p"])
        self.wind = wind or NoWind()

    def state(self, t, h, lat, lon):
        if h < self._h[0] or h > self._h[-1]:
            raise ValueError("altitude {:.1f} m outside the US76 table (0..{:.0f} m)".format(h, self._h[-1]))
        rho = math.exp(np.interp(h, self._h, self._log_rho))
        p = math.exp(np.interp(h, self._h, self._log_p))
        T = float(np.interp(h, self._h, self._T))
        return Freestream(rho, T, p, M_BAR_AIR, self.wind.wind_enu(h))


class NRLMSISE00Atmosphere:
    """NRLMSISE-00 through pymsis (version 0) at the run epoch + t, with the fap-file solar indices."""

    def __init__(self, epoch, solar, wind=None):
        from pymsis import msis        # imported here so the rest of the package works without pymsis
        self._msis = msis
        self.epoch = epoch
        self.solar = solar
        self.wind = wind or NoWind()

    def state(self, t, h, lat, lon):
        when = np.datetime64(self.epoch + timedelta(seconds=float(t)))
        out = self._msis.calculate(np.array([when]), [math.degrees(lon)], [math.degrees(lat)], [h / 1e3],
                                   f107s=[self.solar.f107], f107as=[self.solar.f107a],
                                   aps=[[self.solar.ap] * 7], version=0)
        row = np.asarray(out).reshape(-1, 11)[0]
        rho = float(row[0])
        n = float(np.nansum(row[1:10]))                 # species number densities (NO is NaN for version 0)
        T = float(row[10])
        m_bar = rho / n
        return Freestream(rho, T, n * K_BOLTZMANN * T, m_bar, self.wind.wind_enu(h))


class ReplayAtmosphere:
    """rho(h), m_bar(h), T(h) reconstructed from a SESAM run: rho from its density column; m_bar = air below
    100 km and, above, back-solved from the Knudsen column (lambda = Kn D, n = 1/(sqrt2 pi d^2 lambda));
    T from Ma and V (T = (V/Ma)^2 m_bar / (gamma k_B)) on rows with Ma > 0.3. Interpolated in altitude."""

    def __init__(self, reference, wind=None):
        order = np.argsort(reference.altitude)
        h = reference.altitude[order]
        rho = reference.density[order]
        kn = reference.knudsen[order]
        ma = reference.mach[order]
        V = reference.velocity[order]
        keep = np.concatenate([[True], np.diff(h) > 0])          # strictly increasing altitude nodes
        h, rho, kn, ma, V = h[keep], rho[keep], kn[keep], ma[keep], V[keep]
        m_bar = np.full_like(h, M_BAR_AIR)
        high = (h >= 100e3) & (kn > 0)
        lam = kn[high] * reference.diameter
        n = 1.0 / (math.sqrt(2.0) * math.pi * HARD_SPHERE_DIAMETER ** 2 * lam)
        m_bar[high] = np.clip(rho[high] / n, 16.0 * ATOMIC_MASS_UNIT, 30.0 * ATOMIC_MASS_UNIT)
        ok = ma > 0.3
        T = (V[ok] / ma[ok]) ** 2 * m_bar[ok] / (GAMMA_AIR * K_BOLTZMANN)
        self._h, self._log_rho, self._m_bar = h, np.log(rho), m_bar
        self._hT, self._T = h[ok], T
        self.wind = wind or NoWind()

    def state(self, t, h, lat, lon):
        rho = math.exp(np.interp(h, self._h, self._log_rho))
        m_bar = float(np.interp(h, self._h, self._m_bar))
        T = float(np.interp(h, self._hT, self._T))
        return Freestream(rho, T, rho / m_bar * K_BOLTZMANN * T, m_bar, self.wind.wind_enu(h))


class VacuumAtmosphere:
    """No atmosphere at all (dynamics tests)."""

    def state(self, t, h, lat, lon):
        return Freestream(0.0, 200.0, 0.0, M_BAR_AIR, np.zeros(3))
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_atmosphere.py -q`
Expected: 9 passed. `TestNRLMSISE00::test_reproduces_sesam_density_at_the_reference_starts` is the one informative tolerance: if it fails, print `s.rho / r.density[0]` for both cases and record the ratios in the README verification table (Task 10) — a ratio outside 0.95–1.05 means SESAM's F10.7/Ap convention differs from `fap.solar_indices` (try same-day F10.7 instead of previous-day and report which matches).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/atmosphere.py tests/test_reentry_model_atmosphere.py
git commit -m "reentry_model: US76 table, NRLMSISE-00 via pymsis, SESAM replay atmosphere, static winds

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: The body hook and the trajectory integrator (`body.py`, `trajectory.py`)

**Files:**
- Create: `reentry_model/body.py`, `reentry_model/trajectory.py`, `tests/test_reentry_model_trajectory.py`

**Interfaces:**
- Consumes: `earth`, `aero`, `atmosphere` (any object with `.state(t, h, lat, lon) -> Freestream`), `constants`.
- Produces:
  - `body.Body` protocol (`mass(t)`, `temperature(t)`, `on_step(t, state, freestream, aero)`), `body.ConstantBody(mass_kg, temperature_K=300.0)`
  - `trajectory.InitialState(velocity, altitude, flight_path, heading, lat, lon, epoch)` (SI, rad)
  - `trajectory.Settings(diameter, gravity="j2", rotating_frame=True, rtol=1e-9, atol_position=1e-6, atol_velocity=1e-9, cadence=1.0, t_max=3600.0, ground_altitude=0.0, escape_altitude=150e3)`
  - `trajectory.AeroState(h, lat, lon, freestream, v_rel, V, kn, ma, cd, a_drag, q_dyn)`
  - `trajectory.History(columns: dict[str, np.ndarray], states: np.ndarray (N x 6), end_reason: str, results: dict)`; `trajectory.CSV_COLUMNS`
  - `trajectory.Simulator(initial, body, atmosphere, tables, bridging, settings)` with `.initial_state_vector()`, `.aero_state(t, r, v)`, `.rhs(t, y)`, `.run(extra_times=None) -> History`
  - `trajectory.write_history_csv(history, path)`, `trajectory.read_history_csv(path) -> History`, `trajectory.write_run_json(path, doc)`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_trajectory.py`:
```python
"""body.py and trajectory.py: the rotating-Earth 3-DOF integrator and its outputs."""
import csv
import json
import math
from datetime import datetime

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, earth
from reentry_model import trajectory as tj
from reentry_model.constants import MU_EARTH, OMEGA_EARTH

EPOCH = datetime(2024, 8, 1, 12, 53, 7)
R100 = tj.InitialState(7500.0, 77500.133, math.radians(-0.959331), math.radians(347.168296),
                       math.radians(29.546067), math.radians(-82.134333), EPOCH)
ORBIT = tj.InitialState(7800.0, 300e3, 0.0, math.radians(90.0), 0.0, 0.0, EPOCH)
MASS_100MM = 2813.0 * 4.0 / 3.0 * math.pi * 0.05 ** 3


def make(initial, atm, mass=MASS_100MM, **settings):
    s = tj.Settings(diameter=0.1, **settings)
    return tj.Simulator(initial, body.ConstantBody(mass, 300.0), atm, aero.SphereDragTables.from_json(), aero.SesamTable(), s)


def test_constant_body():
    b = body.ConstantBody(1.5, 310.0)
    assert b.mass(0.0) == 1.5 and b.mass(100.0) == 1.5 and b.temperature(5.0) == 310.0
    assert b.on_step(0.0, None, None, None) is None


def test_initial_state_vector_round_trips_the_inputs():
    y0 = make(R100, atmosphere.VacuumAtmosphere()).initial_state_vector()
    h, lat, lon = earth.ecef_to_geodetic(y0[:3])
    assert h == pytest.approx(77500.133, abs=1e-6) and lat == pytest.approx(R100.lat) and lon == pytest.approx(R100.lon)
    V, gamma, heading = earth.flight_angles(y0[3:], lat, lon)
    assert V == pytest.approx(7500.0) and gamma == pytest.approx(R100.flight_path) and heading == pytest.approx(R100.heading)


def test_vacuum_point_mass_non_rotating_conserves_energy_and_angular_momentum():
    sim = make(ORBIT, atmosphere.VacuumAtmosphere(), gravity="point", rotating_frame=False,
               cadence=60.0, t_max=600.0, escape_altitude=1e9)
    hist = sim.run()
    assert hist.end_reason == "t_max" and hist.states.shape == (11, 6)
    r, v = hist.states[:, :3], hist.states[:, 3:]
    energy = 0.5 * np.sum(v * v, axis=1) - MU_EARTH / np.linalg.norm(r, axis=1)
    momentum = np.linalg.norm(np.cross(r, v), axis=1)
    assert np.max(np.abs(energy - energy[0]) / abs(energy[0])) < 1e-9
    assert np.max(np.abs(momentum - momentum[0]) / momentum[0]) < 1e-9


def test_vacuum_rotating_frame_conserves_the_jacobi_integral():
    sim = make(ORBIT, atmosphere.VacuumAtmosphere(), gravity="point", rotating_frame=True,
               cadence=60.0, t_max=600.0, escape_altitude=1e9)
    hist = sim.run()
    r, v = hist.states[:, :3], hist.states[:, 3:]
    jacobi = 0.5 * np.sum(v * v, axis=1) - MU_EARTH / np.linalg.norm(r, axis=1) \
        - 0.5 * OMEGA_EARTH ** 2 * (r[:, 0] ** 2 + r[:, 1] ** 2)
    assert np.max(np.abs(jacobi - jacobi[0]) / abs(jacobi[0])) < 1e-9


def test_flight_path_angle_rate_matches_sesam_and_vinh():
    # SESAM: -0.00790 deg/s at the R100 start; Vinh's rotating-Earth equations: -0.00792 deg/s (facts note s.7).
    sim = make(R100, atmosphere.VacuumAtmosphere(), cadence=2.5, t_max=2.5, escape_altitude=1e9)
    hist = sim.run()
    gamma = hist.columns["flight_path_deg"]
    rate = (gamma[-1] - gamma[0]) / hist.columns["time_s"][-1]
    assert rate == pytest.approx(-0.0079, abs=0.0002)


def test_us76_flight_ends_on_the_ground_and_samples_as_requested():
    sim = make(R100, atmosphere.US76TableAtmosphere(), cadence=10.0)
    hist = sim.run(extra_times=[3.3, 7.7])
    t = hist.columns["time_s"]
    assert hist.end_reason == "ground"
    assert abs(hist.columns["altitude_km"][-1]) < 1e-6                 # within 1 mm of h = 0
    assert 300.0 < t[-1] < 450.0                                        # SESAM: 366 s on US76
    assert t[0] == 0.0 and np.all(np.diff(t) > 0)
    for wanted in (3.3, 7.7, 10.0, 20.0):
        assert np.any(np.isclose(t, wanted))
    assert hist.results["impact_time_s"] == t[-1] and hist.results["end_reason"] == "ground"
    assert 60.0 < hist.results["knudsen_crossings"]["0.01"] < 76.0     # km; SESAM: 70.1 km
    assert np.all(hist.columns["mass_kg"] == MASS_100MM) and np.all(hist.columns["temperature_K"] == 300.0)
    assert hist.columns["drag"][0] == pytest.approx(aero.drag_coefficient(hist.columns["knudsen"][0], hist.columns["mach"][0],
                                                                            aero.SphereDragTables.from_json(), aero.SesamTable()))
    assert 6.0 < hist.results["max_deceleration_g"] < 10.0 and 35.0 < hist.results["altitude_of_max_deceleration_km"] < 50.0   # SESAM: 8.2 g at 41 km


def test_csv_and_json_round_trip(tmp_path):
    sim = make(R100, atmosphere.US76TableAtmosphere(), cadence=50.0, t_max=100.0)
    hist = sim.run()
    csv_path = tmp_path / "h.csv"
    tj.write_history_csv(hist, str(csv_path))
    with open(csv_path) as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0].keys()) == tj.CSV_COLUMNS and len(rows) == len(hist.columns["time_s"])
    back = tj.read_history_csv(str(csv_path))
    assert np.allclose(back.columns["velocity_kms"], hist.columns["velocity_kms"], rtol=1e-6)
    tj.write_run_json(str(tmp_path / "h.json"), {"results": hist.results, "epoch": EPOCH})
    doc = json.load(open(tmp_path / "h.json"))
    assert doc["results"]["end_reason"] == "t_max" and doc["epoch"].startswith("2024-08-01")
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_trajectory.py -q`
Expected: ERROR `cannot import name 'body'`

- [ ] **Step 3: Implement `body.py`**

```python
"""The body the trajectory carries. Step 1: constant mass and temperature. Step 2 replaces
ConstantBody with a thermal model whose on_step advances the temperature field."""
from dataclasses import dataclass
from typing import Protocol


class Body(Protocol):
    def mass(self, t: float) -> float: ...
    def temperature(self, t: float) -> float: ...
    def on_step(self, t, state, freestream, aero) -> None: ...


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def temperature(self, t):
        return self.temperature_K

    def on_step(self, t, state, freestream, aero):
        return None
```

- [ ] **Step 4: Implement `trajectory.py`**

```python
"""Rotating-Earth 3-DOF trajectory of a sphere, integrated in ECEF Cartesian coordinates.

d r/dt = v ;  d v/dt = g(r) - 2 w x v - w x (w x r) + a_drag,  a_drag = -1/2 rho |v_rel| v_rel C_D A / m
(spec section 6). Positions/velocities are relative to the rotating Earth; v_rel subtracts the wind.
"""
import csv
import json
import math
import time
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
from scipy.integrate import solve_ivp

from . import aero, earth
from .constants import G0, OMEGA_EARTH

CSV_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "temperature_K", "mass_kg", "thick_mm", "lat_deg",
               "lon_deg", "downrange_km", "flight_path_deg", "heading_deg", "drag", "lift", "side", "knudsen",
               "mach", "density_kgm3", "dynamic_pressure_Pa", "load_factor_g"]
KNUDSEN_THRESHOLDS = (10.0, 1.0, 0.1, 0.01)


@dataclass(frozen=True)
class InitialState:
    velocity: float      # m/s relative to the rotating atmosphere
    altitude: float      # m, geodetic
    flight_path: float   # rad, negative = descending
    heading: float       # rad, clockwise from north
    lat: float           # rad, geodetic
    lon: float           # rad
    epoch: datetime


@dataclass
class Settings:
    diameter: float                  # m
    gravity: str = "j2"              # point | j2 | j2j4
    rotating_frame: bool = True
    rtol: float = 1e-9
    atol_position: float = 1e-6      # m
    atol_velocity: float = 1e-9      # m/s
    cadence: float = 1.0             # s between history samples
    t_max: float = 3600.0            # s
    ground_altitude: float = 0.0     # m
    escape_altitude: float = 150e3   # m


@dataclass
class AeroState:
    h: float
    lat: float
    lon: float
    freestream: object
    v_rel: np.ndarray
    V: float
    kn: float
    ma: float
    cd: float
    a_drag: np.ndarray
    q_dyn: float


@dataclass
class History:
    columns: dict
    states: np.ndarray
    end_reason: str
    results: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.columns["time_s"])


class Simulator:
    def __init__(self, initial, body, atmosphere, tables, bridging, settings):
        self.initial, self.body, self.atmosphere = initial, body, atmosphere
        self.tables, self.bridging, self.settings = tables, bridging, settings
        self.area = math.pi * settings.diameter ** 2 / 4.0
        self.omega = np.array([0.0, 0.0, OMEGA_EARTH]) if settings.rotating_frame else np.zeros(3)
        self.r0 = earth.geodetic_to_ecef(initial.altitude, initial.lat, initial.lon)
        self.v0 = earth.velocity_from_flight_angles(initial.velocity, initial.flight_path, initial.heading,
                                                    initial.lat, initial.lon)

    def initial_state_vector(self):
        return np.concatenate([self.r0, self.v0])

    def aero_state(self, t, r, v):
        h, lat, lon = earth.ecef_to_geodetic(r)
        fs = self.atmosphere.state(t, h, lat, lon)
        e, n, u = earth.enu_basis(lat, lon)
        wind = fs.wind_enu[0] * e + fs.wind_enu[1] * n + fs.wind_enu[2] * u
        v_rel = v - wind
        V = float(np.linalg.norm(v_rel))
        ma = aero.mach(V, fs.T, fs.m_bar)
        if fs.rho > 0.0:
            kn = aero.knudsen(fs.rho, fs.m_bar, self.settings.diameter)
            cd = aero.drag_coefficient(kn, ma, self.tables, self.bridging)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * self.area / self.body.mass(t)
        else:
            kn, cd, a_drag = math.inf, self.tables.cd_free_molecular(ma), np.zeros(3)
        return AeroState(h, lat, lon, fs, v_rel, V, kn, ma, cd, a_drag, 0.5 * fs.rho * V * V)

    def rhs(self, t, y):
        r, v = y[:3], y[3:]
        a = earth.gravity(r, self.settings.gravity) - 2.0 * np.cross(self.omega, v) \
            - np.cross(self.omega, np.cross(self.omega, r)) + self.aero_state(t, r, v).a_drag
        return np.concatenate([v, a])

    def _events(self):
        s = self.settings

        def ground(t, y):
            return earth.ecef_to_geodetic(y[:3])[0] - s.ground_altitude

        def escape(t, y):
            return s.escape_altitude - earth.ecef_to_geodetic(y[:3])[0]

        ground.terminal, ground.direction = True, -1
        escape.terminal, escape.direction = True, -1
        return [ground, escape]

    def run(self, extra_times=None):
        s = self.settings
        started = time.perf_counter()
        sol = solve_ivp(self.rhs, (0.0, s.t_max), self.initial_state_vector(), method="DOP853",
                        rtol=s.rtol, atol=[s.atol_position] * 3 + [s.atol_velocity] * 3,
                        dense_output=True, events=self._events())
        if not sol.success:
            raise RuntimeError("integration failed: {}".format(sol.message))
        t_end = float(sol.t[-1])
        if sol.t_events[0].size:
            end_reason = "ground"
        elif sol.t_events[1].size:
            end_reason = "escape"
        else:
            end_reason = "t_max"
        times = np.arange(0.0, t_end, s.cadence)
        if extra_times is not None:
            extra = np.asarray(extra_times, dtype=float)
            times = np.union1d(times, extra[extra <= t_end])
        times = np.union1d(times, [t_end])
        cols = {k: [] for k in CSV_COLUMNS}
        states = []
        for t in times:
            y = sol.sol(t)
            r, v = y[:3], y[3:]
            a = self.aero_state(t, r, v)
            self.body.on_step(t, y, a.freestream, a)
            V, gamma, heading = earth.flight_angles(v, a.lat, a.lon)
            cols["time_s"].append(t)
            cols["altitude_km"].append(a.h / 1e3)
            cols["velocity_kms"].append(V / 1e3)
            cols["temperature_K"].append(self.body.temperature(t))
            cols["mass_kg"].append(self.body.mass(t))
            cols["thick_mm"].append(s.diameter * 500.0)
            cols["lat_deg"].append(math.degrees(a.lat))
            cols["lon_deg"].append(math.degrees(a.lon))
            cols["downrange_km"].append(earth.great_circle_distance(self.initial.lat, self.initial.lon, a.lat, a.lon) / 1e3)
            cols["flight_path_deg"].append(math.degrees(gamma))
            cols["heading_deg"].append(math.degrees(heading))
            cols["drag"].append(a.cd)
            cols["lift"].append(0.0)
            cols["side"].append(0.0)
            cols["knudsen"].append(a.kn)
            cols["mach"].append(a.ma)
            cols["density_kgm3"].append(a.freestream.rho)
            cols["dynamic_pressure_Pa"].append(a.q_dyn)
            cols["load_factor_g"].append(float(np.linalg.norm(a.a_drag)) / G0)     # = SESAM's column (drag / m g0, verified)
            states.append(y)
        columns = {k: np.array(v, dtype=float) for k, v in cols.items()}
        history = History(columns, np.array(states), end_reason)
        history.results = self._results(history, sol, time.perf_counter() - started)
        return history

    def _results(self, history, sol, runtime):
        c = history.columns
        i_dec = int(np.argmax(c["load_factor_g"]))
        crossings = {}
        for thr in KNUDSEN_THRESHOLDS:
            below = np.nonzero(c["knudsen"] < thr)[0]
            crossings["{:g}".format(thr)] = float(c["altitude_km"][below[0]]) if below.size else None
        results = {
            "end_reason": history.end_reason,
            "final_time_s": float(c["time_s"][-1]),
            "impact_time_s": float(c["time_s"][-1]) if history.end_reason == "ground" else None,
            "final_velocity_kms": float(c["velocity_kms"][-1]),
            "final_altitude_km": float(c["altitude_km"][-1]),
            "max_deceleration_g": float(c["load_factor_g"][i_dec]),
            "altitude_of_max_deceleration_km": float(c["altitude_km"][i_dec]),
            "time_of_max_deceleration_s": float(c["time_s"][i_dec]),
            "max_dynamic_pressure_Pa": float(c["dynamic_pressure_Pa"].max()),
            "knudsen_start": float(c["knudsen"][0]),
            "knudsen_crossings": crossings,
            "n_samples": len(history),
            "rhs_evaluations": int(sol.nfev),
            "runtime_s": runtime,
        }
        return results


def write_history_csv(history, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(history.columns[k][i]) for k in CSV_COLUMNS])


def read_history_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    columns = {k: np.array([float(r[k]) for r in rows]) for k in CSV_COLUMNS if k in rows[0]}
    return History(columns, np.zeros((len(rows), 6)), "unknown")


def write_run_json(path, doc):
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
```

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_trajectory.py -q`
Expected: 8 passed, in well under a minute (the US76 flight is ~370 s of simulated time). If `test_flight_path_angle_rate_matches_sesam_and_vinh` fails, check the sign of the Coriolis term and that `velocity_from_flight_angles` uses the geodetic latitude — those are the only two ways to be off by the 20 % the test discriminates.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/body.py reentry_model/trajectory.py tests/test_reentry_model_trajectory.py
git commit -m "reentry_model: rotating-Earth 3-DOF integrator in ECEF with events, sampling, CSV/JSON output

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 8: Comparison metrics and plots (`compare.py`)

**Files:**
- Create: `reentry_model/compare.py`, `tests/test_reentry_model_compare.py`

**Interfaces:**
- Consumes: `trajectory.History`, `sesam_io.Reference`, `analysis/plot_style.py` (optional styling).
- Produces: `align(history, reference) -> dict[str, np.ndarray]` (keys `t, V_model, V_ref, h_model, h_ref, gamma_model, gamma_ref, heading_model, heading_ref, lat_model, lat_ref, lon_model, lon_ref, kn_model, kn_ref, cd_model, cd_ref`; SI, angles deg); `metrics(history, reference) -> dict` with top-level `n_points, model_end_time_s, reference_end_time_s, d_end_time_s, d_end_time_rel, final_velocity_model_ms, final_velocity_reference_ms` and phase dicts `all` and `hypersonic` each with `dV_max_ms, dV_rms_ms, dV_rel_max, dV_rel_rms, dh_max_m, dh_rms_m, n_points`; `plot_all(history, reference, outdir, title) -> list[str]` writing `velocity_time.png, altitude_time.png, altitude_velocity.png, angles_time.png, ground_track.png, regime_drag.png`; `PLOT_NAMES`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_compare.py`:
```python
"""compare.py: model-vs-SESAM metrics and plots on a synthetic history with known offsets."""
import os

import numpy as np
import pytest

from reentry_model import compare, sesam_io
from reentry_model import trajectory as tj

R100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind"


def synthetic_history(ref, dv=5.0, dh=-20.0, truncate=None):
    n = len(ref.time) if truncate is None else truncate
    cols = {k: np.zeros(n) for k in tj.CSV_COLUMNS}
    cols["time_s"] = ref.time[:n].copy()
    cols["velocity_kms"] = (ref.velocity[:n] + dv) / 1e3
    cols["altitude_km"] = (ref.altitude[:n] + dh) / 1e3
    cols["flight_path_deg"] = np.degrees(ref.flight_path[:n])
    cols["heading_deg"] = np.degrees(ref.heading[:n])
    cols["lat_deg"] = np.degrees(ref.lat[:n])
    cols["lon_deg"] = np.degrees(ref.lon[:n])
    cols["knudsen"] = ref.knudsen[:n]
    cols["drag"] = ref.drag[:n]
    cols["mach"] = ref.mach[:n]
    cols["density_kgm3"] = ref.density[:n]
    return tj.History(cols, np.zeros((n, 6)), "ground", {})


@pytest.fixture(scope="module")
def ref():
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, R100 + ".csv"))


def test_metrics_recover_known_offsets(ref):
    m = compare.metrics(synthetic_history(ref), ref)
    assert m["n_points"] == len(ref.time)
    assert m["all"]["dV_max_ms"] == pytest.approx(5.0, abs=1e-6) and m["all"]["dV_rms_ms"] == pytest.approx(5.0, abs=1e-6)
    assert m["all"]["dh_max_m"] == pytest.approx(20.0, abs=1e-6)
    assert m["hypersonic"]["n_points"] < m["all"]["n_points"]
    assert m["hypersonic"]["dV_rel_max"] == pytest.approx(5.0 / 1000.0, abs=3e-4)   # worst at the last V_ref > 1 km/s
    assert m["d_end_time_s"] == 0.0 and m["final_velocity_model_ms"] == pytest.approx(ref.velocity[-1] + 5.0)


def test_alignment_uses_only_reference_times_inside_the_model_flight(ref):
    short = synthetic_history(ref, truncate=100)
    a = compare.align(short, ref)
    assert a["t"][-1] <= short.columns["time_s"][-1] and len(a["t"]) == 100
    m = compare.metrics(short, ref)
    assert m["d_end_time_s"] < 0.0


def test_plots_are_written(ref, tmp_path):
    paths = compare.plot_all(synthetic_history(ref), ref, str(tmp_path), "synthetic")
    assert [os.path.basename(p) for p in paths] == list(compare.PLOT_NAMES)
    for p in paths:
        assert os.path.getsize(p) > 5000
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py -q`
Expected: ERROR `cannot import name 'compare'`

- [ ] **Step 3: Implement `compare.py`**

```python
"""Model history vs a SESAM reference: the model sampled at the reference's own time stamps, error
metrics over the whole flight and the hypersonic phase (V_ref > 1 km/s), and six overlay plots."""
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))
try:
    from plot_style import INK, MUTED, SECOND, apply_rcparams, strip_top_right_spines
except ImportError:                                    # analysis/ not present: fall back to plain matplotlib
    INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"

    def apply_rcparams(plt, font_size=10.5):
        return None

    def strip_top_right_spines(ax):
        return None

MODEL_COLOR, REF_COLOR = "#3987e5", "#0b0b0b"
HYPERSONIC_MS = 1000.0
PLOT_NAMES = ("velocity_time.png", "altitude_time.png", "altitude_velocity.png", "angles_time.png",
              "ground_track.png", "regime_drag.png")


def align(history, reference):
    t_mod = history.columns["time_s"]
    mask = reference.time <= t_mod[-1]
    t = reference.time[mask]

    def at(col, scale=1.0):
        return np.interp(t, t_mod, history.columns[col]) * scale

    return {
        "t": t,
        "V_model": at("velocity_kms", 1e3), "V_ref": reference.velocity[mask],
        "h_model": at("altitude_km", 1e3), "h_ref": reference.altitude[mask],
        "gamma_model": at("flight_path_deg"), "gamma_ref": np.degrees(reference.flight_path[mask]),
        "heading_model": at("heading_deg"), "heading_ref": np.degrees(reference.heading[mask]),
        "lat_model": at("lat_deg"), "lat_ref": np.degrees(reference.lat[mask]),
        "lon_model": at("lon_deg"), "lon_ref": np.degrees(reference.lon[mask]),
        "kn_model": at("knudsen"), "kn_ref": reference.knudsen[mask],
        "cd_model": at("drag"), "cd_ref": reference.drag[mask],
    }


def _phase(a, mask):
    dV = a["V_model"][mask] - a["V_ref"][mask]
    dh = a["h_model"][mask] - a["h_ref"][mask]
    rel = np.abs(dV) / np.maximum(a["V_ref"][mask], 1.0)
    rms = lambda x: float(math.sqrt(np.mean(x * x))) if x.size else float("nan")
    return {"n_points": int(mask.sum()),
            "dV_max_ms": float(np.abs(dV).max()) if dV.size else float("nan"), "dV_rms_ms": rms(dV),
            "dV_rel_max": float(rel.max()) if rel.size else float("nan"), "dV_rel_rms": rms(rel),
            "dh_max_m": float(np.abs(dh).max()) if dh.size else float("nan"), "dh_rms_m": rms(dh)}


def metrics(history, reference):
    a = align(history, reference)
    t_model_end, t_ref_end = float(history.columns["time_s"][-1]), float(reference.time[-1])
    return {
        "n_points": int(a["t"].size),
        "model_end_time_s": t_model_end, "reference_end_time_s": t_ref_end,
        "d_end_time_s": t_model_end - t_ref_end, "d_end_time_rel": (t_model_end - t_ref_end) / t_ref_end,
        "final_velocity_model_ms": float(history.columns["velocity_kms"][-1] * 1e3),
        "final_velocity_reference_ms": float(reference.velocity[-1]),
        "all": _phase(a, np.ones(a["t"].size, dtype=bool)),
        "hypersonic": _phase(a, a["V_ref"] > HYPERSONIC_MS),
    }


def _overlay_with_residual(a, key, ylabel, resid_label, scale, path, title):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(a["t"], a[key + "_ref"] * scale, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(a["t"], a[key + "_model"] * scale, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_ylabel(ylabel); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    rx.plot(a["t"], (a[key + "_model"] - a[key + "_ref"]), color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel(resid_label); rx.set_xlabel("time [s]")
    for x in (ax, rx):
        strip_top_right_spines(x)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_all(history, reference, outdir, title):
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    a = align(history, reference)
    paths = [os.path.join(outdir, n) for n in PLOT_NAMES]
    _overlay_with_residual(a, "V", "velocity [km/s]", "model - SESAM [m/s]", 1e-3, paths[0], title)
    _overlay_with_residual(a, "h", "altitude [km]", "model - SESAM [m]", 1e-3, paths[1], title)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(a["V_ref"] / 1e3, a["h_ref"] / 1e3, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(a["V_model"] / 1e3, a["h_model"] / 1e3, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("velocity [km/s]"); ax.set_ylabel("altitude [km]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    strip_top_right_spines(ax); fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)

    fig, (g, hd) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for ax, key, lab in ((g, "gamma", "flight-path angle [deg]"), (hd, "heading", "heading [deg]")):
        ax.plot(a["t"], a[key + "_ref"], color=REF_COLOR, lw=1.6, label="SESAM")
        ax.plot(a["t"], a[key + "_model"], color=MODEL_COLOR, lw=1.2, ls="--", label="model")
        ax.set_ylabel(lab); strip_top_right_spines(ax)
    g.legend(frameon=False); g.set_title(title, color=SECOND, fontsize=10); hd.set_xlabel("time [s]")
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(a["lon_ref"], a["lat_ref"], color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(a["lon_model"], a["lat_model"], color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("longitude [deg]"); ax.set_ylabel("latitude [deg]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    strip_top_right_spines(ax); fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)

    fig, (kx, cx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    kx.semilogy(a["t"], np.maximum(a["kn_ref"], 1e-7), color=REF_COLOR, lw=1.6, label="SESAM")
    kx.semilogy(a["t"], np.maximum(a["kn_model"], 1e-7), color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    for thr in (10.0, 1.0, 0.1, 0.01):
        kx.axhline(thr, color=MUTED, lw=0.5, ls=":")
    kx.set_ylabel("Knudsen number"); kx.legend(frameon=False); kx.set_title(title, color=SECOND, fontsize=10)
    cx.plot(a["t"], a["cd_ref"], color=REF_COLOR, lw=1.6); cx.plot(a["t"], a["cd_model"], color=MODEL_COLOR, lw=1.2, ls="--")
    cx.set_ylabel("C_D"); cx.set_xlabel("time [s]")
    for x in (kx, cx):
        strip_top_right_spines(x)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    return paths
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py -q`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add reentry_model/compare.py tests/test_reentry_model_compare.py
git commit -m "reentry_model: metrics and overlay plots against SESAM references

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 9: Command line, package entry point, README section (`cli.py`, `__main__.py`)

**Files:**
- Create: `reentry_model/cli.py`, `reentry_model/__main__.py`, `tests/test_reentry_model_cli.py`
- Modify: `README.md` (new section before "## Tests")

**Interfaces:**
- Consumes: everything above.
- Produces: `cli.main(argv=None) -> int`; `cli.build_parser()`; `cli.model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name) -> str`; `cli.make_atmosphere(spec, epoch, wind_name) -> (atmosphere, name, info_dict)`; `cli.DEFAULT_OUTDIR`; the run JSON layout `{schema_version, run_name, inputs, settings, results, comparison, provenance}`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_cli.py`:
```python
"""cli.py: `run` and `compare` end to end on short flights."""
import csv
import json
import os

import pytest

from reentry_model import cli, compare, sesam_io

R100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind.csv")
BASE = ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331"]


def test_run_name():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none") == \
        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"


def test_run_us76_short_flight(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--t-max", "20", "--cadence", "5", "--outdir", str(tmp_path), "--quiet"])
    assert rc == 0
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none")
    rows = list(csv.DictReader(open(tmp_path / (name + ".csv"))))
    assert [float(r["time_s"]) for r in rows] == [0.0, 5.0, 10.0, 15.0, 20.0]
    doc = json.load(open(tmp_path / (name + ".json")))
    assert doc["settings"]["atmosphere"] == "us76" and doc["results"]["end_reason"] == "t_max"
    assert doc["inputs"]["mass_kg"] == pytest.approx(1.4728833557580150) and doc["comparison"] is None
    assert doc["provenance"]["package_version"] and "git_commit" in doc["provenance"]


def test_run_replay_with_reference_writes_comparison_and_plots(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "replay:" + R100, "--reference", R100, "--t-max", "30", "--cadence", "10",
                          "--outdir", str(tmp_path), "--name", "r100_short", "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "r100_short.json"))
    m = doc["comparison"]["metrics"]
    assert doc["settings"]["atmosphere"] == "replay" and doc["settings"]["replay_reference"].endswith("_nowind")
    assert m["all"]["n_points"] >= 25 and m["hypersonic"]["dV_max_ms"] < 50.0     # 30 s of a matched flight
    for plot in compare.PLOT_NAMES:
        assert os.path.isfile(tmp_path / "r100_short" / plot)


def test_compare_subcommand(tmp_path):
    cli.main(BASE + ["--atmosphere", "us76", "--t-max", "30", "--cadence", "1", "--outdir", str(tmp_path), "--name", "m", "--quiet"])
    rc = cli.main(["compare", "--model", str(tmp_path / "m.csv"), "--reference", R100, "--outdir", str(tmp_path / "cmp"), "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "cmp" / "m_vs_reference.json"))
    assert doc["reference"].endswith("_nowind") and "hypersonic" in doc["metrics"]


@pytest.mark.parametrize("argv", [
    BASE + ["--atmosphere", "gram"],
    BASE + ["--bridging", "legge"],
    BASE + ["--gravity", "j6"],
    ["run", "--diameter", "-1", "--velocity", "7.5", "--altitude", "77.5"],
    BASE + ["--epoch", "yesterday"],
])
def test_bad_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
```

- [ ] **Step 2: Run to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: ERROR `cannot import name 'cli'`

- [ ] **Step 3: Implement `cli.py` and `__main__.py`**

`reentry_model/__main__.py`:
```python
import sys

from .cli import main

sys.exit(main())
```

`reentry_model/cli.py`:
```python
"""Command line of the re-entry model.

    python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
        [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
    python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]
"""
import argparse
import math
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from . import __version__, aero, atmosphere, body, compare, fap, sesam_io
from . import trajectory as tj
from .earth import GRAVITY_MODELS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output")
# the reference cases' break-off state (parent satellite, see README "Reference runs")
DEFAULT_HEADING_DEG = 347.168296
DEFAULT_LAT_DEG = 29.546067
DEFAULT_LON_DEG = -82.134333
DEFAULT_EPOCH = "2024-08-01T12:53:07"
DEFAULT_MATERIAL_DENSITY = 2813.0          # drama-AA7075
ATMOSPHERES = ("nrlmsise", "us76")         # plus replay:<path>
WINDS = ("none", "static")


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))


def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name):
    return "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
        diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)


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
    return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0]}


def build_parser():
    p = argparse.ArgumentParser(prog="reentry_model", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="integrate one sphere trajectory")
    r.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    r.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s], relative to the rotating atmosphere")
    r.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    r.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending")
    r.add_argument("--heading", type=float, default=DEFAULT_HEADING_DEG, help="[deg] clockwise from north (default %(default)s)")
    r.add_argument("--lat", type=float, default=DEFAULT_LAT_DEG, help="geodetic latitude [deg] (default %(default)s)")
    r.add_argument("--lon", type=float, default=DEFAULT_LON_DEG, help="longitude [deg] (default %(default)s)")
    r.add_argument("--epoch", type=parse_epoch, default=parse_epoch(DEFAULT_EPOCH), help="UTC YYYY-MM-DDTHH:MM:SS (default {})".format(DEFAULT_EPOCH))
    r.add_argument("--material-density", type=float, default=DEFAULT_MATERIAL_DENSITY, help="[kg/m3] (default %(default)s)")
    r.add_argument("--temperature", type=float, default=300.0, help="initial temperature [K], recorded only (default %(default)s)")
    r.add_argument("--atmosphere", default="nrlmsise", help="nrlmsise (default) | us76 | replay:<sesam.csv>")
    r.add_argument("--wind", choices=WINDS, default="none")
    r.add_argument("--bridging", choices=aero.BRIDGING_NAMES, default="sesam-table")
    r.add_argument("--gravity", choices=GRAVITY_MODELS, default="j2")
    r.add_argument("--rtol", type=float, default=1e-9)
    r.add_argument("--cadence", type=float, default=1.0, help="history sample spacing [s] (default %(default)s)")
    r.add_argument("--t-max", type=float, default=3600.0, help="[s] (default %(default)s)")
    r.add_argument("--reference", default=None, help="SESAM run CSV to compare against (also samples the model at its times)")
    r.add_argument("--outdir", default=DEFAULT_OUTDIR)
    r.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>)")
    r.add_argument("--quiet", action="store_true")

    c = sub.add_parser("compare", help="metrics and plots for an existing model history")
    c.add_argument("--model", required=True, help="model history CSV")
    c.add_argument("--reference", required=True, help="SESAM run CSV")
    c.add_argument("--outdir", default=DEFAULT_OUTDIR)
    c.add_argument("--title", default=None)
    c.add_argument("--quiet", action="store_true")
    return p


def cmd_run(args, parser):
    for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                         ("--cadence", args.cadence), ("--t-max", args.t_max)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.atmosphere not in ATMOSPHERES and not args.atmosphere.startswith("replay:"):
        parser.error("--atmosphere must be one of {} or replay:<sesam.csv>".format(ATMOSPHERES))

    initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                              flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
                              lat=math.radians(args.lat), lon=math.radians(args.lon), epoch=args.epoch)
    settings = tj.Settings(diameter=args.diameter * 1e-3, gravity=args.gravity, rtol=args.rtol,
                           cadence=args.cadence, t_max=args.t_max)
    mass = args.material_density * 4.0 / 3.0 * math.pi * (settings.diameter / 2.0) ** 3
    atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
    reference = sesam_io.load_reference(args.reference) if args.reference else None
    sim = tj.Simulator(initial, body.ConstantBody(mass, args.temperature), atm,
                       aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
    history = sim.run(extra_times=reference.time if reference is not None else None)

    name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind)
    os.makedirs(args.outdir, exist_ok=True)
    csv_path = os.path.join(args.outdir, name + ".csv")
    json_path = os.path.join(args.outdir, name + ".json")
    tj.write_history_csv(history, csv_path)
    doc = {
        "schema_version": 1,
        "run_name": name,
        "inputs": {"diameter_mm": args.diameter, "initial_velocity_kms": args.velocity, "initial_altitude_km": args.altitude,
                   "flight_path_angle_deg": args.flight_path_angle, "heading_deg": args.heading, "latitude_deg": args.lat,
                   "longitude_deg": args.lon, "epoch_utc": args.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
                   "material_density_kgm3": args.material_density, "mass_kg": mass, "initial_temperature_K": args.temperature},
        "settings": {"atmosphere": atm_name, "wind": args.wind, "bridging": args.bridging, "gravity": args.gravity,
                     "rtol": args.rtol, "cadence_s": args.cadence, "t_max_s": args.t_max, **atm_info},
        "results": history.results,
        "comparison": None,
        "provenance": provenance(),
        "files": {"csv": os.path.abspath(csv_path)},
    }
    if reference is not None:
        plots = compare.plot_all(history, reference, os.path.join(args.outdir, name), name)
        doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
                             "reference_sha256": reference.sha256, "metrics": compare.metrics(history, reference),
                             "plots": [os.path.abspath(p) for p in plots]}
    tj.write_run_json(json_path, doc)
    if not args.quiet:
        res = history.results
        print("{}: {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
            name, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
            history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
        if reference is not None:
            hyp = doc["comparison"]["metrics"]["hypersonic"]
            print("  vs {}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m; end time {:+.1f} s".format(
                reference.name, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], doc["comparison"]["metrics"]["d_end_time_s"]))
        print("  csv  -> {}\n  json -> {}".format(os.path.abspath(csv_path), os.path.abspath(json_path)))
    return 0


def cmd_compare(args, parser):
    history = tj.read_history_csv(args.model)
    reference = sesam_io.load_reference(args.reference)
    stem = os.path.splitext(os.path.basename(args.model))[0]
    title = args.title or "{} vs {}".format(stem, reference.name)
    os.makedirs(args.outdir, exist_ok=True)
    plots = compare.plot_all(history, reference, args.outdir, title)
    doc = {"model": os.path.abspath(args.model), "reference": reference.name, "reference_sha256": reference.sha256,
           "metrics": compare.metrics(history, reference), "plots": [os.path.abspath(p) for p in plots]}
    out = os.path.join(args.outdir, stem + "_vs_reference.json")
    tj.write_run_json(out, doc)
    if not args.quiet:
        hyp = doc["metrics"]["hypersonic"]
        print("{}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m -> {}".format(title, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], out))
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            return cmd_run(args, parser)
        return cmd_compare(args, parser)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: 9 passed (5 parametrized bad-argument cases). Note `test_bad_arguments_exit_2` expects `SystemExit(2)` from `parser.error`, while a bad `--atmosphere` string is also caught by `parser.error` before any work is done.

- [ ] **Step 5: Add the README section**

Insert before the `## Tests` heading of `README.md`:
```markdown
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
`--wind none|static`; `--bridging sesam-table|sesam-erf|sin2|textbook`; `--gravity point|j2|j2j4`. Outputs go to
`reentry_model_output/` (git-ignored): `<run>.csv` with the same columns as the SESAM histories, `<run>.json`
(inputs, settings, results, comparison metrics, provenance), and with `--reference` a folder of six plots
(V(t) and h(t) overlays with residuals, h(V), angles, ground track, Knudsen/C_D). `data/reference_runs/` holds the
four committed SESAM references (100 mm from 77.5 km, 50 mm from 115 km, winds on/off).
```

- [ ] **Step 6: Run the whole suite and commit**

Run: `"$PY" -m pytest -m "not drama" -q`
Expected: all green (the wrapper's 207 plus the model's tests).

```bash
git add reentry_model/cli.py reentry_model/__main__.py tests/test_reentry_model_cli.py README.md
git commit -m "reentry_model: run/compare command line, package entry point, README section

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 10: Verification against the four SESAM references

**Files:**
- Create: `analysis/reentry_model_verification.py`, `tests/test_reentry_model_reference.py`
- Modify: `README.md` (results table under the model section)

**Interfaces:**
- Consumes: `cli.main`, the reference runs, `compare.metrics`.
- Produces: `reentry_model_output/verification/summary.md` and `summary.json`; the `reference`-marked integration test with the confirmed thresholds.

- [ ] **Step 1: Write the verification script**

Create `analysis/reentry_model_verification.py`:
```python
#!/usr/bin/env python3
"""Run reentry_model against the four committed SESAM reference runs, twice each (SESAM density replay =
dynamics only; NRLMSISE-00 = full model), and tabulate the hypersonic-phase errors.

    "$PY" analysis/reentry_model_verification.py [--outdir reentry_model_output/verification]
"""
import argparse
import json
import math
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

NAMES = [
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind",
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis",
]
MODES = ("replay", "nrlmsise")


def run_case(name, mode, outdir):
    csv_path = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    ref = sesam_io.load_reference(csv_path)
    run_name = "{}__{}".format(name, mode)
    argv = ["run", "--diameter", "{:.6g}".format(ref.diameter * 1e3), "--velocity", "{:.6g}".format(ref.initial.velocity / 1e3),
            "--altitude", "{:.9g}".format(ref.initial.altitude / 1e3), "--flight-path-angle", "{:.9g}".format(math.degrees(ref.initial.flight_path)),
            "--heading", "{:.9g}".format(math.degrees(ref.initial.heading)), "--lat", "{:.9g}".format(math.degrees(ref.initial.lat)),
            "--lon", "{:.9g}".format(math.degrees(ref.initial.lon)), "--epoch", ref.initial.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
            "--material-density", "{:.6g}".format(ref.material_density),
            "--atmosphere", "replay:" + csv_path if mode == "replay" else "nrlmsise",
            "--reference", csv_path, "--outdir", outdir, "--name", run_name, "--quiet"]
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run failed for {} ({}): exit {}".format(name, mode, rc))
    with open(os.path.join(outdir, run_name + ".json")) as fh:
        doc = json.load(fh)
    m = doc["comparison"]["metrics"]
    return {"case": name, "mode": mode, "winds_in_reference": ref.use_wind,
            "hyp_dV_rel_max": m["hypersonic"]["dV_rel_max"], "hyp_dV_max_ms": m["hypersonic"]["dV_max_ms"],
            "hyp_dh_max_m": m["hypersonic"]["dh_max_m"], "all_dV_max_ms": m["all"]["dV_max_ms"], "all_dh_max_m": m["all"]["dh_max_m"],
            "d_end_time_s": m["d_end_time_s"], "d_end_time_rel": m["d_end_time_rel"],
            "runtime_s": doc["results"]["runtime_s"], "rhs_evaluations": doc["results"]["rhs_evaluations"]}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification"))
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    rows = [run_case(name, mode, args.outdir) for name in NAMES for mode in MODES]
    lines = ["| case | mode | winds in ref | hypersonic max dV | hypersonic max dh | whole-flight max dV / dh | end time | runtime |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append("| {} | {} | {} | {:.1f} m/s ({:.3%}) | {:.0f} m | {:.1f} m/s / {:.0f} m | {:+.1f} s ({:+.2%}) | {:.0f} s, {} evals |".format(
            r["case"].replace("sphere_", "").replace("_T0300.0K_v07.50000kms", "").replace("_mAA7075_nomelt_msis", ""), r["mode"],
            "on" if r["winds_in_reference"] else "off", r["hyp_dV_max_ms"], r["hyp_dV_rel_max"], r["hyp_dh_max_m"],
            r["all_dV_max_ms"], r["all_dh_max_m"], r["d_end_time_s"], r["d_end_time_rel"], r["runtime_s"], r["rhs_evaluations"]))
    table = "\n".join(lines)
    with open(os.path.join(args.outdir, "summary.md"), "w") as fh:
        fh.write(table + "\n")
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    print(table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run it**

Run: `"$PY" analysis/reentry_model_verification.py`
Expected: eight runs, a printed table, and `reentry_model_output/verification/summary.md`. Read the table against the spec's expectations (section 9): replay ≤ 0.2 % / 100 m over the hypersonic phase and ≤ 1 % in end time; NRLMSISE-00 ≤ 1 % / 500 m; winds-on references at most ~8 m/s / 20 m worse than winds-off.

- [ ] **Step 3: If an expectation is missed, find the cause before touching thresholds**

Use the plots in `reentry_model_output/verification/<run>/` and check, in this order (each is a one-line experiment via the CLI flags):
1. A constant velocity offset from t = 0 in replay mode: SESAM's `velocity` may be airspeed vs our ground-relative speed; only possible with winds in the reference — compare the winds-off case first.
2. A residual growing with time in replay mode: gravity (`--gravity point` vs `j2` vs `j2j4`) or the Ma < 5 drag clamp (visible only after the hypersonic phase).
3. A density-shaped residual in NRLMSISE-00 mode: the F10.7 convention (`fap.solar_indices` uses the previous day's F10.7; try the same day) or `m_bar` above 100 km — the atmosphere unit test `TestNRLMSISE00` reports the density ratio at both reference starts.
Record what was found in the README table's notes. Change a threshold only with the measured number and the reason written next to it.

- [ ] **Step 4: Write the integration test with the confirmed thresholds**

Create `tests/test_reentry_model_reference.py`:
```python
"""Full-model verification against the committed SESAM references (marker: reference)."""
import os

import pytest

from reentry_model import aero, atmosphere, body, cli, fap, sesam_io
from reentry_model import compare
from reentry_model import trajectory as tj

NAMES = [
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind",
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis",
]
# Spec section 9 expectations over the hypersonic phase (V_ref > 1 km/s). Task 10 step 3 may revise a
# value only together with the measured number and the reason, recorded in the README verification table.
THRESHOLDS = {
    "replay":   {"dV_rel_max": 0.002, "dh_max_m": 100.0, "d_end_time_rel": 0.01},
    "nrlmsise": {"dV_rel_max": 0.010, "dh_max_m": 500.0, "d_end_time_rel": 0.03},
}


def simulate(ref, mode):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    if mode == "replay":
        atm = atmosphere.ReplayAtmosphere(ref)
    else:
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), ref.initial.epoch.date())
        atm = atmosphere.NRLMSISE00Atmosphere(ref.initial.epoch, solar)
    mass = ref.material_density * 4.0 / 3.0 * 3.141592653589793 * (ref.diameter / 2.0) ** 3
    sim = tj.Simulator(initial, body.ConstantBody(mass), atm, aero.SphereDragTables.from_json(), aero.SesamTable(),
                       tj.Settings(diameter=ref.diameter, cadence=5.0))
    return sim.run(extra_times=ref.time)


@pytest.mark.reference
@pytest.mark.parametrize("mode", ["replay", "nrlmsise"])
@pytest.mark.parametrize("name", NAMES)
def test_model_matches_sesam(name, mode):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, name + ".csv"))
    history = simulate(ref, mode)
    assert history.end_reason == "ground"
    m = compare.metrics(history, ref)
    hyp, thr = m["hypersonic"], THRESHOLDS[mode]
    assert hyp["dV_rel_max"] <= thr["dV_rel_max"], (name, mode, hyp)
    assert hyp["dh_max_m"] <= thr["dh_max_m"], (name, mode, hyp)
    assert abs(m["d_end_time_rel"]) <= thr["d_end_time_rel"], (name, mode, m["d_end_time_s"])
```

- [ ] **Step 5: Run the integration test and the whole suite**

Run: `"$PY" -m pytest -m reference -q` then `"$PY" -m pytest -m "not drama" -q`
Expected: 8 passed for the reference marker; everything green overall.

- [ ] **Step 6: Record the results in the README**

Append under the model section of `README.md`, replacing the placeholders with the numbers from `summary.md`:
```markdown
### Verification (Task 10, `analysis/reentry_model_verification.py`)

Model sampled at SESAM's own time stamps; errors over the hypersonic phase (V > 1 km/s). "replay" feeds SESAM's
density/temperature into the model (dynamics and drag only); "nrlmsise" is the full model. The winds-on
references are compared with the wind-free model, so their rows include SESAM's HWM14 wind effect (≤ 8 m/s).

<paste the table from reentry_model_output/verification/summary.md>

Notes: <what step 3 found, one line each; "none" if every expectation held>.
```

- [ ] **Step 7: Commit**

```bash
git add analysis/reentry_model_verification.py tests/test_reentry_model_reference.py README.md
git commit -m "reentry_model: verification against the four SESAM references with recorded thresholds

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Plan self-review

**Spec coverage.** §2 facts → Task 1 (data), 4 (drag), 6 (atmosphere), 7 (dynamics). §3 references → Task 1 copies them, Task 10 uses all four. §4 scope → Tasks 6–9 implement every "in" item; the deferred items appear only as interfaces (`Matting`, `StaticProfileWind`, `Body`). §5 layout → matches the file structure above. §6.1–6.2 → Task 2 + 7 (ECEF form, J2/J4, dγ/dt test). §6.3 → Task 6 (three atmospheres, winds; fap inputs from Task 3). §6.4–6.5 → Task 4 (Kn/Ma definitions, tables, clamping, bridging objects; reference reproduction test). §6.6 → Task 7 (DOP853, tolerances, events, cadence + reference-time sampling). §7 CLI → Task 9. §8 outputs → Task 7 (CSV columns, precision), Task 9 (JSON layout), Task 8 (metrics, six plots). §9 verification → Task 10 (sequence, expectations, "confirm or justify"). §10 tests → Tasks 1–9 unit tests, Task 10 integration test with the `reference` marker registered in Task 1. §11 environment → Task 1 step 1.

**Placeholders.** None: every step carries code or an exact command; the two "if it fails" notes name the specific checks to make.

**Type consistency.** `Freestream(rho, T, p, m_bar, wind_enu)` is produced by Task 6 and consumed by Task 7's `aero_state`; `Reference` fields used in Tasks 6, 8, 9, 10 all exist in Task 5's dataclass (`time, altitude, velocity, lat, lon, flight_path, heading, knudsen, mach, density, drag, diameter, material_density, initial.*, use_wind, sha256, name, csv_path`); `History(columns, states, end_reason, results)` is constructed the same way in Tasks 7 and 8; `Settings` keyword names in the tests (`gravity, rotating_frame, cadence, t_max, escape_altitude, rtol`) match Task 7; `cli.model_run_name` argument order is identical in Task 9's tests and implementation; `aero.BRIDGING_NAMES` and `earth.GRAVITY_MODELS` are the `choices` used by the CLI.
