# Coupled Trajectory + 3D FEM Heat Transfer of a Sphere (Step 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `reentry_model` so that the Step 1 trajectory and a three-dimensional finite-element temperature field of the sphere are solved together — aerothermal heating on every surface patch (SESAM-equivalent verification mode and a physics mode), radiation, P1-tetrahedra conduction — with VTK/animation output and a verification against SESAM's heat input and lumped temperature on the two no-melt US76 references.

**Architecture:** New modules `material`, `mesh`, `thermal/` (protocol + scikit-fem backend + FEniCSx backend), `gas`, `heating`, `coupled`, `viz`; `body.py` becomes a stepper protocol with `ThermalBody`; `trajectory.Simulator` gains `advance(dt)`; `compare`, `sesam_io` and `cli` grow the heat/temperature columns, metrics, plots and flags. Per 0.5 s macro step: trajectory → aero state → per-patch heating → backward-Euler conduction step with Newton on the radiation term (secant heat capacity, precomputed element matrices, AMG-preconditioned CG) → history row / VTK frame. Every physics module is a small file with a plain-array interface and its own tests; the wrapper (`sphere_reentry.py`) and the Step 1 modules keep their meaning.

**Tech Stack:** Python 3.12 in `drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`); numpy 2.5, scipy 1.18, pymsis 0.13 (Step 1); new: scikit-fem 12.0.2, gmsh 4.15.2, meshio 5.3.5, pyamg 5.3.0, pyvista 0.49.0 (VTK 9.7), imageio 2.37.4 + imageio-ffmpeg 0.6.0, cantera 3.2.0 (all pip wheels, arm64); FEniCSx (dolfinx ≥ 0.9) only in a separate conda env `fenicsx_env`, created only on the user's explicit go-ahead.

**Spec:** `docs/superpowers/specs/2026-09-18-thermal-fem-design.md` — read it first. Facts it relies on: `Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md` §§1–14 (Task 12 appends §15). The measured deviations from the spec listed below are binding; each is recorded in the spec by Task 12.

## Global Constraints

- Interpreter for everything: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` (never the system python). Unit tests: `"$PY" -m pytest -m "not drama and not reference" -q` (Step 1: 278 passed in ~60 s; this plan adds ~120 tests and ~3 min, dominated by the coupled and thermal solver tests). Reference tier: `"$PY" -m pytest -m reference -q` (~35 min after Task 12).
- SI units inside the package (m, s, kg, K, W, rad); the CLI/CSV keep the wrapper's units (km, km/s, deg, mm) plus SI heat columns (W, J, K). Angles in radians; conversion only in `cli.py`, `sesam_io.py` and the CSV writer.
- The SESAM-equivalent heating distribution is uniform over the surface by construction and **not physical**; the class docstring, a comment at the line applying the 0.27471 factor, the CLI help and the README must say so (spec §6.2 — user requirement).
- Never modify DRAMA's databases or the wrapper (`sphere_reentry.py`, `sphere_sweep.py`); the material for the thermal model is a copy of `data/user_materials/AA7075_nomelt.json` inside the package (a test keeps the copy identical).
- `drama_env` gains the packages of `requirements-step2.txt` (Task 1). The FEniCSx environment `fenicsx_env` is **not** created by this plan: Task 11 writes the backend and its tests skip when `dolfinx` is not importable; the CLI exits 2 with a message naming `fenicsx_env`.
- Git-ignored output root `reentry_model_output/` (already in `.gitignore`); meshes are cached under `reentry_model_output/meshes/`; never write elsewhere by default.
- Repo conventions: module docstrings, `argparse`, exit codes 0/1/2, tests under `tests/` named `test_reentry_model_<module>.py`, one commit per task, commit messages in the imperative like the existing history, ending with `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
- Every task's tests are run with the exact command given in the task; a task is done only when the whole unit tier passes (`"$PY" -m pytest -m "not drama and not reference" -q`).
- Code in this plan was executed and its tests passed on 2026-09-18 in a throwaway copy of the package (scikit-fem 12.0.2 etc. as pinned). Transcribe it verbatim; where a task says "replace the file", the block is the complete new file.

## Measured facts and spec amendments (2026-09-18, prototype in a throwaway copy)

These were measured while writing the plan and override the corresponding spec statements; Task 12 writes them into the spec and the facts note.

1. **Cantera `airNASA9.yaml` has no transport data.** Thermodynamics/composition come from it (NASA-9, 200–20 000 K, ions included); viscosity comes from Blottner curve fits (N2, O2, NO, N, O; Gnoffo, Gupta & Shinn, NASA TP-2867, 1989) with Wilke's mixing rule (N2 at 300 K: 1.786e-5 Pa s). Spec §6.3 "which also gives μ_s, μ_w" is amended accordingly.
2. **Fay–Riddell vs Sutton–Graves at the 100 mm start** (ρ 2.727e-5 kg/m³, T∞ 203.6 K, 7.5 km/s, R 0.05 m, T_w 300 K): equilibrium stagnation state T_s 5820 K, p_s 1495 Pa, ρ_s 5.20e-4 kg/m³, h_s 2.803e7 J/kg, μ_s 1.69e-4 Pa s, X_N 0.59, h_D/h_s 0.718; Fay–Riddell 2.22e6 W/m² = 1.29 × Sutton–Graves (1.716e6) and 1.13 × DKR (1.957e6); with the Lewis term off (`catalycity` 0) 1.14 × Sutton–Graves. The spec's "within 15 % of Sutton–Graves" did not hold (72 % of h_s is dissociation enthalpy at 1.5 kPa; Sutton–Graves is a Le = 1 fit); the test asserts the measured band 1.2–1.4. Over the flight the ratio falls to 1.25 (60 km, 7 km/s) and 1.11 (4 km/s).
3. **Mesh sizes.** gmsh Distance/Threshold grading on the 100 mm sphere: 1 mm/8 mm = 76 289 nodes (spec estimated ≈40 k), 2 mm/8 mm = 18 896, 2 mm/16 mm = 12 643, 4 mm/20 mm = 2 934, uniform 4 mm = 6 663. **Defaults become h_surface = 2 mm, h_core = 8 mm** (18.9 k nodes); the convergence run halves h_surface (1 mm/8 mm, 76 k nodes). A node is embedded at the centre (gmsh `embed`) so the centre temperature is a nodal value.
4. **Solver cost** (12.6 k nodes): scikit-fem's generic `asm` 0.23 s per Newton iteration vs 0.02 s for rescaling precomputed element matrices into a fixed CSR pattern (identical to 1e-15, kept as the conformance test); SuperLU 0.41 s per solve vs pyamg SA-CG 0.03 s (76 k nodes: 24 s vs 0.2 s). **`--linear-solver amg` is the default**; `direct` is for tests/small meshes. The AMG hierarchy is rebuilt every 30 solves and the Newton start is extrapolated from the previous step (2.0–2.1 iterations per step at tol 1e-6). Full 100 mm flight on the default mesh: 732 macro steps, **170 s** (SESAM-equivalent) / **300 s** (physics incl. animation); 50 mm from 115 km: 1117 steps, 52 s.
5. **Secant heat capacity.** The mass coefficient is c = [h(T_k) − h(T_old)]/(T_k − T_old) per element, so the discrete energy balance ΔE = (Q_conv − Q_rad)Δt holds to the Newton tolerance for the tabulated c_p(T) as well (measured 1e-8 over full flights), not only for constant c_p; this is also the hook for Step 3's latent heat.
6. **Analytic tests.** Lumped limit and radiative cooling use k × 1e4 (k × 1e6 makes CG indefinite) and agree with the lumped ODE to 8e-5 / 1.7e-4; Carslaw–Jaeger needs a uniform 4 mm test mesh (the graded default's 8–20 mm core puts the centre 5–7 % off): centre 0.5 %, volume mean 0.2 % at 2–4 s.
7. **SESAM's heating has a hot-wall factor.** Q_SESAM/(A · 0.27471 · q_DKR) = max(0, 1 − c_p(T_lumped − T∞)/(V²/2)) with c_p = 1004.5 J/kg/K (best fit 985; rms 0.008, max 0.015 over 189 rows to 2450 K), i.e. it is 1 only while the sphere is cold — the facts note's "no hot-wall correction detectable" was measured below 850 K where the factor is > 0.97. It clamps to exactly 0 (Q = 0 while V < ~1.8 km/s and T ~ 2000 K) and is 0.5 with no hot-wall term below Ma 1 (like SESAM's C_D rule).
8. **SESAM's transitional heating is not the drag blend.** Its factor F_h(Kn) = Q/(A · 0.27471 · q_DKR · hot-wall) is 1.005 in the continuum, 0.94 at Kn 0.04, 0.46 at Kn 0.2, 0.14 at Kn 1 (where (1 − f)q_DKR + f·½ρV³ would give 0.5) and 0.059 at Kn 40 (= 0.78 × the cos θ average of ½ρV³ — so SESAM's free-molecular limit is the textbook value with α ≈ 0.8, not "13× too low"; the facts note §4 decomposition, measured at Kn 0.03, was misattributing the transitional deficit to q_FM). Measured in 0.125-decade bins (`aero.SesamHeatTable`, rms 0.7 %, max 3.4 %). **The SESAM-equivalent mode is therefore q = 0.27471 · q_DKR · F_h(Kn) · hot-wall on every patch** (spec §6.2's `(1 − f)·q_c + f·½αρV³` form gave +12 % / +17 % integrated heat on the 100 / 50 mm references; the measured form gives +0.3 % / +0.03 %).
9. **Verification results with 8** (SESAM-equivalent, default mesh, Δt 0.5 s, US76, winds off): 100 mm — Q_conv max 0.46 % of peak, 2.2 % point-wise (continuum, Q_ref > 10 % of peak), integrated heat +0.31 % (end of hypersonic phase), |ΔT_eq| max 24.5 K (1.2 %), radiated 6.7 % of peak; 50 mm — 0.37 %, 2.4 %, +0.03 %, 28 K (1.2 %), 3.4 %. Physics mode integrated heat / SESAM = 0.744 (100 mm, peak stagnation T 2319 K vs mean 1717 K) and 0.77 (50 mm).
10. **Metric definitions** (spec §9): convective and radiated power errors are normalised by SESAM's peak over the hypersonic phase; the point-wise relative Q_conv error is restricted to Kn_ref < 0.01 **and Q_ref > 10 % of its peak** (SESAM's clamp makes relative errors next to Q = 0 unbounded); the radiated-power threshold is 8 % (= 4 × the 2 % temperature threshold, because the resolved surface radiates at its own hotter temperature); everything else keeps the spec's 3 % / 2 %.
11. **History columns** gain `absorbed_heat_J` (∫(Q_conv − Q_rad)dt) and `T_centre_K` next to the spec's list; `integrated_heat_J` keeps SESAM's meaning (∫Q_conv dt). Run layout: `<outdir>/<name>.csv/.json` as in Step 1, everything else (plots, `vtk/` series, frames, animation, stills) under `<outdir>/<name>/`.
12. Cantera's equilibrium solver is skipped when T∞ + V²/(2c_p) < 1500 K (frozen air; avoids `ChemEquil` warnings below 298 K at subsonic speeds); T∞ is clamped to ≥ 200 K (NASA-9 floor).

---

## File structure

```
requirements-step2.txt                         pinned Step 2 packages for drama_env
reentry_model/material.py                      Material: rho, k(T), c_p(T), emissivity, exact enthalpy h(T) and inverse, latent-heat hook
reentry_model/data/materials/AA7075_nomelt.json  copy of data/user_materials/AA7075_nomelt.json
reentry_model/mesh.py                          gmsh graded sphere (cached), loader, VolumeMesh (mutable points, tets), SurfaceMesh (patches)
reentry_model/thermal/__init__.py              ThermalSolver protocol, StepResult, MissingBackend, thermal_solver(name)
reentry_model/thermal/skfem_backend.py         P1 tets, backward Euler, secant c_p, Newton on radiation, fixed-pattern assembly, AMG/direct
reentry_model/thermal/fenicsx_backend.py       same scheme in UFL/dolfinx (lazy import; stub in Task 3, full in Task 11)
reentry_model/gas.py                           EquilibriumAir (Cantera airNASA9), Blottner/Wilke viscosity, air enthalpy table, DKR, Sutton-Graves, Fay-Riddell
reentry_model/aero.py                          + SesamHeatTable (SESAM's measured heat bridging F_h(Kn))
reentry_model/heating.py                       HeatingResult, lees_shape, cosine_shape, matting_bridge, SesamEquivalentHeating, PhysicsHeating, TabulatedHeating slot
reentry_model/body.py                          Body stepper protocol, ConstantBody, ThermalBody (bookkeeping, surface stats)
reentry_model/trajectory.py                    Simulator.advance(dt), sample_row(), results(); CSV IO of every column
reentry_model/coupled.py                       CoupledSettings, CoupledRun (lockstep loop), THERMAL_COLUMNS, VTK frames + PVD collections
reentry_model/viz.py                           read_series, render_frame, frame_title, animate (MP4/GIF + stills)
reentry_model/compare.py                       + has_thermal, align_thermal, thermal_metrics, plot_thermal, THERMAL_PLOT_NAMES
reentry_model/sesam_io.py                      + convective_heat, rad_cooling, integrated_heat on Reference
reentry_model/cli.py                           + thermal flags, build_thermal(), coupled run path, exit 2 for missing backends
analysis/reentry_model_thermal_verification.py runs both spheres x both heating modes, tabulates the thermal metrics
tests/conftest.py                              + session fixtures coarse_sphere_mesh, uniform_test_mesh
tests/test_reentry_model_{material,mesh,thermal,gas,heating,coupled,viz,fenicsx,reference_thermal}.py   new
tests/test_reentry_model_{trajectory,compare,cli}.py   extended
README.md, docs/superpowers/specs/2026-09-18-thermal-fem-design.md, sesam_verified_facts.md (§15)   Task 12
```

Dependency direction (spec §4): `coupled` → `trajectory`, `heating`, `body`; `heating` → `gas`, `aero`; `thermal/*` → `mesh`, `material` and their own library; `viz` and `compare` read exported files only.

---
### Task 1: Environment and the material model

**Files:**
- Create: `requirements-step2.txt`
- Create: `reentry_model/data/materials/AA7075_nomelt.json` (copy of `data/user_materials/AA7075_nomelt.json`)
- Create: `reentry_model/material.py`
- Test: `tests/test_reentry_model_material.py`

**Interfaces:**
- Consumes: `reentry_model.DATA_DIR` (package data directory, Step 1).
- Produces: `Material` dataclass with `rho`, `emissivity`, `k(T)`, `cp(T)`, `cp_eff(T)`, `liquid_fraction(T)`, `enthalpy(T)` (J/kg above `T_REF` = 293 K, exact integral of the piecewise-linear c_p), `temperature_from_enthalpy(h)`; `Material.from_drama_json(path=None)`; module constants `DEFAULT_MATERIAL`, `T_REF`. All methods accept scalars or arrays.

- [ ] **Step 1: Pin and install the Step 2 packages in drama_env**

Create `requirements-step2.txt`:

```
scikit-fem==12.0.2
gmsh==4.15.2
meshio==5.3.5
pyamg==5.3.0
pyvista==0.49.0
imageio==2.37.4
imageio-ffmpeg==0.6.0
cantera==3.2.0
```

Run:

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
"$PY" -m pip install -r requirements-step2.txt
"$PY" -c "import skfem, gmsh, meshio, pyamg, pyvista, imageio, imageio_ffmpeg, cantera; print(skfem.__version__, gmsh.__version__, pyvista.__version__, cantera.__version__)"
```

Expected: `12.0.2 4.15.2 0.49.0 3.2.0` (vtk 9.7.0 comes with pyvista). If a wheel is missing for this platform, stop and report — do not substitute versions silently.

- [ ] **Step 2: Copy the material file into the package**

```bash
mkdir -p reentry_model/data/materials
cp data/user_materials/AA7075_nomelt.json reentry_model/data/materials/AA7075_nomelt.json
```

- [ ] **Step 3: Write the failing tests**

Create `tests/test_reentry_model_material.py`:

```python
"""material.py: DRAMA material tables, enthalpy and its inverse."""
import filecmp
import os

import numpy as np
import pytest

from reentry_model import material
from helpers import REPO_ROOT


def test_packaged_material_is_the_wrapper_file():
    assert filecmp.cmp(material.DEFAULT_MATERIAL, os.path.join(REPO_ROOT, "data", "user_materials", "AA7075_nomelt.json"), shallow=False)


def test_default_material_tables():
    m = material.Material.from_drama_json()
    assert m.name == "AA7075_nomelt" and m.rho == 2813.0 and m.emissivity == 0.4
    assert m.cp(293.0) == 877.5 and m.cp(850.0) == 1131.6 and m.cp(2000.0) == 1131.6       # held above 850 K
    assert m.k(293.0) == 163.89 and m.k(850.0) == 128.19 and m.k(2000.0) == 128.19
    assert m.cp(303.0) == pytest.approx(881.5) and m.k(100.0) == 163.89                     # linear between, clamped below
    assert np.all(m.liquid_fraction(np.array([300.0, 900.0])) == 0.0) and m.cp_eff(500.0) == m.cp(500.0)


def test_enthalpy_is_the_exact_integral_of_the_piecewise_linear_cp():
    m = material.Material.from_drama_json()
    assert m.enthalpy(293.0) == 0.0
    assert m.enthalpy(313.0) == pytest.approx(20.0 * (877.5 + 885.5) / 2.0)
    assert m.enthalpy(303.0) == pytest.approx(10.0 * (877.5 + 881.5) / 2.0)
    assert m.enthalpy(1850.0) == pytest.approx(m.enthalpy(850.0) + 1000.0 * 1131.6)          # constant c_p above the table
    assert m.enthalpy(273.0) == pytest.approx(-20.0 * 877.5)                                 # constant c_p below it
    T = np.array([250.0, 293.0, 400.0, 733.0, 850.0, 1500.0, 3000.0])
    assert np.all(np.diff(m.enthalpy(T)) > 0)


def test_temperature_from_enthalpy_round_trips():
    m = material.Material.from_drama_json()
    T = np.array([250.0, 293.0, 310.5, 512.0, 733.0, 850.0, 1234.5, 5000.0])
    assert m.temperature_from_enthalpy(m.enthalpy(T)) == pytest.approx(T, rel=1e-12)
    assert float(m.temperature_from_enthalpy(0.0)) == 293.0


def test_constant_property_material():
    m = material.Material("const", 1000.0, 0.5, [200.0, 2000.0], [900.0, 900.0], [200.0, 2000.0], [10.0, 10.0])
    assert m.enthalpy(393.0) == pytest.approx(900.0 * 100.0) and m.k(1500.0) == 10.0
    assert m.temperature_from_enthalpy(900.0 * 1000.0) == pytest.approx(1293.0)
```

- [ ] **Step 4: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py -q`
Expected: ImportError / ModuleNotFoundError for `reentry_model.material` (collection error).

- [ ] **Step 5: Write `reentry_model/material.py`**

```python
"""Solid material for the thermal model: density, k(T), c_p(T), emissivity, and the specific enthalpy h(T) with a
latent-heat hook (zero in Step 2). Tables come from a DRAMA material JSON (the wrapper's `user_materials` format);
np.interp holds the end values outside the tabulated range. Spec section 6.4."""
import json
import os
from dataclasses import dataclass, field

import numpy as np

from . import DATA_DIR

DEFAULT_MATERIAL = os.path.join(DATA_DIR, "materials", "AA7075_nomelt.json")
T_REF = 293.0                      # K, zero of the enthalpy scale (first row of DRAMA's tables)


@dataclass
class Material:
    name: str
    rho: float                      # kg/m3
    emissivity: float
    T_cp: np.ndarray                # K, nodes of the c_p table
    cp_table: np.ndarray            # J/(kg K)
    T_k: np.ndarray                 # K, nodes of the k table
    k_table: np.ndarray             # W/(m K)
    latent_heat: float = 0.0        # J/kg; Step 3 sets it together with a melting range
    T_solidus: float = np.inf
    T_liquidus: float = np.inf
    _T_h: np.ndarray = field(init=False, repr=False)
    _h_nodes: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.T_cp, self.cp_table = np.asarray(self.T_cp, dtype=float), np.asarray(self.cp_table, dtype=float)
        self.T_k, self.k_table = np.asarray(self.T_k, dtype=float), np.asarray(self.k_table, dtype=float)
        # h(T) on the c_p nodes (plus T_REF): the trapezoid rule is exact for the piecewise-linear c_p
        T = np.union1d(self.T_cp, [T_REF])
        cp = np.interp(T, self.T_cp, self.cp_table)
        h = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(T))])
        self._T_h, self._h_nodes = T, h - np.interp(T_REF, T, h)

    @classmethod
    def from_drama_json(cls, path=None):
        with open(path or DEFAULT_MATERIAL) as fh:
            d = json.load(fh)
        cp = np.array(d["specificHeatCapacity"], dtype=float)
        k = np.array(d["heatConductivity"], dtype=float)
        return cls(d["name"], float(d["density"]), float(d["emissivity"][0][1]), cp[:, 0], cp[:, 1], k[:, 0], k[:, 1])

    def k(self, T):
        return np.interp(T, self.T_k, self.k_table)

    def cp(self, T):
        return np.interp(T, self.T_cp, self.cp_table)

    def liquid_fraction(self, T):
        """Melt fraction f_l(T): identically zero in Step 2 (mass and phase are constant); Step 3 replaces it."""
        return np.zeros_like(np.asarray(T, dtype=float))

    def cp_eff(self, T):
        """Effective heat capacity c_p + L_f df_l/dT; equals c_p while f_l is zero."""
        return self.cp(T)

    def enthalpy(self, T):
        """Specific enthalpy above T_REF [J/kg]: the exact integral of the piecewise-linear c_p (quadratic inside each
        table interval, linear beyond the table) plus L_f f_l(T)."""
        T = np.asarray(T, dtype=float)
        i = np.clip(np.searchsorted(self._T_h, T, side="right") - 1, 0, len(self._T_h) - 2)
        T0, T1 = self._T_h[i], self._T_h[i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        x = np.clip(T, T0, T1) - T0
        h = self._h_nodes[i] + cp0 * x + 0.5 * (cp1 - cp0) / (T1 - T0) * x * x
        h = np.where(T < self._T_h[0], self._h_nodes[0] + self.cp_table[0] * (T - self._T_h[0]), h)
        h = np.where(T > self._T_h[-1], self._h_nodes[-1] + self.cp_table[-1] * (T - self._T_h[-1]), h)
        return h + self.latent_heat * self.liquid_fraction(T)

    def temperature_from_enthalpy(self, h):
        """Inverse of enthalpy() (monotonic): the energy-equivalent temperature of a body holding h per kg."""
        h = np.asarray(h, dtype=float)
        i = np.clip(np.searchsorted(self._h_nodes, h, side="right") - 1, 0, len(self._T_h) - 2)
        T0, T1 = self._T_h[i], self._T_h[i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        a = (cp1 - cp0) / (T1 - T0)
        dh = np.clip(h, self._h_nodes[i], self._h_nodes[i + 1]) - self._h_nodes[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(np.abs(a) > 1e-12, (np.sqrt(cp0 * cp0 + 2.0 * a * dh) - cp0) / a, dh / cp0)
        T = T0 + x
        T = np.where(h < self._h_nodes[0], self._T_h[0] + (h - self._h_nodes[0]) / self.cp_table[0], T)
        T = np.where(h > self._h_nodes[-1], self._T_h[-1] + (h - self._h_nodes[-1]) / self.cp_table[-1], T)
        return T
```

- [ ] **Step 6: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py -q`
Expected: `5 passed`.

- [ ] **Step 7: Commit**

```bash
git add requirements-step2.txt reentry_model/data/materials/AA7075_nomelt.json reentry_model/material.py tests/test_reentry_model_material.py
git commit -m "reentry_model: Step 2 environment pins and the AA7075 material model with exact enthalpy

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Sphere mesh generation and surface geometry

**Files:**
- Create: `reentry_model/mesh.py`
- Modify: `tests/conftest.py` (append two session fixtures)
- Test: `tests/test_reentry_model_mesh.py`

**Interfaces:**
- Consumes: gmsh (Python API), meshio.
- Produces: `VolumeMesh(points (n,3) mutable, tets (ne,4), params)` with `n_nodes`, `n_elements`, `element_volumes()`, `volume()`, `surface() -> SurfaceMesh`, `boundary_nodes()`, `centre_node()`; `SurfaceMesh(faces (nf,3), centroids, normals (outward unit), areas)` with `n_patches`, `area`, `angles_to(v_hat) -> theta (nf,)`, `facet_mean(nodal)`, `patch_toward(direction)`; `boundary_faces(tets)`; `mesh_file_name(R, hs, hc)`; `generate_sphere_mesh(radius, h_surface, h_core, mesh_dir) -> path` (cached); `load_mesh(path)`; `sphere_mesh(radius, h_surface, h_core, mesh_dir) -> VolumeMesh`; constants `DEFAULT_H_SURFACE = 2e-3`, `DEFAULT_H_CORE = 8e-3`, `DEFAULT_MESH_DIR = <repo>/reentry_model_output/meshes`.
- Conventions later tasks rely on: `theta = arccos(n · v_hat)` with `v_hat` the direction the body moves in, so the stagnation patch has θ = 0; boundary faces are derived from the tets (any tet mesh works); a node sits exactly at the origin.

- [ ] **Step 1: Add the session fixtures**

Append to `tests/conftest.py`:

```python
@pytest.fixture(scope="session")
def coarse_sphere_mesh(tmp_path_factory):
    """100 mm sphere, 4 mm surface / 20 mm core (~2.9 k nodes): the mesh for the fast solver tests."""
    from reentry_model import mesh
    return mesh.sphere_mesh(0.05, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")))


@pytest.fixture(scope="session")
def uniform_test_mesh(tmp_path_factory):
    """100 mm sphere, uniform 4 mm elements (~6.7 k nodes): resolves the centre for the Carslaw-Jaeger test."""
    from reentry_model import mesh
    return mesh.sphere_mesh(0.05, 4e-3, 4e-3, str(tmp_path_factory.mktemp("meshes")))
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_reentry_model_mesh.py`:

```python
"""mesh.py: gmsh sphere generation and cache, boundary faces, outward normals, patch angles."""
import math
import os

import numpy as np
import pytest

from reentry_model import mesh

R = 0.05


def test_coarse_sphere_geometry(coarse_sphere_mesh):
    m = coarse_sphere_mesh
    assert 2000 < m.n_nodes < 5000 and m.n_elements > 4 * m.n_nodes // 2
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * R ** 3, rel=3e-3)
    s = m.surface()
    assert s.area == pytest.approx(4.0 * math.pi * R ** 2, rel=2e-3)
    assert np.allclose(np.linalg.norm(s.centroids, axis=1), R, atol=3e-4)      # centroids sag by ~h^2/(6R) below the sphere
    assert np.all(np.einsum("ij,ij->i", s.normals, s.centroids) > 0.95 * np.linalg.norm(s.centroids, axis=1))   # outward, within 18 deg of radial
    assert np.allclose(np.linalg.norm(s.normals, axis=1), 1.0)
    assert np.linalg.norm(m.points[m.centre_node()]) < 1e-9
    assert m.boundary_nodes().size == np.unique(s.faces).size and m.params["radius_m"] == R


def test_angles_and_patch_lookup(coarse_sphere_mesh):
    s = coarse_sphere_mesh.surface()
    theta = s.angles_to([1.0, 0.0, 0.0])
    assert theta.min() < math.radians(5.0) and theta.max() > math.radians(175.0)
    windward = theta < math.pi / 2
    assert s.areas[windward].sum() == pytest.approx(0.5 * s.area, rel=2e-2)
    assert s.normals[s.patch_toward([1.0, 0.0, 0.0])][0] > 0.99
    assert s.normals[s.patch_toward([-1.0, 0.0, 0.0])][0] < -0.99
    x = coarse_sphere_mesh.points[:, 0]
    assert np.allclose(s.facet_mean(x), s.centroids[:, 0])


def test_boundary_faces_are_single_use_faces():
    tets = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    faces, opposite = mesh.boundary_faces(tets)
    assert faces.shape == (6, 3) and not any(set(f) == {1, 2, 3} for f in faces)
    assert opposite.shape == (6,)


def test_generation_is_cached(tmp_path):
    p1 = mesh.generate_sphere_mesh(0.01, 2e-3, 4e-3, str(tmp_path))
    mtime = os.path.getmtime(p1)
    p2 = mesh.generate_sphere_mesh(0.01, 2e-3, 4e-3, str(tmp_path))
    assert p1 == p2 and os.path.getmtime(p2) == mtime
    assert os.path.basename(p1) == "sphere_R10.000mm_hs2.000mm_hc4.000mm.msh"
    m = mesh.load_mesh(p1)
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * 0.01 ** 3, rel=3e-2)        # crude 2 mm elements on a 10 mm sphere


def test_points_are_mutable_and_surface_follows(coarse_sphere_mesh):
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets)
    area0 = m.surface().area
    m.points *= 0.5
    assert m.surface().area == pytest.approx(0.25 * area0) and m.volume() == pytest.approx(0.125 * coarse_sphere_mesh.volume())
```

- [ ] **Step 3: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: collection error (`reentry_model.mesh` missing).

- [ ] **Step 4: Write `reentry_model/mesh.py`**

```python
"""Tetrahedral sphere meshes (gmsh, graded) and their surface geometry. Spec section 5.

A VolumeMesh is node coordinates (a mutable array, so Step 3 can recede the surface), tetrahedra, and the boundary
triangles derived from the tetrahedra themselves (faces used by exactly one element) with outward unit normals,
centroids and areas -- so any gmsh volume mesh loads, tagged or not. The generator makes a sphere whose element
size grows linearly from h_surface at the surface to h_core at the centre (gmsh Distance/Threshold field) with a
node embedded at the centre, and caches the .msh by (R, h_surface, h_core).
"""
import os
from dataclasses import dataclass, field

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MESH_DIR = os.path.join(REPO_ROOT, "reentry_model_output", "meshes")
DEFAULT_H_SURFACE = 2.0e-3          # m; 2 mm / 8 mm gives 18.9 k nodes on the 100 mm sphere (measured 2026-09-18)
DEFAULT_H_CORE = 8.0e-3             # m
_FACE_OF_VERTEX = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])     # face opposite each tet vertex


@dataclass
class SurfaceMesh:
    faces: np.ndarray            # (nf, 3) node ids of the boundary triangles ("patches")
    centroids: np.ndarray        # (nf, 3) m
    normals: np.ndarray          # (nf, 3) outward unit normals
    areas: np.ndarray            # (nf,) m2

    @property
    def n_patches(self):
        return len(self.faces)

    @property
    def area(self):
        return float(self.areas.sum())

    def angles_to(self, v_hat):
        """theta per patch [rad]: the angle between the outward normal and v_hat, the direction the body moves in.
        The stagnation patch has its normal along v_hat (theta = 0); theta > pi/2 is leeward."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        return np.arccos(np.clip(self.normals @ v, -1.0, 1.0))

    def facet_mean(self, nodal):
        """Mean of a nodal field over each patch's three nodes."""
        return np.asarray(nodal)[self.faces].mean(axis=1)

    def patch_toward(self, direction):
        """Index of the patch whose outward normal is closest to `direction`."""
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))


@dataclass
class VolumeMesh:
    points: np.ndarray                                   # (n, 3) m, mutable
    tets: np.ndarray                                     # (ne, 4) node ids
    params: dict = field(default_factory=dict)           # generator parameters / source file
    _faces: np.ndarray = field(init=False, repr=False, default=None)
    _opposite: np.ndarray = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.points = np.array(self.points, dtype=float)
        self.tets = np.asarray(self.tets, dtype=np.int64)
        self._faces, self._opposite = boundary_faces(self.tets)

    @property
    def n_nodes(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.tets)

    def element_volumes(self):
        x = self.points[self.tets]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        return np.abs(np.linalg.det(J)) / 6.0

    def volume(self):
        return float(self.element_volumes().sum())

    def surface(self):
        """The boundary patches with their current geometry (recomputed from `points` on every call)."""
        a, b, c = (self.points[self._faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        areas = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * areas)[:, None]
        centroids = (a + b + c) / 3.0
        inward = np.einsum("ij,ij->i", n, centroids - self.points[self._opposite]) < 0.0
        n[inward] *= -1.0
        return SurfaceMesh(self._faces, centroids, n, areas)

    def boundary_nodes(self):
        return np.unique(self._faces)

    def centre_node(self):
        """Node nearest the origin (the embedded centre node of a generated sphere)."""
        return int(np.argmin(np.linalg.norm(self.points, axis=1)))


def boundary_faces(tets):
    """Faces used by exactly one tetrahedron, as (nf, 3) node ids, and the opposite vertex of that tetrahedron."""
    faces = tets[:, _FACE_OF_VERTEX].reshape(-1, 3)
    opposite = np.repeat(tets, 4, axis=0)[np.arange(4 * len(tets)), np.tile(np.arange(4), len(tets))]
    _, first, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    keep = first[counts == 1]
    return faces[keep], opposite[keep]


def mesh_file_name(radius, h_surface, h_core):
    return "sphere_R{:.3f}mm_hs{:.3f}mm_hc{:.3f}mm.msh".format(radius * 1e3, h_surface * 1e3, h_core * 1e3)


def generate_sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR):
    """gmsh sphere of `radius` graded from h_surface (surface) to h_core (centre); returns the cached .msh path."""
    os.makedirs(mesh_dir, exist_ok=True)
    path = os.path.join(mesh_dir, mesh_file_name(radius, h_surface, h_core))
    if os.path.isfile(path):
        return path
    import gmsh                                   # imported here: only mesh generation needs gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("sphere")
        vol = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, radius)
        centre = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h_core)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(0, [centre], 3, vol)
        surfaces = [s[1] for s in gmsh.model.getBoundary([(3, vol)], oriented=False)]
        gmsh.model.addPhysicalGroup(3, [vol], tag=1, name="body")
        gmsh.model.addPhysicalGroup(2, surfaces, tag=2, name="surface")
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", surfaces)
        gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
        thr = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thr, "InField", dist)
        gmsh.model.mesh.field.setNumber(thr, "SizeMin", h_surface)
        gmsh.model.mesh.field.setNumber(thr, "SizeMax", h_core)
        gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", radius)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.model.mesh.generate(3)
        gmsh.write(path)
    finally:
        gmsh.finalize()
    return path


def load_mesh(path):
    """Any gmsh/meshio volume mesh with tetrahedra (other cell types are ignored)."""
    import meshio
    m = meshio.read(path)
    tets = [c.data for c in m.cells if c.type == "tetra"]
    if not tets:
        raise ValueError("no tetrahedra in {}".format(path))
    return VolumeMesh(np.asarray(m.points, dtype=float), np.vstack(tets), {"path": os.path.abspath(path)})


def sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR):
    mesh = load_mesh(generate_sphere_mesh(radius, h_surface, h_core, mesh_dir))
    mesh.params.update({"radius_m": radius, "h_surface_m": h_surface, "h_core_m": h_core})
    return mesh
```

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: `5 passed` in a few seconds (gmsh generates the coarse mesh once per session). The coarse mesh has 2 934 nodes / 10 191 tets; its polyhedral area and volume are 0.14 % / 0.25 % below the sphere's — the analytic solver tests compare against the discrete sphere's own area and volume for that reason.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/mesh.py tests/conftest.py tests/test_reentry_model_mesh.py
git commit -m "reentry_model: gmsh graded sphere mesh with cached generation and patch geometry

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: Thermal solver protocol and the scikit-fem backend

**Files:**
- Create: `reentry_model/thermal/__init__.py`
- Create: `reentry_model/thermal/skfem_backend.py`
- Create: `reentry_model/thermal/fenicsx_backend.py` (stub; Task 11 replaces it)
- Test: `tests/test_reentry_model_thermal.py`

**Interfaces:**
- Consumes: `mesh.VolumeMesh` (`points`, `tets`, `surface()`), `material.Material` (`rho`, `k`, `cp_eff`, `enthalpy`), fixtures `coarse_sphere_mesh`, `uniform_test_mesh`.
- Produces: `thermal.SIGMA_SB`, `thermal.SOLVER_NAMES`, `thermal.MissingBackend`, `thermal.StepResult(T, Q_conv, Q_rad, iterations)`, protocol `ThermalSolver` with `setup(mesh, material, emissivity)`, `set_temperature(T)`, `step(dt, q_conv, T_amb, dirichlet=None) -> StepResult`, `temperature()`, `energy()`, `radiated_power(T_amb)`; `thermal.thermal_solver(name, **options)`; `SkfemThermalSolver(linear_solver="amg"|"direct", lumped_mass=False, newton_tol=1e-6, max_iterations=30, amg_rebuild_every=30, cg_tol=1e-10)` with the extra attributes `faces`, `areas`, `vol`, `tets`, `points`, `T`, `facet_temperature()`, `facet_load(q)`, `operators(T, T_old=None)`, `reference_operators(T)`, `last_cg_iterations`.
- Numerics (spec §7 as amended): backward Euler; per Newton iterate K with k(T̄ₑ), M with the secant heat capacity, radiation from the facet-mean temperature with its exact Jacobian (A_f/9 all-ones block); the system is SPD; CG (pyamg V-cycle preconditioner, hierarchy rebuilt every `amg_rebuild_every` solves) or SuperLU; predictor T_old + (T_old − T_prev).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_thermal.py`:

```python
"""thermal/: the scikit-fem conduction solver against scikit-fem's own assembly and the analytic cases of spec
section 10 (lumped limit, radiative cooling, Carslaw-Jaeger, per-step energy balance)."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from reentry_model import material, thermal
from reentry_model.thermal import SIGMA_SB, skfem_backend

R = 0.05
RHO, CP, K0, EPS = 2813.0, 900.0, 160.0, 0.4


def constant_material(k=K0, cp=CP, rho=RHO, eps=EPS):
    return material.Material("const", rho, eps, [100.0, 20000.0], [cp, cp], [100.0, 20000.0], [k, k])


def solver(mesh, mat, **options):
    s = thermal.thermal_solver("skfem", **options)
    s.setup(mesh, mat, mat.emissivity)
    return s


def lumped_reference(q, A, m, cp, eps, T0, t_end, T_amb=0.0):
    """Exact lumped body: m c_p dT/dt = q A - eps sigma A (T^4 - T_amb^4)."""
    sol = solve_ivp(lambda t, y: [(q * A - eps * SIGMA_SB * A * (y[0] ** 4 - T_amb ** 4)) / (m * cp)], (0.0, t_end), [T0],
                    rtol=1e-12, atol=1e-10)
    return float(sol.y[0, -1])


def test_factory_and_options():
    assert isinstance(thermal.thermal_solver("skfem"), skfem_backend.SkfemThermalSolver)
    with pytest.raises(ValueError):
        thermal.thermal_solver("nope")
    with pytest.raises(ValueError):
        thermal.thermal_solver("skfem", linear_solver="magic")


def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets].mean(axis=1)) * s.vol).sum(), rel=1e-12)
    u = 2.0 * s.points[:, 0] + 3.0 * s.points[:, 1] - s.points[:, 2]           # linear field: u^T K u = k |grad u|^2 V
    K1, _ = solver(coarse_sphere_mesh, constant_material(k=1.0)).operators(T)
    assert u @ (K1 @ u) == pytest.approx(14.0 * s.vol.sum(), rel=1e-10)


def test_lumped_limit_matches_the_lumped_ode(coarse_sphere_mesh):
    """k x 1e4 (Bi -> 0), uniform flux 1e6 W/m2 with radiation, dt 0.5 s for 50 s: T_eq vs the lumped ODE. The lumped
    body uses the discrete sphere's own area and volume. (k x 1e6 makes CG fail on conditioning: keep 1e4.)"""
    s = solver(coarse_sphere_mesh, constant_material(k=K0 * 1e4), newton_tol=1e-8)
    s.set_temperature(300.0)
    A, m = s.areas.sum(), RHO * s.vol.sum()
    q = np.full(len(s.faces), 1e6)
    for _ in range(100):
        res = s.step(0.5, q, 0.0)
    assert res.Q_conv == pytest.approx(1e6 * A, rel=1e-12) and 1 <= res.iterations <= 6
    T_eq = s.energy() / (m * CP) + material.T_REF
    assert T_eq == pytest.approx(lumped_reference(1e6, A, m, CP, EPS, 300.0, 50.0), rel=1e-3)      # spec: 0.1 %
    assert s.T.max() - s.T.min() < 0.1                                                              # isothermal


def test_radiative_cooling_of_an_isothermal_sphere(coarse_sphere_mesh):
    """No flux, k x 1e4, T0 1500 K, dt 0.5 s for 100 s: T(t) = [T0^-3 + 3 eps sigma A t / (m c_p)]^(-1/3), 0.1 %."""
    s = solver(coarse_sphere_mesh, constant_material(k=K0 * 1e4), newton_tol=1e-8)
    s.set_temperature(1500.0)
    A, m = s.areas.sum(), RHO * s.vol.sum()
    for _ in range(200):
        res = s.step(0.5, np.zeros(len(s.faces)), 0.0)
    exact = (1500.0 ** -3 + 3.0 * EPS * SIGMA_SB * A * 100.0 / (m * CP)) ** (-1.0 / 3.0)
    assert s.energy() / (m * CP) + material.T_REF == pytest.approx(exact, rel=1e-3)
    assert res.Q_rad == pytest.approx(EPS * SIGMA_SB * A * exact ** 4, rel=2e-3) and res.Q_conv == 0.0


def test_carslaw_jaeger_step_change_at_the_surface(uniform_test_mesh):
    """Surface stepped from 300 to 800 K and held (Dirichlet), constant properties, dt 0.05 s: centre temperature vs the
    series 2 sum (-1)^(n+1) exp(-n^2 pi^2 alpha t / R^2) within 1 % of the 500 K rise, volume mean within 0.5 %."""
    s = solver(uniform_test_mesh, constant_material(eps=0.0), newton_tol=1e-10)
    s.set_temperature(300.0)
    alpha = K0 / (RHO * CP)
    boundary = uniform_test_mesh.boundary_nodes()
    centre, V = uniform_test_mesh.centre_node(), s.vol.sum()

    def series(t, coefficient):
        return sum(coefficient(n) * math.exp(-n * n * math.pi ** 2 * alpha * t / R ** 2) for n in range(1, 80))

    for i in range(1, 81):
        s.step(0.05, np.zeros(len(s.faces)), 0.0, dirichlet=(boundary, np.full(boundary.size, 800.0)))
        t = 0.05 * i
        if i in (40, 80):
            T_centre_exact = 800.0 - 500.0 * 2.0 * series(t, lambda n: (-1) ** (n + 1))
            T_mean_exact = 800.0 - 500.0 * (6.0 / math.pi ** 2) * series(t, lambda n: 1.0 / n ** 2)
            T_mean = float((s.vol * s.T[s.tets].mean(axis=1)).sum() / V)
            assert abs(s.T[centre] - T_centre_exact) < 0.01 * 500.0, (t, s.T[centre], T_centre_exact)
            assert abs(T_mean - T_mean_exact) < 0.005 * 500.0, (t, T_mean, T_mean_exact)
    assert np.all(s.T[boundary] == 800.0)


@pytest.mark.parametrize("lumped", [False, True])
def test_energy_balance_every_step(coarse_sphere_mesh, lumped):
    """Constant c_p, cos(theta) flux 2e6 W/m2 on the windward side, radiation on: dE = (Q_conv - Q_rad) dt to 1e-6."""
    s = solver(coarse_sphere_mesh, constant_material(), newton_tol=1e-12, lumped_mass=lumped)
    s.set_temperature(300.0)
    theta = np.arccos(np.clip(coarse_sphere_mesh.surface().normals[:, 0], -1.0, 1.0))
    q = 2e6 * np.where(theta < math.pi / 2, np.cos(theta), 0.0)
    for _ in range(20):
        E0 = s.energy()
        res = s.step(0.5, q, 0.0)
        assert abs((s.energy() - E0) - (res.Q_conv - res.Q_rad) * 0.5) < 1e-6 * abs((res.Q_conv - res.Q_rad) * 0.5)
    assert 700.0 < s.T.max() < 800.0 and 305.0 < s.T.min() < 320.0 and res.Q_rad > 0.0


def test_variable_properties_energy_balance_is_exact_too(coarse_sphere_mesh):
    """With the AA7075 tables the secant heat capacity keeps dE = (Q_conv - Q_rad) dt to the Newton tolerance."""
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json(), newton_tol=1e-12)
    s.set_temperature(300.0)
    q = np.full(len(s.faces), 5e5)
    for _ in range(40):
        E0 = s.energy()
        res = s.step(0.5, q, 0.0)
        assert abs((s.energy() - E0) - (res.Q_conv - res.Q_rad) * 0.5) < 1e-6 * abs((res.Q_conv - res.Q_rad) * 0.5)
    assert s.T.max() > 400.0


def test_direct_and_amg_agree(coarse_sphere_mesh):
    q = np.full(len(coarse_sphere_mesh.surface().faces), 1e6)
    results = []
    for linear in ("direct", "amg"):
        s = solver(coarse_sphere_mesh, material.Material.from_drama_json(), linear_solver=linear, newton_tol=1e-10)
        s.set_temperature(300.0)
        for _ in range(5):
            s.step(0.5, q, 0.0)
        results.append(s.T)
    assert np.abs(results[0] - results[1]).max() < 1e-6


def test_step_result_and_temperature_setters(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, constant_material())
    s.set_temperature(np.linspace(300.0, 400.0, coarse_sphere_mesh.n_nodes))
    assert s.temperature()[0] == 300.0 and s.temperature() is not s.T
    res = s.step(0.5, np.zeros(len(s.faces)), 200.0)
    assert isinstance(res, thermal.StepResult) and res.T.shape == (coarse_sphere_mesh.n_nodes,)
    assert res.Q_rad == pytest.approx(s.radiated_power(200.0)) and res.Q_conv == 0.0
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`
Expected: collection error (`reentry_model.thermal` missing).

- [ ] **Step 3: Write `reentry_model/thermal/__init__.py`**

```python
"""Finite-element conduction solvers behind one protocol (spec section 8); `thermal_solver(name)` picks the backend.

Both backends solve rho c_p(T) dT/dt = div(k(T) grad T) with -k grad T . n = eps sigma (T^4 - T_amb^4) - q_conv on the
boundary, backward Euler in time, P1 tetrahedra in space, and report the same quantities: nodal temperatures,
stored enthalpy, radiated power. `dirichlet=(nodes, values)` in `step` is for the analytic tests only."""
from dataclasses import dataclass
from typing import Protocol

import numpy as np

SIGMA_SB = 5.670374419e-8         # Stefan-Boltzmann [W/(m2 K4)]
SOLVER_NAMES = ("skfem", "fenicsx")


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


def thermal_solver(name, **options):
    if name == "skfem":
        from .skfem_backend import SkfemThermalSolver
        return SkfemThermalSolver(**options)
    if name == "fenicsx":
        from .fenicsx_backend import FenicsxThermalSolver      # raises MissingBackend without dolfinx
        return FenicsxThermalSolver(**options)
    raise ValueError("thermal solver must be one of {}, got {!r}".format(SOLVER_NAMES, name))
```

- [ ] **Step 4: Write the FEniCSx stub `reentry_model/thermal/fenicsx_backend.py`** (Task 11 replaces the whole file)

```python
"""FEniCSx (dolfinx) backend: Task 11 fills in the solver; until then the constructor only reports the missing library."""
from . import MissingBackend


class FenicsxThermalSolver:
    def __init__(self, **options):
        try:
            import dolfinx  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        raise NotImplementedError("FenicsxThermalSolver is implemented in Task 11")
```

- [ ] **Step 5: Write `reentry_model/thermal/skfem_backend.py`**

```python
"""scikit-fem backend: P1 tetrahedra, backward Euler, Picard on k(T)/c_p(T) with Newton on the radiation term.

Assembly. The P1 element stiffness and mass matrices depend on T only through one coefficient per element
(k(T_e) and rho c_p(T_e) at the element-mean temperature), so the unit-coefficient element matrices are computed
once and rescaled into a fixed CSR pattern on every iteration (~20 ms for 50 k tets). scikit-fem's generic `asm`
on the same MeshTet/ElementTetP1 costs 0.2 s per iteration at 12 k nodes (measured 2026-09-18) and is kept as the
reference in `reference_operators` for the conformance test.

Radiation. The boundary functional uses the facet-mean temperature: eps sigma (T_f^4 - T_amb^4) A_f/3 to each of
the facet's three nodes, with its exact Jacobian 4 eps sigma T_f^3 A_f/9 on every node pair of the facet. Each
Newton iterate solves the symmetric positive definite system
    (M/dt + K + B_k) T = M T_old/dt + F_conv - F_rad(T_k) + B_k T_k
with M, K frozen at the previous iterate; because the rows of K sum to zero the discrete energy balance
1^T M (T - T_old) = dt (Q_conv - Q_rad) holds to the Newton tolerance.

Linear solver. `direct`: SciPy SuperLU; `amg` (default): CG preconditioned by a pyamg smoothed-aggregation
hierarchy rebuilt every `amg_rebuild_every` solves (the operator changes slowly). Measured on the 100 mm
2 mm/8 mm mesh (18.9 k nodes): AMG 0.03 s per solve vs SuperLU 0.6 s."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import SIGMA_SB, StepResult


class _Pattern:
    """Fixed CSR sparsity for repeated (rows, cols) entries; `assemble` sums values into it with one bincount."""

    def __init__(self, rows, cols, n):
        key = rows.astype(np.int64) * n + cols.astype(np.int64)
        unique, self._map = np.unique(key, return_inverse=True)
        self._indptr = np.searchsorted(unique // n, np.arange(n + 1))
        self._indices = (unique % n).astype(np.int32)
        self.n, self.nnz = n, unique.size

    def assemble(self, values):
        data = np.bincount(self._map, weights=np.asarray(values).ravel(), minlength=self.nnz)
        return sp.csr_matrix((data, self._indices, self._indptr), shape=(self.n, self.n))


def element_matrices(points, tets):
    """Per-element volumes, unit stiffness (grad phi_a . grad phi_b V_e) and consistent mass (V_e/20 (1 + delta_ab))."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    grad_ref = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    grads = np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), grad_ref)
    Ke = np.einsum("eai,ebi->eab", grads, grads) * vol[:, None, None]
    Me = (np.ones((4, 4)) + np.eye(4))[None] * (vol / 20.0)[:, None, None]
    return vol, Ke, Me


class SkfemThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=False, newton_tol=1e-6, max_iterations=30,
                 amg_rebuild_every=30, cg_tol=1e-10):
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        self.linear_solver, self.lumped_mass = linear_solver, lumped_mass
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.amg_rebuild_every, self.cg_tol = amg_rebuild_every, cg_tol
        self._ml, self._solves, self.last_cg_iterations = None, 0, 0

    def setup(self, mesh, material, emissivity):
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        self.points, self.tets = mesh.points, mesh.tets
        surface = mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        self.vol, self.Ke, self.Me = element_matrices(self.points, self.tets)
        n = len(self.points)
        self.pattern = _Pattern(np.repeat(self.tets, 4, axis=1).ravel(), np.tile(self.tets, (1, 4)).ravel(), n)
        self.facet_pattern = _Pattern(np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel(), n)
        self.Bf = np.ones((3, 3))[None] * (self.areas / 9.0)[:, None, None]
        self.T, self._T_prev = np.full(n, 300.0), None

    def set_temperature(self, T):
        self.T = np.full(len(self.points), float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    def facet_temperature(self, T=None):
        return (self.T if T is None else T)[self.faces].mean(axis=1)

    def facet_load(self, q):
        """Nodal load vector of a per-facet flux q [W/m2]: A_f/3 to each of the facet's nodes."""
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=len(self.points))

    def operators(self, T, T_old=None):
        """Stiffness K with k at the element-mean temperature, and mass M with the secant heat capacity
        [h(T_e) - h(T_old,e)] / (T_e - T_old,e) (c_p_eff where the element has not moved), so that at convergence
        1^T M (T - T_old) equals the enthalpy increment exactly, whatever h(T) is (Step 3's latent heat included)."""
        Te = T[self.tets].mean(axis=1)
        K = self.pattern.assemble(self.Ke * self.material.k(Te)[:, None, None])
        c = self.material.cp_eff(Te)
        if T_old is not None:
            Te_old = T_old[self.tets].mean(axis=1)
            dT = Te - Te_old
            moved = np.abs(dT) > 1e-9
            c = np.where(moved, (self.material.enthalpy(Te) - self.material.enthalpy(Te_old)) / np.where(moved, dT, 1.0), c)
        M = self.pattern.assemble(self.Me * (self.material.rho * c)[:, None, None])
        if self.lumped_mass:
            M = sp.diags(np.asarray(M.sum(axis=1)).ravel(), format="csr")
        return K, M

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def energy(self):
        """Stored enthalpy above T_REF: sum_e rho V_e h(T_e); equals 1^T M T for constant c_p (consistent or lumped mass)."""
        Te = self.T[self.tets].mean(axis=1)
        return float((self.material.rho * self.vol * self.material.enthalpy(Te)).sum())

    def step(self, dt, q_conv, T_amb, dirichlet=None):
        T_old = self.T
        # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        T_new, iteration = T_k, 0
        for iteration in range(1, self.max_iterations + 1):
            K, M = self.operators(T_k, T_old)
            Tf = self.facet_temperature(T_k)
            B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            A = (M / dt + K + B).tocsr()
            b = M @ T_old / dt + F_conv - F_rad + B @ T_k
            T_new = self._solve_with_dirichlet(A, b, T_k, dirichlet)
            if np.linalg.norm(T_new - T_k) <= self.newton_tol * np.linalg.norm(T_new):
                break
            T_k = T_new
        self._T_prev, self.T = T_old, T_new
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration)

    def _solve_with_dirichlet(self, A, b, x0, dirichlet):
        if dirichlet is None:
            return self._solve(A, b, x0)
        nodes, values = dirichlet
        free = np.ones(A.shape[0], dtype=bool)
        free[nodes] = False
        x = np.zeros(A.shape[0])
        x[nodes] = values
        x[free] = self._solve(A[free][:, free], b[free] - A[free][:, ~free] @ x[~free], x0[free], fresh=True)
        return x

    def _solve(self, A, b, x0, fresh=False):
        if self.linear_solver == "direct":
            return spla.spsolve(A.tocsc(), b)
        if fresh or self._ml is None or self._solves % self.amg_rebuild_every == 0:
            import pyamg
            self._ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
        self._solves += 1
        counter = []
        x, info = spla.cg(A, b, x0=x0, rtol=self.cg_tol, maxiter=500, M=self._ml.aspreconditioner(cycle="V"),
                          callback=lambda _: counter.append(1))
        self.last_cg_iterations = len(counter)
        if info != 0:
            raise RuntimeError("CG did not converge (info {})".format(info))
        return x

    def reference_operators(self, T):
        """K and M assembled by scikit-fem with the same element-mean coefficients (conformance test only)."""
        from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP1, MeshTet, asm
        from skfem.helpers import dot, grad
        basis = Basis(MeshTet(self.points.T.copy(), self.tets.T.copy()), ElementTetP1())
        basis0 = basis.with_element(ElementTetP0())
        Te = T[self.tets].mean(axis=1)

        @BilinearForm
        def stiffness(u, v, w):
            return w.k * dot(grad(u), grad(v))

        @BilinearForm
        def mass(u, v, w):
            return w.c * u * v

        K = asm(stiffness, basis, k=basis0.interpolate(self.material.k(Te)))
        M = asm(mass, basis, c=basis0.interpolate(self.material.rho * self.material.cp_eff(Te)))
        return K, M
```

- [ ] **Step 6: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`
Expected: `10 passed` in ~30 s (Carslaw–Jaeger 7 s, radiative cooling 5 s, lumped limit 4 s). Measured accuracies on 2026-09-18: lumped limit 8e-5 vs the exact ODE (4e-7 vs backward Euler of the same ODE), radiative cooling 1.7e-4, Carslaw–Jaeger centre 0.5 % / mean 0.2 % of the rise, energy balance 1e-10 (constant c_p) and 1e-9 (AA7075 tables, secant c_p).

- [ ] **Step 7: Commit**

```bash
git add reentry_model/thermal/__init__.py reentry_model/thermal/skfem_backend.py reentry_model/thermal/fenicsx_backend.py tests/test_reentry_model_thermal.py
git commit -m "reentry_model: thermal solver protocol and scikit-fem P1 backward-Euler backend with exact energy balance

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Equilibrium-air gas model and stagnation correlations

**Files:**
- Create: `reentry_model/gas.py`
- Test: `tests/test_reentry_model_gas.py`

**Interfaces:**
- Consumes: Cantera (`airNASA9.yaml`, imported lazily inside `EquilibriumAir.__init__`), scipy `brentq`.
- Produces: `GasState(p, T, rho, h, mu, h_D, X)`; `blottner_viscosity(species, T)`; `wilke_viscosity(T, mole_fractions: dict, molar_masses: dict)`; `EquilibriumAir(mechanism="airNASA9.yaml")` with `air_enthalpy(T)` (frozen-air table, vectorised), `freestream(rho, T)`, `stagnation(rho_inf, T_inf, V)`, `wall(T_w, p)`, `molar_masses`; `dkr(rho, V, radius)`, `sutton_graves(rho, V, radius)`, `fay_riddell(stag, wall, p_inf, radius, catalycity=1.0, Pr=0.71, Le=1.4)`; constants `AIR`, `T_FLOOR`, `CP_AIR`, `T_EQUILIBRATE`, `PRANDTL`, `LEWIS`, `BLOTTNER`.
- Enthalpies are on Cantera's basis (elements at 298.15 K): h(air, 298 K) = 0, h_∞(203.6 K) = −9.55e4 J/kg, h_s = h_∞ + V²/2.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_gas.py`:

```python
"""gas.py: the stagnation correlations, Blottner/Wilke viscosity and the Cantera equilibrium stagnation state at the
100 mm reference start (US76 77.5 km: rho 2.727e-5 kg/m3, T 203.6 K, 7.5 km/s, R 0.05 m)."""
import warnings

import pytest

from reentry_model import gas

RHO, T_INF, V, R = 2.727e-5, 203.6, 7500.0, 0.05


def test_stagnation_correlations_at_the_start_state():
    assert gas.dkr(RHO, V, R) == pytest.approx(1.957e6, rel=2e-3)
    assert gas.sutton_graves(RHO, V, R) == pytest.approx(1.716e6, rel=2e-3)
    assert gas.sutton_graves(RHO, V, R) / gas.dkr(RHO, V, R) == pytest.approx(0.877, abs=2e-3)
    assert gas.blottner_viscosity("N2", 300.0) == pytest.approx(1.786e-5, rel=2e-3)
    assert gas.blottner_viscosity("O", 300.0) == pytest.approx(2.049e-5, rel=2e-3)
    assert gas.wilke_viscosity(300.0, {"N2": 1.0}, {"N2": 28.014}) == gas.blottner_viscosity("N2", 300.0)


@pytest.fixture(scope="module")
def air():
    pytest.importorskip("cantera")
    return gas.EquilibriumAir()


def test_equilibrium_stagnation_state(air):
    s = air.stagnation(RHO, T_INF, V)
    assert s.T == pytest.approx(5820.0, rel=1e-2) and s.p == pytest.approx(1495.0, rel=1e-2)
    assert s.h == pytest.approx(float(air.air_enthalpy(T_INF)) + 0.5 * V * V, rel=1e-6)
    assert s.h_D / s.h == pytest.approx(0.718, abs=0.02) and s.X["N"] > 0.5 and s.mu == pytest.approx(1.69e-4, rel=3e-2)
    w = air.wall(300.0, s.p)
    assert w.h == pytest.approx(float(air.air_enthalpy(300.0)), rel=1e-3) and w.mu == pytest.approx(1.955e-5, rel=1e-2)
    fs = air.freestream(RHO, T_INF)
    assert fs.p == pytest.approx(1.60, rel=1e-2)
    assert float(air.air_enthalpy(150.0)) == float(air.air_enthalpy(200.0))                    # table floor
    with warnings.catch_warnings():
        warnings.simplefilter("error")                                                         # subsonic: frozen path, no equilibrate warning
        slow = air.stagnation(1e-3, 220.0, 300.0)
    assert slow.h_D == 0.0 and 250.0 < slow.T < 300.0


def test_fay_riddell_against_sutton_graves(air):
    s, fs = air.stagnation(RHO, T_INF, V), air.freestream(RHO, T_INF)
    q_fr = gas.fay_riddell(s, air.wall(300.0, s.p), fs.p, R)
    assert q_fr == pytest.approx(2.22e6, rel=3e-2)
    # Fay-Riddell (Le 1.4, fully catalytic) sits above the Sutton-Graves fit at this low-pressure, 72 %-dissociated
    # state (measured 1.29 on 2026-09-18; the spec's 15 % expectation did not hold); with the Lewis term off it is 1.14.
    assert 1.2 < q_fr / gas.sutton_graves(RHO, V, R) < 1.4
    assert 1.05 < gas.fay_riddell(s, air.wall(300.0, s.p), fs.p, R, catalycity=0.0) / gas.sutton_graves(RHO, V, R) < 1.25
    assert gas.fay_riddell(s, air.wall(2000.0, s.p), fs.p, R) < 0.9 * q_fr                     # hot wall
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_gas.py -q`
Expected: collection error (`reentry_model.gas` missing).

- [ ] **Step 3: Write `reentry_model/gas.py`**

```python
"""Equilibrium-air stagnation state and the stagnation-point heating correlations (spec section 6.3).

Thermodynamics and composition: Cantera with `airNASA9.yaml` (NASA-9 polynomials, 200-20000 K, 11 species incl.
ions). That mechanism carries no transport data, so viscosity comes from Blottner's curve fits for N2, O2, NO, N, O
(Blottner, Johnson & Ellis 1971, as tabulated in Gnoffo, Gupta & Shinn, NASA TP-2867, 1989) combined with Wilke's
mixing rule over the equilibrium mole fractions; ions are neglected in the mixture viscosity. Freestream
composition fixed at N2:O2 = 0.79:0.21 by mole; T_inf is clamped to >= 200 K (the polynomial floor; the
enthalpy difference is negligible against V^2/2).
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

AIR = "N2:0.79, O2:0.21"
T_FLOOR = 200.0
CP_AIR = 1004.5              # J/(kg K), perfect-gas air: the stagnation-temperature estimate that decides whether to equilibrate
T_EQUILIBRATE = 1500.0       # K; below this the composition is frozen (no dissociation) and Cantera's equilibrium solver is skipped
PRANDTL = 0.71
LEWIS = 1.4
BLOTTNER = {                       # mu [Pa s] = 0.1 exp((A ln T + B) ln T + C), T in K
    "N2": (0.0268142, 0.3177838, -11.3155513),
    "O2": (0.0449290, -0.0826158, -9.2019475),
    "NO": (0.0436378, -0.0335511, -9.5767430),
    "N": (0.0115572, 0.6031679, -12.4327495),
    "O": (0.0203144, 0.4294404, -11.6031403),
}


@dataclass(frozen=True)
class GasState:
    p: float          # Pa
    T: float          # K
    rho: float        # kg/m3
    h: float          # J/kg (Cantera basis: elements at 298.15 K)
    mu: float         # Pa s
    h_D: float        # J/kg, chemical enthalpy of the free atoms (dissociation enthalpy)
    X: dict           # mole fractions > 1e-6


def blottner_viscosity(species, T):
    A, B, C = BLOTTNER[species]
    lnT = math.log(T)
    return 0.1 * math.exp((A * lnT + B) * lnT + C)


def wilke_viscosity(T, mole_fractions, molar_masses):
    """Mixture viscosity [Pa s] of the Blottner species present in `mole_fractions` (dict name -> X), Wilke's rule."""
    names = [s for s in BLOTTNER if mole_fractions.get(s, 0.0) > 0.0]
    x = np.array([mole_fractions[s] for s in names]); x = x / x.sum()
    mu = np.array([blottner_viscosity(s, T) for s in names])
    m = np.array([molar_masses[s] for s in names])
    total = 0.0
    for i in range(len(names)):
        phi = sum(x[j] * (1.0 + math.sqrt(mu[i] / mu[j]) * (m[j] / m[i]) ** 0.25) ** 2 / math.sqrt(8.0 * (1.0 + m[i] / m[j]))
                  for j in range(len(names)))
        total += x[i] * mu[i] / phi
    return float(total)


class EquilibriumAir:
    """Normal shock + isentropic compression to rest in equilibrium air, and the wall state at (T_w, p_s)."""

    def __init__(self, mechanism="airNASA9.yaml"):
        import cantera as ct                    # imported here: only `--heating physics` needs Cantera
        self.gas = ct.Solution(mechanism)
        self.molar_masses = dict(zip(self.gas.species_names, self.gas.molecular_weights))
        self._h_atoms = {s: self.gas.species(s).thermo.h(298.15) / self.molar_masses[s] for s in ("N", "O")}   # J/kg
        T = np.arange(T_FLOOR, 6000.0 + 1.0, 25.0)
        self._h_table = (T, np.array([self._enthalpy_TP(t, 101325.0) for t in T]))

    def _enthalpy_TP(self, T, p):
        self.gas.TPX = T, p, AIR
        return float(self.gas.enthalpy_mass)

    def air_enthalpy(self, T):
        """Frozen-composition air enthalpy h(T) [J/kg] (ideal gas: pressure-independent), table lookup, vectorised."""
        return np.interp(np.asarray(T, dtype=float), self._h_table[0], self._h_table[1])

    def _state(self):
        g = self.gas
        X = {s: float(x) for s, x in zip(g.species_names, g.X) if x > 1e-6}
        Y = dict(zip(g.species_names, g.Y))
        h_D = sum(Y[s] * self._h_atoms[s] for s in ("N", "O"))
        return GasState(float(g.P), float(g.T), float(g.density), float(g.enthalpy_mass),
                        wilke_viscosity(float(g.T), X, self.molar_masses), float(h_D), X)

    def freestream(self, rho, T):
        self.gas.TDX = max(T, T_FLOOR), rho, AIR
        return self._state()

    def stagnation(self, rho_inf, T_inf, V):
        """Equilibrium stagnation state behind the bow shock for freestream (rho, T, V)."""
        g = self.gas
        T_inf = max(T_inf, T_FLOOR)
        g.TDX = T_inf, rho_inf, AIR
        p_inf, h_inf = float(g.P), float(g.enthalpy_mass)
        equilibrate = T_inf + 0.5 * V * V / CP_AIR > T_EQUILIBRATE     # frozen air below ~1500 K (subsonic/low supersonic states)
        eps = 0.1                                                     # rho_inf / rho_2, Rankine-Hugoniot fixed point
        for _ in range(100):
            u2 = eps * V
            p2 = p_inf + rho_inf * V * (V - u2)
            g.HPX = h_inf + 0.5 * (V * V - u2 * u2), p2, AIR
            if equilibrate:
                g.equilibrate("HP")
            eps_new = rho_inf / float(g.density)
            if abs(eps_new - eps) < 1e-10:
                break
            eps = 0.5 * (eps + eps_new)
        s2, h_s = float(g.entropy_mass), h_inf + 0.5 * V * V             # isentropic compression to rest

        def residual(log_p):
            g.SPX = s2, math.exp(log_p), AIR
            if equilibrate:
                g.equilibrate("SP")
            return float(g.enthalpy_mass) - h_s

        p_s = math.exp(brentq(residual, math.log(p2), math.log(2.0 * p2), xtol=1e-12))
        g.SPX = s2, p_s, AIR
        if equilibrate:
            g.equilibrate("SP")
        return self._state()

    def wall(self, T_w, p):
        """Air at the wall temperature and stagnation pressure (equilibrium; molecular below ~2500 K)."""
        self.gas.TPX = T_w, p, AIR
        if T_w > T_EQUILIBRATE:
            self.gas.equilibrate("TP")
        return self._state()


def dkr(rho, V, radius):
    """Detra-Kemp-Riddell stagnation-point heating [W/m2], SESAM's continuum correlation (facts note s.4)."""
    return 1.1035e8 * radius ** -0.5 * (rho / 1.225) ** 0.5 * (V / 7925.0) ** 3.15


def sutton_graves(rho, V, radius):
    """Sutton-Graves (1971) stagnation-point heating for air [W/m2], cold wall."""
    return 1.7415e-4 * math.sqrt(rho / radius) * V ** 3


def fay_riddell(stag, wall, p_inf, radius, catalycity=1.0, Pr=PRANDTL, Le=LEWIS):
    """Fay-Riddell (1958) equilibrium-boundary-layer stagnation heating [W/m2].

    `catalycity` scales the Lewis-number (atom diffusion + wall recombination) term; 1 is the fully catalytic
    equilibrium form, 0 leaves the Le = 1 value. The frozen-boundary-layer non-catalytic reduction
    (1 - h_D/h_s) is not represented (Step 3 option)."""
    due_dx = math.sqrt(2.0 * (stag.p - p_inf) / stag.rho) / radius
    return 0.763 * Pr ** -0.6 * (wall.rho * wall.mu) ** 0.1 * (stag.rho * stag.mu) ** 0.4 * math.sqrt(due_dx) \
        * (stag.h - wall.h) * (1.0 + (Le ** 0.52 - 1.0) * catalycity * stag.h_D / stag.h)
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_gas.py -q`
Expected: `3 passed` (Cantera init 0.2 s; one stagnation state 23 ms).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/gas.py tests/test_reentry_model_gas.py
git commit -m "reentry_model: equilibrium-air stagnation state (Cantera) with Blottner-Wilke viscosity and the stagnation correlations

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Heating models (SESAM-equivalent verification mode and physics mode)

**Files:**
- Modify: `reentry_model/aero.py` (insert `SesamHeatTable` between `SesamTable` and `SesamErf`)
- Create: `reentry_model/heating.py`
- Test: `tests/test_reentry_model_heating.py`

**Interfaces:**
- Consumes: `gas` (Task 4), `aero.SesamTable`, `trajectory.AeroState` (`freestream.rho/T/p`, `V`, `kn`, `ma`), `constants.GAMMA_AIR`.
- Produces: `HeatingResult(q_conv (nf,), q_stag, q_stag_c, q_stag_fm, blend)` with `total(areas)`; `lees_shape(theta, mach)`, `cosine_shape(theta)`, `matting_bridge(q_fm, q_c, n=1.0)`; `SesamEquivalentHeating(heat_bridging=None, shape_factor=0.27471)`, `PhysicsHeating(stagnation, bridging, matting_n, accommodation, catalycity, air=None)`, `TabulatedHeating` (raises); every model has `.name` and `evaluate(state, theta, T_wall, radius, T_mean=None) -> HeatingResult`; `heating_by_name(name, **options)`; constants `HEATING_NAMES`, `STAGNATION_NAMES`, `BRIDGING_HEAT_NAMES`, `SESAM_SHAPE_FACTOR`, `CP_AIR_SESAM`, `T_WALL_COLD`; `aero.SesamHeatTable()(kn)`.
- The SESAM-equivalent mode is the measured SESAM heat input (amendments 7–8): uniform 0.27471 · q_DKR · F_h(Kn) · hot-wall; **not physical, verification only** — docstrings and the comment at the `np.full(...)` line say so.

- [ ] **Step 1: Insert `SesamHeatTable` into `reentry_model/aero.py`**

Place this class directly after `class SesamTable` (before `class SesamErf`):

```python
class SesamHeatTable:
    """SESAM's convective-heat factor F_h(Kn) = Q / (A x 0.27471 x q_DKR x hot-wall factor), measured on the two no-melt
    US76 reference runs (2026-09-18; 371 rows with Ma >= 5, hot-wall factor > 0.3) in 0.125-decade bins of log10 Kn;
    linear in log10 Kn between bin centres, clamped outside. It is NOT the drag's f(Kn) blend with a free-molecular
    flux: SESAM's transitional heating lies far below both the DKR and the textbook free-molecular values (0.14 x
    continuum at Kn 1, where the drag blend would give 0.5), decaying like Kn^-0.24 in the free-molecular tail
    (0.059 at Kn 40 = 0.78 x the cos-theta average of 1/2 rho V^3). Reproduces the measured rows to rms 0.7 %,
    max 3.4 %. Used only by the SESAM-equivalent verification heating."""
    LOG10_KN = np.array([-2.688, -2.562, -2.438, -2.312, -2.188, -2.062, -1.938, -1.812, -1.688, -1.562, -1.438, -1.312,
                         -1.188, -1.062, -0.938, -0.812, -0.688, -0.562, -0.438, -0.312, -0.188, -0.062, 0.062, 0.188,
                         0.312, 0.438, 0.562, 0.688, 0.812, 0.938, 1.062, 1.188, 1.312, 1.438, 1.562])
    F = np.array([1.0077, 1.0073, 1.0066, 1.0062, 1.0056, 1.0053, 1.0050, 1.0020, 0.9923, 0.9734, 0.9390, 0.8909,
                  0.8228, 0.7405, 0.6548, 0.5639, 0.4645, 0.3673, 0.2868, 0.2176, 0.1693, 0.1442, 0.1371, 0.1334,
                  0.1289, 0.1233, 0.1169, 0.1102, 0.1029, 0.0949, 0.0868, 0.0792, 0.0715, 0.0641, 0.0594])

    def __call__(self, kn):
        if kn <= 0.0:
            return float(self.F[0])
        return float(np.interp(math.log10(kn), self.LOG10_KN, self.F))


```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_reentry_model_heating.py`:

```python
"""heating.py: distributions, bridging, the SESAM-equivalent mode and the physics mode at the 100 mm reference start."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from reentry_model import atmosphere, gas, heating
from reentry_model import trajectory as tj

# 100 mm sphere at the reference start: US76 77.5 km (rho 2.727e-5 kg/m3, T 203.6 K), 7.5 km/s, Kn 0.0298, Ma 26.2
RHO, T_INF, V, R, KN, MA = 2.727e-5, 203.6, 7500.0, 0.05, 0.0298, 26.2


def state(rho=RHO, T=T_INF, V=V, kn=KN, ma=MA):
    fs = atmosphere.Freestream(rho, T, rho / 28.9644e-3 * 8.314462618 * T, 28.9644 * 1.66053906660e-27, np.zeros(3))
    return tj.AeroState(77.5e3, 0.0, 0.0, fs, np.array([V, 0.0, 0.0]), V, kn, ma, 0.9, np.zeros(3), 0.5 * rho * V * V)


def windward_integral(shape, **kw):
    """Integral of the shape over the windward hemisphere divided by the whole sphere area."""
    return 0.5 * quad(lambda th: float(shape(np.array([th]), **kw)[0]) * math.sin(th), 0.0, math.pi / 2)[0]


def test_lees_distribution_limits_and_integrals():
    assert heating.lees_shape(np.array([0.0]), np.inf)[0] == 1.0
    assert heating.lees_shape(np.array([1e-3, 2.0]), 20.0)[0] == pytest.approx(1.0, abs=1e-5)
    assert heating.lees_shape(np.array([2.0, 3.0]), 20.0).tolist() == [0.0, 0.0]                  # leeward
    assert windward_integral(heating.lees_shape, mach=np.inf) == pytest.approx(0.196, abs=1e-3)
    assert windward_integral(heating.lees_shape, mach=10.0) == pytest.approx(0.200, abs=1e-3)
    assert windward_integral(heating.lees_shape, mach=5.0) == pytest.approx(0.209, abs=1e-3)
    assert heating.lees_shape(np.array([math.pi / 4]), np.inf)[0] == pytest.approx(0.5965, abs=1e-3)
    assert windward_integral(heating.cosine_shape) == pytest.approx(0.25, abs=1e-9)
    assert heating.lees_shape(np.array([0.5]), 0.3).tolist() == heating.lees_shape(np.array([0.5]), 1.0).tolist()   # Ma clamp


def test_matting_bridge_limits_and_closed_forms():
    assert heating.matting_bridge(1.0, 1e6) == pytest.approx(1.0, rel=1e-5)               # q_fm << q_c -> q_fm
    assert heating.matting_bridge(1e6, 1.0) == pytest.approx(1.0, rel=1e-9)               # q_fm >> q_c -> q_c
    assert heating.matting_bridge(3.0, 2.0) == pytest.approx(2.0 * (1.0 - math.exp(-1.5)))
    x = math.sqrt(2.0 * 3.0 / 2.0)
    assert heating.matting_bridge(3.0, 2.0, n=2.0) == pytest.approx(2.0 * (1.0 - (1.0 + x) * math.exp(-x)))
    assert heating.matting_bridge(3.0, 2.0, n=1.5) < heating.matting_bridge(3.0, 2.0, n=1.0)
    assert heating.matting_bridge(1.0, 0.0) == 0.0
    ratios = [heating.matting_bridge(q, 1.0) for q in (0.1, 0.5, 1.0, 2.0, 5.0)]
    assert all(b > a for a, b in zip(ratios, ratios[1:])) and ratios[-1] < 1.0


def test_sesam_equivalent_mode_is_uniform_and_totals_to_the_shape_factor(coarse_sphere_mesh):
    surface = coarse_sphere_mesh.surface()
    theta = surface.angles_to([1.0, 0.0, 0.0])
    model = heating.SesamEquivalentHeating()
    res = model.evaluate(state(), theta, np.full(surface.n_patches, 300.0), R)
    F_h = heating.aero.SesamHeatTable()(KN)
    hot_wall_300 = 1.0 - 1004.5 * (300.0 - T_INF) / (0.5 * V * V)                              # 0.9966 at the start
    q_stag = gas.dkr(RHO, V, R) * F_h * hot_wall_300
    assert 0.95 < F_h < 1.0 and res.blend == F_h and res.q_stag == pytest.approx(q_stag) and res.q_stag_c == pytest.approx(gas.dkr(RHO, V, R))
    assert np.all(res.q_conv == res.q_conv[0]) and res.q_conv[0] == pytest.approx(0.27471 * q_stag)
    assert res.total(surface.areas) == pytest.approx(0.27471 * q_stag * surface.area, rel=1e-12)
    assert res.total(surface.areas) == pytest.approx(16244.0, rel=1e-2)                      # SESAM's own t = 0 value, 16243.74 W
    assert model.evaluate(state(rho=0.0), theta, np.full(surface.n_patches, 300.0), R).total(surface.areas) == 0.0
    hot = model.evaluate(state(), theta, np.full(surface.n_patches, 2000.0), R)
    assert np.all(hot.q_conv == pytest.approx(res.q_conv * (1.0 - 1004.5 * (2000.0 - T_INF) / (0.5 * V * V)) / hot_wall_300))   # SESAM's hot-wall factor
    lumped = model.evaluate(state(), theta, np.full(surface.n_patches, 2000.0), R, T_mean=1000.0)
    assert lumped.q_conv[0] == pytest.approx(res.q_conv[0] * (1.0 - 1004.5 * (1000.0 - T_INF) / (0.5 * V * V)) / hot_wall_300)
    slow = model.evaluate(state(rho=1e-3, V=2000.0, kn=1e-5, ma=6.7), theta, np.full(surface.n_patches, 2300.0), R)
    assert np.all(slow.q_conv == 0.0)                                                          # clamped: h_w > h_s
    subsonic = model.evaluate(state(rho=1e-2, V=250.0, kn=1e-6, ma=0.8), theta, np.full(surface.n_patches, 1500.0), R)
    assert subsonic.q_conv[0] == pytest.approx(0.5 * 0.27471 * gas.dkr(1e-2, 250.0, R) * heating.aero.SesamHeatTable()(1e-6))
    fm = model.evaluate(state(rho=4e-8, T=300.0, V=7500.0, kn=35.0, ma=22.0), theta, np.full(surface.n_patches, 300.0), 0.025)
    assert fm.blend == pytest.approx(0.060, abs=3e-3)
    assert 0.6 < fm.total(surface.areas) / (0.25 * 0.5 * 4e-8 * 7500.0 ** 3 * surface.area) < 1.0   # SESAM: 0.78 x the cos-theta average at Kn 35


def test_sesam_heat_table():
    table = heating.aero.SesamHeatTable()
    assert table(1e-4) == pytest.approx(1.0077) and table(100.0) == pytest.approx(0.0594) and table(0.0) == pytest.approx(1.0077)
    assert table(1.0) == pytest.approx(0.14, abs=0.01) and table(0.1) == pytest.approx(0.70, abs=0.05)
    values = [table(10 ** x) for x in np.linspace(-2.7, 1.6, 50)]
    assert all(b <= a + 1e-12 for a, b in zip(values, values[1:]))                            # monotone in Kn


@pytest.fixture(scope="module")
def air():
    pytest.importorskip("cantera")
    return gas.EquilibriumAir()


@pytest.mark.parametrize("stagnation", heating.STAGNATION_NAMES)
def test_physics_mode_distribution_and_totals(air, coarse_sphere_mesh, stagnation):
    surface = coarse_sphere_mesh.surface()
    theta = surface.angles_to([1.0, 0.0, 0.0])
    model = heating.PhysicsHeating(stagnation=stagnation, air=air)
    res = model.evaluate(state(), theta, np.full(surface.n_patches, 300.0), R)
    assert res.q_stag < res.q_stag_c and 0.0 < res.blend < 0.2                               # Matting: q_fm/q_c ~ 3 -> w ~ 0.05
    assert res.blend == pytest.approx(math.exp(-res.q_stag_fm / res.q_stag_c), rel=1e-9)
    expected = res.q_stag * ((1.0 - res.blend) * 0.196 + res.blend * 0.25) * surface.area
    assert res.total(surface.areas) == pytest.approx(expected, rel=2e-2)
    assert np.all(res.q_conv[theta > math.pi / 2] == 0.0) and res.q_conv[int(np.argmin(theta))] == pytest.approx(res.q_stag, rel=1e-2)
    if stagnation == "dkr":
        assert res.q_stag_c == pytest.approx(gas.dkr(RHO, V, R), rel=1e-3)                     # cold wall: factor 1
    hot = model.evaluate(state(), theta, np.full(surface.n_patches, 2000.0), R)
    assert 0.85 < hot.total(surface.areas) / res.total(surface.areas) < 0.97                  # ~7 % hot-wall reduction


def test_physics_mode_free_molecular_limit_and_sesam_bridging(air, coarse_sphere_mesh):
    surface = coarse_sphere_mesh.surface()
    theta = surface.angles_to([1.0, 0.0, 0.0])
    T_wall = np.full(surface.n_patches, 300.0)
    fm = heating.PhysicsHeating(stagnation="sutton-graves", air=air).evaluate(state(rho=1e-10, kn=1e4), theta, T_wall, R)
    assert fm.blend > 0.99 and fm.q_stag == pytest.approx(fm.q_stag_fm, rel=2e-2)
    assert fm.total(surface.areas) == pytest.approx(0.25 * fm.q_stag * surface.area, rel=3e-2)  # cos(theta) shape
    sesam = heating.PhysicsHeating(stagnation="sutton-graves", bridging="sesam-table", air=air).evaluate(state(), theta, T_wall, R)
    assert sesam.blend == pytest.approx(heating.aero.SesamTable()(KN))
    with pytest.raises(ValueError):
        heating.PhysicsHeating(stagnation="magic", air=air)
    with pytest.raises(NotImplementedError):
        heating.TabulatedHeating("table.nc")
    assert isinstance(heating.heating_by_name("sesam"), heating.SesamEquivalentHeating)
```

- [ ] **Step 3: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_heating.py -q`
Expected: collection error (`reentry_model.heating` missing).

- [ ] **Step 4: Write `reentry_model/heating.py`**

```python
"""Convective heating on the sphere's surface patches (spec section 6).

Two models behind one call, `evaluate(state, theta, T_wall, radius) -> HeatingResult`, where `state` is the
trajectory's AeroState (freestream rho/T/p, airspeed V, Kn, Ma), `theta` the angle of every patch's outward normal
to the direction of motion (0 at the stagnation point), `T_wall` every patch's temperature and `radius` the nose
radius. `q_conv` is W/m2 into the body per patch.

SesamEquivalentHeating -- VERIFICATION DEVICE, NOT A PHYSICAL MODEL. It reproduces SESAM's lumped-body heat input as
measured on the two no-melt US76 references (2026-09-18, superseding the (1 - f) q_DKR + f q_FM decomposition of
facts note s.4): q = 0.27471 x q_DKR x F_h(Kn) x max(0, 1 - c_p,air (T_lumped - T_inf) / (V^2/2)) for Ma >= 1, with
SESAM's own heat bridging F_h (aero.SesamHeatTable: 1.005 in the continuum, 0.14 at Kn 1, 0.06 at Kn 40) and its
hot-wall factor (c_p 1004.5 J/kg/K: rms 0.008, max 0.015 of the ratio over 189 rows up to 2450 K; the factor is
1/2 with no hot-wall term below Ma 1), the surface-average factor 0.27471 applied uniformly to EVERY patch, front
and back. Real heating is concentrated on the windward face (Lees: 0.196 of the stagnation flux averaged over the
sphere, zero leeward). The uniform distribution exists so that our conduction, time stepping, material curves and
coupling can be checked against SESAM's lumped temperature and heat totals with no distribution question in between.

PhysicsHeating -- the model proper: continuum stagnation flux from Fay-Riddell (equilibrium air, `gas.py`),
Sutton-Graves or DKR (with the hot-wall factor (h_s - h_w)/(h_s - h_w,300K) where the correlation is cold-wall),
free-molecular stagnation flux A_cq rho V (h_s - h_w) (Matting Eq. 1), Matting's bridging (facts note s.14) or
SESAM's f(Kn), and the distribution (1 - w) Lees(theta, Ma) + w cos(theta) with the free-molecular weight
w = 1 - q_stag/q_c (Matting) or f(Kn) (SESAM table); per patch the flux is scaled by its own hot-wall factor
relative to the stagnation patch. Fixed attitude: theta is whatever the caller passes.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.special import gamma, gammainc

from . import aero, gas
from .constants import GAMMA_AIR

SESAM_SHAPE_FACTOR = 0.27471        # ATDB_SPHERE continuum heat factor: surface average / stagnation flux
CP_AIR_SESAM = 1004.5               # J/(kg K): SESAM's hot-wall enthalpy h_w - h_inf = c_p (T_w - T_inf) (best fit 985; 1004.5 within 1.5 %)
T_WALL_COLD = 300.0                 # K, wall temperature the cold-wall correlations are referred to
HEATING_NAMES = ("sesam", "physics")
STAGNATION_NAMES = ("fay-riddell", "sutton-graves", "dkr")
BRIDGING_HEAT_NAMES = ("matting", "sesam-table")


@dataclass
class HeatingResult:
    q_conv: np.ndarray        # W/m2 per patch, positive into the body
    q_stag: float             # bridged stagnation-point flux
    q_stag_c: float           # continuum stagnation value
    q_stag_fm: float          # free-molecular stagnation value
    blend: float              # free-molecular weight of the distribution: w (physics) or f(Kn) (sesam)

    def total(self, areas):
        return float((self.q_conv * areas).sum())


def lees_shape(theta, mach):
    """Lees (1956) laminar sphere distribution q(theta)/q_stag for theta <= pi/2 (zero leeward), finite Mach through
    eps = 1/(gamma Ma^2); Ma is clamped to >= 1. Integrates over the windward hemisphere to 0.196 (Ma -> inf),
    0.200 (Ma 10), 0.209 (Ma 5) of q_stag x sphere area."""
    theta = np.asarray(theta, dtype=float)
    eps = 0.0 if not np.isfinite(mach) else 1.0 / (GAMMA_AIR * max(mach, 1.0) ** 2)
    D = (1.0 - eps) * (theta ** 2 - 0.5 * theta * np.sin(4.0 * theta) + 0.125 * (1.0 - np.cos(4.0 * theta))) \
        + 4.0 * eps * (theta ** 2 - theta * np.sin(2.0 * theta) + 0.5 * (1.0 - np.cos(2.0 * theta)))
    small = theta < 1e-4                                       # limit theta -> 0 is exactly 1
    with np.errstate(divide="ignore", invalid="ignore"):
        q = 2.0 * theta * np.sin(theta) * ((1.0 - eps) * np.cos(theta) ** 2 + eps) / np.sqrt(D)
    q = np.where(small, 1.0, q)
    return np.where(theta <= 0.5 * math.pi, q, 0.0)


def cosine_shape(theta):
    """Free-molecular distribution cos(theta) on the windward hemisphere (integrates to 0.25 of q_stag x area)."""
    theta = np.asarray(theta, dtype=float)
    return np.where(theta <= 0.5 * math.pi, np.cos(theta), 0.0)


def matting_bridge(q_fm, q_c, n=1.0):
    """Matting (1971) Eq. 18: q = q_c P(n, [Gamma(n+1) q_fm/q_c]^(1/n)) with P the regularised lower incomplete gamma
    function; n = 1 gives q_c [1 - exp(-q_fm/q_c)]. Limits: q -> q_fm when q_fm << q_c, q -> q_c when q_fm >> q_c."""
    if q_c <= 0.0:
        return 0.0
    if n == 1.0:
        return q_c * (1.0 - math.exp(-q_fm / q_c))
    return q_c * float(gammainc(n, (gamma(n + 1.0) * q_fm / q_c) ** (1.0 / n)))


class SesamEquivalentHeating:
    """SESAM's heating on every patch: q = 0.27471 x q_DKR x F_h(Kn) x max(0, 1 - c_p (T_lumped - T_inf) / (V^2/2))
    for Ma >= 1, and 0.5 x 0.27471 x q_DKR x F_h(Kn) with no hot-wall term below Ma 1 (both measured, module docstring).

    NOT PHYSICAL -- a verification device: the uniform factor reproduces SESAM's tumbling-average lumped-body input so
    that the conduction/coupling can be compared with SESAM's lumped temperature. No distribution over theta; the
    hot-wall factor uses the body's lumped (energy-equivalent) temperature `T_mean` as SESAM does, falling back to
    the mean of `T_wall`. HeatingResult: q_stag is the effective stagnation flux q_DKR x F_h x hot-wall, q_stag_c the
    DKR value, q_stag_fm the textbook 1/2 rho V^3 for reference only, blend = F_h(Kn)."""

    name = "sesam"

    def __init__(self, heat_bridging=None, shape_factor=SESAM_SHAPE_FACTOR):
        self.heat_bridging = heat_bridging or aero.SesamHeatTable()
        self.shape_factor = shape_factor

    def evaluate(self, state, theta, T_wall, radius, T_mean=None):
        rho, V = state.freestream.rho, state.V
        n_patches = len(np.asarray(theta))
        if rho <= 0.0 or V <= 0.0:
            return HeatingResult(np.zeros(n_patches), 0.0, 0.0, 0.0, 1.0)
        q_c = gas.dkr(rho, V, radius)
        F_h = self.heat_bridging(state.kn) if np.isfinite(state.kn) else float(self.heat_bridging.F[-1])
        T_lumped = float(np.mean(T_wall)) if T_mean is None else float(T_mean)
        if state.ma >= 1.0:
            hot_wall = max(0.0, 1.0 - CP_AIR_SESAM * (T_lumped - state.freestream.T) / (0.5 * V * V))
        else:
            hot_wall = 0.5                       # SESAM below Ma 1: half the continuum value, no hot-wall factor (measured)
        q_stag = q_c * F_h * hot_wall
        # verification only: SESAM's surface-average factor on EVERY patch, windward and leeward alike (not physical)
        return HeatingResult(np.full(n_patches, self.shape_factor * q_stag), q_stag, q_c, 0.5 * rho * V ** 3, F_h)


class PhysicsHeating:
    """Stagnation correlation + Matting/SESAM bridging + Lees/cos(theta) distribution (module docstring)."""

    name = "physics"

    def __init__(self, stagnation="fay-riddell", bridging="matting", matting_n=1.0, accommodation=0.8,
                 catalycity=1.0, air=None):
        if stagnation not in STAGNATION_NAMES:
            raise ValueError("stagnation must be one of {}, got {!r}".format(STAGNATION_NAMES, stagnation))
        if bridging not in BRIDGING_HEAT_NAMES:
            raise ValueError("bridging must be one of {}, got {!r}".format(BRIDGING_HEAT_NAMES, bridging))
        self.stagnation, self.bridging, self.matting_n = stagnation, bridging, matting_n
        self.accommodation, self.catalycity = accommodation, catalycity
        self.air = air or gas.EquilibriumAir()
        self.sesam_table = aero.SesamTable()

    def stagnation_fluxes(self, state, T_wall_stag, radius):
        """(q_c, q_fm, h_s) at the stagnation patch for the freestream `state` and the stagnation wall temperature."""
        fs = state.freestream
        rho, V = fs.rho, state.V
        h_s = float(self.air.air_enthalpy(fs.T)) + 0.5 * V * V
        h_w = float(self.air.air_enthalpy(T_wall_stag))
        hot_wall = (h_s - h_w) / (h_s - float(self.air.air_enthalpy(T_WALL_COLD)))
        if self.stagnation == "fay-riddell":
            stag = self.air.stagnation(rho, fs.T, V)
            wall = self.air.wall(T_wall_stag, stag.p)
            q_c, h_s = gas.fay_riddell(stag, wall, fs.p, radius, self.catalycity), stag.h
        elif self.stagnation == "sutton-graves":
            q_c = gas.sutton_graves(rho, V, radius) * hot_wall
        else:
            q_c = gas.dkr(rho, V, radius) * hot_wall
        q_fm = self.accommodation * rho * V * (h_s - h_w)
        return q_c, q_fm, h_s

    def evaluate(self, state, theta, T_wall, radius, T_mean=None):
        theta, T_wall = np.asarray(theta, dtype=float), np.asarray(T_wall, dtype=float)
        rho, V = state.freestream.rho, state.V
        if rho <= 0.0 or V <= 0.0:
            return HeatingResult(np.zeros(theta.size), 0.0, 0.0, 0.0, 1.0)
        i_stag = int(np.argmin(theta))
        q_c, q_fm, h_s = self.stagnation_fluxes(state, float(T_wall[i_stag]), radius)
        if self.bridging == "matting":
            q_stag = matting_bridge(q_fm, q_c, self.matting_n)
            w = 1.0 - q_stag / q_c if q_c > 0.0 else 1.0
        else:
            w = self.sesam_table(state.kn) if np.isfinite(state.kn) else 1.0
            q_stag = (1.0 - w) * q_c + w * q_fm
        shape = (1.0 - w) * lees_shape(theta, state.ma) + w * cosine_shape(theta)
        h_w = self.air.air_enthalpy(T_wall)
        per_patch_hot_wall = (h_s - h_w) / (h_s - h_w[i_stag])
        return HeatingResult(q_stag * shape * per_patch_hot_wall, float(q_stag), float(q_c), float(q_fm), float(w))


class TabulatedHeating:
    """Reserved (spec 6.3): q_stag and q(theta)/q_stag interpolated from an offline CFD/DSMC table over (V, rho, T_w).
    Only the interface exists; constructing it raises NotImplementedError."""

    name = "tabulated"

    def __init__(self, table_path):
        raise NotImplementedError("TabulatedHeating is an interface slot; no table format is defined yet")

    def evaluate(self, state, theta, T_wall, radius, T_mean=None):
        raise NotImplementedError


def heating_by_name(name, **options):
    if name == "sesam":
        return SesamEquivalentHeating(**options)
    if name == "physics":
        return PhysicsHeating(**options)
    raise ValueError("heating must be one of {}, got {!r}".format(HEATING_NAMES, name))
```

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_heating.py tests/test_reentry_model_aero.py -q`
Expected: heating `8 passed` (the SESAM-equivalent total at the 100 mm start is 16.2 kW, SESAM's own t = 0 value is 16243.74 W), aero unchanged.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/aero.py reentry_model/heating.py tests/test_reentry_model_heating.py
git commit -m "reentry_model: patch heating models -- measured SESAM-equivalent verification mode and Fay-Riddell/Matting/Lees physics mode

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: Body stepper protocol, ThermalBody and the trajectory stepper

**Files:**
- Modify: `reentry_model/body.py` (replace the file)
- Modify: `reentry_model/trajectory.py` (replace the file)
- Modify: `tests/test_reentry_model_trajectory.py` (two lines in `test_constant_body`, two new tests)

**Interfaces:**
- Consumes: `mesh.VolumeMesh`, `material.Material`, a `ThermalSolver` (Task 3), `heating.HeatingResult` (Task 5, only its `q_conv`).
- Produces: `body.Body` protocol `mass(t)`, `advance(t, dt, loads)`, `surface_temperature()`, `mean_temperature()`, `energy()`, `field()`; `ConstantBody(mass_kg, temperature_K)`; `ThermalBody(mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1,0,0))` with attributes `surface`, `theta`, `i_stag`, `i_back`, `i_centre`, `volume`, `emissivity`, `T_ambient`, `energy0`, `integrated_heat`, `absorbed_heat`, `radiated_heat`, `iterations`, `last` and methods `radiated_power()`, `energy_balance_residual()`, `surface_stats()` (keys `surface_T_max_K`, `surface_T_min_K`, `surface_T_mean_K`, `T_stagnation_K`, `T_back_K`, `T_centre_K`). `trajectory.Simulator` gains `t`, `y`, `end_reason`, `nfev`, `advance(dt) -> advanced time`, `sample_row(t, y, aero) -> dict` and `results(history, nfev, runtime)` (was `_results`); `run()` is unchanged in behaviour; `write_history_csv` writes every column (Step 1 columns first), `read_history_csv` reads every column.
- The old `Body.temperature(t)` / `on_step(...)` are gone; `run()` uses `body.mean_temperature()`.

- [ ] **Step 1: Update the trajectory tests**

In `tests/test_reentry_model_trajectory.py`, replace the last two lines of `test_constant_body`

```python
    assert b.mass(0.0) == 1.5 and b.mass(100.0) == 1.5 and b.temperature(5.0) == 310.0
    assert b.on_step(0.0, None, None, None) is None
```

with

```python
    assert b.mass(0.0) == 1.5 and b.mass(100.0) == 1.5 and b.mean_temperature() == 310.0
    assert b.advance(0.0, 0.5, None) is None and b.field() is None and b.energy() == 0.0
    assert b.surface_temperature().tolist() == [310.0]
```

and append at the end of the file:

```python
def test_advance_reproduces_run_and_truncates_at_the_ground():
    """The stepper (a fresh DOP853 solve per macro step) vs the single dense-output solve of run(): < 0.01 m/s, < 0.1 m."""
    atm = atmosphere.US76TableAtmosphere()
    hist = make(R100, atm, cadence=0.5, t_max=60.0).run()
    sim = make(R100, atm, t_max=60.0)
    assert sim.t == 0.0 and sim.end_reason is None
    rows = [sim.sample_row(sim.t, sim.y, sim.aero_state(sim.t, sim.y[:3], sim.y[3:]))]
    while sim.end_reason is None:
        assert sim.advance(0.5) == pytest.approx(0.5)
        rows.append(sim.sample_row(sim.t, sim.y, sim.aero_state(sim.t, sim.y[:3], sim.y[3:])))
    assert sim.end_reason == "t_max" and sim.t == 60.0 and len(rows) == len(hist) and sim.nfev > 0
    V = np.array([r["velocity_kms"] for r in rows]) * 1e3
    h = np.array([r["altitude_km"] for r in rows]) * 1e3
    assert np.abs(V - hist.columns["velocity_kms"] * 1e3).max() < 0.01
    assert np.abs(h - hist.columns["altitude_km"] * 1e3).max() < 0.1
    with pytest.raises(RuntimeError):
        sim.advance(0.5)
    # ground event truncates the final step
    low = tj.InitialState(300.0, 200.0, math.radians(-60.0), R100.heading, R100.lat, R100.lon, EPOCH)
    sim = make(low, atm)
    advanced = sim.advance(10.0)
    assert sim.end_reason == "ground" and 0.5 < advanced < 1.5 and abs(earth.ecef_to_geodetic(sim.y[:3])[0]) < 1e-3


def test_history_csv_round_trips_extra_columns(tmp_path):
    hist = make(R100, atmosphere.US76TableAtmosphere(), cadence=5.0, t_max=10.0).run()
    hist.columns["convective_heat_W"] = np.array([1.0, 2.0, 3.0])
    tj.write_history_csv(hist, tmp_path / "h.csv")
    with open(tmp_path / "h.csv") as fh:
        header = fh.readline().strip().split(",")
    assert header[:len(tj.CSV_COLUMNS)] == tj.CSV_COLUMNS and header[-1] == "convective_heat_W"
    back = tj.read_history_csv(tmp_path / "h.csv")
    assert back.columns["convective_heat_W"].tolist() == [1.0, 2.0, 3.0] and len(back) == 3
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_trajectory.py -q`
Expected: `test_constant_body` fails (`mean_temperature` missing), the two new tests fail (`advance`/`sample_row` missing, extra CSV column dropped).

- [ ] **Step 3: Replace `reentry_model/body.py`**

```python
"""The body the trajectory carries, as a stepper (spec section 4). ConstantBody keeps mass and one temperature
(Step 1 behaviour, `--thermal none`); ThermalBody wraps a mesh, a material and a conduction solver and advances the
temperature field with each macro step's per-patch convective loads, keeping the heat bookkeeping SESAM reports."""
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


class ThermalBody:
    """Mesh + material + ThermalSolver. `v_hat` is the direction of motion in the body frame (fixed attitude: the
    stagnation patch is the one whose normal is along v_hat, spec 13.1)."""

    def __init__(self, mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        self.mesh, self.material, self.solver, self.mass_kg = mesh, material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient = float(T_ambient)
        self.surface = mesh.surface()
        self.theta = self.surface.angles_to(v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(v_hat), self.surface.patch_toward(-np.asarray(v_hat, dtype=float))
        self.i_centre = mesh.centre_node()
        self.volume = mesh.volume()
        self.solver.setup(mesh, material, self.emissivity)
        self.solver.set_temperature(T0)
        self.energy0 = self.solver.energy()
        self.integrated_heat = 0.0       # J, integral of Q_conv dt (SESAM's integrated_heat_J)
        self.absorbed_heat = 0.0         # J, integral of (Q_conv - Q_rad) dt
        self.radiated_heat = 0.0         # J
        self.iterations = []             # Newton iterations per step
        self.last = None                 # thermal.StepResult of the last step

    def mass(self, t):
        return self.mass_kg

    def advance(self, t, dt, loads):
        res = self.solver.step(dt, loads.q_conv, self.T_ambient)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res

    def surface_temperature(self):
        return self.surface.facet_mean(self.solver.temperature())

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / (self.material.rho * self.volume)))

    def energy(self):
        return self.solver.energy()

    def field(self):
        return self.solver.temperature()

    def radiated_power(self):
        return self.solver.radiated_power(self.T_ambient)

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat) / absorbed heat: zero for an exact discrete balance."""
        return (self.energy() - self.energy0 - self.absorbed_heat) / self.absorbed_heat if self.absorbed_heat else 0.0

    def surface_stats(self):
        T, Tf = self.solver.temperature(), self.surface_temperature()
        return {"surface_T_max_K": float(Tf.max()), "surface_T_min_K": float(Tf.min()),
                "surface_T_mean_K": float((Tf * self.surface.areas).sum() / self.surface.area),
                "T_stagnation_K": float(Tf[self.i_stag]), "T_back_K": float(Tf[self.i_back]), "T_centre_K": float(T[self.i_centre])}
```

- [ ] **Step 4: Replace `reentry_model/trajectory.py`**

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
        if initial.altitude >= settings.escape_altitude:
            raise ValueError("initial altitude {:.0f} m must be below the escape altitude {:.0f} m".format(
                initial.altitude, settings.escape_altitude))
        self.initial, self.body, self.atmosphere = initial, body, atmosphere
        self.tables, self.bridging, self.settings = tables, bridging, settings
        self.area = math.pi * settings.diameter ** 2 / 4.0
        self.omega = np.array([0.0, 0.0, OMEGA_EARTH]) if settings.rotating_frame else np.zeros(3)
        self.r0 = earth.geodetic_to_ecef(initial.altitude, initial.lat, initial.lon)
        self.v0 = earth.velocity_from_flight_angles(initial.velocity, initial.flight_path, initial.heading,
                                                    initial.lat, initial.lon)
        # stepper state (advance): time, ECEF state vector, end reason once an event or t_max is reached
        self.t, self.y, self.end_reason, self.nfev = 0.0, self.initial_state_vector(), None, 0

    def initial_state_vector(self):
        return np.concatenate([self.r0, self.v0])

    def aero_state(self, t, r, v):
        h, lat, lon = earth.ecef_to_geodetic(r)
        # DOP853's adaptive stages can trial-evaluate the RHS a little beyond the ground or the
        # top of a hard-bounded table (e.g. US76's [0, 150000] m) before backing off. Clamp ONLY
        # the altitude used for the atmosphere.state() lookup to [ground_altitude, escape_altitude]
        # so that transient excursion doesn't raise; `h` itself (and thus AeroState.h, used for
        # the reported/history altitude) stays the TRUE, unclamped geodetic altitude.
        # The upper clamp only ever sees this trial-stage overshoot: __init__ rejects any initial
        # altitude at or above escape_altitude, so a real (non-trial) state can't reach here already
        # past the top of the table.
        h_atm = min(max(h, self.settings.ground_altitude), self.settings.escape_altitude)
        fs = self.atmosphere.state(t, h_atm, lat, lon)
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

    def advance(self, dt):
        """One macro step from the stored state: DOP853 over [t, t + dt] at the Step 1 tolerances, a fresh solve_ivp
        per step (spec section 7). A ground/escape event or t_max truncates the step and sets `end_reason`.
        Returns the time actually advanced."""
        if self.end_reason is not None:
            raise RuntimeError("the flight already ended ({})".format(self.end_reason))
        s = self.settings
        t_target = min(self.t + dt, s.t_max)
        sol = solve_ivp(self.rhs, (self.t, t_target), self.y, method="DOP853", rtol=s.rtol,
                        atol=[s.atol_position] * 3 + [s.atol_velocity] * 3, events=self._events())
        if not sol.success:
            raise RuntimeError("integration failed: {}".format(sol.message))
        self.nfev += int(sol.nfev)
        t_new, self.y = float(sol.t[-1]), sol.y[:, -1].copy()
        if sol.t_events[0].size:
            self.end_reason = "ground"
        elif sol.t_events[1].size:
            self.end_reason = "escape"
        elif t_new >= s.t_max:
            self.end_reason = "t_max"
        advanced, self.t = t_new - self.t, t_new
        return advanced

    def sample_row(self, t, y, a):
        """The Step 1 history columns for state y at time t with its AeroState a (the body supplies temperature/mass)."""
        r, v = y[:3], y[3:]
        # velocity_kms/flight_path_deg/heading_deg are kinematic and use the ground-relative
        # (rotating-frame) velocity v -- the quantity SESAM reports relative to the rotating
        # atmosphere -- while mach/knudsen/drag/dynamic_pressure_Pa use aero_state's
        # wind-relative v_rel; the two coincide unless a wind model is active.
        V, gamma, heading = earth.flight_angles(v, a.lat, a.lon)
        return {
            "time_s": t, "altitude_km": a.h / 1e3, "velocity_kms": V / 1e3,
            "temperature_K": self.body.mean_temperature(), "mass_kg": self.body.mass(t),
            "thick_mm": self.settings.diameter * 500.0,
            "lat_deg": math.degrees(a.lat), "lon_deg": math.degrees(a.lon),
            "downrange_km": earth.great_circle_distance(self.initial.lat, self.initial.lon, a.lat, a.lon) / 1e3,
            "flight_path_deg": math.degrees(gamma), "heading_deg": math.degrees(heading),
            "drag": a.cd, "lift": 0.0, "side": 0.0, "knudsen": a.kn, "mach": a.ma,
            "density_kgm3": a.freestream.rho, "dynamic_pressure_Pa": a.q_dyn,
            "load_factor_g": float(np.linalg.norm(a.a_drag)) / G0,     # = SESAM's column (drag / m g0, verified)
        }

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
            row = self.sample_row(t, y, self.aero_state(t, y[:3], y[3:]))
            for k in CSV_COLUMNS:
                cols[k].append(row[k])
            states.append(y)
        columns = {k: np.array(v, dtype=float) for k, v in cols.items()}
        history = History(columns, np.array(states), end_reason)
        history.results = self.results(history, int(sol.nfev), time.perf_counter() - started)
        return history

    def results(self, history, nfev, runtime):
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
            "rhs_evaluations": nfev,
            "runtime_s": runtime,
        }
        return results


def write_history_csv(history, path):
    """All of the history's columns, the Step 1 columns first, then any thermal columns (coupled runs)."""
    names = CSV_COLUMNS + [k for k in history.columns if k not in CSV_COLUMNS]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(names)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(history.columns[k][i]) for k in names])


def read_history_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    columns = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return History(columns, np.zeros((len(rows), 6)), "unknown")


def write_run_json(path, doc):
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
```

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_trajectory.py tests/test_reentry_model_cli.py tests/test_reentry_model_compare.py tests/test_reentry_model_reference.py -q -m "not reference"`
Expected: all pass (the stepper reproduces `run()` to < 0.01 m/s and < 0.1 m over 60 s; the Step 1 CLI and compare tests are unaffected).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/body.py reentry_model/trajectory.py tests/test_reentry_model_trajectory.py
git commit -m "reentry_model: Body stepper protocol with ThermalBody, Simulator.advance(dt), all-column history IO

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: Lockstep coupled driver with VTK output

**Files:**
- Create: `reentry_model/coupled.py`
- Test: `tests/test_reentry_model_coupled.py`

**Interfaces:**
- Consumes: `trajectory.Simulator` stepper (Task 6), `body.ThermalBody` (Task 6), a heating model (Task 5), pyvista (VTK writing, imported inside `write_vtk_frame`).
- Produces: `THERMAL_COLUMNS` (12 names), `PVD_TEMPLATE`, `CoupledSettings(dt=0.5, frames_every=0, output_dir=None, frames=[])`, `CoupledRun(sim, body, heating_model, settings=None)` with `loads_at(t, y)`, `row(t, y, a, loads)`, `run() -> History` (columns = Step 1 columns + `THERMAL_COLUMNS`; `results` = Step 1 results + `peak_surface_T_K`, `time_of_peak_surface_T_s`, `altitude_of_peak_surface_T_km`, `peak_mean_T_K`, `time_of_peak_mean_T_s`, `peak_convective_heat_W`, `time_of_peak_heating_s`, `altitude_of_peak_heating_km`, `integrated_heat_J`, `absorbed_heat_J`, `radiated_heat_J`, `energy_balance_residual`, `n_macro_steps`, `mean_newton_iterations`, `n_frames`), `thermal_results(history)`; `write_vtk_frame(output_dir, k, body, loads)` writes `field_<k>.vtu` (nodal `T`) and `surface_<k>.vtp` (point `T`, cell `q_conv`, `q_rad`, `T_patch`); `run()` writes `field.pvd` and `surface.pvd` collections when frames were written.
- `rad_cooling_W` is negative (SESAM's sign); `integrated_heat_J` = ∫Q_conv dt (SESAM's meaning); `absorbed_heat_J` = ∫(Q_conv − Q_rad)dt; the heating at row t uses the state at t and the wall temperatures at t; `T_mean=body.mean_temperature()` is passed to the heating model.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_coupled.py`:

```python
"""coupled.py: the lockstep loop with an inert body (Step 1 regression) and with the thermal body."""
import math
import os
from datetime import datetime

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, coupled, heating, material, mesh, thermal
from reentry_model import trajectory as tj

EPOCH = datetime(2024, 8, 1, 12, 53, 7)
R100 = tj.InitialState(7500.0, 77500.133, math.radians(-0.959331), math.radians(347.168296),
                       math.radians(29.546067), math.radians(-82.134333), EPOCH)
MASS_100MM = body.sphere_mass(0.1, 2813.0)


def simulator(the_body, **settings):
    return tj.Simulator(R100, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(),
                        tj.Settings(diameter=0.1, **settings))


def thermal_body(the_mesh, T0=300.0, **solver_options):
    mat = material.Material.from_drama_json()
    return body.ThermalBody(the_mesh, mat, thermal.thermal_solver("skfem", **solver_options), MASS_100MM, T0=T0)


def test_coupled_full_flight_trajectory_matches_run(coarse_sphere_mesh):
    """The coupled loop (constant mass) over the whole 100 mm US76 flight, dt 2 s, vs the dense-output run() of Step 1:
    < 0.01 m/s, < 0.1 m (spec section 10, coupling)."""
    reference = simulator(body.ConstantBody(MASS_100MM), cadence=2.0).run()
    b = thermal_body(coarse_sphere_mesh)
    hist = coupled.CoupledRun(simulator(b), b, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=2.0)).run()
    assert hist.end_reason == "ground" and reference.end_reason == "ground"
    assert abs(hist.columns["time_s"][-1] - reference.columns["time_s"][-1]) < 0.01
    t = hist.columns["time_s"]
    V_ref = np.interp(t, reference.columns["time_s"], reference.columns["velocity_kms"]) * 1e3
    h_ref = np.interp(t, reference.columns["time_s"], reference.columns["altitude_km"]) * 1e3
    assert np.abs(hist.columns["velocity_kms"] * 1e3 - V_ref).max() < 0.01
    assert np.abs(hist.columns["altitude_km"] * 1e3 - h_ref).max() < 0.1
    assert set(coupled.THERMAL_COLUMNS) <= set(hist.columns) and hist.columns["temperature_K"].max() > 500.0
    assert hist.results["end_reason"] == "ground" and hist.results["n_macro_steps"] == len(hist) - 1
    assert abs(hist.results["energy_balance_residual"]) < 1e-6 and hist.results["peak_surface_T_K"] > hist.results["peak_mean_T_K"]


def test_thermal_body_bookkeeping(coarse_sphere_mesh):
    b = thermal_body(coarse_sphere_mesh)
    assert b.mass(0.0) == MASS_100MM and b.mean_temperature() == pytest.approx(300.0) and b.field().shape == (coarse_sphere_mesh.n_nodes,)
    assert b.theta[b.i_stag] < 0.1 and b.theta[b.i_back] > 3.0 and b.energy() == pytest.approx(b.energy0)
    q = np.full(b.surface.n_patches, 1e5)
    loads = heating.HeatingResult(q, 1e5 / 0.27471, 0.0, 0.0, 0.0)
    E0 = b.energy()
    b.advance(0.5, 0.5, loads)
    Q = 1e5 * b.surface.area
    assert b.integrated_heat == pytest.approx(0.5 * Q) and b.last.Q_conv == pytest.approx(Q)
    assert b.absorbed_heat == pytest.approx(0.5 * (Q - b.last.Q_rad)) and b.radiated_heat == pytest.approx(0.5 * b.last.Q_rad)
    assert b.energy() - E0 == pytest.approx(b.absorbed_heat, rel=1e-5) and abs(b.energy_balance_residual()) < 1e-5
    stats = b.surface_stats()
    assert stats["surface_T_max_K"] > stats["T_centre_K"] > 299.9 and 300.0 < stats["surface_T_mean_K"] < stats["surface_T_max_K"] + 1e-9
    assert b.mean_temperature() > 300.0 and b.surface_temperature().shape == (b.surface.n_patches,)


def test_coupled_flight_with_the_thermal_body(coarse_sphere_mesh, tmp_path):
    """30 s of the 100 mm flight in SESAM-equivalent mode on the coarse mesh, VTK every 20 steps."""
    b = thermal_body(coarse_sphere_mesh)
    settings = coupled.CoupledSettings(dt=0.5, frames_every=20, output_dir=str(tmp_path / "vtk"))
    hist = coupled.CoupledRun(simulator(b, t_max=30.0), b, heating.SesamEquivalentHeating(), settings).run()
    c = hist.columns
    assert hist.end_reason == "t_max" and len(hist) == 61 and c["time_s"][-1] == 30.0
    assert c["convective_heat_W"][0] == pytest.approx(16244.0, rel=1e-2)        # SESAM's own t = 0 value is 16243.74 W
    assert np.all(np.diff(c["integrated_heat_J"]) > 0) and c["integrated_heat_J"][0] == 0.0
    assert c["integrated_heat_J"][-1] == pytest.approx(np.sum(c["convective_heat_W"][1:] * np.diff(c["time_s"])), rel=1e-9)
    assert np.all(c["rad_cooling_W"] < 0.0) and c["rad_cooling_W"][0] == pytest.approx(-0.4 * 5.670374419e-8 * b.surface.area * 300.0 ** 4, rel=1e-6)
    assert c["temperature_K"][-1] > c["temperature_K"][0] + 20.0 and c["surface_T_max_K"][-1] > c["temperature_K"][-1] > c["T_centre_K"][-1]
    assert np.all(c["heating_blend_f"] == [aero.SesamHeatTable()(kn) for kn in c["knudsen"]])
    assert hist.results["energy_balance_residual"] == pytest.approx(0.0, abs=1e-6)
    assert hist.results["integrated_heat_J"] == c["integrated_heat_J"][-1] and hist.results["n_frames"] == 4
    for k in range(4):
        assert os.path.isfile(tmp_path / "vtk" / "field_{}.vtu".format(k)) and os.path.isfile(tmp_path / "vtk" / "surface_{}.vtp".format(k))
    pvd = open(tmp_path / "vtk" / "field.pvd").read()
    assert 'timestep="10"' in pvd and 'file="field_1.vtu"' in pvd and os.path.isfile(tmp_path / "vtk" / "surface.pvd")


def test_frames_need_an_output_dir(coarse_sphere_mesh):
    with pytest.raises(ValueError):
        coupled.CoupledRun(simulator(body.ConstantBody(MASS_100MM)), body.ConstantBody(MASS_100MM),
                           heating.SesamEquivalentHeating(), coupled.CoupledSettings(frames_every=5))
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py -q`
Expected: collection error (`reentry_model.coupled` missing).

- [ ] **Step 3: Write `reentry_model/coupled.py`**

```python
"""Lockstep coupling of the trajectory stepper with the thermal body (spec section 7).

Per macro step of dt: Simulator.advance(dt) (DOP853, events truncate the last step) -> aero state at the end of the
step -> heating with that state and the wall temperatures at the start of the step -> ThermalBody.advance (radiation
implicit in the solver) -> history row, and every `frames_every` steps a VTK frame (nodal T on the volume mesh,
q_conv / q_rad / T per patch on the surface). First-order operator splitting; the dt-halving test bounds its error.
Mass is constant in Step 2; the loop already carries the body's mass and the mesh so Step 3 can change both."""
import os
import time
from dataclasses import dataclass, field

import numpy as np

from . import trajectory as tj

THERMAL_COLUMNS = ["convective_heat_W", "rad_cooling_W", "integrated_heat_J", "absorbed_heat_J",
                   "surface_T_max_K", "surface_T_min_K", "surface_T_mean_K", "T_stagnation_K", "T_back_K",
                   "T_centre_K", "q_stag_Wm2", "heating_blend_f"]
PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'


@dataclass
class CoupledSettings:
    dt: float = 0.5                  # s, macro step
    frames_every: int = 0            # VTK frame every n macro steps (0: none)
    output_dir: str = None           # directory of the VTK series (required when frames_every > 0)
    frames: list = field(default_factory=list)


class CoupledRun:
    def __init__(self, sim, body, heating_model, settings=None):
        self.sim, self.body, self.heating = sim, body, heating_model
        self.settings = settings or CoupledSettings()
        if self.settings.frames_every and not self.settings.output_dir:
            raise ValueError("frames_every > 0 needs an output_dir")

    def loads_at(self, t, y):
        a = self.sim.aero_state(t, y[:3], y[3:])
        return a, self.heating.evaluate(a, self.body.theta, self.body.surface_temperature(), self.sim.settings.diameter / 2.0,
                                        T_mean=self.body.mean_temperature())

    def row(self, t, y, a, loads):
        body = self.body
        r = self.sim.sample_row(t, y, a)
        r.update({"convective_heat_W": loads.total(body.surface.areas), "rad_cooling_W": -body.radiated_power(),
                  "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat,
                  "q_stag_Wm2": loads.q_stag, "heating_blend_f": loads.blend})
        r.update(body.surface_stats())
        return r

    def run(self):
        started = time.perf_counter()
        sim, body, s = self.sim, self.body, self.settings
        a, loads = self.loads_at(sim.t, sim.y)
        rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
        self._frame(step, sim.t, loads)
        while sim.end_reason is None:
            dt = sim.advance(s.dt)
            a, loads = self.loads_at(sim.t, sim.y)
            body.advance(sim.t, dt, loads)
            step += 1
            rows.append(self.row(sim.t, sim.y, a, loads))
            states.append(sim.y.copy())
            self._frame(step, sim.t, loads)
        self._collection()
        columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
        history = tj.History(columns, np.array(states), sim.end_reason)
        history.results = sim.results(history, sim.nfev, time.perf_counter() - started)
        history.results.update(self.thermal_results(history))
        return history

    def thermal_results(self, history):
        c, body = history.columns, self.body
        i_surf, i_mean, i_q = int(np.argmax(c["surface_T_max_K"])), int(np.argmax(c["temperature_K"])), int(np.argmax(c["convective_heat_W"]))
        return {
            "peak_surface_T_K": float(c["surface_T_max_K"][i_surf]), "time_of_peak_surface_T_s": float(c["time_s"][i_surf]),
            "altitude_of_peak_surface_T_km": float(c["altitude_km"][i_surf]),
            "peak_mean_T_K": float(c["temperature_K"][i_mean]), "time_of_peak_mean_T_s": float(c["time_s"][i_mean]),
            "peak_convective_heat_W": float(c["convective_heat_W"][i_q]), "time_of_peak_heating_s": float(c["time_s"][i_q]),
            "altitude_of_peak_heating_km": float(c["altitude_km"][i_q]),
            "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat, "radiated_heat_J": body.radiated_heat,
            "energy_balance_residual": body.energy_balance_residual(),
            "n_macro_steps": len(body.iterations), "mean_newton_iterations": float(np.mean(body.iterations)) if body.iterations else 0.0,
            "n_frames": len(self.settings.frames),
        }

    def _frame(self, step, t, loads):
        s = self.settings
        if not s.frames_every or step % s.frames_every:
            return
        k = len(s.frames)
        write_vtk_frame(s.output_dir, k, self.body, loads)
        s.frames.append((t, k))

    def _collection(self):
        s = self.settings
        if not s.frames:
            return
        for name in ("field", "surface"):
            entries = "".join('<DataSet timestep="{:.6g}" file="{}_{}.{}"/>\n'.format(t, name, k, "vtu" if name == "field" else "vtp")
                              for t, k in s.frames)
            with open(os.path.join(s.output_dir, name + ".pvd"), "w") as fh:
                fh.write(PVD_TEMPLATE.format(entries))


def write_vtk_frame(output_dir, k, body, loads):
    """field_<k>.vtu: nodal T on the volume mesh; surface_<k>.vtp: q_conv, q_rad, T per patch (PyVista/VTK XML)."""
    import pyvista as pv
    from .thermal import SIGMA_SB
    os.makedirs(output_dir, exist_ok=True)
    mesh, surface, T = body.mesh, body.surface, body.field()
    cells = np.hstack([np.full((mesh.n_elements, 1), 4), mesh.tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(mesh.n_elements, pv.CellType.TETRA), mesh.points)
    grid.point_data["T"] = T
    grid.save(os.path.join(output_dir, "field_{}.vtu".format(k)))
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(mesh.points, faces)
    Tf = surface.facet_mean(T)
    poly.point_data["T"] = T
    poly.cell_data["q_conv"] = loads.q_conv
    poly.cell_data["q_rad"] = body.emissivity * SIGMA_SB * (Tf ** 4 - body.T_ambient ** 4)
    poly.cell_data["T_patch"] = Tf
    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py -q`
Expected: `4 passed` in ~60 s (the full-flight regression with Δt 2 s on the coarse mesh takes most of it).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/coupled.py tests/test_reentry_model_coupled.py
git commit -m "reentry_model: lockstep trajectory-thermal driver with heat bookkeeping and VTK series

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 8: Surface-temperature animation and stills

**Files:**
- Create: `reentry_model/viz.py`
- Test: `tests/test_reentry_model_viz.py`

**Interfaces:**
- Consumes: the `surface.pvd` + `surface_<k>.vtp` series and the run's `History` (Task 7); pyvista off-screen, imageio (+ffmpeg).
- Produces: `read_series(run_dir, name="surface") -> [(t, path)]`, `render_frame(poly, clim, title, path=None, radius=0.05) -> RGB uint8 array`, `frame_title(t, history)`, `animate(run_dir, history, radius, fps=10, animation=True, stills=True) -> {"animation": path|None, "stills": [paths]}` (frames in `<run_dir>/frames/frame_NNNN.png`, `animation.mp4` or `animation.gif`, stills `<run_dir>/stills/{start,peak_heating,peak_surface_T,end}_t<N>s.png`); constants `CAMERA`, `WINDOW`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_viz.py`:

```python
"""viz.py: headless rendering of one frame and the animation/stills from a short coupled run's VTK series."""
import os

import numpy as np
import pytest

from reentry_model import viz

pytest.importorskip("pyvista")


def test_render_one_frame(coarse_sphere_mesh, tmp_path):
    import pyvista as pv
    surface = coarse_sphere_mesh.surface()
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(coarse_sphere_mesh.points, faces)
    poly.point_data["T"] = 300.0 + 500.0 * np.clip(coarse_sphere_mesh.points[:, 0] / 0.05, 0.0, 1.0)
    image = viz.render_frame(poly, (300.0, 800.0), "t = 1.0 s", str(tmp_path / "f.png"))
    assert image.shape == (720, 960, 3) and image.dtype == np.uint8 and os.path.getsize(tmp_path / "f.png") > 1000
    assert image.std() > 10.0                                            # not a blank frame


def test_animation_from_a_short_run(coarse_sphere_mesh, tmp_path):
    from test_reentry_model_coupled import simulator, thermal_body
    from reentry_model import coupled, heating
    b = thermal_body(coarse_sphere_mesh)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=10, output_dir=run_dir)
    hist = coupled.CoupledRun(simulator(b, t_max=15.0), b, heating.SesamEquivalentHeating(), settings).run()
    series = viz.read_series(run_dir)
    assert [t for t, _ in series] == [0.0, 5.0, 10.0, 15.0]
    out = viz.animate(run_dir, hist, 0.05, fps=5)
    assert os.path.isfile(out["animation"]) and os.path.getsize(out["animation"]) > 1000
    assert len(out["stills"]) == 4 and all(os.path.isfile(p) for p in out["stills"])
    assert sorted(os.listdir(os.path.join(run_dir, "frames"))) == ["frame_{:04d}.png".format(k) for k in range(4)]
    assert viz.frame_title(5.0, hist).startswith("t = 5.0 s   h = ")
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_viz.py -q`
Expected: collection error (`reentry_model.viz` missing).

- [ ] **Step 3: Write `reentry_model/viz.py`**

```python
"""Surface-temperature animation and stills from a coupled run's VTK series (spec section 9).

Reads only exported files: `surface.pvd` + `surface_<k>.vtp` (nodal T, per-patch q_conv) written by coupled.py and the
run's history CSV (time, altitude, velocity, surface_T_max_K for the fixed colour scale). PyVista renders off-screen;
frames go to MP4 through imageio-ffmpeg, with a GIF fallback when the MP4 writer is unavailable."""
import os
import xml.etree.ElementTree as ET

import numpy as np

CAMERA = [(0.19, -0.14, 0.11), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]     # three-quarter view of the windward (+x) face, for R = 0.05 m
WINDOW = (960, 720)


def read_series(run_dir, name="surface"):
    """[(time, path)] of a PVD collection written by coupled.py."""
    root = ET.parse(os.path.join(run_dir, name + ".pvd")).getroot()
    return [(float(d.get("timestep")), os.path.join(run_dir, d.get("file"))) for d in root.iter("DataSet")]


def render_frame(poly, clim, title, path=None, radius=0.05):
    """One off-screen frame of the surface coloured by nodal T; returns the RGB array (and writes PNG if `path`)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(poly, scalars="T", cmap="inferno", clim=clim, smooth_shading=True,
                     scalar_bar_args={"title": "surface T [K]", "fmt": "%.0f", "vertical": True, "position_x": 0.86, "position_y": 0.15,
                                      "width": 0.05, "height": 0.6, "title_font_size": 14, "label_font_size": 12})
    plotter.add_text(title, position="upper_left", font_size=11)
    scale = radius / 0.05
    plotter.camera_position = [tuple(scale * c for c in CAMERA[0]), CAMERA[1], CAMERA[2]]
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def frame_title(t, history):
    c = history.columns
    return "t = {:.1f} s   h = {:.1f} km   V = {:.2f} km/s".format(
        t, np.interp(t, c["time_s"], c["altitude_km"]), np.interp(t, c["time_s"], c["velocity_kms"]))


def animate(run_dir, history, radius, fps=10, animation=True, stills=True):
    """MP4 (GIF fallback) of the surface temperature over the run, plus stills at the start, peak heating, peak surface
    temperature and the end (with animation=False only the stills' frames are rendered). Returns {"animation", "stills"}."""
    import imageio.v2 as imageio
    import pyvista as pv
    series = read_series(run_dir)
    c = history.columns
    clim = (float(c["temperature_K"][0]), float(c["surface_T_max_K"].max()))
    times = np.array([t for t, _ in series])
    marks = {"start": 0.0, "peak_heating": float(c["time_s"][int(np.argmax(c["convective_heat_W"]))]),
             "peak_surface_T": float(c["time_s"][int(np.argmax(c["surface_T_max_K"]))]), "end": float(c["time_s"][-1])}
    wanted = {int(np.argmin(np.abs(times - t))) for t in marks.values()} if stills else set()
    frames, out = {}, {"animation": None, "stills": []}
    frames_dir = os.path.join(run_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    for k, (t, path) in enumerate(series):
        if animation or k in wanted:
            frames[k] = render_frame(pv.read(path), clim, frame_title(t, history), os.path.join(frames_dir, "frame_{:04d}.png".format(k)), radius)
    if animation and frames:
        try:
            writer = imageio.get_writer(os.path.join(run_dir, "animation.mp4"), fps=fps, codec="libx264", quality=7, macro_block_size=8)
            for k in sorted(frames):
                writer.append_data(frames[k])
            writer.close()
            out["animation"] = os.path.join(run_dir, "animation.mp4")
        except Exception:                                   # no ffmpeg: GIF
            imageio.mimsave(os.path.join(run_dir, "animation.gif"), [frames[k] for k in sorted(frames)], duration=1.0 / fps)
            out["animation"] = os.path.join(run_dir, "animation.gif")
    if stills and series:
        stills_dir = os.path.join(run_dir, "stills")
        os.makedirs(stills_dir, exist_ok=True)
        for label, t in marks.items():
            k = int(np.argmin(np.abs(times - t)))
            path = os.path.join(stills_dir, "{}_t{:.0f}s.png".format(label, times[k]))
            imageio.imwrite(path, frames[k])
            out["stills"].append(path)
    return out
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_viz.py -q`
Expected: `2 passed` (~10 s; the first VTK frame takes ~4 s to initialise). Open one frame (`tests` write them under pytest's tmp dir, or render one by hand) and check: sphere in a three-quarter view with the windward (+x) face to the right, `inferno` colours, a vertical colour bar labelled "surface T [K]" with integer ticks, the title top-left.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/viz.py tests/test_reentry_model_viz.py
git commit -m "reentry_model: off-screen surface-temperature frames, MP4/GIF animation and stills

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 9: Heat and temperature comparison against SESAM

**Files:**
- Modify: `reentry_model/sesam_io.py` (replace the file)
- Modify: `reentry_model/compare.py` (replace the file)
- Modify: `tests/test_reentry_model_compare.py` (append)

**Interfaces:**
- Consumes: `History` with `THERMAL_COLUMNS` (Task 7), `Reference` heat columns.
- Produces: `Reference.convective_heat`, `.rad_cooling` (negative), `.integrated_heat` (each `None` when the CSV lacks the column); `compare.THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")`, `CONTINUUM_KN = 0.01`, `has_thermal(history, reference=None)`, `align_thermal(history, reference)`, `thermal_metrics(history, reference)` (keys `n_points`, `n_continuum_points`, `Q_conv{max, rms, peak_reference, continuum_rel_max, continuum_rel_rms}`, `integrated_heat{rel_error_end_of_hypersonic, rel_error_end, model_end_J, reference_end_J, ratio_end}`, `temperature{dT_max_K, dT_rel_max, dT_rms_K, model_peak_K, reference_peak_K}`, `radiated{max, rms, peak_reference}`), `plot_thermal(history, reference, outdir, title) -> paths`. Metric definitions per amendment 10.

- [ ] **Step 1: Append the tests**

Append to `tests/test_reentry_model_compare.py`:

```python
def _thermal_history_like(reference, factor=1.0, dT=0.0):
    """A model history built from the reference's own heat columns (scaled), at half its time stamps."""
    from reentry_model import trajectory as tj
    t = reference.time[::2]
    pick = lambda arr: np.interp(t, reference.time, arr)
    columns = {"time_s": t, "velocity_kms": pick(reference.velocity) / 1e3, "altitude_km": pick(reference.altitude) / 1e3,
               "knudsen": pick(reference.knudsen), "convective_heat_W": factor * pick(reference.convective_heat),
               "integrated_heat_J": factor * pick(reference.integrated_heat), "rad_cooling_W": pick(reference.rad_cooling),
               "temperature_K": pick(reference.temperature) + dT, "surface_T_max_K": pick(reference.temperature) + dT + 50.0,
               "surface_T_min_K": pick(reference.temperature) + dT - 20.0}
    return tj.History(columns, np.zeros((t.size, 6)), "ground")


def test_reference_carries_the_heat_columns(ref):
    assert ref.convective_heat[0] == pytest.approx(15819.606) and ref.rad_cooling[0] == pytest.approx(-5.772)
    assert ref.integrated_heat[0] == 0.0 and ref.integrated_heat[-1] > 1e5


def test_thermal_metrics_of_a_scaled_copy_and_plots(ref, tmp_path):
    hist = _thermal_history_like(ref, factor=1.02, dT=3.0)
    assert compare.has_thermal(hist, ref) and not compare.has_thermal(synthetic_history(ref), ref)
    m = compare.thermal_metrics(hist, ref)
    assert m["Q_conv"]["max"] == pytest.approx(0.02, abs=2e-3) and 0.015 < m["Q_conv"]["continuum_rel_max"] < 0.04   # + interpolation at half the stamps
    assert m["integrated_heat"]["rel_error_end"] == pytest.approx(0.02, abs=1e-3) and m["integrated_heat"]["ratio_end"] == pytest.approx(1.02, abs=1e-3)
    assert m["temperature"]["dT_max_K"] == pytest.approx(3.0, abs=0.5) and m["radiated"]["max"] < 0.05
    assert m["n_continuum_points"] > 100 and m["integrated_heat"]["rel_error_end_of_hypersonic"] == pytest.approx(0.02, abs=2e-3)
    paths = compare.plot_thermal(hist, ref, str(tmp_path), "scaled copy")
    assert [os.path.basename(p) for p in paths] == list(compare.THERMAL_PLOT_NAMES) and all(os.path.getsize(p) > 5000 for p in paths)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py -q`
Expected: the two new tests fail (`Reference` has no `convective_heat`; `compare.has_thermal` missing).

- [ ] **Step 3: Replace `reentry_model/sesam_io.py`**

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
    convective_heat: np.ndarray = None      # W (SESAM's convective_heat_W), None if the CSV lacks the column
    rad_cooling: np.ndarray = None          # W, negative (SESAM's rad_cooling_W)
    integrated_heat: np.ndarray = None      # J, integral of the convective heat (SESAM's integrated_heat_J)


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
    optional = lambda key: col(key) if key in rows[0] else None
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
        convective_heat=optional("convective_heat_W"), rad_cooling=optional("rad_cooling_W"), integrated_heat=optional("integrated_heat_J"),
    )
```

- [ ] **Step 4: Replace `reentry_model/compare.py`**

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

MODEL_COLOR = "#3987e5"
REF_COLOR = INK
HYPERSONIC_MS = 1000.0
PLOT_NAMES = ("velocity_time.png", "altitude_time.png", "altitude_velocity.png", "angles_time.png",
              "ground_track.png", "regime_drag.png")
THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")
CONTINUUM_KN = 0.01


def align(history, reference):
    t_mod = history.columns["time_s"]
    # Exclude reference points outside the model's time window to avoid extrapolation
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
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


def _overlay(ax, x_ref, y_ref, x_model, y_model, ylabel, xlabel=None, title=None, legend=True, log_y=False):
    """Draw reference and model overlay lines with consistent styling."""
    plot_fn = ax.semilogy if log_y else ax.plot
    plot_fn(x_ref, y_ref, color=REF_COLOR, lw=1.6, label="SESAM")
    plot_fn(x_model, y_model, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title, color=SECOND, fontsize=10)
    if legend:
        ax.legend(frameon=False)
    strip_top_right_spines(ax)


def metrics(history, reference):
    a = align(history, reference)
    t_model_end, t_ref_end = float(history.columns["time_s"][-1]), float(reference.time[-1])
    v_model_final = float(history.columns["velocity_kms"][-1] * 1e3)
    v_ref_final = float(reference.velocity[-1])
    return {
        "n_points": int(a["t"].size),
        "model_end_time_s": t_model_end, "reference_end_time_s": t_ref_end,
        "d_end_time_s": t_model_end - t_ref_end, "d_end_time_rel": (t_model_end - t_ref_end) / t_ref_end,
        "final_velocity_model_ms": v_model_final,
        "final_velocity_reference_ms": v_ref_final,
        "d_final_velocity_ms": v_model_final - v_ref_final,
        "all": _phase(a, np.ones(a["t"].size, dtype=bool)),
        "hypersonic": _phase(a, a["V_ref"] > HYPERSONIC_MS),
    }


def _overlay_with_residual(a, key, ylabel, resid_label, scale, path, title):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    _overlay(ax, a["t"], a[key + "_ref"] * scale, a["t"], a[key + "_model"] * scale, ylabel, title=title)
    rx.plot(a["t"], (a[key + "_model"] - a[key + "_ref"]), color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel(resid_label); rx.set_xlabel("time [s]")
    strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_all(history, reference, outdir, title):
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    a = align(history, reference)
    paths = [os.path.join(outdir, n) for n in PLOT_NAMES]
    _overlay_with_residual(a, "V", "velocity [km/s]", "model - SESAM [m/s]", 1e-3, paths[0], title)
    _overlay_with_residual(a, "h", "altitude [km]", "model - SESAM [m]", 1e-3, paths[1], title)

    fig, ax = plt.subplots(figsize=(7, 5))
    _overlay(ax, a["V_ref"] / 1e3, a["h_ref"] / 1e3, a["V_model"] / 1e3, a["h_model"] / 1e3,
             ylabel="altitude [km]", xlabel="velocity [km/s]", title=title)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)

    fig, (g, hd) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    _overlay(g, a["t"], a["gamma_ref"], a["t"], a["gamma_model"], ylabel="flight-path angle [deg]", title=title)
    _overlay(hd, a["t"], a["heading_ref"], a["t"], a["heading_model"], ylabel="heading [deg]", xlabel="time [s]", legend=False)
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    _overlay(ax, a["lon_ref"], a["lat_ref"], a["lon_model"], a["lat_model"],
             ylabel="latitude [deg]", xlabel="longitude [deg]", title=title)
    fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)

    fig, (kx, cx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    _overlay(kx, a["t"], np.maximum(a["kn_ref"], 1e-7), a["t"], np.maximum(a["kn_model"], 1e-7),
             ylabel="Knudsen number", title=title, log_y=True)
    for thr in (10.0, 1.0, 0.1, 0.01):
        kx.axhline(thr, color=MUTED, lw=0.5, ls=":")
    _overlay(cx, a["t"], a["cd_ref"], a["t"], a["cd_model"],
             ylabel="C_D", xlabel="time [s]", legend=False)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    return paths


def has_thermal(history, reference=None):
    """True when the model history carries the coupled run's heat columns (and the reference has SESAM's)."""
    ok = "convective_heat_W" in history.columns
    return ok and (reference is None or reference.convective_heat is not None)


def align_thermal(history, reference):
    """Model heat/temperature columns interpolated at the reference's time stamps (within the model's window)."""
    t_mod = history.columns["time_s"]
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
    t = reference.time[mask]
    at = lambda col: np.interp(t, t_mod, history.columns[col])
    return {
        "t": t, "V_ref": reference.velocity[mask], "kn_ref": reference.knudsen[mask],
        "Q_model": at("convective_heat_W"), "Q_ref": reference.convective_heat[mask],
        "H_model": at("integrated_heat_J"), "H_ref": reference.integrated_heat[mask],
        "R_model": -at("rad_cooling_W"), "R_ref": -reference.rad_cooling[mask],
        "T_model": at("temperature_K"), "T_ref": reference.temperature[mask],
        "Tmax_model": at("surface_T_max_K"), "Tmin_model": at("surface_T_min_K"),
    }


def _peak_relative(delta, ref, mask):
    """Errors normalised by the reference's peak over the mask (robust where the reference is small)."""
    peak = float(np.abs(ref[mask]).max()) if mask.any() else float("nan")
    d = np.abs(delta[mask]) / peak if mask.any() else np.array([])
    return {"max": float(d.max()) if d.size else float("nan"), "rms": float(math.sqrt(np.mean(d * d))) if d.size else float("nan"),
            "peak_reference": peak}


def thermal_metrics(history, reference):
    """Spec section 9 metrics. Convective and radiated power: max/rms error relative to SESAM's peak over the hypersonic
    phase (SESAM's free-molecular heating is ~13x below the textbook value, facts note s.4, so a point-wise relative
    error would be dominated by the high-Kn start where the absolute heat is negligible); point-wise relative error
    of Q_conv restricted to the continuum part (Kn_ref < 0.01 and Q_ref above 10 % of its peak); integrated heat at the end of the hypersonic phase and
    at the end; energy-equivalent temperature vs SESAM's lumped temperature over the whole flight."""
    a = align_thermal(history, reference)
    hyp = a["V_ref"] > HYPERSONIC_MS
    peak = float(np.abs(a["Q_ref"][hyp]).max()) if hyp.any() else 0.0
    # point-wise only where SESAM's heat is at least 10 % of its peak: its hot-wall factor clamps the heating to zero
    # late in the flight, and relative errors next to the clamp are unbounded by construction
    cont = hyp & (a["kn_ref"] < CONTINUUM_KN) & (a["Q_ref"] > 0.1 * peak)
    dQ = a["Q_model"] - a["Q_ref"]
    rel_cont = np.abs(dQ[cont]) / a["Q_ref"][cont]
    i_hyp = int(np.nonzero(hyp)[0][-1]) if hyp.any() else len(a["t"]) - 1
    dT = a["T_model"] - a["T_ref"]
    return {
        "n_points": int(a["t"].size), "n_continuum_points": int(cont.sum()),
        "Q_conv": {**_peak_relative(dQ, a["Q_ref"], hyp), "continuum_rel_max": float(rel_cont.max()) if rel_cont.size else float("nan"),
                   "continuum_rel_rms": float(math.sqrt(np.mean(rel_cont * rel_cont))) if rel_cont.size else float("nan")},
        "integrated_heat": {"rel_error_end_of_hypersonic": float(a["H_model"][i_hyp] / a["H_ref"][i_hyp] - 1.0) if a["H_ref"][i_hyp] else float("nan"),
                            "rel_error_end": float(a["H_model"][-1] / a["H_ref"][-1] - 1.0) if a["H_ref"][-1] else float("nan"),
                            "model_end_J": float(a["H_model"][-1]), "reference_end_J": float(a["H_ref"][-1]),
                            "ratio_end": float(a["H_model"][-1] / a["H_ref"][-1]) if a["H_ref"][-1] else float("nan")},
        "temperature": {"dT_max_K": float(np.abs(dT).max()), "dT_rel_max": float((np.abs(dT) / a["T_ref"]).max()),
                        "dT_rms_K": float(math.sqrt(np.mean(dT * dT))),
                        "model_peak_K": float(a["T_model"].max()), "reference_peak_K": float(a["T_ref"].max())},
        "radiated": _peak_relative(a["R_model"] - a["R_ref"], a["R_ref"], hyp),
    }


def plot_thermal(history, reference, outdir, title):
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    a = align_thermal(history, reference)
    paths = [os.path.join(outdir, n) for n in THERMAL_PLOT_NAMES]
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    _overlay(ax, a["t"], a["Q_ref"] / 1e3, a["t"], a["Q_model"] / 1e3, "convective power [kW]", title=title)
    rx.plot(a["t"], (a["Q_model"] - a["Q_ref"]) / 1e3, color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel("model - SESAM [kW]"); rx.set_xlabel("time [s]")
    strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(paths[0], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    _overlay(ax, a["t"], a["T_ref"], a["t"], a["T_model"], "temperature [K]", xlabel="time [s]", title=title)
    ax.plot(a["t"], a["Tmax_model"], color=MODEL_COLOR, lw=0.8, ls=":", label="model surface max")
    ax.plot(a["t"], a["Tmin_model"], color=MUTED, lw=0.8, ls=":", label="model surface min")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(paths[1], dpi=150); plt.close(fig)

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    _overlay(ax, a["t"], a["H_ref"] / 1e3, a["t"], a["H_model"] / 1e3, "integrated convective heat [kJ]", title=title)
    with np.errstate(divide="ignore", invalid="ignore"):
        rx.plot(a["t"], 100.0 * (a["H_model"] / a["H_ref"] - 1.0), color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel("model / SESAM - 1 [%]"); rx.set_xlabel("time [s]"); rx.set_ylim(-20.0, 20.0)
    strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)
    return paths
```

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_sesam_io.py -q`
Expected: all pass (`6 passed` in compare).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/sesam_io.py reentry_model/compare.py tests/test_reentry_model_compare.py
git commit -m "reentry_model: SESAM heat columns on Reference, thermal metrics and the three heat/temperature plots

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 10: Command line for the coupled model

**Files:**
- Modify: `reentry_model/cli.py` (replace the file)
- Modify: `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: everything above (`coupled`, `heating`, `material`, `mesh`, `thermal`, `viz`, `compare`).
- Produces: `run` flags `--thermal none|fem` (default none), `--heating physics|sesam` (default physics), `--stagnation`, `--bridging-heat`, `--matting-n`, `--accommodation`, `--catalycity`, `--material`, `--emissivity`, `--t-ambient`, `--mesh-size`, `--h-surface` (mm), `--h-core` (mm), `--thermal-solver skfem|fenicsx`, `--linear-solver direct|amg`, `--lumped-mass`, `--dt`, `--frames-every`, `--animate`, `--stills`; run name suffix `_fem-<heating>`; `build_thermal(args, settings, mass) -> (ThermalBody, HeatingModel, settings dict)`; JSON `settings` gains `thermal`, `heating`, `material`, `material_file`, `emissivity`, `t_ambient_K`, `mesh_file`, `h_surface_mm`, `h_core_mm`, `n_nodes`, `n_elements`, `n_patches`, `thermal_solver`, `linear_solver`, `lumped_mass`, `macro_step_s`, `frames_every` (+ `stagnation`, `bridging_heat`, `matting_n`, `accommodation`, `catalycity` in physics mode); `comparison.thermal_metrics` and the three extra plots when a reference is given; `files.vtk_dir`, `files.animation`, `files.stills`; `provenance` lists the versions of `skfem`, `gmsh`, `pyamg`, `pyvista`, `cantera`, `dolfinx`; exit 2 on `MissingBackend` and on a missing optional library; `compare` reports thermal metrics/plots when both files carry the columns. The CLI help of `--heating` says the sesam mode is a verification device, not physical.

- [ ] **Step 1: Append the tests**

Append to `tests/test_reentry_model_cli.py`:

```python
US76_100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind.csv")
FEM = BASE + ["--atmosphere", "us76", "--thermal", "fem", "--heating", "sesam", "--h-surface", "4", "--h-core", "20", "--quiet"]


def test_run_name_with_heating():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "sesam").endswith("_none_fem-sesam")


def test_thermal_run_writes_columns_plots_and_json(tmp_path):
    rc = cli.main(FEM + ["--t-max", "20", "--dt", "0.5", "--reference", US76_100, "--outdir", str(tmp_path), "--name", "fem_short",
                         "--frames-every", "10", "--stills"])
    assert rc == 0
    rows = list(csv.DictReader(open(tmp_path / "fem_short.csv")))
    assert len(rows) == 41 and float(rows[-1]["time_s"]) == 20.0 and "surface_T_max_K" in rows[0] and "convective_heat_W" in rows[0]
    assert float(rows[-1]["temperature_K"]) > 300.0 and float(rows[0]["convective_heat_W"]) > 1e4
    doc = json.load(open(tmp_path / "fem_short.json"))
    s = doc["settings"]
    assert s["thermal"] == "fem" and s["heating"] == "sesam" and s["n_nodes"] > 2000 and s["macro_step_s"] == 0.5 and s["frames_every"] == 10
    assert s["material"] == "AA7075_nomelt" and s["emissivity"] == 0.4 and s["thermal_solver"] == "skfem"
    assert doc["results"]["n_macro_steps"] == 40 and abs(doc["results"]["energy_balance_residual"]) < 1e-6
    assert "thermal_metrics" in doc["comparison"] and doc["comparison"]["thermal_metrics"]["Q_conv"]["max"] < 0.3
    for plot in compare.PLOT_NAMES + compare.THERMAL_PLOT_NAMES:
        assert os.path.isfile(tmp_path / "fem_short" / plot)
    assert os.path.isfile(tmp_path / "fem_short" / "vtk" / "field.pvd") and doc["files"]["vtk_dir"].endswith("vtk")
    assert doc["files"]["animation"] is None and len(doc["files"]["stills"]) == 4 and all(os.path.isfile(p) for p in doc["files"]["stills"])
    assert doc["provenance"]["skfem"] and doc["provenance"]["gmsh"]


def test_thermal_none_is_step_one(tmp_path):
    rc = cli.main(BASE + ["--atmosphere", "us76", "--t-max", "10", "--cadence", "5", "--outdir", str(tmp_path), "--name", "plain", "--quiet"])
    assert rc == 0
    rows = list(csv.DictReader(open(tmp_path / "plain.csv")))
    assert "convective_heat_W" not in rows[0] and json.load(open(tmp_path / "plain.json"))["settings"]["thermal"] == "none"
    with pytest.raises(SystemExit) as exc:                                                    # parser.error: needs --thermal fem
        cli.main(BASE + ["--atmosphere", "us76", "--animate", "--outdir", str(tmp_path), "--quiet"])
    assert exc.value.code == 2


def test_missing_fenicsx_exits_2(tmp_path, capsys):
    try:
        import dolfinx  # noqa: F401
        pytest.skip("dolfinx is importable here")
    except ImportError:
        pass
    rc = cli.main(FEM + ["--thermal-solver", "fenicsx", "--t-max", "5", "--outdir", str(tmp_path)])
    assert rc == 2 and "fenicsx_env" in capsys.readouterr().err


def test_compare_subcommand_with_thermal_columns(tmp_path):
    cli.main(FEM + ["--t-max", "10", "--outdir", str(tmp_path), "--name", "m"])
    rc = cli.main(["compare", "--model", str(tmp_path / "m.csv"), "--reference", US76_100, "--outdir", str(tmp_path / "cmp"), "--quiet"])
    assert rc == 0
    doc = json.load(open(tmp_path / "cmp" / "m_vs_reference.json"))
    assert "thermal_metrics" in doc and os.path.isfile(tmp_path / "cmp" / "heating_time.png")
```

- [ ] **Step 2: Run them to verify they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: the five new tests fail (`model_run_name` has no `heating_name` argument, `--thermal` unknown).

- [ ] **Step 3: Replace `reentry_model/cli.py`**

```python
"""Command line of the re-entry model.

    python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
        [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
        [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
    python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]

Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
(cantera for --heating physics, dolfinx for --thermal-solver fenicsx: create the fenicsx_env environment).
"""
import argparse
import math
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from . import __version__, aero, atmosphere, body, compare, coupled, fap, heating, material, mesh, sesam_io, thermal, viz
from . import trajectory as tj
from .earth import GRAVITY_MODELS
from .thermal import MissingBackend

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
THERMAL_MODES = ("none", "fem")


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
    versions = {}
    for name in ("skfem", "gmsh", "pyamg", "pyvista", "cantera", "dolfinx"):
        try:
            versions[name] = __import__(name).__version__
        except Exception:
            versions[name] = None
    return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0], **versions}


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
    r.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")
    r.add_argument("--quiet", action="store_true")
    th = r.add_argument_group("thermal model (Step 2)")
    th.add_argument("--thermal", choices=THERMAL_MODES, default="none", help="none: Step 1 trajectory only (default); fem: coupled 3D conduction")
    th.add_argument("--heating", choices=heating.HEATING_NAMES, default="physics",
                    help="physics (default) or sesam: SESAM's uniform 0.27471 x q_stag on every patch -- a verification device, not physical")
    th.add_argument("--stagnation", choices=heating.STAGNATION_NAMES, default="fay-riddell")
    th.add_argument("--bridging-heat", choices=heating.BRIDGING_HEAT_NAMES, default="matting")
    th.add_argument("--matting-n", type=float, default=1.0, help="Matting exponent n (default %(default)s)")
    th.add_argument("--accommodation", type=float, default=0.8, help="free-molecular energy accommodation A_cq (default %(default)s)")
    th.add_argument("--catalycity", type=float, default=1.0, help="wall catalycity 0..1 in Fay-Riddell (default %(default)s)")
    th.add_argument("--material", default=material.DEFAULT_MATERIAL, help="DRAMA material JSON (default: packaged AA7075_nomelt)")
    th.add_argument("--emissivity", type=float, default=None, help="override the material's emissivity")
    th.add_argument("--t-ambient", type=float, default=0.0, help="radiation background [K] (default %(default)s, SESAM's)")
    th.add_argument("--mesh-size", type=float, default=1.0, help="multiplies --h-surface and --h-core (default %(default)s)")
    th.add_argument("--h-surface", type=float, default=mesh.DEFAULT_H_SURFACE * 1e3, help="surface element size [mm] (default %(default)s)")
    th.add_argument("--h-core", type=float, default=mesh.DEFAULT_H_CORE * 1e3, help="core element size [mm] (default %(default)s)")
    th.add_argument("--thermal-solver", choices=thermal.SOLVER_NAMES, default="skfem")
    th.add_argument("--linear-solver", choices=("direct", "amg"), default="amg")
    th.add_argument("--lumped-mass", action="store_true")
    th.add_argument("--dt", type=float, default=0.5, help="macro step [s] (default %(default)s)")
    th.add_argument("--frames-every", type=int, default=0, help="VTK frame every n macro steps (default 0: none; 10 with --animate/--stills)")
    th.add_argument("--animate", action="store_true", help="MP4/GIF of the surface temperature plus stills")
    th.add_argument("--stills", action="store_true", help="only the four stills (start, peak heating, peak surface T, end)")

    c = sub.add_parser("compare", help="metrics and plots for an existing model history")
    c.add_argument("--model", required=True, help="model history CSV")
    c.add_argument("--reference", required=True, help="SESAM run CSV")
    c.add_argument("--outdir", default=DEFAULT_OUTDIR)
    c.add_argument("--title", default=None)
    c.add_argument("--quiet", action="store_true")
    return p


def build_thermal(args, settings, mass):
    """(ThermalBody, HeatingModel, settings-provenance dict) for --thermal fem."""
    radius = settings.diameter / 2.0
    h_surface, h_core = args.h_surface * 1e-3 * args.mesh_size, args.h_core * 1e-3 * args.mesh_size
    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core)
    mat = material.Material.from_drama_json(args.material)
    solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=args.lumped_mass)
    the_body = body.ThermalBody(the_mesh, mat, solver, mass, T0=args.temperature, emissivity=args.emissivity, T_ambient=args.t_ambient)
    if args.heating == "sesam":
        heating_model = heating.SesamEquivalentHeating()
    else:
        heating_model = heating.PhysicsHeating(stagnation=args.stagnation, bridging=args.bridging_heat, matting_n=args.matting_n,
                                               accommodation=args.accommodation, catalycity=args.catalycity)
    info = {"heating": args.heating, "material": mat.name, "material_file": os.path.abspath(args.material), "emissivity": the_body.emissivity,
            "t_ambient_K": args.t_ambient, "mesh_file": the_mesh.params["path"], "h_surface_mm": h_surface * 1e3, "h_core_mm": h_core * 1e3,
            "n_nodes": the_mesh.n_nodes, "n_elements": the_mesh.n_elements, "n_patches": the_body.surface.n_patches,
            "thermal_solver": args.thermal_solver, "linear_solver": args.linear_solver, "lumped_mass": args.lumped_mass}
    if args.heating == "physics":
        info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                     "accommodation": args.accommodation, "catalycity": args.catalycity})
    return the_body, heating_model, info


def cmd_run(args, parser):
    for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                         ("--cadence", args.cadence), ("--t-max", args.t_max)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.atmosphere not in ATMOSPHERES and not args.atmosphere.startswith("replay:"):
        parser.error("--atmosphere must be one of {} or replay:<sesam.csv>".format(ATMOSPHERES))
    for label, value in (("--dt", args.dt), ("--mesh-size", args.mesh_size), ("--h-surface", args.h_surface), ("--h-core", args.h_core)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if not 0.0 <= args.catalycity <= 1.0:
        parser.error("--catalycity must be within [0, 1]")
    if args.thermal == "none" and (args.animate or args.stills or args.frames_every):
        parser.error("--animate/--stills/--frames-every need --thermal fem")

    initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                              flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
                              lat=math.radians(args.lat), lon=math.radians(args.lon), epoch=args.epoch)
    settings = tj.Settings(diameter=args.diameter * 1e-3, gravity=args.gravity, rtol=args.rtol,
                           cadence=args.cadence, t_max=args.t_max)
    mass = body.sphere_mass(settings.diameter, args.material_density)
    atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
    reference = sesam_io.load_reference(args.reference) if args.reference else None
    name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                       args.heating if args.thermal == "fem" else None)
    run_dir = os.path.join(args.outdir, name)
    os.makedirs(args.outdir, exist_ok=True)
    thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
    if args.thermal == "fem":
        the_body, heating_model, thermal_info = build_thermal(args, settings, mass)
    sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
    if args.thermal == "fem":
        frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
        run = coupled.CoupledRun(sim, the_body, heating_model,
                                 coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk")))
        history = run.run()
        thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
    else:
        history = sim.run(extra_times=reference.time if reference is not None else None)
    csv_path = os.path.join(args.outdir, name + ".csv")
    json_path = os.path.join(args.outdir, name + ".json")
    tj.write_history_csv(history, csv_path)
    prov = provenance()
    # spec section 8: the reference file(s)' SHA-256 also live under provenance, alongside the copies
    # already recorded under "comparison" (--reference) and "settings" (replay_reference/replay_sha256).
    prov["reference_sha256"] = reference.sha256 if reference is not None else None
    prov["replay_sha256"] = atm_info.get("replay_sha256")
    doc = {
        "schema_version": 1,
        "run_name": name,
        "inputs": {"diameter_mm": args.diameter, "initial_velocity_kms": args.velocity, "initial_altitude_km": args.altitude,
                   "flight_path_angle_deg": args.flight_path_angle, "heading_deg": args.heading, "latitude_deg": args.lat,
                   "longitude_deg": args.lon, "epoch_utc": args.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
                   "material_density_kgm3": args.material_density, "mass_kg": mass, "initial_temperature_K": args.temperature},
        "settings": {"atmosphere": atm_name, "wind": args.wind, "bridging": args.bridging, "gravity": args.gravity,
                     "rtol": args.rtol, "atol_position_m": settings.atol_position, "atol_velocity_ms": settings.atol_velocity,
                     "cadence_s": args.cadence, "t_max_s": args.t_max, **atm_info, "thermal": args.thermal, **thermal_info},
        "results": history.results,
        "comparison": None,
        "provenance": prov,
        "files": {"csv": os.path.abspath(csv_path)},
    }
    if reference is not None:
        plots = compare.plot_all(history, reference, run_dir, name)
        doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
                             "reference_sha256": reference.sha256, "metrics": compare.metrics(history, reference),
                             "plots": [os.path.abspath(p) for p in plots]}
        if compare.has_thermal(history, reference):
            doc["comparison"]["thermal_metrics"] = compare.thermal_metrics(history, reference)
            doc["comparison"]["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, run_dir, name)]
    if args.thermal == "fem" and (args.animate or args.stills):
        out = viz.animate(os.path.join(run_dir, "vtk"), history, settings.diameter / 2.0, animation=args.animate)
        doc["files"]["animation"], doc["files"]["stills"] = out.get("animation"), out["stills"]
    if args.thermal == "fem":
        doc["files"]["vtk_dir"] = os.path.abspath(os.path.join(run_dir, "vtk")) if history.results.get("n_frames") else None
    tj.write_run_json(json_path, doc)
    if not args.quiet:
        res = history.results
        print("{}: {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
            name, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
            history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
        if args.thermal == "fem":
            print("  thermal: peak surface T {:.0f} K at t = {:.0f} s, peak mean T {:.0f} K, integrated heat {:.3g} J, "
                  "energy balance residual {:.1e}, {} macro steps, {:.1f} Newton iterations/step".format(
                      res["peak_surface_T_K"], res["time_of_peak_surface_T_s"], res["peak_mean_T_K"], res["integrated_heat_J"],
                      res["energy_balance_residual"], res["n_macro_steps"], res["mean_newton_iterations"]))
        if reference is not None:
            hyp = doc["comparison"]["metrics"]["hypersonic"]
            print("  vs {}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m; end time {:+.1f} s".format(
                reference.name, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], doc["comparison"]["metrics"]["d_end_time_s"]))
            if "thermal_metrics" in doc["comparison"]:
                tm = doc["comparison"]["thermal_metrics"]
                print("  heat vs SESAM: Q_conv max {:.2%} of peak ({:.2%} point-wise, continuum), integrated heat {:+.2%} (hypersonic) "
                      "{:+.2%} (end), |dT_eq| max {:.1f} K ({:.2%}), radiated max {:.2%} of peak".format(
                          tm["Q_conv"]["max"], tm["Q_conv"]["continuum_rel_max"], tm["integrated_heat"]["rel_error_end_of_hypersonic"],
                          tm["integrated_heat"]["rel_error_end"], tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"], tm["radiated"]["max"]))
        print("  csv  -> {}\n  json -> {}".format(os.path.abspath(csv_path), os.path.abspath(json_path)))
    # spec section 7: exit 1 when the flight escaped rather than reaching the ground or t_max;
    # the CSV/JSON are already written above so the escaped trajectory is still available.
    if history.end_reason == "escape":
        return 1
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
    if compare.has_thermal(history, reference):
        doc["thermal_metrics"] = compare.thermal_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, args.outdir, title)]
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
    except KeyError as exc:
        # e.g. fap.solar_indices(): the run epoch (or the day before it) has no record in
        # data/fap_day.dat; could also be a reference CSV missing an expected column.
        key = exc.args[0] if exc.args else exc
        print("ERROR: missing key {!r} (epoch outside the fap file, or a reference CSV without that column)".format(key),
              file=sys.stderr)
        return 2
    except (ValueError, FileNotFoundError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except MissingBackend as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except ModuleNotFoundError as exc:
        print("ERROR: missing optional library {!r}: install requirements-step2.txt into drama_env (dolfinx: the separate "
              "fenicsx_env environment)".format(exc.name), file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
```

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: `17 passed` (~15 s). `test_missing_fenicsx_exits_2` exercises the stub from Task 3 (dolfinx absent in drama_env).

- [ ] **Step 5: Run the whole unit tier**

Run: `"$PY" -m pytest -m "not drama and not reference" -q`
Expected: everything passes (~4 min: Step 1's 278 plus ~120 new tests).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/cli.py tests/test_reentry_model_cli.py
git commit -m "reentry_model: --thermal fem command line with heating, mesh, solver, animation options

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 11: FEniCSx backend (behind the `--thermal-solver fenicsx` switch)

**Files:**
- Modify: `reentry_model/thermal/fenicsx_backend.py` (replace the stub)
- Test: `tests/test_reentry_model_fenicsx.py`

**Interfaces:**
- Consumes: the `ThermalSolver` protocol (Task 3), `mesh.VolumeMesh`, `material.Material`; dolfinx ≥ 0.9, ufl, basix, petsc4py, mpi4py (all lazy).
- Produces: `FenicsxThermalSolver(linear_solver="amg"|"direct", lumped_mass=False, newton_tol=1e-6, max_iterations=30, cg_tol=1e-10)` implementing `setup`, `set_temperature`, `step` (with `dirichlet`), `temperature`, `energy`, `radiated_power`, plus `faces`, `areas`, `cell_volumes`, `cell_dofs`, `node_of_dof`, `dof_of_node`, `facet_load(q)`; raises `MissingBackend` without dolfinx and `ValueError` for `lumped_mass=True`.
- **This backend cannot be executed in `drama_env`** (no dolfinx) and the user has not yet asked for `fenicsx_env`. Write it exactly as below (dolfinx 0.9/0.10 API: `dolfinx.mesh.create_mesh(comm, cells, x, ufl.Mesh(basix.ufl.element(...)))`, `fem.functionspace`, `fem.form`, `dolfinx.fem.petsc.assemble_matrix/assemble_vector/apply_lifting/set_bc`, `Function.x.petsc_vec`), keep the tests skipping (`pytest.importorskip("dolfinx")`), and report in the task report that the backend is untested until `fenicsx_env` exists. When the user creates `fenicsx_env` (`conda create -n fenicsx_env -c conda-forge fenics-dolfinx mpich gmsh python=3.12` plus `pip install -r requirements-step2.txt` there), the tests run with that interpreter: `/path/to/fenicsx_env/bin/python -m pytest tests/test_reentry_model_fenicsx.py -q`.

- [ ] **Step 1: Write the tests** (they skip in drama_env)

Create `tests/test_reentry_model_fenicsx.py`:

```python
"""thermal/fenicsx_backend.py: the same analytic checks as the scikit-fem backend and a cross-check against it.
Skipped unless dolfinx is importable (run with the fenicsx_env interpreter, spec section 12)."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from reentry_model import material, thermal
from reentry_model.thermal import SIGMA_SB

pytest.importorskip("dolfinx")

R = 0.05
RHO, CP, K0, EPS = 2813.0, 900.0, 160.0, 0.4


def constant_material(k=K0, cp=CP, rho=RHO, eps=EPS):
    return material.Material("const", rho, eps, [100.0, 20000.0], [cp, cp], [100.0, 20000.0], [k, k])


def solver(name, mesh, mat, **options):
    s = thermal.thermal_solver(name, **options)
    s.setup(mesh, mat, mat.emissivity)
    return s


def test_missing_backend_message_is_not_raised_here():
    assert thermal.thermal_solver("fenicsx").linear_solver == "amg"
    with pytest.raises(ValueError):
        thermal.thermal_solver("fenicsx", lumped_mass=True)


def test_temperature_round_trip_uses_the_mesh_node_order(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, constant_material())
    T = 300.0 + 100.0 * coarse_sphere_mesh.points[:, 0] / R
    s.set_temperature(T)
    assert np.allclose(s.temperature(), T)
    Te = T[coarse_sphere_mesh.tets].mean(axis=1)                       # the cell set is the same, whatever dolfinx's cell order
    expected = RHO * CP * float((coarse_sphere_mesh.element_volumes() * (Te - material.T_REF)).sum())
    assert s.energy() == pytest.approx(expected, rel=1e-6)


def test_lumped_limit_matches_the_lumped_ode(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, constant_material(k=K0 * 1e4), linear_solver="direct", newton_tol=1e-8)
    s.set_temperature(300.0)
    A, m = s.areas.sum(), RHO * s.cell_volumes.sum()
    q = np.full(len(s.faces), 1e6)
    for _ in range(100):
        res = s.step(0.5, q, 0.0)
    exact = solve_ivp(lambda t, y: [(1e6 * A - EPS * SIGMA_SB * A * y[0] ** 4) / (m * CP)], (0.0, 50.0), [300.0], rtol=1e-12, atol=1e-10).y[0, -1]
    assert s.energy() / (m * CP) + material.T_REF == pytest.approx(exact, rel=1e-3) and res.Q_conv == pytest.approx(1e6 * A, rel=1e-12)


def test_radiative_cooling_of_an_isothermal_sphere(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, constant_material(k=K0 * 1e4), linear_solver="direct", newton_tol=1e-8)
    s.set_temperature(1500.0)
    A, m = s.areas.sum(), RHO * s.cell_volumes.sum()
    for _ in range(200):
        s.step(0.5, np.zeros(len(s.faces)), 0.0)
    exact = (1500.0 ** -3 + 3.0 * EPS * SIGMA_SB * A * 100.0 / (m * CP)) ** (-1.0 / 3.0)
    assert s.energy() / (m * CP) + material.T_REF == pytest.approx(exact, rel=2e-3)


def test_energy_balance_every_step(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, material.Material.from_drama_json(), newton_tol=1e-12)
    s.set_temperature(300.0)
    theta = np.arccos(np.clip(coarse_sphere_mesh.surface().normals[:, 0], -1.0, 1.0))
    q = 2e6 * np.where(theta < math.pi / 2, np.cos(theta), 0.0)
    for _ in range(20):
        E0 = s.energy()
        res = s.step(0.5, q, 0.0)
        assert abs((s.energy() - E0) - (res.Q_conv - res.Q_rad) * 0.5) < 1e-6 * abs((res.Q_conv - res.Q_rad) * 0.5)


def test_carslaw_jaeger_with_dirichlet(uniform_test_mesh):
    s = solver("fenicsx", uniform_test_mesh, constant_material(eps=0.0), newton_tol=1e-10)
    s.set_temperature(300.0)
    alpha = K0 / (RHO * CP)
    boundary, centre = uniform_test_mesh.boundary_nodes(), uniform_test_mesh.centre_node()
    for i in range(1, 81):
        s.step(0.05, np.zeros(len(s.faces)), 0.0, dirichlet=(boundary, np.full(boundary.size, 800.0)))
    exact = 800.0 - 500.0 * 2.0 * sum((-1) ** (n + 1) * math.exp(-n * n * math.pi ** 2 * alpha * 4.0 / R ** 2) for n in range(1, 80))
    assert abs(s.temperature()[centre] - exact) < 0.01 * 500.0 and np.all(s.temperature()[boundary] == 800.0)


def test_cross_check_against_scikit_fem(coarse_sphere_mesh):
    """Same mesh, same AA7075 material, same cos(theta) loads for 20 steps of 0.5 s: energy-equivalent temperature,
    surface extremes and every nodal temperature agree to 0.1 %."""
    mat = material.Material.from_drama_json()
    theta = np.arccos(np.clip(coarse_sphere_mesh.surface().normals[:, 0], -1.0, 1.0))
    q = 1.5e6 * np.where(theta < math.pi / 2, np.cos(theta), 0.0)
    fields = {}
    for name in ("skfem", "fenicsx"):
        s = solver(name, coarse_sphere_mesh, mat, newton_tol=1e-10)
        s.set_temperature(300.0)
        for _ in range(20):
            s.step(0.5, q, 0.0)
        fields[name] = (s.temperature(), s.energy())
    T1, E1 = fields["skfem"]
    T2, E2 = fields["fenicsx"]
    assert abs(E2 / E1 - 1.0) < 1e-3 and np.abs(T2 - T1).max() / T1.max() < 1e-3
    assert abs(T2.max() / T1.max() - 1.0) < 1e-3 and abs(T2.min() / T1.min() - 1.0) < 1e-3
```

- [ ] **Step 2: Run them**

Run: `"$PY" -m pytest tests/test_reentry_model_fenicsx.py -q`
Expected: `1 skipped` (module-level importorskip).

- [ ] **Step 3: Replace `reentry_model/thermal/fenicsx_backend.py`**

```python
"""FEniCSx (dolfinx >= 0.9) backend: the scheme of skfem_backend written in UFL on the same mesh (spec section 8).

Per Newton iterate the linearised, symmetric positive definite system
    a(u, v) = rho c u v / dt + k(T_k) grad u . grad v + 4 eps sigma T_k^3 u v |_Gamma
    L(v)    = rho c T_old v / dt + eps sigma (3 T_k^4 + T_amb^4) v |_Gamma  (+ the nodal convective loads)
is assembled with dolfinx and solved with PETSc (CG + hypre BoomerAMG, or LU); k and the secant heat capacity
c = [h(T_k) - h(T_old)] / (T_k - T_old) are DG0 cell coefficients refreshed every iteration, exactly as in the skfem
backend, so both backends discretise the volume terms identically; the radiation term is integrated by quadrature
on nodal T (skfem uses the facet mean) and the convective load is the same nodal vector A_f/3 per facet node.
dolfinx renumbers vertices: `node_of_dof`/`dof_of_node` map between the VolumeMesh's node ids and the P1 dofs.
Serial by default; every dolfinx/PETSc call is MPI-aware by construction (not exercised here). Lumped mass is not
implemented in this backend. `dolfinx` is imported lazily: the constructor raises MissingBackend without it."""
import numpy as np

from . import SIGMA_SB, MissingBackend, StepResult


class FenicsxThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=False, newton_tol=1e-6, max_iterations=30, cg_tol=1e-10, **_):
        try:
            import dolfinx  # noqa: F401
            import ufl  # noqa: F401
            from mpi4py import MPI  # noqa: F401
            from petsc4py import PETSc  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        if lumped_mass:
            raise ValueError("lumped mass is not implemented in the fenicsx backend")
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        self.linear_solver, self.newton_tol, self.max_iterations, self.cg_tol = linear_solver, newton_tol, max_iterations, cg_tol

    def setup(self, mesh, material, emissivity):
        import basix.ufl
        import ufl
        from dolfinx import fem
        from dolfinx import mesh as dmesh
        from mpi4py import MPI
        from petsc4py import PETSc
        from scipy.spatial import cKDTree
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
        self.msh = dmesh.create_mesh(MPI.COMM_WORLD, mesh.tets.astype(np.int64), mesh.points, domain)
        self.V = fem.functionspace(self.msh, ("Lagrange", 1))
        self.V0 = fem.functionspace(self.msh, ("DG", 0))
        n_local = self.V.dofmap.index_map.size_local
        _, self.node_of_dof = cKDTree(mesh.points).query(self.V.tabulate_dof_coordinates()[:n_local])
        self.dof_of_node = np.empty(n_local, dtype=np.int64)
        self.dof_of_node[self.node_of_dof] = np.arange(n_local)
        self.cell_dofs = np.asarray(self.V.dofmap.list)[:, :4]
        x = self.V.tabulate_dof_coordinates()[self.cell_dofs]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        self.cell_volumes = np.abs(np.linalg.det(J)) / 6.0
        surface = mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        self.T, self.T_old, self.T_k = fem.Function(self.V), fem.Function(self.V), fem.Function(self.V)
        self.k_fun, self.c_fun = fem.Function(self.V0), fem.Function(self.V0)
        self.dt_c, self.T_amb_c = fem.Constant(self.msh, PETSc.ScalarType(1.0)), fem.Constant(self.msh, PETSc.ScalarType(0.0))
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        es = self.emissivity * SIGMA_SB
        self.a = fem.form(self.c_fun / self.dt_c * u * v * ufl.dx + self.k_fun * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
                          + 4.0 * es * self.T_k ** 3 * u * v * ufl.ds)
        self.L = fem.form(self.c_fun / self.dt_c * self.T_old * v * ufl.dx + es * (3.0 * self.T_k ** 4 + self.T_amb_c ** 4) * v * ufl.ds)
        self.rad_form = fem.form(es * (self.T ** 4 - self.T_amb_c ** 4) * ufl.ds)
        self.ksp = PETSc.KSP().create(self.msh.comm)
        if self.linear_solver == "direct":
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.setType("cg")
            self.ksp.getPC().setType("hypre")
            self.ksp.getPC().setHYPREType("boomeramg")
            self.ksp.setTolerances(rtol=self.cg_tol, max_it=500)
        self.T.x.array[:] = 300.0
        self._T_prev = None

    def set_temperature(self, T):
        self.T.x.array[:] = float(T) if np.ndim(T) == 0 else np.asarray(T, dtype=float)[self.node_of_dof]
        self.T.x.scatter_forward()
        self._T_prev = None

    def temperature(self):
        return self.T.x.array[self.dof_of_node].copy()

    def facet_load(self, q):
        """Nodal loads A_f/3 per facet node, in the VolumeMesh's node order."""
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=len(self.mesh.points))

    def _coefficients(self, T_dofs, T_old_dofs):
        Te, Te_old = T_dofs[self.cell_dofs].mean(axis=1), T_old_dofs[self.cell_dofs].mean(axis=1)
        dT = Te - Te_old
        moved = np.abs(dT) > 1e-9
        c = np.where(moved, (self.material.enthalpy(Te) - self.material.enthalpy(Te_old)) / np.where(moved, dT, 1.0), self.material.cp_eff(Te))
        self.k_fun.x.array[:] = self.material.k(Te)
        self.c_fun.x.array[:] = self.material.rho * c

    def energy(self):
        Te = self.T.x.array[self.cell_dofs].mean(axis=1)
        return float((self.material.rho * self.cell_volumes * self.material.enthalpy(Te)).sum())

    def radiated_power(self, T_amb):
        from dolfinx import fem
        from mpi4py import MPI
        self.T_amb_c.value = T_amb
        return float(self.msh.comm.allreduce(fem.assemble_scalar(self.rad_form), op=MPI.SUM))

    def step(self, dt, q_conv, T_amb, dirichlet=None):
        from dolfinx import fem
        from dolfinx.fem.petsc import apply_lifting, assemble_matrix, assemble_vector, set_bc
        from petsc4py import PETSc
        self.dt_c.value, self.T_amb_c.value = dt, T_amb
        self.T_old.x.array[:] = self.T.x.array
        T_old = self.T.x.array.copy()
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))[self.node_of_dof]
        bcs = []
        if dirichlet is not None:
            g = fem.Function(self.V)
            dofs = self.dof_of_node[np.asarray(dirichlet[0])]
            g.x.array[dofs] = dirichlet[1]
            bcs = [fem.dirichletbc(g, dofs.astype(np.int32))]
        T_new, iteration = fem.Function(self.V), 0
        for iteration in range(1, self.max_iterations + 1):
            self.T_k.x.array[:] = T_k
            self._coefficients(T_k, T_old)
            A = assemble_matrix(self.a, bcs=bcs)
            A.assemble()
            b = assemble_vector(self.L)
            b.array[:] += F_conv
            apply_lifting(b, [self.a], bcs=[bcs])
            b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            set_bc(b, bcs)
            self.ksp.setOperators(A)
            self.ksp.solve(b, T_new.x.petsc_vec)
            T_new.x.scatter_forward()
            if np.linalg.norm(T_new.x.array - T_k) <= self.newton_tol * np.linalg.norm(T_new.x.array):
                break
            T_k = T_new.x.array.copy()
        self._T_prev = T_old
        self.T.x.array[:] = T_new.x.array
        return StepResult(self.temperature(), float(F_conv.sum()), self.radiated_power(T_amb), iteration)
```

- [ ] **Step 4: Check the module imports and the missing-backend path**

Run: `"$PY" -c "from reentry_model.thermal import fenicsx_backend; print('import ok')"` and `"$PY" -m pytest tests/test_reentry_model_fenicsx.py tests/test_reentry_model_cli.py -q -k fenicsx`
Expected: `import ok`; `1 passed, 1 skipped` (the CLI still exits 2 with the `fenicsx_env` message).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/thermal/fenicsx_backend.py tests/test_reentry_model_fenicsx.py
git commit -m "reentry_model: FEniCSx conduction backend (dolfinx, lazy) with conformance tests skipped without dolfinx

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 12: Verification against SESAM, README, spec and facts-note amendments

**Files:**
- Create: `tests/test_reentry_model_reference_thermal.py`
- Create: `analysis/reentry_model_thermal_verification.py`
- Modify: `pytest.ini` (marker description), `README.md` (new Step 2 section + Tests section), `docs/superpowers/specs/2026-09-18-thermal-fem-design.md` (amendments), `Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md` (append §15; outside the git repo, no commit)

**Interfaces:**
- Consumes: everything. Produces: the recorded verification (JSON under `reentry_model_output/verification_thermal/`, README table), the reference-tier tests with their thresholds.

- [ ] **Step 1: Write the reference-tier tests**

Create `tests/test_reentry_model_reference_thermal.py`:

```python
"""Coupled model vs the two no-melt US76 SESAM references (marker: reference, ~25 min): SESAM-equivalent heating on
the default mesh against the acceptance thresholds of spec section 10, physics mode reported (not thresholded),
and the mesh / time-step refinement checks on the 100 mm case. Each run's metrics are also written to
reentry_model_output/verification_thermal/ for the README table."""
import json
import os

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, compare, coupled, heating, material, mesh, sesam_io, thermal
from reentry_model import trajectory as tj

NAMES = {
    "d100": "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind",
    "d050": "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_nowind",
}
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "verification_thermal")
# Spec section 10 expectations, except the radiated power: the resolved surface radiates at its own (hotter)
# temperature, so a 2 % temperature margin is a 4 x 2 % = 8 % margin on eps sigma T^4 (measured 6.7 % / 3.4 % on
# 2026-09-18 while T_eq was within 1.2 %). Task 12 may revise a value only together with the measured number and
# the reason, recorded in the README verification table.
THRESHOLDS = {"Q_conv_peak_rel_max": 0.03, "Q_conv_continuum_rel_max": 0.03, "integrated_heat_rel_hypersonic": 0.03,
              "dT_rel_max": 0.02, "radiated_peak_rel_max": 0.08}


def load(key):
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, NAMES[key] + ".csv"))


def coupled_run(ref, heating_name, h_surface=mesh.DEFAULT_H_SURFACE, h_core=mesh.DEFAULT_H_CORE, dt=0.5, t_max=3600.0):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    settings = tj.Settings(diameter=ref.diameter, t_max=t_max)
    the_mesh = mesh.sphere_mesh(ref.diameter / 2.0, h_surface, h_core)
    the_body = body.ThermalBody(the_mesh, material.Material.from_drama_json(), thermal.thermal_solver("skfem"),
                                body.sphere_mass(ref.diameter, ref.material_density))
    model = heating.SesamEquivalentHeating() if heating_name == "sesam" else heating.PhysicsHeating()
    sim = tj.Simulator(initial, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(), settings)
    return coupled.CoupledRun(sim, the_body, model, coupled.CoupledSettings(dt=dt)).run()


def record(label, history, ref):
    os.makedirs(OUTDIR, exist_ok=True)
    doc = {"case": ref.name, "results": history.results, "thermal_metrics": compare.thermal_metrics(history, ref),
           "trajectory_metrics": compare.metrics(history, ref)}
    with open(os.path.join(OUTDIR, label + ".json"), "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
    return doc["thermal_metrics"]


@pytest.mark.reference
@pytest.mark.parametrize("key", list(NAMES))
def test_sesam_equivalent_mode_matches_sesam(key):
    ref = load(key)
    history = coupled_run(ref, "sesam")
    assert history.end_reason == "ground" and abs(history.results["energy_balance_residual"]) < 1e-6
    m = record(key + "__sesam", history, ref)
    assert m["Q_conv"]["max"] <= THRESHOLDS["Q_conv_peak_rel_max"], m["Q_conv"]
    assert m["Q_conv"]["continuum_rel_max"] <= THRESHOLDS["Q_conv_continuum_rel_max"], m["Q_conv"]
    assert abs(m["integrated_heat"]["rel_error_end_of_hypersonic"]) <= THRESHOLDS["integrated_heat_rel_hypersonic"], m["integrated_heat"]
    assert m["temperature"]["dT_rel_max"] <= THRESHOLDS["dT_rel_max"], m["temperature"]
    assert m["radiated"]["max"] <= THRESHOLDS["radiated_peak_rel_max"], m["radiated"]


@pytest.mark.reference
@pytest.mark.parametrize("key", list(NAMES))
def test_physics_mode_is_reported(key):
    """Lees x Fay-Riddell delivers less heat than SESAM's tumbling average: the ratio is recorded, not thresholded
    (expected 0.6-0.9: 0.196/0.2747 x Fay-Riddell/DKR x hot wall)."""
    ref = load(key)
    history = coupled_run(ref, "physics")
    m = record(key + "__physics", history, ref)
    assert history.end_reason == "ground" and 0.4 < m["integrated_heat"]["ratio_end"] < 1.1


@pytest.mark.reference
def test_mesh_and_time_step_refinement():
    """100 mm, SESAM-equivalent, to 200 s (past peak heating and peak surface temperature): halving h_surface
    (1 mm / 8 mm, 76 k nodes) changes the surface-temperature history by < 1 %; halving dt changes the peak surface
    temperature by < 0.5 %."""
    ref = load("d100")
    base = coupled_run(ref, "sesam", t_max=200.0)
    fine = coupled_run(ref, "sesam", h_surface=0.5 * mesh.DEFAULT_H_SURFACE, t_max=200.0)
    t = base.columns["time_s"]
    for col in ("surface_T_max_K", "surface_T_mean_K", "temperature_K"):
        d = np.abs(np.interp(t, fine.columns["time_s"], fine.columns[col]) - base.columns[col]) / base.columns[col]
        assert d.max() < 0.01, (col, d.max())
    half = coupled_run(ref, "sesam", dt=0.25, t_max=200.0)
    assert abs(half.results["peak_surface_T_K"] / base.results["peak_surface_T_K"] - 1.0) < 0.005
    with open(os.path.join(OUTDIR, "d100__refinement.json"), "w") as fh:
        json.dump({"base": base.results, "h_surface_halved": fine.results, "dt_halved": half.results}, fh, indent=2, default=str)
```

- [ ] **Step 2: Write the verification script**

Create `analysis/reentry_model_thermal_verification.py`:

```python
#!/usr/bin/env python3
"""Run the coupled thermal model against the two no-melt US76 SESAM references in both heating modes and tabulate
the heat/temperature errors (Step 2 verification, spec section 10).

    "$PY" analysis/reentry_model_thermal_verification.py [--outdir reentry_model_output/verification_thermal]
        [--cases d100,d050] [--modes sesam,physics] [--animate]

Cases left out by --cases/--modes are read back from an existing <outdir>/<case>__<mode>.json, so summary.md /
summary.json always cover every result available in the output directory.
"""
import argparse
import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

CASES = {
    "d100": ("sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind", ["--diameter", "100", "--altitude", "77.500133"]),
    "d050": ("sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_nowind", ["--diameter", "50", "--altitude", "115"]),
}
MODES = ("sesam", "physics")
COLUMNS = ["case", "mode", "Q_conv max (of peak)", "Q_conv max (continuum, point-wise)", "integrated heat (end of hypersonic)",
           "integrated heat (end)", "max |dT_eq|", "radiated max (of peak)", "peak surface T", "runtime"]


def run_case(key, mode, outdir, animate):
    name, args = CASES[key]
    ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    label = "{}__{}".format(key, mode)
    argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                             "--heating", mode, "--reference", ref, "--outdir", outdir, "--name", label, "--quiet"]
    if animate:
        argv.append("--animate")
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run {} failed with exit code {}".format(label, rc))
    return label


def row_from_doc(key, mode, doc):
    tm, res = doc["comparison"]["thermal_metrics"], doc["results"]
    return {"case": key, "mode": mode, "Q_conv max (of peak)": "{:.2%}".format(tm["Q_conv"]["max"]),
            "Q_conv max (continuum, point-wise)": "{:.2%}".format(tm["Q_conv"]["continuum_rel_max"]),
            "integrated heat (end of hypersonic)": "{:+.2%}".format(tm["integrated_heat"]["rel_error_end_of_hypersonic"]),
            "integrated heat (end)": "{:+.2%} (ratio {:.3f})".format(tm["integrated_heat"]["rel_error_end"], tm["integrated_heat"]["ratio_end"]),
            "max |dT_eq|": "{:.1f} K ({:.2%})".format(tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"]),
            "radiated max (of peak)": "{:.2%}".format(tm["radiated"]["max"]),
            "peak surface T": "{:.0f} K at {:.0f} s".format(res["peak_surface_T_K"], res["time_of_peak_surface_T_s"]),
            "runtime": "{:.0f} s, {} steps, {:.1f} it/step".format(res["runtime_s"], res["n_macro_steps"], res["mean_newton_iterations"])}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_thermal"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--modes", default=",".join(MODES))
    p.add_argument("--animate", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    wanted = {(k, m) for k in args.cases.split(",") for m in args.modes.split(",")}
    rows = []
    for key in CASES:
        for mode in MODES:
            label = "{}__{}".format(key, mode)
            if (key, mode) in wanted:
                run_case(key, mode, args.outdir, args.animate)
            path = os.path.join(args.outdir, label + ".json")
            if os.path.isfile(path):
                with open(path) as fh:
                    rows.append(row_from_doc(key, mode, json.load(fh)))
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    lines += ["| " + " | ".join(str(r[c]) for c in COLUMNS) + " |" for r in rows]
    with open(os.path.join(args.outdir, "summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 3: Update the marker description in `pytest.ini`**

Replace the `reference:` line with:

```
    reference: runs the full model against the committed SESAM reference flights (Step 1: eight trajectories, ~15 min; Step 2: the coupled thermal runs of tests/test_reentry_model_reference_thermal.py, ~35 min); no DRAMA needed; select with -m reference
```

- [ ] **Step 4: Run the verification script (both spheres, both modes, with the animation)**

Run:

```bash
"$PY" analysis/reentry_model_thermal_verification.py --animate
```

Expected (~12 min; measured 2026-09-18 with identical code, default mesh 2 mm/8 mm, Δt 0.5 s): the summary table with, for `sesam`: d100 — Q_conv max 0.46 % of peak, 2.2 % point-wise (continuum), integrated heat +0.31 % (hypersonic) / +0.30 % (end), |ΔT_eq| 24.5 K (1.2 %), radiated 6.7 %, peak surface T 2089 K at 139 s, ~170 s runtime, 732 steps, 2.1 it/step; d050 — 0.37 %, 2.4 %, +0.03 % / +0.02 %, 28 K (1.2 %), 3.4 %, 2442 K at 245 s, ~52 s, 1117 steps; for `physics`: integrated-heat ratio 0.744 (d100, peak stagnation 2319 K vs mean 1717 K) and ≈0.77 (d050). Also `reentry_model_output/verification_thermal/d100__physics/vtk/animation.mp4` and the stills. If any value differs from these by more than the noise of a different machine (runtimes) or 0.1 % (metrics), stop and investigate before writing the README.

- [ ] **Step 5: Run the reference-tier thermal tests**

Run: `"$PY" -m pytest tests/test_reentry_model_reference_thermal.py -m reference -q`
Expected: `5 passed` (~35 min: two SESAM-equivalent flights, two physics flights, and the refinement test — the 1 mm/8 mm mesh run to 200 s takes ~6 min). Measured refinement on 2026-09-18: halving h_surface changes the surface-temperature history by ≤ 3e-4 (centre 1.2e-3), halving Δt changes the peak surface temperature by 3.8e-4. If a threshold fails, revise it only together with the measured value and the reason, in the test's comment and the README table.

- [ ] **Step 6: Write the README section**

In `README.md`, insert the following section after the Step 1 verification notes (i.e. immediately before `## Tests`), replacing the numbers marked "measured" with the values from Steps 4–5 if they differ:

````markdown
## Physics model — `reentry_model` (Step 2: coupled 3D heat transfer)

Step 2 solves the trajectory and the temperature field inside the sphere together: every 0.5 s macro step the
trajectory advances (DOP853, Step 1 tolerances), an aerothermal model turns the freestream state into a convective
flux on each of the ~19 000 surface patches (100 mm sphere, 2 mm surface elements), and a finite-element conduction step (P1 tetrahedra, backward Euler,
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
temperature, the quantity SESAM's lumped model reports); `<run>/vtk/` holds `field.pvd` + `field_<k>.vtu` (nodal T) and
`surface.pvd` + `surface_<k>.vtp` (per-patch q_conv, q_rad, T), `<run>/vtk/animation.mp4` (GIF fallback), `frames/`,
`stills/`; with `--reference`, three more plots (`heating_time`, `temperature_time`, `integrated_heat`) and
`comparison.thermal_metrics` in the JSON. A 100 mm flight takes ~3 min (SESAM-equivalent) / ~5 min (physics with
animation); the FEniCSx backend needs the separate `fenicsx_env` (conda-forge `fenics-dolfinx`) and is untested until
that environment exists (selecting it in `drama_env` exits 2).

**`--heating sesam` is a verification device, not a physical model.** It applies SESAM's tumbling-average heat input —
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
| d100.00mm_h077.500km | sesam | 0.46 % | 2.2 % | +0.31 % / +0.30 % | 24.5 K (1.2 %) | 6.7 % | 2089 K at 139 s (SESAM lumped peak 2103 K) | 170 s, 732 steps |
| d050.00mm_h115.000km | sesam | 0.37 % | 2.4 % | +0.03 % / +0.02 % | 28.0 K (1.2 %) | 3.4 % | 2442 K at 245 s | 52 s, 1117 steps |
| d100.00mm_h077.500km | physics | — | — | ratio to SESAM 0.744 | — | — | 2319 K stagnation, mean 1717 K | 300 s incl. animation |
| d050.00mm_h115.000km | physics | — | — | ratio to SESAM 0.77 | — | — | 2493 K stagnation, mean 2151 K | 70 s |

Thresholds (`tests/test_reentry_model_reference_thermal.py`): Q_conv 3 % of peak and 3 % point-wise, integrated heat
3 %, T_eq 2 %, radiated power 8 % (= 4 × the temperature margin: the resolved surface radiates at its own, hotter
temperature — measured 6.7 % while T_eq was within 1.2 %). Refinement (100 mm, to 200 s): halving `h_surface`
(1 mm / 8 mm, 76 k nodes) changes the surface-temperature history by ≤ 0.03 %, halving Δt changes the peak surface
temperature by 0.04 %. Energy balance closes to 1e-8 over every flight.

Findings recorded while building this step (details in `sesam_verified_facts.md` §15 and the spec's amendments):
- SESAM's convective heating carries a hot-wall factor max(0, 1 − c_p(T − T∞)/(V²/2)) with c_p ≈ 1004.5 J/kg/K
  (invisible below 850 K, where the earlier facts were measured; it makes SESAM's heating vanish below ~1.8 km/s once
  the no-melt sphere is at 2000 K), and it halves the continuum heating below Mach 1.
- SESAM's transitional heating is not its drag blend: relative to 0.27471 × q_DKR it is 0.94 at Kn 0.04, 0.46 at
  Kn 0.2, 0.14 at Kn 1 and 0.059 at Kn 40 — the last being 0.78 × the textbook free-molecular cos θ average, so SESAM's
  free-molecular limit is the ordinary ½ρV³ with α ≈ 0.8 (the "13× too low q_FM" of the Step 1 facts was the
  transitional deficit misread at Kn 0.03). The measured factor F_h(Kn) (`aero.SesamHeatTable`) is what the
  verification mode uses; a first attempt with the textbook blend gave +12 % / +17 % integrated heat.
- Fay–Riddell with Cantera's equilibrium air is 1.29 × Sutton–Graves at the 100 mm start (72 % of the stagnation
  enthalpy is dissociation at 1.5 kPa; Sutton–Graves is a Le = 1 fit), 1.11 × at 4 km/s; `airNASA9.yaml` has no
  transport data, so viscosity comes from Blottner fits with Wilke mixing.
- Physics mode delivers 0.74–0.77 × SESAM's integrated heat (0.196/0.2747 × Fay–Riddell/DKR × hot wall) with a
  600 K stagnation-to-mean temperature difference at peak heating — the reason Step 3 needs the resolved field.
````

Then replace the `## Tests` section with:

````markdown
## Tests

```bash
"$PY" -m pytest -m "not drama and not reference" -q   # unit tests (~4 min; the thermal solver and coupled tests dominate)
"$PY" -m pytest -m "not drama"     # unit tests, no DRAMA needed (real SESAM outputs in tests/fixtures/)
"$PY" -m pytest                    # also the integration tests that run SESAM
"$PY" -m pytest -m reference -q   # the Step 1 reference flights (~15 min) and the Step 2 coupled runs (~35 min)
```

`tests/test_reentry_model_fenicsx.py` runs only with an interpreter that can import `dolfinx` (the `fenicsx_env`
environment); elsewhere it is skipped. To refresh a fixture see `tests/fixtures/README.md`.
````

Also update the Step 1 section's opening line `## Physics model — \`reentry_model\` (Step 1: trajectory)` to mention that Step 2 follows below (one clause), and its note "matting is reserved and not yet implemented — selecting it exits 1" stays true for the *drag* bridging (Matting is implemented for heat only).

- [ ] **Step 7: Amend the spec**

Edit `docs/superpowers/specs/2026-09-18-thermal-fem-design.md`:

1. Status line: `Status: implemented (plan docs/superpowers/plans/2026-09-18-thermal-fem.md); amended with the measurements listed in section 14`.
2. §2 table, row "SESAM heating": replace the cell with: `Total convective power Q = A_sphere · 0.27471 · q_DKR · F_h(Kn) · max(0, 1 − c_p(T − T∞)/(V²/2)) for Ma ≥ 1 (½ · 0.27471 · q_DKR below Ma 1), q_DKR = 1.1035e8 R^-1/2 (ρ/1.225)^1/2 (V/7925)^3.15 W/m², c_p = 1004.5 J/kg/K, F_h measured in 0.125-decade Kn bins (aero.SesamHeatTable: 1.005 continuum, 0.14 at Kn 1, 0.059 at Kn 40 = 0.78 × the free-molecular cos θ average). Supersedes the (1 − f)·q_DKR + f·q_FM decomposition of facts §4 (measured 2026-09-18, facts §15).`
3. §2 row "SESAM shape factor": replace "No hot-wall correction detectable." with "Hot-wall factor 1 − c_p(T − T∞)/(V²/2), clamped at 0 (facts §15)."
4. §5: replace "Defaults for the 100 mm sphere: `h_surface` = 1.0 mm, `h_core` = 8 mm (≈ 40 k nodes)" with "Defaults: `h_surface` = 2.0 mm, `h_core` = 8 mm (18.9 k nodes on the 100 mm sphere; 1 mm/8 mm is 76 k nodes and is the convergence mesh); a node is embedded at the centre".
5. §6.2: replace the definition paragraph's formulas with the amended form of item 2 (uniform 0.27471 · q_DKR · F_h · hot-wall on every patch; the hot-wall factor uses the body's energy-equivalent temperature as SESAM's lumped model does); keep the "verification device, not a physical model" paragraph.
6. §6.3: after "(Cantera, `airNASA9.yaml`)" add "for the thermodynamics and composition; that mechanism has no transport data, so μ_s and μ_w come from Blottner's curve fits (N2, O2, NO, N, O) with Wilke's mixing rule; the equilibrium solver is used only where T∞ + V²/(2c_p) > 1500 K"; add to the Fay–Riddell sentence "γ_cat scales the Lewis-number term of the equilibrium form (γ_cat = 0 gives the Le = 1 value; the frozen non-catalytic reduction 1 − h_D/h_s is a Step 3 option)"; add "each patch's flux is scaled by its own (h_s − h_w)/(h_s − h_w,stag)".
7. §6.4 outputs: add `absorbed_heat_J` = ∫(Q_conv − Q_rad)dt and `T_centre_K`; state that `integrated_heat_J` is ∫Q_conv dt (SESAM's meaning).
8. §7: replace "k(T), c_p(T) lagged from the previous Newton iterate (Picard)" with "k(T) at the element-mean temperature of the previous iterate and the secant heat capacity [h(T_k) − h(T_old)]/(T_k − T_old), so the discrete energy balance is exact to the Newton tolerance for the tabulated c_p"; replace "SciPy sparse direct (`--linear-solver direct`, default up to ~50 k nodes) or CG with a pyamg smoothed-aggregation preconditioner (`amg`)" with "CG with a pyamg smoothed-aggregation preconditioner rebuilt every 30 solves (`amg`, default; 0.03 s per solve at 19 k nodes) or SciPy SuperLU (`direct`, 0.4–24 s per solve at 12–76 k nodes, tests only); Newton starts from the extrapolated previous step (2.0–2.1 iterations per step)"; replace the cost target sentence with "Measured: 100 mm reference, default mesh, 732 steps, 170 s (SESAM-equivalent), 300 s (physics with animation); 50 mm from 115 km, 1117 steps, 52 s"; extend the history column list with `absorbed_heat_J` and `T_centre_K`.
9. §8 `skfem_backend`: replace "assembly per Newton iterate" with "unit element stiffness/mass matrices precomputed once and rescaled into a fixed CSR pattern per iterate (20 ms vs 230 ms for scikit-fem's generic `asm` at 12.6 k nodes, which is kept as the reference operator in the conformance test); radiation from the facet-mean temperature with its exact Jacobian". `fenicsx_backend`: add "written against dolfinx 0.9/0.10; untested until `fenicsx_env` exists; lumped mass not implemented there".
10. §9: run directory: "`<outdir>/<name>.csv/.json` as in Step 1; plots, `vtk/` (field.pvd + field_<k>.vtu, surface.pvd + surface_<k>.vtp), `vtk/animation.mp4`, `vtk/frames/`, `vtk/stills/` under `<outdir>/<name>/`"; metrics: "power errors relative to SESAM's peak over the hypersonic phase; point-wise Q_conv error over Kn_ref < 0.01 and Q_ref > 10 % of peak; integrated heat at the end of the hypersonic phase and at the end; max |T_eq − T_SESAM| and relative".
11. §10 heating checks: replace "Fay–Riddell within 15 % of Sutton–Graves (its fit)" with "Fay–Riddell 2.22e6 W/m² = 1.2–1.4 × Sutton–Graves (measured 1.29; 1.14 with the Lewis term off)"; item 3: "on a uniform 4 mm test mesh (centre 0.5 %, volume mean 0.2 %)"; items 1–2: "k × 1e4"; coupling thresholds: radiated power 8 % (reason as in the README); add the measured values of README's table.
12. §12: pin the versions of `requirements-step2.txt`.
13. Append `## 14. Amendments (2026-09-18)` listing items 1–12 of this plan's "Measured facts and spec amendments" in one line each.

- [ ] **Step 8: Append §15 to the facts note**

Append to `Literature Review/Sphere Demise Model - Planning References/sesam_facts/sesam_verified_facts.md` (outside the repo; no commit):

````markdown
## 15. SESAM's convective heating, re-measured on the no-melt US76 references (2026-09-18, Step 2 verification)
Measured on the two `_mAA7075_nomelt_nowind` references (US76, winds off; 100 mm from 77.5 km, 50 mm from 115 km),
where the temperature reaches 2100–2450 K instead of stopping at the 850 K melting point:

- **Hot-wall factor.** Q / (A · 0.27471 · q_DKR) = max(0, 1 − c_p (T_lumped − T∞) / (V²/2)) with c_p = 1004.5 J/kg/K
  (best fit 985; rms 0.008, max 0.015 over 189 rows with Ma ≥ 5, Kn < 0.01, up to 2450 K). It is > 0.97 below 850 K,
  which is why §4 and §6 found "no hot-wall correction"; it makes SESAM's heating exactly zero once
  c_p (T − T∞) > V²/2 (100 mm: from 1.8 km/s / 178 s at ~1950 K to the transonic phase).
- **Below Ma 1** the ratio is 0.500 ± 0.006 with no hot-wall term (both spheres; the same halving as C_D). In a narrow
  band 1.0 < Ma < 1.1 the factor is ≈ 1.0 (a few watts; not reproduced).
- **Transitional heating is not the drag blend.** With the hot-wall factor divided out, F_h(Kn) = Q / (A · 0.27471 ·
  q_DKR · hot-wall) in 0.125-decade bins of log10 Kn (371 rows, Ma ≥ 5, hot-wall factor > 0.3):
  1.0077 (Kn 0.002), 1.005 (0.012), 0.973 (0.027), 0.891 (0.049), 0.741 (0.087), 0.564 (0.15), 0.367 (0.27), 0.218 (0.49),
  0.144 (0.87), 0.133 (1.5), 0.117 (3.6), 0.095 (8.7), 0.072 (21), 0.059 (37). Interpolation reproduces the rows to
  rms 0.7 %, max 3.4 % (`reentry_model/aero.py: SesamHeatTable`). (1 − f_drag)·q_DKR + f_drag·q_FM with any constant
  q_FM does not fit (rms 23–28 %, up to 60 % at Kn ~0.8): SESAM's transitional heating lies far below both limits.
- **Free-molecular limit.** At Kn 35 (50 mm, 115 km) Q / A = 0.78 × [0.25 × ½ρV³], i.e. the textbook cos θ-averaged
  free-molecular flux with accommodation ≈ 0.8, decaying like Kn^-0.24 toward the transition. §4's "q_FM,avg
  anomalously small (0.075 × textbook)" was measured at Kn 0.03 on the 100 mm run and attributed the transitional
  deficit to the free-molecular coefficient; it is superseded by this section. The Klett coefficient of §10 is
  presumably the Kn-dependent transitional reduction, not a low free-molecular value.
- Consequences for the physics model: the Step 2 SESAM-equivalent verification mode uses the measured F_h(Kn) and the
  hot-wall factor (uniform over the surface, not physical); the physics mode (Fay–Riddell / Matting / Lees) is
  compared against SESAM only as a ratio (0.74–0.77 of SESAM's integrated heat).
````

- [ ] **Step 9: Run the whole unit tier once more**

Run: `"$PY" -m pytest -m "not drama and not reference" -q`
Expected: all pass.

- [ ] **Step 10: Commit**

```bash
git add tests/test_reentry_model_reference_thermal.py analysis/reentry_model_thermal_verification.py pytest.ini README.md docs/superpowers/specs/2026-09-18-thermal-fem-design.md
git commit -m "reentry_model: Step 2 verification against SESAM (thermal reference tests, script, README), spec amendments

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## After the plan

- Deferred to Step 3 (first task): prism-layer surface mesh; melting range in `Material.liquid_fraction` (the secant heat capacity already handles a latent-heat jump); surface recession through the mutable `VolumeMesh.points`; the frozen non-catalytic Fay–Riddell option; a `TabulatedHeating` table format.
- Deferred follow-ups from Step 1 still stand (replay density from 2q/V², smoothed replay T(h), pymsis per-step caching).
- `fenicsx_env` creation and the first run of `tests/test_reentry_model_fenicsx.py` need the user's go-ahead (spec §12).
