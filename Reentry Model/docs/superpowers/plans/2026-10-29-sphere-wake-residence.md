# Sphere Wake Residence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A general parent-plume / follower overlap engine, and the sweep it was built for: how long a solid AA7075
sphere released from the Generic Spacecraft stays inside the spacecraft's reentry plume, for 8 diameters and 6
release altitudes, with the wake gas's forward velocity and low density as a switchable effect on the sphere's drag.

**Architecture:** Pure modules (`laws.py`, `track.py`, `overlap.py`, `wake_gas.py`) hold the plume laws, the regime
selection, the track geometry, the union-of-balls margin, the residence search and the wake gas state; they import
nothing that knows about spheres, droplets, DRAMA or atmospheres. `adapters.py` connects them to the existing
integrators: a twin of the DRAMA parent integrated by a `trajectory.Simulator` subclass, a `WakeAtmosphere` wrapper
that feeds the wake gas to any follower through the atmosphere's `state()` call, and dense-output helpers.
`shock_layer.py` models the parent's shock-layer gas, which the same wrapper feeds to rarefied followers released on
the parent's windward surface (§4.3d).
`scenarios.py` runs one case; `analysis/` runs the sweep, sensitivity, tables and plots. `trajectory.py` and
`droplet_reentry/model.py` are not modified.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy (`solve_ivp`, `CubicHermiteSpline`, `cKDTree`,
`brentq`, `minimize_scalar`), matplotlib, pytest.

**Spec:** `docs/superpowers/specs/2026-10-29-sphere-wake-residence-design.md` (read it; this plan argues from it and
cites its sections as §n).

## Global Constraints

- Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; every command runs from inside `Reentry Model/`.
- Work on branch `wake-residence` (create it from `main`; a worktree via superpowers:using-git-worktrees is fine).
  Never stage `CLAUDE.md`, `reentry_model/droplet_reentry/model.py`,
  `reentry_model/droplet_reentry/tests/test_droplet_reentry_model.py`, the oxidation scripts or `Literature Review/`
  — they hold Asha's uncommitted oxidation work. Stage files by explicit path only, never `git add -A` or `git add .`.
- Do not modify `reentry_model/trajectory.py` or `reentry_model/droplet_reentry/model.py` (§4.3a).
- AFFIx release classifier (§4.3b, amended 2026-10-29): spacecraft releases windward only, theta 90-180 degrees;
  far-wake releases start on the wake axis two parent diameters behind the parent's centroid; escape releases at
  the windward surface point; shock-band cases not integrated; s_w 2 degrees, `mach_min` 5, fitted r range 1/16-1.
- Front plane (§4.3a): the plume exists only behind the parent's base plane; first exits through it are labelled
  `reached_parent` (fell behind, caught up) or `ahead_of_parent` (never fell behind).
- Follower regime (§4.3c, amended 2026-10-29): Kn_f = mean free path of the LOCAL gas at the release point (the
  shock-layer surface state) / follower diameter; `continuum` below 0.01, `free_molecular` at or above
  `kn_free_molecular` 10, `transitional` between. Continuum: AFFIx `governing`. Transitional: AFFIx `indicative`,
  flagged. Free-molecular: AFFIx `not_applicable`. Rarefied (transitional and free-molecular) followers start at the
  windward surface point at theta and are integrated through the shock-layer field (§4.3d) with the composed drag
  law, the pressure-gradient force and the parent as an obstacle (`hit_parent`).
- Shock-layer field (§4.3d) defaults: `shock_standoff` "lighthill", `phi_nose_deg` 40, `phi_blend_deg` 20,
  `shock_thickness_mfp` 3, `pressure_gradient` on, `follower_wall_T_K` 300, `kn_shock_layer_flag` 1e-3, dynamic-pressure
  change flag 0.10; the field is frozen at the release free stream and faded over r0 behind the base plane.
- `overlap.py`, `laws.py`, `track.py`, `wake_gas.py` and `affix.py` must not import `trajectory`, `droplet_reentry`,
  `parent_trajectory`, `body`, `aero` or `atmosphere`, and contain no altitude constants or altitude arguments (§4.0, §4.2).
- All code in SI units; cm Hg-cm appears only inside `pressure_parameter`.
- AA7075 solid density 2813 kg/m3 (`DRAMA Sphere Experiments/data/user_materials/AA7075_nomelt.json`).
- Parent frontal area 6.215 m2 (r0 = 1.407 m, equal-area diameter 2.813 m) (§4.2).
- Thresholds (named, overridable, §4.2): `kn_continuum` 0.01, `re_crit` 5e4, `p_lam` 0.5, `p_turb` 2.0 (cm Hg-cm);
  Taylor K 0.66, exponent 1/3, fitted to 2 <= x/sqrt(C_D A) <= 1200 and 3.7-6.2 km/s; Park C 0.7, exponent 1/3,
  `s_park_limit` 30e3 m, 5.8-7.9 km/s; Feldman f 0.8, Pr 0.71, `T_axis_hot` 3000 K, `laminar_edge_factor` 1;
  `u_min_parent` 1 m/s; `rarefied_fallback` None; `k_E` 1; `tau_pre` 60 s; residence tolerance 1e-3 s.
- Run names encode diameter, release altitude, bound, wake_flow, law and every non-default setting (§3).
- Summary JSON records `git_commit`, `git_dirty` (reuse `wake_column.run._git_commit` / `_git_dirty`) and the
  SHA-256 of the parent file (`ParentTrajectory.sha256`).
- Outputs under `output/wake_residence_output/` (`wake_residence.OUTPUT`), git-ignored by the existing `output/` rule.
- Commit messages: a scope prefix and a sentence describing the behaviour (`wake_residence: ...`), ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A track with non-increasing or duplicate times, or a non-finite position or flow value, must raise `ValueError`
   naming the first bad index rather than produce a silent wrong spline (Task 2).
2. A follower window that starts before the parent track begins, or ends after it ends, must raise `ValueError`
   naming both ranges; a follower that begins later than the parent is normal (Task 3).
3. A parent slower than `u_min_parent` at emission (a sphere parent settling at terminal speed) must give regime
   `"too_slow"`, no radius, and no division by zero in Feldman's x/(f U) (Task 1).
4. A follower that moves ahead of the parent from release must report a small positive residence (time for the
   margin to go negative), not "never inside" (Task 3).
5. A release altitude outside the twin's integrated range must raise `ValueError` naming the range, not
   extrapolate the track (Task 5).

---

## File Structure

```
reentry_model/wake_residence/
  __init__.py        package docstring (module map), OUTPUT path
  laws.py            EmissionState, LawSettings, regime criteria, Park / Taylor / Feldman radii, plume_radius()
  track.py           Track: interpolated path, path length, per-sample flow state, interpolation error
  overlap.py         Plume (radius, margin, wake_state) and residence() with its result records
  wake_gas.py        Ambient, GasState, gas_state(): momentum/energy-conserving Gaussian wake profiles
  adapters.py        DRAMA beta table, TwinSimulator, dense_solve, twin_parent, release_state, WakeAtmosphere,
                     sphere_simulator, solution_track, sphere_parent_track
  affix.py           Register et al. (2020) release classifier: zone, interaction time, validity, release offsets
  shock_layer.py     the parent's shock-layer field (§4.3d): Billig shape, equilibrium oblique jump, Lighthill nose
                     solution, Newtonian surface state, tabulated field with smoothing and pressure gradient
  scenarios.py       ScenarioSettings, run_name, sphere_in_spacecraft_plume, droplet_in_sphere_plume, write_case
  analysis/__init__.py
  analysis/sweep.py        default sweep (resumable, --workers), summary.csv, crossover.csv
  analysis/sensitivity.py  §4.6 set, sensitivity.csv
  analysis/tables.py       plume_size.csv, ballistic_coefficient.csv, twin_vs_drama.csv
  analysis/plots.py        the §2 goal-6 plots into <out>/plots/
  tests/__init__.py, tests/conftest.py, tests/test_laws.py, tests/test_track.py, tests/test_overlap.py,
  tests/test_wake_gas.py, tests/test_twin.py, tests/test_follower.py, tests/test_affix.py, tests/test_shock_layer.py,
  tests/test_rarefied_follower.py, tests/test_scenarios.py,
  tests/test_analysis.py
Modify: pytest.ini (testpaths), docs/model_assumptions.md (new section 10), README.md (root; commands),
        CLAUDE.md (layout + commands; committed on the branch — ruling R3)
```

---

### Task 1: Package skeleton, plume laws and regime selection

**Files:**
- Create: `reentry_model/wake_residence/__init__.py`, `reentry_model/wake_residence/laws.py`,
  `reentry_model/wake_residence/tests/__init__.py`, `reentry_model/wake_residence/tests/conftest.py`,
  `reentry_model/wake_residence/tests/test_laws.py`
- Modify: `pytest.ini` (append `reentry_model/wake_residence/tests` to `testpaths`)

**Interfaces:**
- Produces (in `laws.py`):
  - `@dataclass(frozen=True) EmissionState(rho, T, p, mu, mfp, U, d, cda, r0)` — parent's ambient density
    [kg/m3], temperature [K], pressure [Pa], viscosity [Pa s], mean free path [m], air-relative speed [m/s],
    diameter [m], drag area C_D A [m2], equal-area radius [m] at one emission time.
  - `@dataclass(frozen=True) LawSettings` with the Global-Constraints defaults as fields: `kn_continuum, re_crit,
    p_lam, p_turb, turbulent_law ("taylor"|"park"), taylor_K, taylor_n, taylor_xmax_over_L=1200.0,
    taylor_speed_range=(3.7e3, 6.2e3), park_C, park_n, s_park_limit, park_speed_range=(5.8e3, 7.9e3), feldman_f,
    feldman_Pr, T_axis_hot, laminar_edge_factor, rarefied_fallback (None|"taylor"|"park"|"feldman_low"|"feldman_high"),
    u_min_parent`; method `non_default() -> dict` (fields differing from defaults, for run names).
  - `sutherland_mu(T) -> float` (Feldman Eq. 13: 1.72e-5 Pa s × (T/273)^1.5 × 386/(113+T)).
  - `knudsen(s: EmissionState) -> float` (mfp/d), `reynolds(s) -> float` (rho U d / mu),
    `pressure_parameter(s) -> float` (p × sqrt(cda), converted Pa·m to cm Hg·cm by ×100/1333.22).
  - `classify(s, settings) -> str` in `{"rarefied","too_slow","laminar","transitional","turbulent"}`
    (order: U < u_min_parent → too_slow; Kn > kn_continuum → rarefied; then the Re/P rule of §4.2).
  - `park_radius(s_behind, r0, C, n) -> float`, `taylor_radius(x, cda, r0, K, n) -> float`,
    `feldman_radius(x, s, T_axis, f, Pr, edge_factor) -> float`.
  - `@dataclass(frozen=True) RadiusResult(r: float, regime: str, laws: tuple[str, ...], flags: dict)` — `r` is NaN
    when no law applies; `flags` keys: `speed_out_of_range`, `distance_out_of_range`, `feldman_extrapolation`
    (U outside 4.6-10.7 km/s, i.e. Feldman's 15,000-35,000 ft/s, or d outside his hemisphere-cylinder sizes —
    take both ranges from `citations/laminar_wake_growth_laws.md` §1 and name them as module constants).
  - `plume_radius(x, s, settings, bound) -> RadiusResult` — `bound` is `"low"` or `"high"`; laminar → min/max of
    Feldman at T_axis = s.T and T_axis = settings.T_axis_hot; transitional → min/max of those two and the
    turbulent law; turbulent → the turbulent law; rarefied → fallback law if set, else NaN; too_slow → NaN.
- Produces (in `__init__.py`): `OUTPUT = os.path.join(OUTPUT_DIR, "wake_residence_output")`.
- `tests/conftest.py`: copy the pattern of `reentry_model/wake_column/tests/conftest.py` (project root on
  `sys.path`, a session `parent` fixture).

- [ ] **Step 1: Write the failing tests** in `tests/test_laws.py`:

```python
import math, pytest
from reentry_model.wake_residence.laws import (EmissionState, LawSettings, classify, park_radius, taylor_radius,
    feldman_radius, plume_radius, sutherland_mu, pressure_parameter)

def state(**kw):
    base = dict(rho=1e-4, T=220.0, p=6.3, mu=sutherland_mu(220.0), mfp=1e-3, U=6000.0, d=2.813, cda=5.5, r0=1.4065)
    base.update(kw); return EmissionState(**base)

def test_park_radius_equals_r0_at_parent_and_scales_as_cube_root():
    assert park_radius(0.0, 1.4, 0.7, 1/3) == pytest.approx(1.4, rel=1e-12)
    s0 = (1/0.7) ** 3 * 1.4
    big = park_radius(8000.0 - s0, 1.4, 0.7, 1/3) / park_radius(1000.0 - s0, 1.4, 0.7, 1/3)
    assert big == pytest.approx(2.0, rel=1e-12)
    assert park_radius(500.0, 1.4, 0.9, 1/3) > park_radius(500.0, 1.4, 0.7, 1/3)

def test_taylor_radius_value_and_floor():
    L = math.sqrt(5.5)
    assert taylor_radius(1000.0, 5.5, 1.4, 0.66, 1/3) == pytest.approx(0.5 * 0.66 * L * (1000.0 / L) ** (1/3))
    assert taylor_radius(0.1, 5.5, 1.4, 0.66, 1/3) == pytest.approx(1.4)        # never below r0

def test_feldman_square_grows_linearly_and_bounds_order():
    s = state()
    r1, r2 = (feldman_radius(x, s, 3000.0, 0.8, 0.71, 1.0) for x in (1000.0, 2000.0))
    kappa = sutherland_mu(3000.0) / (s.rho * 0.71)
    assert r2**2 - r1**2 == pytest.approx(4 * kappa * 1000.0 / (0.8 * s.U), rel=1e-12)
    assert feldman_radius(500.0, s, s.T, 0.8, 0.71, 1.0) <= feldman_radius(500.0, s, 3000.0, 0.8, 0.71, 1.0)
    lo = plume_radius(500.0, s, LawSettings(T_axis_hot=s.T), "low")
    hi = plume_radius(500.0, s, LawSettings(T_axis_hot=s.T), "high")
    assert lo.r == pytest.approx(hi.r)

@pytest.mark.parametrize("kw,expected", [
    (dict(mfp=0.1), "rarefied"),                                  # Kn = 0.036 > 0.01
    (dict(rho=1e-5, p=0.5), "laminar"),                           # low Re and low P
    (dict(rho=1e-2, p=300.0), "turbulent"),                       # high Re and high P
    (dict(rho=1e-2, p=0.5), "transitional"),                      # Re says turbulent, P says laminar
    (dict(rho=1e-6, p=300.0), "transitional"),                    # Re says laminar, P says turbulent
    (dict(U=0.5), "too_slow"),
])
def test_classify_uses_only_flow_quantities(kw, expected):
    assert classify(state(**kw), LawSettings()) == expected

def test_small_body_is_rarefied_where_large_body_is_not():
    assert classify(state(d=0.005, r0=0.0025, cda=2e-5), LawSettings()) == "rarefied"
    assert classify(state(), LawSettings()) != "rarefied"

def test_thresholds_are_overridable():
    s = state(rho=1e-2, p=0.5)                                    # transitional by default
    assert classify(s, LawSettings(p_turb=1e-3, p_lam=1e-4)) == "turbulent"
    assert classify(state(mfp=0.1), LawSettings(kn_continuum=1.0)) != "rarefied"

def test_transitional_range_holds_all_three_laws():
    s = state(rho=1e-2, p=0.5)
    lo, hi = (plume_radius(800.0, s, LawSettings(), b) for b in ("low", "high"))
    assert lo.regime == hi.regime == "transitional"
    assert set(hi.laws) == {"feldman_low", "feldman_high", "taylor"} and lo.r <= hi.r

def test_rarefied_has_no_radius_unless_fallback():
    s = state(mfp=0.1)
    assert math.isnan(plume_radius(100.0, s, LawSettings(), "low").r)
    assert plume_radius(100.0, s, LawSettings(rarefied_fallback="taylor"), "low").r > 0

def test_too_slow_parent_gives_no_radius_and_no_error():
    assert math.isnan(plume_radius(100.0, state(U=0.0), LawSettings(), "high").r)

def test_pressure_parameter_units():
    # 1333.22 Pa x 1 m = 1 cm Hg x 100 cm
    assert pressure_parameter(state(p=1333.22, cda=1.0)) == pytest.approx(100.0)

def test_speed_and_distance_flags():
    r = plume_radius(1e5, state(rho=1e-2, p=300.0, U=7000.0), LawSettings(), "low")
    assert r.flags["speed_out_of_range"] and r.flags["distance_out_of_range"]
```

- [ ] **Step 2: Run them to confirm they fail**

Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_laws.py -q`
Expected: collection error (`ModuleNotFoundError: reentry_model.wake_residence`).

- [ ] **Step 3: Implement** `__init__.py`, `laws.py`, `tests/__init__.py`, `tests/conftest.py`, and add the test path
to `pytest.ini`. Each law's coefficient comes from `LawSettings`; the module docstring cites the three sources
(Taylor et al. 1964 Eq. 4-5 p. 1737; Park 1971 Eq. 11; Feldman 1961 Eqs. 6, 8, 11, 14 as assembled in
`citations/laminar_wake_growth_laws.md`) and states that R is density-weighted (§4.2).

- [ ] **Step 4: Run them to confirm they pass**

Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_laws.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add pytest.ini reentry_model/wake_residence/__init__.py reentry_model/wake_residence/laws.py \
  reentry_model/wake_residence/tests/__init__.py reentry_model/wake_residence/tests/conftest.py \
  reentry_model/wake_residence/tests/test_laws.py
git commit -m "wake_residence: plume laws (Taylor, Park, Feldman bounds) and regime selection from Kn, Re_d and Taylor's pressure parameter, never altitude"
```

---

### Task 2: Track

**Files:**
- Create: `reentry_model/wake_residence/track.py`, `reentry_model/wake_residence/tests/test_track.py`

**Interfaces:**
- Consumes: `laws.EmissionState`.
- Produces: `class Track` with
  - `Track.from_arrays(t, r, v=None, speed=None, flow=None) -> Track` — `t` (N,), `r` (N,3), `v` (N,3) optional,
    `speed` (N,) air-relative optional (default |v|, or |dr/dt| of the spline), `flow` optional dict of (N,) arrays
    with keys `rho, T, p, mu, mfp, d, cda, r0`. Raises `ValueError` naming the first index on non-increasing `t`
    or any non-finite value (Review Focus 1).
  - `Track.from_callable(position, t0, t1, velocity=None, speed=None, flow=None, n_samples=2001) -> Track` —
    `position(t)` returns (3,); the callable is used for positions; samples feed the k-d tree and path length.
  - attributes `t0, t1`; methods `position(t)`, `velocity(t)`, `speed(t)` (scalar or vectorised),
    `path_length(t)` (∫ speed dt from t0, monotone), `emission_state(tau) -> EmissionState` (interpolated flow,
    `U = speed(tau)`; raises if the track has no flow), `samples(spacing_m) -> (tau: ndarray, xyz: ndarray)` at
    roughly that path spacing, `interpolation_error() -> float` metres.
- Interpolation: `scipy.interpolate.CubicHermiteSpline` when `v` is given, `CubicSpline` otherwise; path length by
  cumulative integration of `speed` (Hermite/trapezoid on the samples is fine if the test below passes).
- `interpolation_error()`: rebuild the interpolant from every other sample and return the largest distance between
  the two interpolants at the dropped samples.

- [ ] **Step 1: Write the failing tests** in `tests/test_track.py`:

```python
import numpy as np, pytest
from reentry_model.wake_residence.track import Track

def circle(n):
    t = np.linspace(0, 10, n); R, w = 6.4e6, 7500 / 6.4e6
    r = np.c_[R*np.cos(w*t), R*np.sin(w*t), 0*t]; v = np.c_[-R*w*np.sin(w*t), R*w*np.cos(w*t), 0*t]
    return t, r, v

def test_straight_line_is_exact_and_path_length_is_speed_times_time():
    t = np.linspace(0, 2, 5); r = np.c_[7000*t, 0*t, 0*t]
    tr = Track.from_arrays(t, r, v=np.tile([7000, 0, 0], (5, 1)))
    assert np.allclose(tr.position(1.234), [8638.0, 0, 0], atol=1e-9)
    assert tr.path_length(1.5) == pytest.approx(10500.0, rel=1e-9)

def test_hermite_on_orbit_arc_is_millimetre_accurate():
    t, r, v = circle(101)
    tr = Track.from_arrays(t, r, v=v)
    tt = 3.21; R, w = 6.4e6, 7500 / 6.4e6
    assert np.linalg.norm(tr.position(tt) - [R*np.cos(w*tt), R*np.sin(w*tt), 0]) < 1e-3

def test_interpolation_error_flags_coarse_sampling():
    t, r, v = circle(6)
    assert Track.from_arrays(t, r).interpolation_error() > Track.from_arrays(*circle(201)[:2]).interpolation_error()

@pytest.mark.parametrize("bad", ["dup", "nan"])
def test_bad_input_names_index(bad):
    t = np.array([0.0, 1.0, 2.0, 3.0]); r = np.zeros((4, 3)); r[:, 0] = t
    if bad == "dup": t[2] = 1.0
    else: r[3, 1] = np.nan
    with pytest.raises(ValueError, match="index"):
        Track.from_arrays(t, r)

def test_emission_state_interpolates_flow():
    t = np.array([0.0, 1.0]); r = np.c_[[0.0, 1.0], [0, 0], [0, 0]]
    flow = {k: np.array([1.0, 3.0]) for k in ("rho", "T", "p", "mu", "mfp", "d", "cda", "r0")}
    s = Track.from_arrays(t, r, flow=flow).emission_state(0.5)
    assert s.rho == pytest.approx(2.0) and s.U == pytest.approx(1.0)
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_track.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement `track.py`** to the interface above.
- [ ] **Step 4: Run them to confirm they pass.** Same command. Expected: all pass.
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/track.py reentry_model/wake_residence/tests/test_track.py
git commit -m "wake_residence: Track, an interpolated path with path length, per-sample flow state and an interpolation-error estimate"
```

---

### Task 3: Overlap engine (margin and residence)

**Files:**
- Create: `reentry_model/wake_residence/overlap.py`, `reentry_model/wake_residence/tests/test_overlap.py`

**Interfaces:**
- Consumes: `Track`, `LawSettings`, `plume_radius`, `RadiusResult`, `classify`.
- Produces:
  - `@dataclass Plume(parent: Track, settings: LawSettings = LawSettings(), bound: str = "low", tau_pre: float = 60.0,
    r0: float | Callable | None = None)` — `r0` overrides the flow's r0 (a callable takes tau).
  - `Plume.radius(tau, t) -> RadiusResult` with `x = parent.path_length(t) - parent.path_length(tau)`.
  - `@dataclass MarginResult(margin, tau_star, distance, s_behind, r_b, regime, laws, flags, emission)` —
    `margin = r_b - distance` maximised over admissible tau; NaN with `regime` of the nearest point when no
    admissible tau has a radius (rarefied/too_slow).
  - `Plume.margin(x: ndarray(3), t: float) -> MarginResult`; `Plume.margin_bruteforce(x, t, n=20001) -> MarginResult`
    (dense grid over tau in [max(t0, t - tau_pre), t], test reference only).
  - `@dataclass ExitState(t, position, follower_speed, lag_m, distance_m, r_b, s_behind)`;
    `@dataclass ValidityRecord(interval, regimes, laws, kn, re_d, p_param, flags)` (sets/tuples over the interval);
  - `@dataclass ResidenceResult(intervals: list[tuple[float, float | None]], residence_time: float,
    total_inside: float, censored: bool, exits: list[ExitState], validity: list[ValidityRecord],
    interp_error_m: float, reason: str)` — `reason` in `{"ok", "censored", "starts_outside", "rarefied_no_law",
    "too_slow_no_law"}`; `residence_time` is the first interval's length measured from `t_start` (0 if the
    follower starts outside, NaN when the start point has no law).
  - `residence(plume, follower: Track, t_start=None, t_end=None, tol=1e-3, scan_dt=None) -> ResidenceResult`.
- Margin search (the algorithm the tests do not determine): build a `cKDTree` on `parent.samples(spacing)` with
  spacing = min(r0/2, 50 m) once per `Plume`; for a query, take the k = 16 nearest samples with
  max(t0, t - tau_pre) <= tau <= t (double k until at least one qualifies or all are tried), and for each refine
  `-(r_b(tau) - |x - P(tau)|)` with `minimize_scalar(method="bounded")` on the interval between its neighbouring
  samples, clipped to the admissible range; also evaluate tau = t exactly (the parent's own position). Return the
  best. Must agree with `margin_bruteforce` (test below).
- Residence search: sample the margin on a grid from `t_start` to `t_end` with step `scan_dt` (default: the time
  for the follower–parent relative speed at the start, or 1 m/s if smaller, to cover r0/4, clipped to [1e-4, 0.05] s),
  bracket every sign change and refine it with `brentq` to `tol`; record each crossing as entry or exit.
  If the follower is still inside at `t_end`, the last interval's end is None and `censored=True`.
  Raise `ValueError` naming both ranges when [t_start, t_end] is not inside [parent.t0, parent.t1] (Review Focus 2).

- [ ] **Step 1: Write the failing tests** in `tests/test_overlap.py` (helpers build turbulent-flow tracks so the
  radius is Taylor's unless a test says otherwise):

```python
import ast, inspect, math, numpy as np, pytest
from reentry_model.wake_residence import overlap
from reentry_model.wake_residence.overlap import Plume, residence
from reentry_model.wake_residence.track import Track
from reentry_model.wake_residence.laws import LawSettings, taylor_radius

U = 6000.0
TURB = dict(rho=1e-2, T=250.0, p=300.0, mu=1.6e-5, mfp=1e-5, d=2.813, cda=5.5, r0=1.4065)

def line_parent(T=20.0, n=2001, flow=TURB):
    t = np.linspace(0, T, n); r = np.c_[U*t, 0*t, 0*t]
    return Track.from_arrays(t, r, v=np.tile([U, 0, 0], (n, 1)), flow={k: np.full(n, v) for k, v in flow.items()})

def follower(t, r):
    return Track.from_arrays(np.asarray(t), np.asarray(r))

@pytest.mark.parametrize("modname", ["overlap", "laws", "track", "wake_gas"])
def test_engine_imports_nothing_physical(modname):
    mod = __import__("reentry_model.wake_residence." + modname, fromlist=[modname])
    tree = ast.parse(inspect.getsource(mod))
    names = {a.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))
             for a in n.names} | {n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    for bad in ("trajectory", "droplet_reentry", "parent_trajectory", "body", "aero", "atmosphere"):
        assert not any(bad in x.split(".") for x in names), (modname, bad)

def test_regime_functions_take_no_altitude():
    from reentry_model.wake_residence import laws
    for fn in (laws.classify, laws.plume_radius, overlap.Plume.margin):
        assert not any("alt" in p or p in ("h", "h_km") for p in inspect.signature(fn).parameters)

def test_fast_margin_matches_bruteforce():
    p = Plume(line_parent()); rng = np.random.default_rng(1)
    for _ in range(50):
        t = rng.uniform(2, 20); x = np.array([U*t - rng.uniform(0, 5000), rng.uniform(-30, 30), rng.uniform(-30, 30)])
        assert p.margin(x, t).margin == pytest.approx(p.margin_bruteforce(x, t).margin, abs=1e-6)

def test_straight_divergence_matches_hand_computation():
    # follower stays at the release point (x = 0) while the parent moves on; lateral offset grows at 1 m/s
    t = np.linspace(0, 20, 401); f = follower(t, np.c_[0*t, 1.0*t, 0*t])
    res = residence(Plume(line_parent()), f, t_start=0.0)
    from scipy.optimize import brentq
    # the best emission point is the release point itself (tau = 0): radius at s = U t, distance t
    expected = brentq(lambda tt: taylor_radius(U*tt, 5.5, 1.4065, 0.66, 1/3) - tt, 0.01, 19.9)
    assert res.residence_time == pytest.approx(expected, abs=1e-3)

def test_curved_parent_matches_bruteforce_exit():
    R = 6.4e6; w = U / R; t = np.linspace(0, 20, 4001)
    r = np.c_[R*np.cos(w*t), R*np.sin(w*t), 0*t]; v = np.c_[-U*np.sin(w*t), U*np.cos(w*t), 0*t]
    par = Track.from_arrays(t, r, v=v, flow={k: np.full(t.size, x) for k, x in TURB.items()})
    f = follower(t, np.c_[R*np.ones_like(t), 0*t, 2.0*t])           # stays at the start point, rises 2 m/s
    p = Plume(par); res = residence(p, f, t_start=0.0)
    te = res.intervals[0][1]
    assert p.margin_bruteforce(f.position(te), te).margin == pytest.approx(0.0, abs=1e-3)

def test_follower_starting_outside_enters_later_and_reenters():
    t = np.linspace(0, 20, 2001)
    y = 22 + 10*np.cos(0.8*t)        # 32 m at t=0 (outside), 12 m near t=4 (inside), 32 m near t=8 (outside), 12 m near t=12
    res = residence(Plume(line_parent()), follower(t, np.c_[0*t, y, 0*t]), t_start=0.0)
    assert res.residence_time == 0.0 and res.reason == "starts_outside"
    assert len(res.intervals) >= 2

def test_r0_changing_with_emission_time():
    t = np.linspace(0, 20, 401); f = follower(t, np.c_[0*t, 1.0*t, 0*t])
    park = LawSettings(turbulent_law="park")                         # Park's radius scales with r0^(2/3)
    fixed = residence(Plume(line_parent(), park, r0=1.4065), f, t_start=0.0).residence_time
    shrinking = residence(Plume(line_parent(), park, r0=lambda tau: 1.4065 * 0.5), f, t_start=0.0).residence_time
    assert shrinking < fixed

def test_subsampled_follower_gives_the_same_exit_and_a_larger_error_estimate():
    t = np.linspace(0, 20, 801); r = np.c_[0*t, 0.05*t**2 + 0.5*t, 0*t]
    fine = residence(Plume(line_parent()), follower(t, r), t_start=0.0)
    coarse = residence(Plume(line_parent()), follower(t[::2], r[::2]), t_start=0.0)
    assert abs(fine.residence_time - coarse.residence_time) < 1e-3
    assert coarse.interp_error_m >= fine.interp_error_m

def test_follower_identical_to_parent_is_censored():
    par = line_parent(); res = residence(Plume(par), par, t_start=1.0, t_end=19.0)
    assert res.censored and res.intervals[-1][1] is None and res.reason == "censored"

def test_follower_ahead_of_parent_has_small_positive_residence():           # Review Focus 4
    t = np.linspace(0, 20, 2001); f = follower(t, np.c_[1.05*U*t, 0*t, 0*t])
    res = residence(Plume(line_parent()), f, t_start=0.0)
    assert 0.0 < res.residence_time < 0.05 and res.reason == "ok"

def test_window_outside_parent_raises():                                    # Review Focus 2
    t = np.linspace(-5, 20, 101); f = follower(t, np.c_[0*t, 0*t, 0*t])
    with pytest.raises(ValueError, match="parent"):
        residence(Plume(line_parent()), f, t_start=-5.0)

def test_law_switching_along_one_track_tags_intervals():
    n = 3001; t = np.linspace(0, 30, n)
    flow = {k: np.full(n, v) for k, v in TURB.items()}
    flow["mfp"] = np.where(t < 10, 0.1, 1e-5)                         # rarefied, then continuum
    flow["rho"] = np.where(t < 20, 1e-5, 1e-2); flow["p"] = np.where(t < 20, 0.5, 300.0)   # laminar, then turbulent
    par = Track.from_arrays(t, np.c_[U*t, 0*t, 0*t], v=np.tile([U, 0, 0], (n, 1)), flow=flow)
    p = Plume(par)
    regimes = [p.margin(par.position(tt) + [0, 0.5, 0], tt).regime for tt in (5.0, 15.0, 25.0)]
    assert regimes == ["rarefied", "laminar", "turbulent"]

def test_rarefied_start_reports_no_law():
    flow = dict(TURB, mfp=0.1)
    t = np.linspace(0, 20, 401)
    res = residence(Plume(line_parent(flow=flow)), follower(t, np.c_[0*t, 0.1*t, 0*t]), t_start=0.0)
    assert math.isnan(res.residence_time) and res.reason == "rarefied_no_law"

def test_validity_record_is_populated():
    t = np.linspace(0, 20, 401)
    res = residence(Plume(line_parent()), follower(t, np.c_[0*t, 1.0*t, 0*t]), t_start=0.0)
    v = res.validity[0]
    assert v.regimes == {"turbulent"} and v.kn and v.re_d and v.p_param and "taylor" in v.laws
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_overlap.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement `overlap.py`** to the interface and algorithms above.
- [ ] **Step 4: Run them to confirm they pass.** Same command, then the whole package: `"$PY" -m pytest reentry_model/wake_residence -q`. Expected: all pass.
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/overlap.py reentry_model/wake_residence/tests/test_overlap.py
git commit -m "wake_residence: overlap engine -- union-of-balls margin with a k-d tree search checked against brute force, every inside interval, censoring and a per-interval validity record"
```

---

### Task 4: Wake gas state

**Files:**
- Create: `reentry_model/wake_residence/wake_gas.py`, `reentry_model/wake_residence/tests/test_wake_gas.py`
- Modify: `reentry_model/wake_residence/overlap.py` (add `Plume.wake_state`)

**Interfaces:**
- Consumes: `EmissionState`, `MarginResult`.
- Produces:
  - `@dataclass(frozen=True) Ambient(p, T, rho)`; `@dataclass(frozen=True) GasState(rho, T, p, u: ndarray(3))`.
  - `CP_AIR = 1004.5`, `R_AIR = 287.05`, `WAKE_FLOW_MODES = ("off", "velocity", "density", "both")`.
  - `gas_state(r_axis, r_b, emission, direction, ambient, mode, k_E=1.0, cutoff=6.0) -> GasState`.
  - `Plume.wake_state(x, t, ambient, mode, k_E=1.0) -> tuple[GasState, MarginResult]` — direction = unit parent
    velocity at `tau_star`; `r_axis` = `MarginResult.distance`; ambient air when the margin is NaN or mode is off.
- Profile algorithm (§4.3a, the call recorded in spec §8):

```
D    = 0.5 * e.rho * e.U**2 * e.cda                       # drag at emission [N]
R_E  = max(r_b, k_E * sqrt(D / e.p))                      # energy radius: core or line-blast radius
T_st = ambient.T + e.U**2 / (2 * CP_AIR)                  # stagnation-temperature cap
q_c  = min(D * R_AIR / (ambient.p * CP_AIR * pi * R_E**2), 1 - ambient.T / T_st)
q    = q_c * exp(-(r_axis / R_E)**2)                      # q = 1 - T_inf/T, Gaussian; integral = D R/(p c_p)
T    = ambient.T / (1 - q);  rho = ambient.p / (R_AIR * T)
m    = (D / e.U) / (pi * r_b**2) * exp(-(r_axis / r_b)**2) # forward momentum density, integral = D/U
u    = min(m / rho_used, e.U) * direction                 # rho_used = rho (both) or ambient.rho (velocity)
mode "density": u = 0; mode "velocity": T, rho ambient; beyond cutoff * max(r_b, R_E): exactly ambient
```

- [ ] **Step 1: Write the failing tests** in `tests/test_wake_gas.py`:

```python
import math, numpy as np, pytest
from scipy.integrate import quad
from reentry_model.wake_residence.wake_gas import gas_state, Ambient, CP_AIR, R_AIR
from reentry_model.wake_residence.laws import EmissionState, sutherland_mu

E = EmissionState(rho=3e-4, T=240.0, p=20.0, mu=sutherland_mu(240.0), mfp=1e-4, U=4000.0, d=2.813, cda=5.5, r0=1.4065)
AMB = Ambient(p=20.0, T=240.0, rho=20.0 / (R_AIR * 240.0)); X = np.array([1.0, 0, 0])

def test_momentum_and_energy_are_conserved():
    r_b = 40.0; D = 0.5 * E.rho * E.U**2 * E.cda
    mom = quad(lambda r: gas_state(r, r_b, E, X, AMB, "both").rho * gas_state(r, r_b, E, X, AMB, "both").u[0]
               * 2 * math.pi * r, 0, 10 * r_b, limit=200)[0]
    q = lambda r: 1 - AMB.T / gas_state(r, r_b, E, X, AMB, "both").T
    energy = quad(lambda r: q(r) * 2 * math.pi * r, 0, 10 * max(r_b, math.sqrt(D / E.p)), limit=200)[0]
    assert mom == pytest.approx(D / E.U, rel=1e-3)           # only exact where neither cap binds (r_b chosen so)
    assert energy == pytest.approx(D * R_AIR / (AMB.p * CP_AIR), rel=1e-3)

def test_caps_hold_close_to_a_fast_parent():
    s = gas_state(0.0, 0.5, EmissionState(**{**E.__dict__, "U": 7800.0, "rho": 1e-5, "p": 0.6}), X,
                  Ambient(p=0.6, T=200.0, rho=0.6 / (R_AIR * 200.0)), "both")
    assert s.u[0] <= 7800.0 and s.T <= 200.0 + 7800.0**2 / (2 * CP_AIR) + 1e-9

def test_state_is_smooth_across_the_plume_edge():
    r_b = 40.0; f = lambda r: gas_state(r, r_b, E, X, AMB, "both")
    for g in (lambda s: s.T, lambda s: s.u[0]):
        d1 = (g(f(r_b + 1e-3)) - g(f(r_b - 1e-3))) / 2e-3
        d2 = (g(f(r_b + 2e-3)) - g(f(r_b - 2e-3))) / 4e-3
        assert d1 == pytest.approx(d2, rel=1e-3)

def test_off_and_far_field_are_exactly_ambient():
    for s in (gas_state(1.0, 40.0, E, X, AMB, "off"), gas_state(1e5, 40.0, E, X, AMB, "both")):
        assert (s.rho, s.T, s.p) == (AMB.rho, AMB.T, AMB.p) and not s.u.any()

def test_modes_switch_one_effect_each():
    v, d = gas_state(5.0, 40.0, E, X, AMB, "velocity"), gas_state(5.0, 40.0, E, X, AMB, "density")
    assert v.T == AMB.T and v.u[0] > 0 and d.T > AMB.T and not d.u.any()
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_wake_gas.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement `wake_gas.py` and `Plume.wake_state`.** If the conservation test needs r_b larger so no cap
  binds, change the test's r_b, not the tolerance, and say so in the commit message.
- [ ] **Step 4: Run them to confirm they pass**, then the package: `"$PY" -m pytest reentry_model/wake_residence -q`.
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/wake_gas.py reentry_model/wake_residence/overlap.py \
  reentry_model/wake_residence/tests/test_wake_gas.py
git commit -m "wake_residence: wake gas state -- forward gas velocity and hot, thin core from the parent's momentum and drag work, Gaussian profiles, caps and the four wake_flow modes"
```

---

### Task 5: Twin parent from the DRAMA file

**Files:**
- Create: `reentry_model/wake_residence/adapters.py` (first part), `reentry_model/wake_residence/tests/test_twin.py`
- Modify: `reentry_model/wake_residence/tests/conftest.py` (session fixture `twin`)

**Interfaces:**
- Consumes: `parent_trajectory.load_parent`, `ParentTrajectory` (`.rows`, `.branch`, `.path`, `.epoch`,
  `.sha256`), `model.parent_release`, `trajectory.Simulator`, `trajectory.Settings`, `trajectory.InitialState`,
  `earth.ecef_to_geodetic`, `earth.flight_angles(v, lat, lon) -> (V, gamma, heading)`, `aero.mean_free_path`,
  `particle_trajectory.viscosity` (Sutherland), `atmosphere.US76TableAtmosphere`, `constants.G0`, `Track`.
- Produces (in `adapters.py`):
  - `@dataclass BetaTable(h_m, beta, mass)` with methods `beta_at(h)`, `mass_at(h)`, `cda_at(h)` (= mass/beta),
    linear in altitude, clamped at the ends; `drama_beta_table(parent, smooth_window=5) -> BetaTable` — beta =
    `dynamicPressure / (loadFactor * G0)` on the descending branch, rows with loadFactor below 1e-4 dropped,
    smoothed with a centred running median of `smooth_window` rows (document the choice); mass from the same
    object's `_AeroThermalHistory.txt` table (read from `parent.path`), interpolated in altitude.
  - `class TwinSimulator(Simulator)`: `__init__(initial, beta_table, atmosphere, settings, beta_scale=1.0)`;
    overrides `aero_state` to return the parent `AeroState` with `a_drag = -0.5 * rho * V * v_rel / (beta_scale * beta_at(h))`
    (call `super().aero_state` for everything else; it needs a body for `mass(t)` — pass a `ConstantBody(1.0)`
    and a dummy diameter, since only the overridden drag is used).
  - `dense_solve(sim, t_max, events=None, method="DOP853", rtol=None, atol=None) -> DenseRun(sol, t_end, end_reason)`
    — `solve_ivp(sim.rhs, (0, t_max), sim.initial_state_vector(), dense_output=True, events=events or sim._events())`
    with the simulator's tolerances by default.
  - `flow_arrays(t, states, atmosphere, d, cda, r0) -> dict` — per sample rho, T, p, mu (Sutherland), mfp
    (`aero.mean_free_path(rho, m_bar)`), d, cda (scalar or array), r0; plus `speed` (air-relative).
  - `@dataclass TwinParent(track: Track, run: DenseRun, beta: BetaTable, epoch, sha256)` with
    `t_at_altitude(h_km) -> float` (root of altitude(t) - h on the descending track; raises `ValueError` naming the
    track's altitude range if h is outside it — Review Focus 5).
  - `twin_parent(parent=None, h_stop_km=45.0, atmosphere=None, frontal_area=6.215, beta_scale=1.0, cadence=0.02)
    -> TwinParent` — starts from `parent_release(parent, branch-top altitude)`, integrates to `h_stop_km`,
    samples the dense output every `cadence` s into `Track.from_arrays(t, r, v, speed, flow)`.
  - `release_state(track, t, epoch) -> InitialState` — position/velocity of the track at t, converted with
    `ecef_to_geodetic` and `flight_angles`.
  - `twin_vs_drama(twin, parent) -> list[dict]` — one dict per DRAMA row on the branch with keys `time_s`,
    `drama_altitude_km`, `drama_speed_ms`, `drama_downrange_km`, `twin_speed_ms` (twin speed at the DRAMA
    altitude), `twin_altitude_km` (twin altitude at the DRAMA time), `speed_rel_diff`
    (`twin_speed_ms / drama_speed_ms - 1`).

- [ ] **Step 1: Write the failing tests** in `tests/test_twin.py` (`twin` is a session fixture `twin_parent()`):

```python
import numpy as np, pytest
from reentry_model.wake_residence.adapters import drama_beta_table, release_state, twin_vs_drama
from reentry_model.wake_residence.overlap import Plume

def test_beta_table_matches_builder_check(parent):
    b = drama_beta_table(parent)
    assert 120 < b.beta_at(85e3) < 160 and 20 < b.beta_at(60e3) < 40
    assert b.mass_at(90e3) == pytest.approx(740.0, rel=1e-3) and b.mass_at(66e3) < 300

def test_twin_reproduces_drama_speed_within_one_percent_above_60km(twin, parent):
    rows = [r for r in twin_vs_drama(twin, parent) if r["drama_altitude_km"] > 60.0]
    worst = max(abs(r["speed_rel_diff"]) for r in rows)
    print("twin vs DRAMA worst speed difference above 60 km:", worst)
    assert worst <= 0.01

def test_release_starts_exactly_on_the_track(twin):
    t_r = twin.t_at_altitude(80.0)
    st = release_state(twin.track, t_r, twin.epoch)
    from reentry_model import earth
    x = earth.geodetic_to_ecef(st.altitude, st.lat, st.lon)
    m = Plume(twin.track).margin(x, t_r)
    assert m.margin == pytest.approx(m.r_b, abs=1e-3) and m.distance < 1e-3

def test_altitude_outside_track_raises(twin):
    with pytest.raises(ValueError, match="km"):
        twin.t_at_altitude(130.0)
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_twin.py -q -s`. Expected: ImportError.
- [ ] **Step 3: Implement** the adapters above and the `twin` session fixture.
- [ ] **Step 4: Run them.** Same command. Expected: pass, and the printed worst difference is at most 0.01.
  **If the 1 percent test fails**: do not loosen it. Add `DramaColumnAtmosphere(parent)` to `adapters.py`
  (density from the file's `density` column, temperature from `(velocity / mach)^2 * m_bar / (1.4 k_B)`, pressure
  from the ideal gas, interpolated in altitude, no wind), add `twin_parent(..., atmosphere=DramaColumnAtmosphere(parent))`
  as a second variant, rerun, and report both figures to the controller; the controller decides with Asha which
  atmosphere the twin and the followers share (spec §4.3 says a miss is reported, not hidden).
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/adapters.py reentry_model/wake_residence/tests/test_twin.py \
  reentry_model/wake_residence/tests/conftest.py
git commit -m "wake_residence: twin of the DRAMA parent integrated by the repository's own integrator with DRAMA's altitude-dependent ballistic coefficient; reproduces DRAMA's speed to <worst> above 60 km"
```
(Replace `<worst>` with the measured figure.)

---

### Task 6: Followers in the wake (WakeAtmosphere, sphere follower, sphere parent)

**Files:**
- Modify: `reentry_model/wake_residence/adapters.py`
- Create: `reentry_model/wake_residence/tests/test_follower.py`

**Interfaces:**
- Consumes: `Plume.wake_state`, `Ambient`, `atmosphere.Freestream`, `earth.enu_basis`, `earth.geodetic_to_ecef`,
  `body.ConstantBody`, `aero.SphereDragTables.from_json`, `aero.bridging_by_name`, `aero.drag_coefficient`,
  `aero.knudsen`, `aero.mach`, `trajectory.Simulator`, `trajectory.Settings`, `dense_solve`, `flow_arrays`.
- Produces:
  - `class WakeAtmosphere(base, plume, t_offset, mode, k_E=1.0)` with `state(t, h, lat, lon) -> Freestream`: query
    `plume.wake_state(geodetic_to_ecef(h, lat, lon), t + t_offset, Ambient(base.p, base.T, base.rho), mode)`;
    return the base `Freestream` object itself when mode is `"off"` or the point is beyond the cutoff; otherwise
    `Freestream(rho, T, p, base.m_bar, base.wind_enu + ENU(u))`. Counts its calls (`n_calls`) for the timing step.
  - `sphere_simulator(initial, d_m, atmosphere, density=2813.0, bridging="sesam-table", t_max=600.0) -> Simulator`
    — exactly as `cli.py` builds a `--thermal none` run (lines 207-214): `ConstantBody(body.sphere_mass(d_m, density))`,
    `SphereDragTables.from_json()`, `bridging_by_name(bridging)`, `Settings(diameter=d_m, gravity="j2",
    rotating_frame=True, t_max=t_max)`.
  - `solution_track(run, t_offset, flow=None, cadence=None) -> Track` — samples the dense solution (default cadence
    0.002 s, plus the end time) and shifts times by `t_offset` onto the parent clock.
  - `sphere_parent_track(run, sim, t_offset, atmosphere) -> Track` — like `solution_track`, with flow arrays whose
    `d` is the sphere diameter, `r0` = d/2 and `cda` = C_D(Kn, Ma) × area from the simulator's own tables.

- [ ] **Step 1: Write the failing tests** in `tests/test_follower.py`:

```python
import time, numpy as np, pytest
from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.wake_residence.adapters import (WakeAtmosphere, sphere_simulator, dense_solve, release_state)
from reentry_model.wake_residence.overlap import Plume
from reentry_model.wake_residence.laws import LawSettings

def release(twin, h_km=60.0):
    t_r = twin.t_at_altitude(h_km); return t_r, release_state(twin.track, t_r, twin.epoch)

def test_off_returns_the_wrapped_freestream(twin):
    base = US76TableAtmosphere(); w = WakeAtmosphere(base, Plume(twin.track), 0.0, "off")
    a, b = w.state(1.0, 70e3, 0.5, 0.2), base.state(1.0, 70e3, 0.5, 0.2)
    assert (a.rho, a.T, a.p, a.m_bar) == (b.rho, b.T, b.p, b.m_bar) and np.array_equal(a.wind_enu, b.wind_enu)

def test_off_run_reproduces_plain_simulator(twin):
    t_r, st = release(twin)
    plain = dense_solve(sphere_simulator(st, 0.05, US76TableAtmosphere()), 5.0)
    wrapped = dense_solve(sphere_simulator(st, 0.05, WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, "off")), 5.0)
    assert np.allclose(plain.sol(5.0), wrapped.sol(5.0), rtol=0, atol=1e-9)

def test_tailwind_and_thin_gas_reduce_deceleration(twin):
    t_r, st = release(twin)
    speeds = {}
    for mode in ("off", "both"):
        run = dense_solve(sphere_simulator(st, 0.05, WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, mode)), 1.0)
        speeds[mode] = np.linalg.norm(run.sol(1.0)[3:])
    assert speeds["both"] > speeds["off"]

def test_timing_of_one_wake_run_is_recorded(twin, capsys):
    t_r, st = release(twin)
    w = WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, "both")
    t0 = time.perf_counter(); dense_solve(sphere_simulator(st, 0.05, w), 2.0)
    print("wake run: {:.1f} s for 2 s of flight, {} atmosphere calls".format(time.perf_counter() - t0, w.n_calls))
    assert w.n_calls > 0
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_follower.py -q -s`. Expected: ImportError.
- [ ] **Step 3: Implement** the adapters above.
- [ ] **Step 4: Run them.** Same command. Expected: pass; report the printed timing to the controller. If one
  second of flight costs more than about 30 s, profile `Plume.margin` and speed it up (cache the k-d tree query
  across calls at nearly the same point, restrict the candidate tau window) before moving on.
- [ ] **Step 5: Confirm the integrator is still the verified one.** Run the existing SESAM comparison:
  `"$PY" -m pytest tests/test_reentry_model_reference.py -m reference -q -k "<one fast case>"` (pick the shortest
  case listed in that file's parametrisation; its thresholds `dV_rel_max`, `dh_max_m`, `d_end_time_rel` are the
  repository's). Expected: pass. Record the case name in the commit message.
- [ ] **Step 6: Commit**

```bash
git add reentry_model/wake_residence/adapters.py reentry_model/wake_residence/tests/test_follower.py
git commit -m "wake_residence: WakeAtmosphere feeds the wake gas to any follower through the atmosphere's state() call; sphere followers and sphere parents built as cli.py builds a --thermal none run"
```

---

### Task 7: Release classifier (Register et al. 2020, AFFIx)

**Files:**
- Create: `reentry_model/wake_residence/affix.py`, `reentry_model/wake_residence/tests/test_affix.py`
- Modify: `reentry_model/wake_residence/tests/test_overlap.py` (add `affix` to the import-purity parametrisation)

**Interfaces:**
- Consumes: `laws.EmissionState`, `laws.LawSettings.kn_continuum` (pass the value in, do not import settings
  objects you do not need).
- Produces (in `affix.py`, a pure module: no imports of trajectory/droplet_reentry/parent_trajectory/body/aero/
  atmosphere, no altitude argument):
  - `R_STAR = 0.353`, `R_FIT_RANGE = (1/16, 1.0)`, `GAMMA_AIR = 1.4`, `R_AIR = 287.05`, zone names
    `ZONES = ("trapped", "far_wake", "shock_band", "escape")`.
  - `theta_inn(r) -> float` = 122 r + 47 (degrees) for r <= R_STAR, NaN otherwise.
  - `theta_bif(r) -> float` = 100 exp(-1.618 r) + 40 (degrees).
  - `@dataclass(frozen=True) AffixSettings(s_w_deg=2.0, kn_continuum=0.01, mach_min=5.0)`.
  - `@dataclass(frozen=True) AffixResult(zone, r, theta_deg, theta_inn, theta_bif, t_hat, t_fz, flags: dict)` —
    flags keys: `r_below_fit`, `r_above_rstar`, `rarefied_parent`, `rarefied_follower`, `not_hypersonic`,
    `unequal_density`; plus `parent_bulk_density` (NaN when no parent mass is given).
  - `affix_regime(d_follower, d_parent, theta_deg, emission, follower_density, parent_mass=None,
    settings=AffixSettings()) -> AffixResult` — zone: trapped if theta < theta_inn; far_wake if theta_inn <= theta
    < theta_bif - s_w; shock_band if |theta - theta_bif| <= s_w; escape if theta > theta_bif + s_w. For r > R_STAR
    theta_inn is NaN, `r_above_rstar` is set, and the zone is decided by theta_bif alone (far_wake below the band).
    t_hat: 0 (trapped); A0(r) + 0.5 exp[-0.05(theta - theta_inn)] with A0 = -7.2 r^2 + 3 r + 0.1 (far_wake; NaN when
    r > R_STAR); 0.88 r + 0.32 + 0.01 exp[0.125(theta - theta_bif - 10)] (escape); NaN (shock_band).
    t_fz = d_parent * sqrt(follower_density / q_inf) * t_hat with q_inf = 0.5 * emission.rho * emission.U**2.
    Flags: r < 1/16; mfp/d > kn_continuum for parent and for follower; Mach = U / sqrt(GAMMA_AIR R_AIR T) below
    mach_min; unequal_density when parent_mass is given and m/(pi d_parent^3/6) differs from follower_density by more
    than 1 percent. Raises ValueError for theta outside [0, 180] or non-positive diameters.
  - `release_offset(theta_deg, d_parent, v_hat, up) -> ndarray(3)` — the surface point relative to the centroid:
    (d_parent/2)(-cos(theta) v_hat + sin(theta) n_hat), n_hat = unit(up - (up·v_hat) v_hat) (raise if up is
    parallel to v_hat).
  - `far_wake_offset(d_parent, v_hat) -> ndarray(3)` = -2 d_parent v_hat (Asha: two parent diameters behind, on the
    axis).
- Source notes: `/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/citations/fragment_shock_separation.md` §00 (main checkout only; cite in the docstring).

- [ ] **Step 1: Write the failing tests** in `tests/test_affix.py`:

```python
import math, numpy as np, pytest
from reentry_model.wake_residence.affix import (affix_regime, theta_inn, theta_bif, AffixSettings,
    release_offset, far_wake_offset)
from reentry_model.wake_residence.laws import EmissionState, sutherland_mu

E = EmissionState(rho=1.85e-5, T=199.0, p=1.05, mu=sutherland_mu(199.0), mfp=4e-3, U=7600.0, d=2.813,
                  cda=5.5, r0=1.4065)          # roughly the Generic Spacecraft at 80 km

@pytest.mark.parametrize("r", [1/16, 1/4])
def test_thresholds_follow_register(r):
    assert theta_inn(r) == pytest.approx(122 * r + 47)
    assert theta_bif(r) == pytest.approx(100 * math.exp(-1.618 * r) + 40)

def test_r_above_rstar_has_no_inner_threshold():
    assert math.isnan(theta_inn(0.5)) and theta_bif(0.5) == pytest.approx(100 * math.exp(-0.809) + 40)

def test_windward_small_sphere_far_wake_and_stagnation_escape():
    d2 = 0.01
    assert affix_regime(d2, E.d, 90.0, E, 2813.0).zone == "far_wake"
    assert affix_regime(d2, E.d, 180.0, E, 2813.0).zone == "escape"
    tb = theta_bif(d2 / E.d)
    assert affix_regime(d2, E.d, tb + 0.5, E, 2813.0).zone == "shock_band"
    assert affix_regime(d2, E.d, 30.0, E, 2813.0).zone == "trapped"

def test_interaction_times():
    d2 = 0.176; r = d2 / E.d
    q = 0.5 * E.rho * E.U**2; scale = E.d * math.sqrt(2813.0 / q)
    esc = affix_regime(d2, E.d, 180.0, E, 2813.0)
    assert esc.t_fz == pytest.approx(scale * (0.88*r + 0.32 + 0.01*math.exp(0.125*(180 - theta_bif(r) - 10))))
    fw = affix_regime(d2, E.d, 90.0, E, 2813.0)
    A0 = -7.2*r**2 + 3*r + 0.1
    assert fw.t_fz == pytest.approx(scale * (A0 + 0.5*math.exp(-0.05*(90 - theta_inn(r)))))
    assert affix_regime(d2, E.d, 10.0, E, 2813.0).t_fz == 0.0

def test_validity_flags():
    f = affix_regime(0.01, E.d, 90.0, E, 2813.0, parent_mass=740.0).flags
    assert f["r_below_fit"] and f["rarefied_follower"] and not f["rarefied_parent"] and f["unequal_density"]
    slow = EmissionState(**{**E.__dict__, "U": 800.0})
    assert affix_regime(0.3, E.d, 90.0, slow, 2813.0).flags["not_hypersonic"]

def test_theta_range_and_no_altitude_argument():
    import inspect
    with pytest.raises(ValueError):
        affix_regime(0.01, E.d, 181.0, E, 2813.0)
    assert not any("alt" in p or p in ("h", "h_km") for p in inspect.signature(affix_regime).parameters)

def test_release_offsets():
    v = np.array([1.0, 0, 0]); up = np.array([0, 0, 1.0])
    assert np.allclose(far_wake_offset(2.813, v), [-5.626, 0, 0])
    assert np.allclose(release_offset(180.0, 2.813, v, up), [1.4065, 0, 0])
    assert np.allclose(release_offset(90.0, 2.813, v, up), [0, 0, 1.4065], atol=1e-12)
    assert np.linalg.norm(release_offset(135.0, 2.813, v, up)) == pytest.approx(1.4065)
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_affix.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement `affix.py`** and add `"affix"` to the import-purity parametrisation in `test_overlap.py`.
- [ ] **Step 4: Run them to confirm they pass**, then `"$PY" -m pytest reentry_model/wake_residence -q`.
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/affix.py reentry_model/wake_residence/tests/test_affix.py \
  reentry_model/wake_residence/tests/test_overlap.py
git commit -m "wake_residence: Register et al. (2020) release classifier -- trapped, far wake, shock band or escape from the size ratio and release angle, with the interaction time and a validity record"
```

---

### Task 7b: Shock-layer field and the rarefied-follower simulator (§4.3d)

Added 2026-10-29 after Asha chose option 1 (an engineering shock-layer model). Read spec §4.3c and §4.3d first: they
fix every formula below, including Lighthill's conditions and the sign of the r^4 coefficient.

**Files:**
- Create: `reentry_model/wake_residence/shock_layer.py`, `reentry_model/wake_residence/tests/test_shock_layer.py`,
  `reentry_model/wake_residence/tests/test_rarefied_follower.py`
- Modify: `reentry_model/wake_residence/adapters.py` (WakeAtmosphere option, pressure gradient, rarefied simulator,
  collision event, exit classification)

**Interfaces:**
- Consumes: `gas.EquilibriumAir` (its `gas` Cantera object, `stagnation`, `_state`-style GasState), `aero.mean_free_path`,
  `droplet_reentry.particle_trajectory.drag_coefficient(V, rho, T, p, d, temperature_ratio, law="blend")`,
  `overlap.Plume`, `Track.position/velocity`, `WakeAtmosphere` (Task 6), `affix.release_offset`, `earth.enu_basis`,
  `earth.geodetic_to_ecef`, `trajectory.Simulator`, `AeroState`.
- Produces in `shock_layer.py` (may import `reentry_model.gas` and `aero`; no altitude constants or arguments; no
  imports of `trajectory`, `droplet_reentry`, `atmosphere`):
  - `@dataclass(frozen=True) ShockLayerSettings(shock_standoff="lighthill", phi_nose_deg=40.0, phi_blend_deg=20.0,
    shock_thickness_mfp=3.0, n_phi=91, n_eta=21, pressure_gradient=True, follower_wall_T_K=300.0,
    kn_shock_layer_flag=1e-3, q_change_flag=0.10)` with `non_default() -> dict` (as `LawSettings`), validating
    `shock_standoff in ("lighthill", "billig")` and positive numbers.
  - `@dataclass(frozen=True) FreeStream(rho, T, p, U, mfp, m_bar)` (the parent's free stream at release).
  - `@dataclass(frozen=True) PostShock(rho, T, p, h, s, m_bar, w_n, w_t)` and
    `oblique_jump(air, fs, sigma_rad) -> PostShock` — equilibrium Rankine-Hugoniot on the normal component
    `fs.U sin(sigma)` (same damped fixed point as `EquilibriumAir.stagnation`, same `T_EQUILIBRATE` gate),
    `w_t = fs.U cos(sigma)`; the post-shock velocities are in the shock-normal/tangent frame relative to the parent.
  - `@dataclass(frozen=True) Lighthill(eps, U, R, R_s, A, B, K)` with `velocity(r, phi) -> (u_r, u_phi)`;
    `lighthill(eps, U, R) -> Lighthill` (root of g(R) = 0 closest to the body, K = -U (1-eps)^2 / (10 eps R_s^2));
    `lighthill_standoff_approx(eps) -> float` = eps / (1 + sqrt(8 eps / 3)) (as a fraction of R_s).
  - `billig_standoff(R, M)`, `billig_vertex_radius(R, M)`; `@dataclass(frozen=True) ShockShape(R, M, standoff)` with
    `x_s(y)`, `sigma(y)` (local shock angle, radians), `r_s(phi)` (distance from the centroid to the shock along the
    ray at phi; root-find), `inside(x, y) -> bool` and `normal_distance(x, y) -> float` (signed, positive inside,
    first order along the local shock normal).
  - `@dataclass FieldState(rho, T, p, m_bar, w_x, w_y, weight)` — `w_x`, `w_y` the gas velocity relative to the parent
    (x forward along v_hat, y away from the axis), `weight` the field's share in [0, 1] (1 inside the layer, 0 beyond
    the smoothing outside the shock and beyond r0 behind the base plane).
  - `class ShockLayerField` with `ShockLayerField.build(fs: FreeStream, R, r0, settings, air=None)` (computes epsilon
    from the normal jump, the standoff per `settings.shock_standoff`, the shape, the surface states, the (phi, eta) grid
    and bicubic splines of p, T, rho, m_bar, w_x, w_y), attributes `eps`, `standoff`, `standoff_billig`, `shape`,
    `fs`, `build_time_s`; methods `state(x, y) -> FieldState` and `excess_pressure_gradient(x, y) -> (dpx, dpy)` =
    the gradient of weight * (p_field - fs.p) (so it includes the jump across the smoothed shock and is exactly zero
    outside the field); `surface_state(theta_deg) -> FieldState` (eta = 0 at phi = 180 - theta).
- Produces in `adapters.py`:
  - `WakeAtmosphere(base, plume, t_offset, mode, k_E=1.0, shock_layer=None)`: with a field, `state()` blends the field
    (weight w) with the wake-or-ambient gas (1 - w): rho, T, p and m_bar linearly, and the Earth-fixed velocity
    w * (parent velocity at the parent time + field velocity rotated into ECEF) + (1 - w) * wake velocity. Mode `"off"`
    with a field still applies the field. Without a field, behaviour is unchanged (the Task 6 tests must still pass).
  - `WakeAtmosphere.parent_frame(t, x_ecef) -> (x, y, e_x, e_y)` (e_y the unit lateral vector; any perpendicular on
    the axis) and `WakeAtmosphere.pressure_gradient(t, x_ecef) -> ndarray(3)` (ECEF; zeros without a field).
  - `class RarefiedSphereSimulator(Simulator)` and `rarefied_sphere_simulator(initial, d_m, atmosphere,
    density=AA7075_DENSITY, wall_T_K=300.0, pressure_gradient=True, t_max=600.0)`: `aero_state` replaces the drag
    with -(1/2) rho V v_rel C_D A / m, C_D from `particle_trajectory.drag_coefficient` with the local gas and
    `temperature_ratio = wall_T_K / T`, the force taken as zero when V == 0 (no NaN), and adds
    `-atmosphere.pressure_gradient(t, r) / density` when `pressure_gradient` (only if the atmosphere has that method).
    Records the largest |a_p| and |a_drag| seen inside the field (weight > 0.5) as attributes.
  - `parent_collision_event(wake_atmosphere, R)` -> solve_ivp event, distance from the parent's centroid minus R,
    `terminal = True`, `direction = -1`.
  - `@dataclass ShockLayerExit(kind, t, x, y, v_rel)` and `shock_layer_exit(run, wake_atmosphere, field, hit=False,
    tol=1e-3) -> ShockLayerExit` — kind `hit_parent` if the collision event ended the run inside the envelope, else
    the first of `through_shock` (x >= x_s(y)) or `past_base_plane` (x < 0 inside the envelope) found on the dense
    output to `tol` seconds, else `inside_at_end`.

- [ ] **Step 1: Write the failing tests** in `tests/test_shock_layer.py`:

```python
import math, numpy as np, pytest
from reentry_model.gas import EquilibriumAir
from reentry_model.wake_residence.shock_layer import (ShockLayerSettings, FreeStream, oblique_jump, lighthill,
    lighthill_standoff_approx, billig_standoff, billig_vertex_radius, ShockShape, ShockLayerField)

R, R0 = 1.4065, 1.4065
FS = FreeStream(rho=1.85e-5, T=199.0, p=1.05, U=7400.0, mfp=4e-3, m_bar=4.81e-26)   # about 80 km
M = FS.U / math.sqrt(1.4 * 287.05 * FS.T)

@pytest.fixture(scope="module")
def air():
    return EquilibriumAir()

@pytest.fixture(scope="module")
def field(air):
    return ShockLayerField.build(FS, R, R0, ShockLayerSettings(), air=air)

def test_billig_formulas():
    assert billig_standoff(R, M) == pytest.approx(R * 0.143 * math.exp(3.24 / M**2))
    assert billig_vertex_radius(R, M) == pytest.approx(R * 1.143 * math.exp(0.54 / (M - 1) ** 1.2))
    s = ShockShape(R, M, billig_standoff(R, M))
    assert s.x_s(0.0) == pytest.approx(R + billig_standoff(R, M))
    for phi in np.radians([0, 45, 89]):
        assert s.inside(R * math.cos(phi), R * math.sin(phi))

def test_normal_jump_then_compression_gives_stagnation(air):
    stag = air.stagnation(FS.rho, FS.T, FS.U)
    ps = oblique_jump(air, FS, math.pi / 2)
    assert ps.w_t == pytest.approx(0.0, abs=1e-6)
    assert ps.p < stag.p and ps.p == pytest.approx(stag.p, rel=0.05)     # compression to rest adds a few percent
    from reentry_model.gas import AIR
    from scipy.optimize import brentq
    h0 = ps.h + 0.5 * ps.w_n**2
    def res(lp):
        air.gas.SPX = ps.s, math.exp(lp), AIR; air.gas.equilibrate("SP")
        return air.gas.enthalpy_mass - h0
    p0 = math.exp(brentq(res, math.log(ps.p), math.log(2 * ps.p)))
    assert p0 == pytest.approx(stag.p, rel=0.005)

def test_oblique_tangential_preserved_and_weak_limit(air):
    beta = math.asin(1 / M)
    ps = oblique_jump(air, FS, math.radians(60))
    assert ps.w_t == pytest.approx(FS.U * math.cos(math.radians(60)))
    weak = oblique_jump(air, FS, beta + math.radians(0.5))
    assert weak.rho / FS.rho == pytest.approx(1.0, rel=0.05)

@pytest.mark.parametrize("eps", [0.05, 0.1, 0.2])
def test_lighthill_conditions_and_standoff(eps):
    L = lighthill(eps, 1.0, 1.0)
    g = lambda r: L.A * r**2 + L.B / r + L.K * r**4
    dg = lambda r: 2 * L.A * r - L.B / r**2 + 4 * L.K * r**3
    assert g(1.0) == pytest.approx(0.0, abs=1e-12)
    assert g(L.R_s) == pytest.approx(-0.5 * eps * L.R_s**2, rel=1e-10)
    assert dg(L.R_s) == pytest.approx(-L.R_s, rel=1e-10)
    assert (L.R_s - 1.0) / L.R_s == pytest.approx(lighthill_standoff_approx(eps), rel=0.10)
    ur, uphi = L.velocity(1.0, math.radians(30))
    assert ur == pytest.approx(0.0, abs=1e-12)                       # body is a streamline

def test_field_stagnation_body_shock(field, air):
    stag = air.stagnation(FS.rho, FS.T, FS.U)
    s = field.state(R + 1e-9, 0.0)
    assert s.p == pytest.approx(stag.p, rel=0.01) and s.T == pytest.approx(stag.T, rel=0.01)
    assert abs(s.w_x) < 1.0 and abs(s.w_y) < 1.0 and s.weight == pytest.approx(1.0)
    phi = math.radians(50); ex, ey = math.cos(phi), math.sin(phi)
    b = field.state(R * ex * (1 + 1e-9), R * ey * (1 + 1e-9))
    assert abs(b.w_x * ex + b.w_y * ey) < 1e-3 * FS.U                 # no flow through the body
    stag_p = stag.p
    assert b.p == pytest.approx(FS.p + (stag_p - FS.p) * math.cos(phi)**2, rel=0.01)

def test_mass_flux_through_layer(field):
    for phi_deg in (30, 60):
        phi = math.radians(phi_deg); rs = field.shape.r_s(phi)
        r = np.linspace(R, rs, 400); flux = 0.0
        for a, b in zip(r[:-1], r[1:]):
            rm = 0.5 * (a + b); st = field.state(rm * math.cos(phi), rm * math.sin(phi))
            u_phi = -st.w_x * math.sin(phi) + st.w_y * math.cos(phi)
            flux += st.rho * abs(u_phi) * 2 * math.pi * rm * math.sin(phi) * (b - a)
        ys = rs * math.sin(phi)
        assert flux == pytest.approx(math.pi * ys**2 * FS.rho * FS.U, rel=0.25)

def test_outside_is_ambient_and_edges_are_continuous(field):
    x_far = field.shape.x_s(0.0) + 20 * FS.mfp
    assert field.state(x_far, 0.0).weight == 0.0
    assert field.state(-1.01 * R0, 0.5 * R).weight == 0.0
    ws = [field.state(x, 1.2 * R).weight for x in np.linspace(0.05, -0.05, 41)]   # across the base plane
    assert max(abs(a - b) for a, b in zip(ws[:-1], ws[1:])) < 0.05
    xs = field.shape.x_s(0.0); d = FS.mfp * 3
    ws = [field.state(x, 0.0).weight for x in np.linspace(xs - 2 * d, xs + 2 * d, 81)]
    assert max(abs(a - b) for a, b in zip(ws[:-1], ws[1:])) < 0.1

def test_pressure_gradient_matches_finite_difference_and_points_to_body(field):
    x, y, h = R + 0.3 * field.standoff, 0.2, 1e-5
    gx, gy = field.excess_pressure_gradient(x, y)
    pe = lambda x, y: (lambda s: s.weight * (s.p - FS.p))(field.state(x, y))
    assert gx == pytest.approx((pe(x + h, y) - pe(x - h, y)) / (2 * h), rel=1e-3)
    assert gy == pytest.approx((pe(x, y + h) - pe(x, y - h)) / (2 * h), rel=1e-3, abs=1e-6)
    assert field.excess_pressure_gradient(R + 0.5 * field.standoff, 0.0)[0] < 0   # pressure rises towards the body

def test_no_altitude_argument():
    import inspect
    assert not any("alt" in p or p in ("h", "h_km") for p in inspect.signature(ShockLayerField.build).parameters)
```

  and in `tests/test_rarefied_follower.py` (uses the session `twin`; 90 km is inside its range):

```python
import math, numpy as np, pytest
from reentry_model import aero
from reentry_model.droplet_reentry import particle_trajectory as pt
from reentry_model.wake_residence.adapters import (WakeAtmosphere, rarefied_sphere_simulator, dense_solve,
    parent_collision_event, shock_layer_exit, release_state)
from reentry_model.wake_residence.overlap import Plume
from reentry_model.wake_residence.laws import LawSettings
from reentry_model.wake_residence.shock_layer import ShockLayerSettings, ShockLayerField, FreeStream
from reentry_model.wake_residence.affix import release_offset

def _setup(twin, h_km=90.0, theta=180.0, dv=0.0, settings=ShockLayerSettings(), d=1e-4, clearance=0.0):
    # build: t_r = twin.t_at_altitude(h_km); FreeStream from the twin's flow at t_r; field; WakeAtmosphere with the
    # field and mode "off"; initial state at P(t_r) + release_offset(theta, D1, v_hat, up) * (1 + clearance / R)
    # with the parent's velocity plus dv * v_hat; returns (sim, wake_atmosphere, field, t_r)
    ...

def test_drag_is_the_composed_law_and_finite_at_zero_relative_speed(twin):
    sim, wa, field, t_r = _setup(twin, settings=ShockLayerSettings(pressure_gradient=False))
    y0 = sim.initial_state_vector()
    a = sim.aero_state(0.0, y0[:3], y0[3:])
    assert np.all(np.isfinite(a.a_drag))
    fs = a.freestream
    cd = pt.drag_coefficient(a.V, fs.rho, fs.T, fs.p, 1e-4, temperature_ratio=300.0 / fs.T)
    mass = 2813.0 * math.pi * 1e-12 / 6
    assert np.allclose(a.a_drag, -0.5 * fs.rho * a.V * a.v_rel * cd * math.pi * 0.25e-8 / mass, rtol=1e-9)

def test_stagnation_release_moves_forward_relative_to_parent(twin):
    sim, wa, field, t_r = _setup(twin)
    run = dense_solve(sim, 0.2, events=[parent_collision_event(wa, field.shape.R)])
    x0, _, _, _ = wa.parent_frame(0.0, run.sol.sol(0.0)[:3])
    x1, _, _, _ = wa.parent_frame(0.2, run.sol.sol(0.2)[:3])
    assert x1 > x0 + 1e-3                                # the parent decelerates more than the follower

def test_follower_thrown_at_the_parent_hits_it(twin):
    sim, wa, field, t_r = _setup(twin, dv=-50.0, clearance=1e-3)
    run = dense_solve(sim, 1.0, events=[parent_collision_event(wa, field.shape.R)])
    assert shock_layer_exit(run, wa, field, hit=run.end_reason == "event").kind == "hit_parent"

def test_grid_halving_changes_exit_time_by_under_one_percent(twin):
    fine = ShockLayerSettings(n_phi=181, n_eta=41)
    t = []
    for s in (ShockLayerSettings(), fine):
        sim, wa, field, t_r = _setup(twin, theta=150.0, settings=s)
        run = dense_solve(sim, 20.0, events=[parent_collision_event(wa, field.shape.R)])
        t.append(shock_layer_exit(run, wa, field).t)
    assert t[1] == pytest.approx(t[0], rel=0.01)
```

  (The exact `dense_solve` end-reason string for an event follows Task 5's `DenseRun.end_reason`; use it.)

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_shock_layer.py reentry_model/wake_residence/tests/test_rarefied_follower.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement `shock_layer.py`**, then the `adapters.py` additions. Compute the field on the (phi, eta) grid
  with phi from 0 to 90 degrees; at phi = 0 use the axis limit. For y > 0 the frame's e_y is the lateral unit vector;
  on the axis w_y = 0 by symmetry. Report in the commit message and the task report: epsilon, Lighthill's and
  Billig's standoffs at 80 km, the build time of one field, and the wall time of the 20 s follower run.
- [ ] **Step 4: Run them to confirm they pass**, then `"$PY" -m pytest reentry_model/wake_residence -q` (the Task 6
  WakeAtmosphere tests must pass unchanged).
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/shock_layer.py reentry_model/wake_residence/adapters.py \
  reentry_model/wake_residence/tests/test_shock_layer.py reentry_model/wake_residence/tests/test_rarefied_follower.py
git commit -m "wake_residence: the parent's shock-layer flow (Billig shape with a Lighthill standoff, equilibrium oblique-shock and Newtonian surface states, Lighthill nose velocities) fed to rarefied followers through WakeAtmosphere, with the composed drag law, a pressure-gradient force and the parent as an obstacle"
```

---

### Task 8: Scenarios, run names and case output

**Files:**
- Create: `reentry_model/wake_residence/scenarios.py`, `reentry_model/wake_residence/tests/test_scenarios.py`
- Modify: `reentry_model/wake_residence/adapters.py` (add `state_from_ecef`)

**Interfaces:**
- Consumes: everything above (Tasks 1-7b, and the 6b front-plane change: `ResidenceResult.reason` values include
  `reached_parent` and `ahead_of_parent`; `ExitState.through`); `droplet_reentry.model.DropletSimulator`, `weber`,
  `WE_CRIT`, `T_RELEASE`, `VOLATILE`; `wake_column.run._git_commit`, `_git_dirty`; `affix.*`; Task 7b's
  `ShockLayerSettings`, `FreeStream`, `ShockLayerField`, `WakeAtmosphere(..., shock_layer=)`,
  `rarefied_sphere_simulator`, `parent_collision_event`, `shock_layer_exit`.
- Produces:
  - In `adapters.py`: `state_from_ecef(r, v, epoch) -> InitialState` (geodetic position and flight angles from an
    Earth-fixed position and velocity; `release_state` should call it).
  - `@dataclass(frozen=True) ScenarioSettings(law: LawSettings = LawSettings(), bound="low", wake_flow="off",
    density=2813.0, frontal_area=6.215, beta_scale=1.0, k_E=1.0, t_max=600.0, scan_dt=None,
    affix: AffixSettings = AffixSettings(), kn_free_molecular=10.0,
    shock: ShockLayerSettings = ShockLayerSettings())`.
  - `follower_regime(kn_f, kn_continuum, kn_free_molecular) -> str` — `continuum` / `transitional` /
    `free_molecular` (§4.3c), and `AFFIX_STATUS = {"continuum": "governing", "transitional": "indicative",
    "free_molecular": "not_applicable"}`.
  - `run_name(kind, d_mm, h_km, settings, droplet_d_um=None, release="far_wake", theta_deg=None) -> str` — e.g.
    `sphere_d0050.0mm_h080.0km_farwake_low_both_taylor`, or `..._escape_th160.0_...` for escape releases, or
    `..._surface_th160.0_...` for rarefied followers (the far-wake run name has no theta: one run serves every
    far-wake theta), or `droplet_s0050.0mm_d0100.0um_h080.0km_...`; plus `_<field><value>` for each
    `LawSettings.non_default()` entry and each non-default `ScenarioSettings` field (including non-default
    `AffixSettings` and `ShockLayerSettings` fields), in sorted order.
  - `@dataclass CaseResult(name, config: dict, residence: ResidenceResult | None, summary: dict, history: list[dict])`.
  - `classify_release(d_mm, h_km, theta_deg, twin, settings) -> AffixResult` — the parent's EmissionState at the
    twin's release time, the sphere density, the parent mass at that altitude (`twin.beta.mass_at`), D1 = equal-area
    diameter.
  - `sphere_in_spacecraft_plume(d_mm, h_km, theta_deg=90.0, settings=ScenarioSettings(), twin=None) -> CaseResult`.
    Raises ValueError unless 90 <= theta_deg <= 180 (windward only, Asha). Order (§4.3c):
    1. Parent regime from `laws.classify` at t_r: `rarefied` or `too_slow` → no integration, `residence=None`, reason
       `rarefied_no_law` / `too_slow_no_law`.
    2. Build the shock-layer field from the parent's free stream at t_r (`FreeStream` from the twin's flow arrays;
       R = D1/2, r0). Kn_f = `aero.mean_free_path(rho, m_bar)` of `field.surface_state(theta)` / d; follower regime
       from `follower_regime(Kn_f, settings.affix.kn_continuum, settings.kn_free_molecular)`. Classify with AFFIx
       in every case and record `affix_status = AFFIX_STATUS[regime]`.
    3. Rarefied follower (transitional or free_molecular): release kind `surface`, position P(t_r) +
       `release_offset(theta, D1, v_hat, up)`, velocity = parent velocity; integrate
       `rarefied_sphere_simulator(..., wall_T_K=settings.shock.follower_wall_T_K,
       pressure_gradient=settings.shock.pressure_gradient)` through `WakeAtmosphere(..., mode=settings.wake_flow,
       shock_layer=field)` with `events=[parent_collision_event(wa, D1/2)]`; then `shock_layer_exit`, and
       `residence(...)` as for the other kinds (reason `hit_parent` when the collision event ended the run; AFFIx's
       zone does not change what is integrated). Transitional: flag `transitional_follower`.
    4. Continuum follower: the AFFIx branches below, exactly as before (no shock-layer field in its WakeAtmosphere):
    - `far_wake`: release position P(t_r) + `far_wake_offset(D1, v_hat)`, velocity = parent velocity at t_r;
    - `escape`: release position P(t_r) + `release_offset(theta, D1, v_hat, up)` with `up` the local geodetic up at
      P(t_r), velocity = parent velocity;
    - `shock_band`: no integration; `residence=None`, summary reason `"shock_band"`, residence_time NaN;
    - `trapped`: no integration; reason `"near_wake_trapped"`, `censored=True`, residence_time NaN.
    Integrated cases: `dense_solve` of `sphere_simulator(state_from_ecef(...), d, WakeAtmosphere(..., t_offset=t_r,
    mode=settings.wake_flow))` with t_max = min(settings.t_max, wake_atmosphere.follower_t_max) (Task 6 finding: a
    longer run passes the twin's end); `solution_track`; `residence(Plume(twin.track, settings.law, settings.bound),
    follower, t_start=t_r, t_end=follower.t1)`.
    `summary` has: diameter, altitude, theta, release kind (far_wake/escape), the AFFIx zone, theta_inn, theta_bif,
    t_hat, t_fz and every AFFIx flag, `residence_exceeds_t_fz` (escape: residence > 3 t_fz), bound, wake_flow,
    residence_time, reason, censored, total_inside, n_intervals, `follower_t_end` and `follower_end_reason` (the
    integration's end time on the parent clock and its end reason — Task 3 ruling: a censored case reports when it
    ended), exit altitude/speed/lag/distance/r_b/s/through, regimes and laws at release and exit, Kn/Re/P at
    release, every plume flag, interp_error_m, sphere and parent ballistic coefficients at release, parent sha256,
    git_commit, git_dirty. Added by §4.3c-d: `kn_f_release` (shock-layer surface gas), `kn_f_free_stream` (the same
    follower in the parent's free stream, for comparison), `kn_f_exit` (the WakeAtmosphere gas at exit or end),
    `follower_regime`, `affix_status`, `release_kind`, `t_first_entry`, `first_interval`, `total_inside`; and for
    surface releases `shock_layer_exit`, `t_shock_layer_exit`, `x_exit_parent_frame`, `y_exit_parent_frame`,
    `v_rel_shock_layer_exit`, `hit_theta_deg` (NaN unless hit), `epsilon`, `standoff_m`, `standoff_billig_m`, `mach`,
    `a_p_max`, `a_drag_max_shock_layer`, `q_change`, `cd_composed_release`, `cd_sesam_release`,
    `field_build_time_s`, and the flags `mach_outside_billig`, `shock_layer_nonequilibrium`, `transitional_follower`,
    `q_change_large`.
  - Summary keys the tests read: `reason`, `residence_time`, `parent_sha256`, `git_commit`, `git_dirty`,
    `regime_at_release`, `regime_at_exit`, `kn_at_release`, `re_at_release`, `p_at_release`, `affix_zone`, `t_fz`,
    `follower_t_end`, `follower_regime`, `affix_status`, `release_kind`, `kn_f_release`, `kn_f_free_stream`,
    `shock_layer_exit`. History keys: `t`, `follower_x/y/z`, `parent_x/y/z` (parent position at the same t), `lag_m`,
    `distance_m`, `r_b`, `margin`, `regime`.
  - `droplet_in_sphere_plume(sphere_d_mm, h_km, droplet_d_um, settings=ScenarioSettings(), twin=None) -> CaseResult`
    — sphere released from the twin at h_km at the far-wake start, integrated with `wake_flow="off"`; parent Track
    from `sphere_parent_track`; droplet `InitialState` = the sphere's state at its release (no AFFIx offset for the
    droplet: the droplet-from-sphere release geometry is not decided; record `release="centroid"`); Weber pre-check
    as `run_case` does (`weber(fs.rho, V, d0, T_RELEASE) >= WE_CRIT` → summary `reason="skipped_weber"`, no run);
    `DropletSimulator(release, droplet_d_um, atmosphere=WakeAtmosphere(...))` solved with
    `dense_solve(..., method="LSODA", rtol=1e-8, atol=<the list DropletSimulator.run builds>, events=[ground, gone])`
    replicating `DropletSimulator.run`'s two events (`model.py` lines 463-469).
  - `write_case(case, out_dir) -> (csv_path, json_path)` into `<out_dir>/runs/` (no CSV when residence is None).

- [ ] **Step 1: Write the failing tests** in `tests/test_scenarios.py`:

```python
import json, math, numpy as np, pytest
from reentry_model.wake_residence.scenarios import (ScenarioSettings, run_name, sphere_in_spacecraft_plume,
    droplet_in_sphere_plume, write_case, classify_release)
from reentry_model.wake_residence.laws import LawSettings

def test_run_names_differ_for_every_setting():
    base = run_name("sphere", d_mm=50, h_km=80, settings=ScenarioSettings())
    others = [ScenarioSettings(bound="high"), ScenarioSettings(wake_flow="both"),
              ScenarioSettings(law=LawSettings(taylor_K=0.5)), ScenarioSettings(frontal_area=3.1),
              ScenarioSettings(beta_scale=1.1), ScenarioSettings(density=2700.0)]
    names = {run_name("sphere", d_mm=50, h_km=80, settings=s) for s in others}
    assert base not in names and len(names) == len(others)
    assert run_name("sphere", 50, 80, ScenarioSettings(), release="escape", theta_deg=160.0) != base

def test_far_wake_case_runs_and_writes(twin, tmp_path):
    c = sphere_in_spacecraft_plume(10, 60, 90.0, ScenarioSettings(wake_flow="both"), twin=twin)
    csv_path, json_path = write_case(c, tmp_path)
    s = json.load(open(json_path))
    assert s["affix_zone"] == "far_wake" and s["t_fz"] > 0
    assert s["reason"] in ("ok", "censored", "reached_parent") and s["residence_time"] >= 0
    assert s["parent_sha256"] and "git_commit" in s and s["regime_at_release"] and s["follower_t_end"] > 0
    h = c.history                                      # starts on the axis, two parent diameters behind
    gap = np.linalg.norm([h[0]["follower_x"] - h[0]["parent_x"], h[0]["follower_y"] - h[0]["parent_y"],
                          h[0]["follower_z"] - h[0]["parent_z"]])
    assert gap == pytest.approx(2 * 2 * math.sqrt(6.215 / math.pi), rel=1e-6)

def test_follower_regime_thresholds():
    from reentry_model.wake_residence.scenarios import follower_regime, AFFIX_STATUS
    assert [follower_regime(k, 0.01, 10.0) for k in (0.005, 0.5, 10.0)] == ["continuum", "transitional",
                                                                             "free_molecular"]
    assert [AFFIX_STATUS[r] for r in ("continuum", "transitional", "free_molecular")] == [
        "governing", "indicative", "not_applicable"]

def test_rarefied_follower_starts_on_the_surface_and_reports_the_shock_layer(twin):
    c = sphere_in_spacecraft_plume(1, 90, 150.0, twin=twin)
    s = c.summary
    assert s["follower_regime"] in ("transitional", "free_molecular") and s["release_kind"] == "surface"
    assert s["affix_status"] in ("indicative", "not_applicable")
    assert s["kn_f_release"] < s["kn_f_free_stream"]                     # the shock-layer gas is denser
    assert s["shock_layer_exit"] in ("through_shock", "past_base_plane", "hit_parent", "inside_at_end")
    h = c.history
    gap = np.linalg.norm([h[0]["follower_x"] - h[0]["parent_x"], h[0]["follower_y"] - h[0]["parent_y"],
                          h[0]["follower_z"] - h[0]["parent_z"]])
    assert gap == pytest.approx(math.sqrt(6.215 / math.pi), rel=1e-6)   # on the surface, D1/2 from the centroid

def test_escape_and_shock_band_and_windward_limit(twin):
    esc = sphere_in_spacecraft_plume(10, 60, 180.0, twin=twin)
    assert esc.summary["follower_regime"] == "continuum" and esc.summary["affix_status"] == "governing"
    assert esc.summary["affix_zone"] == "escape" and esc.summary["residence_time"] < 0.05
    tb = classify_release(10, 60, 90.0, twin, ScenarioSettings()).theta_bif
    band = sphere_in_spacecraft_plume(10, 60, tb, twin=twin)
    assert band.summary["reason"] == "shock_band" and band.residence is None
    with pytest.raises(ValueError):
        sphere_in_spacecraft_plume(10, 60, 45.0, twin=twin)

def test_droplet_in_sphere_plume_worked_examples(twin):
    # controller ruling R2: (50, 80, 100) may be rarefied (no law) -- assert completion and the recorded regime;
    # (100, 60, 50) is a continuum, laminar sphere parent and must give a finite residence with validity records
    c = droplet_in_sphere_plume(50, 80, 100, twin=twin)
    assert c.summary["regime_at_release"] and c.summary["kn_at_release"] > 0
    c = droplet_in_sphere_plume(100, 60, 50, twin=twin)
    assert c.summary["reason"] != "skipped_weber"
    assert math.isfinite(c.residence.residence_time) and c.residence.residence_time >= 0 and c.residence.validity
    h = c.history
    assert np.allclose([h[0]["follower_x"], h[0]["follower_y"], h[0]["follower_z"]],
                       [h[0]["parent_x"], h[0]["parent_y"], h[0]["parent_z"]], atol=1e-3)
```

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_scenarios.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement `scenarios.py` and `state_from_ecef`.** If the (100, 60, 50) droplet is skipped by the Weber
  check, use a smaller droplet that is not, and note it in the commit message.
- [ ] **Step 4: Run them to confirm they pass**, then `"$PY" -m pytest reentry_model/wake_residence -q`.
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/scenarios.py reentry_model/wake_residence/adapters.py \
  reentry_model/wake_residence/tests/test_scenarios.py
git commit -m "wake_residence: scenarios -- sphere in the spacecraft's plume, released by the follower's own regime (continuum: the AFFIx outcome, far wake on the axis two parent diameters behind or escape at the windward surface point; rarefied: the surface point, integrated through the shock layer) and a molten droplet in a solid sphere's plume, with configuration-encoding run names and case output"
```

---

### Task 9: Sweep, sensitivity and tables

**Files:**
- Create: `reentry_model/wake_residence/analysis/__init__.py`, `analysis/sweep.py`, `analysis/sensitivity.py`,
  `analysis/tables.py`, `reentry_model/wake_residence/tests/test_analysis.py`

**Interfaces:**
- Consumes: `scenarios.*`, `twin_parent`, `twin_vs_drama`, `drama_beta_table`, `plume_radius`, `classify`,
  `affix.*`, `Plume.wake_state`.
- Produces:
  - `sweep.main(argv=None)` — options `--diameters` (default `1 5 10 25 50 100 200 400` mm), `--altitudes`
    (default `100 90 80 70 60 50` km; any outside the twin's range is dropped with a printed message), `--thetas`
    (default `90 120 135 145 160 180` degrees), `--bounds` (default `low high`), `--wake-flow` (default `off both`),
    `--workers` (default CPU count less one), `--out` (default `wake_residence.OUTPUT`). For each (diameter,
    altitude, theta) decide the follower regime and classify; for continuum followers integrate the far-wake run
    once per (diameter, altitude, bound, wake_flow) and reuse it for every far-wake theta, integrate escape runs per
    theta, and do not integrate shock-band cases; integrate rarefied followers per theta from the surface (§4.3d);
    print the count of continuum, transitional and free-molecular followers per altitude. Resumable: a
    run whose JSON exists is skipped. A case that raises is reported at the end without stopping the rest. Writes
    `<out>/summary.csv` — one row per (diameter, altitude, theta, bound, wake_flow), sorted in that order, joining
    the classification with the run that serves it — and `<out>/crossover.csv` (per altitude and wake_flow, from the
    far-wake runs: the largest diameter with positive lag behind the parent at 1 s after release and the smallest
    with negative lag — the crossover bracket).
  - `sensitivity.main(argv=None)` — the §4.6 set at theta = 90 (far wake) (spheres 10, 50, 200 mm at 80 and 60 km;
    each variation one at a time against the baseline), plus the rarefied set of §4.6 (10 mm at 90 km, theta 180 and
    120: `shock_standoff="billig"`, `pressure_gradient=False`, `follower_wall_T_K=1000`, `shock_thickness_mfp` 1 and
    10; if 10 mm is continuum there, use the largest default diameter that is rarefied and say so),
    `<out>/sensitivity.csv` with columns case, parameter, value, bound, residence_time, change_vs_baseline,
    shock_layer_exit, t_shock_layer_exit.
  - `tables.main(argv=None)` — writes:
    - `<out>/plume_size.csv` (for the twin's state at each default release altitude: regime, Kn, Re_d, P, and r_b
      low/high at s = 10, 100, 1000, 3000, 10000, 30000 m with time since passage s/U);
    - `<out>/ballistic_coefficient.csv` (twin beta and C_D A against altitude, and each default sphere's ballistic
      coefficient there);
    - `<out>/twin_vs_drama.csv` — `twin_vs_drama` rows plus `twin_downrange_km` at the DRAMA time (great-circle
      distance from the twin's start point, `earth.great_circle_distance`) and its difference from DRAMA's downrange
      (Task 5 ruling: spec §4.3 asks for altitude AND downrange against time); a header comment notes the twin's
      ~6.5 s time lag behind DRAMA and that the engine uses the twin's own clock;
    - `<out>/affix_regimes.csv` (diameter, altitude, theta, zone, theta_inn, theta_bif, t_hat, t_fz, every flag);
    - `<out>/drag_ratio_check.csv` (per release altitude and bound: our drag ratio at the far-wake start point,
      rho_w (V - u)^2 / (rho_inf V^2) from `Plume.wake_state` with mode "both", next to Register's 0.14);
    - `<out>/wake_gas_sanity.csv` (spec §4.3a sanity checks: the model's core gas speed at 30 km behind the parent
      and core temperature 300 m behind, at each release altitude, next to the assumed values ~10 m/s and
      1,000-2,000 K).
- Commands (from `Reentry Model/`): `"$PY" -m reentry_model.wake_residence.analysis.sweep`,
  `... .analysis.sensitivity`, `... .analysis.tables`.

- [ ] **Step 1: Write the failing tests** in `tests/test_analysis.py`:

```python
import csv, filecmp, os, shutil
from reentry_model.wake_residence.analysis import sweep, tables

ARGS = ["--diameters", "10", "--altitudes", "60", "--thetas", "90", "180", "--bounds", "low",
        "--wake-flow", "off", "--workers", "1"]

def test_short_sweep_is_reproducible_and_resumable(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    sweep.main(ARGS + ["--out", str(a)]); sweep.main(ARGS + ["--out", str(b)])
    assert filecmp.cmp(a / "summary.csv", b / "summary.csv", shallow=False)
    rows = list(csv.DictReader(open(a / "summary.csv")))
    assert {r["affix_zone"] for r in rows} == {"far_wake", "escape"} and len(rows) == 2
    mtime = os.path.getmtime(next((a / "runs").glob("*.json")))
    sweep.main(ARGS + ["--out", str(a)])                                   # resumes: nothing rerun
    assert os.path.getmtime(next((a / "runs").glob("*.json"))) == mtime

def test_tables_report_regimes_as_outputs(tmp_path):
    tables.main(["--out", str(tmp_path)])
    rows = list(csv.DictReader(open(tmp_path / "plume_size.csv")))
    assert {r["regime"] for r in rows} <= {"rarefied", "laminar", "transitional", "turbulent", "too_slow"}
    for name in ("twin_vs_drama.csv", "ballistic_coefficient.csv", "affix_regimes.csv", "drag_ratio_check.csv",
                 "wake_gas_sanity.csv"):
        assert (tmp_path / name).exists(), name
    assert "twin_downrange_km" in open(tmp_path / "twin_vs_drama.csv").read()
```

(The `summary.csv` comparison requires that no wall-clock time is written into it; keep `runtime_s` in the
per-case JSON only.)

- [ ] **Step 2: Run them to confirm they fail.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_analysis.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement** the three scripts.
- [ ] **Step 4: Run them to confirm they pass.** Same command.
- [ ] **Step 5: Run the real thing** (in the background; report wall time and any failed case):
  `"$PY" -m reentry_model.wake_residence.analysis.tables`, then `... .sweep`, then `... .sensitivity`.
  Expected: `summary.csv` with 8 × 6 × 6 × 2 × 2 = 1152 rows less any dropped altitudes; report to the controller
  the twin-vs-DRAMA figures, the regime at each release altitude, the AFFIx zone counts, the follower-regime counts,
  the shock-layer outcome counts for rarefied followers, the drag-ratio check, the wake-gas sanity values and the
  crossover brackets.
- [ ] **Step 6: Commit** (code and tests only; outputs are git-ignored)

```bash
git add reentry_model/wake_residence/analysis/__init__.py reentry_model/wake_residence/analysis/sweep.py \
  reentry_model/wake_residence/analysis/sensitivity.py reentry_model/wake_residence/analysis/tables.py \
  reentry_model/wake_residence/tests/test_analysis.py
git commit -m "wake_residence: resumable sweep over diameter, altitude and release angle (far-wake runs shared across angles), the sensitivity set, and the plume-size, ballistic-coefficient, twin-versus-DRAMA, AFFIx, drag-ratio and wake-gas tables"
```

---

### Task 10: Plots

**Files:**
- Create: `reentry_model/wake_residence/analysis/plots.py`; add a test to `tests/test_analysis.py`

**Interfaces:**
- Consumes: the CSV and run files written by Task 9; `DRAMA Sphere Experiments/analysis/plot_style.py` the way
  `reentry_model/compare.py` imports it (fall back to plain matplotlib without it).
- Produces: `plots.main(argv=None)` with `--out` (default `wake_residence.OUTPUT`), writing to `<out>/plots/`:
  `plume_size.png` (r_b low/high against s and against time since passage, one line per release altitude,
  Taylor's 1200 (C_D A)^(1/2) and Park's 30 km limits marked, regime in the legend), `trajectories_<h>km.png`
  (twin and far-wake spheres, altitude against downrange, for each release altitude present), `margin_<h>km.png`
  (distance to track and r_b against time for every diameter's far-wake run at that altitude, wake flow off and
  on), `residence_map.png` (far-wake residence time over diameter × altitude, one panel per wake_flow, low–high
  shown as two panels or as hatching where they differ, censored, reached-parent, ahead-of-parent and no-law cases
  marked), `affix_zones.png` (zone over diameter × release angle with theta_bif drawn, one panel per altitude only
  if the validity flags differ by altitude, otherwise one panel), and `shock_layer_paths_<h>km.png` for each
  altitude with rarefied followers (their paths in the parent frame, x forward against y from the axis, over the
  first second or to the shock-layer exit if later, with the body circle and the shock curve drawn and the exit kind
  in the legend; skipped with a printed note where there are none). Prints how many runs it found.

- [ ] **Step 1: Write the failing test** (append to `tests/test_analysis.py`):

```python
from reentry_model.wake_residence.analysis import plots

def test_plots_are_written_from_a_short_sweep(tmp_path):
    sweep.main(ARGS + ["--out", str(tmp_path)]); tables.main(["--out", str(tmp_path)])
    plots.main(["--out", str(tmp_path)])
    for name in ("plume_size.png", "residence_map.png", "margin_60km.png", "trajectories_60km.png",
                 "affix_zones.png"):          # ARGS has no rarefied follower, so no shock_layer_paths plot here
        assert (tmp_path / "plots" / name).stat().st_size > 0
```

- [ ] **Step 2: Run it to confirm it fails.** Run: `"$PY" -m pytest reentry_model/wake_residence/tests/test_analysis.py -q -k plots`. Expected: ImportError.
- [ ] **Step 3: Implement `plots.py`.**
- [ ] **Step 4: Run it to confirm it passes**, then run `"$PY" -m reentry_model.wake_residence.analysis.plots`
  on the real outputs and look at each image.
- [ ] **Step 5: Commit**

```bash
git add reentry_model/wake_residence/analysis/plots.py reentry_model/wake_residence/tests/test_analysis.py
git commit -m "wake_residence: plots of plume size, trajectories, distance-to-track against plume radius, the residence-time map and the AFFIx release zones"
```

---

### Task 11: Documentation and final checks

**Files:**
- Modify: `docs/model_assumptions.md` (new section 10, "Wake residence"), root `README.md` (a "Wake residence"
  section with the four commands), `CLAUDE.md` (a "The wake residence" subsection under the droplet and wake column
  ones, and the `pytest.ini` line). In the worktree CLAUDE.md has none of Asha's oxidation edits, so it IS committed
  on the branch (controller ruling R3; merging will need her uncommitted CLAUDE.md edits committed or set aside).

**Interfaces:**
- Consumes: the measured figures from Tasks 5, 6, 6b and 9.

- [ ] **Step 1: Write section 10** of `docs/model_assumptions.md` in the style of section 9: every assumption
  marked **(verified)** with the test that asserts it, or **(assumed)**. Include at least: the three laws with their
  sources, fitted ranges and caveats (Park's 0.7 second-hand from Schapker and Camac; Feldman's combined form
  assembled, not printed; R density-weighted); the regime criteria and thresholds (outputs, not altitude bands);
  the twin and its measured agreement with DRAMA (0.52 % in speed above 60 km; the ~6.5 s time lag and its gravity-
  model sensitivity; downrange agreement from Task 9); the frontal area; the union-of-balls definition and the
  front plane with the `reached_parent` / `ahead_of_parent` labels; the wake gas model, its well-mixed temperature,
  the density relative to the local ambient, and the drag-ratio check against Register's 0.14; the AFFIx release
  classifier, the windward-only release angles, the far-wake start two parent diameters behind (Asha's estimate from
  Register et al.) and the extrapolation below r = 1/16; the follower-regime gate (Kn from the shock-layer surface
  gas); the shock-layer field of §4.3d item by item (Billig shape, Lighthill standoff and nose velocities, equilibrium
  oblique jump, Newtonian surface state, linear interpolation across the layer, frozen free stream, smoothing, the
  composed drag law with the 300 K wall, the pressure-gradient force, the parent as an obstacle) with the measured
  epsilon, standoffs and field build time from Task 7b; the non-goals (no melting of the sphere, no integrated
  shock forces on continuum followers, no plume drift, no shock-layer path for droplets); the measured cost of a
  wake run.
- [ ] **Step 2: Add the README and CLAUDE.md sections.**
- [ ] **Step 3: Run the whole unit tier.** `"$PY" -m pytest -m "not reference" -q`. Expected: all pass (report the
  count and wall time).
- [ ] **Step 4: Confirm the protected files are untouched.**
  Run: `git diff --stat main -- reentry_model/trajectory.py reentry_model/droplet_reentry/model.py`.
  Expected: no change to either file on the branch.
- [ ] **Step 5: Commit**

```bash
git add docs/model_assumptions.md ../README.md CLAUDE.md
git commit -m "docs: wake residence -- assumptions with their verification status, the plume laws' and AFFIx sources and ranges, and the commands"
```
