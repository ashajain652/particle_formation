# Sphere re-entry trajectory model (Step 1 of the physics model) — Design

Date: 2026-09-17
Status: design for review, awaiting implementation plan

## 1. Purpose

First step of a first-principles model of an AA7075 sphere re-entering the atmosphere
(later steps: 3D finite-element heating of the body, then melt spraying). This step is
the **trajectory only**: a solid sphere of constant mass, from an initial velocity,
altitude and flight-path angle, through the free-molecular, transitional and continuum
regimes to the ground, with the physics formulated the way ESA DRAMA's SARA/SESAM does
it, and **verified against SESAM runs of the same sphere**: overlaid plots of velocity
and altitude against time, with error metrics.

The model must be usable on its own (a CLI and an importable package) and must leave
clean attachment points for Step 2 (heating and body temperature) without changing the
trajectory code.

Planning references, the measured SESAM facts this design relies on, and the SESAM
reference runs are described in
`Literature Review/Sphere Demise Model - Planning References/` (`references.md`,
`sesam_facts/sesam_verified_facts.md`) and in this repository's README ("Reference runs").

## 2. Verified facts about SESAM (measured 2026-09-17 on DRAMA 4.1.4)

Implementers rely on these; every item was measured on this machine from
`TOOLS/SARA/REENTRY/data/ATDB_SPHERE.nc` and the run histories under
`sphere_sweep_output/`.

| Topic | Fact |
|---|---|
| Equations of motion | 3-DOF point mass over a rotating Earth. SESAM's dγ/dt at the 100 mm reference start is −0.00790 °/s; Vinh's rotating-Earth equations give −0.00792 °/s, the non-rotating form −0.00650 °/s. Velocity, flight-path angle and heading are relative to the rotating atmosphere. Gravity includes J2, J4, J22 (DRAMA documentation via the Jeanne thesis); Earth shape WGS-84; altitude is geodetic. |
| Atmosphere | Selectable per run through the wrapper (`--atmosphere`). `static` = `StaticEnvironmentData.csv` = the US Standard Atmosphere 1976 (0–150 km every 100 m: p, ρ, T, wind N/E/down, γ, O-fraction), no winds. `nrlmsise` = NRLMSISE-00 with daily F10.7/Ap from `data/fap_day.dat` (`fap_mon.dat`) and HWM14 winds (log lines "Using dynamic atmosphere model NRLMSISI-00 with dynamic solar activity", "Using dynamic wind model HWM14" / "Simulating without wind"). |
| Knudsen number | Column `knudsen` = λ/D with a hard-sphere mean free path: λ = 1/(√2 π d² n), n = ρ/m̄. With d = 3.65e-10 m and m̄ = 28.96 u it reproduces Kn at 77.5 km (2.98 mm); at 115 km (NRLMSISE-00) the reported Kn = 41.0 corresponds to m̄ ≈ 26.7 u, i.e. SESAM uses the local mean molecular mass. |
| Mach number | Ma = V/√(γ R T/m̄) with γ = 1.4 (26.2 at 77.5 km, T = 203 K). |
| Drag, Ma ≥ 5 | `ATDB_SPHERE.nc`: C_D,fm(Ma) = 2.360635, 2.147897, 2.089229, 2.062249, 2.046842, 2.036933 and C_D,c(Ma) = 0.898818, 0.910198, 0.912322, 0.913067, 0.913411, 0.913599 at Ma = 5, 10, 15, 20, 25, 30; linear in Ma; clamped above 30. |
| Drag, Ma < 5 | Continuum: C_D,c = 0.898818 (the Ma-5 value) for 1 ≤ Ma < 5 and exactly half of it, 0.449409, for Ma < 1 (466 samples, no scatter). |
| Drag bridging | C_D = C_D,c + (C_D,fm − C_D,c) f(Kn) with f measured on 764 samples: f = ½[1 + erf((log10 Kn + 0.845)/0.585)], rms 0.006 in f (f = 0.003, 0.05, 0.18, 0.50, 0.78, 0.93, 0.98 at Kn = 0.01, 0.03, 0.06, 0.14, 0.3, 0.6, 1). The textbook sin²[π(0.5 + 0.25 log10 Kn)] does not fit (rms 0.09). |
| Output cadence | Histories at ~1 s (adaptive); trajectory and aerothermal files joined on time by the wrapper. |
| Heating / thermal | Not part of Step 1. (DKR × 0.2747 shape factor, Klett free-molecular, same f(Kn), εσT⁴ to 0 K, lumped mass — recorded in the facts note for Step 2.) |
| pyDRAMA quirk | Booleans reach SESAM as `True`/`False`, which it does not parse; the wrapper now writes `yes`/`no` (commit a687828). |

## 3. Reference cases

Produced with the wrapper (commit a687828), material `AA7075_nomelt` (drama-AA7075 with
T_melt = 1e5 K, so no mass loss), 300 K, 7.5 km/s, γ = −0.959331°, heading 347.168296°,
lat 29.546067°, lon −82.134333°, epoch 2024-08-01T12:53:07, NRLMSISE-00, winds on and off.
Files: `sphere_sweep_output/reference_AA7075_nomelt/runs/<name>.csv|.json` (raw DRAMA
trees under `raw/`); the first implementation task copies the four CSV+JSON pairs into
`data/reference_runs/` and commits them, so the verification does not depend on
git-ignored output.

| Case | Run name (`..._nowind` for winds off) | Coverage |
|---|---|---|
| R100 | `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis` | Kn 0.03 → continuum; impact 368 s |
| R50 | `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis` | Kn 41 → 10 (106.5 km) → 1 (93.8 km) → 0.1 (80.4 km) → 0.01 (65 km); impact 561 s |

Measured sensitivities that bound what the comparison can resolve: NRLMSISE-00 vs US76
changes R100 by ≤ 57 m/s and ≤ 0.4 km (density ratio 0.96–1.04); HWM14 winds change either
case by ≤ 8 m/s and ≤ 20 m.

## 4. Scope

In: constant-mass sphere; rotating-Earth 3-DOF dynamics with J2 gravity (J4/J22 optional);
atmosphere = NRLMSISE-00 (pymsis) with fap-file solar activity, or the US76 table, or a
SESAM density/temperature replay; winds none or the static profile; SESAM's sphere drag
tables and Knudsen bridging (bridging swappable); DOP853 integration with ground/escape
events; CSV/JSON output in the wrapper's column layout; a comparison tool producing the
verification plots and metrics; unit tests; an integration test against the reference cases.

Out (later steps or never): heating and temperature (Step 2), mass loss and spraying
(Step 3), HWM14 winds (measured negligible; interface kept), lift, tumbling dynamics,
non-spherical shapes, fragmentation, Monte Carlo.

## 5. Architecture and layout

```
reentry_model/                      importable package (Python 3.12, drama_env; numpy, scipy, pymsis, matplotlib)
  __init__.py
  constants.py                      WGS-84, μ, J2/J4, ω_E, gas constants, hard-sphere d
  atmosphere.py                     Freestream state at (t, h, lat, lon): Atmosphere protocol + 3 implementations
  aero.py                           Knudsen, Mach, sphere C_D (tables + bridging functions)
  earth.py                          geodetic <-> ECEF, ENU frame, gravity with J2 (J4 optional)
  trajectory.py                     state, right-hand side, integrator, events, history assembly
  sesam_io.py                       read a wrapper run CSV/JSON as a Reference (times, columns, inputs)
  compare.py                        metrics + plots model vs reference
  cli.py                            `python -m reentry_model ...`
  data/atdb_sphere.json             the 6-Mach tables copied from ATDB_SPHERE.nc (provenance in the file)
  data/us76_static_environment.csv  copy of DRAMA's StaticEnvironmentData.csv (provenance in a header line)
data/reference_runs/                the four SESAM reference CSV+JSON pairs (committed by plan task 1)
tests/test_reentry_model_*.py       unit tests; `-m drama_reference` integration test
reentry_model_output/               git-ignored: histories, comparison plots and metrics
docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md   this file
```

Dependency direction: `cli` → `compare`, `trajectory`, `sesam_io`; `trajectory` → `aero`,
`atmosphere`, `earth`; nothing imports `cli` or `compare`. All quantities SI inside the
package; the CLI and CSV use the wrapper's units (km, km/s, deg, mm).

Step 2 attachment point: `trajectory` takes a `Body` object with `mass(t)`,
`temperature(t)` and `on_step(t, state, freestream, aero)`; the right-hand side reads
`body.mass(t)`, and the history writer reads `body.temperature(t)`. Step 1 ships
`ConstantBody` (constant mass, temperature = initial temperature, `on_step` a no-op).
Step 2 supplies a body whose `on_step` advances the thermal solver; `trajectory.py` does
not change.

## 6. Physics

### 6.1 State and frames
Integrate in Earth-centred Earth-fixed (ECEF) Cartesian coordinates: position **r** and
velocity **v** relative to the rotating Earth (= relative to a wind-free atmosphere).
This avoids the polar singularities of the (V, γ, ψ) form and makes J2 exact. Conversions:
- geodetic (h, φ, λ) ↔ ECEF on WGS-84 (a = 6378137 m, f = 1/298.257223563); altitude is
  geodetic, as in SESAM.
- Initial velocity from (V, γ, ψ): in the local East-North-Up frame at the initial point,
  **v** = V (cos γ sin ψ, cos γ cos ψ, sin γ) with ψ the heading clockwise from north (SESAM's
  heading 347.2° = NNW, consistent with the parent ground track).
- Outputs V = |**v**|, γ = asin(v_up/V), ψ = atan2(v_east, v_north), downrange = great-circle
  distance from the initial point on the WGS-84 sphere of radius a.

### 6.2 Equations of motion (rotating frame)
d**r**/dt = **v**
d**v**/dt = **g**(**r**) − 2 **ω**×**v** − **ω**×(**ω**×**r**) + **a**_drag,
**a**_drag = −½ ρ |**v**_rel| **v**_rel C_D A / m, **v**_rel = **v** − **w** (wind, zero by default),
A = π R², m = ρ_material 4/3 π R³ (constant in Step 1).
**g** from the potential U = −μ/r [1 − Σ J_n (a/r)ⁿ P_n(sin φ_c)] with J2 always and J4
optional (J22 omitted: its effect on a ~10 min flight is below the verification
resolution; a test will quantify J2 vs point-mass, expected ≤ 2 m/s in V). Constants:
μ = 3.986004418e14 m³/s², J2 = 1.08262668e-3, J4 = −1.61962e-6, ω = 7.2921159e-5 rad/s.

Check used as a unit test: at the R100 initial state the rate dγ/dt from this formulation
must equal −0.00792 °/s ± 2 % (Vinh's rotating-Earth value that SESAM reproduces).

### 6.3 Atmosphere (`Atmosphere` protocol: `state(t, h, lat, lon) -> Freestream(rho, T, p, m_bar, wind_enu)`)
- `NRLMSISE00Atmosphere` — pymsis (NRLMSISE-00 option, not MSIS 2.x) at the run epoch + t,
  geodetic (lat, lon, h). Solar/geomagnetic inputs from `data/fap_day.dat`: F10.7 of the
  previous day, the file's 81-day mean column (F3M) as F10.7a, daily Ap; species number
  densities give the mean molecular mass m̄ used by Kn and Ma. Cached per call step.
- `US76TableAtmosphere` — `StaticEnvironmentData.csv` (copied into `reentry_model/data/`)
  with log-linear interpolation of ρ and p and linear T; m̄ = 28.96 u; static winds available.
- `ReplayAtmosphere` — ρ(h) and T(h) taken from a reference run CSV (log-linear in altitude),
  m̄ back-solved from the reference Kn column. Diagnostic only: it isolates the dynamics from
  the atmosphere so that a mismatch can be attributed.
Winds: `NoWind` (default) or the static profile's N/E/down columns; HWM14 deferred.

### 6.4 Regime numbers (`aero.py`)
λ = 1/(√2 π d² n), n = ρ/m̄, d = 3.65e-10 m; Kn = λ/D; a = √(γ R T/m̄), Ma = V/a.
Requirement: reproduce SESAM's `knudsen` and `mach` columns of the reference runs to
within 2 % and 1 % respectively with the replay atmosphere (which supplies SESAM's own ρ, T).

### 6.5 Sphere drag (`aero.py`)
C_D,fm(Ma), C_D,c(Ma): linear interpolation of the six-point tables for 5 ≤ Ma ≤ 30, clamped
at 30; for Ma < 5: C_D,c = 0.898818 (1 ≤ Ma < 5) and 0.449409 (Ma < 1); C_D,fm clamped at its
Ma-5 value (irrelevant in practice: Kn is tiny wherever Ma < 5).
C_D = C_D,c + (C_D,fm − C_D,c) f(Kn), with `Bridging` objects:
- `SesamErf` (default): f = ½[1 + erf((log10 Kn + 0.845)/0.585)]
- `Sin2(kn_lo, kn_hi)`: f = sin²[(π/2)(log10 Kn − log10 kn_lo)/(log10 kn_hi − log10 kn_lo)]
  clipped to [0, 1] (the textbook form is kn_lo = 0.01, kn_hi = 1)
- `Matting`: placeholder raising NotImplementedError until the 1971 relation is transcribed
  from the paper (Step 2 uses it for heating; the class exists so the interface is fixed).
Requirement: with the replay atmosphere, reproduce SESAM's `drag` column to within 0.01
(absolute C_D) for Ma ≥ 5 and exactly below.

### 6.6 Integration
`scipy.integrate.solve_ivp`, method DOP853, rtol 1e-9, atol 1e-6 m / 1e-9 m/s, dense
output. Events: geodetic altitude ≤ 0 (ground impact, terminal), altitude ≥ 150 km
(escape, terminal), t ≥ t_max (default 3600 s). The history is sampled at a fixed cadence
(default 1 s) plus the terminal point; when a reference run is given, additionally at the
reference's own time stamps so the comparison needs no interpolation.

## 7. Command line

```
python -m reentry_model run --velocity 7.5 --altitude 115 --flight-path-angle -0.959331 \
    --diameter 50 [--heading 347.168296 --lat 29.546067 --lon -82.134333 --epoch 2024-08-01T12:53:07] \
    [--material-density 2813] [--temperature 300] \
    [--atmosphere nrlmsise|us76|replay:<reference.csv>] [--wind none|static] \
    [--bridging sesam-erf|sin2|textbook] [--gravity j2|point|j2j4] \
    [--reference <reference.csv>] [--outdir reentry_model_output] [--cadence 1.0] [--name <run_name>]
python -m reentry_model compare --model <model.csv> --reference <reference.csv> [--outdir ...]
```
`run` writes `<outdir>/<name>.csv` and `<name>.json`; with `--reference` it also runs
`compare` and writes the metrics into the JSON. Defaults reproduce the reference-case
settings (parent heading/lat/lon/epoch, NRLMSISE-00, no wind, SesamErf, J2). Exit codes:
0 ok, 1 integration failed/escaped, 2 bad arguments.

## 8. Outputs

**History CSV** — the wrapper's columns so the comparison and plotting code is shared:
`time_s, altitude_km, velocity_kms, temperature_K, mass_kg, thick_mm, lat_deg, lon_deg,
downrange_km, flight_path_deg, heading_deg, drag, lift, side, knudsen, mach, density_kgm3,
dynamic_pressure_Pa, load_factor_g` (temperature constant, mass constant, thick = radius in
mm, lift = side = 0; the aerothermal columns are absent in Step 1). Precision: 6
significant digits (more than SESAM's 3-decimal km/s so that residuals are ours, not
rounding's).

**Run JSON** — `inputs` (state, sphere, material density), `settings` (atmosphere, wind,
bridging, gravity, tolerances, cadence, fap-file values used: F10.7, F10.7a, Ap), `results`
(impact time, final velocity, max deceleration and its altitude, max dynamic pressure,
Kn at start, altitudes where Kn crosses 10/1/0.1/0.01), `comparison` (see below, when a
reference was given), `provenance` (git commit of this repo, package versions, reference
file SHA-256).

**Comparison** (`compare.py`) — model sampled at the reference times; metrics: max and RMS
of ΔV (m/s) and of ΔV/V_ref, max and RMS of Δh (m), Δ(time to impact) and Δ(final
velocity), all reported over the whole flight and over the hypersonic phase (V > 1 km/s).
Plots (PNG, `analysis/plot_style.py` palette): (1) V(t) model and SESAM overlaid with a
residual panel; (2) h(t) likewise; (3) h(V); (4) γ(t) and heading; (5) ground track; (6) Kn
and C_D vs time. Each plot titles the case and the settings.

## 9. Verification and acceptance

Sequence, all four reference cases (R100, R50 × winds on/off):
1. Replay atmosphere + SesamErf + J2: isolates dynamics and drag. Expected: |ΔV| ≤ 0.2 % of
   V (≈ 15 m/s at 7.5 km/s) and |Δh| ≤ 100 m over the hypersonic phase; impact time within
   1 %.
2. NRLMSISE-00 (pymsis) instead of replay: adds the atmosphere implementation. Expected:
   |ΔV| ≤ 1 %, |Δh| ≤ 0.5 km over the hypersonic phase (NRLMSISE-00 vs US76 moved R100 by
   57 m/s and 0.4 km, so two implementations of the same model should differ by less).
3. Winds-on references compared with the wind-free model: the residual must be ≤ the measured
   wind effect (≤ 8 m/s, ≤ 20 m) — this is a consistency check, not a target to close.
The thresholds are expectations; the plan's last task runs the comparisons, records the
measured numbers in the run JSONs and in the README, and either confirms the thresholds
or explains a justified change. The integration test asserts the confirmed thresholds.

## 10. Testing

Unit (no DRAMA needed; fixtures = the four reference CSV/JSON pairs in `data/reference_runs/`):
- `earth`: geodetic↔ECEF round trip at the reference initial points to 1e-6 m; ENU velocity
  construction returns the input (V, γ, ψ); gravity magnitude at h = 0, 77.5 km vs g0 (a/r)²
  plus the J2 term; dγ/dt at the R100 initial state = −0.00792 °/s ± 2 %.
- `atmosphere`: US76 table reproduces the CSV rows exactly and interpolates monotonically;
  replay reproduces the reference density/temperature columns; pymsis at the R100 and R50
  initial rows within 3 % of SESAM's `density_kgm3` (a looser, informative tolerance — the
  measured value is recorded); fap-file parsing returns F10.7 = 234, F10.7a = 194, Ap = 19
  for 2024-08-01.
- `aero`: Kn and Ma from replay states match the reference columns (2 %, 1 %); table
  interpolation hits the six ATDB points; bridging limits f(1e-3) < 0.005 and f(10) > 0.995;
  C_D column reproduced within 0.01 for Ma ≥ 5 and exactly for Ma < 5; SesamErf values at the
  seven tabulated Kn.
- `trajectory`: vacuum, point-mass, non-rotating: energy and angular momentum conserved to
  1e-9 over 600 s; vacuum rotating frame: Jacobi integral conserved; ground event fires within
  1 mm of h = 0; cadence and reference-time sampling produce the expected rows.
- `sesam_io` / `compare`: reading the reference JSON/CSV, metrics on synthetic histories with
  known offsets, plot files created.
Integration (`-m drama_reference`, runs in seconds): the four comparisons of §9 with the
confirmed thresholds.

## 11. Environment

`drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`, Python 3.12, numpy 2.5,
matplotlib) plus `pip install pymsis==0.13.0 scipy pytest` (a macOS arm64 wheel for pymsis
0.13.0 exists). No DRAMA install is needed to run the model or the tests; `sphere_reentry.py`
is only needed to regenerate reference runs.
