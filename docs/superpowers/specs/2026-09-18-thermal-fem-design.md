# Coupled trajectory + 3D FEM heat transfer of a sphere (Step 2 of the physics model) — Design

Date: 2026-09-18
Status: design for review, awaiting implementation plan
Builds on: `docs/superpowers/specs/2026-09-17-reentry-trajectory-model-design.md` (Step 1, merged at df8aa65)

## 1. Purpose

Solve the re-entry trajectory of a solid AA7075 sphere and the temperature field inside it
**together**: at every macro time step the trajectory advances, an aerothermal model turns the
freestream state into a convective heat flux on every surface patch (free-molecular, transitional
and continuum regimes), and a three-dimensional finite-element conduction solve advances the
temperature with thermal radiation from every surface patch. Deliverables:

1. the coupled run's history (the Step 1 trajectory columns plus heating and temperature columns),
2. the temperature field over time as a VTK time series,
3. an MP4/GIF animation of the surface temperature on the 3D sphere, with stills,
4. a verification against SESAM's heat input and lumped temperature for the same sphere.

Mass stays constant in this step (the material's latent-heat term is zero); melting, mass loss and
spraying are Step 3. The architecture is chosen so Step 3 adds physics without restructuring:
the coupled loop already carries mass and radius, the material model is enthalpy-capable, the mesh
is an input, and node coordinates are state.

Reference material: `Literature Review/Sphere Demise Model - Planning References/`
(`references.md`; `sesam_facts/sesam_verified_facts.md` §§1–14; `pdfs/matting-2012-...pdf`).

## 2. Facts this design relies on (measured or transcribed, 2026-09-17/18)

| Topic | Fact |
|---|---|
| SESAM heating | Total convective power Q = A_sphere · [(1 − f(Kn))·0.27471·q_DKR + f(Kn)·q_FM]; q_DKR = 1.1035e8 R^-1/2 (ρ/1.225)^1/2 (V/7925)^3.15 W/m² reproduces SESAM within 2 % at Kn ≤ 0.02; f(Kn) is the same measured blend as the drag (`aero.SesamTable`); SESAM's own q_FM is ~13× below the textbook ½ρV³ (Klett coefficient), which matters only where Kn > 0.1 and the absolute heating is negligible. |
| SESAM shape factor | 0.27471 = surface-average ÷ stagnation flux, constant in Mach in the ATDB; applied to the whole surface area of a lumped mass under the "random tumbling" attitude assumption. No hot-wall correction detectable. |
| Lees integral | Lees' laminar sphere distribution over the windward hemisphere, divided by the whole sphere area: 0.196 (M→∞), 0.200 (M=10), 0.209 (M=5). A fixed-attitude, windward-only model delivers ~28 % less total heat than SESAM at the same stagnation flux. |
| SESAM radiation | rad_cooling = −ε σ A T⁴ with T_ambient = 0 K, ε = 0.40 (`drama-AA7075`). |
| SESAM thermal | Lumped (one temperature), c_p(T) and k(T) tables of `drama-AA7075`; the committed `AA7075_nomelt` material holds those curves at their 850 K values above 850 K. |
| Reference runs | `data/reference_runs/sphere_d100.00mm_..._h077.500km_mAA7075_nomelt_nowind.{csv,json}` and `..._d050.00mm_..._h115.000km_mAA7075_nomelt_nowind.{csv,json}`: SESAM on its static US76 table, winds off, no mass loss; columns `convective_heat_W, rad_cooling_W, integrated_heat_J, temperature_K` (lumped). The Step 1 model on the same table matches their trajectories to 0.24 % / 18 m and 0.06 % / 4 m. |
| Matting (1971) | Convective bridging, Eq. (18): q/q_c = Γ(n)⁻¹ γ(n, [Γ(n+1) q_FM/q_c]^1/n), γ = lower incomplete gamma; a function of the ratio q_FM/q_c, not of Kn; correct limits, no overshoot. n = 1 recommended for hypersonic stagnation regions: q = q_c [1 − exp(−q_FM/q_c)]; n = 2: q = q_c [1 − (1+x)e^-x], x = √(2 q_FM/q_c); n = 1.5 fitted subsonic sphere data. q_FM = A_cq ρ∞ V∞ (h_s − h_w), Eq. (1), with A_cq = 0.8 in his validation. |
| Thermal scales | AA7075: α = k/(ρ c_p) ≈ 5e-5 m²/s; a 50 mm radius has a conduction time R²/α ≈ 50 s and Bi ≈ 1 at peak heating; heating varies on ~10 s scales. |

## 3. Scope

In: gmsh graded tetrahedral sphere mesh with tagged surface patches; AA7075 material with
temperature-dependent k, c_p and an enthalpy hook; two heating modes (SESAM-equivalent, physics);
radiation from every patch; 3D FEM conduction with two solver backends (scikit-fem primary,
FEniCSx optional behind a switch); lockstep coupling with the Step 1 trajectory (trajectory
refactored into a stepper); history, VTK series, animation, stills; comparison plots and metrics
against SESAM; analytic solver tests; conformance and cross-check tests for the backends; CLI;
README.

Out (Step 3 or later): melting and mass loss (latent-heat hook present, zero), surface recession,
the prism-layer surface mesh, film/spray, tumbling or attitude change, non-spherical shapes,
ablation blowing, shock-layer radiation, oxidation heating, a CFD/DSMC heating table (interface
slot reserved), MPI runs of the FEniCSx backend (supported by construction, not exercised).

## 4. Architecture

Package additions (SI inside; plain arrays over surface patches at every boundary; the Step 1
modules unchanged in meaning):

```
reentry_model/mesh.py                 gmsh sphere volume mesh (graded), loader for any tagged gmsh mesh, SurfaceMesh of facet patches
reentry_model/material.py             AA7075: rho, k(T), c_p(T), eps; enthalpy H(T) with latent-heat hook (zero here); c_p_eff
reentry_model/heating.py              HeatingModel protocol; SesamEquivalentHeating; PhysicsHeating (stagnation, bridging, distribution); TabulatedHeating slot
reentry_model/gas.py                  equilibrium-air stagnation state and transport for Fay–Riddell (Cantera airNASA9), Sutton–Graves/DKR fits
reentry_model/thermal/__init__.py     ThermalSolver protocol, factory by name
reentry_model/thermal/skfem_backend.py   scikit-fem P1 tets, backward Euler + Newton, SciPy direct or pyamg CG
reentry_model/thermal/fenicsx_backend.py dolfinx implementation of the same protocol (lazy import)
reentry_model/body.py                 Body protocol as a stepper; ConstantBody (inert); ThermalBody (mesh + material + solver + radiation bookkeeping)
reentry_model/coupled.py              lockstep driver: trajectory macro step -> aero state -> heating -> thermal step -> sample; VTK series
reentry_model/trajectory.py           Simulator.advance(dt) stepper; run() kept as the uncoupled wrapper
reentry_model/viz.py                  PyVista off-screen frames -> MP4/GIF + stills
reentry_model/cli.py                  run/compare options for the thermal model
reentry_model/compare.py              heating/temperature metrics and three new plots
```

Dependency direction: `coupled` imports `trajectory`, `heating`, `body`, `mesh`, `material`;
`heating` imports `gas`, `aero`; `thermal/*` import only `mesh`, `material` and their own library;
`viz` and `compare` import nothing from the solvers (they read the exported fields).

The `Body` protocol becomes

```
mass(t) -> float
advance(t, dt, loads) -> None          # loads: per-patch q_conv, freestream, aero state
surface_temperature() -> np.ndarray    # per patch
mean_temperature() -> float            # energy-equivalent (section 6.4)
energy() -> float                      # stored enthalpy above the reference state
field() -> np.ndarray | None           # nodal temperatures (None for ConstantBody)
```

`ConstantBody` keeps mass and a constant temperature (Step 1 behaviour, `--thermal none`).

## 5. Mesh (`mesh.py`)

- gmsh OCC `addSphere(0,0,0,R)`; a Distance-from-boundary/Threshold size field grading element size
  from `h_surface` at the surface to `h_core` at the centre; linear tetrahedra; the outer boundary is
  one physical surface. Defaults for the 100 mm sphere: `h_surface` = 1.0 mm, `h_core` = 8 mm
  (≈ 40 k nodes); `--mesh-size s` multiplies both. Generated in-process, written as `.msh`, cached
  by (R, h_surface, h_core) under `reentry_model_output/meshes/`.
- `SurfaceMesh`: for every boundary triangle its centroid, outward unit normal, area and node ids;
  `angles_to(v_hat)` gives θ per patch each step. Heating and radiation are applied per facet.
- The loader accepts any gmsh volume mesh with one tagged outer surface (Step 3 prism layers and
  other shapes are drop-ins). Node coordinates are stored as a mutable array (later recession).
- Convergence is a test: halving `h_surface` changes the surface-temperature history by < 1 %.

## 6. Physics

### 6.1 Heating interface
`HeatingModel.evaluate(freestream, V, patches, T_wall, radius) -> q_conv` (W/m², positive into the
body, one value per patch). Inputs per macro step: ρ∞, T∞, m̄, airspeed V, Kn, Ma, nose radius R,
the wall temperature of every patch; θ per patch from the velocity vector.

### 6.2 SESAM-equivalent mode (`--heating sesam`)
q_stag,c = 1.1035e8 R^-1/2 (ρ/1.225)^1/2 (V/7925)^3.15; q_stag,fm = ½ α ρ V³ with α = 1;
q_stag = (1 − f)·q_stag,c + f·q_stag,fm with f = `aero.SesamTable()(Kn)`; **q = 0.27471·q_stag on
every patch, front and back**; no hot-wall factor.

*This distribution is a verification device, not a physical model.* It reproduces SESAM's
tumbling-average assumption so that our resolved conduction, time stepping, material curves and
coupling can be checked against SESAM's lumped temperature and heat totals with no distribution
question in between. The class docstring, a comment at the line applying the factor, and the README
must say so (the README paragraph is part of the acceptance of this step).

### 6.3 Physics mode (`--heating physics`)
- Stagnation, continuum: Fay–Riddell
  q_s = 0.763 Pr^-0.6 (ρ_w μ_w)^0.1 (ρ_s μ_s)^0.4 √(du_e/dx) (h_s − h_w) [1 + (Le^0.52 − 1) γ_cat h_D/h_s],
  du_e/dx = (1/R) √(2 (p_s − p∞)/ρ_s), Pr = 0.71, Le = 1.4; the stagnation state from a normal shock
  followed by isentropic compression to rest in equilibrium air (Cantera, `airNASA9.yaml`), which also
  gives μ_s, μ_w (at T_w) and the dissociation enthalpy fraction h_D/h_s from the equilibrium
  composition; γ_cat ∈ [0, 1] scales the recombination term (1 fully catalytic — SESAM's assumption;
  0 non-catalytic). `--stagnation sutton-graves` (1.7415e-4 √(ρ/R) V³) and `dkr` are the cheap
  alternatives; all three carry the hot-wall factor (1 − h_w/h_s) where the correlation is cold-wall.
- Stagnation, free-molecular: q_fm = A_cq ρ V (h_s − h_w) (Matting Eq. 1, ≈ ½ A_cq ρ V³); A_cq default
  0.8, range 0.8–1.0 (`--accommodation`).
- Bridging: Matting Eq. (18) with `--matting-n` (default 1 → q = q_c [1 − exp(−q_fm/q_c)]; general n
  via `scipy.special.gammainc`); `--bridging-heat sesam-table` selects SESAM's f(Kn) instead.
- Distribution: Lees (1956) laminar sphere with finite Mach, q(θ)/q_s = 2θ sinθ [(1−ε)cos²θ + ε]/√D,
  ε = 1/(γ M²), D = (1−ε)(θ² − ½θ sin4θ + ⅛(1 − cos4θ)) + 4ε(θ² − θ sin2θ + ½(1 − cos2θ)), for θ ≤ 90°,
  zero leeward; in the free-molecular limit the shape is cosθ; the applied shape is
  (1 − w)·Lees + w·cosθ with the free-molecular weight w = 1 − q_stag/q_c under Matting bridging
  (w → 0 in continuum where q_stag → q_c, w → 1 in free-molecular flow where q_stag ≪ q_c; for n = 1,
  w = exp(−q_fm/q_c)) and w = f(Kn) under SESAM bridging.
- `TabulatedHeating` (reserved): q_stag and q(θ)/q_stag interpolated from an offline CFD/DSMC table
  over (V, ρ, T_w); interface defined, not implemented.

### 6.4 Radiation and conduction
- q_rad = ε σ (T_w⁴ − T_amb⁴) per patch, ε = 0.40 (`--emissivity`), T_amb default 0 K (`--t-ambient`,
  200 K optional).
- Conduction: ρ c_p(T) ∂T/∂t = ∇·(k(T) ∇T) in the volume; −k ∇T·n = q_rad − q_conv on the surface.
  Weak form on P1 tets: ∫ ρ c_p Ṫ v + ∫ k ∇T·∇v = ∫_Γ (q_conv − ε σ (T⁴ − T_amb⁴)) v.
- Material (`material.py`): ρ = 2813; k(T), c_p(T) from `data/user_materials/AA7075_nomelt.json`
  (DRAMA's `drama-AA7075` curves, held at 850 K values above 850 K); enthalpy H(T) = ∫ρ c_p dT +
  ρ L_f f_l(T); effective c_p = c_p + L_f df_l/dT; f_l ≡ 0 in this step (Step 3 sets the melting
  range). Initial condition: uniform T₀.
- Outputs per step: nodal T; surface T per patch; **energy-equivalent mean temperature** T_eq defined
  by H(T_eq) = (1/V) ∫ H(T) dV — the temperature of a lumped body holding the same enthalpy, the
  quantity compared with SESAM's lumped `temperature_K`; total convective power Σ q_conv A_patch;
  total radiated power; integrated absorbed heat ∫(Q_conv − Q_rad) dt; stagnation-point and
  back-point surface temperatures; surface max/min/mean.

## 7. Coupling and numerics (`coupled.py`, `trajectory.py`)

- Macro step Δt = 0.5 s default (`--dt`). Per step: `Simulator.advance(Δt)` integrates the trajectory
  from t to t + Δt with DOP853 at the Step 1 tolerances (ground/escape events truncate the final step
  and end the run) → aero state at the end of the step → heating with that state and the current
  wall temperatures → radiation → `ThermalBody.advance` → history sample and (every `--frames-every`
  steps) VTK output. First-order operator splitting; the Δt-halving test bounds its error.
- `trajectory.py` refactor: `advance(dt)` performs one macro step from the stored state (a fresh
  `solve_ivp` call per step; continuity of the state; events as before); `run()` becomes the
  uncoupled wrapper (inert body, `cadence` sampling) and must reproduce the current results to
  < 0.01 m/s and < 0.1 m; the CSV/JSON writers and `read_history_csv` are unchanged.
- Thermal step: backward Euler; k(T), c_p(T) lagged from the previous Newton iterate (Picard), the
  radiation term by Newton with its Jacobian 4εσT³ on the boundary mass — the system stays symmetric
  positive definite; 2–4 iterations to 1e-6 relative; consistent mass by default (`--lumped-mass`).
  Linear solver: SciPy sparse direct (`--linear-solver direct`, default up to ~50 k nodes) or CG with
  a pyamg smoothed-aggregation preconditioner (`amg`).
- History columns: the Step 1 columns (unchanged meaning; `mass_kg` constant) plus
  `convective_heat_W, rad_cooling_W, integrated_heat_J, temperature_K` (= T_eq), `surface_T_max_K,
  surface_T_min_K, surface_T_mean_K, T_stagnation_K, T_back_K, q_stag_Wm2, heating_blend_f`;
  sampled every macro step. Existing column names keep the wrapper's meaning so `compare` and the
  plots work unchanged.
- Cost target: 100 mm reference at the default mesh ≈ 750 steps × (assembly ≈ 0.2 s + 2–3 solves
  ≈ 0.1 s) ≈ 4–6 min; pymsis is called once per macro step.

## 8. Solver backends (`thermal/`)

- Protocol `ThermalSolver`: `setup(mesh, material, emissivity)`, `step(dt, q_conv, T_amb) -> T`,
  `energy() -> float`, `temperature() -> np.ndarray`; factory `thermal_solver(name)`.
- `skfem_backend` (default `--thermal-solver skfem`): scikit-fem `MeshTet` + `ElementTetP1`,
  `Basis`/`FacetBasis` on the tagged boundary, assembly per Newton iterate, solvers as in §7.
- `fenicsx_backend` (`--thermal-solver fenicsx`): the same weak form in UFL on the same gmsh mesh
  (`dolfinx.io.gmshio`), backward Euler with `dolfinx.nls.petsc.NewtonSolver`, PETSc LU or CG +
  hypre/gamg; serial by default, MPI-parallel under `mpirun` by construction (not exercised here).
  `dolfinx` is imported lazily; selecting the backend without it exits 2 with a message naming the
  `fenicsx_env` environment.
- Conformance: both backends pass the same analytic tests (§10) and a cross-check on one mesh and
  load history to 0.1 %. FEniCSx tests are skipped automatically when `dolfinx` is not importable.

## 9. Outputs and visualization

- Run directory `reentry_model_output/<name>/`: `<name>.csv`, `<name>.json` (settings incl. heating
  mode and parameters, mesh parameters and node count, solver, Δt; results incl. peak surface and
  mean temperatures with times/altitudes, total absorbed heat, run time), `field.pvd` + `field_<k>.vtu`
  (nodal T; surface q_conv, q_rad, T), `frames/` and `animation.mp4` (GIF fallback), `stills/`.
- Animation (`--animate`): PyVista off-screen, one frame per `--frames-every` macro steps; sphere
  coloured by surface temperature on a fixed scale (T₀ … run maximum); three-quarter front camera;
  title with t, altitude, velocity; colour bar; stills at start, peak heating, peak temperature,
  impact. Headless rendering of one frame is a unit test.
- `compare` with a SESAM reference: the six Step 1 plots plus `heating_time.png` (Q_conv model vs
  SESAM, residual), `temperature_time.png` (T_eq, surface max/min, SESAM lumped), `integrated_heat.png`;
  metrics: max/rms relative error of Q_conv over the hypersonic phase, relative error of integrated
  heat at the end of the hypersonic phase and at impact, max |T_eq − T_SESAM| and its relative value,
  max/rms relative error of the radiated power.

## 10. Verification and acceptance

**Thermal solver, analytic (unit tests, coarse mesh, both backends):**
1. Lumped limit: k × 1e6, uniform flux → T_eq vs the lumped ODE, 0.1 %.
2. Radiative cooling of an isothermal sphere: no flux, k → ∞ → T(t) = [T₀⁻³ + 3εσA t/(m c_p)]⁻¹ᐟ³, 0.1 %
   (constant c_p).
3. Transient conduction: surface temperature stepped and held (Dirichlet variant of the solver used
   only in the test) → centre temperature vs the Carslaw–Jaeger series, 1 %, constant properties.
4. Energy balance each step: ΔH = (Q_conv − Q_rad) Δt to 1e-6 relative.
5. Mesh halving < 1 % on the surface-temperature history; Δt halving < 0.5 % on peak surface T.
6. Backend cross-check 0.1 %.

**Heating model (unit tests):** DKR, Sutton–Graves and Fay–Riddell at the 100 mm start state
(ρ = 2.727e-5 kg/m³, V = 7.5 km/s, R = 0.05 m) against hand-computed values: DKR = 1.96e6 W/m²,
Sutton–Graves = 1.72e6 W/m² (ratio 0.875), Fay–Riddell within 15 % of Sutton–Graves (its fit);
Lees windward integral 0.196 (M→∞) / 0.209 (M=5); Matting limits and n = 1 closed form
vs the incomplete-gamma form; SESAM-equivalent total = 0.27471·q_stag·4πR² to 1e-6; the
free-molecular shape cosθ integrates to 0.25.

**Coupling:** inert-body regression vs Step 1 (< 0.01 m/s, < 0.1 m); `reference`-marked comparisons
of the coupled run (`--heating sesam`, `--atmosphere us76`) against both US76 SESAM references:
- Q_conv within 3 % (max relative error over the hypersonic phase),
- integrated heat within 3 % at the end of the hypersonic phase,
- T_eq within 2 % (in kelvin) of SESAM's lumped temperature over the whole flight,
- radiated power within 3 %,
and the physics-mode run reported (integrated-heat ratio to SESAM, expected ≈ 0.7, with the
0.196/0.2747 and Fay–Riddell/DKR explanation) but not thresholded. As in Step 1, thresholds are
expectations: the plan's last task measures, records the numbers in the README, and changes a
threshold only with the measured value and the reason beside it.

## 11. CLI

`run` adds: `--thermal none|fem` (default `none`: Step 1 behaviour), `--heating sesam|physics`,
`--stagnation fay-riddell|sutton-graves|dkr`, `--bridging-heat matting|sesam-table`, `--matting-n`,
`--accommodation`, `--catalycity`, `--emissivity`, `--t-ambient`, `--mesh-size`, `--h-surface`,
`--h-core`, `--thermal-solver skfem|fenicsx`, `--linear-solver direct|amg`, `--lumped-mass`, `--dt`,
`--frames-every`, `--animate`, `--stills`. Exit codes as in Step 1 (2 for a missing backend).
`compare` accepts the new columns and produces the three new plots when present.

## 12. Environment

`drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`, Python 3.12, arm64) gains
`scikit-fem`, `gmsh`, `pyamg`, `pyvista`, `imageio[ffmpeg]`, `cantera` (all pip wheels available).
FEniCSx: a separate conda-forge environment `fenicsx_env` (`fenics-dolfinx`, `mpich`, `gmsh`,
plus this package's pip dependencies), created only after explicit go-ahead; its tests run with that
interpreter and are skipped elsewhere.

## 13. Assumptions to state in the thesis

1. Fixed attitude (windward face toward the velocity vector) in physics mode; SESAM's tumbling
   average reproduced only in the verification mode.
2. Laminar, equilibrium boundary layer; catalycity as a parameter, fully catalytic by default.
3. No ablation blowing, no oxidation heating, no shock-layer radiation (negligible for R ≤ 0.1 m
   below 9 km/s).
4. Material properties of DRAMA's `drama-AA7075` table; behaviour above 850 K by holding the last
   tabulated value (the no-melt reference material).
5. Radiation to a black background at T_amb (0 K to match SESAM; 200 K optional); ε = 0.40 constant.
6. First-order operator splitting with Δt = 0.5 s; backward Euler in time; P1 tetrahedra.
7. The SESAM-equivalent heating distribution is uniform over the surface by construction and is not
   a physical model (see §6.2).
