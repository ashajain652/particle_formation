# Coupled trajectory + 3D FEM heat transfer of a sphere (Step 2 of the physics model) — Design

Date: 2026-09-18
Status: implemented (plan docs/superpowers/plans/2026-09-18-thermal-fem.md); amended with the measurements listed in section 14
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
| SESAM heating | Total convective power Q = A_sphere · 0.27471 · q_DKR · F_h(Kn) · max(0, 1 − c_p(T − T∞)/(V²/2)) for Ma ≥ 1 (½ · 0.27471 · q_DKR below Ma 1), q_DKR = 1.1035e8 R^-1/2 (ρ/1.225)^1/2 (V/7925)^3.15 W/m², c_p = 1004.5 J/kg/K, F_h measured in 0.125-decade Kn bins (aero.SesamHeatTable: 1.005 continuum, 0.14 at Kn 1, 0.059 at Kn 40 = 0.78 × the free-molecular cos θ average). Supersedes the (1 − f)·q_DKR + f·q_FM decomposition of facts §4 (measured 2026-09-18, facts §15). |
| SESAM shape factor | 0.27471 = surface-average ÷ stagnation flux, constant in Mach in the ATDB; applied to the whole surface area of a lumped mass under the "random tumbling" attitude assumption. Hot-wall factor 1 − c_p(T − T∞)/(V²/2), clamped at 0 (facts §15). |
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
  one physical surface. Defaults: `h_surface` = 2.0 mm, `h_core` = 8 mm (18.9 k nodes on the 100 mm sphere; 1 mm/8 mm is 76 k
  nodes and is the convergence mesh); a node is embedded at the centre; `--mesh-size s` multiplies both. Generated in-process, written as `.msh`, cached
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
q_DKR = 1.1035e8 R^-1/2 (ρ/1.225)^1/2 (V/7925)^3.15; F_h(Kn) the measured bridging (`aero.SesamHeatTable`,
facts §15); hot-wall factor max(0, 1 − c_p(T_eq − T∞)/(V²/2)) for Ma ≥ 1 (½, no hot-wall term, below Ma 1),
c_p = 1004.5 J/kg/K, using the body's energy-equivalent temperature T_eq (`T_mean`) exactly as SESAM's
lumped model does; **q = 0.27471 · q_DKR · F_h(Kn) · hot-wall on every patch, front and back**.

*This distribution is a verification device, not a physical model.* It reproduces SESAM's
tumbling-average assumption so that our resolved conduction, time stepping, material curves and
coupling can be checked against SESAM's lumped temperature and heat totals with no distribution
question in between. The class docstring, a comment at the line applying the factor, and the README
must say so (the README paragraph is part of the acceptance of this step).

### 6.3 Physics mode (`--heating physics`)
- Stagnation, continuum: Fay–Riddell
  q_s = 0.763 Pr^-0.6 (ρ_w μ_w)^0.1 (ρ_s μ_s)^0.4 √(du_e/dx) (h_s − h_w) [1 + (Le^0.52 − 1) γ_cat h_D/h_s],
  du_e/dx = (1/R) √(2 (p_s − p∞)/ρ_s), Pr = 0.71, Le = 1.4; the stagnation state from a normal shock
  followed by isentropic compression to rest in equilibrium air (Cantera, `airNASA9.yaml`) for the
  thermodynamics and composition; that mechanism has no transport data, so μ_s and μ_w come from
  Blottner's curve fits (N2, O2, NO, N, O) with Wilke's mixing rule; the equilibrium solver is used
  only where T∞ + V²/(2c_p) > 1500 K, which also gives μ_s, μ_w (at T_w) and the dissociation
  enthalpy fraction h_D/h_s from the equilibrium composition; γ_cat ∈ [0, 1] scales the recombination
  term (1 fully catalytic — SESAM's assumption; 0 non-catalytic); γ_cat scales the Lewis-number term
  of the equilibrium form (γ_cat = 0 gives the Le = 1 value; the frozen non-catalytic reduction
  1 − h_D/h_s is a Step 3 option). `--stagnation sutton-graves` (1.7415e-4 √(ρ/R) V³) and `dkr` are the
  cheap alternatives; all three carry the hot-wall factor (1 − h_w/h_s) where the correlation is
  cold-wall; each patch's flux is scaled by its own (h_s − h_w)/(h_s − h_w,stag).
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
  total radiated power; `integrated_heat_J` = ∫Q_conv dt (SESAM's meaning: convective input only, not
  net of radiation); `absorbed_heat_J` = ∫(Q_conv − Q_rad) dt; stagnation-point, back-point and
  `T_centre_K` (nodal temperature at the mesh centre) surface/volume temperatures; surface max/min/mean.

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
- Thermal step: backward Euler; k(T) at the element-mean temperature of the previous iterate and the
  secant heat capacity [h(T_k) − h(T_old)]/(T_k − T_old), so the discrete energy balance is exact to
  the Newton tolerance for the tabulated c_p; the radiation term by Newton with its Jacobian 4εσT³ on
  the boundary mass — the system stays symmetric positive definite; 2–4 iterations to 1e-6 relative;
  consistent mass by default (`--lumped-mass`). Linear solver: CG with a pyamg smoothed-aggregation
  preconditioner rebuilt every 30 solves (`amg`, default; 0.03 s per solve at 19 k nodes) or SciPy
  SuperLU (`direct`, 0.4–24 s per solve at 12–76 k nodes, tests only); Newton starts from the
  extrapolated previous step (2.0–2.1 iterations per step).
- History columns: the Step 1 columns (unchanged meaning; `mass_kg` constant) plus
  `convective_heat_W, rad_cooling_W, integrated_heat_J, absorbed_heat_J, temperature_K` (= T_eq),
  `surface_T_max_K, surface_T_min_K, surface_T_mean_K, T_stagnation_K, T_back_K, T_centre_K,
  q_stag_Wm2, heating_blend_f`; sampled every macro step. Existing column names keep the wrapper's
  meaning so `compare` and the plots work unchanged.
- Measured (Task 12, 2026-09-18, this machine): 100 mm reference, default mesh, 732 steps, ~158 s
  (SESAM-equivalent), ~182 s including the surface-temperature animation (physics mode; both modes
  run the same 732 steps to the ground); 50 mm from 115 km, 1117 steps, ~47 s (SESAM-equivalent),
  ~54 s (physics mode, same 1117 steps to the ground). Both heating modes reach the ground on both
  references (facts §16: the equilibrium-air shock iteration behind Fay–Riddell used to stall near
  Mach 1 during descent; fixed by skipping the shock at Ma ≤ 1.1, where none exists).

## 8. Solver backends (`thermal/`)

- Protocol `ThermalSolver`: `setup(mesh, material, emissivity)`, `step(dt, q_conv, T_amb) -> T`,
  `energy() -> float`, `temperature() -> np.ndarray`; factory `thermal_solver(name)`.
- `skfem_backend` (default `--thermal-solver skfem`): scikit-fem `MeshTet` + `ElementTetP1`,
  `Basis`/`FacetBasis` on the tagged boundary, unit element stiffness/mass matrices precomputed once
  and rescaled into a fixed CSR pattern per iterate (20 ms vs 230 ms for scikit-fem's generic `asm`
  at 12.6 k nodes, which is kept as the reference operator in the conformance test); radiation from
  the facet-mean temperature with its exact Jacobian; solvers as in §7.
- `fenicsx_backend` (`--thermal-solver fenicsx`): the same weak form in UFL on the same gmsh mesh
  (`dolfinx.io.gmshio`), backward Euler with `dolfinx.nls.petsc.NewtonSolver`, PETSc LU or CG +
  hypre/gamg; serial by default, MPI-parallel under `mpirun` by construction (not exercised here).
  `dolfinx` is imported lazily; selecting the backend without it exits 2 with a message naming the
  `fenicsx_env` environment. Written against dolfinx 0.9/0.10; untested until `fenicsx_env` exists;
  lumped mass not implemented there.
- Conformance: both backends pass the same analytic tests (§10) and a cross-check on one mesh and
  load history to 0.1 %. FEniCSx tests are skipped automatically when `dolfinx` is not importable.

## 9. Outputs and visualization

- Run directory: `<outdir>/<name>.csv/.json` as in Step 1; plots, `vtk/` (field.pvd + field_<k>.vtu,
  surface.pvd + surface_<k>.vtp), `vtk/animation.mp4`, `vtk/frames/`, `vtk/stills/` under
  `<outdir>/<name>/`.
- Animation (`--animate`): PyVista off-screen, one frame per `--frames-every` macro steps; sphere
  coloured by surface temperature on a fixed scale (T₀ … run maximum); three-quarter front camera;
  title with t, altitude, velocity; colour bar; stills at start, peak heating, peak temperature,
  impact. Headless rendering of one frame is a unit test.
- `compare` with a SESAM reference: the six Step 1 plots plus `heating_time.png` (Q_conv model vs
  SESAM, residual), `temperature_time.png` (T_eq, surface max/min, SESAM lumped), `integrated_heat.png`;
  metrics: power errors relative to SESAM's peak over the hypersonic phase; point-wise Q_conv error
  over Kn_ref < 0.01 and Q_ref > 10 % of peak; integrated heat at the end of the hypersonic phase and
  at the end; max |T_eq − T_SESAM| and relative.

## 10. Verification and acceptance

**Thermal solver, analytic (unit tests, coarse mesh, both backends):**
1. Lumped limit: k × 1e4, uniform flux → T_eq vs the lumped ODE, 0.1 %.
2. Radiative cooling of an isothermal sphere: no flux, k × 1e4 → T(t) = [T₀⁻³ + 3εσA t/(m c_p)]⁻¹ᐟ³, 0.1 %
   (constant c_p).
3. Transient conduction: surface temperature stepped and held (Dirichlet variant of the solver used
   only in the test) → centre temperature vs the Carslaw–Jaeger series on a uniform 4 mm test mesh
   (centre 0.5 %, volume mean 0.2 %), constant properties.
4. Energy balance each step: ΔH = (Q_conv − Q_rad) Δt to 1e-6 relative.
5. Mesh halving < 1 % on the surface-temperature history; Δt halving < 0.5 % on peak surface T.
6. Backend cross-check 0.1 %.

**Heating model (unit tests):** DKR, Sutton–Graves and Fay–Riddell at the 100 mm start state
(ρ = 2.727e-5 kg/m³, V = 7.5 km/s, R = 0.05 m) against hand-computed values: DKR = 1.96e6 W/m²,
Sutton–Graves = 1.72e6 W/m² (ratio 0.875), Fay–Riddell 2.22e6 W/m² = 1.2–1.4 × Sutton–Graves
(measured 1.29; 1.14 with the Lewis term off); Lees windward integral 0.196 (M→∞) / 0.209 (M=5);
Matting limits and n = 1 closed form vs the incomplete-gamma form; SESAM-equivalent total =
0.27471·q_stag·4πR² to 1e-6; the free-molecular shape cosθ integrates to 0.25.

**Coupling:** inert-body regression vs Step 1 (< 0.01 m/s, < 0.1 m); `reference`-marked comparisons
of the coupled run (`--heating sesam`, `--atmosphere us76`) against both US76 SESAM references:
- Q_conv within 3 % of peak over the hypersonic phase (measured 0.46 % / 0.37 %) and 3 % point-wise
  in the continuum (measured 2.24 % / 2.36 %),
- integrated heat within 3 % at the end of the hypersonic phase (measured +0.31 % / +0.03 %),
- T_eq within 2 % (in kelvin) of SESAM's lumped temperature over the whole flight (measured
  1.20 % / 1.22 %),
- radiated power within 8 % (= 4 × the temperature margin: the resolved surface radiates at its own,
  hotter temperature; measured 6.65 % / 3.40 %),
and the physics-mode run reported (integrated-heat ratio to SESAM at the end of the flight, measured
0.744 / 0.770, matching the 0.196/0.2747 × Fay–Riddell/DKR × hot-wall estimate) but not thresholded.
Both references' `--heating physics` runs now reach the ground (facts §16: the Fay–Riddell shock
iteration, `gas.EquilibriumAir.stagnation`, used to stall in the Mach ≈ 1.00–1.01 band both references
cross during descent; fixed by skipping the shock entirely at Ma ≤ 1.1, where none exists). As in
Step 1, thresholds are expectations: the plan's last task measures, records the numbers in the
README, and changes a threshold only with the measured value and the reason beside it.

## 11. CLI

`run` adds: `--thermal none|fem` (default `none`: Step 1 behaviour), `--heating sesam|physics`,
`--stagnation fay-riddell|sutton-graves|dkr`, `--bridging-heat matting|sesam-table`, `--matting-n`,
`--accommodation`, `--catalycity`, `--emissivity`, `--t-ambient`, `--mesh-size`, `--h-surface`,
`--h-core`, `--thermal-solver skfem|fenicsx`, `--linear-solver direct|amg`, `--lumped-mass`, `--dt`,
`--frames-every`, `--animate`, `--stills`. Exit codes as in Step 1 (2 for a missing backend).
`compare` accepts the new columns and produces the three new plots when present.

## 12. Environment

`drama_env` (`/Users/ashajain/miniforge3/envs/drama_env/bin/python`, Python 3.12, arm64) gains, pinned
in `requirements-step2.txt`: `scikit-fem==12.0.2`, `gmsh==4.15.2`, `meshio==5.3.5`, `pyamg==5.3.0`,
`pyvista==0.49.0`, `imageio==2.37.4`, `imageio-ffmpeg==0.6.0`, `cantera==3.2.0` (all pip wheels).
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

## 14. Amendments (2026-09-18)

Measured facts and spec amendments applied by Task 12 (verification against SESAM; facts note §§15–16):

1. Status line updated to "implemented", pointing at this section.
2. §2 SESAM heating: total power is A·0.27471·q_DKR·F_h(Kn)·hot-wall (max(0, 1 − c_p(T − T∞)/(V²/2))
   for Ma ≥ 1, half that with no hot-wall term below Ma 1); supersedes the (1 − f)·q_DKR + f·q_FM
   decomposition of facts §4.
3. §2 SESAM shape factor: the hot-wall factor is present and quantified (c_p = 1004.5 J/kg/K), clamped
   at 0, not "no hot-wall correction detectable".
4. §5 mesh defaults: `h_surface` = 2.0 mm / `h_core` = 8 mm (18.9 k nodes) is the default; 1.0 mm/8 mm
   (76 k nodes) is the convergence mesh, not the default.
5. §6.2 SESAM-equivalent formula rewritten to the measured q_DKR · F_h(Kn) · hot-wall form, the
   hot-wall factor using the body's energy-equivalent temperature T_eq.
6. §6.3 Fay–Riddell: Cantera supplies thermodynamics/composition only (no transport data); μ_s, μ_w
   from Blottner + Wilke; the equilibrium solver activates only above ~1500 K stagnation enthalpy;
   γ_cat's role in the Lewis-number term clarified; each patch scaled by its own hot-wall ratio.
7. §6.4 outputs: `integrated_heat_J` is ∫Q_conv dt (SESAM's meaning, not net of radiation);
   `absorbed_heat_J` and `T_centre_K` added to the documented outputs.
8. §7 numerics: element-mean k(T) with the secant heat capacity (energy-exact), not a Picard-lagged
   c_p; solver performance (amg 0.03 s/solve at 19 k nodes, direct 0.4–24 s at 12–76 k nodes, 2.0–2.1
   Newton iterations/step) and measured step counts/runtimes for both references recorded.
9. §8 backends: `skfem_backend`'s custom assembly measured at 20 ms vs 230 ms for scikit-fem's generic
   `asm` (12.6 k nodes); `fenicsx_backend`'s dolfinx version target and untested/no-lumped-mass status
   noted.
10. §9: run-directory layout and comparison-metric definitions corrected to match the implementation
    (`vtk/` subdirectory; point-wise continuum Q_conv error; heat at end of hypersonic phase and at
    the end).
11. §10 acceptance: Fay–Riddell/Sutton–Graves ratio measured (1.29; 1.14 with the Lewis term off,
    both at 2.22e6 W/m² catalytic); the Carslaw–Jaeger check uses a uniform 4 mm mesh (measured 0.5 %
    centre, 0.2 % volume mean); the lumped/radiative-cooling checks use k × 1e4, not 1e6; the
    radiated-power threshold widened to 8 % (reason: the resolved surface radiates at its own, hotter
    temperature); the measured values of the README's verification table recorded, including the
    physics-mode integrated-heat ratio at the end of the flight (0.744 / 0.770).
12. §12 environment: `requirements-step2.txt` versions pinned and listed explicitly.
13. §6.3/§10 (fix round 1, 2026-09-18): `gas.EquilibriumAir.stagnation` skips the normal-shock fixed
    point entirely at Ma ≤ 1.1 (no shock exists at Ma ≤ 1, and a Ma < 1.1 shock's entropy jump is
    < 0.1 %) and compresses isentropically from the freestream instead, fixing the stall the fixed
    point hit as its contraction factor → 1 near Ma 1 (facts §16); both references' `--heating physics`
    runs now reach the ground.
