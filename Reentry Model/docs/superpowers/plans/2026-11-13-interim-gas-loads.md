# Interim Gas Loads Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the sphere flight's surface pressure, wall shear and convective heating with the interim gas loads of
Part A of the design: a similarity-based boundary layer (shear and heating from one solution, unreacted post-shock
gas) below Kn_body 0.05, Daniels & Schwartzentruber's shear above it, and real loads on the leeward side, all in
`reentry_model/shared/` so both architectures use them.

**Architecture:** Three new shared modules: `daniels_schwartzentruber` (the 2026 sphere model), `similarity` (frozen-air
properties, the compressible similarity solver and the march round the meridian) and `gas_loads` (the replaceable
provider: one `GasLoads.meridian(...)` per step gives pressure, shear and heating on 0-180 degree bins, cached so the
surface flow and the heating share one solve). `gas.py` gains an unreacted (chemically frozen) post-shock gas.
`SurfaceFlow` and a new `SimilarityHeating` read the provider; with no provider (`--gas-loads classic`, and every state
saved before this change) the code paths are today's, byte for byte.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy (`solve_bvp`, `CubicSpline`), Cantera (`airNASA9.yaml`), pytest.

**Spec:** `docs/superpowers/specs/2026-11-12-shoulder-melt-loss-design.md`, Part A (sections A1-A7, with A4b, A6).
Reference calculations, data and plots: `analysis/gas_loads/` (README there); test fixtures: `tests/fixtures/gas_loads/`.

## Global Constraints

- All new physics lives in `reentry_model/shared/` and imports nothing from the rest of `reentry_model` except
  package-root constants; `tests/test_reentry_model_shared.py` must pass, with the three new modules added to its list.
- Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; run tests from `Reentry Model/` with
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest ... -p no:cacheprovider`.
- `--gas-loads classic` must reproduce `main` (08cca55) byte for byte: `analysis/skyscrapers_port/same_results.py
  compare` exit 0 on the cases `step1 skins element nomelt` (Task 8).
- The droplet stage is untouched: `droplet_reentry` keeps building `PhysicsHeating` itself; its reproduction test
  (`reentry_model/droplet_reentry/tests/test_droplet_reentry_model.py`, rel 1e-9) must still pass.
- Existing `SurfaceFlowResult` and `HeatingResult` field names and meanings are kept (the level-set stage-3 branch
  reads `tau`, `dp_ds`, ...). `PhysicsHeating` is not changed.
- Pinned values: shear switch at `Kn_body = 0.05` (`state.kn`, `aero.knudsen`, hard sphere d = 3.65 angstrom);
  similarity Pr 0.71; frozen air `"N2:0.79, O2:0.21"`; separation where f''(0) <= 0.005 or no converged solution;
  first solve by a wall-temperature homotopy of 12 geometric steps from 0.8 T_e down to T_w; meridian bins every
  2 degrees from 0 to 180 by default (1 degree for the reference tests); Matting n = 1, accommodation 0.8; Prandtl-Meyer
  gamma = `SurfaceFlow.gamma_pm` (default 1.15), continued past 90 degrees uncapped and joined to Daniels &
  Schwartzentruber's lee pressure at the first angle above 90 degrees where it is at or below theirs; corners the
  centroid leeward test shelters are evaluated at theta = pi (separated wake).
- Daniels & Schwartzentruber as in the reference (`analysis/gas_loads/reference/ds_model.py`): perfect gas gamma 1.4,
  omega 0.75, Pr 0.72, Sutherland viscosity at T_inf, Kn by their (A 1), fully diffuse wall (sigma_t = sigma_n = 1),
  theta_s capped at 180 degrees, (2.29) as printed, the composite rarefied formula (A 26)-(A 27) at every Kn.
- Cost: one `GasLoads.meridian` call at Kn_body about 0.003 conditions takes at most 2 s on the laptop after its first
  call (target; the test bound is 5 s).
- Run names: `--gas-loads classic` adds `_loads-classic`; the default adds nothing.
- States saved before this change resume on the classic loads (their pickled `SurfaceFlow` has no `gas_loads`
  attribute; the class attribute `gas_loads = None` keeps them classic; their `PhysicsHeating` stays as pickled).
- Commit messages: a sentence in the imperative or a scope prefix and a sentence describing the behaviour that changed,
  ending with the line `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A late, slow flight state (Mach about 1-3, small total enthalpy): `GasLoads.meridian` must return finite loads with
   no exception (the homotopy and Prandtl-Meyer continuation at low Mach). Test in Task 4.
2. Free-molecular (Kn_body >= 10) and degenerate states (rho or V <= 0): no solve is attempted where none is needed, no
   exception, heating tends to Matting's free-molecular limit. Tests in Tasks 5 and 6.
3. Patch angles at exactly 0, pi/2 and pi, unsorted and repeated, with a different wall temperature on every patch:
   loads interpolated at the patch angles, the hot-wall factor per patch. Tests in Tasks 5 and 6.
4. A solve that fails part-way round the meridian is separation, not an exception, and the next step starts clean
   (the warm start is reset). Test in Task 3.
5. Saving and resuming: a run built with `GasLoads` pickles and resumes; a state whose `SurfaceFlow` lacks `gas_loads`
   evaluates exactly as classic. Tests in Tasks 4 and 7.

---

### Task 1: Daniels & Schwartzentruber sphere loads

**Files:**
- Create: `reentry_model/shared/daniels_schwartzentruber.py`
- Test: `tests/test_reentry_model_daniels_schwartzentruber.py`
- Modify: `tests/test_reentry_model_shared.py` (add the module to its list of shared modules)

**Interfaces:**
- Produces:
  - `mu_sutherland(T: float) -> float` [Pa s]
  - `@dataclass(frozen=True) class DSCondition: re_inf, mach_inf, kn_inf, tw_t0, tw_tinf: float; gamma=1.4; omega=0.75; pr=0.72; sigma_t=1.0; sigma_n=1.0`
  - `condition(rho_inf, T_inf, V, diameter, T_wall, gamma=1.4, omega=0.75, R=287.05) -> DSCondition` (Re and Kn by the
    paper's definitions, perfect-gas T0)
  - `parameters(cond: DSCondition) -> dict[str, float]` with the keys of `ds_params.csv` (`Re, Kn, M, T0, Tw_T0, Re_s,
    ..., Br, w_fm`)
  - `coefficients(theta: np.ndarray, cond: DSCondition) -> tuple[np.ndarray, np.ndarray]`: (Cp, Cf) at polar angles in
    radians, Cp = (p - p_inf)/(0.5 rho U^2), Cf = tau/(0.5 rho U^2), composite regime
  - `drag_coefficient(cond: DSCondition) -> float`

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_daniels_schwartzentruber.py`:
  - `test_reproduces_reference_distributions`: for the four Holman-Boyd conditions (rho 1.007e-4, 2.014e-5, 3.987e-6,
    8.057e-7 kg/m3; T_inf 200 K; V 7108 m/s; D 0.3048 m; T_w 1500 K) `coefficients(np.radians(theta), condition(...))`
    equals `tests/fixtures/gas_loads/ds_model.csv` columns `Cp`, `Ctau` at every 1-degree angle, `rtol=1e-9, atol=1e-12`.
  - `test_reproduces_reference_parameters`: `parameters(...)` equals `ds_params.csv` row by row for every column the row
    has except `Cd_HB_DSMC`, `rtol=1e-6`.
  - `test_drag`: `drag_coefficient` is 1.0702, 1.236, 1.534, 1.796 (`rtol=2e-3`) for the four conditions.
  - `test_separation_angle_capped`: `parameters(...)["theta_s"] <= 180` for a Kn 0.5 condition.
- [ ] **Step 2: Run them; expect `ModuleNotFoundError`.**
  `"$PY" -m pytest tests/test_reentry_model_daniels_schwartzentruber.py -q`
- [ ] **Step 3: Implement the module** by porting the formulas of `analysis/gas_loads/reference/ds_model.py`
  (`kn_from`, `post_shock`, `theta_ic_s`, `angle`, `parameters`, `continuum`, `free_molecular`, `distribution`, `drag`,
  `taw_T0`, `mu_sutherland`) behind the interface above, angles in radians at the interface; no file I/O.
- [ ] **Step 4: Run the tests and the shared-import test; expect PASS.**
  `"$PY" -m pytest tests/test_reentry_model_daniels_schwartzentruber.py tests/test_reentry_model_shared.py -q`
- [ ] **Step 5: Commit** (`reentry_model: Daniels & Schwartzentruber (2026) sphere pressure and shear in the shared package ...`).

### Task 2: Unreacted post-shock gas

**Files:**
- Modify: `reentry_model/shared/gas.py` (`EquilibriumAir`, around :66-160)
- Test: `tests/test_reentry_model_gas.py` (create if absent; else add to it)

**Interfaces:**
- Produces: `class UnreactedAir(EquilibriumAir)`: the same `stagnation(rho_inf, T_inf, V)`, `expand(stag, p)`,
  `wall(T_w, p)`, `freestream(rho, T)` and `air_enthalpy(T)`, but the composition stays `AIR` throughout (no
  `equilibrate` call anywhere); `EquilibriumAir` gets an instance attribute `t_equilibrate` (default `T_EQUILIBRATE`,
  used where the module constant is used today) and `UnreactedAir` sets it to `math.inf`.

- [ ] **Step 1: Write the failing tests:**
  - `test_unreacted_stagnation_at_holman_boyd_kn_0_002`: `UnreactedAir().stagnation(1.007e-4, 200.0, 7108.0)` has
    T = 16411 K (`rel=1e-3`), rho = 1.023e-3 kg/m3 (`rel=2e-3`), p = 4837 Pa (`rel=2e-3`), mu = 3.28e-4 Pa s
    (`rel=5e-3`), and its mole fractions are those of `AIR`.
  - `test_equilibrium_unchanged`: `EquilibriumAir().stagnation(1.007e-4, 200.0, 7108.0)` has T 5989 K, rho 1.767e-3,
    p 4947 Pa (`rel=1e-3`), exactly as before this change (also compare to a value computed with `t_equilibrate` left at
    its default).
  - `test_unreacted_expand_keeps_composition`: `expand(stag, 0.5 * stag.p)` keeps `AIR` and has `h < stag.h`.
- [ ] **Step 2: Run; expect FAIL (`UnreactedAir` undefined).**
- [ ] **Step 3: Implement** `t_equilibrate` and `UnreactedAir` in `gas.py`.
- [ ] **Step 4: Run these tests and `tests/test_reentry_model_heating.py tests/test_reentry_model_surface_flow.py`; expect PASS** (classic users unchanged).
- [ ] **Step 5: Commit.**

### Task 3: Similarity boundary-layer solver

**Files:**
- Create: `reentry_model/shared/similarity.py`
- Test: `tests/test_reentry_model_similarity.py`
- Modify: `tests/test_reentry_model_shared.py` (module list)

**Interfaces:**
- Consumes: `gas.UnreactedAir`, `gas.wilke_viscosity` (Task 2).
- Produces:
  - `class FrozenAirProperties(air: gas.EquilibriumAir)`: cubic splines on T = 150 to 20000 K every 50 K of frozen-air
    h(T) (Cantera, `"N2:0.79, O2:0.21"`), T(h), mu(T) (`wilke_viscosity`) and their derivatives; methods `T_of_h(h)`,
    `h_of_T(T)`, `mu_of_T(T)`, `dmu_of_T(T)`, `cp_of_T(T)`.
  - `class ConstantProperties`: the same methods for a calorically perfect gas with mu proportional to T (so C = 1 and
    T/T_e = h/h_e); for the Falkner-Skan checks.
  - `@dataclass class SimilarityPoint: fpp0: float; gp0: float; g_w: float`
  - `class SimilaritySolver(props, pr=0.71)`: `solve(beta, H_e, u_e, T_w) -> SimilarityPoint`, warm-started from its
    last solution, the first call (and the first after `reset()`) by the 12-step wall-temperature homotopy;
    `reset()`. The equations, unknowns, boundary conditions and heat-flux formula are those of
    `analysis/gas_loads/reference/similarity_heat.py` (`SimilarityHeat._solve`), eta in [0, 12] on 500 points,
    `solve_bvp(tol=1e-4, max_nodes=200000)`; a failed solve raises `SimilarityError`.
  - `@dataclass class MeridianSolution: theta_deg, beta, fpp0, gp0, tau_c [Pa], q_c [W/m2]: np.ndarray; separation_deg: float` (nan if attached to the last bin)
  - `march(theta_deg, radius, p_w, u_e, H_e, T_w, rho_w_mu_w, solver) -> MeridianSolution`: xi on a 20-times refined
    grid, beta from `np.gradient(u_e, s)`, the stagnation point by beta = 1/2 with u_e r / sqrt(2 xi) -> sqrt(2 beta0 /
    (rho_w mu_w)), each bin warm-started from the previous; the first bin with f''(0) <= 0.005 or a `SimilarityError`
    is the separation angle, and tau_c = q_c = 0 there and beyond. Same arithmetic as `similarity_heat.py` `main()`.

- [ ] **Step 1: Write the failing tests:**
  - `test_falkner_skan`: with `ConstantProperties`, `g_w = 1` (T_w = T_e) and u_e = 0, `solve(beta, ...)` gives
    f''(0) = 0.9277, 0.4696, 1.2326 for beta = 0.5, 0.0, 1.0 (`abs=5e-4`).
  - `test_cold_wall_stagnation`: `FrozenAirProperties`, beta 0.5, H_e 25162715.7 J/kg, u_e 0, T_w 1500 K gives
    f''(0) = 0.5024 and g'(0) = 0.3313 (`abs=2e-3`); T_w 16395 K gives f''(0) = 0.9273 (`abs=1e-3`).
  - `test_march_reproduces_reference`: for each of the four Holman-Boyd conditions, with `p_w = p_over_pinf * p_inf`
    from `tests/fixtures/gas_loads/similarity_lee.csv` on 1-degree bins, the unreacted stagnation state and edge
    (`UnreactedAir.expand`), T_w 1500 K, `march(...)` gives `tau_c/(0.5 rho U^2)` and `q_c/(0.5 rho U^3)` equal to the
    `Ctau_continuum` and `Cq_continuum` columns of `tests/fixtures/gas_loads/similarity_heat.csv` (`rtol=1e-6` where
    the fixture is non-empty) and the same separation angle (151, 151, 151, 158 degrees).
  - `test_failed_solve_is_separation`: monkeypatch `SimilaritySolver.solve` to raise `SimilarityError` from the bin at
    100 degrees on; `march` returns separation 100 and zeros beyond, raises nothing; after `reset()` a new march solves
    from the homotopy again (Review Focus 4).
- [ ] **Step 2: Run; expect `ModuleNotFoundError`.**
- [ ] **Step 3: Implement `similarity.py`.**
- [ ] **Step 4: Run the tests and the shared-import test; expect PASS.**
- [ ] **Step 5: Commit.**

### Task 4: The gas-loads provider and the leeward pressure

**Files:**
- Create: `reentry_model/shared/gas_loads.py`
- Test: `tests/test_reentry_model_gas_loads.py`
- Modify: `tests/test_reentry_model_shared.py` (module list)

**Interfaces:**
- Consumes: Tasks 1-3; `surface_flow.wall_pressure`, `surface_flow.mach_from_turn`, `surface_flow.GAMMA_PM`.
- Produces:
  - `KN_SHEAR_SWITCH = 0.05`; `BIN_DEG = 2.0`
  - `joined_pressure(theta_deg, p_stag, p_inf, q_inf, gamma_pm, p_ds) -> tuple[np.ndarray, float]`: the model's
    `wall_pressure` to 90 degrees; beyond it the same Prandtl-Meyer expansion uncapped (no p_inf floor); from the first
    angle above 90 degrees where it is <= `p_ds` (Daniels & Schwartzentruber's pressure on the same bins), `p_ds`;
    returns (p, meet angle in degrees, nan if they never meet). Same as `similarity_lee.joined_pressure`.
  - `@dataclass class MeridianLoads`: `theta_deg, p_w, dp_ds, tau, q_c` (n_bins,) arrays over 0-180 degrees;
    `q_c0` (stagnation similarity flux, W/m2), `H_e`, `h_w` (J/kg), `T_wall`, `kn_body`, `shear_source`
    (`"similarity"` or `"daniels-schwartzentruber"`), `meet_deg`, `separation_deg`, `p_stag`, `p_inf`, `q_inf`;
    method `at(theta_rad, field) -> np.ndarray` (linear interpolation in angle; theta clipped to [0, pi]).
    `tau` here is the CONTINUUM shear below the switch (the caller applies slip and blend) and Daniels &
    Schwartzentruber's shear (already rarefied) at and above it.
  - `class GasLoads(gamma_pm=GAMMA_PM, kn_switch=KN_SHEAR_SWITCH, bin_deg=BIN_DEG)`:
    `meridian(state, radius, T_wall) -> MeridianLoads` for a trajectory `AeroState`, cached on
    `(state.freestream.rho, state.freestream.T, state.V, radius, T_wall)` so the surface flow and the heating share one
    solve per step; it holds an `UnreactedAir`, `FrozenAirProperties` and one `SimilaritySolver`; `__getstate__` drops
    the cache and the solver's warm start (rebuilt on first use after loading). The pressure law is evaluated on the
    UNREACTED stagnation pressure (design A4b, last bullet); Daniels & Schwartzentruber from
    `condition(rho, T_inf, V, 2 radius, T_wall)`. The similarity solution (for `q_c`, and for `tau` below the switch) is
    computed at every Kn_body; rho or V <= 0 returns all-zero loads with `p_w = p_inf`.

- [ ] **Step 1: Write the failing tests:**
  - `test_joined_pressure_reproduces_reference`: on 1-degree bins at the four Holman-Boyd conditions (unreacted
    stagnation pressure, `gamma_pm` 1.15, `p_ds` from `ds_model.csv`), `joined_pressure / p_inf` equals
    `similarity_lee.csv` `p_over_pinf` (`rtol=1e-9`) and the meet angle is 150, 150, 150, 157.
  - `test_meridian_matches_reference_on_fine_bins`: `GasLoads(bin_deg=1.0).meridian(state, 0.1524, 1500.0)` for the
    Kn 0.002 Holman-Boyd state reproduces `similarity_heat.csv` (`Ctau_continuum`, `Cq_continuum`, `rtol=1e-6`) and
    `shear_source == "similarity"`.
  - `test_default_bins_close_to_fine`: with the default 2-degree bins, `tau` and `q_c` interpolated at 0-140 degrees are
    within 3 % of the 1-degree result and the separation angle within 2 degrees.
  - `test_switch_to_daniels_schwartzentruber`: for the Kn 0.05 Holman-Boyd state (Kn_body about 0.067) `shear_source ==
    "daniels-schwartzentruber"` and `tau / q_inf` equals `ds_model.csv` `Ctau` at the bins (`rtol=1e-9`).
  - `test_cache_shared`: two `meridian` calls with the same arguments return the same object; a different `T_wall`
    recomputes.
  - `test_low_speed_state`: V 600 m/s at 30 km (Mach about 2): finite `p_w, tau, q_c`, no exception (Review Focus 1).
  - `test_pickle_round_trip`: `pickle.loads(pickle.dumps(loads))` evaluates the same `meridian` (Review Focus 5).
  - `test_meridian_cost`: after one call, a second call at a slightly different state takes < 5 s (print the time).
- [ ] **Step 2: Run; expect `ModuleNotFoundError`.**
- [ ] **Step 3: Implement `gas_loads.py`.**
- [ ] **Step 4: Run; expect PASS; note the measured cost in the commit message.**
- [ ] **Step 5: Commit.**

### Task 5: The surface flow reads the provider

**Files:**
- Modify: `reentry_model/shared/surface_flow.py` (`SurfaceFlow.__init__` :268, `evaluate` :300-375)
- Test: `tests/test_reentry_model_surface_flow.py`

**Interfaces:**
- Consumes: `GasLoads`, `MeridianLoads` (Task 4).
- Produces: `SurfaceFlow.__init__(..., gas_loads=None)`; class attribute `gas_loads = None`. With `gas_loads is None`
  `evaluate` is today's code, unchanged. With a provider, in the shock-layer and merged branches:
  - `p_w` everywhere from `L.at(theta, "p_w")`; `dp_ds` everywhere from `L.at(theta, "dp_ds")`;
    `G = -dp_ds - rho_liquid * W * sin(theta)` everywhere (leeward included).
  - below the switch: `tau_continuum = L.at(theta, "tau")` (both sides); `tau` from it through today's slip factor and
    branch logic (shock layer: free-molecular shear where `kn_local >= KN_LOCAL_SLIP`; merged: `(1-f) tau_slip + f
    tau_fm`, `tau_fm` windward only as today); `kn_local` on the lee from the lee `p_w`.
  - at and above the switch (and in the free-molecular branch): `tau = L.at(theta, "tau")` both sides,
    `tau_continuum = tau`; the free-molecular branch keeps its `p_w = p_inf`.
  - every other field (equilibrium edge state, `delta_a`, `closure`, `flagged`, `u_eff`, `mach_e`, `shear_integral`)
    as today; `T_wall` for the provider is `T_wall[argmin(theta)]` (300 K when `T_wall` is None).
- [ ] **Step 1: Write the failing tests** (existing classic tests stay as they are and must keep passing):
  - `test_provider_leeward_pressure_and_shear`: at the Kn 0.002 Holman-Boyd state, theta 120 degrees: `p_w` equals the
    provider's joined pressure, `tau > 0`; at 170 degrees `tau == 0` (beyond separation) and `p_w` equals Daniels &
    Schwartzentruber's lee pressure; `dp_ds` is non-zero at 120.
  - `test_provider_windward_pressure_is_the_law_on_unreacted_stagnation`: `p_w` at 0-90 equals `wall_pressure(theta,
    p_stag_unreacted, p_inf, q_inf, 1.15)`.
  - `test_provider_rear_stagnation`: theta = pi gives `tau == 0` and `p_w` = the Daniels & Schwartzentruber base pressure.
  - `test_provider_above_switch_uses_ds_shear`: at the Kn 0.25 state `tau / q_inf` equals `ds_model.csv` `Ctau`.
  - `test_unsorted_repeated_angles`: theta `[pi, 0, pi/2, 0.3, 0.3]` gives the same values as the sorted call, element
    by element (Review Focus 3).
  - `test_old_state_is_classic`: a provider-built `SurfaceFlow` with its instance `gas_loads` deleted evaluates exactly
    as `SurfaceFlow()` (Review Focus 5).
  - `test_free_molecular_and_degenerate`: Kn_body 20 and rho = 0 return finite results, no exception (Review Focus 2).
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run `tests/test_reentry_model_surface_flow.py tests/test_reentry_model_spray.py`; expect PASS.**
- [ ] **Step 5: Commit.**

### Task 6: Similarity heating

**Files:**
- Modify: `reentry_model/shared/heating.py` (add a class; `heating_by_name`)
- Test: `tests/test_reentry_model_heating.py`

**Interfaces:**
- Consumes: `GasLoads.meridian` (Task 4), `matting_bridge`, `cosine_shape`, `HeatingResult`.
- Produces: `class SimilarityHeating(gas_loads, bridging="matting", matting_n=1.0, accommodation=0.8)`, `name =
  "physics"`, with `evaluate(state, theta, T_wall, radius, T_mean=None) -> HeatingResult` (same signature as
  `PhysicsHeating.evaluate`; `T_mean` ignored): `L = gas_loads.meridian(state, radius, T_wall[argmin(theta)])`;
  `q_c = L.q_c0`; `q_fm = accommodation * rho * V * (L.H_e - L.h_w)`; `q_stag = matting_bridge(q_fm, q_c, n)`;
  `w = 1 - q_stag/q_c`; `q = q_stag * ((1 - w) * L.at(theta, "q_c")/q_c + w * cosine_shape(theta)) * hot`, with the
  per-patch hot-wall factor `hot = (L.H_e - h_w(T_wall))/(L.H_e - L.h_w)` (frozen-air h_w); rho or V <= 0 returns the
  zero result as `PhysicsHeating` does. `heating_by_name("physics", gas_loads=...)` returns it when `gas_loads` is given.
- [ ] **Step 1: Write the failing tests:**
  - `test_reproduces_reference_heating`: at the four Holman-Boyd states, uniform T_w 1500 K, 1-degree provider,
    `q_conv/(0.5 rho U^3)` equals `similarity_heat.csv` `Cq` (`rtol=1e-6` where non-empty; zero where empty).
  - `test_leeward_heating_until_separation`: at Kn 0.002, q > 0 at 120 degrees and 0 at 160 and at pi.
  - `test_hot_wall_per_patch`: two patches at the same angle with T_w 1000 and 1500 K: ratio = `(H_e - h_w(1000))/(H_e
    - h_w(1500))` (Review Focus 3).
  - `test_free_molecular_limit`: at Kn_body 20 `blend > 0.99` and `q_stag` within 1 % of `q_fm` (Review Focus 2).
  - `PhysicsHeating` tests unchanged and passing.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run `tests/test_reentry_model_heating.py`; expect PASS.**
- [ ] **Step 5: Commit.**

### Task 7: Wiring: sheltered corners, the `--gas-loads` option, saving and resuming

**Files:**
- Modify: `reentry_model/skyscrapers/skinbody.py` (`load_angles` :910-919)
- Modify: `reentry_model/cli.py` (options near :419-425 and :559-587, `build_thermal` :611-716, `model_run_name`
  :254-276, the JSON settings :692-712)
- Test: `tests/test_reentry_model_skinbody.py` (:2442-2583), `tests/test_reentry_model_cli.py` (:1233-1243 and new),
  `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: `GasLoads` (Task 4), `SurfaceFlow(gas_loads=...)` (Task 5), `SimilarityHeating` (Task 6).
- Produces:
  - `SkinBody.load_angles()` under `centroid` returns `np.where(windward_corners(), theta, np.pi)` (decision 2a:
    sheltered corners in the separated wake); under `facing`, `theta` as today.
  - CLI `--gas-loads {similarity,classic}`, default `similarity`. `build_thermal` builds one `GasLoads(gamma_pm=
    args.gamma_pm)` when similarity and passes it to every `SurfaceFlow` it builds and, with `--heating physics`, to
    `SimilarityHeating` (with `--bridging-heat matting`, `--matting-n`, `--accommodation`; `--stagnation` and
    `--catalycity` are classic-only and recorded as unused); `--heating sesam` is unchanged. `classic` builds exactly
    what it builds today.
  - The run JSON's settings record `"gas_loads": "similarity"|"classic"`; `model_run_name` appends `_loads-classic`
    for classic.
- [ ] **Step 1: Write the failing tests:**
  - skinbody: replace the 90 deg + 1e-6 expectations (:2442-2464) by pi; the sheltered-corner tests (:2467-2583) now
    assert `tau == 0`, `dp_ds` finite, and zero heating at a sheltered corner under the provider; under `--gas-loads
    classic` the same corner still gets zero loads.
  - cli: `--gas-loads` in `BASE_RUN_OPTION_STRINGS`; `model_run_name` ends `_loads-classic` only for classic;
    `build_thermal` with similarity gives a `SurfaceFlow` and a `SimilarityHeating` sharing one `GasLoads` object; with
    classic, `PhysicsHeating` and `SurfaceFlow().gas_loads is None`; `--gas-loads bogus` exits 2.
  - restart: a short similarity run saved and resumed continues with `GasLoads` (Review Focus 5); a state whose
    `SurfaceFlow` instance lacks `gas_loads` resumes on classic loads.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run `tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py tests/test_reentry_model_coupled.py tests/test_reentry_model_melting.py -n 4`; expect PASS.** If a qualitative bound in coupled or melting fails because the physics changed, report the old and new values to the controller instead of loosening it silently.
- [ ] **Step 5: Commit.**

### Task 8: Verification and documentation

**Files:**
- Create: `analysis/gas_loads/plots/model_vs_dsmc.py` (the model's own loads, through `SurfaceFlow` and
  `SimilarityHeating` with a `GasLoads`, at the four Holman-Boyd conditions against `reference/hb_dsmc.csv`, written to
  `output/gas_loads/model_vs_dsmc.png`; same layout as `plots/windward_loads_comparison.py`, 0-180 degrees)
- Modify: `docs/model_assumptions.md` (new section "Gas loads (interim, 2026-11-13)"; the statements that heating is
  zero leeward in sections 2 and 8), the root `README.md` (:307-399, "zero leeward" at :351-352), `CLAUDE.md` (the
  `--heating physics` description under "Things that are easy to get wrong", and the shared-module list),
  `analysis/gas_loads/README.md`, the design's status line.

- [ ] **Step 1: Classic is byte-identical to main.** Make a detached checkout of 08cca55 (`git worktree add --detach
  "$SCRATCH/main-08cca55" 08cca55`), then from `Reentry Model/`:
  `"$PY" analysis/skyscrapers_port/same_results.py run --code "$SCRATCH/main-08cca55/Reentry Model" --out "$SCRATCH/sr_main" --cases step1,skins,element,nomelt`
  `"$PY" analysis/skyscrapers_port/same_results.py run --code . --out "$SCRATCH/sr_classic" --cases step1,skins,element,nomelt --append "--gas-loads classic"`
  `"$PY" analysis/skyscrapers_port/same_results.py compare "$SCRATCH/sr_main" "$SCRATCH/sr_classic"` -> exit 0.
  (The run names differ by `_loads-classic`; if `compare` pairs files by name, rename the classic folder's files
  first and say so.) Remove the detached worktree afterwards.
- [ ] **Step 2: The new default flies.** Run the `skins` and `element` cases with no `--gas-loads` into
  `"$SCRATCH/sr_similarity"`; both exit 0. Report the wall time per macro step against the classic runs and the
  history's stagnation heating, shoulder shear and the share of heat landing behind 90 degrees.
- [ ] **Step 3: Regenerate the comparison from the model.** Write and run `model_vs_dsmc.py`; check by eye that the
  model's curves lie on the reference's (`similarity_heat.csv`, `ds_model.csv`) and copy the figure into
  `docs/superpowers/specs/figures/2026-11-12-shoulder-melt-loss/model_vs_dsmc.png`. Also measure the jump at the
  switch (design A4b): `SurfaceFlow` shear at 30, 60 and 90 degrees for the Holman-Boyd freestream with rho set so
  that Kn_body is 0.0495 and 0.0505; report the ratios (not smoothed).
- [ ] **Step 4: The unit tier and the droplet tests.**
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`
  and `"$PY" -m pytest reentry_model/droplet_reentry/tests -q`; all pass (report the counts).
- [ ] **Step 5: The Step 2 reference tier** (about 35 minutes; run in the background):
  `"$PY" -m pytest tests/test_reentry_model_reference_thermal.py -m reference -q`. If the physics-mode ratio to SESAM
  moves (it was 0.744 and 0.770, bounds 0.4-1.1 in `PHYSICS_RATIO`), record the new values in
  `analysis/reentry_model_thermal_verification.py`, `docs/model_assumptions.md` and `README.md` together; if a bound
  fails, stop and report rather than widening it.
- [ ] **Step 6: Documentation.** The new model_assumptions section states, each with its measurement or marked
  assumed: the similarity solution (equations, unreacted gas, frozen-air properties, separation rule), the
  Daniels-Schwartzentruber switch at Kn_body 0.05, the joined leeward pressure, the sheltered-corner rule, the known
  accuracy against Holman-Boyd DSMC from design A6, the deferred nose-heating item A-5, and `--gas-loads classic`.
- [ ] **Step 7: Commit** (`docs: ...`), then report to the controller: the same_results result, the smoke-flight costs
  and loads, the unit-tier counts, the reference-tier ratios.
