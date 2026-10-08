# ---------------------------------------------------------------------------------------------------------- Task 5
emit(r'''### Task 5: The flow-regime gate, the wall pressure and the surface flow per patch

**Files:**
- Modify: `reentry_model/gas.py` (four edits)
- Create: `reentry_model/surface_flow.py`
- Test: `tests/test_reentry_model_surface_flow.py`

**Interfaces:**
- Consumes: `gas.EquilibriumAir.stagnation` (Step 2), `aero.mean_free_path` (SESAM's own lambda, measured fact 17), `aero.SesamTable`, `trajectory.AeroState` (freestream rho/T/p/m_bar, V, kn, ma, a_drag).
- Produces: `GasState.s`, `.a`, `.m_bar` (defaulted fields), `EquilibriumAir.expand(stag, p)`; `surface_flow.SurfaceFlow(air=None, rarefied_shear="slip"|"bridged", bridging=None, sigma_v=1, sigma_t=1, gamma_pm=1.15, kn_body_shock=0.01, kn_body_fm=10.0)` with `.branch_of(state)`, `.edge_table(state, stag, radius)`, `.evaluate(state, theta, radius, rho_liquid, T_wall=None) -> SurfaceFlowResult`, `.last_bins`; `SurfaceFlowResult` carrying the body-scale diagnostics `branch, kn_body, re_shock, mach_inf, p_stag, p_inf, phi_sonic, deceleration, gamma_pm` and the per-patch arrays `p_w, u_e, u_eff, rho_e, T_e, mu_e, mach_e, delta_a, lambda_w, kn_local, closure, flagged, tau, tau_continuum, tau_fm, G`, with `.branch_name` and `.closure_fractions(areas, windward)`; functions `ranger_psi`, `ranger_thickness`, `boundary_layer_thickness`, `prandtl_meyer(mach, gamma)`, `mach_from_turn(turn, gamma)`, `sonic_angle(p_stag, p_inf, q_inf, gamma)`, `wall_pressure(theta, p_stag, p_inf, q_inf, gamma, blend_deg)`, `mean_free_path_maxwell`, `wall_knudsen(p_w, T_wall, mu_wall, length, R_s=R_SPECIFIC_AIR)`; constants `BRANCH_SHOCK_LAYER/MERGED/FREE_MOLECULAR = 0/1/2`, `BRANCH_NAMES`, `CLOSURE_GIRIN/COUETTE = 0/1`, `KN_BODY_SHOCK = 0.01`, `KN_BODY_FM = 10.0`, `KN_LOCAL_CONTINUUM = 0.01`, `KN_LOCAL_SLIP = 0.1`, `RE_SHOCK_MERGED = 100.0`, `R_SPECIFIC_AIR`, `RANGER_C = 58.08`, `THWAITES_C = 0.45`, `GAMMA_PM = 1.15`, `PM_BLEND_DEG = 10.0`, `PM_MAX_DEG = 110.0`, `RAREFIED_SHEAR_NAMES`, `THETA_BINS` (0-90 deg by 1 deg).

The module docstring carries the physics and every threshold's justification -- the standoff argument for 0.01, the Re2 cross-check, the declared conservatism and its measured size, why the Knudsen length is geometric and the state is the wall's, why the Prandtl-Meyer blend must start at phi*, and why no oblique-shock machinery is built. Transcribe it verbatim: it is the spec's section 18 written in the code.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_surface_flow.py`:

''')
emit(file_block("tests/test_reentry_model_surface_flow.py"))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_surface_flow.py -q`
Expected: ImportError (`reentry_model.surface_flow`).

- [ ] **Step 3: Extend `reentry_model/gas.py`**

Five edits. (a) After `AIR = "N2:0.79, O2:0.21"` add `AVOGADRO = 6.02214076e23`. (b) In `GasState`, after `X: dict           # mole fractions > 1e-6` add:

```python
    s: float = 0.0    # J/(kg K), specific entropy (Cantera basis) -- Step 3: the isentropic expansion to the edge state
    a: float = 0.0    # m/s, frozen sound speed
    m_bar: float = 0.0  # kg, mean molecular mass (Step 3: the edge mean free path)
```

(c) In `_state()` replace the `return GasState(...)` statement with:

```python
        return GasState(float(g.P), float(g.T), float(g.density), float(g.enthalpy_mass),
                        wilke_viscosity(float(g.T), X, self.molar_masses), float(h_D), X, float(g.entropy_mass),
                        float(g.sound_speed), float(g.mean_molecular_weight) * 1e-3 / AVOGADRO)
```

(d) Before `def wall(self, T_w, p):` insert:

```python
    def expand(self, stag, p):
        """Isentropic expansion of the stagnation state `stag` to pressure p (equilibrium above T_EQUILIBRATE):
        the boundary-layer edge state (Step 3, spec section 7)."""
        g = self.gas
        g.SPX = stag.s, p, AIR
        if g.T > T_EQUILIBRATE:
            g.equilibrate("SP")
        return self._state()

```

(e) Nothing else changes; `"$PY" -m pytest tests/test_reentry_model_gas.py tests/test_reentry_model_heating.py -q` must still pass.

- [ ] **Step 4: Create `reentry_model/surface_flow.py`**

''')
emit(file_block("reentry_model/surface_flow.py"))
emit(r'''
- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_surface_flow.py tests/test_reentry_model_gas.py tests/test_reentry_model_heating.py -q`
Expected: 7 passed plus the Step 2 gas/heating tests. `evaluate` at 71 km takes ~55 ms (one shock solve, 91 Cantera isentropic expansions; the wall state is analytic). The branch test asserts the measured values there: merged branch, Kn_body 0.0113, Re2 153, every closure Couette and flagged, Kn_local 1.4e-4 at the nose (Kn_body/84), and p_w(89 deg) between 150 and 350 Pa instead of p_inf.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/gas.py reentry_model/surface_flow.py tests/test_reentry_model_surface_flow.py
git commit -m "Add the flow-regime gate, the Newtonian+Prandtl-Meyer wall pressure and the wall-Knudsen melt closure (Step 3 Task 5)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 6
emit(r'''### Task 6: The melt film — lubrication and runoff

**Files:**
- Create: `reentry_model/film.py`
- Test: `tests/test_reentry_model_film.py`

**Interfaces:**
- Consumes: `SurfaceMesh.edges/normals/centroids/areas` (Task 1).
- Produces: `lubrication(tau, G, b, delta_m, mu_l, b_layer=None) -> (V_s, q, shear_rate, thick)` -- `b` is the mobile film, `b_layer` the contiguous liquid depth the branch test uses when given (fact 28); `Runoff(surface, points, windward=None)` with `.i, .j, .length, .n_i, .n_j, .n_patches`, `.edge_coefficients(q, b, t_hat, areas) -> (c_ij, c_ji)`, `.transport(m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4) -> (m_f, n_solves, moved_kg)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_film.py`:

''')
emit(file_block("tests/test_reentry_model_film.py"))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_film.py -q`
Expected: ImportError (`reentry_model.film`).

- [ ] **Step 3: Create `reentry_model/film.py`**

''')
emit(file_block("reentry_model/film.py"))
emit(r'''
- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_film.py -q`
Expected: 3 passed (the strip's interior reaches the closed-form thickness to 1e-3 exactly at the scheme's steady state).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/film.py tests/test_reentry_model_film.py
git commit -m "Add the melt film: lubrication branches and the linearly implicit upwind runoff (Step 3 Task 6)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 7
emit(r'''### Task 7: Spraying — instability branches, release bookkeeping, published reference values

**Files:**
- Create: `reentry_model/spray.py`, `data/reference_values/girin2017_table1.json`, `data/reference_values/girin1994_tables.json`
- Test: `tests/test_reentry_model_spray.py`

**Interfaces:**
- Consumes: `dispersion.DispersionTable` (Task 4), `surface_flow.BRANCH_FREE_MOLECULAR`, `CLOSURE_GIRIN` and `SurfaceFlowResult` fields (Task 5), `material.LiquidProperties` (Task 2).
- Produces: `melt_layer(flow, liquid) -> (delta_m, factor)` -- **nan/0 wherever the patch's closure is not `CLOSURE_GIRIN`**, since Eq. (2) and Ranger's Psi(phi) are continuum constructions and Task 5's gate decides whether they may be used at all; `rayleigh_taylor(deceleration, b, liquid) -> (active, lambda*, tau*)` (the unbounded fastest mode, reported only; `deceleration` is the component **normal to the film**, Girin & Kopyt's W sin(Theta) = W cos(phi) on this body, per patch or scalar); `rayleigh_taylor_bounded(deceleration, h, extent, liquid) -> (active, lambda, tau)`, the same mode in a pool of finite lateral extent, which admits only k >= pi/extent, with `extent` the contiguous **molten region** each patch belongs to (facts 30 and 33; reported only); `wave_fits(lam, extent) -> bool` -- one whole wavelength must fit inside that region, applied on every branch; `thin_film_mode(mach, momentum_flux, liquid) -> (lambda*, tau*)`; `SprayModel(liquid, k_r=0.17, k_t=1.1, we_critical=4.62, table=None).evaluate(flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, b_layer=None, extent=None, deceleration_n=None) -> SprayResult(branch, we_s, unstable, r, mdot, dm, dn, rt_active, delta_m, v_s, growth)` (`b_layer` is the contiguous liquid depth that decides the thick/thin branch, `extent` the molten region's lateral extent for the wave-fits test, `deceleration_n` the surface-normal deceleration for the reported Rayleigh-Taylor flag, and `growth` the selected mode's growth time per patch). Girin's criterion `We_s > We_cr` with We_s = rho_l V_s^2 min(delta_m, layer)/Sigma gates **every** branch, not just the thick one (fact 32); Girin's regime test then runs in **two stages**, in his order (fact 36). Stage one asks whether the *theoretical* melt boundary layer delta_m of his Eq. (2) can form in the liquid present: where it cannot, layer <= delta_m, the rigid core still stabilises the disturbances and the case is his regime 1, the Girin & Kopyt (1994) thin branch. Where it has formed, stage two decides the mechanism on the kinematic viscosities nu_melt = mu_l/rho_l against nu_gas = mu_e/rho_e -- `BRANCH_REGIME2` (regime 2, classical Kelvin-Helmholtz on the tangential-discontinuity profile, reusing Girin & Kopyt's lambda* and tau* with their rigid-wall factor cth(Lambda) -> 1 and the thick branch's torus rate rho_l pi r^2/(lambda* tau*), r = k_r lambda*, which does not reference the film depth) where the melt is the more viscous medium, and the thick branch (regime 3, his gradient instability) where the gas is. Liquid aluminium sits 3.3e4 to 6.9e6 from that threshold, so regime 2 is inert on these flights and only the unit tests exercise it; `source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state) -> list of 22-value rows`; `histogram(r, dn, dm) -> (dn per bin, dM per bin)`; `SprayModel(..., rt_spray=True)` applies the front-surface Rayleigh-Taylor mode on the patches where it passes both criteria and grows faster than the shear mode (fact 35; `rt_spray=False` reports it only, as every run before 2026-09-24 did). Constants `K_R, K_T, B_MIN, N_BINS = 40, R_MIN = 1e-6, R_MAX = 1e-2, BIN_EDGES, WE_BREAKUP = 12, CAPILLARY_TAU, SOURCE_COLUMNS, BRANCH_THICK/THIN/RAREFIED/RT = 0/1/2/3, BRANCH_REGIME2 = 4` (Girin thick / thin-film on the edge state / thin-film on the freestream / front-surface Rayleigh-Taylor / Girin & Kopyt's step profile for a deep film whose melt is the more viscous medium). The `closure` column replaces `regime` in `SOURCE_COLUMNS`. Spraying is never switched off by rarefaction: the gate selects the closure, and hence whether delta_m exists to define a thick film at all.

- [ ] **Step 1: Write the two reference-value files**

`data/reference_values/girin2017_table1.json`:

''')
emit(block(read("data/reference_values/girin2017_table1.json"), "json"))
emit(r'''
`data/reference_values/girin1994_tables.json`:

''')
emit(block(read("data/reference_values/girin1994_tables.json"), "json"))
emit(r'''
- [ ] **Step 2: Write the failing tests**

Create `tests/test_reentry_model_spray.py`:

''')
emit(file_block("tests/test_reentry_model_spray.py"))
emit(r'''
- [ ] **Step 3: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_spray.py -q`
Expected: ImportError (`reentry_model.spray`).

- [ ] **Step 4: Create `reentry_model/spray.py`**

''')
emit(file_block("reentry_model/spray.py"))
emit(r'''
- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_spray.py -q`
Expected: 12 passed (the plan said 4 until 2026-09-25; the file has grown with the gates of fact 32, the applied Rayleigh-Taylor mode of fact 35 and the three regime-2 tests of fact 36).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/spray.py data/reference_values tests/test_reentry_model_spray.py
git commit -m "Add melt spraying: Girin's thick, thin and rarefied branches, release bookkeeping and the published reference values (Step 3 Task 7)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 8
emit(r'''### Task 8: Girin's published cases — the Girin-as-published driver and its analysis script

**Files:**
- Create: `reentry_model/girin_case.py`, `analysis/girin_reference.py`
- Test: `tests/test_reentry_model_girin.py`

**Interfaces:**
- Consumes: `dispersion.DispersionTable` (Task 4), `surface_flow.ranger_psi` (Task 5), `spray.thin_film_mode/rayleigh_taylor/histogram/BIN_EDGES` (Task 7), `compare`'s plot style (Step 2), the two reference-value files (Task 7).
- Produces: `girin_case.GirinVariant(name, V, rho_m, mu_m, sigma, R0=3e-3)`, `dimensionless_numbers(v, rho_a=RHO_A, mu_a=MU_AIR, compression=COMPRESSION) -> (Re, We, GI, alpha, mu, p, t_ch)`, `critical_angle(GI, p, we_cr=4.62)`, `surface_state(...)`, `run_variant(v, k_r, k_t, we_cr, re_density="shock"|"ambient", ...) -> dict` (keys `GI, phi_cr_deg_eq3, t_f_min_us_tau0, r_um_90deg_tau0, N, r_med_um, r_min_um, r_max_um, t_sd_ms, tau_end, radii, counts, times, mass_history, ...`); constants `RHO_A, MU_AIR, COMPRESSION, RE_DENSITY_NAMES`. The script writes `girin2017.json`, `girin1994.json`, `girin_summary.md` and the plots under `reentry_model_output/verification_melt/girin/`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_girin.py`:

''')
emit(file_block("tests/test_reentry_model_girin.py"))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_girin.py -q`
Expected: ImportError (`reentry_model.girin_case`).

- [ ] **Step 3: Create `reentry_model/girin_case.py`**

''')
emit(file_block("reentry_model/girin_case.py"))
emit(r'''
- [ ] **Step 4: Create `analysis/girin_reference.py`**

''')
emit(file_block("analysis/girin_reference.py"))
emit(r'''
- [ ] **Step 5: Run the tests and the script**

```bash
"$PY" -m pytest tests/test_reentry_model_girin.py -q
"$PY" analysis/girin_reference.py
```

Expected: 5 passed; the script prints the summary table (variant I ambient: GI 13.04 (13.00, ×1.00), φ_cr 16.3 (16.1), t_f 7.0 (5.7, ×1.22), N 1.31e6 (1.50e6, ×0.87), r_med 25.8 (26.9, ×0.96), t_s.d. 2.82 (5.90, ×0.48); 1994 Table 1 shock-layer factor 7.77 with all r_d/τ_d/ṁ ratios 0.99–1.01; Table 2 λ*/Eq. (14) 1.00–1.02, λ*/printed 10.0, τ* 1.00) and writes eleven files under `reentry_model_output/verification_melt/girin/`.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/girin_case.py analysis/girin_reference.py tests/test_reentry_model_girin.py
git commit -m "Reproduce Girin's published spraying cases with his own simplifications (Step 3 Task 8)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')
