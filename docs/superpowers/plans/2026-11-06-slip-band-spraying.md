# Slip-Band Spraying Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let Girin's spraying run for body Knudsen numbers from 0.01 up to 0.1, with `delta_m` from a momentum balance, replicating what the code does below 0.01.

**Architecture:** The surface flow gains the integrated shear from the nose. A new `spray_layer` returns the existing Eq. (2) depth, plus the momentum-balance depth on slip-band patches when the band is on. Only the spraying calls use it: film transport and the feed keep `melt_layer`, so they do not change.

**Tech Stack:** Python, numpy, pytest, in the planning copy `prototype/work-2026-10-28-skin-spraying-plan/code`.

**Spec:** No spec file. The decisions are Asha's of 2026-11-06, listed below; the derivation is `I = integral of r tau ds` and `delta_m^3 = 3 mu^2 I / (rho tau^2 r)`.

## Asha's decisions (2026-11-06)

1. The body gate stays the model's rule; the band is body Knudsen 0.01 to 0.1, on the merged branch.
2. `delta_m` in the band comes from the momentum balance, not Eq. (2). Eq. (2) stays below 0.01, so `delta_m` jumps by about 2 to 2.5 times at 0.01; accepted.
3. The entire molten depth may spray, to define `b`: the band reads the flowing depth with liquid properties, as the code does below 0.01. The liquid-slurry viscosity-contrast instability is future work.
4. `--slurry-spray on` also lets the slurry spray in the band, by the existing `slurry_modes`. With it off, the release stays capped at the liquid.
5. Shear is the bridged `flow.tau`. The thin branch keeps its edge-state wavelength (a flagged assumption, those values being the bridge's continuum endpoint). The element model gets the band too.
6. The viscous floor applies on the liquid thin branch. The code already applies it (`viscous_film`, `spray.py`, thin loop), so the band pins this with a test and changes no gate.

## Global Constraints

- With `--slip-band off`, the default, every output is identical to the last bit.
- `KN_BODY_SLIP = 0.1`, in `surface_flow.py`; `--slip-band on|off`, default `off`.
- `pickle`d states from Task 5a must still load: new attributes have class-level defaults.
- Commit only when Asha asks. Trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Never hand-edit `22-skin-spraying-flight.md`: once built, this becomes a new task in `part_22.md`, named Task 5c because Task 5b is already taken by the resume command, with its files added to `make_plan22.py`.
- Do not touch the flights, the DRAMA databases or `23-notes.md`.
- `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. Run tests from `code/`.

## Review Focus

- A patch at the nose (`r = 0` or `tau = 0`): `delta_m` is nan, with no warning or division error, and the patch behaves as it does today.
- Leeward patches: `I` is built from windward patches only, and leeward `delta_m` is nan.
- Two patches at the same angle (many columns share a `theta`): no division by zero and `I` is non-decreasing with angle.
- Body Knudsen exactly 0.01 and exactly 0.1: 0.01 is in the band (it is the merged branch), 0.1 is outside.
- A state saved by Task 5a, or a `FakeFlow` with no `shear_integral`, works with the band off.

---

### Task 1: The integrated shear and the momentum-balance depth

**Files:**
- Modify: `reentry_model/surface_flow.py` (`KN_BODY_SLIP`, `shear_moment`, field `SurfaceFlowResult.shear_integral`, set in `evaluate`'s final return)
- Modify: `reentry_model/spray.py` (`momentum_layer`, `spray_layer`, class attribute `SprayModel.slip_band`)
- Test: `tests/test_reentry_model_surface_flow.py`, `tests/test_reentry_model_spray.py`

**Interfaces:**
- Produces: `surface_flow.shear_moment(theta, tau, radius, windward) -> ndarray`, `I` [Pa m2] per patch, nan on leeward patches.
- Produces: `SurfaceFlowResult.shear_integral: ndarray | None = None`, last field; `None` in the early returns; a property `in_slip_band -> bool`, true when `branch == BRANCH_MERGED and kn_body < KN_BODY_SLIP`.
- Produces: `spray.momentum_layer(tau, I, r, liquid) -> ndarray`, `delta_m` [m].
- Produces: `spray.spray_layer(flow, liquid, band=False, theta=None, radius=None) -> ndarray`, `melt_layer`'s depth where `band` is false or `flow.in_slip_band` is false, otherwise `momentum_layer` on windward patches with `r = radius * sin(theta)`.

- [ ] **Step 1: Write the failing tests**

```python
def test_shear_moment_matches_the_closed_form():
    theta = np.linspace(0.0, np.pi / 2, 2001)
    tau = 40.0 * np.sin(theta)
    I = surface_flow.shear_moment(theta, tau, 0.05, theta <= np.pi / 2)
    exact = 40.0 * 0.05 ** 2 * (theta / 2 - np.sin(2 * theta) / 4)
    assert I[1:] == pytest.approx(exact[1:], rel=1e-3)

def test_shear_moment_keeps_order_ties_and_leeward():
    theta = np.array([0.6, 0.2, 0.2, 2.0])
    I = surface_flow.shear_moment(theta, np.array([3.0, 1.0, 1.0, 9.0]), 0.05, theta <= np.pi / 2)
    assert I[1] == I[2] and I[0] > I[1] and np.isnan(I[3])

def test_momentum_layer_values():
    liq = LiquidProperties(2400.0, 0.0013, 0.8)
    d = spray.momentum_layer(np.array([40.0]), np.array([1e-2]), np.array([0.02]), liq)
    assert d[0] == pytest.approx((3 * 0.0013 ** 2 * 1e-2 / (2400 * 40.0 ** 2 * 0.02)) ** (1 / 3))
    # We_s with V_s = tau delta_m / mu reduces to 3 I / (r sigma), free of viscosity
    assert 2400 * (40.0 * d[0] / 0.0013) ** 2 * d[0] / 0.8 == pytest.approx(3 * 1e-2 / (0.02 * 0.8))

def test_momentum_layer_is_nan_at_the_nose():
    liq = LiquidProperties(2400.0, 0.0013, 0.8)
    d = spray.momentum_layer(np.array([0.0, 40.0]), np.array([0.0, 1e-2]), np.array([0.0, 0.0]), liq)
    assert np.isnan(d).all()

def test_spray_layer_is_melt_layer_outside_the_band():
    # FakeFlow with no shear_integral, band off; and in_slip_band False with band on
    assert_array_equal(spray.spray_layer(flow, liq), spray.melt_layer(flow, liq)[0])
```

Add to the surface-flow tests: with a real `evaluate`, `in_slip_band` is false at Knudsen 0.009 and 0.1, and true at 0.01 and 0.05.

- [ ] **Step 2: Run them and see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_surface_flow.py tests/test_reentry_model_spray.py -q`
Expected: failures for the missing names.

- [ ] **Step 3: Implement**

`shear_moment`: sort windward patches by `theta` (stable), prepend `(theta=0, integrand=0)`, cumulative trapezoid of `radius**2 * sin(theta) * tau` in `theta`, map back. `momentum_layer`: nan where `tau <= 0`, `r <= 0` or `I` is not positive. `evaluate` computes `shear_moment` from the final `tau` on the merged branch and passes `radius` as given. `SprayModel.slip_band = False` is a class attribute, and the constructor takes `slip_band=False`.

- [ ] **Step 4: Run them and see them pass**

Run: the same command. Expected: all pass.

- [ ] **Step 5: Commit, if Asha has asked**

### Task 2: The band in the skins, with the slurry on its switch

**Files:**
- Modify: `reentry_model/skinbody.py` (`girin_plan`, `slurry_modes`' inner `layer`, the sheared-depth history at about line 1080)
- Modify: `reentry_model/cli.py` (`--slip-band`, `spray_model_of`, the settings dict with `slurry_spray`)
- Test: `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_cli.py`

**Interfaces:**
- Consumes: `spray.spray_layer`, `SprayModel.slip_band`, `flow.shear_integral` (Task 1).
- Produces: `--slip-band on|off`; `settings["slip_band"]` in the results. `slurry_modes` gains the argument `theta`, and `radius` remains the body's transverse radius.

`girin_plan` replaces `melt_layer` with `spray_layer(flow, liq, model.slip_band, theta, radius)`. `b`, `regime_layer=F`, `rt_layer=L`, the mass run and `cap` stay as they are, so decision 3 needs no other edit and `on_slurry` stops holding a patch wherever `delta_m` is finite. In `slurry_modes`, `layer` calls `spray_layer` with the slurry's `LiquidProperties` and `sub`.

- [ ] **Step 1: Write the failing tests**

```python
def test_band_off_leaves_a_slipping_film_held(band_case):
    plan = girin_plan(*band_case, slurry=False)                # model.slip_band False
    assert plan.liquid_rate.sum() == 0.0

def test_band_on_sprays_the_liquid_and_reads_the_flowing_depth(band_case_on):
    plan = girin_plan(*band_case_on, slurry=False)
    assert plan.liquid_rate.sum() > 0.0
    assert_array_equal(plan.res.branch == BRANCH_THICK, plan.flowing > plan.res.delta_m)   # on wet windward columns

def test_band_release_is_capped_by_the_liquid_unless_the_slurry_is_on(band_case_on):
    assert_array_equal(girin_plan(*band_case_on, slurry=False).cap, band_case_on.skins.depth_above(T_LIQUIDUS))
    assert girin_plan(*band_case_on, slurry=True).slurry_rate.sum() > 0.0

def test_band_thin_branch_keeps_the_viscous_floor(band_case_on_thin):
    res = girin_plan(*band_case_on_thin, slurry=False).res
    assert res.gamma_tau[thin] >= 1.0 and (res.mdot[thin] == 0.0).all()

def test_band_above_knudsen_0p1_does_nothing(band_case_on_kn02):
    assert girin_plan(*band_case_on_kn02, slurry=False).liquid_rate.sum() == 0.0

def test_cli_band_is_off_by_default_and_recorded():
    assert parser.parse_args(RUN_ARGS).slip_band == "off"
```

The fixtures build a real `SkinSet` with 40 um of liquid over slurry, and a `SurfaceFlow.evaluate` at the stated Knudsen number; use the attribute names `SprayPlan` already has, and rename here to match.

- [ ] **Step 2: Run them and see them fail.** Run: `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -k band`. Expected: failures on the missing option.

- [ ] **Step 3: Implement** the edits named above.

- [ ] **Step 4: Run them and see them pass.** Then the off-switch check: the smoke run with and without this task's code has identical history and droplet files.

- [ ] **Step 5: Commit, if Asha has asked**

### Task 3: The band in the element model

**Files:**
- Modify: `reentry_model/body.py` (`_film_and_spray`, about lines 780-795)
- Test: `tests/test_reentry_model_melting.py`

**Interfaces:**
- Consumes: `spray.spray_layer`, `self.spray.slip_band`.

`delta_m` from `melt_layer` at `body.py:440` still drives the feed, deep runoff and transport. In `_film_and_spray` add a separate `spray_delta = spray_layer(flow, liq, self.spray.slip_band, self.theta, self.nose_radius())`, then a separate `v_spray` from `lubrication`, and pass those two to `self.spray.evaluate`. `held_mass` uses `spray_delta` too. `q` and the runoff keep the original depth, because decision 1's scope is spraying.

- [ ] **Step 1: Write the failing tests**

```python
def test_element_band_sprays_and_leaves_transport_alone(element_band_case):
    off, on = step(element_band_case, slip_band=False), step(element_band_case, slip_band=True)
    assert on.released > off.released
    assert on.feed == off.feed and on.runoff_flux == pytest.approx(off.runoff_flux)
```

- [ ] **Step 2: Run it and see it fail.** Run: `"$PY" -m pytest tests/test_reentry_model_melting.py -q -k band`.

- [ ] **Step 3: Implement** as above.

- [ ] **Step 4: Run it and see it pass.**

- [ ] **Step 5: Commit, if Asha has asked**

### Task 4: Measure, document and fold into sub-plan 22

**Files:**
- Modify: `reentry_model/spray.py` (module docstring: the band, its assumptions)
- Modify: `prototype/work-2026-10-28-skin-spraying-plan/part_22.md`, `make_plan22.py` (Task 5c)
- Create: `prototype/work-2026-10-28-skin-spraying/measure/task5c/` (the run script and notes)

- [ ] **Step 1: Run the full unit tier.** Run: `"$PY" -m pytest tests -q`. Expected: no failures against the baseline of 541 passed, 1 skipped, plus the new tests.

- [ ] **Step 2: Repeat the short 74 km case** (2.34 g handed to the skins), band on with the slurry off, then with it on. Compare with the earlier runs, as flown 0 g, hold lifted 0.73 g, slurry on 1.92 g. Record the mass sprayed, median radius and branch shares.

- [ ] **Step 3: Write Task 5c in `part_22.md`, add its files to `make_plan22.py`, regenerate, and confirm that every replay check passes.**

- [ ] **Step 4: Commit, if Asha has asked**

## Known limits (not tasks)

- The thick branch has no viscous gate, so with `--slurry-spray on` a slurry layer at a Reynolds number near 3 can spray by inviscid theory.
- `I` assumes an axisymmetric sphere; on an eroded body it is an approximation.
- `k_r` was calibrated with the Ranger-based depth, so droplet sizes in the band carry that.
- Records still label band patches as Couette closure, since `closure` is unchanged; `results` carries `slip_band`.
- The deposit folding that stopped the skin flights is not addressed.
