# Sub-plan: Task 7 — Spraying: instability branches, release bookkeeping, published reference values

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 3010–3949). Read `00-shared-context.md` first. **Amended 2026-10-02: the deep runoff and the per-patch conjugate depth** (section below). **Amended 2026-10-03: the molten cascade** (the section after it).

## Amendment of 2026-10-02 — what the spray sees of the deep liquid (no code change)

> Part of the deep-runoff amendment (sub-plan 09's amendment of this date; facts 46–53 in `00-shared-context.md`).

Task 9 now keeps a second liquid account per patch, `m_d`: the contiguous liquid below Girin's conjugate depth δ_m
that the deep runoff has moved there. Asha's three-zone rule of 2026-10-02 says that liquid is not sprayed; only the
skin within δ_m is. Nothing in `spray.py` changes, and these are the consequences, checked rather than assumed:

- **The deep liquid is never offered to a branch.** `SprayModel.evaluate` is called with the film alone (`b` and `m_f`),
  so every release — thick, thin, rarefied, regime 2 and the front-surface Rayleigh–Taylor mode — is drawn from the
  film. Deep liquid reaches the film only from the top, as far as the film is thinner than δ_m (one skin's worth per
  macro step), which Task 9 does before this step. `test_the_deep_liquid_is_never_sprayed_only_the_skin_is`
  (sub-plan 09) sprays a supercritical skin over a millimetre of deep liquid and finds the deep account unchanged to
  the bit.
- **The layer the regime test and the Rayleigh–Taylor criteria read is deeper where deep liquid sits.** `b_layer` is
  now the film, plus the contiguous molten material beneath it, plus the deep account's thickness. That changes no
  branch — a patch holding deep liquid is already deeper than δ_m — but the Rayleigh–Taylor depth criterion
  W cos φ h² ρ_l > 3Σ and the bounded form's tanh(kh) factor see the deeper pool, so the front-surface mode takes more
  of the nose: on the 100 mm flight to 120 s it released 72.5 g against 53.0 g without the deep runoff at the default
  0.5 s step (fact 49), 19.7 g against 17.2 g at 0.25 s and 11.4 g against 8.5 g at 0.125 s, where far less liquid
  lies below the conjugate depth (fact 50). The release itself falls sixfold between those steps with or without the
  deep runoff, because the molten layer the criterion reads is a time-step artefact of the melt step (fact 50).
- **Open, for Asha (fact 53).** The front-surface Rayleigh–Taylor mode is an instability of the whole liquid layer —
  its criterion is evaluated on the full depth — yet it releases only the film. Whether it should also release the
  deep liquid where it is applied is a modelling decision the three-zone rule does not settle: zone 2 is "not
  sprayed", but the rule is written about Girin's shear mechanism, and the front-surface mode is Girin & Kopyt's
  inertial one.

## Amendment of 2026-10-03 — what the spray sees of the molten cascade (no code change)

> Part of the molten-cascade amendment (sub-plan 09's amendment of this date; facts 54–61 in `00-shared-context.md`).

Task 9 now lets the surface recede through fully molten material within the step: an element a death exposes fully
molten (its mean temperature at or above the top of the feed ramp) is fed whole to the film, dies, and exposes the next.
Nothing in `spray.py` changes, and these are the consequences, checked rather than assumed:

- **The cascade's liquid is sprayed one step later.** The cascade runs in step (v), after the spray of step (iii), like
  the death remainders of fact 4, so the liquid it feeds is film when the next step's runoff and spray act on it. Every
  branch draws from that film as from any other; nothing is released twice and nothing skips the critical-Weber gate or
  the wave-fits test.
- **The layer the branch test and the Rayleigh–Taylor criteria read is what the conduction left molten in the step,
  not an accumulated backlog.** Measured on the 100 mm flight at the default step (fact 57): the molten layer's mean
  depth falls from 1.26 to 0.97 mm and the liquid held below the conjugate depth from a median 4.75 to 1.23 g, but the
  share of wet windward patches on the thick branch only from 12.4 % to 11.9 %, because a single molten wall-owning
  element already counts as 0.67–1.16 mm of liquid — more than the conjugate depth on every patch of the production
  mesh — so the test asks whether the owner's mean temperature is above T_feed when the spray runs (fact 58).
- **The release.** Unchanged in kind and nearly in amount at the default step: thin-branch share of the
  sprayed mass 10.5 % against 9.4 %, front-surface Rayleigh–Taylor release 51.9 against 53.0 g (56.3 against 72.5 g
  with the deep runoff on, because the deep liquid that deepened the layer its criterion sees is a tenth of what it
  was), 18.1 against 17.2 million droplets, median radius by number 180 µm either way. The sixfold fall of the
  Rayleigh–Taylor release between the default step and 0.125 s that the amendment of 2026-10-02 recorded is therefore
  not the backlog's doing; it follows the branch test's step dependence (fact 58), which this amendment does not
  change.

---

**Depends on:** Task 2 (material/liquid properties), Task 4 (dispersion table), and Task 5 (flow regime, wall pressure, surface flow).
**Produces, for later tasks:** the thick-film, thin-film, regime-2, and front-surface Rayleigh–Taylor spraying branches, the release bookkeeping, and the two Girin reference-value JSON files that Task 8 checks against.
**Character:** physics-heavy — the second-largest module in the plan.
**Read before implementing:** the master plan's own text notes that Girin's criterion "gates every branch," that the regime test is a two-stage test (depth, then viscosity) applied "in his order," and — importantly — that the "regime 2" branch is "inert on these flights and only the unit tests exercise it." Implement and test regime 2 fully, but do not assume it has been validated against a real physics run; only synthetic tests currently exercise it. Measured facts 30, 32, 33, and 36 in the shared context record two successive corrections to the Rayleigh–Taylor criterion — read them before touching that branch, since it was wrong in two different ways before the plan's current version.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan, and flag the regime-2 coverage gap explicitly in its test plan.

---

### Task 7: Spraying — instability branches, release bookkeeping, published reference values

**Files:**
- Create: `reentry_model/spray.py`, `data/reference_values/girin2017_table1.json`, `data/reference_values/girin1994_tables.json`
- Test: `tests/test_reentry_model_spray.py`

**Interfaces:**
- Consumes: `dispersion.DispersionTable` (Task 4), `surface_flow.BRANCH_FREE_MOLECULAR`, `CLOSURE_GIRIN` and `SurfaceFlowResult` fields (Task 5), `material.LiquidProperties` (Task 2).
- Produces: `melt_layer(flow, liquid) -> (delta_m, factor)` -- **nan/0 wherever the patch's closure is not `CLOSURE_GIRIN`**, since Eq. (2) and Ranger's Psi(phi) are continuum constructions and Task 5's gate decides whether they may be used at all; `rayleigh_taylor(deceleration, b, liquid) -> (active, lambda*, tau*)` (the unbounded fastest mode, reported only; `deceleration` is the component **normal to the film**, Girin & Kopyt's W sin(Theta) = W cos(phi) on this body, per patch or scalar); `rayleigh_taylor_bounded(deceleration, h, extent, liquid) -> (active, lambda, tau)`, the same mode in a pool of finite lateral extent, which admits only k >= pi/extent, with `extent` the contiguous **molten region** each patch belongs to (facts 30 and 33; reported only); `wave_fits(lam, extent) -> bool` -- one whole wavelength must fit inside that region, applied on every branch; `thin_film_mode(mach, momentum_flux, liquid) -> (lambda*, tau*)`; `SprayModel(liquid, k_r=0.17, k_t=1.1, we_critical=4.62, table=None).evaluate(flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, b_layer=None, extent=None, deceleration_n=None) -> SprayResult(branch, we_s, unstable, r, mdot, dm, dn, rt_active, delta_m, v_s, growth)` (`b_layer` is the contiguous liquid depth that decides the thick/thin branch, `extent` the molten region's lateral extent for the wave-fits test, `deceleration_n` the surface-normal deceleration for the reported Rayleigh-Taylor flag, and `growth` the selected mode's growth time per patch). Girin's criterion `We_s > We_cr` with We_s = rho_l V_s^2 min(delta_m, layer)/Sigma gates **every** branch, not just the thick one (fact 32); Girin's regime test then runs in **two stages**, in his order (fact 36). Stage one asks whether the *theoretical* melt boundary layer delta_m of his Eq. (2) can form in the liquid present: where it cannot, layer <= delta_m, the rigid core still stabilises the disturbances and the case is his regime 1, the Girin & Kopyt (1994) thin branch. Where it has formed, stage two decides the mechanism on the kinematic viscosities nu_melt = mu_l/rho_l against nu_gas = mu_e/rho_e -- `BRANCH_REGIME2` (regime 2, classical Kelvin-Helmholtz on the tangential-discontinuity profile, reusing Girin & Kopyt's lambda* and tau* with their rigid-wall factor cth(Lambda) -> 1 and the thick branch's torus rate rho_l pi r^2/(lambda* tau*), r = k_r lambda*, which does not reference the film depth) where the melt is the more viscous medium, and the thick branch (regime 3, his gradient instability) where the gas is. Liquid aluminium sits 3.3e4 to 6.9e6 from that threshold, so regime 2 is inert on these flights and only the unit tests exercise it; `source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state) -> list of 22-value rows`; `histogram(r, dn, dm) -> (dn per bin, dM per bin)`; `SprayModel(..., rt_spray=True)` applies the front-surface Rayleigh-Taylor mode on the patches where it passes both criteria and grows faster than the shear mode (fact 35; `rt_spray=False` reports it only, as every run before 2026-09-24 did). Constants `K_R, K_T, B_MIN, N_BINS = 40, R_MIN = 1e-6, R_MAX = 1e-2, BIN_EDGES, WE_BREAKUP = 12, CAPILLARY_TAU, SOURCE_COLUMNS, BRANCH_THICK/THIN/RAREFIED/RT = 0/1/2/3, BRANCH_REGIME2 = 4` (Girin thick / thin-film on the edge state / thin-film on the freestream / front-surface Rayleigh-Taylor / Girin & Kopyt's step profile for a deep film whose melt is the more viscous medium). The `closure` column replaces `regime` in `SOURCE_COLUMNS`. Spraying is never switched off by rarefaction: the gate selects the closure, and hence whether delta_m exists to define a thick film at all.

- [ ] **Step 1: Write the two reference-value files**

`data/reference_values/girin2017_table1.json`:

```json
{
  "_source": "Girin, O. G. (2017), A&A 606, A63, Table 1 and section 5 (R0 = 0.3 cm meteoroids); media properties from section 5; atmosphere density rho_a = 1e-7 g/cm3 (section 5, 'gas species density around the meteoroid'), rho_inf = rho_a/6 (the sixfold compression, section 5) reproduces the printed GI with We_inf on rho_inf and Re_inf on rho_a (measured 2026-09-20).",
  "constants": {
    "R0_m": 0.003,
    "rho_a_kgm3": 0.0001,
    "rho_inf_kgm3": 1.6666666666666667e-05,
    "mu_air_Pas": 6.8e-05,
    "k_r": 0.17,
    "k_t": 1.1,
    "we_cr": 4.62,
    "we_cr_theory": 3.08
  },
  "variants": {
    "I": {
      "description": "fast iron",
      "V_ms": 60000.0,
      "rho_m": 7800.0,
      "mu_m": 0.0058,
      "sigma": 1.2,
      "Re_inf": 529,
      "GI": 13.0,
      "m_g": 0.9,
      "t_sd_ms": 5.9,
      "sigma_ablation_s2cm2": 1.04e-12,
      "path_m": 350,
      "N": 1500000.0,
      "phi_cr_deg": 16.1,
      "t_f_us": 5.7,
      "z0_km": 105,
      "r_med_um": 26.9,
      "r_min_um": 4.0,
      "r_max_um": 135.0,
      "tau_end": 2.72,
      "t_ch_ms": 2.16
    },
    "II": {
      "description": "slow iron",
      "V_ms": 25000.0,
      "rho_m": 7800.0,
      "mu_m": 0.0058,
      "sigma": 1.2,
      "Re_inf": 222,
      "GI": 3.55,
      "m_g": 0.9,
      "t_sd_ms": 16.0,
      "sigma_ablation_s2cm2": 5.3e-12,
      "path_m": 400,
      "N": 370000.0,
      "phi_cr_deg": 31.4,
      "t_f_us": 31.0,
      "z0_km": 90,
      "r_med_um": 41.6,
      "r_min_um": 16.0,
      "r_max_um": 265.0,
      "tau_end": 3.11,
      "t_ch_ms": 5.15
    },
    "III": {
      "description": "stony",
      "V_ms": 60000.0,
      "rho_m": 3500.0,
      "mu_m": 0.174,
      "sigma": 0.36,
      "Re_inf": 529,
      "GI": 43.5,
      "m_g": 0.4,
      "t_sd_ms": 149.0,
      "sigma_ablation_s2cm2": 1.8e-14,
      "path_m": 9000,
      "N": 832,
      "phi_cr_deg": 8.8,
      "t_f_us": 194.0,
      "z0_km": 120,
      "r_med_um": null,
      "r_min_um": null,
      "r_max_um": null,
      "tau_end": 105.0,
      "t_ch_ms": 1.45
    }
  }
}
```


`data/reference_values/girin1994_tables.json`:

```json
{
  "_source": "Girin, A. G. & Kopyt, N. Kh. (1994), J. Aerosol Sci. 25(7), 1353-1357, Tables 1-2 and Eqs. (9)-(14). Table 1: side-surface KH mode, sigma_t = 1200 dyn/cm, rho_1 = 7.8 g/cm3, M = 3, r_d = lambda*/4, tau_d = tau*, mass rate m per unit area; the printed r_d imply an effective dynamic pressure 7.8 x rho_2 V0^2 (their shock-layer treatment: 'deceleration ~10 times ... acceleration to M = 2-3'), fitted on the first row and checked on the other five (measured 2026-09-20); m = rho_1 r_d / (2 tau_d) reproduces the printed rates. Table 2: front-surface Rayleigh-Taylor mode, Eq. (14); the printed lambda* are 10 x smaller than Eq. (14) gives (a units slip in the table), tau* matches to three digits.",
  "table1": {
    "sigma_Nm": 1.2,
    "rho1_kgm3": 7800.0,
    "mach": 3.0,
    "rho2_kgm3": [
      0.0001,
      0.001,
      0.01
    ],
    "V0_ms": [
      3000.0,
      10000.0
    ],
    "r_d_m": [
      [
        0.000193,
        1.73e-05
      ],
      [
        1.93e-05,
        1.73e-06
      ],
      [
        1.93e-06,
        1.73e-07
      ]
    ],
    "tau_d_s": [
      [
        0.00138,
        3.71e-05
      ],
      [
        4.34e-05,
        1.17e-06
      ],
      [
        1.38e-06,
        3.71e-08
      ]
    ],
    "m_kgm2s": [
      [
        547.0,
        1830.0
      ],
      [
        1740.0,
        5800.0
      ],
      [
        5470.0,
        18300.0
      ]
    ]
  },
  "table2": {
    "sigma_Nm": 1.2,
    "rho1_kgm3": 7800.0,
    "W_ms2": [
      1000.0,
      10000.0,
      100000.0
    ],
    "lambda_star_m_as_printed": [
      0.00042,
      0.000135,
      4.27e-05
    ],
    "lambda_star_m_eq14": [
      0.0042,
      0.00135,
      0.000427
    ],
    "tau_star_s": [
      0.001,
      0.00018,
      3.19e-05
    ]
  }
}
```


- [ ] **Step 2: Write the failing tests**

Create `tests/test_reentry_model_spray.py`:

```python
"""spray.py: delta_m, the thin-film and Rayleigh-Taylor modes against Girin & Kopyt (1994), the thick branch, release
bookkeeping, size caps and histograms."""
import json
import os

import numpy as np
import pytest

from reentry_model import dispersion, material, spray, surface_flow as sf
from helpers import REPO_ROOT

LIQ = material.LiquidProperties(2400.0, 1.3e-3, 0.86)


class Iron:
    rho, mu, sigma = 7800.0, 5.8e-3, 1.2


def test_thin_film_mode_reproduces_the_1994_table_1_scaling():
    """r_d = lambda*/4 and tau_d = tau* for the six (rho_2, V0) cases: one shock-layer factor on rho_2 V0^2 (fitted on
    the first entry, 7.8) reproduces every r_d within 0.5 % and tau_d within 1 %; the printed mass rates are
    rho_1 r_d / (2 tau_d) within 1 %."""
    t = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")))["table1"]
    rho2, V0 = np.array(t["rho2_kgm3"]), np.array(t["V0_ms"])
    lam, tau = spray.thin_film_mode(np.full((3, 2), t["mach"]), rho2[:, None] * V0[None, :] ** 2, Iron)
    r_tab, tau_tab, m_tab = (np.array(t[k]) for k in ("r_d_m", "tau_d_s", "m_kgm2s"))
    f_s = (lam / 4.0)[0, 0] / r_tab[0, 0]
    assert 7.5 < f_s < 8.1
    assert np.abs(lam / 4.0 / f_s / r_tab - 1.0).max() < 5e-3 and np.abs(tau / f_s ** 1.5 / tau_tab - 1.0).max() < 1e-2
    assert np.abs(Iron.rho * (lam / 4.0 / f_s) / (2.0 * tau / f_s ** 1.5) / m_tab - 1.0).max() < 1e-2
    assert spray.CAPILLARY_TAU == pytest.approx(0.798, rel=1e-3)                                  # tau* = 2 x 2 pi / omega_cap
    assert np.isinf(spray.thin_film_mode(np.array([3.0]), np.array([0.0]), Iron)[0][0])


def test_rayleigh_taylor_mode_reproduces_the_1994_table_2():
    t = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")))["table2"]
    for W, lam_ref, tau_ref in zip(t["W_ms2"], t["lambda_star_m_eq14"], t["tau_star_s"]):
        active, lam, tau = spray.rayleigh_taylor(W, np.array([1e-3, 1e-5]), Iron)
        assert lam[0] == pytest.approx(lam_ref, rel=2e-2) and tau[0] == pytest.approx(tau_ref, rel=1e-2)
        assert active[0] == (W * 1e-6 * Iron.rho > 3.0 * Iron.sigma) and not active[1]
    assert not spray.rayleigh_taylor(30.0, np.array([1e-4]), LIQ)[0][0]                          # our deceleration: inactive


class FakeFlow:
    """The fields spray.py reads: a shock-layer branch with the Girin closure unless told otherwise."""

    def __init__(self, n, closure=sf.CLOSURE_GIRIN, branch=sf.BRANCH_SHOCK_LAYER, u=2000.0, rho=7e-4, mu=1.6e-4,
                 delta_a=7e-3, mach=1.0, decel=30.0):
        self.closure = np.full(n, closure)
        self.branch = branch
        self.u_eff = self.u_e = np.full(n, u)
        self.rho_e, self.mu_e, self.delta_a, self.mach_e = np.full(n, rho), np.full(n, mu), np.full(n, delta_a), np.full(n, mach)
        self.deceleration = decel


class FakeState:
    class freestream:
        rho = 8e-5
    V = 7200.0
    ma = 24.0
    h = 71e3


def test_melt_layer_and_branches():
    flow = FakeFlow(3)
    delta_m, factor = spray.melt_layer(flow, LIQ)
    alpha, mu = 7e-4 / 2400.0, 1.6e-4 / 1.3e-3
    assert np.allclose(delta_m, (alpha / mu ** 2) ** (1 / 3) * 7e-3) and np.allclose(factor, (alpha * mu) ** (1 / 3) / (1 + (alpha * mu) ** (1 / 3)))
    assert 5e-5 < delta_m[0] < 2e-4                                                                 # ~0.1 mm (spec estimate)
    model = spray.SprayModel(LIQ)
    b = np.array([1e-3, 2e-5, 0.0])                                                                 # thick, thin, dry
    v_s = np.array([5.0, 20.0, 0.0])          # both films supercritical: Girin's criterion now gates every branch
    areas, m_f = np.full(3, 1e-5), b * LIQ.rho * 1e-5
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True, True]), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN, -1]
    we_s = LIQ.rho * 25.0 * delta_m[0] / LIQ.sigma
    assert res.we_s[0] == pytest.approx(we_s) and res.unstable[0] == (we_s > 4.62)
    d_f, im_f, _ = dispersion.DispersionTable()(we_s)
    lam = 2 * np.pi * delta_m[0] / d_f
    assert res.r[0] == pytest.approx(min(0.17 * lam, (3 * m_f[0] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    assert res.mdot[0] == pytest.approx(LIQ.rho * np.pi * (0.17 * lam) ** 2 / (lam * 1.1 * delta_m[0] / (5.0 * im_f)))
    lam_t, tau_t = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    assert res.r[1] == pytest.approx(min(lam_t[0] / 4.0, (3 * m_f[1] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    assert res.mdot[1] == pytest.approx(LIQ.rho * min(2e-5, lam_t[0] / 8.0) / tau_t[0])
    assert res.dm[1] == pytest.approx(min(res.mdot[1] * 1e-5 * 0.5, m_f[1])) and res.dm[2] == 0.0
    assert res.dn[1] == pytest.approx(res.dm[1] / (4 / 3 * np.pi * LIQ.rho * res.r[1] ** 3))
    assert not res.rt_active
    # Couette closure (the gate declined to certify the patch): no delta_m, so the thin mode with the edge state
    flow_c = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
    dm_c, factor_c = spray.melt_layer(flow_c, LIQ)
    assert np.all(np.isnan(dm_c)) and np.all(factor_c == 0.0)
    res = model.evaluate(flow_c, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]          # a thick film cannot take Girin's branch here
    # free-molecular branch: no edge state exists, the freestream drives the mode
    flow_fm = FakeFlow(2, closure=sf.CLOSURE_COUETTE, branch=sf.BRANCH_FREE_MOLECULAR)
    res = model.evaluate(flow_fm, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
    assert res.branch.tolist() == [spray.BRANCH_RAREFIED, spray.BRANCH_RAREFIED]
    lam_r, _ = spray.thin_film_mode(np.array([24.0]), np.array([8e-5 * 7200.0 ** 2]), LIQ)
    assert res.r[0] == pytest.approx(min(lam_r[0] / 4.0, (3 * (m_f[0] + 1e-9) / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    # leeward or dry: nothing
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([False, False, False]), 0.5, areas, m_f)
    assert not res.unstable.any() and res.dm.sum() == 0.0
    with pytest.raises(ValueError):
        spray.SprayModel(LIQ, k_r=0.0)


def test_source_rows_and_histogram():
    flow = FakeFlow(2)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.array([1e-3, 1e-3])
    m_f = b * LIQ.rho * 1e-5
    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.array([5.0, 5.0]), np.array([True, True]), 0.5, np.full(2, 1e-5), m_f)
    rows = spray.source_rows(10.0, 71e3, 7200.0, np.radians([30.0, 60.0]), np.zeros((2, 3)), np.tile([0.0, 1.0, 0.0], (2, 1)), flow, res, b, LIQ, FakeState())
    assert len(rows) == 2 and len(rows[0]) == len(spray.SOURCE_COLUMNS)
    row = dict(zip(spray.SOURCE_COLUMNS, rows[0]))
    assert row["time_s"] == 10.0 and row["altitude_km"] == 71.0 and row["theta_deg"] == pytest.approx(30.0) and row["branch"] == spray.BRANCH_THICK
    assert row["closure"] == sf.CLOSURE_GIRIN
    assert row["we_d"] == pytest.approx(8e-5 * 7200.0 ** 2 * 2 * row["r_m"] / LIQ.sigma) and row["oh"] == pytest.approx(LIQ.mu / np.sqrt(LIQ.rho * LIQ.sigma * 2 * row["r_m"]))
    assert row["breakup"] == float(row["we_d"] > 12.0) and row["dm_kg"] == res.dm[0] and row["v_s_ms"] == 5.0
    n_hist, m_hist = spray.histogram(np.array([2e-6, 5e-5, 5e-5, np.nan, 2e-2]), np.array([1.0, 2.0, 3.0, 4.0, 5.0]), np.array([1.0, 1.0, 1.0, 1.0, 1.0]))
    assert n_hist.sum() == 6.0 and m_hist.sum() == 3.0 and n_hist.shape == (spray.N_BINS,)          # nan and out-of-range dropped
    assert spray.BIN_EDGES[0] == 1e-6 and spray.BIN_EDGES[-1] == 1e-2


def test_girins_critical_weber_number_gates_the_thin_branch_too():
    """Girin & Kopyt's inviscid side mode is unstable at every wavelength and carries no threshold of its own, so the
    thin branch used to strip film wherever a wavelength existed -- at any angle, however small the shear. Girin's
    (2017) criterion We_s > We_cr now gates it as well, measured over the liquid the shear can actually reach,
    min(delta_m, layer): on a film shallower than the conjugate layer that is the film's own depth."""
    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)                # no delta_m: the thin branch
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.array([2e-5, 2e-5])
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    v_s = np.array([0.5, 20.0])                                   # subcritical, supercritical
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True]), 0.5, areas, m_f)
    we = LIQ.rho * v_s ** 2 * b / LIQ.sigma                       # layer = b here, and delta_m is nan
    assert res.we_s == pytest.approx(we) and we[0] < 4.62 < we[1]
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]          # classified either way
    assert not res.unstable[0] and res.dm[0] == 0.0 and np.isnan(res.r[0])        # but stable: nothing leaves
    assert res.unstable[1] and res.dm[1] > 0.0                                    # and the supercritical one sprays
    lam, tau = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    assert np.isfinite(lam[0]) and np.isfinite(tau[0])            # a wavelength existed for both -- only We decided


def test_a_wave_longer_than_its_melt_region_cannot_form():
    """A mode cannot develop in a puddle narrower than itself, and the domain is the contiguous molten region, not the
    mesh facet. The thin-branch wavelength diverges as the edge dynamic pressure falls, so near the stagnation point it
    is the fit that has to decide (2026-09-24)."""
    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, v_s = np.full(2, 2e-5), np.full(2, 20.0)                   # both supercritical
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    lam, _ = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    tight, loose = 0.5 * lam[0], 2.0 * lam[0]
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True]), 0.5, areas, m_f,
                         extent=np.array([tight, loose]))
    assert not res.unstable[0] and res.dm[0] == 0.0               # the wave is twice as wide as its melt
    assert res.unstable[1] and res.dm[1] > 0.0                    # here it fits
    assert spray.wave_fits(lam, None).all()                       # no extent given disables the test
    assert not spray.wave_fits(lam, np.zeros(1))[0]               # a region of no extent admits nothing


def test_the_rayleigh_taylor_driving_is_the_deceleration_normal_to_the_film():
    """Girin & Kopyt's dispersion relation carries the deceleration as W sin(Theta), Theta being the angle between the
    relative velocity of the media and the deceleration vector; on this body that is W cos(phi). It is the whole
    deceleration at the stagnation point and zero at the equator, where W lies in the surface and cannot lift the
    interface. Passing the whole deceleration instead over-drove every patch off the nose by 1/cos(phi)."""
    W, h = 95.0, np.full(4, 6e-3)          # deep enough that the flank still qualifies at cos 60 = 1/2:
    #                                       W cos(phi) h^2 rho_l > 3 Sigma needs h > 4.76 mm there, 3.37 mm at the nose
    phi = np.radians([0.0, 60.0, 90.0, 150.0])                    # nose, flank, equator, leeward
    w_n = W * np.maximum(np.cos(phi), 0.0)
    on, lam, tau = spray.rayleigh_taylor(w_n, h, LIQ)
    assert on[0] and on[1] and not on[2] and not on[3]            # driven at the nose, dead at and past the equator
    assert lam[0] == pytest.approx(2 * np.pi * np.sqrt(3 * LIQ.sigma / (W * LIQ.rho)))
    assert lam[1] == pytest.approx(lam[0] * np.sqrt(2.0))          # cos 60 = 1/2, and lambda* goes as W^-1/2
    assert lam[2] > 1e3                                           # cos 90 is round-off: kilometres of wavelength, no mode
    assert not np.isfinite(lam[3]) and not np.isfinite(tau[3])     # leeward, clamped to zero: no mode at all
    assert spray.rayleigh_taylor(W, h, LIQ)[0].tolist() == [True] * 4          # a scalar still broadcasts
    # and the bounded form takes the same normal component
    onb, lamb, taub = spray.rayleigh_taylor_bounded(w_n, h, np.full(4, 0.2), LIQ)
    assert onb[0] and not onb[2] and not onb[3] and taub[0] < taub[1] < np.inf


def test_the_front_surface_rayleigh_taylor_mode_is_applied_where_it_is_the_faster_mode():
    """Girin & Kopyt assign the front surface to their aperiodic Rayleigh-Taylor mode, not to the side mode, and with
    Girin's critical angle gating the shear branches the stagnation cap otherwise has no mechanism at all. The mode is
    applied where it passes *both* criteria -- the depth threshold W cos(phi) h^2 rho_l > 3 Sigma, which is what confines
    it to the nose, and an admissible wave in the molten region -- and where its growth time beats the shear mode's,
    since two modes of one interface cannot both break it up."""
    flow = FakeFlow(3, closure=sf.CLOSURE_COUETTE)            # no delta_m: the shear mode here is the weak thin one
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.full(3, 1e-4)
    layer = np.array([5e-3, 1e-3, 5e-3])                      # deep, too shallow for the threshold, deep
    v_s = np.array([0.1, 0.1, 40.0])                          # stable shear, stable shear, strong shear
    w_n = np.full(3, 95.0)                                    # at the stagnation point cos(phi) = 1
    areas, m_f = np.full(3, 1e-5), b * LIQ.rho * 1e-5
    deep, _, _ = spray.rayleigh_taylor(w_n, layer, LIQ)
    assert deep.tolist() == [True, False, True]               # h > lambda*/2 pi only for the deep pair
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f,
                         b_layer=layer, deceleration_n=w_n)
    assert res.branch.tolist() == [spray.BRANCH_RT, spray.BRANCH_THIN, spray.BRANCH_THIN]
    assert res.dm[0] > 0.0 and res.dm[1] == 0.0 and res.dm[2] > 0.0
    _, lam_rt, tau_rt = spray.rayleigh_taylor_bounded(w_n, layer, np.full(3, np.inf), LIQ)
    assert res.growth[0] == pytest.approx(tau_rt[0])          # the winning mode's growth time is the RT one
    assert res.growth[2] < res.growth[0]                      # and the shear mode is the faster one on patch 2
    assert res.mdot[0] == pytest.approx(LIQ.rho * min(b[0], lam_rt[0] / 2.0) / tau_rt[0])
    cap = (3 * m_f[0] / (4 * np.pi * LIQ.rho)) ** (1 / 3)
    assert lam_rt[0] > 1e-2 and res.r[0] == pytest.approx(cap)   # lambda* is centimetres; the film-mass cap binds
    # with the mode reported but not applied -- every run before 2026-09-24 -- the deep stable patch releases nothing
    off = spray.SprayModel(LIQ, rt_spray=False).evaluate(
        flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f, b_layer=layer, deceleration_n=w_n)
    assert off.branch.tolist() == [spray.BRANCH_THIN] * 3 and off.dm[0] == 0.0 and off.dm[2] > 0.0
    assert off.dm[2] == pytest.approx(res.dm[2])              # the shear-dominated patch is untouched either way


def test_the_applied_rayleigh_taylor_mode_needs_the_wave_to_fit_its_melt_region():
    """The lateral bound applies to this mode as it does to the shear ones: a pool narrower than half the marginal
    wavelength admits nothing at all, however deep it is."""
    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, layer, v_s = np.full(2, 1e-4), np.full(2, 5e-3), np.full(2, 0.1)
    w_n = np.full(2, 95.0)
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    lam_c = 2 * np.pi * np.sqrt(LIQ.sigma / (95.0 * LIQ.rho))          # the marginal wavelength
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f,
                         b_layer=layer, extent=np.array([0.4 * lam_c, 20.0 * lam_c]), deceleration_n=w_n)
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_RT]   # too narrow, then wide enough
    assert res.dm[0] == 0.0 and res.dm[1] > 0.0


def test_regime_2_gives_a_deep_film_with_a_viscous_melt_the_step_profile_mode():
    """Regime 2: where delta_m has formed and the melt is the more viscous medium, the patch takes Girin's
    Kelvin-Helmholtz cell -- the 1994 step-profile mode -- and not regime 3's gradient instability.

    delta_m is supplied directly rather than through `melt_layer`, and it has to be: Eq. (2) puts the predicted melt
    boundary layer at metres for a melt this viscous, which is Girin's own observation that for a high-viscosity melt the
    layer "becomes comparable with the radius of the meteoroid remnant". On a body this size such a melt therefore never
    passes stage one and lands in regime 1 instead, which the next test measures. Here the branch itself is exercised on
    the state stage two would hand it."""
    nu_gas = 1.6e-4 / 7e-4                                     # the edge state FakeFlow carries, 0.229 m2/s
    MELT_VISCOUS = material.LiquidProperties(2400.0, 10.0 * nu_gas * 2400.0, 0.86)   # nu_melt = 10 x nu_gas
    assert MELT_VISCOUS.mu / MELT_VISCOUS.rho > nu_gas
    flow = FakeFlow(2)
    model = spray.SprayModel(MELT_VISCOUS)
    delta_m = np.full(2, 1e-4)                                 # a conjugate depth the film can exceed
    b, v_s = np.full(2, 1e-3), np.full(2, 10.0)                # deep film; v_s large enough to clear We_cr
    areas, m_f = np.full(2, 1e-5), np.full(2, 1e-3 * 2400.0 * 1e-5)
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_REGIME2] * 2    # not BRANCH_THICK, and the RT mode does not take it
    assert res.we_s[0] > model.we_critical                     # Girin's gate still applies, as on every branch
    # the wavelength and growth time are Girin & Kopyt's Eqs. (11) and (12) -- the thin branch's pair, same mode
    lam = 1.5 * flow.mach_e[0] * MELT_VISCOUS.sigma / (flow.rho_e[0] * flow.u_eff[0] ** 2)
    tau = spray.CAPILLARY_TAU * lam ** 1.5 * np.sqrt(MELT_VISCOUS.rho / MELT_VISCOUS.sigma)
    assert res.growth[0] == pytest.approx(tau, rel=1e-12)
    assert res.r[0] == pytest.approx(spray.K_R * lam, rel=1e-12)            # r = k_r lambda*, below both caps here
    # the rate is Girin 2017's torus shedding on that wavelength, with no min(b, lambda*/8) factor: one torus of
    # cross-section r = k_r lambda* per wavelength per growth time, mdot = rho_l pi r^2 / (lambda* tau*)
    mdot = MELT_VISCOUS.rho * np.pi * (spray.K_R * lam) ** 2 / (lam * tau)
    assert res.mdot[0] == pytest.approx(mdot, rel=1e-12)
    # "not depth-limited" means the rate does not reference the film depth at all, as on the thick branch: doubling the
    # film leaves it unchanged, where the thin branch's min(b, lambda*/8) factor follows a film shallower than the wave
    deeper = model.evaluate(flow, FakeState(), 2.0 * b, delta_m, v_s, np.full(2, True), 0.5, areas, 2.0 * m_f)
    assert deeper.branch.tolist() == [spray.BRANCH_REGIME2] * 2
    assert deeper.mdot[0] == pytest.approx(res.mdot[0], rel=1e-12)
    # mdot is the instability's demand; the film is what limits the release, exactly as on every other branch
    assert mdot * areas[0] * 0.5 > m_f[0] and res.dm[0] == pytest.approx(m_f[0])


def test_a_kelvin_helmholtz_melt_fails_girins_first_stage_on_this_body():
    """A melt viscous enough for regime 2 cannot satisfy stage one here, so it routes to regime 1 -- correctly.

    delta_m/delta_a = (alpha/mu^2)^(1/3) = ((nu_melt/nu_gas)^2 rho_l/rho_e)^(1/3), so at the viscosity threshold itself
    the predicted melt boundary layer is already (rho_l/rho_e)^(1/3) ~ 151 air boundary-layer thicknesses, of order a
    metre at delta_a = 7 mm, and no film a 100 mm sphere can hold accommodates it. Stage one then says the layer never
    formed, which is Girin's dominant-ablation case and his Girin & Kopyt (1994) mode: the thin branch. That is what his
    own ordering produces and not a gap in it -- recorded so the regime-2 branch is not later mistaken for dead code that
    ought to have fired."""
    flow = FakeFlow(1)
    nu_gas = flow.mu_e[0] / flow.rho_e[0]
    for ratio in (1.0, 10.0, 100.0):
        liq = material.LiquidProperties(2400.0, ratio * nu_gas * 2400.0, 0.86)
        delta_m, _ = spray.melt_layer(flow, liq)
        expected = (ratio ** 2 * liq.rho / flow.rho_e[0]) ** (1.0 / 3.0)            # delta_m / delta_a
        assert delta_m[0] / flow.delta_a[0] == pytest.approx(expected, rel=1e-12)
        assert delta_m[0] > 1.0                                  # metres of liquid would be needed to out-deepen it
    # at the threshold itself the factor is (rho_l/rho_e)^(1/3), and no film the body can hold reaches it
    liq = material.LiquidProperties(2400.0, nu_gas * 2400.0, 0.86)
    delta_m, _ = spray.melt_layer(flow, liq)
    assert delta_m[0] / flow.delta_a[0] == pytest.approx((liq.rho / flow.rho_e[0]) ** (1.0 / 3.0), rel=1e-12)
    b = np.full(1, 1e-2)                                         # 10 mm of melt: a tenth of the sphere's radius
    res = spray.SprayModel(liq).evaluate(flow, FakeState(), b, delta_m, np.full(1, 10.0), np.full(1, True), 0.5,
                                        np.full(1, 1e-5), b * liq.rho * 1e-5)
    assert res.branch.tolist() == [spray.BRANCH_THIN]            # regime 1, the stabilised-core case, as Girin has it


def test_regime_2_never_fires_for_liquid_aluminium():
    """The regression guard: for the real melt the branch is inert, so every earlier run is reproduced exactly.

    Liquid aluminium has nu_melt = 5.42e-7 m2/s against a post-shock edge nu_gas measured at 0.018-3.75 m2/s over the
    two physics flights, a ratio of 3.3e4 to 6.9e6, so the melt is far the less viscous medium everywhere and a formed
    layer is always regime 3 -- Girin's gradient instability, the thick branch -- exactly as before 2026-09-25 (0 of
    533 960 windward wet patch-steps on the 100 mm flight and 0 of 42 268 on the 50 mm one were assigned regime 2)."""
    flow = FakeFlow(2)
    assert LIQ.mu / LIQ.rho < flow.mu_e[0] / flow.rho_e[0]       # aluminium is far the less viscous medium
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, v_s = np.full(2, 1e-3), np.full(2, 10.0)
    assert (b > delta_m).all()                                   # a genuinely deep film by Girin's own conjugate depth
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_THICK] * 2
```


- [ ] **Step 3: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_spray.py -q`
Expected: ImportError (`reentry_model.spray`).

- [ ] **Step 4: Create `reentry_model/spray.py`**

```python
"""Melt spraying: Girin's gradient instability per patch and the droplet release bookkeeping (spec Step 3 section 9).

Branches per windward patch with film (regime from surface_flow, b from the film, delta_m from Girin's Eq. 2):
  thick, continuum/slip (b > delta_m): delta_m = (alpha/mu^2)^(1/3) delta_a, V_s = (alpha mu)^(1/3)/(1 + (alpha
      mu)^(1/3)) u_eff, We_s = rho_l V_s^2 delta_m / Sigma; unstable when We_s > We_cr; lambda_f = 2 pi delta_m /
      Delta_f(We_s), r = k_r lambda_f, t_per = k_t delta_m / (V_s Im Omega_f(We_s)), stripping rate per area
      mdot = rho_l pi r^2 / (lambda_f t_per) (one torus of cross-section r per wavelength per period, Girin 2017
      Eqs. 5-8 and Appendix D).
  thin, continuum/slip (b <= delta_m): Girin & Kopyt 1994 side mode with the edge state: lambda* = 1.5 M_e Sigma /
      (rho_e u_eff^2) (their 1.5 M d / We_d, the film thickness cancels), r = lambda*/4, growth time tau* =
      2 tau_v(lambda*) = 0.798 lambda*^1.5 (rho_l/Sigma)^1/2 (their Eq. 12: twice the capillary-wave period),
      release rate mdot = rho_l min(b, lambda*/8) / tau* (their Table 1's mass rate is rho_1 r_d / (2 tau_d) =
      rho_1 lambda*/(8 tau*), reproduced 2026-09-20; the film supplies at most its thickness). Active when b >= B_MIN.
      (Their dissipation cut-off lambda_t = lambda*/3 is always below lambda*, so it never limits the mode.)
  regime 2, continuum/slip (liquid layer > delta_m, as for the thick branch, *and* nu_melt > nu_gas): the
      Kelvin-Helmholtz corner of Girin's own two-regime classification -- a deep film whose melt is the more viscous
      medium, where the interfacial profile is a tangential discontinuity rather than a resolved gradient. Wavelength and
      growth time are then Girin & Kopyt's pair, lambda* = 1.5 M_e Sigma/(rho_e u_eff^2) and
      tau* = 0.798 lambda*^1.5 (rho_l/Sigma)^1/2 (their Eqs. 11 and 12), with the rigid-wall factor cth(Lambda) of their
      Eq. (5) -> 1 because the film is deeper than the wave sees, so the unmelted body no longer stiffens the response.
      The release is the torus form Girin (2017) uses on the thick branch, mdot = rho_l pi r^2/(lambda* tau*) with
      r = k_r lambda*, and *not* the thin branch's rho_l min(b, lambda*/8)/tau*: that factor is there because a film
      shallower than the wave can give up no more than its own depth, which is not this film's position. Like the thick
      branch's rate it therefore does not reference b at all -- and it is 0.73 x the thin branch's rate at the same
      wavelength, not a larger one, wherever that branch is not itself depth-limited; what bounds the release here is the
      film mass, through dm = min(mdot A dt, m_f) as on every branch.
      PROVENANCE: both tests are Girin's, and the numbering is his. Girin (2017) Sect. 1 sets out "two regimes of
      meteoroid ablation due to instability". In *dominant ablation* the mass loss overtakes fusion, the molten layer is
      thinner than the melt velocity boundary layer, b < delta_m, the rigid core still stabilises the disturbances, and
      "the results obtained in Girin & Kopyt (1994) are thus valid" -- this module's thin branch. In *dominant fusion*,
      b > delta_m, the profile is the completed boundary layer: free conjugated layers, delta_a in air and delta_m in the
      melt, "independent from the solid core". His Sect. 2 then selects the mechanism inside that conjugated picture on
      the kinematic viscosities. The "classical Kelvin-Helmholtz type ... is in action only when the liquid kinematic
      viscosity is greater than that one of the gas: nu_m > nu_a", because then "the liquid boundary layer is much
      thicker than the gas one" and V_s << V_a, so "the velocity profile is close to the discontinuous one (tangential
      discontinuity) which is the base for the Kelvin-Helmholtz mechanism". For the inverse inequality nu_a > nu_m the
      thicknesses become comparable, delta_m = O(delta_a) and V_s = O(V_a), the profile is "inflated", and the mechanism
      is "completely different from the Kelvin-Helmholtz one", being set by his Eq. (1) -- the gradient instability the
      thick branch solves through `dispersion`. So regime 2 is his Kelvin-Helmholtz cell and the thick branch is his
      gradient cell of the same two-by-two; what this model supplies is only the arithmetic for the cell he does not
      develop, by reusing the 1994 mode that belongs to that profile.
      UNREACHABLE FOR THIS PROBLEM, measured 2026-09-25, and consistent with Girin's own description rather than in
      tension with it. (a) Liquid aluminium is nowhere near the Kelvin-Helmholtz cell: nu_melt = 5.42e-7 m2/s against an
      edge nu_gas measured at 0.0178-2.04 m2/s over the 100 mm flight to 110 s and 0.0951-3.75 m2/s over the whole 50 mm
      flight, so nu_gas/nu_melt runs 3.3e4 to 6.9e6 and never approaches one. This body sits firmly in Girin's
      nu_a > nu_m case, which is exactly why his gradient instability is the right model for it. Directly: 0 of 533 960
      windward wet patch-steps (100 mm) and 0 of 42 268 (50 mm) take this branch. (b) The two tests reinforce each other
      rather than crossing, because Eq. (2) carries the same ratio:
      delta_m/delta_a = (alpha/mu^2)^(1/3) = ((nu_melt/nu_gas)^2 rho_l/rho_e)^(1/3), so at nu_melt = nu_gas the conjugate
      depth is already (rho_l/rho_e)^(1/3) ~ 151 air boundary-layer thicknesses, about a metre at delta_a = 7 mm -- which
      is Girin's "much thicker than the gas one" restated numerically. A melt viscous enough to be Kelvin-Helmholtz
      therefore carries a delta_m no film here can exceed, so the gate puts it at b < delta_m, which is Girin's
      dominant-ablation case and his Girin & Kopyt (1994) mode: the thin branch. That is the correct routing and not a
      gap. Reaching regime 2 needs a viscous melt *and* a melt layer deeper than a metre-scale delta_m, which this
      problem never produces; the branch is kept because it is where that cell belongs if the material or the altitude
      band changes, and it is exercised by test with delta_m supplied directly.
      Where it would show up, checked rather than assumed: the source table's `branch` column (`SOURCE_COLUMNS` index 8)
      and nowhere else. The surface VTK carries surface_flow's `closure`, not the spray branch, so a regime-2 patch is
      not distinguishable there; no history column is added; and `thick_branch_fraction` is `film.lubrication`'s own
      thick flag -- isfinite(delta_m) & (layer > delta_m), the whole deep set -- so it would count a regime-2 patch as
      thick. That column measures how deep the liquid is, not which instability took the patch.
Two gates apply to every branch (2026-09-24):
  * Girin's (2017) stability criterion We_s > We_cr, with We_s = rho_l V_s^2 min(delta_m, layer) / Sigma -- his
    closed-form critical angle phi_cr evaluated against the model's own local flow, which holds on an eroded body where
    a single angle does not. On the thick branch this is his We_s = rho_l V_s^2 delta_m / Sigma unchanged; on a film
    shallower than the conjugate layer it is measured over the film's own depth. The thin branch had no threshold
    before: Girin & Kopyt's inviscid side mode has none of its own (it is unstable at every wavelength, their
    "unlimiting potential intensity of high-frequency disturbances"), so it stripped film at any angle however small
    the shear -- 2.2 g inside 10 deg on the 100 mm flight, from a wave 6 x wider than its own facet.
  * `wave_fits`: one whole wavelength must fit inside the *contiguous molten region* the patch belongs to, not on the
    patch (a mesh facet is not a physical boundary and the melt is continuous across it).
KNOWN LIMITATION, deferred: Girin & Kopyt derive the side mode for the side surface (their sin Theta -> 0,
V_0 -> V_inf) and assign the *front* surface (V_0 -> 0, sin Theta -> 1) to a different, aperiodic Rayleigh-Taylor
solution, their Eq. (14). The branch selection here has no angular dependence -- it turns only on layer vs delta_m --
so the side-surface asymptotics are applied from the stagnation point to the equator, including where the paper's own
premise fails. Implementing their front/side split (sin Theta = cos phi on this body) is future work; until then the
near-nose release rests on a solution used outside its stated domain, which is 0.2 % of the sprayed mass.
  free-molecular branch, thin film: the thin mode with the freestream momentum flux, lambda* = 1.5 M_inf Sigma /
      (rho_inf V^2), and the Couette film velocity from tau_fm -- an extrapolation of a continuum film theory
      (spec section 17.3), sensitive to --rarefied-shear.
The branch selection follows the closure gate of surface_flow, not a regime of its own: Girin's dispersion relation
contains no gas parameters (it is a melt-side instability), so spraying is never switched off by rarefaction -- what
the gate decides is which closure supplies the melt velocity scale, and hence whether delta_m exists to define a
"thick" film at all (spec section 18).
The 1994 front-surface Rayleigh-Taylor criterion W b^2 rho_l > 3 Sigma is evaluated and reported (lambda* =
2 pi (3 Sigma/(W rho_l))^1/2, tau* = (27 Sigma/(4 W^3 rho_l))^1/4), never applied. W is the deceleration *normal to the
film* -- their W sin(Theta), which on this body is W cos(phi): the whole deceleration at the stagnation point, zero at
the equator where it lies in the surface (corrected 2026-09-24; the whole deceleration over-drove every patch off the
nose by 1/cos(phi)).

Release per patch and step: dm = min(mdot A dt, m_f) (the film after feed and runoff), dn = dm / (4/3 pi rho_l r^3);
one source row per emitting patch-step; histograms dn(r), dM(r) on N_BINS log bins over [R_MIN, R_MAX]."""
from dataclasses import dataclass, field

import numpy as np

from . import dispersion
from .surface_flow import BRANCH_FREE_MOLECULAR, CLOSURE_GIRIN

K_R, K_T = 0.17, 1.1
B_MIN = 1.0e-6                     # m, films thinner than this do not spray (numerical floor)
N_BINS, R_MIN, R_MAX = 40, 1.0e-6, 1.0e-2
BIN_EDGES = np.logspace(np.log10(R_MIN), np.log10(R_MAX), N_BINS + 1)
WE_BREAKUP = 12.0                  # Pilch-Erdman: secondary breakup expected above this droplet Weber number
CAPILLARY_TAU = 4.0 * np.pi / (2.0 * np.pi) ** 1.5      # 0.798: tau* = 2 x 2 pi / omega_cap(lambda*)
SOURCE_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "theta_deg", "x_m", "y_m", "z_m", "closure", "branch", "b_m",
                  "delta_m_m", "we_s", "r_m", "dn", "dm_kg", "v_s_ms", "tx", "ty", "tz", "we_d", "oh", "breakup"]
BRANCH_THICK, BRANCH_THIN, BRANCH_RAREFIED, BRANCH_RT = 0, 1, 2, 3      # Girin thick / thin-film with the edge state /
#                                                                        thin-film with the freestream / front-surface Rayleigh-Taylor
BRANCH_REGIME2 = 4                 # deep film whose melt is the more viscous medium: Girin & Kopyt's step profile again


def melt_layer(flow, liquid):
    """delta_m per patch from Girin's Eq. (2) and the shear velocity factor (alpha mu)^(1/3)/(1 + (alpha mu)^(1/3)).

    Both are nan/0 wherever the patch's closure is not CLOSURE_GIRIN: Eq. (2) is a conjugate-boundary-layer result and
    Ranger's Psi(phi) for delta_a is a continuum construction, so the wall Knudsen gate of surface_flow decides whether
    they may be used at all (spec section 18). Where they may not, the film takes the Couette closure
    V_s = tau_w b / mu_melt (film.lubrication's thin branch) and the thin-film spray mode."""
    girin = (flow.closure == CLOSURE_GIRIN) & np.isfinite(flow.delta_a)
    alpha = flow.rho_e / liquid.rho
    mu = flow.mu_e / liquid.mu
    with np.errstate(divide="ignore", invalid="ignore"):
        delta_m = np.where(girin & (mu > 0.0), (alpha / np.where(mu > 0.0, mu, 1.0) ** 2) ** (1.0 / 3.0) * flow.delta_a, np.nan)
        am = (alpha * mu) ** (1.0 / 3.0)
        factor = np.where(girin & np.isfinite(am), am / (1.0 + am), 0.0)
    return delta_m, factor


def rayleigh_taylor_bounded(deceleration, h, extent, liquid):
    """The front-surface Rayleigh-Taylor mode in a pool of *finite* lateral extent, which is what a melt pool is.

    `deceleration` is the surface-normal component W cos(phi), as in `rayleigh_taylor`; `extent` is the lateral extent
    of the contiguous molten region the patch belongs to, not of the patch (the mesh facet is not a physical boundary).

    For heavy liquid over light gas the dispersion relation with surface tension is n^2 = (W k - sigma k^3/rho) tanh(k h),
    unstable for k below the cutoff k_c = sqrt(W rho/sigma) -- that is, for wavelengths *longer* than
    lambda_c = 2 pi sqrt(sigma/(W rho)) -- and fastest at k_max = k_c/sqrt(3), which is the lambda* the unbounded form
    reports. A pool of lateral extent L with its rim pinning the interface admits only k >= k_1 = pi/L (half a
    wavelength across the pool), so:
      * L < lambda_c/2  : no admissible wavelength is unstable -- the pool is stable however deep it is;
      * L < lambda*/2   : unstable, but the fastest mode does not fit and the growth rate is the one at k_1, not 1/tau*;
      * otherwise       : the fastest mode fits and the unbounded growth rate applies.
    The tanh(k h) factor is the finite-depth correction and is not negligible here: at h = 2 mm and W = 30 m/s^2,
    k h = 0.58 and tanh = 0.52, so a shallow pool grows at about half the deep-layer rate.

    Comparing lambda* against the nose diameter, as an earlier diagnostic did, is wrong twice over: the mode must fit
    in the *unstable region*, not on the nose, and the yes/no question is decided by the marginal wavelength lambda_c,
    not by the fastest-growing one (measured 2026-09-23). Returns (active, wavelength [m], growth time [s])."""
    W = np.maximum(np.asarray(deceleration, dtype=float), 0.0)    # the surface-normal component, W cos(phi)
    h = np.asarray(h, dtype=float)
    L = np.asarray(extent, dtype=float)
    shape = np.broadcast_shapes(W.shape, h.shape, L.shape)
    live = W > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        k_c = np.where(live, np.sqrt(np.where(live, W, 0.0) * liquid.rho / liquid.sigma), 0.0)
        k_max = k_c / np.sqrt(3.0)
        k_1 = np.where(L > 0.0, np.pi / np.where(L > 0.0, L, 1.0), np.inf)
    k = np.maximum(k_1, k_max)                              # the fastest wavelength the pool actually admits
    active = live & (k < k_c) & (h > 0.0) & np.isfinite(k)
    with np.errstate(invalid="ignore", divide="ignore"):
        n2 = (W * k - liquid.sigma * k ** 3 / liquid.rho) * np.tanh(np.clip(k * h, 0.0, 30.0))
        tau = np.where(active & (n2 > 0.0), 1.0 / np.sqrt(np.where(n2 > 0.0, n2, 1.0)), np.inf)
        lam = np.where(active, 2.0 * np.pi / k, np.inf)
    on = np.broadcast_to(active & (n2 > 0.0), shape).copy()
    return on, np.broadcast_to(lam, shape).copy(), np.broadcast_to(tau, shape).copy()


def rayleigh_taylor(deceleration, b, liquid):
    """(active, lambda*, tau*) of the 1994 front-surface RT mode for the film thickness b.

    `deceleration` is the component of the body deceleration **normal to the film**, which is what drives the mode:
    Girin & Kopyt write it W sin(Theta) with Theta the angle between the relative velocity of the media and the
    deceleration vector, and on a sphere that is W cos(phi) -- full at the stagnation point, zero at the equator where
    the deceleration lies in the surface and cannot lift the interface at all. Passing the whole deceleration instead
    over-drives every patch off the nose by 1/cos(phi) (corrected 2026-09-24). A scalar still broadcasts."""
    W = np.maximum(np.asarray(deceleration, dtype=float), 0.0)
    b = np.asarray(b, dtype=float)
    shape = np.broadcast_shapes(W.shape, b.shape)
    active = np.broadcast_to(W * b * b * liquid.rho > 3.0 * liquid.sigma, shape).copy()
    live = W > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(live, 2.0 * np.pi * np.sqrt(3.0 * liquid.sigma / np.where(live, W * liquid.rho, 1.0)), np.inf)
        tau = np.where(live, (27.0 * liquid.sigma / (4.0 * np.where(live, W, 1.0) ** 3 * liquid.rho)) ** 0.25, np.inf)
    return active, np.broadcast_to(lam, shape).copy(), np.broadcast_to(tau, shape).copy()


def thin_film_mode(mach, momentum_flux, liquid):
    """(lambda*, tau*) of the Girin & Kopyt 1994 side mode for gas Mach number and rho V^2 [Pa]."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(momentum_flux > 0.0, 1.5 * mach * liquid.sigma / momentum_flux, np.inf)
    tau = CAPILLARY_TAU * lam ** 1.5 * np.sqrt(liquid.rho / liquid.sigma)
    return lam, tau


def wave_fits(lam, extent):
    """True where one whole wavelength fits inside the contiguous molten region the patch belongs to.

    A mode cannot form in a puddle narrower than itself. The region, not the patch, is the domain: a mesh facet is a
    bookkeeping boundary, not a physical one, and the melt is continuous across it. `extent` of None disables the test;
    a region of zero extent admits nothing. The Rayleigh-Taylor form of the same argument is in
    `rayleigh_taylor_bounded`, where the rim pins the interface and half a wavelength is enough; a travelling shear
    wave has to complete a period inside the melt, so one whole wavelength is required here (decided 2026-09-24)."""
    if extent is None:
        return np.ones(np.shape(lam), dtype=bool)
    L = np.asarray(extent, dtype=float)
    with np.errstate(invalid="ignore"):
        return (L > 0.0) & np.isfinite(lam) & (lam <= L)


@dataclass
class SprayResult:
    branch: np.ndarray          # 0 thick, 1 thin, 2 rarefied, -1 none
    we_s: np.ndarray
    unstable: np.ndarray        # bool
    r: np.ndarray               # m, droplet radius (nan where none)
    mdot: np.ndarray            # kg/(m2 s) stripping rate
    dm: np.ndarray              # kg released this step per patch
    dn: np.ndarray              # droplets released this step per patch
    rt_active: bool
    delta_m: np.ndarray
    v_s: np.ndarray
    growth: np.ndarray = None   # s, the selected mode's growth time per patch (nan where no mode): what decides
    #                             whether the shear instability outruns the Rayleigh-Taylor mode on the same liquid


class SprayModel:
    def __init__(self, liquid, k_r=K_R, k_t=K_T, we_critical=dispersion.WE_CRITICAL_PRACTICAL, table=None,
                 rt_spray=True):
        if k_r <= 0.0 or k_t <= 0.0 or we_critical <= 0.0:
            raise ValueError("k_r, k_t and the critical Weber number must be > 0")
        self.liquid, self.k_r, self.k_t, self.we_critical = liquid, k_r, k_t, we_critical
        self.rt_spray = bool(rt_spray)               # apply the front-surface mode, or only report it
        self.table = table or dispersion.DispersionTable()

    def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, b_layer=None,
                 extent=None, deceleration_n=None):
        """Per-patch instability and release for the film thickness b, melt-layer thickness delta_m and film surface
        velocity v_s (from film.lubrication), over the step dt; m_f is the film mass available, radius the body's.

        `b_layer`, if given, is the depth of liquid beneath the wall (film plus the contiguous molten material under
        it): the thick/thin test belongs on that, since it asks whether the gas shear reaches the bottom of the liquid,
        while the release rates and the droplet cap stay on the film, which is the mass that can actually leave.
        `extent`, if given, is the lateral extent of the contiguous molten region each patch belongs to, for the
        wave-fits test; `deceleration_n` is the surface-normal deceleration W cos(phi) for the reported
        Rayleigh-Taylor flag, defaulting to the whole deceleration when it is not supplied."""
        liq = self.liquid
        n = b.size
        layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
        branch = np.full(n, -1)
        growth = np.full(n, np.nan)          # the mode's growth time per patch [s], for comparison with other modes
        r = np.full(n, np.nan)
        mdot = np.zeros(n)
        has_film = windward & (b >= B_MIN)
        # Girin's regime test, in his order, two stages (2026-09-25). Stage one: delta_m is a *theoretical* thickness --
        # the melt velocity boundary layer Eq. (2) predicts -- so the question the liquid depth answers is whether that
        # layer can physically form in the melt that is present. If the liquid cannot accommodate it, layer <= delta_m,
        # the layer did not form, the rigid core still stabilises the disturbances and the case is Girin's dominant
        # ablation: regime 1, his Girin & Kopyt (1994) mode, the thin branch below. If it did form, layer > delta_m, the
        # profile is the conjugated pair (delta_a in air, delta_m in melt) and the case is dominant fusion.
        deep = has_film & np.isfinite(delta_m) & (layer > delta_m)         # delta_m fits in the liquid: it formed
        # Stage two, only where it formed: which mechanism the conjugated profile carries, decided on the kinematic
        # viscosities. nu_melt > nu_gas puts the thicker viscous layer in the melt, V_s << V_a and the profile close to a
        # tangential discontinuity -- classical Kelvin-Helmholtz, regime 2. nu_gas > nu_melt makes the two thicknesses
        # comparable and the profile inflated, which is the gradient instability of his Eq. (1) -- regime 3, the thick
        # branch. Both stages are Girin (2017), Sects. 1 and 2; see the module docstring for the quotations and for the
        # measurement that puts liquid aluminium firmly in regime 3.
        nu_melt = liq.mu / liq.rho
        with np.errstate(divide="ignore", invalid="ignore"):
            nu_gas = np.where(flow.rho_e > 0.0, flow.mu_e / np.where(flow.rho_e > 0.0, flow.rho_e, 1.0), np.inf)
        step_profile = np.asarray(nu_melt > nu_gas)                        # the melt carries the thicker viscous layer
        regime2 = deep & step_profile                                     # Kelvin-Helmholtz
        thick = deep & ~step_profile                                      # regime 3: Girin 2017's gradient instability
        free_molecular = flow.branch == BRANCH_FREE_MOLECULAR
        rarefied = has_film & ~deep & free_molecular                      # no edge state exists: the freestream drives the mode
        thin = has_film & ~deep & ~rarefied
        # Girin's (2017) stability criterion, now on *every* branch. His phi_cr comes from We_s > We_cr with the sphere's
        # own edge solution substituted in; evaluated locally instead of through that closed form it holds on an eroded
        # body too, which a single critical angle does not. The shear acts over the liquid it can reach, min(delta_m,
        # layer), so the thick branch is his We_s = rho_l V_s^2 delta_m / Sigma unchanged while a film shallower than the
        # conjugate layer is measured on its own depth. Below We_cr the surface is stable and nothing leaves. The thin
        # branch had no stability threshold before this (Girin & Kopyt's inviscid side mode has none of its own, being
        # unstable at every wavelength), so it stripped film at any angle however small the shear (added 2026-09-24).
        L = None if extent is None else np.asarray(extent, dtype=float)
        with np.errstate(invalid="ignore"):
            d_shear = np.where(np.isfinite(delta_m), np.minimum(delta_m, layer), layer)
            we_s = np.where(has_film, liq.rho * v_s ** 2 * d_shear / liq.sigma, 0.0)
        supercritical = we_s > self.we_critical
        # thick: Girin 2017
        if thick.any():
            vs, dm_ = v_s[thick], delta_m[thick]
            delta_f, im_f, _ = self.table(we_s[thick])
            with np.errstate(divide="ignore", invalid="ignore"):
                lam_f = np.where(delta_f > 0.0, 2.0 * np.pi * dm_ / np.where(delta_f > 0.0, delta_f, 1.0), np.inf)
            unstable = (supercritical[thick] & (im_f > 0.0) & (vs > 0.0)
                        & wave_fits(lam_f, None if L is None else L[thick]))
            lam = np.where(unstable, lam_f, np.nan)
            rr = self.k_r * lam
            with np.errstate(divide="ignore", invalid="ignore"):
                t_per = np.where(unstable, self.k_t * dm_ / (vs * np.where(unstable, im_f, 1.0)), np.inf)
                rate = np.where(unstable, liq.rho * np.pi * rr * rr / (lam * t_per), 0.0)
            r[thick], mdot[thick], branch[thick] = rr, rate, BRANCH_THICK
            growth[thick] = t_per
        # regime 2: the other half of the deep-film set -- Girin & Kopyt's step profile, in their deep-film limit.
        # Wavelength and growth time are their Eqs. (11) and (12), the same pair the thin branch uses, because the mode is
        # the same one; what differs is the supply. The thin branch multiplies by min(b, lambda*/8) because a film
        # shallower than the wave can give up no more than its own depth, and here the film is deeper than the sheared
        # layer, so that factor is dropped and the rate is Girin's (2017) torus shedding on this wavelength,
        # mdot = rho_l pi r^2/(lambda tau*) with r = k_r lambda -- one torus of cross-section r per wavelength per growth
        # time. mdot stays the instability's *demand*; what the film can actually supply is imposed once, below, by
        # dm = min(mdot A dt, m_f), exactly as on every other branch.
        if regime2.any():
            m = regime2
            lam, tau = thin_film_mode(flow.mach_e[m], (flow.rho_e * flow.u_eff ** 2)[m], liq)
            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & supercritical[m]
                  & wave_fits(lam, None if L is None else L[m]))
            rr = np.where(ok, self.k_r * lam, np.nan)
            with np.errstate(divide="ignore", invalid="ignore"):
                rate = np.where(ok, liq.rho * np.pi * self.k_r ** 2 * lam / np.where(ok, tau, 1.0), 0.0)
            r[m], mdot[m], branch[m] = rr, rate, BRANCH_REGIME2
            growth[m] = np.where(ok, tau, np.nan)
        # thin: Girin & Kopyt 1994 with the edge state
        for mask, code, mach, flux in ((thin, BRANCH_THIN, flow.mach_e, flow.rho_e * flow.u_eff ** 2),
                                       (rarefied, BRANCH_RAREFIED, np.full(n, state.ma), np.full(n, state.freestream.rho * state.V ** 2))):
            if not mask.any():
                continue
            lam, tau = thin_film_mode(mach[mask], flux[mask], liq)
            ok = (np.isfinite(lam) & (lam > 0.0) & (tau > 0.0) & supercritical[mask]
                  & wave_fits(lam, None if L is None else L[mask]))
            rr = np.where(ok, lam / 4.0, np.nan)
            rate = np.where(ok, liq.rho * np.minimum(b[mask], lam / 8.0) / np.where(ok, tau, 1.0), 0.0)
            r[mask], mdot[mask], branch[mask] = rr, rate, code
            growth[mask] = np.where(ok, tau, np.nan)
        # (iv) the front-surface Rayleigh-Taylor mode, now *applied* rather than only reported.
        #
        # Girin & Kopyt assign the front surface -- where the flow influence is negligible, V_0 -> 0 and sin Theta -> 1 --
        # to this aperiodic mode rather than to their side mode, and with Girin's critical angle now gating the shear
        # branches the stagnation cap has no other mechanism at all. Three modelling decisions, none of them theirs,
        # because their paper gives no mass-loss rate for the front surface (its Table 2 has no mdot column):
        #   * the gate is *both* criteria, because each supplies something the other does not. The depth criterion
        #     W cos(phi) h^2 rho_l > 3 Sigma -- equivalently h > lambda*/2 pi -- is the one with a threshold, and it is
        #     what confines the mode to the nose: the normal deceleration falls as cos(phi) so the depth needed rises as
        #     1/sqrt(cos phi) while the melt gets shallower. The bounded criterion has no threshold of its own (with a
        #     wide molten region k = k_max is always admissible and tanh(k h) > 0 for any depth, so it answers yes
        #     wherever there is liquid); what it supplies is the admissible wavelength and the growth rate, including
        #     the lateral fit and the finite-depth factor. Requiring only the bounded form made every wet patch
        #     Rayleigh-Taylor-unstable, which is how this was caught;
        #   * the rate carries their side-surface construction across. For the side they give r_d = lambda*/4 and
        #     tau_d = tau*, which this model writes mdot = rho_l min(b, r_d/2)/tau*; for the front they give
        #     r_d ~ lambda* and tau_d ~ tau*, so the same form is mdot = rho_l min(b, lambda*/2)/tau*. Since lambda*
        #     is centimetre-scale here and b millimetre-scale, it reduces to rho_l b / tau*: the film is shaken off
        #     within one growth time, which is what their text describes. No new closure constant is introduced.
        #   * the droplet radius is lambda* (their front-surface r_d), then capped as every branch is, by the film mass
        #     on the patch and a quarter of the body radius -- and here the film-mass cap is what actually binds.
        # The two modes are two modes of one interface, so the one that reaches tearing-off first breaks it up: the
        # shorter growth time takes the patch and the other does not act on the same film. That is also what keeps the
        # shear from being spent twice (the double-counting checklist of the plan's fact 29).
        w_n = flow.deceleration if deceleration_n is None else deceleration_n   # normal component, W cos(phi)
        if self.rt_spray:
            deep, _, _ = rayleigh_taylor(w_n, layer, liq)               # the threshold: h > lambda*/2 pi
            rt_on, lam_rt, tau_rt = rayleigh_taylor_bounded(w_n, layer, np.full(n, np.inf) if L is None else L, liq)
            rt_on = np.asarray(rt_on) & np.asarray(deep) & has_film & np.isfinite(tau_rt)
            with np.errstate(divide="ignore", invalid="ignore"):
                rate_rt = np.where(rt_on, liq.rho * np.minimum(b, lam_rt / 2.0) / np.where(rt_on, tau_rt, 1.0), 0.0)
            faster = rt_on & (rate_rt > 0.0) & (~np.isfinite(growth) | (tau_rt < growth))
            r = np.where(faster, lam_rt, r)
            mdot = np.where(faster, rate_rt, mdot)
            branch = np.where(faster, BRANCH_RT, branch)
            growth = np.where(faster, tau_rt, growth)
        unstable = mdot > 0.0
        dm = np.where(unstable, np.minimum(mdot * areas * dt, m_f), 0.0)
        # a droplet is never larger than the film on its patch nor than a quarter of the body radius
        # (the long near-critical waves Girin excludes; the edge dynamic pressure vanishes toward the equator)
        with np.errstate(invalid="ignore"):
            r = np.minimum(r, np.minimum((3.0 * m_f / (4.0 * np.pi * liq.rho)) ** (1.0 / 3.0), 0.25 * radius))
        with np.errstate(divide="ignore", invalid="ignore"):
            dn = np.where(dm > 0.0, dm / (4.0 / 3.0 * np.pi * liq.rho * np.where(dm > 0.0, r, 1.0) ** 3), 0.0)
        rt_active, _, _ = rayleigh_taylor(w_n, layer, liq)     # the unbounded criterion, still reported for continuity
        return SprayResult(branch, we_s, unstable, r, mdot, dm, dn, bool(np.any(rt_active & has_film)), delta_m, v_s, growth)


def source_rows(t, h, V, theta, centroids, t_hat, flow, res, b, liquid, state):
    """Source-table rows (list of lists in SOURCE_COLUMNS order) for the patches that released mass this step."""
    k = np.flatnonzero(res.dm > 0.0)
    if k.size == 0:
        return []
    r = res.r[k]
    we_d = state.freestream.rho * state.V ** 2 * 2.0 * r / liquid.sigma
    oh = liquid.mu / np.sqrt(liquid.rho * liquid.sigma * 2.0 * r)
    rows = np.column_stack([np.full(k.size, t), np.full(k.size, h / 1e3), np.full(k.size, V / 1e3), np.degrees(theta[k]),
                            centroids[k, 0], centroids[k, 1], centroids[k, 2], flow.closure[k], res.branch[k], b[k],
                            res.delta_m[k], res.we_s[k], r, res.dn[k], res.dm[k], res.v_s[k], t_hat[k, 0], t_hat[k, 1], t_hat[k, 2],
                            we_d, oh, (we_d > WE_BREAKUP).astype(float)])
    return rows.tolist()


def histogram(r, dn, dm):
    """(dn per bin, dM per bin) on BIN_EDGES for droplet radii r with counts dn and masses dm."""
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(r) & (r > 0.0)
    n_hist, _ = np.histogram(r[ok], BIN_EDGES, weights=np.asarray(dn)[ok])
    m_hist, _ = np.histogram(r[ok], BIN_EDGES, weights=np.asarray(dm)[ok])
    return n_hist, m_hist
```


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

