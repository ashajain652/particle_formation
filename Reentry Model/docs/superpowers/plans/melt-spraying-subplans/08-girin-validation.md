# Sub-plan: Task 8 — Girin's published cases: the Girin-as-published driver and its analysis script

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 3950–4348). Read `00-shared-context.md` first.

**Depends on:** Task 4 (dispersion table), Task 5 (flow regime), and Task 7 (spraying).
**Produces:** a standalone driver and analysis/plotting script reproducing Girin's own published test cases — not this project's sphere — as an independent literature-validation track.
**Character:** verification/analysis, deliberately separate from the main coupled model (Tasks 9–14).
**Read before implementing:** none beyond the shared context. Because this never touches the main flight simulation, it can be executed in parallel with Task 9 once Tasks 4/5/7 have landed, rather than waiting behind Task 9.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan.

---

### Task 8: Girin's published cases — the Girin-as-published driver and its analysis script

**Files:**
- Create: `reentry_model/girin_case.py`, `analysis/girin_reference.py`
- Test: `tests/test_reentry_model_girin.py`

**Interfaces:**
- Consumes: `dispersion.DispersionTable` (Task 4), `surface_flow.ranger_psi` (Task 5), `spray.thin_film_mode/rayleigh_taylor/histogram/BIN_EDGES` (Task 7), `compare`'s plot style (Step 2), the two reference-value files (Task 7).
- Produces: `girin_case.GirinVariant(name, V, rho_m, mu_m, sigma, R0=3e-3)`, `dimensionless_numbers(v, rho_a=RHO_A, mu_a=MU_AIR, compression=COMPRESSION) -> (Re, We, GI, alpha, mu, p, t_ch)`, `critical_angle(GI, p, we_cr=4.62)`, `surface_state(...)`, `run_variant(v, k_r, k_t, we_cr, re_density="shock"|"ambient", ...) -> dict` (keys `GI, phi_cr_deg_eq3, t_f_min_us_tau0, r_um_90deg_tau0, N, r_med_um, r_min_um, r_max_um, t_sd_ms, tau_end, radii, counts, times, mass_history, ...`); constants `RHO_A, MU_AIR, COMPRESSION, RE_DENSITY_NAMES`. The script writes `girin2017.json`, `girin1994.json`, `girin_summary.md` and the plots under `reentry_model_output/verification_melt/girin/`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_girin.py`:

```python
"""girin_case.py: Girin (2017) Table 1 -- the exact tier (GI, phi_cr within 2 %), the integrated tier (t_f, N, r_med
within 30 % with the ambient-density Reynolds number) and the reported t_s.d. (spec Step 3 section 13.2 as amended)."""
import json
import os

import numpy as np
import pytest

from reentry_model import girin_case as gc
from helpers import REPO_ROOT

TABLE = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin2017_table1.json")))


def variant(name):
    d = TABLE["variants"][name]
    return gc.GirinVariant(name, d["V_ms"], d["rho_m"], d["mu_m"], d["sigma"]), d


def test_exact_tier_gi_and_critical_angle():
    for name in ("I", "II", "III"):
        v, d = variant(name)
        Re, We, GI, alpha, mu, p, t_ch = gc.dimensionless_numbers(v)
        assert Re == pytest.approx(d["Re_inf"], rel=2e-2) and GI == pytest.approx(d["GI"], rel=2e-2)
        assert t_ch * 1e3 == pytest.approx(d["t_ch_ms"], rel=2e-2)                                # alpha on the ambient density
        assert np.degrees(gc.critical_angle(GI, p)) == pytest.approx(d["phi_cr_deg"], rel=2e-2)   # We_cr = 4.62, his practical value
        assert np.degrees(gc.critical_angle(GI, p, 3.08)) < d["phi_cr_deg"]                        # 3.08 would give ~17 % smaller angles
    assert gc.critical_angle(0.3, 1.0) == np.pi / 2                                                # below GI 0.4: nothing unstable


@pytest.mark.parametrize("name", ["I", "II", "III"])
def test_integrated_tier_with_the_ambient_reynolds_number(name):
    v, d = variant(name)
    o = gc.run_variant(v, re_density="ambient")
    assert o["t_f_min_us_tau0"] == pytest.approx(d["t_f_us"], rel=0.3)          # measured x1.22 / x0.88 / x1.07
    assert o["N"] == pytest.approx(d["N"], rel=0.3)                              # measured x0.87 / x0.97 / x0.76
    if d["r_med_um"]:
        assert o["r_med_um"] == pytest.approx(d["r_med_um"], rel=0.3)            # measured x0.96 / x0.95
        assert o["r_min_um"] < d["r_med_um"] < o["r_max_um"]
    assert o["t_sd_ms"] > 0.0 and o["m_end_kg"] < 1e-3 * o["m0_kg"]
    assert o["radii"].size == o["counts"].size and o["mass_history"][-1, 1] < o["mass_history"][0, 1]
    assert 0.4 < o["t_sd_ms"] / d["t_sd_ms"] < 0.6 or name == "III"              # reported: half his duration (iron), far below for stone


def test_shock_density_option_and_argument_checks():
    v, d = variant("I")
    shock, ambient = gc.run_variant(v, re_density="shock"), gc.run_variant(v, re_density="ambient")
    assert shock["Re_used"] == pytest.approx(6.0 * ambient["Re_used"]) and shock["r_med_um"] < 0.5 * ambient["r_med_um"]
    assert shock["GI"] == ambient["GI"] and shock["phi_cr_deg_eq3"] == ambient["phi_cr_deg_eq3"]
    with pytest.raises(ValueError):
        gc.run_variant(v, re_density="wrong")
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_girin.py -q`
Expected: ImportError (`reentry_model.girin_case`).

- [ ] **Step 3: Create `reentry_model/girin_case.py`**

```python
"""Girin's (2017) published spraying cases run with his own simplifications ("Girin-as-published" mode, spec Step 3
section 13.2): a sphere of radius R0 in a constant gas density, potential-flow edge velocity V_a = 1.5 (V_inf - w)
sin phi, Ranger's boundary layer delta_a = 2.2 R Re_a^-1/2 Psi(phi) with Re_a = Re_inf R~ (1 - W), the conjugated-
layer relations of his Eq. (2), the dispersion table for Delta_f and Im Omega_f, no melt limit and no runoff, spherical
belts of width R dphi, the induction rule tau_ind = int dt / t_per >= 1 releasing R dphi / lambda_f tori of
cross-section radius r = k_r lambda_f, and the deceleration law W = 1 - exp(-C tau), C = 2 alpha^1/2, tau = t / t_ch,
t_ch = 2 R0 / (alpha^1/2 V_inf).

Densities (measured 2026-09-20 against his Table 1): We_inf uses the ambient density rho_inf = rho_a / 6 and Re_inf the
compressed density rho_a = 1e-7 g/cm3 -- that reproduces his GI = We_inf Re_inf^-1/2 (13.0 / 3.55 / 43.5) and his
t_ch (alpha = rho_inf / rho_m); the Reynolds number entering delta_a is a switch (`re_density`): "shock" (rho_a, as the
printed Re_inf = 529) or "ambient" (rho_inf); with "ambient" his N and r_med are reproduced within 15 % for the iron
variants, with "shock" the droplets come out 2.5 x smaller (the paper does not say which density its delta_a used).
The stripping rate per area is independent of delta_a (it scales with V_s only), so t_s.d. is the same either way:
0.5 x his printed value for the iron variants and far below it for the stony one, whose lambda_f (3.5 mm) exceeds the
body radius -- his belt discretisation and induction handling are unstated (see the Step 3 plan's measured facts)."""
from dataclasses import dataclass

import numpy as np

from . import dispersion
from .surface_flow import ranger_psi

RHO_A, MU_AIR = 1.0e-4, 6.8e-5          # kg/m3, Pa s (Girin 2017 section 5)
COMPRESSION = 6.0
RE_DENSITY_NAMES = ("shock", "ambient")


@dataclass
class GirinVariant:
    name: str
    V: float          # m/s
    rho_m: float      # kg/m3
    mu_m: float       # Pa s
    sigma: float      # N/m
    R0: float = 3.0e-3


def dimensionless_numbers(v, rho_a=RHO_A, mu_a=MU_AIR, compression=COMPRESSION):
    """(Re_inf, We_inf, GI, alpha, mu, p, t_ch) with Girin's density conventions."""
    rho_inf = rho_a / compression
    Re = rho_a * 2.0 * v.R0 * v.V / mu_a
    We = rho_inf * 2.0 * v.R0 * v.V ** 2 / v.sigma
    alpha, mu = rho_inf / v.rho_m, mu_a / v.mu_m
    p = 1.0 + (alpha * mu) ** (1.0 / 3.0)
    t_ch = 2.0 * v.R0 / (np.sqrt(alpha) * v.V)
    return Re, We, We / np.sqrt(Re), alpha, mu, p, t_ch


def critical_angle(GI, p, we_cr=dispersion.WE_CRITICAL_PRACTICAL):
    """phi_cr [rad] at tau = 0 from Girin's Eq. (3): 2.475 p^-2 sin^2 phi Psi(phi) GI = We_cr."""
    from scipy.optimize import brentq
    f = lambda phi: 2.475 / p ** 2 * np.sin(phi) ** 2 * ranger_psi(phi) * GI - we_cr
    if f(np.pi / 2 - 1e-9) <= 0.0:
        return np.pi / 2
    return brentq(f, 1e-4, np.pi / 2 - 1e-9)


def surface_state(v, phi, Re_a, R, W, alpha, mu, p, table):
    """Per-belt delta_a, delta_m, V_s, We_s, Delta_f, Im Omega_f, lambda_f, r, t_per, t_f at the current radius R."""
    delta_a = 2.2 * R * ranger_psi(phi) / np.sqrt(Re_a)
    delta_m = (alpha / mu ** 2) ** (1.0 / 3.0) * delta_a
    V_a = 1.5 * v.V * (1.0 - W) * np.sin(phi)
    V_s = (alpha * mu) ** (1.0 / 3.0) / p * V_a
    we_s = v.rho_m * V_s ** 2 * delta_m / v.sigma
    delta_f, im_f, _ = table(np.maximum(we_s, table.we_onset))
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = 2.0 * np.pi * delta_m / delta_f
        t_f = delta_m / (V_s * im_f)
    return delta_a, delta_m, V_s, we_s, delta_f, im_f, lam, t_f


def run_variant(v, k_r=0.17, k_t=1.1, we_cr=dispersion.WE_CRITICAL_PRACTICAL, re_density="shock", dphi_deg=0.5, dt_fraction=0.2,
                mass_fraction_end=1e-4, table=None):
    """Girin's spraying flight of one variant. Returns a dict with the exact-tier numbers (GI, phi_cr, t_f at tau = 0,
    r at tau = 0 and 90 degrees, t_ch) and the integrated ones (N, r_med, r_min, r_max, t_sd, mass-loss history and the
    droplet list)."""
    if re_density not in RE_DENSITY_NAMES:
        raise ValueError("re_density must be one of {}".format(RE_DENSITY_NAMES))
    table = table or dispersion.DispersionTable()
    Re, We, GI, alpha, mu, p, t_ch = dimensionless_numbers(v)
    Re_used = Re if re_density == "shock" else Re / COMPRESSION
    C = 2.0 * np.sqrt(alpha)
    phi = np.radians(np.arange(dphi_deg / 2.0, 90.0, dphi_deg))
    dphi = np.radians(dphi_deg)
    m0 = v.rho_m * 4.0 / 3.0 * np.pi * v.R0 ** 3
    # tau = 0 exact-tier quantities
    _, _, _, we0, _, _, lam0, tf0 = surface_state(v, phi, Re_used, v.R0, 0.0, alpha, mu, p, table)
    unstable0 = we0 > we_cr
    out = {"Re_inf": Re, "We_inf": We, "GI": GI, "alpha": alpha, "mu": mu, "p": p, "t_ch_s": t_ch, "Re_used": Re_used,
           "phi_cr_deg_eq3": float(np.degrees(critical_angle(GI, p, we_cr))),
           "phi_cr_deg_direct": float(np.degrees(phi[unstable0][0])) if unstable0.any() else None,
           "t_f_min_us_tau0": float(tf0[unstable0].min() * 1e6) if unstable0.any() else None,
           "r_um_90deg_tau0": float(k_r * lam0[-1] * 1e6), "we_s_90deg_tau0": float(we0[-1])}
    # the flight
    m, t = m0, 0.0
    tind = np.zeros_like(phi)
    r_rel, n_rel, t_rel, m_hist = [], [], [], [(0.0, m0)]
    while m > mass_fraction_end * m0:
        R = (3.0 * m / (4.0 * np.pi * v.rho_m)) ** (1.0 / 3.0)
        W = 1.0 - np.exp(-C * t / t_ch)
        Re_a = Re_used * (R / v.R0) * (1.0 - W)
        _, delta_m, V_s, we_s, delta_f, im_f, lam, t_f = surface_state(v, phi, Re_a, R, W, alpha, mu, p, table)
        active = we_s > we_cr
        if not active.any():
            break
        t_per = k_t * t_f
        dt = dt_fraction * float(t_per[active].min())
        tind[active] += dt / t_per[active]
        tind[~active] = 0.0
        release = active & (tind >= 1.0)
        if release.any():
            r = k_r * lam[release]
            n_tor = R * dphi / lam[release]
            volume = n_tor * 2.0 * np.pi ** 2 * r ** 2 * R * np.sin(phi[release])
            dn = volume / (4.0 / 3.0 * np.pi * r ** 3)
            r_rel.append(r); n_rel.append(dn); t_rel.append(np.full(r.size, t))
            m -= float((v.rho_m * volume).sum())
            tind[release] -= 1.0
            m_hist.append((t, m))
        t += dt
    r = np.concatenate(r_rel) if r_rel else np.zeros(0)
    n = np.concatenate(n_rel) if n_rel else np.zeros(0)
    if r.size:
        order = np.argsort(r)
        cn = np.cumsum(n[order])
        r_med = float(r[order][np.searchsorted(cn, 0.5 * cn[-1])])
    else:
        r_med = None
    out.update({"N": float(n.sum()), "r_med_um": r_med * 1e6 if r_med else None, "r_min_um": float(r.min() * 1e6) if r.size else None,
                "r_max_um": float(r.max() * 1e6) if r.size else None, "t_sd_ms": t * 1e3, "tau_end": t / t_ch, "m0_kg": m0, "m_end_kg": m,
                "radii": r, "counts": n, "times": np.concatenate(t_rel) if t_rel else np.zeros(0), "mass_history": np.array(m_hist)})
    return out
```


- [ ] **Step 4: Create `analysis/girin_reference.py`**

```python
"""Girin's published spraying cases against the model's transcription of his theory (spec Step 3 sections 13.2-13.3).

    "$PY" analysis/girin_reference.py [--outdir reentry_model_output/verification_melt/girin] [--re-density ambient|shock|both]

Runs the three variants of Girin (2017) Table 1 in Girin-as-published mode (reentry_model.girin_case) and checks the
Girin & Kopyt (1994) Tables 1-2 with the model's thin-film and Rayleigh-Taylor modes (reentry_model.spray). Writes
girin2017.json / girin1994.json (every number next to its published value and ratio), the size distributions dn(r),
dM(r) per variant with the table values marked, the mass-loss law, and log-log plots of the 1994 tables with
residuals, plus girin_summary.md. Thresholds (README, Step 3): exact tier GI and phi_cr within 2 %; integrated tier
t_f, N and r_med within 30 % with the ambient-density Reynolds number; t_s.d. reported."""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
from reentry_model import girin_case, spray  # noqa: E402
from reentry_model.compare import INK, MODEL_COLOR, MUTED, SECOND, apply_rcparams, strip_top_right_spines  # noqa: E402

DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt", "girin")
TABLE_2017 = os.path.join(REPO_ROOT, "data", "reference_values", "girin2017_table1.json")
TABLE_1994 = os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")


def ratio(model, ref):
    return None if (model is None or ref in (None, 0)) else float(model / ref)


def run_2017(outdir, densities):
    tab = json.load(open(TABLE_2017))
    results = {}
    for name, d in tab["variants"].items():
        v = girin_case.GirinVariant(name, d["V_ms"], d["rho_m"], d["mu_m"], d["sigma"])
        results[name] = {}
        for rd in densities:
            o = girin_case.run_variant(v, re_density=rd)
            entry = {"GI": (o["GI"], d["GI"]), "phi_cr_deg": (o["phi_cr_deg_eq3"], d["phi_cr_deg"]), "t_f_us": (o["t_f_min_us_tau0"], d["t_f_us"]),
                     "N": (o["N"], d["N"]), "r_med_um": (o["r_med_um"], d["r_med_um"]), "r_min_um": (o["r_min_um"], d["r_min_um"]),
                     "r_max_um": (o["r_max_um"], d["r_max_um"]), "t_sd_ms": (o["t_sd_ms"], d["t_sd_ms"]), "tau_end": (o["tau_end"], d["tau_end"]),
                     "t_ch_ms": (o["t_ch_s"] * 1e3, d["t_ch_ms"]), "Re_inf": (o["Re_inf"], d["Re_inf"])}
            results[name][rd] = {k: {"model": m, "published": p, "ratio": ratio(m, p)} for k, (m, p) in entry.items()}
            results[name][rd]["r_um_90deg_tau0"] = o["r_um_90deg_tau0"]
            # plots: size distributions and the mass-loss law
            apply_rcparams(plt)
            fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(14, 4.2))
            if o["radii"].size:
                n_hist, m_hist = spray.histogram(o["radii"], o["counts"], o["counts"] * 4.0 / 3.0 * np.pi * v.rho_m * o["radii"] ** 3)
                mid = np.sqrt(spray.BIN_EDGES[:-1] * spray.BIN_EDGES[1:]) * 1e6
                ax.loglog(mid[n_hist > 0], n_hist[n_hist > 0], color=MODEL_COLOR, lw=1.4, label="model")
                bx.loglog(mid[m_hist > 0], m_hist[m_hist > 0] * 1e3, color=MODEL_COLOR, lw=1.4, label="model")
                for value, label in ((d["r_med_um"], "published r_med"), (d["r_min_um"], "published range"), (d["r_max_um"], None)):
                    if value:
                        ax.axvline(value, color=INK, lw=0.8, ls=":", label=label)
                        bx.axvline(value, color=INK, lw=0.8, ls=":")
                mh = o["mass_history"]
                cx.plot(mh[:, 0] * 1e3, mh[:, 1] / mh[0, 1], color=MODEL_COLOR, lw=1.4, label="model")
                t_sd = mh[-1, 0]
                cx.plot(mh[:, 0] * 1e3, (1.0 - np.minimum(mh[:, 0] / t_sd, 1.0)) ** 3, color=MUTED, lw=1.0, ls="--", label="(1 - t/t_s.d.)^3")
                cx.axvline(d["t_sd_ms"], color=INK, lw=0.8, ls=":", label="published t_s.d.")
            ax.set_xlabel("r [um]"); ax.set_ylabel("dn per bin"); ax.legend(frameon=False)
            bx.set_xlabel("r [um]"); bx.set_ylabel("dM per bin [g]")
            cx.set_xlabel("t [ms]"); cx.set_ylabel("m / m0"); cx.legend(frameon=False)
            ax.set_title("Girin 2017 variant {} ({}), Re on {} density".format(name, d["description"], rd), color=SECOND, fontsize=10)
            for a in (ax, bx, cx):
                strip_top_right_spines(a)
            fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin2017_variant_{}_{}.png".format(name, rd)), dpi=150); plt.close(fig)
    return results


def run_1994(outdir):
    tab = json.load(open(TABLE_1994))
    t1, t2 = tab["table1"], tab["table2"]

    class Liquid:
        rho, sigma = t1["rho1_kgm3"], t1["sigma_Nm"]
        mu = 0.0
    rho2, V0 = np.array(t1["rho2_kgm3"]), np.array(t1["V0_ms"])
    r_tab, tau_tab, m_tab = (np.array(t1[k]) for k in ("r_d_m", "tau_d_s", "m_kgm2s"))
    lam, tau = spray.thin_film_mode(np.full((3, 2), t1["mach"]), rho2[:, None] * V0[None, :] ** 2, Liquid)
    f_s = float((lam / 4.0)[0, 0] / r_tab[0, 0])                     # their shock-layer factor, fitted on the first entry
    r_model, tau_model = lam / 4.0 / f_s, tau / f_s ** 1.5
    m_model = Liquid.rho * r_model / (2.0 * tau_model)
    W = np.array(t2["W_ms2"])
    lam_rt, tau_rt = zip(*[spray.rayleigh_taylor(w, np.array([1.0]), Liquid)[1:] for w in W])
    lam_rt, tau_rt = np.array([x[0] for x in lam_rt]), np.array([x[0] for x in tau_rt])
    res = {"table1": {"shock_layer_factor": f_s, "r_d_ratio": (r_model / r_tab).tolist(), "tau_d_ratio": (tau_model / tau_tab).tolist(),
                      "m_ratio": (m_model / m_tab).tolist(), "r_d_model_m": r_model.tolist(), "tau_d_model_s": tau_model.tolist()},
           "table2": {"lambda_star_model_m": lam_rt.tolist(), "lambda_star_ratio_to_eq14_column": (lam_rt / np.array(t2["lambda_star_m_eq14"])).tolist(),
                      "lambda_star_ratio_to_printed": (lam_rt / np.array(t2["lambda_star_m_as_printed"])).tolist(),
                      "tau_star_model_s": tau_rt.tolist(), "tau_star_ratio": (tau_rt / np.array(t2["tau_star_s"])).tolist()}}
    apply_rcparams(plt)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), gridspec_kw={"height_ratios": [3, 1]})
    for j, (name, model, published, unit) in enumerate((("r_d", r_model, r_tab, "m"), ("tau_d", tau_model, tau_tab, "s"), ("m", m_model, m_tab, "kg/m2/s"))):
        ax, rx = axes[0, j], axes[1, j]
        for k, V in enumerate(V0):
            ax.loglog(rho2, published[:, k], "o", color=INK, label="published, V0 = {:.0f} km/s".format(V / 1e3))
            ax.loglog(rho2, model[:, k], "--", color=MODEL_COLOR, label="model x f_s" if k == 0 else None)
            rx.semilogx(rho2, 100.0 * (model[:, k] / published[:, k] - 1.0), "s-", color=MODEL_COLOR if k == 0 else MUTED, lw=0.8)
        ax.set_ylabel("{} [{}]".format(name, unit)); ax.legend(frameon=False, fontsize=8); ax.set_title("Girin & Kopyt 1994 Table 1", color=SECOND, fontsize=10)
        rx.axhline(0.0, color=MUTED, lw=0.6); rx.set_xlabel("rho_2 [kg/m3]"); rx.set_ylabel("model/published - 1 [%]")
        strip_top_right_spines(ax); strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin1994_table1.png"), dpi=150); plt.close(fig)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2))
    ax.loglog(W, t2["lambda_star_m_eq14"], "o", color=INK, label="Eq. (14) (printed x 10)")
    ax.loglog(W, t2["lambda_star_m_as_printed"], "x", color=MUTED, label="Table 2 as printed")
    ax.loglog(W, lam_rt, "--", color=MODEL_COLOR, label="model")
    ax.set_xlabel("W [m/s2]"); ax.set_ylabel("lambda* [m]"); ax.legend(frameon=False, fontsize=8); ax.set_title("Girin & Kopyt 1994 Table 2", color=SECOND, fontsize=10)
    bx.loglog(W, t2["tau_star_s"], "o", color=INK, label="published"); bx.loglog(W, tau_rt, "--", color=MODEL_COLOR, label="model")
    bx.set_xlabel("W [m/s2]"); bx.set_ylabel("tau* [s]"); bx.legend(frameon=False, fontsize=8)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin1994_table2.png"), dpi=150); plt.close(fig)
    return res


def summary_markdown(r2017, r1994, densities):
    lines = ["| variant | Re density | GI | phi_cr [deg] | t_f [us] | N | r_med [um] | range [um] | t_s.d. [ms] |", "|---|---|---|---|---|---|---|---|---|"]
    fmt = lambda e, f="{:.3g}": "{} ({}{})".format(f.format(e["model"]) if e["model"] is not None else "n/a", f.format(e["published"]) if e["published"] is not None else "n/a",
                                                  ", x{:.2f}".format(e["ratio"]) if e["ratio"] else "")
    for name, per in r2017.items():
        for rd in densities:
            e = per[rd]
            rng = "{}-{}".format("{:.1f}".format(e["r_min_um"]["model"]) if e["r_min_um"]["model"] else "n/a", "{:.1f}".format(e["r_max_um"]["model"]) if e["r_max_um"]["model"] else "n/a")
            rng += " ({}-{})".format(e["r_min_um"]["published"] or "n/a", e["r_max_um"]["published"] or "n/a")
            lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(name, rd, fmt(e["GI"], "{:.2f}"), fmt(e["phi_cr_deg"], "{:.1f}"), fmt(e["t_f_us"], "{:.1f}"),
                                                                       fmt(e["N"], "{:.2e}"), fmt(e["r_med_um"], "{:.1f}"), rng, fmt(e["t_sd_ms"], "{:.2f}")))
    lines += ["", "Values in parentheses: published (Girin 2017 Table 1) and model/published.", "",
              "Girin & Kopyt 1994 Table 1: shock-layer factor {:.2f}; r_d ratios {}; tau_d ratios {}; m ratios {}".format(
                  r1994["table1"]["shock_layer_factor"], np.round(r1994["table1"]["r_d_ratio"], 3).tolist(), np.round(r1994["table1"]["tau_d_ratio"], 3).tolist(),
                  np.round(r1994["table1"]["m_ratio"], 3).tolist()),
              "Girin & Kopyt 1994 Table 2: lambda* / Eq. (14) {}; lambda* / printed {}; tau* ratios {}".format(
                  np.round(r1994["table2"]["lambda_star_ratio_to_eq14_column"], 3).tolist(), np.round(r1994["table2"]["lambda_star_ratio_to_printed"], 3).tolist(),
                  np.round(r1994["table2"]["tau_star_ratio"], 3).tolist())]
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--re-density", choices=("ambient", "shock", "both"), default="both")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    densities = ("ambient", "shock") if args.re_density == "both" else (args.re_density,)
    r2017 = run_2017(args.outdir, densities)
    r1994 = run_1994(args.outdir)
    json.dump(r2017, open(os.path.join(args.outdir, "girin2017.json"), "w"), indent=2)
    json.dump(r1994, open(os.path.join(args.outdir, "girin1994.json"), "w"), indent=2)
    md = summary_markdown(r2017, r1994, densities)
    open(os.path.join(args.outdir, "girin_summary.md"), "w").write(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```


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

