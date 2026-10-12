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
