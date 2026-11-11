"""The viscosity of liquid and slurry AA7075 (sub-plan 21 Task 3; the Step 4 plan's section 4): Li et al.'s slurry law,
its shear rate held at 367 1/s, the bridge to the liquid at the liquidus, and rigid metal at half solid."""
import numpy as np
import pytest

from reentry_model import material, rheology

MU_L = 1.3e-3


def scheil():
    return material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])


def test_at_rest_li_et_al_runs_from_0_87_to_5_6_pa_s_670_to_4300_times_the_liquid():
    assert rheology.li_viscosity(0.0, 0.0) == pytest.approx(0.871)
    assert rheology.li_viscosity(0.5, 0.0) == pytest.approx(0.871 * np.exp(0.5 * 3.7311))         # 5.63 Pa s
    assert rheology.li_viscosity(0.0, 0.0) / MU_L == pytest.approx(670.0)
    assert rheology.li_viscosity(0.5, 0.0) / MU_L == pytest.approx(4328.0, rel=1e-3)
    # the published form at a measured point: f_s = 0.1, 61 1/s
    assert rheology.li_viscosity(0.1, 61.0) == pytest.approx((0.871 - 0.00849 * 61.0 ** 0.74924) * np.exp(0.37311))


def test_the_shear_rate_is_held_at_367_per_second_where_the_fit_still_holds():
    f_s = np.linspace(0.1, 0.5, 5)
    np.testing.assert_array_equal(rheology.li_viscosity(f_s, 450.0), rheology.li_viscosity(f_s, 367.0))
    assert np.all(rheology.li_viscosity(f_s, 367.0) > 0.0)                    # the bracket's zero is at 483 1/s
    rates = np.array([0.0, 61.0, 150.0, 367.0])
    assert np.all(np.diff(rheology.li_viscosity(0.3, rates)) < 0.0)           # shear-thinning


def test_the_bridge_runs_log_linearly_from_the_slurry_to_the_liquid_at_the_liquidus():
    f_b = 0.0975
    top = rheology.li_viscosity(f_b, 100.0)
    eta = lambda T, f_s: float(rheology.viscosity(T, f_s, 100.0, MU_L, f_b))
    assert eta(906.4, f_b) == pytest.approx(top)
    assert eta(907.2, 0.05) == pytest.approx(np.sqrt(top * MU_L))             # half way: the geometric mean
    assert eta(908.0, 0.0) == MU_L and eta(950.0, 0.0) == MU_L


def test_metal_half_solid_or_more_is_rigid():
    assert np.isinf(rheology.viscosity(895.0, 0.5, 50.0, MU_L, 0.0975))
    assert np.isinf(rheology.viscosity(850.0, 0.9, 50.0, MU_L, 0.0975))
    assert rheology.viscosity(895.2, 0.4999, 50.0, MU_L, 0.0975) == pytest.approx(rheology.li_viscosity(0.4999, 50.0))


def test_on_aa7075_s_scheil_curve_the_viscosity_falls_without_a_break_from_the_rigid_front_to_the_liquid():
    mat = scheil()
    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    assert f_b == pytest.approx(0.0975, abs=1e-4)                              # Scheil's curve at 906.4 K
    T = np.linspace(mat.T_rigid + 1e-3, 910.0, 4001)
    eta = rheology.viscosity(T, 1.0 - mat.liquid_fraction(T), 100.0, mat.liquid.mu, f_b)
    assert np.all(np.isfinite(eta)) and np.all(np.diff(eta) <= 1e-12 * eta[:-1])
    for T_edge in (rheology.BRIDGE_LOW, rheology.BRIDGE_HIGH):                 # no jump where the bridge begins or ends
        Ts = np.array([T_edge - 1e-9, T_edge + 1e-9])
        lo, hi = rheology.viscosity(Ts, 1.0 - mat.liquid_fraction(Ts), 100.0, mat.liquid.mu, f_b)
        assert hi == pytest.approx(lo, rel=1e-6)
