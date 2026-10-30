import math, pytest
from reentry_model.parent_follower_separation_interactions.laws import (EmissionState, LawSettings, classify, park_radius, taylor_radius,
    feldman_radius, plume_radius, sutherland_mu, pressure_parameter)

def state(**kw):
    base = dict(rho=1e-4, T=220.0, p=6.3, mu=sutherland_mu(220.0), mfp=1e-3, U=6000.0, d=2.813, cda=5.5, r0=1.4065)
    base.update(kw); return EmissionState(**base)

def test_park_radius_equals_r0_at_parent_and_scales_as_cube_root():
    assert park_radius(0.0, 1.4, 0.7, 1/3) == pytest.approx(1.4, rel=1e-12)
    s0 = (1/0.7) ** 3 * 1.4
    big = park_radius(8000.0 - s0, 1.4, 0.7, 1/3) / park_radius(1000.0 - s0, 1.4, 0.7, 1/3)
    assert big == pytest.approx(2.0, rel=1e-12)
    assert park_radius(500.0, 1.4, 0.9, 1/3) > park_radius(500.0, 1.4, 0.7, 1/3)

def test_taylor_radius_value_and_floor():
    L = math.sqrt(5.5)
    assert taylor_radius(1000.0, 5.5, 1.4, 0.66, 1/3) == pytest.approx(0.5 * 0.66 * L * (1000.0 / L) ** (1/3))
    assert taylor_radius(0.1, 5.5, 1.4, 0.66, 1/3) == pytest.approx(1.4)        # never below r0

def test_feldman_square_grows_linearly_and_bounds_order():
    s = state()
    r1, r2 = (feldman_radius(x, s, 3000.0, 0.8, 0.71, 1.0) for x in (1000.0, 2000.0))
    kappa = sutherland_mu(3000.0) / (s.rho * 0.71)
    assert r2**2 - r1**2 == pytest.approx(4 * kappa * 1000.0 / (0.8 * s.U), rel=1e-12)
    assert feldman_radius(500.0, s, s.T, 0.8, 0.71, 1.0) <= feldman_radius(500.0, s, 3000.0, 0.8, 0.71, 1.0)
    # Controller ruling R1: the default state is TRANSITIONAL (Re_d about 1.2e5, P about 1.1), where the low and
    # high bounds also include the Taylor radius and cannot coincide, so the coincidence check uses a laminar state.
    s_lam = state(rho=1e-5, p=0.5)
    lo = plume_radius(500.0, s_lam, LawSettings(T_axis_hot=s_lam.T), "low")
    hi = plume_radius(500.0, s_lam, LawSettings(T_axis_hot=s_lam.T), "high")
    assert lo.regime == hi.regime == "laminar"
    assert lo.r == pytest.approx(hi.r)

@pytest.mark.parametrize("kw,expected", [
    (dict(mfp=0.1), "rarefied"),                                  # Kn = 0.036 > 0.01
    (dict(rho=1e-5, p=0.5), "laminar"),                           # low Re and low P
    (dict(rho=1e-2, p=300.0), "turbulent"),                       # high Re and high P
    (dict(rho=1e-2, p=0.5), "transitional"),                      # Re says turbulent, P says laminar
    (dict(rho=1e-6, p=300.0), "transitional"),                    # Re says laminar, P says turbulent
    (dict(U=0.5), "too_slow"),
])
def test_classify_uses_only_flow_quantities(kw, expected):
    assert classify(state(**kw), LawSettings()) == expected

def test_small_body_is_rarefied_where_large_body_is_not():
    assert classify(state(d=0.005, r0=0.0025, cda=2e-5), LawSettings()) == "rarefied"
    assert classify(state(), LawSettings()) != "rarefied"

def test_thresholds_are_overridable():
    s = state(rho=1e-2, p=0.5)                                    # transitional by default
    assert classify(s, LawSettings(p_turb=1e-3, p_lam=1e-4)) == "turbulent"
    assert classify(state(mfp=0.1), LawSettings(kn_continuum=1.0)) != "rarefied"

def test_transitional_range_holds_all_three_laws():
    s = state(rho=1e-2, p=0.5)
    lo, hi = (plume_radius(800.0, s, LawSettings(), b) for b in ("low", "high"))
    assert lo.regime == hi.regime == "transitional"
    assert set(hi.laws) == {"feldman_low", "feldman_high", "taylor"} and lo.r <= hi.r

def test_rarefied_has_no_radius_unless_fallback():
    s = state(mfp=0.1)
    assert math.isnan(plume_radius(100.0, s, LawSettings(), "low").r)
    assert plume_radius(100.0, s, LawSettings(rarefied_fallback="taylor"), "low").r > 0

def test_too_slow_parent_gives_no_radius_and_no_error():
    assert math.isnan(plume_radius(100.0, state(U=0.0), LawSettings(), "high").r)

def test_pressure_parameter_units():
    # 1333.22 Pa x 1 m = 1 cm Hg x 100 cm
    assert pressure_parameter(state(p=1333.22, cda=1.0)) == pytest.approx(100.0)

def test_speed_and_distance_flags():
    r = plume_radius(1e5, state(rho=1e-2, p=300.0, U=7000.0), LawSettings(), "low")
    assert r.flags["speed_out_of_range"] and r.flags["distance_out_of_range"]
