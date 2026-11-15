import math, numpy as np, pytest
from reentry_model.parent_follower_separation_interactions.affix import (affix_regime, theta_inn, theta_bif, AffixSettings,
    release_offset, far_wake_offset)
from reentry_model.parent_follower_separation_interactions.laws import EmissionState, sutherland_mu

E = EmissionState(rho=1.85e-5, T=199.0, p=1.05, mu=sutherland_mu(199.0), mfp=4e-3, U=7600.0, d=2.813,
                  cda=5.5, r0=1.4065)          # roughly the Generic Spacecraft at 80 km

@pytest.mark.parametrize("r", [1/16, 1/4])
def test_thresholds_follow_register(r):
    assert theta_inn(r) == pytest.approx(122 * r + 47)
    assert theta_bif(r) == pytest.approx(100 * math.exp(-1.618 * r) + 40)

def test_r_above_rstar_has_no_inner_threshold():
    assert math.isnan(theta_inn(0.5)) and theta_bif(0.5) == pytest.approx(100 * math.exp(-0.809) + 40)

def test_windward_small_sphere_far_wake_and_stagnation_escape():
    d2 = 0.01
    assert affix_regime(d2, E.d, 90.0, E, 2813.0).zone == "far_wake"
    assert affix_regime(d2, E.d, 180.0, E, 2813.0).zone == "escape"
    tb = theta_bif(d2 / E.d)
    assert affix_regime(d2, E.d, tb + 0.5, E, 2813.0).zone == "shock_band"
    assert affix_regime(d2, E.d, 30.0, E, 2813.0).zone == "trapped"

def test_interaction_times():
    d2 = 0.176; r = d2 / E.d
    q = 0.5 * E.rho * E.U**2; scale = E.d * math.sqrt(2813.0 / q)
    esc = affix_regime(d2, E.d, 180.0, E, 2813.0)
    assert esc.t_fz == pytest.approx(scale * (0.88*r + 0.32 + 0.01*math.exp(0.125*(180 - theta_bif(r) - 10))))
    fw = affix_regime(d2, E.d, 90.0, E, 2813.0)
    A0 = -7.2*r**2 + 3*r + 0.1
    assert fw.t_fz == pytest.approx(scale * (A0 + 0.5*math.exp(-0.05*(90 - theta_inn(r)))))
    assert affix_regime(d2, E.d, 10.0, E, 2813.0).t_fz == 0.0

def test_validity_flags():
    f = affix_regime(0.01, E.d, 90.0, E, 2813.0, parent_mass=740.0).flags
    assert f["r_below_fit"] and f["rarefied_follower"] and not f["rarefied_parent"] and f["unequal_density"]
    slow = EmissionState(**{**E.__dict__, "U": 800.0})
    assert affix_regime(0.3, E.d, 90.0, slow, 2813.0).flags["not_hypersonic"]

def test_theta_range_and_no_altitude_argument():
    import inspect
    with pytest.raises(ValueError):
        affix_regime(0.01, E.d, 181.0, E, 2813.0)
    assert not any("alt" in p or p in ("h", "h_km") for p in inspect.signature(affix_regime).parameters)

def test_release_offsets():
    v = np.array([1.0, 0, 0]); up = np.array([0, 0, 1.0])
    assert np.allclose(far_wake_offset(2.813, v), [-5.626, 0, 0])
    assert np.allclose(release_offset(180.0, 2.813, v, up), [1.4065, 0, 0])
    assert np.allclose(release_offset(90.0, 2.813, v, up), [0, 0, 1.4065], atol=1e-12)
    assert np.linalg.norm(release_offset(135.0, 2.813, v, up)) == pytest.approx(1.4065)

# --- additions beyond the brief's tests -----------------------------------------------------------------------

def test_band_edges_and_zone_boundaries():
    d2 = 0.176; r = d2 / E.d; tb = theta_bif(r); ti = theta_inn(r); s = AffixSettings().s_w_deg
    z = lambda th: affix_regime(d2, E.d, th, E, 2813.0).zone
    assert z(ti - 1e-9) == "trapped" and z(ti) == "far_wake"
    assert z(tb - s - 1e-6) == "far_wake" and z(tb - s) == "shock_band"
    assert z(tb + s) == "shock_band" and z(tb + s + 1e-6) == "escape"
    assert math.isnan(affix_regime(d2, E.d, tb, E, 2813.0).t_hat)

def test_trapped_has_zero_t_hat_and_result_fields():
    res = affix_regime(0.176, E.d, 10.0, E, 2813.0)
    assert res.t_hat == 0.0 and res.theta_deg == 10.0 and res.r == pytest.approx(0.176 / E.d)
    assert res.theta_inn == pytest.approx(theta_inn(res.r)) and res.theta_bif == pytest.approx(theta_bif(res.r))
    assert math.isnan(res.flags["parent_bulk_density"])

def test_large_follower_decided_by_theta_bif_alone():
    d2 = 0.8 * E.d; r = 0.8; tb = theta_bif(r)
    lo = affix_regime(d2, E.d, tb - 10.0, E, 2813.0)
    assert lo.zone == "far_wake" and lo.flags["r_above_rstar"] and math.isnan(lo.theta_inn) and math.isnan(lo.t_hat)
    assert affix_regime(d2, E.d, 10.0, E, 2813.0).zone == "far_wake"   # no trapped zone above r*
    assert affix_regime(d2, E.d, tb + 10.0, E, 2813.0).zone == "escape"
    assert affix_regime(d2, E.d, tb + 1.0, E, 2813.0).zone == "shock_band"

def test_equal_density_parent_mass_and_fit_range_flags():
    m = 2813.0 * math.pi * E.d**3 / 6
    res = affix_regime(0.3, E.d, 90.0, E, 2813.0, parent_mass=m * 1.005)
    assert not res.flags["unequal_density"] and res.flags["parent_bulk_density"] == pytest.approx(2813.0 * 1.005)
    assert affix_regime(0.3, E.d, 90.0, E, 2813.0, parent_mass=m * 1.02).flags["unequal_density"]
    assert not res.flags["r_below_fit"] and not res.flags["not_hypersonic"]
    assert affix_regime(0.1, E.d, 90.0, E, 2813.0).flags["r_below_fit"]

def test_settings_are_honoured():
    d2 = 0.176; tb = theta_bif(d2 / E.d)
    assert affix_regime(d2, E.d, tb + 3.0, E, 2813.0).zone == "escape"
    assert affix_regime(d2, E.d, tb + 3.0, E, 2813.0, settings=AffixSettings(s_w_deg=3.0)).zone == "shock_band"
    assert affix_regime(0.3, E.d, 90.0, E, 2813.0, settings=AffixSettings(kn_continuum=1e-4)).flags["rarefied_parent"]
    assert affix_regime(0.3, E.d, 90.0, E, 2813.0, settings=AffixSettings(mach_min=40.0)).flags["not_hypersonic"]

@pytest.mark.parametrize("args", [(0.0, 2.813, 90.0), (-1.0, 2.813, 90.0), (0.1, 0.0, 90.0), (0.1, 2.813, -0.1),
                                  (0.1, 2.813, float("nan"))])
def test_bad_inputs_raise(args):
    with pytest.raises(ValueError):
        affix_regime(*args[:2], args[2], E, 2813.0)
    with pytest.raises(ValueError):
        affix_regime(0.1, 2.813, 90.0, E, 0.0)

def test_release_offset_follows_up_vector_and_rejects_parallel():
    v = np.array([0.0, 3.0, 0.0]); up = np.array([1.0, 1.0, 0.0])      # up has a component along v
    off = release_offset(90.0, 2.0, v, up)
    assert np.allclose(off, [1.0, 0, 0]) and abs(off @ v) < 1e-12
    off = release_offset(135.0, 2.0, v, up)
    assert np.allclose(off, [math.sin(math.radians(135)), -math.cos(math.radians(135)), 0.0])
    with pytest.raises(ValueError):
        release_offset(90.0, 2.0, v, 2.0 * v)
    with pytest.raises(ValueError):
        release_offset(90.0, 2.0, np.zeros(3), up)
