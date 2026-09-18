"""heating.py: distributions, bridging, the SESAM-equivalent mode and the physics mode at the 100 mm reference start."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from reentry_model import atmosphere, gas, heating
from reentry_model import trajectory as tj

# 100 mm sphere at the reference start: US76 77.5 km (rho 2.727e-5 kg/m3, T 203.6 K), 7.5 km/s, Kn 0.0298, Ma 26.2
RHO, T_INF, V, R, KN, MA = 2.727e-5, 203.6, 7500.0, 0.05, 0.0298, 26.2


def state(rho=RHO, T=T_INF, V=V, kn=KN, ma=MA):
    fs = atmosphere.Freestream(rho, T, rho / 28.9644e-3 * 8.314462618 * T, 28.9644 * 1.66053906660e-27, np.zeros(3))
    return tj.AeroState(77.5e3, 0.0, 0.0, fs, np.array([V, 0.0, 0.0]), V, kn, ma, 0.9, np.zeros(3), 0.5 * rho * V * V)


def windward_integral(shape, **kw):
    """Integral of the shape over the windward hemisphere divided by the whole sphere area."""
    return 0.5 * quad(lambda th: float(shape(np.array([th]), **kw)[0]) * math.sin(th), 0.0, math.pi / 2)[0]


def test_lees_distribution_limits_and_integrals():
    assert heating.lees_shape(np.array([0.0]), np.inf)[0] == 1.0
    assert heating.lees_shape(np.array([1e-3, 2.0]), 20.0)[0] == pytest.approx(1.0, abs=1e-5)
    assert heating.lees_shape(np.array([2.0, 3.0]), 20.0).tolist() == [0.0, 0.0]                  # leeward
    assert windward_integral(heating.lees_shape, mach=np.inf) == pytest.approx(0.196, abs=1e-3)
    assert windward_integral(heating.lees_shape, mach=10.0) == pytest.approx(0.200, abs=1e-3)
    assert windward_integral(heating.lees_shape, mach=5.0) == pytest.approx(0.209, abs=1e-3)
    assert heating.lees_shape(np.array([math.pi / 4]), np.inf)[0] == pytest.approx(0.5965, abs=1e-3)
    assert windward_integral(heating.cosine_shape) == pytest.approx(0.25, abs=1e-9)
    assert heating.lees_shape(np.array([0.5]), 0.3).tolist() == heating.lees_shape(np.array([0.5]), 1.0).tolist()   # Ma clamp


def test_matting_bridge_limits_and_closed_forms():
    assert heating.matting_bridge(1.0, 1e6) == pytest.approx(1.0, rel=1e-5)               # q_fm << q_c -> q_fm
    assert heating.matting_bridge(1e6, 1.0) == pytest.approx(1.0, rel=1e-9)               # q_fm >> q_c -> q_c
    assert heating.matting_bridge(3.0, 2.0) == pytest.approx(2.0 * (1.0 - math.exp(-1.5)))
    x = math.sqrt(2.0 * 3.0 / 2.0)
    assert heating.matting_bridge(3.0, 2.0, n=2.0) == pytest.approx(2.0 * (1.0 - (1.0 + x) * math.exp(-x)))
    assert heating.matting_bridge(3.0, 2.0, n=1.5) < heating.matting_bridge(3.0, 2.0, n=1.0)
    assert heating.matting_bridge(1.0, 0.0) == 0.0
    ratios = [heating.matting_bridge(q, 1.0) for q in (0.1, 0.5, 1.0, 2.0, 5.0)]
    assert all(b > a for a, b in zip(ratios, ratios[1:])) and ratios[-1] < 1.0


def test_sesam_equivalent_mode_is_uniform_and_totals_to_the_shape_factor(coarse_sphere_mesh):
    surface = coarse_sphere_mesh.surface()
    theta = surface.angles_to([1.0, 0.0, 0.0])
    model = heating.SesamEquivalentHeating()
    res = model.evaluate(state(), theta, np.full(surface.n_patches, 300.0), R)
    F_h = heating.aero.SesamHeatTable()(KN)
    hot_wall_300 = 1.0 - 1004.5 * (300.0 - T_INF) / (0.5 * V * V)                              # 0.9966 at the start
    q_stag = gas.dkr(RHO, V, R) * F_h * hot_wall_300
    assert 0.95 < F_h < 1.0 and res.blend == F_h and res.q_stag == pytest.approx(q_stag) and res.q_stag_c == pytest.approx(gas.dkr(RHO, V, R))
    assert np.all(res.q_conv == res.q_conv[0]) and res.q_conv[0] == pytest.approx(0.27471 * q_stag)
    assert res.total(surface.areas) == pytest.approx(0.27471 * q_stag * surface.area, rel=1e-12)
    assert res.total(surface.areas) == pytest.approx(16244.0, rel=1e-2)                      # SESAM's own t = 0 value, 16243.74 W
    assert model.evaluate(state(rho=0.0), theta, np.full(surface.n_patches, 300.0), R).total(surface.areas) == 0.0
    hot = model.evaluate(state(), theta, np.full(surface.n_patches, 2000.0), R)
    assert np.all(hot.q_conv == pytest.approx(res.q_conv * (1.0 - 1004.5 * (2000.0 - T_INF) / (0.5 * V * V)) / hot_wall_300))   # SESAM's hot-wall factor
    lumped = model.evaluate(state(), theta, np.full(surface.n_patches, 2000.0), R, T_mean=1000.0)
    assert lumped.q_conv[0] == pytest.approx(res.q_conv[0] * (1.0 - 1004.5 * (1000.0 - T_INF) / (0.5 * V * V)) / hot_wall_300)
    slow = model.evaluate(state(rho=1e-3, V=2000.0, kn=1e-5, ma=6.7), theta, np.full(surface.n_patches, 2300.0), R)
    assert np.all(slow.q_conv == 0.0)                                                          # clamped: h_w > h_s
    subsonic = model.evaluate(state(rho=1e-2, V=250.0, kn=1e-6, ma=0.8), theta, np.full(surface.n_patches, 1500.0), R)
    assert subsonic.q_conv[0] == pytest.approx(0.5 * 0.27471 * gas.dkr(1e-2, 250.0, R) * heating.aero.SesamHeatTable()(1e-6))
    fm = model.evaluate(state(rho=4e-8, T=300.0, V=7500.0, kn=35.0, ma=22.0), theta, np.full(surface.n_patches, 300.0), 0.025)
    assert fm.blend == pytest.approx(0.060, abs=3e-3)
    assert 0.6 < fm.total(surface.areas) / (0.25 * 0.5 * 4e-8 * 7500.0 ** 3 * surface.area) < 1.0   # SESAM: 0.78 x the cos-theta average at Kn 35


def test_sesam_heat_table():
    table = heating.aero.SesamHeatTable()
    assert table(1e-4) == pytest.approx(1.0077) and table(100.0) == pytest.approx(0.0594) and table(0.0) == pytest.approx(1.0077)
    assert table(1.0) == pytest.approx(0.14, abs=0.01) and table(0.1) == pytest.approx(0.70, abs=0.05)
    values = [table(10 ** x) for x in np.linspace(-2.7, 1.6, 50)]
    assert all(b <= a + 1e-12 for a, b in zip(values, values[1:]))                            # monotone in Kn


@pytest.fixture(scope="module")
def air():
    pytest.importorskip("cantera")
    return gas.EquilibriumAir()


@pytest.mark.parametrize("stagnation", heating.STAGNATION_NAMES)
def test_physics_mode_distribution_and_totals(air, coarse_sphere_mesh, stagnation):
    surface = coarse_sphere_mesh.surface()
    theta = surface.angles_to([1.0, 0.0, 0.0])
    model = heating.PhysicsHeating(stagnation=stagnation, air=air)
    res = model.evaluate(state(), theta, np.full(surface.n_patches, 300.0), R)
    assert res.q_stag < res.q_stag_c and 0.0 < res.blend < 0.2                               # Matting: q_fm/q_c ~ 3 -> w ~ 0.05
    assert res.blend == pytest.approx(math.exp(-res.q_stag_fm / res.q_stag_c), rel=1e-9)
    expected = res.q_stag * ((1.0 - res.blend) * 0.196 + res.blend * 0.25) * surface.area
    assert res.total(surface.areas) == pytest.approx(expected, rel=2e-2)
    assert np.all(res.q_conv[theta > math.pi / 2] == 0.0) and res.q_conv[int(np.argmin(theta))] == pytest.approx(res.q_stag, rel=1e-2)
    if stagnation == "dkr":
        assert res.q_stag_c == pytest.approx(gas.dkr(RHO, V, R), rel=1e-3)                     # cold wall: factor 1
    hot = model.evaluate(state(), theta, np.full(surface.n_patches, 2000.0), R)
    assert 0.85 < hot.total(surface.areas) / res.total(surface.areas) < 0.97                  # ~7 % hot-wall reduction


def test_physics_mode_free_molecular_limit_and_sesam_bridging(air, coarse_sphere_mesh):
    surface = coarse_sphere_mesh.surface()
    theta = surface.angles_to([1.0, 0.0, 0.0])
    T_wall = np.full(surface.n_patches, 300.0)
    fm = heating.PhysicsHeating(stagnation="sutton-graves", air=air).evaluate(state(rho=1e-10, kn=1e4), theta, T_wall, R)
    assert fm.blend > 0.99 and fm.q_stag == pytest.approx(fm.q_stag_fm, rel=2e-2)
    assert fm.total(surface.areas) == pytest.approx(0.25 * fm.q_stag * surface.area, rel=3e-2)  # cos(theta) shape
    sesam = heating.PhysicsHeating(stagnation="sutton-graves", bridging="sesam-table", air=air).evaluate(state(), theta, T_wall, R)
    assert sesam.blend == pytest.approx(heating.aero.SesamTable()(KN))
    with pytest.raises(ValueError):
        heating.PhysicsHeating(stagnation="magic", air=air)
    with pytest.raises(NotImplementedError):
        heating.TabulatedHeating("table.nc")
    assert isinstance(heating.heating_by_name("sesam"), heating.SesamEquivalentHeating)
