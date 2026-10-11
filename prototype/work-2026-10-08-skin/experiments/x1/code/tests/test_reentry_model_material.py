"""material.py: DRAMA material tables, enthalpy and its inverse."""
import filecmp
import os

import numpy as np
import pytest

from reentry_model import material
from helpers import REPO_ROOT


def test_packaged_material_is_the_wrapper_file():
    assert filecmp.cmp(material.DEFAULT_MATERIAL, os.path.join(REPO_ROOT, "data", "user_materials", "AA7075_nomelt.json"), shallow=False)


def test_default_material_tables():
    m = material.Material.from_drama_json()
    assert m.name == "AA7075_nomelt" and m.rho == 2813.0 and m.emissivity == 0.4
    assert m.cp(293.0) == 877.5 and m.cp(850.0) == 1131.6 and m.cp(2000.0) == 1131.6       # held above 850 K
    assert m.k(293.0) == 163.89 and m.k(850.0) == 128.19 and m.k(2000.0) == 128.19
    assert m.cp(303.0) == pytest.approx(881.5) and m.k(100.0) == 163.89                     # linear between, clamped below
    assert np.all(m.liquid_fraction(np.array([300.0, 900.0])) == 0.0) and m.cp_eff(500.0) == m.cp(500.0)


def test_enthalpy_is_the_exact_integral_of_the_piecewise_linear_cp():
    m = material.Material.from_drama_json()
    assert m.enthalpy(293.0) == 0.0
    assert m.enthalpy(313.0) == pytest.approx(20.0 * (877.5 + 885.5) / 2.0)
    assert m.enthalpy(303.0) == pytest.approx(10.0 * (877.5 + 881.5) / 2.0)
    assert m.enthalpy(1850.0) == pytest.approx(m.enthalpy(850.0) + 1000.0 * 1131.6)          # constant c_p above the table
    assert m.enthalpy(273.0) == pytest.approx(-20.0 * 877.5)                                 # constant c_p below it
    T = np.array([250.0, 293.0, 400.0, 733.0, 850.0, 1500.0, 3000.0])
    assert np.all(np.diff(m.enthalpy(T)) > 0)


def test_temperature_from_enthalpy_round_trips():
    m = material.Material.from_drama_json()
    T = np.array([250.0, 293.0, 310.5, 512.0, 733.0, 850.0, 1234.5, 5000.0])
    assert m.temperature_from_enthalpy(m.enthalpy(T)) == pytest.approx(T, rel=1e-12)
    assert float(m.temperature_from_enthalpy(0.0)) == 293.0


def test_constant_property_material():
    m = material.Material("const", 1000.0, 0.5, [200.0, 2000.0], [900.0, 900.0], [200.0, 2000.0], [10.0, 10.0])
    assert m.enthalpy(393.0) == pytest.approx(900.0 * 100.0) and m.k(1500.0) == 10.0
    assert m.temperature_from_enthalpy(900.0 * 1000.0) == pytest.approx(1293.0)


# ---------------------------------------------------------------------------------------------------------------
# Step 3: latent heat, melting ranges, feed fraction, liquid properties

def test_material_names_and_liquid_properties():
    a, r, n = (material.Material.from_drama_json(k) for k in ("AA7075", "AA7075_range", "AA7075_nomelt"))
    assert a.melts and r.melts and not n.melts and n.latent_heat == 0.0 and n.liquid is None
    assert a.latent_heat == r.latent_heat == 400e3 and a.rho == r.rho == 2813.0 and a.emissivity == 0.4
    assert (a.T_solidus, a.T_liquidus) == (848.0, 852.0) and (r.T_solidus, r.T_liquidus) == (750.0, 908.0)
    assert a.liquid.rho == 2400.0 and a.liquid.mu == 1.3e-3 and a.liquid.sigma == 0.86
    assert np.allclose(a.k(a.T_k), n.k(a.T_k)) and np.allclose(a.cp(a.T_cp), n.cp(a.T_cp))    # DRAMA's tables verbatim
    assert a.T_cp.max() == 850.0 and n.T_cp.max() == 1e5


def test_liquid_and_feed_fractions():
    a, r = material.Material.from_drama_json("AA7075"), material.Material.from_drama_json("AA7075_range")
    assert np.allclose(a.liquid_fraction([840.0, 848.0, 850.0, 852.0, 900.0]), [0.0, 0.0, 0.5, 1.0, 1.0])
    assert np.allclose(a.feed_fraction([840.0, 848.0, 850.0, 852.0, 900.0]), [0.0, 0.0, 0.5, 1.0, 1.0])     # single T: feed = liquid fraction
    assert np.allclose(r.liquid_fraction([750.0, 829.0, 908.0]), [0.0, 0.5, 1.0])
    assert np.allclose(r.feed_fraction([900.0, 906.0, 908.0, 910.0]), [0.0, 0.0, 0.5, 1.0])              # feed only at the liquidus (+-2 K)
    assert r.T_feed == 910.0 and a.T_feed == 852.0
    assert a.cp_eff(850.0) == pytest.approx(a.cp(850.0) + 400e3 / 4.0) and r.cp_eff(800.0) == pytest.approx(r.cp(800.0) + 400e3 / 158.0)
    assert a.cp_eff(700.0) == a.cp(700.0) and r.cp_eff(950.0) == r.cp(950.0)


def test_enthalpy_jump_is_exact_and_invertible():
    for name in ("AA7075", "AA7075_range"):
        m = material.Material.from_drama_json(name)
        sensible = np.trapezoid(m.cp(np.linspace(m.T_solidus, m.T_liquidus, 2001)), np.linspace(m.T_solidus, m.T_liquidus, 2001))
        assert m.enthalpy(m.T_liquidus) - m.enthalpy(m.T_solidus) == pytest.approx(400e3 + sensible, rel=1e-9)
        T = np.linspace(200.0, 2000.0, 7201)
        h = m.enthalpy(T)
        assert np.all(np.diff(h) > 0.0) and np.abs(m.temperature_from_enthalpy(h) - T).max() < 1e-8
        assert m.h_liquid == pytest.approx(m.enthalpy(m.T_feed))
    n = material.Material.from_drama_json("AA7075_nomelt")
    assert n.enthalpy(1000.0) == pytest.approx(material.Material.from_drama_json("AA7075").enthalpy(1000.0) - 400e3, rel=1e-9)


def test_invalid_melting_data():
    with pytest.raises(ValueError):
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], -1.0, 850.0, 850.0)
    with pytest.raises(ValueError):
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 900.0, 850.0)


def test_liquid_and_mixed_enthalpy_invert_exactly(request):
    """The melt film is liquid: h_liquid(T) = h(T) + L_f (1 - f_l) carries the latent heat at every temperature, and a
    node holding a fraction w of film has a latent plateau only (1 - w) as tall. enthalpy_mixed and its inverse are
    each other's exact inverse for every w -- the solver's Newton update inverts the node's own mixture, and a
    mismatch there limit-cycled the iteration between 777 K and 864 K (Step 3)."""
    for name in ("AA7075", "AA7075_range"):
        mat = material.Material.from_drama_json(name)
        T = np.linspace(300.0, 1400.0, 2001)
        solid, liquid = T < mat.T_solidus, T > mat.T_liquidus
        assert mat.enthalpy_liquid(T)[solid] == pytest.approx(mat.enthalpy(T)[solid] + mat.latent_heat)
        assert mat.enthalpy_liquid(T)[liquid] == pytest.approx(mat.enthalpy(T)[liquid])
        assert mat.enthalpy_mixed(T, 0.0) == pytest.approx(mat.enthalpy(T)) and mat.enthalpy_mixed(T, 1.0) == pytest.approx(mat.enthalpy_liquid(T))
        assert mat.cp_mixed(T, 0.0) == pytest.approx(mat.cp_eff(T)) and mat.cp_mixed(T, 1.0) == pytest.approx(mat.cp(T))
        w = np.random.default_rng(0).random(T.size)
        for weights in (np.zeros_like(T), np.full_like(T, 0.3), np.ones_like(T), w):
            assert mat.temperature_from_enthalpy_mixed(mat.enthalpy_mixed(T, weights), weights) == pytest.approx(T, abs=1e-9)
        assert np.all(np.diff(mat.enthalpy_mixed(T, 0.4)) > 0.0)                     # monotone: the inverse is a function

# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-09-27: the Scheil material, a separate variant (AA7075_range stays as it is)

def scheil_fl(T, k=0.4, T_pure=933.0, T_liq=908.0):
    return ((T_pure - np.asarray(T, dtype=float)) / (T_pure - T_liq)) ** (1.0 / (k - 1.0))


def test_scheil_material_is_a_separate_variant():
    s, r = material.Material.from_drama_json("AA7075_scheil"), material.Material.from_drama_json("AA7075_range")
    assert s.scheil and not r.scheil and s.melts and s.partition_coefficient == 0.4 and s.T_pure == 933.0
    assert (r.T_solidus, r.T_liquidus) == (750.0, 908.0)                     # the range material is untouched
    assert (s.T_solidus, s.T_liquidus) == (748.0, 908.0)                     # foot of the eutectic ramp at the solidus
    assert s.latent_heat == 390e3 and r.latent_heat == 400e3 and s.rho == r.rho     # Step 4's latent heat (2026-09-28)
    assert s.liquid.sigma == 0.80 and r.liquid.sigma == 0.86                         # Step 4's surface tension (2026-09-28)
    assert s.liquid.rho == r.liquid.rho and s.liquid.mu == r.liquid.mu
    assert np.array_equal(s.T_cp, r.T_cp) and np.array_equal(s.cp_table, r.cp_table) and np.array_equal(s.k_table, r.k_table)
    assert s.T_feed == r.T_feed == 910.0 and s.h_liquid == pytest.approx(r.h_liquid - 10e3, rel=1e-12)   # fully liquid: L less


def test_scheil_liquid_fraction_and_latent_heat():
    s = material.Material.from_drama_json("AA7075_scheil")
    assert s.liquid_fraction(908.0) == 1.0 and s.liquid_fraction(747.9) == 0.0 and s.liquid_fraction(920.0) == 1.0
    assert scheil_fl(750.0) == pytest.approx(0.0362, abs=1e-4)               # the eutectic remainder at the solidus
    assert s.liquid_fraction(752.0) == pytest.approx(scheil_fl(752.0), rel=1e-12)      # top of the eutectic ramp
    assert s.liquid_fraction(750.0) == pytest.approx(0.5 * scheil_fl(752.0), rel=1e-12)
    T = np.arange(752.0, 908.0, 0.37)
    assert np.abs(s.liquid_fraction(T) - scheil_fl(T)).max() < 1e-3          # the 1 K table against the law
    assert s.liquid_fraction(895.1) == pytest.approx(0.5, abs=1e-3)          # half liquid 12.9 K below the liquidus
    def sens(a, b):                         # exact for the piecewise-linear c_p: the grid holds its kinks
        x = np.union1d(np.linspace(a, b, 2001), s.T_cp[(s.T_cp > a) & (s.T_cp < b)])
        return np.trapezoid(s.cp(x), x)
    latent = lambda a, b: s.enthalpy(b) - s.enthalpy(a) - sens(a, b)
    assert latent(748.0, 908.0) == pytest.approx(390e3, rel=1e-9)            # all of L_f, exactly
    assert latent(895.1, 908.0) / 390e3 == pytest.approx(0.5, abs=1e-3)      # half of it in the top 13 K ...
    r = material.Material.from_drama_json("AA7075_range")
    assert (r.enthalpy(908.0) - r.enthalpy(895.1) - sens(895.1, 908.0)) / 400e3 == pytest.approx(12.9 / 158.0, rel=1e-6)   # ... 8 % linearly


def test_scheil_enthalpy_is_exact_monotone_and_invertible():
    s = material.Material.from_drama_json("AA7075_scheil")
    T = np.linspace(200.0, 2000.0, 7201)
    h = s.enthalpy(T)
    assert np.all(np.diff(h) > 0.0) and np.abs(s.temperature_from_enthalpy(h) - T).max() < 1e-8
    Tm = np.array([749.3, 760.5, 820.5, 890.5, 907.5])                        # cp_eff is dh/dT inside every interval
    assert s.cp_eff(Tm) == pytest.approx((s.enthalpy(Tm + 1e-4) - s.enthalpy(Tm - 1e-4)) / 2e-4, rel=1e-6)
    assert s.cp_eff(700.0) == s.cp(700.0) and s.cp_eff(950.0) == s.cp(950.0)
    solid, liquid = T < s.T_solidus, T > s.T_liquidus
    assert s.enthalpy_liquid(T)[solid] == pytest.approx(s.enthalpy(T)[solid] + s.latent_heat)
    assert s.enthalpy_liquid(T)[liquid] == pytest.approx(s.enthalpy(T)[liquid])
    assert s.cp_mixed(T, 0.0) == pytest.approx(s.cp_eff(T)) and s.cp_mixed(T, 1.0) == pytest.approx(s.cp(T))
    w = np.random.default_rng(0).random(T.size)
    for weights in (np.zeros_like(T), np.full_like(T, 0.3), np.ones_like(T), w):
        assert s.temperature_from_enthalpy_mixed(s.enthalpy_mixed(T, weights), weights) == pytest.approx(T, abs=1e-9)
    assert np.all(np.diff(s.enthalpy_mixed(T, 0.4)) > 0.0)


def test_invalid_scheil_data():
    args = ("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 750.0, 908.0, None)
    for k, T_pure in ((0.0, 933.0), (1.0, 933.0), (0.4, 900.0)):
        with pytest.raises(ValueError):
            material.Material(*args, k, T_pure)
    with pytest.raises(ValueError):                                          # range narrower than the eutectic ramp
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 850.0, 852.0, None, 0.4, 933.0)

# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-09-28: AA7075-empiricaldata, AA7075 with only the empirical surface tension and heat capacity

def test_empirical_material_is_aa7075_with_two_changes():
    import json
    da = json.load(open(material.MATERIAL_NAMES["AA7075"]))
    de = json.load(open(material.MATERIAL_NAMES["AA7075-empiricaldata"]))
    assert {k for k in da if da[k] != de[k]} == {"_comment", "name", "liquid", "specificHeatCapacity"}
    assert {k for k in da["liquid"] if da["liquid"][k] != de["liquid"][k]} == {"surfaceTension", "_sources"}
    assert de["specificHeatCapacity"][:len(da["specificHeatCapacity"])] == da["specificHeatCapacity"]
    assert de["specificHeatCapacity"][len(da["specificHeatCapacity"]):] == [[900.0, 1131.6], [1000.0, 1131.6], [1100.0, 1131.6], [1200.0, 1131.6]]
    e, a = material.Material.from_drama_json("AA7075-empiricaldata"), material.Material.from_drama_json("AA7075")
    assert e.name == "AA7075-empiricaldata" and e.liquid.sigma == 0.80 and a.liquid.sigma == 0.86
    assert (e.latent_heat, e.T_solidus, e.T_liquidus, e.rho, e.emissivity) == (a.latent_heat, a.T_solidus, a.T_liquidus, a.rho, a.emissivity)
    T = np.linspace(200.0, 2000.0, 7201)
    assert np.array_equal(e.cp(T), a.cp(T)) and np.array_equal(e.k(T), a.k(T))   # the held value made explicit: same c_p
    assert e.enthalpy(T) == pytest.approx(a.enthalpy(T), rel=1e-12) and np.abs(e.temperature_from_enthalpy(e.enthalpy(T)) - T).max() < 1e-8


# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-10-06: the rigid substrate -- the temperature at which the material is half liquid

def test_the_rigid_temperature_is_where_half_the_material_is_liquid():
    """T_rigid separates the coherent mush, which carries load, from the slurry, which flows: 50 % liquid, where Chen
    et al.'s (2016) semi-solid strength law ends and Li et al.'s (2014) slurry viscosity data begin. It follows each
    material's own liquid-fraction law."""
    assert material.RIGID_LIQUID_FRACTION == 0.5
    r = material.Material.from_drama_json("AA7075_range")
    assert r.T_rigid == pytest.approx(750.0 + 0.5 * 158.0, abs=1e-9)                     # 829 K on the linear range
    s = material.Material.from_drama_json("AA7075_scheil")
    assert s.T_rigid == pytest.approx(895.1, abs=0.05)                                    # Scheil: 13 K below the liquidus
    a = material.Material.from_drama_json("AA7075")
    assert a.T_rigid == pytest.approx(850.0, abs=1e-9)                                    # single temperature: mid-ramp
    for m in (r, s, a):
        assert float(m.liquid_fraction(m.T_rigid)) == pytest.approx(0.5, abs=1e-9)
    assert np.isinf(material.Material.from_drama_json("AA7075_nomelt").T_rigid)          # never melts: never slurry
