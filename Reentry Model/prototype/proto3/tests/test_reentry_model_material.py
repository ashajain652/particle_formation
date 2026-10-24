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
