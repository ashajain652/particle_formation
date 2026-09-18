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
