"""Tests of the level set's property table (stage 3 Task 4): the shared Material's h, cp_eff, k and liquid fraction
tabulated on a 0.01 K grid and interpolated linearly -- a numerical evaluation device, not new physics -- with an
enthalpy inversion consistent with the table's own h."""
import numpy as np
import pytest

from reentry_model.levelset.properties import PropertyTable
from reentry_model.shared.material import Material

from materials import ConstMaterial, LinearMeltMaterial


def relative(a, b):
    return np.abs(a - b) / np.maximum(np.abs(b), np.finfo(float).tiny)


def test_the_table_matches_the_shared_scheil_material():
    """At 1000 random temperatures over 300-1200 K (through the Scheil range, 748-908 K) the table gives the Material's
    h and liquid fraction within 1e-6 relative and its cp_eff and k within 1e-4 relative. Measured 2026-11-12 (seed
    20261112): h 5.7e-10, f_l 4.6e-15 (absolute; exactly zero where the Material's is), cp_eff 6.5e-11
    (its steps on the grid's nodes kept by the one-sided ends), k 8.8e-16."""
    mat = Material.from_drama_json("AA7075_scheil")
    table = PropertyTable(mat)
    T = np.random.default_rng(20261112).uniform(300.0, 1200.0, 1000)
    assert np.all(relative(table.enthalpy(T), mat.enthalpy(T)) <= 1e-6)
    f, f_mat = table.liquid_fraction(T), mat.liquid_fraction(T)
    assert np.all(np.abs(f - f_mat) <= 1e-6 * np.abs(f_mat))
    assert np.all(relative(table.cp_eff(T), mat.cp_eff(T)) <= 1e-4)
    assert np.all(relative(table.k(T), mat.k(T)) <= 1e-4)
    mush = (T > mat.T_solidus) & (T < mat.T_liquidus)
    assert mush.sum() > 100 and (f[mush] > 0.0).all()                  # the draw did cross the Scheil range


def test_the_enthalpy_inversion_round_trips_on_the_tables_own_h():
    mat = Material.from_drama_json("AA7075_scheil")
    table = PropertyTable(mat)
    rng = np.random.default_rng(5)
    T = np.concatenate([rng.uniform(150.0, 3000.0, 2000), rng.uniform(740.0, 915.0, 2000),
                        150.0 + 0.01 * rng.integers(0, 285001, 200), [150.0, 3000.0, 895.1, 908.0]])
    assert np.max(np.abs(table.temperature_from_enthalpy(table.enthalpy(T)) - T)) <= 1e-9
    # and the Material's own enthalpies come back to their temperatures within the table's interpolation
    assert np.max(np.abs(table.temperature_from_enthalpy(mat.enthalpy(T)) - T)) <= 1e-6


def test_the_table_is_vectorised_and_its_slopes_are_its_own():
    table = PropertyTable(LinearMeltMaterial())
    T = np.array([[300.0, 879.995, 894.0], [908.0, 950.0, 2999.0]])
    for f in (table.enthalpy, table.cp_eff, table.k, table.liquid_fraction, table.dh_dT, table.dk_dT):
        assert f(T).shape == T.shape
    h, dh, k, dk = table.conduction(T)
    assert np.array_equal(h, table.enthalpy(T)) and np.array_equal(k, table.k(T))
    # the slope is the tabulated h's own across the grid interval each temperature falls in
    i = np.floor((T - 150.0) / 0.01)
    lo, hi = 150.0 + 0.01 * i, 150.0 + 0.01 * (i + 1)
    assert dh == pytest.approx((table.enthalpy(hi) - table.enthalpy(lo)) / (hi - lo), rel=1e-9)
    assert np.array_equal(dk, np.zeros_like(T))                          # the test material's k is constant
    assert table.temperature_from_enthalpy(h).shape == T.shape


def test_a_table_needs_a_rising_enthalpy():
    class Flat(ConstMaterial):
        def enthalpy(self, T):
            return np.zeros_like(np.asarray(T, float))
    with pytest.raises(ValueError, match="rise"):
        PropertyTable(Flat())
