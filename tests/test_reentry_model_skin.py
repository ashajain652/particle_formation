"""The layered skins (layered-skin design section 6, sub-plan 20): cells with variable counts, split and merge, the
minimum mass, heat through the depth with liquid removed as it forms, and heat sideways between columns."""
import numpy as np
import pytest

from reentry_model import material as material_mod
from reentry_model import skin


def scheil():
    return material_mod.Material.from_drama_json(material_mod.MATERIAL_NAMES["AA7075_scheil"])


def sharp(k=100.0, rho=2800.0, cp=1000.0, latent=4.0e5, T_melt=900.0):
    """Constant properties, a melting point (with the material's +-2 K ramp) and liquid as dense as solid: the
    material of the steady-ablation solution."""
    return material_mod.Material("sharp", rho, 0.0, np.array([0.0, 5000.0]), np.array([cp, cp]), np.array([0.0, 5000.0]),
                                 np.array([k, k]), latent, T_melt, T_melt, material_mod.LiquidProperties(rho, 1e-3, 0.8))


def filled(mat, n=3, area=4e-6, T=700.0, settings=None):
    s = settings or skin.SkinSettings()
    sk = skin.Skins(mat, np.full(n, area), s)
    sk.fill(skin.minimum_mass(mat, sk.area, s), T)
    return sk


# -- settings and cells (Task 2) ---------------------------------------------------------------------------------------
def test_settings_reject_impossible_values():
    for bad in ({"min_thickness": 0.0}, {"cell": -1e-6}, {"substep": 0.0}, {"growth": 0.9}, {"cell": 0.3e-3}):
        with pytest.raises(ValueError):
            skin.SkinSettings(**bad)


def test_a_filled_skin_holds_its_mass_in_fine_cells():
    mat = scheil()
    sk = filled(mat)
    assert np.all(sk.count == 40)                                   # 0.4 mm of 10 um cells
    assert sk.thickness() == pytest.approx(np.full(3, 0.4e-3), rel=1e-12)
    assert sk.mass() == pytest.approx(3 * mat.rho * 4e-6 * 0.4e-3, rel=1e-12)
    assert sk.energy() == pytest.approx(sk.mass() * float(mat.enthalpy(700.0)), rel=1e-12)
    assert np.allclose(sk.T[sk.valid], 700.0)


def test_split_and_merge_conserve_mass_and_enthalpy_and_keep_cells_near_their_targets():
    mat = scheil()
    sk = filled(mat, n=4)
    rng = np.random.default_rng(1)
    v = sk.valid
    sk.m = np.where(v, sk.m * rng.uniform(0.1, 3.0, sk.m.shape), 0.0)       # cells of every size
    sk.H = sk.m * mat.enthalpy(rng.uniform(600.0, 890.0, sk.m.shape))
    sk._install(sk.m, sk.H, sk.count)
    M0, E0 = sk.mass(), sk.energy()
    sk.regrid()
    assert sk.mass() == pytest.approx(M0, rel=1e-14) and sk.energy() == pytest.approx(E0, rel=1e-14)
    dz, t = sk.dz(), sk.targets()
    v = sk.valid
    assert np.all(dz[v] >= skin.MERGE_BELOW * t[v]) and np.all(dz[v] <= skin.SPLIT_ABOVE * t[v])


def test_the_top_cells_are_untouched_by_a_split_or_merge_deeper_down():
    mat = scheil()
    sk = filled(mat, n=2)
    top_m, top_H = sk.m[:, :10].copy(), sk.H[:, :10].copy()
    sk.add_bottom(np.array([1e-7, 0.0]), float(mat.enthalpy(700.0)))       # a bottom cell 1.9 x its target: split
    sk.take_bottom(np.array([0.0, 0.7 * sk.m[1, 0]]))                       # a thin bottom cell: merged
    sk.regrid()
    assert np.array_equal(sk.m[:, :10], top_m) and np.array_equal(sk.H[:, :10], top_H)
    assert sk.count[0] > 40 and sk.count[1] == 39


def test_cells_below_the_fine_zone_grow_by_the_growth_factor():
    mat = scheil()
    sk = filled(mat, n=1)
    sk.add_bottom(np.array([mat.rho * 4e-6 * 2.6e-3]), float(mat.enthalpy(700.0)))   # 3 mm in all
    sk.regrid()
    dz = sk.dz()[0, :sk.count[0]]
    assert sk.thickness()[0] == pytest.approx(3.0e-3, rel=1e-12)
    # the old bottom cell took the 2.6 mm and was halved until every piece is within its bounds: 39 fine cells above
    # it, then pieces doubling in steps as the target grows (75 cells in all, measured while planning; a continuous
    # 10 % growth would need 40 + 10 ln(1 + 0.1 x 2600 / 10) = 73)
    assert np.allclose(dz[:39], 10e-6)
    assert 60 <= sk.count[0] <= 90
    assert np.all(np.diff(dz) >= -1e-15)                             # never a thick cell above a thinner one


def test_taking_and_exchanging_at_the_bottom_and_top_are_exact():
    mat = scheil()
    sk = filled(mat)
    M0, E0 = sk.mass(), sk.energy()
    dm = np.array([0.0, 0.25, 1.5]) * sk.m[0, 0]
    H = sk.take_bottom(dm)
    assert sk.mass() == pytest.approx(M0 - dm.sum(), rel=1e-14) and sk.energy() == pytest.approx(E0 - H.sum(), rel=1e-14)
    sk.exchange_bottom(np.array([1e-9, -2e-9, 0.0]), np.array([1e-3, -3e-3, 5e-4]))
    assert sk.mass() == pytest.approx(M0 - dm.sum() - 1e-9, rel=1e-14)
    assert sk.energy() == pytest.approx(E0 - H.sum() - 1.5e-3, rel=1e-14)
    m_top, H_top = sk.take_top(np.array([2, 0, 1]))
    assert np.all(sk.count == np.array([38, 40, 38]))                  # the third had lost a cell and a half
    assert m_top[1] == 0.0 and m_top[0] == pytest.approx(2 * sk.m[0, 0], rel=1e-12)
    with pytest.raises(ValueError):
        sk.take_bottom(np.full(3, 1.0))


def test_depths_read_from_the_cells():
    mat = sharp()                                                    # liquid as dense as solid: the cells keep 10 um
    sk = filled(mat, n=1)
    z_mid = np.arange(40) * 10e-6 + 5e-6
    sk.T[0, :40] = 920.0 - 1e5 * z_mid                               # a straight profile, 10 K per 100 um
    sk.H = sk.m * mat.enthalpy(sk.T)
    assert sk.depth_above(900.0)[0] == pytest.approx(200e-6, abs=1e-9)
    assert sk.depth_above(1000.0)[0] == 0.0 and sk.depth_above(0.0)[0] == pytest.approx(0.4e-3)


# -- the minimum mass (Task 3) -----------------------------------------------------------------------------------------
def test_the_minimum_is_400_um_of_solid_unless_girins_layer_needs_more():
    mat = scheil()
    s = skin.SkinSettings()
    area = np.full(4, 4e-6)
    solid = mat.rho * 4e-6 * 0.4e-3
    # 400 um of solid is the liquid of 469 um (2813 / 2400 x 400): delta_m + 50 um decides above 419 um
    m = skin.minimum_mass(mat, area, s, np.array([np.nan, 300e-6, 419e-6, 500e-6]))
    assert m[0] == solid and m[1] == solid
    assert m[2] == pytest.approx(solid, rel=2e-3)
    assert m[3] == pytest.approx(mat.liquid.rho * 4e-6 * 550e-6, rel=1e-12) and m[3] > solid
    assert np.all(skin.minimum_mass(mat, area, s) == solid)


