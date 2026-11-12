"""Tests of the level set's skins (design section 20.6, stage 3): the cells' storage, their grading below the fine zone,
the regrid's exact books, the curved cross-section A(z) = A0 (1 - k1 z)(1 - k2 z), and the shared mixture density."""
import numpy as np
import pytest

from reentry_model.levelset.skin import SKIN_PRESETS, Skins, SkinSettings
from reentry_model.shared.material import Material

from materials import ConstMaterial, LinearMeltMaterial

AREA = 4e-6                              # m2, a 2 mm square column


def column(mat, thick, T, kappa=(0.0, 0.0), settings=None, area=AREA):
    """One flat-estimate column of `thick` metal at T, filled in cells of about the fine size."""
    sk = Skins(mat, np.array([area]), np.array([kappa], dtype=float), settings or SkinSettings())
    sk.fill(mat.rho * area * thick, T)
    return sk


def test_settings_defaults_and_presets():
    s = SkinSettings()
    assert (s.cell, s.fine_zone, s.growth, s.stub, s.substep, s.swell) == (10e-6, 0.4e-3, 1.10, 0.8e-3, 0.010, True)
    assert SKIN_PRESETS == {"production": 10e-6, "dev": 20e-6}
    with pytest.raises(ValueError):
        SkinSettings(growth=0.9)


def test_cells_grade_ten_percent_below_the_fine_zone():
    mat = Material.from_drama_json("AA7075_scheil")
    sk = column(mat, 3e-3, 600.0)
    assert sk.count[0] == 300
    sk.regrid()
    n = sk.count[0]
    dz, t = sk.dz()[0, :n], sk.targets()[0, :n]
    z_top = np.cumsum(dz) - dz
    assert np.sum(z_top + dz <= 0.4e-3 * (1 + 1e-9)) == 40                    # the fine zone keeps its 40 cells
    assert np.allclose(dz[:40], 10e-6, rtol=1e-12)
    deep = z_top >= 0.4e-3 * (1 - 1e-12)
    # each deeper target is 10 % of the cell above it thicker than that cell's: the continuous form of "10 % thicker"
    i = np.flatnonzero(deep[:-1] & deep[1:])
    assert len(i) > 5
    assert t[i + 1] == pytest.approx(t[i] + 0.10 * dz[i], rel=1e-12)
    assert np.all((dz >= 0.5 * t) & (dz <= 1.5 * t))                         # every cell within its regrid band
    assert sk.thickness()[0] == pytest.approx(3e-3, rel=1e-12)
    assert n - 40 < 60                                                        # graded: far fewer than 260 deep cells


def test_regrid_conserves_mass_enthalpy_and_momentum():
    rng = np.random.default_rng(7)
    mat = LinearMeltMaterial()
    n, K = 5, 120
    sk = Skins(mat, rng.uniform(1e-6, 9e-6, n), rng.uniform(-20.0, 30.0, (n, 2)), SkinSettings())
    count = rng.integers(1, K + 1, n)
    m = rng.uniform(0.05, 3.0, (n, K)) * mat.rho * sk.area[:, None] * 10e-6
    T = rng.uniform(600.0, 1000.0, (n, K))
    sk.install(m, m * mat.enthalpy(T), count, rng.normal(0.0, 1.0, (n, K, 3)))
    M0, E0, P0 = sk.totals()
    scale = (sk.m[..., None] * np.abs(sk.U)).sum(axis=1)
    before = sk.count.copy()
    sk.regrid()
    assert not np.array_equal(sk.count, before)
    M1, E1, P1 = sk.totals()
    assert M1 == pytest.approx(M0, rel=1e-14)
    assert E1 == pytest.approx(E0, rel=1e-14)
    assert np.all(np.abs(P1 - P0) <= 1e-14 * scale)
    dz, t = sk.dz(), sk.targets()
    v = sk.valid
    assert np.all(dz[v] <= 1.5 * t[v])
    assert np.all(((dz >= 0.5 * t) | (sk.count[:, None] == 1))[v])


def test_the_curved_cross_section_holds_a_spherical_shell():
    R = 25e-3
    mat = ConstMaterial()
    for t_flat in (0.3e-3, 3e-3, 8e-3):                      # R/3 = 8.33 mm fills the sphere
        sk = column(mat, t_flat, 500.0, kappa=(1 / R, 1 / R))
        sk.regrid()
        t = sk.thickness()[0]
        vol = sk.totals()[0][0] / mat.rho
        shell = (4 * np.pi / 3) * (R ** 3 - (R - t) ** 3) / (4 * np.pi * R ** 2) * AREA   # per 4 pi R^2 of top area
        assert shell == pytest.approx(vol, rel=1e-12)
        assert t > t_flat                                                     # a convex column narrows with depth
    # a concave column widens with depth: thinner than flat, and each cell of one mass is thinner than the one above
    sk = column(mat, 3e-3, 500.0, kappa=(-1 / R, -1 / R))
    assert sk.thickness()[0] < 3e-3
    dz = sk.dz()[0, :sk.count[0]]
    assert np.all(np.diff(dz) < 0.0)
    A = sk.cross_section(np.array([[0.0, 1e-3, 2e-3]]))[0]
    assert np.all(np.diff(A) > 0.0) and A[0] == AREA
    # a skin deeper than its radius of curvature has no positive cross-section: refused, not given negative areas
    with pytest.raises(ValueError, match="radius of curvature"):
        column(mat, 9e-3, 500.0, kappa=(1 / 5e-3, 0.0)).dz()


def test_a_flat_skin_has_the_flat_thickness():
    mat = ConstMaterial()
    m = mat.rho * AREA * 1.7e-3
    sk = Skins(mat, np.array([AREA]), np.zeros((1, 2)), SkinSettings())
    sk.fill(m, 700.0)
    assert sk.thickness()[0] == pytest.approx(m / (mat.rho * AREA), rel=1e-15)
    assert np.all(sk.dz()[0, :sk.count[0]] == sk.m[0, :sk.count[0]] / (mat.rho * AREA))


def test_density_is_the_shared_mixture_law():
    mat = Material.from_drama_json("AA7075_scheil")
    T = np.array([600.0, 800.0, 895.1, 905.0, 1000.0])
    f = mat.liquid_fraction(T)
    sk = Skins(mat, np.array([AREA]), np.zeros((1, 2)), SkinSettings())
    assert np.array_equal(sk.density(T), 1.0 / ((1.0 - f) / mat.rho + f / mat.liquid.rho))
    dry = Skins(mat, np.array([AREA]), np.zeros((1, 2)), SkinSettings(swell=False))
    assert np.array_equal(dry.density(T), np.full(5, mat.rho))


def test_temperature_follows_the_enthalpy_and_phase_is_not_stored():
    mat = Material.from_drama_json("AA7075_scheil")
    sk = column(mat, 0.2e-3, 850.0)
    assert sk.T()[0, :sk.count[0]] == pytest.approx(850.0, abs=1e-9)
    sk.H[0, 0] = sk.m[0, 0] * mat.enthalpy(950.0)
    assert sk.T()[0, 0] == pytest.approx(950.0, abs=1e-9)


def test_add_bottom_take_top_and_depth_above_keep_the_books():
    mat = LinearMeltMaterial()
    sk = Skins(mat, np.full(3, AREA), np.zeros((3, 2)), SkinSettings())
    sk.fill(mat.rho * AREA * 0.4e-3, np.array([950.0, 600.0, 700.0]))
    sk.U[:] = np.array([1.0, 2.0, 0.0])
    M0, E0, P0 = sk.totals()
    c0 = sk.count[0]
    dm, dH = mat.rho * AREA * 50e-6, mat.rho * AREA * 50e-6 * mat.enthalpy(650.0)
    sk.add_bottom(np.array([1, 1]), np.array([dm, dm]), np.array([dH, dH]))   # a repeated index adds twice
    M1, E1, P1 = sk.totals()
    assert M1 == pytest.approx(M0 + np.array([0, 2 * dm, 0]), rel=1e-14)
    assert E1 == pytest.approx(E0 + np.array([0, 2 * dH, 0]), rel=1e-14)
    assert P1 == pytest.approx(P0, rel=1e-14)                                 # metal from below arrives at rest
    mass, H, P = sk.take_top(np.array([0, 2]), np.array([3, 0]))
    assert mass[0] == pytest.approx(3 * M0[0] / c0, rel=1e-12) and mass[1] == 0.0
    M2, E2, P2 = sk.totals()
    assert M2[0] + mass[0] == pytest.approx(M1[0], rel=1e-14) and E2[0] + H[0] == pytest.approx(E1[0], rel=1e-14)
    assert P2[0] + P[0] == pytest.approx(P1[0], rel=1e-14)
    assert sk.count[0] == c0 - 3
    # depth above a level: 0 when the top is below it, the whole skin when all is above
    assert sk.depth_above(1000.0)[0] == 0.0
    assert sk.depth_above(900.0)[0] == pytest.approx(sk.thickness()[0], rel=1e-14)
    sk.H[2, 10:] = sk.m[2, 10:] * mat.enthalpy(500.0)                         # ten cells at 700 K over 500 K
    assert sk.depth_above(600.0)[2] == pytest.approx(sk.dz()[2, :10].sum(), rel=1e-12)
    # an emptied column gets a new cell
    sk.take_top(np.array([2]), np.array([sk.count[2]]))
    assert sk.count[2] == 0
    sk.add_bottom(np.array([2]), np.array([dm]), np.array([dH]))
    assert sk.count[2] == 1 and sk.T()[2, 0] == pytest.approx(650.0, abs=1e-9)
