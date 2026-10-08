"""reentry_model.skin: the one-dimensional skins of the 2026-10-07 design (Step 3 sub-plan 18) -- storage and rezoning,
the top node, depths, hand-over, the implicit sub-step against exact solutions, feed and draw, and the macro step's
trial and real passes."""
import math

import numpy as np
import pytest
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.special import erf, erfc

from reentry_model import material, skin


def constant_material(latent=0.0, T_melt=float("inf")):
    """rho 2800, c_p 1000, k 150, no radiation; a +-2 K melting range at T_melt when latent > 0."""
    return material.Material("const", 2800.0, 0.0, np.array([200.0, 3000.0]), np.array([1000.0, 1000.0]),
                             np.array([200.0, 3000.0]), np.array([150.0, 150.0]), latent,
                             T_melt - 2.0 if latent else float("inf"), T_melt + 2.0 if latent else float("inf"),
                             material.LiquidProperties(2400.0, 1.3e-3, 0.86))


def field(mat=None, n=3, nc=8, thickness=0.4e-3, T=None):
    mat = mat or material.Material.from_drama_json("AA7075_scheil")
    f = skin.SkinField(mat, skin.SkinSettings(thickness=thickness, cell=thickness / nc))
    T = np.linspace(900.0, 820.0, nc)[None, :].repeat(n, axis=0) if T is None else T
    f.add(np.arange(10, 10 + n), np.full(n, 2.0e-6), T, np.full(n, thickness))
    return f


def node_total(f, row, L):
    """The skin's metal enthalpy plus the liquid riding its top node [J]."""
    return float(f.H[row].sum() + L * f.material.enthalpy_liquid(f.T[row, 0]))


# -- storage, rezoning, the top node, depths, hand-over (Task 2) ----------------------------------------------------
def test_settings_are_validated():
    assert skin.SkinSettings().n_cells == 40 and skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"]).n_cells == 20
    for bad in (dict(thickness=0.0), dict(cell=-1e-6), dict(substep=0.0), dict(thickness=1e-4, cell=6e-5)):
        with pytest.raises(ValueError):
            skin.SkinSettings(**bad)


def test_added_skins_hold_their_slab_and_its_enthalpy():
    f = field()
    mat = f.material
    assert f.n == 3 and f.m.shape == (3, 8)
    assert f.mass() == pytest.approx(3 * mat.rho * 2.0e-6 * 0.4e-3, rel=1e-14)
    assert f.energy() == pytest.approx(float((f.m * mat.enthalpy(f.T)).sum()), rel=1e-14)
    assert np.allclose(f.thickness(), 0.4e-3, rtol=1e-14)
    assert list(f.rows_of(np.array([11, 99, 10]))) == [1, -1, 0]


def test_remap_conserves_mass_and_enthalpy_and_keeps_the_order():
    f = field()
    f.m[0] *= np.linspace(0.5, 1.5, 8)                  # unequal cells
    f.H[0] = f.m[0] * f.material.enthalpy(f.T[0])
    M = f.m[0].sum()
    L = np.array([3e-9])
    before = node_total(f, 0, L[0])
    f.remap(np.array([0]), L)
    assert f.m[0].sum() == pytest.approx(M, rel=1e-14) and np.allclose(f.m[0], M / 8, rtol=1e-14)
    assert node_total(f, 0, L[0]) == pytest.approx(before, rel=1e-12)
    assert np.all(np.diff(f.T[0, 1:]) <= 1e-9)          # still hottest at the top


def test_the_top_node_keeps_its_energy_when_liquid_is_booked_or_refreezes():
    f = field()
    rows, L = np.array([1]), np.array([5e-9])
    E0 = f.top_energy(rows, L)
    f.book_top(rows, np.array([1e-4]), L)
    assert f.top_energy(rows, L) == pytest.approx(E0 + 1e-4, rel=1e-12)
    before, T1, M = node_total(f, 1, 5e-9), f.T[1, 0], f.m[1].sum()
    f.refreeze_top(rows, np.array([4e-9]), np.array([1e-9]))     # most of the liquid freezes onto cell 0
    assert f.m[1].sum() == pytest.approx(M + 4e-9, rel=1e-12)
    assert node_total(f, 1, 1e-9) == pytest.approx(before, rel=1e-12)
    assert f.T[1, 0] >= T1 - 1e-9                       # freezing releases latent heat


def test_depth_above_interpolates_between_cell_centres():
    mat = constant_material()
    nc, s = 10, 1.0e-3
    T = np.linspace(1000.0, 900.0, nc)[None, :]         # linear: T = 1000 - 100 (z - dz/2)/(s - dz)
    f = field(mat, n=1, nc=nc, thickness=s, T=T)
    dz = s / nc
    depth, through = f.depth_above(950.0)
    assert depth[0] == pytest.approx(dz / 2 + 0.5 * (s - dz), rel=1e-12) and not through[0]
    depth, through = f.depth_above(800.0)
    assert depth[0] == pytest.approx(s) and through[0]
    depth, _ = f.depth_above(1200.0)
    assert depth[0] == 0.0


def test_hand_over_passes_skins_layer_by_layer():
    f = field()
    M, E = f.mass(), f.energy()
    areas = {100: 1.5e-6, 101: 2.5e-6}
    tgt = np.array([[100, 101], [100, -1]])
    w = np.array([[0.4, 0.6], [1.0, 0.0]])
    _, handed = f.hand_over(np.array([0, 1]), tgt, w, lambda faces: np.array([areas[int(x)] for x in faces]))
    assert handed == pytest.approx(float(f.m[[0, 1]].sum()), rel=1e-14)
    f.keep(np.setdiff1d(np.arange(f.n), [0, 1]))
    assert f.mass() == pytest.approx(M, rel=1e-13) and f.energy() == pytest.approx(E, rel=1e-12)
    r = f.rows_of(np.array([100, 101]))
    assert (r >= 0).all() and np.all(np.diff(f.T[r[0]]) <= 1e-9)          # still hottest at the top
    assert np.allclose(f.m[r[0]], f.m[r[0]].mean(), rtol=1e-12)          # equal cells


def test_hand_over_falls_back_and_drops_empty_skins():
    f = field()
    f.m[2] = 0.0
    f.H[2] = 0.0
    _, handed = f.hand_over(np.array([2]), np.array([[105, -1]]), np.array([[1.0, 0.0]]), lambda faces: np.full(len(faces), 1e-6))
    f.keep(np.flatnonzero(f.m.sum(axis=1) > 0.0))
    assert f.rows_of(np.array([105]))[0] == -1 and handed == 0.0 and f.n == 2


# -- the implicit sub-step (Task 3) ---------------------------------------------------------------------------------
def test_thomas_matches_a_banded_solver():
    rng = np.random.default_rng(0)
    n, k = 5, 9
    sub, sup = -rng.random((n, k - 1)), -rng.random((n, k - 1))
    diag = 3.0 + rng.random((n, k))
    rhs = rng.random((n, k))
    x = skin.thomas(sub, diag, sup, rhs)
    for i in range(n):
        ab = np.zeros((3, k))
        ab[0, 1:], ab[1], ab[2, :-1] = sup[i], diag[i], sub[i]
        assert np.allclose(x[i], solve_banded((1, 1), ab, rhs[i]), rtol=1e-13, atol=0.0)


def test_constant_flux_on_a_semi_infinite_solid():
    """Carslaw and Jaeger 2.9: T - T_i = (2 q / k) sqrt(alpha t) ierfc(z / (2 sqrt(alpha t))), at the top cell's centre."""
    mat = constant_material()
    s, nc = 4.0e-3, 400
    f = field(mat, n=1, nc=nc, thickness=s, T=np.full((1, nc), 300.0))
    q, t_end, dt = 2.0e6, 0.05, 2.5e-4
    for _ in range(int(round(t_end / dt))):
        f.conduct(dt, np.array([q]), 0.0, np.array([0.0]), np.array([300.0]), np.array([0.0]))
    alpha = 150.0 / (2800.0 * 1000.0)
    x = 0.5 * s / nc / (2.0 * math.sqrt(alpha * t_end))
    ierfc = math.exp(-x * x) / math.sqrt(math.pi) - x * math.erfc(x)
    exact = 300.0 + 2.0 * q / 150.0 * math.sqrt(alpha * t_end) * ierfc
    assert f.T[0, 0] == pytest.approx(exact, abs=0.01 * (exact - 300.0))


def test_neumanns_melting_front():
    """Two-phase Neumann problem, equal properties: the front at 2 lambda sqrt(alpha t), lambda from
    St/(exp(l^2) erf l) - St/(exp(l^2) erfc l) = l sqrt(pi) with St = c (T_s - T_m)/L = c (T_m - T_i)/L; the front is
    where the material is half liquid."""
    mat = constant_material(latent=4.0e5, T_melt=900.0)
    s, nc = 8.0e-3, 400
    f = field(mat, n=1, nc=nc, thickness=s, T=np.full((1, nc), 800.0))
    t_end, dt = 0.05, 2.5e-4
    for _ in range(int(round(t_end / dt))):
        f.conduct(dt, np.array([0.0]), 0.0, np.array([0.0]), np.array([800.0]), np.array([0.0]), T_surface=1000.0)
    St = 1000.0 * 100.0 / 4.0e5
    lam = brentq(lambda l: St / (math.exp(l * l) * erf(l)) - St / (math.exp(l * l) * erfc(l)) - l * math.sqrt(math.pi), 1e-3, 3.0)
    alpha = 150.0 / (2800.0 * 1000.0)
    front, _ = f.depth_above(900.0)
    assert front[0] == pytest.approx(2.0 * lam * math.sqrt(alpha * t_end), rel=0.03)


def test_a_substep_closes_the_skin_books_with_the_conducted_base_flux():
    f = field()
    L = np.array([0.0, 2e-9, 0.0])
    E0 = sum(node_total(f, i, L[i]) for i in range(3))
    E_top, E_rad, E_base, _ = f.conduct(0.01, np.full(3, 1.5e6), 0.0, np.full(3, 0.3), np.full(3, 800.0), L)
    E1 = sum(node_total(f, i, L[i]) for i in range(3))
    assert E1 - E0 == pytest.approx(float((E_top - E_rad - E_base).sum()), rel=1e-12)
    dz = f.thickness() / f.nc
    conducted = f.area * f.material.k(f.T[:, -1]) / (0.5 * dz) * (f.T[:, -1] - 800.0) * 0.01
    assert np.allclose(E_base, conducted, rtol=1e-6)
