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


# -- the skin's heat (Task 3) ------------------------------------------------------------------------------------------
from reentry_model.shared.thermal import SIGMA_SB


def layered(mat, layers, area=AREA, kappa=(0.0, 0.0), settings=None, cell=10e-6):
    """One column of layers [(thickness [m], T [K]), ...] from the top down, in cells of `cell` by the solid density."""
    m, H = [], []
    for thick, T in layers:
        k = max(int(round(thick / cell)), 1)
        mc = mat.rho * area * thick / k
        m += [mc] * k
        H += [mc * float(mat.enthalpy(T))] * k
    sk = Skins(mat, np.array([area]), np.array([kappa], dtype=float), settings or SkinSettings(swell=False))
    sk.install(np.array([m]), np.array([H]), np.array([len(m)]))
    return sk


def test_steady_ablation_recedes_at_the_energy_limited_speed_with_the_exponential_profile():
    # a sharp-melting test material (1 K melting range just below T_feed = T_liquidus), constant k, cp, rho; k chosen so
    # the decay length alpha/v (1.0 mm) is short against the 5 mm column (exp(-vL/alpha) = 0.7 %) and the surface
    # gradient (0.6 K per micrometre) costs under 2 % of the profile's span over half a 20 um cell
    mat = LinearMeltMaterial(T_s=907.0, T_l=908.0, rho_l=2800.0, k=97.0)
    q, T0, dt, L = 1e8, 300.0, 0.5e-3, 5e-3
    sk = layered(mat, [(L, T0)], settings=SkinSettings(cell=20e-6, fine_zone=L, swell=False), cell=20e-6)
    m_min = sk.totals()[0].copy()
    dh = float(mat.enthalpy(mat.T_feed) - mat.enthalpy(T0))
    v_ref = q / (mat.rho * dh)
    alpha = mat._k / (mat.rho * mat._cp)
    steps, settle = 400, 200                     # 0.2 s, measured over its second half (alpha/v^2 is 0.027 s)
    gone = 0.0
    for i in range(steps):
        sk.conduct(dt, q, T0, 0.0, T0)
        mass, _ = sk.remove_feed()
        sk.draw(T0, 0.0, m_min)
        sk.regrid()
        if i >= settle:
            gone += mass[0]
    v = gone / (mat.rho * AREA * (steps - settle) * dt)
    assert v == pytest.approx(v_ref, rel=0.03)
    n = sk.count[0]
    T = sk.T()[0, :n]
    dz = sk.dz()[0, :n]
    z = np.cumsum(dz) - 0.5 * dz
    ref = T0 + (mat.T_feed - T0) * np.exp(-v_ref * z / alpha)
    assert np.max(np.abs(T - ref)) <= 0.02 * (mat.T_feed - T0)


def test_conduction_books_close():
    mat = LinearMeltMaterial(rho_l=2800.0)
    R = 20e-3
    sk = Skins(mat, np.array([AREA, 2 * AREA]), np.array([[1 / R, 0.5 / R], [-1 / R, 0.0]]), SkinSettings())
    sk.fill(mat.rho * sk.area * 1.2e-3, np.array([880.0, 700.0]))
    q, T_amb, eps, T_b, dt = np.array([3e6, 1e6]), 300.0, 0.8, np.array([850.0, 690.0]), 0.010
    E_base = np.zeros(2)
    for _ in range(3):
        H0 = sk.totals()[1]
        book = sk.conduct(dt, q, T_amb, eps, T_b, trial=True)
        H1 = sk.totals()[1]
        scale = np.abs(book.E_top) + np.abs(H1 - H0)
        assert np.all(np.abs(book.E_top - book.E_rad - book.E_base + book.E_side - (H1 - H0)) <= 1e-12 * scale)
        assert np.array_equal(book.E_side, np.zeros(2))
        assert book.E_top == pytest.approx(q * sk.area * dt, rel=1e-15)
        T = sk.T()
        assert book.E_rad == pytest.approx(eps * SIGMA_SB * sk.area * (T[:, 0] ** 4 - T_amb ** 4) * dt, rel=1e-12)
        # the residual is the heat the bottom half-cell conducts into the base at the end of the sub-step
        last = sk.count - 1
        r = np.arange(2)
        dz = sk.dz(T)
        A_b = sk.cross_section(dz.sum(axis=1)[:, None])[:, 0]
        k_last = mat.k(T[r, last])
        assert book.E_base == pytest.approx(A_b * k_last / (0.5 * dz[r, last]) * (T[r, last] - T_b) * dt, rel=1e-6)
        E_base += book.E_base
        sk.regrid()
    a, b = sk.base_law(3 * dt)
    assert np.all(b <= 0.0)
    assert a + b * T_b == pytest.approx(E_base / (3 * dt), rel=1e-12)


def test_the_trial_sensitivity_is_the_derivative_of_the_base_heat():
    # constant k and rho_l = rho: the sub-step's geometry and conductances do not move with T_b, so the linearised
    # step's derivative is exact; mushy cells exercise the latent heat's cp_eff, a draw with a gradient the profile
    mat = LinearMeltMaterial(rho_l=2800.0)
    def start():
        sk = layered(mat, [(0.2e-3, 900.0), (0.3e-3, 885.0), (0.4e-3, 800.0)], settings=SkinSettings(swell=False))
        sk.H[0, -1] *= 0.5                                         # the bottom cell half a cell: the draw has room
        sk.m[0, -1] *= 0.5
        return sk
    m_min = start().totals()[0] * 1.02
    q, T_amb, eps, dt, g = 4e6, 300.0, 0.8, 0.010, 3e4

    def run(T_b, trial):
        sk = start()
        E = 0.0
        dE = 0.0
        for _ in range(4):
            book = sk.conduct(dt, q, T_amb, eps, T_b, trial=trial)
            E += book.E_base[0]
            if trial:
                dE += book.dE_base_dTb[0]
            sk.remove_feed()
            sk.draw(T_b, g, m_min)
            sk.regrid()
        return E, dE, sk

    T_b, d = 780.0, 1e-2
    E, dE, sk = run(T_b, True)
    E_hi, E_lo = run(T_b + d, False)[0], run(T_b - d, False)[0]
    assert dE < 0.0
    assert dE == pytest.approx((E_hi - E_lo) / (2 * d), rel=1e-6)
    a, b = sk.base_law(4 * dt)
    assert b[0] == pytest.approx(dE / (4 * dt), rel=1e-14)
    assert a[0] + b[0] * T_b == pytest.approx(E / (4 * dt), rel=1e-12)


def test_draw_and_hand_down_are_exact():
    mat = LinearMeltMaterial()
    sk = Skins(mat, np.array([AREA, AREA, 2 * AREA]), np.array([[30.0, 10.0], [0.0, 0.0], [-20.0, 0.0]]),
               SkinSettings())
    sk.fill(mat.rho * sk.area * np.array([0.3e-3, 2.5e-3, 0.5e-3]), np.array([950.0, 700.0, 850.0]))
    sk.regrid()
    # the minimum: a stub of solid metal over the base's cross-section
    z_b = sk.thickness()
    A_b = sk.cross_section(z_b[:, None])[:, 0]
    assert sk.minimum_mass() == pytest.approx(mat.rho * A_b * sk.settings.stub, rel=1e-14)
    M0, E0, _ = sk.totals()
    m_min = sk.minimum_mass()
    T_b, g = np.array([890.0, 600.0, 920.0]), np.array([2e4, 0.0, 1e4])
    mass, H = sk.draw(T_b, g, m_min)
    M1, E1, _ = sk.totals()
    assert np.all(np.abs(M1 - M0 - mass) <= 1e-14 * M1)
    assert np.all(np.abs(E1 - E0 - H) <= 1e-14 * np.abs(E1))
    assert mass == pytest.approx(np.maximum(m_min - M0, 0.0), rel=1e-14)
    assert mass[1] == 0.0 and mass[0] > 0.0 and mass[2] > 0.0
    # the profile enthalpy: T(s) = T_b - g s capped at T_rigid, averaged over the depth the shortfall takes as solid
    dd = mass / (mat.rho * A_b)
    s = np.linspace(0.0, 1.0, 20001)[None, :] * dd[:, None]
    h = mat.enthalpy(np.minimum(T_b[:, None] - g[:, None] * s, mat.T_rigid))
    h_avg = 0.5 * (h[:, 1:] + h[:, :-1]).mean(axis=1)
    assert H[0] / mass[0] == pytest.approx(h_avg[0], rel=1e-6)
    assert H[2] / mass[2] == pytest.approx(h_avg[2], rel=1e-6)
    sk.regrid()
    # hand-down: column 1 is all rigid, so everything deeper than the stub below its top leaves, whole cells only
    M1, E1, _ = sk.totals()
    dz = sk.dz()
    z_top = np.cumsum(dz, axis=1) - dz
    leaving = sk.valid & (z_top >= sk.settings.stub * (1 - 1e-12))
    rho = sk.density(sk.T())
    m_before, H_before = sk.m.copy(), sk.H.copy()
    mass, H, vol = sk.hand_down(sk.settings.stub)
    M2, E2, _ = sk.totals()
    assert np.all(np.abs(M1 - M2 - mass) <= 1e-14 * M1)
    assert np.all(np.abs(E1 - E2 - H) <= 1e-14 * np.abs(E1))
    assert mass[1] == pytest.approx(m_before[1][leaving[1]].sum(), rel=1e-14)
    assert H[1] == pytest.approx(H_before[1][leaving[1]].sum(), rel=1e-14)
    assert vol[1] == pytest.approx((m_before[1] / rho[1])[leaving[1]].sum(), rel=1e-14)
    assert mass[1] > 0.0 and mass[0] == 0.0          # column 0: its anchored run is only the drawn metal, under the stub
    assert sk.thickness()[1] >= sk.settings.stub * (1 - 1e-12)
    assert sk.thickness()[1] <= sk.settings.stub + sk.dz()[1, :sk.count[1]].max() * (1 + 1e-12)


def test_a_crust_is_never_handed_down():
    mat = LinearMeltMaterial()
    stub = 0.8e-3
    # rigid crust 0.3 mm over 0.2 mm of melt over 2 mm of rigid metal joined to the floor
    sk = layered(mat, [(0.3e-3, 850.0), (0.2e-3, 950.0), (2.0e-3, 800.0)])
    n = sk.count[0]
    crust, anch = sk.crust()[0, :n], sk.anchored()[0, :n]
    assert crust[:30].all() and not crust[30:].any()
    assert anch[50:].all() and not anch[:50].any()
    top = sk.m[0, :50].copy(), sk.H[0, :50].copy()
    mass, H, vol = sk.hand_down(stub)
    assert sk.count[0] == 50 + 80                                # the anchored run keeps its shallowest 0.8 mm
    assert mass[0] == pytest.approx(mat.rho * AREA * 1.2e-3, rel=1e-12)
    assert vol[0] == pytest.approx(AREA * 1.2e-3, rel=1e-12)
    assert np.array_equal(sk.m[0, :50], top[0]) and np.array_equal(sk.H[0, :50], top[1])
    # with the melt down to the floor nothing is anchored: the crust above it stays whatever the depth
    sk = layered(mat, [(0.3e-3, 850.0), (2.0e-3, 950.0)])
    assert not sk.anchored().any() and sk.crust()[0, :30].all()
    M0 = sk.totals()[0]
    mass, H, vol = sk.hand_down(stub)
    assert mass[0] == 0.0 and H[0] == 0.0 and vol[0] == 0.0 and sk.totals()[0] == M0


def test_a_skin_reduced_to_one_cell_keeps_its_books():
    mat = LinearMeltMaterial()
    sk = Skins(mat, np.array([AREA, AREA]), np.zeros((2, 2)), SkinSettings())
    mc = mat.rho * AREA * 4e-6                                   # one cell of 4 um, thinner than a fine cell
    sk.install(np.array([[mc], [mc]]), np.array([[mc * mat.enthalpy(890.0)], [mc * mat.enthalpy(1000.0)]]),
               np.array([1, 1]))
    m_min = sk.minimum_mass()
    dt, q, T_b = 0.010, np.array([2e6, 2e6]), np.array([880.0, 880.0])
    for _ in range(3):
        M0, E0, _ = sk.totals()
        book = sk.conduct(dt, q, 300.0, 0.8, T_b)
        M1, E1, _ = sk.totals()
        scale = np.abs(book.E_top) + np.abs(E1 - E0)
        assert np.all(np.abs(book.E_top - book.E_rad - book.E_base - (E1 - E0)) <= 1e-12 * scale)
        assert np.array_equal(M1, M0)
        gone, gone_H = sk.remove_feed()
        M2, E2, _ = sk.totals()
        assert np.all(M2 + gone == M1) and np.all(np.abs(E2 + gone_H - E1) <= 1e-14 * np.abs(E1))
        drawn, drawn_H = sk.draw(T_b, 0.0, m_min)
        M3, E3, _ = sk.totals()
        assert np.all(np.abs(M3 - M2 - drawn) <= 1e-14 * M3)
        assert np.all(np.abs(E3 - E2 - drawn_H) <= 1e-14 * np.abs(E3))
        sk.regrid()
        assert np.all(sk.m[sk.valid] > 0.0) and np.all(sk.count >= 1)
    # the liquid column lost its one cell on the first removal and was refilled by the draw at the base's enthalpy
    assert sk.totals()[0] == pytest.approx(m_min, rel=1e-14)
    # an empty column passes its heat straight to the base
    sk.take_top(np.array([0]), np.array([sk.count[0]]))
    book = sk.conduct(dt, q, 300.0, 0.8, T_b)
    assert book.E_base[0] == book.E_top[0] and book.E_rad[0] == 0.0
