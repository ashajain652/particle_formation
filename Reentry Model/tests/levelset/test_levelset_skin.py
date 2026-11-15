"""Tests of the level set's skins (design section 20.6, stage 3): the cells' storage, their grading below the fine zone,
the regrid's exact books, the curved cross-section A(z) = A0 (1 - k1 z)(1 - k2 z), and the shared mixture density."""
import numpy as np
import pytest

from reentry_model.levelset.skin import SKIN_PRESETS, Skins, SkinSettings
from reentry_model.shared.material import Material

from materials import ConstMaterial, LinearMeltMaterial, VariableKMaterial

AREA = 4e-6                              # m2, a 2 mm square column


def column(mat, thick, T, kappa=(0.0, 0.0), settings=None, area=AREA):
    """One flat-estimate column of `thick` metal at T, filled in cells of about the fine size."""
    sk = Skins(mat, np.array([area]), np.array([kappa], dtype=float), settings or SkinSettings())
    sk.fill(mat.rho * area * thick, T)
    return sk


def test_settings_defaults_and_presets():
    s = SkinSettings()
    assert (s.cell, s.fine_zone, s.growth, s.stub, s.substep, s.swell, s.depth) == \
        (10e-6, 0.4e-3, 1.10, 0.8e-3, 0.010, True, 1.2e-3)
    assert SKIN_PRESETS == {"production": 10e-6, "dev": 20e-6}
    for cell in SKIN_PRESETS.values():                       # every preset passes the depth checks
        assert SkinSettings(cell=cell).depth == 1.2e-3
    with pytest.raises(ValueError):
        SkinSettings(growth=0.9)
    with pytest.raises(ValueError, match="fine zone"):
        SkinSettings(depth=0.3e-3)                          # a fine zone deeper than the skin
    with pytest.raises(ValueError):
        SkinSettings(depth=0.0, fine_zone=0.0)
    with pytest.raises(ValueError, match="two fine cells"):
        SkinSettings(cell=0.25e-3)                          # fewer than two fine cells in the fine zone
    with pytest.raises(ValueError, match="two fine cells"):
        SkinSettings(cell=0.25e-3, fine_zone=0.6e-3, stub=0.4e-3)   # the stub's own check stays until Task 5


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
    sk = layered(mat, [(L, T0)], settings=SkinSettings(cell=20e-6, fine_zone=L, depth=L, swell=False), cell=20e-6)
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


@pytest.mark.parametrize("mat", [LinearMeltMaterial(rho_l=2800.0), VariableKMaterial(rho_l=2800.0)],
                         ids=["constant_k", "variable_k"])
def test_the_trial_sensitivity_is_the_derivative_of_the_base_heat(mat):
    # rho_l = rho: the sub-step's geometry does not move with T_b, so the linearised step's derivative is exact; mushy
    # cells exercise the latent heat's cp_eff, a draw with a gradient the profile; with k rising with T (Task 4) the
    # conductances move with the temperatures and the Jacobian's dk/dT terms carry it
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


def test_a_hand_down_keeps_the_minimum_mass():
    """Stage 3 Task 9: with `m_keep` a hand-down never takes a skin below that mass -- the bottom cells leave only while
    what stays holds at least m_keep -- so that the hand-down (by depth, at each cell's own density) and the draw (by the
    minimum mass, at the solid's density) cannot hand down and draw back the same metal step after step where the
    rigid metal is less dense than the solid (a Scheil alloy's mushy metal below T_rigid)."""
    mat = LinearMeltMaterial()
    stub = 0.8e-3
    sk = layered(mat, [(2.0e-3, 800.0)])
    sk2 = layered(mat, [(2.0e-3, 800.0)])
    free, _, _ = sk.hand_down(stub)
    keep = mat.rho * AREA * 1.5e-3
    mass, H, vol = sk2.hand_down(stub, m_keep=np.array([keep]))
    assert free[0] == pytest.approx(mat.rho * AREA * 1.2e-3, rel=1e-12)
    assert 0.0 < mass[0] <= mat.rho * AREA * 0.5e-3 * (1 + 1e-12)
    assert sk2.totals()[0][0] >= keep * (1 - 1e-12)
    assert sk2.totals()[0][0] < keep + sk2.m[0, :sk2.count[0]].max()     # whole cells: at most one cell over
    nothing, _, _ = layered(mat, [(2.0e-3, 800.0)]).hand_down(stub, m_keep=np.array([1.0]))
    assert nothing[0] == 0.0


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


# -- the layered column flow (Task 4) ----------------------------------------------------------------------------------
from reentry_model.shared import rheology

STEADY = 1.0e10                   # s: one backward-Euler step this long is the steady state (inertia about 1e-11 of it)
MU = 1.3e-3                       # Pa s, the test materials' liquid viscosity (the shared AA7075 liquid's)
X = np.array([1.0, 0.0, 0.0])


def liquid_mat():
    """Constant viscosity: every flowing cell of these tests is above the shared bridge's top (908 K), so it takes the
    liquid's viscosity whatever its shear rate."""
    return LinearMeltMaterial(T_s=880.0, T_l=908.0, rho_l=2400.0)


def heights(sk, n, row=0):
    """The heights [m] of the first n cells' centres above the bottom of the n-th (the rigid front)."""
    dz = sk.dz()[row, :n]
    return dz.sum() - (np.cumsum(dz) - 0.5 * dz)


def momentum_closes(sk, P0, book, rel=1e-12):
    """The momentum change equals the impulses in (push and shear) less the momentum given to the wall, to `rel` of
    the impulses."""
    P1 = sk.totals()[2]
    scale = np.abs(book.push).max() + np.abs(book.shear).max() + np.abs(P0).max()
    return np.abs(P1 - P0 - book.push - book.shear + book.wall).max() <= rel * scale


def test_flowing_and_liquid_follow_the_shared_temperatures():
    mat = Material.from_drama_json("AA7075_scheil")
    sk = layered(mat, [(30e-6, 950.0), (30e-6, 908.5), (30e-6, 900.0), (30e-6, mat.T_rigid - 0.05), (30e-6, 880.0)])
    flowing, liquid = sk.flowing()[0], sk.liquid()[0]
    assert np.array_equal(flowing, np.arange(15) < 9)                  # above T_rigid (half liquid): flows
    assert np.array_equal(liquid, np.arange(15) < 6)                   # at or above the liquidus: liquid
    eta = sk.viscosity(0.0)[0]
    assert np.all(eta[:6] == mat.liquid.mu) and np.all(np.isinf(eta[9:]))
    f_s = 1.0 - mat.liquid_fraction(np.array([900.0]))
    f_b = 1.0 - mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))
    assert eta[6] == pytest.approx(float(rheology.viscosity(900.0, f_s, 0.0, mat.liquid.mu, f_b)[0]), rel=1e-12)


def test_linear_couette():
    """Gas shear tau on 0.3 mm of liquid over rigid metal, steady: u = tau s / mu at every cell centre (s the height
    above the rigid front), the bottom cell dragging on the rigid metal over its own half-thickness; the surface speed,
    extrapolated from the top cell, tau b / mu within 1e-6."""
    mat = liquid_mat()
    sk = layered(mat, [(0.3e-3, 950.0), (0.2e-3, 850.0)], settings=SkinSettings())
    P0 = sk.totals()[2]
    book = sk.momentum_step(STEADY, 1.0 * X, np.zeros(3))
    s = heights(sk, 30)
    b = sk.dz()[0, :30].sum()
    np.testing.assert_allclose(sk.U[0, :30, 0], s / MU, rtol=1e-9)
    top = sk.U[0, 0, 0] + 0.5 * sk.dz()[0, 0] * (sk.U[0, 0, 0] - sk.U[0, 1, 0]) / (s[0] - s[1])
    assert top == pytest.approx(1.0 * b / MU, rel=1e-6)
    assert np.all(sk.U[0, 30:] == 0.0) and np.all(sk.U[0, :, 1:] == 0.0) and book.passes == 1
    assert book.shear[0] == pytest.approx(1.0 * AREA * STEADY * X, rel=1e-15)
    assert momentum_closes(sk, P0, book)


def test_a_curved_column_carries_the_shear_force_through_every_interface():
    """On a convex column the cross-section narrows with depth, so at steady state the force tau A0 passes through each
    interface, A(z_i) D_i (u_i - u_(i+1)), and the stress rises as A0 / A(z_i)."""
    mat = liquid_mat()
    sk = layered(mat, [(0.5e-3, 950.0), (0.1e-3, 850.0)], kappa=(1 / 2e-3, 1 / 4e-3), settings=SkinSettings())
    sk.momentum_step(STEADY, 1.0 * X, np.zeros(3))
    dz = sk.dz()[0, :50]
    z = np.cumsum(dz)[:-1]
    A = sk.cross_section(z[None, :])[0]
    force = A * MU / (0.5 * dz[:-1] + 0.5 * dz[1:]) * (sk.U[0, :49, 0] - sk.U[0, 1:50, 0])
    np.testing.assert_allclose(force, 1.0 * AREA, rtol=1e-9)
    A_b = sk.cross_section(np.array([[dz.sum()]]))[0, 0]
    assert A_b * MU / (0.5 * dz[-1]) * sk.U[0, 49, 0] == pytest.approx(1.0 * AREA, rel=1e-9)


def test_half_parabola_second_order():
    """A body force G alone on a free liquid layer over rigid metal, steady: u = (G / mu)(h s - s^2 / 2); the error
    falls by four when the cells halve."""
    G = 1e4
    mat = liquid_mat()

    def error(N):
        sk = layered(mat, [(0.3e-3, 950.0), (0.1e-3, 850.0)], cell=0.3e-3 / N, settings=SkinSettings())
        P0 = sk.totals()[2]
        book = sk.momentum_step(STEADY, np.zeros(3), G * X)
        assert momentum_closes(sk, P0, book)
        s, h = heights(sk, N), sk.dz()[0, :N].sum()
        exact = G / MU * (h * s - 0.5 * s * s)
        return np.abs(sk.U[0, :N, 0] - exact).max() / exact.max()
    e15, e30, e60 = error(15), error(30), error(60)
    assert e30 < 1e-3
    assert e15 / e30 == pytest.approx(4.0, rel=0.1) and e30 / e60 == pytest.approx(4.0, rel=0.1)


def steady_slurry_rate(mat, T, tau):
    """The steady shear rate [1/s] of slurry at T [K] under the stress tau [Pa]: tau over its own viscosity at that
    rate, by fixed point on the shared rheology."""
    f_s = 1.0 - float(mat.liquid_fraction(np.array([T]))[0])
    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    gamma = 0.0
    for _ in range(60):
        gamma = tau / float(rheology.viscosity(T, f_s, gamma, mat.liquid.mu, f_b))
    return gamma


def test_liquid_over_slurry_kinks():
    """The shared AA7075 rheology: 1 Pa of gas shear on 15 liquid cells (950 K) over 15 of slurry (900 K, 37 % solid)
    over rigid metal. The stress is 1 Pa all through, so each layer is linear and the profile's slope jumps at the
    interface by the viscosity ratio, the slurry's at its own shear rate (found by the fixed point) over the liquid's."""
    mat = Material.from_drama_json("AA7075_scheil")
    sk = layered(mat, [(0.15e-3, 950.0), (0.15e-3, 900.0), (0.1e-3, 880.0)], settings=SkinSettings())
    P0 = sk.totals()[2]
    book = sk.momentum_step(STEADY, 1.0 * X, np.zeros(3))
    s, u = heights(sk, 30), sk.U[0, :30, 0]
    slope_l = (u[0] - u[14]) / (s[0] - s[14])
    slope_s = (u[15] - u[29]) / (s[15] - s[29])
    gamma = steady_slurry_rate(mat, 900.0, 1.0)
    eta_s = 1.0 / gamma
    assert slope_l / slope_s == pytest.approx(eta_s / mat.liquid.mu, rel=1e-3)
    assert slope_l == pytest.approx(1.0 / mat.liquid.mu, rel=1e-9)
    assert 1 < book.passes <= 5
    assert momentum_closes(sk, P0, book)


def test_stokes_first_problem():
    """Tau suddenly applied to the top of a deep liquid layer at rest: u(z, t) = (2 tau / mu)[sqrt(nu t / pi)
    exp(-z^2 / 4 nu t) - (z / 2) erfc(z / 2 sqrt(nu t))], the surface (2 tau / mu) sqrt(nu t / pi) (read from the top
    cell's centre and the surface gradient tau / mu the shear sets); the surface and the cells within 2 % of the surface
    speed at three times (diffusion lengths sqrt(nu t) of 40, 80 and 120 um in 1.2 mm), and the
    momentum the shear's impulse exactly."""
    from scipy.special import erfc
    mat = liquid_mat()
    nu = MU / mat.liquid.rho
    sk = layered(mat, [(1.2e-3, 950.0), (0.05e-3, 850.0)], settings=SkinSettings())
    N = 120
    dz = sk.dz()[0, :N]
    z = np.cumsum(dz) - 0.5 * dz
    dt, t = 2e-5, 0.0
    for L in (40e-6, 80e-6, 120e-6):
        t_end = L * L / nu
        steps = int(round((t_end - t) / dt))
        for _ in range(steps):
            sk.momentum_step(dt, 1.0 * X, np.zeros(3))
            t += dt
        L = np.sqrt(nu * t)
        exact = 2.0 / MU * (L / np.sqrt(np.pi) * np.exp(-z * z / (4 * L * L)) - 0.5 * z * erfc(z / (2 * L)))
        surface = sk.U[0, 0, 0] + 1.0 * 0.5 * dz[0] / MU          # the top cell's centre plus the shear's gradient
        assert surface == pytest.approx(2.0 / MU * L / np.sqrt(np.pi), rel=0.02)
        assert np.abs(sk.U[0, :N, 0] - exact).max() <= 0.02 * exact[0]
        assert sk.totals()[2][0, 0] == pytest.approx(1.0 * AREA * t, rel=1e-12)


def test_a_floating_crust_moves_as_one_piece():
    """A crust (three rigid cells) on liquid over the rigid metal joined to the floor, every cell moving at the start,
    sheared and pushed: the crust's cells end at one speed, the anchored cells at rest (their momentum to the wall), and
    the momentum change is the impulses less the wall's to 1e-12. Steady under shear alone, the crust rides the
    liquid's Couette profile at its top. A wholly rigid column stops and gives the wall all its momentum."""
    mat = liquid_mat()
    sk = layered(mat, [(30e-6, 850.0), (0.15e-3, 950.0), (50e-6, 850.0)], settings=SkinSettings())
    two = Skins(mat, np.full(2, AREA), np.zeros((2, 2)), SkinSettings())
    rigid = layered(mat, [(0.23e-3, 850.0)])
    two.install(np.vstack([sk.m, rigid.m]), np.vstack([sk.H, rigid.H]), np.array([23, 23]))
    assert two.crust()[0, :3].all() and not two.crust()[0, 3:].any() and two.anchored()[0, 18:].all()
    two.U = np.random.default_rng(3).normal(0.0, 0.05, two.U.shape) * two.valid[..., None]
    P0 = two.totals()[2]
    book = two.momentum_step(1e-3, np.array([[3.0, -1.0, 0.0], [2.0, 0.0, 0.0]]), np.array([2e4, -1e4, 0.0]))
    assert np.ptp(two.U[0, :3], axis=0).max() == 0.0 and np.all(two.U[0, 0] != 0.0)       # one piece
    assert np.all(two.U[0, 18:] == 0.0) and np.all(two.U[1] == 0.0)
    assert book.wall[1] == pytest.approx(P0[1], rel=1e-15) and np.all(book.shear[1] == 0.0)
    assert momentum_closes(two, P0, book)
    # steady under shear alone: the stress is tau all through, the crust at the liquid's speed at its top
    book = sk.momentum_step(STEADY, 1.0 * X, np.zeros(3))
    s = heights(sk, 18)[3:]
    h_l = sk.dz()[0, 3:18].sum()
    np.testing.assert_allclose(sk.U[0, 3:18, 0], s / MU, rtol=1e-6)
    np.testing.assert_allclose(sk.U[0, :3, 0], h_l / MU, rtol=1e-6)


def _three_columns(mat):
    """Three columns of 4, 3 and 5 cells, each cell at its own temperature and velocity (a profile), areas 4, 3, 5 mm2."""
    rng = np.random.default_rng(9)
    sk = Skins(mat, np.array([4e-6, 3e-6, 5e-6]), np.array([[50.0, 0.0], [0.0, 0.0], [-20.0, 10.0]]),
               SkinSettings(cell=20e-6))
    count = np.array([4, 3, 5])
    m = np.where(np.arange(5)[None, :] < count[:, None], rng.uniform(1e-7, 4e-7, (3, 5)), 0.0)
    T = rng.uniform(700.0, 950.0, (3, 5))
    sk.install(m, m * mat.enthalpy(T), count, rng.normal(0.0, 1.0, (3, 5, 3)))
    return sk


def test_cut_gives_fractions_of_a_column_with_its_profile_and_keeps_the_rest():
    """Task 9b (fix round 1, minor 2): `Skins.cut` gives each new column its fraction of every cell of its source --
    mass and enthalpy times the fraction, velocity kept, the same count -- and that fraction of the area, with the
    source's curvatures; the source keeps the rest of every cell and of its area. A row cut twice keeps one less both
    fractions. Mass, enthalpy and momentum add to round-off."""
    sk = _three_columns(LinearMeltMaterial())
    m0, H0, U0, A0, k0 = sk.m.copy(), sk.H.copy(), sk.U.copy(), sk.area.copy(), sk.kappa.copy()
    P0 = (sk.m[..., None] * sk.U).sum(axis=(0, 1))
    new = sk.cut([0, 2, 0], [0.25, 0.4, 0.5])
    assert list(new) == [3, 4, 5] and sk.n == 6
    for row, src, f in ((3, 0, 0.25), (4, 2, 0.4), (5, 0, 0.5)):
        assert sk.count[row] == sk.count[src]
        np.testing.assert_allclose(sk.m[row], f * m0[src], rtol=1e-15, atol=0.0)
        np.testing.assert_allclose(sk.H[row], f * H0[src], rtol=1e-15, atol=0.0)
        assert np.array_equal(sk.U[row], U0[src]) and np.array_equal(sk.kappa[row], k0[src])
        assert sk.area[row] == pytest.approx(f * A0[src], rel=1e-15)
    np.testing.assert_allclose(sk.m[0], 0.25 * m0[0], rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(sk.m[2], 0.6 * m0[2], rtol=1e-14, atol=0.0)
    assert np.array_equal(sk.m[1], m0[1]) and sk.area[0] == pytest.approx(0.25 * A0[0], rel=1e-14)
    P1 = (sk.m[..., None] * sk.U).sum(axis=(0, 1))
    assert sk.m.sum() == pytest.approx(m0.sum(), rel=1e-15) and sk.H.sum() == pytest.approx(H0.sum(), rel=1e-15)
    np.testing.assert_allclose(P1, P0, rtol=0.0, atol=1e-15 * np.abs(m0[..., None] * U0).sum())


def test_cut_of_a_whole_column_empties_it_within_its_slack():
    """Task 9b (fix round 1, minor 2): fractions of one column summing to one -- exactly, or short of it by less than
    CUT_SLACK (here by 1e-16) -- give the whole column away: the source is left empty with no area, and the fractions
    are scaled to sum to one, so the pieces hold all its metal to round-off. Short by more than the slack (1e-12), the
    source keeps the rest."""
    from reentry_model.levelset.skin import CUT_SLACK
    for short in (0.0, 1e-16):
        sk = _three_columns(LinearMeltMaterial())
        m0, H0, A0 = sk.m[1].sum(), sk.H[1].sum(), sk.area[1]
        assert short < CUT_SLACK
        new = sk.cut([1, 1], [0.3, 0.7 - short])
        assert sk.count[1] == 0 and sk.m[1].sum() == 0.0 and sk.area[1] == 0.0
        assert sk.m[new].sum() == pytest.approx(m0, rel=1e-15) and sk.H[new].sum() == pytest.approx(H0, rel=1e-15)
        assert sk.area[new].sum() == pytest.approx(A0, rel=1e-15)
    sk = _three_columns(LinearMeltMaterial())
    m0 = sk.m[1].sum()
    sk.cut([1, 1], [0.3, 0.7 - 1e-12])
    assert sk.count[1] == 3 and sk.m[1].sum() == pytest.approx(1e-12 * m0, rel=1e-3)


def test_cut_refuses_bad_fractions_and_an_open_trial():
    """Task 9b (fix round 1, minor 2): `Skins.cut` refuses fractions that are not positive and finite, a fraction count
    other than the rows', fractions of one column summing past one by more than CUT_SLACK, and a trial pass open (its
    sensitivities would not follow the cut); a refused cut leaves the columns as they were."""
    sk = _three_columns(LinearMeltMaterial())
    m0 = sk.m.copy()
    for rows, fractions in (([0], [0.0]), ([0], [-0.1]), ([0], [np.nan]), ([0, 1], [0.5]), ([0, 0], [0.6, 0.5])):
        with pytest.raises(ValueError):
            sk.cut(rows, fractions)
    assert np.array_equal(sk.m, m0) and sk.n == 3
    sk.conduct(1e-3, 0.0, 0.0, 0.0, 800.0, trial=True)
    with pytest.raises(RuntimeError, match="trial pass"):
        sk.cut([0], [0.5])


def test_append_bottom_adds_whole_cells_and_refuses_a_gap():
    """Stage 3 Task 10 (`Skins.append_bottom`, the volume's metal a top-zone skin takes when it turns full-depth): new
    cells at rest below each column's bottom, top first, the books adding exactly; a row's trailing zeros end it, but a
    zero-mass cell with metal after it is refused (fix round 1: it was silently dropped, with the metal after it)."""
    mat = ConstMaterial()
    sk = Skins(mat, np.full(2, 1e-6), np.zeros((2, 2)), SkinSettings(cell=20e-6))
    sk.fill(mat.rho * 1e-6 * 0.4e-3, 600.0)
    M0, E0 = sk.m.sum(), sk.H.sum()
    c0 = sk.count.copy()
    m = np.array([[1e-9, 2e-9, 0.0], [3e-9, 0.0, 0.0]])
    sk.append_bottom([0, 1], m, 1e3 * m)
    assert np.array_equal(sk.count, c0 + [2, 1])
    assert sk.m.sum() == pytest.approx(M0 + m.sum(), rel=1e-15) and sk.H.sum() == pytest.approx(E0 + 1e3 * m.sum(), rel=1e-15)
    assert np.all(sk.U == 0.0)
    with pytest.raises(ValueError, match="zero"):
        sk.append_bottom([0], [[1e-9, 0.0, 2e-9]], [[1e-6, 0.0, 2e-6]])


# -- one skin depth: measured cross-sections, the base exchange and combining by depth (plan 03b Task 3) ---------------
DEPTH = 1.2e-3                               # m, the skin depth H (SkinSettings.depth)
LEVELS = np.linspace(0.0, DEPTH, 25)         # the division's levels (division.DIVISION_LEVELS = 24 bands)
PINCH = LEVELS[16]                           # m, 0.8 mm: where a pinching column's cross-section closes


def profile_volume(area, z, levels=LEVELS):
    """The test's oracle: the volume [m3] under the cross-section `area` (B+1,), linear between `levels`, from the top
    down to the depth z [m] (at most H): the trapezoid rule on the levels above z and z itself, exact for a
    piecewise-linear integrand."""
    z = min(float(z), float(levels[-1]))
    pts = np.union1d(levels[levels < z], [z])
    return float(np.trapezoid(np.interp(pts, levels, area), pts))


def profiled(mat, area, top=None, settings=None):
    """Empty columns with the measured cross-sections `area` (n, B+1) at LEVELS; `top` the areas the top fluxes act on,
    by default the profile's own top."""
    area = np.atleast_2d(np.asarray(area, dtype=float))
    top = area[:, 0].copy() if top is None else np.asarray(top, dtype=float)
    sk = Skins(mat, top, np.zeros((len(area), 2)), settings or SkinSettings())
    sk.set_profile(LEVELS, area)
    return sk


def stack(sk, cells):
    """Every column's cells installed by depth: cells[i] = (the depths [m] of column i's cell boundaries, top first,
    the cells' temperatures [K][, their velocities (k, 3) [m/s]]); each cell holds the metal its span of the column's
    profile holds at its own density."""
    K = max(len(c[0]) - 1 for c in cells)
    m, H, U = np.zeros((sk.n, K)), np.zeros((sk.n, K)), np.zeros((sk.n, K, 3))
    count = np.zeros(sk.n, dtype=np.int64)
    for i, c in enumerate(cells):
        V = np.array([profile_volume(sk.profile[i], z) for z in c[0]])
        T = np.broadcast_to(np.asarray(c[1], dtype=float), (len(V) - 1,))
        k = len(V) - 1
        m[i, :k] = sk.density(T) * np.diff(V)
        H[i, :k] = m[i, :k] * sk.props.enthalpy(T)
        if len(c) > 2:
            U[i, :k] = c[2]
        count[i] = k
    sk.install(m, H, count, U)


def books(sk):
    """Mass [kg], enthalpy [J] and momentum [kg m/s] summed over every column."""
    M, E, P = sk.totals()
    return M.sum(), E.sum(), P.sum(axis=0)


def test_a_measured_profile_holds_a_spherical_shell():
    """A column under a 25 mm sphere's cross-section A0 ((R - z)/R)^2, given at the division's levels: its volume is the
    shell's to the trapezoid rule's error (h^2 / 6 R^2 = 7e-7 at 50 um levels); filled with that volume, its cells end
    at H and each spans the depths over which the profile holds it. The top fluxes' area stays the one given (ruling
    T1-g) and the curvature is ignored."""
    R, mat = 25e-3, ConstMaterial()
    area = AREA * ((R - LEVELS) / R) ** 2
    sk = profiled(mat, area, top=1.05 * AREA)
    assert sk.area[0] == 1.05 * AREA                         # set_profile leaves the top fluxes' area alone
    V = sk.column_volume()[0]
    assert V == pytest.approx(AREA * (R ** 3 - (R - DEPTH) ** 3) / (3 * R ** 2), rel=1e-6)
    assert V == pytest.approx(profile_volume(area, DEPTH), rel=1e-14)
    z = np.array([[0.0, 0.31e-3, 0.6e-3, DEPTH, 1.5e-3]])
    A = sk.cross_section(z)[0]
    np.testing.assert_allclose(A[:4], np.interp(z[0, :4], LEVELS, area), rtol=1e-14)
    assert A[4] == 0.0                                       # nothing beyond H
    sk.fill(mat.rho * V, 500.0)
    sk.regrid()
    n = sk.count[0]
    assert sk.thickness()[0] == pytest.approx(DEPTH, rel=1e-12)
    assert sk.cells_volume()[0] == pytest.approx(V, rel=1e-14)
    z_bot = np.cumsum(sk.dz()[0, :n])
    np.testing.assert_allclose([profile_volume(area, zb) for zb in z_bot], np.cumsum(sk.m[0, :n]) / mat.rho, rtol=1e-12)
    dz, t = sk.dz()[0, :n], sk.targets()[0, :n]
    assert np.all((dz >= 0.5 * t) & (dz <= 1.5 * t))
    assert sk.base_area()[0] == pytest.approx(area[-1], rel=1e-14)       # a full column's base is the profile's A(H)
    for nudge in (1.0 + 1e-15, 1.0 - 1e-15):                 # whichever side of its volume round-off leaves the cells
        sk.m[0, :n] *= nudge
        sk.H[0, :n] *= nudge
        assert sk.base_area()[0] == pytest.approx(area[-1], rel=1e-14)
    sk.kappa[:] = 1.0 / 2e-3                                 # a curvature that would close the column at 2 mm: ignored
    assert sk.thickness()[0] == pytest.approx(DEPTH, rel=1e-12)


def test_a_column_that_pinches_fills_exactly_to_its_pinch():
    """A column whose cross-section falls linearly to zero at 0.8 mm (beside a convex edge, it meets the column across
    the edge there): zero below the pinch and beyond H; 0.3 mm of metal at its top, filled at its base by the exchange
    to its volume, ends at the pinch to 1e-12 m, holding its volume to 1e-12, and the regrid settles there, where the
    cross-section closes."""
    mat = LinearMeltMaterial()
    area = AREA * np.clip(1.0 - LEVELS / PINCH, 0.0, None)
    sk = profiled(mat, area)
    V = sk.column_volume()[0]
    assert V == pytest.approx(0.5 * AREA * PINCH, rel=1e-14)
    A = sk.cross_section(np.array([[0.4e-3, 0.79e-3, PINCH, 0.9e-3, DEPTH, 1.5e-3]]))[0]
    assert A[0] == pytest.approx(0.5 * AREA, rel=1e-12) and A[1] > 0.0 and np.all(A[2:] == 0.0)
    sk.fill(mat.rho * profile_volume(area, 0.3e-3), 700.0)
    assert sk.thickness()[0] == pytest.approx(0.3e-3, rel=1e-12)
    sk.exchange_base([0], V - sk.cells_volume()[0], 650.0)
    sk.regrid()
    assert abs(sk.thickness()[0] - PINCH) <= 1e-12
    assert sk.cells_volume()[0] == pytest.approx(V, rel=1e-12)
    n = sk.count[0]
    dz, t = sk.dz()[0, :n], sk.targets()[0, :n]
    assert np.all((dz >= 0.5 * t) & (dz <= 1.5 * t))
    assert sk.base_area()[0] == 0.0                          # it has no base


def test_conduction_books_close_on_a_measured_profile():
    """Two columns with measured profiles -- a jagged narrowing one (level-to-level jumps of up to 30 %, as the
    division's edge columns have, ruling T1-h) and one that pinches -- and top areas other than the profiles' tops:
    the heat books close over a trial pass, the top fluxes act on the top areas, the base heat goes through the
    profile's base (none where the column pinches, whose bottom is insulated). At steady state the heat q A_top passes
    through every interface over the profile's cross-section there."""
    mat = LinearMeltMaterial(rho_l=2800.0)
    jag = AREA * (1.0 - 0.5 * LEVELS / DEPTH) * np.random.default_rng(5).uniform(0.7, 1.3, len(LEVELS))
    pinch = AREA * np.clip(1.0 - LEVELS / PINCH, 0.0, None)
    sk = profiled(mat, np.vstack([jag, pinch]), top=np.array([1.1, 0.9]) * AREA)
    sk.fill(mat.rho * sk.column_volume(), np.array([880.0, 700.0]))
    sk.regrid()
    q, T_amb, eps, T_b, dt = np.array([3e6, 1e6]), 300.0, 0.8, np.array([850.0, 690.0]), 0.010
    E_base = np.zeros(2)
    r = np.arange(2)
    for _ in range(3):
        H0 = sk.totals()[1]
        book = sk.conduct(dt, q, T_amb, eps, T_b, trial=True)
        H1 = sk.totals()[1]
        scale = np.abs(book.E_top) + np.abs(H1 - H0)
        assert np.all(np.abs(book.E_top - book.E_rad - book.E_base + book.E_side - (H1 - H0)) <= 1e-12 * scale)
        assert book.E_top == pytest.approx(q * sk.area * dt, rel=1e-15)
        T = sk.T()
        assert book.E_rad == pytest.approx(eps * SIGMA_SB * sk.area * (T[:, 0] ** 4 - T_amb ** 4) * dt, rel=1e-12)
        last, dz = sk.count - 1, sk.dz(T)
        A_b = sk.base_area()
        assert A_b[0] == pytest.approx(jag[-1], rel=1e-9) and A_b[1] == 0.0
        k_last = mat.k(T[r, last])
        assert book.E_base[0] == pytest.approx(A_b[0] * k_last[0] / (0.5 * dz[0, last[0]]) * (T[0, last[0]] - T_b[0])
                                               * dt, rel=1e-6)
        assert abs(book.E_base[1]) <= 1e-6 * book.E_top[1]    # only the iteration's residue: the bottom is insulated
        E_base += book.E_base
        sk.regrid()
    a, b = sk.base_law(3 * dt)
    assert np.all(b <= 0.0) and b[1] == 0.0
    assert a + b * T_b == pytest.approx(E_base / (3 * dt), rel=1e-12)
    # steady: one step of 1e6 s; constant k, so each interface carries q A_top through A(z_i) k dT / (dz_i/2 + dz_j/2)
    const = ConstMaterial()
    st = profiled(const, jag, top=1.1 * AREA, settings=SkinSettings(swell=False))
    st.fill(const.rho * st.column_volume(), 600.0)
    st.regrid()
    st.conduct(1e6, 1e6, 0.0, 0.0, 600.0)
    n = st.count[0]
    T, dz = st.T()[0, :n], st.dz()[0, :n]
    A_i = st.cross_section(np.cumsum(dz)[None, :-1])[0]
    flux = A_i * const._k * (T[:-1] - T[1:]) / (0.5 * dz[:-1] + 0.5 * dz[1:])
    np.testing.assert_allclose(flux, 1e6 * 1.1 * AREA, rtol=1e-6)


def test_the_base_exchange_is_exact():
    """Draw and return at the base. A return takes the bottom volume whatever its phase -- here a molten bottom cell
    (above the liquidus, swollen) whole and half the cell above it -- and gives its mass, enthalpy and momentum; a draw
    adds a new cell at rest below the bottom of the volume asked for, at the solid's density and the drawn temperature.
    What the exchange reports is what crossed the base into each column (signed), so the columns' books change by it
    exactly (1e-12)."""
    from reentry_model.levelset.skin import BaseExchange
    mat = LinearMeltMaterial()
    area = AREA * (1.0 - 0.3 * LEVELS / DEPTH)
    sk = profiled(mat, np.vstack([area, area, area]))
    rng = np.random.default_rng(11)
    edges = np.linspace(0.0, 0.9e-3, 37)
    T0 = np.r_[np.full(34, 700.0), 890.0, 950.0]                 # a mushy cell, then a molten one at the bottom
    stack(sk, [(edges, T0, rng.normal(0.0, 0.1, (36, 3))), (edges, 700.0, rng.normal(0.0, 0.1, (36, 3))),
               (edges, 700.0)])
    assert sk.T()[0, 35] > mat.T_liquidus
    M0, E0, P0 = sk.totals()
    m0, H0, U0 = sk.m.copy(), sk.H.copy(), sk.U.copy()
    vol = sk.m / sk.density(sk.T())
    dV = np.array([-(vol[0, 35] + 0.5 * vol[0, 34]), 0.1e-3 * area[-1], 0.0])
    ex = sk.exchange_base([0, 1, 2], dV, T_draw=np.array([np.nan, 800.0, np.nan]))
    assert isinstance(ex, BaseExchange)
    assert np.array_equal(ex.cells, [-1, 1, 0])
    M1, E1, P1 = sk.totals()
    assert np.all(np.abs(M1 - M0 - ex.mass) <= 1e-12 * M0)
    assert np.all(np.abs(E1 - E0 - ex.H) <= 1e-12 * np.abs(E0))
    assert np.abs(P1 - P0 - ex.P).max() <= 1e-12 * np.abs(m0[..., None] * U0).sum()
    # the return: the molten bottom cell whole, half the one above at its own enthalpy and velocity
    assert -ex.mass[0] == pytest.approx(m0[0, 35] + 0.5 * m0[0, 34], rel=1e-12)
    assert -ex.H[0] == pytest.approx(H0[0, 35] + 0.5 * H0[0, 34], rel=1e-12)
    np.testing.assert_allclose(-ex.P[0], m0[0, 35] * U0[0, 35] + 0.5 * m0[0, 34] * U0[0, 34], rtol=1e-12)
    assert sk.count[0] == 35 and sk.m[0, 34] == pytest.approx(0.5 * m0[0, 34], rel=1e-12)
    assert sk.T()[0, 34] == pytest.approx(T0[34], abs=1e-9) and np.array_equal(sk.U[0, 34], U0[0, 34])
    assert sk.cells_volume()[0] == pytest.approx(vol[0].sum() + dV[0], rel=1e-12)
    # the draw: a new bottom cell at rest, rho_s dV at 800 K
    assert ex.mass[1] == pytest.approx(mat.rho * dV[1], rel=1e-15) and np.all(ex.P[1] == 0.0)
    assert ex.H[1] == pytest.approx(ex.mass[1] * float(sk.props.enthalpy(800.0)), rel=1e-15)
    assert sk.count[1] == 37 and sk.T()[1, 36] == pytest.approx(800.0, abs=1e-9) and np.all(sk.U[1, 36] == 0.0)
    assert sk.cells_volume()[1] == pytest.approx(vol[1].sum() + dV[1], rel=1e-12)
    assert np.array_equal(sk.m[2, :36], m0[2]) and sk.count[2] == 36 and ex.mass[2] == 0.0
    # a return of a whole column empties it; more than a column holds, or a draw with no temperature, is refused
    ex = sk.exchange_base([2], -sk.cells_volume()[2])
    assert sk.count[2] == 0 and ex.cells[0] == -36 and -ex.mass[0] == pytest.approx(m0[2].sum(), rel=1e-14)
    with pytest.raises(ValueError):
        sk.exchange_base([0], -1.01 * sk.cells_volume()[0])
    with pytest.raises(ValueError, match="T_draw"):
        sk.exchange_base([0], 1e-12)


def test_combining_by_depth_keeps_each_depth_s_metal():
    """A pinching column (to 0.8 mm) given to a flat one, both hot (900 K) in their top 0.3 mm over cold (600 K) metal,
    each cell moving: by depth, each of the keeper's cells takes the gone column's metal at its own depths, so the
    keeper's 15 cells above 0.3 mm hold only hot metal and the cells below only cold, its cell boundaries do not move
    (its profile is the sum of the two), and mass, enthalpy and momentum add exactly. By mass (today's rule, the
    default), the gone column's hot metal is spread down to 61 % of the keeper's mass."""
    mat = ConstMaterial()
    flat = np.full(len(LEVELS), AREA)
    pinch = AREA * np.clip(1.0 - LEVELS / PINCH, 0.0, None)
    rng = np.random.default_rng(13)

    def pair():
        sk = profiled(mat, np.vstack([flat, pinch]))
        ek, eg = np.linspace(0.0, DEPTH, 61), np.linspace(0.0, PINCH, 33)
        stack(sk, [(ek, np.where(ek[:-1] < 0.3e-3 - 1e-12, 900.0, 600.0), rng.normal(0.0, 0.1, (60, 3))),
                   (eg, np.where(eg[:-1] < 0.3e-3 - 1e-12, 900.0, 600.0), rng.normal(0.0, 0.1, (32, 3)))])
        return sk

    sk = pair()
    M0, E0, P0 = books(sk)
    dz0 = sk.dz()[0, :60].copy()
    scale = np.abs(sk.m[..., None] * sk.U).sum()
    sk.combine([0], [1], by="depth")
    M1, E1, P1 = books(sk)
    assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
    assert np.abs(P1 - P0).max() <= 1e-14 * scale
    assert sk.count[0] == 60 and sk.count[1] == 0
    np.testing.assert_allclose(sk.profile[0], flat + pinch, rtol=1e-15)
    assert np.all(sk.profile[1] == 0.0) and sk.area[1] == 0.0
    np.testing.assert_allclose(sk.dz()[0, :60], dz0, rtol=1e-12)        # the keeper's boundaries do not move
    T = sk.T()[0, :60]
    np.testing.assert_allclose(T[:15], 900.0, atol=1e-9)                 # the hot top stays above the cold bottom
    np.testing.assert_allclose(T[15:], 600.0, atol=1e-9)
    by_mass = pair()
    by_mass.combine([0], [1])
    T = by_mass.T()[0, :by_mass.count[0]]
    deep = np.cumsum(by_mass.dz()[0, :by_mass.count[0]]) - by_mass.dz()[0, :by_mass.count[0]] > 0.31e-3
    assert T[deep].max() > 650.0                                         # by mass the hot metal goes deeper


def test_taking_metal_between_two_depths_is_exact():
    """`take_between` removes the metal between two depths of a narrowing column, splitting the two cells the depths
    cut: what it takes is the profile's volume between them at the cells' density, mass, enthalpy and momentum are
    exact, the cells above keep their metal and the cells below move up. Two depths inside one cell take a share of
    it."""
    mat = ConstMaterial()
    area = AREA * (1.0 - 0.4 * LEVELS / DEPTH)
    sk = profiled(mat, np.vstack([area, area]))
    rng = np.random.default_rng(17)
    edges = np.linspace(0.0, DEPTH, 61)
    U = rng.normal(0.0, 0.1, (60, 3))
    stack(sk, [(edges, rng.uniform(500.0, 900.0, 60), U), (edges, 700.0, U)])
    M0, E0, P0 = sk.totals()
    m0, H0, U0 = sk.m.copy(), sk.H.copy(), sk.U.copy()
    z_top, z_bot = 0.333e-3, 0.871e-3                       # inside cells 16 and 43 (20 um cells)
    m, H, P = sk.take_between([0], z_top, z_bot)
    assert m[0] == pytest.approx(mat.rho * (profile_volume(area, z_bot) - profile_volume(area, z_top)), rel=1e-12)
    M1, E1, P1 = sk.totals()
    assert M1[0] + m[0] == pytest.approx(M0[0], rel=1e-14) and E1[0] + H[0] == pytest.approx(E0[0], rel=1e-14)
    assert np.abs(P1[0] + P[0] - P0[0]).max() <= 1e-14 * np.abs(m0[0, :, None] * U0[0]).sum()
    assert np.array_equal(sk.m[0, :16], m0[0, :16]) and np.array_equal(sk.H[0, :16], H0[0, :16])
    assert sk.m[0, 16] == pytest.approx(mat.rho * (profile_volume(area, z_top) - profile_volume(area, edges[16])),
                                        rel=1e-10)
    assert sk.count[0] == 60 - 26                           # cells 17 to 42 wholly between the depths
    assert sk.m[0, 17] == pytest.approx(mat.rho * (profile_volume(area, edges[44]) - profile_volume(area, z_bot)),
                                        rel=1e-10)
    assert np.array_equal(sk.m[0, 18:34], m0[0, 44:60]) and np.array_equal(sk.U[0, 18:34], U0[0, 44:60])
    assert np.array_equal(sk.m[1], m0[1])
    # two depths inside one cell take its share between them; the cell stays, at its own enthalpy and velocity
    m, H, P = sk.take_between([1], 0.402e-3, 0.407e-3)
    assert sk.count[1] == 60 and m[0] == pytest.approx(mat.rho * (profile_volume(area, 0.407e-3)
                                                                  - profile_volume(area, 0.402e-3)), rel=1e-10)
    assert sk.m[1, 20] + m[0] == pytest.approx(m0[1, 20], rel=1e-14) and np.array_equal(sk.U[1, 20], U0[1, 20])


def test_profiled_columns_divide_cut_and_combine_with_their_profiles():
    """Ruling T3-a: a column divided, cut or combined takes its share of the profile as it does of the top area and
    the cells (halved, the fraction, the sum), so every column holds exactly its volume afterwards and the cells keep
    their depths; a removed column's profile row goes with it."""
    mat = ConstMaterial()
    jag = AREA * (1.0 - 0.5 * LEVELS / DEPTH) * np.random.default_rng(19).uniform(0.7, 1.3, len(LEVELS))
    pinch = AREA * np.clip(1.0 - LEVELS / PINCH, 0.0, None)
    flat = np.full(len(LEVELS), AREA)
    sk = profiled(mat, np.vstack([jag, pinch, flat]), top=np.array([1.1, 0.9, 1.0]) * AREA)
    sk.fill(mat.rho * sk.column_volume(), np.array([600.0, 700.0, 800.0]))
    sk.regrid()
    p0, V0, dz0 = sk.profile.copy(), sk.column_volume(), sk.dz()
    new = sk.divide([0])
    assert list(new) == [3]
    np.testing.assert_array_equal(sk.profile[0], 0.5 * p0[0])
    np.testing.assert_array_equal(sk.profile[3], 0.5 * p0[0])
    np.testing.assert_allclose(sk.column_volume()[[0, 3]], 0.5 * V0[0], rtol=1e-15)
    np.testing.assert_allclose(sk.dz()[3], dz0[0], rtol=1e-12, atol=0.0)
    pieces = sk.cut([1, 1], [0.3, 0.2])
    np.testing.assert_allclose(sk.profile[pieces], [0.3 * p0[1], 0.2 * p0[1]], rtol=1e-15)
    np.testing.assert_allclose(sk.profile[1], 0.5 * p0[1], rtol=1e-14)
    for row in (1, *pieces):
        np.testing.assert_allclose(sk.dz()[row, :sk.count[row]], dz0[1, :sk.count[1]], rtol=1e-10, atol=0.0)
    sk.combine([2], [pieces[0]])
    np.testing.assert_allclose(sk.profile[2], flat + 0.3 * p0[1], rtol=1e-15)
    np.testing.assert_allclose(sk.cells_volume(), sk.column_volume(), rtol=1e-12)    # every column exactly full
    sk.remove([pieces[0]])
    assert sk.profile.shape == (sk.n, len(LEVELS)) and sk.n == 5
    np.testing.assert_allclose(sk.cells_volume(), sk.column_volume(), rtol=1e-12)


def test_an_overfull_column_continues_below_its_end():
    """Between a gain at its top and the next fill a column may hold more than its volume (here 5 % more). Its cells
    then reach below its end at `cell_section`'s cross-section -- a column with a base at its base's, A(H), continued,
    one that pinches at its mean cross-section -- so every cell keeps a positive area, while the column itself
    (`cross_section`) still ends there and a column that pinches still has no base. The heat books close, and a
    return of the excess at the base leaves both exactly full, ending at H and at the pinch."""
    mat = ConstMaterial()
    narrow = AREA * (1.0 - 0.5 * LEVELS / DEPTH)
    pinch = AREA * np.clip(1.0 - LEVELS / PINCH, 0.0, None)
    sk = profiled(mat, np.vstack([narrow, pinch]))
    V = sk.column_volume()
    extra = 0.05 * V
    sk.fill(mat.rho * (V + extra), 600.0)
    sk.regrid()
    A_ext = np.array([narrow[-1], V[1] / PINCH])
    np.testing.assert_allclose(sk.thickness(), np.array([DEPTH, PINCH]) + extra / A_ext, rtol=1e-12)
    z = np.array([[1.3e-3], [1.0e-3]])
    np.testing.assert_allclose(sk.cell_section(z)[:, 0], A_ext, rtol=1e-14)
    assert np.all(sk.cross_section(z) == 0.0)
    np.testing.assert_allclose(sk.base_area(), [narrow[-1], 0.0], rtol=1e-14)
    n = sk.count
    dz, t = sk.dz(), sk.targets()
    assert all(np.all((dz[i, :n[i]] >= 0.5 * t[i, :n[i]]) & (dz[i, :n[i]] <= 1.5 * t[i, :n[i]])) for i in range(2))
    H0 = sk.totals()[1]
    book = sk.conduct(0.010, 2e6, 300.0, 0.8, 600.0)
    scale = np.abs(book.E_top) + np.abs(sk.totals()[1] - H0)
    assert np.all(np.abs(book.E_top - book.E_rad - book.E_base - (sk.totals()[1] - H0)) <= 1e-12 * scale)
    assert abs(book.E_base[1]) <= 1e-6 * book.E_top[1]          # the pinching column's bottom stays insulated
    sk.exchange_base([0, 1], -extra)
    np.testing.assert_allclose(sk.cells_volume(), V, rtol=1e-12)
    assert abs(sk.thickness()[0] - DEPTH) <= 1e-12 and abs(sk.thickness()[1] - PINCH) <= 1e-12
