"""Tests of the level set's skins (design sections 20.6 and 21, stage 3 and plan 03b): the cells' storage, their grading
below the fine zone, the regrid's exact books, the flat column and the measured cross-section (the curved one was
retired by plan 03b Task 5), and the shared mixture density."""
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
    assert (s.cell, s.fine_zone, s.growth, s.substep, s.swell, s.depth) == (10e-6, 0.4e-3, 1.10, 0.010, True, 1.2e-3)
    assert not hasattr(s, "stub")                            # retired with the full-depth skins (plan 03b Task 5)
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
        short = m_min - sk.totals()[0]                  # the column kept full at its base (plan 03b Task 5's fill)
        if short[0] > 0.0:
            sk.exchange_base([0], short / sk.rho_s, T0)
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
    sk = Skins(mat, np.array([AREA, 2 * AREA]), np.zeros((2, 2)), SkinSettings())
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


TRIAL_FD_RTOL = 1e-5    # the trial derivative against the central difference over +-0.01 K, relative: the heat sub-steps'
                        # Newton tolerance (NEWTON_TOL, 1e-6 K) sets the difference's own error (5.0e-6 measured once,
                        # with constant k, once plan 03b Task 5 took the draw and its exact share of the derivative out;
                        # 1e-6 held with it)


@pytest.mark.parametrize("mat", [LinearMeltMaterial(rho_l=2800.0), VariableKMaterial(rho_l=2800.0)],
                         ids=["constant_k", "variable_k"])
def test_the_trial_sensitivity_is_the_derivative_of_the_base_heat(mat):
    # rho_l = rho: the sub-step's geometry does not move with T_b, so the linearised step's derivative is exact; mushy
    # cells exercise the latent heat's cp_eff; with k rising with T (Task 4) the conductances move with the temperatures
    # and the Jacobian's dk/dT terms carry it. (Plan 03b Task 5 retired the draw within a pass, and with it the drawn
    # metal's sensitivity and the feed removal between sub-steps this test once ran: without the draw's refill the top
    # cells reach T_feed within the four sub-steps, and a removal is not differentiable in T_b.)
    def start():
        return layered(mat, [(0.2e-3, 900.0), (0.3e-3, 885.0), (0.4e-3, 800.0)], settings=SkinSettings(swell=False))
    q, T_amb, eps, dt = 4e6, 300.0, 0.8, 0.010

    def run(T_b, trial):
        sk = start()
        E = 0.0
        dE = 0.0
        for _ in range(4):
            book = sk.conduct(dt, q, T_amb, eps, T_b, trial=trial)
            E += book.E_base[0]
            if trial:
                dE += book.dE_base_dTb[0]
            sk.regrid()
        return E, dE, sk

    T_b, d = 780.0, 1e-2
    E, dE, sk = run(T_b, True)
    E_hi, E_lo = run(T_b + d, False)[0], run(T_b - d, False)[0]
    assert dE < 0.0
    assert dE == pytest.approx((E_hi - E_lo) / (2 * d), rel=TRIAL_FD_RTOL)
    a, b = sk.base_law(4 * dt)
    assert b[0] == pytest.approx(dE / (4 * dt), rel=1e-14)
    assert a[0] + b[0] * T_b == pytest.approx(E / (4 * dt), rel=1e-12)


def test_a_skin_reduced_to_one_cell_keeps_its_books():
    mat = LinearMeltMaterial()
    sk = Skins(mat, np.array([AREA, AREA]), np.zeros((2, 2)), SkinSettings())
    mc = mat.rho * AREA * 4e-6                                   # one cell of 4 um, thinner than a fine cell
    sk.install(np.array([[mc], [mc]]), np.array([[mc * mat.enthalpy(890.0)], [mc * mat.enthalpy(1000.0)]]),
               np.array([1, 1]))
    m_min = sk.totals()[0].copy()                               # each column kept at its mass by the base (the fill)
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
        short = np.maximum(m_min - M2, 0.0)
        rows = np.flatnonzero(short > 0.0)
        ex = sk.exchange_base(rows, short[rows] / sk.rho_s, T_b[rows])
        drawn, drawn_H = np.zeros(2), np.zeros(2)
        drawn[rows], drawn_H[rows] = ex.mass, ex.H
        M3, E3, _ = sk.totals()
        assert np.all(np.abs(M3 - M2 - drawn) <= 1e-14 * M3)
        assert np.all(np.abs(E3 - E2 - drawn_H) <= 1e-14 * np.abs(E3))
        sk.regrid()
        assert np.all(sk.m[sk.valid] > 0.0) and np.all(sk.count >= 1)
    # the liquid column lost its one cell on the first removal and was refilled at its base at the base's enthalpy
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
    sk = layered(mat, [(30e-6, 950.0), (30e-6, 908.5), (30e-6, 900.0), (30e-6, mat.T_rigid - 0.05), (30e-6, 880.0)],
                 settings=SkinSettings(swell=False, rheology="bridge"))      # the test's point is the bridge law's numbers
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


def test_a_narrowing_column_carries_the_shear_force_through_every_interface():
    """On a column whose measured cross-section narrows with depth (beside a convex edge; plan 03b Task 5 retired the
    curved law this test once used, and Task 3's review asked for momentum under a profile) at steady state the force
    tau A0 passes through each interface, A(z_i) D_i (u_i - u_(i+1)), and the stress rises as A0 / A(z_i): 0.5 mm of
    liquid over 0.1 mm of rigid metal, the cross-section falling by 60 % to H."""
    mat = liquid_mat()
    area = AREA * (1.0 - 0.6 * LEVELS / DEPTH)
    sk = profiled(mat, area[None, :])
    z = np.linspace(0.0, 0.6e-3, 61)
    stack(sk, [(z, np.where(np.arange(60) < 50, 950.0, 850.0))])
    sk.momentum_step(STEADY, 1.0 * X, np.zeros(3))
    dz = sk.dz()[0, :50]
    zi = np.cumsum(dz)
    A = sk.cell_section(zi[None, :])[0]
    force = A[:-1] * MU / (0.5 * dz[:-1] + 0.5 * dz[1:]) * (sk.U[0, :49, 0] - sk.U[0, 1:50, 0])
    np.testing.assert_allclose(force, 1.0 * AREA, rtol=1e-9)
    assert A[-1] * MU / (0.5 * dz[-1]) * sk.U[0, 49, 0] == pytest.approx(1.0 * AREA, rel=1e-9)
    assert A[-1] < 0.8 * AREA


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
    # the bridge law, whose slurry viscosity depends on the shear rate: `steady_slurry_rate` iterates its
    # `rheology.viscosity`, and the fixed point takes more than one pass only under a rate-dependent law (the composite
    # law, the default, depends on the solid fraction alone)
    sk = layered(mat, [(0.15e-3, 950.0), (0.15e-3, 900.0), (0.1e-3, 880.0)], settings=SkinSettings(rheology="bridge"))
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


@pytest.mark.parametrize("preset", ["dev", "production"])
def test_a_regridded_column_moved_by_round_off_regrids_to_the_same_cells(preset):
    """A full column under a measured profile (narrowing to half its top's cross-section at H), filled in equal cells of
    the fine size: the fill puts a cell exactly on the merge threshold (the cell whose top is fine_zone + 10 cells deep,
    0.6 mm in the development preset, 0.5 mm in production, is half its target there), so the regrid's decision on it
    is round-off's. Regridded, then its metal scaled by 1 +- 1e-12 (its volume moved by 1e-12 of itself, about what the
    division's closure moves a column by between two regrids, Task 7d), a second regrid changes no cell. (Task 7d round
    2, measured once: before REGRID_RTOL the +1e-12 column regridded from 43 cells to 42 in the development preset and
    from 72 to 71 in production.)"""
    from reentry_model.levelset.skin import MERGE_BELOW
    mat = LinearMeltMaterial()
    area = AREA * (1.0 - 0.5 * LEVELS / DEPTH)
    for nudge in (1.0 + 1e-12, 1.0 - 1e-12):
        sk = profiled(mat, area, settings=SkinSettings(cell=SKIN_PRESETS[preset]))
        sk.fill(sk.density(np.array([890.0])) * sk.column_volume(), 890.0)
        tie = np.abs(sk.dz() / (MERGE_BELOW * sk.targets()) - 1.0)[sk.valid]
        assert tie.min() <= 1e-13                                # the fill leaves a cell on the threshold
        sk.regrid()
        sk.m *= nudge
        sk.H *= nudge
        count, m = sk.count.copy(), sk.m.copy()
        sk.regrid()
        assert np.array_equal(sk.count, count) and np.array_equal(sk.m, m)


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


def test_a_sliver_cell_combined_by_depth_joins_whole_at_its_depth():
    """Plan 03b Task 6 (found when merges took the by-depth rule, the controller's ruling T6-c): a gone column whose
    bottom cell holds 1e-27 kg, its volume some 1e-16 of the column's (left so by earlier cuts and combines on the test
    slab), lies wholly at the round-off of the volumes its overlaps with the keeper's cells are reckoned from, so every
    overlap rounded to nothing and its shares were 0/0 (NaN masses reached the property table). It joins whole the
    keeper's cell at its depth: no floating-point warning, no NaN, and mass, enthalpy and momentum add exactly."""
    mat = ConstMaterial()
    flat = np.full(len(LEVELS), AREA)
    sk = profiled(mat, np.vstack([flat, flat]))
    edges = np.linspace(0.0, DEPTH, 61)
    stack(sk, [(edges, 600.0), (edges, 700.0, np.tile([0.1, 0.0, 0.0], (60, 1)))])
    m, H, U = (np.concatenate([a, np.zeros((2, 1) + a.shape[2:])], axis=1) for a in (sk.m, sk.H, sk.U))
    m[1, 60], H[1, 60], U[1, 60] = 1e-27, 1e-27 * float(mat.enthalpy(700.0)), [0.1, 0.0, 0.0]
    sk.install(m, H, [60, 61], U)
    M0, E0, P0 = books(sk)
    with np.errstate(all="raise"):
        sk.combine([0], [1], by="depth")
    M1, E1, P1 = books(sk)
    assert np.all(np.isfinite(sk.m)) and np.all(np.isfinite(sk.H))
    assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
    assert np.abs(P1 - P0).max() <= 1e-14 * np.abs(P0).max()


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
    then reach below its end at `cell_section`'s cross-section -- the larger of its end's and its mean (plan 03b Task 5:
    the narrowing column's mean, above its base's A(H); a column that pinches at its mean) -- so every cell keeps a
    positive area, while the column itself
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
    A_ext = np.array([V[0] / DEPTH, V[1] / PINCH])
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


def test_a_column_left_with_no_area_combines_on_top_of_its_neighbour():
    """Task 3 fix round 1 (ruling 6: a marker left with no area gives its skin to a neighbour by depth, not to the
    volume). A full 1.2 mm skin, hot (900 K) over its top 0.3 mm and moving, whose profile row is then zeroed -- with no
    top area, or with one -- is given by depth to a full flat keeper at 600 K. Its own cross-section holds no volume, so
    it has no depths: the pair combines by the by-mass rule (bit for bit), mass, enthalpy and momentum add exactly, all
    its heat lands in the keeper's top 0.6 mm (above H), the metal below H is the cold metal the next fill returns, and
    no warning is raised."""
    import warnings
    mat = ConstMaterial()
    flat = np.full(len(LEVELS), AREA)
    rng = np.random.default_rng(23)
    U = rng.normal(0.0, 0.1, (60, 3))
    for top in (0.0, AREA):

        def pair():
            sk = profiled(mat, np.vstack([flat, flat]))
            edges = np.linspace(0.0, DEPTH, 61)
            stack(sk, [(edges, 600.0, U), (edges, np.where(edges[:-1] < 0.3e-3 - 1e-12, 900.0, 600.0), U[::-1])])
            sk.profile[1] = 0.0                                  # its marker left with no area
            sk.area[1] = top
            return sk

        sk, ref = pair(), pair()
        M0, E0, P0 = books(sk)
        scale = np.abs(sk.m[..., None] * sk.U).sum()
        hot = (sk.H - sk.m * sk.props.enthalpy(600.0)).sum()     # the gone column's heat above 600 K
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            sk.combine([0], [1], by="depth")
            ref.combine([0], [1])
            dz = sk.dz()
        for a, b in ((sk.m, ref.m), (sk.H, ref.H), (sk.U, ref.U), (sk.count, ref.count)):
            assert np.array_equal(a, b)                          # the by-mass rule for this pair
        M1, E1, P1 = books(sk)
        assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
        assert np.abs(P1 - P0).max() <= 1e-14 * scale
        assert sk.count[1] == 0 and np.all(np.isfinite(dz))
        n = sk.count[0]
        z_bot = np.cumsum(dz[0, :n])
        excess = sk.H[0, :n] - sk.m[0, :n] * sk.props.enthalpy(600.0)
        assert excess[z_bot <= 0.6e-3 * (1 + 1e-9)].sum() == pytest.approx(hot, rel=1e-12)
        np.testing.assert_allclose(sk.T()[0, :n][z_bot - dz[0, :n] >= DEPTH], 600.0, atol=1e-9)
        assert sk.thickness()[0] == pytest.approx(2 * DEPTH, rel=1e-12)   # overfull by the gone column's metal


def test_metal_in_a_column_with_no_volume_and_no_area_has_no_depth():
    """Task 3 fix round 1: a column holding metal whose profile holds no volume and whose top area is zero has nowhere
    to put it, so its depths raise a clear ValueError instead of turning to NaN; with no metal it has none to place."""
    mat = ConstMaterial()
    sk = profiled(mat, np.vstack([np.full(len(LEVELS), AREA)] * 2))
    sk.fill(mat.rho * sk.column_volume(), 600.0)
    sk.profile[1] = 0.0
    sk.area[1] = 0.0
    with pytest.raises(ValueError, match="neither a profile volume nor a top area"):
        sk.dz()
    with pytest.raises(ValueError, match="neither a profile volume nor a top area"):
        sk.base_area()
    sk.take_top(np.array([1]), np.array([sk.count[1]]))
    assert np.all(sk.dz()[1] == 0.0) and sk.thickness()[0] == pytest.approx(DEPTH, rel=1e-12)


# -- sideways conduction between neighbouring columns (plan 03b Task 4, planning ruling 5) ------------------------------
from scipy.special import erf

from reentry_model.levelset.skin import SideFaces

ROW_SPACING, ROW_WIDTH = 0.25e-3, 1.0e-3      # m: a row's columns' spacing and the width of the boundary they share


def row(mat, T, spacing=ROW_SPACING, width=ROW_WIDTH, settings=None):
    """Columns in a row `spacing` apart, each a prism `spacing` by `width` through the skin depth (a constant measured
    profile), filled full in cells of equal depth at the temperatures T (n,), each sharing with the next a boundary
    `width` long at every level. Returns the skins and their SideFaces, with every column's cross-section and first
    moment (its centroid on the row's axis at every depth: the conductance's centroids, the controller's ruling T6-b of
    plan 03b Task 6, `spacing` apart as the markers are)."""
    T = np.asarray(T, dtype=float)
    n = len(T)
    section = np.full((n, len(LEVELS)), spacing * width)
    sk = profiled(mat, section, settings=settings)
    sk.fill(sk.density(T) * sk.column_volume(), T)
    e = np.arange(n - 1)
    centre = np.zeros((n, len(LEVELS), 3))
    centre[..., 0] = ((np.arange(n) + 0.5) * spacing)[:, None]
    return sk, SideFaces(e, e + 1, np.full((n - 1, len(LEVELS)), width), LEVELS, np.full(n - 1, spacing), section,
                         section[..., None] * centre)


def insulated_step(x, length, x0, s, T_cold, step, images=3):
    """The temperatures [K] at x [m] along a rod `length` [m] long with insulated ends, `step` [K] hotter than T_cold on
    [0, x0] at t = 0, at the time when 2 sqrt(alpha t) = s [m]: the error-function solution summed over the ends'
    images (the hot interval reflected about both ends repeats every 2 length; the controller's ruling P2), exact on the
    finite rod."""
    T = np.full(len(x), float(T_cold))
    for m in range(-images, images + 1):
        a, b = 2.0 * length * m - x0, 2.0 * length * m + x0
        T += 0.5 * step * (erf((x - a) / s) - erf((x - b) / s))
    return T


def test_a_lateral_temperature_step_relaxes_at_the_analytic_rate():
    """Forty columns in a row 0.25 mm apart (10 mm), constant properties, insulated at the top, at the base and at both
    ends, the first twenty 100 K hotter than the rest: after 100 sub-steps to sqrt(alpha t) = 2 mm every cell lies
    within 1 % of the step of the error-function solution summed over the ends' images (the plain error function misses
    the insulated ends by about 3.9 K there, ruling P2). The heat moves only sideways: every sub-step's books close, the
    sideways heat sums to zero and the total enthalpy is kept."""
    mat = ConstMaterial()
    n, T_cold, step = 40, 500.0, 100.0
    sk, side = row(mat, np.where(np.arange(n) < n // 2, T_cold + step, T_cold),
                   settings=SkinSettings(cell=SKIN_PRESETS["dev"]))
    alpha = mat._k / (mat.rho * mat._cp)
    t, N = (2e-3) ** 2 / alpha, 100
    E0 = sk.H.sum()
    for _ in range(N):
        H0 = sk.totals()[1]
        book = sk.conduct(t / N, 0.0, 0.0, 0.0, None, side=side)
        H1 = sk.totals()[1]
        scale = np.abs(book.E_side).sum()
        assert scale > 0.0 and abs(book.E_side.sum()) <= 1e-12 * scale
        assert np.abs(book.E_side - (H1 - H0)).max() <= 1e-9 * scale
    x = (np.arange(n) + 0.5) * ROW_SPACING
    ref = insulated_step(x, n * ROW_SPACING, 0.5 * n * ROW_SPACING, 2.0 * np.sqrt(alpha * t), T_cold, step)
    T = sk.T()
    err = np.abs(np.where(sk.valid, T - ref[:, None], 0.0))
    print("lateral step: largest miss {:.3f} K of the {:.0f} K step ({} cells a column)".format(err.max(), step,
                                                                                            sk.count[0]))
    assert err.max() <= 0.01 * step
    assert sk.H.sum() == pytest.approx(E0, rel=1e-13)


def test_sideways_conduction_closes_its_books():
    """The stage 3 box's 132 columns with the division's measured profiles and shared boundaries (cases.stage3_division:
    boundaries across edges, some with no length at phi), filled full at temperatures rising 160 K along the box, heated
    at the top, radiating and held at the base: over a trial pass of three sub-steps and a real sub-step the heat gained
    sideways sums to zero over the columns to 1e-12 of its magnitude, each column's books close, and the hottest face's
    columns lose heat sideways while the coldest face's gain it."""
    from cases import stage3_box, stage3_division
    box, d = stage3_box(), stage3_division()
    mat = LinearMeltMaterial(rho_l=2800.0)
    sk = Skins(mat, d.area[:, 0], np.zeros((len(d.area), 2)), SkinSettings(cell=SKIN_PRESETS["dev"]))
    sk.set_profile(d.levels, d.area)
    T0 = 600.0 + 2e4 * (box.x[:, 0] - box.x[:, 0].min())
    sk.fill(sk.density(T0) * sk.column_volume(), T0)
    sk.regrid()
    side = SideFaces(d.i, d.j, d.length, d.levels, np.linalg.norm(box.x[d.j] - box.x[d.i], axis=1), d.area, d.moment)
    q, T_amb, eps, dt = 2e6, 300.0, 0.8, 0.010
    T_b = T0 - 5.0
    for trial in (True, True, True, False):
        if not trial and sk.dH is not None:
            sk.base_law(3 * dt)
        H0 = sk.totals()[1]
        book = sk.conduct(dt, q, T_amb, eps, T_b, trial=trial, side=side)
        H1 = sk.totals()[1]
        size = np.abs(book.E_side).sum()
        assert size > 0.0 and abs(book.E_side.sum()) <= 1e-12 * size
        scale = np.abs(book.E_top) + np.abs(H1 - H0) + np.abs(book.E_side)
        assert np.all(np.abs(book.E_top - book.E_rad - book.E_base + book.E_side - (H1 - H0)) <= 1e-12 * scale)
        hottest, coldest = T0 >= T0.max() - 1.0, T0 <= T0.min() + 1.0      # the +x face's columns and the -x face's
        assert book.E_side[hottest].sum() < 0.0 < book.E_side[coldest].sum()
        sk.regrid()


def test_the_sideways_conductance_integrates_the_shared_length_over_the_overlap():
    """Planning ruling 5 on two full columns whose cells do not line up (one in equal cells, the other in cells growing
    with depth and a boundary at 0.6 mm), sharing a boundary whose length falls by 60 % with depth, 2 mm between their
    markers, the conductivity rising with temperature: the first column at 600 K throughout, the second 100 K hotter
    above 0.6 mm and at 600 K below. Over a sub-step of 1 ps (the temperatures move by under 1e-4 K) the heat passed
    sideways is dt 100 K k_h int_0^0.6mm L dz / 2 mm to 1e-6, k_h the harmonic mean of the conductivities at 600 and
    700 K: the shared length integrated over the overlaps of the hot cells alone, whatever the cells."""
    mat = VariableKMaterial(rho_l=2800.0)
    sk = profiled(mat, np.full((2, len(LEVELS)), AREA), settings=SkinSettings(swell=False))
    half = LEVELS[12]
    graded = np.union1d(DEPTH * np.linspace(0.0, 1.0, 23) ** 1.5, [half])
    hot = np.where(graded[:-1] < half, 700.0, 600.0)
    stack(sk, [(np.linspace(0.0, DEPTH, 61), 600.0), (graded, hot)])
    L = 2e-3 * (1.0 - 0.6 * LEVELS / DEPTH)
    side = SideFaces(np.array([0]), np.array([1]), L[None, :], LEVELS, np.array([2e-3]))
    dt = 1e-12
    book = sk.conduct(dt, 0.0, 0.0, 0.0, None, side=side)
    k0, k1 = float(mat.k(600.0)), float(mat.k(700.0))
    expected = dt * 100.0 * 2.0 * k0 * k1 / (k0 + k1) * np.trapezoid(L[:13], LEVELS[:13]) / 2e-3
    assert book.E_side[0] == pytest.approx(expected, rel=1e-6)
    assert book.E_side[1] == -book.E_side[0]
    assert np.abs(sk.T()[1, :sk.count[1]] - hot).max() < 1e-4


def test_the_sideways_conductance_divides_by_the_columns_centroids_over_the_overlap():
    """The controller's ruling T6-b (plan 03b Task 6, amending planning ruling 5): two full columns of one constant
    cross-section in equal cells, their markers 2 mm apart, sharing a boundary 2 mm long at every depth, the first's
    centroid on its marker's axis at every depth and the second's drawing in from 2 mm to 1 mm away with depth (as
    columns across a convex edge do); the second 100 K hotter above 0.6 mm. Over a sub-step of 1 ps the heat passed
    sideways is the sum over the hot cells of dt 100 K k L dz over the distance between the two centroids at the cell's
    mid-depth (both linear in depth, so the centroid over a cell is the one at its middle) to 1e-9, 1.15 times what the
    markers' distance alone passes; with no centroids given, the markers' distance stands."""
    mat = LinearMeltMaterial(rho_l=2800.0)
    sk = profiled(mat, np.full((2, len(LEVELS)), AREA), settings=SkinSettings(swell=False))
    edges = np.linspace(0.0, DEPTH, 61)
    mid = 0.5 * (edges[1:] + edges[:-1])
    hot = np.where(mid < 0.6e-3, 700.0, 600.0)
    stack(sk, [(edges, 600.0), (edges, hot)])
    away = 2e-3 - 1e-3 * LEVELS / DEPTH                            # the second's centroid's distance at each level
    moment = np.zeros((2, len(LEVELS), 3))
    moment[1, :, 0] = AREA * away
    L = np.full((1, len(LEVELS)), 2e-3)
    side = SideFaces(np.array([0]), np.array([1]), L, LEVELS, np.array([2e-3]), np.full((2, len(LEVELS)), AREA),
                     moment)
    dt, k = 1e-12, float(mat.k(650.0))
    expected = dt * 100.0 * k * (2e-3 * np.diff(edges) / np.interp(mid, LEVELS, away))[mid < 0.6e-3].sum()
    marker = dt * 100.0 * k * 2e-3 * 0.6e-3 / 2e-3
    before = sk.state()
    book = sk.conduct(dt, 0.0, 0.0, 0.0, None, side=side)
    print("centroid conductance: {:.6e} J passed, {:.6e} J expected, {:.3f} times the markers' distance's".format(
        -book.E_side[1], expected, -book.E_side[1] / marker))
    assert -book.E_side[1] == pytest.approx(expected, rel=1e-9) and book.E_side[0] == -book.E_side[1]
    assert expected > 1.1 * marker
    sk.set_state(before)
    plain = sk.conduct(dt, 0.0, 0.0, 0.0, None, side=side._replace(section=None, moment=None))
    assert -plain.E_side[1] == pytest.approx(marker, rel=1e-9)


def test_with_sideways_heat_the_trial_sensitivity_is_the_derivative_for_a_common_base_shift():
    """Three columns of different measured profiles in a row (a jagged narrowing one, one that pinches, a constant one)
    at different temperatures, conducting sideways, over a trial pass of four sub-steps with regrids: each column's
    summed dE_base/dT_b is the derivative of its summed base heat when every column's base temperature moves together
    (the trial's derivative with sideways heat), against central differences to 1e-6. rho_l = rho, so the sub-steps'
    geometry does not move with T_b."""
    mat = LinearMeltMaterial(rho_l=2800.0)
    jag = AREA * (1.0 - 0.5 * LEVELS / DEPTH) * np.random.default_rng(7).uniform(0.7, 1.3, len(LEVELS))
    pinch = AREA * np.clip(1.0 - LEVELS / PINCH, 0.0, None)
    area = np.vstack([jag, pinch, np.full(len(LEVELS), AREA)])
    L = np.vstack([np.full(len(LEVELS), 2e-3), 2e-3 * (1.0 - 0.5 * LEVELS / DEPTH)])
    side = SideFaces(np.array([0, 1]), np.array([1, 2]), L, LEVELS, np.array([2e-3, 2e-3]))
    T_start = np.array([890.0, 820.0, 760.0])

    def run(T_b, trial):
        sk = profiled(mat, area, top=np.array([1.1, 0.9, 1.0]) * AREA, settings=SkinSettings(swell=False))
        sk.fill(mat.rho * sk.column_volume(), T_start)
        sk.regrid()
        E, dE = np.zeros(3), np.zeros(3)
        for _ in range(4):
            book = sk.conduct(0.010, 3e6, 300.0, 0.8, T_b, trial=trial, side=side)
            E += book.E_base
            if trial:
                dE += book.dE_base_dTb
            sk.regrid()
        return E, dE

    T_b, d = np.array([850.0, 800.0, 740.0]), 1e-2
    E, dE = run(T_b, True)
    E_hi, E_lo = run(T_b + d, False)[0], run(T_b - d, False)[0]
    fd = (E_hi - E_lo) / (2 * d)
    assert np.all(dE[[0, 2]] < 0.0)
    np.testing.assert_allclose(dE, fd, rtol=1e-6, atol=1e-9 * np.abs(fd).max())


class StalledGMRES:
    """A stand-in for the skin module's `scipy.sparse.linalg` whose `gmres` stalls `share` of the right-hand side's norm
    short of the solution: it returns the exact solution (a sparse LU solve) moved by J^-1 of a residual of that norm in a
    random direction, with info 1 (GMRES's "not converged within its iterations"), and counts its calls; every other
    name is scipy's."""

    def __init__(self, share):
        self.share, self.solves = share, 0

    def gmres(self, A, b, *a, **k):
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        A = sp.csc_matrix(A)
        r = np.random.default_rng(self.solves).standard_normal(len(b))
        self.solves += 1
        return spla.spsolve(A, b) + spla.spsolve(A, self.share * np.linalg.norm(b) * r / np.linalg.norm(r)), 1

    def __getattr__(self, name):
        import scipy.sparse.linalg as spla
        return getattr(spla, name)


SIDE_STALLED_ACCEPTED = 1e-11   # a stalled sideways solve's residual, of its right-hand side, about ten times the round-off
                                # floor measured at 10 um cells (9.9e-13) and under SIDE_ACCEPT_RTOL
SIDE_STALLED_REFUSED = 1e-8     # and one ten times over SIDE_ACCEPT_RTOL


def _stalling_row(mat):
    """Eight columns in a row (`row`, development cells) from 700 to 900 K, across the test material's melting range."""
    return row(mat, np.linspace(700.0, 900.0, 8), settings=SkinSettings(cell=SKIN_PRESETS["dev"]))


def test_a_sideways_solve_stalled_near_round_off_is_accepted_and_the_books_close(monkeypatch):
    """The final review's I3: a sideways heat solve that stalls short of SIDE_RTOL is accepted as an inexact Newton step
    when its residual is at most SIDE_ACCEPT_RTOL of its right-hand side. Eight columns in a row (`_stalling_row`, the
    conductivity rising with temperature), heated at the top, radiating and held 5 K under their temperatures at the
    base, over a trial pass of two sub-steps and a real sub-step, with every sideways solve stalled at
    SIDE_STALLED_ACCEPTED (`StalledGMRES`): every sub-step runs, the skins count every stalled solve
    (`Skins.side_at_floor`), the heat gained sideways sums to zero to 1e-12 of its magnitude and each column's books
    close to 1e-12 (the solve's residual never enters them), and the cells' temperatures are the converged solves'
    within Newton's own tolerance (NEWTON_TOL)."""
    from reentry_model.levelset import skin as skin_mod
    mat = VariableKMaterial()
    q, T_amb, eps, dt = 2e6, 300.0, 0.8, 0.010

    def run():
        sk, side = _stalling_row(mat)
        T_b = sk.T()[np.arange(sk.n), sk.count - 1] - 5.0
        for trial in (True, True, False):
            if not trial:
                sk.base_law(2 * dt)
            H0 = sk.totals()[1]
            book = sk.conduct(dt, q, T_amb, eps, T_b, trial=trial, side=side)
            H1 = sk.totals()[1]
            size = np.abs(book.E_side).sum()
            assert size > 0.0 and abs(book.E_side.sum()) <= 1e-12 * size
            scale = np.abs(book.E_top) + np.abs(H1 - H0) + np.abs(book.E_side)
            assert np.all(np.abs(book.E_top - book.E_rad - book.E_base + book.E_side - (H1 - H0)) <= 1e-12 * scale)
        return sk, book
    converged, book_c = run()
    assert converged.side_at_floor == 0
    stalled = StalledGMRES(SIDE_STALLED_ACCEPTED)
    monkeypatch.setattr(skin_mod, "spla", stalled)
    sk, book = run()
    T, T_c = sk.T(), converged.T()
    print("stalled at {:.0e}: {} solves accepted; temperatures {:.2e} K from the converged solves', base heat {:.2e} of "
          "it".format(SIDE_STALLED_ACCEPTED, sk.side_at_floor, np.abs(np.where(sk.valid, T - T_c, 0.0)).max(),
                      np.abs(book.E_base - book_c.E_base).max() / np.abs(book_c.E_base).max()))
    assert stalled.solves > 0 and sk.side_at_floor == stalled.solves
    assert np.array_equal(sk.count, converged.count)
    assert np.abs(np.where(sk.valid, T - T_c, 0.0)).max() <= skin_mod.NEWTON_TOL


def test_a_sideways_solve_stalled_above_the_floor_raises(monkeypatch):
    """The final review's I3: a sideways heat solve that stalls at SIDE_STALLED_REFUSED of its right-hand side, above
    SIDE_ACCEPT_RTOL, still raises RuntimeError, naming the residual it was left at, and nothing is counted as
    accepted."""
    import re
    from reentry_model.levelset import skin as skin_mod
    sk, side = _stalling_row(VariableKMaterial())
    T_b = sk.T()[np.arange(sk.n), sk.count - 1] - 5.0
    monkeypatch.setattr(skin_mod, "spla", StalledGMRES(SIDE_STALLED_REFUSED))
    with pytest.raises(RuntimeError, match="sideways heat solve did not converge") as stop:
        sk.conduct(0.010, 2e6, 300.0, 0.8, T_b, side=side)
    left = float(re.search(r"residual ([0-9.e+-]+) of the right-hand side", str(stop.value)).group(1))
    assert left == pytest.approx(SIDE_STALLED_REFUSED, rel=0.01)
    assert sk.side_at_floor == 0


# -- the rheology law (plan 03b Task 9) --------------------------------------------------------------------------------

def rheology_column(law):
    """One AA7075 column across the liquidus, the bridge, the slurry and below T_rigid, in the law `law`."""
    mat = Material.from_drama_json("AA7075_scheil")
    layers = [(30e-6, 950.0), (30e-6, 908.5), (30e-6, 907.0), (30e-6, 900.0), (30e-6, 890.0), (30e-6, mat.T_rigid - 0.05),
              (30e-6, 880.0)]
    return mat, layered(mat, layers, settings=SkinSettings(swell=False, rheology=law))


def test_the_default_law_is_the_composite_and_a_bad_name_is_refused():
    assert SkinSettings().rheology == rheology.DEFAULT_LAW == "composite"
    with pytest.raises(ValueError, match="rheology"):
        SkinSettings(rheology="nope")


def test_settings_pickled_before_the_rheology_field_unpickle_as_the_bridge_law():
    import pickle
    old = SkinSettings(rheology="composite")
    state = dict(old.__dict__)
    del state["rheology"]
    fresh = SkinSettings.__new__(SkinSettings)
    fresh.__setstate__(state)
    assert fresh.rheology == "bridge" and fresh.cell == old.cell
    assert pickle.loads(pickle.dumps(SkinSettings(rheology="bridge"))).rheology == "bridge"
    assert pickle.loads(pickle.dumps(old)).rheology == "composite"


def test_composite_viscosity_ignores_the_shear_rate():
    mat, sk = rheology_column("composite")
    T = sk.T()
    flowing = sk.flowing(T)
    f_s = np.minimum(1.0 - sk.props.liquid_fraction(T[flowing]), np.nextafter(rheology.RIGID_SOLID_FRACTION, 0.0))
    want = rheology.law_viscosity("composite", T[flowing], f_s, 0.0, mat.liquid.mu, rheology.bridge_solid_fraction(sk.props))
    base = sk.viscosity(0.0)
    assert np.isfinite(base[flowing]).all() and np.isinf(base[~flowing]).all()
    for gamma in (0.0, 100.0, 5000.0):
        np.testing.assert_allclose(sk.viscosity(gamma), base, rtol=1e-15, atol=0.0)
    np.testing.assert_allclose(base[flowing], want, rtol=1e-15, atol=0.0)
    assert base[flowing].max() > 10.0 * mat.liquid.mu                 # the column does span the slurry


def test_bridge_viscosity_is_unchanged():
    """The bridge law through the shared `law_viscosity` gives today's `rheology.viscosity` values (the skins' own
    property table's solid fraction at the bridge's low end), to 1e-15, and it still depends on the shear rate."""
    mat, sk = rheology_column("bridge")
    T = sk.T()
    flowing = sk.flowing(T)
    f_s = np.minimum(1.0 - sk.props.liquid_fraction(T[flowing]), np.nextafter(rheology.RIGID_SOLID_FRACTION, 0.0))
    f_b = 1.0 - float(sk.props.liquid_fraction(rheology.BRIDGE_LOW))
    for gamma in (0.0, 61.0, 300.0, 5000.0):
        want = rheology.viscosity(T[flowing], f_s, gamma, mat.liquid.mu, f_b)
        got = sk.viscosity(gamma)
        np.testing.assert_allclose(got[flowing], want, rtol=1e-15, atol=0.0)
    assert not np.allclose(sk.viscosity(0.0)[flowing], sk.viscosity(300.0)[flowing], rtol=1e-3)


@pytest.mark.parametrize("law", rheology.LAWS)
def test_the_same_column_gives_the_same_viscosities_in_both_architectures(law):
    from reentry_model.skyscrapers import skin as sky
    mat, ours = rheology_column(law)
    theirs = sky.Skins(mat, ours.area.copy(), sky.SkinSettings(rheology=law))
    theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy())
    flowing = ours.flowing()
    assert flowing.sum() == 12 and (~flowing & ours.valid).sum() == 9     # liquid, bridge and slurry cells, then rigid ones
    T = ours.T()
    f_s = np.minimum(1.0 - ours.props.liquid_fraction(T[flowing]), np.nextafter(rheology.RIGID_SOLID_FRACTION, 0.0))
    for gamma in (0.0, 120.0, 400.0, 3000.0):
        a = ours.viscosity(gamma, T)
        assert np.isinf(a[ours.valid & ~flowing]).all()
        # at the same temperatures and solid fractions the two architectures' laws are the one law: to 1e-15
        np.testing.assert_allclose(a[flowing], sky.viscosity_of(theirs, T[flowing], f_s, gamma), rtol=1e-15, atol=0.0)
        # through each skin's own temperature (each inverts its own enthalpy table, which differ by about 1e-10 in the
        # viscosity on this column, measured once), the same viscosities
        np.testing.assert_allclose(a[flowing], sky.cell_viscosity(theirs, gamma)[flowing], rtol=1e-9, atol=0.0)
