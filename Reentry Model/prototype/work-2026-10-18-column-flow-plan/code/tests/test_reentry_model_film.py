"""film.py: the lubrication branches, the implicit upwind runoff (conservation, positivity, the strip steady state)."""
import numpy as np
import pytest

from reentry_model import film, material, mesh

LIQ = material.LiquidProperties(2400.0, 1.3e-3, 0.86)


def test_lubrication_branches():
    tau, G, mu = 30.0, 5e4, LIQ.mu
    V, q, rate, thick = film.lubrication([tau, tau, tau], [G, G, G], [5e-5, 5e-4, 5e-4], [1e-4, 1e-4, np.nan], mu)
    b, d = 5e-5, 1e-4
    assert not thick[0] and V[0] == pytest.approx(tau * b / mu + G * b * b / (2 * mu)) and q[0] == pytest.approx(tau * b ** 2 / (2 * mu) + G * b ** 3 / (3 * mu))
    assert rate[0] == pytest.approx(V[0] / b)
    b = 5e-4
    assert thick[1] and V[1] == pytest.approx(tau * d / mu) and q[1] == pytest.approx(V[1] * d / 2 + G * b ** 3 / (3 * mu)) and rate[1] == pytest.approx(V[1] / d)
    assert not thick[2] and V[2] == pytest.approx(tau * b / mu + G * b * b / (2 * mu))               # no delta_m: thin branch
    assert film.lubrication([tau], [-1e6], [1e-3], [np.nan], mu)[1][0] < 0.0                         # a strong adverse gradient reverses q


def strip_surface(n=20, w=1e-3):
    """A planar strip of n squares (2 triangles each) along +x, width w, in the plane z = 0 with the normal +z."""
    x = np.arange(n + 1) * w
    pts = np.array([[xi, yi, 0.0] for xi in x for yi in (0.0, w)])
    faces = []
    for k in range(n):
        a, b, c, d = 2 * k, 2 * k + 1, 2 * k + 2, 2 * k + 3
        faces += [[a, c, b], [b, c, d]]
    faces = np.array(faces)
    p = pts[faces]
    normals = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    areas = 0.5 * np.linalg.norm(normals, axis=1)
    normals = normals / (2 * areas)[:, None]
    return pts, mesh.SurfaceMesh(faces, p.mean(axis=1), normals, areas, np.arange(len(faces)), np.arange(len(faces)))


def test_runoff_conserves_mass_and_reaches_the_strip_steady_state():
    """Constant shear on a strip fed at its first patch: in the steady state the thin-film flux through every edge
    equals the feed, so b = sqrt(2 mu S / (rho tau w)) on the interior patches (closed form) within 0.1 %."""
    pts, s = strip_surface()
    ro = film.Runoff(s, pts)
    assert len(ro.i) == 2 * (s.n_patches // 2) - 1                                     # 20 diagonals + 19 shared verticals
    t_hat = np.tile([1.0, 0.0, 0.0], (s.n_patches, 1))
    tau, w = 30.0, 1e-3
    S = 2e-6                                                                            # kg/s fed into patch 0
    q_of_b = lambda b: film.lubrication(np.full(s.n_patches, tau), np.zeros(s.n_patches), b, np.full(s.n_patches, np.nan), LIQ.mu)[1]
    m = np.zeros(s.n_patches)
    dt = 0.05
    for k in range(400):
        m[0] += S * dt
        total = m.sum()
        m, n, moved = ro.transport(m, q_of_b, t_hat, LIQ.rho, s.areas, dt, substeps=1)   # one implicit step per feed pulse
        assert m.sum() == pytest.approx(total, rel=1e-12) and (m >= 0.0).all() and n == 1
        m[-2:] = 0.0                                                                    # the strip's end is stripped (a sink)
    b = m / (LIQ.rho * s.areas)
    b_exact = np.sqrt(2.0 * LIQ.mu * S / (LIQ.rho * tau * w))
    assert np.abs(b[4:-4] / b_exact - 1.0).max() < 1e-3 and moved > 0.0


def test_runoff_on_the_sphere_stops_at_the_equator(coarse_sphere_mesh):
    s = coarse_sphere_mesh.surface()
    v = np.array([1.0, 0.0, 0.0])
    theta, t_hat = s.angles_to(v), s.tangent_from(v)
    windward = theta <= np.pi / 2
    ro = film.Runoff(s, coarse_sphere_mesh.points, windward)
    assert windward[ro.i].all() and windward[ro.j].all()
    m0 = np.where(theta < 0.5, 1e-4 * s.areas * LIQ.rho, 0.0)
    tau = np.where(windward, 30.0 * np.sin(theta), 0.0)
    G = np.where(windward, 5e4 * np.sin(theta) * np.cos(theta), 0.0)
    q_of_b = lambda b: film.lubrication(tau, G, b, np.full(b.size, np.nan), LIQ.mu)[1]
    m, n, moved = ro.transport(m0, q_of_b, t_hat, LIQ.rho, s.areas, 0.5)
    assert m.sum() == pytest.approx(m0.sum(), rel=1e-12) and (m >= 0.0).all() and moved > 0.0
    assert m[~windward].sum() == 0.0 and np.degrees(theta[m > 1e-9 * m.max()].max()) > 35.0        # spread outward (>= 4 patches), never leeward
    assert ro.transport(np.zeros(s.n_patches), q_of_b, t_hat, LIQ.rho, s.areas, 0.5)[1] == 0    # a dry surface costs nothing


# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-10-02: the liquid below the conjugate depth runs off under the pressure gradient and deceleration


def test_deep_flux_is_the_poiseuille_flux_and_completes_fact_29s_two_depth_limits():
    """With no film on top, the liquid beneath it carries the body-force-driven lubrication flux of a layer of depth h,
    q = G h^3 / (3 mu): the half-channel profile u(y) = (G / mu)(h y - y^2 / 2), no slip at the solid and no stress
    from this part at the free surface, integrated over the depth. Under a film of thickness b it carries the rest of
    that flux, so that the thick branch's film flux plus this is plan fact 29's column flux
    tau d^2 / (2 mu) + G h^3 / (3 mu) with d = min(delta_m, h): the shear part over the conjugate depth only, the
    pressure part over the whole liquid depth. Since the amendment of 2026-10-05 that holds where the film fills the
    conjugate layer; a thinner film carries only its own part of the layer's shear flux, so the column is short by the
    shear flux of the rest of the layer, tau (delta_m - b)^2 / (2 mu), which lies in the elements and no account moves."""
    mu, G, h = LIQ.mu, 4.0e4, 1.0e-3
    y = np.linspace(0.0, h, 200001)
    u = G / mu * (h * y - 0.5 * y * y)
    assert film.deep_flux(G, 0.0, h, mu) == pytest.approx(np.trapezoid(u, y), rel=1e-9)
    assert film.deep_flux(G, 0.0, h, mu) == pytest.approx(G * h ** 3 / (3.0 * mu), rel=1e-12)
    tau, delta_m, b = 30.0, 2.5e-4, np.array([1e-5, 1e-4, 4e-4])       # films thinner and thicker than delta_m
    V, q, _, thick = film.lubrication(np.full(3, tau), np.full(3, G), b, np.full(3, delta_m), mu, b_layer=np.full(3, h))
    assert thick.all()                                                  # the column is deeper than delta_m
    column = tau * delta_m ** 2 / (2.0 * mu) + G * h ** 3 / (3.0 * mu)
    held = tau * (delta_m - np.minimum(b, delta_m)) ** 2 / (2.0 * mu)  # the shear flux of the layer below a thin film
    assert q + film.deep_flux(G, b, h, mu) == pytest.approx(np.full(3, column) - held, rel=1e-12)
    assert held[2] == 0.0 and (held[:2] > 0.0).all()                    # exact where the film fills the conjugate layer
    assert film.deep_flux(G, 2e-3, h, mu) == 0.0                        # no liquid beneath a film as deep as the column
    assert film.deep_flux(-G, 0.0, h, mu) == pytest.approx(-G * h ** 3 / (3.0 * mu))   # the deceleration's pull: toward the nose
    # a thin column (h <= delta_m) is all within the shear's reach: lubrication alone is the column, unchanged
    _, q_thin, _, thick_thin = film.lubrication([tau], [G], [5e-5], [delta_m], mu, b_layer=[2e-4])
    assert not thick_thin[0] and q_thin[0] == pytest.approx(tau * 5e-5 ** 2 / (2.0 * mu) + G * 5e-5 ** 3 / (3.0 * mu))


# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-10-05: on a thick patch the film carries only its own part of the conjugate layer's flux


def test_a_thick_patch_moves_only_the_film_present():
    """The thick branch's shear-driven velocity is linear within Girin's conjugate layer -- V at the surface, zero at
    depth delta_m -- and the film is the top b of that layer, so a film thinner than delta_m carries
    V b (1 - b / (2 delta_m)): the integral of the profile over its own depth. That goes to zero with the film, keeps
    q / b below V however thin the film, and equals the whole layer's V delta_m / 2 at b = delta_m. The surface
    velocity, the shear rate and the branch are those of before, and so is everything at b >= delta_m and on thin
    patches, bit for bit. Before the amendment a vanishing film on a thick patch was handed the whole layer's
    V delta_m / 2, so q / b grew like 1/b (the crash of the 0.0125 s time-step study, plan facts 69-70)."""
    mu, tau, G, dm = LIQ.mu, 30.0, 4.0e4, 2.5e-4
    V = tau * dm / mu
    b = np.array([5e-324, 1e-300, 1e-200, 1e-100, 1e-30, 1e-12, 1e-9, 1e-6, 1e-5, 1e-4, 2.4e-4])
    n = b.size
    with np.errstate(over="ignore"):        # the thin branch's V / b, evaluated everywhere and discarded on thick patches
        V_s, q, rate, thick = film.lubrication(np.full(n, tau), np.zeros(n), b, np.full(n, dm), mu, b_layer=np.full(n, 2e-3))
    assert thick.all() and np.array_equal(V_s, np.full(n, V)) and np.array_equal(rate, np.full(n, V / dm))
    assert (q > 0.0).all() and (q <= V * b).all()                       # vanishes with the film; q / b <= V
    assert q[1:] / b[1:] == pytest.approx(V * (1.0 - b[1:] / (2.0 * dm)), rel=1e-12)
    for bb in (1e-5, 1e-4, 2.4e-4):                                     # the profile's integral over the top b
        z = np.linspace(0.0, bb, 100001)
        q1 = film.lubrication([tau], [0.0], [bb], [dm], mu, b_layer=[2e-3])[1][0]
        assert q1 == pytest.approx(np.trapezoid(V * (1.0 - z / dm), z), rel=1e-9)
    # continuous at b = delta_m, with the pressure and deceleration part G b^3 / (3 mu) on top as before
    q_at = [film.lubrication([tau], [G], [x], [dm], mu, b_layer=[2e-3])[1][0] for x in (dm * (1 - 1e-9), dm, dm * (1 + 1e-9))]
    assert q_at[0] == pytest.approx(q_at[1], rel=1e-8) and q_at[2] == pytest.approx(q_at[1], rel=1e-8)
    # unchanged, bit for bit: a film at least delta_m deep on a thick patch, and every thin patch
    bt = np.array([dm, 2.6e-4, 1e-3, 5e-3])
    tt, gt, dt_ = np.full(4, tau), np.full(4, G), np.full(4, dm)
    V_old = tt * dt_ / mu
    assert np.array_equal(film.lubrication(tt, gt, bt, dt_, mu, b_layer=bt + 1e-3)[1], V_old * dt_ / 2.0 + gt * bt ** 3 / (3.0 * mu))
    for bn, dm_thin, extra in ((np.array([1e-9, 5e-5, 2e-4, dm]), dt_, 0.0),               # the layer within delta_m
                               (np.array([1e-9, 5e-5, 1e-3, 5e-3]), np.full(4, np.nan), 1e-3)):  # no conjugate depth
        V_n, q_n, _, thick_n = film.lubrication(tt, gt, bn, dm_thin, mu, b_layer=bn + extra)
        assert not thick_n.any()
        assert np.array_equal(q_n, tt * bn * bn / (2.0 * mu) + gt * bn ** 3 / (3.0 * mu))
        assert np.array_equal(V_n, tt * bn / mu + gt * bn * bn / (2.0 * mu))


def test_runoff_stays_finite_on_nearly_dry_thick_patches():
    """Two thick patches whose surface directions point at each other drain into each other -- as they do in the
    stagnation region and in the craters of the eroded front. With a film of 1e-40 kg on each, the thick branch's flux
    before the amendment gave them emptying rates of 1e37 per second, the linearly implicit system lost its pivots to
    cancellation (1 + c dt - c dt = 0) and spsolve returned NaN for every patch: the failure of the 0.0125 s flight,
    "Matrix is exactly singular" included. Now the rate is bounded by V l / A, and the transport stays finite,
    non-negative and conservative for nearly dry films all the way down to the smallest positive double. The deep
    liquid's flux needed no change: its rate q_deep / (A h_D) <= |G| (b + h_D)^2 / (mu A) stays bounded as h_D -> 0."""
    import warnings
    from scipy.sparse.linalg import MatrixRankWarning
    pts, s = strip_surface()
    ro = film.Runoff(s, pts)
    n = s.n_patches
    t_hat = np.tile([1.0, 0.0, 0.0], (n, 1))
    t_hat[1::2] = [-1.0, 0.0, 0.0]                                      # the two triangles of each square face each other
    tau, dm, G = 30.0, 2.5e-4, np.full(n, -4.0e4)
    molten = np.full(n, 2e-3)                                           # 2 mm of molten material under every patch
    q_of_b = lambda b: film.lubrication(np.full(n, tau), G, b, np.full(n, dm), LIQ.mu, b_layer=b + molten)[1]
    V = tau * dm / LIQ.mu
    bound = (V + 4.0e4 * (1e-5 * 1.001) ** 2 / (3.0 * LIQ.mu)) * ro.length.max() / s.areas.min()
    for tiny in (1e-20, 1e-40, 1e-100, 1e-200, 1e-300, np.finfo(float).tiny, 5e-324):
        m = np.zeros(n)
        m[10:14] = 1e-5 * LIQ.rho * s.areas[10:14]                      # ten microns of film on two squares
        m[6] = m[7] = tiny                                              # and a nearly dry pair [kg]
        b = m / (LIQ.rho * s.areas)
        with warnings.catch_warnings(), np.errstate(over="ignore"):    # (the thin branch's V / b on a 5e-324 film)
            assert np.concatenate(ro.edge_coefficients(q_of_b(b), b, t_hat, s.areas)).max() <= bound
            warnings.simplefilter("error", MatrixRankWarning)           # no singular matrix
            out = m.copy()
            for _ in range(8):                                          # eight macro steps of 0.0125 s, 32 sub-steps
                out, k, moved = ro.transport(out, q_of_b, t_hat, LIQ.rho, s.areas, 0.0125)
        assert np.isfinite(out).all() and (out >= 0.0).all() and k == 4
        assert out.sum() == pytest.approx(m.sum(), rel=1e-12, abs=0.0)
    # the deep liquid beneath a film, nearly gone: a bounded rate, a finite transport
    b_film = np.full(n, 1e-4)
    deep = lambda hd: film.deep_flux(G, b_film, b_film + hd, LIQ.mu)
    for tiny in (1e-40, 1e-200, 5e-324):
        D = np.zeros(n)
        D[6] = D[7] = tiny
        D[10:14] = 1e-4 * LIQ.rho * s.areas[10:14]
        hd = D / (LIQ.rho * s.areas)
        c = np.concatenate(ro.edge_coefficients(deep(hd), hd, t_hat, s.areas))
        assert c.max() <= 4.0e4 * (2e-4 * 1.001) ** 2 / LIQ.mu * ro.length.max() / s.areas.min()
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            out, _, _ = ro.transport(D, deep, t_hat, LIQ.rho, s.areas, 0.0125)
        assert np.isfinite(out).all() and (out >= 0.0).all() and out.sum() == pytest.approx(D.sum(), rel=1e-12, abs=0.0)
