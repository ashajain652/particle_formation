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
    pressure part over the whole liquid depth."""
    mu, G, h = LIQ.mu, 4.0e4, 1.0e-3
    y = np.linspace(0.0, h, 200001)
    u = G / mu * (h * y - 0.5 * y * y)
    assert film.deep_flux(G, 0.0, h, mu) == pytest.approx(np.trapezoid(u, y), rel=1e-9)
    assert film.deep_flux(G, 0.0, h, mu) == pytest.approx(G * h ** 3 / (3.0 * mu), rel=1e-12)
    tau, delta_m, b = 30.0, 2.5e-4, np.array([1e-5, 1e-4, 4e-4])       # films thinner and thicker than delta_m
    V, q, _, thick = film.lubrication(np.full(3, tau), np.full(3, G), b, np.full(3, delta_m), mu, b_layer=np.full(3, h))
    assert thick.all()                                                  # the column is deeper than delta_m
    column = tau * delta_m ** 2 / (2.0 * mu) + G * h ** 3 / (3.0 * mu)
    assert q + film.deep_flux(G, b, h, mu) == pytest.approx(np.full(3, column), rel=1e-12)
    assert film.deep_flux(G, 2e-3, h, mu) == 0.0                        # no liquid beneath a film as deep as the column
    assert film.deep_flux(-G, 0.0, h, mu) == pytest.approx(-G * h ** 3 / (3.0 * mu))   # the deceleration's pull: toward the nose
    # a thin column (h <= delta_m) is all within the shear's reach: lubrication alone is the column, unchanged
    _, q_thin, _, thick_thin = film.lubrication([tau], [G], [5e-5], [delta_m], mu, b_layer=[2e-4])
    assert not thick_thin[0] and q_thin[0] == pytest.approx(tau * 5e-5 ** 2 / (2.0 * mu) + G * 5e-5 ** 3 / (3.0 * mu))


def test_nearly_dry_patches_that_feed_each_other_do_not_make_the_runoff_singular():
    """Post-reconstruction addition (2026-10-07): two adjacent patches holding 1e-39 kg whose tangents point at each
    other, with a flux that does not vanish with the film (the thick branch's V_s delta_m / 2): their emptying numbers
    are ~1e33, the 1 of each diagonal is lost to round-off and, before EMPTYING_MAX, SuperLU found the matrix exactly
    singular and returned NaN (the 100 mm default flight's failure at ~110 s). Capped, the solve is finite, positive
    and conservative, and a third patch's film is untouched."""
    import warnings
    pts, s = strip_surface(4)
    ro = film.Runoff(s, pts)
    t_hat = np.zeros((s.n_patches, 3))
    t_hat[0], t_hat[1] = [1.0, 1.0, 0.0], [-1.0, -1.0, 0.0]                            # across their shared diagonal
    t_hat /= np.maximum(np.linalg.norm(t_hat, axis=1), 1e-300)[:, None]
    m0 = np.zeros(s.n_patches)
    m0[:3] = [1e-39, 3e-39, 1e-6]
    with warnings.catch_warnings():
        warnings.simplefilter("error")                                                  # a MatrixRankWarning fails the test
        m, n, _ = ro.transport(m0, lambda b: np.full(b.size, 1e-4), t_hat, LIQ.rho, s.areas, 0.5)
    assert np.isfinite(m).all() and (m >= 0.0).all() and m.sum() == pytest.approx(m0.sum(), rel=1e-12)
    assert m[0] + m[1] == pytest.approx(4e-39, rel=1e-9) and m[2] == m0[2] and ro.capped_rows > 0
    # below the bound nothing changes: an ordinary film moves exactly as before, with no row capped
    ro2 = film.Runoff(s, pts)
    m1 = np.where(np.arange(s.n_patches) < 2, 1e-6, 0.0)
    ro2.transport(m1, lambda b: np.full(b.size, 1e-9), t_hat, LIQ.rho, s.areas, 0.5)
    assert ro2.capped_rows == 0
