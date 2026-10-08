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
