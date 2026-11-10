import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.transfer import integral, interpolate, locate, project

H0 = 1.0e-3
BASE = BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), H0)
C = np.array([1.4e-3, 1.6e-3, 1.5e-3])


def _near(r, s):
    return lambda x, sz: np.where(np.linalg.norm(x - C, axis=1) < r, s, np.inf)


def test_locate_finds_the_containing_cell_with_nonnegative_barycentrics():
    m = build(BASE, _near(1e-3, 0.25e-3))
    rng = np.random.default_rng(1)
    pts = rng.uniform(0, 3e-3, size=(500, 3))
    cell, b = locate(m, pts)
    assert (b.min(axis=1) > -1e-10).all()
    assert np.allclose(np.einsum("ni,nij->nj", b, m.points[m.tets[cell]]), pts, atol=1e-15)


def test_a_linear_field_is_interpolated_exactly():
    m = build(BASE, _near(1e-3, 0.25e-3))
    f = lambda x: 3.0 + 2e3 * x[:, 0] - 5e3 * x[:, 1] + 7e3 * x[:, 2]
    pts = np.random.default_rng(2).uniform(0, 3e-3, size=(200, 3))
    assert np.allclose(interpolate(m, f(m.points), pts), f(pts), rtol=1e-12)


def test_projection_onto_a_refinement_reproduces_the_field_exactly():
    old = build(BASE, _near(0.8e-3, 0.5e-3))
    new = build(BASE, _near(1.2e-3, 0.25e-3))
    u = np.sin(1e3 * old.points[:, 0]) * np.cos(2e3 * old.points[:, 1])
    v = project(old, new, u)
    assert np.allclose(v, interpolate(old, u, new.points), atol=1e-10)


def test_projection_onto_a_coarsening_conserves_the_integral():
    old = build(BASE, _near(1.2e-3, 0.25e-3))
    new = build(BASE, _near(0.6e-3, 0.5e-3))
    u = 900.0 + 50.0 * np.exp(-((old.points - C) ** 2).sum(axis=1) / (0.5e-3) ** 2)
    v = project(old, new, u)
    assert integral(new, v) == pytest.approx(integral(old, u), rel=1e-12)
