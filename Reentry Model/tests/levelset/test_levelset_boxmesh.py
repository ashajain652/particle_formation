import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, base_mesh, build

H0 = 1.0e-3


def _faces(mesh):
    T = mesh.tets
    f = np.vstack([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]])
    f.sort(axis=1)
    uniq, counts = np.unique(f, axis=0, return_counts=True)
    return uniq, counts


def _area(pts, f):
    p = pts[f]
    return 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)


def _check_conforming(mesh):
    uniq, counts = _faces(mesh)
    assert counts.max() <= 2
    lo, hi = np.asarray(mesh.base.lo), np.asarray(mesh.base.hi)
    box_area = 2 * sum((hi - lo)[i] * (hi - lo)[j] for i, j in ((0, 1), (0, 2), (1, 2)))
    assert _area(mesh.points, uniq[counts == 1]).sum() == pytest.approx(box_area, rel=1e-12)
    assert mesh.volumes.sum() == pytest.approx(np.prod(hi - lo), rel=1e-12)
    assert (mesh.volumes > 0).all()


def test_the_base_grid_is_six_kuhn_tetrahedra_per_cube_and_conforming():
    base = BoxBase.around((0, 0, 0), (3e-3, 2e-3, 2e-3), H0)
    m = base_mesh(base)
    assert len(m.tets) == 6 * 3 * 2 * 2
    assert len(m.points) == 4 * 3 * 3
    _check_conforming(m)


def test_refinement_toward_a_point_stays_conforming_and_meets_the_target():
    base = BoxBase.around((0, 0, 0), (4e-3, 4e-3, 4e-3), H0)
    c = np.array([1.3e-3, 2.2e-3, 2.9e-3])
    target = lambda x, s: np.where(np.linalg.norm(x - c, axis=1) < 0.8e-3, 0.25e-3, np.inf)
    m = build(base, target)
    _check_conforming(m)
    near = np.linalg.norm(m.centroids - c, axis=1) < 0.8e-3
    assert near.any() and (m.size[near] <= 0.25e-3 + 1e-15).all()
    assert m.size.max() == H0


def test_achieved_size_is_the_largest_halving_not_above_the_target():
    base = BoxBase.around((0, 0, 0), (2e-3, 2e-3, 2e-3), H0)
    m = build(base, lambda x, s: np.full(len(x), 0.4e-3))
    assert set(np.unique(m.size)) == {0.25e-3}


def test_tetrahedra_stay_in_a_bounded_number_of_shapes():
    base = BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), H0)
    c = np.array([1.5e-3, 1.5e-3, 1.5e-3])
    m = build(base, lambda x, s: np.where(np.linalg.norm(x - c, axis=1) < 1.0e-3, 0.125e-3, np.inf))
    p = m.points[m.tets]
    edges = np.stack([np.linalg.norm(p[:, i] - p[:, j], axis=1) for i, j in
                      ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))], axis=1)
    quality = m.volumes / edges.max(axis=1) ** 3
    # Kuhn tetrahedra: volume / longest_edge^3 = (1/6) / 3**1.5; bisection keeps it within a fixed factor
    assert quality.min() > 0.25 * (1 / 6) / 3 ** 1.5
    assert len(np.unique(np.round(quality / quality.max(), 6))) <= 3


def test_rebuilding_with_the_same_targets_gives_the_same_mesh():
    base = BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), H0)
    tgt = lambda x, s: np.where(x[:, 0] < 1.2e-3, 0.25e-3, np.inf)
    a, b = build(base, tgt), build(base, tgt)
    assert np.array_equal(a.tets, b.tets) and np.array_equal(a.points, b.points)
