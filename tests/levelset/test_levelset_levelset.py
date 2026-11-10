import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import advect, correct_volumes, piece_volumes, pieces, point_triangle_distance, redistance


C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _sphere_mesh(h=0.125e-3):
    return build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))


def test_point_triangle_distance_in_every_region():
    a, b, c = np.array([[0., 0, 0]]), np.array([[1., 0, 0]]), np.array([[0., 1, 0]])
    p = np.array([[0.2, 0.2, 0.5], [-1, -1, 0], [2, 0, 0], [0.5, -1, 0], [1, 1, 0], [-1, 0.5, 0]])
    want = [0.5, np.sqrt(2), 1.0, 1.0, np.sqrt(2) / 2, 1.0]
    got = point_triangle_distance(p, *(np.repeat(v, len(p), axis=0) for v in (a, b, c)))
    assert np.allclose(got, want)


def test_redistancing_a_scaled_sphere_field_recovers_the_distance():
    m = _sphere_mesh()
    r = np.linalg.norm(m.points - C, axis=1)
    d = redistance(m, 7.0 * (r - 1.1e-3))
    assert np.abs(d - (r - 1.1e-3)).max() < 0.2 * 0.125e-3


def test_a_prescribed_radial_recession_gives_the_analytic_sphere():
    m = _sphere_mesh()
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    v = 1e-3                                     # 1 mm/s inward
    inward = lambda x: -v * (x - C) / np.maximum(np.linalg.norm(x - C, axis=1, keepdims=True), 1e-12)
    for _ in range(6):
        phi = advect(m, phi, inward, 0.05)       # 6 x 0.05 mm
    P = patches(m, phi)
    radius = np.linalg.norm(P.centroid - C, axis=1)
    assert np.abs(radius - 0.8e-3).max() < 0.125e-3


def test_the_volume_correction_is_exact_to_round_off():
    m = _sphere_mesh()
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    label, n = pieces(m, phi)
    target = [4 / 3 * np.pi * 1.1e-3 ** 3]
    phi2, worst = correct_volumes(m, phi, label, target)
    assert n == 1
    assert piece_volumes(m, phi2, label, n)[0] == pytest.approx(target[0], rel=1e-13)
    assert 0 < worst < 0.01


def test_a_dumbbell_pinches_into_two_pieces_with_exact_volumes():
    m = build(BoxBase.around((0, 0, 0), (4e-3, 2e-3, 2e-3), 0.125e-3))
    c1, c2 = np.array([1.1e-3, 1e-3, 1e-3]), np.array([2.9e-3, 1e-3, 1e-3])
    x = m.points
    s1 = np.linalg.norm(x - c1, axis=1) - 0.7e-3
    s2 = np.linalg.norm(x - c2, axis=1) - 0.7e-3
    rod = np.maximum(np.linalg.norm(x[:, 1:] - 1e-3, axis=1) - 0.2e-3, np.abs(x[:, 0] - 2e-3) - 0.9e-3)
    phi = redistance(m, np.minimum(np.minimum(s1, s2), rod))
    assert pieces(m, phi)[1] == 1
    phi = redistance(m, phi + 0.3e-3)            # the whole surface recedes 0.3 mm: the 0.2 mm neck vanishes
    label, n = pieces(m, phi)
    assert n == 2
    targets = [4 / 3 * np.pi * 0.4e-3 ** 3] * 2
    phi, _ = correct_volumes(m, phi, label, targets)
    assert np.allclose(piece_volumes(m, phi, label, n), targets, rtol=1e-13)


def test_two_pieces_three_cells_apart_are_corrected_independently():
    m = build(BoxBase.around((0, 0, 0), (4e-3, 2e-3, 2e-3), 0.125e-3))
    x = m.points
    a = np.linalg.norm(x - [1.2e-3, 1e-3, 1e-3], axis=1) - 0.6e-3
    b = np.linalg.norm(x - [2.8e-3, 1e-3, 1e-3], axis=1) - 0.6e-3        # a 0.4 mm gap: three cells, inside two rings
    # (pieces closer than about two cells share faces of metal-holding cells and count as one piece)
    phi = redistance(m, np.minimum(a, b))
    label, n = pieces(m, phi)
    assert n == 2
    targets = [4 / 3 * np.pi * 0.62e-3 ** 3, 4 / 3 * np.pi * 0.58e-3 ** 3]
    phi, _ = correct_volumes(m, phi, label, targets)
    assert np.allclose(piece_volumes(m, phi, label, n), targets, rtol=1e-13)


def test_a_zero_volume_target_is_refused():
    m = _sphere_mesh(0.25e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    label, n = pieces(m, phi)
    with pytest.raises(ValueError, match="piece 0"):
        correct_volumes(m, phi, label, [0.0])


def test_departure_points_outside_the_box_give_finite_values():
    m = _sphere_mesh(0.25e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    out = advect(m, phi, lambda x: np.tile([2e-2, 0.0, 0.0], (len(x), 1)), 0.05)   # 1 mm shift
    assert np.isfinite(out).all()
