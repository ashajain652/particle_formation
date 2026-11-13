import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import advect, correct_volumes, piece_volumes, pieces, point_triangle_distance, redistance


C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _sphere_mesh(h=0.125e-3):
    return build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))


def _matched_volumes(m, phi, label, n):
    """The volume [m3] of each piece of phi after a correction, the pieces labelled afresh and matched to the labels
    before it (`label`, n pieces) by overlap: each fresh piece goes to the old label holding most of its metal by cell,
    and the matching must be one to one. Measured afresh, so metal a correction grew in cells that held none before
    counts, and against a relative tolerance with no absolute one, so a wrong small piece fails (fix round 5)."""
    from reentry_model.levelset.cutcell import inside_volumes
    fresh, nf = pieces(m, phi)
    assert nf == n
    vol = inside_volumes(m, phi)
    both = (fresh >= 0) & (label >= 0)
    overlap = np.zeros((nf, n))
    np.add.at(overlap, (fresh[both], label[both]), vol[both])
    match = np.argmax(overlap, axis=1)
    assert sorted(match) == list(range(n))
    out = np.zeros(n)
    out[match] = piece_volumes(m, phi, fresh, nf)
    return out


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
    np.testing.assert_allclose(_matched_volumes(m, phi2, label, n), target, rtol=1e-12, atol=0.0)
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
    phi2, _ = correct_volumes(m, phi, label, targets)
    np.testing.assert_allclose(_matched_volumes(m, phi2, label, n), targets, rtol=1e-12, atol=0.0)


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
    phi2, _ = correct_volumes(m, phi, label, targets)
    np.testing.assert_allclose(_matched_volumes(m, phi2, label, n), targets, rtol=1e-12, atol=0.0)


def test_two_pieces_whose_dilations_meet_hold_the_total_volume_exactly():
    """Two spheres 0.25 mm apart on 0.125 mm cells, one grown 5 % and the other shrunk as much: the cells between them
    lie in both pieces' dilations, so the metal a shift grows there was counted for both and the total fell 3e-5 short;
    counted once, with the pieces corrected sweep after sweep, the total and each piece hold their targets (Task 9, fix
    round 4)."""
    from reentry_model.levelset.cutcell import inside_volumes
    m = build(BoxBase.around((0, 0, 0), (4e-3, 2e-3, 2e-3), 0.125e-3))
    x = m.points
    a = np.linalg.norm(x - [1.275e-3, 1e-3, 1e-3], axis=1) - 0.6e-3
    b = np.linalg.norm(x - [2.725e-3, 1e-3, 1e-3], axis=1) - 0.6e-3
    phi = redistance(m, np.minimum(a, b))
    label, n = pieces(m, phi)
    assert n == 2
    targets = piece_volumes(m, phi, label, n) * np.array([1.05, 1.0 / 1.05])
    phi2, _ = correct_volumes(m, phi, label, targets)
    assert inside_volumes(m, phi2).sum() == pytest.approx(targets.sum(), rel=1e-13, abs=0.0)
    np.testing.assert_allclose(_matched_volumes(m, phi2, label, n), targets, rtol=1e-12, atol=0.0)


def test_a_volume_correction_that_does_not_converge_says_where_it_stood():
    """Stage 3 Task 11's review, M6: a piece whose correction does not converge raises with what is needed to diagnose
    it -- its target, the volume it held, its last shift and bracket, the iterations, its claimed nodes and cells."""
    m = _sphere_mesh(0.25e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    label, n = pieces(m, phi)
    with pytest.raises(RuntimeError, match=r"piece 0 did not converge: target .* held .* bracket .* 1 iterations .* "
                                           r"nodes and \d+ cells claimed"):
        correct_volumes(m, phi, label, [1.2 * 4 / 3 * np.pi * 1.1e-3 ** 3], max_iter=1)


def test_a_speck_grown_to_its_target_converges_after_a_first_step_into_saturation():
    """Stage 3 Task 11's fix round 1, the review's concern (b): a speck of metal a base move makes beside the body, with a
    point's hand-down in it, must grow some 1e5 times. From its tiny patch area the first Newton step overshoots far
    beyond the nodes it may shift, where the volume stops responding (every cell it claimed full) while the cut cells at
    the claim's edge keep a patch area, so each Newton step from there moved the shift by only that area's share and the
    correction crawled for 60 iterations without reaching the root (the q6p box at 1.0 to 1.5 s: a 2e-17 m3 speck asked
    to hold 1.7e-12 m3; the shoulder at 1.5 to 2.0 s likewise). Where a Newton step stays in the bracket but fails to
    halve the step before the last (rtsafe's test), the shift is now found by Brent's method on the bracket, so the
    speck reaches its target to round-off and the body keeps its own (a step leaving the bracket is bisected, as before)."""
    m = _sphere_mesh(0.125e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 0.8e-3)
    speck = np.argmin(np.linalg.norm(m.points - (C + [1.25e-3, 0.0, 0.0]), axis=1))
    phi[speck] = -1e-6                                   # one node just inside: a speck of metal 0.45 mm off the body
    label, n = pieces(m, phi)
    assert n == 2
    v = piece_volumes(m, phi, label, n)
    targets = np.array([v[0], 1.0e-12])
    print("speck: {:.2e} m3 asked to hold {:.2e} m3".format(v[1], targets[1]))
    assert targets[1] > 1e5 * v[1]
    phi2, _ = correct_volumes(m, phi, label, targets)
    np.testing.assert_allclose(_matched_volumes(m, phi2, label, n), targets, rtol=1e-12, atol=0.0)


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
