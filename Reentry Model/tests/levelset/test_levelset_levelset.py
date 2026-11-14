import time

import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import inside_volumes, patches
from reentry_model.levelset.levelset import (MOVE_BAND_SHARE, MOVE_RTOL, advect, correct_volumes, move_by_volumes,
                                             piece_volumes, pieces, point_triangle_distance, redistance, swept_volumes)
from reentry_model.levelset.transfer import interpolate

from cases import (BOX_H0, DESIGN_SPACING, body_on, box_corners, box_sdf, edge_length, polygon_sdf, rounding_volume,
                   stage3_box, volume_between)


C = np.array([1.5e-3, 1.5e-3, 1.5e-3])
SHARP = 0.1                # of the rounding volume: an edge is sharp within this (plan 03b, planning ruling 11)
BODY_BAND = 2.0 * (4e-3 + BOX_H0)   # m: the box body's distance-reset band (skinbody.BAND_FACTOR (fine_margin + h0))


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


# The moves of the gas surface by each group's exact volume (design section 21.6; plan 03b Task 2, planning rulings 2 and
# 11, the controller's ruling P1). Each body is a level set phi that is an exact signed distance on the body's mesh, with
# markers seeded as the body seeds them. Each group is asked for the volume it sweeps between phi and an ideal shape
# (`swept_volumes`, by the groups move_by_volumes divides the cells into), so the move must find that shape from the
# groups' volumes alone. The ideal shape's exact signed distance on the same nodes is the reference (ruling P1): the P1
# grid's own cut at every edge is then on both sides. Sharp is judged against rounding (ruling 11): within a tenth of
# what rounding the edges to the move's depth would change, (1 - pi/4) d^2 per unit length, both in volume and in the
# volume between the moved contour and the ideal one.


def _move_to(body, ideal, band=BODY_BAND, factor=1.0):
    """phi moved by the volume each group sweeps between phi and the level set `ideal` (times `factor`, one or per
    group). Returns (the Moved, the volumes asked [m3], the move's seconds)."""
    dV = factor * swept_volumes(body.mesh, body.phi, ideal, body.x, body.normal)
    t0 = time.perf_counter()
    moved = move_by_volumes(body.mesh, body.phi, body.x, body.normal, dV, band=band)
    return moved, dV, time.perf_counter() - t0


def _judge(name, body, moved, took, ideal, analytic, tol):
    """The moved phi against the ideal shape on the same nodes: every group's miss within MOVE_RTOL, the volume and the
    volume between the contours within `tol` [m3], and the groups' swept volumes summing to phi's change."""
    V0 = inside_volumes(body.mesh, body.phi).sum()
    V = inside_volumes(body.mesh, moved.phi).sum()
    ref = inside_volumes(body.mesh, ideal).sum()
    between = volume_between(body.mesh, moved.phi, ideal)
    print("{}: V {:.4f} mm3, grid reference {:.4f}, analytic {:.4f}, tolerance {:.4f}; between {:.4f} mm3; {} round(s), "
          "largest miss {:.2e}, {:.1f} s".format(name, V * 1e9, ref * 1e9, analytic * 1e9, tol * 1e9, between * 1e9,
                                                 moved.rounds, moved.miss.max(), took))
    assert moved.miss.max() <= MOVE_RTOL
    assert abs(V - ref) <= tol
    assert between <= tol
    assert abs(moved.swept.sum() - (V - V0)) <= 1e-12 * V0
    return V


def test_an_equal_recession_keeps_a_box_s_edges_sharp():
    """Design section 21.8 item 2: every group of the stage 3 box loses what a recession of 1.0 mm sweeps in it; the box
    recedes to the box shrunk by 1.0 mm with its convex edges sharp."""
    box, (lo, hi), d = stage3_box(), box_corners(), -1.0e-3
    ideal = box_sdf(lo - d, hi + d)(box.mesh.points)
    moved, _, took = _move_to(box, ideal)
    _judge("equal recession", box, moved, took, ideal, np.prod(hi - lo + 2.0 * d),
           SHARP * rounding_volume(d, edge_length(hi - lo + 2.0 * d)))
    assert moved.rounds == 1


def test_extra_loss_along_one_edge_rounds_it():
    """Design section 21.8 item 3: the groups within one marker spacing of one edge of the stage 3 box (its heated face's
    edge with the +y face) lose twice what a recession of 0.2 mm sweeps in them, the others that once. Every group sweeps
    its volume, and phi at the edge's midpoint rises by more than at the centre of any face; and by more than the
    sqrt(2) times 0.2 mm that the equal recession alone gives a sharp edge's midpoint.

    0.2 mm, because a group sweeps only the metal of the part of the mesh its marker owns: the columns at the box's
    corners narrow to nothing within about a millimetre, and twice what 0.5 mm sweeps in the two corner groups here is
    1.36 and 1.50 times all their parts hold (measured once); twice 0.2 mm's is at most 0.79 times."""
    box, (lo, hi), d = stage3_box(), box_corners(), -0.2e-3
    ideal = box_sdf(lo - d, hi + d)(box.mesh.points)
    near = np.hypot(box.x[:, 0] - lo[0], box.x[:, 1] - hi[1]) < DESIGN_SPACING
    moved, _, took = _move_to(box, ideal, factor=np.where(near, 2.0, 1.0))
    mid = 0.5 * (lo + hi)
    edge = np.array([[lo[0], hi[1], mid[2]]])
    faces = np.array([[lo[0], mid[1], mid[2]], [hi[0], mid[1], mid[2]], [mid[0], lo[1], mid[2]],
                      [mid[0], hi[1], mid[2]], [mid[0], mid[1], lo[2]], [mid[0], mid[1], hi[2]]])
    pts = np.concatenate([edge, faces])
    rise = interpolate(box.mesh, moved.phi, pts) - interpolate(box.mesh, box.phi, pts)
    print("loss at an edge: {} groups doubled; phi rises {:.4f} mm at the edge's midpoint, {} mm at the face centres; "
          "{} round(s), largest miss {:.2e}, {:.1f} s".format(near.sum(), rise[0] * 1e3, np.round(rise[1:] * 1e3, 4),
                                                              moved.rounds, moved.miss.max(), took))
    assert moved.miss.max() <= MOVE_RTOL
    assert rise[0] > rise[1:].max()
    assert rise[0] > np.sqrt(2.0) * abs(d)


def test_a_deposit_on_a_convex_edge_rounds_it_by_the_depth_added():
    """Design section 21.8 item 4: every group of the stage 3 box gains what a deposit 0.5 mm deep sweeps in it; the box
    grows to its offset by 0.5 mm with every edge rounded to 0.5 mm (the box's sum with a ball, its exact signed distance
    less 0.5 mm), short of the sharp offset by the rounding volume and the corners' (8 - 4 pi / 3) d^3."""
    box, (lo, hi), d = stage3_box(), box_corners(), 0.5e-3
    ideal = box.phi - d
    moved, _, took = _move_to(box, ideal)
    size, L = hi - lo, edge_length(hi - lo)
    area = 2.0 * (size[0] * size[1] + size[1] * size[2] + size[2] * size[0])
    rounded = np.prod(size) + area * d + 0.25 * np.pi * d * d * L + 4.0 * np.pi / 3.0 * d ** 3
    V = _judge("deposit", box, moved, took, ideal, rounded, SHARP * rounding_volume(d, L))
    sharp = inside_volumes(box.mesh, box_sdf(lo - d, hi + d)(box.mesh.points)).sum()
    print("  the sharp offset holds {:.4f} mm3 on the mesh ({:.4f} analytic)".format(
        sharp * 1e9, np.prod(size + 2.0 * d) * 1e9))
    assert sharp - V >= (1.0 - SHARP) * rounding_volume(d, L)


def test_a_receding_concave_corner_rounds_by_the_depth_removed():
    """Design section 21.8 item 4: an L-shaped prism (an 8 mm square less its 4 mm quarter, through the mesh along z) whose
    every group loses what a recession of 0.5 mm sweeps in it: its convex edges stay sharp and its concave edge rounds to
    0.5 mm (the L's exact signed distance plus 0.5 mm), so that it loses less than the mitred offset would, by the
    rounding volume along the concave edge."""
    OFF, d = 0.13e-3, -0.5e-3
    corners = OFF + np.array([[0.0, 0.0], [8.0, 0.0], [8.0, 4.0], [4.0, 4.0], [4.0, 8.0], [0.0, 8.0]]) * 1e-3
    base = BoxBase.around((-BOX_H0, -BOX_H0, -4e-3), (8e-3 + BOX_H0, 8e-3 + BOX_H0, 4e-3), BOX_H0)
    body = body_on(base, polygon_sdf(corners))
    ideal = body.phi - d
    moved, _, took = _move_to(body, ideal)
    length = base.n[2] * base.h0                       # the concave edge's length, the mesh's height
    perimeter = 32e-3
    mitre = (perimeter * abs(d) - 4.0 * d * d) * length    # five convex corners overlap, the concave one adds
    removed = (mitre - rounding_volume(d, length))
    V0 = inside_volumes(body.mesh, body.phi).sum()
    V = _judge("concave corner", body, moved, took, ideal, (64e-6 - 16e-6) * length - removed,
               SHARP * rounding_volume(d, length))
    print("  removed {:.4f} mm3 against {:.4f} on the grid, {:.4f} analytic (the mitre's {:.4f})".format(
        (V0 - V) * 1e9, (V0 - inside_volumes(body.mesh, ideal).sum()) * 1e9, removed * 1e9, mitre * 1e9))


def test_a_move_longer_than_the_band_is_taken_in_rounds():
    """Planning ruling 2: a move longer than MOVE_BAND_SHARE of the reset's band is taken in rounds, each within it. A
    flat face (phi = y - 0.13 mm, the metal below it filling the mesh) on 0.25 mm cells through the move, markers 2 mm
    apart on it: with a band of 1.0 mm, every group loses what a recession of 1.5 mm sweeps in it. The move is taken in
    three rounds, every group sweeps its volume, and the face lies where the recession puts it.

    A flat face, because the rounds are equal shares of each group's volume and each round divides the cells afresh,
    the markers left where the move began. On the stage 3 box (a band of 2.0 mm, a recession of 3.0 mm, measured once)
    the first two rounds met their volumes to 6e-9 and 2e-10, but in the third a corner group's part of the mesh held
    0.10 of the 0.70 mm3 asked of it: the corner columns narrow to nothing within about a millimetre. And the band must
    exceed each round's move by more than the cells the contour crosses, whose nodes must hold distances: on the box
    with a band of 1.0 mm on its 0.7 mm cells, nodes held at the band lay in the cells the second round's contour
    crossed, and its groups missed by up to 0.56."""
    off, d, band, fine = 0.13e-3, -1.5e-3, 1.0e-3, 0.25e-3
    base = BoxBase((-5e-3, -3e-3, -5e-3), 0.5e-3, (20, 8, 20))
    mesh = build(base, lambda p, s: np.where((p[:, 1] - off > d - 2.0 * fine - s) & (p[:, 1] - off < 2.0 * fine + s),
                                              fine, np.inf))
    phi = mesh.points[:, 1] - off
    g = np.arange(-5e-3 + 0.5 * DESIGN_SPACING, 5e-3, DESIGN_SPACING) + 0.07e-3     # off the mesh's node planes
    X, Z = np.meshgrid(g, g, indexing="ij")
    x = np.column_stack([X.ravel(), np.full(X.size, off), Z.ravel()])
    normal = np.tile([0.0, 1.0, 0.0], (len(x), 1))
    ideal = phi - d
    dV = swept_volumes(mesh, phi, ideal, x, normal)
    t0 = time.perf_counter()
    moved = move_by_volumes(mesh, phi, x, normal, dV, band=band)
    took = time.perf_counter() - t0
    V, ref = inside_volumes(mesh, moved.phi).sum(), inside_volumes(mesh, ideal).sum()
    print("long move: {} rounds, largest miss {:.2e}, V - V_ideal {:.2e} of it, between {:.2e} mm3, {:.1f} s".format(
        moved.rounds, moved.miss.max(), V / ref - 1.0, volume_between(mesh, moved.phi, ideal) * 1e9, took))
    assert abs(d) > band * MOVE_BAND_SHARE
    assert moved.rounds == 3
    assert moved.miss.max() <= MOVE_RTOL
    assert abs(V / ref - 1.0) <= MOVE_RTOL
    assert volume_between(mesh, moved.phi, ideal) <= MOVE_RTOL * abs(dV.sum())


def test_no_change_of_volume_leaves_phi_as_it_was():
    """Design section 21.2: no construction moves the gas surface where no mass changed."""
    m = _sphere_mesh(0.25e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    x = C + 1.1e-3 * np.eye(3)
    moved = move_by_volumes(m, phi, x, np.eye(3), np.zeros(3))
    assert np.array_equal(moved.phi, phi) and moved.rounds == 0
    assert not moved.swept.any() and not moved.miss.any()
