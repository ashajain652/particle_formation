"""spreading.py: the metal being handed down spread to lower neighbouring corners where the base would fold past a limit
(sub-plan 22 Task 5c; Asha 2026-11-06), on a small plate of near-equilateral 1 mm triangles, six neighbours to a corner as
on the production surface."""
import numpy as np
import pytest

from reentry_model import columns, movemesh, spreading

RHO = 2813.0                    # kg/m3
EDGE = 1e-3                     # m: the plate's triangles' edge
SPECIFIC = 5e5                  # J/kg: the enthalpy of the metal handed down, unless a test says otherwise


def plate(nx=25, ny=25, edge=EDGE, jitter=0.0):
    """A flat plate in z = 0: `ny` rows of `nx` corners, every other row shifted half an edge and each corner up to `jitter`
    edges off its place in x and y, as a `columns.Columns` (its triangles outward, +z), with `row(j)`, the corners of row j."""
    x, y = np.meshgrid(np.arange(nx) * edge, np.arange(ny) * edge * np.sqrt(3.0) / 2.0, indexing="ij")
    x = x + 0.5 * edge * (np.arange(ny) % 2)[None, :]
    points = np.c_[x.ravel(), y.ravel(), np.zeros(nx * ny)]
    points[:, :2] += jitter * edge * np.random.default_rng(1).uniform(-1.0, 1.0, (nx * ny, 2))
    at = lambda i, j: i * ny + j
    faces = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            a, b, c, d = at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1)
            faces += [[a, b, d], [a, d, c]] if j % 2 == 0 else [[a, b, c], [b, d, c]]
    col = columns.Columns(points, np.array(faces))
    col.row = lambda j: np.array([at(i, j) for i in range(nx)])
    return col


def raised(col, rows, heights):
    """The plate with each of the rows in `rows` lifted to its height in `heights` [m]."""
    p = col.points.copy()
    for j, z in zip(rows, heights):
        p[col.row(j), 2] += z
    col.update(p)
    return col


def spread(col, handed, limit=60.0, drawn=None, cap=None, enthalpy=None, passes=spreading.SPREAD_PASSES):
    drawn = np.zeros(col.n) if drawn is None else drawn
    allowance = 0.25 * col.edge_length.mean() * RHO * col.area
    H = SPECIFIC * handed if enthalpy is None else enthalpy
    return spreading.spread_hand_down(col, col.points, col.normal, col.area, RHO, drawn, handed, H,
                                      allowance if cap is None else cap, limit, passes=passes)


def test_the_fold_is_the_surface_s_measure_with_a_ridge_convex_and_a_valley_not():
    """`fold_angles` reads the same angle the moving mesh's `surface_fold` reports, and says which way the edge folds; an
    edge on an open boundary has one triangle and folds nowhere."""
    col = raised(plate(), [12], [0.5 * EDGE])
    angle, convex = spreading.fold_angles(col, col.points)
    assert angle.max() == pytest.approx(movemesh.max_dihedral(col.points, col.faces), rel=1e-12)
    crest, toe = np.isin(col.edges, col.row(12)).all(axis=1), np.isin(col.edges, col.row(11)).all(axis=1)
    assert angle[crest][2:-2].min() == pytest.approx(60.0, abs=1e-9) and convex[crest].all()     # two 30 degree faces
    assert angle[toe][2:-2].min() == pytest.approx(30.0, abs=1e-9) and not convex[toe].any()
    assert (angle[col._tri[:, 1] < 0] == 0.0).all() and not convex[col._tri[:, 1] < 0].any()


def test_a_limit_that_does_not_bind_returns_the_hand_down_itself():
    """On a flat plate nothing folds, so the hand-down and its enthalpy come back as the very arrays given and nothing has
    moved: a limit that does not bind changes nothing to the last bit."""
    col = plate()
    handed = np.full(col.n, 1e-7)
    H = SPECIFIC * handed
    s = spread(col, handed, enthalpy=H)
    assert s.handed is handed and s.enthalpy is H
    assert (s.moved, s.unplaced, s.over) == (0.0, 0.0, 0) and not s.landed.any() and s.landed.shape == (col.n,)


def test_a_ridge_hands_its_metal_to_its_lower_neighbours_and_keeps_the_sums():
    """A ridge 1 mm high along a row of corners, folded about 98 degrees against a 60 degree limit, handing down a quarter
    of the move allowance, with its own enthalpy per column: the sums of the mass and the enthalpy are kept to round-off,
    every crest corner hands down less than it would have, the metal lands on corners lower than the crest and on no other,
    and what lands first on a row that handed nothing carries the mean enthalpy of the crest's metal."""
    col = raised(plate(), [12], [1e-3])
    crest = col.row(12)
    handed = np.zeros(col.n)
    handed[crest] = 0.25 * 0.25 * EDGE * RHO * col.area[crest]
    unit = np.full(col.n, SPECIFIC)
    unit[crest] = SPECIFIC * (1.0 + 0.01 * np.arange(len(crest)) / len(crest))     # no two columns' metal alike
    H = unit * handed
    s = spread(col, handed, enthalpy=H)
    assert s.handed is not handed and s.over > 0 and s.moved > 0.0
    assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12) and s.enthalpy.sum() == pytest.approx(H.sum(), rel=1e-12)
    inner = crest[2:-2]                                            # clear of the plate's open ends
    assert (s.handed[inner] < handed[inner]).all()
    rose = s.handed > handed
    assert rose.any() and (col.points[rose, 2] < col.points[crest, 2].max()).all()
    assert (s.handed[col.row(11)[2:-2]] > 0.0).all() and (s.handed[col.row(13)[2:-2]] > 0.0).all()
    assert (s.handed >= 0.0).all()
    got = s.enthalpy[rose] / s.handed[rose]
    assert got.min() >= SPECIFIC * (1.0 - 1e-12) and got.max() <= unit[crest].max() * (1.0 + 1e-12)
    assert s.moved <= spreading.SPREAD_PASSES * handed.sum()


def test_a_valley_fold_gives_from_the_rim_to_the_floor():
    """A trench along a row, its rims 0.8 mm above its floor: the fold is the floor's, a valley of about 95 degrees (the rims
    fold about half that, within the 70 degree limit), so the rims hand their metal to the floor and the floor's corners,
    which handed down nothing, receive it."""
    col = plate()
    p = col.points.copy()
    p[:, 2] = 0.8e-3
    p[col.row(12), 2] = 0.0
    col.update(p)
    angle, convex = spreading.fold_angles(col, col.points)
    assert not convex[angle > 70.0].any() and (angle > 70.0).sum() > 0
    floor, rims = col.row(12), np.r_[col.row(11), col.row(13)]
    handed = np.zeros(col.n)
    handed[rims] = 0.25 * 0.25 * EDGE * RHO * col.area[rims]
    s = spread(col, handed, limit=70.0)
    assert (s.handed[floor[2:-2]] > 0.0).all() and handed[floor].sum() == 0.0
    assert s.handed[rims].sum() < handed[rims].sum()
    assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12)
    rose = np.flatnonzero(s.handed > handed)
    assert set(rose) <= set(floor)


def test_a_corner_that_receives_and_gives_more_onward_is_marked_as_landed():
    """A two-step cascade down a staircase whose 49 degree toe is past a 40 degree limit: the crest gives its hotter metal to the
    flank below it, and the flank, which holds far more, gives more onward than it took, so it ends handing down less than it
    started with and its total alone would not show that anything arrived. It carries the crest's heat all the same, so `landed`
    marks it, and the toe it gave to; nothing arrives at the crest (Finding 1 of Task 5c's review: comparing the totals missed
    these corners, and their neighbours' heat crossed the base through `h_outward`)."""
    col = raised(plate(), [10, 11, 12, 13, 14], [0.0, 1e-3, 2e-3, 1e-3, 0.0])
    crest, flank, toe = col.row(12)[2:-2], col.row(11)[2:-2], col.row(10)[2:-2]
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    handed = np.zeros(col.n)
    handed[col.row(12)], handed[col.row(11)], handed[col.row(13)] = (0.1 * cap[col.row(12)], 0.8 * cap[col.row(11)],
                                                                     0.8 * cap[col.row(13)])
    unit = np.full(col.n, SPECIFIC)
    unit[col.row(12)] = 1.1 * SPECIFIC
    s = spread(col, handed, limit=40.0, cap=cap, enthalpy=unit * handed)
    assert (s.handed[flank] < handed[flank]).all() and (s.handed[flank] > 0.0).all()      # net donors, with metal left
    assert s.landed[flank].all() and s.landed[toe].all() and not s.landed[crest].any()
    assert s.landed.dtype == bool and s.landed.shape == (col.n,)
    mean = s.enthalpy[flank] / s.handed[flank]
    assert (mean > SPECIFIC * (1.0 + 1e-6)).all() and (mean < 1.1 * SPECIFIC).all()       # the crest's heat came with it
    assert (s.landed[s.handed > handed]).all()                                            # whatever rose, landed


def growing_ridge(limit, parts=200, nx=37, ny=37):
    """A jittered plate on which the middle of the middle row, 21 corners stopping short of the rim, hands down a tenth of
    the move allowance every part: the spreading pass (with `limit`, or none), then the body's move, `SkinBody._move_base`'s
    -- every corner out along its normal by what lands on it, with the tangential smoothing made volume-neutral as
    `MovingMesh.neutral_smoothing` makes it -- the plate's open rim held still, as a closed surface has no rim for the
    smoothing to draw in. Returns the largest fold [deg] after every part and the plate."""
    col = plate(nx, ny, jitter=0.1)
    crest = col.row(ny // 2)[8:-8]
    rim = np.unique(col.edges[col._tri[:, 1] < 0].ravel())
    edges = movemesh.surface_edges(col.faces)

    def held(u):
        u[rim] = 0.0
        return u
    folds = []
    for _ in range(parts):
        p, normal, area = col.points, col.normal, col.area
        allowance = 0.25 * col.edge_length.mean() * RHO * area
        handed = np.zeros(col.n)
        handed[crest] = 0.1 * 0.25 * EDGE * RHO * area[crest]
        landing = handed
        if limit is not None:
            s = spreading.spread_hand_down(col, p, normal, area, RHO, np.zeros(col.n), handed, SPECIFIC * handed, allowance,
                                           limit)
            assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12) and (s.handed <= allowance * (1.0 + 1e-12)).all()
            landing = s.handed
        u = held(movemesh.tangential_smoothing(p, col.faces, normal, edges=edges))
        for _ in range(movemesh.NEUTRAL_PASSES):
            u = held(u + normal * (movemesh.swept_volumes(p, p + u, col.faces, col.faces, col.n) / area)[:, None])
        col.update(p + u + held(normal * (landing / (RHO * area))[:, None]))
        folds.append(float(spreading.fold_angles(col, col.points)[0].max()))
    return np.array(folds), col


def test_repeated_parts_hold_a_growing_ridge_near_the_limit():
    """The test of the pass (Asha's idea, 2026-11-06): a crest handing down 0.025 of an edge every part for 200 parts, moved
    as the body moves its corners. With the pass off it passes 90 degrees, past spec 12's 75-95 degrees, by part 31 and
    goes on climbing; with a 60 degree limit the base stays within 5 degrees of the limit the whole way, not flattened below
    it either, the ridge widening instead. Measured 2026-11-06 on this plate: 60.03 degrees at most and 59.96 at least over
    the last 100 parts. SPREAD_SHARE stayed at the 0.25 the design proposed: from 0.0625 to 1.0 the largest fold is 60.03 to
    61.02, 0.25 the closest, so the share is not what holds the fold; a limit of 45, 75 or 90 degrees holds to within 0.5; one
    pass a part overshoots to 70.6, two to 60.4 and four, SPREAD_PASSES, to 60.03. At twice this crest's rate the fold holds
    at 60.7 for 140 parts, then a valley grows at the foot of the ridge's end, where only 4 of the 16 corners around it hold
    any handed metal and the smoothing has squeezed the triangles to 0.2 mm edges: only handed metal moves."""
    off = growing_ridge(None, parts=60)[0]
    assert off[:40].max() > 90.0 and off.max() > 120.0
    folds, col = growing_ridge(60.0)
    assert folds.max() <= 65.0
    assert folds[-100:].min() >= 55.0
    assert (col.points[:, 2] > 0.1e-3).sum() > 3 * 21             # it has widened: more than three rows of the crest's length


def test_a_receiving_corner_stays_within_its_move_allowance():
    """However much the ridge would give, no corner hands down more than a part's allowance: a corner's own hand-down plus
    what lands on it moves its base by at most the allowance, the cap of 2026-10-27, so a part is never split for it."""
    col = raised(plate(), [12], [2e-3])
    crest = col.row(12)
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    handed = np.zeros(col.n)
    handed[crest] = cap[crest]
    handed[col.row(11)] = 0.9 * cap[col.row(11)]                   # next to nothing left of the allowance
    s = spread(col, handed, cap=cap)
    assert s.moved > 0.0
    assert (s.handed <= cap * (1.0 + 1e-12)).all()
    assert s.handed[col.row(11)].max() > 0.9 * cap[col.row(11)].max()


def test_a_donor_keeps_what_balances_its_draw():
    """A corner drawing metal moves in by its draw less its hand-down, and that move too is within the allowance, or the part
    is split for it: a donor gives only down to its draw less the allowance, and a corner already at or under that gives
    nothing."""
    col = raised(plate(), [12], [2e-3])
    crest = col.row(12)
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    drawn, handed = np.zeros(col.n), np.zeros(col.n)
    drawn[crest], handed[crest] = 1.6 * cap[crest], cap[crest]     # moves in by 0.6 allowances, and may reach 1
    idle = crest[10]
    drawn[idle], handed[idle] = 2.0 * cap[idle], 0.5 * cap[idle]   # over the allowance already: it gives nothing
    s = spread(col, handed, drawn=drawn, cap=cap)
    floor = drawn - cap
    assert s.moved > 0.0
    others = crest[crest != idle]
    assert (s.handed[others] >= floor[others] * (1.0 - 1e-12)).all()
    assert (s.handed[others[2:-2]] < handed[others[2:-2]]).any()
    assert s.handed[idle] >= handed[idle]


def test_what_cannot_be_placed_stays_and_is_counted():
    """Lower corners with no room left: with every neighbour of the crest at its allowance nothing can land, so the crest
    keeps its metal and the mass that wanted to move is reported, as `held_mass` reports what the tolerance holds. With a
    little room, that much lands and the rest is reported."""
    col = raised(plate(), [12], [2e-3])
    crest, flanks = col.row(12), np.r_[col.row(11), col.row(13)]
    cap = 0.25 * col.edge_length.mean() * RHO * col.area
    handed = np.zeros(col.n)
    handed[crest], handed[flanks] = 0.5 * cap[crest], cap[flanks]
    # the flanks' toes are folded too and may give to the next rows, leaving room: a limit only the crest exceeds stops that
    q = col.points + col.normal * (handed / (RHO * col.area))[:, None]
    limit = float(spreading.fold_angles(col, q)[0].max()) - 5.0
    s = spread(col, handed, cap=cap, limit=limit, passes=1)
    assert s.moved == 0.0 and s.unplaced > 0.0 and s.over > 0
    np.testing.assert_array_equal(s.handed, handed)
    handed[flanks] = 0.9 * cap[flanks]
    room = (cap - handed)[flanks].sum()
    s = spread(col, handed, cap=cap, limit=limit, passes=1)
    assert 0.0 < s.moved <= room * (1.0 + 1e-12) and s.unplaced > 0.0


def test_a_donor_that_gives_all_it_has_hands_down_exactly_nothing():
    """A corner with a sliver of a hand-down on a ridge far past the limit gives it all: its hand-down and its enthalpy are
    exactly zero, not round-off, whose quotient would be the mean enthalpy the moving mesh books its sweep at, and the sums
    still hold."""
    col = raised(plate(), [12], [1e-3])
    crest = col.row(12)
    handed = np.zeros(col.n)
    handed[crest] = 1e-3 * 0.25 * EDGE * RHO * col.area[crest]
    H = SPECIFIC * handed
    s = spread(col, handed, enthalpy=H)
    inner = crest[2:-2]
    assert (s.handed[inner] == 0.0).all() and (s.enthalpy[inner] == 0.0).all() and not s.landed[inner].any()
    assert s.handed.sum() == pytest.approx(handed.sum(), rel=1e-12) and s.enthalpy.sum() == pytest.approx(H.sum(), rel=1e-12)
