"""Spreading the re-frozen deposit (sub-plan 22 Task 5c; Asha's idea and decisions of 2026-11-06,
docs/superpowers/specs/2026-11-06-deposit-spreading-design.md).

Both skin flights of 2026-11-05 stopped because a corner's move could no longer be absorbed by the 3D mesh, at the toe of a
ring of re-frozen metal travelling outward ahead of the flow: at the 100 mm sphere's worst edge the base steps about 1.4 mm
across a triangle about 1.5 mm wide. A 2 mm triangle cannot carry a ridge narrower than itself, so where the base as it
would stand after this part's move folds past a limit angle, part of the metal a corner is handing down goes to its lower
neighbouring corners instead: the angle-of-repose relaxation used for loose material on terrain, applied to the deposit
(`SkinBody._move_base` calls `spread_hand_down` after `skin.hand_down`). It moves solid metal sideways with no physical
mechanism, so the deposit's shape at the 2 mm scale is then set by the limit, and the sensitivity to it has to be measured.

Only handed metal moves: draws and absorption stay where they are, and the skins, which have already given the metal up,
do not change. Each pass finds the edges of the base as it would stand that fold past the limit. Across a ridge along an edge
(a, b), a and b give to the corners opposite it, c and e; across a valley, c and e give to a and b. Each of those four pairs
wants the share SPREAD_SHARE of the fold's excess over its angle, of the metal in the donor's dual cell above the receiver.
No donor falls below what balances its draw, so its move stays within the part's allowance, and no receiver rises above that
allowance: the part's move ratio holds, and a hand-down still never splits a part (spec 6.6, amended 2026-10-27). The metal
keeps the donor's mean enthalpy. SPREAD_PASSES passes a part at most, so metal travels at most a ring or two of 2 mm
neighbours in one part. What is left wanting to move after the last pass, and could leave its donor, is reported, as
`held_mass` reports what the tolerance holds. A limit that does not bind returns the very arrays it was given, so a run under
it is the run without it to the last bit."""
from dataclasses import dataclass

import numpy as np

SPREAD_PASSES = 4       # passes a part at most: metal travels at most a ring or two of 2 mm neighbours in one part
SPREAD_SHARE = 0.25     # of the fold's excess over a donor's four pairs
SPENT = 1e-14           # of a donor's hand-down: round-off of a donor that gave all it had, kept as exactly nothing


@dataclass
class Spread:
    handed: np.ndarray    # kg per column, the hand-down after spreading; the argument itself if the limit does not bind
    enthalpy: np.ndarray  # J per column, likewise
    moved: float          # kg moved between columns
    unplaced: float       # kg still wanting to move after the last pass
    over: int             # edges over the limit on the first estimate


def fold_angles(columns, points):
    """Per edge of the columns' triangulation with corners at `points` [m, (n, 3)]: the angle [deg] between its two
    triangles' normals, the measure of `movemesh.max_dihedral` and so of `surface_fold`, zero on an open boundary; and
    whether it is convex, a ridge along the edge, rather than a valley."""
    p = np.asarray(points, dtype=float)
    a, b, c = (p[columns.faces[:, v]] for v in range(3))
    unit = np.cross(b - a, c - a)
    unit /= np.linalg.norm(unit, axis=1)[:, None]
    angle, convex = np.zeros(len(columns.edges)), np.zeros(len(columns.edges), dtype=bool)
    two = np.flatnonzero(columns._tri[:, 1] >= 0)
    first, second = unit[columns._tri[two, 0]], unit[columns._tri[two, 1]]
    angle[two] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", first, second), -1.0, 1.0)))
    beyond = p[columns._opp[two, 1]] - p[columns.edges[two, 0]]         # the second triangle's far corner, under the
    convex[two] = np.einsum("ij,ij->i", beyond, first) < 0.0             # first's plane where the edge is a ridge
    return angle, convex


def _pairs(columns, q, normal, area, rho, angle, convex, over, limit_deg):
    """The donors, the receivers and the mass [kg] each pair wants to move, on the base at `q`, for the edges `over` the
    limit (an open boundary's edge has no fold, so no -1 corner comes in): the four pairs an edge makes, those whose donor
    stands above its receiver along the donor's normal."""
    ends, opposite = columns.edges[over].T, columns._opp[over].T
    high, low = np.where(convex[over], ends, opposite), np.where(convex[over], opposite, ends)
    donor = np.concatenate([high[0], high[0], high[1], high[1]])
    receiver = np.concatenate([low[0], low[1], low[0], low[1]])
    dz = np.einsum("ij,ij->i", q[donor] - q[receiver], normal[donor])
    above = dz > 0.0
    fold = np.tile(angle[over], 4)[above]
    want = SPREAD_SHARE * (fold - limit_deg) / fold * rho * area[donor[above]] * dz[above]
    return donor[above], receiver[above], want


def _within(want, who, room, n):
    """`want` [kg per pair] scaled so that what the pairs of each column `who` add up to stays within its `room` [kg]."""
    asked = np.bincount(who, want, n)
    return want * np.minimum(room / np.where(asked > 0.0, asked, 1.0), 1.0)[who]


def spread_hand_down(columns, points, normal, area, rho, drawn, handed, enthalpy, cap, limit_deg, passes=SPREAD_PASSES):
    """Spread the part's hand-down (`handed` [kg] and `enthalpy` [J] per column) over the corners at `points` [m], outward
    `normal`, dual `area` [m2], where the base as it would stand folds past `limit_deg` (module docstring): a Spread. `drawn`
    [kg] is the part's draws and `cap` [kg] its allowance per column (`SkinBody._hand_cap`). The base's estimate leaves out
    the part's absorption and smoothing. The sums of the hand-down and its enthalpy are kept to round-off."""
    n, first, moved = len(handed), 0, 0.0
    h, H = handed, enthalpy
    floor = np.maximum(drawn - cap, 0.0)                 # a donor keeps what balances its draw: the part's move ratio holds
    for k in range(passes + 1):
        q = points - normal * ((drawn - h) / (rho * area))[:, None]     # the base as it would stand
        angle, convex = fold_angles(columns, q)
        over = np.flatnonzero(angle > limit_deg)
        if k == 0:
            if not over.size:
                return Spread(handed, enthalpy, 0.0, 0.0, 0)
            first = int(over.size)
        if not over.size:
            return Spread(h, H, moved, 0.0, first)
        donor, receiver, want = _pairs(columns, q, normal, area, rho, angle, convex, over, limit_deg)
        want = _within(want, donor, np.maximum(h - floor, 0.0), n)
        if k == passes:
            return Spread(h, H, moved, float(want.sum()), first)
        want = _within(want, receiver, np.maximum(cap - h, 0.0), n)
        going = want > 0.0
        donor, receiver, want = donor[going], receiver[going], want[going]
        unit = H[donor] / h[donor]                       # the metal keeps the donor's mean enthalpy
        gave = np.bincount(donor, want, n)
        rest = 1.0 - gave / np.where(h > 0.0, h, 1.0)
        rest = np.where(rest < SPENT, 0.0, rest)         # all it had: exactly nothing, or the mean enthalpy the moving
        h = h * rest + np.bincount(receiver, want, n)    # mesh books its sweep at would be round-off over round-off
        H = H * rest + np.bincount(receiver, want * unit, n)
        moved += float(want.sum())
