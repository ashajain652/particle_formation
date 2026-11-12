"""Tests of the level set's runoff between markers (design section 20.6, stage 3 Task 6): the pushes, the momentum
through the depth with its response, the thickness and capillary terms implicit in the heights, the implicit upwind
transport across the markers' faces with the arrival rule, the landing conduction and the flow sub-steps' count. The
strip cases are the skyscrapers' (tests/test_reentry_model_skin.py: landing velocity ~1625, runoff onto cold metal
~2051; the pushed plate, tests/test_reentry_model_skinbody.py ~926), built here on straight strips of markers whose faces
are given directly, with no mesh."""
from types import SimpleNamespace

import numpy as np
import pytest

from reentry_model.levelset import runoff
from reentry_model.levelset.markers import Faces, Markers
from reentry_model.levelset.skin import Skins, SkinSettings
from reentry_model.shared.material import Material

from materials import LinearMeltMaterial

LIQUID, RIGID, COLD = 950.0, 880.0, 700.0     # K: all liquid; rigid (AA7075's T_rigid 895.1 K); the cold metal runoff
                                              # lands on (the skyscrapers' strip cases, test_reentry_model_skin.py)
CELL = 1.0e-5                                 # m: each stacked cell holds 10 um of solid metal's mass (their `stacked`)
SETTINGS = SkinSettings(cell=10e-6, fine_zone=0.4e-3, stub=0.4e-3)   # their CELLS_04 geometry
# the skyscrapers' PUSHES["flight"] (analysis/skin_plate.py): the 100 mm sphere's shoulder 45 degrees from the stagnation
# point at 93 s -- dp/ds [Pa/m], gas shear [Pa], deceleration [m/s2] -- copied as data
FLIGHT_ALPHA, FLIGHT_DP_DS, FLIGHT_TAU, FLIGHT_DECEL = 45.0, -3.424e5, 80.29, 29.77


def aa7075():
    return Material.from_drama_json("AA7075_scheil")


def grid_markers(nx, ny, h=1e-3):
    """Markers on a flat regular grid of nx by ny points h apart in the plane z = 0, normals +z, each with a square cell
    of area h^2 and a face of length h to each of its four grid neighbours (the exact Voronoi tessellation), built with
    no mesh. Returns the markers (faces set) and each marker's grid indices (ix, iy)."""
    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    ix, iy = ix.ravel(), iy.ravel()
    x = np.column_stack([ix * h, iy * h, np.zeros(ix.size)])
    n = len(x)
    idx = np.arange(n).reshape(nx, ny)
    i = np.r_[idx[:-1, :].ravel(), idx[:, :-1].ravel()]
    j = np.r_[idx[1:, :].ravel(), idx[:, 1:].ravel()]
    normal = x[j] - x[i]
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    mk = Markers(x, np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, h * h), np.zeros((n, 2)), np.arange(n),
                 np.full(n, -1), h)
    mk.faces = Faces(i, j, np.full(len(i), h), normal)
    return mk, ix, iy


def stacked(mk, temps, mat=None):
    """Skins on the markers, marker i holding a cell of 10 um of solid metal's mass at each temperature of temps[i], top
    first (the skyscrapers' `stacked`)."""
    mat = aa7075() if mat is None else mat
    sk = Skins(mat, mk.area, mk.kappa, SETTINGS)
    count = np.array([len(t) for t in temps])
    K = int(count.max())
    T = np.array([list(t) + [300.0] * (K - len(t)) for t in temps])
    m = np.where(np.arange(K)[None, :] < count[:, None], mat.rho * sk.area[:, None] * CELL, 0.0)
    sk.install(m, m * sk.props.enthalpy(T), count)
    mk.skins = sk
    return sk


def flight_forcing(mk, v_hat, dp_ds=0.0, tau=0.0, decel=0.0, sigma=0.8):
    """A Forcing on the markers: the surface flow's dp/ds and shear (as a SurfaceFlowResult gives them, one per marker),
    the deceleration along the flight direction v_hat (the acceleration the metal feels), surface tension sigma, phi's
    mean curvature from the markers' principal curvatures."""
    flow = SimpleNamespace(dp_ds=np.full(mk.n, dp_ds), tau=np.full(mk.n, tau))
    return runoff.Forcing.from_surface_flow(flow, decel * np.asarray(v_hat), v_hat, sigma, mk.kappa.sum(axis=1))


# the strip's flight direction: windward face +z, the meridian along +x (away from the stagnation point toward +x)
STRIP_V_HAT = np.array([-0.6, 0.0, 0.8])


def momentum_closes(sk, P0, book, rel):
    P1 = sk.totals()[2].sum(axis=0)
    impulse = (book.push + book.shear + book.wall + book.carried).sum(axis=0)
    size = sum(np.linalg.norm(v, axis=1).sum() for v in (book.push, book.shear, book.wall))
    return np.linalg.norm(P1 - P0 - impulse) <= rel * size


# -- the meridian and the pushes ---------------------------------------------------------------------------------------
def test_the_pushes_point_along_the_meridian_away_from_the_stagnation_point():
    """On the strip the meridian runs along +x: a pressure falling along it pushes toward +x, as the gas shear does; the
    deceleration pulls the metal back toward the nose (-x) and out of the face; a thickness gradient under that pull
    draws metal toward the thicker side; the capillary term vanishes on a flat, even layer."""
    mk, ix, _ = grid_markers(12, 1)
    sk = stacked(mk, [[LIQUID] * 20 + [RIGID] * 5] * mk.n)
    f = flight_forcing(mk, STRIP_V_HAT, dp_ds=-2e4, tau=20.0, decel=50.0)
    G, terms = runoff.pushes(mk, sk, f)
    rho_l = sk.rho_l
    np.testing.assert_allclose(terms["pressure"], np.tile([2e4, 0.0, 0.0], (mk.n, 1)), rtol=1e-12)
    np.testing.assert_allclose(terms["deceleration"], np.tile([-0.6 * 50.0 * rho_l, 0.0, 0.0], (mk.n, 1)), rtol=1e-12)
    assert np.abs(terms["thickness"]).max() < 1e-9 and np.abs(terms["capillary"]).max() < 1e-9
    np.testing.assert_allclose(G, sum(terms.values()))
    shear = runoff.along_meridian(mk.normal, STRIP_V_HAT, np.full(mk.n, 20.0))
    np.testing.assert_allclose(shear, np.tile([20.0, 0.0, 0.0], (mk.n, 1)), rtol=1e-12)


# -- the transport's arrival rule --------------------------------------------------------------------------------------
def test_landed_metal_flows_on_at_its_arrival_velocity():
    """The skyscrapers' case (test_reentry_model_skin.py ~1625): liquid on the first five markers of a strip of 1 mm
    markers moving at 0.5 m/s toward rigid metal, one 10 ms transport sub-step (a Courant number of 5 along the strip):
    the metal landing on the first rigid top flows on in the same solve over every rigid marker to the strip's end, each
    holding 5/6 of the last within 1 % -- the one-dimensional implicit upwind scheme's nu / (1 + nu) for a pulse entering
    markers at rest -- at the velocity it arrived with; the books close and no rigid cell changes."""
    mk, ix, _ = grid_markers(20, 1)
    sk = stacked(mk, [[LIQUID] * 10 + [RIGID] * 5 if x < 5 else [RIGID] * 15 for x in ix])
    u = np.array([0.5, 0.0, 0.0])
    sk.U = np.where((sk.valid & (sk.T() > 900.0))[..., None], u, 0.0)
    before = [a.copy() for a in (sk.m, sk.H, sk.count)]
    M0 = sk.m.sum()
    book = runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
    per = book.mass
    assert np.all(per[5:] > 0.0)                                         # 15 markers in one sub-step
    np.testing.assert_allclose(per[6:19] / per[5:18], 5.0 / 6.0, rtol=0.01)
    new = sk.count > before[2]
    assert np.array_equal(new, ix >= 5) and np.abs(sk.U[new, 0] - u).max() < 1e-12
    assert abs(book.mass.sum()) < 1e-16 * M0 and sk.m.sum() == pytest.approx(M0, rel=1e-14)
    for c in np.flatnonzero(new):
        for now, then in ((sk.m, before[0]), (sk.H, before[1])):
            assert np.array_equal(now[c, 1:sk.count[c]], then[c, :before[2][c]])


def test_nothing_passes_through_anchored_metal():
    """Anchored rigid cells neither leave nor receive whatever velocity they hold: a liquid layer running at 1 m/s onto a
    marker that is rigid through its depth lands on top of it as a new top cell; the receiver's anchored cells, and the
    donor's own anchored cells (given a velocity too), are untouched to the last bit. A floating crust's metal lands on
    top wherever it goes, even onto a flowing top."""
    mk, ix, _ = grid_markers(3, 1)
    sk = stacked(mk, [[LIQUID] * 10 + [RIGID] * 5, [RIGID] * 15, [LIQUID] * 10 + [RIGID] * 5])
    sk.U[:] = [1.0, 0.0, 0.0]                                       # anchored cells too: they must not move
    m0, H0 = sk.m.copy(), sk.H.copy()
    anch = sk.anchored()
    book = runoff.transport(sk, runoff.Sheet.of(mk), 1e-3)
    assert sk.count[1] == 16 and book.landed[1] and not book.landed[0]
    assert np.array_equal(sk.m[1, 1:16], m0[1, :15]) and np.array_equal(sk.H[1, 1:16], H0[1, :15])
    assert np.array_equal(sk.m[0, 10:15], m0[0, 10:15]) and np.array_equal(sk.H[0, 10:15], H0[0, 10:15])
    assert np.array_equal(sk.m[2, 10:15], m0[2, 10:15]) and anch[:, 10:].all()
    assert sk.m.sum() == pytest.approx(m0.sum(), rel=1e-14) and abs(book.mass.sum()) < 1e-16 * m0.sum()
    # a crust: rigid cells over flowing metal over the floor, moving at 1 m/s toward a flowing neighbour
    mk2, _, _ = grid_markers(2, 1)
    sk2 = stacked(mk2, [[RIGID] * 3 + [LIQUID] * 7 + [RIGID] * 5, [LIQUID] * 10 + [RIGID] * 5])
    sk2.U[0, :10] = [1.0, 0.0, 0.0]
    crust = sk2.crust()
    assert crust[0, :3].all()
    m1 = sk2.m.copy()
    book2 = runoff.transport(sk2, runoff.Sheet.of(mk2), 1e-3)
    assert book2.landed[1] and sk2.count[1] == 16                     # the crust's metal (and only it) lands on top
    assert sk2.T()[1, 0] == pytest.approx(RIGID, abs=1e-6)          # rigid metal, the crust's own
    assert np.array_equal(sk2.m[1, 11:16], m1[1, 10:15])


# -- the divergence-free flow on a sphere (Task 5's cell-closure ruling) ------------------------------------------------
def test_a_divergence_free_flow_does_not_pile_up_mass():
    """A uniform liquid skin on the 25 mm sphere's markers (2 mm cells), the cells' velocities held at a rigid rotation
    omega x x about an axis -- divergence-free on the sphere -- for 20 transport sub-steps of 1 ms at 0.5 m/s at the
    equator (10 mm, five markers' spacings): the total mass is exact and every skin's mass stays within PILE_BOUND of the
    mean. The bound: the implicit upwind scheme on a divergence-free field moves no mass where the cells close and the
    velocity is uniform across a cell; what remains is the Voronoi cells' non-closure (Task 5: median 1.5 %, max 2.6 %
    of the perimeter) times the Courant number per sub-step, accumulated, plus the upwind error of the velocity's
    variation over a cell (h / R about 0.08); see the report for the measured spread."""
    from reentry_model.levelset.boxmesh import BoxBase, build
    R, h0 = 25e-3, 2e-3
    mesh = build(BoxBase.around((-31e-3,) * 3, (31e-3,) * 3, h0))
    phi = np.linalg.norm(mesh.points, axis=1) - R
    mat = LinearMeltMaterial()
    mk = Markers.create(mesh, phi, None, material=mat, settings=SkinSettings(cell=20e-6))
    mk.tessellate(mesh, phi)
    sk = mk.skins
    sk.area[mk.skin_index] = mk.area
    sk.kappa[mk.skin_index] = mk.kappa
    sk.fill(mat.liquid.rho * sk.area * 0.2e-3, LIQUID)
    sheet = runoff.Sheet.of(mk)
    axis = np.array([0.3, -0.2, 0.93])
    axis /= np.linalg.norm(axis)
    omega = 0.5 / R * axis
    U = np.zeros((sk.n, 3))
    U[mk.skin_index] = np.cross(omega, mk.x)
    M0 = sk.m.sum()
    for _ in range(20):
        sk.U = np.where(sk.valid[..., None], U[:, None, :], 0.0)
        runoff.transport(sk, sheet, 1e-3)
        sk.regrid()
    M = sk.m.sum(axis=1) / sk.area
    spread = np.abs(M / M.mean() - 1.0).max()
    print("divergence-free rotation: largest skin mass per area off the mean by", spread)
    assert sk.m.sum() == pytest.approx(M0, rel=1e-13)
    assert spread < PILE_BOUND


PILE_BOUND = 0.01


# -- the landing conduction and runoff onto cold metal ------------------------------------------------------------------
def test_landing_conduction_keeps_the_books_and_touches_only_the_landed_blocks():
    """A film of 20 um at 950 K on 0.8 mm of rigid metal at 700 K, one landing conduction of 0.5 ms (the skyscrapers'
    case, ~2005): the block of the cells within 0.4 mm of the top is conducted with no heat across its top or bottom; the
    film freezes, the block keeps its energy to round-off, and the cells beneath it and the markers nothing landed on are
    untouched to the last bit."""
    mk, _, _ = grid_markers(4, 2)
    sk = stacked(mk, [[LIQUID] * 2 + [COLD] * 80] * mk.n)
    landed = np.arange(mk.n) < 4
    m0, H0 = sk.m.copy(), sk.H.copy()
    runoff.landing_conduct(sk, landed, 5e-4)
    dz = sk.dz()
    k = int(((np.cumsum(dz, axis=1) - dz) < runoff.LANDING_DEPTH)[0].sum())
    T = sk.T()
    assert 10 < k < 80 and np.all(T[landed, :2] < sk.material.T_rigid) and np.all(T[landed, :k] > COLD)
    np.testing.assert_allclose(sk.H.sum(axis=1), H0.sum(axis=1), rtol=1e-14, atol=0.0)
    assert np.array_equal(sk.H[landed, k + 1:], H0[landed, k + 1:]) and np.array_equal(sk.H[~landed], H0[~landed])
    assert np.array_equal(sk.m, m0)


def cold_runoff(landing=True, monkeypatch=None, cycles=10, tau=20.0, edge=10):
    """The skyscrapers' `cold_runoff` (~2033): liquid 0.2 mm deep (20 cells at 950 K over five at 700 K) on the first
    `edge` markers of a strip of 1 mm markers, rigid metal at 700 K (25 cells) beyond, sheared along the strip by tau Pa
    for one 10 ms heat sub-step of `cycles` flow sub-steps (no heat sub-step conduction), the landing conduction on or
    off. Returns the skins, the flow book, each marker's place along the strip and the start's (mass, enthalpy,
    momentum, per-marker mass)."""
    mk, ix, _ = grid_markers(40, 1)
    sk = stacked(mk, [[LIQUID] * 20 + [COLD] * 5 if x < edge else [COLD] * 25 for x in ix])
    if not landing:
        monkeypatch.setattr(runoff, "landing_conduct", lambda sk_, landed, dt, depth=None: None)
    f = flight_forcing(mk, STRIP_V_HAT, tau=tau)
    f.cycles = cycles
    start = (sk.m.sum(), sk.H.sum(), sk.totals()[2].sum(axis=0), sk.m.sum(axis=1))
    book = runoff.flow_substeps(sk, mk, f, 0.01)
    return sk, book, ix, start


def test_runoff_onto_cold_metal_freezes_where_it_lands_and_the_books_close(monkeypatch):
    """The skyscrapers' case (~2051): the sheared layer runs onto metal at 700 K for one 10 ms heat sub-step of ten 1 ms
    flow sub-steps. Without the landing conduction the landed metal coasts on over the cold tops and stays liquid:
    several markers past the edge hold metal above T_rigid. With it the metal freezes where it lands: only the edge's
    first marker may stay non-rigid, and the deposit's centroid past the edge is under 0.4 of the no-landing case's.
    Mass and enthalpy close to 1e-13 (no heat enters or leaves), the momentum to 1e-12 of the impulses."""
    sk, book, x, (M0, E0, P0, m0) = cold_runoff()
    T_rigid = sk.material.T_rigid
    hot = lambda s: sorted(set(x[(s.valid & (s.T() > T_rigid)).any(axis=1) & (x >= 10)]))
    gained = lambda s, m: np.maximum(s.m.sum(axis=1) - m, 0.0)[10:]
    centroid = lambda g: float((np.arange(len(g)) * g).sum() / g.sum())
    assert book.sub_steps == 10 and book.courant > 0.0
    assert hot(sk) in ([], [10])
    assert sk.m.sum() == pytest.approx(M0, rel=1e-13) and sk.H.sum() == pytest.approx(E0, rel=1e-13)
    assert momentum_closes(sk, P0, book, 1e-12)
    here = centroid(gained(sk, m0))
    free, book_f, _, (M1, E1, P1, m1) = cold_runoff(landing=False, monkeypatch=monkeypatch)
    there = centroid(gained(free, m1))
    print("runoff onto cold metal: hot past the edge", hot(sk), "against", hot(free), "; centroid", here, "against",
          there, "markers")
    assert len(hot(free)) >= 4 and here < 0.4 * there
    assert free.m.sum() == pytest.approx(M1, rel=1e-13) and momentum_closes(free, P1, book_f, 1e-12)


# -- the pushed face ---------------------------------------------------------------------------------------------------
def test_a_pushed_face_runs_downstream():
    """A flat patch of 1 mm markers as a patch of the 100 mm sphere's shoulder 45 degrees from the stagnation point at
    93 s (the skyscrapers' PUSHES["flight"]: dp/ds -3.424e5 Pa/m and 80.29 Pa of shear down the face's meridian, and
    29.77 m/s2 of deceleration along the tilted flight direction), the meridian turned 30 degrees off the markers' grid
    so that the faces do not line up with the push, a disc of liquid 0.2 mm deep at 950 K over rigid metal
    at 880 K in its middle, rigid metal at 880 K around it, five 10 ms heat sub-steps of flow, the flow sub-steps
    chosen: the books close (mass and enthalpy to 1e-13, momentum to 1e-12 of the impulses), the markers downstream of
    the disc's centre (along the meridian) gain more than ten times what those upstream do, and the patch's volume flux points within
    10 degrees of the push."""
    nx = ny = 41                     # the runoff (some 10 m/s) stays 15 mm short of the patch's edge over the 50 ms;
    mk, ix, iy = grid_markers(nx, ny)   # on a 21 x 21 patch it reaches the edge and piles along it, turning the flux
    xc = (ix - 20) * 1e-3               # 16.5 degrees off the push (measured 2026-11-12)
    yc = (iy - 20) * 1e-3
    disc = np.hypot(xc, yc) < 5.5e-3
    sk = stacked(mk, [[LIQUID] * 20 + [RIGID] * 10 if d else [RIGID] * 30 for d in disc])
    a = np.radians(FLIGHT_ALPHA)
    turn = np.radians(30.0)
    meridian = np.array([-np.sin(turn), -np.cos(turn), 0.0])    # away from the stagnation point
    v_hat = -np.sin(a) * meridian + np.array([0.0, 0.0, np.cos(a)])
    f = flight_forcing(mk, v_hat, FLIGHT_DP_DS, FLIGHT_TAU, FLIGHT_DECEL, sigma=sk.material.liquid.sigma)
    G, _ = runoff.pushes(mk, sk, f)
    push = G[disc].mean(axis=0)
    assert push @ meridian > 0.0
    M0, E0, P0, m0 = sk.m.sum(), sk.H.sum(), sk.totals()[2].sum(axis=0), sk.m.sum(axis=1)
    flux = np.zeros(3)
    total = None
    for _ in range(5):
        book = runoff.flow_substeps(sk, mk, f, 0.01)
        flux += (sk.dz()[..., None] * sk.U).sum(axis=(0, 1))
        total = book if total is None else (total.add(book) or total)
    gain = np.maximum(sk.m.sum(axis=1) - m0, 0.0)
    along = xc * meridian[0] + yc * meridian[1]
    down, up = gain[along > 0.0].sum(), gain[along < 0.0].sum()
    angle = np.degrees(np.arccos(flux @ push / (np.linalg.norm(flux) * np.linalg.norm(push))))
    print("pushed face: gained downstream", down, "upstream", up, "kg; flux", angle, "deg from the push; cycles",
          total.cycles, "courant", total.courant)
    assert sk.m.sum() == pytest.approx(M0, rel=1e-13) and sk.H.sum() == pytest.approx(E0, rel=1e-13)
    assert momentum_closes(sk, P0, total, 1e-12)
    assert down > 10.0 * up and down > 0.0
    assert angle < 10.0


# -- the flow sub-steps' count -----------------------------------------------------------------------------------------
def test_the_adaptive_count_meets_its_targets_within_its_floor_and_cap():
    """The skyscrapers' count (`flow_cycles`): from no memory, 10; the next flow sub-step meets the Courant target (6) in
    proportion to the sub-step and the growth target (0.02) in proportion to its square, grows by at most twofold, and
    the count stays within 1 and 40."""
    dt, C, g = 0.01, runoff.FLOW_COURANT, runoff.FLOW_GROWTH
    assert runoff.flow_cycles(np.zeros(3), dt) == runoff.FLOW_CYCLES == 10
    assert runoff.flow_cycles([2.0 * C, 0.0, 1e-3], dt) == 20
    assert runoff.flow_cycles([0.5 * C, 0.0, 1e-3], dt) == 5
    assert runoff.flow_cycles([0.1 * C, 0.0, 1e-3], dt) == 5
    assert runoff.flow_cycles([0.0, 4.0 * g, 1e-3], dt) == 20
    assert runoff.flow_cycles([0.0, 0.0, 1e-2], dt) == runoff.FLOW_CYCLES_MIN == 1
    assert runoff.flow_cycles([1e4 * C, 0.0, 1e-3], dt) == runoff.FLOW_CYCLES_MAX == 40


# -- carried from Task 6 to Task 9: the release in every flow sub-step, the transport's iteration count ------------------
def test_the_flow_sub_steps_release_after_the_landing_and_book_it():
    """A release callable given to `flow_substeps` is called once in every flow sub-step, after the landing conduction
    (carried item (b) of Task 9), with the sub-step's length; what it takes leaves the skins and is booked in the
    FlowBook (mass, enthalpy, momentum), so that the skins' mass falls by exactly the released mass and their momentum
    changes by the impulses less the released momentum."""
    mk, ix, _ = grid_markers(12, 1)
    sk = stacked(mk, [[LIQUID] * 20 + [RIGID] * 5] * mk.n)
    f = flight_forcing(mk, STRIP_V_HAT, dp_ds=-2e4, tau=20.0)
    f.cycles = 4
    calls, order = [], []
    real_landing = runoff.landing_conduct

    def spy_landing(*a, **k):
        order.append("landing")
        return real_landing(*a, **k)

    def release(skins, h):
        calls.append(h)
        order.append("release")
        take = np.where(skins.valid[:, 0], 0.1 * skins.m[:, 0], 0.0)
        mass, H, solid, P, _ = skins.take_top_mass(take, skins.valid)
        return SimpleNamespace(mass=mass, enthalpy=H, momentum=P, solid=solid)

    M0, E0, P0 = sk.m.sum(), sk.H.sum(), sk.totals()[2].sum(axis=0)
    runoff.landing_conduct = spy_landing
    try:
        book = runoff.flow_substeps(sk, mk, f, 0.01, release=release)
    finally:
        runoff.landing_conduct = real_landing
    assert calls == [pytest.approx(0.0025)] * 4
    assert order[-1] == "release" and all(order[k + 1] == "release" for k, o in enumerate(order) if o == "landing")
    assert book.released.sum() > 0.0
    assert sk.m.sum() == pytest.approx(M0 - book.released.sum(), rel=1e-13)
    assert sk.H.sum() == pytest.approx(E0 - book.released_H.sum(), rel=1e-12)
    P1 = sk.totals()[2].sum(axis=0)
    impulse = (book.push + book.shear + book.wall + book.carried - book.released_p).sum(axis=0)
    size = sum(np.linalg.norm(v, axis=1).sum() for v in (book.push, book.shear, book.wall, book.released_p))
    assert np.linalg.norm(P1 - P0 - impulse) <= 1e-12 * size


def test_the_transport_books_its_corrections_and_does_not_fall_back_silently(monkeypatch):
    """Carried item (c) of Task 9: the transport's iterative solve books how many corrections it made and whether it
    converged; where it does not converge within TRANSPORT_CORRECTIONS it raises (no silent direct LU, which hung on
    133 000 cells), unless TRANSPORT_FALLBACK is set, when it finishes with the LU, warns and books converged False."""
    def strip():
        mk, ix, _ = grid_markers(20, 1)
        sk = stacked(mk, [[LIQUID] * 10 + [RIGID] * 5 if x < 5 else [RIGID] * 15 for x in ix])
        sk.U = np.where((sk.valid & (sk.T() > 900.0))[..., None], np.array([0.5, 0.0, 0.0]), 0.0)
        return mk, sk

    mk, sk = strip()
    book = runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
    assert book.converged and book.corrections == 0          # one-way flow: the upwind order solves it at once
    # a flow with a cycle (metal running both ways across one face) needs corrections, which the book counts
    mk, sk = strip()
    sk.U[1::2] *= -1.0
    book = runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
    assert book.converged and book.corrections > 0
    # with none allowed it raises
    mk, sk = strip()
    sk.U[1::2] *= -1.0
    monkeypatch.setattr(runoff, "TRANSPORT_CORRECTIONS", 0)
    with pytest.raises(RuntimeError, match="did not converge"):
        runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
    monkeypatch.setattr(runoff, "TRANSPORT_FALLBACK", True)
    mk, sk = strip()
    sk.U[1::2] *= -1.0
    with pytest.warns(RuntimeWarning, match="did not converge"):
        book = runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
    assert not book.converged
