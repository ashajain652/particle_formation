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


def test_a_new_top_cell_that_receives_nothing_is_not_created():
    """Stage 3 Task 11 (the q6p box's second step on the design's meshes): a new top cell whose arrival is exactly zero
    in floating point -- there a new top fed only through a cycle's back-coupling, which the solve left at zero within
    its tolerance, fed further new tops onward; here a flowing cell of 1e-300 kg moving at 1e-40 m/s onto a rigid top,
    whose arrival underflows -- is not installed as a cell of no mass (its velocity was 0/0 and the next landing
    conduction raised): the receiver keeps its cells exactly, nothing is booked as landed, and the books are unchanged."""
    mk, ix, _ = grid_markers(3, 1)
    sk = stacked(mk, [[RIGID] * 5, [RIGID] * 5, [RIGID] * 5])
    m, H, count = sk.m.copy(), sk.H.copy(), sk.count.copy()
    tiny = 1e-300                                        # kg, a liquid cell on top of the first marker
    m = np.concatenate([np.zeros((3, 1)), m], axis=1)
    H = np.concatenate([np.zeros((3, 1)), H], axis=1)
    m[0, 0], H[0, 0] = tiny, tiny * float(sk.props.enthalpy(np.array([LIQUID]))[0])
    m[1:, :-1], H[1:, :-1] = m[1:, 1:].copy(), H[1:, 1:].copy()
    m[1:, -1] = H[1:, -1] = 0.0
    sk.install(m, H, count + np.array([1, 0, 0]))
    sk.U[0, 0] = [1e-40, 0.0, 0.0]                       # toward the second marker
    before = [a.copy() for a in (sk.m, sk.H, sk.count)]
    book = runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
    assert np.array_equal(sk.count, before[2]) and not book.landed.any()
    assert np.array_equal(sk.m, before[0]) and np.array_equal(sk.H, before[1])
    assert np.all(np.isfinite(sk.U)) and np.all(sk.m[sk.valid] > 0.0)


def test_a_round_off_arrival_joins_the_top_cell_instead_of_landing():
    """Stage 3 Task 11's fix: an arrival on a rigid top lighter than LANDED_SLIVER of the top cell -- on the q6p box's cold
    step 14 972 of them, at most 4e-12 of their cell, from velocities at round-off, landed on 11 of its 13 face skins and
    turned on the landing conduction, which held the skins' fine zones adiabatic each flow sub-step and flattened their
    temperature profiles -- joins the receiver's top cell: no new cell, nothing booked as landed, the receiver's top cell
    holding the arrival's mass, enthalpy and momentum exactly. A cell of liquid moving at 1e-9 m/s onto rigid metal sends
    about 1e-9 of its mass in a 10 ms sub-step; at 1 m/s it lands as a new top cell, as before."""
    for speed, lands in ((1e-9, False), (1.0, True)):
        mk, ix, _ = grid_markers(3, 1)
        sk = stacked(mk, [[LIQUID] + [RIGID] * 4, [RIGID] * 5, [RIGID] * 5])
        sk.U[0, 0] = [speed, 0.0, 0.0]
        before = [a.copy() for a in (sk.m, sk.H, sk.count)]
        M0, E0, P0 = sk.m.sum(), sk.H.sum(), sk.totals()[2].sum(axis=0)
        book = runoff.transport(sk, runoff.Sheet.of(mk), 0.01)
        assert bool(book.landed[1]) == lands and (sk.count[1] > before[2][1]) == lands
        assert sk.m.sum() == pytest.approx(M0, rel=1e-14) and sk.H.sum() == pytest.approx(E0, rel=1e-14)
        P1 = sk.totals()[2].sum(axis=0)
        assert np.linalg.norm(P1 - P0 - (book.momentum + book.wall).sum(axis=0)) <= 1e-12 * np.linalg.norm(P0)
        if not lands:
            got = sk.m[1, 0] - before[0][1, 0]
            assert 0.0 < got <= runoff.LANDED_SLIVER * before[0][1, 0] and book.mass[1] == pytest.approx(got, rel=1e-12)
            assert np.array_equal(sk.m[1, 1:sk.count[1]], before[0][1, 1:before[2][1]])


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


def test_the_transport_reports_its_residual_against_its_tolerance(monkeypatch):
    """Task 9b (minor): a solve that does not converge reports its residual relative to its scale and the tolerance
    separately. The 2 x 2 system [[2, -1], [-1, 2]] x = (1, 0) in the natural order has one back-coupling; with no
    corrections allowed the solve stops at x = M^-1 b = (0.5, 0.25), whose residual (0.25, 0) is 0.25 / 2.25 = 0.1111 of
    the largest |b| + |A| |x| (2.25): 1.111e+12 times the tolerance of 1e-13. (Task 9's step 3 printed "6.5e-10 times its
    tolerance" for a residual of 6.5e-10 against 1e-13.)"""
    import scipy.sparse as sp
    monkeypatch.setattr(runoff, "TRANSPORT_CORRECTIONS", 0)
    A = sp.csr_matrix(np.array([[2.0, -1.0], [-1.0, 2.0]]))
    with pytest.raises(RuntimeError) as err:
        runoff._solve(A, np.array([1.0, 0.0]), np.arange(2))
    message = str(err.value)
    print(message)
    assert runoff.TRANSPORT_TOL == 1e-13
    assert "residual 1.111e-01 of its scale" in message
    assert "tolerance of 1e-13" in message and "1.111e+12 times" in message


# -- depth-resolved faces from the division (plan 03b Task 4, planning ruling 5) -----------------------------------------
import copy
import dataclasses

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.division import divide

from cases import DESIGN_SPACING, stage3_box, stage3_division

SHEET_OFF = 0.13e-3        # m: the flat face off the mesh's node planes (the stage 3 box's offset)
GRID_OFF = 0.07e-3         # m: the markers' grid off the node planes (test_levelset_division's offset) ...
JITTER = 0.05e-3           # m: ... and each marker off it by an irregular offset of up to this (the controller's ruling 2:
                           # a marker equidistant from two planes of nodes drops boundary length at some levels)
DEPTH_SETTINGS = SkinSettings(cell=20e-6)    # the development cells, H = 1.2 mm


def flat_division(n_side, seed=3):
    """A flat face (phi = y - SHEET_OFF, normal +y, the metal below it filling a box mesh of 0.5 mm cells) with n_side by
    n_side markers DESIGN_SPACING apart, off the mesh's planes of nodes by GRID_OFF and an irregular offset of up to
    JITTER, and the layer divided among them at the skin depth. Returns (positions (n, 3), normals (n, 3), Division)."""
    half = 0.5 * n_side * DESIGN_SPACING
    cells = int(round(2.0 * half / 0.5e-3))
    mesh = build(BoxBase((-half, -2.5e-3, -half), 0.5e-3, (cells, 6, cells)))
    phi = mesh.points[:, 1] - SHEET_OFF
    g = -half + (np.arange(n_side) + 0.5) * DESIGN_SPACING + GRID_OFF
    X, Z = np.meshgrid(g, g, indexing="ij")
    rng = np.random.default_rng(seed)
    x = np.column_stack([X.ravel() + rng.uniform(-JITTER, JITTER, X.size), np.full(X.size, SHEET_OFF),
                         Z.ravel() + rng.uniform(-JITTER, JITTER, X.size)])
    normal = np.tile([0.0, 1.0, 0.0], (len(x), 1))
    return x, normal, divide(mesh, phi, DEPTH_SETTINGS.depth, x, normal)


def division_markers(x, normal, d, mat, T, settings=DEPTH_SETTINGS):
    """Markers at x with the normals given, each with a skin whose measured profile is its column of the division d and
    whose top area is the column's top, filled full at T [K] (one or per marker). Returns the markers (skins set)."""
    n = len(x)
    sk = Skins(mat, d.area[:, 0], np.zeros((n, 2)), settings)
    sk.set_profile(d.levels, d.area)
    T = np.broadcast_to(np.asarray(T, dtype=float), (n,))
    sk.fill(sk.density(T) * sk.column_volume(), T)
    return Markers(x, normal, d.area[:, 0], np.zeros((n, 2)), np.arange(n), np.full(n, -1), DESIGN_SPACING, skins=sk)


def _same_cells(a, b, rel):
    assert np.array_equal(a.count, b.count)
    for p, q in ((a.m, b.m), (a.H, b.H), (a.U, b.U)):
        assert np.abs(p - q).max() <= rel * np.abs(q).max()


def test_depth_faces_with_constant_lengths_transport_as_before():
    """On a flat face the division's boundaries have the same length at every level and its columns the same
    cross-section, to 1e-12. With exactly constant lengths and cross-sections the sheet of the division
    (`Sheet.from_division`: each cell's rate (u.n) L(z)/A(z) at its mid-depth, ruling 5) and the sheet of the markers
    given the division's top faces (`Sheet.of`, stage 3's faces) carry the same metal: liquid over rigid metal (and four
    columns rigid at the top, which metal lands on and flows over), running one way at speeds varying with the column
    and the depth, one transport sub-step with a face push leaves the same cells and books to 1e-14, and a heat
    sub-step's flow sub-steps (the pushes, the momentum, the thickness terms, the transport, the landing conduction) to
    1e-12."""
    x, normal, d = flat_division(6)
    B = len(d.levels) - 1
    assert np.all(d.length[:, 0] > 0.0)
    assert np.abs(d.length / d.length[:, :1] - 1.0).max() <= 1e-12
    assert np.abs(d.area / d.area[:, :1] - 1.0).max() <= 1e-12
    const = dataclasses.replace(d, length=np.repeat(d.length[:, :1], B + 1, axis=1),
                                area=np.repeat(d.area[:, :1], B + 1, axis=1))
    mat = LinearMeltMaterial(rho_l=2800.0)            # no swelling: the temperatures below move no cell's depth
    mk = division_markers(x, normal, const, mat, RIGID)
    sk = mk.skins
    n = mk.n
    rng = np.random.default_rng(11)
    dz = sk.dz()
    z = np.cumsum(dz, axis=1) - 0.5 * dz
    liquid_to = rng.uniform(0.2e-3, 0.6e-3, n)
    liquid_to[[3, 14, 21, 32]] = 0.0                  # rigid tops
    T = np.where(z < liquid_to[:, None], LIQUID, RIGID)
    sk.H = np.where(sk.valid, sk.m * sk.props.enthalpy(T), 0.0)
    angle = np.deg2rad(25.0)
    e = np.array([np.cos(angle), 0.0, np.sin(angle)])
    speed = 0.2 * rng.uniform(0.5, 1.5, n)[:, None] * (1.0 + z / DEPTH_SETTINGS.depth)
    sk.U = np.where((sk.valid & sk.flowing())[..., None], speed[..., None] * e, 0.0)
    mk.faces = Faces(d.i, d.j, d.length[:, 0], d.normal)
    old, new = runoff.Sheet.of(mk), runoff.Sheet.from_division(const, mk)
    assert new.length_z.shape == (len(d.i), B + 1) and np.array_equal(new.levels, d.levels)
    np.testing.assert_allclose(new.f, old.f, rtol=1e-15, atol=0.0)
    np.testing.assert_allclose(new.w, old.w, rtol=1e-15)
    E = len(d.i)
    push = runoff.FacePush(np.where(sk.valid & sk.flowing(), 1e-7, 0.0), 0.2 * rng.uniform(-1.0, 1.0, 2 * E)
                           * np.r_[old.w, old.w] * 1e3)
    a, b = copy.deepcopy(sk), copy.deepcopy(sk)
    book_a = runoff.transport(a, old, 4e-3, faces=push)
    book_b = runoff.transport(b, new, 4e-3, faces=push)
    assert book_a.landed[[3, 14, 21, 32]].any() and book_a.courant > 0.1
    _same_cells(b, a, 1e-14)
    for p, q in ((book_b.mass, book_a.mass), (book_b.enthalpy, book_a.enthalpy), (book_b.momentum, book_a.momentum),
                 (book_b.wall, book_a.wall)):
        assert np.abs(p - q).max() <= 1e-14 * np.abs(q).max()
    v_hat = np.array([-0.6 * np.cos(angle), 0.8, -0.6 * np.sin(angle)])     # the meridian along e
    forcing = flight_forcing(mk, v_hat, dp_ds=-2e4, tau=20.0, decel=50.0)
    forcing.cycles = 3
    a, b = copy.deepcopy(sk), copy.deepcopy(sk)
    runoff.flow_substeps(a, old, forcing, 0.01)
    runoff.flow_substeps(b, new, forcing, 0.01)
    _same_cells(b, a, 1e-12)


def test_a_uniform_flow_over_a_flat_face_neither_gathers_nor_spreads():
    """Liquid filling every column of a flat face of 8 by 8 markers (the division's real lengths and cross-sections),
    every cell moving at 0.1 m/s in one direction 30 degrees off the markers' grid, one transport sub-step at a Courant
    number of about 1e-5: every interior column's mass changes by at most 1e-12 of itself. The interior is the columns
    at least two markers inside the face's edge: the columns at the edge, whose cells are open on the outer side, lose
    or gain metal at the Courant number, and the implicit sub-step carries that into their neighbours at its square
    and into the next at its cube (1e-15)."""
    n_side = 8
    x, normal, d = flat_division(n_side, seed=5)
    mk = division_markers(x, normal, d, LinearMeltMaterial(), LIQUID)
    sk = mk.skins
    angle = np.deg2rad(30.0)
    sk.U = np.where(sk.valid[..., None], 0.1 * np.array([np.cos(angle), 0.0, np.sin(angle)]), 0.0)
    sheet = runoff.Sheet.from_division(d, mk)
    m0 = sk.m.sum(axis=1)
    book = runoff.transport(sk, sheet, 2e-7)
    change = np.abs(sk.m.sum(axis=1) / m0 - 1.0)
    g = np.arange(n_side)
    gx, gz = np.meshgrid(g, g, indexing="ij")
    inner = ((gx >= 2) & (gx < n_side - 2) & (gz >= 2) & (gz < n_side - 2)).ravel()
    print("uniform flow: Courant {:.2e}; interior columns' largest change {:.2e}, the edge's {:.2e}".format(
        book.courant, change[inner].max(), change[~inner].max()))
    assert 1e-6 < book.courant < 1e-4 and change[~inner].max() > 1e-7
    assert change[inner].max() <= 1e-12
    assert sk.m.sum() == pytest.approx(m0.sum(), rel=1e-14)


def test_runoff_crosses_a_box_edge():
    """The stage 3 box's division (cases.stage3_division), liquid filling every column, and the pair of columns across
    one of its edges (normals more than 60 degrees apart) sharing the longest boundary at the skin depth: the cells of
    the column on one face run toward the edge at 0.2 m/s for one transport sub-step of 0.1 ms. Each of its cells sends
    at ruling 5's rate, the sum over its boundaries of (u.n)+ L(z)/A(z) at its mid-depth (the division's lengths and the
    column's cross-section, an independent oracle, 1e-12); the column across the edge receives the metal at the depths
    it left (each of its cells gains what the donor's cells sent from the depths it spans, 1e-9); the books close. The
    sheet hands the same boundaries to the sideways heat (`Sheet.side_faces`, in the skins' row order)."""
    box, d = stage3_box(), stage3_division()
    mk = division_markers(box.x, box.normal, d, LinearMeltMaterial(), LIQUID)
    sk = mk.skins
    full = d.area[:, -1] > 0.0
    edge = (np.einsum("ek,ek->e", box.normal[d.i], box.normal[d.j]) < 0.5) & full[d.i] & full[d.j]
    e = int(np.argmax(np.where(edge, d.length[:, -1], -1.0)))
    a, b = int(d.i[e]), int(d.j[e])
    t = d.normal[e] - (d.normal[e] @ box.normal[a]) * box.normal[a]
    u = 0.2 * t / np.linalg.norm(t)
    sk.U[a] = np.where(sk.valid[a][:, None], u, 0.0)
    sheet = runoff.Sheet.from_division(d, mk)
    side = sheet.side_faces()                                       # the same boundaries for the sideways heat
    assert np.array_equal(side.i, d.i) and np.array_equal(side.length, d.length)
    np.testing.assert_allclose(side.distance, np.linalg.norm(box.x[d.j] - box.x[d.i], axis=1), rtol=1e-15)
    M0, E0, P0 = sk.m.sum(), sk.H.sum(), sk.totals()[2].sum(axis=0)
    m_old, count_old = sk.m.copy(), sk.count.copy()
    dz = sk.dz()
    edges = [np.r_[0.0, np.cumsum(dz[c, :sk.count[c]])] for c in (a, b)]
    mid = 0.5 * (edges[0][1:] + edges[0][:-1])
    dt = 1e-4
    book = runoff.transport(sk, sheet, dt)
    # the donor's rate (ruling 5) from the division itself
    faces = np.flatnonzero((d.i == a) | (d.j == a))
    rate = np.zeros(len(mid))
    sends = {}
    A = np.interp(mid, d.levels, d.area[a])
    for f in faces:
        n_out = d.normal[f] if d.i[f] == a else -d.normal[f]
        c = max(u @ n_out, 0.0) * np.interp(mid, d.levels, d.length[f]) / A
        rate += c
        sends[f] = c
    k = count_old[a]
    kept = m_old[a, :k] / (1.0 + dt * rate)
    np.testing.assert_allclose(sk.m[a, :k], kept, rtol=1e-12)
    sent = dt * sends[e] * kept                                     # what each donor cell sent across the edge
    assert np.all(sent > 0.0) and sk.count[b] == count_old[b]       # into the receiver's flowing layer, no new cell
    # each receiving cell's gain: the donor's metal from the depths it spans (both columns full to H)
    za, zb = edges
    overlap = np.maximum(np.minimum(za[1:, None], zb[None, 1:]) - np.maximum(za[:-1, None], zb[None, :-1]), 0.0)
    expected = (sent[:, None] * overlap / np.diff(za)[:, None]).sum(axis=0)
    got = sk.m[b, :count_old[b]] - m_old[b, :count_old[b]]
    print("box edge: columns {} and {} (their normals' dot product {:.2f}), boundary {:.3f} to {:.3f} mm long from phi "
          "to H; the donor's rate across it from {:.1f} to {:.1f} 1/s".format(
              a, b, box.normal[a] @ box.normal[b], d.length[e, 0] * 1e3, d.length[e, -1] * 1e3, sends[e].min(),
              sends[e].max()))
    np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-9 * expected.max())
    assert book.mass[b] == pytest.approx(sent.sum(), rel=1e-12)
    assert sk.m.sum() == pytest.approx(M0, rel=1e-14) and sk.H.sum() == pytest.approx(E0, rel=1e-14)
    assert abs(book.mass.sum()) <= 1e-14 * M0
    P1 = sk.totals()[2].sum(axis=0)
    assert np.linalg.norm(P1 - P0 - (book.momentum + book.wall).sum(axis=0)) <= 1e-12 * np.linalg.norm(P0)


def test_landing_conduction_conducts_on_the_profile():
    """The controller's ruling 4: the landing conduction's blocks carry their columns' measured profile, so that metal
    landed on a profiled skin conducts on the profile's cross-sections, not on the flat top area. A film of 20 um at
    950 K on rigid metal at 700 K in a column narrowing by 60 % to H with a top area 10 % above the profile's top, one
    landing conduction of 0.5 ms: the block's cells end as a profiled column of those cells alone conducted with no heat
    across its top or bottom ends them, to 1e-12, and differently from the flat top area's conduction; the block keeps
    its energy and the column nothing landed on is untouched."""
    mat = LinearMeltMaterial()
    area = 4e-6 * (1.0 - 0.6 * np.linspace(0.0, 1.0, 25))
    levels = np.linspace(0.0, DEPTH_SETTINGS.depth, 25)
    settings = SkinSettings(cell=10e-6)
    sk = Skins(mat, np.full(2, 1.1 * area[0]), np.zeros((2, 2)), settings)
    sk.set_profile(levels, np.vstack([area, area]))
    sk.fill(mat.rho * sk.column_volume() * 0.8, 700.0)
    sk.H[:, :2] = sk.m[:, :2] * sk.props.enthalpy(LIQUID)
    m0, H0 = sk.m.copy(), sk.H.copy()
    dz = sk.dz()
    k = int(max(((np.cumsum(dz, axis=1) - dz)[0] < runoff.LANDING_DEPTH).sum(), runoff.LANDING_CELLS))
    landed = np.array([True, False])
    runoff.landing_conduct(sk, landed, 5e-4)

    def alone(profile):
        ref = Skins(mat, sk.area[:1], np.zeros((1, 2)), settings)
        if profile:
            ref.set_profile(levels, area[None, :])
        ref.install(m0[:1, :k], H0[:1, :k], [k])
        heat = ref.conduct(5e-4, 0.0, 0.0, 0.0, None)
        ref.H[0, k - 1] += heat.E_base[0]
        return ref.H[0, :k]

    on_profile, on_flat = alone(True), alone(False)
    change = np.abs(on_profile - H0[0, :k]).max()
    assert change > 0.0
    np.testing.assert_allclose(sk.H[0, :k], on_profile, rtol=0.0, atol=1e-12 * change)
    assert np.abs(on_flat - on_profile).max() > 1e-3 * change
    assert sk.H[0, :k].sum() == pytest.approx(H0[0, :k].sum(), rel=1e-14)
    assert np.array_equal(sk.H[0, k:], H0[0, k:]) and np.array_equal(sk.H[1], H0[1]) and np.array_equal(sk.m, m0)
