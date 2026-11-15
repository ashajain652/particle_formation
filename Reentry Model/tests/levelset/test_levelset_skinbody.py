"""Tests of the level set's skinned body (design sections 10.1, 12.2, 20 and 21; level-set stage 3 Task 9, plan 03b Task
5): one macro step coupling the skins on the gas-surface markers to the cut-cell 3D solve -- the books of mass and
enthalpy, a cold body left alone with its sharp edges, the gas surface phi moved by each column's volume and the skin
base beta derived one skin depth inside it, every column exactly full, the deferred load reaching the volume, and the
markers on phi after the surface moves.

The cube is the AA7075 Scheil alloy's (the shared Material). The test meshes are coarser than the design's defaults
(fine 1.4 mm instead of 0.7 mm, a 0.7 mm curvature floor instead of 0.25 mm, a 2 mm fine margin instead of 4 mm) so that a
macro step takes seconds; the skins are the development preset's (20 um cells)."""
import time

import numpy as np
import pytest

from reentry_model.levelset import skinbody
from reentry_model.levelset.boxmesh import BoxBase
from reentry_model.levelset.cutcell import inside_volumes
from reentry_model.levelset.division import DIVISION_TOL
from reentry_model.levelset.levelset import pieces
from reentry_model.levelset.skin import SKIN_PRESETS, Skins, SkinSettings
from reentry_model.levelset.surfaces import skin_base
from reentry_model.levelset.transfer import interpolate
from reentry_model.shared.material import Material

from cases import OWNERSHIP_FACTOR, gas_state, nearest_ratio

HALF = 5e-3            # m: the cube's half side (a 10 mm cube)
OFF = 0.13e-3          # m: the cube's centre off the grid, so that its faces cut the cells
SETTINGS = SkinSettings(cell=SKIN_PRESETS["dev"])
H = SETTINGS.depth     # m: the skin depth
MESH = dict(fine_margin=2e-3, fine_size=1.4e-3, cell_floor=0.7e-3)
OWED_SHARE = 0.05      # the heat still owed to the volume at a step's end, summed in absolute value per point, and the last
                       # part's interface mismatch, each at most this share of the step's heat in (fix round 1, Important 2)
OWED_SHARE_HEATED = 0.01    # the same on the heated face, whose step owes little (fix round 5)
PENDING_FLUX_SHARE = 10.0   # the pending heat's flux on a patch at most this times the largest surface flux (fix round 3)
CLOSE_RTOL = 1e-12          # each piece's phi against its metal's volume, relative (plan 03b Task 5)


def owed(body, book):
    """(the heat still owed to the volume, summed in absolute value [J], the last part's sum |mismatch| [J]) as shares of
    the step's heat in, printed with the step's deferred residuals by source; the groups within the guard's range and the
    pending flux bounded."""
    pend = np.abs(body.pending_J).sum() / book.enthalpy_in
    mism = book.interface_mismatch_J / book.enthalpy_in
    lo, hi = skinbody.GROUP_AREA_RANGE
    # a group formed again stops at the first patch that takes it past the base area, so it may overshoot the range's
    # top by that one patch (accepted, fix round 5): less its largest patch, every group is within the range
    assert lo <= book.group_area_ratio_min and book.group_area_ratio_trimmed_max <= hi
    print("  pending flux at most {:.3e} W/m2 against the surface's {:.3e} W/m2".format(book.pending_flux_max,
                                                                                    book.surface_flux_max))
    assert book.pending_flux_max <= PENDING_FLUX_SHARE * book.surface_flux_max
    print("  owed |pending| {:.2e} and mismatch {:.2e} of the heat in; deferred Newton {:.2e} J, moves {:.2e} J, "
          "remesh {:.2e} J; groups {:.2f} to {:.2f} of the base areas".format(
              pend, mism, book.residual_newton_J, book.residual_move_J, book.residual_remesh_J,
              book.group_area_ratio_min, book.group_area_ratio_max))
    return pend, mism


def closure(body):
    """The largest share by which a piece of phi's metal volume misses the skins' cells (m / rho(T)) and the volume's
    metal (beta's) in its cells and rows (plan 03b Task 5: every column exactly full). The carry the step end leaves
    (the volume's metal phi has yet to take) is not counted here: the steps' `density_carry_m3` and
    `density_carry_share_max` bound it (the controller's ruling T5-i)."""
    mesh, mk, sk = body.mesh, body.markers, body.skins
    label, n = pieces(mesh, body.phi)
    v_phi, v_beta = inside_volumes(mesh, body.phi), inside_volumes(mesh, body.beta)
    piece = np.where(body.division.piece >= 0, body.division.piece, 0)
    cells = sk.cells_volume()[mk.skin_index]
    worst = 0.0
    for p in range(n):
        V = v_phi[label == p].sum()
        metal = cells[piece == p].sum() + v_beta[label == p].sum()
        worst = max(worst, abs(V - metal) / V)
    return worst


def aa7075():
    return Material.from_drama_json("AA7075_scheil")


def cube(x):
    return np.max(np.abs(x - OFF) - HALF, axis=1)


ROUNDING = 2e-3        # m: the spraying cube's edge radius (the controller's ruling T6-d; more than the skin depth H)


def rounded_cube(x):
    """The 10 mm cube with its edges and corners rounded to ROUNDING (a signed distance), for the spraying fixture."""
    q = np.abs(x - OFF) - (HALF - ROUNDING)
    return np.linalg.norm(np.maximum(q, 0.0), axis=1) + np.minimum(q.max(axis=1), 0.0) - ROUNDING


def cube_body(T0=850.0, phi0=cube, **kw):
    base = BoxBase.around((-9e-3,) * 3, (9e-3,) * 3, 2.8e-3)
    args = dict(MESH, emissivity=0.0)
    args.update(kw)
    return skinbody.LevelSetSkinBody(aa7075(), base, phi0, T0, SETTINGS, **args)


PLATE_W, PLATE_DEPTH, PLATE_GAS, PLATE_H0 = 5.6e-3, 8e-3, 3.2e-3, 2.8e-3     # m: a column of a plate, walls the mesh's


def plate_body(T0=850.0, **kw):
    """One column of a plate (stage 3 Task 11's plates): metal below z = 0 to the mesh box's floor, its walls the mesh
    box's faces, one marker; the face is flat, so its column's cross-section is its area at every depth."""
    base = BoxBase.around((-PLATE_W / 2, -PLATE_W / 2, -PLATE_DEPTH), (PLATE_W / 2, PLATE_W / 2, PLATE_GAS), PLATE_H0)
    args = dict(spacing=PLATE_W)
    args.update(kw)
    return skinbody.LevelSetSkinBody(aa7075(), base, lambda x: x[:, 2], T0, SETTINGS, **args)


def surface_height(body, x=0.0, y=0.0):
    """The height [m] of phi's zero above the point (x, y) of the plate (linear between samples 1 um apart)."""
    z = np.linspace(-6e-3, 2e-3, 8001)
    v = interpolate(body.mesh, body.phi, np.column_stack([np.full(len(z), x), np.full(len(z), y), z]))
    k = np.flatnonzero((v[:-1] < 0.0) & (v[1:] >= 0.0))[-1]
    return float(z[k] + (z[k + 1] - z[k]) * (-v[k]) / (v[k + 1] - v[k]))


def one_face(q):
    """q [W/m2] on the markers facing +x (the cube's +x face, its edges by their normals' share), none elsewhere."""
    return lambda markers: q * np.clip((markers.normal[:, 0] - 0.5) / 0.5, 0.0, 1.0)


@pytest.fixture(scope="module")
def heated():
    """The 10 mm cube at 850 K heated 2 MW/m2 on its +x face for three 0.5 s macro steps."""
    body = cube_body()
    before = body.books()
    books, closures, seconds = [], [], []
    for _ in range(3):
        started = time.perf_counter()
        books.append(body.advance(0.5, q_prescribed=one_face(2e6)))
        seconds.append(time.perf_counter() - started)
        closures.append(closure(body))
    return body, before, books, closures, seconds


def test_a_macro_step_closes_mass_and_enthalpy(heated):
    """The 10 mm cube at 850 K heated 2 MW/m2 on its +x face for one 0.5 s macro step: the mass closes to 1e-12 and the
    enthalpy to 1e-8 of the stored enthalpy plus the heat in, the volume and the skins both gained heat, the deferred
    interface load is a small part of the heat in, and the heat in is what the face took."""
    body, before, books, _, seconds = heated
    book = books[0]
    print("heated cube: {:.1f} s for the step, {} parts, at most {} passes; heat in {:.4f} J; stored volume {:.4f} J + "
          "skins {:.4f} J; pending {:.3e} J; residuals mass {:.2e}, enthalpy {:.2e}; markers {}".format(
              seconds[0], book.parts, book.passes, book.enthalpy_in, book.enthalpy_volume, book.enthalpy_skins,
              book.pending_J, book.mass_residual, book.enthalpy_residual, book.n_markers))
    assert abs(book.mass_residual) <= 1e-12
    assert abs(book.enthalpy_residual) <= 1e-8
    area = (body.markers.area * np.clip((body.markers.normal[:, 0] - 0.5) / 0.5, 0.0, 1.0)).sum()
    assert book.enthalpy_in == pytest.approx(2e6 * area * 0.5, rel=0.05)
    assert book.enthalpy_volume > before["enthalpy_volume"] and book.enthalpy_skins > before["enthalpy_skins"]
    assert max(owed(body, book)) <= OWED_SHARE_HEATED
    for key in ("n_markers", "marker_splits", "marker_merges", "markers_trapped", "interface_mismatch_J",
                "volume_correction_max", "n_pieces", "frontal_area_m2", "layer_factor_max", "fill_miss_max",
                "move_miss_max", "n_baseless", "density_carry_m3"):
        assert key in book.row()
    for key in ("n_top", "modes_to_top", "mass_absorbed", "enthalpy_modes_J"):        # the retired modes' columns
        assert key not in book.row()
    assert book.n_pieces == 1 and book.frontal_area_m2 == pytest.approx((2 * HALF) ** 2, rel=0.02)
    for _, key in skinbody.SKIN_SUMMARY_COLUMNS:
        assert key in book.row()


CARRY_LEFT = 1.6e-11      # m3: the carry a heated step's end leaves for the next step's first move, |density_carry_m3|
                          # at most this (8.1e-12 m3, 8e-6 of the cube's volume, at most over the three steps)
CARRY_SHARE = 2.8e-3      # and no row's carry more than this share of its column's volume (1.42e-3 at most); each twice
                          # the value measured once, the controller's ruling T5-i


def test_every_step_closes_the_layer_and_the_books(heated):
    """Plan 03b Task 5 (design 21.8 item 5): over the heated cube's three steps, every division closes each piece to
    DIVISION_TOL (`layer_factor_max`), every fill leaves each column within 1e-12 of its volume (`fill_miss_max`), each
    piece's phi encloses its metal -- the skins' cells at their own densities and the volume's metal, beta's -- to 1e-12,
    the carry the step end leaves for the next step's first move is within CARRY_LEFT in all and CARRY_SHARE of any
    column (the controller's ruling T5-i), and the books close within the Global Constraints (mass 1e-12, enthalpy
    1e-8)."""
    body, _, books, closures, seconds = heated
    V = float(inside_volumes(body.mesh, body.phi).sum())
    for k, (book, c) in enumerate(zip(books, closures)):
        print("step {}: {:.1f} s, {} parts; factor {:.2e}, fill {:.1e}, closure {:.1e}, move miss {:.2e}, baseless {}, "
              "carry left {:.3e} m3 ({:.2e} of the body), at most {:.3e} of a column; drawn {:.3e} handed {:.3e} kg; "
              "residuals mass {:.1e} enthalpy {:.1e}".format(
                  k + 1, seconds[k], book.parts, book.layer_factor_max, book.fill_miss_max, c, book.move_miss_max,
                  book.n_baseless, book.density_carry_m3, book.density_carry_m3 / V, book.density_carry_share_max,
                  book.mass_drawn, book.mass_handed, book.mass_residual, book.enthalpy_residual))
        assert book.layer_factor_max <= DIVISION_TOL
        assert book.fill_miss_max <= skinbody.FILL_RTOL
        assert c <= CLOSE_RTOL
        assert abs(book.density_carry_m3) <= CARRY_LEFT and book.density_carry_share_max <= CARRY_SHARE
        assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert books[-1].mass_handed > 0.0                   # the melting face swelled and returned metal at its bases


def test_a_cold_box_keeps_its_edges_and_columns():
    """Design 21.8 item 1: the sharp 10 mm cube at rest, markers 2 mm apart, one 0.5 s step with no heat, no flow and no
    radiation: phi within the band is unchanged to 1e-12 m at every node (no geometric rule moves it, and the edges
    stay as sharp as the grid holds them), every column's volume is unchanged to 1e-12 relative, nothing crosses the
    base (`mass_drawn` and `mass_handed` zero), beta is skin_base(phi), the volume's temperatures and the skins' books
    are unchanged, and the mesh is kept."""
    body = cube_body()
    sk, mk = body.skins, body.markers
    phi0, T0, V0, n0 = body.phi.copy(), body.solver.T.copy(), sk.column_volume().copy(), mk.n
    M0, E0 = sk.m.sum(), sk.H.sum()
    nodes = len(body.mesh.points)
    book = body.advance(0.5)
    band = np.abs(phi0) < body.band
    print("cold box: phi moved {:.1e} m at most within the band; columns {:.1e}; drawn {:.1e} handed {:.1e} kg; "
          "markers {} then {}; residual mass {:.1e}".format(
              np.abs(body.phi - phi0)[band].max() if len(body.phi) == len(phi0) else np.inf,
              np.abs(sk.column_volume() / V0 - 1.0).max() if sk.n == len(V0) else np.inf, book.mass_drawn,
              book.mass_handed, n0, mk.n, book.mass_residual))
    assert len(body.mesh.points) == nodes and mk.n == n0
    assert np.abs(body.phi - phi0)[band].max() <= 1e-12
    np.testing.assert_allclose(sk.column_volume(), V0, rtol=1e-12, atol=0.0)
    assert book.mass_drawn == 0.0 and book.mass_handed == 0.0
    assert np.array_equal(body.beta, skin_base(body.mesh, body.phi, H, band=body.band))
    active = body.solver.active
    assert np.abs(body.solver.T[active] - T0[active]).max() <= 1e-9
    assert sk.m.sum() == pytest.approx(M0, rel=1e-14) and sk.H.sum() == pytest.approx(E0, rel=1e-14)
    assert book.enthalpy_in == 0.0 and abs(book.mass_residual) <= 1e-14


def test_markers_follow_the_recession():
    """The cube's +x face at 6 MW/m2 for 0.25 s with the liquid removed as it forms (remove="feed"): the face recedes
    (metal leaves, phi moves in by the columns' losses, beta follows one skin depth beneath, and the skins draw from the
    volume at their bases; 0.44 mm at its centre, measured once), the books close, and after the step every marker lies
    within 1e-9 m of phi = 0. (For 0.5 s, the step of stage 3's test, the face's centre recedes 1.7 mm while its rim,
    heated by its normals' share, stands; the rim is then thinner than twice the skin depth, and the step end stops on
    `stepend.THIN_METAL`, as design 21.11 asks: measured once.) The step's
    book lists every volume correction's miss before it (`move_corrections`), each entry labelled with what made it --
    a part's move of phi ("part"), the step end's move by the last part's carry ("carry"), the remesh ("remesh") -- so
    that the largest of the list is `volume_correction_max`; each entry gives beside the worst piece's miss the body's,
    the largest piece's."""
    body = cube_body(remove="feed")
    x0 = body.markers.x[body.markers.normal[:, 0] > 0.9, 0].mean()
    book = body.advance(0.25, q_prescribed=one_face(6e6))
    mk = body.markers
    face = mk.normal[:, 0] > 0.9
    off = np.abs(interpolate(body.mesh, body.phi, mk.x))
    worst = [e["worst"] for e in book.move_corrections]
    print("receding face: removed {:.3e} kg, drawn {:.3e} kg; face moved {:.3f} mm; markers off phi at most {:.2e} m; "
          "residuals mass {:.2e}, enthalpy {:.2e}; {} volume corrections ({}) missing by {:.2e} at most (median {:.2e}) "
          "of their piece before the correction; {} parts".format(
              book.mass_removed, book.mass_drawn, (x0 - mk.x[face, 0].mean()) * 1e3, off.max(), book.mass_residual,
              book.enthalpy_residual, len(worst), ", ".join(sorted({e["kind"] for e in book.move_corrections})),
              max(worst, default=np.nan), np.median(worst) if worst else np.nan, book.parts))
    kinds = [e["kind"] for e in book.move_corrections]
    assert set(kinds) <= {"part", "carry", "remesh"} and kinds.count("part") >= 1 and min(worst) >= 0.0
    assert max(worst) == book.volume_correction_max
    body_miss = [e["body"] for e in book.move_corrections]
    assert all(0.0 <= m <= w for m, w in zip(body_miss, worst))
    assert book.mass_removed > 0.0 and book.mass_drawn > 0.0
    assert mk.x[face, 0].mean() < x0
    assert off.max() <= 1e-9
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert max(owed(body, book)) <= OWED_SHARE


def test_a_sprayed_column_draws_from_the_volume_at_its_base():
    """Design 21.5 and 21.6: metal sprayed from a face column's top moves phi in by its volume over the column, beta
    follows one skin depth beneath, and the column, still as deep, draws from the volume at its base exactly what it
    lost. The plate's one flat column (its walls the mesh box's), 0.1 mm of its top taken (as a release takes it): the
    volume it draws at its base, the new cell's m / rho(T), equals the volume taken to 1e-12; beta falls beneath it by
    that volume over the column's area; the column is full; and the body's mass, with what was taken, closes to 1e-12."""
    body = plate_body()
    sk = body.skins
    A, M0 = float(sk.column_volume()[0] / H), body.books()
    beta_z = lambda: zero_height(body, body.beta)
    b0, z0 = beta_z(), surface_height(body)
    m, Hc, _ = (a[0] for a in sk.take_between([0], 0.0, 0.1e-3))       # the release takes the top 0.1 mm
    taken = float(sk.column_volume()[0] - sk.cells_volume()[0])
    tally = skinbody._Tally(sk.n)
    v0 = sk.cells_volume()[0]
    body._layer(sk.cells_volume() - sk.column_volume(), tally, "part")
    drawn = float(sk.cells_volume()[0] - v0)                            # what the base gave, at the cells' density
    M1 = body.books()
    mass = (M1["mass_volume"] + M1["mass_skins"] + m - M0["mass_volume"] - M0["mass_skins"]) / \
        (M0["mass_volume"] + M0["mass_skins"])
    print("sprayed column: took {:.6e} m3, drew {:.6e} m3 ({:.1e} apart); phi fell {:.4f} mm, beta {:.4f} mm, asked "
          "{:.4f} mm; mass {:.1e}".format(taken, drawn, abs(drawn / taken - 1), (z0 - surface_height(body)) * 1e3,
                                         (b0 - beta_z()) * 1e3, taken / A * 1e3, mass))
    assert drawn == pytest.approx(taken, rel=1e-12)
    assert b0 - beta_z() == pytest.approx(taken / A, rel=1e-6)
    assert abs(sk.cells_volume()[0] / sk.column_volume()[0] - 1.0) <= skinbody.FILL_RTOL
    assert abs(mass) <= 1e-12 and tally.drawn > 0.0 and tally.handed == 0.0


DRAW_FIELD = (800.0, 5e3)   # K, K/m: the volume's temperature T0 + g (x + y + z) in the draw test
DRAW_APART_K = 1.0          # K: the field read at the edge column's base point differs from its group's mean by more


def test_an_edge_column_draws_at_its_group_s_base_temperature(monkeypatch):
    """Plan 03b Task 5's fix round 1 (the controller's ruling T5-k): a column with a base draws the volume's metal at
    its own base temperature, the area mean over its group of beta's patch temperatures that its base law reads, not at
    the base point x - H n, which lies outside beta for every edge and corner marker. On the test cube with the volume's
    temperature a linear field T0 + g (x + y + z) (DRAW_FIELD), the edge or corner column whose base point lies farthest
    outside beta, 0.1 mm of its top taken: the fill draws it at exactly its group's mean of the patch temperatures,
    which the field at its base point misses by more than DRAW_APART_K, and the column is full."""
    body = cube_body()
    sk, mk, solver = body.skins, body.markers, body.solver
    T0, g = DRAW_FIELD
    solver.T = T0 + g * body.mesh.points.sum(axis=1)
    rows = mk.skin_index[np.abs(mk.normal).max(axis=1) < 0.9]                        # the edge and corner columns
    rows = rows[~body._baseless()[rows]]
    outside = interpolate(body.mesh, body.beta, body._bases()[rows])
    row = int(rows[np.argmax(outside)])
    sk.take_between([row], 0.0, 0.1e-3)
    S, A_group, _ = body._interface()
    assert A_group[row] > 0.0
    want = float((S @ solver.patch_T())[row])
    at_point = float(interpolate(body.mesh, solver.T, body._bases()[[row]])[0])
    calls = []
    exchange = sk.exchange_base

    def spy(r, dV, T_draw=None):
        calls.append((np.array(r), None if T_draw is None else np.broadcast_to(T_draw, np.shape(r)).copy()))
        return exchange(r, dV, T_draw)
    monkeypatch.setattr(sk, "exchange_base", spy)
    body._fill(skinbody._Tally(sk.n))
    got = [T[list(r).index(row)] for r, T in calls if row in r]
    print("edge column {} (base point {:.3f} mm outside beta): drew at {} K; its group's mean {:.4f} K; the field at "
          "its base point {:.4f} K".format(row, outside.max() * 1e3, np.round(got, 4), want, at_point))
    assert len(got) >= 1 and all(t == pytest.approx(want, rel=1e-12) for t in got)
    assert abs(at_point - want) > DRAW_APART_K
    assert abs(sk.cells_volume()[row] / sk.column_volume()[row] - 1.0) <= skinbody.FILL_RTOL


def zero_height(body, f, x=0.0, y=0.0):
    """The height [m] of the level set f's zero above the point (x, y) of the plate (linear between 1 um samples)."""
    z = np.linspace(-4e-3, 2e-3, 6001)
    v = interpolate(body.mesh, f, np.column_stack([np.full(len(z), x), np.full(len(z), y), z]))
    k = np.flatnonzero((v[:-1] < 0.0) & (v[1:] >= 0.0))[-1]
    return float(z[k] + (z[k + 1] - z[k]) * (-v[k]) / (v[k + 1] - v[k]))


def _edge_column(body):
    """An edge column of the cube: the marker whose normal leans most from every face's (the least largest
    component)."""
    return int(np.argmin(np.abs(body.markers.normal).max(axis=1)))


def _baseless_edge(body, monkeypatch):
    """An edge column of the test cube made baseless (`_edge_column`; its division's cross-section at H given to the
    neighbour it shares its deepest level with, so that the columns still tile the layer; `stepend.divide` wrapped for
    the rest of the test) and the layer divided afresh: (its marker, its row, that neighbour's marker)."""
    from dataclasses import replace
    from reentry_model.levelset import stepend
    m = _edge_column(body)
    d = body.division
    e = np.flatnonzero((d.i == m) | (d.j == m))
    deep = e[np.argmax([np.max(np.flatnonzero(d.length[k] > 0.0), initial=-1) for k in e])]
    other = d.j[deep] if d.i[deep] == m else d.i[deep]
    w = np.full(len(d.levels), d.levels[1]) * np.r_[0.5, np.ones(len(d.levels) - 2), 0.5]     # the trapezoid's
    divide = stepend.divide

    def baseless(*args, **kwargs):                      # the edge column's base given to that neighbour's
        out = divide(*args, **kwargs)
        area = out.area.copy()
        area[other, -1] += area[m, -1]
        area[m, -1] = 0.0
        return replace(out, area=area, volume=area @ w)
    monkeypatch.setattr(stepend, "divide", baseless)
    body._divide(skinbody._Tally(body.skins.n))
    return m, body.markers.skin_index[m], other


def test_a_column_with_no_base_takes_its_shortfall_from_its_neighbours(monkeypatch):
    """Design 21.5: a column that ends before beta, beside a convex edge, has no base: it never draws or returns, and is
    filled from the neighbour across the boundary it shares deepest. On the test cube an edge column made baseless (its
    division's cross-section at H given to that neighbour's, so that the columns still tile the layer), 0.05 mm of its
    top taken (as a release takes it) and the layer accepted (phi moved by it, the division, the fill): the fill never exchanges it at a base, it is full afterwards, the neighbour it
    took from (`_neighbours`: one with a base, sharing its deepest level) is filled at its own base, and the books -- the skins' and the volume's mass with what was taken, the enthalpy with the
    heat owed -- are exact."""
    body = cube_body()
    sk, mk = body.skins, body.markers
    m, row, other = _baseless_edge(body, monkeypatch)
    assert body._baseless()[row] and not body._baseless()[mk.skin_index[other]]
    body._fill(skinbody._Tally(sk.n))                    # the two columns' volumes moved, their sum kept: settle them
    before = body.books()
    giver, _ = body._neighbours(row, ~body._baseless())
    m_out, H_out, _ = (a[0] for a in sk.take_between([row], 0.0, 0.05e-3))
    exchanged = []
    original = Skins.exchange_base

    def watched(self, rows, dV, T_draw=None):
        exchanged.extend(np.asarray(rows).tolist())
        return original(self, rows, dV, T_draw)
    monkeypatch.setattr(Skins, "exchange_base", watched)
    tally = skinbody._Tally(sk.n)
    body._layer(sk.cells_volume() - sk.column_volume(), tally, "part")
    after = body.books()
    M0, M1 = before["mass_volume"] + before["mass_skins"], after["mass_volume"] + after["mass_skins"] + m_out
    E0 = before["enthalpy_volume"] + before["enthalpy_skins"] + before["pending_J"]
    E1 = after["enthalpy_volume"] + after["enthalpy_skins"] + after["pending_J"] + H_out
    miss = abs(sk.cells_volume()[row] / sk.column_volume()[row] - 1.0)
    print("baseless edge column: took {:.3e} kg; exchanged at the bases of {} columns ({} its neighbour), its miss after "
          "{:.1e}; books mass {:.1e} enthalpy {:.1e}".format(m_out, len(set(exchanged)), giver in exchanged, miss,
                                                             M1 / M0 - 1, E1 / E0 - 1))
    assert row not in exchanged and giver in exchanged
    assert miss <= skinbody.FILL_RTOL and tally.fill_miss <= skinbody.FILL_RTOL
    assert abs(M1 / M0 - 1.0) <= 1e-12 and abs(E1 / E0 - 1.0) <= 1e-12


# -- clear stops and the speck rule (plan 03b Task 7b; the controller's rulings T7-f and T7-g) ------------------------
SPECK_DEPTH = 6.5e-5     # m: the parted speck's node below phi = 0 in the speck test (the spraying cube's, diagnosis D2)
SPECK_SHORT = 0.2e-6     # m: how far the move leaves the body shy of its metal there (D2: the correction then lowered
                         # phi 0.24 um on the spraying cube and 3.1 um on the shoulder, far more than PIECE_GONE's 1 nm)


def test_a_speck_the_move_removes_stays_removed_after_the_volume_correction(monkeypatch):
    """The controller's ruling T7-f (diagnosis D2): a speck the move of phi parts from the body (a piece with no column
    and none of the volume's metal) is lifted out by `_move_phi`'s speck rule, its nodes to `stepend.PIECE_GONE` (1 nm
    outside the metal), and must stay out through the per-piece volume correction. Before, the correction's uniform shift
    over the body's claimed nodes (two rings of cells about it) took a lifted node back below zero whenever the body grew
    by more than a nanometre, and the speck came back with no marker: it captured two drifted markers' columns on the
    spraying cube (the split that never ended) and was swallowed by no group's cells on the shoulder (MOVE_SUM_RTOL). The
    plate's column, its move of phi made (`stepend.move_by_volumes` wrapped) to leave the body SPECK_SHORT shy of its
    metal and to part a speck at a node two rings of cells above its face: after `_move_phi` phi is one piece holding its
    metal to 1e-12, the speck counted as used up, its nodes left out of the correction's claims (`fixed`); and the
    correction given those nodes to claim raises, naming the piece it made in cells that held no metal before it, rather
    than make the speck again."""
    from reentry_model.levelset import levelset as ls
    from reentry_model.levelset import stepend
    body = plate_body()
    mesh, sk = body.mesh, body.skins
    label, _ = pieces(mesh, body.phi)
    rings, nodes = [], np.unique(mesh.tets[label == 0])
    for _ in range(2):                                   # the body's claim in the correction: two rings of cells
        nodes = np.unique(mesh.tets[np.isin(mesh.tets, nodes).any(axis=1)])
        rings.append(nodes)
    two = np.setdiff1d(rings[1], rings[0])
    up = two[mesh.points[two, 2] > 0.0]
    node = int(up[np.argmin(np.linalg.norm(mesh.points[up, :2], axis=1))])
    V0 = float(inside_volumes(mesh, body.phi).sum())
    real_move, real_correct = stepend.move_by_volumes, stepend.correct_volumes
    calls = []

    def parted(mesh_, phi, x, normal, dV, band=None):
        moved = real_move(mesh_, phi, x, normal, dV, band=band)
        out = moved.phi + SPECK_SHORT
        out[node] = -SPECK_DEPTH
        gained = float(inside_volumes(mesh_, out).sum() - inside_volumes(mesh_, moved.phi).sum())
        return moved._replace(phi=out, swept=moved.swept + gained / len(moved.swept))

    def correct(mesh_, phi, label_, targets, **kw):
        calls.append((phi.copy(), label_.copy(), np.array(targets, float), kw))
        return real_correct(mesh_, phi, label_, targets, **kw)
    monkeypatch.setattr(stepend, "move_by_volumes", parted)
    monkeypatch.setattr(stepend, "correct_volumes", correct)
    tally = skinbody._Tally(sk.n)
    body._move_phi(np.zeros(sk.n), tally, "part")
    label, n = pieces(mesh, body.phi)
    V1 = float(inside_volumes(mesh, body.phi).sum())
    phi_in, label_in, targets, kw = calls[0]
    print("speck at node {} ({:.3f} mm above the face): pieces after the move {}; phi at the node {:.3e} m; the body's "
          "volume {:.1e} from its metal".format(node, mesh.points[node, 2] * 1e3, n, body.phi[node], V1 / V0 - 1.0))
    assert n == 1 and body.phi[node] > 0.0
    assert tally.pieces_used == 1
    assert V1 == pytest.approx(V0, rel=1e-12)
    assert kw["fixed"][node] and phi_in[node] == stepend.PIECE_GONE
    with pytest.raises(RuntimeError, match="held no metal before"):
        ls.correct_volumes(mesh, phi_in, label_in, targets)


def test_a_column_the_fill_cannot_fill_stops_the_run_naming_it(monkeypatch):
    """The controller's ruling T7-g (diagnosis D2): a column the fill leaves beyond FILL_RTOL of its volume after its
    rounds stops the run with a RuntimeError naming its marker, its position, its piece, its cells' volume and its
    column's volume. Before, the fill skipped a column with no base and no neighbour silently (the spraying cube's
    second step: two columns stranded on a speck kept 0.73 and 0.52 mm3 of metal in columns of 1e-20 m3, the next part's
    move ratio was 4.3e9, and the split never ended). On the test cube an edge column made baseless (`_baseless_edge`)
    and given no neighbour, 0.05 mm of its top taken: the fill raises, naming it."""
    from reentry_model.levelset import stepend
    body = cube_body()
    sk, mk = body.skins, body.markers
    m, row, _ = _baseless_edge(body, monkeypatch)
    body._fill(skinbody._Tally(sk.n))
    sk.take_between([row], 0.0, 0.05e-3)
    monkeypatch.setattr(stepend.StepEnd, "_neighbours", lambda self, r, base: None)
    with pytest.raises(RuntimeError, match="the fill left") as stop:
        body._fill(skinbody._Tally(sk.n))
    print("unfilled column: {}".format(stop.value))
    x = "({:.3f}, {:.3f}, {:.3f}) mm".format(*(mk.x[m] * 1e3))
    held, column = "{:.6e} m3".format(sk.cells_volume()[row]), "{:.6e} m3".format(sk.column_volume()[row])
    for part in ("marker {} ".format(m), x, "piece {}".format(body.division.piece[m]), held, column):
        assert part in str(stop.value), part


def test_a_piece_with_layer_and_no_column_stops_the_run(monkeypatch):
    """The controller's ruling T7-g (diagnosis D2): a division that finds a piece of phi with layer but no column (no
    marker on it, its closure factor infinite, as two of the shoulder's divisions did unseen before its stop) stops the
    run with a RuntimeError naming the piece, and so does a piece whose closure factor departs from one by more than
    `stepend.DIVISION_GROSS` (the speck that took two columns on the spraying cube closed with a factor of 0.0033); the
    body only recorded the factor before. Two spheres on 1 mm cells, 2.2 mm and 1.8 mm in radius (the larger piece 0)
    and 2.8 mm apart, markers on the larger alone (26 directions): `_divide` raises naming the smaller, piece 1; with the
    division's factor for the larger replaced by 1 + 2 DIVISION_GROSS, it raises naming that."""
    from dataclasses import replace
    from types import SimpleNamespace
    from reentry_model.levelset import stepend
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.levelset import redistance
    from reentry_model.levelset.markers import Markers

    class Host(stepend.StepEnd):
        pass
    mesh = build(BoxBase.around((-7e-3, -4e-3, -4e-3), (7e-3, 4e-3, 4e-3), 1e-3))
    a, b, R = np.array([-3.4e-3, OFF, OFF]), np.array([3.4e-3, OFF, OFF]), 2.2e-3
    phi = redistance(mesh, np.minimum(np.linalg.norm(mesh.points - a, axis=1) - R,
                                      np.linalg.norm(mesh.points - b, axis=1) - 1.8e-3))
    assert pieces(mesh, phi)[1] == 2
    u = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)], float)
    u /= np.linalg.norm(u, axis=1)[:, None]
    n = len(u)
    host = Host()
    host.mesh, host.phi, host.settings = mesh, phi, SETTINGS
    host.markers = Markers(a + R * u, u, np.ones(n), np.zeros((n, 2)), np.arange(n), np.full(n, -1), 2e-3)
    host.skins = SimpleNamespace(n=n, set_profile=lambda levels, area: None)
    with pytest.raises(RuntimeError, match="layer but no column") as stop:
        host._divide(skinbody._Tally(n))
    print("a piece with no column: {}".format(stop.value))
    assert "piece 1" in str(stop.value)
    divide = stepend.divide
    monkeypatch.setattr(stepend, "divide", lambda *a_, **k: (lambda d: replace(
        d, factor=np.r_[1.0 + 2.0 * stepend.DIVISION_GROSS, d.factor[1:]]))(divide(*a_, **k)))
    with pytest.raises(RuntimeError, match="closure factor") as stop:
        host._divide(skinbody._Tally(n))
    print("a gross factor: {}".format(stop.value))
    assert "piece 0" in str(stop.value)


def test_a_split_past_PARTS_MAX_stops_naming_the_worst_column(monkeypatch):
    """The controller's ruling T7-g (diagnosis D2): a split never asks for more than PARTS_MAX shorter parts; past that the
    step stops with a RuntimeError naming the column whose move sets the part's ratio (its marker and its ratio), where
    it once asked for ceil(4.3e9) parts at every level to SPLIT_DEPTH_MAX and never ended. The test cube at rest with a
    carry (planning ruling 4's volume change at a column's top) of three PARTS_MAX allowances on one face column, which
    no part's length shrinks: the step raises after its first part, naming that column's marker and the ratio 3 PARTS_MAX
    (a part run again would be the split, and stops the test instead of hanging it)."""
    import re
    body = cube_body()
    sk, mk = body.skins, body.markers
    m = int(np.flatnonzero(mk.normal[:, 0] > 0.9)[0])
    row = mk.skin_index[m]
    allowance = skinbody.MOVE_SHARE * body.spacing * sk.profile[row, 0]          # m3: the volume a part may move it
    body.carry[row] = 3.0 * skinbody.PARTS_MAX * allowance
    runs = []
    part_once = skinbody.LevelSetSkinBody._part_once

    def once(self, *a, **k):
        if runs:
            raise AssertionError("the part was split into {} parts".format(skinbody.PARTS_MAX * 3))
        runs.append(1)
        return part_once(self, *a, **k)
    monkeypatch.setattr(skinbody.LevelSetSkinBody, "_part_once", once)
    with pytest.raises(RuntimeError, match="PARTS_MAX") as stop:
        body.advance(0.05)
    print("split past PARTS_MAX: {}".format(stop.value))
    assert "marker {} ".format(m) in str(stop.value)
    ratio = float(re.search(r"ratio ([0-9.e+-]+)", str(stop.value)).group(1))
    assert ratio == pytest.approx(3.0 * skinbody.PARTS_MAX, rel=1e-6)


BEAM_X = np.array([1.0, 0.0, 0.0])   # a beam arriving along -x, on the cube's +x face
BEAM_Q = 6.0e6                       # W/m2


def _beam(body, q0=BEAM_Q):
    """A prescribed beam q0 along -BEAM_X as each marker takes it on phi as it is: the mean of its patches' projections
    max(n_p . BEAM_X, 0) over the area of phi it owns under the ownership that sets its area (`Markers.owned_patches`;
    none for a marker that owns none), times q0 [W/m2] per marker, a function of the markers (`advance`'s q_prescribed)."""
    from reentry_model.levelset.cutcell import patches

    def q(mk):
        P = patches(body.mesh, body.phi)
        W = mk.owned_patches(body.mesh, body.phi, P)
        A = np.asarray(W.sum(axis=1)).ravel()
        f = q0 * np.maximum(P.normal @ BEAM_X, 0.0)
        return np.divide(np.asarray(W @ f).ravel(), A, out=np.zeros(mk.n), where=A > 0.0)
    return q


def _watch_passes(monkeypatch):
    """Every pass of the body's parts recorded (a list, returned): (phi then, the marker fluxes it takes [W/m2] in the
    markers' order, the heat they take per second, sum q_i A_i [W], the beam's power on phi then [W], and the runoff's
    forcing it takes)."""
    from reentry_model.levelset.cutcell import patches
    seen = []
    real_pass = skinbody.LevelSetSkinBody._pass

    def watched(self, dt, T_b, ctx, trial, sheet):
        mk, sk = self.markers, self.skins
        P = patches(self.mesh, self.phi)
        q = ctx["q"][mk.skin_index].copy()
        seen.append((self.phi, q, float((q * sk.area[mk.skin_index]).sum()),
                     BEAM_Q * float((P.area * np.maximum(P.normal @ BEAM_X, 0.0)).sum()), ctx["forcing"]))
        return real_pass(self, dt, T_b, ctx, trial, sheet)
    monkeypatch.setattr(skinbody.LevelSetSkinBody, "_pass", watched)
    return seen


def _surfaces(seen):
    """The surfaces (phi's arrays) the recorded passes ran on, each once, in order."""
    out = []
    for phi, *_ in seen:
        if not out or out[-1] is not phi:
            out.append(phi)
    return out


def test_each_markers_heat_is_retaken_on_the_moved_surface_every_part(monkeypatch):
    """The controller's ruling T7-i (diagnosis D1; design 21.5, one partition of phi for the areas and the fluxes, and
    21.6, the areas given afresh on the moved phi): each marker's top flux is taken under `Markers.owned_patches`, the
    ownership that sets its area, and taken again on the moved surface whenever the areas are given afresh, at every
    part. Before, the flux was fixed at the step's start while the areas it multiplies followed phi, so that on the
    stage 3 box the heat booked in a part ran 25 % over the beam's power once its face was uneven. The test cube with
    its liquid removed as it forms, one 0.25 s step of a 6 MW/m2 beam on its +x face (`_beam`, a function of the
    markers), in which the face recedes unevenly (its rim stands) over several parts: at every pass of every part the
    heat the skins take, sum q_i A_i, is the beam's power on phi as it is to 1e-9."""
    seen = _watch_passes(monkeypatch)
    body = cube_body(remove="feed")
    book = body.advance(0.25, q_prescribed=_beam(body))
    worst = max(abs(took / beam - 1.0) for _, _, took, beam, _ in seen)
    print("beam: {} parts, {} passes on {} surfaces; the heat taken against the beam's power on phi as it is, {:.1e} at "
          "worst".format(book.parts, len(seen), len(_surfaces(seen)), worst))
    assert len(_surfaces(seen)) >= 2
    assert worst <= 1e-9


def test_a_flight_takes_its_loads_once_per_part_under_the_part_s_ownership(monkeypatch):
    """The controller's ruling T7-i in flights, as ruling T7b-a settles it: the convective flux (`loads.by_patches`,
    from the step's one aero state) is taken again on the moved surface at every part, under the ownership that sets the
    markers' areas there (`Markers.owned_patches`), no new gas state evaluated, while the runoff's forcing is held at the
    step's start (taken again, it folded the spraying fixture's liquid under the cube's bottom face, measured once in
    plan 03b Task 7b). The gas body at the gas state raised to 1e-4 kg/m3, a 0.1 s step started as three parts:
    `by_patches` is called once for each surface a part runs on, in order, with that surface's ownership, every pass on
    it takes that call's fluxes, and every pass takes the step start's forcing."""
    from reentry_model.levelset import loads
    seen = _watch_passes(monkeypatch)
    calls = []
    real = loads.by_patches

    def by_patches(loads_, state, markers, P, *a, W=None, **k):
        apart = abs(W - markers.owned_patches(body.mesh, body.phi, P)).max() if W is not None else np.inf
        q, flow = real(loads_, state, markers, P, *a, W=W, **k)
        calls.append((body.phi, apart, q.copy()))
        return q, flow
    monkeypatch.setattr(loads, "by_patches", by_patches)
    body = gas_body()
    body._hint = 3                                       # the step starts as three parts
    book = body.advance(0.1, state=gas_state(rho=1e-4))
    ran = _surfaces(seen)
    print("flight: {} parts, {} by_patches calls, {} surfaces the parts ran on; ownership apart by {:.1e} m2 at "
          "most".format(book.parts, len(calls), len(ran), max(c[1] for c in calls)))
    assert len(ran) >= 2
    assert len(calls) == len(ran) and all(c[0] is phi for c, phi in zip(calls, ran))    # once per surface, in order
    assert all(c[1] == 0.0 for c in calls)               # under the ownership that sets the markers' areas there
    for phi, q, _, _, _ in seen:
        assert np.array_equal(q, next(c[2] for c in calls if c[0] is phi))
    assert all(f is seen[0][4] for *_, f in seen)          # the forcing held at the step's start (ruling T7b-a)


def test_a_molten_cell_returned_to_the_volume_lowers_phi_above_its_column():
    """Design 21.6 and planning ruling 4: the volume keeps one density, the solid's, so a molten cell returned at a
    column's base shrinks by 1 - rho_l / rho_s, and that difference is carried as a volume change at the top of that
    column into the next move of phi (the next part's, or the step end's). The plate's column given a molten cell of
    0.1 mm at 950 K at its bottom (overfull): its part's move takes phi up by the cell's volume, the fill returns the cell,
    and by the step end phi has fallen back by (1 - rho_l / rho_s) times the cell's volume over the column's area, within
    5 %; the books close to 1e-12 throughout."""
    body = plate_body()
    sk = body.skins
    A = float(sk.column_volume()[0] / H)
    T = 950.0
    rho_l = float(sk.density(np.array([T]))[0])
    v = 0.1e-3 * A
    M0 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    sk._append_bottom([0], rho_l * v, rho_l * v * float(sk.props.enthalpy(np.array([T]))[0]))
    z0 = surface_height(body)
    tally = skinbody._Tally(sk.n)
    body._layer(sk.cells_volume() - sk.column_volume(), tally, "part")
    z1 = surface_height(body)
    M1 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    body._end_step(tally, body._samples())
    z2 = surface_height(body)
    M2 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    want = (1.0 - rho_l / sk.rho_s) * v / A
    added = rho_l * v
    print("molten cell returned: phi rose {:.4f} mm with it, fell {:.4f} mm by the step end, the densities give {:.4f} mm "
          "(1 - rho_l/rho_s = {:.3f}); books {:.1e} {:.1e}".format((z1 - z0) * 1e3, (z1 - z2) * 1e3, want * 1e3,
                                                                  1.0 - rho_l / sk.rho_s, (M1 - M0 - added) / M0,
                                                                  (M2 - M0 - added) / M0))
    assert z1 - z0 == pytest.approx(v / A, rel=1e-6)
    assert z1 - z2 == pytest.approx(want, rel=0.05)
    assert abs(M1 - M0 - added) <= 1e-12 * M0 and abs(M2 - M0 - added) <= 1e-12 * M0
    assert tally.handed > 0.0


def test_a_remesh_carries_phi_and_rebuilds_beta(monkeypatch):
    """The step end's remesh (the controller's ruling T2-b): phi carried to the new mesh -- interpolated, reset to a
    distance keeping its contour and corrected per piece -- holds each piece's volume to 1e-12, beta is skin_base(phi)
    on the new mesh, the columns are full and the books close. The test cube's mesh rebuilt to a target three quarters
    of its own (`_needs_rebuild` held to say so once)."""
    body = cube_body()
    label, n = pieces(body.mesh, body.phi)
    V0 = np.array([inside_volumes(body.mesh, body.phi)[label == p].sum() for p in range(n)])
    M0 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    target = body._target
    monkeypatch.setattr(body, "_target", lambda *a, **k: (lambda f: lambda x, s: 0.75 * np.asarray(f(x, s)))(
        target(*a, **k)))
    monkeypatch.setattr(body, "_needs_rebuild", lambda target, region: True)
    nodes = len(body.mesh.points)
    tally = skinbody._Tally(body.skins.n)
    same = body._rebuild(tally, body._samples())
    label, n1 = pieces(body.mesh, body.phi)
    V1 = np.array([inside_volumes(body.mesh, body.phi)[label == p].sum() for p in range(n1)])
    M1 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    print("remesh: {} nodes then {}; piece volumes {:.1e} apart; books {:.1e}; closure {:.1e}".format(
        nodes, len(body.mesh.points), np.abs(V1 / V0 - 1).max(), M1 / M0 - 1, closure(body)))
    assert not same and len(body.mesh.points) != nodes and n1 == n
    np.testing.assert_allclose(V1, V0, rtol=1e-12, atol=0.0)
    assert np.array_equal(body.beta, skin_base(body.mesh, body.phi, H, band=body.band))
    assert tally.fill_miss <= skinbody.FILL_RTOL and closure(body) <= CLOSE_RTOL
    assert abs(M1 / M0 - 1) <= 1e-12
    assert [e["kind"] for e in tally.corrections] == ["remesh"]


SPECK_AWAY = 7.5e-3       # m: x of the point near which the remesh's speck is made (the cube's +x face is at 5.13 mm)


def test_a_remesh_that_parts_off_a_markerless_speck_lifts_it(monkeypatch):
    """The controller's ruling T7b-b (plan 03b Task 7b, finding B): the speck rule extends to the step end's remesh. phi
    carried to the new mesh (interpolated and reset) can part off a piece no marker stands on: the stage 3 box's second
    step met one of 2.3e-16 m3 in its unreceded rim, which the remesh's rule (the pieces share the old metal by their
    volumes) kept and the division then stopped on (a piece with layer and no column, ruling T7-g). By `_move_phi`'s
    criterion (`_piece_targets`: no column maps to it and it holds none of the volume's metal) such a piece is lifted to
    PIECE_GONE before the remesh's correction, its metal the largest piece's, and its nodes are left out of the
    correction's claims. The test cube's mesh rebuilt to three quarters of its target (as in
    `test_a_remesh_carries_phi_and_rebuilds_beta`), its reset phi given a speck SPECK_DEPTH deep at a node in the gas
    near (SPECK_AWAY, OFF, OFF), a piece of its own: the step end's rebuild divides without stopping, phi is one piece
    holding the old metal to 1e-12, the speck is counted as used up and the books close."""
    from reentry_model.levelset.stepend import StepEnd
    body = cube_body()
    V0 = float(inside_volumes(body.mesh, body.phi).sum())
    M0 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    target = body._target
    monkeypatch.setattr(body, "_target", lambda *a, **k: (lambda f: lambda x, s: 0.75 * np.asarray(f(x, s)))(
        target(*a, **k)))
    monkeypatch.setattr(body, "_needs_rebuild", lambda target, region: True)
    old_mesh, made = body.mesh, []

    def reset(mesh, phi):
        out = StepEnd._reset(body, mesh, phi)
        if mesh is not old_mesh:                         # the remesh's phi: a speck in the gas, a piece of its own
            far = np.flatnonzero(out > 1.5e-3)
            for k in far[np.argsort(np.linalg.norm(mesh.points[far] - [SPECK_AWAY, OFF, OFF], axis=1))]:
                trial = out.copy()
                trial[k] = -SPECK_DEPTH
                if pieces(mesh, trial)[1] == 2:
                    made.append(k)
                    return trial
        return out
    monkeypatch.setattr(body, "_reset", reset)
    tally = skinbody._Tally(body.skins.n)
    body._rebuild(tally, body._samples())
    label, n = pieces(body.mesh, body.phi)
    V1 = float(inside_volumes(body.mesh, body.phi).sum())
    M1 = sum(body.books()[k] for k in ("mass_volume", "mass_skins"))
    print("remesh speck at node {} ({} mm): pieces after {}; phi there {:.1e} m; volume {:.1e}; books {:.1e}".format(
        made, np.round(body.mesh.points[made[0]] * 1e3, 2), n, body.phi[made[0]], V1 / V0 - 1.0, M1 / M0 - 1.0))
    assert made and n == 1 and body.phi[made[0]] > 0.0
    assert tally.pieces_used == 1
    assert V1 == pytest.approx(V0, rel=1e-12) and abs(M1 / M0 - 1.0) <= 1e-12


def test_metal_thinner_than_twice_the_skin_depth_stops_the_run():
    """Design 21.11: where two columns whose normals oppose share a boundary in the layer, the metal there is thinner
    than twice the skin depth, and the run stops with a RuntimeError that names the place. A slab 2.0 mm thick (its
    walls the mesh box's, four markers on each face), against 2H = 2.4 mm: the body raises as it is built."""
    base = BoxBase.around((-PLATE_W / 2, -PLATE_W / 2, -3.2e-3), (PLATE_W / 2, PLATE_W / 2, 3.2e-3), PLATE_H0)
    with pytest.raises(RuntimeError, match="^" + skinbody.THIN_METAL) as stop:
        skinbody.LevelSetSkinBody(aa7075(), base, lambda x: np.abs(x[:, 2] - 0.07e-3) - 1.0e-3, 850.0, SETTINGS,
                                  spacing=PLATE_H0)
    print("thin slab: {}".format(stop.value))
    assert "mm" in str(stop.value)


def test_a_marker_left_with_no_area_gives_its_skin_to_its_neighbour(monkeypatch):
    """Planning ruling 6 (the Review Focus's first item): a marker the step end's tessellation leaves no area is
    released, and its skin goes to the nearest marker that stays by depth (`Skins.combine(..., by="depth")`), not to the
    volume. Here the test cube's step end reports a merge, so that it tessellates a second time, and that tessellation
    leaves one marker no area, its skin moving at 1 m/s: the marker leaves, counted as released; the skins' mass, enthalpy
    and momentum are unchanged to 1e-14 (the volume untouched); and every marker left has an area. Released alone
    (`_release`), a marker's skin adds to the nearest marker's exactly."""
    from reentry_model.levelset.markers import Markers
    body = cube_body()
    mk, sk = body.markers, body.skins
    sk.U[mk.skin_index[0]] = np.where(sk.valid[mk.skin_index[0]][:, None], [1.0, 0.0, 0.0], 0.0)
    calls = {"n": 0}
    tessellate = Markers.tessellate

    def tess(self, mesh, phi):
        faces = tessellate(self, mesh, phi)
        calls["n"] += 1
        if calls["n"] == 2:                               # the tessellation after split and merge
            self.area[0] = 0.0
        return faces
    monkeypatch.setattr(Markers, "tessellate", tess)
    monkeypatch.setattr(Markers, "merge", lambda self, spacing: 1)
    M0, E0, P0 = (np.sum(t, axis=0) for t in sk.totals())
    n0, released0, volume0 = mk.n, mk.trapped, inside_volumes(body.mesh, body.beta).sum()
    body._settle(skinbody._Tally(sk.n))
    M1, E1, P1 = (np.sum(t, axis=0) for t in sk.totals())
    print("a marker with no area: {} markers then {}, released {}; skins' books {:.1e} {:.1e}, momentum {:.1e}".format(
        n0, mk.n, mk.trapped - released0, M1 / M0 - 1, E1 / E0 - 1, np.abs(P1 - P0).max() / np.abs(P0).max()))
    assert mk.n == n0 - 1 and mk.trapped == released0 + 1
    assert np.all(mk.area > 0.0)
    assert abs(M1 / M0 - 1) <= 1e-14 and abs(E1 / E0 - 1) <= 1e-14
    assert np.abs(P1 - P0).max() <= 1e-14 * np.abs(P0).max()
    assert inside_volumes(body.mesh, body.beta).sum() == volume0
    x0 = mk.x.copy()
    keeper = int(np.argsort(np.linalg.norm(x0 - x0[0], axis=1))[1])
    sk.U[mk.skin_index[0]] = np.where(sk.valid[mk.skin_index[0]][:, None], [0.0, 1.0, 0.0], 0.0)
    (M_g, E_g, P_g), (M_k, E_k, P_k) = ([t[mk.skin_index[i]] for t in sk.totals()] for i in (0, keeper))
    body._release(np.array([0]), skinbody._Tally(sk.n))
    k = int(np.argmin(np.linalg.norm(mk.x - x0[keeper], axis=1)))
    M_n, E_n, P_n = (t[mk.skin_index[k]] for t in sk.totals())
    assert M_n == pytest.approx(M_k + M_g, rel=1e-14) and E_n == pytest.approx(E_k + E_g, rel=1e-14)
    np.testing.assert_allclose(P_n, P_k + P_g, rtol=0.0, atol=1e-14 * np.abs(P_g).max())


def test_a_part_that_would_empty_a_column_is_split():
    """The Review Focus's second item: a part whose feed removal would take more than a column's depth is split (the
    part ratio, its move of phi over MOVE_SHARE of the spacing), and never leaves a column with a negative volume. The
    plate at 880 K with four markers 2.8 mm apart (a move of 0.7 mm allowed a part) under 10 MW/m2 with the liquid removed
    as it reaches T_feed, one 0.5 s step: more than the allowance, and more than H, would leave in one part, so the step
    takes more than one part; every column's volume and its cells' are positive, and the books close."""
    body = plate_body(T0=880.0, remove="feed", spacing=PLATE_W / 2)
    book = body.advance(0.5, q_prescribed=lambda mk: np.full(mk.n, 1e7))
    sk = body.skins
    print("emptying part: {} parts (split depth {}), removed {:.3e} kg, the face moved {:.3f} mm; residuals mass {:.1e} "
          "enthalpy {:.1e}".format(book.parts, book.split_depth_max, book.mass_removed, -surface_height(body) * 1e3,
                                   book.mass_residual, book.enthalpy_residual))
    assert book.mass_removed > body.rho * float(sk.column_volume().sum())  # more than a column's depth left
    assert book.parts > 1
    assert np.all(sk.column_volume() > 0.0) and np.all(sk.cells_volume() > 0.0)
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8


EDGE_T = (890.0, 1000.0)   # K: the cube, rigid just below T_rigid (895.1 K), and its +x side's columns, liquid
EDGE_PULL = 0.05           # m/s2: the acceleration the metal feels along +z in the runoff-over-an-edge test
EDGE_STEP = 0.01           # s: its one macro step (one heat sub-step)


def _shuffle_rows(body, seed=3):
    """The body's skins' rows put in a random order, everything kept per row (the cells, the top areas and curvatures,
    the profile, the feed's momentum and the carry), so that `Markers.skin_index` is a real mapping (the body keeps the
    rows in the markers' order otherwise)."""
    sk, mk = body.skins, body.markers
    perm = np.random.default_rng(seed).permutation(sk.n)          # new row r holds old row perm[r]
    for name in ("m", "H", "U", "count", "area", "kappa", "momentum_out", "profile"):
        setattr(sk, name, getattr(sk, name)[perm])
    body.carry = body.carry[perm]
    mk.skin_index = np.argsort(perm)[mk.skin_index]
    assert not np.array_equal(mk.skin_index, np.arange(mk.n))


def test_runoff_crosses_a_box_edge_through_the_body(monkeypatch):
    """Plan 03b Task 6 (design 21.5: runoff crosses the shared boundaries of the one division of the layer): the test
    cube at 890 K (rigid) with its skins' rows shuffled (`_shuffle_rows`), the columns on its +x side (normals within 72
    degrees of +x) made liquid at 1000 K through their depth, and one 0.01 s step with the metal pulled along +z
    (EDGE_PULL, a deceleration with no gas state), toward the +x face's edge with the +z face, its liquid's surface
    tension off (the capillary push of phi's recovered curvature at the sharp cube's edges, some 2.5e5 N/m3 on its face
    markers, folds the liquid face within the step and stops it on THIN_METAL: plan 03b Task 6's diagnosis of the
    spraying fixture, the controller's carry T5-h). Every flow sub-step runs on the body's division (each transport's
    sheet depth-resolved, `runoff.Sheet.from_division`, its rows the shuffled ones); metal arrives in the +z face's
    columns across the edge (their net gain over the real passes' transports is positive); and the step's books close
    (mass 1e-12, enthalpy 1e-8)."""
    from dataclasses import replace
    from reentry_model.levelset import runoff
    body = cube_body(T0=EDGE_T[0])
    body.material.liquid = replace(body.material.liquid, sigma=0.0)
    _shuffle_rows(body)
    sk, mk = body.skins, body.markers
    hot = mk.skin_index[mk.normal[:, 0] > 0.3]
    T = np.where(sk.valid, EDGE_T[1], np.nan)
    sk.H[hot] = np.where(sk.valid[hot], sk.m[hot] * sk.props.enthalpy(T[hot]), 0.0)
    gained, sheets = np.zeros(sk.n), []
    transport = runoff.transport

    def spy(skins, sheet, dt, faces=None):
        book = transport(skins, sheet, dt, faces=faces)
        sheets.append(sheet)
        if skins.dH is None:                             # a real pass (a trial pass carries its sensitivity)
            gained[:] += book.mass
        return book
    monkeypatch.setattr(runoff, "transport", spy)
    rows_z = mk.skin_index[mk.normal[:, 2] > 0.9]          # the +z face's rows within the step (the markers settle after)
    book = body.advance(EDGE_STEP, decel_vec=(0.0, 0.0, EDGE_PULL))
    print("runoff over the edge: {} transports; the +z face's columns gained {:.3e} kg over the real passes, the most "
          "one gained {:.3e} kg; residuals mass {:.1e}, enthalpy {:.1e}; {} parts".format(
              len(sheets), gained[rows_z].sum(), gained[rows_z].max(), book.mass_residual, book.enthalpy_residual,
              book.parts))
    assert sheets and all(s.length_z is not None for s in sheets)
    assert all(not np.array_equal(s.marker, np.arange(s.n)) for s in sheets)   # the rows are not the markers' order
    assert gained[rows_z].sum() > 0.0
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8


def test_pending_load_reaches_the_volume():
    """A deferred load of 5 J at a point under the +x face is delivered to the volume, at the patches of beta nearest the
    point, in the next part (a 0.05 s step): the body's stored enthalpy rises by the 5 J less the load still deferred
    (to 1e-9), the volume receives most of it, the skins above the warmed base take the rest, and the volume's hottest
    node lies near the point; the load is spread over the group of the marker whose base is nearest the point, never on
    one patch alone."""
    body = cube_body()
    point = np.array([HALF + OFF - 1.6e-3, OFF, OFF])
    body.add_pending(point, 5.0)
    before = body.books()
    book = body.advance(0.05)
    gained = book.enthalpy_volume + book.enthalpy_skins - before["enthalpy_volume"] - before["enthalpy_skins"]
    print("pending load: volume gained {:.6f} J, skins {:.6f} J, still deferred {:.3e} J".format(
        book.enthalpy_volume - before["enthalpy_volume"], book.enthalpy_skins - before["enthalpy_skins"],
        book.pending_J))
    assert gained + book.pending_J == pytest.approx(5.0, rel=1e-9)
    assert book.enthalpy_volume - before["enthalpy_volume"] > 2.5
    assert book.enthalpy_skins - before["enthalpy_skins"] >= 0.0
    T = body.solver.T
    hot = np.argmax(np.where(body.solver.active, T, -np.inf))
    assert np.linalg.norm(body.mesh.points[hot] - point) < 3e-3
    body.add_pending(point, 5.0)
    S, A_group, base = body._interface()
    spread = body._pending_per_patch(S, base, A_group > 0.0)
    assert spread.sum() == pytest.approx(body.pending_J.sum(), rel=1e-12)
    assert (spread != 0.0).sum() > 1


def test_a_step_at_a_gas_state_with_the_flow_and_the_spray_on_closes_its_books():
    """The cube at a gas state (77.5 km, 7.5 km/s, the freestream density raised to 1e-4 kg/m3 so that it melts within
    the step) with the shared heating and surface flow (`loads.GasLoads`), the runoff driven by the surface flow and the
    deceleration, radiation, and Girin's plan (the shared SprayModel): one 0.1 s macro step closes mass and enthalpy,
    the plan's books hold exactly what the step booked as sprayed, the history row carries the droplet columns, and
    what a trial pass is undone to restores the cells, the flow memory, the profile and the plan's books."""
    body = gas_body()
    started = time.perf_counter()
    book = body.advance(0.1, state=gas_state(rho=1e-4))
    print("gas state: {:.1f} s, {} parts, {} passes; heat in {:.3f} J, radiated {:.3e} J, sprayed {:.3e} kg; residuals "
          "mass {:.2e}, enthalpy {:.2e}; momentum to the wall {}".format(
              time.perf_counter() - started, book.parts, book.passes, book.enthalpy_in, book.enthalpy_radiated,
              book.mass_sprayed, book.mass_residual, book.enthalpy_residual, book.momentum_to_wall))
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert book.enthalpy_in > 0.0 and book.enthalpy_radiated > 0.0
    assert max(owed(body, book)) <= OWED_SHARE
    plan = body.plan
    assert plan is not None
    assert plan.mass.sum() + plan.mass_slurry.sum() == pytest.approx(book.mass_sprayed, rel=1e-12, abs=1e-30)
    row = book.row()
    assert row["altitude_km"] == pytest.approx(77.5) and np.isfinite(row["flowing_depth_mean_mm"])
    ctx = {"plan": plan}
    saved = body._save(ctx)
    sk = body.skins
    m, memory, mass, profile = sk.m.copy(), sk.flow_memory.copy(), plan.mass.copy(), sk.profile.copy()
    plan.mass += 1.0
    sk.flow_memory = memory + 1.0
    sk.m[0, 0] *= 2.0
    sk.profile[0] *= 2.0
    body._restore(saved, ctx)
    assert np.array_equal(plan.mass, mass) and np.array_equal(sk.flow_memory, memory)
    assert np.array_equal(sk.m, m) and np.array_equal(sk.profile, profile)


def gas_body(phi0=cube):
    """The cube (`phi0`) with the shared heating and surface flow (`loads.GasLoads`) and Girin's spray model, radiating."""
    from reentry_model.levelset import loads
    from reentry_model.shared import heating, spray, surface_flow
    from reentry_model.shared.material import LiquidProperties
    mat = aa7075()
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), mat.liquid.rho)
    model = spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), spray_gates="2026-11-10")
    return cube_body(phi0=phi0, loads=gl, spray_model=model, T_ambient=200.0, emissivity=0.2)


PART_BOUND = 100      # parts a flow-on step of the test cube may take at most: a tenth of PARTS_MAX
CONTINUUM = dict(rho=1e-4, kn=1e-3)   # 7.5 km/s at a density of 1e-4 kg/m3 with a body Knudsen number of 1e-3: the shock
                                      # layer's branch, where Girin's closure gives a conjugate depth


SPURIOUS_PAIRS = ("plan 03b Task 7b (rulings T7-f, T7b-a and T7-k), measured once: the runaway split is gone (the speck "
                  "that stranded two columns stays out), and the second 0.05 s step runs its parts but stops at its end "
                  "on stepend.THIN_METAL with 2 pairs of opposed columns (n_i . n_j -0.56 and -0.51) whose markers are "
                  "2.90 and 3.24 mm apart with metal all the way between them (2.89 and 3.23 mm of it along the "
                  "straight line, reaching 0.67 and 0.75 mm below phi): not metal thinner than twice the skin depth "
                  "(2.4 mm), the spurious pairs of ruling T7-k, Task 7c's subject")


class ThinMetalStop(RuntimeError):
    """A step's THIN_METAL stop (`stepend.THIN_METAL`) and no other, for an expected failure to name."""


@pytest.fixture(scope="module")
def sprayed():
    """The cube with its edges rounded to ROUNDING (`rounded_cube`, the controller's ruling T6-d: the fixture tests the
    spraying and the books, not sharp-edge capillarity; on the sharp cube the capillary push of phi's recovered
    curvature at the edges folded the liquid face and the step stopped on THIN_METAL, plan 03b Task 6's diagnosis) with
    liquid formed on its +x face by 0.25 s of 6 MW/m2 (no flow, nothing removed), then one 0.05 s macro step at the
    continuum gas state with the runoff and Girin's plan on: (body, the melting step's book, the spraying step's book or
    None, its time [s], the RuntimeError that stopped it or None)."""
    body = gas_body(rounded_cube)
    melt = body.advance(0.25, q_prescribed=one_face(6e6))
    started = time.perf_counter()
    try:
        book, stop = body.advance(0.05, state=gas_state(**CONTINUUM)), None
    except RuntimeError as err:
        book, stop = None, err
    return body, melt, book, time.perf_counter() - started, stop


def test_liquid_formed_on_the_face_sprays_and_the_books_close(sprayed):
    """Fix round 1 of Task 9 (ruling b): liquid formed on the +x face sprays in the next 0.05 s macro step at the
    continuum gas state: the released mass is positive and equals the plan's books and the history's, the mass and
    enthalpy close, and the heat still owed and the interface mismatch stay within OWED_SHARE of the heat in."""
    body, melt, book, seconds, stop = sprayed
    assert stop is None, stop
    plan = body.plan
    liquid = melt.summary["liquid_depth_mean_mm"]
    print("spraying face: liquid {:.3f} mm mean before; {:.1f} s, {} parts; sprayed {:.3e} kg from {} markers; heat in "
          "{:.3f} J; residuals mass {:.2e}, enthalpy {:.2e}; droplets {:.0f}, median radius {:.1f} um".format(
              liquid, seconds, book.parts, book.mass_sprayed, int((plan.rate > 0).sum()), book.enthalpy_in,
              book.mass_residual, book.enthalpy_residual, book.summary["n_released"], book.summary["r_median_um"]))
    assert liquid > 0.0 and book.mass_sprayed > 0.0
    assert plan.mass.sum() + plan.mass_slurry.sum() == pytest.approx(book.mass_sprayed, rel=1e-12)
    assert book.summary["released_step_mass_kg"] == pytest.approx(book.mass_sprayed, rel=1e-12)
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert max(owed(body, book)) <= OWED_SHARE


@pytest.mark.xfail(raises=ThinMetalStop, strict=True, reason=SPURIOUS_PAIRS)
def test_two_flow_on_steps_keep_their_part_count_bounded(sprayed):
    """Fix round 1 of Task 9 (ruling a), the smallest case found that reproduced the skyscrapers' runaway: the spraying
    step above and a second 0.05 s step at the same state. A part that meets its allowance only after splits, or never,
    carries the splits' shortness into the skyscrapers' estimate; here a step's hint is at most 1.25 times the parts it
    took: both steps stay within PART_BOUND parts, the second step's hint is at most the skyscrapers' for it, both hints
    are 1.25 times their step's parts at most, and the books close. It fails (SPURIOUS_PAIRS): the second step stops on
    THIN_METAL, which the test raises as ThinMetalStop (any other stop is a failure)."""
    body, _, first, _, stop = sprayed
    assert stop is None, stop
    assert first.parts_next <= max(1, int(np.ceil(skinbody.MOVE_HINT_MARGIN * first.parts)))
    try:
        second = body.advance(0.05, state=gas_state(**CONTINUUM))
    except RuntimeError as err:
        if str(err).startswith(skinbody.THIN_METAL):
            raise ThinMetalStop(str(err)) from err
        raise
    print("two flow-on steps: parts {} then {}; hints {} then {}; the skyscrapers' hints {} then {}".format(
        first.parts, second.parts, first.parts_next, second.parts_next, first.parts_skyscrapers,
        second.parts_skyscrapers))
    assert first.parts <= PART_BOUND and second.parts <= PART_BOUND and second.parts_next <= second.parts_skyscrapers
    assert second.parts_next <= max(1, int(np.ceil(skinbody.MOVE_HINT_MARGIN * second.parts)))
    assert abs(second.mass_residual) <= 1e-12 and abs(second.enthalpy_residual) <= 1e-8
    assert max(owed(body, second)) <= OWED_SHARE


FACE_GROUP_MEDIAN = (0.95, 1.05)   # the test cube's +x face skins' groups over their base areas: the median within this
FACE_GROUP_MIN = 0.7               # and none below this
FACE_GROUP_OWN = 0.8               # at least this share of them left by the guard, each its own patches in proportion


@pytest.fixture(scope="module")
def rest():
    """The test cube at rest, as built (read only)."""
    return cube_body()


def test_a_skin_s_group_is_the_base_beneath_its_own_cell(rest):
    """Design 21.5 (one division serves the base groups): each patch of beta belongs to the marker owning its foot on
    phi (`_beta_owners`, the division's base level), so a skin's base heat enters the volume through the base of its own
    column: the patches' areas are divided whole, the test cube's +x face skins' groups hold their columns' bases (the
    division's A(H)) to a median within FACE_GROUP_MEDIAN and none below FACE_GROUP_MIN, at least FACE_GROUP_OWN of
    them are left by the guard and spread their law over their own patches in proportion to the area they hold of
    each, and the law stays exact: the groups' shares sum to one per marker with a base."""
    body = rest
    P = body.solver.P
    G = body._beta_owners(P)
    np.testing.assert_allclose(np.asarray(G.sum(axis=0)).ravel(), P.area, rtol=1e-12)
    S, A_group, _ = body._interface()
    mk = body.markers
    A = body.skins.base_area()
    face = mk.skin_index[mk.normal[:, 0] > 0.9]
    ratio = A_group[face] / A[face]
    print("face groups {} (guard replaced {})".format(np.round(np.sort(ratio), 3), body.last_guarded))
    lo, hi = FACE_GROUP_MEDIAN
    assert lo <= np.median(ratio) <= hi and ratio.min() >= FACE_GROUP_MIN
    has = A_group > 0.0
    np.testing.assert_allclose(np.asarray(S.sum(axis=1)).ravel()[has], 1.0, rtol=1e-14)
    held = np.asarray(G.sum(axis=1)).ravel()
    g_lo, g_hi = skinbody.GROUP_AREA_RANGE
    own = face[(held[face] >= g_lo * A[face]) & (held[face] <= g_hi * A[face])]     # the groups the guard left
    print("{} of {} face groups left by the guard".format(len(own), len(face)))
    assert len(own) >= FACE_GROUP_OWN * len(face)
    want = G[own].multiply(1.0 / held[own][:, None]).toarray()
    np.testing.assert_allclose(S[own].toarray(), want, rtol=0.0, atol=1e-12)        # their own patches, in proportion


def test_each_patch_of_beta_is_owned_as_the_division_owns_its_base_level(rest):
    """Plan 03b Task 5's fix round 1 (the controller's ruling T5-j): beta keeps phi + H's values in the cells it cuts,
    so its patches are the division's base level, and the groups are the division's ownership of it, found once per
    division. On the test cube at rest, with the division computed afresh on phi as it is (so the body's own is
    current): every patch of beta lies in a cell the base level cuts, with that level's area there to 1e-12, and is held
    by exactly the markers that own that level's patch, each the same area to 1e-12 of the patch's."""
    from reentry_model.levelset.division import divide
    body = rest
    mk, P = body.markers, body.solver.P
    d = divide(body.mesh, body.phi, H, mk.x, mk.normal)
    at = {int(c): k for k, c in enumerate(d.base_cell)}
    assert all(int(c) in at for c in P.cell)
    k = np.array([at[int(c)] for c in P.cell])
    np.testing.assert_allclose(P.area, d.base_patch_area[k], rtol=1e-12, atol=0.0)
    W = d.base_owner[:, k].toarray()
    want = np.zeros((body.skins.n, len(P.area)))
    want[mk.skin_index] = W
    G = body._beta_owners(P).toarray()
    print("{} patches of beta, {} owner entries; largest difference {:.1e} of a patch's area".format(
        len(P.area), int((G > 0.0).sum()), (np.abs(G - want) / P.area).max()))
    assert np.array_equal(G > 0.0, want > 0.0)
    assert (np.abs(G - want) / P.area).max() <= 1e-12


BEAM_Z = np.array([0.0, 0.0, 1.0])   # a beam arriving along -z


def test_the_skins_take_a_beam_over_the_area_of_phi_they_own():
    """Task 11's fix round 1, C1 (the controller's ruling): one partition of phi. Each marker's area is the area of phi
    nearer it than any other marker, the area the gas loads average over, so that on the test cube under a beam along -z
    (each patch's projection max(n . BEAM_Z, 0) as its flux) the heat the skins take, each its own flux times its own
    area, is the beam times phi's projected area to round-off, and the skins' areas sum to phi's."""
    from reentry_model.levelset import loads
    from reentry_model.levelset.cutcell import patches
    body = cube_body()
    mk, sk = body.markers, body.skins
    P = patches(body.mesh, body.phi)
    proj = np.maximum(P.normal @ BEAM_Z, 0.0)
    taken = (sk.area[mk.skin_index] * loads.patch_mean(mk, P, proj)).sum()
    want = (P.area * proj).sum()
    print("the skins take {:.15f} of the beam's projected heat; their areas {:.15f} of phi's".format(
        taken / want, sk.area.sum() / P.area.sum()))
    assert taken == pytest.approx(want, rel=1e-12)
    assert sk.area.sum() == pytest.approx(P.area.sum(), rel=1e-12)


def test_a_marker_whose_nearest_patches_are_few_gets_a_grown_group_and_the_law_stays_exact():
    """Fix round 2 of Task 9 (the ruling on concern 1): a marker's base law is spread over a group of patches of beta
    whose area is within GROUP_AREA_RANGE of its skin's base area. Flat patches of 0.1 x 0.1 mm on a plane carry five
    base points of base area 1 mm2 each, one at the origin ringed by four 0.02 mm away, so that by nearest base it owns
    almost nothing. Every group ends within range, the crowded marker's grown from its nearest patches, and the law each
    patch receives sums back, over the patches, to every marker's a and b exactly (to round-off); a row marked insulated
    (a column with no base, plan 03b Task 5) takes no group."""
    h = 0.1e-3
    ix, iy = np.meshgrid(np.arange(-15, 26), np.arange(-15, 16), indexing="ij")
    centroid = np.column_stack([(ix.ravel() + 0.5) * h, (iy.ravel() + 0.5) * h, np.zeros(ix.size)])
    area = np.full(len(centroid), h * h)
    d = 2e-5
    base = np.array([[0.0, 0.0, 0.0], [d, 0.0, 0.0], [-d, 0.0, 0.0], [0.0, d, 0.0], [0.0, -d, 0.0]])
    A_base = np.full(len(base), 1e-6)
    near = skinbody.patch_groups(centroid, area, base, np.zeros(len(base)))[1]      # nearest base alone, no target
    S, A_group = skinbody.patch_groups(centroid, area, base, A_base)
    ratio = A_group / A_base
    print("patch groups: nearest-base areas {} mm2, grown {} mm2".format(np.round(near * 1e6, 3),
                                                                          np.round(A_group * 1e6, 3)))
    assert near[0] < 0.1 * A_base[0]
    lo, hi = skinbody.GROUP_AREA_RANGE
    assert np.all(ratio >= lo) and np.all(ratio <= hi)
    np.testing.assert_allclose(np.asarray(S.sum(axis=1)).ravel(), 1.0, rtol=1e-14)
    a, b = np.array([3.0, -2.0, 5.0, 1.0, 0.5]), np.array([-16.0, -0.4, -0.1, -2.0, -1.0])
    q0, q1 = (S.T @ a) / area, (S.T @ b) / area
    assert (q0 * area).sum() == pytest.approx(a.sum(), rel=1e-14)
    assert (q1 * area).sum() == pytest.approx(b.sum(), rel=1e-14)
    off = np.array([True, False, False, False, False])
    S_off, A_off = skinbody.patch_groups(centroid, area, base, A_base, insulated=off)
    assert A_off[0] == 0.0 and S_off[0].nnz == 0 and np.all(A_off[1:] > 0.0)


def _hints(move, hint, steps=5, dt=0.5):
    """The part-count hints over consecutive macro steps of length dt [s] driven by the body's own splitting and hint
    rules (`skinbody.split_part`, `skinbody.next_hint`) with a synthetic move: `move(t0, t1)` [-] the largest move
    ratio of a part from t0 to t1 [s] into the step. Returns the hints, the first being the one the steps start with."""
    out = [hint]
    for _ in range(steps):
        clock = {"t": 0.0}

        def run(h):
            t0 = clock["t"]

            def accept():
                clock["t"] = t0 + h
            return move(t0, t0 + h), accept, lambda: None
        counts = [skinbody.split_part(run, dt / hint, dt, 0 if hint == 1 else 1) for _ in range(hint)]
        assert clock["t"] == pytest.approx(dt, rel=1e-12)
        hint = skinbody.next_hint(counts)
        out.append(hint)
    return out


def test_a_move_that_does_not_shrink_with_the_part_does_not_raise_the_part_count():
    """A jump of three allowances at each step's start and a refill of a twentieth of the allowance in every part, neither
    shrinking with the part's length (the full-size cube's hand-down of a backlog, fix round 5): splitting cannot help,
    the parts it produced do not count, and the hint falls to one and stays there."""
    hints = _hints(lambda t0, t1: (3.0 if t0 == 0.0 else 0.0) + 0.05, 4)
    print("fixed jump: hints", hints)
    assert all(b <= a for a, b in zip(hints, hints[1:]))
    assert hints[-1] == 1


def test_a_move_that_scales_with_the_part_sets_the_part_count_and_holds_it():
    """5.5 allowances of move per macro step, spread evenly: from one part the hint grows to ceil(1.25 x 5.5) = 7, what
    the move asks for, and holds there."""
    hints = _hints(lambda t0, t1: 5.5 * (t1 - t0) / 0.5, 1)
    print("scaling move: hints", hints)
    assert hints[1:] == [7] * 5


def test_small_moves_let_the_part_count_fall():
    """A fifth of an allowance per macro step: from 20 parts the hint falls to one and stays."""
    hints = _hints(lambda t0, t1: 0.2 * (t1 - t0) / 0.5, 20)
    print("small moves: hints", hints)
    assert hints[1:] == [1] * 5


def test_the_snap_books_the_markers_its_last_round_brings_from_afar_and_their_slide():
    """Task 9b (minor, left from Task 9): phi = 0.05 (|x| - 2 mm), rising only 0.05 per metre along its normal, and
    three markers 1 mm outside its zero contour: each recovered-normal step removes a twentieth of the distance, so 19
    rounds leave them 1 mm x 0.95^38 = 0.14 mm off, where |phi| = 7e-6 is beyond the Newton steps' reach (1e-6); the last
    round brings all three on by Newton steps from afar. The step's tally counts the three (`snap_fallback`) and holds
    the largest sideways part of those Newton steps (`snap_slide`), as a replay of the rounds measures them, and every
    marker ends on phi."""
    from reentry_model.levelset import markers as M
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.markers import Markers
    from reentry_model.levelset.stepend import SNAP_ROUNDS, SNAP_TOL, StepEnd

    class Host(StepEnd):
        pass
    mesh = build(BoxBase.around((-4e-3,) * 3, (4e-3,) * 3, 1e-3))
    phi = 0.05 * (np.linalg.norm(mesh.points, axis=1) - 2e-3)
    u = np.array([[1.0, 0.2, 0.1], [-0.3, 1.0, 0.4], [0.2, -0.5, -1.0]])
    u /= np.linalg.norm(u, axis=1)[:, None]
    x0 = 3e-3 * u

    def markers():
        return Markers(x0, u, np.ones(3), np.zeros((3, 2)), np.arange(3), np.full(3, -1), 1e-3)
    host = Host()
    host.mesh, host.phi, host.markers = mesh, phi, markers()
    tally = skinbody._Tally(0)
    host._snap(tally)
    mk = markers()
    for _ in range(SNAP_ROUNDS - 1):
        mk.snap(mesh, phi)
    x = mk.x
    for _ in range(M.SNAP_PASSES):
        val, n = M._field(mesh, phi, x)
        x = x - val[:, None] * n
    val, _ = M._field(mesh, phi, x)
    far = int((np.abs(val) > M.NEWTON_REACH).sum())
    after = M._cell_newton(mesh, phi, x)
    _, n = M._field(mesh, phi, after)
    d = after - x
    slide = np.linalg.norm(d - np.einsum("ij,ij->i", d, n)[:, None] * n, axis=1).max()
    off = np.abs(M._field(mesh, phi, host.markers.x)[0]).max()
    print("snap from afar: {} markers by the fallback, largest slide {:.3e} m (replay {:.3e} m), off phi {:.1e}".format(
        tally.snap_fallback, tally.snap_slide, slide, off))
    assert far == 3 and tally.snap_fallback == 3
    assert tally.snap_slide == slide and slide > 0.0
    assert off <= SNAP_TOL


def test_the_test_cube_s_markers_keep_their_neighbours():
    """Task 9b (the ruling on the neighbour test), on the division since plan 03b Task 6: on the test cube as the body
    builds it (its edges sharp since plan 03b Task 5; markers one per 2 mm voxel) every column's top cross-section in
    the body's division (its A(0), which the skins fill and the runoff's cells cross) is within a factor of
    OWNERSHIP_FACTOR (1.5) of the area of phi nearest its marker. (Task 5 widened the tangent-plane cells' bound to 2.5:
    their areas were 0.47 to 1.23 of the nearest-marker areas on the sharp cube, measured once.)"""
    body = cube_body()
    mk = body.markers
    r = nearest_ratio(mk, body.mesh, body.phi, body.division)
    print("test cube: {} markers; the division's top cross-sections {:.3f} to {:.3f} times the nearest-marker "
          "areas".format(mk.n, r.min(), r.max()))
    assert r.min() >= 1.0 / OWNERSHIP_FACTOR and r.max() <= OWNERSHIP_FACTOR


# -- the step end's re-partition (Task 9b) -----------------------------------------------------------------------------
BOOKS_RTOL = 1e-14    # the skins' mass and enthalpy over a step end, relative (Task 9b)
SLAB_HALF = (6e-3, 6e-3, 2e-3)    # m: the flat box's half sides; its top face is the flat surface the markers move on


def _slab_host(mat, T=850.0):
    """A host of the step end (`stepend.StepEnd`) on a closed flat box (a slab 12 x 12 x 4 mm, its centre OFF off the
    grid; 1 mm cells): on its top face 36 markers on a regular grid 2 mm apart (x, y = OFF -5, -3, ..., 5 mm), on its
    other faces markers seeded one per 2 mm voxel as the body seeds them. beta is skin_base(phi); the layer is divided
    and each skin filled with its column of solid metal at T [K]; the host is settled until a step end splits and merges
    none (at most five times). Returns (host, the top grid's markers (6, 6) by x then y)."""
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.cutcell import patches
    from reentry_model.levelset.markers import Markers
    from reentry_model.levelset.stepend import StepEnd

    class Host(StepEnd):
        pass
    h = 2e-3
    mesh = build(BoxBase.around((-8e-3, -8e-3, -4e-3), (8e-3, 8e-3, 4e-3), 1e-3))
    phi = np.max(np.abs(mesh.points - OFF) - np.array(SLAB_HALF), axis=1)
    P = patches(mesh, phi)
    side = np.flatnonzero(P.normal[:, 2] < 0.5)
    key = np.floor((P.centroid[side] + 8e-3) / h).astype(np.int64)
    order = np.lexsort((-P.area[side], key[:, 2], key[:, 1], key[:, 0]))
    k = key[order]
    seeds = side[np.sort(order[np.r_[True, np.any(k[1:] != k[:-1], axis=1)]])]
    g = OFF + (np.arange(6) - 2.5) * h
    X, Y = np.meshgrid(g, g, indexing="ij")
    grid = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, OFF + SLAB_HALF[2])])
    x = np.concatenate([grid, P.centroid[seeds]])
    n = len(x)
    normal = np.concatenate([np.tile([0.0, 0.0, 1.0], (len(grid), 1)), P.normal[seeds]])
    sk = Skins(mat, np.ones(n), np.zeros((n, 2)), SETTINGS)
    mk = Markers(x, normal, np.ones(n), np.zeros((n, 2)), np.arange(n), np.full(n, -1), h, skins=sk)
    mk.tessellate(mesh, phi)
    host = Host()
    host.mesh, host.phi, host.markers, host.skins, host.spacing = mesh, phi, mk, sk, h
    host.band, host.fine_margin, host.material, host.rho, host.settings = 4e-3, 2e-3, mat, mat.rho, SETTINGS
    host.beta = skin_base(mesh, phi, H, band=host.band)
    host.carry, host.division = np.zeros(n), None
    host._sync()
    host._divide(skinbody._Tally(0))
    sk.fill(sk.density(np.full(n, T)) * sk.column_volume(), T)
    for _ in range(5):
        changes = mk.splits + mk.merges
        host._settle(skinbody._Tally(0))
        if mk.splits + mk.merges == changes:
            break
    else:
        raise RuntimeError("the slab's markers did not settle")
    host._divide(skinbody._Tally(0))
    sk.fill(sk.density(np.full(sk.n, T)) * sk.column_volume(), T)
    sk.regrid()
    host.carry = np.zeros(sk.n)
    top = np.array([[np.flatnonzero(np.all(np.isclose(mk.x, p, rtol=0.0, atol=1e-12), axis=1))[0] for p in row]
                    for row in grid.reshape(6, 6, 3)])
    return host, top


def _flat_depth(sk):
    """Each skin row's flat depth [m]: its volume, the cells' m / rho(T), over its top area (what mixing columns at one
    density averages by area)."""
    return sk.cells_volume() / sk.area


def _would_fill(host):
    """What the fill after a step end's settle would move through the base [m3] (n,): each column's volume, its
    division's afresh, less its cells' volume (`stepend.StepEnd._fill`), on the host itself (it is divided)."""
    host._divide(skinbody._Tally(0))
    return host.skins.column_volume() - host.skins.cells_volume()


THROUGH_SHARE = 1e-12  # of the inner columns' volume: what the fill after a re-partition by depth moves through the base
LEVEL_SHARE = 1e-12    # of a column's mass: each depth's metal in a new column against what its donors held there
DEPTH_SHARE = 1e-12    # of H: each new cell's depth span against its donors' (the same in every flat column here)
MUSH = (885.0, 905.0)  # K: the columns' temperatures, inside LinearMeltMaterial's mushy range (880 to 908 K), where the
                       # specific volume is linear in the enthalpy (mixing keeps volume) and the density varies with depth
                       # and column (about 2440 to 2730 kg/m3), so that a mix by mass fraction would move metal in depth


def _profiled(sk, T_of):
    """Every column's cells at the temperatures T_of(rows (n, 1), mid-depths (n, K)) [K], each cell's mass the one that
    keeps its span full (its density times the span's volume under the column's profile, `Skins._volume_at`), so that
    every column stays exactly full; the velocities kept."""
    dz = sk.dz()
    z = np.concatenate([np.zeros((sk.n, 1)), np.cumsum(dz, axis=1)], axis=1)
    span = np.diff(sk._volume_at(z), axis=1)
    T = np.where(sk.valid, T_of(np.arange(sk.n)[:, None], 0.5 * (z[:, 1:] + z[:, :-1])), MUSH[0])
    m = np.where(sk.valid, sk.density(T) * span, 0.0)
    sk.install(m, m * sk.props.enthalpy(T), sk.count, sk.U)


def test_a_repartition_by_depth_moves_no_metal_through_the_base_and_keeps_each_depth():
    """Plan 03b Task 6 (the controller's ruling P4; design 21.5: metal stays where it is, only its owner changes): mushy
    skins on a closed flat box, every column exactly full, its temperature (and so its density) varying with depth and
    from column to column (MUSH); the 16 markers in the middle of its top face shifted along it by up to 0.4 mm each, the
    surface unchanged. The step end's tessellation re-partitions the area (the shifted markers' areas change by 0.82 to
    1.16 times) and the re-partition cuts every old column by the transfer across the new division's boundaries and
    gives each piece to its new column by depth (`Skins.combine(..., by="depth")`). The skins' mass and enthalpy close
    to 1e-14 and their momentum to round-off; the stretch is one, since the surface did not change; on the top face's
    middle, where the columns are flat, the fill that follows moves nothing through the base (THROUGH_SHARE: Task 5
    measured 2.3e-3 of those columns' volume after a re-partition by mass), and each level's metal in each new column
    is what its donors held at that level, in the shares the transfer gave: every cell keeps its depth span (DEPTH_SHARE)
    and holds its donors' metal from that span (LEVEL_SHARE). (The profiles here are not flat at their last two levels,
    1.15 and 1.2 mm deep, where the slab's 1 mm cells smear the distance's kink at its mid-plane, 2 mm deep: the P1
    slope there is 0.74, so the coarea measure is 1.35 times the face's; the check is therefore cell by cell, the
    resolution at which by-depth mixing keeps the metal, not between the levels within a cell.)"""
    from materials import LinearMeltMaterial
    mat = LinearMeltMaterial()
    host, top = _slab_host(mat)
    mk, sk = host.markers, host.skins
    middle = top[1:5, 1:5].ravel()
    inner = top[2:4, 2:4].ravel()
    rng = np.random.default_rng(7)
    rng.uniform(1.5, 3.0, mk.n)                       # (the draws Task 9b's skins of one to three stubs took: the
    rng.random(len(middle)), rng.uniform(1.0, 3.0, len(middle))      # same shifts as its measurement below)
    sk.fill(sk.density(np.full(sk.n, MUSH[1])) * sk.column_volume(), MUSH[1])
    share = np.random.default_rng(8).uniform(0.0, 1.0, sk.n)    # each column's mushy profile: from MUSH[1] at its top
    lo, hi = MUSH                                                # to MUSH[0] at H, bent by its own share
    _profiled(sk, lambda r, z: hi - (hi - lo) * (z / H) ** (0.5 + share[r]))
    sk.regrid()                                       # (as the settle leaves them: its own regrid then changes nothing)
    sk.U = np.where(sk.valid[..., None], np.concatenate([rng.normal(0.0, 1.0, (sk.n, 1, 2)),
                                                         np.zeros((sk.n, 1, 1))], axis=2), 0.0)
    assert np.abs(sk.cells_volume() / sk.column_volume() - 1.0).max() <= 1e-13          # every column full
    M0, E0, P0 = (np.sum(t, axis=0) for t in sk.totals())
    P_size = np.abs(sk.m[..., None] * sk.U).sum()
    held, spans = sk.m[mk.skin_index].copy(), np.cumsum(sk.dz(), axis=1)[mk.skin_index]   # each old column's cells and
    count0 = sk.count[mk.skin_index].copy()                                                # their depths (markers' order)
    A0, n0 = mk.area.copy(), mk.n
    mk.x[middle, :2] += rng.uniform(-0.4e-3, 0.4e-3, (len(middle), 2))
    host._settle(skinbody._Tally(0))
    assert mk.n == n0
    M1, E1, P1 = (np.sum(t, axis=0) for t in sk.totals())
    ratio = mk.area[middle] / A0[middle]
    T, stretch = host.last_settle[0].transfer, host.last_settle[0].stretch
    f = T.multiply(1.0 / np.asarray(T.sum(axis=1)).ravel()[:, None]).tocsc()      # each old column's shares (j, i)
    rows = mk.skin_index[inner]
    donors = np.unique(f[:, inner].nonzero()[0])
    K = int(count0[inner[0]])
    assert np.all(count0[donors] == K) and np.all(sk.count[rows] == K)      # the same cells in every flat column here
    depth_miss = float(np.abs(np.cumsum(sk.dz(), axis=1)[rows, :K] - spans[inner, :K]).max() / H)
    want = np.asarray(f[:, inner].T @ held)[:, :K]
    level_miss = float((np.abs(sk.m[rows, :K] - want).sum(axis=1) / sk.m[rows].sum(axis=1)).max())
    miss = _would_fill(host)
    through = float(np.abs(miss[rows]).sum() / sk.column_volume()[rows].sum())
    print("re-partition by depth: shifted areas {:.3f} to {:.3f} of their old; books mass {:.1e} enthalpy {:.1e} "
          "momentum {:.1e}; stretch {:.1e} from one; the fill would move {:.1e} of the inner columns' volume through "
          "the base; their cells' depths {:.1e} of H and metal {:.1e} off their donors'".format(
              ratio.min(), ratio.max(), M1 / M0 - 1, E1 / E0 - 1, np.abs(P1 - P0).max() / P_size,
              np.abs(stretch - 1.0).max(), through, depth_miss, level_miss))
    assert ratio.min() < 0.85 and ratio.max() > 1.1                    # the shift does re-partition the area
    assert abs(M1 / M0 - 1) <= BOOKS_RTOL and abs(E1 / E0 - 1) <= BOOKS_RTOL
    assert np.abs(P1 - P0).max() <= BOOKS_RTOL * P_size
    assert np.abs(stretch - 1.0).max() <= 1e-14
    assert through <= THROUGH_SHARE
    assert depth_miss <= DEPTH_SHARE and level_miss <= LEVEL_SHARE


def test_a_thick_skin_losing_area_gives_its_metal_to_its_neighbour():
    """Task 9b: on the closed flat box every skin holds its column of solid at 850 K but one in the middle of the top
    face, overfull with 0.8 mm more of liquid at 950 K above its column. Its two neighbours along x close in on it, each
    by half a spacing, so its area halves and a quarter of its old area passes to each. Its metal follows: its own
    column keeps its depth (it does not thicken), and the heat its liquid holds above the stubs' (its enthalpy less that
    of the same mass at 850 K, which adds when columns mix) goes, column by column, exactly where the transfer sends its
    area: most of it to the two neighbours. (Before Task 9b the column kept its metal on half the area and stood twice
    as deep.)"""
    from materials import LinearMeltMaterial
    mat = LinearMeltMaterial()
    host, top = _slab_host(mat)
    mk, sk = host.markers, host.skins
    h = host.spacing
    thick, left, right = top[2, 3], top[1, 3], top[3, 3]
    row = mk.skin_index[thick]
    n_l = 40
    K = int(sk.count.max())
    m, Hc = np.zeros((sk.n, K + n_l)), np.zeros((sk.n, K + n_l))
    m[:, :K], Hc[:, :K] = sk.m, sk.H
    count = sk.count.copy()
    c = count[row]
    m[row, n_l:n_l + c], Hc[row, n_l:n_l + c] = sk.m[row, :c], sk.H[row, :c]
    m[row, :n_l] = mat.liquid.rho * sk.area[row] * 0.8e-3 / n_l
    Hc[row, :n_l] = m[row, :n_l] * mat.enthalpy(950.0)
    count[row] = c + n_l
    sk.install(m, Hc, count)

    def excess():
        """Each marker's enthalpy above the same mass at 850 K [J] (additive when columns mix)."""
        return (sk.H.sum(axis=1) - sk.m.sum(axis=1) * mat.enthalpy(850.0))[mk.skin_index]
    E_ex, M0, depth0, A0, n0 = excess(), sk.m.sum(), _flat_depth(sk)[row], mk.area[thick], mk.n
    mk.x[left, 0] += 0.5 * h
    mk.x[right, 0] -= 0.5 * h
    host._settle(skinbody._Tally(0))
    assert mk.n == n0
    depth1, after = _flat_depth(sk)[row], excess()
    print("thick skin: area {:.3f} of its old, flat depth {:.4f} of its old; its liquid's heat now {:.3f} in itself, "
          "{:.3f} and {:.3f} in its neighbours".format(mk.area[thick] / A0, depth1 / depth0, after[thick] / E_ex[thick],
                                                       after[left] / E_ex[thick], after[right] / E_ex[thick]))
    assert mk.area[thick] / A0 == pytest.approx(0.5, abs=0.02)
    assert depth1 <= depth0 * (1.0 + 1e-12) and depth1 >= 0.9 * depth0
    assert after[left] > 0.2 * E_ex[thick] and after[right] > 0.2 * E_ex[thick]
    assert abs(sk.m.sum() / M0 - 1.0) <= BOOKS_RTOL
    tr = host.last_settle[0].transfer
    share = np.asarray(tr[thick].todense()).ravel() / tr[thick].sum()
    np.testing.assert_allclose(after, E_ex[thick] * share, rtol=0.0, atol=1e-12 * E_ex[thick])
    assert share[left] + share[right] >= 0.9 * (1.0 - share[thick])


# -- heat across a box edge (plan 03b Task 6) ---------------------------------------------------------------------------
BOX_T0, BOX_FLUX, BOX_TIME = 700.0, 1.0e6, 0.5    # K, W/m2, s: the stage 3 box, the flux on its -x face, the one step
EDGE_PROBE = (0.5e-3, 0.6e-3)   # m: the probe 0.5 mm from the heated face's edge along the unheated +y face, 0.6 mm
                                # below that face (at the edge's middle, z = BOX_OFF)
REF_FLOOR = 0.1e-3              # m: the reference's cells near that edge (0.0875 mm, the box mesh's level below it)
REF_STEPS = 100                 # the reference's backward-Euler steps over BOX_TIME
EDGE_RTOL = 0.05                # the probe column's temperature rise against the reference's over that column, relative
REGION_STEP = 0.02e-3           # m: the samples of a column's part of a level surface beside the edge (`_column_region`)
REGION_REACH = (4.99e-3, 4.0e-3)   # m: they reach this far along the edge from its middle and from the edge itself
EDGE_FAILS = ("the 2 mm markers do not resolve the heated edge (Task 7b's diagnosis D1, the controller's rulings T6-e and "
              "T7-j; re-measured in Task 7b, once each): the probe sits in a column that straddles the edge (an "
              "unheated-face column, its marker's normal +y, 0.50 to 1.98 mm from the heated face at the probe's 0.496 mm "
              "depth, 28 % of its part there under the heated face), whose one temperature per depth cannot follow the "
              "steep fall across the edge: it rises 26.30 K, 7.3 % under the reference's 28.37 K mean over that part "
              "(30.71 K at the probe itself); with the markers 1 mm apart the probe falls in a heated-face column, which "
              "rises 28.87 K against 30.73 K, 6.0 % under, still beyond EDGE_RTOL (5 %, kept): design 21.10's finer "
              "markers near edges are the remedy")


def _box_phi0():
    """The stage 3 box's phi0 (analysis/levelset/stage3_cases.py's `box_phi0`): the largest of its six planes'
    distances, an exact distance inside, from `cases.box_corners`."""
    from cases import box_corners
    lo, hi = box_corners()
    return lambda x: np.max(np.concatenate([lo - x, x - hi], axis=1), axis=1)


def _box_flux(P):
    """BOX_FLUX on the box's -x face as a beam along +x puts it [W/m2] on each patch P of phi: BOX_FLUX max(-n_x, 0), so
    that the heat in is BOX_FLUX times phi's shadow along x however the edges' cut cells blend their normals (the face's
    10 x 10 mm less what the P1 contour rounds off the shadow's corners: 1.6e-4 of it on the box's own mesh)."""
    return BOX_FLUX * np.maximum(-P.normal[:, 0], 0.0)


def _layer_temperature(body, p):
    """The body's temperature [K] at the point p (3,) between phi and beta: the column whose marker owns its foot on phi
    (p - phi(p) n, n phi's recovered normal; the nearest marker by the same-side rule, as the division divides the
    layer), at the depth -phi(p), linear between its cells' mid-depths. Returns (T, the marker, the depth [m])."""
    from scipy.spatial import cKDTree
    from reentry_model.levelset import markers as M
    mk, sk = body.markers, body.skins
    val, g = M._field(body.mesh, body.phi, np.asarray(p, float)[None])
    foot = p[None] - val[:, None] * g
    m = int(M._nearest(cKDTree(mk.x), foot, 1, mk.normal, g)[1][0, 0])
    row = mk.skin_index[m]
    dz = sk.dz()[row, :sk.count[row]]
    return float(np.interp(-val[0], np.cumsum(dz) - 0.5 * dz, sk.T()[row, :sk.count[row]])), m, float(-val[0])


def _column_region(body, m, depth):
    """Points sampling marker m's column's part of the level surface `depth` [m] below the box's surface beside the
    heated face's +y edge (the controller's ruling T6-b of plan 03b Task 6): the two planes of that level surface, under
    the heated face (x = lo + depth) and under the +y face (y = hi - depth), REGION_STEP apart, within REGION_REACH of
    the edge's middle along it and of the edge across it; each sample's foot on its face (along that face's normal)
    owned by the nearest marker the same-side rule allows, as the division divides the layer. Returns the column's
    samples (k, 3) [m] (equal areas)."""
    from scipy.spatial import cKDTree
    from cases import box_corners
    from reentry_model.levelset import markers as M
    lo, hi = box_corners()
    along = np.arange(-REGION_REACH[0], REGION_REACH[0], REGION_STEP) + 0.5 * (lo[2] + hi[2])
    across = np.arange(0.0, REGION_REACH[1], REGION_STEP)
    Z, U = (g.ravel() for g in np.meshgrid(along, across, indexing="ij"))
    heated = np.column_stack([np.full(Z.size, lo[0] + depth), hi[1] - depth - U, Z])
    unheated = np.column_stack([lo[0] + depth + U, np.full(Z.size, hi[1] - depth), Z])
    feet_h, feet_u = heated.copy(), unheated.copy()
    feet_h[:, 0], feet_u[:, 1] = lo[0], hi[1]
    points, feet = np.concatenate([heated, unheated]), np.concatenate([feet_h, feet_u])
    n = np.concatenate([np.tile([-1.0, 0.0, 0.0], (len(heated), 1)), np.tile([0.0, 1.0, 0.0], (len(unheated), 1))])
    owner = M._nearest(cKDTree(body.markers.x), feet, 1, body.markers.normal, n)[1][:, 0]
    return points[owner == m]


def _edge_reference(phi0):
    """CutThermal alone on the box (its phi0 exact on the nodes) at BOX_T0, BOX_FLUX on its -x face (`_box_flux`) for
    BOX_TIME in REF_STEPS steps, on the box's base mesh refined to REF_FLOOR within 1 mm of the heated face's +y edge, to
    0.175 mm within 2 mm of it, to 0.35 mm within 4 mm of the heated face and 0.7 mm elsewhere in the metal. Returns
    (the mesh, its nodal temperatures [K])."""
    from cases import box_corners
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.cutcell import patches
    from reentry_model.levelset.cutthermal import CutThermal
    lo, hi = box_corners()
    base = BoxBase.around((-2.8e-3, -5e-3 - 2.8e-3, -5e-3 - 2.8e-3), (8e-3 + 2.8e-3, 5e-3 + 2.8e-3, 5e-3 + 2.8e-3),
                          2.8e-3)

    def target(x, size):
        f = phi0(x)
        metal = f < 0.87 * size
        edge = np.hypot(x[:, 0] - lo[0], x[:, 1] - hi[1])
        want = np.where(edge < 1e-3 + 0.87 * size, REF_FLOOR, np.where(edge < 2e-3 + 0.87 * size, 0.175e-3,
                        np.where(x[:, 0] - lo[0] < 4e-3 + 0.87 * size, 0.35e-3, 0.7e-3)))
        return np.where(metal, want, np.inf)
    mesh = build(base, target)
    phi = phi0(mesh.points)
    solver = CutThermal(aa7075(), emissivity=0.0)
    solver.setup(mesh, phi)
    solver.set_temperature(BOX_T0)
    q = _box_flux(solver.P)
    for _ in range(REF_STEPS):
        solver.step(BOX_TIME / REF_STEPS, q)
    return mesh, solver.T


@pytest.mark.reference
@pytest.mark.xfail(strict=True, reason=EDGE_FAILS)
def test_heat_passes_a_box_edge_between_the_faces_columns():
    """Plan 03b Task 6 (design 21.5 and planning ruling 5 as the controller's ruling T6-b amends it: the skins conduct
    sideways across the division's shared boundaries, over the columns' centroids' distance): the stage 3 box (8 x 10 x
    10 mm, its default meshes, development skins) at 700 K, 1 MW/m2 on its -x face as a beam along +x puts it
    (`_box_flux`) for one 0.5 s step, nothing melting. At the probe 0.6 mm below the unheated +y face, 0.5 mm along it
    from the heated face's edge (EDGE_PROBE), the body's temperature rise -- its column's cell at the probe's depth
    (`_layer_temperature`), a mean over the column's part of the layer there -- is within EDGE_RTOL of the reference's
    mean over that same part (`_column_region`), CutThermal alone on the same box with 0.1 mm cells near that edge
    (`_edge_reference`); the books close. The reference's rise at the probe itself is printed too (Task 6 measured
    30.71 K there against 28.37 K over the column: the column's resolution, which this test does not judge)."""
    from cases import box_corners
    from reentry_model.levelset.cutcell import patches
    lo, hi = box_corners()
    phi0 = _box_phi0()
    probe = np.array([lo[0] + EDGE_PROBE[0], hi[1] - EDGE_PROBE[1], 0.5 * (lo[2] + hi[2])])
    started = time.perf_counter()
    base = BoxBase.around((-2.8e-3, -5e-3 - 2.8e-3, -5e-3 - 2.8e-3), (8e-3 + 2.8e-3, 5e-3 + 2.8e-3, 5e-3 + 2.8e-3),
                          2.8e-3)
    body = skinbody.LevelSetSkinBody(aa7075(), base, phi0, BOX_T0, SETTINGS, remove="feed", v_hat=(-1.0, 0.0, 0.0))
    built = time.perf_counter() - started
    P = patches(body.mesh, body.phi)
    W = body.markers.owned_patches(body.mesh, body.phi, P)
    q = np.asarray(W @ _box_flux(P)).ravel() / body.markers.area          # each marker's mean over what it owns
    heat = float(_box_flux(P) @ P.area) * BOX_TIME                          # the beam over phi's patches [J]
    started = time.perf_counter()
    book = body.advance(BOX_TIME, q_prescribed=q)
    stepped = time.perf_counter() - started
    T_body, m, depth = _layer_temperature(body, probe)
    started = time.perf_counter()
    mesh, T_ref = _edge_reference(phi0)
    referenced = time.perf_counter() - started
    region = _column_region(body, m, depth)
    at_probe = float(interpolate(mesh, T_ref, probe[None])[0]) - BOX_T0
    rise, want = T_body - BOX_T0, float(interpolate(mesh, T_ref, region).mean()) - BOX_T0
    print("box edge: built {:.0f} s, stepped {:.0f} s ({} parts, at most {} passes; heat in {:.4f} J against {:.4f} J; "
          "mismatch {:.2e} of it; residuals mass {:.1e} enthalpy {:.1e}); the probe in marker {}'s column (normal {}) at "
          "{:.3f} mm: rise {:.3f} K against the reference's {:.3f} K over the column's part ({} samples, {:.2f} to {:.2f} "
          "mm from the heated face; {:+.2%}) and {:.3f} K at the probe ({} cells, {:.0f} s)".format(
              built, stepped, book.parts, book.passes, book.enthalpy_in, heat,
              book.interface_mismatch_J / book.enthalpy_in, book.mass_residual, book.enthalpy_residual, m,
              np.round(body.markers.normal[m], 2), depth * 1e3, rise, want, len(region),
              (region[:, 0].min() - lo[0]) * 1e3, (region[:, 0].max() - lo[0]) * 1e3, rise / want - 1.0, at_probe,
              len(mesh.tets), referenced))
    assert np.nanmax(body.skins.T()) < body.material.T_solidus                     # nothing melting
    assert abs(book.mass_residual) <= 1e-12
    assert abs(book.enthalpy_residual) <= 1e-8
    assert book.enthalpy_in == pytest.approx(heat, rel=1e-9)
    assert abs(rise / want - 1.0) <= EDGE_RTOL


def _melt_all(body, T=1000.0):
    """The body's volume and skins all at T [K], above the liquidus."""
    body.solver.set_temperature(T)
    sk = body.skins
    sk.H = np.where(sk.valid, sk.m * sk.props.enthalpy(np.full(sk.m.shape, T)), 0.0)


def test_a_body_with_no_rigid_core_left_stops_with_a_clear_message():
    """The controller's ruling (d) of Task 10: a body whose metal is all above T_rigid has no rigid metal left, nothing for
    its molten metal to stand on; stage 3 holds the volume's melt (no flow solve), so a step that finds it so stops
    with a RuntimeError that says so, not some later failure of the geometry. The test cube, set to 1000 K throughout
    (volume and skins), one 0.01 s step with no heating."""
    body = cube_body()
    _melt_all(body)
    with pytest.raises(RuntimeError, match="no rigid metal left"):
        body.advance(0.01)


def test_a_failed_part_is_reported_as_no_rigid_metal_only_on_restored_state(monkeypatch):
    """Fix round 1 of Task 10 (the controller's ruling (4)): when a part fails, the step asks whether any rigid metal is
    left only once the part has restored the state it started from, and chains the original error. The test cube with
    the skins' conduction made to fail: with rigid metal left the original error comes through unchanged; all molten,
    the no-rigid-metal error comes, caused by the original; and if psi itself cannot be found then, the original comes
    through, not psi's error. The part restores the skins and the volume's temperatures it started from."""
    body = cube_body()

    def boom(self, *args, **kwargs):
        raise RuntimeError("boom: the skins' conduction failed")
    monkeypatch.setattr(Skins, "conduct", boom)
    m0, T0 = body.skins.m.copy(), body.solver.T.copy()
    with pytest.raises(RuntimeError, match="boom") as rigid:
        body.advance(0.01)
    assert rigid.value.__cause__ is None
    assert np.array_equal(body.skins.m, m0) and np.array_equal(body.solver.T, T0)
    _melt_all(body)
    with pytest.raises(RuntimeError, match="no rigid metal left") as molten:
        body.advance(0.01)
    assert "boom" in str(molten.value.__cause__)

    def no_psi(self, phi):
        raise ValueError("psi failed")
    monkeypatch.setattr(skinbody.LevelSetSkinBody, "_psi", no_psi)
    with pytest.raises(RuntimeError, match="boom"):
        body.advance(0.01)


def test_a_body_pickled_with_the_retired_mode_settings_loads_without_them():
    """Plan 03b Task 5: a body pickled before the modes were retired holds `mode_settings`; it loads without it, and one
    with no division or carry (from before the one skin depth) takes none and an empty carry, to be adopted by its next
    step (`_adopt`); so too one whose division has no base level's ownership (from before Task 5's fix round 1)."""
    import pickle
    from dataclasses import replace
    body = cube_body()
    body.__dict__["mode_settings"] = {"hold_thick": 1}
    back = pickle.loads(pickle.dumps(body))
    assert "mode_settings" not in back.__dict__ and back.division is not None
    division = body.division
    body.division = replace(division, base_cell=None, base_patch_area=None, base_owner=None)
    assert pickle.loads(pickle.dumps(body)).division is None
    body.division = division
    del body.__dict__["division"], body.__dict__["carry"], body.__dict__["mode_settings"]
    old = pickle.loads(pickle.dumps(body))
    assert old.division is None and np.array_equal(old.carry, np.zeros(old.skins.n))


def test_a_step_book_pickled_with_the_old_name_of_its_repartitioned_share_loads_under_the_new_one():
    """Plan 03b Task 7 (Task 6's minor item): since Task 6 the step end's re-partitioned share is a share of the layer's
    volume, so the book's field and the history's column are `volume_repartitioned_share`; a book pickled with the old
    name `area_repartitioned_share` (a body's `last_step`, kept by the cases' states) loads with its value under the new
    one, and a book pickled now loads unchanged."""
    import pickle
    from dataclasses import MISSING, fields
    book = skinbody.StepBook(**{f.name: 0 for f in fields(skinbody.StepBook)
                                if f.default is MISSING and f.default_factory is MISSING})
    book.volume_repartitioned_share = 0.25
    assert book.row()["volume_repartitioned_share"] == 0.25 and "area_repartitioned_share" not in book.row()
    assert pickle.loads(pickle.dumps(book)) == book
    state = dict(book.__dict__)
    state["area_repartitioned_share"] = state.pop("volume_repartitioned_share")
    old = skinbody.StepBook.__new__(skinbody.StepBook)
    old.__setstate__(state)
    assert old == book and "area_repartitioned_share" not in old.__dict__
