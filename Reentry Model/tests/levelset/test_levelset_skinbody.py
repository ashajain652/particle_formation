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

from cases import CAPPED_SHARE, gas_state, nearest_ratio

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


def cube_body(**kw):
    base = BoxBase.around((-9e-3,) * 3, (9e-3,) * 3, 2.8e-3)
    args = dict(MESH, emissivity=0.0)
    args.update(kw)
    return skinbody.LevelSetSkinBody(aa7075(), base, cube, 850.0, SETTINGS, **args)


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


def test_a_column_with_no_base_takes_its_shortfall_from_its_neighbours(monkeypatch):
    """Design 21.5: a column that ends before beta, beside a convex edge, has no base: it never draws or returns, and is
    filled from the neighbour across the boundary it shares deepest. On the test cube an edge column made baseless (its
    division's cross-section at H given to that neighbour's, so that the columns still tile the layer), 0.05 mm of its
    top taken (as a release takes it) and the layer accepted (phi moved by it, the division, the fill): the fill never exchanges it at a base, it is full afterwards, the neighbour it
    took from (`_neighbours`: one with a base, sharing its deepest level) is filled at its own base, and the books -- the skins' and the volume's mass with what was taken, the enthalpy with the
    heat owed -- are exact."""
    from dataclasses import replace
    from reentry_model.levelset import stepend
    body = cube_body()
    sk, mk = body.skins, body.markers
    m = _edge_column(body)
    row = mk.skin_index[m]
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
    body._divide(skinbody._Tally(sk.n))
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


def gas_body():
    """The cube with the shared heating and surface flow (`loads.GasLoads`) and Girin's spray model, radiating."""
    from reentry_model.levelset import loads
    from reentry_model.shared import heating, spray, surface_flow
    from reentry_model.shared.material import LiquidProperties
    mat = aa7075()
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), mat.liquid.rho)
    model = spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), spray_gates="2026-11-10")
    return cube_body(loads=gl, spray_model=model, T_ambient=200.0, emissivity=0.2)


PART_BOUND = 100      # parts a flow-on step of the test cube may take at most: a tenth of PARTS_MAX
CONTINUUM = dict(rho=1e-4, kn=1e-3)   # 7.5 km/s at a density of 1e-4 kg/m3 with a body Knudsen number of 1e-3: the shock
                                      # layer's branch, where Girin's closure gives a conjugate depth


FOLDED = ("plan 03b Task 5: at this state the runoff carries the +x face's melt, every skin's 1.2 mm of it now flowing "
          "(the retired top-zone skins held their deep melt in the volume), so far within the 0.05 s step that the face "
          "stands from 1.0 mm below to 1.8 mm above where it began (measured once), its folds' columns face each other, "
          "and the step end stops on stepend.THIN_METAL (19 pairs); for the controller, with Task 6's runoff on the "
          "division's faces")


@pytest.fixture(scope="module")
def sprayed():
    """The cube with liquid formed on its +x face by 0.25 s of 6 MW/m2 (no flow, nothing removed), then one 0.05 s macro
    step at the continuum gas state with the runoff and Girin's plan on: (body, the melting step's book, the spraying
    step's book or None, its time [s], the RuntimeError that stopped it or None; FOLDED)."""
    body = gas_body()
    melt = body.advance(0.25, q_prescribed=one_face(6e6))
    started = time.perf_counter()
    try:
        book, stop = body.advance(0.05, state=gas_state(**CONTINUUM)), None
    except RuntimeError as err:
        book, stop = None, err
    return body, melt, book, time.perf_counter() - started, stop


def _unfolded(stop):
    """xfail the test where the spraying step stopped on THIN_METAL (FOLDED), after checking that it did."""
    if stop is not None:
        assert str(stop).startswith(skinbody.THIN_METAL), stop
        pytest.xfail(FOLDED)


def test_liquid_formed_on_the_face_sprays_and_the_books_close(sprayed):
    """Fix round 1 of Task 9 (ruling b): liquid formed on the +x face sprays in the next 0.05 s macro step at the
    continuum gas state: the released mass is positive and equals the plan's books and the history's, the mass and
    enthalpy close, and the heat still owed and the interface mismatch stay within OWED_SHARE of the heat in."""
    body, melt, book, seconds, stop = sprayed
    _unfolded(stop)
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


def test_two_flow_on_steps_keep_their_part_count_bounded(sprayed):
    """Fix round 1 of Task 9 (ruling a), the smallest case found that reproduced the skyscrapers' runaway: the spraying
    step above and a second 0.05 s step at the same state. A part that meets its allowance only after splits, or never,
    carries the splits' shortness into the skyscrapers' estimate; here a step's hint is at most 1.25 times the parts it
    took: both steps stay within PART_BOUND parts, the second step's hint is at most the skyscrapers' for it, both hints
    are 1.25 times their step's parts at most, and the books close."""
    body, _, first, _, stop = sprayed
    _unfolded(stop)
    assert first.parts_next <= max(1, int(np.ceil(skinbody.MOVE_HINT_MARGIN * first.parts)))
    second = body.advance(0.05, state=gas_state(**CONTINUUM))
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


SHARP_FACTOR = 2.5    # the tangent-plane cells' areas over the area nearest each marker on the sharp test cube, within
                      # this either way (plan 03b Task 5: 0.47 to 1.23 measured once, against OWNERSHIP_FACTOR's 1.5 when
                      # the stub rounded the cube's edges; the cells serve only runoff's faces until Task 6 retires them)


def test_the_test_cube_s_markers_keep_their_neighbours():
    """Task 9b (the ruling on the neighbour test): on the test cube as the body builds it (its edges sharp since plan 03b
    Task 5; markers one per 2 mm voxel) at most CAPPED_SHARE of the cells are bounded by the cap and every tangent-plane
    cell's area is within a factor of SHARP_FACTOR of the area nearest the marker."""
    body = cube_body()
    mk = body.markers
    mk.tessellate(body.mesh, body.phi)
    r = nearest_ratio(mk, body.mesh, body.phi)
    print("test cube: {} markers, {} capped; areas {:.3f} to {:.3f} times the nearest-marker areas".format(
        mk.n, int(mk.capped.sum()), r.min(), r.max()))
    assert mk.capped.mean() < CAPPED_SHARE
    assert r.min() >= 1.0 / SHARP_FACTOR and r.max() <= SHARP_FACTOR


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


THROUGH_SHARE = 1e-2  # of the inner columns' volume: what the fill after a re-partition by mass moves through the base


def test_a_repartition_closes_the_skins_books_and_keeps_their_depths():
    """Task 9b: liquid skins on a closed flat box, the 16 markers in the middle of its top face shifted along it by up to
    0.4 mm each, the surface unchanged: the step end's tessellation re-partitions the area (the shifted markers' areas
    change by 0.82 to 1.16 times, the area of phi each owns) and each column's metal follows its area by the hand-over
    rule. The skins' mass and enthalpy close to 1e-14 and their momentum to round-off; every new column's flat depth lies
    within its donors'; the stretch is one, since the surface did not change; and on the top face's middle, where the
    columns are flat, the fill that follows moves through the base only what the transfer's sampled areas miss of the
    division's (THROUGH_SHARE: 2.3e-3 of those columns' volume, measured once; the re-partition by depth across the
    division's faces, plan 03b Task 6, is to remove it)."""
    from materials import LinearMeltMaterial
    mat = LinearMeltMaterial()
    host, top = _slab_host(mat)
    mk, sk = host.markers, host.skins
    middle = top[1:5, 1:5].ravel()
    rng = np.random.default_rng(7)
    rng.uniform(1.5, 3.0, mk.n)                       # (the draws Task 9b's skins of one to three stubs took: the
    rng.random(len(middle)), rng.uniform(1.0, 3.0, len(middle))      # same shifts as its measurement below)
    sk.fill(sk.density(np.full(sk.n, 950.0)) * sk.column_volume(), 950.0)    # liquid: above the liquidus, 908 K
    sk.U = np.where(sk.valid[..., None], np.concatenate([rng.normal(0.0, 1.0, (sk.n, 1, 2)),
                                                         np.zeros((sk.n, 1, 1))], axis=2), 0.0)
    M0, E0, P0 = (np.sum(t, axis=0) for t in sk.totals())
    P_size = np.abs(sk.m[..., None] * sk.U).sum()
    thick0, A0, n0 = _flat_depth(sk)[mk.skin_index], mk.area.copy(), mk.n
    mk.x[middle, :2] += rng.uniform(-0.4e-3, 0.4e-3, (len(middle), 2))
    host._settle(skinbody._Tally(0))
    assert mk.n == n0
    M1, E1, P1 = (np.sum(t, axis=0) for t in sk.totals())
    ratio = mk.area[middle] / A0[middle]
    inner = mk.skin_index[top[2:4, 2:4].ravel()]
    miss = _would_fill(host)
    through = float(np.abs(miss[inner]).sum() / sk.column_volume()[inner].sum())
    print("re-partition: shifted areas {:.3f} to {:.3f} of their old; books mass {:.1e} enthalpy {:.1e} momentum {:.1e}; "
          "the fill would move {:.1e} of the inner columns' volume through the base".format(
              ratio.min(), ratio.max(), M1 / M0 - 1, E1 / E0 - 1, np.abs(P1 - P0).max() / P_size, through))
    assert ratio.min() < 0.85 and ratio.max() > 1.1                    # the shift does re-partition the area
    assert through <= THROUGH_SHARE
    assert abs(M1 / M0 - 1) <= BOOKS_RTOL and abs(E1 / E0 - 1) <= BOOKS_RTOL
    assert np.abs(P1 - P0).max() <= BOOKS_RTOL * P_size
    T, stretch = host.last_settle[0].transfer, host.last_settle[0].stretch
    assert np.abs(stretch - 1.0).max() <= 1e-14
    thick1 = _flat_depth(sk)[mk.skin_index]
    for i in range(mk.n):
        donors = T[:, i].nonzero()[0]
        assert thick0[donors].min() * (1 - 1e-12) <= thick1[i] <= thick0[donors].max() * (1 + 1e-12)


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
