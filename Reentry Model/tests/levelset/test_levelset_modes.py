"""Tests of the level set's two skin modes (design sections 20.5-20.6; level-set stage 3 Task 10): the thin test per
marker with its hold and its group (`modes.update`, stage 1's `thin.ModeState` keyed by marker), the exact mode changes
(`modes.to_thick`, `modes.to_thin`), a top-zone skin's exchange with the volume at its base (`modes.hand_down`,
`modes.minimum_mass`, `modes.draw_cap`) and Ruling 3, the exchange between a full-depth skin and a top-zone neighbour
through the top 0.4 mm only. The columns are built directly, cell by cell, on markers whose faces are given (no mesh)."""
import numpy as np
import pytest

from reentry_model.levelset import modes, runoff
from reentry_model.levelset.markers import Faces, Markers
from reentry_model.levelset.skin import Skins, SkinSettings
from reentry_model.shared.material import Material

LIQUID, RIGID = 950.0, 860.0                 # K: all liquid; rigid (AA7075's T_rigid 895.1 K)
CELL = 1.0e-5                                # m: each stacked cell holds 10 um of solid metal's mass
SETTINGS = SkinSettings(cell=10e-6, fine_zone=0.4e-3, stub=0.8e-3)
BOOKS_RTOL = 1e-14                           # the brief's tolerance on a mode change's mass and enthalpy


def aa7075():
    return Material.from_drama_json("AA7075_scheil")


def grid_markers(nx, ny, h=1e-3):
    """Markers on a flat regular grid of nx by ny points h apart in z = 0, normals +z, square cells of area h^2 and a
    face of length h to each grid neighbour (built with no mesh), all full-depth with no hold counted."""
    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    x = np.column_stack([ix.ravel() * h, iy.ravel() * h, np.zeros(ix.size)])
    n = len(x)
    idx = np.arange(n).reshape(nx, ny)
    i = np.r_[idx[:-1, :].ravel(), idx[:, :-1].ravel()]
    j = np.r_[idx[1:, :].ravel(), idx[:, 1:].ravel()]
    unit = (x[j] - x[i]) / h
    mk = Markers(x, np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, h * h), np.zeros((n, 2)), np.arange(n),
                 np.full(n, -1), h)
    mk.faces = Faces(i, j, np.full(len(i), h), unit)
    return mk


def stacked(area, temps, mat=None, rng=None, settings=SETTINGS):
    """Skins of the given top areas [m2], column i holding a cell of 10 um of solid metal's mass at each temperature of
    temps[i], top first; every cell moving with a random velocity (a part along the +z normal included) given `rng`."""
    mat = aa7075() if mat is None else mat
    area = np.asarray(area, dtype=float)
    sk = Skins(mat, area, np.zeros((len(area), 2)), settings)
    count = np.array([len(t) for t in temps])
    K = int(count.max())
    T = np.array([list(t) + [300.0] * (K - len(t)) for t in temps])
    valid = np.arange(K)[None, :] < count[:, None]
    m = np.where(valid, mat.rho * area[:, None] * CELL, 0.0)
    U = None if rng is None else np.where(valid[..., None], rng.normal(0.0, 0.3, (len(area), K, 3)), 0.0)
    sk.install(m, m * sk.props.enthalpy(T), count, U)
    return sk


def _split(P, normal):
    """Momentum (k, 3) split into its parts along the surface and along the normal (k, 3)."""
    along = np.einsum("ij,ij->i", P, normal)[:, None] * normal
    return P - along, along


# -- the mode changes are exact ------------------------------------------------------------------------------------------
def test_mode_changes_conserve_mass_enthalpy_and_tangential_momentum():
    """A column of melt over rigid metal, both ways (design 20.6). Thin to thick (`to_thick`): the full-depth skin's cells
    below its fine zone -- the deeper melt and the rigid stub -- go to the volume, what stays keeping at least a top-zone
    skin's minimum mass; mass and enthalpy close to 1e-14, the momentum along the surface exactly (the skin's plus the
    hand-over's) and the normal part, which no skin carries, booked apart; beta steps out by the mass over the volume's
    density. Thick to thin (`to_thin`): the volume's metal beneath the top-zone skin down to a depth -- 1 mm of melt then
    the 0.8 mm stub of rigid metal below psi -- joins the skin, cells filled from the volume's profile along the normal,
    its mass the volume's density times the depth's volume exactly, at rest; beta steps in by it. Two columns, the second
    tilted against its normal so that the normal part is not trivial."""
    mat = aa7075()
    rng = np.random.default_rng(11)
    temps = [[LIQUID] * 150 + [RIGID] * 80, [LIQUID] * 90 + [RIGID] * 100]
    sk = stacked([1e-6, 2e-6], temps, mat, rng)
    normal = np.array([[0.0, 0.0, 1.0], [0.6, 0.0, 0.8]])
    rows = np.array([0, 1])
    M0, E0, P0 = sk.totals()
    top0 = sk.m.copy()
    ho = modes.to_thick(sk, rows, normal=normal, rho=mat.rho)
    M1, E1, P1 = sk.totals()
    print("thin to thick: handed {} kg, kept {} m of skin".format(ho.mass, sk.thickness()))
    np.testing.assert_allclose(M1 + ho.mass, M0, rtol=BOOKS_RTOL, atol=0.0)
    np.testing.assert_allclose(E1 + ho.enthalpy, E0, rtol=BOOKS_RTOL, atol=0.0)
    t0, n0 = _split(P0, normal)
    t1, n1 = _split(P1, normal)
    size = np.abs(sk.m[..., None] * sk.U).sum() + np.abs(P0).sum()
    np.testing.assert_allclose(t1 + ho.momentum, t0, rtol=0.0, atol=1e-15 * size)
    np.testing.assert_allclose(n1 + ho.normal, n0, rtol=0.0, atol=1e-15 * size)
    assert np.abs(ho.normal).max() > 1e-3 * np.abs(ho.momentum).max()          # the normal part is booked, not lost
    np.testing.assert_allclose(np.einsum("ij,ij->i", ho.momentum, normal), 0.0, atol=1e-15 * size)
    np.testing.assert_allclose(ho.volume, ho.mass / mat.rho, rtol=1e-15)
    # the top zone stays as it was, cell for cell; at least a top-zone skin's minimum mass, no more than a cell beyond
    for r in rows:
        c = sk.count[r]
        assert np.array_equal(sk.m[r, :c], top0[r, :c])
    m_keep = mat.rho * sk.area * SETTINGS.fine_zone
    assert np.all(sk.m.sum(axis=1) >= m_keep * (1 - 1e-14))
    assert np.all(sk.m.sum(axis=1) - sk.m[rows, sk.count - 1] <= m_keep * (1 + 1e-12))
    assert np.all(sk.thickness() < SETTINGS.fine_zone * mat.rho / mat.liquid.rho + 2 * CELL)

    # thick to thin: 1 mm of melt under the top-zone skin, then rigid metal, absorbed down to the stub below psi
    melt = 1.0e-3
    depth = np.full(2, melt + SETTINGS.stub)
    profile = lambda k, s: np.where(s < melt, LIQUID, RIGID)
    M1, E1, P1 = sk.totals()
    H1 = sk.thickness()
    hi = modes.to_thin(sk, rows, depth, profile, mat.rho)
    M2, E2, P2 = sk.totals()
    print("thick to thin: took {} kg from the volume, skin now {} m".format(-hi.mass, sk.thickness()))
    np.testing.assert_allclose(-hi.mass, mat.rho * sk.area * depth, rtol=BOOKS_RTOL)
    np.testing.assert_allclose(M2, M1 - hi.mass, rtol=BOOKS_RTOL, atol=0.0)
    np.testing.assert_allclose(E2, E1 - hi.enthalpy, rtol=BOOKS_RTOL, atol=0.0)
    np.testing.assert_allclose(hi.volume, hi.mass / mat.rho, rtol=1e-15)
    assert np.all(hi.momentum == 0.0) and np.all(hi.normal == 0.0)
    np.testing.assert_allclose(P2, P1, rtol=0.0, atol=1e-15 * size)
    # the volume's profile along the normal: melt then the stub of rigid metal, anchored at the bottom
    T = sk.T()
    anchored = sk.anchored()
    dz = sk.dz(T)
    run = np.where(anchored, dz, 0.0).sum(axis=1)
    np.testing.assert_allclose(run, SETTINGS.stub, rtol=0.02)
    assert np.all(sk.thickness() > H1 + melt)


# -- the thin test per marker: the hold and the group -----------------------------------------------------------------
def test_a_new_marker_may_change_mode_only_after_the_hold():
    """Design 20.5: a marker changes mode only after its thin test has given the other answer at `hold` step ends in a
    row. With a hold of three both ways: new markers (as the body creates them, all full-depth with no hold counted)
    whose test wants thick from the first step end do not change at the first nor the second update, and change at the
    third. A marker split between the second and third updates takes its parent's mode and hold, so its new half changes
    with the group (alone it could never make the group a change needs). A marker whose test wavers starts its count
    again."""
    st = modes.ModeSettings(hold_thick=3, hold_thin=3)
    mk = grid_markers(4, 4)
    want = np.zeros(mk.n, bool)                         # thin_candidate's answer: False, thick, everywhere
    want[15] = True                                     # one corner wavers: thick, thick, thin, thick ...
    assert not modes.update(mk, want, st).any() and not mk.top.any()
    assert np.array_equal(mk.hold, np.where(want, 0, 1))
    want[15] = False
    assert not modes.update(mk, want, st).any()
    split = mk.split(0.4 * mk.area[5])                  # every marker is above twice that: all split, halves appended
    assert len(split) == 16 and mk.n == 32
    assert np.array_equal(mk.hold[16:], mk.hold[:16]) and not mk.top.any()
    i, j = np.arange(16), 16 + np.arange(16)            # each half a neighbour of its parent
    old = grid_markers(4, 4).faces
    mk.faces = Faces(np.r_[old.i, i], np.r_[old.j, j], np.r_[old.length, np.full(16, 1e-3)],
                     np.r_[old.normal, np.tile([1.0, 0.0, 0.0], (16, 1))])
    want = np.zeros(mk.n, bool)
    changed = modes.update(mk, want, st)
    late = np.array([15, 31])                           # the wavering corner and its half: counted twice, not three times
    assert changed.sum() == 30 and not changed[late].any()
    assert mk.top[~np.isin(np.arange(32), late)].all() and not mk.top[late].any()
    assert np.all(mk.hold[~np.isin(np.arange(32), late)] == 0) and np.all(mk.hold[late] == 2)
    # the two late ones are a group of two: under the seven a change needs, they stay full-depth whatever the test says
    for _ in range(4):
        assert not modes.update(mk, want, st).any()
    assert not mk.top[late].any()


def test_the_hold_is_one_step_end_to_thick_and_three_back_to_thin():
    """The controller's ruling (2) of Task 10's fix round 1 (a design value for Asha): the defaults change a full-depth
    skin to top-zone after one step end of the test saying thick (`hold_thick`), since deep pooled melt in a full-depth
    skin is the failure path (the full-size cube failed in its second step with a hold of three, before any change), and
    back to full-depth only after three step ends saying thin (`hold_thin`), against flicker; the group stays seven. New
    markers wanting thick change at the first update; wanting thin again, they change back at the third, not before; one
    step end of thick in between starts the count again."""
    st = modes.ModeSettings()
    assert (st.c_kappa, st.c_slope, st.hold_thick, st.hold_thin, st.min_markers) == (0.1, 0.2, 1, 3, 7)
    mk = grid_markers(4, 4)
    want = np.zeros(mk.n, bool)
    assert modes.update(mk, want, st).all() and mk.top.all() and np.all(mk.hold == 0)
    want[:] = True
    assert not modes.update(mk, want, st).any() and not modes.update(mk, want, st).any()
    want[:] = False                                     # the test wavers: thick once
    assert not modes.update(mk, want, st).any() and mk.top.all() and np.all(mk.hold == 0)
    want[:] = True
    for k in range(3):
        changed = modes.update(mk, want, st)
        assert changed.all() == (k == 2) and (not changed.any()) == (k < 2)
    assert not mk.top.any()


def test_a_small_group_does_not_change_mode():
    """Design 20.5: a mode change needs a connected group of at least `min_markers` (7) markers whose holds are full.
    On a row of 20 markers, six neighbours wanting thick stay thin however long their test says so; seven change at the
    first update (the default hold to thick); a lone marker never does."""
    st = modes.ModeSettings()
    mk = grid_markers(20, 1)
    want = np.ones(mk.n, bool)
    want[2:8] = False                                   # six in a row want thick
    want[12:19] = False                                 # seven in a row want thick
    for k in range(6):
        changed = modes.update(mk, want, st)
        if k == 0:
            assert np.array_equal(np.flatnonzero(changed), np.arange(12, 19))
        else:
            assert not changed.any()
    assert not mk.top[2:8].any() and mk.top[12:19].all() and mk.top.sum() == 7
    # and back: the seven, wanting thin again, change back together after the hold to thin (three)
    want[12:19] = True
    for k in range(3):
        changed = modes.update(mk, want, st)
        assert changed.any() == (k == 2)
    assert np.array_equal(np.flatnonzero(changed), np.arange(12, 19)) and not mk.top.any()


# -- a top-zone skin's base and Ruling 3 ------------------------------------------------------------------------------
def test_a_top_zone_skin_keeps_its_fine_zone_against_the_volume():
    """Design 20.6: a top-zone skin keeps the top 0.4 mm (its minimum mass a solid fine zone's, against a full-depth
    skin's stub), returning every cell below it to the volume at the base whatever its phase, with its momentum; a
    full-depth skin returns only the anchored rigid metal below its stub. Drawn short of its minimum, a top-zone skin
    takes the volume's metal at the volume's own temperature (thick melt beneath it), a full-depth skin at most
    T_rigid (rigid metal beneath its stub)."""
    mat = aa7075()
    rng = np.random.default_rng(5)
    temps = [[LIQUID] * 100 + [RIGID] * 120, [LIQUID] * 100 + [RIGID] * 120]
    sk = stacked([1e-6, 1e-6], temps, mat, rng)
    top = np.array([False, True])
    m_min = modes.minimum_mass(sk, top)
    np.testing.assert_allclose(m_min, mat.rho * sk.base_area() * np.array([SETTINGS.stub, SETTINGS.fine_zone]),
                               rtol=1e-15)
    np.testing.assert_array_equal(modes.draw_cap(sk, top), [mat.T_rigid, np.inf])
    M0, E0, P0 = sk.totals()
    z_top = np.cumsum(sk.dz(), axis=1) - sk.dz()
    below_stub = sk.valid[0] & (z_top[0] >= z_top[0, 100] + SETTINGS.stub * (1 - 1e-12))
    mass, H, P = modes.hand_down(sk, top, SETTINGS.stub)
    M1, E1, P1 = sk.totals()
    np.testing.assert_allclose(M1 + mass, M0, rtol=BOOKS_RTOL)
    np.testing.assert_allclose(E1 + H, E0, rtol=BOOKS_RTOL)
    np.testing.assert_allclose(P1 + P, P0, rtol=0.0, atol=1e-15 * np.abs(P0).sum())
    # the full-depth column: its rigid run (860 K, 17 % liquid, a little swollen) is over 1.2 mm deep, and its cells
    # below the stub's depth beneath the run's top go, with whatever momentum they held
    assert 40 <= below_stub.sum() < 50
    assert mass[0] == pytest.approx(below_stub.sum() * mat.rho * 1e-6 * CELL, rel=1e-12)
    # the top-zone column: everything below the fine zone, melt and rigid alike, keeping the minimum mass
    assert sk.m[1].sum() >= m_min[1] * (1 - 1e-14)
    assert sk.thickness()[1] < SETTINGS.fine_zone * mat.rho / mat.liquid.rho + 2 * CELL     # its melt swells
    assert mass[1] > 150 * mat.rho * 1e-6 * CELL and np.abs(P[1]).max() > 0.0
    # draws: a top-zone skin emptied to 0.2 mm takes 0.2 mm of the volume's melt at 950 K, a full-depth one is capped
    sk2 = stacked([1e-6, 1e-6], [[LIQUID] * 20, [LIQUID] * 20], mat)
    drawn, dH = sk2.draw(LIQUID, 0.0, modes.minimum_mass(sk2, np.array([False, True])),
                         cap=modes.draw_cap(sk2, np.array([False, True])))
    h = dH / drawn
    assert h[1] == pytest.approx(float(sk2.props.enthalpy(np.array([LIQUID]))[0]), rel=1e-12)
    assert h[0] == pytest.approx(float(sk2.props.enthalpy(np.array([mat.T_rigid]))[0]), rel=1e-12)


def _row_of_three(top):
    """Three markers 1 mm apart along x: A (1.5 mm of liquid over 0.4 mm of rigid metal), B (0.4 mm of liquid: a top-zone
    skin over held melt) and C (as A), every flowing cell moving at 0.05 m/s along +x; `top` their modes."""
    mk = grid_markers(3, 1)
    mk.top = np.asarray(top, bool)
    temps = [[LIQUID] * 150 + [RIGID] * 40, [LIQUID] * 40, [LIQUID] * 150 + [RIGID] * 40]
    sk = stacked(mk.area, temps, settings=SkinSettings(cell=10e-6, fine_zone=0.4e-3, stub=0.4e-3))
    flowing = sk.flowing()
    sk.U = np.where(flowing[..., None], np.array([0.05, 0.0, 0.0]), 0.0)
    mk.skins = sk
    return mk, sk


def test_thin_and_thick_neighbours_exchange_only_their_top_zone():
    """Ruling 3 (flow off at the thin-thick boundary): a full-depth skin and a neighbouring top-zone skin exchange only
    through the top 0.4 mm, layer by layer; the full-depth skin's deeper cells neither run into the held volume under the
    top-zone skin nor receive from it. On the row A (full) - B (top) - C (full), one transport step: A's cells below its
    fine zone keep their metal, C's cells below its fine zone receive none, the books close; with all three full-depth,
    A's deeper cells do run on into B. Across a mode boundary the thickness term sees no step in the heights (the floor
    of a top-zone skin is held melt, not the surface)."""
    dt = 2e-3
    out = {}
    for top in ((False, True, False), (False, False, False)):
        mk, sk = _row_of_three(top)
        sheet = runoff.Sheet.of(mk)
        m0, H0 = sk.m.copy(), sk.H.copy()
        z_top = np.cumsum(sk.dz(), axis=1) - sk.dz()
        deep = sk.valid & (z_top >= 0.4e-3 * (1 - 1e-9))
        book = runoff.transport(sk, sheet, dt)
        assert abs(sk.m.sum() / m0.sum() - 1.0) <= 1e-14 and abs(sk.H.sum() / H0.sum() - 1.0) <= 1e-14
        out[top] = (m0, sk.m.copy(), deep, book)
    m0, m1, deep, book = out[(False, True, False)]
    print("mixed row: A gave {:.3e} kg, B net {:.3e} kg, C gained {:.3e} kg".format(-book.mass[0], book.mass[1],
                                                                                book.mass[2]))
    assert book.mass[0] < 0.0 and book.mass[2] > 0.0
    np.testing.assert_array_equal(m1[0][deep[0]], m0[0][deep[0]])       # A's deep cells kept their metal
    K = min(m0.shape[1], m1.shape[1])
    shallow_C = int((~deep[2] & (np.arange(m0.shape[1]) < 190)).sum())   # C received only into its top zone (no new top)
    assert m1.shape[1] >= 190
    np.testing.assert_array_equal(m1[2, shallow_C:190], m0[2, shallow_C:190])
    m0f, m1f, deepf, _ = out[(False, False, False)]
    assert np.any(m1f[0][deepf[0]] < m0f[0][deepf[0]])               # all full-depth: the deep cells do run on
    # the thickness term across the mode boundary: none at A (its only neighbour is B, a top-zone skin)
    mk, sk = _row_of_three((False, True, False))
    forcing = runoff.Forcing(np.zeros(3), np.zeros(3), np.array([0.0, 0.0, -30.0]), np.array([0.0, 0.0, 1.0]), 0.8,
                             np.zeros(3))
    _, terms = runoff._push(runoff.Sheet.of(mk), sk, forcing.dp_ds, forcing.g, forcing.v_hat, forcing.sigma,
                            np.zeros(3), None)
    assert np.all(terms["thickness"][[0, 2]] == 0.0)
    mk, sk = _row_of_three((False, False, False))
    _, terms = runoff._push(runoff.Sheet.of(mk), sk, forcing.dp_ds, forcing.g, forcing.v_hat, forcing.sigma,
                            np.zeros(3), None)
    assert np.abs(terms["thickness"][0]).max() > 0.0
