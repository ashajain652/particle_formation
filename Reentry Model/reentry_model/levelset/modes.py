"""The two modes of the level set's skins (design sections 20.5-20.6; level-set stage 3 Task 10): each marker's skin is
full-depth (thin mode: from the gas down through the molten layer to the 0.8 mm stub of rigid metal below psi) or a
top-zone skin (thick mode: the top 0.4 mm, the fine zone, riding on thick melt that the volume mesh holds), and the
changes between the two are exact.

The thin test per marker (`measures`, `candidates`, `update`). Stage 1's test (`thin.patch_measures`, `thin_candidate`)
read at the markers: the molten depth h from phi to psi, the largest principal curvature of phi at the marker and of
psi at its foot, and the slope of the depth along the surface; thin when h is at most c_kappa over the curvature and the
slope at most c_slope. Stage 1's hysteresis (`thin.ModeState`) keyed by marker: a marker changes mode only after the
test has given the other answer at a number of step ends in a row -- `hold_thick` (1) to turn top-zone, `hold_thin` (3)
to turn back full-depth (the controller's ruling (2) of Task 10's fix round 1, a design value for Asha: deep pooled melt
in a full-depth skin is the failure path, and the full-size cube failed in its second step before a hold of three let
any marker change; flicker back to thin is what the hold guards against) -- and only with a connected group (on the
markers' faces) of at least `min_markers` markers ready to change. The mode and the count are the markers' own
(`Markers.top`, `Markers.hold`), so they survive split (a new half takes its parent's), merge (a top-zone marker's, when
the two differ) and release; markers start full-depth with no hold counted.

The mode changes, exact in mass and enthalpy (the volume's side through the body's pending load), each a `Handover` of
what passed to the volume (negative: from it):
- thin to thick (`to_thick`): the full-depth skin's cells lying wholly below its fine zone go to the volume, from the
  bottom up while the skin keeps at least a solid fine zone's mass (as its base exchange then keeps it, so that the next
  draw takes nothing back); beta steps out by their mass over the volume's density, the nodes it brings inside the metal
  taking the handed metal's temperature (`stepend.StepEnd._fill`); their momentum goes to the held volume (`momentum_to_volume`, Ruling 4), its part along the surface exactly and its normal part, which no skin
  carries, booked apart.
- thick to thin (`to_thin`): the volume's metal beneath the marker, from the skin's base down to the stub below psi,
  joins the skin as new cells at rest (the volume flow is held), graded as the skins' cells are, each the volume's
  density times its span's volume under the column's cross-section, at the temperature the volume's field gives at its
  middle along the normal; beta steps in by it.
- A top-zone skin keeps its fine zone against the volume every part (`hand_down`, `minimum_mass`, `draw_cap`): every
  cell below the fine zone returns at the base whatever its phase, and a skin short of a solid fine zone's mass draws the
  volume's metal at the volume's own temperature (thick melt; a full-depth skin draws rigid metal, capped at T_rigid).
  Its base velocity is zero (the skins' no-slip floor; the volume flow is held, Ruling 3), and the momentum its base
  takes goes to the volume.
- The step end's re-partition across modes (the controller's ruling (a) of 2026-11-12): a top-zone skin that received
  metal below its fine zone hands it to the volume by `to_thick`'s rule, beta stepping out; a full-depth skin that
  received top-zone area takes the thinner column (`stepend.StepEnd._modes`).

Every number here is a numerical setting (`ModeSettings`, the thin test's coefficients, design 20.5's starting values);
the physics is the shared package's."""
from collections import namedtuple
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix

from .criteria import BIG, finite_values
from .skin import volume_per_area
from .thin import ModeState, patch_measures, thin_candidate
from .transfer import locate


@dataclass(frozen=True)
class ModeSettings:
    """The thin test's settings (design 20.5, the builder's starting values, recorded with every run): `c_kappa`
    (`--thin-curvature`), `c_slope` (`--thin-slope`), the holds (`--mode-hold`, step ends in a row) `hold_thick` to turn
    top-zone and `hold_thin` to turn back full-depth (module docstring), and `min_markers` (`--mode-min-markers`)."""
    c_kappa: float = 0.1
    c_slope: float = 0.2
    hold_thick: int = 1
    hold_thin: int = 3
    min_markers: int = 7


Handover = namedtuple("Handover", "rows mass enthalpy volume momentum normal cells")
Handover.__doc__ = """What a mode change passed to the volume, per skin row `rows`: `mass` [kg], `enthalpy` [J] (above T_REF),
`volume` [m3] the mass over the volume's density (beta steps out by it; negative: in), `momentum` [kg m/s] (k, 3) its
part along the surface and `normal` (k, 3) its part along the marker's normal; `cells` the cells handed, (mass, enthalpy)
(k, K) each in its place in the column, for the volume to take their profile (None for `to_thin`). Negative mass and
enthalpy: from the volume to the skin."""

Measures = namedtuple("Measures", "h kmax slope")
_Points = namedtuple("_Points", "centroid normal")


# -- the thin test per marker ----------------------------------------------------------------------------------------
def measures(mesh, phi, psi, markers):
    """Stage 1's measures at the markers (`thin.patch_measures` with the markers for the patches): (h [m], kmax [1/m],
    slope [-]) per marker, h the depth from phi to psi (infinite where no rigid metal lies beneath)."""
    return Measures(*patch_measures(mesh, phi, psi, _Points(markers.x, markers.normal)))


def melt_depth(mesh, psi, x):
    """The depth [m] (k,) from phi to psi at the points x (k, 3) on phi, as `measures` reads it: psi's P1 value there,
    at least zero; infinite where psi is (no rigid metal)."""
    cell, b = locate(mesh, x)
    h = np.maximum(np.einsum("ni,ni->n", b, finite_values(psi)[mesh.tets[cell]]), 0.0)
    return np.where(h > 0.5 * BIG, np.inf, h)


def candidates(m, settings):
    """Whether the thin test answers thin (full-depth) at each marker, from its `Measures`."""
    with np.errstate(invalid="ignore"):
        return thin_candidate(m.h, m.kmax, m.slope, settings.c_kappa, settings.c_slope)


def adjacency(markers):
    """The markers' neighbours as a sparse (n, n) matrix, from their faces (`Markers.faces`)."""
    n = markers.n
    if markers.faces is None:
        raise ValueError("the markers have no faces: tessellate them first")
    i, j = np.asarray(markers.faces.i), np.asarray(markers.faces.j)
    return coo_matrix((np.ones(2 * len(i)), (np.r_[i, j], np.r_[j, i])), shape=(n, n)).tocsr()


def update(markers, candidate, settings):
    """One step end's hysteresis (stage 1's `ModeState`, keyed by marker): each marker's count of step ends in a row its
    test (`candidate`, True for thin) gave the other answer, and the change of every connected group of at least
    `settings.min_markers` markers whose counts reached their hold (`settings.hold_thick` for a full-depth marker,
    `settings.hold_thin` for a top-zone one). Updates `markers.top` and `markers.hold`; returns which markers changed
    (n,) bool."""
    hold = np.where(markers.top, settings.hold_thin, settings.hold_thick)
    st = ModeState(~markers.top, hold=hold, min_group=settings.min_markers)
    st.count = np.asarray(markers.hold, dtype=np.int64).copy()
    changed = st.update(np.asarray(candidate, bool), adjacency(markers))
    markers.top, markers.hold = ~st.thin, st.count
    return changed


# -- a skin's base, by mode ------------------------------------------------------------------------------------------
def fine_zone_mass(skins):
    """A solid fine zone's mass [kg] (n,) under each column's cross-section at the fine zone's depth: a top-zone skin's
    minimum, the same whatever the column's depth, so that its draws and its returns at the base meet at one mass."""
    fz = skins.settings.fine_zone
    return skins.rho_s * skins.cross_section(np.full((skins.n, 1), fz))[:, 0] * fz


def minimum_mass(skins, top):
    """Each skin's minimum mass [kg] (n,): a solid stub over its base's cross-section (full-depth,
    `Skins.minimum_mass`) or a solid fine zone (top-zone, `fine_zone_mass`)."""
    return np.where(np.asarray(top, bool), fine_zone_mass(skins), skins.minimum_mass())


def draw_cap(skins, top):
    """The temperature [K] (n,) a draw's profile is capped at: T_rigid beneath a full-depth skin's stub (rigid metal),
    none beneath a top-zone skin (the thick melt, drawn at its own temperature)."""
    return np.where(np.asarray(top, bool), np.inf, skins.material.T_rigid)


def hand_down(skins, top, stub, cells=False):
    """One part's return at the bases (module docstring): a full-depth skin's anchored rigid cells below `stub` [m]
    (`Skins.hand_down`), a top-zone skin's cells below its fine zone whatever their phase (`Skins.take_bottom`), each
    keeping its minimum mass. Returns (mass [kg], enthalpy [J], momentum [kg m/s] (n, 3)) per row: the cells' own
    momentum (anchored cells are at rest after a momentum step, but a cell that froze since may still hold some); with
    `cells` also the top-zone skins' cells handed, (mass, enthalpy) (n, K) in place (`Skins.take_bottom`)."""
    top = np.asarray(top, bool)
    keep = minimum_mass(skins, top)
    mass, H, P, given = skins.take_bottom(np.where(top, skins.settings.fine_zone, np.inf), keep, cells=True)
    P_held = skins.totals()[2]
    m_f, H_f, _ = skins.hand_down(stub, m_keep=np.where(top, np.inf, keep))                   # the full-depth skins
    out = (mass + m_f, H + H_f, P + (P_held - skins.totals()[2]))
    return out + (given,) if cells else out


# -- the mode changes ------------------------------------------------------------------------------------------------
def _split(P, normal):
    along = np.einsum("ij,ij->i", P, normal)[:, None] * normal
    return P - along, along


def to_thick(skins, i, normal=None, rho=None):
    """Thin to thick (module docstring): the skins rows `i` give their cells below the fine zone to the volume, keeping
    at least a solid fine zone's mass under the cross-section at the fine zone's depth. `normal` (k, 3) the markers'
    unit normals, which split the momentum (without them it is all along the surface); `rho` [kg/m3] the volume's
    density (default the skins' solid density). Returns a Handover. Applied at every step end to every top-zone skin, it
    also returns what the re-partition gave a top-zone skin below its fine zone."""
    i = np.asarray(i, dtype=np.int64).reshape(-1)
    rho = skins.rho_s if rho is None else float(rho)
    depth = np.full(skins.n, np.inf)
    depth[i] = skins.settings.fine_zone
    mass, H, P, (m_c, H_c) = skins.take_bottom(depth, fine_zone_mass(skins), cells=True)
    P = P[i]
    if normal is None:
        along, normal_part = P, np.zeros_like(P)
    else:
        along, normal_part = _split(P, np.asarray(normal, dtype=float).reshape(-1, 3))
    return Handover(i, mass[i], H[i], mass[i] / rho, along, normal_part, (m_c[i], H_c[i]))


def _cells(skins, rows, depth):
    """The new cells' spans below each column's base [m]: (z_top, z_bot) (k, P) from its present thickness down
    `depth` (k,), each near the grading's target at its depth (the last joined to its neighbour when under half a
    target), zero-length past the end."""
    s = skins.settings
    z0 = skins.thickness()[rows]
    tops, bots = [], []
    for a, d in zip(z0, depth):
        end, z, t_k, b_k = a + d, a, [], []
        while d > 0.0 and z < end * (1.0 - 1e-15):
            t = s.cell + (s.growth - 1.0) * max(z - s.fine_zone, 0.0)
            nxt = z + t
            if end - nxt < 0.5 * t:
                nxt = end
            t_k.append(z)
            b_k.append(nxt)
            z = nxt
        tops.append(t_k)
        bots.append(b_k)
    P = max([len(t) for t in tops] + [1])
    zt, zb = np.zeros((len(rows), P)), np.zeros((len(rows), P))
    for k, (t, b) in enumerate(zip(tops, bots)):
        zt[k, :len(t)], zb[k, :len(b)] = t, b
        zt[k, len(t):] = zb[k, len(b):] = b[-1] if b else z0[k]
    return zt, zb


def to_thin(skins, i, depth, temperature, rho=None, kappa_share=None):
    """Thick to thin (module docstring): the skins rows `i` absorb the volume's metal from their bases down `depth` [m]
    (k,): new cells at rest, graded as the skins' cells, each holding the volume's density `rho` (default the skins'
    solid density) times its span's volume under the column's cross-section A0 (1 - k1 z)(1 - k2 z), at the temperature
    `temperature(k, s)` [K] gives at its middle, k the entry of `i` and s [m] the depth below the base (the volume's field
    along the marker's normal). With `kappa_share` each column's curvatures are first held within kappa_share over its
    new depth, so that its cross-section stays open to its bottom (the body's `stepend.cap_curvature` rule). Returns a
    Handover, negative (from the volume), with no momentum (the volume is held)."""
    i = np.asarray(i, dtype=np.int64).reshape(-1)
    depth = np.maximum(np.asarray(depth, dtype=float).reshape(-1), 0.0)
    rho = skins.rho_s if rho is None else float(rho)
    if kappa_share is not None and len(i):
        lim = kappa_share / np.maximum(skins.thickness()[i] + depth, skins.settings.stub)
        skins.kappa[i] = np.clip(skins.kappa[i], -lim[:, None], lim[:, None])
    zt, zb = _cells(skins, i, depth)
    k1, k2 = skins.kappa[i, 0][:, None], skins.kappa[i, 1][:, None]
    vol = skins.area[i][:, None] * (volume_per_area(zb, k1, k2) - volume_per_area(zt, k1, k2))
    m = np.where(zb > zt, rho * vol, 0.0)
    mid = 0.5 * (zt + zb) - skins.thickness()[i][:, None]
    kk = np.broadcast_to(np.arange(len(i))[:, None], m.shape)
    some = m > 0.0
    H = np.zeros(m.shape)
    if some.any():
        H[some] = m[some] * skins.props.enthalpy(np.asarray(temperature(kk[some], mid[some]), dtype=float))
    mass, enthalpy = m.sum(axis=1), H.sum(axis=1)
    skins.append_bottom(i, m, H)
    zero = np.zeros((len(i), 3))
    return Handover(i, -mass, -enthalpy, -mass / rho, zero, zero.copy(), None)
