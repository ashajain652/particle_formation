"""The skins' upstream-ray shelter test (`reentry_model/skyscrapers/shelter.py`, design 2026-11-13-fitted-sphere-loads).

The ring test raises the tops of the corners whose base normal is 45 to 52 degrees from the flight direction by
RING_RAISE [m] along their base normals; RING_RAISE is 6 mm, not the brief's 4 mm:
4 mm leaves no ring corner whose top normal faces away from the flow (theta > pi/2) on this 4 mm surface mesh, and 6 mm
leaves 24 of them (8 mm, 42), so 6 mm is the least that tests the exposed crest."""
import numpy as np
import pytest

from reentry_model import columns
from reentry_model.skyscrapers import shelter
from test_reentry_model_skinbody import make

RING_RAISE = 6e-3
VS = [np.array([-1.0, 0, 0]), np.array([1.0, 0, 0]), np.array([0.6, -0.3, 0.74])]


def _unit(v):
    return v / np.linalg.norm(v)


def _sphere(m, v):
    """Tops of the plain sphere (the surface corners moved out 0.8 mm) and the base angle [deg] of each corner."""
    b = make(m)
    v = _unit(v)
    base = b.surface.points.copy()
    n = b.columns.normal
    return b, base, n, np.degrees(np.arccos(np.clip(n @ v, -1, 1)))


def _args(b, pts):
    s = columns.CornerSurface.from_points(pts, b.columns.faces)
    return s.points, s.faces, s.normals, s.areas


def _ahead(b, pts, v):
    return (pts - b.mass_centroid()) @ _unit(v) > 0


@pytest.mark.parametrize("v", VS)
def test_no_corner_of_a_sphere_is_sheltered(coarse_sphere_mesh, v):
    b, pts, n, ang = _sphere(coarse_sphere_mesh, v)
    A = _args(b, pts)
    ahead = _ahead(b, pts, v)
    assert ahead.any()
    assert not shelter.upstream_hits(*A, v, candidates=ahead).any()


def _ringed(m, v):
    b, pts, n, ang = _sphere(m, v)
    ring = (ang >= 45) & (ang <= 52)
    pts = pts.copy()
    pts[ring] += RING_RAISE * n[ring]
    return b, pts, ang, ring


@pytest.mark.parametrize("v", VS[:1] + VS[2:])
def test_a_raised_ring_shelters_the_band_behind_it_and_exposes_its_crest(coarse_sphere_mesh, v):
    b, pts, ang, ring = _ringed(coarse_sphere_mesh, v)
    A = _args(b, pts)
    ahead = _ahead(b, pts, v)
    theta = np.arccos(np.clip(A[2] @ _unit(v), -1, 1))
    hit = shelter.upstream_hits(*A, v, candidates=ahead)
    assert not hit[ang < 30].any()
    assert hit[(ang > 52) & (ang < 65)].any()
    assert (ang[hit] >= 45).all() and (ang[hit] <= 75).all()
    assert (ring & (theta > np.pi / 2) & ~hit).any()


@pytest.mark.parametrize("v", [VS[0], VS[2]])
def test_the_grid_gives_the_plain_search_s_answer(coarse_sphere_mesh, v):
    b, pts, ang, ring = _ringed(coarse_sphere_mesh, v)
    rng = np.random.default_rng(7)
    jit = pts + rng.uniform(-0.3e-3, 0.3e-3, pts.shape)
    for p in (pts, jit):
        A = _args(b, p)
        ahead = _ahead(b, p, v)
        assert np.array_equal(shelter.upstream_hits(*A, v, candidates=ahead),
                              shelter.upstream_hits_plain(*A, v, candidates=ahead))
    assert shelter.upstream_hits(*_args(b, pts), v, candidates=_ahead(b, pts, v)).any()


def test_only_candidates_are_tested(coarse_sphere_mesh):
    b, pts, ang, ring = _ringed(coarse_sphere_mesh, VS[0])
    A = _args(b, pts)
    n = len(pts)
    assert not shelter.upstream_hits(*A, VS[0], candidates=np.zeros(n, bool)).any()
    assert not shelter.upstream_hits_plain(*A, VS[0], candidates=np.zeros(n, bool)).any()
