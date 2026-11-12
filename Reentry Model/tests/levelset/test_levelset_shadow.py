import numpy as np

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.shadow import shadowed


C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _mesh(h=0.0625e-3, box=(3e-3, 3e-3, 3e-3)):
    return build(BoxBase.around((0, 0, 0), box, h))


def test_a_sphere_shadows_exactly_its_lee_side():
    m = _mesh(0.125e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.0e-3)
    P = patches(m, phi)
    hit = shadowed(m, phi, P, (1.0, 0.0, 0.0))
    nx = P.normal[:, 0]
    assert hit[nx < -0.2].all() and not hit[nx > 0.2].any()


def test_an_overhang_shelters_a_forward_facing_ledge_and_one_piece_hides_another():
    m = _mesh(0.0625e-3, box=(3e-3, 3e-3, 3e-3))
    x = m.points
    rr = np.linalg.norm(x[:, 1:] - 1.5e-3, axis=1)
    cyl = lambda r, x0, x1: np.maximum(rr - r, np.maximum(x0 - x[:, 0], x[:, 0] - x1))
    body = np.minimum(np.minimum(cyl(0.3e-3, 0.3e-3, 2.3e-3), cyl(0.8e-3, 0.3e-3, 1.0e-3)),
                      cyl(1.0e-3, 2.0e-3, 2.3e-3))                                  # base step, stem, cap
    phi = redistance(m, body)
    P = patches(m, phi)
    hit = shadowed(m, phi, P, (1.0, 0.0, 0.0))
    pr = np.linalg.norm(P.centroid[:, 1:] - 1.5e-3, axis=1)
    ledge = (P.normal[:, 0] > 0.9) & (np.abs(P.centroid[:, 0] - 1.0e-3) < 0.05e-3) & (pr > 0.4e-3) & (pr < 0.7e-3)
    cap_face = (P.normal[:, 0] > 0.9) & (np.abs(P.centroid[:, 0] - 2.3e-3) < 0.05e-3) & (pr < 0.9e-3)
    assert ledge.sum() > 10 and hit[ledge].all()
    assert cap_face.sum() > 10 and not hit[cap_face].any()


def test_a_ray_grazing_a_long_face_parallel_to_the_flow_finishes_when_one_cell_is_refined_to_a_small_floor():
    base = BoxBase.around((0, 0, 0), (30e-3, 4e-3, 4e-3), 0.25e-3)
    small = lambda x, s: np.where(np.linalg.norm(x - np.array([29e-3, 3.5e-3, 3.5e-3]), axis=1) < 0.2e-3,
                                  0.0625e-3, np.inf)
    m = build(base, small)
    phi = redistance(m, m.points[:, 2] - 2.0e-3)          # metal below a flat face along the flow
    P = patches(m, phi)
    hit = shadowed(m, phi, P, (1.0, 0.0, 0.0))
    assert len(P.cell) > 0 and not hit.any()
