import numpy as np
from scipy.sparse import diags

from lsproto.boxmesh import BoxBase, build
from lsproto.cutcell import patches
from lsproto.levelset import redistance
from lsproto.shadow import shadowed
from lsproto.thin import ModeState, patch_adjacency, patch_measures, thin_candidate

C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _mesh(h=0.0625e-3, box=(3e-3, 3e-3, 3e-3)):
    return build(BoxBase.around((0, 0, 0), box, h))


def test_the_decision_uses_the_largest_principal_curvature_not_the_mean():
    # a saddle (k1 = +1/mm, k2 = -1/mm) has zero mean curvature; 0.3 mm of melt on it is not thin
    assert not thin_candidate(0.3e-3, max(abs(1e3), abs(-1e3)), 0.0)
    assert thin_candidate(0.05e-3, 1e3, 0.0)
    assert not thin_candidate(0.05e-3, 1e3, 0.3)


def test_a_thin_shell_is_thin_and_a_deep_pool_is_thick():
    m = _mesh()
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.1e-3)
    P = patches(m, phi)
    for depth, want in ((0.06e-3, True), (0.5e-3, False)):
        psi = np.maximum(redistance(m, r - (1.1e-3 - depth)), phi)
        h, k, s = patch_measures(m, phi, psi, P)
        assert np.median(h) == np.median(h) and abs(np.median(h) - depth) < 0.02e-3
        assert np.mean(thin_candidate(h, k, s)) > 0.95 if want else np.mean(thin_candidate(h, k, s)) < 0.05


def test_a_floor_tilted_against_the_surface_fails_on_slope():
    m = _mesh(0.125e-3)
    phi = redistance(m, m.points[:, 2] - 2.0e-3)
    psi = np.maximum(redistance(m, m.points[:, 2] - (1.5e-3 + 0.3 * (m.points[:, 0] - 1.5e-3))), phi)
    P = patches(m, phi)
    mid = (np.abs(P.centroid[:, 0] - 1.5e-3) < 0.5e-3) & (np.abs(P.centroid[:, 1] - 1.5e-3) < 1e-3)
    h, k, s = patch_measures(m, phi, psi, P)
    assert np.allclose(s[mid], 0.3 / np.sqrt(1 + 0.09), atol=0.03)


def test_mode_changes_wait_for_the_hold_and_a_large_enough_group():
    n = 20
    adj = (diags(np.ones(n - 1), 1) + diags(np.ones(n - 1), -1)).tocsr()          # a chain of patches
    st = ModeState(np.ones(n, bool), hold=3, min_group=7)
    cand = np.ones(n, bool)
    cand[2:10] = False                                                             # 8 in a row want thick
    cand[15] = False                                                               # a lone one wants thick
    assert not st.update(cand, adj).any()
    assert not st.update(cand, adj).any()
    changed = st.update(cand, adj)
    assert changed[2:10].all() and not changed[15] and changed.sum() == 8
    assert (~st.thin[2:10]).all() and st.thin[15]


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
