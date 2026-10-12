"""spheral_frag.runner.packages.hourglass: the stiffness against SolidSPH's tensile pairing (Spheral M2 Task 10b).

The fit operator A = I - C (numpy and scipy only, no Spheral): (A u)_i = u_i - sum_j c_ij u_j with c_ij the affine
least-squares weights of i's neighbours' displacements at i's reference position. On lattice balls in RZ's half
plane and in 3D it is exact on every affine field (translation, rotation, uniform strain) to round-off, including
the staircase's corners (which take a wider neighbourhood); A 1 = 0, so the forces -A^T (k A u) sum to zero
(momentum); and the alternating plane mode of check 13 (whole lattice planes x = const shifted alternately along
the axis) is seen at order one. The measured tables are tested in test_spheral_frag_runner_task10b.py."""
import numpy as np
import pytest

from spheral_frag.runner.packages.hourglass import fit_operator


def lattice(dim, n, half_plane=False):
    g = np.arange(-n, n + 1, dtype=float)
    X = np.array(np.meshgrid(*[g]*dim, indexing="ij")).reshape(dim, -1).T
    X = X[np.linalg.norm(X, axis=1) <= n - 0.3]
    if half_plane:                                   # RZ: rows at r = (j + 1/2) dx
        X = X[X[:, 1] >= 0.0] + np.array([0.0, 0.5])
    return X


@pytest.mark.parametrize("dim,half", [(2, True), (2, False), (3, False)])
def test_affine_fields_are_reproduced_exactly(dim, half):
    rng = np.random.default_rng(7)
    X = lattice(dim, 9, half)
    A, widened = fit_operator(X, 1.0)
    for _ in range(3):
        u = X @ rng.normal(size=(dim, dim)).T + rng.normal(size=dim)
        assert np.abs(A @ u).max() < 1e-12
    assert np.abs(A @ np.ones(len(X))).max() < 1e-12
    # nearest and diagonal neighbours: 8 (2D) and 18 (3D) in the interior, so 9 and 19 entries per row there (in
    # RZ away from the axis: the axis row's fit weights only its two axial neighbours, the row above at zero)
    r = np.linalg.norm(X, axis=1)
    inner = np.flatnonzero((r < 5.0) & ((X[:, 1] > 1.0) if half else True))
    nnz = np.diff(A.indptr)[inner]
    assert set(nnz.tolist()) == {9 if dim == 2 else 19}
    assert widened >= 0


@pytest.mark.parametrize("dim", [2, 3])
def test_forces_sum_to_zero_and_vanish_on_affine_fields(dim):
    rng = np.random.default_rng(11)
    X = lattice(dim, 8)
    A, _ = fit_operator(X, 1.0)
    k = rng.uniform(0.5, 2.0, len(X))
    u = rng.normal(size=X.shape)*1e-2
    F = -(A.T @ (k[:, None]*(A @ u)))
    assert np.abs(F.sum(axis=0)).max() < 1e-13*np.abs(F).max()*len(X)
    aff = X @ rng.normal(size=(dim, dim)).T
    assert np.abs(A.T @ (k[:, None]*(A @ aff))).max() < 1e-11
    # the energy 1/2 sum k |A u|^2 is positive for a non-affine field
    assert float(np.sum(k*np.einsum("ij,ij->i", A @ u, A @ u))) > 0.0


@pytest.mark.parametrize("dim", [2, 3])
def test_the_alternating_plane_mode_is_seen(dim):
    """Check 13's pairing: planes x = const shifted alternately along x by +-a. |A u| is of order a wherever a
    particle's neighbourhood spans three planes (at the ball's tips two planes make the pattern affine), against
    (dx/L)^2-small for a smooth quadratic field of the same amplitude."""
    X = lattice(dim, 10)
    A, _ = fit_operator(X, 1.0)
    z = np.zeros_like(X)
    z[:, 0] = (-1.0)**np.round(X[:, 0]).astype(int)
    inner = np.linalg.norm(X, axis=1) < 7.0
    d = np.abs((A @ z)[inner, 0])
    assert d.min() > 0.5 and np.mean(d) > 1.0
    q = np.zeros_like(X)
    q[:, 0] = (X[:, 0]/10.0)**2                            # amplitude 1 over the ball's radius
    assert np.abs((A @ q)[inner, 0]).max() < 0.05
