"""spheral_frag.runner.nodes and .checkpoint: node lists from plain arrays, permanent IDs, checkpoints (Spheral M2 plan,
Task 1; decision D3).

A 5 mm lattice sphere of radius 50 mm (3D: 4,169 particles; RZ: 156 rings) is written here as plain arrays in
drama_env's numpy, with shuffled rows, scattered pids and per-particle masses, densities and energies that no
generator would produce, and built under Spheral's Python on 1 and 4 processes by
`tests/fixtures/spheral_frag/runner_nodes_check.py`:

  - every array arrives bitwise (gathered to rank 0 in pid order), in RZ with massRZ = m / (2 pi r) and
    massDensityRZ = rho;
  - the pids follow deleteNodes: each rank deletes its share of every 13th particle and the rest keep their values;
  - a run of 20 steps (SolidSPH, CheapSynchronousRK2, Murnaghan + ConstantStrength under an imposed expansion, so
    that the deviatoric stress and the plastic strain are non-zero) with a checkpoint after step 10 equals, bitwise,
    10 steps continued from that checkpoint on the same process count, in every field, the time and the time step;
  - the checkpoint written on 4 processes, read on 1, rebuilds the same particles in pid order bitwise; and it is a
    plain numpy-1.26 .npz that drama_env's numpy 2 reads without pickles.

The continuation is bitwise because the uninterrupted run, too, restarts at the checkpoint (nodes.rebuild +
nodes.prime_derivatives): CheapSynchronousRK2's first half step uses the derivatives left from the previous step's
midpoint, which no node field carries. Measured 2026-10-09 on 4 processes in 3D at this size: continuing in memory
instead (no restart) differs from the restarted run after 10 more steps by up to 3.7e-3 of max |S| and 9.3e-3 of
max |v| pointwise; re-evaluating the derivatives in place without rebuilding still differs from the resumed run by
1.8e-3 of max |S| (see the fixture's --boundary none | prime)."""
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

from conftest import spheral_launcher
from helpers import REPO_ROOT

from spheral_frag.runner import checkpoint, nodes

SCRIPT = "tests/fixtures/spheral_frag/runner_nodes_check.py"
R, DX, EDOT = 0.05, 0.005, 100.0          # m, m, 1/s
STEPS, CK_AT = 20, 10
PROCS = (1, 4)
FORMS = ("3d", "rz")


def lattice_sphere(form, seed=3):
    """Plain arrays of a lattice sphere: 3D a cubic lattice with a node at the centre; RZ x along the axis and
    r = (j + 1/2) dx. Rows shuffled, pids scattered, masses and densities perturbed per particle."""
    k = int(math.floor(R/DX*(1 + 1e-12)))
    if form == "3d":
        g = np.arange(-k, k + 1)*DX
        X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    else:
        X = np.stack(np.meshgrid(np.arange(-k, k + 1)*DX, (np.arange(0, k + 1) + 0.5)*DX, indexing="ij"),
                     -1).reshape(-1, 2)
    X = X[np.einsum("ij,ij->i", X, X) <= (R*(1 + 1e-9))**2]
    V = np.full(len(X), DX**3) if form == "3d" else 2*math.pi*X[:, 1]*DX**2
    rng = np.random.default_rng(seed)
    n = len(X)
    perm = rng.permutation(n)
    X, V = X[perm], V[perm]
    rho = 2813.0*(1 + 1e-4*np.sin(37.0*X[:, 0]))
    m = rho*V*(1 + 1e-3*rng.standard_normal(n))
    return dict(x_ref=X, m=m, rho=rho, u=1e5*(1 + rng.random(n)), pid=(1000 + 7*rng.permutation(n)).astype(np.int64),
                v=EDOT*X)


def launch(*args, nproc=1, timeout=900):
    cmd = [spheral_launcher()] + (["-n", str(nproc)] if nproc > 1 else []) + [SCRIPT, *args]
    r = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r


def rel(p):
    return os.path.relpath(p, REPO_ROOT)


def sub(npz, prefix):
    return {k[len(prefix) + 1:]: npz[k] for k in npz.files if k.startswith(prefix + ".")}


@pytest.fixture(scope="module")
def work():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")      # inside the repository: the container sees it
    os.makedirs(base, exist_ok=True)
    d = tempfile.mkdtemp(prefix="runner_nodes_", dir=base)
    inputs = {}
    for form in FORMS:
        a = lattice_sphere(form)
        p = os.path.join(d, f"in_{form}.npz")
        np.savez(p, dx=np.array(DX), **a)
        np.savez(os.path.join(d, f"del_{form}.npz"), pid=a["pid"][::13])
        inputs[form] = a
    yield d, inputs
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(scope="module")
def runs(work):
    """(form, nproc) -> (uninterrupted run, continued run, checkpoint path), launched once each."""
    d, _ = work
    cache = {}

    def get(form, nproc):
        if (form, nproc) not in cache:
            ck = os.path.join(d, f"ck_{form}_n{nproc}.npz")
            full, cont = (os.path.join(d, f"{w}_{form}_n{nproc}.npz") for w in ("full", "cont"))
            launch("run", "--form", form, "--input", rel(os.path.join(d, f"in_{form}.npz")), "--steps", str(STEPS),
                   "--checkpoint-at", str(CK_AT), "--checkpoint-out", rel(ck), "--out", rel(full), nproc=nproc)
            launch("run", "--checkpoint", rel(ck), "--steps", str(STEPS - CK_AT), "--out", rel(cont), nproc=nproc)
            cache[(form, nproc)] = (dict(np.load(full)), dict(np.load(cont)), ck)
        return cache[(form, nproc)]
    return get


# ---------------------------------------------------------------------------------------------------------------
# Under Spheral
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.spheral
@pytest.mark.parametrize("nproc", PROCS)
@pytest.mark.parametrize("form", FORMS)
def test_arrays_arrive_bitwise_and_ids_follow_delete(work, form, nproc):
    d, inputs = work
    out = os.path.join(d, f"build_{form}_n{nproc}.npz")
    launch("build", "--form", form, "--input", rel(os.path.join(d, f"in_{form}.npz")),
           "--delete", rel(os.path.join(d, f"del_{form}.npz")), "--out", rel(out), nproc=nproc)
    z = np.load(out)
    a = inputs[form]
    n, dim = len(a["pid"]), nodes.DIM[form]
    order = np.argsort(a["pid"])
    exp = {k: v[order] for k, v in a.items()}

    counts = z["counts"]
    assert int(z["processes"]) == nproc and counts.shape == (nproc, 2)
    assert counts[:, 0].sum() == n and np.array_equal(counts[:, 0], counts[:, 1])   # pid field length = local count
    if nproc > 1:
        assert counts[:, 0].min() > 0

    b = sub(z, "built")
    for k in ("pid", "x_ref", "m", "rho", "u", "v"):
        assert b[k].dtype == exp[k].dtype and np.array_equal(b[k], exp[k]), k
    assert np.array_equal(b["x"], exp["x_ref"])
    H = np.zeros((n, dim, dim))
    H[:, np.arange(dim), np.arange(dim)] = 1.0/(nodes.NPERH*DX)
    assert np.array_equal(b["H"], H)
    assert not b["S"].any() and not b["plastic_strain"].any() and not b["plastic_strain_rate"].any()
    if form == "rz":
        assert np.array_equal(b["m_rz"], exp["m"]/(2.0*math.pi*exp["x_ref"][:, 1]))
        assert np.array_equal(b["rho_rz"], exp["rho"])
    else:
        assert "m_rz" not in b

    dl = np.load(os.path.join(d, f"del_{form}.npz"))["pid"]
    keep = ~np.isin(exp["pid"], dl)
    after = z["counts_after"]
    assert after[:, 0].sum() == n - len(dl) and np.array_equal(after[:, 0], after[:, 1])
    g = sub(z, "deleted")
    for k in ("pid", "x_ref", "m", "rho", "u", "v"):
        assert np.array_equal(g[k], exp[k][keep]), k
    assert np.array_equal(g["x"], exp["x_ref"][keep])


@pytest.mark.spheral
@pytest.mark.parametrize("nproc", PROCS)
@pytest.mark.parametrize("form", FORMS)
def test_continued_from_checkpoint_equals_uninterrupted_bitwise(runs, form, nproc):
    full, cont, _ = runs(form, nproc)
    assert int(full["processes"]) == int(cont["processes"]) == nproc
    assert int(full["steps"]) == int(cont["steps"]) == STEPS
    assert int(full["cycle"]) == int(cont["cycle"]) == STEPS
    for k in ("time", "dt"):
        assert full[k] == cont[k], k
    keys = [k for k in full if k.startswith("final.")]
    assert len(keys) == len(nodes.STATE_KEYS) + (len(nodes.RZ_KEYS) if form == "rz" else 0)
    for k in keys:
        assert np.array_equal(full[k], cont[k]), k
    # The run exercised what the checkpoint carries: stress, plasticity, motion.
    assert np.abs(full["final.S"]).max() > 1e6 and full["final.plastic_strain"].max() > 0.0
    assert not np.array_equal(full["final.x"], full["final.x_ref"])
    if form == "rz":
        # the hoop deviatoric stress (SymTensor2d's zz element) is evolved and carried (checkpoint schema 2, Task 5)
        assert np.abs(full["final.S_zz"]).max() > 1e6


@pytest.mark.spheral
@pytest.mark.parametrize("nproc", PROCS)
@pytest.mark.parametrize("form", FORMS)
def test_checkpoint_is_plain_npz_with_schema(runs, work, form, nproc):
    _, inputs = work
    _, _, ck = runs(form, nproc)
    with np.load(ck, allow_pickle=False) as z:             # numpy 2 reads Spheral's numpy-1.26 file, no pickles
        assert all(z[k].dtype != object for k in z.files)
    c = checkpoint.read(ck)
    assert c["schema"] == checkpoint.SCHEMA and c["schema_version"] == checkpoint.SCHEMA_VERSION
    assert c["run"] == "runner_nodes_check" and c["k"] == CK_AT and c["form"] == form
    assert c["processes"] == nproc and c["dx"] == DX and c["nperh"] == nodes.NPERH and c["scale"] == 1.0
    assert c["extra"]["steps"] == CK_AT and c["extra"]["cycle"] == CK_AT and c["extra"]["time"] > 0.0
    assert np.array_equal(c["pid"], np.sort(inputs[form]["pid"]))
    a = {k: v[np.argsort(inputs[form]["pid"])] for k, v in inputs[form].items()}
    assert np.array_equal(c["x_ref"], a["x_ref"])
    if form == "3d":
        assert np.array_equal(c["m"], a["m"])                     # 3D masses never change
    else:
        # RZ: massRZ is constant bitwise, but Spheral moves some ring masses by the last bit (measured 2026-10-09:
        # 4 of 156 at the controller's start-up, 12 of 156 after 10 steps, at most 2.2e-16 relative; source not
        # located). Threshold: two units in the last place.
        assert np.array_equal(c["m_rz"], a["m"]/(2.0*math.pi*a["x_ref"][:, 1]))
        assert np.abs(c["m"]/a["m"] - 1.0).max() <= 4.5e-16
    assert np.abs(c["S"]).max() > 1e6 and not np.array_equal(c["x"], c["x_ref"])   # a moving, stressed state


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_checkpoint_from_4_processes_rebuilds_on_1_bitwise(runs, work, form):
    d, _ = work
    _, _, ck = runs(form, 4)
    out = os.path.join(d, f"rebuild_{form}.npz")
    launch("rebuild", "--checkpoint", rel(ck), "--out", rel(out), nproc=1)
    z = np.load(out)
    assert int(z["processes"]) == 1
    c = checkpoint.read(ck)
    b = sub(z, "built")
    keys = nodes.STATE_KEYS + (nodes.RZ_KEYS if form == "rz" else ())
    assert set(b) == set(keys)
    for k in keys:
        assert np.array_equal(b[k], c[k]), k


# ---------------------------------------------------------------------------------------------------------------
# Without Spheral (numpy only)
# ---------------------------------------------------------------------------------------------------------------

def test_checkpoint_read_refuses_another_schema(tmp_path):
    p = tmp_path / "x.npz"
    np.savez(p, schema=np.array("something.else"), schema_version=np.array(1))
    with pytest.raises(ValueError, match="not a"):
        checkpoint.read(str(p))
    np.savez(p, schema=np.array(checkpoint.SCHEMA), schema_version=np.array(checkpoint.SCHEMA_VERSION + 1))
    with pytest.raises(ValueError, match="newer"):
        checkpoint.read(str(p))


def test_check_arrays_rejects_bad_input():
    a = lattice_sphere("rz")
    nodes._check_arrays("rz", a)
    with pytest.raises(ValueError, match="unique"):
        nodes._check_arrays("rz", dict(a, pid=np.zeros_like(a["pid"])))
    with pytest.raises(ValueError, match="integers"):
        nodes._check_arrays("rz", dict(a, pid=a["pid"].astype(float)))
    bad = a["x_ref"].copy()
    bad[0, 1] = 0.0
    with pytest.raises(ValueError, match="r > 0"):
        nodes._check_arrays("rz", dict(a, x_ref=bad))
    with pytest.raises(ValueError, match="shape"):
        nodes._check_arrays("rz", dict(a, m=a["m"][:-1]))
    with pytest.raises(ValueError, match="lack"):
        nodes._check_arrays("rz", {k: v for k, v in a.items() if k != "u"})
    with pytest.raises(ValueError, match="RZ-only"):
        nodes._check_arrays("3d", dict(lattice_sphere("3d"), m_rz=np.ones(4169)))


def test_global_id_range_covers_everything_once():
    for n in (0, 1, 7, 156, 4169):
        for procs in (1, 2, 4, 18):
            parts = [nodes.global_id_range(n, r, procs) for r in range(procs)]
            assert parts[0][0] == 0 and parts[-1][1] == n
            assert all(parts[i][1] == parts[i + 1][0] for i in range(procs - 1))
            sizes = [hi - lo for lo, hi in parts]
            assert max(sizes) - min(sizes) <= 1
