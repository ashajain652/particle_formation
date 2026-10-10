"""Solid node lists built from plain arrays, permanent particle IDs and gathers in ID order (Spheral M2 plan, Task 1).

A body arrives as plain arrays (the body stage's output, or a checkpoint): reference positions `x_ref` (n, dim),
masses `m`, densities `rho`, specific thermal energies `u` and permanent IDs `pid` (unique integers), optionally the
current state of a checkpoint (`x`, `v`, `H`, `S`, `plastic_strain`, `plastic_strain_rate`, and in RZ `m_rz`,
`rho_rz`, `S_zz`, `H_zz`: the RZ tensors' hoop element, see RZ_KEYS). `make_body` builds a solid node list from them and distributes it over the processes with Spheral's
`distributeNodes` (Peano-Hilbert), through `ArrayGenerator`, a node generator that serves the arrays the way Spheral's
own generators serve theirs. The pids travel through the redistribution as one of distributeNodesGeneric's "extra"
lists; after it every local particle's values are written again from the arrays, looked up by pid, so that the node
list holds the arrays' values bitwise: no generator's mass is used (plan fact 3), and neither distributeNodesGeneric's
RZ mass (it multiplies the served mass by 2 pi r) nor its clamped H tensors survive.

Forms. "3d": Spheral3d, positions (x, y, z). "rz": SpheralRZ, x the flight axis and y = r > 0; a particle's mass is
its whole ring's (plan fact 2), its `massRZ` is m / (2 pi r) and `massDensityRZ` starts equal to rho, as Spheral's
RZ generators set them (distributeNodesGeneric.py). SolidSPHRZ evolves massDensityRZ and holds massRZ as state, so a
checkpoint carries both. massRZ stays constant bitwise in a run, but Spheral moves some ring masses m by the last bit
(measured 2026-10-09 on the 5 mm sphere, 156 rings: 4 at SpheralController's start-up, 12 after 10 steps, at most
2.2e-16 relative; the source was not located), so RZ accounts that need exact masses take them from a gather, not
from the body's arrays. One form per process: importing SpheralRZ switches Spheral's global coordinate registrar to
RZ (GeometryRegistrar.coords() is RZ after importing both, measured 2026-10-09).

The permanent ID field is a plain IntField named "pid", registered on the node list (plan fact 4: it keeps its length
and order through deleteNodes, and redistribution carries it as it carries distributeNodesGeneric's own extra
fields). The `Body` keeps the Python object, since a Spheral field leaves its node list's registry when it is
destroyed; the reference positions are a VectorField "x_ref" kept the same way.

Frame boundaries are restarts. The node fields are not the whole state of a run: CheapSynchronousRK2 starts each
step from the derivatives the packages kept from the previous step's midpoint. A run continues bitwise from a
checkpoint only when the uninterrupted run restarts there too: `rebuild` (a new body from the gathered state, new
packages, integrator and controller) and `prime_derivatives`. Tested in tests/test_spheral_frag_runner_nodes.py.

Spheral is imported inside the functions (the drama_env import-boundary test inspects this module without it)."""
from __future__ import annotations

import math

import numpy as np

FORMS = ("3d", "rz")
DIM = {"3d": 3, "rz": 2}
NPERH = 2.01                   # nodes per smoothing length (M0's bench.py)
RHO_RANGE_FACTORS = (0.2, 4.0)  # default rhoMin, rhoMax as fractions of the smallest and largest density (bench.py)

REQUIRED = ("x_ref", "m", "rho", "u", "pid")
# Every per-particle array of a checkpoint (plan D3), in the order they are written; m_rz, rho_rz only in RZ.
STATE_KEYS = ("pid", "x_ref", "x", "v", "H", "m", "rho", "u", "S", "plastic_strain", "plastic_strain_rate")
RZ_KEYS = ("m_rz", "rho_rz", "S_zz", "H_zz")
# RZ tensors are Spheral's SymTensor2d, which holds a third diagonal element zz besides xx, xy, yy: SolidSPHRZ
# evolves the deviatoric stress's hoop component in it (SolidSPHRZ.cc: DvDx.zz = v_r / r, the deviatoric
# deformation traceless in 3D), so an RZ state carries S_zz and H_zz beside the (n, 2, 2) tensors.


def spheral_module(form):
    """(Spheral module, distributeNodes) for the form: Spheral3d and distributeNodes3d, or SpheralRZ and
    distributeNodes2d (Peano-Hilbert ordered, as M0's bench.py)."""
    if form == "3d":
        import Spheral3d as S
        from PeanoHilbertDistributeNodes import distributeNodes3d as distribute
    elif form == "rz":
        import SpheralRZ as S
        from PeanoHilbertDistributeNodes import distributeNodes2d as distribute
    else:
        raise ValueError(f"form must be one of {FORMS}, got {form!r}")
    return S, distribute


def _comm():
    from mpi4py import MPI
    return MPI.COMM_WORLD


# ---------------------------------------------------------------------------------------------------------------
# Spheral values <-> numpy rows
# ---------------------------------------------------------------------------------------------------------------

class ZZ:
    """The zz element of an RZ SymTensor field, as a field of its own for `field_to_numpy` and `set_field`."""

    def __init__(self, field):
        self.field = field
        self.name = field.name + ".zz"

    @property
    def numInternalElements(self):
        return self.field.numInternalElements


def _kind(field):
    if isinstance(field, ZZ):
        return "ZZ"
    name = type(field).__name__
    for k in ("SymTensor", "Vector", "Scalar", "Int"):
        if name.startswith(k + "Field"):
            return k
    raise TypeError(f"unsupported field type {name}")


def field_to_numpy(field, dim):
    """The internal values of a Spheral field as numpy: Scalar (n,), Int (n,) int64, Vector (n, dim),
    SymTensor (n, dim, dim)."""
    kind = _kind(field)
    if kind == "ZZ":
        return np.array([t.zz for t in field.field.internalValues()], dtype=np.float64)
    vals = field.internalValues()
    if kind == "Scalar":
        return np.array(list(vals), dtype=np.float64)
    if kind == "Int":
        return np.array(list(vals), dtype=np.int64)
    if kind == "Vector":
        if dim == 3:
            rows = [(v.x, v.y, v.z) for v in vals]
        else:
            rows = [(v.x, v.y) for v in vals]
        return np.array(rows, dtype=np.float64).reshape(len(rows), dim)
    if dim == 3:
        rows = [((t.xx, t.xy, t.xz), (t.yx, t.yy, t.yz), (t.zx, t.zy, t.zz)) for t in vals]
    else:
        rows = [((t.xx, t.xy), (t.yx, t.yy)) for t in vals]
    return np.array(rows, dtype=np.float64).reshape(len(rows), dim, dim)


def to_spheral(S, kind, row):
    """One numpy row as the Spheral value of the field kind (exact: float64 -> Python float)."""
    if kind == "Scalar":
        return float(row)
    if kind == "Int":
        return int(row)
    if kind == "Vector":
        return S.Vector(*[float(c) for c in row])
    return S.SymTensor(*[float(c) for c in np.asarray(row).ravel()])


def set_field(S, field, values):
    """Assign the local internal values of a field from numpy rows (len(values) == numInternalElements)."""
    kind = _kind(field)
    n = field.numInternalElements
    if len(values) != n:
        raise ValueError(f"{field.name}: {len(values)} values for {n} internal nodes")
    if kind == "ZZ":
        f = field.field
        for i in range(n):
            t = f[i]
            f[i] = S.SymTensor(t.xx, t.xy, t.yx, t.yy, float(values[i]))
        return
    for i in range(n):
        field[i] = to_spheral(S, kind, values[i])


def h_tensor_rows(dim, n, dx, nperh):
    """Round H tensors 1/(nperh dx) times the identity, as rows (n, dim, dim)."""
    H = np.zeros((n, dim, dim))
    H[:, np.arange(dim), np.arange(dim)] = 1.0/(nperh*dx)
    return H


# ---------------------------------------------------------------------------------------------------------------
# The generator that serves the arrays
# ---------------------------------------------------------------------------------------------------------------

def global_id_range(n, rank, procs):
    """This rank's slice [lo, hi) of n items split as evenly as possible (NodeGeneratorBase.globalIDRange)."""
    base, rem = divmod(n, procs)
    lo = rank*base + min(rank, rem)
    return lo, lo + base + (1 if rank < rem else 0)


class ArrayGenerator:
    """A node generator over plain arrays, the interface distributeNodesGeneric reads (NodeGeneratorBase at 116c71f):
    localNumNodes, globalNumNodes, localPosition, localMass, localMassDensity, localVelocity, localHtensor,
    localMatE. Each rank serves its contiguous slice of the arrays (which every rank holds in full); `mass` is what
    the generator serves as mass, i.e. massRZ in RZ, which distributeNodesGeneric multiplies by 2 pi r."""

    def __init__(self, S, x, mass, rho, v, H, u):
        self.S = S
        self.x, self.m, self.rho, self.v, self.H, self.u = x, mass, rho, v, H, u

    def localNumNodes(self):
        return len(self.m)

    def globalNumNodes(self):
        return _comm().allreduce(len(self.m))

    def localPosition(self, i):
        return to_spheral(self.S, "Vector", self.x[i])

    def localMass(self, i):
        return float(self.m[i])

    def localMassDensity(self, i):
        return float(self.rho[i])

    def localVelocity(self, i):
        return to_spheral(self.S, "Vector", self.v[i])

    def localHtensor(self, i):
        return to_spheral(self.S, "SymTensor", self.H[i])

    def localMatE(self, i):
        return float(self.u[i])


# ---------------------------------------------------------------------------------------------------------------
# The body
# ---------------------------------------------------------------------------------------------------------------

class Body:
    """A distributed solid node list with its permanent IDs. Attributes: S (the Spheral module), form, dim, nodes,
    db (a DataBase holding the node list), pid (IntField "pid"), x_ref (VectorField "x_ref"), eos, strength and
    kernel (kept alive here: Spheral holds them by reference), dx, nperh."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def state_fields(self):
        """name -> Spheral field for every per-particle array of a checkpoint (STATE_KEYS, + RZ_KEYS in RZ)."""
        n = self.nodes
        f = dict(pid=self.pid, x_ref=self.x_ref, x=n.positions(), v=n.velocity(), H=n.Hfield(), m=n.mass(),
                 rho=n.massDensity(), u=n.specificThermalEnergy(), S=n.deviatoricStress(),
                 plastic_strain=n.plasticStrain(), plastic_strain_rate=n.plasticStrainRate())
        if self.form == "rz":
            f.update(m_rz=n.massRZ(), rho_rz=n.massDensityRZ(), S_zz=ZZ(n.deviatoricStress()), H_zz=ZZ(n.Hfield()))
        return f

    def local_pids(self):
        return field_to_numpy(self.pid, self.dim)

    def delete_by_id(self, pids):
        """Delete the local particles whose pid is in `pids` (each rank deletes its own; collective-free)."""
        local = self.local_pids()
        idx = np.flatnonzero(np.isin(local, np.asarray(pids, dtype=np.int64)))
        if len(idx):
            self.nodes.deleteNodes(self.S.vector_of_ULL([int(i) for i in idx]))
        return len(idx)

    def gather(self, fields=None):
        """gather_by_id of the state fields (or the given names) to rank 0; None on other ranks."""
        sf = self.state_fields()
        sel = sf if fields is None else {k: sf[k] for k in set(fields) | {"pid"}}
        return gather_by_id(self.nodes, sel, self.dim)


def _check_arrays(form, arrays):
    dim = DIM[form]
    for k in REQUIRED:
        if k not in arrays:
            raise ValueError(f"arrays lack {k!r} (need {REQUIRED})")
    n = len(arrays["pid"])
    pid = np.asarray(arrays["pid"])
    if not np.issubdtype(pid.dtype, np.integer):
        raise ValueError("pid must be integers")
    if len(np.unique(pid)) != n:
        raise ValueError("pid must be unique")
    shapes = dict(x_ref=(n, dim), x=(n, dim), v=(n, dim), H=(n, dim, dim), S=(n, dim, dim), m=(n,), rho=(n,),
                  u=(n,), plastic_strain=(n,), plastic_strain_rate=(n,), m_rz=(n,), rho_rz=(n,), S_zz=(n,), H_zz=(n,))
    for k, shp in shapes.items():
        if k in arrays and np.shape(arrays[k]) != shp:
            raise ValueError(f"{k} has shape {np.shape(arrays[k])}, expected {shp}")
    if form == "rz":
        for k in ("x_ref", "x"):
            if k in arrays and not np.all(np.asarray(arrays[k])[:, 1] > 0.0):
                raise ValueError(f"RZ needs y = r > 0 for every particle ({k})")
    elif any(k in arrays for k in RZ_KEYS):
        raise ValueError(f"{RZ_KEYS} are RZ-only")
    return n, dim


def set_ids(S, nodes, pid):
    """The registered IntField "pid" of the node list, filled once with the local particles' permanent IDs (plan
    fact 4). Keep the returned object: the field stays registered (and follows deleteNodes and redistribution)
    only while it lives."""
    f = S.IntField("pid", nodes)
    set_field(S, f, np.asarray(pid, dtype=np.int64))
    return f


def make_body(S, form, arrays, eos, strength, nperh=NPERH, *, dx, kernel=None, name="body", rho_range=None,
              distribute=None):
    """Build and distribute a solid node list from plain arrays (module docstring). Collective: every rank passes
    the same full arrays. dx is the spacing that sets the round H = 1/(nperh dx) when `arrays` has no "H".
    rho_range = (rhoMin, rhoMax), default RHO_RANGE_FACTORS times the smallest and largest density. Returns a
    `Body`."""
    n, dim = _check_arrays(form, arrays)
    if distribute is None:
        _, distribute = spheral_module(form)
    comm = _comm()
    a = {k: np.asarray(v) for k, v in arrays.items()}
    pid = a["pid"].astype(np.int64)
    x_ref = a["x_ref"].astype(np.float64)
    x = a.get("x", x_ref).astype(np.float64)
    v = a.get("v", np.zeros((n, dim))).astype(np.float64)
    H = a.get("H", h_tensor_rows(dim, n, dx, nperh)).astype(np.float64)
    m = a["m"].astype(np.float64)
    rho = a["rho"].astype(np.float64)
    u = a["u"].astype(np.float64)
    zeros = np.zeros(n)
    Sdev = a.get("S", np.zeros((n, dim, dim))).astype(np.float64)
    ps = a.get("plastic_strain", zeros).astype(np.float64)
    psr = a.get("plastic_strain_rate", zeros).astype(np.float64)
    if form == "rz":
        m_rz = a.get("m_rz", m/(2.0*math.pi*x[:, 1])).astype(np.float64)
        rho_rz = a.get("rho_rz", rho).astype(np.float64)
        S_zz = a.get("S_zz", zeros).astype(np.float64)
        H_zz = a.get("H_zz", zeros).astype(np.float64)
    served_mass = m_rz if form == "rz" else m

    if kernel is None:
        kernel = S.TableKernel(S.WendlandC4Kernel(), 1000)
    if rho_range is None:
        rho_range = (RHO_RANGE_FACTORS[0]*float(rho.min()), RHO_RANGE_FACTORS[1]*float(rho.max()))
    extent = 10.0*max(float(np.abs(x).max()), float(np.abs(x_ref).max()), dx)
    box = (S.Vector(*([-extent]*dim)), S.Vector(*([extent]*dim)))
    nodes = S.makeSolidNodeList(name, eos, strength, nPerh=nperh, kernelExtent=kernel.kernelExtent,
                                hmin=1e-3*dx, hmax=100.0*dx, rhoMin=rho_range[0], rhoMax=rho_range[1],
                                xmin=box[0], xmax=box[1])

    lo, hi = global_id_range(n, comm.rank, comm.size)
    sl = slice(lo, hi)
    gen = ArrayGenerator(S, x[sl], served_mass[sl], rho[sl], v[sl], H[sl], u[sl])
    local_pid = [int(p) for p in pid[sl]]
    distribute((nodes, gen, local_pid))        # local_pid now holds this rank's pids after redistribution

    order = np.argsort(pid, kind="stable")
    lp = np.asarray(local_pid, dtype=np.int64)
    rows = order[np.searchsorted(pid, lp, sorter=order)] if len(lp) else np.zeros(0, dtype=np.int64)
    if len(lp) and not np.array_equal(pid[rows], lp):
        raise RuntimeError("a redistributed pid is not among the arrays' pids")
    if nodes.numInternalNodes != len(lp):
        raise RuntimeError(f"{nodes.numInternalNodes} local nodes but {len(lp)} pids after distribution")

    pid_field = set_ids(S, nodes, lp)
    x_ref_field = S.VectorField("x_ref", nodes)
    set_field(S, x_ref_field, x_ref[rows])
    set_field(S, nodes.positions(), x[rows])
    set_field(S, nodes.velocity(), v[rows])
    set_field(S, nodes.Hfield(), H[rows])
    set_field(S, nodes.mass(), m[rows])
    set_field(S, nodes.massDensity(), rho[rows])
    set_field(S, nodes.specificThermalEnergy(), u[rows])
    set_field(S, nodes.deviatoricStress(), Sdev[rows])
    set_field(S, nodes.plasticStrain(), ps[rows])
    set_field(S, nodes.plasticStrainRate(), psr[rows])
    if form == "rz":
        set_field(S, nodes.massRZ(), m_rz[rows])
        set_field(S, nodes.massDensityRZ(), rho_rz[rows])
        set_field(S, ZZ(nodes.deviatoricStress()), S_zz[rows])
        set_field(S, ZZ(nodes.Hfield()), H_zz[rows])
    nodes.neighbor().updateNodes()

    db = S.DataBase()
    db.appendNodeList(nodes)
    return Body(S=S, form=form, dim=dim, nodes=nodes, db=db, pid=pid_field, x_ref=x_ref_field, eos=eos,
                strength=strength, kernel=kernel, dx=dx, nperh=nperh, rho_range=tuple(rho_range), name=name,
                distribute=distribute)


def rebuild(body, name):
    """A new Body from this body's gathered state (what a checkpoint holds), with the same material, kernel,
    spacing and density range: the frame-boundary restart (see `prime_derivatives`). Collective. The old body's
    physics packages, integrator and controller must be rebuilt on the new one; `name` must differ from the old
    node list's while that one lives."""
    state = _comm().bcast(body.gather(), root=0)
    return make_body(body.S, body.form, state, body.eos, body.strength, body.nperh, dx=body.dx, kernel=body.kernel,
                     name=name, rho_range=body.rho_range, distribute=body.distribute)


def prime_derivatives(body, integrator):
    """Evaluate every package's derivatives at the current state, as SpheralController(initializeDerivatives=True)
    does. Call once after building the integrator on a body (fresh, rebuilt or read from a checkpoint).

    Why: CheapSynchronousRK2 takes each step's first half step with the derivatives the packages still hold from
    the previous step's midpoint (CheapSynchronousRK2.cc: no evaluation before `state.update(derivs, hdt, ...)`),
    and SolidSPH keeps them across steps on purpose (SolidSPH.cc registerDerivatives). Those derivatives are hidden
    state that no checkpoint of the node fields carries, so a run continues bitwise from a checkpoint only when the
    uninterrupted run also restarts at that boundary: rebuilt from the gathered state (`rebuild`) and primed here.
    Without the priming a fresh body's first half step uses zero derivatives."""
    pk = integrator.physicsPackages()
    S = body.S
    state = S.State(body.db, pk)
    derivs = S.StateDerivatives(body.db, pk)
    integrator.setGhostNodes()
    integrator.applyGhostBoundaries(state, derivs)
    integrator.preStepInitialize(state, derivs)
    t, dt = integrator.currentTime, integrator.lastDt
    integrator.initializeDerivatives(t, dt, state, derivs)
    derivs.Zero()
    integrator.evaluateDerivatives(t, dt, body.db, state, derivs)
    integrator.finalizeDerivatives(t, dt, body.db, state, derivs)


def gather_by_id(nodes, fields, dim=None):
    """All ranks' internal values of `fields` (name -> Spheral field of `nodes`, including "pid") to rank 0, as
    numpy arrays sorted by pid; None on the other ranks. Collective. dim defaults to the fields' own (the "2d" or
    "3d" ending of their type name; RZ fields are 2d)."""
    if "pid" not in fields:
        raise ValueError('fields must include "pid"')
    if dim is None:
        dim = int(type(fields["pid"]).__name__[-2])
    comm = _comm()
    nloc = nodes.numInternalNodes
    local = {}
    for k, f in fields.items():
        arr = field_to_numpy(f, dim)
        if len(arr) != nloc:
            raise ValueError(f"{k}: {len(arr)} values for {nloc} internal nodes")
        local[k] = arr
    parts = comm.gather(local, root=0)
    if comm.rank != 0:
        return None
    out = {k: np.concatenate([p[k] for p in parts], axis=0) for k in fields}
    order = np.argsort(out["pid"], kind="stable")
    return {k: v[order] for k, v in out.items()}
