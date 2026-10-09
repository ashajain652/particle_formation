"""spheral_frag.runner.packages: the replay's Python physics packages (Spheral M2 plan, Task 4; spec §7, §9.1).

  temperature.ImposedTemperature   u = h(T) with T the body stage's temperature table, linear in time (D1, D2)
  damping.Damping                  a = -2 zeta omega_1 v, mass-proportional (D5)
  loads.SurfaceLoads               surface particles, normals, exposure, inclination, area shares and the
                                   aerodynamic pressure and shear from the prepared load tables (spec §7.2)
  loads.Deceleration               the frame's deceleration sum F / sum m on every particle (spec §7.1)

Each is a Python subclass of Spheral's `Physics` (plan fact 5: a package adds to
`derivs.vectorFields(HydroFieldNames.hydroAcceleration)` in `evaluateDerivatives`, and `dt()` returns a tuple),
built for one Spheral module by a `make_*` function (one form per process, as eos.py). They are appended to the
integrator after the hydro package.

**Built per interval.** A frame boundary is a restart (nodes.py: rebuild, new packages, `prime_derivatives`), so every
package here is built from the current body and plain arrays at the start of each interval and holds nothing that
must survive it; what it records (the deceleration, the temperature monitor) is read out by the driver before the
packages are dropped. Each package caches the local particles' rows at construction and refuses to run if the local
node count changes under it (a redistribution or deletion inside an interval is not supported; the driver deletes
at the boundary, before it rebuilds).

**Time** is Spheral's time, which the replay keeps equal to the flight time (mass scaling makes sound slow enough that
0.5 s of Spheral time is 0.5 s of flight): the driver sets the integrator's `currentTime` to the frame's
`frames_time_s` at the start of the interval.

Fields move between Spheral and numpy through `Field.serialize`/`deserialize` (internal values only), which costs
about as much as a scalar field's `internalValues()` (plan fact 6) also for vectors, instead of 42 ms per vector
field through Python `Vector`s at 49,000 particles."""
from __future__ import annotations

import struct

import numpy as np

F64 = np.dtype("<f8")


def get_vector(field, dim):
    """The internal values of a VectorField as an (n, dim) float64 array (a copy)."""
    b = field.serialize()
    n = field.numInternalElements
    off = 16 + len(field.name.encode())
    return np.frombuffer(b, dtype=F64, offset=off, count=n*dim).reshape(n, dim).copy()


def put_vector(field, values):
    """Write the internal values of a VectorField from an (n, dim) array (ghosts untouched)."""
    name = field.name.encode()
    v = np.ascontiguousarray(values, dtype=F64)
    if len(v) != field.numInternalElements:
        raise ValueError(f"{field.name}: {len(v)} rows for {field.numInternalElements} internal nodes")
    field.deserialize(struct.pack("<Q", len(name)) + name + struct.pack("<Q", len(v)) + v.tobytes())


def add_acceleration(S, derivs, dim, a):
    """Add the (n_local, dim) array `a` to the hydro acceleration of the derivatives (the one node list's field)."""
    fl = derivs.vectorFields(S.HydroFieldNames.hydroAcceleration)
    fields = list(fl)
    if len(fields) != 1:
        raise RuntimeError(f"expected one node list's acceleration field, got {len(fields)}")
    f = fields[0]
    put_vector(f, get_vector(f, dim) + a)


def state_velocity(S, state, dim):
    """The (n_local, dim) velocities of the state (the one node list's field)."""
    fields = list(state.vectorFields(S.HydroFieldNames.velocity))
    if len(fields) != 1:
        raise RuntimeError(f"expected one node list's velocity field, got {len(fields)}")
    return get_vector(fields[0], dim)


def comm():
    from mpi4py import MPI
    return MPI.COMM_WORLD


def allsum(x):
    """The sum over every rank of a float or numpy array (float64)."""
    from mpi4py import MPI
    a = np.asarray(x, dtype=np.float64)
    out = np.empty_like(a)
    comm().Allreduce(np.ascontiguousarray(a), out, op=MPI.SUM)
    return out if out.shape else float(out)


def allmax(x):
    from mpi4py import MPI
    return comm().allreduce(float(x), op=MPI.MAX)


def allmin(x):
    from mpi4py import MPI
    return comm().allreduce(float(x), op=MPI.MIN)


class LocalCount:
    """Guards a package's cached local rows: the node list's internal count must not change within an interval."""

    def __init__(self, nodes):
        self.nodes = nodes
        self.n = nodes.numInternalNodes

    def check(self, who):
        n = self.nodes.numInternalNodes
        if n != self.n:
            raise RuntimeError(f"{who}: the local node count changed from {self.n} to {n} within an interval; "
                               "rebuild the packages after deleting or redistributing particles")
