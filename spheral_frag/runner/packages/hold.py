"""Hold: particles held at rest exactly (Spheral M2 plan, Task 7, check 5: the bar held at both ends).

The held particles' rows of the hydro acceleration are set to zero in `evaluateDerivatives`, after every other
package has added to it (the package is appended last). CheapSynchronousRK2 advances every step twice from the
derivatives: the half step from the previous midpoint's (or, at the start of an interval, `prime_derivatives`'),
the full step from this midpoint's; with a held particle's acceleration zero in both, its velocity stays at its start
value (zero) through both stages, and its position derivative -- the velocity, with XSPH off -- is zero at both, so
neither moves by a single bit.

Why not Spheral's `ConstantVelocityBoundary` (the plan's suggestion, plan fact 11): it resets the velocity field in
`enforceBoundary`, which CheapSynchronousRK2 calls once, at the end of the step (CheapSynchronousRK2.cc), not after
the midpoint update. The midpoint velocity of a held particle is therefore v_mid = a dt/2, and the full step moves it
by a dt^2/2 with its own (reaction) acceleration a: measured on the RZ 2.2 mm bar (2026-10-10, check 5's
development runs, 20 crossings), the held particles had moved by up to 17.6 micrometres (0.008 dx) after the first
interval and 82 micrometres after the third -- not a hold. Its velocity is zero at every step end all the same, so
the drift shows only in the positions. (That bar's ringing turned out to be the RZ 2.2 mm rod's own instability,
present with this package too; check 5's docstring.)

Hold(body, held_local_rows). `dt()` votes no limit. The rows are cached at construction (LocalCount guards them)."""
from __future__ import annotations

import numpy as np

from spheral_frag.runner.packages import LocalCount, get_vector, put_vector

_CLASSES = {}


def hold_class(S):
    if S in _CLASSES:
        return _CLASSES[S]

    class Hold(S.Physics):
        """Zero hydro acceleration on the held local rows (module docstring)."""

        def __init__(self, body, rows):
            S.Physics.__init__(self)
            self.body = body
            self.dim = body.dim
            self.rows = np.asarray(rows, dtype=np.int64)
            self.count = LocalCount(body.nodes)

        def evaluateDerivatives(self, t, dt, db, state, derivs):
            if not len(self.rows):
                return
            self.count.check("Hold")
            fields = list(derivs.vectorFields(S.HydroFieldNames.hydroAcceleration))
            if len(fields) != 1:
                raise RuntimeError(f"expected one node list's acceleration field, got {len(fields)}")
            a = get_vector(fields[0], self.dim)
            a[self.rows] = 0.0
            put_vector(fields[0], a)

        def dt(self, db, state, derivs, t):
            return (1.0e100, "Hold: no limit")

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "Hold"

    _CLASSES[S] = Hold
    return Hold


def make_hold(S, body, rows):
    """A Hold package for the Spheral module S on the body's local rows `rows`."""
    return hold_class(S)(body, rows)
