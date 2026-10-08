"""Diagnostic driver (2026-10-07, not part of the prototype): run the CLI with scipy's spsolve wrapped; the first
solve whose result is not finite dumps the matrix, the right-hand side and the locals of Runoff.transport (and which
body method called it) to <dump.npz>, then lets the run continue to its failure.

    cd prototype/proto3 && python ../rebuild/diag/singular_trap.py <dump.npz> -- <cli run args>
"""
import inspect
import os
import sys
import warnings

sys.path.insert(0, os.getcwd())
import numpy as np
import scipy.sparse.linalg as spla

out = sys.argv[1]
rest = sys.argv[3:] if sys.argv[2] == "--" else sys.argv[2:]
_spsolve = spla.spsolve
state = {"dumped": False, "calls": 0}


def spsolve(A, b, *a, **k):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        x = _spsolve(A, b, *a, **k)
    state["calls"] += 1
    if (w or not np.all(np.isfinite(x))) and not state["dumped"]:
        state["dumped"] = True
        fr = inspect.currentframe().f_back
        loc = fr.f_locals
        caller = fr.f_back.f_code.co_name if fr.f_back else "?"
        A = A.tocoo()
        d = {"rows": A.row, "cols": A.col, "vals": A.data, "b": np.asarray(b), "x": np.asarray(x), "caller": caller,
             "function": fr.f_code.co_name, "warnings": " | ".join(str(i.message) for i in w)}
        for key in ("c_ij", "c_ji", "out", "m", "areas", "dt_s", "n", "t_hat"):
            if key in loc:
                d[key] = np.asarray(loc[key])
        if "self" in loc and hasattr(loc["self"], "i"):
            d["edge_i"], d["edge_j"], d["length"] = loc["self"].i, loc["self"].j, loc["self"].length
        bf = fr.f_back.f_locals if fr.f_back else {}
        body = bf.get("self")
        for key in ("D", "b", "G", "girin"):
            if key in bf:
                d["caller_" + key] = np.asarray(bf[key])
        if body is not None and hasattr(body, "m_f"):
            d["m_f"], d["m_d"], d["windward"] = body.m_f, body.m_d, body.windward
        np.savez(out, **d)
        print("TRAP: non-finite / warned spsolve in {} called from {}: {}".format(fr.f_code.co_name, caller, d["warnings"]), flush=True)
    return x


spla.spsolve = spsolve
from reentry_model import cli  # noqa: E402

sys.exit(cli.main(rest))
