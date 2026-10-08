"""Diagnostic driver (2026-10-08, not part of the prototype): run the CLI and report, at the end, how many donor rows
of the runoff solves film.EMPTYING_MAX held (post/05), split by caller (the film's runoff, the deep runoff), and the
largest emptying number dt_s * out seen. Read-only: it registers each Runoff and wraps edge_coefficients' caller.

    cd <tree> && python ../rebuild/diag/capped_count.py <summary.json> -- <cli run args>
"""
import json
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np

out = sys.argv[1]
rest = sys.argv[3:] if sys.argv[2] == "--" else sys.argv[2:]
from reentry_model import body as body_mod, cli, film  # noqa: E402

stats = {"film": {"capped_rows": 0, "solves": 0, "max_emptying": 0.0, "solves_over_1e8": 0},
         "deep": {"capped_rows": 0, "solves": 0, "max_emptying": 0.0, "solves_over_1e8": 0}}
current = {"who": "film"}
_transport = film.Runoff.transport
_coef = film.Runoff.edge_coefficients


def edge_coefficients(self, q, b, t_hat, areas):
    c_ij, c_ji = _coef(self, q, b, t_hat, areas)
    o = np.bincount(self.i, c_ij, self.n_patches) + np.bincount(self.j, c_ji, self.n_patches)
    dt_s = current.get("dt_s", 0.125)
    e = float(dt_s * o.max()) if o.size else 0.0
    st = stats[current["who"]]
    st["solves"] += 1
    st["max_emptying"] = max(st["max_emptying"], e)
    st["solves_over_1e8"] += int(e > 1e8)
    return c_ij, c_ji


def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
    current["dt_s"] = dt / substeps
    before = self.capped_rows
    r = _transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps)
    stats[current["who"]]["capped_rows"] += self.capped_rows - before
    return r


_deep = body_mod.MeltingBody._deep_runoff


def deep(self, *a, **k):
    current["who"] = "deep"
    try:
        return _deep(self, *a, **k)
    finally:
        current["who"] = "film"


film.Runoff.edge_coefficients = edge_coefficients
film.Runoff.transport = transport
body_mod.MeltingBody._deep_runoff = deep
try:
    code = cli.main(rest)
finally:
    json.dump(stats, open(out, "w"), indent=1)
sys.exit(code)
