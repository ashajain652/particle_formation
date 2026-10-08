"""Diagnostic driver (2026-10-07, not part of the prototype): run the CLI with the MeltingBody methods that move film
or deep liquid wrapped; the first one after which m_f or m_d is not finite (or Runoff.transport returning a
non-finite result) dumps its inputs to <dump.npz> and prints which it was.

    cd <tree> && python ../rebuild/diag/nan_trap.py <dump.npz> -- <cli run args>
"""
import os
import sys
import warnings

sys.path.insert(0, os.getcwd())
import numpy as np

out = sys.argv[1]
rest = sys.argv[3:] if sys.argv[2] == "--" else sys.argv[2:]
from reentry_model import body as body_mod, cli, film  # noqa: E402

state = {"dumped": False}


def dump(tag, **arrays):
    if state["dumped"]:
        return
    state["dumped"] = True
    np.savez(out, tag=tag, **{k: np.asarray(v) for k, v in arrays.items() if v is not None and not callable(v)})
    print("TRAP:", tag, flush=True)


_transport = film.Runoff.transport


def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
    m0 = np.array(m_f, dtype=float)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = _transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps)
    if (w or not np.all(np.isfinite(res[0]))) and not state["dumped"]:
        # replay the sub-steps to find the one that failed, and keep its matrix
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        m, dt_s = m0.copy(), dt / substeps
        info = {}
        for n in range(1, substeps + 1):
            b = m / (rho_l * areas)
            q = q_of_thickness(b)
            c_ij, c_ji = self.edge_coefficients(q, b, t_hat, areas)
            outc = np.bincount(self.i, c_ij, self.n_patches) + np.bincount(self.j, c_ji, self.n_patches)
            rows = np.concatenate([np.arange(self.n_patches), self.j, self.i])
            cols = np.concatenate([np.arange(self.n_patches), self.i, self.j])
            vals = np.concatenate([1.0 + dt_s * outc, -dt_s * c_ij, -dt_s * c_ji])
            A = sp.csc_matrix((vals, (rows, cols)), shape=(self.n_patches, self.n_patches))
            with warnings.catch_warnings(record=True) as w2:
                warnings.simplefilter("always")
                m_new = np.maximum(spla.spsolve(A, m), 0.0)
            if w2 or not np.all(np.isfinite(m_new)):
                info = dict(substep=n, m=m, b=b, q=q, c_ij=c_ij, c_ji=c_ji, out=outc, rows=rows, cols=cols, vals=vals,
                            m_new=m_new, warn=" | ".join(str(x.message) for x in w2))
                break
            m = m_new
        dump("transport", m0=m0, t_hat=t_hat, areas=areas, edge_i=self.i, edge_j=self.j, length=self.length,
             n_i=self.n_i, n_j=self.n_j, warn=" | ".join(str(x.message) for x in w), **info)
    return res


film.Runoff.transport = transport


def wrap(name):
    f = getattr(body_mod.MeltingBody, name)

    def g(self, *a, **k):
        r = f(self, *a, **k)
        if not state["dumped"] and (not np.all(np.isfinite(self.m_f)) or not np.all(np.isfinite(self.m_d))):
            dump("after " + name, m_f=self.m_f, m_d=self.m_d)
        return r
    setattr(body_mod.MeltingBody, name, g)


for name in ("_deep_runoff", "_film_and_spray", "_kill", "_add_to_film", "_freeze_back", "_feed_exposed", "_hand_over"):
    wrap(name)

_deep = body_mod.MeltingBody._deep_runoff


def deep(self, dt, flow, delta_m, reach, T, h_node, h_e, h_p):
    self._trap_t = getattr(self, "_trap_t", 0) + 1
    return _deep(self, dt, flow, delta_m, reach, T, h_node, h_e, h_p)


body_mod.MeltingBody._deep_runoff = deep
sys.exit(cli.main(rest))
