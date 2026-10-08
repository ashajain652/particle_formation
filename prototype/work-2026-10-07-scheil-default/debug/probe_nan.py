"""Throwaway probe for the 0.00625 s crash of 2026-10-07 ("film mass must be one finite non-negative value per node").
Re-runs the dt00625 series run from the frozen copy with the same step switch as harness/dtswitch.py, and after every
stage of every melt step from the switch on checks the film and deep accounts (finite, non-negative), the non-rigid
depth, film.lubrication's outputs and the runoff transport's result. At the first bad value it saves the arrays and the
stage's inputs to bad.npz, prints where, and stops. No package file is edited."""
import json, os, sys, traceback
copy_dir, out = sys.argv[1], sys.argv[2]
KN_SWITCH, DT_FINE = 0.01, float(sys.argv[3])
sys.path.insert(0, copy_dir)
import numpy as np
from reentry_model import body as B, cli, coupled, film

S = {"run": None, "switched": False, "t": None, "step": 0, "log": []}
orig_init = coupled.CoupledRun.__init__
def init(self, *a, **k):
    orig_init(self, *a, **k); S["run"] = self
coupled.CoupledRun.__init__ = init

def bad(x):
    x = np.asarray(x, dtype=float)
    return int((~np.isfinite(x)).sum()), int((x < 0).sum())

def check(b, stage, **extra):
    if not S["switched"]:
        return
    nf, nneg = bad(b.m_f); df, dneg = bad(b.m_d)
    items = {k: bad(v) for k, v in extra.items() if v is not None}
    worst = nf or nneg or df or dneg or any(a or c for a, c in items.values())
    S["log"].append((S["t"], stage, nf, nneg, df, dneg, {k: v for k, v in items.items()}))
    if worst:
        print("BAD at t=%.5f step %d stage %s: m_f nonfinite %d negative %d; m_d nonfinite %d negative %d; %s" % (
            S["t"], S["step"], stage, nf, nneg, df, dneg, items), flush=True)
        save = {"m_f": b.m_f, "m_d": b.m_d, "phi": b.phi, "T": b.solver.temperature(), "areas": b.surface.areas,
                "theta": b.theta, "windward": b.windward}
        for k, v in extra.items():
            if v is not None:
                save["x_" + k] = np.asarray(v, dtype=float)
        np.savez(os.path.join(out, "bad.npz"), **save)
        with open(os.path.join(out, "log.json"), "w") as fh:
            json.dump([list(map(str, r)) for r in S["log"][-200:]], fh)
        traceback.print_stack(limit=6)
        os._exit(7)

orig_adv = B.MeltingBody.advance
def advance(self, t, dt, loads, state=None):
    S["t"] = t; S["step"] += 1
    check(self, "before conduction")
    orig_adv(self, t, dt, loads, state)
    check(self, "end of step")
    run = S["run"]
    if not S["switched"] and state is not None and run is not None and np.isfinite(state.kn) and state.kn < KN_SWITCH:
        run.settings.dt = DT_FINE; S["switched"] = True
        print("switch at", t, flush=True)
B.MeltingBody.advance = advance

def wrap(name):
    orig = getattr(B.MeltingBody, name)
    def f(self, *a, **k):
        r = orig(self, *a, **k)
        check(self, "after " + name)
        return r
    setattr(B.MeltingBody, name, f)
for name in ("_add_to_film", "_deep_runoff", "_film_and_spray", "_freeze_back", "_kill", "_feed_exposed"):
    if hasattr(B.MeltingBody, name):
        wrap(name)

orig_nr = B.MeltingBody.nonrigid_depth
def nonrigid(self, *a, **k):
    d = orig_nr(self, *a, **k)
    if S["switched"] and (~np.isfinite(d)).any():
        print("nonrigid_depth non-finite at t=%.5f: %d of %d (inf %d)" % (S["t"], int((~np.isfinite(d)).sum()), d.size, int(np.isinf(d).sum())), flush=True)
    return d
B.MeltingBody.nonrigid_depth = nonrigid

orig_lub = film.lubrication
def lubrication(tau, G, b, delta_m, mu_l, b_layer=None):
    o = orig_lub(tau, G, b, delta_m, mu_l, b_layer)
    if S["switched"]:
        for i, nm in ((0, "V"), (1, "q")):
            v = np.asarray(o[i], dtype=float)
            if (~np.isfinite(v)).any():
                bl = np.asarray(b_layer, dtype=float) if b_layer is not None else None
                print("lubrication %s non-finite at t=%.5f: %d; b nonfinite %d; b_layer nonfinite %d; tau nonfinite %d; G nonfinite %d" % (
                    nm, S["t"], int((~np.isfinite(v)).sum()), bad(b)[0], bad(bl)[0] if bl is not None else -1, bad(tau)[0], bad(G)[0]), flush=True)
    return o
film.lubrication = lubrication

orig_tr = film.Runoff.transport
def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
    m_in = np.array(m_f, dtype=float)
    o = orig_tr(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps)
    if S["switched"] and (bad(o[0])[0] or bad(o[0])[1]):
        print("transport output bad at t=%.5f: in nonfinite %d negative %d; out nonfinite %d negative %d" % (
            S["t"], *bad(m_in), *bad(o[0])), flush=True)
    return o
film.Runoff.transport = transport

argv = ["run", "--diameter", "100", "--altitude", "77.500133", "--velocity", "7.5", "--flight-path-angle", "-0.959331",
        "--atmosphere", "us76", "--thermal", "fem", "--melt", "on", "--heating", "physics", "--removal", "girin",
        "--outdir", out, "--t-max", "120", "--quiet"] + sys.argv[4:]
code = cli.main(argv)
print("exit", code, "at t", S["t"], flush=True)
