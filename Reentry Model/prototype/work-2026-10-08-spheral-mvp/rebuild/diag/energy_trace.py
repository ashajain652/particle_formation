"""Diagnostic driver (2026-10-07, not part of the prototype): run `python -m reentry_model run ...` from the current
directory's reentry_model with MeltingBody.advance wrapped so that every macro step records the energy-balance
numerator (E - E0 - absorbed + removed + pending + dropped) before the conduction solve, after it, and after
melt_step, plus Newton iterations, the solve's own balance and the nodes on the steep parts of h(T).

    cd prototype/proto3 && python ../rebuild/diag/energy_trace.py <trace.npz> [--newton-tol X] -- <cli run args>
"""
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np

from reentry_model import body as body_mod, cli
from reentry_model.thermal import skfem_backend

out = sys.argv[1]
rest = sys.argv[2:]
tol = None
if rest and rest[0] == "--newton-tol":
    tol = float(rest[1]); rest = rest[2:]
if rest and rest[0] == "--":
    rest = rest[1:]

if tol is not None:
    _init = skfem_backend.SkfemThermalSolver.__init__

    def init(self, *a, **k):
        _init(self, *a, **k)
        self.newton_tol = tol
    skfem_backend.SkfemThermalSolver.__init__ = init

DECOMP = None
if rest and rest[0] == "--decomp":                 # the final iterate's energy error split (step_decomp.py)
    rest = rest[1:]
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import step_decomp
    skfem_backend.SkfemThermalSolver.step = step_decomp.step
    DECOMP = step_decomp.DECOMP
if rest and rest[0] == "--":
    rest = rest[1:]
rec = {k: [] for k in ("t", "num0", "num1", "num2", "absorbed", "iters", "dE_solve", "q_net_dt", "n_eut", "n_liqramp",
                       "cpeff_max", "feed", "cascade", "passes", "n_dead", "frozen", "released", "film", "pend", "dropped")}


def numerator(b):
    return b.energy() - b.energy0 - b.absorbed_heat + b.removed_enthalpy + b.pending_load.sum() + b.total_dropped_load


_adv = body_mod.MeltingBody.advance


def advance(self, t, dt, loads, state=None):
    mat = self.material
    n0 = numerator(self)
    orig_melt = self.melt_step
    box = {}

    def melt_step(t_, dt_, st_):
        box["n1"] = numerator(self)
        T = self.solver.temperature()
        act = np.zeros(len(T), bool)
        act[self.mesh.tets[self.mesh.active].ravel()] = True
        box["n_eut"] = int(((T > mat.T_solidus) & (T < mat.T_solidus + 4.0) & act).sum())
        box["n_liqramp"] = int(((T > mat.T_liquidus - 15.0) & (T < mat.T_liquidus) & act).sum())
        box["cpeff_max"] = float(mat.cp_eff(T[act]).max())
        return orig_melt(t_, dt_, st_)
    self.melt_step = melt_step
    try:
        _adv(self, t, dt, loads, state)
    finally:
        del self.melt_step
    n2 = numerator(self)
    r = self.last
    lm = self.last_melt
    rec["t"].append(t); rec["num0"].append(n0); rec["num1"].append(box["n1"]); rec["num2"].append(n2)
    rec["absorbed"].append(self.absorbed_heat); rec["iters"].append(r.iterations)
    rec["q_net_dt"].append((r.Q_conv - r.Q_rad) * dt)
    rec["dE_solve"].append(box["n1"] - n0)
    rec["n_eut"].append(box["n_eut"]); rec["n_liqramp"].append(box["n_liqramp"]); rec["cpeff_max"].append(box["cpeff_max"])
    rec["feed"].append(lm.get("feed_mass", 0.0)); rec["cascade"].append(lm.get("cascade_mass", 0.0))
    rec["passes"].append(lm.get("cascade_passes", 0)); rec["n_dead"].append(lm.get("n_dead", 0))
    rec["frozen"].append(lm.get("frozen_mass", 0.0)); rec["released"].append(lm.get("released_mass", 0.0))
    rec["film"].append(float(self.m_f.sum())); rec["pend"].append(float(self.pending_load.sum()))
    rec["dropped"].append(self.total_dropped_load)
    if len(rec["t"]) % 20 == 0:
        save()


def save():
    extra = {"decomp": np.array(DECOMP, dtype=float)} if DECOMP is not None else {}
    np.savez(out, **{k: np.array(v, dtype=float) for k, v in rec.items()}, **extra)


body_mod.MeltingBody.advance = advance
try:
    code = cli.main(rest)
finally:
    save()
sys.exit(code)
