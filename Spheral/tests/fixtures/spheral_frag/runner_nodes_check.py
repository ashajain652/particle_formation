"""Spheral-side driver of tests/test_spheral_frag_runner_nodes.py (Spheral M2 plan, Task 1). Run through the launcher:

    spheral_frag/container/spheral [-n N] tests/fixtures/spheral_frag/runner_nodes_check.py MODE --form F ...

Modes (every output a plain .npz written by rank 0, read by the drama_env test):

  build    --input IN [--delete DEL] --out OUT
           make_body from IN's arrays; gather by id ("built.<key>"); the per-rank local counts and pid-field lengths;
           with DEL (an array "pid"), each rank deletes its local particles among them, then gather again
           ("deleted.<key>").
  run      (--input IN | --checkpoint CK) --steps N [--checkpoint-at M --checkpoint-out CKOUT] --out OUT
           build (from IN at t = 0, or from CK at its time, step count and last dt), SolidSPH + CheapSynchronousRK2 +
           SpheralController without the initial H iteration, N steps (each as bench.py's timed loop: integrator
           step, step count, periodic work), a checkpoint after step M when asked, then gather by id ("final.<key>")
           with time, dt and steps. --boundary restart (default) makes step M a frame-boundary restart in the
           uninterrupted run as in the resumed one: rebuilt from the gathered state and primed
           (nodes.rebuild, nodes.prime_derivatives); "prime" and "none" measure what is lost without it.
  rebuild  --checkpoint CK --out OUT
           make_body from a checkpoint (any process count), gather by id ("built.<key>").

The material is a test material (the thermal classes come in Task 3): Murnaghan with bench.py's AA7075 parameters
and ConstantStrength (G 27.6 GPa, Y 20 MPa), so that the deviatoric stress grows and the plastic strain becomes
non-zero within the steps under the imposed expansion the input's velocities carry."""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from spheral_frag.runner import checkpoint, nodes as rn  # noqa: E402

RHO0, ETAMIN, ETAMAX = 2813.0, 0.2, 4.0
K = 71.7e9/(3.0*(1.0 - 2.0*0.33))
G, Y = 27.6e9, 2.0e7
CFL = 0.25


def material(S):
    eos = S.MurnaghanEquationOfState(RHO0, ETAMIN, ETAMAX, 1.0, K, 26.98, S.MKS())
    strength = S.ConstantStrength(G, Y, eos)
    return eos, strength


def load_arrays(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files if k != "dx"}, float(z["dx"])


def build(S, form, arrays, dx, name="body"):
    eos, strength = material(S)
    return rn.make_body(S, form, arrays, eos, strength, dx=dx, rho_range=(ETAMIN*RHO0, ETAMAX*RHO0), name=name)


def setup(S, body, dx, t0, steps0, dt0, cycle0):
    from SpheralController import SpheralController
    hydro = S.SPH(dataBase=body.db, W=body.kernel, cfl=CFL, compatibleEnergyEvolution=True,
                  densityUpdate=S.IntegrateDensity, XSPH=False, ASPH=False)
    integrator = S.CheapSynchronousRK2Integrator(body.db)
    integrator.appendPhysicsPackage(hydro)
    integrator.dtMin = 1e-15
    integrator.dtMax = 1.0
    integrator.verbose = False
    controller = SpheralController(integrator, body.kernel, statsStep=10**9, printStep=10**9,
                                   redistributeStep=None, restartStep=None, vizStep=None, vizTime=None, SPH=True,
                                   iterateInitialH=False, initialTime=t0)
    integrator.lastDt = 0.1*CFL*dx/np.sqrt((K + 4.0/3.0*G)/RHO0) if dt0 is None else dt0
    integrator.currentCycle = cycle0
    controller.totalSteps = steps0
    return integrator, controller, hydro


def save(path, prefix_dicts, **scalars):
    out = {}
    for prefix, d in prefix_dicts.items():
        for k, v in d.items():
            out[f"{prefix}.{k}"] = v
    out.update({k: np.asarray(v) for k, v in scalars.items()})
    np.savez(path, **out)


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("build", "run", "rebuild"))
    p.add_argument("--form", choices=rn.FORMS)
    p.add_argument("--input")
    p.add_argument("--checkpoint")
    p.add_argument("--delete")
    p.add_argument("--steps", type=int, default=0)
    p.add_argument("--checkpoint-at", type=int, default=None)
    p.add_argument("--checkpoint-out")
    p.add_argument("--boundary", choices=("none", "prime", "restart"), default="restart",
                   help="at the checkpoint (and on resuming): restart = rebuild from the gathered state and prime "
                        "the derivatives (the runner's frame boundary); prime = prime in place; none = nothing")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    ck = checkpoint.read(a.checkpoint) if a.checkpoint else None
    form = ck["form"] if ck else a.form
    S, _ = rn.spheral_module(form)
    if ck:
        arrays, dx = checkpoint.arrays(ck), ck["dx"]
    else:
        arrays, dx = load_arrays(a.input)
    body = build(S, form, arrays, dx)

    if a.mode in ("build", "rebuild"):
        nloc = body.nodes.numInternalNodes
        counts = comm.gather((nloc, len(body.pid.internalValues())), root=0)
        built = body.gather()
        deleted = None
        if a.delete:
            with np.load(a.delete) as z:
                body.delete_by_id(z["pid"])
            counts_after = comm.gather((body.nodes.numInternalNodes, len(body.pid.internalValues())), root=0)
            deleted = body.gather()
        if comm.rank == 0:
            d = {"built": built}
            sc = dict(counts=np.array(counts), processes=comm.size)
            if deleted is not None:
                d["deleted"] = deleted
                sc["counts_after"] = np.array(counts_after)
            save(a.out, d, **sc)
        return 0

    # ---------------------------------------------------------------- run
    t0 = float(ck["extra"]["time"]) if ck else 0.0
    steps0 = int(ck["extra"]["steps"]) if ck else 0
    dt0 = float(ck["extra"]["dt"]) if ck else None
    cycle0 = int(ck["extra"]["cycle"]) if ck else 0
    integrator, controller, keep = setup(S, body, dx, t0, steps0, dt0, cycle0)
    if ck and a.boundary in ("prime", "restart"):
        rn.prime_derivatives(body, integrator)
    for i in range(a.steps):
        integrator.step(1e40)
        controller.totalSteps += 1
        controller.doPeriodicWork()
        if a.checkpoint_at is not None and i + 1 == a.checkpoint_at:
            extra = dict(time=controller.time(), dt=integrator.lastDt, steps=controller.totalSteps,
                         cycle=integrator.currentCycle)
            checkpoint.write(a.checkpoint_out, body, run="runner_nodes_check", k=steps0 + i + 1, scale=1.0,
                             extra=extra)
            if a.boundary == "restart":
                del integrator, controller, keep
                body = rn.rebuild(body, name=f"body{i + 1}")
                integrator, controller, keep = setup(S, body, dx, extra["time"], extra["steps"], extra["dt"],
                                                     extra["cycle"])
            if a.boundary in ("prime", "restart"):
                rn.prime_derivatives(body, integrator)
    final = body.gather()
    if comm.rank == 0:
        save(a.out, {"final": final}, time=controller.time(), dt=integrator.lastDt, steps=controller.totalSteps,
             cycle=integrator.currentCycle, processes=comm.size)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
