"""Spheral-side driver of tests/test_spheral_frag_runner_material.py (Spheral M2 plan, Task 3). Run through the
launcher from the repository:

    spheral_frag/container/spheral [-n N] tests/fixtures/spheral_frag/runner_material_check.py MODE --form F ...

Modes (outputs are plain .npz written by rank 0 and read by the drama_env test):

  methods  --table TAB --input IN --scale S --out OUT
           IN holds synthetic states (T, u, rho at scale 1, plastic strain rate). Every method of PyThermalEOS and
           PyThermalStrength (laws flight, elastic, constant:2e7) is called on Spheral fields of a node list with
           that many internal nodes, at scale 1 with rho and at scale S with S rho, and once more on an EOS built at
           1 whose `scale` is then set to S; also setPressure at rho (1 +- 1e-6) for the finite-difference dP/drho.
  rest     --table TAB --input IN --scale S --T T --steps N --out OUT
           A body from IN's lattice (x_ref, V, pid) at uniform T: rho = S rho_free(T), m = rho V, u = h(T), at rest;
           SolidSPH (S.SPH on a solid node list) + CheapSynchronousRK2 for N steps; gather P, S, v, rho by pid.
  bench    --table TAB --body DIR --material {thermal,builtin} --steps N --warmup W [--max-load L] --out OUT
           bench.py's timed loop on a body stage's body (positions, masses, densities, T0), with bench.py's uniform
           expansion (10 1/s): PyThermalEOS + PyThermalStrength (flight law), or Murnaghan (n = 1, K(300 K)) +
           ConstantStrength (G(300 K), 20 MPa). Per-step wall times (barrier to barrier) to OUT."""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from spheral_frag import material  # noqa: E402
from spheral_frag.runner import nodes as rn  # noqa: E402
from spheral_frag.runner.eos import get_internal, make_eos, put_internal  # noqa: E402
from spheral_frag.runner.strength import make_strength  # noqa: E402

CFL = 0.25
EDOT = 10.0          # 1/s, bench.py's imposed expansion
LAWS = ("flight", "elastic", "constant:2e7")


def load(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def methods(S, a, table):
    d = load(a.input)
    n = len(d["T"])
    dim = rn.DIM[a.form]
    keep = []

    def nodelist(eos, strength, label):
        box = (S.Vector(*([-1.0]*dim)), S.Vector(*([1.0]*dim)))
        nl = S.makeSolidNodeList(label, eos, strength, xmin=box[0], xmax=box[1])
        nl.numInternalNodes = n
        keep.extend([eos, strength, nl])
        return nl

    def scalar(nl, name, values):
        f = S.ScalarField(name, nl)
        put_internal(f, values)
        keep.append(f)
        return f

    out = {}
    s = a.scale
    eos1 = make_eos(S, table, 1.0)
    eoss = make_eos(S, table, s)
    eoset = make_eos(S, table, 1.0)
    eoset.scale = s
    out["refdens"] = np.array([eos1.referenceDensity, eoss.referenceDensity, eoset.referenceDensity])
    out["valid"] = np.array([eos1.valid(), eoss.valid(), eoset.valid()])
    for tag, eos, scale in (("s1", eos1, 1.0), ("s", eoss, s), ("set", eoset, s)):
        strengths = {law: make_strength(S, eos, law) for law in LAWS}
        nl = nodelist(eos, strengths["flight"], "m_" + tag)
        rho = scalar(nl, "rho", d["rho"]*scale)
        u = scalar(nl, "u", d["u"])
        T = scalar(nl, "T", d["T"])
        P = scalar(nl, "P", np.zeros(n))
        psr = scalar(nl, "psr", d["psr"])
        ps = scalar(nl, "ps", np.zeros(n))
        D = S.SymTensorField("D", nl)
        keep.append(D)

        def call(m, *args):
            f = S.ScalarField("out", nl)
            getattr(eos, m)(f, *args)
            return get_internal(f).copy()

        for m, args in (("setPressure", (rho, u)), ("setTemperature", (rho, u)),
                        ("setSpecificThermalEnergy", (rho, T)), ("setSpecificHeat", (rho, T)),
                        ("setSoundSpeed", (rho, u)), ("setGammaField", (rho, u)), ("setBulkModulus", (rho, u)),
                        ("setEntropy", (rho, u))):
            out[f"{tag}.{m}"] = call(m, *args)
        fs = [S.ScalarField(f"pd{i}", nl) for i in range(3)]
        eos.setPressureAndDerivs(*fs, rho, u)
        for i, k in enumerate(("P", "dPdu", "dPdrho")):
            out[f"{tag}.setPressureAndDerivs.{k}"] = get_internal(fs[i]).copy()
        for sign, lbl in ((1.0, "plus"), (-1.0, "minus")):
            r2 = scalar(nl, "rho_" + lbl, d["rho"]*scale*(1.0 + sign*1e-6))
            out[f"{tag}.P_{lbl}"] = call("setPressure", r2, u)
        for law, st in strengths.items():
            G = S.ScalarField("G", nl)
            st.shearModulus(G, rho, u, P, D)
            Y = S.ScalarField("Y", nl)
            st.yieldStrength(Y, rho, u, P, ps, psr, D)
            out[f"{tag}.shearModulus.{law}"] = get_internal(G).copy()
            out[f"{tag}.yieldStrength.{law}"] = get_internal(Y).copy()
            out[f"{tag}.provides.{law}"] = np.array([st.providesSoundSpeed(), st.providesBulkModulus()])
    np.savez(a.out, **out)
    return 0


def setup(S, body, dx):
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
                                   iterateInitialH=False)
    return hydro, integrator, controller


def rest(S, a, table):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    d = load(a.input)
    dx = float(d["dx"])
    T = np.full(len(d["pid"]), a.T)
    rho = a.scale*table.rho_free(T)
    arrays = dict(x_ref=d["x_ref"], pid=d["pid"], rho=rho, m=rho*d["V"], u=table.enthalpy(T))
    eos = make_eos(S, table, a.scale)
    strength = make_strength(S, eos, "flight")
    body = rn.make_body(S, a.form, arrays, eos, strength, dx=dx)
    hydro, integrator, controller = setup(S, body, dx)
    K = float(table.bulk_modulus(a.T))
    c = np.sqrt(K/(a.scale*float(table.rho_free(a.T))))
    integrator.lastDt = 0.1*CFL*dx/c
    for _ in range(a.steps):
        integrator.step(1e40)
        controller.totalSteps += 1
        controller.doPeriodicWork()
    g = rn.gather_by_id(body.nodes, dict(pid=body.pid, P=hydro.pressure[0], S=body.nodes.deviatoricStress(),
                                         v=body.nodes.velocity(), rho=body.nodes.massDensity(),
                                         u=body.nodes.specificThermalEnergy(), x=body.nodes.positions(),
                                         cs=hydro.soundSpeed[0]), body.dim)
    if comm.rank == 0:
        np.savez(a.out, time=controller.time(), dt=integrator.lastDt, steps=controller.totalSteps, K=K,
                 rho0=rho[0], u0=arrays["u"][0], processes=comm.size, **g)
    return 0


def bench(S, a, table):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    load_avg = os.getloadavg() if comm.rank == 0 else None
    load_avg = comm.bcast(load_avg, root=0)
    if load_avg[0] > a.max_load:
        if comm.rank == 0:
            print(f"bench: load average {load_avg[0]:.2f} > {a.max_load}; not starting")
        return 2
    b = load(os.path.join(a.body, "body.npz"))
    dx = float(b["dx_m"])
    x = b["x_ref"]
    u = table.enthalpy(b["T0"])
    arrays = dict(x_ref=x, pid=b["pid"], m=b["m"], rho=b["rho"], u=u, v=EDOT*x)
    T300 = material.T_ZERO
    K, G = float(table.bulk_modulus(T300)), float(table.shear_modulus(T300))
    rho_ref = float(table.rho_free(T300))
    if a.material == "thermal":
        eos = make_eos(S, table, 1.0)
        strength = make_strength(S, eos, "flight")
    else:
        eos = S.MurnaghanEquationOfState(rho_ref, 0.2, 4.0, 1.0, K, 26.98, S.MKS())
        strength = S.ConstantStrength(G, 2.0e7, eos)
    body = rn.make_body(S, a.form, arrays, eos, strength, dx=dx)
    hydro, integrator, controller = setup(S, body, dx)
    integrator.lastDt = 0.1*CFL*dx/np.sqrt((K + 4.0/3.0*G)/rho_ref)
    for _ in range(a.warmup):
        integrator.step(1e40)
        controller.totalSteps += 1
        controller.doPeriodicWork()
    walls, dts = [], []
    for _ in range(a.steps):
        comm.Barrier()
        t0 = time.perf_counter()
        integrator.step(1e40)
        controller.totalSteps += 1
        controller.doPeriodicWork()
        comm.Barrier()
        walls.append(time.perf_counter() - t0)
        dts.append(integrator.lastDt)
    n = comm.allreduce(body.nodes.numInternalNodes)
    if comm.rank == 0:
        np.savez(a.out, walls=np.array(walls), dts=np.array(dts), n=n, processes=comm.size,
                 load=np.array(load_avg), material=a.material)
        w = np.array(walls)
        print(f"bench {a.material}: n={n} procs={comm.size} median {np.median(w)*1e3:.1f} ms mean "
              f"{w.mean()*1e3:.1f} ms min {w.min()*1e3:.1f} ms dt {dts[-1]:.3e} load {load_avg[0]:.2f}")
    return 0


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("methods", "rest", "bench"))
    p.add_argument("--form", choices=rn.FORMS, required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--input")
    p.add_argument("--body")
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--T", type=float, default=300.0)
    p.add_argument("--steps", type=int, default=100)
    p.add_argument("--warmup", type=int, default=10)
    p.add_argument("--material", choices=("thermal", "builtin"), default="thermal")
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    S, _ = rn.spheral_module(a.form)
    table = material.MaterialTable.load(a.table)
    return {"methods": methods, "rest": rest, "bench": bench}[a.mode](S, a, table)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
