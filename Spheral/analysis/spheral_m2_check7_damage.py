"""Spheral M2 check 7, Step 3: does Spheral's damage run with each hydro scheme (spec §14.1 condition 2)?

    spheral_frag/container/spheral analysis/spheral_m2_check7_damage.py --scheme sph|crksph --form 3d|rz
        [--out spheral_output/m2/check7/damage] [--steps 400] [--edot 1e4] [--no-damage]

Runs under Spheral's Python (the container), not drama_env. The body and material are M0's benchmark's
(spheral_frag/m0/bench.py: a lattice sphere, Murnaghan with exponent 1, Steinberg-Guinan 6061-T6, the
`ProbabilisticDamageModel` with the library's aluminium Weibull constants and the documented seed), on a 10 mm
radius sphere at 1 mm spacing, one process; the 3D masses are set to rho0 dx^3 (the 3D lattice generator's are wrong,
M2 plan fact 3; RZ's are the rings' already). The hydro is built as the replay builds it (run.make_hydro: SolidSPH
or SolidCRKSPH, integrated density, no XSPH; CRKSPH with mass-over-density volumes, and in RZ without compatible
energy), with CheapSynchronousRK2 and the controller without the initial H iteration (as the replay; under
SolidCRKSPHRZ SpheralController's iterateIdealH segfaults in RKUtilities::evaluateKernelAndGradient before the
corrections exist, measured 2026-10-11). A uniform expansion v = edot x with edot = 1e4 1/s pulls the sphere into
tension well past the Weibull activation strains (400 steps of about 3-5e-8 s: a volumetric strain of 0.4-0.7), so
that a damage model that runs must show damage growing. Recorded per run (`<out>/<scheme>_<form>.json`): whether it
constructed and stepped, the steps done, the error if any, the largest damage eigenvalue and the share of particles
with damage above 0.01 and 0.99 at the end, the strain reached, the cost per step. `--no-damage` runs the same without
the damage model (diagnosis). Not a validation of the damage model (M3's): only whether it runs with the scheme at all."""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

R_M = 0.010
DX_M = 0.001
EDOT = 1.0e4


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--scheme", choices=("sph", "crksph"), required=True)
    p.add_argument("--form", choices=("3d", "rz"), required=True)
    p.add_argument("--steps", type=int, default=400)
    p.add_argument("--edot", type=float, default=EDOT)
    p.add_argument("--no-damage", action="store_true", help="the same run without the damage model (diagnosis)")
    p.add_argument("--out", default=os.path.join("spheral_output", "m2", "check7", "damage"))
    a = p.parse_args(argv)
    out = dict(scheme=a.scheme, form=a.form, steps_requested=a.steps, radius_m=R_M, dx_m=DX_M, edot=a.edot,
               constructed=False, stepped=False, steps_done=0, error=None)
    try:
        run(a, out)
    except Exception as e:                       # the point of the probe: record whatever Spheral says
        out["error"] = f"{type(e).__name__}: {e}"
        out["traceback"] = traceback.format_exc()
    os.makedirs(a.out, exist_ok=True)
    tag = "_nodamage" if a.no_damage else ""
    with open(os.path.join(a.out, f"{a.scheme}_{a.form}{tag}.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "traceback"}, indent=2))
    return 0


def run(a, out):
    from spheral_frag.m0 import bench as b
    from spheral_frag.runner import run as rr
    if a.form == "3d":
        import Spheral3d as S
        from PeanoHilbertDistributeNodes import distributeNodes3d as distribute
    else:
        import SpheralRZ as S
        from PeanoHilbertDistributeNodes import distributeNodes2d as distribute
    from SpheralController import SpheralController
    units = S.MKS()
    eos = S.MurnaghanEquationOfState(b.RHO0, b.ETAMIN, b.ETAMAX, b.MURNAGHAN["n"], b.MURNAGHAN["K"],
                                     b.MURNAGHAN["atomicWeight"], units)
    zero = S.NinthOrderPolynomialFit(*([0.0]*10))
    melt = S.NinthOrderPolynomialFit(*([b.EMELT_PER_VOLUME] + [0.0]*9))
    strength = S.SteinbergGuinanStrength(eos, b.SG["G0"], b.SG["A"], b.SG["B"], b.SG["Y0"], b.SG["Ymax"],
                                         b.SG["Yp"], b.SG["beta"], b.SG["gamma0"], b.SG["nhard"], zero, melt)
    W = S.TableKernel(S.WendlandC4Kernel(), 1000)
    big = 10.0*R_M
    dim = 3 if a.form == "3d" else 2
    box = (S.Vector(*([-big]*dim)), S.Vector(*([big]*dim)))
    nodes = S.makeSolidNodeList("probe", eos, strength, nPerh=b.NPERH, kernelExtent=W.kernelExtent,
                                hmin=1e-3*DX_M, hmax=100.0*DX_M, rhoMin=b.ETAMIN*b.RHO0, rhoMax=b.ETAMAX*b.RHO0,
                                xmin=box[0], xmax=box[1])
    gen, _ = b.lattice_generator(a.form, R_M, DX_M)
    distribute((nodes, gen))
    n = nodes.numInternalNodes
    out["n_particles"] = n
    if a.form == "3d":
        mass = nodes.mass()
        for i in range(n):
            mass[i] = b.RHO0*DX_M**3
    nodes.specificThermalEnergy(S.ScalarField("eps0", nodes, 0.0))
    pos, vel = nodes.positions(), nodes.velocity()
    for i in range(n):
        vel[i] = pos[i]*a.edot
    db = S.DataBase()
    db.appendNodeList(nodes)

    class _Body:                                  # what run.make_hydro reads
        pass
    body = _Body()
    body.db, body.kernel, body.form = db, W, a.form
    hydro = rr.make_hydro(S, a.scheme, body)
    integ = S.CheapSynchronousRK2Integrator(db)
    integ.appendPhysicsPackage(hydro)
    damage = None
    if not a.no_damage:
        damage = S.ProbabilisticDamageModel(materialName="aluminum", units=units, nodeList=nodes, kernel=W,
                                            seed=b.SEED)
        integ.appendPhysicsPackage(damage)
    out["damage_model"] = not a.no_damage
    cs0 = math.sqrt((b.MURNAGHAN["K"] + (4.0/3.0)*b.SG["G0"])/b.RHO0)
    integ.lastDt = 0.1*b.CFL*DX_M/cs0
    integ.dtMin, integ.dtMax, integ.verbose = 1e-15, 1.0, False
    ctl = SpheralController(integ, W, statsStep=10**9, printStep=10**9, redistributeStep=None, restartStep=None,
                            vizStep=None, vizTime=None, SPH=True, iterateInitialH=False,
                            **rr.controller_options(S, a.scheme))
    out["constructed"] = True
    out["packages"] = [type(pk).__name__ for pk in integ.physicsPackages()]
    for attr in ("minFlawStrain", "maxFlawStrain", "kWeibull", "mWeibull", "volumeMultiplier"):
        try:
            if damage is None:
                break
            v = getattr(damage, attr)
            out[attr] = v() if callable(v) else float(v)
        except Exception:                       # noqa: BLE001 -- informational only
            pass
    t0 = time.perf_counter()
    for _ in range(a.steps):
        integ.step(1e40)
        ctl.totalSteps += 1
        out["steps_done"] += 1
        out["stepped"] = True
    wall = time.perf_counter() - t0
    D = nodes.damage()
    dmax = [D[i].eigenValues().maxElement() for i in range(n)]
    t = float(integ.currentTime)
    out.update(time_s=t, volumetric_strain=(1.0 + a.edot*t)**3 - 1.0, damage_max=max(dmax),
               share_damage_gt_0p01=sum(d > 0.01 for d in dmax)/n, share_damage_gt_0p99=sum(d > 0.99 for d in dmax)/n,
               wall_per_step_s=wall/a.steps, finite=all(math.isfinite(d) for d in dmax))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
