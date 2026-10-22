"""Spheral M2 Task 10b's remedy probe (run inside the container by analysis/spheral_m2_task10b.py):

    spheral_frag/container/spheral [-n N] analysis/spheral_m2_task10b_probe.py interval CONFIG SRC_RUN K OUT N1,N2,...
    spheral_frag/container/spheral [-n N] analysis/spheral_m2_task10b_probe.py timing CONFIG SRC_RUN K OUT

`interval` reruns the replay interval K -> K+1 of the run in SRC_RUN from its checkpoint K at each crossings number
N_i with the configuration CONFIG applied, and writes each end state to OUT/k<K>_N<N>.npz (the gathered state of
`Replay.attempt`, pid order, with the interval's scale, steps, wall time, mean step and internal KE ratio), as
analysis/spheral_m2_check13_probe.py does for the runner's defaults. `timing` runs one attempt at 1 crossing over
1 % of the interval (the same steps per crossing, about 330 at RZ 1.1 mm) and writes OUT/timing.json: Spheral's
mean time step and cost per step under CONFIG (the epsTensile diagnosis).

CONFIG is JSON; every key optional, the runner's defaults otherwise:
  kernel    "C2" | "C4" (default) | "C6" | "B" (Spheral's WendlandC2/C4/C6 or cubic B-spline TableKernel)
  nperh     smoothing lengths per spacing (default 2.01): the checkpoint's H is scaled by 2.01 / nperh and the node
            list built with it (IdealH then evolves it); the loads' surface detection keeps its h = 2.01 dx
  eps, ntensile   SPH's Monaghan artificial stress epsTensile and nTensile (default 0 and 4)
  xsph      SPH's XSPH position update (default false)
  Q         "default" (the factory's LimitedMonaghanGingold) | "mg" (unlimited MonaghanGingold, Cl = Cq = 1)
  zeta      the damping ratio (run.py --damping)
  hg, beta  the hourglass stiffness kappa (run.py --hourglass) and its viscous ratio (packages/hourglass.py BETA)

Nothing of the source run is touched: Replay's run directory is OUT, and `attempt` writes no checkpoint."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from spheral_frag.runner import checkpoint, run, nodes as rn  # noqa: E402

KERNELS = {"C2": "WendlandC2Kernel", "C4": "WendlandC4Kernel", "C6": "WendlandC6Kernel", "B": "BSplineKernel"}


def configure(cfg, rp):
    """Apply CONFIG to the runner's module functions and the Replay `rp` (module docstring)."""
    S = rp.S
    rp.kernel = S.TableKernel(getattr(S, KERNELS[cfg.get("kernel", "C4")])(), 1000)
    nperh = float(cfg.get("nperh", rn.NPERH))
    orig_make_body = rn.make_body

    def make_body(S_, form, arrays, eos, strength, nperh_=rn.NPERH, **kw):
        arrays = dict(arrays)
        f = rn.NPERH/nperh
        if f != 1.0:
            arrays["H"] = np.asarray(arrays["H"])*f
            if "H_zz" in arrays:
                arrays["H_zz"] = np.asarray(arrays["H_zz"])*f
        body = orig_make_body(S_, form, arrays, eos, strength, nperh, **kw)
        body.nperh = rn.NPERH                   # the loads' surface detection keeps h = 2.01 dx
        return body
    rn.make_body = make_body

    def make_hydro(S_, scheme, body):
        kw = {}
        if cfg.get("Q", "default") == "mg":
            kw["Q"] = S_.MonaghanGingoldViscosity(Clinear=1.0, Cquadratic=1.0, kernel=body.kernel)
        return S_.SPH(dataBase=body.db, W=body.kernel, cfl=run.CFL, compatibleEnergyEvolution=True,
                      densityUpdate=S_.IntegrateDensity, XSPH=bool(cfg.get("xsph", False)), ASPH=False,
                      epsTensile=float(cfg.get("eps", 0.0)), nTensile=float(cfg.get("ntensile", 4.0)), **kw)
    run.make_hydro = make_hydro
    if "beta" in cfg:
        from spheral_frag.runner.packages import hourglass as phg
        orig = phg.make_hourglass
        phg.make_hourglass = lambda S_, body, kappa, table, beta=None: orig(S_, body, kappa, table,
                                                                           float(cfg["beta"]))


def main(argv):
    mode, cfg = argv[0], json.loads(argv[1])
    src, k, out = argv[2], int(argv[3]), argv[4]
    with open(os.path.join(src, "meta.json")) as fh:
        meta = json.load(fh)
    Ns = [int(x) for x in argv[5].split(",")] if mode == "interval" else [10]
    extra = ["--max-crossings", str(max(Ns))]
    if "zeta" in cfg:
        extra += ["--damping", repr(float(cfg["zeta"]))]
    extra += ["--hourglass", repr(float(cfg.get("hg", 0.0)))]       # the runner's default is 0.01 since Task 10b
    a = run.build_parser().parse_args(["replay", "--body", meta["body_dir"], "--prepared", meta["prepared_dir"],
                                       "--frames", "{}-{}".format(*meta["frames"]), "--out", out] + extra)
    inp = run.Inputs(a)
    os.makedirs(out, exist_ok=True)
    rp = run.Replay(a, inp, out)
    configure(cfg, rp)
    ck = run.read_ck(os.path.join(src, "checkpoints", checkpoint.filename(k)))
    k_next = inp.ks[inp.ks.index(k) + 1]
    if mode == "timing":
        t0 = inp.times[inp.column(k)]
        inp.times = list(inp.times)
        inp.times[inp.column(k_next)] = t0 + 0.01*(inp.times[inp.column(k_next)] - t0)
        st, info = rp.attempt(k, k_next, ck, 1, "t")
        if run.comm().rank == 0:
            doc = dict(config=cfg, steps=info["steps"], wall_s=info["wall_s"], dt_mean_s=info["dt_mean_s"],
                       ms_per_step=1e3*info["wall_per_step_s"], ke_internal_ratio=info["ke_internal_ratio"],
                       processes=run.comm().size)
            run.write_json(os.path.join(out, "timing.json"), doc)
            print(json.dumps(doc), flush=True)
        return 0
    for N in Ns:
        path = os.path.join(out, "k{}_N{}.npz".format(k, N))
        if os.path.isfile(path):
            continue
        st, info = rp.attempt(k, k_next, ck, N, "p{}".format(N))
        if run.comm().rank == 0:
            np.savez(path + ".tmp.npz", scale=info["scale"], steps=info["steps"], wall=info["wall_s"],
                     dt_mean=info["dt_mean_s"], ke_internal_ratio=info["ke_internal_ratio"],
                     hourglass_energy_ratio=info.get("hourglass_energy_ratio", 0.0), k0=k, k1=k_next, crossings=N,
                     config=json.dumps(cfg, sort_keys=True), processes=run.comm().size,
                     **{key: v for key, v in st.items()})
            os.replace(path + ".tmp.npz", path)
            print("k {} -> {}  N {}  steps {}  wall {:.1f} s  internal KE/E {:.3e}".format(
                k, k_next, N, info["steps"], info["wall_s"], info["ke_internal_ratio"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
