"""Spheral M2 check 13's crossings probe (Task 10; run inside the container by analysis/spheral_m2_check13.py
--only probe):

    spheral_frag/container/spheral [-n N] analysis/spheral_m2_check13_probe.py RUN_DIR K OUT N1,N2,... [EPS_TENSILE]

reruns the replay interval K -> K+1 of the run in RUN_DIR from its checkpoint K at each crossings number N_i (every
other setting the run's defaults, as the run itself), and writes each end state to OUT/k<K>_N<N>.npz (the gathered
state of `Replay.attempt`, pid order, with the interval's scale, steps, wall time and internal KE ratio). Nothing of
the run is touched: the probe's Replay writes into OUT only, and `attempt` writes no checkpoint. The states at N and
at the largest N show how far the accepted state at N is from the converged one, which the monitor estimates from N
and 2 N alone (run.py, "The crossings-error monitor").

EPS_TENSILE (optional, default 0: the runner's SolidSPH) builds the hydro with Spheral's Monaghan artificial stress
`epsTensile` (SPH's factory argument; nTensile its default 4) instead, a diagnostic of the axial pairing that RZ
1.1 mm develops under the core's axial tension (module docstring of analysis/spheral_m2_check13.py); the states go to
OUT/k<K>_N<N>_epst<EPS>.npz. The runner itself has no such option."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from spheral_frag.runner import checkpoint, run  # noqa: E402


def main(argv):
    src_run, k, out = argv[0], int(argv[1]), argv[2]
    Ns = [int(x) for x in argv[3].split(",")]
    eps = float(argv[4]) if len(argv) > 4 else 0.0
    suffix = "_epst{:g}".format(eps) if eps > 0.0 else ""
    if eps > 0.0:
        def make_hydro(S, scheme, body):
            assert scheme == "sph"
            return S.SPH(dataBase=body.db, W=body.kernel, cfl=run.CFL, compatibleEnergyEvolution=True,
                         densityUpdate=S.IntegrateDensity, XSPH=False, ASPH=False, epsTensile=eps)
        run.make_hydro = make_hydro
    with open(os.path.join(src_run, "meta.json")) as fh:
        meta = json.load(fh)
    a = run.build_parser().parse_args(["replay", "--body", meta["body_dir"], "--prepared", meta["prepared_dir"],
                                       "--frames", "{}-{}".format(*meta["frames"]), "--out", out,
                                       "--max-crossings", str(max(Ns)),
                                       "--hourglass", repr(float(meta["brackets"].get("hourglass", 0.0)))])
    inp = run.Inputs(a)
    if inp.name != meta["run_name"]:
        raise SystemExit("the probe's settings give {}, not the run's {}".format(inp.name, meta["run_name"]))
    os.makedirs(out, exist_ok=True)
    rp = run.Replay(a, inp, out)
    ck = run.read_ck(os.path.join(src_run, "checkpoints", checkpoint.filename(k)))
    k_next = inp.ks[inp.ks.index(k) + 1]
    for N in Ns:
        path = os.path.join(out, "k{}_N{}{}.npz".format(k, N, suffix))
        if os.path.isfile(path):
            continue
        st, info = rp.attempt(k, k_next, ck, N, "p{}".format(N))
        if run.comm().rank == 0:
            np.savez(path + ".tmp.npz", scale=info["scale"], steps=info["steps"], wall=info["wall_s"],
                     ke_internal_ratio=info["ke_internal_ratio"], k0=k, k1=k_next, crossings=N, eps_tensile=eps,
                     **{key: v for key, v in st.items()})
            os.replace(path + ".tmp.npz", path)
            print("k {} -> {}  N {}  steps {}  wall {:.1f} s  internal KE/E {:.3e}".format(
                k, k_next, N, info["steps"], info["wall_s"], info["ke_internal_ratio"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
