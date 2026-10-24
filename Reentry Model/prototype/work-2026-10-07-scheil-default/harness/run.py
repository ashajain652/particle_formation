"""One plain CLI run of a prototype copy, kept at full precision.

    run.py <copy dir> <out dir> -- <cli run arguments>      (--outdir <out dir> is appended)

No seeding here and no instrumentation: the copy's own cli.main runs exactly as a user would run it, so a copy with the
2026-10-05 amendment seeds itself and an older copy does not. The only addition is a wrapper around
coupled.CoupledRun.run that keeps the returned history and the body, so that after the run every history column and the
final state (nodal temperatures, element fractions, film and deep accounts, source-row count) are written to
<out dir>/full.npz at full precision -- the CSV the CLI writes keeps 9 significant digits, which cannot show a last-bit
difference. Refuses an existing output directory (fact 35: a stale output must never pass for a fresh one).
"""
import json
import os
import sys
import time

copy_dir, out_dir = sys.argv[1], sys.argv[2]
argv = sys.argv[sys.argv.index("--") + 1:] + ["--outdir", out_dir]
if os.path.exists(out_dir):
    sys.exit("refusing: {} exists".format(out_dir))
os.makedirs(out_dir)
sys.path.insert(0, copy_dir)
import numpy as np  # noqa: E402
from reentry_model import cli, coupled  # noqa: E402

KEPT = {}
original_run = coupled.CoupledRun.run


def run_and_keep(self):
    history = original_run(self)
    KEPT["history"], KEPT["body"] = history, self.body
    return history


coupled.CoupledRun.run = run_and_keep
started = time.perf_counter()
code = cli.main(argv)
wall = time.perf_counter() - started
arrays = {}
if "history" in KEPT:
    h, b = KEPT["history"], KEPT["body"]
    arrays.update({"col__" + k: np.asarray(v) for k, v in h.columns.items()})
    arrays["final__T"] = np.asarray(b.solver.temperature())
    for k in ("phi", "m_f", "m_d"):
        if hasattr(b, k):
            arrays["final__" + k] = np.asarray(getattr(b, k))
    if hasattr(b, "mesh") and hasattr(b.mesh, "active"):
        arrays["final__active"] = np.asarray(b.mesh.active)
    np.savez(os.path.join(out_dir, "full.npz"), **arrays)
record = {"argv": argv, "copy": copy_dir, "exit_code": code, "wall_s": wall, "pid": os.getpid(),
          "seed_attr": getattr(cli, "DEFAULT_SEED", None), "started": time.strftime("%Y-%m-%d %H:%M:%S")}
with open(os.path.join(out_dir, "record.json"), "w") as fh:
    json.dump(record, fh, indent=1)
print("run.py: exit {} in {:.0f} s -> {}".format(code, wall, out_dir))
sys.exit(code)
