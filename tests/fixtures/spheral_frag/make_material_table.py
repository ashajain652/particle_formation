"""Write a material table and the finite-element material's own values at sample temperatures (Spheral M1, Task 3).

    "$PY" tests/fixtures/spheral_frag/make_material_table.py --fe-package prototype/work-2026-10-08-spheral-mvp/code --material AA7075_scheil \
        --out tests/fixtures/spheral_frag/material --samples 5000

writes `<out>/material_<name>.npz` + `.json` (exactly what `prepare` writes as `material_table.npz`) and
`<out>/material_<name>_fe_samples.npz` (spheral_frag.fe.fe_samples: T, the FE's h, f_l, c_p, inverse and fully
liquid enthalpy). The committed fixture `tests/fixtures/spheral_frag/material/` is the command above; it lets the
material tests run on machines without the git-ignored prototype. The tests also run this script in a subprocess
(one finite-element package per process) with 100,000 samples, against `SPHERAL_FRAG_FE_PACKAGE` or the repository's
own package. Prints the measured agreement as JSON."""
import argparse
import json
import os
import shlex
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fe-package", default=None, help="directory holding reentry_model (default: the repository)")
    ap.add_argument("--material", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samples", type=int, default=5000, help="uniform random temperatures in 250-1,500 K")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    sys.path.insert(0, REPO_ROOT)
    import numpy as np
    from spheral_frag import fe

    if args.fe_package and not os.path.isdir(args.fe_package):
        print("no finite-element package directory {}".format(args.fe_package), file=sys.stderr)
        return 2
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    command = "tests/fixtures/spheral_frag/make_material_table.py " + " ".join(shlex.quote(a) for a in (argv or sys.argv[1:]))
    table = fe.build_material_table(args.material, os.path.join(out, "material_{}.npz".format(args.material)),
                                    package=args.fe_package, command=command)
    mat = fe.fe_material(args.material, args.fe_package)
    s = fe.fe_samples(mat, fe.sample_temperatures(table, args.samples, args.seed))
    np.savez(os.path.join(out, "material_{}_fe_samples.npz".format(args.material)), **s)
    print(json.dumps(fe.compare_samples(table, s)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
