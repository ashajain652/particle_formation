"""Write the committed synthetic finite-element runs of tests/fixtures/spheral_frag/ (Spheral M1 plan, Task 2).

    "$PY" tests/fixtures/spheral_frag/make_fixtures.py

Run once; rerun to refresh after a change to tests/spheral_frag_synthetic.py or to spheral_frag/contract.py (a renamed
field renames it in the files). Each fixture directory is removed and rewritten. The builders are deterministic
(gmsh's HXT sphere and single-threaded Delaunay dumbbell), and tests/test_spheral_frag_fixtures.py checks that the
regenerated sphere run equals the committed one bitwise. Needs drama_env (gmsh, pyvista)."""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
for path in (REPO_ROOT, os.path.join(REPO_ROOT, "tests")):
    if path not in sys.path:
        sys.path.insert(0, path)

import spheral_frag_synthetic as syn  # noqa: E402

FIXTURES = {
    "sphere_run": lambda: syn.sphere_run_parts()[:3],
    "dumbbell_frame": syn.dumbbell_parts,
    "slab_frame": syn.slab_parts,
}
MAX_TOTAL_BYTES = 3 * 1024 * 1024        # M1 plan, Global constraints: synthetic fixtures below 3 MB in total


def directory_size(path):
    return sum(os.path.getsize(os.path.join(d, f)) for d, _, files in os.walk(path) for f in files)


def main():
    total = 0
    for name, parts in FIXTURES.items():
        out = os.path.join(HERE, name)
        shutil.rmtree(out, ignore_errors=True)
        frames, history, doc = parts()
        syn.write_fe_run(out, frames, history, doc)
        size = directory_size(out)
        total += size
        print("{:15s} {:2d} frame(s), {:6d} tets, {:5d} patches, {:8.1f} kB -> {}".format(
            name, len(frames), len(frames[-1].tets), len(frames[-1].faces), size / 1024.0, out))
    print("total {:.1f} kB".format(total / 1024.0))
    if total >= MAX_TOTAL_BYTES:
        sys.exit("synthetic fixtures are {} bytes, above the 3 MB budget".format(total))


if __name__ == "__main__":
    main()
