"""Extract a few frames' prepared load tables from the MVP flight's prepared input (Spheral M2 plan, Task 4).

    "$PY" tests/fixtures/spheral_frag/make_loads_fixture.py [--prepared spheral_output/prepare/<prepared>] \
        [--frames 0,1,182,183] [--out tests/fixtures/spheral_frag/loads/loads_mvp.npz]

writes the rows of `<prepared>/loads.npz` for the frames (spheral_frag.loads.stack's arrays: k, time_s, theta_deg,
edges_deg, p, tau, area, lee, theta_last_deg, p_stag; windward tables, no lee extension) and each frame's surface
area (`surface_area_m2`, the sum of its patch areas), so that the runner's package tests run without the git-ignored
prepared input. Deterministic: the arrays regenerate bitwise."""
import argparse
import os
import sys

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
PREPARED = os.path.join("spheral_output", "prepare",
                        "prep_model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_"
                        "material-AA7075_scheil_deeprunoff-off_k00000-01254")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--prepared", default=PREPARED)
    ap.add_argument("--frames", default="0,1,182,183")
    ap.add_argument("--out", default=os.path.join("tests", "fixtures", "spheral_frag", "loads", "loads_mvp.npz"))
    a = ap.parse_args(argv)
    sys.path.insert(0, REPO_ROOT)
    from spheral_frag import frames
    ks = [int(k) for k in a.frames.split(",")]
    prepared = os.path.join(REPO_ROOT, a.prepared)
    with np.load(os.path.join(prepared, "loads.npz"), allow_pickle=False) as z:
        full = {k: z[k] for k in z.files}
    rows = [int(np.flatnonzero(full["k"] == k)[0]) for k in ks]
    out = {k: (v if k in ("theta_deg", "edges_deg") else v[rows]) for k, v in full.items()}
    out["surface_area_m2"] = np.array([float(np.sum(frames.load_prepared_frame(
        os.path.join(prepared, "frames", frames.prepared_frame_name(k))).patch["area"])) for k in ks])
    path = os.path.join(REPO_ROOT, a.out)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **out)
    print({k: np.shape(v) for k, v in out.items()}, out["surface_area_m2"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
