"""Fitted-sphere windward loads (2026-11-13), Task 5: the upstream-ray shelter test and the sphere fit on saved skin states.

For each saved state (default: the 100 mm and 50 mm `--hand-down off` states of the task5_compare flights in the subplan20
worktree) this loads the body, sets `windward_loads = "sphere"` and prints the corners ahead of the mass-centroid plane (the
candidates), how many of them the upstream ray finds sheltered and exposed, the sphere fitted to the exposed corners, the flat-face
bound, the radius the loads use, the fitted centre's distance along the flight direction from the mesh's coordinate origin (the
sphere's original centre) [mm], the rms of the exposed corners about the sphere [mm], and the best-of-3 wall time [s] of
`shelter.upstream_hits` (the grid version) and `shelter.upstream_hits_plain` (the reference). Also checks that the two agree.

Usage, from `Reentry Model/`:  python analysis/fitted_sphere_check.py [STATE.pkl ...]
"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reentry_model.skyscrapers import restart, shelter  # noqa: E402

HERE = Path(__file__).resolve()
# the repository root is the folder above `.claude/worktrees/` when this runs in a worktree, else two levels up
REPO = next((p.parent.parent.parent for p in HERE.parents if p.parent.name == "worktrees" and p.parent.parent.name == ".claude"),
            HERE.parents[2])
FOLDER = (REPO / ".claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/"
          "reentry_model_output/task5_compare")
DEFAULTS = [FOLDER / "d100__handdown-off.state.pkl", FOLDER / "d050__handdown-off.state.pkl"]


def best_of(function, n=3):
    times, result = [], None
    for _ in range(n):
        start = time.perf_counter()
        result = function()
        times.append(time.perf_counter() - start)
    return min(times), result


def check(path):
    print(f"\n== {Path(path).name}", flush=True)
    body = restart.load(str(path))["run"].body
    body.windward_loads = "sphere"
    body._sphere_raw = None
    s, v = body.surface, body.v_hat
    ahead = body._leeward()[1]
    t_grid, hit = best_of(lambda: shelter.upstream_hits(s.points, s.faces, s.normals, s.areas, v, candidates=ahead))
    t_plain, hit_plain = best_of(lambda: shelter.upstream_hits_plain(s.points, s.faces, s.normals, s.areas, v, candidates=ahead))
    f = body.sphere_fit
    print(f"corners {len(s.points)}, faces {len(s.faces)}; candidates (ahead of the mass-centroid plane) {int(ahead.sum())}")
    print(f"sheltered {int(f.sheltered.sum())}, exposed {int(f.exposed.sum())} (fitted to {f.points} corners), concave {f.concave}")
    print(f"fitted radius {f.fitted * 1e3:.2f} mm, flat-face bound {f.bound * 1e3:.2f} mm, radius used {f.used * 1e3:.2f} mm")
    print(f"fitted centre along the flight direction from the original centre {float(f.centre @ v) * 1e3:.2f} mm, "
          f"rms {f.rms * 1e3:.3f} mm")
    print(f"upstream_hits {t_grid:.3f} s, upstream_hits_plain {t_plain:.3f} s (best of 3); "
          f"agree {bool(np.array_equal(hit, hit_plain))}")


if __name__ == "__main__":
    for p in (sys.argv[1:] or DEFAULTS):
        check(p)
