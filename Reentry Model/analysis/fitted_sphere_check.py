"""Fitted-sphere windward loads (2026-11-13), Task 5: the upstream-ray shelter test and the sphere fit on saved skin states.

For each saved state (default: the 100 mm and 50 mm `--hand-down off` states of the task5_compare flights in the subplan20
worktree) this loads the body, sets `windward_loads = "sphere"` and prints the corners ahead of the mass-centroid plane (the
candidates), how many of them the upstream ray finds sheltered and exposed, the sphere fitted to the exposed corners, the flat-face
bound, the radius the loads use, the fitted centre's distance along the flight direction from the mesh's coordinate origin (the
sphere's original centre) [mm], the rms of the exposed corners about the sphere [mm], and the best-of-3 wall time [s] of
`shelter.upstream_hits` (the grid version) and `shelter.upstream_hits_plain` (the reference). Also checks that the two agree.

Added in the final review's fix wave (2026-11-13), so that the documents' comparison with the earlier rule has a committed source:
- the push's direction: the fitted sphere's unit meridian at each exposed corner (the direction away from the fitted
  stagnation point along the sphere, as `SkinBody.load_meridian` builds it) projected onto the corner's base tangent plane;
  its length there before `load_meridian` makes it a unit vector again (1 where the base is tangent to the fitted sphere,
  0 where the meridian is normal to the base), the shortest lengths and how many fall below 0.5, over the exposed corners
  whose meridian is defined, and over those whose ray angle is below 90 degrees (the only ones the push acts on under the
  classic loads, on which the saved states flew: their surface flow gives no pressure gradient or shear beyond 90 degrees;
  on the interim gas loads' provider, the default since the merge of 2026-11-15, the lee loads act beyond 90 degrees too,
  model_assumptions section 15.5);
- the earlier rule on the same body, `windward_loads = "normals"` with the `--size-feedback` and `--leeward-test` that a
  `run --windward-loads normals` resolves by default (through `cli`'s own parser and helpers): its 30-degree cap fit, bound and
  radius used, and the ratio of the two rules' radii, which sets the stagnation heating (one over the square root of the
  radius) and the scale of the pressure gradient (one over the radius);
- the corners that change sides: windward under `normals` (its `windward_corners()`) but sheltered under `sphere`, and
  leeward under `normals` (not its `windward_corners()`) but exposed under `sphere`, with how many of the latter are windward
  under `sphere` (their ray angle at most 90 degrees), so heated.

Usage, from `Reentry Model/`:  python analysis/fitted_sphere_check.py [STATE.pkl ...]
"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reentry_model import cli  # noqa: E402
from reentry_model.skyscrapers import restart, shelter  # noqa: E402

HERE = Path(__file__).resolve()
# the repository root is the folder above `.claude/worktrees/` when this runs in a worktree, else two levels up
REPO = next((p.parent.parent.parent for p in HERE.parents if p.parent.name == "worktrees" and p.parent.parent.name == ".claude"),
            HERE.parents[2])
FOLDER = (REPO / ".claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/"
          "reentry_model_output/task5_compare")
DEFAULTS = [FOLDER / "d100__handdown-off.state.pkl", FOLDER / "d050__handdown-off.state.pkl"]


def normals_defaults():
    """(--windward-loads, --size-feedback, --leeward-test) as a `run --windward-loads normals` resolves them by default, through
    the command line's own parser and helpers; the three required flight arguments are placeholders these helpers do not read."""
    args = cli.build_parser().parse_args(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "65",
                                          "--windward-loads", "normals"])
    return cli.skin_windward_loads_of(args), cli.skin_size_feedback_of(args), cli.skin_leeward_test_of(args)


def shortest(x, k=3):
    return ", ".join(f"{y:.3f}" for y in np.sort(x)[:k])


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

    # the push's direction: the fitted unit meridian projected onto each exposed corner's base tangent plane (load_meridian)
    if not f.concave:
        n, k = body.columns.normal, np.flatnonzero(f.exposed)
        u = s.points[k] - f.centre
        u /= np.linalg.norm(u, axis=1)[:, None]
        away = -(v[None, :] - (u @ v)[:, None] * u)          # along the fitted sphere's meridian, away from its stagnation point
        size = np.linalg.norm(away, axis=1)
        ok = size > 1e-12                                     # no meridian at the fitted stagnation point itself
        unit, nk = away[ok] / size[ok, None], n[k][ok]
        proj = np.linalg.norm(unit - np.einsum("ij,ij->i", unit, nk)[:, None] * nk, axis=1)
        front = np.degrees(np.arccos(np.clip(u[ok] @ v, -1.0, 1.0))) < 90.0
        print(f"fitted unit meridian projected onto the base tangent plane, its length: over the {int(ok.sum())} exposed corners "
              f"with a meridian, shortest {shortest(proj)}, {int((proj < 0.5).sum())} below 0.5; over the {int(front.sum())} of "
              f"them with a ray angle below 90 deg, shortest {shortest(proj[front])}, {int((proj[front] < 0.5).sum())} below 0.5")
    else:
        print("no valid sphere: the push keeps each base normal's own meridian (along_meridian)")

    # the earlier rule on the same body, at the command line's defaults for it
    windward_sphere = body.windward_corners()
    saved = {k: vars(body).get(k, "(class default)") for k in ("size_feedback", "leeward_test")}
    rule, size_feedback, leeward_test = normals_defaults()
    body.windward_loads, body.size_feedback, body.leeward_test = rule, size_feedback, leeward_test
    body._nose_raw = None                                     # the cap fit made again under these settings
    g = body.nose_fit
    windward_normals = body.windward_corners()
    ratio = g.used / f.used
    print(f"earlier rule (--windward-loads {rule}, --size-feedback {size_feedback}, --leeward-test {leeward_test}; the state "
          f"saved size_feedback {saved['size_feedback']}, leeward_test {saved['leeward_test']}): 30-degree cap fit "
          f"{g.fitted * 1e3:.2f} mm over {g.points} corners, concave {g.concave}, bound {g.bound * 1e3:.2f} mm, "
          f"radius used {g.used * 1e3:.2f} mm")
    print(f"radius used, normals over sphere {ratio:.4f}: under sphere the stagnation heating is {np.sqrt(ratio):.4f} times and "
          f"the pressure gradient's scale {ratio:.4f} times the earlier rule's")
    print(f"windward under normals and sheltered under sphere {int((windward_normals & f.sheltered).sum())}; leeward under "
          f"normals and exposed under sphere {int((~windward_normals & f.exposed).sum())}, of which windward under sphere "
          f"(heated) {int((~windward_normals & f.exposed & windward_sphere).sum())}")


if __name__ == "__main__":
    for p in (sys.argv[1:] or DEFAULTS):
        check(p)
