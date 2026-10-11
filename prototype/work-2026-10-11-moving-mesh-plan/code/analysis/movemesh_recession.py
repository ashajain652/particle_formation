"""Sub-plan 19 measurement: the moving 3D mesh on the default 100 mm sphere through a prescribed recession of the
front, with heating, to measure rebuild count, rebuild time, energy residuals, element quality and the cost of a
macro step on moving geometry.

Every corner recedes along its inward normal by NOSE_STEP x max(cos theta, 0) per 0.5 s macro step (theta its angle
from the flight direction, -x), so the nose goes 0.5 mm per step, a speed of 1 mm/s, as on the 100 mm flight at
50-60 s (fact: 16.6 mm at 50 s, 29.0 mm at 60 s). The front is heated at Q_STAG x max(cos theta, 0) (no melting
material, `AA7075_nomelt`) so the field the remap and the rebuilds carry has the flight's steep near-surface gradient.
The run stops at NOSE_TOTAL of nose recession.

    python analysis/movemesh_recession.py [--nose-total-mm 30] [--solver skfem|fenicsx] [--out PATH.json]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import material, mesh, movemesh  # noqa: E402
from reentry_model.thermal import thermal_solver  # noqa: E402

DT, NOSE_STEP, Q_STAG = 0.5, 0.5e-3, 1.0e6

parser = argparse.ArgumentParser()
parser.add_argument("--nose-total-mm", type=float, default=30.0)
parser.add_argument("--solver", default="skfem")
parser.add_argument("--out", default="movemesh_recession.json")
parser.add_argument("--no-smoothing", action="store_true", help="corners move along their normals only")
parser.add_argument("--smooth-lambda", type=float, default=None)
parser.add_argument("--smooth-clamp", type=float, default=None)
args = parser.parse_args()

mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_nomelt"])
m = mesh.sphere_mesh(0.05)
solver = thermal_solver(args.solver)
solver.setup(m, mat, mat.emissivity)
solver.set_temperature(300.0)
moving = movemesh.MovingMesh(m, solver, mat)
flight = np.array([-1.0, 0.0, 0.0])
rows, t_total, failure = [], time.perf_counter(), None
steps = int(round(args.nose_total_mm * 1e-3 / NOSE_STEP))
for k in range(steps):
    t0 = time.perf_counter()
    normal, _ = moving.frame()
    cos = np.clip(normal @ flight, 0.0, None)
    u = -NOSE_STEP * cos[:, None] * normal
    if not args.no_smoothing:
        u = u + moving.smoothing(normal, args.smooth_lambda, args.smooth_clamp)
    try:
        result = moving.move_with_policy(u)
    except RuntimeError as exc:                       # a move no rebuild can absorb: record where, and stop
        failure = {"failed_at_mm": (k + 1) * NOSE_STEP * 1e3, "error": str(exc)}
        print(json.dumps(failure), flush=True)
        break
    t_move = time.perf_counter() - t0
    surface = moving.mesh.surface()
    q = Q_STAG * np.clip(surface.normals @ flight, 0.0, None)
    t1 = time.perf_counter()
    step = solver.step(DT, q, 300.0)
    t_step = time.perf_counter() - t1
    ratio = movemesh.radius_ratio(moving.mesh.points, moving.mesh.tets)
    row = {"step": k + 1, "nose_mm": (k + 1) * NOSE_STEP * 1e3, "rebuilt": result.rebuilt,
           "rebuilds": moving.rebuilds, "correction_J": result.correction,
           "swept_J": float(result.swept_enthalpy.sum()), "worst_ratio": float(ratio.min()),
           "worst_ratio_moved": result.worst_ratio,
           "surface_quality": moving.surface_quality(), "surface_fold_deg": moving.surface_fold(), "move_s": t_move, "step_s": t_step,
           "T_max": float(step.T.max()), "n_nodes": moving.mesh.n_nodes, "n_elements": moving.mesh.n_elements}
    if result.rebuilt:
        rb = moving.last_rebuild
        row.update({"rebuild_s": rb.seconds, "rebuild_residual_J": rb.residual,
                    "rebuild_relative_residual": rb.relative_residual, "rebuild_reason": rb.reason})
    rows.append(row)
    if (k + 1) % 10 == 0 or result.rebuilt:
        print(json.dumps({key: (round(v, 6) if isinstance(v, float) else v) for key, v in row.items()}), flush=True)
summary = {"failure": failure, "completed_mm": rows[-1]["nose_mm"] if rows else 0.0, "solver": args.solver, "steps": steps, "nose_total_mm": args.nose_total_mm, "rebuilds": moving.rebuilds,
           "wall_s": time.perf_counter() - t_total,
           "rebuild_s": [r["rebuild_s"] for r in rows if "rebuild_s" in r],
           "rebuild_relative_residuals": [r["rebuild_relative_residual"] for r in rows if "rebuild_s" in r],
           "max_abs_correction_over_swept": max(abs(r["correction_J"]) / r["swept_J"] for r in rows if r["swept_J"] > 0),
           "worst_ratio_min": min(r["worst_ratio"] for r in rows),
           "surface_quality_min": min(r["surface_quality"] for r in rows), "smoothing": not args.no_smoothing,
           "move_s_median": float(np.median([r["move_s"] for r in rows])),
           "step_s_median": float(np.median([r["step_s"] for r in rows]))}
print(json.dumps(summary, indent=1))
json.dump({"summary": summary, "rows": rows}, open(args.out, "w"), indent=1)
