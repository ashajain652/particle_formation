"""Stage 2 measurement: the cost of the cut-cell heat solve (`CutThermal`) on the adaptive mesh around a sphere.

Usage: "$PY" analysis/levelset/measure_thermal.py [diameter_mm] [steps] (default 50 mm, 3 steps). The mesh is the one
of the plan's fact 10: a box of 1.5 diameters a side, 2.8 mm base cells, 0.7 mm cells from 4.6 mm inside the surface
to 1.5 mm outside it. The sphere starts at 300 K; each 0.5 s step heats the windward half (flow along +x) by
2 MW/m2 times the projected-area factor max(-n_x, 0) and radiates with emissivity 0.1 to 250 K. The material is this
branch's AA7075 without melting (`reentry_model.material`): the surface stays below about 500 K over three steps, far
from the 750 K solidus, so the Scheil material of the prototype would give the same solve. Prints one JSON line and
writes it to reentry_model_output/levelset/measure_thermal_<D>mm.json."""
import json
import os
import resource
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO_ROOT)

import numpy as np  # noqa: E402

from reentry_model.levelset.boxmesh import BoxBase, build  # noqa: E402
from reentry_model.levelset.cutthermal import CutThermal  # noqa: E402
from reentry_model.material import Material  # noqa: E402

D = float(sys.argv[1]) * 1e-3 if len(sys.argv) > 1 else 50e-3
STEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
H0, FINE, DT = 2.8e-3, 0.7e-3, 0.5
R, L = D / 2, 0.75 * D
base = BoxBase.around((-L,) * 3, (L,) * 3, H0)
out = {"diameter_mm": D * 1e3, "base_mm": H0 * 1e3, "fine_mm": FINE * 1e3, "box_mm": round(2 * L * 1e3, 3), "dt_s": DT}
t = {}


def tick(name, f):
    t0 = time.perf_counter()
    r = f()
    t[name] = round(time.perf_counter() - t0, 2)
    return r


def target(c, s):
    r = np.linalg.norm(c, axis=1)
    return np.where((r > R - 4.6e-3) & (r < R + 1.5e-3), FINE, np.inf)


m = tick("build_s", lambda: build(base, target))
phi = np.linalg.norm(m.points, axis=1) - R
solver = CutThermal(Material.from_drama_json(), emissivity=0.1, linear_solver="amg", newton_tol=1e-8)
tick("setup_s", lambda: solver.setup(m, phi))
solver.set_temperature(300.0)
E0 = solver.energy()
q0 = 2e6 * np.maximum(-solver.P.normal[:, 0], 0.0)      # windward half, by projected area (plan fact 3)
step_s, iterations, heat_in = [], [], 0.0
for k in range(STEPS):
    t0 = time.perf_counter()
    heat_in += solver.step(DT, q0, rad=1.0, T_amb=250.0)
    step_s.append(round(time.perf_counter() - t0, 2))
    iterations.append(int(solver.iterations))
t["step_s"] = step_s
out.update(times=t, newton_iterations=iterations, cells=int(len(m.tets)), active_cells=int(len(solver.cells)),
           active_nodes=int(solver.active.sum()), patches=int(len(solver.P.area)),
           max_surface_T_K=round(float(solver.patch_T().max()), 1),
           energy_books_relative=float(abs(solver.energy() - E0 - heat_in) / abs(heat_in)),
           peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2))
line = json.dumps(out)
print(line)
dest = os.path.join(REPO_ROOT, "reentry_model_output", "levelset")
os.makedirs(dest, exist_ok=True)
with open(os.path.join(dest, "measure_thermal_{:.0f}mm.json".format(D * 1e3)), "w") as fh:
    fh.write(line + "\n")
