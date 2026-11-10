"""Stage 1 measurement: build the adaptive mesh around a sphere and time each geometric operation.

Usage: "$PY" analysis/levelset/measure_geometry.py <diameter_mm> [base_mm] [fine_mm]. Box: the sphere plus 25 mm each side
(the design's 150 mm box for the 100 mm sphere). Fine region: metal within margin d = 4 mm (+0.6 mm lead) of
the surface (no melt, so psi = phi). Prints one JSON line."""
import json
import os
import resource
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO_ROOT)

import numpy as np  # noqa: E402

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.criteria import kappa_max_nodal, make_target
from reentry_model.levelset.cutcell import inside_volumes, patches
from reentry_model.levelset.levelset import correct_volumes, pieces, redistance
from reentry_model.levelset.shadow import shadowed
from reentry_model.levelset.thin import patch_measures

D = float(sys.argv[1]) * 1e-3
H0 = float(sys.argv[2]) * 1e-3 if len(sys.argv) > 2 else 2.8e-3
FINE = float(sys.argv[3]) * 1e-3 if len(sys.argv) > 3 else 0.7e-3
R = D / 2
c = np.full(3, R + 25e-3)
base = BoxBase.around((0, 0, 0), 2 * c, H0)
out = {"diameter_mm": D * 1e3, "base_mm": H0 * 1e3, "fine_mm": FINE * 1e3, "box_mm": base.hi[0] * 1e3}
t = {}


def tick(name, f):
    t0 = time.perf_counter()
    r = f()
    t[name] = round(time.perf_counter() - t0, 2)
    return r


sdf = lambda x: np.linalg.norm(x - c, axis=1) - R
first = lambda x, s: np.where((sdf(x) < 0.94 * s) & (sdf(x) > -(4.6e-3 + 0.94 * s)), FINE, np.inf)
m = tick("build_from_analytic_s", lambda: build(base, first))
BAND = 4.6e-3 + 2 * H0                  # margin + lead + two base cells
phi = tick("redistance_s", lambda: redistance(m, sdf(m.points), band=BAND))
P = tick("patches_s", lambda: patches(m, phi))
label, n = tick("pieces_s", lambda: pieces(m, phi))
phi2, worst = tick("volume_correction_s", lambda: correct_volumes(m, phi, label, [4 / 3 * np.pi * R ** 3]))
km = tick("curvature_s", lambda: kappa_max_nodal(m, phi, band=2 * 0.94 * H0))
tgt = make_target(m, phi, phi, km, km, fine=FINE)
m2 = tick("rebuild_from_fields_s", lambda: build(base, tgt))
hit = tick("shadow_s", lambda: shadowed(m, phi, P, (1.0, 0.0, 0.0)))
psi = np.maximum(redistance(m, sdf(m.points) + 0.3e-3, band=BAND), phi)
h = tick("thin_measures_s", lambda: patch_measures(m, phi, psi, P))
out.update(times=t, cells=int(len(m.tets)), nodes=int(len(m.points)), cells_rebuilt=int(len(m2.tets)),
           cells_at_fine=int((m.size <= FINE + 1e-12).sum()), patches=int(len(P.cell)), pieces=int(n),
           volume_error_before_correction=float(worst), shadowed_share=float(hit.mean()),
           area_error=float(P.area.sum() / (4 * np.pi * R ** 2) - 1),
           peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2))
print(json.dumps(out))
