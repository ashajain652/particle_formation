import sys, time, resource; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from reentry_model.levelset.boxmesh import BoxBase, build
from lsproto2.cutthermal import CutThermal
import material as mm
D = float(sys.argv[1]); R = D/2
mat = mm.Material.from_drama_json(mm.MATERIAL_NAMES['AA7075_scheil'])
L = 0.75*D
base = BoxBase.around((-L,)*3, (L,)*3, 2.8e-3)
def target(c, sz):
    r = np.linalg.norm(c, axis=1)
    return np.where((r > R - 4.6e-3) & (r < R + 1.5e-3), 0.7e-3, np.inf)
t = time.time(); m = build(base, target); tb = time.time() - t
phi = np.linalg.norm(m.points, axis=1) - R
s = CutThermal(mat, emissivity=0.1, linear_solver='amg', newton_tol=1e-8)
t = time.time(); s.setup(m, phi); ts = time.time() - t
s.set_temperature(300.0)
q0 = 2e6 * np.maximum(-s.P.normal[:, 0], 0)          # flow along +x: windward face heated
times = []
for k in range(3):
    t = time.time(); s.step(0.5, q0, rad=1.0, T_amb=250.0); times.append((time.time() - t, s.iterations))
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
print('D %.0f mm: cells %d active %d nodes %d | build %.1f s setup %.1f s | steps %s | peak %.2f GB' % (
    D*1e3, len(m.tets), len(s.cells), int(s.active.sum()), tb, ts, ['%.1f s/%d it' % x for x in times], rss), flush=True)
