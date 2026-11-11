import sys, time; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.transfer import interpolate
from lsproto2.cutthermal import CutThermal
from lsproto2.simplemat import ConstMaterial
mat = ConstMaterial()
h = 1e-3
m = build(BoxBase.around((-3.37e-3, -3.37e-3, -8e-3), (3.37e-3, 3.37e-3, 3e-3), h))
x, y, z = m.points.T
for ghost in (0.0, 0.1):
    s = CutThermal(mat, ghost=ghost, linear_solver='direct')
    s.setup(m, z - 0.0); s.set_temperature(300.0)
    probe = np.array([[0.2e-3, 0.1e-3, -3.0e-3]])
    Tp, Ts, conds = [], [], []
    top = 0.0
    for k in range(80):                       # recede 0.05 h per step under a steady 1 MW/m2 for 0.02 s steps
        s.step(0.02, 1e6 * np.maximum(s.P.normal[:, 2], 0))
        top -= 0.05 * h
        s.move(z - top)
        Tp.append(interpolate(m, s.T, probe)[0])
    d2 = np.abs(np.diff(Tp, 2))
    dT = np.diff(Tp); print('  per-step rise median %.3f K, total %.1f K, spike steps %s' % (np.median(np.abs(dT)), Tp[-1]-300, list(np.flatnonzero(d2 > 5*np.median(d2))))); print('ghost %.1f: probe T second differences max %.4f K median %.4f K (ratio %.1f)' % (ghost, d2.max(), np.median(d2), d2.max()/np.median(d2)))
# remesh: refine near the surface, carry T, check energy residual before the offset
s = CutThermal(mat, ghost=0.1, linear_solver='direct')
phi = z + 0.3e-3
s.setup(m, phi); s.set_temperature(300 + 2e4 * (z + 8e-3))
E0 = s.energy()
m2 = build(BoxBase.around((-3.37e-3, -3.37e-3, -8e-3), (3.37e-3, 3.37e-3, 3e-3), h), lambda c, sz: np.where(np.abs(c[:, 2] + 0.3e-3) < 1.5e-3, 0.25e-3, np.inf))
phi2 = m2.points[:, 2] + 0.3e-3
r0 = s.remesh(m2, phi2)
print('remesh: cells %d -> %d, energy residual before offset %.2e, after %.2e, metal volume %.6e -> %.6e' % (len(m.tets), len(m2.tets), r0, s.energy()/E0 - 1, 0, s.metal_volume()))
