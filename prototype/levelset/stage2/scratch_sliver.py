import sys; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.transfer import interpolate
from lsproto2.cutthermal import CutThermal
from lsproto2.simplemat import ConstMaterial
h = 1e-3
m = build(BoxBase.around((-3.37e-3, -3.37e-3, -8e-3), (3.37e-3, 3.37e-3, 3e-3), h))
z = m.points[:, 2]
for g in (0.01, 0.03):
    s = CutThermal(ConstMaterial(), ghost=g, linear_solver="direct"); s.setup(m, z); s.set_temperature(300.0)
    probe, Tp, top, Q, A = np.array([[0.2e-3, 0.1e-3, -3.0e-3]]), [], 0.0, [], []
    for k in range(60):
        Q.append(s.step(0.02, 1e6 * np.maximum(s.P.normal[:, 2], 0)))
        A.append(s.P.area[s.P.normal[:,2]>0.99].sum())
        top -= 0.05 * h
        s.move(z - top)
        Tp.append(interpolate(m, s.T, probe)[0])
    d2 = np.abs(np.diff(Tp, 2))
    bad = np.flatnonzero(d2[3:] > 3*np.median(d2[3:])) + 3
    print('ghost', g, 'spikes at steps', bad, 'top then (mm)', [round(-(k+1)*0.05,2) for k in bad], 'd2', np.round(d2[bad],3), 'Q', np.round(np.array(Q)[bad[:3]],5) if len(bad) else '', 'A', np.round(np.array(A)[[0,10,19,20,21]]*1e6,4))
