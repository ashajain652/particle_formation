import sys, time; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from reentry_model.levelset.boxmesh import BoxBase, build
from lsproto2.cutthermal import CutThermal
import material as mm
mat = mm.Material.from_drama_json(mm.MATERIAL_NAMES['AA7075_nomelt'])
R = 10e-3
base = BoxBase.around((-14e-3,)*3, (14e-3,)*3, 2e-3)
m = build(base)
phi = np.linalg.norm(m.points, axis=1) - R
s = CutThermal(mat, emissivity=0.0, linear_solver='direct')
t=time.time(); s.setup(m, phi); print('setup', time.time()-t, 'cells', len(m.tets))
V = 4/3*np.pi*R**3
print('metal vol rel err', s.metal_volume()/V-1, 'mass sum', s.Me.sum()/s.Vin.sum()-1, 'patch area rel', s.P.area.sum()/(4*np.pi*R**2)-1)
s.set_temperature(300.0)
E0 = s.energy(); Q=0
for i in range(5):
    Q += s.step(0.1, 1e6)
print('energy books rel', (s.energy()-E0-Q)/Q, 'iters', s.iterations, 'Tmin/max', s.T[s.active].min(), s.T[s.active].max())
