import sys, time; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from math import erf, sqrt, pi
from reentry_model.levelset.boxmesh import BoxBase, build
from lsproto2.cutthermal import CutThermal
from lsproto2.simplemat import ConstMaterial
mat = ConstMaterial()
alpha = 150/(2800*900); q=2e6; t_end=0.5
def slab(h, ghost, off=0.37, dt=0.01, lumped=False):
    W = 4e-3; top = 0.0; D = 12e-3
    lo = np.array([-W/2-h*off-h, -W/2-h*off-h, -D-h*off-h]); hi = -lo; hi[2] = 2*h
    m = build(BoxBase.around(lo, hi, h))
    x,y,z = m.points.T
    sx = W/2 + h*off
    phi = np.maximum.reduce([z - top, np.abs(x) - sx, np.abs(y) - sx, -z - D])
    s = CutThermal(mat, ghost=ghost, lumped=lumped, linear_solver="direct"); s.setup(m, phi); s.set_temperature(300.)
    top_patch = s.P.normal[:,2] > 0.9
    t=0
    while t < t_end - 1e-9:
        Q = s.step(dt, q * np.maximum(s.P.normal[:, 2], 0.0)); t += dt
    up = s.P.normal[:, 2] > 0.99
    Ts = (s.patch_T()[up] * s.P.area[up]).sum() / s.P.area[up].sum()
    print('   heat in per step %.4f J, expected %.4f J' % (Q, q*dt*(2*sx)**2), flush=True)
    exact = 300 + 2*q/150*sqrt(alpha*t_end/pi)
    # minimum eigen / conditioning proxy: min diagonal of A
    return Ts, exact, len(m.tets)
for h, dt, l in ((1e-3, 0.01, False), (0.5e-3, 0.01, False), (0.25e-3, 0.01, False)):
    Ts, ex, n = slab(h, 0.1, dt=dt, lumped=l)
    print('h %.2f mm dt %.4f lumped %s: rise err %.2f %% cells %d' % (h*1e3, dt, l, (Ts-ex)/(ex-300)*100, n), flush=True)
