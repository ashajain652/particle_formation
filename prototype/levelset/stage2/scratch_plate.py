"""Level-2 plate centre as a column with adiabatic sides: recession and non-rigid depth vs the 1D reference."""
import sys, time, json; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.transfer import interpolate
from lsproto2.cutthermal import CutThermal
import material as mm
import argparse
p = argparse.ArgumentParser(); p.add_argument('--h', type=float, default=0.7e-3); p.add_argument('--dt', type=float, default=0.05)
p.add_argument('--T0', type=float, default=700.); p.add_argument('--q', type=float, default=6e6); p.add_argument('--t-end', type=float, default=2.0)
p.add_argument('--off', type=float, default=0.37); p.add_argument('--top-off', type=float, default=0.0)
a = p.parse_args()
mat = mm.Material.from_drama_json(mm.MATERIAL_NAMES['AA7075_scheil'])
h = a.h; W = 6*h; D = 8e-3; top = a.top_off*h
lo = np.array([-W/2-h*a.off-h, -W/2-h*a.off-h, -D-2*h]); hi = -lo; hi[2] = 2*h
m = build(BoxBase.around(lo, hi, h))
x, y, z = m.points.T; sx = W/2 + h*a.off
phi = np.maximum.reduce([z - top, np.abs(x) - sx, np.abs(y) - sx, -z - D + top])
s = CutThermal(mat, ghost=0.1, linear_solver='direct'); s.setup(m, phi); s.set_temperature(a.T0)
E0 = s.energy(); Qin = 0.0; Eout = 0.0
line = np.column_stack([np.zeros(4001), np.zeros(4001), np.linspace(top, top - D, 4001)])
t = 0.0; out = []; tic = time.time()
nsteps = int(round(a.t_end / a.dt))
for k in range(1, nsteps + 1):
    up = np.maximum(s.P.normal[:, 2], 0.0)
    Qin += s.step(a.dt, a.q * up); t = k * a.dt
    g = redistance(m, s.T - mat.T_feed)
    phin = np.maximum(s.phi, g)
    if (phin > s.phi + 1e-15).any():
        Eout += s.move(phin)
    # front: phi zero along the centre line; non-rigid: T = T_rigid below it
    ph = interpolate(m, s.phi, line); zf = line[np.argmax(ph < 0), 2]
    Tl = interpolate(m, s.T, line); inside = ph < 0
    cold = np.flatnonzero(inside & (Tl <= mat.T_rigid))
    zr = line[cold[0], 2] if len(cold) else zf
    out.append((t, (top - zf) * 1e3, (zf - zr) * 1e3))
Tn = s.T; act = s.active
zz = m.points[:,2]; c = (np.abs(m.points[:,0])<1e-9)&(np.abs(m.points[:,1])<1e-9)
for i in np.argsort(-zz[c])[:14]:
    j=np.flatnonzero(c)[i]; print('  z %.2f mm T %.1f phi %.3f active %s' % (zz[j]*1e3, Tn[j], s.phi[j]*1e3, act[j]))
print('  Qin %.1f J (expected %.1f)  Eout %.1f J  stored %.1f J' % (Qin, a.q*a.t_end*(W+2*h*a.off)**2, Eout, s.energy()-E0))
books = (s.energy() - E0 - Qin + Eout) / Qin
print('h %.2f mm dt %.3f s  books %.1e  wall %.0f s  cells %d' % (h*1e3, a.dt, books, time.time()-tic, len(m.tets)))
for t_ in (0.5, 1.0, 1.5, 2.0):
    r = min(out, key=lambda o: abs(o[0]-t_))
    print('  t %.2f: recession %.4f mm  non-rigid %.4f mm' % r)
json.dump(out, open('plate_h%.2f_dt%.3f_T%d.json' % (h*1e3, a.dt, a.T0), 'w'))
