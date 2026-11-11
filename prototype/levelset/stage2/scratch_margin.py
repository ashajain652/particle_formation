"""Plate column on an adaptive mesh that follows psi: fine (0.7 mm) where psi > -d, coarse (2.8 mm) elsewhere."""
import sys, time, json, argparse; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np
from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.transfer import interpolate
from lsproto2.cutthermal import CutThermal
import material as mm
p = argparse.ArgumentParser(); p.add_argument('--d', type=float, default=4e-3); p.add_argument('--dt', type=float, default=0.05)
p.add_argument('--T0', type=float, default=700.); p.add_argument('--q', type=float, default=6e6); p.add_argument('--t-end', type=float, default=2.0)
a = p.parse_args()
mat = mm.Material.from_drama_json(mm.MATERIAL_NAMES['AA7075_scheil'])
H, h = 2.8e-3, 0.7e-3; W = 2*H; D = 8e-3; off = 0.0
base = BoxBase.around((-W/2-H, -W/2-H, -D-2*H+0.37*h), (W/2+H, W/2+H, H+0.37*h), H)
sx = W/2 + off
def plate_phi(pts, top):
    x, y, z = pts.T
    return np.maximum.reduce([z - top, np.abs(x) - sx, np.abs(y) - sx, -z - D - 0.5*h])
def mesh_for(z_rigid):
    return build(base, lambda c, sz: np.where(c[:, 2] > z_rigid - a.d - 0.5*h, h, np.inf))
line = np.column_stack([np.zeros(4001), np.zeros(4001), np.linspace(0, -D, 4001)])
def fronts(s):
    ph = interpolate(s.mesh, s.phi, line); zf = line[np.argmax(ph < 0), 2]
    Tl = interpolate(s.mesh, s.T, line); cold = np.flatnonzero((ph < 0) & (Tl <= mat.T_rigid))
    return zf, (line[cold[0], 2] if len(cold) else zf)
def zf_exact(s):
    up = s.P.normal[:, 2] > 0.99
    return (s.P.centroid[up, 2] * s.P.area[up]).sum() / s.P.area[up].sum()
m = mesh_for(0.0)
s = CutThermal(mat, ghost=0.1, linear_solver='direct'); s.setup(m, plate_phi(m.points, 0.0)); s.set_temperature(a.T0)
E0 = s.energy(); Qin = Eout = 0.0; out = []; resid = []; tic = time.time()
for k in range(1, int(round(a.t_end/a.dt)) + 1):
    Qin += s.step(a.dt, a.q * np.maximum(s.P.normal[:, 2], 0.0))
    g = redistance(s.mesh, s.T - mat.T_feed)
    phin = np.maximum(s.phi, g)
    if (phin > s.phi + 1e-15).any():
        Eout += s.move(phin)
    zf, zr = fronts(s)
    m2 = mesh_for(zr)
    if len(m2.tets) != len(s.mesh.tets) or not np.array_equal(m2.level, s.mesh.level):
        ph2 = interpolate(s.mesh, s.phi, m2.points)
        resid.append(s.remesh(m2, ph2))
    out.append((k*a.dt, -zf*1e3, (zf - zr)*1e3))
books = (s.energy() - E0 - Qin + Eout) / Qin
print('d %.1f mm dt %.2f T0 %d: books %.1e, remeshes %d, worst energy residual before offset %.1e, wall %.0f s, cells now %d' % (
    a.d*1e3, a.dt, a.T0, books, len(resid), max(np.abs(resid)) if resid else 0, time.time()-tic, len(s.mesh.tets)))
for t_ in (1.0, 1.5, 2.0):
    r = min(out, key=lambda o: abs(o[0]-t_)); print('  t %.2f: recession %.4f mm  non-rigid %.4f mm' % r)
