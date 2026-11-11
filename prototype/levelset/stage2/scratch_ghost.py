import sys; sys.path[:0] = ['.', '../../..', 'ext']
import numpy as np, scipy.sparse.linalg as spla, pyamg
from reentry_model.levelset.boxmesh import BoxBase, build
from lsproto2.cutthermal import CutThermal
from lsproto2.simplemat import ConstMaterial
mat = ConstMaterial()
m = build(BoxBase.around((-14e-3,)*3, (14e-3,)*3, 1e-3))
for shift in (0.0, 0.999e-3*0.5, 1e-3*0.001):
    for R in (10e-3 + shift,):
        phi = np.linalg.norm(m.points, axis=1) - R
        for g in (0.0, 0.01, 0.03, 0.1):
            s = CutThermal(mat, ghost=g); s.setup(m, phi); s.set_temperature(300.)
            T = s.T; K, M = s.operators(T, T); idx = np.flatnonzero(s.active)
            A = (M/1e-3 + K).tocsr()[idx][:, idx]
            vfrac = s.Vin / m.volumes[s.cells]
            d = A.diagonal(); b = np.random.default_rng(0).standard_normal(len(idx))
            ml = pyamg.smoothed_aggregation_solver(A, symmetry='symmetric'); it = []
            spla.cg(A, b, rtol=1e-12, maxiter=2000, M=ml.aspreconditioner(), callback=lambda _: it.append(1))
            it2 = []; spla.cg(A, b, rtol=1e-12, maxiter=5000, callback=lambda _: it2.append(1))
            print('R %.4f mm ghost %.1f: smallest metal fraction %.1e, diag ratio %.1e, AMG-CG iters %d, plain CG iters %d' % (R*1e3, g, vfrac.min(), d.min()/d.max(), len(it), len(it2)))
