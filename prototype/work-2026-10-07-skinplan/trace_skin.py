import numpy as np
from reentry_model import material
mat = material.Material.from_drama_json("AA7075_scheil")
T = np.array([892.299, 891.624, 890.916, 890.18, 889.417, 888.631, 887.823, 886.995, 886.149, 885.287, 884.41, 883.519,
              882.615, 881.7, 880.773, 879.837, 878.891, 877.936, 876.973, 876.003, 875.026, 874.041, 873.051, 872.055,
              871.053, 870.046, 869.034, 868.017, 866.996, 865.971, 864.942, 863.909, 862.872, 861.832, 860.789, 859.743,
              858.694, 857.642, 856.587, 855.53])
src = open('reentry_model/skin.py').read()
src = src.replace("            T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try\n",
                  "            T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try\n"
                  "            j = int(np.abs(R).argmax()); print('lam', lam, 'moved', moved, 'maxR', np.abs(R).max(axis=1), 'cell', j, 'T', round(float(T[0, j]), 4))\n")
src = src.replace("from .thermal import SIGMA_SB", "from reentry_model.thermal import SIGMA_SB")
ns = {}
exec(compile(src, 'skin_dbg', 'exec'), ns)
g = ns['SkinField'](mat, ns['SkinSettings'](cell=10e-6))
g.add(np.array([0]), np.array([3e-6]), T[None, :], np.array([0.4e-3]))
try:
    g.conduct(0.01, np.array([8205094.4]), 0.0, np.array([0.4]), np.array([855.0]), np.array([0.0]))
except RuntimeError as e:
    print(e)
