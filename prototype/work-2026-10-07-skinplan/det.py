import numpy as np
from reentry_model import material, thermal, mesh
m = mesh.sphere_mesh(0.05, 4e-3, 20e-3, "/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/930f32b6-e08a-42f0-a4ce-595da1b5cc42/scratchpad/skinplan/meshes", band=0.0)
mat = material.Material.from_drama_json("AA7075_scheil")
def mk():
    s = thermal.thermal_solver("skfem"); s.setup(m, mat, 0.3); s.set_temperature(700.0); return s
a, b = mk(), mk()
q = np.full(len(a.areas), 2e5)
ra = a.step(0.5, q, 0.0); rb = b.step(0.5, q, 0.0)
print('plain vs plain identical', np.array_equal(ra.T, rb.T), np.abs(ra.T - rb.T).max())
c = mk(); nf = len(c.areas)
rc = c.step(0.5, q, 0.0, interface=(np.zeros(nf, dtype=bool), np.zeros(nf), np.zeros(nf)))
print('plain vs empty mask', np.array_equal(ra.T, rc.T), np.abs(ra.T - rc.T).max())
