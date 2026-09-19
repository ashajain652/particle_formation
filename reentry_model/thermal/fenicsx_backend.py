"""FEniCSx (dolfinx >= 0.9) backend: the scheme of skfem_backend written in UFL on the same mesh (spec section 8).

Per Newton iterate the linearised, symmetric positive definite system
    a(u, v) = rho c u v / dt + k(T_k) grad u . grad v + 4 eps sigma T_k^3 u v |_Gamma
    L(v)    = rho c T_old v / dt + eps sigma (3 T_k^4 + T_amb^4) v |_Gamma  (+ the nodal convective loads)
is assembled with dolfinx and solved with PETSc (CG + hypre BoomerAMG, or LU); k and the secant heat capacity
c = [h(T_k) - h(T_old)] / (T_k - T_old) are DG0 cell coefficients refreshed every iteration, exactly as in the skfem
backend, so both backends discretise the volume terms identically; the radiation term is integrated by quadrature
on nodal T (skfem uses the facet mean) and the convective load is the same nodal vector A_f/3 per facet node.
dolfinx renumbers vertices: `node_of_dof`/`dof_of_node` map between the VolumeMesh's node ids and the P1 dofs.
Serial by default; written for serial runs against dolfinx 0.10 (the nodal-load addition and the dof/node maps
assume one process; an MPI version would scatter the loads and map ghost dofs); untested until fenicsx_env exists.
Lumped mass is not implemented in this backend. `dolfinx` is imported lazily: the constructor raises MissingBackend
without it."""
import numpy as np

from . import SIGMA_SB, MissingBackend, StepResult


class FenicsxThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=False, newton_tol=1e-6, max_iterations=30, cg_tol=1e-10, **_):
        try:
            import dolfinx  # noqa: F401
            import ufl  # noqa: F401
            from mpi4py import MPI  # noqa: F401
            from petsc4py import PETSc  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        if lumped_mass:
            raise ValueError("lumped mass is not implemented in the fenicsx backend")
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.newton_tol, self.max_iterations, self.cg_tol = linear_solver, newton_tol, max_iterations, cg_tol

    def setup(self, mesh, material, emissivity):
        import basix.ufl
        import ufl
        from dolfinx import fem
        from dolfinx import mesh as dmesh
        from mpi4py import MPI
        from petsc4py import PETSc
        from scipy.spatial import cKDTree
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
        self.msh = dmesh.create_mesh(MPI.COMM_WORLD, mesh.tets.astype(np.int64), mesh.points, domain)
        self.V = fem.functionspace(self.msh, ("Lagrange", 1))
        self.V0 = fem.functionspace(self.msh, ("DG", 0))
        n_local = self.V.dofmap.index_map.size_local
        _, self.node_of_dof = cKDTree(mesh.points).query(self.V.tabulate_dof_coordinates()[:n_local])
        self.dof_of_node = np.empty(n_local, dtype=np.int64)
        self.dof_of_node[self.node_of_dof] = np.arange(n_local)
        self.cell_dofs = np.asarray(self.V.dofmap.list)[:, :4]
        x = self.V.tabulate_dof_coordinates()[self.cell_dofs]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        self.cell_volumes = np.abs(np.linalg.det(J)) / 6.0
        surface = mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        self.T, self.T_old, self.T_k = fem.Function(self.V), fem.Function(self.V), fem.Function(self.V)
        self.k_fun, self.c_fun = fem.Function(self.V0), fem.Function(self.V0)
        self.dt_c, self.T_amb_c = fem.Constant(self.msh, PETSc.ScalarType(1.0)), fem.Constant(self.msh, PETSc.ScalarType(0.0))
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        es = self.emissivity * SIGMA_SB
        self.a = fem.form(self.c_fun / self.dt_c * u * v * ufl.dx + self.k_fun * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
                          + 4.0 * es * self.T_k ** 3 * u * v * ufl.ds)
        self.L = fem.form(self.c_fun / self.dt_c * self.T_old * v * ufl.dx + es * (3.0 * self.T_k ** 4 + self.T_amb_c ** 4) * v * ufl.ds)
        self.rad_form = fem.form(es * (self.T ** 4 - self.T_amb_c ** 4) * ufl.ds)
        self.ksp = PETSc.KSP().create(self.msh.comm)
        if self.linear_solver == "direct":
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.setType("cg")
            self.ksp.getPC().setType("hypre")
            self.ksp.getPC().setHYPREType("boomeramg")
            self.ksp.setTolerances(rtol=self.cg_tol, max_it=500)
        self.T.x.array[:] = 300.0
        self._T_prev = None

    def set_temperature(self, T):
        self.T.x.array[:] = float(T) if np.ndim(T) == 0 else np.asarray(T, dtype=float)[self.node_of_dof]
        self.T.x.scatter_forward()
        self._T_prev = None

    def temperature(self):
        return self.T.x.array[self.dof_of_node].copy()

    def facet_load(self, q):
        """Nodal loads A_f/3 per facet node, in the VolumeMesh's node order."""
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=len(self.mesh.points))

    def _coefficients(self, T_dofs, T_old_dofs):
        Te, Te_old = T_dofs[self.cell_dofs].mean(axis=1), T_old_dofs[self.cell_dofs].mean(axis=1)
        dT = Te - Te_old
        moved = np.abs(dT) > 1e-9
        c = np.where(moved, (self.material.enthalpy(Te) - self.material.enthalpy(Te_old)) / np.where(moved, dT, 1.0), self.material.cp_eff(Te))
        self.k_fun.x.array[:] = self.material.k(Te)
        self.c_fun.x.array[:] = self.material.rho * c

    def energy(self):
        Te = self.T.x.array[self.cell_dofs].mean(axis=1)
        return float((self.material.rho * self.cell_volumes * self.material.enthalpy(Te)).sum())

    def radiated_power(self, T_amb):
        from dolfinx import fem
        from mpi4py import MPI
        self.T_amb_c.value = T_amb
        return float(self.msh.comm.allreduce(fem.assemble_scalar(self.rad_form), op=MPI.SUM))

    def step(self, dt, q_conv, T_amb, dirichlet=None):
        from dolfinx import fem
        from dolfinx.fem.petsc import apply_lifting, assemble_matrix, assemble_vector, set_bc
        from petsc4py import PETSc
        self.dt_c.value, self.T_amb_c.value = dt, T_amb
        self.T_old.x.array[:] = self.T.x.array
        T_old = self.T.x.array.copy()
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))[self.node_of_dof]
        bcs = []
        if dirichlet is not None:
            g = fem.Function(self.V)
            dofs = self.dof_of_node[np.asarray(dirichlet[0])]
            g.x.array[dofs] = dirichlet[1]
            bcs = [fem.dirichletbc(g, dofs.astype(np.int32))]
        T_new, iteration = fem.Function(self.V), 0
        converged = False
        last_relative_change = None
        for iteration in range(1, self.max_iterations + 1):
            self.T_k.x.array[:] = T_k
            self._coefficients(T_k, T_old)
            A = assemble_matrix(self.a, bcs=bcs)
            A.assemble()
            b = assemble_vector(self.L)
            b.array[:] += F_conv
            apply_lifting(b, [self.a], bcs=[bcs])
            b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            set_bc(b, bcs)
            self.ksp.setOperators(A)
            self.ksp.solve(b, T_new.x.petsc_vec)
            T_new.x.scatter_forward()
            last_relative_change = np.linalg.norm(T_new.x.array - T_k) / np.linalg.norm(T_new.x.array)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            T_k = T_new.x.array.copy()
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev = T_old
        self.T.x.array[:] = T_new.x.array
        return StepResult(self.temperature(), float(F_conv.sum()), self.radiated_power(T_amb), iteration)
