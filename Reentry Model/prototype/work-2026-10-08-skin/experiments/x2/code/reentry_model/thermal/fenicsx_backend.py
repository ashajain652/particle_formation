"""FEniCSx (dolfinx >= 0.11) backend: the scheme of skfem_backend on the same mesh (spec sections 8 and Step 3 6).

dolfinx assembles the stiffness K = int phi_e k(T_e) grad u . grad v dx (a DG0 coefficient refreshed every Newton
iterate); everything nodal -- the lumped capacity diag(sum_e phi_e V_e/4 rho c_p,eff(T_i)), the enthalpy-rate vector
E = diag(sum_e phi_e V_e/4 rho c_sec,i)(T - T_old), the convective loads A_f/3 per facet node, the radiation
eps sigma (T_f^4 - T_amb^4) A_f/3 with its Jacobian 4 eps sigma T_f^3 A_f/9 on the facet node pairs, and the Step 3
nodal loads -- is built in numpy on the VolumeMesh's numbering and added to the PETSc operator (the facet pairs are
edges of cells, so they lie inside K's pattern). Both backends therefore discretise identically; the radiation
follows the current boundary of the active set (element death), which a UFL `ds` measure could not. Dirichlet
nodes (analytic tests) and pinned material-free nodes are imposed with MatZeroRowsColumns. PETSc CG with hypre
BoomerAMG (hierarchy reused for `amg_rebuild_every` solves) or LU. dolfinx renumbers vertices: `node_of_dof` /
`dof_of_node` map between the VolumeMesh's node ids and the P1 dofs. Serial (the nodal vectors assume one
process). Verified with dolfinx 0.11.0 in fenicsx_env (2026-09-20). `dolfinx` is imported lazily: the constructor
raises MissingBackend without it."""
import numpy as np
import scipy.sparse as sp

from . import SIGMA_SB, MissingBackend, StepResult


class FenicsxThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=True, newton_tol=1e-6, max_iterations=30, amg_rebuild_every=30,
                 cg_tol=1e-10, **_):
        try:
            import dolfinx  # noqa: F401
            import ufl  # noqa: F401
            from mpi4py import MPI  # noqa: F401
            from petsc4py import PETSc  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        if not lumped_mass:
            raise ValueError("the fenicsx backend implements the lumped (nodal-enthalpy) capacity only")
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.newton_tol, self.max_iterations, self.cg_tol = linear_solver, newton_tol, max_iterations, cg_tol
        self.amg_rebuild_every, self._solves = amg_rebuild_every, 0
        self.lumped_mass = True

    def setup(self, mesh, material, emissivity):
        import basix.ufl
        import ufl
        from dolfinx import fem
        from dolfinx import mesh as dmesh
        from dolfinx.fem import petsc
        from mpi4py import MPI
        from petsc4py import PETSc
        from scipy.spatial import cKDTree
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
        self.msh = dmesh.create_mesh(MPI.COMM_WORLD, mesh.tets.astype(np.int64), domain, mesh.points)   # dolfinx >= 0.9: (comm, cells, element, x)
        self.V = fem.functionspace(self.msh, ("Lagrange", 1))
        self.V0 = fem.functionspace(self.msh, ("DG", 0))
        n_local = self.V.dofmap.index_map.size_local
        _, self.node_of_dof = cKDTree(mesh.points).query(self.V.tabulate_dof_coordinates()[:n_local])
        self.dof_of_node = np.empty(n_local, dtype=np.int64)
        self.dof_of_node[self.node_of_dof] = np.arange(n_local)
        self.cell_dofs = np.asarray(self.V.dofmap.list)[:, :4]
        # dolfinx cell k is the k-th cell passed in (create_mesh keeps the order in serial): cell volumes and the map
        # from dolfinx cells to mesh elements
        x = self.V.tabulate_dof_coordinates()[self.cell_dofs]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        self.cell_volumes = np.abs(np.linalg.det(J)) / 6.0
        self.cell_nodes = self.node_of_dof[self.cell_dofs]                   # mesh node ids of each dolfinx cell
        centroids = mesh.points[mesh.tets].mean(axis=1)
        self.element_of_cell = cKDTree(centroids).query(mesh.points[self.cell_nodes].mean(axis=1))[1]
        self.coef = fem.Function(self.V0)                                   # phi_e k(T_e) per cell
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        self.a = fem.form(self.coef * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx)
        self.K = petsc.create_matrix(self.a)
        self.ksp = PETSc.KSP().create(self.msh.comm)
        if self.linear_solver == "direct":
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.setType("cg")
            self.ksp.getPC().setType("hypre")
            self.ksp.getPC().setHYPREType("boomeramg")
            self.ksp.setTolerances(rtol=self.cg_tol, max_it=500)
        self.T = np.full(mesh.n_nodes, 300.0)                               # nodal, mesh numbering
        self._T_prev = None
        self.phi = np.where(mesh.active, 1.0, 0.0)
        self.film_mass = np.zeros(mesh.n_nodes)                   # kg per node (see the skfem backend)
        self.pinned = np.zeros(mesh.n_nodes, dtype=bool)
        self._structure_changed = True
        self._refresh_surface()

    # -- geometry / fractions -----------------------------------------------------------------------------------
    def _refresh_surface(self):
        surface = self.mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        n = self.mesh.n_nodes
        rows, cols = np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel()
        self._facet_rows, self._facet_cols = rows, cols
        pinned = np.ones(n, dtype=bool)
        pinned[np.unique(self.mesh.tets[self.phi > 0.0])] = False
        pinned &= self.film_mass <= 0.0                           # a node carrying film keeps a heat capacity
        if not np.array_equal(pinned, self.pinned):
            self._structure_changed = True
        self.pinned = pinned
        self._radiating = None                            # the facets changed: last step's skin mask no longer applies

    def set_film_mass(self, mass):
        """Nodal mass [kg] of the melt film riding on the boundary (see the skfem backend)."""
        mass = np.asarray(mass, dtype=float)
        if mass.shape != (self.mesh.n_nodes,) or (mass < 0.0).any() or not np.isfinite(mass).all():
            raise ValueError("film mass must be one finite non-negative value per node")
        self.film_mass = mass
        self._refresh_surface()

    def set_fractions(self, phi):
        self.phi = np.asarray(phi, dtype=float)
        if self.phi.shape != (self.mesh.n_elements,) or (self.phi < 0.0).any() or (self.phi > 1.0).any():
            raise ValueError("fractions must be one value in [0, 1] per element")
        self._refresh_surface()

    def set_temperature(self, T):
        self.T = np.full(self.mesh.n_nodes, float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    def state(self):
        """What set_state needs to repeat a step: the temperatures, the predictor's previous field, the solve counter."""
        return (self.T.copy(), None if self._T_prev is None else self._T_prev.copy(), self._solves)

    def set_state(self, state):
        self.T = state[0].copy()
        self._T_prev = None if state[1] is None else state[1].copy()
        self._solves = state[2]

    # -- nodal pieces (mesh numbering) -----------------------------------------------------------------------------
    def facet_load(self, q):
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=self.mesh.n_nodes)

    def facet_temperature(self, T=None):
        return (self.T if T is None else T)[self.faces].mean(axis=1)   # same default as the skfem backend: one protocol

    def lumped(self, c_nodal, c_film=None):
        """diag(sum_e phi_e V_e/4 c_i) on the mesh nodes, plus the film's own (liquid) capacity where it rides."""
        vol = self.mesh.element_volumes()
        return np.bincount(self.mesh.tets.ravel(), weights=np.repeat(self.phi * vol / 4.0, 4) * c_nodal[self.mesh.tets].ravel(),
                           minlength=self.mesh.n_nodes) + self.film_mass * (c_nodal if c_film is None else c_film) / self.material.rho

    def element_energies(self, T=None):
        h = self.material.enthalpy(self.T if T is None else T)
        return self.phi * self.material.rho * self.mesh.element_volumes() * h[self.mesh.tets].mean(axis=1)

    def film_weight(self):
        """Share of each node's mass that is melt film, 0 to 1 (see the skfem backend)."""
        solid = np.bincount(self.mesh.tets.ravel(),
                            weights=np.repeat(self.phi * self.mesh.element_volumes() / 4.0, 4) * self.material.rho,
                            minlength=self.mesh.n_nodes)
        total = solid + self.film_mass
        return np.divide(self.film_mass, total, out=np.zeros_like(total), where=total > 0.0)

    def nodal_capacity(self):
        """Heat capacity carried by each node [J/K], material and film together (see the skfem backend)."""
        T = self.temperature()
        c, vol = self.material.cp_eff(T) * self.material.rho, self.mesh.element_volumes()
        return np.bincount(self.mesh.tets.ravel(), weights=np.repeat(self.phi * vol / 4.0, 4) * c[self.mesh.tets].ravel(),
                           minlength=self.mesh.n_nodes) + self.film_mass * self.material.cp(T)

    def energy(self):
        return float(self.element_energies().sum())

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(self.T if T is None else T)
        q = self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas
        radiating = getattr(self, "_radiating", None)
        if radiating is not None and radiating.shape == q.shape:
            q = np.where(radiating, q, 0.0)                    # a facet under a skin radiates from the skin's top
        return float(q.sum())

    def _assemble_stiffness(self, T):
        from dolfinx.fem.petsc import assemble_matrix
        Te = T[self.mesh.tets].mean(axis=1)
        coef = (self.phi > 0.0) * self.material.k(Te)                     # k unscaled by phi (see skfem_backend.operators)
        self.coef.x.array[:] = coef[self.element_of_cell]
        self.K.zeroEntries()
        assemble_matrix(self.K, self.a)
        self.K.assemble()

    def _petsc_from_csr(self, matrix):
        from petsc4py import PETSc
        m = matrix.tocsr()
        m.sum_duplicates()
        return PETSc.Mat().createAIJ(size=m.shape, csr=(m.indptr.astype(np.int32), m.indices.astype(np.int32), m.data), comm=self.msh.comm)

    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None):
        from petsc4py import PETSc
        mat, n = self.material, self.mesh.n_nodes
        T_old = self.T
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        F_extra = np.zeros(n) if nodal_load is None else np.asarray(nodal_load, dtype=float)
        Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
        fixed = self.pinned.copy()
        fixed_values = T_old.copy()
        if dirichlet is not None:
            fixed[np.asarray(dirichlet[0])] = True
            fixed_values[np.asarray(dirichlet[0])] = dirichlet[1]
        vol = self.mesh.element_volumes()
        film_w = self.film_weight() if self.film_mass.any() else np.zeros(n)
        # the skins' base law (Step 3 sub-plan 18): on the facets under a skin the heat flux a + b T_f enters
        # implicitly, linearised like the radiation, and the facet does not radiate (its skin's top does)
        mask = None
        if interface is not None and np.any(interface[0]):
            mask = np.asarray(interface[0], dtype=bool)
            a_int = np.where(mask, np.asarray(interface[1], dtype=float), 0.0)
            b_int = np.where(mask, np.asarray(interface[2], dtype=float), 0.0)
        eps_f = self.emissivity if mask is None else np.where(mask, 0.0, self.emissivity)
        self._radiating = None if mask is None else ~mask
        converged, last_relative_change, previous_change, damping, decreases = False, None, None, 1.0, 0
        T_new, iteration = T_k, 0
        for iteration in range(1, self.max_iterations + 1):
            self._assemble_stiffness(T_k)
            c_tan = mat.cp_eff(T_k)
            dT = T_k - T_old
            moved = np.abs(dT) > 1e-9
            c_sec = np.where(moved, (mat.enthalpy(T_k) - mat.enthalpy(T_old)) / np.where(moved, dT, 1.0), c_tan)
            c_liq = mat.cp(T_k)                                   # the film is liquid: no latent plateau (skfem backend)
            c_sec_l = np.where(moved, (mat.enthalpy_liquid(T_k) - mat.enthalpy_liquid(T_old)) / np.where(moved, dT, 1.0), c_liq)
            M_tan = self.lumped(mat.rho * c_tan, mat.rho * c_liq)
            E = self.lumped(mat.rho * c_sec, mat.rho * c_sec_l) * dT
            Tf = self.facet_temperature(T_k)
            B = sp.csr_matrix((np.repeat(4.0 * eps_f * SIGMA_SB * Tf ** 3 * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
            F_rad = self.facet_load(eps_f * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            F_int = 0.0
            if mask is not None:
                B = B + sp.csr_matrix((np.repeat(-b_int * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
                F_int = self.facet_load(a_int + b_int * Tf)
            b = M_tan * T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra + F_int
            # operator in dof numbering: K (dolfinx) + diag(M_tan/dt) + B
            A = self.K.copy()
            extra = sp.diags(M_tan / dt) + B
            extra = extra.tocsr()[self.node_of_dof][:, self.node_of_dof]
            A.axpy(1.0, self._petsc_from_csr(extra), structure=PETSc.Mat.Structure.SUBSET_NONZERO_PATTERN)
            rhs = PETSc.Vec().createWithArray(b[self.node_of_dof].copy(), comm=self.msh.comm)
            x = PETSc.Vec().createWithArray(fixed_values[self.node_of_dof].copy(), comm=self.msh.comm)
            fixed_dofs = self.dof_of_node[np.flatnonzero(fixed)].astype(np.int32)
            if fixed_dofs.size:
                A.zeroRowsColumns(fixed_dofs, diag=1.0, x=x, b=rhs)
            self.ksp.setOperators(A)
            reuse = self.linear_solver == "amg" and dirichlet is None and not self._structure_changed and self._solves % self.amg_rebuild_every != 0
            self.ksp.getPC().setReusePreconditioner(reuse)
            self._structure_changed = False
            self._solves += 1
            sol = PETSc.Vec().createWithArray(T_k[self.node_of_dof].copy(), comm=self.msh.comm)
            self.ksp.solve(rhs, sol)
            if self.ksp.getConvergedReason() <= 0:
                raise RuntimeError("PETSc KSP did not converge (reason {})".format(self.ksp.getConvergedReason()))
            T_new = sol.getArray()[self.dof_of_node].copy()
            T_new[fixed] = fixed_values[fixed]
            if mat.melts:
                free, w = ~fixed, film_w[~fixed]
                T_new[free] = mat.temperature_from_enthalpy_mixed(
                    mat.enthalpy_mixed(T_k[free], w) + mat.cp_mixed(T_k[free], w) * (T_new[free] - T_k[free]), w)
            last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            if previous_change is not None and last_relative_change > previous_change:
                damping, decreases = max(0.25, 0.5 * damping), 0
            elif damping < 1.0:
                decreases += 1
                if decreases >= 2:
                    damping, decreases = 1.0, 0
            previous_change = last_relative_change
            T_k = T_k + damping * (T_new - T_k) if damping < 1.0 else T_new
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev, self.T = T_old, T_new
        self.last_damping = damping
        Q_int = float(((a_int + b_int * self.facet_temperature(T_new)) * self.areas).sum()) if mask is not None else 0.0
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()),
                          Q_dropped, Q_int)
