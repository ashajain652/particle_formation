"""Diagnostic (2026-10-07, not part of the prototype): SkfemThermalSolver.step with the final Newton iterate's energy
error split into storage (true nodal enthalpy increment minus the one the linear system enforced), radiation (at the
accepted T minus its linearisation) and the CG residual, appended to DECOMP per step. Generated from proto3's
skfem_backend.step by make_step_decomp (the body of the method is unchanged otherwise)."""
import numpy as np
import scipy.sparse as sp
from reentry_model.thermal import SIGMA_SB, StepResult

DECOMP = []


def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
    """One backward-Euler step with the per-facet convective flux q_conv [W/m2] and an optional nodal load vector
    [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped)."""
    T_old = self.T
    # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
    T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
    F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
    F_extra = np.zeros(len(self.points)) if nodal_load is None else np.asarray(nodal_load, dtype=float)
    Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
    T_new, iteration = T_k, 0
    converged = False
    last_relative_change, previous_change, damping, decreases = None, None, 1.0, 0
    mat = self.material
    film_w = self.film_weight() if self.film_mass.any() else np.zeros(len(self.points))
    for iteration in range(1, self.max_iterations + 1):
        K, M, E = self.operators(T_k, T_old)
        Tf = self.facet_temperature(T_k)
        B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
        F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
        # Newton on R(T) = E(T)/dt + K T - F_conv + F_rad(T) - F_extra with the tangent M/dt + K + B:
        # A T_new = A T_k - R(T_k); with a smooth c_p (c_tan = c_sec) this is the Step 2 iteration unchanged
        A = (M / dt + K + B).tocsr()
        b = M @ T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
        if self.pinned.any():                              # material-free nodes keep their temperature (identity rows)
            A = A + sp.diags(self.pinned.astype(float), format="csr")
            b[self.pinned] = T_old[self.pinned]
        T_new = self._solve_with_dirichlet(A, b, T_k, dirichlet)
        _T_lin, _A, _b, _E, _M, _B, _Frad, _Tk = T_new.copy(), A, b, E, M, B, F_rad, T_k
        if mat.melts:
            # enthalpy-consistent update: the linearised step is an enthalpy increment c_p,eff(T_k) (T_new - T_k)
            # per node; inverting the true h(T) puts a node that would overshoot the melting range where the
            # latent heat actually leaves it (identity where h is linear). Without it nodes jump across the +-2 K
            # ramp of a single-temperature material and the iteration cycles (measured 2026-09-20). The enthalpy
            # inverted is the node's own mixture of material and film (see film_weight), not the material's.
            free, w = ~self.pinned, film_w[~self.pinned]
            T_new[free] = mat.temperature_from_enthalpy_mixed(
                mat.enthalpy_mixed(T_k[free], w) + mat.cp_mixed(T_k[free], w) * (T_new[free] - T_k[free]), w)
        last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
        if last_relative_change <= self.newton_tol:
            converged = True
            break
        # damping when the change grows (a fallback for cycling iterates), released again after two decreases
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
    free = ~self.pinned
    r = _b - _A @ _T_lin
    S_sys = float(_E[free].sum() + (_M @ (_T_lin - _Tk))[free].sum())
    S_true = float(self.operators(T_new, T_old)[2].sum())
    rad_lin = float(_Frad[free].sum() + (_B @ (_T_lin - _Tk))[free].sum())
    rad_true = self.radiated_power(T_amb)
    DECOMP.append((iteration, last_relative_change, S_true - S_sys, dt * (rad_true - rad_lin), -dt * float(r[free].sum()),
                   float(np.abs(T_new - _T_lin).max()), getattr(self, "last_cg_iterations", -1), self.direct_fallbacks))
    self.last_damping = damping
    return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)

