"""Compressible boundary-layer similarity solution (Lees local similarity, total-enthalpy energy equation with
dissipation, real frozen-air properties) and its march along a body meridian, ported from
analysis/gas_loads/reference/similarity_heat.py.

    (C f'')' + f f'' + beta (T/T_e - f'^2) = 0
    (C g'/Pr)' + f g' + (u_e^2/H_e) [ (1 - 1/Pr) C f' f'' ]' = 0,     g = H/H_e
    h = H_e g - u_e^2 f'^2 / 2,  T = T(h),  C = rho mu / (rho_w mu_w) = (T_w / T) mu(T)/mu(T_w),
    f(0) = f'(0) = 0, f'(inf) = 1;  g(0) = h_w/H_e, g(inf) = 1;
    tau_w = rho_w mu_w u_e^2 r f''(0) / sqrt(2 xi),   q_w = rho_w mu_w u_e r H_e g'(0) / (Pr sqrt(2 xi)).
For a calorically perfect or frozen gas rho_e/rho = T/T_e exactly. Past separation (f''(0) <= 0.005, or no converged
solution) nothing is computed and the loads are zero.
"""
import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_bvp
from scipy.interpolate import CubicSpline

from . import gas

PRANDTL = 0.71
ETA = np.linspace(0.0, 12.0, 500)
AIR_X = {"N2": 0.79, "O2": 0.21}
F_SEPARATION = 0.005
HOMOTOPY_STEPS = 12


class SimilarityError(RuntimeError):
    """The boundary-value solver did not converge."""


class FrozenAirProperties:
    """Smooth (cubic-spline) frozen-air h(T), T(h), mu(T) and their derivatives on 150 to 20000 K every 50 K: the
    collocation solver needs a continuous right-hand side, and piecewise-linear tables make it refine without end."""

    def __init__(self, air):
        g = air.gas
        self.T = np.arange(150.0, 20000.0 + 1.0, 50.0)
        h = []
        for t in self.T:
            g.TPX = t, 101325.0, "N2:0.79, O2:0.21"
            h.append(g.enthalpy_mass)
        self.h = np.array(h)
        self.mu = np.array([gas.wilke_viscosity(float(t), AIR_X, air.molar_masses) for t in self.T])
        self._h = CubicSpline(self.T, self.h)
        self._T = CubicSpline(self.h, self.T)
        self._mu = CubicSpline(self.T, self.mu)
        self._cp, self._dmu = self._h.derivative(), self._mu.derivative()

    def T_of_h(self, h):
        return self._T(np.clip(h, self.h[0], self.h[-1]))

    def h_of_T(self, T):
        return float(self._h(T))

    def mu_of_T(self, T):
        return self._mu(T)

    def dmu_of_T(self, T):
        return self._dmu(T)

    def cp_of_T(self, T):
        return self._cp(T)


class ConstantProperties:
    """Calorically perfect gas with viscosity proportional to temperature, so C = 1 and T/T_e = h/h_e (the
    Falkner-Skan problem)."""

    def __init__(self, cp=1005.0, mu_per_K=1.0e-7):
        self.cp, self.k = cp, mu_per_K

    def T_of_h(self, h):
        return np.asarray(h, dtype=float) / self.cp

    def h_of_T(self, T):
        return float(self.cp * T)

    def mu_of_T(self, T):
        return self.k * np.asarray(T, dtype=float)

    def dmu_of_T(self, T):
        return self.k * np.ones_like(np.asarray(T, dtype=float))

    def cp_of_T(self, T):
        return self.cp * np.ones_like(np.asarray(T, dtype=float))


@dataclass
class SimilarityPoint:
    fpp0: float     # f''(0), the dimensionless wall shear
    gp0: float      # g'(0), the dimensionless wall enthalpy gradient
    g_w: float      # g at the wall, h_w / H_e


class SimilaritySolver:
    """One similarity station at a time, warm-started from the previous solution; the first solve (and the first
    after reset()) eases in from a hot wall by a 12-step geometric homotopy in the wall temperature."""

    def __init__(self, props, pr=PRANDTL):
        self.props, self.pr, self._sol = props, pr, None

    def reset(self):
        self._sol = None

    def solve(self, beta, H_e, u_e, T_w):
        if self._sol is None:
            T_e = float(self.props.T_of_h(H_e - 0.5 * u_e * u_e))
            for t in np.geomspace(0.8 * T_e, T_w, HOMOTOPY_STEPS)[:-1]:
                self._solve(beta, H_e, u_e, float(t))
        return self._solve(beta, H_e, u_e, T_w)

    def _solve(self, beta, He, ue, Tw_):
        tab, PR = self.props, self.pr
        hw = tab.h_of_T(Tw_)
        gw = hw / He
        Te = float(tab.T_of_h(He - 0.5 * ue * ue))
        muw = float(tab.mu_of_T(Tw_))
        D = (1.0 - 1.0 / PR) * ue * ue / He

        def props(y):
            f, fp, fpp, g, gp = y
            h = He * g - 0.5 * ue * ue * fp * fp
            T = np.clip(tab.T_of_h(h), 150.0, 20000.0)
            C = (Tw_ / T) * tab.mu_of_T(T) / muw
            dCdT = (-(Tw_ / T**2) * tab.mu_of_T(T) + (Tw_ / T) * tab.dmu_of_T(T)) / muw
            dh = He * gp - ue * ue * fp * fpp
            Cp = dCdT * dh / tab.cp_of_T(T)
            return T, C, Cp

        def rhs(eta, y):
            f, fp, fpp, g, gp = y
            T, C, Cp = props(y)
            fppp = -(Cp * fpp + f * fpp + beta * (T / Te - fp**2)) / C
            gpp = -(PR * (f * gp + D * (Cp * fp * fpp + C * fpp**2 + C * fp * fppp)) + Cp * gp) / C
            return np.vstack([fp, fpp, fppp, gp, gpp])

        bc = lambda a, b: np.array([a[0], a[1], b[1] - 1.0, a[3] - gw, b[3] - 1.0])
        if self._sol is None:
            e = np.exp(-ETA)
            guess = (ETA, np.vstack([ETA - 1 + e, 1 - e, e, 1 - (1 - gw) * e, (1 - gw) * e]))
        else:
            guess = (self._sol.x, self._sol.y)
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            s = solve_bvp(rhs, bc, *guess, tol=1e-4, max_nodes=200000)
        if not s.success:
            raise SimilarityError(s.message)
        self._sol = s
        return SimilarityPoint(float(s.y[2, 0]), float(s.y[4, 0]), gw)


@dataclass
class MeridianSolution:
    theta_deg: np.ndarray
    beta: np.ndarray
    fpp0: np.ndarray
    gp0: np.ndarray
    tau_c: np.ndarray       # Pa, continuum wall shear
    q_c: np.ndarray         # W/m2, continuum wall heat flux
    separation_deg: float   # nan if attached to the last bin


def march(theta_deg, radius, p_w, u_e, H_e, T_w, rho_w_mu_w, solver):
    """Continuum wall shear and heat flux along the meridian bins theta_deg [deg from the stagnation point]; xi is
    integrated on a 20-times refined grid, beta = 2 xi (du_e/ds) / (rho_w mu_w u_e^2 r_body^2), the stagnation point
    uses beta = 1/2 with u_e r / sqrt(2 xi) -> sqrt(2 beta0 / (rho_w mu_w)). The last bin is not marched (as in the
    reference); the first bin with f''(0) <= 0.005 or a failed solve is the separation angle."""
    th = np.asarray(theta_deg, dtype=float)
    n = len(th)
    s = radius * np.radians(th)
    r = radius * np.sin(np.radians(th))
    rmw = np.asarray(rho_w_mu_w, dtype=float)
    u_e = np.asarray(u_e, dtype=float)
    fine = np.linspace(0, s[-1], 20 * (n - 1) + 1)
    xi = np.interp(s, fine, cumulative_trapezoid(np.interp(fine, s, rmw) * np.interp(fine, s, u_e)
                                                 * (radius * np.sin(fine / radius))**2, fine, initial=0.0))
    due = np.gradient(u_e, s)
    beta, fpp0, gp0 = np.zeros(n), np.zeros(n), np.zeros(n)
    tau, q = np.zeros(n), np.zeros(n)
    sep = float("nan")
    beta0 = due[1]
    pt = solver.solve(0.5, H_e, 0.0, T_w)
    beta[0], fpp0[0], gp0[0] = 0.5, pt.fpp0, pt.gp0
    q[0] = math.sqrt(2 * rmw[0] * beta0) * H_e * pt.gp0 / solver.pr
    for i in range(1, n - 1):
        b = 2 * xi[i] * due[i] / (rmw[i] * u_e[i]**2 * r[i]**2)
        try:
            pt = solver.solve(b, H_e, u_e[i], T_w)
            fpp, gp = pt.fpp0, pt.gp0
        except SimilarityError:
            fpp, gp = 0.0, 0.0
        if fpp <= F_SEPARATION:
            sep = float(th[i])
            break
        beta[i], fpp0[i], gp0[i] = b, fpp, gp
        fac = rmw[i] * u_e[i] * r[i] / math.sqrt(2 * xi[i])
        tau[i] = fac * u_e[i] * fpp
        q[i] = fac * H_e * gp / solver.pr
    return MeridianSolution(th, beta, fpp0, gp0, tau, q, sep)
