"""Equilibrium-air stagnation state and the stagnation-point heating correlations (spec section 6.3).

Thermodynamics and composition: Cantera with `airNASA9.yaml` (NASA-9 polynomials, 200-20000 K, 11 species incl.
ions). That mechanism carries no transport data, so viscosity comes from Blottner's curve fits for N2, O2, NO, N, O
(Blottner, Johnson & Ellis 1971, as tabulated in Gnoffo, Gupta & Shinn, NASA TP-2867, 1989) combined with Wilke's
mixing rule over the equilibrium mole fractions; ions are neglected in the mixture viscosity. Freestream
composition fixed at N2:O2 = 0.79:0.21 by mole; T_inf is clamped to >= 200 K (the polynomial floor; the
enthalpy difference is negligible against V^2/2).
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

AIR = "N2:0.79, O2:0.21"
T_FLOOR = 200.0
CP_AIR = 1004.5              # J/(kg K), perfect-gas air: the stagnation-temperature estimate that decides whether to equilibrate
T_EQUILIBRATE = 1500.0       # K; below this the composition is frozen (no dissociation) and Cantera's equilibrium solver is skipped
PRANDTL = 0.71
LEWIS = 1.4
BLOTTNER = {                       # mu [Pa s] = 0.1 exp((A ln T + B) ln T + C), T in K
    "N2": (0.0268142, 0.3177838, -11.3155513),
    "O2": (0.0449290, -0.0826158, -9.2019475),
    "NO": (0.0436378, -0.0335511, -9.5767430),
    "N": (0.0115572, 0.6031679, -12.4327495),
    "O": (0.0203144, 0.4294404, -11.6031403),
}


@dataclass(frozen=True)
class GasState:
    p: float          # Pa
    T: float          # K
    rho: float        # kg/m3
    h: float          # J/kg (Cantera basis: elements at 298.15 K)
    mu: float         # Pa s
    h_D: float        # J/kg, chemical enthalpy of the free atoms (dissociation enthalpy)
    X: dict           # mole fractions > 1e-6


def blottner_viscosity(species, T):
    A, B, C = BLOTTNER[species]
    lnT = math.log(T)
    return 0.1 * math.exp((A * lnT + B) * lnT + C)


def wilke_viscosity(T, mole_fractions, molar_masses):
    """Mixture viscosity [Pa s] of the Blottner species present in `mole_fractions` (dict name -> X), Wilke's rule."""
    names = [s for s in BLOTTNER if mole_fractions.get(s, 0.0) > 0.0]
    x = np.array([mole_fractions[s] for s in names]); x = x / x.sum()
    mu = np.array([blottner_viscosity(s, T) for s in names])
    m = np.array([molar_masses[s] for s in names])
    total = 0.0
    for i in range(len(names)):
        phi = sum(x[j] * (1.0 + math.sqrt(mu[i] / mu[j]) * (m[j] / m[i]) ** 0.25) ** 2 / math.sqrt(8.0 * (1.0 + m[i] / m[j]))
                  for j in range(len(names)))
        total += x[i] * mu[i] / phi
    return float(total)


class EquilibriumAir:
    """Normal shock + isentropic compression to rest in equilibrium air, and the wall state at (T_w, p_s)."""

    def __init__(self, mechanism="airNASA9.yaml"):
        import cantera as ct                    # imported here: only `--heating physics` needs Cantera
        self.gas = ct.Solution(mechanism)
        self.molar_masses = dict(zip(self.gas.species_names, self.gas.molecular_weights))
        self._h_atoms = {s: self.gas.species(s).thermo.h(298.15) / self.molar_masses[s] for s in ("N", "O")}   # J/kg
        T = np.arange(T_FLOOR, 6000.0 + 1.0, 25.0)
        self._h_table = (T, np.array([self._enthalpy_TP(t, 101325.0) for t in T]))

    def _enthalpy_TP(self, T, p):
        self.gas.TPX = T, p, AIR
        return float(self.gas.enthalpy_mass)

    def air_enthalpy(self, T):
        """Frozen-composition air enthalpy h(T) [J/kg] (ideal gas: pressure-independent), table lookup, vectorised."""
        return np.interp(np.asarray(T, dtype=float), self._h_table[0], self._h_table[1])

    def _state(self):
        g = self.gas
        X = {s: float(x) for s, x in zip(g.species_names, g.X) if x > 1e-6}
        Y = dict(zip(g.species_names, g.Y))
        h_D = sum(Y[s] * self._h_atoms[s] for s in ("N", "O"))
        return GasState(float(g.P), float(g.T), float(g.density), float(g.enthalpy_mass),
                        wilke_viscosity(float(g.T), X, self.molar_masses), float(h_D), X)

    def freestream(self, rho, T):
        self.gas.TDX = max(T, T_FLOOR), rho, AIR
        return self._state()

    def stagnation(self, rho_inf, T_inf, V):
        """Equilibrium stagnation state behind the bow shock for freestream (rho, T, V)."""
        g = self.gas
        T_inf = max(T_inf, T_FLOOR)
        g.TDX = T_inf, rho_inf, AIR
        p_inf, h_inf = float(g.P), float(g.enthalpy_mass)
        equilibrate = T_inf + 0.5 * V * V / CP_AIR > T_EQUILIBRATE     # frozen air below ~1500 K (subsonic/low supersonic states)
        eps = 0.1                                                     # rho_inf / rho_2, Rankine-Hugoniot fixed point
        converged = False
        for _ in range(3000):
            u2 = eps * V
            p2 = p_inf + rho_inf * V * (V - u2)
            g.HPX = h_inf + 0.5 * (V * V - u2 * u2), p2, AIR
            if equilibrate:
                g.equilibrate("HP")
            eps_new = rho_inf / float(g.density)
            if abs(eps_new - eps) < 1e-10:
                converged = True
                break
            eps = 0.5 * (eps + eps_new)
        if not converged:
            raise RuntimeError(f"normal shock iteration did not converge in 3000 iterations (rho {rho_inf:.3e} kg/m3, T {T_inf:.1f} K, V {V:.0f} m/s, last change {abs(eps_new - eps):.2e})")
        s2, h_s = float(g.entropy_mass), h_inf + 0.5 * V * V             # isentropic compression to rest

        def residual(log_p):
            g.SPX = s2, math.exp(log_p), AIR
            if equilibrate:
                g.equilibrate("SP")
            return float(g.enthalpy_mass) - h_s

        p_s = math.exp(brentq(residual, math.log(p2), math.log(2.0 * p2), xtol=1e-12))
        g.SPX = s2, p_s, AIR
        if equilibrate:
            g.equilibrate("SP")
        return self._state()

    def wall(self, T_w, p):
        """Air at the wall temperature and stagnation pressure (equilibrium; molecular below ~2500 K)."""
        self.gas.TPX = T_w, p, AIR
        if T_w > T_EQUILIBRATE:
            self.gas.equilibrate("TP")
        return self._state()


def dkr(rho, V, radius):
    """Detra-Kemp-Riddell stagnation-point heating [W/m2], SESAM's continuum correlation (facts note s.4)."""
    return 1.1035e8 * radius ** -0.5 * (rho / 1.225) ** 0.5 * (V / 7925.0) ** 3.15


def sutton_graves(rho, V, radius):
    """Sutton-Graves (1971) stagnation-point heating for air [W/m2], cold wall."""
    return 1.7415e-4 * math.sqrt(rho / radius) * V ** 3


def fay_riddell(stag, wall, p_inf, radius, catalycity=1.0, Pr=PRANDTL, Le=LEWIS):
    """Fay-Riddell (1958) equilibrium-boundary-layer stagnation heating [W/m2].

    `catalycity` scales the Lewis-number (atom diffusion + wall recombination) term; 1 is the fully catalytic
    equilibrium form, 0 leaves the Le = 1 value. The frozen-boundary-layer non-catalytic reduction
    (1 - h_D/h_s) is not represented (Step 3 option)."""
    due_dx = math.sqrt(2.0 * (stag.p - p_inf) / stag.rho) / radius
    return 0.763 * Pr ** -0.6 * (wall.rho * wall.mu) ** 0.1 * (stag.rho * stag.mu) ** 0.4 * math.sqrt(due_dx) \
        * (stag.h - wall.h) * (1.0 + (Le ** 0.52 - 1.0) * catalycity * stag.h_D / stag.h)
