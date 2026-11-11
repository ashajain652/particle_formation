"""Constant-property and simple melting materials for the analytic tests."""
import numpy as np


class ConstMaterial:
    def __init__(self, rho=2800.0, k=150.0, cp=900.0, T_ref=293.0):
        self.rho, self._k, self._cp, self.T_ref = rho, k, cp, T_ref

    def k(self, T):
        return np.full_like(np.asarray(T, float), self._k)

    def cp_eff(self, T):
        return np.full_like(np.asarray(T, float), self._cp)

    def enthalpy(self, T):
        return self._cp * (np.asarray(T, float) - self.T_ref)


class LinearMeltMaterial(ConstMaterial):
    """Constant properties plus a latent heat L released linearly between T_s and T_l (a mushy range)."""

    def __init__(self, L=390e3, T_s=880.0, T_l=908.0, **kw):
        super().__init__(**kw)
        self.L, self.T_s, self.T_l = L, T_s, T_l

    def liquid_fraction(self, T):
        return np.clip((np.asarray(T, float) - self.T_s) / (self.T_l - self.T_s), 0.0, 1.0)

    def cp_eff(self, T):
        T = np.asarray(T, float)
        return self._cp + np.where((T > self.T_s) & (T < self.T_l), self.L / (self.T_l - self.T_s), 0.0)

    def enthalpy(self, T):
        return super().enthalpy(T) + self.L * self.liquid_fraction(T)
