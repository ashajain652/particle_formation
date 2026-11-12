"""Test-only materials for the level-set tests: constant properties (`ConstMaterial`) and constant properties with a
latent heat released linearly over a mushy range (`LinearMeltMaterial`). Moved here from test_levelset_cutthermal.py in
stage 3 Task 2 and extended with what the skins read from the shared `Material` (`liquid` with rho, mu and sigma,
`liquid_fraction`, `temperature_from_enthalpy`, `T_solidus`, `T_liquidus`, `T_feed`, `T_rigid`). Test-only: no physics
of the level set lives here; the package reads the shared Material."""
from types import SimpleNamespace

import numpy as np


class ConstMaterial:
    """Constant density, conductivity and heat capacity; never melts. Its liquid properties (asked for by the skins'
    mixture density) are the solid's density with the shared AA7075 liquid's viscosity and surface tension."""

    def __init__(self, rho=2800.0, k=150.0, cp=900.0, T_ref=293.0, rho_l=None, mu=1.3e-3, sigma=0.8):
        self.rho, self._k, self._cp, self.T_ref = rho, k, cp, T_ref
        self.liquid = SimpleNamespace(rho=rho if rho_l is None else rho_l, mu=mu, sigma=sigma)

    T_solidus = T_liquidus = T_feed = T_rigid = np.inf

    def k(self, T):
        return np.full_like(np.asarray(T, float), self._k)

    def cp_eff(self, T):
        return np.full_like(np.asarray(T, float), self._cp)

    def enthalpy(self, T):
        return self._cp * (np.asarray(T, float) - self.T_ref)

    def liquid_fraction(self, T):
        return np.zeros_like(np.asarray(T, float))

    def temperature_from_enthalpy(self, h):
        return self.T_ref + np.asarray(h, float) / self._cp


class LinearMeltMaterial(ConstMaterial):
    """Constant properties plus a latent heat L released linearly between T_s and T_l (a mushy range)."""

    def __init__(self, L=390e3, T_s=880.0, T_l=908.0, rho_l=2400.0, **kw):
        super().__init__(rho_l=rho_l, **kw)
        self.L, self.T_s, self.T_l = L, T_s, T_l

    @property
    def T_solidus(self):
        return self.T_s

    @property
    def T_liquidus(self):
        return self.T_l

    @property
    def T_feed(self):
        return self.T_l

    @property
    def T_rigid(self):
        return 0.5 * (self.T_s + self.T_l)          # half liquid, as the shared Material's T_rigid

    def liquid_fraction(self, T):
        return np.clip((np.asarray(T, float) - self.T_s) / (self.T_l - self.T_s), 0.0, 1.0)

    def cp_eff(self, T):
        T = np.asarray(T, float)
        return self._cp + np.where((T > self.T_s) & (T < self.T_l), self.L / (self.T_l - self.T_s), 0.0)

    def enthalpy(self, T):
        return super().enthalpy(T) + self.L * self.liquid_fraction(T)

    def temperature_from_enthalpy(self, h):
        h = np.asarray(h, float)
        h_s, h_l = self._cp * (self.T_s - self.T_ref), self._cp * (self.T_l - self.T_ref) + self.L
        mush = self.T_s + (h - h_s) / (self._cp + self.L / (self.T_l - self.T_s))
        return np.where(h <= h_s, self.T_ref + h / self._cp,
                        np.where(h >= h_l, self.T_ref + (h - self.L) / self._cp, mush))


class VariableKMaterial(LinearMeltMaterial):
    """`LinearMeltMaterial` with a conductivity rising linearly with temperature, k = k0 (1 + slope (T - T_ref)), so
    that a heat sub-step's conductances move with its temperatures (stage 3 Task 4: the dk/dT terms of the Jacobian)."""

    def __init__(self, slope=1.0e-3, **kw):
        super().__init__(**kw)
        self.k_slope = slope

    def k(self, T):
        return self._k * (1.0 + self.k_slope * (np.asarray(T, float) - self.T_ref))
