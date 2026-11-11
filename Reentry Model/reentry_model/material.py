"""Material for the thermal and melting model: density, k(T), c_p(T), emissivity, the latent heat with its melting
range, the liquid-phase properties of the film, and the specific enthalpy h(T) with its inverse.

Tables come from a DRAMA material JSON (the wrapper's `user_materials` format); np.interp holds the end values
outside the tabulated range. Melting (Step 3, spec section 6): `meltingHeat` L_f is released between the solidus and
the liquidus (`solidusTemperature`/`liquidusTemperature`, default both = `meltingTemperature`); a single-temperature
material gets a numerical ramp of +-MELT_RAMP around it (the enthalpy jump is exact, only its slope is smoothed).
`liquid_fraction(T)` is the enthalpy method's f_l; `feed_fraction(T)` is the fraction of an element's material that
the film receives -- a +-MELT_RAMP ramp at the liquidus for every material, so that melt and runoff start at the
liquidus and the mushy range counts as solid (decided 2026-09-20, spec sections 8 and 17.2). A `meltingTemperature`
above NO_MELT_ABOVE (the wrapper's 1e5 K device) means no melting at all. Enthalpy is exact: the integral of the
piecewise-linear c_p (quadratic inside each table interval) plus L_f f_l(T), and the inverse solves the same
quadratic, so the FEM's secant heat capacity conserves energy through melting to round-off.

A Scheil material (`solidification: "scheil"`, user decision of 2026-09-27) releases the same L_f along Scheil's
non-equilibrium curve f_l = ((T_pure - T)/(T_pure - T_liquidus))^(1/(k - 1)) instead of linearly: half of it in the
13 K below the liquidus for AA7075 (k = 0.4, T_pure 933 K), with the eutectic remainder f_l(T_solidus) (3.6 %)
melting or freezing over +-MELT_RAMP at the solidus. The curve is tabulated at SCHEIL_DT and its nodes join the
enthalpy table, so the latent slope stays constant inside every interval and the exact inversion is unchanged. The
linear range material (`AA7075_range`) is untouched by this: its code path is the one above, bit for bit."""
import json
import os
from dataclasses import dataclass, field

import numpy as np

from . import DATA_DIR

DEFAULT_MATERIAL = os.path.join(DATA_DIR, "materials", "AA7075_nomelt.json")
MATERIAL_NAMES = {"AA7075_nomelt": DEFAULT_MATERIAL, "AA7075": os.path.join(DATA_DIR, "materials", "AA7075.json"),
                  "AA7075_range": os.path.join(DATA_DIR, "materials", "AA7075_range.json"),
                  "AA7075_scheil": os.path.join(DATA_DIR, "materials", "AA7075_scheil.json"),
                  "AA7075-empiricaldata": os.path.join(DATA_DIR, "materials", "AA7075-empiricaldata.json")}
T_REF = 293.0                      # K, zero of the enthalpy scale (first row of DRAMA's tables)
MELT_RAMP = 2.0                    # K, half-width of the numerical melting/feed ramps
NO_MELT_ABOVE = 5000.0             # K: a melting temperature above this means "never melts"
SCHEIL_DT = 1.0                    # K, spacing of the tabulated Scheil liquid fraction (interpolation error < 1e-3)
RIGID_LIQUID_FRACTION = 0.5        # above this liquid fraction mush is a slurry, not a rigid substrate (amendment of 2026-10-06)


@dataclass
class LiquidProperties:
    rho: float                     # kg/m3
    mu: float                      # Pa s
    sigma: float                   # N/m


@dataclass
class Material:
    name: str
    rho: float                      # kg/m3
    emissivity: float
    T_cp: np.ndarray                # K, nodes of the c_p table
    cp_table: np.ndarray            # J/(kg K)
    T_k: np.ndarray                 # K, nodes of the k table
    k_table: np.ndarray             # W/(m K)
    latent_heat: float = 0.0        # J/kg
    T_solidus: float = np.inf
    T_liquidus: float = np.inf
    liquid: LiquidProperties = None
    partition_coefficient: float = None     # Scheil's k; None for the linear range (solidification "linear")
    T_pure: float = None                    # K, melting point of the pure solvent in Scheil's law
    _T_h: np.ndarray = field(init=False, repr=False)
    _h_nodes: np.ndarray = field(init=False, repr=False)
    _latent_slope: np.ndarray = field(init=False, repr=False)    # L_f df_l/dT inside each enthalpy-table interval
    _fl_T: np.ndarray = field(init=False, repr=False, default=None)  # Scheil: nodes of the tabulated liquid fraction
    _fl: np.ndarray = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.T_cp, self.cp_table = np.asarray(self.T_cp, dtype=float), np.asarray(self.cp_table, dtype=float)
        self.T_k, self.k_table = np.asarray(self.T_k, dtype=float), np.asarray(self.k_table, dtype=float)
        if self.latent_heat < 0.0 or self.T_liquidus < self.T_solidus:
            raise ValueError("latent heat must be >= 0 and the liquidus >= the solidus")
        if self.melts and self.T_solidus == self.T_liquidus:            # single-temperature material: numerical ramp
            self.T_solidus, self.T_liquidus = self.T_liquidus - MELT_RAMP, self.T_liquidus + MELT_RAMP
        if self.scheil:
            k, T_F = self.partition_coefficient, self.T_pure
            if not (0.0 < k < 1.0 and T_F > self.T_liquidus and self.T_liquidus - self.T_solidus > 2.0 * MELT_RAMP):
                raise ValueError("Scheil needs 0 < k < 1, T_pure above the liquidus and a range wider than 2 MELT_RAMP")
            # Scheil's law from the top of the eutectic ramp to the liquidus, tabulated; the eutectic remainder
            # f_l(T_solidus + MELT_RAMP) is released linearly over +-MELT_RAMP at the solidus, whose foot is then the
            # material's T_solidus (as a single-temperature material's ramp foot is)
            T_top = self.T_solidus + MELT_RAMP
            T_nodes = np.union1d(np.arange(np.ceil(T_top), self.T_liquidus, SCHEIL_DT), [T_top, self.T_liquidus])
            f_nodes = ((T_F - T_nodes) / (T_F - self.T_liquidus)) ** (1.0 / (k - 1.0))
            self.T_solidus = self.T_solidus - MELT_RAMP
            self._fl_T, self._fl = np.concatenate([[self.T_solidus], T_nodes]), np.concatenate([[0.0], f_nodes])
        # h(T) on the c_p nodes (plus T_REF and the melting range): the trapezoid rule is exact for the piecewise-linear c_p
        extra = [T_REF] + ([self.T_solidus, self.T_liquidus] if self.melts else []) + (list(self._fl_T) if self.scheil else [])
        T = np.union1d(self.T_cp, extra)
        cp = np.interp(T, self.T_cp, self.cp_table)
        h = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(T))])
        h = h - np.interp(T_REF, T, h)
        slope = np.zeros(len(T) - 1)
        if self.melts:
            inside = (T[:-1] >= self.T_solidus - 1e-9) & (T[1:] <= self.T_liquidus + 1e-9)
            if self.scheil:                     # piecewise linear on the table's own nodes: one slope per interval
                slope[inside] = (self.latent_heat * np.diff(self._liquid_fraction_raw(T)) / np.diff(T))[inside]
            else:
                slope[inside] = self.latent_heat / (self.T_liquidus - self.T_solidus)
            h = h + self.latent_heat * self._liquid_fraction_raw(T)
        self._T_h, self._h_nodes, self._latent_slope = T, h, slope

    @property
    def melts(self):
        return self.latent_heat > 0.0 and np.isfinite(self.T_liquidus)

    @property
    def scheil(self):
        """The latent heat follows Scheil's law (a `partition_coefficient` is given), not the linear range."""
        return self.partition_coefficient is not None

    @property
    def T_feed(self):
        """Top of the feed ramp: material at or above it is fully liquid for the film."""
        return self.T_liquidus + (0.0 if self.T_liquidus - self.T_solidus <= 2.0 * MELT_RAMP + 1e-9 else MELT_RAMP)

    @property
    def h_liquid(self):
        """Specific enthalpy of the film (liquid at the liquidus) [J/kg]."""
        return float(self.enthalpy(self.T_feed))

    @property
    def T_rigid(self):
        """Temperature at which the material is RIGID_LIQUID_FRACTION liquid [K]; inf for a non-melting material. Below it
        the mush is coherent and carries load, a rigid substrate for the melt film; above it the mush is a slurry that
        flows (Chen et al. 2016's semi-solid law ends at 50 % liquid and Li et al. 2014's slurry data begin there;
        amendment of 2026-10-06). Found by bisection on the material's own liquid fraction, which is monotonic."""
        if not self.melts:
            return np.inf
        lo, hi = self.T_solidus, self.T_liquidus
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self._liquid_fraction_raw(mid) < RIGID_LIQUID_FRACTION else (lo, mid)
        return 0.5 * (lo + hi)

    @classmethod
    def from_drama_json(cls, path=None):
        """A DRAMA material file; `path` may also be a name in MATERIAL_NAMES (AA7075_nomelt, AA7075, AA7075_range,
        AA7075_scheil, AA7075-empiricaldata). `solidification: "scheil"` with `partitionCoefficient` and
        `pureMeltingTemperature` selects Scheil's law; anything else keeps the linear range."""
        path = MATERIAL_NAMES.get(path, path) or DEFAULT_MATERIAL
        with open(path) as fh:
            d = json.load(fh)
        cp = np.array(d["specificHeatCapacity"], dtype=float)
        k = np.array(d["heatConductivity"], dtype=float)
        T_melt = float(d.get("meltingTemperature", np.inf))
        latent = float(d.get("meltingHeat", 0.0)) if T_melt < NO_MELT_ABOVE else 0.0
        T_s, T_l = float(d.get("solidusTemperature", T_melt)), float(d.get("liquidusTemperature", T_melt))
        liquid = LiquidProperties(float(d["liquid"]["density"]), float(d["liquid"]["viscosity"]), float(d["liquid"]["surfaceTension"])) if "liquid" in d else None
        scheil = latent > 0.0 and d.get("solidification") == "scheil"
        return cls(d["name"], float(d["density"]), float(d["emissivity"][0][1]), cp[:, 0], cp[:, 1], k[:, 0], k[:, 1],
                   latent, T_s if latent else np.inf, T_l if latent else np.inf, liquid,
                   float(d["partitionCoefficient"]) if scheil else None, float(d["pureMeltingTemperature"]) if scheil else None)

    def k(self, T):
        return np.interp(T, self.T_k, self.k_table)

    def cp(self, T):
        return np.interp(T, self.T_cp, self.cp_table)

    def _liquid_fraction_raw(self, T):
        if self._fl_T is not None:              # Scheil: the tabulated curve (np.interp holds 0 below, 1 above)
            return np.interp(np.asarray(T, dtype=float), self._fl_T, self._fl)
        return np.clip((np.asarray(T, dtype=float) - self.T_solidus) / (self.T_liquidus - self.T_solidus), 0.0, 1.0)

    def liquid_fraction(self, T):
        """Melt fraction f_l(T): 0 below the solidus, 1 above the liquidus, linear between for a range material and
        Scheil's law (tabulated at SCHEIL_DT, eutectic ramp at the solidus) for a Scheil one (zero for a non-melting
        material)."""
        T = np.asarray(T, dtype=float)
        return self._liquid_fraction_raw(T) if self.melts else np.zeros_like(T)

    def feed_fraction(self, T):
        """Fraction of an element's material the film receives: a +-MELT_RAMP ramp ending at T_feed (the liquidus)."""
        T = np.asarray(T, dtype=float)
        if not self.melts:
            return np.zeros_like(T)
        return np.clip((T - (self.T_feed - 2.0 * MELT_RAMP)) / (2.0 * MELT_RAMP), 0.0, 1.0)

    def cp_eff(self, T):
        """Effective heat capacity c_p + L_f df_l/dT."""
        T = np.asarray(T, dtype=float)
        c = self.cp(T)
        if self.melts and self.scheil:          # the slope of the enthalpy-table interval T falls in
            i = np.clip(np.searchsorted(self._T_h, T, side="right") - 1, 0, len(self._T_h) - 2)
            c = c + np.where((T > self.T_solidus) & (T < self.T_liquidus), self._latent_slope[i], 0.0)
        elif self.melts:
            c = c + np.where((T > self.T_solidus) & (T < self.T_liquidus), self.latent_heat / (self.T_liquidus - self.T_solidus), 0.0)
        return c

    def enthalpy(self, T):
        """Specific enthalpy above T_REF [J/kg]: the exact integral of the piecewise-linear c_p (quadratic inside each
        table interval, linear beyond the table) plus L_f f_l(T)."""
        T = np.asarray(T, dtype=float)
        i = np.clip(np.searchsorted(self._T_h, T, side="right") - 1, 0, len(self._T_h) - 2)
        T0, T1 = self._T_h[i], self._T_h[i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        x = np.clip(T, T0, T1) - T0
        h = self._h_nodes[i] + (cp0 + self._latent_slope[i]) * x + 0.5 * (cp1 - cp0) / (T1 - T0) * x * x
        h = np.where(T < self._T_h[0], self._h_nodes[0] + self.cp_table[0] * (T - self._T_h[0]), h)
        h = np.where(T > self._T_h[-1], self._h_nodes[-1] + self.cp_table[-1] * (T - self._T_h[-1]), h)
        return h

    def enthalpy_liquid(self, T):
        """Specific enthalpy of fully liquid material [J/kg]: h(T) with the whole latent heat added, whatever the
        equilibrium liquid fraction at T would be. The melt film is liquid by construction -- that is what makes it a
        film -- so this, not the mixture enthalpy h(T), is what a kilogram of film holds, and feeding the film costs
        the latent heat the mass has not yet paid. Book the film at h(T) instead and melting a body held on the ramp
        is free: it turned a 100 mm sphere entirely to film on a quarter of its latent heat (measured 2026-09-22)."""
        T = np.asarray(T, dtype=float)
        return self.enthalpy(T) + (self.latent_heat * (1.0 - self.liquid_fraction(T)) if self.melts else 0.0)

    def enthalpy_mixed(self, T, w):
        """Specific enthalpy of a node holding a fraction `w` of melt film and 1 - w of ordinary material: the film
        is liquid, so it carries its latent heat at every temperature, and the material carries L_f f_l(T)."""
        w = np.asarray(w, dtype=float)
        return self.enthalpy(T) + (w * self.latent_heat * (1.0 - self.liquid_fraction(T)) if self.melts else 0.0)

    def cp_mixed(self, T, w):
        """d/dT of enthalpy_mixed: the film has no latent plateau, the material does."""
        w = np.asarray(w, dtype=float)
        return (1.0 - w) * self.cp_eff(T) + w * self.cp(T) if self.melts else self.cp(T)

    def temperature_from_enthalpy_mixed(self, h, w):
        """Inverse of enthalpy_mixed in T (monotonic in T for every w). A node the film owns has a shorter latent
        plateau -- only its material part has one -- so inverting the material's h(T) there throws the node clean
        across the ramp and the Newton iteration limit-cycles between 777 K and 864 K (measured 2026-09-22)."""
        h, w = np.asarray(h, dtype=float), np.broadcast_to(np.asarray(w, dtype=float), np.shape(h))
        T = self.temperature_from_enthalpy(h)                          # exact, and what a node with no film gets
        hot = w > 0.0
        if not self.melts or not hot.any():
            return T
        # the same quadratic inversion, against this node's own enthalpy table: the latent plateau is shortened to
        # (1 - w) of its height and the whole table is lifted by w L_f below the solidus
        hw, ww = h[hot], w[hot]
        tab = self._h_nodes[None, :] + ww[:, None] * self.latent_heat * (1.0 - self._liquid_fraction_raw(self._T_h))[None, :]
        i = np.clip((tab <= hw[:, None]).sum(axis=1) - 1, 0, len(self._T_h) - 2)
        rows, T0, T1 = np.arange(len(hw)), self._T_h[i], self._T_h[i + 1]
        h0, h1 = tab[rows, i], tab[rows, i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        b = cp0 + (1.0 - ww) * self._latent_slope[i]
        a = (cp1 - cp0) / (T1 - T0)
        dh = np.clip(hw, h0, h1) - h0
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(np.abs(a) > 1e-12, (np.sqrt(b * b + 2.0 * a * dh) - b) / a, dh / b)
        Tw = np.where(hw < tab[:, 0], self._T_h[0] + (hw - tab[:, 0]) / self.cp_table[0],
                      np.where(hw > tab[:, -1], self._T_h[-1] + (hw - tab[:, -1]) / self.cp_table[-1], T0 + x))
        T = T.copy()
        T[hot] = Tw
        return T

    def temperature_from_enthalpy(self, h):
        """Inverse of enthalpy() (monotonic): the energy-equivalent temperature of a body holding h per kg."""
        h = np.asarray(h, dtype=float)
        i = np.clip(np.searchsorted(self._h_nodes, h, side="right") - 1, 0, len(self._T_h) - 2)
        T0, T1 = self._T_h[i], self._T_h[i + 1]
        cp0, cp1 = np.interp(T0, self.T_cp, self.cp_table), np.interp(T1, self.T_cp, self.cp_table)
        b = cp0 + self._latent_slope[i]
        a = (cp1 - cp0) / (T1 - T0)
        dh = np.clip(h, self._h_nodes[i], self._h_nodes[i + 1]) - self._h_nodes[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(np.abs(a) > 1e-12, (np.sqrt(b * b + 2.0 * a * dh) - b) / a, dh / b)
        T = T0 + x
        T = np.where(h < self._h_nodes[0], self._T_h[0] + (h - self._h_nodes[0]) / self.cp_table[0], T)
        T = np.where(h > self._h_nodes[-1], self._T_h[-1] + (h - self._h_nodes[-1]) / self.cp_table[-1], T)
        return T
