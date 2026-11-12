"""The shared Material's properties tabulated for the level set's skins (stage 3 Task 4): h(T), cp_eff(T), k(T) and the
liquid fraction f_l(T) on a uniform grid of `step` (0.01 K) over T_min to T_max (150-3000 K), read back by linear
interpolation, as the skyscrapers' PropertyTable does (matched, never imported). A numerical evaluation device, not new
physics: every value is the Material's own at the grid's nodes, and a lookup is index arithmetic instead of the
Material's binary searches, which dominated a conduction sub-step (Task 3 measured about 2 s for 20 000 columns of 64
cells on AA7075).

The skins' books close on the table's own h: they store each cell's enthalpy as m h_table(T) and read its temperature by
inverting the same piecewise-linear h (`temperature_from_enthalpy`, exact to round-off, so T(h(T)) returns T), and the
conduction's Newton tangent is the slope of the tabulated h across the grid interval (`dh_dT`), so the Jacobian is that
of the residual it solves. The table's h differs from the Material's by under 1e-9 relative on AA7075 (the Material's h
is quadratic inside each c_p interval). cp_eff is interpolated between one-sided values at each interval's ends, so a
step of the Material's cp_eff that falls on a grid node (Scheil's nodes at whole kelvins) is kept as a step rather than
smeared over 0.01 K. Beyond the range each property continues the end interval's line."""
import numpy as np

TABLE_STEP = 0.01               # K, the grid (the skyscrapers' TABLE_DT)
TABLE_RANGE = (150.0, 3000.0)   # K, the tabulated range (the skyscrapers' TABLE_RANGE)
ONE_SIDED = 1.0e-6              # of a step: how far inside each interval its end values of cp_eff are taken
BUCKETS_MAX = 8                 # the inversion's uniform enthalpy buckets at most this many times the grid's nodes


class PropertyTable:
    """`material`'s enthalpy above its T_REF [J/kg], cp_eff [J/(kg K)], k [W/(m K)] and liquid fraction on the grid
    T_min + i step, i = 0 .. (T_max - T_min) / step (module docstring). Every method takes arrays of any shape."""

    def __init__(self, material, T_min=TABLE_RANGE[0], T_max=TABLE_RANGE[1], step=TABLE_STEP):
        if not (step > 0.0 and T_max > T_min + step):
            raise ValueError("a property table needs step > 0 and at least two intervals")
        self.T0, self.step = float(T_min), float(step)
        self.size = int(round((T_max - T_min) / step)) + 1
        T = self.T0 + self.step * np.arange(self.size)
        self._h = np.asarray(material.enthalpy(T), dtype=float)
        self._k = np.asarray(material.k(T), dtype=float)
        self._f = np.asarray(material.liquid_fraction(T), dtype=float)
        eps = ONE_SIDED * self.step
        self._c0 = np.asarray(material.cp_eff(T[:-1] + eps), dtype=float)    # each interval's cp_eff at its start
        self._c1 = np.asarray(material.cp_eff(T[1:] - eps), dtype=float)     # and at its end
        self._dh = np.diff(self._h)
        if not np.all(self._dh > 0.0):
            raise ValueError("the material's enthalpy must rise with temperature over the table's range")
        self._dk = np.diff(self._k)
        # the inversion's index: buckets of uniform enthalpy no wider than the narrowest interval, so that at most one
        # node lies inside a bucket and a lookup is a division and one correction (searchsorted was 0.27 s for 1.3e6
        # cells, as long as the Material's own inversion); a table whose narrowest interval would need more than
        # BUCKETS_MAX buckets per node falls back to the binary search
        width = self._dh.min()
        n_buckets = int(np.ceil((self._h[-1] - self._h[0]) / width)) + 1
        self._bucket = None
        if n_buckets <= BUCKETS_MAX * self.size:
            self._bucket_width = width
            edges = self._h[0] + width * np.arange(n_buckets)
            self._bucket = np.clip(np.searchsorted(self._h, edges, side="right") - 1, 0, self.size - 2)

    def _locate(self, T):
        """Each temperature's grid interval i and its place x in [0, 1) within it (beyond the range, the end interval
        and x outside [0, 1]; a nan temperature, as a padded cell's, reads nan)."""
        u = (np.asarray(T, dtype=float) - self.T0) / self.step
        i = np.clip(np.floor(np.where(np.isfinite(u), u, 0.0)), 0, self.size - 2).astype(np.int64)
        return i, u - i

    def enthalpy(self, T):
        i, x = self._locate(T)
        return self._h[i] + x * self._dh[i]

    def dh_dT(self, T):
        """The slope of the tabulated h across the interval T falls in [J/(kg K)]: the conduction's Newton tangent."""
        i, _ = self._locate(T)
        return self._dh[i] / self.step

    def cp_eff(self, T):
        i, x = self._locate(T)
        return self._c0[i] + x * (self._c1[i] - self._c0[i])

    def k(self, T):
        i, x = self._locate(T)
        return self._k[i] + x * self._dk[i]

    def dk_dT(self, T):
        """The slope of the tabulated k across the interval T falls in [W/(m K2)]."""
        i, _ = self._locate(T)
        return self._dk[i] / self.step

    def liquid_fraction(self, T):
        i, x = self._locate(T)
        return np.clip(self._f[i] + x * (self._f[i + 1] - self._f[i]), 0.0, 1.0)

    def conduction(self, T):
        """(h, dh/dT, k, dk/dT) at T with one index computation, for the conduction's Newton iteration."""
        i, x = self._locate(T)
        dh, dk = self._dh[i], self._dk[i]
        return self._h[i] + x * dh, dh / self.step, self._k[i] + x * dk, dk / self.step

    def temperature_from_enthalpy(self, h):
        """The temperature [K] whose tabulated enthalpy is h [J/kg]: the inverse of the piecewise-linear `enthalpy`."""
        h = np.asarray(h, dtype=float)
        if self._bucket is None:
            i = np.clip(np.searchsorted(self._h, h, side="right") - 1, 0, self.size - 2)
        else:
            j = np.clip(np.floor((h - self._h[0]) / self._bucket_width), 0, len(self._bucket) - 1).astype(np.int64)
            i = self._bucket[j]
            i = i + ((h >= self._h[np.minimum(i + 1, self.size - 1)]) & (i < self.size - 2))   # a node in the bucket
            i = i - ((h < self._h[i]) & (i > 0))                                              # round-off of j
        return self.T0 + self.step * (i + (h - self._h[i]) / self._dh[i])
