"""Particle microphysics primitives of the wake-column model: size bins, particle size, settling, collision kernels.

Everything here is a small pure function or table so that the box model (box.py) can assemble them.

Bins (decision 8). A particle is characterised by its number of aluminium atoms n. The bins are pivots
n_k = 2 q^k (q = bin_ratio), k = 0..K-1, from a two-atom alumina dimer (0.476 nm) up to the first pivot whose alumina
diameter reaches d_max. Number and aluminium atoms are both conserved by fixed-pivot splitting (Kumar and Ramkrishna
1996): a freshly made entity of n atoms with n_k <= n < n_{k+1} puts the fraction zeta = (n_{k+1} - n)/(n_{k+1} - n_k)
of its number in bin k and 1 - zeta in bin k+1, so zeta n_k + (1 - zeta) n_{k+1} = n.
  - n exactly at pivot k gives zeta = 1 for that pivot.
  - n at or above the last pivot cannot be split; the entity goes to the last bin with its number scaled by
    n / n_last, which conserves atoms (number is then not conserved, as it cannot be: a bigger-than-largest particle
    stands for several largest particles). `split` returns (K-2, 0.0) in that case, `number_scale` returns n / n_last
    and `place` returns the bins and number fractions to add, with the scaling already applied.
  - n below the first pivot is not a particle and raises ValueError.

Size and density (decision 10). Volume per aluminium atom is (1 - f) v_Al2O3 + f v_Al(OH)3 with f = min(H/Al / 3, 1);
the particle is a sphere of that total volume.

Free-molecule kernel and capture rate (decision 7). Hard-sphere free-molecule collisions have rate coefficient
pi (r1 + r2)^2 sqrt(8 k_B T / (pi mu)). For two molecules A and B the formation (capture) rate coefficient is the
larger of the hard-sphere value (sigma_AB = (sigma_A + sigma_B)/2) and the Keesom-averaged dipole-dipole capture rate
k_DD = pi Gamma(2/3) vbar (2 C6 / k_B T)^(1/3), C6 = (2/3)(mu_A mu_B / 4 pi eps0)^2 / (k_B T); the dipole term is used
only when both dipole moments are known and nonzero. The spin factor g_P / (g_A g_B) is the share of encounters on the
attractive surface (decision 7).
"""
import csv
import math
import os
from dataclasses import dataclass
from itertools import combinations_with_replacement

import numpy as np
from scipy.constants import epsilon_0, k as K_B
from scipy.special import gamma

from reentry_model.wake_column import DATA
from reentry_model.wake_column.ambient import _interp_log
from reentry_model.wake_column.species import AL_SPECIES, AMU, ANGSTROM, DEBYE
from reentry_model.wake_column.transport import AIR_MASS_AMU

V_AL2O3 = 2.8218e-29            # m3 per Al atom in alumina Al2O3 (3.0 g/cm3), decision 10
V_ALOH3 = 5.3524e-29            # m3 per Al atom in aluminium hydroxide Al(OH)3 (2.42 g/cm3), decision 10
SMOKE_CSV = os.path.join(DATA, "megner2006_smoke.csv")
SMOKE_DLNR = math.log(200.0) / 34.0     # width of every bin of the smoke table in ln r (its header)


# ---- size bins ------------------------------------------------------------------------------------------------

def volume_per_al(h_per_al):
    """Particle volume per aluminium atom (m3) for a hydrogen-to-aluminium ratio h_per_al (decision 10)."""
    f = min(max(float(h_per_al) / 3.0, 0.0), 1.0)
    return (1.0 - f) * V_AL2O3 + f * V_ALOH3


def diameter(n_al, h_per_al):
    """Diameter (m) of the sphere holding n_al aluminium atoms at hydrogen-to-aluminium ratio h_per_al."""
    return (6.0 * n_al * volume_per_al(h_per_al) / math.pi) ** (1.0 / 3.0)


class Bins:
    """Aluminium-atom pivots n_k = 2 q^k and fixed-pivot splitting (module docstring, decision 8)."""

    def __init__(self, bin_ratio=2.0, d_max_nm=400.0):
        if not bin_ratio > 1.0:
            raise ValueError(f"bin_ratio must exceed 1, got {bin_ratio}")
        self.bin_ratio = float(bin_ratio)
        self.d_max_nm = float(d_max_nm)
        pivots = []
        k = 0
        while True:
            n = 2.0 * self.bin_ratio ** k           # exact for q = 2 (powers of two)
            pivots.append(n)
            if diameter(n, 0.0) >= self.d_max_nm * 1e-9:
                break
            k += 1
        self.n_atoms = np.array(pivots)
        self.K = len(pivots)

    def split(self, n_atoms):
        """(k, zeta): lower bin with n_k <= n < n_{k+1} and the number fraction zeta put there.

        At or above the last pivot returns (K-2, 0.0): everything goes to the last bin (see `number_scale`)."""
        n = float(n_atoms)
        if not n >= self.n_atoms[0]:
            raise ValueError(f"{n} aluminium atoms is below the smallest bin ({self.n_atoms[0]})")
        if n >= self.n_atoms[-1]:
            return self.K - 2, 0.0
        k = int(np.searchsorted(self.n_atoms, n, side="right")) - 1
        lo, hi = self.n_atoms[k], self.n_atoms[k + 1]
        return k, float((hi - n) / (hi - lo))

    def number_scale(self, n_atoms):
        """Factor on the particle number so atoms are conserved: n / n_last above the last pivot, else 1."""
        return max(1.0, float(n_atoms) / self.n_atoms[-1])

    def place(self, n_atoms):
        """[(bin, number fraction), ...] to add per entity of n_atoms atoms; conserves atoms (and number in range)."""
        k, z = self.split(n_atoms)
        if float(n_atoms) >= self.n_atoms[-1]:
            return [(self.K - 1, self.number_scale(n_atoms))]
        out = []
        if z > 0.0:
            out.append((k, z))
        if z < 1.0:
            out.append((k + 1, 1.0 - z))
        return out

    def atom_fraction(self, n_atoms, k, zeta):
        """Share of the entity's atoms that goes to the lower bin k when a number fraction zeta goes there."""
        return float(zeta * self.n_atoms[k] / n_atoms)


# ---- settling and collision kernels ---------------------------------------------------------------------------

def epstein_speed(r_p, rho_p, rho_gas, T, g):
    """Terminal fall speed (m/s) in the free-molecule (Epstein) regime: rho_p g r_p / (rho_gas cbar (1 + pi/8)).

    cbar is the mean thermal speed of air; (1 + pi/8) is the diffuse-reflection-with-accommodation drag factor."""
    cbar = math.sqrt(8.0 * K_B * T / (math.pi * AIR_MASS_AMU * AMU))
    return rho_p * g * r_p / (rho_gas * cbar * (1.0 + math.pi / 8.0))


def _reduced_mass(m1, m2):
    if math.isinf(m1):
        return m2
    if math.isinf(m2):
        return m1
    return m1 * m2 / (m1 + m2)


def fm_kernel(r1, m1, r2, m2, T):
    """Free-molecule collision kernel (m3/s) pi (r1 + r2)^2 sqrt(8 k_B T / (pi mu)); masses in kg, math.inf allowed."""
    return math.pi * (r1 + r2) ** 2 * math.sqrt(8.0 * K_B * T / (math.pi * _reduced_mass(m1, m2)))


def capture_rate(a, b, T):
    """Formation collision rate coefficient (m3/s) for species a and b: max(k_DD, k_HS) (decision 7)."""
    mu = _reduced_mass(a.mass_amu * AMU, b.mass_amu * AMU)
    vbar = math.sqrt(8.0 * K_B * T / (math.pi * mu))
    sigma_ab = 0.5 * (a.sigma_A + b.sigma_A) * ANGSTROM
    k_hs = math.pi * sigma_ab ** 2 * vbar
    da, db = a.dipole_D, b.dipole_D
    if math.isnan(da) or math.isnan(db) or da == 0.0 or db == 0.0:
        return k_hs
    c6 = (2.0 / 3.0) * (da * DEBYE * db * DEBYE / (4.0 * math.pi * epsilon_0)) ** 2 / (K_B * T)
    k_dd = math.pi * gamma(2.0 / 3.0) * vbar * (2.0 * c6 / (K_B * T)) ** (1.0 / 3.0)
    return max(k_dd, k_hs)


def spin_factor(a, b):
    """g_P / (g_A g_B): g_P is 1 if the pair's total electron count is even and 2 if odd (decision 7)."""
    g_p = 1 if (a.electrons + b.electrons) % 2 == 0 else 2
    return g_p / (a.multiplicity * b.multiplicity)


@dataclass(frozen=True)
class Pathway:
    pid: str
    a: str
    b: str


def formation_pathways():
    """All 28 unordered pairs of aluminium-bearing species, self-pairs included, pid "A+B" in AL_SPECIES order."""
    return [Pathway(f"{a}+{b}", a, b) for a, b in combinations_with_replacement(AL_SPECIES, 2)]


# ---- background meteoric smoke --------------------------------------------------------------------------------

class SmokeBackground:
    """Megner et al. (2006) smoke size distribution (data/wake_column/megner2006_smoke.csv) as an uptake sink."""

    def __init__(self, z_km, r_m, dN_dlnr_m3, path="<memory>"):
        self.z_km = np.asarray(z_km, dtype=float)               # (nz,)
        self.r_m = np.asarray(r_m, dtype=float)                 # (nbin,) radii in metres
        self.dN_dlnr = np.asarray(dN_dlnr_m3, dtype=float)      # (nz, nbin) m^-3
        self.path = path
        if self.dN_dlnr.shape != (len(self.z_km), len(self.r_m)):
            raise ValueError(f"smoke table {path} is not a rectangular altitude x radius grid")
        if np.any(np.diff(self.z_km) <= 0):
            raise ValueError(f"altitudes in {path} are not increasing")

    @classmethod
    def load(cls, path=None):
        path = path or SMOKE_CSV
        with open(path, newline="") as f:
            rows = [(float(r["z_km"]), float(r["r_nm"]), float(r["dN_dlnr_m3"]))
                    for r in csv.DictReader(line for line in f if not line.startswith("#"))]
        zs = sorted({z for z, _, _ in rows})
        radii = sorted({r for _, r, _ in rows})
        grid = {(z, r): n for z, r, n in rows}
        if len(grid) != len(zs) * len(radii):
            raise ValueError(f"smoke table {path} does not hold every radius at every altitude")
        return cls(zs, np.array(radii) * 1e-9, [[grid[(z, r)] for r in radii] for z in zs], path)

    def rows(self):
        """(z_km, r_nm, dN_dlnr_m3) for every table row."""
        return [(float(z), float(r) * 1e9, float(self.dN_dlnr[i, j]))
                for i, z in enumerate(self.z_km) for j, r in enumerate(self.r_m)]

    def number_per_bin(self, z_km):
        """Smoke particles per m3 in each radius bin at z_km (log-linear in altitude, bin by bin)."""
        z = float(z_km)
        zmin, zmax = float(self.z_km[0]), float(self.z_km[-1])
        if not (zmin <= z <= zmax):         # also false for nan
            raise ValueError(f"altitude {z_km} km outside {self.path} ({zmin}-{zmax} km)")
        return np.array([_interp_log(z, self.z_km, self.dN_dlnr[:, j]) for j in range(len(self.r_m))]) * SMOKE_DLNR

    def uptake_coefficient(self, sp, z_km, T, density_kgm3):
        """Loss rate (1/s) of one molecule of species sp to the smoke, sticking probability 1.

        Sum over radius bins of n_bin pi (r_bin + r_mol)^2 sqrt(8 k_B T / (pi mu)), with mu the reduced mass of the
        molecule and one smoke particle of that bin: a sphere of radius r_bin and density density_kgm3, mass
        (4/3) pi r_bin^3 density. The smallest bin (0.2 nm) weighs 40-60 amu, as much as AlOH, so the smoke is not
        treated as infinitely heavy. r_mol = sigma_A / 2."""
        if not density_kgm3 > 0.0:
            raise ValueError(f"smoke density must be positive, got {density_kgm3!r} kg/m3")
        n_bin = self.number_per_bin(z_km)
        m_mol = sp.mass_amu * AMU
        r_mol = 0.5 * sp.sigma_A * ANGSTROM
        return float(sum(n * fm_kernel(r, 4.0 / 3.0 * math.pi * r ** 3 * density_kgm3, r_mol, m_mol, T)
                         for n, r in zip(n_bin, self.r_m)))
