"""The species table of the wake-column network, read from data/wake_column/species.csv.

Masses are kept in atomic mass units and diameters in angstrom as in the table; consumers convert with AMU and
ANGSTROM (SI) where they need kilograms and metres.
"""
import csv
import os
from dataclasses import dataclass

from reentry_model.wake_column import DATA

AMU = 1.66053906660e-27     # kg per atomic mass unit
ANGSTROM = 1e-10            # m per angstrom
DEBYE = 3.33564e-30         # C m per debye

AL_SPECIES = ("Al", "AlO", "OAlO", "OAlO2", "AlCO3", "AlOH", "AlOH2")     # AlOH2 is Al(OH)2
AMBIENT_SPECIES = ("N2", "O2", "O", "O3", "H", "H2", "H2O", "CO", "CO2")

LABEL = {"Al": "Al", "AlO": "AlO", "OAlO": "OAlO", "OAlO2": "OAlO$_2$", "AlCO3": "AlCO$_3$", "AlOH": "AlOH",
         "AlOH2": "Al(OH)$_2$", "N2": "N$_2$", "O2": "O$_2$", "O": "O", "O3": "O$_3$", "H": "H", "H2": "H$_2$",
         "H2O": "H$_2$O", "CO": "CO", "CO2": "CO$_2$"}

SPECIES_CSV = os.path.join(DATA, "species.csv")


@dataclass(frozen=True)
class Species:
    name: str
    mass_amu: float
    n_Al: int
    n_O: int
    n_H: int
    deliver_O: int          # oxygen atoms a molecule adds to a particle on sticking (decision 9)
    deliver_H: int          # hydrogen atoms a molecule adds to a particle on sticking
    dipole_D: float         # debye; nan = unknown
    multiplicity: int       # 2S+1 of the ground state
    electrons: int          # total electron count
    sigma_A: float          # hard-sphere collision diameter, angstrom
    ambient: bool
    source: str


def _float(s):
    s = s.strip()
    return float("nan") if s == "" else float(s)


def load_species(path=None):
    """Return {name: Species}, in table order."""
    path = path or SPECIES_CSV
    with open(path, newline="") as f:
        rows = csv.DictReader(line for line in f if not line.startswith("#"))
        out = {}
        for r in rows:
            sp = Species(name=r["name"], mass_amu=float(r["mass_amu"]), n_Al=int(r["n_Al"]), n_O=int(r["n_O"]),
                         n_H=int(r["n_H"]), deliver_O=int(r["deliver_O"]), deliver_H=int(r["deliver_H"]),
                         dipole_D=_float(r["dipole_D"]), multiplicity=int(r["multiplicity"]),
                         electrons=int(r["electrons"]), sigma_A=float(r["sigma_A"]),
                         ambient=bool(int(r["ambient"])), source=r["source"])
            if sp.name in out:
                raise ValueError(f"duplicate species {sp.name!r} in {path}")
            out[sp.name] = sp
    return out
