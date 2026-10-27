"""The reaction table of the wake-column network (data/wake_column/reactions.csv) and its rate evaluation.

Rates are returned in SI: bimolecular m^3/s, termolecular (with the third body M) m^6/s, photolysis 1/s. The
table's literature units (cm^3/s, cm^6/s) are converted once, on load, by folding the factor into A (and into c0
for the log-quadratic form), so Reaction.k never converts.

M stays in the reactant and product tuples so the stoichiometry is visible; `order` counts it. The box model
substitutes the total air number density for M (decision 12).

Each rate is evaluated only inside its validity range [T_min, T_max]: outside it the temperature is clipped to the
nearest limit and, if a RateWarnings collector is passed, one message per reaction is kept (Review Focus 5), so a
run reports a clipped reaction once rather than once per right-hand-side call.
"""
import csv
import math
import os
from dataclasses import dataclass

from reentry_model.wake_column import DATA

REACTIONS_CSV = os.path.join(DATA, "reactions.csv")

_UNIT_FACTOR = {"cm3/s": 1e-6, "cm6/s": 1e-12, "1/s": 1.0, "-": 1.0}
_FORMS = ("arrhenius", "log_quadratic", "photolysis", "threshold")


class RateWarnings:
    """Collects the temperatures at which rates were clipped; keeps the first message per reaction id."""

    def __init__(self):
        self._first = {}        # rid -> message
        self.counts = {}        # rid -> number of clipped evaluations

    def add(self, rid, T, T_min, T_max):
        self.counts[rid] = self.counts.get(rid, 0) + 1
        if rid not in self._first:
            self._first[rid] = (f"{rid}: T = {T:.1f} K is outside the validity range [{T_min:g}, {T_max:g}] K; "
                                f"the rate is evaluated at the nearest limit (reported once per reaction)")

    @property
    def messages(self):
        return list(self._first.values())


@dataclass(frozen=True)
class Reaction:
    rid: str
    reactants: tuple
    products: tuple
    form: str               # arrhenius | log_quadratic | photolysis | threshold
    order: int              # number of reactants, counting M
    A: float                # SI; nan when not used by the form
    n: float
    Ea_K: float
    c0: float               # SI shift already included
    c1: float
    c2: float
    J: float                # photolysis rate, 1/s
    E0_kJmol: float         # threshold energy of a hot-start reaction
    P0: float               # reaction probability above the threshold
    T_min: float
    T_max: float
    stage: str              # box | hot | both
    day_only: bool
    validity_source: str
    source: str

    def k(self, T, warn=None):
        """Rate coefficient in SI at temperature T (K), with T clipped to [T_min, T_max]."""
        if self.form == "threshold":
            raise ValueError(f"{self.rid} is a threshold (hot-start) reaction and has no thermal rate; "
                             "Stage A evaluates it from E0_kJmol and P0")
        if self.form == "photolysis":
            return self.J                       # temperature-independent, so nothing to clip
        Tc = min(max(T, self.T_min), self.T_max)
        if Tc != T and warn is not None:
            warn.add(self.rid, T, self.T_min, self.T_max)
        if self.form == "arrhenius":
            return self.A * (Tc / 300.0) ** self.n * math.exp(-self.Ea_K / Tc)
        lg = math.log10(Tc)                     # log_quadratic
        return 10.0 ** (self.c0 + self.c1 * lg + self.c2 * lg * lg)


def _float(s):
    s = s.strip()
    return float("nan") if s == "" else float(s)


def load_reactions(path=None, stage=("box", "both")):
    """Return the reactions whose stage is in `stage`, in table order, with rates converted to SI."""
    path = path or REACTIONS_CSV
    out, seen = [], set()
    with open(path, newline="") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if r["rid"] in seen:
                raise ValueError(f"duplicate reaction id {r['rid']!r} in {path}")
            seen.add(r["rid"])
            if r["form"] not in _FORMS:
                raise ValueError(f"{r['rid']}: unknown form {r['form']!r}")
            fac = _UNIT_FACTOR[r["units"]]
            reactants, products = tuple(r["reactants"].split("+")), tuple(r["products"].split("+"))
            if int(r["order"]) != len(reactants):
                raise ValueError(f"{r['rid']}: order {r['order']} does not match its {len(reactants)} reactants")
            if r["stage"] not in ("box", "hot", "both"):
                raise ValueError(f"{r['rid']}: unknown stage {r['stage']!r}")
            if r["stage"] not in stage:
                continue
            c0 = _float(r["c0"])
            out.append(Reaction(
                rid=r["rid"], reactants=reactants, products=products, form=r["form"], order=int(r["order"]),
                A=_float(r["A"]) * fac, n=_float(r["n"]), Ea_K=_float(r["Ea_K"]),
                c0=c0 + math.log10(fac), c1=_float(r["c1"]), c2=_float(r["c2"]), J=_float(r["J"]),
                E0_kJmol=_float(r["E0_kJmol"]), P0=_float(r["P0"]), T_min=float(r["T_min"]),
                T_max=float(r["T_max"]), stage=r["stage"], day_only=bool(int(r["day_only"])),
                validity_source=r["validity_source"], source=r["source"]))
    return out
