"""The settings of one wake-column run, the case it runs, and the run-name scheme.

Every efficiency and scale that the literature does not fix is a named field here (spec anti-goal: nothing
unmeasured is hard-coded). The run name encodes the case and every non-default physics field, so a sweep can skip a
run whose summary JSON already exists.
"""
from dataclasses import dataclass

HOT_STARTS = ("slowing_down", "droplet_temperature")
AMBIENTS = ("waccm", "msis_afgl")


def _fmt(v):
    """A value for a run name: shortest exact text (the 6-digit %g form when it is exact, else repr)."""
    if isinstance(v, int):
        return str(v)
    g = f"{v:g}"
    return g if float(g) == v else repr(v)


@dataclass(frozen=True)
class Settings:
    # hot start (Stage A)
    hot_start: str = "slowing_down"         # "slowing_down" Monte Carlo, or "droplet_temperature": the same Stage A with
                                            # the droplet speed set to zero (atoms leave at the droplet temperature
                                            # only; the reactions stay on, but collisions seldom reach a threshold)
    hot_reaction_scale: float = 1.0         # multiplier on the hot-stage reaction probabilities
    thermal_cutoff: float = 2.0             # a test particle is thermal once its energy is below this many k_B T
    n_test_particles: int = 10000
    n_emission_samples: int = 20
    seed: int = 20261027
    # spreading
    l_inner_m: float = 30.0                 # inner eddy scale, m: below it eddies do not mix the plume
    eddy_scale: float = 1.0                 # multiplier on the Megner eddy diffusion coefficient
    # efficiencies
    formation_scale: float = 1.0            # multiplier on every gas-phase dimer formation rate
    formation_eff: tuple = ()               # ((pathway id, value), ...) per-pathway formation efficiencies
    sticking_scale: float = 1.0             # multiplier on every molecule-to-particle sticking efficiency
    sticking: tuple = ()                    # ((species, value), ...) per-species sticking efficiencies
    oxidation_eff: float = 1.0              # efficiency of O, O2 and O3 uptake by particles
    smoke_sticking: float = 1.0             # efficiency of uptake onto background meteoric smoke
    # ambient and bins
    ambient: str = "waccm"                  # "waccm", or "msis_afgl" (NRLMSIS + AFGL cross-check)
    bin_ratio: float = 2.0 ** 0.5           # atom-number ratio q between neighbouring particle bins; 2^(1/2), the
                                            # largest of 2, 2^(1/2), 2^(1/3) that keeps the self-preserving width
                                            # within 1.45 +- 0.03 (q = 2 gives 1.506 by numerical diffusion)
    d_max_nm: float = 400.0                 # alumina diameter at which the bins stop, nm
    # time grid and solver
    t_end_s: float = 43200.0                # 12 hours
    t_first_s: float = 1e-6
    points_per_decade: int = 10
    rtol: float = 1e-6
    # test-only switches (not part of the run name)
    chemistry: bool = True
    particles: bool = True
    coagulation: bool = True
    settling: bool = True
    freeze_radius: bool = False
    reactions: tuple = None                 # tuple of rids to keep, or None for all
    smoke: bool = True

    def __post_init__(self):
        if self.hot_start not in HOT_STARTS:
            raise ValueError(f"hot_start must be one of {HOT_STARTS}, got {self.hot_start!r}")
        if self.ambient not in AMBIENTS:
            raise ValueError(f"ambient must be one of {AMBIENTS}, got {self.ambient!r}")

    def tag(self):
        """'_<short><value>' for every field that differs from its default, in a fixed order ('' when all are
        default): first the physics fields (form, stick, ox, hot, hotdroplet, eddy, lin, cut, smk, q, fe-, st-),
        then the numerical, sampling and test-only fields. `ambient` is not here: run_name puts it in the name."""
        d = Settings()
        parts = []
        for short, field in (("form", "formation_scale"), ("stick", "sticking_scale"), ("ox", "oxidation_eff"),
                             ("hot", "hot_reaction_scale")):
            if getattr(self, field) != getattr(d, field):
                parts.append(f"{short}{_fmt(getattr(self, field))}")
        if self.hot_start != d.hot_start:
            parts.append("hotdroplet")              # the only other hot_start allowed
        for short, field in (("eddy", "eddy_scale"), ("lin", "l_inner_m"), ("cut", "thermal_cutoff"),
                             ("smk", "smoke_sticking"), ("q", "bin_ratio")):
            if getattr(self, field) != getattr(d, field):
                parts.append(f"{short}{_fmt(getattr(self, field))}")
        for pid, value in sorted(self.formation_eff):
            parts.append(f"fe-{pid}{_fmt(value)}")
        for name, value in sorted(self.sticking):
            parts.append(f"st-{name}{_fmt(value)}")
        # sampling, time grid and solver
        for short, field in (("ntp", "n_test_particles"), ("nem", "n_emission_samples"), ("seed", "seed"),
                             ("dmax", "d_max_nm"), ("tend", "t_end_s"), ("tfirst", "t_first_s"),
                             ("ppd", "points_per_decade"), ("rtol", "rtol")):
            if getattr(self, field) != getattr(d, field):
                parts.append(f"{short}{_fmt(getattr(self, field))}")
        # test-only switches
        for short, field in (("nochem", "chemistry"), ("nopart", "particles"), ("nocoag", "coagulation"),
                             ("nosettle", "settling")):
            if getattr(self, field) != getattr(d, field):
                parts.append(short)
        if self.freeze_radius != d.freeze_radius:
            parts.append("freeze")
        if self.reactions != d.reactions:
            parts.append("rx-" + ("-".join(sorted(self.reactions)) if self.reactions else "none"))
        if self.smoke != d.smoke:
            parts.append("nosmoke")
        return "".join("_" + p for p in parts)


@dataclass(frozen=True)
class Case:
    d0_um: float            # initial droplet diameter, micrometres
    h0_km: float            # release altitude, km
    time_of_day: str        # "day" or "night"

    def __post_init__(self):
        if self.time_of_day not in ("day", "night"):
            raise ValueError(f"time_of_day must be 'day' or 'night', got {self.time_of_day!r}")


def run_name(case, settings):
    """The run's name, encoding the case, the ambient source and every non-default physics setting."""
    return (f"wake_d{case.d0_um:06.2f}um_h{case.h0_km:05.1f}km_{case.time_of_day}_{settings.ambient}"
            f"{settings.tag()}")
