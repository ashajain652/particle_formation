"""Plume-radius laws and the choice between them, from the parent's own flow state.

The plume behind a parent body is described by a radius r_b(x) at a distance x behind the point where the parent
emitted it. Three laws are implemented, and which applies is decided per emission point from flow quantities only
(Knudsen number, Reynolds number, Taylor's pressure parameter); nothing here takes or knows an altitude, so the same
code serves any parent-follower pair (design spec section 4.2).

Sources
  Taylor, Melcher and Washburn (1964), AIAA J. 2(10) 1731, equations 4 and 5 (printed page 1737): turbulent wake core,
      w / sqrt(C_D A) = K [x / sqrt(C_D A)]^(1/3), K = 0.66, radius = w / 2. Fitted to 2 <= x / sqrt(C_D A) <= 1200 and
      3.7-6.2 km/s. The regime division by the pressure parameter p sqrt(C_D A) in cm Hg-cm is from p. 1735.
      Notes: citations/wake_radius_growth_laws.md section 0 (in the main checkout, not copied here).
  Park (1971), NASA TM X-62,052, equation 11: turbulent wake, r_b / r0 = C (s' / r0)^n, C = 0.7, n = 1/3, with s' the
      distance from a virtual origin placed so that r_b = r0 at the parent (s' = s + s0, s0 = (1/C)^(1/n) r0).
      Valid to s about 30 km, 5.8-7.9 km/s. C and n were reconstructed from a garbled scan (spec section 4.2).
  Feldman (1961), J. Aerospace Sci. 28:433, equations 6, 8, 11 (p. 442) and 14 (p. 443): laminar conduction-controlled
      trail. The combined width law is not printed in the paper; it is assembled in
      citations/laminar_wake_growth_laws.md section 1:  R^2 = R0^2 + 4 kappa x / (f U),
      kappa = mu(T_axis) / (rho Pr), with Sutherland's viscosity (Feldman Eq. 13).
      Feldman's validity: 15,000-35,000 ft/s and 100,000-250,000 ft altitude for hemisphere-cylinders. We never test
      altitude. The speed range is the named constant FELDMAN_SPEED_RANGE. The notes give no numeric body-size
      range; the only statement (Feldman p. 441: below 3 or 30 cm the expansion and conduction regions must be
      computed together) is used as the lower diameter limit FELDMAN_DIAMETER_RANGE, with no upper limit.

R is Feldman's density-weighted (Howarth) 1/e radius of the temperature excess, in which hot, low-density gas counts
for less than its physical volume, so the physical radius is larger wherever the gas is hot. This first version uses R
as the plume radius r_b anyway (a call recorded in the design spec section 8). Taylor's r_b is the edge of the
turbulent core, so "inside the plume" means slightly different things in the two regimes; `laminar_edge_factor`
multiplies Feldman's R to let that be varied.

Units: SI throughout. Taylor's pressure parameter is in cm Hg-cm; the conversion lives only in `pressure_parameter`.
"""
from dataclasses import dataclass, field, fields
import math

FT = 0.3048                                       # metres per foot
FELDMAN_SPEED_RANGE = (15000.0 * FT, 35000.0 * FT)   # m/s: 4.57-10.67 km/s (Feldman 1961 abstract, p. 433)
FELDMAN_DIAMETER_RANGE = (0.03, math.inf)         # m: lower limit from Feldman p. 441 (see module docstring)
PA_TO_CMHG = 1.0 / 1333.22                        # one cm Hg is 1333.22 Pa
M_TO_CM = 100.0

REGIMES = ("rarefied", "too_slow", "laminar", "transitional", "turbulent")
TURBULENT_LAWS = ("taylor", "park")
FALLBACK_LAWS = (None, "taylor", "park", "feldman_low", "feldman_high")


@dataclass(frozen=True)
class EmissionState:
    """The parent's flow state at one emission time (all SI)."""
    rho: float      # ambient density, kg/m3
    T: float        # ambient temperature, K
    p: float        # ambient pressure, Pa
    mu: float       # ambient dynamic viscosity, Pa s
    mfp: float      # ambient mean free path, m
    U: float        # air-relative speed of the parent, m/s
    d: float        # parent equal-area diameter, m
    cda: float      # drag coefficient times frontal area C_D A, m2
    r0: float       # equal-area radius, m


@dataclass(frozen=True)
class LawSettings:
    """Named, overridable thresholds and law constants (design spec section 4.2)."""
    kn_continuum: float = 0.01                    # Kn above this: rarefied, no law
    re_crit: float = 5e4                          # Re_d threshold (Lees 1964 p. 418)
    p_lam: float = 0.5                            # cm Hg-cm: laminar below this (Taylor 1964 p. 1735)
    p_turb: float = 2.0                           # cm Hg-cm: turbulent above this
    turbulent_law: str = "taylor"                 # "taylor" or "park"
    taylor_K: float = 0.66
    taylor_n: float = 1.0 / 3.0
    taylor_xmax_over_L: float = 1200.0            # fit limit of x / sqrt(C_D A)
    taylor_speed_range: tuple = (3.7e3, 6.2e3)    # m/s
    park_C: float = 0.7
    park_n: float = 1.0 / 3.0
    s_park_limit: float = 30e3                    # m
    park_speed_range: tuple = (5.8e3, 7.9e3)      # m/s
    feldman_f: float = 0.8                        # hot-core speed over flight speed
    feldman_Pr: float = 0.71
    T_axis_hot: float = 3000.0                    # K: axis temperature of the upper Feldman bound
    laminar_edge_factor: float = 1.0              # multiplies Feldman's R
    rarefied_fallback: object = None              # None, "taylor", "park", "feldman_low" or "feldman_high"
    u_min_parent: float = 1.0                     # m/s: a slower parent leaves no plume

    def __post_init__(self):
        if self.turbulent_law not in TURBULENT_LAWS:
            raise ValueError(f"turbulent_law must be one of {TURBULENT_LAWS}, got {self.turbulent_law!r}")
        if self.rarefied_fallback not in FALLBACK_LAWS:
            raise ValueError(f"rarefied_fallback must be one of {FALLBACK_LAWS}, got {self.rarefied_fallback!r}")

    def non_default(self) -> dict:
        """The fields that differ from the defaults, for run names."""
        default = LawSettings()
        return {f.name: getattr(self, f.name) for f in fields(self)
                if getattr(self, f.name) != getattr(default, f.name)}


@dataclass(frozen=True)
class RadiusResult:
    """A plume radius (NaN when no law applies), the regime, the laws that produced it and diagnostics."""
    r: float
    regime: str
    laws: tuple = ()
    flags: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------------------------- flow quantities

def sutherland_mu(T: float) -> float:
    """Dynamic viscosity of air, Pa s (Feldman Eq. 13: 1.72e-5 (T/273)^1.5 386/(113+T))."""
    return 1.72e-5 * (T / 273.0) ** 1.5 * 386.0 / (113.0 + T)


def knudsen(s: EmissionState) -> float:
    """Mean free path over parent diameter."""
    return s.mfp / s.d


def reynolds(s: EmissionState) -> float:
    """Free-stream Reynolds number on the parent diameter, rho U d / mu."""
    return s.rho * s.U * s.d / s.mu


def pressure_parameter(s: EmissionState) -> float:
    """Taylor's pressure parameter p sqrt(C_D A) in cm Hg-cm (p in Pa, sqrt(C_D A) in m)."""
    return s.p * math.sqrt(s.cda) * PA_TO_CMHG * M_TO_CM


def classify(s: EmissionState, settings: LawSettings) -> str:
    """Regime of the parent's plume at this emission point, from flow quantities only (spec section 4.2)."""
    if s.U < settings.u_min_parent:
        return "too_slow"
    if knudsen(s) > settings.kn_continuum:
        return "rarefied"
    re_high = reynolds(s) >= settings.re_crit
    pp = pressure_parameter(s)
    if not re_high and pp < settings.p_lam:
        return "laminar"
    if re_high and pp > settings.p_turb:
        return "turbulent"
    return "transitional"


# ---------------------------------------------------------------------------------------------------- the laws

def park_radius(s_behind: float, r0: float, C: float, n: float) -> float:
    """Park (1971) Eq. 11 with the virtual origin: equals r0 at the parent (s_behind = 0)."""
    s0 = (1.0 / C) ** (1.0 / n) * r0
    return r0 * C * ((max(s_behind, 0.0) + s0) / r0) ** n


def taylor_radius(x: float, cda: float, r0: float, K: float, n: float) -> float:
    """Taylor et al. (1964) Eq. 4: half the turbulent-core width, never below the parent radius r0."""
    L = math.sqrt(cda)
    return max(0.5 * K * L * (max(x, 0.0) / L) ** n, r0)


def feldman_radius(x: float, s: EmissionState, T_axis: float, f: float, Pr: float, edge_factor: float) -> float:
    """Laminar radius sqrt(r0^2 + 4 kappa x / (f U)) times edge_factor; NaN if there is no flow to carry the heat."""
    if s.U <= 0.0 or f <= 0.0 or s.rho <= 0.0:
        return math.nan
    kappa = sutherland_mu(T_axis) / (s.rho * Pr)
    return edge_factor * math.sqrt(s.r0 ** 2 + 4.0 * kappa * max(x, 0.0) / (f * s.U))


def _in(value: float, rng) -> bool:
    return rng[0] <= value <= rng[1]


def _law_radius(name: str, x: float, s: EmissionState, st: LawSettings) -> float:
    if name == "taylor":
        return taylor_radius(x, s.cda, s.r0, st.taylor_K, st.taylor_n)
    if name == "park":
        return park_radius(x, s.r0, st.park_C, st.park_n)
    T_axis = s.T if name == "feldman_low" else st.T_axis_hot
    return feldman_radius(x, s, T_axis, st.feldman_f, st.feldman_Pr, st.laminar_edge_factor)


def _flags(laws, x: float, s: EmissionState, st: LawSettings) -> dict:
    """Diagnostics for the laws used: recorded, never used to refuse a computation."""
    speed_out = distance_out = False
    for name in laws:
        if name == "taylor":
            speed_out |= not _in(s.U, st.taylor_speed_range)
            distance_out |= x > st.taylor_xmax_over_L * math.sqrt(s.cda)
        elif name == "park":
            speed_out |= not _in(s.U, st.park_speed_range)
            distance_out |= x > st.s_park_limit
        else:
            speed_out |= not _in(s.U, FELDMAN_SPEED_RANGE)
    feldman = any(n.startswith("feldman") for n in laws)
    extrapolated = feldman and not (_in(s.U, FELDMAN_SPEED_RANGE) and _in(s.d, FELDMAN_DIAMETER_RANGE))
    return {"speed_out_of_range": bool(speed_out), "distance_out_of_range": bool(distance_out),
            "feldman_extrapolation": bool(extrapolated)}


def plume_radius(x: float, s: EmissionState, settings: LawSettings, bound: str) -> RadiusResult:
    """Plume radius at distance x behind the emission point; bound is "low" or "high" (the minimum or maximum
    over the laws that apply). r is NaN, with no error, when the regime has no law (rarefied without a fallback,
    or too slow)."""
    if bound not in ("low", "high"):
        raise ValueError(f"bound must be 'low' or 'high', got {bound!r}")
    regime = classify(s, settings)
    if regime == "too_slow":
        laws = ()
    elif regime == "rarefied":
        laws = (settings.rarefied_fallback,) if settings.rarefied_fallback else ()
    elif regime == "laminar":
        laws = ("feldman_low", "feldman_high")
    elif regime == "transitional":
        laws = ("feldman_low", "feldman_high", settings.turbulent_law)
    else:
        laws = (settings.turbulent_law,)
    if not laws:
        return RadiusResult(math.nan, regime, (), {"speed_out_of_range": False, "distance_out_of_range": False,
                                                   "feldman_extrapolation": False})
    radii = [_law_radius(n, x, s, settings) for n in laws]
    r = min(radii) if bound == "low" else max(radii)
    return RadiusResult(r, regime, tuple(laws), _flags(laws, x, s, settings))
