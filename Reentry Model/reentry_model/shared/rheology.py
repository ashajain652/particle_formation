"""The viscosity of liquid and slurry AA7075 (the Step 4 plan's section 4, docs/superpowers/plans/2026-09-27-deformation-
ring-shedding.md; the layered-skin design section 6.4; the near-liquidus viscosity design of 2026-11-12, docs/superpowers/
specs/2026-11-12-near-liquidus-viscosity-design.md). Shared physics: it imports neither architecture's modules, so that
both architectures call the same law through `law_viscosity` (Asha, 2026-11-12).

Two laws, chosen by name (LAWS). DEFAULT_LAW = "composite" is the run command's default; "bridge" is the library default
of `SkinSettings` and the law every state saved before the composite law was added keeps.

THE COMPOSITE LAW (`composite_viscosity`; the design's section 3), a function of the solid fraction f_s alone, with no
shear rate and no history:

- Liquid, f_s = 0: the material's own viscosity mu_l (1.3 mPa s for AA7075; Assael et al. 2006, J. Phys. Chem. Ref. Data
  35, 285-300, Eq. 19: 1.42 mPa s for pure aluminium at 908 K, known to 14 %, extrapolated 25 K).
- Below LI_SOLID_MIN = 10 % solid, the grains as a suspension in their envelope fraction: the Mori-Ototake law (Moon's
  form, S = SUSPENSION_S = 3 for spheres), eta = mu_l [1 + S / (1/phi - 1/phi_m)], phi = f_s / f_i, phi_m =
  SUSPENSION_PHI_M = 0.55, with f_i = INTERNAL_SOLID_FRACTION = 0.5 the share of each grain's envelope that is solid.
  Sources: Moon 1990 (MIT thesis, pp. 101-102, Eq. 51, Fig. 5.14) and Moon, Cornie and Flemings 1991 (Mater. Sci. Eng.
  A144, 253-265), measured on silicon carbide in liquid aluminium, read in full; the analogues (Ishibashi and Sato, basalt
  to 31 % crystals; AZ91D Couette; ice slurries), measured on other systems, abstract level; the envelope fraction from
  Beckermann and Wang 1997 (JOM) and Wu and Ludwig 2009 (Acta Mater. 57, 5621-5631), the solidification models'
  construction; f_i = 0.5 from Zheng et al. 2018 (Leoben) and assumed for 7075, bracketed by 1 (compact grains) and 0.25
  in the design's tests. The nearest aluminium values, Hirai et al. 1993 (ISIJ Int. 33, Al-10Cu: 1.08, 1.22, 1.49 times
  the liquid at 2, 5 and 10 % solid at 200 1/s, and 2.5, 5.1 and 10.3 times at 75 1/s), are inferred from fitted curves
  and not used in the law. The law gives
  at most 1.94 times the liquid at f_s just below 0.1 (1.37 with f_i = 1, 5.4 with f_i = 0.25).
- From LI_SOLID_MIN = 0.1 up to RIGID_SOLID_FRACTION = 0.5: Li Yageng et al. (2014, China Foundry 11, 79-84, Table 3,
  pp. 82-83), their fit at the highest rate, 489.88 1/s, where their slurry had its finest and roundest grains,
  eta = LI_BROKEN_A exp(LI_BROKEN_B f_s) = 0.04346 exp(3.23 f_s) Pa s, taken at every shear rate (their method, pp. 79-81,
  and the structures of their Fig. 11). It is fitted to 7075 data
  at 10, 20, 30, 40 and 50 % solid, read in full: 46.2 times the liquid at 10 % solid and 168 times at 50 %. That the
  structure in flight is the most broken-down one rests on Yurko and Flemings (A357), via Kirkwood and Ward 2004, and Liu
  et al. 2003 (Sn-15Pb), for the breakdown, and Florides et al. 2024 (A356, about 520 s) for the rebuild: assumed for
  7075, measured on other alloys, abstract level except Kirkwood and Ward. That the viscosity at a fixed structure does
  not depend on the rate rests on Kirkwood and Ward 2004, Atkinson and Favier 2016 and Koke and Modigell 2003: assumed,
  abstract level except Koke and Modigell. Dropping the history and the falling branch of the flow curve rests on
  Divoux et al. 2016 (the falling branch "not found in experiment"), Yerushalmi, Katz and Shinnar 1970 (the instability,
  via Divoux et al.), Coussot et al. 2002 and Moller, Mewis and Bonn 2006, and on Li et al.'s method, one slurry per
  rate: a consensus in soft matter, with no semi-solid metal measurement of it. Li et al.'s data are placed by solid
  fraction, not temperature (their Table 2, pp. 80-81), their 10 % solid (903.45 K by their DSC) being the model's
  906.36 K on its Scheil curve (906.4 K, BRIDGE_LOW, is where that curve stands at 9.75 % solid): the builder's call.
- Rigid, at or above RIGID_SOLID_FRACTION = 0.5: infinite, as before.
- THE JUMP AT 10 % SOLID IS UNRESOLVED. The law jumps from the suspension value to Li et al.'s: from 1.94 to 46.2 times
  the liquid with f_i = 0.5 (a factor of 24), from 1.37 with f_i = 1 (34), from 5.40 with f_i = 0.25 (8.6). No data say
  how the two join; Li et al.'s 7075 at 10 % solid is about a hundred times Hirai et al.'s aluminium-copper, which the
  reviews attribute to different structures. On the model's Scheil curve it sits at 906.36 K.

What is assumed and what is uncertain (the design's section 3.7, numbered as there; "abstract level" means the review
read only an abstract or a search snippet).

At 10 % solid and above:
  1. Assumed: the structure in flight is Li et al.'s most broken-down one. Breakdown in milliseconds (A357 fell to
     6 Pa s in about 4 ms at 800 1/s, Yurko and Flemings, in Kirkwood and Ward 2004) and rebuilding over about 500 s
     (A356's fitted rebuild rate 0.00191 1/s, Florides et al. 2024) are extrapolated to 7075 and to the gas's rates.
     Metal that has never been sheared (remelting from below, or a re-frozen deposit) would be stiffer, and the law
     does not represent it.
  2. Assumed: no rate dependence at a fixed structure. Structural models are Newtonian at fixed structure (Kirkwood
     and Ward 2004, citing Chen and Fan); rate-jump tests show the stress rising in proportion to the rate or faster
     (Atkinson and Favier 2016, abstract level; Koke and Modigell 2003).
  3. Upper bound above 489.88 1/s. In steady tin-lead data the viscosity was still falling at 2000 1/s (Turng and Wang,
     in Fan 2002, Fig. 4). A 7075 layer sheared faster than 490 1/s may break down further than this law allows.
  4. Li et al.'s lowest values are uncertain. At 10 % solid and 490 1/s their plotted points scatter from 0.05 to
     0.2 Pa s, and their rig was not calibrated with a known liquid.
  5. Different processing. Li et al.'s slurry was solidified while stirred; the sphere remelts wrought metal.
Below 10 % solid (no aluminium measurement exists):
  1. Nothing has been measured. No rheological measurement of any aluminium alloy between the liquid and 10 % solid was
     found; the nearest are fitted curves and analogue systems.
  2. The internal solid fraction is unmeasured for 7075. 0.5 is the solidification models' default. The only measured
     values, 0.1 to 0.55, are for a transparent analogue, succinonitrile-acetone (Badillo, Ceynar and Beckermann). The
     review infers about 1 for remelting wrought grains and 0.2 to 0.3 for re-frozen deposit. The bracket rows of
     section 3.6 show the effect.
  3. No rate dependence is assumed. Hirai et al.'s aluminium-10 % copper fits are strongly rate dependent for branched
     crystals: 5.1 times the liquid at 75 1/s against 1.22 times at 200 1/s, at 5 % solid (values inferred from fitted
     curves, not data points). At the gas's rates, of order 10^3 1/s and above, the sheared values are the relevant ones.
  4. Coherency may come earlier in coarse deposits. Coherency, where grains first form a connected network, ranges from
     below 10 % solid in poorly grain-refined aluminium alloys to above 50 % in fine-grained ones (Arnberg, Chai and
     Backerud 1993; Chai et al. 1995; Malekan and Shabestari 2009, A319, 16-21 %). For 7075 it was measured at 55 to 85 %
     by a thermal method that reads high (Ahmad, Naher and Brabazon 2013). A coarse re-frozen deposit could be coherent
     below 10 % solid, which this law cannot represent.
  5. The jump at 10 % solid is unresolved (above).
  6. The layer's position is uncertain by more than its width. Measured 7075 liquidus temperatures span 903 to 913 K,
     about six times the 1.6 K the model's Scheil curve gives the 0-10 % solid layer.
  7. The liquid itself is uncertain. Its viscosity is known to 14 % for pure aluminium (Assael et al. 2006), alloying
     moves it by up to 25 % (Wang and Overfelt 2002), and no 7xxx liquid has been measured.
  8. S = 3 comes from rigid particles. It was measured on silicon carbide in liquid aluminium; Moon's ratios are taken
     against his rig's 0.020 Pa s "liquid", 15-16 times the true value, so they hold only if that error is proportional.
  9. The law is a steady-state law. The stress just after a sudden jump in rate rises as the rate to the power 2.07
     (Koke and Modigell 2003, Fig. 18), which is not included.

THE BRIDGE LAW (`viscosity`; the present law, unchanged):

- Liquid, at and above the liquidus: the material's own viscosity (1.3 mPa s for AA7075).
- Bridge, BRIDGE_LOW to BRIDGE_HIGH (906.4-908 K): log-linear in the temperature from Li et al.'s value at BRIDGE_LOW,
  where Scheil's curve stands at 10 % solid, to the liquid at the liquidus (labelled: no data).
- Slurry, from the bridge down to half solid: Li Yageng et al. (2014, China Foundry 11, 79-84),
  eta = [0.871 - 0.00849 gamma^0.74924] exp(3.7311 f_s) Pa s, measured at f_s = 0.1-0.5 and shear rates gamma of
  61-490 1/s on continuously sheared slurry; gamma is held at no more than GAMMA_MAX = 367 1/s, the highest rate at
  which the fit still matches the data (its bracket would reach zero at 483 1/s). Between zero and 61 1/s it is
  extrapolated, and there it is the runny side of the data (1.4-1.9 times low at 61 1/s; the Step 4 plan's fact 10).
  From 0.87 Pa s (f_s = 0) to 5.6 Pa s (f_s = 0.5) at rest, 670 to 4 300 times the liquid.
- Rigid, at or above half solid: infinite. A rigid cell does not flow; it is anchored (spec 6.4)."""
import numpy as np


LI_A = 0.871                       # Pa s
LI_B = 0.00849                     # Pa s^(1 + LI_N)
LI_N = 0.74924
LI_C = 3.7311
GAMMA_MAX = 367.0                  # 1/s: Li et al.'s law is held at this shear rate
BRIDGE_LOW, BRIDGE_HIGH = 906.4, 908.0     # K: the log-linear bridge from the slurry to the liquid
RIGID_SOLID_FRACTION = 0.5         # at or above it a cell is rigid (material.RIGID_LIQUID_FRACTION's complement)

LI_BROKEN_RATE = 489.88            # 1/s: Li et al.'s highest rate, whose fit the composite law takes at every rate
LI_BROKEN_A = 0.04346              # Pa s: that fit's eta = A exp(B f_s), Table 3
LI_BROKEN_B = 3.23
LI_SOLID_MIN = 0.1                 # the composite law is the suspension below it and Li et al.'s from it
SUSPENSION_S = 3.0                 # Mori-Ototake's S for spheres (Moon)
SUSPENSION_PHI_M = 0.55            # the suspension's packing fraction, in the envelope fraction
INTERNAL_SOLID_FRACTION = 0.5      # ASSUMED for 7075 (the solidification models' default; no measurement)
LAWS = ("composite", "bridge")
DEFAULT_LAW = "composite"          # the run command's default; SkinSettings' library default is "bridge"


def li_viscosity(f_s, gamma):
    """Li et al.'s slurry viscosity [Pa s] at solid fraction f_s and shear rate gamma [1/s], the rate held within
    [0, GAMMA_MAX]."""
    g = np.clip(np.asarray(gamma, dtype=float), 0.0, GAMMA_MAX)
    return (LI_A - LI_B * g ** LI_N) * np.exp(LI_C * np.asarray(f_s, dtype=float))


def viscosity(T, f_s, gamma, mu_liquid, f_s_bridge):
    """The viscosity [Pa s] of metal at temperature T [K] and solid fraction f_s, sheared at gamma [1/s] (module
    docstring); `mu_liquid` [Pa s] is the liquid's and `f_s_bridge` the material's solid fraction at BRIDGE_LOW (0.098 on
    AA7075's Scheil curve). Infinite where the metal is rigid."""
    T, f_s = np.asarray(T, dtype=float), np.asarray(f_s, dtype=float)
    top = li_viscosity(f_s_bridge, gamma)
    w = np.clip((T - BRIDGE_LOW) / (BRIDGE_HIGH - BRIDGE_LOW), 0.0, 1.0)
    bridge = np.exp((1.0 - w) * np.log(top) + w * np.log(mu_liquid))
    eta = np.where(T >= BRIDGE_HIGH, mu_liquid, np.where(T > BRIDGE_LOW, bridge, li_viscosity(f_s, gamma)))
    return np.where(f_s >= RIGID_SOLID_FRACTION, np.inf, eta)


def _internal(internal_fraction):
    """The internal solid fraction, read from the module at call time when None; it must exceed LI_SOLID_MIN /
    SUSPENSION_PHI_M (so that the envelope fraction stays below the packing fraction under 10 % solid) and not exceed 1."""
    f_i = INTERNAL_SOLID_FRACTION if internal_fraction is None else internal_fraction
    if not LI_SOLID_MIN / SUSPENSION_PHI_M < f_i <= 1.0:
        raise ValueError("internal_fraction must lie in ({:.4f}, 1], got {}".format(LI_SOLID_MIN / SUSPENSION_PHI_M, f_i))
    return float(f_i)


def suspension_viscosity(f_s, mu_liquid, internal_fraction=None):
    """The grains as a suspension in their envelope fraction phi = f_s / f_i [Pa s]: the Mori-Ototake law
    mu_liquid [1 + S / (1/phi - 1/phi_m)], mu_liquid at f_s <= 0 and infinite where phi >= SUSPENSION_PHI_M."""
    f_i = _internal(internal_fraction)
    phi = np.asarray(f_s, dtype=float) / f_i
    with np.errstate(divide="ignore", invalid="ignore"):
        grown = mu_liquid * (1.0 + SUSPENSION_S / (1.0 / phi - 1.0 / SUSPENSION_PHI_M))
    return np.where(phi <= 0.0, float(mu_liquid), np.where(phi >= SUSPENSION_PHI_M, np.inf, grown))


def composite_viscosity(f_s, mu_liquid, internal_fraction=None):
    """The composite law's viscosity [Pa s] at solid fraction f_s (module docstring): the suspension below
    LI_SOLID_MIN, Li et al.'s broken-down fit from there to RIGID_SOLID_FRACTION, infinite at and above it. The jump at
    LI_SOLID_MIN is unresolved."""
    f = np.asarray(f_s, dtype=float)
    sus = suspension_viscosity(f, mu_liquid, internal_fraction)
    li = LI_BROKEN_A * np.exp(LI_BROKEN_B * f)
    return np.where(f >= RIGID_SOLID_FRACTION, np.inf, np.where(f >= LI_SOLID_MIN, li, sus))


def law_viscosity(law, T, f_s, gamma, mu_liquid, f_s_bridge):
    """The viscosity [Pa s] under the named law (LAWS): "bridge" is `viscosity`, unchanged; "composite" is
    `composite_viscosity` of f_s, broadcast over the shapes of T, f_s and gamma (as "bridge" is), and ignores gamma."""
    if law == "bridge":
        return viscosity(T, f_s, gamma, mu_liquid, f_s_bridge)
    if law == "composite":
        eta = composite_viscosity(f_s, mu_liquid)
        return np.array(np.broadcast_to(eta, np.broadcast_shapes(np.shape(T), np.shape(f_s), np.shape(gamma))), dtype=float)
    raise ValueError("unknown viscosity law {!r}; expected one of {}".format(law, LAWS))


def rate_cap(law):
    """The shear rate [1/s] at which a law's caller holds a recorded rate: GAMMA_MAX under "bridge", infinite under
    "composite" (which does not depend on the rate)."""
    if law == "bridge":
        return GAMMA_MAX
    if law == "composite":
        return float("inf")
    raise ValueError("unknown viscosity law {!r}; expected one of {}".format(law, LAWS))


def bridge_solid_fraction(material):
    """The solid fraction [-] at BRIDGE_LOW on the material's own curve (0.0975 for AA7075's Scheil curve); uses only
    the object's `liquid_fraction`."""
    return 1.0 - float(material.liquid_fraction(np.array([BRIDGE_LOW]))[0])
