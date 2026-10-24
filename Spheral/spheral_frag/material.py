"""The material table Spheral and the finite element share (spec §8.1; §8.3's viscosity row; M1 plan, Task 3).

`prepare` builds the table in drama_env by calling the finite-element `Material` (`fe.material_arrays`); this
module evaluates it with numpy alone, so the runner needs no copy of the finite-element package. Two forms live in one file:

- **The exact node table**: the finite-element material's own enthalpy nodes (`h_T`, `h_nodes`, the c_p at each
  interval's ends, the latent slope per interval) and liquid-fraction nodes. `enthalpy`, `temperature` and
  `liquid_fraction` repeat the finite element's arithmetic operation for operation (the piecewise quadratic, its
  quadratic-formula inverse, `np.interp` on the Scheil nodes), so they agree with it to round-off at any temperature,
  not only on a grid -- f_l bitwise, also under a numpy that rounds np.interp differently (`interp`).
- **The 1 K table** spec §8.1 asks for (250-1,500 K): `h300`, `f_l` and the mechanical columns `rho_free`, `E`,
  `nu`, `K`, `G` (from `spheral_frag/data/aa7075_mechanical.json`), for inspection, plots and the free density and
  moduli, which are evaluated by `np.interp` on it.

**Zero of energy.** The finite element's enthalpy is zero at 293 K (`T_ref_fe`); the spec's is zero at 300 K.
`enthalpy` returns h_FE(T) - h_FE(300 K): the finite element's value first, then the stored offset subtracted, so
the shift is exact in differences and `enthalpy(T)` equals the finite element's `enthalpy(T) - enthalpy(300)`
bitwise. `temperature(h)` adds the offset back and inverts exactly as the finite element does.

**Liquid fraction kinds** (`fl_kind`): 1 = Scheil, `np.interp(T, fl_T, fl)` (Step 3's tabulated curve, eutectic ramp
at the solidus); 2 = linear range, `clip((T - T_s)/(T_l - T_s), 0, 1)` written as the finite element writes it
(np.interp on two nodes rounds differently); 0 = no melting, identically zero.

**Viscosity** (spec §8.3; Step 4 plan §4; decision 7 of the M1 review): Li et al. (2014)
eta = [0.871 - 0.00849 g^0.74924] exp(3.7311 f_s) Pa s for 0.1 <= f_s <= 0.5, the shear rate g held at no more than
367 1/s; a log-linear bridge in T from its f_s = 0.1 value at `T_bridge` (where f_l = 0.9 on the table's own curve:
906.37 K for AA7075_scheil, the spec's "906.4 K") to the liquid's viscosity at the liquidus (908 K); the liquid's
1.3 mPa s above; inf below 50 % liquid (no viscous row: the solid and mush rows of §8.3 are M2's).

**Flow stress** (spec §8.3; M2 plan Task 3, fact 10, decisions at review 1 and 4): `flow_stress(table, T, rate, law)`.
Chen G. et al. (2016, J. Alloys Compd. 674, 26-36) Eq. 20, sigma in MPa inside the law and Pa outside, R = 8.314:
plastic branch rate = 6.75872e11 sinh(0.014677 sigma)^6.85385 exp(-186163/RT); semi-solid branch
(1 - 2 f_l)^15.81342 rate = 1.28926e8 sinh(0.03006 sigma)^6.521135 exp(-175131/RT). Both are inverted in logs,
sigma = arcsinh(exp(z))/a with z = (ln rate + m ln(1 - 2 f_l) + Q/RT - ln A)/n, so exp(Q/RT) never overflows (at
1 K it would) and arcsinh(e^z) is taken as z + ln 2 where e^-2z is below round-off. Law "flight": the plastic branch
where f_l = 0 (decision 1: below 623 K it is Chen's branch extended, the labelled `provisional` placeholder of open
input 1, `FLOW_STRESS_PROVISIONAL`; from 623 K to the solidus Chen's own range), the semi-solid branch for
0 < f_l < 0.5, zero at and above 50 % liquid (the stiff replay carries no viscous rows; M4 adds them); the rate is
floored at `rate_min` (D6: 1e-6 1/s, bracket 1e-7 and 1e-5), since the law gives zero stress at zero rate. Law
"elastic": infinite. Law "constant:<Pa>": that value everywhere (checks 4-6). Measured (tests): Chen's own peak
stresses at 600 C, 8.2 MPa at 1e-3 1/s and 22.2 MPa at 1 1/s, are the semi-solid branch at f_l = 0.191, his measured
liquid fraction (his data reach f_l = 0.19); this table's Scheil curve has f_l = 0.233 there. The tearing law of §8.4
is M3's.

**Equation of state and shear modulus** (spec §8.2-8.3; M2 plan Task 3): `thermal_eos` (P = K(T)(rho/(s rho_free(T))
- 1) with its derivatives, sound speed and bulk modulus; u is the enthalpy from 300 K) and `thermal_shear_modulus`
(G(T), zero at 50 % liquid and above) and `thermal_sound_speed` (the longitudinal sqrt(c_f^2 + 4G/(3 rho)) the
strength provides for Spheral's time step): the numpy core of the runner's PyThermalEOS and PyThermalStrength.

**Mechanical columns** (decision 2: labelled M1 placeholders; decision 3: the liquid's free density is the finite
element's liquid density), built by `mechanical_columns`:
E = E_room to the solidus, then linear in T to zero at 50 % liquid (`T_half`), zero above; nu constant;
G = E / (2 (1 + nu)); K = E / (3 (1 - 2 nu)) to the solidus and held at its solidus value above (mush, slurry,
liquid); rho_free from the specific volumes of the expanded solid, rho_s / (1 + alpha (T - 293 K))^3, and of the
finite-element liquid, mixed by liquid fraction: 1/rho_free = (1 - f_l)/rho_solid(T) + f_l/rho_liquid, exactly
rho_liquid where f_l = 1. The 6.5 % volume change on melting is recorded in the input file but does not enter
rho_free (decision 3 replaces it); the change it implies at the liquidus is `implied_melting_volume_change`.
`MaterialTable.provisional` lists every column built from a provisional input.

Numpy and the standard library only (the runner imports this module under Spheral's Python)."""
from __future__ import annotations

import hashlib
import json
import math
import os

import numpy as np

SCHEMA = "spheral_frag.material_table"
SCHEMA_VERSION = 1
T_ZERO = 300.0                     # K, zero of the core's enthalpy scale (spec §8.1)
T_GRID = (250.0, 1500.0, 1.0)      # K, the 1 K table's range and spacing (plan Task 3, Step 1)
T_REF_EXPANSION = 293.0            # K, where rho_solid is the room-temperature density (the FE's T_REF)

FL_NONE, FL_INTERP, FL_LINEAR = 0, 1, 2

# Chen G. et al. (2016) Eq. 20 (M2 plan fact 10; Step 4 plan §4): rate = A sinh(a sigma)^n exp(-Q/RT), sigma in MPa
R_GAS = 8.314                                                # J/(mol K), as Chen's fit
CHEN_PLASTIC = dict(A=6.75872e11, a=0.014677, n=6.85385, Q=186163.0)
CHEN_SEMI_SOLID = dict(A=1.28926e8, a=0.03006, n=6.521135, Q=175131.0, m=15.81342)   # times (1 - 2 f_l)^m
T_CHEN_MIN = 623.0          # K (350 C), the lowest temperature of Chen's tests; the plastic branch below it is
                            # decision 1's provisional placeholder
FL_FLOW_ZERO = 0.5          # liquid fraction at and above which the flight law's flow stress is zero
RATE_MIN = 1.0e-6           # 1/s, the strain-rate floor (D6, decision 4); bracket RATE_MIN_BRACKET
RATE_MIN_BRACKET = (1.0e-7, 1.0e-5)
FLOW_LAWS = ("flight", "elastic", "constant:<Pa>")
FLOW_STRESS_PROVISIONAL = ("yield strength below 623 K: Chen 2016's plastic branch extended (decision 1 of the M2 "
                           "review, 2026-10-09); the handbook curve of 7075-T6 (open input 1) is owed before "
                           "production")

# Li Yageng et al. (2014), China Foundry 11, 79-84: semi-solid 7075 slurry viscosity (Step 4 plan §4; spec §2)
LI_A, LI_B, LI_N, LI_C = 0.871, 0.00849, 0.74924, 3.7311    # eta = [A - B g^N] exp(C f_s) Pa s
LI_FS_MIN, LI_FS_MAX = 0.1, 0.5                              # the fit's solid-fraction range
LI_SHEAR_MAX = 367.0                                         # 1/s, highest rate where the fit matches the data

# np.interp's last step, slope * (x - xp[j]) + fp[j], is contracted into one fused multiply-add by some numpy builds
# and not by others (measured 2026-10-07: drama_env's numpy 2.5.3, clang on macOS arm64, fuses; Spheral's numpy
# 1.26.4, gcc on Linux aarch64, does not), so the two differ by one unit in the last place at about 0.06 % of
# temperatures. The table records whether the numpy that built it (the one that wrote the frames) fuses
# (`interp_fma` in the header); `interp` reproduces that build's result on the other kind: the unfused form with plain
# numpy operations, the fused one in IEEE quad precision where numpy's longdouble is quad (Linux aarch64; the product
# of two doubles is exact there, so only a double rounding, of probability about 2^-60, could differ). Where neither
# applies (an x86-64 longdouble), np.interp is used and `MaterialTable.interp_exact` is False.
_FMA_PROBE = (float.fromhex("0x1.a278fee007243p+9"), (float.fromhex("0x1.a2p+9"), float.fromhex("0x1.a28p+9")),
              (float.fromhex("0x1.ab897af86fe0cp-4"), float.fromhex("0x1.b2fc3c31416fdp-4")),
              float.fromhex("0x1.b293e4de5d9bdp-4"))        # x, xp, fp, the fused result (unfused: ...9bc)
QUAD = np.finfo(np.longdouble).nmant >= 112


def _sinh_law_inverse(log_rate, T, A, a, n, Q):
    """Pa: sigma with rate = A sinh(a sigma[MPa])^n exp(-Q/RT), from ln(rate) (-inf gives 0), in logs (module
    docstring)."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        z = (log_rate + Q/(R_GAS*T) - math.log(A))/n
        big = z > 20.0                                       # e^-2z < 5e-18: arcsinh(e^z) = z + ln 2 to round-off
        asinh = np.where(big, z + math.log(2.0), np.arcsinh(np.exp(np.where(big, 0.0, z))))
        return 1.0e6*asinh/a                                 # MPa -> Pa (inf only as T -> 0)


def _log_rate(rate):
    rate = np.asarray(rate, dtype=float)
    with np.errstate(divide="ignore"):
        return np.log(np.maximum(rate, 0.0))


def chen_plastic(T, rate) -> np.ndarray:
    """Pa, Chen et al. (2016) Eq. 20's plastic branch at temperature T (K) and strain rate (1/s); zero at zero rate."""
    T, lr = np.broadcast_arrays(np.asarray(T, dtype=float), _log_rate(rate))
    c = CHEN_PLASTIC
    return _sinh_law_inverse(lr, T, c["A"], c["a"], c["n"], c["Q"])


def chen_semi_solid(T, f_l, rate) -> np.ndarray:
    """Pa, Chen et al. (2016) Eq. 20's semi-solid branch at T (K), liquid fraction f_l and strain rate (1/s); zero at
    f_l >= 0.5 (where (1 - 2 f_l)^m vanishes) and at zero rate."""
    T, f_l, rate = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(f_l, dtype=float),
                                       np.asarray(rate, dtype=float))
    c = CHEN_SEMI_SOLID
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = _log_rate(rate) + c["m"]*np.log(np.maximum(1.0 - 2.0*f_l, 0.0))
    s = _sinh_law_inverse(lr, T, c["A"], c["a"], c["n"], c["Q"])
    return np.where(f_l >= FL_FLOW_ZERO, 0.0, s)


def parse_flow_law(law):
    """("flight" | "elastic" | "constant", value in Pa or None); ValueError for anything else."""
    if law in ("flight", "elastic"):
        return law, None
    if isinstance(law, str) and law.startswith("constant:"):
        try:
            v = float(law.split(":", 1)[1])
        except ValueError:
            v = float("nan")
        if math.isfinite(v) and v >= 0.0:
            return "constant", v
    raise ValueError("flow-stress law must be one of {}, got {!r}".format(FLOW_LAWS, law))


def flow_stress_provisional(law) -> tuple:
    """The provisional labels a law carries (for prepare and run to print): the flight law's solid yield below 623 K."""
    return (FLOW_STRESS_PROVISIONAL,) if parse_flow_law(law)[0] == "flight" else ()


def flow_stress(table, T, rate, law="flight", rate_min=RATE_MIN) -> np.ndarray:
    """Pa, the flow (yield) stress of spec §8.3 at temperature T and strain rate (module docstring): law "flight"
    (Chen's plastic branch where f_l = 0, provisional below 623 K; the semi-solid branch for 0 < f_l < 0.5; zero
    above), "elastic" (inf) or "constant:<Pa>". The rate is floored at rate_min (D6) for the flight law."""
    kind, value = parse_flow_law(law)
    T = np.asarray(T, dtype=float)
    T, rate = np.broadcast_arrays(T, np.asarray(rate, dtype=float))
    if kind == "elastic":
        return np.full(T.shape, np.inf)
    if kind == "constant":
        return np.full(T.shape, value)
    r = np.maximum(rate, rate_min)
    f_l = table.liquid_fraction(T)
    plastic = chen_plastic(T, r)
    semi = chen_semi_solid(T, f_l, r)
    return np.where(f_l <= 0.0, plastic, np.where(f_l < FL_FLOW_ZERO, semi, 0.0))


def thermal_eos(table, rho, u, scale=1.0) -> dict:
    """The equation of state of spec §8.2 on whole arrays (numpy only; the runner's PyThermalEOS writes these values
    into Spheral's fields). u is the specific thermal energy, J/kg from 300 K: the table's enthalpy, latent heat
    included, so T = table.temperature(u) (the table's own zero is 300 K: enthalpy(300 K) = 0 exactly). With the mass
    scale factor s, eta = rho / (s rho_free(T)):

      P       = K(T) (eta - 1)
      dP/drho = K(T) / (s rho_free(T)),  c = sqrt(dP/drho),  bulk modulus rho dP/drho = K(T) eta
      dP/du   = (dP/dT)_rho / (dh/dT), dP/dT = K'(T)(eta - 1) - K(T) eta rho_free'(T)/rho_free(T): exact for the
                table's piecewise-linear K and rho_free and piecewise-quadratic h (one-sided at their nodes)
      gamma   = 1 + (s/rho) dP/du, the Grueneisen form 1 + (1/rho) dP/du with the physical density rho/s
      c_p     = table.cp(T) (no latent part)

    Returns a dict of arrays: T, P, dPdrho, dPdu, cs, K (the bulk modulus rho dP/drho), cp, gamma. P, s dPdrho,
    sqrt(s) cs, K, dPdu and gamma depend on (rho, s) only through rho/s: the stresses are invariant under mass
    scaling."""
    rho = np.asarray(rho, dtype=float)
    u = np.asarray(u, dtype=float)
    T = table.temperature(u)
    K = table.bulk_modulus(T)
    rf = table.rho_free(T)
    srf = scale*rf
    eta = rho/srf
    dPdrho = K/srf
    dPdT = table.grid_slope("K", T)*(eta - 1.0) - K*eta*table.grid_slope("rho_free", T)/rf
    dPdu = dPdT/table.dh_dT(T)
    return dict(T=T, P=K*(eta - 1.0), dPdrho=dPdrho, dPdu=dPdu, cs=np.sqrt(dPdrho), K=rho*dPdrho,
                cp=table.cp(T), gamma=1.0 + scale*dPdu/rho)


def thermal_shear_modulus(table, T) -> np.ndarray:
    """Pa, the shear modulus of spec §8.3 at T: the table's G(T), zero at and above 50 % liquid."""
    T = np.asarray(T, dtype=float)
    return np.where(table.liquid_fraction(T) >= FL_FLOW_ZERO, 0.0, table.shear_modulus(T))


def thermal_sound_speed(table, T, rho, fluid_sound_speed) -> np.ndarray:
    """m/s, the longitudinal sound speed the strength provides (as Spheral's SteinbergGuinan does):
    sqrt(c_f^2 + 4 G(T) / (3 rho)), c_f the equation of state's sound speed and rho the density Spheral holds (scaled
    by s, so that c scales as 1/sqrt(s) like c_f); G is `thermal_shear_modulus`, zero at 50 % liquid and above."""
    cf = np.asarray(fluid_sound_speed, dtype=float)
    return np.sqrt(cf*cf + (4.0/3.0)*thermal_shear_modulus(table, T)/np.asarray(rho, dtype=float))


def interp_fuses() -> bool:
    """Whether this numpy's np.interp computes its last step as a fused multiply-add."""
    x, xp, fp, fused = _FMA_PROBE
    return float(np.interp(x, xp, fp)) == fused


def _interp_emulated(x, xp, fp, fused):
    """np.interp's algorithm (increasing xp, numpy's end and node rules) with its last step fused or not."""
    x = np.asarray(x, dtype=float)
    n = len(xp)
    j = np.clip(np.searchsorted(xp, x, side="right") - 1, 0, n - 2)
    x0, f0 = xp[j], fp[j]
    slope = (fp[j + 1] - f0) / (xp[j + 1] - x0)
    d = x - x0
    if fused:
        res = (slope.astype(np.longdouble) * d.astype(np.longdouble) + f0.astype(np.longdouble)).astype(np.float64)
    else:
        res = slope * d + f0
    res = np.where(x == x0, f0, res)
    return np.where(x < xp[0], fp[0], np.where(x >= xp[-1], fp[-1], res))


def interp(x, xp, fp, fused=None):
    """np.interp(x, xp, fp) as a numpy whose last step is fused (True) or not (False) computes it; None: this one's."""
    if fused is None or fused == interp_fuses() or len(xp) < 2 or (fused and not QUAD):
        return np.interp(np.asarray(x, dtype=float), xp, fp)
    return _interp_emulated(x, np.asarray(xp, dtype=float), np.asarray(fp, dtype=float), fused)


# the node table: name -> what it is (all float64; written by fe.material_arrays)
NODE_ARRAYS = {
    "h_T": "K, enthalpy-table nodes (the FE's _T_h)",
    "h_nodes": "J/kg from the FE's T_ref_fe, enthalpy at h_T (the FE's _h_nodes)",
    "h_cp0": "J/(kg K), c_p at each interval's lower node, as the FE's enthalpy() computes it",
    "h_cp1": "J/(kg K), c_p at each interval's upper node",
    "h_slope": "J/(kg K), latent slope L df_l/dT inside each interval (the FE's _latent_slope; 0 without one)",
    "cp_end": "J/(kg K), the c_p table's first and last values (linear extrapolation beyond h_T)",
    "T_cp": "K, the c_p table's nodes",
    "cp_table": "J/(kg K), the c_p table",
    "fl_T": "K, liquid-fraction nodes (Scheil: the FE's _fl_T; linear range: [T_s, T_l]; none: empty)",
    "fl": "-, liquid fraction at fl_T",
}
SCALARS = {
    "fl_kind": "0 none, 1 Scheil (np.interp), 2 linear range",
    "T_ref_fe": "K, the FE's zero of enthalpy",
    "h_offset_300": "J/kg, h_FE(300 K), subtracted from every FE enthalpy",
    "T_solidus": "K, the loaded FE material's solidus (Scheil: the foot of the eutectic ramp, 748 K)",
    "T_liquidus": "K",
    "T_feed": "K, the FE's T_feed (the film feed ramp's top); NaN where the FE has none",
    "T_half": "K, 50 % liquid on the table's own f_l curve",
    "T_bridge": "K, f_l = 0.9 (f_s = 0.1) on the table's own curve: start of the viscosity bridge",
    "latent_heat": "J/kg",
    "rho_solid": "kg/m^3, the FE density",
    "rho_liquid": "kg/m^3, the FE liquid density (NaN without liquid properties)",
    "mu_liquid": "Pa s",
    "sigma_liquid": "N/m",
}
GRID_COLUMNS = {
    "T_grid": "K, 250-1,500 K every 1 K",
    "h300": "J/kg from 300 K, = enthalpy(T_grid) exactly",
    "f_l": "-, = liquid_fraction(T_grid) exactly",
    "rho_free": "kg/m^3, free density (spec §8.2)",
    "E": "Pa, Young's modulus",
    "nu": "-, Poisson's ratio",
    "K": "Pa, bulk modulus",
    "G": "Pa, shear modulus",
}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def header_path(path):
    """The JSON header that goes with a table file: `<stem>.json` next to `<stem>.npz`."""
    stem, _ = os.path.splitext(os.fspath(path))
    return stem + ".json"


def load_mechanical(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def mechanical_columns(T, f_l, T_solidus, T_half, rho_solid, rho_liquid, mechanical: dict):
    """The 1 K table's mechanical columns from the input file's entries (module docstring) and the provisional
    columns: (columns dict, tuple of column names built from an entry whose status is not "approved")."""
    e = mechanical["entries"]
    T, f_l = np.asarray(T, dtype=float), np.asarray(f_l, dtype=float)
    E0, nu0, alpha = float(e["E_room"]["value"]), float(e["nu"]["value"]), float(e["alpha"]["value"])
    if not np.isfinite(T_solidus):                       # never melts: E constant, expansion only
        E = np.full_like(T, E0)
    else:
        ramp = np.clip((T_half - T) / (T_half - T_solidus), 0.0, 1.0)
        E = np.where(T <= T_solidus, E0, E0 * ramp)
    nu = np.full_like(T, nu0)
    G = E / (2.0 * (1.0 + nu))
    K_s = E0 / (3.0 * (1.0 - 2.0 * nu0))                  # E is E0 up to the solidus
    K = np.where(T <= T_solidus, E / (3.0 * (1.0 - 2.0 * nu)), K_s)
    rho_sol = rho_solid / (1.0 + alpha * (T - T_REF_EXPANSION)) ** 3
    if np.isfinite(rho_liquid):
        with np.errstate(divide="ignore", invalid="ignore"):
            mixed = 1.0 / ((1.0 - f_l) / rho_sol + f_l / rho_liquid)
        rho_free = np.where(f_l >= 1.0, rho_liquid, np.where(f_l <= 0.0, rho_sol, mixed))
    else:
        rho_free = rho_sol
    cols = {"rho_free": rho_free, "E": E, "nu": nu, "K": K, "G": G}
    depends = {"rho_free": ("rho_solid", "alpha", "rho_liquid"), "E": ("E_room", "E_T"), "nu": ("nu",),
               "K": ("E_room", "nu", "K_T"), "G": ("E_room", "E_T", "nu")}
    prov = tuple(c for c, keys in depends.items() if any(e[k]["status"] != "approved" for k in keys))
    return cols, prov


def write_table(path, arrays: dict, header: dict):
    """Write `arrays` (NODE_ARRAYS, SCALARS as 0-d, GRID_COLUMNS) as a plain numeric npz and `header` (plus the
    npz's sha256) as `<stem>.json`. Returns the header written."""
    missing = sorted((set(NODE_ARRAYS) | set(SCALARS) | set(GRID_COLUMNS)) - set(arrays))
    if missing:
        raise ValueError("material table lacks " + ", ".join(missing))
    out = {k: np.asarray(arrays[k], dtype=np.int64 if k == "fl_kind" else np.float64)
           for k in (*NODE_ARRAYS, *SCALARS, *GRID_COLUMNS)}
    path = os.fspath(path)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    np.savez(path, **out)
    header = dict(header, schema=SCHEMA, schema_version=SCHEMA_VERSION, file=os.path.basename(path),
                  npz_sha256=_sha256(path))
    with open(header_path(path), "w", encoding="utf-8") as fh:
        json.dump(header, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return header


class MaterialTable:
    """The finite-element material, tabulated (module docstring). Attributes: the SCALARS (floats; `fl_kind` an
    int), the node arrays, the 1 K columns, `name` (the FE material), `provisional` (columns from provisional
    inputs) and `header` (the JSON header)."""

    def __init__(self, arrays: dict, header: dict | None = None):
        self.header = dict(header or {})
        for k in NODE_ARRAYS:
            setattr(self, k, np.asarray(arrays[k], dtype=np.float64))
        for k in SCALARS:
            setattr(self, k, int(arrays[k]) if k == "fl_kind" else float(arrays[k]))
        self.grid = {k: np.asarray(arrays[k], dtype=np.float64) for k in GRID_COLUMNS}
        self.T_grid = self.grid["T_grid"]
        self.name = self.header.get("fe_material", "")
        self.provisional = tuple(self.header.get("provisional", ()))
        self.interp_fma = bool(self.header.get("interp_fma", interp_fuses()))
        self.interp_exact = self.interp_fma == interp_fuses() or not self.interp_fma or QUAD
        if self.fl_kind not in (FL_NONE, FL_INTERP, FL_LINEAR):
            raise ValueError("unknown fl_kind {}".format(self.fl_kind))

    @classmethod
    def load(cls, path) -> "MaterialTable":
        """`material_table.npz` (or the directory holding it) with its JSON header; the header's sha256 must match."""
        path = os.fspath(path)
        if os.path.isdir(path):
            path = os.path.join(path, "material_table.npz")
        with open(header_path(path), encoding="utf-8") as fh:
            header = json.load(fh)
        if header.get("schema") != SCHEMA:
            raise ValueError("{} is not a {} header".format(header_path(path), SCHEMA))
        if header.get("npz_sha256") not in (None, _sha256(path)):
            raise ValueError("{} does not match the sha256 in its header".format(path))
        with np.load(path, allow_pickle=False) as z:
            arrays = {k: z[k] for k in z.files}
        return cls(arrays, header)

    @property
    def melts(self) -> bool:
        return self.fl_kind != FL_NONE

    # -- energy ----------------------------------------------------------------------------------------------

    def _enthalpy_fe(self, T):
        """The finite element's h(T) (zero at T_ref_fe), operation for operation."""
        T = np.asarray(T, dtype=float)
        hT = self.h_T
        i = np.clip(np.searchsorted(hT, T, side="right") - 1, 0, len(hT) - 2)
        T0, T1 = hT[i], hT[i + 1]
        cp0, cp1 = self.h_cp0[i], self.h_cp1[i]
        x = np.clip(T, T0, T1) - T0
        h = self.h_nodes[i] + (cp0 + self.h_slope[i]) * x + 0.5 * (cp1 - cp0) / (T1 - T0) * x * x
        h = np.where(T < hT[0], self.h_nodes[0] + self.cp_end[0] * (T - hT[0]), h)
        h = np.where(T > hT[-1], self.h_nodes[-1] + self.cp_end[1] * (T - hT[-1]), h)
        return h

    def enthalpy(self, T) -> np.ndarray:
        """J/kg from 300 K, latent heat included: h_FE(T) - h_FE(300 K)."""
        return self._enthalpy_fe(T) - self.h_offset_300

    def temperature(self, h) -> np.ndarray:
        """Exact inverse of `enthalpy` (the finite element's quadratic solve on h + h_FE(300 K))."""
        h = np.asarray(h, dtype=float) + self.h_offset_300
        hT, hn = self.h_T, self.h_nodes
        i = np.clip(np.searchsorted(hn, h, side="right") - 1, 0, len(hT) - 2)
        T0, T1 = hT[i], hT[i + 1]
        cp0, cp1 = self.h_cp0[i], self.h_cp1[i]
        b = cp0 + self.h_slope[i]
        a = (cp1 - cp0) / (T1 - T0)
        dh = np.clip(h, hn[i], hn[i + 1]) - hn[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(np.abs(a) > 1e-12, (np.sqrt(b * b + 2.0 * a * dh) - b) / a, dh / b)
        T = T0 + x
        T = np.where(h < hn[0], hT[0] + (h - hn[0]) / self.cp_end[0], T)
        T = np.where(h > hn[-1], hT[-1] + (h - hn[-1]) / self.cp_end[1], T)
        return T

    def cp(self, T) -> np.ndarray:
        """J/(kg K), the finite element's c_p(T) (np.interp on its table; no latent part)."""
        return interp(T, self.T_cp, self.cp_table, self.interp_fma)

    def liquid_fraction(self, T) -> np.ndarray:
        """f_l(T), as the finite element computes it (bitwise; module docstring, fl_kind)."""
        T = np.asarray(T, dtype=float)
        if self.fl_kind == FL_INTERP:
            return interp(T, self.fl_T, self.fl, self.interp_fma)
        if self.fl_kind == FL_LINEAR:
            return np.clip((T - self.fl_T[0]) / (self.fl_T[1] - self.fl_T[0]), 0.0, 1.0)
        return np.zeros_like(T)

    def h_full_liquid(self) -> float:
        """J/kg from 300 K: the enthalpy at the liquidus, fully liquid (the FE's enthalpy_liquid(T_liquidus) less
        h_FE(300 K)); inf for a material that never melts."""
        T = np.asarray(self.T_liquidus)
        h = self._enthalpy_fe(T) + (self.latent_heat * (1.0 - self.liquid_fraction(T)) if self.melts else 0.0)
        return float(h - self.h_offset_300)

    # -- mechanical columns (1 K table) ----------------------------------------------------------------------

    def _grid(self, key, T):
        return interp(T, self.T_grid, self.grid[key], self.interp_fma)

    def rho_free(self, T) -> np.ndarray:
        """kg/m^3, the free density (spec §8.2), np.interp on the 1 K table."""
        return self._grid("rho_free", T)

    def young(self, T) -> np.ndarray:
        return self._grid("E", T)

    def bulk_modulus(self, T) -> np.ndarray:
        return self._grid("K", T)

    def shear_modulus(self, T) -> np.ndarray:
        return self._grid("G", T)

    def poisson(self, T) -> np.ndarray:
        return self._grid("nu", T)

    def grid_slope(self, key, T) -> np.ndarray:
        """d(column)/dT of the 1 K column `key` as np.interp evaluates it: the slope of the interval holding T (the
        upper one at a node, as np.interp's search), zero outside the grid (np.interp holds the end values)."""
        T = np.asarray(T, dtype=float)
        xp, fp = self.T_grid, self.grid[key]
        j = np.clip(np.searchsorted(xp, T, side="right") - 1, 0, len(xp) - 2)
        slope = (fp[j + 1] - fp[j])/(xp[j + 1] - xp[j])
        return np.where((T < xp[0]) | (T >= xp[-1]), 0.0, slope)

    def dh_dT(self, T) -> np.ndarray:
        """J/(kg K), the derivative of `enthalpy`, latent heat included (c_p plus L df_l/dT), exactly as the
        piecewise quadratic gives it (the upper interval at a node, cp_end outside the nodes)."""
        T = np.asarray(T, dtype=float)
        hT = self.h_T
        i = np.clip(np.searchsorted(hT, T, side="right") - 1, 0, len(hT) - 2)
        T0, T1 = hT[i], hT[i + 1]
        cp0, cp1 = self.h_cp0[i], self.h_cp1[i]
        x = np.clip(T, T0, T1) - T0
        d = cp0 + self.h_slope[i] + (cp1 - cp0)/(T1 - T0)*x
        d = np.where(T < hT[0], self.cp_end[0], d)
        return np.where(T > hT[-1], self.cp_end[1], d)

    def implied_melting_volume_change(self) -> float:
        """rho_free just below the liquidus' solid branch over the liquid density, minus 1: the volume change on
        melting the table implies at the liquidus (decision 3), against the input file's placeholder."""
        alpha = float(self.header.get("mechanical", {}).get("entries", {}).get("alpha", {}).get("value", np.nan))
        rho_sol = self.rho_solid / (1.0 + alpha * (self.T_liquidus - T_REF_EXPANSION)) ** 3
        return float(rho_sol / self.rho_liquid - 1.0)

    # -- viscosity (spec §8.3 slurry and liquid rows) --------------------------------------------------------

    @staticmethod
    def li_viscosity(f_s, shear_rate=0.0) -> np.ndarray:
        """Li et al. (2014), Pa s; f_s clipped to the fit's 0.1-0.5, the shear rate to 0-367 1/s."""
        g = np.clip(np.asarray(shear_rate, dtype=float), 0.0, LI_SHEAR_MAX)
        f_s = np.clip(np.asarray(f_s, dtype=float), LI_FS_MIN, LI_FS_MAX)
        return (LI_A - LI_B * g ** LI_N) * np.exp(LI_C * f_s)

    def viscosity(self, T, shear_rate=0.0) -> np.ndarray:
        """Pa s: inf below 50 % liquid; Li et al. at f_s = 1 - f_l(T) up to `T_bridge`; log-linear in T to the
        liquid's viscosity at the liquidus; the liquid's above (module docstring)."""
        T = np.asarray(T, dtype=float)
        T, g = np.broadcast_arrays(T, np.asarray(shear_rate, dtype=float))
        f_l = self.liquid_fraction(T)
        slurry = self.li_viscosity(1.0 - f_l, g)
        eta0 = self.li_viscosity(LI_FS_MIN, g)
        with np.errstate(all="ignore"):              # the bridge is evaluated everywhere, used only inside it
            s = (T - self.T_bridge) / (self.T_liquidus - self.T_bridge)
            bridge = np.exp(np.log(eta0) + s * (math.log(self.mu_liquid) - np.log(eta0))) \
                if self.mu_liquid > 0.0 else np.full(T.shape, np.nan)
        eta = np.where(T >= self.T_liquidus, self.mu_liquid, np.where(T > self.T_bridge, bridge, slurry))
        return np.where(f_l < 0.5, np.inf, eta)
