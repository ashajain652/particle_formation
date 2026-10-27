# Wake Column Chemistry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `reentry_model/wake_column/`, which follows the aluminium a demising AA7075 droplet vaporises: it slows
the hot vapour to ambient (Stage A), then runs a 12-hour box model of a diffusing horizontal cylinder with Plane et al.
(2021) gas chemistry, sectional alumina/hydroxide particles tallied by source molecule and route, uptake onto
background meteoric smoke, coagulation and settling, for every release altitude by day and by night.

**Architecture:** Data-driven (species, reactions, ambient columns and Megner profiles are committed tables under
`data/wake_column/`). Small modules each own one piece of physics. `box.py` assembles them into one stiff ODE system
of extensive counts (molecules and atoms in the whole plume), so aluminium conservation is a plain sum. `run.py`
chains droplet run, then trail adapter, then Stage A, then box, then output files.

**Tech Stack:** Python 3.12 in `drama_env` (numpy, scipy `solve_ivp` BDF, scipy.io.netcdf_file, pymsis 0.13,
matplotlib); pytest.

**Spec:** `Reentry Model/docs/superpowers/specs/2026-10-27-wake-column-chemistry-design.md` (approved by Asha
2026-10-27). Read it before any task; section numbers below (§4.2 etc.) refer to it.

## Global Constraints

- Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; every command runs from inside `Reentry Model/`.
- Package `reentry_model/wake_column/`; tests in `reentry_model/wake_column/tests/`, added to `testpaths` in `pytest.ini`.
- Analysis scripts in `reentry_model/wake_column/analysis/`, run with `"$PY" -m reentry_model.wake_column.analysis.<name>`, writing only to `output/wake_column_output/` (git-ignored).
- Committed inputs in `data/wake_column/`; every table carries a `source` column or header comment naming its origin, and marks values that are assumed rather than sourced.
- Run names encode the whole configuration (diameter, release altitude, day/night, ambient source, every non-default setting); a sweep skips a run whose summary JSON exists.
- No unmeasured rate or efficiency is hard-coded: each is a named field of `Settings` or a row of a committed table (spec anti-goal).
- Aluminium is conserved to integration tolerance: gas + particles + smoke + settled = vaporised (relative 1e-6).
- The condensation pathways are never lumped: one labelled pathway per pair of aluminium-bearing species.
- `reentry_model/droplet_reentry/` is not modified.
- SI units inside the code (m, m^-3, m^3/s, m^6/s); literature units appear only in the committed tables, converted on load.
- Commit messages: scope prefix and a sentence describing the physics or behaviour, e.g. `wake_column: ...`; end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Work on branch `wake-column`.

## Decisions this plan makes that the spec left to the builder

These are recorded here once so every task uses the same choice. Task 10 copies them into `docs/model_assumptions.md`.

1. **Extensive state.** The ODE state is total counts in the plume, not concentrations. Ambient species gain `c_ambient · dV/dt` (entrainment); aluminium-bearing species have no ambient source, so their concentration falls only through dilution. This is the spec's §4.2 entrainment equation written for totals.
2. **Box volume** (spec §4.2 "the builder derives the expression"). The box concentration is the centreline value of a uniform line segment of length L broadened by a Gaussian with r² = r0² + 4∫D dt:
   `V(r, L) = π r² L / erf(L / (2 r))`. It tends to π r² L when r ≪ L and to π^1.5 r³ when r ≫ L. (A mean-square volume would be twice as large for second-order rates. The spec fixed V = π r² L, so the centreline form is used and this is an assumption entry.)
3. **Eddy switch** `S(r)`: 0 below `l_inner/√2`, 1 above `√2·l_inner`, cubic smoothstep in ln r between.
4. **Molecular diffusion.** Hard-sphere (first-order) Chapman–Enskog, `D = 3/(8 n σ12²) · sqrt(k_B T / (2π μ))`, with air as the pseudo-species (σ = 3.711 Å, 28.96 amu). Lennard-Jones well depths for aluminium molecules are unknown, so the collision integral is 1. Stage A uses the same diameters.
5. **Collision diameters.** Svehla (1962, NASA TR R-132): N2 3.798, O2 3.467, O 3.050, H 2.708, H2 2.827, H2O 2.641, CO 3.690, CO2 3.941, Al 2.655 Å. Every one-aluminium molecule is 3.5 Å, and O3 is 3.6 Å (both assumed, flagged).
6. **Dipole moments.** AlO 4.45 D (Bai & Steimle 2020, as used by Plane 2021); AlOH 0.97 D (Sikorska & Skurski 2009, as used by Plane 2021). Others are unknown, and pairs involving them use the hard-sphere collision rate (flagged).
7. **Formation rate** for pair (A, B) = `max(k_DD, k_HS) × spin × efficiency × formation_scale`:
   - `k_DD = π Γ(2/3) v̄ (2 C6 / k_B T)^(1/3)`, with the Keesom-averaged `C6 = (2/3)(μ_A μ_B / 4πε0)² / (k_B T)` and `v̄ = sqrt(8 k_B T / π μ)`.
   - `k_HS = π σ_AB² v̄`.
   - For AlO + AlO at 298 K with 4.36 D this gives 1.52e-9 cm³/s. Gómez Martín et al. (2017) quote 1.1e-9; the difference goes in the assumption entry.
   - `spin = g_P / (g_A g_B)`, where g_P = 1 when the pair's total electron count is even and 2 when it is odd.
8. **Bins by aluminium atoms.** Pivots `n_k = 2 · q^k` (`q = Settings.bin_ratio`, default 2^(1/2) = 1.41421: Task 8's coagulation test gave a self-preserving geometric standard deviation of 1.506 at q = 2, outside 1.45 ± 0.03 by numerical diffusion, and 1.476 at 2^(1/2) and 1.468 at 2^(1/3); 60 bins at the default). They run from bin 0 (two Al atoms, a 0.476 nm alumina sphere) until the alumina diameter passes `d_max_nm` (default 400). Number and Al atoms are conserved by fixed-pivot splitting (Kumar & Ramkrishna 1996):
   - An entity of n atoms with `n_k ≤ n < n_{k+1}` puts the fraction `ζ = (n_{k+1} − n)/(n_{k+1} − n_k)` of its number in bin k and `1 − ζ` in bin k+1.
   - Its tallies split in proportion to the atoms each bin receives.
   - Above the last pivot, the entity goes to the last bin with its number scaled to conserve atoms.
9. **Atoms delivered per molecule** (on growth, formation and smoke uptake), as (Al, O, H):
   - Al (1,0,0), AlO (1,1,0), OAlO (1,2,0), OAlO2 (1,3,0), AlOH (1,1,1), AlOH2 (1,2,2).
   - AlCO3 (1,1,0): CO2 is assumed to leave on sticking.
   - The O atom released by AlO + AlO → Al2O + O is not tracked.
   - Particle oxidation does not consume gas-phase O, O2 or O3. The aluminium is a trace.
10. **Particle size and density.** Volume per Al atom = `(1 − f) v_Al2O3 + f v_Al(OH)3`, with `f = min(H/Al / 3, 1)`, `v_Al2O3 = 2.8218e-29 m³` (3.0 g/cm³) and `v_Al(OH)3 = 5.3524e-29 m³` (2.42 g/cm³).
11. **Oxygen uptake.** Rate per bin = Σ over O, O2 and O3 of collision rate × `oxidation_eff` × `g(O/Al)`, adding one O atom per event. `g` is 1 below O/Al = 1.485, 0 at 1.5 and above, and a cubic smoothstep (3s² − 2s³ in s = (1.5 − O/Al)/0.015) between; Task 8 replaced the linear ramp originally written here because its slope jump from −67 to 0 at 1.5 stalled BDF's Newton iteration.
12. **Three-body steps.** [M] is the total air number density, using Plane's N2-measured coefficients for all air (O2 assumed as efficient as N2). The spec's test 4 therefore uses `1/(k6 n_total [O2])`.
13. **Duplicate reaction removed.** The spec lists "AlO + O2 → AlO2 + O and its reverse". The reverse, OAlO + O → AlO + O2, is Plane's R3, so only R3 is kept and the reaction table notes this.
14. **Reaction validity ranges** where the source states none (column `validity_source = assumed`):
    - Plane's (T/300)^(1/6) entries: 100–1000 K.
    - R1: 23–295 K (Le Picard 1997 CRESU and Gómez Martín 2017 measurements).
    - R6 and R7: 100–600 K.
    - Gobrecht/Mangan additions: 100–2500 K.
    - Stage A calibrates each hot reaction at the top of its range.
15. **Stage A emission.** Summed over a uniformly evaporating sphere, cosine-law emission is isotropic, so each atom gets an isotropic thermal velocity with speed drawn from the flux-weighted Maxwellian (pdf ∝ v³ exp(−m v²/2 k_B T_d)), plus V along +x.
    - Only the droplet's speed matters, because the air is still and isotropic, so the adapter passes speeds, not vectors.
    - Emission samples sit at the midpoints of `n_emission_samples` equal-aluminium-mass intervals of the loss history, each with weight 1/n.
16. **Trail length.** L is the path flown from release to the moment 90 % of the total aluminium loss has occurred (the spec's table: 1 µm at 80 km gives 51 m; 10 µm at 60 km gives 40 m). The altitude is the Al-mass-weighted mean over the whole loss.
17. **"Condensed"** means everything that left the gas phase: plume particles + smoke uptake + settled. t50 and t90 are the first output times this reaches 0.5 and 0.9 of the inventory, or null.
18. **Steady state.** The steady time is the earliest output time t_ss such that, for every t ≥ t_ss with 10 t ≤ t_end, every aluminium reservoir (each gas species, plume particles, smoke, settled) changes by less than 0.01 of the inventory between t and 10 t. At least one decade of evidence is required (t_ss ≤ t_end/10); it is null otherwise.
19. **Solver scaling.** Aluminium-bearing counts are divided by the inventory N_Al (or by 1 when N_Al is 0); each ambient count by c_ambient · V0; r² by r0². BDF uses rtol = `Settings.rtol` (1e-6) and atol = 1e-12 on the scaled variables.

## Review Focus

1. **No aluminium vaporised.** A droplet case that vaporises no aluminium (inventory 0) should produce a summary JSON with `skipped: "no aluminium vaporised"`, not a division by zero. The test is in Task 9.
2. **Altitude outside the tables.** A trail whose mean altitude falls outside the ambient or Megner tables (e.g. a 50 km release that sinks below the table bottom) should raise a `ValueError` naming the altitude and the table, never a silent extrapolation. The test is in Task 2.
3. **Extreme length-to-radius ratios.** Very short trails or very wide plumes (L/r from 1e-6 to 1e6) should give a finite, smooth `volume()` and `dlnV_dr2()` with the right limits. The test is in Task 3.
4. **Degenerate Stage A input.** A sample with droplet speed 0, or a column with zero atomic oxygen (night, low altitude), should give finite fractions and r0, with no NaN probabilities or zero-division in the free-flight sampling. The test is in Task 6.
5. **Rates outside their validity range.** A rate evaluated outside its range (WACCM can reach below 150 K, though not in this column) should be clipped and warned about once per reaction per run, not once per right-hand-side call. The test is in Task 1.

---

### Task 1: Package skeleton, species and reaction tables, settings

**Files:**
- Create: `reentry_model/wake_column/__init__.py`, `species.py`, `reactions.py`, `settings.py`
- Create: `data/wake_column/species.csv`, `data/wake_column/reactions.csv`
- Create: `reentry_model/wake_column/tests/__init__.py` (empty), `tests/conftest.py`, `tests/test_wake_tables.py`
- Modify: `pytest.ini` (append `reentry_model/wake_column/tests` to `testpaths`)

**Interfaces:**
- Produces:
  - `wake_column.DATA` (= `data/wake_column`) and `wake_column.OUTPUT` (= `output/wake_column_output`), the same pattern as `droplet_reentry.OUTPUT`.
  - `species.Species(name, mass_amu, n_Al, n_O, n_H, deliver_O, deliver_H, dipole_D, multiplicity, electrons, sigma_A, ambient: bool, source)`, frozen.
  - `species.load_species(path=None) -> dict[str, Species]`.
  - `species.AL_SPECIES = ("Al","AlO","OAlO","OAlO2","AlCO3","AlOH","AlOH2")` and `species.AMBIENT_SPECIES = ("N2","O2","O","O3","H","H2","H2O","CO","CO2")`. `AlOH2` is Al(OH)2; `species.LABEL` maps names to display labels.
  - `reactions.Reaction(rid, reactants: tuple, products: tuple, form, order, A, n, Ea_K, c0, c1, c2, J, E0_kJmol, P0, T_min, T_max, stage, day_only, validity_source, source)`, frozen. It has `k(T, warn: RateWarnings | None = None) -> float` in SI and clips T to [T_min, T_max].
  - `reactions.RateWarnings`: collects `(rid, T)`, keeps one message per rid, and has `.messages -> list[str]`.
  - `reactions.load_reactions(path=None, stage=("box","both")) -> list[Reaction]`.
  - `settings.Settings` (frozen dataclass, fields below) with `.tag() -> str`.
  - `settings.Case(d0_um: float, h0_km: float, time_of_day: str)`.
  - `settings.run_name(case, settings) -> str`.

`species.csv` columns follow the `Species` fields. Values (mass in amu; dipole in debye, blank = unknown; σ in Å):
- **Aluminium species**:

  | name | mass | dipole | multiplicity | σ |
  |---|---|---|---|---|
  | Al | 26.982 | 0 | 2 | 2.655 |
  | AlO | 42.981 | 4.45 | 2 | 3.5 |
  | OAlO | 58.980 | unknown | 2 | 3.5 |
  | OAlO2 | 74.979 | unknown | 2 | 3.5 |
  | AlCO3 | 86.990 | unknown | 2 | 3.5 |
  | AlOH | 43.989 | 0.97 | 1 | 3.5 |
  | AlOH2 | 60.996 | unknown | 2 | 3.5 |

- **Ambient species**: N2 28.014, O2 31.998, O 15.999, O3 47.998, H 1.008, H2 2.016, H2O 18.015, CO 28.010, CO2 44.009 amu, with the Svehla diameters from decision 5.
- `deliver_O` and `deliver_H` follow decision 9.
- `electrons` is the total electron count.
- The `source` column cites decision 5 or 6, or says "assumed".

`reactions.csv` columns: `rid, reactants, products, form, order, units, A, n, Ea_K, c0, c1, c2, J, E0_kJmol, P0, T_min, T_max, stage, day_only, validity_source, source, note`. Rows:
- **R1–R17** exactly as in Plane et al. (2021) Table 1, all `stage=box`:
  - R1: Al+O2→AlO+O, 1.7e-10 (T/300)^(1/6). `stage=both`, 23–295 K.
  - R2: AlO+O3→OAlO+O2, 1.3e-10.
  - R3: OAlO+O→AlO+O2, 1.9e-10. `stage=both`.
  - R4: OAlO+CO→AlO+CO2, 2.6e-11.
  - R5: AlO+CO→Al+CO2, 2.0e-12.
  - R6: AlO+O2+M→OAlO2+M, `log_quadratic` −35.137, 6.1052, −1.4089 (cm6).
  - R7: AlO+CO2+M→AlCO3+M, −38.736, 8.7342, −2.0202 (cm6).
  - R8: OAlO2+O→OAlO+O2, 1.2e-10.
  - R9: OAlO+H→AlOH+O, 2e-10.
  - R10: OAlO2+H→AlOH+O2, 2e-10.
  - R11: OAlO2+H2O→AlOH2+O2, 2e-10.
  - R12: AlOH2+H→AlOH+H2O, 5e-11.
  - R13: AlCO3+O→OAlO+CO2, 4e-11.
  - R14: AlCO3+O2→OAlO2+CO2, 6.6e-11.
  - R15: AlCO3+H→AlOH+CO2, 5e-11.
  - R16: AlCO3+H2O→AlOH2+CO2, 2e-10.
  - R17: AlOH→AlO+H, `photolysis` J = 3.3e-3, `day_only=1`.
  - All the bimolecular rates above have n = 1/6.
- **G42**: AlO+O2→OAlO+O, 3.54e-11 (T/300)^0.19 exp(−13140/T), `stage=both`. Note: "reverse is R3".
- **G33**: AlO+O→Al+O2, 7.15e-12 exp(−945/T), `stage=both`.
- **G26**: Al+O+M→AlO+M, 3.34e-32 (T/300)^−1.79 (cm6).
- **G38**: AlO+O+M→OAlO+M, 6.74e-30 (T/300)^−2.53 (cm6).
- **GH2**: AlO+H2→AlOH+H, 5.37e-13 (T/300)^2.77 exp(−2190/T).
- The Gobrecht entries cite Gobrecht et al. (2022) Table B.1 with the reaction number as given in the research notes; the implementer confirms each number against the PDF in `Literature Review/Alumina Nanoparticles/`.
- **D1**: AlO+M→Al+O+M, `form=threshold`, E0 = 498 kJ/mol (Gobrecht CBS-QB3 AlO bond energy), P0 = 1, `stage=hot`. Note: "upper bound, my call".
- Validity ranges follow decision 14.

`Settings` fields and defaults:
- **Hot start:** `hot_start="slowing_down"`, `hot_reaction_scale=1.0`, `thermal_cutoff=2.0`, `n_test_particles=10000`, `n_emission_samples=20`, `seed=20261027`.
- **Spreading:** `l_inner_m=30.0`, `eddy_scale=1.0`.
- **Efficiencies:** `formation_scale=1.0`, `formation_eff=()` (tuple of (pathway id, value)), `sticking_scale=1.0`, `sticking=()`, `oxidation_eff=1.0`, `smoke_sticking=1.0`.
- **Ambient and bins:** `ambient="waccm"`, `bin_ratio=2.0 ** 0.5` (2.0 until Task 8; see decision 8), `d_max_nm=400.0`.
- **Time grid and solver:** `t_end_s=43200.0`, `t_first_s=1e-6`, `points_per_decade=10`, `rtol=1e-6`.
- **Test-only switches:** `chemistry=True`, `particles=True`, `coagulation=True`, `settling=True`, `freeze_radius=False`, `reactions=None` (a tuple of rids, or None for all), `smoke=True`.

`tag()` appends `_<short><value:g>` for each non-default physics field, in a fixed order. The shorts are:
- form = `formation_scale`
- stick = `sticking_scale`
- ox = `oxidation_eff`
- hot = `hot_reaction_scale`
- hotdroplet = `hot_start="droplet_temperature"` (no value)
- eddy = `eddy_scale`
- lin = `l_inner_m`
- cut = `thermal_cutoff`
- smk = `smoke_sticking`
- q = `bin_ratio`
- fe-<pid> and st-<species> for the per-pathway and per-species entries

`run_name` = `wake_d{d0:06.2f}um_h{h0:05.1f}km_{time_of_day}_{ambient}{tag}`.

- [ ] **Step 1: Write the failing tests** in `tests/test_wake_tables.py`

```python
import math
import warnings
import pytest
from reentry_model.wake_column import species as sp, reactions as rx
from reentry_model.wake_column.settings import Settings, Case, run_name

CM3, CM6 = 1e-6, 1e-12

def test_species_table_has_every_network_species():
    s = sp.load_species()
    assert set(sp.AL_SPECIES) | set(sp.AMBIENT_SPECIES) <= set(s)
    assert s["AlO"].dipole_D == pytest.approx(4.45) and math.isnan(s["OAlO"].dipole_D)
    assert all(s[n].n_Al == 1 for n in sp.AL_SPECIES)

@pytest.mark.parametrize("T", [200.0, 300.0])
def test_plane_table_values_reproduced(T):            # spec test 8
    r = {x.rid: x for x in rx.load_reactions()}
    assert r["R1"].k(min(T, 295.0)) == pytest.approx(1.7e-10 * (min(T, 295.0) / 300) ** (1 / 6) * CM3, rel=1e-9)
    assert r["R11"].k(T) == pytest.approx(2e-10 * (T / 300) ** (1 / 6) * CM3, rel=1e-9)
    lg = math.log10(T)
    assert r["R6"].k(T) == pytest.approx(10 ** (-35.137 + 6.1052 * lg - 1.4089 * lg * lg) * CM6, rel=1e-9)
    assert r["R17"].k(T) == pytest.approx(3.3e-3) and r["R17"].day_only

def test_reverse_of_G42_is_not_duplicated():
    rids = [x.rid for x in rx.load_reactions()]
    assert "R3" in rids and "G42" in rids and len(rids) == len(set(rids))
    assert not any(x.reactants == ("OAlO", "O") for x in rx.load_reactions() if x.rid != "R3")

def test_hot_only_dissociation_not_in_box_network():
    assert "D1" not in [x.rid for x in rx.load_reactions()]
    assert "D1" in [x.rid for x in rx.load_reactions(stage=("hot", "both"))]

def test_clipping_warns_once_per_reaction():          # Review Focus 5
    r = {x.rid: x for x in rx.load_reactions()}["R6"]
    w = rx.RateWarnings()
    for _ in range(1000):
        assert r.k(50.0, w) == r.k(r.T_min)
    assert len(w.messages) == 1

def test_run_name_encodes_non_default_settings():
    c = Case(14.1, 76.0, "day")
    assert run_name(c, Settings()) == "wake_d014.10um_h076.0km_day_waccm"
    assert run_name(c, Settings(formation_scale=0.1, hot_start="droplet_temperature")).endswith("_form0.1_hotdroplet")
    assert run_name(c, Settings(ambient="msis_afgl")) != run_name(c, Settings())
```

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_tables.py -q`. Expected: errors with `ModuleNotFoundError: reentry_model.wake_column`.
- [ ] **Step 3: Write the two CSV tables and the four modules** to the interfaces above. The conftest mirrors `droplet_reentry/tests/conftest.py`: the project root on `sys.path`, plus a session fixture `parent` and `fixtures_dir`.
- [ ] **Step 4: Run the tests to see them pass.** Expected: 7 passed.
- [ ] **Step 5: Commit.** `git add` the files above, then commit with `wake_column: species and reaction tables (Plane 2021 R1-R17 plus Gobrecht/Mangan oxidation steps and the hot-only AlO dissociation) as data, rate evaluation with validity clipping, and the settings and run-name scheme`.

---

### Task 2: Ambient columns (WACCM baseline, NRLMSIS + AFGL cross-check)

**Files:**
- Create: `reentry_model/wake_column/ambient.py`, `analysis/__init__.py`, `analysis/extract_waccm_columns.py`, `analysis/build_msis_afgl_columns.py`
- Create (generated, committed): `data/wake_column/ambient_waccm_40.7N.csv`, `data/wake_column/ambient_msis_afgl_40.7N.csv`, `data/wake_column/afgl_midlatitude_winter.csv`
- Test: `reentry_model/wake_column/tests/test_wake_ambient.py`

**Interfaces:**
- Produces:
  - `ambient.AmbientState(z_km, T_K, p_Pa, n_total_m3, n: dict[str, float])`, frozen; `n` holds every name in `AMBIENT_SPECIES`, in m^-3.
  - `ambient.AmbientColumn.load(source: str, time_of_day: str) -> AmbientColumn`, with source "waccm" or "msis_afgl" and time_of_day "day" or "night".
  - `.at(z_km) -> AmbientState`: linear interpolation of T, log-linear of every density. Outside the table it raises `ValueError(f"altitude {z} km outside {path} ({zmin}-{zmax} km)")`.
- CSV schema (both files): `time_of_day, z_km, p_Pa, T_K, n_total_m3, N2, O2, O, O3, H, H2, H2O, CO, CO2`, levels 45–100 km, with a header comment naming the source file, the column longitude and the date.

The extraction script takes the WACCM path as an argument. The file is outside the repository:
`../../WACCM-Cheyenne/carma-cesm2-inputdata/FWmaCARMAnudged.f19_f19_mg17.carma_trop_strat.2000.001.cam.i.2001-01-01-00000_c220325.nc`, relative to the repository root.

The script generalises the throwaway script whose method is fixed here:
- Latitude nearest 40.7°; day at lon 180, night at lon 0.
- `p = hyam·P0 + hybm·PS`; hypsometric altitude from the interfaces with Rd = 287.05 and g = 9.80665.
- Number density = `q · 28.96/M · p/(k_B T)`.
- The extra variables are H2 and CO; N2 = n_total − Σ(other species).

`build_msis_afgl_columns.py` works as follows:
- `pymsis.calculate(np.datetime64("2001-01-01T00:00"), lon, 40.7, z_km, f107s=[150], f107as=[150], aps=[[4]*7])`, with lon 180 for day and 0 for night (my call: moderate, fixed solar indices; no download of space-weather files).
- The output index order is total mass density, N2, O2, O, He, H, Ar, N, anomalous O, NO, T (m^-3). The implementer checks this against pymsis 0.13's `Variable` enum.
- n_total is the sum of the number densities.
- O3, H2O and CO2 come from the AFGL midlatitude-winter volume mixing ratios times n_total. H2 and CO are copied from the WACCM column, because AFGL lacks them (assumed, noted in the header).
- The AFGL table (Anderson et al. 1986) is committed as `afgl_midlatitude_winter.csv`, taken from libRadtran's `data/atmmod/afglmw.dat` or the `joseki` package. If neither can be obtained, stop and ask Asha rather than inventing values.

- [ ] **Step 1: Run both scripts and commit nothing yet.** Run `"$PY" -m reentry_model.wake_column.analysis.extract_waccm_columns <waccm path>`. Expected: it writes `ambient_waccm_40.7N.csv` with day and night blocks. At about 71 km, O is near 3.7e15 m^-3 by day and 9.7e11 at night, and H2O is near 5.5e15 by day (the values found during design, 2026-10-27).
- [ ] **Step 2: Write the failing tests**

```python
import pytest
from reentry_model.wake_column.ambient import AmbientColumn

@pytest.mark.parametrize("src", ["waccm", "msis_afgl"])
def test_columns_cover_45_to_100_km_and_close(src):
    for tod in ("day", "night"):
        col = AmbientColumn.load(src, tod)
        for z in (45.0, 60.0, 76.0, 90.0, 100.0):
            a = col.at(z)
            assert 130 < a.T_K < 320
            assert sum(a.n.values()) == pytest.approx(a.n_total_m3, rel=0.02)
            assert all(v >= 0 for v in a.n.values())

def test_day_night_atomic_oxygen_contrast_below_80_km():
    day, night = AmbientColumn.load("waccm", "day").at(70.0), AmbientColumn.load("waccm", "night").at(70.0)
    assert day.n["O"] / night.n["O"] > 100
    assert day.n["H2O"] == pytest.approx(night.n["H2O"], rel=0.3)

def test_waccm_and_msis_agree_on_total_density_within_factor_2():
    for z in (60.0, 76.0, 90.0):
        a, b = AmbientColumn.load("waccm", "night").at(z), AmbientColumn.load("msis_afgl", "night").at(z)
        assert 0.5 < a.n_total_m3 / b.n_total_m3 < 2.0

def test_outside_table_raises_with_altitude_in_message():   # Review Focus 2
    with pytest.raises(ValueError, match="44.0 km"):
        AmbientColumn.load("waccm", "day").at(44.0)
```

- [ ] **Step 3: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_ambient.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 4: Implement `ambient.py`.** Run the msis/afgl script, then rerun the tests. Expected: 6 passed.
- [ ] **Step 5: Commit** the scripts, `ambient.py`, the tests and the three CSVs with `wake_column: ambient columns at 40.7 N, WACCM noon and midnight from the full-chemistry starting file, and NRLMSIS 2.1 + AFGL midlatitude winter as the cross-check`.

---

### Task 3: Megner profiles, diffusion and cylinder geometry

**Files:**
- Create: `reentry_model/wake_column/transport.py`
- Create (digitised, committed): `data/wake_column/megner2006_eddy_diffusion.csv` (`z_km, K_m2s`), `data/wake_column/megner2006_smoke.csv` (`z_km, r_nm, dN_dlnr_m3`)
- Test: `reentry_model/wake_column/tests/test_wake_transport.py`

**Interfaces:**
- Produces:
  - `transport.EddyProfile.load() -> EddyProfile`, with `.K(z_km) -> m²/s`: log-linear interpolation, `ValueError` outside the table.
  - `transport.eddy_switch(r_m, l_inner_m) -> float` (decision 3).
  - `transport.D_molecular(sp: Species, amb: AmbientState) -> m²/s` (decision 4).
  - `transport.D_total(r_m, sp, amb, eddy_K, settings) -> m²/s`, which is `D_mol + eddy_scale · K · S(r)`.
  - `transport.volume(r_m, L_m) -> m³` and `transport.dlnV_dr2(r_m, L_m) -> 1/m²` (decision 2).
  - `transport.settling_rate(v_s, r_m) -> 1/s` = `2 v_s/(π r)`.

Digitising:
- Render page 4 of `Literature Review/Alumina Nanoparticles/Megner_2006_Meteoric_Smoke_Distribution_Sensitivity.pdf` at 200 dpi with `/Users/ashajain/miniforge3/bin/python3` (it has PyMuPDF; `drama_env` does not), and read the figures.
- **Eddy diffusion:** the Fig. 3c reference (solid) profile, every 2.5 km over the plotted range.
- **Smoke:** Fig. 2's number and size distribution, at every altitude the figure resolves, as dN/dln r on its radius grid.
- Each CSV's header states "digitised by eye from Megner et al. (2006) ACP 6, 4415, Fig. 3c / Fig. 2, page 4; reading uncertainty about ±10 % in log value", and the date.
- If Fig. 2 does not reach 45 km, the smoke table stops where the figure stops. Task 4's smoke lookup then returns zero smoke below the table, with a logged note (smoke is negligible below the mesosphere in Megner's model).

`dlnV_dr2` is analytic. It must stay finite for L/r from 1e-6 to 1e6: use `erf` and `exp(-x²)` directly with no subtraction of nearly equal numbers, and switch to the series form when L/(2r) < 1e-4.

- [ ] **Step 1: Digitise the two figures** into the CSVs, viewing the rendered page.
- [ ] **Step 2: Write the failing tests**

```python
import math
import numpy as np
import pytest
from reentry_model.wake_column import transport as tr, species as sp
from reentry_model.wake_column.ambient import AmbientColumn

def test_eddy_profile_is_megner_reference():
    e = tr.EddyProfile.load()
    Ks = [e.K(z) for z in np.arange(50, 100, 2.5)]
    assert 5.0 <= max(Ks) <= 15.0 and min(Ks) >= 0.05

def test_molecular_diffusion_matches_spec_magnitudes():      # spec 4.1: ~0.01 m2/s at 50 km, ~3 at 90 km
    s, col = sp.load_species(), AmbientColumn.load("waccm", "night")
    assert 0.003 < tr.D_molecular(s["AlO"], col.at(50.0)) < 0.03
    assert 1.0 < tr.D_molecular(s["AlO"], col.at(90.0)) < 10.0

def test_switch_limits():
    assert tr.eddy_switch(10.0, 30.0) == 0.0 and tr.eddy_switch(60.0, 30.0) == 1.0
    assert 0.0 < tr.eddy_switch(30.0, 30.0) < 1.0

def test_volume_limits():
    assert tr.volume(1.0, 1e3) == pytest.approx(math.pi * 1e3, rel=1e-9)            # cylinder
    assert tr.volume(1e4, 1.0) == pytest.approx(math.pi ** 1.5 * 1e12, rel=1e-3)    # sphere
    assert tr.dlnV_dr2(1.0, 1e6) == pytest.approx(1.0, rel=1e-9)                     # d ln(pi r^2 L)/d r^2 = 1/r^2
    assert tr.dlnV_dr2(1.0, 1e-6) == pytest.approx(1.5, rel=1e-6)                    # d ln(r^3)/d r^2 = 1.5/r^2

def test_volume_finite_and_monotone_over_extreme_ratios():  # Review Focus 3
    r = np.logspace(-6, 6, 241)
    V = np.array([tr.volume(x, 1.0) for x in r]); g = np.array([tr.dlnV_dr2(x, 1.0) for x in r])
    assert np.all(np.isfinite(V)) and np.all(np.diff(V) > 0)
    assert np.all(np.isfinite(g)) and np.all((g * r**2 > 0.99) & (g * r**2 < 1.51))

def test_settling_rate():
    assert tr.settling_rate(0.03, 2.0) == pytest.approx(2 * 0.03 / (math.pi * 2.0))
```

- [ ] **Step 3: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_transport.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 4: Implement `transport.py`**, then run the tests. Expected: 6 passed.
- [ ] **Step 5: Commit** with `wake_column: Megner (2006) eddy-diffusion and smoke profiles digitised, hard-sphere Chapman-Enskog diffusion, the inner-scale eddy switch, and the line-segment plume volume that tends from cylinder to sphere`.

---

### Task 4: Particle microphysics primitives

**Files:**
- Create: `reentry_model/wake_column/particles.py`
- Test: `reentry_model/wake_column/tests/test_wake_particles.py`

**Interfaces:**
- Consumes: `Species` (Task 1), `AmbientState` (Task 2), `megner2006_smoke.csv` (Task 3).
- Produces:
  - `particles.Bins(bin_ratio=2.0, d_max_nm=400.0)`, with `.n_atoms: np.ndarray` (pivots, decision 8) and `.K: int`.
  - `.split(n_atoms) -> (k: int, zeta: float)`: the lower bin and the number fraction placed in it.
  - `.atom_fraction(n_atoms, k, zeta) -> float`: the share of the atoms going to bin k.
  - `particles.volume_per_al(h_per_al) -> m³` (decision 10).
  - `particles.diameter(n_al, h_per_al) -> m`.
  - `particles.epstein_speed(r_p, rho_p, rho_gas, T, g) -> m/s`, equal to `rho_p g r_p / (rho_gas c̄ (1 + π/8))`.
  - `particles.fm_kernel(r1, m1, r2, m2, T) -> m³/s`, equal to `π (r1 + r2)² sqrt(8 k_B T / (π μ))`.
  - `particles.capture_rate(a, b, T) -> m³/s` and `particles.spin_factor(a, b) -> float` (decision 7).
  - `particles.Pathway(pid: str, a: str, b: str)` and `particles.formation_pathways() -> list[Pathway]`: all 28 unordered pairs of `AL_SPECIES` including self-pairs, with pid "A+B" in `AL_SPECIES` order.
  - `particles.SmokeBackground.load()` with `.uptake_coefficient(sp, z_km, T) -> 1/s`: Σ over smoke radii of `n_smoke · fm_kernel(r_smoke, ∞-mass limit, r_mol, m_mol)`, per molecule, with sticking 1.

- [ ] **Step 1: Write the failing tests**

```python
import dataclasses
import math
import pytest
from reentry_model.wake_column import particles as pa, species as sp

def test_first_bin_is_a_half_nanometre_alumina_dimer():
    b = pa.Bins()
    assert b.n_atoms[0] == 2 and b.n_atoms[1] == 4
    assert pa.diameter(2, 0.0) == pytest.approx(0.4759e-9, rel=1e-3)
    assert pa.diameter(b.n_atoms[-1], 0.0) >= 400e-9 > pa.diameter(b.n_atoms[-2], 0.0)

@pytest.mark.parametrize("n", [2.0, 3.0, 5.0, 4095.0])
def test_split_conserves_number_and_atoms(n):
    b = pa.Bins(); k, z = b.split(n)
    assert z * b.n_atoms[k] + (1 - z) * b.n_atoms[k + 1] == pytest.approx(n)
    assert pa.Bins().atom_fraction(n, k, z) == pytest.approx(z * b.n_atoms[k] / n)

def test_hydroxide_is_bulkier_per_aluminium():
    assert pa.volume_per_al(3.0) / pa.volume_per_al(0.0) == pytest.approx(5.3524 / 2.8218, rel=1e-3)

def test_epstein_speed_hand_value():                        # r = 10 nm alumina at about 80 km
    assert pa.epstein_speed(1e-8, 3000.0, 1.846e-5, 198.6, 9.565) == pytest.approx(0.02929, rel=1e-3)

def test_dipole_capture_AlO_AlO_at_298K():                  # decision 7; Gomez Martin 2017 quote 1.1e-9 cm3/s
    s = sp.load_species()
    a = dataclasses.replace(s["AlO"], dipole_D=4.36)
    assert pa.capture_rate(a, a, 298.0) == pytest.approx(1.516e-15, rel=1e-3)

def test_unknown_dipole_falls_back_to_hard_sphere():
    s = sp.load_species(); a, b = s["OAlO"], s["AlO"]
    mu = a.mass_amu * b.mass_amu / (a.mass_amu + b.mass_amu) * 1.66053906660e-27
    hs = math.pi * (3.5e-10) ** 2 * math.sqrt(8 * 1.380649e-23 * 200 / (math.pi * mu))
    assert pa.capture_rate(a, b, 200.0) == pytest.approx(hs, rel=1e-9)

def test_spin_factors():
    s = sp.load_species()
    assert pa.spin_factor(s["AlO"], s["AlO"]) == 0.25
    assert pa.spin_factor(s["AlO"], s["AlOH"]) == 1.0
    assert pa.spin_factor(s["AlOH"], s["AlOH"]) == 1.0

def test_every_pair_is_its_own_pathway():
    p = pa.formation_pathways()
    assert len(p) == 28 and len({x.pid for x in p}) == 28 and "AlO+AlO" in {x.pid for x in p}

def test_smoke_uptake_is_positive_in_the_mesosphere():
    sm, s = pa.SmokeBackground.load(), sp.load_species()
    assert sm.uptake_coefficient(s["AlO"], 85.0, 200.0) > 0
```

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_particles.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 3: Implement `particles.py`.**
- [ ] **Step 4: Run the tests to see them pass.** Expected: 12 passed.
- [ ] **Step 5: Commit** with `wake_column: particle primitives -- aluminium-atom pivot bins with fixed-pivot splitting, composition-weighted size, Epstein settling, free-molecule kernel, dipole-dipole capture with spin factors for all 28 formation pairs, and uptake onto Megner smoke`.

---

### Task 5: Droplet adapter (the trail)

**Files:**
- Create: `reentry_model/wake_column/trail.py`
- Test: `reentry_model/wake_column/tests/test_wake_trail.py`

**Interfaces:**
- Consumes: `reentry_model.droplet_reentry.model.run_case(d0_um, h0_km)`, which returns `(final_state, rows)`. Each row has `t_s, h_km, V_ms, T_K, m_Al_kg, mdot_Al_kgs, m_Zn_kg, m_Mg_kg`.
- Produces:
  - `trail.EmissionSample(t_s, h_km, V_ms, T_d_K, weight)`, frozen.
  - `trail.Trail(d0_um, h0_km, L_m, N_Al, m_Al_kg, m_Zn_kg, m_Mg_kg, z_mean_km, samples: tuple[EmissionSample, ...])`, frozen. `N_Al` is the number of aluminium atoms vaporised.
  - `trail.trail_from_rows(rows, d0_um, h0_km, n_samples=20, loss_fraction=0.9) -> Trail` (decisions 15 and 16).
  - `trail.trail_for(d0_um, h0_km, n_samples=20, parent=None) -> Trail`, which calls `run_case`.
  - A run whose `rows` is None (skipped by the Weber check) or that loses no aluminium returns `N_Al = 0` and an empty `samples`.

- [ ] **Step 1: Write the failing tests**

```python
import pytest
from reentry_model.wake_column.trail import trail_for, trail_from_rows

@pytest.mark.parametrize("d, h, L", [(1.0, 80.0, 51.0), (10.0, 60.0, 40.0)])   # spec 4.2 table (2026-10-27 runs)
def test_trail_length_matches_design_measurements(parent, d, h, L):
    t = trail_for(d, h, parent=parent)
    assert t.L_m == pytest.approx(L, rel=0.05)
    assert abs(t.z_mean_km - h) < 0.05
    assert sum(s.weight for s in t.samples) == pytest.approx(1.0)

def test_inventory_is_the_whole_aluminium_loss(parent):
    from reentry_model.droplet_reentry.model import run_case
    fs, rows = run_case(10.0, 80.0, parent=parent)
    t = trail_from_rows(rows, 10.0, 80.0)
    lost = rows[0]["m_Al_kg"] - rows[-1]["m_Al_kg"]
    assert t.m_Al_kg == pytest.approx(lost, rel=1e-9)
    assert t.N_Al == pytest.approx(lost / (26.9815385e-3 / 6.02214076e23), rel=1e-9)
    assert all(100.0 < s.V_ms < 8000.0 and 900.0 < s.T_d_K < 2000.0 for s in t.samples)

def test_no_rows_gives_empty_trail():
    t = trail_from_rows(None, 1.0, 50.0)
    assert t.N_Al == 0 and t.samples == ()
```

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_trail.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 3: Implement `trail.py`.** The path length is `cumulative ∫ V_ms dt` (trapezoid). Every sample property is interpolated against cumulative aluminium loss.
- [ ] **Step 4: Run the tests to see them pass.** Expected: 4 passed, in a few seconds each.
- [ ] **Step 5: Commit** with `wake_column: the droplet adapter -- trail length to 90 % of the aluminium loss, the whole loss as the inventory, its mass-weighted altitude, and equal-mass emission samples for the slowing-down stage`.

---

### Task 6: Stage A — slowing-down Monte Carlo

**Files:**
- Create: `reentry_model/wake_column/slowing.py`
- Test: `reentry_model/wake_column/tests/test_wake_slowing.py`

**Interfaces:**
- Consumes: `load_reactions(stage=("hot","both"))`, `load_species`, `AmbientColumn.at`, `EmissionSample`, `Settings`.
- Produces:
  - `slowing.HotChannel(rid, reactant, partner, products: tuple, E0_J, P0)` (energies per molecule).
  - `slowing.calibrate(reactions, species) -> list[HotChannel]`:
    - E0 = `Ea_K · k_B` (or `E0_kJmol` for `form=threshold`).
    - P0 = `k(T_ref) exp(Ea/T_ref) / (π σ12² v̄(T_ref))` at T_ref = T_max, capped at 1.
    - Barrierless channels get E0 = 0. Three-body rows are skipped.
    - The expected P0 values are about 0.88 for R1, 0.082 for G42, 0.53 for R3 and 0.010 for G33.
  - `slowing.loc_rate(ch, T, species) -> m³/s`: the analytic thermal average `P0 π σ12² v̄ exp(−E0/k_B T)`.
  - `slowing.sample_rate(ch, T, species, n, rng) -> m³/s`: the same average estimated by sampling Maxwellian relative velocities. It is for test 10.
  - `slowing.HotStartResult(fractions: dict[str, float], r0_m, mean_collisions, mean_stop_time_s, share_above_ionisation, share_above_AlO_bond, per_sample: tuple[dict, ...])`.
  - `slowing.slow_down(samples, column: AmbientColumn, settings, species, channels, partners_at_rest=False, reactions_on=True) -> HotStartResult`.

Algorithm (decision 15 and spec §4.4):
- **Particles.** `n_test_particles` are split equally over the samples and advanced together, vectorised in numpy. Each sample uses `column.at(sample.h_km)`. The droplet velocity is V x̂ (zero when `hot_start == "droplet_temperature"`).
- **Free flight.** Use the null-collision method. The majorant rate is `Σ_j n_j π σ_j² (|v| + 4 v̄_j)`. A candidate partner j is drawn by weight; a candidate partner velocity is drawn from the Maxwellian at T_ambient (or is zero when `partners_at_rest`). The collision is accepted with probability `n_j π σ_j² |v − u| / majorant`.
- **Collisions.** Hard-sphere, isotropic in the centre-of-mass frame.
- **Reactions.** Before scattering, channel probabilities `scale · P0 (1 − E0/E)` (for E > E0) are drawn for the current species and partner. When their sum exceeds 1, they are renormalised to sum to 1. A product keeps the reactant's velocity and takes its own mass.
- **Stop.** A particle stops when `½ m |v|² < thermal_cutoff · 1.5 k_B T_ambient`. Record its species, its perpendicular distance `sqrt(y² + z²)`, its collision count and its elapsed time.
- **Diagnostics.** `share_above_ionisation` is the share of Al-atom collisions with E above 578 kJ/mol (ionisation of Al); `share_above_AlO_bond` is the share of AlO collisions above the AlO bond energy.
- **Outputs.** `r0_m` is the weighted root-mean-square perpendicular distance; `fractions` covers Al, AlO and OAlO, sums to 1 and is weighted over samples.
- **Seed.** `Settings.seed` seeds `numpy.random.default_rng`.

- [ ] **Step 1: Write the failing tests**

```python
import numpy as np
import pytest
from reentry_model.wake_column import slowing as sl, species as sp, reactions as rx
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.settings import Settings
from reentry_model.wake_column.trail import EmissionSample

S = sp.load_species()
CH = sl.calibrate(rx.load_reactions(stage=("hot", "both")), S)
NIGHT = AmbientColumn.load("waccm", "night")

def one(V, T, h=80.0):
    return (EmissionSample(0.0, h, V, T, 1.0),)

def test_calibrated_probabilities():
    P = {c.rid: c.P0 for c in CH}
    assert P["R1"] == pytest.approx(0.88, abs=0.03) and P["G42"] == pytest.approx(0.082, abs=0.01)
    assert P["D1"] == 1.0

def test_mean_energy_loss_per_collision_hard_sphere():        # spec test 9
    st = Settings(n_test_particles=20000)
    r = sl.slow_down(one(7600.0, 300.0), NIGHT, st, S, CH, partners_at_rest=True, reactions_on=False)
    assert r.per_sample[0]["mean_fractional_loss_vs_N2"] == pytest.approx(0.4998, abs=0.01)

def test_stopping_distance_scales_inversely_with_density():   # spec test 9
    st = Settings(n_test_particles=5000)
    a = sl.slow_down(one(7600.0, 1500.0, 70.0), NIGHT, st, S, CH, reactions_on=False)
    b = sl.slow_down(one(7600.0, 1500.0, 80.0), NIGHT, st, S, CH, reactions_on=False)
    ratio = NIGHT.at(70.0).n_total_m3 / NIGHT.at(80.0).n_total_m3
    assert b.r0_m / a.r0_m == pytest.approx(ratio, rel=0.15)

@pytest.mark.parametrize("rid", ["R1", "G42", "R3", "G33"])
def test_line_of_centers_reproduces_thermal_rate(rid):         # spec test 10
    ch = next(c for c in CH if c.rid == rid)
    T = next(x for x in rx.load_reactions(stage=("hot", "both")) if x.rid == rid).T_max
    rng = np.random.default_rng(1)
    assert sl.sample_rate(ch, T, S, 400000, rng) == pytest.approx(sl.loc_rate(ch, T, S), rel=0.05)

def test_cold_start_has_no_threshold_reactions():             # spec test 11
    T_amb = NIGHT.at(80.0).T_K
    r = sl.slow_down(one(0.0, T_amb), NIGHT, Settings(n_test_particles=5000), S, CH)
    assert r.fractions["OAlO"] == 0.0 and r.share_above_AlO_bond == 0.0
    lam = 1.0 / (NIGHT.at(80.0).n_total_m3 * np.pi * (3.711e-10) ** 2)
    assert 0.1 * lam < r.r0_m < 10 * lam

def test_packet_conservation_and_hot_products():               # spec test 12
    r = sl.slow_down(one(7600.0, 1800.0), NIGHT, Settings(n_test_particles=5000), S, CH)
    assert sum(r.fractions.values()) == pytest.approx(1.0, abs=1e-12)
    assert set(r.fractions) == {"Al", "AlO", "OAlO"}
    assert r.fractions["OAlO"] > 0.0 and r.share_above_ionisation < 1e-3

def test_degenerate_inputs_stay_finite():                      # Review Focus 4
    import dataclasses
    a = NIGHT.at(60.0)
    no_O = dataclasses.replace(a, n={**a.n, "O": 0.0})
    col = type("C", (), {"at": lambda self, z: no_O})()
    r = sl.slow_down(one(0.0, 908.0, 60.0), col, Settings(n_test_particles=2000), S, CH)
    assert np.isfinite(r.r0_m) and all(np.isfinite(v) for v in r.fractions.values())
```

`per_sample[i]["mean_fractional_loss_vs_N2"]` is the mean of ΔE/E over N2 collisions, recorded when `partners_at_rest`.

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_slowing.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 3: Implement `slowing.py`.**
- [ ] **Step 4: Run the tests to see them pass.** Expected: 10 passed, under a minute in total.
- [ ] **Step 5: Run the diagnostic and record it.** Run Stage A for 14.1 µm at 60, 76 and 90 km by day, using `trail_for`. Write the fractions, r0 and the two collision shares into the commit message. The design's rough estimate was "half or more AlO, a few percent OAlO at the highest speeds"; say whether it held.
- [ ] **Step 6: Commit** with `wake_column: Stage A slowing-down Monte Carlo -- hard-sphere test particles from the droplet's speed and temperature, line-of-centers reactions calibrated to the thermal rates, products and r.m.s. stopping distance for the cylinder; <the diagnostic numbers>`.

---

### Task 7: The box — gas phase, geometry and entrainment

**Files:**
- Create: `reentry_model/wake_column/box.py`
- Create: `reentry_model/wake_column/tests/box_helpers.py`, imported in tests as `from reentry_model.wake_column.tests.box_helpers import make_box, S`
- Test: `reentry_model/wake_column/tests/test_wake_box_gas.py`

**Interfaces:**
- Consumes: Tasks 1–6.
- Produces:
  - `box.Box(trail: Trail, hot: HotStartResult, amb: AmbientState, time_of_day, settings, *, species=None, reactions=None, bins=None, smoke=None, eddy=None, initial_bins: dict[int, float] | None = None)`.
  - `.layout`: a `box.Layout` with index slices `gas`, `r2`, `N`, `tally` (shape K × 7 sources × 2 routes), `O`, `H`, `smoke` (7 sources), `settled` (7 × 2), and `.unpack(y) -> dict`.
  - `.y0() -> np.ndarray` (scaled, decision 19).
  - `.rhs(t, y) -> np.ndarray`.
  - `.solve() -> box.Solution(t: np.ndarray, Y: np.ndarray, box)`.
  - `Solution.al_budget(i) -> dict[str, float]`: fractions of N_Al for each gas species and for `particles`, `smoke`, `settled`.
  - `Solution.concentration(name, i) -> m^-3`.
  - `Solution.r(i) -> m`, `Solution.V(i) -> m³`.
  - `Solution.bins(i) -> dict` of arrays: N, n_atoms, d_m, H_per_Al, O_per_Al and `tally[source][route]` in atoms.
  - `RateWarnings` messages are exposed as `.warnings`.

This task builds the gas part only:
- **Initial state.** Gas counts start at `fraction · N_Al` for Al, AlO and OAlO (from `hot.fractions`); ambient counts start at `c_amb · V0`, with `V0 = volume(r0, L)`.
- **Chemistry.** Every box reaction runs on concentrations (count / V) at the ambient T; three-body steps use [M] = n_total (decision 12); R17 runs only when `time_of_day == "day"`.
- **Entrainment.** Each ambient species gains `c_amb · dV/dt`.
- **Radius.** `d(r²)/dt = 4 D_total(r)` with D for AlO, which is zero when `freeze_radius`.
- **Switches.** `settings.reactions` restricts the network, and `settings.chemistry=False` removes it.
- **Particles.** The particle slices exist but stay zero until Task 8.

- [ ] **Step 1: Write the helper and the failing tests.** The helper is `make_box(*, z=80.0, time_of_day="night", N_Al=1e12, start="Al", r0=0.05, L=1e9, initial_bins=None, **settings_fields) -> Box`, with `S = load_species()`. It builds a `Trail` by hand (`d0_um=1`, `h0_km=z`, the given L and N_Al, no samples) and a `HotStartResult(fractions={start: 1.0}, r0_m=r0, …)` with zero diagnostics, so the tests run neither the droplet nor Stage A. The ambient is `AmbientColumn.load("waccm", time_of_day).at(z)`. `Box` exposes `.amb`, `.bins`, `.reaction(rid)` and `.settling_rate_of_bin(k)`.

```python
def test_inert_tracer_follows_r0_over_r_squared():            # spec test 2
    b = make_box(chemistry=False, particles=False, eddy_scale=0.0, t_end_s=1e3)
    sol = b.solve()
    D = transport.D_molecular(S["AlO"], b.amb)
    for i, t in enumerate(sol.t):
        r2 = 0.05**2 + 4 * D * t
        assert sol.r(i) ** 2 == pytest.approx(r2, rel=1e-6)
        assert sol.concentration("Al", i) * sol.V(i) == pytest.approx(1e12, rel=1e-6)

def test_ambient_stays_ambient_without_aluminium():             # spec test 3
    b = make_box(N_Al=0.0, t_end_s=1e3); sol = b.solve()
    for name in ("O", "O2", "O3", "H", "H2O", "CO2"):
        assert sol.concentration(name, -1) == pytest.approx(b.amb.n[name], rel=1e-9)

@pytest.mark.parametrize("z", [70.0, 80.0, 90.0])
def test_AlO_lifetime_against_R6(z):                            # spec test 4 (decision 12)
    b = make_box(z=z, start="AlO", reactions=("R6",), freeze_radius=True, particles=False, t_end_s=1e4)
    tau = 1.0 / (b.reaction("R6").k(b.amb.T_K) * b.amb.n_total_m3 * b.amb.n["O2"])
    sol = b.solve()
    i = np.searchsorted(sol.t, tau)
    assert sol.al_budget(i)["AlO"] == pytest.approx(np.exp(-sol.t[i] / tau), rel=1e-3)

def test_Al_lifetime_against_R1():
    b = make_box(reactions=("R1",), freeze_radius=True, particles=False, t_end_s=1.0)
    tau = 1.0 / (b.reaction("R1").k(b.amb.T_K) * b.amb.n["O2"])
    sol = b.solve(); i = np.searchsorted(sol.t, tau)
    assert sol.al_budget(i)["Al"] == pytest.approx(np.exp(-sol.t[i] / tau), rel=1e-3)

def test_gas_aluminium_is_conserved_with_full_chemistry():
    for tod in ("day", "night"):
        sol = make_box(time_of_day=tod, particles=False, t_end_s=3600.0).solve()
        for i in range(len(sol.t)):
            assert sum(sol.al_budget(i).values()) == pytest.approx(1.0, abs=1e-6)
```

In the R6 test, record the three lifetimes in the test output; the spec quotes about 0.08 s, 2 s and 40 s.

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_box_gas.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 3: Implement the gas part of `box.py`.** The output grid is log-spaced from `t_first_s` to `t_end_s` at `points_per_decade`, with t = 0 prepended. Use `solve_ivp(method="BDF", t_eval=grid, rtol=settings.rtol, atol=1e-12)`.
- [ ] **Step 4: Run the tests to see them pass.** Expected: 7 passed.
- [ ] **Step 5: Commit** with `wake_column: the box gas phase -- extensive counts in a diffusing line-segment plume with entrainment of ambient air, the full neutral network at ambient temperature, photolysis by day only; inert-tracer, ambient-stability and single-reaction lifetimes checked (R6: <the three lifetimes>)`.

---

### Task 8: The box — particles, smoke, coagulation, settling and provenance

**Files:**
- Modify: `reentry_model/wake_column/box.py`
- Test: `reentry_model/wake_column/tests/test_wake_box_particles.py`

**Interfaces:**
- Consumes: `Bins`, `capture_rate`, `spin_factor`, `formation_pathways`, `fm_kernel`, `epstein_speed`, `SmokeBackground`, `settling_rate` (Tasks 3–4); `Box` from Task 7.
- Produces: the same `Box` and `Solution` API, now with the particle terms. `Solution.bins(i)` and `al_budget(i)["particles"|"smoke"|"settled"]` become non-trivial. `Solution.provenance(i) -> dict[(source, route), fraction of N_Al]` covers routes "formation", "growth" and "smoke", and includes settled atoms under their original route.

Right-hand-side terms, all on concentrations at the current V and converted back to counts:
- **Formation.** For each pathway (A, B), the rate is `k_form · c_A · c_B · V`, with a self-pair factor ½ in the event rate. Each event removes one A and one B, adds one particle of 2 Al atoms to bin 0, and tallies one Al to (A, formation) and one to (B, formation). It adds O and H per decision 9. `k_form` follows decision 7 × `formation_eff[pid]` × `formation_scale`.
- **Growth.** For each species s and bin k, the rate is `fm_kernel(r_s, m_s, r_k, m_k) · c_s · N_k / V · sticking[s] · sticking_scale` (events per second in the plume). Each event adds one Al atom to a bin-k particle, which is redistributed by `Bins.split(n_k + 1)` and tallied (s, growth). The bin's existing tallies move with the atoms.
- **Oxygen uptake.** Per decision 11.
- **Smoke.** For each species s, the rate is `uptake_coefficient(s, z, T) · smoke_sticking · count_s`. The atoms go to `smoke[s]`. The smoke itself is not a state.
- **Coagulation.** Free-molecule kernel between all bin pairs (½ for self-pairs). The merged entity of `n_i + n_j` atoms is split by `Bins.split`; the tallies, O and H of both partners are carried, and the losses are removed from bins i and j in proportion.
- **Settling.** For each bin, the rate is `settling_rate(epstein_speed(r_k, ρ_k, ρ_gas, T, g), r)`, with g = 9.80665 · (6371/(6371 + z))². The bin's tallies move into `settled[source, route]`.
- **Switches.** `particles`, `coagulation`, `settling` and `smoke` each turn their terms off. `initial_bins` seeds bins with `N` particles whose atoms are tallied to (AlO, formation).

- [ ] **Step 1: Write the failing tests**

```python
def test_aluminium_conserved_in_full_run():                     # spec test 1
    for tod in ("day", "night"):
        sol = make_box(time_of_day=tod, N_Al=1e14, t_end_s=3600.0).solve()
        for i in range(len(sol.t)):
            assert sum(sol.al_budget(i).values()) == pytest.approx(1.0, abs=1e-6)

def test_coagulation_reaches_self_preserving_width():          # spec test 5 (decision 8 fallback below)
    b = make_box(N_Al=0.0, chemistry=False, settling=False, smoke=False, freeze_radius=True,
                 initial_bins={3: 1e15}, r0=1.0, L=1.0, t_end_s=1e6)
    sol = b.solve(); bins = sol.bins(-1)
    assert bins["N"].sum() < 1e-3 * 1e15
    lnd = np.log(bins["d_m"]); w = bins["N"] / bins["N"].sum(); mu = (w * lnd).sum()
    assert np.exp(np.sqrt((w * (lnd - mu) ** 2).sum())) == pytest.approx(1.45, abs=0.03)

def test_settling_only_decays_at_two_vs_over_pi_r():            # spec test 6
    b = make_box(N_Al=0.0, chemistry=False, coagulation=False, smoke=False, freeze_radius=True,
                 initial_bins={10: 1e6}, r0=2.0, t_end_s=1e3)
    sol = b.solve()
    rate = b.settling_rate_of_bin(10)
    for i in range(len(sol.t)):
        assert sol.bins(i)["N"][10] == pytest.approx(1e6 * np.exp(-rate * sol.t[i]), rel=1e-4)

def test_provenance_closes_in_every_bin():                      # spec test 7
    sol = make_box(time_of_day="day", N_Al=1e14, t_end_s=600.0).solve()
    for i in range(len(sol.t)):
        b = sol.bins(i)
        tallied = sum(b["tally"][s][r] for s in b["tally"] for r in b["tally"][s])
        assert np.allclose(tallied, b["N"] * b["n_atoms"], rtol=1e-6, atol=1e-6 * 1e14)
        assert sum(sol.provenance(i).values()) == pytest.approx(
            sum(v for k, v in sol.al_budget(i).items() if k in ("particles", "smoke", "settled")), abs=1e-9)

def test_every_pathway_and_the_smoke_route_are_reported():
    sol = make_box(time_of_day="night", N_Al=1e14, t_end_s=600.0).solve()
    prov = sol.provenance(-1)                                     # every (source, route) key is always present
    assert len(prov) == 7 * 3
    assert sum(v for (_, r), v in prov.items() if r == "formation") > 0
    assert sum(v for (_, r), v in prov.items() if r == "smoke") > 0
```

**Fallback for the coagulation test.** If test 5 fails at `bin_ratio=2` because the sectional scheme diffuses numerically, try `bin_ratio` 2^(1/2) and then 2^(1/3). Make the largest value that passes the default, and say which in the commit message and in decision 8 of this plan.

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_box_particles.py -q`. Expected: FAIL (zero particles, or a missing attribute).
- [ ] **Step 3: Implement the particle terms in `box.py`.** Precompute every kernel and split table once per run, because T and the altitude are fixed.
- [ ] **Step 4: Run the whole wake test suite.** Run `"$PY" -m pytest reentry_model/wake_column/tests -q`. Expected: everything passes, the Task 7 tests included.
- [ ] **Step 5: Commit** with `wake_column: particles in the box -- formation by all 28 pairs, growth, oxygen uptake, uptake onto background smoke, free-molecule coagulation and Epstein settling, with every aluminium atom tallied by source molecule and route; aluminium conserved, self-preserving width <value> at bin ratio <q>`.

---

### Task 9: Running one case and writing its outputs

**Files:**
- Create: `reentry_model/wake_column/run.py`
- Test: `reentry_model/wake_column/tests/test_wake_run.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `run.run(case: Case, settings=Settings(), out_dir=OUTPUT, trail: Trail | None = None, overwrite=False) -> dict`. It returns the summary dict and writes `out_dir/<run_name>/`.
  - `run.summarise(sol, trail, hot, case, settings) -> dict`.
  - `run.steady_state_time(t, reservoirs: dict[str, np.ndarray]) -> float | None` (decision 18).

Behaviour:
- If `<run_name>/summary.json` exists and not `overwrite`, load it and return it.
- Otherwise build the trail (unless given) with `trail_for(case.d0_um, case.h0_km, settings.n_emission_samples)`.
- If `N_Al == 0`, write a summary with `skipped: "no aluminium vaporised"` and the Zn and Mg masses, and return.
- Otherwise build `AmbientColumn.load(settings.ambient, case.time_of_day)`, `amb = column.at(trail.z_mean_km)`, then calibrate, `slow_down`, `Box(...).solve()`.
- Write `history.csv`, `size_distribution.csv` and `summary.json` (spec §6):
  - **history.csv** columns: `t_s, r_m, V_m3, T_K`; `frac_<species>` for the seven gas species; `frac_particles, frac_smoke, frac_settled, frac_condensed`; `frac_<source>_<route>` for every source and route; and `c_<ambient species>_m3` for O, O2, O3, H, H2O.
  - **size_distribution.csv** is written at 1e-3, 1, 60, 600, 3600 and 43200 s (the nearest grid times). Columns: `t_s, bin, n_atoms, d_nm, N, O_per_Al, H_per_Al, alumina_like, Al_<source>_<route>…`. `alumina_like` means H/Al below 0.5.
  - **summary.json** keys:
    - the case and the settings;
    - `run_name`;
    - the trail: `L_m`, `z_mean_km`, `N_Al`, `m_Al_kg`, `m_Zn_kg`, `m_Mg_kg`;
    - the hot start: `hot_fractions`, `r0_m`, `share_above_ionisation`, `share_above_AlO_bond`;
    - the outcome: `condensed_fraction_end`, `t50_s`, `t90_s`, `d_number_mean_nm`, `d_volume_mean_nm`, `alumina_like_share`, `provenance_end` (as "source/route" keys), `settled_fraction_end`, `smoke_fraction_end`, `steady_state_time_s`;
    - the provenance: `rate_warnings`, `git_commit` (from `git rev-parse HEAD`) and `runtime_s`.
  - The mean diameters are over bins at the final time (decision 17 defines "condensed").

- [ ] **Step 1: Write the failing tests**

```python
def test_short_run_writes_all_outputs(tmp_path, parent):
    st = Settings(t_end_s=60.0, n_test_particles=500, n_emission_samples=4)
    s = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    d = tmp_path / s["run_name"]
    assert (d / "history.csv").exists() and (d / "size_distribution.csv").exists()
    for k in ("condensed_fraction_end", "t50_s", "provenance_end", "hot_fractions", "r0_m", "steady_state_time_s"):
        assert k in s
    assert 0.0 <= s["condensed_fraction_end"] <= 1.0

def test_existing_summary_is_reused(tmp_path, parent):
    st = Settings(t_end_s=1.0, n_test_particles=200, n_emission_samples=2)
    a = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    b = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    assert a["runtime_s"] == b["runtime_s"]

def test_zero_inventory_is_skipped_cleanly(tmp_path):           # Review Focus 1
    t = Trail(1.0, 50.0, 0.0, 0.0, 0.0, 1e-15, 1e-16, 50.0, ())
    s = run.run(Case(1.0, 50.0, "day"), Settings(), out_dir=tmp_path, trail=t)
    assert s["skipped"] == "no aluminium vaporised"

def test_steady_state_time():
    t = np.logspace(-6, 4.6, 107)
    flat = {"a": np.where(t < 10, t / 10, 1.0)}
    assert run.steady_state_time(t, flat) == pytest.approx(10.0, rel=0.3)
    assert run.steady_state_time(t, {"a": np.log10(t) / 5}) is None
```

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_run.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 3: Implement `run.py`.**
- [ ] **Step 4: Run the tests to see them pass**, then time one full baseline case. Run `"$PY" -c "from reentry_model.wake_column.run import run; from reentry_model.wake_column.settings import Case; print(run(Case(14.1, 76.0, 'day'))['runtime_s'])"`. Expected: tests pass; record the runtime in the commit message, because Task 10's sweep size depends on it.
- [ ] **Step 5: Commit** with `wake_column: one-case runner -- droplet, trail, Stage A and the 12-hour box chained, with history, size-distribution and summary outputs under a configuration-encoding run name; 14.1 um at 76 km by day takes <runtime>`.

---

### Task 10: Sweeps, sensitivity set, plots and documentation

**Files:**
- Create: `reentry_model/wake_column/analysis/sweep.py`, `analysis/sensitivity.py`, `analysis/plot_sweep.py`, `analysis/plot_sensitivity.py`
- Modify: `docs/model_assumptions.md` (new section "9. Wake column"), `CLAUDE.md` ("The wake column" under the architecture, mirroring "The droplet stage"), `README.md` (one command block and the verification list), `reentry_model/wake_column/__init__.py` (module docstring listing the modules, as `droplet_reentry/__init__.py` does)
- Test: `reentry_model/wake_column/tests/test_wake_analysis.py`

**Interfaces:**
- Consumes: `run.run`, `Settings`, `Case`.
- Produces:
  - `sweep.default_cases() -> list[Case]`: 21 altitudes (90 to 50 km in 2 km steps) × (1.0, 14.1, 200.0) µm × (day, night), 126 cases.
  - `sweep.main(argv)` with `--ambient {waccm,msis_afgl}`, `--workers N` (multiprocessing pool, default os.cpu_count() − 1), `--diameters` and `--altitudes` overrides.
  - `sensitivity.settings_set() -> list[tuple[str, Settings]]`: the 13 one-at-a-time settings of spec §4.7, (a)–(h), labelled "formation x0.1", "formation x0.01", "sticking x0.1", "oxidation x0.1", "hot reactions x0.1", "hot reactions x3", "droplet-temperature start", "formation x80", "eddy x0.1", "eddy x100", "molecular only", "inner scale 10 m", "inner scale 100 m".
  - `sensitivity.cases() -> list[Case]`: 14.1 µm at 60, 76 and 90 km × day and night.
  - The plot scripts read only `summary.json` and `history.csv` files from `OUTPUT` and write PNGs to `OUTPUT/plots/`. They use `DRAMA Sphere Experiments/analysis/plot_style.py` when importable, as `compare.py` does.

Plots (spec §6):
- **plot_sweep:**
  - condensed fraction against time for each altitude, day and night panels, one figure per diameter;
  - provenance stack by altitude (end of run);
  - size distributions by altitude at 12 h;
  - day minus night condensed fraction;
  - WACCM against NRLMSIS + AFGL condensed fraction.
- **plot_sensitivity:** a tornado plot of the condensed fraction and the volume-mean diameter at 12 h, per case.

- [ ] **Step 1: Write the failing tests**

```python
from reentry_model.wake_column.analysis import sweep, sensitivity

def test_default_sweep_is_126_cases():
    c = sweep.default_cases()
    assert len(c) == 126 and len({(x.d0_um, x.h0_km, x.time_of_day) for x in c}) == 126

def test_sensitivity_set_is_13_settings_over_6_cases():
    s = sensitivity.settings_set()
    assert len(s) == 13 and len({run_name(Case(14.1, 76.0, "day"), x) for _, x in s}) == 13
    assert len(sensitivity.cases()) == 6 and {c.h0_km for c in sensitivity.cases()} == {60.0, 76.0, 90.0}
    labels = dict(s)
    assert labels["droplet-temperature start"].hot_start == "droplet_temperature"
    assert labels["molecular only"].eddy_scale == 0.0

def test_plot_scripts_run_on_a_tiny_sweep(tmp_path, parent):
    sweep.main(["--diameters", "1.0", "--altitudes", "80", "--workers", "1", "--out", str(tmp_path),
                "--t-end", "60", "--n-test-particles", "200"])
    from reentry_model.wake_column.analysis import plot_sweep
    plot_sweep.main(["--out", str(tmp_path)])
    assert list((tmp_path / "plots").glob("*.png"))
```

The `--out`, `--t-end` and `--n-test-particles` options exist for testing; the defaults are `OUTPUT` and the `Settings` defaults.

- [ ] **Step 2: Run them to see them fail.** Run `"$PY" -m pytest reentry_model/wake_column/tests/test_wake_analysis.py -q`. Expected: FAIL on the missing module.
- [ ] **Step 3: Implement the four scripts.**
- [ ] **Step 4: Write the documentation.**
  - `docs/model_assumptions.md` §9 lists every assumption of the spec's §10 and every decision 1–19 of this plan, each marked **(verified)** with its test name or **(assumed)**.
  - `CLAUDE.md` and `README.md` gain the wake column's commands and layout.
- [ ] **Step 5: Run the full Reentry Model unit tier.** Run `"$PY" -m pytest -m "not reference" -q`. Expected: the earlier 175 tests plus the wake tests all pass.
- [ ] **Step 6: Commit** with `wake_column: default sweep (126 cases, WACCM or NRLMSIS + AFGL), the 78-run sensitivity set, sweep and tornado plots, and the wake column's assumptions, commands and layout in model_assumptions, CLAUDE.md and the README`.
- [ ] **Step 7: Launch the default sweep in the background and report.** Run `"$PY" -m reentry_model.wake_column.analysis.sweep --ambient waccm` in the background. Report the expected wall time from the Task 9 runtime. Do not start the cross-check or sensitivity sweeps until Asha has seen the baseline plots.
