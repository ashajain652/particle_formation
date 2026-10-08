# Sub-plan: Task 2 — Material with latent heat, melting ranges and liquid properties

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 769–1108), **amended 2026-09-27: Scheil solidification as its own material variant**, and **amended 2026-09-28: Step 4's latent heat and surface tension in the Scheil variant, and the new `AA7075-empiricaldata` file** (sections below; the code blocks are the re-tested code), and **amended 2026-10-06: `Material.T_rigid`, the temperature at which the material is half liquid** (the first section below). Read `00-shared-context.md` first.

**Depends on:** Step 2's existing material loader only.
**Produces, for later tasks:** latent heat, the solidus/liquidus melting range (linear in `AA7075_range`, Scheil's law in the new `AA7075_scheil`), the empirical-data copy of AA7075 (`AA7075-empiricaldata`), feed fraction, and — critically — four separate enthalpy-related functions (mixed-phase enthalpy, liquid-only enthalpy, mixed-phase heat capacity, temperature-from-mixed-enthalpy) that Tasks 3, 6, and 9 all consume by name.
**Character:** physics/data — small in scope, mostly formulas plus four JSON data files (Step 3's two checked materials, the Scheil variant and the empirical-data copy of AA7075).
**Read before implementing:** Measured facts 4, 5, and 25 in the shared context describe two specific wrong implementations that were tried and rejected during prototyping (debiting only the destination element on feed; booking the film at mixture enthalpy instead of liquid-only enthalpy) — both reproduced real bugs (melt reappearing, a "free melting" runaway). Do not collapse the four enthalpy functions into one general-purpose formula; that simplification is exactly what caused the rejected version. Read the amendment of 2026-09-27 below as well: Scheil is a separate material variant, and `AA7075_range` must not be edited to implement it.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan — file list, function signatures, test plan, and acceptance criteria.

## Amendment of 2026-10-06 — `Material.T_rigid`, the temperature at which the material is half liquid

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements). The code below is the tested code, as a diff against this sub-plan's
> `material.py` of 2026-09-28 (the copy the measurements used carried it for the Scheil runs).

`RIGID_LIQUID_FRACTION = 0.5` and the property `Material.T_rigid`: the temperature at which the material's own liquid
fraction is one half, found by bisection (the liquid fraction is monotonic), infinite for a material that does not melt.
Below it mush is coherent and carries load, a rigid substrate for the melt film; above it mush is a slurry that flows —
Chen et al. 2016's semi-solid law ends at 50 % liquid and Li et al. 2014's slurry data begin there, the boundary Step 4
and the large-fragment design use. Values (fact 78): 829.0 K for `AA7075_range`, 895.10 K for `AA7075_scheil`, 850.0 K
for `AA7075`. Nothing else in the material changes; the existing tests are unchanged and the new one is appended after
the empirical-data test.

In `reentry_model/material.py`:

```diff
--- a/reentry_model/material.py
+++ b/reentry_model/material.py
@@ -35,6 +35,7 @@
 MELT_RAMP = 2.0                    # K, half-width of the numerical melting/feed ramps
 NO_MELT_ABOVE = 5000.0             # K: a melting temperature above this means "never melts"
 SCHEIL_DT = 1.0                    # K, spacing of the tabulated Scheil liquid fraction (interpolation error < 1e-3)
+RIGID_LIQUID_FRACTION = 0.5        # above this liquid fraction mush is a slurry, not a rigid substrate (amendment of 2026-10-06)
 
 
 @dataclass
@@ -119,6 +120,20 @@
         """Specific enthalpy of the film (liquid at the liquidus) [J/kg]."""
         return float(self.enthalpy(self.T_feed))
 
+    @property
+    def T_rigid(self):
+        """Temperature at which the material is RIGID_LIQUID_FRACTION liquid [K]; inf for a non-melting material. Below it
+        the mush is coherent and carries load, a rigid substrate for the melt film; above it the mush is a slurry that
+        flows (Chen et al. 2016's semi-solid law ends at 50 % liquid and Li et al. 2014's slurry data begin there;
+        amendment of 2026-10-06). Found by bisection on the material's own liquid fraction, which is monotonic."""
+        if not self.melts:
+            return np.inf
+        lo, hi = self.T_solidus, self.T_liquidus
+        for _ in range(100):
+            mid = 0.5 * (lo + hi)
+            lo, hi = (mid, hi) if self._liquid_fraction_raw(mid) < RIGID_LIQUID_FRACTION else (lo, mid)
+        return 0.5 * (lo + hi)
+
     @classmethod
     def from_drama_json(cls, path=None):
         """A DRAMA material file; `path` may also be a name in MATERIAL_NAMES (AA7075_nomelt, AA7075, AA7075_range,
```

Append to `tests/test_reentry_model_material.py`, after the empirical-data test:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-10-06: the rigid substrate -- the temperature at which the material is half liquid

def test_the_rigid_temperature_is_where_half_the_material_is_liquid():
    """T_rigid separates the coherent mush, which carries load, from the slurry, which flows: 50 % liquid, where Chen
    et al.'s (2016) semi-solid strength law ends and Li et al.'s (2014) slurry viscosity data begin. It follows each
    material's own liquid-fraction law."""
    assert material.RIGID_LIQUID_FRACTION == 0.5
    r = material.Material.from_drama_json("AA7075_range")
    assert r.T_rigid == pytest.approx(750.0 + 0.5 * 158.0, abs=1e-9)                     # 829 K on the linear range
    s = material.Material.from_drama_json("AA7075_scheil")
    assert s.T_rigid == pytest.approx(895.1, abs=0.05)                                    # Scheil: 13 K below the liquidus
    a = material.Material.from_drama_json("AA7075")
    assert a.T_rigid == pytest.approx(850.0, abs=1e-9)                                    # single temperature: mid-ramp
    for m in (r, s, a):
        assert float(m.liquid_fraction(m.T_rigid)) == pytest.approx(0.5, abs=1e-9)
    assert np.isinf(material.Material.from_drama_json("AA7075_nomelt").T_rigid)          # never melts: never slurry
```

## Amendment of 2026-09-28 — Step 4's values in the Scheil variant, and `AA7075-empiricaldata`

**Decisions (user, 2026-09-28).**
- **`AA7075_scheil` takes Step 4's values** (Step 4 plan, `docs/superpowers/plans/2026-09-27-deformation-ring-shedding.md`):
  latent heat **390 kJ/kg** instead of DRAMA's 400 (source: Modulus Metal's AA7075-T6 data sheet, which lists
  384–393 kJ/kg and cites no primary source; alternatives noted and not run: 397 kJ/kg for pure aluminium, NIST-JANAF,
  and 358 kJ/kg, Mills 2002 via ASM Handbook Vol. 15), and liquid surface tension **0.80 N/m** instead of 0.86
  (Bainbridge & Taylor 2013, Metall. Mater. Trans. A 44A, 3901–3909, doi 10.1007/s11661-013-1696-9, Table II:
  commercial 7075 by sessile drop at the liquidus + 50 K, 0.809 ± 0.041 N/m as melted in vacuum, 0.843 ± 0.018 after
  the oxide skin was broken, 0.777 ± 0.061 after exposure to dry air, 0.607 ± 0.083 after breaking it again in dry air).
  Its heat capacity above 850 K stays DRAMA's last value, 1131.6 J/(kg K), which the literature supports to about
  ±4 % for solid, mush and liquid (see the next item). Everything else in the variant is unchanged.
- **A new file, `reentry_model/data/materials/AA7075-empiricaldata.json`: a duplicate of `AA7075.json` whose only
  changes are the new surface tension and heat capacity values** — the liquid surface tension 0.80 N/m (as above) and
  the heat capacity table made explicit above 850 K with rows at 900, 1000, 1100 and 1200 K of 1131.6 J/(kg K), held
  beyond (the compiled liquid 7075 value is 1130 J/(kg K), Mills 2002 via ASM Vol. 15 Table 4; the mass-weighted sum of
  the NIST-JANAF / SGTE liquid element heat capacities gives 1132; liquid aluminium is 1177 in NIST-JANAF against 1127
  measured by pulse heating to 1491 K, Leitner et al. 2017, doi 10.1007/s11661-017-4053-6). Besides these two values
  it differs from `AA7075.json` only in its `name`, which must differ because run names encode the material, and its
  `_comment`. Its latent heat stays DRAMA's 400 kJ/kg and its melting temperature DRAMA's 850 K, as in `AA7075.json`.
  It is registered in `MATERIAL_NAMES` as `AA7075-empiricaldata`.
- **`AA7075.json` and `AA7075_range.json` stay exactly as they are** (byte-identical output of the generator).

**Measured on 2026-09-28** in a throwaway copy of the prototype (the 2026-09-27 copy with this amendment applied; the
code blocks below are the tested code):
- The generator writes `AA7075.json` and `AA7075_range.json` byte-identical to the checked files (`cmp`).
  `AA7075_scheil.json` differs from its 2026-09-27 version only in `meltingHeat`, `liquid.surfaceTension`,
  `liquid._sources` and `_comment`; `AA7075-empiricaldata.json` differs from `AA7075.json` only in `_comment`, `name`,
  `liquid.surfaceTension`, `liquid._sources` and the four heat-capacity rows above 850 K.
- `AA7075-empiricaldata`'s c_p(T) and k(T) are bit-identical to `AA7075`'s at 7201 temperatures from 200 to 2000 K
  (the explicit rows repeat the value `np.interp` already held), and its enthalpy agrees to 1e-12 relative.
- `tests/test_reentry_model_material.py`: 15 passed (the 14 of 2026-09-27, the Scheil ones updated to 390 kJ/kg and
  0.80 N/m, and the new one). Mutation check: writing 390 kJ/kg into the empirical file makes the new test fail. The
  whole unit tier: 206 passed, 1 skipped, 2 failed and 3 errors — the baseline of 205 plus the new test; the failures
  and errors are the known missing melting SESAM references of Task 11.
- 50 mm physics flight, the same command as below with `--material AA7075_scheil`, against the 2026-09-27 values
  (400 kJ/kg, 0.86 N/m): melt onset unchanged at 172.0 s / 79.06 km; demise 208.0 s / 67.98 km, 0.5 s earlier;
  sprayed mass 0.1804 kg (−0.2 %); 3.22e6 droplets (+22 %) with median radius 222 µm (−5 %), the smaller surface
  tension making the drops smaller; peak surface temperature 1304 K; absorbed heat 205.6 kJ (−0.8 %); energy-balance
  residual 1.0e-10 (exact); 2.33 Newton iterations per step; runtime 120 s.

**Follow-ups outside this sub-plan** (not edited here): sub-plan 13's `--material` help text should name
`AA7075_scheil` and `AA7075-empiricaldata`; Task 15 records both files and the 2026-09-28 values in the README,
`docs/model_assumptions.md` §9 and the spec's amendments; the file list in `00-shared-context.md` should name
`AA7075-empiricaldata.json`; the Step 4 plan's open decision on where its material values live is answered here.

## Amendment of 2026-09-27 — Scheil solidification as its own material variant

**Decision (user, 2026-09-27).** The body gets one Scheil material representation: the latent heat, the liquid
fraction and the viscosity law built on it all follow Scheil's non-equilibrium solidification curve. It is a
**separate material variant, `AA7075_scheil`, in a new file `reentry_model/data/materials/AA7075_scheil.json`.**
**`AA7075_range`, Step 3's checked material, stays exactly as it is:** the Step 1 generator writes `AA7075.json` and
`AA7075_range.json` byte-identical to the files it wrote before this amendment, the linear computations in
`material.py` are the same statements as before (they now sit in the `else` branches of new `if self.scheil:` checks), and the
existing tests are unchanged — the Scheil tests are appended after them. Do not implement Scheil by editing
`AA7075_range` or its file.

**The law.** f_l(T) = ((T_pure − T)/(T_pure − T_liquidus))^(1/(k − 1)) with partition coefficient k = 0.4 and
T_pure = 933 K (pure aluminium), and the same 908 K liquidus, 750 K solidus and 400 kJ/kg latent heat as
`AA7075_range` (the latent heat and the liquid surface tension became 390 kJ/kg and 0.80 N/m on 2026-09-28, above). It leaves f_l(750 K) = 3.6 % of eutectic liquid, which melts or freezes linearly over ±`MELT_RAMP`
around the solidus (748–752 K), so the variant's `T_solidus` is the ramp foot, 748 K, exactly as a
single-temperature material's is. The curve is tabulated every `SCHEIL_DT` = 1 K (largest deviation from the formula
below 1e-3) and its nodes join the enthalpy table, so the latent slope stays constant inside every interval and the
exact enthalpy inversion — the film's mixed enthalpy included — works unchanged. Half of the latent heat is released
in the 13 K below the liquidus (f_l = 0.5 at 895.1 K), against 8 % for the linear range; `T_feed` (910 K) and
`h_liquid` are the same as `AA7075_range`'s, because a fully liquid kilogram holds the same enthalpy either way
(since 2026-09-28 `h_liquid` is 10 kJ/kg lower, by the latent heat's difference).

**Known limitation, labelled.** Scheil describes a casting solidifying; the body is wrought 7075 melting, and
wrought 7075 first melts on heating at 769–818 K (Brehm et al. 2022, SAND2022-9908; Gu et al. 2023, Materials 16,
6145), not at 750 K. The variant therefore starts melting 20–70 K early, but with little liquid: 3.6 % at 750 K and
6 % at 800 K, against 32 % at 800 K for the linear range.

**Measured on 2026-09-27** in a throwaway copy of the prototype, with the variant's 400 kJ/kg and 0.86 N/m of that date
(kept as the record of that version; the code blocks below are now the 2026-09-28 code):
- The Step 1 generator writes `AA7075.json` and `AA7075_range.json` byte-identical to the checked files (`cmp`),
  and `AA7075_scheil.json` with the three Scheil keys added.
- The old and the new `material.py` give bit-identical results for `AA7075_nomelt`, `AA7075` and `AA7075_range`: the
  enthalpy tables, `enthalpy`, `cp_eff`, `liquid_fraction`, `temperature_from_enthalpy`, `enthalpy_liquid` and the
  film's mixed inverse compare exactly equal at 200 001 temperatures from 150 to 2500 K.
- `tests/test_reentry_model_material.py`: 14 passed (the 10 existing tests and the 4 new ones). The whole unit tier:
  205 passed, 1 skipped, 2 failed and 3 errors, which is the prototype's baseline of 201 passed plus the four new
  tests; the failures and errors are the known missing melting SESAM references of Task 11. A first version of the
  latent-heat test failed by 3 mJ/kg because its uniform trapezoid grid missed the c_p table's kinks (713, 733,
  753 ... 850 K); the test now integrates on a grid that contains them, and the material itself is exact to 1e-9 J/kg.
- 50 mm physics flight (the model proper, same code, only `--material` changed; `--velocity 7.5 --altitude 115
  --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --removal girin`):
  `AA7075_range` reproduces the unchanged prototype's run to a largest relative difference of 2.2e-8 in any of the 79
  history columns over 416 rows, the same order as fact 36's 1.6e-8 between two runs of one build. `AA7075_scheil`
  against `AA7075_range`: melt onset 172.0 s / 79.06 km against 174.5 s / 78.32 km (0.74 km higher, consistent with
  Scheil taking up only 6 % of the latent heat by 800 K instead of 32 %); demise 208.5 s / 67.83 km against
  207.5 s / 68.08 km; sprayed mass 0.1807 kg against 0.1792 kg (+0.8 %); 2.63e6 droplets against 3.05e6 (−14 %) with
  median radius 234 µm against 224 µm (+4.4 %); peak surface temperature 1303 K against 1411 K; absorbed heat
  207.3 kJ against 214.5 kJ; energy-balance residual 1.5e-10 (exact); 2.32 Newton iterations per step against 2.25;
  runtime 123 s against 126 s.

**Follow-ups outside this sub-plan** (not edited here; sub-plans 00, 13, 14 and 15 are being amended for the
surface-recession design at the same time).
- Task 13 (sub-plan 13) sets the `--melt on` default material, today `AA7075_range`. Making `AA7075_scheil` the
  default for melting runs, as approved on 2026-09-27, is a one-line change to its `--material` default and help text.
- The material-dependent measured facts in `00-shared-context.md` (the physics-mode results of facts 12, 24, 25, 28,
  32 and 35 among them) are re-measured by Task 14's verification runs once the default is `AA7075_scheil`; the
  bookkeeping device (`AA7075`) and the Girin cases do not depend on it.
- Task 15 records the variant in the README, `docs/model_assumptions.md` §9 and the spec's amendments, and the file
  list in `00-shared-context.md` should name `AA7075_scheil.json`.

---

### Task 2: Material with latent heat, melting ranges and liquid properties

**Files:**
- Modify (replace): `reentry_model/material.py`
- Create: `reentry_model/data/materials/AA7075.json`, `reentry_model/data/materials/AA7075_range.json`, `reentry_model/data/materials/AA7075_scheil.json` (the Scheil variant; `AA7075_range.json` is unchanged by it), `reentry_model/data/materials/AA7075-empiricaldata.json` (the empirical-data copy of `AA7075.json`; amendment of 2026-09-28)
- Test: `tests/test_reentry_model_material.py` (append)

**Interfaces:**
- Consumes: `reentry_model/data/materials/AA7075_nomelt.json` (Step 2).
- Produces: `Material` fields `latent_heat`, `T_solidus`, `T_liquidus`, `liquid: LiquidProperties(rho, mu, sigma)`; properties `melts`, `T_feed`, `h_liquid`; methods `liquid_fraction(T)`, `feed_fraction(T)`, `cp_eff(T)`, `enthalpy(T)` (exact, latent slope inside the range), `temperature_from_enthalpy(h)`, and the melt film's four: `enthalpy_liquid(T)` = h(T) + L_f (1 - f_l(T)) (what a kilogram of *liquid* holds at T -- the film carries its latent heat wherever it sits), `enthalpy_mixed(T, w)` and `cp_mixed(T, w)` for a node holding a fraction `w` of film and 1 - w of material, and `temperature_from_enthalpy_mixed(h, w)`, their exact inverse in T for every w (round-trip 4e-12 K, measured); `Material.from_drama_json(path=None)` accepting the names in `MATERIAL_NAMES` (`AA7075_nomelt`, `AA7075`, `AA7075_range`, `AA7075_scheil`, `AA7075-empiricaldata`); constants `MELT_RAMP = 2.0`, `NO_MELT_ABOVE = 5000.0`, `SCHEIL_DT = 1.0`. Scheil variant (amendment of 2026-09-27): optional `Material` fields `partition_coefficient` and `T_pure` (read from `solidification: "scheil"`, `partitionCoefficient`, `pureMeltingTemperature`), property `scheil`; for it `liquid_fraction` is Scheil's law tabulated at `SCHEIL_DT` with the eutectic ramp at the solidus, `cp_eff` is the tabulated slope, and `T_solidus` is the ramp foot (748 K); every other interface is the same for both range materials. Single-temperature materials get a ±MELT_RAMP ramp; `feed_fraction` is a ±MELT_RAMP ramp ending at `T_feed` (the liquidus, +2 K for range materials).

- [ ] **Step 1: Write the four material files**

Generate them from the packaged no-melt file (the DRAMA tables verbatim, the 1e5 K holding row dropped, the liquid properties added). `AA7075_scheil` is built from `AA7075_range` and adds the three Scheil keys, Step 4's latent heat (390 kJ/kg) and liquid surface tension (0.80 N/m) and its own comment; `AA7075-empiricaldata` is `AA7075` with only the liquid surface tension (0.80 N/m) and the heat capacity above 850 K (made explicit to 1200 K) changed, besides its name and comment:

```bash
"$PY" - <<'EOF'
import json
d = json.load(open("reentry_model/data/materials/AA7075_nomelt.json"))
base = {k: v for k, v in d.items() if k != "_comment"}
base["specificHeatCapacity"] = [r for r in d["specificHeatCapacity"] if r[0] <= 850.0]
base["heatConductivity"] = [r for r in d["heatConductivity"] if r[0] <= 850.0]
liquid = {"density": 2400.0, "viscosity": 1.3e-3, "surfaceTension": 0.86,
          "_sources": "pure aluminium near the liquidus: rho_l 2375-2400 kg/m3 (Smithells Metals Reference Book, 8th ed., Table 14.1), mu_l 1.2-1.4 mPa s at 933-1000 K (Assael et al. 2006, J. Phys. Chem. Ref. Data 35, 285), sigma 0.86-0.91 N/m (Smithells; ASM Handbook Vol. 2). Alloy corrections for AA7075 (5.6 % Zn, 2.5 % Mg, 1.6 % Cu) are within 10 %; used as assumptions (spec 2026-09-20 sections 2 and 16)."}
a = dict(base); a["name"] = "AA7075"; a["meltingTemperature"] = 850.0; a["liquid"] = liquid
a["_comment"] = "DRAMA 4.1.4 drama-AA7075 verbatim (TOOLS/material_database.xml: density, cp(T) and k(T) to 850 K, meltingHeat 400 kJ/kg, meltingTemperature 850 K, emissivity 0.4, catalycity 1, oxidation off) plus the liquid-phase properties DRAMA does not carry (`liquid`). Above 850 K the tables hold their last value (cp 1131.6 J/kgK, k 128.19 W/mK). Step 3 material (spec 2026-09-20 section 6): single melting temperature, numerical ramp of +-2 K in the model. DRAMA's database itself is untouched."
r = dict(a); r["name"] = "AA7075_range"; r["solidusTemperature"] = 750.0; r["liquidusTemperature"] = 908.0; r["meltingTemperature"] = 908.0
r["_comment"] = "AA7075 (see AA7075.json) with the alloy's melting range instead of DRAMA's single temperature: solidus 750 K, liquidus 908 K (ASM Handbook Vol. 2, Properties and Selection: Nonferrous Alloys, AA7075 477-635 C), latent heat 400 kJ/kg spread linearly across the range. meltingTemperature is set to the liquidus. Default material of --melt on (spec 2026-09-20 sections 6 and 14); melt and runoff start at the liquidus (section 8)."
# the Step 4 values of 2026-09-28: liquid surface tension, and the heat capacity above 850 K made explicit
SIGMA_7075, SIGMA_SOURCE = 0.80, "sigma 0.80 N/m for fully liquid AA7075 (user decision 2026-09-28): Bainbridge & Taylor 2013, Metall. Mater. Trans. A 44A, 3901-3909, doi 10.1007/s11661-013-1696-9, Table II, commercial 7075 by sessile drop at the liquidus + 50 K: 0.809 +- 0.041 N/m as melted in vacuum (0.843 +- 0.018 after the oxide skin was broken, 0.777 +- 0.061 after exposure to dry air, 0.607 +- 0.083 after breaking it again in dry air); held constant with temperature."
CP_ABOVE_850 = [[900.0, 1131.6], [1000.0, 1131.6], [1100.0, 1131.6], [1200.0, 1131.6]]
CP_SOURCE = "Above 850 K the heat capacity is DRAMA's last value, 1131.6 J/kgK, made explicit to 1200 K and held beyond, for solid, mush and liquid alike, uncertainty about +-4 % (user decision 2026-09-28): the compiled liquid 7075 value is 1130 J/kgK (Mills 2002, reprinted in ASM Handbook Vol. 15, Table 4), the mass-weighted sum of the NIST-JANAF / SGTE liquid heat capacities of Al, Zn, Mg, Cu and Cr gives 1132 J/kgK, and liquid aluminium is 1177 J/kgK in NIST-JANAF against 1127 J/kgK measured by pulse heating to 1491 K (Leitner et al. 2017, doi 10.1007/s11661-017-4053-6). This is the sensible heat capacity; the latent heat is added separately."
liquid_step4 = dict(liquid); liquid_step4["surfaceTension"] = SIGMA_7075
liquid_step4["_sources"] = liquid["_sources"].replace("sigma 0.86-0.91 N/m (Smithells; ASM Handbook Vol. 2)", "sigma see below") + " " + SIGMA_SOURCE
s = dict(r); s["name"] = "AA7075_scheil"; s["solidification"] = "scheil"; s["partitionCoefficient"] = 0.4; s["pureMeltingTemperature"] = 933.0
s["meltingHeat"] = 390000.0; s["liquid"] = liquid_step4
s["_comment"] = "AA7075_range (see AA7075_range.json, which stays as it is) with the latent heat released along Scheil's non-equilibrium solidification curve instead of linearly: f_l = ((T_pure - T)/(T_pure - T_liquidus))^(1/(k - 1)), partition coefficient k 0.4, T_pure 933 K (pure aluminium), liquidus 908 K, solidus 750 K. Half of the latent heat is released in the 13 K below the liquidus (f_l = 0.5 at 895.1 K) against 8 % for the linear range; the eutectic remainder f_l(750 K) = 3.6 % melts or freezes over +-2 K at the solidus. One Scheil material for the whole body -- heat, liquid fraction and the viscosity law that uses it (user decision 2026-09-27). Step 4 values (user decision 2026-09-28): latent heat 390 kJ/kg instead of DRAMA's 400 (Modulus Metal AA7075-T6 data sheet, 384-393 kJ/kg, which cites no primary source; alternatives 397 kJ/kg pure aluminium NIST-JANAF, 358 kJ/kg Mills 2002) and liquid surface tension 0.80 N/m (see liquid._sources); the heat capacity above 850 K is held at 1131.6 J/kgK, which the literature supports (see AA7075-empiricaldata.json). Known limitation, labelled: wrought 7075 first melts on heating at 769-818 K (Brehm et al. 2022, SAND2022-9908; Gu et al. 2023, Materials 16, 6145), not at the 750 K solidus of a casting, so melting starts 20-70 K early -- though with only 3.6 % liquid at 750 K and 6 % at 800 K, against 32 % at 800 K for the linear range."
e = dict(a); e["name"] = "AA7075-empiricaldata"; e["liquid"] = liquid_step4
e["specificHeatCapacity"] = a["specificHeatCapacity"] + CP_ABOVE_850
e["_comment"] = "AA7075 (see AA7075.json) duplicated with only two changes, both empirical values adopted on 2026-09-28 (user decision): the liquid surface tension, 0.80 N/m instead of 0.86 (see liquid._sources), and the heat capacity above 850 K. " + CP_SOURCE + " Everything else -- DRAMA's tables to 850 K, the 400 kJ/kg latent heat, the 850 K melting temperature, the conductivity, emissivity and liquid density and viscosity -- is AA7075's."
for name, doc in (("AA7075", a), ("AA7075_range", r), ("AA7075_scheil", s), ("AA7075-empiricaldata", e)):
    keys = ["_comment", "name", "materialType", "catalycity", "density", "meltingHeat", "meltingTemperature"] + (["solidusTemperature", "liquidusTemperature"] if "solidusTemperature" in doc else []) + (["solidification", "partitionCoefficient", "pureMeltingTemperature"] if "solidification" in doc else []) + ["specificHeatCapacity", "heatConductivity", "emissivity", "oxideActivationTemperature", "oxideEmissivity", "oxideHeatOfFormation", "oxideReactionProbability", "liquid"]
    json.dump({k: doc[k] for k in keys}, open("reentry_model/data/materials/%s.json" % name, "w"), indent=2)
print("written")
EOF
```

Expected: all four files exist; `AA7075.json` has 29 c_p rows ending at 850 K, `meltingHeat` 400000, `meltingTemperature` 850; `AA7075_range.json` adds `solidusTemperature` 750 and `liquidusTemperature` 908 (`meltingTemperature` 908); `AA7075_scheil.json` is `AA7075_range.json` plus `solidification` "scheil", `partitionCoefficient` 0.4 and `pureMeltingTemperature` 933, with `meltingHeat` 390000 and `liquid.surfaceTension` 0.80; `AA7075-empiricaldata.json` is `AA7075.json` with `liquid.surfaceTension` 0.80 and 33 c_p rows ending at 1200 K. `AA7075.json` and `AA7075_range.json` are byte-identical to what this script wrote before the Scheil amendment (checked with `cmp` on 2026-09-27 and again on 2026-09-28).

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_reentry_model_material.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: latent heat, melting ranges, feed fraction, liquid properties

def test_material_names_and_liquid_properties():
    a, r, n = (material.Material.from_drama_json(k) for k in ("AA7075", "AA7075_range", "AA7075_nomelt"))
    assert a.melts and r.melts and not n.melts and n.latent_heat == 0.0 and n.liquid is None
    assert a.latent_heat == r.latent_heat == 400e3 and a.rho == r.rho == 2813.0 and a.emissivity == 0.4
    assert (a.T_solidus, a.T_liquidus) == (848.0, 852.0) and (r.T_solidus, r.T_liquidus) == (750.0, 908.0)
    assert a.liquid.rho == 2400.0 and a.liquid.mu == 1.3e-3 and a.liquid.sigma == 0.86
    assert np.allclose(a.k(a.T_k), n.k(a.T_k)) and np.allclose(a.cp(a.T_cp), n.cp(a.T_cp))    # DRAMA's tables verbatim
    assert a.T_cp.max() == 850.0 and n.T_cp.max() == 1e5


def test_liquid_and_feed_fractions():
    a, r = material.Material.from_drama_json("AA7075"), material.Material.from_drama_json("AA7075_range")
    assert np.allclose(a.liquid_fraction([840.0, 848.0, 850.0, 852.0, 900.0]), [0.0, 0.0, 0.5, 1.0, 1.0])
    assert np.allclose(a.feed_fraction([840.0, 848.0, 850.0, 852.0, 900.0]), [0.0, 0.0, 0.5, 1.0, 1.0])     # single T: feed = liquid fraction
    assert np.allclose(r.liquid_fraction([750.0, 829.0, 908.0]), [0.0, 0.5, 1.0])
    assert np.allclose(r.feed_fraction([900.0, 906.0, 908.0, 910.0]), [0.0, 0.0, 0.5, 1.0])              # feed only at the liquidus (+-2 K)
    assert r.T_feed == 910.0 and a.T_feed == 852.0
    assert a.cp_eff(850.0) == pytest.approx(a.cp(850.0) + 400e3 / 4.0) and r.cp_eff(800.0) == pytest.approx(r.cp(800.0) + 400e3 / 158.0)
    assert a.cp_eff(700.0) == a.cp(700.0) and r.cp_eff(950.0) == r.cp(950.0)


def test_enthalpy_jump_is_exact_and_invertible():
    for name in ("AA7075", "AA7075_range"):
        m = material.Material.from_drama_json(name)
        sensible = np.trapezoid(m.cp(np.linspace(m.T_solidus, m.T_liquidus, 2001)), np.linspace(m.T_solidus, m.T_liquidus, 2001))
        assert m.enthalpy(m.T_liquidus) - m.enthalpy(m.T_solidus) == pytest.approx(400e3 + sensible, rel=1e-9)
        T = np.linspace(200.0, 2000.0, 7201)
        h = m.enthalpy(T)
        assert np.all(np.diff(h) > 0.0) and np.abs(m.temperature_from_enthalpy(h) - T).max() < 1e-8
        assert m.h_liquid == pytest.approx(m.enthalpy(m.T_feed))
    n = material.Material.from_drama_json("AA7075_nomelt")
    assert n.enthalpy(1000.0) == pytest.approx(material.Material.from_drama_json("AA7075").enthalpy(1000.0) - 400e3, rel=1e-9)


def test_invalid_melting_data():
    with pytest.raises(ValueError):
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], -1.0, 850.0, 850.0)
    with pytest.raises(ValueError):
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 900.0, 850.0)


def test_liquid_and_mixed_enthalpy_invert_exactly(request):
    """The melt film is liquid: h_liquid(T) = h(T) + L_f (1 - f_l) carries the latent heat at every temperature, and a
    node holding a fraction w of film has a latent plateau only (1 - w) as tall. enthalpy_mixed and its inverse are
    each other's exact inverse for every w -- the solver's Newton update inverts the node's own mixture, and a
    mismatch there limit-cycled the iteration between 777 K and 864 K (Step 3)."""
    for name in ("AA7075", "AA7075_range"):
        mat = material.Material.from_drama_json(name)
        T = np.linspace(300.0, 1400.0, 2001)
        solid, liquid = T < mat.T_solidus, T > mat.T_liquidus
        assert mat.enthalpy_liquid(T)[solid] == pytest.approx(mat.enthalpy(T)[solid] + mat.latent_heat)
        assert mat.enthalpy_liquid(T)[liquid] == pytest.approx(mat.enthalpy(T)[liquid])
        assert mat.enthalpy_mixed(T, 0.0) == pytest.approx(mat.enthalpy(T)) and mat.enthalpy_mixed(T, 1.0) == pytest.approx(mat.enthalpy_liquid(T))
        assert mat.cp_mixed(T, 0.0) == pytest.approx(mat.cp_eff(T)) and mat.cp_mixed(T, 1.0) == pytest.approx(mat.cp(T))
        w = np.random.default_rng(0).random(T.size)
        for weights in (np.zeros_like(T), np.full_like(T, 0.3), np.ones_like(T), w):
            assert mat.temperature_from_enthalpy_mixed(mat.enthalpy_mixed(T, weights), weights) == pytest.approx(T, abs=1e-9)
        assert np.all(np.diff(mat.enthalpy_mixed(T, 0.4)) > 0.0)                     # monotone: the inverse is a function
```

Then append the Scheil variant's tests (amendment of 2026-09-27, values updated 2026-09-28) after them, leaving the tests above unchanged:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-09-27: the Scheil material, a separate variant (AA7075_range stays as it is)

def scheil_fl(T, k=0.4, T_pure=933.0, T_liq=908.0):
    return ((T_pure - np.asarray(T, dtype=float)) / (T_pure - T_liq)) ** (1.0 / (k - 1.0))


def test_scheil_material_is_a_separate_variant():
    s, r = material.Material.from_drama_json("AA7075_scheil"), material.Material.from_drama_json("AA7075_range")
    assert s.scheil and not r.scheil and s.melts and s.partition_coefficient == 0.4 and s.T_pure == 933.0
    assert (r.T_solidus, r.T_liquidus) == (750.0, 908.0)                     # the range material is untouched
    assert (s.T_solidus, s.T_liquidus) == (748.0, 908.0)                     # foot of the eutectic ramp at the solidus
    assert s.latent_heat == 390e3 and r.latent_heat == 400e3 and s.rho == r.rho     # Step 4's latent heat (2026-09-28)
    assert s.liquid.sigma == 0.80 and r.liquid.sigma == 0.86                         # Step 4's surface tension (2026-09-28)
    assert s.liquid.rho == r.liquid.rho and s.liquid.mu == r.liquid.mu
    assert np.array_equal(s.T_cp, r.T_cp) and np.array_equal(s.cp_table, r.cp_table) and np.array_equal(s.k_table, r.k_table)
    assert s.T_feed == r.T_feed == 910.0 and s.h_liquid == pytest.approx(r.h_liquid - 10e3, rel=1e-12)   # fully liquid: L less


def test_scheil_liquid_fraction_and_latent_heat():
    s = material.Material.from_drama_json("AA7075_scheil")
    assert s.liquid_fraction(908.0) == 1.0 and s.liquid_fraction(747.9) == 0.0 and s.liquid_fraction(920.0) == 1.0
    assert scheil_fl(750.0) == pytest.approx(0.0362, abs=1e-4)               # the eutectic remainder at the solidus
    assert s.liquid_fraction(752.0) == pytest.approx(scheil_fl(752.0), rel=1e-12)      # top of the eutectic ramp
    assert s.liquid_fraction(750.0) == pytest.approx(0.5 * scheil_fl(752.0), rel=1e-12)
    T = np.arange(752.0, 908.0, 0.37)
    assert np.abs(s.liquid_fraction(T) - scheil_fl(T)).max() < 1e-3          # the 1 K table against the law
    assert s.liquid_fraction(895.1) == pytest.approx(0.5, abs=1e-3)          # half liquid 12.9 K below the liquidus
    def sens(a, b):                         # exact for the piecewise-linear c_p: the grid holds its kinks
        x = np.union1d(np.linspace(a, b, 2001), s.T_cp[(s.T_cp > a) & (s.T_cp < b)])
        return np.trapezoid(s.cp(x), x)
    latent = lambda a, b: s.enthalpy(b) - s.enthalpy(a) - sens(a, b)
    assert latent(748.0, 908.0) == pytest.approx(390e3, rel=1e-9)            # all of L_f, exactly
    assert latent(895.1, 908.0) / 390e3 == pytest.approx(0.5, abs=1e-3)      # half of it in the top 13 K ...
    r = material.Material.from_drama_json("AA7075_range")
    assert (r.enthalpy(908.0) - r.enthalpy(895.1) - sens(895.1, 908.0)) / 400e3 == pytest.approx(12.9 / 158.0, rel=1e-6)   # ... 8 % linearly


def test_scheil_enthalpy_is_exact_monotone_and_invertible():
    s = material.Material.from_drama_json("AA7075_scheil")
    T = np.linspace(200.0, 2000.0, 7201)
    h = s.enthalpy(T)
    assert np.all(np.diff(h) > 0.0) and np.abs(s.temperature_from_enthalpy(h) - T).max() < 1e-8
    Tm = np.array([749.3, 760.5, 820.5, 890.5, 907.5])                        # cp_eff is dh/dT inside every interval
    assert s.cp_eff(Tm) == pytest.approx((s.enthalpy(Tm + 1e-4) - s.enthalpy(Tm - 1e-4)) / 2e-4, rel=1e-6)
    assert s.cp_eff(700.0) == s.cp(700.0) and s.cp_eff(950.0) == s.cp(950.0)
    solid, liquid = T < s.T_solidus, T > s.T_liquidus
    assert s.enthalpy_liquid(T)[solid] == pytest.approx(s.enthalpy(T)[solid] + s.latent_heat)
    assert s.enthalpy_liquid(T)[liquid] == pytest.approx(s.enthalpy(T)[liquid])
    assert s.cp_mixed(T, 0.0) == pytest.approx(s.cp_eff(T)) and s.cp_mixed(T, 1.0) == pytest.approx(s.cp(T))
    w = np.random.default_rng(0).random(T.size)
    for weights in (np.zeros_like(T), np.full_like(T, 0.3), np.ones_like(T), w):
        assert s.temperature_from_enthalpy_mixed(s.enthalpy_mixed(T, weights), weights) == pytest.approx(T, abs=1e-9)
    assert np.all(np.diff(s.enthalpy_mixed(T, 0.4)) > 0.0)


def test_invalid_scheil_data():
    args = ("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 750.0, 908.0, None)
    for k, T_pure in ((0.0, 933.0), (1.0, 933.0), (0.4, 900.0)):
        with pytest.raises(ValueError):
            material.Material(*args, k, T_pure)
    with pytest.raises(ValueError):                                          # range narrower than the eutectic ramp
        material.Material("bad", 2813.0, 0.4, [300.0, 900.0], [900.0, 900.0], [300.0, 900.0], [150.0, 150.0], 4e5, 850.0, 852.0, None, 0.4, 933.0)
```

Then append the empirical-data file's test (amendment of 2026-09-28):

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-09-28: AA7075-empiricaldata, AA7075 with only the empirical surface tension and heat capacity

def test_empirical_material_is_aa7075_with_two_changes():
    import json
    da = json.load(open(material.MATERIAL_NAMES["AA7075"]))
    de = json.load(open(material.MATERIAL_NAMES["AA7075-empiricaldata"]))
    assert {k for k in da if da[k] != de[k]} == {"_comment", "name", "liquid", "specificHeatCapacity"}
    assert {k for k in da["liquid"] if da["liquid"][k] != de["liquid"][k]} == {"surfaceTension", "_sources"}
    assert de["specificHeatCapacity"][:len(da["specificHeatCapacity"])] == da["specificHeatCapacity"]
    assert de["specificHeatCapacity"][len(da["specificHeatCapacity"]):] == [[900.0, 1131.6], [1000.0, 1131.6], [1100.0, 1131.6], [1200.0, 1131.6]]
    e, a = material.Material.from_drama_json("AA7075-empiricaldata"), material.Material.from_drama_json("AA7075")
    assert e.name == "AA7075-empiricaldata" and e.liquid.sigma == 0.80 and a.liquid.sigma == 0.86
    assert (e.latent_heat, e.T_solidus, e.T_liquidus, e.rho, e.emissivity) == (a.latent_heat, a.T_solidus, a.T_liquidus, a.rho, a.emissivity)
    T = np.linspace(200.0, 2000.0, 7201)
    assert np.array_equal(e.cp(T), a.cp(T)) and np.array_equal(e.k(T), a.k(T))   # the held value made explicit: same c_p
    assert e.enthalpy(T) == pytest.approx(a.enthalpy(T), rel=1e-12) and np.abs(e.temperature_from_enthalpy(e.enthalpy(T)) - T).max() < 1e-8
```


- [ ] **Step 3: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py -q`
Expected: the ten new tests fail (`from_drama_json("AA7075")`, `from_drama_json("AA7075_scheil")` and `MATERIAL_NAMES["AA7075-empiricaldata"]` are bad paths or names; no `melts`, `feed_fraction`, `liquid`, `enthalpy_liquid`; `Material` takes no Scheil arguments).

- [ ] **Step 4: Replace `reentry_model/material.py`**

```python
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
```


- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py tests/test_reentry_model_thermal.py tests/test_reentry_model_cli.py -q`
Expected: all pass (the no-melt material behaves exactly as in Step 2: `latent_heat` 0, `melts` False, enthalpy unchanged; `AA7075_range` behaves exactly as before the Scheil amendment).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/material.py reentry_model/data/materials/AA7075.json reentry_model/data/materials/AA7075_range.json reentry_model/data/materials/AA7075_scheil.json reentry_model/data/materials/AA7075-empiricaldata.json tests/test_reentry_model_material.py
git commit -m "Add the latent heat, melting ranges (linear and Scheil), feed fraction, liquid properties and the empirical-data copy of AA7075 to the material (Step 3 Task 2)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

