# Gas loads: reference calculations and comparison plots

Evidence behind Part A of `docs/superpowers/specs/2026-11-12-shoulder-melt-loss-design.md` (interim gas loads),
made on 2026-11-12 before the loads were built into the package. Everything here is a stand-alone script; nothing in
`reentry_model/` imports it. The package implementation is tested against the CSV copies in
`tests/fixtures/gas_loads/`.

## `reference/`

Run each with the model's environment from this folder, for example
`"$PY" -I similarity_heat.py`; each writes its CSV next to itself.

- `hb_dsmc.csv`, `hb_dsmc_flags.csv` — Holman & Boyd (2011, Phys. Fluids 23, 027101) DSMC surface pressure, shear and
  heat-flux coefficients digitized from their Figs 5-6, 8-9, 11-12 and 14-15 at 1-degree spacing for Kn 0.002, 0.01,
  0.05 and 0.25 (0.3048 m sphere, Mach 25 reacting air at 7108 m/s, 200 K, 1500 K wall). Checks: their tabulated
  drag within 2.4 %, their tabulated stagnation heating within 3 %. Their rear "Cp" is p/(0.5 rho U^2) (design A1).
- `ds_model.py` -> `ds_model.csv`, `ds_params.csv` — Daniels & Schwartzentruber (2026, arXiv:2609.28658) Appendix A,
  all regimes, perfect gas, fully diffuse wall; `verify_fig7/9/10.png` overlay it on the paper's own figures.
- `model_loads.py` -> `model_loads*.csv`, `model_loads_summary.json` — the re-entry model's loads before this change
  (`SurfaceFlow(rarefied_shear="slip", gamma_pm=1.15)`, `PhysicsHeating(fay-riddell, matting)`) at the same conditions.
- `similarity_shear.py` — Lees local-similarity shear (shear-only form) on the model's edge state; the class
  `Similarity` passes the Falkner-Skan checks 0.9277, 0.4696, 1.2326.
- `similarity_hb_pressure.py` (`frozen` argument for unreacted gas) — the same on Holman & Boyd's pressure.
- `similarity_lee.py` -> `similarity_lee.csv` — round the whole body on the decided pressure (model windward pressure,
  Prandtl-Meyer continued past 90 deg until it meets Daniels & Schwartzentruber's lee pressure), unreacted gas.
- `similarity_heat.py` -> `similarity_heat.csv` — the decided single solution: total-enthalpy energy equation with
  dissipation and real frozen-air properties, giving shear and heating together, to separation; the model's Matting
  bridge on the heating.

## `plots/`

`windward_loads_comparison.py`, `leeward_loads_comparison.py`, `leeward_shear_options.py`, `heating_comparison.py`
draw the four figures in the design (copies in `docs/superpowers/specs/figures/2026-11-12-shoulder-melt-loss/`) from
`reference/`, into `output/gas_loads/` (git-ignored).
