# Sub-plan 22 — spraying from the skins, and the flight (stage 4 of the layered-skin design) — scope

> Read `00-shared-context.md` first, then the spec. **This sub-plan is a scoped outline, not yet a code-level plan.**
> Its tasks are written out with tested code after sub-plan 21 is built and verified (spec §11.2).

**Goal:** Girin's spraying acts on each column's top cells, with the film thickness b the whole flowing depth; the old
film, deep account and their rules are removed; and the model flies the 100 mm and 50 mm cases to the spec's acceptance.

**Spec:** §8 (spraying), §9 (what is retired), §10.2 (outputs), §11 (build stage 4; level 3: convergence and the
acceptance criterion; the comparison with the element model), §14 (the assumptions the thesis states).

**Builds on:** sub-plan 21 — the flowing layer, the droplet momentum sink, the resolved sheared depth.

## Tasks (to be written in full)

1. **Spraying per column**: b the flowing depth, δ_m from Girin's Eq. (2), the rigid-substrate regime test reading the
   skin's depths, the thin, thick and front-surface Rayleigh–Taylor modes; the release leaves the top cells, liquid
   first; the droplets' momentum leaves with them. `--slurry-spray` carries over off. Tests: each mode's release from a
   prescribed column; books with spraying.
2. **Droplet records**: combined per column, per mode, per macro step; `solid_fraction` and `h_kg` added to the source
   table. Tests: combining changes no radius; histograms exact.
3. **Retirement**: element fill fractions, `PHI_DEATH`, deaths and hand-overs, the film and deep accounts, the feed rule
   and gate, the deep runoff, deep surfacing, the molten cascade, freeze-back and the stage-1 return device, removed from
   the skin path; `--surface-model elements` keeps the current model bit for bit.
4. **Outputs**: the history columns, frame fields and results of spec §10.2.
5. **Flights**: the 100 mm flight to 120 s and the 50 mm flight; convergence under halving of the macro step, 10 against
   20 µm cells, and a skin surface at 2 against 1 mm spacing; the acceptance criterion (droplet count, thick share,
   medians by number and mass, sprayed mass and mass at 120 s within the scatter between seeds under halving and
   refinement); the comparison with the element model's flight; cost and memory.
6. **Documentation**: the README, `docs/model_assumptions.md` and the facts, as the repository's rules require.

## Inputs to settle before the code-level plan is written

- Sub-plans 19–21's measurements: the fold limit's resolution, the rebuild transfer, the skins' and the sideways solve's
  costs, which together decide the default macro step and whether the flights fit the laptop or need the cluster.
