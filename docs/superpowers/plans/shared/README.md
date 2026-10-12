# Shared physics sub-plans (stage S)

These sub-plans promote Step 3's architecture-independent physics from the prototype work copies into
`reentry_model/shared/`, so that both model architectures — **skyscrapers** and **level set** — import one copy.
Design of record: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4.2 and §14.
They are built after stage 0 (the shared package, the re-exports and the `--architecture` flag) and before either
architecture's melting code lands in the package. A change to a shared module must keep both architectures' tests
passing. The skins are never shared.

| Sub-plan | Shared module(s) | Source | Depends on | Status |
|---|---|---|---|---|
| S1 material and rheology | `shared/material.py`, `shared/rheology.py` | skyscrapers sub-plan 02; the Step 4 plan's slurry law | stage 0 | to be written |
| S2 dispersion table | `shared/dispersion.py` | skyscrapers sub-plan 04 | stage 0 | to be written |
| S3 surface flow | `shared/surface_flow.py`, additions to `shared/gas.py` | skyscrapers sub-plan 05 | stage 0 | to be written |
| S4 spray physics | `shared/spray.py` (modes, `SprayModel`, `SOURCE_COLUMNS`, `source_rows`, `histogram`) | skyscrapers sub-plan 07, core only | S1–S3 | to be written |
| S5 Girin's published cases | `shared/girin_case.py`, `analysis/girin_reference.py` | skyscrapers sub-plan 08 | S2–S4 | to be written |
| S6 droplet outputs | `shared/droplets.py` (particle files, size distributions, plots) | skyscrapers sub-plans 10 (`write_particles`) and 12 | S4 | to be written |

Each is prototyped by moving the prototype module into a scratch copy of the package with `shared/` and running both
the prototype's tests and the package's; it is accepted when the module reproduces its prototype's outputs bit for bit.
