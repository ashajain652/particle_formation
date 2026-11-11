# Sub-plan 20 — working design notes (planning copy, 2026-10-12)

Base: `prototype/work-2026-10-12-moving-mesh/code` at 70ecbdd (sub-plan 19 built, sub-plans 11 and 12's reference
data and loader). Planning copy: `code/` here, tracked in the outer repository; the generator `make_plan20.py` cuts the
plan's blocks from it and checks they rebuild it from the base exactly.

## Rulings (each to be stated in the plan with its reason)

1. **Stage 1's skin is renamed, not deleted.** `reentry_model/skin.py` -> `skin1d.py` (and its tests); MeltingBody's
   stage-1 path is selected by `surface_model="skin1d"` and keeps its behaviour; sub-plan 22 removes it. The new
   `skin.py` is the layered skin (spec 10.3). `--surface-model skin` builds the new `SkinBody`. The default stays
   `elements` until sub-plan 22: stage 2 has neither flow nor spraying, so a skin run is a verification device.
2. **Columns.** Dual cells are barycentric (a third of each touching triangle, spec 5.1). Sideways conduction uses the
   cotangent weight of each surface edge, w = (cot a + cot b)/2, which is what the flux of a piecewise-linear field
   through the dual cell's two segments gives (exact for laterally linear fields); clipped at zero for obtuse pairs.
   The surface gradient and curvature operators are deferred to sub-plan 21, which is the first to use them.
3. **Cell targets** depend on the depth z of the cell's top below the column's top: dz_f for z < Z_f (the fine zone,
   the minimum thickness), dz_f + (g - 1)(z - Z_f) below (the continuous form of 10 % growth per cell). Split above
   1.5 x target into two equal halves; merge below 0.5 x target into the cell beneath (the bottom cell into the one
   above). Repeated until no cell qualifies.
4. **Draw happens inside the sub-steps** (a skin below its minimum after a sub-step's removal draws back to it at the
   base temperature), so a long macro step never empties a column; absorb and hand-down act after the sub-steps.
   Draw + absorb - hand-down become one inward displacement per corner, plus the volume-neutral smoothing; the move's
   actual swept volume and enthalpy are then reconciled on the skin's bottom cell, so mass and energy are exact.
5. **Liquid removed as it forms** (stage 2 only): after each sub-step the contiguous cells from the top at or above
   T_feed leave with their own enthalpy.
6. **3D model under skins**: no convective load and no radiation on its boundary (emissivity 0); the base law is a
   nodal load a_i + b_i T_i [W] at each corner node, implicit, in both backends.
7. **Heating and radiation** use the top surface (the base corners moved out by the skin thickness): angle from its
   area-weighted normals, power over its dual areas.
8. **SkinBody** is a ThermalBody with a corner surface (`columns.CornerSurface`: n_patches = corners), so CoupledRun's
   loop works unchanged; skin history columns through `skin_stats()`.

## State at the pause for Asha's decisions (2026-10-12)

Written and tested in `code/` (not yet assembled into the plan document; no generator yet):
Task 0 `movemesh.neutral_smoothing` (+ own-sweep `swept_volumes`, pending approval); Task 1 `columns.py`; Tasks 2-5 and 7
`skin.py` (cells, minimum, heat through depth with removal and draw, sideways, hand-down, PropertyTable); Task 6 nodal
law in both backends; Task 7 `skinbody.py`; Task 8 CLI (`--surface-model skin`, flags, run names), coupled-run columns,
results, frames; stage 1 renamed `skin1d` (value "skin1d", Python only). Harnesses: `analysis/ablation_1d.py`,
`analysis/skin_box.py`, `../measure/timing_100mm.py`. Reference data: `data/reference_values/box_q6p_fine.json`.

Rulings made while writing (to state in the plan): draw at min(T_b, T_rigid) (3D node overshoot to 930 K on 4 mm
elements made a liquid draw-remove loop); property tables at 0.01 K with books on stored enthalpy; regrid only rows that
need it; Newton tangent = table slope; gradient/curvature deferred to 21; velocity deferred to 21; SkinBody in its own
module; CLI default stays elements; skin + melt needs --removal instant; --animate refused with skins.

Findings needing Asha:
1. Third-of-prism booking (spec 5.2) credits non-moving neighbours with a moving corner's sweep; at the box's edges the
   side walls bulged out 0.33 mm and folded the surface (fail at 1.3 s). Own-sweep booking (A_t u_k.n_t / 3 per corner,
   second-order rest shared in proportion) fixes the bulge, keeps all 33 movemesh tests, and makes the neutral
   smoothing converge ~20x per pass (worst corner 0.46 um plain, 3e-11 m after 3 passes, 2 mm sphere).
2. The level-2 box cannot reach 2 s: the 90-degree edges hit the fold limit (2 mm mesh: 1.4 s; 1 mm: 0.9 s). Before
   failing, removed-mass recession on 1 mm matched the 1D reference at 0.75 s (0.101 vs 0.101 mm); centre lags.
   1D reference (1 um, 0.5 ms): 0.789, 2.712, 5.219 mm at 1.0, 1.5, 2.0 s; fine box 0.8045, 2.729, 5.206.
3. Run time, 100 mm default sphere, production: 51-57 s per 0.5 s macro step solid (375k cells), 60-115 s melting
   (450k cells); 3D under 1 s; interface repeat every step (doubles cost).

## After Asha's decisions (2026-10-12)

Superseded by `part_20.md` (the plan's prose) and `make_plan20.py` (the generator). Decisions: own-sweep booking, the
plate for level 2 (remeshing the back-up), the repeat-trigger measurement. Further ruling under the spec's fallback: the
rebuild's L2 projection. The predictor was tried and dropped. Facts 108-114 in 00-shared-context.md.
