# Re-frozen metal stays in the skins: design draft

**Status: draft, 2026-11-08. Nothing built.** Three changes, to be built and flown together: the part-count fix, a switch
that keeps re-frozen metal in the skins, and sideways heat between cells that touch.

## Why

The spreading flights of 2026-11-07 (limit 60 degrees, resumed from the states saved at 170.5 s and 22.0 s) passed the
old stops, 179.04 s on the 50 mm sphere and 33.0 s on the 100 mm sphere, with the largest fold 68.0 and 86.5 degrees
against 144 and 138. Both then stalled and were stopped by Asha at 16:55 on 2026-11-08 (179.66 s and 36.18 s): the
parts of each 12.5 ms macro step grew by a quarter every step, from 9 to 9 765 and from 1 to 1 537.

Every fold that stopped a flight sat at the outer edge of the ring of re-frozen metal, where the hand-down pushes the
base outward; the receded nose had none. Asha (2026-11-08): the re-frozen metal need not pile up, because it can reheat,
flow on to other skins or spray, so it may belong in the skins. The earlier flights agree that it melts again: on the
100 mm sphere 64.1 g was handed down over the flight but about 26 g stood above the starting surface at the stop, on the
50 mm sphere 27.7 g against 7.8 g. Today that metal re-melts only by a round trip, handed down and then absorbed again.

## 1. The part count (fix)

`SkinBody.advance` sets the next macro step's parts to `ceil(MOVE_HINT_MARGIN * max(ratio * dt_macro / dt))`, where the
ratio is a part's largest move over the allowance (a quarter of the mean surface edge). Scaling by the number of parts
assumes a move shrinks with the part. Drawing does; a hand-down capped at the allowance (2026-10-27, because a held cell
released whole does not shrink with the part) does not, so one capped hand-down in a step makes the next 25 % finer,
without end. The 50 mm sequence 837, 1 047, 1 309, ... 9 765 is that rule applied exactly (inferred, not instrumented).

**Change:** only the draws set the next step's count. The split test is unchanged; since a capped hand-down is always just
under the allowance, only a draw can split a part. Cost: held metal drains at one allowance a part (at most 1.5 mm and
1.1 mm deep on the two flights, about three parts).

**Guard (Claude's proposal; Asha may reject):** a flight stops with a clear message when a macro step needs more than a
set number of parts, provisionally 1 000.

## 2. Keeping re-frozen metal in the skins (switch)

A setting turns the hand-down off; the default keeps today's behaviour. With it off, a skin returns toward its minimum
only by re-melting and flowing on, or by spray. The base then moves only inward (draws and absorption), so the deposit no
longer folds the 3D mesh, and spreading and the capped hand-down have nothing to act on.

Expected costs, untested: deeper columns (about 90 cells for a 10 mm column against 59 at the minimum, estimated from
§6.2's 10 % grading); the fold moves to the top of the metal (tops reached 178 and 179 degrees in the earlier flights, and
the flow's push and Girin's spraying read the top's slope); straight columns misplace volume on a curved body as they
deepen (the 50 mm pools are already 22-45 mm deep on a 25 mm radius). A middle option, not chosen: keep metal in the
skin up to a set depth and hand down below it.

## 3. Sideways heat between cells that touch (Asha, 2026-11-08)

Spec §6.3 pairs each cell with the neighbour's cells at the same depth below the top. Where the tops differ in height
this pairs cells that do not touch: at the ring's outer edge the top steps about 1.4 mm across a 1.5 mm triangle, so the
cell 1 mm below the taller top, which is beside gas, is paired with metal 1.4 mm lower. While the deposit is handed to
the 3D body the mesh carries this correctly; kept in the skins, it would not.

**Change:** depth ranges are measured upward from the base, the continuous surface of the 3D body beneath the skins,
plus the offset between the two columns' base positions along their direction, and heat is shared by the overlap of
those ranges times the shared face's length, as now. Asha's alternative, tracking every cell's Cartesian position, gives
the same answer at 2 mm spacing at more cost and is not chosen (Claude's call).

What it brings out, stated: the part of a taller column above its neighbour's top is a side wall facing the gas; it
correctly exchanges no heat with the neighbour, and the model does not heat or cool it from the gas.

Not changed: the flow's sideways transport, which matches cells by fraction of the flowing depth and concerns only the
liquid.

## Tests

The part count holds steady on a small body that hands down a capped cell every step (and grows by a quarter today).
Mass, enthalpy and momentum exact to round-off with the hand-down off. Heat exchange between two columns of different top
heights pairs cells by height above the base, and equals today's where tops are level. Defaults reproduce existing runs
to the last bit.

## The check

Resume both spheres from the states saved at 170.5 s and 22.0 s with the hand-down off, the part-count fix and
height-matched sideways heat, beside the spreading flights: do they pass the old stops, how deep the skins get, how far
the tops fold, how much is sprayed, and whether the books close. Run on mains power.

## Open

The guard's limit; whether the middle option is needed if columns grow too deep; whether `--slurry-spray on` flights
are added.
