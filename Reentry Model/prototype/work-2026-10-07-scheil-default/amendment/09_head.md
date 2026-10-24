## Amendment of 2026-10-07 — freeze-back never leaves a negative film

> Found by the time-step series of 2026-10-07 (facts 88–96 in `00-shared-context.md`). The code below is the tested
> code, as a diff against the copy of the rigid-substrate amendment of 2026-10-06 with sub-plan 13's Scheil default of
> this date applied.

**What failed.** The 0.00625 s run of the series stopped at its first fine step, at 49.50625 s, with the thermal
solver's input check "film mass must be one finite non-negative value per node" (exit code 2, fact 77 (c)). A probe that
checks the film and deep accounts after every stage of the melt step, re-running the same seeded flight, found the first
bad value right after `_freeze_back`: a film of −1.29e-25 kg on one windward patch, 51° from the stagnation point, whose
deep liquid freeze-back had just taken whole. Nothing was non-finite; the non-rigid depth, `film.lubrication`'s outputs
and the runoff transport were all clean up to that point.

**Why.** The deep-runoff amendment of 2026-10-02 caps what freeze-back takes from a patch at its film plus its deep
liquid, takes the deep liquid first and debits the film by the rest, computed as the capped amount less the deep part:
the film is left with m_f − ((m_f + m_d) − m_d). In floating point (m_f + m_d) − m_d is m_f rounded to the precision of
m_d, so when the deep liquid dwarfs the film, the film's part comes out one rounding unit above the film for about half
of all such pairs (measured on 200 000 random pairs: 99 860), and the film is left at about −1e-25 kg. The spray stage's
clean-up, which zeroes films below 1e-30 kg, runs before freeze-back, so the negative reaches the solver. It needs a
patch holding both film and deep liquid whose liquid all freezes back in one step. No earlier run stopped on it, and
none of the 2026-10-06 flights is changed by the fix: re-run on the fixed code, the 100 mm and 50 mm flights on the
linear range with the rule on and off and the 100 mm flight on Scheil's curve with the rule off are bit-identical to the
originals in every history column, result and source-table column (fact 89). The Scheil series met it at its first fine
step.

**The change.** Each account gives up at most what it holds, and the film's part is taken directly: `from_deep =
min(m_d, wanted)`, `from_film = min(m_f, wanted - from_deep)`, `taken = from_deep + from_film`. Each subtraction then
takes from a number no smaller than what it takes, which correctly rounded arithmetic never rounds below zero. Where the
cap does not bind, the amounts are those of before; where it binds they can differ in the last bit.

**Test.** `test_freeze_back_of_film_and_deep_liquid_leaves_neither_negative`: one patch with 2.34e-14 kg of film over
1.88e-9 kg of deep liquid — a pair the old arithmetic gets wrong, which the test asserts — all of it freezing back into
an owner element with room. Both accounts must be exactly zero afterwards and the owner must gain their mass. It failed
before the change with the film negative and passes after it. Unit tier 252 passed, with Task 11's five known reference
failures; FEniCSx tier 12 passed.

