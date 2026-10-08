## Amendment of 2026-10-07 (continuum step) — the drivers switch by default, and the harness retires

> Asha's decision of 2026-10-07 on fact 96 (a) (facts 97–100 in `00-shared-context.md`). The code below is the tested
> code, as a diff against the copy of this date's earlier amendment. These items stack on this sub-plan's amendments
> above.

1. **Every physics and sensitivity run now switches.** `melt_verification.py`'s `physics` mode and every
   `melt_sensitivity.py` variant run with `--removal girin`, so from the first step under the continuum boundary they
   take 0.0125 s steps (sub-plan 13's default). The bookkeeping and resolved modes run `--removal instant` and Girin
   removal with `--size-feedback initial` respectively: the first keeps the default step throughout (the flag's default
   is off for instant removal), the second switches like the physics mode. The 50 mm flight never crosses the continuum
   boundary, so it is unchanged bit for bit (fact 99); the 100 mm flight switches at 49.5 s and to 120 s costs about twelve times
   as much (fact 99 gives the measured times).
2. **A `dtcoff` variant** (`--dt-continuum off`): the default 0.5 s step throughout, as every row before this amendment
   was run, so the table carries the step's effect against `base` like every other switch's. `dt025` now halves the step
   before the switch only.
3. **The switched-step harness retires.** `dtswitch.py` and `compare.py` (2026-10-05) wrapped the model from outside the
   package to switch the step; `--dt-continuum <s>` does the same from the model's own command line, and reproduces the
   harness's 0.0125 s and 0.5 s runs of the 2026-10-07 series bit for bit (fact 99). A step series is now a set of
   ordinary runs with `--dt-continuum 0.05`, `0.025`, `0.0125`, `0.00625` and `off`; the harness's runoff-by-latitude
   instrumentation stays in the throwaway copy, as it never belonged in `analysis/`.

