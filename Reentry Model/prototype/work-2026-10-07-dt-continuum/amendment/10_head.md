## Amendment of 2026-10-07 (continuum step) — the macro step shortens where Girin's closure begins

> Asha's decision of 2026-10-07 on fact 96 (a): make the switched step a model option at 0.0125 s (facts 97–100 in
> `00-shared-context.md`). The code below is the tested code, as a diff against the copy of this date's earlier
> amendment (the Scheil default and the freeze-back fix). Sub-plan 13's amendment of this date holds the flag.

**The change.** `CoupledSettings` gains `dt_continuum` (s; `None`, the default, for no switch) and `kn_switch` (the body
Knudsen number of the switch, default `surface_flow.KN_BODY_SHOCK`, 0.01). After each macro step, on the aero state the
body was just given — before any re-evaluation on a surface the step's deaths changed — `CoupledRun._switch` sets the
macro step to `dt_continuum` the first time the trajectory's body Knudsen number (`AeroState.kn`, the mean free path
over the body's current reference length, the history's `knudsen` column) is below `kn_switch`, and latches it. That is
the rule of the 2026-10-05 harness (`dtswitch.py`), which applied it from outside the package, so the option reproduces
the harness's runs bit for bit (fact 98). `dt_continuum` must be greater than zero and no longer than `dt`.

**Frames.** `frames_every` counts macro steps until the switch and flight time after it — a frame every `frames_every`
first steps — so a fine step does not write a frame per step; that is also the harness's rule, which wrote frames every
10 s of flight time.

**Results.** With `dt_continuum` set the run's results gain `dt_continuum_s`, `dt_switch_time_s`, `dt_switch_kn` and
`dt_switch_altitude_km` (the three `None` if the flight never crosses). History rows stay one per macro step, so a run
with the switch has many more of them after it (5 740 rows to 120 s at 0.0125 s on the 100 mm flight, against 241).

**Tests.** `test_the_continuum_step_reproduces_the_switched_step_harness`: a hot body on the coarse mesh from 48 s, at
0.5 s to the switch (48.5 s there: the body's Knudsen number reads its shrinking size) and 0.05 s to 50.5 s, once with
the option and once with the harness's rule applied by hand: every history column identical, the switch recorded where
the harness switched, the steps 0.5 s before and 0.05 s after it, and frames at 48.0, 48.5, 49.0, 49.5, 50.0 and 50.5 s
— by step, then by flight time. `test_the_continuum_step_is_off_by_default_and_rejects_a_coarser_step`. Both failed
before the change (no such setting) and pass after it.

