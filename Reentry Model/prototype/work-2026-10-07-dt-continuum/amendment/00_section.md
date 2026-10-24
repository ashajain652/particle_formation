## Amendment of 2026-10-07 (continuum step) — the macro step shortens where Girin's closure begins (facts 97–100)

Asha's decision of 2026-10-07 on fact 96 (a): the switched step becomes a model option, `--dt-continuum`, at 0.0125 s,
option (2). It amends sub-plans 10 (the switch in the coupled loop), 13 (the flag), 14 (the drivers) and 15 (the README,
the assumptions and the spec). Facts 1–96 stand except where these say otherwise; **fact 96 (a) is answered, and fact 77
(b)'s question — should the switch become a model option — with it.** Measured in a throwaway copy
(`prototype/work-2026-10-07-dt-continuum/`, git-ignored) made from the copy of facts 88–96, which that amendment's seven
diff blocks rebuild exactly from the earlier copy; the measurement copy was frozen with a checksum manifest before the
runs and never edited while one was in flight; every run was made on mains power under `caffeinate -i`, the power source
logged at the start and the end.

97. **The option.** `CoupledSettings.dt_continuum` and `kn_switch`, and `--dt-continuum off|<s>`. After each macro step,
    on the aero state the body was given, the step becomes `dt_continuum` the first time the trajectory's body Knudsen
    number — the mean free path over the body's current reference length, the history's `knudsen` column — is below
    `--kn-body-shock` (0.01, the surface flow's own continuum boundary, where Girin's closure begins), and stays there.
    That is the harness's rule (fact 69). Its default depends on the removal mode, as `--size-feedback`'s does: 0.0125 s
    for `--removal girin`, the model proper, where the series converges (facts 92–94); off for `--removal instant`, the
    bookkeeping device, which has no film and whose SESAM thresholds were measured at the default step. Frames count
    steps until the switch and flight time after it. The results record the value and the switch's time, Knudsen number
    and altitude. Only a non-default value changes the run name (`_dtcontinuum-off`, `_dtcontinuum-<s>`); a value of
    zero or below, one longer than `--dt`, or the flag without `--melt on` exits 2.

98. **Tests.** A coupled test runs a hot body on the coarse mesh from 48 s through the switch twice — with the option
    and with the harness's rule applied by hand — and requires every history column to be identical, the switch recorded
    where the harness switched (48.5 s there, the body's Knudsen number reading its shrinking size), the steps 0.5 s
    before it and 0.05 s after, and frames by step then by flight time. The CLI tests check the default per removal
    mode, `off`, the run names and three bad values. `test_melting_run_writes_columns_files_and_json` crosses the
    continuum boundary at about 9.5 s of its 15 s from 71 km, so with the default it ran 460 macro steps instead of 30;
    it now passes `--dt-continuum 0.1` and checks the switch it records. Unit tier 258 passed, 1 skipped, with only Task
    11's five known reference failures; FEniCSx tier 12 passed. Writing the coupled test showed one thing that is not
    the option's: where fine steps sum to just short of `t_max`, the trajectory stepper ends the run with a sliver step
    (1.1e-13 s in the test; the series' 0.05 and 0.00625 s runs end the same way), which adds one history row and
    changes nothing measurable (fact 100 (b)).

99. **The option reproduces the harness's runs bit for bit.** Run with the defaults — `AA7075_scheil`, the rule on,
    `--dt-continuum` at its default — the 100 mm flight to 120 s switches at 49.5 s, at 69.93 km and a body Knudsen
    number of 0.009905, and is bit-identical to the series' 0.0125 s harness run (fact 91): all 5 740 rows of the 90
    history columns, every result field the harness run has but the run time and the frame count, all 22 source-table
    columns, and its 13 frames at the same flight times. Its energy balance closes to 3.1e-9. With `--dt-continuum off`
    the same flight is bit-identical to the series' 0.5 s run — all 241 rows of the 90 history columns, every result
    field the older run has but the run time and the frame count, all 22 source-table columns, and its 13 frames at the
    same flight times. The 50 mm flight on the linear range never crosses the continuum boundary — it records no switch
    — and is bit-identical to the run of fact 89 with the option at its default. So the harness can retire: a step
    series is a set of ordinary runs. The cost: the default run took 8 207 s (2.3 h) of model time for 5 739 macro
    steps, 5 640 of them fine, about 1.4 s each on average and 1.9 s early in the window, while the body is large
    (measured from its frames: 10 s of flight in 26 min between 50 and 60 s); the run with the switch off took 643 s
    (10.7 min), but it shared the machine with two other runs for most of its time, so the ratio between them, 12.8, is
    approximate.

100. **Not done, and for Asha to decide.** (a) **The cost of a whole flight — decide this first.** A whole 100 mm flight
     now switches at 49.5 s and keeps the fine step to the ground. The committed non-melting references reach the ground
     at 366–368 s, so a whole flight is about 25 000 fine steps; at the 1 to 1.9 s a fine step costs (fact 99; the
     larger figure early in the window, while the body is large) that is about 7 to 13 hours, against about 40 minutes
     at 0.5 s (fact 23), and every 100 mm row of Task 14's tables costs it. The spray has nearly stopped by 120 s — at
     0.0125 s it releases about 20 g/s from 60 to 100 s, 7.5 g/s over 100–110 s and 1.5 g/s over 110–120 s, the hottest
     surface at 910 K, the top of the feed ramp — so most of those fine steps fall where nothing sprays. Options: (1)
     keep it and budget the runs; (2) return to the default step once the melting is over — for example from the first
     step after the switch at which nothing is fed or sprayed and no surface node is above the half-liquid temperature —
     as a change of its own, with a measurement that the remnant's state at the ground does not depend on the step; (3)
     end Task 14's physics runs at 120 s, which loses the survivor's state at the ground. Recommendation: (2), designed
     and measured on its own; until then the step-sensitive droplet quantities, which are settled by 120 s, can be
     quoted from runs to 120 s. (b) **The sliver last step** (fact 98): snap the trajectory stepper's last step to
     `t_max` when what remains is below a small fraction of the step; it would change the last history row of every run
     that ends this way, so it belongs in a change of its own with the runs re-measured. (c) Facts 96 (b) — the thick
     branch's drift at the finest step — and 96 (c) — the run name does not carry the material — stand.

