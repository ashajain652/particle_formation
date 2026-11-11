# Test fixtures

Real SESAM 2.3.0 (DRAMA 4.1.4) outputs captured on 2026-09-13 during the
design spike, one directory per case. Each holds the four files the parsers
read: `*_AeroThermalHistory.txt`, `*_Trajectory.txt`,
`PySara.ImpactingFragments.xml`, `sesam.log`.

| Case | Sphere | Start state (alt km, lat, lon, v km/s, fpa, heading, epoch) | energyThreshold | Ends |
|---|---|---|---|---|
| `T1_demised_50mm` | 50 mm, 300 K | 101.247, -3.749, -74.640, 7.907, -0.190, 347.960, 2024-08-01T12:44:58 | 15 J | melts to 0 kg, `uncritical` at 80.14 km |
| `T3_survivor_50mm` | 50 mm, 300 K | 39.940, 35.911, -83.878, 0.500, -32.790, 347.122, 2024-08-01T12:55:51 | 15 J | `ground impact`, XML mass 1.8411041946975187e-01 kg |
| `T5_5mm_750K` | 5 mm, 750 K | 77.500, 29.546, -82.134, 7.500, -0.959, 347.168, 2024-08-01T12:53:07 | 15 J | `uncritical` after 8.3 s, mass column all 0.000 |
| `E1_ballooning_5mm` | 5 mm, 750 K | same as T5 | 1e-6 J | `ballooning` after 8.5 s |

To regenerate a case once `sphere_reentry.py` exists: run it with the state
above and `--keep-raw`, then copy the four files from
`sphere_sweep_output/raw/<run_name>/run_0/reentry/`. Expected values in the
tests were taken from these exact files and may need refreshing afterwards.
