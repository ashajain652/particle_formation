| case | mode | surface | Q_conv max (of peak) | Q_conv max (continuum, point-wise) | integrated heat (end of hypersonic) | integrated heat (end) | max |dT_eq| | radiated max (of peak) | peak surface T | runtime | verdict | skins |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| d100 | sesam | elements | 0.46% | 2.23% | +0.31% | +0.30% (ratio 1.003) | 24.6 K (1.20%) | 6.66% | 2089 K at 139 s | 32 s, 733 steps, 2.0 it/step | pass |  |
| d100 | physics | elements | 25.32% | 30.35% | -25.60% | -25.63% (ratio 0.744) | 410.0 K (20.57%) | 52.14% | 2319 K at 118 s | 39 s, 733 steps, 2.1 it/step | reported |  |
| d100 | physics | skin | 25.26% | 28.49% | -25.49% | -25.52% (ratio 0.745) | 410.0 K (20.61%) | 51.94% | 2322 K at 118 s | 78815 s, 733 steps, 2.4 it/step | reported | 9417 columns, 0 rebuilds, 335 repeats, 0 base faults; 66832 s in the skins, 1498 s in the 3D solve, 107.5 s a step |
| d050 | sesam | elements | 0.37% | 2.36% | +0.03% | +0.02% (ratio 1.000) | 28.0 K (1.22%) | 3.40% | 2442 K at 245 s | 8 s, 1118 steps, 2.0 it/step | pass |  |
| d050 | sesam | skin | 0.51% | 0.83% | +0.15% | +0.14% (ratio 1.001) | 31.1 K (1.38%) | 3.07% | 2443 K at 245 s | 56749 s, 1118 steps, 2.1 it/step | pass | 2367 columns, 0 rebuilds, 792 repeats, 0 base faults; 52138 s in the skins, 612 s in the 3D solve, 50.8 s a step |
| d050 | physics | elements | 25.96% | 35.42% | -22.94% | -22.98% (ratio 0.770) | 351.5 K (15.96%) | 40.60% | 2494 K at 239 s | 12 s, 1118 steps, 2.0 it/step | reported |  |
| d050 | physics | skin | 25.82% | 33.69% | -22.83% | -22.86% (ratio 0.771) | 353.4 K (15.73%) | 40.44% | 2497 K at 240 s | 57206 s, 1118 steps, 2.2 it/step | reported | 2367 columns, 0 rebuilds, 835 repeats, 0 base faults; 52492 s in the skins, 648 s in the 3D solve, 51.2 s a step |

physics-mode power errors are differences from SESAM's surface-averaged model, not thresholds
