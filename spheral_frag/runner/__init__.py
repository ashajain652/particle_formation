"""spheral_frag.runner: the Spheral side of spheral_frag (spec §4.3; Spheral M2 plan).

Runs under Spheral's Python only (numpy 1.26, scipy 1.13, mpi4py 4.0), through the launcher
`spheral_frag/container/spheral`. Its modules import only the standard library, numpy, scipy, mpi4py, Spheral's
modules and the core's numpy-only modules; the core never imports this package (tests/test_spheral_frag_imports.py).
Spheral itself is imported inside the functions that need it, so that the modules can be inspected without it.

One form per process: importing SpheralRZ switches Spheral's global coordinate registrar to RZ, so a process builds
either 3D or RZ bodies, never both."""
