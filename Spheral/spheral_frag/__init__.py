"""spheral_frag: large fragments of the re-entering sphere with Spheral (spec 2026-10-02-spheral-large-fragments-design).

The core (every module at this level) imports only the standard library and numpy at module level, so that the
runner can import it under Spheral's Python; functions that need scipy, pyvista, vtk or gmsh import them inside the
function ("prepare side only"). Only `fe.py` touches the finite-element package, and only inside its functions.
`runner/` (later milestones), `m0/` and `container/` need Spheral and run under its Python only."""

__version__ = "0.1.0"
