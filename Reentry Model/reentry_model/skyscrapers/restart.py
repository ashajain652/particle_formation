"""A flight's whole state on disk, and back (sub-plan 22 Task 5a; Asha 2026-11-01).

The convergence repeats of Task 5 -- a halved macro step, another seed, 20 instead of 10 um skin cells -- are compared
over a window of about 10 s around peak heating, 93.5 s into the 100 mm flight; flown from the start each would cost
nearly a whole flight. So the base flight saves its state once on its way (`coupled.CoupledSettings.save_at`, the
command line's --save-state-at) and each repeat resumes from it (`coupled.CoupledRun.run(resume=...)`).

The state is the coupled run itself -- the trajectory, the body with its 3D mesh, heat solver, skins, moving mesh and
spray model -- with numpy's global random state, whose only draws in a run are pyamg's (cli.cmd_run), pickled whole.
Measured on the 100 mm sphere's production mesh (2026-11-01): 1.32 GB, written in 0.4 s, and everything in it pickles but
pyamg's coarse-grid solver, a class made inside `pyamg.multilevel.coarse_grid_solver`, which pickle cannot name. It is
saved as the name of its solver and the factor it has cached (`P`, `LU` or `L`), and comes back as a new solver of that
name holding that factor, so a hierarchy built before the save solves after it exactly as it would have: a flight
resumed from its state in a fresh process repeated the unbroken one to the last bit. The objects the flight holds are
never altered by a save."""
import ast
import io
import os
import pickle

FORMAT = 1                      # of a flight's saved state (`coupled.CoupledRun.state`); a resume refuses any other
PART = ".part"                  # the state is written beside its path under this suffix, then renamed into place
COARSE_SOLVERS = ("pinv", "pinv2", "lu", "cholesky")    # pyamg's coarse solvers whose caches are plain arrays


def _is_coarse_solver(obj):
    cls = type(obj)
    return cls.__module__ == "pyamg.multilevel" and cls.__qualname__.endswith("coarse_grid_solver.<locals>.GenericSolver")


def _coarse_name(obj):
    """The solver a pyamg coarse-grid solver was made with, which must be one of COARSE_SOLVERS given no arguments: the
    arguments live in the closure that made it, where a save cannot carry them, so a solver with any is refused."""
    try:
        name = ast.literal_eval(type(obj).name())
        cells = dict(zip(type(obj).__call__.__code__.co_freevars, type(obj).__call__.__closure__ or ()))
        solve = cells["solve"].cell_contents
        kwargs = dict(zip(solve.__code__.co_freevars, solve.__closure__ or ())).get("kwargs")
        kwargs = {} if kwargs is None else kwargs.cell_contents
    except (KeyError, ValueError, SyntaxError, AttributeError, TypeError):
        raise TypeError("a pyamg coarse solver whose make cannot be read cannot be saved: {!r}".format(obj))
    if name not in COARSE_SOLVERS or kwargs:
        raise TypeError("a pyamg coarse solver {!r} with arguments {!r} cannot be saved (saved: {} with none)".format(
            name, kwargs, ", ".join(COARSE_SOLVERS)))
    return name


def _coarse_solver(name, cache):
    """A pyamg coarse-grid solver of `name` holding `cache`, the attributes the saved one had cached."""
    from pyamg.multilevel import coarse_grid_solver
    solver = coarse_grid_solver(name)
    vars(solver).update(cache)
    return solver


class _Pickler(pickle.Pickler):
    def reducer_override(self, obj):
        if _is_coarse_solver(obj):
            return _coarse_solver, (_coarse_name(obj), dict(vars(obj)))
        return NotImplemented


def dumps(state):
    """`state` pickled (protocol 5), the coarse solvers as above, as bytes. These are the pickling buffer itself, not a
    copy: `BytesIO.getvalue` hands over the buffer it wrote into, trimmed to the pickle's length, when nothing else holds
    it, as here. So a save or a copy of a completed step holds the pickle once (final fix wave of the geometry stop,
    2026-11-12, measured with tracemalloc: at its peak 1.125 times the pickle, the buffer's growth margin, and the pickle
    alone once returned; peak resident memory did not rise at the `getvalue` of an 800 MB pickle). `getbuffer` would keep
    the growth margin, up to an eighth of the pickle, for as long as the bytes are held."""
    buf = io.BytesIO()
    _Pickler(buf, protocol=5).dump(state)
    return buf.getvalue()


def loads(data):
    """The state `dumps` made `data` of, as a new object graph."""
    return pickle.loads(data)


def write_bytes(path, data):
    """Write `data` to `path` through path + PART and a rename, so a write cut short leaves no file at `path` (and its
    part removed when the failure is Python's)."""
    part = path + PART
    try:
        with open(part, "wb") as fh:
            fh.write(data)
        os.replace(part, path)
    finally:
        if os.path.exists(part):
            os.remove(part)


def save(path, state):
    """Write `state` to `path`: `dumps` of it, by `write_bytes`."""
    write_bytes(path, dumps(state))


def load(path):
    """The state `save` wrote at `path`."""
    with open(path, "rb") as fh:
        return pickle.load(fh)
