"""The runner's own checkpoints (Spheral M2 plan, decision D3; Asha accepted 2026-10-09).

At a frame boundary the runner writes `checkpoint_<k>.npz`: every particle's state gathered to rank 0 in pid order,
so that the file is independent of the domain decomposition, and the run resumes or branches from it (a window
converts the state anyway, D2). Spheral's restart files are not used.

Format: a plain `.npz` without pickles (numpy 1.26 writes it under Spheral's Python, numpy 2 reads it in drama_env,
M1 fact 10). Per-particle arrays, named as `nodes.STATE_KEYS` (+ `nodes.RZ_KEYS` in RZ):

  pid (n,) int64; x_ref, x, v (n, dim); H, S (n, dim, dim) (S the deviatoric stress); m, rho, u, plastic_strain,
  plastic_strain_rate (n,); RZ only: m_rz, rho_rz (n,) (SolidSPHRZ's massRZ and massDensityRZ state) and S_zz, H_zz
  (n,), the zz elements of the RZ tensors (the deviatoric stress's hoop component; schema version 2, Task 5: version 1
  dropped it, so an RZ restart zeroed the hoop deviatoric stress)

Metadata as 0-d arrays: schema "spheral_frag.checkpoint", schema_version, run (the run name), k (the frame), form,
dx, nperh, scale (the mass-scaling factor s in force), processes (that wrote it). Driver quantities (time, last time
step, step count, ...) go in `extra`, stored as "extra.<name>" and returned by `read` under "extra".

`read` needs only numpy; `write` needs Spheral (through the body) and MPI and is collective."""
from __future__ import annotations

import os

import numpy as np

from . import nodes as _nodes

SCHEMA = "spheral_frag.checkpoint"
SCHEMA_VERSION = 2
META_KEYS = ("schema", "schema_version", "run", "k", "form", "dx", "nperh", "scale", "processes")
EXTRA_PREFIX = "extra."


def filename(k):
    return f"checkpoint_{int(k):05d}.npz"


def write(path, body, extra=None, *, run, k, scale=1.0):
    """Gather the body's state to rank 0 in pid order and write it to `path` (atomically: a temporary file then a
    rename). Collective; returns the path on rank 0 and None elsewhere. `extra` maps names to numbers or arrays."""
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    state = body.gather()
    if comm.rank != 0:
        comm.barrier()
        return None
    write_state(path, state, form=body.form, dx=body.dx, nperh=body.nperh, run=run, k=k, scale=scale,
                processes=comm.size, extra=extra)
    comm.barrier()
    return path


def write_state(path, state, *, form, dx, nperh, run, k, scale, processes, extra=None):
    """Write a gathered state (per-particle arrays named as nodes.STATE_KEYS + RZ_KEYS, sorted by pid) as a
    checkpoint, atomically. Not collective: the caller (rank 0) holds the state, e.g. a gather with particles
    removed at a frame boundary (runner/run.py)."""
    keys = _nodes.STATE_KEYS + (_nodes.RZ_KEYS if form == "rz" else ())
    missing = [key for key in keys if key not in state]
    if missing:
        raise ValueError(f"state lacks {missing}")
    if np.any(np.diff(np.asarray(state["pid"])) <= 0):
        raise ValueError("state rows must be sorted by unique pid")
    out = {key: np.asarray(state[key]) for key in keys}
    out.update(schema=np.array(SCHEMA), schema_version=np.array(SCHEMA_VERSION, dtype=np.int64),
               run=np.array(str(run)), k=np.array(int(k), dtype=np.int64), form=np.array(form),
               dx=np.array(float(dx)), nperh=np.array(float(nperh)), scale=np.array(float(scale)),
               processes=np.array(int(processes), dtype=np.int64))
    for name, val in (extra or {}).items():
        arr = np.asarray(val)
        if arr.dtype == object:
            raise TypeError(f"extra {name!r} is not numeric")
        out[EXTRA_PREFIX + name] = arr
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        np.savez(f, **out)
    os.replace(tmp, path)
    return path


def read(path):
    """The checkpoint as a dict: the per-particle arrays, the metadata as Python scalars and `extra` as a dict
    (0-d entries as Python scalars). Plain numpy; refuses another schema or a newer version."""
    with np.load(path, allow_pickle=False) as z:
        raw = {k: z[k] for k in z.files}
    if str(raw.get("schema", "")) != SCHEMA:
        raise ValueError(f"{path}: not a {SCHEMA} file")
    if int(raw["schema_version"]) > SCHEMA_VERSION:
        raise ValueError(f"{path}: schema version {int(raw['schema_version'])} is newer than {SCHEMA_VERSION}")
    out, extra = {}, {}
    for k, v in raw.items():
        if k.startswith(EXTRA_PREFIX):
            extra[k[len(EXTRA_PREFIX):]] = v.item() if v.ndim == 0 else v
        elif k in META_KEYS:
            out[k] = v.item()
        else:
            out[k] = v
    out["extra"] = extra
    keys = _nodes.STATE_KEYS + (_nodes.RZ_KEYS if out["form"] == "rz" else ())
    missing = [k for k in keys if k not in out]
    if missing:
        raise ValueError(f"{path}: missing {missing}")
    return out


def arrays(ckpt):
    """The per-particle arrays of a read checkpoint, ready for nodes.make_body."""
    keys = _nodes.STATE_KEYS + (_nodes.RZ_KEYS if ckpt["form"] == "rz" else ())
    return {k: ckpt[k] for k in keys}
