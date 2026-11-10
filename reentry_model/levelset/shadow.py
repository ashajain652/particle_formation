"""The ray-cast shadow test of design §20.8: from each patch, offset a twentieth of its cell outward, a ray is marched
upstream (along the direction of motion) through phi by sphere tracing, with steps of phi but at least a quarter of the
smallest cell; a ray that re-enters metal before leaving the box marks the patch as shadowed. (A whole-cell offset was
tried first: near the outline of a sphere it carried leeward rays outside the body, so they missed it.)"""
import numpy as np

from .transfer import locate


def shadowed(mesh, phi, P, direction, max_steps=400):
    d = np.asarray(direction, float)
    d = d / np.linalg.norm(d)
    cell, _ = locate(mesh, P.centroid)
    size = mesh.size[cell]
    x = P.centroid + P.normal * (0.05 * size)[:, None]
    min_step = 0.25 * mesh.size.min()
    lo, hi = np.asarray(mesh.base.lo), np.asarray(mesh.base.hi)
    hit = np.zeros(len(x), bool)
    active = np.ones(len(x), bool)
    for _ in range(max_steps):
        idx = np.flatnonzero(active)
        if len(idx) == 0:
            break
        c, b = locate(mesh, x[idx])
        f = np.einsum("ni,ni->n", b, phi[mesh.tets[c]])
        inside = f < 0
        hit[idx[inside]] = True
        active[idx[inside]] = False
        go = idx[~inside]
        x[go] += np.maximum(f[~inside], min_step)[:, None] * d
        out = ((x[go] < lo) | (x[go] > hi)).any(axis=1)
        active[go[out]] = False
    else:
        raise RuntimeError("shadow rays did not finish in {} steps".format(max_steps))
    return hit
