"""The step's end of the level set's skinned body (design sections 20.3-20.6; level-set stage 3 Task 9, split from
`skinbody` in its fix round): the geometry that follows the skins -- the skin base beta moved by the skins' exchange with
the volume, the outer surface phi and the molten-rigid surface psi rebuilt, the mesh rebuilt where it no longer meets its
refinement targets, and the markers settled on phi -- as `StepEnd`, a mixin of `skinbody.LevelSetSkinBody`, whose
attributes it reads (mesh, beta, phi, psi, solver, markers, skins, material, base, settings, spacing, band,
fine_margin, fine_size, cell_floor, rho). The physics is the geometry modules' (`levelset`, `surfaces`, `criteria`,
`boxmesh`, `markers`) and the shared package's; every number here is a numerical setting.

Moving beta (`_move_beta`): volume changes dV at points, each a displacement along the normal, extended to the nodes
within fine_margin of beta from the nearest point; beta - delta reset to a distance and corrected per piece to the
exact volume target (`levelset.correct_volumes`, the targets summed piece by piece with pairwise sums); then
`CutThermal.move`, whose enthalpy against what the skins took and gave the caller books as pending heat. A piece other
than the largest whose target is not positive -- the draws of the markers standing on it ask more than it holds -- is
used up: its nodes are lifted just outside the metal (PIECE_GONE), so that its old metal and the draws at its points
fall to the largest piece's target, and the pieces are found again (fix round 5 of Task 9: on the full-size cube's third
step three markers over a 1.9e-10 m3 sliver of beta drew 2.0e-9 m3). The total stays exact; where the metal is drawn
from is then the largest piece, not the sliver beneath the markers. Draws that ask more than the largest piece then
holds raise RuntimeError: nothing is left to give them.

At the macro step's end (`_rebuild`): phi = `outer_surface(beta, H)` with every node taking the thickness of the marker
whose base point is nearest; psi = `molten_rigid_surface(phi, T)` with T the volume's field inside beta, each skin's
cell at the node's depth between beta and phi, and the top cell outside phi; the mesh rebuilt from
`criteria.make_target` only where it no longer meets the targets (`_needs_rebuild`, with slack: a rebuild at every
step does not settle), with `CutThermal.remesh` (beta carried, reset and corrected to the exact volume; the remesh's
enthalpy residual booked as pending heat); phi and psi again; the markers settled (`_settle`): snapped onto phi to
SNAP_TOL, the trapped and the patchless released (their metal to the volume by a beta move, their momentum to the
volume), tessellated, their skins re-partitioned (below), split and merged, snapped, tessellated and re-partitioned
again, the skins regridded and given their markers' areas and capped curvatures.

The re-partition (Task 9b, the controller's ruling of 2026-11-12). A tessellation hands the surface out afresh: most of
a marker's area change at a step end is area passing between markers, not the surface stretching (measured on the test
cube and the full-size cube, Task 9b's report), and a skin that kept its metal on its new area thickened or thinned and
moved metal through its base for no reason (a backlog handed down 10 mm at once on the full-size cube). Each
tessellation is followed by `_repartition`: from each marker's old owners (its position and normal at the tessellation
before, `Markers.x_tess`), its skin's area and the surface then (the step's old phi for the first, phi as it is for the
one after split and merge), and the new tessellation, `markers.transfer_matrix` gives every old column's stretch and
the transfer T[j, i], the area of old column j that new column i receives (`repartition`). Each old column is scaled by its stretch
(its top area becomes its old area times its stretch, its metal unchanged: section 9.2's draw or return at its base
follows in the next part, by the draw to the minimum mass and the hand-down, as for any stretched or gathered skin),
then cut into the fractions T[j, i] / (its row's sum) of its profile (`Skins.cut`) and each new column is the hand-over
mix (`Skins.combine`) of what it receives. Mass, enthalpy and momentum add to round-off; the volume is untouched. Each
re-partition is kept in `last_settle` (the transfer, the stretch, the old and new areas) for the step's books.

A skin deeper than its radius of curvature (`Skins.dz` raises for all columns if one is): each column's principal
curvatures are held, either way, at KAPPA_DEPTH_SHARE over its flat depth (`cap_curvature`), so the cross-section
A0 (1 - k1 z)(1 - k2 z) stays positive to about three times the column's depth; only the columns over the cap change.
A surface bending more tightly than its skin's depth is sub-grid (design 20.6, known limitation 1); the markers keep
phi's curvature for the capillary push."""
from collections import namedtuple

import numpy as np
from scipy.spatial import cKDTree

from .boxmesh import build
from .criteria import finite_values, kappa_max_nodal, make_target
from .cutcell import inside_volumes
from .levelset import correct_volumes, pieces, redistance
from .markers import NEWTON_REACH
from .repartition import SAMPLE_SHARE, apply_transfer, surface_samples, transfer_matrix
from .surfaces import molten_rigid_surface, outer_surface
from .transfer import interpolate, locate

KAPPA_DEPTH_SHARE = 0.1           # a column's principal curvature at most this over its flat depth (module docstring)
SNAP_TOL = 1.0e-10                # m: the markers are snapped onto phi until none is farther from its zero contour
SNAP_ROUNDS = 20                  # rounds of `Markers.snap` at most
PIECE_ROUNDS = 8                  # rounds at most in which a move's draws use up pieces of beta (`_move_beta`)
PIECE_GONE = 1.0e-9               # m: the level set a used-up piece's nodes are lifted to, just outside the metal
LEAD = 0.6e-3                     # m: the fine region's lead beyond the margin (criteria.make_target's default)
REGION_SLACK = 0.5                # fine cells: a rebuild is needed once a leaf is too coarse for a lead this much shorter
CURVATURE_SLACK = 2.0             # or once a leaf is more than this times its curvature's target size
COARSEN_SHARE = 0.1               # of the metal's volume: leaves that could be COARSEN_SIZES times coarser before the mesh
MODES = ("full", "top")           # a marker's skin: full-depth (Task 9, everywhere) or top-zone (Task 10)
COARSEN_SIZES = 4.0               # is rebuilt (two sizes: the closure's grading keeps one-size steps finer than asked)

BARE_ROUNDS = 4                   # tessellations at most, each releasing the markers the last left no area (`_tessellate`)

Repartition = namedtuple("Repartition", "transfer stretch old_area area")   # one re-partition of a step end


def cap_curvature(skins, kappa, depth=None):
    """Each column's principal curvatures (n, 2) [1/m] held within KAPPA_DEPTH_SHARE over the column's flat depth either
    way (module docstring), `depth` [m] (n,) or the skins' flat depth sum m / (rho(T) A0), at least the stub (fix round
    3: a column runoff had all but emptied kept phi's sharp curvature and the draw then refilled it to the stub past
    1/kappa); columns within it are unchanged. Convex curvatures are capped so that the cross-section stays positive; concave ones too, since a column
    widening with depth raises its own minimum mass, rho_s A(z_bottom) stub, faster than its mass where |k| stub > 1, so
    that its draws would run away (fix round 1 of Task 9)."""
    kappa = np.asarray(kappa, dtype=float)
    if depth is None:
        T, valid = skins.T(), skins.valid
        vol = np.zeros(skins.m.shape)
        vol[valid] = skins.m[valid] / skins.density(T[valid])
        depth = vol.sum(axis=1) / skins.area
    depth = np.maximum(np.asarray(depth, dtype=float), skins.settings.stub)   # the draw refills an emptied column
    limit = KAPPA_DEPTH_SHARE / depth
    return np.clip(kappa, -limit[:, None], limit[:, None])


class StepEnd:
    """The step's end of `LevelSetSkinBody` (module docstring). Every method that defers heat to the volume adds its
    size to the step's tally (`tally.move` for the moves, `tally.remesh` for the remesh), and a release adds the released
    skins' momentum to `tally.to_volume`."""

    def _target(self, mesh, phi, psi, curvature=True, lead=LEAD):
        band = 2.0 * self.base.h0
        if curvature:
            kp, ks = kappa_max_nodal(mesh, phi, band=band), kappa_max_nodal(mesh, finite_values(psi), band=band)
        else:
            kp = ks = np.zeros(len(phi))
        return make_target(mesh, phi, psi, kp, ks, fine=self.fine_size, margin=self.fine_margin, lead=lead,
                           floor=self.cell_floor)

    def _row_marker(self):
        """Each skin row's marker (n,)."""
        out = np.empty(self.skins.n, dtype=np.int64)
        out[self.markers.skin_index] = np.arange(self.markers.n)
        return out

    def _bases(self):
        """Each skin row's base point [m] (n, 3), x - H n (Ruling 1)."""
        mk, sk = self.markers, self.skins
        r = self._row_marker()
        return mk.x[r] - sk.thickness()[:, None] * mk.normal[r]

    def _thickness_nodal(self):
        """The skins' thickness [m] at every node, from the marker whose base point is nearest (design 20.3)."""
        _, near = cKDTree(self._bases()).query(self.mesh.points)
        return self.skins.thickness()[near]

    def _sync(self):
        """The skins take their markers' areas and (capped) curvatures."""
        mk, sk = self.markers, self.skins
        sk.area[mk.skin_index] = mk.area
        sk.kappa[mk.skin_index] = mk.kappa
        sk.kappa = cap_curvature(sk, sk.kappa)

    def _cap(self):
        self.skins.kappa = cap_curvature(self.skins, self.skins.kappa)

    def _temperature_field(self, phi):
        """The temperature at every node: the volume's inside beta; between beta and phi each skin's cell at the node's
        depth below phi (the nearest marker's), and outside phi its top cell (so that psi sees no crossing there)."""
        T = self.solver.T.copy()
        out = self.beta >= 0.0
        if not out.any():
            return T
        sk = self.skins
        Tc = sk.T()
        dz = sk.dz(Tc)
        z_bot = np.where(sk.valid, np.cumsum(dz, axis=1), np.inf)
        nodes = np.flatnonzero(out)
        _, near = cKDTree(self.markers.x).query(self.mesh.points[nodes])
        r = self.markers.skin_index[near]
        depth = np.maximum(-phi[nodes], 0.0)
        cell = np.minimum((z_bot[r] < depth[:, None]).sum(axis=1), np.maximum(sk.count[r] - 1, 0))
        T[nodes] = np.where(sk.count[r] > 0, Tc[r, cell], T[nodes])
        return T

    def _psi(self, phi):
        psi, _ = molten_rigid_surface(self.mesh, phi, self._temperature_field(phi), self.material.T_rigid,
                                      band=self.band)
        return psi

    def _book_move(self, points, weight, J):
        """Book J [J] owed to the volume after a move, at the points in proportion to `weight`."""
        if J == 0.0 or not len(points):
            return
        w = weight / weight.sum() if weight.sum() > 0.0 else np.full(len(weight), 1.0 / len(weight))
        keep = w > 0.0
        self.add_pending(points[keep], J * w[keep])

    def _move_beta(self, points, dV, disp, tally):
        """Move beta by the volume changes dV [m3] at the points [m] (k, 3), each displaced along its normal by disp [m]
        (outward positive), extended to the nodes within fine_margin of beta from the nearest point; reset, corrected per
        piece to the exact volume target, and `CutThermal.move`. Returns the enthalpy the move removed from the volume
        [J]."""
        mesh, beta = self.mesh, self.beta
        new = beta.copy()
        if np.any(disp != 0.0):
            sel = np.flatnonzero(np.abs(beta) < self.fine_margin)
            _, near = cKDTree(points).query(mesh.points[sel])
            new[sel] = beta[sel] - disp[near]
            new = redistance(mesh, new, band=self.band)
        for _ in range(PIECE_ROUNDS):
            label, n = pieces(mesh, new)
            if n == 0:
                raise RuntimeError("the skin base has no metal left after a move")
            targets = self._piece_targets(beta, new, label, n, points, dV)
            used = np.flatnonzero(targets[1:] <= 0.0) + 1
            if not len(used):
                break
            # a piece the draws use up leaves, and what they ask beyond it the largest piece gives (fix round 5)
            gone = np.unique(mesh.tets[np.isin(label, used)])
            new[gone] = np.maximum(new[gone], PIECE_GONE)
            tally.pieces_used += len(used)
        else:
            raise RuntimeError("pieces kept being used up by a move after {} rounds".format(PIECE_ROUNDS))
        if targets[0] <= 0.0:                              # no piece is left to give what the draws ask
            raise RuntimeError("the draws ask {:.3e} m3 more than the skin base's largest piece holds".format(-targets[0]))
        new, worst = correct_volumes(mesh, new, label, targets)
        tally.correction = max(tally.correction, worst)
        removed = self.solver.move(new)
        self.beta = new
        return removed

    def _piece_targets(self, old, new, label, n, points, dV):
        """The exact volume of each piece of `new` [m3]: the old metal in its cells (the old metal left in no piece goes
        to the largest) plus the changes dV at the points that lie in it (the largest's where outside every piece)."""
        mesh = self.mesh
        v_old = inside_volumes(mesh, old)
        cell, _ = locate(mesh, points)
        where = np.maximum(label[cell], 0)
        # summed piece by piece with numpy's pairwise sums: bincount's running sum lost 1e-13 of the volume a move
        t = np.array([v_old[label == p].sum() + dV[where == p].sum() for p in range(n)])
        t[0] += v_old[label < 0].sum()
        return t

    def _rebuild(self, tally):
        """phi and psi, the mesh rebuilt from the refinement targets with the solve carried to it, then phi and psi on
        it and the markers settled (module docstring). Returns whether the mesh stayed as it was."""
        mesh = self.mesh
        surface = self._samples()                          # the surface the markers' areas were given on
        phi = outer_surface(mesh, self.beta, self._thickness_nodal(), band=self.band)
        psi = self._psi(phi)
        target = self._target(mesh, phi, psi)
        same = not self._needs_rebuild(target, self._target(mesh, phi, psi, False, LEAD - REGION_SLACK * self.fine_size))
        if not same:
            new = build(self.base, target)
            beta = redistance(new, interpolate(mesh, self.beta, new.points), band=self.band)
            v_old = inside_volumes(mesh, self.beta)
            label, n = pieces(new, beta)
            if n == 0:
                raise RuntimeError("the skin base has no metal on the rebuilt mesh")
            targets = np.zeros(n)
            targets[0] = v_old.sum()
            if n > 1:                                      # the pieces share the old metal by their volumes here
                v_new = inside_volumes(new, beta)
                share = np.array([v_new[label == p].sum() for p in range(n)])
                targets = share * (v_old.sum() / share.sum())
            beta, worst = correct_volumes(new, beta, label, targets)
            tally.correction = max(tally.correction, worst)
            E0 = self.solver.energy()
            self.solver.remesh(new, beta)
            self.mesh, self.beta = new, beta
            owed = E0 - self.solver.energy()
            tally.remesh += abs(owed)
            self._book_move(self._bases(), self.skins.area, owed)
            phi = outer_surface(new, self.beta, self._thickness_nodal(), band=self.band)
        self.phi = phi
        self.psi = self._psi(phi)
        self._settle(tally, surface)
        return same

    def _needs_rebuild(self, target, region):
        """Whether the mesh must be rebuilt (a rebuild builds every leaf to `target`): some leaf coarser than the fine
        region asks with its lead shortened by REGION_SLACK fine cells (`region`, no curvature), some leaf more than
        CURVATURE_SLACK times its full target, or more than COARSEN_SHARE of the metal (by volume) in leaves that could be
        COARSEN_SIZES times coarser. A rebuild at every step would not settle: the skin base carried to a new mesh is a
        slightly different contour, and phi's recovered curvature a slightly different field, which move leaves across
        the targets' thresholds and so change the next mesh (measured on the test cube: four rebuilds in a row, each
        with some 36 leaves over their curvature target)."""
        mesh = self.mesh
        size, x = mesh.size, mesh.centroids
        want = np.asarray(target(x, size), float)
        if np.any(size > np.asarray(region(x, size), float) * (1.0 + 1e-9)):
            return True
        if np.any(size > CURVATURE_SLACK * want):
            return True
        v = inside_volumes(mesh, self.beta)                 # the metal's leaves; closure grades the rest anyway
        coarser = COARSEN_SIZES * size <= want
        return float(v[coarser].sum() / v.sum()) > COARSEN_SHARE

    def _snap(self, tally):
        """The markers onto phi, `Markers.snap` round after round until none is farther than SNAP_TOL from its zero
        contour, at most SNAP_ROUNDS rounds, the last finishing every marker by its cell's Newton steps whatever its
        distance (`Markers.snap`'s reach; fix round 4); the step's `tally` counts that round's fallback markers and largest
        sideways slide (`snap_fallback`, `snap_slide`; fix round 5)."""
        mk = self.markers
        for k in range(SNAP_ROUNDS):
            if np.abs(interpolate(self.mesh, self.phi, mk.x)).max() <= SNAP_TOL:
                return
            last = k == SNAP_ROUNDS - 1
            fallback, slide = mk.snap(self.mesh, self.phi, reach=np.inf if last else NEWTON_REACH)
            if last and tally is not None:
                tally.snap_fallback += fallback
                tally.snap_slide = max(tally.snap_slide, slide)
        if np.abs(interpolate(self.mesh, self.phi, mk.x)).max() <= SNAP_TOL:
            return
        raise RuntimeError("the markers did not snap onto phi within {} m in {} rounds".format(SNAP_TOL, SNAP_ROUNDS))

    def _release(self, idx, tally):
        """Release the markers idx: their skins' metal and enthalpy return to the volume (beta moved by their volume,
        the move's enthalpy against theirs booked as pending at their points), their momentum to the volume
        (`tally.to_volume`)."""
        mk = self.markers
        x, mass, H, P = mk.release(idx)
        tally.to_volume += P.sum(axis=0)
        removed = self._move_beta(x, mass / self.rho, np.zeros(len(x)), tally)
        owed = float(H.sum()) + removed
        tally.move += abs(owed)
        self._book_move(x, mass, owed)
        self.phi = outer_surface(self.mesh, self.beta, self._thickness_nodal(), band=self.band)
        self.psi = self._psi(self.phi)

    def _samples(self):
        """phi's surface samples for the re-partition, SAMPLE_SHARE of the target spacing across (`surface_samples`)."""
        return surface_samples(self.mesh, self.phi, SAMPLE_SHARE * self.spacing)

    def _owners(self, surface):
        """The markers' old owners for a re-partition: [their positions and normals at the last tessellation, their
        skins' areas, the samples `surface` of the surface those areas were given on]."""
        mk = self.markers
        return [mk.x_tess.copy(), mk.normal_tess.copy(), self.skins.area[mk.skin_index].copy(), surface]

    def _repartition(self, old, surface):
        """The skins follow the area the last tessellation re-partitioned (module docstring): every old column scaled by
        its stretch, cut by the transfer and mixed into its receivers by the hand-over rule; the re-partition kept in
        `last_settle`. `old` from `_owners`, `surface` the samples of phi now."""
        mk, sk = self.markers, self.skins
        old_x, old_normal, old_area, old_surface = old
        T, stretch = transfer_matrix(old_x, old_normal, old_area, old_surface, mk.x, mk.area, mk.faces, surface)
        apply_transfer(sk, mk.skin_index, old_area, T, stretch)
        self.last_settle.append(Repartition(T, stretch, old_area, mk.area.copy()))

    def _drop_owners(self, old, idx):
        """`old` (from `_owners`) without the rows of the markers idx, which are about to be released."""
        stay = np.setdiff1d(np.arange(self.markers.n), idx)
        return [a[stay] for a in old[:3]] + [old[3]]

    def _tessellate(self, old, tally):
        """The markers tessellated on phi; those it leaves no area (the patchless) released and the rest snapped and
        tessellated again until none is left, at most BARE_ROUNDS times, their rows leaving `old` (fix round 1 of Task 9b:
        after every tessellation of the step end, not only the first). Returns `old`."""
        mk = self.markers
        for _ in range(BARE_ROUNDS):
            mk.tessellate(self.mesh, self.phi)
            bare = np.flatnonzero(mk.area <= 0.0)
            if not len(bare):
                return old
            old = self._drop_owners(old, bare)
            self._release(bare, tally)
            self._snap(tally)
        raise RuntimeError("markers were still left no area after {} tessellations".format(BARE_ROUNDS))

    def _settle(self, tally, surface=None):
        """The markers onto phi, the trapped (and the patchless) released, tessellated and their skins re-partitioned,
        split, merged, tessellated and re-partitioned again, and the skins regridded with their markers' areas and
        curvatures. `surface`: the samples of the surface the markers' areas were last given on (`_rebuild` passes the
        step's old phi), phi as it is if None."""
        mk = self.markers
        self.last_settle = []
        old = self._owners(self._samples() if surface is None else surface)
        self._snap(tally)
        cell, _ = locate(self.mesh, mk.x)
        trapped = mk.trapped_markers(self.mesh, self.phi, self.mesh.size[cell])
        if len(trapped):
            old = self._drop_owners(old, trapped)
            self._release(trapped, tally)
            self._snap(tally)
        old = self._tessellate(old, tally)
        now = self._samples()
        self._repartition(old, now)
        self._sync()
        big = mk.split(self.spacing ** 2)
        merged = mk.merge(self.spacing)
        if len(big) or merged:
            old = self._owners(now)
            self._snap(tally)
            old = self._tessellate(old, tally)
            self._repartition(old, self._samples())
        self._sync()                                       # a merged column holds both columns' metal: cap it first
        self.skins.regrid()
        self._sync()
        self.mode = np.full(mk.n, MODES[0])          # every skin full-depth in Task 9; Task 10 sets and carries the modes
