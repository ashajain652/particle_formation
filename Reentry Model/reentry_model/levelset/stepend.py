"""The layer and the step's end of the level set's skinned body (design sections 20.3-20.6 and 21; level-set stage 3 Task
9, split from `skinbody` in its fix round; plan 03b Task 5): the geometry that follows the skins -- the gas surface phi
moved by each column's volume, the skin base beta derived one skin depth inside it, the division of the layer between
them, every column filled exactly, the molten-rigid surface psi rebuilt, the mesh rebuilt where it no longer meets its
refinement targets, and the markers settled on phi -- as `StepEnd`, a mixin of `skinbody.LevelSetSkinBody`, whose
attributes it reads (mesh, phi, beta, psi, division, carry, solver, markers, skins, material, base, settings, spacing,
band, fine_margin, fine_size, cell_floor, rho). The physics is the geometry modules' (`levelset`, `surfaces`,
`division`, `criteria`, `boxmesh`, `markers`) and the shared package's; every number here is a numerical setting.

The two surfaces (design 21.4). phi is the evolved level set; beta is always `surfaces.skin_base(phi, H)`, the zero
contour of phi + H (H = `SkinSettings.depth`), the boundary of the volume's metal (the CutThermal solve runs on
{beta < 0}). The layer between them is divided among the markers (`division.divide`): column i, marker i's, has the
cross-section A_i(z) at depth z, its volume per unit depth, and the skins take it as their profile
(`Skins.set_profile`). A column whose base A_i(H) is at most BASELESS_SHARE of its top A_i(0) has no base (beside a
convex edge it pinches before H): its bottom is insulated, and it never exchanges metal with the volume.

Moving phi (`_move_phi`, planning ruling 3 steps 2-5): each column's volume change dV_i -- its cells' change since the
last fill, plus the density carry below -- is swept by its marker's group (`levelset.move_by_volumes`, each group's
swept volume exact to MOVE_RTOL; a group asked more than its region holds misses, and the per-piece correction makes up
what it missed: the controller's ruling T2-e); the groups' swept volumes are checked against phi's change of volume
(MOVE_SUM_RTOL of the body's: cells no marker may own are in no group); phi is corrected per piece to its metal's
volume, the skins' cells at their own densities, the volume's metal (beta's) and the carry (`levelset.correct_volumes`);
beta is derived afresh, and `CutThermal.move` fills the volume's new nodes (an island the old metal's extension cannot
reach takes the skins' temperature at its nodes, `_temperature_field`), the enthalpy it removed being the caller's to
book as heat owed to the volume. A piece of phi the move parts from the body with no column and none of the volume's
metal in it (a speck) is lifted out, its metal the largest piece's, and its nodes are left out of the correction's
claims, so that it stays out (plan 03b Task 7b, the controller's ruling T7-f: the body's shift over its two rings of
cells took a lifted node back below zero whenever the body grew by more than PIECE_GONE, and the speck came back with no
target: on the spraying cube it captured two drifted markers' columns and the split never ended, and on the shoulder no
group counted the cells that swallowed it; diagnosis D2). The markers stay where the step began until its end, but
their areas, the top fluxes', are given afresh on the moved phi at every acceptance (`_own`), and the top fluxes are
taken again under the same ownership before the next part (`skinbody.LevelSetSkinBody._take_loads`).

Clear stops (plan 03b Task 7b, the controller's ruling T7-g): a division that finds a piece of phi with layer but no
column, or a piece whose closure factor departs from one by more than DIVISION_GROSS, raises (`_divide`); so does a fill
that leaves a column beyond FILL_RTOL of its volume after its rounds (`_fill`), naming it. Every completed step of every
recorded case filled to 4e-16 of a column's volume at worst and closed within 4e-3, so neither fires on a good step.

Filling the columns (`_fill`, planning ruling 3 step 7). After the division (`_divide`), every column whose cells'
volume m / rho(T) differs from its own by more than FILL_RTOL is filled, in at most FILL_ROUNDS rounds: first the columns
with no base, from the neighbour across the boundary they share deepest (preferring a neighbour with a base): a short
one takes the metal its neighbour holds just above that depth (`Skins.take_between`), as a new bottom cell; an overfull
one gives its bottom excess into the neighbour's cell at that depth; then every column with a base at its base
(`Skins.exchange_base`): a short one draws the volume's metal at its base temperature, the area mean over its group of
beta's patch temperatures (the base law's T_b, `skinbody.LevelSetSkinBody._interface`; the controller's ruling T5-k:
read at the base point x - H n instead, it fell outside beta at every edge and corner marker of the test cube, by up to
0.55 mm at rest, where the solve holds its cut cells' extensions or stale values), whatever its phase (design 21.3, the
volume's own temperature), and an overfull one returns its bottom cells whatever their phase. (Task 5 first capped the
draw at T_rigid, as stage 3's full-depth draw had; with the markers' areas given afresh within a step, `_own`, and the
draw at the group's temperature, the test cube at 6 MW/m2 with the liquid fed off took 10 parts in a 0.5 s step
without the cap against 9 with it, and 3 against 3 then 8 against 8 in two 0.25 s steps, measured once each, so the
controller's ruling T5-d removed it.) What a fill's draw takes of the volume's enthalpy is owed to the volume at the base, as any
fill's difference is. Metal changes owner, never place: the volume's one density is the
solid's (design 21.6), so what crosses the base at another density changes the metal's true volume against the
geometry: that difference, the cells' volume change less the exchanged mass over rho_s, is each column's carry (planning
ruling 4), a volume change at the top of that column taken into the next move of phi (the next part's, or the step
end's); so too the cells' volume changes no move of phi took, the mixing of columns at a step end (a Scheil alloy's
volumes do not add) and of metal passed between columns by a fill. The volume's mass is rho_s (V(beta) + the carry)
(`skinbody.LevelSetSkinBody.books`), so the books stay exact while phi settles one move late.

At the macro step's end (`_end_step`): the last part's carry moves phi (ruling 4); psi from phi (`_psi`); the mesh
rebuilt from `criteria.make_target` only where it no longer meets the targets (`_needs_rebuild`, with slack), phi
carried to it (interpolated, reset keeping its contour, corrected to its pieces' volumes; the controller's ruling T2-b;
a piece the carried phi parts off with no column and none of the volume's metal is lifted out first by the speck rule of
the moves, `_lift_specks`, its nodes left out of the correction's claims: the controller's ruling T7b-b, after the stage
3 box's second step met one of 2.3e-16 m3 in its rim, which the division then stopped on), beta derived, and
`CutThermal.remesh` (its enthalpy residual booked as heat owed); the markers settled (`_settle`): snapped onto phi to
SNAP_TOL, the trapped and the patchless released (each giving its skin to the nearest marker that
may own its part of phi, by depth, `Markers.release`; planning ruling 6), given their areas of phi (`_tessellate`), the
layer divided among them and their skins re-partitioned across that division (below), split and merged, snapped, given
their areas, divided and re-partitioned again, the skins regridded and given their markers' areas; then given the last
division's profile and filled, and the thin-metal stop (THIN_METAL, design 21.11) checked. The step end's own settle and
fill leave a carry of their own (the re-partition's mixing and the drawn metal's density difference), which waits for
the next step's first move (the controller's ruling T5-i, reading planning ruling 4's "before the closure check" as a
bound on the carry left): until then phi misses its metal's true volume by it, the step's `density_carry_m3` (largest
share of a column, `density_carry_share_max`), which the tests bound.

The re-partition (Task 9b, the controller's ruling of 2026-11-12; by depth across the division since plan 03b Task 6,
the controller's ruling P4). A settle hands the layer out afresh: most of a column's change at a step end is layer
passing between markers. Each giving of the markers' areas (`_tessellate`, which releases the markers it leaves no area
and gives them again) is followed by `_repartition`: the layer divided among the markers as they now stand
(`division.divide`), and from each marker's old owners (its position and normal at the last giving, `Markers.x_tess`),
its old column's volume (its skin's profile, the last division's) and the surface then, and the new division,
`repartition.transfer_matrix` gives every old column's stretch and the transfer T[j, i], the volume of old column j that
new column i receives, across the new division's boundaries; `repartition.apply_transfer` stretches each old column
(its profile with it), cuts it into the fractions T[j, i] / (its row's sum) of itself (`Skins.cut`) and gives each new
column what it receives by depth (`Skins.combine(..., by="depth")`). Mass, enthalpy and momentum add to round-off; the
volume is untouched; the step end then gives the skins the last division's profile and the fill settles each column.
Each re-partition is kept in `last_settle`.

No rigid metal left (the controller's ruling (d) of stage 3 Task 10): a step end whose psi finds no rigid metal raises
RuntimeError saying so (`_check_core`): the molten metal has nothing to stand on, and stage 3 holds the volume's melt
(no flow solve). What is checked is that no metal below T_rigid is left anywhere; a rigid crust over a molten core
passes (psi then takes the crust for the core)."""
from collections import namedtuple

import numpy as np
from scipy.spatial import cKDTree

from .boxmesh import build
from .criteria import finite_values, kappa_max_nodal, make_target
from .cutcell import inside_volumes, patches
from .division import _marker_pieces, divide
from .levelset import correct_volumes, move_by_volumes, pieces, redistance
from .markers import NEWTON_REACH, OPPOSED_NORMALS
from .repartition import SAMPLE_SHARE, apply_transfer, surface_samples, transfer_matrix
from .surfaces import molten_rigid_surface, skin_base
from .transfer import interpolate, locate

SNAP_TOL = 1.0e-10                # m: the markers are snapped onto phi until none is farther from its zero contour
SNAP_ROUNDS = 20                  # rounds of `Markers.snap` at most
LEAD = 0.6e-3                     # m: the fine region's lead beyond the margin (criteria.make_target's default)
REGION_SLACK = 0.5                # fine cells: a rebuild is needed once a leaf is too coarse for a lead this much shorter
CURVATURE_SLACK = 2.0             # or once a leaf is more than this times its curvature's target size
COARSEN_SHARE = 0.1               # of the metal's volume: leaves that could be COARSEN_SIZES times coarser before the mesh
COARSEN_SIZES = 4.0               # is rebuilt (two sizes: the closure's grading keeps one-size steps finer than asked)
BARE_ROUNDS = 4                   # tessellations at most, each releasing the markers the last left no area
                                  # (`_tessellate`)
FILL_RTOL = 1e-12                 # of a column's volume: a column whose cells' volume is within this of it is full
                                  # (plan 03b Task 5; the books' mass tolerance, so that what a full column leaves is
                                  # within the books)
FILL_ROUNDS = 8                   # rounds of the fill at most (each a pass over the columns with no base and one over
                                  # those with one; the second settles the property table's round-off of the first)
BASELESS_SHARE = 1e-6             # a column whose base A(H) is at most this share of its top A(0) has no base (plan 03b
                                  # Task 5: its column pinches before H beside a convex edge)
PIECE_ROUNDS = 8                  # rounds at most in which a move leaves pieces of phi with no metal (`_move_phi`)
PIECE_GONE = 1.0e-9               # m: the level set such a piece's nodes are lifted to, just outside the metal
AREA_FLOOR = 1e-6                 # of the target spacing squared: the least top area a marker is given within a step
                                  # (`_own`; one whose part of phi the move took keeps this sliver until the step end
                                  # releases it, so that the runoff's operators, over the top areas, stay finite)
MOVE_SUM_RTOL = 1e-12             # of the body's volume: the groups' swept volumes must sum to phi's change of volume to
                                  # this (the controller's ruling 3 of Task 5's dispatch)
DIVISION_GROSS = 0.1              # a piece's closure factor more than this from one stops the run (`_divide`; the
                                  # controller's ruling T7-g: the speck that took two columns on the spraying cube closed
                                  # at 0.0033, every recorded step within 4e-3; DIVISION_TOL, 1e-3, is the tests' bound)
NO_CORE = ("the body has no rigid metal left: psi finds no metal below T_rigid ({:.1f} K) anywhere in it, so its molten "
           "metal has nothing to stand on (it has melted through); that is beyond level-set stage 3, which holds the "
           "volume's melt (no flow solve)")
THIN_METAL = "metal thinner than twice the skin depth"   # the message's start (design section 21.11)

Repartition = namedtuple("Repartition", "transfer stretch old_volume volume")   # one re-partition of a step end [m3]


class StepEnd:
    """The layer and the step's end of `LevelSetSkinBody` (module docstring). Every method that defers heat to the
    volume adds its size to the step's tally (`tally.move` for the moves and fills, `tally.remesh` for the remesh)."""

    # -- checks -------------------------------------------------------------------------------------------------------
    def _check_core(self, psi):
        """RuntimeError (NO_CORE) when psi finds no rigid metal: no node of the metal below T_rigid anywhere (a rigid crust
        over a molten core passes: `surfaces.molten_rigid_surface` then takes it for the core)."""
        if not np.any(np.asarray(psi) < 0.0):
            raise RuntimeError(NO_CORE.format(self.material.T_rigid))

    def _check_thin(self):
        """RuntimeError starting with THIN_METAL where two columns whose normals oppose share a boundary in the layer
        (`Division.opposed`, design section 21.11: the metal there is thinner than twice the skin depth, or an edge is
        sharper than about 60 degrees), naming the place; the thin case is designed later."""
        d = self.division
        if d is None or not d.opposed.any():
            return
        e = np.flatnonzero(d.opposed)
        a, b = self.markers.x[d.i[e[0]]] * 1e3, self.markers.x[d.j[e[0]]] * 1e3
        raise RuntimeError("{}: {} pair(s) of columns whose normals oppose (n_i . n_j < {}) share a boundary in the "
                           "layer, the first between the markers at ({:.3f}, {:.3f}, {:.3f}) and ({:.3f}, {:.3f}, {:.3f}) "
                           "mm (design section 21.11; the skin depth is {:.3g} mm, and the thin case is designed "
                           "later)".format(THIN_METAL, len(e), OPPOSED_NORMALS, *a, *b, self.settings.depth * 1e3))

    # -- fields ---------------------------------------------------------------------------------------------------------
    def _target(self, mesh, phi, psi, curvature=True, lead=LEAD):
        band = 2.0 * self.base.h0
        if curvature:
            kp, ks = kappa_max_nodal(mesh, phi, band=band), kappa_max_nodal(mesh, finite_values(psi), band=band)
        else:
            kp = ks = np.zeros(len(phi))
        return make_target(mesh, phi, psi, kp, ks, fine=self.fine_size, margin=self.fine_margin, lead=lead,
                           floor=self.cell_floor)

    def _reset(self, mesh, phi):
        """phi reset to a signed distance within the band, keeping its zero contour (`redistance(..., keep=True)`; the
        controller's ruling T2-b: the plain reset rounds convex edges at every call)."""
        return redistance(mesh, phi, band=self.band, keep=True)

    def _row_marker(self):
        """Each skin row's marker (n,)."""
        out = np.empty(self.skins.n, dtype=np.int64)
        out[self.markers.skin_index] = np.arange(self.markers.n)
        return out

    def _row_normals(self):
        """Each skin row's marker's unit normal (n, 3)."""
        out = np.zeros((self.skins.n, 3))
        out[self.markers.skin_index] = self.markers.normal
        return out

    def _bases(self):
        """Each skin row's base point [m] (n, 3), x - H n: one skin depth below its marker along its normal."""
        mk = self.markers
        r = self._row_marker()
        return mk.x[r] - self.settings.depth * mk.normal[r]

    def _sync(self):
        """The skins take their markers' areas (the top fluxes' areas, the controller's ruling T1-g) and curvatures."""
        mk, sk = self.markers, self.skins
        sk.area[mk.skin_index] = mk.area
        sk.kappa[mk.skin_index] = mk.kappa

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

    # -- the layer: phi moved, beta derived, the division and the fill ----------------------------------------------
    def _needs_move(self, dV):
        """Whether the volume changes dV [m3] (n,) per row ask phi to move: some column's change is beyond FILL_RTOL of
        its volume (a column with no volume, beyond none)."""
        V = self.skins.column_volume()
        return bool(np.any(np.abs(dV) > FILL_RTOL * V))

    def _cells_change(self):
        """Each row's volume change [m3] (n,) since its column was last filled: its cells' volume m / rho(T) less its
        column's volume (planning ruling 3, step 1), plus its carry (ruling 4)."""
        sk = self.skins
        return sk.cells_volume() - sk.column_volume() + self.carry

    def _piece_targets(self, phi, beta_old, phi_old, mesh=None):
        """Each piece of phi's metal volume [m3], label and count: the cells' volume (m / rho(T)) and the carry of the
        rows whose columns it holds, and beta_old's metal in its cells (where outside every piece, the largest's). A
        row's piece is the piece of phi that holds most of the metal its column's piece of `phi_old` held (the division's
        pieces, `Division.piece`, numbered as `levelset.pieces(phi_old)`): within a step the markers stand where it
        began, off the moved phi, and placed by position one could fall to a speck the move parted off, which would then
        be filled with its column. phi_old and beta_old are on the body's mesh; phi is too, or on `mesh` (the remesh's
        new mesh, `_rebuild`; the controller's ruling T7b-b), where each old cell's metal is placed by its centroid
        (`transfer.locate`) and an old piece's new piece is the one holding most of its old metal."""
        old_mesh, mk, sk = self.mesh, self.markers, self.skins
        mesh = old_mesh if mesh is None else mesh
        label, n = pieces(mesh, phi)
        if n == 0:
            raise RuntimeError("the gas surface holds no metal after a move")
        old, n_old = pieces(old_mesh, phi_old)
        if mesh is old_mesh:                                # each old cell's new piece, and what weighs it
            at_new, weight = label, inside_volumes(mesh, phi)
        else:
            at_new, weight = label[locate(mesh, old_mesh.centroids)[0]], inside_volumes(old_mesh, phi_old)
        image = np.zeros(max(n_old, 1), dtype=np.int64)          # each old piece's new piece: where most of it lies
        for q in range(n_old):
            at = (old == q) & (at_new >= 0)
            if at.any():
                image[q] = int(np.argmax(np.bincount(at_new[at], weight[at], minlength=n)))
        if self.division is not None and len(self.division.piece) == mk.n:
            was = self.division.piece
        else:
            was = _marker_pieces(old_mesh, mk.x, old, patches(old_mesh, phi_old))
        piece = image[np.clip(was, 0, max(n_old - 1, 0))]
        held = (sk.cells_volume() + self.carry)[mk.skin_index]
        v_beta = inside_volumes(old_mesh, beta_old)
        where = np.where(at_new >= 0, at_new, 0)
        # summed piece by piece with numpy's pairwise sums (bincount's running sum lost 1e-13 of a move's volume)
        t = np.array([held[piece == p].sum() + v_beta[where == p].sum() for p in range(n)])
        return t, label, n

    def _lift_specks(self, mesh, phi, beta_old, phi_old, tally, what):
        """phi (on `mesh`) with every piece but the largest that is owed no metal (`_piece_targets`: no column maps to
        it and none of beta_old's metal lies in it) lifted to PIECE_GONE, round after round (at most PIECE_ROUNDS), its
        metal the largest piece's; the lifted pieces counted in `tally.pieces_used`. Returns (phi, the pieces' targets,
        label, count, the lifted nodes (bool per node), to be left out of the correction's claims: ruling T7-f).
        `what` names the caller in the error raised if pieces keep being left with no metal."""
        lifted = np.zeros(len(phi), dtype=bool)
        for _ in range(PIECE_ROUNDS):
            targets, label, n = self._piece_targets(phi, beta_old, phi_old, mesh)
            gone = np.flatnonzero(targets[1:] <= 0.0) + 1
            if not len(gone):
                return phi, targets, label, n, lifted
            nodes = np.unique(mesh.tets[np.isin(label, gone)])
            phi = phi.copy()
            phi[nodes] = np.maximum(phi[nodes], PIECE_GONE)
            lifted[nodes] = True
            tally.pieces_used += len(gone)
        raise RuntimeError("pieces of phi kept being left with no metal by {} after {} rounds".format(what, PIECE_ROUNDS))

    def _move_phi(self, dC, tally, kind):
        """phi moved by the rows' volume changes dC [m3] (n,) and their carry, corrected per piece, beta derived and the
        volume moved to it (module docstring; planning ruling 3, steps 2-5); the carry is taken. `kind` labels the
        correction's entry in `tally.corrections` ("part", "carry"). Returns the enthalpy the volume's move removed
        [J], for the caller to book."""
        mesh, mk, sk = self.mesh, self.markers, self.skins
        phi0, beta0 = self.phi, self.beta
        dV = np.asarray(dC, dtype=float) + self.carry
        V0 = inside_volumes(mesh, phi0)
        moved = move_by_volumes(mesh, phi0, mk.x, mk.normal, dV[mk.skin_index], band=self.band)
        V1 = inside_volumes(mesh, moved.phi)
        change, swept = (V1 - V0).sum(), moved.swept.sum()
        if abs(swept - change) > MOVE_SUM_RTOL * V0.sum():
            raise RuntimeError("the groups of phi's move swept {:.6e} m3 but phi's volume changed by {:.6e} m3, more "
                               "than {:.0e} of the body's volume apart: cells no marker may own changed".format(
                                   swept, change, MOVE_SUM_RTOL))
        tally.move_miss = max(tally.move_miss, float(moved.miss.max(initial=0.0)))
        # a speck the move parted from the body (no marker on it and none of the volume's metal) leaves, its metal the
        # largest piece's (stage 3 Task 9's rule for used-up pieces of beta), its nodes out of the correction's claims
        phi, targets, label, _, lifted = self._lift_specks(mesh, moved.phi, beta0, phi0, tally, "a move")
        if targets[0] <= 0.0:
            raise RuntimeError("the gas surface's largest piece has no metal left after a move ({:.3e} m3)".format(
                targets[0]))
        misses = []
        phi, worst = correct_volumes(mesh, phi, label, targets, misses=misses, fixed=lifted)
        tally.correction = max(tally.correction, worst)
        tally.corrections.append({"kind": kind, "worst": float(worst), "body": float(misses[0])})
        self.phi = phi
        self.carry = np.zeros(sk.n)
        beta = skin_base(mesh, phi, self.settings.depth, band=self.band)
        removed = self.solver.move(beta, seed=lambda nodes: self._temperature_field(phi0)[nodes])
        tally.islands += getattr(self.solver, "seeded", 0)
        self.beta = beta
        return removed

    def _own(self):
        """The markers' areas, and so the skins' top areas (`_sync`), afresh on phi as it now is: the area of phi each
        marker owns (`Markers.owned_patches`, as `Markers.tessellate` gives them: the same-side and same-piece rules),
        the markers where they stand. Within a step phi moves by its columns' volumes while the markers stay where the
        step began (they settle at its end); the top fluxes follow phi (plan 03b Task 5: on the test cube at 6 MW/m2
        with the liquid fed off, the edges receded up to 1.9 mm within one 0.5 s step, and an edge marker's step-start
        area came to be 40 times its column's cross-section at phi, its heat poured into a column that was no longer
        there). At least AREA_FLOOR of the spacing squared each."""
        mk = self.markers
        area = np.asarray(mk.owned_patches(self.mesh, self.phi).sum(axis=1)).ravel()
        mk.area = np.maximum(area, AREA_FLOOR * mk.spacing ** 2)
        self._sync()

    def _division(self):
        """The layer between phi and beta divided among the markers as they stand (`division.divide`)."""
        mk = self.markers
        return divide(self.mesh, self.phi, self.settings.depth, mk.x, mk.normal)

    def _divide(self, tally, d=None):
        """The layer divided among the markers (`_division`; or `d`, a division of them as they stand on phi as it is,
        such as the step end's settle leaves) and the skins given its profile; the largest departure of a piece's
        closure factor from one into `tally.layer_factor`. RuntimeError, before the skins take the profile, where a
        piece of phi has layer but no column (no marker on it: its factor infinite) or a factor more than DIVISION_GROSS
        from one (module docstring; the controller's ruling T7-g)."""
        mk, sk = self.markers, self.skins
        d = self._division() if d is None else d
        if len(d.area) != mk.n:
            raise ValueError("the division must be of the markers as they stand: one column per marker")
        if len(d.factor):
            tally.layer_factor = max(tally.layer_factor, float(np.abs(d.factor - 1.0).max()))
            self._check_division(d)
        area = np.zeros((sk.n, len(d.levels)))
        area[mk.skin_index] = d.area
        sk.set_profile(d.levels, area)
        self.division = d

    def _check_division(self, d):
        """RuntimeError naming the first piece of phi that the division `d` left with layer but no column, or whose
        closure factor departs from one by more than DIVISION_GROSS (`_divide`)."""
        off = np.abs(d.factor - 1.0) > DIVISION_GROSS
        if not off.any():
            return
        p = int(np.flatnonzero(off)[0])
        label, _ = pieces(self.mesh, self.phi)
        layer = inside_volumes(self.mesh, self.phi) - inside_volumes(self.mesh, self.phi + self.settings.depth)
        cells = label == p
        at = (self.mesh.centroids[cells] * layer[cells][:, None]).sum(axis=0) / max(layer[cells].sum(), 1e-300) * 1e3
        where = "its layer {:.3e} m3 about ({:.3f}, {:.3f}, {:.3f}) mm".format(layer[cells].sum(), *at)
        if np.isinf(d.factor[p]):
            raise RuntimeError("piece {} of phi has layer but no column: no marker stands on it ({}; the controller's "
                               "ruling T7-g)".format(p, where))
        raise RuntimeError("the division's closure factor of piece {} of phi is {:.6g}, more than DIVISION_GROSS ({}) from "
                           "one, over its {} column(s) ({}; the controller's ruling T7-g)".format(
                               p, d.factor[p], DIVISION_GROSS, int((d.piece == p).sum()), where))

    def _baseless(self):
        """Each row's column has no base (n,) bool: its division base at most BASELESS_SHARE of its top."""
        d, mk = self.division, self.markers
        out = np.zeros(self.skins.n, dtype=bool)
        out[mk.skin_index] = d.base <= BASELESS_SHARE * d.area[:, 0]
        return out

    def _neighbours(self, row, base):
        """The neighbour a column with no base fills from (`_fill`): (its row, the depth [m] of the deepest level at
        which they share a boundary), preferring one with a base (`base` (n,) bool), then the deepest, then the
        longest there; None for a column that shares no boundary."""
        d, mk = self.division, self.markers
        m = self._row_marker()[row]
        e = np.flatnonzero((d.i == m) | (d.j == m))
        if not len(e):
            return None
        L = d.length[e]
        some = L > 0.0
        deepest = np.where(some.any(axis=1), L.shape[1] - 1 - np.argmax(some[:, ::-1], axis=1), -1)
        ok = deepest >= 0
        if not ok.any():
            return None
        e, deepest = e[ok], deepest[ok]
        other = np.where(d.i[e] == m, d.j[e], d.i[e])
        rows = mk.skin_index[other]
        best = np.lexsort((-d.length[e, deepest], -deepest, ~base[rows]))[0]
        return int(rows[best]), float(d.levels[deepest[best]])

    def _give_at_depth(self, row, z, m, H, P):
        """Metal (m [kg], H [J], P [kg m/s] (3,)) joins the cell of column `row` that spans the depth z [m] (its bottom
        cell below its cells, a new cell in an empty column), mixing with it: books exact."""
        sk = self.skins
        if sk.count[row] == 0:
            sk._append_bottom([row], m, H, P / m)
            return
        z_bot = np.cumsum(sk.dz()[row, :sk.count[row]])
        k = int(min(np.searchsorted(z_bot, z), sk.count[row] - 1))
        P0 = sk.m[row, k] * sk.U[row, k]
        sk.m[row, k] += m
        sk.H[row, k] += H
        sk.U[row, k] = (P0 + P) / sk.m[row, k]

    def _fill_baseless(self, rows, miss, base):
        """The columns `rows` with no base filled from their neighbours (module docstring): each column's miss [m3]
        (its volume less its cells', positive short) taken from, or given to, the neighbour at the depth of their
        deepest shared level. Returns each row's volume change no exchange with the volume made (n,) [m3] (the mixing
        the metal moved between columns makes; the carry's)."""
        sk = self.skins
        mixed = np.zeros(sk.n)
        for r, d in zip(rows, miss):
            near = self._neighbours(r, base)
            if near is None:
                continue
            j, z = near
            v0 = sk.cells_volume()
            if d > 0.0:                                     # short: the neighbour's metal just above that depth
                Vz = float(sk._volume_at(np.array([[z]]), [j])[0, 0])
                top = max(Vz - d, 0.0)
                za, zb = sk._depth_of(np.array([[top, top + d]]), [j])[0]
                mass, H, P = (a[0] for a in sk.take_between([j], za, zb))
                if mass <= 0.0:
                    continue
                sk._append_bottom([r], mass, H, P / mass)
            else:                                           # overfull: its bottom excess into the neighbour's cell
                V = float(v0[r])
                za = float(sk._depth_of(np.array([[max(V + d, 0.0)]]), [r])[0, 0])
                below = 2.0 * float(sk.dz()[r].sum()) + self.settings.depth       # beyond the cells' bottom
                mass, H, P = (a[0] for a in sk.take_between([r], za, below))
                if mass <= 0.0:
                    continue
                self._give_at_depth(j, z, mass, H, P)
            v1 = sk.cells_volume()
            mixed[j] += (v1[r] - v0[r]) + (v1[j] - v0[j])   # what the pair's cells gained between them (mixing)
        return mixed

    def _fill(self, tally):
        """Every column filled to its volume (module docstring; planning ruling 3, step 7): the columns with no base
        from their neighbours, then those with one at their bases, in at most FILL_ROUNDS rounds; each row's carry
        takes its cells' volume change less the mass the volume gave it over rho_s. A column with a base draws at its
        base temperature, its group's area mean of beta's patch temperatures (`_interface`, once; the controller's ruling
        T5-k); RuntimeError if a short one has no group, and RuntimeError naming the worst (its marker, position and
        piece, its cells' volume and its column's) where a column is still beyond FILL_RTOL of its volume after the
        rounds (the controller's ruling T7-g: a column with no base and no neighbour was skipped silently). Counts in
        `tally`: `drawn`, `handed` [kg] the base exchange each way, `to_volume` [kg m/s] the returned cells' momentum (the
        held volume takes it), `fill_miss` the largest miss left as a share of its column's volume, `baseless` the
        columns with no base. Returns the enthalpy [J] the skins took from the volume (negative where they returned
        more), for the caller to book as heat owed."""
        sk = self.skins
        V_col = sk.column_volume()
        base = ~self._baseless()
        tally.baseless = int((~base).sum())
        S, A_group, _ = self._interface()
        T_draw = np.full(sk.n, np.nan)                      # each column's base temperature (NaN: no group)
        has = A_group > 0.0
        T_draw[has] = (S @ self.solver.patch_T())[has]
        taken = 0.0
        for _ in range(FILL_ROUNDS):
            miss = V_col - sk.cells_volume()
            need = np.abs(miss) > FILL_RTOL * V_col
            if not need.any():
                break
            lone = np.flatnonzero(need & ~base)
            if len(lone):
                self.carry += self._fill_baseless(lone, miss[lone], base)
            miss = V_col - sk.cells_volume()
            rows = np.flatnonzero((np.abs(miss) > FILL_RTOL * V_col) & base)
            if not len(rows):
                continue
            v0 = sk.cells_volume()[rows]
            d = miss[rows]
            short = d > 0.0
            if np.isnan(T_draw[rows[short]]).any():
                raise RuntimeError("{} short column(s) with a base have no group of beta to draw from".format(
                    int(np.isnan(T_draw[rows[short]]).sum())))
            dV = d.copy()
            dV[short] = d[short] * sk.density(T_draw[rows[short]]) / sk.rho_s      # a draw fills d at T_draw's density
            ex = sk.exchange_base(rows, dV, T_draw[rows])
            self.carry[rows] += (sk.cells_volume()[rows] - v0) - ex.mass / sk.rho_s
            tally.drawn += float(ex.mass[ex.mass > 0.0].sum())
            tally.handed -= float(ex.mass[ex.mass < 0.0].sum())
            tally.to_volume -= ex.P.sum(axis=0)
            taken += float(ex.H.sum())
        held = sk.cells_volume()
        miss = np.abs(V_col - held)
        share = np.divide(miss, V_col, out=np.where(miss > 0.0, np.inf, 0.0), where=V_col > 0.0)
        tally.fill_miss = max(tally.fill_miss, float(share.max(initial=0.0)))
        left = miss > FILL_RTOL * V_col
        if left.any():
            r = int(np.argmax(np.where(left, share, -1.0)))
            m = int(self._row_marker()[r])
            raise RuntimeError("the fill left {} column(s) beyond FILL_RTOL ({:.0e}) of their volume after {} rounds; the "
                               "worst, marker {} at ({:.3f}, {:.3f}, {:.3f}) mm on piece {} ({}), its cells hold {:.6e} m3 "
                               "against its column's {:.6e} m3, {:.3e} m3 apart ({:.3e} of its volume; the controller's "
                               "ruling T7-g)".format(int(left.sum()), FILL_RTOL, FILL_ROUNDS, m, *(self.markers.x[m] * 1e3),
                                                     int(self.division.piece[m]), "with a base" if base[r] else "no base",
                                                     held[r], V_col[r], miss[r], share[r]))
        return taken

    def _layer(self, dC, tally, kind):
        """The layer's acceptance (planning ruling 3, steps 2-7): phi moved by the rows' volume changes dC [m3] and
        their carry (`_move_phi`), the layer divided and the columns filled; the heat the move and the fill leave owed to
        the volume booked at the bases, in proportion to each row's volume change."""
        dV = np.asarray(dC, dtype=float) + self.carry
        removed = self._move_phi(dC, tally, kind)
        self._own()
        self._divide(tally)
        owed = removed - self._fill(tally)
        tally.move += abs(owed)
        self._book_move(self._bases(), np.abs(dV), owed)

    # -- the step's end -----------------------------------------------------------------------------------------------
    def _end_step(self, tally, surface=None, check_core=True):
        """The macro step's end (module docstring): the last part's carry taken by phi (planning ruling 4), then
        `_rebuild` from the surface the markers' areas were last given on (`surface`, the samples of phi as it is by
        default: the parts give the areas afresh as phi moves, `_own`). The carry `_rebuild`'s settle and fill leave
        waits for the next step's first move (the controller's ruling T5-i)."""
        if self._needs_move(self.carry):
            weight = np.abs(self.carry)
            owed = self._move_phi(np.zeros(self.skins.n), tally, "carry")
            tally.move += abs(owed)
            self._book_move(self._bases(), weight, owed)
            self._own()
        return self._rebuild(tally, self._samples() if surface is None else surface, check_core=check_core)

    def _rebuild(self, tally, surface, check_core=False):
        """psi, the mesh rebuilt from the refinement targets with phi carried to it and the solve remeshed, the markers
        settled, the layer divided and filled, the thin-metal stop and psi again (module docstring). `surface`: the
        samples of the surface the markers' areas were last given on. With `check_core`, a psi with no rigid metal
        raises (`_check_core`). Returns whether the mesh stayed as it was."""
        mesh = self.mesh
        psi = self._psi(self.phi)
        if check_core:
            self._check_core(psi)
        target = self._target(mesh, self.phi, psi)
        same = not self._needs_rebuild(target, self._target(mesh, self.phi, psi, False, LEAD - REGION_SLACK *
                                                            self.fine_size))
        if not same:
            new = build(self.base, target)
            phi = self._reset(new, interpolate(mesh, self.phi, new.points))
            v_old = inside_volumes(mesh, self.phi)
            label, n = pieces(new, phi)
            if n == 0:
                raise RuntimeError("the gas surface has no metal on the rebuilt mesh")
            lifted = np.zeros(len(phi), dtype=bool)
            if n > 1:                                      # a speck the carried phi parted off leaves, as a move's does
                phi, _, label, n, lifted = self._lift_specks(new, phi, self.beta, self.phi, tally, "the remesh")
            targets = np.zeros(n)
            targets[0] = v_old.sum()
            if n > 1:                                      # the pieces share the old metal by their volumes here
                v_new = inside_volumes(new, phi)
                share = np.array([v_new[label == p].sum() for p in range(n)])
                targets = share * (v_old.sum() / share.sum())
            misses = []
            phi, worst = correct_volumes(new, phi, label, targets, misses=misses, fixed=lifted)
            tally.correction = max(tally.correction, worst)
            tally.corrections.append({"kind": "remesh", "worst": float(worst), "body": float(misses[0])})
            beta = skin_base(new, phi, self.settings.depth, band=self.band)
            E0 = self.solver.energy()
            self.solver.remesh(new, beta)
            self.mesh, self.phi, self.beta = new, phi, beta
            owed = E0 - self.solver.energy()
            tally.remesh += abs(owed)
            self._book_move(self._bases(), self.skins.area, owed)
        self._divide(tally, self._settle(tally, surface))
        owed = -self._fill(tally)
        tally.move += abs(owed)
        self._book_move(self._bases(), self.skins.area, owed)
        self._check_thin()
        self.psi = self._psi(self.phi)
        return same

    def _needs_rebuild(self, target, region):
        """Whether the mesh must be rebuilt (a rebuild builds every leaf to `target`): some leaf coarser than the fine
        region asks with its lead shortened by REGION_SLACK fine cells (`region`, no curvature), some leaf more than
        CURVATURE_SLACK times its full target, or more than COARSEN_SHARE of the metal (by volume, phi's) in leaves that
        could be COARSEN_SIZES times coarser. A rebuild at every step would not settle: phi carried to a new mesh is a
        slightly different contour, and its recovered curvature a slightly different field, which move leaves across the
        targets' thresholds and so change the next mesh (measured on the test cube: four rebuilds in a row, each with
        some 36 leaves over their curvature target)."""
        mesh = self.mesh
        size, x = mesh.size, mesh.centroids
        want = np.asarray(target(x, size), float)
        if np.any(size > np.asarray(region(x, size), float) * (1.0 + 1e-9)):
            return True
        if np.any(size > CURVATURE_SLACK * want):
            return True
        v = inside_volumes(mesh, self.phi)                  # the metal's leaves; closure grades the rest anyway
        coarser = COARSEN_SIZES * size <= want
        return float(v[coarser].sum() / v.sum()) > COARSEN_SHARE

    # -- the markers ----------------------------------------------------------------------------------------------------
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
        """Release the markers idx: each gives its skin to the nearest marker that stays and may own its part of phi
        (the same-side and same-piece rules on phi as it is), by depth (`Markers.release`, planning ruling 6), before
        any column's profile is replaced; the next fill settles both."""
        self.markers.release(idx, self.mesh, self.phi)

    def _samples(self):
        """phi's surface samples for the re-partition, SAMPLE_SHARE of the target spacing across (`surface_samples`)."""
        return surface_samples(self.mesh, self.phi, SAMPLE_SHARE * self.spacing)

    def _owners(self, surface):
        """The markers' old owners for a re-partition: [their positions and normals at the last giving of their areas
        (`Markers.x_tess`, `normal_tess`), the samples `surface` of the surface those areas were given on]."""
        mk = self.markers
        return [mk.x_tess.copy(), mk.normal_tess.copy(), surface]

    def _repartition(self, old, surface):
        """The skins follow the layer the markers' last giving of areas re-partitioned (module docstring): the layer
        divided among the markers as they stand, every old column (its volume its profile's) cut by the transfer across
        that division and given to its receivers by depth; the re-partition kept in `last_settle`. `old` from
        `_owners`, `surface` the samples of phi now. Returns the division (the skins keep their pieces' profiles until
        the caller gives them its own, `_divide`)."""
        mk, sk = self.markers, self.skins
        old_x, old_normal, old_surface = old
        d = self._division()
        old_volume = sk.column_volume()[mk.skin_index]
        T, stretch = transfer_matrix(old_x, old_normal, old_volume, old_surface, mk.x, d.volume, d, surface)
        apply_transfer(sk, mk.skin_index, T, stretch)
        self.last_settle.append(Repartition(T, stretch, old_volume, d.volume.copy()))
        return d

    def _drop_owners(self, old, idx):
        """`old` (from `_owners`) without the rows of the markers idx, which are about to be released."""
        stay = np.setdiff1d(np.arange(self.markers.n), idx)
        return [a[stay] for a in old[:2]] + [old[2]]

    def _tessellate(self, old, tally):
        """The markers given their areas of phi (`Markers.tessellate`: their patches, normals, curvatures, owned areas
        and longest directions); those it leaves no area (the patchless) released and the rest snapped and given their
        areas again until none is left, at most BARE_ROUNDS times, their rows leaving `old` (fix round 1 of Task 9b:
        after every giving of the step end, not only the first). Returns `old` (None stays None)."""
        mk = self.markers
        for _ in range(BARE_ROUNDS):
            mk.tessellate(self.mesh, self.phi)
            bare = np.flatnonzero(mk.area <= 0.0)
            if not len(bare):
                return old
            old = None if old is None else self._drop_owners(old, bare)
            self._release(bare, tally)
            self._snap(tally)
        raise RuntimeError("markers were still left no area after {} tessellations".format(BARE_ROUNDS))

    def _settle(self, tally, surface=None, repartition=True):
        """The markers onto phi, the trapped (and the patchless) released, given their areas, the layer divided and
        their skins re-partitioned across the division, split, merged, given their areas, divided and re-partitioned
        again, and the skins regridded with their markers' areas. `surface`: the samples of the surface the markers'
        areas were last given on (the step's start's phi), phi as it is if None. Without `repartition` (construction,
        whose skins are filled afresh after it) the skins follow the markers' changes alone and nothing is divided. The
        cells' volume change this makes (columns mixed at another temperature, a Scheil alloy's volumes not adding; no
        metal crosses beta) is carried, spread over the rows by their cells' volume (module docstring). Returns the
        division of the settled markers (None without `repartition`), for the caller to give the skins (`_divide`)."""
        mk, sk = self.markers, self.skins
        self.last_settle = []
        V0, C0 = float(sk.cells_volume().sum()), float(self.carry.sum())
        old = self._owners(self._samples() if surface is None else surface) if repartition else None
        self._snap(tally)
        cell, _ = locate(self.mesh, mk.x)
        trapped = mk.trapped_markers(self.mesh, self.phi, self.mesh.size[cell])
        if len(trapped):
            old = None if old is None else self._drop_owners(old, trapped)
            self._release(trapped, tally)
            self._snap(tally)
        old = self._tessellate(old, tally)
        d = None
        if repartition:
            now = self._samples()
            d = self._repartition(old, now)
        self._sync()
        big = mk.split(self.spacing ** 2)
        merged = mk.merge(self.spacing)
        if len(big) or merged:
            old = self._owners(now) if repartition else None
            self._snap(tally)
            old = self._tessellate(old, tally)
            d = self._repartition(old, self._samples()) if repartition else None
        self._sync()
        self.skins.regrid()
        self._sync()
        v = sk.cells_volume()
        self.carry = np.zeros(sk.n)
        if v.sum() > 0.0:                                  # the carry before it and the cells' change, spread
            self.carry += (C0 + float(v.sum()) - V0) * v / v.sum()
        return d
