"""The body the trajectory carries, as a stepper (spec section 4). ConstantBody keeps mass and one temperature
(Step 1 behaviour, `--thermal none`); ThermalBody wraps a mesh, a material and a conduction solver and advances the
temperature field with each macro step's per-patch convective loads, keeping the heat bookkeeping SESAM reports."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def reference_area(self): ...                        # m2 drag reference area, or None for the fixed pi D^2/4
    def reference_length(self): ...                      # m length of the body Knudsen number, or None for the fixed D
    def drag_shape_factor(self): ...                     # continuum C_D of the current shape / a sphere's (1 for a sphere)
    def surface_temperature(self) -> np.ndarray: ...     # K per patch
    def mean_temperature(self) -> float: ...             # K, energy-equivalent (spec 6.4)
    def energy(self) -> float: ...                       # J stored above the material's reference temperature
    def field(self): ...                                 # nodal temperatures, or None


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def reference_length(self):
        return None

    def drag_shape_factor(self):
        return 1.0

    def advance(self, t, dt, loads):
        return None

    def surface_temperature(self):
        return np.array([self.temperature_K])

    def mean_temperature(self):
        return self.temperature_K

    def energy(self):
        return 0.0

    def field(self):
        return None


class ThermalBody:
    """Mesh + material + ThermalSolver. `v_hat` is the direction of motion in the body frame (fixed attitude: the
    stagnation patch is the one whose normal is along v_hat, spec 13.1)."""

    def __init__(self, mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        self.mesh, self.material, self.solver, self.mass_kg = mesh, material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient = float(T_ambient)
        self.surface = mesh.surface()
        self.theta = self.surface.angles_to(v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(v_hat), self.surface.patch_toward(-np.asarray(v_hat, dtype=float))
        self.i_centre = mesh.centre_node()
        self.volume = mesh.volume()
        self.solver.setup(mesh, material, self.emissivity)
        self.solver.set_temperature(T0)
        self.energy0 = self.solver.energy()
        self.integrated_heat = 0.0       # J, integral of Q_conv dt (SESAM's integrated_heat_J)
        self.absorbed_heat = 0.0         # J, integral of (Q_conv - Q_rad) dt
        self.radiated_heat = 0.0         # J
        self.iterations = []             # Newton iterations per step
        self.last = None                 # thermal.StepResult of the last step

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def reference_length(self):
        """Length scale of the body Knudsen number, or None for the fixed initial diameter."""
        return None

    def drag_shape_factor(self):
        """Continuum drag of the body's shape relative to a sphere's: 1 for the sphere of Steps 1-2."""
        return 1.0

    def nose_radius(self):
        """Stagnation-point radius of curvature [m] for the heating and the surface flow: the sphere's."""
        return float(self.mesh.params.get("radius_m", 0.05))

    def advance(self, t, dt, loads, state=None):
        res = self.solver.step(dt, loads.q_conv, self.T_ambient)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res

    def surface_temperature(self):
        return self.surface.facet_mean(self.solver.temperature())

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / (self.material.rho * self.volume)))

    def energy(self):
        return self.solver.energy()

    def field(self):
        return self.solver.temperature()

    def radiated_power(self):
        return self.solver.radiated_power(self.T_ambient)

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat) / absorbed heat: zero for an exact discrete balance."""
        return (self.energy() - self.energy0 - self.absorbed_heat) / self.absorbed_heat if self.absorbed_heat else 0.0

    def surface_stats(self):
        T, Tf = self.solver.temperature(), self.surface_temperature()
        return {"surface_T_max_K": float(Tf.max()), "surface_T_min_K": float(Tf.min()),
                "surface_T_mean_K": float((Tf * self.surface.areas).sum() / self.surface.area),
                "T_stagnation_K": float(Tf[self.i_stag]), "T_back_K": float(Tf[self.i_back]), "T_centre_K": float(T[self.i_centre])}


SIZE_FEEDBACK_NAMES = ("current", "initial")
FEED_DEPTH_NAMES = ("conjugate", "all")     # how deep melt may leave an element for the film (Step 3, 2026-09-22)
CP_MAX_NEWTONIAN = 1.84          # only the ratio to the meshed sphere's own value is used, so this cancels
NOSE_CAP_ANGLE = 30.0            # deg: the windward cap fitted for the nose radius (depth (1 - cos 30 deg) R_t behind the front)
NOSE_CAP_FACTOR = 1.67           # a flat face of radius R_t heats like a sphere of 1.67 R_t (its stagnation velocity gradient is
                                 # ~0.6 x a sphere's of the same radius, Boison & Curtiss 1959): the cap on the fitted radius
PHI_MIN = 1.0e-3                 # element fraction kept by elements that do not own a patch (they cannot die: no cavities)
NEAREST_PATCHES = 4              # patches that receive an interior element's liquid (area-weighted)
LOAD_DT_MAX = 1000.0             # K, the most a deferred melt load may move one node in one macro step
PHI_DEATH = 0.05                 # a patch owner below this fraction dies (its remainder goes to the film): keeps the surface
                                 # nodes' thermal mass above 5 % of an element's, which the Newton iteration needs (measured
                                 # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
MAX_CASCADE_PASSES = 32          # safety cap on the molten cascade's passes in one macro step; what it leaves waits for the
                                 # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
NONRIGID_MAX_CROSSINGS = 64      # elements the non-rigid depth's line may cross before it stops (amendment of 2026-10-06)
REMOVAL_NAMES = ("girin", "instant")


@dataclass
class MeltSettings:
    removal: str = "girin"           # girin: film + runoff + spraying; instant: liquid removed as it forms (verification device)
    runoff: bool = True
    demise_fraction: float = 0.01    # the run ends when the mass falls below this fraction of the initial mass
    particles: bool = True           # keep the source-table rows
    size_feedback: str = "current"    # current: Kn on the equivalent diameter of the remaining mass and the nose radius fitted
                                      # to the windward cap; initial: D0 and R0 throughout (SESAM's convention, for the devices)
    feed_depth: str = "conjugate"     # conjugate: only liquid the gas shear reaches leaves its element (fact 28); all: any depth
    deep_runoff: bool = True          # the liquid below the conjugate depth runs off, never sprayed (amendment of 2026-10-02)
    molten_cascade: bool = True       # a death that exposes a fully molten element (mean T >= T_feed) feeds it within the
                                      # step, so the surface recedes through molten material by more than one element per
                                      # macro step (amendment of 2026-10-03)
    rigid_substrate: bool = True      # the thin branch needs a rigid substrate: the regime test reads the film plus the slurry
                                      # (above 50 % liquid) beneath it, and off Girin's closure a film on slurry does not
                                      # spray (amendment of 2026-10-06)

    def __post_init__(self):
        if self.removal not in REMOVAL_NAMES:
            raise ValueError("removal must be one of {}, got {!r}".format(REMOVAL_NAMES, self.removal))
        if self.size_feedback not in SIZE_FEEDBACK_NAMES:
            raise ValueError("size_feedback must be one of {}, got {!r}".format(SIZE_FEEDBACK_NAMES, self.size_feedback))
        if self.feed_depth not in FEED_DEPTH_NAMES:
            raise ValueError("feed_depth must be one of {}, got {!r}".format(FEED_DEPTH_NAMES, self.feed_depth))
        if not 0.0 < self.demise_fraction < 1.0:
            raise ValueError("demise_fraction must be within (0, 1)")


def fit_sphere(points):
    """Algebraic least-squares sphere through `points` (n, 3): (centre, radius)."""
    x = np.asarray(points, dtype=float)
    A = np.column_stack([2.0 * x, np.ones(len(x))])
    p = np.linalg.lstsq(A, (x * x).sum(axis=1), rcond=None)[0]
    c = p[:3]
    return c, float(math.sqrt(max(p[3] + c @ c, 0.0)))


class MeltingBody(ThermalBody):
    """ThermalBody with melting (spec Step 3 sections 6, 8-11): element fractions phi_e, the film per patch, the
    gas-side surface flow, spraying, element death and the mass/area/energy accounting.

    Per macro step (`advance`): the conduction step with the deferred melt loads of the previous step -> melt step:
    (i) every element's liquid inventory f_feed(T_e) phi_e rho V_e becomes film on its patches (owners) or the
    nearest patch (interior elements keep PHI_MIN so no cavity can open), netted against (iv) below; what leaves is
    the element's molten part, at the enthalpy that part carries, so the melt front moves at the energy-limited rate
    whatever the element size; (ii-a) the deep runoff (`_deep_runoff`, amendment of 2026-10-02): the contiguous liquid
    below the conjugate depth, which (i) holds in its elements, moves along the surface under the pressure gradient and
    the deceleration into a per-patch deep account `m_d` that is never sprayed and becomes film only from the top, as
    far as the film is thinner than the conjugate depth; (ii) surface flow, delta_m, lubrication, runoff transport;
    (iii) spraying and release;
    (iv) re-solidification: the fraction 1 - f_feed(T_patch) of each patch's film returns to its owner element
    (raising phi_e, capped at 1) once the patch falls back through the ramp -- the mirror of the feed rule, netted
    with it so that no element both melts and freezes in one step; (v) patch owners at phi <= PHI_DEATH die: the
    mesh's active set, the surface, the film (handed to the nearest surviving patch) and the solver's fractions are
    refreshed; with the molten cascade (`_feed_exposed`, amendment of 2026-10-03) an element a death exposes fully
    molten (its mean temperature at or above T_feed, `molten_depth`'s test) is fed whole within the step, so that it
    dies in turn and the surface recedes through molten material by as many elements as are molten rather than by one
    per macro step. With `removal = "instant"`
    the liquid leaves the body at h_liquid instead, and the difference to the element's own enthalpy is a nodal load on
    its nodes over the next step.

    The film has the surface's own temperature, because it is thermally thin: q b / k_l = 0.22 K across a 10 um film
    at 2 MW/m2 (22 K even at 1 mm) and b^2/alpha = 0.3 ms against the 0.5 s macro step. Rather than give it an energy
    equation, its *mass* is handed to the solver (`set_film_mass`), which weighs it with the liquid heat capacity, so
    the film rides the boundary nodes at the surface temperature by construction; it heats, cools and re-solidifies
    with the surface, and the droplets carry away whatever superheat the surface has. A film is liquid by
    construction, so it holds the *liquid* enthalpy h(T) + L_f (1 - f_l) wherever it sits (`film_enthalpy`): feeding
    it costs the latent heat its mass has not yet paid, which is what keeps a body that merely sits on the melting
    ramp from turning into film for free. `film_energy` sums that enthalpy over the patches -- identical to the nodal
    sum the solver's capacity matrix carries -- and every transfer -- feed, freeze-back, runoff between patches at
    different temperatures, the hand-over when a patch dies -- books the enthalpy difference it carries, at the
    patch it arrives on (`melt_step`), so the coupled balance stays exact (measured 1e-10). `removal = "instant"`
    removes the liquid inventory of every element as it forms, without film, runoff or spraying: the lumped-melting
    device that reproduces SESAM's Q/L_f law (SESAM hollows the sphere at fixed outer geometry, measured
    2026-09-20, so the geometry is kept until elements die).
    mass(t) = sum phi rho V + sum m_f; the drag reference area is the current surface's projection on the flight
    direction (pi R^2 while intact). Size feedback (decided 2026-09-21, replacing spec section 10's fixed R0): the body
    Knudsen number uses the equivalent diameter of the remaining mass, and the stagnation radius for the heating and
    the surface flow is fitted to the current windward cap -- a least-squares sphere through the patch centroids within
    (1 - cos 30 deg) R_t of the front-most point (R_t the transverse radius about the mass centre), bounded to
    [0.1 R_t, NOSE_CAP_FACTOR R_t] (the front erodes fastest and flattens, which lowers the stagnation heating as R^-1/2;
    a flat face heats like a sphere of 1.67 x its radius). The sphere drag tables are kept. `size_feedback = "initial"`
    keeps D0 and R0 (SESAM's convention, used by the verification devices)."""

    def __init__(self, mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0=300.0,
                 emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        if not material.melts or material.liquid is None:
            raise ValueError("MeltingBody needs a material with a latent heat and liquid properties (AA7075, AA7075_range or AA7075_scheil)")
        super().__init__(mesh, material, solver, mass_kg, T0, emissivity, T_ambient, v_hat)
        from . import film as film_mod, spray as spray_mod, surface_flow
        self.settings = settings or MeltSettings()
        self.flow = flow or surface_flow.SurfaceFlow()
        self.spray = spray_model or spray_mod.SprayModel(material.liquid)
        self._film_mod = film_mod
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.liquid = material.liquid
        self.element_mass = material.rho * mesh.element_volumes()          # kg at phi = 1
        self.phi = np.ones(mesh.n_elements)
        self.mass0 = float(self.element_mass.sum())
        self.m_f = np.zeros(self.surface.n_patches)
        self.m_d = np.zeros(self.surface.n_patches)                         # kg, mobile liquid below the conjugate depth
        self.pending_load = np.zeros(mesh.n_nodes)                          # J, deferred melt energy for the next step
        self.sprayed_mass = self.runoff_mass = self.removed_mass = self.removed_enthalpy = self.frozen_mass = 0.0
        self.deep_runoff_mass = 0.0                                         # kg the deep runoff took out of the elements
        self.deep_surfaced_mass = 0.0                                       # kg of deep liquid that became film from the top
        self.cascade_mass = 0.0                                             # kg the molten cascade fed (2026-10-03)
        self.cascade_capped_steps = 0                                       # steps whose cascade stopped at MAX_CASCADE_PASSES
        self.n_released = 0.0
        self.source_rows = []
        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
        self.last_flow = self.last_spray = None
        self.last_delta_m = None                                           # the step's conjugate depth per patch [m]
        self.last_nonrigid = None                                          # the step's non-rigid depth per patch [m]
        self.last_face_ids = None                                          # face ids of the surface they were evaluated on
        self.last_melt = {"n_dead": 0, "runoff_substeps": 0, "released_mass": 0.0, "n_released": 0.0, "feed_mass": 0.0}
        self.melt_onset = self.spray_onset = None
        self.consumed = False
        self._refresh_geometry()
        self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation

    # -- geometry -----------------------------------------------------------------------------------------------
    def _refresh_geometry(self):
        self.surface = self.mesh.surface()
        self.theta = self.surface.angles_to(self.v_hat)
        self.t_hat = self.surface.tangent_from(self.v_hat)
        self.windward = self.theta <= 0.5 * np.pi
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)
        self.runoff = self._film_mod.Runoff(self.surface, self.mesh.points, self.windward) if self.settings.runoff else None
        owners = self.surface.owner
        self.owner_area = np.bincount(owners, self.surface.areas, self.mesh.n_elements)      # patch area per owner element
        from scipy.spatial import cKDTree
        self._patch_tree = cKDTree(self.surface.centroids)
        self.patch_of_face = np.full(len(self.mesh._face_nodes), -1, dtype=np.int64)     # face id -> patch index (-1: not a patch)
        self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
        self._fit_nose()

    def _fit_nose(self):
        """Mass centre, transverse radius and the windward-cap radius of the current surface (module docstring)."""
        m, s, v = self.mesh, self.surface, self.v_hat
        w = self.phi * m.element_volumes()
        self.mass_centre = (w[:, None] * m.points[m.tets].mean(axis=1)).sum(axis=0) / max(w.sum(), 1e-300)
        X = m.points[np.unique(s.faces)] - self.mass_centre
        self.transverse_radius = float(np.sqrt(np.maximum(np.linalg.norm(X, axis=1) ** 2 - (X @ v) ** 2, 0.0)).max())
        x = (s.centroids - self.mass_centre) @ v
        band = x >= x.max() - (1.0 - math.cos(math.radians(NOSE_CAP_ANGLE))) * self.transverse_radius
        r0 = float(m.params.get("radius_m", 0.05))
        if band.sum() >= 4 and self.transverse_radius > 0.0:
            _, r = fit_sphere(s.centroids[band])
            self.fitted_nose_radius = r
            self.cap_nose_radius = float(min(max(r, 0.1 * self.transverse_radius), NOSE_CAP_FACTOR * self.transverse_radius))
        else:
            self.fitted_nose_radius = self.cap_nose_radius = r0

    def nose_radius(self):
        if self.settings.size_feedback == "initial":
            return float(self.mesh.params.get("radius_m", 0.05))
        return self.cap_nose_radius

    def newtonian_drag(self):
        """Modified-Newtonian drag coefficient of the current windward silhouette: sum C_p,max (n.v)^3 A / sum (n.v) A
        over the convex hull of the surface (0.92 for a sphere, 1.84 for a flat face, with C_p,max = 1.84).

        The hull, not the raw facets: the staircase left by element death scatters the normals and biases the integral
        low and non-monotonically (0.79-1.07 against the hull's smooth 0.92-1.73, measured 2026-09-21), and the flow
        sees the silhouette -- a shallow cavity recovers roughly the stagnation pressure at its mouth. It is therefore
        an upper bound where the face is cratered."""
        from scipy.spatial import ConvexHull
        points = self.mesh.points[np.unique(self.surface.faces)]
        try:
            hull = ConvexHull(points)
        except Exception:                                        # degenerate remnant: fall back to the sphere
            return CP_MAX_NEWTONIAN / 2.0, self.surface.projected_area(self.v_hat)
        tri = points[hull.simplices]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        area = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * np.maximum(area, 1e-300))[:, None]
        inward = np.einsum("ij,ij->i", n, tri.mean(axis=1) - points.mean(axis=0)) < 0.0
        n[inward] *= -1.0
        c = n @ self.v_hat
        w = c > 0.0
        projected = float((area[w] * c[w]).sum())
        if projected <= 0.0:
            return CP_MAX_NEWTONIAN / 2.0, self.surface.projected_area(self.v_hat)
        return float((CP_MAX_NEWTONIAN * c[w] ** 3 * area[w]).sum() / projected), projected

    def drag_shape_factor(self):
        """The continuum drag of the current shape relative to a sphere's, for aero.drag_coefficient: exactly 1 while
        the body is intact (it is normalised by the meshed sphere's own Newtonian value, so the mesh discretisation
        error cancels and the ATDB sphere table is reproduced bit for bit) and 2.00 at the flat-face limit, where the
        ATDB disc entry sits. `size_feedback = "initial"` pins it to 1 (SESAM's convention)."""
        if self.settings.size_feedback == "initial":
            return 1.0
        return self.newtonian_drag()[0] / self.newtonian_drag_sphere

    def reference_length(self):
        if self.settings.size_feedback == "initial":
            return None
        return 2.0 * self.equivalent_radius()

    def mass(self, t):
        return float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())

    def reference_area(self):
        """Drag reference area: the hull's projected area, paired with the C_D of drag_shape_factor (the raw surface's
        projection counts forward-facing patches inside a crater and is 2-4 % larger)."""
        if self.settings.size_feedback == "initial":
            return self.surface.projected_area(self.v_hat)
        return self.newtonian_drag()[1]

    def equivalent_radius(self):
        return (3.0 * self.mass(0.0) / (4.0 * math.pi * self.material.rho)) ** (1.0 / 3.0)

    def energy(self):
        return self.solver.energy() + self.film_energy()

    def mean_temperature(self):
        m = self.mass(0.0)
        return float(self.material.temperature_from_enthalpy(self.energy() / m)) if m > 0.0 else 0.0

    def melt_front_depth(self):
        """Deepest molten point [m]: the distance from the original surface of the innermost element with f_l > 0."""
        f_l = self.material.liquid_fraction(self.solver.temperature())[self.mesh.tets].max(axis=1)
        molten = self.mesh.active & (f_l > 0.0)
        if not molten.any():
            return 0.0
        r = np.linalg.norm(self.mesh.points[self.mesh.tets[molten]].mean(axis=1), axis=1)
        return float(self.mesh.params.get("radius_m", r.max()) - r.min())

    # -- the step --------------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
        # A deferred melt load is energy the transferred mass delivered to the nodes it arrived at; a node that has
        # since melted away has nothing to heat with it, so no node is asked to move more than LOAD_DT_MAX in one
        # step and the remainder waits (it is applied as soon as the node has the capacity, dropped into Q_dropped
        # if the node dies first, and counted either way -- the balance sees `pending_load` and `Q_dropped` alike).
        room = LOAD_DT_MAX * self.solver.nodal_capacity()
        applied = np.clip(self.pending_load, -room, room)
        self.pending_load = self.pending_load - applied
        res = self.solver.step(dt, loads.q_conv, self.T_ambient, nodal_load=applied / dt)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res
        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
        self.melt_step(t, dt, state)

    def melt_step(self, t, dt, state):
        """Feed, film transport, spraying, re-solidification and element death (class docstring).

        Energy is moved between three places -- the finite-element solid, the film, and the outside world -- and every
        transfer is booked by the same rule: mass that moves at one temperature books nothing, and mass that arrives
        somewhere colder or hotter than it left books the difference it carries, as a deferred load on the nodes it
        arrives at. So a transfer of m kilograms from enthalpy h_src to h_dst changes the accounted energy by
        m (h_dst - h_src) and applies m (h_src - h_dst) to the destination: the balance closes exactly, and nothing
        larger than the difference ever touches the temperature field. Booking the two halves separately -- the
        solid's m h_src on its element's nodes, the film's m h_dst on its patch's nodes -- closes the balance just as
        exactly and wrecks the field: the same mass is spread by 1/4 over four nodes on one side and by 1/3 over
        three on the other, which left +-m h/12 on every surface node and swung the body to 1618 K and -1202 K in two
        macro steps (measured 2026-09-22). The enthalpies are the ones the two energy functionals actually use: the
        element's is the mean of h(T_i) over its four nodes, the film's the mean of the *liquid* h over its patch's
        three, and what leaves an element is its molten part, at the enthalpy that part carries."""
        mat, liq, s = self.material, self.liquid, self.settings
        from . import spray as spray_mod
        T = self.solver.temperature()
        tets = self.mesh.tets
        # The surface flow is evaluated once per step, here rather than inside the film step, because the feed needs
        # Girin's conjugate melt-layer depth: the gas shear penetrates the liquid only that far, so only liquid within
        # that depth of the wall can be carried away, and material deeper than it keeps its melt until the surface has
        # receded to it. Feeding from any depth let the interior drain through an intact skin -- by 82 s of the 100 mm
        # flight the core was being consumed faster than the outermost shell (measured 2026-09-22).
        flow = delta_m = None
        if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
            delta_m, _ = spray_mod.melt_layer(flow, liq)
            self.last_face_ids = self.surface.face_ids.copy()      # the surface delta_m belongs to (the deaths come later)
        self.last_delta_m = delta_m                                # per patch, nan off Girin's closure; for the VTK frame
        h_node = mat.enthalpy(T)
        h_e = h_node[tets].mean(axis=1)                          # as solver.energy() weighs the solid
        h_liq = mat.enthalpy_liquid(T)
        h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
        # (i) feed and (iv) re-solidification, netted. The feed hands the molten fraction of what each element still
        # holds to the film (f_feed of phi_e, so an element that is a quarter molten hands over a quarter of its
        # remainder); the mirror rule hands back the fraction of the film that has fallen below the ramp. What stops
        # the feed from eating a body that merely sits on the ramp is not the rule but the energy: film is liquid and
        # carries L_f wherever it sits (`Material.enthalpy_liquid`), so every kilogram fed debits the latent heat it
        # has not paid and the surface falls back to the solidus. Book the film at the mixture enthalpy instead and
        # melting is free on the ramp -- a 100 mm sphere turned entirely to film on a quarter of its latent heat
        # (measured 2026-09-22). The two directions are netted per element because they are one equilibrium seen
        # from opposite sides: run separately they cycled 3 % of the body's mass through the film every step with no
        # net effect, pinning the surface at T_feed and paying Newton iterations for it.
        fn = mat.feed_fraction(T)[tets]
        f = fn.mean(axis=1) * self.mesh.active
        reach = self._shear_reaches(delta_m) if s.feed_depth == "conjugate" and s.removal != "instant" else None
        if reach is not None:
            f = f * reach
        h_hot = (fn * h_node[tets]).mean(axis=1) * self.mesh.active     # what the molten part carries, per kg of element
        cap = np.where(self.owner_area > 0.0, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
        gross = np.minimum(f * self.phi, cap) * self.element_mass
        # the liquid on each patch below the ramp re-solidifies: the film and the deep liquid beneath it alike
        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * (self.m_f + self.m_d) if s.removal != "instant" else np.zeros(0)
        want = np.bincount(self.surface.owner, solid, self.mesh.n_elements) if solid.size else np.zeros(self.mesh.n_elements)
        net = gross - want
        fed, want = np.maximum(net, 0.0), np.maximum(-net, 0.0)
        # the enthalpy the fed mass carries: the *molten* part of the element, not its mean. The mass that leaves sits
        # at the nodes that are above the ramp, so it carries h(T_i) there -- and what stays behind is the colder
        # rest, which is why the element must be debited the difference. Debit only the destination and the element
        # keeps a melt fraction it no longer has and feeds it again next step: the surface then melted three times
        # faster than the heat allowed and the body was gone in ten steps (measured 2026-09-22).
        with np.errstate(divide="ignore", invalid="ignore"):
            fed_h = np.where(f > 0.0, fed * h_hot / np.where(f > 0.0, f, 1.0), 0.0)
        self.phi = self.phi - fed / self.element_mass
        feed_mass = float(fed.sum())
        if feed_mass > 0.0 and self.melt_onset is None:
            self.melt_onset = t
        deep_liquid = 0.0
        if s.removal == "instant":
            self.removed_mass += feed_mass
            self.removed_enthalpy += feed_mass * mat.h_liquid
            self._defer_to_elements(fed * (h_e - mat.h_liquid))  # it leaves the body at the liquidus, not at h_e
            released = feed_mass
        else:
            self._defer_to_elements(fed * h_e - fed_h)           # the element keeps only what the melt left behind
            delta, carried = self._add_to_film(fed, fed_h)
            self._defer_to_patches(carried - delta * h_p)        # ... and the melt arrives at its patch's temperature
            if s.deep_runoff and flow is not None:
                deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
            released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
        frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
        # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed --
        # and, with the molten cascade (amendment of 2026-10-03), are fed within the step if the death exposed them
        # fully molten (`_feed_exposed`), so that they die in the next pass and expose the next
        n_dead, passes, cascade, capped = 0, 0, 0.0, False
        # fully molten as `molten_depth` has it -- the element's mean nodal temperature at or above T_feed -- so the cascade
        # removes exactly the contiguous molten layer the branch test reads, wherever a death exposes it
        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" else None
        while not self.consumed:
            dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
            if dead.size == 0:
                break
            rest = self.phi[dead] * self.element_mass[dead]
            if s.removal == "instant":
                self.removed_mass += float(rest.sum())
                self.removed_enthalpy += float(rest.sum()) * mat.h_liquid
                self._defer_to_elements(rest * (h_e[dead] - mat.h_liquid), dead)
            else:
                extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
                extra[dead], extra_h[dead] = rest, rest * h_e[dead]
                delta, carried = self._add_to_film(extra, extra_h)
                self._defer_to_patches(carried - delta * h_p)
            self.phi[dead] = 0.0
            before = float(((self.m_f + self.m_d) * h_p).sum())
            was_owner = self.owner_area > 0.0
            self._kill(dead)                                     # the liquid of a dead patch moves to another patch...
            if self.m_f.size:                                    # ... which is at its own temperature
                h_p = self.film_enthalpy(h_liq)
                self._spread_to_patches(before - float(((self.m_f + self.m_d) * h_p).sum()), self.m_f + self.m_d)
            n_dead += int(dead.size)
            if full is not None and not self.consumed:          # the molten cascade: walls the deaths exposed fully molten
                exposed = np.flatnonzero(full & (self.owner_area > 0.0) & ~was_owner & (self.phi > PHI_DEATH))
                if exposed.size and passes < MAX_CASCADE_PASSES:
                    cascade += self._feed_exposed(exposed, h_e, h_p)
                    passes += 1
                elif exposed.size:                               # the cap: these wait for the next step's feed
                    capped = True
        self.solver.set_fractions(self.phi)
        self.solver.set_film_mass(self._film_nodal())
        self.frozen_mass += frozen
        self.cascade_mass += cascade
        self.cascade_capped_steps += int(capped)
        self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass, "frozen_mass": frozen,
                               "deep_liquid": deep_liquid, "cascade_mass": cascade, "cascade_passes": passes,
                               "cascade_capped": capped})

    def _freeze_back(self, want, solid, h_e, h_p):
        """Give each element back the `want` kilograms of film that have fallen below the feed ramp, drawn from its
        own patches in proportion to what each wants to give (`solid`). An element can only recover film that is
        still there -- what ran off or sprayed away is gone -- and phi_e is capped at 1, because a patch may hold
        what its neighbours' runoff delivered and the mesh cannot grow a layer outside itself, so film with nowhere
        to go simply stays film (`film_frozen_fraction` records how much). The deep liquid beneath the film (amendment
        of 2026-10-02) lies against the solid, so it is drawn first; both hold the same liquid enthalpy on the same
        nodes, so which of the two gives up the mass books nothing."""
        if not self.m_f.size or not want.any():
            return 0.0
        owner = self.surface.owner
        have = np.bincount(owner, solid, self.mesh.n_elements)
        room = np.maximum(1.0 - self.phi, 0.0) * self.element_mass
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(have > 0.0, np.minimum(want, room) / np.where(have > 0.0, have, 1.0), 0.0)
        taken = np.minimum(solid * share[owner], self.m_f + self.m_d)
        if not taken.any():
            return 0.0
        per_element = np.bincount(owner, taken, self.mesh.n_elements)
        self.phi = np.clip(self.phi + per_element / self.element_mass, 0.0, 1.0)
        from_deep = np.minimum(self.m_d, taken)
        self.m_d = self.m_d - from_deep
        self.m_f = self.m_f - (taken - from_deep)
        # the film arrives in the element it froze onto, which is colder than the surface it left
        self._defer_to_elements(np.bincount(owner, taken * h_p, self.mesh.n_elements) - per_element * h_e)
        return float(taken.sum())

    def _feed_exposed(self, exposed, h_e, h_p):
        """(v-a) The molten cascade (amendment of 2026-10-03): feed the fully molten elements a death has just exposed,
        within the step. Returns the mass fed [kg].

        The melt step runs after the conduction, so between two melt steps the heat flux melts the material beneath
        and beside the wall-owning element as well as the element itself. The feed hands the owner's liquid to the film
        and the owner dies once consumed, but the element its death exposes was buried a moment earlier, so the feed
        gate of fact 28(b) held its melt in place, and it used to wait for the next step's feed: the surface receded
        through molten material by one element per macro step, and whatever melted below the owner waited as a backlog
        that shrank with the step (plan fact 50). Exposed, the element is a wall owner -- in the sheared layer by
        definition, which is what the gate admits -- so it is fed now.

        "Fully molten" is `molten_depth`'s definition, the element's mean nodal temperature at or above T_feed, so the
        cascade removes exactly the contiguous molten layer the branch test reads. It is not "every node above T_feed":
        the feed's own energy debit keeps every element it is feeding on the feed ramp, so at the melt front no element
        has all four nodes above it -- on the 100 mm flight not one of the 1700-1900 elements of the chains below the
        wall-owning elements does, their coldest node sitting at a median 897 K, 93 % liquid by the enthalpy, and a
        cascade on that test fed 32 g in 120 s (deep runoff off) and left the backlog as it was (measured 2026-10-02).
        The element is fed whole: its remainder phi rho V leaves at its mean nodal enthalpy h_e and nothing is booked on
        it, and it arrives on the faces it owns at their liquid enthalpy with the difference released there -- the owner
        rule of step (i) where every node is above T_feed (the molten part is then all of it), and the death rule of fact
        4 for the rest, whose residual latent heat the faces it lands on pay (facts 5 and 25). Left at phi = 0, it dies
        in the next pass of the death loop and hands that liquid down to the faces it exposes, like any dying patch's
        film; the cascade ends where the elements a death exposes are no longer fully molten -- the contiguity rule of
        `molten_depth` -- so molten material that cooler material separates from the wall is never fed. Only elements
        that would survive the pass are fed: one at or below PHI_DEATH dies anyway and leaves its remainder by the death
        rule, as before."""
        fed = np.zeros(self.mesh.n_elements)
        fed[exposed] = self.phi[exposed] * self.element_mass[exposed]
        self.phi[exposed] = 0.0
        delta, carried = self._add_to_film(fed, fed * h_e)     # the whole element, at its mean nodal enthalpy
        self._defer_to_patches(carried - delta * h_p)          # ... and the melt arrives at its patch's temperature
        return float(fed.sum())

    def _defer_to_elements(self, energy, elements=None):
        """Book a per-element energy [J] on the elements' nodes, to be applied over the next step."""
        e = np.asarray(energy, dtype=float)
        if not e.any():
            return
        nodes = self.mesh.tets if elements is None else self.mesh.tets[elements]
        np.add.at(self.pending_load, nodes.ravel(), np.repeat(e / 4.0, 4))

    def _spread_to_patches(self, energy, weights):
        """Book one energy [J] over the patches that received mass, in proportion to how much each received. Used
        where the transfer's per-edge detail is not carried back (runoff, and the hand-over when a patch dies): the
        total is exact and it lands on the arriving liquid, which is where the difference is released."""
        total = float(weights.sum())
        if total > 0.0 and energy != 0.0:
            self._defer_to_patches(energy * weights / total)

    def _defer_to_patches(self, energy):
        """Book a per-patch energy [J] on the patches' nodes, to be applied over the next step."""
        e = np.asarray(energy, dtype=float)
        if e.any():
            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))

    def _film_nodal(self):
        """The film's mass per node [kg]: each patch's film shared over its three nodes. The film's thermal state
        lives on the nodes because that is where the solver's capacity and the solver's enthalpy live. The deep liquid
        beneath the film rides the same nodes in the same way (amendment of 2026-10-02)."""
        if not self.m_f.size:
            return np.zeros(self.mesh.n_nodes)
        return np.bincount(self.surface.faces.ravel(), np.repeat((self.m_f + self.m_d) / 3.0, 3), self.mesh.n_nodes)

    def _deep_runoff(self, dt, flow, delta_m, reach, T, h_node, h_e, h_p):
        """(ii-a) The deep runoff (amendment of 2026-10-02): the contiguous liquid below the conjugate depth runs off
        under the pressure gradient and the deceleration, and is never sprayed. Returns the deep liquid it saw [kg].

        Plan fact 28(b) keeps melt in its element unless the gas shear can reach it, and fact 29 asked that the
        pressure-gradient and deceleration-driven runoff integrate over the whole liquid depth rather than the film's.
        Both hold here. The deep liquid of a patch is (a) the fully liquid inventory f_feed phi rho V of the elements of
        its contiguous molten chain (`_molten_chain`, the march of `molten_depth`) that the shear does not reach -- the
        ones the feed gate held back -- shared among the patches whose chains pass through them by patch area, plus
        (b) the mobile deep liquid `m_d` the runoff has already delivered to it. With h_D that liquid's thickness and b
        the film's, it moves by `film.deep_flux` G ((b + h_D)^3 - b^3) / (3 mu_l) on the same linearly implicit upwind
        transport as the film, on windward patches under Girin's closure only: elsewhere there is no conjugate depth,
        so there is no layer below the shear's reach. A patch that loses deep liquid gives it up in proportion from
        its chain's elements and its `m_d`; what arrives goes to `m_d`. In situ the liquid stays where it is: an element
        is debited only by what actually flows out of it, so the melt pool keeps its place in the mesh and is still
        delivered to the film by the owner rule as the surface recedes to it.
        Only the skin sprays. `m_d` is never offered to the spray: it becomes film only from the top, as far as the
        film is thinner than the conjugate depth (the free surface has receded into it), and all at once where there
        is no conjugate depth or on a leeward patch. Every transfer is booked by the rule of facts 5 and 25: liquid
        leaves an element at the enthalpy its molten part carries, the element keeps the difference to its mean, the
        liquid arrives at its destination patch's liquid enthalpy and the difference is released there; liquid moving
        between `m_d` and the film of one patch moves at one temperature and books nothing. `deep_runoff_mass` counts
        what left the elements, `deep_surfaced_mass` what became film from the top."""
        liq, mat, mesh = self.liquid, self.material, self.mesh
        n, areas = self.surface.n_patches, self.surface.areas
        if not n:
            return 0.0
        girin = self.windward & np.isfinite(delta_m)          # a conjugate depth exists: liquid can lie below it
        rows = e = share = None
        D = self.m_d.copy()
        if reach is not None:                                  # feed_depth "all" feeds every element: nothing is held back
            _, chain = self._molten_chain()
            rows, level = np.nonzero(chain >= 0)
            e = chain[rows, level]
            held = girin[rows] & (reach[e] <= 0.0)             # in a chain under a Girin patch, and the feed gate held it
            rows, e = rows[held], e[held]
            f = mat.feed_fraction(T)[mesh.tets].mean(axis=1) * mesh.active          # its fully liquid part (facts 5, 28)
            liquid = np.minimum(f * self.phi, np.maximum(self.phi - PHI_MIN, 0.0)) * self.element_mass   # never below PHI_MIN
            weight = np.bincount(e, areas[rows], mesh.n_elements)
            share = liquid[e] * areas[rows] / np.where(weight[e] > 0.0, weight[e], 1.0)
            D = D + np.bincount(rows, share, n)
        seen = float(D.sum())
        if self.runoff is not None and seen > 0.0:
            b = self.m_f / (liq.rho * areas)
            G = np.where(girin, flow.G, 0.0)
            D_new, _, _ = self.runoff.transport(D, lambda hd: self._film_mod.deep_flux(G, b, b + hd, liq.mu),
                                                    self.t_hat, liq.rho, areas, dt)
            change = D_new - D
            gain, loss = np.maximum(change, 0.0), np.maximum(-change, 0.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(D > 0.0, np.minimum(loss / np.where(D > 0.0, D, 1.0), 1.0), 0.0)
            from_m_d = r * self.m_d
            carried, left = float((from_m_d * h_p).sum()), float(from_m_d.sum())
            if e is not None and e.size:
                out = np.bincount(e, r[rows] * share, mesh.n_elements)          # kg leaving each held element
                fn = mat.feed_fraction(T)[mesh.tets]
                fl = fn.mean(axis=1)
                with np.errstate(divide="ignore", invalid="ignore"):            # per kg: what its molten part carries
                    h_out = np.where(fl > 0.0, (fn * h_node[mesh.tets]).mean(axis=1) / np.where(fl > 0.0, fl, 1.0), h_e)
                self.phi = self.phi - out / self.element_mass
                self._defer_to_elements(out * h_e - out * h_out)               # the element keeps what the liquid left
                carried += float((out * h_out).sum())
                left += float(out.sum())
                self.deep_runoff_mass += float(out.sum())                     # mobilised: what left the elements
            # what arrives is exactly what left: the direct solve conserves only to its conditioning, 1e-12 to 2e-11 of
            # the mass moved under these stiff coefficients (measured 2026-10-02), so the arrivals are scaled to the
            # departures -- a correction of that size, which makes the books exact by construction
            if gain.sum() > 0.0:
                gain = gain * (left / float(gain.sum()))
            self.m_d = self.m_d - from_m_d + gain
            self._spread_to_patches(carried - float((gain * h_p).sum()), gain)  # ... which arrives at its patch's enthalpy
        # the skin: deep liquid becomes film only from the top, as far as the film is thinner than the conjugate depth
        with np.errstate(invalid="ignore"):
            room = np.where(girin, np.maximum(liq.rho * areas * np.where(girin, delta_m, 0.0) - self.m_f, 0.0), np.inf)
        up = np.minimum(self.m_d, room)
        if up.any():
            self.m_f = self.m_f + up
            self.m_d = self.m_d - up
            self.deep_surfaced_mass += float(up.sum())
        return seen



    def _add_to_film(self, fed, carried):
        """Distribute each element's freed liquid to the film: owners to their own patches by area, interior elements
        to the NEAREST_PATCHES nearest patches by area. `carried` is the enthalpy that liquid takes with it [J per
        element]; it rides the same weights, so the caller knows what arrived on each patch. Returns the per-patch
        increment [kg] and the per-patch enthalpy carried in [J]."""
        before, got = self.m_f.copy(), np.zeros_like(self.m_f)
        k = np.flatnonzero(fed > 0.0)
        if k.size == 0:
            return np.zeros_like(self.m_f), got
        owner = self.owner_area[k] > 0.0
        # owners: shared by patch area; interior elements: the nearest patch
        if owner.any():
            ko = k[owner]
            share, share_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
            share[ko], share_h[ko] = fed[ko] / self.owner_area[ko], carried[ko] / self.owner_area[ko]
            self.m_f += share[self.surface.owner] * self.surface.areas
            got += share_h[self.surface.owner] * self.surface.areas
        if (~owner).any():                                    # interior elements: the NEAREST_PATCHES nearest patches, by area
            ki = k[~owner]
            kk = min(NEAREST_PATCHES, self.surface.n_patches)
            _, near = self._patch_tree.query(self.mesh.points[self.mesh.tets[ki]].mean(axis=1), k=kk)
            near = near.reshape(len(ki), kk)
            w = self.surface.areas[near] / self.surface.areas[near].sum(axis=1, keepdims=True)
            np.add.at(self.m_f, near.ravel(), (fed[ki][:, None] * w).ravel())
            np.add.at(got, near.ravel(), (carried[ki][:, None] * w).ravel())
        return self.m_f - before, got

    def _film_and_spray(self, t, dt, state, h_p, flow=None, delta_m=None):
        liq, s, mat = self.liquid, self.settings, self.material
        from . import spray as spray_mod
        if state is None or self.m_f.sum() <= 0.0 or flow is None:
            self.last_flow = self.last_spray = self.last_nonrigid = None
            return 0.0
        areas = self.surface.areas
        # the depth of liquid under each patch: its film plus the contiguous molten material beneath it. The branch
        # tests below ask whether the gas shear reaches the bottom of the liquid, so they belong on this; the fluxes
        # and the release stay on the film, which is the mass that can actually move (decided 2026-09-22).
        molten = self.molten_depth()
        deep = self.m_d / (liq.rho * areas)          # the deep liquid the runoff delivered: part of the layer, never film
        # The thin branch needs a rigid substrate (amendment of 2026-10-06): the regime test -- lubrication's branch and
        # the spray's -- reads what lies under the film down to rigid material, the slurry above 50 % liquid included,
        # instead of the fully molten material alone. The Rayleigh-Taylor criteria and the molten region the wave-fits
        # test measures against keep the liquid.
        nonrigid = self.nonrigid_depth() if s.rigid_substrate else None
        under = molten if nonrigid is None else nonrigid
        # (ii) lubrication and runoff
        n_sub, moved = 0, 0.0
        if self.runoff is not None:
            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + under + deep)[1]
            before = self.m_f
            self.m_f, n_sub, moved = self.runoff.transport(self.m_f, q_of_b, self.t_hat, liq.rho, areas, dt)
            self.runoff_mass += moved                                              # mass that arrived on another patch
            # film that runs to a colder patch takes its enthalpy with it and arrives at that patch's temperature;
            # the difference is released where it lands (the per-edge detail is not carried back, so it is shared
            # over the patches that gained film, which is exact in total and second order in the attribution)
            d = self.m_f - before
            self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
        b = self.m_f / (liq.rho * areas)
        layer = b + molten + deep
        regime_layer = layer if nonrigid is None else b + nonrigid + deep
        v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=regime_layer)
        # the molten surface: the area over which a wave could form at all, either wetted by the film or molten in its
        # own right. Its contiguous extent is what every mode's wavelength is measured against (2026-09-24).
        wetted = (b >= spray_mod.B_MIN) | (molten > 0.0)
        extent = self.region_extent(wetted)
        # Girin & Kopyt's W sin(Theta): the deceleration normal to the film, which on this body is W cos(phi). It is the
        # whole deceleration at the stagnation point and vanishes at the equator, where W lies in the surface.
        w_n = flow.deceleration * np.maximum(np.cos(self.theta), 0.0)
        # (iii) spraying. Off Girin's closure a film whose base is slurry takes no shear mode; what it holds is recorded
        on_slurry = None if nonrigid is None else nonrigid > 0.0
        held = None if nonrigid is None else (self.windward & (b >= spray_mod.B_MIN) & ~np.isfinite(delta_m) & on_slurry)
        held_mass = float("nan") if held is None else float(self.m_f[held].sum())
        res = self.spray.evaluate(flow, state, b, delta_m, v_s, self.windward, dt, areas, self.m_f,
                                  self.transverse_radius, b_layer=layer, extent=extent, deceleration_n=w_n,
                                  regime_layer=None if nonrigid is None else regime_layer, on_slurry=on_slurry)
        released = float(res.dm.sum())
        if released > 0.0:
            if self.spray_onset is None:
                self.spray_onset = t
            if s.particles:
                self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state)
            hn, hm = spray_mod.histogram(res.r, res.dn, res.dm)
            self.hist_n += hn
            self.hist_m += hm
            self.removed_enthalpy += float((res.dm * h_p).sum())   # the droplets keep the surface's superheat, and go
            self.m_f = self.m_f - res.dm
            self.sprayed_mass += released
            self.removed_mass += released
            self.n_released += float(res.dn.sum())
        self.m_f[self.m_f < 1e-30] = 0.0                       # no denormal films (they made 0/0 coefficients in the runoff)
        self.last_flow, self.last_spray, self.last_b, self.last_nonrigid = flow, res, b, nonrigid
        self.last_face_ids = self.surface.face_ids.copy()          # deaths later in this step rebuild the surface
        # what the rigid-substrate test did this step, on the wet windward patches (thick_branch_fraction's set); the two
        # fractions are over those under Girin's closure, where the depth decides: made thick only by the slurry, and made
        # thin by the line's depth where the whole-element liquid layer said thick (fact 58's count, removed)
        wet_w = self.windward & (self.m_f > 0.0)
        girin_wet = wet_w & np.isfinite(delta_m)
        if nonrigid is None or not girin_wet.any():
            slurry_thick = rigid_thin = float("nan")
        else:
            with np.errstate(invalid="ignore"):
                liquid_thick = layer > delta_m
            slurry_thick = float((thick & ~liquid_thick)[girin_wet].mean())
            rigid_thin = float((liquid_thick & ~thick)[girin_wet].mean())
        nonrigid_mean = float(nonrigid[wet_w].mean() * 1e3) if nonrigid is not None and wet_w.any() else float("nan")
        # The Rayleigh-Taylor criterion is reported, never applied, so report it usefully: a body-level "any patch"
        # boolean says nothing about how much melt is involved or whether the unstable wave even fits on the nose.
        # `rt_mass_fraction` is the share of the film on unstable patches; `rt_wavelength_over_nose` is the shortest
        # unstable wavelength against the nose diameter -- above 1 the mode cannot develop on the cap at all.
        rt_on, rt_lam, rt_tau = spray_mod.rayleigh_taylor(w_n, layer, liq)
        # the depth criterion above only says the pool is deep enough for the fastest mode to see it as deep
        # (W h^2 rho > 3 Sigma is h > lambda*/2pi); whether a wave fits across the pool is a second, lateral question,
        # answered against the molten region each patch belongs to rather than the body's largest one
        rt_bound, rt_bl, rt_bt = spray_mod.rayleigh_taylor_bounded(w_n, layer, extent, liq)
        rt_extent = float(extent[np.asarray(rt_on) & (self.m_f > 0.0)].max()) if (np.asarray(rt_on) & (self.m_f > 0.0)).any() else 0.0
        wet = self.m_f > 0.0
        rt_sel = np.asarray(rt_on) & wet
        rt_mass = float(self.m_f[rt_sel].sum() / self.m_f.sum()) if self.m_f.sum() > 0.0 else 0.0
        with np.errstate(invalid="ignore"):
            rt_fit = float(np.nanmin(rt_lam[rt_sel]) / (2.0 * self.nose_radius())) if rt_sel.any() else float("nan")
        r = res.r[res.dm > 0.0]
        self.last_melt.update({
            "runoff_substeps": n_sub, "n_released": float(res.dn.sum()), "regime_fractions": flow.closure_fractions(areas, self.windward),
            "kn_body": flow.kn_body, "re_shock": flow.re_shock, "flow_branch": float(flow.branch),
            "kn_local_stag": float(flow.kn_local[self.i_stag]), "p_w_stag": float(flow.p_w[self.i_stag]),
            "phi_sonic_deg": math.degrees(flow.phi_sonic) if np.isfinite(flow.phi_sonic) else float("nan"),
            "theta_cr_deg": float(np.degrees(self.theta[res.unstable].min())) if res.unstable.any() else float("nan"),
            "spraying_area_m2": float(areas[res.unstable].sum()), "rt_active": float(res.rt_active),
            "r_median_um": float(np.median(r) * 1e6) if r.size else float("nan"), "r_max_um": float(r.max() * 1e6) if r.size else float("nan"),
            "film_thickness_max_mm": float(b.max() * 1e3),
            "rt_mass_fraction": rt_mass, "rt_wavelength_over_nose": rt_fit,
            "rt_growth_ms": float(np.nanmedian(np.asarray(rt_tau)[rt_sel]) * 1e3) if rt_sel.any() else float("nan"),
            "rt_region_mm": rt_extent * 1e3,
            "rt_bounded_fraction": float(self.m_f[np.asarray(rt_bound) & (self.m_f > 0.0)].sum() / self.m_f.sum())
            if self.m_f.sum() > 0.0 else 0.0,
            "rt_bounded_growth_ms": float(np.nanmedian(np.asarray(rt_bt)[np.asarray(rt_bound) & (self.m_f > 0.0)]) * 1e3)
            if (np.asarray(rt_bound) & (self.m_f > 0.0)).any() else float("nan"),
            "spray_growth_ms": float(np.nanmedian(res.growth[res.dm > 0.0]) * 1e3)
            if res.growth is not None and (res.dm > 0.0).any() and np.isfinite(res.growth[res.dm > 0.0]).any() else float("nan"),
            "molten_depth_max_mm": float(molten.max() * 1e3),
            "molten_depth_mean_mm": float((self.m_f * molten).sum() / self.m_f.sum() * 1e3) if self.m_f.sum() > 0.0 else 0.0,
            "delta_m_mean_um": self._mean_conjugate_um(delta_m),
            "thick_branch_fraction": float(thick[self.windward & (self.m_f > 0.0)].mean())
            if (self.windward & (self.m_f > 0.0)).any() else float("nan"),
            "nonrigid_depth_mean_mm": nonrigid_mean, "slurry_thick_fraction": slurry_thick,
            "rigid_thin_fraction": rigid_thin, "slurry_held_mass_kg": held_mass})
        return released

    def _kill(self, dead):
        old_surface, old_m_f, old_m_d = self.surface, self.m_f, self.m_d
        gone, new = self.mesh.deactivate(dead)
        if self.mesh.n_active == 0:                              # nothing left: the run ends (demise) without a surface
            self.consumed = True                                 # the liquid still on it leaves with the body
            self.removed_mass += float(self.m_f.sum() + self.m_d.sum())
            self.removed_enthalpy += float(((self.m_f + self.m_d) * self.film_enthalpy()).sum())
            self.m_f, self.m_d = np.zeros(0), np.zeros(0)
            return
        self._refresh_geometry()
        keep = self.patch_of_face[old_surface.face_ids]
        kept = keep >= 0
        self.m_f = self._hand_over(old_surface, old_m_f, keep, kept)
        self.m_d = self._hand_over(old_surface, old_m_d, keep, kept)   # the deep liquid beneath it goes the same way

    def _hand_over(self, old_surface, old, keep, kept):
        """A per-patch liquid mass of the surface before the deaths, carried onto the current surface: kept patches keep
        theirs; a vanished patch's goes to the faces its own element exposed (below)."""
        m_f = np.zeros(self.surface.n_patches)
        m_f[keep[kept]] += old[kept]
        lost = ~kept & (old > 0.0)
        if lost.any():
            # The film on a vanished patch has not moved: the surface receded *into* it, so it now rests on the faces
            # of that same element which have just become boundary -- the face on its inward side, whose neighbour is
            # the next layer in, and the walls of the pit its lateral neighbours expose -- shared by their areas.
            # Handing it to the nearest surviving centroid instead sends it sideways to whatever triangle happens to
            # be closest, which put a third of a collapsing body's film onto one 0.6 mm2 face (measured 2026-09-22),
            # and is a rule about the mesh's numbering rather than about where the liquid is.
            idx = np.flatnonzero(lost)
            targets = self.patch_of_face[self.mesh.element_faces(old_surface.owner[idx])]   # (n, 4), -1: not a patch
            ok = targets >= 0
            safe = np.where(ok, targets, 0)
            # weight by each exposed face's area *projected on the vanished patch's own outward normal*, not by its raw
            # area: the face the recession uncovered underneath lies parallel to the patch that vanished (its new
            # outward normal points the same way), while the pit walls stand perpendicular to it. A prism element's
            # side faces can out-area its inward face -- 0.25 x 1.5 mm slivers against a 1 mm2 floor -- so raw areas
            # route most of the film sideways into walls no shear is pushing it towards, and concentrate it by the
            # element's aspect ratio (3 to 8 here). The projection is the footprint each face offers the liquid.
            cosine = np.clip(np.einsum("nij,nj->ni", self.surface.normals[safe], old_surface.normals[idx]), 0.0, None)
            area = np.where(ok, self.surface.areas[safe] * cosine, 0.0)
            total = area.sum(axis=1)
            flat = ok & (total[:, None] <= 0.0)                  # no face faces the right way: fall back to raw areas
            if flat.any():
                area = np.where(flat, self.surface.areas[safe], area)
                total = area.sum(axis=1)
            has = total > 0.0
            if has.any():
                share = old[idx[has]][:, None] * area[has] / total[has][:, None]
                np.add.at(m_f, targets[has][ok[has]], share[ok[has]])
            if (~has).any():                                     # the element exposed nothing: its death opened a
                orphan = idx[~has]                               # hole right through, so fall back to the nearest patch
                np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old[orphan])
        return m_f

    # -- reporting -----------------------------------------------------------------------------------------------
    def on_current_surface(self, values, fill):
        """A per-patch array from the last flow and spray evaluation (`last_flow`, `last_spray`), carried onto the
        current surface. Element deaths later in the same melt step rebuild the surface, so the arrays no longer line
        up with it by index; patches are matched by face id instead, and a patch that was not on the surface at the
        evaluation (a face the deaths exposed) gets `fill`. Matching by index instead left every frame written after
        a step with a death -- every melting frame -- at the defaults: 0 of 104 frames of the 50 mm physics flight
        carried a wall pressure (measured 2026-09-27)."""
        n = self.surface.n_patches
        if values is None or self.last_face_ids is None or not self.last_face_ids.size:
            return np.full(n, fill, dtype=float)
        index = np.full(len(self.patch_of_face), -1, dtype=np.int64)
        index[self.last_face_ids] = np.arange(self.last_face_ids.size)
        j = index[self.surface.face_ids]
        return np.where(j >= 0, np.asarray(values, dtype=float)[np.maximum(j, 0)], fill)

    def film_temperature(self):
        """The film's temperature per patch: the surface's own (the film is thermally thin, class docstring)."""
        return self.surface_temperature()

    def film_enthalpy(self, h_node=None):
        """Specific enthalpy of each patch's film [J/kg]: the mean over the patch's three nodes of the *liquid*
        enthalpy h(T_i) + L_f (1 - f_l) -- the film is liquid by construction, so it carries its latent heat wherever
        it sits. Two things this must not be: the enthalpy of the mean temperature (the mean of the enthalpies and
        the enthalpy of the mean differ by the whole latent heat across the ramp, which left a -0.4 residual in the
        coupled balance), and the mixture enthalpy h(T) (which makes melting free on the ramp)."""
        if h_node is None:
            h_node = self.material.enthalpy_liquid(self.solver.temperature())
        return h_node[self.surface.faces].mean(axis=1) if self.m_f.size else np.zeros(0)

    def film_energy(self):
        """Enthalpy held by the film [J]. Identical to the nodal sum the solver carries (`set_film_mass` gives each
        node sum_p m_p/3), so the film's sensible heat is inside 1^T M dT and the balance closes on it exactly. The
        deep liquid beneath the film holds the same liquid enthalpy on the same nodes (amendment of 2026-10-02)."""
        return float(((self.m_f + self.m_d) * self.film_enthalpy()).sum()) if self.m_f.size else 0.0

    def film_frozen_fraction(self):
        """Fraction of the film sitting on patches below the feed ramp: mass the enthalpy calls solid but the model
        still treats as liquid, because the mesh cannot grow a crust outside itself and its own element may be gone
        (section 9.5). It is the honest measure of that limitation, so it is recorded every step."""
        total = float(self.m_f.sum())
        if total <= 0.0:
            return 0.0
        return float(self.m_f[self.material.feed_fraction(self.film_temperature()) <= 0.0].sum() / total)

    @staticmethod
    def _mean_conjugate_um(delta_m):
        """Mean conjugate melt-layer depth [um] over the patches that have one; nan where the gate denies them all."""
        ok = np.isfinite(delta_m)
        return float(delta_m[ok].mean() * 1e6) if ok.any() else float("nan")

    def _shear_reaches(self, delta_m):
        """Per element, may its melt leave for the film this step? An element that owns part of the wall is in the
        sheared layer by definition and always may. A buried element may only if its centroid lies within Girin's
        conjugate melt-layer depth of the nearest patch, because that is how far the gas shear penetrates the liquid;
        deeper melt is stagnant and stays where it is until the surface reaches it. Where the wall Knudsen gate denies
        the Girin closure there is no conjugate depth, and then only wall-owning elements may feed, which is the
        conservative reading of the same statement."""
        owner = self.owner_area > 0.0
        if delta_m is None or not self.surface.n_patches:
            return owner.astype(float)
        dist, near = self._patch_tree.query(self.mesh.points[self.mesh.tets].mean(axis=1))
        d = delta_m[near]
        return (owner | (np.isfinite(d) & (dist <= d))).astype(float)

    def region_extent(self, selected):
        """Lateral extent [m] of *each* patch's own contiguous selected region, as that region's equivalent diameter
        2 sqrt(A/pi); zero where the patch is not selected.

        This is the domain a surface wave has to fit inside, and the region -- not the mesh facet -- is the physical
        boundary of it: melt is continuous across a facet edge, so a facet width is a bookkeeping artefact and using it
        would make the limit mesh-dependent (decided 2026-09-24). Contiguity is taken on the patch adjacency graph, so
        melt separated by solid does not add up to one wide pool. Each patch gets its own region rather than the
        body's largest, since a patch in an isolated puddle must not be credited with the width of the main pool."""
        n = self.surface.n_patches
        sel = np.asarray(selected, dtype=bool)
        out = np.zeros(n)
        if not n or not sel.any():
            return out
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        _, i, j = self.surface.edges()
        keep = sel[i] & sel[j]
        i, j = i[keep], j[keep]
        g = coo_matrix((np.ones(i.size * 2), (np.concatenate([i, j]), np.concatenate([j, i]))), shape=(n, n))
        count, label = connected_components(g, directed=False)
        area = np.bincount(label[sel], self.surface.areas[sel], count)
        out[sel] = 2.0 * np.sqrt(area[label[sel]] / np.pi)
        return out

    def unstable_region_extent(self, unstable):
        """Lateral extent [m] of the largest contiguous flagged region: the maximum of `region_extent`, kept for the
        reported diagnostics. A Rayleigh-Taylor wave has to fit inside the region where the melt is deep enough to
        support it, not on the nose as a whole: comparing its wavelength against the nose diameter asks whether the
        wave would fit on a surface most of which carries no pool (corrected 2026-09-23)."""
        e = self.region_extent(unstable)
        return float(e.max()) if e.size else 0.0

    def molten_depth(self, Te=None, max_levels=8):
        """Depth of *contiguous* fully molten material under each patch [m], measured inward from the patch.

        Marches inward through the mesh from the patch's owner element, hopping to the face-adjacent element furthest
        along the inward normal, and stops at the first element that is not fully molten (below `Material.T_feed`) or
        at the mesh boundary. Only contiguous liquid counts: molten material separated from the patch by solid is
        blocked from it, cannot be part of the layer the gas shear sees, and must not be credited to it -- accumulating
        every nearby molten element instead would count melt that cannot reach the wall (decided 2026-09-22).

        This is the thickness that belongs in the thick/thin instability test against Girin's conjugate melt-layer
        depth delta_m: the film mass on a patch is the mobile inventory, not the depth of liquid beneath the wall.
        Measured on the 100 mm physics flight, the two differ by an order of magnitude and invert the branch choice."""
        return self._molten_chain(Te, max_levels)[0]

    def _molten_chain(self, Te=None, max_levels=8):
        """`molten_depth`'s march, also returning the elements it passed through: (depth per patch [m], chain), the
        chain an (n_patches, max_levels) array of element ids, level 0 the patch's owner and -1 past the chain's end.
        The deep runoff needs the elements (amendment of 2026-10-02); the branch test needs only the depth."""
        mat, mesh = self.material, self.mesh
        chain = np.full((self.surface.n_patches, max_levels), -1, dtype=np.int64)
        if not self.surface.n_patches or not mat.melts:
            return np.zeros(self.surface.n_patches), chain
        if Te is None:
            Te = self.solver.temperature()[mesh.tets].mean(axis=1)
        molten = mesh.active & (Te >= mat.T_feed)
        centroid = mesh.points[mesh.tets].mean(axis=1)
        inward = -self.surface.normals
        cur = self.surface.owner.copy()
        alive = molten[cur]
        depth = np.zeros(self.surface.n_patches)
        rows = np.arange(len(cur))
        for level in range(max_levels):
            if not alive.any():
                break
            chain[alive, level] = cur[alive]
            reach = np.einsum("ij,ij->i", centroid[cur] - self.surface.centroids, inward)      # projected depth so far
            depth = np.where(alive, np.maximum(depth, reach), depth)
            pair = mesh._face_elements[mesh.element_faces(cur)]                               # (n, 4, 2)
            nb = np.where(pair[:, :, 0] == cur[:, None], pair[:, :, 1], pair[:, :, 0])
            ok = nb >= 0
            step = np.einsum("nkj,nj->nk", centroid[np.where(ok, nb, 0)] - centroid[cur][:, None, :], inward)
            pick = np.where(ok, step, -np.inf).argmax(axis=1)
            nxt = nb[rows, pick]
            good = (nxt >= 0) & (step[rows, pick] > 0.0)
            alive = alive & good & molten[np.where(good, nxt, 0)]
            cur = np.where(good, nxt, cur)
        # the chain's last centroid sits half an element short of the far wall of that element: add that half
        return np.maximum(depth, 0.0) * 1.5, chain

    def liquid_layer_depth(self, Te=None):
        """Depth of liquid at each patch [m]: its film plus the contiguous molten material beneath it, plus the deep
        liquid the runoff has delivered there (amendment of 2026-10-02). What the gas shear sees, and so what decides
        whether the film is thick or thin against the conjugate melt-layer depth."""
        if not self.m_f.size:
            return np.zeros(0)
        a = self.liquid.rho * self.surface.areas
        return self.m_f / a + self.molten_depth(Te) + self.m_d / a

    def nonrigid_depth(self, T=None, max_crossings=NONRIGID_MAX_CROSSINGS):
        """Depth of non-rigid material under each patch [m]: how far inward from the wall the material stays more than
        half liquid (above `Material.T_rigid`), each element's stretch weighed by phi_e (amendment of 2026-10-06).

        A line runs from the patch centre along the inward normal through the elements it actually crosses, leaving each
        by the face its barycentric coordinates reach first. The P1 field is linear along the line inside an element, so
        the point where it falls to T_rigid is found exactly instead of being counted in whole elements, and plan fact
        58's element-size artefact does not arise. Only material continuous with the wall counts: the line stops at the
        first point at or below T_rigid, at the active mesh's boundary, or after `max_crossings` elements. Each
        element's stretch counts phi_e of its length, because what it has already fed to the film is in the film
        account the layer adds on top (fact 58's double count). The thin branch is Girin's dominant ablation, where the
        rigid core stabilises the disturbances, and slurry is no rigid core: this is the depth the regime test reads."""
        mat, mesh, surf = self.material, self.mesh, self.surface
        n = surf.n_patches
        depth = np.zeros(n)
        if not n or not mat.melts:
            return depth
        T = self.solver.temperature() if T is None else np.asarray(T, dtype=float)
        T_r = mat.T_rigid
        c, d = surf.centroids, -surf.normals
        elem, s_in, live = surf.owner.copy(), np.zeros(n), np.arange(n)
        for _ in range(max_crossings):
            if live.size == 0:
                break
            e, s0 = elem[live], s_in[live]
            x = mesh.points[mesh.tets[e]]                                                       # (k, 4, 3)
            inv = np.linalg.inv(np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2))
            grad = np.concatenate([-inv.sum(axis=1, keepdims=True), inv], axis=1)               # grad of lambda_0..3
            lam = np.einsum("kij,kj->ki", inv, c[live] - x[:, 0])
            lam = np.concatenate([1.0 - lam.sum(axis=1, keepdims=True), lam], axis=1)           # barycentrics of the centre
            rate = np.einsum("kij,kj->ki", grad, d[live])                                       # d lambda / ds on the line
            s_face = np.where(rate < 0.0, -lam / np.where(rate < 0.0, rate, -1.0), np.inf)       # where each face is reached
            exit_face = s_face.argmin(axis=1)
            rows = np.arange(len(e))
            s1 = np.maximum(s_face[rows, exit_face], s0)
            Tn = T[mesh.tets[e]]
            T0, dT = (lam * Tn).sum(axis=1), (rate * Tn).sum(axis=1)                            # T on the line: T0 + s dT
            Ta, Tb = T0 + s0 * dT, T0 + s1 * dT
            hot = Ta > T_r
            ends = hot & (Tb <= T_r)
            with np.errstate(divide="ignore", invalid="ignore"):
                s_end = np.where(ends, s0 + (Ta - T_r) / np.where(ends, Ta - Tb, 1.0) * (s1 - s0), s1)
            depth[live] += np.where(hot, self.phi[e] * (s_end - s0), 0.0)
            pair = mesh._face_elements[mesh.element_faces(e)[rows, exit_face]]                  # across the exit face
            nxt = np.where(pair[:, 0] == e, pair[:, 1], pair[:, 0])
            go = hot & ~ends & np.isfinite(s1) & (nxt >= 0)
            go &= mesh.active[np.maximum(nxt, 0)]
            elem[live[go]], s_in[live[go]] = nxt[go], s1[go]
            live = live[go]
        return depth

    def film_thickness_max(self):
        """Thickest film [m] over the patches where a thickness means anything: at least a tenth of the median patch
        area (slivers excluded) and no deeper than the patch is wide (b <= sqrt(A)). A film deeper than its patch is
        wide is not a film -- m_f/(rho_l A) is then a volume per area, not a depth, and the lubrication picture the
        number belongs to has already failed. Those patches are reported separately by `film_blob_fraction` rather
        than quietly averaged in: without the second test this diagnostic read 2431 mm on a 50 mm sphere in the last
        step of its flight, where thousands of dying patches handed their film to a few 0.6 mm2 survivors (measured
        2026-09-22; the mass was real, the depth was not)."""
        if not self.m_f.size:
            return 0.0
        a = self.surface.areas
        b = self.m_f / (self.liquid.rho * a)
        ok = (a >= 0.1 * np.median(a)) & (b <= np.sqrt(a))
        return float(b[ok].max()) if ok.any() else 0.0

    def deep_blob_fraction(self):
        """Fraction of the deep liquid (`m_d`) sitting deeper than its patch is wide (b > sqrt(A)), the deep account's
        counterpart of `film_blob_fraction`. The deep runoff carries the liquid below the conjugate depth to where the
        flow converges and nothing strips it there, so it piles up on the patches the patch graph cannot drain -- on
        the 100 mm flight almost all of it (amendment of 2026-10-02). The mass is exact; the depth is not a film's."""
        total = float(self.m_d.sum())
        if total <= 0.0:
            return 0.0
        a = self.surface.areas
        return float(self.m_d[self.m_d / (self.liquid.rho * a) > np.sqrt(a)].sum() / total)

    def film_blob_fraction(self):
        """Fraction of the film sitting deeper than its patch is wide (b > sqrt(A)) -- melt the lubrication film
        model cannot describe, because it no longer fits on the facet it is booked to. It is a geometry limit of the
        patch graph, not a mass error: the mass and the energy are still accounted exactly. Measured on the 50 mm
        physics flight it is zero until the last 5 % of the flight, where the body collapses from 3026 patches to 56."""
        total = float(self.m_f.sum())
        if total <= 0.0:
            return 0.0
        a = self.surface.areas
        return float(self.m_f[self.m_f / (self.liquid.rho * a) > np.sqrt(a)].sum() / total)

    def film_thickness_mean(self):
        """Film mass over the area it meaningfully wets [m]: patches at least B_MIN deep, the same floor below which
        spraying ignores a film. Counting every patch holding a nonzero number instead makes this diagnostic
        hypersensitive -- a patch with 1e-20 kg on it adds its whole area to the denominator, so a last-bit difference
        in the temperature field that nudges one element across the feed ramp moves the reported mean by 1-2 % while
        the film mass is bit-identical. That is how a reporting artefact was mistaken for a physics change
        (2026-09-23); the mass-based diagnostics (`film_blob_fraction`, `film_frozen_fraction`) never had the problem,
        because a patch holding nothing contributes nothing to a mass fraction."""
        from .spray import B_MIN
        if not self.m_f.size:                                    # no film account at all (e.g. --removal instant)
            return 0.0
        wet = self.m_f > self.liquid.rho * self.surface.areas * B_MIN
        return float(self.m_f[wet].sum() / (self.liquid.rho * self.surface.areas[wet].sum())) if wet.any() else 0.0

    def demised(self):
        """True once the body's material (film excluded) is below the demise fraction of the initial mass."""
        return self.consumed or float((self.phi * self.element_mass).sum()) < self.settings.demise_fraction * self.mass0

    def energy_balance_residual(self):
        """(E_body + E_film - E0 - absorbed heat + removed enthalpy + pending load + dropped load) / absorbed heat:
        zero for an exact discrete balance (the deferred loads cancel between the FEM and the removed material)."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy + self.pending_load.sum()
                + self.total_dropped_load) / self.absorbed_heat

    @property
    def total_applied_load(self):
        return getattr(self, "_applied_total", 0.0)

    @property
    def total_dropped_load(self):
        return getattr(self, "_dropped_total", 0.0)

    def melt_stats(self):
        lm = self.last_melt
        fractions = lm.get("regime_fractions", [0.0, 0.0, 0.0])
        return {"mass_kg": self.mass(0.0), "film_mass_kg": float(self.m_f.sum()), "sprayed_mass_kg": self.sprayed_mass,
                "runoff_mass_kg": self.runoff_mass, "removed_mass_kg": self.removed_mass,
                "melt_front_depth_max_mm": self.melt_front_depth() * 1e3, "equivalent_radius_mm": self.equivalent_radius() * 1e3,
                "n_active_elements": float(self.mesh.n_active), "spraying_area_m2": lm.get("spraying_area_m2", 0.0),
                "theta_cr_deg": lm.get("theta_cr_deg", float("nan")), "n_released": self.n_released,
                "released_mass_kg": lm.get("released_mass", 0.0), "r_median_um": lm.get("r_median_um", float("nan")),
                "r_max_um": lm.get("r_max_um", float("nan")), "closure_fraction_girin": fractions[0],
                "closure_fraction_couette_slip": fractions[1], "closure_fraction_couette_fm": fractions[2], "rt_active": lm.get("rt_active", 0.0),
                "kn_body": lm.get("kn_body", float("nan")), "kn_local_stag": lm.get("kn_local_stag", float("nan")),
                "re_shock": lm.get("re_shock", float("nan")), "flow_branch": lm.get("flow_branch", float("nan")),
                "p_w_stag_Pa": lm.get("p_w_stag", float("nan")), "phi_sonic_deg": lm.get("phi_sonic_deg", float("nan")),
                "drag_shape_factor": self.drag_shape_factor(), "frozen_mass_kg": self.frozen_mass,
                "film_T_max_K": float(self.film_temperature().max()) if self.m_f.size else float("nan"),
                "film_T_mean_K": float((self.m_f * self.film_temperature()).sum() / self.m_f.sum()) if self.m_f.sum() > 0.0 else float("nan"),
                "film_frozen_fraction": self.film_frozen_fraction(), "unapplied_load_J": float(self.pending_load.sum()),
                "film_blob_fraction": self.film_blob_fraction(),
                "rt_mass_fraction": lm.get("rt_mass_fraction", float("nan")),
                "rt_growth_ms": lm.get("rt_growth_ms", float("nan")),
                "rt_region_mm": lm.get("rt_region_mm", float("nan")),
                "rt_bounded_fraction": lm.get("rt_bounded_fraction", float("nan")),
                "rt_bounded_growth_ms": lm.get("rt_bounded_growth_ms", float("nan")),
                "spray_growth_ms": lm.get("spray_growth_ms", float("nan")),
                "rt_wavelength_over_nose": lm.get("rt_wavelength_over_nose", float("nan")),
                "molten_depth_max_mm": lm.get("molten_depth_max_mm", float("nan")),
                "molten_depth_mean_mm": lm.get("molten_depth_mean_mm", float("nan")),
                "delta_m_mean_um": lm.get("delta_m_mean_um", float("nan")),
                "thick_branch_fraction": lm.get("thick_branch_fraction", float("nan")),
                "nonrigid_depth_mean_mm": lm.get("nonrigid_depth_mean_mm", float("nan")),
                "slurry_thick_fraction": lm.get("slurry_thick_fraction", float("nan")),
                "rigid_thin_fraction": lm.get("rigid_thin_fraction", float("nan")),
                "slurry_held_mass_kg": lm.get("slurry_held_mass_kg", float("nan")),
                "deep_liquid_kg": lm.get("deep_liquid", 0.0), "deep_mass_kg": float(self.m_d.sum()),
                "deep_runoff_mass_kg": self.deep_runoff_mass, "deep_surfaced_mass_kg": self.deep_surfaced_mass,
                "deep_blob_fraction": self.deep_blob_fraction(),
                "cascade_passes": float(lm.get("cascade_passes", 0)), "cascade_mass_kg": self.cascade_mass,
                "removed_enthalpy_J": self.removed_enthalpy,
                "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
                "fitted_nose_radius_mm": self.fitted_nose_radius * 1e3,
                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active)}
