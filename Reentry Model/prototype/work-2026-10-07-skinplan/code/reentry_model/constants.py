"""Physical, geodetic and gas constants (SI)."""
# WGS-84 ellipsoid
WGS84_A = 6378137.0                       # semi-major axis [m]
WGS84_F = 1.0 / 298.257223563             # flattening
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)      # first eccentricity squared
# gravity field
MU_EARTH = 3.986004418e14                 # [m^3/s^2]
J2 = 1.08262668e-3
J4 = -1.61962159e-6
OMEGA_EARTH = 7.2921159e-5                # Earth rotation rate [rad/s]
G0 = 9.80665                              # standard gravity [m/s^2], for load factors
# gas
K_BOLTZMANN = 1.380649e-23                # [J/K]
ATOMIC_MASS_UNIT = 1.66053906660e-27      # [kg]
M_BAR_AIR = 28.9644 * ATOMIC_MASS_UNIT    # mean molecular mass of air below ~100 km [kg]
GAMMA_AIR = 1.4
HARD_SPHERE_DIAMETER = 3.65e-10           # [m]; reproduces SESAM's Knudsen column (facts note s.2)
