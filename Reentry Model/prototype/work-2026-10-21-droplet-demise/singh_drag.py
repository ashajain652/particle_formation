"""Sphere drag coefficient of Singh, Kroells, Li, Ching, Ihme, Hogan & Schwartzentruber, "General Drag Coefficient
for Flow over Spherical Particles", AIAA Journal (2021), DOI 10.2514/1.J060648, Appendix C, Eqs. (C1)-(C12).
Transcribed 2026-10-22 from the published PDF (author copy, Stanford) and checked against the arXiv preprint
2012.04813 (Appendix C, Eqs. 50-61); see citations/sphere_drag_rarefied_compressible.md.

    C_d = C_d^c(M, Re) f_{Kn,Wr} / (1 + Br^eta) + C_d,fm Br^eta / (1 + Br^eta)                         (C1)

  * C_d^c, the continuum drag (C8): Abraham's boundary-layer form C0 (1 + delta0 / Re~^1/2)^2 (Stokes' 24/Re as
    Re -> 0), corrected for compressibility through the Stewartson-Howarth function Theta(M) (C9) and, above Mach 1,
    evaluated behind a curved bow shock (normal-shock relations mapped by alpha(M), C12) plus a shock-pressure term
    C1 (1 - alpha U_s/U) tuned so that C_d^c -> 0.9 in the hypersonic limit (C11).
  * f_{Kn,Wr} (C5): the Cunningham-Millikan-Davies slip correction (A1, A2, A3 = 2.514, 0.8, 0.55 with Kn = lambda/D)
    times a high-speed rarefaction correction 1/(1 + alpha_hoc Wr^T), Wr^T = M^(2 omega)/Re (1 + T_p/T_s)^omega.
  * C_d,fm (C6): the Schaaf & Chambre / Patterson free-molecule sphere, full diffuse accommodation.
  * Br (C2): the bridging parameter Wr^T (M^(2 omega - 1) + 1) / M^(2 omega - 1), which tends to Wr^T at high Mach
    and to Kn (times a constant of order 1) at low Mach.

Kn = lambda / D = (M / Re) sqrt(pi gamma / 2) (Eq. 2) and s = M sqrt(gamma / 2) (C7).

One interpretation is ours: the paper defines the post-shock state (T_s, U_s, M_s, alpha) for M >= 1 only. Below
Mach 1 we take the free-stream values (T_s = T, U_s = U, M_s = M, alpha = 1), which is what the paper's subsonic form
(its Eq. 5) is, and which makes every term continuous at M = 1, where the normal-shock relations are the identity.
The paper writes Wr^T with the condition "M > 1" in its Eq. (15) but applies it without condition in Appendix C; we
follow Appendix C, using T_s = T below Mach 1.
"""
import math

from scipy.special import erf

GAMMA = 1.4
DELTA0 = 9.4                  # boundary-layer thickness scaling (Table 1)
C0 = 24.0 / DELTA0 ** 2       # C0 delta0^2 = 24 recovers Stokes (Sec. II.A.1)
ALPHA0 = 0.356                # shock-curvature parameter (Table 1)
A1, A2, A3 = 2.514, 0.8, 0.55 # low-speed rarefaction correction, Kn = lambda/D (Table 1)
ETA = 1.8                     # bridging-function modulator (Table 1)
ALPHA_HOC = 1.27              # high-speed rarefaction correction (Table 1)
OMEGA = 0.74                  # viscosity power-law exponent mu ~ T^omega for air (Sec. II.A.3)
CD_HYPERSONIC = 0.9           # hypersonic continuum limit (Hornung et al. 2019, Sec. II.A.4)


def theta(M, g=GAMMA):
    """Stewartson-Howarth compressibility function (C9)."""
    return (1.0 + 0.5 * (g - 1.0) * M * M) ** (g / (g - 1.0))


def post_shock(M, g=GAMMA):
    """(T_s/T, U_s/U, M_s, alpha): normal-shock relations (Eqs. A2.x) and the curvature mapping (C12) for M >= 1;
    the free stream unchanged below Mach 1."""
    if M < 1.0:
        return 1.0, 1.0, M, 1.0
    M2 = M * M
    Ts = ((g - 1.0) * M2 + 2.0) * (2.0 * g * M2 - (g - 1.0)) / ((g + 1.0) ** 2 * M2)
    Us = ((g - 1.0) * M2 + 2.0) / ((g + 1.0) * M2)
    Ms = math.sqrt(((g - 1.0) * M2 + 2.0) / (2.0 * g * M2 - (g - 1.0)))
    alpha = 1.0 / (ALPHA0 * M + 1.0 - ALPHA0)
    return Ts, Us, Ms, alpha


def cd_continuum(M, Re, g=GAMMA, omega=OMEGA):
    """Continuum drag C_d^c(M, Re) (C8), with C1 from (C11)."""
    Ts, Us, Ms, alpha = post_shock(M, g)
    th_s = theta(Ms, g)
    re_s = Re * ((1.0 / alpha ** 2) / Ts) ** omega * th_s ** ((g + 1.0 - 2.0 * (g - 1.0) * omega) / (2.0 * g))
    viscous = C0 * th_s * (1.0 + DELTA0 / math.sqrt(re_s)) ** 2
    if M < 1.0:
        return viscous
    c1 = (CD_HYPERSONIC - C0 * (1.0 + (g - 1.0) ** 2 / (4.0 * g)) ** (g / (g - 1.0))) \
        / (1.0 - (1.0 / (ALPHA0 * M)) * (g - 1.0) / (g + 1.0))
    return c1 * (1.0 - alpha * Us) + viscous


def cd_free_molecular(s, Tp_over_T):
    """Patterson / Schaaf & Chambre sphere, full diffuse accommodation (C6)."""
    if s < 1e-3:
        return (16.0 / (3.0 * math.sqrt(math.pi)) + 2.0 * math.sqrt(math.pi * Tp_over_T) / 3.0) / s
    return ((1.0 + 2.0 * s * s) * math.exp(-s * s) / (s ** 3 * math.sqrt(math.pi))
            + (4.0 * s ** 4 + 4.0 * s * s - 1.0) * erf(s) / (2.0 * s ** 4)
            + 2.0 * math.sqrt(math.pi * Tp_over_T) / (3.0 * s))


def drag_coefficient(M, Re, Tp_over_T, g=GAMMA, omega=OMEGA):
    """Singh et al. (2021) C_d (C1) from the free-stream Mach and Reynolds numbers (Re on the diameter) and the
    particle-to-gas temperature ratio. Returns (C_d, diagnostics)."""
    Kn = (M / Re) * math.sqrt(math.pi * g / 2.0)                 # lambda / D (Eq. 2)
    s = M * math.sqrt(g / 2.0)                                   # (C7)
    Ts, _, _, _ = post_shock(M, g)
    hot = (1.0 + Tp_over_T / Ts) ** omega                        # (1 + T_p/T_s)^omega, T_s in units of T
    m = M ** (2.0 * omega - 1.0)
    Wr_T = Kn * m * math.sqrt(2.0 / (math.pi * g)) * hot         # (C3)-(C4)
    Br = Kn * math.sqrt(2.0 / (math.pi * g)) * hot * (m + 1.0)   # (C2), written without the 0/0 at M -> 0
    f = 1.0 / (1.0 + Kn * (A1 + A2 * math.exp(-A3 / Kn))) / (1.0 + ALPHA_HOC * Wr_T)   # (C5)
    w = Br ** ETA / (1.0 + Br ** ETA)
    cdc = cd_continuum(M, Re, g, omega)
    cdfm = cd_free_molecular(s, Tp_over_T)
    cd = cdc * f * (1.0 - w) + cdfm * w
    return cd, dict(Kn=Kn, s=s, Br=Br, Wr_T=Wr_T, w_fm=w, cd_continuum=cdc, f_kn=f, cd_fm=cdfm)
