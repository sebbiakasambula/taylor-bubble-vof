"""Reference results for the Taylor bubble. Each function says how confident its source is.

INERTIAL LIMIT, confident. For high Eo and high Nf, surface tension and viscosity are negligible and
Fr tends to a constant. Dumitrescu (1943) obtained 0.351 from potential flow; Davies and Taylor (1950)
obtained 0.328 by an approximate method. Viana et al. (2003, eq 16) find Fr = 0.34 from measurements
for Eo > 40. 0.351 is context only; the validation reference is fr_viana. Sources:
  Dumitrescu 1943, ZAMM 23, 139-149, doi 10.1002/zamm.19430230303
  Davies and Taylor 1950, Proc. R. Soc. A 200, 375-390, doi 10.1098/rspa.1950.0023

VIANA ET AL. (2003), implemented as fr_viana. J. Fluid Mech. 494, 379-398, doi 10.1017/S0022112003006165.
Eq (23), valid for Eo > 6. Coefficients checked against the JFM version on 10 Nov 2026 (docs/rise_velocity_gas_bubble.pdf),
identical to the preprint (docs/UnivCorrelation2.pdf).

LLEWELLIN ET AL. (2012), Proc. R. Soc. A 468, 1041-1064, doi 10.1098/rspa.2011.0476. Available as
fr_llewellin (eq 2.5, Eo > 40) and film_llewellin (eq 4.2, film over pipe radius, Eo > 40).

FILM THICKNESS, first principles, for mesh design only. Volume conservation in the bubble frame plus a
laminar Nusselt falling film. This is an estimate to size the mesh, not a validation reference. Use
film_llewellin when a reference film thickness is needed.
"""
import math

from scipy.optimize import brentq

FR_DUMITRESCU = 0.351
FR_DAVIES_TAYLOR = 0.328


def fr_inertial():
    """The inviscid-limit Froude number, with the spread between the two classical values. Context only."""
    return dict(value=FR_DUMITRESCU, low=FR_DAVIES_TAYLOR, high=FR_DUMITRESCU)


def film_thickness(Fr, Nf, R=0.5):
    """Developed film thickness delta (in units of D) from volume conservation and a Nusselt film.

    Bubble frame: liquid far ahead moves down at U_b over the whole pipe, and through the film annulus
    at U_b + U_f, where U_f is the lab-frame mean film velocity. Conservation gives
        U_b (R - delta)^2 = U_f (2 R delta - delta^2).
    A laminar film falling under gravity with a free surface has U_f = g delta^2 / (3 nu).
    With g = 1 and nu = 1/Nf (drho close to rho_l), and U_b = Fr in these units:
        Fr (R - delta)^2 = (Nf / 3) delta^3 (2 R - delta).
    """
    f = lambda d: Fr * (R - d) ** 2 - (Nf / 3.0) * d ** 3 * (2 * R - d)
    return brentq(f, 1e-6, R - 1e-9)


def film_reynolds(delta, Nf):
    """Re_f = U_f delta / nu = Nf^2 delta^3 / 3. Laminar below about 20, wavy laminar to about 400."""
    return Nf ** 2 * delta ** 3 / 3.0


def fr_viana(Eo, Nf):
    """Froude number from Viana et al. (2003), eq (23), valid for Eo > 6.

    Source: J. Fluid Mech. 494, 379-398, doi 10.1017/S0022112003006165. The coefficients match the JFM
    abstract and p. 393 (checked 10 Nov 2026), so the values are final. Eo < 40 is where the
    full form matters. Viana report most of their data within 10% of this fit.
    Their R is sqrt(D^3 g drho rho_l) / mu_l, which is this project's Nf. At large Eo the form reduces
    to Llewellin (2.5).
    """
    if Eo <= 6:
        raise ValueError(f"Eo = {Eo} is outside the data behind Viana eq (23), which needs Eo > 6.")
    num = 0.34 / (1 + 3805 / Eo ** 3.06) ** 0.58
    k = 1 + 7.22e13 / Eo ** 9.93
    inner = (Nf / (31.08 * (1 + 778.76 / Eo ** 1.96) ** 0.49)) ** (-1.45 * k ** 0.094)
    return num / (1 + inner) ** (0.71 * k ** (-0.094))


def fr_llewellin(Nf):
    """Froude number from Llewellin et al. (2012), eq (2.5), valid for Eo > 40."""
    return 0.34 * (1 + (31.08 / Nf) ** 1.45) ** (-0.71)


def film_development_length(Fr, Nf):
    """Length below the nose, in D, before the film reaches equilibrium: Llewellin et al. (2012), eq (2.27).

    Z / D = 0.0852 (Fr^2 Nf)^(2/3). Araujo et al. (2012), eq (7), use a stricter test and get more, about 2.9 D
    at Nf 1000, so report both when the film window is short.
    """
    return 0.0852 * (Fr ** 2 * Nf) ** (2.0 / 3.0)


def film_llewellin(Nf):
    """Film thickness over pipe RADIUS (not diameter) from Llewellin et al. (2012), eq (4.2).

    Valid for Eo > 40 and 0.1 < Nf < 1e5. The film is laminar below Nf about 3000 and the tail shape
    is unstable above Nf about 600.
    """
    if not 0.1 < Nf < 1e5:
        raise ValueError(f"Nf = {Nf} is outside 0.1 < Nf < 1e5, the range of Llewellin eq (4.2).")
    return 0.204 + 0.123 * math.tanh(2.66 - 1.15 * math.log10(Nf))


if __name__ == "__main__":
    print("inertial limit:", fr_inertial())
    print(f"{'Nf':>6} {'delta/D':>8} {'Re_film':>8}")
    for Nf in (10, 50, 200, 500, 1000, 2000):
        d = film_thickness(FR_DUMITRESCU, Nf)
        print(f"{Nf:6g} {d:8.4f} {film_reynolds(d, Nf):8.1f}")
