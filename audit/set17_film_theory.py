"""Set 17. Measured Fr and film against laminar film theory, every run at T1, T2 and T3 (17 Nov).

The film below a long bubble is taken as fully developed and laminar with no shear from the gas. Its axial velocity is
Brown's profile, eq. (3.5) of Abubakar and Matar,
    w(r) = (G/4)(R^2 - r^2) - (G a^2/2) ln(R/r),   G = drho g / mu_l,  a = R - delta,
and the liquid it carries down must equal the volume the bubble sweeps up (Llewellin et al., Q_f = -Q_b), so
    Fr_film = Q / (pi a^2) = G / (8 a^2) [R^4 - 4 R^2 a^2 + 3 a^4 + 4 a^4 ln(R/a)].
Brown's own relation (Llewellin eq. 2.14, the cubic of their eq. 2.25 and Abubakar and Matar eq. 3.3) uses the flux of a
flat film over the interface circumference, Fr_B = (2G/3) delta^3 / (R - delta). Both are printed.

Fr is the straight-line fit of the nose over the window, the film the mean over the snapshots of the window, from the
per-snapshot tables in audit/m (same numbers as sets 2, 5, 6, 7, 8, 9 and 14). R is the wall the film is measured from,
the wedge mid-plane 0.499524 in OpenFOAM and 0.5 in Fluent.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from tb_params import fluid  # noqa: E402

WALL_OF = 0.5 * math.cos(math.radians(2.5))
CASE = {"T1": (300, 1000), "T2": (100, 50), "T3": (20, 500)}


def load(p):
    rows = list(csv.DictReader(open(p)))
    def num(v):
        try:
            return float(v)
        except ValueError:
            return float("nan")
    return {k: np.array([num(r[k]) for r in rows]) for k in ("t", "nose", "film")}


def window(M, lo, hi):
    return (M["t"] >= lo - 1e-6) & (M["t"] <= hi + 1e-3)


def fr_fit(M, lo, hi):
    w = window(M, lo, hi) & np.isfinite(M["nose"])
    t, y = M["t"][w], M["nose"][w]
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    se = math.sqrt(np.sum(r ** 2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2)) if len(t) > 2 else float("nan")
    return float(p[0]), se, int(w.sum())


def film_mean(M, lo, hi):
    w = window(M, lo, hi) & np.isfinite(M["film"])
    return float(np.mean(M["film"][w])), float(np.std(M["film"][w])), int(w.sum())


def G_of(case):
    f = fluid(*CASE[case], 1000)
    return (f["rho_l"] - f["rho_g"]) * 1.0 / f["mu_l"]


def fr_exact(delta, R, G):
    a = R - delta
    return G / (8 * a * a) * (R ** 4 - 4 * R * R * a * a + 3 * a ** 4 + 4 * a ** 4 * math.log(R / a))


def fr_brown(delta, R, G):
    return 2 * G / 3 * delta ** 3 / (R - delta)


def delta_exact(fr, R, G):
    return brentq(lambda d: fr_exact(d, R, G) - fr, 1e-4, R - 1e-4)


def flux_check(R, G, delta):
    """Integrate the profile numerically as a check of the closed form."""
    a = R - delta
    r = np.linspace(a, R, 20001)
    w = G / 4 * (R * R - r * r) - G * a * a / 2 * np.log(R / r)
    return np.trapezoid(w * 2 * np.pi * r, r) / (np.pi * a * a)


AUD = HERE
RUNS = [
    # case, label, file, wall, window
    ("T2", "OF baseline", "m/T2_prod.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF baseline, half step", "m/T2_A9_halfdt.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF baseline, 10 deg wedge", "m/T2_B1_wedge10.csv", 0.5 * math.cos(math.radians(5.0)), (5, 8.4)),
    ("T2", "OF baseline, smoothing", "m/T2_B3_smooth2.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF baseline, density weight", "m/T2_rhoST.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF predictor", "m/T2_mompred.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF predictor, 3 outer", "m/T2_A3_outer3.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF iso isoAlpha", "m/T2_A1_iso.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF iso plicRDF", "m/T2_A2_plic.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF baseline coarse", "m/T2_C1_coarse.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF baseline fine", "m/T2_C2_fine.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF predictor coarse", "m/T2_M1c_mompred_coarse.csv", WALL_OF, (5, 8.4)),
    ("T2", "OF iso plicRDF coarse", "m/T2_A2c_coarse.csv", WALL_OF, (5, 8.4)),
    ("T2", "Fluent", "m/F_T2.csv", 0.5, (5, 8.4)),
    ("T2", "Fluent coarse", "flcopy/m_F_T2_coarse.csv", 0.5, (5, 8.4)),
    ("T2", "Fluent fine", "flcopy/m_F_T2_fine.csv", 0.5, (5, 8.4)),
    ("T2", "Fluent node gradient", "flcopy/m_F_T2_A6.csv", 0.5, (5, 8.4)),
    ("T2", "Fluent Modified HRIC", "flcopy/m_F_T2_A7.csv", 0.5, (5, 8.4)),
    ("T2", "Fluent Courant 0.25", "flcopy/m_F_T2_A8.csv", 0.5, (5, 8.4)),
    ("T1", "OF baseline", "m/V4_medium.csv", WALL_OF, (5, 10)),
    ("T1", "OF predictor", "m/T1med_mompred.csv", WALL_OF, (5, 10)),
    ("T1", "OF iso plicRDF", "m/T1med_A2plic.csv", WALL_OF, (5, 10)),
    ("T1", "Fluent", "m/F_T1.csv", 0.5, (5, 10)),
    ("T1", "OF baseline 5-20", "m/V4_medium.csv", WALL_OF, (5, 20)),
    ("T1", "OF fine 5-20", "m/V4_fine.csv", WALL_OF, (5, 20)),
    ("T1", "Fluent 5-20", "m/F_T1.csv", 0.5, (5, 20)),
    ("T3", "OF baseline", "m/T3_prod.csv", WALL_OF, (5, 8.4)),
    ("T3", "OF iso plicRDF", "m/T3_A5_plic.csv", WALL_OF, (5, 8.4)),
    ("T3", "Fluent", "m/F_T3.csv", 0.5, (5, 8.4)),
    ("T3", "OF baseline 5-20", "m/T3_prod.csv", WALL_OF, (5, 20)),
    ("T3", "Fluent 5-20", "m/F_T3.csv", 0.5, (5, 20)),
]

if __name__ == "__main__":
    print("Closed form against numerical integration, T2, delta 0.134: %.6f %.6f" % (
        fr_exact(0.134, 0.5, G_of("T2")), flux_check(0.5, G_of("T2"), 0.134)))
    for c in CASE:
        f = fluid(*CASE[c], 1000)
        print("%s: mu_l %.6g, nu_l %.6g, G %.6g" % (c, f["mu_l"], f["nu_l"], G_of(c)))
    print()
    print("| case | run | window | Fr (s.e., n) | film (sd, n) | Fr_film exact | Fr / Fr_film - 1 (%) | Fr_film Brown | Fr / Fr_B - 1 (%) | film for the measured Fr (exact) | film / that - 1 (%) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    out = []
    for case, lab, f, R, (lo, hi) in RUNS:
        M = load(AUD / f)
        fr, se, n = fr_fit(M, lo, hi)
        d, sd, nf = film_mean(M, lo, hi)
        G = G_of(case)
        fe = fr_exact(d, R, G)
        fb = fr_brown(d, R, G)
        de = delta_exact(fr, R, G)
        out.append((case, lab, fr, d, fe))
        print("| %s | %s | %g to %g | %.4f (%.4f, %d) | %.4f (%.4f, %d) | %.4f | %+.2f | %.4f | %+.2f | %.4f | %+.2f |" % (
            case, lab, lo, hi, fr, se, n, d, sd, nf, fe, 100 * (fr / fe - 1), fb, 100 * (fr / fb - 1), de,
            100 * (d / de - 1)))
    # Abubakar and Matar at T2: Fr 0.2531 (Fig. 6a), film 0.1403 (Fig. 7b), both interpolated to Nf 50 at Eo 100.
    # Their gas is passive, so their Nf has rho_l, not drho: G = rho_l g / mu = Nf.
    for lab, fr, d in (("A&M (Nf 50, interpolated)", 0.2531, 0.1403),):
        G = 50.0
        print("| T2 | %s | steady | %.4f | %.4f | %.4f | %+.2f | %.4f | %+.2f | %.4f | %+.2f |" % (
            lab, fr, d, fr_exact(d, 0.5, G), 100 * (fr / fr_exact(d, 0.5, G) - 1), fr_brown(d, 0.5, G),
            100 * (fr / fr_brown(d, 0.5, G) - 1), delta_exact(fr, 0.5, G), 100 * (d / delta_exact(fr, 0.5, G) - 1)))
    # their markers at Nf 40, 60, 80 (Eo 100): Fr from Fig. 6a, film from Fig. 7b
    for nf, fr, d in ((40, 0.2336, 0.1458), (60, 0.2692, 0.1351), (80, 0.2879, 0.1266)):
        G = float(nf)
        print("| T2 | A&M Nf %d | steady | %.4f | %.4f | %.4f | %+.2f | %.4f | %+.2f | %.4f | %+.2f |" % (
            nf, fr, d, fr_exact(d, 0.5, G), 100 * (fr / fr_exact(d, 0.5, G) - 1), fr_brown(d, 0.5, G),
            100 * (fr / fr_brown(d, 0.5, G) - 1), delta_exact(fr, 0.5, G), 100 * (d / delta_exact(fr, 0.5, G) - 1)))
    # references at T2 and T1
    for case, lab, fr, d in (("T2", "Viana with Llewellin", 0.2504, 0.1394), ("T1", "Viana with Llewellin", 0.3384, 0.0615)):
        G = G_of(case)
        print("| %s | %s | correlations | %.4f | %.4f | %.4f | %+.2f | %.4f | %+.2f | %.4f | %+.2f |" % (
            case, lab, fr, d, fr_exact(d, 0.5, G), 100 * (fr / fr_exact(d, 0.5, G) - 1), fr_brown(d, 0.5, G),
            100 * (fr / fr_brown(d, 0.5, G) - 1), delta_exact(fr, 0.5, G), 100 * (d / delta_exact(fr, 0.5, G) - 1)))
