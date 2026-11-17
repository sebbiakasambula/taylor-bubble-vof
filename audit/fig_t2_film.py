"""Paper figure fig_t2_film.eps (17 Nov): the film at T2 against laminar film theory.

(a) Downward axial speed of the liquid cells across the film, in the middle column of the film window of the last
snapshot of each run (set17_profiles.py, profiles.npz), with Brown's profile for the film of that column (eq. 3.5 of
Abubakar and Matar). (b) Fr against the film thickness over 5 to 8.4 s for every T2 run, with the Fr that a developed
laminar film of that thickness carries (set17_film_theory.fr_exact), the steady results of Abubakar and Matar
interpolated to Nf 50, and the correlations of Viana et al. (Fr) and Llewellin et al. (film).
"""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(HERE.parent / "scripts"))  # shared modules and config.py
from config import WORK as FIG  # noqa: E402  (the EPS files land in the work folder; copy them to paper/figures)
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "scripts"))
import set17_film_theory as T  # noqa: E402
from fig_t2_ladder2 import rows as LADDER  # noqa: E402  (same values as the ladder figure)

OF, OF2, OF3, FL, REF, GRID = '#2a78d6', '#86b6ec', '#123f73', '#eb6834', '#6b6b66', '#e4e3dd'
plt.rcParams.update({'font.size': 8.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'savefig.dpi': 300,
                     'legend.frameon': False, 'ps.fonttype': 42})
G = T.G_of("T2")
P = np.load(HERE / "profiles.npz")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.5, 3.15), gridspec_kw=dict(width_ratios=[1, 1.08], wspace=0.32))
runs = [("OF_baseline", "OpenFOAM baseline", OF, "o"), ("OF_predictor", "OpenFOAM predictor", OF2, "s"),
        ("OF_geometric_and_predictor", "OpenFOAM geometric and predictor", OF3, "^"), ("Fluent", "Fluent", FL, "D")]
handles = []
for key, lab, c, mk in runs:
    t, xr, delta, R = P[key + "_meta"]
    r, w = P[key + "_r"], P[key + "_w"]
    a = R - delta
    liq = r > a
    rr = np.linspace(a, R, 200)
    axA.plot(R - rr, G / 4 * (R * R - rr * rr) - G * a * a / 2 * np.log(R / rr), color=c, lw=1.0, zorder=2)
    axA.plot(R - r[liq], w[liq], mk, ms=3.2, color=c, mec="white", mew=0.4, zorder=3)
    handles.append(Line2D([], [], marker=mk, ls="", color=c, mec="white", ms=5, label=lab))
handles.append(Line2D([], [], color="#333", lw=1.0, label="Brown's profile, same film"))
axA.set_xlabel("distance from the wall (D)")
axA.set_ylabel("downward axial speed (m/s)")
axA.set_xlim(0, 0.15)
axA.set_ylim(0, 0.47)

# (b) Fr against film
d = np.linspace(0.115, 0.148, 300)
axB.plot(d, [T.fr_exact(x, 0.5, G) for x in d], color="#333", lw=1.0, zorder=1)
for lab, ((fr, se), (fi, sf)), c in LADDER:
    axB.plot(fi, fr, "o", ms=4.6, color=c, mec="white", mew=0.6, zorder=3)
for f, c in (("m/T2_C1_coarse.csv", OF), ("m/T2_C2_fine.csv", OF), ("m/T2_M1c_mompred_coarse.csv", OF),
             ("m/T2_A2c_coarse.csv", OF), ("flcopy/m_F_T2_coarse.csv", FL), ("flcopy/m_F_T2_fine.csv", FL)):
    M = T.load(HERE / f)
    fr, se, n = T.fr_fit(M, 5, 8.4)
    fi, sd, nf = T.film_mean(M, 5, 8.4)
    axB.plot(fi, fr, "o", ms=4.6, mfc="white", color=c, mew=1.0, zorder=3)
axB.plot(0.1403, 0.2531, "D", ms=5, color=REF, mec="white", mew=0.6, zorder=4)
axB.plot(0.1394, 0.2504, "x", ms=5.5, color="#111", mew=1.2, zorder=4)
axB.set_xlabel(r"film thickness $\delta/D$")
axB.set_ylabel("Fr")
axB.set_xlim(0.118, 0.146)
axB.set_ylim(0.145, 0.285)
handles = [h for h in handles if h.get_label() != "Brown's profile, same film"]
for h in handles:
    h.set_label(h.get_label() + " (a)")
hB = [Line2D([], [], marker="o", ls="", color=OF, mec="white", ms=5, label="OpenFOAM, medium mesh (b)"),
      Line2D([], [], marker="o", ls="", color=FL, mec="white", ms=5, label="Fluent, medium mesh (b)"),
      Line2D([], [], marker="o", ls="", mfc="white", color="#555", mew=1.0, ms=5, label="coarse and fine meshes (b)"),
      Line2D([], [], marker="D", ls="", color=REF, mec="white", ms=5, label="Abubakar and Matar, Nf 50 (b)"),
      Line2D([], [], marker="x", ls="", color="#111", mew=1.2, ms=5.5, label="Viana et al. with Llewellin et al. (b)"),
      Line2D([], [], color="#333", lw=1.0, label="developed laminar film (a, b)")]
for ax, lab in ((axA, "(a)"), (axB, "(b)")):
    ax.text(-0.02, 1.02, lab, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5)
fig.legend(handles=handles + hB, loc="upper left", bbox_to_anchor=(0.005, 0.995), ncol=3, fontsize=6.7,
           handlelength=1.6, columnspacing=1.0, borderaxespad=0.2)
fig.subplots_adjust(left=0.085, right=0.985, bottom=0.13, top=0.72)
fig.savefig(FIG / "fig_t2_film.eps")
fig.savefig(HERE / "fig_t2_film_preview.png", dpi=130)
for key, lab, c, mk in runs:
    t, xr, delta, R = P[key + "_meta"]
    print("%s: t %.4f, %.4f D behind the nose, film %.4f, wall %.6f" % (lab, t, -xr, delta, R))
