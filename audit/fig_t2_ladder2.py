"""Paper figure fig_t2_tests.eps, two panels (17 Nov): (a) Fr and (b) film thickness at T2 over 5 to 8.4 s for every
configuration and test, with the standard error of each fit or mean as a bar.

Same rows and sampling as fig_tests.py: OpenFOAM runs on all their snapshots in the window, Fluent tests on the
snapshots nearest the five Fluent baseline times. References: Abubakar and Matar at Nf 50 (Fr 0.2531, range 0.2514 to
0.2544, their Fig. 6a; film 0.1403, range 0.1401 to 0.1405, their Fig. 7b, abu/fig7/digit7b.py), Viana et al. (Fr),
White and Beardmore at Eo 128 (Fr), Llewellin et al. (film). EPS holds no transparency, so bands are opaque.
"""
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(HERE.parent / "scripts"))  # shared modules and config.py
from config import WORK as FIG  # noqa: E402  (the EPS files land in the work folder; copy them to paper/figures)
sys.path.insert(0, str(HERE.parent / "scripts"))
from tb_references import fr_viana, film_llewellin  # noqa: E402

OF, FL, REF, GRID, BAND = '#2a78d6', '#eb6834', '#6b6b66', '#e4e3dd', '#dededd'
plt.rcParams.update({'font.size': 8.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'savefig.dpi': 300,
                     'legend.frameon': False, 'ps.fonttype': 42})
LO, HI = 5.0 - 1e-6, 8.4 + 1e-3


def load(p):
    rows = list(csv.DictReader(open(p)))
    def num(v):
        try:
            return float(v)
        except ValueError:
            return float("nan")
    return {k: np.array([num(r[k]) for r in rows]) for k in ("t", "nose", "film")}


def fit(t, y):
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    return float(p[0]), math.sqrt(np.sum(r ** 2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2))


def mean_se(v):
    v = v[np.isfinite(v)]
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))


def all_rows(M):
    w = (M["t"] >= LO) & (M["t"] <= HI) & np.isfinite(M["nose"])
    return fit(M["t"][w], M["nose"][w]), mean_se(M["film"][w])


def matched(M, tb):
    i = [int(np.argmin(np.abs(M["t"] - x))) for x in tb]
    return fit(M["t"][i], M["nose"][i]), mean_se(M["film"][i])


def M(n):
    return load(HERE / "m" / ("%s.csv" % n))


FT2 = M("F_T2")
tbF = FT2["t"][(FT2["t"] >= LO) & (FT2["t"] <= HI)]
rows = [("OpenFOAM baseline", all_rows(M("T2_prod")), OF),
        ("baseline, half the time step", all_rows(M("T2_A9_halfdt")), OF),
        (u"baseline, 10° wedge", all_rows(M("T2_B1_wedge10")), OF),
        ("baseline, curvature smoothing 2", all_rows(M("T2_B3_smooth2")), OF),
        ("baseline, density-weighted surface tension", all_rows(M("T2_rhoST")), OF),
        ("momentum predictor", all_rows(M("T2_mompred")), OF),
        ("predictor and 3 outer loops", all_rows(M("T2_A3_outer3")), OF),
        ("predictor and isoAdvector, isoAlpha", all_rows(M("T2_A1_iso")), OF),
        ("predictor and isoAdvector, plicRDF", all_rows(M("T2_A2_plic")), OF),
        ("Fluent baseline", all_rows(FT2), FL),
        ("Fluent, Green-Gauss node based gradient", matched(load(HERE / "flcopy/m_F_T2_A6.csv"), tbF), FL),
        ("Fluent, Modified HRIC", matched(load(HERE / "flcopy/m_F_T2_A7.csv"), tbF), FL),
        ("Fluent, global Courant number 0.25", matched(load(HERE / "flcopy/m_F_T2_A8.csv"), tbF), FL)]
SEPS = (5, 9)


def panel(ax, vals, xlim, fmt, step):
    n = len(rows)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    for i, ((v, se), c) in enumerate(zip(vals, [r[2] for r in rows])):
        y = n - 1 - i
        ax.text(1.02, y, fmt % v, transform=tr, va="center", ha="left", fontsize=7, color="#333")
        if v < xlim[0]:
            ax.annotate("", xy=(xlim[0], y), xytext=(xlim[0] + step, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.1))
            ax.text(xlim[0] + 1.15 * step, y, "off scale", va="center", ha="left", fontsize=7, color="#333")
            continue
        ax.errorbar(v, y, xerr=se, fmt="o", ms=5, color=c, mec="white", mew=0.9, ecolor=c, elinewidth=1.0, capsize=2,
                    zorder=3)
    for s in SEPS:
        ax.axhline(n - 1 - s + 0.5, color='#bbb', lw=0.6, zorder=1)
    ax.set_yticks(range(n))
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, n - 0.4)
    ax.grid(axis="y", visible=False)


fig = plt.figure(figsize=(6.6, 4.15))
axA = fig.add_axes([0.405, 0.105, 0.228, 0.745])
axB = fig.add_axes([0.718, 0.105, 0.196, 0.745])
fv = fr_viana(100, 50)
fl = 0.5 * film_llewellin(50)
axA.axvspan(0.2514, 0.2544, color=BAND, lw=0, zorder=0)
axA.axvline(0.2531, color=REF, lw=1.0, zorder=1)
axA.axvline(fv, color=REF, lw=1.0, ls="--", zorder=1)
axA.axvline(0.2671, color=REF, lw=1.2, ls=":", zorder=1)
axB.axvspan(0.1401, 0.1405, color=BAND, lw=0, zorder=0)
axB.axvline(0.1403, color=REF, lw=1.0, zorder=1)
axB.axvline(fl, color=REF, lw=1.0, ls="-.", zorder=1)
panel(axA, [r[1][0] for r in rows], (0.205, 0.278), "%.4f", 0.0045)
panel(axB, [r[1][1] for r in rows], (0.1305, 0.1445), "%.4f", 0.0009)
axA.set_yticklabels([r[0] for r in rows][::-1], fontsize=7.4)
axB.set_yticklabels([])
axA.set_xlabel("Fr")
axB.set_xlabel(r"film thickness $\delta/D$")
axA.set_xticks([0.21, 0.23, 0.25, 0.27])
axB.set_xticks([0.132, 0.136, 0.140, 0.144])
for ax, lab in ((axA, "(a)"), (axB, "(b)")):
    ax.text(0.0, 1.015, lab, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5)
handles = [Line2D([], [], marker="o", ls="", color=OF, mec="white", ms=5.5, label="OpenFOAM"),
           Line2D([], [], marker="o", ls="", color=FL, mec="white", ms=5.5, label="Fluent"),
           Patch(facecolor=BAND, edgecolor=REF, lw=0.8, label="Abubakar and Matar, Nf 50"),
           Line2D([], [], color=REF, lw=1.0, ls="--", label="Viana et al. (a)"),
           Line2D([], [], color=REF, lw=1.2, ls=":", label="White and Beardmore, Eo 128 (a)"),
           Line2D([], [], color=REF, lw=1.0, ls="-.", label="Llewellin et al. (b)")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=3, fontsize=7.3, handlelength=1.9,
           columnspacing=1.3, borderaxespad=0.2)
fig.savefig(FIG / "fig_t2_tests.eps")
fig.savefig(HERE / "fig_t2_ladder2_preview.png", dpi=130)
for lab, ((fr, se), (fi, sf)), _ in rows:
    print("%-44s Fr %.4f +- %.4f   film %.4f +- %.5f" % (lab, fr, se, fi, sf))
print("Viana %.4f, Llewellin film %.4f" % (fv, fl))
