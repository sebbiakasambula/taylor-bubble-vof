"""Paper figure fig_transfer.eps (17 Nov): the same settings at T1, T2 and T3, Fr with the standard error of each fit.

Replaces fig_t3_tests.eps (its T3 rows are the right panel here, same values). Windows: T1 5 to 10 s (the length of
the predictor and geometric advection runs at T1, set 14), T2 and T3 5 to 8.4 s. Sampling as in the paper: OpenFOAM
runs on all snapshots in the window, except at T3, where both OpenFOAM runs are measured at the highest point of the
largest connected gas region (a5/*.npy) and the geometric advection run at the snapshots nearest the five baseline times;
Fluent tests at the snapshots nearest the Fluent baseline times.
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

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(HERE.parent / "scripts"))  # shared modules and config.py
from config import WORK as FIG  # noqa: E402  (the EPS files land in the work folder; copy them to paper/figures)
sys.path.insert(0, str(HERE.parent / "scripts"))
from tb_references import fr_viana  # noqa: E402

OF, FL, REF, GRID, BAND = '#2a78d6', '#eb6834', '#6b6b66', '#e4e3dd', '#dededd'
plt.rcParams.update({'font.size': 8.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'savefig.dpi': 300,
                     'legend.frameon': False, 'ps.fonttype': 42})


def load(p):
    rows = list(csv.DictReader(open(p)))
    def num(v):
        try:
            return float(v)
        except ValueError:
            return float("nan")
    return {k: np.array([num(r[k]) for r in rows]) for k in ("t", "nose")}


def fit(t, y):
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    return float(p[0]), math.sqrt(np.sum(r ** 2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2))


def win(M, lo, hi):
    return (M["t"] >= lo - 1e-6) & (M["t"] <= hi + 1e-3) & np.isfinite(M["nose"])


def fr_all(M, lo, hi):
    w = win(M, lo, hi)
    return fit(M["t"][w], M["nose"][w])


def fr_matched(M, tb):
    i = [int(np.argmin(np.abs(M["t"] - x))) for x in tb]
    return fit(M["t"][i], M["nose"][i])


def m(n):
    return load(HERE / "m" / ("%s.csv" % n))


ROWS = ["OpenFOAM baseline", "momentum predictor", "geometric advection and predictor", "Fluent baseline",
        "Fluent, Green-Gauss node based gradient", "Fluent, Modified HRIC"]
COL = [OF, OF, OF, FL, FL, FL]

# T1, 5 to 10 s
t1 = {0: fr_all(m("V4_medium"), 5, 10), 1: fr_all(m("T1med_mompred"), 5, 10), 2: fr_all(m("T1med_A2plic"), 5, 10),
      3: fr_all(m("F_T1"), 5, 10)}
# T2, 5 to 8.4 s
FT2 = m("F_T2")
tb2 = FT2["t"][win(FT2, 5, 8.4)]
t2 = {0: fr_all(m("T2_prod"), 5, 8.4), 1: fr_all(m("T2_mompred"), 5, 8.4), 2: fr_all(m("T2_A2_plic"), 5, 8.4),
      3: fr_all(FT2, 5, 8.4), 4: fr_matched(load(HERE / "flcopy/m_F_T2_A6.csv"), tb2),
      5: fr_matched(load(HERE / "flcopy/m_F_T2_A7.csv"), tb2)}
# T3, 5 to 8.4 s
P = np.load(HERE / "a5" / "T3_prod.npy")
A5 = np.load(HERE / "a5" / "T3_A5_plic.npy")
i5 = [int(np.argmin(np.abs(A5[:, 0] - x))) for x in P[:, 0]]
FT3 = m("F_T3")
tb3 = FT3["t"][win(FT3, 5, 8.4)]
t3 = {0: fit(P[:, 0], P[:, 1]), 2: fit(A5[i5, 0], A5[i5, 1]), 3: fr_all(FT3, 5, 8.4),
      4: fr_matched(load(HERE / "flcopy/m_F_T3_D1.csv"), tb3), 5: fr_matched(load(HERE / "flcopy/m_F_T3_D2.csv"), tb3)}

fig = plt.figure(figsize=(6.6, 2.75))
x0, w, gap = 0.372, 0.178, 0.034
axes = [fig.add_axes([x0 + k * (w + gap), 0.16, w, 0.60]) for k in range(3)]
cases = [("(a) T1, 5 to 10 s", t1, (0.300, 0.350), fr_viana(300, 1000)),
         ("(b) T2, 5 to 8.4 s", t2, (0.205, 0.275), fr_viana(100, 50)),
         ("(c) T3, 5 to 8.4 s", t3, (0.262, 0.340), fr_viana(20, 500))]
n = len(ROWS)
for ax, (title, vals, xlim, fv) in zip(axes, cases):
    if "T2" in title:
        ax.axvspan(0.2514, 0.2544, color=BAND, lw=0, zorder=0)
        ax.axvline(0.2531, color=REF, lw=1.0, zorder=1)
        ax.axvline(0.2671, color=REF, lw=1.2, ls=":", zorder=1)
    ax.axvline(fv, color=REF, lw=1.0, ls="--", zorder=1)
    for i in range(n):
        y = n - 1 - i
        if i in vals:
            v, se = vals[i]
            ax.errorbar(v, y, xerr=se, fmt="o", ms=5, color=COL[i], mec="white", mew=0.9, ecolor=COL[i],
                        elinewidth=1.0, capsize=2, zorder=3)
        else:
            ax.text(0.5, y, "not run", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=6.8,
                    color="#888")
    ax.axhline(n - 1 - 3 + 0.5, color='#bbb', lw=0.6, zorder=1)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(range(n))
    ax.set_yticklabels([])
    ax.grid(axis="y", visible=False)
    ax.set_title(title, fontsize=8.3, loc="left", pad=3)
    ax.set_xlabel("Fr")
    ax.tick_params(axis="x", labelsize=7.4)
axes[0].set_yticklabels(ROWS[::-1], fontsize=7.4)
axes[0].set_xticks([0.30, 0.32, 0.34])
axes[1].set_xticks([0.21, 0.23, 0.25, 0.27])
axes[2].set_xticks([0.27, 0.29, 0.31, 0.33])
handles = [Line2D([], [], marker="o", ls="", color=OF, mec="white", ms=5.5, label="OpenFOAM"),
           Line2D([], [], marker="o", ls="", color=FL, mec="white", ms=5.5, label="Fluent"),
           Line2D([], [], color=REF, lw=1.0, ls="--", label="Viana et al."),
           Patch(facecolor=BAND, edgecolor=REF, lw=0.8, label="Abubakar and Matar, Nf 50 (T2)"),
           Line2D([], [], color=REF, lw=1.2, ls=":", label="White and Beardmore, Eo 128 (T2)")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=3, fontsize=7.3, handlelength=1.9,
           columnspacing=1.3, borderaxespad=0.2)
fig.savefig(FIG / "fig_transfer.eps")
fig.savefig(HERE / "fig_transfer_preview.png", dpi=130)
for (title, vals, xlim, fv) in cases:
    base = vals[0][0]
    print(title, "Viana %.4f" % fv)
    for i, (v, se) in sorted(vals.items()):
        print("   %-42s %.4f +- %.4f  (%+.2f%% vs OF baseline, %+.2f%% vs Viana)" % (ROWS[i], v, se, 100 * (v / base - 1), 100 * (v / fv - 1)))
