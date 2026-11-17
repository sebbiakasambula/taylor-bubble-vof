"""The code comparison on the medium mesh over 5 to 20 s as two figures (the reference panels side by side,
the rest as a group of four, every marker explained in a legend above the plots).

fig_cases.eps: (a) Fr with Viana et al. and their stated scatter, (b) film with Llewellin et al.
fig_cases2.eps: (a) bubble length with its range, (b) gas outside the bubble and the part within 2 D below the tail,
(c) smeared gas, (d) largest speed in all cells and in the liquid.
Values from audit/cmp/cmp_table.csv, as fig_cmp.py (the six-panel version, kept for reverting). No alpha (EPS).
"""
import csv
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
OF, FL, REF, GRID, INK = '#2a78d6', '#eb6834', '#6b6b66', '#e4e3dd', '#333333'
plt.rcParams.update({'font.size': 8.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'savefig.dpi': 300,
                     'legend.frameon': False, 'ps.fonttype': 42})

V = {}
for r in csv.DictReader(open(HERE / "cmp" / "cmp_table.csv")):
    V[(r["case"], r["code"], r["statistic"])] = float(r["value"])
CASES = ["T1", "T2", "T3"]
X = np.arange(3)


def get(code, stat):
    return np.array([V.get((c, code, stat), np.nan) for c in CASES])


VIANA = [0.3384, 0.2504, 0.2770]
LLEW = [0.0615, 0.1394, np.nan]
D = 0.12
D2 = (0.20, 0.07)
MS = 6
CODES = [Line2D([], [], ls='', marker='o', c=OF, ms=MS, label='OpenFOAM'),
         Line2D([], [], ls='', marker='s', c=FL, ms=MS, label='Fluent')]


def finish(ax, letter):
    ax.set_xticks(X, CASES)
    ax.set_xlim(-0.5, 2.5)
    ax.grid(axis='x', visible=False)
    ax.set_title("(%s)" % letter, loc='left', fontsize=8.5, pad=3)


def key(ax, handles):
    ax.legend(handles=handles, loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=len(handles), fontsize=7.3,
              handletextpad=0.3, columnspacing=0.9, borderaxespad=0.2)


# Figure 7: the two quantities with reference values
fig, (a, b) = plt.subplots(1, 2, figsize=(6.3, 2.35))
a.errorbar(X, VIANA, yerr=0.05, fmt='_', c=REF, ms=12, capsize=3.5, lw=1.1)
a.plot(X - D, get("OF", "Fr"), 'o', c=OF, ms=MS)
a.plot(X + D, get("F", "Fr"), 's', c=FL, ms=MS)
a.set_ylabel("Fr")
a.set_ylim(0.15, 0.40)
b.plot(X[:2], LLEW[:2], '_', c=REF, ms=18, mew=2, zorder=0)
b.plot(X - D, get("OF", "film_mean"), 'o', c=OF, ms=MS)
b.plot(X + D, get("F", "film_mean"), 's', c=FL, ms=MS)
b.set_ylabel(r"Film thickness, $\delta/D$")
b.set_ylim(0, 0.16)
finish(a, "a")
finish(b, "b")
h = CODES + [Line2D([], [], ls='', marker='_', c=REF, ms=12, mew=1.1, label=r'Viana et al., $\pm0.05$'),
             Line2D([], [], ls='', marker='_', c=REF, ms=14, mew=2, label='Llewellin et al.')]
fig.legend(handles=h, loc='upper center', ncol=4, fontsize=7.3, handletextpad=0.4, columnspacing=2.0,
           bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.0)
fig.savefig(FIG / "fig_cases.eps")
fig.savefig(HERE / "fig_cases_preview.png", dpi=110)
plt.close(fig)

# Figure 8: the four quantities without a reference value
fig, axs = plt.subplots(2, 2, figsize=(6.3, 4.5))
c, d, e, f = axs.ravel()
for code, dx, mk, col in (("OF", -D, 'o', OF), ("F", D, 's', FL)):
    c.vlines(X + dx, get(code, "L_b_min"), get(code, "L_b_max"), color=col, lw=1.2)
    c.plot(X + dx, get(code, "L_b_mean"), mk, c=col, ms=MS)
c.set_ylabel(r"Bubble length, $L/D$")
c.set_ylim(1.8, 4.0)
key(c, [Line2D([], [], ls='', marker='o', c=INK, ms=MS - 1, label='mean'),
        Line2D([], [], ls='', marker='|', c=INK, ms=10, mew=1.2, label='range over the window')])

for code, s, mk, col in (("OF", -1, 'o', OF), ("F", 1, 's', FL)):
    d.plot(X + s * D2[0], get(code, "out_b_mean"), mk, c=col, ms=MS)
    d.plot(X + s * D2[1], get(code, "wake_b_mean"), mk, c=col, mfc='white', mew=1.2, ms=MS - 0.5)
d.set_ylabel("Gas outside the bubble (%)")
d.set_ylim(-0.6, 16)
key(d, [Line2D([], [], ls='', marker='o', c=INK, ms=MS - 1, label='all'),
        Line2D([], [], ls='', marker='o', c=INK, mfc='white', mew=1.1, ms=MS - 1, label=r'within 2 $D$ below the tail')])

e.plot(X - D, get("OF", "smear_mean"), 'o', c=OF, ms=MS)
e.plot(X + D, get("F", "smear_mean"), 's', c=FL, ms=MS)
e.set_ylabel("Smeared gas (%)")
e.set_ylim(-0.4, 10)

for code, s, mk, col in (("OF", -1, 'o', OF), ("F", 1, 's', FL)):
    f.plot(X + s * D2[0], get(code, "pk_all_median"), mk, c=col, ms=MS)
    f.plot(X + s * D2[1], get(code, "pk_liq_median"), mk, c=col, mfc='white', mew=1.2, ms=MS - 0.5)
f.set_ylabel("Largest speed (m/s)")
f.set_ylim(0, 5.6)
key(f, [Line2D([], [], ls='', marker='o', c=INK, ms=MS - 1, label='all cells'),
        Line2D([], [], ls='', marker='o', c=INK, mfc='white', mew=1.1, ms=MS - 1, label='liquid only')])

for k, ax in enumerate(axs.ravel()):
    finish(ax, "abcd"[k])
fig.legend(handles=CODES, loc='upper center', ncol=2, fontsize=7.3, handletextpad=0.3, columnspacing=1.1,
           bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=1.6, w_pad=2.0)
fig.savefig(FIG / "fig_cases2.eps")
fig.savefig(HERE / "fig_cases2_preview.png", dpi=110)
for st in ("Fr", "film_mean", "L_b_mean", "L_b_min", "L_b_max", "out_b_mean", "wake_b_mean", "smear_mean",
           "pk_all_median", "pk_liq_median"):
    print("%-14s OF %s  F %s" % (st, np.round(get("OF", st), 4), np.round(get("F", st), 4)))
