"""Fr and film thickness by case (T1 to T3), OpenFOAM and Fluent: compact version of fig_cases.pdf, 6.3 in wide, at most 2.05 in tall.

The numbers are computed as in figs.py (header and "numbers" section: Fr as the slope fit of the front position over
5 to 20 s, film thickness averaged over 5 to 20 s, OpenFOAM T1 film re-measured by of_film with the wedge-wall shift
applied), the plot calls are those of the section "Fig 2 Fr and film by case": same markers, sizes, colours, error bars
(Viana et al., +-0.05), Llewellin marks, axis limits, tick labels and legend entries. Only the height changes (2.6 in to
2.05 in). figs.py is not run and nothing is written into the paper folder. Output in the work folder ($P2_WORK):
    fig_cases_compact.pdf, fig_cases_compact_preview.png (110 dpi)

Run from any folder:  python fig_cases_compact.py
"""
import sys
import re
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from config import RUNS, WORK
S = P3 = WORK
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])      # scal, vec, tb_post, as in figs.py
RUN = RUNS
OF, FL, REF, GRID = '#2a78d6', '#eb6834', '#6b6b66', '#e4e3dd'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 1.6,
                     'savefig.dpi': 300, 'legend.frameon': False})
W, H = 6.3, 2.05     # inches: width as before, height cap
# Legends: same entries and font size (7.5 pt) as the original, but shorter handles and less padding, because in the
# 1.5 in tall axes the original left legend (2 columns, 174 pt) is wider than the axes (165 pt) and the right one
# comes within 0.4 pt of the T2 Fluent marker. Gaps from the legend (box and handles) to the nearest data item: see checks/A13b_fig.md.
LEG_L = dict(loc='lower left', ncol=2, handlelength=1.5, columnspacing=1.0, handletextpad=0.5, borderaxespad=0.2, borderpad=0.2)
LEG_R = dict(loc='upper right', labelspacing=0.3, borderaxespad=0.2, borderpad=0.2)       # default handle length: the Llewellin mark is 22 pt wide and needs it


def fit(t, y, a, b):
    m = (t >= a) & (t <= b + 1e-3); return np.polyfit(t[m], y[m], 1)[0]


def film(t, f, a, b):
    m = (t >= a) & (t <= b + 1e-3); return np.nanmean(f[m])


def of_film(case, a0, b0, nc):
    C = vec(RUN / case / '0/C'); x, y = C[:, 0], np.round(C[:, 1], 9); v = []
    for d in sorted([p for p in (RUN / case).iterdir() if re.fullmatch(r'[0-9.]+', p.name)], key=lambda p: float(p.name)):
        t = float(d.name)
        if a0 <= t <= b0 + 1e-3:
            a = scal(d / 'alpha.water'); n = tb_post.nose_position(x, y, a); tl = tb_post.tail_position(x, y, a)
            xi, ri = tb_post.interface_radius(x, y, a); w = (xi > tl + 0.5) & (xi < n - nc); v.append((0.5 - ri[w]).mean())
    return np.mean(v)


def numbers():
    o = {k: np.load(S / (k + '.npy')) for k in ['V4_medium', 'T2_prod', 'T3_prod']}
    f = {k: np.load(S / (k + '.npy')) for k in ['f_t1', 'f_t2', 'f_t3']}
    ofk = ['V4_medium', 'T2_prod', 'T3_prod']; flk = ['f_t1', 'f_t2', 'f_t3']
    frO = [fit(o[k][:, 0], o[k][:, 1], 5, 20) for k in ofk]; frF = [fit(f[k][:, 1], f[k][:, 2], 5, 20) for k in flk]
    fiO = [film(o[k][:, 0], o[k][:, 4], 5, 20) for k in ofk]; fiF = [film(f[k][:, 1], f[k][:, 4], 5, 20) for k in flk]
    fiO[0] = of_film('V4_medium', 5, 20, 2.0)
    # OpenFOAM films are measured from the flat wedge wall, which lies at 0.5 cos(2.5 deg) in the mid-plane, not 0.5
    WALL_SHIFT = 0.5 - 0.5 * np.cos(np.radians(2.5)); fiO = [v - WALL_SHIFT for v in fiO]
    return frO, frF, fiO, fiF


def draw(frO, frF, fiO, fiF, legend_kw=None):
    cases = ['T1', 'T2', 'T3']; X = np.arange(3); vi = [0.3384, 0.2504, 0.2770]; ll = [0.0615, 0.1394, np.nan]
    fig, ax = plt.subplots(1, 2, figsize=(W, H))
    ax[0].errorbar(X, vi, yerr=0.05, fmt='_', c=REF, ms=14, capsize=4, lw=1.2, label=r'Viana et al., $\pm0.05$')
    ax[0].plot(X - 0.12, frO, 'o', c=OF, ms=7, label='OpenFOAM'); ax[0].plot(X + 0.12, frF, 's', c=FL, ms=7, label='Fluent')
    ax[0].set_ylabel('Fr, fit over 5 to 20 s')
    ax[0].set_yticks(np.linspace(0.10, 0.40, 7)); ax[0].set_ylim(0.10, 0.40)   # the seven tick values of the original (the auto locator would drop to four in the shorter axes); ticks first, then limits
    ax[1].plot(X[:2], ll[:2], '_', c=REF, ms=22, mew=2, zorder=0, label='Llewellin et al.')
    ax[1].plot(X - 0.12, fiO, 'o', c=OF, ms=7, label='OpenFOAM'); ax[1].plot(X + 0.12, fiF, 's', c=FL, ms=7, label='Fluent')
    ax[1].set_ylabel(r'Film thickness, $\delta/D$'); ax[1].set_ylim(0, 0.24)
    for a in ax:
        a.set_xticks(X, cases); a.set_xlim(-0.5, 2.5); a.grid(axis='x', visible=False)
    kw0, kw1 = legend_kw if legend_kw else (LEG_L, LEG_R)
    ax[0].legend(fontsize=7.5, **kw0); ax[1].legend(fontsize=7.5, **kw1)
    return fig, ax


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    frO, frF, fiO, fiF = numbers()
    print('Fr OF', np.round(frO, 4), 'Fl', np.round(frF, 4)); print('film OF', np.round(fiO, 4), 'Fl', np.round(fiF, 4))
    fig, ax = draw(frO, frF, fiO, fiF)
    fig.tight_layout()
    fig.savefig(P3 / 'fig_cases_compact.pdf'); fig.savefig(P3 / 'fig_cases_compact_preview.png', dpi=110)
    print('figure %.3f x %.3f in' % tuple(fig.get_size_inches()))
    plt.close(fig)
