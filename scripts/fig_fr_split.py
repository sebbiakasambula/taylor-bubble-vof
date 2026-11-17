"""Fr over 5 s windows, T1: the two panels of fig_fr_t1 as two single-panel figures, 3.4 in wide.

Data, colours, markers, line styles, legend entries, labels, tick values and the dotted Viana line (Fr = 0.3384) are
those of the two-panel figure made in figs.py (section "Fig 3"). Fonts are the same point sizes (9 pt title, label and
tick labels, 7.5 pt legend). Nothing is written into the paper folder; the output goes to the work folder ($P2_WORK).

Run from any folder:  python fig_fr_split.py
"""
import sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from config import WORK
S = P3 = WORK
OF, FL, REF, GRID = '#2a78d6', '#eb6834', '#6b6b66', '#e4e3dd'
BL = ['#86b6ef', '#2a78d6', '#104281']   # coarse, medium, fine
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 1.6,
                     'savefig.dpi': 300, 'legend.frameon': False})
W_NEW = 3.4          # inches, width of each new figure
FIGW_OLD, FIGH_OLD = 6.3, 2.8   # the original two-panel figure


def fit(t, y, a, b):
    m = (t >= a) & (t <= b + 1e-3); return np.polyfit(t[m], y[m], 1)[0]


def win(t, y, end):
    W = np.arange(5, end, 2.5)
    return np.array([(w + 2.5, fit(t, y, w, w + 5)) for w in W if ((t >= w) & (t <= w + 5)).sum() >= 3])


o = {k: np.load(S / (k + '.npy')) for k in ['V4_coarse', 'V4_medium', 'V4_fine', 'VA_square']}
f = {k: np.load(S / (k + '.npy')) for k in ['f_t1', 'f_t1sq']}


def draw_left(a, legend_kw):
    for k, c, l in zip(['V4_coarse', 'V4_medium', 'V4_fine'], BL, ['Coarse', 'Medium', 'Fine']):
        w = win(o[k][:, 0], o[k][:, 1], 37); a.plot(w[:, 0], w[:, 1], '-o', c=c, ms=3.5, label=l)
    a.set_title('OpenFOAM, stretched meshes', fontsize=9)
    a.set_xlabel('Window centre, $t$ (s)'); a.axhline(0.3384, c=REF, lw=1, ls=':')
    a.set_ylabel('Fr over a 5 s window'); a.legend(fontsize=7.5, **legend_kw)


def draw_right(a, legend_kw):
    for d, t, y, c, ls, l in [(o['V4_medium'], 0, 1, OF, '-', 'OF, stretched'), (o['VA_square'], 0, 1, OF, '--', 'OF, square'),
                              (f['f_t1'], 1, 2, FL, '-', 'F, stretched'), (f['f_t1sq'], 1, 2, FL, '--', 'F, square')]:
        w = win(d[:, t], d[:, y], 17.5); sq = (ls == '--')
        a.plot(w[:, 0], w[:, 1], ls, marker='s' if sq else 'o', mfc='white' if sq else c, c=c, ms=4.5 if sq else 3.5, label=l)
    a.set_title('Medium stretched and square meshes', fontsize=9)
    a.set_xlabel('Window centre, $t$ (s)'); a.axhline(0.3384, c=REF, lw=1, ls=':')
    a.set_ylabel('Fr over a 5 s window')      # the original has no y label here (sharey with the left panel)
    a.legend(fontsize=7.5, handlelength=3, **legend_kw)


# Legend placement per panel (found by the overlap test in checks/A12_check.py).
# left: free upper right corner, as in the original. right: the original 2-column legend (2.3 in wide) does not fit
# into the 2.4 in wide axes, so it becomes one column in the free upper right corner, pushed to 0.1 of the font size
# from the axes corner to keep clear of the F, square markers at x = 12.5 s.
LEG_L = dict(loc='upper right')
LEG_R = dict(loc='upper right', ncol=1, borderaxespad=0.1)
# the legend calls of the original two-panel figure, used only to rebuild its layout (axes size) for the reference
ORIG_L = {}
ORIG_R = dict(ncol=2, loc='center right', bbox_to_anchor=(1, 0.37))


# ---- the shared y range and the tick values of the original pair (two panels, sharey), built only to read them
def reference_pair():
    fig, ax = plt.subplots(1, 2, figsize=(FIGW_OLD, FIGH_OLD), sharey=True)
    draw_left(ax[0], ORIG_L); draw_right(ax[1], ORIG_R)
    ax[1].set_ylabel('')                     # the original right panel has no y label
    fig.tight_layout(); fig.canvas.draw()
    ref = dict(ylim=ax[0].get_ylim(), yticks=ax[0].get_yticks(),
               xlim=[a.get_xlim() for a in ax], xticks=[a.get_xticks() for a in ax],
               ratio=[(a.get_position().height * FIGH_OLD) / (a.get_position().width * FIGW_OLD) for a in ax])
    plt.close(fig)
    return ref


REF_PAIR = reference_pair()


def make(draw, legend_kw, idx):
    """One panel in a 3.4 in wide figure; height chosen so that the axes (plot area) h/w equals the original panel's."""
    ratio = REF_PAIR['ratio'][idx]
    H = 2.8
    for _ in range(4):
        fig, a = plt.subplots(figsize=(W_NEW, H))
        draw(a, legend_kw)
        # ticks first, then limits (set_xticks would widen the limits to include a tick outside them)
        a.set_xticks(REF_PAIR['xticks'][idx]); a.set_xlim(REF_PAIR['xlim'][idx])
        a.set_yticks(REF_PAIR['yticks']); a.set_ylim(REF_PAIR['ylim'])
        fig.tight_layout(); fig.canvas.draw()
        pos = a.get_position()
        w_ax = pos.width * W_NEW; h_ax = pos.height * H
        margins = H - h_ax                     # top + bottom, in inches (independent of H for fixed text)
        H_new = ratio * w_ax + margins
        if abs(H_new - H) < 1e-6:
            break
        plt.close(fig); H = H_new
    return fig, a


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    for draw, leg, idx, name in [(draw_left, LEG_L, 0, 'fig_fr_t1_left'), (draw_right, LEG_R, 1, 'fig_fr_sq')]:
        fig, a = make(draw, leg, idx)
        fig.savefig(P3 / (name + '.pdf'))
        fig.savefig(P3 / (name + '_preview.png'), dpi=110)
        pos = a.get_position(); Hf = fig.get_size_inches()[1]
        print(name, 'fig %.4f x %.4f in' % (W_NEW, Hf), 'axes %.4f x %.4f in' % (pos.width * W_NEW, pos.height * Hf),
              'h/w %.5f (original %.5f)' % (pos.height * Hf / (pos.width * W_NEW), REF_PAIR['ratio'][idx]))
        plt.close(fig)
