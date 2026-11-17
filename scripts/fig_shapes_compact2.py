"""Interface shapes at gas fraction 0.5 (T1, T2, T3), OpenFOAM and Fluent: compact version 2 of fig_shapes.pdf, at most 2.70 in tall.

Same data, colours, line widths, wall lines, legend entries, panel titles, tick labels and y range (-5 to 0.2) as
fig_shapes_compact.py (which holds the data code copied from figs.py, section "Fig 1 shapes"); only the height is
2.70 in instead of 3.30 in, so with equal scale (1 D in r = 1 D in x) the three panels get narrower. The titles stay on
one line (9 pt); the width is the smallest at which neighbouring titles keep a 12 pt gap, which is below the 4.2 in limit.
Nothing is written into the paper folder; the output goes to the work folder ($P2_WORK):
    fig_shapes_compact2.pdf, fig_shapes_compact2_preview.png (110 dpi)

Run from any folder:  python fig_shapes_compact2.py
"""
import sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fig_shapes_compact as base      # same folder; the layout functions read base.W and base.H at call time

H2 = 2.70            # inches, height cap (used in full)
W_MAX = 4.2          # inches, width limit
TITLE_GAP = 12.0     # pt, smallest gap between neighbouring titles


def make2(sh):
    base.H = H2
    base.W = W_MAX
    fig, ax, leg = base.make(sh)          # first pass at the width limit: gives c1 and the title widths
    r = fig.canvas.get_renderer()
    tb = [a.title.get_window_extent(r) for a in ax]
    tw = [b.width / fig.dpi for b in tb]
    cs = [0.5 * (a.get_window_extent(r).x0 + a.get_window_extent(r).x1) / fig.dpi for a in ax]
    c1 = cs[0]
    pitch = max(0.5 * (tw[0] + tw[1]), 0.5 * (tw[1] + tw[2])) + TITLE_GAP / 72
    W_new = c1 + 2 * pitch + tw[2] / 2 + base.EDGE
    plt.close(fig)
    base.W = min(W_new, W_MAX)
    return base.make(sh)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    sh = base.load_shapes()
    fig, ax, leg = make2(sh)
    fig.savefig(base.P3 / 'fig_shapes_compact2.pdf')
    fig.savefig(base.P3 / 'fig_shapes_compact2_preview.png', dpi=110)
    Wf, Hf = fig.get_size_inches()
    for i, a in enumerate(ax):
        bb = a.get_position()
        print('panel', i + 1, 'axes %.4f x %.4f in' % (bb.width * Wf, bb.height * Hf), 'centre x %.4f in' % ((bb.x0 + bb.width / 2) * Wf))
    print('figure %.4f x %.4f in' % (Wf, Hf))
    plt.close(fig)
