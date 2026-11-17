"""Interface shapes at gas fraction 0.5 (T1, T2, T3), OpenFOAM and Fluent: compact version of fig_shapes.pdf.

Same data as the figure made in figs.py (section "Fig 1 shapes"): same snapshots, same lines (each interface mirrored
about the axis), same colours, line widths, wall lines, legend entries, panel titles and tick labels. Only the layout
changes: 4.2 in wide, at most 3.3 in tall, so it is printed at 100 % without scaling in LaTeX. The three panels share
the y axis (tick labels on the first panel only), 1 D in r is 1 D in x in every panel, and the axes are placed in
inches so that the plot area, the text and the empty margins come out as measured, not as tight_layout guesses.

Nothing is written into the paper folder (figs.py is not run); the output goes to the work folder ($P2_WORK):
    fig_shapes_compact.pdf, fig_shapes_compact_preview.png (110 dpi)

Run from any folder:  python fig_shapes_compact.py
"""
import sys
import re
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from config import RUNS, FLUENT, WORK
S = P3 = WORK
# helper functions of snap.py (scal, vec, tb_post), the part above "case=Path" as figs.py does
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])
FD = FLUENT
RUN = RUNS
OF, FL, REF, GRID = '#2a78d6', '#eb6834', '#6b6b66', '#e4e3dd'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 1.6,
                     'savefig.dpi': 300, 'legend.frameon': False})

# ---- the figure
W, H = 4.2, 3.3            # inches: width target, height cap (the height is used in full, so the plot area is as large as it can be)
XLIM = (-0.55, 0.55)       # r/D: the data lie in -0.4834 .. 0.4834 (walls at +-0.5)
YLIM = (-5.0, 0.2)         # x/D: see the note at the end of this block
# The lowest plotted point of each line (x/D): T1 OpenFOAM -6.50, T1 Fluent -10.06 (gas shed into the wake),
# T2 -3.26 / -3.43, T3 -3.42 / -3.56. The original figure cuts every panel at -5, which keeps the T1 tail (-3.1 to -4.6)
# and the first wake fragments the caption talks about; the lines are clipped there exactly as in the original.
# A lower limit of -3.7 would show T2 and T3 whole but would cut the T1 tail off, so -5 stays. The upper limit goes
# from 0.3 to 0.2 (the noses are at -0.002 .. -0.022, the line is about 0.05 D thick).
FS_TITLE, FS_LABEL, FS_XTICK, FS_LEG = 9, 9, 8, 8      # pt, as in the original figure (y tick labels: default 9 pt)
GAP_T = 0.03               # in, empty strip between the legend and the panel titles
EDGE = 0.02                # in, empty strip between the outermost text and the canvas edge


def of_shape(case, t):
    C = vec(RUN / case / '0/C'); a = scal(RUN / case / t / 'alpha.water'); x, y = C[:, 0], np.round(C[:, 1], 9)
    xi, ri = tb_post.interface_radius(x, y, a); return xi - tb_post.nose_position(x, y, a), ri


def fl_shape(folder, name):
    d = np.loadtxt(FD / folder / name, delimiter=',', skiprows=1); x, y, g = d[:, 1], np.round(d[:, 2], 9), d[:, 9]
    xi, ri = tb_post.interface_radius(x, y, 1 - g); return xi - tb_post.nose_position(x, y, 1 - g), ri


def load_shapes():
    return [(of_shape('V4_medium', '20.315988'), fl_shape('T1_Case', 's_snap-21645'), 'T1, Eo 300, Nf 1000'),
            (of_shape('T2_prod', '19.967232'), fl_shape('T2_case', 'T2_Taylor_bubble-23760'), 'T2, Eo 100, Nf 50'),
            (of_shape('T3_prod', '20.305026'), fl_shape('T3_Case', 'T3_Taylor_bubble-78000'), 'T3, Eo 20, Nf 500')]


def extent_in(fig, artist):
    """Bounding box of a text/legend artist in inches from the lower left canvas corner: (x0, y0, x1, y1)."""
    bb = artist.get_window_extent(fig.canvas.get_renderer())
    return bb.x0 / fig.dpi, bb.y0 / fig.dpi, bb.x1 / fig.dpi, bb.y1 / fig.dpi


def draw_panels(sh):
    """Create the figure with three axes (placed later by place()); returns fig, axes, legend."""
    fig = plt.figure(figsize=(W, H))
    ax = []
    for i, (so, sf, tt) in enumerate(sh):
        a = fig.add_axes([0.1 + 0.3 * i, 0.1, 0.2, 0.6], sharey=ax[0] if ax else None)
        for (xx, rr), c, l in [(so, OF, 'OpenFOAM'), (sf, FL, 'Fluent')]:
            k = np.argsort(xx); xx, rr = xx[k], rr[k]
            a.plot(rr, xx, c=c, label=l); a.plot(-rr, xx, c=c)
        a.axvline(-0.5, c=REF, lw=1.2); a.axvline(0.5, c=REF, lw=1.2)
        a.set_xlim(*XLIM); a.set_ylim(*YLIM)
        a.set_title(tt, fontsize=FS_TITLE, pad=4)
        a.set_xlabel('$r/D$', fontsize=FS_LABEL, labelpad=2)
        a.set_xticks([-0.5, 0, 0.5], ['-0.5', '0', '0.5'], fontsize=FS_XTICK)
        a.tick_params(axis='x', pad=5)                    # 1.5 pt more than the default, keeps -0.5 clear of the -5 of the y axis
        a.set_aspect('equal'); a.grid(False)
        if i > 0:
            a.tick_params(axis='y', labelleft=False)      # y tick labels on the first panel only
        ax.append(a)
    ax[0].set_ylabel('Distance from the nose, $x/D$', fontsize=FS_LABEL, labelpad=2)
    h, l = ax[0].get_legend_handles_labels()
    leg = fig.legend(h, l, loc='upper center', ncol=2, fontsize=FS_LEG, borderaxespad=0, borderpad=0.2,
                     bbox_to_anchor=(0.5, 1 - EDGE / H))
    return fig, ax, leg


def place(fig, ax, leg, left, bottom, h, c1, c3):
    """Axes placed in inches: bottom edge at `bottom`, height h, width h * 1.1 / 5.2 (equal scale), centres from c1 to c3."""
    w = h * (XLIM[1] - XLIM[0]) / (YLIM[1] - YLIM[0])
    cs = [c1, 0.5 * (c1 + c3), c3]
    for a, c in zip(ax, cs):
        a.set_position([(c - w / 2) / W, bottom / H, w / W, h / H])
    leg.set_bbox_to_anchor((cs[1] / W, 1 - EDGE / H), transform=fig.transFigure)
    return w, cs


def make(sh):
    fig, ax, leg = draw_panels(sh)
    # provisional layout, then three measure-and-place passes (text sizes do not depend on the layout)
    h, bottom, c1, c3 = 2.3, 0.5, 0.9, 3.4
    for it in range(4):
        w, cs = place(fig, ax, leg, 0, bottom, h, c1, c3)
        fig.canvas.draw()
        ab = [extent_in(fig, a) for a in ax]                      # axes boxes in inches
        # left: y label (and y tick labels, which it clears) against the left canvas edge
        yl = extent_in(fig, ax[0].yaxis.label)
        left_need = ab[0][0] - yl[0] + EDGE                       # distance from the first axes' left edge to the canvas edge
        # bottom: x label of the lowest of the three
        xl_low = min(extent_in(fig, a.xaxis.label)[1] for a in ax)
        bottom_need = ab[0][1] - xl_low + EDGE
        # top: legend height, strip, tallest title (measured from the axes top edge)
        lg = extent_in(fig, leg)
        title_top = max(extent_in(fig, a.title)[3] for a in ax) - ab[0][3]
        top_need = EDGE + (lg[3] - lg[1]) + GAP_T + title_top
        # title widths (T1 first, T3 last) set the outer panel centres
        tw = [extent_in(fig, a.title)[2] - extent_in(fig, a.title)[0] for a in ax]
        h_new = H - top_need - bottom_need
        w_new = h_new * (XLIM[1] - XLIM[0]) / (YLIM[1] - YLIM[0])
        c1_new = max(left_need + w_new / 2, tw[0] / 2 + EDGE)
        c3_new = W - EDGE - tw[2] / 2
        if max(abs(h_new - h), abs(bottom_need - bottom), abs(c1_new - c1), abs(c3_new - c3)) < 1e-6:
            break
        h, bottom, c1, c3 = h_new, bottom_need, c1_new, c3_new
    place(fig, ax, leg, 0, bottom, h, c1, c3)
    fig.canvas.draw()
    return fig, ax, leg


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    sh = load_shapes()
    fig, ax, leg = make(sh)
    fig.savefig(P3 / 'fig_shapes_compact.pdf')
    fig.savefig(P3 / 'fig_shapes_compact_preview.png', dpi=110)
    Wf, Hf = fig.get_size_inches()
    for i, a in enumerate(ax):
        bb = a.get_position()
        print('panel', i + 1, 'axes %.4f x %.4f in' % (bb.width * Wf, bb.height * Hf), 'centre x %.4f in' % ((bb.x0 + bb.width / 2) * Wf))
    print('figure %.4f x %.4f in' % (Wf, Hf))
    plt.close(fig)
