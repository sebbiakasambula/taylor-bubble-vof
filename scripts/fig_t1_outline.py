"""T1 bubble outline, stretched against square mesh, both codes, in one strip (alpha = 0.5 contour of the main bubble)."""
import math, sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import contourpy
from pypdf import PdfReader

from config import RUNS, FLUENT, WORK
SP = P3 = WORK
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])      # scal, vec, tb_post
RUN = RUNS
FD = FLUENT
W = 0.5 - 0.5 * math.cos(math.radians(2.5))                  # OpenFOAM wedge: moves the wall to r = 0.5

# name, code, case folder, snapshot, flow time in s
fl_str_t = float(np.load(SP / 'f_t1.npy')[:, :2][np.load(SP / 'f_t1.npy')[:, 0] == 21645][0, 1])
_g = np.loadtxt(FD / 'T1_square_w_completion/gas-centroid-rfile_6001.out', skiprows=3)
fl_sq_t = float(_g[_g[:, 0] == 8190][0, 2])
SNAP = [('OpenFOAM stretched', 'of', 'V4_medium', '19.700352', 19.700352),
        ('OpenFOAM square', 'of', 'VA_square', '19.700352', 19.700352),
        ('Fluent stretched', 'fl', 'T1_Case', 's_snap-21645', fl_str_t),
        ('Fluent square', 'fl', 'T1_square_w_completion', 's_snap_completion-8190', fl_sq_t)]
COL = {'of': '#2a78d6', 'fl': '#eb6834'}
LW = 1.4


def load(kind, case, snap):
    if kind == 'of':
        C = vec(RUN / case / '0/C'); a = scal(RUN / case / snap / 'alpha.water')
        return C[:, 0], C[:, 1], a, W
    d = np.loadtxt(FD / case / snap, delimiter=',', skiprows=1)
    return d[:, 1], d[:, 2], 1 - d[:, 9], 0.0


def area(p):
    return 0.5 * np.sum(p[:-1, 0] * p[1:, 1] - p[1:, 0] * p[:-1, 1])


def outline(name, kind, case, snap, t):
    x, r, a, w = load(kind, case, snap)
    nose = tb_post.nose_position(x, np.round(r, 9), a)
    X = np.unique(np.round(x, 6)); R = np.unique(np.round(r, 6)); nx, nr = len(X), len(R)
    assert nx * nr == len(x), (name, nx, nr, len(x))
    ix = np.searchsorted(X, np.round(x, 6)); ir = np.searchsorted(R, np.round(r, 6))
    A = np.zeros((nr, nx)); A[ir, ix] = a
    Rs = R + w
    Rm = np.r_[-Rs[::-1], Rs]; Am = np.vstack([A[::-1], A])      # mirrored about the axis
    lines = contourpy.contour_generator(X, Rm, Am).lines(0.5)
    closed = [p for p in lines if len(p) > 3 and np.allclose(p[0], p[-1])]
    ar = np.array([abs(area(p)) for p in closed])
    order = np.argsort(ar)[::-1]
    loop = closed[order[0]].copy(); loop[:, 0] -= nose
    k = np.searchsorted(X, nose) - 1
    return dict(name=name, kind=kind, t=t, nose=nose, loop=loop, area=ar[order[0]], area2=ar[order[1]] if len(ar) > 1 else 0.0,
                nclosed=len(closed), nopen=len(lines) - len(closed), dx_nose=float(X[k + 1] - X[k]), nx=nx, nr=nr,
                is_closed=bool(np.allclose(loop[0], loop[-1])))


def half(loop, upper):
    """Part of a closed loop with r >= 0 (upper) or r <= 0, as polylines that start and end on the axis."""
    p = loop[:-1]; n = len(p); pts = []
    for i in range(n):
        a, b = p[i], p[(i + 1) % n]
        pts.append(a)
        if (a[1] < 0) != (b[1] < 0) and a[1] != 0 and b[1] != 0:
            s = a[1] / (a[1] - b[1]); pts.append(np.array([a[0] + s * (b[0] - a[0]), 0.0]))
    pts = np.array(pts)
    inside = pts[:, 1] >= 0 if upper else pts[:, 1] <= 0
    if inside.all(): return [np.vstack([pts, pts[:1]])]
    first_out = np.argmin(inside)                      # rotate to start on an outside point
    pts = np.roll(pts, -first_out, axis=0); inside = np.roll(inside, -first_out)
    runs, cur = [], []
    for q, ok in zip(pts, inside):
        if ok: cur.append(q)
        elif cur: runs.append(np.array(cur)); cur = []
    if cur: runs.append(np.array(cur))
    return runs


res = [outline(*s) for s in SNAP]
for o in res:
    L = o['loop']
    o['pieces'] = len(half(L, o['kind'] == 'of')); o['xmin'] = L[:, 0].min(); o['xmax'] = L[:, 0].max(); o['rmax'] = np.abs(L[:, 1]).max()
    print(f"{o['name']:20s} t {o['t']:.4f} nose {o['nose']:.4f} loop x {o['xmin']:.3f}..{o['xmax']:.4f} rmax {o['rmax']:.4f} "
          f"area {o['area']:.3f} (2nd {o['area2']:.4f}) closed {o['nclosed']} open {o['nopen']} dx_nose {o['dx_nose']:.4f} grid {o['nx']}x{o['nr']}")

with open(P3 / 't1_outline_snap.csv', 'w') as f:
    f.write('name,time_s,nose_x,loop_min_x_minus_nose,loop_max_r\n')
    for o in res: f.write(f"{o['name']},{o['t']:.4f},{o['nose']:.4f},{o['xmin']:.4f},{o['rmax']:.4f}\n")

# ---------------- figure ----------------
plt.rcParams.update({'font.size': 8.5})
loop_min = min(o['xmin'] for o in res)
XL = (-3.9 if loop_min >= -3.75 else math.floor((loop_min - 0.15) * 10) / 10 - 0.0, 0.3)
if len(sys.argv) > 1: XL = (float(sys.argv[1]), 0.3)
YL = (-0.55, 0.55)
ASP = (YL[1] - YL[0]) / (XL[1] - XL[0])
TARGET_W, PAD, WCORR = 6.3, 0.02, 0.006
LEG_H = 0.30                                       # space above the axes for the legend row


def build(wax, wc, hc, ox, oy):
    fig = plt.figure(figsize=(wc, hc)); hax = wax * ASP
    ax = fig.add_axes([ox / wc, oy / hc, wax / wc, hax / hc])
    ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect('equal')
    ax.axhline(0, color='#6b6b66', lw=0.6, ls=(0, (1, 2)), zorder=1)
    for s in (-0.5, 0.5): ax.axhline(s, color='#6b6b66', lw=1.0, zorder=1)
    for o in res:
        upper = o['kind'] == 'of'
        sq = 'square' in o['name']
        for seg in half(o['loop'], upper):
            ln, = ax.plot(seg[:, 0], seg[:, 1], color=COL[o['kind']], lw=LW, zorder=3, solid_capstyle='butt',
                          **({'dashes': (4, 2)} if sq else {}))
    ax.set_xlabel('Distance from the nose, $x/D$'); ax.set_ylabel('$r/D$')
    ax.set_yticks([-0.5, 0, 0.5]); ax.set_xticks(np.arange(math.ceil(XL[0]), 1))
    xt = ax.text(XL[0] + 0.06, 0.28, 'OpenFOAM', color=COL['of'], ha='left', va='center', fontsize=8.5)
    ft = ax.text(XL[0] + 0.06, -0.28, 'Fluent', color=COL['fl'], ha='left', va='center', fontsize=8.5)
    h = [Line2D([], [], color='k', lw=LW), Line2D([], [], color='k', lw=LW, dashes=(4, 2))]
    leg = ax.legend(h, ['stretched', 'square'], loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2, frameon=False,
                    fontsize=8.5, handlelength=2.6, columnspacing=1.6, borderaxespad=0.15, borderpad=0.2)
    return fig, ax, leg, (xt, ft)


def tight(fig):
    fig.canvas.draw()
    return fig.get_tightbbox(fig.canvas.get_renderer())


wax, ox, oy = 5.4, 1.0, 1.0
for it in range(8):
    fig, ax, leg, tx = build(wax, 9, 6, ox, oy)
    tb = tight(fig); plt.close(fig)
    dw = (TARGET_W + WCORR - 2 * PAD) - tb.width
    ox, oy = ox - tb.x0 + PAD, oy - tb.y0 + PAD
    wax += dw
    print(f'iter {it}: content {tb.width:.4f} x {tb.height:.4f} -> wax {wax:.4f}')
    if abs(dw) < 1e-4: break
# the final canvas is the tight box (height from the last pass)
fig, ax, leg, tx = build(wax, 9, 6, ox, oy)
tb = tight(fig); plt.close(fig)
fig, ax, leg, tx = build(wax, tb.width + 2 * PAD, tb.height + 2 * PAD, ox - tb.x0 + PAD, oy - tb.y0 + PAD)
fig.canvas.draw(); R_ = fig.canvas.get_renderer()

# ---------------- checks ----------------
texts = [('xlabel', ax.xaxis.label), ('ylabel', ax.yaxis.label), ('text OpenFOAM', tx[0]), ('text Fluent', tx[1])]
texts += [(f'xtick {t.get_text()}', t) for t in ax.xaxis.get_majorticklabels()] + [(f'ytick {t.get_text()}', t) for t in ax.yaxis.get_majorticklabels()]
texts += [(f'legend {t.get_text()}', t) for t in leg.get_texts()]
texts = [(n, t) for n, t in texts if t.get_visible() and t.get_text().strip()]
boxes = [(n, t, t.get_window_extent(R_)) for n, t in texts]


def sep(a, b):
    dx = max(b.x0 - a.x1, a.x0 - b.x1); dy = max(b.y0 - a.y1, a.y0 - b.y1)
    return math.hypot(dx, dy) if dx > 0 and dy > 0 else max(dx, dy)


pairs = sorted((sep(boxes[i][2], boxes[j][2]), boxes[i][0], boxes[j][0]) for i in range(len(boxes)) for j in range(i + 1, len(boxes)))
overlaps = [p for p in pairs if p[0] < 0]
Wf, Hf = fig.bbox.width, fig.bbox.height
outside = [n for n, t, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > Wf + 0.5 or b.y1 > Hf + 0.5]
sizes = sorted({round(t.get_fontsize(), 2) for n, t, b in boxes})
dpi = fig.dpi
# plotted line vertices inside the two in-axes texts (display coordinates, 2 px margin)
hits = {}
for n, t, b in boxes:
    if n.startswith('text'):
        c = 0
        for ln in ax.lines:
            xy = ax.transData.transform(np.c_[ln.get_xdata(), ln.get_ydata()]) if len(ln.get_xdata()) > 2 else np.empty((0, 2))
            c += int(np.sum((xy[:, 0] > b.x0 - 2) & (xy[:, 0] < b.x1 + 2) & (xy[:, 1] > b.y0 - 2) & (xy[:, 1] < b.y1 + 2)))
        hits[n] = c
bb = ax.get_window_extent(R_); xr = ax.get_xlim(); yr = ax.get_ylim()
aspect = ((bb.height / abs(yr[1] - yr[0])) / (bb.width / abs(xr[1] - xr[0])))

fig.savefig(P3 / 'fig_t1_outline.pdf', bbox_inches='tight', pad_inches=PAD, dpi=300)
fig.savefig(P3 / 'fig_t1_outline_preview.png', bbox_inches='tight', pad_inches=PAD, dpi=110)
mb = PdfReader(str(P3 / 'fig_t1_outline.pdf')).pages[0].mediabox
pw, ph = float(mb.width) / 72, float(mb.height) / 72

Lm = ['# B3 check: fig_t1_outline', '',
      f'Script: `fig_t1_outline.py`. Texts checked: {len(boxes)} visible text artists (x label, y label, the two in-axes labels, tick labels, legend entries).', '',
      '## 1. Text overlap test (after fig.canvas.draw())', '',
      f'- Intersecting pairs: **{len(overlaps)}**',
      f'- Texts extending outside the figure canvas: **{len(outside)}**' + ('' if not outside else ' ' + ', '.join(outside)),
      '- Closest non-intersecting pairs (separation in pt): ' + '; '.join(f'{p[1]} / {p[2]} = {p[0] / dpi * 72:.1f}' for p in pairs[:5]),
      '- Plotted outline vertices inside the in-axes texts (2 px margin): ' + ', '.join(f'{k} = {v}' for k, v in hits.items()),
      '', '## 2. Font sizes', '',
      f'- Smallest font size used: **{min(sizes):.1f} pt**; all sizes present: {", ".join(f"{s:.1f}" for s in sizes)} pt',
      '', '## 3. Aspect', '', f'- Axes set_aspect(\'equal\'); data aspect measured from the axes box: {aspect:.4f} ({(aspect - 1) * 100:+.2f} %)',
      f'- x limits {XL[0]:.2f} to {XL[1]:.2f}, r limits {YL[0]} to {YL[1]}; leftmost point of any kept loop {loop_min:.3f}',
      '', '## 4. Kept loops (main bubble = closed alpha = 0.5 loop of largest enclosed area on the mirrored grid)', '',
      '| snapshot | time (s) | nose x | grid (x by r) | closed loops | open lines | kept area | next area | closed | min x - nose | max x - nose | axial cell at nose | within one cell | max r | pieces drawn |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
ok_all = True
for o in res:
    within = abs(o['xmax']) <= o['dx_nose']; ok_all &= within and o['is_closed']
    Lm.append(f"| {o['name']} | {o['t']:.4f} | {o['nose']:.4f} | {o['nx']} by {o['nr']} | {o['nclosed']} | {o['nopen']} | {o['area']:.4f} | {o['area2']:.4f} | "
              f"{'yes' if o['is_closed'] else 'NO'} | {o['xmin']:.4f} | {o['xmax']:.4f} | {o['dx_nose']:.4f} | {'yes' if within else 'NO'} | {o['rmax']:.4f} | {o['pieces']} |")
Lm += ['', f'- All four kept loops closed and rightmost point within one axial cell of the nose: **{"yes" if ok_all else "NO"}**',
       '', '## 5. Final figure size from the PDF (MediaBox)', '',
       f'- **{pw:.3f} in wide by {ph:.3f} in tall** ({float(mb.width):.1f} x {float(mb.height):.1f} pt); axes {wax:.3f} in wide by {wax * ASP:.3f} in tall']
(P3 / 'checks/B3_fig.md').write_text('\n'.join(Lm) + '\n', encoding='utf-8')
print('\n'.join(Lm))
