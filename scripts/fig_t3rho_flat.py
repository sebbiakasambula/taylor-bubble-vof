"""T3 speed field and interface in the pipe, x/D from 0 to 9: four times of the density-weighted-surface-tension run and the production run at 8.27 s."""
import math, sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import Normalize; import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from pypdf import PdfReader

from config import PV as PVDIR, WORK
PV = str(PVDIR) + '/'   # VTK inputs and ParaView strips
P3 = str(WORK) + '/'
EXT = (0.0, 9.0, -0.55, 0.55)
RHO = 'T3, OpenFOAM, density-weighted surface tension, '
PANELS = [('T3R_1',   '(a) ' + RHO + '$t$ = 1.0 s',  'T3_rhoST', '1'),
          ('T3R_24',  '(b) ' + RHO + '$t$ = 2.4 s',  'T3_rhoST', '2.4'),
          ('T3R_5',   '(c) ' + RHO + '$t$ = 5.0 s',  'T3_rhoST', '5'),
          ('T3R_84',  '(d) ' + RHO + '$t$ = 8.4 s',  'T3_rhoST', '8.4'),
          ('T3P_827', '(e) T3, OpenFOAM, $t$ = 8.27 s', 'T3_prod', '8.272418')]
N = len(PANELS); TAGS = 'abcde'
TARGET_W, PAD = 6.3, 0.02                   # inches, final PDF width and tight-bbox pad
GAP = 0.31                                  # inches between strips (same as the T1 figure)
WCORR = 0.03                                # PDF text extents are a little narrower than Agg ones, so build this much wider
CB_GAP, CB_W = 0.12, 0.14                   # inches, gap axes to colour bar and its width

plt.rcParams.update({'font.size': 8.5})
IMGS = {p[0]: mpimg.imread(PV + 'strips/' + p[0] + '.png') for p in PANELS}   # strips from paraview/render_t3rho_flat.py
ASP = (EXT[3] - EXT[2]) / (EXT[1] - EXT[0])   # strip height over width, so the data aspect is 1


def build(wax, gap, wc, hc, ox, oy):
    """Strips placed in absolute inches: ox = left edge of the axes, oy = bottom edge of the lowest strip."""
    fig = plt.figure(figsize=(wc, hc))
    hax = wax * ASP
    axs = []
    for i, (n, t, *_) in enumerate(PANELS):
        yb = oy + (N - 1 - i) * (hax + gap)
        a = fig.add_axes([ox / wc, yb / hc, wax / wc, hax / hc], sharex=axs[0] if axs else None)
        a.imshow(IMGS[n], extent=EXT, aspect='equal')
        a.set_title(t, fontsize=8, loc='left')
        a.set_ylabel('$r/D$'); a.set_yticks([-0.5, 0, 0.5]); a.set_xticks(range(0, 10))
        if i < N - 1: a.tick_params(labelbottom=False)
        axs.append(a)
    axs[-1].set_xlabel('$x/D$')
    c = fig.add_axes([(ox + wax + CB_GAP) / wc, oy / hc, CB_W / wc, ((N - 1) * (hax + gap) + hax) / hc])
    cb = fig.colorbar(cm.ScalarMappable(Normalize(0, 2), 'viridis'), cax=c, extend='max')
    cb.set_ticks([0, 0.5, 1, 1.5, 2]); c.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    cb.set_label('Speed (m/s)')
    return fig, axs, c


def tight(fig):
    fig.canvas.draw()
    return fig.get_tightbbox(fig.canvas.get_renderer())


# solve for the axes width so the cropped figure is TARGET_W wide; the height follows from the strip aspect
wax, ox, oy = 4.85, 1.5, 1.5
for it in range(8):
    fig, axs, c = build(wax, GAP, 9, 12, ox, oy)
    tb = tight(fig); plt.close(fig)
    dw = (TARGET_W + WCORR - 2 * PAD) - tb.width
    ox, oy = ox - tb.x0 + PAD, oy - tb.y0 + PAD
    wax += dw
    print(f'iter {it}: content {tb.width:.4f} x {tb.height:.4f}  -> wax {wax:.4f}')
    if abs(dw) < 1e-4: break
fig, axs, c = build(wax, GAP, 9, 12, ox, oy)
tb = tight(fig); plt.close(fig)
hc = tb.height + 2 * PAD
ox, oy = ox - tb.x0 + PAD, oy - tb.y0 + PAD

fig, axs, cax = build(wax, GAP, TARGET_W + WCORR, hc, ox, oy)
fig.canvas.draw(); R = fig.canvas.get_renderer()

# ---------------- checks ----------------
texts = []   # (label, Text artist)
for k, a in enumerate(axs):
    tag = TAGS[k]
    texts.append((f'title ({tag})', a._left_title))
    texts.append((f'ylabel ({tag})', a.yaxis.label))
    texts.append((f'xlabel ({tag})', a.xaxis.label))
    for t in a.xaxis.get_majorticklabels() + a.xaxis.get_minorticklabels(): texts.append((f'xtick ({tag}) {t.get_text()}', t))
    for t in a.yaxis.get_majorticklabels() + a.yaxis.get_minorticklabels(): texts.append((f'ytick ({tag}) {t.get_text()}', t))
    texts.append((f'xoffset ({tag})', a.xaxis.offsetText)); texts.append((f'yoffset ({tag})', a.yaxis.offsetText))
texts.append(('colourbar label', cax.yaxis.label))
for t in cax.yaxis.get_majorticklabels() + cax.yaxis.get_minorticklabels(): texts.append((f'colourbar tick {t.get_text()}', t))
texts.append(('colourbar xlabel', cax.xaxis.label))
texts = [(n, t) for n, t in texts if t.get_visible() and t.get_text().strip()]
boxes = [(n, t, t.get_window_extent(R)) for n, t in texts]


def sep(a, b):   # separation of two boxes in display px, negative if they intersect
    dx = max(b.x0 - a.x1, a.x0 - b.x1); dy = max(b.y0 - a.y1, a.y0 - b.y1)
    if dx > 0 and dy > 0: return math.hypot(dx, dy)
    return max(dx, dy)


pairs = []
for i in range(len(boxes)):
    for j in range(i + 1, len(boxes)):
        pairs.append((sep(boxes[i][2], boxes[j][2]), boxes[i][0], boxes[j][0]))
overlaps = [p for p in pairs if p[0] < 0]
dpi = fig.dpi
closest = sorted(pairs)[:6]
title_gap = {}
for k in range(1, N):
    tb_ = [b for n, t, b in boxes if n == f'title ({TAGS[k]})'][0]
    above = [(sep(tb_, b), n) for n, t, b in boxes if f'({TAGS[k-1]})' in n and not n.startswith('title')]
    title_gap[k] = min(above)
W, H = fig.bbox.width, fig.bbox.height
outside = [(n, b.x0, b.y0, b.x1, b.y1) for n, t, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5]
sizes = sorted({round(t.get_fontsize(), 2) for n, t, b in boxes})
aspects = []
for k, a in enumerate(axs):
    bb = a.get_window_extent(R); x0, x1 = a.get_xlim(); y0, y1 = a.get_ylim()
    aspects.append(abs(((bb.height / abs(y1 - y0)) / (bb.width / abs(x1 - x0)))))
# rendered image: data per pixel in y over data per pixel in x (ParaView window against the 9 D by 1.1 D extent)
img_asp = {n: ((EXT[3] - EXT[2]) / im.shape[0]) / ((EXT[1] - EXT[0]) / im.shape[1]) for n, im in IMGS.items()}

fig.savefig(P3 + 'fig_t3rho.pdf', bbox_inches='tight', pad_inches=PAD, dpi=300)
fig.savefig(P3 + 'fig_t3rho_preview.png', bbox_inches='tight', pad_inches=PAD, dpi=110)
mb = PdfReader(P3 + 'fig_t3rho.pdf').pages[0].mediabox
pw, ph = float(mb.width) / 72, float(mb.height) / 72


# ---------------- numbers from the VTK input arrays ----------------
def vtk(name):
    with open(PV + name + '.vtk') as f:
        for _ in range(5): f.readline()
        n = int(f.readline().split()[1])
        pts = np.loadtxt(f, max_rows=n)
        for _ in range(3): f.readline()
        sp = np.loadtxt(f, max_rows=n)
        for _ in range(2): f.readline()
        al = np.loadtxt(f, max_rows=n)
    return pts, sp, al


rows = []
for n, t, case, folder in PANELS:
    pts, sp, al = vtk(n)
    g = al < 0.5
    rows.append((n, case, folder, len(sp) // 2, sp.max(), sp[al > 0.99].max(), pts[g, 0].min(), pts[g, 0].max(), pts[:, 0].min(), pts[:, 0].max(),
                 int((sp > 2).sum()) // 2, 100 * (sp > 2).mean()))

L = ['# A11 check: fig_t3rho', '',
    f'Script: `fig_t3rho.py` (images from `pv/render_t3rho.py`, VTK from `tovtk_t3rho.py`). Texts checked: {len(boxes)} visible non-empty text artists '
    f'({N} titles, y labels, the x label, all tick labels, colour-bar label and tick labels).', '',
    'Snapshots: T3_rhoST folders `1`, `2.4`, `5`, `8.4`; T3_prod folder `8.272418` (the folder nearest 8.27 s in `runs_p2/T3_prod`; its neighbours are 7.52038 and 9.024456; there is no folder `8.272406`). Panel (e) title uses 8.27 s.', '',
    '## 1. Text overlap test (after fig.canvas.draw())', '',
    f'- Intersecting pairs: **{len(overlaps)}**' + ('' if not overlaps else ' ' + str(overlaps)),
    f'- Texts extending outside the figure canvas: **{len(outside)}**' + ('' if not outside else ' ' + str(outside)),
    '- Closest non-intersecting pairs (separation in pt): ' + '; '.join(f'{p[1]} / {p[2]} = {p[0] / dpi * 72:.1f}' for p in closest),
    '- Title against the text of the strip above (separation in pt): ' + '; '.join(
        f'({TAGS[k]}) title vs {title_gap[k][1]} = {title_gap[k][0] / dpi * 72:.1f}' for k in range(1, N)),
    '', '## 2. Font sizes', '',
    f'- Smallest font size used: **{min(sizes):.1f} pt**; all sizes present: {", ".join(f"{s:.1f}" for s in sizes)} pt',
    '', '## 3. Data aspect of each strip (y data per px over x data per px, 1 = equal)', '']
L += [f'- strip ({TAGS[k]}): {aspects[k]:.4f} ({(aspects[k] - 1) * 100:+.2f} %)' for k in range(N)]
L += [f'- rendered ParaView images, {IMGS["T3R_1"].shape[1]} x {IMGS["T3R_1"].shape[0]} px against the 9 D x 1.1 D extent: ' +
      ', '.join(f'{n} {(v - 1) * 100:+.3f} %' for n, v in img_asp.items()),
      '', '## 4. Final figure size from the PDF (MediaBox)', '',
      f'- **{pw:.3f} in wide by {ph:.3f} in tall** ({float(mb.width):.1f} x {float(mb.height):.1f} pt)',
      f'- Axes width {wax:.3f} in, strip height {wax * ASP:.3f} in, gap between strips {GAP:.3f} in; colour bar {CB_W:.2f} in wide, spans from the bottom of (e) to the top of (a)',
      '', '## 5. Numbers from the VTK input arrays (cell-centre values; the half domain is mirrored, so mirrored copies do not change any max or min)', '',
      '| panel | case | folder | cells (half) | max speed, all cells (m/s) | max speed, alpha > 0.99 (m/s) | min x of cells with alpha < 0.5 | max x of cells with alpha < 0.5 | x range of grid | cells with speed > 2 m/s |',
      '|---|---|---|---|---|---|---|---|---|---|']
for (n, case, folder, nc, mx, ml, xmin, xmax, gx0, gx1, n2, p2), tag in zip(rows, TAGS):
    L.append(f'| ({tag}) {n} | {case} | {folder} | {nc} | {mx:.3f} | {ml:.3f} | {xmin:.4f} | {xmax:.4f} | {gx0:.4f} to {gx1:.4f} | {n2} ({p2:.1f} %) |')
open(P3 + 'fig_t3rho_autocheck.md', 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L))
