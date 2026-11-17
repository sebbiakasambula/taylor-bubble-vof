import math, sys
sys.stdout.reconfigure(encoding='utf-8')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import Normalize; import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from pypdf import PdfReader

from config import PV as PVDIR, WORK
PV = str(PVDIR / 'strips') + '/'   # strips (images) from paraview/render_t1_flat.py
P3 = str(WORK) + '/'
EXT = (-5.325, 0.725, -0.55, 0.55)
PANELS = [('T1_OF_str', '(a) T1, OpenFOAM, stretched, $t$ = 19.7 s'),
          ('T1_OF_sq',  '(b) T1, OpenFOAM, square, $t$ = 19.7 s'),
          ('T1_F_str',  '(c) T1, Fluent, stretched, $t$ = 19.7 s'),
          ('T1_F_sq',   '(d) T1, Fluent, square, $t$ = 19.2 s')]
TARGET_W, TARGET_H, PAD = 6.3, 5.15, 0.02   # inches, final PDF size and tight-bbox pad
WCORR = 0.03                                # PDF text extents are a little narrower than Agg ones, so build this much wider
CB_GAP, CB_W = 0.12, 0.14                   # inches, gap axes to colour bar and its width

plt.rcParams.update({'font.size': 8.5})
IMGS = {n: mpimg.imread(PV + n + '.png') for n, _ in PANELS}
ASP = (EXT[3] - EXT[2]) / (EXT[1] - EXT[0])   # strip height over width, so the data aspect is 1


def build(wax, gap, wc, hc, ox, oy):
    """Strips placed in absolute inches: ox = left edge of the axes, oy = bottom edge of strip (d)."""
    fig = plt.figure(figsize=(wc, hc))
    hax = wax * ASP
    axs = []
    for i, (n, t) in enumerate(PANELS):
        yb = oy + (3 - i) * (hax + gap)
        a = fig.add_axes([ox / wc, yb / hc, wax / wc, hax / hc], sharex=axs[0] if axs else None)
        a.imshow(IMGS[n], extent=EXT, aspect='equal')
        a.set_title(t, fontsize=8, loc='left')
        a.set_ylabel('$r/D$'); a.set_yticks([-0.5, 0, 0.5]); a.set_xticks(range(-5, 1))
        if i < 3: a.tick_params(labelbottom=False)
        axs.append(a)
    axs[-1].set_xlabel('Distance from the nose, $x/D$')
    c = fig.add_axes([(ox + wax + CB_GAP) / wc, oy / hc, CB_W / wc, (3 * (hax + gap) + hax) / hc])
    cb = fig.colorbar(cm.ScalarMappable(Normalize(0, 2), 'viridis'), cax=c, extend='max')
    cb.set_ticks([0, 0.5, 1, 1.5, 2]); c.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    cb.set_label('Speed (m/s)')
    return fig, axs, c


def tight(fig):
    fig.canvas.draw()
    return fig.get_tightbbox(fig.canvas.get_renderer())


# solve for the axes width and the strip gap so the cropped figure is TARGET_W by TARGET_H
wax, gap, ox, oy = 4.85, 0.34, 1.5, 1.5
for it in range(8):
    fig, axs, c = build(wax, gap, 9, 9, ox, oy)
    tb = tight(fig); plt.close(fig)
    cw, chh = tb.width, tb.height
    dw = (TARGET_W + WCORR - 2 * PAD) - cw
    dh = (TARGET_H - 2 * PAD) - chh - 4 * ASP * dw
    ox, oy = ox - tb.x0 + PAD, oy - tb.y0 + PAD
    wax += dw; gap += dh / 3
    print(f'iter {it}: content {cw:.4f} x {chh:.4f}  -> wax {wax:.4f} gap {gap:.4f}')
    if abs(dw) < 1e-4 and abs(dh) < 1e-4: break

fig, axs, cax = build(wax, gap, TARGET_W + WCORR, TARGET_H, ox, oy)
fig.canvas.draw(); R = fig.canvas.get_renderer()

# ---------------- checks ----------------
texts = []   # (label, Text artist)
for k, a in enumerate(axs):
    tag = 'abcd'[k]
    texts.append((f'title ({tag})', a._left_title))
    texts.append((f'ylabel ({tag})', a.yaxis.label))
    texts.append((f'xlabel ({tag})', a.xaxis.label))
    for j, t in enumerate(a.xaxis.get_majorticklabels() + a.xaxis.get_minorticklabels()): texts.append((f'xtick ({tag}) {t.get_text()}', t))
    for j, t in enumerate(a.yaxis.get_majorticklabels() + a.yaxis.get_minorticklabels()): texts.append((f'ytick ({tag}) {t.get_text()}', t))
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
# each title against the texts of the strip above (a title only has strips above for b, c, d)
title_gap = {}
for k in (1, 2, 3):
    tb_ = [b for n, t, b in boxes if n == f'title ({"abcd"[k]})'][0]
    above = [(sep(tb_, b), n) for n, t, b in boxes if f'({"abcd"[k-1]})' in n and not n.startswith('title')]
    title_gap[k] = min(above)
W, H = fig.bbox.width, fig.bbox.height
outside = [(n, b.x0, b.y0, b.x1, b.y1) for n, t, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5]
sizes = sorted({round(t.get_fontsize(), 2) for n, t, b in boxes})
size_by = {}
for n, t, b in boxes: size_by.setdefault(round(t.get_fontsize(), 2), []).append(n.split(' ')[0] + ' ' + n.split(' ')[1] if ' ' in n else n)
aspects = []
for k, a in enumerate(axs):
    bb = a.get_window_extent(R); x0, x1 = a.get_xlim(); y0, y1 = a.get_ylim()
    aspects.append(abs(((bb.height / abs(y1 - y0)) / (bb.width / abs(x1 - x0)))))

fig.savefig(P3 + 'fig_t1_fields.pdf', bbox_inches='tight', pad_inches=PAD, dpi=300)
fig.savefig(P3 + 'fig_t1_fields_preview.png', bbox_inches='tight', pad_inches=PAD, dpi=110)
mb = PdfReader(P3 + 'fig_t1_fields.pdf').pages[0].mediabox
pw, ph = float(mb.width) / 72, float(mb.height) / 72

L = ['# B2 check: fig_t1_fields', '',
     f'Script: `fig_t1_fields.py`. Texts checked: {len(boxes)} visible non-empty text artists '
     '(4 titles, y labels, the x label, all tick labels, colour-bar label and tick labels).', '',
     '## 1. Text overlap test (after fig.canvas.draw())', '',
     f'- Intersecting pairs: **{len(overlaps)}**',
     f'- Texts extending outside the figure canvas: **{len(outside)}**' + ('' if not outside else ' ' + str(outside)),
     f'- Closest non-intersecting pairs (separation in pt): ' + '; '.join(f'{p[1]} / {p[2]} = {p[0] / dpi * 72:.1f}' for p in closest),
     '- Title against the text of the strip above (separation in pt): ' + '; '.join(
         f'({"abcd"[k]}) title vs {title_gap[k][1]} = {title_gap[k][0] / dpi * 72:.1f}' for k in (1, 2, 3)),
     '', '## 2. Font sizes', '',
     f'- Smallest font size used: **{min(sizes):.1f} pt**; all sizes present: {", ".join(f"{s:.1f}" for s in sizes)} pt',
     '', '## 3. Data aspect of each strip (y data per px over x data per px, 1 = equal)', '']
L += [f'- strip ({"abcd"[k]}): {aspects[k]:.4f} ({(aspects[k] - 1) * 100:+.2f} %)' for k in range(4)]
L += ['', '## 4. Final figure size from the PDF (MediaBox)', '',
      f'- **{pw:.3f} in wide by {ph:.3f} in tall** ({float(mb.width):.1f} x {float(mb.height):.1f} pt)',
      f'- Axes width {wax:.3f} in, strip height {wax * ASP:.3f} in, gap between strips {gap:.3f} in; colour bar {CB_W:.2f} in wide, spans from the bottom of (d) to the top of (a)']
open(P3 + 'fig_t1_fields_autocheck.md', 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L))
