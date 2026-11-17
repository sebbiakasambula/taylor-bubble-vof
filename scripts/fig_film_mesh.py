"""B10: mesh across the liquid film at T1 (t = 19.7 s), stretched against square, one figure of two panels.
Field rendering (cell colours, cell edges, alpha = 0.5 line) comes from ParaView: paraview/render_film_mesh.py, fed by
prep_film_mesh.py (front wedge face of each OpenFOAM mesh, read from polyMesh, read-only).
This script adds the axes, wall line, correlation line, legend, measures the mesh and film numbers and runs the checks."""
import math, sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter
from pypdf import PdfReader

from config import WORK
S = str(WORK) + '/'
PV, P3 = S + 'pv/', S
XW, RW = (-2.56, -2.44), (0.42, 0.505)
CASES = [('V4_medium', '(a) Stretched, 96 radial cells'), ('VA_square', '(b) Square, 48 radial cells')]
A = math.radians(2.5); COSA = math.cos(A); Y_WALL = 0.5 * COSA          # wall in the wedge mid-plane
DELTA = 0.0615                                                          # film thickness from the correlation, in D
PX_PER_UNIT = 20000
WPX, HPX = round((XW[1] - XW[0]) * PX_PER_UNIT), round((RW[1] - RW[0]) * PX_PER_UNIT)   # render size in px
TARGET_W, PAD = 6.3, 0.02
GAS, LIQ = (0.90, 0.95, 1.00), (0.76, 0.76, 0.76)
IFACE, CORR = '#D55E00', '#0072B2'

# ---------------------------------------------------------------- measurements (mid-plane, from the mesh and alpha field)
M = {}
for case, _ in CASES:
    d = np.load(PV + f'film_mesh_{case}.npz')
    ux, ur, a2 = d['ux'], d['ur'], d['a2']                  # column centres (x - nose), row centres (mid-plane y), alpha[col, row]
    ylo = np.zeros(len(ur)); yhi = np.zeros(len(ur)); xlo = np.zeros(len(ux)); xhi = np.zeros(len(ux))
    for ir in range(len(ur)):
        m = d['ir'] == ir; ylo[ir] = d['y_lo'][m][0]; yhi[ir] = d['y_hi'][m][0]
        assert np.allclose(d['y_lo'][m], ylo[ir]) and np.allclose(d['y_hi'][m], yhi[ir])
    for ic in range(len(ux)):
        m = d['ix'] == ic; xlo[ic] = d['x_lo'][m][0]; xhi[ic] = d['x_hi'][m][0]
    assert abs(yhi[-1] - Y_WALL) < 2e-6, (yhi[-1], Y_WALL)   # top row ends on the mid-plane wall
    # (1) cells whose centre lies within DELTA of the wall, and size of the wall cell, axial size
    n_in = int(np.sum(Y_WALL - ur <= DELTA))
    r_mid = 0.5 * (ylo + yhi) / COSA                         # same count with the front-face radius (wall at r = 0.5), arithmetic mid-radius
    n_in_r = int(np.sum(r_mid > 0.5 - DELTA))
    dy_wall = yhi[-1] - ylo[-1]; dx = xhi - xlo
    # (2) film thickness per column: first inward crossing of alpha = 0.5 from the wall, linear between cell centres
    in_win = np.where((ux >= XW[0]) & (ux <= XW[1]))[0]
    th, nc, nfrac, ncross = [], [], [], []
    for ic in in_win:
        a = a2[ic]
        assert a[-1] >= 0.5, 'wall cell is not liquid'
        k = len(a) - 1
        while k > 0 and a[k - 1] >= 0.5: k -= 1              # k = innermost row still liquid (alpha >= 0.5) in the unbroken run from the wall
        assert k > 0, 'no crossing inside the kept rows'
        a1, a0 = a[k], a[k - 1]; y1, y0 = ur[k], ur[k - 1]   # a1 >= 0.5 > a0, y1 > y0
        yc = y1 + (0.5 - a1) / (a0 - a1) * (y0 - y1)
        ncross.append(int(np.sum(np.diff(np.sign(a - 0.5)) != 0)))
        th.append(Y_WALL - yc)
        nc.append(int(np.sum(ur >= yc)))                     # cell centres between the wall and the crossing
        edges = np.concatenate([yhi[::-1], [ylo[0]]])        # y edges from the wall inward
        j = int(np.searchsorted(-edges, -yc) - 1)            # edges[j] >= yc > edges[j+1]
        nfrac.append(j + (edges[j] - yc) / (edges[j] - edges[j + 1]))
    M[case] = dict(ux=ux, ur=ur, a2=a2, ylo=ylo, yhi=yhi, xlo=xlo, xhi=xhi, n_in=n_in, n_in_r=n_in_r, dy_wall=dy_wall,
                   dx=dx.mean(), dx_spread=dx.max() - dx.min(), th=np.array(th), nc=np.array(nc), nfrac=np.array(nfrac),
                   ncross=ncross, ncols=len(in_win), in_win=in_win, nose=float(d['nose']))

# ---------------------------------------------------------------- render alignment: pixel rows of cell edges and of the orange line
IMGS = {c: mpimg.imread(PV + f'film_{c}.png')[:, :, :3] for c, _ in CASES}
align = {}
for c, _ in CASES:
    img = IMGS[c]; H, W = img.shape[:2]
    assert (W, H) == (round((XW[1] - XW[0]) * PX_PER_UNIT), round((RW[1] - RW[0]) * PX_PER_UNIT))
    m = M[c]
    exp_rows = np.unique(np.round((RW[1] - np.concatenate([m['ylo'], m['yhi']]) / COSA) * PX_PER_UNIT, 3))   # expected edge rows (px from top)
    exp_rows = exp_rows[(exp_rows > 3) & (exp_rows < H - 3)]
    # a column strip in the middle of a cell, away from vertical edges and from the orange line
    cols = np.arange(W)
    xedge_px = np.round((np.concatenate([m['xlo'], m['xhi']]) - XW[0]) * PX_PER_UNIT)
    ok = np.ones(W, bool)
    for xe in xedge_px: ok &= np.abs(cols - xe) > 12
    ior = (np.abs(img[:, :, 0] - 0.835) < 0.06) & (np.abs(img[:, :, 1] - 0.369) < 0.06) & (np.abs(img[:, :, 2] - 0.0) < 0.08)
    orows = np.where(ior.mean(axis=1) > 0.2)[0]               # rows covered by the orange line (it spans the whole panel width)
    prof = img[:, ok, :].mean(axis=(1, 2))                    # mean intensity per row over the strip columns
    dark = prof < 0.5                                         # edge lines are about 0.15 grey, fills are 0.76 / 0.95
    rows = np.where(dark)[0]
    groups = np.split(rows, np.where(np.diff(rows) > 1)[0] + 1)
    cen = np.array([g.mean() + 0.5 for g in groups if len(g) <= 8])   # px-centre of each dark row group (the 12 px orange line is longer)
    far = np.array([len(orows) == 0 or np.min(np.abs(orows - e)) > 14 for e in exp_rows])   # skip edges hidden under the orange line
    off = [np.min(np.abs(cen - e)) for e in exp_rows[far]]
    offo = []
    for ic in m['in_win']:
        xpx = int(round((m['ux'][ic] - XW[0]) * PX_PER_UNIT))
        a = m['a2'][ic]; k = len(a) - 1
        while k > 0 and a[k - 1] >= 0.5: k -= 1
        yc = m['ur'][k] + (0.5 - a[k]) / (a[k - 1] - a[k]) * (m['ur'][k - 1] - m['ur'][k])
        exp_px = (RW[1] - yc / COSA) * PX_PER_UNIT
        r_ = np.where(ior[:, xpx])[0]
        offo.append(r_.mean() + 0.5 - exp_px if len(r_) else np.nan)
    align[c] = dict(n_exp=int(far.sum()), n_hidden=int((~far).sum()), n_found=len(cen), edge_max=max(off), edge_mean=float(np.mean(off)), iface_max=float(np.nanmax(np.abs(offo))),
                    iface_n=len(offo), iface_nan=int(np.isnan(offo).sum()))

# ---------------------------------------------------------------- figure
plt.rcParams.update({'font.size': 8.5})
ASP = (RW[1] - RW[0]) / (XW[1] - XW[0])
GAP, LEG_DROP = 0.16, 0.42        # inches: gap between panels; legend top below the panel bottom (clears tick labels and x label)


def build(wax, wc, hc, ox, oy):
    fig = plt.figure(figsize=(wc, hc))
    hax = wax * ASP
    axs = []
    for i, (c, title) in enumerate(CASES):
        a = fig.add_axes([(ox + i * (wax + GAP)) / wc, oy / hc, wax / wc, hax / hc])
        a.imshow(IMGS[c], extent=(XW[0], XW[1], RW[0], RW[1]), aspect='equal', interpolation='antialiased')
        a.set_xlim(XW); a.set_ylim(RW)
        a.axhline(0.5, color='k', lw=1.1, solid_capstyle='butt')                                 # wall
        a.axhline(0.5 - DELTA, color=CORR, lw=1.0, ls=(0, (4, 2)))                             # film thickness from the correlation
        a.set_title(title, fontsize=8, loc='left')
        a.set_xlabel(r'$(x - x_\mathrm{nose})/D$')
        a.set_xticks([-2.55, -2.50, -2.45]); a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        a.set_yticks([0.42, 0.44, 0.46, 0.48, 0.50]); a.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        if i == 0: a.set_ylabel(r'$r/D$')
        else: a.tick_params(labelleft=False)
        axs.append(a)
    hs = [Patch(facecolor=LIQ, edgecolor='0.15', lw=0.5, label=r'liquid ($\alpha$ = 1)'),
          Patch(facecolor=GAS, edgecolor='0.15', lw=0.5, label=r'gas ($\alpha$ = 0)'),
          Line2D([], [], color=IFACE, lw=1.0, label=r'$\alpha$ = 0.5'),
          Line2D([], [], color='k', lw=1.1, label='wall'),
          Line2D([], [], color=CORR, lw=1.0, ls=(0, (4, 2)), label=r'$\delta$ = 0.0615 $D$ (correlation)')]
    leg = fig.legend(handles=hs, loc='upper center', bbox_to_anchor=((ox + (2 * wax + GAP) / 2) / wc, (oy - LEG_DROP) / hc),
                     ncol=5, frameon=False, fontsize=8, handlelength=1.5, columnspacing=1.0, handletextpad=0.5, borderpad=0)
    return fig, axs, leg


def tight(fig):
    fig.canvas.draw()
    return fig.get_tightbbox(fig.canvas.get_renderer())


WCORR = 0.0
wax, ox, oy = 2.8, 1.0, 1.0
for outer in range(6):
    for it in range(8):
        fig, axs, leg = build(wax, 9, 9, ox, oy)
        tb = tight(fig); plt.close(fig)
        dw = (TARGET_W + WCORR - 2 * PAD) - tb.width
        ox, oy = ox - tb.x0 + PAD, oy - tb.y0 + PAD
        wax += dw / 2
        if abs(dw) < 1e-4: break
    hc = tb.height + 2 * PAD
    fig, axs, leg = build(wax, TARGET_W + WCORR, hc, ox, oy)
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    fig.savefig(P3 + 'fig_film_mesh.pdf', bbox_inches='tight', pad_inches=PAD, dpi=864)   # dpi sets the embedded image resolution: 864 = 2400 px over 2.78 in, no downsampling
    mb = PdfReader(P3 + 'fig_film_mesh.pdf').pages[0].mediabox
    pw, ph = float(mb.width) / 72, float(mb.height) / 72
    print(f'outer {outer}: wax {wax:.4f} WCORR {WCORR:.4f} -> PDF {pw:.4f} x {ph:.4f} in')
    if abs(pw - TARGET_W) < 0.002: break
    WCORR += TARGET_W - pw
    plt.close(fig)
fig.savefig(P3 + 'fig_film_mesh_preview.png', bbox_inches='tight', pad_inches=PAD, dpi=110)
fig.savefig(S + 'pv/film_mesh_300.png', bbox_inches='tight', pad_inches=PAD, dpi=300)   # for my own look at print scale (not a deliverable)

# ---------------------------------------------------------------- checks
fig.canvas.draw(); R = fig.canvas.get_renderer()   # savefig at other dpi moved the cached text positions; redraw at the figure dpi
texts = []
for k, a in enumerate(axs):
    tag = 'ab'[k]
    texts += [(f'title ({tag})', a._left_title), (f'ylabel ({tag})', a.yaxis.label), (f'xlabel ({tag})', a.xaxis.label)]
    for t in a.xaxis.get_majorticklabels() + a.xaxis.get_minorticklabels(): texts.append((f'xtick ({tag}) {t.get_text()}', t))
    for t in a.yaxis.get_majorticklabels() + a.yaxis.get_minorticklabels(): texts.append((f'ytick ({tag}) {t.get_text()}', t))
    texts += [(f'xoffset ({tag})', a.xaxis.offsetText), (f'yoffset ({tag})', a.yaxis.offsetText)]
for t in leg.get_texts(): texts.append((f'legend "{t.get_text()}"', t))
texts = [(n, t) for n, t in texts if t.get_visible() and t.get_text().strip()]
boxes = [(n, t, t.get_window_extent(R)) for n, t in texts]


def sep(a, b):
    dx = max(b.x0 - a.x1, a.x0 - b.x1); dy = max(b.y0 - a.y1, a.y0 - b.y1)
    if dx > 0 and dy > 0: return math.hypot(dx, dy)
    return max(dx, dy)


pairs = [(sep(boxes[i][2], boxes[j][2]), boxes[i][0], boxes[j][0]) for i in range(len(boxes)) for j in range(i + 1, len(boxes))]
overlaps = [p for p in pairs if p[0] < 0]
closest = sorted(pairs)[:6]
W_, H_ = fig.bbox.width, fig.bbox.height
outside = [(n, b.x0, b.y0, b.x1, b.y1) for n, t, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W_ + 0.5 or b.y1 > H_ + 0.5]
sizes = sorted({round(t.get_fontsize(), 2) for n, t, b in boxes})
# legend box against the axes boxes (the legend must not sit on an axes)
lb = leg.get_window_extent(R); axb = [a.get_window_extent(R) for a in axs]
leg_vs_ax = [min(0, sep(lb, b)) for b in axb]
leg_gap_ax = [sep(lb, b) / fig.dpi * 72 for b in axb]
aspects = []
for a in axs:
    bb = a.get_window_extent(R); x0, x1 = a.get_xlim(); y0, y1 = a.get_ylim()
    aspects.append(abs((bb.height / abs(y1 - y0)) / (bb.width / abs(x1 - x0))))
dpi = fig.dpi
import os, json
_tb = fig.get_tightbbox(R)
json.dump(dict(a_left_in=ox - (_tb.x0 - PAD), wax_in=wax, hax_in=wax * ASP, total_h_in=_tb.height + 2 * PAD, top_axes_from_top_in=(_tb.y1 + PAD) - (oy + wax * ASP)),
          open(S + 'pv/film_mesh_geom.json', 'w'))
pdf_kb = os.path.getsize(P3 + 'fig_film_mesh.pdf') / 1024
# ---------------------------------------------------------------- write the check file (visual section is added by hand after looking at the preview)
f4 = lambda v: f'{v:.4f}'
f5 = lambda v: f'{v:.5f}'
L = ['# B10 check: fig_film_mesh', '',
     'Script: `p3/fig_film_mesh.py` (axes and assembly, measurements, checks); `pv/render_film_mesh.py` (ParaView, cell colours, edges, alpha = 0.5 line); '
     '`pv/prep_film_mesh.py` (builds the VTK files from polyMesh). Case folders were only read.', '',
     '## 0. How the mesh was read',
     '',
     '- The ParaView OpenFOAM reader needs `case.foam` inside the case folder (tested: with `case.foam` in a new folder elsewhere it fails with "Cannot open directory .../constant/"), so it was not used.',
     '- Instead the cell polygons are the front-wedge-face quads of each mesh, built in Python from `constant/polyMesh` (points, faces, owner: faces of patch `front`, owner cell of each face) with `alpha.water` of time 19.700352 as cell data, written as legacy VTK (UNSTRUCTURED_GRID, quads) to `pv/film_mesh_<case>.vtk`; ParaView draws them as Surface With Edges.',
     '- Plot coordinates: x - x_nose (nose 11.3199 for V4_medium, 11.5606 for VA_square) and r = sqrt(y^2 + z^2) of the front-face vertices, so the wall is at r = 0.5 in the figure.',
     '- The alpha = 0.5 line is the 0.5 contour of alpha at the cell centres (cell-centre coordinates from the solver file `0/C`, checked against the analytic wedge-cell centroid: max difference 5e-9), i.e. the same interpolation as the film-thickness measurement below.',
     '- Numbers below are in the wedge mid-plane (cell centres at z = 0, wall at y = 0.5 cos 2.5 deg = ' + f'{Y_WALL:.6f}' + '), as requested.', '']
L += ['## 1. Mesh numbers (from the mesh)', '',
      '| quantity | V4_medium (stretched) | VA_square (square) |', '|---|---|---|']
a, b = M['V4_medium'], M['VA_square']
L.append(f'| radial cells | 96 | 48 |')
L.append(f'| radial cells whose centre lies within {DELTA} of the wall (mid-plane, wall {Y_WALL:.6f}) | **{a["n_in"]}** | **{b["n_in"]}** |')
L.append(f'| same count with front-face radius (wall 0.5, mid-radius of each cell) | {a["n_in_r"]} | {b["n_in_r"]} |')
L.append(f'| radial size of the wall cell, mid-plane dy | {f5(a["dy_wall"])} | {f5(b["dy_wall"])} |')
L.append(f'| radial size of the wall cell, front face dr = dy / cos 2.5 deg | {f5(a["dy_wall"] / COSA)} | {f5(b["dy_wall"] / COSA)} |')
L.append(f'| axial cell size dx | {f5(a["dx"])} | {f5(b["dx"])} |')
L.append(f'| wall cell aspect dx / dy | {a["dx"] / a["dy_wall"]:.2f} | {b["dx"] / b["dy_wall"]:.2f} |')
L.append(f'| axial size spread in the window (max - min) | {a["dx_spread"]:.1e} | {b["dx_spread"]:.1e} |')
L += ['', f'Wall-cell radial size and axial size are exact from the polyMesh vertices (not rounded to the figure). dx = 20/480 and 20/1920 resp.: '
      f'{20 / 480:.5f}, {20 / 1920:.5f}.', '']
L += [f'## 2. Film thickness from the alpha field in the panel window (x - x_nose from {XW[0]} to {XW[1]})', '',
      'Method: for each axial column whose cell centre lies in the window, take alpha.water at the cell centres along r, '
      'walk in from the wall (alpha >= 0.5) to the first cell centre with alpha < 0.5, interpolate linearly between those two cell centres to alpha = 0.5, '
      f'thickness = {Y_WALL:.6f} - y(crossing). Cells across the film: (i) cell centres between the wall and the crossing; '
      '(ii) the wall-to-interface distance in cell units (full cells plus the fraction of the cell holding the crossing).', '',
      '| quantity | V4_medium (stretched) | VA_square (square) |', '|---|---|---|',
      f'| columns in the window | {a["ncols"]} | {b["ncols"]} |',
      f'| crossings of alpha = 0.5 per column along the whole kept profile (all columns) | {sorted(set(a["ncross"]))} | {sorted(set(b["ncross"]))} |',
      f'| mean film thickness, mid-plane | **{f4(a["th"].mean())}** | **{f4(b["th"].mean())}** |',
      f'| min / max over columns | {f4(a["th"].min())} / {f4(a["th"].max())} | {f4(b["th"].min())} / {f4(b["th"].max())} |',
      f'| mean thickness / correlation 0.0615 | {a["th"].mean() / DELTA:.3f} | {b["th"].mean() / DELTA:.3f} |',
      f'| cells across the film (i): cell centres between wall and crossing, mean (min-max) | **{a["nc"].mean():.2f}** ({a["nc"].min()}-{a["nc"].max()}) | **{b["nc"].mean():.2f}** ({b["nc"].min()}-{b["nc"].max()}) |',
      f'| cells across the film (ii): thickness in cell units, mean (min-max) | **{a["nfrac"].mean():.2f}** ({a["nfrac"].min():.2f}-{a["nfrac"].max():.2f}) | **{b["nfrac"].mean():.2f}** ({b["nfrac"].min():.2f}-{b["nfrac"].max():.2f}) |',
      f'| thickness / wall-cell radial size | {a["th"].mean() / a["dy_wall"]:.1f} | {b["th"].mean() / b["dy_wall"]:.1f} |', '']
L += ['## 3. Render alignment (ParaView PNG against the mesh numbers)', '',
      f'Cell-edge rows found in the render versus the expected pixel row of every radial cell edge ({PX_PER_UNIT} px per unit of D), and the orange line against the alpha = 0.5 crossing above:', '']
for c, _ in CASES:
    q = align[c]
    L.append(f'- {c}: expected edge rows {q["n_exp"]} (plus {q["n_hidden"]} under the orange line, skipped), dark row groups found {q["n_found"]}; largest offset of an expected edge from the nearest dark row {q["edge_max"]:.2f} px '
             f'(mean {q["edge_mean"]:.2f} px); orange line centre vs expected crossing, largest offset over {q["iface_n"]} window columns {q["iface_max"]:.2f} px '
             f'({q["iface_nan"]} columns without orange pixels). 1 px = {1 / PX_PER_UNIT:.5f} D.')
L += ['', '## 4. Text overlap test (after fig.canvas.draw())', '',
      f'- Texts checked: {len(boxes)} (2 titles, 2 x labels, y label, all tick labels, 5 legend labels).',
      f'- Intersecting pairs: **{len(overlaps)}**' + ('' if not overlaps else ' ' + str(overlaps)),
      f'- Texts extending outside the figure canvas: **{len(outside)}**' + ('' if not outside else ' ' + str(outside)),
      '- Closest non-intersecting pairs (separation in pt): ' + '; '.join(f'{p[1]} / {p[2]} = {p[0] / dpi * 72:.1f}' for p in closest),
      f'- Legend box against axes (a), (b): separation {leg_gap_ax[0]:.1f} pt, {leg_gap_ax[1]:.1f} pt ({"no overlap" if min(leg_vs_ax) == 0 else "OVERLAP"}).', '',
      '## 5. Font sizes', '',
      f'- Smallest font size used: **{min(sizes):.1f} pt**; all sizes present: {", ".join(f"{s:.1f}" for s in sizes)} pt', '',
      '## 6. Data aspect (y data per px over x data per px, 1 = equal)', '',
      f'- panel (a): {aspects[0]:.4f} ({(aspects[0] - 1) * 100:+.2f} %); panel (b): {aspects[1]:.4f} ({(aspects[1] - 1) * 100:+.2f} %)', '',
      '## 7. Final figure size from the PDF', '',
      f'- **{pw:.3f} in wide by {ph:.3f} in tall** ({float(mb.width):.1f} x {float(mb.height):.1f} pt); file size {pdf_kb:.0f} kB',
      f'- Axes width {wax:.3f} in each, panel height {wax * ASP:.3f} in, gap between panels {GAP:.2f} in; the panel images are {WPX} x {HPX} px each '
      f'({WPX / wax:.0f} px per inch in the PDF); in the PDF, one unit of r/D is {wax / (XW[1] - XW[0]):.2f} in, so the wall cell of (a) is {a["dy_wall"] / COSA * wax / (XW[1] - XW[0]) * 72:.2f} pt tall and a cell of (b) is {b["dy_wall"] / COSA * wax / (XW[1] - XW[0]) * 72:.2f} pt tall.']
open(P3 + 'checks/B10_fig.md', 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L))
