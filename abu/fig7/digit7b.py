"""Read the markers and the lines of Abubakar and Matar (arXiv 2109.09091) Fig. 7b, film thickness against Eo.

Same method as ../digit.py for Fig. 6a. The axes are found from the spines and tick marks, the markers by colour after
an opening that removes the thin lines. The lines are Brown's relation, their eq. (3.3), with Ub from eq. (3.1), so they
are read at Eo 80, 120 and 200 and checked against the same relation computed here. That gives the reading error.
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as nd
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
from tb_references import fr_viana  # noqa: E402

im = np.array(Image.open(HERE / 'fig7_img1_x556.png').convert('RGB')).astype(int)
H, W, _ = im.shape
dark = im.sum(2) < 150
colsum = dark.sum(0)
rowsum = dark.sum(1)
xl = int(np.argmax(colsum[:W // 2]))
xr = W // 2 + int(np.argmax(colsum[W // 2:]))
yb = H // 2 + int(np.argmax(rowsum[H // 2:]))
yt = int(np.argmax(rowsum[:H // 2]))
print("spines x", xl, xr, "y", yt, yb)
xt = [x for x in range(xl, xr + 1) if dark[yb + 1:yb + 5, x].all()]
ytk = [y for y in range(yt, yb + 1) if dark[y, xl - 4:xl].all()]


def groups(v):
    g = [[v[0]]]
    for a in v[1:]:
        (g[-1].append(a) if a - g[-1][-1] <= 1 else g.append([a]))
    return [float(np.mean(x)) for x in g]


xt = groups(xt)
ytk = groups(ytk)
print("x tick px", np.round(xt, 2))
print("y tick px", np.round(ytk, 2))
xv = np.arange(0, 301, 50)[:len(xt)]
yv = np.round(np.arange(0.15, 0.0749, -0.01), 2)[:len(ytk)]
print("tick values x", xv, "y", yv)
ax = np.polyfit(xt, xv, 1)
ay = np.polyfit(ytk, yv, 1)
print("x fit resid", np.round(np.polyval(ax, xt) - xv, 4), "y fit resid", np.round(np.polyval(ay, ytk) - yv, 6))

cols = {40: (31, 119, 180), 60: (255, 127, 14), 80: (44, 160, 44), 100: (214, 39, 40)}
legend_box_y = None
marks = {}
for nf, c in cols.items():
    m = (np.abs(im - np.array(c)).sum(2) < 60)
    m[:, :xl + 1] = False
    m[yb:, :] = False
    line = m.copy()
    m = nd.binary_opening(m, structure=np.ones((5, 5)))
    lab, n = nd.label(m)
    pts = []
    for i in range(1, n + 1):
        yy, xx = np.nonzero(lab == i)
        if len(yy) < 20:
            continue
        pts.append((float(np.polyval(ax, xx.mean())), float(np.polyval(ay, yy.mean())), len(yy)))
    pts.sort()
    marks[nf] = pts
    print("Nf %d markers:" % nf, " ".join("Eo %.1f d %.4f (%d px)" % p for p in pts))
    # line: thin pixels of this colour left after removing the markers, read at three Eo values
    line &= ~nd.binary_dilation(m, structure=np.ones((5, 5)))
    for eo in (80.0, 120.0, 200.0):
        xpix = (eo - ax[1]) / ax[0]
        col = int(round(xpix))
        ys = np.nonzero(line[:, col - 1:col + 2].any(1))[0]
        ys = ys[(ys > yt + 2) & (ys < yb - 2)]
        if len(ys) == 0:
            print("   line at Eo %.0f: none" % eo)
            continue
        dl = float(np.polyval(ay, ys.mean()))
        ub = fr_viana(eo, nf)
        th = brentq(lambda d: 4 * nf / (3 * ub) * d ** 3 + 2 * d - 1, 1e-4, 0.5)
        print("   line at Eo %.0f: read %.4f, Brown with Viana %.4f, diff %+.4f" % (eo, dl, th, dl - th))

print()
at = {}
for nf, pts in marks.items():
    near = [p for p in pts if abs(p[0] - 100) < 3]
    if near:
        at[nf] = near[0][1]
        print("Eo 100, Nf %d: %.4f" % (nf, at[nf]))
nfs = np.array(sorted(k for k in at if k <= 80), float)
dv = np.array([at[int(k)] for k in nfs])
lin = float(np.interp(50, nfs[:2], dv[:2]))
quad = float(np.polyval(np.polyfit(nfs, dv, 2), 50))
# ratio to Brown with Viana, interpolated, times Brown with Viana at Nf 50
def brown(nf, eo=100.0):
    ub = fr_viana(eo, nf)
    return brentq(lambda d: 4 * nf / (3 * ub) * d ** 3 + 2 * d - 1, 1e-4, 0.5)
ratios = dv / np.array([brown(k) for k in nfs])
rat = float(np.interp(50, nfs[:2], ratios[:2])) * brown(50)
print("Nf 50: linear %.4f quadratic %.4f ratio-to-Brown %.4f" % (lin, quad, rat))
print("Brown with Viana at Nf 40 60 80 50:", ["%.4f" % brown(k) for k in (40, 60, 80, 50)])
print("ratios A&M/Brown at 40 60 80:", np.round(ratios, 4))
