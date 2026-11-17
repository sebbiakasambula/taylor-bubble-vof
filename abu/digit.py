"""Digitise Abubakar & Matar (arXiv 2109.09091) Fig. 6a markers from the embedded raster (585 x 440)."""
import numpy as np
from PIL import Image
from scipy import ndimage as nd
im = np.array(Image.open('fig6a_raster.png').convert('RGB')).astype(int)
H, W, _ = im.shape
dark = im.sum(2) < 150
# spines: longest dark column / row
colsum = dark.sum(0); rowsum = dark.sum(1)
xl = int(np.argmax(colsum[:W//2])); xr = W//2 + int(np.argmax(colsum[W//2:]))
yb = H//2 + int(np.argmax(rowsum[H//2:])); yt = int(np.argmax(rowsum[:H//2]))
print("spines x", xl, xr, "y", yt, yb)
# x ticks: dark pixels just below the bottom spine
xt = [x for x in range(xl, xr+1) if dark[yb+1:yb+5, x].all()]
# y ticks: dark pixels just left of the left spine
ytk = [y for y in range(yt, yb+1) if dark[y, xl-4:xl].all()]
def groups(v):
    g=[[v[0]]]
    for a in v[1:]:
        (g[-1].append(a) if a-g[-1][-1]<=1 else g.append([a]))
    return [np.mean(x) for x in g]
xt = groups(xt); ytk = groups(ytk)
print("x tick px", np.round(xt,2)); print("y tick px", np.round(ytk,2))
xv = np.arange(0, 301, 50)[:len(xt)]; yv = np.round(np.arange(0.30, 0.049, -0.05), 2)[:len(ytk)]
ax = np.polyfit(xt, xv, 1); ay = np.polyfit(ytk, yv, 1)
print("x fit resid", np.round(np.polyval(ax, xt)-xv, 4), "y fit resid", np.round(np.polyval(ay, ytk)-yv, 5))
cols = {'Nf40': (31,119,180), 'Nf60': (255,127,14), 'Nf80': (44,160,44)}
for name, c in cols.items():
    m = (np.abs(im - np.array(c)).sum(2) < 60); m[:, :xl+1] = False; m[yb:, :] = False
    m = nd.binary_opening(m, structure=np.ones((5,5)))   # removes the 1.5 px line, keeps discs
    lab, n = nd.label(m)
    pts = []
    for i in range(1, n+1):
        yy, xx = np.nonzero(lab == i)
        if len(yy) < 20: continue
        pts.append((np.polyval(ax, xx.mean()), np.polyval(ay, yy.mean()), len(yy)))
    pts.sort()
    print(name, " ".join("Eo %.1f Ub %.4f (%d px)" % p for p in pts))
