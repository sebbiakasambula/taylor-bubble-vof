import sys, numpy as np
from PIL import Image
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "scripts"))
from tb_references import fr_viana
im = np.array(Image.open('fig6a_raster.png').convert('RGB')).astype(int)
ax = np.polyfit([86,163,240,316,393,470,546], [0,50,100,150,200,250,300], 1)
ay = np.polyfit([56,118,180,242,304,366], [0.30,0.25,0.20,0.15,0.10,0.05], 1)
for eo in (80, 120, 200):
    x = int(round((eo - ax[1]) / ax[0]))
    for nf, c in ((40,(31,119,180)), (60,(255,127,14)), (80,(44,160,44))):
        w = []; 
        for xx in (x-1, x, x+1):
            col = np.abs(im[:, xx] - np.array(c)).sum(1); ys = np.nonzero(col < 80)[0]
            if len(ys): w.append(np.average(ys, weights=80 - col[ys]))
        y = np.mean(w); print("Eo %3d Nf %d line %.4f  Viana code %.4f  diff %+.4f" % (eo, nf, np.polyval(ay, y), fr_viana(eo, nf), np.polyval(ay, y) - fr_viana(eo, nf)))
