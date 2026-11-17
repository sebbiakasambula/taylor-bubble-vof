"""Hysing et al. (2009) test case 1: benchmark quantities for the nine H1 runs, against Table XII p. 1275.

V_c is the gas-weighted mean velocity (bubbleU function object, every step). Its maximum and the time of the
maximum come from a parabola through the five samples around the largest sample. y_c(3) is the gas centroid at
t = 3 from the written fields. Circularity is pi d_a / P, with d_a the diameter of the circle of the gas area and P
the length of the alpha = 0.5 contour drawn through the cell centres of the uniform n x 2n grid. Its minimum comes
from a parabola through the three write times around the smallest value (writes every 0.05 s).
"""
import re
import sys
from pathlib import Path

import contourpy
import numpy as np

exec(open(Path(__file__).with_name('snap.py')).read().split("case=Path")[0])   # scal(), vec()
RUN = Path.home() / 'cfdml/runs_p2'
REF = dict(vmax=0.2419, tvmax=(0.921, 0.932), cmin=0.9012, yc3=1.081)


def series(p):
    t, v = [], []
    for l in open(p):
        if l.startswith('#'):
            continue
        a = l.replace('(', ' ').replace(')', ' ').split()
        t.append(float(a[0])); v.append([float(x) for x in a[1:]])
    return np.array(t), np.array(v)


def peak(t, y, k, half):
    i = int(k)
    s = slice(max(i - half, 0), i + half + 1)
    c = np.polyfit(t[s] - t[i], y[s], 2)
    tv = -c[1] / (2 * c[0])
    return t[i] + tv, np.polyval(c, tv)


def case(name):
    c = RUN / name
    n = int(name.rsplit('_', 1)[1]); h = 1.0 / n
    t, U = series(c / 'postProcessing/bubbleU/0/volFieldValue.dat')
    tv, vmax = peak(t, U[:, 1], np.argmax(U[:, 1]), 2)
    tg, G = series(c / 'postProcessing/gasVolume/0/volFieldValue.dat')
    drift = G[-1, 0] / G[0, 0] - 1
    C = vec(c / '0/C'); V = scal(c / '0/V')
    ix = np.floor(C[:, 0] * n).astype(int); iy = np.floor(C[:, 1] * n).astype(int)
    xg = (np.arange(n) + 0.5) * h; yg = (np.arange(2 * n) + 0.5) * h
    times = sorted([d for d in c.iterdir() if d.is_dir() and re.fullmatch(r'[0-9.]+', d.name) and d.name != '0'],
                   key=lambda d: float(d.name))
    ct, cv, yc3 = [], [], None
    for d in times:
        a = scal(d / 'alpha.water'); g = 1 - a
        area = (g * V).sum() / h
        grid = np.empty((2 * n, n)); grid[iy, ix] = a
        lines = contourpy.contour_generator(xg, yg, grid).lines(0.5)
        P = sum(np.sqrt((np.diff(L, axis=0) ** 2).sum(1)).sum() for L in lines)
        ct.append(float(d.name)); cv.append(np.pi * 2 * np.sqrt(area / np.pi) / P)
        if abs(float(d.name) - 3.0) < 1e-6:
            yc3 = (g * V * C[:, 1]).sum() / (g * V).sum()
    ct, cv = np.array(ct), np.array(cv)
    k = np.argmin(cv)
    tc, cmin = peak(ct, cv, k, 1)
    return dict(n=n, vmax=vmax, tv=tv, cmin=cmin, tc=tc, yc3=yc3, drift=drift, pieces=len(lines))


rows = {}
for cfg in ('mules', 'mulespred', 'plic'):
    for n in (40, 80, 160):
        rows[cfg, n] = r = case(f'H1_{cfg}_{n}')
        print(f"{cfg:9s} 1/{n:<3d} Vmax {r['vmax']:.5f} at {r['tv']:.3f} ({100*(r['vmax']/REF['vmax']-1):+.2f}%)"
              f" | circ min {r['cmin']:.5f} at {r['tc']:.2f} ({100*(r['cmin']/REF['cmin']-1):+.3f}%)"
              f" | yc(3) {r['yc3']:.5f} ({100*(r['yc3']/REF['yc3']-1):+.3f}%) | gas drift {100*r['drift']:+.4f}%"
              f" | contour pieces at t=3: {r['pieces']}")
np.save(Path(__file__).with_name('hysing_rows.npy'), rows, allow_pickle=True)
