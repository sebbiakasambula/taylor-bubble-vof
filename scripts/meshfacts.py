"""Mesh facts for the paper mesh table from polyMesh points (radial and axial edges), and gas drift of the medium-mesh runs."""
import math, re
import numpy as np
from pathlib import Path
from config import RUNS as R
def scal(f):
    s = Path(f).read_text(); m = re.search(r'internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(', s)
    return np.fromstring(s[m.end():s.index(')', m.end())], sep=' ', count=int(m[1]))
def points(case):
    s = (R / case / 'constant/polyMesh/points').read_text()
    m = re.search(r'\n(\d+)\s*\n\(', s); n = int(m[1]); b = s[m.end():]; b = b[:b.index('\n)')].replace('(', ' ').replace(')', ' ')
    return np.fromstring(b, sep=' ', count=3 * n).reshape(n, 3)
c25 = math.cos(math.radians(2.5))
for case in ['V4_coarse', 'V4_medium', 'V4_fine', 'VA_square', 'V1a_static', 'D_coarseAR1', 'D_square']:
    P = points(case)
    front = P[P[:, 2] > 1e-9]
    r = np.unique(np.round(np.hypot(front[:, 1], front[:, 2]), 7))   # radius of each vertex ring
    x = np.unique(np.round(P[:, 0], 7))
    dr = np.diff(r); dx = np.diff(x)
    ar = dx.mean() / dr
    infilm = int(np.sum(r[1:] - dr / 2 > 0.5 - 0.0615))   # cell centres (mid-radius) within 0.0615 of the wall at r = 0.5
    print(f"{case:12s} radial {len(r)-1:4d} axial {len(x)-1:5d} dx {dx.mean():.5f} (spread {dx.max()-dx.min():.1e}) "
          f"dr axis {dr[0]:.5f} wall {dr[-1]:.5f} ratio wall/axis {dr[-1]/dr[0]:.3f} AR max {ar.max():.2f} min {ar.min():.2f} r_wall {r[-1]:.6f} in film {infilm}")
for case in ['V4_medium', 'T2_prod', 'T3_prod']:
    c = R / case
    ts = sorted([d.name for d in c.iterdir() if re.fullmatch(r'[0-9.]+', d.name)], key=float)
    V = scal(c / '0/V')
    G = {t: ((1 - scal(c / t / 'alpha.water')) * V).sum() for t in (ts[0], ts[1], ts[-1])}
    print(case, 'drift from', ts[0], 'to', ts[-1], '%.5f %%' % (100 * (G[ts[-1]] / G[ts[0]] - 1)),
          '| from', ts[1], '%.5f %%' % (100 * (G[ts[-1]] / G[ts[1]] - 1)))
