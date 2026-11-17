import re, numpy as np
from pathlib import Path
from config import RUNS as R
def scal(f):
    s = Path(f).read_text(); m = re.search(r'internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(', s)
    return np.fromstring(s[m.end():s.index(')', m.end())], sep=' ', count=int(m[1]))
for case in ['V4_medium', 'T2_prod', 'T3_prod']:
    c = R / case; V = scal(c / '0/V')
    ts = sorted([d.name for d in c.iterdir() if re.fullmatch(r'[0-9.]+', d.name)], key=float)
    G = np.array([((1 - scal(c / t / 'alpha.water')) * V).sum() for t in ts])
    d = 100 * (G / G[0] - 1); i = np.argmax(np.abs(d))
    print(case, 'end %.5f%%  max |d| %.5f%% at t=%s  (n %d)' % (d[-1], d[i], ts[i], len(ts)))
