"""Fr (nose fit) and film over 5 to 8.4 s for the T2 runs, to compare the density-weighted surface tension run
with the baseline and the predictor test. Usage: python t2rho.py CASE [CASE ...]"""
import re, sys
from pathlib import Path
import numpy as np
exec(open(Path(__file__).with_name('snap.py')).read().split("case=Path")[0])
from config import RUNS, WORK
RUN = RUNS
for name in sys.argv[1:]:
    case = RUN / name
    C = vec(case / '0/C'); x, y = C[:, 0], np.round(C[:, 1], 9)
    rows = []
    for d in sorted([p for p in case.iterdir() if re.fullmatch(r'[0-9.]+', p.name) and p.name != '0'], key=lambda p: float(p.name)):
        t = float(d.name)
        if t < 4.99 or t > 8.41 or not (d / 'alpha.water').exists():
            continue
        a = scal(d / 'alpha.water'); n = tb_post.nose_position(x, y, a); tl = tb_post.tail_position(x, y, a)
        xi, ri = tb_post.interface_radius(x, y, a); w = (xi > tl + 0.5) & (xi < n - 1.5)
        rows.append((t, n, (0.5 - ri[w]).mean() if w.sum() >= 3 else np.nan))
    R = np.array(rows); m = (R[:, 0] >= 5 - 1e-6) & (R[:, 0] <= 8.4 + 1e-3)
    fr = np.polyfit(R[m, 0], R[m, 1], 1)[0]
    print(f"{name:14s} n={m.sum():2d} t {R[m,0].min():.2f}-{R[m,0].max():.2f}  Fr nose {fr:.4f}  film {np.nanmean(R[m,2]):.4f}")
    np.save(WORK / (name + '_5to84.npy'), R)
