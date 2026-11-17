"""Grid-study films from the raw fields: per-snapshot film (tail+0.5 to nose-2.0), then mean and median over
the 5 s to end window, with the wall at 0.5 (as used) and at the wedge mid-plane 0.499524."""
import re, math
import numpy as np
from pathlib import Path
exec(open(Path(__file__).with_name('snap.py')).read().split("case=Path")[0])
from config import RUNS, WORK
RUN = RUNS
WALL = 0.5 * math.cos(math.radians(2.5))
out = {}
for case in ('V4_coarse', 'V4_medium', 'V4_fine', 'VA_square'):
    C = vec(RUN / case / '0/C'); x, y = C[:, 0], np.round(C[:, 1], 9)
    print(case, 'max cell-centre y %.6f' % y.max())
    T, F = [], []
    for d in sorted([p for p in (RUN / case).iterdir() if re.fullmatch(r'[0-9.]+', p.name) and p.name != '0'], key=lambda p: float(p.name)):
        t = float(d.name)
        if t < 5:
            continue
        a = scal(d / 'alpha.water'); n = tb_post.nose_position(x, y, a); tl = tb_post.tail_position(x, y, a)
        if n is None or tl is None:
            T.append(t); F.append(np.nan); continue
        xi, ri = tb_post.interface_radius(x, y, a); w = (xi > tl + 0.5) & (xi < n - 2.0)
        T.append(t); F.append((0.5 - ri[w]).mean() if w.sum() >= 3 else np.nan)
    T, F = np.array(T), np.array(F); out[case] = (T, F)
    for lab, m in (('5-end', T >= 5), ('5-20', (T >= 5) & (T <= 20.001))):
        f = F[m]; ok = ~np.isnan(f)
        print('  %-6s n %2d/%2d  mean %.5f median %.5f | wall-corrected mean %.5f median %.5f' % (
            lab, ok.sum(), m.sum(), np.nanmean(f), np.nanmedian(f), np.nanmean(f) - (0.5 - WALL), np.nanmedian(f) - (0.5 - WALL)))
np.save(WORK / 'v_film.npy', np.array([out], dtype=object), allow_pickle=True)
def order(f3, f2, f1, N=(23667, 46080, 90048)):
    r21, r32 = math.sqrt(N[2] / N[1]), math.sqrt(N[1] / N[0]); e21, e32 = f2 - f1, f3 - f2; s = np.sign(e32 / e21); p = 2.0
    for _ in range(200):
        p = abs(math.log(abs(e32 / e21)) + math.log((r21 ** p - s) / (r32 ** p - s))) / math.log(r21)
    ea = abs(e21 / f1); return p, 100 * 1.25 * ea / (r21 ** p - 1)
for stat in (np.nanmean, np.nanmedian):
    for lab, lo, hi in (('5-end', 5, 1e9), ('5-20', 5, 20.001)):
        v = [stat(out[k][1][(out[k][0] >= lo) & (out[k][0] <= hi)]) for k in ('V4_coarse', 'V4_medium', 'V4_fine')]
        for wall_shift in (0.0, 0.5 - WALL):
            vv = [z - wall_shift for z in v]; p, g = order(*vv)
            print('%-9s %-6s wall %s  films %.5f %.5f %.5f  p %.3f  GCI fine %.3f%%' % (stat.__name__, lab, '0.5' if wall_shift == 0 else 'mid', *vv, p, g))
