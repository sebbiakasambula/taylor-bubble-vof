"""T1 stretched against square, both codes, near 19.7 s: VTK files for ParaView and the numbers the text will quote."""
import math
from config import WORK
import numpy as np
from pathlib import Path
exec(open(Path(__file__).with_name('tovtk.py')).read().split("of('V1a_static'")[0])
W = 0.5 - 0.5 * math.cos(math.radians(2.5))
SNAP = {
    'T1_OF_str': ('of', 'V4_medium', '19.700352'),
    'T1_OF_sq': ('of', 'VA_square', '19.700352'),
    'T1_F_str': ('fl', 'T1_Case', 's_snap-21645'),
    'T1_F_sq': ('fl', 'T1_square_w_completion', 's_snap_completion-8190'),
}
def load(kind, case, snap):
    if kind == 'of':
        C = vec(RUN / case / '0/C'); V = scal(RUN / case / '0/V'); a = scal(RUN / case / snap / 'alpha.water'); U = vec(RUN / case / snap / 'U')
        return C[:, 0], C[:, 1], np.linalg.norm(U, axis=1), a, V, W
    d = np.loadtxt(FD / case / snap, delimiter=',', skiprows=1)
    return d[:, 1], d[:, 2], np.hypot(d[:, 4], d[:, 5]), 1 - d[:, 9], d[:, 8], 0.0
shapes = {}
for name, (kind, case, snap) in SNAP.items():
    x, y, sp, a, V, w = load(kind, case, snap)
    yr = np.round(y, 9); n = tb_post.nose_position(x, yr, a); tl = tb_post.tail_position(x, yr, a)
    write(name, x, y, sp, a, n)
    xi, ri = tb_post.interface_radius(x, yr, a); m = (xi > tl + 0.5) & (xi < n - 2.0)
    film = (0.5 - ri[m]).mean() - w
    gas = (1 - a) * V; wake = (x < tl - 0.05) & (x > tl - 2.0)
    print(f"{name:10s} nose {n:.3f} tail {tl:.3f} length {n-tl:.2f} film {film:.4f} | peak all {sp.max():.2f} liquid {sp[a>0.99].max():.2f} gas {sp[a<0.01].max():.2f} | wake gas {100*gas[wake].sum()/gas.sum():.2f}%  wake liquid median {np.median(sp[wake & (a>0.999)]):.3f}")
    shapes[name] = (xi - n, ri + (w if kind == 'of' else 0.0))   # OpenFOAM radius moved to a wall at 0.5 for the overlay
np.save(WORK / 't1_shapes.npy', np.array([shapes], dtype=object), allow_pickle=True)
