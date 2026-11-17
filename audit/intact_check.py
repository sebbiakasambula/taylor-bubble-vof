"""Did the T2 bubble stay in one piece? For each run and each snapshot up to 8.4 s, label the cells with alpha < 0.5 that
touch through faces, and give the number of regions, the share of the gas in the largest and second largest region, and
the axial extent and lowest radius of the largest region. Gas share uses (1 - alpha) V over the region cells against all
the gas in the domain. Axial centres closer than 1e-5 are merged first (write-precision noise on the coarse meshes).

python -I intact_check.py RUN [RUN ...]   writes audit/intact/<RUN>.csv and prints a summary line per run
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from of_metrics import RUNS, scal, vec  # noqa: E402

OUT = HERE / "intact"
OUT.mkdir(exist_ok=True)


def merged_index(v, tol=1e-5):
    u = np.unique(v)
    groups = np.concatenate([[0], np.cumsum(np.diff(u) > tol)])
    centres = np.array([u[groups == g].mean() for g in range(groups[-1] + 1)])
    return centres, groups[np.searchsorted(u, v)]


def run(name, tmax=8.4 + 1e-3):
    case = RUNS / name
    C = vec(case / "0/C")
    V = scal(case / "0/V")
    ux, ix = merged_index(C[:, 0])
    uy, iy = merged_index(np.round(C[:, 1], 9))
    assert len(ux) * len(uy) == len(V), (len(ux), len(uy), len(V))
    ts = sorted([d.name for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.]+", d.name)], key=float)
    ts = [t for t in ts if float(t) <= tmax and (case / t / "alpha.water").exists()]
    rows = []
    for t in ts:
        a = scal(case / t / "alpha.water")
        A = np.ones((len(ux), len(uy)))
        A[ix, iy] = a
        G = np.zeros_like(A)
        G[ix, iy] = (1 - a) * V
        lab, n = ndimage.label(A < 0.5)
        tot = G.sum()
        if n == 0:
            rows.append((float(t), 0, 0.0, 0.0, 100.0, np.nan, np.nan, np.nan))
            continue
        s = np.sort(ndimage.sum(G, lab, range(1, n + 1)))[::-1]
        k = 1 + int(np.argmax(ndimage.sum(G, lab, range(1, n + 1))))
        body = lab == k
        rr = np.where(body.any(axis=1))[0]
        rows.append((float(t), n, 100 * s[0] / tot, 100 * (s[1] if n > 1 else 0.0) / tot, 100 * (1 - s[0] / tot),
                     ux[rr.max()], ux[rr.max()] - ux[rr.min()], tot))
    with open(OUT / ("%s.csv" % name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "regions", "largest_pct", "second_pct", "outside_largest_pct", "top_x", "extent", "gas_total"])
        for r in rows:
            w.writerow(["%.6g" % r[0], r[1]] + ["%.6f" % v for v in r[2:]])
    R = np.array(rows, dtype=float)
    late = R[R[:, 0] >= 1.0]
    print("%-22s snaps %3d  max regions %3d  max second %.4f%%  max outside largest %.4f%%  extent at ~1,2,3,5,8.4 s %s"
          % (name, len(R), int(np.nanmax(R[:, 1])), np.nanmax(R[:, 3]), np.nanmax(late[:, 4]),
             " ".join("%.2f" % R[np.argmin(np.abs(R[:, 0] - x)), 6] for x in (1, 2, 3, 5, 8.4))))


if __name__ == "__main__":
    for nm in sys.argv[1:]:
        run(nm)
