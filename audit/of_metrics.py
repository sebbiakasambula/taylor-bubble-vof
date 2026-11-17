"""Data audit: per-snapshot metrics of one OpenFOAM rising run, recomputed from the raw fields with the frozen
p3/metrics.py, exactly as p3/run_of.py does (same field readers, y rounded to 9 decimals, wall at the wedge mid-plane).

python -I of_metrics.py RUN CLEAR [RUN CLEAR ...]   writes audit/m/<RUN>.csv with every written snapshot, t = 0 included
(t = 0 gives the starting gas volume; its other metrics are those of the initial field).
Skips a run whose CSV already holds every time folder.
"""
import csv
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
import metrics  # noqa: E402

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
WALL = 0.5 * math.cos(math.radians(2.5))  # wedge mid-plane wall, 0.499524111, as run_of.py
OUT = HERE / "m"


def scal(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(", s)
    n = int(m[1])
    return np.fromstring(s[m.end():s.index(")", m.end())], sep=" ", count=n)


def vec(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<vector>\s*(\d+)\s*\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    return np.fromstring(body, sep=" ", count=3 * n).reshape(n, 3)


def fmt(v):
    if isinstance(v, (int, np.integer)):
        return "%d" % v
    if isinstance(v, str):
        return v
    return repr(float(v))


def half_angle(case):
    s = (case / "constant/polyMesh/points").read_text()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    P = np.fromstring(body, sep=" ", count=3 * n).reshape(n, 3)
    P = P[P[:, 1] > 1e-6]
    return math.degrees(math.atan(np.max(np.abs(P[:, 2]) / P[:, 1])))


def run(name, clear):
    case = RUNS / name
    ha = half_angle(case)
    wall = WALL if abs(ha - 2.5) < 1e-6 else 0.5 * math.cos(math.radians(ha))
    folders = sorted([(float(d.name), d) for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.]+", d.name)],
                     key=lambda p: p[0])
    out = OUT / ("%s.csv" % name)
    if out.exists():
        with open(out) as f:
            have = sum(1 for _ in f) - 1
        if have == len(folders):
            print("%s: cached, %d snapshots" % (name, have))
            return
    C = vec(case / "0/C")
    V = scal(case / "0/V")
    x, y = C[:, 0], np.round(C[:, 1], 9)
    try:
        metrics.structured_index(x, y)
        snapped = 0.0
    except ValueError:
        # write-precision noise (about 1e-7) splits an axial column in two on some meshes (V4_coarse): merge axial
        # centres closer than 1e-5 into one column at their mean, and report the largest shift
        ux = np.unique(x)
        grp = np.concatenate([[0], np.cumsum(np.diff(ux) > 1e-5)])
        means = np.array([ux[grp == g].mean() for g in range(grp[-1] + 1)])
        xs = means[grp[np.searchsorted(ux, x)]]
        snapped = float(np.max(np.abs(xs - x)))
        x = xs
        metrics.structured_index(x, y)
        print("%s: axial centres merged, largest shift %.3g" % (name, snapped), flush=True)
    t0 = time.perf_counter()
    rows = []
    for t, d in folders:
        alpha = scal(d / "alpha.water")
        Uf = d / "U"
        U = vec(Uf) if Uf.exists() and "nonuniform" in Uf.read_text()[:3000] else np.zeros((len(x), 3))
        m = metrics.snapshot_metrics(x, y, alpha, np.linalg.norm(U, axis=1), V, wall, clear)
        rows.append(dict(t=t, **m))
    OUT.mkdir(exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t"] + list(metrics.KEYS))
        for r in rows:
            w.writerow([repr(float(r["t"]))] + [fmt(r[k]) for k in metrics.KEYS])
    print("%s: %d snapshots, t %.6f to %.6f, %d cells, clear %.1f, half-angle %.4f, wall %.9f, %.1f s" % (
        name, len(rows), rows[0]["t"], rows[-1]["t"], len(x), clear, ha, wall, time.perf_counter() - t0), flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    for i in range(0, len(a), 2):
        run(a[i], float(a[i + 1]))
