"""Per-snapshot metrics for one OpenFOAM run in $P2_RUNS/<case>.

python run_of.py --case V4_medium --t0 5 --t1 20 --clear 2.0 --out m_V4_medium.csv [--extra 20.3]

Snapshots with t0 - 1e-6 <= t <= t1 + 1e-3 are read, plus any --extra time folder outside the window.
--out is a file name inside the work folder ($P2_WORK), or an absolute path.
One CSV row per snapshot, t first.
"""
import argparse
import csv
import math
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import metrics

exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])

from config import RUNS, WORK
WALL = 0.5 * math.cos(math.radians(2.5))  # wedge mid-plane wall, 0.499524111


def fmt(v):
    """Counts as integers, text as text, everything else as a full-precision float."""
    if isinstance(v, (int, np.integer)):
        return "%d" % v
    if isinstance(v, str):
        return v
    return repr(float(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--t0", type=float, required=True)
    ap.add_argument("--t1", type=float, required=True)
    ap.add_argument("--clear", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--extra", type=float, action="append", default=[])
    a = ap.parse_args()

    case = RUNS / a.case
    folders = sorted([(float(d.name), d) for d in case.iterdir()
                      if d.is_dir() and re.fullmatch(r"[0-9.]+", d.name)], key=lambda p: p[0])
    pick = [(t, d) for t, d in folders if a.t0 - 1e-6 <= t <= a.t1 + 1e-3]
    for te in a.extra:
        m = [(t, d) for t, d in folders if abs(t - te) < 1e-6]
        if not m:
            sys.exit("no time folder at t = %g in %s" % (te, case))
        if m[0] not in pick:
            pick.append(m[0])
    pick.sort(key=lambda p: p[0])

    C = vec(case / "0/C")
    V = scal(case / "0/V")
    x, y = C[:, 0], np.round(C[:, 1], 9)

    t_start = time.perf_counter()
    rows = []
    for t, d in pick:
        alpha = scal(d / "alpha.water")
        U = vec(d / "U")
        m = metrics.snapshot_metrics(x, y, alpha, np.linalg.norm(U, axis=1), V, WALL, a.clear)
        rows.append(dict(t=t, **m))
    with open(WORK / a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t"] + list(metrics.KEYS))
        for r in rows:
            w.writerow([repr(float(r["t"]))] + [fmt(r[k]) for k in metrics.KEYS])
    print("%s: %d snapshots, t %.6f to %.6f, %d cells, %.1f s" % (
        a.case, len(rows), rows[0]["t"], rows[-1]["t"], len(x), time.perf_counter() - t_start))


if __name__ == "__main__":
    main()
