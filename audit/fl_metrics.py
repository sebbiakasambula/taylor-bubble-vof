"""Data audit: per-snapshot metrics of the Fluent baseline runs, recomputed from the raw cell-centre exports with the
frozen p3/metrics.py, as p3/run_fl.py does (columns 1 x, 2 y, 4 and 5 velocity, 8 volume, 9 gas fraction; wall at
r = 0.5). Snapshot files are chosen as run_fl.sources() does.

Times are rebuilt from the run's own report files, independently of the f_t*.npy arrays run_fl.py used: the .out row
of the step where there is one, otherwise the step and flow-time pairs in the .dat.h5 files, with the mean step between
known points checked against the step cap to 1e-9 (fl2/run_new.time_map). The cap is the step size of the last .out rows.

python -I fl_metrics.py T1 T1sq T1_40it T2 T3      writes audit/m/F_<name>.csv
"""
import csv
import glob
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
# run_new.py sits in this folder
import metrics  # noqa: E402
import run_new  # noqa: E402  (time_map, out_rows, dat_anchors)

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import FLUENT as FD  # noqa: E402
WALL = 0.5


def files_in(folder, prefix, max_step=None):
    out = {}
    for p in glob.glob(str(FD / folder / (prefix + "-*"))):
        tail = Path(p).name.rsplit("-", 1)[-1]
        if tail.isdigit() and Path(p).name.startswith(prefix + "-") and (max_step is None or int(tail) <= max_step):
            assert int(tail) not in out
            out[int(tail)] = Path(p)
    return out


def cap_of(folder):
    rows = run_new.out_rows(folder)
    s = sorted(rows)[-200:]
    d = np.diff([rows[k] for k in s]) / np.diff(s)
    return float(np.median(d))


# name: list of (folder, prefix, max_step), clear
RUNS = {
    "T1": ([("T1_Case", "s_snap", None)], 2.0),
    "T1sq": ([("T1_square", "s_snap", 6000), ("T1_square_w_completion", "s_snap_completion", None)], 2.0),
    "T1_40it": ([("40_iterations", "T1_Taylor_bubble", None)], 2.0),
    "T2": ([("T2_case", "T2_Taylor_bubble", None)], 1.5),
    "T3": ([("T3_Case", "T3_Taylor_bubble", None)], 1.5),
}


def fmt(v):
    if isinstance(v, (int, np.integer)):
        return "%d" % v
    if isinstance(v, str):
        return v
    return repr(float(v))


SP = HERE.parent / "data"   # the stored f_t*.npy arrays
ARR = {"T1": "f_t1.npy", "T1sq": "f_t1sq.npy", "T2": "f_t2.npy", "T3": "f_t3.npy"}


def run(name):
    """Times: the original arrays (step, t) where they exist, each checked against the .out row, a .dat.h5 anchor or,
    off those, the cap stepping from the nearest known point (runs held at the cap). Without an array (T1_40it) the
    times come from run_new.time_map."""
    parts, clear = RUNS[name]
    t0 = time.perf_counter()
    A = {int(a): float(b) for a, b in np.load(SP / ARR[name])[:, :2]} if name in ARR else None
    rows, chk = [], dict(out=0, anchor=0, cap=0, unchecked=0, worst=0.0)
    for folder, prefix, mx in parts:
        fdir = FD / folder
        src = files_in(folder, prefix, mx)
        out = run_new.out_rows(fdir)
        anc = run_new.dat_anchors(fdir)
        cap = cap_of(fdir)
        known = dict(anc); known.update(out)
        ks = np.array(sorted(known))
        tmap = None
        for s in sorted(src):
            if A is None:
                if tmap is None:
                    tmap = run_new.time_map(fdir, cap)[0]
                t = out[s] if s in out else tmap(s)[0]
            else:
                if s not in A:
                    raise SystemExit("%s: step %d has a file but no array time" % (name, s))
                t = A[s]
            if s in out:
                d = abs(t - out[s]); chk["out"] += 1
            elif s in anc:
                d = abs(t - anc[s]); chk["anchor"] += 1
            else:
                k = int(ks[np.argmin(np.abs(ks - s))])
                pred = known[k] + (s - k) * cap
                d = abs(t - pred)
                if d < 1e-6:
                    chk["cap"] += 1
                else:
                    chk["unchecked"] += 1
                    print("  %s step %d: array t %.6f, cap stepping from step %d gives %.6f (step not at the cap there)" % (name, s, t, k, pred))
                    d = 0.0
            chk["worst"] = max(chk["worst"], d)
            dd = np.loadtxt(src[s], delimiter=",", skiprows=1)
            m = metrics.snapshot_metrics(dd[:, 1], dd[:, 2], 1.0 - dd[:, 9], np.hypot(dd[:, 4], dd[:, 5]), dd[:, 8], WALL, clear)
            rows.append((s, t, m))
    rows.sort(key=lambda r: r[1])
    with open(HERE / "m" / ("F_%s.csv" % name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "t"] + list(metrics.KEYS))
        for s, t, m in rows:
            w.writerow(["%d" % s, repr(float(t))] + [fmt(m[k]) for k in metrics.KEYS])
    print("F_%s: %d snapshots, t %.4f to %.4f; times checked on .out %d, anchors %d, cap stepping %d, not checkable %d, largest difference %.3g s; %.1f s" % (
        name, len(rows), rows[0][1], rows[-1][1], chk["out"], chk["anchor"], chk["cap"], chk["unchecked"], chk["worst"], time.perf_counter() - t0), flush=True)


if __name__ == "__main__":
    for n in sys.argv[1:]:
        run(n)
