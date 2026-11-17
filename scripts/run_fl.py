"""Per-snapshot metrics for the Fluent runs (cell-centre CSV snapshots), with the frozen metrics.snapshot_metrics.

python run_fl.py [T1] [T1sq] [T2] [T3]      (no argument: all four, one after the other)

Fluent snapshot columns: 1 x, 2 y, 4 axial velocity, 5 radial velocity, 8 cell volume, 9 gas fraction (phase-2 vof).
Step and time are taken from columns 0 and 1 of the arrays f_t1.npy, f_t1sq.npy, f_t2.npy, f_t3.npy in the work folder ($P2_WORK).
One CSV row per array row, step and t first. The wall is at r = 0.5.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"          # single-thread BLAS, the laptop stays usable
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import csv
import ctypes
import glob
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
from config import WORK, FLUENT
sys.path.insert(0, str(HERE))
import numpy as np
import metrics

try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x00004000)  # BELOW_NORMAL
except Exception:
    pass

FD = FLUENT
WALL = 0.5


def step_of(p):
    """Step number if the file name ends in -<digits>, else None (skips .cas, .dat, .out files)."""
    tail = Path(p).name.rsplit("-", 1)[-1]
    return int(tail) if tail.isdigit() else None


def files_in(folder, prefix, max_step=None):
    out = {}
    for p in glob.glob(str(FD / folder / (prefix + "-*"))):
        s = step_of(p)
        if s is not None and Path(p).name.startswith(prefix + "-") and (max_step is None or s <= max_step):
            assert s not in out, (folder, s)
            out[s] = Path(p)
    return out


def sources(name):
    if name == "T1":
        return files_in("T1_Case", "s_snap")
    if name == "T1sq":
        a = files_in("T1_square", "s_snap", max_step=6000)
        b = files_in("T1_square_w_completion", "s_snap_completion")
        assert not set(a) & set(b), sorted(set(a) & set(b))
        return {**a, **b}
    if name == "T2":
        return files_in("T2_case", "T2_Taylor_bubble")
    if name == "T3":
        return files_in("T3_Case", "T3_Taylor_bubble")
    raise KeyError(name)


RUNS = {  # name: (array, clear, output)
    "T1": ("f_t1.npy", 2.0, "m_F_T1.csv"),
    "T1sq": ("f_t1sq.npy", 2.0, "m_F_T1sq.csv"),
    "T2": ("f_t2.npy", 1.5, "m_F_T2.csv"),
    "T3": ("f_t3.npy", 1.5, "m_F_T3.csv"),
}


def fmt(v):
    """Counts as integers, text as text, everything else as a full-precision float."""
    if isinstance(v, (int, np.integer)):
        return "%d" % v
    if isinstance(v, str):
        return v
    return repr(float(v))


def run(name):
    arr, clear, out = RUNS[name]
    A = np.load(WORK / arr)
    src = sources(name)
    missing = [int(s) for s in A[:, 0] if int(s) not in src]
    if missing:
        sys.exit("%s: steps in the array without a snapshot file: %s" % (name, missing))
    extra = sorted(set(src) - {int(s) for s in A[:, 0]})
    t_start = time.perf_counter()
    rows = []
    for s, t in zip(A[:, 0].astype(int), A[:, 1]):
        d = np.loadtxt(src[s], delimiter=",", skiprows=1)
        m = metrics.snapshot_metrics(d[:, 1], d[:, 2], 1.0 - d[:, 9], np.hypot(d[:, 4], d[:, 5]), d[:, 8], WALL, clear)
        rows.append((int(s), float(t), m))
    with open(WORK / out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "t"] + list(metrics.KEYS))
        for s, t, m in rows:
            w.writerow(["%d" % s, repr(t)] + [fmt(m[k]) for k in metrics.KEYS])
    print("%s: %d snapshots (steps %d to %d, t %.4f to %.4f), %d cells, clear %.1f, %.1f s; files not in the array: %s" % (
        name, len(rows), rows[0][0], rows[-1][0], rows[0][1], rows[-1][1], len(d), clear, time.perf_counter() - t_start, extra or "none"))


if __name__ == "__main__":
    for n in (sys.argv[1:] or list(RUNS)):
        run(n)
