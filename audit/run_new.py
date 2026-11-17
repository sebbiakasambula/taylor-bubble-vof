"""Per-snapshot metrics for the Fluent T2 runs of 12 Nov (A6, A7, A8, coarse, fine), with the frozen
metrics.snapshot_metrics, exactly as run_fl.py does for the medium baseline.

Step to flow time: the .out report rows where they exist (they stop early in every run except the coarse one), and
past them the step and flow-time pairs stored in the .dat.h5 files (current step and the autosave history). Between
consecutive known points the mean step must equal the cap to 1e-9 relative, so t = t_known + (step - step_known) * cap
is exact there. The script stops if that check fails.

python -I run_new.py [A6 A7 A8 coarse fine]
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import csv
import ctypes
import glob
import re
import sys
from pathlib import Path

import h5py
import numpy as np

P3 = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(P3))
import metrics  # noqa: E402

try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x00004000)  # below normal
except Exception:  # not Windows
    pass

sys.path.insert(0, str(P3))
from config import FLUENT as FD  # noqa: E402  (the Fluent exports of the test runs sit here too)
OUT = Path(__file__).resolve().parent / "flcopy"
WALL, CLEAR = 0.5, 1.5
RUNS = {  # name: (folder, snapshot prefix, dt cap)
    "A6": ("T2_A6", "A6_s-snap", 0.00085),
    "A7": ("T2_A7", "A7_s-snap", 0.00085),
    "A8": ("T2_A8", "A8_s-snap", 0.00085),
    "coarse": ("T2_Grid_study/Coarse", "Coarse_s-snap", 0.001407),
    "fine": ("T2_Grid_study/Fine", "T2_Taylor_bubble_Fine_Mesh", 0.000516),
    "T3D1": ("T3_D1", "D1_s-snap", 0.00026),
    "T3D2": ("T3_D2", "D2_s-snap", 0.00026),
}


def out_rows(folder):
    rows = {}
    for line in open(folder / "maximum-velocity-rfile.out"):
        s = line.split()
        try:
            rows[int(s[0])] = float(s[2])
        except (ValueError, IndexError):
            pass
    return rows


def dat_anchors(folder):
    anc = {}
    for p in glob.glob(str(folder / "*.dat.h5")):
        dv = h5py.File(p, "r")["settings/Data Variables"][0]
        dv = dv.decode() if isinstance(dv, bytes) else str(dv)
        now = int(re.search(r"\(time-step (\d+)\)", dv).group(1))
        for s, t in zip(re.findall(r"\(time-step \. (\d+)\)", dv), re.findall(r"\(flow-time \. ([0-9.eE+-]+)\)", dv)):
            if int(s) > now:  # a history entry later than the file itself is left from another run (T3 D2: 23000 in the -7 files)
                print("  skipped stale anchor %s (t %s) in %s, file at step %d" % (s, t, Path(p).name, now))
                continue
            anc[int(s)] = float(t)
        anc[now] = float(re.search(r"\(flow-time ([0-9.eE+-]+)\)", dv).group(1))
    return anc


def time_map(folder, cap):
    out = out_rows(folder)
    anc = dat_anchors(folder)
    for s, t in anc.items():  # anchors inside the .out range must agree with it
        if s in out:
            assert abs(out[s] - t) < 1e-9, (folder, s, out[s], t)
    known = dict(anc)
    last_out = max(out)
    known[last_out] = out[last_out]
    ks = sorted(k for k in known if k >= last_out)
    for a, b in zip(ks, ks[1:]):
        mean = (known[b] - known[a]) / (b - a)
        assert abs(mean / cap - 1) < 1e-9, ("step not at the cap", folder, a, b, mean)

    allk = dict(out)
    allk.update(anc)
    sk = sorted(allk)

    def t_of(step):
        if step in out:
            return out[step], "out"
        if step < last_out:  # inside a gap of the .out rows (T3 D2): bracket by the known points either side
            if step in allk:
                return allk[step], "gap"
            a = max(k for k in sk if k <= step)
            b = min(k for k in sk if k >= step)
            mean = (allk[b] - allk[a]) / (b - a)
            assert abs(mean / cap - 1) < 1e-9, ("step not at the cap in .out gap", folder, a, b, mean)
            return allk[a] + (step - a) * cap, "gap"
        a = max(k for k in ks if k <= step)
        return known[a] + (step - a) * cap, "dat"
    return t_of, last_out, ks


def run(name):
    folder, prefix, cap = RUNS[name]
    folder = FD / folder
    t_of, last_out, ks = time_map(folder, cap)
    snaps = {}
    for p in glob.glob(str(folder / (prefix + "-*"))):
        tail = p.rsplit("-", 1)[-1]
        if tail.isdigit():
            snaps[int(tail)] = Path(p)
    rows, src = [], {"out": 0, "dat": 0, "gap": 0}
    for s in sorted(snaps):
        t, how = t_of(s)
        src[how] += 1
        d = np.loadtxt(snaps[s], delimiter=",", skiprows=1)
        m = metrics.snapshot_metrics(d[:, 1], d[:, 2], 1.0 - d[:, 9], np.hypot(d[:, 4], d[:, 5]), d[:, 8], WALL, CLEAR)
        rows.append((s, t, m))
    with open(OUT / (("m_F_T3_%s.csv" % name[2:]) if name.startswith("T3") else ("m_F_T2_%s.csv" % name)), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "t"] + list(metrics.KEYS))
        for s, t, m in rows:
            w.writerow(["%d" % s, repr(t)] + [("%d" % m[k]) if isinstance(m[k], (int, np.integer)) else
                                               (m[k] if isinstance(m[k], str) else repr(float(m[k]))) for k in metrics.KEYS])
    print("%-6s %d snapshots, steps %d..%d, t %.4f..%.4f, cells %d; time from .out %d, from .dat anchors %d "
          "(.out ends at step %d; known points past it %s; in .out gaps %d)" % (name, len(rows), rows[0][0], rows[-1][0], rows[0][1],
                                                               rows[-1][1], len(d), src["out"], src["dat"], last_out, ks, src["gap"]))


if __name__ == "__main__":
    for n in (sys.argv[1:] or list(RUNS)):
        run(n)
