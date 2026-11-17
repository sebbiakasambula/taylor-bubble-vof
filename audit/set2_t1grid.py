"""Data audit, set 2: T1 grid study (V4_coarse, V4_medium, V4_fine) and the time-step test (V3_halfdt), OpenFOAM.

Per-snapshot metrics from audit/m/<run>.csv (of_metrics.py, frozen p3/metrics.py on the raw fields).
  Fr      least-squares slope of the nose over 5 - 1e-6 <= t <= end + 1e-3 (the paper's "5 s to the end")
  film    mean of the finite film values over 5 - 1e-6 <= t <= 20.001
  GCI     Celik et al. (2008), Fs 1.25, r = square root of the cell-count ratio, as fl2/results.py
  windows 5 s windows starting every 2.5 s from 5 s, at least 3 snapshots, as p3/fig_fr_split.py
  mesh    p3/meshfacts.py rules (vertex rings, cells in film = radial centres within 0.0615 of r = 0.5)
Everything to 4 decimal places. Writes set2_t1grid_out.md.
"""
import csv
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
LO = 5.0 - 1e-6


def load(name):
    rows = list(csv.DictReader(open(HERE / "m" / ("%s.csv" % name))))
    M = {}
    for k in rows[0]:
        try:
            M[k] = np.array([float(r[k]) for r in rows])
        except ValueError:
            M[k] = np.array([r[k] for r in rows])
    return M


def fit(t, y):
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    n = len(t)
    se = math.sqrt(np.sum(r ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)) if n > 2 else float("nan")
    return float(p[0]), se, n


def points(case):
    s = (RUNS / case / "constant/polyMesh/points").read_text()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    n = int(m[1])
    b = s[m.end():]
    b = b[:b.index("\n)")].replace("(", " ").replace(")", " ")
    return np.fromstring(b, sep=" ", count=3 * n).reshape(n, 3)


def meshfacts(case):
    P = points(case)
    front = P[P[:, 2] > 1e-9]
    r = np.unique(np.concatenate([[0.0], np.round(np.hypot(front[:, 1], front[:, 2]), 7)]))  # axis ring added (z = 0 there)
    x = np.unique(np.round(P[:, 0], 7))
    dr, dx = np.diff(r), np.diff(x)
    return dict(radial=len(r) - 1, axial=len(x) - 1, cells=(len(r) - 1) * (len(x) - 1), dx=dx.mean(),
                ar_wall=dx.mean() / dr[-1], infilm=int(np.sum(r[1:] - dr / 2 > 0.5 - 0.0615)),
                infilm_mid=int(np.sum(r[1:] - dr / 2 > 0.5 * math.cos(math.radians(2.5)) - 0.0615)),
                length=x[-1] - x[0], x_top=x[-1])


def celik(f3, f2, f1, r21, r32):
    e21, e32 = f2 - f1, f3 - f2
    s = np.sign(e32 / e21)
    p = 1.0
    for _ in range(500):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        pn = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(pn - p) < 1e-13:
            break
        p = pn
    ea21, ea32 = abs((f1 - f2) / f1), abs((f2 - f3) / f2)
    return dict(p=p, s=s, gci_fine=1.25 * ea21 / (r21 ** p - 1), gci_med_e32=1.25 * ea32 / (r32 ** p - 1),
                gci_med_pair21=1.25 * abs((f1 - f2) / f2) * r21 ** p / (r21 ** p - 1),
                ext=(r21 ** p * f1 - f2) / (r21 ** p - 1))


def f4(v):
    return "%.4f" % v


def wall_time(case):
    txt = (RUNS / case / "log.interFoam").read_text(errors="replace")
    ct = re.findall(r"ClockTime = ([0-9.]+) s", txt)
    et = re.findall(r"ExecutionTime = ([0-9.]+) s", txt)
    return float(ct[-1]), float(et[-1])


def max_dt(case):
    return float(re.search(r"maxDeltaT\s+([0-9.eE+-]+)", (RUNS / case / "system/controlDict").read_text())[1])


def courant(case, cap, limit=0.2):
    """Per time step: deltaT and the largest interface Courant number. Returns the share of steps at the cap and
    limit / (interface Courant) over the capped steps."""
    dts, cos = [], []
    co = None
    for line in open(RUNS / case / "log.interFoam", errors="replace"):
        if line.startswith("Interface Courant Number"):
            co = float(line.split("max:")[1])
        elif line.startswith("deltaT = ") and co is not None:
            dts.append(float(line.split("=")[1]))
            cos.append(co)
            co = None
    dts, cos = np.array(dts), np.array(cos)
    at_cap = dts / cap > 0.999  # the log prints the step rounded (0.00091205 against 0.000912233)
    ratio = limit / cos[at_cap & (cos > 0)]
    return dict(steps=len(dts), share_cap=100 * at_cap.mean(), ratio_median=np.median(ratio), ratio_p10=np.percentile(ratio, 10),
                ratio_min=ratio.min(), share_below=100 * (~at_cap).mean())


NAMES = ["V4_coarse", "V4_medium", "V4_fine", "V3_halfdt"]
D = {n: load(n) for n in NAMES}
MF = {n: meshfacts(n) for n in NAMES[:3]}
out = ["## Set 2. T1 grid study and time step (OpenFOAM), recomputed 14 Nov\n",
       "Script `code/analysis/audit/set2_t1grid.py`, per-snapshot tables `code/analysis/audit/m/<run>.csv`.\n",
       "### Meshes\n", "| mesh | radial | axial | cells | dx | dx/dr at wall | cells in film (r = 0.5) | same, wall at 0.4995 | domain length |",
       "|---|---|---|---|---|---|---|---|---|"]
for n in NAMES[:3]:
    m = MF[n]
    out.append("| %s | %d | %d | %d | %s | %s | %d | %d | %s |" % (n, m["radial"], m["axial"], m["cells"], f4(m["dx"]), f4(m["ar_wall"]),
                                                             m["infilm"], m["infilm_mid"], f4(m["length"])))

out += ["\n### Fr, film and gas per run\n",
        "| run | snapshots from 5 s | t last | Fr 5 s to end | s.e. | film 5 to 20 s (D) | film rows 5 to 20 s | snapshots after 20 s | with a film | gas drift t = 0 to end (%) | largest drift (%) | nose at end | outlet minus nose (D) | maxDeltaT | wall time (h) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
R = {}
for n in NAMES:
    M = D[n]
    t = M["t"]
    w = (t >= LO) & np.isfinite(M["nose"])
    fr, se, k = fit(t[w], M["nose"][w])
    wf = (t >= LO) & (t <= 20.001) & np.isfinite(M["film"])
    after = t > 20.001
    G = M["G"]
    drift = 100 * (G[-1] / G[0] - 1)
    big = 100 * np.max(np.abs(G / G[0] - 1))
    top = MF.get(n, MF["V4_medium"])["x_top"]
    wt = wall_time(n)[0] / 3600
    R[n] = dict(fr=fr, film=np.mean(M["film"][wf]))
    out.append("| %s | %d | %s | %s | %s | %s | %d | %d | %d | %s | %s | %s | %s | %.6g | %s |" % (
        n, k, f4(t[-1]), f4(fr), "%.6f" % se, f4(R[n]["film"]), int(wf.sum()), int(after.sum()), int(np.isfinite(M["film"][after]).sum()),
        f4(drift), f4(big), f4(M["nose"][-1]), f4(top - M["nose"][-1]), max_dt(n), f4(wt)))

r21 = math.sqrt(MF["V4_fine"]["cells"] / MF["V4_medium"]["cells"])
r32 = math.sqrt(MF["V4_medium"]["cells"] / MF["V4_coarse"]["cells"])
c = celik(R["V4_coarse"]["fr"], R["V4_medium"]["fr"], R["V4_fine"]["fr"], r21, r32)
cf = celik(R["V4_coarse"]["film"], R["V4_medium"]["film"], R["V4_fine"]["film"], r21, r32)
out += ["\n### Grid convergence (Celik et al. 2008, Fs 1.25)\n",
        "- r21 = %s, r32 = %s" % (f4(r21), f4(r32)),
        "- Fr: coarse %s, medium %s, fine %s; coarse to medium %s%%, medium to fine %s%%" % (
            f4(R["V4_coarse"]["fr"]), f4(R["V4_medium"]["fr"]), f4(R["V4_fine"]["fr"]),
            f4(100 * (R["V4_medium"]["fr"] / R["V4_coarse"]["fr"] - 1)), f4(100 * (R["V4_fine"]["fr"] / R["V4_medium"]["fr"] - 1))),
        "- Fr: %s, order %s, GCI fine %s%%, GCI medium from e32 %s%%, GCI medium from the fine-medium pair %s%%, extrapolated %s" % (
            "monotone" if c["s"] > 0 else "oscillatory", f4(c["p"]), f4(100 * c["gci_fine"]), f4(100 * c["gci_med_e32"]),
            f4(100 * c["gci_med_pair21"]), f4(c["ext"])),
        "- Film: coarse %s, medium %s, fine %s; coarse to medium %s%%, medium to fine %s%%, order %s" % (
            f4(R["V4_coarse"]["film"]), f4(R["V4_medium"]["film"]), f4(R["V4_fine"]["film"]),
            f4(100 * (R["V4_medium"]["film"] / R["V4_coarse"]["film"] - 1)), f4(100 * (R["V4_fine"]["film"] / R["V4_medium"]["film"] - 1)),
            f4(cf["p"]))]

out += ["\n### Fr over 5 s windows (start every 2.5 s, at least 3 snapshots)\n"]
for n in NAMES[:3]:
    M = D[n]
    t, y = M["t"], M["nose"]
    ws = []
    for a in np.arange(5, t[-1], 2.5):
        m = (t >= a) & (t <= a + 5 + 1e-3) & np.isfinite(y)
        if m.sum() >= 3:
            ws.append((a + 2.5, np.polyfit(t[m], y[m], 1)[0], int(m.sum())))
    dev = [100 * (v / R[n]["fr"] - 1) for _, v, _ in ws]
    out.append("- %s: %s. Largest deviation from the 5 s to end Fr: %s%%; last window (centre %s, %d snapshots) %s%%" % (
        n, ", ".join("%s: %s" % (f4(cw), f4(v)) for cw, v, _ in ws), f4(max(dev, key=abs)), f4(ws[-1][0]), ws[-1][2], f4(dev[-1])))

cap = max_dt("V4_medium")
cm = courant("V4_medium", cap)
wt_m, wt_h = wall_time("V4_medium")[0] / 3600, wall_time("V3_halfdt")[0] / 3600
out += ["\n### Time step\n",
        "- V4_medium maxDeltaT %.6g s, V3_halfdt maxDeltaT %.6g s (ratio %s)" % (cap, max_dt("V3_halfdt"), f4(max_dt("V3_halfdt") / cap)),
        "- V4_medium: %d steps, %s%% at the cap. At the capped steps, 0.2 / largest interface Courant number: median %s, 10th percentile %s, smallest %s" % (
            cm["steps"], f4(cm["share_cap"]), f4(cm["ratio_median"]), f4(cm["ratio_p10"]), f4(cm["ratio_min"])),
        "- Wall time on 2 processes (ClockTime): V4_medium %s h, V3_halfdt %s h" % (f4(wt_m), f4(wt_h)),
        "- Fr V3_halfdt %s against V4_medium %s: %s%%" % (f4(R["V3_halfdt"]["fr"]), f4(R["V4_medium"]["fr"]),
                                                       f4(100 * (R["V3_halfdt"]["fr"] / R["V4_medium"]["fr"] - 1)))]

(HERE / "set2_t1grid_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
