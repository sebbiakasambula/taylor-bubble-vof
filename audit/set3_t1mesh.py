"""Data audit, set 3: T1 mesh shape. VA_square (OpenFOAM square mesh) and Fluent T1 square, against V4_medium and
Fluent T1 on the stretched mesh.

1. The regenerated t1mesh_table.csv (cmp/cmp_table.py run on the recomputed per-snapshot tables, identical to the stored
   table) printed to 4 decimal places.
2. The numbers of the section that the table does not hold, computed here from audit/m/*.csv and the raw files.
Writes set3_t1mesh_out.md.
"""
import csv
import glob
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# run_new.py sits in this folder
import run_new  # noqa: E402

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import FLUENT as FD  # noqa: E402
LO = 5.0 - 1e-6
WALL = 0.5 * math.cos(math.radians(2.5))


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
    return float(p[0]), math.sqrt(np.sum(r ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)), n


def fr(M, lo, hi=1e9):
    t = M["t"]
    w = (t >= lo) & (t <= hi + 1e-3) & np.isfinite(M["nose"])
    return fit(t[w], M["nose"][w])


def f4(v):
    return "%.4f" % v


out = ["## Set 3. T1 mesh shape (OpenFOAM and Fluent), recomputed 14 Nov\n",
       "Per-snapshot tables recomputed from the raw files (audit/m), identical to the stored ones. `cmp_table.py` rerun on them "
       "(`code/analysis/audit/cmp/`) gives t1mesh_table.csv, cmp_table.csv, base_5to84.csv and wake_t3.csv identical to the stored "
       "tables (largest difference 0). Window 5 to 20 s unless stated.\n",
       "### t1mesh_table (regenerated), 4 decimal places\n",
       "| id | statistic | OF stretched | OF square | F stretched | F square |", "|---|---|---|---|---|---|"]
T = list(csv.DictReader(open(HERE / "cmp" / "t1mesh_table.csv")))
order, cells = [], {}
for r in T:
    key = (r["id"], r["statistic"])
    if key not in cells:
        order.append(key)
        cells[key] = {}
    v = r["value"]
    try:
        v = "%.4f" % float(v)
    except ValueError:
        pass
    cells[key][r["run"]] = v + ("" if r.get("n") in (None, "") else " (n %s)" % r["n"])
for key in order:
    c = cells[key]
    out.append("| %s | %s | %s | %s | %s | %s |" % (key[0], key[1], c.get("OF stretched", ""), c.get("OF square", ""),
                                                   c.get("F stretched", ""), c.get("F square", "")))

M = {n: load(n) for n in ("V4_medium", "VA_square", "V4_fine", "F_T1", "F_T1sq")}
sq5, sq10, fi5, fi10 = fr(M["VA_square"], LO), fr(M["VA_square"], 10 - 1e-6), fr(M["V4_fine"], LO), fr(M["V4_fine"], 10 - 1e-6)
GCI_FINE = 1.0577  # set 2
out += ["\n### Square mesh in OpenFOAM, whole run\n",
        "- VA_square Fr 5 s to end %s (n %d); 10 s to end %s (n %d)" % (f4(sq5[0]), sq5[2], f4(sq10[0]), sq10[2]),
        "- V4_fine Fr 5 s to end %s; 10 s to end %s" % (f4(fi5[0]), f4(fi10[0])),
        "- Square against fine, 5 s to end: %s%%, which is %s times the fine-grid GCI of %s%%" % (
            f4(100 * (sq5[0] / fi5[0] - 1)), f4(100 * (sq5[0] / fi5[0] - 1) / GCI_FINE), f4(GCI_FINE)),
        "- Square against fine, 10 s to end: %s%%" % f4(100 * (sq10[0] / fi10[0] - 1))]

# mesh facts of the square mesh, meshfacts rules with the axis ring
def points(case):
    s = (RUNS / case / "constant/polyMesh/points").read_text()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    n = int(m[1])
    b = s[m.end():]
    b = b[:b.index("\n)")].replace("(", " ").replace(")", " ")
    return np.fromstring(b, sep=" ", count=3 * n).reshape(n, 3)


mf = {}
for case in ("V4_medium", "VA_square"):
    P = points(case)
    front = P[P[:, 2] > 1e-9]
    r = np.unique(np.concatenate([[0.0], np.round(np.hypot(front[:, 1], front[:, 2]), 7)]))
    x = np.unique(np.round(P[:, 0], 7))
    dr, dx = np.diff(r), np.diff(x)
    mf[case] = dict(radial=len(r) - 1, axial=len(x) - 1, dx=dx.mean(), ar=dx.mean() / dr[-1],
                    infilm=int(np.sum(r[1:] - dr / 2 > 0.5 - 0.0615)), faces=r)
out += ["- VA_square mesh: %d radial by %d axial, %d cells, dx/dr at wall %s, cells in film %d. Axial cell of V4_medium / VA_square = %s; cells in film V4_medium %d" % (
    mf["VA_square"]["radial"], mf["VA_square"]["axial"], mf["VA_square"]["radial"] * mf["VA_square"]["axial"], f4(mf["VA_square"]["ar"]),
    mf["VA_square"]["infilm"], f4(mf["V4_medium"]["dx"] / mf["VA_square"]["dx"]), mf["V4_medium"]["infilm"])]

# film across the cells 2.5 D behind the nose at 19.7 s (B10 method)
def scal(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(", s)
    return np.fromstring(s[m.end():s.index(")", m.end())], sep=" ", count=int(m[1]))


def vec(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<vector>\s*(\d+)\s*\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    return np.fromstring(body, sep=" ", count=3 * n).reshape(n, 3)


for case in ("V4_medium", "VA_square"):
    C = vec(RUNS / case / "0/C")
    a = scal(RUNS / case / "19.700352/alpha.water")
    Mx = M[case]
    nose = float(Mx["nose"][np.argmin(np.abs(Mx["t"] - 19.700352))])
    xc = nose - 2.5
    cols = np.unique(np.round(C[:, 0], 6))
    cols = cols[np.abs(cols - xc) <= 0.0625 + 1e-9]
    th, cu = [], []
    faces = mf[case]["faces"] * math.cos(math.radians(2.5))  # mid-plane radii of the vertex rings
    for c in cols:
        k = np.abs(np.round(C[:, 0], 6) - c) < 1e-9
        y, al = C[k, 1], a[k]
        o = np.argsort(-y)
        y, al = y[o], al[o]
        j = int(np.argmax(al < 0.5))
        yc = y[j - 1] + (0.5 - al[j - 1]) * (y[j] - y[j - 1]) / (al[j] - al[j - 1])
        t_ = WALL - yc
        th.append(t_)
        # cells across the film: whole cells between the wall and the crossing plus the fraction of the cut cell
        n_full = np.sum(faces[:-1] >= yc) - 0  # rings at or above the crossing, inside the wall
        lo_ring = faces[faces < yc].max()
        hi_ring = faces[faces >= yc].min()
        cu.append(np.sum(faces[:-1] >= yc) + (hi_ring - yc) / (hi_ring - lo_ring))  # full cells inside the film plus the cut cell
    out.append("- %s at 19.700352 s: nose %s, window centre %s, %d columns, film %s D (min %s, max %s), cells across %s" % (
        case, f4(nose), f4(xc), len(cols), f4(np.mean(th)), f4(min(th)), f4(max(th)), f4(np.mean(cu))))

# Fluent: starting gas, time step on the square mesh, code differences
G1, G2 = M["F_T1"]["G"][0], M["F_T1sq"]["G"][0]
out += ["\n### Fluent square mesh\n",
        "- Gas volume, first snapshot: stretched %.9f, square %.9f, difference %s%%" % (G1, G2, f4(100 * (G2 / G1 - 1)))]
for folder, fname in (("T1_square", "maximum-velocity-rfile.out"), ("T1_square_w_completion", "maximum-velocity-rfile_6001.out")):
    rows = {}
    for line in open(FD / folder / fname):
        s = line.split()
        try:
            rows[int(s[0])] = float(s[2])
        except (ValueError, IndexError):
            pass
    st = np.array(sorted(rows))
    tt = np.array([rows[k] for k in st])
    k = (np.diff(st) == 1) & (st[:-1] > 200)
    dt = np.diff(tt)[k]
    out.append("- %s/%s: step size over %d consecutive rows after step 200: median %.6g s, 10th to 90th percentile %.6g to %.6g s" % (
        folder, fname, k.sum(), np.median(dt), np.percentile(dt, 10), np.percentile(dt, 90)))
v = {r["run"] + "|" + r["statistic"]: float(r["value"]) for r in T if r["id"] == "M1" and r["statistic"] == "Fr"}
of_st, of_sq, f_st, f_sq = v["OF stretched|Fr"], v["OF square|Fr"], v["F stretched|Fr"], v["F square|Fr"]
out += ["- Fluent against OpenFOAM, 5 to 20 s: stretched %s%%, square %s%%" % (f4(100 * (f_st / of_st - 1)), f4(100 * (f_sq / of_sq - 1))),
        "- Square against stretched, 5 to 20 s: OpenFOAM %s%%, Fluent %s%%" % (f4(100 * (of_sq / of_st - 1)), f4(100 * (f_sq / f_st - 1)))]
L = {r["run"]: float(r["value"]) for r in T if r["statistic"] == "L_b_mean"}
out.append("- Fluent bubble against OpenFOAM: stretched %s%%, square %s%%" % (f4(100 * (L["F stretched"] / L["OF stretched"] - 1)),
                                                                          f4(100 * (L["F square"] / L["OF square"] - 1))))
g = {r["run"]: float(r["value"]) for r in T if r["statistic"] in ("gap05_mean",)}
if g:
    out.append("- gap05 mean, square against stretched: OpenFOAM %s%%, Fluent %s%%" % (
        f4(100 * (g["OF square"] / g["OF stretched"] - 1)), f4(100 * (g["F square"] / g["F stretched"] - 1))))
ob = {r["run"]: float(r["value"]) for r in T if r["statistic"] == "out_b_mean"}
out.append("- out_b mean OpenFOAM / Fluent: stretched %s, square %s; square minus stretched: OpenFOAM %s, Fluent %s points" % (
    f4(ob["OF stretched"] / ob["F stretched"]), f4(ob["OF square"] / ob["F square"]), f4(ob["OF square"] - ob["OF stretched"]),
    f4(ob["F square"] - ob["F stretched"])))

(HERE / "set3_t1mesh_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out[-16:]))
