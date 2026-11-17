"""Data audit, set 5: T2 grid study in both codes. OpenFOAM T2_C1_coarse, T2_prod, T2_C2_fine (audit/m, recomputed);
Fluent coarse and fine (fl2/m_F_T2_coarse.csv, m_F_T2_fine.csv, computed from the raw exports by fl2/run_new.py) and
Fluent medium (audit/m/F_T2.csv). Window 5 to 8.4 s. 'matched' = the snapshot nearest each medium-run time of the same
code; 'all' = every snapshot in the window. Celik et al. (2008) as fl2/results.py. Writes set5_t2grid_out.md.
"""
import csv
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FL2 = HERE / "flcopy"
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
LO, HI = 5.0 - 1e-6, 8.4 + 1e-3
CELLS = {"coarse": 23667, "medium": 46080, "fine": 90048}
R21, R32 = math.sqrt(CELLS["fine"] / CELLS["medium"]), math.sqrt(CELLS["medium"] / CELLS["coarse"])


def load(p):
    rows = list(csv.DictReader(open(p)))
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
    return float(p[0]), (math.sqrt(np.sum(r ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)) if n > 2 else float("nan")), n


def celik(f3, f2, f1):
    e21, e32 = f2 - f1, f3 - f2
    s = np.sign(e32 / e21)
    p = 1.0
    for _ in range(500):
        q = math.log((R21 ** p - s) / (R32 ** p - s))
        pn = abs(math.log(abs(e32 / e21)) + q) / math.log(R21)
        if abs(pn - p) < 1e-13:
            break
        p = pn
    ea21, ea32 = abs((f1 - f2) / f1), abs((f2 - f3) / f2)
    return dict(p=p, s=s, gf=100 * 1.25 * ea21 / (R21 ** p - 1), gm=100 * 1.25 * ea32 / (R32 ** p - 1), ext=(R21 ** p * f1 - f2) / (R21 ** p - 1))


def f4(v):
    return "%.4f" % v


def window(M):
    return (M["t"] >= LO) & (M["t"] <= HI)


out = ["## Set 5. T2 grid study, OpenFOAM and Fluent, recomputed 14 Nov\n",
       "r21 = %s, r32 = %s (square root of the cell-count ratio).\n" % (f4(R21), f4(R32))]
G = {}
for code, files in (("OpenFOAM", {"coarse": HERE / "m/T2_C1_coarse.csv", "medium": HERE / "m/T2_prod.csv", "fine": HERE / "m/T2_C2_fine.csv"}),
                    ("Fluent", {"coarse": FL2 / "m_F_T2_coarse.csv", "medium": HERE / "m/F_T2.csv", "fine": FL2 / "m_F_T2_fine.csv"})):
    D = {k: load(v) for k, v in files.items()}
    med = D["medium"]
    tb = med["t"][window(med)]
    res = {}
    out += ["### %s\n" % code, "| mesh | rows in window | Fr all | s.e. | Fr matched | largest offset (s) | film all | film matched | gas at first row | gas against medium (%) |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for k in ("coarse", "medium", "fine"):
        M = D[k]
        w = window(M) & np.isfinite(M["nose"])
        fa = fit(M["t"][w], M["nose"][w])
        idx = [int(np.argmin(np.abs(M["t"] - x))) for x in tb]
        off = max(abs(M["t"][i] - x) for i, x in zip(idx, tb))
        fm = fit(M["t"][idx], M["nose"][idx])
        fla = np.nanmean(M["film"][window(M)])
        flm = np.nanmean(M["film"][idx])
        res[k] = dict(fa=fa[0], fm=fm[0], fla=fla, flm=flm, G=M["G"][0])
        out.append("| %s | %d | %s | %s | %s | %s | %s | %s | %.9f | %s |" % (k, int(w.sum()), f4(fa[0]), "%.6f" % fa[1], f4(fm[0]), f4(off),
                                                                            f4(fla), f4(flm), M["G"][0], f4(100 * (M["G"][0] / D["medium"]["G"][0] - 1))))
    G[code] = {k: res[k]["G"] for k in res}
    for key, lab in (("fm", "matched"), ("fa", "all rows")):
        c = celik(res["coarse"][key], res["medium"][key], res["fine"][key])
        out.append("- Fr %s: coarse to medium %s%%, medium to fine %s%%, %s, order %s, GCI fine %s%%, GCI medium %s%%, extrapolated %s" % (
            lab, f4(100 * (res["medium"][key] / res["coarse"][key] - 1)), f4(100 * (res["fine"][key] / res["medium"][key] - 1)),
            "monotone" if c["s"] > 0 else "oscillatory", f4(c["p"]), f4(c["gf"]), f4(c["gm"]), f4(c["ext"])))
    for key, lab in (("flm", "matched"), ("fla", "all rows")):
        fl = [res[k][key] for k in ("coarse", "medium", "fine")]
        out.append("- film %s: %s, %s, %s; changes %s%% and %s%%" % (lab, f4(fl[0]), f4(fl[1]), f4(fl[2]), f4(100 * (fl[1] / fl[0] - 1)), f4(100 * (fl[2] / fl[1] - 1))))
    if code == "OpenFOAM":
        long = load(HERE / "m/T2_prod.csv")
        w = (long["t"] >= LO) & np.isfinite(long["nose"])
        fl_ = fit(long["t"][w], long["nose"][w])
        c = celik(res["coarse"]["fm"], fl_[0], res["fine"]["fm"])
        out.append("- with the medium Fr from 5 s to the end (%s, n %d): order %s, GCI fine %s%%, GCI medium %s%%" % (
            f4(fl_[0]), fl_[2], f4(c["p"]), f4(c["gf"]), f4(c["gm"])))
        OFGCI = dict(matched=celik(res["coarse"]["fm"], res["medium"]["fm"], res["fine"]["fm"])["gm"], long=c["gm"])
    else:
        FGCI = celik(res["coarse"]["fm"], res["medium"]["fm"], res["fine"]["fm"])["gm"]
    out.append("")

out.append("- OpenFOAM GCI medium / Fluent GCI medium: matched %s, with the 5 s to end medium fit %s" % (f4(OFGCI["matched"] / FGCI), f4(OFGCI["long"] / FGCI)))
out.append("- OpenFOAM starting gas against Fluent medium: coarse %s%%, medium %s%%, fine %s%%" % tuple(
    f4(100 * (G["OpenFOAM"][k] / G["Fluent"]["medium"] - 1)) for k in ("coarse", "medium", "fine")))
out.append("- Fluent starting gas against Fluent medium: coarse %s%%, fine %s%%" % (
    f4(100 * (G["Fluent"]["coarse"] / G["Fluent"]["medium"] - 1)), f4(100 * (G["Fluent"]["fine"] / G["Fluent"]["medium"] - 1))))
for n in ("T2_C1_coarse", "T2_prod", "T2_C2_fine"):
    out.append("- %s maxDeltaT %s" % (n, re.search(r"maxDeltaT\s+([0-9.eE+-]+)", (RUNS / n / "system/controlDict").read_text())[1]))
out.append("- Fluent caps (fl2/run_new.py RUNS): coarse 0.001407, medium 0.00085, fine 0.000516")
(HERE / "set5_t2grid_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
