"""Data audit, set 7: the OpenFOAM configuration ladder at T2 and the related runs: T2_mompred, T2_A1_iso, T2_A2_plic,
T2_A3_outer3, T2_A2c_coarse (A2 on the coarse mesh), T2_B1_wedge10 (10 degree wedge), T2_A9_halfdt (half the time step,
never read before 14 Nov). Window 5 to 8.4 s. 'all' = every snapshot in the window; 'matched' = the snapshot nearest
each of the four T2_prod window times. References: Fluent 0.263163 (set 5), Abubakar and Matar 0.2531 (0.2514 to 0.2544,
code/analysis/abu), White and Beardmore 0.2671, Viana 0.250410. Writes set7_ladder_out.md.
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
LO, HI = 5.0 - 1e-6, 8.4 + 1e-3
FL, ABU, ABU_LO, ABU_HI, WB, VIANA = 0.263163, 0.2531, 0.2514, 0.2544, 0.2671, 0.250410


def f4(v):
    return "%.4f" % v


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
    return float(p[0]), (math.sqrt(np.sum(r ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)) if n > 2 else float("nan")), n


def wall(name):
    ct = re.findall(r"ClockTime = ([0-9.]+) s", (RUNS / name / "log.interFoam").read_text(errors="replace"))
    return float(ct[-1]) / 3600 if ct else float("nan")


def share0(name):
    from of_metrics import scal  # noqa: E402
    c = RUNS / name
    V = scal(c / "0/V")
    a = scal(c / "0/alpha.water")
    return float(((1 - a) * V).sum() / V.sum())


import sys  # noqa: E402
sys.path.insert(0, str(HERE))

B = load("T2_prod")
wb = (B["t"] >= LO) & (B["t"] <= HI)
TB = B["t"][wb]
base = fit(B["t"][wb], B["nose"][wb])[0]
base_film = float(np.mean(B["film"][wb]))
out = ["## Set 7. OpenFOAM configuration ladder at T2, recomputed 14 Nov\n",
       "Baseline T2_prod over 5 to 8.4 s: Fr %s from its 4 snapshots (%s), film %s.\n" % (f4(base), ", ".join(f4(x) for x in TB), f4(base_film)),
       "| run | change | Fr all (n) | Fr matched | gap to Fluent all (%) | phi all | vs Abubakar all (%) | vs White and Beardmore all (%) | vs Viana all (%) | gap matched (%) | vs Abubakar matched (%) | film all | largest gas drift (%) | starting gas share vs T2_prod (%) | wall time (h) |",
       "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
RES = {}
for name, ch in (("T2_prod", "baseline"), ("T2_mompred", "momentum predictor"), ("T2_A3_outer3", "predictor + 3 outer loops"),
                 ("T2_A1_iso", "interIsoFoam isoAlpha + predictor"), ("T2_A2_plic", "interIsoFoam plicRDF + predictor"),
                 ("T2_A2c_coarse", "A2 on the coarse mesh"), ("T2_B1_wedge10", "baseline on a 10 degree wedge"),
                 ("T2_A9_halfdt", "baseline at half the time step")):
    M = load(name)
    w = (M["t"] >= LO) & (M["t"] <= HI) & np.isfinite(M["nose"])
    fa = fit(M["t"][w], M["nose"][w])
    idx = [int(np.argmin(np.abs(M["t"] - x))) for x in TB]
    fm = fit(M["t"][idx], M["nose"][idx])[0]
    RES[name] = (fa[0], fm)
    out.append("| %s | %s | %s (%d) | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        name, ch, f4(fa[0]), fa[2], f4(fm), f4(100 * (FL / fa[0] - 1)), f4((fa[0] - base) / (FL - base)), f4(100 * (fa[0] / ABU - 1)),
        f4(100 * (fa[0] / WB - 1)), f4(100 * (fa[0] / VIANA - 1)), f4(100 * (FL / fm - 1)), f4(100 * (fm / ABU - 1)),
        f4(float(np.nanmean(M["film"][w]))), f4(100 * np.max(np.abs(M["G"] / M["G"][0] - 1))), f4(100 * (share0(name) / share0("T2_prod") - 1)),
        f4(wall(name))))

mp, a1, a2, a3 = (RES[k][0] for k in ("T2_mompred", "T2_A1_iso", "T2_A2_plic", "T2_A3_outer3"))
out += ["\n### Derived\n",
        "- A3 against the predictor run: %s%%; A1 against A2 (all rows): %s%%; A2 against the predictor run: phi %s" % (
            f4(100 * (a3 / mp - 1)), f4(100 * (a1 / a2 - 1)), f4((a2 - mp) / (FL - mp))),
        "- Fluent against Abubakar %s%%, against White and Beardmore %s%%; Abubakar band %s to %s" % (
            f4(100 * (FL / ABU - 1)), f4(100 * (FL / WB - 1)), f4(ABU_LO), f4(ABU_HI)),
        "- B1 against the baseline: all %s%%, matched %s%%; A9 against the baseline: all %s%%, matched %s%%" % (
            f4(100 * (RES["T2_B1_wedge10"][0] / base - 1)), f4(100 * (RES["T2_B1_wedge10"][1] / base - 1)),
            f4(100 * (RES["T2_A9_halfdt"][0] / base - 1)), f4(100 * (RES["T2_A9_halfdt"][1] / base - 1))),
        "- A9 against T2_prod over the same 18 snapshot times is not possible (T2_prod wrote 4 in the window); A9 maxDeltaT %s against %s" % (
            re.search(r"maxDeltaT\s+([0-9.eE+-]+)", (RUNS / "T2_A9_halfdt/system/controlDict").read_text())[1],
            re.search(r"maxDeltaT\s+([0-9.eE+-]+)", (RUNS / "T2_prod/system/controlDict").read_text())[1])]

# A3 residuals: medians of the initial residual of the first solve of each field in each outer loop, 5 to 8.4 s
steps, cur, t = [], None, None
for line in open(RUNS / "T2_A3_outer3/log.interFoam", errors="replace"):
    if line.startswith("Time = "):
        if cur:
            steps.append((t, cur))
        t = float(line.split("=")[1])
        cur = []
    elif line.startswith("PIMPLE: iteration") or line.startswith("PIMPLE: Iteration"):
        cur.append({})
    elif "Solving for" in line and cur:
        f = re.search(r"Solving for ([\w.]+), Initial residual = ([0-9.eE+-]+)", line)
        if f and f[1] not in cur[-1]:
            cur[-1][f[1]] = float(f[2])
if cur:
    steps.append((t, cur))
inwin = [s for tt, s in steps if LO <= tt <= HI and len(s) == 3]
out.append("- A3 log: %d steps with 3 outer loops in 5 to 8.4 s" % len(inwin))
for fld in ("p_rgh", "Uy", "Ux", "alpha.water"):
    med = [np.median([s[k][fld] for s in inwin if fld in s[k]]) if any(fld in s[k] for s in inwin) else float("nan") for k in range(3)]
    out.append("  %s median initial residual, loops 1 to 3: %.4g, %.4g, %.4g (first over third %.4f)" % (fld, med[0], med[1], med[2], med[0] / med[2]))
(HERE / "set7_ladder_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
