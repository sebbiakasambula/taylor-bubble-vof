"""Data audit, set 16: T2_M2f_mompred_fine, the momentum predictor (MULES, 1 outer loop) on the fine T2 mesh, the last
run of the predictor grid study (hint 2). The bubble splits, so no Fr is given; the run is compared snapshot by snapshot
with the fine baseline (T2_C2_fine) and A2 on the fine mesh (T2_A2f_fine). Reads audit/m/*.csv and audit/intact/*.csv.
Writes set16_m2f_out.md."""
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(folder, name):
    with open(HERE / folder / ("%s.csv" % name)) as f:
        return {round(float(r["t"]), 3): r for r in csv.DictReader(f)}


def num(r, k):
    try:
        v = float(r[k])
        return v if math.isfinite(v) else float("nan")
    except (KeyError, ValueError, TypeError):
        return float("nan")


runs = [("T2_C2_fine", "baseline, fine"), ("T2_M2f_mompred_fine", "predictor, fine"), ("T2_A2f_fine", "A2, fine")]
M = {n: load("m", n) for n, _ in runs}
I = {n: load("intact", n) for n, _ in runs}

out = ["## Set 16. T2_M2f_mompred_fine, the momentum predictor on the fine T2 mesh (hint 2)\n",
       "Ran 14 Nov 11:56 to 21:59, exit 0, to 8.4 s, 43 snapshots. Case files identical to T2_C2_fine except "
       "momentumPredictor yes (checked 14 Nov 19:14). Metrics: of_metrics.py (clear 1.5); regions: intact_check.py.\n"]

g0 = {n: num(M[n][0.0], "G") for n, _ in runs}
gE = {n: num(M[n][8.4], "G") for n, _ in runs}
out.append("Starting gas G: " + ", ".join("%s %.6f" % (lab, g0[n]) for n, lab in runs) +
           ". M2f against the fine baseline %+.4f%%. Gas drift 0 to 8.4 s: " % (100 * (g0["T2_M2f_mompred_fine"] / g0["T2_C2_fine"] - 1)) +
           ", ".join("%s %+.4f%%" % (lab, 100 * (gE[n] / g0[n] - 1)) for n, lab in runs) + ".\n")

out.append("| t | " + " | ".join("%s: regions, largest %%, axis nose, L_b" % lab for _, lab in runs) + " |")
out.append("|---|" + "---|" * len(runs))
for t in [0.2, 1.0, 2.0, 2.4, 2.6, 2.8, 3.0, 3.6, 3.8, 4.0, 5.0, 5.2, 6.0, 7.0, 8.4]:
    cells = []
    for n, _ in runs:
        m, i = M[n].get(t), I[n].get(t)
        if m is None or i is None:
            cells.append("-")
            continue
        cells.append("%s, %.2f, %.3f, %.3f" % (i["regions"], float(i["largest_pct"]), num(m, "nose"), num(m, "L_b")))
    out.append("| %.1f | " % t + " | ".join(cells) + " |")

m2, c2 = M["T2_M2f_mompred_fine"], M["T2_C2_fine"]
out.append("")
out.append("Axis nose speed before the first split (2.0 to 2.6 s): predictor %.3f D/s, fine baseline %.3f D/s." % (
    (num(m2[2.6], "nose") - num(m2[2.0], "nose")) / 0.6, (num(c2[2.6], "nose") - num(c2[2.0], "nose")) / 0.6))
fl = [(t, num(m2[t], "film"), num(c2[t], "film")) for t in (1.0, 1.6, 2.0) if t in m2 and t in c2]
out.append("Film before the split (D), predictor against baseline: " +
           ", ".join("t %.1f: %.4f / %.4f" % x for x in fl) + ".")
pk = [(t, num(m2[t], "pk_gas"), num(c2[t], "pk_gas")) for t in (1.0, 2.0, 2.4, 2.6) if t in m2 and t in c2]
out.append("Largest gas speed (m/s), predictor against baseline: " +
           ", ".join("t %.1f: %.3f / %.3f" % x for x in pk) + ".")
out.append("")
out.append("Reading (facts only): on the fine mesh the predictor run splits between 2.6 and 2.8 s and again between 3.6 "
           "and 3.8 s; the fine baseline stays one bubble; A2 splits once between 5.0 and 5.2 s. No fine Fr exists for "
           "the predictor configuration, so its grid study has two points (coarse, medium). Why is not known.")
text = "\n".join(out) + "\n"
(HERE / "set16_m2f_out.md").write_text(text, encoding="utf-8")
print(text)
