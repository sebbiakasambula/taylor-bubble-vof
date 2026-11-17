"""Set 18. Validation uncertainty at T2 in the form of ASME V&V 20, for Fr and the film (17 Nov).

E = S - D for each configuration S and reference D, and u_val = sqrt(u_num^2 + u_fit^2 + u_wedge^2 + u_D^2), all in
per cent of the quantity. u_num is the grid convergence index of the medium mesh as it stands: Celik et al. with
Fs 1.25 where three meshes give an order (OpenFOAM baseline, Fluent), and Roache's two-grid index with Fs 3 for an
assumed order p = 2 and p = 1 where only the coarse and medium meshes exist (predictor, geometric advection with the
predictor). The time step scales with the cell size in these grid studies, so u_num covers both.

Values: Fr and film from the per-snapshot tables, 5 to 8.4 s. Three-grid studies use the snapshots nearest the medium
run's times (matched), as set 5 and the paper. Two-grid studies use all rows (18 on both meshes), as set 15.
"""
import csv
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LO, HI = 5.0 - 1e-6, 8.4 + 1e-3
R21, R32 = math.sqrt(90048 / 46080), math.sqrt(46080 / 23667)  # as set 5


def load(p):
    rows = list(csv.DictReader(open(p)))
    def num(v):
        try:
            return float(v)
        except ValueError:
            return float("nan")
    return {k: np.array([num(r[k]) for r in rows]) for k in ("t", "nose", "film")}


def fit(t, y):
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    return float(p[0]), math.sqrt(np.sum(r ** 2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2))


def win(M):
    return (M["t"] >= LO) & (M["t"] <= HI)


def all_rows(M):
    w = win(M) & np.isfinite(M["nose"])
    fr, se = fit(M["t"][w], M["nose"][w])
    return fr, se, float(np.nanmean(M["film"][w]))


def matched(M, tb):
    i = [int(np.argmin(np.abs(M["t"] - x))) for x in tb]
    fr, se = fit(M["t"][i], M["nose"][i])
    return fr, float(np.nanmean(M["film"][i]))


def celik(f3, f2, f1):
    e21, e32 = f2 - f1, f3 - f2
    Rr = e21 / e32
    s = np.sign(e32 / e21)
    p = 1.0
    for _ in range(500):
        q = math.log((R21 ** p - s) / (R32 ** p - s))
        pn = abs(math.log(abs(e32 / e21)) + q) / math.log(R21)
        if abs(pn - p) < 1e-13:
            break
        p = pn
    ea32 = abs((f2 - f3) / f2)
    return dict(p=p, R=Rr, gm=100 * 1.25 * ea32 / (R32 ** p - 1), ext=(R21 ** p * f1 - f2) / (R21 ** p - 1))


def roache2(fc, fm, p, r=R32):
    return 100 * 3.0 * abs((fc - fm) / fm) / (r ** p - 1)


M = {k: load(HERE / v) for k, v in {
    "base_c": "m/T2_C1_coarse.csv", "base_m": "m/T2_prod.csv", "base_f": "m/T2_C2_fine.csv",
    "pred_c": "m/T2_M1c_mompred_coarse.csv", "pred_m": "m/T2_mompred.csv",
    "geo_c": "m/T2_A2c_coarse.csv", "geo_m": "m/T2_A2_plic.csv",
    "fl_c": "flcopy/m_F_T2_coarse.csv", "fl_m": "m/F_T2.csv", "fl_f": "flcopy/m_F_T2_fine.csv",
    "wedge": "m/T2_B1_wedge10.csv", "A8": "flcopy/m_F_T2_A8.csv"}.items()}

REF_FR = {"Abubakar and Matar": (0.2531, math.hypot(100 * 0.0015 / 0.2531, 100 * 0.0005 / 0.2531)),
          "White and Beardmore": (0.2671, 0.9)}
REF_FILM = {"Abubakar and Matar": (0.1403, math.hypot(100 * 0.0002 / 0.1403, 100 * 0.0002 / 0.1403)),
            "Llewellin et al.": (0.1394, float("nan"))}

out = ["# Set 18. Validation uncertainty at T2 (ASME V&V 20 form), 17 Nov", "",
       "All in per cent of the quantity. 5 to 8.4 s, medium mesh unless stated.", ""]

# medium values (all rows, as the paper) and fit errors
res = {}
for key in ("base", "pred", "geo", "fl"):
    fr, se, film = all_rows(M[key + "_m"])
    res[key] = dict(fr=fr, se=se, film=film, u_fit=100 * se / fr)

# grid studies
tb_of = M["base_m"]["t"][win(M["base_m"])]
tb_fl = M["fl_m"]["t"][win(M["fl_m"])]
g = {}
for key, tb in (("base", tb_of), ("fl", tb_fl)):
    c3 = [matched(M[key + s], tb) for s in ("_c", "_m", "_f")]
    cfr = celik(c3[0][0], c3[1][0], c3[2][0])
    cfi = celik(c3[0][1], c3[1][1], c3[2][1])
    g[key] = dict(fr3=[x[0] for x in c3], film3=[x[1] for x in c3], cfr=cfr, cfi=cfi)
    out.append("- %s three meshes, matched: Fr %s, order %.4f, R %.4f, GCI medium %.4f%%, extrapolated %.4f; film %s, order %.4f, R %.4f, GCI medium %.4f%%, extrapolated %.4f" % (
        key, ", ".join("%.4f" % v for v in g[key]["fr3"]), cfr["p"], cfr["R"], cfr["gm"], cfr["ext"],
        ", ".join("%.4f" % v for v in g[key]["film3"]), cfi["p"], cfi["R"], cfi["gm"], cfi["ext"]))
# Fluent film: e21/e32 > 1 means the change grows under refinement (Celik: monotone divergence), so no order. Two-grid
# index of Roache from the larger of the two changes, the medium-fine pair, with r21.
fl_f3 = g["fl"]["film3"]
g["fl"]["film_u"] = {p: roache2(fl_f3[2], fl_f3[1], p, R21) for p in (2, 1)}
out.append("- Fluent film: changes coarse to medium %.4f%%, medium to fine %.4f%%; Roache two-grid from the medium-fine pair: p 2 %.4f%%, p 1 %.4f%%" % (
    100 * (fl_f3[1] / fl_f3[0] - 1), 100 * (fl_f3[2] / fl_f3[1] - 1), g["fl"]["film_u"][2], g["fl"]["film_u"][1]))
for key in ("pred", "geo"):
    frc, sec, filmc = all_rows(M[key + "_c"])
    g[key] = dict(frc=frc, filmc=filmc,
                  fr_u={p: roache2(frc, res[key]["fr"], p) for p in (2, 1)},
                  film_u={p: roache2(filmc, res[key]["film"], p) for p in (2, 1)})
    out.append("- %s two meshes, all rows: Fr coarse %.4f, medium %.4f (%+.4f%%), GCI Fs 3 p 2 %.4f%%, p 1 %.4f%%; film coarse %.4f, medium %.4f (%+.4f%%), GCI p 2 %.4f%%, p 1 %.4f%%" % (
        key, frc, res[key]["fr"], 100 * (frc / res[key]["fr"] - 1), g[key]["fr_u"][2], g[key]["fr_u"][1], filmc,
        res[key]["film"], 100 * (filmc / res[key]["film"] - 1), g[key]["film_u"][2], g[key]["film_u"][1]))

# other components
wfr, wse, wfilm = all_rows(M["wedge"])
u_wedge_fr = 100 * abs(wfr / res["base"]["fr"] - 1)
u_wedge_film = 100 * abs(wfilm / res["base"]["film"] - 1)
a8 = matched(M["A8"], tb_fl)
rep = 100 * abs(a8[0] / matched(M["fl_m"], tb_fl)[0] - 1)
out.append("- 10 degree wedge against the baseline: Fr %.4f%%, film %.4f%% (used as u_wedge for OpenFOAM)" % (u_wedge_fr, u_wedge_film))
out.append("- Fluent Courant 0.25 run (the step stayed at its cap, so a repeat), matched: %.4f%% (not added; inside u_fit)" % rep)
out.append("- references, Fr: " + "; ".join("%s %.4f (u_D %.2f%%)" % (k, *v) for k, v in REF_FR.items()))
out.append("- references, film: " + "; ".join("%s %.4f (u_D %.2f%%)" % (k, *v) for k, v in REF_FILM.items()))
out.append("")

labels = {"base": "OpenFOAM baseline", "pred": "OpenFOAM predictor", "geo": "OpenFOAM geometric and predictor", "fl": "Fluent"}
out += ["## Fr", "", "| configuration | Fr | u_num | u_fit | u_wedge | E vs A&M | u_val vs A&M | E vs W&B | u_val vs W&B |",
        "|---|---|---|---|---|---|---|---|---|"]
rows_fr = {}
for key in ("base", "pred", "geo", "fl"):
    r = res[key]
    uw = u_wedge_fr if key != "fl" else 0.0
    if key in ("base", "fl"):
        un = {"": g[key]["cfr"]["gm"]}
    else:
        un = {" (p 2)": g[key]["fr_u"][2], " (p 1)": g[key]["fr_u"][1]}
    for tag, u in un.items():
        cells = []
        for ref, (D, uD) in REF_FR.items():
            E = 100 * (r["fr"] / D - 1)
            uv = math.sqrt(u ** 2 + r["u_fit"] ** 2 + uw ** 2 + uD ** 2)
            cells += ["%+.2f" % E, "%.2f" % uv]
            rows_fr[(key, tag, ref)] = (E, uv)
        out.append("| %s%s | %.4f | %.2f | %.2f | %.2f | %s |" % (labels[key], tag, r["fr"], u, r["u_fit"], uw, " | ".join(cells)))

out += ["", "## Film", "", "| configuration | film | u_num | u_wedge | E vs A&M | u_val vs A&M | E vs Llewellin |", "|---|---|---|---|---|---|---|"]
for key in ("base", "pred", "geo", "fl"):
    r = res[key]
    uw = u_wedge_film if key != "fl" else 0.0
    if key == "base":
        un = {"": g[key]["cfi"]["gm"]}
    elif key == "fl":
        un = {" (p 2)": g[key]["film_u"][2], " (p 1)": g[key]["film_u"][1]}
    else:
        un = {" (p 2)": g[key]["film_u"][2], " (p 1)": g[key]["film_u"][1]}
    for tag, u in un.items():
        D, uD = REF_FILM["Abubakar and Matar"]
        E = 100 * (r["film"] / D - 1)
        uv = math.sqrt(u ** 2 + uw ** 2 + uD ** 2)
        El = 100 * (r["film"] / REF_FILM["Llewellin et al."][0] - 1)
        out.append("| %s%s | %.4f | %.2f | %.2f | %+.2f | %.2f | %+.2f |" % (labels[key], tag, r["film"], u, uw, E, uv, El))

out += ["", "Extrapolated (three meshes): OpenFOAM baseline Fr %.4f (%+.2f%% against A&M), film %.4f (%+.2f%%); Fluent Fr %.4f (%+.2f%%)" % (
    g["base"]["cfr"]["ext"], 100 * (g["base"]["cfr"]["ext"] / 0.2531 - 1), g["base"]["cfi"]["ext"],
    100 * (g["base"]["cfi"]["ext"] / 0.1403 - 1), g["fl"]["cfr"]["ext"], 100 * (g["fl"]["cfr"]["ext"] / 0.2531 - 1))]
(HERE / "set18_budget_out.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
