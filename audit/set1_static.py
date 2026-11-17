"""Data audit, set 1: the static bubble runs (case S), recomputed from the raw OpenFOAM files.

Runs: V1a_static, D_coarseAR1, D_square, D_tol, D_visc10, D_T1prop, D_T3prop, S_rhoST.
Definitions are the ones behind the paper's numbers (tb_post, b11_srho.py, b2chk.py):
  mean maximum speed  sample mean of max|U| from postProcessing/maxSpeed over 0.3 <= t <= 1.0001 s
  liquid peak         max|U| over cells with alpha > 0.99, per written snapshot
  gas peak            max|U| over cells with alpha < 0.01, per written snapshot
  jump                tb_post.laplace_jump (mean p inside minus outside, centre x = 2, radius 0.25, margin 0.3)
  Laplace value       2 sigma / R_eff, R_eff from the gas volume (tb_post.effective_radius, wedge half-angle from the mesh)
  smeared share       gas in cells with alpha > 0.5, as a share of all gas (metrics.py smear)
"to 1 s" = the written snapshots with 0 < t <= 1.0001 s. Everything printed to 4 decimal places.

python -I set1_static.py   (writes set1_static_out.md and set1_static_snap.csv next to this file)
"""
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(CODE))
import tb_post  # noqa: E402

RUNS = Path.home() / "cfdml" / "runs_p2"
HERE = Path(__file__).resolve().parent
NAMES = ["V1a_static", "D_coarseAR1", "D_square", "D_tol", "D_visc10", "D_T1prop", "D_T3prop", "S_rhoST"]


def scal(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(", s)
    return np.array(s[m.end():s.index(")", m.end())].split(), dtype=float)


def vec(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<vector>\s*(\d+)\s*\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    return np.array(body.split(), dtype=float).reshape(n, 3)


def points(case):
    s = (case / "constant/polyMesh/points").read_text()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    return np.array(body.split(), dtype=float).reshape(n, 3)


def props(case):
    t = (case / "constant/transportProperties").read_text()
    sigma = float(re.search(r"sigma\s+([0-9.eE+-]+)", t)[1])
    w = re.search(r"water\s*\{[^}]*nu\s+([0-9.eE+-]+);[^}]*rho\s+([0-9.eE+-]+);", t)
    a = re.search(r"air\s*\{[^}]*nu\s+([0-9.eE+-]+);[^}]*rho\s+([0-9.eE+-]+);", t)
    nul, rhol, nug, rhog = float(w[1]), float(w[2]), float(a[1]), float(a[2])
    g = (case / "constant/g").read_text()
    gv = [float(x) for x in re.search(r"value\s*\(([^)]*)\)", g)[1].split()]
    return dict(sigma=sigma, nu_l=nul, rho_l=rhol, nu_g=nug, rho_g=rhog, g=gv)


def mesh_facts(case, C):
    P = points(case)
    r = np.round(np.hypot(P[:, 1], P[:, 2]), 6)
    x = np.round(P[:, 0], 6)
    rf = np.unique(r)
    xf = np.unique(x)
    yz = P[P[:, 1] > 1e-9]
    half = math.degrees(math.atan(np.max(np.abs(yz[:, 2]) / yz[:, 1])))
    dr = np.diff(rf)
    dx = np.diff(xf)
    return dict(n_cells=len(C), n_rad=len(rf) - 1, n_ax=len(xf) - 1, length=xf[-1] - xf[0], r_wall=rf[-1],
                dr_wall=dr[-1], dr_min=dr.min(), dr_mean=dr.mean(), dx=dx.mean(), dx_spread=dx.max() - dx.min(),
                ar_wall=dx.mean() / dr[-1], half_angle=half, faces=rf)


def monitor(case):
    t, mx = [], []
    for line in open(case / "postProcessing/maxSpeed/0/fieldMinMax.dat"):
        g = re.match(r"\s*([0-9.e+-]+)\s+mag\(U\)\s+([0-9.e+-]+)\s+\([^)]*\)\s+([0-9.e+-]+)", line)
        if g:
            t.append(float(g[1]))
            mx.append(float(g[3]))
    return np.array(t), np.array(mx)


def snapshots(case):
    ds = [d for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.]+", d.name)]
    return sorted(ds, key=lambda d: float(d.name))


def f4(v):
    return "nan" if v is None or (isinstance(v, float) and math.isnan(v)) else "%.4f" % v


out = []
snaprows = []
R = {}
for name in NAMES:
    case = RUNS / name
    C = vec(case / "0/C")
    V = scal(case / "0/V")
    pr = props(case)
    mf = mesh_facts(case, C)
    t_m, m_m = monitor(case)
    rows = []
    for d in snapshots(case):
        t = float(d.name)
        a = scal(d / "alpha.water")
        G = tb_post.gas_volume(a, V)
        if t == 0:
            rows.append(dict(t=0.0, G=G, xc=tb_post.centroid_x(C[:, 0], a, V)))
            continue
        U = vec(d / "U")
        p = scal(d / "p")
        sp = np.linalg.norm(U, axis=1)
        gas, liq = a < 0.01, a > 0.99
        itf = ~gas & ~liq
        reff = tb_post.effective_radius(a, V, mf["half_angle"])
        j = tb_post.laplace_jump(C[:, 0], C[:, 1], a, p, 2.0, 0.25)
        gas_all = np.sum((1 - a) * V)
        rows.append(dict(t=t, G=G, xc=tb_post.centroid_x(C[:, 0], a, V), liq=sp[liq].max(), gas=sp[gas].max(),
                         jump=(np.nan if j is None else j), reff=reff, target=2 * pr["sigma"] / reff,
                         smear=100 * np.sum(((1 - a) * V)[a > 0.5]) / gas_all))
    R[name] = dict(pr=pr, mf=mf, t_m=t_m, m_m=m_m, rows=rows)
    for r in rows:
        snaprows.append([name] + ["%.6f" % r.get(k, np.nan) if k in r else "" for k in
                                  ("t", "liq", "gas", "jump", "reff", "target", "G", "xc", "smear")])
    print("read", name, len(rows), "snapshots", file=sys.stderr)

with open(HERE / "set1_static_snap.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run", "t", "liquid_peak", "gas_peak", "jump", "R_eff", "laplace", "gas_volume", "centroid_x", "smeared_pct"])
    w.writerows(snaprows)


def window(name, lo, hi):
    return [r for r in R[name]["rows"] if r["t"] > 0 and lo <= r["t"] <= hi]


def mon_mean(name, lo=0.3, hi=1.0001):
    t, m = R[name]["t_m"], R[name]["m_m"]
    k = (t >= lo) & (t <= hi)
    return float(m[k].mean()), float(np.trapezoid(m[k], t[k]) / (t[k][-1] - t[k][0])), int(k.sum()), float(m[k].max()), float(m.max()), float(t[np.argmax(m)])


out.append("## Set 1. Static bubble S (OpenFOAM), recomputed 14 Nov from the raw run files\n")
out.append("Script `code/analysis/audit/set1_static.py`, per-snapshot values in `set1_static_snap.csv`.\n")
out.append("### Runs, properties and meshes\n")
out.append("| run | sigma | nu_l | rho_l | rho_g | g | Eo (g = 1) | Nf (g = 1) | cells | radial | axial | length | dx | dr at wall | dx/dr at wall | mean dr | half-angle (deg) | maxDeltaT | half capillary limit |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for name in NAMES:
    pr, mf = R[name]["pr"], R[name]["mf"]
    mu_l = pr["nu_l"] * pr["rho_l"]
    drho = pr["rho_l"] - pr["rho_g"]
    Eo = drho * 1.0 / pr["sigma"]
    Nf = math.sqrt(drho * pr["rho_l"]) / mu_l
    cd = (RUNS / name / "system/controlDict").read_text()
    mdt = float(re.search(r"maxDeltaT\s+([0-9.eE+-]+)", cd)[1])
    faces = mf["faces"]
    near = np.diff(faces)[faces[:-1] < 0.25 * 1.2]
    dxc = min(near.min(), mf["dx"])
    cap = math.sqrt(0.5 * (pr["rho_l"] + pr["rho_g"]) * dxc ** 3 / (2 * math.pi * pr["sigma"]))
    out.append("| %s | %.5g | %.6g | %g | %g | %s | %s | %s | %d | %d | %d | %s | %s | %s | %s | %s | %s | %.6g | %.6g |" % (
        name, pr["sigma"], pr["nu_l"], pr["rho_l"], pr["rho_g"], "(%g %g %g)" % tuple(pr["g"]), f4(Eo), f4(Nf),
        mf["n_cells"], mf["n_rad"], mf["n_ax"], f4(mf["length"]), f4(mf["dx"]), f4(mf["dr_wall"]), f4(mf["ar_wall"]),
        f4(mf["dr_mean"]), f4(mf["half_angle"]), mdt, 0.5 * cap))

out.append("\n### To t = 1 s (the paper's static table window)\n")
out.append("| run | monitor samples 0.3-1 s | mean max speed (m/s) | trapezoid mean | max in 0.3-1 s | peak of whole monitor (t) | snapshots | liquid peak min | liquid peak max | liquid peak mean | gas peak min | gas peak max | jump min (Pa) | jump max | jump mean | Laplace mean | jump mean error (%) |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
S = {}
for name in ["V1a_static", "D_coarseAR1", "D_square", "D_tol", "D_visc10", "S_rhoST"]:
    mm, mt, n, mx, pk, tpk = mon_mean(name)
    w = window(name, 0.0, 1.0001)
    liq = np.array([r["liq"] for r in w])
    gas = np.array([r["gas"] for r in w])
    jp = np.array([r["jump"] for r in w])
    lap = np.mean([r["target"] for r in w])
    S[name] = dict(mm=mm, liq=liq, jp=jp, lap=lap, gas=gas)
    out.append("| %s | %d | %s | %s | %s | %s (%s) | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        name, n, f4(mm), f4(mt), f4(mx), f4(pk), f4(tpk), len(w), f4(liq.min()), f4(liq.max()), f4(liq.mean()),
        f4(gas.min()), f4(gas.max()), f4(np.nanmin(jp)), f4(np.nanmax(jp)), f4(np.nanmean(jp)), f4(lap),
        f4(100 * (np.nanmean(jp) / lap - 1))))

out.append("\n### Whole runs\n")
out.append("| run | snapshots | t first to last | jump mean (Pa) | jump sd (population) | jump sd (sample) | jump min | jump max | n jump defined | Laplace mean | error (%) | liquid peak min | liquid peak max | liquid peak slope over t (m/s per s) | gas peak min | gas peak max | gas volume drift t=0 to end (%) | centroid x min | centroid x max | centroid x end | smeared % end |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for name in NAMES:
    w = window(name, 0.0, 99.0)
    jp = np.array([r["jump"] for r in w])
    liq = np.array([r["liq"] for r in w])
    tt = np.array([r["t"] for r in w])
    gas = np.array([r["gas"] for r in w])
    lap = np.mean([r["target"] for r in w])
    G0 = R[name]["rows"][0]["G"]
    G1 = R[name]["rows"][-1]["G"]
    xc = np.array([r["xc"] for r in w])
    slope = np.polyfit(tt, liq, 1)[0] if len(tt) > 2 else float("nan")
    out.append("| %s | %d | %s to %s | %s | %s | %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        name, len(w), f4(tt[0]), f4(tt[-1]), f4(np.nanmean(jp)), f4(np.nanstd(jp)), f4(np.nanstd(jp, ddof=1)),
        f4(np.nanmin(jp)), f4(np.nanmax(jp)), int(np.isfinite(jp).sum()), f4(lap), f4(100 * (np.nanmean(jp) / lap - 1)),
        f4(liq.min()), f4(liq.max()), f4(slope), f4(gas.min()), f4(gas.max()), f4(100 * (G1 / G0 - 1)),
        f4(xc.min()), f4(xc.max()), f4(xc[-1]), f4(w[-1]["smear"])))

out.append("\n### Ratios the paper uses\n")
v, c, s = S["V1a_static"], S["D_coarseAR1"], S["D_square"]
mm_full = {n: float(R[n]["m_m"][(R[n]["t_m"] >= 0.3)].mean()) for n in ("V1a_static", "S_rhoST")}
lines = [
    ("Mean max speed, stretched / 96 by 12 uniform", v["mm"] / c["mm"]),
    ("Mean max speed, stretched / 384 by 48 uniform", v["mm"] / s["mm"]),
    ("Liquid peak, stretched / 96 by 12: min over min", v["liq"].min() / c["liq"].min()),
    ("Liquid peak, stretched / 96 by 12: max over max", v["liq"].max() / c["liq"].max()),
    ("Liquid peak, stretched / 96 by 12: mean over mean", v["liq"].mean() / c["liq"].mean()),
    ("Cells, 96 by 12 / stretched", R["D_coarseAR1"]["mf"]["n_cells"] / R["V1a_static"]["mf"]["n_cells"]),
    ("Mean radial cell, 96 by 12 / stretched", R["D_coarseAR1"]["mf"]["dr_mean"] / R["V1a_static"]["mf"]["dr_mean"]),
    ("Axial cell, stretched / 384 by 48", R["V1a_static"]["mf"]["dx"] / R["D_square"]["mf"]["dx"]),
    ("Axial cell, stretched / 96 by 12", R["V1a_static"]["mf"]["dx"] / R["D_coarseAR1"]["mf"]["dx"]),
    ("D_tol against V1a, mean max speed change (%)", 100 * (S["D_tol"]["mm"] / v["mm"] - 1)),
    ("D_tol against V1a, mean liquid peak change (%)", 100 * (S["D_tol"]["liq"].mean() / v["liq"].mean() - 1)),
    ("D_tol against V1a, mean jump change (%)", 100 * (np.nanmean(S["D_tol"]["jp"]) / np.nanmean(v["jp"]) - 1)),
    ("Liquid peak, V1a / D_visc10: mean over mean", v["liq"].mean() / S["D_visc10"]["liq"].mean()),
    ("Liquid peak, V1a / D_visc10: min over min", v["liq"].min() / S["D_visc10"]["liq"].min()),
    ("Liquid peak, V1a / D_visc10: max over max", v["liq"].max() / S["D_visc10"]["liq"].max()),
    ("Mean max speed, V1a / S_rhoST, 0.3 to 1 s", v["mm"] / S["S_rhoST"]["mm"]),
    ("Mean max speed, V1a / S_rhoST, 0.3 s to the end", mm_full["V1a_static"] / mm_full["S_rhoST"]),
]
for lab, val in lines:
    out.append("- %s: %s" % (lab, f4(val)))

out.append("\n### Gas centroid of V1a_static and S_rhoST per snapshot (start x = 2)\n")
out.append("| t | V1a xc | S_rhoST xc | V1a jump | S_rhoST jump | V1a liquid peak | S_rhoST liquid peak | V1a smeared % | S_rhoST smeared % |")
out.append("|---|---|---|---|---|---|---|---|---|")
for a, b in zip(R["V1a_static"]["rows"], R["S_rhoST"]["rows"]):
    if a["t"] == 0:
        continue
    out.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (f4(a["t"]), f4(a["xc"]), f4(b["xc"]), f4(a["jump"]), f4(b["jump"]),
                                                             f4(a["liq"]), f4(b["liq"]), f4(a["smear"]), f4(b["smear"])))

(HERE / "set1_static_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
