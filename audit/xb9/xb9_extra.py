"""B9 supplement (not requested): where is the gas in T3_rhoST versus T3_prod?  Appends to checks/B9_t3check.md.
Reuses the functions of xb9.py (own code) by executing its definitions only."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
from pathlib import Path
HERE_ = Path(__file__).resolve().parent
src = (HERE_ / "xb9.py").read_text(encoding="utf-8").split("gs, rs = run(")[0]
exec(src)


def rows_for(case, lo, hi):
    g = Grid(case)
    out = []
    for t, d in time_folders(RUNS / case, lo, hi):
        a = read_field(d / "alpha.water", "scalar")
        A2 = a[g.idx2d]; gas2 = A2 < 0.5
        gv = (1.0 - A2) * g.V2; G = gv.sum()
        xc = float((gv * g.X2).sum() / G)                      # gas-volume weighted centroid in x
        cols = np.flatnonzero(gas2.any(axis=1))
        lab, n = ndimage.label(gas2, structure=S4)
        cg = ndimage.sum(gv * gas2, lab, index=np.arange(1, n + 1))
        nose, tail, _ = nose_and_tail(g.x, g.y, a)
        out.append(dict(t=t, name=d.name, xc=xc, xmax=float(g.ux[cols[-1]]), xmin=float(g.ux[cols[0]]), naxis=int(gas2[:, 0].sum()), ncomp=int(n),
                        big_pct=100.0 * float(cg.max()) / G, nose=(np.nan if nose is None else nose)))
    return out


R = {"T3_rhoST": rows_for("T3_rhoST", -1.0, 8.4 + 1e-3), "T3_prod": rows_for("T3_prod", -1.0, 9.0)}
L = ["", "## Supplement (not requested): where is the gas?", "",
     "Gas-volume centroid xc = sum((1-alpha) V x) / sum((1-alpha) V); xmax / xmin = largest / smallest x of any cell with alpha < 0.5; axis cells = number of axis-row cells with alpha < 0.5; big comp. = share of gas volume in the largest 4-connected alpha < 0.5 region (cells with alpha < 0.5 only, so below 100 % even for one bubble).", ""]
for case in ("T3_rhoST", "T3_prod"):
    rr = R[case]
    L += ["### %s" % case, "", "| t (s) | nose (xb5 rule) | centroid xc | xmax | xmin | axis cells | components | big comp. (% of gas) |", "|---|---|---|---|---|---|---|---|"]
    for r in rr:
        L.append("| %s | %s | %s | %s | %s | %d | %d | %s |" % (r["name"], f6(r["nose"]), f6(r["xc"]), f6(r["xmax"]), f6(r["xmin"]), r["naxis"], r["ncomp"], f6(r["big_pct"])))
    L.append("")
    for lo, hi, lab_ in ((5.0 - 1e-6, 8.4 + 1e-3, "window 5 to 8.4"), (-1.0, 2.6 + 1e-3, "t = 0 to 2.6")):
        sel = [r for r in rr if lo <= r["t"] <= hi]
        if len(sel) < 3:
            continue
        t = [r["t"] for r in sel]
        okn = [r for r in sel if np.isfinite(r["nose"])]
        fc = slope_se(t, [r["xc"] for r in sel]); fx = slope_se(t, [r["xmax"] for r in sel])
        fn = slope_se([r["t"] for r in okn], [r["nose"] for r in okn])
        L.append("- %s, %s: slope of centroid %s (s.e. %s, n %d); slope of xmax %s (s.e. %s); slope of nose %s (s.e. %s, n %d)." % (case, lab_, f6(fc[0]), f6(fc[1]), fc[2], f6(fx[0]), f6(fx[1]), f6(fn[0]), f6(fn[1]), fn[2]))
    L.append("")
with open(HERE_ / "checks/B9_t3check.md", "a", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
