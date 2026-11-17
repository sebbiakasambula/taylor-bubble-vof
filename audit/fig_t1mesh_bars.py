"""Draft figure: change on the square mesh against the medium stretched mesh at T1 over 5 to 20 s, OpenFOAM and Fluent.
Values from the regenerated t1mesh_table.csv (audit/cmp), relative changes 100 (square / stretched - 1) from unrounded values.
Writes paper/figures/fig_t1mesh.eps (the paper uses it) and a preview png next to this script. No alpha (EPS).
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OF, FL, GRID = '#2a78d6', '#eb6834', '#e4e3dd'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'savefig.dpi': 300,
                     'legend.frameon': False, 'ps.fonttype': 42})
V = {}
for r in csv.DictReader(open(HERE / "cmp" / "t1mesh_table.csv")):
    try:
        V[(r["run"], r["statistic"])] = float(r["value"])
    except ValueError:
        pass
rows = [("Fr", "Fr"), ("Film thickness", "film_mean"), ("Gap 0.5 $D$ behind the nose", "gap05_mean"),
        ("Bubble length", "L_b_mean"), ("Wall layer below the tail", "skirt_b_mean"),
        ("Gas outside the bubble", "out_b_mean"), ("Largest speed (median)", "pk_all_median"),
        ("Largest liquid speed (median)", "pk_liq_median"), ("Mean gas speed in the bubble", "u_core_median")]
ch = {c: [100 * (V[(c + " square", k)] / V[(c + " stretched", k)] - 1) for _, k in rows] for c in ("OF", "F")}
fig, ax = plt.subplots(figsize=(6.3, 3.2))
y = np.arange(len(rows))[::-1]
h = 0.36
ax.barh(y + h / 2, ch["OF"], height=h, color=OF, label="OpenFOAM", zorder=3)
ax.barh(y - h / 2, ch["F"], height=h, color=FL, label="Fluent", zorder=3)
for yy, a, b in zip(y, ch["OF"], ch["F"]):
    for v, dy in ((a, h / 2), (b, -h / 2)):
        ax.text(v + (0.6 if v >= 0 else -0.6), yy + dy, "%+.1f" % v, va="center", ha="left" if v >= 0 else "right", fontsize=7, color="#333")
ax.axvline(0, color="#555", lw=0.8, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in rows], fontsize=8)
ax.set_xlabel("Change on the square mesh against the stretched mesh (%)")
lo, hi = min(min(ch["OF"]), min(ch["F"])), max(max(ch["OF"]), max(ch["F"]))
ax.set_xlim(lo - 6, hi + 6)
ax.grid(axis="y", visible=False)
ax.legend(loc="lower right", fontsize=7.5)
fig.tight_layout()
import sys  # noqa: E402
sys.path.insert(0, str(HERE.parent / "scripts"))  # config.py
from config import WORK  # noqa: E402
fig.savefig(WORK / "fig_t1mesh.eps")  # copy it to paper/figures
fig.savefig(HERE / "fig_t1mesh_bars_preview.png", dpi=110)
for (lab, _), a, b in zip(rows, ch["OF"], ch["F"]):
    print("%-32s OF %+8.4f%%  F %+8.4f%%" % (lab, a, b))
