"""Set 19. Regime numbers per case (17 Nov): Re, Mo and the viscous diffusion number of the wall cell.

Re = Fr Nf with Fr of Viana et al., Mo = Eo^3 / Nf^4 (equal to g mu_l^4 drho / (rho_l^2 sigma^3) with these
definitions), and nu dt / dr^2 for the liquid and the gas at the wall cell of the medium mesh, with dt the time step
cap of the medium run (maxDeltaT in its controlDict, held in 87% of the T1 steps and by Fluent throughout).
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from tb_params import fluid, morton  # noqa: E402
from tb_references import fr_viana  # noqa: E402

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402


def radii(case):
    s = (RUNS / case / "constant/polyMesh/points").read_text()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    P = np.fromstring(body, sep=" ", count=3 * n).reshape(n, 3)
    return np.unique(np.round(np.hypot(P[:, 1], P[:, 2]), 9))


def max_dt(case):
    return float(re.search(r"maxDeltaT\s+([0-9.eE+-]+);", (RUNS / case / "system/controlDict").read_text())[1])


r = radii("T2_prod")
dr_wall, dr_axis = r[-1] - r[-2], r[1] - r[0]
lines = ["# Set 19. Regime numbers", "", "Medium mesh: wall cell %.7f D, axis cell %.7f D." % (dr_wall, dr_axis), "",
         "| case | Eo | Nf | Fr Viana | Re | Mo (Eo^3/Nf^4) | Mo (tb_params) | dt cap | nu_l | nu_g | nu_l dt/dr^2 wall | nu_g dt/dr^2 wall | nu_g dt/dr^2 axis |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for case, run, Eo, Nf in (("T1", "V4_medium", 300, 1000), ("T2", "T2_prod", 100, 50), ("T3", "T3_prod", 20, 500)):
    f = fluid(Eo, Nf, 1000)
    fr = fr_viana(Eo, Nf)
    dt = max_dt(run)
    lines.append("| %s | %d | %d | %.4f | %.2f | %.4g | %.4g | %.6g | %.6g | %.6g | %.4f | %.4f | %.4f |" % (
        case, Eo, Nf, fr, fr * Nf, Eo ** 3 / Nf ** 4, morton(f), dt, f["nu_l"], f["nu_g"],
        f["nu_l"] * dt / dr_wall ** 2, f["nu_g"] * dt / dr_wall ** 2, f["nu_g"] * dt / dr_axis ** 2))
f = fluid(10, 2000, 1000)
lines.append("| S | 10 | 2000 | | | %.4g | %.4g | | %.6g | %.6g | | | |" % (10 ** 3 / 2000 ** 4, morton(f), f["nu_l"], f["nu_g"]))
lines.append("")
lines.append("Half the time step at T2 (T2_A9_halfdt, maxDeltaT %.7g): liquid %.4f, gas %.4f at the wall." % (
    max_dt("T2_A9_halfdt"), fluid(100, 50, 1000)["nu_l"] * max_dt("T2_A9_halfdt") / dr_wall ** 2,
    fluid(100, 50, 1000)["nu_g"] * max_dt("T2_A9_halfdt") / dr_wall ** 2))
(HERE / "set19_regime_out.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
