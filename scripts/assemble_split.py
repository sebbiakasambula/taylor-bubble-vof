"""Assemble the fine-mesh split figure from the ParaView strips (render_split.py). EPS without transparency."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from config import PV, WORK  # noqa: E402
S = str(PV / "split_strips") + "/"
OUT = str(WORK / "fig_t2_split")  # copy the EPS to paper/figures
plt.rcParams.update({"font.size": 8.5, "ps.fonttype": 42})
fig = plt.figure(figsize=(6.3, 3.3))
g = fig.add_gridspec(3, 2, width_ratios=[1, 0.022], hspace=0.55, wspace=0.03)
panels = [("split_base_2.6", "(a) Momentum predictor off, $t$ = 2.6 s"),
          ("split_pred_2.6", "(b) Momentum predictor on, $t$ = 2.6 s"),
          ("split_pred_2.8", "(c) Momentum predictor on, $t$ = 2.8 s")]
for k, (n, t) in enumerate(panels):
    a = fig.add_subplot(g[k, 0])
    a.imshow(mpimg.imread(S + n + ".png"), extent=(1.5, 9.0, -0.55, 0.55), aspect="equal")
    a.set_title(t, fontsize=8.5, loc="left")
    a.set_ylabel("$r/D$")
    a.set_yticks([-0.5, 0, 0.5])
    if k == 2:
        a.set_xlabel("Height above the bottom of the pipe, $x/D$")
    else:
        a.set_xticklabels([])
c = fig.add_subplot(g[:, 1])
fig.colorbar(cm.ScalarMappable(LogNorm(1e-2, 1e2), "viridis"), cax=c, extend="both").set_label("Speed (m/s)")
fig.savefig(OUT + ".eps", bbox_inches="tight", dpi=300)
fig.savefig(S + "fig_t2_split_preview.png", bbox_inches="tight", dpi=110)
print("ok")
