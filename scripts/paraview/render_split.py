"""pvbatch: render the three strips of the fine-mesh split figure (unlit, log speed, interface at alpha 0.5 in white),
same style as rerender/render_flat.py. Window x 1.5 to 9.0 D, r -0.55 to 0.55 D."""
import os

from paraview.simple import *  # noqa: F401,F403

import sys  # noqa: E402
from pathlib import Path as _Path  # noqa: E402
sys.path.insert(0, str(_Path(globals().get('__file__', sys.argv[0])).resolve().parent.parent))  # scripts/ holds config.py
from config import PV as _PV  # noqa: E402
D = str(_PV) + "/"
O = D + "split_strips/"
os.makedirs(O, exist_ok=True)


def panel(name, rng, centre, half_h, res):
    r = LegacyVTKReader(FileNames=[D + name + ".vtk"])
    v = CreateRenderView()
    v.OrientationAxesVisibility = 0
    v.Background = [1, 1, 1]
    try:
        v.UseColorPaletteForBackground = 0
    except Exception:
        pass
    v.ViewSize = res
    d = Show(r, v)
    ColorBy(d, ("POINTS", "speed"))
    d.Ambient = 1.0
    d.Diffuse = 0.0
    lut = GetColorTransferFunction("speed_" + name.replace(".", "_"))
    lut.ApplyPreset("Viridis", True)
    lut.UseLogScale = 1
    lut.RescaleTransferFunction(*rng)
    lut.AutomaticRescaleRangeMode = "Never"
    d.LookupTable = lut
    d.SetScalarBarVisibility(v, False)
    c = Contour(Input=r, ContourBy=["POINTS", "alpha"], Isosurfaces=[0.5])
    tr = Transform(Input=c)
    tr.Transform.Translate = [0.0, 0.0, 0.05]
    cd = Show(tr, v)
    cd.ColorArrayName = ["POINTS", ""]
    cd.Ambient = 1.0
    cd.Diffuse = 0.0
    cd.AmbientColor = [1, 1, 1]
    cd.DiffuseColor = [1, 1, 1]
    cd.LineWidth = 2.5
    Render(v)
    v.CameraPosition = [centre[0], centre[1], 10]
    v.CameraFocalPoint = [centre[0], centre[1], 0]
    v.CameraViewUp = [0, 1, 0]
    v.CameraParallelProjection = 1
    v.CameraParallelScale = half_h
    Render(v)
    SaveScreenshot(O + name + ".png", v, ImageResolution=res, TransparentBackground=0)


# the first view of a session renders without the contour line, so render one throwaway panel first
for n in ["split_pred_2.8", "split_base_2.6", "split_pred_2.6", "split_pred_2.8"]:
    panel(n, (1e-2, 1e2), (5.25, 0.0), 0.55, [2400, 352])
