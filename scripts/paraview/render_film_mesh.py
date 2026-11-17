"""B10 render: the front-wedge-face cells of two meshes in a window behind the T1 nose, coloured by alpha.water, with all
cell edges and the alpha = 0.5 line. Input: film_mesh_<case>.vtk and film_ctr_<case>.vtk from prep_film_mesh.py.
Output: film_<case>.png, exactly the window x - x_nose in [-2.56, -2.44], r in [0.42, 0.505] (parallel projection, 1:1).
The wall line, the dashed correlation line and the axes are added by fig_film_mesh.py (matplotlib)."""
import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(globals().get('__file__', sys.argv[0])).resolve().parent.parent))   # scripts/ holds config.py
from config import PV as _PV
from paraview.simple import *

D = str(_PV) + '/'
XW = (-2.56, -2.44); RW = (0.42, 0.505)
PX_PER_UNIT = 20000                                  # 0.12 -> 2400 px, 0.085 -> 1700 px
RES = [int(round((XW[1] - XW[0]) * PX_PER_UNIT)), int(round((RW[1] - RW[0]) * PX_PER_UNIT))]
GAS, LIQ = [0.90, 0.95, 1.00], [0.76, 0.76, 0.76]    # light blue gas, light grey liquid
EDGE = [0.15, 0.15, 0.15]
EDGE_PX = 3.0                                        # about 0.25 pt in the final 2.8 in wide panel
IFACE = [0.835, 0.369, 0.0]                          # vermilion
TUBE_R = 0.0003                                      # data units; 2 R = 0.0006 = about 1 pt in the final panel (23 in per unit of D)


def panel(case):
    r = LegacyVTKReader(registrationName='mesh_' + case, FileNames=[D + 'film_mesh_' + case + '.vtk'])
    v = CreateRenderView(); v.OrientationAxesVisibility = 0; v.Background = [1, 1, 1]
    try: v.UseColorPaletteForBackground = 0
    except Exception: pass
    v.ViewSize = RES
    d = Show(r, v); d.SetRepresentationType('Surface With Edges')
    ColorBy(d, ('CELLS', 'alpha.water'))
    lut = GetColorTransferFunction('alpha.water')
    lut.RGBPoints = [0.0] + GAS + [1.0] + LIQ
    lut.ColorSpace = 'RGB'; lut.AutomaticRescaleRangeMode = 'Never'; lut.RescaleTransferFunction(0.0, 1.0)
    d.LookupTable = lut; d.SetScalarBarVisibility(v, False)
    d.EdgeColor = EDGE; d.LineWidth = EDGE_PX; d.Ambient = 1.0; d.Diffuse = 0.0; d.Interpolation = 'Flat'
    # alpha = 0.5 through the cell-centre values, drawn as a flat-coloured tube so the width is not limited by the GL line width
    g = LegacyVTKReader(registrationName='ctr_' + case, FileNames=[D + 'film_ctr_' + case + '.vtk'])
    c = Contour(Input=g, ContourBy=['POINTS', 'alpha.water'], Isosurfaces=[0.5])
    t = Tube(Input=c); t.Radius = TUBE_R; t.NumberofSides = 20; t.Capping = 1
    td = Show(t, v); ColorBy(td, None)
    td.AmbientColor = IFACE; td.DiffuseColor = IFACE; td.Ambient = 1.0; td.Diffuse = 0.0; td.Specular = 0.0
    Render(v)
    cx, cy = 0.5 * (XW[0] + XW[1]), 0.5 * (RW[0] + RW[1])
    v.CameraPosition = [cx, cy, 10]; v.CameraFocalPoint = [cx, cy, 0]; v.CameraViewUp = [0, 1, 0]
    v.CameraParallelProjection = 1; v.CameraParallelScale = 0.5 * (RW[1] - RW[0])
    Render(v)
    SaveScreenshot(D + 'film_' + case + '.png', v, ImageResolution=RES, TransparentBackground=0)
    print(case, 'saved', RES, 'cells', r.GetDataInformation().GetNumberOfCells(), 'contour pts', c.GetDataInformation().GetNumberOfPoints())


for case in ['V4_medium', 'VA_square']:
    panel(case)
print('done')
