import sys, os
import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(globals().get('__file__', sys.argv[0])).resolve().parent.parent))   # scripts/ holds config.py
from config import PV as _PV
from paraview.simple import *
D=str(_PV)+'/'
O=str(_PV)+'/strips/'   # strips go here; the VTK inputs are read from the same pv folder
os.makedirs(O,exist_ok=True)
def panel(name,rng,log,centre,up,half_h,res):
    r=LegacyVTKReader(FileNames=[D+name+'.vtk'])
    v=CreateRenderView(); v.OrientationAxesVisibility=0; v.Background=[1,1,1]
    try: v.UseColorPaletteForBackground=0
    except: pass
    v.ViewSize=res
    d=Show(r,v); ColorBy(d,('POINTS','speed'))
    d.Ambient=1.0; d.Diffuse=0.0     # unlit, so the strip colours equal the colour bar (default lighting darkens them by about 22 %)
    lut=GetColorTransferFunction('speed_'+name); lut.ApplyPreset('Viridis',True)
    if log: lut.UseLogScale=1
    lut.RescaleTransferFunction(*rng); lut.AutomaticRescaleRangeMode='Never'; d.LookupTable=lut; d.SetScalarBarVisibility(v,False)
    c=Contour(Input=r,ContourBy=['POINTS','alpha'],Isosurfaces=[0.5])
    cd=Show(c,v); ColorBy(cd,None); cd.AmbientColor=[1,1,1]; cd.DiffuseColor=[1,1,1]; cd.LineWidth=2.5
    Render(v); v.CameraPosition=[centre[0],centre[1],10]; v.CameraFocalPoint=[centre[0],centre[1],0]; v.CameraViewUp=up
    v.CameraParallelProjection=1; v.CameraParallelScale=half_h
    Render(v); SaveScreenshot(O+name+'.png',v,ImageResolution=res,TransparentBackground=0)
# static bubble, x horizontal, window x 1.4 to 2.6, r -0.5 to 0.5
for n in ['S_stretched','S_uniform']: panel(n,(1e-3,10),True,(2.0,0.0),[0,1,0],0.5,[1200,1000])
# T3, x up, nose at 0, window x -5 to 0.5
for n in ['T3_OF','T3_F']: panel(n,(0,2),False,(-2.3,0.0),[0,1,0],0.55,[2200,400])
print('done')
