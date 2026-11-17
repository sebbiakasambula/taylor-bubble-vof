import numpy as np
from pathlib import Path
from config import RUNS, FLUENT, PV, WORK
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])
RUN=RUNS; FD=FLUENT
OUT=PV
def write(name,x,r,sp,al,shift=0.0):
    # structured grid of cell centres, mirrored about the axis
    X=np.unique(np.round(x,6)); R=np.unique(np.round(r,6)); nx,nr=len(X),len(R)
    assert nx*nr==len(x),(name,nx,nr,len(x))
    ix=np.searchsorted(X,np.round(x,6)); ir=np.searchsorted(R,np.round(r,6))
    S=np.zeros((nr,nx)); A=np.zeros((nr,nx)); S[ir,ix]=sp; A[ir,ix]=al
    Rm=np.r_[-R[::-1],R]; S=np.vstack([S[::-1],S]); A=np.vstack([A[::-1],A])
    XX,RR=np.meshgrid(X-shift,Rm)
    with open(OUT/(name+'.vtk'),'w') as f:
        f.write('# vtk DataFile Version 3.0\n%s\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS %d %d 1\nPOINTS %d float\n'%(name,nx,2*nr,nx*2*nr))
        np.savetxt(f,np.c_[XX.ravel(),RR.ravel(),np.zeros(XX.size)],fmt='%.6g')
        f.write('POINT_DATA %d\nSCALARS speed float 1\nLOOKUP_TABLE default\n'%(XX.size)); np.savetxt(f,S.ravel(),fmt='%.6g')
        f.write('SCALARS alpha float 1\nLOOKUP_TABLE default\n'); np.savetxt(f,A.ravel(),fmt='%.6g')
    print(name,nx,nr,'max speed %.3f'%sp.max(),'liquid max %.3f'%sp[al>0.99].max())
def of(case,t,name,nose=False):
    C=vec(RUN/case/'0/C'); a=scal(RUN/case/t/'alpha.water'); U=vec(RUN/case/t/'U')
    x,y=C[:,0],C[:,1]; sh=tb_post.nose_position(x,np.round(y,9),a) if nose else 0.0
    write(name,x,y,np.linalg.norm(U,axis=1),a,sh); return sh
def fl(folder,fn,name):
    d=np.loadtxt(FD/folder/fn,delimiter=',',skiprows=1); x,y=d[:,1],d[:,2]; a=1-d[:,9]
    sh=tb_post.nose_position(x,np.round(y,9),a)
    write(name,x,y,np.hypot(d[:,4],d[:,5]),a,sh); return sh
of('V1a_static','0.988705','S_stretched'); of('D_coarseAR1','0.999955','S_uniform')
print('nose OF',of('T3_prod','20.305026','T3_OF',True)); print('nose F',fl('T3_Case','T3_Taylor_bubble-78000','T3_F'))
