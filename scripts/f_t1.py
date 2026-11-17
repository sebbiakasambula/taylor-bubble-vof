import sys,glob,os,re
import numpy as np
from config import WORK, FLUENT
import tb_post
D=str(FLUENT/"T1_Case")
rep=np.loadtxt(D+"/volume-gas-rfile.out",skiprows=3)
step2t=dict(zip(rep[:,0].astype(int),rep[:,2]))
def tstep(s): return step2t[s] if s in step2t else rep[-1,2]+(s-int(rep[-1,0]))*9.1e-4
rows=[]
for f in sorted(glob.glob(D+"/s_snap-*"),key=lambda p:int(p.split('-')[-1])):
    s=int(f.split('-')[-1]); d=np.loadtxt(f,delimiter=',',skiprows=1)
    x,y,p,ux,ur,V,g=d[:,1],d[:,2],d[:,3],d[:,4],d[:,5],d[:,8],d[:,9]
    y=np.round(y,9); a=1-g
    n=tb_post.nose_position(x,y,a); tl=tb_post.tail_position(x,y,a)
    xi,ri=tb_post.interface_radius(x,y,a)
    w=(xi>tl+0.5)&(xi<n-2.0); fl=(0.5-ri[w]).mean() if w.sum()>=3 else np.nan
    sp=np.hypot(ux,ur); vg=(g*V).sum()*2*np.pi
    rows.append((s,tstep(s),n,tl,fl,sp.max(),sp[g<0.01].max(),sp[g>0.99].max(),vg,(g*V*x).sum()/(g*V).sum()))
    print("%6d t %7.3f nose %7.3f tail %7.3f film %s umax %.2f liq %.2f gas %.2f Vgas2pi %.5f cen %.3f"%(rows[-1][0],rows[-1][1],n,tl,"%.4f"%fl,rows[-1][5],rows[-1][6],rows[-1][7],vg,rows[-1][9]))
R=np.array(rows); print("cells",len(x),"radial",len(np.unique(y)),"axial",len(np.unique(np.round(x,9))))
for lo,hi in((5,20.01),):
    m=(R[:,1]>=lo)&(R[:,1]<=hi)
    print("fit %g-%g n=%d Fr %.4f  film mean %.4f n %d sd %.4f"%(lo,hi,m.sum(),np.polyfit(R[m,1],R[m,2],1)[0],np.nanmean(R[m,4]),np.isfinite(R[m,4]).sum(),np.nanstd(R[m,4])))
np.save(os.path.join(WORK,"f_t1.npy"),R)
