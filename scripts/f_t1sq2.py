# Fluent T1 square: original run to step 6000, restart from the step-6000 data file after it (10 Nov).
# Every snapshot time is read from a report file, no interpolation.
import sys,glob
import numpy as np
from config import WORK, FLUENT
import tb_post
B=str(FLUENT)+"/"
o=np.loadtxt(B+"T1_square/volume-gas-rfile.out",skiprows=3); o=o[o[:,0]<=6734]
n=np.loadtxt(B+"T1_square_w_completion/volume-gas-rfile_6001.out",skiprows=3)
files=[(int(f.split('-')[-1]),f,o) for f in glob.glob(B+"T1_square/s_snap-*") if int(f.split('-')[-1])<=6000]
files+=[(int(f.split('-')[-1]),f,n) for f in glob.glob(B+"T1_square_w_completion/s_snap_completion-*")]
rows=[]
for s,f,rep in sorted(files):
    k=np.where(rep[:,0]==s)[0]; assert len(k)==1,(s,f); t=rep[k[0],2]
    d=np.loadtxt(f,delimiter=',',skiprows=1)
    x,y,ux,ur,V,g=d[:,1],np.round(d[:,2],9),d[:,4],d[:,5],d[:,8],d[:,9]; al=1-g
    nz=tb_post.nose_position(x,y,al); tl=tb_post.tail_position(x,y,al)
    xi,ri=tb_post.interface_radius(x,y,al)
    w=(xi>tl+0.5)&(xi<nz-2.0); fl=(0.5-ri[w]).mean() if w.sum()>=3 else np.nan
    sp=np.hypot(ux,ur); smear=(g*V)[g<0.5].sum()/(g*V).sum()
    rows.append((s,t,nz,tl,fl,sp.max(),sp[g<0.01].max(),sp[g>0.99].max(),(g*V).sum()*2*np.pi,(g*V*x).sum()/(g*V).sum(),smear,0))
    r=rows[-1]; print("%5d %s t %7.4f nose %6.3f tail %6.3f film %.4f umax %5.2f liq %.2f cen %.3f smear %.3f"%(s,'orig' if rep is o else 'rest',t,nz,tl,fl,r[5],r[6],r[9],smear))
R=np.array(rows); np.save(WORK/"f_t1sq.npy",R)
print("interval speeds",np.round(np.diff(R[:,2])/np.diff(R[:,1]),4))
for lo,hi in((5,20.01),(5,15.61),(10,20.01),(5,21)):
    m=(R[:,1]>=lo)&(R[:,1]<=hi); p,cov=np.polyfit(R[m,1],R[m,2],1,cov=True)
    print("fit %g-%g n=%d Fr nose %.4f (se %.4f) cen %.4f film %.4f (n %d)"%(lo,hi,m.sum(),p[0],np.sqrt(cov[0,0]),np.polyfit(R[m,1],R[m,9],1)[0],np.nanmean(R[m,4]),np.isfinite(R[m,4]).sum()))
m=(R[:,1]>=5)&(R[:,1]<=20.01); print("median peak 5-20: all %.2f liq %.2f"%(np.median(R[m,5]),np.median(R[m,6])))
print("gas volume drift %.2e %%"%(100*abs(np.r_[o[:,1],n[:,1]]/o[0,1]-1).max()))
