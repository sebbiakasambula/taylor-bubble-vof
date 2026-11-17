import sys,glob,os
import numpy as np
from config import WORK, FLUENT
import tb_post
D=str(FLUENT/"T2_case")
a=np.loadtxt(D+"/volume-gas-rfile.out",skiprows=3); b=np.loadtxt(D+"/volume-gas-rfile_15001.out",skiprows=3)
rep=np.vstack([a[a[:,0]<15001],b])   # later run wins
assert np.all(np.diff(rep[:,0])>0)
def tstep(s): return np.interp(s,rep[:,0],rep[:,2])
print("report drift max %.2e %%"%(100*abs(rep[:,1]/rep[0,1]-1).max()))
rows=[]
for f in sorted([p for p in glob.glob(D+"/T2_Taylor_bubble-[0-9]*") if p.split("-")[-1].isdigit()],key=lambda p:int(p.split('-')[-1])):
    s=int(f.split('-')[-1]); d=np.loadtxt(f,delimiter=',',skiprows=1)
    x,y,ux,ur,V,g=d[:,1],np.round(d[:,2],9),d[:,4],d[:,5],d[:,8],d[:,9]; al=1-g
    n=tb_post.nose_position(x,y,al); tl=tb_post.tail_position(x,y,al)
    xi,ri=tb_post.interface_radius(x,y,al)
    w=(xi>tl+0.5)&(xi<n-1.5); fl=(0.5-ri[w]).mean() if w.sum()>=3 else np.nan
    sp=np.hypot(ux,ur); vg=(g*V).sum()*2*np.pi
    smear=(g*V)[g<0.5].sum()/(g*V).sum()
    rows.append((s,tstep(s),n,tl,fl,sp.max(),sp[g<0.01].max(),sp[g>0.99].max(),vg,(g*V*x).sum()/(g*V).sum(),smear))
    r=rows[-1]; print("%6d t %6.2f nose %6.3f tail %6.3f film %.4f umax %.2f liq %.2f gas %.2f cen %.3f smear %.3f"%(s,r[1],n,tl,fl,r[5],r[6],r[7],r[9],smear))
R=np.array(rows); print("cells",len(x),"radial",len(np.unique(y)),"axial",len(np.unique(np.round(x,9))))
for lo,hi in((5,20.0),(5,99),(10,20.0)):
    m=(R[:,1]>=lo)&(R[:,1]<=hi)
    print("fit %g-%g n=%d Fr nose %.4f cen %.4f film %.4f sd %.4f"%(lo,hi,m.sum(),np.polyfit(R[m,1],R[m,2],1)[0],np.polyfit(R[m,1],R[m,9],1)[0],np.nanmean(R[m,4]),np.nanstd(R[m,4])))
T=R[:,1]
for lo in range(0,21,5):
    k=(T>=lo)&(T<lo+5)
    if k.sum()>=3: print("  %2d-%2d Fr nose %.4f"%(lo,lo+5,np.polyfit(T[k],R[k,2],1)[0]))
np.save(WORK/"f_t2.npy",R)
O=np.load(WORK/"T2_prod.npy"); m=(O[:,0]>=5)&(O[:,0]<=20)
print("OpenFOAM T2 5-20: n=%d Fr nose %.4f cen %.4f film %.4f  peak median liq %.2f all %.2f"%(m.sum(),np.polyfit(O[m,0],O[m,1],1)[0],np.polyfit(O[m,0],O[m,3],1)[0],np.nanmean(O[m,4]),np.median(O[m,6]),np.median(O[m,5])))
m=(R[:,1]>=5)&(R[:,1]<=20)
print("Fluent T2 5-20 peak median liq %.2f all %.2f"%(np.median(R[m,6]),np.median(R[m,5])))
