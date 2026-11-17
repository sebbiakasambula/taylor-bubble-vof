import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import LogNorm; import matplotlib.cm as cm
from config import PV, WORK
R=str(WORK)+'/'; P=str(PV)+'/strips/'; OUT=R+'fig_s_fields.pdf'   # strips from paraview/render_flat.py
plt.rcParams.update({'font.size':8.5})
fig=plt.figure(figsize=(5.3,2.3))
g=fig.add_gridspec(1,3,width_ratios=[1,1,0.05],wspace=0.12)
for k,(n,t) in enumerate([('S_stretched','(a) S, stretched, $t$ = 0.99 s'),('S_uniform','(b) S, 96 by 12 uniform, $t$ = 1.00 s')]):
    a=fig.add_subplot(g[0,k]); a.imshow(mpimg.imread(P+n+'.png'),extent=(1.4,2.6,-0.5,0.5)); a.set_title(t,fontsize=8); a.set_xlabel('$x/D$'); a.set_xticks([1.5,2.0,2.5])
    if k==0: a.set_ylabel('$r/D$')
    else: a.set_yticklabels([])
c=fig.add_subplot(g[0,2]); fig.colorbar(cm.ScalarMappable(LogNorm(1e-3,10),'viridis'),cax=c).set_label('Speed (m/s)')
fig.savefig(OUT,bbox_inches='tight',dpi=300); fig.savefig(R+'fig_s_fields_preview.png',bbox_inches='tight',dpi=110)
