"""Fig. 1 (GT phantom and VOIs), Fig. 2 (blurred images) and Supplementary Fig. S1 (tumour TAC fits).
Outputs: results/figs/Fig1.png, Fig2.png, FigS1.png
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results')
import numpy as np, pandas as pd, json, glob, pydicom, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from matplotlib.colors import PowerNorm
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica','DejaVu Sans'],'font.size':8,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
                     'xtick.major.size':2.5,'ytick.major.size':2.5,'legend.frameon':False})
MM=1/25.4; W2=174*MM; W1=84*MM
C={'Tumour':'#D9731E','Spleen':'#1F8A94','Right kidney':'#2F4F9E','Left kidney':'#7C99D4','Liver':'#8A7650'}
CONT={1:'#C9B27C',2:'#F08A24',3:'#3FB6C1',4:'#5B8CFF',5:'#A9C4FF'}
MAP={'Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney','Liver':'Liver'}
GT={'Tumour':28.1,'Spleen':5.39,'Right kidney':3.6,'Left kidney':3.61,'Liver':1.70}
FW=np.array([0,5,7.5,10,12.5,15.])
S=pd.read_csv(f'{ROOT}/results/summary_recommended_models.csv'); S['V']=S.VOI.map(MAP)
AD={v:S[S.V==v].sort_values('FWHM_mm').AD_Gy.values for v in C}; SD={v:S[S.V==v].sort_values('FWHM_mm').AD_SD_Gy.values for v in C}
seg=np.load(f'{ROOT}/derived/seg.npy')
VX=4.419644; x0,y0,z0=-279.947394,-279.047394,284.897394; ox=oy=-250.0; sxy=0.976562; oz=-167.0; sz=2.5
def load(t,fw):
    p=f'{ROOT}/input/truth_{t}_reg.dcm' if fw==0 else f'{ROOT}/results/blurred_dicom/FWHM_{fw:04.1f}mm/truth_{t}_reg_G{fw:g}mm.dcm'
    d=pydicom.dcmread(p); return d.pixel_array*float(d.RescaleSlope)
def axial(ax,img,k,xl,yl,vmax,gamma=0.8,lw=0.8):
    zc=z0-k*VX; ext=(x0-VX/2,x0+127.5*VX,y0+127.5*VX,y0-VX/2)
    ax.imshow(img[k],cmap='gray',norm=PowerNorm(gamma,0,vmax),extent=ext,interpolation='nearest')
    j0=int(np.floor((zc-VX/2-oz)/sz)); j1=int(np.ceil((zc+VX/2-oz)/sz)); xs=ox+np.arange(512)*sxy; ys=oy+np.arange(512)*sxy
    for l in (1,2,3,4,5):
        m=np.mean([(seg[j]==l) for j in range(max(j0,0),min(j1,seg.shape[0]))],axis=0)
        if m.max()>0: ax.contour(xs,ys,m,levels=[0.5],colors=[CONT[l]],linewidths=lw)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)
def coronal(ax,img,r,xl,zl,vmax,gamma=0.8,lw=0.8):
    yc=y0+r*VX; ext=(x0-VX/2,x0+127.5*VX,z0-127.5*VX,z0+VX/2)
    ax.imshow(img[:,r,:],cmap='gray',norm=PowerNorm(gamma,0,vmax),extent=ext,interpolation='nearest',aspect='equal')
    jr=int(round((yc-oy)/sxy)); xs=ox+np.arange(512)*sxy; zs=oz+np.arange(seg.shape[0])*sz
    for l in (1,2,3,4,5):
        m=np.mean([(seg[:,j,:]==l) for j in range(jr-2,jr+3)],axis=0)
        if m.max()>0: ax.contour(xs,zs,m,levels=[0.5],colors=[CONT[l]],linewidths=lw)
    ax.set_xlim(*xl); ax.set_ylim(*zl); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)
def lab(ax,s,x=0.02,y=0.98,color='white'): ax.text(x,y,s,transform=ax.transAxes,fontsize=9,fontweight='bold',color=color,va='top',ha='left')
legh=[Line2D([0],[0],color=CONT[i],lw=1.5) for i in (1,2,3,4,5)]; legl=['Liver','Tumour','Spleen','Right kidney','Left kidney']

from resdos.manifest import record
# ---------------- Fig 1
img=load('24h',0.0); vmax=img.max()*0.4
fig=plt.figure(figsize=(W2,78*MM),facecolor='white')
a=fig.add_axes([0.0,0.1,0.36,0.9]); coronal(a,img,62,(-140,140),(-150,140),vmax); lab(a,'(a)')
b=fig.add_axes([0.37,0.55,0.315,0.45]); axial(b,img,55,(-122,113),(72,-78),vmax); lab(b,'(b)')
c=fig.add_axes([0.37,0.1,0.315,0.43]); axial(c,img,76,(-122,113),(88,-62),vmax*0.75); lab(c,'(c)')
d=fig.add_axes([0.7,0.1,0.3,0.9]); d.axis('off')
txt=['Ground-truth activity map','(24 h, 128$^3$ voxels of 4.42 mm)','','CT-based VOIs:']
d.text(0.02,0.95,'\n'.join(txt),fontsize=7.5,va='top')
d.legend(legh,legl,loc='upper left',bbox_to_anchor=(0.0,0.80),fontsize=7.5,handlelength=1.5)
d.text(0.02,0.36,'Partial-volume representation:\nboundary voxels hold the\nvolume-weighted mean of\ntheir tissues',fontsize=7,color='#4A6272',va='top')
a.plot([-120,-70],[-140,-140],color='white',lw=1.5); a.text(-95,-136,'50 mm',color='white',fontsize=6.5,ha='center',va='bottom')
fig.savefig('figs/Fig1.png',dpi=600,facecolor='white'); plt.close(fig)
# ---------------- Fig 2
fws=[0.0,5.0,7.5,10.0,12.5,15.0]; ref=load('24h',0.0); vmax=ref.max()*0.4
H2=52*MM
fig=plt.figure(figsize=(W2,H2),facecolor='white')
lm=0.055; cw=(1-lm-0.005-5*0.006)/6; ih=cw*W2*150/235/H2
for j,fw in enumerate(fws):
    im=load('24h',fw); x=lm+j*(cw+0.006)
    fig.text(x+cw/2,0.935,f'{fw:g} mm',ha='center',fontsize=8,fontweight='bold')
    a1=fig.add_axes([x,0.90-ih,cw,ih]); axial(a1,im,55,(-122,113),(72,-78),vmax,lw=0.6)
    a2=fig.add_axes([x,0.90-2*ih-0.015,cw,ih]); axial(a2,im,76,(-122,113),(88,-62),vmax*0.75,lw=0.6)
fig.text(0.005,0.90-ih/2,'Tumour\nlevel',fontsize=7,va='center'); fig.text(0.005,0.90-1.5*ih-0.015,'Kidney\nlevel',fontsize=7,va='center')
fig.legend(legh,legl,loc='lower center',ncol=5,fontsize=7,bbox_to_anchor=(0.53,-0.02),handlelength=1.5)
fig.savefig('figs/Fig2.png',dpi=600,facecolor='white'); plt.close(fig)
# ---------------- Fig 3: tumour TAC fits
T=np.array([1,4,24,40,72,144.]); lp=4.345e-3; lbc=41.59
f4b=lambda t,A1,l1,A2,l2: A1*np.exp(-(l1+lp)*t)-A2*np.exp(-(l2+lp)*t)-(A1-A2)*np.exp(-(lbc+lp)*t)
D=pd.concat([pd.read_csv(f,sep='\t').set_index('Segment')['Mean'] for f in sorted(glob.glob(f'{ROOT}/derived/opendose3d/fw00/*.tsv'),key=lambda s:float(s.split('_')[-1].replace('HR.tsv','')))],axis=1)
y=D.loc['Tumour'].values
P=[float(v) for v in pd.read_csv(f'{ROOT}/results/tac_fit_results_0mm.csv').set_index(['VOI','model']).loc[('Tumour','f4b'),'params'].split(';')]
fig,ax=plt.subplots(1,2,figsize=(W2,62*MM))
for a_,(xm,log) in zip(ax,[(160,False),(600,True)]):
    tt=np.linspace(0,xm,1201)
    a_.plot(tt,f4b(tt,*P),color=C['Tumour'],lw=1.3,label='f4b, converged (T$_{eff}$ 95 h; 24.1 Gy)')
    a_.plot(tt,f4b(tt,439.9,3.735e-4,287.5,2.644e-21),'--',color='#4A6272',lw=1.1,label='f4b, $\\lambda_2\\rightarrow$0 (T$_{eff}$ 126 h; 27.1 Gy)')
    a_.plot(T,y,'o',color='k',ms=3.5,label='VOI mean dose rate')
    a_.set_xlim(0,xm); a_.set_xlabel('Time after administration (h)'); a_.set_ylabel('Tumour dose rate (mGy h$^{-1}$)')
    for s_ in ['top','right']: a_.spines[s_].set_visible(False)
    if log: a_.set_yscale('log'); a_.set_ylim(3,250)
    else: a_.set_ylim(0,170)
ax[0].legend(fontsize=6.5,loc='lower left'); lab(ax[0],'(a)',x=-0.17,y=1.05,color='k'); lab(ax[1],'(b)',x=-0.17,y=1.05,color='k')
fig.tight_layout(); fig.savefig('figs/FigS1.png',dpi=600,facecolor='white'); plt.close(fig)

record('figs/Fig1.png','05_figures_tables/f01_fig1_fig2_figS1.py',3); record('figs/Fig2.png','05_figures_tables/f01_fig1_fig2_figS1.py',1); record('figs/FigS1.png','05_figures_tables/f01_fig1_fig2_figS1.py',2)
print('Fig 1, Fig 2, Fig S1 done')
