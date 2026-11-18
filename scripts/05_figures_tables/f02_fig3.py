"""Fig. 3: recovery coefficients (AD and TIA) versus FWHM.
Outputs: results/figs/Fig3.png
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
tia=json.load(open('tia_rc.json')); sph={'validation':json.load(open('sphere_model/tumour_validation.json'))}
KEY={'Tumour':'Tumour','Spleen':'spleen','Right kidney':'Right Kidney','Left kidney':'Left Kidney','Liver':'Liver'}
st=json.load(open('stats.json'))
# Fig 4: RC (AD and TIA) + sphere model
tia=json.load(open('tia_rc.json')); sph={'validation':json.load(open('sphere_model/tumour_validation.json'))}
fig,ax=plt.subplots(figsize=(120*MM,80*MM))
KEY={'Tumour':'Tumour','Spleen':'spleen','Right kidney':'Right Kidney','Left kidney':'Left Kidney','Liver':'Liver'}
for k in C:
    rc=AD[k]/AD[k][0]; cv=SD[k]/AD[k]; err=np.r_[0,rc[1:]*np.sqrt(cv[1:]**2+cv[0]**2)]
    ax.errorbar(FW,rc,yerr=err,color=C[k],marker='o',ms=3.2,lw=1.2,capsize=1.8,elinewidth=0.6,label=f'{k} (AD)')
    rt=[tia[KEY[k]][str(f)]['RC_TIA'] for f in FW]
    ax.plot(FW,rt,'s',mfc='white',mec=C[k],ms=3.6,mew=0.8)
sv=sph['validation']; ax.plot(sv['FWHM'],sv['RC'],':',color='k',lw=1.0,label='Sphere model (15.9 mL, contrast 16.8)')
ax.plot([],[],'s',mfc='white',mec='k',ms=3.6,mew=0.8,label='TIA (open squares)')
ax.set_xlabel('3D Gaussian FWHM (mm)'); ax.set_ylabel('RC = AD(FWHM) / AD(0 mm)'); ax.set_xticks(FW); ax.set_xticklabels(['0','5','7.5','10','12.5','15'])
ax.set_ylim(0.55,1.05); ax.grid(axis='y',color='#E3E8EB',lw=0.5); ax.set_axisbelow(True)
for s_ in ['top','right']: ax.spines[s_].set_visible(False)
ax.legend(fontsize=6.3,loc='lower left',ncol=1)
fig.tight_layout(); fig.savefig('figs/Fig3.png',dpi=600,facecolor='white'); plt.close(fig)

from resdos.manifest import record
record('figs/Fig3.png','05_figures_tables/f02_fig3.py',1)
print('Fig 3 done')
