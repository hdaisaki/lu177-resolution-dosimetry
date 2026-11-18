"""Fig. 4 (RE versus FWHM) and Fig. 6 (measured SPECT: profile, local model held-out, proportional AD estimates).
Outputs: results/figs/Fig4.png, Fig6.png
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results')
import numpy as np, pandas as pd, json, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy import ndimage as ndi
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica','DejaVu Sans'],'font.size':8,'axes.linewidth':0.6,'legend.frameon':False,'xtick.major.size':2.5,'ytick.major.size':2.5})
MM=1/25.4
C={'Tumour':'#D9731E','spleen':'#1F8A94','Right Kidney':'#2F4F9E','Left Kidney':'#7C99D4','Liver':'#8A7650'}
NM={'Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney','Liver':'Liver'}
cats=['Tumour','spleen','Right Kidney','Left Kidney','Liver']
GT={'Tumour':28.1,'spleen':5.39,'Right Kidney':3.6,'Left Kidney':3.61,'Liver':1.70}
S=pd.read_csv(f'{ROOT}/results/summary_recommended_models.csv')
FW=np.array([0,5,7.5,10,12.5,15.])
def clean(a):
    for s_ in ['top','right']: a.spines[s_].set_visible(False)
# ---------------- Fig 5: RE vs FWHM (no regression)
fig,ax=plt.subplots(figsize=(120*MM,78*MM))
for v in cats:
    d=S[S.VOI==v].sort_values('FWHM_mm'); re=d.RE_vs_GT_pct.values; e=d.AD_SD_Gy.values/GT[v]*100
    ax.errorbar(FW,re,yerr=e,color=C[v],marker='o',ms=3.2,lw=1.1,elinewidth=0.6,capsize=1.8,label=NM[v])
ax.set_xlabel('3D Gaussian FWHM (mm)'); ax.set_ylabel('Relative error vs ground truth (%)'); ax.set_xticks(FW); ax.set_xticklabels(['0','5','7.5','10','12.5','15'])
ax.set_ylim(-50,0); ax.grid(axis='y',color='#E3E8EB',lw=0.5); ax.set_axisbelow(True); clean(ax); ax.legend(fontsize=6.5,loc='lower left')
fig.tight_layout(); fig.savefig('figs/Fig4.png',dpi=600,facecolor='white'); plt.close(fig)
# (Fig. 5 of the manuscript is produced by sphere_model.py and must not be regenerated here)
# ---------------- Fig 8
Sp=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy'); k=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['k_abs']
MG={t:json.load(open(f'{ROOT}/results/spect/model_grid_{t}.json')) for t in TS}; LF={t:json.load(open(f'{ROOT}/results/spect/local_fit_{t}.json')) for t in TS}
FG=np.array(sorted(float(f) for f in MG['1h']))
tab={'Tumour':[1916.5,1809.2,1864.5,1713.2,1396.6,825.9],'spleen':[374.5,428.0,500.4,433.7,314.4,151.7],'Right Kidney':[494.9,466.2,403.5,326.1,214.4,82.7],
     'Left Kidney':[490.0,456.5,404.2,325.1,212.4,82.4],'Liver':[113.5,162.4,137.9,111.2,84.2,46.8]}
L={'Liver':1,'Tumour':2,'spleen':3,'Right Kidney':4,'Left Kidney':5}
fig,ax=plt.subplots(1,3,figsize=(174*MM,62*MM),gridspec_kw=dict(width_ratios=[1,1,1.1]))
sl=float(pydicom.dcmread(f'{ROOT}/input/truth_24h_reg.dcm',stop_before_pixels=True).RescaleSlope); Ct=32000*sl
tr=truth('24h'); c=np.round(ndi.center_of_mass(lab==2)).astype(int); xs=(np.arange(-7,8))*VX
prof=lambda v: v[c[0],c[1],c[2]-7:c[2]+8]/Ct
ax[0].step(xs,prof(to_spect(tr,order=0)),where='mid',color='k',lw=1.0,label='Ground truth')
ax[0].plot(xs,prof(to_spect(blur(tr,9.5),order=3)),'--',color='#4A6272',lw=1.1,label='Gaussian 9.5 mm')
ax[0].plot(xs,prof(Sp[2]/k),'-o',color=C['Tumour'],ms=2.5,lw=1.2,label='Measured SPECT')
ax[0].set_xlabel('Position (mm)'); ax[0].set_ylabel('Concentration / true tumour'); ax[0].legend(fontsize=6,loc='upper left'); ax[0].set_ylim(0,2.0)
lo,hi=0.6,1.05
YH=1.16
ax[1].fill_between([lo,hi],[lo*0.95,hi*0.95],[lo*1.05,hi*1.05],color='#EEEEEE',zorder=0); ax[1].plot([lo,hi],[lo,hi],'k-',lw=0.6)
for v,l in L.items():
    for i,t in enumerate(TS):
        true=tab[v][i]/0.5388; meas=Sp[i][lab==l].mean()/k/true
        oth=[x for x in TS if x!=t]; f_=np.median([LF[x][v]['fwhm'] for x in oth]); s_=np.median([LF[x][v]['scale'] for x in oth])
        mod=lambda f: np.interp(f,FG,[MG[t][str(g)][v] for g in FG])/true
        ax[1].plot(mod(9.5),meas,'o',mfc='white',mec=C[v],ms=3.2,mew=0.7)
        ax[1].plot(mod(f_)*s_,meas,'o',color=C[v],mec='k',mew=0.3,ms=3.4,label=NM[v] if i==0 else None)
ax[1].plot([],[],'o',mfc='white',mec='#888888',ms=3.2,label='Global 9.5 mm, no scale'); ax[1].plot([],[],'o',color='#888888',mec='k',mew=0.3,ms=3.4,label='Local model, held-out')
ax[1].set_xlim(lo,hi); ax[1].set_ylim(lo,YH); ax[1].set_xlabel('Predicted VOI / true concentration'); ax[1].set_ylabel('Measured VOI / true concentration')
ax[1].legend(fontsize=5.0,loc='upper left',ncol=2,columnspacing=0.6,handletextpad=0.3)
P=json.load(open('pvc_summary.json'))
x=np.arange(5); w=0.27
un=[P[v]['ad']['unc'] for v in cats]; co=[P[v]['ad']['gtm'] for v in cats]; rc=[P[v]['ad']['rc'] for v in cats]
lo_=[min(P[v]['gtm_ad_rng']) for v in cats]; hi_=[max(P[v]['gtm_ad_rng']) for v in cats]
ax[2].bar(x-w,un,w,color='#9AA9B2',label='Uncorrected')
ax[2].bar(x,co,w,color='#0177A2',label='GTM (PSF 9.5 mm; whiskers 8–12 mm)')
ax[2].errorbar(x,co,yerr=[np.clip(np.array(co)-np.array(lo_),0,None),np.clip(np.array(hi_)-np.array(co),0,None)],fmt='none',ecolor='k',elinewidth=0.7,capsize=2)
ax[2].bar(x+w,rc,w,color='#E97132',label='RC correction (NEMA spheres)')
ax[2].plot(x+w,[P[v]['rc_ct_ad'] for v in cats],'_',color='k',ms=6,mew=1.0)
ax[2].axhline(0,color='k',lw=0.6); ax[2].set_xticks(x); ax[2].set_xticklabels(['Tumour','Spleen','R kidney','L kidney','Liver'],rotation=30,ha='right')
ax[2].set_ylabel('Proportional AD estimate vs GT (%)'); ax[2].set_ylim(-25,52); ax[2].legend(fontsize=5.4,loc='upper right')
for a_,l_ in zip(ax,['(a)','(b)','(c)']): clean(a_); a_.text(-0.22,1.04,l_,transform=a_.transAxes,fontsize=9,fontweight='bold')
fig.tight_layout(); fig.savefig('figs/Fig6.png',dpi=600,facecolor='white'); plt.close(fig)

from resdos.manifest import record
record('figs/Fig4.png','05_figures_tables/f03_fig4_fig6.py',1); record('figs/Fig6.png','05_figures_tables/f03_fig4_fig6.py',3)
print('Fig 4, Fig 6 done')
