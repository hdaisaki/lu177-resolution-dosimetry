"""Supplementary Fig. S3: NEMA recovery coefficients, RC curves, Gaussian model and anthropomorphic recoveries.
Outputs: results/figs/FigS3.png, results/nema/gauss_nema.json
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results')
# Supplementary Fig. S3: NEMA IEC recovery coefficients and anthropomorphic-phantom recoveries
import numpy as np, json, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
import resdos.sphere_model as sm
from resdos.manifest import record
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica','DejaVu Sans'],'font.size':8,'axes.linewidth':0.6,'legend.frameon':False})
R=json.load(open(f'{ROOT}/results/nema/nema_rc.json')); Q=json.load(open(f'{ROOT}/results/spect/organ_results_final.json'))['organs']
rc_fun=lambda V,al,be: 1.0/(1.0+(al/V)**be)
geo=sorted(R['results']['CT per sphere | ICF 42.82 cps/MBq'],key=lambda r:r['d'])
od3d=[(0.52,0.28),(1.21,0.44),(2.5,0.57),(5.61,0.60),(11.57,0.63),(26.68,0.69)]
rng=np.random.default_rng(0); gauss=[]
for r in geo:
    vals=[]
    for k in range(16):
        n=sm.box_size(r['V'],4.419644,20); occ,_=sm.sphere_occupancy(r['V'],4.419644,n,rng.uniform(-0.5,0.5,3)); voi=occ>=0.5
        vals.append(sm.blur(occ,9.5,4.419644)[voi].mean())
    gauss.append(float(np.mean(vals)))
json.dump(dict(gauss=gauss),open(f'{ROOT}/results/nema/gauss_nema.json','w'))
Vg=np.logspace(np.log10(0.3),np.log10(2000),300)
fig,ax=plt.subplots(figsize=(120/25.4,85/25.4))
ax.plot(Vg,rc_fun(Vg,2.53,0.41),'-',color='#0177A2',lw=1.3,label='RC curve, OpenDose3D (α 2.53 mL, β 0.41)')
ax.plot(Vg,rc_fun(Vg,2.16,0.39),'--',color='#888888',lw=1.0,label='RC curve, CT-positioned spheres (α 2.16 mL, β 0.39)')
ax.plot([v for v,_ in od3d],[r for _,r in od3d],'o',color='#0177A2',ms=4.2,label='NEMA spheres, volume threshold (OpenDose3D)')
ax.plot([r['V'] for r in geo],[r['rc_conc'] for r in geo],'s',mfc='white',mec='#555555',ms=3.8,label='NEMA spheres, CT-positioned VOI')
ax.plot([r['V'] for r in geo],gauss,'x',color='#7C3AED',ms=4.5,mew=1.0,label='Gaussian 9.5 mm model, same spheres')
col={'Tumour':'#D9731E','spleen':'#1F8A94','Right Kidney':'#2F4F9E','Left Kidney':'#7C99D4','Liver':'#8A7650'}
NM={'Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney','Liver':'Liver'}
VOL={'Tumour':15.8,'spleen':121.8,'Right Kidney':110.4,'Left Kidney':133.3,'Liver':1291.9}
for v in VOL: ax.plot(VOL[v],Q[v]['Q_spect'],'D',color=col[v],ms=5,mec='k',mew=0.4,label=f'{NM[v]}, anthropomorphic phantom')
ax.set_xscale('log'); ax.set_xlabel('Volume (mL)'); ax.set_ylabel('Recovery (VOI / true concentration)'); ax.set_ylim(0,1.05)
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.legend(fontsize=5.5,loc='lower right')
fig.tight_layout(); fig.savefig('figs/FigS3.png',dpi=600,facecolor='white'); plt.close(fig)
record('figs/FigS3.png','05_figures_tables/f05_figS3.py',1); print('Fig S3 done')
