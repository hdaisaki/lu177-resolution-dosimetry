"""Supplementary Fig. S2: components of the 0-mm relative error.
Outputs: results/figs/FigS2.png
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results')
# Supplementary Fig. S2: components of the 0-mm relative error (SAAM II-based decomposition)
import json, numpy as np, matplotlib, sys
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from resdos.manifest import record
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica','DejaVu Sans'],'font.size':8,'axes.linewidth':0.6,'legend.frameon':False})
comp=json.load(open('components_0mm.json'))
order=['Tumour','spleen','Right Kidney','Left Kidney','Liver']; NM={'Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney','Liver':'Liver'}
MM=1/25.4; fig,ax=plt.subplots(figsize=(150*MM,88*MM)); x=np.arange(5); w=0.16
ser=[('Voxelisation (GT-compartment VOI)','vox','#0177A2'),('VOI–compartment difference','delin','#7FB8D4'),('Time integration vs SAAM II AUC','integ','#E97132'),
     ('Dose conversion (OpenDose3D vs GT)','conv','#F2B880'),('Observed RE at 0 mm','obs','#0B2E44')]
for i,(lab,k,col) in enumerate(ser):
    vals=[comp[v][k] for v in order]; ax.bar(x+(i-2)*w,vals,w,color=col,label=lab)
    for xi,vi in zip(x+(i-2)*w,vals): ax.text(xi,vi-0.3 if vi<0 else vi+0.3,f'{vi:.1f}',ha='center',va='top' if vi<0 else 'bottom',fontsize=5,rotation=90)
ax.axhline(0,color='k',lw=0.6); ax.set_xticks(x); ax.set_xticklabels([NM[v] for v in order]); ax.set_ylabel('Contribution at 0 mm (%)'); ax.set_ylim(-19,7)
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.legend(fontsize=5.8,loc='upper center',bbox_to_anchor=(0.5,-0.10),ncol=3)
fig.tight_layout(); fig.savefig('figs/FigS2.png',dpi=600,facecolor='white'); plt.close(fig)
record('figs/FigS2.png','05_figures_tables/f04_figS2.py',1); print('Fig S2 done')
