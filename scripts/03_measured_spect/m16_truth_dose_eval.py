"""Voxel-based dose evaluation of the exact (unblurred) GT image: binary VOI, occupancy-weighted and interior voxels.
Outputs: results/truth_dose_eval.json
Manuscript: Supplementary Methods S4 (background for Table S12).
"""
from resdos.paths import ROOT
import sys, json, numpy as np
from resdos.common import truth, TS
from scipy.signal import fftconvolve
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); K=np.load(f'{ROOT}/derived/vsv_lu177_4.42mm.npy'); VV=0.4419644**3
F=np.load(f'{ROOT}/derived/frac_tp_avg.npy'); lab=np.load(f'{ROOT}/derived/seglab.npy')
COMP={'Liver':F[...,3],'Tumour':F[...,9],'spleen':F[...,4],'Right Kidney':F[...,7]+F[...,8],'Left Kidney':F[...,5]+F[...,6]}
L={'Liver':1,'Tumour':2,'spleen':3,'Right Kidney':4,'Left Kidney':5}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
GT={'Tumour':28.12285,'spleen':5.382068,'Right Kidney':3.598351,'Left Kidney':3.599789,'Liver':1.702361}
drb={v:[] for v in L}; drc={v:[] for v in L}; drpure={v:[] for v in L}
for t in TS:
    a=truth(t); D=fftconvolve(a*0.5388*VV/1000,K,mode='same')*3600
    for v in L:
        drb[v].append(D[lab==L[v]].mean()); w=COMP[v]; drc[v].append((D*w).sum()/w.sum())
        core=w>0.99; drpure[v].append(D[core].mean() if core.sum()>0 else np.nan)
I=lambda m,y: ft.integral(m,ft.fit(m,T,np.asarray(y,float))[0])
out={}
for v in L:
    b=I(rec[v],drb[v])/1000; c=I(rec[v],drc[v])/1000; p=I(rec[v],drpure[v])/1000
    out[v]=dict(binary=(b/GT[v]-1)*100,comp=(c/GT[v]-1)*100,interior=(p/GT[v]-1)*100)
    print(f"{v:12s} {out[v]['binary']:+6.1f} | {out[v]['comp']:+6.1f} | {out[v]['interior']:+6.1f}")
json.dump(out,open(f'{ROOT}/results/truth_dose_eval.json','w'),indent=1)
