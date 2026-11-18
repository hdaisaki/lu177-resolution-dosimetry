"""Components of the 0-mm relative error: voxelisation, VOI-compartment difference, time integration and dose conversion.
Inputs : results/gt_auc_validation.json, gt_closure.json, delineation.json, osirix_vs_slicer.json, holdout_144.json, summary
Outputs: results/components_0mm.json
Manuscript: Table S5, Fig. S2.
"""
from resdos.paths import ROOT
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
V=json.load(open(f'{ROOT}/results/gt_auc_validation.json'))
S=pd.read_csv(f'{ROOT}/results/summary_recommended_models.csv'); R0=pd.read_csv(f'{ROOT}/results/tac_fit_results_0mm.csv')
ho=json.load(open(f'{ROOT}/results/holdout_144.json')); gc=json.load(open(f'{ROOT}/results/gt_closure.json')); de=json.load(open(f'{ROOT}/results/delineation.json')); ovs=json.load(open(f'{ROOT}/results/osirix_vs_slicer.json'))
NM={'Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney','Liver':'Liver'}; order=list(NM)
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
sg=lambda x: (f'{x:+.1f}').replace('-','−')
s1=[]
for v in order:
    r=R0[(R0.VOI==v)&(R0.model==rec[v])].iloc[0]; te=S[S.VOI==v].Teff_72_144_h; e=V[v]['integration']
    oth=[e[m]*100 for m in ['f2','f3','f4b'] if m!=rec[v]]
    s1.append([NM[v],rec[v],f'{te.min():.0f}–{te.max():.0f}' if round(te.min())!=round(te.max()) else f'{te.min():.0f}',f'{S[S.VOI==v].max_rel_residual_pct.max():.1f}',
               f'{r.AD_CV_pct:.1f}',f'{r.extrapolation_pct:.1f}',sg(ho[v][rec[v]]*100),sg(e[rec[v]]*100),sg(e['model-free']*100),f'{sg(min(oth))} to {sg(max(oth))}',sg(V[v]['dose_conversion']*100)])
nom={'Liver':1306.7,'Tumour':15.9,'spleen':124.95,'Right Kidney':111.4,'Left Kidney':133.8}
vol={'Liver':1291.9,'Tumour':15.8,'spleen':121.8,'Right Kidney':110.4,'Left Kidney':133.3}
okey={'Liver':'Liver','Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney'}
s5=[]; comp={}
for v in order:
    o=ovs[okey[v]]; indep=v in ('Liver','Tumour','spleen')
    voi=gc[v]['voi_term']; integ=V[v]['integration'][rec[v]]; conv=V[v]['dose_conversion']
    prod=((1+voi)*(1+integ)*(1+conv)-1)*100
    comp[v]=dict(voi=voi*100,vox=de[v]['vox_gtvoi']*100,delin=de[v]['delineation']*100,integ=integ*100,conv=conv*100,prod=prod,obs=gc[v]['obs'])
    s5.append([NM[v],f"{vol[v]:.1f} ({sg((vol[v]/nom[v]-1)*100)})",(f"{o['dice']:.3f} / {o['msd']:.2f} / {o['hd95']:.1f}" if indep else 'GT compartment'),
               sg(voi*100),sg(de[v]['vox_gtvoi']*100),sg(de[v]['delineation']*100),sg(integ*100),sg(conv*100),sg(prod),sg(gc[v]['obs'])])
json.dump(comp,open(f'{ROOT}/results/components_0mm.json','w'),indent=1)
for r in s1: print(r)
for r in s5: print(r)
