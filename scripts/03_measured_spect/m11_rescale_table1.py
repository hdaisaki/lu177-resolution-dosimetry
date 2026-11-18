"""Express recoveries relative to the Table 1 concentrations (the GT maps use time-invariant occupancy).
Outputs: results/spect/organ_results_final.json, gtm_final.json, local_decomp_final.json
Manuscript: Table 5 (Q, local FWHM, scale), Fig. S3.
"""
from resdos.paths import ROOT
import numpy as np, json, pydicom, importlib.util
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.])
KT=1916.5/(32000*float(pydicom.dcmread(f'{ROOT}/input/truth_1h_reg.dcm',stop_before_pixels=True).RescaleSlope))
tab={'Tumour':[1916.5,1809.2,1864.5,1713.2,1396.6,825.9],'spleen':[374.5,428.0,500.4,433.7,314.4,151.7],'Right Kidney':[494.9,466.2,403.5,326.1,214.4,82.7],
     'Left Kidney':[490.0,456.5,404.2,325.1,212.4,82.4],'Liver':[113.5,162.4,137.9,111.2,84.2,46.8]}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
R=json.load(open(f'{ROOT}/results/spect/organ_results_k_abs.json')); G=json.load(open(f'{ROOT}/results/spect/gtm_k_abs.json')); L=json.load(open(f'{ROOT}/results/spect/local_decomp.json'))
I=lambda m,y: ft.integral(m,ft.fit(m,T,np.asarray(y,float))[0])
for v in tab:
    f=I(rec[v],R['organs'][v]['tac_true'])/I(rec[v],np.array(tab[v])/KT)
    o=R['organs'][v]
    o['Q_spect']*=f; o['Q_model']={k:x*f for k,x in o['Q_model'].items()}
    for kk in G[v]: G[v][kk]*=f
    L[v]['Q_pred']*=f; L[v]['Q_meas']*=f; L[v]['Q_model_local']*=f
    o['table1_factor']=f
    print(f"{v:12s} factor {f:.4f} | Q_SPECT {o['Q_spect']:.3f} | Q_model 0/9.5 {o['Q_model']['0.0']:.3f}/{o['Q_model']['9.5']:.3f} | local pred {L[v]['Q_pred']:.3f} | GTM 8/9.5/11/12 {G[v]['8.0']:.3f}/{G[v]['9.5']:.3f}/{G[v]['11.0']:.3f}/{G[v]['12.0']:.3f}")
json.dump(R,open(f'{ROOT}/results/spect/organ_results_final.json','w'),indent=1); json.dump(G,open(f'{ROOT}/results/spect/gtm_final.json','w'),indent=1); json.dump(L,open(f'{ROOT}/results/spect/local_decomp_final.json','w'),indent=1)
