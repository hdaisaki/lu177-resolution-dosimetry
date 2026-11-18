"""Hold-out test of the integration functions on the Table 1 curves (fit 1-72 h, predict 144 h).
Outputs: results/holdout_144.json
Manuscript: Table S1.
"""
from resdos.paths import ROOT
import numpy as np, json, importlib.util
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.])
tab={'Tumour':[1916.5,1809.2,1864.5,1713.2,1396.6,825.9],'spleen':[374.5,428.0,500.4,433.7,314.4,151.7],'Right Kidney':[494.9,466.2,403.5,326.1,214.4,82.7],
     'Left Kidney':[490.0,456.5,404.2,325.1,212.4,82.4],'Liver':[113.5,162.4,137.9,111.2,84.2,46.8]}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
lp=4.345e-3
out={}
print('Hold-out test on the GT (Table 1) curves: fit 1-72 h, predict 144 h (error %)')
for v,y in tab.items():
    y=np.array(y); r={}
    for m in ['f2','f3','f4b']:
        p=ft.fit(m,T[:5],y[:5])[0]; r[m]=ft.MODELS[m][0](144.0,*p)/y[5]-1
    lam=max(-np.polyfit(T[2:5],np.log(y[2:5]),1)[0],lp); r['mono-exp tail']=y[4]*np.exp(-lam*72)/y[5]-1
    out[v]=r
    print(f"{v:12s} "+'  '.join(f"{m}{'*' if m==rec[v] else ''} {e*100:+6.1f}" for m,e in r.items()))
json.dump(out,open(f'{ROOT}/results/holdout_144.json','w'),indent=1)
