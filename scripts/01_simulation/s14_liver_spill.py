"""Decomposition of the liver-VOI time-integrated concentration by source compartment (liver, tumour, background,
kidneys, spleen) at 0 and 15 mm, using one-compartment maps and a common linear (trapezoidal) integral.
Inputs : GT maps, derived/frac_tp_avg.npy, derived/seglab.npy
Outputs: results/liver_spill_full.json
Manuscript: Results 3.2, Table S10.
"""
from resdos.paths import ROOT
import sys, json, numpy as np
from resdos.common import truth, blur, TS
T=np.array([1,4,24,40,72,144.])
F=np.load(f'{ROOT}/derived/frac_tp_avg.npy'); CN=['Air','BG','Lung','Liver','Spleen','KidLC','KidLM','KidRC','KidRM','Tumor']
lab=np.load(f'{ROOT}/derived/seglab.npy'); liv=lab==1
tr=[truth(t) for t in TS]
conc={}
for k,c in enumerate(CN):
    pure=F[...,k]>0.95
    conc[c]=[float(tr[i][pure].mean()) if pure.sum()>5 else 0.0 for i in range(6)]
groups={'Liver':['Liver'],'Tumour':['Tumor'],'Background (incl. lung)':['BG','Lung'],'Spleen':['Spleen'],'Kidneys':['KidLC','KidLM','KidRC','KidRM']}
def comp_map(g,i): return sum(conc[c][i]*F[...,CN.index(c)] for c in groups[g])
y0=np.array([tr[i][liv].mean() for i in range(6)]); lam=max(-np.polyfit(T[-3:],np.log(y0[-3:]),1)[0],4.345e-3)
def integ(y): y=np.asarray(y); return 0.5*y[0]*T[0]+np.trapezoid(y,T)+y[-1]/lam
res={}
for fw in [0.0,15.0]:
    tot=[blur(tr[i],fw)[liv].mean() for i in range(6)]
    res[fw]={'Total':integ(tot)}
    for g in groups:
        res[fw][g]=integ([blur(comp_map(g,i),fw)[liv].mean() for i in range(6)])
ref=res[0.0]['Total']
for g in list(groups)+['Total']: print(f"{g:26s} {res[0.0][g]/ref*100:8.2f} {res[15.0][g]/ref*100:8.2f}")
json.dump({str(k):v for k,v in res.items()},open(f'{ROOT}/results/liver_spill_full.json','w'),indent=1)
