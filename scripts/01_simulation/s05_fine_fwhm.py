"""Phantom FWHM limits for RC >= 0.95/0.90/0.80 on a 1.25-mm FWHM grid (TIA-based RC), and the own-VSV AD at 0 mm.
Inputs : GT maps, derived/seglab.npy, results/same_mask.json
Outputs: results/fine_fwhm.json, results/ownvsv_ad0.json
Manuscript: Table 4 (phantom rows), Fig. 5c.
"""
from resdos.paths import ROOT
import sys, json, numpy as np, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import truth, blur, TS
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); lab=np.load(f'{ROOT}/derived/seglab.npy')
V={'Liver':(1,'f3'),'Tumour':(2,'f4b'),'spleen':(3,'f4b'),'Right Kidney':(4,'f2'),'Left Kidney':(5,'f2')}
tr={t:truth(t) for t in TS}
FWS=[0.0,2.5,3.75,5.0,6.25,7.5,8.75,10.0,11.25,12.5,13.75,15.0]
act={fw:{v:[float(blur(tr[t],fw)[lab==l].mean()) for t in TS] for v,(l,m) in V.items()} for fw in FWS}
I=lambda m,y: ft.integral(m,ft.fit(m,T,np.array(y))[0])
out={}
for v,(l,m) in V.items():
    rc=np.array([I(m,act[f][v]) for f in FWS]); rc=rc/rc[0]
    thr={str(t):(float(np.interp(-t,-rc,FWS)) if rc.min()<=t else None) for t in (0.95,0.90,0.80)}
    out[v]=dict(FW=FWS,RC=rc.tolist(),thr=thr)
    print(f"{v:12s} thresholds RC 0.95/0.90/0.80: {[None if x is None else round(x,1) for x in thr.values()]}")
json.dump(out,open(f'{ROOT}/results/fine_fwhm.json','w'),indent=1)
# own-VSV 0-mm AD (same mask) for reference
R=json.load(open(f'{ROOT}/results/same_mask.json'))
GT={'Tumour':28.1,'spleen':5.39,'Right Kidney':3.6,'Left Kidney':3.61,'Liver':1.70}
ref={}
for v,(l,m) in V.items():
    a0=I(m,R['0.0']['dr'][v])/1000; ref[v]=a0; print(f'{v:12s} own-VSV AD at 0 mm {a0:.2f} Gy ({(a0/GT[v]-1)*100:+.1f}% vs GT)')
json.dump(ref,open(f'{ROOT}/results/ownvsv_ad0.json','w'))
