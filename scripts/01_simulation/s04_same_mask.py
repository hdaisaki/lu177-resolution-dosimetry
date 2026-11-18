"""TIA- and AD-based RC with identical masks: dose rates recalculated from the (blurred) GT maps by voxel S-value
convolution, VOI means of activity and dose rate, time integration with the fixed functions.
Inputs : GT maps, derived/seglab.npy, derived/vsv_lu177_4.42mm.npy
Outputs: results/same_mask.json, results/same_mask_rc.json
Manuscript: Table 3 (last column), Results 3.2. Usage: python s04_same_mask.py
"""
from resdos.paths import ROOT
import sys, json, numpy as np; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import truth, TS
from scipy.signal import fftconvolve
K=np.load(f'{ROOT}/derived/vsv_lu177_4.42mm.npy'); lab=np.load(f'{ROOT}/derived/seglab.npy'); VV=0.4419644**3
V={'Liver':1,'Tumour':2,'spleen':3,'Right Kidney':4,'Left Kidney':5}
res={}
for fw in [0.0,5.0,7.5,10.0,12.5,15.0]:
    act={v:[] for v in V}; dr={v:[] for v in V}
    for t in TS:
        a=truth(t,fw); D=fftconvolve(a*0.5388*VV/1000.0,K,mode='same')*3600
        for v,l in V.items(): act[v].append(float(a[lab==l].mean())); dr[v].append(float(D[lab==l].mean()))
    res[str(fw)]=dict(act=act,dr=dr); print(fw,'done')
json.dump(res,open(f'{ROOT}/results/same_mask.json','w'))

# ---- RC from identical masks
import json, numpy as np, importlib.util
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); R=json.load(open(f'{ROOT}/results/same_mask.json'))
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
I=lambda m,y: ft.integral(m,ft.fit(m,T,np.array(y))[0])
FW=['0.0','5.0','7.5','10.0','12.5','15.0']; out={}
for v,m in rec.items():
    ia=[I(m,R[f]['act'][v]) for f in FW]; idr=[I(m,R[f]['dr'][v]) for f in FW]
    rca=np.array(ia)/ia[0]; rcd=np.array(idr)/idr[0]
    out[v]=dict(RC_TIA=rca.tolist(),RC_AD=rcd.tolist())
    print(f'{v:12s} RC_TIA {np.round(rca[1:],3)} | RC_AD {np.round(rcd[1:],3)} | max|diff| {np.max(abs(rca-rcd)):.3f}')
json.dump(out,open(f'{ROOT}/results/same_mask_rc.json','w'),indent=1)
