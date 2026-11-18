"""Per-time-point concentrations for uncorrected SPECT and GTM (separate cortex and medulla; argument: PSF FWHM),
voxel-based dose rates of uncorrected and GTM-corrected images, self-dose fraction.
Usage: python m14_unified_eval.py 8.0 | 9.5 | 12.0
Outputs: results/spect/unified_{FW}.json
"""
from resdos.paths import ROOT
import sys, json, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.signal import fftconvolve
from scipy.ndimage import binary_erosion
import resdos.fit_tac as ft
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy'); frac=np.load(f'{ROOT}/results/spect/segfrac_spect.npy')
K=np.load(f'{ROOT}/derived/vsv_lu177_4.42mm.npy'); kabs=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['k_abs']
KT=0.5388; VV=(VX/10)**3; FW=float(sys.argv[1]) if len(sys.argv)>1 else 9.5
F=np.load(f'{ROOT}/derived/frac_tp_avg.npy')
cR=to_spect(F[...,7],order=1); mR=to_spect(F[...,8],order=1); cL=to_spect(F[...,5],order=1); mL=to_spect(F[...,6],order=1)
body=to_spect(truth('24h'),order=0)>0
tpl=[frac[1],frac[2],frac[3],cR,mR,cL,mL]
rows=[lab==1,lab==2,lab==3,(cR>=0.5)&(lab==4),(mR>=0.5)&(lab==4),(cL>=0.5)&(lab==5),(mL>=0.5)&(lab==5)]
bg=np.clip(body.astype(float)-sum(tpl),0,1); tpl=tpl+[bg]; rows=rows+[binary_erosion(body,iterations=2)&(lab==0)]
B=[blur(t,FW) for t in tpl]; W=np.array([[b[r].mean() for b in B] for r in rows])
V={'Liver':1,'Tumour':2,'spleen':3,'Right Kidney':4,'Left Kidney':5}
out={v:dict(conc_unc=[],conc_gtm=[],dr_unc=[],dr_gtm=[],dr_self=[],dr_gtm_comp=[]) for v in V}
COMP={'Liver':tpl[0],'Tumour':tpl[1],'spleen':tpl[2],'Right Kidney':tpl[3]+tpl[4],'Left Kidney':tpl[5]+tpl[6]}
for i in range(6):
    m=np.array([S[i][r].mean() for r in rows])/kabs
    c=np.linalg.solve(W,m)*KT                                    # kBq/mL per compartment
    img=sum(ci/KT*t for ci,t in zip(c,tpl))                      # truth units
    Dc=fftconvolve(img*KT*VV/1000.0,K,mode='same')*3600
    Au=S[i]/(10.705*1800.0); Du=fftconvolve(Au,K,mode='same')*3600
    gtm_c={'Liver':c[0],'Tumour':c[1],'spleen':c[2],'Right Kidney':0.679*c[3]+0.321*c[4],'Left Kidney':0.679*c[5]+0.321*c[6]}
    for v,l in V.items():
        msk=lab==l
        out[v]['conc_unc'].append(float(S[i][msk].mean()/kabs*KT)); out[v]['conc_gtm'].append(float(gtm_c[v]))
        out[v]['dr_unc'].append(float(Du[msk].mean())); out[v]['dr_gtm'].append(float(Dc[msk].mean()))
        w=COMP[v]; out[v]['dr_gtm_comp'].append(float((Dc*w).sum()/w.sum()))
        Ds=fftconvolve(Au*msk,K,mode='same')*3600; out[v]['dr_self'].append(float(Ds[msk].mean()))
json.dump(out,open(f'{ROOT}/results/spect/unified_{FW}.json','w'))
print('done',FW)
