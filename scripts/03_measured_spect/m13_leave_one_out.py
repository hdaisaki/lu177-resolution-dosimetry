"""Leave-one-time-point-out prediction of the VOI concentrations with the local model.
Outputs: results/spect/loo_local.json
Manuscript: Table 5, Table S6, Fig. 6b.
"""
from resdos.paths import ROOT
import json, numpy as np
TS=['1h','4h','24h','40h','72h','144h']
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy')
k=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['k_abs']
MG={t:json.load(open(f'{ROOT}/results/spect/model_grid_{t}.json')) for t in TS}
LF={t:json.load(open(f'{ROOT}/results/spect/local_fit_{t}.json')) for t in TS}
V={'Liver':1,'Tumour':2,'spleen':3,'Right Kidney':4,'Left Kidney':5}
FG=np.array(sorted(float(f) for f in MG['1h']))
res={}
print('Leave-one-time-point-out prediction of VOI mean (SPECT/model) — errors in %')
for v,l in V.items():
    e_loo=[]; e_glob=[]; e_fit=[]
    for i,t in enumerate(TS):
        meas=S[i][lab==l].mean()/k
        others=[x for x in TS if x!=t]
        f_loo=np.median([LF[x][v]['fwhm'] for x in others]); s_loo=np.median([LF[x][v]['scale'] for x in others])
        mod=lambda f: np.interp(f,FG,[MG[t][str(g)][v] for g in FG])
        pred_loo=mod(f_loo)*s_loo; pred_glob=mod(9.5)*1.0; pred_fit=mod(LF[t][v]['fwhm'])*LF[t][v]['scale']
        e_loo.append(meas/pred_loo-1); e_glob.append(meas/pred_glob-1); e_fit.append(meas/pred_fit-1)
    e_loo=np.array(e_loo)*100; e_glob=np.array(e_glob)*100; e_fit=np.array(e_fit)*100
    res[v]=dict(loo=e_loo.tolist(),glob=e_glob.tolist(),insample=e_fit.tolist())
    print(f"{v:12s} LOO local: {np.round(e_loo,1)} (MAE {np.mean(abs(e_loo)):.1f}) | global 9.5 mm no scale: MAE {np.mean(abs(e_glob)):.1f} | in-sample: MAE {np.mean(abs(e_fit)):.1f}")
json.dump(res,open(f'{ROOT}/results/spect/loo_local.json','w'),indent=1)
