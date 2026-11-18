"""Medians of local FWHM and scale over the six time points and predicted vs measured recovery.
Outputs: results/spect/local_decomp.json
"""
from resdos.paths import ROOT
import json, numpy as np
TS=['1h','4h','24h','40h','72h','144h']
L={t:json.load(open(f'{ROOT}/results/spect/local_fit_{t}.json')) for t in TS}
R=json.load(open(f'{ROOT}/results/spect/organ_results_k_abs.json'))['organs']
out={}
for v in ['Tumour','spleen','Right Kidney','Left Kidney','Liver']:
    f=np.array([L[t][v]['fwhm'] for t in TS]); s=np.array([L[t][v]['scale'] for t in TS])
    qm=R[v]['Q_model']; xs=np.array(sorted(float(x) for x in qm)); ys=np.array([qm[str(x)] if str(x) in qm else qm[f'{x}'] for x in xs])
    fm=float(np.median(f)); sm=float(np.median(s)); qloc=float(np.interp(fm,xs,ys)); pred=qloc*sm; meas=R[v]['Q_spect']
    out[v]=dict(fwhm_local_median=fm,fwhm_range=[float(f.min()),float(f.max())],scale_median=sm,scale_range=[float(s.min()),float(s.max())],Q_model_local=qloc,Q_pred=pred,Q_meas=meas)
    print(f"{v:12s} local FWHM {fm:4.1f} ({f.min():.1f}-{f.max():.1f}) | scale {sm:.3f} ({s.min():.3f}-{s.max():.3f}) | pred {pred:.3f} vs meas {meas:.3f} ({(meas/pred-1)*100:+.1f}%)")
json.dump(out,open(f'{ROOT}/results/spect/local_decomp.json','w'),indent=1)
