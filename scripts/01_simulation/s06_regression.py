"""Descriptive regression of RE(FWHM) (linear and quadratic) and FWHM thresholds of the six-point series.
Inputs : results/summary_recommended_models.csv
Outputs: results/stats.json, results/regression_desc.json
Manuscript: Table S3, Fig. 4.
"""
from resdos.paths import ROOT
import numpy as np, pandas as pd, json
from scipy import stats
S=pd.read_csv(f'{ROOT}/results/summary_recommended_models.csv')
FW=np.array([0,5,7.5,10,12.5,15.])
out={}
for v in ['Tumour','spleen','Right Kidney','Left Kidney','Liver']:
    d=S[S.VOI==v].sort_values('FWHM_mm'); re=d.RE_vs_GT_pct.values; rc=d.RC_vs_0mm.values
    r={}
    for name,y in [('RE',re),('RC',rc)]:
        pl=np.polyfit(FW,y,1); pq=np.polyfit(FW,y,2)
        rssl=((y-np.polyval(pl,FW))**2).sum(); rssq=((y-np.polyval(pq,FW))**2).sum()
        F=(rssl-rssq)/(rssq/(6-3)); p=1-stats.f.cdf(F,1,3)
        aic=lambda rss,k: 6*np.log(rss/6)+2*k+2*k*(k+1)/(6-k-1)
        r[name]=dict(quad=pq.tolist(),lin=pl.tolist(),F=F,p=p,dAICc=aic(rssq,3)-aic(rssl,2),R2q=1-rssq/((y-y.mean())**2).sum(),R2l=1-rssl/((y-y.mean())**2).sum())
    # FWHM thresholds for RC 0.95 / 0.90 / 0.80 by linear interpolation of measured RC
    thr={}
    for t in (0.95,0.90,0.80):
        thr[t]=float(np.interp(-t,-rc,FW)) if rc.min()<=t else None
    r['thr']=thr; out[v]=r
    print(f"{v:12s} RE: quad vs lin F={r['RE']['F']:.1f} p={r['RE']['p']:.4f} dAICc={r['RE']['dAICc']:+.1f} | RC: p={r['RC']['p']:.4f} | FWHM for RC 0.95/0.90/0.80: {[None if x is None else round(x,1) for x in thr.values()]}")
json.dump(out,open(f'{ROOT}/results/stats.json','w'),indent=1)
import json, numpy as np, pandas as pd
S=pd.read_csv(f'{ROOT}/results/summary_recommended_models.csv'); FW=np.array([0,5,7.5,10,12.5,15.])
out={}
for v in ['Tumour','spleen','Right Kidney','Left Kidney','Liver']:
    y=S[S.VOI==v].sort_values('FWHM_mm').RE_vs_GT_pct.values
    pl=np.polyfit(FW,y,1); pq=np.polyfit(FW,y,2)
    rl=y-np.polyval(pl,FW); rq=y-np.polyval(pq,FW)
    r2=lambda r: 1-(r**2).sum()/((y-y.mean())**2).sum()
    out[v]=dict(lin=pl.tolist(),quad=pq.tolist(),R2l=r2(rl),R2q=r2(rq),maxres_l=float(abs(rl).max()),maxres_q=float(abs(rq).max()))
    print(f"{v:12s} linear R2 {r2(rl):.3f} max|res| {abs(rl).max():.2f} pt | quadratic R2 {r2(rq):.3f} max|res| {abs(rq).max():.2f} pt | a {pq[0]:.4f} b {pq[1]:.3f} c {pq[2]:.1f}")
json.dump(out,open(f'{ROOT}/results/regression_desc.json','w'),indent=1)
