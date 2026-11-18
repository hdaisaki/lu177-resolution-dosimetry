"""Sensitivity of RC and RE to the choice of the integration function (f2, f3, f4b).
Inputs : results/tac_fit_results_{fw}mm.csv
Outputs: results/model_sensitivity.json
Manuscript: Table S4.
"""
from resdos.paths import ROOT
import pandas as pd, numpy as np, json
FWS=[0,5,7.5,10,12.5,15]
GT={'Tumour':28.1,'spleen':5.39,'Right Kidney':3.6,'Left Kidney':3.61,'Liver':1.70}
D=pd.concat([pd.read_csv(f'{ROOT}/results/tac_fit_results_{f:g}mm.csv').assign(FWHM=f) for f in FWS])
out={}
for v in GT:
    ref=D[(D.VOI==v)].drop_duplicates('FWHM').sort_values('FWHM').reference_AD_Gy.values; rcr=ref/ref[0]
    print(f'== {v}  model-free ref RC15 {rcr[-1]:.3f}, RE0 {(ref[0]/GT[v]-1)*100:+.1f}%, RE15 {(ref[-1]/GT[v]-1)*100:+.1f}%')
    for m in ['f2','f3','f4b']:
        d=D[(D.VOI==v)&(D.model==m)].sort_values('FWHM'); ad=d.AD_Gy.values; rc=ad/ad[0]; re=(ad/GT[v]-1)*100
        out.setdefault(v,{})[m]=dict(AD=ad.tolist(),RC=rc.tolist(),RE=re.tolist(),maxres=float(d.max_rel_residual_pct.max()),Teff=d.Teff_72_144_h.tolist())
        print(f'  {m:3s} RC15 {rc[-1]:.3f} (vs ref {rc[-1]-rcr[-1]:+.3f}) | RE0 {re[0]:+6.1f}% RE15 {re[-1]:+6.1f}% | maxres {d.max_rel_residual_pct.max():5.1f}% | Teff {d.Teff_72_144_h.iloc[0]:.0f}->{d.Teff_72_144_h.iloc[-1]:.0f} h')
    out[v]['ref']=dict(RC=rcr.tolist(),AD=ref.tolist(),RE=((ref/GT[v]-1)*100).tolist())
    rcs=np.array([out[v][m]['RC'] for m in ['f2','f3','f4b']]); res=np.array([out[v][m]['RE'] for m in ['f2','f3','f4b']])
    print(f'  spread across functions: RC {np.max(rcs.max(0)-rcs.min(0)):.3f} | RE {np.max(res.max(0)-res.min(0)):.1f} points')
json.dump(out,open(f'{ROOT}/results/model_sensitivity.json','w'),indent=1)
