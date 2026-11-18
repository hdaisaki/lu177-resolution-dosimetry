"""Sensitivity of uncorrected and corrected errors to the absolute calibration over assumed ranges
(Jaszczak-based factor x0.965/1/1.035 with refitted NEMA RC curve; factor four x0.98/1/1.04).
Outputs: results/calib_sensitivity.json
Manuscript: Methods 2.6, Results 3.5, Supplementary Methods S7, Table S13.
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results')
# Sensitivity of uncorrected and corrected TIAC / proportional AD errors to the absolute calibration.
# c: error of the Jaszczak-based factor (common to NEMA and anthropomorphic images; NEMA RCs scale by 1/c and the RC curve is refitted)
# s: error of the factor four (anthropomorphic images only)
import json, numpy as np, itertools
from scipy.optimize import curve_fit
P=json.load(open('pvc_summary.json')); V=json.load(open('gt_auc_validation.json'))
order=['Tumour','spleen','Right Kidney','Left Kidney','Liver']
VOL={'Tumour':15.8,'spleen':121.8,'Right Kidney':110.4,'Left Kidney':133.3,'Liver':1291.9}
od3d_V=np.array([0.52,1.21,2.5,5.61,11.57,26.68]); od3d_rc=np.array([0.28,0.44,0.57,0.60,0.63,0.69])
rc_fun=lambda v,a,b: 1/(1+(a/v)**b)
fit=lambda y: curve_fit(rc_fun,od3d_V,y,p0=[2.5,0.4],bounds=([1e-3,0.05],[100,5]))[0]
a0,b0=fit(od3d_rc)
fit_implied=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))
implied=np.array([r[2]/r[3] for r in fit_implied['abs_rows']][:5])     # measured/true total activity at 1-72 h with 10.705
C=[0.965,1.0,1.035]; Sx=[0.98,1.0,1.04]
out={'basis':{'c':'Jaszczak activity 441.2 +/- 13.8 MBq (3.1%, README) and NEMA total-activity difference +3.5%',
              's':f'factor implied by the known total activity at 1-72 h: {10.705*implied.min():.2f}-{10.705*implied.max():.2f} cps/MBq ({(implied.min()-1)*100:+.1f}% to {(implied.max()-1)*100:+.1f}%)'},
     'refit':{},'scen':{}}
for c in C:
    a,b=fit(od3d_rc/c); out['refit'][str(c)]=[float(a),float(b)]
for v in order:
    t=P[v]['tiac']; conv=V[v]['dose_conversion']; rows=[]
    for c,s in itertools.product(C,Sx):
        fA=1/(c*s); a,b=out['refit'][str(c)]
        unc=(1+t['unc']/100)*fA-1
        gtm=(1+t['gtm']/100)*fA-1
        rc=(1+t['unc']/100)*fA/rc_fun(VOL[v],a,b)-1
        ad=lambda e: ((1+e)*(1+conv)-1)*100
        rows.append(dict(c=c,s=s,tiac=dict(unc=unc*100,gtm=gtm*100,rc=rc*100),ad=dict(unc=ad(unc),gtm=ad(gtm),rc=ad(rc)),
                         closer_tiac='GTM' if abs(gtm)<abs(rc) else 'RC',closer_ad='GTM' if abs(ad(gtm))<abs(ad(rc)) else 'RC'))
    out['scen'][v]=rows
json.dump(out,open('calib_sensitivity.json','w'),indent=1)
print('RC refit (alpha, beta):',{k:np.round(v,2).tolist() for k,v in out['refit'].items()}, ' base', round(a0,2), round(b0,2))
print(out['basis'])
for v in order:
    R=out['scen'][v]
    rng=lambda k,m: (min(r[k][m] for r in R),max(r[k][m] for r in R))
    base=[r for r in R if r['c']==1 and r['s']==1][0]
    print(f"{v:12s} TIAC GTM {base['tiac']['gtm']:+5.1f} [{rng('tiac','gtm')[0]:+5.1f},{rng('tiac','gtm')[1]:+5.1f}]  RC {base['tiac']['rc']:+5.1f} [{rng('tiac','rc')[0]:+5.1f},{rng('tiac','rc')[1]:+5.1f}]  "
          f"| AD GTM [{rng('ad','gtm')[0]:+5.1f},{rng('ad','gtm')[1]:+5.1f}] RC [{rng('ad','rc')[0]:+5.1f},{rng('ad','rc')[1]:+5.1f}] unc [{rng('ad','unc')[0]:+5.1f},{rng('ad','unc')[1]:+5.1f}]"
          f" | closer (AD): {sorted(set(r['closer_ad'] for r in R))} ({sum(r['closer_ad']=='GTM' for r in R)}/9 GTM)")
