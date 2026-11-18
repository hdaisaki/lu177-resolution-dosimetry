"""Collect the selected (fixed) function per VOI over FWHM into one table with RC (relative to 0 mm) and RE (relative to GT).
Inputs : results/tac_fit_results_{fw}mm.csv
Outputs: results/summary_recommended_models.csv
Manuscript: Tables 2, 3, S2; Figs 3, 4.
"""
from resdos.paths import ROOT
import pandas as pd, os
GT={'Tumour':28.1,'spleen':5.39,'Right Kidney':3.6,'Left Kidney':3.61,'Liver':1.70}
NK={'Tumour':[28.41,26.16,24.12,21.81,19.65,17.56],'spleen':[4.95,4.71,4.47,4.22,3.96,3.71],'Right Kidney':[3.36,3.16,2.98,2.80,2.63,2.47],
    'Left Kidney':[3.37,3.18,3.01,2.84,2.67,2.51],'Liver':[1.53,1.51,1.50,1.48,1.47,1.46]}
FWS=[0,5,7.5,10,12.5,15]
rows=[]
for fw in FWS:
    f=f'{ROOT}/results/tac_fit_results_{fw:g}mm.csv'
    if not os.path.exists(f): continue
    d=pd.read_csv(f); d=d[d.recommended]
    for _,r in d.iterrows():
        rows.append(dict(FWHM_mm=fw,VOI=r.VOI,model=r.model,AD_Gy=r.AD_Gy,AD_SD_Gy=r.AD_SD_Gy,Teff_72_144_h=r.Teff_72_144_h,
                         max_rel_residual_pct=r.max_rel_residual_pct,parameter_at_bound=r.parameter_at_bound))
S=pd.DataFrame(rows)
S['RC_vs_0mm']=S.apply(lambda r:r.AD_Gy/S[(S.FWHM_mm==0)&(S.VOI==r.VOI)].AD_Gy.iloc[0],axis=1)
S['RE_vs_GT_pct']=S.apply(lambda r:(r.AD_Gy/GT[r.VOI]-1)*100,axis=1)
S['NUKFIT_AD_Gy']=S.apply(lambda r:NK[r.VOI][FWS.index(r.FWHM_mm)],axis=1)
S['NUKFIT_RC']=S.apply(lambda r:NK[r.VOI][FWS.index(r.FWHM_mm)]/NK[r.VOI][0],axis=1)
S.to_csv(f'{ROOT}/results/summary_recommended_models.csv',index=False)
for v in GT:
    a=S[S.VOI==v]
    print(f'{v:12s} '+' | '.join(f'{fw:g}mm {ad:.3f}±{sd:.3f} RC {rc:.3f} (NK {nrc:.3f}) RE {re:+.1f}%' for fw,ad,sd,rc,nrc,re in zip(a.FWHM_mm,a.AD_Gy,a.AD_SD_Gy,a.RC_vs_0mm,a.NUKFIT_RC,a.RE_vs_GT_pct)))
print('Teff by FWHM:'); print(S.pivot(index='VOI',columns='FWHM_mm',values='Teff_72_144_h').round(0))
