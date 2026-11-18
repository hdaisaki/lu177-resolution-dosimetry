"""Three-stage evaluation of uncorrected, GTM-corrected and RC-corrected data under identical conditions:
concentration per time point (vs Table 1), time-integrated concentration (vs SAAM II) and proportional AD estimate
(common OpenDose3D dose-rate-to-concentration ratio).
Inputs : results/spect/unified_{8.0,9.5,12.0}.json, results/gt_auc_validation.json, results/spect/gtm_final.json
Outputs: results/unified_eval.json, results/unified_stage3.json, results/pvc_summary.json
Manuscript: Table 5, Fig. 6c, Tables S8, S11.
"""
from resdos.paths import ROOT
import json, numpy as np
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); I=lambda m,y: ft.integral(m,ft.fit(m,T,np.asarray(y,float))[0])
tab={'Tumour':[1916.4519,1809.2385,1864.4781,1713.1579,1396.5609,825.8863],'spleen':[374.4555,428.0064,500.3910,433.7019,314.4332,151.7456],
     'Right Kidney':[494.9367,466.2140,403.4738,326.1331,214.4026,82.6633],'Left Kidney':[489.9914,456.4592,404.1572,325.1385,212.3693,82.3926],
     'Liver':[113.4567,162.4198,137.9349,111.2321,84.2231,46.8034]}
GT={'Tumour':dict(V=15.901,AUC=5267.948,AD=28.12285),'spleen':dict(V=124.948,AUC=7814.396,AD=5.382068),'Right Kidney':dict(V=111.4,AUC=4704.988,AD=3.598351),
    'Left Kidney':dict(V=133.768,AUC=5649.918,AD=3.599789),'Liver':dict(V=1306.698,AUC=24937.645,AD=1.702361)}
VOL={'Tumour':15.8,'spleen':121.8,'Right Kidney':110.4,'Left Kidney':133.3,'Liver':1291.9}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
R=lambda V: 1/(1+(2.53/V)**0.41)          # NEMA RC curve (OpenDose3D tool, volume threshold)
Rct=lambda V: 1/(1+(2.16/V)**0.39)        # NEMA RC curve, CT-positioned spheres
U={fw:json.load(open(f'{ROOT}/results/spect/unified_{fw}.json')) for fw in ['8.0','9.5','12.0']}
Vv=json.load(open(f'{ROOT}/results/gt_auc_validation.json')); G=json.load(open(f'{ROOT}/results/spect/gtm_final.json'))
E={}; S3={}; P={}
for v in tab:
    u=U['9.5'][v]; tru=np.array(tab[v]); cu=np.array(u['conc_unc']); cg=np.array(u['conc_gtm']); r=R(VOL[v]); cr=cu/r
    Itrue=GT[v]['AUC']/GT[v]['V']*1000
    tu,tg,tr_=I(rec[v],cu),I(rec[v],cg),I(rec[v],cr)
    au=I(rec[v],u['dr_unc'])/1000; ag=I(rec[v],u['dr_gtm'])/1000; ac=I(rec[v],u['dr_gtm_comp'])/1000; ar=au/r
    E[v]=dict(RC=r,conc_err=dict(unc=((cu/tru-1)*100).tolist(),gtm=((cg/tru-1)*100).tolist(),rc=((cr/tru-1)*100).tolist()),
              tiac_err=dict(unc=(tu/Itrue-1)*100,gtm=(tg/Itrue-1)*100,rc=(tr_/Itrue-1)*100),
              ad_err=dict(unc=(au/GT[v]['AD']-1)*100,gtm=(ag/GT[v]['AD']-1)*100,rc=(ar/GT[v]['AD']-1)*100,gtm_comp=(ac/GT[v]['AD']-1)*100),
              scale_check=tr_/(tu/r),self_dose_fraction=float(np.mean(np.array(u['dr_self'])/np.array(u['dr_unc']))))
    conv=Vv[v]['dose_conversion']; t=E[v]['tiac_err']
    ad=lambda e: ((1+e/100)*(1+conv)-1)*100
    g8,g12=[(I(rec[v],U[fw][v]['conc_gtm'])/Itrue-1)*100 for fw in ('8.0','12.0')]
    S3[v]=dict(conv=conv*100,tiac=t,ad={m:ad(t[m]) for m in ('unc','gtm','rc')},gtm_tiac_psf=[g8,g12])
    tun=1+t['unc']/100; rc_ct_t=(tun/Rct(VOL[v])-1)*100
    ads=[ad(x) for x in (g8,t['gtm'],g12)]
    P[v]=dict(conv=conv*100,tiac=t,ad=S3[v]['ad'],rc_ct_t=rc_ct_t,rc_ct_ad=ad(rc_ct_t),gtm_t_rng=[g8,g12],gtm_ad_rng=[min(ads),max(ads)],
              gtm1_t=(G[v]['9.5']*(1+Vv[v]['integration'][rec[v]])-1)*100)
    print(f"{v:12s} TIAC unc {t['unc']:+6.1f} GTM {t['gtm']:+6.1f} RC {t['rc']:+6.1f} | AD GTM {S3[v]['ad']['gtm']:+6.1f} RC {S3[v]['ad']['rc']:+6.1f}")
json.dump(E,open(f'{ROOT}/results/unified_eval.json','w'),indent=1)
json.dump(S3,open(f'{ROOT}/results/unified_stage3.json','w'),indent=1)
json.dump(P,open(f'{ROOT}/results/pvc_summary.json','w'),indent=1)
