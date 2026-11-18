"""Time integration of the Table 1 concentration curves against the published SAAM II AUC/volume, and the dose-conversion
difference between OpenDose3D (0 mm) and the GT.
Inputs : derived/opendose3d/fw00, results/activity_tacs.json
Outputs: results/gt_auc_validation.json
Manuscript: Results 3.1, Tables S1, S5, Fig. S2.
"""
from resdos.paths import ROOT
import numpy as np, json, pandas as pd, glob, importlib.util, pydicom
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); lp=4.345e-3
GT={'Tumour':dict(V=15.901,AUC=5267.94767803804,AD=28.12285014),'spleen':dict(V=124.948,AUC=7814.39634647587,AD=5.382067629),
    'Right Kidney':dict(V=111.4,AUC=4704.987759768963,AD=3.598350646),'Left Kidney':dict(V=133.768,AUC=5649.9177447367665,AD=3.599789479),
    'Liver':dict(V=1306.698,AUC=24937.6450682712,AD=1.702360556)}
tab={'Tumour':[1916.4519257596373,1809.2385388446949,1864.4780811467326,1713.1578719436197,1396.5608888720014,825.8863336149234],
     'spleen':[374.45548945516043,428.00639191272677,500.3909791332273,433.7018936548291,314.4331900422871,151.74564482079342],
     'Right Kidney':[494.93665904060174,466.21397870140015,403.47382124167785,326.1331220001673,214.40258445779244,82.6633472963747],
     'Left Kidney':[489.991394221854,456.45922754403887,404.15724040893633,325.13849437430326,212.3692621937584,82.39258618278674],
     'Liver':[113.45671940757998,162.41981709249492,137.93488624036078,111.23210386733963,84.22313036128853,46.80336440731882]}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
KT=1916.4519257596373/(32000*float(pydicom.dcmread(f'{ROOT}/input/truth_1h_reg.dcm',stop_before_pixels=True).RescaleSlope))
key=lambda s: float(s.split('_')[-1].replace('HR.tsv',''))
DR=pd.concat([pd.read_csv(f,sep='\t').set_index('Segment')['Mean'].rename(key(f)) for f in sorted(glob.glob(f'{ROOT}/derived/opendose3d/fw00/*.tsv'),key=key)],axis=1)
act=json.load(open(f'{ROOT}/results/activity_tacs.json'))
def trap_tail(y):
    lam=max(-np.polyfit(T[-3:],np.log(y[-3:]),1)[0],lp); return 0.5*y[0]+np.trapezoid(y,T)+y[-1]/lam
out={}
print('Time integration of the GT concentration curves vs SAAM II AUC/volume (independent of dose conversion), %')
for v,g in GT.items():
    y=np.array(tab[v]); I_GT=g['AUC']/g['V']*1000          # kBq*h/mL
    # sanity: activity at 1 h = conc * V
    e={m:ft.integral(m,ft.fit(m,T,y)[0])/I_GT-1 for m in ['f2','f3','f4b']}; e['model-free']=trap_tail(y)/I_GT-1
    # dose conversion: OpenDose3D dose rate per concentration (k, per truth unit) vs GT dose per TIAC
    k=np.mean(DR.loc[v].values/np.array(act[v]['0.0']))       # mGy/h per truth unit
    k_kbq=k/KT                                                # mGy/h per kBq/mL
    d_gt=g['AD']*1000/I_GT                                    # mGy per (kBq*h/mL)
    conv=k_kbq/d_gt-1
    out[v]=dict(I_GT=I_GT,integration={m:float(x) for m,x in e.items()},k_mGy_h_per_kBqml=k_kbq,dgt_mGy_per_kBqhml=d_gt,dose_conversion=conv)
    print(f"{v:12s} "+'  '.join(f"{m}{'*' if m==rec[v] else ''} {x*100:+5.1f}" for m,x in e.items())+f" | dose conversion OD3D/GT {conv*100:+5.1f}% (OD3D {k_kbq:.4f}, GT {d_gt:.4f} mGy per kBq h/mL)")
json.dump(out,open(f'{ROOT}/results/gt_auc_validation.json','w'),indent=1)
