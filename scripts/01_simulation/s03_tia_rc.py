"""TIA-based recovery coefficients with the CT-based VOIs from the (blurred) GT activity maps.
Inputs : input/truth_{t}_reg.dcm, results/blurred_dicom, derived/seglab.npy
Outputs: results/tia_rc.json, results/activity_tacs.json
Manuscript: Table 3, Fig. 3.
"""
from resdos.paths import ROOT
import os; os.chdir(ROOT)
import numpy as np, pydicom, json, importlib.util, pandas as pd
import resdos.fit_tac as ft
TP=['1h','4h','24h','40h','72h','144h']; T=np.array([1,4,24,40,72,144.])
FWS=[0.0,5.0,7.5,10.0,12.5,15.0]
lab=np.load('derived/seglab.npy')
L={'Tumour':(2,'f4b'),'spleen':(3,'f4b'),'Right Kidney':(4,'f2'),'Left Kidney':(5,'f2'),'Liver':(1,'f3')}
def img(t,fw):
    p=f'{ROOT}/input/truth_{t}_reg.dcm' if fw==0 else f'{ROOT}/results/blurred_dicom/FWHM_{fw:04.1f}mm/truth_{t}_reg_G{fw:g}mm.dcm'
    d=pydicom.dcmread(p); return d.pixel_array*float(d.RescaleSlope)
res={}
for fw in FWS:
    ims=[img(t,fw) for t in TP]
    for v,(l,m) in L.items():
        y=np.array([im[lab==l].mean() for im in ims])
        p,rss,aicc,cv,ab,adcv=ft.fit(m,T,y)
        res.setdefault(v,{})[fw]=dict(TIAC=ft.integral(m,p),cv=adcv,teff=72*np.log(2)/np.log(ft.MODELS[m][0](72.,*p)/ft.MODELS[m][0](144.,*p)),tac=y.tolist())
S=pd.read_csv(f'{ROOT}/results/summary_recommended_models.csv')
out={}
print('VOI          FWHM  RC_TIA  RC_AD   diff')
for v in L:
    t0=res[v][0.0]['TIAC']
    for fw in FWS:
        rct=res[v][fw]['TIAC']/t0
        rca=float(S[(S.VOI==v)&(S.FWHM_mm==fw)].RC_vs_0mm.iloc[0])
        out.setdefault(v,{})[fw]=dict(RC_TIA=rct,RC_AD=rca,TIA_CV=res[v][fw]['cv'],Teff=res[v][fw]['teff'])
        print(f'{v:12s} {fw:5.1f}  {rct:.4f}  {rca:.4f}  {rct-rca:+.4f}   Teff {res[v][fw]["teff"]:.0f} h')
json.dump(out,open('results/tia_rc.json','w'),indent=1)
json.dump({v:{str(f):res[v][f]['tac'] for f in FWS} for v in L},open('results/activity_tacs.json','w'))
