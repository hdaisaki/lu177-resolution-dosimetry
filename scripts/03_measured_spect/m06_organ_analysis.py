"""TIA recovery Q of measured SPECT and of the Gaussian model on the SPECT grid. Usage: python m06_organ_analysis.py k_abs
Outputs: results/spect/organ_results_k_abs.json
"""
from resdos.paths import ROOT
import sys, json, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
import resdos.fit_tac as ft
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy')
cal=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['cal']
kt=np.array([cal[t] for t in TS]); kg=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))[sys.argv[1]]
V={'Tumour':(2,'f4b'),'spleen':(3,'f4b'),'Right Kidney':(4,'f2'),'Left Kidney':(5,'f2'),'Liver':(1,'f3')}
# true concentration TACs (image units)
C={'Liver':[1895,2872,2367,2077,1929,1813],'Spleen':[6253,7570,8588,8101,7204,5878],'KidLC':[10386,10262,8789,7733,6193,4068],'KidLM':[3511,3435,3012,2559,2053,1341],
   'KidRC':[10556,10519,8849,7772,6267,4080],'KidRM':[3410,3431,2849,2533,2044,1344],'Tumor':[32000]*6}
F=np.load(f'{ROOT}/derived/frac_tp_avg.npy'); CN=['Air','BG','Lung','Liver','Spleen','KidLC','KidLM','KidRC','KidRM','Tumor']
sl=np.array([float(pydicom.dcmread(f'{ROOT}/input/truth_{t}_reg.dcm',stop_before_pixels=True).RescaleSlope) for t in TS])
parts={'Tumour':['Tumor'],'spleen':['Spleen'],'Right Kidney':['KidRC','KidRM'],'Left Kidney':['KidLC','KidLM'],'Liver':['Liver']}
ytrue={}
for v,p in parts.items():
    w=np.array([F[...,CN.index(q)].sum() for q in p]); ytrue[v]=sum(wi*np.array(C[q])*sl for wi,q in zip(w,p))/w.sum()
def tiac(m,y): p=ft.fit(m,T,np.asarray(y,float))[0]; return ft.integral(m,p)
out={'kt':kt.tolist(),'kg':float(kg)}
# model curves on SPECT grid
FWM=[5.0,7.5,9.5,10.0,12.5,15.0,17.5,20.0]
tr=[truth(t) for t in TS]
model={v:{} for v in V}
for f in FWM:
    ims=[to_spect(blur(x,f),order=3) for x in tr]
    for v,(l,m) in V.items():
        y=[im[lab==l].mean() for im in ims]; model[v][f]=y
act=json.load(open(f'{ROOT}/results/activity_tacs.json'))
for v,(l,m) in V.items():
    model[v][0.0]=act[v]['0.0']   # truth grid, same CT VOI rule
spect_pt=[]; 
res={}
for v,(l,m) in V.items():
    yS=np.array([S[i][lab==l].mean() for i in range(6)])
    yS_t=yS/kt; yS_g=yS/kg
    It=tiac(m,ytrue[v])
    q={f:tiac(m,model[v][f])/It for f in sorted(model[v])}
    qs_t=tiac(m,yS_t)/It; qs_g=tiac(m,yS_g)/It
    # equivalent FWHM by interpolation of Q_model(F) (monotone decreasing)
    fs=np.array(sorted(q)); qq=np.array([q[f] for f in fs])
    feq=float(np.interp(-qs_g,-qq,fs)) if qq.min()<=qs_g<=qq.max() else None
    res[v]=dict(Q_model={str(k):float(x) for k,x in q.items()},Q_spect_total=float(qs_t),Q_spect=float(qs_g),FWHM_eq=feq,
                tac_spect=yS_g.tolist(),tac_true=ytrue[v].tolist(),tac_ratio=(yS_g/ytrue[v]).tolist())
    print(f'{v:12s} Q(model) 0mm {q[0.0]:.3f} 5 {q[5.0]:.3f} 9.5 {q[9.5]:.3f} 10 {q[10.0]:.3f} 15 {q[15.0]:.3f} 20 {q[20.0]:.3f} | SPECT Q (const cal) {qs_g:.3f} (total-count cal {qs_t:.3f}) | FWHM_eq {feq if feq is None else round(feq,1)} mm')
    print('      SPECT/true per time point:',np.round(yS_g/ytrue[v],3),' model 9.5 mm:',np.round(np.array(model[v][9.5])/ytrue[v],3))
out['organs']=res
json.dump(out,open(f'{ROOT}/results/spect/organ_results_{sys.argv[1]}.json','w'),indent=1)
