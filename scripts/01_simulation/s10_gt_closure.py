"""VOI term of the 0-mm error from the images alone (article: "VOI term"; voxelisation plus VOI-compartment difference).
Inputs : derived/opendose3d/fw00, results/activity_tacs.json
Outputs: results/gt_closure.json
Manuscript: Table S5, Fig. S2.
"""
from resdos.paths import ROOT
import numpy as np, json, pandas as pd, glob, pydicom, importlib.util
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); lp=4.345e-3
KT=1916.5/(32000*float(pydicom.dcmread(f'{ROOT}/input/truth_1h_reg.dcm',stop_before_pixels=True).RescaleSlope))  # kBq/mL per image unit
tab={'Tumour':[1916.5,1809.2,1864.5,1713.2,1396.6,825.9],'spleen':[374.5,428.0,500.4,433.7,314.4,151.7],'Right Kidney':[494.9,466.2,403.5,326.1,214.4,82.7],
     'Left Kidney':[490.0,456.5,404.2,325.1,212.4,82.4],'Liver':[113.5,162.4,137.9,111.2,84.2,46.8]}
key=lambda s: float(s.split('_')[-1].replace('HR.tsv',''))
DR=pd.concat([pd.read_csv(f,sep='\t').set_index('Segment')['Mean'].rename(key(f)) for f in sorted(glob.glob(f'{ROOT}/derived/opendose3d/fw00/*.tsv'),key=key)],axis=1)
act=json.load(open(f'{ROOT}/results/activity_tacs.json'))
GT={'Tumour':28.1,'spleen':5.39,'Right Kidney':3.6,'Left Kidney':3.61,'Liver':1.70}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
OBS={'Tumour':-14.1,'spleen':-7.5,'Right Kidney':-11.8,'Left Kidney':-11.9,'Liver':-6.0}
def trap_tail(y):
    lam=max(-np.polyfit(T[-3:],np.log(y[-3:]),1)[0],lp); return 0.5*y[0]+np.trapezoid(y,T)+y[-1]/lam
def integ(m,y): return ft.integral(m,ft.fit(m,T,np.asarray(y,float))[0])
out={}
print('Integration of the GT (Table 1) concentration curves vs GT-implied TIAC (%):')
for v in tab:
    k=np.mean(DR.loc[v].values/np.array(act[v]['0.0']))
    ytrue=np.array(tab[v])/KT                     # image units
    tgt=GT[v]*1000/k
    e={m:integ(m,ytrue)/tgt-1 for m in ['f2','f3','f4b']}; e['ref']=trap_tail(ytrue)/tgt-1
    # voxelisation term: VOI TAC vs true TAC (time-integrated, same fixed function)
    vox=integ(rec[v],act[v]['0.0'])/integ(rec[v],ytrue)-1
    vox_ref=trap_tail(np.array(act[v]['0.0']))/trap_tail(ytrue)-1
    pred=((1+vox)*(1+e[rec[v]])-1)*100
    out[v]=dict(k=k,integration=e,voi_term=vox,voi_term_ref=vox_ref,pred=pred,obs=OBS[v])
    print(f"{v:12s} f2 {e['f2']*100:+5.1f} f3 {e['f3']*100:+5.1f} f4b {e['f4b']*100:+5.1f} ref {e['ref']*100:+5.1f} | voxelisation {vox*100:+5.1f} (ref {vox_ref*100:+5.1f}) | predicted RE0 {pred:+5.1f} vs observed {OBS[v]:+5.1f}")
json.dump(out,open(f'{ROOT}/results/gt_closure.json','w'),indent=1)
