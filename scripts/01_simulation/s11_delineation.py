"""Split of the VOI term into voxelisation (GT-compartment VOIs) and the VOI-compartment difference (CT VOIs).
Inputs : GT maps, derived/frac_tp_avg.npy, derived/seglab.npy
Outputs: results/delineation.json
Manuscript: Table S5, Fig. S2.
"""
from resdos.paths import ROOT
import numpy as np, json, pydicom, importlib.util
import resdos.fit_tac as ft
T=np.array([1,4,24,40,72,144.]); TP=['1h','4h','24h','40h','72h','144h']
KT=1916.5/(32000*float(pydicom.dcmread(f'{ROOT}/input/truth_1h_reg.dcm',stop_before_pixels=True).RescaleSlope))
print('kBq/mL per unit from 11.94 cps/MBq, 1800 s: %.4f | from Table 1: %.4f'%(1000/(11.94*1800)/(0.4419644**3),KT))
F=np.load(f'{ROOT}/derived/frac_tp_avg.npy')
grp=np.stack([F[...,0],F[...,1],F[...,2],F[...,3],F[...,4],F[...,5]+F[...,6],F[...,7]+F[...,8],F[...,9]],-1)
labgt=(grp+1e-3*np.arange(8)).argmax(-1)   # GT-compartment majority VOIs: 3 liver,4 spleen,5 kidL,6 kidR,7 tumour
labct=np.load(f'{ROOT}/derived/seglab.npy')    # CT (3D Slicer) VOIs: 1 liver,2 tumour,3 spleen,4 R kidney,5 L kidney
ims=[]
for t in TP:
    d=pydicom.dcmread(f'{ROOT}/input/truth_{t}_reg.dcm'); ims.append(d.pixel_array*float(d.RescaleSlope))
tab={'Tumour':[1916.5,1809.2,1864.5,1713.2,1396.6,825.9],'spleen':[374.5,428.0,500.4,433.7,314.4,151.7],'Right Kidney':[494.9,466.2,403.5,326.1,214.4,82.7],
     'Left Kidney':[490.0,456.5,404.2,325.1,212.4,82.4],'Liver':[113.5,162.4,137.9,111.2,84.2,46.8]}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
mp={'Tumour':(7,2),'spleen':(4,3),'Right Kidney':(6,4),'Left Kidney':(5,5),'Liver':(3,1)}
I=lambda m,y: ft.integral(m,ft.fit(m,T,np.asarray(y,float))[0])
out={}
for v,(lg,lc) in mp.items():
    ytrue=np.array(tab[v])/KT; It=I(rec[v],ytrue)
    ygt=[im[labgt==lg].mean() for im in ims]; yct=[im[labct==lc].mean() for im in ims]
    vg=I(rec[v],ygt)/It-1; vc=I(rec[v],yct)/It-1
    m1=labgt==lg; m2=labct==lc; dice=2*(m1&m2).sum()/(m1.sum()+m2.sum())
    out[v]=dict(vox_gtvoi=vg,vox_ctvoi=vc,delineation=vc-vg,dice=float(dice))
    print(f'{v:12s} VOI term: GT-compartment VOI {vg*100:+5.1f}% | CT (3D Slicer) VOI {vc*100:+5.1f}% | delineation part {(vc-vg)*100:+4.1f} pt | Dice {dice:.3f}')
json.dump(out,open(f'{ROOT}/results/delineation.json','w'),indent=1)
