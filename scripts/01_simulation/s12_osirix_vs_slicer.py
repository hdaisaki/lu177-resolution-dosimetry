"""Agreement between the OsiriX GT segmentation and the 3D Slicer VOIs (Dice, mean surface distance, HD95, volumes).
Inputs : derived/osirix_masks_+0.5.npz, derived/seg.npy
Outputs: results/osirix_vs_slicer.json
Manuscript: Methods 2.1, Table S5.
"""
from resdos.paths import ROOT
import numpy as np, json
from scipy import ndimage as ndi
O=np.load(f'{ROOT}/derived/osirix_masks_+0.5.npz'); seg=np.load(f'{ROOT}/derived/seg.npy')[:124]
sp=np.array([2.5,0.976562,0.976562]); vv=np.prod(sp)/1000; tum=O['tumor']
osx={'Liver':O['liver']&~tum,'Tumour':tum,'Spleen':O['spleen'],'Right kidney':O['R_kidney'],'Left kidney':O['L_kidney']}
slc={'Liver':seg==1,'Tumour':seg==2,'Spleen':seg==3,'Right kidney':seg==4,'Left kidney':seg==5}
GTm={'Liver':1306.7,'Tumour':15.9,'Spleen':124.9,'Right kidney':111.4,'Left kidney':133.8}
res={}
for k in osx:
    a,b=osx[k],slc[k]
    sa=a&~ndi.binary_erosion(a); sb=b&~ndi.binary_erosion(b)
    d=np.r_[ndi.distance_transform_edt(~sb,sampling=sp)[sa],ndi.distance_transform_edt(~sa,sampling=sp)[sb]]
    dice=2*(a&b).sum()/(a.sum()+b.sum())
    res[k]=dict(dice=float(dice),msd=float(d.mean()),hd95=float(np.percentile(d,95)),vol_osirix=float(a.sum()*vv),vol_slicer=float(b.sum()*vv))
    print(f"{k:12s} Dice {dice:.3f} | MSD {d.mean():.2f} mm | HD95 {np.percentile(d,95):.2f} mm | volume OsiriX {a.sum()*vv:7.1f} / Slicer {b.sum()*vv:7.1f} mL (dataset mass {GTm[k]} g)")
json.dump(res,open(f'{ROOT}/results/osirix_vs_slicer.json','w'),indent=1)
