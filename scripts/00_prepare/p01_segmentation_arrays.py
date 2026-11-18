"""Rasterise the 3D Slicer segmentation (liver, tumour, spleen, kidneys) and map it to the 4.42-mm GT grid.
Inputs : input/MRT_Segmentation.nrrd (3D Slicer, CT grid)
Outputs: derived/seg.npy (CT grid labels); derived/segfrac.npy (occupancy on the GT grid); derived/seglab.npy (binary VOIs, >= 50%)
Manuscript: Methods 2.1 (VOIs); used by all VOI analyses.
Requires SimpleITK.
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/derived')
import SimpleITK as sitk, numpy as np
im=sitk.ReadImage(f'{ROOT}/input/MRT_Segmentation.nrrd')
a=sitk.GetArrayFromImage(im)  # z,y,x
print('size',im.GetSize(),'spacing',im.GetSpacing(),'origin',im.GetOrigin(),'dir',im.GetDirection())
r=sitk.ImageFileReader(); r.SetFileName(f'{ROOT}/input/MRT_Segmentation.nrrd'); r.ReadImageInformation()
names={}
for k in r.GetMetaDataKeys():
    if k.endswith('_Name') or k.endswith('_LabelValue'): names[k]=r.GetMetaData(k)
for i in range(10):
    if f'Segment{i}_Name' in names: print(i,names[f'Segment{i}_Name'],'label',names[f'Segment{i}_LabelValue'])
vv=np.prod(im.GetSpacing())/1000
u,c=np.unique(a,return_counts=True)
for uu,cc in zip(u,c): print('label',uu,'voxels',cc,'vol %.1f mL'%(cc*vv))
np.save('seg.npy',a)

# ---- map to the GT grid (formerly segmap.py)
import numpy as np
seg=np.load('seg.npy')            # (z,y,x) CT grid, LPS
sx=sy=0.976562; sz=2.5; ox,oy,oz=-250.0,-250.0,-167.0
VX=4.419644; x0,y0,z0=-279.947394,-279.047394,284.897394; dzs=-4.41964
N=128; S=3
offs=(np.arange(S)+0.5)/S-0.5
frac=np.zeros((6,N,N,N),np.float32)
# bounding box of labels in truth index space
zz,yy,xx=np.nonzero(seg)
xmin,xmax=ox+xx.min()*sx,ox+xx.max()*sx; ymin,ymax=oy+yy.min()*sy,oy+yy.max()*sy; zmin,zmax=oz+zz.min()*sz,oz+zz.max()*sz
c0,c1=int((xmin-x0)/VX)-2,int((xmax-x0)/VX)+3; r0,r1=int((ymin-y0)/VX)-2,int((ymax-y0)/VX)+3
i0,i1=int((z0-zmax)/-dzs)-2,int((z0-zmin)/-dzs)+3
print('truth index bbox z',i0,i1,'y',r0,r1,'x',c0,c1)
I,R,Cc=np.meshgrid(np.arange(i0,i1),np.arange(r0,r1),np.arange(c0,c1),indexing='ij')
for a in offs:
    for b in offs:
        for c in offs:
            X=x0+(Cc+c)*VX; Y=y0+(R+b)*VX; Z=z0+(I+a)*dzs
            ix=np.round((X-ox)/sx).astype(int); iy=np.round((Y-oy)/sy).astype(int); iz=np.round((Z-oz)/sz).astype(int)
            ok=(ix>=0)&(ix<512)&(iy>=0)&(iy<512)&(iz>=0)&(iz<seg.shape[0])
            lab=np.zeros(ix.shape,np.uint8); lab[ok]=seg[iz[ok],iy[ok],ix[ok]]
            for l in range(6):
                frac[l,i0:i1,r0:r1,c0:c1]+=(lab==l)
frac/=S**3
frac[0]=1-frac[1:].sum(0)
np.save('segfrac.npy',frac)
VV=(VX/10)**3
names=['bg','Liver','Tumour','Spleen','Right Kidney','Left Kidney']
lab=frac.argmax(0)
np.save('seglab.npy',lab.astype(np.uint8))
for l in range(1,6): print(f'{names[l]:12s} frac-sum {frac[l].sum()*VV:7.1f} mL   majority VOI {(lab==l).sum()*VV:7.1f} mL')
