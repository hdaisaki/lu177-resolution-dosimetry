"""CT-based VOIs on the SPECT grid (occupancy and binary >= 50%).
Outputs: results/spect/segfrac_spect.npy, results/spect/seglab_spect.npy
"""
from resdos.paths import ROOT
import numpy as np
seg=np.load(f'{ROOT}/derived/seg.npy'); sxy=0.976562; sz=2.5; ox=oy=-250.0; oz=-167.0
VX=4.419644; x0,y0,z0=-279.947394,-280.047394,280.647394   # SPECT grid
N=128; S=3; offs=(np.arange(S)+0.5)/S-0.5
frac=np.zeros((6,N,N,N),np.float32)
zz,yy,xx=np.nonzero(seg)
c0=int((ox+xx.min()*sxy-x0)/VX)-2; c1=int((ox+xx.max()*sxy-x0)/VX)+3
r0=int((oy+yy.min()*sy if False else oy+yy.min()*sxy - y0)/VX)-2; r1=int((oy+yy.max()*sxy-y0)/VX)+3
i0=int((z0-(oz+zz.max()*sz))/VX)-2; i1=int((z0-(oz+zz.min()*sz))/VX)+3
I,R,C=np.meshgrid(np.arange(i0,i1),np.arange(r0,r1),np.arange(c0,c1),indexing='ij')
for a in offs:
  for b in offs:
    for c in offs:
      X=x0+(C+c)*VX; Y=y0+(R+b)*VX; Z=z0-(I+a)*VX
      ix=np.round((X-ox)/sxy).astype(int); iy=np.round((Y-oy)/sxy).astype(int); iz=np.round((Z-oz)/sz).astype(int)
      ok=(ix>=0)&(ix<512)&(iy>=0)&(iy<512)&(iz>=0)&(iz<seg.shape[0])
      lab=np.zeros(ix.shape,np.uint8); lab[ok]=seg[iz[ok],iy[ok],ix[ok]]
      for l in range(6): frac[l,i0:i1,r0:r1,c0:c1]+=(lab==l)
frac/=S**3; frac[0]=1-frac[1:].sum(0)
lab=frac.argmax(0).astype(np.uint8)
np.save(f'{ROOT}/results/spect/segfrac_spect.npy',frac); np.save(f'{ROOT}/results/spect/seglab_spect.npy',lab)
print({l:int((lab==l).sum()) for l in range(1,6)})
