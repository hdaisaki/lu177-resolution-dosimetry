"""Resample the 177Lu soft-tissue voxel S-value kernel (DOODLE kernel, 4.7952 mm) to the 4.42-mm grid,
conserving the electron energy deposited near the source voxel.
Inputs : derived/kernels/DK_softtissueICRP_mGyperdecay.img, the unmodified kernel of PyTheranostics 0.1 (Uribe C et al.,
         https://github.com/qurit/PyTheranostics; Apache License 2.0), included in the derived-data record; if missing,
         the original wheel is downloaded (see resdos/kernel_source.py and THIRD_PARTY.md). SHA-256 checked.
Outputs: derived/vsv_lu177_4.42mm.npy
Manuscript: Supplementary Methods S3.
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/derived')
from resdos.kernel_source import ensure_kernel
KERNEL = ensure_kernel(os.path.join(ROOT, 'derived'))
import numpy as np
L_old=4.7952; L_new=4.419644
k=np.fromfile(KERNEL,dtype='<f4').reshape(51,51,51).astype(float)*1e6  # mGy/(MBq s)
c=25
idx=np.indices(k.shape)-c
ri=np.sqrt((idx**2).sum(0))          # radius in voxels (old)
cheb=np.abs(idx).max(0)
# photon/long-range part: voxels not touching the source voxel (cheb>=2)
far=cheb>=2
r_mm=ri*L_old
rs=np.unique(np.round(r_mm[far],3))
g=np.array([k[far & (np.round(r_mm,3)==r)].mean() for r in rs])
# smooth radial function: interpolate log(r^2 S) vs r
y=np.log(rs**2*g)
def S_ph(r):
    r=np.asarray(r,float); rr=np.clip(r,rs[0],rs[-1])
    out=np.exp(np.interp(rr,rs,y))/np.maximum(r,1e-9)**2
    # extrapolate inward for r<rs[0] using r^2S constant
    return out
# near field (center, face, edge, corner) electron energies (proportional to S*L^3)
near_types={'c':(0,0,0),'f':(0,0,1),'e':(0,1,1),'k':(1,1,1)}
Sold={t:k[c+v[0],c+v[1],c+v[2]] for t,v in near_types.items()}
rad={t:np.sqrt(sum(np.array(v)**2)) for t,v in near_types.items()}
mult={'c':1,'f':6,'e':12,'k':8}
Ee={t:(Sold[t]-(S_ph(rad[t]*L_old) if t!='c' else 0.0))*L_old**3 for t in near_types}
Etot=sum(mult[t]*Ee[t] for t in Ee)
esc_old=1-Ee['c']/Etot
kk=L_old/L_new
esc_new=esc_old*kk
print('electron escape fraction old %.4f -> new %.4f'%(esc_old,esc_new))
knew=S_ph(ri*L_new)  # long-range for all
for t,v in near_types.items():
    m=(np.sort(np.abs(idx),axis=0)==np.array(sorted(v))[:,None,None,None]).all(0)
    if t=='c': Ee_new=Etot*(1-esc_new)
    else: Ee_new=Ee[t]*kk
    knew[m]=Ee_new/L_new**3 + (S_ph(rad[t]*L_new) if t!='c' else 0.0)
np.save('vsv_lu177_4.42mm.npy',knew)
LED=lambda L:0.1484*1.602e-13*1e6/((L/10)**3*1e-3)*1e3
print('old center %.4f (LED %.4f) | new center %.4f (LED %.4f)'%(k[c,c,c],LED(L_old),knew[c,c,c],LED(L_new)))
print('energy within box  old %.4f  new %.4f  (S*L^3 sum, mGy/MBqs*cm3)'%((k*L_old**3).sum()/1e3,(knew*L_new**3).sum()/1e3))
print('face old %.5f new %.5f ; r=2 old %.3e new %.3e'%(k[c,c,c+1],knew[c,c,c+1],k[c,c,c+2],knew[c,c,c+2]))
