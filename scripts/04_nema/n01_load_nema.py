"""Load the NEMA IEC SPECT (Xeleris, OSEM 4x10, RR) and CT of the public dataset; sphere centroids on SPECT.
Inputs : input/GE_TEW_RR_4X10.dcm, input/nema_ct/*.dcm (CT.zip of the dataset)
Outputs: results/nema/spect.npy, cents.npy, ct.npy, ct_z.npy
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results/nema')
import pydicom, numpy as np
from scipy import ndimage as ndi
d=pydicom.dcmread(f'{ROOT}/input/GE_TEW_RR_4X10.dcm'); a=d.pixel_array.astype(float)
lab,n=ndi.label(a>0.15*a.max()); sizes=ndi.sum(np.ones_like(a),lab,range(1,n+1))
print('components',n,'sizes',sorted(sizes.astype(int),reverse=True)[:10])
VX=4.419644; x0,y0,z0=-279.947394,-280.047394,280.647394
order=np.argsort(sizes)[::-1][:6]
cents=[]
for k in order:
    m=lab==(k+1); c=np.array(ndi.center_of_mass(a*m))
    cents.append(c); print(f'blob {k+1}: voxels {int(sizes[k])}, max {a[m].max():.0f}, centre (i,r,c) {np.round(c,2)} -> (x,y,z) mm ({x0+c[2]*VX:.1f}, {y0+c[1]*VX:.1f}, {z0-c[0]*VX:.1f})')
# total counts
print('total counts FOV %.4e'%a.sum())
# decay-corrected concentration and total activity
lam=np.log(2)/(6.647*24); dt=(16+49/60+16/3600)-(10+50/60)
C=2.051*np.exp(-lam*dt); vols=np.pi/6*np.array([10,13,17,22,28,37])**3/1000
print('dt %.3f h, conc at scan %.4f MBq/mL, sphere volumes %s, total %.3f mL, total activity %.2f MBq'%(dt,C,np.round(vols,4),vols.sum(),C*vols.sum()))
print('implied sensitivity with 3600 s: %.2f cps/MBq; with 1800 s: %.2f'%(a.sum()/(3600*C*vols.sum()),a.sum()/(1800*C*vols.sum())))
np.save('spect.npy',a); np.save('cents.npy',np.array(cents))
# ---- CT
import pydicom, numpy as np, glob
fs=glob.glob(f'{ROOT}/input/nema_ct/*.dcm')
sl=[pydicom.dcmread(f) for f in fs]
sl=[s for s in sl if hasattr(s,'ImagePositionPatient')]
sl.sort(key=lambda s: float(s.ImagePositionPatient[2]))
z=np.array([float(s.ImagePositionPatient[2]) for s in sl]); print('n slices',len(sl),'z range',z.min(),z.max(),'dz',np.round(np.diff(z)[:3],3))
vol=np.stack([s.pixel_array*float(s.RescaleSlope)+float(s.RescaleIntercept) for s in sl])
print('CT volume',vol.shape,'HU range',vol.min(),vol.max())
np.save('ct.npy',vol.astype(np.float32)); np.save('ct_z.npy',z)
ox,oy=[float(v) for v in sl[0].ImagePositionPatient[:2]]; ps=float(sl[0].PixelSpacing[0]); print('origin',ox,oy,'ps',ps)
# slice nearest z=32.5
k=int(np.argmin(abs(z-32.5))); img=vol[k]
# profile across the largest sphere centre (x=27.4, y=51.4)
cx=int(round((27.4-ox)/ps)); cy=int(round((51.4-oy)/ps))
print('slice z',z[k],' profile HU along x through sphere 37 mm:',np.round(img[cy,cx-25:cx+26:2]).astype(int))
