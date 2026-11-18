from resdos.paths import ROOT
import numpy as np, pydicom
from scipy.ndimage import map_coordinates, convolve1d
from scipy.special import ive
VX=4.419644; TS=['1h','4h','24h','40h','72h','144h']; T=np.array([1,4,24,40,72,144.])
# truth grid origin vs SPECT grid origin (physical, LPS)
T0=np.array([284.897394,-279.047394,-279.947394])   # z(first frame), y, x  (truth_reg)
S0=np.array([280.647394,-280.047394,-279.947394])   # SPECT
# truth index of SPECT voxel centre:  i_t = i_s + (T0z - S0z)/VX ; r_t = r_s + (S0y - T0y)/VX ; c_t = c_s
DZ=(T0[0]-S0[0])/VX; DY=(S0[1]-T0[1])/VX; DX=(S0[2]-T0[2])/VX
def truth(t,fw=0.0):
    p=f'{ROOT}/input/truth_{t}_reg.dcm' if fw==0 else f'{ROOT}/results/blurred_dicom/FWHM_{fw:04.1f}mm/truth_{t}_reg_G{fw:g}mm.dcm'
    d=pydicom.dcmread(p); return d.pixel_array.astype(np.float64)*float(d.RescaleSlope)
def g1(s):
    t=s*s; R=int(np.ceil(5*s))+3; n=np.arange(-R,R+1); w=ive(np.abs(n),t); return w/w.sum()
def blur(img,fw):
    if fw<=0: return img
    out=img
    for ax in range(3): out=convolve1d(out,g1(fw/2.3548/VX),axis=ax,mode='constant')
    return out
N=128
_I,_R,_C=np.meshgrid(np.arange(N),np.arange(N),np.arange(N),indexing='ij')
def to_spect(vol,shift=(0,0,0),order=1):
    """sample a truth-grid volume at SPECT voxel centres (+optional extra shift in voxels)"""
    coords=np.array([_I+DZ+shift[0],_R+DY+shift[1],_C+DX+shift[2]])
    return map_coordinates(vol,coords,order=order,mode='constant',cval=0.0)
