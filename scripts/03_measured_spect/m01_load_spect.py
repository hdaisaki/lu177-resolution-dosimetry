"""Load the six reconstructed SPECT images of the public dataset.
Inputs : input/{t}_SPECT_GE.dcm (MRTDosimetry, OSF)
Outputs: results/spect/spect.npy
"""
from resdos.paths import ROOT
import pydicom, numpy as np
d=pydicom.dcmread(f'{ROOT}/input/1h_SPECT_GE.dcm')
print(repr(d.ImageComments))
ts=['1h','4h','24h','40h','72h','144h']
S=np.stack([pydicom.dcmread(f'{ROOT}/input/{t}_SPECT_GE.dcm').pixel_array.astype(np.float32)*float(getattr(pydicom.dcmread(f'{ROOT}/input/{t}_SPECT_GE.dcm',stop_before_pixels=True),'RescaleSlope',1)) for t in ts])
np.save(f'{ROOT}/results/spect/spect.npy',S)
print('min',S.min(),'neg voxels',(S<0).sum())
