"""Stack the GT activity maps (six time points) for the occupancy estimation.
Inputs : input/truth_{t}_reg.dcm
Outputs: derived/raw.npy
"""
from resdos.paths import ROOT
import numpy as np, pydicom
TS=['1h','4h','24h','40h','72h','144h']
A=np.stack([pydicom.dcmread(f'{ROOT}/input/truth_{t}_reg.dcm').pixel_array.astype(np.float32) for t in TS])
np.save(f'{ROOT}/derived/raw.npy',A)
