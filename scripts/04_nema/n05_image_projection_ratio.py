"""Ratio of reconstructed image counts to acquired projection counts (DICOM Counts Accumulated) for the NEMA and
anthropomorphic images (calibration-independent check of the count scaling).
Outputs: results/nema/img_proj_ratio.json
Manuscript: Supplementary Methods S6.
"""
from resdos.paths import ROOT
import pydicom, numpy as np, json
def ratio(p):
    d=pydicom.dcmread(p); img=d.pixel_array.astype(float)*float(getattr(d,'RescaleSlope',1) or 1)+float(getattr(d,'RescaleIntercept',0) or 0)
    return img.sum()/float(d.CountsAccumulated)
nema=ratio(f'{ROOT}/input/GE_TEW_RR_4X10.dcm')
anth=[ratio(f'{ROOT}/input/{t}_SPECT_GE.dcm') for t in ['1h','4h','24h','40h','72h','144h']]
print(f'NEMA {nema:.2f}; anthropomorphic {np.round(anth,2)}; ratio {np.round(nema/np.array(anth),2)}')
json.dump(dict(nema=nema,anth=anth),open(f'{ROOT}/results/nema/img_proj_ratio.json','w'))
