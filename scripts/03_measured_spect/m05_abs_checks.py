"""Background, liver-interior and tumour-peak checks with the absolute calibration.
Outputs: results/spect/fwhm_fit.json (abs_bg, abs_lc, abs_tumour_peak)
Manuscript: Table S7, Results 3.5.
"""
from resdos.paths import ROOT
import sys, json; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.ndimage import binary_erosion, binary_dilation
from scipy import ndimage as ndi
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy')
d=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json')); k=d['k_abs']
tr0=to_spect(truth('24h'),order=0); body=tr0>0
bg=binary_erosion(body,iterations=4)&~binary_dilation(lab>0,iterations=4); lc=binary_erosion(lab==1,iterations=3)
out={'bg':[],'lc':[]}
for i,t in enumerate(TS):
    m=to_spect(blur(truth(t),9.5),order=3)
    out['bg'].append(float(S[i][bg].mean()/k/m[bg].mean())); out['lc'].append(float(S[i][lc].mean()/k/m[lc].mean()))
print('background (abs) measured/model:',np.round(out['bg'],3)); print('liver interior (abs):',np.round(out['lc'],3))
sl=float(pydicom.dcmread(f'{ROOT}/input/truth_24h_reg.dcm',stop_before_pixels=True).RescaleSlope); Ct=32000*sl
tm=lab==2; peaks=[]
for i,t in enumerate(TS):
    slt=float(pydicom.dcmread(f'{ROOT}/input/truth_{t}_reg.dcm',stop_before_pixels=True).RescaleSlope)
    peaks.append(float(S[i][tm].max()/k/(32000*slt)))
print('tumour peak / true (abs) per time point:',np.round(peaks,2))
d['abs_bg']=out['bg']; d['abs_lc']=out['lc']; d['abs_tumour_peak']=peaks; json.dump(d,open(f'{ROOT}/results/spect/fwhm_fit.json','w'),indent=1)
