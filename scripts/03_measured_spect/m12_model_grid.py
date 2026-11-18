"""Blurred-GT VOI means on the SPECT grid for FWHM 6-14 mm. Usage: python m12_model_grid.py 0 1 2 3 4 5
Outputs: results/spect/model_grid_{t}.json
"""
from resdos.paths import ROOT
import sys, json; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy')
FG=[6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
V={'Liver':1,'Tumour':2,'spleen':3,'Right Kidney':4,'Left Kidney':5}
for ti in map(int,sys.argv[1:]):
    t=TS[ti]; tr=truth(t); out={}
    for f in FG:
        m=to_spect(blur(tr,f),order=3)
        out[str(f)]={v:float(m[lab==l].mean()) for v,l in V.items()}
    json.dump(out,open(f'{ROOT}/results/spect/model_grid_{t}.json','w'))
    print(t,'done')
