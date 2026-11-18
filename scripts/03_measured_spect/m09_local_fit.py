"""Local FWHM and local quantitative scale per organ neighbourhood. Usage: python m09_local_fit.py <time index 0-5>
Outputs: results/spect/local_fit_{t}.json
"""
from resdos.paths import ROOT
import sys, json; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.ndimage import binary_dilation
ti=int(sys.argv[1]); t=TS[ti]
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy')
k_abs=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['k_abs']
regions={v:binary_dilation(lab==l,iterations=4) for v,l in [('Liver',1),('Tumour',2),('spleen',3),('Right Kidney',4),('Left Kidney',5)]}
regions['Liver']=regions['Liver']&~binary_dilation(lab==2,iterations=5)
FG=[5,6,7,8,9,10,11,12,13,14]
tr=truth(t); rss={v:[] for v in regions}
for f in FG:
    m=to_spect(blur(tr,float(f)),order=3)
    for v,reg in regions.items():
        y=S[ti][reg]; x=m[reg]; kk=(x*y).sum()/(x*x).sum(); rss[v].append(((y-kk*x)**2).sum()); rss[v][-1]=(rss[v][-1],kk)
out={}
for v in regions:
    r=np.array([a for a,b in rss[v]]); j=int(np.argmin(r))
    if 0<j<len(FG)-1:
        p=np.polyfit(FG[j-1:j+2],r[j-1:j+2],2); fb=float(-p[1]/(2*p[0]))
    else: fb=float(FG[j])
    out[v]=dict(fwhm=fb,scale=float(rss[v][j][1]/k_abs))
json.dump(out,open(f'{ROOT}/results/spect/local_fit_{t}.json','w'))
print(t,{v:(round(o['fwhm'],1),round(o['scale'],3)) for v,o in out.items()})
