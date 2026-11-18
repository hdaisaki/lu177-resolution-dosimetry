"""Global effective Gaussian FWHM by least-squares fit of the blurred GT to each SPECT image.
Outputs: results/spect/fwhm_fit.json
Manuscript: Results 3.5.
"""
from resdos.paths import ROOT
import sys; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.ndimage import binary_dilation
import json
S=np.load(f'{ROOT}/results/spect/spect.npy')
res={'cal':{},'fwhm':{}}
FG=[8,10,12,13,14,15,16,17,18,20,22,25]
for i,t in enumerate(TS):
    tr=truth(t)
    body=binary_dilation(to_spect(blur(tr,10.0),order=3)>0.02*tr.max()*0.05,iterations=3)
    # calibration on total activity in a generous body mask
    tot_true=tr.sum(); k_tot=S[i][body].sum()/tot_true
    rss=[]
    for f in FG:
        m=to_spect(blur(tr,float(f)),order=3)
        k=(S[i][body]*m[body]).sum()/(m[body]**2).sum()
        r=S[i][body]-k*m[body]; rss.append(((r**2).sum(),k))
    j=int(np.argmin([x[0] for x in rss]))
    # parabolic refinement
    if 0<j<len(FG)-1:
        x=np.array(FG[j-1:j+2],float); y=np.array([rss[j-1][0],rss[j][0],rss[j+1][0]])
        p=np.polyfit(x,y,2); fb=-p[1]/(2*p[0])
    else: fb=FG[j]
    res['cal'][t]=k_tot; res['fwhm'][t]=dict(best=float(fb),k_fit=rss[j][1])
    print(f'{t:4s}: counts per truth-unit (total) {k_tot:.4f} | best-fit FWHM {fb:.1f} mm (k {rss[j][1]:.4f}) | rel RSS at 10/15/20 mm: '+' '.join(f'{rss[FG.index(f)][0]/rss[j][0]:.2f}' for f in (10,15,20)))
json.dump(res,open(f'{ROOT}/results/spect/fwhm_fit.json','w'),indent=1)
