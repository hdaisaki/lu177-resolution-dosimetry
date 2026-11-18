"""GTM correction with one compartment per kidney at PSF 8, 9.5 and 11 mm (sensitivity analysis). Usage: python m07_gtm_one_compartment.py k_abs
Outputs: results/spect/gtm_k_abs.json
"""
from resdos.paths import ROOT
import sys, json, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.ndimage import binary_erosion, binary_dilation
import resdos.fit_tac as ft
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy'); frac=np.load(f'{ROOT}/results/spect/segfrac_spect.npy')
k=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))[sys.argv[1]]
R=json.load(open(f'{ROOT}/results/spect/organ_results_{sys.argv[1]}.json'))['organs']
V={'Liver':(1,'f3'),'Tumour':(2,'f4b'),'spleen':(3,'f4b'),'Right Kidney':(4,'f2'),'Left Kidney':(5,'f2')}
body=to_spect(truth('24h'),order=0)>0          # body outline (in patients: from CT)
bgfrac=np.clip(body.astype(float)-frac[1:].sum(0),0,1)
templates=[frac[1],frac[2],frac[3],frac[4],frac[5],bgfrac]
vois=[lab==1,lab==2,lab==3,lab==4,lab==5,binary_erosion(body,iterations=2)&(lab==0)]
def gtm_matrix(fw):
    W=np.zeros((6,6))
    for j,tpl in enumerate(templates):
        b=blur(tpl,fw)
        for i,vm in enumerate(vois): W[i,j]=b[vm].mean()
    return W
def tiac(m,y): p=ft.fit(m,T,np.asarray(y,float))[0]; return ft.integral(m,p)
out={}
for fw in [8.0,9.5,11.0]:
    W=gtm_matrix(fw)
    Cg=[]
    for i in range(6):
        m=np.array([S[i][vm].mean() for vm in vois])/k
        Cg.append(np.linalg.solve(W,m))
    Cg=np.array(Cg)
    for j,(v,(l,mod)) in enumerate(V.items()):
        It=tiac(mod,R[v]['tac_true'])
        q=tiac(mod,Cg[:,j])/It
        out.setdefault(v,{})[fw]=q
for v in V:
    q=R[v]; qs=q['Q_spect']; qm=q['Q_model']
    rc_corr=qs/qm['9.5']                          # model-based RC correction (to true concentration)
    print(f"{v:12s} uncorrected Q {qs:.3f} | RC-correction (Gaussian 9.5) {rc_corr:.3f} | GTM 8/9.5/11 mm {out[v][8.0]:.3f}/{out[v][9.5]:.3f}/{out[v][11.0]:.3f}")
json.dump({v:{str(k2):x for k2,x in d.items()} for v,d in out.items()},open(f'{ROOT}/results/spect/gtm_{sys.argv[1]}.json','w'),indent=1)
