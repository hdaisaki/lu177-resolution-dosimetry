"""As m07 for PSF 12 mm. Usage: python m08_gtm_one_compartment_12mm.py 12
Outputs: results/spect/gtm_k_abs.json (adds 12 mm)
"""
from resdos.paths import ROOT
import sys, json, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.ndimage import binary_erosion
import resdos.fit_tac as ft
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy'); frac=np.load(f'{ROOT}/results/spect/segfrac_spect.npy')
k=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['k_abs']
R=json.load(open(f'{ROOT}/results/spect/organ_results_k_abs.json'))['organs']
V={'Liver':(1,'f3'),'Tumour':(2,'f4b'),'spleen':(3,'f4b'),'Right Kidney':(4,'f2'),'Left Kidney':(5,'f2')}
body=to_spect(truth('24h'),order=0)>0
bgfrac=np.clip(body.astype(float)-frac[1:].sum(0),0,1)
templates=[frac[1],frac[2],frac[3],frac[4],frac[5],bgfrac]
vois=[lab==1,lab==2,lab==3,lab==4,lab==5,binary_erosion(body,iterations=2)&(lab==0)]
fw=float(sys.argv[1])
W=np.zeros((6,6))
for j,tpl in enumerate(templates):
    b=blur(tpl,fw)
    for i,vm in enumerate(vois): W[i,j]=b[vm].mean()
Cg=np.array([np.linalg.solve(W,np.array([S[i][vm].mean() for vm in vois])/k) for i in range(6)])
G=json.load(open(f'{ROOT}/results/spect/gtm_k_abs.json'))
for j,(v,(l,mod)) in enumerate(V.items()):
    It=ft.integral(mod,ft.fit(mod,T,np.array(R[v]['tac_true']))[0])
    G[v][str(fw)]=ft.integral(mod,ft.fit(mod,T,Cg[:,j])[0])/It
json.dump(G,open(f'{ROOT}/results/spect/gtm_k_abs.json','w'),indent=1)
print({v:round(G[v][str(fw)],3) for v in V})
