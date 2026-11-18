"""Doses from time-integrated activity images (SAAM II values, fitted Table 1 curves, GTM) on the SPECT grid with
identical templates; evaluation in binary VOIs, occupancy-weighted and interior voxels.
Outputs: results/tia_image_decomp.json
Manuscript: Table S12, Results 3.5.
"""
from resdos.paths import ROOT
import sys, json, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.signal import fftconvolve
from scipy.ndimage import binary_erosion
import resdos.fit_tac as ft
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy'); frac=np.load(f'{ROOT}/results/spect/segfrac_spect.npy')
K=np.load(f'{ROOT}/derived/vsv_lu177_4.42mm.npy'); kabs=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json'))['k_abs']; KT=0.5388; VV=(VX/10)**3
F=np.load(f'{ROOT}/derived/frac_tp_avg.npy')
cR=to_spect(F[...,7],order=1); mR=to_spect(F[...,8],order=1); cL=to_spect(F[...,5],order=1); mL=to_spect(F[...,6],order=1)
body=to_spect(truth('24h'),order=0)>0
names=['Liver','Tumour','Spleen','R cortex','R medulla','L cortex','L medulla']
tpl=[frac[1],frac[2],frac[3],cR,mR,cL,mL]; bg=np.clip(body.astype(float)-sum(tpl),0,1); tpl=tpl+[bg]; names=names+['Background']
# --- true concentration curves (kBq/mL) and SAAM II time-integrated concentrations (kBq h/mL)
C={'Liver':[113.4567,162.4198,137.9349,111.2321,84.2231,46.8034],'Tumour':[1916.4519,1809.2385,1864.4781,1713.1579,1396.5609,825.8863],
   'Spleen':[374.4555,428.0064,500.3910,433.7019,314.4332,151.7456],'R cortex':[632.1880,594.7442,515.5828,416.1147,273.5195,105.3095],
   'R medulla':[204.2500,193.9979,166.0365,135.5598,89.1979,34.7007],'L cortex':[622.0184,580.2426,512.1216,413.9566,270.3382,104.9750],
   'L medulla':[210.3115,194.2422,175.4506,136.9903,89.5705,34.5551],'Background':[101.0,67.4,12.7,9.5,6.9,3.5]}
AUC={'Liver':(24937.645,1306.698),'Tumour':(5267.948,15.901),'Spleen':(7814.396,124.948),'R cortex':(4065.179,75.671),'R medulla':(639.808,35.729),
     'L cortex':(4881.750,90.871),'L medulla':(768.167,42.897)}
fn={'Liver':'f3','Tumour':'f4b','Spleen':'f4b','R cortex':'f2','R medulla':'f2','L cortex':'f2','L medulla':'f2'}
T_=np.array([1,4,24,40,72,144.]); lp=4.345e-3
def trap_tail(y):
    y=np.asarray(y); lam=max(-np.polyfit(T_[-3:],np.log(y[-3:]),1)[0],lp); return 0.5*y[0]+np.trapezoid(y,T_)+y[-1]/lam
I=lambda m,y: ft.integral(m,ft.fit(m,T_,np.asarray(y,float))[0])
tiac_saam={n:(AUC[n][0]/AUC[n][1]*1000 if n in AUC else trap_tail(C[n])) for n in names}
tiac_fit={n:(I(fn[n],C[n]) if n in fn else trap_tail(C[n])) for n in names}
# --- GTM compartment concentrations (PSF 9.5 mm)
rows=[lab==1,lab==2,lab==3,(cR>=0.5)&(lab==4),(mR>=0.5)&(lab==4),(cL>=0.5)&(lab==5),(mL>=0.5)&(lab==5),binary_erosion(body,iterations=2)&(lab==0)]
B=[blur(t,9.5) for t in tpl]; W=np.array([[b[r].mean() for b in B] for r in rows])
cg=np.array([np.linalg.solve(W,np.array([S[i][r].mean() for r in rows])/kabs)*KT for i in range(6)])   # (6, 8) kBq/mL
tiac_gtm={n:(I(fn[n],cg[:,j]) if n in fn else trap_tail(np.clip(cg[:,j],1e-3,None))) for j,n in enumerate(names)}
def dose(tiac):  # mGy
    img=sum(tiac[n]*t for n,t in zip(names,tpl))*VV/1000.0      # MBq h per voxel
    return fftconvolve(img,K,mode='same')*3600.0
D={k:dose(v) for k,v in [('saam',tiac_saam),('fit',tiac_fit),('gtm',tiac_gtm)]}
organs={'Liver':(1,[0]),'Tumour':(2,[1]),'Spleen':(3,[2]),'Right kidney':(4,[3,4]),'Left kidney':(5,[5,6])}
GT={'Liver':1.702361,'Tumour':28.12285,'Spleen':5.382068,'Right kidney':3.598351,'Left kidney':3.599789}
out={}
for o,(l,idx) in organs.items():
    w=sum(tpl[i] for i in idx); reg={'binary VOI':lab==l,'occupancy-weighted':w,'interior (>0.99)':w>0.99}
    out[o]={}
    for rn,r in reg.items():
        vals={}
        for k,d in D.items():
            vals[k]=float((d*r).sum()/r.sum()/1000) if rn=='occupancy-weighted' else float(d[r].mean()/1000)
        out[o][rn]=vals
json.dump(dict(out=out,tiac_saam=tiac_saam,tiac_fit=tiac_fit,tiac_gtm=tiac_gtm),open(f'{ROOT}/results/tia_image_decomp.json','w'),indent=1)
print('AD vs GT (%)  [SAAM II TIA image | fitted Table-1 TIA image | GTM TIA image]')
for o in organs:
    for rn,v in out[o].items():
        print(f'{o:12s} {rn:20s} {(v["saam"]/GT[o]-1)*100:+6.1f} | {(v["fit"]/GT[o]-1)*100:+6.1f} | {(v["gtm"]/GT[o]-1)*100:+6.1f}')
