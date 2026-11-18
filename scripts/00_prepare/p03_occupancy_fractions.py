"""Estimate time-invariant compartment occupancy fractions of each GT voxel (background, lung, liver, spleen,
renal cortex/medulla, tumour) from the six GT maps and the Table 1 concentrations.
Inputs : derived/raw.npy
Outputs: derived/frac_tp_avg.npy
Manuscript: Methods 2.1; used for GTM templates, delineation and spill-in analyses.
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/derived')
import numpy as np, pydicom
from scipy import ndimage as ndi
UP=f'{ROOT}/input'; TP=['1h','4h','24h','40h','72h','144h']
raw=np.load('raw.npy')
CLASSES={ # 生画素値(6時点)
 'Air':[0]*6,
 'BG':[1686,1192,218,177,158,136],
 'Lung':[422,298,54,44,40,34],
 'Liver':[1895,2872,2367,2077,1929,1813],
 'Spleen':[6253,7570,8588,8101,7204,5878],
 'KidLC':[10386,10262,8789,7733,6193,4068],
 'KidLM':[3511,3435,3012,2559,2053,1341],
 'KidRC':[10556,10519,8849,7772,6267,4080],
 'KidRM':[3410,3431,2849,2533,2044,1344],
 'Tumor':[32000]*6,
}
cn=list(CLASSES); K=len(cn)
Fsum=np.zeros(raw.shape[1:]+(K,))
for ti in range(6):
    img=raw[ti]
    vals=np.array([CLASSES[c][ti] for c in cn])
    P=np.full(img.shape,-1,int)
    for k,v in enumerate(vals): P[img==v]=k
    # 近傍の純粋ラベル頻度 (3x3x3, 足りなければ 5x5x5)
    cnt3=np.stack([ndi.uniform_filter((P==k).astype(float),3,mode='constant')*27 for k in range(K)],-1)
    cnt5=np.stack([ndi.uniform_filter((P==k).astype(float),5,mode='constant')*125 for k in range(K)],-1)
    Fr=np.zeros(img.shape+(K,))
    pure=P>=0
    Fr[pure,P[pure]]=1
    imp=np.argwhere(~pure)
    for (z,y,x) in imp:
        v=img[z,y,x]
        c=cnt3[z,y,x]+1e-3*cnt5[z,y,x]
        cand=np.argsort(-c)[:4]; cand=[k for k in cand if c[k]>0]
        best=None
        # 候補ペアのうち、v を [0,1] の分率で説明でき、1/8刻みに最も近いもの（同程度なら近傍頻度の高いペア）
        for i in range(len(cand)):
            for j in range(i+1,len(cand)):
                a,b=cand[i],cand[j]
                if vals[a]==vals[b]: continue
                fb=(v-vals[a])/(vals[b]-vals[a])
                if -0.02<=fb<=1.02:
                    err=abs(fb*8-round(fb*8))/8 - 0.001*(c[a]+c[b])
                    if best is None or err<best[0]: best=(err,a,b,min(max(fb,0),1))
        if best is None:  # 説明できない場合は最近傍値のクラス
            k=int(np.argmin(np.abs(vals-v))); Fr[z,y,x,k]=1
        else:
            _,a,b,fb=best; Fr[z,y,x,a]=1-fb; Fr[z,y,x,b]=fb
    Fsum+=Fr
    print(TP[ti],'impure voxels',len(imp))
Fm=Fsum/6
np.save('frac_tp_avg.npy',Fm.astype(np.float32))
json_c=cn
vv=(4.419644/10)**3
for k,c in enumerate(cn):
    print(f'{c:7s} true vol {Fm[...,k].sum()*vv:8.1f} mL   f>=0.5: {(Fm[...,k]>=0.5).sum()*vv:8.1f} mL')
print('kidL',(Fm[...,5]+Fm[...,6]).sum()*vv,'kidR',(Fm[...,7]+Fm[...,8]).sum()*vv)
