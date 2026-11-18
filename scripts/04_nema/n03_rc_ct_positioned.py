"""Recovery coefficients with CT-positioned spheres of known diameter (>= 50% occupancy; sensitivity analysis)
and calibration variants.
Outputs: results/nema/nema_rc.json
Manuscript: Supplementary Methods S6, Table S9, Fig. S3.
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results/nema')
import numpy as np, json, sys
from scipy.optimize import curve_fit
sys.path.insert(0,f'{ROOT}/derived/sphere_pkg'); import resdos.sphere_model as sm
a=np.load('spect.npy'); cm=json.load(open('ct_match.json'))
VX=4.419644; VV=(VX/10)**3; sx0,sy0,sz0=-279.947394,-280.047394,280.647394
T=3600.0; ICF=42.82
lam=np.log(2)/(6.647*24); dt=(16+49/60+16/3600)-(10+50/60); Ctrue=2.051*np.exp(-lam*dt)   # MBq/mL
D=np.array([c['d'] for c in cm]); Vnom=np.pi/6*D**3/1000
total_true=Ctrue*Vnom.sum()
cal={'ICF 42.82 cps/MBq':1.0/(T*ICF), 'total-activity':total_true/a.sum()}     # MBq per count
off=np.mean([c['shift'] for c in cm],0)
def centre(c,mode):
    x,y,z=c['spect']
    if mode=='SPECT centroid': d=np.zeros(3)
    elif mode=='CT per sphere': d=np.array(c['shift'])
    else: d=off
    x,y,z=x+d[0],y+d[1],z+d[2]
    return np.array([(sz0-z)/VX,(y-sy0)/VX,(x-sx0)/VX])   # (i,r,c) voxel coordinates (voxel centres at integers)
def occupancy(cvox,d,sub=8):
    r=d/2/VX; lo=np.floor(cvox-r-1).astype(int); hi=np.ceil(cvox+r+2).astype(int)
    s=(np.arange(sub)+0.5)/sub-0.5
    occ=np.zeros_like(a)
    I=np.arange(lo[0],hi[0]); R=np.arange(lo[1],hi[1]); Cc=np.arange(lo[2],hi[2])
    gi=(I[:,None]+s).ravel(); gr=(R[:,None]+s).ravel(); gc=(Cc[:,None]+s).ravel()
    X,Y,Z=np.meshgrid(gi,gr,gc,indexing='ij')
    ins=((X-cvox[0])**2+(Y-cvox[1])**2+(Z-cvox[2])**2)<=r*r
    occ[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]]=ins.reshape(len(I),sub,len(R),sub,len(Cc),sub).mean(axis=(1,3,5))
    return occ
def vol_threshold_voi(cvox,V):
    # OpenDose3D-style: within a local box, threshold so that the connected region containing the centre has the nominal volume
    from scipy import ndimage as ndi
    n=max(int(round(V/VV)),1); ci=np.round(cvox).astype(int); h=8
    box=a[ci[0]-h:ci[0]+h+1,ci[1]-h:ci[1]+h+1,ci[2]-h:ci[2]+h+1]
    vals=np.sort(box.ravel())[::-1]
    for k in range(n,len(vals)):
        thr=vals[k]; lab,_=ndi.label(box>=thr); l=lab[h,h,h] if lab[h,h,h]>0 else lab.flat[np.argmax(box)]
        if (lab==l).sum()>=n: m=(lab==l); break
    out=np.zeros_like(a,bool); out[ci[0]-h:ci[0]+h+1,ci[1]-h:ci[1]+h+1,ci[2]-h:ci[2]+h+1]=m; return out
rc_fun=lambda V,al,be: 1.0/(1.0+(al/V)**be)
results={}
for mode in ['SPECT centroid','CT per sphere','SPECT + mean CT offset']:
    for cname,f in cal.items():
        rows=[]
        for c,d,V in zip(cm,D,Vnom):
            cv=centre(c,mode); occ=occupancy(cv,d); voi=occ>=0.5
            conc_rc=a[voi].mean()*f/VV/Ctrue
            act_rc=a[voi].sum()*f/(Ctrue*V)
            rows.append(dict(d=float(d),V=float(V),nvox=int(voi.sum()),Vvoi=float(voi.sum()*VV),rc_conc=float(conc_rc),rc_act=float(act_rc)))
        results[f'{mode} | {cname}']=rows
# volume-threshold VOI (SPECT-based; centre independent of CT)
for cname,f in cal.items():
    rows=[]
    for c,d,V in zip(cm,D,Vnom):
        voi=vol_threshold_voi(centre(c,'SPECT centroid'),V)
        rows.append(dict(d=float(d),V=float(V),nvox=int(voi.sum()),Vvoi=float(voi.sum()*VV),rc_conc=float(a[voi].mean()*f/VV/Ctrue),rc_act=float(a[voi].sum()*f/(Ctrue*V))))
    results[f'volume threshold (SPECT) | {cname}']=rows
fits={}
print(f'concentration at scan {Ctrue:.4f} MBq/mL; total activity {total_true:.2f} MBq; ICF-based total {a.sum()/(T*ICF):.2f} MBq ({(a.sum()/(T*ICF)/total_true-1)*100:+.1f}%)')
print('mean CT offset (mm):',np.round(off,2))
for k,rows in results.items():
    V=np.array([r['V'] for r in rows])
    line=' '.join(f"{r['d']:.0f}:{r['rc_conc']:.3f}/{r['rc_act']:.3f}" for r in sorted(rows,key=lambda r:r['d']))
    out={}
    for which in ['rc_conc','rc_act']:
        y=np.array([r[which] for r in rows])
        p,cov=curve_fit(rc_fun,V,y,p0=[2.5,0.4],bounds=([1e-3,0.05],[100,5]),maxfev=20000)
        out[which]=dict(alpha=float(p[0]),beta=float(p[1]),se=[float(x) for x in np.sqrt(np.diag(cov))])
    fits[k]=out
    print(f'{k:55s} RC conc/act by diameter: {line}\n{"":55s} fit conc: alpha {out["rc_conc"]["alpha"]:.2f} beta {out["rc_conc"]["beta"]:.2f} | act: alpha {out["rc_act"]["alpha"]:.2f} beta {out["rc_act"]["beta"]:.2f}')
json.dump(dict(results=results,fits=fits,Ctrue=Ctrue,offset=off.tolist()),open('nema_rc.json','w'),indent=1)
