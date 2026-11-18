"""Refine the sphere centres on the CT by matching the sphere walls (shell template).
Outputs: results/nema/ct_match.json
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results/nema')
import numpy as np, json
from scipy import ndimage as ndi
ct=np.load('ct.npy'); z=np.load('ct_z.npy'); ox=oy=-250.0; ps=0.976562; dz=2.5
cents=np.load('cents.npy')   # SPECT (i,r,c)
VX=4.419644; sx0,sy0,sz0=-279.947394,-280.047394,280.647394
D=np.array([37,28,22,17,13,10.])  # blobs ordered by size
spect_mm=np.array([[sx0+c[2]*VX, sy0+c[1]*VX, sz0-c[0]*VX] for c in cents])
# CT shell template matching around each SPECT centroid (smooth CT first)
cts=ndi.gaussian_filter(ct,(0.5,1.0,1.0))
zi=lambda zz: (zz-z[0])/dz
res=[]
for k,(d,(x,y,zc)) in enumerate(zip(D,spect_mm)):
    best=None
    R=d/2
    for dzz in np.arange(-4,4.01,0.5):
        for dy in np.arange(-4,4.01,0.5):
            for dx in np.arange(-4,4.01,0.5):
                cx,cy,cz=x+dx,y+dy,zc+dzz
                # sample shell (r = R+0.5..R+1.5 mm wall) minus inside/outside water
                th=np.linspace(0,np.pi,9)[1:-1]; ph=np.linspace(0,2*np.pi,16,endpoint=False)
                T,P=np.meshgrid(th,ph); u=np.stack([np.sin(T)*np.cos(P),np.sin(T)*np.sin(P),np.cos(T)],-1).reshape(-1,3)
                def samp(r):
                    p=np.array([cx,cy,cz])+r*u
                    return ndi.map_coordinates(cts,[zi(p[:,2]),(p[:,1]-oy)/ps,(p[:,0]-ox)/ps],order=1).mean()
                score=samp(R+0.8)-0.5*(samp(R-2.0)+samp(R+3.0))
                if best is None or score>best[0]: best=(score,dx,dy,dzz)
    res.append(dict(d=float(d),spect=[float(v) for v in (x,y,zc)],shift=[float(v) for v in best[1:]],score=float(best[0])))
    print(f'{d:4.0f} mm sphere: SPECT centroid ({x:6.1f},{y:6.1f},{zc:5.1f}) | CT shell offset dx {best[1]:+.1f} dy {best[2]:+.1f} dz {best[3]:+.1f} mm | wall contrast {best[0]:.0f} HU')
json.dump(res,open('ct_match.json','w'),indent=1)
