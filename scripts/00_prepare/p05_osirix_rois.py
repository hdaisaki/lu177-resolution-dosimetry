"""Read the OsiriX MD ROI series used to build the GT maps and rasterise it at CT resolution.
Inputs : input/CTStd_2_5mmTomo_Day0.rois_series
Outputs: derived/osirix_masks_+0.5.npz (argument 0.5 = OsiriX pixel-centre convention)
Requires pytypedstream (pip install pytypedstream). Usage: python p05_osirix_rois.py 0.5
"""
from resdos.paths import ROOT
import typedstream, re, numpy as np, json
from matplotlib.path import Path
root=typedstream.unarchive_from_file(f'{ROOT}/input/CTStd_2_5mmTomo_Day0.rois_series')
slices=list(list(root.elements)[0].elements)
NZ=len(slices); N=512
masks={k:np.zeros((NZ,N,N),bool) for k in ['outline','liver','spleen','R-kidney','L-kidney','inner-R-kidney','inner-L-kidney','tumor']}
def raster(shift):
    yy,xx=np.mgrid[0:N,0:N]; P=np.c_[xx.ravel()+shift,yy.ravel()+shift]
    M={k:np.zeros((NZ,N,N),bool) for k in masks}
    for zi,sl in enumerate(slices):
        for r in sl.elements:
            c=r.contents; typ=c[2].value.value; name=str(c[10].value.value)
            if typ==11.0:
                pts=[]
                for p in c[0].value.elements:
                    m=re.findall(r'[-\d.eE]+',str(p.contents[0].value.value)); pts.append((float(m[0]),float(m[1])))
                if len(pts)<3: continue
                inside=Path(pts).contains_points(P).reshape(N,N)
            else:
                m=[float(x) for x in re.findall(r'[-\d.eE]+',str(c[1].value.value))]
                cx,cy,rx,ry=m
                inside=((xx-cx+shift)/rx)**2+((yy-cy+shift)/ry)**2<=1.0
            M[name][zi]|=inside
    return M
if __name__=='__main__':
    import sys
    shift=float(sys.argv[1]) if len(sys.argv)>1 else 0.0
    M=raster(shift)
    np.savez_compressed(f'{ROOT}/derived/osirix_masks_{shift:+.1f}.npz',**{k.replace('-','_'):v for k,v in M.items()})
    print({k:int(v.sum()) for k,v in M.items()})
