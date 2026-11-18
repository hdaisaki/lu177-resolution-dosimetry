"""OpenDose3D-style volume-threshold sphere VOIs (region grown from a seed at each sphere centre to the nominal volume),
RC = VOI concentration / true concentration (42.82 cps/MBq, 3600 s), and fit of R(V) = 1/[1 + (alpha/V)^beta].
Outputs: results/nema/rc_compare.json
Manuscript: Methods 2.6, Supplementary Methods S6, Table S9.
"""
from resdos.paths import ROOT
"""
Recovery coefficients of the NEMA IEC spheres (MRTDosimetry dataset, GE Discovery 670, Xeleris OSEM 4x10 RR, TEW)
with the OpenDose3D-style volume-based threshold (seeded region growing on SPECT), curve fit and application to
the SPECT-based absorbed doses of the anthropomorphic phantom.
Inputs : GE_TEW_RR_4X10.dcm (NEMA SPECT); sphere centres (SPECT centroids, see nema_rc.py)
Settings: ICF 42.82 cps/MBq, acquisition duration 3600 s, 2.051 MBq/mL at 2019-02-22 10:50 UTC, scan 16:49:16
"""
import numpy as np, pydicom
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
VX = 4.419644; VV = (VX / 10) ** 3
T, ICF = 3600.0, 42.82
lam = np.log(2) / (6.647 * 24); dt = (16 + 49 / 60 + 16 / 3600) - (10 + 50 / 60); CTRUE = 2.051 * np.exp(-lam * dt)
def load(path=f'{ROOT}/input/GE_TEW_RR_4X10.dcm'):
    return pydicom.dcmread(path).pixel_array.astype(float)
def sphere_centres(a, n=6):
    lab, m = ndi.label(a > 0.15 * a.max()); sizes = ndi.sum(np.ones_like(a), lab, range(1, m + 1))
    order = np.argsort(sizes)[::-1][:n]
    return [np.array(ndi.center_of_mass(a * (lab == k + 1))) for k in order]   # largest to smallest
def volume_threshold_voi(a, seed, volume_ml, diameter_mm):
    """grow the region connected to the local maximum near the seed by lowering the threshold until it reaches the nominal volume"""
    n = max(int(round(volume_ml / VV)), 1); h = int(np.ceil(diameter_mm / 2 / VX)) + 4; ci = np.round(seed).astype(int)
    box = a[ci[0]-h:ci[0]+h+1, ci[1]-h:ci[1]+h+1, ci[2]-h:ci[2]+h+1]
    nb = box[h-1:h+2, h-1:h+2, h-1:h+2]; m0 = tuple(np.array(np.unravel_index(np.argmax(nb), nb.shape)) + h - 1)
    prev = None
    for thr in np.sort(np.unique(box[box <= box[m0]]))[::-1]:
        lab, _ = ndi.label(box >= thr); reg = lab == lab[m0]
        if reg.sum() >= n:
            if prev is not None and abs(prev.sum() - n) < abs(reg.sum() - n): reg = prev
            break
        prev = reg
    out = np.zeros_like(a, bool); out[ci[0]-h:ci[0]+h+1, ci[1]-h:ci[1]+h+1, ci[2]-h:ci[2]+h+1] = reg
    return out
rc_fun = lambda V, al, be: 1.0 / (1.0 + (al / V) ** be)
if __name__ == '__main__':
    a = load(); diam = [37, 28, 22, 17, 13, 10]; cents = sphere_centres(a)
    V, RC = [], []
    for d, c in zip(diam, cents):
        vn = np.pi / 6 * d ** 3 / 1000; voi = volume_threshold_voi(a, c, vn, d)
        rc = a[voi].mean() / (T * ICF) / VV / CTRUE; V.append(voi.sum() * VV); RC.append(rc)
        print(f'{d:2d} mm: {voi.sum():3d} voxels ({voi.sum()*VV:5.2f} mL)  RC {rc:.3f}')
    p, _ = curve_fit(rc_fun, np.array(V), np.array(RC), p0=[2.5, 0.4], bounds=([1e-3, 0.05], [100, 5]))
    print(f'fit: alpha {p[0]:.2f} mL, beta {p[1]:.2f}   (OpenDose3D: 2.53, 0.41)')
    # ---- archive the values used in Table S9 and the RC-curve fits
    import json
    geo=json.load(open(f'{ROOT}/results/nema/nema_rc.json'))['results']['CT per sphere | ICF 42.82 cps/MBq']
    g=sorted(geo,key=lambda r:r['d']); pg,_=curve_fit(rc_fun,np.array([r['V'] for r in g]),np.array([r['rc_conc'] for r in g]),p0=[2.5,0.4],bounds=([1e-3,0.05],[100,5]))
    od3dV=np.array([0.52,1.21,2.5,5.61,11.57,26.68]); od3dRC=np.array([0.28,0.44,0.57,0.60,0.63,0.69])   # OpenDose3D tool output (this NEMA image)
    po,_=curve_fit(rc_fun,od3dV,od3dRC,p0=[2.5,0.4],bounds=([1e-3,0.05],[100,5]))
    rows=[dict(d=float(d),V=float(np.pi/6*d**3/1000),nvox=int(round(v/VV)),Vvoi=float(v),rc=float(r)) for d,v,r in zip(diam,V,RC)]
    json.dump(dict(vol_threshold=sorted(rows,key=lambda r:r['d']),fits={'OpenDose3D tool (this NEMA)':[float(x) for x in po],'Volume threshold, this work':[float(x) for x in p],
              'CT-positioned sphere VOI, this work':[float(x) for x in pg]}),open(f'{ROOT}/results/nema/rc_compare.json','w'),indent=1)
