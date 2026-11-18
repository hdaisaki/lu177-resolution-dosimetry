"""Absolute calibration (10.705 cps/MBq for 1800 s) and total activity in the field of view.
Outputs: results/spect/fwhm_fit.json (k_abs, abs_rows)
Manuscript: Methods 2.6, Table S7.
"""
from resdos.paths import ROOT
import sys, json, importlib.util; sys.path.insert(0,f'{ROOT}/derived/spect')
from resdos.common import *
from scipy.ndimage import binary_erosion
S=np.load(f'{ROOT}/results/spect/spect.npy'); lab=np.load(f'{ROOT}/results/spect/seglab_spect.npy')
VV=(VX/10)**3
KBQML_PER_COUNT=1000.0/(10.705*1800.0)/VV           # kBq/mL per reconstructed count
sl=np.array([float(pydicom.dcmread(f'{ROOT}/input/truth_{t}_reg.dcm',stop_before_pixels=True).RescaleSlope) for t in TS])
# kBq/mL per truth unit, from Table 1 (tumour 1916.5 kBq/mL at 1 h = raw 32000 * slope)
KBQML_PER_TRUTH=1916.5/(32000*sl[0])
k_abs=KBQML_PER_TRUTH/KBQML_PER_COUNT                 # expected counts per truth unit
print('kBq/mL per count %.5f | kBq/mL per truth unit %.5f | expected counts per truth unit %.4f'%(KBQML_PER_COUNT,KBQML_PER_TRUTH,k_abs))
# check Table-1 factor consistency across organs/time
tab={'Tumour':[1916.5,1809.2,1864.5,1713.2,1396.6,825.9],'Liver':[113.5,162.4,137.9,111.2,84.2,46.8],'Spleen':[374.5,428.0,500.4,433.7,314.4,151.7],'Background':[101.0,67.4,12.7,9.5,6.9,3.5]}
raw={'Tumour':[32000]*6,'Liver':[1895,2872,2367,2077,1929,1813],'Spleen':[6253,7570,8588,8101,7204,5878],'Background':[1686,1192,218,177,158,136]}
for k_,v in tab.items(): print(' table/truth',k_,np.round(np.array(v)/(np.array(raw[k_])*sl),4))
lc=binary_erosion(lab==1,iterations=3)
rows=[]
for i,t in enumerate(TS):
    tr=truth(t); m=to_spect(blur(tr,9.5),order=3)
    spect_kbq=S[i]*KBQML_PER_COUNT; tr_kbq=m*KBQML_PER_TRUTH
    r_lc=spect_kbq[lc].mean()/tr_kbq[lc].mean()
    tot_spect=S[i].sum()/(10.705*1800.0); tot_true=tr.sum()*KBQML_PER_TRUTH*VV/1000.0
    rows.append((t,float(r_lc),float(tot_spect),float(tot_true)))
    print(f'{t:4s} liver-interior SPECT/true {r_lc:.3f} | total activity in FOV: SPECT {tot_spect:7.1f} MBq vs true {tot_true:7.1f} MBq ({(tot_spect/tot_true-1)*100:+.1f}%)')
r=np.array([x[1] for x in rows]); print('median liver-interior bias %+.1f%%'%((np.median(r)-1)*100))
d=json.load(open(f'{ROOT}/results/spect/fwhm_fit.json')); d['k_abs']=float(k_abs); d['abs_rows']=rows; d['kbqml_per_count']=KBQML_PER_COUNT; d['kbqml_per_truth']=KBQML_PER_TRUTH
json.dump(d,open(f'{ROOT}/results/spect/fwhm_fit.json','w'),indent=1)
