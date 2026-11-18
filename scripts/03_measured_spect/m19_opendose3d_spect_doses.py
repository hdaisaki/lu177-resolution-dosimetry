"""OpenDose3D absorbed doses of the measured SPECT (selected functions) and the agreement of the authors' voxel
S-value dose rates with the OpenDose3D dose rates.
Inputs : results/tac_fit_results_SPECT_GE670.csv, derived/opendose3d/fwspect, results/spect/unified_9.5.json
Outputs: results/spect/spect_od3d.json
Manuscript: Results 3.5 (OpenDose3D ADs; VSV agreement).
"""
from resdos.paths import ROOT
import json, glob, numpy as np, pandas as pd
d=pd.read_csv(f'{ROOT}/results/tac_fit_results_SPECT_GE670.csv')
GT={'Tumour':28.12285,'spleen':5.382068,'Right Kidney':3.598351,'Left Kidney':3.599789,'Liver':1.702361}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
key=lambda s: float(s.split('_')[-1].replace('HR.tsv',''))
DR=pd.concat([pd.read_csv(f,sep='\t').set_index('Segment')['Mean'].rename(key(f)) for f in sorted(glob.glob(f'{ROOT}/derived/opendose3d/fwspect/*.tsv'),key=key)],axis=1)
U=json.load(open(f'{ROOT}/results/spect/unified_9.5.json')); out={}
for v in GT:
    r=d[(d.VOI==v)&(d.model==rec[v])].iloc[0]
    vsv=float(np.mean(np.array(U[v]['dr_unc'])/DR.loc[v].values)-1)*100
    out[v]=dict(AD=float(r.AD_Gy),SD=float(r.AD_SD_Gy),RE=float((r.AD_Gy/GT[v]-1)*100),own_vsv_vs_opendose3d_pct=vsv)
    print(f"{v:12s} OpenDose3D AD {r.AD_Gy:.2f} Gy ({out[v]['RE']:+.1f}%) | own VSV vs OpenDose3D dose rate {vsv:+.1f}%")
json.dump(out,open(f'{ROOT}/results/spect/spect_od3d.json','w'),indent=1)
