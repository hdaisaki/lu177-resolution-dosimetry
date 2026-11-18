"""
Sphere model: RC versus FWHM for uniform spheres (2-128 mL) on the 4.42-mm and 2.21-mm grids, FWHM limits for
RC >= 0.95/0.90/0.80, validation against the phantom tumour and the three-panel figure.
Inputs : results/tia_rc.json (observed tumour RC), results/fine_fwhm.json (phantom targets)
Outputs: results/sphere_model/{rc_grid_4.42mm.csv, rc_grid_2.21mm.csv, fwhm_limits.csv, observed_tumour_rc.csv,
         phantom_targets.csv, tumour_validation.json}, results/figs/Fig5.png
Manuscript: Methods 2.5, Table 4 (spheres), Fig. 5. Run time about 1 min.
"""
import json, os, subprocess, sys, shutil
from resdos.paths import ROOT
out=os.path.join(ROOT,'results','sphere_model'); os.makedirs(out,exist_ok=True); os.chdir(out)
tia=json.load(open(os.path.join(ROOT,'results','tia_rc.json')))['Tumour']
with open('observed_tumour_rc.csv','w') as f:
    f.write('fwhm_mm,rc\n'); [f.write(f"{float(k)},{round(v['RC_TIA'],4)}\n") for k,v in sorted(tia.items(),key=lambda x: float(x[0]))]
thr=json.load(open(os.path.join(ROOT,'results','fine_fwhm.json')))
with open('phantom_targets.csv','w') as f:
    f.write('target,volume_ml,fwhm_rc_0_90\n')
    for k,n,v in [('Tumour','Tumour',15.9),('spleen','Spleen',124.95),('Right Kidney','Right kidney',111.4),('Left Kidney','Left kidney',133.8)]:
        f.write(f"{n},{v},{round(thr[k]['thr']['0.9'],2)}\n")
sm=[sys.executable,'-m','resdos.sphere_model']
subprocess.run(sm+['simulate','--voxel','4.419644','--out','rc_grid_4.42mm.csv'],check=True)
subprocess.run(sm+['simulate','--voxel','2.209822','--out','rc_grid_2.21mm.csv'],check=True)
subprocess.run(sm+['thresholds','rc_grid_4.42mm.csv','rc_grid_2.21mm.csv','--out','fwhm_limits.csv'],check=True)
subprocess.run(sm+['figure','--coarse','rc_grid_4.42mm.csv','--fine','rc_grid_2.21mm.csv','--observed','observed_tumour_rc.csv',
                   '--targets','phantom_targets.csv','--out','sphere_model_figure.png'],check=True)
# validation curve of the 15.9-mL sphere (tumour-to-liver ratio of the GT) for Fig. 3 and Results 3.4
import numpy as np, resdos.sphere_model as smod
fw=[0.0,2.5,5.0,7.5,10.0,12.5,15.0]
rc,sd,base,n=smod.rc_curve(15.9,4.419644,fw,16.812082254899174,32,6,np.random.default_rng(0))
json.dump(dict(V=15.9,contrast=16.812082254899174,FWHM=fw,RC=[float(x) for x in rc]),open('tumour_validation.json','w'),indent=1)
shutil.copy('sphere_model_figure.png',os.path.join(ROOT,'results','figs','Fig5.png'))
from resdos.manifest import record; record('figs/Fig5.png','02_sphere_model/run_sphere_model.py',3)
