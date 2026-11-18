"""Time integration of the OpenDose3D VOI dose-rate tables for every FWHM (fixed NUKFIT-type functions,
proportional weighting, 250 starts, analytic integration, delta-method SD); all three functions per VOI.
Inputs : derived/opendose3d/fw00, fw05, fw075, fw10, fw125, fw15 (OpenDose3D Segment-statistics .tsv), fwspect (measured SPECT)
Outputs: results/tac_fit_results_{fw}mm.csv, results/tac_fit_results_SPECT_GE670.csv
Manuscript: Methods 2.3; Tables 2, S1, S2, S4.
"""
from resdos.paths import ROOT
import subprocess, sys
DIRS={'0':'00','5':'05','7.5':'075','10':'10','12.5':'125','15':'15'}
for fw,d in DIRS.items():
    subprocess.run([sys.executable,'-m','resdos.fit_tac','--dir',f'{ROOT}/derived/opendose3d/fw{d}','--all-models','--out',f'{ROOT}/results/tac_fit_results_{fw}mm.csv'],check=True)
subprocess.run([sys.executable,'-m','resdos.fit_tac','--dir',f'{ROOT}/derived/opendose3d/fwspect','--all-models','--out',f'{ROOT}/results/tac_fit_results_SPECT_GE670.csv'],check=True)
