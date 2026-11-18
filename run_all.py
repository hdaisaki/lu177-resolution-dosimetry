"""Run the complete analysis in order.

    python run_all.py                 # complete analysis, including the preparation of derived inputs (00_prepare)
    python run_all.py --skip-prepare  # skip 00_prepare p01-p05 if derived/ already contains the arrays
    python run_all.py --from m01      # start at a given step
    python run_all.py --only t01      # run one step

The OpenDose3D dose-rate calculation (3D Slicer GUI) is not scripted; its outputs are part of the derived data
(derived/opendose3d). Set RESDOS_ROOT or config.json to the data root first (see README.md).
"""
import argparse, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
S = lambda p: os.path.join(HERE, 'scripts', p)
PREPARE = [('p01', '00_prepare/p01_segmentation_arrays.py', []), ('p02', '00_prepare/p02_truth_stack.py', []),
           ('p03', '00_prepare/p03_occupancy_fractions.py', []), ('p04', '00_prepare/p04_vsv_kernel.py', []),
           ('p05', '00_prepare/p05_osirix_rois.py', ['0.5'])]
STEPS = [('p06', '00_prepare/p06_blur_gt.py', []),
         ('s01', '01_simulation/s01_integrate_dose_rates.py', []), ('s02', '01_simulation/s02_summary.py', []),
         ('s03', '01_simulation/s03_tia_rc.py', []), ('s04', '01_simulation/s04_same_mask.py', []),
         ('s05', '01_simulation/s05_fine_fwhm.py', []), ('s06', '01_simulation/s06_regression.py', []),
         ('s07', '01_simulation/s07_model_sensitivity.py', []), ('s08', '01_simulation/s08_holdout.py', []),
         ('s09', '01_simulation/s09_gt_auc_validation.py', []), ('s10', '01_simulation/s10_gt_closure.py', []),
         ('s11', '01_simulation/s11_delineation.py', []), ('s12', '01_simulation/s12_osirix_vs_slicer.py', []),
         ('s13', '01_simulation/s13_components_0mm.py', []), ('s14', '01_simulation/s14_liver_spill.py', []),
         ('ss1', '02_sphere_model/run_sphere_model.py', []),
         ('m01', '03_measured_spect/m01_load_spect.py', []), ('m02', '03_measured_spect/m02_vois.py', []),
         ('m03', '03_measured_spect/m03_fwhm_fit.py', []), ('m04', '03_measured_spect/m04_abs_cal.py', []),
         ('m05', '03_measured_spect/m05_abs_checks.py', []), ('m06', '03_measured_spect/m06_organ_analysis.py', ['k_abs']),
         ('m07', '03_measured_spect/m07_gtm_one_compartment.py', ['k_abs']), ('m08', '03_measured_spect/m08_gtm_one_compartment_12mm.py', ['12']),
         *[(f'm09.{i}', '03_measured_spect/m09_local_fit.py', [str(i)]) for i in range(6)],
         ('m10', '03_measured_spect/m10_local_decomposition.py', []), ('m11', '03_measured_spect/m11_rescale_table1.py', []),
         ('m12', '03_measured_spect/m12_model_grid.py', ['0', '1', '2', '3', '4', '5']), ('m13', '03_measured_spect/m13_leave_one_out.py', []),
         *[(f'm14.{fw}', '03_measured_spect/m14_unified_eval.py', [fw]) for fw in ('8.0', '9.5', '12.0')],
         ('m15', '03_measured_spect/m15_pvc_summary.py', []), ('m16', '03_measured_spect/m16_truth_dose_eval.py', []),
         ('m17', '03_measured_spect/m17_tia_image_decomposition.py', []), ('m19', '03_measured_spect/m19_opendose3d_spect_doses.py', []),
         ('n01', '04_nema/n01_load_nema.py', []), ('n02', '04_nema/n02_ct_match.py', []), ('n03', '04_nema/n03_rc_ct_positioned.py', []),
         ('n04', '04_nema/n04_rc_volume_threshold.py', []), ('n05', '04_nema/n05_image_projection_ratio.py', []),
         ('m18', '03_measured_spect/m18_calibration_sensitivity.py', []),
         ('f01', '05_figures_tables/f01_fig1_fig2_figS1.py', []), ('f02', '05_figures_tables/f02_fig3.py', []),
         ('f03', '05_figures_tables/f03_fig4_fig6.py', []), ('f04', '05_figures_tables/f04_figS2.py', []),
         ('f05', '05_figures_tables/f05_figS3.py', []), ('t01', '05_figures_tables/t01_make_tables.py', ['--check'])]
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--skip-prepare', action='store_true'); ap.add_argument('--from', dest='start'); ap.add_argument('--only')
    a = ap.parse_args()
    steps = ([] if a.skip_prepare else PREPARE) + STEPS
    if a.only: steps = [s for s in steps if s[0] == a.only or s[0].startswith(a.only + '.')]
    if a.start: steps = steps[[i for i, s in enumerate(steps) if s[0].startswith(a.start)][0]:]
    env = dict(os.environ, PYTHONPATH=HERE + os.pathsep + os.environ.get('PYTHONPATH', ''))
    for key, script, args in steps:
        t0 = time.time(); print(f'== {key}: {script} {" ".join(args)}', flush=True)
        subprocess.run([sys.executable, S(script), *args], check=True, env=env, cwd=HERE)
        print(f'   done in {time.time()-t0:.0f} s', flush=True)
if __name__ == '__main__':
    main()
