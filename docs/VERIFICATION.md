# Verification record (version 1.0.0)

The complete analysis was run with this repository from an empty `results/` folder (Python 3.12.3, numpy 2.4.4,
scipy 1.17.1, pandas 3.0.2, matplotlib 3.10.8, pydicom 3.0.2) and compared with the outputs of the original
analysis that underlie the article.

| Check | Result |
|---|---|
| Preparation (p01-p05): VOI arrays, occupancy maps, GT stack, VSV kernel, OsiriX masks (7 files) | bit-identical |
| Blurred GT maps (p06; 30 DICOM files) | identical pixel data |
| OpenDose3D time integration (s01; 7 files, 105 rows, all three functions) | identical |
| Simulation, sphere model, measured SPECT and NEMA outputs (49 JSON files, summary table, sphere-model CSVs, intermediate arrays) | identical (max relative difference < 1e-6) |
| All tables of the article and supplement (689 values; `t01_make_tables.py --check`) | identical to the published values |
| Figures 1, 2, 4, 5, 6, S1, S2, S3 | pixel-identical |
| Figure 3 | 99.92% of pixels identical (see below) |
| `tests/run_tests.py` (5 tests) | all passed |
| Run time | about 25-40 min on Linux (time integration about 6 min); about 50 min on a Windows laptop |

## Differences from the analysis session

* **Interactive steps converted to scripts.** About twenty results were computed interactively during the analysis
  (for example the SAAM II validation, the error components, the leave-one-out test, the three-stage evaluation and
  the liver decomposition). Their code was recovered from the session records, placed in scripts with unchanged logic
  and verified by the reproduction above.
* **Figure 3, sphere-model curve.** The dotted 15.9-mL sphere curve is now computed with the stand-alone sphere model
  (`resdos/sphere_model.py`, 32 sub-voxel positions) that also produces Table 4 and Fig. 5, instead of an earlier
  prototype (6 positions). The values agree to three decimals (RC at 15 mm 0.6225 in both; maximum difference from
  the phantom tumour 0.008, as reported).
* **`spect/spect_od3d.json`.** Relative errors of the OpenDose3D doses of the measured SPECT are computed against the
  full-precision GT doses (previously rounded values); the doses are unchanged and the relative errors change by at
  most 0.1 percentage point. Fields that are no longer used were removed.
* **`nema/rc_compare.json`.** The obsolete correction by direct division of the absorbed dose by R(V) was removed;
  the RC-corrected results of the article are the proportional estimates in `pvc_summary.json`.
* **Paths.** All absolute paths of the analysis environment were replaced by paths below `RESDOS_ROOT`.
* **Not included.** Exploratory and superseded scripts of the analysis session, and the scripts that assembled the
  manuscript documents.

## Independent verification (pre-release, Windows)

A student who had not been involved in the analysis followed README and data/README alone on Windows 11 (Python
3.12.16, conda, the pinned package versions): the complete run finished without error in about 49 min; all values
of the abstract, the main text, Tables 2-5 and Supplementary Table S13 (about 250 values), all 689 table values
(`t01_make_tables.py --check`) and the content of all figures agreed with the article. Repeating the integration
step on the same computer gave byte-identical files. Differences to the Linux reference: (i) the SD/CV of the
non-identifiable liver f4b fit (parameter CV up to 8e7%) differed by up to 4.7% and its AD by 0.064%; (ii) other
fit-based outputs differed by 1e-6 to 5e-4 (relative); (iii) figures differed only in the font (Liberation Sans not
installed). The pre-release testing led to the following measures, all included in version 1.0.0:

* Non-identifiable fits (parameter CV > 1000%) are flagged (`identifiable` column) and their SD/CV of the time
  integral set to NaN (liver f4b at all FWHM; left kidney f3 at 15 mm, whose SD was already undefined).
* `compare_with_published.py` uses a cross-platform tolerance (rtol 1e-3; values of the non-identifiable liver f4b
  fit excluded, see `results/published/COMPARISON.json`), lists every value outside the tolerance, offers
  `--strict` (rtol 1e-6) and runs without PYTHONPATH; the tests report the failing values.
* Figures fall back to the metric-compatible Arial when Liberation Sans is missing, without font warnings. Fig. S2:
  legend below the axes (no overlap with the bar labels).
* `gt_closure.json` names the VOI term of the article `voi_term`.
* The voxel S-value kernel is the soft-tissue kernel of PyTheranostics 0.1 (Apache License 2.0); it is part of the
  derived data, because version 0.1 was removed from the PyPI index and later versions do not contain it (the original
  wheel remains a checksum-verified fallback). It reproduces the values logged in the original analysis (electron
  escape fraction 0.0624 -> 0.0677, central voxel 0.2066 -> 0.2624). A pre-release wrongly described it as part of
  OpenDose3D.
* The OsiriX contours are also provided as binary NRRD files (`scripts/export_osirix_masks_nrrd.py`), each checked to
  be identical to the mask used in the analysis.
* `pyproject.toml` (`pip install -e .`), `scripts/verify_checksums.py`, and instructions for Windows and for the OSF
  folders were added.

With these measures, the complete run was repeated from the distribution files on Linux: all outputs were identical to
those of the pre-release (strict comparison, rtol 1e-6) apart from the intended changes above, all tables matched
the published values, and all five tests passed. `results/published` was generated with version 1.0.0.
