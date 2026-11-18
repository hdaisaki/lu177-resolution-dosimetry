# lu177-resolution-dosimetry

[![DOI](https://zenodo.org/badge/1463427312.svg)](https://doi.org/10.5281/zenodo.24355540)

Analysis code for

> Ogino T, Ito S, Onuma K, Sasaki F, Kita T, Sagara H, Daisaki H. **Isolating the effect of spatial resolution on
> absorbed-dose estimation in ¹⁷⁷Lu-DOTATATE dosimetry: a noise-free digital phantom study with comparison to measured
> SPECT/CT.** EJNMMI Physics (submitted). DOI: *to be added*

The code reproduces every table and figure of the article and its supplementary material from

1. the public **MRTDosimetry** dataset (measured SPECT/CT of the four-organ phantom, NEMA IEC phantom; OSF,
   [doi:10.17605/OSF.IO/69NGE](https://doi.org/10.17605/OSF.IO/69NGE)), and
2. the authors' **derived data** (ground-truth activity maps, 3D Slicer segmentation, OsiriX ROIs and OpenDose3D
   dose-rate tables; Zenodo, [doi:10.5281/zenodo.24343110](https://doi.org/10.5281/zenodo.24343110)).

Starting from an empty `results/` folder, the complete run reproduces all tables of the article and supplement exactly
and eight of nine figures pixel-identically (see [docs/VERIFICATION.md](docs/VERIFICATION.md)).

**If you use this software in your research, please cite the associated publication** (see [Citation](#citation)).

## Contents

```
resdos/                 Python package: data-root handling, GT/SPECT geometry, time integration, sphere model
scripts/00_prepare/     derived inputs: VOI arrays, occupancy maps, VSV kernel, OsiriX ROIs, Gaussian-blurred GT maps
scripts/01_simulation/  digital-phantom analysis (time integration, RC/RE, same-mask test, thresholds, error components)
scripts/02_sphere_model/sphere model (Table 4, Fig. 5)
scripts/03_measured_spect/  measured SPECT: calibration, effective resolution, local model, GTM, three-stage evaluation,
                        calibration sensitivity
scripts/04_nema/        NEMA IEC recovery coefficients (volume threshold, CT-positioned spheres), count ratio
scripts/05_figures_tables/  figures and all tables (with a check against the published values)
results/published/      values published in the article (tables and analysis outputs) for comparison
tests/                  reproduction tests
docs/                   manuscript map, verification record
data/README.md          how to obtain and arrange the input data
run_all.py              runs the complete analysis in order
```

## Installation

Python 3.12 was used. Create an environment, install the dependencies and the package `resdos`:

```bash
# Linux / macOS
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e .
```

```bat
:: Windows (Anaconda Prompt); conda guarantees Python 3.12 even if another Python is installed
conda create -n resdos python=3.12 -y
conda activate resdos
pip install -r requirements.txt
pip install -e .
```

`pip install -e .` makes `resdos` importable, so that single scripts and tools can be run directly; `run_all.py`
works without it. SimpleITK and pytypedstream are needed only for the preparation steps p01 and p05. The ¹⁷⁷Lu
soft-tissue voxel S-value kernel (from PyTheranostics 0.1, Apache 2.0) is part of the derived-data record; see
[THIRD_PARTY.md](THIRD_PARTY.md).

## Data

Arrange the input data below one data root as described in [data/README.md](data/README.md):

```
<root>/input/      public MRTDosimetry images + ground-truth maps, segmentation and OsiriX ROIs (derived-data record)
<root>/derived/    OpenDose3D dose-rate tables and the voxel S-value kernel (derived-data record)
<root>/results/    created by the scripts
```

Point the code to the root (alternatively copy `config.example.json` to `config.json` and edit `"root"`):

| Shell | Command |
|---|---|
| Linux / macOS | `export RESDOS_ROOT=/path/to/root` |
| Windows cmd | `set RESDOS_ROOT=C:\path\to\root` (lost when the window is closed) |
| Windows PowerShell | `$env:RESDOS_ROOT = "C:\path\to\root"` |
| conda (any OS, persistent) | `conda env config vars set RESDOS_ROOT=C:\path\to\root PYTHONUTF8=1`, then `conda activate resdos` again |

On Windows, `PYTHONUTF8=1` is recommended. Check the derived data with `python scripts/verify_checksums.py <folder>`
(works on all systems; replaces `sha256sum -c`).

## Running

```bash
python run_all.py                 # complete analysis: about 25-40 min (Linux), about 50 min (Windows laptop)
python run_all.py --only m15      # one step;  --from m01 to start at a step;  --skip-prepare to skip p01-p05
python tests/run_tests.py         # reproduction tests (before the data are in place: 2 PASS, 3 SKIP)
python scripts/compare_with_published.py           # all outputs vs results/published, cross-platform tolerance
python scripts/compare_with_published.py --strict  # same-platform check (rtol 1e-6)
```

**What to expect.** All tables of the article are reproduced exactly on every platform (`t01_make_tables.py --check`).
Results of non-linear fits can differ between operating systems by up to a few 1e-4 (relative), because the optimiser
converges to slightly different points; the default comparison therefore uses a relative tolerance of 1e-3. The SD
and CV of non-identifiable fits (parameter CV > 1000%; only fits not used in the article) are reported as NaN.
Figures use Liberation Sans, or the metric-compatible Arial where it is not installed (e.g. Windows).

The absorbed-dose rates of the digital phantom and of the measured SPECT were calculated in **OpenDose3D**
(3D Slicer, convolution in water, no density correction) from the blurred GT maps written by `p06_blur_gt.py`; this
GUI step is not scripted, and its outputs are provided in the derived-data record (`derived/opendose3d/`).

## Manuscript map

| Item | Scripts (in order) | Main outputs |
|---|---|---|
| Table 1 | GT concentrations of the dataset (constants in the scripts) | – |
| Table 2 | s01, s02, t01 | `summary_recommended_models.csv` |
| Table 3 | s02, s03, s04, t01 | `tia_rc.json`, `same_mask_rc.json` |
| Table 4 | s05 (phantom), ss1 (spheres), t01 | `fine_fwhm.json`, `sphere_model/fwhm_limits.csv` |
| Table 5 | m01–m15, t01 | `spect/organ_results_final.json`, `spect/local_decomp_final.json`, `spect/loo_local.json`, `pvc_summary.json` |
| Fig. 1, 2, S1 | p01, p06, s01, f01 | `figs/Fig1.png`, `Fig2.png`, `FigS1.png` |
| Fig. 3 | s02, s03, ss1, f02 | `figs/Fig3.png` |
| Fig. 4 | s02, f03 | `figs/Fig4.png` |
| Fig. 5 | ss1 | `figs/Fig5.png` |
| Fig. 6 | m01–m15, f03 | `figs/Fig6.png` |
| Table S1 | s01, s02, s08, s09 | `holdout_144.json`, `gt_auc_validation.json` |
| Table S2 / S3 / S4 | s02 / s06 / s07 | `summary_recommended_models.csv`, `regression_desc.json`, `model_sensitivity.json` |
| Table S5, Fig. S2 | s09–s13, f04 | `components_0mm.json` |
| Table S6 | m13 | `spect/loo_local.json` |
| Table S7 | m04, m05, m14 | `spect/fwhm_fit.json`, `spect/unified_9.5.json` |
| Table S8, S11 | m14, m15 | `pvc_summary.json`, `unified_eval.json` |
| Table S9, Fig. S3 | n01–n04, m11, f05 | `nema/rc_compare.json`, `nema/nema_rc.json`, `nema/gauss_nema.json` |
| Table S10 | s14 | `liver_spill_full.json` |
| Table S12 | m17, m15, s09 | `tia_image_decomp.json` |
| Table S13 | m18 | `calib_sensitivity.json` |
| Text values | m19 (OpenDose3D SPECT doses, VSV agreement), n05 (count ratio), m16 (GT-image dose evaluation), s12 (segmentation agreement) | |

A detailed list with the inputs and outputs of every script is given in each script's header and in
[docs/manuscript_map.md](docs/manuscript_map.md).

## Notes

* **Ground-truth values.** Table 1 concentrations, SAAM II time-integrated activities and GT absorbed doses of the
  MRTDosimetry phantom are used as constants in the scripts. They are published with the permission of
  Dr Andrew P. Robinson and colleagues at the National Physical Laboratory (NPL), UK, the custodians of the dataset.
* **Calibration.** The measured SPECT is converted with 10.705 cps/MBq (1800 s), one quarter of the Jaszczak-based
  factor of 42.82 cps/MBq; see the article (Methods 2.6, Supplementary Methods S6, S7).
* **Differences from the analysis session** are listed in [docs/VERIFICATION.md](docs/VERIFICATION.md).
* **Use of generative AI.** The code was written with the assistance of a large language model (Claude, Anthropic).
  All code was reviewed and tested by the authors, including against independent implementations (OpenDose3D dose
  rates and the OpenDose3D recovery-coefficient tool) and by the complete reproduction documented in
  docs/VERIFICATION.md.

## Citation

If you use this software in your research, please cite the associated publication:

> Ogino T, Ito S, Onuma K, Sasaki F, Kita T, Sagara H, Daisaki H. Isolating the effect of spatial resolution on
> absorbed-dose estimation in ¹⁷⁷Lu-DOTATATE dosimetry: a noise-free digital phantom study with comparison to measured
> SPECT/CT. EJNMMI Physics (submitted). DOI: *to be added*

and, where appropriate, this software (Zenodo; version 1.0.0: [doi:10.5281/zenodo.24355541](https://doi.org/10.5281/zenodo.24355541);
all versions: [doi:10.5281/zenodo.24355540](https://doi.org/10.5281/zenodo.24355540); see [CITATION.cff](CITATION.cff)). If you
use the MRTDosimetry images, please also cite Robinson AP, Bardiès M, Denis-Bacelar AM, Lassmann M, Tipping J,
Tran-Gia J. MRTDosimetry. OSF; 2021. doi:10.17605/OSF.IO/69NGE.

## Licence

Code: MIT licence (see [LICENSE](LICENSE)); third-party components: [THIRD_PARTY.md](THIRD_PARTY.md). The
MRTDosimetry dataset is licensed under CC BY 4.0.

Contact (corresponding author): Hiromitsu Daisaki, Department of Radiological Sciences, Gunma Prefectural College of Health Sciences
(daisaki@gchs.ac.jp).
