"""3D discrete Gaussian blurring (modified-Bessel kernel) of the GT maps for FWHM 0-15 mm; DICOM output that is
imported into OpenDose3D (dose rates: see data/README.md, step "OpenDose3D").
Inputs : input/truth_{t}_reg.dcm
Outputs: results/blurred_dicom/FWHM_xx.xmm/truth_{t}_reg_G{fw}mm.dcm
Manuscript: Methods 2.2.
"""
from resdos.paths import ROOT
# -*- coding: utf-8 -*-
"""
真値画像 (truth_*_reg.dcm) に 3D Gaussian (FWHM 5, 7.5, 10, 12.5, 15 mm) をかけて DICOM 保存
  - ぼけ: 離散ガウス核（変形ベッセル関数, ITK DiscreteGaussian と同型）, 分散 σ² を厳密保持, ゼロ境界
  - 保存: 元ファイルのヘッダ・幾何情報をそのまま継承し、画素は int16 + RescaleSlope（最大値 -> 32000, 元と同じ規約）
  - UID : SOPInstanceUID / SeriesInstanceUID を新規発行（StudyInstanceUID, FrameOfReferenceUID は元のまま）
"""
import os
import numpy as np, pydicom
from pydicom.uid import generate_uid
from scipy.special import ive
from scipy.ndimage import convolve1d

IN_DIR = f'{ROOT}/input'
OUT_DIR = f'{ROOT}/results/blurred_dicom'
TP = ['1h', '4h', '24h', '40h', '72h', '144h']
FWHMS = [5.0, 7.5, 10.0, 12.5, 15.0]
MAX_STORED = 32000


def gauss1d(sigma_vox):
    t = sigma_vox ** 2
    R = int(np.ceil(5 * sigma_vox)) + 3
    n = np.arange(-R, R + 1)
    w = ive(np.abs(n), t)            # exp(-t) I_n(t)
    return w / w.sum()


def blur(img, fwhm_mm, spacing_mm):
    out = img.astype(np.float64)
    for ax, sp in enumerate(spacing_mm):
        s = fwhm_mm / (2 * np.sqrt(2 * np.log(2))) / sp
        out = convolve1d(out, gauss1d(s), axis=ax, mode='constant', cval=0.0)
    return out


def main():
    for fw in FWHMS:
        sub = os.path.join(OUT_DIR, f'FWHM_{fw:04.1f}mm')
        os.makedirs(sub, exist_ok=True)
        tag = f'G{fw:g}mm'
        for t in TP:
            ds = pydicom.dcmread(os.path.join(IN_DIR, f'truth_{t}_reg.dcm'))
            slope = float(ds.RescaleSlope); icpt = float(getattr(ds, 'RescaleIntercept', 0.0))
            phys = ds.pixel_array.astype(np.float64) * slope + icpt          # (frames, rows, cols)
            # 軸順 (z, y, x) の画素間隔。z はスライス間隔の絶対値
            sz = abs(float(getattr(ds, 'SpacingBetweenSlices', ds.SliceThickness)))
            sy, sx = map(float, ds.PixelSpacing)
            b = blur(phys, fw, (sz, sy, sx))
            b[b < 0] = 0.0
            new_slope = b.max() / MAX_STORED
            stored = np.round(b / new_slope).astype(np.int16)

            ds.PixelData = stored.tobytes()
            ds.RescaleSlope = f'{new_slope:.10g}'
            ds.RescaleIntercept = '0'
            ds.SmallestImagePixelValue = int(stored.min())
            ds.LargestImagePixelValue = int(stored.max())
            ds.WindowCenter = f'{stored.max() / 2:.1f}'
            ds.WindowWidth = f'{float(stored.max()):.1f}'
            ds['SmallestImagePixelValue'].VR = 'SS'; ds['LargestImagePixelValue'].VR = 'SS'

            ds.SOPInstanceUID = generate_uid()
            ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
            ds.SeriesInstanceUID = generate_uid()
            ds.SeriesDescription = f'truth_{t}_reg_{tag}'
            ds.SeriesNumber = 950 + int(round(fw * 10))           # 1000=5mm, 1025=7.5mm, ...
            ds.ImageType = ['DERIVED', 'SECONDARY', 'RECON TOMO', 'EMISSION']
            ds.ImageComments = (str(getattr(ds, 'ImageComments', '')) +
                                f'$3D Gaussian FWHM {fw:g} mm (discrete Bessel kernel, var=sigma^2)')[:10240]
            ds.DerivationDescription = f'3D Gaussian blur FWHM {fw:g} mm applied to truth_{t}_reg'

            ds.save_as(os.path.join(sub, f'truth_{t}_reg_{tag}.dcm'), enforce_file_format=True)
        print('saved', sub)


if __name__ == '__main__':
    main()
