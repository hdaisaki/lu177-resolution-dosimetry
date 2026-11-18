"""Export the OsiriX MD compartment contours (rasterised by p05) as binary NRRD files in an open format.

The masks of derived/osirix_masks_+0.5.npz (CT resolution, OsiriX pixel-centre convention) are written as one
binary NRRD per structure, with the geometry of the 3D Slicer segmentation of the same CT (input/MRT_Segmentation.nrrd);
slices without OsiriX contours are zero. Each file is read back and compared with the mask used in the analysis.
Usage: python scripts/export_osirix_masks_nrrd.py   (after step p05)
Output: <RESDOS_ROOT>/input/osirix_masks/osirix_<structure>.nrrd (uint8, 1 = inside, gzip-compressed)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, SimpleITK as sitk
from resdos.paths import ROOT

NAMES = {'liver': 'liver (including the tumour insert)', 'tumor': 'spherical tumour insert', 'spleen': 'spleen',
         'R_kidney': 'right kidney (outer contour)', 'L_kidney': 'left kidney (outer contour)',
         'inner_R_kidney': 'right kidney, inner contour', 'inner_L_kidney': 'left kidney, inner contour',
         'outline': 'body outline of the phantom'}

def main():
    ref = sitk.ReadImage(os.path.join(ROOT, 'input', 'MRT_Segmentation.nrrd'))
    nz, ny, nx = sitk.GetArrayViewFromImage(ref).shape
    M = np.load(os.path.join(ROOT, 'derived', 'osirix_masks_+0.5.npz'))
    out = os.path.join(ROOT, 'input', 'osirix_masks'); os.makedirs(out, exist_ok=True)
    for key in NAMES:
        m = M[key]
        assert m.shape[1:] == (ny, nx) and m.shape[0] <= nz, (key, m.shape)
        full = np.zeros((nz, ny, nx), np.uint8); full[:m.shape[0]] = m
        img = sitk.GetImageFromArray(full); img.CopyInformation(ref)
        path = os.path.join(out, f'osirix_{key}.nrrd'); sitk.WriteImage(img, path, useCompression=True)
        back = sitk.GetArrayFromImage(sitk.ReadImage(path))
        assert np.array_equal(back[:m.shape[0]].astype(bool), m) and not back[m.shape[0]:].any(), key
        print(f'{path}: {int(m.sum())} voxels, identical to the analysis mask')

if __name__ == '__main__':
    main()
