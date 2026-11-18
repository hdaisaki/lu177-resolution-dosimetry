"""Locate the 177Lu soft-tissue voxel S-value kernel used in the article.

DK_softtissueICRP_mGyperdecay.img (51 x 51 x 51 voxels, 4.7952 mm) was distributed in the wheel of PyTheranostics 0.1
(doodle/kernels/; Apache License 2.0). Version 0.1 is no longer listed on PyPI and later versions do not contain the
kernel, so the unmodified file is redistributed, with its licence, in the derived-data record of this work
(derived/kernels/). Search order: derived/kernels/ -> derived/opendose3d_kernels/ (layout of a pre-release) ->
download of the original wheel from the PyPI file server. The SHA-256 is checked in every case.
"""
import hashlib, io, os, shutil, urllib.request, zipfile

NAME = 'DK_softtissueICRP_mGyperdecay.img'
SHA256 = '533be3195abccfcdd4ebe8fe0bb8e969539c537f52df2b3e62233d89dbba1022'
WHEEL_MEMBER = 'doodle/kernels/' + NAME
WHEEL_URL = ('https://files.pythonhosted.org/packages/10/03/bc2f48835ef7a6eeb545b4c111943fac98b2ec33988b04a694826f1d3013/'
             'PyTheranostics-0.1-py3-none-any.whl')

MANUAL = f"""
The voxel S-value kernel {NAME} was not found.
It is part of the derived-data record of this work: extract the record into the data root so that
    <RESDOS_ROOT>/derived/kernels/{NAME}
exists (see data/README.md). Expected SHA-256: {SHA256}
"""

def _sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()

def ensure_kernel(derived_dir):
    """Return the path of the kernel, from the derived data or, as a fallback, from the original wheel."""
    target = os.path.join(derived_dir, 'kernels', NAME)
    legacy = os.path.join(derived_dir, 'opendose3d_kernels', NAME)
    for p in (target, legacy):
        if os.path.exists(p):
            if _sha(p) != SHA256:
                raise SystemExit(f'{p} has an unexpected SHA-256.' + MANUAL)
            return p
    try:
        print(f'kernel not in the derived data; trying the original PyTheranostics 0.1 wheel: {WHEEL_URL}')
        whl = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(WHEEL_URL, timeout=120).read()))
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with whl.open(WHEEL_MEMBER) as src, open(target, 'wb') as dst:
            shutil.copyfileobj(src, dst)
    except Exception as e:
        raise SystemExit(f'Download failed ({e}).' + MANUAL)
    if _sha(target) != SHA256:
        os.remove(target)
        raise SystemExit('Downloaded kernel has an unexpected SHA-256.' + MANUAL)
    print(f'kernel saved to {target} (SHA-256 verified)')
    return target
