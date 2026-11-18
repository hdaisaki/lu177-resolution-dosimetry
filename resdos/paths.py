"""Location of the data root.

All scripts read and write below one root directory (RESDOS_ROOT) with the layout
    input/    public MRTDosimetry images and the authors' ground-truth (GT) activity maps
    derived/  author-made inputs (VOI arrays, occupancy maps, VSV kernel, OpenDose3D tables, OsiriX ROIs)
    results/  everything produced by the scripts
The root is taken from the environment variable RESDOS_ROOT, else from config.json next to this package,
else the current working directory.
"""
import json, os
from pathlib import Path

def _root():
    if os.environ.get('RESDOS_ROOT'):
        return Path(os.environ['RESDOS_ROOT']).expanduser().resolve()
    cfg = Path(__file__).resolve().parent.parent / 'config.json'
    if cfg.exists():
        return Path(json.load(open(cfg))['root']).expanduser().resolve()
    return Path.cwd()

ROOT = str(_root())
for sub in ('derived', 'results', 'results/spect', 'results/nema', 'results/figs'):
    os.makedirs(os.path.join(ROOT, sub), exist_ok=True)
