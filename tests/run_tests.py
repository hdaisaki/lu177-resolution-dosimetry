"""Reproduction tests.

    python tests/run_tests.py            # runs every test whose inputs are available, skips the others

1. time integration of the Table 1 curves against the SAAM II values (no input files needed)
2. sphere-model FWHM limits (Table 4, spheres; no input files needed, about 30 s)
3. NEMA recovery coefficients with the volume threshold (needs input/GE_TEW_RR_4X10.dcm; skipped before the data are in place)
4. all results against results/published (needs a completed run_all.py)
5. all tables against the published tables (needs a completed run_all.py)
The tests also run under pytest (pytest tests/run_tests.py).
"""
import json, os, subprocess, sys, tempfile
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE); sys.path.insert(0, REPO)
from resdos.paths import ROOT
PUB = os.path.join(REPO, 'results', 'published')
ENV = dict(os.environ, PYTHONPATH=REPO + os.pathsep + os.environ.get('PYTHONPATH', ''))

try:
    import pytest
    Skip = pytest.skip.Exception
except ImportError:
    class Skip(Exception):
        pass

def test_integration_vs_saam():
    import resdos.fit_tac as ft
    T = np.array([1, 4, 24, 40, 72, 144.])
    C = {'Tumour': [1916.4519, 1809.2385, 1864.4781, 1713.1579, 1396.5609, 825.8863], 'spleen': [374.4555, 428.0064, 500.3910, 433.7019, 314.4332, 151.7456],
         'Right Kidney': [494.9367, 466.2140, 403.4738, 326.1331, 214.4026, 82.6633], 'Left Kidney': [489.9914, 456.4592, 404.1572, 325.1385, 212.3693, 82.3926],
         'Liver': [113.4567, 162.4198, 137.9349, 111.2321, 84.2231, 46.8034]}
    AUC = {'Tumour': (5267.948, 15.901), 'spleen': (7814.396, 124.948), 'Right Kidney': (4704.988, 111.4), 'Left Kidney': (5649.918, 133.768), 'Liver': (24937.645, 1306.698)}
    pub = json.load(open(os.path.join(PUB, 'gt_auc_validation.json')))
    worst = 0.0
    for v in C:
        ref = AUC[v][0] / AUC[v][1] * 1000
        for m in ('f2', 'f3', 'f4b'):
            e = ft.integral(m, ft.fit(m, T, np.array(C[v]))[0]) / ref - 1
            worst = max(worst, abs(e - pub[v]['integration'][m]))
    assert worst < 5e-4, f'integration error differs by {worst*100:.3f} points'
    return f'max difference {worst*100:.4f} percentage points'

def test_sphere_thresholds():
    with tempfile.TemporaryDirectory() as d:
        run = lambda *a: subprocess.run([sys.executable, '-m', 'resdos.sphere_model', *a], cwd=d, env=ENV, check=True, capture_output=True)
        run('simulate', '--voxel', '4.419644', '--out', 'c.csv'); run('simulate', '--voxel', '2.209822', '--out', 'f.csv')
        run('thresholds', 'c.csv', 'f.csv', '--out', 'lim.csv')
        a = pd.read_csv(os.path.join(d, 'lim.csv')); b = pd.read_csv(os.path.join(PUB, 'sphere_model', 'fwhm_limits.csv'))
    cols = [c for c in a.columns if c.startswith('fwhm_rc')]
    diff = np.nanmax(np.abs(a[cols].values - b[cols].values))
    assert diff < 1e-9, f'FWHM limits differ by {diff} mm'
    return f'{len(a)} rows identical'

def test_nema_rc():
    path = os.path.join(ROOT, 'input', 'GE_TEW_RR_4X10.dcm')
    if not os.path.exists(path): raise Skip('input/GE_TEW_RR_4X10.dcm not found')
    sys.path.insert(0, os.path.join(REPO, 'scripts', '04_nema'))
    import importlib.util
    spec = importlib.util.spec_from_file_location('vt', os.path.join(REPO, 'scripts', '04_nema', 'n04_rc_volume_threshold.py'))
    vt = importlib.util.module_from_spec(spec); spec.loader.exec_module(vt)          # functions only (main guarded)
    a = vt.load(path); pub = {int(r['d']): r['rc'] for r in json.load(open(os.path.join(PUB, 'nema', 'rc_compare.json')))['vol_threshold']}
    worst = 0.0
    for d, c in zip([37, 28, 22, 17, 13, 10], vt.sphere_centres(a)):
        vn = np.pi / 6 * d ** 3 / 1000; voi = vt.volume_threshold_voi(a, c, vn, d)
        rc = a[voi].mean() / (vt.T * vt.ICF) / vt.VV / vt.CTRUE; worst = max(worst, abs(rc - pub[d]))
    assert worst < 1e-9, f'RC differs by {worst}'
    return 'six sphere RCs identical (0.280-0.695)'

def test_results_vs_published():
    if not os.path.exists(os.path.join(ROOT, 'results', 'pvc_summary.json')): raise Skip('no completed run in RESDOS_ROOT/results')
    r = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'compare_with_published.py')], env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        lines = [l for l in r.stdout.splitlines() if l.startswith(('DIFF', 'MISSING', '  '))]
        raise AssertionError('\n' + '\n'.join(lines[:30]) + '\n' + r.stdout.strip().splitlines()[-1])
    return r.stdout.strip().splitlines()[-1]

def test_tables():
    if not os.path.exists(os.path.join(ROOT, 'results', 'pvc_summary.json')): raise Skip('no completed run in RESDOS_ROOT/results')
    r = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', '05_figures_tables', 't01_make_tables.py'), '--check'], env=ENV, capture_output=True, text=True)
    assert r.returncode == 0, '\n' + '\n'.join(r.stdout.splitlines()[-30:])
    return r.stdout.strip().splitlines()[-1]

if __name__ == '__main__':
    fails = 0
    for name, fn in [(k, v) for k, v in list(globals().items()) if k.startswith('test_')]:
        try:
            print(f'PASS  {name}: {fn()}')
        except Skip as e:
            print(f'SKIP  {name}: {e}')
        except AssertionError as e:
            fails += 1; print(f'FAIL  {name}: {e}')
    sys.exit(1 if fails else 0)
