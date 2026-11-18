"""Compare the results of a run (RESDOS_ROOT/results) with the values published in the article (results/published).

Every numeric value of every JSON/CSV output is compared with a relative tolerance:
    python scripts/compare_with_published.py            # cross-platform check (rtol 1e-3, see below)
    python scripts/compare_with_published.py --strict   # same-platform check (rtol 1e-6, nothing excluded)
Non-linear fits converge to slightly different optima on different operating systems and numerical libraries
(relative differences up to a few 1e-4 were observed between Linux and Windows). Values of fits flagged as
non-identifiable (only the unused liver f4b fit; listed in results/published/COMPARISON.json) are not compared in
the cross-platform check. The tables of the article are checked separately and exactly by t01_make_tables.py --check.
The script lists every value outside the tolerance and exits with code 1 if there is any.
"""
import argparse, glob, json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, pandas as pd
from resdos.paths import ROOT
PUB = os.path.join(os.path.dirname(HERE), 'results', 'published')

def leaves(x, p=''):
    if isinstance(x, dict):
        for k, v in x.items(): yield from leaves(v, f'{p}/{k}')
    elif isinstance(x, list):
        for i, v in enumerate(x): yield from leaves(v, f'{p}[{i}]')
    elif isinstance(x, (int, float)) and not isinstance(x, bool): yield p, float(x)

def load(path, skip_nonidentifiable):
    if path.endswith('.json'):
        return dict(leaves(json.load(open(path, encoding='utf-8'))))
    d = pd.read_csv(path)
    if skip_nonidentifiable and 'identifiable' in d:
        d = d[d['identifiable'].astype(str) == 'True']
    key = d[['VOI', 'model']].astype(str).agg('/'.join, axis=1) if {'VOI', 'model'} <= set(d) else d.index.astype(str)
    out = {}
    for c in d.select_dtypes('number'):
        for k, v in zip(key, d[c]): out[f'{k}/{c}'] = float(v)
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--strict', action='store_true', help='rtol 1e-6 and no exclusions (same platform)')
    ap.add_argument('--tol', type=float, help='override the relative tolerance')
    a = ap.parse_args()
    tol = a.tol if a.tol is not None else (1e-6 if a.strict else 1e-3)
    excl = {} if a.strict else json.load(open(os.path.join(PUB, 'COMPARISON.json'), encoding='utf-8'))['exclude']
    bad_files = 0; n = 0; bad_values = []
    for pub in sorted(glob.glob(os.path.join(PUB, '**', '*.*'), recursive=True)):
        rel = os.path.relpath(pub, PUB).replace(os.sep, '/')
        if rel in ('SHA256SUMS', 'tables_published.json', 'COMPARISON.json') or not rel.endswith(('.json', '.csv')):
            continue
        new = os.path.join(ROOT, 'results', rel)
        if not os.path.exists(new):
            print(f'MISSING  {rel}'); bad_files += 1; continue
        A, B = load(new, not a.strict), load(pub, not a.strict)
        pats = [re.compile(p) for p in excl.get(rel, [])]
        worst = 0.0; fails = []
        for k, y in B.items():
            if any(p.search(k) for p in pats): continue
            if k not in A: fails.append((k, np.nan, y, np.inf)); continue
            x = A[k]
            if np.isnan(x) and np.isnan(y): continue
            d = abs(x - y) / abs(y) if abs(y) > 1e-9 else abs(x - y)
            worst = max(worst, d)
            if not d <= tol: fails.append((k, x, y, d))
        n += 1
        print(f'{"OK  " if not fails else "DIFF"} {rel:40s} max rel diff {worst:.2e}' + (f'  ({len(fails)} value(s) > {tol:g})' if fails else ''))
        if fails:
            bad_files += 1; bad_values += [(rel,) + f for f in fails]
    if bad_values:
        print(f'\nValues outside the tolerance (rtol {tol:g}):')
        for rel, k, x, y, d in sorted(bad_values, key=lambda r: -r[4])[:40]:
            print(f'  {rel}{k}: run {x:.6g} | published {y:.6g} | rel diff {d:.2e}')
    print(f'\n{n} files compared (rtol {tol:g}{", strict" if a.strict else ", cross-platform"}); '
          + ('all within tolerance' if bad_files == 0 else f'{bad_files} file(s) differ or missing'))
    sys.exit(1 if bad_files else 0)

if __name__ == '__main__':
    main()
