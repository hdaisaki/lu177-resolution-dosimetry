#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time-integration of VOI absorbed-dose-rate curves with NUKFIT-type sum-of-exponential functions.

Input
  One 3D Slicer "Segment statistics" table (.tsv) per time point, exported from the OpenDose3D
  dose-rate map (FFT convolution).  The column "Mean" (VOI mean absorbed-dose rate, mGy/h) is used.
  Files are sorted by the number before "HR" in the file name and assigned, in that order, to the
  acquisition times given with --times (h after administration).  The mapping is printed.

Functions (x = time in h, lambda_p = physical decay constant of Lu-177)
  f2 : y = A1*exp(-(l1+lp)x)
  f3 : y = A1*exp(-(l1+lp)x) - A1*exp(-(l2+lp)x)
  f4b: y = A1*exp(-(l1+lp)x) - A2*exp(-(l2+lp)x) - (A1-A2)*exp(-(lbc+lp)x),  lbc fixed (default 41.59 /h)

Fitting
  Weighted least squares with a proportional, data-based error model (sigma_i proportional to y_i),
  i.e. minimise sum(((f(t_i)-y_i)/y_i)^2).  Bounded trust-region-reflective optimiser
  (scipy.optimize.least_squares, method 'trf').  To avoid local minima, the fit is repeated from
  N_STARTS starting points (the nominal start values multiplied by log-normal factors drawn with a
  fixed seed per start) and the solution with the lowest cost is kept.  Data are scaled to max = 1
  for fitting; amplitudes are scaled back afterwards.

Integration
  Analytical integral from 0 to infinity (mGy -> Gy).  Extrapolation (%) = share of the integral
  after the last time point.  A model-free reference (trapezoid, linear from 0 to t1, mono-exponential
  tail fitted to the last three points with lambda >= lambda_p) is reported as a check.

Diagnostics
  max relative residual, effective half-life between 72 and 144 h, AICc (n ln(RSS/n)+2k+2k(k+1)/(n-k-1)),
  approximate parameter CV from the Jacobian (residual variance RSS/(n-k)), CV of the integral by
  the delta method, and a flag when a parameter sits on a bound.

Usage
  python fit_tac_nukfit_like.py --dir ./FWHM_00mm --times 1 4 24 40 72 144 --out results_00mm.csv
  (requires numpy, scipy, pandas)

Fits with a parameter coefficient of variation above CV_LIMIT (1000%) are flagged identifiable=False; their SD and CV
of the time integral are reported as NaN, because the covariance matrix is nearly singular and these values then
depend on floating-point details of the platform (in the article, only unused fits are affected).
"""
import argparse, glob, os, re, sys, platform
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

LN2 = np.log(2.0)
T_HALF_LU177_H = 6.647 * 24.0            # 159.53 h  -> lambda_p = 4.345e-3 /h (as in NUKFIT output)
LP = LN2 / T_HALF_LU177_H
LBC = 41.59                               # fixed fast rate constant of f4b (value used by NUKFIT)
N_STARTS = 250
SPREAD = 1.2                              # SD of the log-normal perturbation of start values
CV_LIMIT = 1000.0                         # % ; a parameter CV above this marks a non-identifiable (degenerate) fit
DEFAULT_MODELS = {'Tumour': 'f4b', 'spleen': 'f4b', 'Right Kidney': 'f2', 'Left Kidney': 'f2', 'Liver': 'f3'}


# ---------------------------------------------------------------- functions
def f2(t, A1, l1):
    return A1 * np.exp(-(l1 + LP) * t)


def f3(t, A1, l1, l2):
    return A1 * (np.exp(-(l1 + LP) * t) - np.exp(-(l2 + LP) * t))


def f4b(t, A1, l1, A2, l2):
    return A1 * np.exp(-(l1 + LP) * t) - A2 * np.exp(-(l2 + LP) * t) - (A1 - A2) * np.exp(-(LBC + LP) * t)


# name: (function, start values, lower bounds, upper bounds, indices of amplitude parameters)
MODELS = {
    'f2':  (f2,  [1.0, 5e-3],             [0, 0],       [1e4, 1],            [0]),
    'f3':  (f3,  [1.1, 2e-3, 2.0],        [0, 0, 0],    [1e4, 1, 100],       [0]),
    'f4b': (f4b, [1.3, 3e-3, 0.3, 0.04],  [0, 0, 0, 0], [1e5, 1, 1e5, 10],   [0, 2]),
}


def integral(model, p, t0=0.0):
    """Analytical integral of the model from t0 to infinity (same units as y * h)."""
    e = lambda A, k: A * np.exp(-k * t0) / k
    if model == 'f2':
        A1, l1 = p; return e(A1, l1 + LP)
    if model == 'f3':
        A1, l1, l2 = p; return e(A1, l1 + LP) - e(A1, l2 + LP)
    A1, l1, A2, l2 = p
    return e(A1, l1 + LP) - e(A2, l2 + LP) - e(A1 - A2, LBC + LP)


def equation(model, p):
    s = lambda v: f'{v:.3e}'
    if model == 'f2':
        return f'f2: y = {s(p[0])}*exp(-({s(p[1])}+{LP:.3e})*x)'
    if model == 'f3':
        return f'f3: y = {s(p[0])}*exp(-({s(p[1])}+{LP:.3e})*x) - {s(p[0])}*exp(-({s(p[2])}+{LP:.3e})*x)'
    return (f'f4b: y = {s(p[0])}*exp(-({s(p[1])}+{LP:.3e})*x) - {s(p[2])}*exp(-({s(p[3])}+{LP:.3e})*x)'
            f' - ({s(p[0])}-{s(p[2])})*exp(-({LBC:.3e}+{LP:.3e})*x)')


# ---------------------------------------------------------------- fitting
def fit(model, t, y):
    f, p0, lb, ub, amp_idx = MODELS[model]
    lb, ub = np.array(lb, float), np.array(ub, float)
    scale = y.max()
    yy = y / scale
    resid = lambda q: (f(t, *q) - yy) / yy              # proportional, data-based weighting
    best = None
    for k in range(N_STARTS):
        if k == 0:
            start = np.array(p0, float)
        else:
            rng = np.random.default_rng(k)               # fixed seed per start -> reproducible
            start = np.array(p0, float) * np.exp(rng.normal(0.0, SPREAD, len(p0)))
        start = np.clip(start, lb + 1e-12, ub - 1e-9)
        r = least_squares(resid, start, bounds=(lb, ub), method='trf', max_nfev=8000)
        if best is None or r.cost < best.cost:
            best = r
    p = best.x.copy()
    n, kpar = len(t), len(p)
    rss = 2.0 * best.cost
    # approximate parameter covariance (relative weighting: residual variance estimated from RSS)
    cv = np.full(kpar, np.nan)
    ad_cv = np.nan
    if n > kpar:
        J = best.jac
        try:
            cov = np.linalg.inv(J.T @ J) * rss / (n - kpar)
            with np.errstate(divide='ignore', invalid='ignore'):
                cv = np.sqrt(np.diag(cov)) / np.abs(p) * 100.0
            # uncertainty of the integral (delta method, numerical gradient in the scaled parameter space)
            I0 = integral(model, p)
            g = np.zeros(kpar)
            for i in range(kpar):
                h = 1e-6 * max(abs(p[i]), 1e-12)
                q = p.copy(); q[i] += h
                g[i] = (integral(model, q) - I0) / h
            ad_cv = float(np.sqrt(max(g @ cov @ g, 0.0)) / abs(I0) * 100.0)
        except np.linalg.LinAlgError:
            pass
    at_bound = bool(np.any(np.isclose(p, lb, rtol=0, atol=1e-10) | np.isclose(p, ub, rtol=1e-6)))
    for i in amp_idx:
        p[i] *= scale
    aicc = n * np.log(rss / n) + 2 * kpar + (2 * kpar * (kpar + 1) / (n - kpar - 1) if n - kpar - 1 > 0 else np.inf)
    return p, rss, aicc, cv, at_bound, ad_cv


def reference_integral(t, y):
    lam = max(-np.polyfit(t[-3:], np.log(y[-3:]), 1)[0], LP)
    return 0.5 * y[0] * t[0] + np.trapezoid(y, t) + y[-1] / lam


# ---------------------------------------------------------------- I/O
def read_tables(folder, times):
    files = glob.glob(os.path.join(folder, '*.tsv'))
    keyed = []
    for fpath in files:
        m = re.search(r'(\d+(?:\.\d+)?)\s*HR', os.path.basename(fpath), flags=re.I)
        if m:
            keyed.append((float(m.group(1)), fpath))
    keyed.sort()
    if len(keyed) != len(times):
        sys.exit(f'Found {len(keyed)} tables with "<n>HR" in the name in {folder}, but {len(times)} times were given.')
    cols = []
    print('File -> acquisition time used')
    for (label, fpath), t in zip(keyed, times):
        print(f'  {os.path.basename(fpath):45s} (label {label:g} h) -> {t:g} h')
        d = pd.read_csv(fpath, sep='\t').set_index('Segment')['Mean'].astype(float)
        cols.append(d.rename(t))
    return pd.concat(cols, axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True, help='folder with the six Segment-statistics .tsv files')
    ap.add_argument('--times', nargs='+', type=float, default=[1, 4, 24, 40, 72, 144], help='acquisition times (h)')
    ap.add_argument('--out', default='tac_fit_results.csv')
    ap.add_argument('--all-models', action='store_true', help='fit f2, f3 and f4b for every VOI')
    a = ap.parse_args()

    print(f'Python {platform.python_version()} | numpy {np.__version__} | scipy {scipy.__version__} | pandas {pd.__version__}')
    print(f'lambda_p = {LP:.4e} /h (T1/2 = {T_HALF_LU177_H:.2f} h), lambda_bc = {LBC} /h, starts = {N_STARTS}')
    t = np.array(a.times, float)
    D = read_tables(a.dir, t)
    rows = []
    for voi in D.index:
        y = D.loc[voi].values
        ref = reference_integral(t, y) / 1000.0
        models = ['f2', 'f3', 'f4b'] if a.all_models else [DEFAULT_MODELS.get(voi, 'f4b')]
        for m in models:
            p, rss, aicc, cv, at_bound, ad_cv = fit(m, t, y)
            identifiable = bool(np.all(np.isfinite(cv)) and np.max(cv) <= CV_LIMIT)
            if not identifiable:
                ad_cv = np.nan     # covariance nearly singular: SD/CV undefined and platform-dependent
            f = MODELS[m][0]
            ad = integral(m, p) / 1000.0
            ext = integral(m, p, t0=t[-1]) / integral(m, p) * 100.0
            teff = 72 * LN2 / np.log(f(72.0, *p) / f(144.0, *p)) if (72 in t and 144 in t) else np.nan
            maxres = np.max(np.abs(f(t, *p) / y - 1.0)) * 100.0
            rows.append(dict(VOI=voi, model=m, recommended=(m == DEFAULT_MODELS.get(voi)), AD_Gy=ad, AD_SD_Gy=ad*ad_cv/100.0, AD_CV_pct=ad_cv,
                             extrapolation_pct=ext, Teff_72_144_h=teff, max_rel_residual_pct=maxres,
                             AICc=aicc, param_CV_pct=';'.join(f'{c:.1f}' for c in cv), identifiable=identifiable, parameter_at_bound=at_bound,
                             reference_AD_Gy=ref, equation=equation(m, p),
                             params=';'.join(f'{v:.6e}' for v in p)))
            print(f'{voi:13s} {m:3s} AD {ad:8.4f} ± {ad*ad_cv/100:.4f} Gy (CV {ad_cv:4.1f}%) | extrap {ext:5.1f}% | Teff {teff:5.0f} h | maxres {maxres:4.1f}% '
                  f'| bound {at_bound} | {equation(m, p)}')
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f'saved {a.out}')


if __name__ == '__main__':
    main()
