#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sphere model for the partial-volume effect (PVE) of a uniform sphere in a uniform background
==========================================================================================

Reproduces the sphere model of:
  "Isolating the effect of spatial resolution on absorbed-dose estimation in 177Lu-DOTATATE
   dosimetry: a noise-free digital phantom study with comparison to measured SPECT/CT"

What the model does
-------------------
1. A uniform sphere (volume V) is placed in a uniform background; concentration ratio
   sphere : background = RATIO (default 16.8; tumour-to-liver ratio of the time-integrated
   concentrations of the ground truth).
2. The object is voxelised on a cubic grid (voxel size VX) in a partial-volume representation:
   each voxel value = background + (RATIO - 1) * occupancy, where the occupancy is estimated
   with SUB x SUB x SUB sub-samples per voxel.
3. The image is blurred with an isotropic 3D discrete Gaussian kernel (modified Bessel functions,
   variance preserved exactly; Lindeberg 1990) for each FWHM. The background is treated as
   infinite (edge replication), so the box size does not influence the result.
4. The VOI is the set of voxels with occupancy >= 0.5 (binary VOI).
5. RC(FWHM) = mean(blurred image in VOI) / mean(unblurred image in VOI).
6. Steps 2-5 are repeated for N_OFFSETS random sub-voxel positions of the sphere centre
   (fixed random seed) and averaged.
7. The largest FWHM compatible with RC >= 0.95, 0.90 and 0.80 is obtained by linear interpolation.

Sub-commands
------------
  simulate    RC versus FWHM for a set of sphere volumes on one grid  -> CSV
  thresholds  largest FWHM for RC >= 0.95 / 0.90 / 0.80 from one or more CSV files -> CSV
  validate    compare a sphere (e.g. 15.9 mL) with observed RC values (e.g. phantom tumour)
  figure      three-panel figure: (a) example slices, (b) RC curves, (c) FWHM limits

Example (reproduces the manuscript)
-----------------------------------
  python sphere_model.py simulate   --voxel 4.419644 --out rc_grid_4.42mm.csv
  python sphere_model.py simulate   --voxel 2.209822 --out rc_grid_2.21mm.csv
  python sphere_model.py thresholds rc_grid_4.42mm.csv rc_grid_2.21mm.csv --out fwhm_limits.csv
  python sphere_model.py validate   --voxel 4.419644 --volume 15.9 --observed observed_tumour_rc.csv
  python sphere_model.py figure     --coarse rc_grid_4.42mm.csv --fine rc_grid_2.21mm.csv \
                                    --observed observed_tumour_rc.csv --targets phantom_targets.csv \
                                    --out sphere_model_figure.png

Requirements: Python >= 3.9, numpy, scipy, pandas, matplotlib
"""
import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.special import ive
from scipy.ndimage import convolve1d

DEFAULT_VOLUMES = [2, 4, 8, 16, 32, 64, 128]                       # mL
DEFAULT_FWHM = [0.0, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0, 11.25, 12.5, 13.75, 15.0]  # mm
DEFAULT_RATIO = 16.8
DEFAULT_SUB = 6
DEFAULT_OFFSETS = 32
DEFAULT_SEED = 0
LEVELS = (0.95, 0.90, 0.80)
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


# ----------------------------------------------------------------------------- core functions
def discrete_gaussian_kernel(sigma_vox):
    """1D discrete Gaussian T(n; t) = exp(-t) I_n(t), t = sigma^2 (in voxels). Sums to 1."""
    t = sigma_vox ** 2
    half = int(np.ceil(5.0 * sigma_vox)) + 3
    n = np.arange(-half, half + 1)
    w = ive(np.abs(n), t)                      # exponentially scaled Bessel: exp(-t) I_n(t)
    return w / w.sum()


def blur(img, fwhm_mm, voxel_mm):
    """Separable 3D discrete Gaussian blur; background replicated at the box edges."""
    if fwhm_mm <= 0:
        return img
    k = discrete_gaussian_kernel(fwhm_mm * FWHM_TO_SIGMA / voxel_mm)
    out = img
    for ax in range(3):
        out = convolve1d(out, k, axis=ax, mode='nearest')
    return out


def sphere_occupancy(volume_ml, voxel_mm, n_box, centre_offset_vox, sub=DEFAULT_SUB):
    """Occupancy (0..1) of each voxel by a sphere of the given volume.
    Voxel i covers [i, i+1) in voxel units; the sphere centre is at n_box/2 + offset."""
    r = (3.0 * volume_ml * 1000.0 / (4.0 * np.pi)) ** (1.0 / 3.0) / voxel_mm   # radius in voxels
    c = n_box / 2.0 + np.asarray(centre_offset_vox, float)
    occ = np.zeros((n_box, n_box, n_box))
    lo = int(np.floor(c.min() - r - 1)); hi = int(np.ceil(c.max() + r + 2))
    lo = max(lo, 0); hi = min(hi, n_box)
    idx = np.arange(lo, hi)
    s = (np.arange(sub) + 0.5) / sub
    g = (idx[:, None] + s[None, :]).ravel()
    X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
    inside = ((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2) <= r * r
    m = len(idx)
    occ[lo:hi, lo:hi, lo:hi] = inside.reshape(m, sub, m, sub, m, sub).mean(axis=(1, 3, 5))
    return occ, r


def box_size(volume_ml, voxel_mm, fwhm_max):
    """Box large enough for the sphere plus the kernel support (edge replication is used anyway)."""
    r_mm = (3.0 * volume_ml * 1000.0 / (4.0 * np.pi)) ** (1.0 / 3.0)
    margin = 5.0 * fwhm_max * FWHM_TO_SIGMA + 3 * voxel_mm
    return int(np.ceil(2 * (r_mm + margin) / voxel_mm)) + 4


def rc_curve(volume_ml, voxel_mm, fwhms, ratio=DEFAULT_RATIO, n_offsets=DEFAULT_OFFSETS, sub=DEFAULT_SUB, rng=None):
    """Mean RC over random sub-voxel positions; also returns the 0-mm VOI bias (voxelisation)."""
    rng = np.random.default_rng(DEFAULT_SEED) if rng is None else rng
    n = box_size(volume_ml, voxel_mm, max(fwhms))
    rcs, base = [], []
    for _ in range(n_offsets):
        occ, _r = sphere_occupancy(volume_ml, voxel_mm, n, rng.uniform(-0.5, 0.5, 3), sub)
        img = 1.0 + (ratio - 1.0) * occ                 # background = 1, sphere = ratio
        voi = occ >= 0.5
        m0 = img[voi].mean()
        base.append(m0 / ratio - 1.0)                   # VOI mean vs true sphere concentration
        rcs.append([blur(img, f, voxel_mm)[voi].mean() / m0 for f in fwhms])
    rcs = np.array(rcs)
    return rcs.mean(0), rcs.std(0), float(np.mean(base)), n


def fwhm_limit(fwhms, rc, level):
    """Largest FWHM with RC >= level (linear interpolation); None if never reached."""
    fwhms = np.asarray(fwhms, float); rc = np.asarray(rc, float)
    if rc.min() > level:
        return None
    return float(np.interp(-level, -rc, fwhms))         # rc decreases monotonically


# ----------------------------------------------------------------------------- sub-commands
def cmd_simulate(a):
    fwhms = [float(x) for x in a.fwhm]
    rng = np.random.default_rng(a.seed)
    rows = []
    print(f'voxel {a.voxel} mm | ratio {a.ratio} | sub-samples {a.sub}^3 | offsets {a.offsets} | seed {a.seed}')
    for v in a.volumes:
        t0 = time.time()
        rc, sd, base, n = rc_curve(float(v), a.voxel, fwhms, a.ratio, a.offsets, a.sub, rng)
        for f, r, s in zip(fwhms, rc, sd):
            rows.append(dict(voxel_mm=a.voxel, volume_ml=float(v), fwhm_mm=f, rc=r, rc_sd_offsets=s,
                             baseline_0mm=base, ratio=a.ratio, sub=a.sub, n_offsets=a.offsets, seed=a.seed, box=n))
        print(f'  {v:>6} mL  box {n:3d}^3  0-mm VOI bias {base * 100:+5.1f}%  RC(15 mm) {rc[-1]:.3f}  ({time.time() - t0:.1f} s)')
    df = pd.DataFrame(rows)
    if a.append and os.path.exists(a.out):
        df = pd.concat([pd.read_csv(a.out), df], ignore_index=True)
    df.to_csv(a.out, index=False)
    print(f'saved {a.out}')


def cmd_thresholds(a):
    out = []
    for f in a.csv:
        d = pd.read_csv(f)
        for (vx, v), g in d.groupby(['voxel_mm', 'volume_ml']):
            g = g.sort_values('fwhm_mm')
            row = dict(voxel_mm=vx, volume_ml=v)
            for lv in LEVELS:
                row[f'fwhm_rc_{lv:.2f}'] = fwhm_limit(g.fwhm_mm, g.rc, lv)
            out.append(row)
    df = pd.DataFrame(out).sort_values(['voxel_mm', 'volume_ml'], ascending=[False, True])
    df.to_csv(a.out, index=False)
    with pd.option_context('display.float_format', '{:.1f}'.format):
        print(df.to_string(index=False))
    print(f'saved {a.out}')


def cmd_validate(a):
    obs = pd.read_csv(a.observed).sort_values('fwhm_mm')
    fwhms = obs.fwhm_mm.astype(float).tolist()
    rc, sd, base, n = rc_curve(a.volume, a.voxel, fwhms, a.ratio, a.offsets, a.sub, np.random.default_rng(a.seed))
    diff = rc - obs.rc.values
    print('FWHM (mm)  sphere RC  observed RC  difference')
    for f, r, o, d in zip(fwhms, rc, obs.rc.values, diff):
        print(f'{f:9.2f}  {r:9.3f}  {o:11.3f}  {d:+.3f}')
    print(f'max |difference| = {np.abs(diff).max():.4f}  (sphere {a.volume} mL, ratio {a.ratio}, voxel {a.voxel} mm)')


def cmd_figure(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.6, 'legend.frameon': False})
    mm = 1 / 25.4
    co = pd.read_csv(a.coarse); fi = pd.read_csv(a.fine) if a.fine else None
    vx = float(co.voxel_mm.iloc[0]); ratio = float(co.ratio.iloc[0])
    fig = plt.figure(figsize=(174 * mm, 108 * mm))
    gs = GridSpec(2, 2, height_ratios=[0.42, 1.0], hspace=0.30, wspace=0.28, left=0.07, right=0.99, top=0.90, bottom=0.09)
    # (a) example slices
    gsa = gs[0, :].subgridspec(1, 9, wspace=0.06)
    vols, fws, half = [4, 16, 64], [0.0, 7.5, 15.0], 8
    k = 0
    for V in vols:
        n = box_size(V, vx, 15.0)
        occ, r = sphere_occupancy(V, vx, n, np.zeros(3))
        img = 1 + (ratio - 1) * occ; voi = occ >= 0.5; c0 = n // 2
        sl = np.s_[c0, c0 - half:c0 + half + 1, c0 - half:c0 + half + 1]
        for f in fws:
            ax = fig.add_subplot(gsa[0, k]); k += 1
            ext = [-half * vx - vx / 2, half * vx + vx / 2, half * vx + vx / 2, -half * vx - vx / 2]
            ax.imshow(blur(img, f, vx)[sl], cmap='gray', vmin=0, vmax=ratio, extent=ext, interpolation='nearest')
            gx = np.arange(-half, half + 1) * vx
            ax.contour(gx, gx, voi[sl].astype(float), levels=[0.5], colors=['#F08A24'], linewidths=0.9)
            off = (n / 2 - (c0 + 0.5)) * vx; th = np.linspace(0, 2 * np.pi, 200)
            ax.plot(r * vx * np.cos(th) + off, r * vx * np.sin(th) + off, '--', color='#7FD3FF', lw=0.7)
            ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
            ax.set_title(f'{f:g} mm', fontsize=7, pad=2)
    fig.text(0.02, 0.965, '(a)', fontsize=9, fontweight='bold')
    for gi, V in enumerate(vols):
        fig.text(0.07 + (gi * 3 + 1.5) / 9 * 0.92, 0.965, f'{V} mL sphere', ha='center', fontsize=7.5, fontweight='bold')
    # (b) RC curves
    ax = fig.add_subplot(gs[1, 0]); cm = plt.cm.viridis
    vv = sorted(co.volume_ml.unique())
    for i, V in enumerate(vv):
        g = co[co.volume_ml == V].sort_values('fwhm_mm')
        ax.plot(g.fwhm_mm, g.rc, color=cm(i / max(len(vv) - 1, 1)), lw=1.1, label=f'{V:g} mL')
    if a.observed:
        ob = pd.read_csv(a.observed)
        ax.plot(ob.fwhm_mm, ob.rc, 'o', color='#D9731E', ms=3.5, label=a.observed_label)
    ax.set_xlabel('3D Gaussian FWHM (mm)'); ax.set_ylabel('RC (sphere in uniform background)')
    ax.set_ylim(0.3, 1.02); ax.legend(fontsize=6, loc='lower left', ncol=2)
    ax.text(-0.17, 1.04, '(b)', transform=ax.transAxes, fontsize=9, fontweight='bold')
    # (c) FWHM limits
    ax = fig.add_subplot(gs[1, 1])
    def lim(d, lv):
        return [fwhm_limit(g.sort_values('fwhm_mm').fwhm_mm, g.sort_values('fwhm_mm').rc, lv) for _, g in d.groupby('volume_ml')]
    x = np.array(sorted(co.volume_ml.unique()))
    for lv, ls in [(0.95, ':'), (0.90, '-'), (0.80, '--')]:
        y = lim(co, lv); m = [i for i, t in enumerate(y) if t is not None]
        ax.plot(x[m], [y[i] for i in m], ls, color='k', lw=1.1, marker='o', ms=2.5, label=f'RC ≥ {lv:.2f} ({vx:.2f}-mm grid)')
    if fi is not None:
        xf = np.array(sorted(fi.volume_ml.unique())); y = lim(fi, 0.90)
        ax.plot(xf, y, '-', color='#9AA9B2', lw=1.1, marker='s', ms=2.5, label=f'RC ≥ 0.90 ({float(fi.voxel_mm.iloc[0]):.2f}-mm grid)')
    if a.targets:
        tg = pd.read_csv(a.targets)
        cols = ['#D9731E', '#1F8A94', '#2F4F9E', '#7C99D4', '#8A7650']
        for i, r_ in tg.iterrows():
            ax.plot(r_.volume_ml, r_.fwhm_rc_0_90, 'D', color=cols[i % len(cols)], ms=4, label=f'{r_.target} (RC ≥ 0.90)')
    ax.set_xscale('log'); ax.set_xticks(x); ax.set_xticklabels([f'{v:g}' for v in x])
    ax.set_xlabel('Sphere volume (mL)'); ax.set_ylabel('Largest FWHM (mm)'); ax.set_ylim(0, 16); ax.legend(fontsize=5.6, loc='upper left')
    ax.text(-0.17, 1.04, '(c)', transform=ax.transAxes, fontsize=9, fontweight='bold')
    for axx in fig.axes[9:]:
        for s in ['top', 'right']: axx.spines[s].set_visible(False)
    fig.savefig(a.out, dpi=a.dpi, facecolor='white')
    print(f'saved {a.out}')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('simulate', help='RC versus FWHM for sphere volumes on one grid')
    s.add_argument('--voxel', type=float, required=True, help='voxel size (mm), e.g. 4.419644 or 2.209822')
    s.add_argument('--volumes', type=float, nargs='+', default=DEFAULT_VOLUMES)
    s.add_argument('--fwhm', type=float, nargs='+', default=DEFAULT_FWHM)
    s.add_argument('--ratio', type=float, default=DEFAULT_RATIO)
    s.add_argument('--sub', type=int, default=DEFAULT_SUB, help='sub-samples per voxel and axis')
    s.add_argument('--offsets', type=int, default=DEFAULT_OFFSETS, help='random sub-voxel positions')
    s.add_argument('--seed', type=int, default=DEFAULT_SEED)
    s.add_argument('--out', required=True)
    s.add_argument('--append', action='store_true', help='append to an existing CSV (to split long runs)')
    t = sp.add_parser('thresholds', help='largest FWHM for RC >= 0.95/0.90/0.80')
    t.add_argument('csv', nargs='+'); t.add_argument('--out', default='fwhm_limits.csv')
    v = sp.add_parser('validate', help='compare one sphere with observed RC values')
    v.add_argument('--voxel', type=float, default=4.419644); v.add_argument('--volume', type=float, required=True)
    v.add_argument('--observed', required=True, help='CSV with columns fwhm_mm, rc')
    v.add_argument('--ratio', type=float, default=DEFAULT_RATIO); v.add_argument('--sub', type=int, default=DEFAULT_SUB)
    v.add_argument('--offsets', type=int, default=DEFAULT_OFFSETS); v.add_argument('--seed', type=int, default=DEFAULT_SEED)
    f = sp.add_parser('figure', help='three-panel figure')
    f.add_argument('--coarse', required=True); f.add_argument('--fine')
    f.add_argument('--observed', help='CSV fwhm_mm, rc (e.g. phantom tumour)'); f.add_argument('--observed-label', default='Phantom tumour (15.9 mL)')
    f.add_argument('--targets', help='CSV target, volume_ml, fwhm_rc_0_90'); f.add_argument('--out', default='sphere_model_figure.png')
    f.add_argument('--dpi', type=int, default=600)
    a = p.parse_args()
    print(f'Python {sys.version.split()[0]} | numpy {np.__version__} | pandas {pd.__version__}')
    {'simulate': cmd_simulate, 'thresholds': cmd_thresholds, 'validate': cmd_validate, 'figure': cmd_figure}[a.cmd](a)


if __name__ == '__main__':
    main()
