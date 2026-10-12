#!/usr/bin/env python
'''
Electron-scattering polarization of the winds of Sirocco models (Step 4),
run on a directory of completed models such as the examples/regress set.


Command line usage (if any):

    usage::

        step4_regress.py [-h] [-n NPHOT] [-n_bm NPHOT] [-seed SEED] [-no_wind]
                         [-plot FILE] [-out FILE] [-from_table FILE]
                         [directory] [root ...]

    directory   where the models were run (default: the current directory)
    root ...    the models to use (default: every root with a .wind_save file)
    -n          photons per model with the model's own sources (default 1e5)
    -n_bm       photons per model for the point-source check (default 2e5)
    -seed       seed for the random number generator
    -no_wind    leave out the emission of the wind
    -plot       write a summary plot to FILE
    -out        write the results table to FILE (default step4_results.txt)
    -from_table only make the plot (-plot FILE), from a results table saved earlier

Description:

    Each model needs root.stdout.txt (for the source luminosities), root.pf
    or root.out.pf, and the windsave2table master and heat tables (run
    windsave2table on each root first; the heat tables must have the
    kpkt_lum and matom_lum columns for macro-atom models).

    For each model two calculations are made with sirocco_grid:

    1. The model's own sources (central object, boundary layer, AGN, disk
       and wind) weighted by their luminosities in the band of the spectral
       cycles, peeled off toward the model's observers.  Only electron
       scattering is included.  The wind emits unpolarized light, which is
       then scattered like the rest.
    2. A check: an isotropic point source at the center with no disk or
       central object to block it, compared with Brown & McLean (1977)
       evaluated directly from the grid cells.  The two should agree when the
       wind is optically thin to electron scattering; spherical winds must
       give no net polarization.

Primary routines:

    main, do_model, plot

Notes:

    Requires numpy and astropy, and matplotlib for -plot.

History:

    261010 ksl Coding begun, with Claude Code
    261011 ksl Wind emission added; source weights from the band (f_) luminosities

'''

import argparse
import os
import sys
import time
from glob import glob

import numpy as np

# Make the polarization package importable when this is run as a script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from polarization.sirocco_grid import SiroccoModel, run, brown_mclean_grid
from polarization.stokes import pol_degree_angle

BM_ANGLES = (30, 60, 90)


def do_model(root, directory, nphot, nphot_bm, rng, wind=True):
    '''Run one model; returns a dict of results.'''
    t0 = time.time()
    model = SiroccoModel(root, directory, wind=wind)
    out = {'root': root, 'coords': [d.coord for d in model.domains], 'sources': model.sources,
           'obs': [], 'bm': []}

    for s in run(model, nphot, rng):
        q, u, sq, su = s.result()
        P, pa = pol_degree_angle(1.0, q, u)
        out['obs'].append((s.angle, s.I, s.I_scat / s.I, q, sq, u, su, P, pa))

    check = SiroccoModel(root, directory, point_source=True)
    for s in run(check, nphot_bm, rng, angles=BM_ANGLES):
        q, u, sq, su = s.result()
        qb, tau0, gamma = brown_mclean_grid(check, s.angle)
        out['bm'].append((s.angle, q, sq, u, su, qb, tau0, gamma))
    out['time'] = time.time() - t0
    return out


def report(res, f):
    def both(line=''):
        print(line)
        f.write(line + '\n')

    r = res
    srcs = ', '.join(f'{k} {v:.3f}' for k, v in r['sources'].items())
    both(f"\n{r['root']}  (grid: {', '.join(r['coords'])}; source weights: {srcs}; "
         f"{r['time']:.0f} s)")
    both('   angle   I        scattered     q (%)            u (%)          P (%)    PA (deg)')
    for ang, I, fs, q, sq, u, su, P, pa in r['obs']:
        both(f'   {ang:5.1f}  {I:.4f}   {fs:6.3f}    {100 * q:+7.3f}±{100 * sq:.3f}   {100 * u:+7.3f}±{100 * su:.3f}   '
             f'{100 * P:6.3f}   {pa:6.1f}')
    tau0, gamma = r['bm'][0][6], r['bm'][0][7]
    both(f'   point-source check (Brown & McLean from the grid: tau0 = {tau0:.3g}, gamma = {gamma:.3f}):')
    for ang, q, sq, u, su, qb, _, _ in r['bm']:
        ratio = f'{q / qb:6.3f}' if abs(qb) > 1e-12 else '   -  '
        both(f'      i = {ang:3.0f}:  q = {100 * q:+7.3f}±{100 * sq:.3f}%   B&M {100 * qb:+7.3f}%   ratio {ratio}'
             f'   u/sig {u / su:+5.1f}')


def plot(results, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    colors = [plt.cm.tab10(k) for k in range(10)] + ['black', 'gold']
    markers = 'osD^v<>ph*'
    for k, r in enumerate(results):
        o = np.array(r['obs'])
        order = np.argsort(o[:, 0])
        ax[0].errorbar(o[order, 0] + 0.4 * (k - len(results) / 2), 100 * o[order, 3], 100 * o[order, 4],
                       fmt='-', marker=markers[k % len(markers)], color=colors[k % len(colors)], ms=5,
                       label=r['root'])
    ax[0].axhline(0, color='k', lw=0.8)
    ax[0].set_yscale('symlog', linthresh=0.1)
    ax[0].set_xlabel('inclination (deg)')
    ax[0].set_ylabel('q = Q/I (%)   [Q > 0: E-vector along the projected axis]')
    ax[0].set_title('Electron-scattering polarization, model sources')
    ax[0].legend(fontsize=7, ncol=2)

    for k, r in enumerate(results):
        b = np.array(r['bm'])
        tau0 = b[0, 6]
        ok = np.abs(b[:, 5]) > 1e-12
        if ok.any():
            ax[1].errorbar(np.full(ok.sum(), tau0) * (1 + 0.03 * (b[ok, 0] - 60) / 30), b[ok, 1] / b[ok, 5],
                           b[ok, 2] / np.abs(b[ok, 5]), fmt=markers[k % len(markers)],
                           color=colors[k % len(colors)], ms=5, label=r['root'])
    ax[1].axhline(1, color='k', lw=0.8)
    ax[1].set_xscale('log')
    ax[1].set_xlabel('tau0 (Brown & McLean, from the grid)')
    ax[1].set_ylabel('q (point source) / q (Brown & McLean)')
    ax[1].set_title('Point-source check (i = 30, 60, 90)')
    ax[1].legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    print(f'\nPlot written to {fname}')


def read_table(fname):
    '''Read the results back from a table written by report(), for replotting.'''
    results = []
    r = None
    for line in open(fname):
        w = line.split()
        if not w or line.startswith('#'):
            continue
        if line[0] not in ' ' and '(grid:' in line:
            r = {'root': w[0], 'obs': [], 'bm': []}
            results.append(r)
        elif r is not None and w[0].replace('.', '', 1).isdigit() and len(w) >= 7:
            ang, I, fs = float(w[0]), float(w[1]), float(w[2])
            q, sq = (float(v) / 100 for v in w[3].split('±'))
            u, su = (float(v) / 100 for v in w[4].split('±'))
            r['obs'].append((ang, I, fs, q, sq, u, su, float(w[5]) / 100, float(w[6])))
        elif r is not None and line.strip().startswith('point-source check'):
            tau0 = float(line.split('tau0 =')[1].split(',')[0])
            gamma = float(line.split('gamma =')[1].split(')')[0])
        elif r is not None and line.strip().startswith('i ='):
            ang = float(line.split('i =')[1].split(':')[0])
            q, sq = (float(v.strip('%')) / 100 for v in line.split('q =')[1].split()[0].split('±'))
            qb = float(line.split('B&M')[1].split('%')[0]) / 100
            entry = (ang, q, sq, 0.0, 1.0, qb, tau0, gamma)
            if entry not in r['bm']:
                r['bm'].append(entry)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('directory', nargs='?', default='.')
    parser.add_argument('roots', nargs='*')
    parser.add_argument('-n', type=float, default=1e5)
    parser.add_argument('-n_bm', type=float, default=2e5)
    parser.add_argument('-seed', type=int, default=20261010)
    parser.add_argument('-no_wind', action='store_true', help='leave out the emission of the wind')
    parser.add_argument('-plot', default=None)
    parser.add_argument('-out', default='step4_results.txt')
    parser.add_argument('-from_table', default=None, help='only make the plot, from a saved results table')
    args = parser.parse_args()

    if args.from_table:
        plot(read_table(args.from_table), args.plot or 'step4_regress.png')
        return

    roots = args.roots or sorted(os.path.basename(f)[:-10] for f in glob(os.path.join(args.directory, '*.wind_save')))
    rng = np.random.default_rng(args.seed)
    results = []
    with open(args.out, 'w') as f:
        f.write(f'# step4_regress.py on {os.path.abspath(args.directory)}: n = {int(args.n)}, n_bm = {int(args.n_bm)}, '
                f'wind emission {"excluded" if args.no_wind else "included"}\n')
        for root in roots:
            try:
                res = do_model(root, args.directory, int(args.n), int(args.n_bm), rng, wind=not args.no_wind)
            except Exception as e:
                print(f'\n{root}: failed ({e})')
                f.write(f'\n{root}: failed ({e})\n')
                continue
            report(res, f)
            f.flush()
            results.append(res)
    if args.plot and results:
        plot(results, args.plot)


if __name__ == '__main__':
    main()
