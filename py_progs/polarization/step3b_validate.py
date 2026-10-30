#!/usr/bin/env python
'''
Test the treatments of Sobolev trapping when the escape probability
depends on direction (Step 3b).


Command line usage (if any):

    usage::

        step3b_validate.py [-h] [-n SCALE] [-seed SEED] [-plot FILE]

    -n      multiplies the default photon numbers (default 1; ~6 minutes)
    -seed   seed for the random number generator
    -plot   write a summary plot to FILE

Description:

    Uses RadialFlowEllipsoid: an oblate ellipsoid (a = 1, c = 0.3, in units
    where v0 = vmax = 0.001 c at r = a) in a radial outflow v ~ r^alpha, where
    the Sobolev optical depth is tau(mu) ~ 1/(1 - (1 - alpha) mu^2).
    alpha = 0.25 traps photons moving radially, alpha = 3 photons moving
    tangentially.

    1. alpha = 1 must reproduce the homologous flow exactly (same random
       numbers), which checks the bisection and the direction-dependent
       code paths.
    2. Scattered continuum: q of the line-scattered light for the explicit,
       single (Sirocco's current scheme) and effective (W2 F with the
       direction-averaged beta) treatments.
    3. Pure line emission: explicit trapping with anisotropic escape makes
       unpolarized line emission polarized; the single and effective
       treatments, which emit line photons unpolarized, cannot.  The size
       of the effect against tau.
    4. A sphere with pure line emission must give no net polarization.

Primary routines:

    main, case_equivalence, case_scattering, case_emission, case_sphere, plot

Notes:

    Requires numpy, and matplotlib for -plot.

History:

    261010 ksl Coding begun, with Claude Code

'''

import argparse
import os
import sys
import time

import numpy as np

# Make the polarization package importable when this is run as a script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from polarization.sobolev_flow import HomologousEllipsoid, RadialFlowEllipsoid, run

VMAX = 0.001


def geom(alpha, tau, trap, line_frac=0.0, c=0.3, w2=1.0):
    return RadialFlowEllipsoid(VMAX, c * VMAX, VMAX, alpha, 0.0, tau, w2=w2, line_frac=line_frac, trap=trap)


def case_equivalence(seed, nphot):
    print(f'\n1. alpha = 1 radial flow vs homologous flow (explicit trapping, tau = 3), {nphot} photons')
    out = []
    for g in (HomologousEllipsoid(VMAX, 0.3 * VMAX, 0.0, 3.0, trap='explicit'),
              RadialFlowEllipsoid(VMAX, 0.3 * VMAX, VMAX, 1.0, 0.0, 3.0, r_in=1e-9 * VMAX, trap='explicit')):
        specs, _ = run(g, nphot, [90], np.random.default_rng(seed))
        out.append(specs[0].integrated()[4])
    print(f'   homologous q_line = {out[0]:+.6f}   radial q_line = {out[1]:+.6f}   difference {out[1] - out[0]:+.1e}')


def case_scattering(rng, nphot, taus=(1.0, 3.0, 10.0)):
    print(f'\n2. Scattered continuum, i = 90, W2 = 1, {nphot} photons: q of the line-scattered light')
    print('   alpha   tau    explicit            single              effective          single/expl  eff/expl')
    rows = []
    for alpha in (0.25, 3.0):
        for tau in taus:
            res = {}
            for mode in ('explicit', 'single', 'effective'):
                specs, _ = run(geom(alpha, tau, mode), nphot, [90], rng)
                res[mode] = specs[0].integrated()[4:]
            (qe, se), (qs, ss), (qf, sf) = res['explicit'], res['single'], res['effective']
            rows.append((alpha, tau, qe, se, qs, ss, qf, sf))
            print(f'   {alpha:5.2f}  {tau:4.0f}   {qe:+.5f}±{se:.5f}  {qs:+.5f}±{ss:.5f}  {qf:+.5f}±{sf:.5f}'
                  f'   {qs / qe:6.2f}     {qf / qe:.3f}±{np.hypot(sf / qe, se * qf / qe ** 2):.3f}')
    return rows


def case_emission(rng, nphot, taus=(0.3, 1.0, 3.0, 10.0, 30.0)):
    print(f'\n3. Pure line emission (no continuum), i = 90, W2 = 1, explicit trapping: q vs tau')
    print('   (single and effective give exactly 0, since they emit line photons unpolarized)')
    print('   alpha   tau     q (%)           q tau (%)')
    rows = []
    for alpha in (0.25, 3.0):
        for tau in taus:
            n = int(nphot * min(1.0, 3.0 / tau))
            specs, _ = run(geom(alpha, tau, 'explicit', line_frac=1.0), n, [90], rng)
            q, u, sq, su, _, _ = specs[0].integrated()
            rows.append((alpha, tau, q, sq))
            print(f'   {alpha:5.2f}  {tau:5.1f}   {100 * q:+.3f}±{100 * sq:.3f}   {100 * q * tau:+7.2f}')
    return rows


def case_sphere(rng, nphot):
    print(f'\n4. Sphere, pure line emission, alpha = 0.25, tau = 3, explicit, {nphot} photons: must be 0')
    specs, _ = run(geom(0.25, 3.0, 'explicit', line_frac=1.0, c=1.0), nphot, [30, 90], rng)
    for s in specs:
        q, u, sq, su, _, _ = s.integrated()
        print(f'   i={s.angle:3d}  q = {q:+.5f} ± {sq:.5f} ({q / sq:+.1f} sig)   u = {u:+.5f} ± {su:.5f} ({u / su:+.1f} sig)')


def plot(results, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    r = np.array(results['scat'])
    for alpha, color in ((0.25, 'C0'), (3.0, 'C1')):
        k = r[:, 0] == alpha
        ax[0].errorbar(r[k, 1], r[k, 2], r[k, 3], fmt='o', color=color, label=f'explicit, alpha={alpha}')
        ax[0].plot(r[k, 1], r[k, 4], 's--', mfc='none', color=color, label=f'single, alpha={alpha}')
        ax[0].plot(r[k, 1], r[k, 6], 'x:', color=color, label=f'effective, alpha={alpha}')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Sobolev tau (at r = a)')
    ax[0].set_ylabel('q of line-scattered light')
    ax[0].set_title('Scattered continuum, anisotropic escape')
    ax[0].legend(fontsize=7)

    r = np.array(results['emit'])
    for alpha, color in ((0.25, 'C0'), (3.0, 'C1')):
        k = r[:, 0] == alpha
        ax[1].errorbar(r[k, 1], 100 * r[k, 2], 100 * r[k, 3], fmt='o-', color=color, label=f'explicit, alpha={alpha}')
    ax[1].axhline(0, color='k', lw=0.8, label='single, effective')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Sobolev tau (at r = a)')
    ax[1].set_ylabel('q (%)')
    ax[1].set_title('Pure line emission: polarization from trapping')
    ax[1].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    print(f'\nPlot written to {fname}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', type=float, default=1.0, help='scale factor for photon numbers')
    parser.add_argument('-seed', type=int, default=20261010)
    parser.add_argument('-plot', default=None, help='write a summary plot to this file')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    results = {}
    case_equivalence(args.seed, int(2e5 * args.n))
    results['scat'] = case_scattering(rng, int(6e5 * args.n))
    results['emit'] = case_emission(rng, int(4e5 * args.n))
    case_sphere(rng, int(4e5 * args.n))
    print(f'\nTotal time {time.time() - t0:.0f} s')
    if args.plot:
        plot(results, args.plot)


if __name__ == '__main__':
    main()
