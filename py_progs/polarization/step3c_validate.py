#!/usr/bin/env python
'''
Test the capped explicit trapping chain proposed for the C version (Step 3c).


Command line usage (if any):

    usage::

        step3c_validate.py [-h] [-n SCALE] [-seed SEED] [-plot FILE]

    -n      multiplies the default photon numbers (default 1; ~15 minutes)
    -seed   seed for the random number generator
    -plot   write a summary plot to FILE

Description:

    The capped mode follows at most nmax polarized local scatters in a
    Sobolev resonance region, sums their peel-off contributions locally and
    attenuates them once per observer, and lets photons still trapped after
    nmax scatters escape unpolarized, as Sirocco's thermal trapping does now.

    1. With nmax effectively infinite, capped must reproduce explicit
       exactly (same random numbers), which checks the local summation of
       the peel-off contributions.
    2. For RadialFlowEllipsoid (oblate; alpha = 0.25 and 3) and tau = 3, 10
       and 30: q of the line-scattered continuum, and q of pure line
       emission, against nmax (fixed caps, and the automatic cap
       ceil(3 / beta_mean)), relative to the uncapped chain.
    3. The cost: local scatters plus fallback tries per line event, against
       nmax = 0, where every photon escapes as in Sirocco now (rejection with
       beta / beta_max, i.e. p_norm).

Primary routines:

    main, case_equivalence, case_cap, plot

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
NMAX_INF = 10 ** 9
NMAX_AUTO = -1
NMAXES = (0, 1, 3, 10, 30, NMAX_AUTO, NMAX_INF)
LABELS = {NMAX_AUTO: 'auto', NMAX_INF: 'inf'}


def geom(alpha, tau, nmax, line_frac, trap='capped'):
    nmax = 'auto' if nmax == NMAX_AUTO else nmax
    return RadialFlowEllipsoid(VMAX, 0.3 * VMAX, VMAX, alpha, 0.0, tau, line_frac=line_frac, trap=trap, nmax=nmax)


def case_equivalence(seed, nphot):
    print(f'\n1. capped with nmax = infinity vs explicit (same random numbers), {nphot} photons')
    worst = 0.0
    for alpha in (0.25, 1.0):
        for lf in (0.0, 0.5):
            out = []
            for trap in ('explicit', 'capped'):
                g = RadialFlowEllipsoid(VMAX, 0.3 * VMAX, VMAX, alpha, 200.0, 3.0, line_frac=lf, trap=trap,
                                        nmax=NMAX_INF, r_in=1e-9 * VMAX if alpha == 1.0 else None)
                specs, _ = run(g, nphot, [60, 90], np.random.default_rng(seed))
                out.append(np.array([[s.I.sum(), s.Q.sum(), s.U.sum()] for s in specs]))
            rel = np.max(np.abs(out[1] - out[0]) / np.maximum(np.abs(out[0]), 1e-30))
            worst = max(worst, rel)
            print(f'   alpha={alpha} line_frac={lf}: largest relative difference in I, Q, U = {rel:.1e}')
    print(f'   worst {worst:.1e}')
    return worst


def case_cap(rng, nphot, taus=(3.0, 10.0, 30.0)):
    print(f'\n2. q against nmax (i = 90), relative to the uncapped chain')
    print('   cost = (local scatters + fallback tries) per line event; nmax = 0 is Sirocco\'s present')
    print('   thermal trapping (rejection with beta / beta_max), so cost/cost(0) is the extra work')
    rows = []
    for kind, lf in (('scattered continuum', 0.0), ('pure line emission', 1.0)):
        for alpha in (0.25, 3.0):
            for tau in taus:
                n = int(nphot * min(1.0, 10.0 / tau))
                res = {}
                for nmax in NMAXES:
                    g = geom(alpha, tau, nmax, lf)
                    specs, _ = run(g, n, [90], rng)
                    q, u, sq, su, ql, sql = specs[0].integrated()
                    qq, ss = (ql, sql) if lf == 0.0 else (q, sq)
                    st = g.stats
                    cost = (st['local'] + st['tries']) / max(st['events'], 1)
                    res[nmax] = (qq, ss, cost, st['fallback'] / max(st['events'], 1))
                q_ref, s_ref = res[NMAX_INF][:2]
                cost0 = res[0][2]
                print(f'   {kind}, alpha={alpha}, tau={tau:g} ({n} photons): uncapped q = {q_ref:+.5f} ± {s_ref:.5f}')
                for nmax in NMAXES:
                    qq, ss, cost, ffb = res[nmax]
                    label = LABELS.get(nmax, str(nmax))
                    print(f'      nmax={label:>4}: q = {qq:+.5f} ± {ss:.5f}   q/q_inf = {qq / q_ref:6.3f}'
                          f'   fallback fraction {ffb:.3f}   cost/event {cost:6.1f} ({cost / cost0:5.1f} x Sirocco)')
                    rows.append((lf, alpha, tau, nmax, qq, ss, q_ref, s_ref, cost, ffb))
    return rows


def plot(rows, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    r = np.array(rows, dtype=float)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    xs = {n: k for k, n in enumerate(NMAXES)}
    labels = [LABELS.get(n, str(n)) for n in NMAXES]
    for panel, lf, title in ((0, 0.0, 'Scattered continuum'), (1, 1.0, 'Pure line emission')):
        for (alpha, tau), color in zip([(a, t) for a in (0.25, 3.0) for t in (3.0, 10.0, 30.0)],
                                       ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']):
            k = (r[:, 0] == lf) & (r[:, 1] == alpha) & (r[:, 2] == tau)
            x = np.array([xs[int(n)] for n in r[k, 3]]) + (0.05 * (color[-1] == '1'))
            ax[panel].errorbar(x, r[k, 4] / r[k, 6], r[k, 5] / np.abs(r[k, 6]), fmt='o-', color=color, ms=3,
                               label=f'alpha={alpha}, tau={tau:g}')
        ax[panel].axhline(1.0, color='k', lw=0.8)
        ax[panel].set_xticks(range(len(NMAXES)))
        ax[panel].set_xticklabels(labels)
        ax[panel].set_xlabel('nmax (local scatters followed)')
        ax[panel].set_ylabel('q / q(uncapped)')
        ax[panel].set_title(title)
        ax[panel].legend(fontsize=7)
    ax[0].set_ylim(-0.2, 1.4)
    ax[1].set_ylim(-0.4, 1.6)

    for (alpha, tau), color in zip([(a, t) for a in (0.25, 3.0) for t in (3.0, 10.0, 30.0)],
                                   ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']):
        k = (r[:, 0] == 0.0) & (r[:, 1] == alpha) & (r[:, 2] == tau)
        x = [xs[int(n)] for n in r[k, 3]]
        ax[2].plot(x, r[k, 8] / r[k, 8][0], 'o-', color=color, ms=3, label=f'alpha={alpha}, tau={tau:g}')
    ax[2].axhline(1.0, color='k', lw=0.8)
    ax[2].set_xticks(range(len(NMAXES)))
    ax[2].set_xticklabels(labels)
    ax[2].set_xlabel('nmax')
    ax[2].set_yscale('log')
    ax[2].set_ylabel('cost per line event / cost at nmax = 0')
    ax[2].set_title('Cost relative to Sirocco\'s present trapping')
    ax[2].legend(fontsize=7)

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
    case_equivalence(args.seed, int(5e4 * args.n))
    rows = case_cap(rng, int(6e5 * args.n))
    print(f'\nTotal time {time.time() - t0:.0f} s')
    if args.plot:
        plot(rows, args.plot)


if __name__ == '__main__':
    main()
