#!/usr/bin/env python
'''
Validate polarized resonance-line scattering in homologous flow (Step 3).


Command line usage (if any):

    usage::

        step3_validate.py [-h] [-n SCALE] [-seed SEED] [-plot FILE]

    -n      multiplies the default photon numbers (default 1; ~3 minutes)
    -seed   seed for the random number generator
    -plot   write a summary plot to FILE

Description:

    Runs five cases in a homologously expanding ellipsoid (v_max = 0.001 c):

    1. The polarizability W2(J_l, J_u) against the textbook values, and the
       fraction (7/10) W2 of the polarization kept per local line scatter.
    2. An optically thin line (pure line scattering of a central continuum):
       the band-integrated polarization must equal W2 times the Brown &
       McLean result for the equivalent electron-scattering opacity.
    3. Trapping in thick lines: q of the line-scattered light for the
       explicit (exact Sobolev), single (Sirocco's current single scatter)
       and effective (single scatter with W2_eff = W2 F) treatments, against
       the analytic factor F = beta / (1 - (1 - beta)(7/10) W2).
    4. A uniform sphere with electron scattering and a thick line: the net
       polarization must vanish.
    5. Peel-off against live-or-die, with a thick line and electron
       scattering.

    It then computes line-effect spectra (an electron-scattering envelope
    with unpolarized line emission and a thick resonance line) for the three
    trapping treatments.

Primary routines:

    main, case_w2, case_thin, case_trapping, case_sphere, case_live,
    case_line_effect, plot

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

from polarization import photon
from polarization.stokes import observer_direction, sky_north
from polarization.line_scatter import w2_factor, scatter_line, peel_line
from polarization.sobolev_flow import HomologousEllipsoid, run, thin_line_brown_mclean, trap_factor, p_escape

VMAX = 0.001


def case_w2():
    print('\n1. Polarizability W2(J_l -> J_u)')
    table = [((0, 1), 1.0), ((0.5, 0.5), 0.0), ((0.5, 1.5), 0.5), ((1, 1), 0.25), ((1, 2), 0.35),
             ((1.5, 1.5), 0.32), ((1.5, 2.5), 0.28), ((2, 1), 0.01), ((2, 2), 0.35), ((2, 3), 0.24)]
    worst = 0.0
    for (jl, ju), expected in table:
        w = w2_factor(jl, ju)
        worst = max(worst, abs(w - expected))
        print(f'   {jl:>3} -> {ju:<3}  W2 = {w:.4f}  (expected {expected})')
    print(f'   largest difference {worst:.1e}')


def case_retention(rng, nphot, ngen=6):
    '''
    Follow a beam through repeated line scatters at one point and measure
    the Q it would send to an edge-on observer after each scatter: the
    ratio between successive generations should be (7/10) W2.
    '''
    print(f'\n1b. Polarization kept per local scatter, {nphot} photons: Q(k+1)/Q(k) vs (7/10) W2')
    n_obs, N = observer_direction(90, 0), sky_north(90, 0)
    out = {}
    for w2 in (1.0, 0.5, 0.25):
        p = photon.make_photons(nphot)
        p['lmn'] = [0.0, 0.0, 1.0]
        photon.set_unpolarized(p)
        p['w'] = 1.0
        Q, sQ = [], []
        for _ in range(ngen):
            f, q, u = peel_line(p, n_obs, N, w2)
            Q.append(np.mean(f * q))
            sQ.append(np.std(f * q) / np.sqrt(nphot))
            scatter_line(p, rng, w2)
        Q, sQ = np.array(Q), np.array(sQ)
        out[w2] = (Q, sQ)
        print(f'   W2 = {w2}: ratios ' + ' '.join(f'{r:.3f}' for r in Q[1:4] / Q[:3]) + f'   (7/10) W2 = {0.7 * w2:.3f}')
    return out


def case_thin(rng, nphot):
    print(f'\n2. Optically thin line (tau = 0.02), {nphot} photons: q vs W2 x Brown & McLean equivalent')
    print('   (expected ratio ~0.996: O(tau) corrections from 1 - exp(-tau) and line absorption of the continuum)')
    rows = []
    for a, c in [(1.0, 0.3), (0.3, 1.0)]:
        for w2 in (1.0, 0.5):
            g = HomologousEllipsoid(a * VMAX, c * VMAX, 0.0, 0.02, w2=w2, trap='single')
            specs, _ = run(g, nphot, [90], rng)
            q, u, sq, su, _, _ = specs[0].integrated()
            pred = thin_line_brown_mclean(g, 90)
            rows.append((a, c, w2, q, sq, pred))
            print(f'   a={a} c={c} W2={w2}:  q = {q:+.4e} ± {sq:.1e}   predicted {pred:+.4e}   ratio {q / pred:.3f} ± {abs(sq / pred):.3f}')
    return rows


def case_trapping(rng, nphot, taus=(0.1, 0.3, 1, 3, 10, 30)):
    print(f'\n3. Trapping in the resonance region (oblate, pure line scattering, i = 90), {nphot} photons')
    print('   q of the line-scattered light; F = beta / (1 - (1 - beta)(7/10) W2)')
    out = {}
    for w2 in (1.0, 0.5):
        print(f'   W2 = {w2}')
        print('     tau    beta     F      explicit            single     effective           expl/single  (eff-expl)/sig')
        rows = []
        for tau in taus:
            res = {}
            for mode in ('explicit', 'single', 'effective'):
                g = HomologousEllipsoid(VMAX, 0.3 * VMAX, 0.0, tau, w2=w2, trap=mode)
                specs, _ = run(g, nphot, [90], rng)
                res[mode] = specs[0].integrated()[4:]
            beta = float(p_escape(tau))
            F = trap_factor(beta, w2)
            (qe, se), (qs, ss), (qf, sf) = res['explicit'], res['single'], res['effective']
            rows.append((tau, beta, F, qe, se, qs, ss, qf, sf))
            print(f'   {tau:5.1f}  {beta:.3f}  {F:.3f}  {qe:+.5f}±{se:.5f}  {qs:+.5f}  {qf:+.5f}±{sf:.5f}'
                  f'   {qe / qs:.3f}±{se / qs:.3f}   {(qf - qe) / np.hypot(se, sf):+5.1f}')
        out[w2] = rows
    return out


def case_sphere(rng, nphot):
    g = HomologousEllipsoid(VMAX, VMAX, 500.0, 3.0, w2=1.0, trap='explicit')
    specs, _ = run(g, nphot, [20, 60, 90], rng)
    print(f'\n4. Uniform sphere, tau_es = 0.5, tau_line = 3, {nphot} photons: net polarization must be 0')
    for s in specs:
        q, u, sq, su, _, _ = s.integrated()
        print(f'   i={s.angle:3d}  q = {q:+.5f} ± {sq:.5f} ({q / sq:+.1f} sig)   u = {u:+.5f} ± {su:.5f} ({u / su:+.1f} sig)')


def case_live(rng, nphot, nbins=5):
    g = HomologousEllipsoid(VMAX, 0.3 * VMAX, 300.0, 3.0, w2=1.0, trap='explicit')
    mu_c = (np.arange(nbins) + 0.5) / nbins
    specs, live = run(g, nphot, np.degrees(np.arccos(mu_c)), rng, live_bins=2 * nbins)
    mu, I, q, u, sq, su = live.result()
    print(f'\n5. Oblate, tau_es = 0.3, tau_line = 3, explicit trapping, {nphot} photons: peel-off vs live-or-die')
    print('     mu    q (peel)            q (live)          diff/sig')
    rows = []
    for k, s in zip(np.where(mu > 0)[0], specs):
        kk = len(mu) - 1 - k
        qp, up, sqp, sup, _, _ = s.integrated()
        ql = 0.5 * (q[k] + q[kk])
        sql = 0.5 * np.hypot(sq[k], sq[kk])
        rows.append((mu[k], qp, sqp, ql, sql))
        print(f'   {mu[k]:.2f}  {qp:+.5f}±{sqp:.5f}  {ql:+.5f}±{sql:.5f}  {(ql - qp) / np.hypot(sqp, sql):+5.1f}')
    return rows


def case_line_effect(rng, nphot):
    print(f'\n6. Line effect: oblate, tau_es = 0.3, line tau = 3, W2 = 0.5, 30% of the light in unpolarized line emission,')
    print(f'   i = 90, {nphot} photons')
    out = {}
    for mode in ('explicit', 'single', 'effective'):
        g = HomologousEllipsoid(VMAX, 0.3 * VMAX, 300.0, 3.0, w2=0.5, line_frac=0.3, trap=mode)
        specs, _ = run(g, nphot, [90], rng, nbins=40)
        y, I, q, u, sq, su = specs[0].binned()
        out[mode] = (y, I, q, sq)
        cont = np.abs(y) > 1.2
        core = np.abs(y) < 0.3
        print(f'   {mode:9s}: q(continuum, |y|>1.2) = {np.average(q[cont], weights=I[cont]):+.5f}'
              f'   q(line core, |y|<0.3) = {np.average(q[core], weights=I[core]):+.5f}')
    return out


def plot(results, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    ax = axes[0]

    tt = np.geomspace(0.05, 50, 200)
    for (w2, rows), color in zip(results['trap'].items(), ('C0', 'C1')):
        r = np.array(rows)
        ax[0].plot(tt, trap_factor(p_escape(tt), w2), '-', lw=1, color=color, label=f'F, W2={w2}')
        ax[0].errorbar(r[:, 0], r[:, 3] / r[:, 5], r[:, 4] / r[:, 5], fmt='o', color=color,
                       label=f'explicit/single, W2={w2}')
        ax[0].plot(r[:, 0], r[:, 7] / r[:, 5], 'x', color='k')
    ax[0].plot(tt, p_escape(tt), 'k:', lw=1, label='beta')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Sobolev tau')
    ax[0].set_ylabel('q(line) / q(single scatter)')
    ax[0].set_title('Depolarization by trapping (x = effective)')
    ax[0].legend(fontsize=7)

    le = results['line_effect']
    for mode, (y, I, q, sq) in le.items():
        ax[1].plot(y, I, drawstyle='steps-mid', label=mode)
        ax[2].errorbar(y, 100 * q, 100 * sq, drawstyle='steps-mid', label=mode, lw=1)
    ax[1].set_xlabel('(nu - nu0) / (nu0 vmax/c)')
    ax[1].set_ylabel('I')
    ax[1].set_title('Line effect: intensity, i = 90')
    ax[1].legend(fontsize=8)
    ax[2].set_xlabel('(nu - nu0) / (nu0 vmax/c)')
    ax[2].set_ylabel('q = Q/I (%)')
    ax[2].set_title('Line effect: polarization, i = 90')
    ax[2].legend(fontsize=8)

    for (w2, (Q, sQ)), color in zip(results['retention'].items(), ('C0', 'C1', 'C2')):
        k = np.arange(len(Q))
        ax2 = axes[1]
        ax2[0].errorbar(k, np.abs(Q), sQ, fmt='o', color=color, label=f'W2 = {w2}')
        ax2[0].plot(k, np.abs(Q[0]) * (0.7 * w2) ** k, '-', color=color, lw=1)
    axes[1][0].set_yscale('log')
    axes[1][0].set_xlabel('number of local scatters k')
    axes[1][0].set_ylabel('|Q| toward an edge-on observer')
    axes[1][0].set_title('Polarization kept per scatter (lines: (0.7 W2)^k)')
    axes[1][0].legend(fontsize=8)

    r = np.array(results['thin'])
    labels = [f'{"oblate" if a > c else "prolate"}\nW2={w2}' for a, c, w2, *_ in results['thin']]
    x = np.arange(len(r))
    axes[1][1].errorbar(x, r[:, 3] / r[:, 5], np.abs(r[:, 4] / r[:, 5]), fmt='o')
    axes[1][1].axhline(1.0, color='k', lw=0.8)
    axes[1][1].axhline(0.996, color='k', lw=0.8, ls=':')
    axes[1][1].set_xticks(x)
    axes[1][1].set_xticklabels(labels, fontsize=8)
    axes[1][1].set_ylabel('q / (W2 x Brown & McLean)')
    axes[1][1].set_title('Optically thin line (dotted: expected O(tau) value)')
    axes[1][1].set_ylim(0.85, 1.15)

    axes[1][2].axis('off')

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
    case_w2()
    results['retention'] = case_retention(rng, int(2e6 * args.n))
    results['thin'] = case_thin(rng, int(1.2e7 * args.n))
    results['trap'] = case_trapping(rng, int(1e6 * args.n))
    case_sphere(rng, int(1e6 * args.n))
    results['live'] = case_live(rng, int(1e6 * args.n))
    results['line_effect'] = case_line_effect(rng, int(1e6 * args.n))
    print(f'\nTotal time {time.time() - t0:.0f} s')
    if args.plot:
        plot(results, args.plot)


if __name__ == '__main__':
    main()
