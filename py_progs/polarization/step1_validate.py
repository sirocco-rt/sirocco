#!/usr/bin/env python
'''
Check the polarized Thomson-scattering routines against their analytic
distributions (Step 1), and plot the comparisons.


Command line usage (if any):

    usage::

        step1_validate.py [-h] [-n SCALE] [-seed SEED] [-plot FILE]

    -n      multiplies the default photon numbers (default 1; ~30 s)
    -seed   seed for the random number generator
    -plot   write a summary plot to FILE

Description:

    1. cos(Theta) drawn by sample_thomson_cos against (3/8)(1 + mu^2).
    2. The azimuth phi drawn by sample_thomson_phi for a polarized beam
       against p(phi | Theta) ~ 1 + a (q cos 2phi + u sin 2phi),
       a = sin^2/(1 + cos^2).
    3. A polarized beam is scattered with scatter_es and the photons leaving
       within a small cone around each of many observer directions are
       collected.  Their number, and their mean q and u on the sky, are
       compared with the peel-off prediction of peel_es.  For one incoming
       beam q and u are deterministic functions of direction, so their only
       difference from peel-off comes from averaging over the cone.

Primary routines:

    main, case_cos, case_phi, case_cones, plot

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
from polarization.stokes import normalize, dot, sky_north, to_sky
from polarization.es_scatter import sample_thomson_cos, sample_thomson_phi, scatter_es, peel_es


def case_cos(rng, n):
    mu = sample_thomson_cos(rng, n)
    edges = np.linspace(-1, 1, 41)
    h, _ = np.histogram(mu, edges, density=True)
    c = 0.5 * (edges[1:] + edges[:-1])
    # the exact average of (3/8)(1 + mu^2) over each bin
    cdf = 0.375 * (edges + edges ** 3 / 3.0)
    expected = np.diff(cdf) / np.diff(edges)
    err = np.sqrt(expected / (np.diff(edges) * n))
    chi2 = np.sum(((h - expected) / err) ** 2)
    print(f'\n1. cos(Theta): chi^2 = {chi2:.1f} for {len(c)} bins')
    return c, h, err, expected


def case_phi(rng, n, q=0.8, u=0.0):
    out = []
    print(f'\n2. phi for a beam with q = {q}, u = {u}:')
    for c in (0.0, 0.6):
        phi = sample_thomson_phi(rng, np.full(n, c), np.full(n, q), np.full(n, u))
        edges = np.linspace(0, 2 * np.pi, 37)
        h, _ = np.histogram(phi, edges, density=True)
        x = 0.5 * (edges[1:] + edges[:-1])
        a = (1 - c * c) / (1 + c * c)
        # the exact average of the analytic density over each bin
        cdf = (edges + 0.5 * a * (q * np.sin(2 * edges) - u * np.cos(2 * edges))) / (2 * np.pi)
        expected = np.diff(cdf) / np.diff(edges)
        err = np.sqrt(expected / (np.diff(edges) * n))
        chi2 = np.sum(((h - expected) / err) ** 2)
        print(f'   cos(Theta) = {c}: chi^2 = {chi2:.1f} for {len(x)} bins')
        out.append((c, x, h, err, expected))
    return out


def case_cones(rng, nphot, nobs=40, q=0.5, u=-0.4, cone_deg=3.0):
    n_in = normalize(np.array([0.3, -0.5, 0.8]))
    p = photon.make_photons(nphot)
    p['lmn'] = n_in
    photon.set_unpolarized(p)
    e1 = p['pol_ref'][0].copy()
    p['stokes'][:, 0] = q
    p['stokes'][:, 1] = u
    p['w'] = 1.0
    p0 = p[:1].copy()
    scatter_es(p, rng)

    # observer directions spread over the sphere
    angles = np.degrees(np.arccos(1 - 2 * (np.arange(nobs) + 0.5) / nobs))
    phases = (np.arange(nobs) * 0.618034) % 1.0
    cone = np.radians(cone_deg)
    rows = []
    for ang, ph in zip(angles, phases):
        psi = -2 * np.pi * ph
        i = np.radians(ang)
        n_obs = np.array([np.sin(i) * np.cos(psi), np.sin(i) * np.sin(psi), np.cos(i)])
        N = sky_north(ang, ph)
        f, qp, up = peel_es(p0, n_obs, N)
        sel = dot(p['lmn'], n_obs) > np.cos(cone)
        ns = sel.sum()
        frac = ns / nphot / ((1 - np.cos(cone)) / 2)
        qs, us = to_sky(p['stokes'][sel, 0], p['stokes'][sel, 1], p['pol_ref'][sel], p['lmn'][sel], N)
        rows.append((f[0], frac, np.sqrt(ns) / nphot / ((1 - np.cos(cone)) / 2),
                     qp[0], qs.mean(), qs.std() / np.sqrt(ns), up[0], us.mean(), us.std() / np.sqrt(ns)))
    r = np.array(rows)
    print(f'\n3. Sampled scatters vs peel-off, {nobs} observer directions, {nphot} photons, {cone_deg} deg cones')
    d = (r[:, 1] - r[:, 0]) / r[:, 2]
    print(f'   weight f: rms (sampled - peel)/sigma = {np.sqrt(np.mean(d * d)):.2f}, max |.| = {np.abs(d).max():.1f}')
    # For a single incoming beam, q and u are deterministic functions of the
    # outgoing direction, so the only difference is from averaging over the
    # finite cone; report it directly.
    for name, j in (('q', 3), ('u', 6)):
        print(f'   {name:8s}: max |sampled - peel| = {np.abs(r[:, j + 1] - r[:, j]).max():.4f}')
    return r


def plot(results, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 4, figsize=(18, 4.3))
    c, h, err, expected = results['cos']
    ax[0].errorbar(c, h, err, fmt='.', label='sample_thomson_cos')
    ax[0].plot(c, expected, 'k-', lw=1, label='(3/8)(1 + mu^2), bin averaged')
    ax[0].set_xlabel('cos Theta')
    ax[0].set_ylabel('probability density')
    ax[0].set_title('Scattering angle')
    ax[0].legend(fontsize=8)

    for (cc, x, hh, ee, ex), color in zip(results['phi'], ('C0', 'C1')):
        ax[1].errorbar(x, hh, ee, fmt='.', color=color, label=f'sampled, cos Theta = {cc}')
        ax[1].plot(x, ex, '-', color=color, lw=1)
    ax[1].set_xlabel('phi (rad)')
    ax[1].set_ylabel('probability density')
    ax[1].set_title('Azimuth for a beam with q = 0.8 (lines: analytic)')
    ax[1].legend(fontsize=8)

    r = results['cones']
    ax[2].errorbar(r[:, 0], r[:, 1], r[:, 2], fmt='o', ms=3)
    lim = [0.25, 1.6]
    ax[2].plot(lim, lim, 'k-', lw=0.8)
    ax[2].set_xlabel('peel-off weight 4 pi p')
    ax[2].set_ylabel('sampled (fraction in cone / solid angle)')
    ax[2].set_title('Direction distribution: sampled vs peel-off')

    ax[3].errorbar(r[:, 3], r[:, 4], r[:, 5], fmt='o', ms=3, label='q')
    ax[3].errorbar(r[:, 6], r[:, 7], r[:, 8], fmt='s', ms=3, mfc='none', label='u')
    lim = [-1, 1]
    ax[3].plot(lim, lim, 'k-', lw=0.8)
    ax[3].set_xlabel('peel-off prediction (sky frame)')
    ax[3].set_ylabel('sampled mean in cone')
    ax[3].set_title('Stokes parameters: sampled vs peel-off')
    ax[3].legend(fontsize=8)

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
    results = {'cos': case_cos(rng, int(2e6 * args.n)),
               'phi': case_phi(rng, int(2e6 * args.n)),
               'cones': case_cones(rng, int(4e6 * args.n))}
    print(f'\nTotal time {time.time() - t0:.0f} s')
    if args.plot:
        plot(results, args.plot)


if __name__ == '__main__':
    main()
