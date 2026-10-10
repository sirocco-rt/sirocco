#!/usr/bin/env python
'''
Validate the toy polarized transport against known results (Step 2).

Cases:

1. Uniform sphere around a point source: the net polarization must vanish
   at every inclination.
2. Optically thin oblate and prolate ellipsoids: compare with Brown &
   McLean (1977).  Oblate gives Q > 0, prolate Q < 0 (IAU convention).
3. A moderately thick oblate ellipsoid: no analytic answer, but peel-off
   and live-or-die must agree.
4. A thick slab lit from below: compare with the deterministic solution in
   milne.py, whose limb polarization is Chandrasekhar's 11.71%.

Usage:

    python step2_validate.py [-n SCALE] [-seed SEED] [-plot FILE]

SCALE multiplies the default photon numbers (default 1; a few minutes).
'''

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toy_transport import Ellipsoid, Slab, run, brown_mclean
from milne import solve_slab


def nsig(x, x0, sig):
    return (x - x0) / sig if sig > 0 else np.inf


def case_sphere(rng, nphot, kappa=2.0, angles=(10, 45, 70, 90)):
    geom = Ellipsoid(1.0, 1.0, kappa)
    specs, _ = run(geom, nphot, angles, rng)
    print(f'\n1. Uniform sphere, tau = {kappa}, {nphot} photons: net polarization must be 0')
    print('   angle        q            u        q/sig  u/sig')
    rows = []
    for s in specs:
        q, u, sq, su = s.result()
        rows.append((s.angle, q, u, sq, su))
        print(f'   {s.angle:5.1f}  {q:+.5f}±{sq:.5f}  {u:+.5f}±{su:.5f}  {q / sq:+5.1f}  {u / su:+5.1f}')
    return rows


def case_brown_mclean(rng, nphot, kappa=0.05, angles=(10, 30, 50, 70, 90)):
    print(f'\n2. Optically thin ellipsoids (kappa = {kappa}), {nphot} photons, vs Brown & McLean (1977)')
    out = {}
    for label, a, c in [('oblate  a=1, c=0.3', 1.0, 0.3), ('prolate a=0.3, c=1', 0.3, 1.0)]:
        geom = Ellipsoid(a, c, kappa)
        specs, _ = run(geom, nphot, angles, rng)
        _, tau0, gamma = brown_mclean(geom, 90)
        print(f'   {label}: tau0 = {tau0:.4f}, gamma = {gamma:.3f}')
        print('   angle      q (MC)            q (B&M)    ratio    u/sig')
        rows = []
        for s in specs:
            q, u, sq, su = s.result()
            qbm = brown_mclean(geom, s.angle)[0]
            rows.append((s.angle, q, sq, qbm, u, su))
            ratio = q / qbm if abs(qbm) > 3 * sq else np.nan
            print(f'   {s.angle:5.1f}  {q:+.6f}±{sq:.6f}  {qbm:+.6f}  {ratio:6.3f}  {u / su:+5.1f}')
        out[label] = rows
    return out


def case_live_vs_peel(rng, nphot, a=1.0, c=0.3, kappa=1.0, nbins=10):
    geom = Ellipsoid(a, c, kappa)
    mu_c = (np.arange(nbins) + 0.5) / nbins
    angles = np.degrees(np.arccos(mu_c))
    specs, live = run(geom, nphot, angles, rng, live_bins=2 * nbins)
    mu, I_l, q_l, u_l, sq_l, su_l = live.result()
    # Fold the two hemispheres (the ellipsoid is symmetric about z = 0).
    hi = mu > 0
    print(f'\n3. Oblate ellipsoid a={a}, c={c}, kappa={kappa}, {nphot} photons: live-or-die vs peel-off')
    print('     mu      I(peel)  I(live)     q (peel)          q (live)        diff/sig')
    rows = []
    for k, s in zip(np.where(hi)[0], specs):
        kk = len(mu) - 1 - k                       # mirror bin, -mu
        Il = 0.5 * (I_l[k] + I_l[kk])
        ql = 0.5 * (q_l[k] + q_l[kk])
        sql = 0.5 * np.hypot(sq_l[k], sq_l[kk])
        q, u, sq, su = s.result()
        d = nsig(ql, q, np.hypot(sq, sql))
        rows.append((mu[k], s.I, Il, q, sq, ql, sql))
        print(f'   {mu[k]:.3f}  {s.I:.4f}  {Il:.4f}  {q:+.5f}±{sq:.5f}  {ql:+.5f}±{sql:.5f}  {d:+5.1f}')
    return rows


def case_slab(rng, nphot, tau_max=10.0, angles=(0, 30, 60, 75, 84, 87, 89)):
    geom = Slab(tau_max)
    specs, live = run(geom, nphot, angles, rng, live_bins=40)
    emergent = solve_slab(tau_max)
    print(f'\n4. Slab, tau_max = {tau_max}, {nphot} photons, vs deterministic solution')
    print('   angle    mu      I/mu (MC)         I (exact)   P% (MC)        P% (exact)   q diff/sig')
    rows = []
    for s in specs:
        m = np.cos(np.radians(s.angle))
        I0, q0 = emergent(m)
        q, u, sq, su = s.result()
        sI = s.errors()[0]
        rows.append((m, s.I / m, sI / m, I0[0], q, sq, q0[0]))
        print(f'   {s.angle:5.1f}  {m:.4f}  {s.I / m:.4f}±{sI / m:.4f}  {I0[0]:.4f}   '
              f'{-100 * q:6.2f}±{100 * sq:.2f}   {-100 * q0[0]:6.2f}   {nsig(q, q0[0], sq):+5.1f}')
    I_lim, q_lim = emergent(1e-4)
    print(f'   exact limb polarization (mu -> 0): {-100 * q_lim[0]:.2f}%  (Chandrasekhar: 11.71%)')
    return rows, live, emergent


def plot(results, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

    for label, rows in results['bm'].items():
        r = np.array(rows)
        ax[0].errorbar(r[:, 0], 100 * r[:, 1], 100 * r[:, 2], fmt='o', label=f'MC {label}')
        ax[0].plot(r[:, 0], 100 * r[:, 3], 'k:', lw=1)
    ax[0].axhline(0, color='0.7', lw=0.5)
    ax[0].set_xlabel('inclination (deg)')
    ax[0].set_ylabel('q = Q/I (%)')
    ax[0].set_title('Thin ellipsoids vs Brown & McLean (dotted)')
    ax[0].legend(fontsize=8)

    r = np.array(results['live'])
    ax[1].errorbar(r[:, 0], 100 * r[:, 3], 100 * r[:, 4], fmt='o', label='peel-off')
    ax[1].errorbar(r[:, 0] + 0.01, 100 * r[:, 5], 100 * r[:, 6], fmt='s', mfc='none', label='live-or-die')
    ax[1].set_xlabel('mu = cos i')
    ax[1].set_ylabel('q = Q/I (%)')
    ax[1].set_title('Oblate ellipsoid, kappa = 1')
    ax[1].legend(fontsize=8)

    rows, live, emergent = results['slab']
    r = np.array(rows)
    mm = np.linspace(1e-3, 1, 300)
    ax[2].plot(mm, -100 * emergent(mm)[1], 'k-', lw=1, label='deterministic')
    ax[2].errorbar(r[:, 0], -100 * r[:, 4], 100 * r[:, 5], fmt='o', label='MC peel-off')
    mu, I, q, u, sq, su = live.result()
    ok = mu > 0
    ax[2].errorbar(mu[ok], -100 * q[ok], 100 * sq[ok], fmt='s', mfc='none', ms=4, label='MC live-or-die')
    ax[2].set_xlabel('mu = cos i')
    ax[2].set_ylabel('P = -q (%)')
    ax[2].set_title('Slab (Milne problem)')
    ax[2].legend(fontsize=8)

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
    results['sphere'] = case_sphere(rng, int(1e6 * args.n))
    results['bm'] = case_brown_mclean(rng, int(4e6 * args.n))
    results['live'] = case_live_vs_peel(rng, int(1e6 * args.n))
    results['slab'] = case_slab(rng, int(1e6 * args.n))
    print(f'\nTotal time {time.time() - t0:.0f} s')
    if args.plot:
        plot(results, args.plot)


if __name__ == '__main__':
    main()
