#!/usr/bin/env python
'''
Compare q from the Sirocco runs made by make_bm.py with single scattering
(single_scatter.py), which is Brown & McLean (1977) with attenuation, over
4500-6500 A (away from the edges of the band, where the hot electrons
Doppler-shift scattered light in or out).  Prints a table, with u as a gauge
of the noise, and writes compare_bm.png:

    left:  q / tau_eq against inclination, for Sirocco (points; errors from
           |u|), single scattering (solid) and Brown & McLean (dotted)
    right: q (Sirocco) / q (single scattering) at 50-90 deg, the ratio that
           testdef.py checks, with the allowed range shaded

usage: compare_bm.py root:tau_eq [root:tau_eq ...]      e.g. compare_bm.py bm01:0.1 bm003:0.03
'''
import os
import sys
sys.path.insert(0, os.path.join(os.environ['SIROCCO'], 'py_progs'))
import numpy as np
import plot_pol
from single_scatter import predict, total_q

WMIN, WMAX = 4500, 6500
INCS = (10, 30, 50, 70, 90)


def compare(root, tau):
    '''
    For each inclination, Sirocco's q and |u| (%) over WMIN-WMAX, the single
    scattering and Brown & McLean predictions (%), and the ratio of q to single
    scattering.
    '''
    r = {round(x['angle']): x for x in plot_pol.band_pol(root + '.spec', WMIN, WMAX)}
    pred = predict(tau, INCS)
    rows = []
    for inc in INCS:
        Q, I_sc, bm = pred[inc]
        ss = total_q(tau, inc, Q, I_sc)
        rows.append(dict(inc=inc, q=r[inc]['q'], u=abs(r[inc]['u']), ss=ss, bm=100 * bm, ratio=r[inc]['q'] / ss))
    return rows


def doit(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for k, (root, tau) in enumerate(cases):
        rows = compare(root, tau)
        print('%s (tau_eq %g), %d-%d A' % (root, tau, WMIN, WMAX))
        for r in rows:
            print('   i=%3d  q = %+8.5f%%  |u| = %.5f%%  single scattering %+8.5f%%  B&M %+8.5f%%  ratio %.4f'
                  % (r['inc'], r['q'], r['u'], r['ss'], r['bm'], r['ratio']))
        c = 'C%d' % k
        incs = np.array([r['inc'] for r in rows])
        q, u, ss = (np.array([r[x] for r in rows]) for x in ('q', 'u', 'ss'))
        ii = np.linspace(0, 90, 91)
        ax[0].plot(ii, 100 * tau / 4 * 0.4 * np.sin(np.radians(ii)) ** 2 / tau, ':', color=c, lw=1)
        ax[0].plot(incs, ss / tau, '-', color=c, lw=1)
        ax[0].errorbar(incs, q / tau, yerr=u / tau, fmt='o', color=c, ms=5, label='tau_eq = %g' % tau)
        m = incs >= 50
        ax[1].errorbar(incs[m], (q / ss)[m], yerr=(u / ss)[m], fmt='o-', color=c, ms=5, label='tau_eq = %g' % tau)
        ax[1].axhspan(0.96, 1 + tau + 0.04, color=c, alpha=0.07)
    ax[0].set_xlabel('inclination (deg)')
    ax[0].set_ylabel('q / tau_eq (%)')
    ax[0].set_title('Sirocco (points; error from |u|), single scattering with attenuation (solid),\n'
                    'Brown & McLean (dotted); %d-%d A' % (WMIN, WMAX), fontsize=9)
    ax[0].legend(fontsize=8)
    ax[1].axhline(1, color='0.5', lw=0.5)
    ax[1].set_xlabel('inclination (deg)')
    ax[1].set_ylabel('q (Sirocco) / q (single scattering)')
    ax[1].set_title('Ratio checked by the test (shaded: allowed range for each tau_eq)', fontsize=9)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig('compare_bm.png', dpi=110)
    print('compare_bm.png')


def steer(argv):
    '''
    Process the command line
    '''

    words = []

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            words.append(argv[i])
        i += 1

    if len(words) == 0:
        print(__doc__)
        return

    doit([(a.split(':')[0], float(a.split(':')[1])) for a in words])


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
