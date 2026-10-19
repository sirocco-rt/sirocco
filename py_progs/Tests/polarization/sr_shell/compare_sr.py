#!/usr/bin/env python
'''
Compare the bolometric single-scattered I and Q from Sirocco with
sr_single.py, relative to the (attenuated) direct light at the same angle.
-nonrel runs (root ending in n) are compared with the static prediction.
Writes compare_sr.txt and compare_sr.png.

usage: compare_sr.py [root ...]        (default: sb0 sb0.1 sb0.3 sb0.3n ss0.3)
'''
import os
import sys
sys.path.insert(0, os.path.join(os.environ['SIROCCO'], 'py_progs'))
import numpy as np
import plot_pol as cs
from sr_single import predict, ANGLES

RUNS = ["sb0", "sb0.1", "sb0.3", "sb0.3n", "ss0.3"]


def measure(root):
    units, freq, lam, cols = cs.read_spec(root + '.spec')
    x = cs.lnu_dnu(units, freq, False)
    out = {}
    for inc in ANGLES:
        name = f'A{inc:02.0f}P0.50'
        d = {}
        for s in (0, 1):
            c = f'{name}_sc:{s}'
            d[s] = [(cols[c + t] * x).sum() for t in ('', '_Q', '_U')]
        # fraction of the scattered light in the outer 2% of bins at each end: is the band wide enough?
        sc = cols[f'{name}_sc:1'] * x
        nb = len(sc) // 50
        edge = (sc[:nb].sum() + sc[-nb:].sum()) / sc.sum()
        out[inc] = (d[0][0], d[1][0], d[1][1], d[1][2], edge)
    return out


def compare(root):
    '''
    Sirocco against the prediction for one run.  Returns beta, whether it is a
    -nonrel run, and for each angle a dict with the measured and predicted
    I_sc/I_0, Q_sc/I_0 and q_sc (%), U_sc/I_0, and the ratios (ratio_q is nan
    where the predicted Q is negligible).
    '''
    cells = np.load(root + '.cells.npz')
    beta = float(cells['beta'])
    nonrel = root.endswith('n')
    m = measure(root)
    rows = []
    for inc in ANGLES:
        # -nonrel: static scattering, and tau without gamma (1 - beta), which is what predict gives for beta = 0
        I, Q, U, tau = predict(cells, inc, 0.0 if nonrel else beta)
        I, Q = I * np.exp(tau), Q * np.exp(tau)
        I0, Is, Qs, Us, edge = m[inc]
        rows.append(dict(inc=inc, I=Is / I0, I_pred=I, Q=Qs / I0, Q_pred=Q, U=Us / I0,
                         q=100 * Qs / Is, q_pred=100 * Q / I if I > 0 else np.nan,
                         ratio_i=Is / I0 / I, ratio_q=(Qs / I0) / Q if abs(Q) > 1e-3 * I else np.nan))
    return beta, nonrel, rows


def main(runs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    lines = []
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    for root in runs:
        beta, nonrel, rows = compare(root)
        lines.append(f'{root}: beta = {beta}' + ('  (-nonrel: compared with the static prediction)' if nonrel else ''))
        lines.append('   i     I_sc/I_0 Sirocco  pred    ratio    Q_sc/I_0 (%) Sirocco  pred   ratio    U_sc/I_0 (%)  q_sc (%) Sir  pred')
        for r in rows:
            lines.append(f'  {r["inc"]:5.1f}  {r["I"]:.4e}  {r["I_pred"]:.4e}  {r["ratio_i"]:6.3f}    {100*r["Q"]:+8.4f}  {100*r["Q_pred"]:+8.4f}  '
                         + (f'{r["ratio_q"]:6.3f}' if np.isfinite(r['ratio_q']) else '   -- ') + f'    {100*r["U"]:+8.4f}    '
                         + f'{r["q"]:+7.2f}  ' + (f'{r["q_pred"]:+7.2f}' if np.isfinite(r['q_pred']) else ''))
        xs = [r['inc'] for r in rows]
        ax[0].plot(xs, [r['ratio_i'] for r in rows], 'o-', label=root)
        if root.startswith('sb'):
            ax[1].plot(xs, [r['ratio_q'] for r in rows], 'o-', label=root)
    ax[0].set_title('I_sc / I_0: Sirocco / exact single scattering')
    ax[1].set_title('Q_sc / I_0: Sirocco / exact single scattering')
    for a in ax:
        a.axhline(1, color='0.6', lw=0.5)
        a.set_xlabel('inclination (deg)')
        a.set_ylim(0.94, 1.06)
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig('compare_sr.png', dpi=110)
    open('compare_sr.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(sys.argv[1:] if len(sys.argv) > 1 else RUNS)
