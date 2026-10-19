#!/usr/bin/env python
'''
Sirocco against single scattering for pbl_es and eq_es: q of the scattered light (top: photons
scattered exactly once, solid, and at least once, dotted) and of the total light (bottom) against
wavelength.  Also prints the averages over the windows between the lines.  The predictions are read
from root.pred.txt, which geom_pol.py root theta_min theta_max -central writes.
Writes es_vs_prediction.png

usage: es_vs_prediction.py [root ...]        (default: pbl_es eq_es)
'''
import os
import sys
sys.path.insert(0, os.path.join(os.environ['SIROCCO'], 'py_progs'))
import plot_pol as p
import numpy as np

INCS = (10, 45, 80)
WINDOWS = ((1700, 1850), (2000, 2700))     # between the resonance lines

title = {'pbl_es': 'pbl_es: polar cone 0-20$^\\circ$ (E parallel to disk: q < 0)',
         'eq_es': 'eq_es: equatorial wedge 60-88$^\\circ$ (E along axis: q > 0)'}


def read_pred(root):
    '''{inc: (q of the light scattered once, q of the total light)} in % from root.pred.txt'''
    pred = {}
    for line in open(root + '.pred.txt'):
        words = line.split()
        pred[int(words[0])] = (float(words[4].rstrip('%')), float(words[6].rstrip('%')))
    return pred


def compare(root):
    '''
    For each inclination, q (%) in the windows between the lines of the light
    scattered once, of the light scattered at least once and of the total
    light, the predictions, and the ratio of light scattered once to light
    scattered more than once.
    '''
    pred = read_pred(root)
    u, f, w, c = p.read_spec(root + '.spec')
    fac = p.lnu_dnu(u, f, False)
    m = np.zeros(len(w), dtype=bool)
    for lo, hi in WINDOWS:
        m |= (w > lo) & (w < hi)
    rows = []
    for inc in INCS:
        a = 'A%dP0.50' % inc
        s1, sm = a + '_sc:1', a + '_sc:>1'
        q1, qm, qt = [100 * (c[x + '_Q'] * fac)[m].sum() / (c[x] * fac)[m].sum() for x in (s1, sm, a)]
        f1 = (c[s1] * fac)[m].sum() / (c[sm] * fac)[m].sum()
        rows.append(dict(inc=inc, q1=q1, qm=qm, qtot=qt, pred1=pred[inc][0], predtot=pred[inc][1], f1=f1))
    return rows


def main(roots):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, len(roots), figsize=(6.5 * len(roots), 8), sharex=True, squeeze=False)
    for j, r in enumerate(roots):
        for row in compare(r):
            print('%-6s %2d  q(1 scatter) %6.1f%%  q(>=1) %6.1f%%  pred %6.1f%% | q_tot %6.2f%%  pred %6.2f%% | '
                  'single/all scattered %.2f' % (r, row['inc'], row['q1'], row['qm'], row['pred1'], row['qtot'],
                                                 row['predtot'], row['f1']))
        pred = read_pred(r)
        u, f, w, c = p.read_spec(r + '.spec')
        fac = p.lnu_dnu(u, f, False)
        nb = 50
        n = (len(w) // nb) * nb
        cb = lambda x: (x * fac)[:n].reshape(-1, nb).sum(1)
        wb = w[:n].reshape(-1, nb).mean(1)
        for k, inc in enumerate(INCS):
            col = 'C%d' % k
            a = 'A%dP0.50' % inc
            s1, sm = a + '_sc:1', a + '_sc:>1'
            ax[0, j].step(wb, 100 * cb(c[s1 + '_Q']) / cb(c[s1]), where='mid', color=col, lw=0.9, label='%d$^\\circ$' % inc)
            ax[0, j].step(wb, 100 * cb(c[sm + '_Q']) / cb(c[sm]), where='mid', color=col, lw=0.8, ls=':', alpha=0.8)
            ax[0, j].axhline(pred[inc][0], color=col, ls='--', lw=1.3)
            ax[1, j].step(wb, 100 * cb(c[a + '_Q']) / cb(c[a]), where='mid', color=col, lw=0.8)
            ax[1, j].axhline(pred[inc][1], color=col, ls='--', lw=1.3)
        ax[0, j].set_title(title.get(r, r), fontsize=10)
        ax[0, j].set_ylabel('q of the scattered light (%)')
        ax[1, j].set_ylabel('q of the total light (%)')
        ax[1, j].set_xlabel('Wavelength ($\\AA$)')
        for a_ in ax[:, j]:
            a_.axhline(0, color='0.5', lw=0.5)
            for lo, hi in WINDOWS:
                a_.axvspan(lo, hi, color='0.92', zorder=0)
    ax[0, 0].legend(title='Sirocco: scattered once (solid), at least once (dotted);\nsingle-scattering prediction, attenuated (dashed)', fontsize=8, title_fontsize=8)
    fig.suptitle('Central source L = 1e37 erg/s, T = 1e5 K; constant density, tau_es = 1; simple atoms, no wind emission '
                 '(shaded: windows between lines)', fontsize=10)
    fig.tight_layout()
    fig.savefig('es_vs_prediction.png', dpi=100)
    print('es_vs_prediction.png')


if __name__ == '__main__':
    main(sys.argv[1:] if len(sys.argv) > 1 else ['pbl_es', 'eq_es'])
