#!/usr/bin/env python
'''
Band-integrated polarization for the Hillier (1994) tests made by make_hillier.py.

usage:
    hillier_pol.py root [root ...]          q (%) at each observer angle
    hillier_pol.py -albedo root             q (%) vs albedo from the nscat-selected
                                            observers: I(a) = sum_j a^j I_j (Hillier App. A)

Band 4500-6500 A (away from the edges, where Doppler broadening by the 1e6 K
electrons moves flux in or out).  Errors: half the difference between the two
halves of the band is NOT independent, so for the angle tables the error is
taken from |u|, which should be 0 for an axisymmetric model.
'''
import os
import sys
sys.path.insert(0, os.path.join(os.environ['SIROCCO'], 'py_progs'))
import numpy as np
import plot_pol as cs

WMIN, WMAX = 4500, 6500


def band(fname):
    units, freq, lam, cols = cs.read_spec(fname)
    x = cs.lnu_dnu(units, freq, False)
    m = (lam > WMIN) & (lam < WMAX)
    out = []
    for o in cs.observers(cols):
        out.append((o, (cols[o] * x)[m].sum(), (cols[o + '_Q'] * x)[m].sum(), (cols[o + '_U'] * x)[m].sum()))
    return out


def angles(roots):
    for root in roots:
        print(root)
        for o, I, Q, U in band(root + '.spec'):
            print(f'  {o:10s} q = {100*Q/I:+7.3f}%  u = {100*U/I:+7.3f}%')


def albedo(root):
    b = band(root + '.spec')
    I = np.array([r[1] for r in b])
    Q = np.array([r[2] for r in b])
    U = np.array([r[3] for r in b])
    nj = len(I) - 1                       # last observer: >= nj scatters
    print(f'{root}: {b[0][0]}; fraction of I in each j (a = 1):')
    print('  ' + ' '.join(f'{v:.3f}' for v in I / I.sum()))
    # tail beyond nj: geometric continuation with the ratio of the last few exact orders
    rI = np.mean(I[nj - 3:nj] / I[nj - 4:nj - 1])
    print(f'  I_j ratio at j = {nj-3}..{nj-1}: {rI:.3f}')
    print('  albedo    q (%)    u (%)   [tail treated as geometric]')
    res = []
    for a in np.round(np.arange(0, 1.0001, 0.05), 3):
        w = a ** np.arange(nj)
        Ia = (w * I[:nj]).sum()
        Qa = (w * Q[:nj]).sum()
        Ua = (w * U[:nj]).sum()
        # >= nj bin: at a = 1 it is exact; otherwise weight each order by a^j with I_j ~ I_nj r^(j-nj)
        if a > 0:
            s = (a ** nj) * (1 - rI) / (1 - a * rI)     # sum_k a^(nj+k) r^k (1-r), normalised to the >=nj total
            Ia += s * I[nj]
            Qa += s * Q[nj]
            Ua += s * U[nj]
        res.append((a, 100 * Qa / Ia, 100 * Ua / Ia))
        print(f'  {a:5.2f}  {100*Qa/Ia:+7.3f}  {100*Ua/Ia:+7.3f}')
    return res


def steer(argv):
    '''
    Process the command line
    '''

    words = []
    xalbedo = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-albedo':
            xalbedo = True
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            words.append(argv[i])
        i += 1

    if len(words) == 0:
        print(__doc__)
        return

    if xalbedo:
        for r in words:
            albedo(r)
    else:
        angles(words)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
