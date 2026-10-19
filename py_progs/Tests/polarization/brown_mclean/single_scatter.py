'''
Exact single scattering (with attenuation in and out) of an isotropic point
source by n_e = n0 (r_in/r)^2 sin^2(theta), r_in < r < r_out, by quadrature.
Returns Q/I0 and I_sc/I0 for an observer at inclination i, where I0 = L/4pi,
and the q of the total light that this implies.  Prints a line for each of
10, 30, 50, 70 and 90 deg; testdef.py calls predict() and total_q().

usage: single_scatter.py tau_eq
'''
import sys
import numpy as np
SIG = 6.6524587e-25
r_in, r_out = 1e11, 1e12


def predict(tau_eq, incs=(10, 30, 50, 70, 90)):
    '''
    Single scattering for each inclination: {inc: (Q/I0, I_sc/I0, Brown & McLean Q/I0)}
    '''
    n0 = tau_eq / (SIG * r_in * (1 - r_in / r_out))

    def ne(x):
        r = np.linalg.norm(x, axis=-1)
        s2 = 1 - (x[..., 2] / np.maximum(r, 1e-30)) ** 2
        return np.where((r >= r_in) & (r <= r_out), n0 * r_in ** 2 / np.maximum(r, 1e-30) ** 2 * s2, 0.0)

    # quadrature in r (Gauss-Legendre in 1/r), mu, phi
    nr, nmu, nphi = 60, 64, 64
    y, wy = np.polynomial.legendre.leggauss(nr)          # in 1/r
    a, b = 1 / r_out, 1 / r_in
    inv = 0.5 * (b - a) * y + 0.5 * (b + a); r = 1 / inv
    wr = 0.5 * (b - a) * wy * r ** 2                       # dr = r^2 d(1/r)
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi; wphi = 2 * np.pi / nphi
    R, MU, PHI = np.meshgrid(r, mu, phi, indexing='ij')
    W = (wr[:, None, None] * wmu[None, :, None] * wphi) * R ** 2      # dV
    ST = np.sqrt(1 - MU ** 2)
    X = np.stack([R * ST * np.cos(PHI), R * ST * np.sin(PHI), R * MU], axis=-1)
    nin = X / R[..., None]
    n_here = ne(X)
    tau_in = SIG * n0 * r_in ** 2 * (1 - MU ** 2) * (1 / r_in - 1 / R)
    # path for tau_out: s from 0 to 2 r_out, Gauss-Legendre in segments
    s_nodes, s_w = np.polynomial.legendre.leggauss(400)
    out = {}
    for inc in incs:
        ii = np.radians(inc)
        o = np.array([np.sin(ii), 0, np.cos(ii)])
        north = np.array([-np.cos(ii), 0, np.sin(ii)])
        east = np.cross(north, o)
        tau_out = np.zeros(R.shape)
        smax = 2.2 * r_out
        s = 0.5 * smax * (s_nodes + 1); ws = 0.5 * smax * s_w
        for k in range(len(s)):
            tau_out += ws[k] * ne(X + s[k] * o)
        tau_out *= SIG
        cosx = nin @ o
        # E-vector of scattered light is along n_in x o (perpendicular to the scattering plane)
        e = np.cross(nin, o)
        en, ee = e @ north, e @ east
        cos2psi = (en ** 2 - ee ** 2) / np.maximum(en ** 2 + ee ** 2, 1e-300)
        f = n_here * SIG * 3 / (16 * np.pi) / R ** 2 * np.exp(-tau_in - tau_out) * W
        Q = np.sum(f * (1 - cosx ** 2) * cos2psi)
        I_sc = np.sum(f * (1 + cosx ** 2))
        bm = tau_eq / 4 * 0.4 * np.sin(ii) ** 2
        out[inc] = (Q, I_sc, bm)
    return out


def total_q(tau_eq, inc, Q, I_sc):
    '''
    q (%) of the total light: the attenuated direct light, e^(-tau_eq sin^2 i)
    times I0, plus the scattered light
    '''
    return 100 * Q / (np.exp(-tau_eq * np.sin(np.radians(inc)) ** 2) + I_sc)


def doit(tau_eq):
    '''
    Print the single-scattering prediction for an equatorial optical
    depth tau_eq
    '''

    for inc, (Q, I_sc, bm) in predict(tau_eq).items():
        print(f'tau_eq {tau_eq}: i={inc:2d}  single-scatter Q/I0 = {100*Q:+.6f}%  B&M {100*bm:+.6f}%  ratio {Q/bm:.4f}   '
              f'I_sc/I0 = {I_sc:.5f}  q_tot = {total_q(tau_eq, inc, Q, I_sc):+.6f}%')


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

    doit(float(words[0]))


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
