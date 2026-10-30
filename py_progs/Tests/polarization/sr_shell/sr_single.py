#!/usr/bin/env python
'''
Single scattering by an optically thin shell with a constant radial
bulk velocity, for a point source at the centre, with exact special relativity.
Independent of Sirocco; uses the cell densities written by make_srshell.py.
These are lab-frame densities n: Sirocco takes the comoving density as
n' = n / gamma (define_wind.c), and so does this calculation.

For a steady flow, the lab-frame energy scattered per unit time into dOmega is

    dE/dt dOmega = sum_cells n' sigma dr  int dOmega_r  (L / 4 pi)
                   gamma^2 (1 - beta)^2  D^3  (3 / 16 pi) (1 + mu'^2)

where gamma (1 - beta) appears twice, once in the scattering rate
(dtau = n' sigma gamma (1 - beta.n) ds) and once for the photon energy in the
comoving frame, D = 1 / (gamma (1 - beta mu)) with mu = r_hat.k, and mu' =
(mu - beta) / (1 - beta mu) is the comoving scattering angle (the incident light
is radial, along the velocity, so it is not aberrated).  The polarization is
P = (1 - mu'^2) / (1 + mu'^2), with E along r_hat x k: in the comoving frame E is
perpendicular to the (r_hat, k') plane, which contains the velocity, and a boost
along a vector in that plane leaves E perpendicular to it.

The direct light is attenuated by tau = sigma sum n' gamma (1 - beta) dr along
the line of sight, as Sirocco's peel-off is.  Everything is relative to the
attenuated direct light, which is what the nscat = 0 observer measures.

usage: sr_single.py root [root ...]     (prints q and I_sc/I_0 at 22.5 ... 90 deg)
'''
import sys
import numpy as np

SIGMA_T = 6.6524587e-25
ANGLES = [22.5, 45.0, 67.5, 90.0]


def predict(cells, inc, beta, nmu=8, nphi=256):
    re, te, ne = cells['re'], cells['te'], cells['ne']
    dr = np.diff(re)
    gam = 1 / np.sqrt(1 - beta ** 2)
    # Sirocco reads a model density as the lab-frame value and uses n' = n / gamma (define_wind.c)
    col = (ne * dr[:, None]).sum(axis=0) * SIGMA_T / gam    # sigma int n' dr, per theta cell
    k = np.array([np.sin(np.radians(inc)), 0.0, np.cos(np.radians(inc))])
    north = np.array([0.0, 0.0, 1.0]) - k[2] * k
    north /= np.linalg.norm(north)

    xg, wg = np.polynomial.legendre.leggauss(nmu)
    phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    I = Q = U = 0.0
    east = np.cross(north, k)       # only to measure U, which should be 0
    for hemi in (1, -1):            # the grid covers 0-90 deg and is mirrored
        for j in range(len(te) - 1):
            mu_a, mu_b = np.cos(np.radians(te[j + 1])), np.cos(np.radians(te[j]))
            mus = 0.5 * (mu_b - mu_a) * xg + 0.5 * (mu_b + mu_a)
            wmu = 0.5 * (mu_b - mu_a) * wg
            ct = hemi * mus[:, None]
            st = np.sqrt(1 - ct ** 2)
            rhat = np.stack([st * np.cos(phi), st * np.sin(phi), ct * np.ones_like(phi)], axis=-1)
            mu = rhat @ k
            mup = (mu - beta) / (1 - beta * mu)
            D = 1 / (gam * (1 - beta * mu))
            g = gam ** 2 * (1 - beta) ** 2 * D ** 3 * 3 / (16 * np.pi) * (1 + mup ** 2)
            P = (1 - mup ** 2) / (1 + mup ** 2)
            m = np.cross(rhat, k)
            m /= np.linalg.norm(m, axis=-1)[..., None]
            c2 = 2 * (m @ north) ** 2 - 1
            s2 = 2 * (m @ north) * (m @ east)
            w = col[j] * wmu[:, None] * (2 * np.pi / nphi)
            I += (w * g).sum()
            Q += (w * g * P * c2).sum()
            U += (w * g * P * s2).sum()

    # attenuation of the direct light along k (theta = inc, folded about the equator)
    th = inc if inc <= 90 else 180 - inc
    j = min(int(np.searchsorted(te, th, side='right')) - 1, len(te) - 2)
    tau = col[j] * gam * (1 - beta)
    return I, Q, U, tau


def doit(roots):
    '''
    Print the prediction for each run in roots (needs root.cells.npz)
    '''

    for root in roots:
        cells = np.load(root + '.cells.npz')
        beta = float(cells['beta'])
        print(f'{root}: beta = {beta}')
        for inc in ANGLES:
            I, Q, U, tau = predict(cells, inc, beta)
            print(f'  i = {inc:5.1f}  I_sc/I_0 = {I:.5e}  Q_sc/I_0 = {Q:+.5e}  U = {U:+.1e}  q_sc = {100*Q/I:+7.3f}%  tau_dir = {tau:.4f}')


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

    doit(words)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
