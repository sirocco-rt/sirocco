#!/usr/bin/env python
'''
Single-scattering estimate of P(lambda) for the make_cvcone.py models,
independent of Sirocco, to understand the wavelength dependence of P.

Rays leave the white dwarf (the hemisphere above the disk) and each annulus
of a standard blackbody disk with Eddington limb darkening, I ~ 1 + 1.5 mu,
as in Sirocco.  Along each ray, light is Thomson scattered by the wind
electrons (n_e from root.master.txt, from windsave2table) toward an observer
at inclination i, with intensity (3/16 pi)(1 + cos^2 Theta) and polarized
fraction sin^2 Theta / (1 + cos^2 Theta) perpendicular to the scattering
plane; the Stokes parameters are referred to the projected axis.  The wind
below the disk is the mirror image: it is included as the wind above the
disk seen from 180 - i, with the light blocked where its path to the
observer crosses the disk.  By default the light is attenuated by electron
scattering on the way in (exp(-tau) along the ray to the scattering point)
and on the way out (along the path to the observer, through the wind on
both sides of the disk plane), and the direct light of the central source is
attenuated along the radial line of sight; -no_atten turns this off.
Multiple scattering, and unpolarized emission by the wind, are neglected, so
the prediction is for photons scattered exactly once (Sirocco's _sc:1).

For each source the result is the scattered I and Q per unit source
luminosity per steradian, and so

    P(lambda) = sum L_lambda Q  /  (direct + sum L_lambda I)

where the direct light is the disk, L_face mu (1 + 1.5 mu) / 2 pi with mu =
cos i, plus the white dwarf, L / 4 pi times the fraction not hidden by the
disk, (1 + cos i) / 2.

usage: geom_pol.py root theta_min theta_max [-nray N] [-incs i1 i2 ...] [-no_atten] [-central]
       (default incs 10 45 80; -central skips the disk annuli, for models whose disk does not radiate)
       writes root_geompol.png, comparing with root.spec, and root.pred.txt, the
       central-source prediction that es_vs_prediction.py reads
'''
from types import SimpleNamespace
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii

sys.path.insert(0, os.path.join(os.environ['SIROCCO'], 'py_progs'))
import plot_pol

G, MSOL, YR, SIGMA_SB = 6.674e-8, 1.989e33, 3.15576e7, 5.6704e-5
H, C, K = 6.62607e-27, 2.99792458e10, 1.380649e-16
SIGMA_T = 6.6524587e-25
R_WD, M_WD, T_WD = 7e8, 0.8 * MSOL, 4e4
MDOT, R_DISK = 1e-8 * MSOL / YR, 2.4e10
INCS = [10, 45, 80]                     # the observers of pbl_es and eq_es


def doit(root, theta_min, theta_max, nray=3000, incs=INCS, no_atten=False, rdisk=2.4e10, seed=2, central=False):
    '''
    Make the single-scattering prediction for the wind of root, between
    theta_min and theta_max (deg from the pole); the parameters are those
    of the command line, described above
    '''

    a = SimpleNamespace(**locals())         # the code below refers to the parameters as a.<name>
    R_DISK = a.rdisk
    INCS = [int(x) for x in a.incs]
    rng = np.random.default_rng(a.seed)

    m = ascii.read(a.root + '.master.txt')
    r_edge = np.unique(np.asarray(m['r']))
    t_edge = np.unique(np.asarray(m['theta']))
    ne_grid = np.zeros((len(r_edge), len(t_edge)))
    for row in m:
        if row['inwind'] == 0:
            ne_grid[np.searchsorted(r_edge, row['r']), np.searchsorted(t_edge, row['theta'])] = row['ne']
    r_in, r_out = r_edge[0], r_edge[-2]
    tmin, tmax = np.radians(a.theta_min), np.radians(a.theta_max)

    s = np.concatenate([[0], np.geomspace(1e7, 4 * r_out, 2500)])
    ds = np.diff(s)
    sm = 0.5 * (s[1:] + s[:-1])


    so = np.concatenate([[0], np.geomspace(1e7, 4 * r_out, 400)])     # steps for the path out to the observer
    dso = np.diff(so)
    smo = 0.5 * (so[1:] + so[:-1])


    def dtau(x):
        '''electron-scattering dtau per cm at positions x (..., 3); the wind is mirror-symmetric about the disk plane'''
        r = np.linalg.norm(x, axis=-1)
        th = np.arccos(np.clip(np.abs(x[..., 2]) / np.maximum(r, 1e-30), 0, 1))
        inside = (r > r_in) & (r < r_out) & (th > tmin) & (th < tmax)
        ii = np.clip(np.searchsorted(r_edge, r) - 1, 0, ne_grid.shape[0] - 1)
        jj = np.clip(np.searchsorted(t_edge, np.degrees(th)) - 1, 0, ne_grid.shape[1] - 1)
        return np.where(inside, ne_grid[ii, jj], 0.0) * SIGMA_T


    def tau_out(xs, k, chunk=20000):
        '''tau from each point xs to infinity along k'''
        out = np.empty(len(xs))
        for i0 in range(0, len(xs), chunk):
            p = xs[i0:i0 + chunk]
            out[i0:i0 + chunk] = (dtau(p[:, None, :] + smo[None, :, None] * k[None, None, :]) * dso[None, :]).sum(1)
        return out


    def tau_los(inc):
        '''tau from the centre to the observer at inclination inc (deg), for the direct light'''
        ti = np.radians(inc)
        return (dtau(smo[:, None] * np.array([np.sin(ti), 0.0, np.cos(ti)])[None, :]) * dso).sum()


    def limb_mu(n):
        '''mu drawn from p(mu) ~ mu (1 + 1.5 mu), as Sirocco's randvcos'''
        out = np.empty(0)
        while len(out) < n:
            x = rng.random(2 * n)
            y = rng.random(2 * n) * 2.5
            out = np.concatenate([out, x[y < x * (1 + 1.5 * x)]])
        return out[:n]


    def observers():
        '''unit vectors k and sky North for observers at INCS (above) and 180 - INCS (the mirrored wind)'''
        obs = []
        for inc in INCS:
            for i in (inc, 180 - inc):
                ti = np.radians(i)
                k = np.array([np.sin(ti), 0.0, np.cos(ti)])
                nth = np.array([0.0, 0.0, 1.0]) - k[2] * k
                nth /= np.linalg.norm(nth)
                obs.append((k, nth))
        return obs


    OBS = observers()


    def scatter(p0, normals, chunk=2000):
        '''scattered (I, Q) per unit luminosity per sr for each observer, from rays from p0 about normals'''
        n = len(p0)
        return sum(len(p0[i:i + chunk]) * scatter1(p0[i:i + chunk], normals[i:i + chunk])
                   for i in range(0, n, chunk)) / n


    def scatter1(p0, normals):
        '''scatter() for one chunk of rays'''
        n = len(p0)
        mu = limb_mu(n)
        phi = 2 * np.pi * rng.random(n)
        t1 = np.cross(normals, [0.0, 1.0, 0.0])
        bad = np.linalg.norm(t1, axis=1) < 1e-6
        t1[bad] = np.cross(normals[bad], [1.0, 0.0, 0.0])
        t1 /= np.linalg.norm(t1, axis=1)[:, None]
        t2 = np.cross(normals, t1)
        st = np.sqrt(1 - mu ** 2)
        d = mu[:, None] * normals + (st * np.cos(phi))[:, None] * t1 + (st * np.sin(phi))[:, None] * t2
        x = p0[:, None, :] + sm[None, :, None] * d[:, None, :]
        r = np.linalg.norm(x, axis=2)
        z = x[:, :, 2]
        th = np.arccos(np.clip(z / np.maximum(r, 1e-30), -1, 1))
        inside = (r > r_in) & (r < r_out) & (th > tmin) & (th < tmax) & (z > 0)
        ii = np.clip(np.searchsorted(r_edge, r) - 1, 0, ne_grid.shape[0] - 1)
        jj = np.clip(np.searchsorted(t_edge, np.degrees(th)) - 1, 0, ne_grid.shape[1] - 1)
        w = np.where(inside, ne_grid[ii, jj], 0.0) * SIGMA_T * ds[None, :]      # dtau at each step
        sel = w > 0
        if not a.no_atten:
            w = w * np.exp(-(np.cumsum(w, axis=1) - 0.5 * w))      # what reaches the step, to the middle of it
        w, xs = w[sel], x[sel]
        dd = np.broadcast_to(d[:, None, :], x.shape)[sel]
        res = []
        for k, nth in OBS:
            c = dd @ k
            e = np.cross(dd, k)                         # E of the scattered light, normal to the scattering plane
            en = np.linalg.norm(e, axis=1)
            e = e / np.maximum(en, 1e-30)[:, None]
            cos2chi = 2 * (e @ nth) ** 2 - 1
            # the path to the observer from the mirrored (lower) wind crosses the disk plane
            if k[2] < 0:
                tcross = xs[:, 2] / -k[2]
                rc = np.hypot(xs[:, 0] + tcross * k[0], xs[:, 1] + tcross * k[1])
                vis = rc > R_DISK
            else:
                vis = np.ones(len(w), dtype=bool)
            if not a.no_atten:
                vis = vis * np.exp(-tau_out(xs, k))
            wI = w * 3 / (16 * np.pi) * (1 + c ** 2) * vis
            wQ = w * 3 / (16 * np.pi) * (1 - c ** 2) * cos2chi * vis
            res.append((wI.sum() / n, wQ.sum() / n))
        return np.array(res)          # (nobs, 2)


    edges = np.geomspace(R_WD, R_DISK, 41)
    Rc = np.sqrt(edges[1:] * edges[:-1])
    tstar = (3 * G * M_WD * MDOT / (8 * np.pi * SIGMA_SB * R_WD ** 3)) ** 0.25
    T_ann = tstar * (R_WD / Rc) ** 0.75 * (1 - np.sqrt(R_WD / Rc)) ** 0.25
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)

    ann = []
    for R in ([] if a.central else Rc):
        ph = 2 * np.pi * rng.random(a.nray)            # all azimuths around the disk
        p0 = np.stack([R * np.cos(ph), R * np.sin(ph), np.ones(a.nray)], axis=1)
        ann.append(scatter(p0, np.tile([0, 0, 1.0], (a.nray, 1))))
    ann = np.array(ann) if ann else np.zeros((len(Rc), len(OBS), 2))     # (nann, nobs, 2)
    if a.central:
        np.seterr(all='ignore')
    cth = rng.random(a.nray)
    phw = 2 * np.pi * rng.random(a.nray)
    nrm = np.stack([np.sqrt(1 - cth ** 2) * np.cos(phw), np.sqrt(1 - cth ** 2) * np.sin(phw), cth], axis=1)
    wd = scatter(R_WD * nrm * 1.0001, nrm)

    # the wavelength range, and the bands used in the printed tables, follow the spectrum file
    units, freq, wave, cols = plot_pol.read_spec(a.root + '.spec')
    wlo, whi = wave.min(), wave.max()
    lam = np.linspace(wlo, whi, 200)
    lc = lam * 1e-8
    wedges = np.linspace(wlo, whi, 4)
    BANDS = [(wedges[i], wedges[i + 1]) for i in range(3)]          # thirds of the range
    B1, B2 = (wlo, wlo + 0.25 * (whi - wlo)), (whi - 0.25 * (whi - wlo), whi)   # lowest and highest quarters


    def bb(T):
        return 2 * H * C ** 2 / lc ** 5 / np.expm1(H * C / (lc * K * T))


    L_wd_up = 0.5 * 4 * np.pi * R_WD ** 2 * np.pi * bb(T_WD)          # the hemisphere above the disk
    L_face = np.array([area[k] * np.pi * bb(T_ann[k]) for k in range(len(Rc))])

    fac = plot_pol.lnu_dnu(units, freq, False)
    nb = 100
    n = (len(wave) // nb) * nb
    wb = wave[:n].reshape(-1, nb).mean(1)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    print('  i   P_pred(%.0f-%.0f)  P_pred(%.0f-%.0f)   Sirocco(%.0f-%.0f)  Sirocco(%.0f-%.0f)' % (B1 + B2 + B1 + B2))
    for k, inc in enumerate(INCS):
        mu = np.cos(np.radians(inc))
        I_sc = np.zeros_like(lam)
        Q_sc = np.zeros_like(lam)
        for o in (2 * k, 2 * k + 1):                     # wind above, mirrored wind below
            I_sc += L_wd_up * wd[o, 0] + (L_face * ann[:, o, 0][:, None]).sum(0)
            Q_sc += L_wd_up * wd[o, 1] + (L_face * ann[:, o, 1][:, None]).sum(0)
        att = 1.0 if a.no_atten else np.exp(-tau_los(inc))       # the disk's direct light is not attenuated
        direct = att * 2 * L_wd_up / (4 * np.pi) * (1 + mu) / 2 + L_face.sum(0) * mu * (1 + 1.5 * mu) / (2 * np.pi)
        P = 100 * np.abs(Q_sc) / (direct + I_sc)
        s = []
        for lo, hi in BANDS:
            m = (lam > lo) & (lam < hi)
            s.append('%.3f %6.2f' % (I_sc[m].sum() / (direct[m] + I_sc[m]).sum(), 100 * Q_sc[m].sum() / I_sc[m].sum()))
        print('%3d  predicted frac_sc, P_sc(%%) by band:  %s' % (inc, '  |  '.join(s)))
        name = 'A%02dP0.50' % inc
        col = 'C%d' % k
        ax[0].plot(lam, P, color=col, lw=2, ls='--', label='single scattering %d$^\\circ$' % inc)
        if name in cols:
            I = (cols[name] * fac)[:n].reshape(-1, nb).sum(1)
            Q = (cols[name + '_Q'] * fac)[:n].reshape(-1, nb).sum(1)
            U = (cols[name + '_U'] * fac)[:n].reshape(-1, nb).sum(1)
            Ps = 100 * np.hypot(Q, U) / I
            ax[0].step(wb, Ps, where='mid', color=col, lw=0.8, label='Sirocco %d$^\\circ$' % inc)
            b1, b2 = (wb > B1[0]) & (wb < B1[1]), (wb > B2[0]) & (wb < B2[1])
            l1, l2 = (lam >= B1[0]) & (lam <= B1[1]), (lam >= B2[0]) & (lam <= B2[1])
            print('%3d   %8.2f           %8.2f            %8.2f            %8.2f' % (
                inc, P[l1].mean(), P[l2].mean(), 100 * np.hypot(Q[b1].sum(), U[b1].sum()) / I[b1].sum(),
                100 * np.hypot(Q[b2].sum(), U[b2].sum()) / I[b2].sum()))
        if inc == INCS[-1]:
            # how polarized is the scattered light from each kind of source at this inclination?
            o1, o2 = 2 * k, 2 * k + 1
            pwd = (wd[o1, 1] + wd[o2, 1]) / (wd[o1, 0] + wd[o2, 0])
            pann = (ann[:, o1, 1] + ann[:, o2, 1]) / (ann[:, o1, 0] + ann[:, o2, 0])
            ax[1].semilogx(Rc / R_WD, -100 * pann, 'o-', ms=3, label='disk annulus at R')
            ax[1].axhline(-100 * pwd, color='C3', ls='--', label='white dwarf')
            fsc = (ann[:, o1, 0] + ann[:, o2, 0])
            ax2 = ax[1].twinx()
            ax2.semilogx(Rc / R_WD, fsc / fsc.max(), color='0.5', lw=0.8)
            ax2.set_ylabel('scattered I per unit L (relative)', color='0.5')
    print('central source only (as polar_bl)%s:' % ('' if a.no_atten else ', with attenuation'))
    pred = open(a.root + '.pred.txt', 'w')
    for k, inc in enumerate(INCS):
        o1, o2 = 2 * k, 2 * k + 1
        # per unit total L: each hemisphere (half the light) illuminates its own half of the wind
        Isc, Qsc = 0.5 * (wd[o1, 0] + wd[o2, 0]), 0.5 * (wd[o1, 1] + wd[o2, 1])
        att = 1.0 if a.no_atten else np.exp(-tau_los(inc))
        direct = att / (4 * np.pi) * (1 + np.cos(np.radians(inc))) / 2
        line = '%3d  frac_sc %.3f  P_sc %6.2f%%  P %6.3f%%' % (inc, Isc / (direct + Isc), 100 * Qsc / Isc,
                                                             100 * Qsc / (direct + Isc))
        print(line)
        pred.write(line + '\n')
    pred.close()
    print('wrote %s.pred.txt' % a.root)
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Wavelength ($\\AA$)')
    ax[0].set_ylabel('P (%)')
    ax[0].legend(fontsize=7, ncol=2)
    ax[0].set_title('%s: Sirocco vs single scattering of WD + disk light (no fitting)' % a.root, fontsize=9)
    ax[1].set_xlabel('disk radius (R$_{wd}$)')
    ax[1].set_ylabel('polarization of the scattered light at %d$^\\circ$ (%%)' % INCS[-1])
    ax[1].legend(fontsize=8, loc='center left')
    ax[1].set_title('polarization of the scattered light, by source', fontsize=9)
    fig.tight_layout()
    fig.savefig(a.root + '_geompol.png', dpi=110)
    print(a.root + '_geompol.png')


def steer(argv):
    '''
    Process the command line
    '''

    words = []
    x = {}

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] in ('-nray', '-seed'):
            x[argv[i][1:]] = int(argv[i + 1])
            i += 1
        elif argv[i] == '-rdisk':
            x['rdisk'] = float(argv[i + 1])
            i += 1
        elif argv[i] == '-incs':
            incs = []
            while i + 1 < len(argv) and is_number(argv[i + 1]):
                i += 1
                incs.append(float(argv[i]))
            x['incs'] = incs
        elif argv[i] in ('-no_atten', '-central'):
            x[argv[i][1:]] = True
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            words.append(argv[i])
        i += 1

    if len(words) != 3:
        print(__doc__)
        return

    doit(words[0], float(words[1]), float(words[2]), **x)


def is_number(word):
    '''True if word can be read as a number'''
    try:
        float(word)
        return True
    except ValueError:
        return False


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
