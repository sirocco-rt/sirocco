#!/usr/bin/env python
'''
Toy transport of polarized photons through a homologously expanding
ellipsoid with electron scattering and one Sobolev resonance line (Step 3
of the polarization project).


Command line usage (if any):

    None; this is a module.  step3_validate.py runs it.

Description:

    Lengths are in units of c t, so the flow velocity in units of c is just
    the position, v = x, and frequencies are in units of the line rest
    frequency nu_0.  Doppler shifts are first order in v/c (as in Sirocco's
    REL_MODE_LINEAR), and aberration is neglected:

        nu_cmf = nu (1 - n.x),      nu = nu_cmf / (1 - n.x)

    A photon moving along n from x reaches the resonance (nu_cmf = 1) after
    a distance s_res = 1 - 1/nu - n.x, and, since nu_cmf only decreases
    along a straight path, it meets the line at most once until it next
    scatters.  The Sobolev optical depth tau_line is the same in every
    direction in homologous flow, so the escape probability from the
    resonance region, beta = (1 - exp(-tau))/tau, is isotropic.

    Sources: a central point source of continuum (unpolarized, isotropic,
    flat in nu) and, optionally, unpolarized line emission spread uniformly
    through the ellipsoid.  The continuum extends 3 vmax beyond the spectral
    window on each side, so that Doppler shifts in electron scattering do not
    remove light from the window near its edges.

    The treatment of trapping in the resonance region is selected by
    ``trap``:

    explicit
        Every local scatter is followed: the photon scatters with the line
        phase matrix and escapes with probability beta, otherwise it is
        reabsorbed on the spot and scatters again.  This is the exact
        Sobolev result.  Peel-off happens at every local scatter, with the
        factor beta.
    single
        One polarized scatter per line interaction, as Sirocco's single
        scatter with thermal trapping does now.  This ignores the
        depolarization caused by the local scatters.
    effective
        One scatter per line interaction, but with the effective
        polarizability W2_eff = W2 F, where

            F = beta / (1 - (1 - beta) (7/10) W2)

        Each local scatter keeps a fraction (7/10) W2 of the K = 2
        (alignment) part of the radiation, which carries both the
        anisotropy and the linear polarization; an escaping photon has
        scattered k times with probability beta (1 - beta)^(k-1), and the
        sum over k gives F.  This is exact for homologous flow, where beta
        is isotropic, and cheap enough for C.

Primary routines:

    run, trans_phot, and the classes HomologousEllipsoid, SpectrumNu

Notes:

    Line emission is treated as the first event of a trapping sequence in
    the explicit mode; in the other modes it escapes directly, unpolarized.

History:

    261010 ksl Coding begun, with Claude Code

'''

import numpy as np

from . import photon
from .photon import P_INWIND, P_ESCAPE, P_TOO_MANY_SCATTERS, NRES_ES
from .stokes import normalize, dot, observer_direction, sky_north
from .es_scatter import scatter_es, peel_es
from .line_scatter import scatter_line, peel_line
from .toy_transport import Ellipsoid, LiveOrDie

RES_EPS = 1e-12      # s_res below this means the photon is at or past resonance
TRAP_MODES = ('explicit', 'single', 'effective')
MAX_LOCAL = 100000   # cap on local scatters in one resonance region


def p_escape(tau):
    '''Sobolev escape probability (1 - exp(-tau))/tau.'''
    tau = np.asarray(tau, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(tau > 1e-6, -np.expm1(-tau) / tau, 1.0 - 0.5 * tau)


def trap_factor(beta, w2):
    '''
    The fraction F = beta / (1 - (1 - beta)(7/10) W2) of the single-scatter
    polarization that survives trapping in a Sobolev resonance region with
    isotropic escape probability beta.
    '''
    return beta / (1.0 - (1.0 - beta) * 0.7 * w2)


class HomologousEllipsoid(Ellipsoid):
    '''
    A homologously expanding ellipsoid (semi-axes a, c in units of c t, which
    are also the maximum speeds in units of c) with uniform electron
    scattering opacity kappa (per unit c t) and a line of uniform Sobolev
    optical depth tau_line and polarizability w2.
    '''

    def __init__(self, a, c, kappa, tau_line, w2=1.0, line_frac=0.0, band=1.5, trap='explicit'):
        super().__init__(a, c, kappa)
        if trap not in TRAP_MODES:
            raise ValueError(f'trap must be one of {TRAP_MODES}')
        self.tau_line = float(tau_line)
        self.w2 = float(w2)
        self.beta = float(p_escape(tau_line)) if tau_line > 0 else 1.0
        self.line_frac = float(line_frac)
        self.vmax = max(a, c)
        self.band = band * self.vmax          # the spectra cover 1 - band .. 1 + band
        # The continuum extends 3 vmax beyond the spectra, so that electron
        # scattering (shifts of up to 2 vmax) moves as much light into the
        # spectral window as out of it.
        self.cont_half = self.band + 3.0 * self.vmax
        self.trap = trap
        # polarizability used for each line scatter
        if trap == 'effective':
            self.w2_scatter = self.w2 * trap_factor(self.beta, self.w2)
        else:
            self.w2_scatter = self.w2

    def s_res(self, x, n, nu):
        '''Distance to the line resonance, or inf if there is none ahead.'''
        s = 1.0 - 1.0 / nu - dot(n, x)
        return np.where(s > RES_EPS, s, np.inf)

    def tau_to_edge(self, x, n, nu):
        '''Optical depth (electron scattering plus line) from x to the edge along n.'''
        s_edge, _ = self.ds_to_edge(x, n)
        tau = self.kappa * s_edge
        if self.tau_line > 0:
            tau = tau + np.where(self.s_res(x, n, nu) < s_edge, self.tau_line, 0.0)
        return tau

    def generate(self, nphot, rng):
        '''
        Continuum photons from the center (flat in nu) and, if line_frac > 0,
        unpolarized line photons emitted uniformly through the ellipsoid.  The
        weights sum to 1.
        '''
        p = photon.make_photons(nphot)
        p['lmn'] = normalize(rng.normal(size=(nphot, 3)))
        nline = rng.binomial(nphot, self.line_frac) if self.line_frac > 0 else 0
        ncont = nphot - nline

        p['x'] = 0.0
        p['freq'][:ncont] = rng.uniform(1.0 - self.cont_half, 1.0 + self.cont_half, ncont)
        p['origin'][:ncont] = photon.PTYPE_STAR

        if nline:
            r = rng.random(nline) ** (1.0 / 3.0)
            pos = normalize(rng.normal(size=(nline, 3))) * r[:, None] * self.scale
            p['x'][ncont:] = pos
            p['freq'][ncont:] = 1.0 / (1.0 - dot(p['lmn'][ncont:], pos))
            p['origin'][ncont:] = photon.PTYPE_WIND
            p['nres'][ncont:] = 0              # born in the line

        p['freq_orig'] = p['freq']
        p['origin_orig'] = p['origin']
        p['w'] = p['w_orig'] = 1.0 / nphot
        photon.set_unpolarized(p)
        return p


class SpectrumNu:
    '''
    Peel-off spectrum for one observer, binned in y = (nu - 1)/vmax, with
    batch errors for the band-integrated q and u.  The analog of one xxspec[].
    '''

    def __init__(self, angle, vmax, nbins=60, ymax=1.5, phase=0.0):
        self.angle = angle
        self.lmn = observer_direction(angle, phase)
        self.north = sky_north(angle, phase)
        self.vmax = vmax
        self.edges = np.linspace(-ymax, ymax, nbins + 1)
        self.I = np.zeros(nbins)
        self.Q = np.zeros(nbins)
        self.U = np.zeros(nbins)
        self.I_line = np.zeros(nbins)    # light whose last interaction was the line
        self.Q_line = np.zeros(nbins)
        self.batches = []
        self._b = np.zeros(5)
        self.bin_batches = []                 # (I, Q, U) per bin for each batch
        self._bb = np.zeros((3, nbins))

    def tally(self, nu, wt, q, u, from_line=False):
        y = (nu - 1.0) / self.vmax
        k = np.searchsorted(self.edges, y) - 1
        ok = (k >= 0) & (k < len(self.I))
        k, wt, q, u = k[ok], wt[ok], q[ok], u[ok]
        np.add.at(self.I, k, wt)
        np.add.at(self.Q, k, wt * q)
        np.add.at(self.U, k, wt * u)
        np.add.at(self._bb[0], k, wt)
        np.add.at(self._bb[1], k, wt * q)
        np.add.at(self._bb[2], k, wt * u)
        if from_line:
            np.add.at(self.I_line, k, wt)
            np.add.at(self.Q_line, k, wt * q)
        self._b += (wt.sum(), (wt * q).sum(), (wt * u).sum(),
                    wt.sum() if from_line else 0.0, (wt * q).sum() if from_line else 0.0)

    def end_batch(self):
        self.batches.append(self._b.copy())
        self._b[:] = 0.0
        self.bin_batches.append(self._bb.copy())
        self._bb[:] = 0.0

    def integrated(self):
        '''
        Band-integrated (q, u, sq, su), and (q_line, sq_line) for the light
        whose last interaction was the line, with errors from the batches.
        '''
        b = np.array(self.batches)
        nb = len(b)
        tot = b.sum(axis=0)
        q, u = tot[1] / tot[0], tot[2] / tot[0]
        ql = tot[4] / tot[3] if tot[3] > 0 else np.nan

        def err(num, den, ratio, total):
            return np.std(b[:, num] * nb - ratio * b[:, den] * nb, ddof=1) / np.sqrt(nb) / total

        sq = err(1, 0, q, tot[0])
        su = err(2, 0, u, tot[0])
        sql = err(4, 3, ql, tot[3]) if tot[3] > 0 else np.nan
        return q, u, sq, su, ql, sql

    def binned(self):
        '''Bin centers in y, I, q, u, and the batch errors of q and u, for each bin.'''
        y = 0.5 * (self.edges[1:] + self.edges[:-1])
        b = np.array(self.bin_batches)
        nb = len(b)
        with np.errstate(divide='ignore', invalid='ignore'):
            q = self.Q / self.I
            u = self.U / self.I
            sq = np.std(b[:, 1] * nb - q * b[:, 0] * nb, axis=0, ddof=1) / np.sqrt(nb) / self.I
            su = np.std(b[:, 2] * nb - u * b[:, 0] * nb, axis=0, ddof=1) / np.sqrt(nb) / self.I
        return y, self.I, q, u, sq, su


def _observed_freq(x, n, nu_cmf):
    return nu_cmf / (1.0 - dot(n, x))


def extract(geom, p, spectra, kind):
    '''
    Peel photons p toward every spectrum.  kind is 'emit' (unpolarized,
    isotropic emission), 'es' (electron scatter), 'line' (a line scatter with
    polarizability geom.w2_scatter, multiplied by beta in the explicit mode),
    or 'line_emit' (line emission in the explicit mode, multiplied by beta).
    '''
    nphot = len(p)
    if nphot == 0:
        return
    x = p['x']
    nu_cmf = p['freq'] * (1.0 - dot(p['lmn'], x))
    zero = np.zeros(nphot)
    for spec in spectra:
        n_obs = np.broadcast_to(spec.lmn, x.shape)
        nu_obs = _observed_freq(x, n_obs, nu_cmf)
        if kind in ('line', 'line_emit'):
            nu_obs = _observed_freq(x, n_obs, 1.0)
        tau = geom.tau_to_edge(x, n_obs, nu_obs)
        if kind == 'emit':
            f, q, u = np.ones(nphot), zero, zero
        elif kind == 'es':
            f, q, u = peel_es(p, spec.lmn, spec.north)
        elif kind == 'line':
            f, q, u = peel_line(p, spec.lmn, spec.north, geom.w2_scatter)
            if geom.trap == 'explicit':
                f = f * geom.beta
        elif kind == 'line_emit':
            f, q, u = geom.beta * np.ones(nphot), zero, zero
        wt = p['w'] * f * np.exp(-tau)
        spec.tally(nu_obs, wt, q, u, from_line=kind != 'es' and kind != 'emit')


def _line_event(geom, p, rng, spectra):
    '''
    Process photons p that have just been absorbed in the line (they are at
    resonance, nu_cmf = 1).  On return they have left the resonance region
    with their new direction, Stokes parameters and observer-frame frequency.
    '''
    if geom.trap in ('single', 'effective'):
        extract(geom, p, spectra, 'line')
        scatter_line(p, rng, geom.w2_scatter)
    else:
        todo = np.arange(len(p))
        nloc = 0
        while todo.size and nloc < MAX_LOCAL:
            pt = p[todo]
            extract(geom, pt, spectra, 'line')
            scatter_line(pt, rng, geom.w2)
            pt['nrscat'] += 1
            p[todo] = pt
            todo = todo[rng.random(todo.size) >= geom.beta]
            nloc += 1
    p['freq'] = _observed_freq(p['x'], p['lmn'], 1.0)
    p['nres'] = 0


def _line_birth(geom, p, rng, spectra):
    '''Line photons at birth: direct escape, or (explicit mode) trapping.'''
    if geom.trap != 'explicit':
        extract(geom, p, spectra, 'emit')
        return
    extract(geom, p, spectra, 'line_emit')
    trapped = rng.random(len(p)) >= geom.beta
    pt = p[trapped]
    if len(pt):
        _line_event(geom, pt, rng, spectra)
        p[trapped] = pt


def trans_phot(geom, p, rng, spectra, live=None, maxscat=100000):
    '''
    Transport photons p through geom until they escape, peeling off toward
    every SpectrumNu in spectra at birth and at every interaction.  Escaping
    photons inside the spectral window, 1 - band < nu < 1 + band, are added to
    live (a LiveOrDie) if one is given.
    '''
    born_line = p['nres'] == 0
    pc = p[~born_line]
    extract(geom, pc, spectra, 'emit')
    pl = p[born_line]
    if len(pl):
        _line_birth(geom, pl, rng, spectra)
        p[born_line] = pl

    active = np.arange(len(p))
    while active.size:
        pa = p[active]
        x, n, nu = pa['x'], pa['lmn'], pa['freq']
        s_edge, istat = geom.ds_to_edge(x, n)
        s_es = -np.log(rng.random(active.size)) / geom.kappa if geom.kappa > 0 else np.full(active.size, np.inf)
        s_res = geom.s_res(x, n, nu) if geom.tau_line > 0 else np.full(active.size, np.inf)

        leave = (s_edge <= s_es) & (s_edge <= s_res)
        es = ~leave & (s_es < s_res)
        res = ~leave & ~es

        pa['x'][leave] += s_edge[leave, None] * n[leave]
        pa['istat'][leave] = istat[leave]

        if es.any():
            pe = pa[es]
            pe['x'] += s_es[es, None] * pe['lmn']
            nu_cmf = pe['freq'] * (1.0 - dot(pe['lmn'], pe['x']))
            pe['nscat'] += 1
            pe['nres'] = NRES_ES
            extract(geom, pe, spectra, 'es')
            scatter_es(pe, rng)
            pe['freq'] = _observed_freq(pe['x'], pe['lmn'], nu_cmf)
            pa[es] = pe

        if res.any():
            pr = pa[res]
            pr['x'] += s_res[res, None] * pr['lmn']
            absorbed = rng.random(len(pr)) < -np.expm1(-geom.tau_line)
            pl = pr[absorbed]
            if len(pl):
                pl['nscat'] += 1
                _line_event(geom, pl, rng, spectra)
                pr[absorbed] = pl
            pa[res] = pr

        pa['istat'][pa['nscat'] >= maxscat] = P_TOO_MANY_SCATTERS
        p[active] = pa
        if live is not None:
            esc = pa[pa['istat'] == P_ESCAPE]
            live.tally(esc[np.abs(esc['freq'] - 1.0) < geom.band])
        active = active[pa['istat'] == P_INWIND]
    return p


def run(geom, nphot, angles, rng, nbatches=20, max_batch=200000, nbins=60, live_bins=20):
    '''
    Run nphot photons through geom, peeling off toward each inclination in
    angles.  Returns (spectra, live).
    '''
    spectra = [SpectrumNu(a, geom.vmax, nbins=nbins, ymax=geom.band / geom.vmax) for a in angles]
    live = LiveOrDie(live_bins)
    nb = max(nbatches, int(np.ceil(nphot / max_batch)))
    sizes = np.full(nb, nphot // nb)
    sizes[:nphot % nb] += 1
    for size in sizes:
        p = geom.generate(size, rng)
        p['w'] *= size / nphot
        trans_phot(geom, p, rng, spectra, live)
        for spec in spectra:
            spec.end_batch()
    return spectra, live


def thin_line_brown_mclean(geom, angle):
    '''
    The optically thin prediction for line scattering of the continuum,
    integrated over the spectral window: q = W2 q_BM(kappa_eff), where
    kappa_eff = tau_line / (2 band) is the electron-scattering opacity that
    produces the same number of scatterings per unit radius relative to the
    continuum in the window.  (Valid to O(vmax), since the resonance radius
    changes with frequency as dr/dnu = 1/nu^2, and to O(tau_line).)
    '''
    from .toy_transport import brown_mclean
    a, c = geom.scale[0], geom.scale[2]
    q, tau0, gamma = brown_mclean(Ellipsoid(a, c, geom.tau_line / (2.0 * geom.band)), angle)
    return geom.w2 * q
