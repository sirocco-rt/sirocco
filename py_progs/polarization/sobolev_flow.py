#!/usr/bin/env python
'''
Toy transport of polarized photons through an expanding ellipsoid with
electron scattering and one Sobolev resonance line (Step 3 of the
polarization project).


Command line usage (if any):

    None; this is a module.  step3_validate.py and step3b_validate.py run it.

Description:

    Lengths are in units of c t (homologous flow) or of an arbitrary
    reference length (radial flow); velocities are in units of c and
    frequencies in units of the line rest frequency nu_0.  Doppler shifts are
    first order in v/c (as in Sirocco's REL_MODE_LINEAR), and aberration is
    neglected:

        nu_cmf = nu (1 - n.v(x)),      nu = nu_cmf / (1 - n.v(x))

    Both flows accelerate in every direction (n.grad(n.v) > 0), so nu_cmf
    falls monotonically along a straight path and a photon meets the line at
    most once until it next scatters.

    Two geometries, both with a central point source of continuum
    (unpolarized, isotropic and flat in nu) and, optionally, unpolarized line
    emission spread uniformly through the line-forming region:

    HomologousEllipsoid
        v = x.  The Sobolev optical depth is the same in every direction, so
        the escape probability beta = (1 - exp(-tau))/tau is isotropic.  The
        resonance point is analytic: s_res = 1 - 1/nu - n.x.
    RadialFlowEllipsoid
        v = v0 (r/a)^alpha in the radial direction, with the line confined to
        r > r_in.  The Sobolev optical depth depends on direction,

            tau(x, n) = tau_line (r/a)^(1 - alpha) / (1 - (1 - alpha) mu^2),

        with mu = n.r_hat, so the escape probability is anisotropic: alpha < 1
        traps photons moving radially, and alpha > 1 photons moving
        tangentially.  alpha = 1 is homologous.  The resonance point is found
        by bisection.

    The treatment of trapping in the resonance region is selected by
    ``trap``:

    explicit
        Every local scatter is followed: the photon scatters with the line
        phase matrix and escapes in its new direction n with probability
        beta(x, n), otherwise it is reabsorbed on the spot and scatters
        again.  This is the exact Sobolev result.  Peel-off happens at every
        local scatter, with the factor beta(x, n_obs).
    single
        One polarized scatter per line interaction, as Sirocco's line
        scatter with thermal trapping does now: the direction is drawn from
        beta(x, n) times the phase function by rejection, and the peel-off
        weight is multiplied by beta(x, n_obs) and by the number of tries
        (cf. nnscat/p_norm in extract.c).  This ignores the depolarization
        caused by the local scatters.
    effective
        As single, but with the effective polarizability W2_eff = W2 F, where

            F = beta_bar / (1 - (1 - beta_bar) (7/10) W2)

        and beta_bar is the escape probability averaged over all directions.
        Each local scatter keeps a fraction (7/10) W2 of the K = 2
        (alignment) part of the radiation, which carries both the anisotropy
        and the linear polarization; an escaping photon has scattered k
        times with probability beta (1 - beta)^(k-1), and the sum over k
        gives F.  This is exact when beta is isotropic.

Primary routines:

    run, trans_phot, extract, trap_factor, thin_line_brown_mclean, and the
    classes HomologousEllipsoid, RadialFlowEllipsoid, SpectrumNu

Notes:

    Line emission is treated as the first event of a trapping sequence in
    the explicit mode; in the other modes it escapes directly, unpolarized,
    in a direction drawn from beta(x, n).

History:

    261010 ksl Coding begun, with Claude Code

'''

import numpy as np

from . import photon
from .photon import P_INWIND, P_ESCAPE, P_TOO_MANY_SCATTERS, NRES_ES
from .stokes import normalize, dot, observer_direction, sky_north
from .es_scatter import scatter_es, peel_es
from .line_scatter import peel_line, sample_line_direction, scattered_state
from .toy_transport import Ellipsoid, LiveOrDie

RES_EPS = 1e-12      # nu_cmf - 1 or s_res below this: at or past resonance
TRAP_MODES = ('explicit', 'single', 'effective')
MAX_LOCAL = 100000   # cap on local scatters in one resonance region
NMU_BETA = 16        # Gauss-Legendre points for the angle-averaged escape probability


def p_escape(tau):
    '''Sobolev escape probability (1 - exp(-tau))/tau.'''
    tau = np.asarray(tau, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(tau > 1e-6, -np.expm1(-tau) / tau, 1.0 - 0.5 * tau)


def trap_factor(beta, w2):
    '''
    The fraction F = beta / (1 - (1 - beta)(7/10) W2) of the single-scatter
    polarization that survives trapping in a Sobolev resonance region with
    escape probability beta.
    '''
    return beta / (1.0 - (1.0 - beta) * 0.7 * w2)


class _SobolevEllipsoid(Ellipsoid):
    '''Common parts of the two flows; see the module documentation.'''

    isotropic_escape = False

    def _setup(self, tau_line, w2, line_frac, band, trap):
        if trap not in TRAP_MODES:
            raise ValueError(f'trap must be one of {TRAP_MODES}')
        self.tau_line = float(tau_line)
        self.w2 = float(w2)
        self.line_frac = float(line_frac)
        self.band = band * self.vmax          # the spectra cover 1 - band .. 1 + band
        # The continuum extends 3 vmax beyond the spectra, so that electron
        # scattering (shifts of up to 2 vmax) moves as much light into the
        # spectral window as out of it.
        self.cont_half = self.band + 3.0 * self.vmax
        self.trap = trap

    def nu_cmf(self, x, n, nu):
        return nu * (1.0 - dot(n, self.velocity(x)))

    def nu_obs(self, x, n, nu_cmf):
        return nu_cmf / (1.0 - dot(n, self.velocity(x)))

    def beta(self, x, n):
        '''Sobolev escape probability at x in direction n.'''
        return p_escape(self.tau_sobolev(x, n))

    def w2_scatter(self, x):
        '''The polarizability used for a single line scatter at x.'''
        if self.trap == 'effective':
            return self.w2 * trap_factor(self.beta_mean(x), self.w2)
        return np.full(len(x), self.w2)

    def tau_to_edge(self, x, n, nu):
        '''Optical depth (electron scattering plus line) from x to the edge along n.'''
        s_edge, _ = self.ds_to_edge(x, n)
        tau = self.kappa * s_edge
        if self.tau_line > 0:
            s = self.s_res(x, n, nu)
            hit = s < s_edge
            if hit.any():
                xr = x[hit] + s[hit, None] * n[hit]
                tau = tau.copy()
                tau[hit] += self.tau_sobolev(xr, n[hit])
        return tau

    def in_line_region(self, x):
        return np.ones(len(x), dtype=bool)

    def generate(self, nphot, rng):
        '''
        Continuum photons from the center (flat in nu) and, if line_frac > 0,
        unpolarized line photons emitted uniformly through the line-forming
        region.  The weights sum to 1.
        '''
        p = photon.make_photons(nphot)
        p['lmn'] = normalize(rng.normal(size=(nphot, 3)))
        nline = rng.binomial(nphot, self.line_frac) if self.line_frac > 0 else 0
        ncont = nphot - nline

        p['x'] = 0.0
        p['freq'][:ncont] = rng.uniform(1.0 - self.cont_half, 1.0 + self.cont_half, ncont)
        p['origin'][:ncont] = photon.PTYPE_STAR

        if nline:
            pos = np.empty((0, 3))
            while len(pos) < nline:
                m = 2 * (nline - len(pos)) + 10
                r = rng.random(m) ** (1.0 / 3.0)
                trial = normalize(rng.normal(size=(m, 3))) * r[:, None] * self.scale
                pos = np.vstack([pos, trial[self.in_line_region(trial)]])
            pos = pos[:nline]
            p['x'][ncont:] = pos
            p['freq'][ncont:] = self.nu_obs(pos, p['lmn'][ncont:], 1.0)
            p['origin'][ncont:] = photon.PTYPE_WIND
            p['nres'][ncont:] = 0              # born in the line

        p['freq_orig'] = p['freq']
        p['origin_orig'] = p['origin']
        p['w'] = p['w_orig'] = 1.0 / nphot
        photon.set_unpolarized(p)
        return p


class HomologousEllipsoid(_SobolevEllipsoid):
    '''
    A homologously expanding ellipsoid (semi-axes a, c in units of c t, which
    are also the maximum speeds in units of c) with uniform electron
    scattering opacity kappa (per unit c t) and a line of uniform Sobolev
    optical depth tau_line and polarizability w2.
    '''

    isotropic_escape = True

    def __init__(self, a, c, kappa, tau_line, w2=1.0, line_frac=0.0, band=1.5, trap='explicit'):
        super().__init__(a, c, kappa)
        self.vmax = max(a, c)
        self._setup(tau_line, w2, line_frac, band, trap)
        self.beta0 = float(p_escape(tau_line)) if tau_line > 0 else 1.0

    def velocity(self, x):
        return x

    def s_res(self, x, n, nu):
        '''Distance to the line resonance, or inf if there is none ahead.'''
        s = 1.0 - 1.0 / nu - dot(n, x)
        return np.where(s > RES_EPS, s, np.inf)

    def tau_sobolev(self, x, n):
        return np.full(len(x), self.tau_line)

    def beta_mean(self, x):
        return np.full(len(x), self.beta0)


class RadialFlowEllipsoid(_SobolevEllipsoid):
    '''
    An ellipsoid (semi-axes a, c) in a radial outflow v = v0 (r/a)^alpha,
    with uniform electron scattering opacity kappa and a line whose Sobolev
    optical depth is tau(x, n) = tau_line (r/a)^(1-alpha) / (1 - (1-alpha) mu^2)
    for r > r_in (and zero inside r_in).
    '''

    def __init__(self, a, c, v0, alpha, kappa, tau_line, r_in=None, w2=1.0, line_frac=0.0,
                 band=1.5, trap='explicit'):
        super().__init__(a, c, kappa)
        if alpha <= 0:
            raise ValueError('alpha must be > 0, so that the flow accelerates in every direction')
        self.v0 = float(v0)
        self.alpha = float(alpha)
        self.a_ref = float(a)
        self.r_in = 0.2 * min(a, c) if r_in is None else float(r_in)
        self.vmax = self.v0 * (max(a, c) / a) ** self.alpha
        self._setup(tau_line, w2, line_frac, band, trap)
        self._mu, self._wmu = np.polynomial.legendre.leggauss(NMU_BETA)

    def velocity(self, x):
        r = np.linalg.norm(x, axis=-1)
        with np.errstate(divide='ignore', invalid='ignore'):
            vr = np.where(r > 0, self.v0 * (r / self.a_ref) ** self.alpha / r, 0.0)
        return vr[..., None] * x

    def in_line_region(self, x):
        return np.linalg.norm(x, axis=-1) > self.r_in

    def tau_sobolev(self, x, n):
        r = np.linalg.norm(x, axis=-1)
        mu = dot(n, x) / np.maximum(r, 1e-300)
        tau = self.tau_line * (r / self.a_ref) ** (1.0 - self.alpha) / (1.0 - (1.0 - self.alpha) * mu * mu)
        return np.where(r > self.r_in, tau, 0.0)

    def beta_mean(self, x):
        '''The escape probability averaged over all directions at x.'''
        r = np.linalg.norm(x, axis=-1)
        t0 = self.tau_line * (r / self.a_ref) ** (1.0 - self.alpha)
        tau = t0[:, None] / (1.0 - (1.0 - self.alpha) * self._mu[None, :] ** 2)
        return 0.5 * np.sum(self._wmu * p_escape(tau), axis=1)

    def s_res(self, x, n, nu):
        '''
        Distance to the line resonance by bisection, or inf if there is none
        ahead inside the ellipsoid, or if it lies inside r_in.
        '''
        s_edge, _ = self.ds_to_edge(x, n)
        f0 = self.nu_cmf(x, n, nu) - 1.0
        f1 = self.nu_cmf(x + s_edge[:, None] * n, n, nu) - 1.0
        out = np.full(len(x), np.inf)
        k = np.where((f0 > RES_EPS) & (f1 < 0))[0]
        if k.size:
            lo = np.zeros(k.size)
            hi = s_edge[k].copy()
            xk, nk, nuk = x[k], n[k], nu[k]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                above = self.nu_cmf(xk + mid[:, None] * nk, nk, nuk) > 1.0
                lo = np.where(above, mid, lo)
                hi = np.where(above, hi, mid)
            s = 0.5 * (lo + hi)
            xr = xk + s[:, None] * nk
            out[k] = np.where(self.in_line_region(xr), s, np.inf)
        return out


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


def extract(geom, p, spectra, kind, w2=None, ntries=None):
    '''
    Peel photons p toward every spectrum.  kind is

    'emit'       unpolarized, isotropic emission (continuum at birth)
    'es'         an electron scatter
    'line'       a line scatter with polarizability w2 (one value per photon)
    'line_emit'  unpolarized line emission

    For 'line' and 'line_emit' the weight is multiplied by the escape
    probability toward the observer, beta(x, n_obs).  In the single and
    effective modes, where one event stands for a whole trapping sequence,
    it is also multiplied by ntries, the number of directions tried before
    one escaped (cf. nnscat/p_norm in extract.c); ntries = 1/beta exactly
    when the escape probability is isotropic.
    '''
    nphot = len(p)
    if nphot == 0:
        return
    x = p['x']
    zero = np.zeros(nphot)
    line_kind = kind in ('line', 'line_emit')
    nu_cmf = np.ones(nphot) if line_kind else geom.nu_cmf(x, p['lmn'], p['freq'])
    for spec in spectra:
        n_obs = np.broadcast_to(spec.lmn, x.shape)
        nu_obs = geom.nu_obs(x, n_obs, nu_cmf)
        tau = geom.tau_to_edge(x, n_obs, nu_obs)
        if kind == 'emit':
            f, q, u = np.ones(nphot), zero, zero
        elif kind == 'es':
            f, q, u = peel_es(p, spec.lmn, spec.north)
        elif kind == 'line':
            f, q, u = peel_line(p, spec.lmn, spec.north, w2)
        else:
            f, q, u = np.ones(nphot), zero, zero
        if line_kind:
            f = f * geom.beta(x, n_obs)
            if ntries is not None:
                f = f * ntries
        wt = p['w'] * f * np.exp(-tau)
        spec.tally(nu_obs, wt, q, u, from_line=line_kind)


def _draw_escaping(geom, rng, x, n, e1, stokes, w2):
    '''
    Draw a direction from beta(x, n') times the line phase function (w2 = 0
    for isotropic emission) by rejection.  Returns the directions and the
    number of tries for each photon.
    '''
    nphot = len(x)
    n_new = np.empty((nphot, 3))
    tries = np.zeros(nphot)
    todo = np.arange(nphot)
    while todo.size:
        trial = sample_line_direction(rng, n[todo], e1[todo], stokes[todo], w2[todo])
        tries[todo] += 1
        ok = rng.random(todo.size) < geom.beta(x[todo], trial)
        n_new[todo[ok]] = trial[ok]
        todo = todo[~ok]
    return n_new, tries


def _apply_line_scatter(p, n_new, w2):
    _, e_perp, st = scattered_state(p['lmn'], p['pol_ref'], p['stokes'], n_new, w2)
    p['lmn'] = n_new
    p['pol_ref'] = e_perp
    p['stokes'] = st


def _line_event(geom, p, rng, spectra):
    '''
    Process photons p that have just been absorbed in the line (they are at
    resonance, nu_cmf = 1).  On return they have left the resonance region
    with their new direction, Stokes parameters and observer-frame frequency.
    '''
    x = p['x']
    if geom.trap in ('single', 'effective'):
        w2 = geom.w2_scatter(x)
        if geom.isotropic_escape:
            n_new = sample_line_direction(rng, p['lmn'], p['pol_ref'], p['stokes'], w2)
            ntries = 1.0 / geom.beta_mean(x)
        else:
            n_new, ntries = _draw_escaping(geom, rng, x, p['lmn'], p['pol_ref'], p['stokes'], w2)
        extract(geom, p, spectra, 'line', w2=w2, ntries=ntries)
        _apply_line_scatter(p, n_new, w2)
    else:
        w2 = np.full(len(p), geom.w2)
        todo = np.arange(len(p))
        nloc = 0
        while todo.size and nloc < MAX_LOCAL:
            pt = p[todo]
            extract(geom, pt, spectra, 'line', w2=w2[todo])
            n_new = sample_line_direction(rng, pt['lmn'], pt['pol_ref'], pt['stokes'], w2[todo])
            _apply_line_scatter(pt, n_new, w2[todo])
            pt['nrscat'] += 1
            escaped = rng.random(todo.size) < geom.beta(pt['x'], n_new)
            p[todo] = pt
            todo = todo[~escaped]
            nloc += 1
    p['freq'] = geom.nu_obs(p['x'], p['lmn'], 1.0)
    p['nres'] = 0


def _line_birth(geom, p, rng, spectra):
    '''
    Line photons at birth (unpolarized, isotropic emission at resonance).
    explicit: escape in the emitted direction with probability beta(x, n),
    otherwise start a trapping sequence.  single/effective: escape directly
    in a direction drawn from beta(x, n).
    '''
    x = p['x']
    if geom.trap != 'explicit':
        if geom.isotropic_escape:
            ntries = 1.0 / geom.beta_mean(x)
        else:
            zero_w2 = np.zeros(len(p))
            n_new, ntries = _draw_escaping(geom, rng, x, p['lmn'], p['pol_ref'], p['stokes'], zero_w2)
            p['lmn'] = n_new
            photon.set_unpolarized(p)
        extract(geom, p, spectra, 'line_emit', ntries=ntries)
    else:
        extract(geom, p, spectra, 'line_emit')
        trapped = rng.random(len(p)) >= geom.beta(x, p['lmn'])
        pt = p[trapped]
        if len(pt):
            _line_event(geom, pt, rng, spectra)
            p[trapped] = pt
    p['freq'] = geom.nu_obs(p['x'], p['lmn'], 1.0)


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
            nu_cmf = geom.nu_cmf(pe['x'], pe['lmn'], pe['freq'])
            pe['nscat'] += 1
            pe['nres'] = NRES_ES
            extract(geom, pe, spectra, 'es')
            scatter_es(pe, rng)
            pe['freq'] = geom.nu_obs(pe['x'], pe['lmn'], nu_cmf)
            pa[es] = pe

        if res.any():
            pr = pa[res]
            pr['x'] += s_res[res, None] * pr['lmn']
            absorbed = rng.random(len(pr)) < -np.expm1(-geom.tau_sobolev(pr['x'], pr['lmn']))
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
    The optically thin prediction for line scattering of the continuum in a
    HomologousEllipsoid, integrated over the spectral window: q = W2
    q_BM(kappa_eff), where kappa_eff = tau_line / (2 band) is the electron-
    scattering opacity that produces the same number of scatterings per unit
    radius relative to the continuum in the window.  (Valid to O(vmax),
    since the resonance radius changes with frequency as dr/dnu = 1/nu^2, and
    to O(tau_line).)
    '''
    from .toy_transport import brown_mclean
    a, c = geom.scale[0], geom.scale[2]
    q, tau0, gamma = brown_mclean(Ellipsoid(a, c, geom.tau_line / (2.0 * geom.band)), angle)
    return geom.w2 * q
