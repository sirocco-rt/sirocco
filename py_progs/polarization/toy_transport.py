#!/usr/bin/env python
'''
Toy Monte Carlo transport of polarized photons through pure electron-
scattering media with analytic geometries (Step 2 of the polarization
project).


Command line usage (if any):

    None; this is a module.  step2_validate.py runs it.

Description:

    The structure follows Sirocco's spectral cycles: photons are generated,
    moved to their next interaction (cf. translate/trans_phot), peeled off
    toward each observer at every scatter (cf. extract/extract_one), and
    then scattered (cf. scatter).  Escaping photons are also binned by
    direction, giving a live-or-die spectrum to check the peel-off against.

    Geometries (the opacity kappa = n_e sigma_T is uniform inside each):

    * Ellipsoid -- an axisymmetric ellipsoid with equatorial semi-axis a and
      polar semi-axis c, around a central, isotropic, unpolarized point
      source.  a = c is a uniform sphere.
    * Slab -- a plane-parallel slab, measured in optical depth, with
      -tau_max < z < 0 and the surface at z = 0 facing +z.  Photons are
      injected upward at the base with p(mu) ~ mu, and photons that leave
      through the base are lost.  For tau_max >> 1 the emergent radiation is
      that of Chandrasekhar's Milne problem.

Primary routines:

    run, trans_phot, extract_one, brown_mclean, and the classes
    Ellipsoid, Slab, Spectrum, LiveOrDie

Notes:

    The tallies are energy per unit solid angle (as in Sirocco's extracted
    spectra), normalized so that an isotropic, unscattered source of total
    weight 1 gives I = 1 in every direction.  For the slab this quantity is
    mu times the specific intensity at the surface.

History:

    261010 ksl Coding begun, with Claude Code

'''

import numpy as np

from . import photon
from .photon import P_INWIND, P_ESCAPE, P_ABSORB, P_TOO_MANY_SCATTERS, NRES_ES
from .stokes import normalize, dot, observer_direction, sky_north, to_sky
from .es_scatter import scatter_es, peel_es

ZAXIS = np.array([0.0, 0.0, 1.0])


class Ellipsoid:
    '''Uniform electron-scattering ellipsoid around a point source at the origin.'''

    def __init__(self, a, c, kappa):
        self.scale = np.array([a, a, c], dtype=float)
        self.kappa = float(kappa)

    def tau_radial(self, mu):
        '''Optical depth from the center to the surface along direction cosine mu.'''
        a, c = self.scale[0], self.scale[2]
        return self.kappa / np.sqrt((1.0 - mu * mu) / a ** 2 + mu * mu / c ** 2)

    def ds_to_edge(self, x, n):
        '''
        Distance from x (inside) to the surface along n, and the status a
        photon has when it reaches it (always P_ESCAPE here).
        '''
        X = x / self.scale
        D = n / self.scale
        A = dot(D, D)
        B = 2.0 * dot(X, D)
        Cc = dot(X, X) - 1.0
        s = (-B + np.sqrt(np.maximum(B * B - 4.0 * A * Cc, 0.0))) / (2.0 * A)
        return np.maximum(s, 0.0), np.full(len(s), P_ESCAPE)

    def generate(self, nphot, rng):
        '''Isotropic, unpolarized photons from the origin, total weight 1.'''
        p = photon.make_photons(nphot)
        p['lmn'] = normalize(rng.normal(size=(nphot, 3)))
        p['x'] = 0.0
        p['w'] = p['w_orig'] = 1.0 / nphot
        p['origin'] = p['origin_orig'] = photon.PTYPE_STAR
        photon.set_unpolarized(p)
        return p

    def emission_factor(self, p, n_obs):
        '''Peel-off factor for direct emission relative to isotropic: 1.'''
        return np.ones(len(p))


class Slab:
    '''Pure-scattering slab in optical-depth units (kappa = 1), lit from below.'''

    def __init__(self, tau_max):
        self.tau_max = float(tau_max)
        self.kappa = 1.0

    def ds_to_edge(self, x, n):
        '''
        Distance to the top (status P_ESCAPE) or the base (status P_ABSORB,
        i.e. lost) along n.
        '''
        z = x[:, 2]
        nz = n[:, 2]
        with np.errstate(divide='ignore', invalid='ignore'):
            s = np.where(nz > 0, -z / nz, np.where(nz < 0, (z + self.tau_max) / -nz, np.inf))
        istat = np.where(nz > 0, P_ESCAPE, P_ABSORB)
        return np.maximum(s, 0.0), istat

    def generate(self, nphot, rng):
        '''Unpolarized photons injected upward at the base, p(mu) ~ mu.'''
        p = photon.make_photons(nphot)
        mu = np.sqrt(rng.random(nphot))
        phi = rng.uniform(0.0, 2.0 * np.pi, nphot)
        st = np.sqrt(1.0 - mu * mu)
        p['lmn'] = np.column_stack([st * np.cos(phi), st * np.sin(phi), mu])
        p['x'] = 0.0
        p['x'][:, 2] = -self.tau_max
        p['w'] = p['w_orig'] = 1.0 / nphot
        p['origin'] = p['origin_orig'] = photon.PTYPE_DISK
        photon.set_unpolarized(p)
        return p

    def emission_factor(self, p, n_obs):
        '''Peel-off factor for direct emission relative to isotropic: 4 mu.'''
        return 4.0 * np.maximum(n_obs[2], 0.0) * np.ones(len(p))


class Spectrum:
    '''
    Peel-off tallies for one observer, the analog of one xxspec[].  Because
    the toy problems are grey, there is a single frequency bin.

    Errors come from the scatter between batches of photons.  The
    contributions of one photon at its successive scatters are correlated,
    so a sum of squared contributions would underestimate them.
    '''

    def __init__(self, angle, phase=0.0):
        self.angle = angle
        self.phase = phase
        self.lmn = observer_direction(angle, phase)
        self.north = sky_north(angle, phase)
        self.I = self.Q = self.U = 0.0
        self.I_scat = 0.0      # scattered light only (cf. f_wind)
        self.batches = []      # (I, Q, U) of each completed batch
        self._b = np.zeros(3)

    def tally(self, wt, q, u, scattered):
        dI, dQ, dU = wt.sum(), (wt * q).sum(), (wt * u).sum()
        self.I += dI
        self.Q += dQ
        self.U += dU
        if scattered:
            self.I_scat += dI
        self._b += (dI, dQ, dU)

    def end_batch(self):
        self.batches.append(self._b.copy())
        self._b[:] = 0.0

    def errors(self):
        '''
        Return the 1-sigma errors of (I, q, u) from the batch-to-batch
        scatter (NaN if there are fewer than 2 batches).
        '''
        b = np.array(self.batches)
        nb = len(b)
        if nb < 2:
            return np.nan, np.nan, np.nan
        I_b = b[:, 0] * nb               # each batch scaled to the full run
        sI = np.std(I_b, ddof=1) / np.sqrt(nb)
        # linearized error of the ratio estimators Q/I and U/I
        q, u = self.Q / self.I, self.U / self.I
        sq = np.std(b[:, 1] * nb - q * I_b, ddof=1) / np.sqrt(nb) / self.I
        su = np.std(b[:, 2] * nb - u * I_b, ddof=1) / np.sqrt(nb) / self.I
        return sI, sq, su

    def result(self):
        '''Return q = Q/I, u = U/I and their 1-sigma errors.'''
        _, sq, su = self.errors()
        return self.Q / self.I, self.U / self.I, sq, su


class LiveOrDie:
    '''
    Escaping photons binned in mu = cos(inclination), averaged over
    azimuth (appropriate for axisymmetric problems).  Each photon's Stokes
    parameters are referred to its own sky North, the projected z axis.
    '''

    def __init__(self, nbins=20):
        self.edges = np.linspace(-1.0, 1.0, nbins + 1)
        self.I = np.zeros(nbins)
        self.Q = np.zeros(nbins)
        self.U = np.zeros(nbins)
        self.n = np.zeros(nbins, dtype=int)
        self.wq2 = np.zeros(nbins)
        self.wu2 = np.zeros(nbins)

    def tally(self, p):
        n = p['lmn']
        north = n[:, 2:3] * n - ZAXIS                    # -(z - (z.n) n): sign is irrelevant
        bad = np.linalg.norm(north, axis=1) < 1e-9
        north[bad] = [1.0, 0.0, 0.0]
        north = normalize(north)
        q, u = to_sky(p['stokes'][:, 0], p['stokes'][:, 1], p['pol_ref'], n, north)
        k = np.clip(np.searchsorted(self.edges, n[:, 2]) - 1, 0, len(self.I) - 1)
        np.add.at(self.I, k, p['w'])
        np.add.at(self.Q, k, p['w'] * q)
        np.add.at(self.U, k, p['w'] * u)
        np.add.at(self.n, k, 1)
        np.add.at(self.wq2, k, (p['w'] * q) ** 2)
        np.add.at(self.wu2, k, (p['w'] * u) ** 2)

    def result(self):
        '''
        Return bin centers in mu, I (energy per unit solid angle, in the
        same units as Spectrum.I), q, u, and the 1-sigma errors of q and u.
        '''
        dmu = np.diff(self.edges)
        mu = 0.5 * (self.edges[1:] + self.edges[:-1])
        with np.errstate(divide='ignore', invalid='ignore'):
            return (mu, self.I / (dmu / 2.0), self.Q / self.I, self.U / self.I,
                    np.sqrt(self.wq2) / self.I, np.sqrt(self.wu2) / self.I)


def extract_one(geom, p, spec, scattered):
    '''
    Peel photons p (at their current position, after any scatter has been
    chosen) toward spec: attenuate by exp(-tau) to the edge, and tally only
    photons whose path toward the observer escapes.
    '''
    n_obs = np.broadcast_to(spec.lmn, p['lmn'].shape)
    s, istat = geom.ds_to_edge(p['x'], n_obs)
    tau = geom.kappa * s
    if scattered:
        f, q, u = peel_es(p, spec.lmn, spec.north)
    else:
        f = geom.emission_factor(p, spec.lmn)
        q = u = np.zeros(len(p))
    wt = p['w'] * f * np.exp(-tau) * (istat == P_ESCAPE)
    spec.tally(wt, q, u, scattered)


def trans_phot(geom, p, rng, spectra, live=None, maxscat=100000):
    '''
    Transport photons p through geom until they escape, are lost, or
    exceed maxscat scatters.  At birth and at every scatter, each photon is
    peeled off toward every Spectrum in spectra.  Escaping photons are
    added to live (a LiveOrDie) if one is given.  p is changed in place;
    on return p['istat'] holds each photon's fate.
    '''
    for spec in spectra:
        extract_one(geom, p, spec, scattered=False)

    active = np.arange(len(p))
    while active.size:
        pa = p[active]
        tau_scat = -np.log(rng.random(active.size))
        s_edge, istat = geom.ds_to_edge(pa['x'], pa['lmn'])
        s = tau_scat / geom.kappa

        leave = s >= s_edge
        pa['x'][leave] += s_edge[leave, None] * pa['lmn'][leave]
        pa['istat'][leave] = istat[leave]

        stay = ~leave
        pa['x'][stay] += s[stay, None] * pa['lmn'][stay]
        ps = pa[stay]
        ps['nscat'] += 1
        ps['nres'] = NRES_ES
        for spec in spectra:
            extract_one(geom, ps, spec, scattered=True)
        scatter_es(ps, rng)
        too_many = ps['nscat'] >= maxscat
        ps['istat'][too_many] = P_TOO_MANY_SCATTERS
        pa[stay] = ps

        p[active] = pa
        if live is not None:
            live.tally(pa[pa['istat'] == P_ESCAPE])
        active = active[pa['istat'] == P_INWIND]
    return p


def brown_mclean(geom, angle, nmu=2001):
    '''
    The optically thin, point-source prediction of Brown & McLean (1977)
    for an Ellipsoid:  q = tau0 (1 - 3 gamma) sin^2(i), with
      tau0  = (3/16) integral_{-1}^{1} tau_r(mu) dmu
      gamma = integral tau_r mu^2 dmu / integral tau_r dmu
    where tau_r(mu) is the radial optical depth.  q > 0 (E-vector along the
    projected axis) for a flattened envelope, in the IAU sign convention.
    Returns (q, tau0, gamma).
    '''
    mu, wmu = np.polynomial.legendre.leggauss(64)
    t = geom.tau_radial(mu)
    tau0 = 3.0 / 16.0 * np.sum(wmu * t)
    gamma = np.sum(wmu * t * mu * mu) / np.sum(wmu * t)
    return tau0 * (1.0 - 3.0 * gamma) * np.sin(np.radians(angle)) ** 2, tau0, gamma


def run(geom, nphot, angles, rng, nbatches=20, live_bins=20, phase=0.0, max_batch=200000):
    '''
    Run nphot photons through geom, peeling off toward each inclination in
    angles.  The photons are run in at least nbatches equal batches (more if
    a batch would exceed max_batch photons), which give the errors.
    Returns (spectra, live).
    '''
    spectra = [Spectrum(a, phase) for a in angles]
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
