#!/usr/bin/env python
'''
Polarized electron-scattering transport through the wind of a Sirocco model
(Step 4 of the polarization project).


Command line usage (if any):

    None; this is a module.  step4_regress.py runs it on Sirocco models.

Description:

    The electron density is read, cell by cell, from the master tables that
    windsave2table writes (root.master.txt, or root.0.master.txt,
    root.1.master.txt, ... for several domains), and photons are traced
    through Sirocco's own grid: cylindrical (x, z), polar (r, theta) or
    spherical (r).  The wind is axisymmetric and mirror symmetric about the
    disk plane, as in Sirocco.  Rays are walked from one cell boundary to the
    next, so the optical depth of the piecewise-constant electron density is
    integrated exactly.

    Only electron scattering (Thomson) is included, so the problem is grey
    and only the relative luminosities of the sources matter.  These are
    taken from the last "xdefine_phot:   f_tot ..." line in Sirocco's
    standard output, i.e. the luminosities in the band of the spectral
    cycles.  (The "xdefine_phot: lum_tot ..." line gives bolometric
    luminosities.)  The sources are

    * the central object, the boundary layer and the AGN, emitted from the
      surface of a sphere of the central-object radius;
    * the disk, emitted from both faces of a flat disk with the radial
      distribution of a standard steady-state disk.

    * the wind, from the cells in proportion to their emission in the band
      (windsave2table's heat table: lum_tot for simple atoms, kpkt_lum +
      matom_lum for macro atoms), uniformly within each cell, isotropically
      and unpolarized.

    Surfaces emit with Eddington limb darkening, p(mu) ~ mu (2 + 3 mu), the
    same law extract.c uses.  The disk and the central object are opaque.

Primary routines:

    SiroccoModel, run, brown_mclean_grid

Notes:

    A vertically extended disk is treated as flat.  Cells that are only
    partly in the wind are given their electron density and emission over
    the whole cell.  Wind photons are emitted unpolarized and isotropically
    in the rest frame of the observer; the polarization that line trapping
    can give to line emission (Step 3b) is not included.

History:

    261010 ksl Coding begun, with Claude Code
    261011 ksl Wind emission added; source weights from the band (f_) luminosities

'''

import os
from glob import glob

import numpy as np

from . import photon
from .stokes import normalize, dot
from .es_scatter import scatter_es, peel_es
from .toy_transport import Spectrum

SIGMA_T = 6.6524587e-25

STAT_INWIND = 0
STAT_ESCAPE = 1
STAT_STAR = 2
STAT_DISK = 3
STAT_INTERACT = 4

EPS_REL = 1e-10       # nudge across a boundary, relative to the size of the model
MAX_STEPS = 100000


# ------------------------------------------------------------------ ray-surface distances

def _smallest_positive(r1, r2, eps):
    r1 = np.where(r1 > eps, r1, np.inf)
    r2 = np.where(r2 > eps, r2, np.inf)
    return np.minimum(r1, r2)


def _quadratic(a, b, c, eps):
    '''Smallest root > eps of a s^2 + b s + c = 0 (inf if none).'''
    disc = b * b - 4 * a * c
    ok = disc >= 0
    sq = np.sqrt(np.where(ok, disc, 0.0))
    with np.errstate(divide='ignore', invalid='ignore'):
        # numerically stable roots
        qq = -0.5 * (b + np.where(b >= 0, sq, -sq))
        r1 = np.where(qq != 0, c / qq, np.inf)
        r2 = np.where(a != 0, qq / a, np.where(b != 0, -c / b, np.inf))
    out = _smallest_positive(r1, r2, eps)
    return np.where(ok, out, np.inf)


def ds_to_plane(x, n, z0, eps):
    with np.errstate(divide='ignore', invalid='ignore'):
        s = (z0 - x[:, 2]) / n[:, 2]
    return np.where(s > eps, s, np.inf)


def ds_to_cylinder(x, n, R, eps):
    a = n[:, 0] ** 2 + n[:, 1] ** 2
    b = 2 * (x[:, 0] * n[:, 0] + x[:, 1] * n[:, 1])
    c = x[:, 0] ** 2 + x[:, 1] ** 2 - R * R
    return _quadratic(a, b, c, eps)


def ds_to_sphere(x, n, R, eps):
    return _quadratic(np.ones(len(x)), 2 * dot(x, n), dot(x, x) - R * R, eps)


def ds_to_cone(x, n, theta, eps):
    '''
    Distance to the double cone at polar angle theta (radians, scalar or one
    per ray) from the z axis; theta = pi/2 is the plane z = 0.  The cone is
    written as cos^2(theta) rho^2 - sin^2(theta) z^2 = 0, which stays finite
    at theta = pi/2.
    '''
    c2 = np.cos(theta) ** 2
    s2 = np.sin(theta) ** 2
    a = c2 * (n[:, 0] ** 2 + n[:, 1] ** 2) - s2 * n[:, 2] ** 2
    b = 2 * (c2 * (x[:, 0] * n[:, 0] + x[:, 1] * n[:, 1]) - s2 * x[:, 2] * n[:, 2])
    c = c2 * (x[:, 0] ** 2 + x[:, 1] ** 2) - s2 * x[:, 2] ** 2
    return _quadratic(a, b, c, eps)


# ------------------------------------------------------------------ domains

class Domain:
    '''
    One Sirocco domain: cell edges and the electron density per cell.

    coord is 'spherical' (edges1 = r), 'cylindrical' (edges1 = x, edges2 = z)
    or 'polar' (edges1 = r, edges2 = theta in radians).  The last edge in each
    direction is the inner edge of Sirocco's guard cells, beyond which the
    domain has no wind.
    '''

    def __init__(self, coord, edges1, edges2, ne, emis=None):
        self.coord = coord
        self.e1 = np.asarray(edges1, dtype=float)
        self.e2 = None if edges2 is None else np.asarray(edges2, dtype=float)
        self.ne = np.asarray(ne, dtype=float)
        self.emis = np.zeros_like(self.ne) if emis is None else np.asarray(emis, dtype=float)

    @classmethod
    def from_master(cls, fname, heat=None, simple=True):
        '''
        Read a master table and, if heat is given, the wind emission in the
        band of the spectral cycles from the matching heat table: lum_tot if
        simple (simple-atom emission is used by the model), plus kpkt_lum and
        matom_lum when the table has them.
        '''
        dm = cls._from_master(fname)
        if heat is not None:
            dm.emis = read_emission(heat, dm.ne.shape, simple)
        return dm

    @classmethod
    def _from_master(cls, fname):
        from astropy.io import ascii
        t = ascii.read(fname)
        ne = np.where(t['inwind'] >= 0, t['ne'], 0.0)
        if 'theta' in t.colnames:
            r = np.unique(np.asarray(t['r']))
            th = np.unique(np.asarray(t['theta']))
            th_rad = np.radians(th) if th.max() > 2 * np.pi else th
            grid = np.zeros((len(r), len(th)))
            grid[np.searchsorted(r, t['r']), np.searchsorted(th, t['theta'])] = ne
            return cls('polar', r, th_rad, grid)
        if 'z' in t.colnames:
            x = np.unique(np.asarray(t['x']))
            z = np.unique(np.asarray(t['z']))
            grid = np.zeros((len(x), len(z)))
            grid[np.searchsorted(x, t['x']), np.searchsorted(z, t['z'])] = ne
            return cls('cylindrical', x, z, grid)
        r = np.unique(np.asarray(t['r']))
        grid = np.zeros(len(r))
        grid[np.searchsorted(r, t['r'])] = ne
        return cls('spherical', r, None, grid)

    def _coords(self, x):
        if self.coord == 'cylindrical':
            return np.hypot(x[:, 0], x[:, 1]), np.abs(x[:, 2])
        r = np.linalg.norm(x, axis=1)
        if self.coord == 'spherical':
            return r, None
        with np.errstate(invalid='ignore', divide='ignore'):
            th = np.arccos(np.clip(np.abs(x[:, 2]) / np.maximum(r, 1e-300), 0, 1))
        return r, th

    def _index(self, x):
        c1, c2 = self._coords(x)
        i = np.searchsorted(self.e1, c1, side='right') - 1
        j = None if c2 is None else np.searchsorted(self.e2, c2, side='right') - 1
        if self.coord == 'polar' and self.e2[-1] >= 0.5 * np.pi - 1e-9:
            # theta = 90 deg exactly (a ray in the disk plane) belongs to the
            # last cell above the plane, not to the guard cells beyond it
            j = np.where(j >= len(self.e2) - 1, len(self.e2) - 2, j)
        return i, j

    def ne_at(self, x):
        i, j = self._index(x)
        n1 = len(self.e1)
        if j is None:
            ok = (i >= 0) & (i < n1)
            return np.where(ok, self.ne[np.clip(i, 0, n1 - 1)], 0.0)
        n2 = len(self.e2)
        ok = (i >= 0) & (i < n1) & (j >= 0) & (j < n2)
        return np.where(ok, self.ne[np.clip(i, 0, n1 - 1), np.clip(j, 0, n2 - 1)], 0.0)

    def ds_to_boundary(self, x, n, eps):
        '''Distance to the nearest boundary of the cell containing x.'''
        i, j = self._index(x)
        s = np.full(len(x), np.inf)
        n1 = len(self.e1)
        for k in (i, i + 1):
            ok = (k >= 0) & (k < n1)
            if ok.any():
                R = self.e1[np.clip(k, 0, n1 - 1)]
                if self.coord == 'cylindrical':
                    d = ds_to_cylinder(x, n, R, eps)
                else:
                    d = ds_to_sphere(x, n, R, eps)
                s = np.minimum(s, np.where(ok, d, np.inf))
        if j is not None:
            n2 = len(self.e2)
            for k in (j, j + 1):
                ok = (k >= 0) & (k < n2)
                if ok.any():
                    v = self.e2[np.clip(k, 0, n2 - 1)]
                    if self.coord == 'cylindrical':
                        d = np.minimum(ds_to_plane(x, n, v, eps), ds_to_plane(x, n, -v, eps))
                    else:
                        d = ds_to_cone(x, n, v, eps)
                    s = np.minimum(s, np.where(ok, d, np.inf))
        return s

    def random_points(self, cell, rng):
        '''Points distributed uniformly in volume within the given cells (flat indices), both hemispheres.'''
        m = len(cell)
        if self.coord == 'spherical':
            r0, r1 = self.e1[cell], self.e1[cell + 1]
            r = np.cbrt(r0 ** 3 + rng.random(m) * (r1 ** 3 - r0 ** 3))
            return normalize(rng.normal(size=(m, 3))) * r[:, None]
        i, j = np.unravel_index(cell, self.ne.shape)
        a0, a1 = self.e1[i], self.e1[i + 1]
        b0, b1 = self.e2[j], self.e2[j + 1]
        phi = rng.uniform(0, 2 * np.pi, m)
        side = np.where(rng.random(m) < 0.5, 1.0, -1.0)
        if self.coord == 'cylindrical':
            rho = np.sqrt(a0 ** 2 + rng.random(m) * (a1 ** 2 - a0 ** 2))
            z = b0 + rng.random(m) * (b1 - b0)
        else:
            r = np.cbrt(a0 ** 3 + rng.random(m) * (a1 ** 3 - a0 ** 3))
            mu = np.cos(b1) + rng.random(m) * (np.cos(b0) - np.cos(b1))
            rho, z = r * np.sqrt(1 - mu * mu), r * mu
        return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), side * z])

    def extent(self):
        '''A radius that encloses the whole grid.'''
        if self.coord == 'cylindrical':
            return np.hypot(self.e1[-1], self.e2[-1])
        return self.e1[-1]

    def cells(self):
        '''Cell centres (r, mu), volumes (both hemispheres) and ne, for cells with outer edges.'''
        if self.coord == 'spherical':
            r0, r1 = self.e1[:-1], self.e1[1:]
            vol = 4 / 3 * np.pi * (r1 ** 3 - r0 ** 3)
            rc = 0.5 * (r0 + r1)
            # a spherical shell: <mu^2> = 1/3
            return rc, np.full(len(rc), np.sqrt(1 / 3)), vol, self.ne[:-1]
        a0, a1 = np.meshgrid(self.e1[:-1], self.e2[:-1], indexing='ij')
        b0, b1 = np.meshgrid(self.e1[1:], self.e2[1:], indexing='ij')
        if self.coord == 'cylindrical':
            vol = 2 * np.pi * (b0 ** 2 - a0 ** 2) * (b1 - a1)
            xc, zc = 0.5 * (a0 + b0), 0.5 * (a1 + b1)
            rc = np.hypot(xc, zc)
            mu = zc / rc
        else:
            vol = 2 * 2 * np.pi / 3 * (b0 ** 3 - a0 ** 3) * (np.cos(a1) - np.cos(b1))
            rc = 0.5 * (a0 + b0)
            mu = np.cos(0.5 * (a1 + b1))
        return rc.ravel(), mu.ravel(), vol.ravel(), self.ne[:-1, :-1].ravel()


# ------------------------------------------------------------------ the model

def read_pf(fname):
    '''Read a parameter file into {keyword without units: [values]}.'''
    out = {}
    with open(fname) as f:
        for line in f:
            w = line.split()
            if len(w) < 2 or w[0].startswith('#'):
                continue
            key = w[0].split('(')[0]
            out.setdefault(key, []).append(w[1])
    return out


def read_luminosities(stdout, band=True):
    '''
    The last xdefine_phot luminosities in Sirocco's output: in the band of the
    spectral cycles (f_tot, f_star, f_bl, f_bh, f_disk, f_wind, f_matom,
    f_kpkt) if band, otherwise bolometric (lum_tot, lum_star, ...).
    '''
    tag, prefix = ('xdefine_phot:   f_tot', 'f_') if band else ('xdefine_phot: lum_tot', 'lum_')
    lum = None
    with open(stdout) as f:
        for line in f:
            if tag in line:
                w = line.split()
                lum = {w[k]: float(w[k + 1]) for k in range(len(w) - 1) if w[k].startswith(prefix)}
    return lum


def read_emission(fname, shape, simple=True):
    '''Wind emission per cell in the band of the spectral cycles, from a windsave2table heat table.'''
    from astropy.io import ascii
    t = ascii.read(fname)
    inwind = np.asarray(t['inwind']) >= 0
    e = np.zeros(len(t))
    if simple:
        e += np.asarray(t['lum_tot'])
    for col in ('kpkt_lum', 'matom_lum'):
        if col in t.colnames:
            e += np.asarray(t[col])
    e = np.where(inwind, e, 0.0)
    grid = np.zeros(shape)
    if len(shape) == 1:
        grid[np.asarray(t['i'])] = e
    else:
        grid[np.asarray(t['i']), np.asarray(t['j'])] = e
    # the guard cells (no outer edge) cannot emit
    if len(shape) == 1:
        grid[-1] = 0.0
    else:
        grid[-1, :] = grid[:, -1] = 0.0
    return grid


class SiroccoModel:
    '''
    The geometry, electron density and sources of a Sirocco run, read from
    root.pf (or root.out.pf), the windsave2table master tables, and the
    standard output (root.stdout.txt, for the source luminosities).
    '''

    def __init__(self, root, directory='.', point_source=False, wind=True):
        d = directory
        pf_name = os.path.join(d, root + '.out.pf')
        if not os.path.exists(pf_name):
            pf_name = os.path.join(d, root + '.pf')
        self.pf = read_pf(pf_name)
        lum = read_luminosities(os.path.join(d, root + '.stdout.txt'))
        self.lum = lum
        masters = sorted(glob(os.path.join(d, root + '.[0-9]*.master.txt')))
        if not masters:
            masters = [os.path.join(d, root + '.master.txt')]
        use_wind = wind and not point_source
        # simple-atom wind emission is in f_wind; macro-atom emission in f_kpkt and f_matom
        simple = lum.get('f_wind', 0) > 0
        self.domains = [Domain.from_master(m, m.replace('.master.txt', '.heat.txt') if use_wind else None, simple)
                        for m in masters]

        def get(key, default=None):
            v = self.pf.get(key)
            return default if v is None else v[0]

        self.r_star = float(get('Central_object.radius', 0.0))
        disk_type = get('Disk.type', 'none')
        self.disk = disk_type != 'none'
        self.r_disk = float(get('Disk.radmax', 0.0)) if self.disk else 0.0
        radmax = [float(v) for v in self.pf.get('Wind.radmax', [])]
        self.r_max = max([dm.extent() for dm in self.domains] + radmax)
        self.eps = EPS_REL * self.r_max
        nobs = int(get('Spectrum.no_observers', 0))
        # .out.pf files can repeat the observer block; take the first nobs
        self.angles = [float(v) for v in self.pf.get('Spectrum.angle', [])][:nobs]
        self.phase = float(get('Spectrum.orbit_phase', 0.5))

        self.point_source = point_source
        if point_source:
            # an isotropic point source, no occultation: the Brown & McLean case
            self.sources = {'point': 1.0}
            self.r_star = 0.0
            self.disk = False
        else:
            central = lum.get('f_star', 0) + lum.get('f_bl', 0) + lum.get('f_bh', 0)
            disk = lum.get('f_disk', 0) if self.disk else 0.0
            wnd = sum(dm.emis.sum() for dm in self.domains)
            tot = central + disk + wnd
            self.sources = {'central': central / tot, 'disk': disk / tot, 'wind': wnd / tot}
        self.disk_r_in = max(self.r_star, 1e-6 * self.r_max)

    @classmethod
    def from_domains(cls, domains, r_max, angles=(30, 60, 90), phase=0.0):
        '''
        A model with an isotropic point source at the center and no disk or
        central object, built directly from Domain objects (for tests).
        '''
        m = cls.__new__(cls)
        m.pf = {}
        m.domains = list(domains)
        m.r_star = 0.0
        m.disk = False
        m.r_disk = 0.0
        m.r_max = float(r_max)
        m.eps = EPS_REL * m.r_max
        m.angles = list(angles)
        m.phase = phase
        m.lum = None
        m.point_source = True
        m.sources = {'point': 1.0}
        m.disk_r_in = 0.0
        return m

    # -------------------------------------------------- opacity and walking

    def kappa(self, x):
        ne = np.zeros(len(x))
        for dm in self.domains:
            ne = ne + dm.ne_at(x)
        return SIGMA_T * ne

    def walk(self, x, n, tau_target):
        '''
        Walk rays from x along n until the accumulated electron-scattering
        optical depth reaches tau_target, or the ray escapes beyond r_max, hits
        the central object, or hits the disk.  Returns (x_end, status, tau).
        '''
        x = x.copy()
        nray = len(x)
        tau = np.zeros(nray)
        status = np.full(nray, STAT_INWIND)
        todo = np.arange(nray)
        steps = 0
        while todo.size and steps < MAX_STEPS:
            xt, nt = x[todo], n[todo]
            s = np.full(todo.size, np.inf)
            for dm in self.domains:
                s = np.minimum(s, dm.ds_to_boundary(xt, nt, self.eps))
            s_out = ds_to_sphere(xt, nt, self.r_max, self.eps)
            s_star = ds_to_sphere(xt, nt, self.r_star, self.eps) if self.r_star > 0 else np.full(todo.size, np.inf)
            s_plane = ds_to_plane(xt, nt, 0.0, self.eps)
            if self.disk:
                xp = xt + np.where(np.isfinite(s_plane), s_plane, 0.0)[:, None] * nt
                rp = np.hypot(xp[:, 0], xp[:, 1])
                s_disk = np.where((rp >= self.disk_r_in) & (rp <= self.r_disk), s_plane, np.inf)
            else:
                s_disk = np.full(todo.size, np.inf)
            outside = dot(xt, xt) >= self.r_max ** 2
            s_out = np.where(outside, 0.0, s_out)       # already beyond r_max: escape now
            s_stop = np.minimum.reduce([s_out, s_star, s_disk])
            step = np.minimum(np.minimum(s, s_stop), s_plane)
            step = np.where(np.isfinite(step), step, 0.0)
            s_out = np.where(np.isfinite(s_out), s_out, 0.0)
            s_stop = np.minimum.reduce([s_out, s_star, s_disk])
            k = self.kappa(xt + 0.5 * step[:, None] * nt)
            dtau = k * step
            hit = tau[todo] + dtau >= tau_target[todo]
            # interaction inside this step
            with np.errstate(divide='ignore', invalid='ignore'):
                s_int = np.where(hit, (tau_target[todo] - tau[todo]) / k, 0.0)
            new_x = xt + np.where(hit, s_int, step)[:, None] * nt
            tau[todo] = np.where(hit, tau_target[todo], tau[todo] + dtau)
            st = np.where(hit, STAT_INTERACT, STAT_INWIND)
            done_stop = ~hit & (step >= s_stop)
            st = np.where(done_stop & (s_out <= step), STAT_ESCAPE, st)
            st = np.where(done_stop & (s_star <= step) & (s_out > step), STAT_STAR, st)
            st = np.where(done_stop & (s_disk <= step) & (s_out > step) & (s_star > step), STAT_DISK, st)
            moving = st == STAT_INWIND
            new_x[moving] += self.eps * nt[moving]          # step into the next cell
            x[todo] = new_x
            status[todo] = st
            todo = todo[moving]
            steps += 1
        status[todo] = STAT_ESCAPE
        return x, status, tau

    def tau_to_observer(self, x, n_obs):
        '''Electron-scattering optical depth to escape along n_obs (inf if blocked).'''
        n = np.broadcast_to(n_obs, x.shape).copy()
        _, status, tau = self.walk(x, n, np.full(len(x), np.inf))
        return np.where(status == STAT_ESCAPE, tau, np.inf)

    # -------------------------------------------------- sources

    def generate(self, nphot, rng):
        '''Photons from the sources, with total weight 1 and the sources' weights.'''
        p = photon.make_photons(nphot)
        names = list(self.sources)
        frac = np.array([self.sources[k] for k in names])
        which = rng.choice(len(names), size=nphot, p=frac / frac.sum())
        p['w'] = p['w_orig'] = 1.0 / nphot
        p['origin'] = -1
        for k, name in enumerate(names):
            idx = np.where(which == k)[0]
            m = idx.size
            if m == 0:
                continue
            if name == 'point':
                p['x'][idx] = 0.0
                p['lmn'][idx] = normalize(rng.normal(size=(m, 3)))
                p['origin'][idx] = photon.PTYPE_STAR
            elif name == 'wind':
                p['x'][idx] = self._wind_positions(rng, m)
                p['lmn'][idx] = normalize(rng.normal(size=(m, 3)))
                p['origin'][idx] = photon.PTYPE_WIND
            elif name == 'central':
                rhat = normalize(rng.normal(size=(m, 3)))
                p['x'][idx] = rhat * (self.r_star + self.eps)
                p['lmn'][idx] = _limb_darkened(rng, rhat)
                p['origin'][idx] = photon.PTYPE_STAR
            else:
                r = self._disk_radii(rng, m)
                phi = rng.uniform(0, 2 * np.pi, m)
                side = np.where(rng.random(m) < 0.5, 1.0, -1.0)
                normal = np.zeros((m, 3))
                normal[:, 2] = side
                p['x'][idx] = np.column_stack([r * np.cos(phi), r * np.sin(phi), side * self.eps])
                p['lmn'][idx] = _limb_darkened(rng, normal)
                p['origin'][idx] = photon.PTYPE_DISK
        photon.set_unpolarized(p)
        return p

    def _wind_positions(self, rng, m):
        '''Positions of wind photons: cells chosen by their emission, uniform within a cell.'''
        emis = np.concatenate([dm.emis.ravel() for dm in self.domains])
        bounds = np.cumsum([0] + [dm.emis.size for dm in self.domains])
        cell = rng.choice(emis.size, size=m, p=emis / emis.sum())
        x = np.empty((m, 3))
        for k, dm in enumerate(self.domains):
            sel = (cell >= bounds[k]) & (cell < bounds[k + 1])
            if sel.any():
                x[sel] = dm.random_points(cell[sel] - bounds[k], rng)
        return x

    def _disk_radii(self, rng, m):
        '''Radii from a standard disk, dL/dr ~ r^-2 (1 - sqrt(r_in/r)).'''
        r_in = self.disk_r_in
        out = np.empty(0)
        while out.size < m:
            k = 2 * (m - out.size) + 10
            # draw from dL/dr ~ r^-2 between r_in and r_disk, then reject
            u = rng.random(k)
            r = 1.0 / (1.0 / r_in - u * (1.0 / r_in - 1.0 / self.r_disk))
            keep = rng.random(k) < 1.0 - np.sqrt(r_in / r)
            out = np.concatenate([out, r[keep]])
        return out[:m]

    def emission_factor(self, p, n_obs):
        '''Peel-off factor at birth relative to isotropic emission: mu (2 + 3 mu), or 1.'''
        f = np.ones(len(p))
        x = p['x']
        cen = p['origin'] == photon.PTYPE_STAR
        dsk = p['origin'] == photon.PTYPE_DISK
        if not self.point_source:
            mu = dot(normalize(x[cen]), n_obs)
            f[cen] = np.where(mu > 0, mu * (2 + 3 * mu), 0.0)
            mu = np.sign(x[dsk, 2]) * n_obs[2]
            f[dsk] = np.where(mu > 0, mu * (2 + 3 * mu), 0.0)
        return f


def _limb_darkened(rng, normal):
    '''Directions about the normals with p(mu) ~ mu (2 + 3 mu).'''
    m = len(normal)
    mu = np.empty(m)
    todo = np.arange(m)
    while todo.size:
        t = rng.random(todo.size)
        ok = rng.random(todo.size) * 5.0 < t * (2 + 3 * t)
        mu[todo[ok]] = t[ok]
        todo = todo[~ok]
    phi = rng.uniform(0, 2 * np.pi, m)
    a = normalize(np.cross(normal, np.where(np.abs(normal[:, 2:3]) < 0.9, [0, 0, 1.0], [1.0, 0, 0])))
    b = np.cross(normal, a)
    st = np.sqrt(1 - mu * mu)
    return normalize(mu[:, None] * normal + st[:, None] * (np.cos(phi)[:, None] * a + np.sin(phi)[:, None] * b))


# ------------------------------------------------------------------ transport

def extract(model, p, spectra, scattered):
    for spec in spectra:
        tau = model.tau_to_observer(p['x'], spec.lmn)
        if scattered:
            f, q, u = peel_es(p, spec.lmn, spec.north)
        else:
            f = model.emission_factor(p, spec.lmn)
            q = u = np.zeros(len(p))
        spec.tally(p['w'] * f * np.exp(-tau), q, u, scattered)


def trans_phot(model, p, rng, spectra, maxscat=1000):
    extract(model, p, spectra, scattered=False)
    active = np.arange(len(p))
    while active.size:
        pa = p[active]
        tau_target = -np.log(rng.random(active.size))
        x, status, _ = model.walk(pa['x'], pa['lmn'], tau_target)
        pa['x'] = x
        inter = status == STAT_INTERACT
        pa['istat'] = np.where(inter, photon.P_INWIND,
                               np.where(status == STAT_ESCAPE, photon.P_ESCAPE,
                                        np.where(status == STAT_STAR, photon.P_HIT_STAR, photon.P_HIT_DISK)))
        if inter.any():
            ps = pa[inter]
            ps['nscat'] += 1
            ps['nres'] = photon.NRES_ES
            extract(model, ps, spectra, scattered=True)
            scatter_es(ps, rng)
            ps['istat'][ps['nscat'] >= maxscat] = photon.P_TOO_MANY_SCATTERS
            pa[inter] = ps
        p[active] = pa
        active = active[pa['istat'] == photon.P_INWIND]
    return p


def run(model, nphot, rng, angles=None, nbatches=20, max_batch=100000):
    '''
    Run nphot photons through the model, peeling off toward each inclination
    (default: the model's Spectrum.angle values, at its orbital phase).
    Returns the list of toy_transport.Spectrum objects.
    '''
    angles = model.angles if angles is None else angles
    spectra = [Spectrum(a, model.phase) for a in angles]
    nb = max(nbatches, int(np.ceil(nphot / max_batch)))
    sizes = np.full(nb, nphot // nb)
    sizes[:nphot % nb] += 1
    for size in sizes:
        p = model.generate(size, rng)
        p['w'] *= size / nphot
        trans_phot(model, p, rng, spectra)
        for spec in spectra:
            spec.end_batch()
    return spectra


def brown_mclean_grid(model, angle):
    '''
    The optically thin, point-source prediction of Brown & McLean (1977)
    evaluated directly from the grid cells:

        tau0  = (3/16) sigma_T sum n_e V / (2 pi r^2)
        gamma = sum n_e V mu^2 / (2 pi r^2) / sum n_e V / (2 pi r^2)
        q     = tau0 (1 - 3 gamma) sin^2 i

    (integral n dr dmu = integral n dV / (2 pi r^2) for an axisymmetric
    density).  Returns (q, tau0, gamma).
    '''
    S = 0.0
    S2 = 0.0
    for dm in model.domains:
        rc, mu, vol, ne = dm.cells()
        w = ne * vol / (2 * np.pi * rc ** 2)
        S += w.sum()
        S2 += (w * mu * mu).sum()
    tau0 = 3.0 / 16.0 * SIGMA_T * S
    gamma = S2 / S
    return tau0 * (1 - 3 * gamma) * np.sin(np.radians(angle)) ** 2, tau0, gamma
