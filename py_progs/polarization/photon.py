#!/usr/bin/env python
'''
A Python mirror of the Sirocco photon structure, for prototyping polarization.


Command line usage (if any):

    None; this is a module, imported as ``from polarization import photon``

Description:

    The field names, enum values and helper routines follow ``struct photon``
    in ``source/sirocco.h`` and ``source/phot_util.c``, so that code written
    against this module can be ported to C more or less line by line.  Two
    fields are new:

    * ``stokes`` -- (q, u, v) = (Q/I, U/I, V/I); I is the photon weight ``w``
    * ``pol_ref`` -- e1, a unit vector perpendicular to ``lmn`` which defines
      the axis along which Q > 0

    A bundle of photons is a NumPy structured array, the analog of
    ``photmain``.  ``p['lmn']`` is then an (N, 3) array, so routines elsewhere
    in this package operate on whole arrays (or a slice, ``p[i:i+1]``, for a
    single photon).

Primary routines:

    make_photons, init_dummy_phot, set_unpolarized, stuff_phot,
    default_pol_ref

Notes:

    The conventions are those in the Obsidian note
    Sirocco-2026-Polarization-Conventions.

History:

    261010 ksl Coding begun, with Claude Code

'''

import numpy as np

# enum istat_enum
P_INWIND = 0
P_SCAT = 1
P_ESCAPE = 2
P_HIT_STAR = 3
P_TOO_MANY_SCATTERS = 4
P_ERROR = 5
P_ABSORB = 6
P_HIT_DISK = 7
P_SEC = 8
P_ADIABATIC = 9
P_ERROR_MATOM = 10
P_LOFREQ_FF = 11
P_REPOSITION_ERROR = 12

# enum frame
F_LOCAL = 0
F_OBSERVER = 1

# enum origin_enum
PTYPE_STAR = 0
PTYPE_BL = 1
PTYPE_DISK = 2
PTYPE_WIND = 3
PTYPE_AGN = 4
PTYPE_STAR_MATOM = 10
PTYPE_BL_MATOM = 11
PTYPE_DISK_MATOM = 12
PTYPE_WIND_MATOM = 13
PTYPE_AGN_MATOM = 14
PTYPE_DUMMY = -1

# nres values for continuum processes
NRES_ES = -1
NRES_FF = -2
NRES_NOT_SET = -3

photon_dtype = np.dtype([
    ('x', 'f8', 3),
    ('lmn', 'f8', 3),
    ('freq', 'f8'),
    ('freq_orig', 'f8'),
    ('w', 'f8'),
    ('w_orig', 'f8'),
    ('tau', 'f8'),
    ('istat', 'i4'),
    ('frame', 'i4'),
    ('nscat', 'i4'),
    ('nrscat', 'i4'),
    ('nmacro', 'i4'),
    ('nres', 'i4'),
    ('line_res', 'i4'),
    ('nnscat', 'i4'),
    ('grid', 'i4'),
    ('origin', 'i4'),
    ('origin_orig', 'i4'),
    ('np', 'i4'),
    ('path', 'f8'),
    ('ds', 'f8'),
    ('stokes', 'f8', 3),    # new: q, u, v
    ('pol_ref', 'f8', 3),   # new: e1
])


def default_pol_ref(lmn):
    '''
    Return the reference axis used for an unpolarized photon moving along lmn.

    e1 = normalize(z x n), or normalize(x x n) if n is (anti)parallel to z.
    The choice is arbitrary for an unpolarized photon, but it is made
    deterministically so that runs are reproducible.

    lmn may be a single vector (3,) or an array (N, 3).
    '''
    n = np.atleast_2d(lmn)
    e1 = np.cross([0.0, 0.0, 1.0], n)
    norm = np.linalg.norm(e1, axis=1)
    bad = norm < 1e-6
    if np.any(bad):
        e1[bad] = np.cross([1.0, 0.0, 0.0], n[bad])
        norm[bad] = np.linalg.norm(e1[bad], axis=1)
    e1 /= norm[:, None]
    return e1.reshape(np.shape(lmn))


def set_unpolarized(p):
    '''
    Make the photons in p unpolarized: q = u = v = 0, and pol_ref set to
    the default axis for the current direction.  This is what every source
    and every non-electron-scattering re-emission does.
    '''
    p['stokes'] = 0.0
    p['pol_ref'] = default_pol_ref(p['lmn'])


def init_dummy_phot(p):
    '''
    Mirror of init_dummy_phot in phot_util.c: put photons in a defined
    state.  Here the photons are also made unpolarized.
    '''
    p['x'] = 0.0
    p['lmn'] = [0.0, 0.0, 1.0]
    p['freq'] = p['freq_orig'] = 0.0
    p['w'] = p['w_orig'] = 0.0
    p['tau'] = 0.0
    p['istat'] = P_INWIND
    p['frame'] = F_OBSERVER
    p['nscat'] = p['nrscat'] = p['nmacro'] = 0
    p['nres'] = NRES_NOT_SET
    p['line_res'] = NRES_NOT_SET
    p['nnscat'] = 1
    p['grid'] = 0
    p['origin'] = p['origin_orig'] = PTYPE_DUMMY
    p['path'] = 0.0
    p['ds'] = 0.0
    set_unpolarized(p)


def make_photons(nphot):
    '''
    Allocate nphot photons (the analog of allocating photmain) and put
    them in a defined state.  Photon numbers np run from 0 to nphot-1.
    '''
    p = np.zeros(nphot, dtype=photon_dtype)
    init_dummy_phot(p)
    p['np'] = np.arange(nphot)
    return p


def stuff_phot(pin, pout):
    '''
    Mirror of stuff_phot in phot_util.c: copy photon(s) pin into pout.
    '''
    pout[...] = pin
    return 0
