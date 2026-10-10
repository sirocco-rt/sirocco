#!/usr/bin/env python
'''
Stokes-parameter bookkeeping for polarized Monte Carlo transport.

A photon moving along n carries (q, u, v) relative to a reference axis
e1 (``pol_ref``), with e2 = n x e1, so (e1, e2, n) is right handed.
Q > 0 means the E-vector lies along e1; U > 0 means along (e1 + e2)/sqrt(2).

Observer (sky) frame: the astronomical (IAU 1974) convention, with North
the projection of the system +z axis on the sky and East = n_obs x North,
so that Q > 0 is a North-South E-vector and position angle increases from
North through East.

See the Obsidian note Sirocco-2026-Polarization-Conventions, sections 2-4
and 7.  All routines work on single vectors (3,) or arrays (N, 3).
'''

import numpy as np

C = 2.99792458e10    # VLIGHT in Sirocco, cm/s

REL_MODE_LINEAR = 0
REL_MODE_FULL = 1


def normalize(v):
    '''Return v / |v| along the last axis.'''
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def dot(a, b):
    '''Dot product along the last axis.'''
    return np.sum(np.asarray(a) * np.asarray(b), axis=-1)


def rotate_stokes(q, u, alpha):
    '''
    Return (q', u') after the reference axis is turned by alpha (radians)
    from e1 toward e2, i.e. e1' = cos(alpha) e1 + sin(alpha) e2.
    v is unchanged by a rotation.
    '''
    c2 = np.cos(2.0 * alpha)
    s2 = np.sin(2.0 * alpha)
    return q * c2 + u * s2, -q * s2 + u * c2


def rotation_angle(e1, n, target):
    '''
    Return the angle alpha that turns the reference axis e1 of a photon
    moving along n onto the axis target (which must be perpendicular to n).
    '''
    e2 = np.cross(n, e1)
    return np.arctan2(dot(target, e2), dot(target, e1))


def reorthogonalize(e1, n):
    '''
    Return e1 with any component along n removed, renormalized.  This is
    the guard applied after a frame transformation.
    '''
    return normalize(e1 - dot(e1, n)[..., None] * n)


def observer_direction(angle, phase=0.0):
    '''
    Direction toward an observer at inclination angle (degrees) and orbital
    phase, exactly as set up for xxspec[].lmn in spectra.c.
    '''
    i = np.radians(angle)
    psi = -2.0 * np.pi * phase
    return np.array([np.sin(i) * np.cos(psi), np.sin(i) * np.sin(psi), np.cos(i)])


def sky_north(angle, phase=0.0):
    '''
    Unit vector toward North on the sky of an observer at (angle, phase):
    the projection of +z, written as -d(n_obs)/d(angle) so that it is still
    defined for a pole-on (angle = 0) observer.
    '''
    i = np.radians(angle)
    psi = -2.0 * np.pi * phase
    return np.array([-np.cos(i) * np.cos(psi), -np.cos(i) * np.sin(psi), np.sin(i)])


def sky_east(angle, phase=0.0):
    '''East = n_obs x North: to the left when North is up, facing the source.'''
    return np.cross(observer_direction(angle, phase), sky_north(angle, phase))


def to_sky(q, u, e1, n_obs, north):
    '''
    Return (q, u) referred to North for photons travelling along n_obs
    whose current reference axis is e1.
    '''
    return rotate_stokes(q, u, rotation_angle(e1, n_obs, north))


def pol_degree_angle(I, Q, U):
    '''
    Return the degree of polarization P = sqrt(Q^2 + U^2)/I and the position
    angle in degrees, measured from North through East, in [0, 180).
    '''
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.sqrt(Q * Q + U * U) / I
    pa = np.degrees(0.5 * np.arctan2(U, Q)) % 180.0
    return P, pa


def boost_pol_ref(e1, n, v):
    '''
    Transform the reference axis e1 of a photon moving along n into a frame
    moving with velocity v (cm/s), treating e1 as the direction of the
    wave's electric field:

        E' = gamma (E + beta x B) - gamma^2/(gamma+1) beta (beta . E),  B = n x E

    The result is perpendicular to the aberrated direction that
    lorentz_transform in frame.c gives for lmn.  Aberration is conformal,
    so q, u and v are unchanged; only the reference axis changes.
    '''
    beta = np.asarray(v, dtype=float) / C
    g = 1.0 / np.sqrt(1.0 - dot(beta, beta))
    g = np.asarray(g)[..., None]
    b = np.cross(n, e1)
    e = g * (e1 + np.cross(beta, b)) - g * g / (g + 1.0) * beta * dot(beta, e1)[..., None]
    return normalize(e)


def lorentz_transform(p, v, rel_mode=REL_MODE_FULL):
    '''
    Mirror of lorentz_transform in frame.c, applied in place to the photons
    in p, with the reference axis transformed as well.

    v is the velocity of the frame being transformed into (cm/s), either
    one vector (3,) or one per photon (N, 3).  In REL_MODE_LINEAR only the
    frequency changes, as in frame.c, so pol_ref is left alone.
    '''
    v = np.asarray(v, dtype=float)
    n = p['lmn']
    vel = dot(n, v)
    f_in = p['freq']

    if rel_mode == REL_MODE_LINEAR:
        local = p['frame'] == 0
        p['freq'] = np.where(local, f_in * (1.0 - vel / C), f_in / (1.0 + vel / C))
        return 0

    gamma = 1.0 / np.sqrt(1.0 - dot(v, v) / (C * C))
    f_out = f_in * gamma * (1.0 - vel / C)
    x = gamma / C * (1.0 - gamma * vel / ((gamma + 1.0) * C))

    e1_new = boost_pol_ref(p['pol_ref'], n, v)
    lmn_new = (f_in / f_out)[..., None] * (n - np.asarray(x)[..., None] * v)

    p['freq'] = f_out
    p['lmn'] = lmn_new
    p['w'] *= f_out / f_in
    p['pol_ref'] = reorthogonalize(e1_new, lmn_new)
    return 0
