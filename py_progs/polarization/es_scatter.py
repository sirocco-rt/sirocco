#!/usr/bin/env python
'''
Polarized electron (Thomson) scattering: the phase matrix, sampling of a
new direction, and peel-off toward an observer.


Command line usage (if any):

    None; this is a module, imported as ``from polarization import es_scatter``

Description:

    The matrix is written in the scattering-plane basis, e1 = e_perp =
    normalize(n x n'), which is perpendicular to both the incoming direction
    n and the outgoing direction n'.  With c = cos(Theta), s2 = 1 - c^2:

                  3 | 1+c^2  s2     0   0 |
        M     =   - | s2     1+c^2  0   0 |
                  4 | 0      0     2c   0 |
                    | 0      0      0  2c |

    normalized so that the intensity row averages to 1 over the sphere for
    unpolarized light.  (Chandrasekhar's parallel/perpendicular basis flips
    the sign of the off-diagonal s2 terms.)  Because e_perp is perpendicular
    to n' as well as n, it serves directly as the reference axis of the
    scattered photon and no rotation is needed afterwards.

Primary routines:

    scatter_es, peel_es, thomson_phase, sample_thomson_cos,
    sample_thomson_phi

Notes:

    Only the Thomson limit is implemented.  The Klein-Nishina (Fano) matrix
    is planned; it should replace thomson_matrix / thomson_phase without
    changing the bookkeeping.  See the Obsidian note
    Sirocco-2026-Polarization-Conventions, sections 5a and 6.

History:

    261010 ksl Coding begun, with Claude Code

'''

import numpy as np

from .stokes import normalize, dot, rotate_stokes, rotation_angle

PERP_EPS = 1e-6    # |n x n'| below which a scatter counts as forward/back


def thomson_matrix(c):
    '''Return the 4x4 Thomson phase matrix for cos(Theta) = c (scalar).'''
    s2 = 1.0 - c * c
    return 0.75 * np.array([[1 + c * c, s2, 0, 0],
                            [s2, 1 + c * c, 0, 0],
                            [0, 0, 2 * c, 0],
                            [0, 0, 0, 2 * c]])


def thomson_phase(c, q_perp, u_perp, v_perp=0.0):
    '''
    Apply the Thomson matrix to a normalized Stokes vector (1, q, u, v)
    that is already referred to e_perp.

    Returns (f, q_out, u_out, v_out), where f = 4 pi p(Theta, phi) is the
    ratio of the scattered intensity in this direction to that of isotropic
    scattering, and (q_out, u_out, v_out) are referred to e_perp.
    '''
    c2 = c * c
    s2 = 1.0 - c2
    D = (1.0 + c2) + s2 * q_perp
    q_out = (s2 + (1.0 + c2) * q_perp) / D
    u_out = 2.0 * c * u_perp / D
    v_out = 2.0 * c * v_perp / D
    return 0.75 * D, q_out, u_out, v_out


def sample_thomson_cos(rng, n):
    '''
    Draw n values of cos(Theta) from the unpolarized Thomson distribution,
    p(mu) = (3/8)(1 + mu^2), by inverting the CDF exactly.  The CDF gives
    mu^3 + 3 mu - (8 xi - 4) = 0, which has one real root:

        mu = cbrt(A + sqrt(A^2 + 1)) + cbrt(A - sqrt(A^2 + 1)),  A = 4 xi - 2

    This is also the exact marginal distribution of Theta for polarized
    light, because the terms in phi average to zero.
    '''
    A = 4.0 * rng.random(n) - 2.0
    r = np.sqrt(A * A + 1.0)
    return np.cbrt(A + r) + np.cbrt(A - r)


def sample_thomson_phi(rng, c, q, u):
    '''
    Draw phi, the angle of e_perp measured from the photon's e1 toward e2,
    from the conditional distribution given cos(Theta) = c:

        p(phi | Theta) ~ 1 + a (q cos 2phi + u sin 2phi),  a = s2/(1 + c^2)

    by rejection with envelope 1 + a P, P = sqrt(q^2 + u^2) <= 1.
    c, q, u are arrays of the same length.
    '''
    a = (1.0 - c * c) / (1.0 + c * c)
    pmax = 1.0 + a * np.sqrt(q * q + u * u)
    phi = np.empty_like(c)
    todo = np.arange(len(c))
    while todo.size:
        trial = rng.uniform(0.0, 2.0 * np.pi, todo.size)
        f = 1.0 + a[todo] * (q[todo] * np.cos(2 * trial) + u[todo] * np.sin(2 * trial))
        ok = rng.random(todo.size) * pmax[todo] < f
        phi[todo[ok]] = trial[ok]
        todo = todo[~ok]
    return phi


def scatter_es(p, rng):
    '''
    Electron-scatter the photons in p, in place: choose new directions from
    the polarized Thomson phase function and update stokes and pol_ref.
    The weights do not change, because the directions are drawn from the
    exact phase function.

    This is the polarized replacement for the direction choice that
    scatter() in resonate.c makes for nres == NRES_ES.  As in scatter(),
    the photons are not moved.
    '''
    nphot = len(p)
    n = p['lmn']
    e1 = p['pol_ref']
    e2 = np.cross(n, e1)
    q = p['stokes'][:, 0]
    u = p['stokes'][:, 1]
    v = p['stokes'][:, 2]

    c = sample_thomson_cos(rng, nphot)
    phi = sample_thomson_phi(rng, c, q, u)
    s = np.sqrt(np.maximum(0.0, 1.0 - c * c))

    e_perp = np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
    n_new = c[:, None] * n + s[:, None] * np.cross(n, e_perp)

    q_perp, u_perp = rotate_stokes(q, u, phi)
    f, q_out, u_out, v_out = thomson_phase(c, q_perp, u_perp, v)

    p['lmn'] = normalize(n_new)
    p['pol_ref'] = e_perp
    p['stokes'] = np.column_stack([q_out, u_out, v_out])
    return 0


def peel_es(p, n_obs, north):
    '''
    Peel-off of electron-scattered photons toward an observer.

    p holds photons at their scattering points, with lmn, stokes and
    pol_ref as they were just before the scatter.  n_obs is the direction
    to the observer and north the sky North vector for that observer.

    Returns (f, q_sky, u_sky): f = 4 pi p(Theta, phi) is the factor by which
    the photon weight must be multiplied (for unpolarized light it reduces
    to the 0.75(1 + cos^2) now used in extract.c), and (q_sky, u_sky) are
    referred to North.  The tallies are then
    I += w f exp(-tau),  Q += w f exp(-tau) q_sky,  U += w f exp(-tau) u_sky.
    '''
    n = p['lmn']
    e1 = p['pol_ref']
    q = p['stokes'][:, 0]
    u = p['stokes'][:, 1]
    v = p['stokes'][:, 2]
    n_obs = np.broadcast_to(n_obs, n.shape)

    c = np.clip(dot(n, n_obs), -1.0, 1.0)
    cr = np.cross(n, n_obs)
    crn = np.linalg.norm(cr, axis=1)
    e_perp = np.where((crn > PERP_EPS)[:, None], cr / np.maximum(crn, PERP_EPS)[:, None], e1)

    phi = rotation_angle(e1, n, e_perp)
    q_perp, u_perp = rotate_stokes(q, u, phi)
    f, q_out, u_out, _ = thomson_phase(c, q_perp, u_perp, v)

    q_sky, u_sky = rotate_stokes(q_out, u_out, rotation_angle(e_perp, n_obs, north))
    return f, q_sky, u_sky
