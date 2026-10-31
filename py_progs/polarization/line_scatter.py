#!/usr/bin/env python
'''
Polarized resonance-line scattering: the polarizability factor W2, the line
phase matrix, sampling of a new direction, and peel-off toward an observer.


Command line usage (if any):

    None; this is a module, imported as ``from polarization import line_scatter``

Description:

    Resonance scattering in a line J_l -> J_u -> J_l is a mixture of a
    Rayleigh (dipole) part with weight W2 and an isotropic, depolarizing part
    with weight 1 - W2 (Hamilton 1947).  In the scattering-plane basis used
    for electron scattering (es_scatter.py), and with q_perp, u_perp the
    incoming Stokes parameters referred to e_perp,

        I = W2 (3/4) [(1 + c^2) + s2 q_perp] + (1 - W2)
        Q = W2 (3/4) [s2 + (1 + c^2) q_perp]
        U = W2 (3/4) 2c u_perp

    W2 depends only on the angular momenta of the two levels,

        W2 = 3 (2 J_u + 1) {1 1 2; J_u J_u J_l}^2

    which gives 1 for J 0 -> 1, 0 for 1/2 -> 1/2 and 1/2 for 1/2 -> 3/2.

Primary routines:

    w2_factor, line_phase, scatter_line, peel_line, scattered_state

Notes:

    The direction of a scattered photon is drawn from the mixture (the
    Thomson distribution with probability W2, isotropic otherwise), and its
    Stokes parameters are then those of the full mixture matrix in that
    direction.  V is carried with the Rayleigh weight; it stays 0 for the
    unpolarized sources used here.

History:

    261010 ksl Coding begun, with Claude Code

'''

from fractions import Fraction
from math import factorial, sqrt

import numpy as np

from .stokes import normalize, dot, rotate_stokes, rotation_angle
from .es_scatter import sample_thomson_cos, sample_thomson_phi, PERP_EPS


def _triangle(a, b, c):
    return (factorial(int(a + b - c)) * factorial(int(a - b + c)) * factorial(int(-a + b + c))
            / factorial(int(a + b + c + 1)))


def wigner_6j(j1, j2, j3, j4, j5, j6):
    '''
    The Wigner 6j symbol {j1 j2 j3; j4 j5 j6}, by the Racah formula.
    Arguments may be integers or half integers (floats or Fractions).
    '''
    j = [Fraction(x).limit_denominator(2) for x in (j1, j2, j3, j4, j5, j6)]
    for a, b, c in [(j[0], j[1], j[2]), (j[0], j[4], j[5]), (j[3], j[1], j[5]), (j[3], j[4], j[2])]:
        if c > a + b or c < abs(a - b) or (a + b + c).denominator != 1:
            return 0.0
    pre = sqrt(_triangle(j[0], j[1], j[2]) * _triangle(j[0], j[4], j[5])
               * _triangle(j[3], j[1], j[5]) * _triangle(j[3], j[4], j[2]))
    total = 0.0
    for t in range(0, int(sum(j)) + 2):
        args = [t - j[0] - j[1] - j[2], t - j[0] - j[4] - j[5], t - j[3] - j[1] - j[5],
                t - j[3] - j[4] - j[2], j[0] + j[1] + j[3] + j[4] - t,
                j[1] + j[2] + j[4] + j[5] - t, j[2] + j[0] + j[5] + j[3] - t]
        if min(args) < 0:
            continue
        d = 1
        for a in args:
            d *= factorial(int(a))
        total += (-1) ** t * factorial(t + 1) / d
    return pre * total


def w2_factor(J_l, J_u):
    '''
    The polarizability W2 of resonance scattering J_l -> J_u -> J_l.
    J = (g - 1)/2 in terms of the statistical weights in the atomic data.
    '''
    return 3.0 * (2 * J_u + 1) * wigner_6j(1, 1, 2, J_u, J_u, J_l) ** 2


def line_phase(c, q_perp, u_perp, v_perp, w2):
    '''
    Apply the line phase matrix to a normalized Stokes vector (1, q, u, v)
    already referred to e_perp.  Returns (f, q_out, u_out, v_out): f = 4 pi p
    is the scattered intensity relative to isotropic scattering, and the
    outgoing Stokes parameters are referred to e_perp.  w2 = 1 gives Thomson
    scattering; w2 = 0 gives isotropic, unpolarized scattering.
    '''
    c2 = c * c
    s2 = 1.0 - c2
    I = w2 * 0.75 * ((1.0 + c2) + s2 * q_perp) + (1.0 - w2)
    Q = w2 * 0.75 * (s2 + (1.0 + c2) * q_perp)
    U = w2 * 0.75 * 2.0 * c * u_perp
    V = w2 * 0.75 * 2.0 * c * v_perp
    return I, Q / I, U / I, V / I


def scattered_state(n, e1, stokes, n_new, w2):
    '''
    For photons moving along n with reference axis e1 and Stokes (q, u, v),
    scattered into directions n_new, return (f, e_perp, stokes_out), where
    stokes_out is referred to e_perp.  For forward or backward scatters,
    e_perp = e1.
    '''
    cr = np.cross(n, n_new)
    crn = np.linalg.norm(cr, axis=1)
    e_perp = np.where((crn > PERP_EPS)[:, None], cr / np.maximum(crn, PERP_EPS)[:, None], e1)
    c = np.clip(dot(n, n_new), -1.0, 1.0)
    phi = rotation_angle(e1, n, e_perp)
    q_perp, u_perp = rotate_stokes(stokes[:, 0], stokes[:, 1], phi)
    f, q, u, v = line_phase(c, q_perp, u_perp, stokes[:, 2], w2)
    return f, e_perp, np.column_stack([q, u, v])


def sample_line_direction(rng, n, e1, stokes, w2):
    '''
    Draw new directions for photons (moving along n, with reference axis e1
    and Stokes parameters stokes) from the line phase function: Thomson with
    probability w2, isotropic otherwise.
    '''
    nphot = len(n)
    w2 = np.broadcast_to(w2, (nphot,))
    rayleigh = rng.random(nphot) < w2
    n_new = normalize(rng.normal(size=(nphot, 3)))
    k = np.where(rayleigh)[0]
    if k.size:
        q = stokes[k, 0]
        u = stokes[k, 1]
        c = sample_thomson_cos(rng, k.size)
        phi = sample_thomson_phi(rng, c, q, u)
        e2 = np.cross(n[k], e1[k])
        e_perp = np.cos(phi)[:, None] * e1[k] + np.sin(phi)[:, None] * e2
        s = np.sqrt(np.maximum(0.0, 1.0 - c * c))
        n_new[k] = normalize(c[:, None] * n[k] + s[:, None] * np.cross(n[k], e_perp))
    return n_new


def scatter_line(p, rng, w2):
    '''
    Resonance-scatter the photons in p, in place: choose new directions from
    the line phase function and update stokes and pol_ref.  The weights do
    not change.  w2 may be a scalar or one value per photon.
    '''
    n_new = sample_line_direction(rng, p['lmn'], p['pol_ref'], p['stokes'], w2)
    _, e_perp, st = scattered_state(p['lmn'], p['pol_ref'], p['stokes'], n_new, w2)
    p['lmn'] = n_new
    p['pol_ref'] = e_perp
    p['stokes'] = st
    return 0


def peel_line(p, n_obs, north, w2):
    '''
    Peel-off of resonance-scattered photons toward an observer.  Returns
    (f, q_sky, u_sky), as peel_es does.  Any Sobolev escape probability must
    be applied separately.
    '''
    n_obs_all = np.broadcast_to(n_obs, p['lmn'].shape)
    f, e_perp, st = scattered_state(p['lmn'], p['pol_ref'], p['stokes'], n_obs_all, w2)
    q_sky, u_sky = rotate_stokes(st[:, 0], st[:, 1], rotation_angle(e_perp, n_obs_all, north))
    return f, q_sky, u_sky
