#!/usr/bin/env python
'''
A deterministic solution for polarized electron scattering in a plane-
parallel slab, used as an independent check of the Monte Carlo slab in
toy_transport.py.

The slab has optical depth tau_max, nothing is incident on the top
(tau = 0), and unpolarized isotropic radiation of intensity I_inc is
incident on the base.  The radiation field is axisymmetric, so U = 0 and
the Stokes parameters are I_l (E-vector in the meridian plane, i.e. along
the projected normal) and I_r (perpendicular to it).  Averaging the dipole
(Thomson) redistribution over azimuth gives the source functions

    S_r(mu) = 3/8 C
    S_l(mu) = 3/8 [mu^2 C + 2 (1 - mu^2) B]

with the moments

    C = integral (I_r + mu'^2 I_l) dmu'
    B = integral (1 - mu'^2) I_l dmu'

The transfer equation is solved by short characteristics with a source
function linear in tau, and the resulting linear system for B and C on the
tau grid is solved exactly (no Lambda iteration).

In the sky convention of this package (North = projected normal),
q = (I_l - I_r)/(I_l + I_r).  For tau_max >> 1 this is Chandrasekhar's
Milne problem, whose limb polarization is q(mu -> 0) = -11.71%.
'''

import numpy as np


def _tau_grid(tau_max, ntau, tau_min=1e-5):
    t = np.concatenate([[0.0], np.geomspace(tau_min, tau_max, ntau - 1)])
    return t


def _weights(dtau, mu):
    '''Short-characteristic weights for a source linear in tau.'''
    d = dtau / abs(mu)
    ed = np.exp(-d)
    with np.errstate(divide='ignore', invalid='ignore'):
        b = np.where(d > 1e-6, (1.0 - ed * (1.0 + d)) / d, 0.5 * d)
    a = (1.0 - ed) - b
    return ed, a, b


def _formal_matrix(tau, mu):
    '''
    Return (L, g) such that I(tau_i, mu) = sum_k L[i, k] S(tau_k) + g[i] I_bc,
    where I_bc is the boundary intensity (at the base for mu > 0, at the
    top for mu < 0).  mu > 0 means travelling upward, toward tau = 0.
    '''
    n = len(tau)
    L = np.zeros((n, n))
    g = np.zeros(n)
    dt = np.diff(tau)
    if mu > 0:
        g[-1] = 1.0
        for i in range(n - 2, -1, -1):
            ed, a, b = _weights(dt[i], mu)
            L[i] = ed * L[i + 1]
            L[i, i] += a
            L[i, i + 1] += b
            g[i] = ed * g[i + 1]
    else:
        g[0] = 1.0
        for i in range(1, n):
            ed, a, b = _weights(dt[i - 1], mu)
            L[i] = ed * L[i - 1]
            L[i, i] += a
            L[i, i - 1] += b
            g[i] = ed * g[i - 1]
    return L, g


def solve_slab(tau_max=10.0, I_inc=4.0, ntau=300, nmu=24):
    '''
    Solve the slab problem.  Returns a function emergent(mu) giving
    (I, q) = (I_l + I_r, (I_l - I_r)/(I_l + I_r)) at the top surface for an
    array of mu > 0.

    I_inc = 4 matches the normalization of toy_transport.Slab, in which a
    source of total weight 1 is injected through the base and an isotropic
    source of weight 1 would give I = 1.
    '''
    tau = _tau_grid(tau_max, ntau)
    n = len(tau)
    x, wx = np.polynomial.legendre.leggauss(nmu)
    mu_up = 0.5 * (x + 1.0)              # double-Gauss on (0, 1)
    w_up = 0.5 * wx
    mus = np.concatenate([mu_up, -mu_up])
    ws = np.concatenate([w_up, w_up])

    # I_l(tau_i, mu_j) = L_j (S_l,j) + g_j bc ; S_l,j = 3/8 [mu^2 C + 2(1-mu^2) B]
    # Assemble C = sum_j w_j (I_r + mu^2 I_l), B = sum_j w_j (1 - mu^2) I_l
    # as a linear system in the unknown vector [B, C].
    M = np.zeros((2 * n, 2 * n))
    rhs = np.zeros(2 * n)
    for mu, w in zip(mus, ws):
        L, g = _formal_matrix(tau, mu)
        m2 = mu * mu
        bc_l = bc_r = 0.5 * I_inc if mu > 0 else 0.0
        # I_l = L (3/8)(m2 C + 2(1-m2) B) + g bc_l ;  I_r = L (3/8) C + g bc_r
        Il_B = L * (0.75 * (1 - m2))
        Il_C = L * (0.375 * m2)
        Ir_C = L * 0.375
        # B row: w (1 - m2) I_l
        M[:n, :n] += w * (1 - m2) * Il_B
        M[:n, n:] += w * (1 - m2) * Il_C
        rhs[:n] += w * (1 - m2) * g * bc_l
        # C row: w (I_r + m2 I_l)
        M[n:, :n] += w * m2 * Il_B
        M[n:, n:] += w * (Ir_C + m2 * Il_C)
        rhs[n:] += w * (g * bc_r + m2 * g * bc_l)
    sol = np.linalg.solve(np.eye(2 * n) - M, rhs)
    B, C = sol[:n], sol[n:]

    def emergent(mu):
        mu = np.atleast_1d(np.asarray(mu, dtype=float))
        I_l = np.empty_like(mu)
        I_r = np.empty_like(mu)
        for k, m in enumerate(mu):
            L, g = _formal_matrix(tau, m)
            S_l = 0.375 * (m * m * C + 2 * (1 - m * m) * B)
            S_r = 0.375 * C
            I_l[k] = L[0] @ S_l + g[0] * 0.5 * I_inc
            I_r[k] = L[0] @ S_r + g[0] * 0.5 * I_inc
        I = I_l + I_r
        return I, (I_l - I_r) / I

    return emergent
