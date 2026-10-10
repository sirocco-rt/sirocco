'''
Step 3 tests: polarized resonance-line scattering, Sobolev trapping in
homologous flow, and the effective-polarizability approximation.  Photon
numbers are kept modest; step3_validate.py runs the same cases with more
photons.

Run with:  python -m pytest py_progs/polarization/tests
'''

import numpy as np
import pytest

from polarization import photon
from polarization.stokes import normalize, dot, rotate_stokes, observer_direction, sky_north, to_sky
from polarization.es_scatter import thomson_phase
from polarization.line_scatter import w2_factor, line_phase, scatter_line, peel_line
from polarization.homologous import (HomologousEllipsoid, run, thin_line_brown_mclean, trap_factor,
                                     p_escape)

NSIG = 4.0
VMAX = 0.001


@pytest.fixture
def rng():
    return np.random.default_rng(31)


def beam(n, lmn, q=0.0, u=0.0):
    p = photon.make_photons(n)
    p['lmn'] = lmn
    photon.set_unpolarized(p)
    p['stokes'][:, 0] = q
    p['stokes'][:, 1] = u
    p['w'] = 1.0
    return p


# ---------------------------------------------------------------- W2 and the line matrix

@pytest.mark.parametrize('jl,ju,w2', [(0, 1, 1.0), (0.5, 0.5, 0.0), (0.5, 1.5, 0.5), (1, 1, 0.25),
                                      (1, 0, 0.0), (1, 2, 0.35), (1.5, 1.5, 0.32), (1.5, 0.5, 0.0),
                                      (1.5, 2.5, 0.28), (2, 1, 0.01), (2, 2, 0.35), (2, 3, 0.24)])
def test_w2_factor(jl, ju, w2):
    assert np.isclose(w2_factor(jl, ju), w2, atol=1e-12)


def test_line_phase_limits(rng):
    c = rng.uniform(-1, 1, 100)
    q, u = rng.uniform(-0.6, 0.6, (2, 100))
    v = np.zeros(100)
    assert np.allclose(np.array(line_phase(c, q, u, v, 1.0)[:3]), np.array(thomson_phase(c, q, u)[:3]))
    f, qo, uo, vo = line_phase(c, q, u, v, 0.0)
    assert np.allclose(f, 1.0) and np.allclose(qo, 0.0) and np.allclose(uo, 0.0)


@pytest.mark.parametrize('w2,q,u', [(0.5, 0.0, 0.0), (0.25, 0.7, -0.2), (0.8, -0.3, 0.6)])
def test_line_phase_normalized(w2, q, u):
    mu, wmu = np.polynomial.legendre.leggauss(64)
    phi = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    M, PHI = np.meshgrid(mu, phi, indexing='ij')
    qp, up = rotate_stokes(q, u, PHI)
    f, *_ = line_phase(M, qp, up, 0.0, w2)
    assert np.isclose(np.sum(wmu[:, None] * f) / len(phi) / 2.0, 1.0)


def test_sampled_line_scatters_agree_with_peel_off(rng):
    '''The scatter_line / peel_line analog of the Step 1 cone test.'''
    nphot = 2000000
    w2 = 0.5
    n_in = normalize(np.array([0.3, -0.5, 0.8]))
    for angle, phase in [(90, 0.0), (40, 0.37)]:
        n_obs = observer_direction(angle, phase)
        N = sky_north(angle, phase)
        f_peel, q_peel, u_peel = peel_line(beam(1, n_in, 0.6, 0.3), n_obs, N, w2)
        p = beam(nphot, n_in, 0.6, 0.3)
        scatter_line(p, rng, w2)
        cone = np.radians(4.0)
        sel = dot(p['lmn'], n_obs) > np.cos(cone)
        frac = f_peel[0] * (1 - np.cos(cone)) / 2
        assert abs(sel.sum() / nphot - frac) < 5 * np.sqrt(frac / nphot)
        qs, us = to_sky(p['stokes'][sel, 0], p['stokes'][sel, 1], p['pol_ref'][sel], p['lmn'][sel], N)
        tol = 5 / np.sqrt(sel.sum()) + 0.01
        assert abs(np.mean(qs) - q_peel[0]) < tol and abs(np.mean(us) - u_peel[0]) < tol


@pytest.mark.parametrize('w2', [1.0, 0.5])
def test_each_local_scatter_keeps_seven_tenths_w2(rng, w2):
    '''
    The K = 2 part of the radiation, seen as Q toward an edge-on observer,
    falls by (7/10) W2 per scatter.  This is the basis of trap_factor.
    '''
    n_obs, N = observer_direction(90, 0), sky_north(90, 0)
    p = beam(1000000, [0, 0, 1.0])
    Q = []
    for _ in range(3):
        f, q, u = peel_line(p, n_obs, N, w2)
        Q.append(np.mean(f * q))
        scatter_line(p, rng, w2)
    assert abs(Q[1] / Q[0] - 0.7 * w2) < 0.02
    assert abs(Q[2] / Q[1] - 0.7 * w2) < 0.04


# ---------------------------------------------------------------- homologous transport

def test_resonance_distance_gives_rest_frequency(rng):
    g = HomologousEllipsoid(VMAX, VMAX, 0.0, 1.0)
    n = 1000
    x = normalize(rng.normal(size=(n, 3))) * rng.uniform(0, 0.5 * VMAX, n)[:, None]
    lmn = normalize(rng.normal(size=(n, 3)))
    nu = 1.0 + rng.uniform(0.0, VMAX, n)
    s = g.s_res(x, lmn, nu)
    ok = np.isfinite(s)
    x_res = x[ok] + s[ok, None] * lmn[ok]
    assert np.allclose(nu[ok] * (1.0 - dot(lmn[ok], x_res)), 1.0, atol=1e-14)


def test_trap_factor_limits():
    assert np.isclose(trap_factor(1.0, 1.0), 1.0)          # thin line: no trapping
    assert np.isclose(trap_factor(0.5, 0.0), 0.5)          # W2 = 0: only first scatters count
    assert np.isclose(trap_factor(1e-6, 1.0), 1e-6 / 0.3, rtol=1e-3)


def test_thin_line_matches_brown_mclean_equivalent(rng):
    g = HomologousEllipsoid(VMAX, 0.3 * VMAX, 0.0, 0.02, w2=1.0, trap='single')
    specs, _ = run(g, 4000000, [90], rng)
    q, u, sq, su, _, _ = specs[0].integrated()
    pred = thin_line_brown_mclean(g, 90)
    assert abs(q - pred) < NSIG * sq + 0.02 * abs(pred)
    assert abs(u) < NSIG * su


def test_effective_trapping_matches_explicit(rng):
    '''
    Thick line (tau = 3): the effective polarizability reproduces the exact
    (explicit) trapping, while a single polarized scatter is far too high.
    '''
    res = {}
    for mode in ('explicit', 'single', 'effective'):
        g = HomologousEllipsoid(VMAX, 0.3 * VMAX, 0.0, 3.0, w2=1.0, trap=mode)
        specs, _ = run(g, 500000, [90], rng)
        res[mode] = specs[0].integrated()[4:]
    (qe, se), (qs, ss), (qf, sf) = res['explicit'], res['single'], res['effective']
    assert abs(qf - qe) < NSIG * np.hypot(se, sf)
    assert (qs - qe) > 10 * np.hypot(se, ss)
    F = trap_factor(float(p_escape(3.0)), 1.0)
    assert abs(qe / qs - F) < NSIG * se / qs


def test_homologous_sphere_has_no_net_polarization(rng):
    g = HomologousEllipsoid(VMAX, VMAX, 500.0, 3.0, w2=1.0, trap='explicit')
    specs, _ = run(g, 300000, [30, 90], rng)
    for s in specs:
        q, u, sq, su, _, _ = s.integrated()
        assert abs(q) < NSIG * sq and abs(u) < NSIG * su


def test_homologous_peel_off_agrees_with_live_or_die(rng):
    nbins = 3
    mu_c = (np.arange(nbins) + 0.5) / nbins
    g = HomologousEllipsoid(VMAX, 0.3 * VMAX, 300.0, 3.0, w2=1.0, trap='explicit')
    specs, live = run(g, 400000, np.degrees(np.arccos(mu_c)), rng, live_bins=2 * nbins)
    mu, I, q, u, sq, su = live.result()
    for k, s in zip(np.where(mu > 0)[0], specs):
        kk = len(mu) - 1 - k
        qp, up, sqp, sup, _, _ = s.integrated()
        ql = 0.5 * (q[k] + q[kk])
        sql = 0.5 * np.hypot(sq[k], sq[kk])
        assert abs(ql - qp) < NSIG * np.hypot(sqp, sql)
