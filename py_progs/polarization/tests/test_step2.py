'''
Step 2 tests: toy polarized transport against known results.  Photon
numbers are kept modest so these run in about a minute; step2_validate.py
runs the same cases with more photons.

Run with:  python -m pytest py_progs/polarization/tests
'''

import numpy as np
import pytest

from polarization.toy_transport import Ellipsoid, Slab, run, brown_mclean
from polarization.milne import solve_slab

NSIG = 4.0


@pytest.fixture
def rng():
    return np.random.default_rng(1010)


@pytest.fixture(scope='module')
def milne10():
    return solve_slab(10.0)


# ---------------------------------------------------------------- deterministic reference

def test_milne_limb_polarization_is_chandrasekhar(milne10):
    '''Chandrasekhar (1960): P(mu = 0) = 11.71% for the Milne problem.'''
    _, q = milne10(1e-4)
    assert abs(-100 * q[0] - 11.71) < 0.03


def test_milne_converged_and_independent_of_thickness(milne10):
    mu = np.array([0.01, 0.1, 0.5, 0.9])
    _, q10 = milne10(mu)
    _, q40 = solve_slab(40.0, ntau=400, nmu=32)(mu)
    assert np.allclose(q10, q40, atol=2e-5)
    assert np.all(np.diff(-q10) < 0)            # P falls monotonically from the limb


# ---------------------------------------------------------------- Monte Carlo cases

def test_sphere_has_no_net_polarization(rng):
    specs, _ = run(Ellipsoid(1.0, 1.0, 2.0), 300000, [20, 60, 90], rng)
    for s in specs:
        q, u, sq, su = s.result()
        assert abs(q) < NSIG * sq and abs(u) < NSIG * su


@pytest.mark.parametrize('a,c,sign', [(1.0, 0.3, +1), (0.3, 1.0, -1)])
def test_thin_ellipsoid_matches_brown_mclean(rng, a, c, sign):
    '''Oblate gives Q > 0 and prolate Q < 0, with the Brown & McLean amplitude.'''
    geom = Ellipsoid(a, c, 0.05)
    specs, _ = run(geom, 1500000, [45, 90], rng)
    for s in specs:
        q, u, sq, su = s.result()
        qbm = brown_mclean(geom, s.angle)[0]
        assert np.sign(q) == sign
        # allow O(tau) departures from the first-order formula
        assert abs(q - qbm) < NSIG * sq + 0.04 * abs(qbm)
        assert abs(u) < NSIG * su


def test_live_or_die_agrees_with_peel_off(rng):
    nbins = 5
    mu_c = (np.arange(nbins) + 0.5) / nbins
    specs, live = run(Ellipsoid(1.0, 0.3, 1.0), 400000, np.degrees(np.arccos(mu_c)), rng,
                      live_bins=2 * nbins)
    mu, I_l, q_l, u_l, sq_l, su_l = live.result()
    for k, s in zip(np.where(mu > 0)[0], specs):
        kk = len(mu) - 1 - k                        # fold the two hemispheres
        q, u, sq, su = s.result()
        ql = 0.5 * (q_l[k] + q_l[kk])
        sql = 0.5 * np.hypot(sq_l[k], sq_l[kk])
        assert abs(ql - q) < NSIG * np.hypot(sq, sql)
        # within a bin I varies by a few per cent; compare loosely
        assert abs(0.5 * (I_l[k] + I_l[kk]) / s.I - 1) < 0.02


def test_slab_matches_deterministic_solution(rng, milne10):
    angles = [0, 60, 84]
    specs, _ = run(Slab(10.0), 200000, angles, rng)
    for s in specs:
        m = np.cos(np.radians(s.angle))
        I0, q0 = milne10(m)
        q, u, sq, su = s.result()
        sI = s.errors()[0]
        assert abs(s.I / m - I0[0]) < NSIG * sI / m
        assert abs(q - q0[0]) < NSIG * sq
        assert abs(u) < NSIG * su
