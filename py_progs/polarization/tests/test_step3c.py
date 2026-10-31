'''
Step 3c tests: the capped explicit trapping chain proposed for the C
version.  Photon numbers are kept modest; step3c_validate.py runs the full
study.

Run with:  python -m pytest py_progs/polarization/tests
'''

import numpy as np
import pytest

from polarization.sobolev_flow import HomologousEllipsoid, RadialFlowEllipsoid, run

NSIG = 4.0
VMAX = 0.001
NMAX_INF = 10 ** 9


def radial(trap, nmax=30, alpha=0.25, tau=3.0, kappa=0.0, line_frac=0.0):
    return RadialFlowEllipsoid(VMAX, 0.3 * VMAX, VMAX, alpha, kappa, tau, line_frac=line_frac, trap=trap, nmax=nmax)


@pytest.mark.parametrize('make', [
    lambda trap, nmax: radial(trap, nmax, kappa=200.0, line_frac=0.5),
    lambda trap, nmax: HomologousEllipsoid(VMAX, 0.3 * VMAX, 200.0, 3.0, line_frac=0.5, trap=trap, nmax=nmax),
])
def test_uncapped_chain_reproduces_explicit_exactly(make):
    '''
    Summing the peel-off contributions locally and attenuating once per
    observer is the same as peeling at every local scatter.
    '''
    out = []
    for trap in ('explicit', 'capped'):
        specs, _ = run(make(trap, NMAX_INF), 30000, [60, 90], np.random.default_rng(8))
        out.append(np.array([[s.I.sum(), s.Q.sum(), s.U.sum()] for s in specs]))
    assert np.allclose(out[0], out[1], rtol=1e-12, atol=1e-18)


def test_nmax_zero_is_unpolarized_escape():
    '''With nmax = 0 every line photon escapes unpolarized, as in Sirocco now.'''
    g = radial('capped', nmax=0)
    specs, _ = run(g, 20000, [90], np.random.default_rng(9))
    q_line = specs[0].integrated()[4]
    assert q_line == 0.0
    assert g.stats['fallback'] == g.stats['events'] and g.stats['local'] == 0


@pytest.mark.parametrize('line_frac', [0.0, 1.0])
def test_cap_of_thirty_matches_uncapped(line_frac):
    '''
    At tau = 3, nmax = 30 leaves almost no photons trapped, so q agrees with
    the uncapped chain (scattered continuum, and pure line emission).
    '''
    res = []
    for nmax in (30, NMAX_INF):
        g = radial('capped', nmax=nmax, line_frac=line_frac)
        specs, _ = run(g, 200000, [90], np.random.default_rng(10 + nmax % 7))
        q, u, sq, su, ql, sql = specs[0].integrated()
        res.append((ql, sql) if line_frac == 0.0 else (q, sq))
        if nmax == 30:
            assert g.stats['fallback'] < 0.02 * g.stats['events']
    (q1, s1), (q2, s2) = res
    assert abs(q1 - q2) < NSIG * np.hypot(s1, s2)


def test_fallback_fraction_falls_with_nmax():
    frac = []
    for nmax in (1, 3, 10):
        g = radial('capped', nmax=nmax, tau=10.0)
        run(g, 20000, [90], np.random.default_rng(11))
        frac.append(g.stats['fallback'] / g.stats['events'])
    assert frac[0] > frac[1] > frac[2]
