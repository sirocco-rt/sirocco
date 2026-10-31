'''
Step 3b tests: Sobolev trapping when the escape probability depends on
direction (RadialFlowEllipsoid).  Photon numbers are kept modest;
step3b_validate.py runs the same cases with more photons.

Run with:  python -m pytest py_progs/polarization/tests
'''

import numpy as np
import pytest

from polarization.stokes import normalize
from polarization.sobolev_flow import HomologousEllipsoid, RadialFlowEllipsoid, run, p_escape

NSIG = 4.0
VMAX = 0.001


@pytest.fixture
def rng():
    return np.random.default_rng(47)


def geom(alpha, tau, trap, line_frac=0.0, c=0.3):
    return RadialFlowEllipsoid(VMAX, c * VMAX, VMAX, alpha, 0.0, tau, line_frac=line_frac, trap=trap)


def test_bisection_finds_resonance(rng):
    g = geom(0.5, 3.0, 'explicit')
    n = 2000
    x = normalize(rng.normal(size=(n, 3))) * rng.uniform(0, 0.25 * VMAX, n)[:, None]
    lmn = normalize(rng.normal(size=(n, 3)))
    nu = 1.0 + rng.uniform(0.0, 3 * VMAX, n)
    s = g.s_res(x, lmn, nu)
    ok = np.isfinite(s)
    assert ok.sum() > 100
    assert np.allclose(g.nu_cmf(x[ok] + s[ok, None] * lmn[ok], lmn[ok], nu[ok]), 1.0, atol=1e-13)


def test_sobolev_tau_and_mean_escape():
    g = geom(0.25, 2.0, 'explicit')
    x = np.array([[0.5 * VMAX, 0.0, 0.0]])
    radial = g.tau_sobolev(x, np.array([[1.0, 0, 0]]))
    tangential = g.tau_sobolev(x, np.array([[0, 1.0, 0]]))
    assert np.isclose(radial[0] / tangential[0], 1.0 / 0.25)      # 1/(1 - (1 - alpha))
    mu = np.linspace(-1, 1, 20001)
    tau = 2.0 * 0.5 ** 0.75 / (1 - 0.75 * mu * mu)
    assert np.isclose(g.beta_mean(x)[0], np.trapz(p_escape(tau), mu) / 2, rtol=1e-4)


def test_alpha_one_reproduces_homologous_flow():
    out = []
    for g in (HomologousEllipsoid(VMAX, 0.3 * VMAX, 0.0, 3.0, trap='explicit'),
              RadialFlowEllipsoid(VMAX, 0.3 * VMAX, VMAX, 1.0, 0.0, 3.0, r_in=1e-9 * VMAX, trap='explicit')):
        specs, _ = run(g, 40000, [90], np.random.default_rng(3))
        out.append(specs[0].integrated()[4])
    assert np.isclose(out[0], out[1], rtol=1e-9)


def test_effective_tracks_explicit_with_anisotropic_escape(rng):
    '''Scattered continuum, alpha = 3, tau = 3: effective agrees, single is far too high.'''
    res = {}
    for mode in ('explicit', 'single', 'effective'):
        specs, _ = run(geom(3.0, 3.0, mode), 300000, [90], rng)
        res[mode] = specs[0].integrated()[4:]
    (qe, se), (qs, ss), (qf, sf) = res['explicit'], res['single'], res['effective']
    assert abs(qf - qe) < NSIG * np.hypot(se, sf) + 0.1 * abs(qe)
    assert qs - qe > 8 * np.hypot(se, ss)


@pytest.mark.parametrize('alpha,sign', [(0.25, +1), (3.0, -1)])
def test_trapping_polarizes_line_emission(rng, alpha, sign):
    '''
    With anisotropic escape, explicit trapping polarizes unpolarized line
    emission: tangential pattern (Q > 0 for an oblate envelope) when radial
    directions are trapped, radial pattern (Q < 0) when tangential
    directions are trapped.  The single and effective treatments give 0.
    '''
    specs, _ = run(geom(alpha, 3.0, 'explicit', line_frac=1.0), 200000, [90], rng)
    q, u, sq, su, _, _ = specs[0].integrated()
    assert sign * q > 6 * sq
    assert abs(u) < NSIG * su
    specs, _ = run(geom(alpha, 3.0, 'effective', line_frac=1.0), 50000, [90], rng)
    assert specs[0].integrated()[0] == 0.0


def test_sphere_line_emission_has_no_net_polarization(rng):
    specs, _ = run(geom(0.25, 3.0, 'explicit', line_frac=1.0, c=1.0), 200000, [30, 90], rng)
    for s in specs:
        q, u, sq, su, _, _ = s.integrated()
        assert abs(q) < NSIG * sq and abs(u) < NSIG * su
