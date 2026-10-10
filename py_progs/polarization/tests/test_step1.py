'''
Step 1 tests: Stokes bookkeeping, the Thomson matrix, direction sampling,
peel-off, the sky frame and the Lorentz transform of the reference axis.

Run with:  python -m pytest py_progs/polarization/tests
'''

import numpy as np
import pytest
from scipy import stats

from polarization import photon
from polarization.stokes import (rotate_stokes, rotation_angle, observer_direction, sky_north, sky_east,
                    to_sky, pol_degree_angle, boost_pol_ref, lorentz_transform, normalize, dot, C,
                    REL_MODE_LINEAR)
from polarization.es_scatter import (thomson_matrix, thomson_phase, sample_thomson_cos, sample_thomson_phi,
                        scatter_es, peel_es)


@pytest.fixture
def rng():
    return np.random.default_rng(20261010)


def random_unit(rng, n):
    return normalize(rng.normal(size=(n, 3)))


def beam(n, lmn, q=0.0, u=0.0, pol_ref=None):
    '''n photons along lmn with the given polarization.'''
    p = photon.make_photons(n)
    p['lmn'] = lmn
    photon.set_unpolarized(p)
    if pol_ref is not None:
        p['pol_ref'] = pol_ref
    p['stokes'][:, 0] = q
    p['stokes'][:, 1] = u
    p['w'] = 1.0
    return p


# ---------------------------------------------------------------- photons

def test_default_pol_ref_is_perpendicular_unit(rng):
    n = np.vstack([random_unit(rng, 1000), [[0, 0, 1], [0, 0, -1]]])
    e1 = photon.default_pol_ref(n)
    assert np.allclose(np.linalg.norm(e1, axis=1), 1.0)
    assert np.allclose(dot(e1, n), 0.0, atol=1e-12)


def test_make_photons_unpolarized():
    p = photon.make_photons(10)
    assert np.all(p['stokes'] == 0.0)
    assert np.allclose(dot(p['pol_ref'], p['lmn']), 0.0)
    assert list(p['np']) == list(range(10))


# ---------------------------------------------------------------- rotations

def test_rotation_by_pi_is_identity(rng):
    q, u = rng.uniform(-0.7, 0.7, (2, 100))
    q2, u2 = rotate_stokes(q, u, np.pi)
    assert np.allclose(q2, q) and np.allclose(u2, u)


def test_rotations_compose_and_preserve_p(rng):
    q, u = rng.uniform(-0.7, 0.7, (2, 100))
    a, b = rng.uniform(-3, 3, (2, 100))
    q1, u1 = rotate_stokes(*rotate_stokes(q, u, a), b)
    q2, u2 = rotate_stokes(q, u, a + b)
    assert np.allclose(q1, q2) and np.allclose(u1, u2)
    assert np.allclose(np.hypot(q1, u1), np.hypot(q, u))


def test_rotation_angle_maps_e1_to_target(rng):
    n = random_unit(rng, 200)
    e1 = photon.default_pol_ref(n)
    target = normalize(np.cross(n, random_unit(rng, 200)))
    alpha = rotation_angle(e1, n, target)
    e2 = np.cross(n, e1)
    new = np.cos(alpha)[:, None] * e1 + np.sin(alpha)[:, None] * e2
    assert np.allclose(new, target)


def test_fully_polarized_along_target_gives_q_one(rng):
    '''A photon with q=1 along e1 must have q=1 in any basis whose e1 is the same axis.'''
    n = random_unit(rng, 50)
    e1 = photon.default_pol_ref(n)
    q, u = to_sky(np.ones(50), np.zeros(50), e1, n, -e1)   # sign of the axis is irrelevant
    assert np.allclose(q, 1.0) and np.allclose(u, 0.0, atol=1e-12)


# ---------------------------------------------------------------- sky frame (IAU)

@pytest.mark.parametrize('angle,phase', [(0, 0), (0, 0.3), (30, 0.1), (62, 0.5), (90, 0), (135, 0.75)])
def test_sky_frame_is_orthonormal_right_handed(angle, phase):
    n = observer_direction(angle, phase)
    N = sky_north(angle, phase)
    E = sky_east(angle, phase)
    assert np.isclose(np.linalg.norm(N), 1) and np.isclose(dot(N, n), 0)
    assert np.isclose(dot(N, E), 0) and np.isclose(np.linalg.norm(E), 1)
    assert np.allclose(np.cross(N, E), n)     # (N, E, n) right handed


def test_sky_north_is_projected_z_axis():
    for angle in (10, 45, 80, 90):
        n = observer_direction(angle, 0.2)
        proj = normalize(np.array([0, 0, 1.0]) - n[2] * n)
        assert np.allclose(sky_north(angle, 0.2), proj)


def test_iau_east_is_on_the_left():
    '''Observer on +x looking back at the origin with z up: left is -y.'''
    assert np.allclose(observer_direction(90, 0), [1, 0, 0])
    assert np.allclose(sky_east(90, 0), [0, -1, 0])


@pytest.mark.parametrize('e_vec,pa', [('N', 0.0), ('E', 90.0), ('NE', 45.0)])
def test_position_angle_north_through_east(e_vec, pa):
    angle, phase = 60.0, 0.3
    n = observer_direction(angle, phase)
    N, E = sky_north(angle, phase), sky_east(angle, phase)
    axis = {'N': N, 'E': E, 'NE': normalize(N + E)}[e_vec]
    e1 = photon.default_pol_ref(n)
    alpha = rotation_angle(e1, n, axis)
    q, u = rotate_stokes(1.0, 0.0, -alpha)     # fully polarized along axis, in the e1 basis
    q_sky, u_sky = to_sky(q, u, e1, n, N)
    P, theta = pol_degree_angle(1.0, q_sky, u_sky)
    assert np.isclose(P, 1.0)
    assert np.isclose(theta, pa, atol=1e-9) or np.isclose(theta, pa + 180.0, atol=1e-9)


# ---------------------------------------------------------------- Thomson matrix

def test_matrix_and_phase_agree(rng):
    for _ in range(20):
        c = rng.uniform(-1, 1)
        q, u, v = rng.uniform(-0.5, 0.5, 3)
        S = thomson_matrix(c) @ np.array([1.0, q, u, v])
        f, qo, uo, vo = thomson_phase(c, q, u, v)
        assert np.isclose(S[0], f)
        assert np.allclose(S[1:] / S[0], [qo, uo, vo])


def test_unpolarized_90_degree_scatter_is_fully_polarized():
    f, q, u, v = thomson_phase(0.0, 0.0, 0.0)
    assert np.isclose(f, 0.75) and np.isclose(q, 1.0) and u == 0.0


def test_forward_scatter_preserves_polarization():
    f, q, u, v = thomson_phase(1.0, 0.3, -0.4)
    assert np.isclose(f, 1.5) and np.isclose(q, 0.3) and np.isclose(u, -0.4)


def test_back_scatter_mirrors_u():
    f, q, u, v = thomson_phase(-1.0, 0.3, -0.4)
    assert np.isclose(q, 0.3) and np.isclose(u, 0.4)


@pytest.mark.parametrize('q,u', [(0, 0), (1, 0), (0, 1), (0.6, -0.8), (-0.3, 0.2)])
def test_phase_function_normalized_for_any_polarization(q, u):
    '''<4 pi p> over the sphere must be 1 for any incoming polarization.'''
    mu, wmu = np.polynomial.legendre.leggauss(64)
    phi = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    M, PHI = np.meshgrid(mu, phi, indexing='ij')
    qp, up = rotate_stokes(q, u, PHI)
    f, *_ = thomson_phase(M, qp, up)
    avg = np.sum(wmu[:, None] * f) / len(phi) / 2.0
    assert np.isclose(avg, 1.0)


def test_scattered_p_never_exceeds_one(rng):
    c = rng.uniform(-1, 1, 10000)
    ang = rng.uniform(0, 2 * np.pi, 10000)
    P = np.sqrt(rng.random(10000))
    _, q, u, _ = thomson_phase(c, P * np.cos(ang), P * np.sin(ang))
    assert np.all(np.hypot(q, u) <= 1.0 + 1e-12)


# ---------------------------------------------------------------- sampling

def test_cos_theta_follows_thomson(rng):
    mu = sample_thomson_cos(rng, 200000)
    cdf = lambda x: 0.375 * (x + x ** 3 / 3.0 + 4.0 / 3.0)
    assert stats.kstest(mu, cdf).pvalue > 0.01


@pytest.mark.parametrize('c,q,u', [(0.0, 1.0, 0.0), (0.3, 0.0, 0.8), (-0.6, -0.5, 0.5), (0.9, 0.7, 0.0)])
def test_phi_moments(rng, c, q, u):
    '''For p(phi) ~ 1 + a(q cos2phi + u sin2phi): <cos2phi> = a q/2, <sin2phi> = a u/2.'''
    n = 400000
    phi = sample_thomson_phi(rng, np.full(n, c), np.full(n, q), np.full(n, u))
    a = (1 - c * c) / (1 + c * c)
    tol = 4.0 / np.sqrt(2 * n)
    assert abs(np.mean(np.cos(2 * phi)) - a * q / 2) < tol
    assert abs(np.mean(np.sin(2 * phi)) - a * u / 2) < tol


def test_scatter_es_postconditions(rng):
    n = 20000
    p = beam(n, random_unit(rng, n), q=0.4, u=-0.3)
    w0 = p['w'].copy()
    scatter_es(p, rng)
    assert np.allclose(np.linalg.norm(p['lmn'], axis=1), 1.0)
    assert np.allclose(np.linalg.norm(p['pol_ref'], axis=1), 1.0)
    assert np.allclose(dot(p['pol_ref'], p['lmn']), 0.0, atol=1e-12)
    assert np.all(np.hypot(p['stokes'][:, 0], p['stokes'][:, 1]) <= 1.0 + 1e-12)
    assert np.array_equal(p['w'], w0)


def test_unpolarized_scatter_is_azimuthally_symmetric(rng):
    '''Unpolarized light along +z: scattered directions are symmetric about z.'''
    n = 200000
    p = beam(n, [0, 0, 1.0])
    scatter_es(p, rng)
    az = np.arctan2(p['lmn'][:, 1], p['lmn'][:, 0])
    assert abs(np.mean(np.cos(2 * az))) < 4 / np.sqrt(n)
    assert abs(np.mean(np.sin(2 * az))) < 4 / np.sqrt(n)


# ---------------------------------------------------------------- sampling vs. peel-off

@pytest.mark.parametrize('q,u', [(0.0, 0.0), (0.8, 0.0), (0.3, -0.6)])
def test_sampled_scatters_agree_with_peel_off(rng, q, u):
    '''
    Scatter a polarized beam many times and collect the photons that leave
    within a small cone around n_obs.  The fraction in the cone must match
    f * (cone solid angle)/(4 pi), and their mean sky Stokes parameters must
    match the peel-off prediction.  This ties scatter_es and peel_es together.
    '''
    nphot = 2000000
    n_in = normalize(np.array([0.3, -0.5, 0.8]))
    e1 = photon.default_pol_ref(n_in)
    for angle, phase in [(90, 0.0), (40, 0.37), (150, 0.8)]:
        n_obs = observer_direction(angle, phase)
        N = sky_north(angle, phase)

        p = beam(1, n_in, q, u, e1)
        f_peel, q_peel, u_peel = peel_es(p, n_obs, N)

        p = beam(nphot, n_in, q, u, e1)
        scatter_es(p, rng)
        cone = np.radians(4.0)
        sel = dot(p['lmn'], n_obs) > np.cos(cone)
        frac_expected = f_peel[0] * (1 - np.cos(cone)) / 2
        nsel = sel.sum()
        assert abs(nsel / nphot - frac_expected) < 5 * np.sqrt(frac_expected / nphot)

        qs, us = to_sky(p['stokes'][sel, 0], p['stokes'][sel, 1], p['pol_ref'][sel], p['lmn'][sel], N)
        tol = 5 / np.sqrt(nsel) + 0.01        # noise plus variation across the cone
        assert abs(np.mean(qs) - q_peel[0]) < tol
        assert abs(np.mean(us) - u_peel[0]) < tol


def test_peel_unpolarized_matches_extract_weight(rng):
    '''For unpolarized photons the peel factor is the 0.75(1+cos^2) used in extract.c.'''
    n = 1000
    p = beam(n, random_unit(rng, n))
    n_obs = observer_direction(37, 0.2)
    f, _, _ = peel_es(p, n_obs, sky_north(37, 0.2))
    c = dot(p['lmn'], n_obs)
    assert np.allclose(f, 0.75 * (1 + c * c))


def test_peel_forward_and_back_scatter_are_finite():
    n_obs = observer_direction(90, 0)
    N = sky_north(90, 0)
    for n_in in (n_obs, -n_obs):
        p = beam(1, n_in, 0.5, 0.2)
        f, q, u = peel_es(p, n_obs, N)
        assert np.all(np.isfinite([f, q, u]))
        assert np.isclose(np.hypot(q, u), np.hypot(0.5, 0.2))


def test_sign_axis_scatter_gives_negative_q():
    '''
    Unpolarized light travelling along the symmetry axis, scattered to an
    edge-on observer: E-vector perpendicular to the projected axis, Q = -1
    (the bipolar sign).  Light travelling in the equatorial plane gives
    Q = +1 (the flattened-envelope sign).
    '''
    n_obs, N = observer_direction(90, 0), sky_north(90, 0)
    _, q, u = peel_es(beam(1, [0, 0, 1.0]), n_obs, N)
    assert np.isclose(q[0], -1.0) and np.isclose(u[0], 0.0, atol=1e-12)
    _, q, u = peel_es(beam(1, [0, 1.0, 0]), n_obs, N)
    assert np.isclose(q[0], 1.0) and np.isclose(u[0], 0.0, atol=1e-12)


# ---------------------------------------------------------------- Lorentz transform

def test_boost_keeps_reference_perpendicular_and_handed(rng):
    nphot = 500
    p = beam(nphot, random_unit(rng, nphot))
    p['freq'] = 1e15
    p['frame'] = photon.F_OBSERVER
    v = random_unit(rng, nphot) * rng.uniform(0.01, 0.6, nphot)[:, None] * C
    e2_old = np.cross(p['lmn'], p['pol_ref'])
    e2_new = boost_pol_ref(e2_old, p['lmn'], v)
    lorentz_transform(p, v)
    assert np.allclose(np.linalg.norm(p['lmn'], axis=1), 1.0)
    assert np.allclose(dot(p['pol_ref'], p['lmn']), 0.0, atol=1e-12)
    # the transformed second axis is still n' x e1': the basis is conformally carried
    assert np.allclose(np.cross(p['lmn'], p['pol_ref']), e2_new, atol=1e-10)


def test_boost_round_trip(rng):
    nphot = 200
    p = beam(nphot, random_unit(rng, nphot))
    p['freq'] = 1e15
    p0 = p.copy()
    v = random_unit(rng, nphot) * 0.3 * C
    lorentz_transform(p, v)
    lorentz_transform(p, -v)
    assert np.allclose(p['lmn'], p0['lmn'])
    assert np.allclose(np.abs(dot(p['pol_ref'], p0['pol_ref'])), 1.0)
    assert np.allclose(p['freq'], p0['freq'])


def test_linear_mode_leaves_reference_alone(rng):
    p = beam(10, random_unit(rng, 10))
    p['freq'] = 1e15
    ref = p['pol_ref'].copy()
    lmn = p['lmn'].copy()
    lorentz_transform(p, [1e8, 0, 0], rel_mode=REL_MODE_LINEAR)
    assert np.array_equal(p['pol_ref'], ref) and np.array_equal(p['lmn'], lmn)
