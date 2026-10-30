'''
Step 4 tests: ray tracing through Sirocco-style grids (sirocco_grid.py),
using synthetic grids so that no Sirocco run is needed.

Run with:  python -m pytest py_progs/polarization/tests
'''

import numpy as np
import pytest

from polarization.stokes import normalize
from polarization.sirocco_grid import (Domain, SiroccoModel, run, brown_mclean_grid, ds_to_cone, ds_to_plane,
                                       SIGMA_T, STAT_ESCAPE, STAT_INTERACT)
from polarization.toy_transport import Ellipsoid, run as toy_run

NSIG = 4.0
A, C = 1.0, 0.3


def cyl_ellipsoid(tau_eq, n=60):
    x = np.linspace(0, 1.4 * A, n)
    z = np.linspace(0, 1.4 * A, n)
    X, Z = np.meshgrid(0.5 * (x[:-1] + x[1:]), 0.5 * (z[:-1] + z[1:]), indexing='ij')
    ne = np.zeros((n, n))
    ne[:-1, :-1] = np.where((X / A) ** 2 + (Z / C) ** 2 < 1, tau_eq / (SIGMA_T * A), 0)
    return Domain('cylindrical', x, z, ne)


def polar_ellipsoid(tau_eq, nr=60, nth=60):
    r = np.linspace(0, 1.4 * A, nr)
    th = np.linspace(0, np.pi / 2, nth)
    R, T = np.meshgrid(0.5 * (r[:-1] + r[1:]), 0.5 * (th[:-1] + th[1:]), indexing='ij')
    ne = np.zeros((nr, nth))
    ne[:-1, :-1] = np.where((R * np.sin(T) / A) ** 2 + (R * np.cos(T) / C) ** 2 < 1, tau_eq / (SIGMA_T * A), 0)
    return Domain('polar', r, th, ne)


@pytest.fixture
def rng():
    return np.random.default_rng(53)


def test_cone_distances():
    x = np.array([[1.0, 0, 1.0], [1.0, 0, -2.0]])
    n = np.array([[0, 0, -1.0], [0, 0, 1.0]])
    assert np.allclose(ds_to_cone(x, n, np.pi / 2, 1e-12), ds_to_plane(x, n, 0.0, 1e-12))
    # from a point on the 45 degree cone, straight down to the other nappe
    assert np.isclose(ds_to_cone(x[:1], n[:1], np.pi / 4, 1e-12)[0], 2.0)


def test_read_master_tables(tmp_path):
    '''Parse small cylindrical and polar master tables in windsave2table format.'''
    rows = ['x z xcen zcen i j inwind ne']
    for i, xv in enumerate([0.0, 1.0, 2.0]):
        for j, zv in enumerate([0.0, 0.5, 1.0]):
            rows.append(f'{xv} {zv} {xv + 0.5} {zv + 0.25} {i} {j} {0 if i < 2 and j < 2 else -1} {10 * i + j + 1}')
    f = tmp_path / 'c.master.txt'
    f.write_text('\n'.join(rows) + '\n')
    d = Domain.from_master(str(f))
    assert d.coord == 'cylindrical'
    assert np.allclose(d.e1, [0, 1, 2]) and np.allclose(d.e2, [0, 0.5, 1])
    assert d.ne[1, 1] == 12 and d.ne[2, 2] == 0           # not in the wind -> 0

    rows = ['r theta r_cen theta_cen x z xcen zcen i j inwind ne']
    for i, rv in enumerate([0.0, 1.0]):
        for j, tv in enumerate([0.0, 45.0, 90.0]):
            rows.append(f'{rv} {tv} 0 0 0 0 0 0 {i} {j} 0 {i + j + 1}')
    f = tmp_path / 'p.master.txt'
    f.write_text('\n'.join(rows) + '\n')
    d = Domain.from_master(str(f))
    assert d.coord == 'polar'
    assert np.allclose(d.e2, np.radians([0, 45, 90]))       # degrees converted


def test_polar_and_cylindrical_grids_agree_on_optical_depth(rng):
    mc = SiroccoModel.from_domains([cyl_ellipsoid(1.0)], 2.0)
    mp = SiroccoModel.from_domains([polar_ellipsoid(1.0)], 2.0)
    x = normalize(rng.normal(size=(2000, 3))) * rng.uniform(0, 0.3, 2000)[:, None]
    n = normalize(rng.normal(size=(2000, 3)))
    t = [m.walk(x, n, np.full(2000, np.inf))[2] for m in (mc, mp)]
    assert abs(np.mean(t[1] - t[0])) < 0.005 and np.std(t[1] - t[0]) < 0.03


def test_ray_in_the_disk_plane_on_a_polar_grid():
    '''theta = 90 deg exactly must use the last cell above the plane (the bug fixed in step 4).'''
    mp = SiroccoModel.from_domains([polar_ellipsoid(1.0)], 2.0)
    tau = mp.tau_to_observer(np.zeros((1, 3)), np.array([1.0, 0, 0]))[0]
    assert abs(tau - 1.0) < 0.03


def test_interaction_probability(rng):
    m = SiroccoModel.from_domains([cyl_ellipsoid(1.0)], 2.0)
    x = np.zeros((20000, 3))
    n = normalize(rng.normal(size=(20000, 3)))
    tau = m.walk(x, n, np.full(20000, np.inf))[2]
    _, st, _ = m.walk(x, n, -np.log(rng.random(20000)))
    assert abs(np.mean(st == STAT_INTERACT) - np.mean(1 - np.exp(-tau))) < 0.015


@pytest.mark.parametrize('make', [cyl_ellipsoid, polar_ellipsoid])
def test_grid_ellipsoid_matches_smooth_ellipsoid(rng, make):
    m = SiroccoModel.from_domains([make(1.0)], 2.0, angles=[90])
    q, u, sq, su = run(m, 100000, rng)[0].result()
    q2, _, sq2, _ = toy_run(Ellipsoid(A, C, 1.0 / A), 100000, [90], rng)[0][0].result()
    assert abs(q - q2) < NSIG * np.hypot(sq, sq2) + 0.03 * abs(q2)
    assert abs(u) < NSIG * su


def test_spherical_grid_has_no_net_polarization(rng):
    r = np.linspace(0.2, 1.0, 30)
    ne = np.full(30, 1.0 / (SIGMA_T * 0.8))
    ne[-1] = 0
    m = SiroccoModel.from_domains([Domain('spherical', r, None, ne)], 1.5, angles=[45, 90])
    for s in run(m, 50000, rng):
        q, u, sq, su = s.result()
        assert abs(q) < NSIG * sq and abs(u) < NSIG * su
    assert np.isclose(brown_mclean_grid(m, 90)[2], 1.0 / 3.0)


def test_wind_emission_positions_and_polarization(rng):
    '''
    Wind photons are placed uniformly in volume in proportion to the cell
    emission, and an emitting sphere seen from any angle is unpolarized.
    '''
    for make in (cyl_ellipsoid, polar_ellipsoid):
        dm = make(0.5)
        # uniform emissivity in the ellipsoid: the emission of a cell is its volume
        vol = dm.cells()[2].reshape(dm.ne[:-1, :-1].shape)
        dm.emis[:-1, :-1] = np.where(dm.ne[:-1, :-1] > 0, vol, 0.0)
        m = SiroccoModel.from_domains([dm], 2.0, angles=[90])
        m.sources = {'wind': 1.0}
        m.point_source = False
        x = m._wind_positions(rng, 40000)
        # uniform within the ellipsoid of uniform emissivity: <(rho/a)^2 + (z/c)^2> = 3/5
        e2 = (np.hypot(x[:, 0], x[:, 1]) / A) ** 2 + (x[:, 2] / C) ** 2
        assert abs(np.mean(e2) - 0.6) < 0.02
        assert abs(np.mean(x[:, 2] > 0) - 0.5) < 0.01

    r = np.linspace(0.2, 1.0, 30)
    ne = np.full(30, 1.0 / (SIGMA_T * 0.8))
    ne[-1] = 0
    emis = np.append(np.diff(r ** 3), 0.0)
    m = SiroccoModel.from_domains([Domain('spherical', r, None, ne, emis)], 1.5, angles=[45, 90])
    m.sources = {'wind': 1.0}
    m.point_source = False
    for s in run(m, 50000, rng):
        q, u, sq, su = s.result()
        assert abs(q) < NSIG * sq and abs(u) < NSIG * su
