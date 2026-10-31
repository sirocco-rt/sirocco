#!/usr/bin/env python
'''
Create the density structure of an AGN torus, on a polar grid that
Sirocco can import.


Synopsis:

    Write root.import.txt, a polar (r-theta) grid containing a torus
    of the kind used by Vander Meulen et al. (2024, A&A 689;
    arXiv:2409.02986) for their X-ray polarization models.

Command line usage:

    usage::

        make_torus.py root.model.pf

    where root.model.pf is a parameter file (read with kpar) that
    describes the torus, and the grid is written to root.import.txt.
    If the file does not exist, one is asked for each parameter, with a
    suggested value, and the file is created for future use.  A name
    given without .pf has .model.pf added.  The .model.pf extension keeps
    the file apart from the Sirocco parameter file of the run, root.pf.

    The parameters, in the order they are read, are

        Model               the density structure; at present only wedge
        R_min, R_max        the inner and outer radius of the torus (cm)
        Theta_t             the half-thickness of the torus, in degrees
                            above and below the equatorial plane
        N_H                 the hydrogen column density (cm**-2) along a
                            line of sight through the torus from the
                            central source
        No.Radial.Bins      the number of radial cells across the torus
        No.Theta.Bins       the number of theta cells across the torus,
                            from the equator to Theta_t (Sirocco mirrors
                            the grid below the equatorial plane)
        Radial.Spacing      log or linear
        M_bh                the mass of the black hole (Msun), which sets
                            the rotation of the torus
        T_init              the initial electron and radiation
                            temperature (K)

Description:

    Version 1 (Model wedge) is the "wedge torus" of Vander Meulen et al.:
    gas of constant density between rmin and rmax, within theta_t of the
    equatorial plane, around a central point source.  The covering factor
    of the torus, as seen from the source, is sin(theta_t).  Vander Meulen
    et al. write the covering factor as cos(theta), where theta =
    90 - theta_t is the half-opening angle of the torus measured from the
    pole.

    Because the density is constant and the torus is radial, every line
    of sight from the source through the torus has the same column
    density, so

        n_H = nh / (rmax - rmin)
        rho = 1.4 m_p n_H

    The factor 1.4 is approximate, for gas of solar composition (H and
    He).  Sirocco calculates its own conversion between rho and n_H from
    the abundances in the atomic data it is given, so the column density
    in a run will differ from nh by a few per cent.  The electron
    scattering optical depth of the torus is roughly 1.2 nh sigma_T if
    the gas is fully ionized.

    The torus rotates.  A wind that is exactly static causes Sirocco to
    abort (issue #1208), so it is given the circular speed about the
    black hole, v_phi = sqrt(G M / r), at every point.  This has little
    effect on the continuum.

    The grid is polar, with r and theta given at the inner edges of the
    cells, and the velocities in Cartesian coordinates, as described in
    the Sirocco documentation on importing models.  In theta there are
    No.Theta.Bins equal cells between 90 - theta_t and 90 degrees (theta
    is measured from the pole), and the empty region between the pole and
    the torus is a single cell.  In r there are No.Radial.Bins cells
    between rmin and rmax.  A row and a column of ghost cells, not in the wind, are
    added beyond rmax and beyond 90 degrees.

    The routine creates only the density structure.  A Sirocco .pf file
    that imports it, with Wind.coord_system polar, has to be made
    separately.

Primary routines:

    doit - creates the grid and writes the import file
    get_inputs - reads (or asks for) the parameters
    steer - processes the command line

Notes:

    Other versions of the torus can be added as further density
    functions, selected with Model.

    doit can also be called directly from python, with the parameters
    as keywords.

History:

    261019 ksl Coding begun

'''

import os
import sys
import numpy as np
from astropy.table import Table
sys.path.insert(0, os.path.join(os.environ.get('SIROCCO', ''), 'py_progs'))
import kpar


MPROT = 1.672621923e-24
SIGMA_T = 6.6524587e-25
G = 6.674e-8
MSOL = 1.989e33
RHO_PER_NH = 1.4 * MPROT        # approximate, solar composition
NE_PER_NH = 1.2                 # fully ionized gas, solar composition

MODELS = ['wedge']


def make_grid(rmin, rmax, theta_t, nr=30, ntheta=30, linear=False):
    '''
    Return the radial and angular cell edges (cm, degrees), including
    the ghost row and column, and the number of theta cells up to the
    equator
    '''

    if linear:
        r_edges = np.linspace(rmin, rmax, nr + 1)
        r_edges = np.append(r_edges, 2 * r_edges[-1] - r_edges[-2])
    else:
        r_edges = np.geomspace(rmin, rmax, nr + 1)
        r_edges = np.append(r_edges, r_edges[-1] * r_edges[-1] / r_edges[-2])

    t_edges = np.linspace(90. - theta_t, 90., ntheta + 1)
    if theta_t < 90.:
        t_edges = np.insert(t_edges, 0, 0.0)
    nth = len(t_edges) - 1
    t_edges = np.append(t_edges, t_edges[-1] + theta_t / ntheta)

    return r_edges, t_edges, nth


def wedge(r, theta, rmin, rmax, theta_t, nh):
    '''
    The density (g/cm**3) of a wedge torus at the centre of a cell at
    r (cm) and theta (degrees from the pole)
    '''

    if rmin < r < rmax and theta > 90. - theta_t and theta < 90.:
        return RHO_PER_NH * nh / (rmax - rmin)
    return 0.0


def doit(root='torus', rmin=1e17, rmax=1e18, theta_t=30., nh=1e24, nr=30, ntheta=30, linear=False,
         mbh=1e8, t_init=1e5, model='wedge'):
    '''
    Create the grid for a torus and write it to root.import.txt
    '''

    if model not in MODELS:
        print('Error: unknown model %s; the choices are %s' % (model, MODELS))
        return

    if rmin >= rmax or theta_t <= 0 or theta_t > 90:
        print('Error: need rmin < rmax and 0 < theta_t <= 90')
        return

    r_edges, t_edges, nth = make_grid(rmin, rmax, theta_t, nr, ntheta, linear)

    ii = []
    jj = []
    inwind = []
    r = []
    theta = []
    v_x = []
    v_y = []
    v_z = []
    rho = []

    for i in range(len(r_edges)):
        for j in range(len(t_edges)):
            if i < nr and j < nth:
                r_cen = 0.5 * (r_edges[i] + r_edges[i + 1])
                t_cen = 0.5 * (t_edges[j] + t_edges[j + 1])
                x = wedge(r_cen, t_cen, rmin, rmax, theta_t, nh)
            else:
                x = 0.0
            ii.append(i)
            jj.append(j)
            inwind.append(0 if x > 0 else -1)
            r.append(r_edges[i])
            theta.append(t_edges[j])
            v_x.append(0.0)
            v_y.append(np.sqrt(G * mbh * MSOL / r_edges[i]))
            v_z.append(0.0)
            rho.append(x)

    xtab = Table([ii, jj, inwind, r, theta, v_x, v_y, v_z, rho],
                 names=['i', 'j', 'inwind', 'r', 'theta', 'v_x', 'v_y', 'v_z', 'rho'])
    xtab['t_e'] = t_init
    xtab['t_r'] = t_init

    xtab['r'].format = '.6e'
    xtab['theta'].format = '.4f'
    xtab['v_x'].format = xtab['v_y'].format = xtab['v_z'].format = '.4e'
    xtab['rho'].format = '.4e'
    xtab['t_e'].format = xtab['t_r'].format = '.2e'

    outfile = root + '.import.txt'
    xtab.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)

    n_h = nh / (rmax - rmin)
    print('%s: %s torus, r = %.3e - %.3e cm, within %g deg of the equator (covering factor %.3f)'
          % (outfile, model, rmin, rmax, theta_t, np.sin(np.radians(theta_t))))
    print('  N_H = %.3e cm**-2, n_H = %.3e cm**-3, rho = %.3e g/cm**3, tau_es ~ %.3f'
          % (nh, n_h, RHO_PER_NH * n_h, NE_PER_NH * nh * SIGMA_T))
    print('  grid: %d x %d cells in the torus (%s in r), %d x %d with the empty and ghost cells'
          % (nr, ntheta, 'linear' if linear else 'logarithmic', len(r_edges), len(t_edges)))

    return xtab


def get_par(name, hint):
    '''
    Read one parameter with kpar, printing a hint first if the
    parameter is not in the file, so it will be asked for
    '''

    if name not in [one[0].split('(')[0] for one in kpar.xxx]:
        print(hint)
    return kpar.rdpar(name)


def get_inputs(pf):
    '''
    Read the parameters of the torus from pf, asking for any that are
    missing, and write them all back to pf
    '''

    kpar.opar(pf)

    x = {}
    x['model'] = str(get_par('Model', 'The density structure, one of %s, e.g. wedge' % MODELS))
    x['rmin'] = float(get_par('R_min', 'The inner radius of the torus (cm), e.g. 1e17'))
    x['rmax'] = float(get_par('R_max', 'The outer radius of the torus (cm), e.g. 1e18'))
    x['theta_t'] = float(get_par('Theta_t', 'The half-thickness of the torus, degrees above and below the equator, '
                                 'e.g. 30'))
    x['nh'] = float(get_par('N_H', 'The H column density through the torus (cm**-2), e.g. 1e24'))
    x['nr'] = int(get_par('No.Radial.Bins', 'The number of radial cells across the torus, e.g. 30'))
    x['ntheta'] = int(get_par('No.Theta.Bins', 'The number of theta cells from the equator to Theta_t, e.g. 30'))
    spacing = str(get_par('Radial.Spacing', 'The spacing of the radial cells, log or linear, e.g. log'))
    x['mbh'] = float(get_par('M_bh', 'The mass of the black hole (Msun), e.g. 1e8'))
    x['t_init'] = float(get_par('T_init', 'The initial temperature (K), e.g. 1e5'))

    kpar.cpar(pf)

    x['linear'] = spacing.lower().startswith('lin')

    return x


def steer(argv):
    '''
    Process the command line
    '''

    pf = argv[1]
    if not pf.endswith('.pf'):
        pf = pf + '.model.pf'

    root = os.path.basename(pf)
    if root.endswith('.model.pf'):
        root = root[:-len('.model.pf')]
    else:
        root = root[:-len('.pf')]

    x = get_inputs(pf)
    doit(root, **x)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
