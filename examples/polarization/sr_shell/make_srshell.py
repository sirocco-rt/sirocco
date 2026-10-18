#!/usr/bin/env python
'''
The Hillier (1994) prolate electron-scattering shell, made optically
thin and given a constant radial bulk velocity beta c, to test the special
relativity of polarized electron scattering in Sirocco.

    n_e = n0 (R_min/r)^4 (1 + 10 cos^2 theta)    (lab frame; Sirocco uses n_e' = n_e / gamma)
    R_min < r < 30 R_min,  v = beta c r_hat

The star is made small (R_* = 0.01 R_min), so the incident light is radial and
is not aberrated in the comoving frame.  For each observer angle there are two
observers: photons scattered 0 times (the direct light) and exactly once.  The
spectral range is wide, so the bolometric I, Q of the single-scattered light can
be compared with sr_single.py.

Writes root.pf, root.import.txt, and root.cells.npz (the grid and n_e, which
sr_single.py and compare_sr.py read).

usage: make_srshell.py root chi0 beta [-sphere] [-nphot N] [-ncycles C]
       (-sphere drops the (1 + 10 cos^2) factor; default 4000000 photons, 2 cycles)
'''
import argparse
import numpy as np

SIGMA_T = 6.6524587e-25
M_U = 1.66053907e-24
MU_H = 1.00794
C = 2.99792458e10

ap = argparse.ArgumentParser()
ap.add_argument('root')
ap.add_argument('chi0', type=float)
ap.add_argument('beta', type=float)
ap.add_argument('-sphere', action='store_true', help='drop the (1 + 10 cos^2) factor')
ap.add_argument('-nphot', type=int, default=4000000)
ap.add_argument('-ncycles', type=int, default=2)
args = ap.parse_args()
root, chi0, beta = args.root, args.chi0, args.beta

r_min = 1e13
rstar, r_out = 0.01 * r_min, 30 * r_min
L = r_min / 2                                   # Hillier's unit of length
nr, nth = 80, 45
re = np.geomspace(r_min, r_out, nr + 1)
re = np.append(re, re[-1] * re[-1] / re[-2])    # guard row
te = np.linspace(0, 90, nth + 1)
te = np.append(te, 90 + 90 / nth)               # guard column
n0 = chi0 / (SIGMA_T * L)
aniso = 0.0 if args.sphere else 10.0


def ne_cell(i, j):
    '''cell average: exact for r^-4 radially, (1 + 10 mu^2) averaged over the cell in mu'''
    r1, r2 = re[i], re[i + 1]
    rad = r_min ** 4 * (r1 ** -3 - r2 ** -3) / (3 * (r2 - r1))
    mu1, mu2 = np.cos(np.radians(te[j + 1])), np.cos(np.radians(te[j]))
    ang = 1 + aniso * (mu2 ** 3 - mu1 ** 3) / (3 * (mu2 - mu1))
    return n0 * rad * ang


# the cell densities are also written out for sr_single.py
ne = np.array([[ne_cell(i, j) for j in range(nth)] for i in range(nr)])
np.savez(root + '.cells.npz', re=re[:nr + 1], te=te[:nth + 1], ne=ne, beta=beta)

v = beta * C
with open(root + '.import.txt', 'w') as f:
    f.write('# i j inwind r theta v_x v_y v_z rho t_e t_r\n')
    for i in range(len(re)):
        for j in range(len(te)):
            inwind = 0 if (i < nr and j < nth) else -1
            rho = ne[i, j] * MU_H * M_U if inwind == 0 else 0.0
            th = np.radians(te[j])
            f.write(f'{i:3d} {j:3d} {inwind:3d} {re[i]:.8e} {te[j]:.6f} {v*np.sin(th):.6e} 0.0 {v*np.cos(th):.6e} {rho:.6e} 1e5 1e5\n')

angles = [22.5, 45.0, 67.5, 90.0]
obs = ''.join(f'Spectrum.angle(0=pole)                     {a}\n' for a in angles for s in (0, 1))
scat = '@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)                   yes\n'
scat += ''.join(f'@Spectrum.select_scatters                   {s}\n' for a in angles for s in (0, 1))
print(f'{root}: chi0 = {chi0}, beta = {beta}, sphere = {args.sphere}, n0 = {n0:.4e}')

pf = f'''System_type(star,cv,bh,agn,previous)                 star
Central_object.mass(msol)                  1
Central_object.radius(cm)                  {rstar:.3e}
Central_object.radiation(yes,no)                  yes
Central_object.rad_type_to_make_wind(bb,models)                   bb
Central_object.temp                        10000
Disk.type(none,flat,vertically.extended,rmin>central.obj.rad)                 none
Boundary_layer.radiation(yes,no)                   no
Wind.number_of_components                  1
Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)             imported
Wind.coord_system(spherical,cylindrical,polar,cyl_var)                polar
Wind.model2import                          {root}.import.txt
@Diag.adjust_grid(yes,no)                        no
Photons_per_cycle                          {args.nphot}
Ionization_cycles                          0
Spectrum_cycles                            {args.ncycles}
Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,LTE_iterate,fixed,matrix_bb,matrix_pow,matrix_est,matrix_multishot)               lte_te
Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping)          escape_prob
Wind.radiation(yes,no)                           no
Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)               absorb
Wind_heating.extra_processes(none,adiabatic,nonthermal,both)                 none
Atomic_data                                h_only_master.dat
@Diag.write_atomicdata(yes,no)                   no
Wind.filling_factor(1=smooth,<1=clumped)   1
Central_object.rad_type_in_final_spectrum(bb,models,uniform,mono)                   bb
Spectrum.nwave                             10000
Spectrum.wavemin(Angstroms)                300
Spectrum.wavemax(Angstroms)                200000
Spectrum.no_observers                      {2*len(angles)}
{obs}Spectrum.live_or_die(live.or.die,extract)              extract
{scat}@Spectrum.select_photons_by_position(yes,no)                   no
@Spectrum.polarization(yes,no)                   yes
Spectrum.type(flambda,fnu,basic)              flambda
@Diag.use_standard_care_factors(yes,no)                  yes
Reverb.type(none,photon,wind,matom)                 none
Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)               t_star
@Diag.extra(yes,no)                              no
'''
open(root + '.pf', 'w').write(pf)
