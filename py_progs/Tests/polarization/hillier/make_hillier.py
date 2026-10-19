#!/usr/bin/env python
'''
Make the Hillier (1994, A&A 289, 492, Sect. 9) test used by Bulla+15 Sect. 3.1:
a point-like star inside a detached prolate electron-scattering shell,

    chi_es = chi0 (R_min/r)^4 (1 + 10 cos^2 beta),   R_min < r < 30 R_min,

with R_* = 0.2 R_min.  Hillier's lengths are in units where R_min = 2, so
chi0 is an opacity per (R_min/2).  tau_av = 2.888 chi0, gamma = 7/13 and the
optically thin (Brown & McLean) result is q = -66.6 chi0 sin^2 i per cent.

Imported polar grid (0-90 deg, mirrored by Sirocco), pure H at 1e6 K so
electron scattering is the only significant opacity (ff/es < 1e-4 for
R_min = 1e13 cm), no wind emission, 0 ionization cycles.

Writes root.pf and root.import.txt.

usage: make_hillier.py root chi0 [-nscat] [-angle A] [-nphot N] [-ncycles C]

    without -nscat: observers at 11.25, 22.5, ... 90 deg, all photons
    with -nscat:    20 observers at the angle given by -angle (default 90),
                    selecting photons scattered exactly 0..18 times and
                    >= 19 times, for I(a) = sum a^j I_j (Hillier App. A)
    -nphot, -ncycles: photons per cycle and spectral cycles (default 8000000, 2)
'''
import sys
import numpy as np

SIGMA_T = 6.6524587e-25
M_U = 1.66053907e-24
MU_H = 1.00794


def doit(root, chi0, nscat=False, angle=90.0, nphot=8000000, ncycles=2):
    '''
    Write root.import.txt and root.pf for a Hillier shell with opacity
    chi0; with nscat, 20 observers at angle selected by the number of
    scatters
    '''

    r_min = 1e13
    rstar, r_out = 0.2 * r_min, 30 * r_min
    L = r_min / 2                                   # Hillier's unit of length
    nr, nth = 80, 45
    re = np.geomspace(r_min, r_out, nr + 1)
    re = np.append(re, re[-1] * re[-1] / re[-2])    # guard row
    te = np.linspace(0, 90, nth + 1)
    te = np.append(te, 90 + 90 / nth)               # guard column
    n0 = chi0 / (SIGMA_T * L)                       # n_e at r_min, beta = 90 deg (equator)


    def ne_cell(i, j):
        '''cell average: exact for r^-4 radially, (1 + 10 mu^2) averaged over the cell in mu'''
        r1, r2 = re[i], re[i + 1]
        rad = r_min ** 4 * (r1 ** -3 - r2 ** -3) / (3 * (r2 - r1))
        mu1, mu2 = np.cos(np.radians(te[j + 1])), np.cos(np.radians(te[j]))
        ang = 1 + 10 * (mu2 ** 3 - mu1 ** 3) / (3 * (mu2 - mu1))
        return n0 * rad * ang


    v = 1e7
    with open(root + '.import.txt', 'w') as f:
        f.write('# i j inwind r theta v_x v_y v_z rho t_e t_r\n')
        for i in range(len(re)):
            for j in range(len(te)):
                inwind = 0 if (i < nr and j < nth) else -1
                rho = ne_cell(i, j) * MU_H * M_U if inwind == 0 else 0.0
                th = np.radians(te[j])
                f.write(f'{i:3d} {j:3d} {inwind:3d} {re[i]:.8e} {te[j]:.6f} {v*np.sin(th):.6e} 0.0 {v*np.cos(th):.6e} {rho:.6e} 1e6 1e6\n')

    # tau_av and gamma from the cells, as Hillier eqs (2), (3)
    S = S2 = 0.0
    for j in range(nth):
        mu1, mu2 = np.cos(np.radians(te[j + 1])), np.cos(np.radians(te[j]))
        taur = sum(ne_cell(i, j) * (re[i + 1] - re[i]) for i in range(nr)) * SIGMA_T
        S += 2 * taur * (mu2 - mu1)
        S2 += 2 * taur * (mu2 ** 3 - mu1 ** 3) / 3
    tau_av, gamma = 0.5 * S, S2 / S
    print(f'{root}: chi0 = {chi0}, n0 = {n0:.4e}, tau_av = {tau_av:.4f} (2.888 chi0 = {2.888*chi0:.4f}), '
          f'gamma = {gamma:.4f} (7/13 = {7/13:.4f}), thin q(90) = {100*0.375*tau_av*(1-3*gamma):+.3f}%')

    if nscat:
        angles = [angle] * 20
        sel = list(range(19)) + [-19]
    else:
        angles = [11.25 * k for k in range(1, 9)]
        sel = None

    obs = ''.join(f'Spectrum.angle(0=pole)                     {a}\n' for a in angles)
    if sel is None:
        scat = '@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)                   no\n'
    else:
        scat = '@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)                   yes\n'
        scat += ''.join(f'@Spectrum.select_scatters                   {s}\n' for s in sel)

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
Photons_per_cycle                          {nphot}
Ionization_cycles                          0
Spectrum_cycles                            {ncycles}
Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,LTE_iterate,fixed,matrix_bb,matrix_pow,matrix_est,matrix_multishot)               lte_te
Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping)          escape_prob
Wind.radiation(yes,no)                           no
Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)               absorb
Wind_heating.extra_processes(none,adiabatic,nonthermal,both)                 none
Atomic_data                                h_only_master.dat
@Diag.write_atomicdata(yes,no)                   no
Wind.filling_factor(1=smooth,<1=clumped)   1
Central_object.rad_type_in_final_spectrum(bb,models,uniform,mono)                   bb
Spectrum.nwave                             1000
Spectrum.wavemin(Angstroms)                3500
Spectrum.wavemax(Angstroms)                7500
Spectrum.no_observers                      {len(angles)}
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


def steer(argv):
    '''
    Process the command line
    '''

    words = []
    nscat = False
    angle = 90.0
    nphot = 8000000
    ncycles = 2

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-nscat':
            nscat = True
        elif argv[i] == '-angle':
            i += 1
            angle = float(argv[i])
        elif argv[i] == '-nphot':
            i += 1
            nphot = int(argv[i])
        elif argv[i] == '-ncycles':
            i += 1
            ncycles = int(argv[i])
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            words.append(argv[i])
        i += 1

    if len(words) != 2:
        print(__doc__)
        return

    doit(words[0], float(words[1]), nscat, angle, nphot, ncycles)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
