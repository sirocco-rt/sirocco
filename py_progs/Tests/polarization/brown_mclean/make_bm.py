#!/usr/bin/env python
'''
Make a point-source Brown & McLean test for the C polarization: an imported
polar grid with n_e = n0 (r_in/r)^2 sin^2(theta) between r_in and r_out,
pure H at 1e6 K (so electron scattering is the only opacity), a small
central star and no wind emission.  Writes ROOT.import.txt and ROOT.pf and
prints the Brown & McLean prediction computed from the same cells.

usage: make_bm.py root tau_eq [-nphot N] [-ncycles C]     (default 8000000 photons, 2 spectral cycles)
'''
import sys
import numpy as np

SIGMA_T = 6.6524587e-25
M_U = 1.66053907e-24
MU_H = 1.00794


def doit(root, tau_eq, nphot=8000000, ncycles=2):
    '''
    Write root.import.txt and root.pf for an envelope with an equatorial
    electron-scattering optical depth tau_eq, and print the Brown & McLean
    prediction for it
    '''

    rstar, r_in, r_out = 1e10, 1e11, 1e12
    nr, nth = 40, 30
    re = np.geomspace(r_in, r_out, nr + 1)
    re = np.append(re, re[-1] * re[-1] / re[-2])          # guard row
    te = np.linspace(0, 90, nth + 1)
    te = np.append(te, 90 + 90 / nth)                       # guard column
    n0 = tau_eq / (SIGMA_T * r_in * (1 - r_in / r_out))     # tau along the equator = tau_eq

    def ne_cell(i, j):
        # radial integral exact for r^-2: use n(r_c) with r_c^2 = r_i r_{i+1}; sin^2 averaged over the cell
        t1, t2 = np.radians(te[j]), np.radians(te[j + 1])
        mu1, mu2 = np.cos(t2), np.cos(t1)
        sin2 = 1 - (mu2 ** 3 - mu1 ** 3) / (3 * (mu2 - mu1))
        return n0 * r_in ** 2 / (re[i] * re[i + 1]) * sin2

    v = 1e7
    with open(root + '.import.txt', 'w') as f:
        f.write('# i j inwind r theta v_x v_y v_z rho t_e t_r\n')
        for i in range(len(re)):
            for j in range(len(te)):
                inwind = 0 if (i < nr and j < nth) else -1
                rho = ne_cell(i, j) * MU_H * M_U if inwind == 0 else 0.0
                th = np.radians(te[j])
                f.write(f'{i:3d} {j:3d} {inwind:3d} {re[i]:.8e} {te[j]:.6f} {v*np.sin(th):.6e} 0.0 {v*np.cos(th):.6e} {rho:.6e} 1e6 1e6\n')

    # Brown & McLean from the same cells (both hemispheres)
    S = S2 = 0.0
    for j in range(nth):
        t1, t2 = np.radians(te[j]), np.radians(te[j + 1])
        mu1, mu2 = np.cos(t2), np.cos(t1)
        taur = sum(ne_cell(i, j) * (re[i + 1] - re[i]) for i in range(nr)) * SIGMA_T
        S += 2 * taur * (mu2 - mu1)
        S2 += 2 * taur * (mu2 ** 3 - mu1 ** 3) / 3
    tau0, gamma = 3 / 16 * S, S2 / S
    print(f'{root}: tau_eq = {tau_eq}, n0 = {n0:.4e}, tau0 = {tau0:.5f}, gamma = {gamma:.5f}')
    for a in (10, 30, 50, 70, 90):
        print(f'  i = {a:2d}: Brown & McLean q = {100*tau0*(1-3*gamma)*np.sin(np.radians(a))**2:+.4f}%')

    pf = f'''System_type(star,cv,bh,agn,previous)                 star

### Parameters for the Central Object
Central_object.mass(msol)                  1
Central_object.radius(cm)                  {rstar:.3e}
Central_object.radiation(yes,no)                  yes
Central_object.rad_type_to_make_wind(bb,models)                   bb
Central_object.temp                        10000

### Parameters for the Disk (if there is one)
Disk.type(none,flat,vertically.extended,rmin>central.obj.rad)                 none

### Parameters for Boundary Layer or the compact object in an X-ray Binary or AGN
Boundary_layer.radiation(yes,no)                   no

### Parameters describing the various winds or coronae in the system
Wind.number_of_components                  1
Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)             imported
Wind.coord_system(spherical,cylindrical,polar,cyl_var)                polar
Wind.model2import                          {root}.import.txt
@Diag.adjust_grid(yes,no)                        no

### Parameters associated with photon number, cycles,ionization and radiative transfer options
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

### Parameters for Domain 0
Wind.filling_factor(1=smooth,<1=clumped)   1

### Parameters defining the spectra seen by observers

Central_object.rad_type_in_final_spectrum(bb,models,uniform,mono)                   bb

### The minimum and maximum wavelengths in the final spectra and the number of wavelength bins
Spectrum.nwave                             1000
Spectrum.wavemin(Angstroms)                3500
Spectrum.wavemax(Angstroms)                7500

### The observers and their location relative to the system
Spectrum.no_observers                      5
Spectrum.angle(0=pole)                     10
Spectrum.angle(0=pole)                     30
Spectrum.angle(0=pole)                     50
Spectrum.angle(0=pole)                     70
Spectrum.angle(0=pole)                     90
Spectrum.live_or_die(live.or.die,extract)              extract
@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)                   no
@Spectrum.select_photons_by_position(yes,no)                   no
@Spectrum.polarization(yes,no)                   yes
Spectrum.type(flambda,fnu,basic)              flambda
@Diag.use_standard_care_factors(yes,no)                  yes

### Parameters for Reverberation Modeling (if needed)
Reverb.type(none,photon,wind,matom)                 none

### Other parameters
Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)               t_star
@Diag.extra(yes,no)                              no
'''
    open(root + '.pf', 'w').write(pf)


def steer(argv):
    '''
    Process the command line
    '''

    words = []
    nphot = 8000000
    ncycles = 2

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
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

    doit(words[0], float(words[1]), nphot, ncycles)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
