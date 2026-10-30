#!/usr/bin/env python
'''
Semi-real polarization examples: a CV-like system with an imported radial
outflow confined to a range of polar angles, e.g. a polar cone (0-20 deg,
pbl_es) or a wedge close to the disk (60-88 deg, eq_es).

The defaults are those of pbl_es and eq_es (2026-10-17), a pure
electron-scattering test that converges:

    make_cvcone.py pbl_es 0 20
    make_cvcone.py eq_es 60 88
    sirocco -d -p pbl_es          (-p: photons per ionization cycle rise
                                   logarithmically to Photons_per_cycle)

i.e. constant density, radial tau_es = 1, a point-like central source
(boundary layer, L = 1e37 erg/s, T = 1e5 K; the disk is present but does
not radiate, and absorbs), simple atoms and no wind emission, 5e6 photons,
20 ionization and 3 spectral cycles, observers at 10, 45 and 80 deg, each
split into all/unscattered/scattered photons, 1000-3000 A, t_init 40 kK.
(A CV-like source, L ~ 5e34, cannot photoionize a tau_es ~ 1 wind out to
1e11 cm; the outer wind then collapses to the temperature floor.)

The system is a CV, so Sirocco includes a secondary star, which hides
from the observers any light whose path to them crosses its Roche lobe.
With the usual period of a few hours the lobe lies inside the wind
(L1 at 4e10 cm for 5.57 hr) and changes the polarization, which a
calculation without a secondary does not reproduce (seen 261019; the
wind itself is not affected).  The period is therefore set to 300 hr,
which puts the secondary at about 1.1e12 cm, well outside the wind.

    rho(r) = A (r_in / r)^n          (-rho_index n; default 0, constant density)
    rho(r) = A / (r^2 v(r))          (-rho_index wind; mass conserving)
    v(r)   = v0 + (v_inf - v0) (1 - r_in / r)^beta,   radial

inside theta_min < theta < theta_max, 0 outside.  A is set so that the
electron-scattering optical depth along a radial path through the wind,
from r_in to r_out, is tau_es, assuming the gas is fully ionized (n_e =
1.19 rho / m_p for Sirocco's standard abundances; Sirocco uses its own n_e
once the ionization is calculated).  Free-free and bound-free opacity go
as n_e^2, electron scattering as n_e, so for a given tau_es a flat profile
keeps the wind dominated by scattering; with 1/(r^2 v) and tau_es ~ 1 the
dense base makes the continuum optically thick.

The grid is polar: nr logarithmic radial cells from r_in to r_out, and
dtheta (2-degree) cells in theta across the wind; each empty region (e.g.
20-90 deg for a polar cone) is a single cell, since the grid must still
cover 0-90 deg (-full_theta: dtheta cells from 0 to 90 deg).  Sirocco
mirrors the grid below the disk.  Velocities are given at cell corners and
densities at cell centres (cell averages), which is Sirocco's convention.

Writes root.import.txt and root.pf.

usage: make_cvcone.py root theta_min theta_max [-tau TAU] [-rho_index n|wind]
                      [-bl L T | -disk] [-macro] [-no_split] [-nphot N]
                      [-ion_cycles C] [-spec_cycles C] [-angles a1 a2 ...]
                      [-wmin W] [-wmax W] [-nwave N] [-t_init T] ...

-bl L T: the central source is a boundary layer of luminosity L and
temperature T; the disk does not radiate, but is still there (and absorbs).
-disk: instead, a radiating standard disk (mdot 1e-8) and white dwarf, as in
examples/regress/cv.pf.
-macro: macro atoms for H and He (as in examples/regress/agn_macro.pf) and
wind emission, instead of simple atoms and Wind.radiation no (which is only
honoured without macro atoms).
-no_split: one observer per angle, instead of four: all photons,
unscattered photons (direct light, and light emitted by the wind), photons
scattered exactly once, and photons scattered at least once (electron or
resonance scattering).  The single-scatter spectra (_sc:1) are the ones to
compare with a single-scattering calculation.
'''
import sys
from types import SimpleNamespace
import numpy as np

SIGMA_T = 6.6524587e-25
MPROT = 1.672621923e-24
NE_PER_RHO = 1.19 / MPROT       # fully ionized, Sirocco's standard abundances
R_WD = 7e8


def doit(root, theta_min, theta_max, tau=1.0, r_in=4 * R_WD, r_out=1e11, v0=1e8, v_inf=3e8, beta=1.0,
         rho_index='0', nr=50, dtheta=2.0, full_theta=False, t_init=4e4, nphot=5000000, ion_cycles=20,
         spec_cycles=3, bl=(1e37, 1e5), disk=False, macro=False, angles=(10, 45, 80), wmin=1000, wmax=3000,
         nwave=10000, no_split=False):
    '''
    Write root.import.txt and root.pf for a wind between theta_min and
    theta_max (deg from the pole); the parameters are those of the command
    line, described above
    '''

    a = SimpleNamespace(**locals())         # the code below refers to the parameters as a.<name>
    a.rho_index = None if a.rho_index == 'wind' else float(a.rho_index)
    a.bl = None if a.disk else a.bl
    a.simple = not a.macro
    a.split = not a.no_split

    re = np.geomspace(a.r_in, a.r_out, a.nr + 1)
    re = np.append(re, re[-1] * re[-1] / re[-2])            # guard row
    if a.full_theta:
        te = np.linspace(0, 90, int(round(90 / a.dtheta)) + 1)
    else:
        # dtheta cells only across the wind; the grid must still run from 0 to 90 deg
        # (rtheta_where_in_grid clamps theta to the grid), so each empty region is one cell
        nw = max(1, int(round((a.theta_max - a.theta_min) / a.dtheta)))
        te = np.linspace(a.theta_min, a.theta_max, nw + 1)
        if a.theta_min > 0:
            te = np.insert(te, 0, 0.0)
        if a.theta_max < 90:
            te = np.append(te, 90.0)
    nth = len(te) - 1
    te = np.append(te, te[-1] + a.dtheta)                   # guard column


    def vel(r):
        return a.v0 + (a.v_inf - a.v0) * np.clip(1 - a.r_in / r, 0, None) ** a.beta


    def shape(r):
        '''rho / A'''
        if a.rho_index is not None:
            return (a.r_in / r) ** a.rho_index
        return 1.0 / (r ** 2 * vel(r))


    # normalization: tau_es = sigma_T * NE_PER_RHO * A * int shape dr
    rr = np.geomspace(a.r_in, a.r_out, 200001)
    integral = np.trapz(shape(rr), rr)
    A = a.tau / (SIGMA_T * NE_PER_RHO * integral)

    # cell averages of rho in r (theta does not matter inside the wind)
    rho_cell = np.zeros(a.nr)
    for i in range(a.nr):
        x = np.geomspace(re[i], re[i + 1], 201)
        rho_cell[i] = A * np.trapz(shape(x), x) / (re[i + 1] - re[i])
    tau_grid = SIGMA_T * NE_PER_RHO * np.sum(rho_cell * np.diff(re[:a.nr + 1]))

    th_c = 0.5 * (te[:-1] + te[1:])
    in_theta = (th_c[:nth] > a.theta_min) & (th_c[:nth] < a.theta_max)
    if not in_theta.any():
        raise SystemExit('no theta cells between %g and %g' % (a.theta_min, a.theta_max))

    # mass-loss rate (both hemispheres) for reference, at r_in
    mu1, mu2 = np.cos(np.radians(te[np.where(in_theta)[0][-1] + 1])), np.cos(np.radians(te[np.where(in_theta)[0][0]]))
    omega = 2 * 2 * np.pi * (mu2 - mu1)
    mdot = rho_cell[0] * a.r_in ** 2 * vel(a.r_in) * omega * 3.15576e7 / 1.989e33

    with open(a.root + '.import.txt', 'w') as f:
        f.write('# i j inwind r theta v_x v_y v_z rho t_e t_r\n')
        for i in range(len(re)):
            for j in range(len(te)):
                inwind = 0 if (i < a.nr and j < nth and in_theta[j]) else -1
                rho = rho_cell[i] if inwind == 0 else 0.0
                v = vel(re[i])
                th = np.radians(te[j])
                f.write('%3d %3d %3d %.8e %.6f %.6e 0.0 %.6e %.6e %.3e %.3e\n'
                        % (i, j, inwind, re[i], te[j], v * np.sin(th), v * np.cos(th), rho, a.t_init, a.t_init))

    print('%s: wind at %g-%g deg (%d theta cells), r = %.2e-%.2e cm' % (a.root, a.theta_min, a.theta_max,
          in_theta.sum(), a.r_in, a.r_out))
    print('  radial tau_es = %.3f (grid %.3f), rho(r_in) = %.3e g/cm3, n_e(r_in) = %.3e, Mdot = %.2e Msun/yr'
          % (a.tau, tau_grid, rho_cell[0], rho_cell[0] * NE_PER_RHO, mdot))

    angles = a.angles
    nsel = [20000, 0, 1, -1] if a.split else [20000]   # >= MAXSCAT (20000) selects all photons
    obs = ''.join('Spectrum.angle(0=pole)                     %g\n' % x for x in angles for n in nsel)
    phase = ''.join('Spectrum.orbit_phase(0=inferior_conjunction)   0.5\n' for x in angles for n in nsel)
    if a.split:
        scat = '@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)                   yes\n'
        scat += ''.join('@Spectrum.select_scatters                   %d\n' % n for x in angles for n in nsel)
    else:
        scat = '@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)                   no\n'
    if a.bl:
        disk = 'Disk.radiation(yes,no)                          no\n'
        bl = ('Boundary_layer.radiation(yes,no)                   yes\n'
              'Boundary_layer.rad_type_to_make_wind(bb,models,power)                   bb\n'
              'Boundary_layer.luminosity(ergs/s)          %.4e\n'
              'Boundary_layer.temp(K)                     %.4e\n' % (a.bl[0], a.bl[1]))
        bl_spec = 'Boundary_layer.rad_type_in_final_spectrum(bb,models,uniform)                   bb\n'
        disk_spec = ''
    else:
        disk = ('Disk.radiation(yes,no)                          yes\n'
                'Disk.rad_type_to_make_wind(bb,models)                   bb\n')
        bl = 'Boundary_layer.radiation(yes,no)                   no\n'
        bl_spec = ''
        disk_spec = 'Disk.rad_type_in_final_spectrum(bb,models,uniform)                   bb\n'

    if a.simple:
        wrad, ltrans, matom, atomic = 'no', 'thermal_trapping', '', 'data/standard80.dat'
    else:
        wrad, ltrans, atomic = 'yes', 'macro_atoms_thermal_trapping', 'data/h20_hetop_standard80.dat'
        matom = 'Matom_transition_mode(mc_jumps,matrix)             mc_jumps\n'

    pf = f'''System_type(star,cv,bh,agn,previous)                   cv

### Parameters for the Central Object
Central_object.mass(msol)                  0.8
Central_object.radius(cm)                  {R_WD:.1e}
Central_object.radiation(yes,no)                  yes
Central_object.rad_type_to_make_wind(bb,models)                   bb
Central_object.temp                        40000
Binary.mass_sec(msol)                      0.6
Binary.period(hr)                          300

### Parameters for the Disk (if there is one)
Disk.type(none,flat,vertically.extended)                 flat
{disk}Disk.temperature.profile(standard,readin)             standard
Disk.mdot(msol/yr)                         1e-8
Disk.radmax(cm)                            2.4e+10

### Parameters for Boundary Layer or the compact object in an X-ray Binary or AGN
{bl}
### Parameters describing the various winds or coronae in the system
Wind.radiation(yes,no)                          {wrad}
Wind.number_of_components                  1
Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)             imported
Wind.coord_system(spherical,cylindrical,polar,cyl_var)          polar
Wind.model2import                          {a.root}.import.txt
@Diag.adjust_grid(yes,no)                        no

### Parameters associated with photon number, cycles,ionization and radiative transfer options
Photons_per_cycle                          {a.nphot}
Ionization_cycles                          {a.ion_cycles}
Spectrum_cycles                            {a.spec_cycles}
Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow,matrix_est)            matrix_bb
Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping) {ltrans}
{matom}Atomic_data                                {atomic}
@Diag.write_atomicdata(yes,no)                   no
Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)               absorb
Wind_heating.extra_processes(none,adiabatic,nonthermal,both)            adiabatic

### Parameters for Domain 0
Wind.radmax(cm)                            {a.r_out:.3e}
Wind.t.init                                {a.t_init:g}
Wind.filling_factor(1=smooth,<1=clumped)   1

### Parameters defining the spectra seen by observers
Central_object.rad_type_in_final_spectrum(bb,models,uniform)                   bb
{disk_spec}{bl_spec}
### The minimum and maximum wavelengths in the final spectra
Spectrum.nwave                    {a.nwave}
Spectrum.wavemin(Angstroms)                {a.wmin:g}
Spectrum.wavemax(Angstroms)                {a.wmax:g}

### The observers and their location relative to the system
Spectrum.no_observers                      {len(angles) * len(nsel)}
{obs}{phase}Spectrum.live_or_die(live.or.die,extract)              extract
{scat}@Spectrum.select_photons_by_position(yes,no)                   no
@Spectrum.polarization(yes,no)                   yes
Spectrum.type(flambda,fnu,basic)              flambda

### Parameters for Reverberation Modeling (if needed)
Reverb.type(none,photon,wind,matom)                 none

### Other parameters
Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)                    cv
@Diag.use_standard_care_factors(yes,no)                  yes
@Diag.extra(yes,no)                              no
'''
    open(a.root + '.pf', 'w').write(pf)
    print('  wrote %s.pf and %s.import.txt; run with  sirocco -d -p %s' % (a.root, a.root, a.root))


def steer(argv):
    '''
    Process the command line
    '''

    words = []
    x = {}

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] in ('-tau', '-r_in', '-r_out', '-v0', '-v_inf', '-beta', '-dtheta', '-t_init', '-wmin',
                         '-wmax'):
            x[argv[i][1:]] = float(argv[i + 1])
            i += 1
        elif argv[i] in ('-nr', '-nphot', '-ion_cycles', '-spec_cycles', '-nwave'):
            x[argv[i][1:]] = int(argv[i + 1])
            i += 1
        elif argv[i] == '-rho_index':
            x['rho_index'] = argv[i + 1]
            i += 1
        elif argv[i] == '-bl':
            x['bl'] = (float(argv[i + 1]), float(argv[i + 2]))
            i += 2
        elif argv[i] == '-angles':
            angles = []
            while i + 1 < len(argv) and is_number(argv[i + 1]):
                i += 1
                angles.append(float(argv[i]))
            x['angles'] = angles
        elif argv[i] in ('-full_theta', '-disk', '-macro', '-no_split'):
            x[argv[i][1:]] = True
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            words.append(argv[i])
        i += 1

    if len(words) != 3:
        print(__doc__)
        return

    doit(words[0], float(words[1]), float(words[2]), **x)


def is_number(word):
    '''True if word can be read as a number'''
    try:
        float(word)
        return True
    except ValueError:
        return False


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
