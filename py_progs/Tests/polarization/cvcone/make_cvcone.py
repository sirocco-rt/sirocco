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
import argparse
import numpy as np

SIGMA_T = 6.6524587e-25
MPROT = 1.672621923e-24
NE_PER_RHO = 1.19 / MPROT       # fully ionized, Sirocco's standard abundances
R_WD = 7e8

ap = argparse.ArgumentParser()
ap.add_argument('root')
ap.add_argument('theta_min', type=float, help='deg from the pole')
ap.add_argument('theta_max', type=float, help='deg from the pole')
ap.add_argument('-tau', type=float, default=1.0, help='radial electron-scattering depth')
ap.add_argument('-r_in', type=float, default=4 * R_WD)
ap.add_argument('-r_out', type=float, default=1e11)
ap.add_argument('-v0', type=float, default=1e8, help='cm/s at r_in')
ap.add_argument('-v_inf', type=float, default=3e8, help='cm/s')
ap.add_argument('-beta', type=float, default=1.0)
ap.add_argument('-rho_index', default='0',
                help='rho ~ r^-n (default 0, constant density); wind: mass conserving, 1/(r^2 v)')
ap.add_argument('-nr', type=int, default=50)
ap.add_argument('-dtheta', type=float, default=2.0, help='deg')
ap.add_argument('-full_theta', action='store_true',
                help='dtheta cells from 0 to 90 deg (default: dtheta cells only across the wind)')
ap.add_argument('-t_init', type=float, default=4e4)
ap.add_argument('-nphot', type=int, default=5000000)
ap.add_argument('-ion_cycles', type=int, default=20)
ap.add_argument('-spec_cycles', type=int, default=3)
ap.add_argument('-bl', type=float, nargs=2, metavar=('L', 'T'), default=[1e37, 1e5],
                help='central source: boundary layer of luminosity L, temperature T (default 1e37 1e5); the disk does not radiate')
ap.add_argument('-disk', action='store_true', help='a radiating disk and white dwarf instead of the boundary layer')
ap.add_argument('-macro', action='store_true',
                help='macro atoms and wind emission (default: simple atoms, no wind emission, a pure scattering test)')
ap.add_argument('-angles', type=float, nargs='+', default=[10, 45, 80], help='observer inclinations (deg)')
ap.add_argument('-wmin', type=float, default=1000, help='shortest wavelength of the spectra (A)')
ap.add_argument('-wmax', type=float, default=3000, help='longest wavelength of the spectra (A)')
ap.add_argument('-nwave', type=int, default=10000)
ap.add_argument('-no_split', action='store_true',
                help='one observer per angle (default four: all photons, unscattered (sc:0), scattered once (sc:1), scattered (sc:>1))')
args = ap.parse_args()
a = args
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
Binary.period(hr)                          5.57

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
