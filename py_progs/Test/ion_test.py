#!/usr/bin/env python
"""
Ionization test: run thin shell models with Sirocco and compare the mean ionization
state of H, He, C, N, O and Fe with Cloudy.

Each model is a single cell, a shell 1 cm thick at 1e11 cm with n_H = 1e7 cm^-3 and T_e fixed at
1e4 K, so the ionization depends only on the SED, U, n_H, T_e and the atomic data, and can be
compared with Cloudy. There are two sets of three models:
    simple atoms (data/standard80.dat):
        ion_test_lowU         1e6 K blackbody, log U = -3.5 (H and He partly ionized)
        ion_test_highU        1e6 K blackbody, log U = -1.25 (C, N and O partly ionized)
        ion_test_pl           Cloudy table power law, log U = -3, with cloudy_test banding, which gives
                              fewer ionization bands than photon generation bands
    H and He as macro atoms (data/h20_hetop_standard80.dat):
        ion_test_matom_lowU   as ion_test_lowU, but log U = -3.25
        ion_test_matom_highU  as ion_test_highU
        ion_test_matom_pl     as ion_test_pl
The reference values are from Cloudy 23.01 with the same SED, shell, n_H, T_e and abundances, with
charge exchange for the simple-atom models and without for the macro-atom ones, matching the
atomic data. The mean ionization state chi, which is 0 for all neutrals
and 1 for fully stripped ions, must be within the tolerance Cloudy for all elements in ELEMENTS.

The script writes the parameter files, runs sirocco -f and windsave2table on each model in the
current folder (making a data link to $SIROCCO/xdata if there is none), and checks the results:

    python ion_test.py                         all the models
    python ion_test.py ion_test_lowU ...       just these models
    python ion_test.py -check_only             check existing windsave2table outputs

It exits with status 1 if any element is outside the tolerance, or a model fails to run.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

ELEMENTS = ["H", "He", "C", "N", "O", "Fe"]

PF_TEMPLATE = """System_type(star,cv,bh,agn,previous)                   bh
Central_object.mass(msol)                  0.8
Central_object.radius(cm)                  {rstar:.6e}
Binary.mass_sec(msol)                      0.4
Binary.period(hr)                          1000
Disk.type(none,flat,vertically.extended,rmin>central.obj.rad)                 none
Central_object.radiation(yes,no)                  yes
{sed}
Central_object.geometry_for_source(sphere,lamp_post,bubble)               sphere
Wind.number_of_components                  1
Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)                shell
Photons_per_cycle                          100000
Ionization_cycles                          2
Spectrum_cycles                            0
Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow,matrix_est)           matrix_pow
Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping)     {line_transfer}
Matom_transition_mode(mc_jumps,matrix) matrix
Wind.radiation(yes,no)                           no
Atomic_data                                {masterfile}
Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)              reflect
Wind_heating.extra_processes(none,adiabatic,nonthermal,both)                 none
Shell.wind_mdot(msol/yr)                   4.717415e-20
Shell.wind.radmin(cm)                      1.000000e+11
Shell.wind.radmax(cm)                      1.000000000010e+11
Shell.wind_v_at_rmin(cm)                   1
Shell.wind.v_at_rmax(cm)                   1.000010
Shell.wind.acceleration_exponent           1
Wind.t.init                                1.000000e+04
Wind.filling_factor(1=smooth,<1=clumped)   1
Reverb.type(none,photon,wind,matom)                 none
Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)          {banding}
Photon_sampling.nbands      10
Photon_sampling.low_energy_limit(eV)    0.1
Photon_sampling.high_energy_limit(eV)  100000.0
"""

# the SED, and the photon banding used with it. The blackbody's U is set by the central object radius
SEDS = {
    "bb": ("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems,mono)       bb\n"
           "Central_object.blackbody_temp(K)           1.000000e+06", "logarithmic"),
    "pl": ("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems,mono)       cloudy\n"
           "Central_object.luminosity(ergs/s)          2.125857e+27\n"
           "Central_object.power_law_index             -0.9\n"
           "Central_object.cloudy.low_energy_break(ev)  0.124\n"
           "Central_object.cloudy.high_energy_break(ev) 50000.0", "cloudy_test"),
}

# (line transfer, masterfile) for simple atoms and for H and He as macro atoms
ATOMS = {
    False: ("thermal_trapping", "data/standard80.dat"),
    True: ("macro_atoms_thermal_trapping", "data/h20_hetop_standard80.dat"),
}

# model parameters for the
# root: (SED, central object radius (cm), macro atoms, Cloudy chi for H, He, C, N, O, Fe)
MODELS = {
    "ion_test_lowU": ("bb", 2.503332e+03, False, [0.6490, 0.6629, 0.2438, 0.1313, 0.0784, 0.0462]),
    "ion_test_highU": ("bb", 3.338246e+04, False, [0.9965, 0.9941, 0.7631, 0.7065, 0.6968, 0.2914]),
    "ion_test_pl": ("pl", 1.000000e+10, False, [0.9896, 0.7844, 0.3284, 0.2782, 0.2538, 0.1217]),
    "ion_test_matom_lowU": ("bb", 3.338246e+03, True, [0.7280, 0.7388, 0.3855, 0.3182, 0.3177, 0.1843]),
    "ion_test_matom_highU": ("bb", 3.338246e+04, True, [0.9952, 0.9941, 0.8383, 0.7542, 0.7322, 0.2933]),
    "ion_test_matom_pl": ("pl", 1.000000e+10, True, [0.9900, 0.7844, 0.3603, 0.3207, 0.3123, 0.1590]),
}


def write_pf(root):
    """
    write the parameter file root.pf using PF_TEMPLATE and the various parameters in MODELS and SEDS
    """
    sed, rstar, macro_atoms, chi = MODELS[root]
    line_transfer, masterfile = ATOMS[macro_atoms]
    with open(root + ".pf", "w") as f:
        f.write(PF_TEMPLATE.format(rstar=rstar, sed=SEDS[sed][0], banding=SEDS[sed][1],
                                   line_transfer=line_transfer, masterfile=masterfile))


def run_model(root, sirocco, windsave2table):
    """
    run sirocco -f and windsave2table on root.pf, with their output in root.log; True on succes
    s"""
    with open(root + ".log", "w") as log:
        for cmd in ([sirocco, "-f", root], [windsave2table, root]):
            if subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode != 0:
                print("{}: {} failed, see {}.log".format(root, os.path.basename(cmd[0]), root))
                return False
    return True


def mean_ion(root, el):
    """
    Calculate the mean ionization level in the (single) wind cell, 
    from the windsave2table file root.el.frac.txt
    """
    tab = np.genfromtxt("{}.{}.frac.txt".format(root, el), names=True)
    fracs = np.array([tab[c] for c in tab.dtype.names if c[0] == "i" and c[1:].isdigit()]).T[tab["inwind"] == 0][0]
    return (np.sum(np.arange(1, len(fracs) + 1) * fracs) - 1) / (len(fracs) - 1)


def check(roots, tol):
    """
    print chi against the reference for each model and element; 
    Returns True if all are within tolerance (tol)
    """

    ok = True
    print("{:22s}".format("chi Sirocco/Cloudy") + "".join("{:>14s}".format(el) for el in ELEMENTS))
    for root in roots:
        if not os.path.exists(root + ".H.frac.txt"):
            print("{}: no windsave2table output".format(root))
            ok = False
            continue
        line = "{:22s}".format(root)
        for el, chi_ref in zip(ELEMENTS, MODELS[root][3]):
            chi = mean_ion(root, el)
            bad = abs(chi - chi_ref) > tol
            ok = ok and not bad
            line += "{:>14s}".format("{:.3f}/{:.3f}{}".format(chi, chi_ref, " X" if bad else ""))
        print(line)
    return ok


def main():

    # see if SIROCCO env var is set 
    if "SIROCCO" in os.environ:
        sirocco_bin_dir = os.path.join(os.environ["SIROCCO"], "bin")
    else:
        sirocco_bin_dir = ""

    # parse the command line
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("roots", nargs="*", default=list(MODELS), help="models to run (default: all)")
    parser.add_argument("-sirocco", default=os.path.join(sirocco_bin_dir, "sirocco"),
                        help="the sirocco executable (default: %(default)s)")
    parser.add_argument("-windsave2table", default=os.path.join(sirocco_bin_dir, "windsave2table"),
                        help="the windsave2table executable (default: %(default)s)")
    parser.add_argument("-tol", type=float, default=0.05, help="tolerance in chi (default: %(default)s)")
    parser.add_argument("-check_only", action="store_true", help="only check existing windsave2table outputs")
    args = parser.parse_args()

    unknown = [root for root in args.roots if root not in MODELS]
    if unknown:
        parser.error("unknown models {}; the models are {}".format(unknown, list(MODELS)))

    if not args.check_only: 
        if not os.path.exists("data"):
            if not sirocco_bin_dir:
                sys.exit("ERROR: no data link here and $SIROCCO is not set")

            os.symlink(os.path.join(os.environ["SIROCCO"], "xdata"), "data")

        for root in args.roots:
            print("running {}".format(root), flush=True)
            write_pf(root)
            run_model(root, args.sirocco, args.windsave2table)

    if check(args.roots, args.tol):
        print("Ionization test passed: every chi within {} of Cloudy".format(args.tol))
    else:
        print("ERROR: ionization test failed: chi outside {} of Cloudy (marked X), or a model did not run".format(args.tol))
        sys.exit(1)


if __name__ == "__main__":
    main()
