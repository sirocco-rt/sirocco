#!/usr/bin/env python

"""
Uses Chianti and Topbase to create the base for macro-atom models.

Use Chianti and Topbase to create a set of data files that
can be used as the basis for creating macro-atom models of various ions


Command line usage (if any):

    usage::

        MakeHMacro.py ion_name nlevels [True]

    where the ion name is in Chianti notation, e.g c_4 for C IV, fe_25 for Fe XXV and
    nlevels is the number of energy levels to include in the model

    and the optional True implies this is the top ion to include.  
    
    * Changes as per 27/08/20
    
    usage (in terminal window)::

        MakeMacro.py ion_name nlevels True/False

Description:  

    This routine uses ChiantiPy to access the Chianti database for level line 
    and collision information about an ion.  It attempts to create MacroAtom versions level and line 
    and collision files from the Chianti data  Then the routine retrieves (if necessary) the TopBase 
    photoinization data from Vizier, and constructs photoinization x-sections files that match
    the level files.  

    As currently constructed the routine (attaches a line for the next ion up to the level file) so that
    if one trying to use these data with multiple ions, one needs to make sure to remove all but the
    highest ion+1 line.  As an example suppose you want to construct a model for CIII, CIV, and CV.
    Then one must remove the last line from the level 5 for c_3 and c_4 but not c_5.  

Primary routines:

    doit is the routine run from the command line, and creates all of the data needed for a single 
    ion

Notes:

    There is no guarantee that the file that is produced will be a file that one wants from a physical
    perspective so this routine should be used with great care.  

    In order for this program to work, one must have Chianti and ChiantiPy properly installed on your
    machine.
                                       
History:

191227 ksl
    Coding begun
221227 ksl
    Relooked at routine, verified it seemed to work and cleaned up some of the comments.
    The functionality is unchanged.

"""

import sys
from astropy.io import ascii
import numpy
import os

# Do not call this when we're on ReadTheDocs
if not os.environ.get('READTHEDOCS'):
    import ChiantiPy.core as ch

import numpy as np
from astropy.table import Table, join
from astropy.io import ascii
import RedoPhot


# The levels format looks like
#
# ```
# #         z ion lvl ion_pot   ex_energy   g  rad_rate
# LevMacro  1  1  1 -13.59843    0.000000   2  1.00e+21 () n=1
# LevMacro  1  1  2  -3.39960    10.19883   8  1.00e-09 () n=2
# LevMacro  1  1  3  -1.51093    12.08750  18  1.00e-09 () n=3
# LevMacro  1  1  4  -0.84989    12.74854  32  1.00e-09 () n=4
# LevMacro  1  1  5  -0.54393    13.05450  50  1.00e-09 () n=5
# ...
# LevMacro  1  1 19  -0.03766    13.56077 722  1.00e-09 () n=19
# LevMacro  1  1 20  -0.03399    13.56444 800  1.00e-09 () n=20
# LevMacro  1  2  1   0.00000    13.59843   1  1.00e+21 () cnt
#
# ```
#
# Note, it appears that for our current H macro atom model we have combined some states, e.g all of the n=2 levles into one.  This is why the g is not the same as above

# Before we write out the file we have to add a row for H2.  One could do this cleverly, and make it general,
# but for now we simple add the row that is needed


def get_levels(ion="h_1", nlevels=10):
    """
    Extract the information needed to write out the
    the levels file, and return it as an astropy table.
    
    Note:  In Keara's notes, she derives the multiplicity
    from J, but the multiplicity can be accessed directly
    
    This returns the information needed to write out a level
    file.  Depending on whether one plans to include more
    ions, one needs to add a dummy line for the next ion up
    for the purposes of created a full data set.

    The excitation energy has to be conistent across ions
    
    Sirocco can read the astropy table for levels directly,
    so normally one would write this to a file
    """
    print('get_levels: Getting levels for %s and nlevels %s' % (ion,nlevels))
    try:
        x = ch.ion(ion, temperature=1e5)
    except:
        print('Error: get_levels: failed for %s' % ion)
        raise IOError
    First_Ion_Pot = x.Ip

    # Find the ionization potential of lower ionization
    # states as well

    ex_offset = 0
    i = 1
    while i < x.Ion:
        words = ion.split("_")
        element_string = "%s_%d" % (words[0], i)
        xx = ch.ion(element_string, 1e5)
        ex_offset += xx.Ip
        i += 1

    n = 0
    lvl = []
    rad_rate = []
    ex_energy = []
    g = []
    ion_pot = []
    islp = []
    configuration = []
    if len(x.Elvlc["lvl"]) < nlevels:
        print(
            "There are only %d levels, not %d as requested"
            % (len(x.Elvlc["lvl"]), nlevels)
        )
        nlevels = len(x.Elvlc["lvl"])
    while n < nlevels:
        lv = x.Elvlc["lvl"][n]
        lvl.append(lv)
        if lv == 1:
            rad_rate.append(1e21)
        else:
            rad_rate.append(1e-9)

        #  Calculate the parity.  Note that one cannot calculate the parity from
        #  L in LS coupling.  It is calculated from the the indidual l's of the
        #  electrons in unfilled shells.  It is given by the number of electrons
        #  in the p, f, or h orbitals, so for a term of 2s2p3, which means you 1
        #  electorn in the 2s state and 3 in the 2p state, you count to 3, and
        #  this has odd parity, which is 1 in our notation.  This information is
        #  contained in the term

        term = x.Elvlc["term"][n]
        xterm = list(term)

        j = 0
        num = 0
        while j < len(xterm):
            if xterm[j] == "p" or xterm[j] == "f" or xterm[j] == "h":
                try:
                    delta = int(xterm[j + 1])
                except:
                    delta = 1
                num += delta
            j += 1

        # print('%14s %2d' % (term,num))

        parity = num % 2
        # if x.Elvlc["l"][n] % 2 == 1:
        #     p = 1
        # else:
        #     p = 0

        xislp = x.Elvlc["spin"][n] * 100 + x.Elvlc["l"][n] * 10 + parity
        islp.append(xislp)
        ex = x.Elvlc["ecm"][n] / 8065.5
        ex_energy.append(ex)
        # g.append(x.Elvlc['j'][n]*2+1)
        g.append(x.Elvlc["mult"][n])
        ion_pot.append(-(First_Ion_Pot - ex))
        pretty = x.Elvlc["pretty"][n]
        pretty = pretty.replace(" ", "_")
        configuration.append(pretty)
        n += 1

    # Before setting up the table add the offset to ex
    ex_energy = np.array(ex_energy) + ex_offset

    xtab = Table(
        [lvl, ion_pot, ex_energy, g, rad_rate, configuration, islp],
        names=["lvl", "ion_pot", "ex", "g", "rad_rate", "config", "islp"],
    )
    xtab["Element"] = x.Z
    xtab["Ion"] = x.Ion
    xtab["Dtype"] = "LevMacro"
    xxtab = xtab[
        "Dtype",
        "Element",
        "Ion",
        "lvl",
        "ion_pot",
        "ex",
        "g",
        "rad_rate",
        "config",
        "islp",
    ]
    xxtab["ion_pot"].format = ".6f"
    xxtab["ex"].format = ".6f"

    # Now try to set up the TopBass Ilv for this file
    # The approach here is to use the absolute value of the
    # difference in the energy, but this may not be ideal.

    qslp = np.unique(xxtab["islp"])

    xxtab["ilv"] = -99
    #    for q in qslp:
    #        i = 0
    #        n = 0
    #        e_old = 0
    #        while i < len(xxtab):
    #            if xxtab["islp"][i] == q:
    #                if n == 0:
    #                    e_old = xxtab["ex"][i]
    #                    n += 1
    #                elif abs(xxtab["ex"][i] - e_old) < 0.1:
    #                    pass
    #                else:
    #                    e_old = xxtab["ex"][i]
    #                    n += 1
    #                xxtab["ilv"][i] = n
    #            i += 1

    # This is an alternative, which simply assumes that
    # All of the sublevels are displayed at once

    #    xxtab["ilv"] = -99
    #    for q in qslp:
    #        i = 0
    #        n = 0
    #        while i < len(xxtab):
    #            if xxtab["islp"][i] == q:
    #                if n==0:
    #                    i_old=i
    #                    n+=1
    #                elif i_old+1==i:
    #                    i_old=i
    #                else:
    #                    i_old=i
    #                    n+=1
    #                xxtab["ilv"][i] = n
    #            i+=1

    # Use the term to idenfify all of the sublevels that are associated with
    # one topbase level. This approach assumes that all fo the sublevels assoicated
    # with a single level are contiguous in the Elvlc dictionary.
    for q in qslp:
        i = 0
        n = 0
        term_old = ""
        while i < len(xxtab):
            if xxtab["islp"][i] == q:
                if n == 0:
                    term_old = x.Elvlc["term"][i]
                    n += 1
                elif x.Elvlc["term"][i] == term_old:
                    pass
                else:
                    term_old = x.Elvlc["term"][i]
                    n += 1
                xxtab["ilv"][i] = n
            i += 1

    print('The number of macro levels is %d' % (len(xxtab)))
    return xxtab



MELEC=9.10956e-28
C=2.997925e10
ECHARGE=4.8035e-10 
def get_f(one_line):
    '''
    Calculate an absorption oscillator strength
    from an Einstein A
    '''
    # print(one_line)
    
    factor=MELEC*C/(8.* np.pi*ECHARGE**2)
    xwave=(one_line['Wave']*1e-8)**2
    a=one_line['a']
    f=factor*xwave*a
    f*=one_line['gu']/one_line['gl']
    # f=factor
    return f

def calculate_oscillator_strength(A_ij, wavelength_AA):
    # Constants in CGS units
    m_e = 9.109e-28  # Electron mass in grams
    c = 2.997e10  # Speed of light in cm/s
    h = 6.626e-27  # Planck's constant in erg s

    wavelength_cm=wavelength_AA*1e8

    # Calculate the absorption oscillator strength (f)
    f = (4 * np.pi**2 * m_e * c**3 / (3 * h)) * (A_ij / wavelength_cm**2)
    
    return f

# # Lines

# The format of the lines file is a follows:
#
# ```
# # z = element, ion= ionstage, f = osc. str., gl(gu) = stat. we. lower(upper) level
# # el(eu) = energy lower(upper) level (eV), ll(lu) = lvl index lower(upper) level
# #        z ion       lambda      f         gl  gu    el          eu        ll   lu
# LinMacro    1   1    1215.33907  0.41620     2     8    0.00000  10.19883     1     2
# LinMacro    1   1    1025.44253  0.07910     2    18    0.00000  12.08750     1     3
# LinMacro    1   1     972.27104  0.02899     2    32    0.00000  12.74854     1     4
# LinMacro    1   1     949.48382  0.01394     2    50    0.00000  13.05450     1     5
# LinMacro    1   1     937.54762  0.00780     2    72    0.00000  13.22070     1     6
# LinMacro    1   1     930.49400  0.00471     2    98    0.00000  13.32092     1     7
# LinMacro    1   1     925.97292  0.00318     2   128    0.00000  13.38596     1     8
# LinMacro    1   1     922.89865  0.00222     2   162    0.00000  13.43055     1     9
# LinMacro    1   1     920.71177  0.00161     2   200    0.00000  13.46245     1    10
# LinMacro    1   1     919.10059  0.00120     2   242    0.00000  13.48605     1    11
# LinMacro    1   1     917.87888  0.00092     2   288    0.00000  13.50400     1    12
# ...
# LinMacro    1   1    6562.83858  0.64080     8    18   10.19883  12.08750     2     3
# LinMacro    1   1    4861.35082  0.11930     8    32   10.19883  12.74854     2     4
#
#
# ```


def get_lines(ion="h_4", nlevels=10):
    """
    Given an astropy table that contains all of the levels, 
    get the associated line information.
    
    Notes:
    This routine calls get_lev to decide which lines to retrieve
    so that it only produces a table with lines for which levels
    are known

    There are some lines in the Chianti database that have 0 wavelength.
    These are all strictly forbidden lines, but Sirocco still needs a
    wavelength to calculate the collision x-section in the absence 
    of Burgess-style collision data, so we fill in the wavelength
    from the upper and lower level energies in this case, and issue
    a warning.
    """

    lev = get_levels(ion, nlevels)
    # Next bit is duplicated if I decide to put into a single routine
    x = ch.ion(ion, temperature=1e6)
    # End of duplication

    wavelength = x.Wgfa["wvl"]
    lower = x.Wgfa["lvl1"]
    upper = x.Wgfa["lvl2"]
    gf = x.Wgfa["gf"]
    a = x.Wgfa["avalue"]
    xtab = Table(
        [wavelength, lower, upper, gf, a], names=["Wave", "ll", "ul", "gf", "a"]
    )
    # xtab.info()
    select = []
    gl = []
    gu = []
    el = []
    eu = []
    i = 0
    for one in xtab:
        got_upper = False
        got_lower = False
        for one_level in lev:
            if one_level["lvl"] == one["ll"]:
                got_lower = True
                xgl = one_level["g"]
                xel = one_level["ion_pot"]
            if one_level["lvl"] == one["ul"]:
                got_upper = True
                xgu = one_level["g"]
                xeu = one_level["ion_pot"]
            if got_lower and got_upper:
                select.append(i)
                gl.append(xgl)
                gu.append(xgu)
                el.append(xel)
                eu.append(xeu)
                break
        # print(i,len(xtab))
        i += 1
    # print(select)
    xxtab = xtab[select]
    xxtab["gl"] = gl
    xxtab["gu"] = gu
    xxtab["f"] = xxtab["gf"] / xxtab["gl"]
    xxtab["el"] = x.Ip + el
    xxtab["eu"] = x.Ip + eu
    xxtab["el"].format = "10.6f"
    xxtab["eu"].format = "10.6f"
    xxtab["f"].format = "9.6f"

    xxtab["Dtype"] = "LinMacro"
    xxtab["z"] = lev["Element"][0]
    xxtab["ion"] = lev["Ion"][0]

    # Check for lines that have a Wavelength of 0 and if tis happens fix the
    # wavelegnth using el and eu, but issue a warning when you do this
    # Also assure that the f values are not zero.  Note that this
    # value may be higher than we want, as the f values for forbidden lines
    # may be even lower than this.  Howevever, if one decides to change this
    # One needs to be that output formats for various routines capture the full range, as the
    # f value is compared to the value for collisions.

    f_min=0.000001
    n_fmissing=0

    for one in xxtab:
        if one["Wave"] == 0:
            one["Wave"] = 12396.1914 / (one["eu"] - one["el"])
            # print(
            #     "Line with ll %d and ul %d is missing wavelength.  Correcting to %.3f using el and eu"
            #     % (one["ll"], one["ul"], one["Wave"])
            # )
        if one["f"]<f_min:
            fff=get_f(one)
            # print('Toast %f %g' % (one['a'],fff))
            one['f']=f_min
            n_fmissing+=1
            # print("Line with ll %d and ul %d is missing f.  Changing to %.6f so line has way out" 
            #       % (one["ll"], one["ul"],one["f"]))

    print('There were %d lies with missing f' % n_fmissing)


    # Convert theoretical wavelengths (which are negative) to positive values
    i=0
    nnn=0
    while i<len(xxtab):
        if xxtab['Wave'][i]<0:
            xxtab['Wave'][i]=-xxtab['Wave'][i]
            nnn+=1
        i+=1

    print('There were %d theoretical wavelengths, converted to positive values' % nnn)
    xxtab.write('xline.txt',format='ascii.fixed_width_two_line',overwrite=True)

    xxxtab = xxtab["Dtype", "z", "ion", "Wave", "f", "gl", "gu", "el", "eu", "ll", "ul"]
    return xxxtab


# # Photoionization data

# We have to first get the data from TopBase and then we have to format it propoerly
#
# There is a search form on the web page, but you can also retrieve the files directly with wget.  
# The formats are slightly different, and so converting them to Sirocco inputs is also going to be a
# bit different, but this allows one to retrieve things without going through a web page.



def get_phot(ion="c_4"):
    """
    Obtain the photoionization x-sectons from TopBase. This routine uses 
    normal astronomical terminology, that is CIV would be element 6 and 
    ion 4
    
    First check if we already have the data, and if not retrieve it
    
    Note that the TopBase file names are based on the element number and the 
    number of electrons that the ion of interest has, but we convert to
    astronomical notation in this routine
    """

    os.makedirs('./Phot',exist_ok=True)
    x = ch.ion(ion, temperature=1e5)
    nelec = x.Z - x.Ion + 1
    fileroot = "p%02d.%02d" % (x.Z, nelec)
    outroot = "./Phot/p%02d.%02d" % (x.Z, x.Ion)
    print("Looking for ", fileroot)
    if os.path.isfile("%s.txt" % outroot):
        print("TopBase Phot file for element %d and ion %d exists" % (x.Z, x.Ion))
        return
    os.system("wget cdsweb.u-strasbg.fr/topbase/p/%s.gz" % (fileroot))
    os.system("gunzip %s.gz" % (fileroot))
    os.system("mv %s %s.txt" % (fileroot, outroot))
    print(
        "TopBase Phot file for element %s and ion %d has now been retrieved"
        % (x.Z, x.Ion)
    )
    return


# The retrieved files begin like
#
# ```
#     6    4    P
#     1    0    0    1
#     393    457
#   3.511710E+00    0.0200
#   3.421708E+00 1.603E+00
#   3.429049E+00 1.468E+00
# ```
#
# and then continue with x-sections for a while and then you see
#
# ```
#   1.000000E+02 1.763E-02
#     1    0    0    2
#     777    820
#   1.805500E+00    0.0200
#   1.715504E+00 1.074E+00
# ```
#
# so the first line is unique, but then you need to find the remainder
#
# The second line ends up being  ISLP and ILV and the second quantity
# in the second line is the number of points I believe
#
# The line after this gives the number of x-sections (as the second number).
# I do not know what the third number is
#
# The line after this which has two entries, gives the excitation energy,
# not a x-section,
#
# There seems to be an empty transition at the end of the file, containg 0 0 0 0
#
# ## There is more work to be done on the routine below, mainly to get
# things into the right units, but also to remove non-physical x-sections.


def make_phot(ion="c_4",macro=True):
    """
    Read a retrieved TopBase file and make a Sirocco-Photon file
    
    Note that althought the TopBase data on Vizier does not use
    astronomical notation, we changed the naming convention
    when we retrieved the files
    """

    x = ch.ion(ion, temperature=1e5)
    z = x.Z
    xion = x.Ion

    fileroot = "p%02d.%02d" % (z, xion)
    try:
        x = open("./Phot/%s.txt" % fileroot)
        lines = x.readlines()
    except:
        print("Error: make_phot %s.txt not found  %s" % fileroot)
        return
    records = []
    for line in lines:
        words = line.split()
        records.append(words)
    f = open("./Phot/%s.dat" % fileroot, "w")
    # z=int(records[0][0])
    # xion=int(records[0][1])
    i = 1
    num = 0
    xsection = 1
    while i < len(records):
        one_record = records[i]
        # print('test5',one_record)
        if len(one_record) == 4:
            # This is a new x-section
            # if i>1:
            #     print('Number of x-sections %d' % num)
            h = records[i]
            # print(h)
            islp = h[0] + h[1] + h[2]
            if islp == "000":
                break
            ilv = h[3]
            i += 1
            # print(records[i])
            h = records[i]
            npts = int(h[1])
            i += 1
            h = records[i]
            ex = eval(h[0]) * 13.605693009
            if ex == 0:
                xh = records[i + 1]
                ex = eval(xh[0]) * 13.605693009
                print(
                    "Warning: ex=0  for  %d %d %s %s %d - changing to %10.6f"
                    % (z, xion, islp, ilv, npts, ex)
                )
            new_xsection=True
            # The string differs depending on whether Macro or NOT.  
            # string = "PhotTopS %d %d %s %s %10.6f %d" % (z, xion, islp, ilv, ex, npts)
            # string='PhotTop %d %d  %d %d %d %d' % (z,xion,islp,ilv,ex,npts)
            # print(string)
            # f.write("%s\n" % string)
            num = 0
            xsection += 1
        else:
            energy = eval(records[i][0]) * 13.605693009
            xsection = eval(records[i][1]) * 1e-18
            if energy >= ex and new_xsection==True:
                if macro==True:
                    string = "PhotTopS %d %d %s %s %10.6f %d" % (z, xion, islp, ilv, ex, npts)
                else:
                    string='PhotTop %d %d  %d %d %d %d' % (z,xion,islp,ilv,ex,npts)
                f.write("%s\n" % string)
                new_xsection=False
            else:
                npts-=1

            if energy >= ex: #added so that only x-sections with energies greater than threshold appear (RG)
                string = "Phot %10.6f %10.6e" % (energy, xsection)
                f.write("%s\n" % string)
            num += 1
        i += 1
        if i == len(records):
            print("Number of x-sections %d" % num)
    f.close()
    return


# **Now we have to generate the Photometry file**
#
# The retrieved files begin like
#
# ```
#     6    4    P
#     1    0    0    1
#     393    457
#   3.511710E+00    0.0200
#   3.421708E+00 1.603E+00
#   3.429049E+00 1.468E+00
# ```
#
# and then continue with x-sections for a while and then you see
#
# ```
#   1.000000E+02 1.763E-02
#     1    0    0    2
#     777    820
#   1.805500E+00    0.0200
#   1.715504E+00 1.074E+00
# ```
#
# so the first line is unique, but then you need to find the remainder
#
# The second line ends up being ISLP and ILV and the second quantity in the second line is the number of points
#
# The line after this gives the number of x-sections (as the second number). I do not know what the third number is
#
# The line after this which has two entries, gives the excitation energy, not a x-section,
#
# There seems to be an empty transition at the end of the file, containg 0 0 0 0
#
#
# We need to convert this into a file with the following format
#
# ```
# PhotMacS       1       1       1       1       13.598430     100
# PhotMac       13.598430   6.3039999e-18
# PhotMac       13.942675   5.8969998e-18
# PhotMac       14.295634   5.5159998e-18
# ```
# or looking at just the headers for H
#
# ```
# Pescado:macroatom : grep PhotMacS data/atomic_macro/h20_phot.dat
# PhotMacS       1       1       1       1       13.598430     100
# PhotMacS       1       1       2       1       3.3996000     100
# PhotMacS       1       1       3       1       1.5109299     100
# PhotMacS       1       1       4       1      0.84988999     100
# PhotMacS       1       1       5       1      0.54392999     100
# PhotMacS       1       1       6       1      0.37773001      99
# PhotMacS       1       1       7       1      0.27750999      99
# PhotMacS       1       1       8       1      0.21247000     100
# PhotMacS       1       1       9       1      0.16788000      99
# ```
# or
# ```
# PhotMacS        6       4       1       1       64.432072       82
# PhotMac 64.432113       6.599E-19
# PhotMac 67.775753       6.068E-19
# PhotMac 71.222782       5.580E-19
# PhotMac 74.777800       5.132E-19
# PhotMac 78.445351       4.720E-19
# ```
#
# or looking just at the headers
#
# ```
# PhotMacS        6       4       1       1       64.432072       82
# PhotMacS        6       4       2       1       56.434374       68
# PhotMacS        6       4       3       1       56.434374       68
# PhotMacS        6       4       4       1       26.927027       71
# PhotMacS        6       4       5       1       24.791477       71
# PhotMacS        6       4       6       1       24.791477       71
# PhotMacS        6       4       7       1       24.210786       77
# PhotMacS        6       4       8       1       24.210786       77
# PhotMacS        6       4       9       1       14.725714       76
# PhotMacS        6       4       10      1       13.860936       77
# ```
#
# The headers are:
# PhotMacS element # ion # level # level up (1) energy threshold [eV] # of data lines
#
# PhotMac energy[eV] cross sections [cm2]
#
# The hard part here is getting the levels to be right.  If I understand Keara's notes, there are fewer sets of
# photoinzation x-sections in top base than there are in Chianti, and so often one needs to have the same set
# of x-section data for multiple Chianti levels, and when you do this you have to right the data out multiple times.
# because Sirocco expects a photionzation x-section for each level
#
# Alternatively, I suspect you can combine some of the levels into a single simpler level, if you know how to modify
# various factors, e.g g, and you make some kind of assumption about what you need to do interms of adding f together.
# This simplfies the model that you have


def write_phot(ion="c_4",outdir='./Adata'):

    x = ch.ion(ion, temperature=1e5)
    xion = x.Ion
    z = x.Z

    fileroot = "./Phot/p%02d.%02d" % (x.Z, x.Ion)

    # First strip of the headers of the phot_file
    try:
        phot = open(fileroot + ".dat")
        lines = xx = phot.readlines()
    except:
        print("Error: write_phot: Could not open %s.dat" % fileroot)

    headers = []
    i = 0
    for one in xx:
        if one.count("PhotTopS"):
            word = one.split()
            word[0] = i
            word[1] = int(word[1])
            word[2] = int(word[2])
            word[3] = int(word[3])
            word[4] = int(word[4])
            word[5] = eval(word[5])
            word[6] = int(word[6])
            headers.append(word)
        i += 1
    phot.close()
    # Now we have a set of records we can make into a table if we are clever
    headers = np.array(headers)
    xxhead = Table(
        headers,
        names=["Row", "Element", "Ion", "islp", "ilv", "e_thresh", "np"],
        dtype=["i", "i", "i", "i", "i", "f", "i"],
    )

    # Write what we have so far for diagnostic purposes. This file is not used.
    # xxhead.write("head.txt", format="ascii.fixed_width_two_line", overwrite=True)

    # Now find, if we can the parts of the photon file that match what we need

    # lev_file = ion + "_levels.dat"
    lev_file = '%s/%s_levels.dat' %(outdir,ion)
    try:
        lev = ascii.read(lev_file)
    except:
        print("Error: Could not open level file %s" % lev_file)

    # lev.info()

    # Join the levels to the photoionization data                                 

    try:
        foo = join(
            lev, xxhead, keys=["Element", "Ion", "islp", "ilv"], join_type="left"
        )
    except:
        print("Error: join failed")
        xxhead.info()
        lev.info()
        # xxhead.write("head.txt", format="ascii.fixed_width_two_line", overwrite=True)
        lev.write("levels.txt", format="ascii.fixed_width_two_line", overwrite=True)



    # Reorder the table so that the photoionization x-sections can be written in 
    # level order

    foo.sort(["Ion", "lvl"])

    foo["delta_e"] = foo["ion_pot"] + foo["e_thresh"]
    foo["delta_e"].format = "8.2f"

    # This write is diagnostic and lev2phot is not used anywhere.  The Row columin 
    # refers to the row in the input file, and has nothing to do with 
    # The location of a particular x-section in the ouput PhotFile, which
    # is reordered to be in energy level order

    foo.write(
        ion + "_lev2phot.txt", format="ascii.fixed_width_two_line", overwrite=True
    )

    # Now write out the PhotFile. We use the row in the original file
    # to understand which lines to write out.  

    # output_file = ion + "_phot.dat"
    output_file = '%s/%s_phot.dat' % (outdir,ion)
    f = open(output_file, "w")
    for one in foo:
        # print(one)
        if one["e_thresh"] > 0:
            xstring = "PhotMacS  %2d %2d %2d %2d  %10.6f %3d" % (
                z,
                xion,
                one["lvl"],
                1,
                one["e_thresh"],
                one["np"],
            )
            f.write("%s\n" % xstring)
            i = int(one["Row"])
            nphot = int(one["np"])
            j = i + 1
            jstop = j + nphot
            # print('%d %d %s'% (j,jstop,lines[j]))
            while j < jstop:
                try:
                    xstring = "PhotMac       13.598430   6.3039999e-18"
                    xstring = lines[j]
                    xstring = xstring.replace("Phot", "PhotMac")
                    f.write(xstring)
                except IndexError:
                    print('Failed on j %d for jstop %d and nphot %d' % (j,jstop,nphot))
                    print(one)
                    break

                j += 1
    f.close()


# # Get the collision data from Chianti

# Generating the collision data files is not entirely straightforward, because we have to match the collision data to the levels, and the source of the level information may not have been Chianti.  The collision data contained in Chianti (and in fact the only type of collision data that Sirocco understnds) uses the so called Burgess approximation, basically a spline fit to the x-sections which is valid up to a certain maximum temperature.
#
# A valid file begins like
#
# ````
# CSTREN Line 1 1 1215.673  0.139 2 2 0.000000 0.200121 0  1  1 3 7.500e-01 2.772e-01 1.478e+00  5 1 1.700e+00
# SCT   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
# SCUPS    1.132e-01   2.708e-01   5.017e-01   8.519e-01   1.478e+00
# CSTREN Line 1 1 1215.668 0.277 2 4 0.000000 10.200166 0 2  1 4 7.500e-01 5.552e-01  2.961e+00 5 1 1.700e+00
# SCT   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
# SCUPS    2.265e-01   5.424e-01   1.005e+00   1.706e+00   2.961e+00
# CSTREN Line 1 1 1025.722  0.026 2 2  0.000000 12.089051 0 3 1 6 8.890e-01 5.268e-02  2.370e-01 5 1 1.600e+00
# SCT   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
# ````
#
# The first 10 values come directly from the line file, and the rest have to be derived from the collision data.
# These values are
#
# ***Note***
# In exploring the collisions in Chianti, I discovered that for H at least here are collisions strengths for between levels that do not have lines assoicated with them.  I am not sure why that is, but it must affect macro atoms, and it did it might be needed.  A simple join between the lines file and the collisions files reveals those that have both, but a left join easily shows you the ones that are just collision data


def get_collisions(ion="h_1", nlev=20,outdir='./Adata'):
    """
    Given a set of levels, get the collision information associatred with these levels.  Note
    that as far as I can determine the scups file does not contain collisionl information 
    for collisional ionization
    
    Currently in Sirocco collisions are tightly tied to the lines that are going to be used
    in a dataset.  The information about line information is used as a mechanism to identify
    the collision x-section.  (This probably dates back to simple atoms, where lines are used
    for radiative transfer and levels are used in calculating densities of upper level states,
    but we don't really have a concept of an integrated atom)
    """
    # Get the collision data

    x = ch.ion(ion, temperature=1e5)
    lower = x.Scups["lvl1"]
    upper = x.Scups["lvl2"]
    de = x.Scups["de"]
    gf = x.Scups["gf"]
    ntemp = x.Scups["ntemp"]
    btemp = x.Scups["btemp"]
    bscups = x.Scups["bscups"]
    lim = x.Scups["lim"]
    ttype = x.Scups["ttype"]
    cups = x.Scups["cups"]

    xtab = Table(
        [lower, upper, de, lim, ntemp, btemp, bscups, ttype, cups],
        names=["ll", "ul", "de", "lim", "ntemp", "btemp", "bscups", "ttype", "cups"],
    )

    print(xtab)
    xtab=xtab[xtab['ttype']<5]

    npossible = 0
    for one in xtab:
        if one["ul"] <= nlev:
            npossible += 1

    # xtab.write("T_cups.txt", format="ascii.fixed_width_two_line", overwrite=True)
    # Get the lines

    linetab = get_lines(ion=ion, nlevels=nlev)
    # linetab.info()
    linetab.write("T_lines.txt", format="ascii.fixed_width_two_line", overwrite=True)

    xxtab = join(xtab, linetab)
    # xxtab.write("T_all.txt", format="ascii.fixed_width_two_line", overwrite=True)
    xxtab["gf"] = xxtab["gl"] * xxtab["f"]

    print(
        "There were %d (of %d) collision x-sections with nlev <= %d, and %d that matched %d lines"
        % (npossible, len(xtab), nlev, len(xxtab), len(linetab))
    )

    # xtab=Table([lower,upper,gf,ntemp,btemp,bscups],names=['ll','de','ul','gf','ntemp','btemp','bscups'])
    upsfile='%s/%s__upsilon.dat' % (outdir,ion)
    # xxtab.info()
    xout = open(upsfile, "w")
    for one in xxtab:
        # print(one)
        xstring = "CSTREN Line %3d %3d %10.6f %9.6f %2d %2d  %10.6f %10.6f " % (
            one["z"],
            one["ion"],
            one["Wave"],
            one["f"],
            one["gl"],
            one["gu"],
            one["el"],
            one["eu"],
        )
        xstring = xstring + "%3d %3d %3d %3d %10.3e %10.3e %10.3e %3d %3d %10.3e" % (
            one["ll"],
            one["ul"],
            one["ll"],
            one["ul"],
            one["de"],
            one["gf"],
            one["lim"],
            one["ntemp"],
            one["ttype"],
            one["cups"],
        )

        # print(xstring)
        xout.write("%s\n" % xstring)

        xstring = "SCT    "
        for one_temp in one["btemp"]:
            xstring = xstring + (" %10.3e" % one_temp)
        xout.write("%s\n" % xstring)

        xstring = "SCUPS  "
        for one_bscup in one["bscups"]:
            xstring = xstring + (" %10.3e" % one_bscup)
        xout.write("%s\n" % xstring)

    xout.close()
    return xxtab


def print_elvlc(ion="c_4"):
    """
    Print out values of information in the Elvlc file that might be used 
    """
    x = ch.ion(ion, temperature=1e5)
    imax = len(x.Elvlc["lvl"])
    print(imax)
    i = 0
    while i < imax:
        if x.Elvlc["l"][i] % 2 == 1:
            p = 1
        else:
            p = 0
        islp = x.Elvlc["spin"][i] * 100 + x.Elvlc["l"][i] * 10 + p
        e = x.Ip - x.Elvlc["erydth"][i] * 13.605693009
        pretty = x.Elvlc["pretty"][i]
        pretty = pretty.replace(" ", "_")
        print(
            "%2d %10s %d %s %d %5.1f %5.1f %15s  %3d %10.5f"
            % (
                x.Elvlc["lvl"][i],
                x.Elvlc["term"][i],
                x.Elvlc["spin"][i],
                x.Elvlc["spd"][i],
                x.Elvlc["l"][i],
                x.Elvlc["j"][i],
                x.Elvlc["mult"][i],
                pretty,
                islp,
                e,
            )
        )
        i += 1


def doit(atom="h_1", nlev=10, next_ion = False,outdir='./Adata'): 
    """
    Create all of the necessary files for 
    a given atom

    Description:

    Notes:

    History:
    
    2020: RG - added feature to add in first level of next ion (only use for adding in fully ionised state to single electron state) (RG)


    """


    os.makedirs(outdir,exist_ok=True)

    nlev = int(nlev)

    xion = ch.ion(atom, temperature=1e5)
    xlevels = get_levels(atom, nlev)

    # If we could not find the number of levels
    # requested, reduce the number so that the
    # rest of the routines will not generate more comments on this
    if len(xlevels) < nlev:
        nlev = len(xlevels)

    # Find the ionization potential of lower ionization
    # states as well

    ex_offset = 0
    i = 1
    while i < xion.Ion:
        words = atom.split("_")
        element_string = "%s_%d" % (words[0], i)
        xx = ch.ion(element_string, 1e5)
        ex_offset += xx.Ip
        i += 1
    
    if next_ion == True: #if True, we add in the fully ionised state level (RG)  
        xlevels.add_row(
            [
                "LevMacro",
                xion.Z,
                xion.Ion + 1,
                1,
                0.00000,
                xion.Ip + ex_offset,
                1,
                1.00e21,
                "Next",
                0,
                0,
            ]
        )
    xlevels.write('%s/%s_levels.dat' % (outdir,atom), format="ascii.fixed_width_two_line", overwrite=True
    )
    xlines = get_lines(atom, nlev)
    xlines.write('%s/%s_lines.dat' % (outdir,atom) , format="ascii.fixed_width_two_line", overwrite=True
    )

    get_phot(atom)
    make_phot(atom)
    write_phot(atom,outdir)
    RedoPhot.redo_one('Adata/%s_phot.dat' % atom, 'Adata/%s' % atom)

    xcol = get_collisions(atom, nlev)
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys

    if len(sys.argv) == 3:
        doit(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        if sys.argv[3]=='True':
            sys.argv[3]=True
        doit(sys.argv[1], sys.argv[2],sys.argv[3])
    else:
        print(__doc__)
