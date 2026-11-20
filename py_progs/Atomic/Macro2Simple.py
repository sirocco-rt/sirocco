#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create versions of a MacroAtom imputs that can be used in 
Simple atom calculations.

This must be run in the same directory as MakeMacro.py


Command line usage (if any):

    usage: Macro2Simple.py [-nlev xxx] ion1 ion2 ...

    where:

        -nlev xxx Resets the maxium number of levels to a specific value
        ion1 ion2 .... etc  names of ion in the format we are using

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240127 ksl Coding begun
240522 ksl Mods to make the program work with current versions of MakeMacro

'''

from MakeMacro import get_levels
from astropy.io import ascii
import RedoPhot
import ChiantiPy.core as ch



def write_simple_levels(ion="h_1", nlevels=10):

    try:
        xtab=get_levels(ion,nlevels)
    except:
        print('Could not get macro levels for ion %s %d' % (ion,nlevels))
        print('Was MakeMacroRun for this ion?')
        return ''

    xtab['Dtype']='LevTop'
    xtab['eqn']=-99.0
    
    ztab=xtab['Dtype','Element','Ion','islp','ilv','ion_pot','ex','g','eqn','rad_rate','config']
    outfile='Adata/%s_levels_simple.dat' % ion
    ztab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return outfile



def write_simple_lines(ion="h_1"):
    infile='Adata/%s_lines.dat' % ion
    outfile='Adata/%s_lines_simple.dat' % ion
    try:
        x=ascii.read(infile)
    except:
        print('Error: write_simple_lines: Could not read %s' % infile)
        print('Was MakeMacroRun for this ion?')
        return []
    x['Dtype']='Line'
    x.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return outfile
    



def doit(ion='h_1',nlevels=10):
    lev_file=write_simple_levels(ion,nlevels)
    if len(lev_file)==0:
        print('Failed for ion %s with levels %d' % (ion,nlevels))
        return

    print(lev_file)
    line_file=write_simple_lines(ion)
    print(line_file)
    x = ch.ion(ion, temperature=1e5)
    filename = "Phot/p%02d.%02d.dat" % (x.Z,x.Ion)
    outfile = "Adata/%s_phot_simple.dat" % (ion)  
    phot_tab,xtab=RedoPhot.extrap(filename,1e5)
    if len(phot_tab)==0:
        print('Error: Exiting because of previous problems')
        return    
    RedoPhot.write_phot_tab(outfile,phot_tab,xtab)
    print(outfile)
    return
    

def steer(argv):
    '''
    Set up the various inputs, and call the actual routine
    to convert an ion that exists to a simple ion
    '''

    i=1
    ion=[]
    nlevels=1000
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-nlev':
            i+=1
            nlevels=int(argv[i])
        else:
            ion.append(argv[i])
        i+=1

    if len(ion)==0:
        print('Error: No ins to convert: ',argv)
        return 

    for one_ion in ion:
        if one_ion.isdigit():
            print('Improperly formatted  command line:',argv)
            print(__doc__)
            return

    for one_ion in ion:
        doit(one_ion,nlevels)

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)



