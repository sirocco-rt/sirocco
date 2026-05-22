#!/usr/bin/env python
# coding: utf-8

'''
Summarize an atomic data file in order to see how large arrays defined in atomic.h must be


Command line usage (if any):

    usage::

        CheckAtomic.py masterfile


Description:  

    The routine reads a master file and
    and then all of the data that would
    be read into Sirocco using this master
    file.  It then parses the concatenated
    list of data and attempts to summarize
    how many inputs are read in.

Primary routines:

    doit

Notes:

    At present it does not work out
    the number of photioniztion x-sections
    or the number of collisional x-sections

    The routine must b run from a directory
    in which the variious files in the 
    master file are properly linked.

                                       
History:

    240511 ksl Coding begun

'''

# # Read and atomic data file and assess the number of macro levels etc



from astropy.table import Table, vstack
import numpy as np
import os


def xread(filename='data/fe_23to27.dat',split=False):
    try:
        x=open(filename,'r')
    except:
        print('Error: Could not locate %s' % filename)
        return None
    lines=x.readlines()
    records=[]
    for one in lines:
        z=one.strip()
        if len(z)>0 and z[0]!='#':
            if split:
                records.append(z.split())
            else:
                records.append(z)
    return records


def get_element_table(lines):
    '''
    Get the elements from the concatenated input file
    '''
    name=[]
    z=[]
    for one in lines:
        # print(one)
        if one[0]=='Element':
            name.append(one[2])
            z.append(int(one[1]))
    print(name)
    # print(z)
    element_tab=Table([name,z],names=['element','z'])
    return element_tab
                                
            


def get_ion_tab(lines,elements):
    '''
    Get the ions from the concatenated
    nput file, given the elements that
    have been already determined.
    '''
    z=[]
    istate=[]
    macro=[]
    name=[]
    max_lev=[]
    max_macro=[]
    for one in lines:
        if one[0].count('Ion'):
            name.append(one[1])
            z.append(int(one[2]))
            istate.append(int(one[3]))
            max_lev.append(int(one[6]))
            max_macro.append(int(one[7]))
            if one[0].count('IonM'):
                macro.append('yes')
            else:
                macro.append('no')
    ion_tab=Table([name, z,istate,macro,max_lev,max_macro],
                  names=['element','z','ion','macro','max_lev','max_macro'])
    
    records=[]
    for one_z in elements['z']:
        xx=ion_tab[ion_tab['z']==one_z]
        if len(xx)>0:
            records.append(xx)
    final_ion_tab=vstack(records)
    return final_ion_tab



def get_level_tab(lines,ions):
    '''
    Get the leverls from the concatenated
    inputs given  the ions that are 
    active
    '''
    z=[]
    istate=[]
    lvl=[]
    macro=[]
    for one in lines:
        # print(one)
        if one[0].count('Lev'):
            z.append(int(one[1]))
            istate.append(int(one[2]))
            lvl.append(int(one[3]))
            if one[0].count('LevMacro'):
                macro.append('yes')
            else:
                macro.append('no')
    # print(z)
    lev_tab=Table([z,istate,lvl,macro],names=['z','ion','level','macro'])
    # print('Hello :',len(lev_tab))

    records=[]
    xions=np.unique(ions['z'])
    for one_z in xions:
        xx=lev_tab[lev_tab['z']==one_z]
        if len(xx)>0:
            records.append(xx)
    final_lev_tab=vstack(records)
    # print('Goodbye',len(final_lev_tab))
    return final_lev_tab



def get_line_tab(lines,lev):
    '''
    Get the lines from the concatenated 
    inputs given a table contaiing the levels
    '''
    z=[]
    istate=[]
    ll=[]
    ul=[]
    macro=[]
    for one in lines:
        # print(one)
        if one[0].count('Lin'):
            z.append(int(one[1]))
            istate.append(int(one[2]))
            ll.append(int(one[9]))
            ul.append(int(one[10]))
            if one[0].count('LinMacro')>0:
                macro.append('yes')
            else:
                macro.append('no')
    # print(z)
    line_tab=Table([z,istate,ll,ul,macro],names=['z','ion','ll','ul','macro'])
    # print(np.unique(line_tab['macro'],return_counts=True))

    records=[]
    zz=np.unique(lev['z'])
    for one_z in zz:
        xx=line_tab[line_tab['z']==one_z]
        # print(len(xx))
        if len(xx)>0:
            records.append(xx)
    final_line_tab=vstack(records)
    # print(np.unique(final_line_tab['macro'],return_counts=True))
    return final_line_tab





def analyze_lines(line_tab):
    '''
    Analyze the lines
    '''
    print('There are %5d lines' % (len(line_tab)))
    macro_lines=line_tab[line_tab['macro']=='yes']
    print('There are %5d macro lines' % (len(macro_lines)))
    zz=np.unique(line_tab['z'])
    # print(np.unique(line_tab['macro'],return_counts=True))
    
    for one_z in zz:
        xtab=line_tab[line_tab['z']==one_z]
        xmacro=macro_lines[macro_lines['z']==one_z]
        zion=np.unique(line_tab['ion'])
        print('\nElement z nlines nmacro: %3d     %10d %10d' % (one_z,len(xtab),len(xmacro)))
        for one_ion in zion:
            xxtab=xtab[xtab['ion']==one_ion]
            xxmacro=xmacro[xmacro['ion']==one_ion]
            if len(xxtab)>0:
                print(' z istate nlines nmacro: %3d %3d %10d %10d' % (one_z, one_ion, len(xxtab),len(xxmacro)))
    return



def doit(master='data/fe_23to27.dat'):
    '''
    Anaylyse a set of atomicdata for use
    with Sirocco, given a masterifle listing
    all of the files that consitute the inputs
    '''
    xmaster=xread(master)
    if xmaster==None:
        return
    zlines=[]
    for one in xmaster:
        try:
            zlines+=xread(one,split=True)
        except:
            print('Error: Could not read ',one)
            return
    # At this point we have all of the records, split into words
    
    element_tab=get_element_table(zlines)
    print('There are %d elements' % (len(element_tab)))
    
    ion_tab=get_ion_tab(zlines,element_tab)
    print('There are %d ions' % (len(ion_tab)))
    
    
    lev_tab=get_level_tab(lines=zlines,ions=ion_tab)
    lev_macro=lev_tab[lev_tab['macro']=='yes']
    
    print('There are %5d levels' % (len(lev_tab)))
    print('There are %5d macro_levels\n' % (len(lev_macro)))

    if len(lev_macro)>0:
       zz,nn=np.unique(lev_macro['z'],return_counts=True)
       j=0
       while j<len(zz):
           print('macro_levels z nlevels: %5d %10d' % (zz[j],nn[j]))
           j+=1
    print('\n')
       
    
    line_tab=get_line_tab(lines=zlines,lev=lev_tab)
    
    
    analyze_lines(line_tab=line_tab)
    
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__ )



