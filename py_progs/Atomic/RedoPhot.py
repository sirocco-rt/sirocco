#!/usr/bin/env python 

'''
Extend phot tables retrieved from Topbase to higher energies and produce a plot file which shows the extended x-section file

Command line usage (if any):

    usage: RedoPhot.py [-h] [-outroot whatever] filename

Description:  

    where:

    -h prints the the help information and then the routine
    quits

    -outroot whatever is the root name for the output files
    If this optional parameter is not given the ouputs
    will be written to test_phot, and the figure will 
    have the name test_phot.png

    filename is the name of a macro atom phot file

Primary routines:

    redo_one

Notes:

    The input file must be a standard photfile usable by Python.

    At present this is only set up for handling phot files for
    MacroAtoms but it would be easy to modify it to work 
    with photfiles intended to work with simple ions
                                       
History:

231102 ksl Coding begun

'''

import sys
from astropy.io import ascii
import os
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


def read_phot(photfile='o_2_phot.dat'):
    '''
    Read a valid Python photometry file
    and split the inputs up into a
    table that contains the header lines
    for each set of x-sections and 
    a list of tables, one table for
    each transition.
    '''
    
    try: 
        f=open(photfile)
    except:
        print('Error: RedoPhot.read_phot: %s not found' % photfile)
        return [],[]

    rlines=f.readlines()
    
    words=[]
    for line in rlines:
        xline=line.strip()
        word=xline.split()
        words.append(word)
        
    # At this point we have everything in words
    
    lineno=[]
    label=[]
    z=[]
    state=[]
    lower_lev=[]
    upper_lev=[]
    e_thresh=[]
    npt=[]
    i=0
    
    xhead=[]
    xtab=[] # list of tables with the xsections
    for word in words:
        # print(word)
        if word[0].count('PhotMacS') or word[0].count('PhotTopS'):
            # print('Got header')
            lineno.append(i)
            label.append(word[0])
            z.append(int(word[1]))
            state.append(int(word[2]))
            lower_lev.append(int(word[3]))
            upper_lev.append(int(word[4]))
            e_thresh.append(eval(word[5]))
            npt.append(int(word[6]))
            
            if len(xhead)>0:
                cross_tab=Table([xhead,xev,xsec],names=['Label','e','sigma'])
                xtab.append(cross_tab)
            
            xhead=[]
            xev=[]
            xsec=[]
        elif word[0]=='PhotMac' or word[0]=='Phot':
            # print('Got data')
            xhead.append(word[0])
            xev.append(eval(word[1]))
            xsec.append(eval(word[2]))
        i+=1
    phot_tab=Table([label,z,state,lower_lev,upper_lev,e_thresh,npt,lineno],names=['Label','z','ion','ll','ul','e_thresh','npt','LineNo'])
    if len(xhead)>0:
            cross_tab=Table([xhead,xev,xsec],names=['Label','e','sigma'])
            
            xtab.append(cross_tab)
    if len(phot_tab)==0:
        print('Error: %s exists but had no x-section headers' % photfile)
        return [], []
    return phot_tab,xtab    


def plot_phot(phot_sum,xcross):
    '''
    Plot each of the x-sections in xcross
    '''
    plt.figure(1,(8,8))
    i=0
    while i<len(xcross):
        one_tab=xcross[i]
        plt.loglog(one_tab['e'],one_tab['sigma'],label='%2d' % (i+1))
        i+=1
    plt.xlabel('Energy (ev)',size=16)
    plt.ylabel(r'$\sigma$',size=16)
    plt.legend()
    plt.tight_layout()



def extrap(photfile='o_6_phot.dat',emax=1e5):
    '''
    This routine adds points to the various
    x-sections so that the maxium energy
    is given by emax

    It returns the results in terms of one
    table that contains header information
    for the x-sections, and a list of tables
    one for each transition.
    '''
    phot_tab,xtab=read_phot(photfile)
    if len(phot_tab)==0:
           return [],[]

    k=0
    for one in xtab:
        k+=1                 
        if one['e'][-1]>=emax:
            print('The x-section for the %3d transition already extends to %e >= %e' % (k,one['e'][-1],emax))
            continue


        sigma_start=sigma_jump=one['sigma'][0]
        jjump=0
        j=1
        while j<len(one):
            if one['sigma'][j]>sigma_start:
                sigma_jump=one['sigma'][j]

                jjump=j
            sigma_start=one['sigma'][j]
            j+=1
        
        jjjump=jjump
        if len(one)-jjjump>10:
            jjjump=len(one)-10
            
        a=(np.log10(one['sigma'][-1])-np.log10(one['sigma'][jjjump]))/(np.log10(one['e'][-1])-np.log10(one['e'][jjjump]))
        b=np.log10(one['sigma'][-1])
        
        # print(a,b)
        
        e_more=np.logspace(np.log10(one['e'][-1]),np.log10(emax),10)
        # print(e_more)
        log_sigma=a*(np.log10(e_more)-np.log10(one['e'][-1]))+b
        sigma_more=10**log_sigma
        
        i=1  # We don't want to add the first value
        while i<len(e_more):
            one.add_row([one['Label'][0],e_more[i],sigma_more[i]])
            i+=1
                   
        # print('%.6e %5d %5d %e %e %e %e' % (sigma_jump,jjump,len(one),np.log10(one['sigma'][jjump]),np.log10(one['sigma'][-1]),a,b))

    return phot_tab,xtab
 

def write_phot_tab(out_name,xsum,xcross):
    '''
    Write the tables that have been used as port of the
    analysis out into a valid Phot file for use with
    python.  

    Note that this routine is not prevented from
    writing data to an existing file.
    '''
    print('RedoPhot: Starting to write phot file to %s' % (out_name))

    f=open(out_name,'w')
    i=0
    while i<len(xsum):
        one=xsum[i]
        string='%s %2d %2d %2d %d %12.6f %4d' % (one['Label'],one['z'],one['ion'],one['ll'],one['ul'],one['e_thresh'],len(xcross[i]))
        # print(string)
        f.write('%s\n' % string)
        for one_x in xcross[i]:
            # print(one_x)
            string='%s %12.4f %12.4e' % (one_x['Label'],one_x['e'],one_x['sigma'])
            # print(string)
            f.write('%s\n' % string)
        i+=1

    f.close()
    print('RedoPhot: Write extended phot file to %s' % (out_name))




def redo_one(phot_file='o_2_phot.dat',outroot=''):
    phot_tab,xtab=extrap(phot_file,1e5)
    if len(phot_tab)==0:
        print('Error: RedoPhot - Exiting because of previous problems')
        return

    plot_phot(phot_sum=phot_tab,xcross=xtab)
    if outroot=='':
        word=phot_file.split('.')
        if len(word)>1:
            outroot=word[-2]
        else:
            outroot=phot_file
    
    
    plt.savefig(outroot+'_phot.png')
    write_phot_tab(outroot+'_phot.dat',phot_tab,xtab)
    print('Finished: The maximum number of x-section values for any transition is %3d' % (np.max(phot_tab['npt'])))
    return

def steer(argv):
    '''
    This is just a steering routine
    '''

    outroot=''
    infile=''
    i=1
    while i<len(argv):
        print(argv)
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-outroot':
            i+=1
            outroot=argv[i]
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s' % argv[i])
            return
        elif infile=='':
            infile=argv[i]
        else:
            print('Error: Extra argument %s' % argv[i])
        i+=1

    if infile=='':
        print('Error: not enough arguments: ',argv)
        return

    redo_one(infile,outroot)
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
