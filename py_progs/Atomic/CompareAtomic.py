#!/usr/bin/env python
# coding: utf-8

'''
Simple routine to determine differences between files read in two sets of atomic data.


Command line usage (if any):

    usage: CompareAtomic.py masterfile1 mastefile2

Description:  

    Routine parses each of two master files, eliminating
    all the comments and blanck lines and simply looks
    to see what files are read in by each masterfile.

    It checks what files would be read in from the
    first file, that are not in the second file,
    and then caries out the chck in the opposite 
    direction.

    The results, which include the line no of
    a file in the original masterfiles 
    are printed to the screen.

Primary routines:

    doit

Notes:

The routine does not tell one if the files that are
read in are actually used, just that they will be
read in if Python uses this atomic data file

                                       
History:

240510 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np
from astropy.table import Table,join


def xread(filename='data/fe_23to27.dat',split=False):
    '''
    Read a set of atomic data, ignoing comments
    and blank lines.  Optionally  split the lines
    into words
    '''
    try:
        x=open(filename,'r')
    except:
        print('Error: could not open ',filename)
        return []
    lines=x.readlines()
    records=[]
    line_no=[]
    i=1
    for one in lines:
        z=one.strip()
        if len(z)>0 and z[0]!='#':
            line_no.append(i)
            if split:
                records.append(z.split())
            else:
                records.append(z)
        i+=1

    xtab=Table([records,line_no],names=['File','LineNo'])
    return xtab



def compare_files(one='xdata/h10_hetop_standard80.dat',two='zdata/master_hhe.dat'):
    
    master_one=xread(one)
    master_two=xread(two)

    if len(master_one)==0  or len(master_two) == 0:
        return
    
    print('\n---------------\n')
    print('Files in %s that are not in %s' % (one,two))
    print('\n---------------\n')

    x1=np.setdiff1d(master_one['File'],master_two['File'])
    x1=Table([x1],names=['File'])
    x1=join(x1,master_one,join_type='left')
    x1.sort(['LineNo'])

    x1=join(x1,master_one,join_type='left')

    for one_file in x1:
        print('%5d %s' % (one_file['LineNo'],one_file['File']))
    
    print('\n---------------\n')
    print('Files in %s that are not in %s' % (two,one))
    print('\n---------------\n')

    x2=np.setdiff1d(master_two['File'],master_one['File'])
    x2=Table([x2],names=['File'])
    x2=join(x2,master_two,join_type='left')
    x2.sort(['LineNo'])


    for one_file in x2:
        print('%5d %s' % (one_file['LineNo'],one_file['File']))
    
    print('\n---------------\n')
    return


def doit(one,two):
    compare_files(one,two)
    return
        
    

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)==3:
        doit(sys.argv[1],sys.argv[2])
    else:
        print(len(sys.argv))
        print (__doc__)
