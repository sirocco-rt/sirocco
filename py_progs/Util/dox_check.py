#!/usr/bin/env python 

'''
Determine which files lack doxygen commments


Command line usage (if any):

    usage::

        dox_check.py filename

Description:  

Primary routines:

    `doit`

Notes:
                                       
History:

    180404 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy

from glob import glob

def doit():
    '''
    Simply check whether .c files appear
    to have doxygen comments or not.

    Description:

    Notes:

    History:


    '''

    files=glob('*.c')

    print(files)

    missing=[]
    ngood=0

    for one in files:
        f=open(one)
        # print(one)
        lines=f.readlines()
        n=0
        for line in lines:
            n+=line.count('@brief')
            n+=line.count('@param')
        print(one,n)
        if n==0:
            missing.append(one)
        else:
            ngood+=1
           

    print('These files do not have doxygen comments:')
    x=open('dox_missing.txt','w')
    for one in missing:
        print(one)
        x.write('%s\n' % one)
    x.close()

    print('\nOf the %d .c files in this directory, \n %3d have doxygen comments, \n %3d do not'
            % (len(files),ngood,len(missing)))
    
    



    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    doit()
