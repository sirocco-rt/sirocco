#!/usr/bin/env python 

'''
Create a  polar coordinate `.pf` file from a `windsave2table` file.

Read the master file produced by windsave2table for a
rtheta (polar-coordinate mode model and produce 
a file to be read by import in python


Command line usage (if any):

    usage: import_rtheta.py rootname

    e.g rtheta to read a the master.txt table

    rtheta.0.master.txt

Description:  

    This routine allows one to create examples of files
    that can be imported into Python.  It is intended
    to be useful for creating examples, and for regression
    testing

Primary routines:

    doit

Notes:

    Windsave2table saves the values of rho in the CMF frame,
    but Python expects models to be in the observer frame
    so this routine corrects for this.
                                       
History:

171106 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
        110729    ksl
            Added optional delimiters
        141209    ksl
            Reinstalled in my standard startup
            script so there was flexibility to
            read any ascii file
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print ("The file %s does not exist" % filename)
        return []   
    
    lines=[]
    
    i=0
    while i<len(xlines):
        z=xlines[i].strip()
        if char=='':
            z=z.split()
        else:
            z=z.split(char)
        if len(z)>0:
            if z[0][0]!='#':
                lines=lines+[z]
        i=i+1
    return lines




def read_table(filename='foo.txt',format=''):
    '''
    Read a file using astropy.io.ascii and 
    return this 

    Description:

    Notes:

    History:


    '''
    try:
        if format=='':
            data=ascii.read(filename)
        else:
            data=ascii.read(filename,format=format)
        for col in data.itercols():
            if col.dtype.kind in 'SU':
                data.replace_column(col.name,col.astype('object'))
    except IOError:
        print ('Error: file %s does not appear to exist' % filename)
        return

    print ('Here are the column names:')
    
    print (data.colnames)

    return data


def doit(root='rtheta',outputfile=''):
    '''
    Read a master.txt file produced by windsave2table for
    an rtheda coordinate model and produce a file
    which can be imported into  python
    

    Description:

    Notes:

    History:

    '''

    filename=root+'.master.txt'
    if outputfile=='':
        outputfile=root+'.import.txt'
    


    data=read_table(filename)


    xdata=data['i','j','inwind','r','theta','v_x','v_y','v_z','rho','t_r']

    C=2.997925e10

    v=numpy.sqrt(xdata['v_x']**2+xdata['v_y']**2+xdata['v_z']**2)

    gamma=1./numpy.sqrt(1-(v/C)**2)
    xdata['rho']*=gamma
    xdata['rho'].format='.2e'
    
    print (xdata)


    # This format is the easy to read back automatically
    ascii.write(xdata,outputfile,format='fixed_width_two_line',overwrite=True)

    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__)
