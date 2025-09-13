#!/usr/bin/env python 

'''
Create a 1-d spherical `.pf` file from a `windsave2table` file.

Read the master file produced by windsave2table for a 1d spherical
model and produce a file which can be imported into Sirocco


Command line usage (if any):

    usage: import_1d.py root     

    where root is the rootname of the mastertable

Description:  

   The routine is intended to provide examples of the types
   of files Sirocco can read, and to allow regression tests
   for importing models

Primary routines:

    doit

Notes:

    Windsave2table saves the values of rho in the CMF frame,
    but Sirocco  expects models to be in the observer frame
    so this routine corrects for this.
                                       

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


def doit(root='star',outputfile=''):
    '''
    Read the table produced by running wind
    save2table on a model for a spherical wind
    and produce a table which python can import


    Description:

    Notes:

    History:

    '''

    filename=root+'.master.txt'
    if outputfile=='':
        outputfile=root+'.import.txt'
    


    data=read_table(filename)

    v=data['v_x']**2+data['v_y']**2+data['v_z']**2
    v=numpy.sqrt(v)

    data['v']=v

    xdata=data['i','inwind','r','v','rho','t_e','t_r']

    C=2.997925e10

    gamma=1./numpy.sqrt(1-(xdata['v']/C)**2)
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
