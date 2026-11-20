#!/usr/bin/env python
'''
The purpose of this routine is to add a new parameter to an existing parameter file.

The routine uses a paremeter already existing in the parameter file to
find where to place the new parameter


Command line usage (if any):

    usage::

        add_param.py new_param value old_param

    where
        new_param
            is the new variable including what is normally included in
            parentheses.  This string will be written to the new .pf file.
        value
            is the default value to assign to the new parameter, which can be a string
        old_param
            is the the variable after which the new parameter is to
            be placed.  Note that this is a minimum match, which is to say
            that one does not need to include any more of the old parameter
            than what is required to make a unique identification

    The various paremeters may need to be enclosed in quotes.

Description:

    The routine first searches the parameter file to see whether the new
    parameter exists, if it does then it does not create a new .pf file

    If that is not the case then the routine puts the new parameter with its
    value after the old_param.

    If a new paremeter file is merited it is written to a file with a prefix
    of `new_`

    If new parameter files are written a command file MoveEm is created, which
    can be run by "source MoveEm" to replace the old .pf files with the new ones.

Primary routines:

    doit- the master routine for doing all of the pf files in a directory
    do_one - the routine that processes each pf file

Notes:

    This is not intended to be a bulled-proof routine so one should use it
    with care, inspecting at least some of the new files before replacing the
    old parameter files is strongly advised.

History:

    190803 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy
from glob import glob


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments

    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:

    110729    ksl    Added optional delimiters
    141209    ksl    Reinstalled in my standard startup script so there was flexibility to read any ascii file
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

def do_one(filename,new_param,value,old_param):
    '''
    Precess a single file


    The routine returns True if a new .pf file has been written, False
    otherwise
    '''

    try:
        f=open(filename)
        lines=f.readlines()
    except:
        print('Error do_one: %s does not appear to exist' % filename)
        return False

    # Check that the new parameter does not exist in the file a already
    # Also check that the parameter we want to put the line behind exits

    words=new_param.split('(')
    basic=words[0]

    old_exists=False
    for line in lines:
        if line.count(basic):
            print('Found %s in file %s already. Nothing to do so returning' % (basic,filename))
            return False
        if line.count(old_param):
            old_exists=True

    if old_exists==False:
        print('Warning: Could not find the old_parameter name in the file %s,  OK if not expected' % filename)
        return False


    g=open('new_'+filename,'w')

    write_it=True
    for line in lines:
        g.write(line)
        if line.count(old_param) and write_it==True:
            g.write('%40s   %s\n' % (new_param,value))
            write_it=False # Prevent writing the parameter multiple times

    g.close()
    return True







def doit(new_param,value,old_param):
    '''
    Do something magnificent

    Description:

    This routine calls do_one many times.  If a new file has been created
    if creates a file which can be used to move the new files onto the
    old ones

    Notes:

    History:


    '''

    files=glob('*.pf')

    command_file=False

    for one in files:
        ireturn=do_one(one,new_param,value,old_param)
        if ireturn==True:
            if command_file==False:
                g=open('MoveEm','w')
                command_file=True
            g.write ('mv new_%s %s\n' % (one,one))

    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)==4:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print (__doc__)
