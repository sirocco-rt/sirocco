#!/usr/bin/env python 

'''
Checks .sig files to see which runs in a grid have completed.

When running a grid of models, it may not be clear whether
all of the runs have been completed, especially if time limits
are placed on the individual runs.  This routine checks the
.sig files to see which of a list of .pf files have been
completed.


Command line usage (if any):

    usage::

        grid_check.py filename

    where filename contains a list of the .pf files that need to be checked

Description:  

    The routine reads a file containing a list of .pf files, excluding those
    with out.pf extensions and then checks the .sig files to see whether
    the run is complete.  
    
    It generates a table Status.txt that contians the completion 
    status of each .pf file and writes to the screen some information about
    those that do not appear to be complete.

    Additionally, if there are runs that look to be incomplete it writes out
    a file XRunRest that contains a proto run file for the incompleted runs.  
    This will need to be edited for each specific situation. 


Primary routines:

    doit

Notes:
                                       
History:

190516 ksl Coding begun

'''

import sys
from astropy.io import ascii
from astropy.table import Table
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



def doit(filename='all.txt',outputfile='status.txt'):
    '''
    Read a file containing a list of .pf files and 
    determine whether a run is complete.

    Description:

    Notes:

    History:


    '''

    records=read_file(filename)
    if len(records)==0:
        print('Exiting because file containing list of .pf files does not exist')
        return

    # read_file returns the file names words for each record,
    # We want to collect the original.pf files and make sure
    # .out.pf are not included

    pf_files=[]
    for one in records:
        if one[0].count('.out.pf')==0 and one[0].count('pf'):
            pf_files.append(one[0])

    print('There are %d runs to look at' % (len(pf_files)))

    x=Table([pf_files],names=['pf_file'])

    x['Complete']='    Unknown'
    x['RunTime']=-99.
    x['RunTime'].format='.1f'

    for one in x:
        try:
            xfile=one['pf_file'].replace('.pf','.sig')
            f=open(xfile)
            lines=f.readlines()
            last=lines[len(lines)-1]
            words=last.split()
            if words[6]=='COMPLETE':
                one['Complete']='yes'
            else:
                one['Complete']='no'
            one['RunTime']=eval(words[5])
        except:
            print('%s was not found' % xfile)
            one['Complete']='NotBegun'

    x.write(outputfile,format='ascii.fixed_width_two_line',overwrite=True)
    print('Evaluations details are reported in %s' % outputfile)

    xbad=[]
    i=0
    while i<len(x):
        if x['Complete'][i] != 'yes':
            xbad.append(i)
        i+=1

    if len(xbad)==0:
        print('All of the %d runs listed in %s appear to be complete!' % (len(x),filename))
        return

    xx=x[xbad]


    print('These runs are incomplete')

    g=open('XRunRest.txt','w')

    for one in xx:
        print('%20s %15s  %.1f' %(one['pf_file'],one['Complete'],one['RunTime']))
        g.write('py %s' % one['pf_file'])

    print('In summary, of %d runs, %d have not been completed' % (len(x),len(xx)))
    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print (__doc__)
