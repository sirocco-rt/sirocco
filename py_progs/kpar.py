#!/usr/bin/env python 

'''
A python version of kpar


Command line usage (if any):

    usage: kpar.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

    160704 ksl
        Coding begun

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
        print("The file %s does not exist" % filename)
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
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        return

    print('Here are the column names:')
    
    print(data.colnames)

    return data


def opar(filename):
    '''
    Read a parameter file
    '''
    global xxx
    global ixxx

    xxx=read_file(filename)
    ixxx=0
    return



def rdpar(string):
    '''
    read a parameter
    '''
    global ixxx
    global xxx

    while ixxx<len(xxx):
        word=xxx[ixxx][0].split('(')
        if word[0]==string:
            try:
                z=eval(xxx[ixxx][1])
            except NameError:
                z=xxx[ixxx][1]
            ixxx=ixxx+1
            return z
        ixxx=ixxx+1
    if ixxx==len(xxx):
        ixxx=0

    while ixxx<len(xxx):
        word=xxx[ixxx][0].split('(')
        if word[0]=='string':
            try:
                z=eval(xxx[ixxx][1])
            except NameError:
                z=xxx[ixxx][1]
            ixxx=ixxx+1
            return z
        ixxx=ixxx+1
    

    answer=input(string+': ')

    try:
        z=eval(answer)
    except NameError:
        z=answer


    xxx.append([string,answer])

    return z






def cpar(filename):
    '''
    Write a parameter file
    '''
    g=open(filename,'w')
    for one in xxx:
        g.write('%-40s %-20s\n' % (one[0],one[1]))
    g.close()
    



    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print('usage: kpar.py filename')
