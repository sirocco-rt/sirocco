#!/usr/bin/env python 

'''
Check whether c routines have proper dox headers, and if not write a new file with a dummy dox header for the user to update


Command line usage:

    usage:
    * `dox.py whatever.c` to process a single file
    * `dox.py *.c`  or `-all` to process all c files in a directory
    * `dox.py -h` to get this help message

Description:  

    This script looks for a doxyegen header for each function in a .c file, and 
    reports on the results.  
    
    If one or more headers are missing, for a file named for exampe foo.c,  it will write out
    a new file named `new_foo.c` with a basic header installed.

    This basic headers in this routine should be edited and the file, in this case `new_foo.c`
    copied back to foo.c

Primary routines:

    `doit` processes a single file
    `do_many` processes all of the files in a directory, by calling doit multiple times
    steer processes the command line

Notes:

    Files that begin `new_` will not be checked

    dox.py can give incorrect results if the routine is not properly indented

    This routine does not check for old style headers used prior to the introduction of doxygen
    to Python, since that conversion should be complete.
                                       
History:

    180912 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy
import os
import subprocess
import os
from glob import glob


def is_installed(program):
    """
    Tests to see if a program has been installed. Code taken from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Args:
        program (str): Program name to test
    Returns:
        Either returns True or raises an OSError
    """
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    raise(OSError("Executable '{}' could not be found!".format(program)))


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

file_header_string='''

/***********************************************************/
/** @file   %s
 * @author ??? initials of primary author or authors, last 
 * should be current maintainer
 * @date   ??? e.g January, 2018
 * @brief  ???? brief description of the purpose of the routine
 * in this file
 *
 * ???  Any other relevant discussion
***********************************************************/

'''


def gen_file_header(filename):
    '''
    Craete a new header for the file since one was not found
    '''

    header=file_header_string % (filename)

    return header

header_start='''
/**********************************************************/
/**
 * @brief   ????? (Does not include function name)
'''

header_end='''*
 * @return ????
 * 
 * Explain here what is returned, not only by the return 
 * statement but also to passed variables.
 *
 * @details     
 * More extended explantion of what the routine does, and how it 
 * is used.  This section is mandatory, 
 *
 * Include here an references to other work, e.g a paper, or
 * someone's thesis.
 *
 * ### Notes ###
 * 
 * This section is optional, and is intended primarily
 * to help maintainers or developers.
 * 
 * This can contain a discussion of how the routine might
 * be improved.  If so, leave initials and date to identify
 * who made the comment.
 *
 * This is not intended to be a history of the routine, but
 * can include some references to the latest changes
 *
 * Remove the section title if it is not used at all.
**********************************************************/
'''

def gen_header(xtype='int',xname='get_domain_params',xvar=['int ndom','double *foo']):
    '''
    generate a doxygen header for function


    '''

    header=header_start
    for one in xvar:
        if one=='void':
            pass
        else:
            string='* @param  [in,out] %s  ?????\n' % one
            header+=string
    header+=header_end
    return header




def doit(filename='setup_domains.c',outputfile=''):
    '''
    Do something magnificent

    Description:

    Notes:

    History:


    '''

    try:
        is_installed('cproto')
    except:
        print('cproto is not installed')
        return



    try:
        f=open(filename)
        lines=f.readlines()
    except IOError :
        print ("The file %s does not exist" % filename)
        return    

    print('\nAnalyzing: %s' % filename)

    proc = subprocess.Popen('cproto %s ' % filename, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split('\n')
    # print('Knox1\n',stdout)


    # Parse the cproto outputs       

    xtype=[]
    xname=[]
    xvar=[]

    iprog=1
    while iprog<len(stdout)-1:
        one_prog=stdout[iprog]
        # First split out the function definition from the variables
        open_paren=one_prog.index('(')
        # print('Prog: ',one_prog[:open_paren])
        # print('Var : ',one_prog[open_paren:])
        word=one_prog[:open_paren]
        words=word.split()
        xtype.append(words[0])
        xname.append(words[1])
        word=one_prog[open_paren:]
        word=word.replace('(',' ')
        word=word.replace(')',' ')
        word=word.replace(';',' ')
        words=word.split(',')
        k=0
        while k<len(words):
            words[k]=words[k].strip()
            k+=1
        xvar.append(words)
        iprog+=1
    # print(xtype)
    # print(xname)
    # print(xvar)


    # Find the start of each routine

    i=0
    start=[]
    while i<len(xname):
        print('Checking:',xname[i])
        line_no=-1
        icheck=0
        j=1
        while j<len(lines):
            long_string=lines[j-1]+lines[j]
            words=long_string.split()
            k=1
            while k<len(words):
                if words[k-1]==xtype[i] and words[k]==xname[i]:
                        icheck+=1
                        if icheck>1:
                            print('Error: Found multiple lines for %s %s' % (xtype[i],xname[i]))
                        line_no=j-1
                k+=1
            j+=1
        start.append(line_no)
        i+=1

    print('Starting lines:',start)

    # Find the end of each routine, assuming that the last } line before the next open module

    j=0
    end_brackets=[]
    while j<len(lines):
        if lines[j].count('}')>0:
            end_brackets.append(j)
        j+=1

    # print(end_brackets)

    #  Assume that the close bracket is the last one before the begining of another routine

    i=0
    end=[]
    while i<len(start)-1:
        for one in end_brackets:
            if one<start[i+1]:
                xend=one
            else:
                break
        end.append(xend)
        i+=1
    # Now add the last one
    end.append(end_brackets[len(end_brackets)-1])

    print('Ending lines  :',end)


    # Now check for the existence of a doxygen header
    # First try to find the file header

    print('Checking for file header')
    i=0
    header_end=0
    file_header=False
    while i<end[0]:
        if lines[i].count('@file'):
            file_header=True
        if file_header==True and lines[i].count('*****************/'):
            print('Results: Found file header for %s' % filename)
            header_end=i
            break
        i+=1
    if file_header==False:
        print('Results: Could not locate overall file header for %s' % filename)

    header=[]

    i=0
    while i<len(xtype):
        brief=False
        if i==0:
            jstart=header_end
        else:
            jstart=end[i-1]
        jstop=start[i]

        j=jstart
        while j<jstop:
            if lines[j].count('@brief'):
                brief=True
            j+=1
        header.append(brief)
        i+=1

    # print(header)

    all_ok=True
    if file_header==False:
        all_ok=False
    i=0
    while i<len(header):
        if header[i]==False:
            all_ok=False
            print('\nResults: routine %s failed dox check' % xname[i])
        i+=1

    if all_ok==True:
        print('Results: %s passed all doxchecks' % filename)
    else:
        if outputfile=='':
            outputfile='new_'+filename
        print('Results: %s has dox issues; writing new file %s' % (filename,outputfile))
        f=open(outputfile,'w')

        if file_header==False:
            q=gen_file_header(filename)
            f.write(q)

        while i<len(lines):
            k=0
            while k<len(xname):
                if start[k]==i and header[k]==False:
                    q=gen_header(xtype[k],xname[k],xvar[k])
                    f.write(q)
                k+=1
            f.write(lines[i])
            i+=1

    return




def do_many():
    ''' 
    Check all of the .c files in a directory for dox commens
    '''

    files=glob('*.c')

    for filename in files:
        if 'new_'!=filename[:4]:
            doit(filename)

    return


def steer(argv):
    '''
    Process the command line
    '''

    files=[]
    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        if argv[i]=='-all':
            do_many()
            return
        else:
            files.append(argv[i])
        i+=1

    todo=[]
    for one in files:
        if one[:4]!='new_' and one[len(one)-2:]=='.c':
            todo.append(one)

    if len(todo)==0:
        print('Error: Read command line but could find nothing to do')
        return

    for one in todo:
        doit(one)

    return





    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        # doit(sys.argv[1])
        steer(sys.argv)
    else:
        print ('usage: dox.py filename')
