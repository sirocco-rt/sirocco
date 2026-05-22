#!/usr/bin/env python 

'''
This is a program that is designed to expand a .pf file into a grid of .pf files.


Usage:

    Use::

        pf_grid.py base.pf

    where base.pf is a correctly formatted parameter file, with variables that one
    wishes to grid indicated by ?

Description:  

    This program was written for creating grids of models with  
    python, although it should work with any .pf file in the case where one
    wants to run a grid of models.  As an example this should be a reasonable way
    to do much of the setup work for runing a bunch of ISM calculations.

    To use the program take a .pf file which is correctly formatted to run with python.
    This is going to be the basis of all of the .pf files.  Replace the values of those
    variables you want to grid with $

    Run the program. It will ask for the .pf file, herafter grid.pf, to use, 
    and a basename for all of the output files.  It will read through grid.pf and 
    will ask the user to give the values of each of the variables to grid. 

    The format of the inputs is very specific, and not particularly user-friendly! 
    (The lines are read with the input or using the raw_input command in python, 
    the scripting language'.  
    
    The two possibities that are supported by the program are:
    
    #. a list format which means one must enclose everything in [], e.g. [14,27,32]
       which would say that for this variable you want values of 14, 27, and 32.
    
    #. a special format to create a logarithminc grid of variables, e.g. log_ints(xmin,xmax,n)
       where xmin and xmax are the minimum and maximum values of the variable in
       question and n is the number of variables.  Note that you enter the
       actual values of xmin and xmax and not the log of this.  For example::

            log_ints(1.e16.,1.e18.,3)

       will end up produciing a list with variables [1.e16,1.e17,1.e18]

       Note also that this input is treated within the routine as a function,
       and so here the inputs have to be enclosed in parentheses.


Returns:

    The program produces a lot of files:

    a file which lists what was varied and the values for them.  The name of this file
    is the same as basename, with a .ls extensition attached, e.g vwhyi.ls
    This is intended for use with a fitting program::

        # Variable Disk.mdot(msol/yr)
        # Variable Wind.mdot(msol/yr)
        PYGRID/vwhyi_0000.spec     1e-08     1e-09 
        PYGRID/vwhyi_0001.spec     1e-08   2.2e-09 
        PYGRID/vwhyi_0002.spec     1e-08   4.6e-09 
        ... 
    
    a file which can run python for the grid, one model after another.  This has the prefix
    `Run_` and the remainder is the basename, e.g something like, `Run_vwhyi`.  It begins
    something like::

        #!/usr/bin/env bash
        py vwhyi_0000
        py vwhyi_0001
        py vwhyi_0002
        py vwhyi_0003
        py vwhyi_0100
        ...

    A number of .pf files, that are like the original pf files, but now the $ have been
    replaced with the variables associated with the grid.

Notes:

    If the program finds .pf files that seems to have been created earlier, it will ask 
    you if you want to delete them
                                       
History:

    0709    ksl
        Coded and debugged
    1701    ksl
        Updated for python3.  It should be backward compatible.  I have partially
        but not completely updated the style to the way I would write this program
        today.

'''

import sys
import math
import os
import fnmatch
from numpy import *


# Next routine is from kslio, just to make the installation slightly easier 
def get_input(question,answer):
    """ 
    normal usage answer=get_input(question,answer) 
    
    This is a simple command line io routine
    """

    string=question+' ('+answer+'):  '

    # The next if statements are needed so that the program will funciton in python 2.7 or 3
    if sys.version_info[0]<3:
        ss=raw_input(string)
    else:
        ss=input(string)
    if ss=='':
        ss=answer
    return ss
#End of code from kslio

def log_ints(xmin,xmax,nints):
    '''
    This generates a set of logarithmic intervals beween xmin and xmax inclusive
    '''

    lmin=math.log10(xmin)
    lmax=math.log10(xmax)
    dx=(lmax-lmin)/(nints-1)
    x=[]
    i=0
    while i < nints:
        q=lmin+dx*i
        q=pow(10.,q)
        x=x+[q]
        i=i+1
    return array(x)

def expand_array(arrayin,newvar,names_old):
    '''
    Expand the the arrays containing the items to ve
    varied to individual runs
    '''

    print(arrayin)
    print(newvar)
    print(names_old)

    names=[]
    if len(arrayin)==0:
        # Then we are just beginning
        outarray=newvar.reshape(len(newvar),1)
        nrows=len(newvar)
        n=0
        while n<nrows:
            names=names+["%02d"%(n)]
            n=n+1
    else:
        inshape=arrayin.shape
        rows_in=inshape[0]
        nmax=newvar.shape[0]
        rows_out=rows_in*nmax
        cols_in=inshape[1] 
        cols_out=cols_in+1 
        outarray=zeros([rows_out,cols_out])
        row=0
        nold=0
        while row<rows_out:
            nold=row//nmax
            fname=names_old[nold]
            fname="%s%02d"% (fname,row%nmax)
            names=names+[fname]
            col=0
            while col < cols_out:
                if col < cols_in:
                    outarray[row,col]=arrayin[nold, col]
                else:
                    outarray[row,col]=newvar[row % nmax]
                col=col+1
            row=row+1
    return outarray, names


# Here names and var_array are parallel arrays.
def export_results(lines, basename,names,var_array):
    '''
    Write out the results
    '''

    outfile=basename+".ls"
    gout=open(outfile,"w")

    # Write header

    n=0;
    while n<len(lines):
        z=lines[n].split()
        gout.write("# Variable %s\n" % z[0])
        n=n+1

    nrows=var_array.shape[0]
    ncols=var_array.shape[1]
    n=0
    while n<nrows:
        filename="PYGRID/%s_%s.spec" % (basename,names[n])
        gout.write("%s " % filename)
        m=0
        while m<ncols:


# Conditions to eliinate specific models ought to go here

# Finish the conditions
            gout.write("%9.2g "% var_array[n,m])
            m=m+1
        gout.write("\n")
        n=n+1


def create_parameter_files(basename,names,lines,z):
    '''
    Actually create the paameter files

    where:
    names are the number portions of the names
    lines contains the original parmeter fiels
    z contains the values

    '''
    n=0
    while n<len(names):
        pffilename=basename+'_'+names[n]+'.pf'
        pffile=open(pffilename,'w')
        m=0
        k=0
        while m<len(lines):
            words=lines[m].split()
            if len(words)>1 and words[1]=='$':
                pffile.write('%s        %s\n' % (words[0],z[n,k]))
                k=k+1
            else :
                pffile.write(lines[m])
            m=m+1
        pffile.close()
        n=n+1


def create_runfile(basename,names):
    '''
    Create a file which can be sourced to run all the
    models
    '''

    n=0
    runfile='Run_'+basename
    rfile=open(runfile,'w')
    rfile.write("#!/usr/bin/env bash\n")
    while n<len(names):
        pffilename=basename+'_'+names[n]
        rfile.write('py %s\n' % pffilename)
        n=n+1
    rfile.close()
    os.system('chmod +x '+runfile)


def cleanse(basename):
    '''
    Remove a preexistigng set of files with the 
    same name
    '''

    files=os.listdir('.')
    pattern=basename+'_*.pf'
    pffiles=fnmatch.filter(files,pattern)
    if len(pffiles)> 0:
        print("There are .pf files with this basename.")
        print(pffiles)
        i=int(get_input("Delete them -1, Quit 0, Continue 1","1"))
        if i==1:
            return
        if i==0:
            sys.exit()
        # Delete them
        os.system("rm "+pattern)


def doit(parameter_file=''):
    '''
    Generate a set of parameter files from a pre-existing one
    with at least one entry  ($) for interactive input

    '''


    if parameter_file=='':
        parameter_file=get_input("Starting.parameter.file","grid.pf")

    try:
        file_in=open(parameter_file,"r") 
        lines=file_in.readlines()
    except:
        print('Error: Could not read %s' % parameter_file)
        return


    # Define a basename for the grid
    basename=get_input("basename for grid","vwhyi")

    cleanse(basename)

    sumfile=basename+'.sum'
    sumf=open(sumfile,'w')

    print("Enter values for each variable, in bracketed format, e.g. [1, 4, 7,8]")
    print("Alternatively for logarithmic intervals, enter, something like: log_ints(1e-9, 2e-8, 5)")
    print("Be careful to use proper brackets or parentheses, as program will fail if you don't")

    i=0

    names=[]
    olines=[]
    z=[]
    while i<len(lines):
        line=lines[i]
        words=line.split()
        if len(words)>1 and words[1]=='$':
            olines=olines+[words[0]]
            string='%s [values]? ' % words[0]

            if sys.version_info[0]<3:
                zz=input(string)
            else:
                zz=eval(input(string))

            zz=array(zz)

        # Summarize the variables
            sumf.write('%s'%words[0])
            kk=0;
            while kk<len(zz):
                sumf.write('  %8.3g'%zz[kk])
                kk=kk+1
            sumf.write('\n')

            z,names=expand_array(z,zz,names)
        i=i+1


    # This creates a summary of all of the variables that are changed.  It is to be used in
    # fitting programs
    export_results(olines,basename,names,z)

    # This creates the individual .pf files
    create_parameter_files(basename,names,lines,z)

    # This creates the file that will run the grid of models
    create_runfile(basename,names)

    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__)

