#!/usr/bin/env python 

'''
Compare two regression runs, highligting the differences between them


Command line usage (if any):

    usage: regression_check.py [-h] run1 [run2]

    where run1 and run2 are the names of two diretories where
    a standard set of regression tests have been run, and where
    -h prints out this message.

Description:  

    The routine checks various files in two regression direct to
    see how different the files are.  The actual diffs (or the
    beginning of them anyway) are written to a file (check.txt by
    default).

    If run1 and run2 are provided then these will be the two
    directories that are compared.

    If only run1 is provided then the comparison will be makde between
    the run1 directory and the last run of the regression test, excluding
    run1.  

    The routine looks for all of the .out.pf files in the two directories,
    and compares them, and then makes comparisons of the .spec and .log_spec_tot
    files.  

    

Primary routines:

    doit

Notes:
                                       
History:

180810 ksl Coding begun

'''

import sys
from astropy.io import ascii
from astropy.table import Table,join
import numpy
import os
from glob import glob

import difflib  # To compare text files
import regression_plot # To make plots


def diff_two_files(file1,file2):
    '''
    Report on the differences between two files
    '''


    try: 
        x1=open(file1)
        text1=x1.readlines()
    except IOError:
        print('Could not open %s' % (file1))
        raise ValueError('Missing file')


    try: 
        x2=open(file2)
        text2=x2.readlines()
    except IOError:
        print('Could not open %s' % (file2))
        raise ValueError('Missing file')
        return

    ## d = difflib.Differ()
    ## diff = d.compare(text1,text2)
    # diff=difflib.ndiff(text1,text2)
    # diff = difflib.Differ.compare(text1,text2)
    diff=difflib.unified_diff(text1,text2) # this is supposed to be faster

    keep=[]
    for e in diff:
        if e.startswith('+'):
            keep.append(e)
        if e.startswith('-'):
            keep.append(e)

    return keep


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729    ksl    Added optional delimiters
    141209    ksl    Reinstalled in my standard startup
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


def doit(run1='py_180809',run2='',outputfile='check.txt'):
    '''
    Run the regression checks on two directories.

    Description:

    If run2 =='', then find the second directory  to compare 2

    Notes:

    History:


    '''

    if run2=='':
        run2=get_other_directory(run1)

    out=open(outputfile,'w')

    if os.path.isdir(run1)==False:
        print('Error: %s is not a directory' % run1)
    if os.path.isdir(run2)==False:
        print('Error: %s is not a directory' % run2)

    pf1=glob('%s/*.out.pf' % run1)
    pf2=glob('%s/*.out.pf' % run2)

    # print(pf1)
    # print(pf2)

    name1=[]
    root1=[]
    for one in pf1:
        x=one.replace('.out.pf','')
        root1.append(x)
        x=x.split('/')
        name1.append(x[1])
    # print(name1)

    name2=[]
    root2=[]
    for one in pf2:
        x=one.replace('.out.pf','')
        root2.append(x)
        x=x.split('/')
        name2.append(x[1])
    # print(name2)


    table1=Table([name1,root1],names=['name','root1'])
    table2=Table([name2,root2],names=['name','root2'])

    try:
        combined=join(table1,table2,join_type='left',keys=['name'])
    except ValueError:
        print('Error: regression_check: run1 %s len %d run2 %s len %d\n' % (run1,len(table1),run2,len(table2)))
        return

    # print(combined)

    # Check various files that have been produces

    pf_count=[]
    log_spec_tot_count=[]
    spec_count=[]
    for one in combined:
        out.write('\nCOMPARING %s\n' % one['name'])
        print('COMPARING %s' % one['name'])

        ext='.out.pf'
        # x1=one['root1']+ext
        # x2=one['root2']+ext

        try:
            x1=one['root1']+ext
            x2=one['root2']+ext
            result=diff_two_files(x1,x2)
            # print('There were differences in %5d lines in %s' % (len(result)//2,ext))
            out.write('There were differences in %5d lines in %s\n' % (len(result)//2,ext))
            if len(result)<50:
                out.write(''.join(result))
            pf_count.append(len(result)//2)
        except ValueError:
            pf_count.append(-99)
        except TypeError:
            pf_count.append(-99)


        ext='.log_spec_tot'
        # x1=one['root1']+ext
        # x2=one['root2']+ext

        try:
            x1=one['root1']+ext
            x2=one['root2']+ext
            result=diff_two_files(x1,x2)
            # print('There were differences in %5d lines in %s' % (len(result)//2,ext))
            out.write('There were differences in %5d lines in %s\n' % (len(result)//2,ext))
            log_spec_tot_count.append(len(result)//2)
            if len(result)<50:
                out.write(''.join(result))
        except ValueError:
            log_spec_tot_count.append(-99)
        except TypeError:
            log_spec_tot_count.append(-99)



        ext='.spec'
        # x1=one['root1']+ext
        # x2=one['root2']+ext

        try:
            x1=one['root1']+ext
            x2=one['root2']+ext
            result=diff_two_files(x1,x2)
            # print('There were differences in %5d lines in %s' % (len(result)//2,ext))
            out.write('There were differences in %5d lines in %s\n' % (len(result)//2,ext))
            spec_count.append(len(result)//2)
            if len(result)<50:
                out.write(''.join(result))
        except ValueError:
            spec_count.append(-99)
        except TypeError:
            spec_count.append(-99)


    combined['pf_count']=pf_count
    combined['log_spec_tot_count']=log_spec_tot_count
    combined['spec_count']=spec_count

    print('\nThe table below shows differences in various files in regression directories %s and %s.' % (run1,run2))
    print('More details can be found in %s.\n' % outputfile )


    print(combined)

    print('\n Make plots of the spectra which will be stored in Xcompare\n')
    regression_plot.do_all(run1,run2)
    print('Plotting completed')



    return


def get_other_directory(run1):
    '''
    Find a second directory.  The second directory is assumed to be
    the regression directory that was most recently created that is not
    run1
    '''
    x=glob('py*')
    dirs=[]
    modtime=[]
    for one in x:
        one=one.rstrip('/')
        if os.path.isdir(one) and one!=run1:
            dirs.append(one)
            modtime.append(os.path.getmtime(one))


    dirs=numpy.array(dirs)
    modtime=numpy.array(modtime)
    iorder=numpy.argsort(modtime)


    if len(dirs)>1:
        run2=dirs[iorder[len(dirs)-1]]
        return run2
    else:
        print('Could not find a second set of tests for a comparison')
        return ''

        

def steer(argv):
    '''
    Parse the command line and set everything up
    '''
    run1=''
    run2=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif run1=='' and os.path.isdir(argv[i]):
            run1=argv[i]
        elif run2=='' and os.path.isdir(argv[i]):
            run2=argv[i]
        else:
            print('Error: Uninterpretable parameter, check dir names')
        i+=1

    # Now deal with the case where we only got one entry

    if run2=='':
        run2=get_other_directory(run1)
        if run2=='':
            return

    # Check that we have everthing we need and then run the check

    if os.path.isdir(run1) and os.path.isdir(run2):
        print('Comparing %s to %s' % (run1,run2))
        doit(run1,run2)
    else:
        print('Not able to parse input line:',''.join(argv))

    return






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
