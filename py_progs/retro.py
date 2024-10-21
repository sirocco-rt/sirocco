#!/usr/bin/env python 

'''
Run multiple commits on the same `.pf` file to see where changes arose.


Synopsis:  

    The functions here are intended to ease running a collection
    of commits on the same parameter file so that one can locate
    where a change in Python occured.


Command line usage (if any):

    None at present. These are intend to be run as part 
    as a Jupyter script or from a python window

Description:  

Primary routines:

    log2table  -
        create an ascii table that summarizes the log and contains
        entries for compiling and running python many times
    compile_many -
        create at lot of Python executables in a local directory
    run_many  -
        runs a set of python executables
    plot_many  -
        plot consecutive sets of the outputs


Notes:

    To use this, create a directory where one wants to run

    Create a text file by running::

        git log > commit.txt (or whatever

    Edit this file so it includes only the commit range in which you are interested.

    Run the function log2table to create a file which contains a summary of the logfile.
    This file, hereafter the masterfile, has additional columnns needed for 
    the remainder of the effort

    Run the function compile_many to compile every nth commit.  The various versions of 
    python will be called py_x001, etc, where 001 corresponds to the nubmer of the commit
    going from most recent to oldest.

    Put the .pf file in your local directory.  The runs are going to occur in a subdirectories
    name r_py_x001 etc.  

    If the run requires special files, such as seds that are not part of the normal python 
    distribution, one needs to adjust the parameter files to work in the subdirectories.
    (The above sentence is worth repeating, because the system can easily hang up if
    this is not correct.  python will exit but mpirun will not catch this problem, and
    so the subprocess routine will not exist until it times out.)

    Run the function run_many to run all of the model that you have compiled.  

    Run the function plot_many to examine what the various spectra look like

                                       
History:

    210220 ksl Coding begun

'''

import sys
import os
import subprocess
import shutil


from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt



from scipy.signal.windows import boxcar
from scipy.signal import convolve


# Define the directory from which the
# program is called.
home=os.getcwd()


FAILURE=True
SUCCESS=False


def read_table(filename='foo.txt',format=''):
    '''
    Read a file using astropy.io.ascii and
    return this

    Description:

    Notes:

    This version handles changing string
    lengths better than a simply reading
    the table with ascii.read

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
        return []

    return data


def log2table(logfile='commits.txt',masterfile='retro_master.txt'):
    '''
    Create a ascii table from a file that contains 
    a normal git log (or the portion of a log that 
    has been edited to be in the time period of interest

    This is normall produced by 

        git log > commits.txt

    One should ensure that the edited file contains complete
    records.
    '''
    x=open(logfile)
    line=x.readlines()
    commit=[]
    author=[]
    date=[]
    info=[]
    xinfo=''
    i=0
    while i < len(line):
        q=line[i].strip()
        word=q.split()
        if len(q)==0:
            pass
        elif word[0]=='commit':
            if len(xinfo):
                if len(xinfo)>100:
                    info.append(xinfo[0:100])
                else:
                    info.append(xinfo)
            xinfo=''
            commit.append(word[1])
            xnext='author'
        elif word[0]=='Author:' and xnext=='author':
            author.append(word[1])
            xnext='date'
        elif word[0]=='Date:' and xnext=='date':
            xdate='%s %s %s %s' %(word[2],word[3],word[4],word[5])
            date.append(xdate)
            xnext='info'
        elif xnext=='info':
            xinfo='%s %s' % (xinfo,q)
        i+=1
    if len(xinfo)>100:
        info.append(xinfo[0:100])
    else:
        info.append(xinfo)
    
    print(len(commit),len(author),len(date),len(info))
    
    xnum=np.linspace(0,len(commit)-1,len(commit),dtype=int)
    
    result=Table([xnum,commit,author,date,info], names=['Number','Commit','Authur','Date','Info'])

    result['Compiled']='No'
    result['Changed']='Unknown'
    result['Ran']='No'
    
    result.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)

    return SUCCESS
        
SIROCCO_SOURCE_DIRECTORY=''

def get_python_source_directory():
    global SIROCCO_SOURCE_DIRECTORY

    result=subprocess.run(['which','py'],capture_output=True,text=True)
    location=str(result.stdout)
    x=location.strip()
    xx=x.replace('bin/py','source/')

    SIROCCO_SOURCE_DIRECTORY=xx
    print('Set Python source directory to :',SIROCCO_SOURCE_DIRECTORY)

    return  SIROCCO_SOURCE_DIRECTORY 


def compile_one(commit='a77ec180c017244cb56f41b50178daf81541748a',number=35,print_output=False):
    '''
    Compile a single commit of Python and move the executable to a local directory.
    The name of the executable is determined  by the number
    '''

    # cwd=os.getcwd()

    if SIROCCO_SOURCE_DIRECTORY=='':
        get_python_source_directory()
    os.chdir(SIROCCO_SOURCE_DIRECTORY)

    result=subprocess.run(['git', 'checkout', commit],capture_output=True,text=True)
    if result.returncode:
        print ('Error: Could not checkout new branch, commit all changes before rerunning')
        print("stderr:", result.stderr)
        return FAILURE
    if print_output:
        print("stdout:", result.stdout)

    result=subprocess.run(['make', 'clean'],capture_output=True,text=True)

    version='VERSION=_x%03d' % number
    
    result=subprocess.run(['make', version, 'python'],capture_output=True,text=True)
    if result.returncode:
        print('Error: Could not compile commit %s',commit)
        print("stderr:", result.stderr)
        return FAILURE

    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    

    if os.path.isfile(home+'/py_x%03d' % number):
        os.remove(home+'/py_x%03d' % number)
    # if os.path.isfile(cwd+'/py_x%03d' % number):
    #     os.remove(cwd+'/py_x%03d' % number)
    
    success=SUCCESS
    try:
        shutil.move('../bin/py_x%03d' % number,home)
        # shutil.move('../bin/py_x%03d' % number,cwd)
    except:
        print('Could not find ../bin/py_x%03d for commit %s' %(number,commit))
        success=FAILURE

    # return the repository to a known state
    result=subprocess.run(['git', 'checkout', 'dev'],capture_output=True,text=True)

    # os.chdir(cwd)
    os.chdir(home)
    return success



def check_inputs(command='py',parameter_file='cv.pf',print_output=True):
    '''
    Check that a parameter file has the appropriate parameters to run with 
    a version of python
    '''
    try:
        result=subprocess.run([command, '-i', parameter_file],capture_output=True,text=True,timeout=5)
    except:
        print('Executable %s appears missing' % command)
        return FAILURE
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    if result.stdout.count('Exiting'):
        print('Parameter File problem with %s %s' % (command,parameter_file))
        return FAILURE
    return SUCCESS



def check_completion(parameter_file='cv.pf'):
    '''
    Check the .sig file to see if a run has completed
    '''
    root=parameter_file.replace('.pf','')
    xsig=root+'.sig'
    f=open(xsig)
    lines=f.readlines()
    print(lines[-1].strip())
    if lines[-1].count('COMPLETE'):
        return SUCCESS
    return  FAILURE

def run_one(command='py',parameter_file='cv.pf',np=8, extra_inputs='',print_output=False,timeout=3000):
    '''
    Run a single model, with a timeout, in the current working directory

    Here the timeout does not use the timeout for Python itself but is on the
    python subprocess call.
    '''
    
    result=subprocess.run(['Setup_Py_Dir'],capture_output=True,text=True)

    if check_inputs(command,parameter_file,print_output)==FAILURE:
        print('Inputs failed')
        return FAILURE
    
    try:
        result=subprocess.run(['mpirun', '-np',str(np),command, extra_inputs,parameter_file],capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        print('Failed for timeout')
        return FAILURE

    if result.returncode:
        print ('Error: could not run command ',command)
        print("stderr:", result.stderr)
        return FAILURE


    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    if result.stdout.count('Exiting'):
        print('Parameter File problem with %s %s' % (command,parameter_file))
        return FAILURE
    return check_completion(parameter_file)


def check4differences(masterfile='retro_master.txt',print_output=True):
    '''
    Check for changes between compiled versions
    '''
    os.chdir(home)
    xsum=read_table(masterfile)

    n_compiled=0
    n_changed=0
    n_unchanged=0
    n_error=0


    i=0
    old=-1
    while i<len(xsum):
        if xsum['Compiled'][i]=='Yes':
            n_compiled+=1
            if old== -1:
                xsum['Changed'][i]='Yes'
                n_changed+=1
            else:
                try:
                    if os.path.getsize('py_x%03d' % xsum['Number'][old])==os.path.getsize('py_x%03d' % xsum['Number'][i]):
                        xsum['Changed'][i]='No'
                        if print_output:
                            print('No Change for py_x%03d vs py_x%03d' % (xsum['Number'][i],xsum['Number'][old]) )
                        n_unchanged+=1
                    else:
                        if print_output:
                            print('   Change for py_x%03d vs py_x%03d' % (xsum['Number'][i],xsum['Number'][old]) )
                        xsum['Changed'][i]='Yes'
                        n_changed+=1
                except:
                    xsum['Changed'][i]='Unknown'
                    n_error+=1
            old=i
        i+=1
    xsum.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)

    print('Of %d compiled, %d changed %d were not changed and there were %d problems' % 
            (n_compiled, n_changed, n_unchanged, n_error))
    return
            


def compile_many(masterfile='retro_master.txt',increment=10, nmin=0,nmax=0,rerun=False):
    '''
    This routine compiles a number of versions of python and checks
    to see if the code has changed between the time it has changed.
    It adds a flag about this to the table

    '''
    xsum=read_table(masterfile)
    i=0
    old=-1
    
    if  nmin == nmax:
        imin=0
        imax=len(xsum)
    else:
        if nmax<len(xsum):
            imax=nmax
        else:
            imax=len(xsum)
        imin=0
        if nmin>0 and nmin<nmax:
            imin=nmin

    i=imin
    while i<=imax:
        # print(xsum['Commit'][i],xsum['Number'][i])
        if rerun == False and  xsum['Compiled'][i] == 'Yes':
            pass
        else:
            success=compile_one(xsum['Commit'][i],xsum['Number'][i],print_output=False)
            print('Tried to compile commit %s or no %d, which returned %d  where SUCCESS== %d' % (xsum['Commit'][i],xsum['Number'][i],success,SUCCESS))
            if success == SUCCESS:
                # print('Yes it was successful')
                xsum['Compiled'][i]='Yes'
            else:
                xsum['Compiled'][i]='Failed'
                xsum['Changed'][i]='Failed'
        # print(xsum[i])
        xsum.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)
        i+=increment

    check4differences(masterfile,print_output=False)
    return 


def run_many(master='retro_master.txt',parameter_file='cv.pf',extra_inputs='',np=8,nmin=0,nmax=0,rerun=False,run_changed_only=False):
    '''
    Run many versions of python on a single parameter file.
    
    The routines are run in a subdirectory that is defined by the
    version number, so the parameter files must be able to run from there.
    In particular if any special files are called then the path to the
    files must be defined correctly, and any set of model files must
    be defined so that they will be found from the subdirectory as well.

    nmin and nmax indicate the range of models to run

    rerun indicates that you want to rerun models that have been rerun previously. (If
    you want to run a differenct parameter file, without resetting the master file up,
    you should use this option.)

    if changed_only is True, then you will only run models whose binary file length
    has changed from the previous commit.  Changes of large scope generally change the
    size of the binary file, but for small changes the file lennght can stat the same,
    one might want to use this option as a way to coarsely isolate a problem in one of the
    commits, but one should be aware that small canges to not necessarily change the
    binary file length.

    Note that if one is sure there has been a large change between commits and the
    file lenght does not change, be sure to chack that the commit number is what 
    one expects between the two runs.

    '''

    start_dir=os.getcwd()
    x=read_table(master)
    # Check that the file contains the column Ran, and if not add it
    
    if  nmin == nmax:
        imin=0
        imax=len(x)
    else:
        if nmax<len(x):
            imax=nmax
        else:
            imax=len(x)
        imin=0
        if nmin>0 and nmin<nmax:
            imin=nmin

    print('Running for versions of Python from %d to %d' % (nmin,nmax))
    if run_changed_only:
        print('Only running models for which the binary file length changes between models')
        print('Warning: Not all code changes result in changes to the binary file length')
    if rerun:
        print('Rerunning all models, not just those that had not been run previously')

    i=imin
    while i<imax:
        one=x[i]

  

        if rerun:
            one['Ran']='No'

        if one['Compiled']!='Yes':
            pass
        elif one['Ran']=='Yes':
            pass
        elif one['Changed']=='Yes' or run_changed_only == False:
            xcommand='py_x%03d' % one['Number']
            print('Starting %s' % xcommand)
            if os.path.isdir('r_'+xcommand)==False:
                os.mkdir('r_'+xcommand)
            shutil.copy(parameter_file,'r_'+xcommand)
            os.chdir('r_'+xcommand)

            success=run_one(command='../'+xcommand,parameter_file=parameter_file,
                    extra_inputs=extra_inputs,np=np, print_output=False,timeout=8000)
            if success==SUCCESS:
                one['Ran']='Yes'
                print('%s completed successfully' % (xcommand))
            else:
                one['Ran']='Failed'
                print('%s failed to complete' % (xcommand))
            os.chdir(start_dir)
            x.write(master,format='ascii.fixed_width_two_line',overwrite=True)
        i+=1
    return 



def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux

    This is just a little utility to 
    allow me  to boxcar smoothe thins
    without having to call the appropirate
    routines from scipy
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)



def plot_two(dir1='r_py_x001',dir2='r_py_x026',file='agn.spec',colname='A01P0.50',smooth=1,xlim=[0,0],ylim=[0,0]):
    '''
    Create plots that compare two Python runs
    '''
    
    try:
        one=ascii.read('%s/%s' %(dir1,file))
    except:
        print('Could not read %s/%s' %(dir1,file))
        return 1
    try:
        two=ascii.read('%s/%s' %(dir2,file))
    except:
        print('Could not read %s/%s' %(dir2,file))   
        return 2
    
    got_col=False
    for col in one.colnames:
        if colname==col:
            got_col=True
    if got_col==False:
        print('File does not have col %s' % colname)
        print('The possible column names are :\n',one.colnames)
        return 1

    root1=dir1.replace('r_py_','')
    root2=dir2.replace('r_py_','')

    if smooth>1:
        y1=xsmooth(one[colname],smooth)
        y2=xsmooth(two[colname],smooth)
    else:
        y1=one[colname]
        y2=two[colname]


    
    plt.figure(1)
    plt.clf()
    plt.plot(one['Lambda'],y1,label=root1)
    plt.plot(two['Lambda'],y2,label=root2)
    if xlim[0]!=xlim[1]:
        plt.xlim(xlim[0],xlim[1])
    if ylim[0]!=ylim[1]:
        plt.ylim(ylim[0],ylim[1])
    plt.legend()
    
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.tight_layout()
    plt.savefig('%s_%s.png' % (root1,root2))
    return 0




def plot_many(masterfile='retro_master.txt',file='agn.spec',colname='A01P0.50',smooth=11,xlim=[0,0],ylim=[0,0]):
    '''
    Plot spectra from a set of models all of which should in various sub-dirctories.

    Here:

    file is the name of the file one of whose columns is to be plotted
    colname is the column to plot
    smooth 
    '''
    x=ascii.read(masterfile)
    xx=x[x['Ran']=='Yes']
    print(len(xx))
    i=1
    while i<len(xx):
        dir1='r_py_x%03d' % (xx['Number'][i-1])
        dir2='r_py_x%03d' % (xx['Number'][i])
        plot_two(dir1,dir2,file,colname,smooth,xlim,ylim)
        i+=1
    return
        

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: retro.py filename')
