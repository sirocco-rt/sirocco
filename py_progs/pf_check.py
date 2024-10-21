#!/usr/bin/env python 

'''
Checks if `.pf` files in a directory need modification.

Synopsis:  

    Check all of the pf files in a directory to see if
    they require modification, alternatively run pieces of travis
    to make sure the parameter files will run.


Command line usage (if any):

    usage::

        pf_check.py directory_name

Description:  

    The program checks all the .pf files in a directory by
    running the file with Python with the -i switch turned 
    on so to check whether the pf file will run, querying
    the user for any missing inputs.  

    In order that ther original directory not be affected,
    it first copies all of the .pf files to a directory 
    with the name::

        pf_check_directory_date

    The full path names to directories in the `$SIROCCO/examples`
    do not need to be given, so for example::

        pf_check.py beta

    will test the parameter files in `$SIROCCO/examples/beta`.

    After running Python in -i mode for each of the files, 
    the routine diffs the .out.pf file with the original.pf
    file to see if changes should be made to the .pf file

    A special procedure is used for the travis directory.  

    There the .travis.yaml file is parsed to see what options
    have been selected, e.g whether -d option has been included
    for that file.  For the travis directory a -i option is added
    if it is missing, but then the travis commands are run.

    It is up to the user to edit the .pf files and put them back
    in the original directory.

    To help with this though, a file DoCommands is place in the
    output directory.  One can source this file to re-execute
    the python runs 

Primary routines:

    doit
        the main routine which runs everything
    travis
        the special routine for parsing the .travis.yaml file to
        get the commands to run
    steer
        interpret the command line

Notes:

    At present, there is no way general way (except for the travis directory)
    to handle .pf files which are supposed to be run with the -d switch.  
    What could be done is to check for diagnostic commands and then add the
    -d switch.
                                       
History:

    190304 ksl
        Coding begun

'''


import sys
from astropy.io import ascii
import numpy
from glob import glob
import time
import subprocess
import os
import shutil
import yaml


def travis():
    '''
    Return a special set of commands when dealing with the travis directory
    '''

    SIROCCO=os.environ['SIROCCO']

    x=open('%s/.travis.yml' % (SIROCCO))
    z=yaml.load(x)
    print(z['script'])
    commands=[]
    for one in z['script']:
        print(one)
        if one[0:3]=='py ':
            commands.append(one)
    print(commands)

    # Wow have to deal with the fact that occassionally we are not running in 
    # interactive mode

    i=0
    while i<len(commands):
        one=commands[i]
        if one.count('-i ')==0:
            commands[i]=one.replace('py ','py -i ')
        i+=1
    
    print(commands)
    return commands




def doit(directory):
    '''
    '''


    xdir=''
    if os.path.isdir(directory):
        xdir=directory
    else:
        # Get the SIROCCO environment variable
        SIROCCO=os.environ['SIROCCO']
        ydir='%s/examples/%s' % (SIROCCO,directory)
        if os.path.isdir(ydir):
            xdir=ydir
        else:
            print('Error: Could not find %s/n' % directory)
            return
    # So at this point we know the directory we want to process

    date=time.strftime("%y%m%d", time.gmtime())

    word=directory.split('/')
    out_dir='pf_check_%s_%s' % (word[len(word)-1],date)

    if os.path.exists(out_dir)==False:
        os.mkdir(out_dir)

    # Now copy the pf files to this directory, avoiding the out.pf files if they exist

    pf_files=glob(xdir+'/*.pf')
    txt_files=glob(xdir+'/*.txt')

    # Screen out the .out.pf files
    select=[]
    for one in pf_files:
        if one.count('.out.pf')==0:
            select.append(one)
    pf_files=select


    if len(pf_files)==0:
        print ('No input files found for %s search' % (pf_dir+'/*pf'))
        return
    else:
        print ('Executing models for the following models:')
        for one in pf_files:
            print(one)



    # Get the rootnames and copy the .pf files
    root_names=[]
    for one in pf_files:
        print('Copy',one)
        shutil.copy(one,out_dir)
        words=one.split('/')
        pf=(words[len(words)-1])
        root_name=pf.replace('.pf','')
        root_names.append(root_name)

    commands=[]
    for one in root_names:
        words=one.split('/')
        pf=(words[len(words)-1])
        root_name=pf.replace('.pf','')
        command='py -i %s.pf' % (pf)
        commands.append(command)

    if directory=='travis':
        commands=travis()

    # get any text files if any

    for one in txt_files:
        shutil.copy(one,out_dir)
    			
    # Switch to the work directory
    cwd=os.getcwd()
    os.chdir(out_dir)


    print('The commands that will be executed will be:')
    f=open('DoCommands','w')
    g=open('Cycle_files','w')
    for one in commands:
        print(one)
        f.write('%s >foo.txt\n' % one)

    proc=subprocess.Popen('Setup_Py_Dir',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)


    for one in commands:
        print('\nRunning: %s\n' % one)
        proc=subprocess.call(one,shell=True,stdout=subprocess.PIPE)



    # outputfile='test'
    # run_cmds(commands,root_names,outputfile)


    # Now check what files have changed

    print('Summary of results for directory %s' % directory)
    print('Look in the output direcory %s for more information, and a command file to rerun the tests' % out_dir) 

    f.write("echo 'Checking for Problems'\n")

    for one in root_names:
        cmd='diff -b  %s.pf %s.out.pf' % (one,one)  #ignore differences in white space length
        f.write("echo 'Check %s.pf'\n" % one)
        f.write('%s\n' % cmd)
        g.write('mv %s.out.pf %s.pf\n' % (one,one))
        proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()
        # print('stdout: ',stdout)
        # print('stderr: ',stderr)
        if len(stdout):
            string="%30s.pf has changed - important to investigate" % one
            print(string)
        else:
            string="%30s.pf unchanged   - test passed"  % one
            print(string)


    f.close()

    os.chdir(cwd)
    return out_dir
	

    
def steer(argv):
    '''
    Parse the command line and set everything up
    '''

    directory=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif directory=='':
            directory=argv[i]
        else:
            print('Error: Uninterpretable inputs')
            return
        i+=1

    doit(directory)





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
