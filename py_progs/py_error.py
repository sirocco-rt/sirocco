#!/usr/bin/env python 

'''
Post-run parser for error logs from parallel simulation runs.
    
This is the post-processing code to deal with the error logs in 
SIROCCO's parallel mode. It basically sums the number of errors
in the various diagnostic files

usage and arguments:

    python py_error.py [-f folder] root
    root is the root name of parameter file
    -f folder is the folder the diagfiles are in, defaults to diag_root

returns:

    Error summary printed to screen, and an astropy table root_error_sum.txt. 

History:

    1307    JM  Coding began -- initial tests conducted successfully
    1802    ksl Updated to be python3 compatable, to write the results to an astropy table, and to be callable from another routine

'''

    
## import modules
import sys, subprocess
from astropy.table import Table
from glob import glob




def doit(root, diagfolder_name=''):
    '''
    Given a root name for a parameter file and a location
    for the diag_folder then summarizes the errors for a
    python run and write the outputs to a file
    '''

    if diagfolder_name=='':
        diagfolder_name = 'diag_' + root


    string ='%s/%s_*.diag' % (diagfolder_name,root)
    files=glob(string)

    if len(files)==0:
        print('Error: No diagnostic files found with %s' % string)
        return(0)

    # number of processors and therefore diag files provided as argument
    nfiles = len(files)
    nprocessors = nfiles



    # create some arrays to store the error logs and increment the counters
    # note that we also count the number of threads an error occurs in
    error_log = []
    error_count = []
    thread_count = []



    #everything now initialised, can collate errors

    for i in range(nfiles):

        filepath = files [i]

        # grep for error summary in the diag file
        command = "grep -A 100 'Error summary: End of program' " + filepath +"&"
        output = subprocess.check_output ( command , shell=True )

        # split the grep output up into lines
        output = output.decode()
        data = output.split("\n")
        if len(data) <=1:
            print('No errors for thread %3d, which is unusual' % i)
        # now cycle through the lines in output
        for j in range(len(data)):
            No_log = True    # variable which is True if this is a new error

            # check we have a non blank entry to avoid indexing error
            if len(data[j])>0:
                error = data[j].split()

                # this next line is to check we have an actual error summary and not a string line
                if error[0].isdigit():

                    # now get the err log string, and its count
                    error = data[j].split()
                    
                    count = int(error[0])
                    log_string = data[j][13:]

                    # search for a match in the already fpund errors
                    for k in range(len(error_log)):
                        if log_string == error_log[k]:

                            No_log = False    # we've found a match, increment the counts
                            error_count[k] += count
                            thread_count[k] += 1
                    
                    # this is a new error, so append to the log and count arrays
                    if No_log:
                        error_log.append ( log_string )
                        error_count.append ( count )    
                        thread_count.append (1)


    n_logs=len(error_log)

    # print the output to screen
    print("Error summary: "+root)
    print("Collated errors for " + str(nprocessors) + " processors")


    try:
        x=Table([error_count,thread_count,error_log],names=['ErrorCount','ThreadCount','Error'])
        x.sort('ErrorCount')
        x.reverse()
        x.write('%s_error_sum.txt' % root, format='ascii.fixed_width_two_line',overwrite=True) 
        x.sort(keys=['ErrorCount'])
        x.reverse()  # put in order of the most errors first
    except:
        print("Problem creating error counts")
        if n_logs==0:
           print('No errors were parsed from the diag files, suggesting Sirocco failed to run to completion')
        return

    print("Recurrences --  number of threads with error -- Description")
    for one in x:
        print("\t%5d -- %5d -- %s" % (one['ErrorCount'],one['ThreadCount'],one['Error']))
    # all done.
    return
                    

if __name__ == "__main__":        # allows one to run from command line without running automatically with write_docs.py
    # read root pf file from command line
    if len(sys.argv) > 1:
        root = sys.argv[-1]
    else:
        print('Error: no arguments provided. Exiting.')
        print(__doc__)
        exit (0)



    # set defaults
    nprocessors = 1
    diagfolder_name = ''  # Use the default

    # check for args
    for i in range(len(sys.argv)):
        if sys.argv[i]=="-f":
            diagfolder_name = sys.argv[i+1]
        elif sys.argv[i]=='-h':
            print(__doc__)
            exit (0)



    doit(root, diagfolder_name)


    

