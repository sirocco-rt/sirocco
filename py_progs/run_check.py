#!/usr/bin/env python 
'''
Synopsis:  
    Sumarize a model run with python, ultimately generating
    an html file with various plots, etc.


Command line usage (if any):
    run_check.py root1 [root2 ...]

    run_check.py root1.pf [root2.pf ...]

    run_check -all

    run_check -h
            

Description:  
    This routine performs basic checks on one or more python runs
    and creates an html file for each that is intended to provide
    a quick summary of a run.  

    The user can enter the runs to be tested from the command line,
    either in the form of a set of root names or .pf names. Wildcarding,
    e.g ``*.pf`` can be used.  ``*.out.pf`` files will be ignored.  

    Alternatively to process all of the files in a directory, one can use
    the switch -all (which supercedes anything else).  

    In all cases the routine checks to see if the appropriate wind_save file
    exists before attempting to run.

    -h delivers this documentation

Primary routines:
    doit - processes a single file
    steer  - processes the command line and calls doit for each file.

Notes:
'''

import sys
from glob import glob
import os
from astropy.io import ascii
import numpy
import subprocess
import matplotlib.pyplot as pylab
import xhtml
import plot_wind
import plot_wind_1d
import plot_spec
import plot_tot


def read_diag(root):
    '''
    Get convergence and possibly other information from the diag file
    '''

    filename='diag_%s/%s_0.diag' % (root,root)

    command="grep 'Check_convergence' %s" % filename

    proc=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

    if len(stderr):
        print("Error - grep")
        print(stderr.decode())
        return [],[]
    else:
        # print(stdout.decode())
        x=stdout.decode()
        # print(len(x))
        lines=x.split('\n')
        # print(len(lines))
        converging=[]
        converged=[]
        t_r=[]
        t_e=[]
        hc=[]
        for line in lines:
            if line.count('converged'):
                xline=line.replace('(',' ')
                xline=xline.replace(')',' ')
                words=xline.split()
                converged.append(eval(words[2]))
                converging.append(eval(words[6]))
                ncells=eval(words[9])
            elif line.count('hc(real'):
                words=line.split()
                t_r.append(eval(words[2]))
                t_e.append(eval(words[4]))
                hc.append(eval(words[8]))

        if len(t_r)>0:
            t_r=numpy.array(t_r)
            t_e=numpy.array(t_e)
            hc=numpy.array(hc)


            t_r=t_r/ncells
            t_e=t_e/ncells
            hc=hc/ncells
        else:
            t_r=[]
            t_e=[]
            hc=[]
            print('Read diag file, but there is not evidence of ionization cycles')


        return converged,converging,t_r,t_e,hc


def windsave2table(root):
    '''
    Run windsave2table
    '''

    command='windsave2table %s' % root

    proc=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    # print(stdout.decode())
    # print('Error')
    # print(stderr.decode())
    # print('return code')
    # print (proc.returncode)
    if proc.returncode:
        print('Error: running windsave2table',proc.returncode)
        return True
    elif len(stderr):
        print('Error: running windsave2table')
        print(stderr.decode())
        return True 
    else:
        return False

def py_error(root):
    '''
    Run py_error.py and capture the output to the 
    screen

    Note:
        py_error could be refactored so that it did not
        need to be run from the command line, but this
        is the simplest way to capture the outputs at
        present
    '''

    command='py_error.py %s' % root

    proc=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    if len(stderr):
        print("Error - py_error.py")
        print(stderr.decode())
        return []
    else:
        # print(stdout.decode())
        x=stdout.decode()
        # print(len(x))
        lines=x.split('\n')
        return(lines)


def plot_converged(root,converged,converging,t_r,t_e,hc):
    '''
    Make a plot of the convergence statistics in the diag directroy
    '''


    pylab.figure(1,(6,6))
    pylab.clf()
    pylab.plot(t_r,'--',label='t_r')
    pylab.plot(t_e,'--',label='t_e')
    pylab.plot(hc,'--',label='heating:cooling')
    pylab.plot(converged,lw=2,label='Converged')
    pylab.plot(converging,lw=2,label='Converging')
    pylab.legend(loc='best')
    pylab.xlabel('Ionization Cycle')
    pylab.ylabel('Fraction of Cells in Wind')
    pylab.savefig('diag_%s/convergence.png' % root)

    return

def check_completion(root):
    '''
    Verify that the run actually completed and
    provide information about the timeing, from
    the .sig file
    '''

    try:
        x=open(root+'.sig')
        lines=x.readlines()
    except:
        print('Error: %s.sig does not appear to exist' % root)
        return []

    ion_string=[]
    spec_string=[]
    restart=False
    tot=0

    i=0
    for line in lines:
        if line.count('Finished'):
            if line.count('ionization'):
                ion_string=line
            elif line.count('spectrum'):
                spec_string=line
        elif line.count('RESTART'):
            restart=True
            ion_string=[]
            spec_string=[]
            word=lines[i-1].split()
            print(word)
            tot+=eval(word[5])
        i+=1

    complete_string=lines[len(lines)-1]

    complete_message=[]

    ion_time=0 # To handle the case where there were no ionization cycles

    if restart:
        message='%s was restarted. The earlier run(s) took %.1f s. ' % (root,tot)
        complete_message.append(message)
        message='Here we report times and progress since the (last) restart.'
        complete_message.append(message)

    if line.count('COMPLETE'):
        word=complete_string.split()
        message='%s ran to completion in %s s.' % (root,word[5])
        complete_message.append(message)
        try:
            word=ion_string.split()
            message='%s ionization cycles were completed in %s s.' % (word[8],word[5])
            ion_time=eval(word[5])
        except:
            message='There were no ionization cycles for this run'
        complete_message.append(message)

        try:
            word=spec_string.split()
            spec_time=eval(word[5])
            message='%s   spectrum cycles were completed in %s s' % (word[8],spec_time-ion_time)
        except:
            message='There were no spectrum cycles for this run'
        complete_message.append(message)

    else:
        message='WARNING: RUN %s HAS NOT COMPLETED SUCCESSFULLY' % root
        complete_message.append(message)
        word=complete_string.split()
        try:
            message='If not running, it stopped after about %s s in  %s of %s %s cycles' % (word[5],word[8],word[10],word[11])
        except IndexError: # add catch for when the above can't be read 
            message='could not read cycle information'
            pass 
        complete_message.append(message)

    return complete_message

convergence_message='''
The plot shows the fraction of cells which satisfy various convergence criteria, comparing the radiation tempererature, t_r, 
the electron temperature, t_e between cycles, or the balance between heating and cooling in a given cycle.  If the fractional
variation in, for example, the electron temperature has changed by less than 5%, then the electron temperature criterion has been passed.
If a cell passes all three tests, it is said to have "converged".  The plot also indicates the number of cells for which the electron 
temperature is evolving to a higher or lower value, rather than oscillating up and down in electron temperature. Such cells are said to be
"converging".
'''

convergence2='''
The plot below indicates which cells have converged.  The values given for each cell indicate  the number of convergence tests 
which have failed, so a value of 0 indicates that a cell is "converged".
'''

error_description='''
Python accumulates a fairly large number of errors and warnings.  Most are benign, especially if they only occur a few
times, but one should beware if an error message occurs many times, or if messages that look unusual start to appear.
A summary of the errors for this run of the program is shown below.  More information about where the errors occured
can be found in the diag files directory.  This summary presents the total number of times an error message of a particularly 
type was generated; many times the same error message will will occur in each of the threads.  
'''

tot_plot_description='''
This plot is based on the .spec_tot file.  It shows the photons that were created and emitted in the last ionization cycle.  
The photons created by a central object or disk (the external sources) are plotted separately from the wind (the internal source).
For a simple atom model, the wind photons are generated based on physical conditions in  the wind at the beginning of the cycle.
Ror a macro-atom models, the wind photons are generated as photons pass through the wind and interact with macro-atoms.  The observed
flux from photons that escape to infinity and a subset of this, the emitted flux of photons arising in the wind are also shown.
Finally, the flux of the photons that hit the disk or star is also shown (as they would be seen at a distance of 100 pc).
'''

spec_plot_description='''
This plot is from the .spec file, and shows the expected spectra (at 100 pc) as a function of the various inclination angles requested. Note
that the y-axis scale varies for different inclination angles.
'''
    
def make_html(root,converge_plot,te_plot,tr_plot,spec_tot_plot,spec_plot,nspectra=3,complete_message=['test'],errors=['test','test2']):
    '''
    Make an html file that collates all the results
    '''

    string=xhtml.begin('%s: How well did the Python run go?' % root)

    string+=xhtml.paragraph('Provide an overview of whether the run of %s has succeeded' % root)

    string+=xhtml.add_list(complete_message)

    string+=xhtml.hline()
    string+=xhtml.h2('Did the run converge?')

    if os.path.isfile('./diag_%s/convergence.png' % root):
        string+=xhtml.image('file:./diag_%s/convergence.png' % root)
    else:
        string+=xhtml.paragraph('There is no convergence plot. OK if no ionization cycles')

    string+=xhtml.paragraph(convergence_message)

    string+=xhtml.paragraph(convergence2)
    # print('Test ',converge_plot)
    string+=xhtml.image('file:%s' % (converge_plot))
    string+=xhtml.hline()

    string+=xhtml.h2('What do the temperatures look like?')
    string+=xhtml.image('file:%s' % (te_plot))
    string+=xhtml.image('file:%s' % (tr_plot))
    # string+=xhtml.image('file:./diag_%s/%s_t_e.png' % (root,root))
    # string+=xhtml.image('file:./diag_%s/%s_t_r.png' % (root,root))
    string+=xhtml.hline()
    string+=xhtml.h2('What do the total spectra look like (somewhat smoothed)?')

    if os.path.isfile(spec_tot_plot):
        string+=xhtml.image('file:%s' % (spec_tot_plot),width=800)
        string+=xhtml.paragraph(tot_plot_description)
    else:
        string+=xhtml.paragraph('There is no total spectrum plot. OK if no ionization cycles')

    string+=xhtml.hline()
    string+=xhtml.h2('What do the final spectra look like (somewhat smoothed)?')
    if spec_plot != 'None':
        if nspectra==1:
            string+=xhtml.image('file:%s' % (spec_plot),width=900,height=600)
        else:
            string+=xhtml.image('file:%s' % (spec_plot),width=900,height=nspectra*300)
        string+=xhtml.paragraph(spec_plot_description)
    else:
        string+=xhtml.paragraph('There is no plot of a detailed spectrum, probably because detailed spectra were not created')
    string+=xhtml.hline()

    string+=xhtml.h2('Errors and Warnings')

    string+=xhtml.paragraph(error_description)

    string+=xhtml.preformat(errors)

    string+=xhtml.hline()

    # Now add the parameter file so we can see what we have

    string+=xhtml.h2('The parameter file used:')

    x=open(root+'.pf')
    lines=x.readlines()

    string+=xhtml.add_list(lines)
    string+=xhtml.hline()



    string+=xhtml.end()

    # print(string)

    g=open(root+'.html','w')
    g.write(string)

def how_many_dimensions(filename):
    '''
    Check whether a windsave file is one or two dimenaions
    '''

    x=ascii.read(filename)
    for one in x.colnames:
        if one=='j':
            return 2
    
    return 1






def doit(root='ixvel',outputfile='out.txt'):
    '''
    Create a summary of a Python run, which will enough information that one can assess
    whether the run was successful

    Description:

    Notes:

    History:


    '''

    print('\nEvaluating %s\n' % root)
    if windsave2table(root):
        print('Exiting becase windsave2table failed')
        return

    complete_message=check_completion(root)
    if len(complete_message)==0:
        print('Exiting because could not parse %s.sig file' % root)

    for one in complete_message:
        print(one)

    xdim=how_many_dimensions('%s.master.txt' % root)

    if xdim==2:
        converge_plot=plot_wind.doit('%s.master.txt' % root,'converge',plot_dir='./diag_%s' % root)
        te_plot=plot_wind.doit('%s.master.txt' % root,'t_e',plot_dir='./diag_%s' % root)
        tr_plot=plot_wind.doit('%s.master.txt' % root,'t_r',plot_dir='./diag_%s' % root)
    else:
        converge_plot=plot_wind_1d.doit('%s.master.txt' % root,'converge',plot_dir='./diag_%s' % root)
        te_plot=plot_wind_1d.doit('%s.master.txt' % root,'t_e',plot_dir='./diag_%s' % root)
        tr_plot=plot_wind_1d.doit('%s.master.txt' % root,'t_r',plot_dir='./diag_%s' % root)


    converged,converging,t_r,t_e,hc=read_diag(root)

    if len(converged)>1:
        plot_converged(root,converged,converging,t_r,t_e,hc)
    else:
        print('There were not enough cycles to plot the convergence by cycle')

    plot_tot.doit(root)
    spec_tot_plot=root+'.spec_tot.png'

    try:
        spec_plot,nspectra=plot_spec.do_mosaic(root,wmin=0,wmax=0)
    except:
        print('Could not construct detailed spectrum plot')
        spec_plot='none'
        nspectra=0

    # spec_plot=root+'.png'

    errors=py_error(root)
     

    make_html(root,converge_plot,te_plot,tr_plot,spec_tot_plot,spec_plot,nspectra,complete_message,errors)


    return

def steer(argv):
    '''
    Process the command line
    '''

    xall=False
    files=[]
    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
        elif argv[i]=='-all':
            xall=True
            break
        elif argv[i].count('.out.pf'):
            pass
        elif argv[i].count('.pf'):
            string=argv[i]
            string=string.replace('.pf','')
            files.append(string)
        elif argv[i][0]=='-':
            print('Unknown switch: %s' % argv[i])
            return
        else:
            files.append(argv[i])
        i+=1

    if xall==True:
        files=glob('*.wind_save')
        for one in files:
            doit(one.replace('.wind_save',''))
    else:
        for one in files:
            if os.path.isfile(one+'.wind_save'):
                    doit(one)
            else:
                print('Error: No windsave file for %s' % one)




    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(sys.argv[1])
        steer(sys.argv)
    else:
        print (__doc__)
