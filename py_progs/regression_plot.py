#!/usr/bin/env python 

'''
Create plots which compare spectra (and other properties) of two runs of Python.

This may become a basis for some kind of ipython notebook


Command line usage (if any):

    usage::

        regression_plot.py run1 run2 [model]

    where run1 and run 2 are in directories containing the runs to compare.
    If model is given, then only that model is compared.  
    If model is not given all the models in the dirctory Xcompare

Description:  

Primary routines:

    doit_two  Creates a plot comparing the spectra from one model
    do_all    Runs doit_two for all of the models in the named directories

Notes:
                                       
History:

181127 ksl Coding begun

'''

import os
import sys
from astropy.io import ascii
import numpy as np
import pylab
from glob import glob

from scipy.signal.windows import boxcar
from scipy.signal import convolve
import time

fig_num=1

def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)

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

def plot_spec(spec_name1,spec_name2,out):
    '''
    Compare two detailed output
    spectra
    '''
    global fig_num
    spec1=ascii.read(spec_name1)
    spec2=ascii.read(spec_name2)

    # spec1.info()

    # pylab.figure(2,(8,8))
    pylab.figure(fig_num,(8,8))
    pylab.clf()

    pylab.subplot(311)
    pylab.title(spec_name1)
    flux1=xsmooth(spec1['Created'],21)
    flux2=xsmooth(spec2['Created'],21)
    diff_created=flux2-flux1
    pylab.plot(spec1['Lambda'],flux1,label='Created1')
    pylab.plot(spec2['Lambda'],flux2,label='Created2')
    pylab.legend(loc='best')
    pylab.ylabel('Flux')
    # pylab.xlabel('Wavelength')

    pylab.subplot(312)
    flux1=xsmooth(spec1['Emitted'],21)
    flux2=xsmooth(spec2['Emitted'],21)
    diff_emitted=flux2-flux1
    pylab.plot(spec1['Lambda'],flux1,label='Emitted1')
    pylab.plot(spec2['Lambda'],flux2,label='Emitted2')
    pylab.legend(loc='best')
    pylab.ylabel('Flux')
    # pylab.xlabel('Wavelength')

    pylab.subplot(313)
    pylab.plot(spec1['Lambda'],diff_created,label='Created')
    pylab.plot(spec1['Lambda'],diff_emitted,label='Emitted')
    pylab.legend(loc='best')
    pylab.ylabel('Flux Difference')
    pylab.xlabel('Wavelength')

    fig_num+=1

    return


def plot_tot(spec_name1,spec_name2,out):
    '''
    Compare two detailed output
    spectra
    '''

    global fig_num

    spec1=ascii.read(spec_name1)
    spec2=ascii.read(spec_name2)

    # spec1.info()


    # pylab.figure(1,(8,8))
    pylab.figure(fig_num,(8,8))
    pylab.clf()

    pylab.subplot(311)
    flux1=xsmooth(spec1['Created'],21)
    flux2=xsmooth(spec2['Created'],21)
    diff_created=flux2-flux1
    pylab.loglog(spec1['Freq.'],flux1,label='Created')
    pylab.loglog(spec2['Freq.'],flux2,label='Created')
    y=pylab.ylim()
    pylab.ylim(y[1]/1e5,y[1])


    flux1=xsmooth(spec1['Emitted'],21)
    flux2=xsmooth(spec2['Emitted'],21)
    diff_emitted=flux2-flux1
    pylab.loglog(spec1['Freq.'],flux1,label='Emitted')
    pylab.loglog(spec2['Freq.'],flux2,label='Emitted')
    y=pylab.ylim()
    pylab.ylim(y[1]/1e5,y[1])
    pylab.ylabel('Flux')
    pylab.legend(loc='best')


    pylab.subplot(312)

    flux1=xsmooth(spec1['CenSrc'],21)
    flux2=xsmooth(spec2['CenSrc'],21)
    diff_CenSrc=flux2-flux1
    pylab.loglog(spec1['Freq.'],flux1,label='CenSrc')
    pylab.loglog(spec2['Freq.'],flux2,label='CenSrc')
    y=pylab.ylim()
    # pylab.ylim(y[1]/1e5,y[1])


    flux1=xsmooth(spec1['Disk'],21)
    flux2=xsmooth(spec2['Disk'],21)
    diff_disk=flux2-flux1
    pylab.loglog(spec1['Freq.'],flux1,label='Disk')
    pylab.loglog(spec2['Freq.'],flux2,label='Disk')
    y=pylab.ylim()
    # pylab.ylim(y[1]/1e5,y[1])


    flux1=xsmooth(spec1['Wind'],21)
    flux2=xsmooth(spec2['Wind'],21)
    diff_wind=flux2-flux1
    pylab.loglog(spec1['Freq.'],flux1,label='Wind')
    pylab.loglog(spec2['Freq.'],flux2,label='Wind')
    y=pylab.ylim()
    pylab.ylim(y[1]/1e5,y[1])
    pylab.ylabel('Flux')

    pylab.legend(loc='best')

    pylab.subplot(313)
    pylab.semilogx(spec1['Freq.'],diff_created,label='Created')
    pylab.semilogx(spec1['Freq.'],diff_emitted,label='Emitted')
    pylab.semilogx(spec1['Freq.'],diff_CenSrc,label='CenSrc')
    pylab.semilogx(spec1['Freq.'],diff_disk,label='Disk')
    pylab.semilogx(spec1['Freq.'],diff_wind,label='Wind')
    pylab.xlabel('Freq')
    pylab.ylabel('Flux difference')
    pylab.legend(loc='best')

    fig_num+=1

    return

def doit_two(run1='py82i_181127',run2='py82i_181126',model='cv_kur',outdir=''):
    '''
    In a single plot summarize it all
    '''

    global fig_num

    # print(run1,run2,model)

    if outdir=='':
        outdir='Xcompare'
    if os.path.isdir(outdir)==False:
        print('making % s to store regression plots' % outdir)
        os.mkdir(outdir)

    pylab.figure(fig_num,(8,16))
    # pylab.figure(1,(8,16))
    pylab.clf()

    # Now compare the spectra made during the ionization cycles

    spec_name1='%s/%s.log_spec_tot' % (run1,model)
    spec_name2='%s/%s.log_spec_tot' % (run2,model)

    if os.path.exists(spec_name1) and os.path.exists(spec_name2):

        spec1=ascii.read(spec_name1)
        spec2=ascii.read(spec_name2)
        # print('Found ',spec_name1,spec_name2)




        pylab.subplot(611)
        pylab.title('%s (1) vs %s (2) for %s' % (run1,run2,model))
        flux1=xsmooth(spec1['Created'],21)
        flux2=xsmooth(spec2['Created'],21)

        try:
            wind=xsmooth(spec1['WCreated'],21)
            flux1=flux1+wind
        except:
            pass


        try:
            wind=xsmooth(spec2['WCreated'],21)
            flux2=flux2+wind
        except:
            pass


        diff_created=flux2-flux1
        pylab.loglog(spec1['Freq.'],flux1,label='Created1')
        pylab.loglog(spec2['Freq.'],flux2,label='Created2')

        
#        iboth=True
#        try:
#            flux1=xsmooth(spec1['WCreated'],21)
#            pylab.loglog(spec1['Freq.'],flux1,label='WCreated1')
#        except:
#            iboth=False

#        try:
#            flux2=xsmooth(spec2['WCreated'],21)
#            pylab.loglog(spec2['Freq.'],flux2,label='WCreated2')
#        except:
#            iboth=False

#        if iboth:
#            diff_wcreated=flux2-flux1

        y=pylab.ylim()
        pylab.ylim(y[1]/1e5,y[1])


        flux1=xsmooth(spec1['Emitted'],21)
        flux2=xsmooth(spec2['Emitted'],21)
        diff_emitted=flux2-flux1
        pylab.loglog(spec1['Freq.'],flux1,label='Emitted1')
        pylab.loglog(spec2['Freq.'],flux2,label='Emitted2')
        y=pylab.ylim()
        pylab.ylim(y[1]/1e5,y[1])
        pylab.ylabel('Flux')
        pylab.legend(loc='best')
        pylab.tight_layout()


        pylab.subplot(612)

        flux1=xsmooth(spec1['CenSrc'],21)
        flux2=xsmooth(spec2['CenSrc'],21)
        diff_CenSrc=flux2-flux1
        pylab.loglog(spec1['Freq.'],flux1,label='CenSrc1')
        pylab.loglog(spec2['Freq.'],flux2,label='CenSrc2')
        y=pylab.ylim()
    # pylab.ylim(y[1]/1e5,y[1])


        flux1=xsmooth(spec1['Disk'],21)
        flux2=xsmooth(spec2['Disk'],21)
        diff_disk=flux2-flux1
        pylab.loglog(spec1['Freq.'],flux1,label='Disk1')
        pylab.loglog(spec2['Freq.'],flux2,label='Disk2')
        y=pylab.ylim()
    # pylab.ylim(y[1]/1e5,y[1])


        flux1=xsmooth(spec1['Wind'],21)
        flux2=xsmooth(spec2['Wind'],21)
        diff_wind=flux2-flux1
        pylab.loglog(spec1['Freq.'],flux1,label='Wind1')
        pylab.loglog(spec2['Freq.'],flux2,label='Wind2')
        y=pylab.ylim()
        pylab.ylim(y[1]/1e5,y[1])
        pylab.ylabel('Flux')

        pylab.legend(loc='best')
        pylab.tight_layout()

        pylab.subplot(613)
        pylab.title('Old (2) - New (1)')
        pylab.semilogx(spec1['Freq.'],diff_created,label='Created')
#        if iboth:
#            pylab.semilogx(spec1['Freq.'],diff_wcreated,label='WCreated')
        pylab.semilogx(spec1['Freq.'],diff_emitted,label='Emitted')
        pylab.semilogx(spec1['Freq.'],diff_CenSrc,label='CenSrc')
        pylab.semilogx(spec1['Freq.'],diff_disk,label='Disk')
        pylab.semilogx(spec1['Freq.'],diff_wind,label='Wind')
        pylab.xlabel('Freq')
        pylab.ylabel('Flux difference')
        pylab.legend(loc='best')
        pylab.tight_layout()

    else:
        print('Error: either %s or %s are missing' % (spec_name1,spec_name2))

    # Now compare the spectra from the spectral cycles
    spec_name1='%s/%s.spec' % (run1,model)
    spec_name2='%s/%s.spec' % (run2,model)


    if os.path.exists(spec_name1) and os.path.exists(spec_name2):
        spec1=ascii.read(spec_name1)
        spec2=ascii.read(spec_name2)

        pylab.subplot(614)
        flux1=xsmooth(spec1['Created'],21)
        try:
            wind=xsmooth(spec1['WCreated'],21)
            flux1=flux1+wind
        except:
            pass

        flux2=xsmooth(spec2['Created'],21)
        try:
            wind=xsmooth(spec2['WCreated'],21)
            flux2=flux2+wind
        except:
            pass




        diff_created=flux2-flux1
        pylab.plot(spec1['Lambda'],flux1,label='Created1')
        pylab.plot(spec2['Lambda'],flux2,label='Created2')

        
#        iboth=True
#        try:
#            flux1=xsmooth(spec1['WCreated'],21)
#            pylab.plot(spec1['Lambda'],flux1,label='WCreated1')
#        except:
#            iboth=False

#        try:
#            flux2=xsmooth(spec2['WCreated'],21)
#            pylab.plot(spec2['Lambda'],flux2,label='WCreated2')
#        except:
#            iboth=False

#        if iboth:
#            diff_wcreated=flux2-flux1

        flux1=xsmooth(spec1['Emitted'],21)
        flux2=xsmooth(spec2['Emitted'],21)
        diff_emitted=flux2-flux1


        xmean=np.average(flux1)
        xmed=np.median(flux1)
        if xmed >0.01*xmean :
            pylab.plot(spec1['Lambda'],flux1,label='Emitted1')
            pylab.plot(spec2['Lambda'],flux2,label='Emitted2')
            ymax=np.max(flux1)
            ymin=np.min(flux1)
            pylab.ylim(0.5*ymin,1.5*ymax)
        else:
            pylab.semilogy(spec1['Lambda'],flux1,label='Emitted1')
            pylab.semilogy(spec2['Lambda'],flux2,label='Emitted2')
            y=pylab.ylim()
            pylab.ylim(y[1]/1e5,y[1])
        


        pylab.legend(loc='best')
        pylab.ylabel('Flux')
        pylab.tight_layout()
    # pylab.xlabel('Wavelength')

        pylab.subplot(615)
        spec_list=spec1.colnames[9:]
        i=len(spec_list)//2
        # print(i,spec_list[i])
        name=spec_list[i]
        # print(spec1.colnames[9:])

        flux1=xsmooth(spec1[name],21)
        flux2=xsmooth(spec2[name],21)
        diff_spec=flux2-flux1

        xmean=np.average(flux1)
        xmed=np.median(flux1)
        if xmed >0.01*xmean :
            pylab.plot(spec1['Lambda'],flux1,label=name)
            pylab.plot(spec2['Lambda'],flux2,label=name)
        else:
            pylab.semilogy(spec1['Lambda'],flux1,label=name)
            pylab.semilogy(spec2['Lambda'],flux2,label=name)
            y=pylab.ylim()
            pylab.ylim(y[1]/1e5,y[1])
        


        pylab.legend(loc='best')
        pylab.ylabel('Flux')
        pylab.tight_layout()
    # pylab.xlabel('Wavelength')

        pylab.subplot(616)
        pylab.title('Old (2) - New (1)')
        pylab.plot(spec1['Lambda'],diff_created,label='Created')
#        if iboth:
#            pylab.semilogx(spec1['Lambda'],diff_wcreated,label='WCreated')
        pylab.plot(spec1['Lambda'],diff_emitted,label='Emitted')
        pylab.plot(spec1['Lambda'],diff_spec,label=name)
        pylab.legend(loc='best')
        pylab.ylabel('Flux Difference')
        pylab.xlabel('Wavelength')
        pylab.tight_layout()




    else:
        print('Error: either %s or %s are missing' % (spec_name1,spec_name2))





    pylab.draw()
    pylab.savefig('%s/%s.png' % (outdir,model))
    print('finished:%s/%s.png'% (outdir,model))

    fig_num+=1









def doit(run1='py82i_181127',run2='py82i_181126',model='cv_kur',outdir=''):
    '''
    Do something magnificent

    Description:

    Notes:

    History:


    '''

    if outdir=='':
        outdir='Xcompare'
    if os.path.isdir(outdir)==False:
        print('making',outdir)
        os.mkdir(outdir)

    # Now compare the spectra made during the ionization cycles

    spec_name1='%s/%s.log_spec_tot' % (run1,model)
    spec_name2='%s/%s.log_spec_tot' % (run2,model)

    if os.path.exists(spec_name1) and os.path.exists(spec_name2):
        plot_tot(spec_name1,spec_name2,outdir)
    else:
        print('Error: either %s or %s are missing' % (spec_name1,spec_name2))

    # Now compare the spectra from the spectral cycles
    spec_name1='%s/%s.spec' % (run1,model)
    spec_name2='%s/%s.spec' % (run2,model)


    if os.path.exists(spec_name1) and os.path.exists(spec_name2):
        plot_spec(spec_name1,spec_name2,outdir)
    else:
        print('Error: either %s or %s are missing' % (spec_name1,spec_name2))


    
    



    return


def do_all(run1='py82i_181127',run2='py82i_181126',outdir=''): 
    '''
    Compare all of the files spectrum files in a two different regression 
    directorys
    '''
    global fig_num

    fig_num=1

    files=glob('%s/*.out.pf' % run1)
    if len(files)==0:
        files=glob('%s/*.spec' % run1)
        ext='.spec'
    else:
        ext='.out.pf'

    for one in files:
        word=one.split('/')
        model=word[1].replace(ext,'')
        try:
            doit_two(run1,run2,model,outdir)
        except Exception as e:
            print('Failed to plot %s' % model)
            print(f"Failed with  error: {e}")
            






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)==3:
        do_all(sys.argv[1],sys.argv[2])
        print('Plots are in the director Xcompare')
    elif len(sys.argv)==4:
        doit_two(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print ('usage: regression_plot.py run1 run2 [model]')
