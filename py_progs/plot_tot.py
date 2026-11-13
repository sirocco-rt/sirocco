#!/usr/bin/env python 
'''
Plot the emergent spectrum from the log_spec_tot file


Command line usage (if any):

    usage::

        plot_tot.py [-smooth 11] rootname


    rootname:
        is the rootname of the files in the run
    `-smooth`:
        is an opticnal parameter indicating how
        much smooth of the orignnal spectrum is to
        be done.  The default is 11 bins


Description:  
    This routines plots to a file a boxcar smoothed version of the
    log spectrum. 
    

Primary routines:

Notes:
                                       
History:

    130620    ksl
        Coding begun
    141125    ksl
        Updated to use astropy.io

'''

import sys
import os
import numpy
import pylab
from astropy.io import ascii
from scipy.signal.windows import boxcar
from scipy.signal import convolve


C=2.997e18 # c in Angstroms/s
Lyman=C/912.
HeII=C/912.*4
LymanA=C/1216
HB=C/4863
HA=C/6563

def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux
    '''

    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)




def doit(rootname='sphere',smooth=21,fig_no=2,outdir=''):
    '''
    Plot the spectra contained in the spec_tot 
    file

    141125    ksl
        Updated for new formats which use astropy
    191210  ksl
        Modified so that what is ploted is nuL_nu, and
        did a better job at setting limits for the plot
    261018  ksl
        The plot is written to outdir, if given
    '''
    # Make sure we only have the rootname

    rootname=rootname.replace('.log_spec_tot','')
    filename=rootname+'.log_spec_tot'

    try:
        data=ascii.read(filename)
    except IOError:
        print('Error: Could not find %s' % filename)
        return

    # print(data.colnames)

    
    pylab.figure(fig_no,(12,6))
    pylab.clf()

    created=data['Freq.']*xsmooth(data['Created'])
    pylab.loglog(data['Freq.'],created,label='Created(Ext)')

    try:
        wind=data['Freq.']*xsmooth(data['WCreated'])
        pylab.loglog(data['Freq.'],wind,label='Created(Wind)')
    except:
        pass



    emitted=data['Freq.']*xsmooth(data['Emitted'])
    hit_surf=data['Freq.']*xsmooth(data['HitSurf'])
    wind=data['Freq.']*xsmooth(data['Wind'])

    pylab.loglog(data['Freq.'],emitted,label='Observed(Tot)')

    pylab.loglog(data['Freq.'],wind,label='Observed(Wind)')

    pylab.loglog(data['Freq.'],hit_surf,label='Hit Surface')

    zz=pylab.axis()

    # Find the maximum values of all of these arrays

    q=numpy.fmax(created,emitted)
    qq=numpy.fmax(q,wind)
    ymax=3*numpy.max(qq)
    pylab.ylim(ymax/1e8, ymax)

    # make a mask of from the values of qq

    test=data['Freq.'][qq>ymax/1e8]
    pylab.xlim(test[0],test[len(test)-1])

    # pylab.xlim(zz[0],zz[1])


    # pylab.axis((zz[0],zz[1],zz[3]/1e8,zz[3]))

    pylab.text(Lyman,0.9*ymax,'H',horizontalalignment='center',verticalalignment='top',size=14)
    pylab.text(HeII,0.9*ymax,'HeII',horizontalalignment='center',verticalalignment='top',size=14)
    pylab.plot([LymanA,LymanA],[0.5*ymax,0.9*ymax],'-k')
    pylab.plot([HB,HB],[0.5*ymax,0.9*ymax],'-k')
    pylab.plot([HA,HA],[0.5*ymax,0.9*ymax],'-k')
    pylab.legend(loc='best')
    pylab.title(rootname,size=16)
    pylab.xlabel(r'$ {\nu} $',size=16)
    pylab.ylabel(r'$\nu$L$_{\nu}$',size=16)


    pylab.draw()
    pylab.savefig(os.path.join(outdir,rootname+'.spec_tot.png'))

    # Since the total spectra are basically written out as L_nu we need to 
    # integrate over frequency

    freq=numpy.array(data['Freq.'])

    dfreq=[]
    i=0
    while i<len(freq):
        if i==0:
            dfreq.append(freq[1]-freq[0])
        elif i==len(freq)-1:
            dfreq.append(freq[i]-freq[i-1])
        else:
            dfreq.append(0.5*(freq[i+1]-freq[i-1]))

        i+=1

    dfreq=numpy.array(dfreq)

    lum_created=numpy.sum(dfreq*data['Created'])
    lum=numpy.sum(dfreq*data['Emitted'])
    
    #dfreq=freq[1]-freq[0]  # For a linear scale the frequencies are all equally spaced

    # lum_created=dfreq*numpy.sum(numpy.array(data['Created']))
    # lum=dfreq*numpy.sum(data['Emitted'])
    # print(data['Created'])
    print('The Created luminosity was ',lum_created)
    print('The emitted luminosity was ',lum)

    return


def steer(argv):
    '''
    Parse the command line
    '''

    smooth=11

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-smooth':
            i+=1
            smooth=eval(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            rootname=argv[i]
        i+=1

    doit(rootname,smooth)

    return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print(__doc__)
