#!/usr/bin/env python 

'''
Perform a standardized comparison between two Python runs

Command line usage (if any):

    usage: xcompare.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

200629 ksl Coding begun

'''

import os
import matplotlib.pyplot as plt
from scipy.signal.windows import boxcar
from scipy.signal import convolve
import numpy as np
from astropy.io import ascii
import subprocess



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

def singlet(name=r'Ly$\alpha',wavelength=1216,bot=0.8,top=0.9):
    limits=plt.axis()
    # print(limits)
    if limits[0]<wavelength and wavelength<limits[1]:
        dy=limits[3]-limits[2]
        xbot=limits[2]+bot*dy
        xtop=limits[2]+top*dy
        plt.plot([wavelength,wavelength],[xbot,xtop],'k')
        # print(xbot,xtop)
        plt.text(wavelength,xtop,name,va='bottom',ha='center')
    return

        
def doublet(name='PV',w1=1118,w2=1128,bot=0.8,top=0.9):
    limits=plt.axis()
    # print(limits)
    if limits[0]<w1 and w2<limits[1]:
        dy=limits[3]-limits[2]
        xbot=limits[2]+bot*dy
        xtop=limits[2]+top*dy
        plt.plot([w1,w1],[xbot,xtop],'k')
        plt.plot([w2,w2],[xbot,xtop],'k')
        # print(xbot,xtop)
        plt.text(0.5*(w1+w2),xtop,name,va='bottom',ha='center')
    return
        
def add_lines():
    singlet(name=r'Ly$\alpha$',wavelength=1215.33,bot=0.8,top=0.9)
    singlet(name=r'Ly$\beta$',wavelength=1025.4,bot=0.8,top=0.9)
    singlet(name=r'Ly$\gamma$',wavelength=972.27,bot=0.8,top=0.85)
    singlet(name=r'Ly$\delta$',wavelength=949.48,bot=0.8,top=0.9)
    singlet(name=r'CIV',wavelength=1549,bot=0.8,top=0.9)
    singlet(name=r'CIII',wavelength=977.02008,bot=0.8,top=0.9)
    singlet(name=r'NIII',wavelength=989.7990118,bot=0.8,top=0.95)
    singlet(name=r'CIII$^{*}$',wavelength=1175,bot=0.8,top=0.9)
    singlet(name=r'NIV$^{*}$',wavelength=1718,bot=0.8,top=0.9)
    singlet(name=r'HeII',wavelength=1640,bot=0.8,top=0.9)
    singlet(name=r'NIV$^{*}$',wavelength=923,bot=0.8,top=0.9)
    singlet(name=r'OV$^{*}$',wavelength=1371.3,bot=0.8,top=0.9)
    singlet(name=r'NIV$^{*}$',wavelength=955.334,bot=0.7,top=0.75)
    
    
    doublet(name='PV',w1=1118,w2=1128,bot=0.8,top=0.9)
    doublet(name='NV',w1=1238.8,w2=1242.8,bot=0.8,top=0.9)
    doublet(name='SiIV',w1=1393.7,w2=1402.7,bot=0.8,top=0.9)
    doublet(name='SVI',w1=933.38,w2=945.55,bot=0.8,top=.85)
    doublet(name='OVI',w1=1031.9,w2=1037.6,bot=0.8,top=.85)
    doublet(name='SIV',w1=1062.6,w2=1076,bot=0.8,top=.85)

    return
    


def observer45(table):
    '''
    The column of the observer at 45 deg, whatever its azimuth
    (A45P0.50, or A45P180.0 as written by x3d)
    '''
    for one in table.colnames:
        if one.startswith('A45P') and not one.endswith(('_Q','_U')):
            return one
    raise KeyError('no observer at 45 deg in this spectrum')

def plot_spec(root='star',root2='star5',wmin=850,wmax=1850,smooth=21):
    '''
    Multi-panel plot comparing a CMFGen and Python model
    in the UV
    '''
    star=ascii.read(root+'.spec')
    star['nuFnu']=star['Lambda']*star[observer45(star)]
    star_nufnu=convolve(star['nuFnu'],boxcar(smooth)/float(smooth),mode='same')
    
    star2=ascii.read(root2+'.spec')
    star2['nuFnu']=star2['Lambda']*star2[observer45(star2)]
    star2_nufnu=convolve(star2['nuFnu'],boxcar(smooth)/float(smooth),mode='same')
    plt.figure(1,(8,12))
    
    
    plt.subplot(311)

    plt.plot(star['Lambda'],star_nufnu,label=root)
    plt.plot(star2['Lambda'],star2_nufnu,label=root2)
    plt.legend(loc='best')
    plt.xlim(850,1200)
    # plt.ylabel(r'$\nu F_{\nu}$ (ergs cm$^{-1}$s$^{-1}$)')
    # plt.xlabel(r'Wavelength ($\AA$)')
    add_lines()
    
    plt.subplot(312)

    plt.plot(star['Lambda'],star_nufnu,label=root)
    plt.plot(star2['Lambda'],star2_nufnu,label=root2)
    plt.legend(loc='best')
    plt.xlim(1150,1500)
    plt.ylabel(r'$\nu F_{\nu}$ (ergs cm$^{-1}$s$^{-1}$)',size=16)
    # plt.xlabel(r'Wavelength ($\AA$)')
    add_lines()
    
    plt.subplot(313)

    plt.plot(star['Lambda'],star_nufnu,label=root)
    plt.plot(star2['Lambda'],star2_nufnu,label=root2)
    plt.legend(loc='best')
    plt.xlim(1450,1800)
    # plt.ylabel(r'$\nu F_{\nu}$ (ergs cm$^{-1}$s$^{-1}$)')
    plt.xlabel(r'Wavelength ($\AA$)',size=16)
    add_lines()
    
    plt.savefig('%s_%s.png' % (root,root2))
    return

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


def plot_tot(xlist,plot_root=''):

    xtab=[]
    for one in xlist:
        x=ascii.read(one+'.log_spec_tot')
        xtab.append(x)

    fig=plt.figure(3,(8,12))
    plt.clf()



    plt.subplot(311)
    qmax=0
    i=0
    for data in xtab:
        root=xlist[i]
        created=data['Freq.']*xsmooth(data['Created'])
        plt.loglog(data['Freq.'],created,label='%s: Ext' %(root))
        wind=data['Freq.']*xsmooth(data['WCreated'])
        plt.loglog(data['Freq.'],wind,label='%s: Wind' % (root))
        q=np.fmax(created,wind)
        qmax=np.fmax(q,qmax)
        i+=1
    ymax=3*np.max(qmax)
    plt.ylim(ymax/1e8, ymax)
    test=data['Freq.'][qmax>ymax/1e8]
    plt.xlim(test[0],test[len(test)-1])

    plt.text(Lyman,0.9*ymax,'H',horizontalalignment='center',verticalalignment='top',size=14)
    plt.text(HeII,0.9*ymax,'HeII',horizontalalignment='center',verticalalignment='top',size=14)
    plt.plot([LymanA,LymanA],[0.5*ymax,0.9*ymax],'-k')
    plt.plot([HB,HB],[0.5*ymax,0.9*ymax],'-k')
    plt.plot([HA,HA],[0.5*ymax,0.9*ymax],'-k')
    plt.legend(loc='best')
    plt.title('Created')
    plt.ylabel(r'$\nu$L$_{\nu}$',size=16)

    plt.subplot(312)

    i=0
    for data in xtab:
        root=xlist[i]
        emitted=data['Freq.']*xsmooth(data['Emitted'])
        plt.loglog(data['Freq.'],emitted,label='%s: Total ' % root )
        wind=data['Freq.']*xsmooth(data['Wind'])
        plt.loglog(data['Freq.'],wind,label='%s: Wind' % root)
        i+=1

    plt.xlim(test[0],test[len(test)-1])
    plt.ylim(ymax/1e8, ymax)
    plt.legend(loc='best')
    plt.title('Observed')
    plt.ylabel(r'$\nu$L$_{\nu}$',size=16)

    plt.subplot(313)

    i=0
    for data in xtab:
        root=xlist[i]
        if np.average(data['CenSrc'])>0:
            censrc=data['Freq.']*xsmooth(data['CenSrc'])
            plt.loglog(data['Freq.'],censrc,label='%s: Star' % root,alpha=1.0)
        if np.average(data['Disk'])>0:
            disk=data['Freq.']*xsmooth(data['Disk'])
            plt.loglog(data['Freq.'],wind,label='%s: Disk' % root,alpha=1.0)
        hit=data['Freq.']*xsmooth(data['HitSurf'])
        plt.loglog(data['Freq.'],hit,label='%s: HitSurf' % root,alpha=1.0)
        i+=1

    plt.xlim(test[0],test[len(test)-1])
    plt.ylim(ymax/1e8, ymax)
    plt.legend(loc='best')
    plt.title('Observed')

    plt.xlabel(r'$ {\nu} $',size=16)
    plt.ylabel(r'$\nu$L$_{\nu}$',size=16)


    if plot_root=='':
        for one in xlist:
            if plot_root=='':
                plot_root='%s' % one
            else:
                plot_root='%s_%s' % (plot_root,one)

    plt.savefig(plot_root+'.tot_spec.png')

    return



    

def xplot(xlist,var='t_e',file='.master.txt'):
    plt.figure(2,(6,6))
    plt.clf()
    for one in xlist:
        filename=one+file
        x=ascii.read(filename)
        z=x[x['inwind']==0]
        plt.loglog(z['r'],z[var],label=one)
    plt.title(var)
    plt.ylabel(var)
    plt.xlabel('R')
    plt.legend(loc='best')
    plt.savefig('%s.png' % var)

    return

def doit (root1,root2):
    '''
    Run some standard plots that compare two-1d (for now) models
    '''


    windsave2table(root1)
    windsave2table(root2)

    plot_spec(root1,root2)

    xlist=[root1,root2]
    plot_tot(xlist)
    xplot(xlist,var='t_e',file='.master.txt')
    xplot(xlist,var='t_r',file='.master.txt')
    xplot(xlist,var='ne',file='.master.txt')

    return





        




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>2:
        doit(sys.argv[1],sys.argv[2])
    else:
        print ('usage: xcompare.py root1 root2')

