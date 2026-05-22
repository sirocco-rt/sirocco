#!/usr/bin/env python 
'''
Synopsis:  
    These are routines for plotting various parameters in of the wind
    after these parameters have been saved to an astropy-compatible
    ascii table


Command line usage 
    plot_wind filename var   

    to make a plot of a single variable from the command line

Description:  


Primary routines:
    doit : Create a plot of a single variable in a file made with 
            windsave2table.  This is the routine called from
            the command line. Additional options are available
            when called from a python script.
    compare_separate    : Compare a single variable in two
            different runs of python and produce 3 separate plots
            one for each run and one containing the difference
    compare:     Similar to compare_separate but produces a 
            single file
'''

import sys
from astropy.io import ascii
import numpy
import matplotlib.pyplot as pylab
from scipy.signal.windows import boxcar
from scipy.signal import convolve
import subprocess
from matplotlib import tri
from astropy.io import ascii
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy



def compare_separate(f1='fiducial_agn_master.txt',f2='fiducial_agn_master.txt', var='t_r',grid='ij',inwind='',scale='guess',zmin=-1e50,zmax=1e50,fig_no=5):
    '''
    This routine compares the same variable from two different runs of python, and produces
    Three separate plots. The plots represent the variable in the first file, the variable in the second file and
    the difference between the two
    '''

    x1,y1,xvar1,title,xlabel,ylabel=get_data(f1,var,grid,inwind,scale,zmin,zmax)
    x2,y2,xvar2,title,xlabel,ylabel=get_data(f2,var,grid,inwind,scale,zmin,zmax)
    xdiff=xvar2-xvar1

    root=f1.split('.')
    root=root[0]
    just_plot(x1,y1,xvar1,root,title,xlabel,ylabel,fig_no)

    root=f2.split('.')
    root=root[0]
    just_plot(x2,y2,xvar2,root,title,xlabel,ylabel,fig_no+1)

    root='diff'
    just_plot(x1,y1,xdiff,root,title,xlabel,ylabel,fig_no+2)
    return

def compare(f1='sv_master.txt',f2='sv_master.txt',var='t_r',grid='ij',inwind='',scale='guess',zmin=-1e50,zmax=1e50,fig_no=5):
    '''
    Compare results of two variables within a single plot
    

    The three plots are of from the first file, the second file, and the difference 
    respectively
    '''


    
    x1,y1,xvar1,title,xlabel,ylabel=get_data(f1,var,grid,inwind,scale,zmin,zmax)
    x2,y2,xvar2,title,xlabel,ylabel=get_data(f2,var,grid,inwind,scale,zmin,zmax)
    xdiff=xvar2-xvar1

    

    pylab.figure(fig_no,(18,6))
    pylab.clf()
    pylab.subplot(131)

    ax=pylab.gca()
    cmap=copy.copy(pylab.cm.jet)
    cmap.set_bad(color='white')
    im=ax.pcolormesh(x1,y1,xvar1,cmap=cmap)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pylab.colorbar(im, cax=cax)
    pylab.tight_layout()


    pylab.subplot(132)

    ax=pylab.gca()
    cmap=copy.copy(pylab.cm.jet)
    cmap.set_bad(color='white')
    im=ax.pcolormesh(x2,y2,xvar2,cmap=cmap)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pylab.colorbar(im, cax=cax)
    pylab.tight_layout()

    pylab.subplot(133)

    ax=pylab.gca()
    cmap=copy.copy(pylab.cm.jet)
    cmap.set_bad(color='white')
    im=ax.pcolormesh(x2,y2,xdiff,cmap=cmap)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pylab.colorbar(im, cax=cax)
    pylab.tight_layout()


    pylab.draw()
    
    outroot='Difference_%s.png' % var
    pylab.savefig(outroot)


    return
    

def get_data(filename='fiducial_agn_master.txt', var='t_r',grid='ij',inwind='',scale='guess',zmin=-1e50,zmax=1e50):
    '''
    This routine reads and scales the data from a single variable that is read from an ascii table 
    representation of one or more of the parameters in a windsave file

    The grid obptions:

    ij
    log
    lin

    '''

    try:
        data=ascii.read(filename)
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        return

    # data.info()

    # print('options: xyscale %s varscale %s min %g max %g' % (grid,scale,zmin,zmax))


    title=var


    ii=data['i']
    jj=data['j']
    iimax=numpy.max(ii)
    jjmax=numpy.max(jj)
    if grid=='ij':
        x=numpy.array(data['i'])
        y=numpy.array(data['j'])
        xlabel='i'
        ylabel='j'
    elif grid=='log':
        x=numpy.array(data['x']) 
        y=numpy.array(data['z'])
        # find the minimum value in x and y that is greater than 0
        xmin=1e50
        ymin=1e50
        i=0
        while i<len(x):
            if x[i]>1 and x[i]<xmin:
                xmin= x[i]
            if y[i]>1 and y[i]<ymin:
                ymin= y[i]
            i+=1
        xlogmin=numpy.log10(xmin/10)
        ylogmin=numpy.log10(xmin/10)
        x=numpy.select([x>1],[numpy.log10(x)],default=xlogmin)
        y=numpy.select([y>1],[numpy.log10(y)],default=ylogmin)
        xlabel='log(x)'
        ylabel='log(z)'
    else:
        x=numpy.array(data['x'])
        y=numpy.array(data['z'])
        xlabel='x'
        ylabel='z'

    xvar=numpy.array(data[var])
    xwind=numpy.array(data['inwind'])
    x=x.reshape(iimax+1,jjmax+1)
    y=y.reshape(iimax+1,jjmax+1)
    xvar=xvar.reshape(iimax+1,jjmax+1)
    xwind=xwind.reshape(iimax+1,jjmax+1)


    # now mask out the appropriate data

    if inwind=='all':
        xwind.fill(0) # Keep everything
    elif inwind.count('partial'):
        xwind=numpy.select([xwind<0],[xwind],default=0)

    # Now make the mask regardless
    mask=numpy.ma.make_mask(xwind)





    xvar=numpy.ma.array(xvar,mask=mask)


    vmin=numpy.ma.min(xvar)
    vmax=numpy.ma.max(xvar)



    if vmin<zmin:
        vmin=zmin
    if vmax>zmax:
        vmax=vmax



    if scale=='guess':
        if (vmax-vmin)>1e4:
            scale='log'
    
    # Need this to go from 
    if scale=='log':
        if vmin<vmax/1e10:
            vmin=vmax/1e10
        xvar=numpy.select([xvar>vmin],[xvar],default=vmin)
        xvar=numpy.ma.log10(xvar)
        vmax=numpy.log10(vmax)
        vmin=numpy.log10(vmin)
        title='log '+title

    # Reapply the mask; need to make bad color work
    xvar=numpy.ma.array(xvar,mask=mask)



    


    return x,y,xvar,title,xlabel,ylabel


def just_plot(x,y,xvar,root,title,xlabel,ylabel,fig_no=1,vmin=0,vmax=0):
    '''
    This routine simply is produces a plot of a variable from
    that has been printed to an astropy table with a routine like
    windsave2table.    This function is simply a plotting routine


    '''

    pylab.figure(fig_no,(6,6))
    pylab.clf()
    ax=pylab.gca()
    cmap=copy.copy(pylab.cm.jet)
    cmap.set_bad(color='white')

    ax.tick_params(labelsize=14)
    if vmin==0 and vmax==0:
        im=ax.pcolormesh(x,y,xvar,cmap=cmap,shading='auto')
    else:
        im=ax.pcolormesh(x,y,xvar,cmap=cmap,shading='auto',vmin=vmin,vmax=vmax)

    pylab.title(title,size=16)
    pylab.xlabel(xlabel,size=16)
    pylab.ylabel(ylabel,size=16)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pylab.colorbar(im, cax=cax)

    pylab.draw()
    title=title.replace(' ','_')
    if root=='':
        filename=title+'.png'
    else:
        filename='%s_%s.png' %(root,title)
    pylab.tight_layout()
    pylab.savefig(filename)
    return filename



def doit(filename='fiducial_agn.master.txt', var='t_r',grid='ij',inwind='',scale='guess',zmin=-1e50,zmax=1e50,
        plot_dir='',root=''):
    '''
    Plot a single variable from an astropy table (normally created with windsave2table, with various
    options

    where var is the variable to plot where grid can be ij, log, or anything else.  If ij then the plot will 
    be in grid coordinates, if log the plot will be in on a log scale in physical coordiantes.  If anything else, 
    the plot will be on a linear scale in physical coordiantes where scale indicates how the variable should be 
    plotted.  guess tells the routine to make a sensible choice linear implies the scale should be linear and log 
    implies a log scale should be used where zmin and zmax overide the max and mimimum in the array (assuming these 
    limits are with the range of the variable)
    '''

    if root=='':
        root=filename.split('.')
        root=root[0]

    x,y,xvar,title,xlabel,ylabel=get_data(filename,var,grid,inwind,scale,zmin,zmax)
    if plot_dir!='':
        root='%s/%s' % (plot_dir,root)
    # print(len(x),len(y),len(xvar))

    if zmin==-1e50 and zmax==1e50:
        autoscale=True
        vmin=vmax=0
    else:
        autoscale=False
        vmin=zmin
        vmax=zmax
        if title.count('log'):
            vmin=numpy.log10(vmin)
            vmax=numpy.log10(vmax)

    plotfile=just_plot(x,y,xvar,root,title,xlabel,ylabel,vmin=vmin,vmax=vmax)

    return plotfile

    
def steer(argv):
    '''
    Controlling routine for py_wind
    '''

    xfile=''
    xvar=''
    xgrid='ij'
    xmin=-1e50
    xmax=1e50
    xscale='guess'

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-log':
            xgrid='log'
        elif argv[i][0:4]=='-lin':
            xgrid='linear'
        elif argv[i][0:4]=='-min':
            i+=1
            xmin=eval(argv[i])
        elif argv[i][0:4]=='-max':
            i+=1
            xmax=eval(argv[i])
        elif argv[i][0]=='-':
            print('Unknown options :',argv)
            return
        elif xfile=='':
            xfile=argv[i]
        elif xvar=='':
            xvar=argv[i]
        i+=1


    doit(filename=xfile, var=xvar,grid=xgrid,inwind='',scale='log',zmin=xmin,zmax=xmax,
        plot_dir='',root='')


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>2:
         steer(sys.argv)
    else:
        print(__doc__)
