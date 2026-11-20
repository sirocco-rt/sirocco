#!/usr/bin/env python

from numpy import *
from pylab import *



def read_python_spec(filename,column):
    """Read the a spectrum file created with python"""
    try: 
        specfile=open(filename,'r')
    except  IOError:
        print((filename, " does not exist"))
        return(-1)
    w=[]
    f=[]
    line=specfile.readline()
    while line!='':
        z=line.split()
        word1=z[0]
        if word1[0]=='#':
            pass
        else:
            w=w+[float(z[1])]
            f=f+[float(z[column+1])]
        line=specfile.readline()
    return w, f
    
    
def read_data_spec(filename,wave,flux,error):
    """Read the ascie data spectrum file"""
    try: 
        specfile=open(filename,'r')
    except  IOError:
        print(filename," does not exist")
        return(-1)
    n=0
    z=specfile.readline()
    while (z != ""):
        q=z.split()
        word1=q[0]
        if word1[0] != '#':
            wave.extend([float(q[0])])
            flux.extend([float(q[1])])
            if len(q) > 2:
                error.extend([float(q[2])])
                ierr=1
            else:
                error.extend([1.0])
                ierr=0
        n=n+1
        z=specfile.readline()
    return filename,n,ierr
    
    
def ave (x,y,xmin,xmax):
    """average a restricted portion of y between xmin and xmax"""
    lx=len(x)
    num=0
    ytot=0
    l=0
    while l < lx :
        if(xmin <= x[l] <= xmax):
            ytot=ytot+y[l]
            num=num+1
        l=l+1
    if (num==0): return(-99)
    return(ytot/num)

def rescale (zin,zout,scale):
    """rescale a list by a scalefactor"""
    lin=len(zin)
    i=0
    while i < lin:
        zout.extend([scale*zin[i]])
        i=i+1
#        print 'len zout', len(zout)
    return len(zout)


def scale_model (wdat,fdat,wmod,fmod,wmin,wmax,fout):
    """Produce a rescaled model"""
    data_ave=ave(wdat,fdat,wmin,wmax)
#    print data_ave
    mod_ave=ave(wmod,fmod,wmin,wmax);
#    print mod_ave
    scale=data_ave/mod_ave
    rescale(fmod,fout,scale)
#    print 'len fout', len(fout)
    return scale


def plot_data_model(data,model,column,wmin,wmax):
    """plot a rescaled model against the data"""
    wave=[]
    flux=[]
    error=[]
    if read_data_spec(data,wave,flux,error) == -1:
        return -1
#        print 'Read data',len(wave),len(flux)
    wmod=[]
    fmod=[]
    if read_python_spec(model,column,wmod,fmod) == -1:
        return -1
#    print 'Read model',len(wmod),len(fmod)
    fout=[]
    scale_model (wave,flux,wmod,fmod,wmin,wmax,fout)
#    print 'scale_model', len(fout)
    plot(wave,flux,'b')
    plot(wmod,fout,'r')
    reset=axis()
    xreset=list(reset)
    xreset[0]=wmin
    xreset[1]=wmax
    xreset[2]=0
    axis(xreset)
    draw()
    

