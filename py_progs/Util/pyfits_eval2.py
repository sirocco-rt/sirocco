#!/usr/bin/env python

import tkinter
from pylab import *
ion() # Turn interactive graphics on



def read_python_spec(filename,column,wave,flux):
    """Read the a spectrum file created with python"""
    try: 
        specfile=open(filename,'r')
    except  IOError:
        print(filename, " does not exist")
        return(-1)
    n=0
    z=specfile.readline()
    q=z.split()
    specname=q[column]
    z=specfile.readline()
    while (z != ""):
        q=z.split()
#        wave[n]=float(q[1])
        wave.extend([float(q[1])])
        if(len(q)< (column-1)):
            print("Not enough columns in line",n)
            return(-1)
#        flux[n]=float(q[column])
        flux.extend([float(q[column])])
        n=n+1
        z=specfile.readline()
    return(specname)
    
    
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
#    clf()
    plot(wave,flux,'b')
    plot(wmod,fout,'r')
    reset=axis()
    xreset=list(reset)
    xreset[0]=wmin
    xreset[1]=wmax
    xreset[2]=0
    axis(xreset)
    

# This routine determines which of the models are good enough
# and shows where in the grid they lie
def analyze_pyfit3_results(filename,col,cutoff,col2,cutoff2,nx,ny):
    """Identify which models are good enough"""    
    # Assume columns in input files are numbered from 1
    print("Analyze",col,cutoff,col2,cutoff)
    try: 
        resultsfile=open(filename,'r')
    except  IOError:
        print(filename, " does not exist")
        return(-1)
    z=resultsfile.readline()
    value=[]
    value2=[]
    x=[]
    y=[]
    xgood=[]
    ygood=[]
    goodmodels=[]
    goodspectrum=[]
    n=0
    ngood=0
    if col > 0:
        col=col-1
    if col2 > 0:
        col2=col2-1
    xaxislabel='Col '+str(nx)
    yaxislabel='Col '+str(ny)
    nx=nx-1
    ny=ny-1
    while(z!=''):
        q=z.split()
        firstword=q[0]
        if firstword[0]== '#':
            print(z.strip())
            # It's a comment line
            secondword=q[1]
            if secondword == 'colmin':
                colmin=int(q[2])
                colmax=int(q[4])
                print('colmin,colmax',colmin,colmax)
            if secondword == 'col.with.modelname':
                modelcolumn=int(q[2])-1
        else:
            value.extend([q[col]])
            x.extend([float(q[nx])])
            y.extend([float(q[ny])])
            value2.extend([q[col2]])
            if float(value[n]) < float(cutoff) and float(value2[n])<float(cutoff2):
                xgood.extend([float(q[nx])])
                ygood.extend([float(q[ny])])
                goodmodels.extend([q[modelcolumn]])
                goodspectrum.extend([int(q[modelcolumn+1])])
                ngood=ngood+1
            # It's a data line
            n=n+1
        z=resultsfile.readline()
    print('For a cutoff of ',cutoff,'there are ', ngood,'good models')
    figure(1)
    clf()
    plot(x,y,'ro')
    plot(xgood,ygood,'bo')
    aaa=axis()
    print("trying to rock and roll", aaa)
    bbb=list(aaa)
    bbb[0]=0.9*bbb[0]
    bbb[1]=1.1*bbb[1]
    bbb[2]=0.9*bbb[2]
    bbb[3]=1.1*bbb[3]
#    print aaa
#    axis(aaa)
    axis(bbb)
    font={'fontname'   : 'Helvetica', 'color'      : 'b', 'fontweight' : 'bold', 'fontsize'   : 16}
    xlabel(xaxislabel,font)
    ylabel(yaxislabel,font)
    return [goodmodels,goodspectrum]
    
    
# Plot the good models
def plot_good(data,models,wmin,wmax):
    """Compare the good models to the data"""
    n=0;
    modelnames=models[0]
    modelcolumns=models[1]
    lin=len(modelcolumns)

    print("plot_good",data,wmin,wmax,"Nmods to plot:  ",lin)
#        print models

    figure(2)
    clf()
    while n < lin:
        model=modelnames[n]
        col=modelcolumns[n]+8
        print('model: ',model,col)
        plot_data_model(data,model,col,wmin,wmax)
        n=n+1    
        font={'fontname'   : 'Helvetica', 'color'      : 'b', 'fontweight' : 'bold', 'fontsize'   : 16}
    xlabel('Wavelength',font)
    ylabel('Flux',font)
    return(lin)    
# OK, here is the gui that runs it all

def gui():
    global v1,v2,v3,vv2,vv3,v4,v5,v6,v7,v8
    global xmodels

    xmodels=[]



    def process1():
        global v1,v2,v3,vv2,vv3,v4,v5,v6,v7,v8
        global xmodels
        xmodels=[13,14,16]
        filename=v1.get()
        col=v2.get()
        cutoff=v3.get()
        col2=vv2.get()
        cutoff2=vv3.get()
        nx=v4.get()
        ny=v5.get()
        print("filename: ",filename)
        print("Col containing chi**2: ",col)
        xmodels=analyze_pyfit3_results(filename,col,cutoff,col2,cutoff2,nx,ny)
#        print "Number of models returned",len(xmodels)
    

# plot the good models
    def process2():
        global v1,v2,v3,vv2,vv3,v4,v5,v6,v7,v8
        global xmodels
        if len(xmodels)==0:
            print("Error: no models to compare. Run select best")
            return(0)
        data=v6.get()
        wmin=v7.get()
        wmax=v8.get()
#        print data,wmin,wmax
#        print xmodels
        plot_good(data,xmodels,wmin,wmax)
    

    def destroy():
        root.destroy()   # Close the gui
        print("on destroy", len(xmodels))
        close("all")   # Close all the plot windows


    root = tkinter.Tk()
    frame=tkinter.Frame(root)
    frame.pack()
    b1 = tkinter.Button(frame,text="Select best",command=process1)
    b2 = tkinter.Button(frame,text="Plot   best",command=process2)
    b3 = tkinter.Button(frame,text="Quit",command=destroy)

    l1=tkinter.Label(frame,text="Results from pyfit3: ")

    l2=tkinter.Label(frame,text="Column to use for chi**2: ")
    l3=tkinter.Label(frame,text="Value of scaled chi**2 to keep: ")

    ll2=tkinter.Label(frame,text="Column to use for chi**2: ")
    ll3=tkinter.Label(frame,text="Value of scaled chi**2 to keep: ")

    l4=tkinter.Label(frame,text="Col x to plot: ")
    l5=tkinter.Label(frame,text="Col y to plot: ")
    l6=tkinter.Label(frame,text="Spectrum: ")
    l7=tkinter.Label(frame,text="Wave min: ")
    l8=tkinter.Label(frame,text="Wave max: ")

    b1.grid(row=0, column=0)
    b2.grid(row=0, column=1)
    b3.grid(row=12, column=1)

    l1.grid(row=2,column=0)

    l2.grid(row=3,column=0)
    l3.grid(row=4,column=0)

    ll2.grid(row=5,column=0)
    ll3.grid(row=6,column=0)

    l4.grid(row=7,column=0)
    l5.grid(row=8,column=0)
    l6.grid(row=9,column=0)
    l7.grid(row=10,column=0)
    l8.grid(row=11,column=0)

    # Data set to analyze
    # Pyfit3 results file
    v1=tkinter.StringVar()
    v1.set("v3885sgr_hst_ixvel.results")
    e1=tkinter.Entry(frame,textvariable=v1)
    e1.grid(row=2,column=1)

    
        # Column to use
    v2=tkinter.IntVar()
    v2.set(5)
    e2=tkinter.Entry(frame,textvariable=v2)
    e2.grid(row=3,column=1)

        # Ratio of chi**2 to keep
    v3=tkinter.DoubleVar()
    v3.set(1.1)
    e3=tkinter.Entry(frame,textvariable=v3)
    e3.grid(row=4,column=1)

        # Column to use
    vv2=tkinter.IntVar()
    vv2.set(6)
    ee2=tkinter.Entry(frame,textvariable=vv2)
    ee2.grid(row=5,column=1)

        # Ratio of chi**2 to keep
    vv3=tkinter.DoubleVar()
    vv3.set(1.1)
    ee3=tkinter.Entry(frame,textvariable=vv3)
    ee3.grid(row=6,column=1)


        # Column x to plot       
    v4=tkinter.IntVar()
    v4.set(7)
    e4=tkinter.Entry(frame,textvariable=v4)
    e4.grid(row=7,column=1)

        # Column y to plot        
    v5=tkinter.IntVar()
    v5.set(8)
    e5=tkinter.Entry(frame,textvariable=v5)
    e5.grid(row=8,column=1)

    v6=tkinter.StringVar()
    v6.set("v3885sgr_hst1")
    e6=tkinter.Entry(frame,textvariable=v6)
    e6.grid(row=9,column=1)


    v7=tkinter.DoubleVar()
    v7.set(1500)
    e7=tkinter.Entry(frame,textvariable=v7)
    e7.grid(row=10,column=1)

    v8=tkinter.DoubleVar()
    v8.set(1600)
    e8=tkinter.Entry(frame,textvariable=v8)
    e8.grid(row=11,column=1)






    root.title("Pyfit3 Eval")
    root.mainloop()

if __name__ == "__main__":        # allows one to run from command line without running automatically with write_docs.py

    gui()
