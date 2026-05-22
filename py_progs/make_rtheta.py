#!/usr/bin/env python 

'''
Create a simple rtheta import file that can be used
with Sirocco.

Synopsis:  

    Create an import file of an rtheta (polar) grid)
    with either angles in theta or angles
    set so that the size of cells at a single
    radius is the same


Command line usage:

    usage::

        make_rtheta.py file.pf

    where file.pf is a parameter file that gives 
    the information needed to create the grid. 
    
    If this file does not exist then one will be 
    asked for the necessary information, and the
    parameter file will be created (for future use).

Description:  

    The routine produces a spherically
    expanding wind on an r-theta grid, that can
    be imported into Sirocco.  

    The model is intended for testing various aspects of
    Sirocco, rather than being a model that is
    physically realistic.  

    Specifically, although the model has a
    velocity that increases linearly from an 
    inner to an outer radius. 

    The file that is created has a uniform density

    The rtheta grid can be one that has equal 
    angle bins, or one where the angular bins
    are adjusted so that the volume of the 
    cells at a particular radius are essentiall
    constant.

    Warning: The routine does not read or create 
    a Sirocco .pf file

    To make use of the file produced here, one
    needs to create an Sirocco pf file that uses
    the import option and has a coordinate sytem
    type of polar


Primary routines:

    doit

Notes:

    This routine was developed for testing of
    the convergence properties of rad/hydro 
    calculations 

    The routine makes use of a very simplified
    python implementation of rdpar.  

    It would be straightfoward to adapt this
    routine to produce a real stellar wind
    model with parameters that are consistent
    with what is used in Sirocco.
                                       
History:

    230601 ksl
        Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np
from kpar import *
from astropy.table import Table


def gen_model(x):
    '''
    Generate a model given a set of
    parameters
    '''

    nr=x['nr']
    ntheta=x['ntheta']
    r_min=x['r_min']
    r_max=x['r_max']
    v_min=x['v_min']
    v_max=x['v_max']
    rho=x['rho']
    t=x['t']

    imax=nr+2
    jmax=ntheta+2

    dlogr = (np.log10 (r_max / r_min)) / (nr - 3)

    ii=[]
    jj=[]

    i=0
    while i<imax:
        j=0
        while j<jmax:
            ii.append(i)
            jj.append(j)
            j+=1
        i+=1

    xtab=Table([ii,jj],names=['i','j'])

    xtab['inwind']=0
    n=0
    while n<len(xtab):
        if xtab['i'][n]==0 or xtab['i'][n]>=nr:
            xtab['inwind'][n]=-1
        if xtab['j'][n]==0 or xtab['j'][n]>=ntheta :
            xtab['inwind'][n]=-1
        n+=1
   

    xtab['r'] = r_min * 10.**(dlogr * (xtab['i'] - 1))
    xtab['r'].format='.4e'

    # For equal angles                                 

    if x['equal_angle']:

        dtheta=90./ntheta
        xtab['theta']=xtab['j']*dtheta

    else:
        theta=z=np.linspace(0.0,np.pi/2.,10000)
        z=np.sin(z)
        z=np.cumsum(z)/np.sum(z)

        # print('theta',theta)
        # print('z',z)
        # So we effectively have a cdf
        values=np.linspace(0.0,1.0,ntheta)
        # print('values',values)
        xtab['theta']=0.0
        kk=0
        jj=1
        xtheta=[-0.001,0]
        while kk<len(z) and jj<len(values):
            if z[kk]>values[jj]:
                # print('gotcha',kk,jj,theta[kk])
                xtheta.append(theta[kk])
                jj+=1
            kk+=1

        xtheta=np.array(xtheta)

        xtheta=xtheta*90./(0.5*np.pi)
        xtheta=np.append(xtheta,[90.,90.5])

        # print('final',len(xtheta),xtheta,jmax)

        n=0
        while n<len(xtab):
            xtab['theta'][n]=xtheta[xtab['j'][n]]
            n+=1



        print(kk,jj)

    xtab['v_x']=xtab['i']*(v_max-v_min)/nr+v_min
    xtab['v_y']=0.0
    xtab['v_z']=xtab['v_x']

    theta=0.5*np.pi/90.*xtab['theta']

    xtab['v_x']*=np.sin(theta)    
    xtab['v_z']*=np.cos(theta)

    xtab['theta'].format='.2f'
    xtab['v_x'].format=xtab['v_y'].format=xtab['v_z'].format='.2e'

    xtab['rho']=rho
    xtab['T']=t
    xtab['rho'].format=xtab['T'].format='.2e'


    outroot=x['root']
    outfile=outroot+'_import.txt'


    xtab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)

    return
    




def get_inputs(pf):
    '''
    Get the inputs
    '''

    opar(pf)

    nr=rdpar('No.Radial.Bins')
    ntheta=rdpar('No.Theta.Bins')
    print('yes - for standard equal angles, no for theta scaled so volumes at each radii are similar')
    equal_angle=rdpar('Equal.Angle')
    print('Remaining parameters should have cgs units')
    r_min=rdpar('R_min')
    r_max=rdpar('R_max')
    v_min=rdpar('v_min')
    v_max=rdpar('v_max')
    rho=rdpar('Rho')
    t=rdpar('T')
    outroot=rdpar('Outroot')

    if equal_angle.count('y') or equal_angle.count('Y'):
        equal_angle=True
    else:
        equal_angle=False

    cpar('foo.txt')

    x={'nr':nr,
            'ntheta':ntheta,
            'equal_angle':equal_angle,
            'r_min':r_min,
            'r_max':r_max,
            'v_min':v_min,
            'v_max':v_max,
            'rho':rho,
            't':t,
            'root':outroot}

    cpar(pf)

    return x




def doit(pf):
    '''
    Run the program
    '''

    x=get_inputs(pf)


    gen_model(x)

    print(x)





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__)
