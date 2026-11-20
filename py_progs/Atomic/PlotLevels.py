#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Make a plot of the levels for a ion as a function of excitation energy and 
ISLP


Command line usage (if any):

    usage: PlotLevels.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

231012 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii
from astropy.table import join
import matplotlib.pyplot as plt
import numpy as np



def plot_config(ion='he_1',annotate=True):
    '''
    where 

    ion is a ion name, eg he_1
    and
    annotate menast to label the levels
    '''
    
    lev_file=ion+'_levels.dat'
    try:
        levels=ascii.read(lev_file)
    except:
        print('Could not read levels from %s' % lev_file)
        return
    
    line_file=ion+'_lines.dat'
    try:
        lines=ascii.read(line_file)
    except:
        print('Could not read levels from %s' % line_file)
        return    
    
    xup=levels['lvl','islp']
    xlow=xup.copy()

    xup.rename_column('lvl','ul')
    xup.rename_column('islp','u_islp')
    xlow.rename_column('lvl','ll')
    xlow.rename_column('islp','l_islp')
    
    lines=join(lines,xlow,join_type='left')
    lines=join(lines,xup,join_type='left')
    
    plt.figure(1,(6,10))
    plt.plot(levels['islp'],levels['ex']-levels['ex'][0],'.')
    if annotate:
        for one in levels:
            plt.annotate(one['lvl'],[one['islp'],one['ex']-levels['ex'][0]])

            
    for one in lines:
        xx=[one['l_islp'],one['u_islp']]
        yy=[one['el'],one['eu']]
        xalpha=one['f']
        if xalpha>1:
            xalpha=1
        if xalpha>0.01:
            plt.plot(xx,yy,'k',alpha=xalpha)
        elif xalpha>0:
            plt.plot(xx,yy,':k',alpha=0.1)
                          
        
        plt.xlabel('ISLP',size=16)
        plt.ylabel('Excitation Energy (eV)',size=16)
        plt.title(ion,size=16)
        plt.tight_layout()  
        plt.savefig(ion+'_config.png')



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        plot_config(sys.argv[1])
    else:
        print (__doc__)




