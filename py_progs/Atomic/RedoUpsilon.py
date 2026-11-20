#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Reformat collision files  so that one simply reads in
the collision strength as a functions of kT/hu


Command line usage (if any):

    usage: RedoUpslon.py old_upsilon file

Description:  

    The routine reads the collision strength files
    and reworks all of the collision strenths to
    be simple linear interpolations 

    This type of collision strength is given
    the label type 5 to indicate that this
    is the format.

    The routine writes a new file, appending
    _X to the root

    A plot file is also generated that shows
    how the collision strength varies.


Primary routines:

    doit

Notes:

    The routine presumes that the input 
    files are "pure", in the sense
    that they do not have extra 
    embedded comments.
                                       
History:

240815 ksl Coding begun

'''

import os
from astropy.io import ascii
from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
import sys

# # Reformat Collision Files
# 
# Create a routine to reformat the current upsilon files to a new format that just uses interplolation



# ```
#    "%*s %*s %d %2d %le %le %le %le %le %le %d %d %d %d %le %le %le %d %d %le",
#     &z, &istate, &wave, &f, &gl, &gu, &el, &eu, &levl, &levu, &c_l, &c_u, &en, &gf, &hlt, &np, &type, &sp));
# ```


def read_collisions(col='he_1__upsilon.dat',return_original=False):
    '''
    Read the collision file and put everything into table for retrieval

    Note that this table cannot be written in ascii.fixed_width_two_Line, one must used something like ascii.ecv

    If return_original = True than the routine returns the lines it has read in as well as the table that
    one would use going forward so that the data can be subsetted easily
    '''

    x=open(col)
    lines=x.readlines()
    print(len(lines))
    cstren_line=[]
    sct_line=[]
    scups_line=[]
    for one_line in lines:
        if one_line.count('CSTREN'):
            cstren_line.append(one_line.strip())
        elif one_line.count('SCT'):
            sct_line.append(one_line.strip())
        elif one_line.count('SCUPS'):
            scups_line.append(one_line.strip())
        else:
            print('Could not interpret ',one_line)

    print(len(cstren_line),len(sct_line),len(scups_line))

    element=[]
    ion=[]
    f=[]
    gl=[]
    gu=[]
    el=[]
    eu=[]
    lower=[]
    upper=[]
    cl=[]
    cu=[]
    xtype=[]
    ntemp=[]
    scaling_param=[]
    tlim=[]
    energy=[]
    gf=[]
    xlam=[]
    sct=[]
    scups=[]
    for one in cstren_line:
        word=one.split()
        # print(len(word),word)
        element.append(int(word[2]))
        ion.append(int(word[3]))
        xlam.append(float(word[4]))
        f.append(float(word[5]))
        gl.append(int(word[6]))
        gu.append(int(word[7]))
        el.append(float(word[8]))
        eu.append(float(word[9]))
                  
                        
        lower.append(int(word[10]))
        upper.append(int(word[11]))

        cl.append(int(word[12]))
        cu.append(int(word[13]))
        
        energy.append(float(word[14]))
        gf.append(float(word[15]))
        tlim.append(float(word[16]))
        ntemp.append(int(word[17]))
        xtype.append(int(word[18]))
        scaling_param.append(float(word[-1]))

    coltab=Table([element,ion,xlam,f,gl,gu,el,eu,lower,upper,cl,cu, energy,gf, xtype,ntemp,tlim,scaling_param],
                 names=['Element','Ion','Wave','f','gl','gu','el','eu','ll','ul','cl','cu','E_ryd','gf','type','ntemp','temp_lim','scale_par'])




    sct=[]
    for one in sct_line:
        word=one.split()
        j=1
        record=[]
        while j<len(word):
            record.append(float(word[j]))
            j+=1
        record=np.array(record)
        sct.append(record)


    scups=[]
    for one in scups_line:
        word=one.split()
        j=1
        record=[]
        while j<len(word):
            record.append(float(word[j]))
            j+=1
        record=np.array(record)
        scups.append(record)



    coltab['sct']=sct
    coltab['scups']=scups


    coltab['Number']=range(len(coltab))

    coltab.write('goo.txt',format='ascii.ecsv',overwrite=True)

    if return_original==True:
        return coltab,cstren_line,sct_line,scups_line
    else:
        return coltab


def convert_ups(xtable,xratio=30,npts=20):
    '''
    Process a single row, as a function of temperature
    '''
    ytable=xtable.copy()
    ytable.remove_column('sct')
    ytable.remove_column('scups')


    sct=[]
    scups=[]
    u0=np.linspace(0,xratio,npts)
    i=0
    while i <len(xtable):
        # print(i)
        row=xtable[i]
     
        if row['type']==5:
            # The format has alredy been converted
            print('what')
            continue
        elif row['type']==1 or row['type']==4:
            x = 1. - (np.log (row['scale_par']) / np.log (u0 + row['scale_par']))
        elif row['type'] == 2 or row['type'] == 3:
            x = u0 / (u0 + row['scale_par'])
        else:
            print('Error - Unknown collision type')
            print (row)
            return

        # print(x)

        z=np.interp(x,row['sct'],row['scups'])
        # plt.plot(u0,z/z[0])

        sct.append(u0)
        scups.append(z)
        i+=1

    sct=np.array(sct)
    scups=np.array(scups)
    
    # print(sct.shape) 
    # print(scups.shape)

    # print(len(ytable))
    
    ytable['sct']=sct


    ytable['scups']=scups
    ytable['type']=5
    ytable['ntemp']=npts
    ytable['scale_par']=-99.
    return ytable
    



#    "%*s %*s %d %2d %le %le %le %le %le %le %d %d %d %d %le %le %le %d %d %le",
#     &z, &istate, &wave, &f, &gl, &gu, &el, &eu, &levl, &levu, &c_l, &c_u, &en, &gf, &hlt, &np, &type, &sp));



def plot_col(xtable):
    plt.figure(1,(6,6))
    for one  in xtable:
        try:
            plt.plot(one['sct'],one['scups']/one['scups'][0])
        except:
            print('Could not plot ',one)
    plt.xlabel('kT/E')
    plt.ylabel('C(kT/E):C(kT=0)')
    plt.tight_layout()

        

def rewrite(xtable,outname='toot.txt'):
    f=open(outname,'w')
    for row in xtable:
        string='CSTREN Line %5d %5d  %12.6f '% (row['Element'],row['Ion'],row['Wave'])
        string+='%8.6f %4d %4d %10.4e %10.4f %4d %4d ' % (row['f'],row['gl'],row['gu'], row['el'], row['eu'], row['ll'], row['ul'])
        string+='%5d %5d %10.4f %6.4f ' % (row['cl'],row['cu'],row['E_ryd'],row['gf'])
        string+='%8.4e %5d %5d %8.3f ' % (row['temp_lim'],row['ntemp'],row['type'],row['scale_par'])
        f.write('%s\n' % string)
        
        string= 'SCT   '
        for one_val in row['sct']:
            string+='%9.4e ' % one_val
        f.write('%s\n' % string)
        
        string= 'SCUPS '
        for one_val in row['scups']:
            string+='%9.4e ' % one_val
        f.write('%s\n' % string) 
    f.close()



def doit(input_file='he_1__upsilon.dat',kt_e_ratio=30,npoints=20, outroot=''):
    ytab=read_collisions(col=input_file,return_original=False)
    final=convert_ups(xtable=ytab,xratio=kt_e_ratio,npts=npoints)
    print(len(final))
    if outroot=='':
        word=input_file.split('.')
        outroot=word[0]
        outroot='%s_X' % outroot

    plot_col(final)
    plt.savefig(outroot+'.png')
    rewrite(xtable=final,outname=outroot+'.dat')

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__)
