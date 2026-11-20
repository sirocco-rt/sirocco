#!/usr/bin/env python
# coding: utf-8
'''
Create a summary of the differences between various .pf files in a directory


Command line usage (if any):

    usage: ModSum.py 

Description:  

    This routine simply looks at all of the .pf files in 
    the current working directory and constructs a summary
    of the runs that are there.  

    It prodcues a file, hardwired to be ModSum.txt that
    is an astropy table that indicates whether the variious
    .pf files have been done.  The remaining columnts are
    show what is different about the various models

Primary routines:

    doit

Notes:
                                       
History:

220812 ksl Converted from and ipynb file to an executable

'''

import sys
import os
from astropy.io import ascii
import numpy as np
from glob import glob
from astropy.table import Table,join,vstack




def get_status(run='big'):
    if run.count('.pf'):
        run=run.replace('.pf','')
    if run.count('.sig')==0:
        run=run+'.sig'
    if os.path.exists(run)==False:
        return 0,'NotRun'
    f=open(run)
    lines=f.readlines()
    last_line=lines[-1]
    word=last_line.split()
    time=eval(word[5])
    status=word[6]
    return time, status



def get_models(xlim=''):
    if xlim=='':
        files=glob('*.pf')
    elif xlim.count('pf'):
        files=glob(xlim)
    else:
        xlim=xlim+'.pf'
        files=glob(xlim)
    
    xfiles=[]
    for file in files:
        if file.count('out'):
            pass
        else:
            xfiles.append(file)
    xtimes=[]
    xstatus=[]
    for one in xfiles:
        time,status=get_status(one)
        xtimes.append(time)
        xstatus.append(status)
    

    xtab=Table([xfiles,xstatus,xtimes],names=['Root','Status','ExecTime'])
    xtab.sort('Root')
        
    return xtab



def get_pf(pf_file='big'):
    if pf_file.count('.pf')==0:
        pf_file=pf_file+'.pf'
    try:
        f=open(pf_file)
    except:
        print('get_pf: %s does not appear to exist' % pf_file)
        raise IOError
    xkeys=[]
    xvals=[]
    lines=f.readlines()
    for one in lines:
        one=one.strip()
        word=one.split()
        if len(word)==0 or word[0].count('#'):
            pass
        else:
            key=word[0]
            j=key.find('(')
            key=key[0:j]
            value=word[1]
            # print(key,value)
            xkeys.append(key)
            xvals.append(value)
    xkeys=np.array(xkeys)
    xvals=np.array(xvals)
    names,counts=np.unique(xkeys,return_counts=True)
    i=0
    while i<len(names):
        if counts[i]>1:
            for one in xkeys:
                j=1
                if one==names[i]:
                    one='%s_%d' % (one,j)
                    j+=1
        i+=1
    xtab=Table([[pf_file]],names=['Root'])
    i=0
    while i <len(xkeys):
        # print(xkeys[i],xvals[i])
        xtab[xkeys[i]]=xvals[i]
        i+=1
    return xtab       
           




def make_master(glob_root=''):
    xmodels=get_models(glob_root)
    xdetail=[]
    for one in xmodels:
        one_detail=get_pf(one['Root'])
        if len(xdetail)==0:
            xdetail=one_detail.copy()
        else:
            xdetail=vstack([xdetail,one_detail])
    print(len(xdetail))
    # xdetail.info()
    xall=join(xmodels,xdetail,join_type='left')
    
    for one_name in xall.colnames:
        if len(np.unique(xall[one_name]))==1:
            xall.remove_column(one_name)
    return xall
        
        

def doit():
    '''
    Run the script

    Description:

    Notes:

    History:


    '''
    xall=make_master()     

    xall.write('ModSum.txt',format='ascii.fixed_width_two_line',overwrite=True)
    

    return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        print (__doc__)
    else:
        doit()
