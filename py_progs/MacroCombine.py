#!/usr/bin/env python 

'''
Combine levels, etc in a MacroAtoms (generated with MakeMacro)  
to to create a more succinct MacroAtom model

Command line usage:

    usage: MacroCombine.py [-guess] level_file line_file phot_file upsilon_file

    with
        -guess implies that the program will try to guess
            what levels should be combined.  In the absence
            of this, the program will check if the 
            levels one wants to combine have beeen
            set already.  

    The order of the input files is IMPORTANT!

Description:  

    This routine attempts to reconfigure the level, line
    and phot and collision files  to allow 
    for the fact that one would often like
    to combine levels to make MacroAtoms smaller

    There are two basically possibilities:

    With the -guess option, the program takes the original level file, in 
    the format produced by MakeMacro, and guesses what levels to combine, 
    producing an 'intermediate level file' that contains information 
    extra columns indicating what to combine.

    It then continues to produce final level fileline and
    photoionization files that reflect this choice.

    Alternatively, without the -guess option the routine will read
    epects to read a level file which contains the extra columns
    that indcicate what needs to combine,  and uses this to produce
    final level, line and photoionization files.

    Generally, speaking one should start with the -guess option
    and then modify the xlev column to reflect a final choice
    of levels to combine, and then one should rerun the routine
    without the -guess option to synchronize the rest of the files.

Primary routines:


Notes:
                                       
History:

231013 ksl Coding begun
240910 ksl Updated to add collisions to the process and to get the
            correct format for the output photometry files

'''

import sys
from astropy.io import ascii
from astropy.table import join, Table
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from RedoUpsilon import convert_ups, rewrite


def read_phot(file='h_1_phot.dat'):
    '''
    Read a Sirocco photoionization file and return 
    the information in a form the is eacy to manipulate
    
    The routine returns a table that summarizes all of the
    different x-sections, and a list of tables, one table for
    each x-section.
    '''
    print('Beginning read_phot ', file)
    f=open(file)
    lines=f.readlines()
    i=0
    summary=[]
    xsections=[] # This is going to have one table for each set of xsections
    i=0
    elem=[]
    ion=[]
    level=[]
    level_up=[]
    ethresh=[]
    nlines=[]
    line_no=[]
    
    for line in lines:
        word=line.split()
        if word[0][0]=='#':
            pass
        elif word[0].count('PhotMacS'):
            elem.append(int(word[1]))
            ion.append(int(word[2]))
            level.append(int(word[3]))
            level_up.append(int(word[4]))
            ethresh.append(float(word[5]))
            nlines.append(int(word[6]))
            line_no.append(i)
            # one_record=[elem,ion,level,level_up,ethresh,nlines,line_no]
            # summary.append(one_record)
            if len(line_no)>1:
                # The we have already accumulated xsections
                x=Table(np.array(one_xsection),names=['e','sigma'])
                xsections.append(x)
            one_xsection=[]
        else:
            one_xsection.append([eval(word[1]),eval(word[2])])
        i+=1
    x=Table(np.array(one_xsection),names=['e','sigma'])
    xsections.append(x)
    
    tab_sum=Table([elem,ion,level,level_up,ethresh,nlines,line_no],names=['z','ion','ll','ul','ethresh','nlines','line_no'])

    ncross=range(len(tab_sum))
    tab_sum['Number']=ncross
    
    # Getting the right format would need to be done earlier, because summary contins floats
    tab_sum.write('phot_sum.txt',format='ascii.fixed_width_two_line',overwrite=True)

    print('Finished  read_phot %s %d %d' % (file,len(tab_sum),len(xsections)))
    
    return tab_sum,xsections



def write_phot(outfile,xtab,xsec):
    '''
    Write a valid Macro photometry xsection file from what amounts
    to an internal format for this routine.  The
    internal form consists of 

    * a table (xtab) which summarizes each x-section and 
    * a list (xsec) that contains tables of the xsections
    for each x-section
    '''

    if len(xtab)!=len(xsec):
        print('Warning: Writing a phot file in which only %d of the %d levels have xsections'  % (len(xsec),len(xtab)))
    else:
        print('All levels have photx-sections')
    
    print(xtab)
    f=open(outfile,'w')
    i=0
    qtab=xtab[xtab['xsection']=='yes']
    while i<len(qtab):
        one=qtab[i]
        string='PhotMacS  %3d %3d %3d %3d %10.6f %3d' %  (one['z'],one['ion'],one['ll'],one['ul'],one['ethresh'],one['nlines'])
        # print(string)
        f.write('%s\n' %string)
        one_x=xsec[i]
        # print('foo \n',one_x)
        j=0
        while j<len(one_x):
            f.write('PhotMac %10.6f %10.6e\n' % (one_x['e'][j],one_x['sigma'][j]))
            j+=1
        i+=1
    


def plot_xsec(selection,xtab,style='-',clear=True):
    '''
    Here 

    * xtab is a list of tables containing the xsections and 
    * selection is a list of the x-sections to plot

    Note: section is 0-based and does not necessarily
    correspond to a specific level

    Aslo not that if one wants to save the plot, one 
    must currently do that separately
    '''
    isel=np.array(selection)
    plt.figure(1,(6,6))
    if clear:
        plt.clf()
    for one in isel:
        xcross=xtab[one]
        plt.semilogy(xcross['e'],xcross['sigma'],style,label='%d' % one, alpha=0.5)
    
    plt.xlabel('E (ev)')
    plt.ylabel(r'$\sigma$ (cm$^{-2}$)')
    plt.legend()
    plt.tight_layout()
    return



def reweight(levels,g,xtab):
    '''
    Reweight a group of of photionsation xsections, based on g

    xtab contains the entire list of xsections, each
    element of which corresponds to the x=sections
    for a single transition.

    level contains a list, e.g [1,2,3]  of the xsections to consider

    and g is a listi, e.g [2,2,4] containign the g factors associated
    with the various x-sections.

    Note: As written this is zero based and this may
    not be what we want.

    '''
    # first check that the energies are the same
    
    xsig=[]
    if len(levels)<1:
        print('Houston we have a problem: the number of levels is ',len(levels),len(g),len(xtab))
    else:
        print('Houston we should be OK: the number of levels is ',len(levels),len(g),len(xtab))

    i=0
    while i<len(g):
        # for one in levels:
        # print('test:',one)
        # print('test\n',xtab[one])
        if i==0:
            xenergy=xtab[i]['e']
            xsig=g[i]*np.array(xtab[i]['sigma'])
        else:
            energy=xtab[i]['e']
            sigma=xtab[i]['sigma']
            interpolation_func = interp1d(energy, sigma, kind='linear', bounds_error=False, fill_value='extrapolate')
            xxsigma=interpolation_func(xenergy)
            # print(i,len(xsig),len(xtab[i]['sigma']),len(xxsigma))
            xsig+=g[i]*xxsigma
            # print(len(xsig))
        i+=1
    xsig/=np.sum(g)
    e=np.array(xtab[0]['e'])
    qtab=Table([e,xsig],names=['e','sigma'])
    qtab['e'].format='.6f'
    qtab['sigma'].format='.7e'
    # print(qtab)
    return qtab


def redo_phot(xlevel_file='h_1_lev2phot_guess.txt',phot_file_orig='h_1_phot.dat',root_out='foo'):
    '''
    Read a file containing the levels we want to compress amd  the associted phot_file
    and produce a new phot file
    
    The photometry file reads data into a summary table (tab_sum) and a list of tables with the photoionization data (xsections)
    We need to construct a new summary tabel (tab_sum_out) and a new lst of talbes (xsections)

    This routine alos writes out a new configuration file.
    '''

    if root_out=='':
        word=phot_file_orig.split('_')
        if len(word)>2:
            root_out='%s_%s' % (word[0],word[1])
        else:
            word=xlevel_file.split('_')
            if len(word)>2:
                root_out='%s_%s' % (word[0],word[1])
            else:
                root_out=foo
                print('No root out provided and no obvious split for inputs files so using %s' % root_out)
         
    
    try:
        xlevel=ascii.read(xlevel_file)
    except:
        print('Could not read %s' % level_file)
        return
   
    tab_sum,xsections=read_phot(phot_file_orig)
    
    print('The original num of levels  is ',len(xlevel))
    print('The number of xsections     is ',len(tab_sum))
    
    # Strip off the next ion if it is there
    xlevel=xlevel[xlevel['config']!='Next']
    nlev=xlevel['xlev'][-1]

    print('There will be  %d new levels' % (nlev))



    
    #create a place to put a revised level_file
    lev_out=xlevel[0:nlev]
    tab_sum_out=tab_sum[0:nlev]
    xsections_out=[]

    tab_sum.rename_column('ll','lvl')

    # print('Tell me\n',tab_sum_out.info)
    j=0
    i=1
    # i refers to the output level

    good=[]
    while i<=nlev:
        print('Start processing ouput level  %d ' % i)
        # cur_leve is a table containing all the levels to map to i
        cur_lev=xlevel[xlevel['xlev']==i]
        xcur=join(cur_lev,tab_sum,keys='lvl',join_type='left')

        zcur=xcur[xcur['Number']>=0]

        # zcur contains the data one wants to combine

        if hasattr(xcur['Number'],'mask'):
            xmask=xcur['Number'].mask
            zcur=xcur[~xmask]
            print('Levels to combine %d out of %d originally' % (len(zcur),len(xcur)))
            print('Not all original levels had  x-sections')

            # print(xcur[xmask])
        else:
            print('All %d levels had x-sections' % (len(xcur)))
            zcur=xcur
       


        if len(zcur)==0:
            print('Level %d had no x-sections' % (i))
            good.append('no')
        if len(zcur)==1:
            fin_sections=xsections[zcur[0]['Number']]
            xsections_out.append(fin_sections)
            tab_sum_out['ethresh'][i-1]=zcur['ethresh'][0]
            tab_sum_out['nlines'][i-1]=len(fin_sections)
            good.append('yes')
        elif len(zcur)>1:
            xlev=np.array(zcur['lvl'])
            glev=np.array(zcur['g'])
            n=np.array(zcur['Number'])
            qsections=[]
            for one in n:
                qsections.append(xsections[one])
            # qsections=xsections[n]

            xlev=xlev-1

            vsections=reweight(xlev,glev,qsections)
            xsections_out.append(vsections)
            tab_sum_out['ethresh'][i-1]=np.min(zcur['ethresh'])
            tab_sum_out['nlines'][i-1]=len(vsections)
            good.append('yes')
                          
        i+=1
    
    # print('xfffff',i,len(xsections_out))
    lev_out['g']=lev_out['G']

    # print('nada\n',tab_sum_out.info())
    i=0
    while i<nlev:

        lev_out['lvl'][i]=i+1
        tab_sum_out['ll'][i]=i+1
        zz=xlevel[xlevel['xlev']==i+1]
        lev_out['ion_pot'][i]=np.average(zz['ion_pot'])
        lev_out['ex'][i]=np.average(zz['ex'])
        if len(zz)>1:
            lev_out['config'][i]='combined'
        else:
            lev_out['config'][i]=xlevel['config'][0]
        lev_out['rad_rate'][i]=np.max(zz['rad_rate'])
        lev_out['islp'][i]=np.median(zz['islp'])
        i+=1        
        
    
    phot_file='%s_phot_final.dat' %   root_out
    lev_file='%s_lev_final.dat' % root_out
    tab_sum_out['xsection']=good
    write_phot(phot_file,tab_sum_out,xsections_out)
    lev_out['ion_pot'].format='.4f'
    lev_out['ex'].format='.4f'

    lev_out.write(lev_file,format='ascii.fixed_width_two_line',overwrite=True)



def add_gtot(filename='he_2_levels_comb.dat'):
    ''' 
    Add the final G factors to the for
    the combined levels to the level
    file
    
    The table that is produced still has 
    the same number of levels as the original 
    table.  It simply has an additional
    column labelled G that will be the 
    multiplicty associated with the final level
    
    '''
    # print('test',filename)
    
    x=ascii.read(filename)
    # x.info()
    x['G']=-1
    try:
        levels=np.array(np.unique(x['xlev']))
    except:
        print('Error: the levels file did not have the column xlev')
        print('This means it does not know what levels to combine')
        print('Rerun the routine with the -guess optio and then, modify as necessary')
        return ''

    # print(levels)
    gg=np.zeros(len(levels))
    xion=x['Ion'][0]  # This is to allow us to avoid the next level up
    for one in x:
        if one['Ion']==xion:
            # print('ok ',one)
            gg[one['xlev']-1]+=one['g']
    
    #print(gg)
    
    # Now we have the sums, we just need to add them to the file
    
    for one in x:
        if one['Ion']==xion:
            # print('ok ',one)
            one['G']=gg[one['xlev']-1]

    x.write(filename,format='ascii.fixed_width_two_line',overwrite=True)
    return filename
  
    


def redo_lines(master='h_1_levels_comb2.dat',lines='h_1_lines.dat'):
    '''
    Write out a new line file, based on the intermidate
    level file
    '''
    xmaster=ascii.read(master)
    xmaster.rename_column('Element','z')
    xmaster.rename_column('Ion','ion')
    xlines=ascii.read(lines)
    # xlines['ll_final']=-99
    # xlines['ul_final']=-99
    # print(len(xlines))
    
    xlow=xmaster['z','ion','lvl','xlev','G']
    xlow.rename_column('lvl','ll')
    xlow.rename_column('xlev','xll')
    xlow.rename_column('G','Gl')

    xlow.write('low.txt',format='ascii.fixed_width_two_line',overwrite=True)
    
    xup=xmaster['z','ion','lvl','xlev','G']
    xup.rename_column('lvl','ul')
    xup.rename_column('xlev','xul')
    xup.rename_column('G','Gu')
    
    xlines=join(xlines,xlow,join_type='left')
    # print(len(xlines))
        
    xlines=join(xlines,xup,join_type='left')
    
    # print(len(xlines))
    
    xlines['A']=xlines['f']*xlines['gl']/xlines['gu']
    xlines['A'].format='.6f'
    
    xlines['x']=xlines['xll']*1000+xlines['xul']
    
    
    xlines.sort('ul')
    # xlines.info()
    xlines.write('foo_lines.dat',format='ascii.fixed_width_two_line',overwrite=True)
    
    z,xindex=np.unique(np.array(xlines['x']),return_index=True)
    
    xfinal=xlines[xindex]

    
    # print(xfinal)
    
    # Now make a new table that contains only one row for each value in z

        
    for one_row in xfinal:
        all_rows=xlines[xlines['x']==one_row['x']]
        one_row['Wave']=np.average(all_rows['Wave'])
        one_row['A']=np.sum(all_rows['A']*all_rows['gu'])/one_row['Gu']
        one_row['f']=one_row['Gu']*one_row['A']/one_row['Gl']
        # one_row['f']=np.sum(all_rows['f']*all_rows['gu'])/np.sum(all_rows['gu'])
        one_row['gu']=np.sum(all_rows['gu'])
        one_row['ul']=np.sum(all_rows['ul'])
        one_row['eu']=np.average(all_rows['eu'])
        
    xfinal['f'].format=xfinal['eu'].format='.5f'
    xfinal['Wave'].format='.5f'


    words=lines.split('.')
    out_name='%s_final.%s' % (words[-2],words[-1])

    
    
    xxx=xfinal['Dtype','z','ion','Wave','f','Gl','Gu','el','eu','xll','xul']
    xxx.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    return xxx



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
    # print(len(lines))
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

    # print(len(cstren_line),len(sct_line),len(scups_line))

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


def combine_collisions(coltab):
    result=convert_ups(xtable=coltab,xratio=30,npts=20)
    result[0]['sct']=np.sum(result['sct'],axis=0)
    # print('Hello')
    # print(result.info)
    return result[0]


def redo_collisions(master='he_1_levels_guess.dat',col='he_1__upsilon.dat'):
    '''
    Write out a new line file, based on the intermidate
    level file
    '''
    print('***Processing collisions***')
    xmaster=ascii.read(master)
    # xmaster.info()
    # print(xmaster['xlev'])
    collisions,cstren,sct,scups=read_collisions(col,return_original=True)
    final_collisions=collisions.copy()
    final_collisions['sct']=np.empty(len(final_collisions),dtype=object)
    final_collisions['scups']=np.empty(len(final_collisions),dtype=object)
    
    # final_collisions['sct']=np.array(final_collisions['sct'],dtype=object)
    # final_collisions['scups']=np.array(final_collisions['scups'],dtype=object)
    xcol=collisions['Number','ll','ul']
    xmaster=xmaster['lvl','xlev']
    xmaster['ll']=xmaster['lvl']
    xmaster['ul']=xmaster['lvl']
    ll=join(xmaster,xcol,keys=['ll'],join_type='inner')
    ul=join(xmaster,xcol,keys=['ul'],join_type='inner')
    ll['xlev_lower']=ll['xlev']
    ul['xlev_upper']=ul['xlev']
    q=join(ll['Number','xlev_lower'],ul['Number','xlev_upper'],join_type='left')
    # At this point q contains what final upper level and what final lower level to associate
    # with each collision strength, and so we need to work out the unique combinations

    lval,counts=np.unique(q['xlev_lower'],return_counts=True)
    uval,ucounts=np.unique(q['xlev_upper'],return_counts=True)
    i=0
    nfinal=0
    while i <len(lval):
        xlev_lower=q[q['xlev_lower']==lval[i]]
        # print(len(xlev_lower))
        j=0
        while j<len(uval):
            xlev_upper=xlev_lower[xlev_lower['xlev_upper']==uval[j]]
            # print(lval[i],uval[j],len(xlev_upper))
            if len(xlev_upper)==0:
                print ('No x-section for  %d %d' % (lval[i],uval[j]))
            elif len(xlev_upper)==1:
                print('Copy x-section for %d %d' % (lval[i],uval[j]))
                final_collisions[nfinal]=collisions[xlev_upper['Number'][0]]
                final_collisions['ll'][nfinal]=final_collisions['cl'][nfinal]=lval[i]
                final_collisions['ul'][nfinal]=final_collisions['cu'][nfinal]=uval[j]
                nfinal+=1
            else:
                print('Combine x-section for %d %d' % (lval[i],uval[j]))
                # print(xlev_upper)
                zcollisions=collisions[xlev_upper['Number']]
                xcombine=combine_collisions(zcollisions)
                # print('testy',nfinal,len(final_collisions))
                # print(final_collisions.info)
                for one_col in final_collisions.colnames:
                    final_collisions[one_col][nfinal]=xcombine[one_col]
                final_collisions['ll'][nfinal]=final_collisions['cl'][nfinal]=lval[i]
                final_collisions['ul'][nfinal]=final_collisions['cu'][nfinal]=uval[j]
                nfinal+=1


            j+=1
        i+=1

    final_collisions=final_collisions[:nfinal]
    outfile=col.replace('.dat','_final.dat')
    rewrite(final_collisions,outfile)

    
    print('***Finished collisions***')
    return 


    
    
def xguess(tab_name='h_1_lev2phot.txt',out_name='',frac=0.01,redo=True):
    '''
    Guess how levels for an atom should be combined and write out
    the intermediate level file that contains the various suggestions

    '''
    xtab=ascii.read(tab_name)
    cols=xtab.colnames
    for one in cols:
        if one=='xlev'and redo!=True:
            print('This file %s is already has levels to combine, so returning' % tab_name )
            print('To force a repeat, set redo to True')
            return
    xtab['xlev']=-99
    
    i=1
    ilev=1
    xpot=xtab['ion_pot'][0]
    xtab['xlev'][0]=1
    while i<len(xtab):
        cpot=xtab['ion_pot'][i]
        delta=np.fabs((cpot-xpot)/cpot)
        if delta<frac:
            xtab['xlev'][i]=ilev
        else:
            ilev+=1
            xtab['xlev'][i]=ilev
            xpot=cpot
        i+=1
        
    if out_name=='':        
        words=tab_name.split('.')
        out_name='%s_guess.%s' % (words[-2],words[-1])
        
    xtab.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    
    return out_name


def steer(argv):
    '''
    This is just steering routine to allow operations from the command line
    '''

    level_file=''
    line_file=''
    phot_file=''
    col_file=''
    iguess=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        if argv[i]=='-guess':
            iguess=True
        elif argv[i][0]=='-':
            print(__doc__)
            print('Unknown option %s' % argv[i])
        elif level_file=='':
            level_file=argv[i]
        elif line_file=='':
            line_file=argv[i]
        elif phot_file=='':
            phot_file=argv[i]
        elif col_file=='':
            col_file=argv[i]
        else:
            print(__doc__)
            print ('Too many arguments' % argv)
            return

        i+=1

    if level_file=='':
        print('No level file provided')
        return
    else:
        if os.path.isfile(level_file)==False:
            print('The level file does not appear to exist: ',level_file)
            return

    if iguess==False:
        x=ascii.read(level_file)
        colnames=x.colnames
        check=False
        for one in colnames:
            if one=='xlev':
                check=True
        if check==False:
            print('Error: levels file does not have columns to show what levels to combine')
            print('Run with -guess option and then modfify the input levels file as necessary')
            return


    if iguess:
        level_file=xguess(level_file)
        add_gtot(level_file)

    if line_file!='':
        if os.path.isfile(line_file)==False:
            print('The line file does not appear to exist: ',line_file)
            return
        add_gtot(level_file)
        redo_lines(level_file,line_file)
    else:
        print('No line file provided')
        return

    if phot_file!='':
        if os.path.isfile(phot_file)==False:
            print('The phot file does ont appear to exist: ',phot_file)
            return
        redo_phot(xlevel_file=level_file,phot_file_orig=phot_file,root_out='')
    else:
        print('No phot file provided')

    if col_file!='':
        if os.path.isfile(col_file)==False:
            print('The line file does not appear to exist: ',col_file)
            return
        add_gtot(level_file)
        redo_collisions(level_file,col_file)
    else:
        print('No col file provided')
        return


                



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)

