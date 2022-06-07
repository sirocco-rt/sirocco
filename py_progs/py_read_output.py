#!/usr/bin/env python 
'''
Synopsis:
    This program enables one to read outputs from the Python radiative transfer code.
    Where possible, we use the astropy.io module to read outputs.
    There are also a number of routines for processing and reshaping various 
    data formats

    see 
    https://github.com/agnwinds/python/wiki/Useful-python-commands-for-reading-and-processing-outputs 
    for usage

Usage:
	
Arguments:
'''


# we need the classes and numpy modules 
import py_classes as cls
import numpy as np
import subprocess
import py_plot_util as util


has_astropy = True 
try:
    from astropy.io import ascii
except ImportError:
    has_astropy = False


def read_spectrum(filename):
    '''
    Load data from a spectrum output file from the radiative
    transfer code Python 

    Parameters:
        filename : file or str
        
        File, filename, or generator to read.  If the filename extension is
        ``.gz`` or ``.bz2``, the file is first decompressed. Note that
        generators should return byte strings for Python 3k.
    
    Returns         
        Success: 
        spectrum
        returns a Table of class astropy.table.table.Table

        Failure returns 1
    '''

    if not '.spec' in filename: 
        if not '.log_spec_tot' in filename:
            if not '.spec_tot' in filename:
                if not ".log_spec" in filename:
                    filename = filename + '.spec' # assume user wants the spectrum file if no suffix

    if has_astropy:
        spectrum = ascii.read(filename)
        return spectrum 

    else:
        print("Please install astropy. returning 1")
        return 1



def read_spectrum_to_class (filename, new=True):
    
    '''
    reads a Python .spec file and places in specclass array,
    which is returned

    Parameters        
        filename : file or str
            File, filename to read.  

        new:
            True means the Created column exists in the file 
    
    Returns         
        Success: 
        spectrum
        returns a spectrum class cls.specclass

        Failure returns 1
    '''
    
    if not '.spec' in filename: 
        if not '.log_spec_tot' in filename:
            if not '.spec_tot' in filename:
                filename = filename + '.spec' # assume user wants the spectrum file if no suffix

    # this deals with whether the Created column exists or not
    if new:
        add = 0
    else:
        add = 1
        
    
    # initialise the spectrum array with blank arrays
    spectrum = cls.specclass ([],[],[],[],[],[],[], [], [], []) 
    
    # first read the file into a temporary storage array
    if has_astropy:
        '''
        astropy is present, so we'll use the ascii.read 
        function to read in the file then put in classes 
        we can possibly deprecate this function but I have
        scripts which use a class format for the spectrum 
        so would like to to retain in short term 
        '''
        s = ascii.read(filename)
        spectrum.freq = s["Freq."]
        spectrum.wavelength = spectrum_temp["Wavelength"]

        if new:
            spectrum.created = s['Created']

        spectrum.emitted = s['Emitted']
        spectrum.censrc = s['CenSrc']
        spectrum.disk = s['Disk']
        spectrum.wind = spectrum_temp['Wind'] 
        spectrum.scattered = spectrum_temp['Scattered']
        spectrum.hitsurf = spectrum_temp['HitSurf']

        nangles = len(s.dtype.names) - 9 + add
        spectrum.spec = np.zeros(nangles)

        for i in range(nangles):
            spectrum.spec[i] = s[s.dtype.names[9 - add + i]]


    else:
        '''
        astropy is not present
        '''
        print("Please install astropy. returning 1")
        return 1
    


    #finally, return the spectrum class which is a series of named arrays      
    return spectrum


def read_pywind_summary(filename, return_inwind=False, mode="2d"):

    '''
    read a py_wind output file using np array reshaping and manipulation

    Parameters             
        filename : file or str
            File, filename to read, e.g. root.ne.dat  

        return_inwind: Bool
            return the array which tells you whether you
            are partly, fully or not inwind.

        mode: string 
            can be used to control different coord systems 
        
    Returns             
        d: astropy.Table.table.table object
            value is the quantity you are concerned with, e.g. ne
    '''

    if has_astropy == False:
        print("Please install astropy. returning 1")
        return 1

    
    if not ".complete" in filename:
        filename = filename + ".complete"

    # first, simply load the filename 
    #d = np.loadtxt(filename, comments="#", dtype = "float", unpack = True)
    d = ascii.read(filename)

    return d



def read_pywind(filename, return_inwind=False, mode="2d", complete=True):

    '''
    read a py_wind output file using np array reshaping and manipulation

    Parameters            
        filename : file or str
            File, filename to read, e.g. root.ne.dat  

        return_inwind: Bool
            return the array which tells you whether you
            are partly, fully or not inwind.

        mode: string 
            can be used to control different coord systems 
    
    Returns          
        x, z, value: masked arrays
            value is the quantity you are concerned with, e.g. ne
    '''

    if has_astropy == False:
        print("Please install astropy. returning 1")
        return 1


    # first, simply load the filename 
    #d = np.loadtxt(filename, comments="#", dtype = "float", unpack = True)
    d = ascii.read(filename)


    return util.wind_to_masked(d, "var", return_inwind=return_inwind, mode=mode)




def read_pf(root):

    '''
    reads a Python .pf file and returns a dictionary

    Parameters
        root : file or str
            File, filename to read.  

        new:
            True means the Created column exists in the file 
        
    Returns
        pf_dict
            Dictionary object containing parameters in pf file
    '''
    OrderedDict_present=True
    try:
        from collections import OrderedDict
    except ImportError:
        OrderedDict_present=False

    if not ".pf" in root:
        root = root + ".pf"

    params, vals = np.loadtxt(root, dtype=str, unpack=True)

    if OrderedDict_present:
        pf_dict = OrderedDict()
    else:
        pf_dict = dict()    # should work with all version of python, but bad for writing
        print("Warning, your dictionary object is not ordered.")

    old_param = None 
    old_val = None

    for i in range(len(params)):


        # convert if it is a float
        try:
            val = float(vals[i])

        except ValueError:
            val = vals[i]

        if params[i] == old_param:

            if isinstance(pf_dict[params[i]], list):
                pf_dict[params[i]].append(val)

            else:
                pf_dict[params[i]] = [old_val, val]

        else:
            pf_dict[params[i]] = val

        old_param = params[i]
        old_val = val


    return pf_dict


def write_pf(root, pf_dict):

    '''
    writes a Python .pf file from a dictionary

    Parameters        
        root : file or str
            File, filename to write.  

        pf_dict:
            dictionary to write
    
    Returns          
        pf_dict
            Dictionary object containing parameters in pf file
    '''

    if not ".pf" in root:
        root = root + ".pf"

    OrderedDict_present=True
    try:
        from collections import OrderedDict
    except ImportError:
        OrderedDict_present=False

    if (isinstance(pf_dict, OrderedDict) == False):
        print("Warning, your dictionary object is not ordered. Output file will be wrong, writing anyway.")


    f = open(root, "w")

    for key,val in pf_dict.items():

        # convert if it is a float
        if isinstance(val, list):           
            for i in range(len(val)):
                f.write("%s    %s\n" % (key, val[i]))

        #elif isinstance(val, float): 
        #    if "photons_per_cycle" not in key:
        #        if val > 1e5:
        #            f.write("%s    %e\n" % (key, val))
        else:
            f.write("%s    %s\n" % (key, val))

    f.close()

    return (0)


def setpars():
    '''
    set some standard parameters for plotting
    '''
    import matplotlib.pyplot as plt
    
    print('Setting plot parameters for matplotlib.')
    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams['axes.linewidth'] = 1.3
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams['text.usetex']='True'
    
    return 0



def read_emissivity ( root ):
    
    '''Read macro atom emissivities from a root diag file. 
       Returns two arrays, kpkt_emiss and matom_emiss.'''
    
    if not "_0.diag" in root:
        root = root + "_0.diag"
    
    # search file for smieeisivities string    
    matom_emiss, matom_abs = [], []
    kpkt_emiss = 0
    with open(root, 'r') as searchfile:
        for line in searchfile:
            
            # check if we have a matom_diagnostics line reporting level emissivities
            if 'emissivities' in line:
                
                data = line.split()
                for i in range(len(data)):
                                       
                    # we now put the appropriate values in corresponding arrays 
                    if data[i] == 'matom_abs':
                        matom_abs.append(float(data[i+1]))
                        
                    if data[i] == 'kpkt_abs':
                        kpkt_abs = float(data[i+1])
                        
                    if data[i] == 'matom_emiss':
                        matom_emiss.append(float(data[i+1]))
                        
                    if data[i] == 'kpkt_emiss':
                        kpkt_emiss = float(data[i+1])
                        

    # convert to numpy arrays
    matom_emiss = np.array(matom_emiss)
    matom_abs = np.array(matom_abs)
    
    return matom_emiss, kpkt_emiss



def thinshell_read ( root ):
    
    '''
    Read py_wind output filename for thin shell models with one cell
    '''
    
    inp = open(root, 'r')
    
    for line in inp:
        data = line.split()
        if data[0]!="#":
            if data[2]=="0":
                value = float(data[1])

    return value


def read_convergence (root ):

    ''' 
    check convergence in a diag file
    '''

    if not "_0.diag" in root:
        if not ".out" in root:
            root = root + "_0.diag"
    
    conv_fraction = []
    
    with open(root, 'r') as searchfile:
        for line in searchfile:

            # check if we have a matom_diagnostics line reporting level emissivities
            if 'Summary  convergence' in line:

                data = line.split()
                conv_fraction.append(float (data[3]))
                
    #print conv_fraction

    final_conv = conv_fraction [-1]
    return final_conv
		
