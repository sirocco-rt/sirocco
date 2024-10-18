#!/usr/bin/env python 
'''
various utilities for processing Python outputs and plotting spectra and wind properties

Usage:
    
Arguments:
'''

import py_read_output as r 
import numpy as np 
import os, sys
import time
 

standard_cmds = np.array(["1", "n","t", "r","v","1","2","3","-1",\
                         "I", "i", "1","1","1","2","0","i","0", \
                         "1","1","1","2","2","2","2","1","2","3",
                         "6", "3","6","4","6","5","0","q"])


ion_standard_variables = ["ionh1", "ionhe1", "ionhe2", "ionc4", "ionc5"]

description_dict = {
    "x": "left-hand lower cell corner x-coordinate, cm", 
    "z": "left-hand lower cell corner z-coordinate, cm", 
    "xcen": "cell centre x-coordinate, cm", 
    "zcen": "cell centre z-coordinate, cm", 
    "i": "cell index (column)", 
    "j": "cell index (row)", 
    "inwind": "is the cell in wind (0), partially in wind (1) or out of wind (<0)", 
    "converge": "how many convergence criteria is the cell failing?", 
    "v_x": "x-velocity, cm/s", 
    "v_y": "y-velocity, cm/s", 
    "v_z": "z-velocity, cm/s",  
    "vol": "volume in cm^3", 
    "rho": "density in g/cm^3", 
    "ne": "electron density in cm^-3", 
    "t_e": "electron temperature in K", 
    "t_r": "radiation temperature in K",  
    "h1": "H1 ion fraction", 
    "he2": "He2 ion fraction", 
    "c4": "C4 ion fraction",  
    "n5": "N5 ion fraction", 
    "o6": "O6 ion fraction", 
    "dmo_dt_x": "momentum rate, x-direction", 
    "dmo_dt_y": "momentum rate, y-direction", 
    "dmo_dt_z": "momentum rate, z-direction", 
    "ip": "U ionization parameter", 
    "xi": "xi ionization parameter", 
    "ntot": "total photons passing through cell", 
    "nrad": "total wind photons produced in cell", 
    "nioniz": "total ionizing photons passing through cell"
}

def get_pywind_summary (fname, vers="", den_or_frac=0):

    '''
    run version vers of py_wind on file fname.wind_save
    and generate the complete wind summary as output

    produce the output fname.complete to read

    if den_or_frac is 1, return fractions, otherwise densities
    '''

    cmds = ["1", "1", den_or_frac, "q"] # these commands create onefile summary

    isys = run_py_wind(fname, vers=vers, cmds=cmds)

    return isys



def run_py_wind (fname, vers="", cmds=None, ilv=None, py_wind_cmd = "py_wind", return_output = False):
    '''
    run version vers of py_wind on file fname.wind_save
    '''

    if cmds == None:
        cmds = standard_cmds

    x = cmds
    np.savetxt("_tempcmd.txt", x, fmt = "%s")
    logfilename = "tempfile"

    print ("Running py_wind...")
    print ("commands = {}".format(cmds))
    cmd_to_run = "{}{} {} < _tempcmd.txt > {}".format(py_wind_cmd, vers, fname, logfilename)
    isys = os.system(cmd_to_run)
    time.sleep(3)

    # return log file to use if required
    if return_output:
        logfile = open(logfilename)
        logfile_contents = logfile.read()
        

    # remove temporary file
    os.system("rm -f _tempcmd.txt {}".format(logfilename))
    if return_output:
        return (isys, logfile_contents)
    else:
        return isys



def read_pywind_smart(filename, return_inwind=False):
    '''
    read a py_wind file using np array reshaping and manipulation

    DEPRECATED
    '''

    # first, simply load the filename 
    d = np.loadtxt(filename, comments="#", dtype = "float", unpack = True)

    # our indicies are already stored in the file- we will reshape them in a sec
    zindices = d[-1]
    xindices = d[-2]

    # we get the grid size by finding the maximum in the indicies list 99 => 100 size grid
    zshape = int(np.max(zindices) + 1)
    xshape = int(np.max(zindices) + 1)


    # reshape our indices arrays
    xindices = xindices.reshape(xshape, zshape)
    zindices = zindices.reshape(xshape, zshape)

    # now reshape our x,z and value arrays
    x = d[0].reshape(xshape, zshape)
    z = d[1].reshape(xshape, zshape)

    values = d[2].reshape(xshape, zshape)

    # these are the values of inwind SIROCCO spits out
    inwind = d[3].reshape(xshape, zshape)

    # create an inwind boolean to use to create mask
    inwind_bool = (inwind >= 0)
    mask = (inwind < 0)

    # finally we have our mask, so create the masked array
    masked_values = np.ma.masked_where ( mask, values )

    #print xshape, zshape, masked_values.shape

    #return the transpose for contour plots.
    if return_inwind:
        return x, z, masked_values.T, inwind_bool.T
    else:
        return x, z, masked_values.T



def wind_to_masked(d, value_string, return_inwind=False, mode="2d", ignore_partial = True):

    '''
    turn a table, one of whose colnames is value_string,
    into a masked array based on values of inwind 

    Parameters:
        d: astropy.table.table.Table object 
            data, probably read from .complete wind data 

        value_string: str 
            the variable you want in the array, e.g. "ne"

        return_inwind: Bool
            return the array which tells you whether you
            are partly, fully or not inwind.
    Returns:
        x, z, value: Floats 
            value is the quantity you are concerned with, e.g. ne
    '''
    # this tuple helpd us decide whether partial cells are in or out of the wind
    if ignore_partial:
        inwind_crit = (0,1)
    else:
        inwind_crit = (0,2)

    if mode == "1d":
        inwind = d["inwind"]
        x = d["r"]
        values = d[value_string]

        # create an inwind boolean to use to create mask
        inwind_bool = (inwind >= inwind_crit[0]) * (inwind < inwind_crit[1])
        mask = ~inwind_bool

    # finally we have our mask, so create the masked array
        masked_values = np.ma.masked_where ( mask, values )

    #return the arrays later, z is None for 1d
        z = None


    elif mode == "2d":
        # our indicies are already stored in the file- we will reshape them in a sec
        zindices = d["j"]
        xindices = d["i"]

        # we get the grid size by finding the maximum in the indicies list 99 => 100 size grid
        zshape = int(np.max(zindices) + 1)
        xshape = int(np.max(xindices) + 1)

        # now reshape our x,z and value arrays
        x = d["x"].reshape(xshape, zshape)
        z = d["z"].reshape(xshape, zshape)

        values = d[value_string].reshape(xshape, zshape)

        # these are the values of inwind SIROCCO spits out
        inwind = d["inwind"].reshape(xshape, zshape)

        # create an inwind boolean to use to create mask
        inwind_bool = (inwind >= inwind_crit[0]) * (inwind < inwind_crit[1])
        mask = ~inwind_bool

        # finally we have our mask, so create the masked array
        masked_values = np.ma.masked_where ( mask, values )


    else:
        print ("Error: mode {} not understood!".format(mode))

    #return the transpose for contour plots.
    if return_inwind:
        return x, z, masked_values, inwind_bool
    else:
        return x, z, masked_values


    





def smooth(x,window_len=20,window='hanning'):

    '''smooth data x by a factor with window of length window_len'''

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]


    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:  
        w = eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')

    return y[window_len:-window_len+1]


def parse_rcparams(fname = "params.rc"):

    '''
    parse the file params.rc and set values in matplotlib.rcparams

    file should be of format 

        font.family             :   serif
        mathtext.fontset        :   custom

    '''

    import matplotlib as mpl

    f = open(fname, "r")

    for line in f:
        data = line.split()


        if len(data) > 0:
            if data[0] != "#":  # comments

                if data[1] != ":":
                    print("parse_rcparams: warning: unexpected format for filename %s" % (fname))


                mpl.rcParams[data[0]] = data[2]


    return 0


def get_windsave_descriptions(data, descr_dict=description_dict):

    for name in data.colnames:
        try:
            print ("{:10s} --  {}".format(name, descr_dict[name]))
        except KeyError:
            print ("no description for column {}".format(name))

    return (descr_dict)

def get_flux_at_wavelength(lambda_array, flux_array, w):

    '''
    Find the flux at wavelength w

    Parameters:
        lambda_array: array-like    
            array of wavelengths in angstroms. 1d 

        flux_array: array-like 
            array of fluxes same shape as lambda_array 

        w: float 
            wavelength in angstroms to find
    
    Returns:
        f: float 
            flux at point w
    '''

    i = np.abs(lambda_array - w).argmin()

    return flux_array[i]



