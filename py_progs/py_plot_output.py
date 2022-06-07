#!/usr/bin/env python 
'''
Synopsis:
    various plotting routines for making standard plots
    from Python outputs 

Usage:
    Either import as a module in a python session e.g.
    import py_plot_output as p 

    or run from the command line e.g. 

    py_plot_output root mode [roots to compare]


    
Arguments:
    root 
        root filename to analyse

    mode 
        mode of plotting 
        wind        plot of common wind quantites
        ions        plot of common ions 
        spec        spectrum for different viewing angles
        spec_comps  components contributing to total spectrum e.g. disk, wind
        compare  compare the root to other roots to compare
        all         make all the above plots
'''
import pylab as p
import py_read_output as r 
import numpy as np 
import os, sys
import py_plot_util as util

has_astropy = True 
try:
    import astropy
except ImportError:
    has_astropy = False

use_pretty = True
try: 
    import brewer2mpl
except ImportError:
    use_pretty = False

def make_spec_plot(s, fname, smooth_factor = 10, angles = True, components = False, with_composite=False):

    '''
    make a spectrum plot from astropy.table.table.Table object. Saves output as "spectrum_%s.png" % (fname)

    Parameters:
        s: astropy.table.table.Table
            table containing spectrum data outputted from Python 

        fname: str 
            filename to save as e.g. sv

        smooth_factor: int 
            factor you would like to smooth by, default 10

        angles: Bool 
            Would you like to plot the viewing angle spectra? 

        components: Bool 
            would you like to plot the individual components e.g. Disk Wind 
    
    Returns:
        Success returns 0
        Failure returns 1
    '''

    if type(s) != astropy.table.table.Table:
        raise TypeError("make_spec_plot takes astropy.table.table.Table object as first arg")
        return 1

    if s.colnames[8] != "Scattered":
        print ("Warning- colnames are not in expected order! {} != Scattered".format(s.colnames[8]))

    ncomponents = 9


    if angles:

        # first make viewing angle plot
        p.figure(figsize=(8,12))

        nspecs = len(s.dtype.names) - ncomponents

        nx = 1
        ny = nspecs
        if nspecs > 4:
            nx = 2
            ny = (1 + nspecs) / nx 


        print ("Making a {} by {} plot, {} spectra".format(nx, ny, nspecs))

        if with_composite:
            lambda_composite, f_composite, errors = np.loadtxt("%s/examples/telfer_qso_composite_hst.asc" % (os.environ["PYTHON"]), unpack=True, comments="#")

        for i in range(nspecs):

            p.subplot(ny, nx, i+1)

            if with_composite:
                f_1100 = util.get_flux_at_wavelength(s["Lambda"],s[s.dtype.names[ncomponents + i]], 1100.0)
                p.plot(s["Lambda"], util.smooth(s[s.dtype.names[ncomponents + i]]/f_1100, window_len = smooth_factor), label="Model")
                p.plot(lambda_composite, f_composite, label="HST composite")
                p.legend()

            else:
                p.plot(s["Lambda"], util.smooth(s[s.dtype.names[ncomponents + i]], window_len = smooth_factor))


            p.title(s.dtype.names[ncomponents + i])
            p.xlabel("Wavelength")
            p.ylabel("Flux")

        p.savefig("spectrum_%s.png" % (fname), dpi=300)
        p.clf()

    if components:

        p.figure(figsize=(8,12))
        p.subplot(211)
        p.plot(s["Lambda"], util.smooth(s["Created"], window_len = smooth_factor), label="Created")
        p.plot(s["Lambda"], util.smooth(s["Emitted"], window_len = smooth_factor), label="Emitted")
        p.legend()

        p.subplot(212)

        for i in range(4,9):
            p.plot(s["Lambda"], util.smooth(s[s.dtype.names[i]], window_len = smooth_factor), label=s.dtype.names[i])

        p.xlabel("Wavelength")
        p.ylabel("Flux")
        p.legend()

        p.savefig("spec_components_%s.png" % (fname), dpi=300)
        p.clf()
        
    return 0


def make_spec_plot_from_class(s, fname, smooth_factor = 10, angles = True, components = False):

    '''
    make a spectrum plot from py_classes.specclass object. Saves output as "spectrum_%s.png" % (fname)

    Parameters:
        s: specclass object
            table containing spectrum data outputted from Python 

        fname: str 
            filename to save as e.g. sv

        smooth_factor: int 
            factor you would like to smooth by, default 10

        angles: Bool 
            Would you like to plot the viewing angle spectra? 

        components: Bool 
            would you like to plot the individual components e.g. Disk Wind 
    
    Returns:
        Success returns 0
        Failure returns 1
    '''


    if angles:

        # first make viewing angle plot
        p.figure(figsize=(8,12))
        nspecs = len(s.spec)

        nx = 1
        ny = nspecs
        if nspecs > 4:
            nx = 2
            ny = (1 + nspecs) / nx 


        print ("Making a {} by {} plot, {} spectra".format(nx, ny, nspecs))

        for i in range(nspecs):


            p.subplot(ny, nx, i+1)

            p.plot(s.wavelength, util.smooth(s.spec[i], window_len = smooth_factor))

            p.xlabel("Wavelength")
            p.ylabel("Flux")

        p.savefig("spectrum_%s.png" % (fname))
        p.clf()

    if components:

        p.figure(figsize=(8,12))
        p.subplot(211)
        p.plot(s.wavelength, util.smooth(s.created, window_len = smooth_factor), label="Created")
        p.plot(s.wavelength, util.smooth(s.emitted, window_len = smooth_factor), label="Emitted")

        p.subplot(212)
        p.plot(s.wavelength, util.smooth(s.censrc, window_len = smooth_factor), label="CenSrc")
        p.plot(s.wavelength, util.smooth(s.disk, window_len = smooth_factor), label="Disk")
        p.plot(s.wavelength, util.smooth(s.wind, window_len = smooth_factor), label="Wind")
        p.plot(s.wavelength, util.smooth(s.hitsurf, window_len = smooth_factor), label="HitSurf")
        p.plot(s.wavelength, util.smooth(s.scattered, window_len = smooth_factor), label="Scattered")
        p.xlabel("Wavelength")
        p.ylabel("Flux")
        p.legend()

        p.savefig("spec_components_%s.png" % (fname))
        p.clf()

    return 0




def make_wind_plot(d, fname, var=None, shape=(4,2), axes="log", den_or_frac=0, fname_prefix="wind", lims=None):
    '''
    make a wind plot from astropy.table.table.Table object. Saves output as "spectrum_%s.png" % (fname)

    Parameters:
        d: astropy.table.table.Table
            table containing wind data outputted from Python 
            if == None then this routine will get the data for you

        fname: str 
            filename to save as e.g. sv

        var: array type 
            array of string colnames to plot 

        angles: Bool 
            Would you like to plot the viewing angle spectra? 

        components: Bool 
            would you like to plot the individual components e.g. Disk Wind 

        axes: str 
            lin or log axes

        den_or_frac: int
            0 calculate ion densities
            1 calculate ion fractions

        lims: array-like
            limits of plot, specified as ((xmin,xmax), (ymin, tmax))
            can be array or tuple. Default is Nonetype.

    Returns:
        Success returns 0
        Failure returns 1
    '''
    
    if d == None:
        util.get_pywind_summary(fname, den_or_frac=den_or_frac)
        d = r.read_pywind_summary(fname)

    if var == None:
        var = ["ne", "te", "tr", "IP", "nphot", "v", "w", "ionC4"]

    if axes != "lin" and axes != "log":
        print ("Error: didn't understand mode {}, defaulting to log".format(axes))
        axes = "log"

    nplots = len(var)

    # check shape is ok for variables required
    if shape[0] * shape[1] < nplots:
        print ("Error: shape is less than length of var array")
        return 1

    if shape[0] * shape[1] > nplots:
        print ("Warning: shape is more than length of var array")

    p.figure(figsize=(8,12))

    # cycle over variables and make plot
    for i in range(nplots):

        p.subplot(shape[0], shape[1], i+1)

        value_string = var[i]

        x,z,v = util.wind_to_masked(d, value_string)

        print (np.mean(v))

        if "ion" in value_string and den_or_frac==1:
            p.pcolormesh(x,z,np.log10(v), vmin=-5,vmax=0.1)
        else:
            p.pcolormesh(x,z,np.log10(v))
        p.colorbar()
        p.title("Log(%s)" % value_string)

        # if lims != None:
        #     p.xlim(lims[0][0], lims[0][1])
        #     p.ylim(lims[1][0], lims[1][1])

        if axes == "log":
            # log axes
            p.loglog()
            p.loglog()

        # Figure out the "best" values for the lower and upper limits
        lower_lim = x[0][0] if x[0][0] != 0 else x[1][0]
        upper_lim = x[-1][0]

        p.xlim(lower_lim, upper_lim)
        p.ylim(lower_lim, upper_lim)

    p.tight_layout()
    p.savefig("%s_%s.png" % (fname_prefix, fname))

    return 0


def make_spec_comparison_plot (s_array, labels, fname="comparison", smooth_factor = 10, angles = True, components = False):

    '''
    make a spectrum comparison plot from array of astropy.table.table.Table objects. Saves output as "spectrum_%s.png" % (fname)

    Parameters:
        s_array: array-like of astropy.table.table.Table objects
            table containing spectrum data outputted from Python 

        labels: array-like
            strings of labels for each spectrum 

        fname: str
            filename to save as e.g. sv

        smooth_factor: int 
            factor you would like to smooth by, default 10

        angles: Bool 
            Would you like to plot the viewing angle spectra? 

        components: Bool 
            would you like to plot the individual components e.g. Disk Wind 
    
    Returns:
        Success returns 0
        Failure returns 1
    '''

    ncomponents = 9



    if angles:

        # first make viewing angle plot
        p.figure(figsize=(8,12))


        nspecs = len(s_array[0].dtype.names) - ncomponents
        nx = 1
        ny = nspecs
        if nspecs > 4:
            nx = 2
            ny = (1 + nspecs) / nx 


        print ("Making a {} by {} comparison plot, {} spectra".format(nx, ny, nspecs))

        for j in range(len(s_array)):
            for i in range(nspecs):

                p.subplot(ny, nx, i+1)

                p.plot(s_array[j]["Lambda"], util.smooth(s_array[j][s_array[j].dtype.names[ncomponents + i]], window_len = smooth_factor), label=labels[j])

                p.xlabel("Wavelength")
                p.ylabel("Flux")


                if i == 0 and j == (len(s_array) - 1):
                    p.legend()

        p.savefig("spectrum_%s.png" % (fname), dpi=300)
        p.clf()

    if components:

        p.figure(figsize=(8,12))

        n = len(s_array)

        for j in range(n):

            p.subplot(j+1,1,1)

            s = s_array[j]

            p.plot(s["Lambda"], util.smooth(s["Created"], window_len = smooth_factor), label="Created")
            p.plot(s["Lambda"], util.smooth(s["Emitted"], window_len = smooth_factor), label="Emitted")

            for i in range(4,9):
                p.plot(s["Lambda"], util.smooth(s[s.dtype.names[i]], window_len = smooth_factor), label=s.dtype.names[i])

            p.xlabel("Wavelength")
            p.ylabel("Flux")
            p.legend()

        p.savefig("spec_components_%s.png" % (fname), dpi=300)
        p.clf()

    return 0




# Next lines permit one to run the routine from the command line with various options -- see docstring
if __name__ == "__main__":
    import sys

    ishow='yes'
    if len(sys.argv)>2:
        fname = sys.argv[1]
        mode = sys.argv[2]
    else:
        print (__doc__)
        sys.exit(1)

    # try to read a parms file in the directory
    io_print = True
    try:
        util.parse_rcparams()
    except IOError:
        print ("Tried to read parameters from params.rc in local directory, but none found- Continuing.")
        io_print = False

    if io_print:
        print ("Read parameters from params.rc")

    if mode == "spec":
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname)

    elif mode == "specc":   # compare to the HST composite QSO spectrum
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname, with_composite=True)

    elif mode == "wind":
        make_wind_plot(None, fname)

    elif mode == "spec_comps":
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname, angles = False, components = True)

    elif mode == "ions":
        make_wind_plot(None, fname, 
                       var = ["ionH1", "ionH2", "ionC3", "ionC4", "ionC5", "ionSi4", "ionN5", "ionO6"],
                       fname_prefix="ions", den_or_frac = 1)

    elif mode == "compare":     # comapre 2 or more spectra
        if len(sys.argv) <= 3:
            print (__doc__)
            sys.exit(1)

        s_array = [r.read_spectrum(sys.argv[1])]
        labels=[fname]
        for i in range(3,len(sys.argv)):
            s = r.read_spectrum(sys.argv[i])
            labels.append(sys.argv[i])
            s_array.append(s)

        make_spec_comparison_plot(s_array, labels)

    elif mode == "all":
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname, components = True)
        make_wind_plot(None, fname)
        make_wind_plot(None, fname, 
                       var = ["ionH1", "ionH2", "ionC3", "ionC4", "ionC5", "ionSi4", "ionN5", "ionO6"],
                       fname_prefix="ions", den_or_frac = 1)

    else:
        print ("didn't understand mode {}".format(mode) )
        print (__doc__)










