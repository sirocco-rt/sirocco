Spectra Files
#############

SIROCCO is intended to produce simulated spectra.  These spectra are all ascii tables intended to be accessible with software packages such as astropy.

All of the ascii begin with commented headers that contain all of the parameters of associated with a run,
along with the date of the run and the specific version of SIROCCO used to make the run.
In principle, if one still has access to any of the spectra, one can reproduce the entire run again.

Broad band spectra are created from the last ionization cycle. (More accurately the broad band spectra are written out at the end of each ionization cycle, so one the program is finished one has the 
broad band spectrum of the last cycle)  

Detailed are calculated from all of the spectral cycles. (Properly normalized spectra are written out at the 
end of each spectral cycle, and with each cycle the photon statistics improves.)

The units in which the spectra are written out is also indicated in the header to the file.  


For a model with root name *cv*, the following broadband spectra will be created:

* **cv.log_spec_tot**
* **cv.spec_tot_wind**
* **cv.log_spec_tot_wind**
* **cv.spec**
* **cv.log_spec**

File types
==========

.log_spec_tot
  An ascii file that contains various spectra from the ionization-calculation phase of the program on a log frequency scale.
  The first few lines of the file (omitting the header) are as follows:

  .. code::

    Freq.        Lambda     Created    WCreated   Emitted    CenSrc     Disk       Wind       HitSurf
    2.524334e+14 11876.102  3.5244e+18          0 3.5244e+18          0 3.5244e+18          0 1.1547e+16
    2.550397e+14 11754.737  3.4721e+18          0 3.4721e+18          0 3.4721e+18          0 2.8761e+15
    2.576461e+14 11635.827  3.4433e+18          0 3.4433e+18          0 3.4433e+18          0 2.8835e+15
    2.602524e+14 11519.299  3.6858e+18          0 3.6858e+18          0 3.6858e+18          0 2.8706e+15
    2.628587e+14 11405.082  3.6711e+18          0 3.6711e+18          0 3.6711e+18          0 1.1528e+16


The first column describes the frequency of a particular bin in :math:`Hz`. The frequency stated is placed on the left edge of the bin, i.e the smallest (minimum) frequency value within that bin range. 
The second column is the wavelength in Angstroms converted from the frequency in the first column.

The remainder indicate the luminosity, that is :math:`L_{\nu}` of the system for specific types of photons. The units are :math:`{\rm erg\: s}^{-1} {\rm Hz}^{-1}`. 

The remaining columns are:

* **Created** is the total spectrum of all of the photons paakets as created, that is before having been translated through the wind
* **WCreated** is the spectrum of the photons that are created in the wind before translation
* **Emitted** is the emergent spectrum after the photons have been translated through the wind
* **CenSrc** is the emergent spectrum from photon bundles originating from the Star or BL, 
* **Disk** is the emergent spectrum from photon bundles originating from the disk,
* **Wind** is the emergent spectrum from photon bundles that have been reprocessed by the wind,
* **HitSurf** represents photons that did not escape the system but ran into a boundary 



.log_spec_tot_wind
  Identical to *.log_spec_tot* but just including photons that were generated in the wind or scattered by the wind



.spec
  an ascii file that contains the final detailed spectra for the wavelengths of interest at a distance of **100 pc** on a linear frequency scale.  
  The units for the detailed spectra are determined by the input parameter Spectrum.type.

  Photons bundles are generated in cycles in SIROCCO and the *.spec* file is actually written out at the end of each cycle
  as the program is running in the spectrum-generation phase of the program. So one can inspect the spectrum as it is building up.

  The beginning of the file (omitting the header) is as follows:

  .. code::


    Freq.        Lambda     Created    WCreated   Emitted    CenSrc     Disk       Wind       HitSurf    Scattered  A01P0.50   A30P0.50   A60P0.50   A80P0.50
    5.998778e+14  4997.560  5.5824e-13          0 5.5824e-13          0 5.5824e-13          0 1.0097e-15          0 1.9797e-12  1.141e-12 4.0282e-13  1.068e-13
    6.001705e+14  4995.122  6.4224e-13          0 6.4224e-13          0 6.4224e-13          0 1.3472e-15          0 2.0123e-12 1.2369e-12 5.1482e-13 1.0398e-13
    6.004632e+14  4992.687  7.2239e-13          0 7.2239e-13          0 7.2239e-13          0          0          0 1.8656e-12 1.2165e-12 4.9179e-13 1.3359e-13
    6.007560e+14  4990.254  7.4183e-13          0 7.4183e-13          0 7.4183e-13          0 6.7702e-16          0 1.7185e-12 1.4226e-12 5.9175e-13 1.6808e-13
    6.010487e+14  4987.824  7.9709e-13          0 7.9709e-13          0 7.9709e-13          0 3.3825e-16          0  2.262e-12 1.6291e-12 7.2959e-13 1.4697e-13



Where the first set columns are as follows:

* Frequency in Hz
* Wavelength in Angstorms
* The spectrum of photons which are created  (before passing through the wind)
* The spectrum of all photons which are created in the wind (before processing by the wind)
* The spectrum of all photons which escape the wind (after passing through the wind)
* The spectrum of all photons created by the star or BH (after passing through the wind)
* The spectrum of all photons created by the wind (after passing though the wind)
* The spectrum of all photons that are scattered by the wind (after passing through the wind)

These data in the first set of columns do not reflect the angular dependence of the emission. They are effectively an angle averaged spectrum. Except for the fact that the units are different and the wavelength range is limited these should resemble the spectra in the output files (such as .spec_tot) that record the spectra constructed in the ionization cycles.  


The remaining columns are the spectra at various inclination angles and binary phases.  The label A30P0.50 means the spectrum is viewed at an inclination angle of 30 degrees and at a phase of 0.5 -- for a binary system this is when the secondary was located behind the primary. 

**An important note to the user:** Due to the nature of modelling with wavelength bounds on the SEDs, the edges (minimum and maximum wavelengths) of observers emergent spectra drop in flux artificially. This is expected behaviour and is not a bug. As a result, if one is interested with a certain spectral region or a particular line. Be sure to widen the pre-determined wavelength minimum and maximum bounds by 50 or so Angstroms extra. This will ensure that your region of interest is appropriately modelled.

.log_spec
 Identical to the spectrum .spec file except with logarithmic intervals.  


.. note::
    Theree is no .spec_tot output spectrum.  The reason for this is that the wavelength/frequency range is so large in the tot spectrum that plots made 
    on a linear scale do not represtent the spectrum well.  For the detailed spectra, which are over a more limited spectral range both types of 
    spectra are produced.
