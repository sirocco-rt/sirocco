Spectra Files
#############

Python is intended to produce simulated spectra.  These spectra are all ascii tables intended to be accessible with software packages such as astropy.

All of the ascii begin with commented headers that contain all of the parameters of associated with a run,
along with the date of the run and the specific version of Python used to make the run.
In principle, if one still has access to any of the spectra, one can reproduce the entire run again.

Broad band spectra are created from the last ionization cycle. (More accurately the broad band spectra are written out at the end of each ionization cycle, so one the program is finished one has the 
broad band spectrum of the last cycle)  

Detailed are calculated from all of the spectral cycles. (Properly normalized spectra are written out at the 
end of each spectral cycle, and with each cycle the photon statistics improves.)

The units in which the spectra are written out is also indicated in the header to the file.  


For a model with root name *cv*, the following broadband spectra will be created:

* **cv.spec_tot**
* **cv.log_spec_tot**
* **cv.spec_tot_wind**
* **cv.log_spec_tot_wind**

File types
==========

.spec_tot
  An ascii file that contains various spectra from the ionization-calculation phase of the program on a linear frequency scale.
  The first few lines of the file (omitting the header) are as follows:

  .. code::

    Freq.        Lambda     Created    WCreated   Emitted    CenSrc     Disk       Wind       HitSurf
    2.524334e+14 11876.102  3.5244e+18          0 3.5244e+18          0 3.5244e+18          0 1.1547e+16
    2.550397e+14 11754.737  3.4721e+18          0 3.4721e+18          0 3.4721e+18          0 2.8761e+15
    2.576461e+14 11635.827  3.4433e+18          0 3.4433e+18          0 3.4433e+18          0 2.8835e+15
    2.602524e+14 11519.299  3.6858e+18          0 3.6858e+18          0 3.6858e+18          0 2.8706e+15
    2.628587e+14 11405.082  3.6711e+18          0 3.6711e+18          0 3.6711e+18          0 1.1528e+16


The first two columns are fairly obvious. Lambda is in Angstroms. 

The remainder indicate the luminosity, that is :math:`L_{\nu}` of the system for specific types of photons. The units are :math:`ergs\: s^{-1} Hz^{-1}`. 

The remaining columns are:

* Created is the total spectrum of all of the photons paakets as created, that is before having been translated through the wind
* WCreated is the spectrum of the photons that are created in the wind before translation
* Emitted is the emergent spectrun after the photons have been translated through the wind
* CenSrc is the emergent spectrum from photons bundles originating on the Star or BL, 
* Disk photons the spectrum due to photons starting in the disk
* and Wind are the same for photons originating in the disk and wind respectively. 
* HitSurf represents photons that did not escape the system but ran into a boundary 


.log_spec_tot
  An ascii file which contains the same information as *.spec_tot*, but with a logarithmically space frequency intervals.
  This gives better sampling of the SED in a lot of cases and is much better for plotting things such as the input spectrum.

.spec_tot_wind
  Identical to *.spec_tot* but just including photons that were generated in the wind or scattered by the wind

.log_spec_tot_wind
  A logarithmic version of *.spec_tot_wind*




.spec
  an ascii file that contains the final detailed spectra for the wavelengths of interest at a distance of **100 pc**.  The units for the detailed spectra are determined by the input parameter Spectrum.type.

  Photons bundles are generated in cycles in Python and the *.spec* file is actually written out at the end of each cycle
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
* The spectrum of all photoons created by the wind (after passing though the wind)
* The spectrum of all photons that are scattered by the win (after passing through the wind)

These data in the first set of columns do not reflect the angular dependence of the emission. They are effectively an angle averaged spectrum.  
Except for the fact that the units are different and the wavelength range is limited these should resemble the spectra in various ouput files showing 
the spectra constructed in the ionization cycles.  


The remaining columns are the spectra at various inclination angles and binary phase.  The name A30P0.50 implies the spectrum at an incination angle 
of 30 degrees and an for a binary system one obtained when the orbital period was 0.5 that is when the secondary was located behind the primary.

.log_spec
 Identical to the spectrum .spec file except with logarithmic intervals.  
