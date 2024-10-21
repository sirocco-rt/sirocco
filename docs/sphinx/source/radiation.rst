Radiation Sources
#################

.. todo::
  Fill in. Add description of how to use your own input spectrum. Finish links to keywords.

External Radiation Sources
==========================

In generic terms, there are two main external radiation sources for any SIROCCO calculation: a **Central Source** which can be a normal star, a WD, or a BH, and a **disk**.  Even though SIROCCO supports the existence of a secondary star for the purposes of calculating when light from a disk system is occulted, the secondary star does not radiate.

Photons for radiation from the central object emerge uniformly over its surface, except when a lamp-post geometry is
specified for the ``bh`` or ``agn`` system types. In this lamp-post case, radiation originates from a point source above
and below the central object, with a specified height.

Emission from a boundary layer can also be defined when this is relevant, from which radiation also emerges uniformly over the surface of the central object.

The Wind as a radiation source
==============================

In macro-atom calculations the wind is NOT a radition source.  All of the photons in a macro-atom
calculation are generated externally, and with minor exceptions photons preserve their weight
throughout their passage through the wind.  (The minor exceptions have to do with processes like
adiabiatic cooling, which result in the loss of photons).

In the simple-atom approach, various processes cause photons passing through the wind to lose energy
as they pass through the wind. This energy heats the plasma.  To account for this, photons are
generated from the wind at the beginning of each cycle.  Processes include, free-free emission, free-bound emission and line emission.

In non-macro-atom calculations wind radiation can be turned on and off using the
:doc:`/input/parameters/wind/Wind/Wind.radiation` keyword.

(In various files that contain the spectra there is a column WCreated that in the simple atom mode gives the spectrum of photons that were created in the wind.  This column, also exists in the macro-atom case, where it records the spectrum of pbotons that have interacted with the wind and been re-emitted.)


Spectra of the external radiation sources
=========================================

For the most part, the various radiation sources can radiate by any of the following process, as appropriate)

1. Blackbody radiation, specified in terms of a temperature.  Depending on the nature of the source, the luminosity
   specified either by the size of the object, or directly as the total luminosity.

2. Bremsstrahlung radiation, specified in terms of a temperature and a luminosity between 2 and 10 keV

3. Power law radiation, specified in terms of a spectral index, and a luminosity between 2 and 10 keV

4. One or more spectral models read from a series of files.  The models must specified in terms of two
   parameters, usually T and log g, each model consists of an ascii file containing the spectra.
   An example of the ascii files
   that can be read in is contained in the xdata folder that is part of the distribution (See below).


In the ionization cycles, the spectra of the central source, boundary layer (if present) and disk are determined by these three keywords:

* :doc:`/input/parameters/central_object/Central_object/Central_object.rad_type_to_make_wind`
* :doc:`/input/parameters/central_object/Boundary_layer/Boundary_layer.rad_type_to_make_wind`
* :doc:`/input/parameters/Disk/Disk.rad_type_to_make_wind`

It is possible to choose different input spectra for the ionization and spectral cycles, so a corresponding keyword of
the form :doc:`/input/parameters/Disk/Disk.rad_type_in_final_spectrum` is also needed.


Spectra from a model grid (details)
===================================

SIROCCO was initially written to model the winds of cataclysmic variables (CVs).  Although the spectra of the disks of cataclymic
variables are often modelled in terms of blackbodies, the spectra of CVs show clear evidence of features that arise from the i
disk (as well as the wind).   The features that arise from the disk resemble in many respects those that arise from
an appropriately weighted set of stellar atmospheres.  To allow for this possibility, SIROCCO can be configured to read
a set of
models characterized by a temperature and log g, and produce spectra of either the central object or the disk by interpolating on t
and log g.  The data that must read in consists of a file that associates a temperature and log g with the indvidual spectra.

For example, as part of the standard distruction there is a file kurucz91.ls, which starts as follows



::

    data/kurucz91/fp00t3500g00k2c125.txt         3500          0.0
    data/kurucz91/fp00t3500g05k2c125.txt         3500          0.5
    data/kurucz91/fp00t3500g10k2c125.txt         3500          1.0
    data/kurucz91/fp00t3500g15k2c125.txt         3500          1.5
    data/kurucz91/fp00t3500g20k2c125.txt         3500          2.0
    data/kurucz91/fp00t3500g25k2c125.txt         3500          2.5
    data/kurucz91/fp00t3500g30k2c125.txt         3500          3.0
    data/kurucz91/fp00t3500g35k2c125.txt         3500          3.5
    data/kurucz91/fp00t3500g40k2c125.txt         3500          4.0
    data/kurucz91/fp00t3500g45k2c125.txt         3500          4.5
    data/kurucz91/fp00t3500g50k2c125.txt         3500          5.0
    data/kurucz91/fp00t3750g00k2c125.txt         3750          0.0
    ...

In this case we have spectra at a temperature of 3500, for 11 different values of log g,
before going on to temperature of 3750 K.  Each spectrum is one of the Kurucz models, and these contain entries which contain a
set of wavelengths and a quantity that is understood to be proportional to :math:`F_{\lambda}`.

The 3 column format above is required.  If one wants to use a set of models that have only a T parameter one should
simply
choose a value for the second column.  The use case here is fairly specific, especially with regard to the first parameter T.
If the disk or central object temperature outside the
temperatures in the grid, then SIROCCO will "adjust" the spectrum assuming that the overall spectrum changes as a BB
would, but
the features in the spectrum are unchanged.  If the gravity goes outside the range of the grid, the closest value is chosen.

One need not use Kurucz models, of course.  Any set of models can be used, as long as the files contain two
columns, a wavelength in Angstroms and something that is proportional to :math:`F_{\lambda}`.  The normalization of the fluxes
does not matter, because the models are only used to establish the shape of the spectrum.  The normalization is
determined by the total luminosity of the component.

.. toctree::
   :glob:

   radiation/*

