Parameters
##########

A top-level parameter overview when initialising a .pf file. These parameters are universally required for all .pf files 99% of the time. The boolean parameters for this list are assumed to be no. If yes, the user is asked about further parameters.

.. toctree::
   :glob:
   :hidden:

   parameters/System_type
   parameters/central_object
   parameters/Disk
   parameters/wind
   parameters/rt_ion
   parameters/Spectrum
   parameters/other
   parameters/*

Primary .pf Parameters
***************************

* :ref:`System_type`

Parameters for the Central Object
---------------------------------

* :ref:`Central_object.mass`
* :ref:`Central_object.radius`
* :ref:`Central_object.radiation`

Parameters for the Disk (if there is one)
-----------------------------------------
* :ref:`Disk.type`
* :ref:`Disk.radiation`
* :ref:`Disk.temperature.profile`
* :ref:`Disk.mdot`
* :ref:`Disk.radmax`

Parameters for Boundary Layer or the Compact Object
----------------------------------------------------
* :ref:`Boundary_layer.radiation`

Parameters describing the various winds or coronae in the system
-----------------------------------------------------------------
* :ref:`Wind Model Parameters`

Parameters associated with photon number, cycles,ionization and radiative transfer options
--------------------------------------------------------------------------------------------
* :ref:`Photons_per_cycle`
* :ref:`Ionization_cycles`
* :ref:`Spectrum_cycles`
* :ref:`Wind.ionization`
* :ref:`Line_transfer`
* :ref:`Wind.radiation`
* :ref:`Surface.reflection.or.absorption`
* :ref:`Wind_heating.extra_processes`
* :ref:`Atomic Data`

Parameters for Domain 0
-----------------------
* :ref:`winds`

The minimum and maximum wavelengths in the final spectra and the number of wavelength bins
------------------------------------------------------------------------------------------

* :ref:`Spectrum.nwave`
* :ref:`Spectrum.wavemin`
* :ref:`Spectrum.wavemax`

The observers and their location relative to the system
-------------------------------------------------------
* :ref:`Spectrum.no_observers`
* :ref:`Spectrum.angle`
* :ref:`Spectrum.orbit_phase`
* :ref:`Spectrum.live_or_die`
* :ref:`Spectrum.type`

Parameters for Reverberation Modeling (if needed)
-------------------------------------------------

* :ref:`Reverb.type`

Others 
------
* :ref:`Photon_sampling.approach`


Diagnostic switches
-------------------
* :ref:`Diag.extra`
* :ref:`Diag.use_standard_care_factors`
* :ref:`Diag.write_atomicdata`
