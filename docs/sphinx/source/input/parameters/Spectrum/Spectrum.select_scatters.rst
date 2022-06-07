Spectrum.select_scatters
========================
Advaned command that allows one to extract photons that
have undergone a certain number of scatters.  

* If n > MAXSCAT,that is to say a very large number, then all photons that could contribue to a spectrum are selected.
* If n lies between 0 and MAXSCAT then only contributions arising from photons that have
  scattered exactly n times will contribute to the spectrum. 
* If n is < 0 then contributions from photons with n or greater scattters will be extracted.

To understand this a little more clearly, consider a photon which undergoes a total on 
n scatters.  In extract mode, this photon will have made n+1 contributions to the total
spectrum, one when it was first emitted, one when it scattered the first time, one
when it scattered the second time, etc.  If one chooses to construct a spectrum from photons that
have one scatter, the contribution of this photon to the total spectra, at it's first scatter
will be reported.  



Type
  Integer

Values
  Greater than 0

File
  `setup.c <https://github.com/agnwinds/python/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Spectrum.select_specific_no_of_scatters_in_spectra`: ``True``


