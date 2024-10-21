Bound Free (TopBase)
####################


Source
======

Obtained from The `Opacity project <http://cdsweb.u-strasbg.fr/topbase/topbase.html>`_. See also `Cunto et at 1993, A&A, 275, L5 <http://articles.adsabs.harvard.edu/full/1993A%26A...275L...5C>`_


Translation to SIROCCO format
============================

ksl - It's not clear that we are now making use of the topbase data in this way but my original attempt to incorporate topbase photoinization data into SIROCCO is contained in the directory topbase. Processing of these files was done by py_top_phot. My feeling is that we can replace these remarks with those that are more up to date, once Nick and James discuss this section, and delete any mention of my original attempt on this in the data-gen archive.



Data format
===========
Our original photoionization cross sectiions came from a combination of \cite{verner96}, supplemented by a set of older values from \cite{verner95}

The TopBase photoionization files have the following format::


  PhotTopS  1  1 200    1    13.605698  50
  PhotTop    13.605698 6.304e-18
  PhotTop    16.627193 3.679e-18

whereas the Macro Atoms look like::

  PhotMacS       1       1       1       1       13.598430     100
  PhotMac       13.598430   6.3039999e-18
  PhotMac       13.942675   5.8969998e-18

The meaning of the columns is the same in both cases here. It may be simply a historical accident that we have both formats, or probably we were worried we would need to change. Topbase is generally the source for this information.

The initial line contains (a) the element, (b) the ion, (c) the level, (d) the level up, (e), the energy threshold in eV, and (f)  the number of x-sections to be read in.
The following lines gives the photon energy in eV, and the cross-section in cm$^{2}$.  To be a valid file the photon x-section  must not appear below the energy threshold

"Level up" corresponds to how many levels the electron is moving in the transition: this is simply 1
 
For the simple atom case, the header line can be parsed as follows
 
* z:  the atomic number
* ionization stage, in the astronomical  convention of the lower level ion, the one that gets ionized
* ISLP
* ll: the lower level index in the ion that gets ionized (with that ISLP)
* e: the energy threshold in eV (only photons above this energy ionize)
* npts: the number of points where the cross-section is measured
 
For the simple atom case the combination ISLP and level is unique
 
For the macro-atom case the entries in the header line are
 
* z:  the atomic number
* ionization stage, in the astronomical  convention of the lower level ion, the one that gets ionized
* ll: the lower level index in the ion that gets ionized (with that ISLP)
* ul: the upper level index in the ion after ionization (usually 1)
* e: the energy threshold in eV (only photons above this energy ionize)
* npts: the number of points where the cross-ection is measured
 
For the macro atom case, the indices relate to the energy levels, that is these must be connected up properly
 
The TopBase energies are inaccurate and so generally adjustments are made using Chianti or some other source to fix up the energy levels.

SIROCCO structure
================

The data are stored in the Topbase_phot stucture which can be found in atomic.h

Criteria for usage in SIROCCO run
================================

Data has to be read into SIROCCO in a logical order.  For a set of  phototionization x-sections to be accepted, the energy levels (or configuratios) must already have been defiend.  See :doc:`levels`

The items that must match are:

- the element (z) 
- the ion (istate)
- the upper level, which will be a level in the next ion (ilv)
- the lower level, which will be in the ion that is being photoionized


Comments
========

**The upper level in the MacroAtom case**

A common error that creates problems in reading in photoionization x-sections in the MacroAtom case is not to include the next ion up, partiulary the bare ion. If one encounters errors where the upper level is
not found, one should check the level file to verify that that the upper level ion is present, and that the inputs allow for the existence of at least the first level
of that ion.

For example, if one wishes to read in photoionization x-sections for N VII (hydrogenic), the levels file should include lines like::

    IonM    N   7   7   2  667.05100 1000   5     1s(2S_{1/2})
    IonM    N   7   8   1 1.0000e+20   1   1     Bare

The following is incorect::

    IonM    N   7   7   2  667.05100 1000   5     1s(2S_{1/2})
    IonM    N   7   8   1 1.0000e+20   0   0     Bare

because although the bare ion is present, the maximum number of levels is set to 0.   This is not an issue for the simple atom case.


**Extrapolation to higher energies**

Some topbase cross-sections do not extend to very high energies, for reasons that 
are not obvious.  This can cause non-physical edges to appear in spectra.  Therefore,
is is important to inspect any additions to the atomic data based on x-sections
retrieved from TopBas

Some tools have been developed To address this probllem.  In particularly,  JM wrote a 
script to extrapolate the cross-section to higher energies, by  
calculating the gradient in log-space at the maximum energy and extrapolating
to 100 keV. A number of cross-sections had unrealistic gradients at the original 
maximum energy, and were identified by eye and then forced to have a :math:`\nu^{-3}` shape.
This is the shape of a hydrogenic cross-section and whilst it is not accurate 
for non-hydrogenic ions, it is more realistic (and conservative) than some of 
the unphysically shallow gradients that were being found.
This is also briefly described in section~3.7.2 of Matthews PhD thesis.
The python scripts can be found in the `data-gen <https://github.com/agnwinds/data-gen>`_ repository progs/extrapolate\_xs/ 
with docstrings describing their use.

