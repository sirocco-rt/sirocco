.. python documentation master file, created by
   sphinx-quickstart on Sun Jan 14 18:04:35 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

########
*python*
########
--------------------------------------
Radiative transfer and ionisation code
--------------------------------------

Python is a Monte-Carlo radiative transfer code designed to simulate the spectrum of biconical (or spherical)
winds in disk systems.  It was origianally written by
`Long and Knigge (2002) <https://ui.adsabs.harvard.edu/abs/2002ApJ...579..725L/abstract>`_ and
was intended for simulating the spectra of winds in cataclysmic variables. Since then, it has
also been used to simulate the spectra of systems ranging from young stellar objects to AGN.

The name Python is today unfortunate, and changing the name is an ongoing debate within the development team.
The program is written in C and can be compiled on systems runining various flavors of linux, including OSX on Macs.


The code is is available on `github <https://github.com/agnwinds/python>`_  Issues regarding the code and suggestions for improvement the code regarding the should be reported there.  We actively 
encourage other to make use of the code for their own science.  If anyone has questions about 
whether the code might be useful for a project, we encourage you to contact one of the authors of
the code.


-------------
Documentation
-------------

Various documentation exists:

* A :doc:`Quick Guide <quick>` describing how to install and run Python (in a fairly mechanistic fashion).

For more information on how this page was generated and how to create documentation for *python*,
look at the page for :doc:`documentation on the documentation <meta>`.

-------
Authors
-------
The authors of the *python* code and their institutions are:

Knox Long
  Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
  Eureka Scientific, Inc., 2452 Delmer St., Suite 100, Oakland, CA 94602-3017, USA

Christian Knigge
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

Stuart Sim
  School of Mathematics and Physics, Queen's University Belfast, University Road, Belfast, BT7 1NN, UK

Nick Higginbottom
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

James Matthews
  Institute of Astronomy, University of Cambridge, Cambridge, CB3 0HA, UK

Sam Mangham
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

Edward Parkinson
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

Mandy Hewitt
  School of Mathematics and Physics, Queen's University Belfast, University Road, Belfast, BT7 1NN, UK

----------------------------------------

.. toctree::
   :titlesonly:
   :glob:
   :hidden:
   :caption: Documentation

   quick
   installation
   running_python
   input
   output
   operation
   radiation
   wind_models
   coordinate
   examples
   physics
   atomic
   meta
   *
