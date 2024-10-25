Plotting \& Processing Outputs
#####################################

SIROCCO produces a large number of files in both binary and ascii format. Tools exist to examine the binary files.

PySi 
==============
PySi, short for PySIROCCO, is a python package designed to interface with the SIROCCO code.
It can be installed from its `own Github repository <https://github.com/sirocco-rt/pysi>`_ 
under the sirocco-rt organisation. We recommend installing this package for reading code outputs.:
isntructions to do so can be found in the README of that repository.

.. admonition :: Warning to users

These tutorials are not complete, and our approach to interfacing 
with the code has not been uniform across the collaboration. Nevertheless, these 
tutorials should give you the basic tools needed to look at spectra and wind properties, 
with more detail on python scripts and packages provided in :ref:`the API documentation <py_progs>`. 

.. toctree::
   :glob:
   
.. nbgallery::
   plotting/plot_spectrum
   plotting/plot_wind