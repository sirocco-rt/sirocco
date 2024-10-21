Bound Bound
###########

This is the data for computing bound bound, or line interactions in simple/macro  atoms. 

Source
======

The Kurucz data used to create simple lines was taken from the `Kurucz website <http://kurucz.harvard.edu/linelists.html>`_.  
The file gfall.dat was used, though a few extra lines known to have been missing were likely added.



There are two main sources of data currently used in SIROCCO.

* Kurucz
* Verner
* Chianti

Kurucz is normally used for simple atoms whereas Chianti is the most common source for information about lines used in macro-atom versions
We have also used a line list from Verner in the past 


  
Translation to SIROCCO format
============================
There are several steps to creating the data used in SIROCCO from that in gfall.dat, that are carried out by py_read_kurucz and py_link. The first routine reads the gfall.dat file and creates two output files, a file containing the lines and the associated such as the effective oscillatory strength and a file which contains information about the ion levels.  py_read_kurucz chooses only a portion of the Kurucz lines, namely those associated with ions with ionization potentials in a certain range and lines with gf factors exceeding a certain value. The second program py_link attempts to create a model ion with links between the levels and the ions.  Both of these routines are driven by .pf files, similar to what are used in python.  Examples of the .pf files are in the directory py_kurucz

In practice we have not used these data for any SIROCCO publications. At some point early in the AGN project, NSH increased the number of lines, and generated lines\_linked\_ver\_2.dat and levels\_ver\_2.dat. I think this was because there was a small bug which meant the oscillator strength cut that was stated was not that which was applied.

Data format
===========



The lines have the following format

For lines, we did not create a specific topbase format, but most of the recent sets of 
data use a format that is similar to what is used  for macro atoms::

  Line  1  1  926.226013  0.003184   2   4     0.000000    13.387685    0    9
  Line  1  1  930.747986  0.004819   2   4     0.000000    13.322634    0    8
  Line  1  1  937.802979  0.007798   2   4     0.000000    13.222406    0    7
  Line  1  1  949.742981  0.013931   2   4     0.000000    13.056183    0    6

whereas for MacroAtoms::

  # z = element, ion= ionstage, f = osc. str., gl(gu) = stat. we. lower(upper) level
  # el(eu) = energy lower(upper) level (eV), ll(lu) = lvl index lower(upper) level
  #        z ion       lambda      f         gl  gu    el          eu        ll   lu
  LinMacro    1   1          1215.33907             0.41620     2     8             0.00000            10.19883     1     2
  LinMacro    1   1          1025.44253             0.07910     2    18             0.00000            12.08750     1     3
  LinMacro    1   1           972.27104             0.02899     2    32             0.00000            12.74854     1     4

For LinMacro the columns are 

* an identifier, 
* the element z, 
* the ion number, 
* the wavelength of the line in A, 
* the absorption oscillator strength, 
* the lower and upper level multiplicities, 
* the energy of the lower level and upper level. 
* the lower and upper level indicies (matched back to the energy levels)

The ultimate source for this information is usually NIST . The main issue with all of this is that 
one needs to index everything self-consistentl



SIROCCO structure
================

Line data is stored in the data structure **lines**

Comments
========
The version of gfall.dat that is used in SIROCCO is out of date, though whether this affects any of the lines actually used in python is unclear.  The website listed above has a link to newer versions of gfall.dat.


