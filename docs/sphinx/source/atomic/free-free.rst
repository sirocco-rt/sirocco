Free-Free Emission
##################



Source
======
The free-free Gaunt factors are taken from  `Sutherland 1998, MNRAS, 300, 321. <http://articles.adsabs.harvard.edu/full/1998MNRAS.300..321S>`_ The data is available for  download `here <http://www.mso.anu.edu.au/~ralph/data/freefree/>`_ where three files exist

- gffew.dat : Free-Free Emission Gaunt factors as a function of scaled electron and photon energy.
- gffgu.dat : Free-Free Emission Gaunt factors for Maxwellian electrons as a function of scaled temperature and scaled frequency.
- gffint.dat : Total Free-Free Emission Gaunt factors for Maxwellian electrons.


The last file is the one we use to calculate free-free emission, since this in integrated gaunt factor over a population of electrons  with a Boltzmann distribution of velocities.  The other two files could be of use in the future should we wish to have gaunt factor corrections for the heating rates,in which case we should use the gffgu.dat data file. However generally speaking free-free heating is never important and there would be significant overhead in calculating a gaunt factor for each photon.

Translation to python
=====================
The file is simply modified by hand to put a label "FF\_GAUNT" at the start of each data line and a hash at the start of each comment line.

Datafile - gffint.dat:
======================
The format of the data file to be read into python is as follows:

+----------+------------------------+---------------------------+-----------+-----------+-----------+
|Label     | :math:`\log(\gamma^2)` |:math:`<g_{ff}(\gamma^2)>` |:math:`s_1`|:math:`s_2`|:math:`s_3`|
+----------+------------------------+---------------------------+-----------+-----------+-----------+
|FF_GAUNT  |-4.00                   |   1.113E+00               | 1.348E-02 | 1.047E-02 |-1.855E-03 |
+----------+------------------------+---------------------------+-----------+-----------+-----------+
|FF_GAUNT  |-3.80                   |1.117E+00                  | 1.744E-02 | 9.358E-03 |5.565E-03  |
+----------+------------------------+---------------------------+-----------+-----------+-----------+
|FF_GAUNT  |-3.60                   |1.121E+00                  | 2.186E-02 | 1.270E-02 |4.970E-03  |
+----------+------------------------+---------------------------+-----------+-----------+-----------+



where  :math:`\gamma^2` is  the "scaled inverse  temperature experienced by an ion"
and the other four numbers allow the free-free Gaunt factor to be computed at any scaled inverse temperature 


:math:`x=\frac{Z^2}{k_BT_e}\frac{ 2\pi^2e^4m_e}{h^2}`

through spline interpolation between the two bracketing values of :math:`\log(\gamma^2)` 

:math:`<g_{ff}(x)>=<g_{ff}(\gamma^2)>+\Delta\left[s_1+\Delta\left[s_2+\Delta s_3\right]\right]`

where

:math:`\Delta=\log(x)-\log(\gamma^2)`

SIROCCO structure
================
This data is held internally in SIROCCO in the structure **gaunt_total** which has members

- log_gsqrd
- gff
- s1, s2, s3


Comments
========
We currently just use the total free free emission gaunt factor as a function of temperature for a Maxwellian distribution of electrons. This is OK for cooling, however we should really use the frequency dependant gaunt factor for free free heating. If we ever have a model where free-free heating is dominant, this should be looked into.
