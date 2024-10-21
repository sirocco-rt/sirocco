Direct Ionization
#################

This is the data to compute ionization rates from collisions between ions and hot electrons.


Source
======

The data comes directly from `Dere 2006, A&A, 466, 771 <https://www.aanda.org/articles/aa/pdf/2007/17/aa6728-06.pdf>`_ .  This paper gives  direct ionization and excitation-autoionization rate coefficients for many ions as a function  of temperature for Maxwellian electron distributions. 


Translation to SIROCCO format
============================


The data table is downloaded in its entirety  from the data table associated with the paper. All that happens is that the table is saved to a text file, and the keyword DI_DERE is just prepended to each row.


Data format
===========

Each line starts with the label DI_DERE and then follows


- Nuclear Charge - z - used to identify the ion
- Ion - state in our normal notation, so 1=neutral
- Number of splines N- the number of spline points for the fit of rate coefficients vs scaled temperature
- Scaled temperatures - there are N of these
- Scaled Rate coefficients - N of these

The scaled temperatures are  given by

:math:`x=1-\frac{\log{f}}{\log(t+f)}`

where :math:`t=kT/I`. :math:`I` is the ionization potential, and :math:`f=2.0`.
The rate coefficient R(T) is recovered from the scaled rate coefficient in the table, :math:`\rho` using

:math:`\rho=t^{1/2}I^{3/2}R(T)/E_{1}(1/t)`

where :math:`E_{1}` is the first exponential integral. In python we use the  gsl_sf_expint_E1 routine in gsl.

SIROCCO structure
================

This data is stored in the  dere_di_rate structure with members


- int nion- Internal cross reference to the ion that this refers to
- int nspline - the number of spline points that the fit is evaluated over
- double temps[DERE_DI_PARAMS]-  temperatures at which the rate is tabulated
- double rates[DERE_DI_PARAMS]-  rates corresponding to those temperatures
- double xi - the ionization energy of this ion
- double min_temp -the minimum temperature to which this rate should be applied


Comments
========
This data is also in Chianti , although in a different form. So we could potentially use this data as part of a push to just use Chianti for all our data uses. 
An updated set of DI data is available `here <https://arxiv.org/pdf/1702.06007.pdf>`_



