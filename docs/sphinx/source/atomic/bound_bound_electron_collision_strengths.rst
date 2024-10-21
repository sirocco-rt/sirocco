Bound-bound electron collision strengths
########################################

Source
======
We use the Chianti atomic database, specifically the \*.scups files. These "contain the effective electron collision strengths 
scaled according to the rules formulated by 
`Burgess & Tully 1992, A&A, 254, 436 <https://ui.adsabs.harvard.edu/abs/1992A%26A...254..436B/abstract>`_
The values in the file are functional fits to :math:`\Upsilon(T)` which we referred to as :math:`\Omega` in our calculations for collisional de-excitation rate coefficient


:math:`q_{21}=\Omega\frac{8.629\times10^{-6}}{g_u\sqrt{T}}`

In the g-bar formulation

:math:`\Omega=4.77\times10^{16}g_l\overline{g}\frac{f_{abs}}{\nu}`

These values of :math:`\Upsilon` simply replace :math:`\Omega`.

In the asbsence of data in this format, the Van Regemorter approximation is utilized.

Translation to SIROCCO format
============================

It is necessary to link each line in our line list with the relevant electron collision strength. This is achieved using the python script "coll_stren_lookup.py" which first reads in the  "lines_linked_ver_2.py" line list, then attempts to work out which lines are which by comparing the energy and the oscillator strength of the line. If these match to within a factor of 10% then the code logs this as a possible match. If better matches come along, then the code adopts those instead.

Each matched line get a line in the data file which is basically all of the line data for the matched line. This is to give SIROCCO the best chance of linking it up with the line internally.

Data format
===========

The collision strength data has the following format::

  CSTREN Line  1  1 1215.673584  0.139000   2   2     0.000000    10.200121    0    1       1      3   7.500e-01   2.772e-01   1.478e+00    5    1   1.700e+00
  SCT   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
  SCUPS    1.132e-01   2.708e-01   5.017e-01   8.519e-01   1.478e+00
  CSTREN Line  1  1 1215.668213  0.277000   2   4     0.000000    10.200166    0    2       1      4   7.500e-01   5.552e-01   2.961e+00    5    1   1.700e+00
  SCT   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
  SCUPS    2.265e-01   5.424e-01   1.005e+00   1.706e+00   2.961e+00
  CSTREN Line  1  1 1025.722900  0.026300   2   2     0.000000    12.089051    0    3       1      6   8.890e-01   5.268e-02   2.370e-01    5    1   1.600e+00



Each record has three lines. The first line has a keyword CSTREN and this contains all the line data for the line to which this collision strength refers as the first 10 fields. These fields are identical to the 10 fields that appear in a standard line file. The next 8 fields are

- 1-2 Upper and lower level of transition - Chianti nomenclature
- 3 Energy of transition - Rydberg
- 4 Oscillator strength x lower level multiplicity (GF)
- 5 High temperature limit value
- 6 Number of scaled temperatures - 5 or 9
- 7 Transition type \cite{1992A&A...254..436B} nomenclature
- 8 Scaling parameter (C) (Burgess & Tully 1992) nomenclature

The next two lines, with labels SCT and SCTUPS are the 5 or 9 point spline fits to  :math:`\Upsilon` vs T
in reduced units y,x.

There are four different types of transitions, each with their own scaling between temperature of the transition and :math:`\Upsilon`.

For example, for type 1 (the most common)

:math:`x=1-\frac{lnC}{ln\left(\frac{kT}{E_ij}+C\right)}`

and

:math:`y(x)=\frac{\Upsilon}{ln\left(\frac{kT}{E_{ij}}+e\right)}`

So, to get :math:`\Upsilon` for a given T, one converts T to x via the correct equation, then linearly interpolate between values of :math:`y(x)`, then convert back to :math:`\Upsilon`.

SIROCCO structure
================

The data is stored in SIROCCO in the Coll\_stren structure which has memebers


- int n - internal index
- int lower, upper - the Chianti levels, not currently used
- double energy - the energy of the transition
- double gf - the effective oscillator strength - just oscillator strength x multiplicity
- double hi_t_lim - the high temperature limit of y
- double n_points -The number of points in the spline fit
- int type - The type of fit, this defines how one computes the scaled temperature and scaled coll strength
- float scaling_param - The scaling parameter C used in the Burgess and Tully calculations
- double sct[N_COLL_STREN_PTS] -The scaled temperature points in the fit
- double scups[N_COLL_STREN_PTS]- The sclaed coll sttengths in ythe fit


There is also a member in the line structure (coll_index) which points to the relevant record

Comments
========


There are currenly 4 types of transitions that are read from the Chianti data

- 1 - Allowed transition
- 2 - Forbiddent transitions
- 3 - Intercombination trantions
- 4 - Allowed transitions with 

which correspond to the transition types idenitfied by Burgess & Tully

There are addtional transition types in the Chianti database

- 5 - Dielectronic capture tranisitions
- 6 - Proton transitions


The latter are not currently used in *SIROCCO*

Discussion of how Chianti handles transitions can be found in 
`The CHIANTI upsilon files (ups and scups) <http://www.chiantidatabase.org/tech_reports/13_scups/chianti_report_13.pdf>`_

