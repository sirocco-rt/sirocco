Stellar_wind.v_infinity
=======================
A parameter controlling the velocity of a stellar wind at large distances from the central object. Described in terms
of the Casters and Larmers equation,

.. math:: 
  v(r) = v_0 + (v_\infty - v_0) * (1 - R_s/r)^\beta.

Type
  Double

Unit
  cm/s

Values
  Greater than 0

File
  `stellar_wind.c <https://github.com/sirocco-rt/sirocco/blob/master/source/stellar_wind.c>`_


Parent(s)
  * :ref:`Wind.type`: ``star``


