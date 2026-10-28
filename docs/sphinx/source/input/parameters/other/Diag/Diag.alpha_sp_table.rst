Diag.alpha_sp_table
===================

Whether to tabulate the spontaneous recombination integral (alpha_sp) for
macro atoms as a function of electron temperature, instead of calculating it
afresh for every cell and trial temperature. The table covers 10\ :sup:`3` to
10\ :sup:`9` K with 40 points per decade, is made for each continuum the first
time it is needed, and is interpolated linearly in log alpha_sp versus log T.
With the ``segments`` method of :doc:`Diag.bf_integration` the interpolation
error is less than 1%, and usually about 0.01%. See issue #1219.

Type
  Boolean (yes/no)

Values
  no
    Calculate the integral directly (the default)

  yes
    Use the table


File
  `diag.c <https://github.com/sirocco-rt/sirocco/blob/master/source/diag.c>`_


Parent(s)
  * Only read when Sirocco is run with the ``-d`` flag
