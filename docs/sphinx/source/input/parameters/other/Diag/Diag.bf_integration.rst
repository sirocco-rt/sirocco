Diag.bf_integration
===================

How the integrals over bound-free cross sections are done for macro atoms.
These are the spontaneous recombination rate (alpha_sp, used in the
temperature calculation), the band-limited version of it used in the
spectral cycles, and the dilute black-body photoionization and stimulated
recombination estimators used when the macro-atom level populations cannot
be calculated for a cell.

The TOPbase cross sections used in most macro-atom data sets have many narrow
resonances. For these the ``segments`` method is much faster than the others,
and more accurate. See issues #1146 and #1219.

Type
  Enumerator

Values
  sum
    A 10,000-point sum for alpha_sp; GSL integration for the band-limited
    alpha_sp and the dilute black-body estimators. This is the default, and
    is what was used before this option was added.

  gsl
    GSL integration for all of the integrals (what was obtained before by
    compiling with ``MATOM_VER=2``).

  segments
    Gauss-Legendre integration over the segments between the tabulated
    points of each cross section, using the same interpolation as is used
    for the cross section elsewhere in the code.


File
  `diag.c <https://github.com/sirocco-rt/sirocco/blob/master/source/diag.c>`_


Parent(s)
  * Only read when Sirocco is run with the ``-d`` flag
