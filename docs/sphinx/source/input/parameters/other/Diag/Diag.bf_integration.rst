Diag.bf_integration
===================

How the integrals over bound-free cross sections are done for macro atoms.
These are the spontaneous recombination rate (alpha_sp, used in the
temperature calculation), the band-limited version of it used in the
spectral cycles, and the dilute black-body photoionization and stimulated
recombination estimators used when the macro-atom level populations cannot
be calculated for a cell.

The TOPbase cross sections used in most macro-atom data sets have many narrow
resonances. For these the ``exact`` method is much faster than the others,
and more accurate. See issues #1146 and #1219, and PR #1220.

Type
  Enumerator

Values
  exact
    Gauss-Legendre integration over each interval between the tabulated
    points of the cross section, using the same interpolation as is used for
    the cross section elsewhere in the code, to machine precision. The
    integrals for the last two temperatures of each continuum are
    remembered. This is the default.

  sum
    A 10,000-point sum for alpha_sp; GSL integration for the band-limited
    alpha_sp and the dilute black-body estimators. This is what was used
    before this option was added. (The constant in the rates is now computed
    from the physical constants, so the rates are 0.094% lower than in
    earlier versions.)

  gsl
    GSL integration for all of the integrals (what was obtained before by
    compiling with ``MATOM_VER=2``).


File
  `diag.c <https://github.com/sirocco-rt/sirocco/blob/master/source/diag.c>`_


Parent(s)
  * Only read when Sirocco is run with the ``-d`` flag
