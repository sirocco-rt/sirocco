Wind.ionization.band_model
==========================
How the radiation field in each cell is described in the matrix ionization modes
(matrix_pow, matrix_est and matrix_multi).

SIROCCO does not normally use a detailed spectrum of the radiation in each cell.  Instead, as photons
pass through a cell during an ionization cycle, it records them in a small number (at most
20) of frequency bands, and at the end of the cycle fits a simple model of the mean
intensity J_nu, a power law or an exponential, in each band.  These band models are then
used wherever the code needs J_nu in the cell:

* the photoionization rates of the simple ions, which set the ionization (matrix_pow and
  matrix_multi);
* the populations of the levels of simple-atom lines, and induced Compton scattering and
  Compton cooling (matrix_pow and matrix_est).

A photoionization rate depends mostly on J_nu just above the ion's threshold, so the
results are much more accurate when the band boundaries lie on the thresholds of the
important ions.  This parameter chooses where the boundaries go, or replaces the band
models by the detailed cell spectrum.

These bands are not the bands used to generate photons, which are set by
:ref:`Photon_sampling.approach` (see :ref:`Photon Banding Strategies`); the two are the
same only with photon_bands.

Type
  Enumerator

Values
  hhe_cno_edges
    Band boundaries at the ground-state photoionization thresholds of H, He, O, C and N,
    taken element by element in that order and, within each element, from the neutral ion
    up; thresholds outside the photon frequency range, or very close to a boundary already
    placed, are skipped, and at most 10 are used.  The widest bands are then split until
    there are 20 bands.  (This is the default.)

  photon_bands
    The bands used to generate photons, set by Photon_sampling.approach.  If there are
    fewer than 7 of them, boundaries are added at the main photoionization edges.  This
    is how the bands were chosen before 2026.

  cell_spectrum
    No band models for the photoionization rates: J_nu is taken as constant in each bin of
    the cell spectrum, which then has 1000 log bins (unless -cell_spec_dim is given) plus a
    boundary at every photoionization threshold.  This is the most detailed description
    but needs many more photons per cell, and the cell spectra take much more memory.
    The band models (with the hhe_cno_edges bands) are still made for the other uses.

  user_bands
    Band boundaries chosen by the user: the number of bands (at most 20), then the
    boundaries between them in eV, increasing (Wind.ionization.band_model.nbands and
    Wind.ionization.band_model.band_boundary).  The outer limits are those of the photon
    generation bands, so that every photon falls in a band; a boundary outside them is
    dropped with a warning.  The widest bands are not split, so the bands are exactly
    those asked for.  This is intended for testing other choices of edges, for example
    for X-ray spectra.  For example, 4 bands with boundaries at the H I, He I and He II
    edges::

      Wind.ionization.band_model(hhe_cno_edges,photon_bands,cell_spectrum,user_bands)  user_bands
      Wind.ionization.band_model.nbands                 4
      Wind.ionization.band_model.band_boundary(eV)      13.6
      Wind.ionization.band_model.band_boundary(eV)      24.59
      Wind.ionization.band_model.band_boundary(eV)      54.42

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Wind.ionization`: matrix_pow, matrix_est, matrix_multi
