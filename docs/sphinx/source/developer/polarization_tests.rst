Polarization Tests
##################

The tests in :code:`$SIROCCO/py_progs/Tests/polarization` check that the polarization computed by SIROCCO (see
:doc:`../physics/polarization`) is right.  Each test is a full SIROCCO run, made with
:code:`@Spectrum.polarization yes`.  The result is compared with a calculation made independently of SIROCCO:
an analytic result, a published one, or a separate numerical calculation.  They complement the other kinds of
testing:

* the unit tests (:doc:`tests`), including :code:`test_polarization.c`, check individual routines;
* the regression tests in :code:`examples/regress` check that results do not change between versions, and
  include :code:`pol_cone`, a short version of one of the tests here;
* the polarization tests check that the answer is correct.

They should be run after any change to electron scattering, the frame transformations, the extraction
(peel-off) of photons, or the spectra.

Running the tests
=================

The tests are run by :code:`pol_tests.py`, in a working directory outside the repository:

.. code:: bash

    $ python $SIROCCO/py_progs/Tests/polarization/pol_tests.py ~/pol_test_runs          # quick tests, ~10 minutes
    $ python $SIROCCO/py_progs/Tests/polarization/pol_tests.py -full ~/pol_test_runs    # full tests, several hours

The main options are :code:`-np` (the number of MPI processes, 8 by default), :code:`-exe` (the SIROCCO
executable to test) and :code:`-w2t` (a :code:`windsave2table` built from the same code, which one test needs).
:code:`-tests` selects some of the tests, and :code:`-check_only` repeats the checks on results already in the
working directory.

For each test, :code:`pol_tests.py` copies the files it needs into a subdirectory of the working directory,
makes and runs the models there, runs the analysis, and compares the results with tolerances.  Nothing is written
in the repository.  The results go to :code:`Summary.txt` in the working directory, which records the
executable and its version, and gives the value, the allowed range and PASS or FAIL for each check.  The exit
status is 0 only if every check passes.

The tests
=========

All the tests use electron scattering only, the process whose polarization SIROCCO calculates.

.. list-table::
   :header-rows: 1
   :widths: 15 40 45

   * - Test
     - What it is
     - What is checked
   * - :code:`brown_mclean`
     - Optically thin, flattened envelopes around a point source
     - The polarization against exact single scattering (Brown & McLean 1977, with attenuation), allowing for
       multiple scattering
   * - :code:`hillier`
     - Hillier's (1994) prolate shell, the test used by Bulla et al. (2015), for optical depths from 0.03 to 3
     - The polarization against Hillier's results, as a function of inclination, optical depth and albedo, to
       0.1-0.3%
   * - :code:`sr_shell`
     - A thin shell moving radially at up to :math:`0.3c`, and the same with :code:`-nonrel`
     - The scattered I and Q against an exact special-relativistic single-scattering calculation, to a few per
       cent
   * - :code:`cvcone`
     - A CV-like wind in a polar cone or an equatorial wedge, lit by a central source, with the ionization
       calculated by SIROCCO
     - The polarization of the light scattered once against an independent single-scattering calculation, to 1
       point

Each test has a quick version, of a few minutes, for checking a change, and a full version, which reproduces
the validation described in :ref:`pol_validation`.  The README in :code:`py_progs/Tests/polarization` describes each
test in more detail, and the :code:`testdef.py` in each test directory documents its tolerances and how they
were set.

Known issues
============

* For the wedge model of :code:`cvcone` at inclinations of 10 and 45 degrees, SIROCCO gives a polarization
  about 2 points more positive than the single-scattering prediction.  This is not understood, and these two
  angles are reported but not checked.
* For the thinnest envelope in :code:`brown_mclean` (:math:`\tau_{eq} = 0.01`), SIROCCO's polarization has been a
  few per cent below single scattering, where multiple scattering should make it slightly higher.  It passes,
  but is close to the limit.
