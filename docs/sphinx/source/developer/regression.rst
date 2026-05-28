Regression Testing
##################

Regression tests verify that changes to SIROCCO do not introduce unexpected
differences in model outputs.  Developers are strongly encouraged to run these
tests before merging into any major branch.

Test Suites
===========

Standard regression tests (``examples/regress/``)
--------------------------------------------------

The standard suite covers a broad range of model types — CV, AGN, XRB, stellar
winds, 2D spherical, cylindrical, and polar (RTHETA) grids, macro-atom and
simple-atom line transfer, and imported hydro models.  These are the same
models run by the GitHub Actions CI workflow on every push to ``dev`` and
``main``.

3D regression tests (``examples/regress3d/``)
----------------------------------------------

A separate suite exercises the CYLIND3D and SPH3D coordinate systems introduced
on the x3d branch.  The models are:

- ``agn_cyl3d.pf`` — AGN wind on a 3D cylindrical grid (CYLIND3D)
- ``cv_cyl3d.pf`` — CV wind on a 3D cylindrical grid (CYLIND3D)
- ``agn_sph3d.pf`` — AGN wind on a 3D spherical-polar grid (SPH3D)
- ``cv_sph3d.pf`` — CV wind on a 3D spherical-polar grid (SPH3D)

Like the standard tests, these are not run to convergence — the goal is to
catch gross regressions in 3D-specific code paths (grid initialisation, photon
transport, wind-save I/O).

Polar regression tests (``examples/regress_polar/``)
-----------------------------------------------------

Models using the RTHETA (polar) coordinate system with full 0–180° coverage,
including lower-hemisphere observers, are kept in ``examples/regress_polar/``.

Running the Tests
=================

SIROCCO should be compiled with ``mpicc`` before running regression tests.

Standard and 3D tests
---------------------

The primary script is ``regression.py``.  A typical workflow::

    cd ~/Regression
    regression.py sirocco

This creates a dated run directory (e.g. ``sirocco_260526``), copies the
parameter files from ``examples/regress/``, runs each model sequentially with
the default number of MPI processes (currently 3), and — if a previous run
exists — calls ``regression_check.py`` and ``regression_plot.py`` to compare
spectra.  Comparison plots are written to an ``Xcompare/`` subdirectory.

To run the 3D suite instead::

    regression.py sirocco --indir $SIROCCO/examples/regress3d

To compare any two existing run directories directly::

    regression_check.py run1 run2

All scripts accept ``-h`` for a full list of options.

Inspecting Results with run_check3d.py
---------------------------------------

``run_check3d.py`` is an interactive diagnostic tool that handles 1D, 2D, and
3D wind models.  It generates a self-contained HTML summary with Plotly figures
(wind structure, spectra, convergence) for one or more models::

    run_check3d.py agn_cyl3d          # single model
    run_check3d.py -all               # every .wind_save in current directory

This is intended as the successor to ``run_check.py`` for 3D models, where
static matplotlib slices are insufficient to visualise the full wind geometry.

Interpreting Results
====================

The regression models are deliberately short runs (not converged).  Differences
between two versions may arise from:

- **Code changes** — genuine physics or numerical differences (investigate).
- **RNG reseeding** — a different execution path changes the photon sequence,
  producing noise-level differences that are not a concern.

There is no automatic way to distinguish these; inspection of the comparison
plots is required.

API Reference
=============

Full API documentation for the regression scripts is in the
:doc:`Python Scripts section <../py_progs/regression>`.
