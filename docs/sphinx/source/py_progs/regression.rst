Regression
----------

Primarily to verify that changes made to SIROCCO do not inadvertently cause unexpected changes
if models, several routines exist to run a fixed set of (relatively fast) models that are
**nominally**  contained in Examples/regress.

Developers are encouraged to use this routines, before they merge anything into one of the
major branches of SIROCCO.

The routines involved are

 * `regression.py`
 * `regression_check.py`
 * `regression_plot.py`

The primary routine is regression.py.  All of the routines can be run with a -h switch
to obtain information of the full set of command line options.

Setup
=====

Typically one should set up a directory, e.g Regression to run the routines, and, if for example,
sirocco87f, is the version of SIROCCO being used when you set up the directory, being run.

SIROCCO should be compiled with mpicc before running the regression program

Basic Procedure
===============

Run::

    regression.py sirocco87f

This will create a directory sirocco87f_231108 where 231108 is the current date.  The pf files from
the regression directory as well as various ancillary files will be copied into this directory,
and all of the models contained therein will run sequentially.
In the absence of command line
switches the routines will be run with a default number of processors (currently 3).
Assuming this is the first time the program is run, no comparison to previous runs will be made

The models have been selected
to test a variety of types of models and/or to highlight areas of concern. As a result, the models that are run are likely
to change occassionaly.

Once changes have been made to python, one reruns the program, e.g::

    regression.py sirocco

This will create a directory sirocco_2311108 (assuming it this is the same day) and repead the previous
precedured.

**If the program is run on the same day with the same version of python, the older models
will be overwritten.  Typically one can avoid this by using py one time and py with the version number
a second time.  But there is also a command line switch to specify the name of the run time directory**

Assuming all of the various models run to complesion, regression.py will call subroutines in regression_check.py
and regression_plot.py to make comparasions between the model just run and the previous one.  The plots (one for each model)
will be contained in a directory called Xcompare.


Interpretation of the results
==============================

The models that are run are simple models, and to allow one to proceed quickly, none of the models is run to convergence.
The outputs compare the spectra that were produced in the two runs of the program, both by looking to see how many lines in
the ionization and detailed spectra have changed, and by generating plots that show comparisions of the spectra.

Many times the results will be identical, but if a change between two versions of the program results in a different
sequence of random numbers, then the spectra will change simply as a result of random noise, which is not a concern.
There is no easy way to quantify changes that are due to this effect or something else, and
so one simply through experience has to gauge the results by inspecting the plots that are produced..


Comments and additions
======================

Although regression.py generally produces a comparison betwen the set of models being run and the last set of models that wre run, one can use
regression_check.py to compare any two sets of runs::

    regression_check.py run1 run2

where one gives the names of the two directories to be compared.

While the regression procedure described here is generally set up to run on the models that are contained in the Examples/regress directory,
regression.py has switches that allow one to do tests on models that are in any input directory.  This can be useful, if one wishes to test
different models in order to solve specific problems, or to run a set of models sequentially.

API Documentation
=================

.. autosummary::
    :toctree: regression

    regression
    regression_check
    regression_plot
    regression_nsh
