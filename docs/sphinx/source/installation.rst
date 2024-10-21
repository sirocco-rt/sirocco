Getting Started
###############

What machines will python run on? SIROCCO is written C.  We have regularly run SIROCCO on various linux systems as well as a variity of Mac machines.
It is compiled using mpicc/gcc, but it can be compiled simply with gcc.

It uses the Gnu Scientific Libraries (gsl)

(Developers should also have cproto in their path in order to create new prototypes, and access to indent to insure that routines are formatted in a standard fashion. They will also want to make sure the py_progs routines are properly installed, as indicated below).

Installation
============

SIROCCO and the various routines associated are set up in a self-contained directory structure.
The basic directory structure and the data files that one needs to run SIROCCO need to be retrieved and compiled.

If you want to obtain a stable (!) release, go to the `Releases <https://github.com/sirocco-rt/sirocco/releases/>`_ page.

If you want to download the latest dev version, you can zip up the git repository by clicking on the zip icon to the right of the GitHub page.
Alternatively, clone it directly as

.. code:: bash

    $ git clone git@github.com:sirocco-rt/sirocco.git

If you anticipate contributing to development we suggest Forking the repository and submitting pull requests with any proposed changes.

Once you have the files, you need to cd to the new directory and set your environment variables

.. code:: bash

    $ export SIROCCO=/path/to/sirocco/
    $ cd $SIROCCO
    $ ./configure
    $ make install
    $ make clean

One can run a more rigorous clean of GSL with :code:`make distclean`, or remove the compiled GSL libraries altogether with :code:`make rm_lib`.

note that export syntax is for bash- for csh use

.. code:: console

    $ setenv SIROCCO /path/to/sirocco/


The atomic data needed to run SIROCCO is included in the distribution.


The source code for SIROCCO is under actively development and is updated fairly often. Normally, one does not need to redo the entire installation process, since this includes GSL setup.
Instead, one can pull in changes and recompile the source code by running

.. code:: bash

    $ cd $SIROCCO/source
    $ make sirocco

which will compile the main program. The program plus full set of auxiliary programs (such as windsave2table and py_wind, see below) can be compiled using :code:`make all`.

Running python
==============

To run python you need to add the following to your $PATH variable:

.. code:: bash

    $SIROCCO/bin

You can then setup your symbolic links by running

.. code:: bash

    $ Setup_SIROCCO_Dir

and run the code by typing, e.g.

.. code:: bash

    $ sirocco root.pf


Running in parallel mode
------------------------

While SIROCCO can be run in single processor mode, it is generally more efficient to run on multiple processors. in multiprocessor mode,
When multiprocessing is invoked, SIROCCO uses mulitple threads for photon transfer and in calcuation ionization equilibrium.  As these
comprise the bulk of the computational load the total time to run is  roughly an inverse of the number of threads.  SIROCCO uses `MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ for parallel processing and so software libraries that implement this must be on the machine that is
being used.  For Macs, mpi can installed with HomeBrew or Fink.  For linux machines, two common libraries are `Open-MPI <https://www.open-mpi.org/>`_ and `MPICH <https://www.mpich.org/>`_  If not already installed, one should
install them.

With mpi installed (and after recompiling with mpicc, which is the default) one would simply run the above program with

.. code:: bash

    $ mpirun -np 8  sirocco root.pf

where  -np followed by a number designates the number of threads assigned.


Auxiliary programs
------------------

There are two programs that are useful for extracting information about models

* windsave2table generates a series of astropy tables that can be used to inspect elements of the various models, including densities of specific ions
* sirocco_wind is a mainly interactive routine that prints similar infomation to the screen.

The two files are run as follows

.. code:: bash

    $ windsave2table root
    $ sirocco_wind root

Brief descriptions of command line options for running these routines can obtained using a -h switch

SIROCCO scripts
--------------

There are a number of Python, the programming language scripts, that can be used to plot results
from a SIROCCO run.  These are not particularly well documented and many have been developed
for looking at various aspects of the code.  A few may require Python packages to be installed.
However, a number are likely to be useful.

To make use of these scripts one should add

$SIROCCO/py_progs both to the PATH and PYTHONPATH variables

One script that is particularly useful is run_check.py, which is run as follows

.. code:: bash

    $ run_check.py root


This should create an html file that contains a summary set of information about a run, with plots that
indicate how much of the wind has converged as a function of cycle, which cells have converged at the end, what
the electron and temperature structrue of the wind is, as well as quick plots of the spectra that were produced.

Directory structure
-------------------

The python directory structure is fairly simple:

source
  Location of source code

bin
  Location of executables

docs
  Location of documentation, including sphinx docs, doxygen, parameters and documentation for the python programs in py_progs.

data
  Location for all datafiles. Files that are mainly for reference should be gzipped to save space. Such files are not recreated in

bin
  The location of the executables. (It is a good idea to put this directory in your path)

software
  This directory contains libraries which are used in in python that must be recompiled when creating an installation on a new machine, primarily Bill Pence's cfitsio package and the GNU scientific library gsl

py_progs
  python programs for helping analyse the code. We recommend adding this directory to your PATH and PYTHON_PATH environment variables.

examples
  A directory with a few examples of python runs. (Note that the input files will have changed and so one may not be able to run these examples without some changes in the input files.)

Please help by reporting bugs in installation
---------------------------------------------

This can be done by submitting a bug under the `Issues <https://github.com/sicorro-rt/sirocco/issues/>`_ page
