Programming Notes
#################

SIROCCO is written in C and is normally tested on linux machines and on Macs, where the compiler usually turns out to be clang. It is also regularly compiled with gcc as part of the travis-CI tests. Certain portions of the code are parallelized using the Message Parsing Interface (MPI).

Version control is (obviously) managed through **git**.  The stable version is on the `master` branch; the main development is carried out on the `dev` branch. This is generally the branch to start with in developing new code. If possible, a developer should use the so-called Fork and Pull model for their version control workflow. See e.g. `this gist post <https://gist.github.com/Chaser324/ce0505fbed06b947d962>`_.

If one modifies the code, a developer needs to be sure to have ``$SIROCCO/py_progs`` both in ``PYTHONPATH`` and ``PATH``.  One should also have a version of indent installed, preferably but not necessarily gnu_indent installed.  This is because, the Makefile will call a script run_indent.py on files that the developer has changed, which enforces a specific indent style on the code.

In addition to indent, one should have cproto or something equivalent installed. cproto is used to prototypes for all of the subroutines (through the make command

.. code:: bash

    make prototypes

(The many warnings that appear when cproto is run on OSX can so far be ignored. cproto for macs is available with brew)

All new routines should have Doxygen headers.

printf statements should be avoided, particular in the main code.  SIROCCO has its own replacements for these commands, e.g Log and Error which standardize the outputsand allow for managing what is printed to the screen in multprocessor mode.  There is alsoa command line switch that contorls the amount of information that is printed to the screen.  Specific errors are only logged for a limited number of times, after which theyare merely counted.  A log of the number of times each error has occurred is printed outat the end of each run and for each thread.  (Additional detailes can be found in the Doxygenheader for xlog.c)

Structures
==========

In order to understand the code, one needs to understand the data structures.

The main header files  are:

* atomic.h - This contains all of the structures that hold atomic data, e.g oscillator
  strengths, photoionization cross-sections, elemental abundances, etc.  These data are
  read in at the beginning of the program (see atomicdata.c and other similarly named
  routines)
* python.h - This contains the structures and other data that comprise the wind as well
  as the parameters of the model.  (This is fairly well-documented, or should be)


Program Flow
============

Basically, (as decribed from a more scientific perspective elsewhere), the program consists
of a number stages

* Data gathering and initialization: This consiss of reading in all of the parameters
  for the model to be calculated, reading in the associated atomic data, and setting up
  program to run.  This procuess involves allocating space for many of the data structures.
* Ionization cycles: During this portion of the program fleets of photons are generated
  and propogated through the wind, interacting with it in various ways. These photons are
  generated over a large range of frequencies, because their purpose is to allow the program
  to determine the ionization state of the wind.  During this
  process various estimators are accumulated that describe the interaction of the photons
  with the wind.  Once all of the photons have propagated through the wind the various
  estimators are used to calculate a new estimate of the ionization state of the various
  wind cells that constitute the wind.  This process is repeated for a number of cycles,
  by which time, hopefully, the wind will have reached a "steady state".
* Spectral cycles: Once the ionization cycles have been completed, the ionization state
  of the wind is fixed, and more detailed spectra are calculated. For this, photons are generated
  in a limited spectral range, depending on the interests of the user.  In contrast to
  the ionization state, where "cycles" are  crucial concept, the only reason to have spectrl
  cycles in the "Spectral cycle" phase is to allow one to gradually improve the statistics
  of the spectrum.  At the end of each spectral cycle, the detailed spectra are written out,
  so one can see the spectra building up.


Parallel Operation
==================

SIROCCO uses MPI to parallelize the most compute intensive portions of the routine.  It has
been run on large machines with 100s of cores without problem.

The portions of the routine that are parallelize are:

* Photon generation and transfer: When run in multiprocesser mode, each thread creates only a
  fraction of the total number of photons.  The weight of the photons in each thread is such
  that the sum of the weights is the total energy expected to be produced in one observer frame second.
  These photons are propagated through the wind, and estimators based on these photons are accumulated.
  At the end of photon transfer by all threads, the various quantities, including the spectra,  that
  have been accumulated in the separate threads are gathered together and averaged or summed as
  appropriate.  For ionization cycles, this means that all of the data needed to calculate the
  ionization in any cell is available on each of the threads.
* Ionization calculation:  Although all of the threads have all of the data needed to calculate
  the ionization in any cell, in practice what happens is that the program assigns a different set of
  cells to each thread to calculate the ionization.  After the thread calculates the new ionization
  state for its assigned cells, the ionization states are then gathered back and broadcast to all
  of the threads, in preparation for the next cycle.
* Preparation for detailed radiative transfer in the macro-atom mode.  When photons go through the
  grid in the simple-atom mode, photon frequencies do not change a great deal, however in macro-atom
  mode the frequencies can shift by large amounts. This creates a problem during the detailed spectral
  generation stage, because one does not know before hand how many photons that started out of the
  desired band end up in the desired band.  (In macro-atom mode, a photon bundle that starts out at
  8 keV photoionizes an Fe ion can turn (for example) into an Hbeta photon).  To handle this, one
  needs to estimate how often this happens and include this (effectively as a source function) in
  radiative transfer involving macro-atoms. This is parallelized, in the same manner as the ionization
  calculation by assigning various cells to various threads and gathering the results back before
  the radiative transfer step in the detailed spectrum phase.


MPI requires intialization. For python this is carried out in python.c.  Various subroutines make
use of MPI, and as a result, programmers need to be aware of this fact when they write auxiliary
routines that use the various subroutines called by SIROCCO.

Input naming conventions
========================

As is evident from an inspection of a typical input file, we have adopted a somewhat hierarchical scheme
for the naming of the input variables, which groups variables associated with the same part of the system
together.  So for example, all of the variables associated with the central object have names like::

    ### Parameters for the Central Object
    Central_object.mass(msol)                  0.8
    Central_object.radius(cm)                  7e+08
    Central_object.radiation(yes,no)                  yes
    Central_object.rad_type_to_make_wind(bb,models)                   bb
    Central_object.temp                        40000


that is, they all begin with Central_object.  This convention should be followed.


External variables
==================

SIROCCO uses lots (and likely too many), what are properly know as  external variables.   (In C, a global
variable is a variable whose scope is all of the routines in a speciric file.  An external varriable
is one that is shared across multiple files.)

In the latest generations of gcc,  the standards for extenral variiables have been tightened.

If one wishes to define an external variable, one must first declare it as eternal, and then one
must initialize it outside a specific routine exactly in one place.

The standard convention is that the variables are declared as external in a header file, e.g python.h,
and then intialized in a separate .c file, e.g python_extern_init.c.   Unless, a variable is actually
initialized, no space will be allocated for the variable.

So if variables are added (or subtracted), one must make a change both in the relavant .h file.

Currently has three.c files atomic_extern_init.c, models_extern_init.c, python_extern_init.c
corresponding to the three main .h files, atommic.h, models.h and python.h


