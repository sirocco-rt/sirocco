Running SIROCCO
##############

The normal way to run SIROCCO is simply to enter

.. code :: bash

    sirocco xxx

where ``xxx`` is the root name of a parameter file.  (The full name ``xxx.pf`` can also
be entered).

However SIROCCO features a number of command line options which can be used
to modify it's operation.  These include the following:

-h
  Causes SIROCCO to print out a brief help message and quit. The help message
  principally describes the command line options.

-i (or --dry-run)
  Causes SIROCCO to read and verify the inputs, writing a clean version of the input
  file ``xxx.pf`` to the output file ``xxx.out.pf``, and then stop. This option is useful
  for setting up a proper ``.pf`` file.  (Often one will want to copy ``xxx.out.pf`` back
  to ``xxx.pf`` before proceeding.

-t time_max
  Limits a run of SIROCCO to approximately time_max in sec.  This switch is
  used in situations where one would like to check whether the routine is operating
  properly be continuing, or where one needs to checkpoint the program after a certain
  period of time (due for example to time limits placed on jobs in a Beowulf cluster).
  The time is checked at the end of ionization and spectral cycles, immediately after
  saving the binary files that describe a model, and so one needs to leave a cushion
  between time_max and the maximum time one wants the program to run

-r
  Restarts a run that has been interrupted or halted, by reading a the ``xxx.wind_save``
  and ``xxx.spec_save`` file (if it exists).  Note that very few values in the .pf
  file are read when this options is used, as most of the information there has
  already been utilized in setting up and executing the run. The main ones that
  can be changed are the numbers of cycles for either ionizaion or detailed spectral
  cycles.  Parameters that will be ignored include those assoicated with the wavelength
  range and extraction angles of the detailed spectra.  The way to make changes to
  the detailed spectra is usually to use the option of setting the System\_type to previous,
  which will allow one to set all of the detailed spectral parameters anew.

-v n
  Changes the amount of information printed to the screen by SIROCCO during a run.
  The default is 4.  Larger numbers increase this. Smaller numbers decrease it.
  The log files are not affected.

--rseed
  Causes SIROCCO to use a random number seed that is time-based, rather than fixed.
  In most cases, a fixed seed is preferred so that problems can be replicated, but if
  is repeating the same calculation multiple times, then one may want a random seed.

--rng
  Save or load the RNG state to file, to allow persistent RNG states between restarts

--version
  Causes SIROCCO to print out the version number and commit hash (and whether
  uncommitted files exist, and then stop.

-p n_steps
  Changes the number of photons generated during ionization cycles so that the
  number increases logarithmically to the maximum value. This can help speed up SIROCCO
  simulations but check the covergence of the wind. The number ``n_steps`` is optional,
  and specifies the number of decades over which the increase takes place.



Special switches
================

SIROCCO has a number of other switches that are not intended for the general user, but
which may be useful in certain special cases.  These include:

-d
  Enables a variety of specialized diagnostic inputs which have been implemented
  to help with solving various problems, and were regarded (by someone) as useful
  enough to maintain in the program.  The user is then queried regarding which
  of these diagnostics to enable for a specific run.  These diagnostic queries all start
  with @ (and can co-exist in the ``.pf`` file, with normal commands. These options accessible
  with this flag are described further in :ref:`Diag`.

-e n
  Where ``n`` is a number, changes the number of errors of a specific type that
  are allowed to occur before the program gives up.  For a variety of reasons,
  errors are expected during SIROCCO runs.
  Most of these errors are harmless in the sense that they occur rarely.
  But if an error occurs too often, something is seriously and so SIROCCO halts at that point.
  The default is :math:`10^{5}` (per process).

-e_write n
  Changes the number of times an error message of a specific type is written
  to a diagnostic file.  When errors occur, a line describing the error is written
  to the diagnostic file the first ``n`` times the error occurs. After that statistics
  are maintained as to the number of times the error occurred, but it is not printed
  to the diagnostic file. The default is 100 (per process)

-classic
  Reverts to using v/c corrections for special relativity and eliminates work done to treat
  co-moving frames properly.  This is for testing, and is likely to be removed in the not
  too distant future.

-srclassic
   Use SIROCCO with full special relativity for Doppler shits, etc., but do not include any co-moving frame effects.

-no-matrix-storage
   Do not store macro-atom transition matrices if using the macro-atom line transfer and the matrix matom_transition_mode.

-ignore_partial_cells
   Ignore wind cells that are only partially filled by the wind (This is now the default)

-include_partial_cells
   Include wind cells that are only partially filled by the wind

Running Different Versions of SIROCCO
=================================

Once you have SIROCCO up and running, you can also install and run different versions of SIROCCO. This is particularly useful if you want to run and compare an older model from a previous paper or how the outputs have evolved.

You can store multiple older versions of SIROCCO by recompiling a newer version.

* Pull in the version of the program you want using git.
* Then navigate with the terminal into SIROCCO's :code:`source` folder.
* Run "make all" to recompile all programs with the new updates.
* The process will put the new binaries into the :code:`bin/` directory and not delete what was already there.

You can then run a specific installed version by replacing the SIROCCO executable, eg :code:`sirocco root.pf`, with the version you desire, eg.

  .. code :: bash

    sirocco87a root.pf