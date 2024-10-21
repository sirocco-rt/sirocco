Evaluation
###########

Determining whether SIROCCO has run successfully from a a scientific point of view depends very specifically on one's goals. Did the spectra turn out to be what one expected? Here by evaluation we mean, did my run complete without significant errors and did the ionization structure converge to a "steady state" given the number of ionization cycles, the number of photons, and the frequency distributions of the photons we chose.

Convergence
============

:doc:`Ionization cycles <../operation/ionization_cycles>` in SIROCCO are intended to establish the correct ion densities and temperature of the various cells in the wind.  The degree to which this happens for a given number of ionization cycles depends on how far the initial guess of electrons temperatures in various portions of the wind and the number of photons generated during each photoionization cycles.  Furthermore, the accuracy of the final model depends on the number of photons that pass through each cell.  As a result, the accuracy with which ion abundances and temperature are determined  will differ on a cell by cell basis. In a typical model with a biconical outflow, a small cells at the outer edge of the accretion disk will record fewer photon passages than one in the middle of the grid that is exposed to large numbers of photons from the disk.

A very basic question about a particular run is, has it reached a "steady state" and if it is in a steady state are cells stable in the sense that fluctuations are small. Hopefully, each ionization cycle brings one closer and closer to to a solution after which the ionization structure no longer evolves. Of course, since SIROCCO is a Monte Carlo code,  the degree to which the solution stays constant from cycle to cycle is limited by counting statistics.  To check convergence in SIROCCO, we monitor the the radiation temperature :math:`T_r`, the electron temperature :math:`T_e`, and the total heating :math:`\mathrm{Heat_{tot}}` and cooling :math:`\mathrm{Cooling_{tot}}` in each cell as a function of the ionization cycle :math:`n`. 

To estimate whether a a model calculation has reached a steady state, SIROCCO carries out three tests, one comparing the difference in the radiation temperature between the current and the preceding ionization cycle, one comparing the electron temperature in the same manner and once comparing heating and cooling rates in the current cycle. If a cell satisfies the following 3 tests,

.. math::
    \left | \frac{T_e^n-T_e^{n-1}}{T_e^n+T_e^{n-1}} \right | < \epsilon

.. math::
    \left | \frac{T_r^n-T_r^{n-1}}{T_r^n+T_r^{n-1}} \right | < \epsilon

.. math::
    \left | 
    \frac{\mathrm{Heat_{tot}}^n- \mathrm{Cooling_{tot}}^{n}}
    {\mathrm{Heat_{tot}}^n + \mathrm{Cooling_{tot}}^{n}} 
    \right | <\epsilon

where :math:`\epsilon = 0.05`, it is flagged as converged. 

It is rare that all of the cells in a model will satisfy all of these criteria.  That is is because the number photons passing that pass through a cell vary considerably and the statistical  noise depends on the the number of photons. It is important to note that the photons that contribute most to the spectra of an object will be those which have the most photons passing through them.  Typically, we consider a model to be well converged if 90% of the cells are converged by this criterion.

The routine ``run_check.py`` (see :doc:`SIROCCO Script documentation <../py_progs>`) produces two plots related to convergence, one showing the fraction of cells that have passed each of the tests above as a function of cycle, and the other showing the number of failed checks for each cell in the wind.

Note that it is not always important that all cells be converged. The Monte Carlo process preferentially picks out the cells which affect the emergent radiation field. Portions of the grid which do not get many photons are typically the ones that are "not converged", but since they don't contribute much to the emergent radiation, one does not need them to be converged (except if one wants to make nice plots of the temperature as a function of position in the wind or of the density of a particular ion species). On the other hand, if one is using SIROCCO in conjunction with a hydrodynamical code one wants all the cells to be converged.

Errors
============

SIROCCO is designed to continue to run unless something catastrophic happens.
As it runs, it logs error messages that can be found in the .diag files.
These messages are a combinations of warnings, and/or unusual occurrences,
that if they start occurring often suggest a real problem.

These error messages are all of the form:

.. code::

   Error: wind2d: Cell   0 ( 0, 0) in domain 0 has 1 corners in wind, but zero volume

that is they begin with the word :code:`Error`. followed by the subroutine in the code where the error occurred followed by what is hopefully a helpful.
If one is concerned about a particular message, one can hopefully determine what is happening by looking for that message in the log files.

SIROCCO keeps a count of the number of times a particular message has occurred and at the end of the program, and the very end of the
diag files contain a listing of how many times a particular error has occurred.

.. code::

   Error summary: End of program, Thread 2 only
   Recurrences --  Description
        7 -- getatomic_data: line input incomplete: %s
      128 -- get_atomicdata: Could not interpret line %d in file %s: %s
        1 -- error_count: This error will no longer be logged: %s
        1 -- get_wind_params: zdom[ndom].rmax 0 for wind type %d
        1 -- wind2d: Cell %3d (%2d,%2d) in domain %d has %d corners in wind, but zero volume
        1 -- check_grid: velocity changes by >1,000 km/s in %i cells
        1 -- check_grid: some cells have large changes. Consider modifying zlog_scale or grid dims

As indicated here, these are the errors for only thread 2 of a program.
In order to get a summary of all the threads, there is a script py_error.py that be run as :code:`py_error.py rootname` from the main run directory.
Note that in many cases, the summary will be the number times an error occurred in one thread times the number of threads, but not always.

One should be aware of these errors, and watch out for situations where the number of errors  of a particular type is much larger than usual.
