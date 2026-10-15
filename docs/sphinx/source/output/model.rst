Model
#####

As SIROCCO is run, it repeatedly writes out two binary files that contain essentially all information about the wind as calculated in the ionization phase of the program,
along with status of the program at the last point where the file was written.
These files along with the parameter file are sufficient to restart the program,
if for example, one wants to check point the program after a certain time, and restart where one left off,
or to add spectral cycles to get better spectra.

.wind_save
  A binary file that contains essentially all information about the wind including ion densities,
  temperatures, and velocities in each cell, along with status of the program at the last point where the file was written.

.spec_save
  A binary file that contains all of the information about the spectra that have created.  This file is not of interest to users directly.  It is used when restarting

Several standalone C programs are included in the distribution for inspecting and
post-processing wind save files.  See :doc:`/c_executables` for full
documentation of all programs.  A brief summary of the most commonly used ones:

windsave2table
  Executed from the command line with :code:`windsave2table rootname`.

  Produces a standard set of ASCII tables showing per-cell quantities such as
  wind velocity, :math:`n_e`, temperatures, and ion densities.  Run
  :code:`windsave2table -h` for a full list of options.

swind
  Executed from the command line with :code:`swind rootname`.

  Allows the user to query for information about the model interactively.
  Results can be written to ASCII files.  Run :code:`swind -h` for options.
