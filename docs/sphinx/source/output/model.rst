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

Several routines in the SIROCCO distribution allow the user to examine the wind model after a run.

windsave2table
  Executed from the command line with :code:`windsave2table rootname`.

  Produces a standard set of ASCII tables showing, for each grid cell, quantities such as wind velocity,
  :math:`n_e`, temperatures, and densities of prominent ions.  Each table has a header row identifying
  columns including ``x``, ``z``, ``xcen``, ``zcen``, ``i``, ``j`` (cylindrical coordinates and grid
  indices), followed by the quantity of interest.

  Various options control how much data is printed.  A summary can be obtained with :code:`windsave2table -h`.
  See :doc:`fits_tools` for the FITS equivalent.

swind
  Executed from the command line with :code:`swind rootname`

  Allows the user to query for information about the model interactively.  The results can be written to ASCII files for future reference.

  Various options are available; a summary can be obtained with :code:`swind -h`.

FITS tools
  Three programs export wind-save data to FITS binary tables that can be read with Python/astropy or TOPCAT:

  * :code:`windsave2fits` — per-cell spectral models (``root_cellspec.fits``)
  * :code:`windstruct2fits` — full wind structure and plasma state including both hemispheres (``root_windstruct.fits``)

  See :doc:`fits_tools` for full column documentation.
