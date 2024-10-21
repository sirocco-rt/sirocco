Bound Free (Verner)
###################

This is data for bound free or photoionization data. There is information for both inner shell (auger) and outer shell PI.


Source
======

There are three sources for this data


- `Verner & Yakovlev 1995 <http://adsabs.harvard.edu/abs/1995A\%26AS..109..125V>`_ : Inner and Outer Shell Cross-sections
- `Verner et al. 1996 <http://adsabs.harvard.edu/abs/1996ApJ...465..487V>`_ :Improved Outer Shell Cross-sections
- `Kaastra \& Mewe 1993 <http://adsabs.harvard.edu/abs/1993A\%26AS...97..443K>`_ :Electron and photon yield data


Translation to SIROCCO format
============================

**Tabulation Process**

The raw VFKY data comes in a series of fit parameters. We decided, circa SIROCCO 78, to tabulate this data instead. Partly, I think I because the on the fly method was time consuming (yes, computing all the pow() commands to commute the cross sections on the fly took a huge amount of time) and we decided that tabulating pre program was better than doing it in the program, so that everything was of the same format.

The script which does this is progs/tabulate\_xs/photo\_xs.py -- it creates a file like photo\_vfky\_tabulated.data.

**Inner and Outer Shells**

For the ground states, we split the cross sections up into outer shell and inner shell cross sections. This allows us to calculate possible auger ionization as ions relax after an inner shell ionization. This is done using the python script "verner_2_py.py. This script takes the normal verner cross sections, which truncate at the first inner shell edge and firstly appends the outer shell data from VY to that to make a full outer shell cross section. These are written out into "vfky_outershell_tab.data"
It then writes out the inner shell cross sections into "vfky_innershell_tab.data". There is a lot of complicated machinery to try and work out the exact shell that is being ionized to allow these rates to be linked up to the relevant electron yield (and flourescent) records.




Data format
===========

Explain the ascii format of the file which is read into SIROCCO

**VFKY_outershell_tab.data**

+----------+--+------+------+------+-----------------+---------+
|Label     |z |state |islp  |level |threshold_energy |n_points |
+----------+--+------+------+------+-----------------+---------+
|PhotVfkyS | 1| 1    | -999 | -999 | 1.360e+01       | 100     |
+----------+--+------+------+------+-----------------+---------+



This data is linked to the relevant ion via z and state, islp and level are not used. the last number n_points, says how many points are used for the fit, and the next n_points lines in the file, each preceeded by the label PhotVfky are pairs of energy (in eV) vs cross section which make up that fit.

**VY_innershell_tab.data**

+---------+--+------+--------+------------+------------------+----------+
|label    |z |state |n_shell | l_subshell | threshold_energy | n_points |
+---------+--+------+--------+------------+------------------+----------+
|InnerVYS |3 |1     |1       |0           |  6.439e+01       | 100      |
+---------+--+------+--------+------------+------------------+----------+



This data is linked to the relevant ion via z and state. the n_shell and l_subshell numbers are used to cross reference to the electron yield records. As above, the last record shows how many points are in the fit, and the data pairs making up the fit follow with the keyword InnerVY.

SIROCCO structure
================

Where the data is stored internally in SIROCCO


Comments
========

The manner in which this data is read into SIROCCO is a bit labyrinthine at the moment. The intention is to use a combination of VFKY and VY for all ground states, an
