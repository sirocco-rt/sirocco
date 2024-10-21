Auger Electron Yields
#####################

This data is linked with the inner shell photoionization data. It gives probabilities for different numbers of electrons to be ejected
following inner shell ionizations.

Source
======

This data comes from `Kaastra and Mewe 1993, A&A, 97, 443 <http://articles.adsabs.harvard.edu/full/1993A%26AS...97..443K>`_ . The data is downloaded from the vizier site linked and put into a file called "electron_yield.data"

Translation to SIROCCO format
============================

The translation takes place using the python script "kaastra_2_py.py" which takes the saved raw data file "electron_yield.data" and compares it line by line to the inner shell cross section data in "vy_innershell_tab.data"(see above). The n shell and l subshell to which each record applies is coded in the KM data and needs to be decoded. This is what the script does, and all the script then does is output the yield data into a new file "kaastra_electron_yield.data" which contains the n and l cross reference.


Data format
===========

This is the data format of the electron yield data

+-----------+---+-------+---+--+------------+------------+-----------+--------+
|Label      | z | state | n |l |IP(eV)      | <E>(eV)    | P(1e)     | P(2e)  |
+-----------+---+-------+---+--+------------+------------+-----------+--------+
|Kelecyield | 4 |1      |1  |0 |1.15e+02    | 9.280e+01  | 0         | 10000  | 
+-----------+---+-------+---+--+------------+------------+-----------+--------+
|Kelecyield | 5 |1      |1  |0 |1.92e+02    | 1.639e+02  | 6         | 9994   |
+-----------+---+-------+---+--+------------+------------+-----------+--------+
|Kelecyield | 5 |2      |1  |0 |2.06e+02    | 1.499e+02  | 0         | 10000  |
+-----------+---+-------+---+--+------------+------------+-----------+--------+



The data is linked to the correct inner shell photoionization cross section (and hence rate) via z, state, n shell and l subshell. The IP is not used. <E>  is the mean electron energy ejected, used to calculate the PI heating rate in radiation.c. The last ten columns in the file (2 shown in the table above) show the chance of various numbers of electrons being ejected in units of 1/10000. 


SIROCCO structure
================

The data is stored in python in the inner_elec_yield structure which contains

- int nion - Index to the ion which was the parent of the inner shell ionization
- int z, istate - element and ionization state of parent ion
- int n, l - Quantum numbers of shell
- double prob[10] - probability for between 1 and 10 electrons being ejected 
- double I - Ionization energy
- double Ea - Average electron energy
 
Comments
========

