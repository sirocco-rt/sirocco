Auger Photon Yields
###################

When inner shell (or Auger) ionization takes place - there is a chance of photons being eected as the inner shells are re-filled. This data
provies the information to compute the photons thus made. It is currently not used.

Source
======
This data comes from `Kaastra and Mewe 1993, A&A, 97, 443 <http://articles.adsabs.harvard.edu/full/1993A%26AS...97..443K>`_ . The data is downloaded from the vizier site linked and put into a file called "fluorescent\_yield.data"

Translation to SIROCCO format
============================

The translation takes place using the python script "kaastra_2_py.py". All identical to electron yield, but input file is "fluorescent_yield.data" and output is "kaastra_fluorescent_yield.data"


Data format
===========

This is the data format of the electron yield data

+-----------+--+------+---+--+-------------------+-----------+
|Label      |z |state | n |l | photon_energy(eV) |yield      |
+-----------+--+------+---+--+-------------------+-----------+
|Kphotyield |5 | 1    | 1 |0 | 1.837e+02         | 6.000e-04 |
+-----------+--+------+---+--+-------------------+-----------+
|Kphotyield |5 |1     |1  |0 | 1.690e+01         | 7.129e-01 |
+-----------+--+------+---+--+-------------------+-----------+
|Kphotyield |6 |1     |1  |0 |2.768e+02          | 2.600e-03 |
+-----------+--+------+---+--+-------------------+-----------+



The data is linked to the correct inner shell photoionization cross section (and hence rate) via z, state, n shell and l subshell. The photon energy field is thew energy of the fluorescent photon in eV, and yield is the number of said photons emitted per ionization multiplied by :math:`10^4`.


SIROCCO structure
================

The data is stored in python in the inner_fluor_yield structure which contains


- int nion - Index to the ion which was the parent of the inner shell ionization
- int z, istate - element and ionization state of parent ion
- int n, l - Quantum numbers of shell
- double freq - the rest frequency of the photon emitted 
- double yield - number of photons per ionization x :math:`10^4`


Comments
========
This data is not currently used


