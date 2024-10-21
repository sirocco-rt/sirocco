Levels
######


Once the element and ion data has been read into \textsc{SIROCCO}, the next step is to read in the level information.

Source:
=======

Level information can be derived from a variety of sources, by:

   * deriving it from a line list such as that of Kurucz, or
   * more commonly, for data abstracted from TopBase or Chianti


Translation to python:
======================



Various Formats of Level Files
==============================

The original format::

 Comment-- Ion H  1
 # There are 17 unique levels for 1 1
 Level   1   1   0    2   0.000000
 Level   1   1   1    2  10.200121
 Level   1   1   2    4  10.200166
 Level   1   1   3    2  12.089051

where the colums are z, istate (in conventional notation), a unique level no,
the multiplicity of the level and the excitation energy of the level in eV

Level files with this type of format are used in  files such as levels\_kur.dat, when the
levels are derived from linelists such as Kurucz.  It is only allowable when dealing
with simple atoms.  

The current format::

  # Maximum excitation above ground (in Ryd)  for inclusion  4.000000
  # Miniumum excitation below continuum (in Ryd) for inclusion -0.100000
  # Topbase levels: Order changed to move config to end of line
  #LevTop  z ion  iSLP iLV iCONF   E(eV) Te(eV) gi RL(s) eqn RL(s)
  # ======================================================================================
  #      i NZ NE iSLP iLV iCONF                 E(RYD)      TE(RYD)   gi     EQN    RL(NS)
  # ======================================================================================
  # ======================================================================================
  LevTop  1  1  200  1 -13.605698   0.000000  2  1.0000 1.00e+21 () 1s
  LevTop  1  1  200  2  -3.401425  10.204273  2  2.0000 1.00e+21 () 2s
  LevTop  1  1  211  1  -3.401425  10.204273  6  2.0000 1.60e-09 () 2p

whereas the for Macro Atoms we have::

  #         z ion lvl ion_pot   ex_energy   g  rad_rate
  LevMacro  1  1  1 -13.59843    0.000000   2  1.00e+21 () n=1
  LevMacro  1  1  2  -3.39960    10.19883   8  1.00e-09 () n=2
  LevMacro  1  1  3  -1.51093    12.08750  18  1.00e-09 () n=3


The columns are similar in the two cases.

Each level is described by an element number and ion number and a level number.  
In the macro-atom case the level number is unique; in the simple atom case the combination of iSLP and the level number are unique.  
 
For the Topbase case for simple atoms, the columns are:
 
 * the atomic number
 * the ion number in the usual astronomical convention
 * iSLP 
 * the level number
 * the energy in eV relative to the continuum
 * the energy in eV relative to the ground state
 * the multiplicity (g) of the level, usually 2J+1
 * the equivalent quantum number
 * the radiative lifetime
 * the configuration
 
There are some specific differences. 
In particular, for LevMacro levels, the excitation energy needs need to be on an absolute scale 
between ions, and so it includes the ionization energy of the lower level ionization states. 
Note that the radiative rates are not used. The original intention was to use this to define the 
difference between metastable and normal levels, with the expectation that if the level was metastable it 
would be put in Boltzmann equilibrium with the ground state. 
Right now python uses :math:`10^{15}` seconds, essentially a Hubble time to do this, but this portion of the 
code is not, according to ss, tested. 

The primary source for this is usually the NIST database, although similar information is usually available in Chianti. 
One normally wants text output, and eV to describe the levels, and then one needs to put things in energy order. 
Since they quote J, one converts to g = 2J+1


The ionization potential is not used, as it is redundant with the excitation energy which is, and the last column giving the configuration is also for information only.

SIROCCO structure:
=================
This data is held in SIROCCO in various fields in structure **config**.

Comments:
=========

