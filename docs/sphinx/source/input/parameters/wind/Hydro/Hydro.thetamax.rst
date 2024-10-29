Hydro.thetamax
==============
The maximum theta value to be read in from a hydrodynamic snapshot.
This is typically used to excise a dense disk from the midplane of
such a snapshot. Use a negative value to tell the code to use all
the data.

Type
  Double

Unit
  Degrees

Values
  -1
    use all data

  X
    use up to that angle (typically less than 90)


File
  `hydro_import.c <https://github.com/sirocco-rt/sirocco/blob/master/source/hydro_import.c>`_


Parent(s)
  * :ref:`Wind.type`: ``hydro``


