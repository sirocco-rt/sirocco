Wind.fixed_concentrations_file
==============================
The filename for the fixed ion concentrations if you have
set Wind_ionization to 2 (fixed). This file has format
[atomic_number  ionizationstage   ion fraction].


And example of a fixed concentrations file is below::

    1   1  0
    1   2  1



In the example, the only element is H, and H is completely ionized.

Note that if one wants electron densities to agree with that
expected from the ions in the in the fixed concentrations
file, then for imported models, one make sure the the elements_ions data file
only activates these elements.

Type
  String

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Wind.ionization`: fixed


