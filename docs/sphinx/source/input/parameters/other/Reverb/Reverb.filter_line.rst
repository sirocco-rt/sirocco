Reverb.filter_line
==================
Line number of one line to include in the output ``.delay_dump`` file. This is
the sirocco internal line number. It can be found using either the macro-atom
mode (which prints out the line number once it's found one) or by doing an
exploratory run with :ref:`reverb.filter_lines` = -1, then looking through the delay
dump file for photons of the right wavelength to see what their line is. This
should almost certainly be changed to be specified using a species and
wavelength!

Currently, the easiest way to find the line is to use the `Atomix <https://github.com/saultyevil/atomix>`_ tool to
search for the line number.

Type
  Integer

Values
  Any valid line index

File
  `setup_reverb.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_reverb.c>`_


Parent(s)
  * :ref:`Reverb.filter_lines`: Greater than 0, once per filer line.


