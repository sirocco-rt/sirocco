Diag.make_ioncycle_tables
=========================
This diagnostic switch causes a series of astropy-compatable tables being 
written for each ion cycle into the main diag directory for a run.  The tables are identical to those that can be generated with the separate routine windsave2table, and are intended to provide an overview of the ionization structure and the importance of various physical processes.

Additionally, this option also causes information about the disk into a
disk diag file after each ionization cycle concerning how much
energy is hitting the disk as a function of disk radius.  

Type
  Boolean (yes/no)

File
  `diag.c <https://github.com/agnwinds/python/blob/master/source/diag.c>`_


Parent(s)
  * :ref:`Diag.extra`: ``True``


